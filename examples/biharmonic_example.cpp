/** @file biharmonic_example.cpp

    @brief A Biharmonic example.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

# include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot = false;

    index_t numRefine  = 3;
    index_t discreteDegree = 3;
    index_t discreteRegularity = 2;
    bool last = false;
    bool neumann = false;

    std::string fn;
    index_t geometry = 1000;

    gsCmdLine cmd("Example for solving the biharmonic problem.");
    cmd.addInt( "p", "discreteDegree","Which discrete degree?", discreteDegree );
    cmd.addInt( "r", "discreteRegularity", "Number of discreteRegularity",  discreteRegularity );
    cmd.addInt( "l", "refinementLoop", "Number of refinementLoop",  numRefine );
    cmd.addString( "f", "file", "Input geometry file", fn );
    cmd.addInt( "g", "geometry", "Which geometry",  geometry );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    cmd.addSwitch("neumann", "Neumann", neumann);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mp;
    gsBoundaryConditions<> bc;
    gsFunctionExpr<> f, ms;
    gsOptionList optionList;

    //! [Read geometry]
    std::string string_geo;
    if (fn.empty())
        string_geo = "planar/geometries/g" + util::to_string(geometry) + ".xml";
    else
        string_geo = fn;

    gsInfo << "Filedata: " << string_geo << "\n";
    gsReadFile<>(string_geo, mp);
    mp.clearTopology();
    mp.computeTopology();

    gsFunctionExpr<>source("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    f.swap(source);
    gsInfo << "Source function " << f << "\n";

    gsFunctionExpr<> solution("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
    ms.swap(solution);
    gsInfo << "Exact function " << ms << "\n";

    //! [Refinement]
    gsMultiBasis<> dbasis(mp, false);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( discreteDegree); // preserve smoothness

    // h-refine each basis
    if (last)
    {
        for (int r =0; r < numRefine; ++r)
            dbasis.uniformRefine(1, discreteDegree-discreteRegularity);
        numRefine = 0;
    }

    //! [Boundary condition]
    for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
        // Laplace
        gsFunctionExpr<> laplace ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);

        // Neumann
        gsFunctionExpr<> sol1der("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                                 "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)", 2);

        bc.addCondition(*bit, condition_type::dirichlet, ms);
        if (neumann)
            bc.addCondition(*bit, condition_type::neumann, sol1der);
        else
            bc.addCondition(*bit, condition_type::laplace, laplace);
    }
    bc.setGeoMap(mp);
    gsInfo << "Boundary conditions:\n" << bc << "\n";
    //! [Boundary condition]

    //! [Problem setup]
    gsExprAssembler<real_t> A(1,1);
    //gsInfo<<"Active options:\n"<< A.options() <<"\n";

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<real_t> ev(A);

    // Set the geometry map
    auto G = A.getMap(mp);

    // Set the source term
    auto ff = A.getCoeff(f, G); // Laplace example

    // Set the discretization space
    gsMappedBasis<2,real_t> bb2;
    auto u = A.getSpace(dbasis);

    // Solution vector and solution variable
    gsMatrix<real_t> solVector;
    auto u_sol = A.getSolution(u, solVector);

    // Recover manufactured solution
    auto u_ex = ev.getVariable(ms, G);
    //! [Problem setup]

    //! [Solver loop]
    gsSparseSolver<real_t>::SimplicialLDLT solver;

    gsVector<real_t> l2err(numRefine+1), h1err(numRefine+1), h2err(numRefine+1),
    IFaceErr(numRefine+1), dofs(numRefine+1);
    gsInfo<< "(dot1=approxC1construction, dot2=assembled, dot3=solved, dot4=got_error)\n"
        "\nDoFs: ";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);
    gsStopwatch timer;
    for (int r=0; r<=numRefine; ++r)
    {
        dbasis.uniformRefine(1,discreteDegree -discreteRegularity);

        gsInfo << dbasis.basis(0).getMaxCellLength() << " ";

        // Setup the system
        u.setup(bc, dirichlet::l2Projection, 0);

        // Initialize the system
        A.initSystem();
        setup_time += timer.stop();

        gsInfo<< A.numDofs() <<std::flush;

        timer.restart();
        // Compute the system matrix and right-hand side
        A.assemble(ilapl(u, G) * ilapl(u, G).tr() * meas(G),u * ff * meas(G));

        // Enforce Laplace conditions to right-hand side
        auto g_L = A.getBdrFunction(G); // Set the laplace bdy value
        //auto g_L = A.getCoeff(laplace, G);
        A.assembleBdr(bc.get("Laplace"), (igrad(u, G) * nv(G)) * g_L.tr() );

        dofs[r] = A.numDofs();
        ma_time += timer.stop();
        gsInfo<< "." <<std::flush;// Assemblying done

        timer.restart();
        solver.compute( A.matrix() );
        solVector = solver.solve(A.rhs());

        slv_time += timer.stop();
        gsInfo<< "." <<std::flush; // Linear solving done

        timer.restart();
        //linferr[r] = ev.max( f-s ) / ev.max(f);

        l2err[r]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) ); // / ev.integral(f.sqNorm()*meas(G)) );
        h1err[r]= l2err[r] +
            math::sqrt(ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) )); // /ev.integral( igrad(f).sqNorm()*meas(G) ) );

        h2err[r]= h1err[r] +
                 math::sqrt(ev.integral( ( ihess(u_ex) - ihess(u_sol,G) ).sqNorm() * meas(G) )); // /ev.integral( ihess(f).sqNorm()*meas(G) )

        err_time += timer.stop();
        gsInfo<< ". " <<std::flush; // Error computations done
    } //for loop

    //! [Error and convergence rates]
    gsInfo<< "\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";
    gsInfo<< "H2 error: "<<std::scientific<<h2err.transpose()<<"\n";
    gsInfo<< "Deriv Interface error: "<<std::scientific<<IFaceErr.transpose()<<"\n";
    gsInfo<< "Dofs: "<<std::fixed<<dofs.transpose()<<"\n";

    if (numRefine>0)
    {
        gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
              <<  ( l2err.head(numRefine).array()  /
                   l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0)
                   <<"\n";

        gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
              <<( h1err.head(numRefine).array() /
                  h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";

        gsInfo<<   "EoC (H2): "<< std::fixed<<std::setprecision(2)
              <<( h2err.head(numRefine).array() /
                  h2err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";

        gsInfo<<   "EoC (Iface): "<< std::fixed<<std::setprecision(2)
              <<( IFaceErr.head(numRefine).array() /
                  IFaceErr.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    }
    //! [Error and convergence rates]

    // Plot solution in paraview
    if (plot)
    {
        // Write approximate and exact solution to paraview files
        gsInfo<<"Plotting in ParaView...\n";
        gsInfo<<"Plotting in Paraview...\n";
        ev.options().setSwitch("plot.elements", true);
        ev.writeParaview( u_sol   , G, "solution");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";

    return  0;
}
