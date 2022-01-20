/** @file biharmonic2_example.cpp

    @brief Tutorial on how to use expression assembler and the (approx.) C1 basis function
                to solve the Biharmonic equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

//! [Include namespace]
#include <gismo.h>

#include <gsUnstructuredSplines/gsApproxC1Spline.h>
//#include <gsUnstructuredSplines/gsDPatch.h>

#include <gsIO/gsCSVWriter.h>

using namespace gismo;
//! [Include namespace]

void setMapperForBiharmonic(gsBoundaryConditions<> & bc, gsMappedBasis<2,real_t> & bb2, gsDofMapper & mapper)
{
    mapper.setIdentity(bb2.nPatches(), bb2.size(), 1);

    gsMatrix<index_t> bnd;
    for (typename gsBoundaryConditions<real_t>::const_iterator
                 it = bc.begin("Dirichlet"); it != bc.end("Dirichlet"); ++it)
    {
        bnd = bb2.basis(it->ps.patch).boundary(it->ps.side());
        mapper.markBoundary(it->ps.patch, bnd, 0);
    }

    for (typename gsBoundaryConditions<real_t>::const_iterator
             it = bc.begin("Neumann"); it != bc.end("Neumann"); ++it)
    {
        bnd = bb2.basis(it->ps.patch).boundaryOffset(it->ps.side(),1);
        mapper.markBoundary(it->ps.patch, bnd, 0);
    }

    mapper.finalize();
}

void gsDirichletNeumannValuesL2Projection2(gsMultiPatch<> & mp, gsMultiBasis<> & dbasis, gsBoundaryConditions<> & bc,
                                           gsMappedBasis<2,real_t> & bb2, const expr::gsFeSpace<real_t> & u)
{
    gsDofMapper mapper = u.mapper();

    gsMatrix<index_t> bnd = mapper.findFree(mapper.numPatches()-1);
    gsDofMapper mapperBdy;
    mapperBdy.setIdentity(bb2.nPatches(), bb2.size(), 1);
    mapperBdy.markBoundary(0, bnd, 0);
    mapperBdy.finalize();

    gsExprAssembler<> A(1,1);
    A.setIntegrationElements(dbasis);

    auto G = A.getMap(mp);
    auto uu = A.getSpace(bb2);
    auto g_bdy = A.getBdrFunction(G);

    uu.setupMapper(mapperBdy);
    gsMatrix<real_t> & fixedDofs_A = const_cast<expr::gsFeSpace<real_t>&>(uu).fixedPart();
    fixedDofs_A.setZero( uu.mapper().boundarySize(), 1 );

    real_t lambda = 1e-5;

    A.initSystem();
    A.assembleBdr(bc.get("Dirichlet"), uu * uu.tr() * meas(G));
    A.assembleBdr(bc.get("Dirichlet"), uu * g_bdy * meas(G));
    A.assembleBdr(bc.get("Neumann"),
                  lambda * (igrad(uu, G) * nv(G).normalized()) * (igrad(uu, G) * nv(G).normalized()).tr() * meas(G));
    A.assembleBdr(bc.get("Neumann"),
                  lambda *  (igrad(uu, G) * nv(G).normalized()) * (g_bdy.tr() * nv(G).normalized()) * meas(G));

    gsSparseSolver<>::SimplicialLDLT solver;
    solver.compute( A.matrix() );
    gsMatrix<real_t> & fixedDofs = const_cast<expr::gsFeSpace<real_t>& >(u).fixedPart();
    fixedDofs = solver.solve(A.rhs());
}


int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t smoothing = 2;

    index_t numRefine  = 3;
    index_t discreteDegree = 3;
    index_t discreteRegularity = 2;
    bool last = false;
    bool info = false;
    bool neumann = false;
    bool nitsche = false;

    std::string xml;
    bool output = false;

    std::string fn;

    index_t geometry = 1000;

    gsCmdLine cmd("Tutorial on solving a Biharmonic problem.");
    cmd.addInt( "s", "smoothing","Smoothing", smoothing );
    cmd.addInt( "p", "discreteDegree","Which discrete degree?", discreteDegree );
    cmd.addInt( "r", "discreteRegularity", "Number of discreteRegularity",  discreteRegularity );
    cmd.addInt( "l", "refinementLoop", "Number of refinementLoop",  numRefine );
    cmd.addString( "f", "file", "Input geometry file", fn );
    cmd.addInt( "g", "geometry", "Which geometry",  geometry );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("info", "Getting the information inside of Approximate C1 basis functions", info);

    cmd.addSwitch("neumann", "Neumann", neumann);
    cmd.addSwitch("nitsche", "Nitsche", nitsche);

    cmd.addString("x", "xml", "Use the information from the xml file", xml);
    cmd.addSwitch("output", "Output in xml (for python)", output);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Initialize data]
    gsMultiPatch<> mp;
    gsBoundaryConditions<> bc;
    gsFunctionExpr<> f, ms;
    gsOptionList optionList;
    //! [Initialize data]

    //! [Read Argument inputs]
    if (xml.empty()) {
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

        optionList = cmd;
        gsInfo << "OptionList: " << optionList << "\n";
        gsInfo << "Finished\n";
    }
    //! [Read Argument inputs]

    //! [Read XML file]
    else
    {
        // id=0 Boundary
        // id=1 Source function
        // id=2 Optionlist
        // id=3 Exact solution
        // id=X Geometry (should be last!)
        gsFileData<> fd(xml); // "planar/biharmonic_pde/bvp1.xml"

        // Geometry
        fd.getAnyFirst(mp);
        mp.computeTopology();
        gsInfo << "Multipatch " << mp << "\n";

        // Functions
        fd.getId(1, f); // Source solution
        gsInfo << "Source function " << f << "\n";

        fd.getId(3, ms); // Exact solution
        gsInfo << "Exact function " << ms << "\n";

        // Boundary condition
        fd.getId(0, bc); // id=2: boundary conditions
        bc.setGeoMap(mp);
        gsInfo << "Boundary conditions:\n" << bc << "\n";

        // Option list
        fd.getId(2, optionList); // id=100: assembler options
        gsInfo << "OptionList: " << optionList << "\n";
    }
    //! [Read XML file]

    //! [Read option list]
    discreteDegree = optionList.getInt("discreteDegree");
    discreteRegularity = optionList.getInt("discreteRegularity");
    numRefine = optionList.getInt("refinementLoop");

    nitsche = optionList.getSwitch("nitsche");

    plot = optionList.getSwitch("plot");
    info = optionList.getSwitch("info");
    //! [Read option list]

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

    // Assume that the condition holds for each patch TODO
    // Refine once
    if (dbasis.basis(0).numElements() < 4)
        dbasis.uniformRefine(1, discreteDegree-discreteRegularity);


    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    //gsInfo<<"Active options:\n"<< A.options() <<"\n";

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    auto G = A.getMap(mp);

    // Set the source term
    auto ff = A.getCoeff(f, G); // Laplace example

    // Set the discretization space
    gsMappedBasis<2,real_t> bb2;
    auto u = nitsche ? A.getSpace(dbasis) : A.getSpace(bb2);

    // The approx. C1 space
    gsSparseMatrix<real_t> global2local;
    gsApproxC1Spline<2,real_t> approxC1(mp,dbasis);
    approxC1.options().setSwitch("info",info);
    approxC1.options().setSwitch("plot",plot);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    auto u_sol = A.getSolution(u, solVector);

    // Recover manufactured solution
    auto u_ex = ev.getVariable(ms, G);
    //! [Problem setup]

    //! [Solver loop]
    gsSparseSolver<>::SimplicialLDLT solver;

    gsVector<> l2err(numRefine+1), h1err(numRefine+1), h2err(numRefine+1),
    IFaceErr(numRefine+1);
    gsInfo<< "(dot1=approxC1construction, dot2=assembled, dot3=solved, dot4=got_error)\n"
        "\nDoFs: ";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);
    gsStopwatch timer;
    for (int r=0; r<=numRefine; ++r)
    {
        dbasis.uniformRefine(1,discreteDegree -discreteRegularity);

        if (!nitsche)
            approxC1.update(bb2);
        gsInfo<< "." <<std::flush; // Approx C1 construction done

        // Setup the mapper
        if (!nitsche) // MappedBasis
        {
            gsDofMapper map;
            setMapperForBiharmonic(bc, bb2,map);

            // Setup the system
            u.setupMapper(map);
            gsDirichletNeumannValuesL2Projection2(mp, dbasis, bc, bb2, u);
        }
        else // Nitsche
        {
            // Setup the system
            u.setup(bc, dirichlet::user, 0);
        }

        //gsMatrix<real_t> u_fixed_new = u.fixedPart();
        //gsDirichletNeumannValuesL2Projection(u, bc); // Maybe too much memory
        //gsDirichletValuesL2Projection(u, bc);
        //gsInfo<< "Error " << (u.fixedPart()-u_fixed_new).norm() <<"\n";

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

        // Enforce Neumann conditions to right-hand side
        //A.assembleRhsBc(ilapl(u, G) * (igrad(u_ex) * nv(G)).tr(), bc.neumannSides() );

        real_t penalty  = 4 * ( dbasis.maxCwiseDegree() + dbasis.dim() ) * ( dbasis.maxCwiseDegree() + 1 );
        real_t m_h      = dbasis.basis(0).getMinCellLength(); //*dbasis.basis(0).getMinCellLength();
        real_t mu       = 2 * penalty / m_h;

        if (nitsche)
        {
            gsInfo << "mu: " << mu << "\n";
            /*
            A.assembleIfc(mp.interfaces(), - 0.5 * ((igrad(u.left(), G.left()) - igrad(u.right(), G.right())) * nv(G)) *
                    (ilapl(u.left(), G.left()) + ilapl(u.right(), G.right()))
                    - 0.5 * ((igrad(u.right(), G.right()) - igrad(u.left(), G.left())) * nv(G)) *
                    (ilapl(u.right(), G.right()) + ilapl(u.left(), G.left()))
                    * meas(G)
            );
             */
            real_t alpha = 1;
            //A.assembleIfc(mp.interfaces(), - 0.5 * ((igrad(u.left(), G.left()) * nv(G)) * ilapl(u.left(), G.left()).tr()) * meas(G),
            //A.assembleIfc(mp.interfaces(), + 0.5 * ((igrad(u.right(), G.right()) * nv(G)) * ilapl(u.right(), G.right()).tr()) * meas(G));

            //A.assembleIfc(mp.interfaces(), - 0.5 * ((igrad(u.left(), G.left()) * nv(G.left())) * ilapl(u.right(), G.right()).tr()) * meas(G));
            //A.assembleIfc(mp.interfaces(), + 0.5 * ((igrad(u.right(), G.right()) * nv(G.right())) * ilapl(u.left(), G.left()).tr()) * meas(G));

            //A.assembleIfc(mp.interfaces(), alpha * ((igrad(u.left(), G.left()) * nv(G)) * (igrad(u.left(), G.left()) * nv(G)).tr()) * meas(G));
            //A.assembleIfc(mp.interfaces(), alpha * ((igrad(u.right(), G.right()) * nv(G)) * (igrad(u.right(), G.right()) * nv(G)).tr()) * meas(G));

            //A.assembleIfc(mp.interfaces(), - alpha * ((igrad(u.left(), G.left()) * nv(G)) * (igrad(u.right(), G.right()) * nv(G)).tr()) * meas(G));
            //A.assembleIfc(mp.interfaces(), - alpha * ((igrad(u.right(), G.right()) * nv(G)) * (igrad(u.left(), G.left()) * nv(G)).tr()) * meas(G));

            A.assembleIfc(mp.interfaces(),
                     //B11
                     -alpha*0.5*igrad( u.left() , G.left()) * nv(G).normalized() * (ilapl(u.left(), G.left())).tr()  * meas(G),
                     //B12
                     -alpha*0.5*igrad( u.left()  , G.left()) * nv(G).normalized() * (ilapl(u.right(), G.right())).tr() * meas(G),
                     //B21
                     alpha*0.5*igrad( u.right(), G.right()) * nv(G).normalized() * (ilapl(u.left(), G.left())).tr() * meas(G),
                     //B22
                     alpha*0.5*igrad( u.right() , G.right()) * nv(G).normalized() * (ilapl(u.right(), G.right())).tr() * meas(G),
/*
                     // symmetry
                     beta *0.5*igrad( u_dg.right(), G) * nv(G).normalized() * u_dg.right().tr() * meas(G),
                     -beta *0.5*igrad(u_dg.right(), G) * nv(G).normalized() * u_dg.left() .tr() * meas(G),
                     beta *0.5*igrad( u_dg.left() , G) * nv(G).normalized() * u_dg.right().tr() * meas(G),
                     -beta *0.5*igrad(u_dg.left() , G) * nv(G).normalized() * u_dg.left() .tr() * meas(G),
*/
                     // E11
                      mu * igrad(u.left(), G.left()) * nv(G).normalized() * (igrad(u.left(), G.left()) * nv(G).normalized()).tr() * meas(G),
                     //-E12
                      -mu * (igrad(u.left(), G.left()) * nv(G).normalized()) * (igrad(u.right(), G.right()) * nv(G).normalized()).tr() * meas(G),
                     //-E21
                      -mu * (igrad(u.right(), G.right()) * nv(G).normalized()) * (igrad(u.left(), G.left()) * nv(G).normalized()).tr() * meas(G),
                     // E22
                      mu * igrad(u.right(), G.right()) * nv(G).normalized() * (igrad(u.right(), G.right()) * nv(G).normalized()).tr() * meas(G)
                );
        }

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

        if (!nitsche)
        {
            gsMatrix<real_t> solFull;
            u_sol.extractFull(solFull);
            gsMappedSpline<2, real_t> mappedSpline(bb2, solFull);

            auto ms_sol = A.getCoeff(mappedSpline);

            IFaceErr[r] = math::sqrt(ev.integralInterface(((igrad(ms_sol.left(), G.left()) -
                                                            igrad(ms_sol.right(), G.right())) *
                                                           nv(G).normalized()).sqNorm() * meas(G)));
        }

        //gsMatrix<> points(2,1);
        //points << 1, 0.5;
        //std::vector<gsMatrix<>> result;
        //bb2.evalAllDers_into(0, points, 1, result);
        //gsDebugVar(result[1].reshape(2,result[1].rows()/2));
        //points << 0, 0.5;
        //bb2.evalAllDers_into(1, points, 1, result);
        //gsDebugVar(result[1].reshape(2,result[1].rows()/2));

        //gsDebugVar(ev.integralInterface((igrad(u_sol.left(),G) * nv(G).normalized()).sqNorm() * meas(G)));
        //gsDebugVar(ev.integralInterface((igrad(u_sol.right(),G) * nv(G).normalized()).sqNorm() * meas(G)));

        //gsDebugVar(ev.integralInterface((igrad(u_sol.left(),G) * nv(G).normalized()
        //- igrad(u_sol.right(),G) * nv(G).normalized()).sqNorm() * meas(G)));

        //gsDebugVar(ev.integralInterface(((igrad(u_sol.left(),G) - igrad(u_sol.right(),G)) * nv(G).normalized()).sqNorm() * meas(G)));

        //gsDebugVar(ev.integralInterface(((igrad(u_sol,G.left()) - igrad(u_sol,G.right())) * nv(G).normalized()).sqNorm() * meas(G)));
        //gsDebugVar(ev.integralInterface(((igrad(u_sol.left(),G.left()) - igrad(u_sol.right(),G.right())) * nv(G).normalized()).sqNorm() * meas(G)));
        //gsDebugVar(ev.integralInterface(((igrad(u_sol.left(),G.left()) - igrad(u_sol.right(),G.right())) * nv(G.left()).normalized()).sqNorm() * meas(G)));

        err_time += timer.stop();
        gsInfo<< ". " <<std::flush; // Error computations done
    } //for loop

    //! [Solver loop]

    timer.stop();
    gsInfo<<"\n\nTotal time: "<< setup_time+ma_time+slv_time+err_time <<"\n";
    gsInfo<<"     Setup: "<< setup_time <<"\n";
    gsInfo<<"  Assembly: "<< ma_time    <<"\n";
    gsInfo<<"   Solving: "<< slv_time   <<"\n";
    gsInfo<<"     Norms: "<< err_time   <<"\n";

    //! [Error and convergence rates]
    gsInfo<< "\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";
    gsInfo<< "H2 error: "<<std::scientific<<h2err.transpose()<<"\n";
    gsInfo<< "Deriv Interface error: "<<std::scientific<<IFaceErr.transpose()<<"\n";

    if (!last && numRefine>0)
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

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        ev.options().setSwitch("plot.elements", false);
        ev.options().setInt   ("plot.npts"    , 1000);
        ev.writeParaview( u_sol   , G, "solution");
        //ev.writeParaview( u_ex    , G, "solution_ex");
        //ev.writeParaview( grad(s), G, "solution_grad");
        //ev.writeParaview( grad(f), G, "solution_ex_grad");
        //ev.writeParaview( (f-s), G, "error_pointwise");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    //! [Export data to xml]
    if (output)
    {
        std::string cmdName = "testtest";

        gsMatrix<> error_collection(l2err.rows(), 4);
        error_collection.col(0) = l2err;
        error_collection.col(1) = h1err;
        error_collection.col(2) = h2err;
        error_collection.col(3) = IFaceErr;

        gsFileData<> xml_out;
        xml_out << error_collection;
        // Add solution
        // [...]
        xml_out.save(cmdName);
    }
    //! [Export data to xml]

    return EXIT_SUCCESS;

}// end main
