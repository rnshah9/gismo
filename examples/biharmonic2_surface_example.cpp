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


using namespace gismo;
//! [Include namespace]


int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;

    index_t smoothing = 0;

    index_t numRefine  = 3;
    index_t discreteDegree = 3;
    index_t discreteRegularity = 2;

    index_t gluingDataDegree = -1;
    index_t gluingDataRegularity = -1;

    bool last = false;
    bool info = false;
    bool neumann = false;
    bool cond = false;
    real_t penalty_init = -1.0;
    std::string xml;
    std::string output;

    std::string fn = "surfaces/triangularThreePatchDisturbed.xml";

    index_t geometry = 1000;

    gsCmdLine cmd("Tutorial on solving a Biharmonic problem.");
    cmd.addInt( "s", "smoothing","Smoothing", smoothing );

    cmd.addInt( "p", "discreteDegree","Which discrete degree?", discreteDegree );
    cmd.addInt( "r", "discreteRegularity", "Number of discreteRegularity",  discreteRegularity );
    cmd.addInt( "l", "refinementLoop", "Number of refinementLoop",  numRefine );

    cmd.addInt( "P", "gluingDataDegree","Which degree for gluing data?", gluingDataDegree );
    cmd.addInt( "R", "gluingDataRegularity", "Which regularity for gluing data?",  gluingDataRegularity );

    cmd.addString( "f", "file", "Input geometry file", fn );
    cmd.addInt( "g", "geometry", "Which geometry",  geometry );

    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("info", "Getting the information inside of Approximate C1 basis functions", info);
    cmd.addSwitch("cond","Estimate condition number (slow!)", cond);

    cmd.addSwitch("neumann", "Neumann", neumann);

    cmd.addReal( "y", "penalty", "Fixed Penalty value for Nitsche's method",  penalty_init);

    cmd.addString("x", "xml", "Use the information from the xml file", xml);
    cmd.addString("o", "output", "Output in xml (for python)", output);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Initialize data]
    gsMultiPatch<real_t> mp;
    gsBoundaryConditions<> bc;
    gsFunctionExpr<real_t> f, ms;
    gsOptionList optionList;
    //! [Initialize data]

    //! [Read Argument inputs]
    if (xml.empty()) {
        //! [Read geometry]
        std::string string_geo;
        string_geo = fn;

        gsInfo << "Filedata: " << string_geo << "\n";
        gsReadFile<>(string_geo, mp);
        mp.clearTopology();
        mp.computeTopology();

        gsWriteParaview(mp, "geometry", 2000);

        optionList = cmd;
        gsInfo << "OptionList: " << optionList << "\n";
        gsInfo << "Finished\n";
    }
    //! [Read Argument inputs]


    //! [Read option list]
    discreteDegree = optionList.getInt("discreteDegree");
    discreteRegularity = optionList.getInt("discreteRegularity");
    numRefine = optionList.getInt("refinementLoop");

    gluingDataDegree = optionList.getInt("gluingDataDegree");
    gluingDataRegularity = optionList.getInt("gluingDataRegularity");

    smoothing = optionList.getInt("smoothing");

    penalty_init = optionList.getReal("penalty");

    cond = optionList.getSwitch("cond");
    plot = optionList.getSwitch("plot");
    info = optionList.getSwitch("info");
    //! [Read option list]

    //! [Refinement]
    gsMultiBasis<real_t> dbasis(mp, true);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( discreteDegree); // preserve smoothness
    //dbasis.degreeElevate(discreteDegree- mp.patch(0).degree(0));


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
    {
        dbasis.uniformRefine(1, discreteDegree-discreteRegularity);
    }

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif
    //! [Refinement]

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
    auto u = A.getSpace(bb2);

    // The approx. C1 space
    gsApproxC1Spline<2,real_t> approxC1(mp,dbasis);
    approxC1.options().setSwitch("info",info);
    approxC1.options().setSwitch("plot",plot);
    approxC1.options().setInt("gluingDataDegree",gluingDataDegree);
    approxC1.options().setInt("gluingDataRegularity",gluingDataRegularity);

    // Solution vector and solution variable
    gsMatrix<real_t> solVector;
    auto u_sol = A.getSolution(u, solVector);

    // Recover manufactured solution
    auto u_ex = ev.getVariable(ms, G);
    //! [Problem setup]

    //! [Solver loop]
    gsVector<real_t> l2err(numRefine+1), h1err(numRefine+1), h2err(numRefine+1),
            IFaceErr(numRefine+1), meshsize(numRefine+1), dofs(numRefine+1),
            cond_num(numRefine+1), penalty(numRefine+1);
    gsInfo<< "(dot1=approxC1construction, dot2=assembled, dot3=solved, dot4=got_error)\n"
        "\nDoFs: ";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);
    gsStopwatch timer;
    for (int r=0; r<=numRefine; ++r)
    {
        dbasis.uniformRefine(1,discreteDegree -discreteRegularity);
        meshsize[r] = dbasis.basis(0).getMinCellLength();
        approxC1.update(bb2);

        gsInfo<< "." <<std::flush; // Approx C1 construction done

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

    gsInfo<< "\nMesh-size: " << meshsize.transpose() << "\n";
    gsInfo<< "\nCondition-number: " << cond_num.transpose() << "\n";

    //! [Error and convergence rates]
    gsInfo<< "\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";
    gsInfo<< "H2 error: "<<std::scientific<<h2err.transpose()<<"\n";
    gsInfo<< "Deriv Interface error: "<<std::scientific<<IFaceErr.transpose()<<"\n";

    if (!last && numRefine>0)
    {
        gsInfo<< "EoC (L2): " << std::fixed<<std::setprecision(2)
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

        gsInfo<<   "EoC (Cnum): "<< std::fixed<<std::setprecision(2)
              <<( cond_num.tail(numRefine).array() /
                      cond_num.head(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    }
    //! [Error and convergence rates]

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        ev.options().setSwitch("plot.elements", true);
        ev.options().setInt   ("plot.npts"    , 1000);
        ev.writeParaview( u_sol   , G, "solution");
        //ev.writeParaview( u_ex    , G, "solution_ex");
        //ev.writeParaview( grad(s), G, "solution_grad");
        //ev.writeParaview( grad(f), G, "solution_ex_grad");
        ev.writeParaview( (u_ex-u_sol), G, "error_pointwise");
        gsWriteParaview( mp, "geom",100,true);
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;
}// end main
