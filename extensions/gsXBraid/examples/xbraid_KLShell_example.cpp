/** @file xbraid_example.cpp

    @brief XBraid integration

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, M. Moeller
*/

#include <gismo.h>
#include <gsXBraid/gsXBraid.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>
// #include <gsThinShell/gsArcLengthIterator.h>
#include <gsStructuralAnalysis/gsArcLengthIterator.h>


using namespace gismo;

#ifdef GISMO_WITH_XBRAID

namespace gismo {

/**
   \brief Derived class implementing the XBraid wrapper for the heat equation
*/
template<typename T>
class gsXBraid_app : public gsXBraid< gsMatrix<T> >
{
private:
  // Spatial discretisation parameters
  index_t numRefine, numElevate;

  // Temporal discretisation parameters
  index_t numSteps, typeMethod;
  T tstart, tstop, tstep;

  // Spatial discretizations
  gsMultiPatch<T> m_patches;
  gsMultiBasis<T> m_basis;
  
  // Boundary conditions
  gsBoundaryConditions<T> m_bc;
  gsPointLoads<T> m_pLoads;

  // Assembler options
  gsOptionList Aopt, Topt; //Sopt
  
  // Shell m_assembler
  gsMaterialMatrixBase<T> * m_materialMatrix;
  mutable gsThinShellAssemblerBase<real_t> * m_assembler;
  gsArcLengthIterator<T> * m_arcLength;
  
  // Solution
  gsMatrix<T> sol;

  // Force vector
  gsVector<> m_Force;
  
  // Make Jacobian and residual objects
  typedef typename std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>                                Jacobian_t;
  typedef typename std::function<gsVector<real_t> (gsVector<real_t> const &, real_t, gsVector<real_t> const &) >   ALResidual_t;

  Jacobian_t m_Jacobian;
  ALResidual_t m_ALResidual;

  // Benchmark options
  index_t testCase;

  // Shell parameters
  T thickness;
  T Exx, Eyy, Gxy;
  T PoissonRatio = 0.3;
  T Density    = 1e0;

  // Material options
  bool m_composite;

  // ALM fixed options
  real_t tol        = 1e-6;
  real_t tolU       = 1e-6;
  real_t tolF       = 1e-3;
  index_t maxit     = 20;

  gsFunction<T> * force;
  gsFunction<T> * t;
  gsFunction<T> * E;
  gsFunction<T> * nu;
  gsFunction<T> * rho;

 public:
  /// Contructor
  gsXBraid_app(const gsMpiComm& comm,
               const T&         tstart,
               const T&         tstop,
               index_t          method,
               index_t          numSteps,
               index_t          numRefine,
               index_t          numElevate,
               bool             composite,
               index_t          testCase,
               std::string&     fn)
    : gsXBraid< gsMatrix<T> >::gsXBraid(comm, tstart, tstop, (int)numSteps),
      numRefine(numRefine),
      numElevate(numElevate),
      numSteps(numSteps),
      testCase(testCase),
      m_composite(composite),
      tstart(tstart),
      tstop(tstop),
      tstep( (tstop-tstart)/numSteps )
  {

    // Read m_assembler options
    gsFileData<T> fd(fn);
    if (this->id() == 0) gsInfo << "Loaded file " << fd.lastPath() << "\n";

    fd.getId(0, Aopt); // id=5: m_assembler options
    if (this->id() == 0) gsInfo << "Assembler options:\n" << Aopt << "\n";

    fd.getId(1, Topt); // id=6: multigrid-in-time options
    if (this->id() == 0) gsInfo << "Multigrid-in-time options:\n" << Topt << "\n";

    // Load multipatch
    gsReadFile<>("surface/scordelis_lo_roof_shallow.xml", m_patches);

    // Set material parameters
    if (m_composite)
    {
      Exx = 3300;
      Eyy = 1100;
      Gxy = 660;
      PoissonRatio = 0.25;
    }
    else
    {
      Exx  = 3102.75;
      PoissonRatio = 0.3;
    }

    if (testCase==1)      thickness = 6.35;
    else if (testCase==2) thickness = 12.7;
    else if (testCase==3) thickness = 16.75;

    // Apply refinement and elevation and build basis
    for(index_t i = 0; i< numElevate; ++i)
      m_patches.patch(0).degreeElevate();    // Elevate the degree

    // h-refine
    for(index_t i = 0; i< numRefine; ++i)
      m_patches.patch(0).uniformRefine();

    m_basis = gsMultiBasis<T>(m_patches);

    // Make boundary conditions and loads
    m_pLoads = gsPointLoads<T>();

    gsVector<T> tmp(m_patches.targetDim());
    tmp.setZero();

    GISMO_ASSERT(m_patches.targetDim()==3,"Geometry must be surface (targetDim=3)!");
    // Diaphragm conditions
    m_bc.addCondition(boundary::north, condition_type::dirichlet, nullptr, 0, false, 0 ); // unknown 0 - x
    m_bc.addCondition(boundary::north, condition_type::dirichlet, nullptr, 0, false, 1 ); // unknown 1 - y
    m_bc.addCondition(boundary::north, condition_type::dirichlet, nullptr, 0, false, 2 ); // unknown 2 - z
    // m_bc.addCornerValue(0,boundary::southwest, 0.0, 0, 0); // (corner,value, patch, unknown)
    m_bc.addCondition(boundary::south, condition_type::dirichlet, nullptr, 0, false, 0 ); // unknown 0 - x
    m_bc.addCondition(boundary::south, condition_type::dirichlet, nullptr, 0, false, 1 ); // unknown 1 - y
    m_bc.addCondition(boundary::south, condition_type::dirichlet, nullptr, 0, false, 2 ); // unknown 2 - z

    T Load = -1e1;
    // Point loads
    gsVector<T> point(2);
    gsVector<T> load (3);
    point<< 0.5, 0.5 ;
    load << 0.0, 0.0, Load ;
    m_pLoads.addLoad(point, load, 0 );

    // Output directories
    /*
      std::string output = "solution";
      std::string dirname = "ArcLengthResults";

      dirname = dirname + "/" +  "Roof_t="+ std::to_string(thickness) + "-r=" + std::to_string(numRefine) + "-e" + std::to_string(numElevate) +"_solution";
      output =  "solution";
      wn = output + "data.txt";

      std::string commands = "mkdir -p " + dirname;
      const char *command = commands.c_str();
      system(command);

      // plot geometry
      if (plot)
        gsWriteParaview(mp,dirname + "/" + "mp",1000,mesh);

      if (write)
        initStepOutput(dirname + "/" + wn, writePoints);

    */

    // Define material object
    gsMaterialMatrixBase<real_t>* m_materialMatrix;



    std::vector<gsFunctionSet<> * > Gs;
    std::vector<gsFunctionSet<> * > Ts;
    std::vector<gsFunctionSet<> * > Phis;

    if (m_composite)
    {
      real_t pi = math::atan(1)*4;
      index_t kmax = 3;
      Gs.resize(kmax);
      Ts.resize(kmax);
      Phis.resize(kmax);

      gsMatrix<> Gmat = gsCompositeMatrix(Exx,Eyy,Gxy,PoissonRatio,PoissonRatio*Eyy/Exx);
      Gmat.resize(Gmat.rows()*Gmat.cols(),1);
      gsConstantFunction<> Gfun(Gmat,3);
      Gs[0] = Gs[1] = Gs[2] = &Gfun;

      gsConstantFunction<> phi1,phi2,phi3;
      phi1.setValue(pi/2,3);
      phi2.setValue(0,3);
      phi3.setValue(pi/2,3);

      Phis[0] = &phi1;
      Phis[1] = &phi2;
      Phis[2] = &phi3;

      gsConstantFunction<> thicks(thickness/kmax,3);
      Ts[0] = Ts[1] = Ts[2] = &thicks;

    }


    force = new gsConstantFunction<>(tmp,3);
    t = new gsFunctionExpr<>(std::to_string(thickness), 3);
    E = new gsFunctionExpr<>(std::to_string(Exx),3);
    nu = new gsFunctionExpr<>(std::to_string(PoissonRatio),3);
    rho = new gsFunctionExpr<>(std::to_string(Density),3);

    std::vector<gsFunction<>*> parameters;
    parameters.resize(2);
    parameters[0] = E;
    parameters[1] = nu;

    gsOptionList options;
    if (m_composite)
    {
        m_materialMatrix = new gsMaterialMatrixComposite<3,real_t>(m_patches,Ts,Gs,Phis);
    }
    else
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
        m_materialMatrix = getMaterialMatrix<3,real_t>(m_patches,*t,parameters,*rho,options);
    }

    // Make shell m_assembler
    m_assembler = new gsThinShellAssembler<3, real_t, true >(m_patches,m_basis,m_bc,*force,m_materialMatrix);
    // Construct m_assembler object
    m_assembler->setOptions(Aopt);
    m_assembler->setPointLoads(m_pLoads);

    // Function for the Jacobian
    m_Jacobian = [this](gsVector<real_t> const &x)
    {
      gsMultiPatch<> mp_def;
      m_assembler->constructSolution(x,mp_def);
      m_assembler->assembleMatrix(mp_def);
      gsSparseMatrix<real_t> m = m_assembler->matrix();
      return m;
    };
    // Function for the Residual
    m_ALResidual = [this](gsVector<real_t> const &x, real_t lam, gsVector<real_t> const &force)
    {
      gsMultiPatch<> mp_def;
      m_assembler->constructSolution(x,mp_def);
      m_assembler->assembleVector(mp_def);
      gsVector<real_t> Fint = -(m_assembler->rhs() - force);
      gsVector<real_t> result = Fint - lam * force;
      return result; // - lam * force;
    };

    // Assemble linear system to obtain the force vector
    m_assembler->assemble();
    m_Force = m_assembler->rhs();

    // Define ALM object
    m_arcLength = new gsArcLengthIterator<real_t>(m_Jacobian, m_ALResidual, m_Force);

    m_arcLength->options().setInt("Solver",0); // LDLT solver
    m_arcLength->options().setInt("BifurcationMethod",0); // 0: determinant, 1: eigenvalue
    m_arcLength->options().setInt("Method",method);
    m_arcLength->options().setInt("AngleMethod",0); // 0: step, 1: iteration
    // m_arcLength->options().setReal("Scaling",0.0);
    m_arcLength->options().setReal("Tol",tol);
    m_arcLength->options().setReal("TolU",tolU);
    m_arcLength->options().setReal("TolF",tolF);
    m_arcLength->options().setInt("MaxIter",maxit);
    m_arcLength->options().setSwitch("Verbose",true);
    // m_arcLength->options().setReal("Relaxation",relax);
    // if (quasiNewtonInt>0)
    // {
    //   quasiNewton = true;
    //   m_arcLength->options().setInt("QuasiIterations",quasiNewtonInt);
    // }
    // m_arcLength->options().setSwitch("Quasi",quasiNewton);

    if (this->id() == 0) gsInfo << "ALM options:\n" << m_arcLength->options() << "\n";

    m_arcLength->applyOptions();
    m_arcLength->initialize();

    /*
     gsParaviewCollection collection(dirname + "/" + output);
     gsParaviewCollection Smembrane(dirname + "/" + "membrane");
     gsParaviewCollection Sflexural(dirname + "/" + "flexural");
     gsParaviewCollection Smembrane_p(dirname + "/" + "membrane_p");
     */
    gsMultiPatch<> deformation = m_patches;

    // Set XBraid options
    this->SetCFactor(Topt.getInt("CFactor"));
    this->SetMaxIter(Topt.getInt("maxIter"));
    this->SetMaxLevels(Topt.getInt("maxLevel"));
    this->SetMaxRefinements(Topt.getInt("numMaxRef"));
    this->SetMinCoarse(Topt.getInt("minCLevel"));
    this->SetNFMG(Topt.getInt("numFMG"));
    this->SetNFMGVcyc(Topt.getInt("numFMGVcyc"));
    this->SetNRelax(Topt.getInt("numRelax"));
    this->SetAccessLevel(Topt.getInt("access"));
    this->SetPrintLevel(Topt.getInt("print"));
    this->SetStorage(Topt.getInt("numStorage"));
    this->SetTemporalNorm(Topt.getInt("norm"));

    if (Topt.getSwitch("fmg"))           this->SetFMG();
    if (Topt.getSwitch("incrMaxLevels")) this->SetIncrMaxLevels();
    if (Topt.getSwitch("periodic"))      this->SetPeriodic(1); else this->SetPeriodic(0);
    if (Topt.getSwitch("refine"))        this->SetRefine(1);   else this->SetRefine(0);
    if (Topt.getSwitch("sequential"))    this->SetSeqSoln(1);  else this->SetSeqSoln(0);
    if (Topt.getSwitch("skip"))          this->SetSkip(1);     else this->SetSkip(0);
    if (Topt.getSwitch("spatial"))       this->SetSpatialCoarsenAndRefine();
    if (Topt.getSwitch("tol"))           this->SetAbsTol(Topt.getReal("absTol"));
    else                                 this->SetRelTol(Topt.getReal("relTol"));

    real_t dL = tstep;

    if (this->id() == 0)
    {

      gsInfo<<"ID = "<<this->id()<<"\n";

      m_arcLength->setLength(dL);

      // Make objects for previous solutions
      real_t Lold = 0;
      gsMatrix<> Uold = m_Force;
      Uold.setZero();

      gsMatrix<> solVector;
      real_t indicator = 0.0;
      m_arcLength->setIndicator(indicator); // RESET INDICATOR
      bool bisected = false;
      real_t dL0 = dL;
      for (index_t k=0; k<numSteps; k++)
      {
        gsInfo<<"Load step "<< k<<"\n";
        gsInfo<<"dL = "<<dL<<"\n";
        // assembler->constructSolution(solVector,solution);
        m_arcLength->step();

        // gsInfo<<"m_U = "<<arcLength.solutionU()<<"\n";
        if (!(m_arcLength->converged()))
        {
          gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
          dL = dL / 2.;
          m_arcLength->setLength(dL);
          m_arcLength->setSolution(Uold,Lold);
          bisected = true;
          k -= 1;
          continue;
        }

        // if (SingularPoint)
        // {
        //   arcLength.computeStability(arcLength.solutionU(),quasiNewton);
        //   if (arcLength.stabilityChange())
        //   {
        //     gsInfo<<"Bifurcation spotted!"<<"\n";
        //     arcLength.computeSingularPoint(1e-4, 5, Uold, Lold, 1e-10, 0, false);
        //     arcLength.switchBranch();
        //     dL0 = dL = dL;
        //     arcLength.setLength(dL);
        //   }
        // }

        indicator = m_arcLength->indicator();

        solVector = m_arcLength->solutionU();
        Uold = solVector;
        Lold = m_arcLength->solutionL();

        if (!bisected)
        {
          dL = dL0;
          m_arcLength->setLength(dL);
        }
        bisected = false;

      }

    }


    /////////////////////////////////////////////////////////////////////////////////////////////
    //                           Code for heat equation starts here                            //
    /////////////////////////////////////////////////////////////////////////////////////////////
    

    // if (this->id() == 0) {
 
    //   gsStopwatch clock;
    //   clock.restart();
      
    //   sol.setZero(M.numDofs());
 
    //   switch((gsXBraid_typeMethod)typeMethod) {
    //   case gsXBraid_typeMethod::FE_FE:
    //   case gsXBraid_typeMethod::FE_BE:
    //     // Forward Euler method
        
    //     for ( int i = 1; i<=numSteps; ++i) // for all timesteps
    //       // Compute the system for the timestep i (rhs is assumed constant wrt time)
    //       sol = m_solver[0]->solveWithGuess(tstep*K.rhs() +
    //                                                (M.matrix()-tstep*K.matrix())*sol,
    //                                                sol);
    //     break;
        
    //   case gsXBraid_typeMethod::BE_BE:
    //     // Backward Euler method
        
    //    for ( int i = 1; i<=numSteps; ++i) // for all timesteps
    //       // Compute the system for the timestep i (rhs is assumed constant wrt time)
    //       sol = m_solver[0]->solveWithGuess(tstep*K.rhs() + (M.matrix())*sol, sol);
    //    break;
        
    //   case gsXBraid_typeMethod::CN_CN:
    //   case gsXBraid_typeMethod::CN_BE:
    //     // Crank-Nicholson method
    //     for ( int i = 1; i<=numSteps; ++i) // for all timesteps
    //       // Compute the system for the timestep i (rhs is assumed constant wrt time)
    //       sol = m_solver[0]->solveWithGuess(tstep*K.rhs() +
    //                                                (M.matrix()-tstep*0.5*K.matrix())*sol,
    //                                                sol);
    //     break;
        
    //   default:
    //     throw std::runtime_error("Unsupported time-stepping method");
    //   }
      
    //   gsInfo << "wall time = " << clock.stop() << "\n"
    //          << "L2 norm of the solution  = " << sol.norm() << "\n";
      
    //   // gsExprEvaluator<T> ev(M);
    //   // solution u_sol = M.getSolution(u_M, sol);
    //   // variable u_ex  = ev.getVariable(ms, G_M);
    //   // T l2err = math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G_M) ) );
    //   // T h1err = l2err +
    //   //   math::sqrt(ev.integral( ( igrad(u_ex) - grad(u_sol)*jac(G_M).inv() ).sqNorm() * meas(G_M) ));

    //   // gsInfo << "L2 error of the solution = " << l2err << "\n"
    //   //        << "H1 error of the solution = " << h1err << std::flush;
    // }
  }

  /// Destructor
  virtual ~gsXBraid_app()
  {

  }
  
  /// Creates instance from command line argument
  static inline gsXBraid_app create(const gsMpiComm& comm,
                                    int              argc,
                                    char**           argv)
  {
    // Input options
    int numElevate    = 1;
    int numRefine       = 1;
    bool plot         = false;
    bool mesh         = false;
    bool stress       = false;
    bool membrane     = false;
    bool quasiNewton  = false;
    int quasiNewtonInt= -1;
    bool adaptive     = false;
    int step          = 10;
    int method        = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method; 4: extended iterations)
    bool deformed     = false;

    bool composite = false;

    real_t relax      = 1.0;

    int testCase      = 1;

    int result        = 0;

    bool write        = false;

    // Arc length method options
    real_t dL        = 0.5; // Ard length to find bifurcation


    std::string wn("data.csv");

    std::string opts("options/MGRITsolver_options.xml");

    index_t numSteps      = 10;
    T       tfinal        = 0.1;
        
    gsCmdLine cmd("Arc-length analysis for thin shells.");

    cmd.addString( "f", "file", "Input XML file for m_assembler options", opts );

    cmd.addInt("t", "testcase", "Test case: 0: clamped-clamped, 1: pinned-pinned, 2: clamped-free", testCase);

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numRefine);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
    cmd.addSwitch("composite", "Composite material", composite);

    cmd.addInt("m","Method", "Arc length method; 1: Crisfield's method; 2: RIks' method.", method);
    cmd.addReal("L","dL", "arc length", dL);
    // cmd.addReal("A","relaxation", "Relaxation factor for arc length method", relax);

    // cmd.addInt("q","QuasiNewtonInt","Use the Quasi Newton method every INT iterations",quasiNewtonInt);
    cmd.addInt("N", "maxsteps", "Maximum number of steps", step);

    cmd.addSwitch("adaptive", "Adaptive length ", adaptive);
    // cmd.addSwitch("quasi", "Use the Quasi Newton method", quasiNewton);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("mesh", "Plot mesh?", mesh);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
    cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);
    cmd.addSwitch("deformed", "plot on deformed shape", deformed);
    cmd.addSwitch("write", "write to file", write);
    
    cmd.getValues(argc,argv);

    // Create instance
    gsXBraid_app<T> app(comm, 0.0, tfinal, method, numSteps, numRefine, numElevate, composite,testCase, opts);
    
    return app;
  }

  /// Initializes a vector
  braid_Int Init(braid_Real    t,
                 braid_Vector *u_ptr)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
  {
    gsMatrix<T>* u = new gsMatrix<T>(m_Force.rows() + 1, 1);
    
    if (t != tstart) {
      // Intermediate solution
      u->setZero();
    } else {
      // Initial solution
      u->setZero();
    }

    *u_ptr = (braid_Vector) u;
    return braid_Int(0);
  }
  
  /// Performs a single step of the parallel-in-time multigrid
  braid_Int Step(braid_Vector    u,                   //// current solution
                 braid_Vector    ustop,               //// initial guess
                 braid_Vector    fstop,               //// forcing
                 BraidStepStatus &status)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
  {
    gsMatrix<T>* u_ptr = (gsMatrix<T>*) u;
    gsMatrix<T>* ustop_ptr = (gsMatrix<T>*) ustop;

    // XBraid forcing
    if (fstop != NULL) {
      gsMatrix<T>* fstop_ptr = (gsMatrix<T>*) fstop;
      *u_ptr += *fstop_ptr;
    }
    
    // Get load step information
    std::pair<braid_Real, braid_Real> load =
      static_cast<gsXBraidStepStatus&>(status).timeInterval();
    T dL(load.second - load.first);

    gsMatrix<T> Uold = u_ptr->block(0,0,m_Force.rows(),1);
    T Lold = u_ptr->at(m_Force.rows());

    m_arcLength->setLength(dL);
    m_arcLength->setSolution(Uold,Lold);
    m_arcLength->resetStep();

    m_arcLength->step();

    // IF WE WANT TO DO SOMETHING SPECIAL FOR LEVEL 0
    // if (static_cast<gsXBraidStepStatus&>(status).level() == 0) {
    // } else {
    // }

      
    // Carry out adaptive refinement in time
    if (static_cast<gsXBraidStepStatus&>(status).level() == 0) {
      braid_Real error = static_cast<gsXBraidStepStatus&>(status).error();
      if (error != braid_Real(-1.0)) {
        braid_Int rfactor = (braid_Int) std::ceil( std::sqrt( error / 1e-3) );
        status.SetRFactor(rfactor);
      } else
        status.SetRFactor(1);
    }
    
    return braid_Int(0);
  }

  /*
      Little bit CHANGED
   */
  /// Sets the size of the MPI communication buffer
  braid_Int BufSize(braid_Int         *size_ptr,
                    BraidBufferStatus &status)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
  {
    *size_ptr = sizeof(T)*(m_Force.size()+2);
    return braid_Int(0);
  }

  /*
      NOT CHANGED

      Writes to disk --> later
   */
  /// Handles access for input/output
  braid_Int Access(braid_Vector       u,
                   BraidAccessStatus &status)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
  {
    if (static_cast<gsXBraidAccessStatus&>(status).done() &&
        static_cast<gsXBraidAccessStatus&>(status).timeIndex() ==
        static_cast<gsXBraidAccessStatus&>(status).times()) {
      gsMatrix<T>* u_ptr = (gsMatrix<T>*) u;
      gsInfo << "norm of the solution = " << u_ptr->norm() << std::endl;    
    }
    return braid_Int(0);
  }

  /*
      NOT CHANGED
   */
  /// Performs spatial coarsening
  braid_Int Coarsen(braid_Vector           fu,
                    braid_Vector          *cu_ptr,
                    BraidCoarsenRefStatus &status)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
  {
    // gsInfo << "Coarsen on level = "
    //        << static_cast<gsXBraidCoarsenRefStatus&>(status).level()
    //        << " of "
    //        << static_cast<gsXBraidCoarsenRefStatus&>(status).levels()
    //        << "\n";
    gsMatrix<T> *fu_ptr = (gsMatrix<T>*) fu;    
    gsMatrix<T>* cu     = new gsMatrix<T>();
    *cu = *fu_ptr;
    *cu_ptr = (braid_Vector) cu;
    return braid_Int(0);
  }
  
  /*
      NOT CHANGED
   */
  // Performs spatial refinement
  braid_Int Refine(braid_Vector           cu,
                   braid_Vector          *fu_ptr,
                   BraidCoarsenRefStatus &status)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
  {
    // gsInfo << "Refine on level = "
    //        << static_cast<gsXBraidCoarsenRefStatus&>(status).level()
    //        << " of "
    //        << static_cast<gsXBraidCoarsenRefStatus&>(status).levels()
    //        << "\n";
    gsMatrix<T> *cu_ptr = (gsMatrix<T>*) cu;    
    gsMatrix<T>* fu     = new gsMatrix<T>();
    *fu = *cu_ptr;
    *fu_ptr = (braid_Vector) fu;
    return braid_Int(0);
  }
};
  
} // ending namespace gismo

#endif

int main(int argc, char**argv)
{
#ifdef GISMO_WITH_XBRAID
  
  // Initialize the MPI environment and obtain the world communicator
  gsMpiComm comm = gsMpi::init(argc, argv).worldComm();

  // Set up app structure
  gsXBraid_app<real_t> app = gsXBraid_app<real_t>::create(comm, argc, argv);

  // Perform parallel-in-time multigrid
  app.solve();

#else

  gsInfo << "\n";
 
#endif

  return 0;
  
}


template <class T>
void initStepOutput(const std::string name, const gsMatrix<T> & points)
{
  std::ofstream file;
  file.open(name,std::ofstream::out);
  file  << std::setprecision(20)
        << "Deformation norm" << ",";
        for (index_t k=0; k!=points.cols(); k++)
        {
          file<< "point "<<k<<" - x" << ","
              << "point "<<k<<" - y" << ","
              << "point "<<k<<" - z" << ",";
        }

  file  << "Lambda" << ","
        << "Indicator"
        << "\n";
  file.close();

  gsInfo<<"Step results will be written in file: "<<name<<"\n";
}

template <class T>
void writeStepOutput(const gsArcLengthIterator<T> & arcLength, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const index_t extreme, const index_t kmax) // extreme: the column of point indices to compute the extreme over (default -1)
{
  gsMatrix<T> P(2,1), Q(2,1);
  gsMatrix<T> out(3,points.cols());
  gsMatrix<T> tmp;

  for (index_t p=0; p!=points.cols(); p++)
  {
    P<<points.col(p);
    deformation.patch(0).eval_into(P,tmp);
    out.col(p) = tmp;
  }

  std::ofstream file;
  file.open(name,std::ofstream::out | std::ofstream::app);
  if (extreme==-1)
  {
    file  << std::setprecision(6)
          << arcLength.solutionU().norm() << ",";
          for (index_t p=0; p!=points.cols(); p++)
          {
            file<< out(0,p) << ","
                << out(1,p) << ","
                << out(2,p) << ",";
          }

    file  << arcLength.solutionL() << ","
          << arcLength.indicator() << ","
          << "\n";
  }
  else if (extreme==0 || extreme==1)
  {
    gsMatrix<T> out2(kmax,points.cols()); // evaluation points in the rows, output (per coordinate) in columns
    for (int p = 0; p != points.cols(); p ++)
    {
      Q.at(1-extreme) = points(1-extreme,p);
      for (int k = 0; k != kmax; k ++)
      {
        Q.at(extreme) = 1.0*k/(kmax-1);
        deformation.patch(0).eval_into(Q,tmp);
        out2(k,p) = tmp.at(2); // z coordinate
      }
    }

    file  << std::setprecision(6)
          << arcLength.solutionU().norm() << ",";
          for (index_t p=0; p!=points.cols(); p++)
          {
            file<< out(0,p) << ","
                << out(1,p) << ","
                << std::max(abs(out2.col(p).maxCoeff()),abs(out2.col(p).minCoeff())) << ",";
          }

    file  << arcLength.solutionL() << ","
          << arcLength.indicator() << ","
          << "\n";
  }
  else
    GISMO_ERROR("Extremes setting unknown");

  file.close();
}
