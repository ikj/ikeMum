#ifndef UGPWITHCVCOMPFLOW_H
#define UGPWITHCVCOMPFLOW_H

#include <iomanip>      // std::setprecision

#include "MiscUtils.h"
#include "tc_vec3d.h"

#include "UgpWithCv2.h"
#include "CdpFilter.h"
#include "MshFilter.h"
#include "Param.h"

#include "AtmosphericCond.h"

#include "HypreSolver.h"
#include "PetscSolver.h"

#include "Logging.h"
using namespace logging;

#define MAX_WARNING_MESSAGES_UGP 50

//#define USE_ARTIF_VISC

//#define temp_reconstruction   // determine whether temperature or pressure reconstruction is used at faces
//#define alpha_limiter         // determine if rho*Phi is limited based on rho or rho*Phi


/**
 * Class contains the main routines for unstructured Euler for now
 *
 * \author Frank Ham, Rene Pecnik and Vincent Terrapon
 * \date March, 2008
 */


class UgpWithCvCompFlow : public UgpWithCv2Op, public ParamMap
{
public:   // constructors/destructors

  /**
   * standard constructor needed for Models that virtual inherits from UgpWithCvCompFlow
   */
  UgpWithCvCompFlow() {    init();  }

  /**
   * constructor, pass ParamMap
   */
  UgpWithCvCompFlow(ParamMap &p) : ParamMap(p) {    init();  }

  /**
   * constructor, pass name of joe's input file
   */
  UgpWithCvCompFlow(char *name) : ParamMap(name) {    init();  }


  /**
   * registers scalars and vectors for main equations (Navier-Stokes)
   * and set main parameters (solver: nsteps; NS-eq: gravity...)
   */
  void init()
  {
	// For IKE
	firstCall_JOEcalcViscousFluxScalar = true;

    Param *param;
    if (getParam(param, "LOGGING_LEVEL"))
    {
      string name = param->getString();
      if      (name == "EVERYTHING") setLoggingLevel(EVERYTHING);
      else if (name == "DEBUG_HI"  ) setLoggingLevel(DEBUG_HI  );
      else if (name == "DEBUG_LO"  ) setLoggingLevel(DEBUG_LO  );
      else if (name == "INFO_HI"   ) setLoggingLevel(INFO_HI   );
      else if (name == "INFO_LO"   ) setLoggingLevel(INFO_LO   );
      else if (name == "WARNING"   ) setLoggingLevel(WARNING   );
      else if (name == "ERROR"     ) setLoggingLevel(ERROR     );
      else if (name == "CRITICAL"  ) setLoggingLevel(CRITICAL  );
      else if (name == "FATAL"     ) setLoggingLevel(FATAL     );
      else if (name == "SILENT"    ) setLoggingLevel(SILENT    );
      else {
        lout(WARNING) << "Warning: Unknown logging level, using default.\n";
      }
    }


    if (mpi_rank == 0)    cout << "UgpWithCvCompFlow()"<< endl;

    // ----------------------------------------------------------------------------------------
    // register data...
    // ----------------------------------------------------------------------------------------

    // conservative variables
    rho   = NULL;         registerScalar(rho,  "RHO",  CV_DATA);
    rhou  = NULL;         registerVector(rhou, "RHOU", CV_DATA);
    rhoE  = NULL;         registerScalar(rhoE, "RHOE", CV_DATA);

    // primitive variables
    vel      = NULL;      registerVector(vel,      "vel",      CV_DATA);
    press    = NULL;      registerScalar(press,    "press",    CV_DATA);
    temp     = NULL;      registerScalar(temp,     "temp",     CV_DATA);
    enthalpy = NULL;      registerScalar(enthalpy, "enthalpy", CV_DATA);

    // material properties
    gamma    = NULL;      registerScalar(gamma,    "gamma",    CV_DATA);
    RoM      = NULL;      registerScalar(RoM,      "RoM",      CV_DATA);
    sos      = NULL;      registerScalar(sos,      "sos",      CV_DATA);
    
    // diffusion coefficient for momentum and energy at cell faces
    mul_fa = NULL;
    lamOcp_fa = NULL;
    
    // turbulence 
    mut_fa = NULL;          // defined at faces, allocated in initializeFromRestartFile
    kine   = NULL;          // registerScalar(kine, "kine", CV_DATA);

    // just to take a look at the residual field of energy
    residField  = NULL;     registerScalar(residField,  "residField",  CV_DATA);  // Residual of RhoE
    residField2 = NULL;     registerScalar(residField2, "residField2", CV_DATA);  // Residual of omega, nuSA, or eps
    residField3 = NULL;     registerScalar(residField3, "residField3", CV_DATA);  // Residual of kine
    residField4 = NULL;     registerScalar(residField4, "residField4", CV_DATA);  // Residual of ZMean
    residField5 = NULL;     registerScalar(residField5, "residField5", CV_DATA);  // Residual of ZVar
    residField6 = NULL;     registerScalar(residField6, "residField6", CV_DATA);  // Residual of CMean

    local_dt = NULL;        registerScalar(local_dt, "local_dt", CV_DATA);


    grad_rho = NULL;        // register only for second order euler flux or if viscous
    grad_u = NULL;          // register only for second order euler flux or if viscous
    grad_p = NULL;          // register only for second order euler flux, if pressure reconstruction is used
    grad_temp = NULL;       // register only for second order euler flux, if temperature reconstruction is used
    grad_enthalpy = NULL;   // register only if viscous
    
    strMag = NULL;
    vortMag = NULL;
    diverg = NULL;


    // ----------------------------------------------------------------------------------------
    // register flow parameters
    // ----------------------------------------------------------------------------------------

    gravity[0] = getDoubleParam("GX", "0.0");
    gravity[1] = getDoubleParam("GY", "0.0");
    gravity[2] = getDoubleParam("GZ", "0.0");


    turbModel = NONE;

    viscMode      = getStringParam("MU_MODE", "POWERLAW");
    mu_ref        = getDoubleParam("MU_REF", "0.0");           // NOTE: if MU_MODE == SUTHERLAND then set MU_REF = 1.716e-5

//IKJ
    d_scalar_comb_ref = getDoubleParam("D_SCALAR_COMB_REF", "0.0");

    mu_power_law  = getDoubleParam("MU_POWER_LAW", "0.76");    // note: specify 0 here for constant visc = mu_ref
    SL_Tref       = getDoubleParam("SL_TREF", "273.15");       // Sutherland law, Tref should be 273.15
    SL_Sref       = getDoubleParam("SL_SREF", "110.4");        // Sutherland law, S value should be 110.4

    // if ALTITUDE is specified in the input file T_ref, p_ref, rho_ref and mu_ref are recalculated using getAthmospericConditions
    if (checkParam("ALTITUDE"))
    {
      if (!getAthmospericConditions(T_ref, p_ref, rho_ref, getDoubleParam("ALTITUDE"), 9.81, 1.4, 287.0))
      {
        cerr << "ERROR in getAthmospericConditions occurred" << endl;
        throw(-1);
      }

      viscMode = "SUTHERLAND";    // change to Sutherland law for viscosity if height is specified 
      mu_ref = 1.716e-5;          // hard coded for reference viscosity at conditions given below
      SL_Tref = 273.15;           // hard coded for reference temperature in Sutherland law
      SL_Sref = 110.4;            // hard coded for constant in Sutherland law
    }
    else
    {
      T_ref    = getDoubleParam("T_REF");
      p_ref    = getDoubleParam("P_REF");
      rho_ref  = getDoubleParam("RHO_REF");
    }

    // ideal gas law is P / rho = R * T ...
    GAMMA  = getDoubleParam("GAMMA", "1.4");
    R_gas  = p_ref/rho_ref/T_ref;
    Pr     = getDoubleParam("PR", "0.72");
    PrTurb = getDoubleParam("PRTURB", "0.9");

    // ----------------------------------------------------------------------------------------
    // register solver parameters
    // ----------------------------------------------------------------------------------------

    // check if second order for Navier-Stokes

    if (checkParam("SPATIAL_SECOND_ORDER"))       sndOrder = true;
    else                                          sndOrder = false;

     if (checkParam("SCALAR_FIRST_ORDER"))        stOrderScalar = true;
    else                                          stOrderScalar = false;

    
    // register type of gradient reconstruction and limiters

    if (!checkParam("GRAD_RECONSTRUCTION"))
    {
      ParamMap::add("GRAD_RECONSTRUCTION  SCHEME=GREENGAUSS  LIMITER=YES  EPS=0.001");    // add default values
      if (mpi_rank == 0)
        cout << "WARNING: added keyword \"GRAD_RECONSTRUCTION  SCHEME=GREENGAUSS  LIMITER=YES  EPS=0.001\"" <<
                " to parameter map!" << endl;
    }
    string gradReconstr = getParam("GRAD_RECONSTRUCTION")->getString("SCHEME");

    int GRAD_GG = 1;
    int GRAD_LS = 2;
    int GRAD_SDWLS = 3;

    if      (gradReconstr == "GREENGAUSS")       gradreconstruction = GRAD_GG;
    else if (gradReconstr == "LEASTSQUARES")     gradreconstruction = GRAD_LS;
    else if (gradReconstr == "SDWLS")            gradreconstruction = GRAD_SDWLS;
    else
    {
      if (mpi_rank == 0)
        cerr << "ERROR: gradient reconstruction type "
             << "\"" << gradReconstr << "\"" << " not recognized, choose GRAD_RECONSTRUCTION = \"GREENGAUSS\", \"LEASTSQUARE\" or \"SDWLS\"" << endl;
      throw(-1);
    }

    string limiterstr;
    if ((gradreconstruction == GRAD_GG) || (gradreconstruction == GRAD_LS))
    {
      limiterstr   = getParam("GRAD_RECONSTRUCTION")->getString("LIMITER");

      if      ((limiterstr == "NOLIMITER") || (limiterstr == "NO"))              limiterNavierS = NOLIMITER;
      else if ((limiterstr == "BARTH_JESPERSEN_MOD") || (limiterstr == "YES"))   limiterNavierS = BARTH_JESPERSEN_MOD;
      else 
      {
        if (mpi_rank == 0)
          cerr << "ERROR: Grandient limiter (" << limiterstr <<
              ") not recognized, choose GRAD_RECONSTRUCTION = SCHEME=<YES, NO, NOLIMITER or BARTH_JESPERSEN_MOD>" << endl;
        throw(-1);
      }

      epsilonSDWLS = getParam("GRAD_RECONSTRUCTION")->getDouble("EPS");
    }
    else if (gradreconstruction == GRAD_SDWLS)
    {
      limiterstr = "VANLEER";
      limiterNavierS = 0;
      epsilonSDWLS = getParam("GRAD_RECONSTRUCTION")->getDouble("EPS");
    }

    // register CFL condition
    
    if (checkParam("CFL"))
    {
      timeStepMode = "CFL";
      cfl = getDoubleParam("CFL");
    }
    else if (checkParam("CFL_MIN"))
    {
      timeStepMode = "CFL_MIN";
      cfl = getDoubleParam("CFL_MIN");
    }
    else if (checkParam("DT"))
    {
      timeStepMode = "DT";
      cfl = 0.0;
    }
    else
    {
      cout << "ERROR: DT, CFL or CFL_MIN not specified" << endl;
      throw(-1);
    }


    nsteps = getIntParam("NSTEPS", "100");
    registerValue(step,"STEP");
    registerValue(time,"TIME");



    write_restart = getIntParam("WRITE_RESTART", "-1");     // default is no writing
    check_interval = getIntParam("CHECK_INTERVAL", "10");


    // ----------------------------------------------------------------------------------------
    // set Navier-Stokes solver
    // ----------------------------------------------------------------------------------------

    Param * p;
    string tIntName = getStringParam("TIME_INTEGRATION", "BACKWARD_EULER");




    // ----------------------------------------------------------------------------------------
    // pick linear solver for Navier-Stokes
    // ----------------------------------------------------------------------------------------

    petscSolver = NULL;
    nbocv_v_global = NULL;

    if (getParam(p,"LINEAR_SOLVER_NS"))
    {
      string name = p->getString();

      if (name == "PETSC_GMRES")
      {
        if (mpi_rank == 0)        cout << "LINEAR_SOLVER_NS: PETSC_GMRES" << endl;
        linearSolverNS = PETSC_GMRES;
      }
      else if (name == "BCGSTAB")
      {
        if (mpi_rank == 0)        cout << "LINEAR_SOLVER_NS: BCGSTAB" << endl;
        linearSolverNS = BCGSTAB;
      }
      else if (name == "PETSC_LU_MUMPS")
      {
        if (mpi_rank == 0)        cout << "LINEAR_SOLVER_NS: PETSC_LU_MUMPS" << endl;
        linearSolverNS = PETSC_LU_MUMPS;
      }
      else
      {
        if (mpi_rank == 0)        cerr << "Error: unrecognized LINEAR_SOLVER_NS: " << name << endl;
        throw(-1);
      }
    }
    else
    {
      if (mpi_rank == 0)          cout << "LINEAR_SOLVER_NS: BCGSTAB (default)" << endl;
      linearSolverNS = BCGSTAB;
    }


    // ----------------------------------------------------------------------------------------
    // pick linear solver for scalar
    // ----------------------------------------------------------------------------------------
    petscSolverScalars = NULL;

    if (getParam(p,"LINEAR_SOLVER_SCALARS"))
    {
      string name = p->getString();

      if (name == "PETSC_GMRES")
      {
        if (mpi_rank == 0)        cout << "LINEAR_SOLVER_SCALARS: PETSC_GMRES" << endl;
        linearSolverScalars = PETSC_GMRES;
      }
      else if (name == "BCGSTAB")
      {
        if (mpi_rank == 0)        cout << "LINEAR_SOLVER_SCALARS: BCGSTAB" << endl;
        linearSolverScalars = BCGSTAB;
      }
      else if (name == "PETSC_LU_MUMPS")
      {
        if (mpi_rank == 0)        cout << "LINEAR_SOLVER_SCALARS: PETSC_LU_MUMPS" << endl;
        linearSolverScalars = PETSC_LU_MUMPS;
      }
      else
      {
        if (mpi_rank == 0)        cerr << "Error: unrecognized LINEAR_SOLVER_NS: " << name << endl;
        throw(-1);
      }
    }
    else
    {
      if (linearSolverNS == PETSC_GMRES)
      {
        linearSolverScalars = PETSC_GMRES;
        if (mpi_rank == 0)          cout << "LINEAR_SOLVER_SCALARS: PETSC_GMRES (default)" << endl;
      }
      if (linearSolverNS == PETSC_LU_MUMPS)
      {
    	  linearSolverScalars = PETSC_LU_MUMPS;
    	  if (mpi_rank == 0)          cout << "LINEAR_SOLVER_SCALARS: PETSC_LU_MUMPS (default)" << endl;
      }
      else
      {
        linearSolverScalars = BCGSTAB;
        if (mpi_rank == 0)          cout << "LINEAR_SOLVER_SCALARS: BCGSTAB (default)" << endl;
      }
    }

    initial_flowfield_output = getStringParam("INITIAL_FLOWFIELD_OUTPUT", "YES");

    // ----------------------------------------------------------------------------------------
    // write some parameters on the screen
    // ----------------------------------------------------------------------------------------

    if (mpi_rank == 0)
    {
    	cout << "Gas properties         " << endl;
    	cout << "    GAMMA            : " << GAMMA << endl;
    	cout << "    R_GAS            : " << R_gas << endl;
    	cout << "    P_REF            : " << p_ref << endl;
    	cout << "    RHO_REF          : " << rho_ref << endl;
    	cout << "    T_REF            : " << T_ref << endl;
    	cout << "    SOS_REF          : " << sqrt(GAMMA*R_gas*T_ref) << endl;
    	cout << "Material properties    " << endl;
    	cout << "    MU_MODE          : " << viscMode << endl;
    	cout << "      mu_ref         : " << mu_ref << endl;
    	cout << "      mu_power_law   : " << mu_power_law << endl;
    	cout << "      SL_Tref        : " << SL_Tref << endl;
    	cout << "      SL_Sref        : " << SL_Sref << endl;
    	cout << "    Pr               : " << Pr << endl;
    	cout << "    PrTurb           : " << PrTurb << endl;
    	cout << "Body Force           : " << endl;
    	cout << "    gravity          : " << gravity[0] << " " << gravity[1] << " " << gravity[2] << endl;
    	cout << "Solver settings        " << endl;
    	cout << "    second order     : " << sndOrder << endl;
    	cout << "    1st order scalars: " << stOrderScalar << endl;
    	cout << "    gradient recon   : " << gradReconstr << ", value = " << gradreconstruction << endl;
    	cout << "    gradient limiter : " << limiterstr << ", value = ";
    	//IKJ
    	switch(limiterNavierS) {
    	case NOLIMITER:            cout << "NOLIMITER"<<endl;            break;
    	case BARTH_JESPERSEN_MOD:  cout << "BARTH_JESPERSEN_MOD"<<endl;  break;
    	default:
    		cout << limiterNavierS << endl;  break;
    	}
    	if (gradreconstruction == GRAD_SDWLS)
    		cout << "                     : " << "epsilonSDWLS = " << epsilonSDWLS << endl;
    	cout << "TIME_INTEGRATION     : " << tIntName << endl;
    	cout << "    nsteps           : " << nsteps << endl;
    	cout << "    timeStepMode     : " << timeStepMode << endl;
    	cout << "--------------------------------------------" << endl << endl;
    }

#ifdef USE_ARTIF_VISC
    // ----------------------------------------------------------------------------------------
    // Artificial viscosity for shock-capturing
    // ----------------------------------------------------------------------------------------
    turnOnArtifVisc = false;  // Default = false

    // Check if the user wants to use artificial viscosity for shock-capturing
    string tempStr;

    Param *pmy;
    if (getParam(pmy, "ARTIFICIAL_VISCOSITY")) {
    	// Example: ARTIFICIAL_VISCOSITY	 TURN_ON= YES  BULKVISC_ONLY= YES  SHOCK_ONLY= YES  SMOOTHSTEP_DIVERG_THRESH= 1.0e4  TYPE= STRAIN_BASE  COEFF= 3.25
    	tempStr = getParam("ARTIFICIAL_VISCOSITY")->getString("TURN_ON");
    	std::transform(tempStr.begin(), tempStr.end(), tempStr.begin(), ::toupper);
    	if(tempStr.compare("YES")==0 || tempStr.compare("Y")==0)
    		turnOnArtifVisc = true;
    }

    if(mpi_rank==0) {
    	cout<<"ARTIFICIAL_VISCOSITY: TURN_ON = "<<turnOnArtifVisc;
    	if(turnOnArtifVisc)
    		cout<<" (TRUE)"<<endl;
    	else
    		cout<<" (FALSE)"<<endl;
    }

    // Default values
    artifVisc_bulkViscOnly     = false;
    artifVisc_shockOnly        = false;
    artifVisc_smoothstepThresh = 1.0;
    artifVisc_type             = "NONE";
    artifVisc_coeff            = 0.0;

    artifVisc_mag = NULL;

    // Get the parameters for artificial viscosity
    if(turnOnArtifVisc) {
    	tempStr = getParam("ARTIFICIAL_VISCOSITY")->getString("BULKVISC_ONLY");
    	std::transform(tempStr.begin(), tempStr.end(), tempStr.begin(), ::toupper);
    	if(tempStr.compare("YES")==0 || tempStr.compare("Y")==0)
    		artifVisc_bulkViscOnly = true;

    	tempStr = getParam("ARTIFICIAL_VISCOSITY")->getString("SHOCK_ONLY");
    	std::transform(tempStr.begin(), tempStr.end(), tempStr.begin(), ::toupper);
    	if(tempStr.compare("YES")==0 || tempStr.compare("Y")==0)
    		artifVisc_shockOnly = true;
    	artifVisc_smoothstepThresh = getParam("ARTIFICIAL_VISCOSITY")->getDouble("SMOOTHSTEP_DIVERG_THRESH");

    	artifVisc_type = getParam("ARTIFICIAL_VISCOSITY")->getString("TYPE");

    	artifVisc_coeff = getParam("ARTIFICIAL_VISCOSITY")->getDouble("COEFF");

    	if(mpi_rank==0) {
    		cout<<"                      BULKVISC_ONLY = "<< artifVisc_bulkViscOnly << endl
    			<<"                      SHOCK_ONLY    = "<< artifVisc_shockOnly << "  (SMOOTHSTEP_DIVERG_THRESH = "<< artifVisc_smoothstepThresh << ")" <<endl
    			<<"                      TYPE          = "<< artifVisc_type <<endl
    			<<"                      COEFF         = "<< artifVisc_coeff <<endl
    			<<endl;
    	}

    	assert(artifVisc_mag == NULL);  	registerScalar(artifVisc_mag, "artifVisc_mag", CV_DATA);
    	assert(strMag == NULL);         	registerScalar(strMag,        "strMag",        CV_DATA);
    	assert(diverg == NULL);        		registerScalar(diverg,        "diverg",        CV_DATA);
    }

    // Statistics for artificial viscosity
    artifViscMagMin =  2.22e22;
    artifViscMagMax = -2.22e22;
#endif

    // ----------------------------------------------------------------------------------------
    // parse any WRITE_DATA parameters...
    // ----------------------------------------------------------------------------------------

    // init writedata
    initWriteData(this);

    // init probes //Moved to JoeWithModels.cpp
    //initProbes(this);
  }

  virtual void initializeFromRestartFile(const string &name)
  {
    // ----------------------------------------------------------------------------------------
    // read file
    // ----------------------------------------------------------------------------------------
    
    string ext = getFileNameExtension(name, ".");

    if (ext == "cdp")
    {
      CdpFilter cf(name);
      cf.minInitUgp(this);
    }
    else if ((ext == "msh") || (ext == "cas"))
    {
      MshFilter mf(name);
      mf.initUgpMin(this);
    }
    else
    {
      // assume a generic restart...
      readRestart(name);
    }

    UgpWithCv2Op::init(GG_GRADIENTS);
    
    // ----------------------------------------------------------------------------------------
    // init memory for face-based data
    // ----------------------------------------------------------------------------------------
    mul_fa    = new double[nfa_b2gg];
    lamOcp_fa = new double[nfa_b2gg];
    mut_fa    = new double[nfa_b2gg];

    // ----------------------------------------------------------------------------------------
    // init memory for scalar diffusion coeff
    // ----------------------------------------------------------------------------------------
    for (ScalarTranspEqIterator data = scalarTranspEqVector.begin(); data < scalarTranspEqVector.end(); data++)
    {
      data->diff     = new double[nfa];
      data->grad_phi = new double[ncv_g][3];
      if (data->reconstruction == "CONSERVATIVE")
      {
        data->rhophi      = new double[ncv_ggff];
        data->grad_rhophi = new double[ncv_g][3];
      }
    }


    // ----------------------------------------------------------------------------------------
    // gradients
    // ----------------------------------------------------------------------------------------
    if (sndOrder == true)
    {
      grad_rho = new double[ncv_g][3];
#ifdef temp_reconstruction
      grad_temp = new double[ncv_g][3];
#else
      grad_p = new double[ncv_g][3];
#endif
    }

    grad_u = new double[ncv_g][3][3];             // always allocate

    if ((mu_ref > 0.0) || (sndOrder == true))     // for now, will change
      grad_enthalpy = new double[ncv_g][3];
  }

  virtual ~UgpWithCvCompFlow()
  {
	if(mpi_rank==0)
		cout<<"~UgpWithCvCompFlow()"<<endl;
    if (mul_fa != NULL)       delete [] mul_fa;
    if (lamOcp_fa != NULL)    delete [] lamOcp_fa;
    if (mut_fa != NULL)       delete [] mut_fa;

//    for (ScalarTranspEqIterator data = scalarTranspEqVector.begin(); data < scalarTranspEqVector.end(); data++)
//    {
//      if (data->diff != NULL)
//        delete [] data->diff;
//    }

    if (grad_enthalpy != NULL)      delete [] grad_enthalpy;
    if (grad_rho != NULL)           delete [] grad_rho;
    if (grad_p != NULL)             delete [] grad_p;
    if (grad_u != NULL)             delete [] grad_u;

    if (petscSolver != NULL)        	delete petscSolver;
    if (petscSolverScalars != NULL) 	delete petscSolverScalars;

    if (nbocv_v_global != NULL) { delete [] nbocv_v_global; 	nbocv_v_global = NULL; }  // I don't know reason but this line will generate a memory crash in the destructor of Ugp (!)
  }

public:   // member variables
  // --------
  // For IKE
  // --------
  bool firstCall_JOEcalcViscousFluxScalar;

#ifdef USE_ARTIF_VISC
  // --------------------
  // Artificial viscosity
  // --------------------
  // Artificial viscosity for shock-capturing: Parameters
  bool turnOnArtifVisc;

  bool artifVisc_bulkViscOnly;       // Bulk-viscosity only, or Add to the molecular viscosity
  bool artifVisc_shockOnly;          // Do not apply to expansion wave (i.e. no artificial viscosity for positive divergence).
  double artifVisc_smoothstepThresh; // If artifVisc_shockOnly is active, a smooth-step function is applied between (-artifVisc_smoothstepThresh < diverg <0.0).
  string artifVisc_type;             // strain-magnitude based, divergence based
  double artifVisc_coeff;            // This is the base for artificial viscosity: the artif-visc is multiplied by this.

  // Artificial viscosity for shock-capturing: data for visualization
  double* artifVisc_mag;

  // Statistics for artificial viscosity
  double artifViscMagMin;
  double artifViscMagMax;
  double maxRatioArtifToReal;
#endif

  // ---------------------
  // Normal JOE variables
  // ---------------------
  double  *rho;             ///< density \f$\rho\f$
  double (*rhou)[3];        ///< momentum \f$\rho u\f$
  double  *rhoE;            ///< energy \f$\rho E\f$

  double (*vel)[3];         ///< velocity
  double *press;            ///< pressure
  double *temp;             ///< temperature
  double *enthalpy;         ///< chemical + sensible enthalpy (no kinetic energy!)

  double *RoM;              ///< gas constant at cell center
  double *gamma;            ///< ratio of specific heat at cell center
  double *sos;              ///< speed of sound at cell center
  
  double *mul_fa;           ///< laminar viscosity at cell faces
  double *lamOcp_fa;        ///< heat conductivity at cell faces
  
  double *mut_fa;           ///< turbulent viscosity at cell faces
  double *kine;             ///< turbulent kinetic energy

  double *residField;
  double *residField2;
  double *residField3;
  double *residField4;
  double *residField5;
  double *residField6;
  double resid_energ_th;

  double *strMag;           ///< magnitude strain rate
  double *vortMag;          ///< magnitude vorticity
  double *diverg;           ///< divergence

  // ----------------------------------------------------------------------------------------
  // boundary faces, allocate memory, TODO: should be incorporated in to the associated arrays ???
  // ----------------------------------------------------------------------------------------
//  double *rho_bfa, *T_bfa, (*vel_bfa)[3], *p_bfa, *h_bfa, *gam_bfa, *RoM_bfa;

  double (*grad_enthalpy)[3]; ///< enthalpy gradient, registered only if 2nd order euler flux calc turned on
  double (*grad_rho)[3];      ///< density gradient, registered only if 2nd order euler flux calc turned on
  double (*grad_u)[3][3];     ///< velocity gradient tensor, registered only 2nd order euler flux calc turned on or viscous calc
  double (*grad_p)[3];        ///< pressure gradient, registered only if 2nd order euler flux calc turned on and temp_reconstruction NOT defined
  double (*grad_temp)[3];     ///< temperature gradient, registered only if 2nd order euler flux calc turned on and temp_reconstruction defined
  
//  double *alpha_rho;          ///< limiter factor for gradients of rho, used to limit scalars in computation of Euler flux

  double p_ref, rho_ref, T_ref, R_gas;
  double GAMMA, Pr, PrTurb;

  string viscMode;
  double mu_ref;
// IKJ
  double d_scalar_comb_ref;

  double mu_power_law;
  double SL_Tref, SL_Sref;    ///< Sutherland law, Tref and S value, should be 273.15 and 110.4

  double gravity[3];

  double cfl, *local_dt;

  string timeStepMode;
  int nsteps, step;
  double time;


  int write_restart, check_interval;

  int gradreconstruction;
  bool sndOrder;
  bool stOrderScalar;
  int limiterNavierS;
  double epsilonSDWLS;
  
  int turbModel;
  enum TurbModel{NONE, SA, KEPS, KOM, KOMSST, V2F};

protected: 

  int linearSolverNS;         ///< linear solver for Navier-Stokes equations
  int linearSolverScalars;    ///< linear solver for scalar equations

  PetscSolver *petscSolver, *petscSolverScalars;   ///< petsc's linear solvers
  int *cvora5;
  int *nbocv_v_global;
  
  string initial_flowfield_output;       ///< tag (YES or NO) for output of initial field

public:   // member functions

  /**
   *  solve coupled linear system for Navier-Stokes equations
   */
  void solveCoupledLinSysNS(double (*phi)[5],
                            double (*Ap)[5][5],
                            double (*rhs)[5],
                            const double zeroAbs,
                            const double zeroRel,
                            const int maxIter)
  {
    switch (linearSolverNS)
    {
      case PETSC_GMRES:
        // on the first time, instantiate the petsc solver...
        if (petscSolver == NULL)
        {
          if (nbocv_v_global == NULL)
          {
            nbocv_v_global = new int[ncv_g];
            for (int icv = 0; icv < ncv; icv++)
              nbocv_v_global[icv] = cvora[mpi_rank] + icv;
            updateCvDataG1(nbocv_v_global, REPLACE_DATA);
          }
          petscSolver = new PetscSolver(cvora, nbocv_i, nbocv_v, 5);
          petscSolver->setTresholds(zeroAbs, zeroRel, maxIter);
        }
        petscSolver->solveGMRES(Ap, phi, rhs, cvora, nbocv_i, nbocv_v, nbocv_v_global, 5);

        break;

      case BCGSTAB:

        solveCvVectorR5Bcgstab(phi, Ap, rhs, zeroAbs, zeroRel, maxIter, "NS-eq");   // solve the system
        break;

      default:
        if (mpi_rank == 0)
          cerr << "Error: unrecognized solver: " << linearSolverNS << endl;
        throw(-1);
        break;
    }
  }

  void solveLinSysScalar(double *phi, double *Ap, double *rhs,
      const double zeroAbs, const double zeroRel, const int maxIter, char *scalarName)
  {
	  switch (linearSolverScalars)
	  {
	  case PETSC_GMRES:
		  // on the first time, instantiate the petsc solver...
		  if (petscSolverScalars == NULL)
		  {
			  if (nbocv_v_global == NULL)
			  {

				  nbocv_v_global = new int[ncv_g];
				  for (int icv = 0; icv < ncv; icv++)
					  nbocv_v_global[icv] = cvora[mpi_rank] + icv;
				  updateCvDataG1(nbocv_v_global, REPLACE_DATA);

			  }
			  petscSolverScalars = new PetscSolver(cvora, nbocv_i, nbocv_v, 1);
		  }

		  petscSolverScalars->setTresholds(zeroAbs, zeroRel, maxIter);
		  petscSolverScalars->solveGMRES(Ap, phi, rhs, cvora, nbocv_i, nbocv_v, nbocv_v_global);

		  break;

	  case BCGSTAB:
		  solveCvScalarBcgstab(phi, Ap, rhs, zeroAbs, zeroRel, maxIter, scalarName);
		  break;

	  default:
		  if (mpi_rank == 0)
			  cerr << "Error: unrecognized solver for scalars: " << linearSolverNS << endl;
		  throw(-1);
		  break;
	  }


  }
     /**
   *  solve coupled linear system for Navier-Stokes equations and scalars together
   */
  void solveCoupledLinSysNSCoupled(double **phi,
                                   double ***Ap,
                                   double **rhs,
                                   const double zeroAbs,
                                   const double zeroRel,
                                   const int maxIter,
                                   int nScal)
  {
    switch (linearSolverNS)
    {
      case PETSC_GMRES:

        // on the first time, instantiate the petsc solver...
        if (petscSolver == NULL)
        {
          if (nbocv_v_global == NULL)
          {
            nbocv_v_global = new int[ncv_g];
            for (int icv = 0; icv < ncv; icv++)
              nbocv_v_global[icv] = cvora[mpi_rank] + icv;
            updateCvDataG1(nbocv_v_global, REPLACE_DATA);
          }

          petscSolver = new PetscSolver(cvora, nbocv_i, nbocv_v, 5+nScal);
          petscSolver->setTresholds(zeroAbs, zeroRel, maxIter);
        }

        petscSolver->solveGMRESCoupled(Ap, phi, rhs, cvora, nbocv_i, nbocv_v, nbocv_v_global, 5+nScal);

        break;

      default:
        if (mpi_rank == 0)
          cerr << "Error: unrecognized solver: " << linearSolverNS << endl;
        throw(-1);
        break;
    }
  }

    void solveLinSysFull(double *phi, double *Ap, double *rhs, int *row_start, int *col_idx, const double zeroAbs, const double zeroRel, const int maxIter, int nVars)
  {
    switch (linearSolverScalars)
    {
      case PETSC_GMRES:
         {
          int *cvora_full, *index_full, *index;

          cvora_full =  new int[mpi_size+1];
          index_full =  new int[ncv_gg*nVars];
          index      =  new int[ncv_gg];

            for (int icv = 0; icv < ncv; icv++)
              index[icv] = cvora[mpi_rank] + icv;

            updateCvDataG1G2(index, REPLACE_DATA);

          for (int i=0; i<=mpi_size; i++)
                cvora_full[i] = cvora[i]*nVars;

	  int count =0;
          for (int icv=0; icv<ncv_gg; icv++) {
             for (int i=0; i<nVars; i++) {
                index_full[count]= index[icv]*nVars+i;
		count++;
             }
	  }

        // on the first time, instantiate the petsc solver...
          if (petscSolverScalars == NULL)
             petscSolverScalars = new PetscSolver(cvora_full, row_start, col_idx, 1);

          petscSolverScalars->setTresholds(zeroAbs, zeroRel, maxIter);
          petscSolverScalars->solveGMRES(Ap, phi, rhs, cvora_full, row_start, col_idx, index_full);
        }
        break;

      case BCGSTAB:
        solveCvScalarBcgstab(phi, Ap, rhs, zeroAbs, zeroRel, maxIter, "Full");
        break;

      default:
        if (mpi_rank == 0)
          cerr << "Error: unrecognized solver for scalars: " << linearSolverNS << endl;
        throw(-1);
        break;
    }
  }



  /**
   *    calc time step
   */
  virtual double calcDt(double cfl = 0.0);


  // IKJ
  /*
   * Method: calcDtForCFLone
   * -----------------------
   * Original code = UgpWithCvCompFlow::calcDt()
   * Calculate Dt corresponding to CFL=1.0 using the "MAX_BASED" formulation
   * This will not affect the "local_dt" variables, so it is safe to call this function anytime
   */
  double calcDtForCFLone();

  /**
   *
   *  explicit euler flux HLLC
   *
   */
  virtual int calcEulerFlux_HLLC(double &Frho, double *Frhou, double &FrhoE, double *FrhoScal,
      const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *ScalL, const double kL,
      const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *ScalR, const double kR,
      const double area, const double *nVec, const int nScal, const double surfVeloc);


  /**
   *
   *  implicit euler flux HLLC supersonic
   *
   */
  virtual int calcEulerFluxMatrices_HLLC(double (*A_L)[5], double (*A_R)[5], double (*A_L_Scal)[6], double (*A_R_Scal)[6],
      const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *scalL, const double kL,
      const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *scalR, const double kR,
      const double area, const double *nVec, const int nScal, const double surfVeloc);
  

  /**
   *
   *  implicit euler flux for HLLC subsonic
   *
   */
  void calcSubSonicJacobeanHLLC(double (*AL)[5], double (*AR)[5],
      double rhoL, const double *uL, double pL, double eL, double qL, double psiL, double SL,
      double rhoSL, double *rhouSL, double eSL,
      double *dSMdUL, double *dSMdUR, double *dpsdUL, double *dpsdUR, double SM, double pS,
      double gamma, const double *nV);

  // nV is not normalized
  void calcJacobianA(double (*A)[5], const double *vel, double pp, double rrho,
      const double *nV, double gamma, double surfVeloc);
  

  /**
   *
   * viscous flux routine for implicit
   * needs lsg_coeff0 and lsg_coeff1 ---> LS gradient coefficients
   *
   */
#ifdef USE_ARTIF_VISC
  virtual void addViscFlux(double *Frhou, double &FrhoE, double (*A0)[5], double (*A1)[5],
      const double rho0, const double *u0, const double (&grad_u0)[3][3], const double h0, const double *grad_h0, const double T0, const double R0, const double gam0, const double kine0,
      const double rho1, const double *u1, const double (&grad_u1)[3][3], const double h1, const double *grad_h1, const double T1, const double R1, const double gam1, const double kine1,
      const double mul, const double mut, const double lambdaOverCp, const double kine_fa, const double *u_fa,
      const double area, const double *nVec, const double smag, const double *sVec, const double artifBulkViscosity = 0.0);
#else
  virtual void addViscFlux(double *Frhou, double &FrhoE, double (*A0)[5], double (*A1)[5],
      const double rho0, const double *u0, const double (&grad_u0)[3][3], const double h0, const double *grad_h0, const double T0, const double R0, const double gam0, const double kine0,
      const double rho1, const double *u1, const double (&grad_u1)[3][3], const double h1, const double *grad_h1, const double T1, const double R1, const double gam1, const double kine1,
      const double mul, const double mut, const double lambdaOverCp, const double kine_fa, const double *u_fa, 
      const double area, const double *nVec, const double smag, const double *sVec);
#endif
  
  /**
   * loop over scalars and solve them
   */
  void setScalarBC(FaZone *zone);
  
  virtual void UserDefinedScalarClipping(const string &name)  {/*empty*/}
  
  void calcViscousFluxScalar(double *rhs_rhoScal, double *Ascal, ScalarTranspEq &transpScal, int flagImplicit);


// *************** COUPLED SOLUTIONS (UNDER DEBUGGING) ******************************


 
  /**
   *  explicit euler flux HLLC, coupled NSE and scalars
   */
  virtual void calcEulerFluxCoupled_HLLC(double *EulerFlux,
      const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *ScalL, const double kL,
      const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *ScalR, const double kR,
      const double area, const double *nVec, const int nScal, const double *ConvTerm, const double surfVeloc, const string BC_treatment);

  /**
   *  matrices for implicit euler flux HLLC, coupled NSE and scalars
   */
  virtual void calcEulerFluxMatricesCoupled_HLLC(double **A_L, double **A_R, 
      const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *scalL, const double *dpress_dscalL, const double kL,
      const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *scalR, const double *dpress_dscalR, const double kR,
      const double area, const double *nVec, const int nScal, const double *ConvTerm, const double surfVeloc, const string BC_treatment);

  /**
   *  supersonic Jacobian matrix for Euler flux, coupled NSE and scalars
   */
  void calcJacobianCoupled_HLLC(double **A, const double *vel, const double pp, const double rrho, const double hh, const double *scal, const double *dpress_dscal,
      const double *nV, const double gamma, const double surfVeloc, const int nScal);

  /**
   *  derivative of contact surface speed
   */  
  void calcdSMdU_HLLC(double *dSMdU, const double *vel, const double pp, const double rrho, const double hh, const double *dpress_dscal,
                 const double *nV, const double gamma, const double SM, const double S, const double invertrho, const double Factor, const int nScal);

  /**
   *  subsonic Jacobian for Euler flux, coupled NSE and scalars
   */  
  void calcSubSonicJacobianCoupled_HLLC(double **A, const double *dSMdU, const double rhoS, const double *rhouS, const double rhoeS, const double *scal,
                                        const double pS, const double Omega, const double sM, const double rhoSmq, const double *nV, const double area, const int nScal);
  
  /**
   *  explicit and implicit viscous flux, coupled NSE and scalars
   */
  virtual void calcViscousFluxCoupled(double *ViscousFlux, double **A0, double **A1,
      const double rho0, const double *u0, const double (&grad_u0)[3][3], const double h0, const double *grad_h0, const double T0, const double R0, const double gam0, const double *Scal0, const double (*gradScal0)[3], const double *dpress_dscal0, const double kine0,
      const double rho1, const double *u1, const double (&grad_u1)[3][3], const double h1, const double *grad_h1, const double T1, const double R1, const double gam1, const double *Scal1, const double (*gradScal1)[3], const double *dpress_dscal1, const double kine1,
      const double mul, const double mut, const double lambdaOverCp, const double kine_fa, const double *u_fa, const double *diff, const double *DiffTerm,
      const double area, const double *nVec, const double smag, const double *sVec, const double alpha, const int nScal);

  

  
  

// *************** COUPLED SOLUTIONS (UNDER DEBUGGING) ******************************
  
  virtual void calcGradVel()
  {
    calcCv2Grad(grad_u, vel, limiterNavierS, sos, epsilonSDWLS);
  }

  virtual void calcStrainRateAndDivergence()
  {
    if ((strMag == NULL) || (diverg == NULL))
    {
      cout << "you have tried to calculate the magnitude of the strain rate and the divergence of U, but one of them has not been registered!" << endl;
      throw(-1);
    }

    for (int icv=0; icv<ncv_g; icv++)
    {
      diverg[icv] = grad_u[icv][0][0] + grad_u[icv][1][1] + grad_u[icv][2][2];

      strMag[icv] = 0.0;
      for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        if (i == j)  strMag[icv] += pow(0.5*(grad_u[icv][i][j] + grad_u[icv][j][i])-1./3.*diverg[icv], 2.0);
        else         strMag[icv] += pow(0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]), 2.0);

      strMag[icv] = sqrt(2.0*strMag[icv]);
    }
  }

  virtual void calcDivergence()
  {
    if (diverg == NULL)
    {
      cout << "you have tried to calculate the divergence of U, but it has not been registered!" << endl;
      throw(-1);
    }

    for (int icv=0; icv<ncv_g; icv++)
      diverg[icv] = grad_u[icv][0][0] + grad_u[icv][1][1] + grad_u[icv][2][2];
  }
  
  virtual void calcStrainRate()
  {
    if (strMag == NULL)
    {
      cout << "you have tried to calculate the magnitude of the strain rate, but it has not been registered!" << endl;
      throw(-1);
    }

    for (int icv=0; icv<ncv_g; icv++)
    {
      strMag[icv] = 0.0;
      double diverg = grad_u[icv][0][0] + grad_u[icv][1][1] + grad_u[icv][2][2];

      for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        if (i == j)  strMag[icv] += pow(0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]) - 1./3.*diverg, 2.0);
        else         strMag[icv] += pow(0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]), 2.0);

      strMag[icv] = sqrt(2.0*strMag[icv]);
    }
  }

  virtual void calcVorticity()
  {
    if (vortMag == NULL)
    {
      cout << "you have tried to calculate the magnitude of the vorticity, but it has not been registered!" << endl;
      throw(-1);
    }

    for (int icv=0; icv<ncv_g; icv++)
    {
      vortMag[icv] = 0.0;
      for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        vortMag[icv] += pow(0.5*(grad_u[icv][i][j] - grad_u[icv][j][i]), 2.0);

      vortMag[icv] = sqrt(2.0*vortMag[icv]);
    }
  }


  /** 
   * calculate wall distance: minimum distance of cv to wall face by looping over all cv's and faces 
   */
  void calcWallDistance(double *phi, double *wd)
  {
    if (mpi_rank == 0)
      cout << "calc wall distance" << endl;

    double wtime, wtime0;
    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0)
      wtime = MPI_Wtime();

    for (int icv=0; icv<ncv; icv++)
      wd[icv] = 1.0e20;

    // count the walls for each rank
    int my_wall_faces = 0;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        Param *param;
        if (getParam(param, zone->getName()))
          if (param->getString() == "WALL")
          {
            for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
              my_wall_faces++;
          }
      }

    int tot_wall_faces = 0;
    MPI_Allreduce(&my_wall_faces, &tot_wall_faces, 1, MPI_INT, MPI_SUM, mpi_comm);

    if (tot_wall_faces == 0)
    {
      if (mpi_rank == 0)
        cerr << "ERROR: you have been trying to calculate the wall distance, but there are no walls in you domain ?!?" << endl;
      throw(-1);
    }

    // set up send side
    int send_count[mpi_size], send_displ[mpi_size];
    for (int r=0; r<mpi_size; r++)
    {
      send_count[r] = my_wall_faces*3;
      send_displ[r] = 0;
    }

    // set up receive side
    int recv_count[mpi_size], recv_displ[mpi_size];
    int wallcoord = my_wall_faces*3;
    MPI_Allgather(&wallcoord, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
    
    recv_displ[0] = 0;
    for (int i = 1; i < mpi_size; i++)
      recv_displ[i] = recv_count[i-1] + recv_displ[i-1];

    // fill up send buffer
    int count = 0;
    double *my_wall_fa        = new double[my_wall_faces*3];
    double *my_wall_fa_normal = new double[my_wall_faces*3];
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        Param *param;
        if (getParam(param, zone->getName()))
          if (param->getString() == "WALL")
            for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
            {
              for (int i=0; i<3; i++)
              {
                my_wall_fa[3*count + i] = x_fa[ifa][i];
                my_wall_fa_normal[3*count + i] = fa_normal[ifa][i];
              }
  
              count++;
            }
      }
    
    // receive all wall face coordinates
    double *wall_fa   = new double[tot_wall_faces*3];
    double *wall_fa_n = new double[tot_wall_faces*3];
    
    MPI_Alltoallv(my_wall_fa, send_count, send_displ, MPI_DOUBLE, 
                  wall_fa,    recv_count, recv_displ, MPI_DOUBLE, mpi_comm);

    MPI_Alltoallv(my_wall_fa_normal, send_count, send_displ, MPI_DOUBLE, 
                  wall_fa_n,         recv_count, recv_displ, MPI_DOUBLE, mpi_comm);

    delete [] my_wall_fa;
    delete [] my_wall_fa_normal;


    // find minimum wall distance 
    for (int icv=0; icv<ncv; icv++)
    {
      double minDist = 1.0e20;
      int ifaSave;

      for (int ifa=0; ifa<tot_wall_faces; ifa++)
      {
        double dist = sqrt(   pow(x_cv[icv][0]-wall_fa[3*ifa  ], 2.0)
                            + pow(x_cv[icv][1]-wall_fa[3*ifa+1], 2.0)
                            + pow(x_cv[icv][2]-wall_fa[3*ifa+2], 2.0));
        
        if (dist < minDist)
        {
          ifaSave = ifa;
          minDist = dist;
        }
      }
      
/*      double fa_x[3] = {wall_fa[3*ifaSave], wall_fa[3*ifaSave+1], wall_fa[3*ifaSave+2]};
      double fa_n[3] = {wall_fa_n[3*ifaSave], wall_fa_n[3*ifaSave+1], wall_fa_n[3*ifaSave+2]};
      double n[3], s_half[3];
      
      normVec3d(n, fa_n);
      vecMinVec3d(s_half, fa_x, x_cv[icv]);
      wd[icv] = fabs(vecDotVec3d(s_half, n));*/
      
      double fa_x[3] = {wall_fa[3*ifaSave], wall_fa[3*ifaSave+1], wall_fa[3*ifaSave+2]};
      double fa_n[3] = {wall_fa_n[3*ifaSave], wall_fa_n[3*ifaSave+1], wall_fa_n[3*ifaSave+2]};
      double n[3], s_half[3];
      
      double nmag = normVec3d(n, fa_n);
      vecMinVec3d(s_half, fa_x, x_cv[icv]);

      double sMag2 = vecDotVec3d(s_half, s_half);
      double dNorm2 = vecDotVec3d(s_half, n)*vecDotVec3d(s_half, n);
      double dTang2 = sMag2-dNorm2;

      double specRadius2 = 0.5*nmag;
      

      if (dTang2 > specRadius2)   wd[icv] = minDist;
      else                        wd[icv] = fabs(vecDotVec3d(s_half, n));


//      wd[icv] = minDist;
    }
    /*
    // correct closest wall
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      Param *param;
      if (getParam(param, zone->getName()))
        if (param->getString() == "WALL")
        for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
        {
          int icv0 = cvofa[ifa][0];
          
          double n[3], s_half[3];
          double nmag = normVec3d(n, fa_normal[ifa]);
          vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
          wd[icv0] = fabs(vecDotVec3d(s_half, n));  //normVec3d(s_half);          
        }
    }*/
    
    MPI_Barrier(mpi_comm);
    
    // update wall distance into ghosts
    updateCvDataG1(wd, REPLACE_DATA);


    if (mpi_rank == 0) 
    {
      double wtime0 = wtime;
      wtime = MPI_Wtime();
      cout << " > time to compute wall distance [s]: " << wtime - wtime0 << endl;
    }
    
    delete [] wall_fa;
    delete [] wall_fa_n;    
  }
  
  // Interpolate variable at the cell center from values at the faces (e.g., for viscosity)
  double InterpolateAtCellCenterFromFaceValues(double *phi, int icv)
  {
    double phiC = 0.0;

    int foc_f = faocv_i[icv];
    int foc_l = faocv_i[icv + 1] - 1;
    for (int foc=foc_f; foc<=foc_l; foc++)
    {
      int ifa = faocv_v[foc];
      phiC += phi[ifa];
    }
    phiC /= (double)(foc_l-foc_f+1);    
    
    return phiC;
  }
  
  
  
  
  
  
  
  // ###########################################################################
  // ---------------------------------------------------------------------------
  //
  // turbmodels turbmodels turbmodels turbmodels turbmodels turbmodels 
  //
  // ---------------------------------------------------------------------------
  // ###########################################################################
  
  virtual void calcRansTurbViscMuet()
  {
    static int flag = 1;
    
    // provide zero mut for laminar calculations
    if (flag == 1)
    {
      flag = 0;
      for (int ifa=0; ifa<nfa; ifa++)
        mut_fa[ifa] = 0.0;
    }
  }

  virtual void initialHookScalarRansTurbModel() {/*empty*/}
  virtual void diffusivityHookScalarRansTurb(const string &name)  {/*empty*/}
  virtual void sourceHookScalarRansTurb(double *rhs, double *A, const string &name, int flagImplicit) {/*empty*/}
  virtual void boundaryHookScalarRansTurb(double *phi_ph, FaZone *zone, const string &name)  {/*empty*/}  
  virtual void sourceHookRansTurb(double *rhsRho, double (*rhsRhou)[3], double *rhsRhoE, double (*A)[5][5]) {/*empty*/}

  virtual void sourceHookRansTurbCoupled(double **rhs, double ***A, int nScal, int flagImplicit) {/*empty*/}

  virtual void pressureDerivativeHookScalarRansTurb() {/*empty*/}
  
  

  
  
  // ###########################################################################
  // ---------------------------------------------------------------------------
  //
  // thermodynamics thermodynamics thermodynamics thermodynamics 
  //
  // ---------------------------------------------------------------------------  
  // ###########################################################################
  
  virtual double calcMuLam(double temp)
  {
    double muLam = 0.0;
    
    if (viscMode == "SUTHERLAND")      muLam = mu_ref*pow(temp/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(temp + SL_Sref); 
    else if (viscMode == "POWERLAW")   muLam = mu_ref*pow(temp/T_ref, mu_power_law);
    else muLam = 0.0;

    return muLam;
  }
  
  virtual double calcMuLam(int icv)
  {
    double muLam = 0.0;
    
    if (viscMode == "SUTHERLAND")      muLam = mu_ref*pow(temp[icv]/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(temp[icv] + SL_Sref); 
    else if (viscMode == "POWERLAW")   muLam = mu_ref*pow(temp[icv]/T_ref, mu_power_law);
    else muLam = 0.0;

    return muLam;
  }
  
  
  /*! \brief Compute constants gas constant R_gas and heat capacity ratio gamma.
   *
   *  Compute R_gas and gamma and initialize the arrays RoM and gam.
   */
  virtual void initialHookScalarRansCombModel()
  {
    for (int icv = 0; icv < ncv_gg; icv++)
    {
      RoM[icv] = R_gas;
      gamma[icv] = GAMMA;
    }
  }



  /*! \brief Compute for each cell the states of variables and material properties based on constant gamma and R/M.
   *
   *  Compute temperature temp, viscosity mul, heat coefficient lambda,
   *  enthalpy enthalp and pressure press.
   */
  virtual void calcRansStateVarAndMaterialProperties(double *rho, double(*rhou)[3], double *rhoE)
  {
    // Compute velocity, pressure, temperature, enthalpy and speed of sound at cell centers
    calcStateVariables();
    
    // Compute viscosity and thermal diffusivity at cell faces
    calcMaterialProperties();
  }

  /*! \brief Compute for each cell the states of variables based on constant gamma and R/M.
   *
   *  Compute velocity vel, temperature temp, enthalpy enthalpy, pressure press, specific heat ratio gamma, speed of sounds sos and gas constant R.
   *  R and gamma not actually computed since constant
   *  
   */
  virtual void calcStateVariables()
  {
	vector<std::string> warningsVectorUgp; // If the simulation breaks down for some reasons,
			                               // it will dump out a massive amount of warning/error messages.
			                               // Thus, this object stores the (only) warning messages so that
			                               // the user can manipulate the amount of messages that should be
			                               // printed out on the screen.
    for (int icv = 0; icv < ncv_gg; icv++)
    {
      if (rho[icv] <= 0.0) {
//        cout << "negative density at xcv: " << x_cv[icv][0] << ", " << x_cv[icv][1] << ", " << x_cv[icv][2] << endl;
    	  std::stringstream ss;
    	  ss<<"WARNING! UgpWithCvCompFlow::calcStateVariables(): Negative density at xcv: "
    		<< x_cv[icv][0] << ", "<< x_cv[icv][1] << ", " << x_cv[icv][2] << ": "
    		<< rho[icv] << endl;

    	  warningsVectorUgp.push_back(ss.str());
      }

      for (int i=0; i<3; i++)
        vel[icv][i] = rhou[icv][i]/rho[icv];

      double kinecv = 0.0;
      if (kine != NULL)
        kinecv = kine[icv];

      double pr = (gamma[icv] - 1.0) * (rhoE[icv] - 0.5 * vecDotVec3d(rhou[icv], rhou[icv]) / rho[icv] - rho[icv] * kinecv);
      if (pr <= 0.0)
      {
    	  std::stringstream ss;
    	  ss<<"WARNING! UgpWithCvCompFlow::calcStateVariables(): Negative pressure at xcv: "
    			  << x_cv[icv][0] << ", " << x_cv[icv][1] << ", " << x_cv[icv][2] << ": " << pr<<endl;
    	  ss<<std::scientific;
    	  ss<<std::setprecision(3);
    	  ss<<"                                                  rhoE="<<rhoE[icv]
    	    <<", rho="<<rho[icv]
    	    <<", rhou=("<<rhou[icv][0]<<","<<rhou[icv][1]<<","<<rhou[icv][2]<<")"
    	    <<", gamma="<<gamma[icv]
    	    <<", kinecv="<<kinecv<<endl;

    	  warningsVectorUgp.push_back(ss.str());
      }
      else
        press[icv] = pr;

      temp[icv] = press[icv] / (rho[icv] * RoM[icv]);
      enthalpy[icv] = gamma[icv] * RoM[icv] / (gamma[icv] - 1.0) * temp[icv];
      sos[icv] = sqrt(fabs(gamma[icv] * press[icv] / rho[icv]));
    }

    int myCountWM = warningsVectorUgp.size();
    int totCountWM;
    MPI_Allreduce(&myCountWM, &totCountWM, 1, MPI_INT, MPI_SUM, mpi_comm);

    if(totCountWM <= MAX_WARNING_MESSAGES_UGP) {
        if(mpi_rank != 0) { MPI_Status status; int dummy; MPI_Recv(&dummy, 1, MPI_INT, mpi_rank-1, 1117, mpi_comm, &status); }

        if(!warningsVectorUgp.empty()) {
            for(size_t i=0; i<warningsVectorUgp.size(); ++i)
                cout<<warningsVectorUgp[i];
        }

        if(mpi_rank < mpi_size-1) { int dummy; MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 1117, mpi_comm); }
    } else {
        // Get the number of warning messages per each cpu
        int *countWMarray = new int [mpi_size];
        MPI_Gather(&myCountWM, 1,  MPI_INT, countWMarray, 1, MPI_INT, 0, mpi_comm); 
       
        MPI_Barrier(mpi_comm);
        for(int i=0; i<54321; ++i) {}
        MPI_Barrier(mpi_comm);

        if(mpi_rank==0) {
            // Get the number of cpus that have more than one error message, minimum nubmer of error messages, and maximum
            int countNonzero = 0;
            int minNonzero = 100*ncv;
            int maxNonzero = 0;
            for(int impi=0; impi<mpi_size; ++impi) {
                if(countWMarray[impi]>0) {
                    ++countNonzero;
                    minNonzero = min(countWMarray[impi], minNonzero);
                    maxNonzero = max(countWMarray[impi], maxNonzero);
                }
            }
            assert(countNonzero>0);
            
            // Number of messages per ecah cpu: a + b[mpi_rank]
            int a = min(int(MAX_WARNING_MESSAGES_UGP/countNonzero), minNonzero);
            int left = MAX_WARNING_MESSAGES_UGP - a*countNonzero;
            for(int impi=0; impi<mpi_size; ++impi) {
                if(countWMarray[impi]>0) {
                    if(left > 0) {
                        countWMarray[impi] = a + 1;
                        --left;
                    } else
                        countWMarray[impi] = a;
                }
            }

            cout<<"TOO MANY WARNINGS ("<< totCountWM <<") FROM UgpWithCvCompFlow::calcStateVariables() at TOTAL "<<countNonzero<<" CPUS"<<endl
                <<"SHOW ONLY FIRST " << MAX_WARNING_MESSAGES_UGP << " WARNING MESSGAES: "<<endl;
        }

        MPI_Bcast(countWMarray, mpi_size, MPI_INT, 0, mpi_comm);

        // Show warning messages
        if(mpi_rank != 0) { MPI_Status status; int dummy; MPI_Recv(&dummy, 1, MPI_INT, mpi_rank-1, 1222, mpi_comm, &status); }

        if(!warningsVectorUgp.empty()) {
            for(size_t i=0; i<max(0, min(countWMarray[mpi_rank], (int) warningsVectorUgp.size())); ++i)
                cout<<warningsVectorUgp[i];
        }

        if(mpi_rank < mpi_size-1) { int dummy; MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 1222, mpi_comm); }

        // Clear memory
        delete [] countWMarray;
    }
    MPI_Barrier(mpi_comm);
  }

  virtual void calcStateVariables(bool showWarningErrorMessages)
  {
	vector<std::string> warningsVectorUgp; // If the simulation breaks down for some reasons,
			                               // it will dump out a massive amount of warning/error messages.
			                               // Thus, this object stores the (only) warning messages so that
			                               // the user can manipulate the amount of messages that should be
			                               // printed out on the screen.
    for (int icv = 0; icv < ncv_gg; icv++)
    {
      if (rho[icv] <= 0.0 && showWarningErrorMessages) {
    	  std::stringstream ss;
    	  ss<<"WARNING! UgpWithCvCompFlow::calcStateVariables(): Negative density at xcv: "
    		<< x_cv[icv][0] << ", "<< x_cv[icv][1] << ", " << x_cv[icv][2] << ": "
    		<< rho[icv] << endl;

    	  warningsVectorUgp.push_back(ss.str());
      }

      for (int i=0; i<3; i++)
        vel[icv][i] = rhou[icv][i]/rho[icv];

      double kinecv = 0.0;
      if (kine != NULL)
        kinecv = kine[icv];

      double pr = (gamma[icv] - 1.0) * (rhoE[icv] - 0.5 * vecDotVec3d(rhou[icv], rhou[icv]) / rho[icv] - rho[icv] * kinecv);
      if (pr <= 0.0)
        if(showWarningErrorMessages) {
    	  std::stringstream ss;
    	  ss<<"WARNING! UgpWithCvCompFlow::calcStateVariables(): Negative pressure at xcv: "
    			  << x_cv[icv][0] << ", " << x_cv[icv][1] << ", " << x_cv[icv][2] << ": " << pr<<endl;
    	  ss<<std::scientific;
    	  ss<<std::setprecision(3);
    	  ss<<"                                                  rhoE="<<rhoE[icv]
    	    <<", rho="<<rho[icv]
    	    <<", rhou=("<<rhou[icv][0]<<","<<rhou[icv][1]<<","<<rhou[icv][2]<<")"
    	    <<", gamma="<<gamma[icv]
    	    <<", kinecv="<<kinecv<<endl;

    	  warningsVectorUgp.push_back(ss.str());
        }
      else
        press[icv] = pr;

      temp[icv]     = press[icv] / (rho[icv] * RoM[icv]);
      enthalpy[icv] = gamma[icv] * RoM[icv] / (gamma[icv] - 1.0) * temp[icv];
      sos[icv]      = sqrt(fabs(gamma[icv] * press[icv] / rho[icv]));
    }

    if(showWarningErrorMessages) {
        int myCountWM = warningsVectorUgp.size();
        int totCountWM;
        MPI_Allreduce(&myCountWM, &totCountWM, 1, MPI_INT, MPI_SUM, mpi_comm);

        if(totCountWM <= MAX_WARNING_MESSAGES_UGP) {
            if(mpi_rank != 0) { MPI_Status status; int dummy; MPI_Recv(&dummy, 1, MPI_INT, mpi_rank-1, 1117, mpi_comm, &status); }

            if(!warningsVectorUgp.empty()) {
                for(size_t i=0; i<warningsVectorUgp.size(); ++i)
                    cout<<warningsVectorUgp[i];
            }

            if(mpi_rank < mpi_size-1) { int dummy; MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 1117, mpi_comm); }
        } else {
            // Get the number of warning messages per each cpu
            int *countWMarray = new int [mpi_size];
            MPI_Gather(&myCountWM, 1,  MPI_INT, countWMarray, 1, MPI_INT, 0, mpi_comm); 
           
            MPI_Barrier(mpi_comm);
            for(int i=0; i<54321; ++i) {}
            MPI_Barrier(mpi_comm);

            if(mpi_rank==0) {
                // Get the number of cpus that have more than one error message, minimum nubmer of error messages, and maximum
                int countNonzero = 0;
                int minNonzero = 100*ncv;
                int maxNonzero = 0;
                for(int impi=0; impi<mpi_size; ++impi) {
                    if(countWMarray[impi]>0) {
                        ++countNonzero;
                        minNonzero = min(countWMarray[impi], minNonzero);
                        maxNonzero = max(countWMarray[impi], maxNonzero);
                    }
                }
                assert(countNonzero>0);
                
                // Number of messages per ecah cpu: a + b[mpi_rank]
                int a = min(int(MAX_WARNING_MESSAGES_UGP/countNonzero), minNonzero);
                int left = MAX_WARNING_MESSAGES_UGP - a*countNonzero;
                for(int impi=0; impi<mpi_size; ++impi) {
                    if(countWMarray[impi]>0) {
                        if(left > 0) {
                            countWMarray[impi] = a + 1;
                            --left;
                        } else
                            countWMarray[impi] = a;
                    }
                }

                cout<<"TOO MANY WARNINGS ("<< totCountWM <<") FROM UgpWithCvCompFlow::calcStateVariables() at TOTAL "<<countNonzero<<" CPUS"<<endl
                    <<"SHOW ONLY FIRST " << MAX_WARNING_MESSAGES_UGP << " WARNING MESSGAES: "<<endl;
            }

            MPI_Bcast(countWMarray, mpi_size, MPI_INT, 0, mpi_comm);

            // Show warning messages
            if(mpi_rank != 0) { MPI_Status status; int dummy; MPI_Recv(&dummy, 1, MPI_INT, mpi_rank-1, 1222, mpi_comm, &status); }

            if(!warningsVectorUgp.empty()) {
                for(size_t i=0; i<max(0, min(countWMarray[mpi_rank], (int) warningsVectorUgp.size())); ++i)
                    cout<<warningsVectorUgp[i];
            }

            if(mpi_rank < mpi_size-1) { int dummy; MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 1222, mpi_comm); }

            // Clear memory
            delete [] countWMarray;
        }
    }
    MPI_Barrier(mpi_comm);
  }


  /*! \brief Compute for each face the material properties.
   *
   *  Compute laminar viscosity mul_fa and heat conductivity lamOcp_fa at the faces.
   *  
   */
  virtual void calcMaterialProperties() {
	  if (mu_ref > 0.0) {
		  if (viscMode == "SUTHERLAND") {
			  // internal faces
			  for (int ifa = nfa_b; ifa < nfa; ifa++) {
				  int icv0 = cvofa[ifa][0];
				  int icv1 = cvofa[ifa][1];

				  double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
				  vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
				  vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
				  double w0 = sqrt(vecDotVec3d(dx0, dx0));
				  double w1 = sqrt(vecDotVec3d(dx1, dx1));

				  double temperature  = (w1*temp[icv0] + w0*temp[icv1])/(w0+w1);
				  mul_fa[ifa] = mu_ref*pow(temperature/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(temperature + SL_Sref);
				  lamOcp_fa[ifa] = mul_fa[ifa] / Pr;
			  }

			  // boundary faces computed next in setBC
		  }
		  else if (viscMode == "POWERLAW") {
			  // internal faces
			  for (int ifa = nfa_b; ifa < nfa; ifa++)
			  {
				  int icv0 = cvofa[ifa][0];
				  int icv1 = cvofa[ifa][1];

				  double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
				  vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
				  vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
				  double w0 = sqrt(vecDotVec3d(dx0, dx0));
				  double w1 = sqrt(vecDotVec3d(dx1, dx1));

				  double temperature  = (w1*temp[icv0] + w0*temp[icv1])/(w0+w1);
				  mul_fa[ifa] = mu_ref*pow(temperature/T_ref, mu_power_law);
				  lamOcp_fa[ifa] = mul_fa[ifa] / Pr;
			  }

			  // boundary faces computed next in setBC
		  }
		  else {
			  cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
			  throw(-1);
		  }
	  }

	  // The part for calculating the artificial viscosity is moved to JoeWithModels.cpp -- Here, grad_u has not been calculated yet.
  }

#ifdef USE_ARTIF_VISC
  // \ Calculation of artificial viscosity.
  // \ It takes the parameters only related to the artificial viscosity, but it can access to the all the member variables of UgpWithCvCompFlows
  virtual void calcArtifVisc(double* artifVisc_mag,
		  const bool artifVisc_bulkViscOnly, const bool artifVisc_shockOnly,
		  const double artifVisc_smoothstepThresh, const string &artifVisc_type, const double artifVisc_coeff) {
	  static bool firstCall = true;

	  assert(diverg != NULL && strMag != NULL);
	  assert(grad_rho != NULL); // grad_rho is not defined in a first-order calculation
	  assert(artifVisc_mag != NULL);

	  calcStrainRateAndDivergence();
	  if(firstCall || !sndOrder)
		  calcCv2Grad(grad_rho, rho, limiterNavierS, rho, epsilonSDWLS);
	  MPI_Barrier(mpi_comm);

	  double myArtifViscMagMin =  2.2e22;
	  double myArtifViscMagMax = -2.2e22;

	  // Calculate artifVisc_mag on the CV cells
	  for(int icv=0; icv<ncv_g; ++icv) {
		  // If SHOCK_ONLY, clip for expansion waves and weak shocks
		  double clipFunc = 1.0;
		  if(artifVisc_shockOnly) {
			  if(diverg[icv] > 0.0) // expansion waves
				  clipFunc = 0.0;
			  else {
				  if(diverg[icv] > artifVisc_smoothstepThresh) {
					  double divergNormalized = -diverg[icv] / artifVisc_smoothstepThresh;
					  clipFunc = 3.0*pow(divergNormalized, 2.0) - 2.0*pow(divergNormalized, 3.0);
					  	  // A smooth-step function whose domain is between 0.0 and 1.0
					  	  // (For details, search Wikipedia with "smoothstep")
				  }
			  }
		  }
		  if(clipFunc<0.0 || clipFunc>1.0)
			  cout<<"WARNING UgpWithCvCompFlow::calcArtifVisc(): The clip function is not in the correct range (0~1) = "<<clipFunc<<endl;

		  // Calculate grid-size (currently MAX_BASE... You can also try the averages of the six grid-scales).
		  // We are only considering the grid spacing projected in the shock normal direction (grad_rho direction)
		  // in order to remove the span-wise grid spacing in a 2D calculation.
		  double delta = 0.0;

//		  if(grad_rho != NULL) {
//			  double unitGradRho[3];
//			  double gradRhoMag = normVec3d(unitGradRho, grad_rho[icv]);
//			  if(gradRhoMag < 1.0e-12)  // Sometimes unitGradRho becomes NaN when gradRhoMag==0.0
//				  for(int i=0; i<3; ++i)
//					  unitGradRho[i] = 0.0;
//
//			  int foc_f = faocv_i[icv];
//			  int foc_l = faocv_i[icv+1]-1;
//			  for (int foc=foc_f; foc<=foc_l; foc++) {
//				  int ifa = faocv_v[foc];
//
//				  double nVec[3];
//				  double area = normVec3d(nVec, fa_normal[ifa]);
//				  double dx = cv_volume[icv] / area; // Approximated grid-spacing in the face-normal direction: works only for a hex mesh
//
//				  double dxProjected = fabs(dx * vecDotVec3d(unitGradRho, nVec));
//
//				  if( dxProjected<0.0 || isnan(dxProjected) ) {
//					  cerr<<"UgpWithCvCompFlow::calcMaterialProperties(): wrong dxProjected is calculated = "<<dxProjected<<endl;
//					  throw(-1);
//				  }
//
////					  delta = max(delta, dxProjected);  // MAX_BASED
//				  delta += dxProjected;  // SUM_BASED
//			  }
//
//			  delta /= double(foc_l - foc_f + 1);
//		  } else { // If grad_rho is not defined, then you must find another (but inaccurate) method
			  double unitRhou[3];
			  double rhouMag = normVec3d(unitRhou, rhou[icv]);
			  assert(!isnan(rhouMag));
			  int countFaces = 0;

			  int foc_f = faocv_i[icv];
			  int foc_l = faocv_i[icv+1]-1;
			  for (int foc=foc_f; foc<=foc_l; foc++) {
				  int ifa = faocv_v[foc];

				  double nVec[3];
				  double area = normVec3d(nVec, fa_normal[ifa]);

				  double flowdirectInnerprod = fabs(vecDotVec3d(unitRhou, nVec));

				  if(fabs(rhouMag) < 1.0e-8 || flowdirectInnerprod > 2.0e-3) { // Note: cos(0.1 degrees) = 0.0017
					  double dx = cv_volume[icv] / area; // Approximated grid-spacing in the face-normal direction: works only for a hex mesh
					  delta += dx;  // Note: MAX_BASE cannot work since it will only take the span-wise grid-spacing in most of the 2D calculations
					  ++countFaces;
				  }
			  }

			  delta /= double(countFaces);
//		  }

		  // Calculate CV-based artificial viscosity
		  if(artifVisc_type.compare("STRAIN_BASE") == 0)
			  artifVisc_mag[icv] = artifVisc_coeff * rho[icv] * pow(delta, 2.0) * fabs(strMag[icv]) * clipFunc;
		  else if(artifVisc_type.compare("DIVERG_BASE") == 0)
			  artifVisc_mag[icv] = artifVisc_coeff * rho[icv] * pow(delta, 2.0) * fabs(diverg[icv]) * clipFunc;
		  else {
			  if(mpi_rank==0) cerr<<"UgpWithCvCompFlow::calcArtifVisc(): artifVisc_type is not supported = "<<artifVisc_type<<endl;
			  throw(-1);
		  }

		  // Check possible errors
		  if( artifVisc_mag[icv] < 0.0 || isnan(artifVisc_mag[icv]) ) {
			  cerr<<"UgpWithCvCompFlow::calcArtifVisc(): artifVisc_mag becomes negative = "<<artifVisc_mag[icv]<<endl;
			  throw(-1);
		  }

		  // Statistics
		  myArtifViscMagMin = min(myArtifViscMagMin, artifVisc_mag[icv]);
		  myArtifViscMagMax = max(myArtifViscMagMax, artifVisc_mag[icv]);
	  }
	  MPI_Allreduce(&myArtifViscMagMin, &artifViscMagMin, 1, MPI_DOUBLE, MPI_MIN, mpi_comm);
	  MPI_Allreduce(&myArtifViscMagMax, &artifViscMagMax, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);

	  // Interpolate artifVisc_mag on the faces
	  double myMaxRatioArtifToReal = 0.0;
	  if(!artifVisc_bulkViscOnly) { // If the user DOES NOT want to add the artificial viscosity only on the bulk viscosity
		  // internal faces
		  for (int ifa = nfa_b; ifa < nfa; ifa++) {
			  int icv0 = cvofa[ifa][0];
			  int icv1 = cvofa[ifa][1];

			  double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
			  vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
			  vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
			  double w0 = sqrt(vecDotVec3d(dx0, dx0));
			  double w1 = sqrt(vecDotVec3d(dx1, dx1));

			  if(mu_ref == 0)
				  mul_fa[ifa] = 0.0;  // Make sure that the baseline viscosity is equal to zero for inviscid calculations

			  double artifVisc_fa = (w1*artifVisc_mag[icv0] + w0*artifVisc_mag[icv1]) / (w0+w1);
			  myMaxRatioArtifToReal = max( myMaxRatioArtifToReal, calcRatioArtifAndRealVisc(mul_fa[ifa], artifVisc_fa, ifa) );
			  mul_fa[ifa] += artifVisc_fa;
		  }

		  // boundary faces computed next in setBC
	  } else { // Just calculate the ratio between the aritifical bulk visc. to molecular+turb visc.
		  // internal faces
		  for (int ifa = nfa_b; ifa < nfa; ifa++) {
			  int icv0 = cvofa[ifa][0];
			  int icv1 = cvofa[ifa][1];

			  double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
			  vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
			  vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
			  double w0 = sqrt(vecDotVec3d(dx0, dx0));
			  double w1 = sqrt(vecDotVec3d(dx1, dx1));

			  double artifVisc_fa = (w1*artifVisc_mag[icv0] + w0*artifVisc_mag[icv1]) / (w0+w1);
			  myMaxRatioArtifToReal = max( myMaxRatioArtifToReal, calcRatioArtifAndRealVisc(mul_fa[ifa], artifVisc_fa, ifa) );
		  }
	  }
	  MPI_Allreduce(&myMaxRatioArtifToReal, &maxRatioArtifToReal, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);

	  if(artifViscMagMin < 0.0) {
		  if(mpi_rank == 0)
			  cout<<"ERROR calcArtifVisc(): artifViscMagMin = "<<artifViscMagMin<<" is less than zero"<<endl;
		  throw(-1);
	  }

	  firstCall = false;
  }

  // \ Calculation of the ratio of the local artificial viscosity to the local viscosity (molecular + turbulent).
  virtual double calcRatioArtifAndRealVisc(const double local_mul_fa, const double local_artif_visc, const int ifa) {
	  double local_muT = 0.0;
	  return local_artif_visc / (local_mul_fa + local_muT + 1.0e-15);
  }
#endif
  
  // \brief Compute for a given temperature, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_T(double &p, double &h, double &R, double &gam, double &rho, double &T, double *Scal, int nScal)
  {
    R = R_gas;
    gam = GAMMA;
    p = rho * R * T;
    h = gam * R / (gam - 1.0) * T;
  }
 
  // \brief Compute for a given pressure, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_p(double &T, double &h, double &R, double &gam, double &rho, double &p, double *Scal, int nScal)
  {
    R = R_gas;
    gam = GAMMA;
    T = p / (rho * R);
    h = gam * R / (gam - 1.0) * T;
  }

  // \brief Compute for a given temperature the properties of the mixture at a face.
  virtual void ComputeBCProperties_T(FaZone *zone)
  {
	  //  for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
	  for (int index = 0; index < zone->faVec.size(); ++index) {
		  int ifa = zone->faVec[index];
		  int icv1 = cvofa[ifa][1];
		  gamma[icv1] = GAMMA;
		  RoM[icv1] = R_gas;
		  enthalpy[icv1] = GAMMA*R_gas/(GAMMA-1.0)*temp[icv1]; //TODO: why is enthalpy here?
	  }

	  if (mu_ref > 0.0)
	  {
		  if (viscMode == "SUTHERLAND")
		  {
			  for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
			  {
				  //for (int index = 0; index < zone->faVec.size(); ++index) {
				  // int ifa = zone->faVec[index];
				  int icv1 = cvofa[ifa][1];
				  mul_fa[ifa] = mu_ref*pow(temp[icv1]/SL_Tref, 1.5)*(SL_Tref+SL_Sref)/(temp[icv1]+SL_Sref);
			  }
		  }
		  else if (viscMode == "POWERLAW")
		  {
			  for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
			  {
				  //for (int index = 0; index < zone->faVec.size(); ++index) {
				  //  int ifa = zone->faVec[index];
				  int icv1 = cvofa[ifa][1];
				  mul_fa[ifa] = mu_ref*pow(temp[icv1]/T_ref, mu_power_law);
			  }
		  }
		  else
		  {
			  cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
			  throw(-1);
		  }

		  for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
		  {
			  //for (int index = 0; index < zone->faVec.size(); ++index) {
			  //  int ifa = zone->faVec[index];
			  lamOcp_fa[ifa] = mul_fa[ifa]/Pr;
		  }
	  }
  }  
  
  /*
   * Method: ComputeBCProperties1D_T
   * -------------------------------
   * brief Compute for a given temperature the properties of the mixture at a face.
   * Original code: ComputeBCProperties1D_T_AD() in IkeWithModel.cpp and ComputeBCProperties_T() UgpWithCvCompFlow.h
   *                In an old version, this method used to be in IkeUgpWithCvCompFlow.h
   *                but moved to here so that CombModel_FPVA_Coeff.h can recognize this method.
   * In some situation, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
   */
  virtual void ComputeBCProperties1D_T(const int ifa) {
	  double *gamma = UgpWithCvCompFlow::gamma;
	  double *RoM = UgpWithCvCompFlow::RoM;
	  double *enthalpy = UgpWithCvCompFlow::enthalpy;
	  double *temp = UgpWithCvCompFlow::temp;
	  double *mul_fa = UgpWithCvCompFlow::mul_fa;
	  double *lamOcp_fa = UgpWithCvCompFlow::lamOcp_fa;

	  int icv1=cvofa[ifa][1];
	  gamma[icv1] = GAMMA;
	  RoM[icv1] = R_gas;
	  enthalpy[icv1] = GAMMA * R_gas / (GAMMA - 1.0) * temp[icv1];

	  if (mu_ref > 0.0) {
		  if (viscMode == "SUTHERLAND") {
			  int icv1=cvofa[ifa][1];
			  mul_fa[ifa] = mu_ref*pow(temp[icv1]/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(temp[icv1] + SL_Sref);
		  } else if (viscMode == "POWERLAW") {
			  int icv1=cvofa[ifa][1];
			  mul_fa[ifa] = mu_ref*pow(temp[icv1]/T_ref, mu_power_law);
		  } else {
			  cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
			  throw(-1);
		  }

		  lamOcp_fa[ifa] = mul_fa[ifa] / Pr;
	  }
  }

  // \brief Compute for a given enthalpy the properties of the mixture at a face.
  virtual void ComputeBCProperties_H(FaZone *zone)
  {
	  //for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
	  for (int index = 0; index < zone->faVec.size(); ++index) {
		  int ifa = zone->faVec[index];
		  int icv1 = cvofa[ifa][1];
		  gamma[icv1] = GAMMA;
		  RoM[icv1] = R_gas;
		  temp[icv1] = enthalpy[icv1]*(GAMMA-1.0)/(GAMMA*R_gas);
	  }

	  if (mu_ref > 0.0)
	  {
		  if (viscMode == "SUTHERLAND")
		  {
			  for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
			  {
				  //for (int index = 0; index < zone->faVec.size(); ++index) {
				  //int ifa = zone->faVec[index];
				  int icv1 = cvofa[ifa][1];
				  mul_fa[ifa] = mu_ref*pow(temp[icv1]/SL_Tref, 1.5)*(SL_Tref+SL_Sref)/(temp[icv1]+SL_Sref);
			  }
		  }
		  else if (viscMode == "POWERLAW")
		  {
			  for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
			  {
				  //for (int index = 0; index < zone->faVec.size(); ++index) {
				  //  int ifa = zone->faVec[index];
				  int icv1 = cvofa[ifa][1];
				  mul_fa[ifa] = mu_ref*pow(temp[icv1]/T_ref, mu_power_law);
			  }
		  }
		  else
		  {
			  cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
			  throw(-1);
		  }
		  for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
		  {
			  //for (int index = 0; index < zone->faVec.size(); ++index) {
			  // int ifa = zone->faVec[index];
			  lamOcp_fa[ifa] = mul_fa[ifa] / Pr;
		  }
	  }
  }

  /*
   * Method: ComputeBCProperties1D_H
   * -------------------------------
   * brief Compute for a given enthalpy the properties of the mixture at a face.
   * Original code: ComputeBCProperties1D_H_AD() in UgpWithCvCompFlowAD.h
   *                In an old version, this method used to be in IkeUgpWithCvCompFlow.h
   *                but moved to here so that CombModel_FPVA_Coeff.h can recognize this method.
   */
  virtual void ComputeBCProperties1D_H(const int ifa) {
	  double *gamma = UgpWithCvCompFlow::gamma;
	  double *RoM = UgpWithCvCompFlow::RoM;
	  double *enthalpy = UgpWithCvCompFlow::enthalpy;
	  double *temp = UgpWithCvCompFlow::temp;
	  double *mul_fa = UgpWithCvCompFlow::mul_fa;
	  double *lamOcp_fa = UgpWithCvCompFlow::lamOcp_fa;

	  int icv1 = cvofa[ifa][1];
	  gamma[icv1] = GAMMA;
	  RoM[icv1] = R_gas;
	  temp[icv1] = enthalpy[icv1]*(GAMMA-1.0)/(GAMMA*R_gas);

	  if (mu_ref > 0.0) {
		  if (viscMode == "SUTHERLAND") {
			  int icv1 = cvofa[ifa][1];
			  mul_fa[ifa] = mu_ref*pow(temp[icv1]/SL_Tref, 1.5)*(SL_Tref+SL_Sref)/(temp[icv1]+SL_Sref);
		  } else if (viscMode == "POWERLAW") {
			  int icv1 = cvofa[ifa][1];
			  mul_fa[ifa] = mu_ref*pow(temp[icv1]/T_ref, mu_power_law);
		  } else {
			  cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
			  throw(-1);
		  }

		  lamOcp_fa[ifa] = mul_fa[ifa] / Pr;
	  }
  }



  // ###########################################################################
  // ---------------------------------------------------------------------------
  //
  // combustion combustion combustion combustion combustion combustion
  //
  // ---------------------------------------------------------------------------  
  // ###########################################################################

  virtual void diffusivityHookScalarRansComb(const string &name)  {/*empty*/}
  virtual void sourceHookScalarRansCombExplicit(double *rhs, const string &name)  {/*empty*/}
  virtual void sourceHookScalarRansComb(double *rhs, double *A, const string &name, int flagImplicit)  {/*empty*/}
///IKJ
  virtual void sourceHookInSourceHookScalarRansComb(double *rhs, double *A, const string &name, int flagImplicit)  {/*empty*/}
///IKJ
  virtual void sourceHookScalarRansComb_ikj(double **rhs_rhoScal, double ***Ascal, vector<ScalarTranspEq> &scalarTranspEqVector, int flagImplicit) {/*empty*/}

  virtual void boundaryHookScalarRansComb(double *phi_ph, FaZone *zone, const string &name)  {/*empty*/}  
  virtual void sourceHookRansComb(double *rhsRho, double(*rhsRhou)[3], double *rhsRhoE, double(*A)[5][5])  {/*empty*/}  

  virtual void sourceHookRansCombCoupled(double **rhs, double ***A, int nScal, int flagImplicit)  {/*empty*/}

  virtual void pressureDerivativeHookScalarRansComb() {/*empty*/}


  
  
  
  
  
  
  
  

  int zoneIsWall(const string & faName) 
  {
    Param *p;
    string paraName;
    int retVal = 0;

    if (getParam(p, faName)) {
      paraName = p->getString();
      retVal = (paraName == "WALL") || (paraName == "Wall") || (paraName == "wall") ||
               (paraName == "WALL.T") || (paraName == "Wall.T") || (paraName == "wall.T");
    }
    return retVal;
  }

  // #########################################################################################
  //
  //    scalar boundaries / scalar boundaries / scalar boundaries
  //
  // #########################################################################################
  int scalarZoneIsHook(const string &faName, const string &scalName) {

    Param *p;
    int retVal = 0;

    if (getParam(p, faName + "." + scalName)) {
      string name = p->getString(1);
      retVal = (name == "HOOK") || (name == "Hook")|| (name == "hook");
    } else
      retVal = 0;

    return retVal;
  }

  int scalarZoneIsDirichlet(double &phiBC, const string &faName, const string &scalName) {

    Param *p;
    int retVal = 0;

    if (getParam(p, faName + "." + scalName)) {
      phiBC = p->getDouble(1);
      retVal = 1;
    } else
      retVal = 0;

    return retVal;
  }

  int scalarZoneIsFlux(double phiBCflux, const string &faName, const string &scalName) {

    Param *p;
    string paraName;
    int retVal = 0;

    if (getParam(p, faName + "." + scalName + ".flux")) {
      phiBCflux = p->getDouble(1);
      retVal = 1;
    } else
      retVal = 0;

    return retVal;
  }


// Utilities developed for IKE but also useful in Joe
#include "../JoeUgpWithCvCompFlow.h"

};


#endif
