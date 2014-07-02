#ifndef UGPWITHCVCOMPFLOWAD_H
#define UGPWITHCVCOMPFLOWAD_H

#define USE_MEM_SAVING_ADVAR
#ifdef USE_MEM_SAVING_ADVAR
	#include "../../ADvar.h"
//	#define USE_CONDASSIGN
#endif

//#define ADOLC_TAPELESS // ADOL-C tapeless mode (faster calculation, less flexibility)
#ifdef ADOLC_TAPELESS
	#define NUMBER_DIRECTIONS 10
#endif

#include "MiscUtils.h"
#include "tc_vec3d.h"

#include "UgpWithCvCompFlow.h"
#include "adolc.h"

#ifdef ADOLC_TAPELESS
	typedef adtl::adouble adouble;
	//#define adouble adtl::adouble
#endif

#ifndef MACHINE_EPS
#define MACHINE_EPS 2.0e-16 // Note: machine epsilon for double precision = 1.11e-16
#endif

#define REALX   double
#define REALQ  adouble
#define REALQS adouble
#define REALA  adouble
#define REALAS adouble

/**
 * Class contains the main routines for joe adjoints
 *
 * \author Frank Ham, Rene Pecnik, Karthik Duraisamy
 * \date June, 2010
 */

enum BOUNDARY_TYPE {WRONG, INTERNAL, HOOK, CBC, CBC_SUBSONIC_INLET, CBC_SUBSONIC_OUTLET, SYMMETRY, NEUMANN, WALL};
enum BOUNDARY_TYPE_SCALAR {INTERNAL_SCALAR, HOOK_SCALAR, DIRICHLET_SCALAR, FLUX_SCALAR, OTHER_SCALAR};

class UgpWithCvCompFlow_AD : virtual public UgpWithCvCompFlow
{
public:   // constructors/destructors
	bool turnOnAD;

  /*
   * constructor, pass ParamMap
   */
   UgpWithCvCompFlow_AD(ParamMap &p) : UgpWithCvCompFlow(p) {    init();  }
  /*
   * constructor, pass name of joe's input file
   */
   UgpWithCvCompFlow_AD(char *name) : UgpWithCvCompFlow(name) {    init();  }

   UgpWithCvCompFlow_AD() {init();}


  /**
   * registers scalars and vectors for main equations (Navier-Stokes)
   * and set main parameters (solver: nsteps; NS-eq: gravity...)
   */
  void init()
  {
    if (mpi_rank == 0)    cout << "UgpWithCvCompFlow_AD()"<< endl;

    scalarTranspEqVector_AD = NULL;

    // ----------------------------------------------------------------------------------------
    // register data...
    // ----------------------------------------------------------------------------------------

    // adjoint variables
    psi_rho   = NULL;         registerScalar(psi_rho,  "PSI_RHO",  CV_DATA);
    psi_rhou  = NULL;         registerVector(psi_rhou, "PSI_RHOU", CV_DATA);
    psi_rhoE  = NULL;         registerScalar(psi_rhoE, "PSI_RHOE", CV_DATA);

    if (mpi_rank == 0)
    {
      cout << endl << "--------------------------------------------" << endl;
      cout<<"Adjoint Solver initialized"<<endl;
      cout << "--------------------------------------------" << endl << endl;
    }

    UgpWithCvCompFlow_AD_init();
  }

  virtual ~UgpWithCvCompFlow_AD()
  {
	  if(mpi_rank==0)
		  cout<<"~UgpWithCvCompFlow_AD()"<<endl;
	  UgpWithCvCompFlow_AD_clear();
  }

public:   // member variables

   double *adj_vars;   ///< Adjoint variable vector.  This is redundant. Fix later.

   double *psi_rho;         ///< Adjoint density
   double (*psi_rhou)[3];   ///< Adjoint momentum
   double *psi_rhoE;        ///< Adjoint energy

   double *psi_nuSA;        ///< Adjoint nuSA
   double *psi_kine;        ///< Adjoint TKE
   double *psi_omega;       ///< Adjoint  omega
   double *psi_ZMean;       ///< Adjoint Mixture fraction
   double *psi_ZVar;        ///< Adjoint Mixture fraction variance
   double *psi_CMean;       ///< Adjoint Progress variable
   double *psi_phi;         ///< Adjoint transfer (TBD_AD: Get rid of this soon - right now this is in for transitioning )


  // Other variables will be needed when we go for turbulence and combustion

  adouble functional;        ///<   J

  double *functional_der;      ///< dJ/dU
  double *residual_der;        ///< (dR/dU)^T V

#ifndef USE_MEM_SAVING_ADVAR
  adouble  *rho_AD;       ///< density derivative
  adouble (*rhou_AD)[3];  ///< momentum derivative
  adouble  *rhoE_AD;      ///< energy derivative

  adouble (*vel)[3];   ///< velocity
  adouble *press;      ///< pressure e
  adouble *temp;       ///< temperature
  adouble *enthalpy;   ///< chemical + sensible enthalpy (no kinetic energy!) derivative

  adouble *gamma;      ///< gamma derivative
  adouble *RoM;        ///< Gas constant derivative
  adouble *sos;        ///< speed of sound derivative

  adouble *mul_fa;           ///< laminar viscosity at cell faces
  adouble *lamOcp_fa;        ///< heat conductivity at cell faces

  REALQS *mut_fa;     ///< turbulent viscosity derivative
  REALQS *kine;      ///< turbulent kinetic energy derivative

  adouble *strMag, *vortMag, *diverg;  ///< magnitude strain rate, divergence derivatives
  adouble *blendFuncF1_AD, *blendFuncF2_AD, *crossDiff_AD; ///<SST model functions

  adouble *CmeanSource_AD,*muLam_AD,*LambdaOverCp_AD,*chi_AD; ///<FPVA model functions

  // ----------------------------------------------------------------------------------------
  // boundary faces, allocate memory, TODO: should be incorporated in to the associated arrays ???
  // ----------------------------------------------------------------------------------------

  adouble (*grad_rho)[3];      ///< density  gradient der
  adouble (*grad_u)[3][3];     ///< velocity gradient tensor der
  adouble (*grad_p)[3];        ///< pressure gradient der
  adouble (*grad_temp)[3];     ///< temperature gradient der
  adouble (*grad_enthalpy)[3]; ///< enthalpy gradient der
#endif

  int nsteps_AD;                  ///< Adjoint solver time steps

#ifndef USE_MEM_SAVING_ADVAR
  class ScalarTranspEq_AD
{
  public:
	  ScalarTranspEq_AD() {
		  phi = NULL;
		  diff = NULL;
		  grad_phi = NULL;
		  rhophi = NULL;
		  grad_rhophi = NULL;
		  dpress_dphi = NULL;
	  }
	  ~ScalarTranspEq_AD() {
		  if(phi != NULL) {
			  delete [] phi; 			phi = NULL;
		  }
		  if(diff != NULL) {
			  delete [] diff; 			diff = NULL;
		  }
		  if(grad_phi != NULL) {
			  delete [] grad_phi; 		grad_phi = NULL;
		  }
		  if(rhophi != NULL) {
			  delete [] rhophi; 		rhophi = NULL;
		  }
		  if(grad_rhophi != NULL) {
			  delete [] grad_rhophi; 	grad_rhophi = NULL;
		  }
		  if(dpress_dphi != NULL) {
			  delete [] dpress_dphi; 	dpress_dphi = NULL;
		  }
	  }

  adouble *phi;
  adouble *diff;
  adouble (*grad_phi)[3];
  adouble *rhophi;
  adouble (*grad_rhophi)[3];
  adouble *dpress_dphi;
  char name[DATA_NAME_LEN];
};
#endif
  
#include "../../IkeUgpWithCvCompFlow.h"

    ScalarTranspEq_AD *scalarTranspEqVector_AD;

  class psi_AD
{
  public:

  double *phi;
};

  psi_AD  *scalarTranspEqVector_psi;

protected: 

public:   // member functions
   void initializeAdjoint()
   {
      int nScal = scalarTranspEqVector.size();

      adj_vars       = new double[(5+nScal)*ncv_ggff]; // TBD_AD We don't really need the fakes here
      functional_der = new double[(5+nScal)*ncv_ggff];     
      residual_der   = new double[(5+nScal)*ncv_ggff];     

      if(nScal > 0) {
    	  scalarTranspEqVector_psi = new psi_AD[nScal];

    	  for(int i=0; i<nScal ; i++)
    		  scalarTranspEqVector_psi[i].phi      = new double[ncv_ggff];
      }
   }

#ifndef USE_MEM_SAVING_ADVAR
   virtual void initialize_adoubles()
   {
    rho_AD        = new adouble [ncv_ggff];
    rhou_AD       = new adouble[ncv_ggff][3];
    rhoE_AD       = new adouble[ncv_ggff];

    vel        = new adouble[ncv_ggff][3];
    press      = new adouble[ncv_ggff];
    temp       = new adouble[ncv_ggff];
    enthalpy   = new adouble[ncv_ggff];

    gamma          = new adouble[ncv_ggff];
    RoM            = new adouble[ncv_ggff];
    sos            = new adouble[ncv_ggff];

    strMag      = new adouble[ncv_gg];
    vortMag     = new adouble[ncv_gg];
    diverg      = new adouble[ncv_gg];

    //TBD_AD MOVE TO KOMSST
   blendFuncF1_AD = new adouble[ncv_gg];
   blendFuncF2_AD = new adouble[ncv_gg];
   crossDiff_AD   = new adouble[ncv_gg];
    //TBD_AD MOVE TO FPVA
   CmeanSource_AD   = new adouble[ncv_gg];
   muLam_AD   	    = new adouble[ncv_gg];
   LambdaOverCp_AD  = new adouble[ncv_gg];
   chi_AD   	    = new adouble[ncv_gg];

   kine    = NULL;

    // ----------------------------------------------------------------------------------------
    // init memory for face-based data
    // ----------------------------------------------------------------------------------------

    mul_fa      = new adouble[nfa_b2gg];
    lamOcp_fa   = new adouble[nfa_b2gg];
    mut_fa      = new adouble[nfa_b2gg];


    // ----------------------------------------------------------------------------------------
    // gradients
    // ----------------------------------------------------------------------------------------
    if (sndOrder == true)
    {
      grad_rho  = new adouble[ncv_gg][3];
#ifdef temp_reconstruction
      grad_temp = new adouble[ncv_gg][3];
#else
      grad_p    = new adouble[ncv_gg][3];
#endif
    }

    grad_u = new adouble[ncv_gg][3][3];         // allocate always

    if ((sndOrder == true) || (mu_ref > 0.0))
      grad_enthalpy = new adouble[ncv_gg][3];

    int nScal = scalarTranspEqVector.size();

    if(nScal>0)
    	scalarTranspEqVector_AD = new ScalarTranspEq_AD[nScal];
    else
    	scalarTranspEqVector_AD = NULL;

    for(int i=0; i<nScal ; i++){
        scalarTranspEqVector_AD[i].phi      = new adouble[ncv_ggff];
        scalarTranspEqVector_AD[i].diff     = new adouble[nfa_b2gg];
        scalarTranspEqVector_AD[i].rhophi   = new adouble[ncv_ggff];
        scalarTranspEqVector_AD[i].grad_phi = new adouble[ncv_gg][3];
        scalarTranspEqVector_AD[i].grad_rhophi = new adouble[ncv_gg][3];
        scalarTranspEqVector_AD[i].dpress_dphi = NULL;
        strcpy(scalarTranspEqVector_AD[i].name, scalarTranspEqVector[i].getName());
//        if(mpi_rank==0) cout<<"Allocated Scalar "<<scalarTranspEqVector_AD[i].name<<endl;
    }
   }


   virtual void destroy_adoubles()
   {
    //TBD_AD DESTROY SCALAR ADOUBLES
    if(scalarTranspEqVector_AD!= NULL) delete [] scalarTranspEqVector_AD;

    if ((sndOrder == true) || (mu_ref > 0.0))
       delete [] grad_enthalpy ;
     delete [] grad_u  ;
    if (sndOrder == true)
    {
#ifdef temp_reconstruction
      delete [] grad_temp;
#else
      delete [] grad_p;
#endif
      delete [] grad_rho;
    }
     delete [] mut_fa ;
     delete [] lamOcp_fa ;
     delete [] mul_fa ;

     delete [] chi_AD ;
     delete [] LambdaOverCp_AD ;
     delete [] muLam_AD ;
     delete [] CmeanSource_AD ;

     delete [] crossDiff_AD ;
     delete [] blendFuncF2_AD ;
     delete [] blendFuncF1_AD ;

     delete [] diverg ;
     delete [] vortMag ;
     delete [] strMag ;

     delete [] sos      ;
     delete [] RoM      ;
     delete [] gamma    ;

     delete [] enthalpy ;
     delete [] temp     ;
     delete [] press    ;
     delete [] vel      ;

     delete []  rhoE_AD  ;
     delete []  rhou_AD  ;
     delete []  rho_AD   ;
   }
#endif

   virtual void initialHookScalarRansCombModel_AD(int loadtable)
  {
    for (int icv = 0; icv < ncv_ggff; icv++)
    {
      RoM[icv] = R_gas;
      gamma[icv] = GAMMA;
    }

  }

   virtual void initialHookScalarRansTurbModel_AD()
  {
  }

   virtual void calcStateVariables_AD(REALQ *rho, REALQ(*rhou)[3], REALQ *rhoE)
  {
    for (int icv = 0; icv < ncv_gg; icv++)
    {
      if (rho[icv] <= 0.0)
        cout << "Negative density at xcv: " << x_cv[icv][0] << ", " 
                                            << x_cv[icv][1] << ", " << x_cv[icv][2] << endl;

      for (int i=0; i<3; i++)
        vel[icv][i] = rhou[icv][i]/rho[icv];

      REALQS kinecv = 0.0;

#ifdef USE_MEM_SAVING_ADVAR
      if (kine != NULL)
        kinecv = (*kine)[icv];
#else
      if (kine != NULL)
        kinecv = kine[icv];
#endif

      REALQ pr = (gamma[icv]-1.0)*(rhoE[icv] 
      - 0.5*(rhou[icv][0]*rhou[icv][0]+rhou[icv][1]*rhou[icv][1]+rhou[icv][2]*rhou[icv][2])/rho[icv] 
      - rho[icv]*kinecv);
      if (pr <= 0.0)
      {
        cout << "Calc RANS negative pressure at xcv: " << x_cv[icv][0] << ", " << x_cv[icv][1] << ", " << x_cv[icv][2] << gamma[icv]<<endl;
      }
      else
        press[icv] = pr;

      temp[icv] = press[icv]/(rho[icv]*RoM[icv]);
      enthalpy[icv] = gamma[icv]*RoM[icv]/(gamma[icv]-1.0)*temp[icv];
      sos[icv] = sqrt(gamma[icv]*press[icv]/rho[icv]);
    }

  }

   virtual void calcMaterialProperties_AD(REALQ *rho, REALQ(*rhou)[3], REALQ *rhoE)
  {

        if (mu_ref > 0.0) 
    {    
      if (viscMode == "SUTHERLAND")
      {    
        // internal faces
       //for (int ifa = nfa_b; ifa < nfa; ifa++)
  	for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) // Ensure we sweep through internal faces of ghosts too
  	{
   	if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg))
   	{

          int icv0 = cvofa[ifa][0];
          int icv1 = cvofa[ifa][1];

          REALX dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
          vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
          vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
          REALX w0 = sqrt(vecDotVec3d(dx0, dx0));
          REALX w1 = sqrt(vecDotVec3d(dx1, dx1));

          REALQ temperature  = (w1*temp[icv0] + w0*temp[icv1])/(w0+w1);
          mul_fa[ifa] = mu_ref*pow(temperature/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(temperature + SL_Sref);
          lamOcp_fa[ifa] = mul_fa[ifa] / Pr;
        }    
        }    
          // boundary faces computed next in setBC
      }    
      else if (viscMode == "POWERLAW")
      {    
        // internal faces
         //for (int ifa = nfa_b; ifa < nfa; ifa++)
  for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) // Ensure we sweep through internal faces of ghosts too
  {
   if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg))
   {

          int icv0 = cvofa[ifa][0];
          int icv1 = cvofa[ifa][1];

          REALX dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
          vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
          vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
          REALX w0 = sqrt(vecDotVec3d(dx0, dx0));
          REALX w1 = sqrt(vecDotVec3d(dx1, dx1));
   
          REALQ temperature  = (w1*temp[icv0] + w0*temp[icv1])/(w0+w1);
          mul_fa[ifa] = mu_ref*pow(temperature/T_ref, mu_power_law);
          lamOcp_fa[ifa] = mul_fa[ifa] / Pr;
        }    
   } 
        // boundary faces computed next in setBC
      }    
      else 
      {    
        cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
        throw(-1);
      }    
    }    


  }












   // \brief Compute for a given temperature, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_T_AD(REALQ &p, REALQ &h, REALQ &R, REALQ &gam, REALQ &rho, REALQ &T, REALQS *Scal, int nScal)
  {
    R = R_gas;
    gam = GAMMA;
    p = rho * R * T;
    h = gam * R / (gam - 1.0) * T;
  }

  // \brief Compute for a given pressure, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_p_AD(REALQ &T, REALQ &h, REALQ &R, REALQ &gam, REALQ &rho, REALQ &p, REALQS *Scal, int nScal)
  {
    R = R_gas;
    gam = GAMMA;
    T = p / (rho * R);
    h = gam * R / (gam - 1.0) * T;
  }


 // \brief Compute for a given temperature the properties of the mixture at a face.
  virtual void ComputeBCProperties_T_AD(FaZone *zone)
  {
    //for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
    //{
      for (int index = 0; index < zone->faVec.size(); ++index) {
      int ifa = zone->faVec[index];
      int icv1=cvofa[ifa][1];
      gamma[icv1] = GAMMA;
      RoM[icv1] = R_gas;
      enthalpy[icv1] = GAMMA * R_gas / (GAMMA - 1.0) * temp[icv1];
    }

    if (mu_ref > 0.0)
    {
      if (viscMode == "SUTHERLAND")
        //for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++){
        for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
          int icv1=cvofa[ifa][1];
          mul_fa[ifa] = mu_ref*pow(temp[icv1]/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(temp[icv1] + SL_Sref);
        }
      else if (viscMode == "POWERLAW")
        //for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++){
	  for (int index = 0; index < zone->faVec.size(); ++index) {
      	  int ifa = zone->faVec[index];
          int icv1=cvofa[ifa][1];
          mul_fa[ifa] = mu_ref*pow(temp[icv1]/T_ref, mu_power_law);
        }
      else
      {
        cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
        throw(-1);
      }

      //for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      for (int index = 0; index < zone->faVec.size(); ++index) {
      int ifa = zone->faVec[index];
        lamOcp_fa[ifa] = mul_fa[ifa] / Pr;
      }
    }
  }


  // \brief Compute for a given enthalpy the properties of the mixture at a face.
  virtual void ComputeBCProperties_H_AD(FaZone *zone)
  {
    //for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
    //{
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
       // for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
        //{
	  for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
          int icv1 = cvofa[ifa][1];
          mul_fa[ifa] = mu_ref*pow(temp[icv1]/SL_Tref, 1.5)*(SL_Tref+SL_Sref)/(temp[icv1]+SL_Sref);
        }
      }
      else if (viscMode == "POWERLAW")
      { 
        //for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
        //{ 
	for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
          int icv1 = cvofa[ifa][1];
          mul_fa[ifa] = mu_ref*pow(temp[icv1]/T_ref, mu_power_law);
        }
      }
      else
      { 
        cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
        throw(-1); 
      }   
      //for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
	for (int index = 0; index < zone->faVec.size(); ++index) {
        int ifa = zone->faVec[index];
        lamOcp_fa[ifa] = mul_fa[ifa] / Pr;
	}
    } 
  }   

    virtual adouble calcMuLam_AD(adouble temp)
  {
    adouble muLam = 0.0;

    if (viscMode == "SUTHERLAND")      muLam = mu_ref*pow(temp/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(temp + SL_Sref);
    else if (viscMode == "POWERLAW")   muLam = mu_ref*pow(temp/T_ref, mu_power_law);
    else muLam = 0.0;

    return muLam;
  }

     virtual adouble calcMuLam_AD(int icv)
  {
    adouble muLam = 0.0;

    if (viscMode == "SUTHERLAND")      muLam = mu_ref*pow(temp[icv]/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(temp[icv] + SL_Sref);
    else if (viscMode == "POWERLAW")   muLam = mu_ref*pow(temp[icv]/T_ref, mu_power_law);
    else muLam = 0.0;

    return muLam;
  }

   virtual void calcRansTurbViscMuet_AD(REALQ *rho, REALQ (*rhou)[3])
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
   


   /**
   *
   *  explicit euler flux HLLC
   *
   */
  int calcEulerFlux_HLLC_new_AD(REALQ &Frho, REALQ *Frhou, REALQ &FrhoE, REALQS *FrhoScal,
      REALQ rhoL, REALQ *uL, REALQ pL, REALQ TL, REALQ h0, REALQ RL, REALQ gammaL, REALQS *ScalL, REALQS kL,
      REALQ rhoR, REALQ *uR, REALQ pR, REALQ TR, REALQ h1, REALQ RR, REALQ gammaR, REALQS *ScalR, REALQS kR,
      REALX area, REALX *nVec, int nScal,   double surfVeloc);

  int calcEulerFlux_HLLC_AD(REALQ &Frho, REALQ *Frhou, REALQ &FrhoE, REALQS *FrhoScal,
      REALQ rhoL, REALQ *uL, REALQ pL, REALQ TL, REALQ h0, REALQ RL, REALQ gammaL, REALQS *ScalL, REALQS kL,
      REALQ rhoR, REALQ *uR, REALQ pR, REALQ TR, REALQ h1, REALQ RR, REALQ gammaR, REALQS *ScalR, REALQS kR,
      REALX area, REALX *nVec, int nScal,   double surfVeloc);

  int calcEulerFlux_HLLC_orig_AD(REALQ &Frho, REALQ *Frhou, REALQ &FrhoE, REALQS *FrhoScal,
      REALQ rhoL, REALQ *uL, REALQ pL, REALQ TL, REALQ h0, REALQ RL, REALQ gammaL, REALQS *ScalL, REALQS kL,
      REALQ rhoR, REALQ *uR, REALQ pR, REALQ TR, REALQ h1, REALQ RR, REALQ gammaR, REALQS *ScalR, REALQS kR,
      REALX area, REALX *nVec, int nScal,   double surfVeloc);

  int calcEulerFluxMatrices_HLLC_AD(double (*A_L)[5], double (*A_R)[5], double (*A_L_Scal)[6], double (*A_R_Scal)[6],
        REALQ rhoL_AD, REALQ *uL_AD, REALQ pL_AD, REALQ TL_AD, REALQ h0_AD, REALQ RL_AD, REALQ gammaL_AD, REALQS *scalL, REALQS kL,
        REALQ rhoR_AD, REALQ *uR_AD, REALQ pR_AD, REALQ TR_AD, REALQ h1_AD, REALQ RR_AD, REALQ gammaR_AD, REALQS *scalR, REALQS kR,
        REALX area, REALX *nVec, int nScal, double surfVeloc);

  /*
   * Explicit viscous fluxes
   */

   /**
   *
   * viscous flux routine for implicit
   * needs lsg_coeff0 and lsg_coeff1 ---> LS gradient coefficients
   *
   */
  virtual void addViscFlux_AD(REALQ *Frhou, REALQ &FrhoE, double (*A0)[5], double (*A1)[5],
      REALQ rho0, REALQ *u0, REALQ (&grad_u0)[3][3], REALQ h0, REALQ *grad_h0, REALQ T0, REALQ R0, REALQ gam0, REALQS kine0,
      REALQ rho1, REALQ *u1, REALQ (&grad_u1)[3][3], REALQ h1, REALQ *grad_h1, REALQ T1, REALQ R1, REALQ gam1, REALQS kine1,
      REALQ mul, REALQS mut, REALQ lambdaOverCp, REALQS kine_fa, REALQ *u_fa,
      REALX area, REALX *nVec, REALX smag, REALX *sVec);

   virtual void addViscFluxJacobians(double (*A0)[5], double (*A1)[5],
    REALQ rho0_AD, REALQ *u0_AD, REALQ (&grad_u0_AD)[3][3], REALQ h0_AD, REALQ *grad_h0_AD, REALQ T0_AD, REALQ R0_AD, REALQ gam0_AD, REALQS kine0_AD,
    REALQ rho1_AD, REALQ *u1_AD, REALQ (&grad_u1_AD)[3][3], REALQ h1_AD, REALQ *grad_h1_AD, REALQ T1_AD, REALQ R1_AD, REALQ gam1_AD, REALQS kine1_AD,
    REALQ mul_AD,  REALQS mut_AD, REALQ lambdaOverCp_AD, REALQS kine_fa_AD, REALQ *u_fa_AD,
    REALX area, REALX *nVec, REALX smag, REALX *sVec);

    void setScalarBC_AD(FaZone *zone);

    void calcViscousFluxScalar_new_AD(adouble *rhs_rhoScal, double *Ascal, ScalarTranspEq &transpScal, ScalarTranspEq_AD &transpScal_AD, int flagImplicit);

    virtual void boundaryHookScalarRansTurb_AD(adouble *phi_ph, FaZone *zone, const string &name)  {/*empty*/}
    virtual void diffusivityHookScalarRansTurb_AD(const string &name)  {/*empty*/}
    virtual void sourceHookScalarRansTurb_new_AD(adouble *rhs, double *A, const string &name, int flagImplicit) {/*empty*/}

    virtual void boundaryHookScalarRansComb_AD(adouble *phi_ph, FaZone *zone, const string &name)  {/*empty*/}
    virtual void diffusivityHookScalarRansComb_AD(const string &name)  {/*empty*/}
    virtual void sourceHookScalarRansComb_new_AD(adouble *rhs, double *A, const string &name, int flagImplicit)  {/*empty*/}

    virtual void sourceHookCoupled_AD(REALA **rhs,          double ***A,  int flagImplicit) {/* empty */}
    virtual void sourceHookRansTurbCoupled_AD(REALA **rhs,  double ***A,  int flagImplicit) {/* empty */}
    virtual void sourceHookRansCombCoupled_AD(REALA **rhs,  double ***A,  int flagImplicit) {/* empty */}
  
    virtual void pressureDerivativeHookScalarRansTurb_AD() {/*empty */}
    virtual void pressureDerivativeHookScalarRansComb_AD() {/*empty */}

    virtual void initialize_turb_adjoint() {/* empty */}
    virtual void copy_turb_adjoint() {/* empty */}
    virtual void initialize_comb_adjoint() {/* empty */}
    virtual void copy_comb_adjoint() {/* empty */}


    void calcEulerFluxCoupled_HLLC_AD(adouble *EulerFlux,
         adouble rhoL, adouble *uL, adouble pL, adouble TL, adouble h0, adouble RL, adouble gammaL, adouble *ScalL, adouble kL,
         adouble rhoR, adouble *uR, adouble pR, adouble TR, adouble h1, adouble RR, adouble gammaR, adouble *ScalR, adouble kR,
         double area, double *nVec, const int nScal, double *ConvTerm, adouble surfVeloc, const string BC_treatment);

    int calcEulerFluxMatricesCoupled_HLLC_AD(double **A_L, double **A_R,
        REALQ rhoL_AD, REALQ *uL_AD, REALQ pL_AD, REALQ TL_AD, REALQ h0_AD, REALQ RL_AD, REALQ gammaL_AD, REALQS *scalL_AD, REALQ *dpress_dscalL_AD, REALQS kL_AD,
        REALQ rhoR_AD, REALQ *uR_AD, REALQ pR_AD, REALQ TR_AD, REALQ h1_AD, REALQ RR_AD, REALQ gammaR_AD, REALQS *scalR_AD, REALQ *dpress_dscalR_AD, REALQS kR_AD,
        double area, double *nVec, int nScal, double *ConvTerm, double surfVeloc, const string BC_treatment);

    void calcViscousFluxCoupled_AD(REALQ *ViscousFlux, double **A0, double **A1,
         REALQ rho0, REALQ *u0, REALQ (&grad_u0)[3][3], REALQ h0, REALQ *grad_h0, REALQ T0, REALQ R0, REALQ gam0, REALQS *Scal0, REALQS (*gradScal0)[3], REALQ *dpress_dscal0, REALQS kine0,
         REALQ rho1, REALQ *u1, REALQ (&grad_u1)[3][3], REALQ h1, REALQ *grad_h1, REALQ T1, REALQ R1, REALQ gam1, REALQS *Scal1, REALQS (*gradScal1)[3], REALQ *dpress_dscal1, REALQS kine1,
         REALQ  mul, REALQ mut, REALQ lambdaOverCp, REALQS kine_fa, REALQ *u_fa, REALQS *diff, double *DiffTerm,
         double area, double *nVec, double smag, double *sVec, double alpha, const int nScal);

    void addViscFluxJacobiansCoupled(double **A0, double **A1,
    REALQ rho0_AD, REALQ *u0_AD, REALQ (&grad_u0_AD)[3][3], REALQ h0_AD, REALQ *grad_h0_AD, REALQ T0_AD, REALQ R0_AD, REALQ gam0_AD, REALQS *Scal0_AD, REALQS *dpress_dscal0_AD, REALQS kine0_AD,
    REALQ rho1_AD, REALQ *u1_AD, REALQ (&grad_u1_AD)[3][3], REALQ h1_AD, REALQ *grad_h1_AD, REALQ T1_AD, REALQ R1_AD, REALQ gam1_AD, REALQS *Scal1_AD, REALQS *dpress_dscal1_AD, REALQS kine1_AD,
    REALQ mul_AD,  REALQS mut_AD, REALQ lambdaOverCp_AD, REALQS kine_fa_AD, REALQ *u_fa_AD, REALQS *diff_AD, double *DiffTerm,
    double area, double *nVec, double smag, double *sVec, double alpha, const int nScal);


  // These should really be in UgpWithCv, but are currently here 
  // for convenience and to promote minimal code interference

  /**
   * 
   * 
   *  
   * GRADIENT ROUTINES, GRADIENT ROUTINES, GRADIENT ROUTINES, GRADIENT ROUTINES
   *  
   * 
   *  
   */



  void calcCv2Grad_AD(adouble (*gradPhi)[3], adouble *phi,
        const int limiter, adouble *refValue, const double epsilon)
  {
    calcCv2ScalGrad_AD(gradPhi, phi);

    switch (limiter)
    {
    case BARTH_JESPERSEN_MOD:
      limitCv2GradBJ_AD(gradPhi, phi, refValue, epsilon);
      break;

    default:
      break;
    }

   // Make sure we have gradients in ncv_g:ncv_gg-1

   double (*buf) = new double[ncv_ggff];
   FOR_I3 {
    FOR_ICV 
	buf[icv] = gradPhi[icv][i].value();

        updateCvDataG1G2(buf, REPLACE_DATA);

    FOR_ICV_G2_ONLY
	gradPhi[icv][i] = buf[icv];
   }

   delete [] buf;

  }

  void calcCv2Grad_AD(adouble (*gradPhi)[3][3], adouble (*phi)[3],
          const int limiter, adouble *refValue, const double epsilon)
  {
    calcCv2VecGrad_AD(gradPhi, phi);

    switch (limiter)
    {
    case BARTH_JESPERSEN_MOD:
      limitCv2GradBJ_AD(gradPhi, phi, refValue, epsilon);
      break;

    default:
     break;
    }

   // Make sure we have gradients in ncv_g:ncv_gg-1

   double (*buf) = new double[ncv_ggff];
   FOR_I3 {
    FOR_J3 {
     FOR_ICV 
	buf[icv] = gradPhi[icv][i][j].value();

        updateCvDataG1G2(buf, REPLACE_DATA);

    FOR_ICV_G2_ONLY
	gradPhi[icv][i][j] = buf[icv];
    }
   }

   delete [] buf;
  }


    void calcCv2ScalGrad_AD(adouble (*dphidxi)[3], adouble * phi)
   {
    // populates the gradient in cv range [0:ncv_g-1]. Note that this requires
    // phi to be defined over all cv's [0:ncv_ggff-1]...

    FOR_ICV_G {
      int noc_f = nbocv_all_i[icv];
      FOR_I3 dphidxi[icv][i] = nbocv_all_grad_coeff[noc_f][i]*phi[icv];

      int noc_l = nbocv_all_i[icv+1]-1;
      for (int noc = noc_f+1; noc <= noc_l; ++noc) {
        int icv_nbr = nbocv_all_v[noc];
        FOR_I3 dphidxi[icv][i] += nbocv_all_grad_coeff[noc][i]*phi[icv_nbr];
      }
    }
  }

  /**
   * populates the gradient in cv range [0:ncv_g-1]. Note that this requires
   * phi to be defined over all cv's [0:ncv_ggff-1]...
   */
  void calcCv2VecGrad_AD(adouble(*duidxj)[3][3], adouble(*u)[3])
  { 
    FOR_ICV_G
    { 
      int noc_f = nbocv_all_i[icv];
      FOR_I3
        FOR_J3
          duidxj[icv][i][j] = nbocv_all_grad_coeff[noc_f][j]*u[icv][i];
      
      int noc_l = nbocv_all_i[icv+1]-1;
      for (int noc = noc_f+1; noc<=noc_l; ++noc)
      { 
        int icv_nbr = nbocv_all_v[noc];
        FOR_I3
          FOR_J3
            duidxj[icv][i][j] += nbocv_all_grad_coeff[noc][j]*u[icv_nbr][i];
      }
    }
  }



  /**
   * apply smooth barth jespersen limiters
   */
  void limitCv2GradBJ_AD(adouble (*grad_p)[3], adouble *phi, adouble *refValue, const double epsilonSDWLS)
  {
    for (int icv = 0; icv < ncv_g; icv++)
    {
      adouble phiMax = phi[icv];
      adouble phiMin = phiMax;

      int noc_f = nbocv_all_i[icv];
      int noc_l = nbocv_all_i[icv + 1] - 1;
      // skip the diagonal in this loop...
      for (int noc = noc_f + 1; noc <= noc_l; noc++)
      {
        int icv_nbr = nbocv_all_v[noc];
        phiMax = max(phiMax, phi[icv_nbr]);
        phiMin = min(phiMin, phi[icv_nbr]);
      }

      adouble alfa = 1.0;

      if ((phiMax-phiMin) > 1.0e-12)
      {
        int foc_f = faocv_i[icv];
        int foc_l = faocv_i[icv + 1] - 1;
        for (int foc = foc_f; foc <= foc_l; foc++)
        {
          int ifa = faocv_v[foc];

          //if (ifa >= nfa_b)
          {
            double r0[3] = {0.0, 0.0, 0.0};
            vecMinVec3d(r0, x_fa[ifa], x_cv[icv]);
            adouble phifa = phi[icv] + vecDotVec3d_AD(r0, grad_p[icv]);
#ifdef BJ
            if ((phifa - phi[icv]) > 0.0)               alfa = min(alfa, (phiMax-phi[icv])/(phifa-phi[icv]));
            else if ((phifa - phi[icv]) < 0.0)          alfa = min(alfa, (phiMin-phi[icv])/(phifa-phi[icv]));
#else
            adouble dp;
            adouble eps = max(epsilonSDWLS*fabs(refValue[icv]), SDWLSClippingScalar);
            adouble eps2 = eps*eps;
            adouble dm = phifa - phi[icv];

            if   ((phifa - phi[icv]) >= 0.0)    dp = phiMax - phi[icv];
            else                                dp = phiMin - phi[icv];

            alfa = min(alfa, (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2));
#endif
          }
        }
      }
      else alfa = 0.0;

      for (int i=0; i<3; i++)
        grad_p[icv][i] *= alfa;
    }
  }

  /**
   * apply smooth barth jespersen limiters
   */
  void limitCv2GradBJ_AD(adouble (*grad_p)[3][3], adouble (*phi)[3], adouble *refValue, const double epsilonSDWLS)
  {
    for (int icv = 0; icv < ncv_g; icv++)
    {
      for (int i=0; i<3; i++)
      {
        adouble phiMax = phi[icv][i];
        adouble phiMin = phiMax;

        int noc_f = nbocv_all_i[icv];
        int noc_l = nbocv_all_i[icv + 1] - 1;
        // skip the diagonal in this loop...
        for (int noc = noc_f + 1; noc <= noc_l; noc++)
        {
          int icv_nbr = nbocv_all_v[noc];
          phiMax = max(phiMax, phi[icv_nbr][i]);
          phiMin = min(phiMin, phi[icv_nbr][i]);
        }

        adouble alfa = 1.0;

        if ((phiMax-phiMin) > 1.0e-8)
        {
          int foc_f = faocv_i[icv];
          int foc_l = faocv_i[icv + 1] - 1;
          for (int foc = foc_f; foc <= foc_l; foc++)
          {
            int ifa = faocv_v[foc];

            //if (ifa >= nfa_b)
            {
              double r0[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(r0, x_fa[ifa], x_cv[icv]);
              adouble phifa = phi[icv][i] + vecDotVec3d_AD(r0, grad_p[icv][i]);

#ifdef BJ
              if ((phifa - phi[icv][i]) > 0.0)               alfa = min(alfa, (phiMax-phi[icv][i])/(phifa-phi[icv][i]));
              else if ((phifa - phi[icv][i]) < 0.0)          alfa = min(alfa, (phiMin-phi[icv][i])/(phifa-phi[icv][i]));
#else
              adouble dp;
              adouble eps  = epsilonSDWLS*refValue[icv];
              adouble eps2 = eps*eps;
              adouble dm = phifa - phi[icv][i];

              if ((phifa - phi[icv][i]) >= 0.0)   dp = phiMax - phi[icv][i];
              else                                dp = phiMin - phi[icv][i];

              alfa = min(alfa, (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2));
#endif
            }
          }
        }
        else alfa = 0.0;

        for (int j=0; j<3; j++)
          grad_p[icv][i][j] *= alfa;
      }
    }
  }




//#ifdef USE_CONDASSIGN
//  adouble max(adouble a, adouble b) {
//	  return fmax(a, b);
//  }
//
//  adouble min(adouble a, adouble b) {
//	  return fmin(a, b);
//  }
//#else
  adouble max(adouble a, adouble b)
  {
    if (a >= b)    return a;
    else          return b;
  }

  adouble min(adouble a, adouble b)
  {
    if (a <= b)    return a;
    else          return b;
  }
//#endif

  double max(double a, double b)
  {
    if (a >= b)    return a;
    else          return b;
  }

  double min(double a, double b)
  {
    if (a <= b)    return a;
    else          return b;
  }
  
  


  /**
   * 
   * 
   *  
   * GRADIENT ROUTINES, GRADIENT ROUTINES, GRADIENT ROUTINES, GRADIENT ROUTINES
   *  
   * 
   *  
   */



adouble vecDotVec3d_AD(adouble * v1, adouble * v2) {

  adouble  dotproduct = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];

  return(dotproduct);
}

adouble vecDotVec3d_AD(double * v1, adouble * v2) {

  adouble  dotproduct = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];

  return(dotproduct);
}

adouble vecDotVec3d_AD(adouble * v1, double * v2) {

  adouble  dotproduct = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];

  return(dotproduct);
}

double vecDotVec3d_AD(double * v1, double * v2) {

   double dotproduct = 0.0;

  for (int i = 0; i < 3; i++)
    dotproduct += v1[i]*v2[i];

  return(dotproduct);
}

void vecMinVec3d_AD(adouble *vRes, adouble *v1, adouble *v2) {

  for (int i = 0; i < 3; i++)
    vRes[i] = v1[i] - v2[i];
}

 // Interpolate variable at the cell center from values at the faces (e.g., for viscosity)
  adouble InterpolateAtCellCenterFromFaceValues(adouble *phi, int icv)
  {
    adouble phiC1 = 0.0;

    int foc_f = faocv_i[icv];
    int foc_l = faocv_i[icv + 1] - 1;
    for (int foc=foc_f; foc<=foc_l; foc++)
    {
      int ifa = faocv_v[foc];
      phiC1 += phi[ifa];
    }
    adouble phiC = phiC1/(double)(foc_l-foc_f+1);

    return phiC;
  }

  virtual void calcStrainRate_AD()
  {
#ifdef USE_MEM_SAVING_ADVAR
    if (strMag.empty())
    {
      cout << "you have tried to calculate the magnitude of the strain rate, but it has not been registered!" << endl;
      throw(-1);
    }
#else
    if (strMag == NULL)
    {
      cout << "you have tried to calculate the magnitude of the strain rate, but it has not been registered!" << endl;
      throw(-1);
    }
#endif

    for (int icv=0; icv<ncv_gg; icv++)
    {
      strMag[icv] = 0.0;
      adouble diverg = grad_u[icv][0][0] + grad_u[icv][1][1] + grad_u[icv][2][2];

      for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        if (i == j)  strMag[icv] += pow(0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]) - 1./3.*diverg, 2.0);
        else         strMag[icv] += pow(0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]), 2.0);

      strMag[icv] = sqrt(2.0*strMag[icv]);
    }
  }

   virtual void calcStrainRateAndDivergence_AD()
  {
#ifdef USE_MEM_SAVING_ADVAR
	if (strMag.empty() || diverg.empty())
	{
	  cout << "you have tried to calculate the magnitude of the strain rate and the divergence of U, but one of them has not been registered!" << endl;
	  throw(-1);
	}
#else
    if ((strMag == NULL) || (diverg == NULL))
    {
      cout << "you have tried to calculate the magnitude of the strain rate and the divergence of U, but one of them has not been registered!" << endl;
      throw(-1);
    }
#endif

    for (int icv=0; icv<ncv_gg; icv++)
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


  virtual void calcVorticity_AD()
  {
#ifdef USE_MEM_SAVING_ADVAR
    if (vortMag.empty())
    {
      cout << "you have tried to calculate the magnitude of the vorticity, but it has not been registered!" << endl;
      throw(-1);
    }
#else
    if (vortMag == NULL)
    {
      cout << "you have tried to calculate the magnitude of the vorticity, but it has not been registered!" << endl;
      throw(-1);
    }
#endif

    for (int icv=0; icv<ncv_gg; icv++)
    {
      vortMag[icv] = 0.0;
      for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        vortMag[icv] += pow(0.5*(grad_u[icv][i][j] - grad_u[icv][j][i]), 2.0);

      vortMag[icv] = sqrt(2.0*vortMag[icv]);
    }
  }


};


#endif
