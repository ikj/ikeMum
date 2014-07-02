#ifndef JOEWITHMODELS_H
#define JOEWITHMODELS_H

#include "MiscUtils.h"
#include "tc_vec3d.h"

#include "UgpWithCvCompFlow.h"


class JoeWithModels: virtual public UgpWithCvCompFlow
{
public:
	double *Residual;

	// For output display (Metric to see the convergence of the simulation)
	double *log10_resid_rhoE;
	double *log10_resid_scalar0;

 private:
  
  int icv_debug;
  
  void setDebugCv() {
    
    icv_debug = -1;
    double x_debug[3] = {46.2617,0.5,0.376642};
    
    for (int icv = 0; icv < ncv_ggff; ++icv) {
      double d2 = 0.0;
      FOR_I3 {
	double dx = x_cv[icv][i] - x_debug[i];
	d2 += dx*dx;
      }
      if (sqrt(d2) < 1.0E-12) {
	cout << "found it on rank: " << mpi_rank << " " << x_cv[icv][0] << " " << x_cv[icv][1] << " " << x_cv[icv][2]<< " "<<icv<<endl;
	if (icv < ncv)
	  cout << " internal cv" << endl;
	else if (icv < ncv_g)
	  cout << " ghost level 1 cv" << endl;
	else if (icv < ncv_gg)
	  cout << " ghost level 2 cv" << endl;
	else if (icv < ncv_ggf)
	  cout << " fake level 1 cv" << endl;
	else
	  cout << " fake level 2 cv" << endl;
	icv_debug = icv;
      }
    }
    MPI_Pause("did we find one?");
    
  }

  void dumpDebugState(char * message) {
    
    if (icv_debug >= 0) 
      cout << "debug state: " << rho[icv_debug] << " " << rhou[icv_debug][0] << " " << rhou[icv_debug][1] << " " << rhoE[icv_debug] << endl;
    MPI_Pause(message);
  }
  
  void dumpDebugData(char * message,double * drho,double (*drhou)[3],double * drhoE) {
    if (icv_debug >= 0) 
      cout << "debug data: " << drho[icv_debug] << " " << drhou[icv_debug][0] << " " << drhou[icv_debug][1] << " " << drhoE[icv_debug] << endl;
    MPI_Pause(message);
  }
  

public:

  /**
   * constructor, pass ParamMap
   */
  JoeWithModels(ParamMap &p) : UgpWithCvCompFlow(p) {
	  nonVirtualInit();
	  init();
  }
  JoeWithModels(ParamMap &p, int s) : UgpWithCvCompFlow(p) {
	  nonVirtualInit();
	  init();
  }
  /**
   * constructor, pass name of joe's input file
   */
  JoeWithModels(char *name) : UgpWithCvCompFlow(name) {
	  nonVirtualInit();
	  init();
  }
  JoeWithModels(char *name, int s) : UgpWithCvCompFlow(name) {
	  nonVirtualInit();
	  init();
  }

  void nonVirtualInit() {
		// For output display (Metric to see the convergence of the simulation)
		log10_resid_rhoE  = NULL;  		registerScalar(log10_resid_rhoE,  "LOG10_RESID_RHOE",  CV_DATA);
		log10_resid_scalar0 = NULL;  	registerScalar(log10_resid_scalar0, "LOG10_RESID_SCALAR0", CV_DATA);
  }

  virtual void init()
  {
    if (mpi_rank == 0)
      cout << "JoeWithModels virtual init()"<< endl;
  }

  virtual ~JoeWithModels() {}

public:

  int doneSolver(const double residEnergy)
  {
    // ---------------------------------------------
    // returns 1 if we are done, 0 otherwise...
    // ---------------------------------------------
    int done = 0;

    if (mpi_rank==0)
    {

      // if the current completed step is greater or equal to the user
      // specified nsteps, then we are done...

      if ((nsteps >= 0) && (step >= nsteps))
      {
        cout << "step: " << step << endl;
        cout << "nsteps: " << nsteps << endl;
        cout << " > step from restart file is greater than NSTEP in Joe.in, you can reset the step by using RESET_STEP in Joe.in. stopping."<<endl;
        done = 1;
      }

      // check energy residual
      if (residEnergy < resid_energ_th)
      {
        cout<<" > reached threshold for energy residual. stopping."<<endl;
        done = 1;
      }

      // check if energy residual is NAN (not a number)
      if (residEnergy != residEnergy)
      {
        cout<<" > nan in energy residual detected. stopping."<<endl;
        done = 1;
      }

      // RUNTIME...
	/*
      if ((runtime_flag) && (MPI_Wtime() > runtime))
      {
        cout<<" > reached specified RUNTIME. stopping."<<endl;
        done = 1;
      }*/

      // killjoe file...
      if (fileExists("killjoe"))
      {
        MPI_File_delete("killjoe", MPI_INFO_NULL);
        cout<<" > Found file killjoe. stopping."<<endl;
        done = 1;
      }
    }
    MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

       return (done);
  }

  /**
   * hook for initializing Navier-Stokes equations
   */
  virtual void initialHook()
  {
    if (mpi_rank == 0)
      cout << "setInitialCondition()"<< endl;

    double rho_init, p_init, u_init[3] = {0.0, 0.0, 0.0};

    Param *p;
    if (getParam(p, "RHO_INITIAL"))       rho_init = p->getDouble(1);
    else                                  rho_init = rho_ref;

    if (getParam(p, "P_INITIAL"))         p_init = p->getDouble(1);
    else                                  p_init = p_ref;

    if (getParam(p, "U_INITIAL")) 
    {
      u_init[0] = p->getDouble(1);
      u_init[1] = p->getDouble(2);
      u_init[2] = p->getDouble(3);
    }

    if (!checkDataFlag(rho))
      for (int icv=0; icv<ncv; icv++)
        rho[icv] = rho_init;

    if (!checkDataFlag(rhou))
      for (int icv=0; icv<ncv; icv++)
        for (int i=0; i<3; i++)
          rhou[icv][i] = rho_init*u_init[i];

    if (!checkDataFlag(rhoE))
      for (int icv=0; icv<ncv; icv++)
        rhoE[icv] = p_init/(gamma[icv]-1.0)+ 0.5*rho_init*vecDotVec3d(u_init, u_init);

//    writeData(0);
  }


  /**
   * boundary HOOK
   */
  virtual void boundaryHook(double *rho, double (*vel)[3], double *press, FaZone *zone)  { /*empty*/ }


  /**
   * source HOOK, rhs already contains fluxes, don't overwrite !!!!!
   * e.g.: gravity
   *
   *  for (int icv = 0; icv < ncv; icv++)
   *  {
   *    for (int i = 0; i < 3; i++)
   *    rhs_rhou[icv][i] += cv_volume[icv]*rho[icv]*gravity[i];
   *    rhs_rhoE[icv] += cv_volume[icv]*(gravity[0]*rhou[icv][0]+gravity[1]*rhou[icv][1]+gravity[2]*rhou[icv][2]);
   *  }
   */
  virtual void sourceHook(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5])  { /*empty*/ }


  /**
   * minimum stuff to run a compressible Navier-Stokes simulation
   */
  virtual void run();

  /**
   * explicit backward Euler
   */
  void runForwardEuler();

  /**
   * explicit 3-step Runge-Kutta
   */
  void runRK();

  /**
   * global implicit fully-coupled relaxation 
   */
  void runBackwardEuler();
  
  /**
   * global implicit fully-coupled relaxation for unsteady calcs 
   */
  void runBDF2();

  /**
   * set boundary conditions for Navier-Stokes and scalars 
   */
  void setBC();
  
  /**
   * calculate RHS, routine used for explicit and implicit calculations,
   * flagImplicit = 1 ... implicit otherwise set to zero
   */
  int calcRhs(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double **rhs_rhoScal, double (*A)[5][5], double ***AScal, int flagImplicit);
  
  /**
   * calculate Euler flux for both NSE and scalars, routine used for explicit and implicit calculations,
   * flagImplicit = 1 ... implicit otherwise set to zero
   */
  int calcEulerFlux(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double **rhs_rhoScal, double (*A)[5][5], double ***AScal, int flagImplicit);

  /**
   * calculate viscous flux for NSE only, routine used for explicit and implicit calculations,
   * flagImplicit = 1 ... implicit otherwise set to zero
   */
  void calcViscousFluxNS(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5], int flagImplicit);
  
  /**
   * show residual
   */
  void showResidue(double *rhsResid);
  

  /*
   * add your own temporal output here
   */
  virtual void temporalHook() { /*empty*/ }

  /*
   * add your own finalization here
   */
  virtual void finalHook() { /*empty*/ }


      void runBackwardEulerCoupled();

      int calcRhsCoupled(double **rhs, double ***A, int nScal, int flagImplicit);

      int calcFluxCoupled(double **rhs, double ***A, int nScal, int flagImplicit);

      virtual void sourceHookCoupled(double **rhs, double ***A, int nScal, int flagImplicit)  { /*empty*/ }

};


#endif



