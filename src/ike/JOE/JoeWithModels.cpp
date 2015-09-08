#include "JoeWithModels.h"
#include "myMem.h"

//IKE
#include <set>

void JoeWithModels::run()
{
  // read mesh or restart file
  initializeFromRestartFile(getStringParam("RESTART"));

  if (checkParam("RESET_STEP")) {
	  step = 0;
	  if (mpi_rank == 0)
		  cout << "RESET_STEP: " << step << endl;

	  time = 0.0;
  }

  // initialize models
  initialHookScalarRansTurbModel();
  initialHookScalarRansCombModel();
  for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
    updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);

  // initialize Navier-Stokes equations
  initialHook(); 

  // init probes 
  initProbes(this);


  updateCvDataG1G2(rho, REPLACE_DATA);
  updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
  updateCvDataG1G2(rhoE, REPLACE_DATA);

  string tIntName = getStringParam("TIME_INTEGRATION");

  if (tIntName == "FORWARD_EULER")                 runForwardEuler();
  else if (tIntName == "RK")                       runRK();
  else if (tIntName == "BACKWARD_EULER")           runBackwardEuler();
  else if (tIntName == "BACKWARD_EULER_COUPLED")   runBackwardEulerCoupled();
  else if (tIntName == "BDF2")                     runBDF2();
  else
    if (mpi_rank == 0)
    {
      cerr << "ERROR: wrong time integration scheme specified !" << endl;
      cerr << "available integration schemes are: FORWARD_EULER, RK, BACKWARD_EULER, BDF2" << endl;
    }
}


void JoeWithModels::runForwardEuler()
{

  //setDebugCv();

  int nScal = scalarTranspEqVector.size();

  double *myResidual = new double[5+nScal];
  Residual = new double[5+nScal];

  double *drho = new double[ncv];
  double (*drhou)[3] = new double[ncv][3];
  double *drhoE = new double[ncv];
  double **drhoScal = NULL;
  if (nScal > 0) getMem2D(&drhoScal, 0, nScal-1, 0, ncv-1, "drhoScal");

  double (*dummyA)[5][5] = NULL;       // provide dummy pointer for explicit mode!
  double ***dummyAScal   = NULL;       // provide dummy pointer for explicit mode!

  // -------------------------------------------------------------------------------------------
  // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
  // -------------------------------------------------------------------------------------------
  calcStateVariables();

  // -------------------------------------------------------------------------------------------
  // update material properties: laminar viscosity and heat conductivity
  // -------------------------------------------------------------------------------------------
  calcMaterialProperties();

  // -------------------------------------------------------------------------------------------
  // set BC's for NS and scalars
  // -------------------------------------------------------------------------------------------
  setBC();

  // -------------------------------------------------------------------------------------------
  // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
  // -------------------------------------------------------------------------------------------
  calcRansTurbViscMuet();


  // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  // Loop over time steps
  // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  int done = 0;
  if ((nsteps >= 0)&&(step >= nsteps))  done = 1;

  if (initial_flowfield_output == "YES")
    writeData(0);

  while (done != 1)
  {
    step++;
    double dt = calcDt(cfl);
    time += dt;


    // =========================================================================================
    // EXPLICIT FORWARD STEP
    // =========================================================================================

    //dumpDebugState("1");
    
    calcRhs(drho, drhou, drhoE, drhoScal, dummyA, dummyAScal, 0);

    //dumpDebugData("Deltas: ",drho,drhou,drhoE);

	/*
    if(mpi_rank==0) cout<<"NOT UPDATING FLOW"<<endl;
	*/
    for (int icv = 0; icv < ncv; icv++)
    {
      double tmp = local_dt[icv]/cv_volume[icv];
      rho[icv] += tmp * drho[icv];
      rhoE[icv] += tmp * drhoE[icv];
      for (int i = 0; i < 3; i++)
        rhou[icv][i] += tmp * drhou[icv][i];
    }

    for (int iScal = 0; iScal < nScal; iScal++)
    {
      double *phi = scalarTranspEqVector[iScal].phi;

      for (int icv = 0; icv < ncv; icv++)
      {
        double tmp = local_dt[icv]/cv_volume[icv];
        phi[icv] = ((rho[icv] - tmp * drho[icv]) * phi[icv] + tmp * drhoScal[iScal][icv]) / rho[icv];
      }
    }

    // -------------------------------------------------------------------------------------------
    // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
    // -------------------------------------------------------------------------------------------
    updateCvDataG1G2(rho, REPLACE_DATA);
    updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
    updateCvDataG1G2(rhoE, REPLACE_DATA);

    for (int iScal = 0; iScal < nScal; iScal++)
      updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);

    // -------------------------------------------------------------------------------------------
    // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
    // -------------------------------------------------------------------------------------------
    calcStateVariables();

    // -------------------------------------------------------------------------------------------
    // update material properties: laminar viscosity and heat conductivity
    // -------------------------------------------------------------------------------------------
    calcMaterialProperties();

    // -------------------------------------------------------------------------------------------
    // set BC's for NS and scalars
    // -------------------------------------------------------------------------------------------
    setBC();

    // -------------------------------------------------------------------------------------------
    // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
    // -------------------------------------------------------------------------------------------
    calcRansTurbViscMuet();


    // =========================================================================================
    // calculate and show residual
    // =========================================================================================
    for (int i = 0; i < 5+nScal; i++)
    {
      myResidual[i] = 0.0;
      Residual[i] = 0.0;
    }

    for (int icv = 0; icv < ncv; icv++)
    {
      myResidual[0] += fabs(drho[icv]);
      for (int i=0; i<3; i++)
        myResidual[i+1] += fabs(drhou[icv][i]);
      myResidual[4] += fabs(drhoE[icv]);
    }
    for (int iScal = 0; iScal < nScal; iScal++)
      for (int icv = 0; icv < ncv; icv++)
        myResidual[5+iScal] += fabs(drhoScal[iScal][icv]);

    MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (step%check_interval == 0)
    {
      if ((mpi_rank == 0) && (step%(check_interval*10) == 0))
        cout << "\ndone step: "<< step << ", cfl: " << cfl << ", min. dt: " << dt << " time: " << time << endl;

      showResidue(Residual);
    }

    // IKJ: For output display (Metric to see the convergence of the simulation)
    for (int icv = 0; icv < ncv; icv++) {
    	residField[icv]       = drhoE[icv];
    	log10_resid_rhoE[icv] = log10(fabs(drhoE[icv]) + 1.0e-15);
    }
    for (int iScal = 0; iScal < nScal; iScal++)
    	for (int icv = 0; icv < ncv; icv++)
    		log10_resid_scalar0[icv] = log10(fabs(drhoScal[iScal][icv]) + 1.0e-15);
    updateCvDataG1G2(residField,       REPLACE_DATA);
    updateCvDataG1G2(log10_resid_rhoE, REPLACE_DATA);
    if(nScal>0)
    	updateCvDataG1G2(log10_resid_scalar0, REPLACE_DATA);

    // IKJ
    double myTotResid_rhs = 0.0;
    double myTotResid_dq  = 0.0;
    for (int icv = 0; icv < ncv; icv++) {
    	double localTotResid_rhs = 0.0;

    	localTotResid_rhs += fabs(drho[icv]);
    	for (int i=0; i<3; i++)
    		localTotResid_rhs += fabs(drhou[icv][i]);
    	localTotResid_rhs += fabs(drhoE[icv]);

    	double tmp = local_dt[icv]/cv_volume[icv];

    	myTotResid_rhs +=       localTotResid_rhs;
    	myTotResid_dq  += tmp * localTotResid_rhs;
    }
	for (int iScal = 0; iScal < nScal; iScal++)
		for (int icv = 0; icv < ncv; icv++) {
			double tmp = local_dt[icv]/cv_volume[icv];

			myTotResid_rhs +=       fabs(drhoScal[iScal][icv]);
			myTotResid_dq  += tmp * fabs(drhoScal[iScal][icv]);
		}
    MPI_Allreduce(&myTotResid_dq,  &totResid_dq,  1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myTotResid_rhs, &totResid_rhs, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    temporalHook();

    writeData(step);

    if ((write_restart > 0) && (step % write_restart == 0))
      writeRestart(step);

    if ((nsteps != -1) && (step >= nsteps))
      done = 1;

  }

  writeRestart();

  finalHook();

  delete [] drho;
  delete [] drhou;
  delete [] drhoE;
  if (nScal > 0) freeMem2D(drhoScal, 0, nScal-1, 0, ncv-1);

  delete [] myResidual;
  delete [] Residual; 	Residual = NULL;
}

void JoeWithModels::runRK()
{
  int nScal = scalarTranspEqVector.size();

  double *myResidual = new double[5+nScal];
  Residual = new double[5+nScal];

  double *drho1 = new double[ncv];
  double (*drhou1)[3] = new double[ncv][3];
  double *drhoE1 = new double[ncv];
  double **drhoScal1 = NULL;
  if (nScal > 0) getMem2D(&drhoScal1, 0, nScal-1, 0, ncv-1, "drhoScal1");

  double *drho2 = new double[ncv];
  double (*drhou2)[3] = new double[ncv][3];
  double *drhoE2 = new double[ncv];
  double **drhoScal2 = NULL;
  if (nScal > 0) getMem2D(&drhoScal2, 0, nScal-1, 0, ncv-1, "drhoScal2");

  double *drho3 = new double[ncv];
  double (*drhou3)[3] = new double[ncv][3];
  double *drhoE3 = new double[ncv];
  double **drhoScal3 = NULL;
  if (nScal > 0) getMem2D(&drhoScal3, 0, nScal-1, 0, ncv-1, "drhoScal3");

  double *rho0 = new double[ncv];
  double (*rhou0)[3] = new double[ncv][3];
  double *rhoE0 = new double[ncv];
  double **rhoScal0 = NULL;
  if (nScal > 0) getMem2D(&rhoScal0, 0, nScal-1, 0, ncv-1, "rhoScal0");

  double (*dummyA)[5][5] = NULL;       // provide dummy pointer for explicit mode!
  double ***dummyAScal   = NULL;       // provide dummy pointer for explicit mode!
  // -------------------------------------------------------------------------------------------
  // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
  // -------------------------------------------------------------------------------------------
  calcStateVariables();

  // -------------------------------------------------------------------------------------------
  // update material properties: laminar viscosity and heat conductivity
  // -------------------------------------------------------------------------------------------
  calcMaterialProperties();

  // -------------------------------------------------------------------------------------------
  // set BC's for NS and scalars
  // -------------------------------------------------------------------------------------------
  setBC();

  // -------------------------------------------------------------------------------------------
  // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
  // -------------------------------------------------------------------------------------------
  calcRansTurbViscMuet();


  // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  // Loop over time steps
  // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  int done = 0;
  if ((nsteps >= 0)&&(step >= nsteps))  done = 1;

  if (initial_flowfield_output == "YES")
    writeData(0);

  while (done != 1)
  {
    step++;
    double dt = calcDt(cfl);
    time += dt;


    // =========================================================================================
    // copy the current solution into rho0, rhou0, etc...
    // =========================================================================================
    for (int icv = 0; icv < ncv; icv++)
    {
      rho0[icv] = rho[icv];
      rhoE0[icv] = rhoE[icv];
      for (int i = 0; i < 3; i++)
        rhou0[icv][i] = rhou[icv][i];
    }
    for (int iScal = 0; iScal < nScal; iScal++)
    {
      double *phi = scalarTranspEqVector[iScal].phi;
      for (int icv = 0; icv < ncv; icv++)
        rhoScal0[iScal][icv] = rho[icv] * phi[icv];
    }


    // =========================================================================================
    // RK STEP 1
    // =========================================================================================
    calcRhs(drho1, drhou1, drhoE1, drhoScal1, dummyA, dummyAScal, 0);

    for (int icv = 0; icv < ncv; icv++)
    {
      double tmp = local_dt[icv]/cv_volume[icv];
      drho1[icv] *= tmp;
      drhoE1[icv] *= tmp;
      for (int i = 0; i < 3; i++)
        drhou1[icv][i] *= tmp;
    }
    for (int iScal = 0; iScal < nScal; iScal++)
      for (int icv = 0; icv < ncv; icv++)
      {
        double tmp = local_dt[icv]/cv_volume[icv];
        drhoScal1[iScal][icv] *= tmp;
      }

    for (int icv = 0; icv < ncv; icv++)
    {
      rho[icv] = rho0[icv] + 0.5*drho1[icv];
      rhoE[icv] = rhoE0[icv] + 0.5*drhoE1[icv];
      for (int i = 0; i < 3; i++)
        rhou[icv][i] = rhou0[icv][i] + 0.5*drhou1[icv][i];
    }

    // -------------------------------------------------------------------------------------------
    // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
    // -------------------------------------------------------------------------------------------
    updateCvDataG1G2(rho, REPLACE_DATA);
    updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
    updateCvDataG1G2(rhoE, REPLACE_DATA);

    for (int iScal = 0; iScal < nScal; iScal++)
      updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);

    // -------------------------------------------------------------------------------------------
    // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
    // -------------------------------------------------------------------------------------------
    calcStateVariables();

    // -------------------------------------------------------------------------------------------
    // update material properties: laminar viscosity and heat conductivity
    // -------------------------------------------------------------------------------------------
    calcMaterialProperties();

    // -------------------------------------------------------------------------------------------
    // set BC's for NS and scalars
    // -------------------------------------------------------------------------------------------
    setBC();

    // -------------------------------------------------------------------------------------------
    // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
    // -------------------------------------------------------------------------------------------
    calcRansTurbViscMuet();


    // =========================================================================================
    // RK STEP 2
    // =========================================================================================
    calcRhs(drho2, drhou2, drhoE2, drhoScal2, dummyA, dummyAScal, 0);

    for (int icv = 0; icv < ncv; icv++) {
      double tmp = local_dt[icv]/cv_volume[icv];
      drho2[icv] *= tmp;
      drhoE2[icv] *= tmp;
      for (int i = 0; i < 3; i++)
        drhou2[icv][i] *= tmp;
    }
    for (int iScal = 0; iScal < nScal; iScal++)
      for (int icv = 0; icv < ncv; icv++)
      {
        double tmp = local_dt[icv]/cv_volume[icv];
        drhoScal2[iScal][icv] *= tmp;
      }

    for (int icv = 0; icv < ncv; icv++) {
      rho[icv] = rho0[icv] - drho1[icv] + 2.0*drho2[icv];
      rhoE[icv] = rhoE0[icv] - drhoE1[icv] + 2.0*drhoE2[icv];
      for (int i = 0; i < 3; i++)
        rhou[icv][i] = rhou0[icv][i] - drhou1[icv][i] + 2.0*drhou2[icv][i];
    }

    for (int iScal = 0; iScal < nScal; iScal++)
    {
      double *phi = scalarTranspEqVector[iScal].phi;
      for (int icv = 0; icv < ncv; icv++)
        phi[icv] = (rhoScal0[iScal][icv] - drhoScal1[iScal][icv] + 2.0*drhoScal2[iScal][icv]) / rho[icv];
    }

    // -------------------------------------------------------------------------------------------
    // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
    // -------------------------------------------------------------------------------------------
    updateCvDataG1G2(rho, REPLACE_DATA);
    updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
    updateCvDataG1G2(rhoE, REPLACE_DATA);

    for (int iScal = 0; iScal < nScal; iScal++)
      updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);

    // -------------------------------------------------------------------------------------------
    // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
    // -------------------------------------------------------------------------------------------
    calcStateVariables();

    // -------------------------------------------------------------------------------------------
    // update material properties: laminar viscosity and heat conductivity
    // -------------------------------------------------------------------------------------------
    calcMaterialProperties();

    // -------------------------------------------------------------------------------------------
    // set BC's for NS and scalars
    // -------------------------------------------------------------------------------------------
    setBC();

    // -------------------------------------------------------------------------------------------
    // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
    // -------------------------------------------------------------------------------------------
    calcRansTurbViscMuet();


    // =========================================================================================
    // RK STEP 3
    // =========================================================================================
    calcRhs(drho3, drhou3, drhoE3, drhoScal3, dummyA, dummyAScal, 0);
    for (int icv = 0; icv < ncv; icv++) {
      double tmp = local_dt[icv]/cv_volume[icv];
      drho3[icv] *= tmp;
      drhoE3[icv] *= tmp;
      for (int i = 0; i < 3; i++)
        drhou3[icv][i] *= tmp;
    }
    for (int iScal = 0; iScal < nScal; iScal++)
      for (int icv = 0; icv < ncv; icv++)
      {
        double tmp = local_dt[icv]/cv_volume[icv];
        drhoScal3[iScal][icv] *= tmp;
      }

    for (int icv = 0; icv < ncv; icv++) {
      rho[icv] = rho0[icv] + (drho1[icv] + 4.0*drho2[icv] + drho3[icv])/6.0;
      rhoE[icv] = rhoE0[icv] + (drhoE1[icv] + 4.0*drhoE2[icv] + drhoE3[icv])/6.0;
      for (int i = 0; i < 3; i++)
        rhou[icv][i] = rhou0[icv][i] + (drhou1[icv][i] + 4.0*drhou2[icv][i] + drhou3[icv][i])/6.0;
    }

    for (int iScal = 0; iScal < nScal; iScal++)
    {
      double *phi = scalarTranspEqVector[iScal].phi;
      for (int icv = 0; icv < ncv; icv++)
        phi[icv] = (rhoScal0[iScal][icv] + (drhoScal1[iScal][icv] + 4.0*drhoScal2[iScal][icv] + drhoScal3[iScal][icv]) / 6.0) / rho[icv];
    }

    // -------------------------------------------------------------------------------------------
    // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
    // -------------------------------------------------------------------------------------------
    updateCvDataG1G2(rho, REPLACE_DATA);
    updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
    updateCvDataG1G2(rhoE, REPLACE_DATA);

    for (int iScal = 0; iScal < nScal; iScal++)
      updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);

    // -------------------------------------------------------------------------------------------
    // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
    // -------------------------------------------------------------------------------------------
    calcStateVariables();

    // -------------------------------------------------------------------------------------------
    // update material properties: laminar viscosity and heat conductivity
    // -------------------------------------------------------------------------------------------
    calcMaterialProperties();

    // -------------------------------------------------------------------------------------------
    // set BC's for NS and scalars
    // -------------------------------------------------------------------------------------------
    setBC();

    // -------------------------------------------------------------------------------------------
    // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
    // -------------------------------------------------------------------------------------------
    calcRansTurbViscMuet();


    // =========================================================================================
    // calculate and show residual
    // =========================================================================================
    for (int i = 0; i < 5+nScal; i++)
    {
      myResidual[i] = 0.0;
      Residual[i] = 0.0;
    }

    for (int icv = 0; icv < ncv; icv++)
    {
      myResidual[0] += fabs(drho1[icv] + 4.0*drho2[icv] + drho3[icv]) / 6.0;
      for (int i=0; i<3; i++)
        myResidual[i+1] += fabs(drhou1[icv][i] + 4.0*drhou2[icv][i] + drhou3[icv][i]) / 6.0;
      myResidual[4] += fabs(drhoE1[icv] + 4.0*drhoE2[icv] + drhoE3[icv]) / 6.0;
    }
    for (int iScal = 0; iScal < nScal; iScal++)
      for (int icv = 0; icv < ncv; icv++)
        myResidual[5+iScal] += fabs(drhoScal1[iScal][icv] + 4.0*drhoScal2[iScal][icv] + drhoScal3[iScal][icv]) / 6.0;

    MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (step%check_interval == 0)
    {
      if ((mpi_rank == 0) && (step%(check_interval*10) == 0))
        cout << "\ndone step: "<< step << ", cfl: " << cfl << ", min. dt: " << dt << " time: " << time << endl;

      showResidue(Residual);
    }

    // IKJ: For output display (Metric to see the convergence of the simulation)
    for (int icv = 0; icv < ncv; icv++) {
    	double tmpInverse = cv_volume[icv] / local_dt[icv];
    	residField[icv]       = tmpInverse * (drhoE1[icv] + 4.0*drhoE2[icv] + drhoE3[icv]) / 6.0;
    	log10_resid_rhoE[icv] = log10(fabs(residField[icv]) + 1.0e-15);
    }
    for (int iScal = 0; iScal < nScal; iScal++)
    	for (int icv = 0; icv < ncv; icv++) {
    		double tmpInverse = cv_volume[icv] / local_dt[icv];
    		log10_resid_scalar0[icv] = log10(tmpInverse * (drhoScal1[iScal][icv] + 4.0*drhoScal2[iScal][icv] + drhoScal3[iScal][icv]) / 6.0 + 1.0e-15);
    	}
    updateCvDataG1G2(residField,       REPLACE_DATA);
    updateCvDataG1G2(log10_resid_rhoE, REPLACE_DATA);
    if(nScal>0)
    	updateCvDataG1G2(log10_resid_scalar0, REPLACE_DATA);

    // IKJ: Note -- Unlike in ForwardEuler, drho stuff was already scaled by tmp in RK
    double myTotResid_dq  = 0.0;
    double myTotResid_rhs = 0.0;
    for (int icv = 0; icv < ncv; icv++) {
    	double localTotResid_dq = 0.0;

    	localTotResid_dq += fabs(drho1[icv] + 4.0*drho2[icv] + drho3[icv]) / 6.0;
    	for (int i=0; i<3; i++)
    		localTotResid_dq += fabs(drhou1[icv][i] + 4.0*drhou2[icv][i] + drhou3[icv][i]) / 6.0;
    	localTotResid_dq += fabs(drhoE1[icv] + 4.0*drhoE2[icv] + drhoE3[icv]) / 6.0;

    	double tmpInverse = cv_volume[icv] / local_dt[icv];

    	myTotResid_rhs += tmpInverse * localTotResid_dq;
    	myTotResid_dq  +=              localTotResid_dq;
    }
	for (int iScal = 0; iScal < nScal; iScal++)
		for (int icv = 0; icv < ncv; icv++) {
			double tmpInverse = cv_volume[icv] / local_dt[icv];

			myTotResid_rhs += tmpInverse * fabs(drhoScal1[iScal][icv] + 4.0*drhoScal2[iScal][icv] + drhoScal3[iScal][icv]) / 6.0;
			myTotResid_dq  +=              fabs(drhoScal1[iScal][icv] + 4.0*drhoScal2[iScal][icv] + drhoScal3[iScal][icv]) / 6.0;
		}
    MPI_Allreduce(&myTotResid_dq,  &totResid_dq,  1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myTotResid_rhs, &totResid_rhs, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    temporalHook();

    writeData(step);

    if ((write_restart > 0) && (step % write_restart == 0))
      writeRestart(step);

    if ((nsteps != -1) && (step >= nsteps))
      done = 1;

  }

  writeRestart();

  finalHook();

  delete [] drho1;
  delete [] drhou1;
  delete [] drhoE1;
  if (nScal > 0) freeMem2D(drhoScal1, 0, nScal-1, 0, ncv-1);

  delete [] drho2;
  delete [] drhou2;
  delete [] drhoE2;
  if (nScal > 0) freeMem2D(drhoScal2, 0, nScal-1, 0, ncv-1);

  delete [] drho3;
  delete [] drhou3;
  delete [] drhoE3;
  if (nScal > 0) freeMem2D(drhoScal3, 0, nScal-1, 0, ncv-1);

  delete [] rho0;
  delete [] rhou0;
  delete [] rhoE0;
  if (nScal > 0) freeMem2D(rhoScal0, 0, nScal-1, 0, ncv-1);

  delete [] myResidual;
  delete [] Residual; 	Residual = NULL;
}




void JoeWithModels::runBackwardEuler()
{
  int nScal = scalarTranspEqVector.size();

  double *myResidual = new double[5+nScal];
  Residual           = new double[5+nScal];

  double *RHSrho = new double[ncv];
  double (*RHSrhou)[3] = new double[ncv][3];
  double *RHSrhoE = new double[ncv];

  double (*A)[5][5] = new double[nbocv_s][5][5];
  double (*dq)[5]   = new double[ncv_g][5];
  double (*rhs)[5]  = new double[ncv][5];

  double ***AScal  = NULL;  if (nScal > 0) getMem3D(&AScal,   0, nScal-1, 0, 5, 0, nbocv_s-1, "AScal");
  double **dScal   = NULL;  if (nScal > 0) getMem2D(&dScal,   0, nScal-1, 0, ncv_g-1, "dScal");
  double **rhsScal = NULL;  if (nScal > 0) getMem2D(&rhsScal, 0, nScal-1, 0, ncv-1, "rhsScal");


  //------------------------------------
  // some parameters that are only relevant for backward euler
  //------------------------------------
  double underRelax = getDoubleParam("UNDER_RELAXATION", "0.3");

  if (!checkParam("LINEAR_SOLVER_NS_TRESHOLDS"))
  {
    ParamMap::add("LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2");    // add default values
    if (mpi_rank == 0)
      cout << "WARNING: added keyword \"LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2\"" <<
              " to parameter map!" << endl;
  }
  int maxIterLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getInt("MAX_ITER");
  double zeroAbsLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("ABS_RESID");
  double zeroRelLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("REL_RESID");

  if (!checkParam("CFL_RAMP"))
  {
    ParamMap::add("CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0");    // add default: no increase of CFL number!
    if (mpi_rank == 0)
      cout << "WARNING: added \"CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0\" to parameter map" << endl;
  }

  int startIncCFL = getParam("CFL_RAMP")->getInt("AFTER_ITER");
  int intervalIncCFL = getParam("CFL_RAMP")->getInt("INTERVAL_ITER");
  double incCFL = getParam("CFL_RAMP")->getDouble("FACTOR_CFL");
  double maxCFL = getParam("CFL_RAMP")->getDouble("MAX_CFL");

  // -------------------------------------------------------------------------------------------
  // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
  // -------------------------------------------------------------------------------------------
  calcStateVariables();

  // -------------------------------------------------------------------------------------------
  // update material properties: laminar viscosity and heat conductivity
  // -------------------------------------------------------------------------------------------
  calcMaterialProperties();

  // -------------------------------------------------------------------------------------------
  // set BC's for NS and scalars
  // -------------------------------------------------------------------------------------------
  setBC();

  // -------------------------------------------------------------------------------------------
  // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
  // -------------------------------------------------------------------------------------------
  calcRansTurbViscMuet();

  // -------------------------------------------------------------------------------------------
  //
  //   Loop over time steps
  //
  // -------------------------------------------------------------------------------------------

  int done = 0;

  if ((nsteps >= 0)&&(step >= nsteps))  done = 1;

  if (initial_flowfield_output == "YES")
    writeData(0);

  // provide total runtime
  double wtime, wtime0;
  MPI_Barrier(mpi_comm);
  if (mpi_rank == 0)
    wtime = MPI_Wtime();

  while (done != 1)
  {
    step++;
    if ((step >= startIncCFL) && (step%intervalIncCFL == 0) && (cfl < maxCFL))      cfl *= incCFL;
    double dt_min = calcDt(cfl);

    // ---------------------------------------------------------------------------------
    // Compute RHS for both NSE and scalars
    // ---------------------------------------------------------------------------------
    for (int noc = 0; noc < nbocv_s; noc++)                             // set A, dq to zero! rhs is set zero in calcRHS
      for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
          A[noc][i][j] = 0.0;
    for (int icv = 0; icv < ncv_g; icv++)
      for (int i = 0; i < 5; i++)
        dq[icv][i] = 0.0;

    for (int iScal = 0; iScal < nScal; iScal++)                         // set AScal, dScal to zero! rhs is set zero in calcRHS
    {
      for (int i = 0; i <= 5; i++)
        for (int noc = 0; noc < nbocv_s; noc++)
          AScal[iScal][i][noc] = 0.0;
      for (int icv = 0; icv < ncv_g; icv++)
        dScal[iScal][icv] = 0.0;
    }


    // ---------------------------------------------------------------------------------
    // calculate RHS for both NSE and scalars
    // ---------------------------------------------------------------------------------
    calcRhs(RHSrho, RHSrhou, RHSrhoE, rhsScal, A, AScal, true);

    // ---------------------------------------------------------------------------------
    // solve linear system for the NSE
    // ---------------------------------------------------------------------------------

    for (int icv=0; icv<ncv; ++icv)                                 // prepare rhs and A
    {
      rhs[icv][0] = underRelax*RHSrho[icv];
      rhs[icv][1] = underRelax*RHSrhou[icv][0];
      rhs[icv][2] = underRelax*RHSrhou[icv][1];
      rhs[icv][3] = underRelax*RHSrhou[icv][2];
      rhs[icv][4] = underRelax*RHSrhoE[icv];

      residField[icv] = RHSrhoE[icv];

      double tmp = cv_volume[icv]/(local_dt[icv]);
      for (int i = 0; i < 5; i++)
        A[nbocv_i[icv]][i][i] += tmp;                               // diagonal part ( vol/dt + A )
    }

    solveCoupledLinSysNS(dq, A, rhs, zeroAbsLS, zeroRelLS, maxIterLS);  // solve linear system

	/*
    if(mpi_rank==0) cout<<"NOT UPDATING FLOW"<<endl;
	*/
    for (int icv=0; icv<ncv; icv++)                                 // update solution: Q_new = Q_old + Delta_Q
    {
      rho[icv]     += dq[icv][0];
      rhou[icv][0] += dq[icv][1];
      rhou[icv][1] += dq[icv][2];
      rhou[icv][2] += dq[icv][3];
      rhoE[icv]    += dq[icv][4];
    }
    

    UpdateCvDataStateVec(dq);                                       // update dq since neighbors needed to compute RHS of scalars


    // ---------------------------------------------------------------------------------
    // solve linear system for the scalars
    // ---------------------------------------------------------------------------------

    // the scalars are solved separately from the NSE but in order to ensure consistency with
    // the continuity equation, the implicit terms (i.e., dependence of the scalars on rho, rhou)
    // are used on the RHS of the equations. This means that AScal[iScal][4] is the only implicit
    // term on the LHS, while AScal[iScal][0-3] are put back to the right side.
	
    for (int iScal = 0; iScal < nScal; iScal++)                               // prepare rhs and A
    {
    	string scalname = scalarTranspEqVector[iScal].getName();
    	for (int icv = 0; icv < ncv; ++icv) {
    		if (iScal == getScalarTransportIndex("omega") || iScal == getScalarTransportIndex("nuSA")
    		    || iScal == getScalarTransportIndex("eps")) residField2[icv] = rhsScal[iScal][icv];
    		if (iScal == getScalarTransportIndex("kine"))  residField3[icv] = rhsScal[iScal][icv];
    		if (iScal == getScalarTransportIndex("ZMean") || iScal == getScalarTransportIndex("v2")) residField4[icv] = rhsScal[iScal][icv];
    		if (iScal == getScalarTransportIndex("ZVar")  || iScal == getScalarTransportIndex("f"))  residField5[icv] = rhsScal[iScal][icv];
    		if (iScal == getScalarTransportIndex("CMean")) residField6[icv] = rhsScal[iScal][icv];

    		//rhsScal[iScal][icv] *= underRelax;
    		rhsScal[iScal][icv] *= scalarTranspEqVector[iScal].relax;

    		int noc_f = nbocv_i[icv];
    		int noc_l = nbocv_i[icv + 1] - 1;

    		double tmp = cv_volume[icv] / local_dt[icv];

    		if (iScal == getScalarTransportIndex("f")) {
    			// do nothing
    		} else {
    			AScal[iScal][5][noc_f] += tmp;                                 // diagonal part ( vol/dt + A )
    		}

    		// move the other implicit terms to the RHS
    		for (int noc = noc_f; noc <= noc_l; noc++)
    			rhsScal[iScal][icv] = rhsScal[iScal][icv]
    			                      - AScal[iScal][0][noc] * dq[nbocv_v[noc]][0]
    			                      - AScal[iScal][1][noc] * dq[nbocv_v[noc]][1]
    			                      - AScal[iScal][2][noc] * dq[nbocv_v[noc]][2]
    			                      - AScal[iScal][3][noc] * dq[nbocv_v[noc]][3]
    			                      - AScal[iScal][4][noc] * dq[nbocv_v[noc]][4];
    	}

    	solveLinSysScalar(dScal[iScal], AScal[iScal][5], rhsScal[iScal],
    			scalarTranspEqVector[iScal].phiZero,
    			scalarTranspEqVector[iScal].phiZeroRel,
    			scalarTranspEqVector[iScal].phiMaxiter,
    			scalarTranspEqVector[iScal].getName());

    	// update scalars and clip
    	double *phi = scalarTranspEqVector[iScal].phi;
    	for (int icv = 0; icv < ncv; icv++) {
    		//if(iScal == getScalarTransportIndex("f"))
    		//phi[icv] = min(max((phi[icv] + dScal[iScal][icv]), scalarTranspEqVector[iScal].lowerBound), scalarTranspEqVector[iScal].upperBound);
    		//else
    		phi[icv] = min(max((phi[icv]*(rho[icv]-dq[icv][0]) + dScal[iScal][icv])/rho[icv], scalarTranspEqVector[iScal].lowerBound), scalarTranspEqVector[iScal].upperBound);
    	}
    }

    // -------------------------------------------------------------------------------------------
    // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
    // -------------------------------------------------------------------------------------------
    updateCvDataG1G2(rho, REPLACE_DATA);
    updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
    updateCvDataG1G2(rhoE, REPLACE_DATA);

    for (int iScal = 0; iScal < nScal; iScal++)
      updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);



    // -------------------------------------------------------------------------------------------
    // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
    // -------------------------------------------------------------------------------------------
    calcStateVariables();

    // -------------------------------------------------------------------------------------------
    // update material properties: laminar viscosity and heat conductivity
    // -------------------------------------------------------------------------------------------
    calcMaterialProperties();

    // -------------------------------------------------------------------------------------------
    // set BC's for NS and scalars
    // -------------------------------------------------------------------------------------------
    setBC();

    // -------------------------------------------------------------------------------------------
    // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
    // -------------------------------------------------------------------------------------------
    calcRansTurbViscMuet();


    // =========================================================================================
    // calculate and show residual
    // =========================================================================================
    for (int i = 0; i < 5+nScal; i++)
    {
      myResidual[i] = 0.0;
      Residual[i] = 0.0;
    }

    for (int icv = 0; icv < ncv; icv++)
    {
      myResidual[0] += fabs(RHSrho[icv]);
      for (int i=0; i<3; i++)
        myResidual[i+1] += fabs(RHSrhou[icv][i]);
      myResidual[4] += fabs(RHSrhoE[icv]);
    }
    for (int iScal = 0; iScal < nScal; iScal++)
      for (int icv = 0; icv < ncv; icv++)
        myResidual[5+iScal] += fabs(rhsScal[iScal][icv]/underRelax);

    MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (step%check_interval == 0)
    {
      if ((mpi_rank == 0) && (step%(check_interval*10) == 0))
        cout << "\ndone step: "<< step << ", cfl: " << cfl << ", min. dt: " << dt_min << " time: " << time << endl;

      showResidue(Residual);
    }

    // IKJ: For output display (Metric to see the convergence of the simulation)
    for (int icv = 0; icv < ncv; icv++) {
    	log10_resid_rhoE[icv] = log10(fabs(RHSrhoE[icv]) + 1.0e-15);
    }
    for (int iScal = 0; iScal < nScal; iScal++)
    	for (int icv = 0; icv < ncv; icv++)
    		log10_resid_scalar0[icv] = log10(fabs(rhsScal[iScal][icv]/underRelax) + 1.0e-15);
    updateCvDataG1G2(residField,       REPLACE_DATA);
    updateCvDataG1G2(log10_resid_rhoE, REPLACE_DATA);
    if(nScal>0)
    	updateCvDataG1G2(log10_resid_scalar0, REPLACE_DATA);

    // IKJ
    double myTotResid_dq  = 0.0;
    double myTotResid_rhs = 0.0;
    for (int icv = 0; icv < ncv; icv++) {
    	for(int i=0; i<5; ++i)
    		myTotResid_dq += dq[icv][i];

    	myTotResid_rhs += fabs(RHSrho[icv]);
    	for (int i=0; i<3; i++)
    		myTotResid_rhs += fabs(RHSrhou[icv][i]);
    	myTotResid_rhs += fabs(RHSrhoE[icv]);
    }
	for (int iScal = 0; iScal < nScal; iScal++) {
		double *phi = scalarTranspEqVector[iScal].phi;
		for (int icv = 0; icv < ncv; icv++) {
			myTotResid_dq  += fabs( (dScal[iScal][icv] - phi[icv]*dq[icv][0]) / (rho[icv] - dq[icv][0]) );
			myTotResid_rhs += fabs(rhsScal[iScal][icv]);
		}
	}
    MPI_Allreduce(&myTotResid_dq,  &totResid_dq,  1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myTotResid_rhs, &totResid_rhs, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);


    temporalHook();
    dumpProbes(step, 0.0);
    writeData(step);

    if ((write_restart > 0) && (step % write_restart == 0))
      writeRestart(step);

    if ((step >= nsteps) || (Residual[4] <= resid_energ_th))   done = 1;
  }

  MPI_Barrier(mpi_comm);
  if (mpi_rank == 0)
  {
    double wtime0 = wtime;
    wtime = MPI_Wtime();
    cout << " > runtime for iterations[s]: " << wtime - wtime0 << endl;
  }


  // ---------------------------------------------------------------------------------
  // output
  // ---------------------------------------------------------------------------------

  temporalHook();
  finalHook();

  if( checkParam("STORE_JACOBIAN") ) {
	  // Write vectors
	  char filenameXcv[50];
	  sprintf(filenameXcv, "VectorXCVcoord.bin");
	  int dimension = 3;
	  double **xcvNew = NULL;
	  getMem2D(&xcvNew, 0, ncv-1, 0, dimension-1, "JoeWithModels::runBackwardEuler():xcvNew", true);
	  for(int icv=0; icv<ncv; ++icv) {
		  for(int i=0; i<dimension; ++i)
			  xcvNew[icv][i] = x_cv[icv][i];
	  }
	  MPI_Barrier(mpi_comm);
	  writeJoe2DvectorBinaryParallel(filenameXcv, xcvNew, dimension, cvora, ncv);
	  freeMem2D(xcvNew, 0, ncv-1, 0, dimension-1); 	xcvNew = NULL;

	  char filenameRhs[50];
	  sprintf(filenameRhs, "VectorRhsJoeBE.bin");
	  writeJoeUncoupledNSrhsBinaryParallel(filenameRhs, RHSrho, RHSrhou, RHSrhoE, cvora, ncv);

	  char filenameDiag[50];
	  sprintf(filenameDiag, "VectorDiagonalJoeBE.bin");
	  int nRowsPerBlock = 5;
	  double **diagPart = NULL;
	  getMem2D(&diagPart, 0, ncv-1, 0, nRowsPerBlock-1, "JoeWithModels::runBackwardEuler():diagPart", true);
	  for(int icv=0; icv<ncv; ++icv) {
		  double tmp = cv_volume[icv]/(local_dt[icv]);
		  for(int i=0; i<nRowsPerBlock; ++i)
			  diagPart[icv][i] = tmp;
	  }
	  MPI_Barrier(mpi_comm);
	  writeJoe2DvectorBinaryParallel(filenameDiag, diagPart, nRowsPerBlock, cvora, ncv);
	  freeMem2D(diagPart, 0, ncv-1, 0, nRowsPerBlock-1); 	diagPart = NULL;

	  // Write matrices
	  char filenameMatrix[50];
	  sprintf(filenameMatrix, "JacobianJoeBE.bin");
	  writeJoeUncoupledNSMatrixBinaryParallel(filenameMatrix, A, cvora, ncv, nbocv_i, nbocv_v_global);

	  char filenameMatrixClearDiag[50];
	  sprintf(filenameMatrixClearDiag, "JacobianClearDiagJoeBE.bin");
	  writeJoeUncoupledNSMatrixClearDiagBinaryParallel(filenameMatrixClearDiag, A, cvora, ncv, nbocv_i, nbocv_v_global, cv_volume, local_dt);

	  char filenameMatrixWeightedDiag[50];
	  sprintf(filenameMatrixWeightedDiag, "JacobianWeightedDiagJoeBE.bin");
	  writeJoeUncoupledNSMatrixWeightedDiagBinaryParallel(filenameMatrixWeightedDiag, A, cvora, ncv, nbocv_i, nbocv_v_global, cv_volume, local_dt);
  }

//  writeRestart();
  writeData(step);



  // ---------------------------------------------------------------------------------
  // delete memory
  // ---------------------------------------------------------------------------------

  delete [] A;
  delete [] rhs;

  delete [] RHSrho;
  delete [] RHSrhou;
  delete [] RHSrhoE;

  delete [] dq;

  delete [] myResidual;
  delete [] Residual; 	Residual = NULL;

  if (nScal > 0) freeMem3D(AScal,   0, nScal-1, 0, 5, 0, nbocv_s-1);
  if (nScal > 0) freeMem2D(dScal,   0, nScal-1, 0, ncv_g-1);
  if (nScal > 0) freeMem2D(rhsScal, 0, nScal-1, 0, ncv-1);
}


void JoeWithModels::runBDF2()
{
	int nScal = scalarTranspEqVector.size();

	double *myResidual = new double[5+nScal];
	Residual           = new double[5+nScal];
	double *fstResidual = new double[5+nScal];
	double *relResidual = new double[5+nScal];

	double *RHSrho = new double[ncv];
	double (*RHSrhou)[3] = new double[ncv][3];
	double *RHSrhoE = new double[ncv];

	double (*A)[5][5] = new double[nbocv_s][5][5];
	double (*dq)[5] = new double[ncv_g][5];
	double (*rhs)[5] = new double[ncv][5];
	double (*qn)[5] = new double[ncv][5];
	double (*qnm1)[5] = new double[ncv][5];

	double ***AScal   = NULL;  if (nScal > 0) getMem3D(&AScal,    0, nScal-1, 0, 5, 0, nbocv_s-1, "AScal");
	double **dScal    = NULL;  if (nScal > 0) getMem2D(&dScal,    0, nScal-1, 0, ncv_g-1, "dScal");
	double **rhsScal  = NULL;  if (nScal > 0) getMem2D(&rhsScal,  0, nScal-1, 0, ncv-1, "rhsScal");
	double **qnScal   = NULL;  if (nScal > 0) getMem2D(&qnScal,   0, nScal-1, 0, ncv-1, "qnScal");
	double **qnm1Scal = NULL;  if (nScal > 0) getMem2D(&qnm1Scal, 0, nScal-1, 0, ncv-1, "qnm1Scal");

	//------------------------------------
	// some parameters
	//------------------------------------
	double underRelax = getDoubleParam("UNDER_RELAXATION", "0.9");
	int bdf2_TVD_limiter = getIntParam("BDF2_TVD_LIMITER", "1");

	if (!checkParam("NEWTON_SOLVER_TRESHOLDS"))
	{
		ParamMap::add("NEWTON_SOLVER_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-6");    // add default values
		if (mpi_rank == 0)
			cout << "WARNING: added keyword \"NEWTON_SOLVER_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-6\"" <<
			" to parameter map!" << endl;
	}
	int    mnNewton       = getParam("NEWTON_SOLVER_TRESHOLDS")->getInt("MAX_ITER");
	double absResidNewton = getParam("NEWTON_SOLVER_TRESHOLDS")->getDouble("ABS_RESID");
	double relResidNewton = getParam("NEWTON_SOLVER_TRESHOLDS")->getDouble("REL_RESID");

	if (!checkParam("LINEAR_SOLVER_NS_TRESHOLDS"))
	{
		ParamMap::add("LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2");    // add default values
		if (mpi_rank == 0)
			cout << "WARNING: added keyword \"LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2\"" <<
			" to parameter map!" << endl;
	}
	int maxIterLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getInt("MAX_ITER");
	double zeroAbsLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("ABS_RESID");
	double zeroRelLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("REL_RESID");

	if (!checkParam("CFL_RAMP"))
	{
		ParamMap::add("CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0");    // add default: no increase of CFL number!
		if (mpi_rank == 0)
			cout << "WARNING: added \"CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0\" to parameter map" << endl;
	}

	int startIncCFL = getParam("CFL_RAMP")->getInt("AFTER_ITER");
	int intervalIncCFL = getParam("CFL_RAMP")->getInt("INTERVAL_ITER");
	double incCFL = getParam("CFL_RAMP")->getDouble("FACTOR_CFL");
	double maxCFL = getParam("CFL_RAMP")->getDouble("MAX_CFL");

	// -------------------------------------------------------------------------------------------
	// update state properties: velocity, pressure, temperature, enthalpy, gamma and R
	// -------------------------------------------------------------------------------------------
	calcStateVariables();

	// -------------------------------------------------------------------------------------------
	// update material properties: laminar viscosity and heat conductivity
	// -------------------------------------------------------------------------------------------
	calcMaterialProperties();

	// -------------------------------------------------------------------------------------------
	// set BC's for NS and scalars
	// -------------------------------------------------------------------------------------------
	setBC();

	// -------------------------------------------------------------------------------------------
	// update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
	// -------------------------------------------------------------------------------------------
	calcRansTurbViscMuet();


	// -------------------------------------------------------------------------------------------
	//
	//   Loop over time steps
	//
	// -------------------------------------------------------------------------------------------

	int done = 0;
	if ((nsteps >= 0)&&(step >= nsteps))  done = 1;

	if (initial_flowfield_output == "YES")
		writeData(0);

	// provide total runtime
	double wtime, wtime0;
	MPI_Barrier(mpi_comm);
	if (mpi_rank == 0)
		wtime = MPI_Wtime();



	//double bdf2Alfa = 1.5, bdf2Beta = 2.0, bdf2Gamma = 0.5;
	double bdf2Alfa = 0.5, bdf2Beta = 1.0, bdf2Gamma = 0.5;

	// store n and n-1 and solve n+1 using second order
	for (int icv=0; icv<ncv; icv++)
	{
		qn[icv][0] = rho[icv];
		qn[icv][1] = rhou[icv][0];
		qn[icv][2] = rhou[icv][1];
		qn[icv][3] = rhou[icv][2];
		qn[icv][4] = rhoE[icv];
	}

	for (int iScal = 0; iScal < nScal; iScal++)
	{
		double *phi = scalarTranspEqVector[iScal].phi;
		for (int icv=0; icv<ncv; icv++)
			qnScal[iScal][icv] = rho[icv]*phi[icv];
	}

	for(int i=0; i<5+nScal; ++i)
		Residual[i] = 1.0e20;

	while (done != 1)
	{
		step++;
		if ((step >= startIncCFL) && (step%intervalIncCFL == 0) && (cfl < maxCFL))      cfl *= incCFL;
		double dt_min = calcDt(cfl);
		time += dt_min;


		// store n and n-1 and solve n+1 using second order
		for (int icv=0; icv<ncv; icv++)
		{
			qnm1[icv][0] = qn[icv][0];
			qnm1[icv][1] = qn[icv][1];
			qnm1[icv][2] = qn[icv][2];
			qnm1[icv][3] = qn[icv][3];
			qnm1[icv][4] = qn[icv][4];

			qn[icv][0] = rho[icv];
			qn[icv][1] = rhou[icv][0];
			qn[icv][2] = rhou[icv][1];
			qn[icv][3] = rhou[icv][2];
			qn[icv][4] = rhoE[icv];
		}

		for (int iScal = 0; iScal < nScal; iScal++)
		{
			double *phi = scalarTranspEqVector[iScal].phi;

			for (int icv=0; icv<ncv; icv++)
			{
				qnm1Scal[iScal][icv] = qnScal[iScal][icv];
				qnScal[iScal][icv] = rho[icv]*phi[icv];
			}
		}


		// ---------------------------------------------------------------------------------
		// start Newton iterations
		// ---------------------------------------------------------------------------------

		if ((mpi_rank == 0) && (step%check_interval == 0))
			printf("dt: %.6le\n", dt_min);

		relResidual[4] = 1.0e20;

		int pN = 0;
		while ((pN < mnNewton) && (relResidual[4] > relResidNewton) && (pN == 0 || Residual[4] > absResidNewton))              // newton steps!!!
		{
			for (int i = 0; i < 5+nScal; i++)
				myResidual[i] = 0.0;

			for (int noc = 0; noc < nbocv_s; noc++)                             // set A, dq to zero! rhs is set zero in calcRHS
				for (int i = 0; i < 5; i++)
					for (int j = 0; j < 5; j++)
						A[noc][i][j] = 0.0;
			for (int icv = 0; icv < ncv_g; icv++)
				for (int i = 0; i < 5; i++)
					dq[icv][i] = 0.0;

			for (int iScal = 0; iScal < nScal; iScal++)                         // set AScal, dScal to zero! rhs is set zero in calcRHS
			{
				for (int i = 0; i <= 5; i++)
					for (int noc = 0; noc < nbocv_s; noc++)
						AScal[iScal][i][noc] = 0.0;
				for (int icv = 0; icv < ncv_g; icv++)
					dScal[iScal][icv] = 0.0;
			}


			// ---------------------------------------------------------------------------------
			// calculate RHS for both NSE and scalars
			// ---------------------------------------------------------------------------------

			calcRhs(RHSrho, RHSrhou, RHSrhoE, rhsScal, A, AScal, true);

			// ---------------------------------------------------------------------------------
			// solve linear system for the NSE
			// ---------------------------------------------------------------------------------

			for (int icv=0; icv<ncv; ++icv)                                 // prepare rhs and A
			{
				double tmp = cv_volume[icv]/(local_dt[icv]);

				double psi[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
				//double psi[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
				if (bdf2_TVD_limiter == 1) {
					bool tvdRzeroDemoninator = false;
					for(int i=0; i<5; ++i)
						if(fabs(qn[icv][i]-qnm1[icv][i]) < 1.0e-10) {
							tvdRzeroDemoninator = true;
							break;
						}

					if( !tvdRzeroDemoninator ) {
						double tvdR[5], minmod[5];

						tvdR[0] = (rho[icv]-    qn[icv][0])/(qn[icv][0]-qnm1[icv][0]);
						tvdR[1] = (rhou[icv][0]-qn[icv][1])/(qn[icv][1]-qnm1[icv][1]);
						tvdR[2] = (rhou[icv][1]-qn[icv][2])/(qn[icv][2]-qnm1[icv][2]);
						tvdR[3] = (rhou[icv][2]-qn[icv][3])/(qn[icv][3]-qnm1[icv][3]);
						tvdR[4] = (rhoE[icv]-   qn[icv][4])/(qn[icv][4]-qnm1[icv][4]);

						for (int i=0; i<5; i++) {
							if (tvdR[i] > 0.0) {
								if (1.0 < fabs(tvdR[i]))    minmod[i] = 1.0;
								else                        minmod[i] = tvdR[i];
							}
							else {
								minmod[i] = 0.0;
							}
							psi[i] = pow(minmod[i]/max(1.0, fabs(tvdR[i])), 0.5);
							if(isnan(psi[i]) || isinf(psi[i])) psi[i] = 1.0;
						}
					}
				}

				rhs[icv][0] = underRelax*(RHSrho[icv]     - ((1.0+psi[0]*bdf2Alfa)*rho[icv]     - (1.0+psi[0]*bdf2Beta)*qn[icv][0] + psi[0]*bdf2Gamma*qnm1[icv][0])*tmp);
				rhs[icv][1] = underRelax*(RHSrhou[icv][0] - ((1.0+psi[1]*bdf2Alfa)*rhou[icv][0] - (1.0+psi[1]*bdf2Beta)*qn[icv][1] + psi[1]*bdf2Gamma*qnm1[icv][1])*tmp);
				rhs[icv][2] = underRelax*(RHSrhou[icv][1] - ((1.0+psi[2]*bdf2Alfa)*rhou[icv][1] - (1.0+psi[2]*bdf2Beta)*qn[icv][2] + psi[2]*bdf2Gamma*qnm1[icv][2])*tmp);
				rhs[icv][3] = underRelax*(RHSrhou[icv][2] - ((1.0+psi[3]*bdf2Alfa)*rhou[icv][2] - (1.0+psi[3]*bdf2Beta)*qn[icv][3] + psi[3]*bdf2Gamma*qnm1[icv][3])*tmp);
				rhs[icv][4] = underRelax*(RHSrhoE[icv]    - ((1.0+psi[4]*bdf2Alfa)*rhoE[icv]    - (1.0+psi[4]*bdf2Beta)*qn[icv][4] + psi[4]*bdf2Gamma*qnm1[icv][4])*tmp);

				residField[icv] = RHSrhoE[icv];

				for (int i=0; i<5; i++)
					myResidual[i] += fabs(rhs[icv][i]);

				for (int i = 0; i < 5; i++)
					A[nbocv_i[icv]][i][i] += (1.0+psi[0]*bdf2Alfa)*tmp;                               // diagonal part ( vol/dt + A )
			}

			solveCoupledLinSysNS(dq, A, rhs, zeroAbsLS, zeroRelLS, maxIterLS);  // solve linear system

			for (int icv=0; icv<ncv; icv++)                                 // update solution: Q_new = Q_old + Delta_Q
			{
				rho[icv]     += dq[icv][0];
				rhou[icv][0] += dq[icv][1];
				rhou[icv][1] += dq[icv][2];
				rhou[icv][2] += dq[icv][3];
				rhoE[icv]    += dq[icv][4];
			}

			UpdateCvDataStateVec(dq);                                       // update dq since neighbors needed to compute RHS of scalars

			// ---------------------------------------------------------------------------------
			// solve linear system for the scalars
			// ---------------------------------------------------------------------------------

			// the scalars are solved separately from the NSE but in order to ensure consistency with
			// the continuity equation, the implicit terms (i.e., dependence of the scalars on rho, rhou)
			// are used on the RHS of the equations. This means that AScal[iScal][4] is the only implicit
			// term on the LHS, while AScal[iScal][0-3] are put back to the right side.

			for (int iScal = 0; iScal < nScal; iScal++)                               // prepare rhs and A
			{
				double *phi = scalarTranspEqVector[iScal].phi;

				for (int icv = 0; icv < ncv; ++icv)
				{
					double tmp = cv_volume[icv]/(local_dt[icv]);
					double rhoPhi = (rho[icv]-dq[icv][0])*phi[icv];

					double psi = 1.0;
					if (bdf2_TVD_limiter == 1)
					{
						// IKJ: We should prevent zero denominators!!
						bool tvdRzeroDemoninator = false;
						if(fabs(qnScal[iScal][icv]-qnm1Scal[iScal][icv]) < 1.0e-10)
							tvdRzeroDemoninator = true;

						if(!tvdRzeroDemoninator) {
							double minmod, tvdR = (rhoPhi-qnScal[iScal][icv])/(qnScal[iScal][icv]-qnm1Scal[iScal][icv]);
							if (tvdR > 0.0)
							{
								if (1.0 < fabs(tvdR))    minmod= 1.0;
								else                     minmod = tvdR;
							}
							else                       minmod = 0.0;

							psi = pow(minmod/max(1.0, fabs(tvdR)), 0.5); // Note: "psi" often becomes NaN
							if(isnan(psi) || isinf(psi)) psi = 1.0;
						}
					}

					rhsScal[iScal][icv] = underRelax*( rhsScal[iScal][icv]
					                                                  -( (1.0+psi*bdf2Alfa)*rhoPhi
					                                                		  - (1.0+psi*bdf2Beta)*qnScal[iScal][icv]
					                                                		                                     + psi*bdf2Gamma*qnm1Scal[iScal][icv]) * tmp );

					int noc_f = nbocv_i[icv];
					int noc_l = nbocv_i[icv + 1] - 1;

					AScal[iScal][5][noc_f] += (1.0+psi*bdf2Alfa)*tmp;                                 // diagonal part ( vol/dt + A )

					// move the other implicit terms to the RHS
					for (int noc = noc_f; noc <= noc_l; noc++)
						rhsScal[iScal][icv] = rhsScal[iScal][icv]
						                                     - AScal[iScal][0][noc] * dq[nbocv_v[noc]][0]
						                                                                               - AScal[iScal][1][noc] * dq[nbocv_v[noc]][1]
						                                                                                                                         - AScal[iScal][2][noc] * dq[nbocv_v[noc]][2]
						                                                                                                                                                                   - AScal[iScal][3][noc] * dq[nbocv_v[noc]][3]
						                                                                                                                                                                                                             - AScal[iScal][4][noc] * dq[nbocv_v[noc]][4];

					myResidual[5+iScal] += fabs(rhsScal[iScal][icv]);
				}

				solveLinSysScalar(dScal[iScal], AScal[iScal][5], rhsScal[iScal],
						scalarTranspEqVector[iScal].phiZero,
						scalarTranspEqVector[iScal].phiZeroRel,
						scalarTranspEqVector[iScal].phiMaxiter,
						scalarTranspEqVector[iScal].getName());

				// update scalars and clip
				for (int icv = 0; icv < ncv; icv++)
					phi[icv] = min(max((phi[icv]*(rho[icv]-dq[icv][0]) + dScal[iScal][icv])/rho[icv], scalarTranspEqVector[iScal].lowerBound), scalarTranspEqVector[iScal].upperBound);

			}

			// -------------------------------------------------------------------------------------------
			// update state properties: velocity, pressure, temperature, enthalpy, gamma and R
			// -------------------------------------------------------------------------------------------
			updateCvDataG1G2(rho, REPLACE_DATA);
			updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
			updateCvDataG1G2(rhoE, REPLACE_DATA);

			for (int iScal = 0; iScal < nScal; iScal++)
				updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);

			calcStateVariables();

			calcMaterialProperties();

			setBC();

			calcRansTurbViscMuet();

			MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);

			if (pN == 0)
				for (int i=0; i<5+nScal; i++)
					fstResidual[i] = Residual[i];

			for (int i=0; i<5+nScal; i++)
				relResidual[i] = Residual[i]/(fstResidual[i]+1.0e-10);

			if (mpi_rank == 0)
			{
				if (step%check_interval == 0)
				{
					Param *param;
					if (getParam(param, "LOGGING_LEVEL")) {
						string name = param->getString();

						bool isLowLogging = ( name.compare("EVERYTHING")==0 || name.compare("DEBUG_HI")==0 || name.compare("DEBUG_LO")==0 ||
								name.compare("INFO_HI")==0 || name.compare("INFO")==0 || name.compare("INFO_LO")==0 );
						        // Note: You can find LOGGING_LEVEL info is UgpWithCvCompFlow.h

						if (isLowLogging || (pN>0 && pN%std::max<int>(int((mnNewton-1)/2), 1)==0) || pN == mnNewton-1) {
							printf("Newton step residual: %d:\t", pN+1);
							for (int i=0; i<5+nScal; i++)
								printf("%.4le\t", relResidual[i]);
							printf("\n");
						}
					}
				}
			}

			pN++;
		}

		// =========================================================================================
		// calculate and show residual
		// =========================================================================================

		if (step%check_interval == 0)
		{
			if ((mpi_rank == 0) && (step%(check_interval*10) == 0))
				cout << "\ndone step: "<< step << ", cfl: " << cfl << ", min. dt: " << dt_min << " time: " << time << endl;

			showResidue(Residual);
		}

		// IKJ: For output display (Metric to see the convergence of the simulation)
		for (int icv = 0; icv < ncv; icv++) {
			residField[icv] = rhs[icv][4];
			log10_resid_rhoE[icv] = log10(fabs(RHSrhoE[icv]) + 1.0e-15);
		}
		for (int iScal = 0; iScal < nScal; iScal++)
			for (int icv = 0; icv < ncv; icv++)
				log10_resid_scalar0[icv] = log10(fabs(rhsScal[iScal][icv]/underRelax) + 1.0e-15);
		updateCvDataG1G2(residField,       REPLACE_DATA);
		updateCvDataG1G2(log10_resid_rhoE, REPLACE_DATA);
		if(nScal>0)
			updateCvDataG1G2(log10_resid_scalar0, REPLACE_DATA);

		// IKJ
		double myTotResid_rhs = 0.0;
		double myTotResid_dq  = 0.0;
		for (int icv = 0; icv < ncv; icv++) {
			myTotResid_rhs += fabs(RHSrho[icv]);
			for (int i=0; i<3; i++)
				myTotResid_rhs += fabs(RHSrhou[icv][i]);
			myTotResid_rhs += fabs(RHSrhoE[icv]);

			myTotResid_dq += fabs(rho[icv]- qn[icv][0]);
			for (int i=0; i<3; i++)
				myTotResid_dq += fabs(rhou[icv][i] - qn[icv][1+i]);
			myTotResid_dq += fabs(rhoE[icv] - qn[icv][4]);
		}
		for (int iScal = 0; iScal < nScal; iScal++) {
			double *phi = scalarTranspEqVector[iScal].phi;
			for (int icv = 0; icv < ncv; icv++) {
				myTotResid_rhs += fabs(rhsScal[iScal][icv]);
				myTotResid_dq  += fabs(phi[icv] - qnScal[iScal][icv]/qn[icv][0]);
			}
		}
		MPI_Allreduce(&myTotResid_dq,  &totResid_dq,  1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(&myTotResid_rhs, &totResid_rhs, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

		temporalHook();
		dumpProbes(step, 0.0);
		writeData(step);

		if ((write_restart > 0) && (step % write_restart == 0))
			writeRestart(step);

		if ((step >= nsteps) || (Residual[4] <= resid_energ_th))   done = 1;
	}

	MPI_Barrier(mpi_comm);
	if (mpi_rank == 0)
	{
		double wtime0 = wtime;
		wtime = MPI_Wtime();
		cout << " > runtime for iterations[s]: " << wtime - wtime0 << endl;
	}


	// ---------------------------------------------------------------------------------
	// output
	// ---------------------------------------------------------------------------------

	temporalHook();
	finalHook();

	writeRestart();


	// ---------------------------------------------------------------------------------
	// delete memory
	// ---------------------------------------------------------------------------------
	delete [] RHSrho;
	delete [] RHSrhou;
	delete [] RHSrhoE;

	delete [] myResidual;
	delete [] Residual; 	Residual = NULL;

	delete [] A;
	delete [] dq;
	delete [] rhs;
	delete [] qn;
	delete [] qnm1;

	if (nScal > 0) freeMem3D(AScal,   0, nScal-1, 0, 5, 0, nbocv_s-1);
	if (nScal > 0) freeMem2D(dScal,   0, nScal-1, 0, ncv_g-1);
	if (nScal > 0) freeMem2D(rhsScal, 0, nScal-1, 0, ncv-1);
	if (nScal > 0) freeMem2D(qnScal, 0, nScal-1, 0, ncv-1);
	if (nScal > 0) freeMem2D(qnm1Scal, 0, nScal-1, 0, ncv-1);

}


int JoeWithModels::calcRhs(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double **rhs_rhoScal, double (*A)[5][5], double ***AScal, int flagImplicit)
{
	static bool firstCall = true;

	firstCall_JOEcalcViscousFluxScalar = firstCall;

	// set RHS to zero
	for (int icv = 0; icv < ncv; icv++)
	{
		rhs_rho[icv] = 0.0;
		for (int i = 0; i < 3; i++)
			rhs_rhou[icv][i] = 0.0;
		rhs_rhoE[icv] = 0.0;
	}

	// set scalars RHS to zero
	for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
		for (int icv = 0; icv < ncv; icv++)
			rhs_rhoScal[iScal][icv] = 0.0;

	// =======================================================================================
	// NAVIER-STOKES
	// =======================================================================================
	// compute Euler Flux for NS and scalars
	int CountReducedOrder = calcEulerFlux(rhs_rho, rhs_rhou, rhs_rhoE, rhs_rhoScal, A, AScal, flagImplicit);

	// compute viscous Flux for NS
#ifdef USE_ARTIF_VISC
	if(UgpWithCvCompFlow::turnOnArtifVisc) {
		UgpWithCvCompFlow::calcArtifVisc(artifVisc_mag,
				artifVisc_bulkViscOnly, artifVisc_shockOnly,
				artifVisc_smoothstepThresh, artifVisc_type, artifVisc_coeff);
	}
	if (mu_ref > 0.0 || UgpWithCvCompFlow::turnOnArtifVisc)
		calcViscousFluxNS(rhs_rho, rhs_rhou, rhs_rhoE, A, flagImplicit);
#else
	if (mu_ref > 0.0)
		calcViscousFluxNS(rhs_rho, rhs_rhou, rhs_rhoE, A, flagImplicit);
#endif

	// add source terms to RHS of Navier-Stokes equations
	sourceHook(rhs_rho, rhs_rhou, rhs_rhoE, A);
	sourceHookRansTurb(rhs_rho, rhs_rhou, rhs_rhoE, A);
	sourceHookRansComb(rhs_rho, rhs_rhou, rhs_rhoE, A);

	// =======================================================================================
	// SCALARS
	// =======================================================================================
	for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
	{
		// compute viscous Flux for scalars and add source terms to RHS of scalar equations
		if (AScal == NULL) {
			if (mu_ref > 0.0 || d_scalar_comb_ref > 0.0)
				calcViscousFluxScalar(rhs_rhoScal[iScal], NULL, scalarTranspEqVector[iScal], flagImplicit);

			sourceHookScalarRansTurb(rhs_rhoScal[iScal], NULL, scalarTranspEqVector[iScal].getName(), flagImplicit);
			sourceHookScalarRansComb(rhs_rhoScal[iScal], NULL, scalarTranspEqVector[iScal].getName(), flagImplicit);
		} else {
			if (mu_ref > 0.0 || d_scalar_comb_ref > 0.0)
				calcViscousFluxScalar(rhs_rhoScal[iScal], AScal[iScal][5], scalarTranspEqVector[iScal], flagImplicit);

			sourceHookScalarRansTurb(rhs_rhoScal[iScal], AScal[iScal][5], scalarTranspEqVector[iScal].getName(), flagImplicit);
			sourceHookScalarRansComb(rhs_rhoScal[iScal], AScal[iScal][5], scalarTranspEqVector[iScal].getName(), flagImplicit);
		}
	}

	firstCall = false;

	return CountReducedOrder;
}

int JoeWithModels::calcEulerFlux(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double **rhs_rhoScal, double (*A)[5][5], double ***AScal, int flagImplicit)
{
	int nScal = scalarTranspEqVector.size();

	double (*Apl)[5] = NULL;
	double (*Ami)[5] = NULL;

	if (flagImplicit)
	{
		Apl = new double[5][5];
		Ami = new double[5][5];

		for (int i = 0; i < 5; i++)
			for (int j = 0; j < 5; j++)
				Apl[i][j] = Ami[i][j] = 0.0;
	}

	// Implicit matrix for scalars: the definition is not easy, so ask an expert (e.g., R. Pecnik)
	double (*AplScal)[6] = NULL;
	double (*AmiScal)[6] = NULL;

	if (flagImplicit)
	{
		if (nScal > 0) AplScal = new double[nScal][6];
		if (nScal > 0) AmiScal = new double[nScal][6];

		for (int iScal = 0; iScal < nScal; iScal++)
			for (int i = 0; i <= 5; i++)
				AplScal[iScal][i] = AmiScal[iScal][i] = 0.0;
	}

	double Frho, Frhou[3], FrhoE;
	double *FrhoScal     = NULL;         if (nScal > 0) FrhoScal     = new double[nScal];
	double *Scalar0      = NULL;         if (nScal > 0) Scalar0      = new double[nScal];            // cell face if second order
	double *Scalar1      = NULL;         if (nScal > 0) Scalar1      = new double[nScal];
	double *ScalCV0      = NULL;         if (nScal > 0) ScalCV0      = new double[nScal];            // cell face if second order
	double *ScalCV1      = NULL;         if (nScal > 0) ScalCV1      = new double[nScal];
	double *ScalConvTerm = NULL;         if (nScal > 0) ScalConvTerm = new double[nScal];            // 0 if convective term no considered, otherwise 1


	// count how many cells switched back to first order due to extrapolated negative density or pressure/temperature at the faces
	int CountReducedOrder = 0;
	int myCountReducedOrder = 0;

	// save the index of kine if defined and save convTerm for speed
	int kine_Index = getScalarTransportIndex("kine");

	for (int iScal = 0; iScal < nScal; iScal++)
		ScalConvTerm[iScal] = double(scalarTranspEqVector[iScal].convTerm);

	if (sndOrder == true || mu_ref>0.0)
		calcCv2Grad(grad_u,   vel, limiterNavierS, sos, epsilonSDWLS);
	// =============================================================================================
	// compute gradients, with boundary values
	// =============================================================================================
	if (sndOrder == true)
	{
		calcCv2Grad(grad_rho, rho, limiterNavierS, rho, epsilonSDWLS);
#ifdef temp_reconstruction
		calcCv2Grad(grad_temp, temp, limiterNavierS, temp, epsilonSDWLS);
#else
		calcCv2Grad(grad_p, press, limiterNavierS, press, epsilonSDWLS);
#endif

		for (int iScal = 0; iScal < nScal; iScal++)
		{
			if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
			{
				// For the scalar reconstruction, Grad(rho*Phi) is required
				// Gradients of rho*Phi are saved in grad_phi temporarily
				// Gradients of rho*Phi are also limited like rho with alpha_rho
				// Boundary face values rho_bfa*Phi_fa are saved in Phi_fa temporarily
				double *rhoPhi = new double[ncv_ggff];
				double *phi = scalarTranspEqVector[iScal].phi;
				double (*grad_phi)[3] = scalarTranspEqVector[iScal].grad_phi;

				// Compute rho*Phi
				for (int icv = 0; icv < ncv_ggff; icv++)
					rhoPhi[icv] = rho[icv] * phi[icv];

				// Compute gradients of rho*Phi and limit based on rho*Phi
				calcCv2Grad(grad_phi, rhoPhi, limiterNavierS, rhoPhi, epsilonSDWLS);

				delete [] rhoPhi;
			}
			else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
			{
				double *phi = scalarTranspEqVector[iScal].phi;
				double (*grad_phi)[3] = scalarTranspEqVector[iScal].grad_phi;
				calcCv2Grad(grad_phi, phi, limiterNavierS, phi, epsilonSDWLS);
			}

			else
			{
				cerr << "### JoeWithModels::calcEulerFlux_new => Wrong reconstruction type for scalars! ###" << endl;
				throw(-1);
			}
		}
	}

	// ===============================================================================================
	// cycle through internal faces, assembling flux to both sides
	// ===============================================================================================
	for (int ifa = nfa_b; ifa < nfa; ifa++)
	{
		int icv0 = cvofa[ifa][0];
		int icv1 = cvofa[ifa][1];
		assert( icv0 >= 0 );
		assert( icv1 >= 0 );

		int noc00, noc01, noc11, noc10;
		if (flagImplicit)
			getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

		// face unit normal and area...
		double nVec[3] = {0.0, 0.0, 0.0};
		double area = normVec3d(nVec, fa_normal[ifa]);

		// .............................................................................................
		// reconstruction of variables at faces: rho, u, T or P, scalars
		// .............................................................................................
		double rho0 = rho[icv0];
		double u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
		double p0 = press[icv0];
		double T0 = temp[icv0];
		double h0 = enthalpy[icv0];
		double gam0 = gamma[icv0];
		double R0 = RoM[icv0];
		double kineCV0 = 0.0;          // cell center
		double kineFA0 = 0.0;          // cell face if second order

		double rho1 = rho[icv1];
		double u1[3] = {vel[icv1][0], vel[icv1][1], vel[icv1][2]};
		double p1 = press[icv1];
		double T1 = temp[icv1];
		double h1 = enthalpy[icv1];
		double gam1 = gamma[icv1];
		double R1 = RoM[icv1];
		double kineCV1 = 0.0;
		double kineFA1 = 0.0;

		for (int iScal = 0; iScal < nScal; iScal++)
		{
			double *phi = scalarTranspEqVector[iScal].phi;
			ScalCV0[iScal] = Scalar0[iScal] = phi[icv0];
			ScalCV1[iScal] = Scalar1[iScal] = phi[icv1];
		}

		if (sndOrder == true)
		{
			double r0[3] = {0.0, 0.0, 0.0}, r1[3] = {0.0, 0.0, 0.0};
			vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);
			vecMinVec3d(r1, x_fa[ifa], x_cv[icv1]);

			// ----------------------------------------
			// left side
			// ----------------------------------------
			rho0 += vecDotVec3d(r0, grad_rho[icv0]);
#ifdef temp_reconstruction
			T0 += vecDotVec3d(r0, grad_temp[icv0]);
			if ((T0 <= 0.0) || (rho0 <= 0.0))
			{
				T0 = temp[icv0];
				rho0 = rho[icv0];
				myCountReducedOrder++;
			}
#else
			p0 += vecDotVec3d(r0, grad_p[icv0]);
			if ((p0 <= 0.0) || (rho0 <= 0.0))
			{
				p0 = press[icv0];
				rho0 = rho[icv0];
				myCountReducedOrder++;
			}
#endif
			else
			{
				for (int i = 0; i < 3; i++)
					u0[i] += vecDotVec3d(r0, grad_u[icv0][i]);
				for (int iScal = 0; iScal < nScal; iScal++)
				{
					double (*grad_phi)[3] = scalarTranspEqVector[iScal].grad_phi;

					if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
						Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d(r0, grad_phi[icv0])) / rho0;
					else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
						Scalar0[iScal] += vecDotVec3d(r0, grad_phi[icv0]);
				}
			}

			// ----------------------------------------
			// right side
			// ----------------------------------------
			rho1 += vecDotVec3d(r1, grad_rho[icv1]);
			if(isNaN0(rho1)) {
				printf("ERROR in JoeWithModels::calcEulerFlux(): NaN detected in rho[icv1] at INTERNAL ifa=%d(%.2e,%.2e,%.2e)\n", ifa, x_fa[ifa][0], x_fa[ifa][1], x_fa[ifa][2]);
				printf("         Details: mpi_rank=%d, icv0=%d(%.2e,%.2e,%.2e), icv1=%d(%.2e,%.2e,%.2e)\n", mpi_rank, icv0, x_cv[icv0][0],x_cv[icv0][1],x_cv[icv0][2], icv1, x_cv[icv1][0],x_cv[icv1][1],x_cv[icv1][2]);
				printf("                  grad_rho[icv1]=(%.3e, %.3e, %.3e)\n",grad_rho[icv1][0],grad_rho[icv1][1],grad_rho[icv1][2]);
				throw(-11);
			}
#ifdef temp_reconstruction
			T1 += vecDotVec3d(r1, grad_temp[icv1]);
			if ((T1 <= 0.0) || (rho1 <= 0.0))
			{
				T1 = temp[icv1];
				rho1 = rho[icv1];
				myCountReducedOrder++;
			}
#else
			p1 += vecDotVec3d(r1, grad_p[icv1]);
			if ((p1 <= 0.0) || (rho1 <= 0.0))
			{
				p1 = press[icv1];
				rho1 = rho[icv1];
				myCountReducedOrder++;
			}
#endif
			else
			{
				for (int i = 0; i < 3; i++)
					u1[i] += vecDotVec3d(r1, grad_u[icv1][i]);
				for (int iScal = 0; iScal < nScal; iScal++)
				{
					double (*grad_phi)[3] = scalarTranspEqVector[iScal].grad_phi;

					if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
						Scalar1[iScal] = (rho[icv1] * Scalar1[iScal] + vecDotVec3d(r1, grad_phi[icv1])) / rho1;

					else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
						Scalar1[iScal] += vecDotVec3d(r1, grad_phi[icv1]);


				}
			}


			// .............................................................................................
			// calculation of other variables at faces: p/T, h, R, gam
			// .............................................................................................
#ifdef temp_reconstruction
			calcThermoProp_T(p0, h0, R0, gam0, rho0, T0, Scalar0, nScal);
			calcThermoProp_T(p1, h1, R1, gam1, rho1, T1, Scalar1, nScal);
#else
			calcThermoProp_p(T0, h0, R0, gam0, rho0, p0, Scalar0, nScal);
			calcThermoProp_p(T1, h1, R1, gam1, rho1, p1, Scalar1, nScal);
#endif
		}

		if (kine_Index > -1)   // save kine if defined
		{
			kineCV0 = ScalCV0[kine_Index];         // cell center left for implicit side (Jacobi is computed with cell center)
			kineFA0 = Scalar0[kine_Index];         // cell face left, if second order
			kineCV1 = ScalCV1[kine_Index];         // cell center left for implicit side (Jacobi is computed with cell center)
			kineFA1 = Scalar1[kine_Index];         // cell face right,
		}

		// .............................................................................................
		// calculation of Euler Flux explicit using HLLC
		// .............................................................................................
		calcEulerFlux_HLLC(Frho, Frhou, FrhoE, FrhoScal,
				rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
				rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
				area, nVec, nScal, 0.0);

		if(isNaN0(Frho) || isNaN0(Frhou[0]) || isNaN0(Frhou[1]) || isNaN0(Frhou[2]) || isNaN0(FrhoE)) {
			printf("ERROR in JoeWithModels::calcEulerFlux(): NaN detected in flux at INTERNAL ifa=%d(%e,%e,%e)\n", ifa, x_fa[ifa][0], x_fa[ifa][1], x_fa[ifa][2]);
			printf("         Details: mpi_rank=%d, icv0=%d(%.4e,%.4e,%.4e), icv1=%d(%.4e,%.4e,%.4e)\n", mpi_rank, icv0, x_cv[icv0][0],x_cv[icv0][1],x_cv[icv0][2], icv1, x_cv[icv1][0],x_cv[icv1][1],x_cv[icv1][2]);
			printf("                  Frho=%.3e, Frhou=(%.3e,%.3e,%.3e), FrhoE=%3e\n", Frho, Frhou[0],Frhou[1],Frhou[2], FrhoE);
			printf("                  rho0=%.3e, u0=(%.3e, %.3e, %.3e), p0=%.3e, T0=%.3e, h0=%.3e, R0=%.3e, gam0=%.3e, kineFA0=%.3e", rho0, u0[0], u0[1], u0[2], p0, T0, h0, R0, gam0, kineFA0);
			for (ScalarTranspEqIterator data = scalarTranspEqVector.begin(); data < scalarTranspEqVector.end(); data++)
				printf(", %s=%.3e", data->getName(), (data->phi)[icv0]);
			printf("\n");
			printf("                  rho1=%.3e, u1=(%.3e, %.3e, %.3e), p1=%.3e, T1=%.3e, h1=%.3e, R1=%.3e, gam1=%.3e, kineFA1=%.3e", rho1, u1[0], u1[1], u1[2], p1, T1, h1, R1, gam1, kineFA1);
			for (ScalarTranspEqIterator data = scalarTranspEqVector.begin(); data < scalarTranspEqVector.end(); data++)
				printf(", %s=%.3e", data->getName(), (data->phi)[icv1]);
			printf("\n\n");
			throw(-11);
		}

		// icv0 is always valid...
		rhs_rho[icv0] -= Frho;
		for (int i = 0; i < 3; i++)
			rhs_rhou[icv0][i] -= Frhou[i];
		rhs_rhoE[icv0] -= FrhoE;

		for (int iScal = 0; iScal < nScal; iScal++)
			rhs_rhoScal[iScal][icv0] -= ScalConvTerm[iScal] * FrhoScal[iScal];

		// icv1 can be ghost...
		if (icv1 < ncv)
		{
			rhs_rho[icv1] += Frho;
			for (int i = 0; i < 3; i++)
				rhs_rhou[icv1][i] += Frhou[i];
			rhs_rhoE[icv1] += FrhoE;
			for (int iScal = 0; iScal < nScal; iScal++)
				rhs_rhoScal[iScal][icv1] += ScalConvTerm[iScal] * FrhoScal[iScal];
		}

		// .............................................................................................
		// calculate implicit matrix using HLLC
		// .............................................................................................
		if (flagImplicit)
		{
			calcEulerFluxMatrices_HLLC(Apl, Ami, AplScal, AmiScal,
					rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, kineCV0,
					rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, kineCV1,
					area, nVec, nScal, 0.0);

			for (int i = 0; i < 5; i++)
				for (int j = 0; j < 5; j++)
				{
					A[noc00][i][j] += Apl[i][j];
					A[noc01][i][j] += Ami[i][j];
				}

			if (icv1 < ncv)  // if icv1 is internal...
					{
				for (int i = 0; i < 5; i++)
					for (int j = 0; j < 5; j++)
					{
						A[noc11][i][j] -= Ami[i][j];
						A[noc10][i][j] -= Apl[i][j];
					}
					}

			for (int iScal = 0; iScal < nScal; iScal++)
				for (int i = 0; i <= 5; i++)
				{
					AScal[iScal][i][noc00] += ScalConvTerm[iScal] * AplScal[iScal][i];
					AScal[iScal][i][noc01] += ScalConvTerm[iScal] * AmiScal[iScal][i];

					if (icv1 < ncv)
					{
						AScal[iScal][i][noc11] -= ScalConvTerm[iScal] * AmiScal[iScal][i];
						AScal[iScal][i][noc10] -= ScalConvTerm[iScal] * AplScal[iScal][i];
					}
				}
		}
	}

	// ===============================================================================================
	// cycle through boundary faces, assembling flux
	// ===============================================================================================
	for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
	{
		if (zone->getKind() == FA_ZONE_BOUNDARY)
		{
			Param *param;

			if (getParam(param, zone->getName()))
			{
				// .............................................................................................
				// SYMMETRY BOUNDARY CONDITION OR WALL BOUNDARY CONDITION, ATTENTION: ONLY FIRST ORDER!!!
				// .............................................................................................
				if ((param->getString() == "SYMMETRY") || (param->getString() == "WALL"))
				{
					for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ++ifa) {
						/*
	    for (int index = 0; index < zone->faVec.size(); ++index) {
	      int ifa = zone->faVec[index];
						 */
						int icv0 = cvofa[ifa][0];
						int icv1 = cvofa[ifa][1];

						double nVec[3] = {0.0, 0.0, 0.0};
						double area = normVec3d(nVec, fa_normal[ifa]);

						for (int iScal = 0; iScal < nScal; iScal++)
							Scalar0[iScal] = scalarTranspEqVector[iScal].phi[icv1];

						double kineFA = 0.0;
						if (kine_Index > -1)
							kineFA = scalarTranspEqVector[kine_Index].phi[icv1];

						calcEulerFlux_HLLC(Frho, Frhou, FrhoE, FrhoScal,
								rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
								rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
								area, nVec, nScal, 0.0);

						rhs_rho[icv0] -= Frho;
						for (int i = 0; i < 3; i++)
							rhs_rhou[icv0][i] -= Frhou[i];
						rhs_rhoE[icv0] -= FrhoE;
						for (int iScal = 0; iScal < nScal; iScal++)
							rhs_rhoScal[iScal][icv0] -= ScalConvTerm[iScal] * FrhoScal[iScal];

						if(isNaN0(Frho) || isNaN0(Frhou[0]) || isNaN0(Frhou[1]) || isNaN0(Frhou[2]) || isNaN0(FrhoE)) {
							printf("ERROR in JoeWithModels::calcEulerFlux(): NaN detected in flux at BOUNDARY (SYMM or WALL) ifa=%d(%.2e,%.2e,%.2e)\n", ifa, x_fa[ifa][0], x_fa[ifa][1], x_fa[ifa][2]);
							printf("         Details: mpi_rank=%d, icv0=%d(%.2e,%.2e,%.2e), icv1=%d(%.2e,%.2e,%.2e)\n", mpi_rank, icv0, x_cv[icv0][0],x_cv[icv0][1],x_cv[icv0][2], icv1, x_cv[icv1][0],x_cv[icv1][1],x_cv[icv1][2]);
							printf("                  Frho=%.3e, Frhou=(%.3e,%.3e,%.3e), FrhoE=%.3e\n", Frho, Frhou[0],Frhou[1],Frhou[2], FrhoE);
							printf("                  rho=%.3e, vel=(%.3e, %.3e, %.3e), press=%.3e, temp=%.3e, enthalpy=%.3e, RoM=%.3e, gamma=%.3e, kineFA=%.3e\n", rho[icv1], vel[icv1][0], vel[icv1][1], vel[icv1][2], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], kineFA);
							throw(-11);
						}

						if (flagImplicit)
						{
							calcEulerFluxMatrices_HLLC(Apl, NULL, AplScal, AmiScal,
									rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
									rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
									area, nVec, nScal, 0.0);

							int noc00 = nbocv_i[icv0]; // icv0's diagonal
							for (int i = 0; i < 5; i++)
								for (int j = 0; j < 5; j++)
									A[noc00][i][j] += Apl[i][j];

							for (int iScal = 0; iScal < nScal; iScal++)
								for (int i = 0; i <= 5; i++)
									AScal[iScal][i][noc00] += ScalConvTerm[iScal] * AplScal[iScal][i];
						}
					}
				}
				// .............................................................................................
				// OTHER BOUNDARY CONDITIONS (HOOK, DIRICHLET(CBC), NEUMANN, ...)
				// .............................................................................................
				else
				{

					for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ++ifa) {
						/*
	      for (int index = 0; index < zone->faVec.size(); ++index) {
	      int ifa = zone->faVec[index];
						 */
						int icv0 = cvofa[ifa][0];
						int icv1 = cvofa[ifa][1];

						double nVec[3] = {0.0, 0.0, 0.0};
						double area = normVec3d(nVec, fa_normal[ifa]);

						double rho0 = rho[icv0];
						double u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
						double p0 = press[icv0];
						double T0 = temp[icv0];
						double h0 = enthalpy[icv0];
						double gam0 = gamma[icv0];
						double R0 = RoM[icv0];
						double kineCV0 = 0.0;           // cell center
						double kineFA0 = 0.0;           // cell face

						double kineFA1 = 0.0;

						for (int iScal = 0; iScal < nScal; iScal++)
						{
							ScalCV0[iScal] = Scalar0[iScal] = scalarTranspEqVector[iScal].phi[icv0];
							Scalar1[iScal] = scalarTranspEqVector[iScal].phi[icv1];
						}

						if (sndOrder == true)
						{
							double r0[3] = {0.0, 0.0, 0.0};
							vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);

							// left side
							rho0 += vecDotVec3d(r0, grad_rho[icv0]);
#ifdef temp_reconstruction
							T0 += vecDotVec3d(r0, grad_temp[icv0]);
							if ((T0 <= 0.0) || (rho0 <= 0.0))
							{
								T0 = temp[icv0];
								rho0 = rho[icv0];
								myCountReducedOrder++;
							}
#else
							p0 += vecDotVec3d(r0, grad_p[icv0]);
							if ((p0 <= 0.0) || (rho0 <= 0.0))
							{
								p0 = press[icv0];
								rho0 = rho[icv0];
								myCountReducedOrder++;
							}
#endif
							else
							{
								for (int i = 0; i < 3; i++)
									u0[i] += vecDotVec3d(r0, grad_u[icv0][i]);
								for (int iScal = 0; iScal < nScal; iScal++)
								{
									double (*grad_phi)[3] = scalarTranspEqVector[iScal].grad_phi;
									if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
										Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d(r0, grad_phi[icv0])) / rho0;
									else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
										Scalar0[iScal] += vecDotVec3d(r0, grad_phi[icv0]);
								}
							}

							// calculation of other variables at faces: p/T, h, R, gam
#ifdef temp_reconstruction
							calcThermoProp_T(p0, h0, R0, gam0, rho0, T0, Scalar0, nScal);
#else
							calcThermoProp_p(T0, h0, R0, gam0, rho0, p0, Scalar0, nScal);
#endif
						}

						if (kine_Index > -1)   // save kine if defined
						{
							kineCV0 = ScalCV0[kine_Index];          // cell center
							kineFA0 = Scalar0[kine_Index];                                 // cell face
							kineFA1 = Scalar1[kine_Index];
						}

						if (rho0 <= 0.0) {
							cout << "got negative rho0: " << rho0 << endl;
							throw(-1);
						}
						if (rho[icv1] <= 0.0) {
							cout << "got negative rho[icv1]: " << rho[icv1] << " at icv1 = "<<icv1<<" ("<<x_cv[icv1][0]<<", "<<x_cv[icv1][1]<<", "<<x_cv[icv1][2]<<")"<<endl;
							throw(-1);
						}

						calcEulerFlux_HLLC(Frho, Frhou, FrhoE, FrhoScal,
								rho0,      u0,        p0,          T0,         h0,             R0,        gam0,        Scalar0, kineFA0,
								rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar1, kineFA1,
								area, nVec, nScal, 0.0);

						rhs_rho[icv0] -= Frho;
						for (int i = 0; i < 3; i++)
							rhs_rhou[icv0][i] -= Frhou[i];
						rhs_rhoE[icv0] -= FrhoE;
						for (int iScal = 0; iScal < nScal; iScal++)
							rhs_rhoScal[iScal][icv0] -= ScalConvTerm[iScal] * FrhoScal[iScal];

						if(isNaN0(Frho) || isNaN0(Frhou[0]) || isNaN0(Frhou[1]) || isNaN0(Frhou[2]) || isNaN0(FrhoE)) {
							printf("ERROR in JoeWithModels::calcEulerFlux(): NaN detected in flux at BOUNDARY (HOOK, DIRI(CBC), NEUM) ifa=%d(%.2e,%.2e,%.2e)\n", ifa, x_fa[ifa][0], x_fa[ifa][1], x_fa[ifa][2]);
							printf("         Details: mpi_rank=%d, icv0=%d(%.2e,%.2e,%.2e), icv1=%d(%.2e,%.2e,%.2e)\n", mpi_rank, icv0, x_cv[icv0][0],x_cv[icv0][1],x_cv[icv0][2], icv1, x_cv[icv1][0],x_cv[icv1][1],x_cv[icv1][2]);
							printf("                  Frho=%.3e, Frhou=(%.3e,%.3e,%.3e), FrhoE=%.3e\n", Frho, Frhou[0],Frhou[1],Frhou[2], FrhoE);
							printf("                  rho0=%.3e, u0=(%.3e, %.3e, %.3e), p0=%.3e, T0=%.3e, h0=%.3e, R0=%.3e, gam0=%.3e, kineFA0=%.3e\n", rho0, u0[0], u0[1], u0[2], p0, T0, h0, R0, gam0, kineFA0);
							printf("                  rho=%.3e, vel=(%.3e, %.3e, %.3e), press=%.3e, temp=%.3e, enthalpy=%.3e, RoM=%.3e, gamma=%.3e, kineFA=%.3e\n", rho[icv1], vel[icv1][0], vel[icv1][1], vel[icv1][2], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], kineFA1);
							throw(-11);
						}

						if (flagImplicit)
						{
							calcEulerFluxMatrices_HLLC(Apl, NULL, AplScal, AmiScal,
									rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, kineCV0,
									rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar1, kineFA1,
									area, nVec, nScal, 0.0);

							int noc00 = nbocv_i[icv0]; // icv0's diagonal
							for (int i = 0; i < 5; i++)
								for (int j = 0; j < 5; j++)
									A[noc00][i][j] += Apl[i][j];

							for (int iScal = 0; iScal < nScal; iScal++)
								for (int i = 0; i <= 5; i++)
									AScal[iScal][i][noc00] += ScalConvTerm[iScal] * AplScal[iScal][i];
						}
					}
				}
			}
		}
	}

	// output the number of times switched back to first order at faces
	MPI_Allreduce(&myCountReducedOrder, &CountReducedOrder, 1, MPI_INT, MPI_SUM, mpi_comm);
	if ((CountReducedOrder > 0) && (mpi_rank == 0))
		cout << "Switched back to first order at " << CountReducedOrder << " face(s)" << endl;

	if (Apl != NULL)  delete [] Apl;
	if (Ami != NULL)  delete [] Ami;

	if (AplScal != NULL) delete [] AplScal;
	if (AmiScal != NULL) delete [] AmiScal;

	if (nScal > 0) delete [] FrhoScal;
	if (nScal > 0) delete [] Scalar0;
	if (nScal > 0) delete [] Scalar1;
	if (nScal > 0) delete [] ScalCV0;
	if (nScal > 0) delete [] ScalCV1;
	if (nScal > 0) delete [] ScalConvTerm;

	return CountReducedOrder;
}

void JoeWithModels::calcViscousFluxNS(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5], int flagImplicit)
{
	double (*A0)[5];
	double (*A1)[5];

	if (flagImplicit)
	{
		A0 = new double[5][5];
		A1 = new double[5][5];

		for (int i = 0; i < 5; i++)
			for (int j = 0; j < 5; j++)
				A0[i][j] = A1[i][j] = 0.0;
	}
	else
		A0 = A1 = NULL;

	double Frhou[3] = {0.0, 0.0, 0.0}, FrhoE = 0.0;

	// save the index of kine if defined
	int kine_index = -1;
	for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
		if (strcmp(scalarTranspEqVector[iScal].getName(),"kine")==0)
			kine_index = iScal;

	// ====================================================================
	// compute gradients
	// ====================================================================
	if (sndOrder != true)
		calcCv2Grad(grad_u, vel, limiterNavierS, sos, epsilonSDWLS);

	calcCv2Grad(grad_enthalpy, enthalpy, limiterNavierS, enthalpy, epsilonSDWLS);

	// ====================================================================
	// cycle through internal faces, assembling flux and matrix
	// ====================================================================
	for (int ifa = nfa_b; ifa < nfa; ifa++)
	{
		int icv0 = cvofa[ifa][0];
		int icv1 = cvofa[ifa][1];

		int noc00, noc01, noc11, noc10;
		if (flagImplicit)
			getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

		// face unit normal and area...
		double nVec[3] = {0.0, 0.0, 0.0};
		double area = normVec3d(nVec, fa_normal[ifa]);
		double sVec[3] = {0.0, 0.0, 0.0};
		vecMinVec3d(sVec, x_cv[icv1], x_cv[icv0]);
		double smag = normVec3d(sVec);

		double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
		vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
		vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
		double w0 = sqrt(vecDotVec3d(dx0, dx0));
		double w1 = sqrt(vecDotVec3d(dx1, dx1));
		double ws = w0 + w1;
		w0 /= ws;
		w1 /= ws;
		double uAux_fa[3] = { w1*vel[icv0][0]+ w0*vel[icv1][0],
				w1*vel[icv0][1]+ w0*vel[icv1][1],
				w1*vel[icv0][2]+ w0*vel[icv1][2]};

		// kinetic energy, if defined
		double kine0 = 0.0;
		double kine1 = 0.0;
		double kine_fa = 0.0;
		if (kine_index > -1)
		{
			kine0 = scalarTranspEqVector[kine_index].phi[icv0];
			kine1 = scalarTranspEqVector[kine_index].phi[icv1];
			kine_fa = w1*kine0 + w0*kine1;
		}

		// calculate viscous flux
#ifdef USE_ARTIF_VISC
		// If the user wants to have artificial bulk viscosity, pass the information to the addViscFlux() method
		// (If the user wants to use artificial viscosity (not bulk only), it was already reflected in mul_fa)
		double artifBulkViscosity = 0.0;
		if(turnOnArtifVisc && artifVisc_bulkViscOnly)
			artifBulkViscosity = w1*artifVisc_mag[icv0] + w0*artifVisc_mag[icv1];

		addViscFlux(Frhou, FrhoE, A0, A1,
				rho[icv0], vel[icv0], grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0], gamma[icv0], kine0,
				rho[icv1], vel[icv1], grad_u[icv1], enthalpy[icv1], grad_enthalpy[icv1], temp[icv1], RoM[icv1], gamma[icv1], kine1,
				mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kine_fa, uAux_fa,
				area, nVec, smag, sVec, artifBulkViscosity);
#else
		addViscFlux(Frhou, FrhoE, A0, A1,
				rho[icv0], vel[icv0], grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0], gamma[icv0], kine0,
				rho[icv1], vel[icv1], grad_u[icv1], enthalpy[icv1], grad_enthalpy[icv1], temp[icv1], RoM[icv1], gamma[icv1], kine1,
				mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kine_fa, uAux_fa,
				area, nVec, smag, sVec);
#endif

		if (flagImplicit)
		{
			for (int i=1; i<5; i++)
				for (int j=0; j<5; j++)
				{
					A[noc00][i][j] -= A0[i][j];
					A[noc01][i][j] -= A1[i][j];
				}

			if (icv1 < ncv)  // if icv1 is internal...
			{
				for (int i=1; i<5; i++)
					for (int j=0; j<5; j++)
					{
						A[noc11][i][j] += A1[i][j];
						A[noc10][i][j] += A0[i][j];
					}
			}
		}

		// icv0 is always valid...
		for (int i = 0; i < 3; i++)
			rhs_rhou[icv0][i] -= Frhou[i];
		rhs_rhoE[icv0] -= FrhoE;

		// icv1 can be ghost...
		if (icv1 < ncv)
		{
			for (int i = 0; i < 3; i++)
				rhs_rhou[icv1][i] += Frhou[i];
			rhs_rhoE[icv1] += FrhoE;
		}
	}


	// ====================================================================
	// cycle through boundary faces, assembling flux and matrix
	// ====================================================================
	for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
	{
		if (zone->getKind() == FA_ZONE_BOUNDARY)
		{
			Param *param;
			if (getParam(param, zone->getName()))
			{
				// .............................................................................................
				// SYMMETRY BOUNDARY CONDITION
				// .............................................................................................
				if ((param->getString() == "SYMMETRY") || (param->getString() == "NOTHING"))
				{
					if (kine_index > -1)
					{
						for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
						{
							int icv0 = cvofa[ifa][0];
							int icv1 = cvofa[ifa][1];
							assert( icv0 >= 0 );
							double nVec[3] = {0.0, 0.0, 0.0};
							double area = normVec3d(nVec, fa_normal[ifa]);

							double *phi = scalarTranspEqVector[kine_index].phi;
							double kine_fa  = phi[icv1];

							double tmp = 1.0/3.0*(rho[icv0] + rho[icv1])*kine_fa;

							for (int i = 0; i < 3; i++)
								rhs_rhou[icv0][i] -= tmp*fa_normal[ifa][i];
						}
					}
				}


				// .............................................................................................
				// WALL BOUNDARY CONDITION
				// .............................................................................................
				else if (param->getString() == "WALL")
				{
					for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ++ifa) {
						//for (int index = 0; index < zone->faVec.size(); ++index) {
						//int ifa = zone->faVec[index];
						int icv0 = cvofa[ifa][0];
						int icv1 = cvofa[ifa][1];

						double nVec[3] = {0.0, 0.0, 0.0};
						double area = normVec3d(nVec, fa_normal[ifa]);
						double sVec[3] = {0.0, 0.0, 0.0};
						vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
						double smag = fabs(vecDotVec3d(sVec, nVec));  // project sVec to wall face normal

						double kine0 = 0.0;
						double kine1 = 0.0;
						double kine_fa = 0.0;
						if (kine_index > -1)
							kine0 = scalarTranspEqVector[kine_index].phi[icv0];

						// calculate viscous flux
						addViscFlux(Frhou, FrhoE, A0, NULL,
								rho[icv0],    vel[icv0],    grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0],    gamma[icv0],  kine0,
								rho[icv1],    vel[icv1],    grad_u[icv0], enthalpy[icv1], grad_enthalpy[icv0], temp[icv1], RoM[icv1],    gamma[icv1],  kine1,
								mul_fa[ifa], 0.0, lamOcp_fa[ifa], kine_fa, vel[icv1],
								area, nVec, smag, nVec);  /* <- use nVec here instead of sVec, to avoid inaccurate correction*/

						if (flagImplicit)
						{
							int noc00 = nbocv_i[icv0]; // icv0's diagonal

							for (int i=1; i<5; i++)
								for (int j=0; j<5; j++)
									A[noc00][i][j] -= A0[i][j];
						}

						for (int i = 0; i < 3; i++)
							rhs_rhou[icv0][i] -= Frhou[i];
						rhs_rhoE[icv0] -= FrhoE;
					}
				}
				// .............................................................................................
				// OTHER BOUNDARY CONDITIONS
				// .............................................................................................
				else
				{
					for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ++ifa) {
						//for (int index = 0; index < zone->faVec.size(); ++index) {
						//int ifa = zone->faVec[index];
						int icv0 = cvofa[ifa][0];
						int icv1 = cvofa[ifa][1];

						double nVec[3] = {0.0, 0.0, 0.0};
						double area = normVec3d(nVec, fa_normal[ifa]);
						double sVec[3] = {0.0, 0.0, 0.0};
						vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
						double smag = normVec3d(sVec);

						double kine0 = 0.0;
						double kine1 = 0.0;
						double kine_fa = 0.0;
						if (kine_index > -1)
						{
							kine0 = scalarTranspEqVector[kine_index].phi[icv0];
							kine1 = scalarTranspEqVector[kine_index].phi[icv1];
							kine_fa = kine1;
						}

						// calculate viscous flux
						addViscFlux(Frhou, FrhoE, A0, NULL,
								rho[icv0], vel[icv0], grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0], gamma[icv0],  kine0,
								rho[icv1], vel[icv1], grad_u[icv0], enthalpy[icv1], grad_enthalpy[icv0], temp[icv1], RoM[icv1], gamma[icv1],  kine1,
								mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kine_fa, vel[icv1],
								area, nVec, smag, sVec);

						if (flagImplicit)
						{
							int noc00 = nbocv_i[icv0]; // icv0's diagonal

							for (int i=1; i<5; i++)
								for (int j=0; j<5; j++)
									A[noc00][i][j] -= A0[i][j];
						}

						for (int i = 0; i < 3; i++)
							rhs_rhou[icv0][i] -= Frhou[i];
						rhs_rhoE[icv0] -= FrhoE;
					}
				}
			}
		}
	}

	if (A0  != NULL)  delete [] A0;
	if (A1  != NULL)  delete [] A1;
}

void JoeWithModels::setBC()
{
  static int first = 1;
  int bc_err = 0;

  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      Param *param;

      if (getParam(param, zone->getName()))
      {
        // .............................................................................................
        // HOOK BOUNDARY CONDITION
        // .............................................................................................
        if (param->getString() == "HOOK")
        {
          if ((first) && (mpi_rank == 0))
            cout << "Applying HOOK                to zone: "<< zone->getName() << endl;

          boundaryHook(temp, vel, press, &(*zone));

          setScalarBC(&(*zone));
          ComputeBCProperties_T(&(*zone));
        }
        // .............................................................................................
        // CBC BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "CBC")
        {
          double u_bc[3], T_bc, p_bc;

          for (int i=0; i<3; i++)
            u_bc[i] = param->getDouble(i+2);
          T_bc = param->getDouble(5);
          p_bc = param->getDouble(6);

          if ((first)&&(mpi_rank == 0))
            cout << "Applying CBC                 to zone: "<< zone->getName() <<"\t u_bc: "<< u_bc[0]<< " "<< u_bc[1]<< " "
                 << u_bc[2] <<" T_bc: " << T_bc << " p_bc: " << p_bc << endl;
	  
          for (int index = 0; index < zone->faVec.size(); ++index) {
	    int ifa = zone->faVec[index];
	    
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];
	    
	    // depending on where the face is, the cvs should be in different ranges...
	    if ( (ifa >= 0)&&(ifa < nfa_b) ) {
	      // old boundary face...
	      assert((icv0 >= 0)&&(icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));
	    }
	    else {
	      // b2 boundary face, between ghost 1 and fake 2...
	      assert( (ifa >= nfa)&&(ifa < nfa_b2) );
	      assert((icv0 >= ncv)&&(icv0 < ncv_g)&&(icv1>=ncv_ggf)&&(icv1<ncv_ggff));
	    }

            double nVec[3];
            double area = normVec3d(nVec, fa_normal[ifa]);

            if (vecDotVec3d(u_bc, nVec) > 0.0)   // outlet
            {
              double velMagN = vecDotVec3d(vel[icv0], nVec);
              double mach = fabs(velMagN)/sos[icv0];

              temp[icv1] = temp[icv0];
              for (int i=0; i<3; i++)
                vel[icv1][i] = vel[icv0][i];
              if (mach >= 1.0)  press[icv1] = press[icv0]; 
              else              press[icv1] = p_bc;
            }
            else      // inlet
            {
              temp[icv1] = T_bc;
              for (int i=0; i<3; i++)
                vel[icv1][i] = u_bc[i];
              press[icv1] = p_bc;
            }
          }

          setScalarBC(&(*zone));
          ComputeBCProperties_T(&(*zone));
       }
        // .............................................................................................
        // CBC SUBSONIC INLET BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "CBC_SUBSONIC_INLET")
        {
          double angleU[3], Ttot, htot, ptot;

          for (int i=0; i<3; i++)
            angleU[i] = param->getDouble(i+2);
          Ttot = param->getDouble(5);
          ptot = param->getDouble(6);

          if ((first)&&(mpi_rank == 0))
            cout << "Applying CBC_SUBSONIC        to zone: "<< zone->getName() <<"\t angleU: "<<angleU[0]<<" "<<angleU[1]<<" "<<angleU[2]
                 << " Ttot: " << Ttot << " Ptot: "<< ptot << endl;

	   for (int index = 0; index < zone->faVec.size(); ++index) {
            int ifa = zone->faVec[index];

            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];
            //assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

            double u[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
            double wPow2 = vecDotVec3d(u, u);               // velocity squared
            double velMag = sqrt(wPow2);                    // apply angle to extrapolated velocity
            for (int i=0; i<3; i++)
              vel[icv1][i] = angleU[i]*velMag;
            temp[icv1] = Ttot;                              // save temporary total temperature (to compute then total enthalpy)
          }

          setScalarBC(&(*zone));
          ComputeBCProperties_T(&(*zone));                  // total enthalpy from total temperature

	   for (int index = 0; index < zone->faVec.size(); ++index) {
            int ifa = zone->faVec[index];
            int icv1 = cvofa[ifa][1];
            double wPow2 = vecDotVec3d(vel[icv1], vel[icv1]);         // velocity squared
            enthalpy[icv1] -= 0.5* wPow2;                                   // static enthalpy
          }

          ComputeBCProperties_H(&(*zone));                  // static temperature and thermo properties from static enthalpy

          // Assumes isentropic relations to determine static pressure (= constant cp)
          // At first approximation ok, but could be improved; should for now be considered in defining p_bc
          for (int index = 0; index < zone->faVec.size(); ++index) {
	    int ifa = zone->faVec[index];
            int icv1 = cvofa[ifa][1];
            press[icv1] = ptot*pow(temp[icv1]/Ttot, gamma[icv1]/(gamma[icv1]-1.0));
          }
        }
        // .............................................................................................
        // CBC SUBSONIC OUTLET BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "CBC_SUBSONIC_OUTLET")
        {
          double p_bc = param->getDouble(2);

          if ((first)&&(mpi_rank == 0))
            cout << "Applying CBC_SUBSONIC_OUTLET to zone: "<< zone->getName() << "\t pOut: "<< p_bc << endl;

          for (int index = 0; index < zone->faVec.size(); ++index) {
	    int ifa = zone->faVec[index];
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];
            //assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

            // Assumes that the temperature is constant across outlet boundary (extrapolation) while p_bc is prescribed
            press[icv1] = p_bc;
            for (int i=0; i<3; i++)
              vel[icv1][i] = vel[icv0][i];
	    
            temp[icv1] = temp[icv0];
          }
	  
          setScalarBC(&(*zone));
          ComputeBCProperties_T(&(*zone));
        }
        // .............................................................................................
        // SYMMETRY BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "SYMMETRY")
        {
        	if ((first)&&(mpi_rank == 0))
        		cout << "Applying SYMMETRY            to zone: "<< zone->getName() << endl;

        	for (int index = 0; index < zone->faVec.size(); ++index) {
        		int ifa = zone->faVec[index];
        		//cout << mpi_rank<< " "<<ifa << endl;
        		int icv0 = cvofa[ifa][0];
        		int icv1 = cvofa[ifa][1];

        		//assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

        		double nVec[3];
        		double area = normVec3d(nVec, fa_normal[ifa]);

        		// flip u, APPROXIMATION ---> take velocity at the cell center
        		double u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};

        		double un = vecDotVec3d(nVec, u0);
        		for (int i = 0; i < 3; i++)
        			vel[icv1][i] = u0[i] - 1.0*un*nVec[i];
        		//assert(fabs(vecDotVec3d(vel_bfa[ifa], nVec)) < 1.0e-10);

        		temp[icv1] = temp[icv0];

        		press[icv1] = press[icv0];
        	}

        	setScalarBC(&(*zone));
        	ComputeBCProperties_T(&(*zone));
        }
        // .............................................................................................
        // NEUMANN BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "NEUMANN")
        {
          if ((first)&&(mpi_rank == 0))
            cout << "Applying NEUMANN             to zone: "<< zone->getName() << endl;

          for (int index = 0; index < zone->faVec.size(); ++index) {
	    int ifa = zone->faVec[index];
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];

            for (int i = 0; i < 3; i++)
              vel[icv1][i] = vel[icv0][i];
	    
            temp[icv1] = temp[icv0];
            press[icv1] = press[icv0];
          }

          setScalarBC(&(*zone));
          ComputeBCProperties_T(&(*zone));
        }
        // .............................................................................................
        // WALL BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "WALL")
        {
          int i=0;
          double T_bc = 0.0;
          if ((i = param->findString("TEMP")) != 0)
            T_bc = param->getDouble(i+1);

          if ((first)&&(mpi_rank == 0))
          {
            if (T_bc > 0.0)    cout << "Applying WALL isothermal     to zone: "<< zone->getName() << "\t Temp_BC: " << T_bc << endl;
            else               cout << "Applying WALL adiabatic      to zone: "<< zone->getName() << endl;
          }

          for (int index = 0; index < zone->faVec.size(); ++index) {
	    int ifa = zone->faVec[index];
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];
            //assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

            for (int i = 0; i < 3; i++)
              vel[icv1][i] = 0.0;

            if (T_bc > 0.0)   temp[icv1] = T_bc;           // wall temperature
            else              temp[icv1] = temp[icv0];     // adiabatic wall

            press[icv1] = press[icv0];                      // Assumes zero pressure gradient at the wall
          }

          setScalarBC(&(*zone));
          ComputeBCProperties_T(&(*zone));
        }
        // .............................................................................................
        // OTHER BOUNDARY CONDITIONS
        // .............................................................................................
        else
	  {
	    if (mpi_rank == 0)
            cerr << "Error: cannot recognize bc for zone: " << zone->getName() << " bc: " << param->getString() << endl;
          bc_err = 1;
        }
      }
      else
      {
        if (mpi_rank == 0)
          cerr << "Error: no bc set for: "<< zone->getName() << endl;
        bc_err = 1;
      }
    }

  // update density at boundary using EOS
  for (int ifa = 0; ifa < nfa_b; ifa++) {
    int icv1 = cvofa[ifa][1];
    rho[icv1] = press[icv1]/(RoM[icv1]*temp[icv1]);
    sos[icv1] = sqrt(gamma[icv1] * press[icv1] / rho[icv1]);
  }
  for (int ifa = nfa; ifa < nfa_b2; ifa++) {
    int icv1 = cvofa[ifa][1];
    rho[icv1] = press[icv1]/(RoM[icv1]*temp[icv1]);
    sos[icv1] = sqrt(gamma[icv1] * press[icv1] / rho[icv1]);
  }
  
  if (bc_err != 0)
    throw(-1);

  first = 0;

}

void JoeWithModels::showResidue(double *rhsResid)
{
  int nScal = scalarTranspEqVector.size();

  // residual label at every 10 output steps
  if ((mpi_rank == 0) && ((step%(check_interval*10) == 0) || (step == 1)))
  {
    printf("                rho          rhou-X       rhou-Y       rhou-Z       rhoE      ");
    for (int iScal = 0; iScal < nScal; iScal++)
      printf("%12s", scalarTranspEqVector[iScal].getName());
    cout << endl;
  }

  // residual value at each output step
  if (mpi_rank == 0)
  {
    printf("RESID: %6d %12.4e %12.4e %12.4e %12.4e %12.4e", step, rhsResid[0], rhsResid[1], rhsResid[2], rhsResid[3], rhsResid[4]);
    for (int iScal = 0; iScal < nScal; iScal++)
      printf("%12.4e", rhsResid[5+iScal]);
    cout << endl;
  }
}

// *************** COUPLED SOLUTIONS (UNDER DEBUGGING) ******************************

void JoeWithModels::runBackwardEulerCoupled()
{
  int nScal = scalarTranspEqVector.size();
  for (int iScal = 0; iScal < nScal; iScal++)
    scalarTranspEqVector[iScal].coupling = "COUPLED";

  double *myResidual = new double[5+nScal];
  double *Residual   = new double[5+nScal];
  
  double ***A;           getMem3D(&A,   0, nbocv_s-1, 0, 5+nScal-1, 0, 5+nScal-1, "JoeWithModels::runBackwardEulerCoupled -> A",   true);
  double **rhs;          getMem2D(&rhs, 0, ncv-1,     0, 5+nScal-1,               "JoeWithModels::runBackwardEulerCoupled -> rhs", true);
  double **dq;           getMem2D(&dq,  0, ncv_g-1,   0, 5+nScal-1,               "JoeWithModels::runBackwardEulerCoupled -> dq",  true);

  // ERROR ESTIMATE STUFF
  double *RHSrho ;
  double (*RHSrhou)[3] ;
  double *RHSrhoE ;        
  RHSrho   = NULL;         registerScalar(RHSrho,  "RHSRHO",  CV_DATA);
  RHSrhou  = NULL;         registerVector(RHSrhou, "RHSRHOU", CV_DATA);
  RHSrhoE  = NULL;         registerScalar(RHSrhoE, "RHSRHOE", CV_DATA);
  // ERROR ESTIMATE STUFF

  
  //------------------------------------
  // some parameters
  //------------------------------------
  double underRelax = getDoubleParam("UNDER_RELAXATION", "0.3");

  if (!checkParam("LINEAR_SOLVER_NS_TRESHOLDS"))
  {
    ParamMap::add("LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2");    // add default values
    if (mpi_rank == 0)
      cout << "WARNING: added keyword \"LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2\"" <<
              " to parameter map!" << endl;
  }
  int maxIterLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getInt("MAX_ITER");
  double zeroAbsLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("ABS_RESID");
  double zeroRelLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("REL_RESID");

  if (!checkParam("CFL_RAMP"))
  {
    ParamMap::add("CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0");    // add default: no increase of CFL number!
    if (mpi_rank == 0)
      cout << "WARNING: added \"CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0\" to parameter map" << endl;
  }

  int startIncCFL = getParam("CFL_RAMP")->getInt("AFTER_ITER");
  int intervalIncCFL = getParam("CFL_RAMP")->getInt("INTERVAL_ITER");
  double incCFL = getParam("CFL_RAMP")->getDouble("FACTOR_CFL");
  double maxCFL = getParam("CFL_RAMP")->getDouble("MAX_CFL");


  // -------------------------------------------------------------------------------------------
  // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
  // -------------------------------------------------------------------------------------------
  calcStateVariables();

  // -------------------------------------------------------------------------------------------
  // update material properties: laminar viscosity and heat conductivity
  // -------------------------------------------------------------------------------------------
  calcMaterialProperties();

  // -------------------------------------------------------------------------------------------
  // set BC's for NS and scalars
  // -------------------------------------------------------------------------------------------
  setBC();

  // -------------------------------------------------------------------------------------------
  // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
  // -------------------------------------------------------------------------------------------
  calcRansTurbViscMuet();


  // -------------------------------------------------------------------------------------------
  //
  //   Loop over time steps
  //
  // -------------------------------------------------------------------------------------------
  int done = doneSolver(1.0e20);

  if (initial_flowfield_output == "YES")     writeData(0);

  // provide total runtime
  double wtime, wtime0;
  MPI_Barrier(mpi_comm);
  if (mpi_rank == 0)
    wtime = MPI_Wtime();

  while (done != 1)
  {
    step++;
    if ((step >= startIncCFL) && (step%intervalIncCFL == 0) && (cfl < maxCFL))      cfl *= incCFL;
    double dt_min = calcDt(cfl);

    // ---------------------------------------------------------------------------------
    // Compute RHS for both NSE and scalars
    // ---------------------------------------------------------------------------------
    for (int noc = 0; noc < nbocv_s; noc++)                             // set A, dq to zero! rhs is set zero in calcRHS
      for (int i = 0; i < 5+nScal; i++)
        for (int j = 0; j < 5+nScal; j++)
          A[noc][i][j] = 0.0;
    
    for (int icv = 0; icv < ncv_g; icv++)
      for (int i = 0; i < 5+nScal; i++)
        dq[icv][i] = 0.0;

    // ---------------------------------------------------------------------------------
    // calculate RHS for both NSE and scalars
    // ---------------------------------------------------------------------------------
    
    calcRhsCoupled(rhs, A, nScal, true);
	
     //ERROR ESTIMATE STUFF
    for (int icv = 0; icv < ncv; icv++)                                     // update solution for NSE: Q_new = Q_old + Delta_Q
    { 
      RHSrho[icv]        = rhs[icv][0];
      RHSrhou[icv][0]    = rhs[icv][1];
      RHSrhou[icv][1]    = rhs[icv][2];
      RHSrhou[icv][2]    = rhs[icv][3];
      RHSrhoE[icv]       = rhs[icv][4];
      residField[icv]    = rhs[icv][4];
    }

    for (int iScal = 0; iScal < nScal; iScal++)                          
     {
      string scalname = scalarTranspEqVector[iScal].getName();
      for (int icv = 0; icv < ncv; ++icv)
      {
        if (iScal == getScalarTransportIndex("omega") || iScal == getScalarTransportIndex("nuSA")) residField2[icv] = rhs[icv][5+iScal];
        if (iScal == getScalarTransportIndex("kine"))  residField3[icv] = rhs[icv][5+iScal];
        if (iScal == getScalarTransportIndex("ZMean")) residField4[icv] = rhs[icv][5+iScal];
        if (iScal == getScalarTransportIndex("ZVar"))  residField5[icv] = rhs[icv][5+iScal];
        if (iScal == getScalarTransportIndex("CMean")) residField6[icv] = rhs[icv][5+iScal];
      }
     }
     //ERROR ESTIMATE STUFF
    
    // ---------------------------------------------------------------------------------
    // calculate residual
    // ---------------------------------------------------------------------------------

    for (int i = 0; i < 5+nScal; i++)
    {
      myResidual[i] = 0.0;
      Residual[i] = 0.0;
    }
    
    for (int icv = 0; icv < ncv; icv++)
      for (int i = 0; i < 5+nScal; i++)
        myResidual[i] += fabs(rhs[icv][i]);
      
    MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);


    // ---------------------------------------------------------------------------------
    // solve linear system for the NSE
    // ---------------------------------------------------------------------------------
    
    for (int icv = 0; icv < ncv; icv++)                                     // prepare rhs and A
    {
      for (int i = 0; i < 5+nScal; i++)                                     // relaxation
        rhs[icv][i] *= underRelax;

      double tmp = cv_volume[icv]/(local_dt[icv]);
      for (int i = 0; i < 5+nScal; i++)
        A[nbocv_i[icv]][i][i] += tmp;                                   // diagonal part ( vol/dt + A )
    }      
    
    solveCoupledLinSysNSCoupled(dq, A, rhs, zeroAbsLS, zeroRelLS, maxIterLS, nScal);  // solve linear system
             
    for (int icv = 0; icv < ncv; icv++)                                     // update solution for NSE: Q_new = Q_old + Delta_Q
    {
      rho[icv]     += dq[icv][0];
      rhou[icv][0] += dq[icv][1];
      rhou[icv][1] += dq[icv][2];
      rhou[icv][2] += dq[icv][3];
      rhoE[icv]    += dq[icv][4];
    }

    updateCvDataG1G2(rho,  REPLACE_DATA);
    updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
    updateCvDataG1G2(rhoE, REPLACE_DATA);

    for (int iScal = 0; iScal < nScal; iScal++)                         // update + clip solution for scalars: Q_new = (rho_old * Q_old + Delta_Q ) / rho_new
    {
      double *phi = scalarTranspEqVector[iScal].phi;
      for (int icv = 0; icv < ncv; icv++)
        phi[icv] = min(max((phi[icv] * (rho[icv] - dq[icv][0]) + dq[icv][5+iScal]) / rho[icv], scalarTranspEqVector[iScal].lowerBound), scalarTranspEqVector[iScal].upperBound);
      updateCvDataG1G2(phi, REPLACE_DATA);
    }


    // -------------------------------------------------------------------------------------------
    // update state properties: velocity, pressure, temperature, enthalpy, gamma and R
    // -------------------------------------------------------------------------------------------
    calcStateVariables();

    // -------------------------------------------------------------------------------------------
    // update material properties: laminar viscosity and heat conductivity
    // -------------------------------------------------------------------------------------------
    calcMaterialProperties();

    // -------------------------------------------------------------------------------------------
    // set BC's for NS and scalars
    // -------------------------------------------------------------------------------------------
    setBC();

    // -------------------------------------------------------------------------------------------
    // update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
    // -------------------------------------------------------------------------------------------
    calcRansTurbViscMuet();


    // =========================================================================================
    // show residual
    // =========================================================================================
    if (step%check_interval == 0)
    {
      if ((mpi_rank == 0) && (step%(check_interval*10) == 0))
        cout << "\ndone step: "<< step << ", cfl: " << cfl << ", min. dt: " << dt_min << " time: " << time << endl;

      showResidue(Residual);
    }

    // IKJ: For output display (Metric to see the convergence of the simulation)
    for (int icv = 0; icv < ncv; icv++) {
    	log10_resid_rhoE[icv] = log10(fabs(RHSrhoE[icv]) + 1.0e-15);
    }
    for (int iScal = 0; iScal < nScal; iScal++)
    	for (int icv = 0; icv < ncv; icv++)
    		log10_resid_scalar0[icv] = log10(fabs(rhs[icv][5+iScal]/underRelax) + 1.0e-15);
    updateCvDataG1G2(residField,       REPLACE_DATA);
    updateCvDataG1G2(log10_resid_rhoE, REPLACE_DATA);
    if(nScal>0)
    	updateCvDataG1G2(log10_resid_scalar0, REPLACE_DATA);

    // IKJ
    double myTotResid_dq  = 0.0;
    double myTotResid_rhs = 0.0;
    for (int icv = 0; icv < ncv; icv++) {
    	for(int i=0; i<5; ++i)
    		myTotResid_dq += dq[icv][i];

    	myTotResid_rhs += fabs(RHSrho[icv]);
    	for (int i=0; i<3; i++)
    		myTotResid_rhs += fabs(RHSrhou[icv][i]);
    	myTotResid_rhs += fabs(RHSrhoE[icv]);
    }
    for (int iScal = 0; iScal < nScal; iScal++) {
    	double *phi = scalarTranspEqVector[iScal].phi;
    	for (int icv = 0; icv < ncv; icv++) {
    		myTotResid_dq  += fabs( (dq[icv][5+iScal] - phi[icv]*dq[icv][0]) / (rho[icv] - dq[icv][0]) );
    		myTotResid_rhs += fabs(rhs[icv][5+iScal]);
    	}
    }
    MPI_Allreduce(&myTotResid_dq,  &totResid_dq,  1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&myTotResid_rhs, &totResid_rhs, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);


    temporalHook();
    dumpProbes(step, 0.0);
    writeData(step);

    if ((write_restart > 0) && (step % write_restart == 0))
      writeRestart(step);

    done = doneSolver(Residual[4]);   // pass the energy residual to determine job cancellation
  }

  MPI_Barrier(mpi_comm);
  if (mpi_rank == 0)
  {
    double wtime0 = wtime;
    wtime = MPI_Wtime();
    cout << " > runtime for iterations[s]: " << wtime - wtime0 << endl;
  }


  // ---------------------------------------------------------------------------------
  // output
  // ---------------------------------------------------------------------------------

  temporalHook();
  finalHook();

  writeRestart();



  // ---------------------------------------------------------------------------------
  // delete memory
  // ---------------------------------------------------------------------------------
  delete [] myResidual;
  delete [] Residual;
  
  freeMem3D(A,   0, nbocv_s-1, 0, 5+nScal-1, 0, 5+nScal-1);         A   = NULL;
  freeMem2D(rhs, 0, ncv-1,     0, 5+nScal-1);                       rhs = NULL;
  freeMem2D(dq,  0, ncv_g-1,   0, 5+nScal-1);                       dq  = NULL;

}


int JoeWithModels::calcRhsCoupled(double **rhs, double ***A, int nScal, int flagImplicit)
{
  // set RHS to zero
  for (int icv = 0; icv < ncv; icv++)
    for (int i = 0; i < 5+nScal; i++)
      rhs[icv][i] = 0.0;
  
  // compute Euler Flux for NS and scalars
  int CountReducedOrder = calcFluxCoupled(rhs, A, nScal, flagImplicit);
  
  // add source terms to RHS of Navier-Stokes equations
  sourceHookCoupled(rhs, A, nScal, flagImplicit);
  sourceHookRansTurbCoupled(rhs, A, nScal, flagImplicit);
  sourceHookRansCombCoupled(rhs, A, nScal, flagImplicit);

  return CountReducedOrder;
}

int JoeWithModels::calcFluxCoupled(double **rhs, double ***A, int nScal, int flagImplicit)
{
  double **Apl = NULL;
  double **Ami = NULL;
  double **A0  = NULL;
  double **A1  = NULL;
  
  if (flagImplicit)
  {
    getMem2D(&Apl, 0, 5+nScal-1, 0, 5+nScal-1, "JoeWithModels::calcFluxCoupled -> Apl", true);
    getMem2D(&Ami, 0, 5+nScal-1, 0, 5+nScal-1, "JoeWithModels::calcFluxCoupled -> Ami", true);
    for (int i = 0; i < 5+nScal; i++)
      for (int j = 0; j < 5+nScal; j++)
      {
        Apl[i][j] = 0.0;
        Ami[i][j] = 0.0;
      }
    if (mu_ref > 0.0)
    {
      getMem2D(&A0,  0, 5+nScal-1, 0, 5+nScal-1, "JoeWithModels::calcFluxCoupled -> A0 ", true);
      getMem2D(&A1,  0, 5+nScal-1, 0, 5+nScal-1, "JoeWithModels::calcFluxCoupled -> A1 ", true);
      for (int i = 0; i < 5+nScal; i++)
        for (int j = 0; j < 5+nScal; j++)
        {
          A0[i][j] = 0.0;
          A1[i][j] = 0.0;
        }
    }
  }

  double *EulerFlux   = new double[5+nScal];
  double *ViscousFlux = new double[5+nScal];
  
  for (int i = 0; i < 5+nScal; i++)
  {
    EulerFlux[i]   = 0.0;
    ViscousFlux[i] = 0.0;
  }
  
  double *Scalar0        = NULL;         if (nScal > 0) Scalar0       = new double[nScal];            // cell face if 2nd order, cell center if 1st order
  double *Scalar1        = NULL;         if (nScal > 0) Scalar1       = new double[nScal];            // cell face if 2nd order, cell center if 1st order
  double *ScalCV0        = NULL;         if (nScal > 0) ScalCV0       = new double[nScal];            // cell center
  double *ScalCV1        = NULL;         if (nScal > 0) ScalCV1       = new double[nScal];            // cell center
  double (*gradScal0)[3] = NULL;         if (nScal > 0) gradScal0     = new double[nScal][3];         // gradient of scalars
  double (*gradScal1)[3] = NULL;         if (nScal > 0) gradScal1     = new double[nScal][3];         // gradient of scalars
  double *dpress_dscal0  = NULL;         if (nScal > 0) dpress_dscal0 = new double[nScal];            // derivative of pressure with respect to scalars
  double *dpress_dscal1  = NULL;         if (nScal > 0) dpress_dscal1 = new double[nScal];            // derivative of pressure with respect to scalars
  double *diffScal       = NULL;         if (nScal > 0) diffScal      = new double[nScal];            // diffusivity of scalars
  double *ConvTerm       = NULL;         if (nScal > 0) ConvTerm      = new double[nScal];            // 0 if convective term not considered, otherwise 1
  double *DiffTerm       = NULL;         if (nScal > 0) DiffTerm      = new double[nScal];            // 0 if diffusive  term not considered, otherwise 1
  
  // count how many cells switched back to first order due to extrapolated negative density or pressure/temperature at the faces
  int CountReducedOrder = 0;
  int myCountReducedOrder = 0;

  // save the index of kine if defined and save convTerm for speed
  int kine_Index = getScalarTransportIndex("kine");
    
  for (int iScal = 0; iScal < nScal; iScal++)
  {
    ConvTerm[iScal] = (double)scalarTranspEqVector[iScal].convTerm;
    DiffTerm[iScal] = (double)scalarTranspEqVector[iScal].diffTerm;
    dpress_dscal0[iScal] = 0.0;
    dpress_dscal1[iScal] = 0.0;
  }
  
  // =============================================================================================
  // compute gradients, with boundary values
  // =============================================================================================
  if ((sndOrder == true) || (mu_ref > 0.0))
    calcCv2Grad(grad_u, vel, limiterNavierS, sos, epsilonSDWLS);                      // velocity gradients
  
  if (sndOrder == true)
  {
   if (stOrderScalar == false) 
   {
    calcCv2Grad(grad_rho, rho, limiterNavierS, rho, epsilonSDWLS);                    // density gradients
#ifdef temp_reconstruction
    calcCv2Grad(grad_temp, temp, limiterNavierS, temp, epsilonSDWLS);                   // temperature gradients
#else
    calcCv2Grad(grad_p, press, limiterNavierS, press, epsilonSDWLS);                    // pressure gradients
#endif
   }else{
	for(int icv=0; icv<ncv_g; icv++){
	   for(int i=0; i<3;i++){
		grad_rho[icv][i]=0.0;
#ifdef temp_reconstruction
		grad_temp[icv][i]=0.0;
#else
		grad_p[icv][i]=0.0;
#endif
	   }
	}
  }
  }
  
  if (mu_ref > 0.0)
    calcCv2Grad(grad_enthalpy, enthalpy, limiterNavierS, enthalpy, epsilonSDWLS);       // enthalpy gradients
  
  for (int iScal = 0; iScal < nScal; iScal++)
  {
    string scalName(scalarTranspEqVector[iScal].getName());

    if ((scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE") && (sndOrder == true)  && (stOrderScalar == false))
    {
      // For the scalar reconstruction, Grad(rho * phi) is required
      double *phi              = scalarTranspEqVector[iScal].phi;
      double *rhoPhi           = scalarTranspEqVector[iScal].rhophi;
      double (*grad_rhophi)[3] = scalarTranspEqVector[iScal].grad_rhophi;

      // Compute rho * phi
      for (int icv = 0; icv < ncv_ggff; icv++)
        rhoPhi[icv] = rho[icv] * phi[icv];

      // Compute gradients of rho*Phi and limit based on rho*Phi
      calcCv2Grad(grad_rhophi, rhoPhi, limiterNavierS, rhoPhi, epsilonSDWLS);     // scalars gradients
    }

    if ( ((scalarTranspEqVector[iScal].reconstruction == "STANDARD") && (sndOrder == true) && (stOrderScalar == false)) || ((mu_ref > 0.0) /*&& (scalarTranspEqVector[iScal].diffTerm == 1.0) */) )
    {
      double *phi           = scalarTranspEqVector[iScal].phi;
      double (*grad_phi)[3] = scalarTranspEqVector[iScal].grad_phi;
      calcCv2Grad(grad_phi, phi, limiterNavierS, phi, epsilonSDWLS);                 // scalars gradients
    }

    if ((mu_ref > 0.0) /*&& (scalarTranspEqVector[iScal].diffTerm == 1.0)*/)
    {
      diffusivityHookScalarRansTurb(scalName);
      diffusivityHookScalarRansComb(scalName);
    }
// IKJ
    else if (d_scalar_comb_ref > 0.0)
    	diffusivityHookScalarRansComb(scalName);
  }
  
  if (flagImplicit)
  {
    pressureDerivativeHookScalarRansTurb();      // compute pressure derivative dP/dScal
    pressureDerivativeHookScalarRansComb();
  }

  // ===============================================================================================
  // cycle through internal faces, assembling flux to both sides
  // ===============================================================================================
  for (int ifa = nfa_b; ifa < nfa; ifa++)
  {
    int icv0 = cvofa[ifa][0];
    int icv1 = cvofa[ifa][1];
    assert( icv0 >= 0 );
    assert( icv1 >= 0 );

    int noc00, noc01, noc11, noc10;
    if (flagImplicit)
      getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

    // face unit normal and area...
    double nVec[3] = {0.0, 0.0, 0.0};
    double area = normVec3d(nVec, fa_normal[ifa]);

    // .............................................................................................
    // reconstruction of variables at faces: rho, u, T or P, scalars
    // .............................................................................................
    double rho0 = rho[icv0];
    double u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
    double p0 = press[icv0];
    double T0 = temp[icv0];
    double h0 = enthalpy[icv0];
    double gam0 = gamma[icv0];
    double R0 = RoM[icv0];
    double kineCV0 = 0.0;          // cell center
    double kineFA0 = 0.0;          // cell face if second order

    double rho1 = rho[icv1];
    double u1[3] = {vel[icv1][0], vel[icv1][1], vel[icv1][2]};
    double p1 = press[icv1];
    double T1 = temp[icv1];
    double h1 = enthalpy[icv1];
    double gam1 = gamma[icv1];
    double R1 = RoM[icv1];
    double kineCV1 = 0.0;
    double kineFA1 = 0.0;
    
    double kine_fa = 0.0;          // interpolated kinetic energy at face

    for (int iScal = 0; iScal < nScal; iScal++)
    {
      double *phi = scalarTranspEqVector[iScal].phi;
      ScalCV0[iScal] = Scalar0[iScal] = phi[icv0];
      ScalCV1[iScal] = Scalar1[iScal] = phi[icv1];
    }

    if (sndOrder == true)
    {
      double r0[3] = {0.0, 0.0, 0.0};
      double r1[3] = {0.0, 0.0, 0.0};
      vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(r1, x_fa[ifa], x_cv[icv1]);

      // ----------------------------------------
      // left side
      // ----------------------------------------
      rho0 += vecDotVec3d(r0, grad_rho[icv0]);
#ifdef temp_reconstruction
      T0 += vecDotVec3d(r0, grad_temp[icv0]);
      if ((T0 <= 0.0) || (rho0 <= 0.0))
      {
        T0 = temp[icv0];
        rho0 = rho[icv0];
        myCountReducedOrder++;
      }
#else
      p0 += vecDotVec3d(r0, grad_p[icv0]);
      if ((p0 <= 0.0) || (rho0 <= 0.0))
      {
        p0 = press[icv0];
        rho0 = rho[icv0];
        myCountReducedOrder++;
      }
#endif
      else
      {
        for (int i = 0; i < 3; i++)
          u0[i] += vecDotVec3d(r0, grad_u[icv0][i]); 
       if (stOrderScalar == false) 
        for (int iScal = 0; iScal < nScal; iScal++)
        {
          if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
            Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d(r0, scalarTranspEqVector[iScal].grad_rhophi[icv0])) / rho0;
          else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
            Scalar0[iScal] += vecDotVec3d(r0, scalarTranspEqVector[iScal].grad_phi[icv0]);
        }
      }

      // ----------------------------------------
      // right side
      // ----------------------------------------
      rho1 += vecDotVec3d(r1, grad_rho[icv1]);
#ifdef temp_reconstruction
      T1 += vecDotVec3d(r1, grad_temp[icv1]);
      if ((T1 <= 0.0) || (rho1 <= 0.0))
      {
        T1 = temp[icv1];
        rho1 = rho[icv1];
        myCountReducedOrder++;
      }
#else
      p1 += vecDotVec3d(r1, grad_p[icv1]);
      if ((p1 <= 0.0) || (rho1 <= 0.0))
      {
        p1 = press[icv1];
        rho1 = rho[icv1];
        myCountReducedOrder++;
      }
#endif
      else
      {
        for (int i = 0; i < 3; i++)
          u1[i] += vecDotVec3d(r1, grad_u[icv1][i]);
       if (stOrderScalar == false) 
        for (int iScal = 0; iScal < nScal; iScal++)
        {
          if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
            Scalar1[iScal] = (rho[icv1] * Scalar1[iScal] + vecDotVec3d(r1, scalarTranspEqVector[iScal].grad_rhophi[icv1])) / rho1;
          else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
            Scalar1[iScal] += vecDotVec3d(r1, scalarTranspEqVector[iScal].grad_phi[icv1]);
        }
      }
     
      
      // .............................................................................................
      // calculation of other variables at faces: p/T, h, R, gam
      // .............................................................................................
#ifdef temp_reconstruction
      calcThermoProp_T(p0, h0, R0, gam0, rho0, T0, Scalar0, nScal);
      calcThermoProp_T(p1, h1, R1, gam1, rho1, T1, Scalar1, nScal);
#else
      calcThermoProp_p(T0, h0, R0, gam0, rho0, p0, Scalar0, nScal);
      calcThermoProp_p(T1, h1, R1, gam1, rho1, p1, Scalar1, nScal);
#endif

    }


    if (kine_Index > -1)   // save kine if defined
    {
      kineCV0 = ScalCV0[kine_Index];                     // cell center left for implicit side (Jacobi is computed with cell center)
      kineFA0 = Scalar0[kine_Index];                     // cell face left, if second order
      kineCV1 = ScalCV1[kine_Index];                     // cell center left for implicit side (Jacobi is computed with cell center)
      kineFA1 = Scalar1[kine_Index];                     // cell face right,
      
      double dx0[3] = {0.0, 0.0, 0.0};
      double dx1[3] = {0.0, 0.0, 0.0};
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      double w0 = sqrt(vecDotVec3d(dx0, dx0));
      double w1 = sqrt(vecDotVec3d(dx1, dx1));
      kine_fa  = (w1 * kineCV0 + w0 * kineCV1) / (w0 + w1);  // cell face interpolated
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
    // .............................................................................................
    // calculate Euler Flux explicit using HLLC
    // .............................................................................................
    calcEulerFluxCoupled_HLLC(EulerFlux,
             rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
             rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
             area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

    // icv0 is always valid...
    for (int i = 0; i < 5+nScal; i++)
      rhs[icv0][i] -= EulerFlux[i];
    
    // icv1 can be ghost...
    if (icv1 < ncv)
      for (int i = 0; i < 5+nScal; i++)
        rhs[icv1][i] += EulerFlux[i];

    // .............................................................................................
    // calculate Euler implicit matrix using HLLC
    // .............................................................................................
    if (flagImplicit)
    {
      for (int iScal = 0; iScal < nScal; iScal++)
        if (scalarTranspEqVector[iScal].dpress_dphi != NULL)  
        {
          dpress_dscal0[iScal] = scalarTranspEqVector[iScal].dpress_dphi[icv0];
          dpress_dscal1[iScal] = scalarTranspEqVector[iScal].dpress_dphi[icv1];
        }

      calcEulerFluxMatricesCoupled_HLLC(Apl, Ami,
               rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, dpress_dscal0, kineCV0,
               rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, dpress_dscal1, kineCV1,
               area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
//      calcEulerFluxMatricesCoupled_HLLC(Apl, Ami,
//               rho0, u0, p0, T0, h0, R0, gam0, Scalar0, dpress_dscal0, kineFA0,
//               rho1, u1, p1, T1, h1, R1, gam1, Scalar1, dpress_dscal1, kineFA1,
//               area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

      for (int i = 0; i < 5+nScal; i++)
      for (int j = 0; j < 5+nScal; j++)
      {
        A[noc00][i][j] += Apl[i][j];
        A[noc01][i][j] += Ami[i][j];
      }

      if (icv1 < ncv)  // if icv1 is internal...
        for (int i = 0; i < 5+nScal; i++)
        for (int j = 0; j < 5+nScal; j++)
        {
          A[noc11][i][j] -= Ami[i][j];
          A[noc10][i][j] -= Apl[i][j];
        }
    }

    // .............................................................................................
    // calculate viscous Flux explicit and implicit matrix
    // .............................................................................................
    if (mu_ref > 0.0)
    {
      double sVec[3] = {0.0, 0.0, 0.0};
      vecMinVec3d(sVec, x_cv[icv1], x_cv[icv0]);
      double smag = normVec3d(sVec);
      double alpha = vecDotVec3d(nVec, sVec);
      assert((alpha > 0.0) && (alpha < 1.0+1.0E-12));             // alpha should now contain s dot n...
      double uAux_fa[3] = { 0.5 * (vel[icv0][0] + vel[icv1][0]),
                            0.5 * (vel[icv0][1] + vel[icv1][1]),
                            0.5 * (vel[icv0][2] + vel[icv1][2])};
      
      for (int iScal = 0; iScal < nScal; iScal++)
      {
        for (int i = 0; i < 3; i++)
        {
          gradScal0[iScal][i] = scalarTranspEqVector[iScal].grad_phi[icv0][i];
          gradScal1[iScal][i] = scalarTranspEqVector[iScal].grad_phi[icv1][i];
        }
        diffScal[iScal] = scalarTranspEqVector[iScal].diff[ifa];
      }      
      
      calcViscousFluxCoupled(ViscousFlux, A0, A1,
                 rho[icv0], vel[icv0], grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0], gamma[icv0], ScalCV0, gradScal0, dpress_dscal0, kineCV0,
                 rho[icv1], vel[icv1], grad_u[icv1], enthalpy[icv1], grad_enthalpy[icv1], temp[icv1], RoM[icv1], gamma[icv1], ScalCV1, gradScal1, dpress_dscal1, kineCV1,
                 mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kine_fa, uAux_fa, diffScal, DiffTerm, 
                 area, nVec, smag, sVec, alpha, nScal);

      // icv0 is always valid...
      for (int i = 0; i < 5+nScal; i++)
        rhs[icv0][i] -= ViscousFlux[i];
      
      // icv1 can be ghost...
      if (icv1 < ncv)
        for (int i = 0; i < 5+nScal; i++)
          rhs[icv1][i] += ViscousFlux[i];
      
      if (flagImplicit)
      {
        for (int i = 0; i < 5+nScal; i++) 
        for (int j = 0; j < 5+nScal; j++)
        {
          A[noc00][i][j] += A0[i][j];
          A[noc01][i][j] += A1[i][j];
        }
  
        if (icv1 < ncv)  // if icv1 is internal...
          for (int i = 0; i < 5+nScal; i++) 
          for (int j = 0; j < 5+nScal; j++)
          {
            A[noc11][i][j] -= A1[i][j];
            A[noc10][i][j] -= A0[i][j];
          }
      }
    }

  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ===============================================================================================
  // cycle through boundary faces, assembling flux
  // ===============================================================================================
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
  {
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      Param *param;

      if (getParam(param, zone->getName()))
      {
        // .............................................................................................
        // SYMMETRY BOUNDARY CONDITION, ATTENTION: ONLY FIRST ORDER!!!
        // .............................................................................................
        if (param->getString() == "SYMMETRY")
        {
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];
            assert( icv0 >= 0 );
            int noc00 = nbocv_i[icv0];          // icv0's diagonal

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);

            for (int iScal = 0; iScal < nScal; iScal++)
              Scalar0[iScal] = scalarTranspEqVector[iScal].phi[icv0];

            double kineFA = 0.0;
            if (kine_Index > -1)
              kineFA = scalarTranspEqVector[kine_Index].phi[icv0];

            // Euler flux
            calcEulerFluxCoupled_HLLC(EulerFlux,
                rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
                rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
//                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
                     area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");          // more "correct" but apparently less stable...

            for (int i = 0; i < 5+nScal; i++)
              rhs[icv0][i] -= EulerFlux[i];

            if (flagImplicit)
            {
              for (int iScal = 0; iScal < nScal; iScal++)
                if (scalarTranspEqVector[iScal].dpress_dphi != NULL)  
                  dpress_dscal0[iScal] = scalarTranspEqVector[iScal].dpress_dphi[icv0];

              calcEulerFluxMatricesCoupled_HLLC(Apl, NULL,
                       rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, dpress_dscal0, kineFA,
                       rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, dpress_dscal0, kineFA,
//                       area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
                       area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");          // more "correct" but apparently less stable...

              for (int i = 0; i < 5+nScal; i++)
              for (int j = 0; j < 5+nScal; j++)
                A[noc00][i][j] += Apl[i][j];
            }

            // Viscous flux: only 2/3 viscosity times trace of strain rate tensor and 2/3 * rho * kine
            // Viscous flux
            if (mu_ref > 0.0)
            {
              double tmp = 2.0 / 3.0 * ((mul_fa[ifa] + mut_fa[ifa]) * (grad_u[icv0][0][0] + grad_u[icv0][1][1] + grad_u[icv0][2][2]) + rho[icv1] * kineFA);
              rhs[icv0][1] -= area * tmp * nVec[0];
              rhs[icv0][2] -= area * tmp * nVec[1];
              rhs[icv0][3] -= area * tmp * nVec[2];
              
              if (flagImplicit)
              {
                // No implicit term considered here!
              }
            }
            
          }
        }
        
        // .............................................................................................
        // NEUMANN BOUNDARY CONDITION, ATTENTION: ONLY FIRST ORDER!!!
        // ASSUMES THAT IF NSE HAS NEUMANN BC, ALL SCALARS HAVE NEUMANN BC!!!
        // .............................................................................................
        else if (param->getString() == "NEUMANN")
        {
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];
            assert( icv0 >= 0 );
            int noc00 = nbocv_i[icv0];          // icv0's diagonal

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);

            for (int iScal = 0; iScal < nScal; iScal++)
              Scalar0[iScal] = scalarTranspEqVector[iScal].phi[icv1];

            double kineFA = 0.0;
            if (kine_Index > -1)
              kineFA = scalarTranspEqVector[kine_Index].phi[icv1];

            // Euler flux
            calcEulerFluxCoupled_HLLC(EulerFlux,
                     rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
                     rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

            for (int i = 0; i < 5+nScal; i++)
              rhs[icv0][i] -= EulerFlux[i];

            if (flagImplicit)
            {
              for (int iScal = 0; iScal < nScal; iScal++)
                if (scalarTranspEqVector[iScal].dpress_dphi != NULL)  
                  dpress_dscal0[iScal] = scalarTranspEqVector[iScal].dpress_dphi[icv0];

              calcEulerFluxMatricesCoupled_HLLC(Apl, Ami,
                       rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, dpress_dscal0, kineFA,
                       rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, dpress_dscal0, kineFA,
                       area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

              for (int i = 0; i < 5+nScal; i++)
              for (int j = 0; j < 5+nScal; j++)
                A[noc00][i][j] += Apl[i][j]; // + Ami[i][j];   // !!! Apl + Ami for Neumann boundary condition: Qr = Ql -> dF/dQ = dF/dQr + dF/dQl but more unstable
            }

            // Viscous flux
            if (mu_ref > 0.0)
            {
              double sVec[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
              double smag = fabs(vecDotVec3d(sVec, nVec));                // project sVec to wall face normal
              double alpha = 1.0;
              
              for (int iScal = 0; iScal < nScal; iScal++)
              {
                diffScal[iScal] = scalarTranspEqVector[iScal].diff[ifa];
                for (int i = 0; i < 3; i++)
                  gradScal0[iScal][i] = 0.0;
              }
              
              calcViscousFluxCoupled(ViscousFlux, A0, A1,
                         rho[icv1],    vel[icv1],    grad_u[icv0], enthalpy[icv1], grad_enthalpy[icv0], temp[icv1], RoM[icv1],    gamma[icv1], Scalar0, gradScal0, dpress_dscal0, kineFA,
                         rho[icv1],    vel[icv1],    grad_u[icv0], enthalpy[icv1], grad_enthalpy[icv0], temp[icv1], RoM[icv1],    gamma[icv1], Scalar0, gradScal0, dpress_dscal0, kineFA,
                         mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kineFA, vel[icv1], diffScal, DiffTerm, 
                         area, nVec, smag, nVec, alpha, nScal);  /* <- use nVec here instead of sVec, to avoid inaccurate correction */

              for (int i = 0; i < 5+nScal; i++)
                rhs[icv0][i] -= ViscousFlux[i];
              
              if (flagImplicit)
              {
                for (int i = 0; i < 5+nScal; i++)
                for (int j = 0; j < 5+nScal; j++)
                  A[noc00][i][j] += A0[i][j]; // + A1[i][j];   // !!!! Apl + Ami for Neumann boundary condition: Qr = Ql -> dF/dQ = dF/dQr + dF/dQl
              }
            }            
          }
        }
        
        // .............................................................................................
        // WALL BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "WALL")
        {
          // Set DiffTerm flag for scalars to 0 if Neumann BC for viscous flux
          for (int iScal = 0; iScal < nScal; iScal++)
          {
            double dummy;
            string scalName(scalarTranspEqVector[iScal].getName());
            if (!(scalarZoneIsHook(zone->getName(), scalName) || scalarZoneIsDirichlet(dummy, zone->getName(), scalName)))
              DiffTerm[iScal] = 0.0;
          }

          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];
            assert( icv0 >= 0 );
            int noc00 = nbocv_i[icv0];            // icv0's diagonal

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);

            for (int iScal = 0; iScal < nScal; iScal++)
            {
              ScalCV0[iScal]  = scalarTranspEqVector[iScal].phi[icv0];
              Scalar0[iScal]  = scalarTranspEqVector[iScal].phi[icv1];
            }
            
            double kineCV0 = 0.0;
            double kineFA = 0.0;
            if (kine_Index > -1)
              kineCV0 = scalarTranspEqVector[kine_Index].phi[icv0];
                        
            // Euler flux
            calcEulerFluxCoupled_HLLC(EulerFlux,
                     rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
                     rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
                     area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");
//                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

            for (int i = 0; i < 5+nScal; i++)
              rhs[icv0][i] -= EulerFlux[i];

            if (flagImplicit)
            {
              for (int iScal = 0; iScal < nScal; iScal++)
                if (scalarTranspEqVector[iScal].dpress_dphi != NULL)  
                  dpress_dscal0[iScal] = scalarTranspEqVector[iScal].dpress_dphi[icv0];

              calcEulerFluxMatricesCoupled_HLLC(Apl, NULL,
                       rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, dpress_dscal0, kineFA,
                       rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, dpress_dscal0, kineFA,
                       area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");
//                       area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

              for (int i = 0; i < 5+nScal; i++)
              for (int j = 0; j < 5+nScal; j++)
                A[noc00][i][j] += Apl[i][j];
            }
            
            // Viscous flux
            if (mu_ref > 0.0)
            {
              double sVec[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
              double smag = fabs(vecDotVec3d(sVec, nVec));                // project sVec to wall face normal
              double alpha = 1.0;
              
              for (int iScal = 0; iScal < nScal; iScal++)
              {
                diffScal[iScal] = scalarTranspEqVector[iScal].diff[ifa];
                for (int i = 0; i < 3; i++)
                  gradScal0[iScal][i] = scalarTranspEqVector[iScal].grad_phi[icv0][i];
              }
              
              // Neumann BC for scalars is enforced by setting DiffTerm to 0
              calcViscousFluxCoupled(ViscousFlux, A0, NULL,
                         rho[icv0],    vel[icv0],    grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0],    gamma[icv0], ScalCV0, gradScal0, dpress_dscal0, kineCV0,
                         rho[icv1],    vel[icv1],    grad_u[icv0], enthalpy[icv1], grad_enthalpy[icv0], temp[icv1], RoM[icv1],    gamma[icv1], Scalar0, gradScal0, dpress_dscal0, kineFA,
                         mul_fa[ifa], 0.0, lamOcp_fa[ifa], kineFA, vel[icv1], diffScal, DiffTerm, 
                         area, nVec, smag, nVec, alpha, nScal);  /* <- use nVec here instead of sVec, to avoid inaccurate correction */

              for (int i = 0; i < 5+nScal; i++)
                rhs[icv0][i] -= ViscousFlux[i];
              
              if (flagImplicit)
              {
                for (int i = 0; i < 5+nScal; i++)
                for (int j = 0; j < 5+nScal; j++)
                  A[noc00][i][j] += A0[i][j];
              }
            }
          }
          
          // Set back DiffTerm flag for scalars to original setting
          for (int iScal = 0; iScal < nScal; iScal++)
            DiffTerm[iScal] = (double)scalarTranspEqVector[iScal].diffTerm;
        }
        
        // .............................................................................................
        // OTHER BOUNDARY CONDITIONS (HOOK, DIRICHLET(CBC), ...)
        // .............................................................................................
        else
        {
          // Set DiffTerm flag for scalars to 0 if Neumann BC for viscous flux
          for (int iScal = 0; iScal < nScal; iScal++)
          {
            double dummy;
            string scalname(scalarTranspEqVector[iScal].getName());
            if (!(scalarZoneIsHook(zone->getName(), scalname) || scalarZoneIsDirichlet(dummy, zone->getName(), scalname)))
              DiffTerm[iScal] = 0.0;
          }

          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];
            assert( icv0 >= 0 );
            int noc00 = nbocv_i[icv0]; // icv0's diagonal

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);

            double rho0 = rho[icv0];
            double u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
            double p0 = press[icv0];
            double T0 = temp[icv0];
            double h0 = enthalpy[icv0];
            double gam0 = gamma[icv0];
            double R0 = RoM[icv0];
            double kineCV0 = 0.0;           // cell center
            double kineFA0 = 0.0;           // cell face

            double kineFA1 = 0.0;

            for (int iScal = 0; iScal < nScal; iScal++)
            {
              ScalCV0[iScal] = Scalar0[iScal] = scalarTranspEqVector[iScal].phi[icv0];
              Scalar1[iScal] = scalarTranspEqVector[iScal].phi[icv1];
            }

            if (sndOrder == true)
            {
              double r0[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);

              // left side
              rho0 += vecDotVec3d(r0, grad_rho[icv0]);
#ifdef temp_reconstruction
              T0 += vecDotVec3d(r0, grad_temp[icv0]);
              if ((T0 <= 0.0) || (rho0 <= 0.0))
              {
                T0 = temp[icv0];
                rho0 = rho[icv0];
                myCountReducedOrder++;
              }
#else
              p0 += vecDotVec3d(r0, grad_p[icv0]);
              if ((p0 <= 0.0) || (rho0 <= 0.0))
              {
                p0 = press[icv0];
                rho0 = rho[icv0];
                myCountReducedOrder++;
              }
#endif
              else
              {
                for (int i = 0; i < 3; i++)
                  u0[i] += vecDotVec3d(r0, grad_u[icv0][i]);
       if (stOrderScalar == false) 
                for (int iScal = 0; iScal < nScal; iScal++)
                {
                  if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
                    Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d(r0, scalarTranspEqVector[iScal].grad_rhophi[icv0])) / rho0;
                  else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
                    Scalar0[iScal] += vecDotVec3d(r0, scalarTranspEqVector[iScal].grad_phi[icv0]);
                }
              }

              // calculation of other variables at faces: p/T, h, R, gam
#ifdef temp_reconstruction
              calcThermoProp_T(p0, h0, R0, gam0, rho0, T0, Scalar0, nScal);
#else
              calcThermoProp_p(T0, h0, R0, gam0, rho0, p0, Scalar0, nScal);  
#endif
            }

            if (kine_Index > -1)   // save kine if defined
            {
              kineCV0 = ScalCV0[kine_Index];          // cell center
              kineFA0 = Scalar0[kine_Index];          // cell face
              kineFA1 = Scalar1[kine_Index];
            }

            calcEulerFluxCoupled_HLLC(EulerFlux,
                     rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         Scalar0, kineFA0,
                     rho[icv1],    vel[icv1], press[icv1],  temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1],  Scalar1, kineFA1,
                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

            for (int i = 0; i < 5+nScal; i++)
              rhs[icv0][i] -= EulerFlux[i];
            
            if (flagImplicit)
            {
              for (int iScal = 0; iScal < nScal; iScal++)
                if (scalarTranspEqVector[iScal].dpress_dphi != NULL)  
                  dpress_dscal0[iScal] = scalarTranspEqVector[iScal].dpress_dphi[icv0];

              // HACK:: Q_bfa might depend on Q0 so that the Jacobi matrix Apl should be adapted, not yet implemented
              calcEulerFluxMatricesCoupled_HLLC(Apl, NULL,
                        rho[icv0],    vel[icv0],    press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, dpress_dscal0, kineCV0,
                        rho[icv1],    vel[icv1],    press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar1, dpress_dscal0, kineFA1,
                        area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
//              calcEulerFluxMatricesCoupled_HLLC(Apl, NULL,
//                        rho0,         u0,           p0,          T0,         h0,         R0,           gam0,         Scalar0, dpress_dscal0, kineFA0,
//                        rho[icv1], vel[icv1], press[icv1],  temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1],      Scalar1, dpress_dscal0, kineFA1,
//                        area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

              for (int i = 0; i < 5+nScal; i++)
              for (int j = 0; j < 5+nScal; j++)
                A[noc00][i][j] += Apl[i][j];
            }

            
            // Viscous flux
            if (mu_ref > 0.0)
            {
              double sVec[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
//              double smag = fabs(vecDotVec3d(sVec, nVec));                // project sVec to wall face normal
//              double alpha = 1.0;
              double smag = normVec3d(sVec);
              double alpha = vecDotVec3d(nVec, sVec);
              
              for (int iScal = 0; iScal < nScal; iScal++)
              {
                diffScal[iScal] = scalarTranspEqVector[iScal].diff[ifa];
                for (int i = 0; i < 3; i++)
                  gradScal0[iScal][i] = scalarTranspEqVector[iScal].grad_phi[icv0][i];
              }
              
              // Neumann BC for scalars is enforced by setting DiffTerm to 0
              calcViscousFluxCoupled(ViscousFlux, A0, NULL,
                         rho[icv0],    vel[icv0],    grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0], gamma[icv0], ScalCV0, gradScal0, dpress_dscal0, kineCV0,
                         rho[icv1],    vel[icv1],    grad_u[icv0], enthalpy[icv1], grad_enthalpy[icv0], temp[icv1], RoM[icv1], gamma[icv1], Scalar0, gradScal0, dpress_dscal0, kineFA1,
                         mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kineFA1, vel[icv1], diffScal, DiffTerm, 
//                         area, nVec, smag, nVec, alpha, nScal);  /* <- use nVec here instead of sVec, to avoid inaccurate correction */
                         area, nVec, smag, sVec, alpha, nScal);  

              for (int i = 0; i < 5+nScal; i++)
                rhs[icv0][i] -= ViscousFlux[i];
              
              if (flagImplicit)
              {
                for (int i = 0; i < 5+nScal; i++)
                for (int j = 0; j < 5+nScal; j++)
                  A[noc00][i][j] += A0[i][j];
              }
            }
          }
          
          // Set back DiffTerm flag for scalars to original setting
          for (int iScal = 0; iScal < nScal; iScal++)
            DiffTerm[iScal] = (double)scalarTranspEqVector[iScal].diffTerm;
        }
      }
    }
  }

  // output the number of times switched back to first order at faces
  MPI_Allreduce(&myCountReducedOrder, &CountReducedOrder, 1, MPI_INTEGER, MPI_SUM, mpi_comm);
  if ((CountReducedOrder > 0) && (mpi_rank == 0))
    cout << "Switched back to first order at " << CountReducedOrder << " face(s)" << endl;

  if (Apl != NULL) {freeMem2D(Apl, 0, 5+nScal-1, 0, 5+nScal-1);   Apl = NULL;}
  if (Ami != NULL) {freeMem2D(Ami, 0, 5+nScal-1, 0, 5+nScal-1);   Ami = NULL;}
  if (A0  != NULL) {freeMem2D(A0,  0, 5+nScal-1, 0, 5+nScal-1);   A0  = NULL;}
  if (A1  != NULL) {freeMem2D(A1,  0, 5+nScal-1, 0, 5+nScal-1);   A1  = NULL;}
  
  if (EulerFlux     != NULL) {delete [] EulerFlux;        EulerFlux     = NULL;}
  if (ViscousFlux   != NULL) {delete [] ViscousFlux;      ViscousFlux   = NULL;}
  
  if (Scalar0       != NULL) {delete [] Scalar0;          Scalar0       = NULL;}
  if (Scalar1       != NULL) {delete [] Scalar1;          Scalar1       = NULL;}
  if (ScalCV0       != NULL) {delete [] ScalCV0;          ScalCV0       = NULL;}
  if (ScalCV1       != NULL) {delete [] ScalCV1;          ScalCV1       = NULL;}
  if (gradScal0     != NULL) {delete [] gradScal0;        gradScal0     = NULL;}
  if (gradScal1     != NULL) {delete [] gradScal1;        gradScal1     = NULL;}
  if (dpress_dscal0 != NULL) {delete [] dpress_dscal0;    dpress_dscal0 = NULL;}
  if (dpress_dscal1 != NULL) {delete [] dpress_dscal1;    dpress_dscal1 = NULL;}
  if (diffScal      != NULL) {delete [] diffScal;         diffScal      = NULL;}
  if (ConvTerm      != NULL) {delete [] ConvTerm;         ConvTerm      = NULL;}
  if (DiffTerm      != NULL) {delete [] DiffTerm;         DiffTerm      = NULL;}

  return CountReducedOrder;
}





