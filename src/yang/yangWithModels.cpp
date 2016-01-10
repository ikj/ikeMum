#include "yangWithModels.h"

/*
 * Method: init
 * ------------
 * Initialize YangWithModels member variables
 */
void YangWithModels::init() {
//	qVec = NULL;

	RHSrho   = NULL; registerScalar(RHSrho,  "RHSRHO",  CV_DATA);
	RHSrhou  = NULL; registerVector(RHSrhou, "RHSRHOU", CV_DATA);
	RHSrhoE  = NULL; registerScalar(RHSrhoE, "RHSRHOE", CV_DATA);

	RHSrhoScal = NULL;
}

/*
 * Method: clear()
 * ---------------
 *
 */
void YangWithModels::clear() {
	int nScal = scalarTranspEqVector.size();

	if(lambda != NULL) {
		delete [] lambda; 	lambda = NULL;
	}

//	if(qVec != NULL) {
//		delete [] qVec; 	qVec = NULL;
//	}

	if (RHSrhoScal != NULL && nScal > 0)
		freeMem2D(RHSrhoScal, 0, nScal-1, 0, ncv-1);
}

/*=******************=*
 *  PUBLIC FUNCTIONS  *
 *=******************=*/
/*
 * Method: run
 * -----------
 * Original code = run() in JoeWithModels.h
 */
void YangWithModels::run() {
	if(mpi_rank == 0)
		cout<<endl
			<<classID<<"::run() "<<endl;

	// read mesh or restart file
	initializeFromRestartFile(getStringParam("RESTART"));

	// initialize models
	initialHookScalarRansTurbModel();
	initialHookScalarRansCombModel();
	for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
		updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);

	// initialize Navier-Stokes equations
	initialHook();

	updateCvDataG1G2(rho, REPLACE_DATA);
	updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
	updateCvDataG1G2(rhoE, REPLACE_DATA);

	calcStateVariables();

	// Initialize probe
	initProbes(this);

//	// Update xcvMin and xcvMax
//	for(int icv=0; icv<ncv; ++icv) {
//		for(int i=0; i<3; ++i) {
//			xcvMin[i] = min(xcvMin[i], x_cv[icv][i]);
//			xcvMax[i] = max(xcvMax[i], x_cv[icv][i]);
//		}
//	}

	int nScal = scalarTranspEqVector.size();
	int nVars = 5 + nScal;

	if(mpi_rank==0)
		cout<<endl
		    <<"=========================================="<<endl
		    <<"=  YANG:                                 ="<<endl
		    <<"=    simple stable branch calculator     ="<<endl
		    <<"=    with Natural parameter continuation ="<<endl
		    <<"=  ikj AT stanford DOT edu               ="<<endl
		    <<"=========================================="<<endl;
	MPI_Barrier(mpi_comm);

	// ==========================
	// Read the point information
	// ==========================
	init_pt  = getIntParam("STARTING_PT"); // the point number in the bifurcation curve
	assert(init_pt >= 0);
	final_pt = getIntParam("FINAL_PT");
	assert(final_pt >= init_pt);

	assert(lambda != NULL);
	char paramName[40];
	for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
		sprintf(paramName, "CONTROL_PARAM%d", iEqn);
		lambda[iEqn] = getParam("LAMBDA_INITIAL_0") -> getDouble(paramName);
	}

	dLambda = getDoubleParam("DELTA_LAMBDA");

	resid_energ_th = getDoubleParam("RESID_ENERG_TH");
	assert(resid_energ_th > 0.0);

	if(mpi_rank==0)
		cout<<endl
		    <<"  STARTING_PT = "<<init_pt<<endl
		    <<"  FINAL_PT    = "<<final_pt<<endl
		    <<"  LAMBDA_INITIAL[0] = "<<lambda[0]<<endl
		    <<"  DELTA_LAMBDA      = "<<dLambda<<endl
		    <<"  RESID_ENERG_TH = "<<resid_energ_th<<endl;

	MPI_Barrier(mpi_comm);

	// ===============
	// Run simulations
	// ===============
//	assert(qVec == NULL);
//	qVec = new double [ncv*nVars];
//	for(int i=0; i<ncv*nVars; ++i)
//		qVec[i] = 0.0;
//
//	MPI_Barrier(mpi_comm);

	// Allocate memory
	assert(NcontrolEqns > 0);
//	double* lambda1FromFile = new double [NcontrolEqns];
//	for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
//		lambda1FromFile[iEqn] = 0.0;

	ipt = init_pt;
	while(ipt <= final_pt) {
	    double wtime0 = MPI_Wtime();

	    if(ipt > init_pt)
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
				lambda[iEqn] += dLambda;

		if(mpi_rank==0)
			cout<<endl
			    <<"---------------------------------"<<endl
			    <<" YANG CALC FOR PT = "<<ipt<<endl
			    <<"               LAMBDA = "<<std::setprecision(14)<<lambda[0]<<endl
			    <<"---------------------------------"<<endl
			    <<endl;
		step = 0;
		if(ipt == init_pt)
			step = getIntParam("INIT_STEP_FIRST_PT", "0");
		time = 0.0;
		
		string booleanString = getStringParam("RUN_INIT_PT", "YES");
		std::transform(booleanString.begin(), booleanString.end(), booleanString.begin(), ::tolower);
		bool runInitPt = (booleanString.compare("true")==0 || booleanString.compare("yes")==0
		                   || booleanString.compare("t")==0 || booleanString.compare("y")==0);

		if((ipt == init_pt && runInitPt) || ipt > init_pt) {
			// -------------------------------------------------
			// calculate steady-state soln with time stepping
			// -------------------------------------------------
			string tIntName = getStringParam("TIME_INTEGRATION");

			if (tIntName == "FORWARD_EULER")                 runForwardEuler();
			else if (tIntName == "RK")                       runRK();
			else if (tIntName == "BACKWARD_EULER")           runBackwardEuler();
			else if (tIntName == "BACKWARD_EULER_COUPLED")   runBackwardEulerCoupled();
			else if (tIntName == "BDF2")                     runBDF2();
			else {
				if (mpi_rank == 0) {
					cerr << "ERROR: wrong time integration scheme specified !" << endl;
					cerr << "available integration schemes are: FORWARD_EULER, RK, BACKWARD_EULER, BDF2" << endl;
				}
				throw(YANG_ERROR_CODE);
			}

			// -------------------------------------------------
			// post-processing: write QoI's and tecplot output
			// -------------------------------------------------
			if(mpi_rank==0)
				cout<<"> Writing QoI"<<endl;
			writeQoIOnFile(ipt, filenameQoI, (ipt==init_pt));

			writeData2(ipt, step, true);

			char restartFilename[50];
			sprintf(restartFilename, "restart.pt%04d.step%05d.out", ipt, step);
			writeRestart(restartFilename);
		} else {
			if(mpi_rank==0)
				cout<<"> Writing QoI"<<endl;
			writeQoIOnFile(ipt, filenameQoI, (ipt==init_pt));

			if (initial_flowfield_output.compare("YES"))
				writeData2(ipt, step, true);
		}

		double wtime1 = MPI_Wtime();
		if(mpi_rank==0)
			cout<<endl
			    <<"> PT = "<<ipt<<" is done: wall time [sec] = "<<wtime1-wtime0<<endl;

		++ipt;

		if(lastResidualRhoE > resid_energ_th) {
			if(mpi_rank == 0)
				cerr<<"ERROR in "<<classID<<"::run(): Simulation cannot converged in "<<nsteps<<" steps: residual[4]="
				    <<Residual[4]<<" > RESID_ENERG_TH="<<resid_energ_th<<endl;
			throw(YANG_ERROR_CODE);
		}
		if(lastResidualRhoE == 0.0 && mpi_rank == 0)
			cout<<"WARNING in "<<classID<<"::run(): lastResidualRhoE is zero. Are you using an integration sheme other than Backward Euler?"<<endl
			    <<"Please check if it is updated at the end of the integration scheme"<<endl;
	}
}

/*
 * Method: setBC1D
 * ---------------
 * Original code = setBC1D in IkeWithModels_AD.cpp and setBC() in JoeWithModels.cpp
 * In some situations, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
 *
 * Note: Currently empty!!!!
 */
void YangWithModels::setBC1D(const int ifa, FaZone* zone) {
	/* empty */
}
//void YangWithModels::setBC1D(const int ifa, FaZone* zone) {
//	int bc_err = 0;
//
//	Param* param;
//
//	double *temp     = UgpWithCvCompFlow::temp;
//	double (*velUgp)[3] =  UgpWithCvCompFlow::vel;  // Note: if the name is just "vel", the compiler is confused with the adouble version vel
//	double *press    = UgpWithCvCompFlow::press;
//	double *sos      = UgpWithCvCompFlow::sos;
//	double *enthalpy = UgpWithCvCompFlow::enthalpy;
//	double *gamma    = UgpWithCvCompFlow::gamma;
//	double *RoM      = UgpWithCvCompFlow::RoM;
//
//	for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
//		if (zone->getKind() == FA_ZONE_BOUNDARY) {
//			Param *param;
//
//			if (getParam(param, zone->getName())) {
//				// .............................................................................................
//				// HOOK BOUNDARY CONDITION
//				// .............................................................................................
//				if (param->getString() == "HOOK") {
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//						if(ifa>=nfa && ifa<nfa_b2) { // Special treatment for ggf faces
//							boundaryHook1D(ifa, temp, velUgp, press, zone);
//
//							setScalarBC1D(ifa);
//							ComputeBCProperties1D_T(ifa);
//
//							// update density at boundary using EOS
//							if((ifa>= 0 && ifa<nfa_b) || (ifa>=nfa && ifa<nfa_b2)) {
//								int icv1 = cvofa[ifa][1];
//								rho[icv1] = press[icv1]/(RoM[icv1]*temp[icv1]);
//								sos[icv1] = sqrt(gamma[icv1]*press[icv1]/RoM[icv1]);
//							}
//						}
//					}
//				}
//				// .............................................................................................
//				// CBC BOUNDARY CONDITION
//				// .............................................................................................
//				else if (param->getString() == "CBC") {
//					double u_bc[3], T_bc, p_bc;
//
//					for (int i=0; i<3; i++)
//						u_bc[i] = param->getDouble(i+2);
//					T_bc = param->getDouble(5);
//					p_bc = param->getDouble(6);
//
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//
//						if(ifa>=nfa && ifa<nfa_b2) { // Special treatment for ggf faces
//							int icv0 = cvofa[ifa][0];
//							int icv1 = cvofa[ifa][1];
//
//							// b2 boundary face, between ghost 1 and fake 2...
//							assert( (ifa >= nfa)&&(ifa < nfa_b2) );
//							assert((icv0 >= ncv)&&(icv0 < ncv_g)&&(icv1>=ncv_ggf)&&(icv1<ncv_ggff));
//
//							double nVec[3];
//							double area = normVec3d(nVec, fa_normal[ifa]);
//
//							if (vecDotVec3d(u_bc, nVec) > 0.0) {  // outlet
//								double velMagN = vecDotVec3d(velUgp[icv0], nVec);
//								double mach = fabs(velMagN)/sos[icv0];
//
//								temp[icv1] = temp[icv0];
//								for (int i=0; i<3; i++)
//									velUgp[icv1][i] = velUgp[icv0][i];
//
//								if (mach >= 1.0)  press[icv1] = press[icv0];
//								else              press[icv1] = p_bc;
//							} else {     // inlet
//								temp[icv1] = T_bc;
//								for (int i=0; i<3; i++)
//									velUgp[icv1][i] = u_bc[i];
//								press[icv1] = p_bc;
//							}
//
//							setScalarBC1D(ifa);
//							ComputeBCProperties1D_T(ifa);
//
//							// update density at boundary using EOS
//							if((ifa>= 0 && ifa<nfa_b) || (ifa>=nfa && ifa<nfa_b2)) {
//								int icv1 = cvofa[ifa][1];
//								rho[icv1] = press[icv1]/(RoM[icv1]*temp[icv1]);
//								sos[icv1] = sqrt(gamma[icv1]*press[icv1]/RoM[icv1]);
//							}
//						}
//					}
//				}
//				// .............................................................................................
//				// CBC SUBSONIC INLET BOUNDARY CONDITION
//				// .............................................................................................
//				else if (param->getString() == "CBC_SUBSONIC_INLET") {
//					double angleU[3], Ttot, htot, ptot;
//
//					for (int i=0; i<3; i++)
//						angleU[i] = param->getDouble(i+2);
//					Ttot = param->getDouble(5);
//					ptot = param->getDouble(6);
//
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//
//						if(ifa>=nfa && ifa<nfa_b2) { // Special treatment for ggf faces
//							int icv0 = cvofa[ifa][0];
//							int icv1 = cvofa[ifa][1];
//							//assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));
//
//							{
//								double u[3] = {velUgp[icv0][0], velUgp[icv0][1], velUgp[icv0][2]};
//								double wPow2 = vecDotVec3d(u, u);            // velocity squared
//								double velMag = sqrt(wPow2);                   // apply angle to extrapolated velocity
//								for (int i=0; i<3; i++)
//									velUgp[icv1][i] = angleU[i]*velMag;
//								temp[icv1] = Ttot;                              // save temporary total temperature (to compute then total enthalpy)
//							}
//
//							setScalarBC1D(ifa);
//							ComputeBCProperties1D_T(ifa);                  // total enthalpy from total temperature ?
//
//							{
//								double wPow2 = vecDotVec3d(velUgp[icv1], velUgp[icv1]);         // velocity squared
//								enthalpy[icv1] -= 0.5* wPow2;                               // static enthalpy
//							}
//
//							ComputeBCProperties1D_H(ifa);                  // static temperature and thermo properties from static enthalpy ?
//
//							// Assumes isentropic relations to determine static pressure (= constant cp)
//							// At first approximation ok, but could be improved; should for now be considered in defining p_bc
//							press[icv1] = ptot*pow(temp[icv1]/Ttot, gamma[icv1]/(gamma[icv1]-1.0));
//
//							// update density at boundary using EOS
//							if((ifa>= 0 && ifa<nfa_b) || (ifa>=nfa && ifa<nfa_b2)) {
//								int icv1 = cvofa[ifa][1];
//								rho[icv1] = press[icv1]/(RoM[icv1]*temp[icv1]);
//								sos[icv1] = sqrt(gamma[icv1]*press[icv1]/RoM[icv1]);
//							}
//						}
//					}
//				}
//				// .............................................................................................
//				// CBC SUBSONIC OUTLET BOUNDARY CONDITION
//				// .............................................................................................
//				else if (param->getString() == "CBC_SUBSONIC_OUTLET") {
//					double p_bc = param->getDouble(2);
//
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//
//						if(ifa>=nfa && ifa<nfa_b2) { // Special treatment for ggf faces
//							int icv0 = cvofa[ifa][0];
//							int icv1 = cvofa[ifa][1];
//							//assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));
//
//							// Assumes that the temperature is constant across outlet boundary (extrapolation) while p_bc is prescribed
//							press[icv1] = p_bc;
//							for (int i=0; i<3; i++)
//								velUgp[icv1][i] = velUgp[icv0][i];
//
//							temp[icv1] = temp[icv0];
//
//							setScalarBC1D(ifa);
//							ComputeBCProperties1D_T(ifa);
//
//							// update density at boundary using EOS
//							if((ifa>= 0 && ifa<nfa_b) || (ifa>=nfa && ifa<nfa_b2)) {
//								int icv1 = cvofa[ifa][1];
//								rho[icv1] = press[icv1]/(RoM[icv1]*temp[icv1]);
//								sos[icv1] = sqrt(gamma[icv1]*press[icv1]/RoM[icv1]);
//							}
//						}
//					}
//				}
//				// .............................................................................................
//				// SYMMETRY BOUNDARY CONDITION
//				// .............................................................................................
//				else if (param->getString() == "SYMMETRY") {
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//
//						if(ifa>=nfa && ifa<nfa_b2) { // Special treatment for ggf faces
//							int icv0 = cvofa[ifa][0];
//							int icv1 = cvofa[ifa][1];
//
//							double nVec[3];
//							double area = normVec3d(nVec, fa_normal[ifa]);
//
//							// flip u, APPROXIMATION ---> take velocity at the cell center
//							double u0[3] = {velUgp[icv0][0], velUgp[icv0][1], velUgp[icv0][2]};
//							double un = vecDotVec3d_AD(nVec, u0);
//							for (int i = 0; i < 3; i++)
//								velUgp[icv1][i] = u0[i] - 1.0*un*nVec[i];
//							//assert(fabs(vecDotVec3d(vel_bfa[ifa], nVec)) < 1.0e-10);
//
//							temp[icv1]  = temp[icv0];
//							press[icv1] = press[icv0];
//
//							setScalarBC1D(ifa);
//							ComputeBCProperties1D_T(ifa);
//
//							// update density at boundary using EOS
//							if((ifa>= 0 && ifa<nfa_b) || (ifa>=nfa && ifa<nfa_b2)) {
//								int icv1 = cvofa[ifa][1];
//								rho[icv1] = press[icv1]/(RoM[icv1]*temp[icv1]);
//								sos[icv1] = sqrt(gamma[icv1]*press[icv1]/RoM[icv1]);
//							}
//						}
//					}
//				}
//				// .............................................................................................
//				// NEUMANN BOUNDARY CONDITION
//				// .............................................................................................
//				else if (param->getString() == "NEUMANN") {
//					if ((first)&&(mpi_rank == 0))
//						cout << "Applying NEUMANN             to zone: "<< zone->getName() << endl;
//
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//
//						if(ifa>=nfa && ifa<nfa_b2) { // Special treatment for ggf faces
//							int icv0 = cvofa[ifa][0];
//							int icv1 = cvofa[ifa][1];
//
//							for (int i = 0; i < 3; i++)
//								velUgp[icv1][i] = velUgp[icv0][i];
//
//							temp[icv1] = temp[icv0];
//							press[icv1] = press[icv0];
//
//							setScalarBC1D(ifa);
//							ComputeBCProperties1D_T(ifa);
//
//							// update density at boundary using EOS
//							if((ifa>= 0 && ifa<nfa_b) || (ifa>=nfa && ifa<nfa_b2)) {
//								int icv1 = cvofa[ifa][1];
//								rho[icv1] = press[icv1]/(RoM[icv1]*temp[icv1]);
//								sos[icv1] = sqrt(gamma[icv1]*press[icv1]/RoM[icv1]);
//							}
//						}
//					}
//				}
//				// .............................................................................................
//				// WALL BOUNDARY CONDITION
//				// .............................................................................................
//				else if (param->getString() == "WALL") {
//					int i=0;
//					double T_bc = 0.0;
//					if ((i = param->findString("TEMP")) != 0)
//						T_bc = param->getDouble(i+1);
//
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//
//						if(ifa>=nfa && ifa<nfa_b2) { // Special treatment for ggf faces
//							int icv0 = cvofa[ifa][0];
//							int icv1 = cvofa[ifa][1];
//							//assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));
//
//							for (int i = 0; i < 3; i++)
//								velUgp[icv1][i] = 0.0;
//
//							if (T_bc > 0.0)   temp[icv1] = T_bc;           // wall temperature
//							else              temp[icv1] = temp[icv0];     // adiabatic wall
//
//							press[icv1] = press[icv0];                      // Assumes zero pressure gradient at the wall
//
//							setScalarBC1D(ifa);
//							ComputeBCProperties1D_T(ifa);
//
//							// update density at boundary using EOS
//							if((ifa>= 0 && ifa<nfa_b) || (ifa>=nfa && ifa<nfa_b2)) {
//								int icv1 = cvofa[ifa][1];
//								rho[icv1] = press[icv1]/(RoM[icv1]*temp[icv1]);
//								sos[icv1] = sqrt(gamma[icv1]*press[icv1]/RoM[icv1]);
//							}
//						}
//					}
//				}
//				// .............................................................................................
//				// OTHER BOUNDARY CONDITIONS
//				// .............................................................................................
//				else {
//					if (mpi_rank == 0)
//						cerr << "Error: cannot recognize bc for zone: " << zone->getName() << " bc: " << param->getString() << endl;
//					bc_err = 1;
//				}
//			} else {
//				if (mpi_rank == 0)
//					cerr << "Error: no bc set for: "<< zone->getName() << endl;
//				bc_err = 1;
//			}
//		}
//
//	if (bc_err != 0)
//		throw(YANG_ERROR_CODE);
//}

/****************************
 * RUBENS SPECIFIC METHODS
 ****************************/
///*
// * Method: findMatchedIndex
// * ------------------------
// *
// */
//int YangWithModels::findMatchedIndex(int& foundIndex, const int icv, const int nCVs, const double (*x_cv_eigen)[3], const double epsil) {
//	int count = 0;
//	foundIndex = -1;
//	for(int index =0; index<nCVs; ++index) {
//		if( fabs(x_cv_eigen[index][0]-x_cv[icv][0])<epsil && fabs(x_cv_eigen[index][1]-x_cv[icv][1])<epsil && fabs(x_cv_eigen[index][2]-x_cv[icv][2])<epsil ) {
//			foundIndex = index;
//			count++;
//		}
//	}
//	return count;
//}
//
///*
// * Method: readPsALCdumpedDataSerial
// * ---------------------------------
// * Read field data from file
// * Original code: IkeWithModels::readPsALCdumpedDataSerial()
// */
//void YangWithModels::readPsALCdumpedDataSerial(const char filename[], int &stepFromFile, double* lambda1FromFile, double* qVec, const int NcontrolEqns) {
//	// Check some possible errors
//	for(int i=0; i<3; ++i)
//		if(xcvMax[i] < xcvMin[i]) {
//			cerr<<"ERROR in YangWithModels::readPsALCdumpedDataSerial(): mpi_rank == "<<mpi_rank<<" : xcvMax["<<i<<"] < xcvMin["<<i<<"]"<<endl;
//			throw(RUBENS_ERROR_CODE);
//		}
//	assert(qVec != NULL);
//	assert(lambda1FromFile != NULL);
//
//	// number of variables
//	int nScal = scalarTranspEqVector.size();
//	int m = nScal+5; // number of variables (rho, rhou, rhov, rhow, rhoE, and scalars)
//
//	// variables related to the Mesh
//	const double gridTol = 1.0e-9;
//
//	int *countFound = new int [ncv];
//	for(int icv=0; icv<ncv; ++icv)
//		countFound[icv] = 0;
//
//	// read the file
//	ifstream infile(filename, ios::in | ios::binary);
//	if(infile.is_open()) {
//		int NcontrolEqns_file;
//		double *lambda0 = NULL;
//		double dsTry;
//
//		int mpi_size_file;
//		int *cvora_file = NULL;
//		double *xcvMinArray_file = NULL;
//		double *xcvMaxArray_file = NULL;
//
//		// 1. Read the header part of the file
//		//      Structure of the header part:
//		//              1. step (int)                      2. NcontrolEqns (int)
//		//              3. lambda0 (double*NcontrolEqns)   4. lambda1 (double*NcontrolEqns)
//		//              5. ds (double)
//		//              6. mpi_size (int)                  7. cvora (double*(mpi_size+1)
//		//              8. xMinArray (double*3*mpi_size)   9. xMaxArray (double*3*mpi_size)
//		infile.read(reinterpret_cast<char*>(&stepFromFile), sizeof(int)); // Note: if you simply use the ">>" operator, it will cause some problem due to bugs in a g++ library on linux
//		infile.read(reinterpret_cast<char*>(&NcontrolEqns_file), sizeof(int));
//		if(NcontrolEqns_file != NcontrolEqns) {
//			if(mpi_rank==0)
//				printf("ERROR in YangWithModels::readPsALCdumpedData: NcontrolEqns from file=%d is not equal to given NcontrolEqns=%d\n", NcontrolEqns_file, NcontrolEqns);
//			delete [] countFound;
//			throw(RUBENS_ERROR_CODE);
//		}
//		if(NcontrolEqns_file > 0) {
//			lambda0 = new double [NcontrolEqns_file];
//		}
//		double dummyDouble;
//		for(int iEqn=0; iEqn<NcontrolEqns_file; ++iEqn) {
//			infile.read(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
//			lambda0[iEqn] = dummyDouble;
//		}
//		for(int iEqn=0; iEqn<NcontrolEqns_file; ++iEqn) {
//			infile.read(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
//			lambda1FromFile[iEqn] = dummyDouble;
//		}
//		infile.read(reinterpret_cast<char*>(&dsTry), sizeof(double));
//		infile.read(reinterpret_cast<char*>(&mpi_size_file), sizeof(int));
//		cvora_file = new int [mpi_size_file+1];
//		infile.read(reinterpret_cast<char*>(cvora_file), sizeof(int)*(mpi_size_file+1));
//		xcvMinArray_file = new double [mpi_size_file*3];
//		infile.read(reinterpret_cast<char*>(xcvMinArray_file), sizeof(double)*(mpi_size_file*3));
//		xcvMaxArray_file = new double [mpi_size_file*3];
//		infile.read(reinterpret_cast<char*>(xcvMaxArray_file), sizeof(double)*(mpi_size_file*3));
//
//		if(mpi_rank==0) {
//			printf("\nReading \'%s\': step=%d, lambda1=", filename, stepFromFile);
//			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
//				printf("%.8f, ", lambda1FromFile[iEqn]);
//			printf("ds=%.8f, total number of CVs=%d, mpi_size of the file=%d\n", dsTry, cvora_file[mpi_size_file], mpi_size_file);
//			cout<<"  Grid tolerance = "<<gridTol<<endl;
//			if(RUBENS_DEBUG_LEVEL > 1) {
//				cout<<"  XCV_BOX = "<<endl;
//				for(int impi=0; impi<mpi_size_file; ++impi)
//					printf("    %3d: min=(%.3f, %.3f, %.3f)   max=(%.3f, %.3f, %.3f)\n", impi,
//							xcvMinArray_file[impi*3],xcvMinArray_file[impi*3+1],xcvMinArray_file[impi*3+2],
//							xcvMaxArray_file[impi*3],xcvMaxArray_file[impi*3+1],xcvMaxArray_file[impi*3+2]);
//			}
//		}
//
//		assert(stepFromFile >= 0);
//		assert(cvora_file[mpi_size_file] == cvora[mpi_size]);
//		for(int mpi_file=0; mpi_file<mpi_size_file; ++mpi_file)
//			assert(cvora_file[mpi_file+1] > 0);
//
//		// 2. Read the body
//		//       Structure of the body part:
//		//            For each mpi,
//		//              1. x_cv (double*3*ncv)
//		//              2. qVec (double*m*ncv)
//		int ncv_file;
//		double (*x_cv_file)[3] = NULL;
//		double *qVec_file = NULL;
//
//		for(int mpi_file=0; mpi_file<mpi_size_file; ++mpi_file) {
//			// read data dumped by a CPU
//			ncv_file = cvora_file[mpi_file+1]-cvora_file[mpi_file];
//
//			assert(x_cv_file==NULL);	x_cv_file = new double [ncv_file][3];
//			assert(qVec_file==NULL);	qVec_file = new double [ncv_file*m];
//
//			for(int icv=0; icv<ncv_file; ++icv) {
//				for(int i=0; i<3; ++i)
//					infile.read(reinterpret_cast<char*>(&x_cv_file[icv][i]), sizeof(double));
//			}
//			infile.read(reinterpret_cast<char*>(qVec_file), sizeof(double)*(ncv_file*m));
//
//			// find the matched data
//			bool matchedBox = true;
//
//			double xcvMin_file[3] = {xcvMinArray_file[mpi_file*3], xcvMinArray_file[mpi_file*3+1], xcvMinArray_file[mpi_file*3+2]};
//			double xcvMax_file[3] = {xcvMaxArray_file[mpi_file*3], xcvMaxArray_file[mpi_file*3+1], xcvMaxArray_file[mpi_file*3+2]};
//
//			if(xcvMin_file[0]>xcvMax[0] || xcvMin_file[1]>xcvMax[1] || xcvMin_file[2]>xcvMax[2])
//				matchedBox = false;
//			if(xcvMax_file[0]<xcvMin[0] || xcvMax_file[1]<xcvMin[1] || xcvMax_file[2]<xcvMin[2])
//				matchedBox = false;
//
//			if(matchedBox){
//				for(int icv=0; icv<ncv; ++icv) {
//					int icv_file;
//					int count = findMatchedIndex(icv_file, icv, ncv_file, x_cv_file, gridTol);
//					if(count>0) {
//						for(int i=0; i<m; ++i)
//							qVec[m*icv+i] = qVec_file[m*icv_file+i];
//						countFound[icv] += count;
//					}
//				}
//			}
//			delete [] x_cv_file;	x_cv_file = NULL;
//			delete [] qVec_file;	qVec_file = NULL;
//		}
//
//		// 3. Read the foot
//		int dummyInt;
//		infile.read(reinterpret_cast<char*>(&dummyInt), sizeof(int));
//		if(dummyInt != EOF_ERROR_CHECK_CODE) {
//			if(mpi_rank==0)
//				cout<<"ERROR! YangWithModels::readPsALCdumpedDataSerial(): file does not end with EOF_ERROR_CHECK_CODE="<<EOF_ERROR_CHECK_CODE<<endl;
//			assert(false);
//		}
//
//		// 4. check if each CV found a matched CV
//		for (int icv=0; icv<ncv; ++icv) {
//			if(countFound[icv]==0) {
//				printf("  Cannot find a matched point for icv=%d (%.6f,%.6f,%.6f) at mpi_rank=%d \n", icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2],mpi_rank);
//				throw(-1);
//			} else if(countFound[icv]>1) {
//				printf("  Too many matched points(%d) for icv=%d (%.6f,%.6f,%.6f) at mpi_rank=%d \n", countFound[icv],icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2],mpi_rank);
//				throw(-1);
//			}
//		}
//
//		delete [] lambda0;
//
//		delete [] cvora_file;
//		delete [] xcvMinArray_file;
//		delete [] xcvMaxArray_file;
//	}  else {
//		if(mpi_rank==0)
//			printf("YangWithModels::readPsALCdumpedDataSerial(): Cannot read \'%s\' \n", filename);
//		throw(-1);
//	}
//	infile.close();
//
//	delete [] countFound; 	countFound = NULL;
//
//	MPI_Barrier(mpi_comm);
//}

