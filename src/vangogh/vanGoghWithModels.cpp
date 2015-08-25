#include "vanGoghWithModels.h"

/*
 * Method: init
 * ------------
 * Initialize VanGoghWithModels member variables
 */
void VanGoghWithModels::init() {
	// initial field (from a restart file)
	rho_init = NULL; 	registerScalar(rho_init, "INIT_RHO", CV_DATA);
	rhou_init = NULL; 	registerVector(rhou_init, "INIT_RHOU", CV_DATA);
	rhoE_init = NULL; 	registerScalar(rhoE_init, "INIT_RHOE", CV_DATA);
	// NOTE: scalars should be registered in the child class

	// data on the unstable branch (from Q0_PT000**.bin)
	rho_unstable = NULL; 	registerScalar(rho_unstable, "UNSTABLE_RHO", CV_DATA);
	rhou_unstable = NULL;	registerVector(rhou_unstable, "UNSTABLE_RHOU", CV_DATA);
	rhoE_unstable = NULL; 	registerScalar(rhoE_unstable, "UNSTABLE_RHOE", CV_DATA);
	// NOTE: scalars should be registered in the child class

	// -------------------------
	// Perturbations
	// -------------------------
	// perturbation
	array_perturb = NULL;
	wdFace = NULL;

	// -------------------------
	// Eigenvecs: Tecplot output
	// -------------------------
	// least-stable global modes
	rho_Dct_1stReal  = NULL; 	registerScalar(rho_Dct_1stReal,  "DCT1ST_RHO", CV_DATA);
	rhou_Dct_1stReal = NULL;	registerVector(rhou_Dct_1stReal, "DCT1ST_RHOU", CV_DATA);
	rhoE_Dct_1stReal = NULL;	registerScalar(rhoE_Dct_1stReal, "DCT1ST_RHOE", CV_DATA);
	// NOTE: scalars should be registered in the child class

	// first adjoint global modes
	rho_Adj_1stReal  = NULL;	registerScalar(rho_Adj_1stReal,  "ADJ1ST_RHO", CV_DATA);
	rhou_Adj_1stReal = NULL;	registerVector(rhou_Adj_1stReal, "ADJ1ST_RHOU", CV_DATA);
	rhoE_Adj_1stReal = NULL;	registerScalar(rhoE_Adj_1stReal, "ADJ1ST_RHOE", CV_DATA);
	// NOTE: scalars should be registered in the child class
}

/*
 * Method: clear()
 * ---------------
 * Clear VanGoghWithModels member arrays
 */
void VanGoghWithModels::clear() {
	if(lambda != NULL) {
		delete [] lambda;	lambda = NULL;
	}

	if(array_perturb != NULL) {
		delete [] array_perturb;  array_perturb = NULL;
	}

	if(wdFace != NULL) {
		delete [] wdFace;         wdFace = NULL;
	}
}

/*=******************=*
 *  PUBLIC FUNCTIONS  *
 *=******************=*/

/*
 * Method: run
 * -----------
 * Original code = run() in JoeWithModels.h
 */
void VanGoghWithModels::run() {
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

	initProbes(this);

	updateCvDataG1G2(rho, REPLACE_DATA);
	updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
	updateCvDataG1G2(rhoE, REPLACE_DATA);

	calcStateVariables();

	// Update xcvMin and xcvMax
	for(int icv=0; icv<ncv; ++icv) {
		for(int i=0; i<3; ++i) {
			xcvMin[i] = min(xcvMin[i], x_cv[icv][i]);
			xcvMax[i] = max(xcvMax[i], x_cv[icv][i]);
		}
	}

	int nScal = scalarTranspEqVector.size();

	if(mpi_rank==0)
		cout<<endl
		    <<"============================="<<endl
		    <<"=  VAN GOGH                 ="<<endl
		    <<"=  ikj AT stanford DOT edu  ="<<endl
		    <<"============================="<<endl;
	MPI_Barrier(mpi_comm);

	// ===================================
	// Read direct modes and adjoint modes
	// ===================================
	int ptOnSCurve;

	string tempString = getStringParam("EIGENPAIRS_SOURCE", "REMBRANDT");
	if(mpi_rank==0)
		cout<<"> EIGENPAIRS_SOURCE = "<<tempString<<endl;
	if(tempString.compare("MATLAB") == 0) {
		ptOnSCurve = getIntParam("PT_ON_S_CURVE", "-1");
//		readLinSystemMATLAB(ptOnSCurve, rho_1stMode, rhou_1stMode, rhoE_1stMode, kine_1stMode, omega_1stMode,
//				rho_adj, rhou_adj, rhoE_adj, kine_adj, omega_adj);
	} else if(tempString.compare("REMBRANDT") == 0) {
		assert(NcontrolEqns > 0);
		double* lambda_file = new double [NcontrolEqns];
		string filename = getStringParam("EIGENPAIRS_FILENAME", "NONE");

		assert(DirectEvecs.empty()); 	DirectEvecs.resize(1);  // Note: We will get only the first vector
		assert(AdjointEvecs.empty()); 	AdjointEvecs.resize(1);  // Note: We will get only the first vector

		if(mpi_rank==0)
			cout<<"> Reading EIGENPAIRS_FILENAME = "<<filename<<endl;
		readEigenPairsRembrandtSerial(filename, ptOnSCurve, lambda_file, NcontrolEqns);
		assert(!DirectEvecs[0].empty());
		assert(!AdjointEvecs[0].empty());

		if(fabs((lambda_file[0]-lambda[0]) / max(lambda[0], MACHINE_EPS)) > 1.0e-6) {
			if(mpi_rank==0)
				cerr<<"ERROR in "<<classID<<"::run(): lambda[0]="<<lambda[0]<<" is different from lambda_file[0]="<<lambda_file[0]<<endl;
			throw(VANGOGH_ERROR_CODE);
		}
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
			lambda[iEqn] = lambda_file[iEqn]; // Update lambda from EIGENPAIRS_SOURCE

		delete [] lambda_file;
	} else {
		if(mpi_rank==0)
			cerr<<"ERROR in "<<classID<<"::run(): EIGENPAIRS_SOURCE="<<tempString<<" is not supported"<<endl;
		throw(VANGOGH_ERROR_CODE);
	}
	MPI_Barrier(mpi_comm);

	// Check these out again!
	assert(!DirectEvecs[0].empty());
	assert(!AdjointEvecs[0].empty());

	updateEigenVecTecplotNS();
	updateEigenVecTecplotScalars();

	// =====================================
	// Read the point on the unstable branch
	// =====================================
	string filenameQ1 = getStringParam("UNSTABLE_Q1_FILENAME", "NONE");
	double* qVecTemp = new double [ncv*(5+nScal)];
	if(mpi_rank==0) cout<<"> UNSTABLE_Q1_FILENAME = "<<filenameQ1<<endl;
	readPsALCdumpedDataSerial(filenameQ1.c_str(), qVecTemp, NcontrolEqns);

	if(mpi_rank==0)
		cout<<"  UPDATE UNSTABLE VECTORS..."<<endl;
	updateUnstableVecNS(qVecTemp, 5+nScal);
	updateUnstableVecScalars(qVecTemp, 5+nScal);

	MPI_Barrier(mpi_comm);
	delete [] qVecTemp; 	qVecTemp = NULL;

	// =====================
	// Perturbation analysis
	// =====================
	// Save the initial data from the restart file if the user wants to reinitialize the flow field by random perturbation
	string boolString = getStringParam("ALREADY_PERTURBED_RESTART", "NO");
	if(mpi_rank==0) cout<<"> ALREADY_PERTURBED_RESTART = "<<boolString<<endl;
	std::transform(boolString.begin(), boolString.end(), boolString.begin(), ::toupper);
	if(boolString.compare("YES")==0 || boolString.compare("Y")==0 || boolString.compare("TRUE")==0 || boolString.compare("T")==0) {
		itest = getIntParam("CASE_NUM", "0");
		ntests = itest+1;
	} else {
		itest = 0;
		ntests = getIntParam("NTESTS", "1");
	}

	storeInitRestartNS();
	storeInitRestartScalars();

	// Get the parameters for the perturbations
	getPerturbParams(perturbParams);

	// Allocate memory: array_perturb & wdFace
	assert(array_perturb == NULL);  	array_perturb = new double [ncv_ggff];
	if(perturbParams.useFiltering) {
		if(mpi_rank==0)
			cout<<"> Allocate memory for wdFace"<<endl;
		assert(wdFace == NULL);  	wdFace = new double [nfa_b2gg]; // wall distance at faces

		calcWallDistanceFace(wdFace);
	}

	bool firstItest = true;
	// run simulations
	while(itest < ntests) {
		if(mpi_rank==0)
			cout<<endl
				<<">> RUNNING "<<itest<<"TH SIMULATION"<<endl
				<<endl;

		step = 0;

		// initialize models
		initialHookScalarRansTurbModel();
		initialHookScalarRansCombModel();
		for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
			updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);

		// initialize Navier-Stokes equations
		initialHook();

		// -------------------------------------------------
		// random perturbation
		// -------------------------------------------------
		// Perturb the field with filtering
		perturbFieldScalars();
		perturbFieldNS();

		updateCvDataG1G2(rho,  REPLACE_DATA);
		updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
		updateCvDataG1G2(rhoE, REPLACE_DATA);
		for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
			updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);

		// -------------------------------------------------
		// run the solver
		// -------------------------------------------------
		// calculate state variables: press, etc.
		calcStateVariables();

		char filenameQoI[30];
		sprintf(filenameQoI, "QoI.txt");
		writeQoIOnFile(itest, filenameQoI, firstItest);

		char filenameResidual[30];
		sprintf(filenameResidual, "residual.case%04d.csv", itest);
		writeResidualOnFile(itest, filenameResidual, firstItest);

		writeData(itest, 0); // Save each case on a different file

		// write restart file
		char restartInitFilename[30];
		sprintf(restartInitFilename, "Restart.case%04d.init", itest);
		writeRestart(restartInitFilename);

		// run solver
		string tIntName = getStringParam("TIME_INTEGRATION");

		if (tIntName == "FORWARD_EULER")                 runForwardEuler();
		else if (tIntName == "RK")                       runRK();
		else if (tIntName == "BACKWARD_EULER")           runBackwardEuler();
		else if (tIntName == "BACKWARD_EULER_COUPLED")   runBackwardEulerCoupled();
		else if (tIntName == "BDF2")                     runBDF2();
		else {
			cerr << "ERROR: wrong time integration scheme specified !" << endl;
			cerr << "available integration schemes are: FORWARD_EULER, RK, BACKWARD_EULER, BDF2" << endl;
		}

		char restartFilename[32];
		sprintf(restartFilename, "Restart.case%04d.%06d.out", itest, step);
		writeRestart(restartFilename);

		++itest;
		firstItest = false;
	}
}

/****************************
 * UTILITY METHODS
 ****************************/

/*
 * Method: calcRMS
 * -------------------
 * Return the root mean square of a field variable (Option: 1. simple  2. volume averaged)
 */
double VanGoghWithModels::calcRMS(const double *vec, bool volumeAvg) {
	double myNormSq = 0.0;
	double NormSq = 0.0;
	double myVolume = 0.0;
	double Volume = 0.0;

	/* calculate 2-Norm */
	for (int i = 0; i < ncv; i++) {
		if(volumeAvg) {
			myNormSq += fabs(vec[i])*fabs(vec[i])*cv_volume[i];
			myVolume += cv_volume[i];
		} else
			myNormSq += fabs(vec[i])*fabs(vec[i]);
	}

	MPI_Allreduce(&myNormSq, &NormSq, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
	MPI_Allreduce(&myVolume, &Volume, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

	double rmsVal;
	if(volumeAvg)
		rmsVal = sqrt(NormSq/Volume);
	else
		rmsVal = sqrt(NormSq/cvora[mpi_size]);

	return rmsVal;
}

/*
 * Method: calcRMS3D
 * -----------------
 * Return the root mean square of a field variable (Option: 1. simple  2. volume averaged)
 */
double VanGoghWithModels::calcRMS3D(const double (*vec)[3], const int dim, bool volumeAvg) {
	double myNormSq = 0.0;
	double NormSq = 0.0;
	double myVolume = 0.0;
	double Volume = 0.0;

	/* calculate 2-Norm */
	for (int icv = 0; icv < ncv; icv++) {
		if(volumeAvg) {
			myNormSq += fabs(vec[icv][dim])*fabs(vec[icv][dim])*cv_volume[icv];
			myVolume += cv_volume[icv];
		} else
			myNormSq += fabs(vec[icv][dim])*fabs(vec[icv][dim]);
	}

	MPI_Allreduce(&myNormSq, &NormSq, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
	MPI_Allreduce(&myVolume, &Volume, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

	double rmsVal;
	if(volumeAvg)
		rmsVal = sqrt(NormSq/Volume);
	else
		rmsVal = sqrt(NormSq/cvora[mpi_size]);

	return rmsVal;
}

/*
 * Method: findMinValue
 * --------------------
 * Find the minimun value of an array and return it
 */
double VanGoghWithModels::findMinValue(double* array, const int size) {
	double myMinVal = 1.0e15;
	for(int i=0; i<size; ++i) {
		if(array[i] < myMinVal) {
			myMinVal = array[i];
		}
	}

	double minVal;
	MPI_Allreduce(&myMinVal, &minVal, 1, MPI_DOUBLE, MPI_MIN, mpi_comm);

	return minVal;
}

/*
 * Method: findMaxalue
 * --------------------
 * Find the maximum value of an array and return it
 */
double VanGoghWithModels::findMaxValue(double* array, const int size) {
	double myMaxVal = -1.0e15;
	for(int i=0; i<size; ++i) {
		if(array[i] > myMaxVal) {
			myMaxVal = array[i];
		}
	}

	double maxVal;
	MPI_Allreduce(&myMaxVal, &maxVal, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);

	return maxVal;
}

/*
 * Method: calcWallDistanceFace
 * ----------------------------
 * calculate wall distance: minimum distance of cv to wall face by looping over all cv's and faces
 */
void VanGoghWithModels::calcWallDistanceFace(double *wd) {
	if (mpi_rank == 0)
		cout << "VanGoghKOm::calcWallDistanceFace() - calc wall distance at the faces" << endl;
	assert(wd!=NULL);

	double wtime, wtime0;
	MPI_Barrier(mpi_comm);
	if (mpi_rank == 0)
		wtime = MPI_Wtime();

	for (int ifa=0; ifa<nfa_b2gg; ifa++)
		wd[ifa] = 1.0e20;

	// count the walls for each rank
	int my_wall_faces = 0;
	for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
		if (zone->getKind() == FA_ZONE_BOUNDARY) {
			Param *param;
			if (getParam(param, zone->getName()))
				if (param->getString() == "WALL") {
					for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
						my_wall_faces++;
				}
		}

	int tot_wall_faces = 0;
	MPI_Allreduce(&my_wall_faces, &tot_wall_faces, 1, MPI_INT, MPI_SUM, mpi_comm);

	if (tot_wall_faces == 0) {
		if (mpi_rank == 0)
			cerr << "ERROR: you have been trying to calculate the wall distance, but there are no walls in you domain ?!?" << endl;
		throw(-1);
	}

	// set up send side
	int send_count[mpi_size], send_displ[mpi_size];
	for (int r=0; r<mpi_size; r++) {
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
		if (zone->getKind() == FA_ZONE_BOUNDARY) {
			Param *param;
			if (getParam(param, zone->getName()))
				if (param->getString() == "WALL")
					for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++) {
						for (int i=0; i<3; i++) {
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
			wall_fa,   recv_count, recv_displ, MPI_DOUBLE, mpi_comm);

	MPI_Alltoallv(my_wall_fa_normal, send_count, send_displ, MPI_DOUBLE,
			wall_fa_n, recv_count, recv_displ, MPI_DOUBLE, mpi_comm);

	delete [] my_wall_fa;
	delete [] my_wall_fa_normal;


	// find minimum wall distance
	for (int ifa0=0; ifa0<nfa_b2gg; ifa0++) {
		double minDist = 1.0e20;
		int ifaSave;

		for (int ifa=0; ifa<tot_wall_faces; ifa++) {
			double dist = sqrt(pow(x_fa[ifa0][0]-wall_fa[3*ifa], 2.0)
					+ pow(x_fa[ifa0][1]-wall_fa[3*ifa+1], 2.0)
					+ pow(x_fa[ifa0][2]-wall_fa[3*ifa+2], 2.0));

			if (dist < minDist) {
				ifaSave = ifa;
				minDist = dist;
			}
		}

		double fa_x[3] = {wall_fa[3*ifaSave], wall_fa[3*ifaSave+1], wall_fa[3*ifaSave+2]};
		double fa_n[3] = {wall_fa_n[3*ifaSave], wall_fa_n[3*ifaSave+1], wall_fa_n[3*ifaSave+2]};
		double n[3], s_half[3];

		double nmag = normVec3d(n, fa_n);
		vecMinVec3d(s_half, fa_x, x_fa[ifa0]);

		double sMag2 = vecDotVec3d(s_half, s_half);
		double dNorm2 = vecDotVec3d(s_half, n)*vecDotVec3d(s_half, n);
		double dTang2 = sMag2-dNorm2;

		double specRadius2 = 0.5*nmag;

		if (dTang2 > specRadius2)   wd[ifa0] = minDist;
		else                        wd[ifa0] = fabs(vecDotVec3d(s_half, n));
	}
	MPI_Barrier(mpi_comm);

	if (mpi_rank == 0) {
		double wtime0 = wtime;
		wtime = MPI_Wtime();
		cout << " > time to compute wall distance [s]: " << wtime - wtime0 << endl;
	}

	delete [] wall_fa;
	delete [] wall_fa_n;
}

/****************************
 * VANGOGH SPECIFIC METHODS
 ****************************/
/*
 * Method: findMatchedIndex
 * ------------------------
 *
 */
int VanGoghWithModels::findMatchedIndex(int& foundIndex, const int icv, const int nCVs, const double (*x_cv_eigen)[3], const double epsil) {
	int count = 0;
	foundIndex = -1;
	for(int index =0; index<nCVs; ++index) {
		if( fabs(x_cv_eigen[index][0]-x_cv[icv][0])<epsil && fabs(x_cv_eigen[index][1]-x_cv[icv][1])<epsil && fabs(x_cv_eigen[index][2]-x_cv[icv][2])<epsil ) {
			foundIndex = index;
			count++;
		}
	}
	return count;
}

/*
 * Method: readLinSystemMATLAB
 * ---------------------------
 * Read eigenmodes and adjoint modes from "dumpedEigen_pt*****.bin" -- generated from a MATLAB script.
 * This methods has less flexibility: it reads only "dumpedEigen_pt*****.bin" and thus ptOnSCurve must be given.
 * This method is modified from the old version (readLinSystem()), but it has a significant flaw!
 */
void VanGoghWithModels::readLinSystemMATLAB(const int ptOnSCurve) {
	if(mpi_rank == 0)
		cout<<"VanGoghWithModels::readLinSystemMATLAB() "<<endl;

	stringstream ss;
	ss<<"dumpedEigen_pt"<<ptOnSCurve<<".bin";

	if(mpi_rank == 0)
		cout<<"ERROR! Not yet developed!"<<endl;

	throw(VANGOGH_ERROR_CODE);
}

/*
 * Method: readEigenPairsRembrandtSerial
 * -------------------------------------
 * Read data dumped by RembrandtWithModels::writeEigenPairsParallel()
 */
void VanGoghWithModels::readEigenPairsRembrandtSerial(string &filename, int &ptOnSCurve_file, double* lambda_file, const int NcontrolEqns) {
	char *filenameArray = new char[filename.size() + 1]; // Since MPI_File_open() doesn't accept "const char *", you cannot use std::string.c_str()
	std::copy(filename.begin(), filename.end(), filenameArray);
	filenameArray[filename.size()] = '\0';
	readEigenPairsRembrandtSerial(filenameArray, ptOnSCurve_file, lambda_file, NcontrolEqns);

	delete [] filenameArray;
}

void VanGoghWithModels::readEigenPairsRembrandtSerial(const char filename[], int &ptOnSCurve_file, double* lambda_file, const int NcontrolEqns) {
	int nScal = scalarTranspEqVector.size();

	vector<complex<double> > Evals;  // Note: 'std::complex' class

	assert(DirectEvecs.size() == 1); // Note: We will get only the first vector
	DirectEvecs[0].resize(5+nScal);
	for(size_t ivar=0; ivar<5+nScal; ++ivar)
		DirectEvecs[0][ivar].resize(ncv);

	assert(AdjointEvecs.size() == 1); // Note: We will get only the first vector
	AdjointEvecs[0].resize(5+nScal);
	for(size_t ivar=0; ivar<5+nScal; ++ivar)
		AdjointEvecs[0][ivar].resize(ncv);

	// Check some possible errors
	for(int i=0; i<3; ++i)
		if(xcvMax[i] < xcvMin[i]) {
			cerr<<"ERROR in VanGoghWithModels::readEigenPairsRembrandtSerial(): mpi_rank == "<<mpi_rank<<" : xcvMax["<<i<<"] < xcvMin["<<i<<"]"<<endl;
			throw(VANGOGH_ERROR_CODE);
		}

	// variables related to the Mesh
	const double gridTol = 1.0e-9;

	int *countFound = new int [ncv];
	for(int icv=0; icv<ncv; ++icv)
		countFound[icv] = 0;

	// Read the file
	int nVars; // number of variables (rho, rhou, rhov, rhow, rhoE, and scalars)
	int nev;

	int mpi_size_file;
	int *cvora_file = NULL;
	double *xcvMinArray_file = NULL;
	double *xcvMaxArray_file = NULL;

	ifstream infile(filename, ios::in | ios::binary);
	if(infile.is_open()) {
		int NcontrolEqns_file;
		// 1. Read the header part of the file
		//      Structure of the header part:
		//              1. step (int)                      2. NcontrolEqns (int)
		//              3. lambda (double*NcontrolEqns)    4. nVars (int)
		//              5. nev (int)                       6. eigenvalues (double*2*nev)
		//              7. mpi_size (int)                  9. cvora (int*(mpi_size+1))
		//              9. xMinArray (double*3*mpi_size)   10. xMaxArray (double*3*mpi_size)
		infile.read(reinterpret_cast<char*>(&ptOnSCurve_file), sizeof(int)); // Note: if you simply use the ">>" operator, it will cause some problem due to bugs in a g++ library on linux
		infile.read(reinterpret_cast<char*>(&NcontrolEqns_file), sizeof(int));
		if(NcontrolEqns_file != NcontrolEqns) {
			if(mpi_rank==0)
				printf("ERROR in VanGoghWithModels::readEigenPairsRembrandtSerial: NcontrolEqns from file=%d is not equal to given NcontrolEqns=%d\n", NcontrolEqns_file, NcontrolEqns);
			delete [] countFound;
			throw(VANGOGH_ERROR_CODE);
		}
		double dummyDouble;
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
			infile.read(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
			lambda_file[iEqn] = dummyDouble;
		}
		infile.read(reinterpret_cast<char*>(&nVars), sizeof(int));
		if(nVars != 5+nScal) {
			if(mpi_rank==0)
				printf("ERROR in VanGoghWithModels::readEigenPairsRembrandtSerial: nVars from file (=%d) is not equal to the number of variables in VanGogh (=%d)\n", nVars, 5+nScal);
			delete [] countFound;
			throw(VANGOGH_ERROR_CODE);
		}
		infile.read(reinterpret_cast<char*>(&nev), sizeof(int));
		Evals.resize(nev);
		double dummyReal, dummyImag;
		for(int iev=0; iev<nev; ++iev) {
			infile.read(reinterpret_cast<char*>(&dummyReal), sizeof(double));
			infile.read(reinterpret_cast<char*>(&dummyImag), sizeof(double));
			Evals[nev] = complex<double>(dummyReal, dummyImag);
		}
		infile.read(reinterpret_cast<char*>(&mpi_size_file), sizeof(int));
		cvora_file = new int [mpi_size_file+1];
		infile.read(reinterpret_cast<char*>(cvora_file), sizeof(int)*(mpi_size_file+1));
		xcvMinArray_file = new double [mpi_size_file*3];
		infile.read(reinterpret_cast<char*>(xcvMinArray_file), sizeof(double)*(mpi_size_file*3));
		xcvMaxArray_file = new double [mpi_size_file*3];
		infile.read(reinterpret_cast<char*>(xcvMaxArray_file), sizeof(double)*(mpi_size_file*3));

		if(mpi_rank==0) {
			printf("\nReading \'%s\': PT=%d, lambda0=", filename, ptOnSCurve_file);
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
				printf("%.8f, ", lambda_file[iEqn]);
			printf(", total number of CVs=%d \n", cvora_file[mpi_size_file]);
			cout<<"  Grid tolerance = "<<gridTol<<endl;

			if(VANGOGH_DEBUG_LEVEL > 1) {
				cout<<endl;
				for(int mpi_file=0; mpi_file<mpi_size_file; ++mpi_file)
					printf("  RANK in the file = %d:  NCV = %d,  X = %e~%e, Y = %e~%e, Z = %e~%e\n", mpi_file, cvora_file[mpi_file+1]-cvora_file[mpi_file],
							xcvMinArray_file[mpi_file*3], xcvMaxArray_file[mpi_file*3],  xcvMinArray_file[mpi_file*3+1], xcvMaxArray_file[mpi_file*3+1],  xcvMinArray_file[mpi_file*3+2], xcvMaxArray_file[mpi_file*3+2]);
				cout<<endl;
			}

			cout<<"  EIGEN-VALUES in the file : Total "<<nev<<" eigenvalues"<<endl;
			for(size_t iev=0; iev<nev; ++iev)
				printf("    [%2d] %g + %gi\n", iev, Evals[iev].real(), Evals[iev].imag());
			cout<<endl;
		}

		assert(ptOnSCurve_file >= 0);
		assert(cvora_file[mpi_size_file] == cvora[mpi_size]);
		for(int mpi_file=0; mpi_file<mpi_size_file; ++mpi_file)
			assert(cvora_file[mpi_file+1] > 0);

		// 2. Read the body
		//       Structure of the body part:
		//            For each mpi,
		//              1. x_cv (double*3*ncv)
		//              2. direct global modes (double*nev*ncv*nVars*2)
		//              3. adjoint global modes (double*nev*ncv*nVars*2)
		int ncv_file;
		double (*x_cv_file)[3] = NULL;
		double **directEvecs_file = NULL;
		double **adjointEvecs_file = NULL;

		for(int mpi_file=0; mpi_file<mpi_size_file; ++mpi_file) {
			// read data dumped by a CPU
			ncv_file = cvora_file[mpi_file+1]-cvora_file[mpi_file];

			assert(x_cv_file==NULL);			x_cv_file = new double [ncv_file][3];
			assert(directEvecs_file==NULL);		directEvecs_file = new double* [nev];
			for(int iev=0; iev<nev; ++iev)
				directEvecs_file[iev] = new double [ncv_file*nVars*2];
			assert(adjointEvecs_file==NULL);	adjointEvecs_file = new double* [nev];
			for(int iev=0; iev<nev; ++iev)
				adjointEvecs_file[iev] = new double [ncv_file*nVars*2];

			for(int icv=0; icv<ncv_file; ++icv) {
				for(int i=0; i<3; ++i)
					infile.read(reinterpret_cast<char*>(&x_cv_file[icv][i]), sizeof(double));
			}
			for(int iev=0; iev<nev; ++iev)
				infile.read(reinterpret_cast<char*>(directEvecs_file[iev]),  sizeof(double)*(ncv_file*nVars*2));
			for(int iev=0; iev<nev; ++iev)
				infile.read(reinterpret_cast<char*>(adjointEvecs_file[iev]), sizeof(double)*(ncv_file*nVars*2));

			// find the matched data
			bool matchedBox = true;

			double xcvMin_file[3] = {xcvMinArray_file[mpi_file*3], xcvMinArray_file[mpi_file*3+1], xcvMinArray_file[mpi_file*3+2]};
			double xcvMax_file[3] = {xcvMaxArray_file[mpi_file*3], xcvMaxArray_file[mpi_file*3+1], xcvMaxArray_file[mpi_file*3+2]};

			if(xcvMin_file[0]>xcvMax[0] || xcvMin_file[1]>xcvMax[1] || xcvMin_file[2]>xcvMax[2])
				matchedBox = false;
			if(xcvMax_file[0]<xcvMin[0] || xcvMax_file[1]<xcvMin[1] || xcvMax_file[2]<xcvMin[2])
				matchedBox = false;

			if(matchedBox) {
				size_t iv = 0;  // Note: We will get only the first vector associated with the first eigen-value

				for(int icv=0; icv<ncv; ++icv) {
					int icv_file;
					int count = findMatchedIndex(icv_file, icv, ncv_file, x_cv_file, gridTol);
					if(count>0) {
						for(size_t ivar=0; ivar<5+nScal; ++ivar) {
							int index = icv_file*nVars*2 + ivar*2;
							DirectEvecs[iv][ivar][icv]  = complex<double>(directEvecs_file[iv][index],  directEvecs_file[iv][index+1]);
							AdjointEvecs[iv][ivar][icv] = complex<double>(adjointEvecs_file[iv][index], adjointEvecs_file[iv][index+1]);
						}

						countFound[icv] += count;
					}
				}
			}

			delete [] x_cv_file;			x_cv_file = NULL;
			for(int iev=0; iev<nev; ++iev)
				delete [] directEvecs_file[iev];
			delete [] directEvecs_file; 	directEvecs_file = NULL;
			for(int iev=0; iev<nev; ++iev)
				delete [] adjointEvecs_file[iev];
			delete [] adjointEvecs_file; 	adjointEvecs_file = NULL;
		}

		// 3. Read the foot
		int dummyInt;
		infile.read(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		if(dummyInt != EOF_ERROR_CHECK_CODE) {
			if(mpi_rank==0)
				cout<<"ERROR! VanGoghWithModels::readEigenPairsRembrandtSerial(): file does not end with EOF_ERROR_CHECK_CODE="<<EOF_ERROR_CHECK_CODE<<endl;
			assert(false);
		}

		// 4. Check if each CV finds a matched CV
		for (int icv=0; icv<ncv; ++icv) {
			if(countFound[icv]==0) {
				printf("  Cannot find a matched point for icv=%d (%.6f,%.6f,%.6f) at mpi_rank=%d \n", icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2],mpi_rank);
				assert(false);
			} else if(countFound[icv]>1) {
				printf("  Too many matched points(%d) for icv=%d (%.6f,%.6f,%.6f) at mpi_rank=%d \n", countFound[icv],icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2],mpi_rank);
				assert(false);
			}
		}
	} else {
		if(mpi_rank==0)
			printf("Cannot read \'%s\' \n", filename);
		assert(false);
	}

	infile.close();

	delete [] countFound;
	delete [] cvora_file;
	delete [] xcvMinArray_file;
	delete [] xcvMaxArray_file;

	MPI_Barrier(mpi_comm);


	double myMinDirectEvec[5+nScal];  	for(int i=0; i<5+nScal; ++i) myMinDirectEvec[i] =  ABSURDLY_BIG_NUMBER;
	double myMaxDirectEvec[5+nScal];  	for(int i=0; i<5+nScal; ++i) myMaxDirectEvec[i] = -ABSURDLY_BIG_NUMBER;
	double myMinAdjointEvec[5+nScal]; 	for(int i=0; i<5+nScal; ++i) myMinAdjointEvec[i] =  ABSURDLY_BIG_NUMBER;
	double myMaxAdjointEvec[5+nScal]; 	for(int i=0; i<5+nScal; ++i) myMaxAdjointEvec[i] = -ABSURDLY_BIG_NUMBER;
	for(int icv=0; icv<ncv; ++icv) {
		for(size_t ivar=0; ivar<5+nScal; ++ivar) {
			double mag = sqrt( pow(DirectEvecs[0][ivar][icv].real(), 2.0) + pow(DirectEvecs[0][ivar][icv].imag(), 2.0) );
			myMinDirectEvec[ivar] = min(myMinDirectEvec[ivar], mag);
			myMaxDirectEvec[ivar] = max(myMaxDirectEvec[ivar], mag);

			mag = sqrt( pow(AdjointEvecs[0][ivar][icv].real(), 2.0) + pow(AdjointEvecs[0][ivar][icv].imag(), 2.0) );
			myMinAdjointEvec[ivar] = min(myMinAdjointEvec[ivar], mag);
			myMaxAdjointEvec[ivar] = max(myMaxAdjointEvec[ivar], mag);
		}
	}

	double minDirectEvec[5+nScal];
	double maxDirectEvec[5+nScal];
	double minAdjointEvec[5+nScal];
	double maxAdjointEvec[5+nScal];
	MPI_Allreduce(myMinDirectEvec,  minDirectEvec,  5+nScal, MPI_DOUBLE, MPI_MIN, mpi_comm);
	MPI_Allreduce(myMaxDirectEvec,  maxDirectEvec,  5+nScal, MPI_DOUBLE, MPI_MAX, mpi_comm);
	MPI_Allreduce(myMinAdjointEvec, minAdjointEvec, 5+nScal, MPI_DOUBLE, MPI_MIN, mpi_comm);
	MPI_Allreduce(myMaxAdjointEvec, maxAdjointEvec, 5+nScal, MPI_DOUBLE, MPI_MAX, mpi_comm);

	if(mpi_rank == 0) {
		cout<<"  RANGE (MAGNITUDE) OF 1ST DIRECT-VECTORS  = ";
		for(int ivar=0; ivar<5+nScal-1; ++ivar)
			printf("%.3e~%.3e | ", minDirectEvec[ivar], maxDirectEvec[ivar]);
		printf("%.3e~%.3e\n", minDirectEvec[5+nScal-1], maxDirectEvec[5+nScal-1]);

		cout<<"  RANGE (MAGNITUDE) OF 1ST ADJOINT-VECTORS = ";
		for(int ivar=0; ivar<5+nScal-1; ++ivar)
			printf("%.3e~%.3e | ", minAdjointEvec[ivar], maxAdjointEvec[ivar]);
		printf("%.3e~%.3e\n", minAdjointEvec[5+nScal-1], maxAdjointEvec[5+nScal-1]);

		cout<<endl;
	}

	MPI_Barrier(mpi_comm);
}

/*
 * Method: readFieldFromQ0bin
 * --------------------------
 * Read field data from Q0_PT000**.bin
 * Original code: readQ0() in joeWithPsALcont.h
 */
//void VanGogh::readFieldFromQ0bin(const int ptOnSCurve, double* rho_unstable, double (*rhou_unstable)[3], double* rhoE_unstable, double* kine_unstable, double* omega_unstable) {
//	// filename
//	stringstream ss;
//	ss<<"Q0_PT"<<setfill('0')<<setw(5)<<ptOnSCurve<<".bin";
//	string filename(ss.str());
//
//	// number of variables
//	int nScal = scalarTranspEqVector.size();
//	int m = nScal+5; // number of variables (rho, rhou, rhov, rhow, rhoE, and scalars)
//
//	// variables related to the Mesh
//	const double gridTol = 1.0e-9;
//	int *countFound = new int [ncv]; 	for(int icv=0; icv<ncv; ++icv) countFound[icv] = 0;
//
//	// read the file
//	ifstream infile(filename.c_str(), ios::in | ios::binary);
//	if(infile.is_open()) {
//		int step0;
//		double lambda0, lambda1, dsTry;
//
//		int mpi_size_file, totalNcv_file;
//		// read some parameters from the file
//		infile.read(reinterpret_cast<char*>(&step0), sizeof(int)); // Note: if you simply use the ">>" operator, it will cause some problems due to bugs in a g++ library on linux
//		infile.read(reinterpret_cast<char*>(&lambda0), sizeof(double));
//		infile.read(reinterpret_cast<char*>(&lambda1), sizeof(double));
//		infile.read(reinterpret_cast<char*>(&dsTry), sizeof(double));
//		infile.read(reinterpret_cast<char*>(&mpi_size_file), sizeof(int));
//		infile.read(reinterpret_cast<char*>(&totalNcv_file), sizeof(int));
//
//		if(mpi_rank==0)
//			printf("VanGoghKOm::readFieldFromQ0bin(): Reading \'%s\' -- pt=%d, lambda=%.8f, arclength=%.8f, total number of CVs=%d \n\n",
//					filename.c_str(), step0, lambda0, dsTry, totalNcv_file);
//		assert(step0 >= 0);
//		assert(totalNcv_file == cvora[mpi_size]);
//
//		// variables
//		int ncv_file;
//		double (*x_cv_file)[3] = NULL;
//		double *q0_file = NULL;
//
//		// read each data
//		for(int mpi_file=0; mpi_file<mpi_size_file; ++mpi_file) {
//			// read data dumped by a CPU
//			infile.read(reinterpret_cast<char*>(&ncv_file), sizeof(int));
//
//			assert(x_cv_file==NULL);	x_cv_file = new double [ncv_file][3];
//			assert(q0_file==NULL);		q0_file = new double [ncv_file*m];
//
//			for(int icv_file=0; icv_file<ncv_file; ++icv_file) {
//				for(int i=0; i<3; ++i)
//					infile.read(reinterpret_cast<char*>(&x_cv_file[icv_file][i]), sizeof(double));
//				for(int i=0; i<m; ++i)
//					infile.read(reinterpret_cast<char*>(&q0_file[m*icv_file+i]), sizeof(double));
//			}
//
//			// find the matached data
//			for(int icv=0; icv<ncv; ++icv) {
//				int icv_file;
//				int count = findMatchedIndex(icv_file, icv, ncv_file, x_cv_file, gridTol);
//				if(count>0) {
//					rho_unstable[icv] = q0_file[m*icv_file];
//					rhou_unstable[icv][0] = q0_file[m*icv_file+1];
//					rhou_unstable[icv][1] = q0_file[m*icv_file+2];
//					rhou_unstable[icv][2] = q0_file[m*icv_file+3];
//					rhoE_unstable[icv] = q0_file[m*icv_file+4];
////					kine_unstable[icv] = q0_file[m*icv_file+5];
////					omega_unstable[icv] = q0_file[m*icv_file+6];
//
//					countFound[icv] += count;
//				}
//			}
//
//			delete [] x_cv_file;	x_cv_file = NULL;
//			delete [] q0_file;		q0_file = NULL;
//		}
//
//		// check if each CV found a matched CV
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
//		updateCvDataG1G2(rho_unstable,  REPLACE_DATA);
//		updateCvDataG1G2(rhou_unstable, REPLACE_ROTATE_DATA);
//		updateCvDataG1G2(rhoE_unstable, REPLACE_DATA);
////		updateCvDataG1G2(kine_unstable, REPLACE_DATA);
////		updateCvDataG1G2(omega_unstable, REPLACE_DATA);
//	}  else {
//		if(mpi_rank==0)
//			printf("VanGoghKOm::readFieldFromQ0bin(): Cannot read \'%s\' \n", filename.c_str());
//		throw(-1);
//	}
//
//	delete [] countFound; 	countFound = NULL;
//
//	infile.close();
//}

/*
 * Method: readPsALCdumpedDataSerial
 * ---------------------------------
 * Read field data from file
 * Original code: IkeWithModels::readPsALCdumpedDataSerial()
 */
void VanGoghWithModels::readPsALCdumpedDataSerial(const char filename[], double* qVec, const int NcontrolEqns) {
MPI_Barrier(mpi_comm); for(int i=0; i<777777; ++i) {}; MPI_Barrier(mpi_comm);
if(mpi_rank==0) cout<<"VanGoghWithModels::readPsALCdumpedDataSerial(): start to read the file"<<endl;
MPI_Barrier(mpi_comm); for(int i=0; i<777777; ++i) {}; MPI_Barrier(mpi_comm);


	// Check some possible errors
	for(int i=0; i<3; ++i)
		if(xcvMax[i] < xcvMin[i]) {
			cerr<<"ERROR in VanGoghWithModels::readPsALCdumpedDataSerial(): mpi_rank == "<<mpi_rank<<" : xcvMax["<<i<<"] < xcvMin["<<i<<"]"<<endl;
			throw(VANGOGH_ERROR_CODE);
		}
	assert(qVec != NULL);

	// number of variables
	int nScal = scalarTranspEqVector.size();
	int m = nScal+5; // number of variables (rho, rhou, rhov, rhow, rhoE, and scalars)

	// variables related to the Mesh
	const double gridTol = 1.0e-9;

	int *countFound = new int [ncv];
	for(int icv=0; icv<ncv; ++icv)
		countFound[icv] = 0;

	// read the file
	ifstream infile(filename, ios::in | ios::binary);
	if(infile.is_open()) {
		int NcontrolEqns_file;
		int step0;
		double *lambda0 = NULL;
		double *lambda1 = NULL;
		double dsTry;

		int mpi_size_file;
		int *cvora_file = NULL;
		double *xcvMinArray_file = NULL;
		double *xcvMaxArray_file = NULL;

		// 1. Read the header part of the file
		//      Structure of the header part:
		//              1. step (int)                      2. NcontrolEqns (int)
		//              3. lambda0 (double*NcontrolEqns)   4. lambda1 (double*NcontrolEqns)
		//              5. ds (double)
		//              6. mpi_size (int)                  7. cvora (double*(mpi_size+1)
		//              8. xMinArray (double*3*mpi_size)   9. xMaxArray (double*3*mpi_size)
		infile.read(reinterpret_cast<char*>(&step0), sizeof(int)); // Note: if you simply use the ">>" operator, it will cause some problem due to bugs in a g++ library on linux
		infile.read(reinterpret_cast<char*>(&NcontrolEqns_file), sizeof(int));
		if(NcontrolEqns_file != NcontrolEqns) {
			if(mpi_rank==0)
				printf("ERROR in VanGoghWithModels::readPsALCdumpedData: NcontrolEqns from file=%d is not equal to given NcontrolEqns=%d\n", NcontrolEqns_file, NcontrolEqns);
			delete [] countFound;
			throw(VANGOGH_ERROR_CODE);
		}
		if(NcontrolEqns_file > 0) {
			lambda0 = new double [NcontrolEqns_file];
			lambda1 = new double [NcontrolEqns_file];
		}
		double dummyDouble;
		for(int iEqn=0; iEqn<NcontrolEqns_file; ++iEqn) {
			infile.read(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
			lambda0[iEqn] = dummyDouble;
		}
		for(int iEqn=0; iEqn<NcontrolEqns_file; ++iEqn) {
			infile.read(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
			lambda1[iEqn] = dummyDouble;
		}
		infile.read(reinterpret_cast<char*>(&dsTry), sizeof(double));
		infile.read(reinterpret_cast<char*>(&mpi_size_file), sizeof(int));
		cvora_file = new int [mpi_size_file+1];
		infile.read(reinterpret_cast<char*>(cvora_file), sizeof(int)*(mpi_size_file+1));
		xcvMinArray_file = new double [mpi_size_file*3];
		infile.read(reinterpret_cast<char*>(xcvMinArray_file), sizeof(double)*(mpi_size_file*3));
		xcvMaxArray_file = new double [mpi_size_file*3];
		infile.read(reinterpret_cast<char*>(xcvMaxArray_file), sizeof(double)*(mpi_size_file*3));

MPI_Barrier(mpi_comm); for(int i=0; i<777777; ++i) {}; MPI_Barrier(mpi_comm);
if(mpi_rank==0) cout<<"reading the header is done!"<<endl;
MPI_Barrier(mpi_comm); for(int i=0; i<777777; ++i) {}; MPI_Barrier(mpi_comm);

		if(mpi_rank==0) {
			printf("\nReading \'%s\': step=%d, lambda1=", filename, step0);
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
				printf("%.8f, ", lambda1[iEqn]);
			printf("ds=%.8f, total number of CVs=%d, mpi_size of the file=%d\n", dsTry, cvora_file[mpi_size_file], mpi_size_file);
			cout<<"  Grid tolerance = "<<gridTol<<endl;
			cout<<"  XCV_BOX = "<<endl;
			for(int impi=0; impi<mpi_size_file; ++impi)
				printf("    %3d: min=(%.3f, %.3f, %.3f)   max=(%.3f, %.3f, %.3f)\n", impi,
						xcvMinArray_file[impi*3],xcvMinArray_file[impi*3+1],xcvMinArray_file[impi*3+2],
						xcvMaxArray_file[impi*3],xcvMaxArray_file[impi*3+1],xcvMaxArray_file[impi*3+2]);
		}

MPI_Barrier(mpi_comm); for(int i=0; i<77777; ++i) {}; MPI_Barrier(mpi_comm);
if(mpi_rank==0) cout<<"cvora_file[mpi_size_file]="<<cvora_file[mpi_size_file]<<endl;
if(mpi_rank==0) cout<<"cvora[mpi_size]="<<cvora[mpi_size]<<endl;
if(mpi_rank==0)
	for(int mpi_file=0; mpi_file<mpi_size_file; ++mpi_file)
		cout<<"cvora_file["<<mpi_file+1<<"]="<<cvora_file[mpi_file+1]<<endl;
MPI_Barrier(mpi_comm); for(int i=0; i<77777; ++i) {}; MPI_Barrier(mpi_comm);

		assert(step0 >= 0);
		assert(cvora_file[mpi_size_file] == cvora[mpi_size]);
		for(int mpi_file=0; mpi_file<mpi_size_file; ++mpi_file)
			assert(cvora_file[mpi_file+1] > 0);

		// 2. Read the body
		//       Structure of the body part:
		//            For each mpi,
		//              1. x_cv (double*3*ncv)
		//              2. qVec (double*m*ncv)
		int ncv_file;
		double (*x_cv_file)[3] = NULL;
		double *qVec_file = NULL;

		for(int mpi_file=0; mpi_file<mpi_size_file; ++mpi_file) {
			// read data dumped by a CPU
			ncv_file = cvora_file[mpi_file+1]-cvora_file[mpi_file];

if(mpi_rank==1) cout<<"reading mpi_file = "<<mpi_file<<endl;

			assert(x_cv_file==NULL);	x_cv_file = new double [ncv_file][3];
			assert(qVec_file==NULL);	qVec_file = new double [ncv_file*m];

			for(int icv=0; icv<ncv_file; ++icv) {
				for(int i=0; i<3; ++i)
					infile.read(reinterpret_cast<char*>(&x_cv_file[icv][i]), sizeof(double));
			}
			infile.read(reinterpret_cast<char*>(qVec_file), sizeof(double)*(ncv_file*m));

if(mpi_rank==1) cout<<"  reading x_cv and qvec are done"<<endl;

			// find the matched data
			bool matchedBox = true;

			double xcvMin_file[3] = {xcvMinArray_file[mpi_file*3], xcvMinArray_file[mpi_file*3+1], xcvMinArray_file[mpi_file*3+2]};
			double xcvMax_file[3] = {xcvMaxArray_file[mpi_file*3], xcvMaxArray_file[mpi_file*3+1], xcvMaxArray_file[mpi_file*3+2]};

			if(xcvMin_file[0]>xcvMax[0] || xcvMin_file[1]>xcvMax[1] || xcvMin_file[2]>xcvMax[2])
				matchedBox = false;
			if(xcvMax_file[0]<xcvMin[0] || xcvMax_file[1]<xcvMin[1] || xcvMax_file[2]<xcvMin[2])
				matchedBox = false;

			if(matchedBox){
if(mpi_rank==1) cout<<"  start to find matched box"<<endl;
				for(int icv=0; icv<ncv; ++icv) {
					int icv_file;
					int count = findMatchedIndex(icv_file, icv, ncv_file, x_cv_file, gridTol);
					if(count>0) {
						for(int i=0; i<m; ++i)
							qVec[m*icv+i] = qVec_file[m*icv_file+i];
						countFound[icv] += count;
					}
				}
			}
if(mpi_rank==1) cout<<"  finding matched box is done"<<endl;
			delete [] x_cv_file;	x_cv_file = NULL;
			delete [] qVec_file;	qVec_file = NULL;
if(mpi_rank==1) cout<<"  done with mpi_file="<<mpi_file<<endl;
		}

MPI_Barrier(mpi_comm); for(int i=0; i<77777; ++i) {}; MPI_Barrier(mpi_comm);
if(mpi_rank==0) cout<<"Reading the body is done!"<<endl;
MPI_Barrier(mpi_comm); for(int i=0; i<77777; ++i) {}; MPI_Barrier(mpi_comm);

		// 3. Read the foot
		int dummyInt;
		infile.read(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		if(dummyInt != EOF_ERROR_CHECK_CODE) {
			if(mpi_rank==0)
				cout<<"ERROR! VanGoghWithModels::readPsALCdumpedDataSerial(): file does not end with EOF_ERROR_CHECK_CODE="<<EOF_ERROR_CHECK_CODE<<endl;
			assert(false);
		}

		// 4. check if each CV found a matched CV
		for (int icv=0; icv<ncv; ++icv) {
			if(countFound[icv]==0) {
				printf("  Cannot find a matched point for icv=%d (%.6f,%.6f,%.6f) at mpi_rank=%d \n", icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2],mpi_rank);
				throw(-1);
			} else if(countFound[icv]>1) {
				printf("  Too many matched points(%d) for icv=%d (%.6f,%.6f,%.6f) at mpi_rank=%d \n", countFound[icv],icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2],mpi_rank);
				throw(-1);
			}
		}

		delete [] lambda0;
		delete [] lambda1;

		delete [] cvora_file;
		delete [] xcvMinArray_file;
		delete [] xcvMaxArray_file;
	}  else {
		if(mpi_rank==0)
			printf("VanGoghWithModels::readPsALCdumpedDataSerial(): Cannot read \'%s\' \n", filename);
		throw(-1);
	}
	infile.close();

	delete [] countFound; 	countFound = NULL;

	MPI_Barrier(mpi_comm);
}

/*
 * Method: updateEigenVecTecplotNS
 * -------------------------------
 * Update NS variabls (rho_Dct_1stReal, rho_Adj_1stReal, etc.) from DirectEvecs and AdjointEvecs.
 * This is just for Tecplot output
 */
void VanGoghWithModels::updateEigenVecTecplotNS() {
	// least-stable direct global modes
	for(int icv=0; icv<ncv; ++icv) {
		rho_Dct_1stReal[icv] = DirectEvecs[0][0][icv].real();
		for(int i=0; i<3; ++i)
			rhou_Dct_1stReal[icv][i] = DirectEvecs[0][1+i][icv].real();
		rhoE_Dct_1stReal[icv] = DirectEvecs[0][4][icv].real();
	}

	// first adjoint global modes
	for(int icv=0; icv<ncv; ++icv) {
		rho_Adj_1stReal[icv] = AdjointEvecs[0][0][icv].real();
		for(int i=0; i<3; ++i)
			rhou_Adj_1stReal[icv][i] = AdjointEvecs[0][1+i][icv].real();
		rhoE_Adj_1stReal[icv] = AdjointEvecs[0][4][icv].real();
	}
}

/*
 * Method: updateUnstableVecNS
 * ---------------------------
 * Update NS variabls (rho_unstable, etc.) from a qVec that has been generated by reading a IKE binary file.
 * This is NOT ONLY for Tecplot output BUT ALSO linear analysis
 */
void VanGoghWithModels::updateUnstableVecNS(double* qVecTemp, const int nVars) {
	assert(qVecTemp!=NULL && rho_unstable!=NULL && rhou_unstable!=NULL && rhoE_unstable!=NULL);

	for(int icv=0; icv<ncv; ++icv) {
		int index = icv*nVars;
		rho_unstable[icv] = qVecTemp[index];
		for(int i=0; i<3; ++i)
			rhou_unstable[icv][i] = qVecTemp[index+1+i];
		rhoE_unstable[icv] = qVecTemp[index+4];
	}
}

/*
 * Method: storeInitRestartNS
 * --------------------------
 * Store initial NS data in arrays (e.g. rho_init, et.c)
 */
void VanGoghWithModels::storeInitRestartNS() {
	if(mpi_rank==0)
		cout<<"VanGoghWithModels::storeInitRestartNS()"<<endl;

	assert(rho_init!=NULL && rhou_init!=NULL && rhoE_init!=NULL);

	for (int icv=0; icv<ncv; ++icv) {
		rho_init[icv] = rho[icv];
		for(int i=0; i<3; ++i)
			rhou_init[icv][i] = rhou[icv][i];
		rhoE_init[icv] = rhoE[icv];
	}
	updateCvDataG1G2(rho_init,  REPLACE_DATA);
	updateCvDataG1G2(rhou_init, REPLACE_ROTATE_DATA);
	updateCvDataG1G2(rhoE_init, REPLACE_DATA);
}

/****************************
 * FUNCTIONS FOR PERTURBATION
 ****************************/
/*
 * Method: getPerturbParams
 * ----------------------
 * Get the parameters for perturbations
 */
void VanGoghWithModels::getPerturbParams(PerturbParams &perturbParams) {
	//  Random perturbation function
	if (!checkParam("RANDOM_PERTURB_DISTRIB")) {
		ParamMap::add("RANDOM_PERTURB_DISTRIB  FUNC=NORMAL  MAGNITUDE=0.05  CLIP=0.001"); // add default values
		if (mpi_rank == 0)
			cout << "WARNING: added keyword \"RANDOM_PERTURB_DISTRIB  FUNC=NORMAL  MAGNITUDE=0.05  CLIP=0.001\" to parameter map!"<<endl;
	}
	string randomDistribFunc = getParam("RANDOM_PERTURB_DISTRIB")->getString("FUNC");
	if(randomDistribFunc.compare("UNIFORM")==0)
		perturbParams.randDistribFunc = UNIFORM_DISTRIB;
	else if(randomDistribFunc.compare("NORMAL")==0)
		perturbParams.randDistribFunc = NORMAL_DISTRIB;
	else {
		if(mpi_rank==0) cout<<"ERROR in VanGoghWithModels::getPerturbParams(): RANDOM_DISTRIB_FUNC = "<<randomDistribFunc<<" is NOT supported"<<endl;
		throw(VANGOGH_ERROR_CODE);
	}
	perturbParams.disturbMag  = getParam("RANDOM_PERTURB_DISTRIB")->getDouble("MAGNITUDE");
	perturbParams.disturbClip = getParam("RANDOM_PERTURB_DISTRIB")->getDouble("CLIP");

	//  Filtering
	if (!checkParam("FILTER")) {
		ParamMap::add("FILTER  APPLY_FILTER=NO  MIN_FILTER_WIDTH=0.01  MAX_FILTER_WIDTH=1.0  FILTER_BOUNDARY_LENGTH=0.0"); // add default values
		if (mpi_rank == 0)
			cout << "WARNING: added keyword \"FILTER  APPLY_FILTER=NO  MIN_FILTER_WIDTH=0.01  MAX_FILTER_WIDTH=1.0  FILTER_BOUNDARY_LENGTH=0.0\" to parameter map!"<<endl;
	}
	string boolStringF = getParam("FILTER")->getString("APPLY_FILTER");
	std::transform(boolStringF.begin(), boolStringF.end(), boolStringF.begin(), ::tolower);
	if(boolStringF.compare("yes")==0 || boolStringF.compare("y")==0 || boolStringF.compare("true")==0 || boolStringF.compare("t")==0) {
		perturbParams.useFiltering = true;

		// get filter size
		perturbParams.minFilterWidth = getParam("FILTER")->getDouble("MIN_FILTER_WIDTH");
		perturbParams.maxFilterWidth = getParam("FILTER")->getDouble("MAX_FILTER_WIDTH");
		perturbParams.filterBoundaryLength = getParam("FILTER")->getDouble("FILTER_BOUNDARY_LENGTH");
	}

	//  Smoothing
	if (!checkParam("SMOOTHING")) {
		ParamMap::add("SMOOTHING  APPLY_SMOOTHING=NO  SMOOTH_XMIN=-1.0e10  SMOOTH_XMAX=1.0e10  SMOOTH_EDGE_SIZE=1.0e-5"); // add default values
		if (mpi_rank == 0)
			cout << "WARNING: added keyword \"SMOOTHING  APPLY_SMOOTHING=NO  SMOOTH_XMIN=-1.0e10  SMOOTH_XMAX=1.0e10  SMOOTH_EDGE_SIZE=1.0e-5\" to parameter map!"<<endl;
	}

	string boolStringSM = getParam("SMOOTHING")->getString("APPLY_SMOOTHING");
	std::transform(boolStringSM.begin(), boolStringSM.end(), boolStringSM.begin(), ::tolower);
	if(boolStringSM.compare("yes")==0 || boolStringSM.compare("y")==0 || boolStringSM.compare("true")==0 || boolStringSM.compare("t")==0) {
		perturbParams.useSmoothing = true;

		// Some parameters for the hat-shaped smoothing function: f(x) = 1.0 - 0.5*(exp(-a*(x-xMin)) + exp(a*(x-xMax))), where a = ln(100)/edgeSize
		// Note: At xMin & xMax, f(x) = 0.0. At the "edges", f(x) = 0.99. At the center of the domain, f(x) ~= 1.0
		//       The smoothing occurs only if useSmoothing is TRUE.
		perturbParams.disturbSmoothXmin = getDoubleParam("DISTURB_SMOOTH_XMIN");
		perturbParams.disturbSmoothXmax = getDoubleParam("DISTURB_SMOOTH_XMAX");
		perturbParams.disturbSmoothEdgeSize = getDoubleParam("DISTURB_SMOOTH_EDGE_SIZE");
	}
}

/*
 * Method: perturbFieldNS
 * ----------------------
 * Perturb the initial field with random variables and apply filter to the perturbations
 */
void VanGoghWithModels::perturbFieldNS() {
	if(mpi_rank==0)
		cout<<classID<<"::perturbFieldNS()"<<endl;

	// ---------------
	// Perturb rho
	// ---------------
	bool applyClipping = true;
	perturbScalar(rho, array_perturb, "rho", applyClipping);

	// ---------------
	// Perturb rhou-X
	// ---------------
	int coord = 0;
	applyClipping = false;
	perturbVector(rhou, coord, array_perturb, "rhou-X", applyClipping);

	// ---------------
	// Perturb rhou-Y
	// ---------------
	coord = 1;
	applyClipping = false;
	perturbVector(rhou, coord, array_perturb, "rhou-Y", applyClipping);

	// ---------------
	// Perturb rhoE
	// ---------------
	applyClipping = true;
	perturbScalar(rhoE, array_perturb, "rhoE", applyClipping);


	if(mpi_rank==0)
		cout<<endl;
	MPI_Barrier(mpi_comm);
}

/*
 * Method: perturbScalar
 * ---------------------
 * Perturb the initial field with random variables and apply filter to the perturbations
 * Note: The array "array_perturb" is used only in this method, but sometimes the user may want to take a look at it.
 *       Thus, it is given as an argument.
 *
 * Return: The coefficient used to construct the perturbation --
 *           coeff = rmsVal * perturbParams.disturbMag ,  where rmsVal is the RMS value of the given scalar array.
 *           perturb[icv] = coeff *array_perturb[icv]
 */
double VanGoghWithModels::perturbScalar(double* scalarArray, double* array_perturb, const char varName[], const bool applyClipping) {
	assert(array_perturb != NULL);

	// Smoothing: Some parameters for the hat-shaped function
	//   f(x) = 1.0 - 0.5*(exp(-a*(x-xMin)) + exp(a*(x-xMax))), where a = ln(100)/edgeSize
	// Note: At xMin & xMax, f(x) = 0.0. At the "edges", f(x) = 0.99. At the center of the domain, f(x) ~= 1.0
	//       The smoothing occurs only if useSmoothing is TRUE.
	double a = log(100.0)/perturbParams.disturbSmoothEdgeSize;
	assert(perturbParams.disturbSmoothXmin < perturbParams.disturbSmoothXmax);
	assert(a>0);

	// -----------------------
	// Generate a random array
	// -----------------------
	if(perturbParams.randDistribFunc == UNIFORM_DISTRIB) {
		randArrayUniform(array_perturb, ncv); // note: randArrayUniform() generates random variables in [0,1]
		for(int icv=0; icv<ncv; ++icv)
			array_perturb[icv] -= 0.5;        // shifting
	} else
		randArrayNormal(array_perturb, ncv);  // Note: default = normal distribution
//	updateCvDataG1G2(array_perturb, REPLACE_DATA);

	// -----------------------
	// Filter the random array
	// -----------------------
	// note: parameters for the filtering such as maxFilterWidth and filterBoundaryLength were already set
	if(perturbParams.useFiltering) {
		double *tempFiltered = new double [ncv_ggff];

		// Check if filter has been constructed.
		if(!filterConsturcted())
			buildCvDifferentialFilter();
		assert(perturbParams.maxFilterWidth > 0.0);

		applyDiffFilterG1G2(tempFiltered, array_perturb);
		for(int icv=0; icv<ncv; ++icv)
			array_perturb[icv] = tempFiltered[icv];
//		updateCvDataG1G2(array_perturb, REPLACE_DATA);

		delete [] tempFiltered;
	}

	// -----------------------------
	// Update the desinated variable
	// -----------------------------
	double rmsVal = calcRMS(scalarArray, false); // get RMS
	double coeff = rmsVal * perturbParams.disturbMag; // calculate the coefficient based on the RMS

	int myCountClip = 0;
	for(int icv=0; icv<ncv; ++icv) {
		double disturb = coeff*array_perturb[icv];
		if(perturbParams.useSmoothing) {
			double x = x_cv[icv][0];
			if(x<perturbParams.disturbSmoothXmin || x>perturbParams.disturbSmoothXmax)
				disturb = 0.0;
			else
				disturb *= 1.0 - 0.5*(exp(-a*(x-perturbParams.disturbSmoothXmin)) + exp(a*(x-perturbParams.disturbSmoothXmax)));
		}
		if(applyClipping && scalarArray[icv]+disturb < perturbParams.disturbClip*scalarArray[icv]) {
			scalarArray[icv] = perturbParams.disturbClip*scalarArray[icv]; // clipping
			myCountClip++;
		} else {
			scalarArray[icv] += disturb;
		}
	}
	updateCvDataG1G2(scalarArray,  REPLACE_DATA);

	int countClip; 	MPI_Allreduce(&myCountClip, &countClip, 1, MPI_INT, MPI_SUM, mpi_comm);
	if(mpi_rank == 0) {
		if(applyClipping) {
			if(perturbParams.randDistribFunc == UNIFORM_DISTRIB)
				printf("    Perturb the %s field with Uniform:  min disturb = %.3e, max disturb = %.3e, # clip = %d \n", varName, -0.5*coeff, 0.5*coeff, countClip);
			else
				printf("    Perturb the %s field with Gaussian:  std = %.3e, # clip = %d \n", varName, coeff, countClip);
		} else {
			if(perturbParams.randDistribFunc == UNIFORM_DISTRIB)
				printf("    Perturb the %s field with Uniform:  min disturb = %.3e, max disturb = %.3e\n", varName, -0.5*coeff, 0.5*coeff);
			else
				printf("    Perturb the %s field with Gaussian:  std = %.3e\n", varName, coeff);
		}
	}

	return coeff;
}

/*
 * Method: perturbVector
 * ---------------------
 * Perturb the initial field with random variables and apply filter to the perturbations
 * Note: The array "array_perturb" is used only in this method, but sometimes the user may want to take a look at it.
 *       Thus, it is given as an argument.
 *
 * Return: The coefficient used to construct the perturbation --
 *           coeff = rmsVal * perturbParams.disturbMag ,  where rmsVal is the RMS value of the given scalar array.
 *           perturb[icv] = coeff *array_perturb[icv]
 */
double VanGoghWithModels::perturbVector(double (*vectorArray)[3], const int coord, double* array_perturb, const char varName[], const bool applyClipping) {
	assert(array_perturb != NULL);
	assert(coord>=0 && coord<3);

	// Smoothing: Some parameters for the hat-shaped function
	//   f(x) = 1.0 - 0.5*(exp(-a*(x-xMin)) + exp(a*(x-xMax))), where a = ln(100)/edgeSize
	// Note: At xMin & xMax, f(x) = 0.0. At the "edges", f(x) = 0.99. At the center of the domain, f(x) ~= 1.0
	//       The smoothing occurs only if useSmoothing is TRUE.
	double a = log(100.0)/perturbParams.disturbSmoothEdgeSize;
	assert(perturbParams.disturbSmoothXmin < perturbParams.disturbSmoothXmax);
	assert(a>0);

	// -----------------------
	// Generate a random array
	// -----------------------
	if(perturbParams.randDistribFunc == UNIFORM_DISTRIB) {
		randArrayUniform(array_perturb, ncv); // note: randArrayUniform() generates random variables in [0,1]
		for(int icv=0; icv<ncv; ++icv)
			array_perturb[icv] -= 0.5;        // shifting
	} else
		randArrayNormal(array_perturb, ncv);  // Note: default = normal distribution
//	updateCvDataG1G2(array_perturb, REPLACE_DATA);

	// -----------------------
	// Filter the random array
	// -----------------------
	// note: parameters for the filtering such as maxFilterWidth and filterBoundaryLength were already set
	if(perturbParams.useFiltering) {
		double *tempFiltered = new double [ncv_ggff];

		// Check if filter has been constructed.
		if(!filterConsturcted())
			buildCvDifferentialFilter();
		assert(perturbParams.maxFilterWidth > 0.0);

		applyDiffFilterG1G2(tempFiltered, array_perturb);
		for(int icv=0; icv<ncv; ++icv)
			array_perturb[icv] = tempFiltered[icv];
//		updateCvDataG1G2(array_perturb, REPLACE_DATA);

		delete [] tempFiltered;
	}

	// -----------------------------
	// Update the desinated variable
	// -----------------------------
	double rmsVal = calcRMS3D(vectorArray, coord, false); // get RMS
	double coeff = rmsVal * perturbParams.disturbMag; // calculate the coefficient based on the RMS

	int myCountClip = 0;
	for(int icv=0; icv<ncv; ++icv) {
		double disturb = coeff*array_perturb[icv];
		if(perturbParams.useSmoothing) {
			double x = x_cv[icv][0];
			if(x<perturbParams.disturbSmoothXmin || x>perturbParams.disturbSmoothXmax)
				disturb = 0.0;
			else
				disturb *= 1.0 - 0.5*(exp(-a*(x-perturbParams.disturbSmoothXmin)) + exp(a*(x-perturbParams.disturbSmoothXmax)));
		}

		if(applyClipping && vectorArray[icv][coord]+disturb < perturbParams.disturbClip*vectorArray[icv][coord]) {
			vectorArray[icv][coord] = perturbParams.disturbClip * vectorArray[icv][coord]; // clipping
			myCountClip++;
		} else {
			vectorArray[icv][coord] += disturb;
		}
	}
	updateCvDataG1G2(vectorArray,  REPLACE_ROTATE_DATA); // Update all the three coordinates!!

	int countClip; 	MPI_Allreduce(&myCountClip, &countClip, 1, MPI_INT, MPI_SUM, mpi_comm);
	if(mpi_rank == 0) {
		if(applyClipping) {
			if(perturbParams.randDistribFunc == UNIFORM_DISTRIB)
				printf("    Perturb the %s field with Uniform:  min disturb = %.3e, max disturb = %.3e, # clip = %d \n", varName, -0.5*coeff, 0.5*coeff, countClip);
			else
				printf("    Perturb the %s field with Gaussian:  std = %.3e, # clip = %d \n", varName, coeff, countClip);
		} else {
			if(perturbParams.randDistribFunc == UNIFORM_DISTRIB)
				printf("    Perturb the %s field with Uniform:  min disturb = %.3e, max disturb = %.3e\n", varName, -0.5*coeff, 0.5*coeff);
			else
				printf("    Perturb the %s field with Gaussian:  std = %.3e\n", varName, coeff);
		}
	}

	return coeff;
}

/****************************
 * OPTIMAL METRIC
 ****************************/

///*
// * Method: calcMetricsAdj
// * ----------------------
// * Calculate metrics using adjoint: Let phi=flow field, phi0=flow field on the unstable branch, y=least stable eigenvector, a=adjoint vector,
// *                                           (phi-phi0)*a
// *                                  metric = ------------
// *                                                y*a
// * Return: metricNumer = (phi-phi0)*a
// *         metricDenom = y*a
// */
//void VanGoghKOm::calcMetricsAdj(double* metricNumer, double* metricDenom,
//		const double *rho_unstable, const double (*rhou_unstable)[3], const double *rhoE_unstable, const double *kine_unstable, const double *omega_unstable,
//		const double *rho_1stMode, const double (*rhou_1stMode)[3], const double *rhoE_1stMode, const double *kine_1stMode, const double *omega_1stMode,
//		const double *rho_adj, const double (*rhou_adj)[3], const double *rhoE_adj, const double *kine_adj, const double *omega_adj) {
//	// define variables
//	int m = 7; // number of variables
//	double myTemp[m];
//
//	// calculate numerator = (phi-phi0)*a
//	for(int i=0; i<m; ++i)
//		myTemp[i] = 0.0;
//	for(int icv=0; icv<ncv; ++icv) {
//		myTemp[0] += (rho[icv]-rho_unstable[icv])*rho_adj[icv];
//		myTemp[1] += (rhou[icv][0]-rhou_unstable[icv][0])*rhou_adj[icv][0];
//		myTemp[2] += (rhou[icv][1]-rhou_unstable[icv][1])*rhou_adj[icv][1];
//		myTemp[3] += (rhou[icv][2]-rhou_unstable[icv][2])*rhou_adj[icv][2];
//		myTemp[4] += (rhoE[icv]-rhoE_unstable[icv])*rhoE_adj[icv];
//		myTemp[5] += (kine[icv]-kine_unstable[icv])*kine_adj[icv];
//		myTemp[6] += (omega[icv]-omega_unstable[icv])*omega_adj[icv];
//	}
//	MPI_Allreduce(myTemp, metricNumer, m, MPI_DOUBLE, MPI_SUM, mpi_comm);
//
//	// calculate denominator = y*a
//	for(int i=0; i<m; ++i)
//		myTemp[i] = 0.0;
//	for(int icv=0; icv<ncv; ++icv) {
//		myTemp[0] += rho_1stMode[icv]*rho_adj[icv];
//		myTemp[1] += rhou_1stMode[icv][0]*rhou_adj[icv][0];
//		myTemp[2] += rhou_1stMode[icv][1]*rhou_adj[icv][1];
//		myTemp[3] += rhou_1stMode[icv][2]*rhou_adj[icv][2];
//		myTemp[4] += rhoE_1stMode[icv]*rhoE_adj[icv];
//		myTemp[5] += kine_1stMode[icv]*kine_adj[icv];
//		myTemp[6] += omega_1stMode[icv]*omega_adj[icv];
//	}
//	MPI_Allreduce(myTemp, metricDenom, m, MPI_DOUBLE, MPI_SUM, mpi_comm);
//}
