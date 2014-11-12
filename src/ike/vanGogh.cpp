#include "vanGogh.h"

void calcPTM(double &p, double &T, double &M, const double gamma, const double Rgas, const double rho, const double rhou, const double rhoE) {
	p = (gamma-1) * (rhoE - 0.5*pow(rhou, 2.0)/rho);
	T = p/rho/Rgas;
	M = rhou/rho / sqrt(fabs(gamma*p/rho));
}

// ###########################################################################################
// ------                                                                               ------
// ------                       Perturbation tests called VanGogh                       ------
// ------                                                                               ------
// ###########################################################################################
/**********************
 * OVERLOADED METHODS
 **********************/
/*
 * Method: run
 * -----------
 * Original code = run() in JoeWithModels.h
 */
void VanGogh::run() {
	if(mpi_rank == 0)
		cout<<endl
			<<"VanGogh::run() "<<endl;

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

	for(int i=0; i<3; ++i) {
		xcvMin[i] =  1.0e22;
		xcvMax[i] = -1.0e22;
	}
	for(int icv=0; icv<ncv; ++icv) {
		for(int i=0; i<3; ++i) {
			xcvMin[i] = min(xcvMin[i], x_cv[icv][i]);
			xcvMax[i] = max(xcvMax[i], x_cv[icv][i]);
		}
	}

	calcStateVariables();

//	// Read eigenmodes and adjoint modes
	int ptOnSCurve = getIntParam("PT_ON_S_CURVE", "-1");
//	readLinSystem(ptOnSCurve, rho_1stMode, rhou_1stMode, rhoE_1stMode, kine_1stMode, omega_1stMode,
//			rho_adj, rhou_adj, rhoE_adj, kine_adj, omega_adj);

//	// Read the point on the unstable branch (Note: the user should set ptOnSCurve to be on the unstable branch)
//	char filename[20];
//	sprintf(filename, "Q1_PT%05d.bin", ptOnSCurve);
//	readPsALCdumpedDataSerial(filename, rho_unstable, rhou_unstable, rhoE_unstable);

	// Save the initial data from the restart file if the user wants to reinitialize the flow field by random perturbation
	string reinitField = getStringParam("REINIT_FIELD", "YES");
	if(reinitField.compare("YES")==0 || reinitField.compare("yes")==0 || reinitField.compare("Y")==0 || reinitField.compare("y")==0) {
		storeInitRestart();
		itest = 0;
		ntests = getIntParam("NTESTS", "1");
	} else {
		itest = getIntParam("CASE_NUM", "0");
		ntests = itest+1;
	}

	// run simulations
	char restartFilename[32];
	while(itest < ntests) {
		if(mpi_rank==0)
			cout<<endl
				<<">> RUNNING "<<itest<<"TH SIMULATION"<<endl
				<<endl;

		step = 0;

//		// initialize models
//		initialHookScalarRansTurbModel();
//		initialHookScalarRansCombModel();
//		for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
//			updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);

		// initialize Navier-Stokes equations
		initialHook();
		updateCvDataG1G2(rho, REPLACE_DATA);
		updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
		updateCvDataG1G2(rhoE, REPLACE_DATA);

		// write restart file
		sprintf(restartFilename, "Restart.case%04d.init", itest);
		writeRestart(restartFilename);

		// run solver
		string tIntName = getStringParam("TIME_INTEGRATION");

		if (tIntName == "FORWARD_EULER")                 runForwardEuler();
		else if (tIntName == "RK")                       runRK();
		else if (tIntName == "BACKWARD_EULER")           runBackwardEuler();
		else if (tIntName == "BACKWARD_EULER_COUPLED")   runBackwardEulerCoupled();
		else if (tIntName == "BDF2")                     runBDF2();
		else
			if (mpi_rank == 0) {
				cerr << "ERROR: wrong time integration scheme specified !" << endl;
				cerr << "available integration schemes are: FORWARD_EULER, RK, BACKWARD_EULER, BDF2" << endl;
			}

		sprintf(restartFilename, "Restart.case%04d.%06d.out", itest, step);
		writeRestart(restartFilename);

		++itest;
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
double VanGogh::calcRMS(const double *vec, bool volumeAvg) {
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
double VanGogh::calcRMS3D(const double (*vec)[3], const int dim, bool volumeAvg) {
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
double VanGogh::findMinValue(double* array, const int size) {
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
double VanGogh::findMaxValue(double* array, const int size) {
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
void VanGogh::calcWallDistanceFace(double *wd) {
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
int VanGogh::findMatchedIndex(int& foundIndex, const int icv, const int nCVs, const double (*x_cv_eigen)[3], const double epsil) {
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
 * Method: readEigenBinaryFileHeader()
 * -----------------------------------
 * Read the header part of a given file.
 */
void VanGogh::readEigenBinaryFileHeader(const string &filename, int &nCVs, int &dimEigen, int &nEigens) {
	ifstream file(filename.c_str(), ios::in | ios::binary);
	if(file.is_open()) {
		// read nCVs and nEigens from the file
		file.read(reinterpret_cast<char*>(&nCVs), sizeof(int)); // if you simply use the ">>" operator, it will cause some problems
		file.read(reinterpret_cast<char*>(&dimEigen), sizeof(int)); // if you simply use the ">>" operator, it will cause some problems
		file.read(reinterpret_cast<char*>(&nEigens), sizeof(int)); // if you simply use the ">>" operator, it will cause some problems

		if(mpi_rank==0)
			printf("VanGoghKOm::readEigenBinaryFileHeader(): \'%s\' -- nCVs=%d, dimEigen=%d, nEigens=%d \n", filename.c_str(), nCVs, dimEigen, nEigens);
	} else {
		if(mpi_rank==0)
			printf("Cannot read \'%s\' \n", filename.c_str());
	}
}

/*
 * Method: readEigenBinaryFileBody()
 * ------------------------
 * Read the body part of a given file.
 */
void VanGogh::readEigenBinaryFileBody(const string &filename, double (*x_cv_eigen)[3], complex<double> (*EVecDirect)[NEIGENS], complex<double> (*EVecAdjoint)[NEIGENS]) {
	assert(x_cv_eigen != NULL && EVecDirect != NULL && EVecAdjoint != NULL);

	ifstream file(filename.c_str(), ios::in | ios::binary);
	if(file.is_open()) {
		int nCVs, dimEigen, nEigens;
		// read nCVs and nEigens from the file
		file.read(reinterpret_cast<char*>(&nCVs), sizeof(int)); // if you simply use the ">>" operator, it will cause some problems
		file.read(reinterpret_cast<char*>(&dimEigen), sizeof(int)); // if you simply use the ">>" operator, it will cause some problems
		file.read(reinterpret_cast<char*>(&nEigens), sizeof(int)); // if you simply use the ">>" operator, it will cause some problems

		if(mpi_rank==0)
			printf("VanGoghKOm::readEigenBinaryFileBody(): \'%s\' \n", filename.c_str());
		assert(nEigens == NEIGENS);
		assert(dimEigen == NVARS_EIGEN*nCVs);

		// read data
		// 1. x,y,z coordinates
		for (int i=0; i<nCVs; ++i)
			file.read(reinterpret_cast<char*>(&x_cv_eigen[i][0]), sizeof(double)); // if you simply use the ">>" operator, it will cause some problems
		for (int i=0; i<nCVs; ++i)
			file.read(reinterpret_cast<char*>(&x_cv_eigen[i][1]), sizeof(double)); // if you simply use the ">>" operator, it will cause some problems
		for (int i=0; i<nCVs; ++i)
			file.read(reinterpret_cast<char*>(&x_cv_eigen[i][2]), sizeof(double)); // if you simply use the ">>" operator, it will cause some problems
		// 2. eigenvectors of the direct Jacobian matrix
		double dummyReal, dummyImag;
		for (int i=0; i<dimEigen; ++i) {
			for (int col=0; col<NEIGENS; ++col) {
				file.read(reinterpret_cast<char*>(&dummyReal), sizeof(double)); // if you simply use the ">>" operator, it will cause some problems
				file.read(reinterpret_cast<char*>(&dummyImag), sizeof(double)); // if you simply use the ">>" operator, it will cause some problems
				EVecDirect[i][col] = complex<double>(dummyReal, dummyImag);
			}
		}
		// 3. eigenvectors of the adjoint Jacobian matrix
		for (int i=0; i<dimEigen; ++i) {
			for (int col=0; col<NEIGENS; ++col) {
				file.read(reinterpret_cast<char*>(&dummyReal), sizeof(double)); // if you simply use the ">>" operator, it will cause some problems
				file.read(reinterpret_cast<char*>(&dummyImag), sizeof(double)); // if you simply use the ">>" operator, it will cause some problems
				EVecAdjoint[i][col] = complex<double>(dummyReal, dummyImag);
			}
		}

	} else {
		if(mpi_rank==0)
			printf("Cannot read \'%s\' \n", filename.c_str());
	}
}

/*
 * Method: readLinSystem
 * ---------------------
 * Read eigenmodes and adjoint modes from "dumpedEigen_pt*****.bin"
 */
void VanGogh::readLinSystem(const int ptOnSCurve,
		double *rho_1stMode, double (*rhou_1stMode)[3], double *rhoE_1stMode, double *kine_1stMode, double *omega_1stMode,
		double *rho_adj, double (*rhou_adj)[3], double *rhoE_adj, double *kine_adj, double *omega_adj) {
	if(mpi_rank == 0)
		cout<<"VanGoghKOm::readLinSystem() "<<endl;

	// ------------------------------
	// read the data file
	// ------------------------------
	// define variables
	int nCVs, dimEigen, nEigens;
	double (*x_cv_eigen)[3]                 = NULL;
	complex<double> (*EVecDirect)[NEIGENS]  = NULL;  // note: 'std::complex' class
	complex<double> (*EVecAdjoint)[NEIGENS] = NULL; // note: 'std::complex' class

	// read the header of the binary file
	stringstream ss;
	ss<<"dumpedEigen_pt"<<ptOnSCurve<<".bin";
	readEigenBinaryFileHeader(ss.str(), nCVs, dimEigen, nEigens);

	// allocate memory
	x_cv_eigen = new double[nCVs][3];
	EVecDirect = new complex<double>[dimEigen][NEIGENS];
	EVecAdjoint = new complex<double>[dimEigen][NEIGENS];

	// read the body part of the binary file
	readEigenBinaryFileBody(ss.str(), x_cv_eigen, EVecDirect, EVecAdjoint);

	// ------------------------------
	// store the data on given arrays
	// ------------------------------
//	double *kine = getR1("kine");
//	double *omega = getR1("omega");

	double epsil = 1.0e-8;
	for (int icv=0; icv<ncv; ++icv) {
		int foundIndex;
		int count=findMatchedIndex(foundIndex, icv, nCVs, x_cv_eigen, epsil);

		if(count==0) {
			printf("  Cannot find a matched point for icv=%d (%.6f,%.6f,%.6f) at mpi_rank=%d \n", icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2],mpi_rank);
			printf("    Current epsil=%.5e... Increase epsil \n",epsil);
			throw(-1);
		} else if(count>1) {
			printf("  Too many matched points(%d) for icv=%d (%.6f,%.6f,%.6f) at mpi_rank=%d \n", count,icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2],mpi_rank);
			printf("    Current epsil=%.5e... Decrease epsil \n",epsil);
			throw(-1);
		}

		rho_1stMode[icv] =     EVecDirect[foundIndex*NVARS_EIGEN][0].real();
		rhou_1stMode[icv][0] = EVecDirect[foundIndex*NVARS_EIGEN+1][0].real();
		rhou_1stMode[icv][1] = EVecDirect[foundIndex*NVARS_EIGEN+2][0].real();
		rhou_1stMode[icv][2] = EVecDirect[foundIndex*NVARS_EIGEN+3][0].real();
		rhoE_1stMode[icv] =    EVecDirect[foundIndex*NVARS_EIGEN+4][0].real();
//		kine_1stMode[icv] =    EVecDirect[foundIndex*NVARS_EIGEN+5][0].real();
//		omega_1stMode[icv] =   EVecDirect[foundIndex*NVARS_EIGEN+6][0].real();

		rho_adj[icv] =     EVecAdjoint[foundIndex*NVARS_EIGEN][0].real();
		rhou_adj[icv][0] = EVecAdjoint[foundIndex*NVARS_EIGEN+1][0].real();
		rhou_adj[icv][1] = EVecAdjoint[foundIndex*NVARS_EIGEN+2][0].real();
		rhou_adj[icv][2] = EVecAdjoint[foundIndex*NVARS_EIGEN+3][0].real();
		rhoE_adj[icv] =    EVecAdjoint[foundIndex*NVARS_EIGEN+4][0].real();
//		kine_adj[icv] =    EVecAdjoint[foundIndex*NVARS_EIGEN+5][0].real();
//		omega_adj[icv] =   EVecAdjoint[foundIndex*NVARS_EIGEN+6][0].real();
	}
	updateCvDataG1G2(rho_1stMode,  REPLACE_DATA);
	updateCvDataG1G2(rhou_1stMode, REPLACE_ROTATE_DATA);
	updateCvDataG1G2(rhoE_1stMode, REPLACE_DATA);
//	updateCvDataG1G2(kine_1stMode, REPLACE_DATA);
//	updateCvDataG1G2(omega_1stMode, REPLACE_DATA);

	updateCvDataG1G2(rho_adj,  REPLACE_DATA);
	updateCvDataG1G2(rhou_adj, REPLACE_ROTATE_DATA);
	updateCvDataG1G2(rhoE_adj, REPLACE_DATA);
//	updateCvDataG1G2(kine_adj, REPLACE_DATA);
//	updateCvDataG1G2(omega_adj, REPLACE_DATA);

	// ------------------------------
	// Release memory
	// ------------------------------
	delete [] x_cv_eigen;		x_cv_eigen = NULL;
	delete [] EVecDirect;		EVecDirect = NULL;
	delete [] EVecAdjoint;		EVecAdjoint = NULL;
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
void VanGogh::readPsALCdumpedDataSerial(const char filename[], double* rho_unstable, double (*rhou_unstable)[3], double* rhoE_unstable) {
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
		double *lambda0 = new double [NcontrolEqns_file];
		double *lambda1 = new double [NcontrolEqns_file];
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
		double dummyDouble;
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
			infile.read(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
			lambda0[iEqn] = dummyDouble;
		}
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
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

		if(mpi_rank==0) {
			printf("\nReading \'%s\': step=%d, lambda1=", filename, step0);
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
				printf("%.8f, ", lambda0[iEqn]);
			printf("ds=%.8f, total number of CVs=%d \n", dsTry, cvora_file[mpi_size_file]);
			cout<<"  Grid tolerance = "<<gridTol<<endl;
			cout<<"  XCV_BOX = "<<endl;
			for(int impi=0; impi<mpi_size_file; ++impi)
				printf("    %3d: min=(%.3f, %.3f, %.3f)   max=(%.3f, %.3f, %.3f)\n", impi,
						xcvMinArray_file[impi*3],xcvMinArray_file[impi*3+1],xcvMinArray_file[impi*3+2],
						xcvMaxArray_file[impi*3],xcvMaxArray_file[impi*3+1],xcvMaxArray_file[impi*3+2]);
		}

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

			assert(x_cv_file==NULL);	x_cv_file = new double [ncv_file][3];
			assert(qVec_file==NULL);	qVec_file = new double [ncv_file*m];

			for(int icv=0; icv<ncv_file; ++icv) {
				for(int i=0; i<3; ++i)
					infile.read(reinterpret_cast<char*>(&x_cv_file[icv][i]), sizeof(double));
			}
			infile.read(reinterpret_cast<char*>(qVec_file), sizeof(double)*(ncv_file*m));

			// find the matched data
			bool matchedBox = true;

			double xcvMin_file[3] = {xcvMinArray_file[mpi_file*3], xcvMinArray_file[mpi_file*3+1], xcvMinArray_file[mpi_file*3+2]};
			double xcvMax_file[3] = {xcvMaxArray_file[mpi_file*3], xcvMaxArray_file[mpi_file*3+1], xcvMaxArray_file[mpi_file*3+2]};

			if(xcvMin_file[0]>xcvMax[0] || xcvMin_file[1]>xcvMax[1] || xcvMin_file[2]>xcvMax[2])
				matchedBox = false;
			if(xcvMax_file[0]<xcvMin[0] || xcvMax_file[1]<xcvMin[1] || xcvMax_file[2]<xcvMin[2])
				matchedBox = false;

			if(matchedBox){
				for(int icv=0; icv<ncv; ++icv) {
					int icv_file;
					int count = findMatchedIndex(icv_file, icv, ncv_file, x_cv_file, gridTol);
					if(count>0) {
						rho_unstable[icv]     = qVec_file[m*icv_file];
						rhou_unstable[icv][0] = qVec_file[m*icv_file+1];
						rhou_unstable[icv][1] = qVec_file[m*icv_file+2];
						rhou_unstable[icv][2] = qVec_file[m*icv_file+3];
						rhoE_unstable[icv]    = qVec_file[m*icv_file+4];
//					kine_unstable[icv] = qVec_file[m*icv_file+5];
//					omega_unstable[icv] = qVec_file[m*icv_file+6];

						countFound[icv] += count;
					}
				}
			}
			delete [] x_cv_file;	x_cv_file = NULL;
			delete [] qVec_file;	qVec_file = NULL;
		}

		// 3. Read the foot
		int dummyInt;
		infile.read(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		if(dummyInt != EOF_ERROR_CHECK_CODE) {
			if(mpi_rank==0)
				cout<<"ERROR! VanGogh::readPsALCdumpedDataSerial(): file does not end with EOF_ERROR_CHECK_CODE="<<EOF_ERROR_CHECK_CODE<<endl;
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

		updateCvDataG1G2(rho_unstable,  REPLACE_DATA);
		updateCvDataG1G2(rhou_unstable, REPLACE_ROTATE_DATA);
		updateCvDataG1G2(rhoE_unstable, REPLACE_DATA);
//		updateCvDataG1G2(kine_unstable, REPLACE_DATA);
//		updateCvDataG1G2(omega_unstable, REPLACE_DATA);

		delete [] lambda0;
		delete [] lambda1;

		delete [] cvora_file;
		delete [] xcvMinArray_file;
		delete [] xcvMaxArray_file;
	}  else {
		if(mpi_rank==0)
			printf("VanGogh::readPsALCdumpedDataSerial(): Cannot read \'%s\' \n", filename);
		throw(-1);
	}
	infile.close();

	delete [] countFound; 	countFound = NULL;

	MPI_Barrier(mpi_comm);
}

/*
 * Method: storeInitRestart
 * ------------------------
 * Store initial data in arrays and calculate the norms of them for future uses
 */
void VanGogh::storeInitRestart() {
	if(mpi_rank==0)
		cout<<"VanGoghKOm::storeInitRestart()"<<endl;

//	assert(rho_init!=NULL && rhou_init!=NULL && rhoE_init!=NULL && kine_init!=NULL && omega_init!=NULL);
	assert(rho_init!=NULL && rhou_init!=NULL && rhoE_init!=NULL);

//	double *kine = getR1("kine");
//	double *omega = getR1("omega");

	if(mpi_rank == 0)
		cout<<"VanGoghKOm::storeInitRestart() ";

	for (int icv=0; icv<ncv; ++icv) {
		rho_init[icv] = rho[icv];
		for(int i=0; i<3; ++i)
			rhou_init[icv][i] = rhou[icv][i];
		rhoE_init[icv] = rhoE[icv];

//		kine_init[icv] = kine[icv];
//		omega_init[icv] = omega[icv];
	}
	updateCvDataG1G2(rho_init,  REPLACE_DATA);
	updateCvDataG1G2(rhou_init, REPLACE_ROTATE_DATA);
	updateCvDataG1G2(rhoE_init, REPLACE_DATA);
//	updateCvDataG1G2(kine_init, REPLACE_DATA);
//	updateCvDataG1G2(omega_init, REPLACE_DATA);
}

/*
 * Method: reinitializeNS
 * ----------------------
 *
 */
void VanGogh::reinitializeNS(const double* rho_init, const double (*rhou_init)[3], const double* rhoE_init) {
	if(mpi_rank==0)
		cout<<"VanGoghKOm::reinitializeNS()"<<endl;

	assert(rho_init!=NULL && rhou_init!=NULL && rhoE_init!=NULL);

	for (int icv=0; icv<ncv; ++icv) {
		rho[icv] = rho_init[icv];
		for(int i=0; i<3; ++i)
			rhou[icv][i] = rhou_init[icv][i];
		rhoE[icv] = rhoE_init[icv];
	}
	updateCvDataG1G2(rho,  REPLACE_DATA);
	updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
	updateCvDataG1G2(rhoE, REPLACE_DATA);

	MPI_Barrier(mpi_comm);
}

///*
// * Method: reinitializeScalars
// * ---------------------------
// *
// */
//void VanGoghKOm::reinitializeScalars(const double* kine_init, const double* omega_init) {
//	if(mpi_rank==0)
//		cout<<"VanGoghKOm::reinitializeScalars()"<<endl;
//
//	assert(kine_init!=NULL && omega_init!=NULL);
//
//	double *kine = getR1("kine");
//	double *omega = getR1("omega");
//
//	for (int icv=0; icv<ncv; ++icv) {
//		kine[icv] = kine_init[icv];
//		omega[icv] = omega_init[icv];
//	}
//	updateCvDataG1G2(kine, REPLACE_DATA);
//	updateCvDataG1G2(omega, REPLACE_DATA);
//}

///*
// * Method: perturbFieldScalars
// * ---------------------------
// * Perturb the initial field with random variables and apply filter to the perturbations
// */
//void VanGoghKOm::perturbFieldScalars(const RAND_DISTRIB_FUNC randDistribFunc, const bool useSmoothing, const bool useFiltering) {
//	if(mpi_rank==0)
//		cout<<"VanGoghKOm::perturbFieldScalars()"<<endl;
//
//	double *kine = getR1("kine");
//	double *omega = getR1("omega");
//
//	// Random variable array and filter array
//	assert(array_perturb == NULL);
//	array_perturb = new double [ncv_ggff];
//
//	double *tempFiltered = NULL;
//	if(useFiltering)
//		tempFiltered = new double [ncv_ggff];
//
//	// Some parameters for the disturbances
//	double disturbMag = getDoubleParam("DISTURB_MAGNITUDE", "0.05");
//	double disturbClip = getDoubleParam("DISTURB_CLIP", "0.001");
//
//	// Some parameters for the hat-shaped smoothing function: f(x) = 1.0 - 0.5*(exp(-a*(x-xMin)) + exp(a*(x-xMax))), where a = ln(100)/edgeSize
//	// Note: At xMin & xMax, f(x) = 0.0. At the "edges", f(x) = 0.99. At the center of the domain, f(x) ~= 1.0
//	//       The smoothing occurs only if useSmoothing is TRUE.
//	double disturbSmoothXmin = getDoubleParam("DISTURB_SMOOTH_XMIN", "-1.0e10");
//	double disturbSmoothXmax = getDoubleParam("DISTURB_SMOOTH_XMAX", "1.0e10");
//	double disturbSmoothEdgeSize = getDoubleParam("DISTURB_SMOOTH_EDGE_SIZE", "1.0e-5");
//	double a = log(100.0)/disturbSmoothEdgeSize;
//	assert(disturbSmoothXmin<disturbSmoothXmax);
//	assert(a>0);
//
//	// -------------------
//	// Generate kine_perturb
//	// -------------------
//	if(randDistribFunc == UNIFORM_DISTRIB) {
//		randArrayUniform(array_perturb, ncv); // note: randArrayUniform() generates random variables in [0,1]
//		for(int icv=0; icv<ncv; ++icv)
//			array_perturb[icv] -= 0.5;        // shifting
//	} else
//		randArrayNormal(array_perturb, ncv);  // default = normal distribution
//
//	updateCvDataG1G2(array_perturb, REPLACE_DATA);
//
//	// -------------------
//	// Filter kine_perturb
//	// -------------------
//	// note: parameters for the filtering such as maxFilterWidth and filterBoundaryLength were already set
//	if(useFiltering) {
//		// Check if filter has been constructed.
//		if(!filterConsturcted())
//			buildCvDifferentialFilter();
//		assert(maxFilterWidth > 0.0);
//
//		applyDiffFilterG1G2(tempFiltered, array_perturb);
//		for(int icv=0; icv<ncv; ++icv)
//			array_perturb[icv] = tempFiltered[icv];
//		updateCvDataG1G2(array_perturb, REPLACE_DATA);
//	}
//
//	// -------------------
//	// Update kine
//	// -------------------
//	double rmsVal = calcRMS(kine, false); // get RMS
//	double coeff = rmsVal*disturbMag; // calculate the coefficient based on the RMS
//	int myCountClip = 0;
//	for(int icv=0; icv<ncv; ++icv) {
//		double disturb = coeff*array_perturb[icv];
//		if(useSmoothing) {
//			double x = x_cv[icv][0];
//			if(x<disturbSmoothXmin || x>disturbSmoothXmax)
//				disturb = 0.0;
//			else
//				disturb *= 1.0 - 0.5*(exp(-a*(x-disturbSmoothXmin)) + exp(a*(x-disturbSmoothXmax)));
//		}
//		if(kine[icv]+disturb < disturbClip*kine[icv]) {
//			kine[icv] = disturbClip*kine[icv];
//			++myCountClip;
//		} else
//			kine[icv] += disturb;
//	}
//	updateCvDataG1G2(kine,  REPLACE_DATA);
//	int countClip; 	MPI_Allreduce(&myCountClip, &countClip, 1, MPI_INT, MPI_SUM, mpi_comm);
//	if(mpi_rank == 0) {
//		if(randDistribFunc == UNIFORM_DISTRIB)
//			printf("    Perturb the kine field with Uniform:  min disturb = %.3e, max disturb = %.3e, # clip = %d \n", -0.5*coeff, 0.5*coeff, countClip);
//		else
//			printf("    Perturb the kine field with Gaussian:  std = %.3e, # clip = %d \n", coeff, countClip);
//	}
//
//	// -------------------
//	// Generate omega_perturb
//	// -------------------
//	if(randDistribFunc == UNIFORM_DISTRIB) {
//		randArrayUniform(array_perturb, ncv); // note: randArrayUniform() generates random variables in [0,1]
//		for(int icv=0; icv<ncv; ++icv)
//			array_perturb[icv] -= 0.5;        // shifting
//	} else
//		randArrayNormal(array_perturb, ncv);  // default = normal distribution
//	updateCvDataG1G2(array_perturb, REPLACE_DATA);
//
//	// -------------------
//	// Filter omega_perturb
//	// -------------------
//	if(useFiltering) {
//		applyDiffFilterG1G2(tempFiltered, array_perturb);
//		for(int icv=0; icv<ncv; ++icv)
//			array_perturb[icv] = tempFiltered[icv];
//		updateCvDataG1G2(array_perturb, REPLACE_DATA);
//	}
//
//	// -------------------
//	// Update omega
//	// -------------------
//	rmsVal = calcRMS(omega, false); // get RMS
//	coeff = rmsVal*disturbMag; // calculate the coefficient based on the RMS
//
//	myCountClip = 0;
//	for(int icv=0; icv<ncv; ++icv) {
//		double disturb = coeff*array_perturb[icv];
//		if(useSmoothing) {
//			double x = x_cv[icv][0];
//			if(x<disturbSmoothXmin || x>disturbSmoothXmax)
//				disturb = 0.0;
//			else
//				disturb *= 1.0 - 0.5*(exp(-a*(x-disturbSmoothXmin)) + exp(a*(x-disturbSmoothXmax)));
//		}
//		if(omega[icv]+disturb < disturbClip*omega[icv]) {
//			omega[icv] = disturbClip*omega[icv];
//			myCountClip++;
//		} else
//			omega[icv] += disturb;
//	}
//	updateCvDataG1G2(omega,  REPLACE_DATA);
//	MPI_Allreduce(&myCountClip, &countClip, 1, MPI_INT, MPI_SUM, mpi_comm);
//	if(mpi_rank == 0) {
//		if(randDistribFunc == UNIFORM_DISTRIB)
//			printf("    Perturb the omega field with Uniform: min disturb = %.3e, max disturb = %.3e, # clip = %d \n", -0.5*coeff, 0.5*coeff, countClip);
//		else
//			printf("    Perturb the omega field with Gaussian: std = %.3e, # clip = %d \n", coeff, countClip);
//	}
//
//	// -------------
//	// Clear memory
//	// -------------
//	delete [] array_perturb; array_perturb = NULL;
//	if(useFiltering)
//		delete [] tempFiltered;
//
//	if(mpi_rank==0)
//		cout<<endl;
//	MPI_Barrier(mpi_comm);
//}

/*
 * Method: perturbFieldNS
 * ----------------------
 * Perturb the initial field with random variables and apply filter to the perturbations
 */
void VanGogh::perturbFieldNS(const RAND_DISTRIB_FUNC randDistribFunc, const bool useSmoothing, const bool useFiltering) {
	if(mpi_rank==0)
		cout<<"VanGoghKOm::perturbFieldNS()"<<endl;

	// Random variable array and filter array
	assert(array_perturb == NULL);
	array_perturb = new double [ncv_ggff];

	double *tempFiltered = NULL;
	if(useFiltering)
		tempFiltered = new double [ncv_ggff];

	// Some parameters for the disturbances
	double disturbMag = getDoubleParam("DISTURB_MAGNITUDE", "0.05");
	double disturbClip = getDoubleParam("DISTURB_CLIP", "0.001");

	// Some parameters for the hat-shaped smoothing function: f(x) = 1.0 - 0.5*(exp(-a*(x-xMin)) + exp(a*(x-xMax))), where a = ln(100)/edgeSize
	// Note: At xMin & xMax, f(x) = 0.0. At the "edges", f(x) = 0.99. At the center of the domain, f(x) ~= 1.0
	//       The smoothing occurs only if useSmoothing is TRUE.
	double disturbSmoothXmin = getDoubleParam("DISTURB_SMOOTH_XMIN", "-1.0e10");
	double disturbSmoothXmax = getDoubleParam("DISTURB_SMOOTH_XMAX", "1.0e10");
	double disturbSmoothEdgeSize = getDoubleParam("DISTURB_SMOOTH_EDGE_SIZE", "1.0e-5");
	double a = log(100.0)/disturbSmoothEdgeSize;
	assert(disturbSmoothXmin<disturbSmoothXmax);
	assert(a>0);

	// -------------------
	// Generate rho_perturb
	// -------------------
	if(randDistribFunc == UNIFORM_DISTRIB) {
		randArrayUniform(array_perturb, ncv); // note: randArrayUniform() generates random variables in [0,1]
		for(int icv=0; icv<ncv; ++icv)
			array_perturb[icv] -= 0.5;        // shifting
	} else
		randArrayNormal(array_perturb, ncv);  // default = normal distribution

//	updateCvDataG1G2(array_perturb, REPLACE_DATA);

	// -------------------
	// Filter rho_perturb
	// -------------------
	// note: parameters for the filtering such as maxFilterWidth and filterBoundaryLength were already set
	if(useFiltering) {
		// Check if filter has been constructed.
		if(!filterConsturcted())
			buildCvDifferentialFilter();
		assert(maxFilterWidth > 0.0);

		applyDiffFilterG1G2(tempFiltered, array_perturb);
		for(int icv=0; icv<ncv; ++icv)
			array_perturb[icv] = tempFiltered[icv];
//		updateCvDataG1G2(array_perturb, REPLACE_DATA);
	}

	// -------------------
	// Update rho
	// -------------------
	double rmsVal = calcRMS(rho, false); // get RMS
	double coeff = rmsVal*disturbMag; // calculate the coefficient based on the RMS
	int myCountClip = 0;
	for(int icv=0; icv<ncv; ++icv) {
		double disturb = coeff*array_perturb[icv];
		if(useSmoothing) {
			double x = x_cv[icv][0];
			if(x<disturbSmoothXmin || x>disturbSmoothXmax)
				disturb = 0.0;
			else
				disturb *= 1.0 - 0.5*(exp(-a*(x-disturbSmoothXmin)) + exp(a*(x-disturbSmoothXmax)));
		}
		if(rho[icv]+disturb < disturbClip*rho[icv]) {
			rho[icv] = disturbClip*rho[icv]; // clipping
			myCountClip++;
		} else {
			rho[icv] += disturb;
		}
	}
	updateCvDataG1G2(rho,  REPLACE_DATA);
	int countClip; 	MPI_Allreduce(&myCountClip, &countClip, 1, MPI_INT, MPI_SUM, mpi_comm);
	if(mpi_rank == 0) {
		if(randDistribFunc == UNIFORM_DISTRIB)
			printf("    Perturb the density field with Uniform:  min disturb = %.3e, max disturb = %.3e, # clip = %d \n", -0.5*coeff, 0.5*coeff, countClip);
		else
			printf("    Perturb the density field with Gaussian:  std = %.3e, # clip = %d \n", coeff, countClip);
	}

	// -------------------
	// Generate rhouX_perturb
	// -------------------
	if(randDistribFunc == UNIFORM_DISTRIB) {
		randArrayUniform(array_perturb, ncv); // note: randArrayUniform() generates random variables in [0,1]
		for(int icv=0; icv<ncv; ++icv)
			array_perturb[icv] -= 0.5;        // shifting
	} else
		randArrayNormal(array_perturb, ncv);  // default = normal distribution

//	updateCvDataG1G2(array_perturb, REPLACE_DATA);

	// -------------------
	// Filter rhouX_perturb
	// -------------------
	// note: parameters for the filtering such as maxFilterWidth and filterBoundaryLength were already set
	if(useFiltering) {
		applyDiffFilterG1G2(tempFiltered, array_perturb);
		for(int icv=0; icv<ncv; ++icv)
			array_perturb[icv] = tempFiltered[icv];
//		updateCvDataG1G2(array_perturb, REPLACE_DATA);
	}

	// -------------------
	// Update x-momentum
	// -------------------
	rmsVal = calcRMS3D(rhou, 0, false); // get RMS
	coeff = rmsVal*disturbMag; // calculate the coefficient based on the RMS
	for(int icv=0; icv<ncv; ++icv) {
		double disturb = coeff*array_perturb[icv];
		if(useSmoothing) {
			double x = x_cv[icv][0];
			if(x<disturbSmoothXmin || x>disturbSmoothXmax)
				disturb = 0.0;
			else
				disturb *= 1.0 - 0.5*(exp(-a*(x-disturbSmoothXmin)) + exp(a*(x-disturbSmoothXmax)));
		}
		rhou[icv][0] += disturb;
	}
	if(mpi_rank == 0) {
		if(randDistribFunc == UNIFORM_DISTRIB)
			printf("    Perturb the momentum field with Uniform: min disturb = %.3e, max disturb = %.3e \n", -0.5*coeff, 0.5*coeff);
		else
			printf("    Perturb the momentum field with Gaussian: x-mom: std = %.3e \n", coeff);
	}

//	// -------------------
//	// Generate rhouY_perturb
//	// -------------------
//	if(randDistribFunc == UNIFORM_DISTRIB) {
//		randArrayUniform(array_perturb, ncv); // note: randArrayUniform() generates random variables in [0,1]
//		for(int icv=0; icv<ncv; ++icv)
//			array_perturb[icv] -= 0.5;        // shifting
//	} else
//		randArrayNormal(array_perturb, ncv);  // default = normal distribution
//
//	updateCvDataG1G2(array_perturb, REPLACE_DATA);
//
//	// -------------------
//	// Filter rhouY_perturb
//	// -------------------
//	// note: parameters for the filtering such as maxFilterWidth and filterBoundaryLength were already set
//	if(useFiltering) {
//		applyDiffFilterG1G2(tempFiltered, array_perturb);
//		for(int icv=0; icv<ncv; ++icv)
//			array_perturb[icv] = tempFiltered[icv];
//		updateCvDataG1G2(array_perturb, REPLACE_DATA);
//	}
//
//	// -------------------
//	// Update y-momentum
//	// -------------------
//	rmsVal = calcRMS3D(rhou, 1, false); // get RMS
//	coeff = rmsVal*disturbMag; // calculate the coefficient based on the RMS
//	for(int icv=0; icv<ncv; ++icv) {
//		double disturb = coeff*array_perturb[icv];
//		if(useSmoothing) {
//			double x = x_cv[icv][0];
//			if(x<disturbSmoothXmin || x>disturbSmoothXmax)
//				disturb = 0.0;
//			else
//				disturb *= 1.0 - 0.5*(exp(-a*(x-disturbSmoothXmin)) + exp(a*(x-disturbSmoothXmax)));
//		}
//		rhou[icv][0] += disturb;
//	}
//	if(mpi_rank == 0) {
//		if(randDistribFunc == UNIFORM_DISTRIB)
//			printf("    Perturb the momentum field with Uniform: min disturb = %.3e, max disturb = %.3e \n", -0.5*coeff, 0.5*coeff);
//		else
//			printf("    Perturb the momentum field with Gaussian: y-mom: std = %.3e \n", coeff);
//	}
//
//	updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);

	// -------------------
	// Generate rhoE_perturb
	// -------------------
	if(randDistribFunc == UNIFORM_DISTRIB) {
		randArrayUniform(array_perturb, ncv); // note: randArrayUniform() generates random variables in [0,1]
		for(int icv=0; icv<ncv; ++icv)
			array_perturb[icv] -= 0.5;        // shifting
	} else
		randArrayNormal(array_perturb, ncv);  // default = normal distribution
//	updateCvDataG1G2(array_perturb, REPLACE_DATA);

	// -------------------
	// Filter rhoE_perturb
	// -------------------
	if(useFiltering) {
		applyDiffFilterG1G2(tempFiltered, array_perturb);
		for(int icv=0; icv<ncv; ++icv)
			array_perturb[icv] = tempFiltered[icv];
//		updateCvDataG1G2(array_perturb, REPLACE_DATA);
	}

	// -------------------
	// Update rhoE
	// -------------------
	int kine_Index = getScalarTransportIndex("kine");

	rmsVal = calcRMS(rhoE, false); // get RMS
	coeff = rmsVal*disturbMag; // calculate the coefficient based on the RMS

	myCountClip = 0;
	for(int icv=0; icv<ncv; ++icv) {
		double interEnergy0 = rhoE[icv] - vecDotVec3d(rhou[icv], rhou[icv])/rho[icv];
		if(kine_Index > -1)
			interEnergy0 -= rho[icv]*scalarTranspEqVector[kine_Index].phi[icv];
		double disturb = coeff*array_perturb[icv];
		if(useSmoothing) {
			double x = x_cv[icv][0];
			if(x<disturbSmoothXmin || x>disturbSmoothXmax)
				disturb = 0.0;
			else
				disturb *= 1.0 - 0.5*(exp(-a*(x-disturbSmoothXmin)) + exp(a*(x-disturbSmoothXmax)));
		}

		double interEnergy = interEnergy0 + disturb;
		if(interEnergy < disturbClip*interEnergy0) {
			rhoE[icv] = disturbClip*interEnergy0 + vecDotVec3d(rhou[icv], rhou[icv])/rho[icv];
			if(kine_Index > -1)
				rhoE[icv] += rho[icv]*scalarTranspEqVector[kine_Index].phi[icv];
			myCountClip++;
		} else {
			rhoE[icv] += disturb;
		}
	}
	updateCvDataG1G2(rhoE,  REPLACE_DATA);
	MPI_Allreduce(&myCountClip, &countClip, 1, MPI_INT, MPI_SUM, mpi_comm);
	if(mpi_rank == 0) {
		if(randDistribFunc == UNIFORM_DISTRIB)
			printf("    Perturb the energy field with Uniform:   min disturb = %.3e, max disturb = %.3e, # clip = %d \n", -0.5*coeff, 0.5*coeff, countClip);
		else
			printf("    Perturb the energy field with Gaussian:   std = %.3e, # clip = %d \n", coeff, countClip);
	}

	// -------------
	// Clear memory
	// -------------
	delete [] array_perturb; array_perturb = NULL;
	if(useFiltering)
		delete [] tempFiltered;

	if(mpi_rank==0)
		cout<<endl;
	MPI_Barrier(mpi_comm);
}

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

// ###########################################################################################
// ------                                                                               ------
// ------                                      MAIN                                     ------
// ------                                                                               ------
// ###########################################################################################

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);

	initMpiStuff();

	// the run number specifies the class which is going to be instantiated
	int run = 0;

	// set input name to default "Ike.in"
	char inputFileName[50];
	sprintf(inputFileName, "VanGogh.in");

	for (int i=1; i<argc; i++) {
		string str(argv[i]);
		if (from_string<int>(run, str, std::dec)) {
			if (mpi_rank == 0)
				cout << "You have specified run number = " << run << endl;
		}
		else
			strcpy(inputFileName, argv[i]);
	}

	if (mpi_rank == 0) {
		cout << "SPECIFIED INPUT NAME = " << inputFileName << endl;
		cout << "SPECIFIED RUN = " << run << endl;
	}


	try {
		// provide total runtime
		double wtime, wtime0;
		MPI_Barrier(mpi_comm);
		if (mpi_rank == 0)
			wtime = MPI_Wtime();

		// declare pointer to JoeWithModels_AD
		JoeWithModels *vanGogh;

		switch (run) {
			case 1:  vanGogh = new VanGogh(inputFileName); 	break;
			default:
				if (mpi_rank == 0)
					cerr << "WARNING: run number not available!" << endl;
				vanGogh = new VanGogh(inputFileName);
		}

		// run vanGogh
		vanGogh->run();

		MPI_Barrier(mpi_comm);
		if (mpi_rank == 0) {
			double wtime0 = wtime;
			wtime = MPI_Wtime();
			cout << " > total runtime [s]: " << wtime - wtime0 << endl;
		}

		// delete vanGogh (make sure memory is deallocated in destructors
		delete vanGogh;
	}
	catch (int e) {
		cerr << "Exception: " << e << endl;
		MPI_Finalize();
		return(-1);
	}
	catch(...) {
		cerr << "unhandled Exception.\n" << endl;
		MPI_Finalize();
		return(-1);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return (0);
}

