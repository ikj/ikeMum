#include "rembrandtWithModels.h"

// ###########################################################################################
// ------                                                                               ------
// ------   Calculate eigen-decomposition with SLEPSC and Export the data for TECPLOT   ------
// ------                                                                               ------
// ###########################################################################################

/*
 * Method: init
 * ------------
 *
 */
void RembrandtWithModels::init() {
	// -------------------------
	// Global-modes calculations
	// -------------------------
	directEvalsReal = NULL; 	directEvalsImag = NULL;
	directEvecsReal = NULL; 	directEvecsImag = NULL;

	adjointEvalsReal = NULL; 	adjointEvalsImag = NULL;
	adjointEvecsReal = NULL; 	adjointEvecsImag = NULL;

	// -------------------------
	// Tecplot output
	// -------------------------
	FirstDirectMode_real_rho  = NULL; 	registerScalar(FirstDirectMode_real_rho,  "FirstDirectMode_REAL_RHO",  CV_DATA);
	FirstDirectMode_real_rhou = NULL;	registerVector(FirstDirectMode_real_rhou, "FirstDirectMode_REAL_RHOU", CV_DATA);
	FirstDirectMode_real_vel  = NULL; 	registerVector(FirstDirectMode_real_vel,  "FirstDirectMode_REAL_VEL",  CV_DATA);
	FirstDirectMode_real_rhoE = NULL; 	registerScalar(FirstDirectMode_real_rhoE, "FirstDirectMode_REAL_RHOE", CV_DATA);
	// NOTE: scalars should be registered in the child class

	FirstAdjointMode_real_rho  = NULL; 	registerScalar(FirstAdjointMode_real_rho,  "FirstAdjointMode_REAL_RHO",  CV_DATA);
	FirstAdjointMode_real_rhou = NULL; 	registerVector(FirstAdjointMode_real_rhou, "FirstAdjointMode_REAL_RHOU", CV_DATA);
	FirstAdjointMode_real_vel  = NULL; 	registerVector(FirstAdjointMode_real_vel,  "FirstAdjointMode_REAL_VEL",  CV_DATA);
	FirstAdjointMode_real_rhoE = NULL; 	registerScalar(FirstAdjointMode_real_rhoE, "FirstAdjointMode_REAL_RHOE", CV_DATA);
	// NOTE: scalars should be registered in the child class

	rho_sens  = NULL;	registerScalar(rho_sens,  "SENSITIVITY_RHO",  CV_DATA);
	rhou_sens = NULL;	registerVector(rhou_sens, "SENSITIVITY_RHOU", CV_DATA);
	vel_sens  = NULL;	registerVector(vel_sens,  "SENSITIVITY_VEL",  CV_DATA);
	rhoE_sens = NULL;	registerScalar(rhoE_sens, "SENSITIVITY_RHOE", CV_DATA);
	// NOTE: scalars should be registered in the child class
}

/*
 * Method: clear
 * -------------
 *
 */
void RembrandtWithModels::clear() {
	// -------------------------
	// Global-modes calculations
	// -------------------------
	if(directEvalsReal != NULL) {
		delete [] directEvalsReal; 	directEvalsReal = NULL;
	}
	if(directEvalsImag != NULL) {
		delete [] directEvalsImag; 	directEvalsImag = NULL;
	}
	if(directEvecsReal != NULL) {
		for(int i=0; i<nev; ++i)
			delete [] directEvecsReal[i];
		delete [] directEvecsReal; 	directEvecsReal = NULL;
	}
	if(directEvecsImag != NULL) {
		for(int i=0; i<nev; ++i)
			delete [] directEvecsImag[i];
		delete [] directEvecsImag; 	directEvecsImag = NULL;
	}

	if(adjointEvalsReal != NULL) {
		delete [] adjointEvalsReal; 	adjointEvalsReal = NULL;
	}
	if(adjointEvalsImag != NULL) {
		delete [] adjointEvalsImag; 	adjointEvalsImag = NULL;
	}
	if(adjointEvecsReal != NULL) {
		for(int i=0; i<nev; ++i)
			delete [] adjointEvecsReal[i];
		delete [] adjointEvecsReal; 	adjointEvecsReal = NULL;
	}
	if(adjointEvecsImag != NULL) {
		for(int i=0; i<nev; ++i)
			delete [] adjointEvecsImag[i];
		delete [] adjointEvecsImag; 	adjointEvecsImag = NULL;
	}

	// -------------------------
	// SLEPc context
	// -------------------------
	if (petscSolver2 != NULL) {
		delete petscSolver2;  	petscSolver2 = NULL;
	}
	if (slepcSolver2 != NULL) {
		delete slepcSolver2;  	slepcSolver2 = NULL;
	}
}

/*
 * Method: run
 * -----------
 *
 */
void RembrandtWithModels::run() {
	nScal = scalarTranspEqVector.size();

	debugLevel = getDebugLevel();

	// +++++++++++++++++++++++++++++++++++++
	// Load the mesh or the restart file
	// +++++++++++++++++++++++++++++++++++++
	if(mpi_rank==0){
		cout<<"======================="<<endl;
		cout<<"LOAD THE INITIAL FIELDS"<<endl;
		cout<<"======================="<<endl;
	}
	// Read mesh or restart file
	initializeFromRestartFile(getStringParam("RESTART"));
	// initialize models: Register scalars if not started from a restart file
	initialHookScalarRansTurbModel();
	initialHookScalarRansCombModel();

	// Make sure ghost cells are populated
	updateCvDataG1G2(rho, REPLACE_DATA);
	updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
	updateCvDataG1G2(rhoE, REPLACE_DATA);
	for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
		updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);

	// Construct (cvora &) nbocv_v_global
	if(mpi_rank==0)
		cout<<endl;
	if(cvora == NULL) {
		if(mpi_rank==0)
			cout<<"CAUTION!! cvora is NULL -- Remember that you should initialize it later"<<endl;
	}
	if (nbocv_v_global == NULL) {
		nbocv_v_global = new int[ncv_gg];
		for (int icv = 0; icv < ncv; icv++)
			nbocv_v_global[icv] = cvora[mpi_rank] + icv;
		updateCvDataG1G2(nbocv_v_global, REPLACE_DATA);

		if(mpi_rank==0)
			cout<<"nbocv_v_global is allocated and initialized"<<endl<<endl;
	}
	// Build the 2-layered CSR structure (nbocv2_i & nbocv2_v)
	if(mpi_rank==0){
		cout<<"======================="<<endl;
		cout<<"BUILD THE 2-CSR STRUCT."<<endl;
		cout<<"======================="<<endl;
	}
	bool addFakeCells = false;
	build2layerCSRstruct(addFakeCells);
	MPI_Barrier(mpi_comm);
	if(mpi_rank==0)
		cout<<"nbocv2_v_global, nbocv2_i, and nbocv2_v are successfully constructed"<<endl<<endl;

	// Load the reference flow values
	getReferenceParams(RefFlowParams);
	if(mpi_rank==0)
		RefFlowParams.showOnScreen();

#ifdef USE_VOLUME_WEIGHTED_INNER_PROD
	// Get the maximum and minimum cv_volume: it will be used for the weights in the inner product
	double myMaxVolume = -ABSURDLY_BIG_NUMBER;
	double myMinVolume = ABSURDLY_BIG_NUMBER;
	for(int icv=0; icv<ncv; ++icv) {
		myMaxVolume = max(myMaxVolume, cv_volume[icv]);
		myMinVolume = min(myMinVolume, cv_volume[icv]);
	}
	MPI_Allreduce(&myMaxVolume, &totMaxCvVol, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
	MPI_Allreduce(&myMinVolume, &totMinCvVol, 1, MPI_DOUBLE, MPI_MIN, mpi_comm);
#endif

	step = 0;

	string tempString = getStringParam("INITIAL_FLOWFIELD_OUTPUT", "NO");
	std::transform(tempString.begin(), tempString.end(), tempString.begin(), ::tolower);
	if(tempString.compare("yes")==0 || tempString.compare("y")==0)
		writeData(0);

	// +++++++++++++++++++++++++++++++++++++
	// Obtain the Jacobian matrix
	// +++++++++++++++++++++++++++++++++++++
	++step;

	NcontrolEqns = getNumControlVar(); 		assert(NcontrolEqns >= 0);

	if(mpi_rank==0){
		cout<<endl;
		cout<<"======================="<<endl;
		cout<<"JACOBIAN MATRIX"<<endl;
		cout<<"======================="<<endl;
	}

	// Obtain the Jacobian matrix
	obtainJacobianMat();

	int myNrows = jacMatrixSTL.get_nRows();
	int myNnz   = jacMatrixSTL.get_nnz();
	int totNrows, totNnz;
	MPI_Allreduce(&myNrows, &totNrows, 1, MPI_INT, MPI_SUM, mpi_comm);  // Note: Getting totNcols doesn't say anything meaningful.
	MPI_Allreduce(&myNnz,   &totNnz,   1, MPI_INT, MPI_SUM, mpi_comm);
	if(mpi_rank==0)
		cout<<"> "<<classID<<"::run(): MatComprsedSTL total Nrows = "<<totNrows<<", total nnz = "<<totNnz<<endl
		    <<endl;

	// +++++++++++++++++++++++++++++++++++++
	// Launch eigen-analysis
	// +++++++++++++++++++++++++++++++++++++
	if(mpi_rank==0){
		cout<<endl;
		cout<<"======================="<<endl;
		cout<<"EIGEN-DECOMPOSITION"<<endl;
		cout<<"======================="<<endl;
	}

	// Perform the eigen-analysis
	runEigenAnalysis();

	// Write the result on a tecplot file
	writeData(step);

	// +++++++++++++++++++++++++++++++++++++
	// End of Rembrandt
	// +++++++++++++++++++++++++++++++++++++
	// Pause for a while
	MPI_Barrier(mpi_comm);
	for(int i=0; i<100000; ++i) {};
	MPI_Barrier(mpi_comm);

	if(mpi_rank==0){
		cout<<endl;
		cout<<"======================="<<endl;
		cout<<"CLOSING REMBRANDT"<<endl;
		cout<<"======================="<<endl;
	}
}

/*
 * Method: obtainJacobianMat
 * -------------------------
 * Obtain the Jacobian matrix: There are two possible options --
 *   1. Calculate the Jacobian matrix using ADOL-C
 *   2. Read the Jacobian matrix from a file (still in development)
 */
void RembrandtWithModels::obtainJacobianMat() {
	int nVars = 5+nScal;

	string tempString = getStringParam("HOW_TO_GET_JACOBIAN", "NONE");
	std::transform(tempString.begin(), tempString.end(), tempString.begin(), ::toupper);

	HOW_TO_GET_JACOBIAN HowToGetJacobian;
	if(tempString.compare("FROM_MAT_BINARY_FILE") == 0)    HowToGetJacobian = FROM_MAT_BINARY_FILE;
	else if(tempString.compare("CALC_USING_ADOLC") == 0)   HowToGetJacobian = CALC_USING_ADOLC;
	else {
		if(mpi_rank==0) cerr<<"ERROR in "<<classID<<"::obtainJacobianMat(): HOW_TO_GET_JACOBIAN = "<<tempString<<" is not supported"<<endl;
		throw(REMBRANDT_ERROR_CODE);
	}

	// Variable  for "FROM_MAT_BINARY_FILE" :
	string Q1Filename;
	string JacFilename;

	switch(HowToGetJacobian) {
	case FROM_MAT_BINARY_FILE:
		if(mpi_rank == 0)
			cerr<<"ERROR in "<<classID<<"::obtainJacobianMat(): HOW_TO_GET_JACOBIAN = "<<tempString<<" has not been implemented yet!"<<endl;
		throw(REMBRANDT_ERROR_CODE);

//		// Read the matrix: Note that the number of CPU cores must match
//		Q1Filename  = getStringParam("Q1_FILENAME", "DEFAULT_FILENAME_NONE");
//		JacFilename = getStringParam("JACOBIAN_MATRIX_FILENAME", "DEFAULT_FILENAME_NONE");
//		readMatrixBinaryRembrandt<MatComprsedSTL>(Q1Filename, JacFilename, jacMatrixSTL, nScal, cvora, nbocv_v_global);

		break;
	case CALC_USING_ADOLC:
		initialHook();
		calcMatrixFromField();
		finalHook();

		break;
	default:
		break;
	}
}

/*
 * Method: runEigenAnalysis
 * ------------------------
 * Calculate both direct and adjoint global-modes and perform the required post-processing.
 * Return = the number of converged eigenpairs (min in Direct and Adjont)
 */
int RembrandtWithModels::runEigenAnalysis() {
	int nVars = 5+nScal;

	// ==============================================
	// Calculate both direct and adjoint global-modes
	// ==============================================
	nev = getIntParam("NUM_EIGEN_PAIRS", "3");

	directEvalsReal = new double [nev];
	directEvalsImag = new double [nev];
	directEvecsReal = new double* [nev];
	for(int i=0; i<nev; ++i)
		directEvecsReal[i] = new double [ncv*(5+nScal)];
	directEvecsImag = new double* [nev];
	for(int i=0; i<nev; ++i)
		directEvecsImag[i] = new double [ncv*(5+nScal)];

	adjointEvalsReal = new double [nev];
	adjointEvalsImag = new double [nev];
	adjointEvecsReal = new double* [nev];
	for(int i=0; i<nev; ++i)
		adjointEvecsReal[i] = new double [ncv*(5+nScal)];
	adjointEvecsImag = new double* [nev];
	for(int i=0; i<nev; ++i)
		adjointEvecsImag[i] = new double [ncv*(5+nScal)];

	// Cacluate the global-modes by calling calcGlobalModes()
	int nconv = calcGlobalModes(nev);  // Calculate both direct and adjoint eigen-pairs.
	                                   // This method also writes the eigenvalues on a file.

	// ===============
	// Post-processing
	// ===============
	if(nconv > 0) {
		int directLeastStableIndex = -1;
		double directLeastStableEigenReal = -ABSURDLY_BIG_NUMBER;
		for(int i=0; i<nconv; ++i)
			if(directLeastStableEigenReal < directEvalsReal[i]) {
				directLeastStableIndex     = i;
				directLeastStableEigenReal = directEvalsReal[i];
			}

		int adjointLeastStableIndex = -1;
		double adjointLeastStableEigenReal = -ABSURDLY_BIG_NUMBER;
		for(int i=0; i<nconv; ++i)
			if(adjointLeastStableEigenReal < adjointEvalsReal[i]) {
				adjointLeastStableIndex     = i;
				adjointLeastStableEigenReal = adjointEvalsReal[i];
			}

		assert(directLeastStableIndex >= 0   && adjointLeastStableIndex >= 0);
		assert(directLeastStableIndex <  nev && adjointLeastStableIndex <  nev);

		if(mpi_rank == 0)
			cout<<">> Start post-processing: least stable directmode : index = "<<directLeastStableIndex <<", eigenvalue (Real part) = "<<directLeastStableEigenReal<<endl
				<<"                          least stable adjointmode: index = "<<adjointLeastStableIndex<<", eigenvalue (Real part) = "<<adjointLeastStableEigenReal<<endl;

		if(mpi_rank == 0)
			cout<<endl
			    <<"WARNING! You may want to check the current formulation for sensitivity in Rembrandt!"<<endl
			    <<endl;

		// Post-processing for rho, rhou, and rhoE
		for(int icv=0; icv<ncv; ++icv) {
			int indexTemp = icv*nVars;

			// -------------------------------------------------
			// Tecplot output for the least-stable direct modes
			// -------------------------------------------------
			FirstDirectMode_real_rho[icv]   = directEvecsReal[directLeastStableIndex][indexTemp];
			for(int i=0; i<3; ++i) {
				FirstDirectMode_real_rhou[icv][i] = directEvecsReal[directLeastStableIndex][indexTemp+1+i];
				FirstDirectMode_real_vel[icv][i]  = directEvecsReal[directLeastStableIndex][indexTemp+1+i]/FirstDirectMode_real_rho[icv];
			}
			FirstDirectMode_real_rhoE[icv]  = directEvecsReal[directLeastStableIndex][indexTemp+4];

			// -------------------------------------------------
			// Tecplot output for the least-stable adjoint modes
			// -------------------------------------------------
			FirstAdjointMode_real_rho[icv]   = adjointEvecsReal[adjointLeastStableIndex][indexTemp];
			for(int i=0; i<3; ++i) {
				FirstAdjointMode_real_rhou[icv][i] = adjointEvecsReal[adjointLeastStableIndex][indexTemp+1+i];
				FirstAdjointMode_real_vel[icv][i]  = adjointEvecsReal[adjointLeastStableIndex][indexTemp+1+i]/FirstAdjointMode_real_rho[icv];
			}
			FirstAdjointMode_real_rhoE[icv]  = adjointEvecsReal[adjointLeastStableIndex][indexTemp+4];

			// -------------------------------------------------
			// Structural sensitivity
			// -------------------------------------------------
//			rho_sens[icv] = calcComplexMag(FirstDirectMode_real_rho[icv], directEvecsImag[directLeastStableIndex][indexTemp])
//					* calcComplexMag(FirstAdjointMode_real_rho[icv], adjointEvecsImag[directLeastStableIndex][indexTemp]) ;
//			for(int i=0; i<3; ++i) {
//				rhou_sens[icv][i] = calcComplexMag(FirstDirectMode_real_rhou[icv][i], directEvecsImag[directLeastStableIndex][indexTemp+1+i])
//						* calcComplexMag(FirstAdjointMode_real_rhou[icv][i], adjointEvecsImag[directLeastStableIndex][indexTemp+1+i]) ;
//
//				double firstDirectMode_imag_vel = directEvecsImag[directLeastStableIndex][indexTemp+1+i] / directEvecsImag[directLeastStableIndex][indexTemp];
//				double firstAdjointMode_imag_vel = adjointEvecsImag[directLeastStableIndex][indexTemp+1+i] / adjointEvecsImag[directLeastStableIndex][indexTemp];
//				vel_sens[icv][i]  = calcComplexMag(FirstDirectMode_real_vel[icv][i], firstDirectMode_imag_vel)
//						* calcComplexMag(FirstAdjointMode_real_vel[icv][i], firstAdjointMode_imag_vel);
//			}
//			rhoE_sens[icv] = calcComplexMag(FirstDirectMode_real_rhoE[icv], directEvecsImag[directLeastStableIndex][indexTemp+4])
//							* calcComplexMag(FirstAdjointMode_real_rhoE[icv], adjointEvecsImag[directLeastStableIndex][indexTemp+4]) ;

			// Because the least stable modes are real anyways, we don't need to worry about the imaginary parts anyways
			rho_sens[icv] = fabs(COEFF_BOOST_SENS * FirstDirectMode_real_rho[icv] * FirstAdjointMode_real_rho[icv]);
			for(int i=0; i<3; ++i) {
				rhou_sens[icv][i] = fabs(COEFF_BOOST_SENS * FirstDirectMode_real_rhou[icv][i] * FirstAdjointMode_real_rhou[icv][i]);
				vel_sens[icv][i]  = fabs(COEFF_BOOST_SENS * FirstDirectMode_real_vel[icv][i] * FirstAdjointMode_real_vel[icv][i]);
			}
			rhoE_sens[icv] = fabs(COEFF_BOOST_SENS * FirstDirectMode_real_rhoE[icv] * FirstAdjointMode_real_rhoE[icv]);
		}

		updateCvDataG1G2(FirstDirectMode_real_rho,  REPLACE_DATA);
		updateCvDataG1G2(FirstDirectMode_real_rhou, REPLACE_ROTATE_DATA);
		updateCvDataG1G2(FirstDirectMode_real_vel,  REPLACE_ROTATE_DATA);
		updateCvDataG1G2(FirstDirectMode_real_rhoE, REPLACE_DATA);

		updateCvDataG1G2(FirstAdjointMode_real_rho,  REPLACE_DATA);
		updateCvDataG1G2(FirstAdjointMode_real_rhou, REPLACE_ROTATE_DATA);
		updateCvDataG1G2(FirstAdjointMode_real_vel,  REPLACE_ROTATE_DATA);
		updateCvDataG1G2(FirstAdjointMode_real_rhoE, REPLACE_DATA);

		updateCvDataG1G2(rho_sens,  REPLACE_DATA);
		updateCvDataG1G2(rhou_sens, REPLACE_ROTATE_DATA);
		updateCvDataG1G2(vel_sens,  REPLACE_ROTATE_DATA);
		updateCvDataG1G2(rhoE_sens, REPLACE_DATA);

		// Post-processing for scalars in TurbModels and CombModels
		postEigenDecompScalarRansTurbModel(directEvecsReal, adjointEvecsReal, directLeastStableIndex, adjointLeastStableIndex);
		postEigenDecompScalarRansCombModel(directEvecsReal, adjointEvecsReal, directLeastStableIndex, adjointLeastStableIndex);

		// Check the maximum and minimum values
		showDataRangeEigenAnalysis();
	}

	// =================
	// WRITE BINARY FILE
	// =================
	if(nconv > 0) {
		stringstream ss;
		ss<<"EigenPairs.pt"<<step<<".bin";
		int nevWrite = getIntParam("NEV_WRITE_BINDARY", "1"); 	assert(nevWrite > 0);
		writeEigenPairsParallel(ss.str(), step, nevWrite);
	}

	return nconv;
}

/*
 * Method: showDataRangeEigenAnalysis
 * ----------------------------------
 * Show data range of the least-stable direct-mode, the least-stable ajoint-mode, and the sensitivity
 */
void RembrandtWithModels::showDataRangeEigenAnalysis(){
	double my_sens_max[5];
	double my_sens_min[5];
	for(int icv=0; icv<ncv; ++icv) {
		my_sens_max[0] = max(my_sens_max[0], rho_sens[icv]);
		for(int i=0; i<3; ++i)
			my_sens_max[1+i] = max(my_sens_max[1+i], rhou_sens[icv][i]);
		my_sens_max[4] = max(my_sens_max[4], rhoE_sens[icv]);

		my_sens_min[0] = min(my_sens_min[0], rho_sens[icv]);
		for(int i=0; i<3; ++i)
			my_sens_min[1+i] = min(my_sens_min[1+i], rhou_sens[icv][i]);
		my_sens_min[4] = min(my_sens_min[4], rhoE_sens[icv]);
	}

	double sens_max[5];
	double sens_min[5];
	MPI_Allreduce(my_sens_max, sens_max, 5, MPI_DOUBLE, MPI_MAX, mpi_comm);
	MPI_Allreduce(my_sens_min, sens_min, 5, MPI_DOUBLE, MPI_MIN, mpi_comm);

	if(mpi_rank == 0) {
		cout<<endl
		    <<">> Data range  rho_sens  = "<<sens_min[0]<<" : "<<sens_max[0]<<endl
		    <<"               rhou_sens = "<<"("<<sens_min[1]<<","<<sens_min[2]<<","<<sens_min[3]<<") : ("<<sens_max[1]<<","<<sens_max[2]<<","<<sens_max[3]<<")"<<endl
		    <<"               rhoE_sens = "<<sens_min[4]<<" : "<<sens_max[4]<<endl;
	}
}

/*
 * Method: readMatrixBinaryRembrandt
 * ---------------------------------
 * Original function = IkeWithPsALC_AD::readMatrixBinary()
 * For some reason, "undefined reference to" error occurs by the linker. Thus, define the same function in Rembrandt.
 */
template <class MatT>
void RembrandtWithModels::readMatrixBinaryRembrandt(string &Q1Filename, string &JacFilename, MatT &jacMatrix,
		const int nScal, const int *cvora, const int *cv_gl) {
	// You must update xcvMin and xcvMax before reading the Q1 file and the Jacobian file
	for(int i=0; i<3; ++i) {
		xcvMin[i] = ABSURDLY_BIG_NUMBER;
		xcvMax[i] = -ABSURDLY_BIG_NUMBER;
	}
	for(int icv=0; icv<ncv; ++icv) {
		for(int i=0; i<3; ++i) {
			if(x_cv[icv][i] < xcvMin[i])
				xcvMin[i] = x_cv[icv][i];
			if(x_cv[icv][i] > xcvMax[i])
				xcvMax[i] = x_cv[icv][i];
		}
	}
	for(int i=0; i<3; ++i) { // Add or subtract the grid tolerance for safety
		xcvMin[i] -= GRID_TOL;
		xcvMax[i] += GRID_TOL;
	}

	/***********
	 ** Read the Q1 file first to get the mesh information of the IKE simulation
	 ***********/
	int* mpi_rank_inFile = new int [ncv];
	int* icv_inFile      = new int [ncv];
	set<int> Set_mpi_rank_inFile;

	if(mpi_rank==0) cout<<">> readMatrixBinaryRembrandt(): Get mesh information from  "<<Q1Filename<<endl;

	readMatchedPtsFromPsALCdumpedDataSerial(Q1Filename, mpi_rank_inFile, icv_inFile, Set_mpi_rank_inFile, NcontrolEqns);


	/***********
	 ** Read the Jacobian matrix
	 ***********/
	if(mpi_rank==0) cout<<">> readMatrixBinaryRembrandt(): Initialize Jacobian matrix from  "<<JacFilename<<endl;

	int nVars = 5 + nScal;

	// Read the data
	int nnz_local;
	vector<unsigned int> rind_local_unsorted;
	vector<unsigned int> cind_local_unsorted;
	vector<double> values_unsorted;

	// To Do: The following function has not been developed yet!
//		readMatrixBinaryForRembrandtParallel(nnz_local, rind_local_unsorted, cind_local_unsorted, values_unsorted, JacFilename, mpi_rank_inFile, icv_inFile, Set_mpi_rank_inFile);

	/***********
	 ** Store the matrix in the data container
	 ***********/
	assert(jacMatrix.empty());

	if(debugLevel>0 && mpi_rank==0)
		cout<<"    > Storing the matrix in the data container ... ";

	jacMatrix.setMatSize(ncv_gg*nVars, ncv*nVars);

	jacMatrix.matCopy(nnz_local, rind_local_unsorted, cind_local_unsorted, values_unsorted);

	if(debugLevel>0 && mpi_rank==0)
		cout<<"DONE!"<<endl;

	/***********
	 ** Clear the memory
	 ***********/
	rind_local_unsorted.clear();
	cind_local_unsorted.clear();
	values_unsorted.clear();

	delete [] mpi_rank_inFile;
	delete [] icv_inFile;

	MPI_Barrier(mpi_comm);
}

/*
 * Method: readMatchedPtsFromPsALCdumpedDataSerial
 * -----------------------------------------------
 * Read data from previous IKE simulation
 */
void RembrandtWithModels::readMatchedPtsFromPsALCdumpedDataSerial(const string &filename,
		int* mpi_rank_inFile, int* icv_inFile, set<int>& Set_mpi_rank_inFile, const int NcontrolEqns) {
	char *filenameArray = new char[filename.size() + 1]; // Since MPI_File_open() doesn't accept "const char *", you cannot use std::string.c_str()
	std::copy(filename.begin(), filename.end(), filenameArray);
	filenameArray[filename.size()] = '\0';
	readMatchedPtsFromPsALCdumpedDataSerial(filenameArray, mpi_rank_inFile, icv_inFile, Set_mpi_rank_inFile, NcontrolEqns);

	delete [] filenameArray;
}

void RembrandtWithModels::readMatchedPtsFromPsALCdumpedDataSerial(const char filename[],
		int* mpi_rank_inFile, int* icv_inFile, set<int>& Set_mpi_rank_inFile, const int NcontrolEqns) {
	// Check some possible errors
	for(int i=0; i<3; ++i)
		if(xcvMax[i] < xcvMin[i]) {
			cerr<<"ERROR in "<<classID<<"::readMatchedPtsFromPsALCdumpedDataSerial(): mpi_rank == "<<mpi_rank<<" : xcvMax["<<i<<"] < xcvMin["<<i<<"]"<<endl;
			throw(REMBRANDT_ERROR_CODE);
		}
	assert(mpi_rank_inFile != NULL);
	assert(icv_inFile      != NULL);

	// Initialize found_mpi_rank & found_icv with a negative number
	for(int icv=0; icv<ncv; ++icv) {
		mpi_rank_inFile[icv] = -1;
		icv_inFile[icv]      = -1;
	}

	// number of variables
	int nScal = scalarTranspEqVector.size();
	int m = nScal+5; // number of variables (rho, rhou, rhov, rhow, rhoE, and scalars)

	// dummy variables for the file information
	int step0;
	double* lambda0 = NULL;
	double* lambda1 = NULL;
	if(NcontrolEqns > 0) {
		lambda0 = new double [NcontrolEqns];
		lambda1 = new double [NcontrolEqns];
	}
	double dsTry;

	// variables related to the Mesh
	const double gridTol = 1.0e-9;

	int *countFound = new int [ncv];
	for(int icv=0; icv<ncv; ++icv)
		countFound[icv] = 0;

	// Read the file
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
		//              3. lambda0 (double*NcontrolEqns)   4. lambda1 (double*NcontrolEqns)
		//              5. ds (double)
		//              6. mpi_size (int)                  7. cvora (double*(mpi_size+1)
		//              8. xMinArray (double*3*mpi_size)   9. xMaxArray (double*3*mpi_size)
		infile.read(reinterpret_cast<char*>(&step0), sizeof(int)); // Note: if you simply use the ">>" operator, it will cause some problem due to bugs in a g++ library on linux
		infile.read(reinterpret_cast<char*>(&NcontrolEqns_file), sizeof(int));
		if(NcontrolEqns_file != NcontrolEqns) {
			if(mpi_rank==0)
				printf("ERROR in %s::readMatchedPtsFromPsALCdumpedDataSerial(): NcontrolEqns from file=%d is not equal to given NcontrolEqns=%d\n", classID.c_str(), NcontrolEqns_file, NcontrolEqns);
			delete [] countFound;
			throw(REMBRANDT_ERROR_CODE);
		}
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
			printf("\nReading \'%s\': step=%d, lambda0=", filename, step0);
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
				printf("%.8f, ", lambda0[iEqn]);
			printf("lambda1=");
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
				printf("%.8f, ", lambda1[iEqn]);
			printf("ds=%.8f, total number of CVs=%d \n", dsTry, cvora_file[mpi_size_file]);
			cout<<"  Grid tolerance = "<<gridTol<<endl;

			if(debugLevel > 1) {
				cout<<endl;
				for(int mpi_file=0; mpi_file<mpi_size_file; ++mpi_file)
					printf("  RANK in the file = %d:  NCV = %d,  X = %e~%e, Y = %e~%e, Z = %e~%e\n", mpi_file, cvora_file[mpi_file+1]-cvora_file[mpi_file],
							xcvMinArray_file[mpi_file*3], xcvMaxArray_file[mpi_file*3],  xcvMinArray_file[mpi_file*3+1], xcvMaxArray_file[mpi_file*3+1],  xcvMinArray_file[mpi_file*3+2], xcvMaxArray_file[mpi_file*3+2]);
				cout<<endl;
			}
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
						mpi_rank_inFile[icv] = mpi_file;
						icv_inFile[icv]      = icv_file;
						if(Set_mpi_rank_inFile.count(mpi_file) == 0)
							Set_mpi_rank_inFile.insert(mpi_file);

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
				cout<<"ERROR! "<<classID<<"::readMatchedPtsFromPsALCdumpedDataSerial(): file does not end with EOF_ERROR_CHECK_CODE="<<EOF_ERROR_CHECK_CODE<<endl;
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

		// 5. Check if there is still a negative number for mpi_rank_inFile & icv_inFile
		for (int icv=0; icv<ncv; ++icv) {
			if(mpi_rank_inFile[icv] < 0 || icv_inFile[icv] < 0) {
				printf("  Cannot find a matched mpi_rank and/or icv from the file for icv=%d (%.6f,%.6f,%.6f) at mpi_rank=%d \n", icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2],mpi_rank);
				assert(false);
			}
		}
	}  else {
		if(mpi_rank==0)
			printf("Cannot read \'%s\' \n", filename);
		assert(false);
	}

	infile.close();

	delete [] countFound;
	delete [] cvora_file;
	delete [] xcvMinArray_file;
	delete [] xcvMaxArray_file;

	if(lambda0 != NULL) 	delete [] lambda0;
	if(lambda1 != NULL) 	delete [] lambda1;

	MPI_Barrier(mpi_comm);
}

//	/*
//	 * Method: readMatrixBinaryForRembrandtParallel
//	 * --------------------------------------------
//	 * Read a matrix from a binary file using MPI-2 parallel I/O
//	 * Output: myNnz_fromFile = local number of non-zeros (NOT the total nnz -- it is the nnz for a core)
//	 *         rind           = GLOBAL row indices array
//	 *         cind           = GLOBAL column indices array
//	 *         values         = values array
//	 */
//	void RembrandtWithModels::readMatrixBinaryForRembrandtParallel(int &nnz_local,
//			vector<unsigned int> &rind_local_unsorted, vector<unsigned int> &cind_local_unsorted, vector<double> &values_unsorted,
//			string &JacFilename, const int* mpi_rank_inFile, const int* icv_inFile, set<int>& Set_mpi_rank_inFile,
//			const int nVars, const int *cvora_inFile, const int *cv_gl_inFile) {
//		/***********
//		 ** Read the data
//		 ***********/
//		char *filenameArray = new char[JacFilename.size() + 1]; // Since MPI_File_open() doesn't accept "const char *", you cannot use std::string.c_str()
//		std::copy(JacFilename.begin(), JacFilename.end(), filenameArray);
//		filenameArray[JacFilename.size()] = '\0';
//
//		/*
//		 * 1. Read the header part
//		 */
//		int nRows_fromFile = -1;
//		int nCols_fromFile = -1;
//		int nnz_fromFile = -1;
//		int NcontrolEqns_fromFile = -1;
//		int mpi_size_fromFile = -1;
//		int *localNnzArray_fromFile = NULL;
//
//		ifstream infile(filename.c_str(), ios::in | ios::binary);
//		if(infile.is_open()) {
//			infile.read(reinterpret_cast<char*>(&nRows_fromFile), sizeof(int)); // Note: if you simply use the ">>" operator, it will cause some problem due to bugs in a g++ library on linux
//			infile.read(reinterpret_cast<char*>(&nCols_fromFile), sizeof(int));
//			infile.read(reinterpret_cast<char*>(&nnz_fromFile),   sizeof(int));
//			infile.read(reinterpret_cast<char*>(&NcontrolEqns_fromFile), sizeof(int));
//			if(NcontrolEqns_fromFile != NcontrolEqns) {
//				if(mpi_rank==0)
//					printf("ERROR in %s::readMatrixBinaryForRembrandtParallel: NcontrolEqns from file=%d is not equal to given NcontrolEqns=%d\n",
//							classID.c_str(), NcontrolEqns_fromFile, NcontrolEqns);
//				delete [] filenameArray;
//				throw(REMBRANDT_ERROR_CODE);
//			}
//
//			infile.read(reinterpret_cast<char*>(&mpi_size_fromFile), sizeof(int));
//			assert(mpi_size_fromFile>0);
//
//			localNnzArray_fromFile = new int [mpi_size_fromFile];
//			infile.read(reinterpret_cast<char*>(localNnzArray_fromFile), sizeof(int)*mpi_size_fromFile);
//			if(debugLevel>0 && mpi_rank==0)
//				cout<<"    > Reading the Jacobian matrix from "<<filename<<": nRows="<<nRows_fromFile<<", nCols="<<nCols_fromFile<<", nnz="<<nnz_fromFile<<endl;
//			if(debugLevel>0 && mpi_rank==0) {
//				cout<<"      > local nnz =";
//				for(int i=0; i<mpi_size_fromFile; ++i)
//					cout<<" "<<localNnzArray_fromFile[i]<<" ";
//				cout<<endl;
//			}
//
//			assert(nRows_fromFile>0 && nCols_fromFile>0 && nnz_fromFile>0);
//			assert(localNnzArray_fromFile[mpi_rank] > 0);
//		} else {
//			if(mpi_rank==0)
//				printf("Cannot read \'%s\' \n", JacFilename.c_str());
//			assert(false);
//		}
//		infile.close();
//		MPI_Barrier(mpi_comm);
//
////		myNnz_fromFile = localNnzArray_fromFile[mpi_rank];
//		nnz_local = 0;
//		int initDisp = sizeof(int)*(5 + mpi_size_fromFile); // initial displacement
//
//		/*
//		 * 2. Read the body (matrix)
//		 */
//		MPI_Status status;
//		MPI_File fh;
//		MPI_Offset displacement;
//
//		// Open the file
//		if (MPI_File_open(mpi_comm, filenameArray, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh)!=0) {
//				// Note: MPI_MODE_RDONLY = read only
//			cerr<<"ERROR! IkeWithPsALC_AD::readMatrixBinary(): Cannot open "<<filenameArray<<endl;
//			delete [] filenameArray;
//			delete [] localNnzArray_fromFile;
//			throw(REMBRANDT_ERROR_CODE);
//		}
//
//		for (set<int>::iterator it=Set_mpi_rank_inFile.begin(); it!=Set_mpi_rank_inFile.end(); ++it) {
//			int mpi_rank_file = *it;
//			int myNnz_fromFile = localNnzArray_fromFile[mpi_rank_file];
//
//			// 2-1. Read the row indices
//			int myOffsetNnz = 0;
//			for(int i=0; i<mpi_rank_file; ++i)
//				myOffsetNnz += localNnzArray_fromFile[i];
//			displacement = MPI_Offset(initDisp + sizeof(int)*myOffsetNnz);
//			if(MPI_File_set_view(fh, displacement, MPI_INT, MPI_INT, "native", MPI_INFO_NULL)!=0) {
//				// Note: native     = Data in this representation are stored in a file exactly as it is in memory
//				//       Internal   = Data in this representation can be used for I/O operations in a homogeneous or heterogeneous environment
//				//       External32 = This data representation states that read and write operations convert all data from and to the "external32" representation
//				cerr<<"ERROR! "<<classID<<"::readMatrixBinaryForRembrandtParallel(): Cannot set the first MPI_File_set_view -- offset="<<displacement<<endl;
//				delete [] filenameArray;
//				delete [] localNnzArray_fromFile;
//				throw(REMBRANDT_ERROR_CODE);
//			}
//
//			int *bufferRind = new int [myNnz_fromFile];
//			MPI_File_read(fh, bufferRind, myNnz_fromFile, MPI_INT, &status);
//
//			// 2-2. Read the column indices
//			displacement += MPI_Offset(nnz_fromFile*sizeof(int));
//			if(MPI_File_set_view(fh, displacement, MPI_INT, MPI_INT, "native", MPI_INFO_NULL)!=0) {
//				cerr<<"ERROR! "<<classID<<"::readMatrixBinaryForRembrandtParallel(): Cannot set the first MPI_File_set_view -- offset="<<displacement<<endl;
//				delete [] filenameArray;
//				delete [] localNnzArray_fromFile;
//				delete [] bufferRind;
//				throw(REMBRANDT_ERROR_CODE);
//			}
//
//			int *bufferCind = new int [myNnz_fromFile];
//			MPI_File_read(fh, bufferCind, myNnz_fromFile, MPI_INT, &status);
//
//			// 2-3. Read the values
//			displacement += MPI_Offset((nnz_fromFile-myOffsetNnz)*sizeof(int) + myOffsetNnz*sizeof(double));
//			if(MPI_File_set_view(fh, displacement, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL)!=0) {
//				cerr<<"ERROR! "<<classID<<"::readMatrixBinaryForRembrandtParallel(): Cannot set the first MPI_File_set_view -- offset="<<displacement<<endl;
//				delete [] filenameArray;
//				delete [] localNnzArray_fromFile;
//				delete [] bufferRind;
//				delete [] bufferCind;
//				throw(REMBRANDT_ERROR_CODE);
//			}
//
//			double *bufferValues = new double [myNnz_fromFile];
//			MPI_File_read(fh, bufferValues, myNnz_fromFile, MPI_DOUBLE, &status);
//
//			delete [] bufferRind;
//			delete [] bufferCind;
//			delete [] bufferValues;
//
//		}
//
//		/*
//		 * 3. Read the foot
//		 */
//		displacement = (MPI_Offset) (initDisp + sizeof(int)*nnz_fromFile*2 + + sizeof(double)*nnz_fromFile);
//		if(MPI_File_set_view(fh, displacement, MPI_INT, MPI_INT, "native", MPI_INFO_NULL)!=0) {
//			cerr<<"ERROR! "<<classID<<"::readMatrixBinaryForRembrandtParallel(): Cannot set the first MPI_File_set_view -- offset="<<displacement<<endl;
//			delete [] filenameArray;
//			delete [] localNnzArray_fromFile;
//			throw(REMBRANDT_ERROR_CODE);
//		}
//
//		int dummyInt;
//		MPI_File_read(fh, &dummyInt, 1, MPI_INT, &status);
//		if(dummyInt != EOF_ERROR_CHECK_CODE) {
//			cerr<<"ERROR! "<<classID<<"::readMatrixBinaryForRembrandtParallel(): The foot of the file (="<<dummyInt<<") is not equal to "<<EOF_ERROR_CHECK_CODE<<endl;
//			delete [] filenameArray;
//			delete [] localNnzArray_fromFile;
//			throw(REMBRANDT_ERROR_CODE);
//		}
//
//		/*
//		 * Close the file and clear some arrays
//		 */
//		MPI_File_close(&fh);
//
//		rind.resize(myNnz_fromFile);
//		cind.resize(myNnz_fromFile);
//		values.resize(myNnz_fromFile);
//		for(int i=0; i<myNnz_fromFile; ++i) {
//			rind[i] = bufferRind[i];
//			cind[i] = bufferCind[i];
//			values[i] = bufferValues[i];
//		}
//
//		delete [] filenameArray;
//		delete [] localNnzArray_fromFile;
//
//		MPI_Barrier(mpi_comm);
//	}

/*
 * Method: calcMatrixFromField
 * ---------------------------
 * Wrapper function for runEigenAnalysis
 */
void RembrandtWithModels::calcMatrixFromField() {
	int nVars = 5 + nScal;

	// +++++++++++++++++++++++++++++++++++++++++
	// Initialize the flow field
	// +++++++++++++++++++++++++++++++++++++++++
	// Update the data field if the Q1 file is specified
	string Q1filename = getStringParam("Q1_FILENAME", "DEFAULT_FILENAME_NONE");

	if(Q1filename.compare("DEFAULT_FILENAME_NONE") != 0 ) {
		double *qVec   = new double[ncv*nVars];

		double *lambda0Temp = NULL;
		if(NcontrolEqns > 0)
			lambda0Temp = new double [NcontrolEqns];

		double arclengthTemp;

		if(mpi_rank==0) cout<<">> calcMatrixFromField(): Initialize flow field from = "<<Q1filename<<endl;

		// You must update xcvMin and xcvMax before calling readPsALCdumpedDataSerial()
		for(int i=0; i<3; ++i) {
			xcvMin[i] =  ABSURDLY_BIG_NUMBER;
			xcvMax[i] = -ABSURDLY_BIG_NUMBER;
		}
		for(int icv=0; icv<ncv; ++icv) {
			for(int i=0; i<3; ++i) {
				if(x_cv[icv][i] < xcvMin[i])
					xcvMin[i] = x_cv[icv][i];
				if(x_cv[icv][i] > xcvMax[i])
					xcvMax[i] = x_cv[icv][i];
			}
		}
		for(int i=0; i<3; ++i) { // Add or subtract the grid tolerance for safety
			xcvMin[i] -= GRID_TOL;
			xcvMax[i] += GRID_TOL;
		}

		// Read the binary file
		readPsALCdumpedDataSerial(Q1filename, qVec, step, lambda0Temp, lambda, NcontrolEqns, arclengthTemp);
		if(mpi_rank==0) {
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
				printf("  lambda%d=%.3e", iEqn, lambda[iEqn]);
			printf("  arclength=%.4e \n", arclengthTemp);
		}

		// Update the flow by the data in the binary file
		for(int icv=0; icv<ncv; ++icv) {
			int tempInt = icv*nVars;

			rho[icv] = qVec[tempInt];
			for(int i=0; i<3; ++i)
				rhou[icv][i] = qVec[tempInt+1+i];
			rhoE[icv] = qVec[tempInt+4];

			if(nScal > 0) {
				for(int iScal=0; iScal<nScal; ++iScal) {
					double *scalArray = scalarTranspEqVector[iScal].phi;
					scalArray[icv] = qVec[tempInt+5+iScal];
				}
			}
		}

		// Clear memory
		delete [] qVec;
		delete [] lambda0Temp;

		MPI_Barrier(mpi_comm);
	} else {
		if(mpi_rank==0) cout<<">> calcMatrixFromField(): Initialize flow from the restart file"<<endl;
	}

	// +++++++++++++++++++++++++++++++++++++++++
	// Allocate memory related to RHS evaluation
	// +++++++++++++++++++++++++++++++++++++++++
	double *rhsVec = new double[ncv*nVars];

	if(weightRhs == NULL)
		weightRhs = new double [ncv*nVars];
	for(int i=0; i<ncv*nVars; ++i)
		weightRhs[i] = 1.0;

	assert(RHSrhoScal == NULL);
	if (nScal > 0)
		getMem2D(&RHSrhoScal,  0, nScal-1, 0, ncv_g-1, "rhsScal");

	// The followings are also required during the RHS evaluation
	if(nScal==1) {
		assert(RHSsa == NULL);
		registerScalar(RHSsa, "RHSSA", CV_DATA);
	} else if(nScal==2) {
		assert(RHSkine == NULL && RHSomega == NULL);
		registerScalar(RHSkine,  "RHSKINE",  CV_DATA);
		registerScalar(RHSomega, "RHSOMEGA", CV_DATA);
	} else if(nScal==3) {
		assert(RHSZMean == NULL && RHSZVar == NULL && RHSCMean == NULL);
		registerScalar(RHSZMean, "RHSZMEAN", CV_DATA);
		registerScalar(RHSZVar,  "RHSZVAR",  CV_DATA);
		registerScalar(RHSCMean, "RHSCMEAN", CV_DATA);
	} else if(nScal==5) {
		assert(RHSkine == NULL && RHSomega == NULL && RHSZMean == NULL && RHSZVar == NULL && RHSCMean == NULL);
		registerScalar(RHSkine,  "RHSKINE",  CV_DATA);
		registerScalar(RHSomega, "RHSOMEGA", CV_DATA);
		registerScalar(RHSZMean, "RHSZMEAN", CV_DATA);
		registerScalar(RHSZVar,  "RHSZVAR",  CV_DATA);
		registerScalar(RHSCMean, "RHSCMEAN", CV_DATA);
	}

	// +++++++++++++++++++++++++++++++++++++++++
	// Perform the basic setting
	// +++++++++++++++++++++++++++++++++++++++++
	calcStateVariables();     // Compute velocity, pressure, temperature, enthalpy and speed of sound at cell centers
//		calcMaterialProperties(); // Compute viscosity and thermal diffusivity at cell faces
	setBC();
	for(int ifa=nfa; ifa<nfa_b2; ++ifa)
		setBC1D(ifa, UgpWithCvCompFlow::rho, UgpWithCvCompFlow::rhou, UgpWithCvCompFlow::rhoE); // Note: ncvggf's cannot be not updated, so the user should take care of it
	calcRansTurbViscMuet();

	// +++++++++++++++++++++++++++++++++++++++++
	// Compatibility check if required
	// +++++++++++++++++++++++++++++++++++++++++
	string temp = getStringParam("SKIP_COMPATIBILITY_CHECK", "YES");
	if(temp=="YES" || temp=="yes" || temp=="Y" || temp=="y") {
		if(mpi_rank==0)
			cout<<endl
				<<"> CAUTION!! SKIP_COMPATIBILITY_CHECK = "<<temp<<": skip the compatibility-check between IKE and JOE"<<endl;
	} else {
		compCheckOn = true;  // We want to allow lambda_AD only during the compatibility check

		bool negligibleError = compatibilityCheck();

		if(!negligibleError) {
			if(mpi_rank==0)
				cout<<"ERROR in "<<classID<<"::run(): compatibility check by compatibilityCheck() FAILs"<<endl;
			throw(REMBRANDT_ERROR_CODE);
		}

		compCheckOn = false;  // We want to allow lambda_AD only during the compatibility check
		MPI_Barrier(mpi_comm);
	}

	// +++++++++++++++++++++++++++++++++++++++++
	// Calculate the RHS and the Jacobian matrix
	// +++++++++++++++++++++++++++++++++++++++++
	string booleanString = getStringParam("PERFORM_NEWTON_METHOD_ONCE", "FALSE");
	if(mpi_rank == 0)
		cout << "> PERFORM_NEWTON_METHOD_ONCE = " << booleanString << endl
		     << endl;
	std::transform(booleanString.begin(), booleanString.end(), booleanString.begin(), ::tolower);

	if(booleanString.compare("true")==0 || booleanString.compare("yes")==0) {
		if(mpi_rank==0)
			cerr<<"ERROR! PERFORM_NEWTON_METHOD_ONCE has not been developed yet"<<endl;
		// Note: Currently the calcSteadySolnByNewtonRembrandt() method has an error --
		//       The PETSc package is frozen while solving the linear system using KSP.
		throw(REMBRANDT_ERROR_CODE);

//		if (!checkParam("NEWTON_PARAMETERS_2NDPT")) {
//			ParamMap::add("NEWTON_PARAMETERS_2NDPT  MAX_ITER=100  ABS_RESID=1.0e-14  REL_RESID=1.0e-8"); // add default values
//			if (mpi_rank == 0)
//				cout<< "WARNING: added keyword \"NEWTON_PARAMETERS_2NDPT  MAX_ITER=100  ABS_RESID=1.0e-14  REL_RESID=1.0e-8\""<< " to parameter map!" << endl;
//		}
//		int maxIterNewton   = getParam("NEWTON_PARAMETERS_2NDPT")->getInt("MAX_ITER");
//		double absTolNewton = getParam("NEWTON_PARAMETERS_2NDPT")->getDouble("ABS_RESID");
//		double relTolNewton = getParam("NEWTON_PARAMETERS_2NDPT")->getDouble("REL_RESID");
//
//		double *qVec = new double[ncv*(5+nScal)];
//
//		// Since IkeWithPsALC_AD::getSteadySolnByNewton() deletes jacMatrix in the routine,
//		// you should call calcSteadySolnByNewtonRembrandt() instead.
//		calcSteadySolnByNewtonRembrandt(qVec, rhsVec, maxIterNewton, absTolNewton, relTolNewton, nScal);
//
//		delete [] qVec;
	} else {
		calcRhsWithJacMatrixRembrandt<MatComprsedSTL>(jacMatrixSTL, rhsVec, nScal);
	}

	// +++++++++++++++++++++++++++++++++++++++++
	// Clear memory
	// +++++++++++++++++++++++++++++++++++++++++
	if(rhsVec != NULL) {
		delete [] rhsVec; 	rhsVec = NULL;
	}
	if(weightRhs != NULL) {
		delete [] weightRhs;  	weightRhs = NULL;
	}
	if(RHSrhoScal != NULL && nScal > 0)
		freeMem2D(RHSrhoScal, 0, nScal-1, 0, ncv-1);
}

/*
 * Method: calcRhsWithJacMatrixRembrandt
 * -------------------------------------
 * Calculate the Jacobian matrix as IKE
 */
template <class MatT>
void RembrandtWithModels::calcRhsWithJacMatrixRembrandt(MatT &jacMatrix, double* rhsVec, const int nScal) {
	int nVars = 5 + nScal;
	assert(rhsVec != NULL);

	if(NcontrolEqns > 0) {
		if(lambda == NULL) {
			if(mpi_rank == 0)
				cout<<"WARNING in "<<classID<<"::calcRhsWithJacMatrixRembrandt(): lambda == NULL. Assign memory."<<endl
				    <<endl;
			lambda = new double [NcontrolEqns];
		}
	}

	// Calculate the Jacobian matrix
	if(mpi_rank==0) cout<<">> calcRhsWithJacMatrixRembrandt(): Calculate Jacobian: "<<endl;

	calcJacobian1DAD(jacMatrixSTL, rhsVec, debugLevel, 0); // Note: arguments = matrix, rhs, debugLevel, number of AD system parameters

	// Calculate the norm of RHS and show it on the screen
	int whichNorm = 1; // one-norm
	if(mpi_rank==0)
		cout << endl
			 << "> NORM = " << whichNorm << "-norm" << endl;

	double* residNormVecFlow = new double [5+nScal];
	calcResidualsFrom1Drhs(residNormVecFlow, rhsVec, whichNorm);
	double residNormTotFlowOnly = calcSumResidual(residNormVecFlow, whichNorm);

	if(mpi_rank == 0) {
		cout << "> Total residual of the flow RHS = " << residNormTotFlowOnly << endl;

		printf("         rho          rhou-X       rhou-Y       rhou-Z       rhoE      ");
		for (int iScal = 0; iScal < nScal; iScal++)
			printf("%12s", scalarTranspEqVector[iScal].getName());
		cout << endl;

		printf("RESID: %12.4e %12.4e %12.4e %12.4e %12.4e", residNormVecFlow[0], residNormVecFlow[1], residNormVecFlow[2], residNormVecFlow[3], residNormVecFlow[4]);
		for (int iScal = 0; iScal < nScal; iScal++)
			printf("%12.4e", residNormVecFlow[5+iScal]);
		cout << endl;

		if(residNormTotFlowOnly > A_LARGE_TOT_FLOW_RESID)
			cout << endl
			     << "WARNING! The total residual (=" << residNormTotFlowOnly << ") seems to be too large!" << endl
			     << "  -- You may want to check if the residual of the original IKE gives the same order of residual." <<endl
			     << "     If not, please run IKE with the same number of CPU cores to produce new restart and Q1 files." << endl
			     << endl;
	}

	// Store the RHS errors
	if(mpi_rank == 0) cout<<"> fabs(RHS error) is stored in RHSRHO, RHSRHOU, RHSRHOE";

	assert(RHSrho  != NULL);
	assert(RHSrhou != NULL);
	assert(RHSrhoE != NULL);
	if(nScal == 2) {
		if(mpi_rank == 0) cout<<", RHSKINE, and RHSOMEGA"<<endl
							  <<endl;

		assert(RHSkine  != NULL);
		assert(RHSomega != NULL);
	}
	for(int icv=0; icv<ncv; ++icv) {
		int indexStart = icv*nVars;
		RHSrho[icv] = fabs(rhsVec[indexStart]);
		for(int i=0; i<3; ++i)
			RHSrhou[icv][i] = fabs(rhsVec[indexStart + 1+i]);
		RHSrhoE[icv] = fabs(rhsVec[indexStart + 4]);

		if(nScal == 2) {
			RHSkine[icv]  = fabs(rhsVec[indexStart + 5]);
			RHSomega[icv] = fabs(rhsVec[indexStart + 6]);
		}
	}

	string tempString = getStringParam("DUMP_FLOWFIELD_AFTER_RHS_CALC", "NO");
	std::transform(tempString.begin(), tempString.end(), tempString.begin(), ::tolower);
	if(tempString.compare("yes")==0 || tempString.compare("y")==0)
		writeData(0, 0);

	// Clear memory
	delete [] residNormVecFlow;
}

/*
 * Method: calcRhsWithJacMatrixRembrandt
 * -------------------------------------
 * Original code = IkeWithPsALC_AD::getSteadySolnByNewton
 * Since IkeWithPsALC_AD::getSteadySolnByNewton() deletes jacMatrix in the routine,
 * this function must be called in order to get the jacMatrix.
 */
void RembrandtWithModels::calcSteadySolnByNewtonRembrandt(double* q, double* rhs,
		const int maxIterNewton, const double absTolNewton, const double relTolNewton, const int nScal) {
	string funcID = "calcSteadySolnByNewtonRembrandt";

	/***********************/
	/*  Initialize lambda  *
	 ***********************/
	if(this->NcontrolEqns > 0) {
		if(lambda == NULL) {
			if(mpi_rank == 0)
				cout<<"WARNING in "<<classID<<"::"<<funcID<<"(): lambda == NULL. Assign memory."<<endl
				    <<endl;
			lambda = new double [this->NcontrolEqns];
		}
	} else {
		if(mpi_rank == 0)
			cout<<"WARNING in "<<classID<<"::"<<funcID<<"(): NcontrolEqns == 0."<<endl
			    <<endl;
	}

	/**************************************************/
	/*  Variables that must be passed to this method  *
	 *  but not required if NcontrolEqns = 0.         *
	 **************************************************/
	const bool writeJacOnFile = false;
	const int NcontrolEqns = 0;  // This variable is zero only in this method!
	const double *q1 = NULL;
	double** q_tangent      = NULL;
	double*  lambda_tangent = NULL;
	const double* lambda0 = NULL;
	const double* lambda1 = NULL;
	const double weightLambda = 0.0;
	const double arcLength = 0.0;

	/****************************************************************************/
	/*  Start the contents in IkeWithPsALC_AD::getSteadySolnByNewton from here  *
	 ****************************************************************************/
	static bool firstCall = true;

	assert(q!=NULL && rhs!=NULL);
	if(NcontrolEqns==0)
		assert(q1==NULL && q_tangent==NULL && lambda_tangent==NULL && lambda1==NULL);
	else
		assert(q1!=NULL && q_tangent!=NULL && lambda_tangent!=NULL && lambda1!=NULL && arcLength!=0.0);

	int m = nScal+5; // number of variables (rho, rhou, rhov, rhow, rhoE, and scalars)

	int countNegative = 0; // number of negative rho or press. If this is greater than 0, the backtracking line search becomes active

	/******************************/
	/*  Parameter initialization  *
	 ******************************/
	/* Read parameters from input file */
	// Read basic Newton parameters
	getNewtonParam();

	// How to calculate the Jacobian matrix (ROW_1D, ORDINARY_2D): For a better maintenance of the code, call a separated function
	HOW_TO_CALC_JAC howToCalcJac = getHowToCalcJac();

	// Modified Newton method
	MODIFIED_NEWTON modifiedNewtonMethod; // "MODIFIED_NEWTON_NS"
	double modifNewtonRelResid; // "MODIFIED_NEWTON_NS_PARAMETERS"-->"REL_RESID": This value will not be used if modifiedNewtonMethod==BASIC
	int modifNewtonFreq;        // "MODIFIED_NEWTON_NS_PARAMETERS"-->"FREQ": This value will not be used if modifiedNewtonMethod==BASIC
	getParamsModifiedNewtonMethod(modifiedNewtonMethod, modifNewtonRelResid, modifNewtonFreq);

	// RHS weight
	string tempString = getStringParam("WEIGHT_RHS_METHOD", "NO_RHSWEIGHT");
	if(tempString.compare("REF_VALUES")==0)                      weightRhsMethod = REF_VALUES;
	else if(tempString.compare("LOCAL_VALUES")==0)               weightRhsMethod = LOCAL_VALUES;
	else if(tempString.compare("REF_VALUES_AND_DT_OVER_VOL")==0) weightRhsMethod = REF_VALUES_AND_DT_OVER_VOL;
	else                                                         weightRhsMethod = NO_RHSWEIGHT;
	if(firstCall && mpi_rank==0)
		cout<<"> WEIGHT_RHS_METHOD == "<<weightRhsMethod<<" ("<<tempString<<")"<<endl;

	// A-priori regularization
	//      Function: If x<=x_2, w = 10^(-digit) + (1.0-10^(-digit)) * exp(-beta*(x-x2)^2)
	//                Otherwise, w = 1.0
	//      Note: 1. The function is uni-directional Gaussian.
	//            2. Maximum function value is 1.0, and Minimum function value is 10^(-digit)
	//      If beta = digit*ln(10) / (x_2-x_1)^2, w(x_1) = 2*10^(-digit)
	//                                            w(x_2) = 1.0
	if (!checkParam("A_PRIORI_WEIGHT_REGULARIZ")) {
		ParamMap::add("A_PRIORI_WEIGHT_REGULARIZ  USE=NO  DIGIT=3.0  X1=0.35  X2=0.4"); // add default values
		if (mpi_rank == 0)
			cout << "WARNING: added keyword \"A_PRIORI_WEIGHT_REGULARIZ  USE=NO  DIGIT=3.0  X1=0.35  X2=0.4\" to parameter map!"<<endl;
	}

	tempString = getParam("A_PRIORI_WEIGHT_REGULARIZ")->getString("USE");
	bool   useAprioriWeightFunction = false;
	if(tempString.compare("YES")==0 || tempString.compare("yes")==0 || tempString.compare("Y")==0 || tempString.compare("y")==0)
		useAprioriWeightFunction = true;

	double digitAprioriWeight = getParam("A_PRIORI_WEIGHT_REGULARIZ")->getDouble("DIGIT"); // If 3.0, Want to have minimum weight of 1.0e-3
	double x1AprioriWeight    = getParam("A_PRIORI_WEIGHT_REGULARIZ")->getDouble("X1"); ;  // Inlet location of Scramjet
	double x2AprioriWeight    = getParam("A_PRIORI_WEIGHT_REGULARIZ")->getDouble("X2"); ;  // Slightly upstream of combustion start
	double betaAprioriWeight  = digitAprioriWeight*2.302585093/pow(x2AprioriWeight - x1AprioriWeight, 2.0); // Note: ln(10) = 2.302585093

	// Run the iterations for few more steps even after they converges
	int moreNTstepsBelowConv_moreSteps;
	double moreNTstepsBelowConv_resid, moreNTstepsBelowConv_deltaQ;
	getParamMoreNewtonSteps(moreNTstepsBelowConv_moreSteps, moreNTstepsBelowConv_resid, moreNTstepsBelowConv_deltaQ, firstCall);

	// Note: q has been initialized before calling this method

	/* Initialize rhs */
	for(int i=0; i<ncv*m; ++i)
		rhs[i] = 0.0;
	if(mpi_rank==mpi_size-1)
		for(int iParam=0; iParam<NcontrolEqns; ++iParam)
			rhs[ncv*(5+nScal)+iParam] = 0.0;

	// Set iterNewton to zero (Some messages will be printed in some methods if iterNewton == 0. Also the barrier works correctly for the first time only if iterNewton is reset)
	iterNewton = 0;

	int startingNewtonIter = 0;
	if(initThirdFromQ1) {
		if(step-startingStep==2) {
			startingNewtonIter = getIntParam("INIT_THIRD_NEWTON_ITERS", "0");
			iterNewton += startingNewtonIter;

			MPI_Barrier(mpi_comm);

			if(firstCall && mpi_rank==0)
				cout<<"> For the THIRD point, iterNewton starts from "<<startingNewtonIter<<endl;
		}
	}

	// InitialHookNewton() & initialHookNewton_firstCall()
	initialHookNewton();
	if(firstCall)
		initialHookNewton_firstCall();

	/*********************************************************/
	/*  Check if RHS's by ADOLC and normal JOE are the same  *
	 *********************************************************/
	if(firstCall) {
		string temp = getStringParam("SKIP_COMPATIBILITY_CHECK", "NO");
		if(temp=="YES" || temp=="yes" || temp=="Y" || temp=="y") {
			if(mpi_rank==0)
				cout<<endl
				    <<"> CAUTION!! SKIP_COMPATIBILITY_CHECK = "<<temp<<": skip the compatibility-check between IKE and JOE"<<endl;
		} else {
			bool negligibleError = true;
			if(firstCall)
				negligibleError = compatibilityCheck();
			if(!negligibleError) {
				if(mpi_rank==0)
					cout<<"ERROR in "<<classID<<"::"<<funcID<<"(): compatibility check by compatibilityCheck() FAILs"<<endl;
				throw(PSALC_ERROR_CODE);
			}
		}
	}

	/******************************************/
	/*  Newton iteration:                     *
	 *    Solve Q_{n+1} = Q_{n} - L^{-1}*RES  *
	 ******************************************/
	if(NcontrolEqns > 0 && mpi_rank == mpi_size-1) {
		if(fabs(lambda_tangent[0]) < 1.0e-6)
			cout<<"> CAUTION!! lambda_tangent (="<<lambda_tangent[0]<<") is TOO small, which can cause a convergence problem"<<endl;
	}

	// Settings for the modified Newton's method
	bool useExactJac = true, useExactJacTempo = true;
	int countModifiedNewton = 0;
	int myNnzJac = 0, nnzJac = 0;

	// Allocate and initialize q_tangent_inJacMat and lambda_tangent_inJacMat.
	//   Note: In the while loop below, they will be initialized again in every iteration because they can be modified to stabilize the Newton's method.
	double **q_tangent_inJacMat     = NULL; // Note: q_tangent may need to be modified: minus sign, stabilization for singular Jacobian, ...
	double *lambda_tangent_inJacMat = NULL; // Note: lambda_tangent may need to be modified: minus sign, stabilization for singular Jacobian, ...
	if(NcontrolEqns > 0) {
		getMem2D(&q_tangent_inJacMat, 0, NcontrolEqns-1, 0, (ncv*m)-1, "q_tangent_inJacMat");
		// Note: The size of q_tangent is the same for both ORDINARY_2D and ROW_1D: In the matrix used for Newton's method, you don't need to have ncv_gg even for ORDINARY_2D

		lambda_tangent_inJacMat = new double [NcontrolEqns];
	}

	// Allocate and initialize phi
	double *phi; // note: this array will be used while solving the linear system
	if(mpi_rank==mpi_size-1)
		phi = new double [ncv*m + NcontrolEqns];
	else
		phi = new double [ncv*m];

	if(q1 != NULL) { // Actually this step is not important
		for(int i=0; i<ncv*(5+nScal) + ((mpi_rank==mpi_size-1) ? NcontrolEqns : 0); ++i)
			phi[i] = q[i] - q1[i];
	} else {
		for(int i=0; i<ncv*(5+nScal) + ((mpi_rank==mpi_size-1) ? NcontrolEqns : 0); ++i)
			phi[i] = 0.0;
	}

	// Residuals
	int whichNorm = 1; // one-norm
	if(firstCall && mpi_rank==0)
		cout<<endl
		    <<"> NORM = "<<whichNorm<<"-norm"<<endl;

	double* residNormVecFlow    = new double [5+nScal];
	double* residNormVecFlowOld = new double [5+nScal];
	double residNormTotFlowOnly = ABSURDLY_BIG_NUMBER;
	residNormTot                = ABSURDLY_BIG_NUMBER;
	residNormTotOld             = ABSURDLY_BIG_NUMBER;

	// The residuals of the tangential conditions (For the details, google "pseudo-arclength continuation method")
	double* Nres = NULL; // Note: 1. Nres = q_tangent*(q_guess-q1) + lambda_tangent*(lambda_guess-lambda1) - ds
	                     //       2. Nres must be zero if the tangential condition is satisfied.
	if(NcontrolEqns>0) {
		Nres = new double [NcontrolEqns];
		for(int iParam=0; iParam<NcontrolEqns>0; ++iParam)
			Nres[iParam] = -ABSURDLY_BIG_NUMBER; // Initialize with a very small number for a debugging purpose
	}

	/*
	 * Calculate RHS and Jacobian based on the initial guess before the Newton iterations
	 */
	try {
		if(howToCalcJac == ROW_1D) {
			countNegative = calcJacobian1DAD(jacMatrixSTL, rhs, debugLevel, NcontrolEqns);
			myNnzJac      = jacMatrixSTL.get_nnz();
		} else if(howToCalcJac == ORDINARY_2D) {
			countNegative = calcJacobianAD(jacMatrix, rhs, debugLevel, NcontrolEqns);
			myNnzJac      = jacMatrix.get_nnz();
		} else { // ERROR
			freeMemGetSteadySolnByNewton(phi, residNormVecFlow, residNormVecFlowOld, Nres, jacMatrix, jacMatrixSTL);
			assert(true);
		}
		MPI_Allreduce(&myNnzJac, &nnzJac, 1, MPI_INT, MPI_SUM, mpi_comm);

		calcResidualsFrom1Drhs(residNormVecFlow, rhs, whichNorm);
		residNormTotFlowOnly = calcSumResidual(residNormVecFlow, whichNorm);
	}
	catch(int e) { // Catch and re-throw
		freeMemGetSteadySolnByNewton(phi, residNormVecFlow, residNormVecFlowOld, Nres, jacMatrix, jacMatrixSTL);
		throw(e);
	}
	catch(...) { // Catch and re-throw
		freeMemGetSteadySolnByNewton(phi, residNormVecFlow, residNormVecFlowOld, Nres, jacMatrix, jacMatrixSTL);
		throw;
	}

	// Fill up the last element of rhs ( Nres == ds - q_tangent*(q_guess-q1) - weightLambda*lambda_tangent*(lambda_guess-lambda1) ) for mpi_rank == mpi_size-1
	// Note: 1. Even though this part of the code is designed to handle multiple parameters (or lambdas), the actual formulation here is only for one parameters.
	//       2. You don't have to pass weighted_q_tangent instead of q_tangent since it will be weighted in getFlowWeightsForInnerProduct() called by vecDotVecWithWeight().
//	if(NcontrolEqns>0) {
//		calcNres(Nres, q, rhs, m, NcontrolEqns, q1, q_tangent, lambda_tangent, lambda1, weightLambda, arcLength);
//		MPI_Barrier(mpi_comm);
//
//		if(mpi_rank==0) {
//			for(int iParam=0; iParam<NcontrolEqns; ++iParam) {
//				if(fabs(Nres[iParam]/arcLength) > 1.0e4*double(NcontrolEqns*cvora[mpi_size])*MACHINE_EPS) { // Note: Nres must be very close to zero
//					cout<<"WARNING in "<<classID<<"::"<<funcID<<"(): iterNewton=="<<iterNewton<<", The tangential condition["<<iParam<<"] shows a big residual."<<endl
//					    <<"                                    Tangential-residual = "<<Nres[iParam]<<" (target arclength="<<arcLength<<")"<<endl;
//				}
//				if(isNaN(Nres[iParam])) {
//					cout<<"ERROR in "<<classID<<"::"<<funcID<<"(): iterNewton=="<<iterNewton<<", The residual of the tangential condition becomes NaN for iParam="<<iParam<<endl;
//					throw(PSALC_ERROR_CODE);
//				}
//			}
//		}
//
//		if(mpi_rank == mpi_size-1) {
//			for(int iParam=0; iParam<NcontrolEqns; ++iParam)
//				rhs[ncv*m+iParam] = Nres[iParam];
//		}
//	}

	// Update residNormTot with Nres and Show the total residual calculated with the flow variables
	if(NcontrolEqns>0) {
		residNormTot = updateTotResidualWithNres(residNormTotFlowOnly, Nres, NcontrolEqns, whichNorm);

		if(debugLevel>0 && mpi_rank==0)
			printf("           >> Residual: flow RHS = %.5e, tangential RHS[0] = %.5e --> Total residual = %.5e\n", residNormTotFlowOnly, fabs(Nres[0]), residNormTot);
	} else {
		residNormTot = residNormTotFlowOnly;

		if(debugLevel>0 && mpi_rank==0)
			printf("           >> Residual of the flow RHS = %.5e\n", residNormTot);
	}

	// Update residNormOld and residNormTotOld (It will be used for backtracking)
	for(int i=0; i<5+nScal; ++i)
		residNormVecFlowOld[i] = residNormVecFlow[i];
	residNormTotOld = residNormTot;

	// Calculate deltaQ from the previous point (q1) to the initial guess (q)
	double deltaQnorm = 0.0;
//	if(NcontrolEqns > 0) { // If NcontrolEqns==0, q1 was not passed (i.e. q1 is NULL if NcontrolEqns==0)
//        double mySumPhiSq = 0.0, sumPhiSq;
//
//		for(int i=0; i<ncv*m; ++i)
//				mySumPhiSq += pow(q[i]-q1[i], 2.0);
//
//		if(mpi_rank==mpi_size-1 && NcontrolEqns>0)
//			for(int i=0; i<NcontrolEqns; ++i)
//				mySumPhiSq += weightLambda * pow(q[ncv*(5+nScal)+i]-q1[ncv*(5+nScal)+i], 2.0);
//
//		MPI_Allreduce(&mySumPhiSq, &sumPhiSq, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
//
//		deltaQnorm = sqrt(sumPhiSq);
//	}
//
//	if(newtonParam.trustRegionSize > 8.0*deltaQnorm && NcontrolEqns > 0) {
//		if(mpi_rank==0)
//			cout<<"WARNING in "<<classID<<"::"<<funcID<<"(): TRUST_REGION_SIZE="<<newtonParam.trustRegionSize<<" is larger than 8 * DELTA_Q of the init. guess="<<deltaQnorm<<endl;
//		newtonParam.trustRegionSize = 8.0*deltaQnorm;
//		if(mpi_rank==0)
//			cout<<"                                    Reduce TRUST_REGION_SIZE to "<<newtonParam.trustRegionSize<<endl
//			    <<endl;
//	}

	// Show the result on the screen
	if(isNaN(residNormTot)) {
		if(mpi_rank==0) cerr<<classID<<"::"<<funcID<<"(): The residual of the initial field is NaN"<<endl;
		if(debugLevel>1)
			for(int i=0; i<ncv*(5+nScal); ++i)
				if(isNaN(rhs[i]))
					printf("     mpi_rank=%d:  rhs[%d]=NaN \n", mpi_rank, i);
		throw(PSALC_ERROR_CODE);
	} else if(debugLevel > 0) {
		if(mpi_rank==0) {
			// Header
			printf("\n");
			printf("        ITER         RESIDUAL     DELTA_Q       NNZ_JAC  LS_ITER       LS_RESID  NEWTON(JAC/RESID)");
			if(NcontrolEqns > 0)
				printf("          LAMBDA\n");
			else
				printf("\n");

			// Residual information
			if(NcontrolEqns > 0)
				printf("   InitGuess  %15.8e  %10.3e  %12d                                            %15.8e\n", residNormTot, deltaQnorm, nnzJac, lambda[0]);
			else
				printf("   InitGuess  %15.8e              %12d\n", residNormTot, nnzJac);
		}
	}

	// Also save the result on a file
	if(mpi_rank==0) {
		FILE *fp;
		if(firstCall) {
			fp = fopen(NEWTON_STATUS_FILENAME, "w");
			fprintf(fp, "STEP,  ITER,      JAC_NNZ,    RESID_INIT,  RELAXATION,   RESID_FINAL,      NRES,        LAMBDA \n");
		} else
			fp = fopen(NEWTON_STATUS_FILENAME, "a");

		fprintf(fp, "%4d, %5d,             ,              ,            , %13.6e,          , %13.6e\n", step, 0, residNormTot, lambda[0]);
		fclose(fp);
	}

	// Gives a warning message if the residual is too high even though the initial guess comes from a binary
	if(step - startingStep==2 && initThirdFromQ1) {
		double startingResidual = getDoubleParam("INIT_THIRD_RESID", "0.0");
		if(residNormTot > startingResidual)
			if(mpi_rank==0)
				cout<<endl
				    <<"WARNING in "<<classID<<"::"<<funcID<<"(): The initial residual is greater than INIT_THIRD_RESID="<<startingResidual<<endl
				    <<"                                    Possible reasons - 1. Check if REDUCE_ARCLENTH_THIRD is correct"<<endl
				    <<"                                                       2. Check if LAMBDA_INITIAL_0 and LAMBDA_INITIAL_1 are correct"<<endl
				    <<endl;
	}

	// Open a file to save the convergence history of the linear solver
	if(mpi_rank==0) {
		FILE *fp = fopen(LS_CONVERGENCE_FILENAME, "w");
		fprintf(fp, "STEP=%d\n", step);
		fprintf(fp, "NT_ITER,   ITER,   RESIDUAL\n");
		fclose(fp);
	}

	/*
	 * Calculate tolNewton
	 */
	initResidNorm = residNormTot;
	double tolNewton = absTolNewton + relTolNewton * initResidNorm;

	MPI_Barrier(mpi_comm);

	/*
	 * Iterations
	 */
	++iterNewton;
	int succeessiveBacktrack = 0; // Count the number of successive backtracking

	bool done = false;
	int  countConvergedSteps = 0; // In the case that the user wants to run the iterations even after the simulations converges.
	int  newtonDumpQ1Interval = getIntParam("NEWTON_DUMP_Q1_INTERVAL", "100");

	double relaxation = getDoubleParam("UNDER_RELAXATION", "1.0");

	double alpha_stab_newton = 1.0; // This will be only used for "CFL_BASED_DIAG" stabilized-Newton
	if(newtonParam.stabilizationMethodType != NO_STAB_NEWTON) // Apply diagonal scaled by CFL and local volume
		alpha_stab_newton = newtonParam.stabilizationAlpha;

	while(!done) {
		/******
		 ** Check iterations and report error
		 ******/
		if(iterNewton > maxIterNewton+startingNewtonIter) { // If the Newton procedure fails to converge (Note: startingNewtonIter=0 unless initial guess for the thrid point)
			if(mpi_rank==0) {
				cout<<"ERROR in "<<classID<<"::"<<funcID<<"(): Newton iter has NOT been converged until the "<<iterNewton-1<<"th iter"<<": residual="<<residNormTot;
				printf("  lambda[0]=%.15f\n", lambda[0]);
			}

			if(NcontrolEqns>0) {
				char filename[50];
				sprintf(filename, "Q1_PT%05d.dumped.%05d.bin", step, iterNewton-1);
				if(mpi_rank==0) {
					cout<<"                                  Dump the current Q vector on "<<filename;
					printf(", lambda[0]=%.15f\n", lambda[0]);
				}
				writePsALCdumpedDataParallel(filename, step, arcLength, lambda1, lambda, NcontrolEqns, q);
			}

			writeData(step, iterNewton-1);

			freeMemGetSteadySolnByNewton(phi, residNormVecFlow, residNormVecFlowOld, Nres, jacMatrix, jacMatrixSTL);
			throw (PSALC_ERROR_CODE);
		} else {
			if(iterNewton>startingNewtonIter+1 && newtonDumpQ1Interval>0) { // If the user wants to dump the data at every "newtonDumpQ1Interval" Newton step
				if((iterNewton-1)%newtonDumpQ1Interval == 0) {
					if(NcontrolEqns>0) {
						char filename[50];
						sprintf(filename, "Q1_PT%05d.dumped.%05d.bin", step, iterNewton-1);
						if(mpi_rank==0)
							cout<<"   > Dump the current Q vector on "<<filename<<": lambda[0]="<<std::setprecision(15)<<lambda[0]<<endl;
						writePsALCdumpedDataParallel(filename, step, arcLength, lambda1, lambda, NcontrolEqns, q);
					}

					writeData(step, iterNewton-1);
				}
			}
		}

		/******
		 ** Ramp maximum iterations and residual (For iterative solver)
		 ******/
		if(newtonParam.intervalDecLS > 0)
			if ((iterNewton > startingNewtonIter+1) && (iterNewton >= newtonParam.startDecLS)) {
				if((iterNewton-newtonParam.startDecLS) % newtonParam.intervalDecLS == 0) {
					if(newtonParam.maxIterLS < newtonParam.maxFinalIterLS)
						newtonParam.maxIterLS += newtonParam.incIterLS;
					if(newtonParam.zeroAbsLS > newtonParam.minZeroAbsLS+MACHINE_EPS || newtonParam.zeroRelLS > newtonParam.minZeroRelLS+MACHINE_EPS) {
						newtonParam.zeroAbsLS *= newtonParam.decZeroLS;
						newtonParam.zeroRelLS *= newtonParam.decZeroLS;
					}
					petscSolver2->setTresholds(newtonParam.zeroAbsLS, newtonParam.zeroRelLS, newtonParam.maxIterLS);

					if(debugLevel>1 && mpi_rank==0) {
						printf("RAMP_LSNT_PARAMETERS: MAX_ITER=%d  MIN_ABS_RESID=%.2e  MIN_REL_RESID=%.2e\n", newtonParam.maxIterLS, newtonParam.zeroAbsLS, newtonParam.zeroRelLS);
					}
				}
			}

//		/******
//		 ** Get q_tangent_inJacMat and lambda_tangent_inJacMat that will be actually used in the Jacobian matrix
//		 ******/
//		if(NcontrolEqns > 0) {
//			assert(q_tangent_inJacMat != NULL && q_tangent != NULL && lambda_tangent_inJacMat != NULL && lambda_tangent != NULL);
//
//			if(howToCalcJac == ROW_1D) {
//				for(int iParam=0; iParam<NcontrolEqns; ++iParam)
//					for(int i=0; i<ncv*m; ++i)
////						q_tangent_inJacMat[iParam][i] = -q_tangent[iParam][i];
//						q_tangent_inJacMat[iParam][i] = -q_tangent[iParam][i] * AreaLambda;
//				for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
////					lambda_tangent_inJacMat[iEqn] = -lambda_tangent[iEqn];
//					lambda_tangent_inJacMat[iEqn] = -lambda_tangent[iEqn] * AreaLambda;
//			} else if(howToCalcJac == ORDINARY_2D) { // Note: ORDINARY_2D is still using old formulation unlike ROW_1D
//				for(int iParam=0; iParam<NcontrolEqns; ++iParam)
//					for(int i=0; i<ncv*m; ++i)
//						q_tangent_inJacMat[iParam][i] = q_tangent[iParam][i];
//				for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
//					lambda_tangent_inJacMat[iEqn] = lambda_tangent[iEqn];
//			} else {
//				throw(PSALC_ERROR_CODE);
//			}
//		}

		/******
		 ** Stabilized Newton's method for ill-conditioned (or almost singular) Jacobians: This will modify the Jacobian matrix
		 ******/
		if(newtonParam.stabilizationMethodType != NO_STAB_NEWTON) {
			MPI_Barrier(mpi_comm);

			if(newtonParam.stabilizationMethodType == CONST_DIAG) { // Apply constant diagonal scaled by local volume
				double mySumVolume = 0.0;
				double myMinDiagValue = ABSURDLY_BIG_NUMBER, myMaxDiagValue = 0.0;
				for(int icv=0; icv<ncv; ++icv) {
					double diagValue = cv_volume[icv] / alpha_stab_newton;

					mySumVolume += cv_volume[icv];

					myMinDiagValue = min(myMinDiagValue, diagValue);
					myMaxDiagValue = max(myMaxDiagValue, diagValue);

					for(int ivar=0; ivar<m; ++ivar) {
						int row = icv*m + ivar;
						int diagIndex = jacMatrixSTL.get_diag_index(row);
						double oldValue = jacMatrixSTL.get_values(diagIndex);

						jacMatrixSTL.set_values(diagIndex, oldValue - diagValue); // Note that the Jacobian matrix is negative (semi-) definite
						++myNonDiagAdd_count; // Just for Statistics
					}
				}
				double sumVolume;
				MPI_Allreduce(&mySumVolume, &sumVolume, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
				double minDiagValue, maxDiagValue;
				MPI_Allreduce(&myMinDiagValue, &minDiagValue, 1, MPI_DOUBLE, MPI_MIN, mpi_comm);
				MPI_Allreduce(&myMaxDiagValue, &maxDiagValue, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);

				if(NcontrolEqns > 0 && mpi_rank == mpi_size-1) {
					double VolumeLambda = 1.0e-2 * AreaLambda * sumVolume;
					double diagValueLambda = VolumeLambda / alpha_stab_newton; // Use total volume for lambda
					for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
						lambda_tangent_inJacMat[iEqn] -= diagValueLambda;
						++myNonDiagAdd_count; // Just for Statistics
					}

					cout<<"           >> CONST_DIAG: DT = "<<alpha_stab_newton
					    <<" --> max diag (flow)="<<maxDiagValue<<", min diag (flow)="<<minDiagValue
					    <<"; diag for lambda="<<diagValueLambda<<endl;
				}

				MPI_Barrier(mpi_comm);
			} else if(newtonParam.stabilizationMethodType == CFL_BASED_DIAG) { // Apply diagonal scaled by CFL and local volume
				// Calculate dt for each CV
				double dt_min = calcDt(alpha_stab_newton);

				double my_dt_max = 0.0, my_vol_sum = 0.0;
				double myMinDiagValue = ABSURDLY_BIG_NUMBER, myMaxDiagValue = 0.0;
				for(int icv=0; icv<ncv; ++icv) {
					my_dt_max = max(my_dt_max, local_dt[icv]);
					my_vol_sum += cv_volume[icv];

					double diagValue = cv_volume[icv] / local_dt[icv];
					myMinDiagValue = min(myMinDiagValue, diagValue);
					myMaxDiagValue = max(myMaxDiagValue, diagValue);

					for(int ivar=0; ivar<m; ++ivar) {
						int row = icv*m + ivar;
						int diagIndex = jacMatrixSTL.get_diag_index(row);
						double oldValue = jacMatrixSTL.get_values(diagIndex);

						jacMatrixSTL.set_values(diagIndex, oldValue - diagValue); // Note that the Jacobian matrix is negative (semi-) definite
						++myNonDiagAdd_count; // Just for Statistics
					}
				}
				double dt_max, vol_sum;
				MPI_Allreduce(&my_dt_max,  &dt_max,  1, MPI_DOUBLE, MPI_MAX, mpi_comm);
				MPI_Allreduce(&my_vol_sum, &vol_sum, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
				double minDiagValue, maxDiagValue;
				MPI_Allreduce(&myMinDiagValue, &minDiagValue, 1, MPI_DOUBLE, MPI_MIN, mpi_comm);
				MPI_Allreduce(&myMaxDiagValue, &maxDiagValue, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);

				if(NcontrolEqns > 0 && mpi_rank == mpi_size-1) {
					double VolumeLambda = 1.0-2 * AreaLambda;  // V = 1 / |lambda_tangent|.
					double diagValueLambda = VolumeLambda / alpha_stab_newton;
					for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
						lambda_tangent_inJacMat[iEqn] -= diagValueLambda; // Note that the Jacobian matrix is negative (semi-) definite
						++myNonDiagAdd_count; // Just for Statistics
					}

					cout<<"           >> CFL_BASED_DIAG: CFL = "<<alpha_stab_newton
					    <<" --> max diag (flow)="<<maxDiagValue<<", min diag (flow)="<<minDiagValue
					    <<"; diag for lambda="<<diagValueLambda<<endl;
				}

				MPI_Barrier(mpi_comm);
			} else if(newtonParam.stabilizationMethodType == GENERAL_NEWTON_HU) {
				                                      // Apply a scaled diagonal matrix - scaled with rhs: See X.Wu, Applied Mathematics and Computation, vol.189, 2007.
				                                      // Basically finding x that minimize "0.5 * || diag(exp(w*(x-x*))) * f(x) ||_2^2"
				                                      //           or finding the root of "diag(exp(w*(x-x*))) * f(x) = 0"
				// Calculate the actual alpha
				double alpha = max( min(residNormTot, newtonParam.stabilizationAlpha), newtonParam.stabilizationAlphaEps);

				// Matrix
				for(int icv=0; icv<ncv; ++icv) {
					for(int ivar=0; ivar<m; ++ivar) {
						int row = icv*m + ivar;
						int diagIndex = jacMatrixSTL.get_diag_index(row);

						double oldValue = jacMatrixSTL.get_values(diagIndex);
						double addVal   = - max(MACHINE_EPS, alpha*fabs(rhs[row])) - newtonParam.stabilizationAlphaEps;
						jacMatrixSTL.set_values(diagIndex, oldValue + addVal); // Note that the Jacobian matrix is negative (semi-) definite

						++myNonDiagAdd_count; // Just for Statistics
					}
				}
			} else if(newtonParam.stabilizationMethodType == GENERAL_NEWTON_HUSEO) {
				                                      // Apply both scaled diagonal matrix and scaled rhs: See J.L.Hueso, J. Comp. and Applied Math., vol.224, 2009.
				                                      // More generalized version than X.Wu, 2007.
				                                      // Note: 1. We will use the same weight in the A-priori regularization for the exponent (m in the formula)
				                                      //       2. If this version is active, you must deactivate the A-priori regularization in the next part of Newton iteration.
				                                      // Basically finding x that minimize "0.5 * || diag(exp(w*(x-x*))) * f(x)^(1/m) ||_2^2"
				                                      //           or finding the root of "diag(exp(w*m*(x-x*))) * diag(m) * f(x) = 0"
				// Calculate the actual alpha
				double alpha = max( min(residNormTot, newtonParam.stabilizationAlpha), newtonParam.stabilizationAlphaEps);

				// The coefficients to calculate the exponent (m in the formulation) were already set-up at the beginning of this method (search for A-priori regularization)
				// If a region should not be modfied based on a priori knowledge (e.g. upstream of scramjet inlet), apply a weight function that is close to zero:
				//     Function: If x<=x_2, w = 10^(-digit) + (1.0-10^(-digit)) * exp(-beta*(x-x2)^2)
				//               Otherwise, w = 1.0
				//     i.e., The function is uni-directional Gaussian
				double minWeight = pow(10.0, -digitAprioriWeight); // If digitAprioriWeight==3, minWeight=0.001

				// Matrix
				for(int icv=0; icv<ncv; ++icv) {
					double xCoord = x_cv[icv][0];

					double exponent = 1.0;
//						if(xCoord <= x2AprioriWeight) {
//							exponent = minWeight + max( (1-minWeight)*exp(-betaAprioriWeight*pow(xCoord-x2AprioriWeight, 2.0)), MACHINE_EPS );
//							exponent = max(MACHINE_EPS, min(exponent, 1.0)); // Make sure that 0 < exponent <= 1.0
//						}

					for(int ivar=0; ivar<m; ++ivar) {
						int row = icv*m + ivar;
						int diagIndex = jacMatrixSTL.get_diag_index(row);

						double oldValue = jacMatrixSTL.get_values(diagIndex);
//							double addVal   = - max(MACHINE_EPS, alpha*exponent*fabs(rhs[row])) - stabilizationAlphaEps;
						double addVal   = -alpha*exponent*rhs[row] - newtonParam.stabilizationAlphaEps;
						jacMatrixSTL.set_values(diagIndex, oldValue + addVal); // Note that the Jacobian matrix is negative (semi-) definite

						++myNonDiagAdd_count; // Just for Statistics
					}
				}
			}

			MPI_Barrier(mpi_comm);
		}

		/******
		 ** solve the linear (NEWTON) system :
		 ** 1. Normal Newton:    q_{n+1} = q_{n} - J^{-1}*RHS
		 ** 2. Pseudo-arclength: Q_{n+1} = Q_{n} - L^{-1}*RESID
		 **                          Where L = [           J                  dF/dlambda
		 **                                      -weighted_q_tangent'  -weighted_lambda_tangent ]
		 **                                RESID = [ RHS(q, lambda)
		 **                                          Nres           ]
		 **                          c.f.) Nres = ds - [q_tangent', lambda_tangent']*W*[dq; dlambda],  where W is the (diagonal) weight matrix
		 ******/
		// solve the sub-system, PHI = L^(-1)*RHS
		int nIter;
		double absResid;
		vector<pair<int, double> > kspMonitorHistory; 	kspMonitorHistory.clear();

		double wtime0, wtimeLS;

// TO DO: The following part has a significant error!
//        The PETSc package is frozen while solving the linear system using KSP --
//        For details, see the solveGMRES_solver() method of the PetscSolver2 class in PetscSolver2.h
		if(howToCalcJac == ROW_1D) {
			if(firstCall && iterNewton==startingNewtonIter+1 && debugLevel>0 && mpi_rank==0)
				wtime0 = MPI_Wtime();

			// Note: If you want to print out the matrix stored in PETSc on the screen, set the ShowPetscMatrixMatlab variable as true
			//       ShowPetscMatrixMatlab is defined as a static variable in the PetscSolver2.h file
			//       e.g. if(iterNewton==1 && NcontrolEqns>0) ShowPetscMatrixMatlab = true;
			bool useOldSolnAsInitGuess = false; // If true, use the previous solution vector as the initial guess.
			                                    // Otherwise, apply PC(pre-conditioner) to get the inital guess (the Knoll trick) -- For Bifurcation, this works better.

			solveLinSysNSCoupled2<MatComprsedSTL>(phi, jacMatrixSTL, rhs, nIter, absResid, kspMonitorHistory,
					useOldSolnAsInitGuess, newtonParam.zeroAbsLS, newtonParam.zeroRelLS, newtonParam.maxIterLS, nScal, monitorConvergInterval,
					NcontrolEqns, ncv_gg, q_tangent_inJacMat, lambda_tangent_inJacMat, step, iterNewton);
			// Note that you should pass nScal to this function instead of m(=5+nScal)

			if(firstCall && iterNewton==startingNewtonIter+1 && debugLevel>0 && mpi_rank==0)
				wtimeLS = MPI_Wtime();
			if(firstCall && iterNewton==startingNewtonIter+1 && debugLevel>0 && mpi_rank==0)
				cout<<"   > "<<classID<<"::"<<funcID<<"(): Runtime for the LS solver   [sec] = "<<wtimeLS - wtime0<<endl;
		} else if(howToCalcJac == ORDINARY_2D) {
			// Note: If you want to print out the matrix stored in PETSc on the screen, set the ShowPetscMatrixMatlab variable as true
			//       ShowPetscMatrixMatlab is defined as a static variable in the PetscSolver2.h file
			//       e.g. if(iterNewton==1 && NcontrolEqns>0) ShowPetscMatrixMatlab = true;
			bool useOldSolnAsInitGuess = false; // If true, use the previous solution vector as the initial guess.
			                                    // Otherwise, apply PC(pre-conditioner) to get the inital guess (the Knoll trick) -- For Bifurcation, this works better.
			int saveConvergInterval = int(max(monitorConvergInterval/10, 1.01));

			solveLinSysNSCoupled2<MatComprsed>(phi, jacMatrix, rhs, nIter, absResid, kspMonitorHistory,
					useOldSolnAsInitGuess, newtonParam.zeroAbsLS, newtonParam.zeroRelLS, newtonParam.maxIterLS, nScal, monitorConvergInterval,
					NcontrolEqns, ncv_gg, q_tangent_inJacMat, lambda_tangent_inJacMat, step, iterNewton);
			// Note that you should pass nScal to this function instead of m(=5+nScal)
		} else {
			throw(PSALC_ERROR_CODE);
		}

		if(nIter == newtonParam.maxIterLS  &&  absResid > newtonParam.zeroAbsLS) {
			if(debugLevel>0 || (newtonParam.incIterLS==0 && newtonParam.incIterLS==1.0)) { // TO DO
				// Show the warning message on the screen
				double trueAbsResid = petscSolver2->calcTrueResidualNorm(); // Since "absResid" contains the last approximate & left-preconditioned residual norm,
				                                                            // we want to calculate the true residual here.

				if(mpi_rank==0)
					cout<<"           >> WARNING in "<<classID<<"::"<<funcID<<"(): Linear solver reaches MAX ITER (iterNewton="<<iterNewton<<"): nIter="<<nIter<<", absResid="<<trueAbsResid<<endl; // absResid<<endl;

				// Write the residual history on a file
				if(mpi_rank==0) {
					if(!kspMonitorHistory.empty()) {
						FILE *fp = fopen(LS_CONVERGENCE_FILENAME, "a");

						size_t minIter  = kspMonitorHistory.size()-1;
						double minResid = ABSURDLY_BIG_NUMBER;
						for(size_t i=0; i<kspMonitorHistory.size(); ++i) {
							// Save the convergence history on a file
							fprintf(fp, "%7d, %6d, %10.4e\n", iterNewton, kspMonitorHistory[i].first, kspMonitorHistory[i].second);

							// Find the minimum
							if(kspMonitorHistory[i].second < minResid) {
								minIter  = kspMonitorHistory[i].first;
								minResid = kspMonitorHistory[i].second;
							}
						}
						fclose(fp);

						cout<<"           >> KSP linear solver residual history: "<<endl;
						for(size_t i=0; i<kspMonitorHistory.size(); ++i) {
							if(kspMonitorHistory[i].first % (monitorConvergInterval*10) == 0) {
								if(kspMonitorHistory[i].first == minIter)
									printf("             %5d\t %11.4e (MINIMUM) \n", kspMonitorHistory[i].first, kspMonitorHistory[i].second);
								else
									printf("             %5d\t %11.4e\n", kspMonitorHistory[i].first, kspMonitorHistory[i].second);
							} else {
								if(kspMonitorHistory[i].first == minIter) // Always show the minimul residual even if kspMonitorHistory[i].first % monitorConvergInterval != 0.
									printf("               %5d\t %11.4e (MINIMUM) \n", kspMonitorHistory[i].first, kspMonitorHistory[i].second);
							}
						}
					}
				}
			}
		}

        // (Re)Initialize relaxation
		relaxation = getDoubleParam("UNDER_RELAXATION", "1.0");
		MPI_Barrier(mpi_comm);

		/******
		 ** Trust-region regularization
		 ******/
        double mySumPhiSq= 0.0, sumPhiSq;

        double *weightVec = new double [5+nScal];
		for(int icv=0; icv<ncv; ++icv) {
			getFlowWeightsForInnerProduct(weightVec, 1.0, icv, 5+nScal);  // 1.0 --> without volume scaling
			int startIndex = icv*(5+nScal);
			for(int i=0; i<5+nScal; ++i)
				mySumPhiSq += pow(weightVec[i] * phi[startIndex+i], 2.0);
		}
		delete [] weightVec;
		if(mpi_rank==mpi_size-1 && NcontrolEqns>0)
			for(int i=0; i<NcontrolEqns; ++i)
				mySumPhiSq += weightLambda * pow(phi[ncv*(5+nScal)+i], 2.0);

		MPI_Allreduce(&mySumPhiSq, &sumPhiSq, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

		double deltaQnorm_beforeRelax = sqrt(sumPhiSq);
        deltaQnorm = relaxation * deltaQnorm_beforeRelax;

        if(deltaQnorm >= newtonParam.trustRegionSize) {
        	if(mpi_rank==0)
        		cout<<"WARNING in "<<classID<<"::"<<funcID<<"(): After solving the linear system,"<<endl
        		    <<"                                    deltaQ="<<deltaQnorm<<" is greater than trustRegionSize="<<newtonParam.trustRegionSize<<endl
        		    <<"                                    Reduce relaxation from "<<relaxation;
        	relaxation *= newtonParam.trustRegionSize / deltaQnorm;
        	if(mpi_rank==0)
        		cout<<" to "<<relaxation<<endl;

            updateFlow_primVars(q, phi, -relaxation, nScal); // note: You should not update q here! (because of backtracking)
                                                             //       Only rho, rhou, rhoE, and scalars should be updated here
            if(NcontrolEqns > 0)
				updateLambda(q, phi, -relaxation, NcontrolEqns, nScal);
        }
        MPI_Barrier(mpi_comm);

		/******
		 ** Check if negative rho or press occurs: If occurs, update both the flow field and the "relaxation"
		 ******/
		int kine_index = getScalarTransportIndex("kine");

        int maxIterFAcheck = 5;
        int negativeValCount_CV, negativeValCount_FA;
        relaxation = backtrackForNegativeVals(negativeValCount_CV, negativeValCount_FA, newtonParam.clipParameter, newtonParam.safeParameter,
        		relaxation, kine_index, nScal, q, phi, maxIterFAcheck);

        int negativeValCount = negativeValCount_CV + negativeValCount_FA;

        if(negativeValCount > 0) {
            if(mpi_rank==0) {
                printf("           >> WARNING! %s::%s(): Checking Negative rho, press, or kine - Negative occurs at total %d CVs / %d FAs\n", classID.c_str(), funcID.c_str(), negativeValCount_CV, negativeValCount_FA);
                printf("                                                                                        Reduce RELAXATION to %.3e\n", relaxation);
            }

            updateFlow_primVars(q, phi, -relaxation, nScal); // note: You should not update q here! (because of backtracking)
                                                             //       Only rho, rhou, rhoE, and scalars should be updated here

            if(NcontrolEqns > 0)
				updateLambda(q, phi, -relaxation, NcontrolEqns, nScal);
        }
        MPI_Barrier(mpi_comm);

//        /******
//         ** Check if Delta(Lambda) is too large
//         ******/
//        double relaxBeforeDlambda = relaxation;
//
//        // If delta_lambda is not that big, then the backtracking algorithm can take care of it.
//        // However, if delta_lambda is extremely large, we need to do something here.
//		if(NcontrolEqns>0) {
//			if(mpi_rank==mpi_size-1) {
//                for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
//                    double delta_lambda      = relaxBeforeDlambda * phi[ncv*(5+nScal)+iEqn];
//                    double delta_lambda_prev = lambda1[iEqn] - lambda0[iEqn];
//                    if(fabs(delta_lambda) > 7.77*fabs(delta_lambda_prev)) {
//                        cout<<"WARNING in "<<classID<<"::"<<funcID<<"(): delta_lambda["<<iEqn<<"]="<< delta_lambda
//                        	<<" is greater than 777% of lambda["<<iEqn<<"]-lambda1["<<iEqn<<"]="<<delta_lambda_prev<<")"<<endl;
//
//                        if(fabs(delta_lambda_prev)>1.0e-10) {
//                        	cout<<"                                    reduce relaxation from "<<relaxation;
//                            relaxation = min(max(RELAXATION_EPS, 7.77/fabs(delta_lambda/delta_lambda_prev)*relaxation), relaxBeforeDlambda);
//                            cout<<" to "<<relaxation<<endl;
//                        }
//                    }
//                }
//			}
//		}
//        MPI_Bcast(&relaxation, 1, MPI_DOUBLE, mpi_size-1, mpi_comm);

		// update the flow field: u_{n+1} = q_{n} - phi
		updateFlow_primVars(q, phi, -relaxation, nScal); // note: You should not update q here! (because of backtracking)
		                                                 //       Only rho, rhou, rhoE, and scalars should be updated here

        if(NcontrolEqns>0)
			updateLambda(q, phi, -relaxation, NcontrolEqns, nScal);
        MPI_Barrier(mpi_comm);

		/******
		 ** Calculate RHS and residuals using JOE: This is to launch backtracking and other stabilization algorithm
		 ** Note: It looks a redundant step, but remember the fact that the RHS calculation by JOE is very fast
		 ******/

		bool useBarrier = true;
		countNegative = calcRhsWithBarrier(rhs, useBarrier);

		calcResidualsFrom1Drhs(residNormVecFlow, rhs, whichNorm);
		residNormTotFlowOnly = calcSumResidual(residNormVecFlow, whichNorm);

#ifdef USE_TOT_NORM_WITH_LAMBDA
		if(NcontrolEqns>0) {
			// Note: 1. Even though this part of the code is designed to handle multiple parameters (or lambdas), the actual formulation here is only for one parameters.
			//       2. You don't have to pass weighted_q_tangent instead of q_tangent since it will be weighted in getFlowWeightsForInnerProduct() called by vecDotVecWithWeight().
//			calcNres(Nres, q, rhs, m, NcontrolEqns, q1, q_tangent, lambda_tangent, lambda1, weightLambda, arcLength);
			calcNresForJOE(Nres, m, NcontrolEqns, q1, q_tangent, lambda_tangent, lambda, lambda1, weightLambda, arcLength);
											// Use calcNresForJOE() instead of calcNres() since q has not been updated yet!

			for(int iParam=0; iParam<NcontrolEqns; ++iParam) {
				if(isNaN(Nres[iParam])) {
					if(mpi_rank==0)
						cerr<<"ERROR in "<<classID<<"::"<<funcID<<"(): iterNewton=="<<iterNewton<<", The residual of the tangential condition becomes NaN for iParam="<<iParam<<endl;
					throw(PSALC_ERROR_CODE);
				}
			}

			MPI_Barrier(mpi_comm);
		}

		if(NcontrolEqns>0) {
			residNormTot = updateTotResidualWithNres(residNormTotFlowOnly, Nres, NcontrolEqns, whichNorm);

			if(debugLevel>0 && mpi_rank==0)
				printf("           >> Residual: flow RHS = %.5e, tangential RHS[0] = %.5e --> Total residual = %.5e\n", residNormTotFlowOnly, fabs(Nres[0]), residNormTot);
		} else {
			residNormTot = residNormTotFlowOnly;

			if(debugLevel>0 && mpi_rank==0)
				printf("           >> Residual of the flow RHS = %.5e\n", residNormTot);
		}
#else
		residNormTot = residNormTotFlowOnly;

		if(debugLevel>0 && mpi_rank==0)
			printf("           >> Residual of the flow field RHS = %.5e \n", residNormTot);
#endif

		/******
		 ** Backtracking (a.k.a. "line search")
		 ** 1. The backtracking will be triggered if there is negative rho/press or residuals increase.
		 ** 2. First, calculate the relaxation size by calling backtrackWithJOE_calcRelaxAndRHS() - The primary variables (rho, rhou, rhoE) will be also updated, but Q won't
		 ** 3. Second, re-calculate RHS by using ADOLC (note that in the backtracking, not be used ADOLC but normal JOE)
		 ******/

		double relaxBeforeBacktrack = relaxation;
		double residNormTotFresh = residNormTot; // This is just to provide information to the user: check the residual before all the treatments (e.g. backtracking) in the following

		bool btNotConverged   = false;
		bool forgiveBacktrack = false;
		if(newtonParam.backtrackMaxIter > 0) {
			bool needBacktrackingDueToResid = checkBacktrackFromResidInc(residNormTotOld, residNormTot, min(relaxation, newtonParam.backtrackRelax_UpperBound), newtonParam.backtrackBarrierCoeff);

			if(countNegative>0 || needBacktrackingDueToResid) { // The primary reason to call the backtracking algorithm is negative scalars (e.g. rho, press, etc.)
				++succeessiveBacktrack;

				// Check if backtracking is required: If the conditions are satisfied, backtracking will be omitted
				if(newtonParam.skipBT_freq>0)
					forgiveBacktrack = ( (newtonParam.skipBT_firstITer && iterNewton==startingNewtonIter+1) || succeessiveBacktrack%newtonParam.skipBT_freq==0 );
				else
					forgiveBacktrack = ( newtonParam.skipBT_firstITer && iterNewton==startingNewtonIter+1 );

				// Check if the user want to skip backtracking for first few times
				if(forgiveBacktrack && countNegative==0) { // Note: If countNegative>0, backtracking will be active anyway
					relaxation = 0.2 * ( newtonParam.backtrackRelax_LowerBound + min(relaxation, newtonParam.backtrackRelax_UpperBound) );

					updateFlow_primVars(q, phi, -relaxation, nScal);
					if(NcontrolEqns > 0)
						updateLambda(q, phi, -relaxation, NcontrolEqns, nScal);

					if(debugLevel > 0 && mpi_rank == 0)
						cout<<"                  Forgive BACKTRACKING to escape from the current point: Instead, set relaxation="<<relaxation<<" (lambda="<<lambda[0]<<")"<<endl;
				} else {
					// Show the reason to use the backtracking algorithm on the screen if the debug level is high
					if(debugLevel > 0) {
						if(mpi_rank==0) {
							cout<<"                  Calling BACKTRACKING ";
							if(debugLevel > 1) {
								cout<<"due to";
								if(countNegative > 0)
									printf(" negative vars=%d", countNegative);
								if(needBacktrackingDueToResid)
									printf(" residual=%.8e", residNormTot);
							}
							cout<<endl;
						}
					}
					MPI_Barrier(mpi_comm);

					// Launch the backtracking algorithm -------
					// First, calculate the relaxation size by calling backtrackWithJOE_calcRelaxAndRHS()
					relaxation = backtrackWithJOE_calcRelaxAndRHS(rhs, q, phi, btNotConverged,
							newtonParam.backtrackMaxIter, relaxation, newtonParam.backtrackRelax_LowerBound, newtonParam.backtrackRelax_UpperBound,
							residNormVecFlowOld, residNormVecFlow, whichNorm, newtonParam.backtrackBarrierCoeff,
							NcontrolEqns, q_tangent, lambda_tangent, weightLambda, q1, Nres, lambda1, arcLength);
							// Since the reduction algorithm based on delta_lambda is heuristics, use "relaxBeforeDlambda" instead of "relaxation"

					calcResidualsFrom1Drhs(residNormVecFlow, rhs, whichNorm);
					residNormTotFlowOnly = calcSumResidual(residNormVecFlow, whichNorm);

					// Second, update the primary variables and lambda
					updateFlow_primVars(q, phi, -relaxation, nScal);
					if(NcontrolEqns>0) {
#ifdef UPDATE_LAMBDA_IN_BACKTRACK
						updateLambda(q, phi, -relaxation, NcontrolEqns, nScal);
#endif
					}

#ifdef USE_TOT_NORM_WITH_LAMBDA
					if(NcontrolEqns>0) {
						// Finally, update the total residual using the tangential condition
						calcNresForJOE(Nres, m, NcontrolEqns, q1, q_tangent, lambda_tangent, lambda, lambda1, weightLambda, arcLength);
								// Use calcNresForJOE() instead of calcNres() since q has not been updated yet!
						MPI_Barrier(mpi_comm);
						residNormTot = updateTotResidualWithNres(residNormTotFlowOnly, Nres, NcontrolEqns, whichNorm);
					}
#endif
				}
			} else {
				succeessiveBacktrack = 0;
			}
		}

		/******
		 ** Update the q vector based on relaxation
		 **     Note that the primary variables (rho, rhou, rhoE) were already updated. q is being updated now because of the backtracking
		 ******/
		if(mpi_rank == mpi_size-1) {
#ifdef UPDATE_LAMBDA_IN_BACKTRACK
			updateFlow_1Dvec(q, phi, -relaxation, ncv*m+NcontrolEqns);
#else
			updateFlow_1Dvec(q, phi, -relaxation, ncv*m);
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
				q[ncv*m+iEqn] = lambda[iEqn];
#endif
		} else
			updateFlow_1Dvec(q, phi, -relaxation, ncv*m);

		/******
		 ** Give the user a chance to apply an ad-hoc modification to the solution
		 ******/
		bool solnModified = temporalHookNewton(q, lambda, phi, relaxation);
		if(solnModified) {
			updateCvDataG1G2(rho, REPLACE_DATA);
			updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
			updateCvDataG1G2(rhoE, REPLACE_DATA);

			for (int iScal = 0; iScal < nScal; iScal++)
				updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);

			calcStateVariables();
		}

		/******
		 ** Check if a modified Newton's method should be used in the next iteration and clear the Jacobian matrix if required
		 ** 1. If next iter % freq == 0                        : Force BASIC
		 ** 2. If residual > modifNewtonRelResid*initResidNorm : Force BASIC
		 ******/
		if (modifiedNewtonMethod == BASIC) {
			useExactJac = true;
		} else {
			if ((iterNewton+1)%modifNewtonFreq == 0) {
				useExactJac = true;
				petscSolver2->setPCreuse(false);
				jacMatrix.clear();
			} else if (residNormTot < modifNewtonRelResid*initResidNorm) {
				useExactJac = true;
				petscSolver2->setPCreuse(false);
				jacMatrix.clear();
			} else {
				useExactJac = false;
				petscSolver2->setPCreuse(true);
			}
		}

		/******
		 ** Calculate RHS and Jacobian for the next iteration
		 ******/
		if(useExactJac) { // Calculate rhs & Jacobian by using IKE
			try {
				if(howToCalcJac == ROW_1D) {
					jacMatrixSTL.clear();
					countNegative = calcJacobian1DAD(jacMatrixSTL, rhs, debugLevel, NcontrolEqns);
					myNnzJac = jacMatrixSTL.get_nnz();
				} else if(howToCalcJac == ORDINARY_2D) {
					jacMatrix.clear();
					countNegative = calcJacobianAD(jacMatrix, rhs, debugLevel, NcontrolEqns);
					myNnzJac = jacMatrix.get_nnz();
				} else { // ERROR
					freeMemGetSteadySolnByNewton(phi, residNormVecFlow, residNormVecFlowOld, Nres, jacMatrix, jacMatrixSTL);
					throw(PSALC_ERROR_CODE);
				}
				MPI_Allreduce(&myNnzJac, &nnzJac, 1, MPI_INT, MPI_SUM, mpi_comm);
			}
			catch(int e) { // Catch and re-throw
				freeMemGetSteadySolnByNewton(phi, residNormVecFlow, residNormVecFlowOld, Nres, jacMatrix, jacMatrixSTL);
				throw(e);
			}
			catch(...) { // Catch and re-throw
				freeMemGetSteadySolnByNewton(phi, residNormVecFlow, residNormVecFlowOld, Nres, jacMatrix, jacMatrixSTL);
				throw;
			}

			if(countNegative>0) {
				freeMemGetSteadySolnByNewton(phi, residNormVecFlow, residNormVecFlowOld, Nres, jacMatrix, jacMatrixSTL);
				if(mpi_rank==0) cout<<"ERROR! "<<classID<<"::"<<funcID<<"(): After negative vars = "<<countNegative<<endl;
				throw(PSALC_ERROR_CODE);
			}
		} else { // Calculate rhs, residNorm, and residNormTotFlowOnly by using normal JOE
			bool useBarrier = true;
			int countNegative = calcRhsWithBarrier(rhs, useBarrier);
			calcResidualsFrom1Drhs(residNormVecFlow, rhs, whichNorm);
			residNormTotFlowOnly = calcSumResidual(residNormVecFlow, whichNorm);
		}

		calcResidualsFrom1Drhs(residNormVecFlow, rhs, whichNorm);

		/******
		 ** fill up the last element of rhs (Nres == q_tangent*(q_guess-q1) + weightLambda*lambda_tangent*(lambda_guess-lambda1) - ds) for mpi_rank == mpi_size-1
		 ** Note: Even though this part of the code is designed to handle multiple parameters (or lambdas), the actual formulation here is implemented for one parameter now.
		 ******/
//		if(NcontrolEqns>0) {
//			calcNres(Nres, q, rhs, m, NcontrolEqns, q1, q_tangent, lambda_tangent, lambda1, weightLambda, arcLength);
//
//			for(int iParam=0; iParam<NcontrolEqns; ++iParam)
//				if(fabs(Nres[iParam]/arcLength) > 1.0e-10 * sqrt(double(NcontrolEqns*cvora[mpi_size]))) {
//					if(mpi_rank==0) {
//						double arcLengthNew = arcLength - Nres[iParam];
//						double lambdaPart = weightLambda * lambda_tangent[iParam] * (lambda[iParam] - lambda1[iParam]);
//						cout<<"WARNING in "<<classID<<"::"<<funcID<<"(): The tangential condition is NOT satisfied: Nres["<<iParam<<"] = "<<Nres[iParam]<<" (target arclength="<<arcLength<<")"<<endl
//							<<"                                    (lambdaPart + flowPart = "<<lambdaPart<<" + "<<arcLengthNew-lambdaPart<<" = "<<arcLengthNew<<")"<<endl;
//					}
//				}
//
//			if(mpi_rank == mpi_size-1) {
//				for(int iParam=0; iParam<NcontrolEqns; ++iParam)
//					rhs[ncv*m+iParam] = Nres[iParam];
//			}
//		}

		if(NcontrolEqns>0) {
			residNormTot = updateTotResidualWithNres(residNormTotFlowOnly, Nres, NcontrolEqns, whichNorm);

			if(debugLevel>0 && mpi_rank==0)
				printf("           >> Residual: flow RHS = %.5e, tangential RHS[0] = %.5e --> Total residual = %.5e\n", residNormTotFlowOnly, fabs(Nres[0]), residNormTot);
		} else {
			residNormTot = residNormTotFlowOnly;

			if(debugLevel>0 && mpi_rank==0)
				printf("           >> Residual of the flow RHS = %.5e\n", residNormTot);
		}

		/******
		 ** Update residNormVecFlowOld and residNormTotOld:
		 ** Note: "residNormTotOld" is used only to compare residuals from RHS calculation without the tangential condition,
		 **       thus it should be updated before the tangential condition.
		 ******/
		if(residNormTot < residNormTotOld) {
			for(int i=0; i<5+nScal; ++i)
				residNormVecFlowOld[i] = residNormVecFlow[i];
			residNormTotOld = residNormTot;
		}

		/******
		 ** Calculate the 2-norm of the final delta_Q
		 ******/
//		deltaQnorm = relaxation*sqrt(sumPhiSq);
		deltaQnorm = relaxation*deltaQnorm_beforeRelax;

		/******
		 ** If debug level is high, then show the current status
		 ******/
		if(debugLevel > 0 && mpi_rank == 0) {
			if(iterNewton>1 && (iterNewton-1)%10==0) {
				printf("        ITER         RESIDUAL     DELTA_Q       NNZ_JAC  LS_ITER       LS_RESID  NEWTON(JAC/RESID)");
				if(NcontrolEqns > 0) {
					printf("          LAMBDA\n");
				} else {
					printf("\n");
				}
			}
			printf("       %5d  %15.8e  %10.3e  %12d  %7d  %13.6e  ", iterNewton, residNormTot, deltaQnorm, nnzJac, nIter, absResid);
			if(useExactJacTempo) cout<<"(Exact/";
			else                 cout<<"("<<modifiedNewtonMethod<<"/";
			if(useExactJac) cout<<"Exact)";
			else            cout<<modifiedNewtonMethod<<")";
			if(NcontrolEqns > 0) {
				printf("     %15.8e\n",lambda[0]);
			} else {
				printf("\n");
			}

			if(debugLevel > 0) {
				printf("           RESID: ");
				for(int i=0; i<5+nScal; ++i)
					printf("%10.3e", residNormVecFlow[i]);
				for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
					printf(" | %10.3e", Nres[iEqn]);
				printf("\n");
			}

			printf("\n");
		}

		/******
		 ** Print out the residual history on a file
		 ******/
		if(mpi_rank==0) {
			FILE *fp = fopen(NEWTON_STATUS_FILENAME, "a");

			double nresNorm = 0.0;
			for(int iParam=0; iParam<NcontrolEqns; ++iParam)
				nresNorm += fabs(Nres[iParam]);

			if(NcontrolEqns > 0)
				fprintf(fp, "%4d, %5d, %12d, %13.6e, %11.4e, %13.6e, %9.2e", step, iterNewton, nnzJac, residNormTotFresh, relaxation, residNormTot, nresNorm);
			else
				fprintf(fp, "%4d, %5d, %12d, %13.6e, %11.4e, %13.6e,          ", step, iterNewton, nnzJac, residNormTotFresh, relaxation, residNormTot);

			if(NcontrolEqns>0)
				fprintf(fp, ", %13.6e\n", lambda[0]);
			else
				fprintf(fp, "\n");
			fclose(fp);
		}

		/******
		 ** update currIter & useExactJacTempo
		 ******/
		++iterNewton;

		useExactJacTempo = useExactJac; // Currently, useExactJacTempo is not used at all

		MPI_Barrier(mpi_comm);

		/******
		 ** Check if the iterations are done
		 ******/
		if(residNormTot <= tolNewton) {
			++countConvergedSteps;

			if(moreNTstepsBelowConv_moreSteps <= 0)
				done = true;
			else {
				if(countConvergedSteps > moreNTstepsBelowConv_moreSteps)
					done = true;
				else {
					if(deltaQnorm <= moreNTstepsBelowConv_deltaQ)
						done = true;
				}
			}

			if(residNormTot <= moreNTstepsBelowConv_resid)
				done=true;

			if(iterNewton > startingNewtonIter+maxIterNewton)
				done = true;

			if(mpi_rank==0 && countConvergedSteps==1) {
				if(!done)
					cout<<"Newton-solver converged but keep running the simulation since MORE_STEPS_BELOW_CONV="<<moreNTstepsBelowConv_moreSteps<<", DELTA_Q="<<moreNTstepsBelowConv_deltaQ<<endl;
				else
					cout<<"Even though MORE_STEPS_BELOW_CONV="<<moreNTstepsBelowConv_moreSteps<<", do not launch more Newton steps since other conditions have met"<<endl;
			}
		}

		if(newtonParam.stabilizationMethodType != NO_STAB_NEWTON) { // Apply diagonal scaled by CFL and local volume
			// Check ramping
			if (!checkParam("RAMP_STAB_NEWTON_ALPHA")) {
				ParamMap::add("RAMP_STAB_NEWTON_ALPHA AFTER_ITER=10  INTERVAL_ITER=10  FACTOR_INC=1.0  TARGET_ALPHA=1.0");    // add default: no increase of CFL number!
				if (mpi_rank == 0)
					cout << "WARNING: added \"RAMP_STAB_NEWTON_ALPHA AFTER_ITER=10  INTERVAL_ITER=10  FACTOR_INC=1.0  TARGET_ALPHA=1.0\" to parameter map" << endl;
			}

			int startIncAlpha    = getParam("RAMP_STAB_NEWTON_ALPHA")->getInt("AFTER_ITER");
			int intervalIncAlpha = getParam("RAMP_STAB_NEWTON_ALPHA")->getInt("INTERVAL_ITER");
			double incAlpha    = getParam("RAMP_STAB_NEWTON_ALPHA")->getDouble("FACTOR_INC");
			double targetAlpha = getParam("RAMP_STAB_NEWTON_ALPHA")->getDouble("TARGET_ALPHA");

			if ((iterNewton >= startIncAlpha) && (iterNewton%intervalIncAlpha == 0)) {
				if(incAlpha>1.0 && alpha_stab_newton<targetAlpha)
					if(relaxation > 1.0/incAlpha)
						alpha_stab_newton *= incAlpha;
				if(incAlpha<1.0 && alpha_stab_newton>targetAlpha)
					if(relaxation > incAlpha)
						alpha_stab_newton *= incAlpha;
			}
		}
////IKJ
//writeData(step, iterNewton);
	}

	if(debugLevel>0 && mpi_rank == 0)
		printf("Newton solver converged after %5dth outer iter: residual = %12.5e \n", iterNewton-1, residNormTot);

	// Write the Jacobian matrix on a file for post-processing(eigen-decomposition)
	if(writeJacOnFile) {
		if(newtonParam.stabilizationMethodType != NO_STAB_NEWTON) { // i.e., If the Jacobian matrix was modified in order not to have a singular matrix
			if(newtonParam.stabilizationMethodType == CONST_DIAG) {
				for(int icv=0; icv<ncv; ++icv) {
					double diagValue = cv_volume[icv] * newtonParam.stabilizationAlpha;

					for(int ivar=0; ivar<m; ++ivar) {
						int row = icv*m + ivar;
						int diagIndex = jacMatrixSTL.get_diag_index(row);
						double oldValue = jacMatrixSTL.get_values(diagIndex);

						jacMatrixSTL.set_values(diagIndex, oldValue + diagValue);
					}
				}
				// Note: You don't need to worry about the tangential equation because it is not a part of jacMatrixSTL
			} else if(newtonParam.stabilizationMethodType == CFL_BASED_DIAG) { // Apply diagonal scaled by CFL and local volume
				double dt_min = calcDt(newtonParam.stabilizationAlpha);  // Take CFL = stabilizationAlpha

				for(int icv=0; icv<ncv; ++icv) {
					double diagValue = cv_volume[icv] / local_dt[icv];

					for(int ivar=0; ivar<m; ++ivar) {
						int row = icv*m + ivar;
						int diagIndex = jacMatrixSTL.get_diag_index(row);
						double oldValue = jacMatrixSTL.get_values(diagIndex);

						jacMatrixSTL.set_values(diagIndex, oldValue + diagValue); // Note that the Jacobian matrix is negative (semi-) definite
					}
				}
				// Note: You don't need to worry about the tangential equation because it is not a part of jacMatrixSTL
			} else {
				if(howToCalcJac == ROW_1D) {
					jacMatrixSTL.clear();
					countNegative = calcJacobian1DAD(jacMatrixSTL, rhs, debugLevel, NcontrolEqns);
				} else if(howToCalcJac == ORDINARY_2D) {
					jacMatrix.clear();
					countNegative = calcJacobianAD(jacMatrix, rhs, debugLevel, NcontrolEqns);
				}
			}
		}

		if(!useExactJacTempo) { // i.e., If an approximated Jacobian is used
			jacMatrix.clear();

			if(howToCalcJac == ROW_1D) {
				jacMatrixSTL.clear();
				countNegative = calcJacobian1DAD(jacMatrixSTL, rhs, debugLevel, NcontrolEqns);
			} else if(howToCalcJac == ORDINARY_2D) {
				jacMatrix.clear();
				countNegative = calcJacobianAD(jacMatrix, rhs, debugLevel, NcontrolEqns);
			}
		}

		stringstream ss;
		ss<<"jacMatrix"<<step<<".bin";
		string filename = ss.str();

		if(jacMatrixSTL.empty()) {
			if(jacMatrix.empty()) { // This is an error!! Both jacMatrixSTL and jacMatrix are empty!
				if(mpi_rank==0) {
					cout<<"ERROR! "<<classID<<"::"<<funcID<<"(): Both jacMatrixSTL and jacMatrix are empty!"<<endl;
					cout<<"                                                 Cannot write the Jacobian matrix on a file"<<endl;
				}
				throw(PSALC_ERROR_CODE);
			}
			writeMatrixBinaryParallel<MatComprsed>(filename, jacMatrix, nScal, cvora, nbocv_v_global);
		} else {
			writeMatrixBinaryParallel<MatComprsedSTL>(filename, jacMatrixSTL, nScal, cvora, nbocv_v_global);
		}
	}

	// Write a summary of the Bifurcation code on a file
	if(mpi_rank==0) {
		FILE *fp = fopen(BIFUR_SUMMARY_FILENAME, "a");
		if(lambda1 != NULL)
			fprintf(fp, "%4d, %11.4e, %11.4e, %11.4e, %8d, %11.4e, %.15f\n", step, arcLength, initResidNorm, residNormTot, iterNewton, lambda[0]-lambda1[0], lambda[0]);
		else
			fprintf(fp, "%4d,            , %11.4e, %11.4e, %8d,            , %.15f\n", step, initResidNorm, residNormTot, iterNewton, lambda[0]);
		fclose(fp);
	}

	/* free memory */
	freeMemGetSteadySolnByNewton(phi, residNormVecFlow, residNormVecFlowOld, Nres, jacMatrix, jacMatrixSTL);

	if(q_tangent_inJacMat != NULL)       freeMem2D(q_tangent_inJacMat, 0, NcontrolEqns-1, 0, (ncv*m)-1);
	if(lambda_tangent_inJacMat != NULL)  delete [] lambda_tangent_inJacMat;

	firstCall = false;

	MPI_Barrier(mpi_comm);
}

/*
 * Method: writeEigenPairsParallel
 * -------------------------------
 * Write eigen-pairs for an van Gogh analysis
 * Format:
 *   Header:
 *       step (int), NcontrolEqns (int), lambda (double[NcontrolEqns]), nVars (int), nev (int), eigenvalues (double[nev*2])
 *       mpi_size (int), cvora (int[mpi_size+1]), xcvMinArray(double[mpi_size*3]), xcvMaxArray(double[mpi_size*3])
 *   Body: For each mpi_rank,
 *       x_cv
 *       eigen-vectors   (variable names = directEvecsReal, directEvecsImag)
 *       adjoint vectors (variable names = adjointEvecsReal, adjointEvecsImag)
 *   Foot:
 *       EOF_ERROR_CHECK_CODE
 */
void RembrandtWithModels::writeEigenPairsParallel(const string &filename, const int step, const int nev) {
	char *filenameArray = new char[filename.size() + 1]; // Since MPI_File_open() doesn't accept "const char *", you cannot use std::string.c_str()
	std::copy(filename.begin(), filename.end(), filenameArray);
	filenameArray[filename.size()] = '\0';

	writeEigenPairsParallel(filenameArray, step, nev);

	delete [] filenameArray;
}

void RembrandtWithModels::writeEigenPairsParallel(char filename[], const int step, const int nev) {
	// Number of variables
	int nScal = scalarTranspEqVector.size();
	int nVars = 5 + nScal; // number of variables (rho, rhou, rhov, rhow, rhoE, and scalars)

	// xcvMin and xcvMax
	double *xcvMinArray = new double [mpi_size*3];
	double *xcvMaxArray = new double [mpi_size*3];
	MPI_Gather(xcvMin, 3, MPI_DOUBLE, xcvMinArray, 3, MPI_DOUBLE, 0, mpi_comm);
	MPI_Gather(xcvMax, 3, MPI_DOUBLE, xcvMaxArray, 3, MPI_DOUBLE, 0, mpi_comm);

	// Error message: remains empty if no error occurs
	string errorMessage;

	/***********
	 ** Write the data
	 ***********/
	double wtime0, wtimeF;
	if(mpi_rank==0)
		wtime0 = MPI_Wtime();

	// 1. Header
	//      Structure of the header part:
	//              1. step (int)                      2. NcontrolEqns (int)
	//              3. lambda (double*NcontrolEqns)    4. nVars (int)
	//              5. nev (int)                       6. eigenvalues (double*2*nev)
	//              7. mpi_size (int)                  9. cvora (int*(mpi_size+1))
	//              9. xMinArray (double*3*mpi_size)   10. xMaxArray (double*3*mpi_size)
	if(mpi_rank==0) {
		// Open the file
		ofstream ofile;
		ofile.open(filename, ios_base::out | ios_base::trunc | ios_base::binary);

		// Write on the file
		int dummyInt;
		double dummyDouble;

		dummyInt=step;				ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int)); // Note: if you simply use the "<<" operator, it will cause some problem due to bugs in a g++ library on linux
		dummyInt=NcontrolEqns;		ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
			dummyDouble=lambda[iEqn];		ofile.write(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
		}
		dummyInt=nVars;		    	ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		dummyInt=nev;		    	ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		for(int ieig=0; ieig<nev; ++ieig) {
			dummyDouble=directEvalsReal[ieig]; 	ofile.write(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
			dummyDouble=directEvalsImag[ieig]; 	ofile.write(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
		}
		dummyInt=mpi_size; 			ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		for(int irank=0; irank<mpi_size+1; ++irank) {
			dummyInt=cvora[irank]; 		ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		}
		for(int irank=0; irank<mpi_size; ++irank) {
			for(int i=0; i<3; ++i) {
				dummyDouble=xcvMinArray[irank*3+i]; 	ofile.write(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
			}
		}
		for(int irank=0; irank<mpi_size; ++irank) {
			for(int i=0; i<3; ++i) {
				dummyDouble=xcvMaxArray[irank*3+i]; 	ofile.write(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
			}
		}

		// Close the file
		ofile.close();
	}
	MPI_Barrier(mpi_comm);
	int initDisp = sizeof(int)*(5 + mpi_size+1) + sizeof(double)*(NcontrolEqns + nev*2 + 2*mpi_size*3);

	// 2. Body
	//       Structure of the body part:
	//            For each mpi,
	//              1. x_cv (double*3*ncv)
	//              2. direct global modes (double*nev*ncv*nVars*2)
	//              3. adjoint global modes (double*nev*ncv*nVars*2)
	MPI_Status status;
	MPI_File fh;
	MPI_Offset displacement;

	// Open the file
	if (MPI_File_open(mpi_comm, filename, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh)!=0) {
			// Note: MPI_MODE_WRONLY = write only
			//       MPI_MODE_CREATE = create the file if it does not exist.
		stringstream ss;
		ss<<"ERROR! RembrandtWithModels::writeEigenPairsParallel(): Cannot open "<<filename<<endl;
		errorMessage = ss.str();
		showMessageParallel(errorMessage);

		delete [] xcvMinArray; 	delete [] xcvMaxArray;
		throw(REMBRANDT_ERROR_CODE);
	}

	// Write the CV coordinates first
	int myOffsetNcv;
	MPI_Scan(&ncv, &myOffsetNcv, 1, MPI_INT, MPI_SUM, mpi_comm);
	myOffsetNcv -= ncv;
	displacement = MPI_Offset(initDisp + sizeof(double)*myOffsetNcv*(3+nVars*nev*4));
	if(MPI_File_set_view(fh, displacement, MPI_INT, MPI_INT, "native", MPI_INFO_NULL)!=0) {
		// Note: native     = Data in this representation are stored in a file exactly as it is in memory
		//       Internal   = Data in this representation can be used for I/O operations in a homogeneous or heterogeneous environment
		//       External32 = This data representation states that read and write operations convert all data from and to the "external32" representation
		stringstream ss;
		ss<<"ERROR! RembrandtWithModels::writeEigenPairsParallel(): Cannot set the first MPI_File_set_view -- offset="<<displacement<<endl;
		errorMessage = ss.str();
		showMessageParallel(errorMessage);

		delete [] xcvMinArray; 	delete [] xcvMaxArray;
		throw(REMBRANDT_ERROR_CODE);
	}
	double* bufferDouble = new double [ncv*3];
	for(int icv=0; icv<ncv; ++icv) {
		int indexStart = icv*3;
		for(int i=0; i<3; ++i)
			bufferDouble[indexStart+i] = x_cv[icv][i];
	}
	MPI_File_write(fh, bufferDouble, ncv*3, MPI_DOUBLE, &status);
	delete [] bufferDouble; 	bufferDouble = NULL;

	// Write the direct eigen vectors
	displacement += MPI_Offset(3*ncv*sizeof(double));

	bufferDouble = new double [ncv*nVars*2];
	for(int ieig=0; ieig<nev; ++ieig) {
		for(int i=0; i<ncv*nVars; ++i) {
			int index = i*2;
			bufferDouble[index]   = directEvecsReal[ieig][i];
			bufferDouble[index+1] = directEvecsImag[ieig][i];
		}

		MPI_File_write(fh, bufferDouble, ncv*nVars*2, MPI_DOUBLE, &status);
	}

	// Write the adjoint vectors
	displacement += MPI_Offset(nev*ncv*nVars*2*sizeof(double));
	for(int ieig=0; ieig<nev; ++ieig) {
		for(int i=0; i<ncv*nVars; ++i) {
			int index = i*2;
			bufferDouble[index]   = adjointEvecsReal[ieig][i];
			bufferDouble[index+1] = adjointEvecsImag[ieig][i];
		}

		MPI_File_write(fh, bufferDouble, ncv*nVars*2, MPI_DOUBLE, &status);
	}

	delete [] bufferDouble; 	bufferDouble = NULL;

	// 2-4. Close the file
	MPI_File_close(&fh);
	MPI_Barrier(mpi_comm);

	// 3. Foot
	if(mpi_rank==mpi_size-1) {
		ofstream ofile;
		ofile.open(filename, ios_base::out | ios_base::app | ios_base::binary);
		int dummyInt = EOF_ERROR_CHECK_CODE; 	ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		ofile.close();
	}
	MPI_Barrier(mpi_comm);

	if(mpi_rank==0)
		wtimeF = MPI_Wtime();

	/***********
	 ** Show the summary on the screen if the debugging level is high
	 ***********/
	if(debugLevel>0 && mpi_rank==0)
		cout<<endl
		    <<">> "<<nev<<" eigen-pairs are written on "<<filename<<" (RUN TIME FOR WRITING = "<<wtimeF-wtime0<<" [sec]) \n"<<endl;

	/***********
	 ** Free the memory
	 ***********/
	delete [] xcvMinArray;
	delete [] xcvMaxArray;
}

//
// PROTECTED METHODS!!
//

/*
 * Method: calcGlobalModes
 * -----------------------
 * Calculate both direct and adjoint global-modes
 */
int RembrandtWithModels::calcGlobalModes(const int nevSlepc) {
	// Get the Slepc EPS (eigen problem solver) parameters
	SlepcEPSParams slepcEPSParams;
	getSlepcEPSParams(slepcEPSParams, nevSlepc);

	// Get the ST (spectral transformation) parameters for Slepc
	SlepcSTparams slepcSTparams; // By default, SlepcSTparams.useST = false
	getSlepcSTparams(slepcSTparams, slepcEPSParams); // This method gets the parameters only if "SLEPC_ST_PARAM" exists in the input file.
	                                                 // Also, it gives ERROR if SLEPC_ST_PARAM is not specified even though EPSsolverType is EPSGD or EPSJD.

	// Allocate variables
	double* relError = new double [nevSlepc];
	int numIter;

	// On the first time, initialize the petsc solver...
	if (slepcSolver2 == NULL) {
		if (nbocv_v_global == NULL) {
			nbocv_v_global = new int[ncv_gg];
			for (int icv = 0; icv < ncv; icv++)
				nbocv_v_global[icv] = cvora[mpi_rank] + icv;
			updateCvDataG1G2(nbocv_v_global, REPLACE_DATA);
		}

		slepcSolver2 = new SlepcSolver2(cvora, nbocv2_i, nbocv2_v, 5 + nScal, "PCNONE", 0, true);
	}

	// -----------------------------
	// Solve for direct global-modes
	// -----------------------------
	string epsConvHistoryFilename = "SLEPC_HISTORY_DIRECT_MODE.csv";
	slepcSolver2->setConvHistoryFileName(epsConvHistoryFilename);
	bool transposeBeforeSolve = false;
	int nconvDirect = slepcSolver2->solveEigenProblemSlepc<MatComprsedSTL>(directEvalsReal, directEvalsImag, directEvecsReal, directEvecsImag, relError, numIter,
										jacMatrixSTL, cvora, nbocv_v_global, nScal, ncv_gg,
										nevSlepc, slepcEPSParams, slepcSTparams, transposeBeforeSolve);

	// Check the ranges of the eigenvectors
	double minVal_Real[nevSlepc], maxVal_Real[nevSlepc], minVal_Imag[nevSlepc], maxVal_Imag[nevSlepc];
	getEvecsRange(minVal_Real, maxVal_Real, minVal_Imag, maxVal_Imag, directEvecsReal, directEvecsImag, nevSlepc);

	// Report the results (both on the screen and on a file)
	if(mpi_rank==0) {
		cout<<endl
			<<">> "<<nconvDirect<<" converged DIRECT eigen-pairs after "<<numIter<<" iterations"<<endl;
		cout<<"   Result for each eigen-pair:"<<endl;
		for(int i=0; i< std::min<int>(nevSlepc, nconvDirect); ++i)
			printf("    %2d   REL_ERROR = %.5e   EVAL = %.5e + %.5ei  EVEC-RANGE = (real=%.3e~%.3e, imag=%.3e~%.3e)\n",
					i, relError[i], directEvalsReal[i], directEvalsImag[i], minVal_Real[i], maxVal_Real[i], minVal_Imag[i], maxVal_Imag[i]);
		if(nconvDirect == 0)
			cout<<"    NOT AVAILABLE because # of converged eigenpairs = "<<nconvDirect<<endl;
		cout<<endl;
	}

	if(mpi_rank==0 && nconvDirect > 0) {
		FILE *fp = fopen(FILENAME_DIRECT_EVALS, "w");

		fprintf(fp, "NO, EVAL_REAL, EVAL_IMAG\n");

		for(int i=0; i< std::min<int>(nevSlepc, nconvDirect); ++i)
			fprintf(fp, "%d, %.8e, %.8e\n", i, directEvalsReal[i], directEvalsImag[i]);

		// close the file
		fclose(fp);
	}

	MPI_Barrier(mpi_comm);
	delete slepcSolver2; 	slepcSolver2 = NULL;

	// ------------------------------
	// Solve for adjoint global-modes
	// ------------------------------
	slepcSolver2 = new SlepcSolver2(cvora, nbocv2_i, nbocv2_v, 5 + nScal, "PCNONE", 0, true);

//	slepcSolver2->transpose();

	epsConvHistoryFilename = "SLEPC_HISTORY_ADJOINT_MODE.csv";
	slepcSolver2->setConvHistoryFileName(epsConvHistoryFilename);
	transposeBeforeSolve = true;
	int nconvAdjoint = slepcSolver2->solveEigenProblemSlepc<MatComprsedSTL>(adjointEvalsReal, adjointEvalsImag, adjointEvecsReal, adjointEvecsImag, relError, numIter,
			jacMatrixSTL, cvora, nbocv_v_global, nScal, ncv_gg,
			nevSlepc, slepcEPSParams, slepcSTparams, transposeBeforeSolve);

	// We don't need to write the adjoint eigenvalues on a file because they must be the same as the direct eigenvalues.
	// Instead, check if the two set of eigenvalues are the same:
	if(nconvDirect != nconvAdjoint) {
		if(mpi_rank==0) cerr<<"ERROR! nconvDirect(="<<nconvDirect<<") is not equal to nconvAdjoint(="<<nconvAdjoint<<")"<<endl;
		throw(REMBRANDT_ERROR_CODE);
	}

	vector<double> diff_real, diff_imag;
	diff_real.resize(nconvAdjoint, ABSURDLY_BIG_NUMBER);
	diff_imag.resize(nconvAdjoint, ABSURDLY_BIG_NUMBER);
	vector<int> correspDirectIndices;
	correspDirectIndices.resize(nconvDirect, -1);
	for(int i=0; i< nconvAdjoint; ++i) {
		double diff_magSq = ABSURDLY_BIG_NUMBER;

		for(int j=0; j< nconvDirect; ++j) {
			double diff_real_temp = adjointEvalsReal[i] - directEvalsReal[j];
			double diff_imag_temp = adjointEvalsImag[i] - directEvalsImag[j];
			double diff_magSq_temp = pow(diff_real_temp, 2.0) + pow(diff_imag_temp, 2.0);

			if(diff_magSq_temp < diff_magSq) {
				diff_magSq = diff_magSq_temp;
				diff_real[i] = diff_real_temp;
				diff_imag[i] = diff_imag_temp;
				correspDirectIndices[i] = j;
			}
		}
	}

	int equalEvalsCount = 0;
	if(mpi_rank == 0) { cout<<">> Comparison between Direct and Adjoint eigenvalues:"<<endl; };
	for(int i=0; i< nconvAdjoint; ++i) {
		if(fabs(diff_real[i]) < THRESHOLD_EVAL_DIFF || fabs(diff_imag[i]) < THRESHOLD_EVAL_DIFF) {
			++equalEvalsCount;
			if(mpi_rank == 0) { printf("     [%2d] Matched adjoint e-val with the %2dth direct e-val: error = %g + %gi\n", i, correspDirectIndices[i], diff_real[i], diff_imag[i]); }
		} else {
			if(mpi_rank == 0) { printf("     [%2d] NO matched e-val! Closest direct e-val = %2d ----- error = %g + %gi\n", i, correspDirectIndices[i], diff_real[i], diff_imag[i]); }
		}
	}
	if(mpi_rank == 0) { cout<<"   Total "<<equalEvalsCount<<" matched e-vals"<<endl; }

	if(equalEvalsCount == 0) {
		if(mpi_rank==0) cerr<<"ERROR! NO matched e-vals between Direct and Adjoint"<<endl;
		throw(REMBRANDT_ERROR_CODE);
	}

	// Check the ranges of the eigenvectors
	getEvecsRange(minVal_Real, maxVal_Real, minVal_Imag, maxVal_Imag, adjointEvecsReal, adjointEvecsImag, nevSlepc);

	// Report the results (only on the screen)
	if(mpi_rank==0) {
		cout<<endl
		    <<">> "<<nconvAdjoint<<" converged ADJOINT eigen-pairs after "<<numIter<<" iterations"<<endl;
		cout<<"   Result for each eigen-pair:"<<endl;
		for(int i=0; i<std::min<int>(nevSlepc, nconvAdjoint); ++i)
			printf("    %2d   REL_ERROR = %.5e   EVAL = %.5e + %.5ei  EVEC-RANGE = (real=%.3e~%.3e, imag=%.3e~%.3e)\n",
					i, relError[i], adjointEvalsReal[i], adjointEvalsImag[i], minVal_Real[i], maxVal_Real[i], minVal_Imag[i], maxVal_Imag[i]);
		if(nconvAdjoint == 0)
			cout<<"    NOT AVAILABLE because # of converged eigenpairs = "<<nconvAdjoint<<endl;
		cout<<endl;
	}

	// Check the difference of the direct global modes and the adjoint global modes
	if(debugLevel > 1) {
		if(nconvDirect > 0 && nconvAdjoint > 0) {
			double diff1Norm_real[nevSlepc], diff1Norm_imag[nevSlepc], diffInfNorm_real[nevSlepc], diffInfNorm_imag[nevSlepc];
			checkTwoEvecsDifferenceRange(diff1Norm_real, diff1Norm_imag, diffInfNorm_real, diffInfNorm_imag,
					directEvecsReal, directEvecsImag, adjointEvecsReal, adjointEvecsImag, nevSlepc);
			if(mpi_rank==0) {
				cout<<endl
						<<">> Difference between the direct and the adjoint modes:"<<endl;
				for(int i=0; i<std::min<int>(nevSlepc, nconvAdjoint); ++i)
					printf("    %2d   DIFF_REAL = (1-norm=%.5e, inf-norm=%.5e)   DIFF_IMAG = (1-norm=%.5e, inf-norm=%.5e)\n",
							i, diff1Norm_real[i], diffInfNorm_real[i], diff1Norm_imag[i], diffInfNorm_imag[i]);
				cout<<endl;
			}
		}
	}

	MPI_Barrier(mpi_comm);

	delete slepcSolver2; 	slepcSolver2 = NULL;

	// Free memory
	delete [] relError;

	// Return
//		return std::min<int>(std::min<int>(nconvDirect, nconvAdjoint), nevSlepc); // Note: Many times nconv >= nevSlepc, which will generate segmentation falut
	return std::min<int>(nconvDirect, nevSlepc); // Note: Many times nconv >= nevSlepc, which will generate segmentation falut
}

/*
 * Method: getSlepcEPSParams
 * -------------------------
 * Get SLEPc EPS (eigen problem solver parameters) from the input file
 *
 * Arguments:
 *   nevSlepc = This is required only for printing output
 */
void RembrandtWithModels::getSlepcEPSParams(SlepcEPSParams &slepcEPSParams, const int nevSlepc) {
	// Get the slepc parameters
	if (!checkParam("SLEPC_EIGEN_PARAM")) {
		ParamMap::add("SLEPC_EIGEN_PARAM  NCV=15  MPD=15  TOL=1.0e-7  MAX_ITER=1000"); // add default values
		if (mpi_rank == 0)
			cout<< "WARNING: added keyword \"SLEPC_EIGEN_PARAM  NCV=15  MPD=15  TOL=1.0e-7  MAX_ITER=1000\""<< " to parameter map!" << endl;
	}
	slepcEPSParams.ncvSlepc = getParam("SLEPC_EIGEN_PARAM")->getInt("NCV");
	slepcEPSParams.mpdSlepc = getParam("SLEPC_EIGEN_PARAM")->getInt("MPD");

	slepcEPSParams.stringWhichEigenOfInterest = getParam("SLEPC_EIGEN_PARAM")->getString("EIGEN_OF_INTEREST");
	slepcEPSParams.WhichEigenOfInterest       = getEigenOfInterestFromString(slepcEPSParams.stringWhichEigenOfInterest);

	slepcEPSParams.stringEPSsolverType = getParam("SLEPC_EIGEN_PARAM")->getString("EPS_SOLVER_TYPE");
	slepcEPSParams.EPSsolverType       = getEPSsolverTypeFromString(slepcEPSParams.stringEPSsolverType);

	slepcEPSParams.tolSlepc      = getParam("SLEPC_EIGEN_PARAM")->getDouble("TOL");
	slepcEPSParams.max_iterSlepc = getParam("SLEPC_EIGEN_PARAM")->getInt("MAX_ITER");

	slepcEPSParams.EPSmonitorInterv = getIntParam("SLEPC_EPS_CHECK_INTERVAL", "100");

	slepcEPSParams.targetValue = 0.0; // If "EPS_TARGET_MAGNITUDE", "EPS_TARGET_REAL", or "EPS_TARGET_IMAGINARY" is used, set the target-value.
	if(slepcEPSParams.WhichEigenOfInterest == EPS_TARGET_MAGNITUDE || slepcEPSParams.WhichEigenOfInterest == EPS_TARGET_REAL
			|| slepcEPSParams.WhichEigenOfInterest == EPS_TARGET_IMAGINARY)
		slepcEPSParams.targetValue = getParam("SLEPC_EIGEN_PARAM")->getDouble("TARGET_VALUE");

	// Show the result on the screen
	if(mpi_rank==0) {
		slepcEPSParams.showParamOnScreen(nevSlepc);
	}
}

/*
 * Method: getEigenOfInterestFromString
 * ------------------------------------
 *
 */
EPSWhich RembrandtWithModels::getEigenOfInterestFromString(string& tempString) {
	// Available selection: EPS_LARGEST_MAGNITUDE   Largest |lambda|
	//                      EPS_SMALLEST_MAGNITUDE  Smallest |lambda|
	//                      EPS_LARGEST_REAL        Largest Re(lambda)
	//                      EPS_SMALLEST_REAL       Smallest Re(lambda)
	//                      EPS_LARGEST_IMAGINARY   Largest Im(lambda)
	//                      EPS_SMALLEST_IMAGINARY  Smallest Im(lambda)
	//                      EPS_TARGET_MAGNITUDE    Smallest |lambda - tau|
	//                      EPS_TARGET_REAL         Smallest |Re(lambda - tau)|
	//                      EPS_TARGET_IMAGINARY    Smallest |Im(lambda - tau)|
	//                      EPS_ALL                 All lambda in [a, b]
	//                      EPS_WHICH_USER          user-defined

	if(tempString.compare("EPS_LARGEST_MAGNITUDE") == 0)
		return EPS_LARGEST_MAGNITUDE;
	else if(tempString.compare("EPS_SMALLEST_MAGNITUDE") == 0)
		return EPS_SMALLEST_MAGNITUDE;
	else if(tempString.compare("EPS_LARGEST_REAL") == 0)
		return EPS_LARGEST_REAL;
	else if(tempString.compare("EPS_SMALLEST_REAL") == 0)
		return EPS_SMALLEST_REAL;
	else if(tempString.compare("EPS_LARGEST_IMAGINARY") == 0)
		return EPS_LARGEST_IMAGINARY;
	else if(tempString.compare("EPS_SMALLEST_IMAGINARY") == 0)
		return EPS_SMALLEST_IMAGINARY;
	else if(tempString.compare("EPS_TARGET_MAGNITUDE") == 0)
		return EPS_TARGET_MAGNITUDE;
	else if(tempString.compare("EPS_TARGET_REAL") == 0)
		return EPS_TARGET_REAL;
	else if(tempString.compare("EPS_TARGET_IMAGINARY") == 0)
		return EPS_TARGET_IMAGINARY;
	else if(tempString.compare("EPS_ALL") == 0)
		return EPS_ALL;
	else {
		if(mpi_rank == 0) {
			cerr<<"ERROR getEigenOfInterestFromString(): EIGEN_OF_INTEREST = "<<tempString<<" is not in the list"<<endl;
		}
		throw(REMBRANDT_ERROR_CODE);
	}
}

/*
 * Method: getStringFromEigenOfInterest
 * ------------------------------------
 *
 */
string RembrandtWithModels::getStringFromEigenOfInterest(const EPSWhich eigenOfInterest) {
	// Available selection: EPS_LARGEST_MAGNITUDE   Largest |lambda|
	//                      EPS_SMALLEST_MAGNITUDE  Smallest |lambda|
	//                      EPS_LARGEST_REAL        Largest Re(lambda)
	//                      EPS_SMALLEST_REAL       Smallest Re(lambda)
	//                      EPS_LARGEST_IMAGINARY   Largest Im(lambda)
	//                      EPS_SMALLEST_IMAGINARY  Smallest Im(lambda)
	//                      EPS_TARGET_MAGNITUDE    Smallest |lambda - tau|
	//                      EPS_TARGET_REAL         Smallest |Re(lambda - tau)|
	//                      EPS_TARGET_IMAGINARY    Smallest |Im(lambda - tau)|
	//                      EPS_ALL                 All lambda in [a, b]
	//                      EPS_WHICH_USER          user-defined

	string tempString;

	if(eigenOfInterest == EPS_LARGEST_MAGNITUDE)
		tempString = "EPS_LARGEST_MAGNITUDE";
	else if(eigenOfInterest == EPS_SMALLEST_MAGNITUDE)
		tempString = "EPS_SMALLEST_MAGNITUDE";
	else if(eigenOfInterest == EPS_LARGEST_REAL)
		tempString = "EPS_LARGEST_REAL";
	else if(eigenOfInterest == EPS_SMALLEST_REAL)
		tempString = "EPS_SMALLEST_REAL";
	else if(eigenOfInterest == EPS_LARGEST_IMAGINARY)
		tempString = "EPS_LARGEST_IMAGINARY";
	else if(eigenOfInterest == EPS_SMALLEST_IMAGINARY)
		tempString = "EPS_SMALLEST_IMAGINARY";
	else if(eigenOfInterest == EPS_TARGET_MAGNITUDE)
		tempString = "EPS_TARGET_MAGNITUDE";
	else if(eigenOfInterest == EPS_TARGET_REAL)
		tempString = "EPS_TARGET_REAL";
	else if(eigenOfInterest == EPS_TARGET_IMAGINARY)
		tempString = "EPS_TARGET_IMAGINARY";
	else if(eigenOfInterest == EPS_ALL)
		tempString = "EPS_ALL";
	else {
		if(mpi_rank == 0) {
			cerr<<"ERROR getStringFromEigenOfInterest(): EIGEN_OF_INTEREST = "<<eigenOfInterest<<" is not in the list"<<endl;
		}
		throw(REMBRANDT_ERROR_CODE);
	}

	return tempString;
}

/*
 * Method: getEPSsolverTypeFromString
 * ----------------------------------
 *
 */
EPSType RembrandtWithModels::getEPSsolverTypeFromString(string &tempString) {
	// Available solvers: Power/Inverse/RQI   :  EPSPOWER       (largest |lambda|;    Any problem type; Complex)
	//                    Subspace Iteration  :  EPSSUBSPACE    (largest |lambda|;    Any problem type; Complex)
	//                    Arnoldi             :  EPSARNOLDI     (Any spectrum;        Any problem type; Complex)
	//                    Lanczos             :  EPSLANCZOS     (Any spectrum;        EPS_HEP,EPS_GHEP; Complex)
	//                    Krylov-Schur        :  EPSKRYLOVSCHUR (Any spectrum;        Any problem type; Complex)
	//                    Generalized Davidson:  EPSGD          (Any spectrum;        Any problem type; Complex)
	//                    Jacobi-Davidson     :  EPSJD          (Any spectrum;        Any problem type; Complex)
	//                    Rayleigh quotient CG:  EPSRQCG        (smallest Re(lambda); EPS_HEP,EPS_GHEP; Complex)
	//                    --------------------------------------
	//                    lapack solver     :  EPSLAPACK    (Any spectrum;                  Any problem type; Complex)
	//                    Wrapper to arpack :  EPSARPACK    (Any spectrum;                  Any problem type; Complex)
	//                    Wrapper to primme :  EPSPRIMME    (largest & smallest Re(lambda); EPS_HEP;          Complex)
	//                    Wrapper to blzpack:  EPSBLZPACK   (smallest Re(lambda);           EPS_HEP,EPS_GHEP; no)
	//                    Wrapper to trlan  :  EPSTRLAN     (largest & smallest Re(lambda); EPS_HEP;          no)
	//                    Wrapper to blopex :  EPSBLOPEX    (smallest Re(lambda);           EPS_HEP,EPS_GHEP; Complex)

	if(tempString.compare("EPSPOWER") == 0)
		return EPSPOWER;
	else if(tempString.compare("EPSSUBSPACE") == 0)
		return EPSSUBSPACE;
	else if(tempString.compare("EPSARNOLDI") == 0)
		return EPSARNOLDI;
	else if(tempString.compare("EPSLANCZOS") == 0)
		return EPSLANCZOS;
	else if(tempString.compare("EPSKRYLOVSCHUR") == 0)
		return EPSKRYLOVSCHUR;
	else if(tempString.compare("EPSGD") == 0)
		return EPSGD;
	else if(tempString.compare("EPSJD") == 0)
		return EPSJD;
//		else if(tempString.compare("EPSRQCG") == 0)  // For some reason, EPSRQCG is not supported by SLEPc
//			return EPSRQCG;
	else if(tempString.compare("EPSLAPACK") == 0)
		return EPSLAPACK;
	else if(tempString.compare("EPSARPACK") == 0)
		return EPSARPACK;
	else if(tempString.compare("EPSPRIMME") == 0)
		return EPSPRIMME;
	else if(tempString.compare("EPSBLZPACK") == 0)
		return EPSBLZPACK;
	else if(tempString.compare("EPSTRLAN") == 0)
		return EPSTRLAN;
	else if(tempString.compare("EPSBLOPEX") == 0)
		return EPSBLOPEX;
	else {
		if(mpi_rank == 0) {
			cerr<<"ERROR getEPSsolverTypeFromString(): EPS_SOLVER_TYPE = "<<tempString<<" is not in the list"<<endl;
		}
		throw(REMBRANDT_ERROR_CODE);
	}
}

/*
 * Method: getSlepcSTparams
 * ------------------------
 * get the spectral transformation (ST) parameters for Slepc from the input file
 *
 * Arguments:
 *   EPSsolverType is required for error checks: For EPSGD or EPSJD, you must use ST.
 *
 * Note: targetValue is not set here -- You must set it separately.
 */
void RembrandtWithModels::getSlepcSTparams(SlepcSTparams& slepcSTparams, const SlepcEPSParams &slepcEPSParams) {
	if (!checkParam("SLEPC_ST_PARAM")) {
		if(slepcEPSParams.EPSsolverType == EPSGD || slepcEPSParams.EPSsolverType == EPSJD) {
			if(mpi_rank == 0)
				cerr<<"ERROR! For EPSGD or EPSJD, you must use ST, but \"SLEPC_ST_PARAM\" cannot be found in the input file."<<endl
				<<"       It should be in the form of \"SLEPC_ST_PARAM  USE_ST=YES  ST_TYPE=STPRECOND  KSP_TYPE=GMRES  PC_METHOD=ASM  ST_SHIFT=0.0\""<<endl;
			throw(REMBRANDT_ERROR_CODE);
		}

		if(mpi_rank == 0)
			cerr<<"WARNING! \"SLEPC_ST_PARAM\" cannot be found in the input file"<<endl
			<<"         If you want to use ST, add \"SLEPC_ST_PARAM  USE_ST=YES  ST_TYPE=STSINVERT  KSP_TYPE=PREONLY  PC_METHOD=LU  ST_SHIFT=0.0\""<<endl
			<<endl;
	} else {
		string tempStringST = getParam("SLEPC_ST_PARAM")->getString("USE_ST");
		std::transform(tempStringST.begin(), tempStringST.end(), tempStringST.begin(), ::tolower);
		if(tempStringST.compare("true")==0 || tempStringST.compare("yes")==0) {
			slepcSTparams.useST = true;

			// sttype
			string tempStringSTParam = getParam("SLEPC_ST_PARAM")->getString("ST_TYPE");
			// Available ST types:
			//   Spectral Transformation | STType    | Operator
			//   -----------------------------------------------------------------------
			//   Shift of Origin         | STSHIFT   | B^(-1)*A + \sigma*I
			//   Spectrum Folding        | STFOLD    | (B^(-1)*A - \sigma*I)^2
			//   Shift-and-invert        | STSINVERT | (A - \sigma*B)^(-1) * B
			//   Generalized Cayley      | STCAYLEY  | (A - \sigma*B)^(-1) * (A + \nu*B)
			//   Preconditioner          | STPRECOND | K^(-1) ~ (A - \sigma*B)^(-1)
			//   -----------------------------------------------------------------------
			//   Shell Transformation    | STSHELL   | user-defined
			// Default = STPRECOND
			if(tempStringSTParam.compare("SHIFT") == 0)
				slepcSTparams.sttype_ = STSHIFT;
			else if(tempStringSTParam.compare("FOLD") == 0)
				slepcSTparams.sttype_ = STFOLD;
			else if(tempStringSTParam.compare("SINVERT") == 0)
				slepcSTparams.sttype_ = STSINVERT;
			else if(tempStringSTParam.compare("CAYLEY") == 0)
				slepcSTparams.sttype_ = STCAYLEY;
			else
				slepcSTparams.sttype_ = STPRECOND; // This is the default

			// ksptype
			tempStringSTParam = getParam("SLEPC_ST_PARAM")->getString("KSP_TYPE");

			if(tempStringSTParam.compare("GMRES") == 0)
				slepcSTparams.ksptype_ = KSPGMRES;
			else if(tempStringSTParam.compare("PREONLY") == 0)  // You need this to use LU decompostion
				slepcSTparams.ksptype_ = KSPPREONLY;
			else if(tempStringSTParam.compare("BICG") == 0)
				slepcSTparams.ksptype_ = KSPBICG;
			else {
				if(mpi_rank == 0)
					cerr<<"ERROR! SLEPC_ST_PARAM->KSP_TYPE = "<<tempStringSTParam<<" is not supported in REMBRANDT"<<endl;
				throw(REMBRANDT_ERROR_CODE);
			}

			// pcMethod
			tempStringSTParam = getParam("SLEPC_ST_PARAM")->getString("PC_METHOD");

			if(tempStringSTParam.compare("LU") == 0)
				slepcSTparams.pcMethod = PCLU;
			else if(tempStringSTParam.compare("BJACOBI") == 0) // Block Jacobi
				slepcSTparams.pcMethod = PCBJACOBI;
			else if(tempStringSTParam.compare("ASM") == 0) // Additive Schwarz
				slepcSTparams.pcMethod = PCASM;
			else if(tempStringSTParam.compare("ILU") == 0) {
				if(mpi_size==1)
					slepcSTparams.pcMethod = PCILU;   // Incomplete LU -- The ILU preconditioner in Petsc only works for serial -- For parallel, call an external package
				else
					slepcSTparams.pcMethod = PCHYPRE; // Incomplete LU -- For parallel, use PCHYPRE with the PCHYPRESetType(pc, "euclid") command
			}
			else {
				if(mpi_rank == 0)
					cerr<<"ERROR! SLEPC_ST_PARAM->PC_METHOD = "<<tempStringSTParam<<" is not supported in REMBRANDT"<<endl;
				throw(REMBRANDT_ERROR_CODE);
			}

			// STshift and nu
			slepcSTparams.STshift = getParam("SLEPC_ST_PARAM")->getDouble("ST_SHIFT");  // Sometimes, STshift is called "sigma".
			// Often, people use the same value as EPS_TARGET for STshift.
			if(slepcSTparams.sttype_ == STCAYLEY)
				slepcSTparams.nu = getParam("SLEPC_ST_PARAM")->getDouble("ST_NU"); // Note: "nu" is only required for Generalized Cayley

			//			if(mpi_rank == 0)
			//				slepcSTparams.showParamOnScreen();
		} else {
			if(slepcEPSParams.EPSsolverType == EPSGD || slepcEPSParams.EPSsolverType == EPSJD) {
				cerr<<"ERROR! For EPSGD or EPSJD, you must use ST"<<endl;
				throw(REMBRANDT_ERROR_CODE);
			}
		}
	}
}

/*
 * Method: getEvecRange
 * --------------------
 * Get the ranges of eigen-vectors
 * Return = minVal_Real[], maxVal_Real[], minVal_Imag[], maxVal_Imag[]
 */
void RembrandtWithModels::getEvecsRange(double* minVal_Real, double* maxVal_Real, double* minVal_Imag, double* maxVal_Imag,
		double** EvecsReal, double** EvecsImag, const int nevSlepc) {
	double myMinVal_Real[nevSlepc], myMaxVal_Real[nevSlepc], myMinVal_Imag[nevSlepc], myMaxVal_Imag[nevSlepc];
	for(int iev=0; iev<nevSlepc; ++iev) {
		myMinVal_Real[iev] =  ABSURDLY_BIG_NUMBER;
		myMaxVal_Real[iev] = -ABSURDLY_BIG_NUMBER;
		myMinVal_Imag[iev] =  ABSURDLY_BIG_NUMBER;
		myMaxVal_Imag[iev] = -ABSURDLY_BIG_NUMBER;
	}
	for(int iev=0; iev<nevSlepc; ++iev) {
		for(int localIndex=0; localIndex<(cvora[mpi_rank+1]-cvora[mpi_rank])*(5+nScal); ++localIndex) {
			myMinVal_Real[iev] = min(myMinVal_Real[iev], EvecsReal[iev][localIndex]);
			myMaxVal_Real[iev] = max(myMaxVal_Real[iev], EvecsReal[iev][localIndex]);
			myMinVal_Imag[iev] = min(myMinVal_Imag[iev], EvecsImag[iev][localIndex]);
			myMaxVal_Imag[iev] = max(myMaxVal_Imag[iev], EvecsImag[iev][localIndex]);
		}
	}
	MPI_Allreduce(myMinVal_Real, minVal_Real, nevSlepc, MPI_DOUBLE, MPI_MIN, mpi_comm);
	MPI_Allreduce(myMaxVal_Real, maxVal_Real, nevSlepc, MPI_DOUBLE, MPI_MAX, mpi_comm);
	MPI_Allreduce(myMinVal_Imag, minVal_Imag, nevSlepc, MPI_DOUBLE, MPI_MIN, mpi_comm);
	MPI_Allreduce(myMaxVal_Imag, maxVal_Imag, nevSlepc, MPI_DOUBLE, MPI_MAX, mpi_comm);
}

/*
 * Method: checkTwoEvecsDifferenceRange
 * ------------------------------------
 * Get the ranges of the difference between two eigen-vectors
 * Return = 1-norm of real and imaginary parts
 */
void RembrandtWithModels::checkTwoEvecsDifferenceRange(double* diff1Norm_real, double* diff1Norm_imag, double* diffInfNorm_real, double* diffInfNorm_imag,
		double** EvecsReal_1st, double** EvecsImag_1st, double** EvecsReal_2nd, double** EvecsImag_2nd,
		const int nevSlepc) {
	double myDiff1Norm_real[nevSlepc], myDiff1Norm_imag[nevSlepc];
	double myDiffInfNorm_real[nevSlepc], myDiffInfNorm_imag[nevSlepc];
	for(int iev=0; iev<nevSlepc; ++iev) {
		myDiff1Norm_real[iev] = 0.0;
		myDiff1Norm_imag[iev] = 0.0;
		myDiffInfNorm_real[iev] = -ABSURDLY_BIG_NUMBER;
		myDiffInfNorm_imag[iev] = -ABSURDLY_BIG_NUMBER;
	}
	for(int iev=0; iev<nevSlepc; ++iev) {
		for(int localIndex=0; localIndex<(cvora[mpi_rank+1]-cvora[mpi_rank])*(5+nScal); ++localIndex) {
			double abs_diff_real = fabs(EvecsReal_1st[iev][localIndex] - EvecsReal_2nd[iev][localIndex]);
			double abs_diff_imag = fabs(EvecsImag_1st[iev][localIndex] - EvecsImag_2nd[iev][localIndex]);

			myDiff1Norm_real[iev] += abs_diff_real;
			myDiff1Norm_imag[iev] += abs_diff_imag;
			myDiffInfNorm_real[iev] = max(myDiffInfNorm_real[iev], abs_diff_real);
			myDiffInfNorm_imag[iev] = max(myDiffInfNorm_imag[iev], abs_diff_imag);
		}
	}

	MPI_Allreduce(myDiff1Norm_real, diff1Norm_real, nevSlepc, MPI_DOUBLE, MPI_SUM, mpi_comm);
	MPI_Allreduce(myDiff1Norm_imag, diff1Norm_imag, nevSlepc, MPI_DOUBLE, MPI_SUM, mpi_comm);
	MPI_Allreduce(myDiffInfNorm_real, diffInfNorm_real, nevSlepc, MPI_DOUBLE, MPI_MAX, mpi_comm);
	MPI_Allreduce(myDiffInfNorm_imag, diffInfNorm_imag, nevSlepc, MPI_DOUBLE, MPI_MAX, mpi_comm);
}

