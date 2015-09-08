#include "vanGoghWithModels.h"

/*
 * Method: init
 * ------------
 * Initialize VanGoghWithModels member variables
 */
void VanGoghWithModels::init() {
	// initial field (from a restart file)
	rho_init = NULL; 	registerScalar(rho_init,  "INIT_RHO", CV_DATA);
	rhou_init = NULL; 	registerVector(rhou_init, "INIT_RHOU", CV_DATA);
	rhoE_init = NULL; 	registerScalar(rhoE_init, "INIT_RHOE", CV_DATA);
	// NOTE: scalars should be registered in the child class

	// data on the unstable branch (from Q0_PT000**.bin)
	rho_unstable = NULL; 	registerScalar(rho_unstable,  "UNSTABLE_RHO", CV_DATA);
	rhou_unstable = NULL;	registerVector(rhou_unstable, "UNSTABLE_RHOU", CV_DATA);
	rhoE_unstable = NULL; 	registerScalar(rhoE_unstable, "UNSTABLE_RHOE", CV_DATA);
	// NOTE: scalars should be registered in the child class

	// -------------------------
	// Perturbations
	// -------------------------
	// perturbation
	array_perturb = NULL; 	registerScalar(array_perturb, "ARRAY_PERTURB", CV_DATA);

	rho_perturb  = NULL; 	registerScalar(rho_perturb,  "PERTURB_RHO", CV_DATA);
	rhou_perturb = NULL;	registerVector(rhou_perturb, "PERTURB_RHOU", CV_DATA);
	rhoE_perturb = NULL; 	registerScalar(rhoE_perturb, "PERTURB_RHOE", CV_DATA);

	wdFace = NULL;

	// -------------------------
	// Eigenvecs: Tecplot output
	// -------------------------
	// least-stable global modes
	rho_Dct_1stReal  = NULL; 	registerScalar(rho_Dct_1stReal,  "DCT1ST_REAL_RHO", CV_DATA);
	rhou_Dct_1stReal = NULL;	registerVector(rhou_Dct_1stReal, "DCT1ST_REAL_RHOU", CV_DATA);
	rhoE_Dct_1stReal = NULL;	registerScalar(rhoE_Dct_1stReal, "DCT1ST_REAL_RHOE", CV_DATA);
	// NOTE: scalars should be registered in the child class

	// first adjoint global modes
	rho_Adj_1stReal  = NULL;	registerScalar(rho_Adj_1stReal,  "ADJ1ST_REAL_RHO", CV_DATA);
	rhou_Adj_1stReal = NULL;	registerVector(rhou_Adj_1stReal, "ADJ1ST_REAL_RHOU", CV_DATA);
	rhoE_Adj_1stReal = NULL;	registerScalar(rhoE_Adj_1stReal, "ADJ1ST_REAL_RHOE", CV_DATA);
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

	if(wdFace != NULL) {
		delete [] wdFace;         wdFace = NULL;
	}

	Evals.clear();
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

	// Get the parameters for the perturbations
	getPerturbParams(perturbParams);

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

	// Update xcvMin and xcvMax
	for(int icv=0; icv<ncv; ++icv) {
		for(int i=0; i<3; ++i) {
			xcvMin[i] = min(xcvMin[i], x_cv[icv][i]);
			xcvMax[i] = max(xcvMax[i], x_cv[icv][i]);
		}
	}

	int nScal = scalarTranspEqVector.size();
	
	double* metricNumer = new double [5+nScal];
	double* metricDenom = new double [5+nScal];

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
		readEigenPairsRembrandtSerial(filename, ptOnSCurve, lambda_file, NcontrolEqns); // Read Rembrandt-dump file and update Evals, DirectEvecs, AdjointEvecs, etc.

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
	updateEigenVecTecplotScalarRansTurbModel();
	updateEigenVecTecplotScalarRansCombModel();

	// =====================================
	// Read the point on the unstable branch
	// =====================================
	string filenameQ1 = getStringParam("UNSTABLE_Q1_FILENAME", "NONE");
	double* qVecTemp = new double [ncv*(5+nScal)];
	if(mpi_rank==0) cout<<"> UNSTABLE_Q1_FILENAME = "<<filenameQ1<<endl;
	readPsALCdumpedDataSerial(filenameQ1.c_str(), qVecTemp, NcontrolEqns);

	updateUnstableVecNS(qVecTemp, 5+nScal);
	updateUnstableVecRansTurbModel(qVecTemp, 5+nScal);
	updateUnstableVecRansCombModel(qVecTemp, 5+nScal);

	MPI_Barrier(mpi_comm);
	delete [] qVecTemp; 	qVecTemp = NULL;

	// ===================================
	// Check the initial unperturbed field
	// ===================================
	// It becomes a bit tricky if ALREADY_PERTURBED_RESTART == YES!
	string boolStringInit = getStringParam("ALREADY_PERTURBED_RESTART", "NO");
	if(mpi_rank==0) cout<<"> ALREADY_PERTURBED_RESTART = "<<boolStringInit<<endl;
	std::transform(boolStringInit.begin(), boolStringInit.end(), boolStringInit.begin(), ::toupper);
	bool isAlreadyPerturbeRestart = false;
	if(boolStringInit.compare("YES")==0 || boolStringInit.compare("Y")==0 || boolStringInit.compare("TRUE")==0 || boolStringInit.compare("T")==0) 
		isAlreadyPerturbeRestart = true;
	
	// Calculate QoI and optimal metric
	if(!isAlreadyPerturbeRestart) {
		itest = 0; // In general, "itest=0" means the initial unperturbed field, whereas "itest=1" is the first perturbed field.
		
		char filenameQoI[15];
		sprintf(filenameQoI, "QoI.txt");
		writeQoIOnFile(itest, filenameQoI, true);
		
		if(perturbParams.calcOptimalMetric) {
			calcOptMetrics(metricNumer, metricDenom);
			
			if(mpi_rank == 0) {
				cout<<"> Optimal metric of the init unperturbed field = ";
				for(int i=0; i<5+nScal-1; ++i) {
					if(i==3) // rhow can be zero
						cout<<metricNumer[i]/(metricDenom[i] + MACHINE_EPS)<<" | ";
					else 
						cout<<metricNumer[i]/metricDenom[i]<<" | ";
				}
				cout<<metricNumer[5+nScal-1]/metricDenom[5+nScal-1]<<endl;
			}
			
			writeOptimalMeticOnFile(itest, OPTIMAL_METRIC_INIT_FILENAME, true, metricNumer, metricDenom, nScal);
		}
	} 
	
	// Save the initial data from the restart file if the user wants to reinitialize the flow field by random perturbation
	if(mpi_rank==0)
		cout<<"> Store the initial field from restart"<<endl;
	storeInitRestartNS();
	storeInitRestartScalarRansTurbModel();
	storeInitRestartScalarRansCombModel();
	
	// =====================
	// Perturbation analysis
	// =====================
	if(isAlreadyPerturbeRestart) {
		itest = getIntParam("TEST_NUM"); 	assert(itest>0);
		ntests = itest;  // You can launch not multiple perturbations but a single perturbation 
		                 // because the initial unperturbed field is not provided.  
	} else {
		itest = 1; // In general, "itest=0" means the initial unperturbed field, whereas "itest=1" is the first perturbed field. 
		ntests = getIntParam("NTESTS", "1");
	}

	// Allocate memory: wdFace
	if(perturbParams.useFiltering) {
		if(mpi_rank==0)
			cout<<"> Allocate memory for wdFace"<<endl;
		assert(wdFace == NULL);  	wdFace = new double [nfa_b2gg]; // wall distance at faces

		calcWallDistanceFace(wdFace);
	}

	if(mpi_rank==0)
		cout<<endl
		    <<"-------------------"<<endl
		    <<" Perturbation test "<<endl
		    <<"-------------------"<<endl;

	// set the stopping criteria for time-advancement : Note that "resid_energ_th" is a member variable in UgpWithCvCompFlow.h
	resid_energ_th = getDoubleParam("RESID_ENERG_TH", 1.0e-20);

	bool firstItest = true;
	// run simulations
	while(itest <= ntests) {
		if(mpi_rank==0)
			cout<<endl
				<<"> RUNNING THE "<<itest<<"TH PERTURBATION TEST"<<endl
				<<endl;

		if(!isAlreadyPerturbeRestart)
			step = 0;

		if(!firstItest) { // Note: We don't need to reinitialize the flow field for the first itest.
			if(mpi_rank==0)
				cout<<"> Note: Reinitialize the flow field with the stored initial field"<<endl;

			// reinitialize models
			reinitialHookScalarRansTurbModel();
			reinitialHookScalarRansCombModel();
			for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
				updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);

			// reinitialize Navier-Stokes equations
			reinitialHook();

			updateCvDataG1G2(rho, REPLACE_DATA);
			updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
			updateCvDataG1G2(rhoE, REPLACE_DATA);
		}

		// -------------------------------------------------
		// random perturbation
		// -------------------------------------------------
		// Perturb the field with filtering
		perturbFieldScalarRansTurbModel();
		perturbFieldScalarRansCombModel();
		perturbFieldNS();

		updateCvDataG1G2(rho,  REPLACE_DATA);
		updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
		updateCvDataG1G2(rhoE, REPLACE_DATA);
		for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
			updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);

		// -------------------------------------------------
		// calculate the optimal metric
		// -------------------------------------------------
		if(perturbParams.calcOptimalMetric) {
			calcOptMetrics(metricNumer, metricDenom);
			if(mpi_rank == 0) {
				cout<<"> Optimal metric of the perturbed field = ";
				for(int i=0; i<5+nScal-1; ++i) {
					if(i==3) // rhow can be zero
						cout<<metricNumer[i]/(metricDenom[i] + MACHINE_EPS)<<" | ";
					else 
						cout<<metricNumer[i]/metricDenom[i]<<" | ";
				}
				cout<<metricNumer[5+nScal-1]/metricDenom[5+nScal-1]<<endl;
			}
			
			if(isAlreadyPerturbeRestart)
				writeOptimalMeticOnFile(itest, OPTIMAL_METRIC_INIT_FILENAME, firstItest, metricNumer, metricDenom, nScal);
			else
				writeOptimalMeticOnFile(itest, OPTIMAL_METRIC_INIT_FILENAME, false, metricNumer, metricDenom, nScal);
		}
		
		// -------------------------------------------------
		// run the solver
		// -------------------------------------------------
		// calculate state variables: press, etc.
		calcStateVariables();

		char filenameResidual[40];
		sprintf(filenameResidual, "residHist.case%04d.csv", itest);
		writeResidualOnFile(itest, filenameResidual, true);

		if(!isAlreadyPerturbeRestart)
			writeData2(itest, 0, true); // Save each case on a different file

		// write restart file
		char restartInitFilename[40];
		sprintf(restartInitFilename, "Restart.case%04d.init", itest);
		writeRestart(restartInitFilename);

		// run solver
		string tIntName = getStringParam("TIME_INTEGRATION");
		{
			if (tIntName == "FORWARD_EULER")                 runForwardEuler();
			else if (tIntName == "RK")                       runRK();
			else if (tIntName == "BACKWARD_EULER")           runBackwardEuler();
			else if (tIntName == "BACKWARD_EULER_COUPLED")   runBackwardEulerCoupled();
			else if (tIntName == "BDF2")                     runBDF2();
			else {
				cerr << "ERROR: wrong time integration scheme specified !" << endl;
				cerr << "available integration schemes are: FORWARD_EULER, RK, BACKWARD_EULER, BDF2" << endl;
			}
		}

		// Post-processing
		char filenameQoI[15];
		sprintf(filenameQoI, "QoI.txt");
		if(isAlreadyPerturbeRestart)
			writeQoIOnFile(itest, filenameQoI, firstItest);
		else
			writeQoIOnFile(itest, filenameQoI, false);

		if(perturbParams.calcOptimalMetric) {
			calcOptMetrics(metricNumer, metricDenom);
			if(mpi_rank == 0) {
				cout<<"> Optimal metric of the perturbed field = ";
				for(int i=0; i<5+nScal-1; ++i) {
					if(i==3) // rhow can be zero
						cout<<metricNumer[i]/(metricDenom[i] + MACHINE_EPS)<<" | ";
					else 
						cout<<metricNumer[i]/metricDenom[i]<<" | ";
				}
				cout<<metricNumer[5+nScal-1]/metricDenom[5+nScal-1]<<endl;
			}
			
			writeOptimalMeticOnFile(itest, OPTIMAL_METRIC_FINAL_FILENAME, firstItest, metricNumer, metricDenom, nScal);
		}
		
		writeData2(itest, step, true); // Save each case on a different file

		char restartFilename[40];
		sprintf(restartFilename, "Restart.test%04d.%06d.out", itest, step);
		writeRestart(restartFilename);

		++itest;
		firstItest = false;
	}
	
	if(metricNumer != NULL) { delete [] metricNumer; 	metricNumer = NULL; }
	if(metricDenom != NULL) { delete [] metricDenom; 	metricDenom = NULL; }
}

/*
 * Method: temporalHook
 * --------------------
 *
 */
void VanGoghWithModels::temporalHook() {
	// Write residual history on a file in every "CHECK_INTERVAL"
	int check_interval = getIntParam("CHECK_INTERVAL", "1");
	if(step > 0) {
		if(step == 1 || step % check_interval == 0) {
			char filenameResidual[40];
			sprintf(filenameResidual, "residHist.case%04d.csv", itest);
			writeResidualOnFile(itest, filenameResidual, false);
		}
	}

	// Write Tecplot file in every "interval"
	writeData2(itest, step);

	// Write Flow history in every "WRITE_FLOWHISTORY_INTERVAL"
	int writeFlowHistroy_interval = getIntParam("WRITE_FLOWHISTORY_INTERVAL");
	if(step > 0) {
		int nScal = scalarTranspEqVector.size();
		double metricNumer[5+nScal];
		double metricDenom[5+nScal];
		for(int i=0; i<5+nScal; ++i) { metricNumer[i] = 0.0; 	metricDenom[i] = 0.0; }

		if(perturbParams.calcOptimalMetric) {
			if(step == 1 || step%writeFlowHistroy_interval == 0)
				calcOptMetrics(metricNumer, metricDenom);
		}

		if(step == 1 || step%writeFlowHistroy_interval == 0) {
			char flowHistoryFilename[60];
			sprintf(flowHistoryFilename, "FlowHistory_test%04d.txt", itest);
			bool rewrite = (step == 1);

			writeFlowHistory(flowHistoryFilename, metricNumer, metricDenom, rewrite);
		}
	}

	// Dump out restart file (and Tecplot file too) before the simulation is forced to quit due to timeout
	double autosaveHours = getDoubleParam("AUTOSAVE_HOURS", "23.8");
	static bool turnOnAutosave = true;
	if(turnOnAutosave && step > 0) {
		double wtime = MPI_Wtime();
		double wtimeHour = (wtime - vanGogh_startWtime) / 3600;
		if(wtimeHour > autosaveHours) {
			if(mpi_rank == 0)
				cout<<"> Autosave current result: AUTOSAVE_HOURS = "<<autosaveHours<<endl;

			char restartFilename[40];
			sprintf(restartFilename, "Restart.test%04d.%06d.out", itest, step);
			writeRestart(restartFilename);

			writeData2(itest, step, true);
		}
	}
	turnOnAutosave = false;
}

/*
 * Method: finalHook
 * -----------------
 *
 */
void VanGoghWithModels::finalHook() {
	char filename[30];
	sprintf(filename, "Q1_pt%04d.bin", itest);
	writeJOEDataParallel(filename);
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

	assert(Evals.empty());

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
		Evals.allocate(nev);  // Note: A data container that has the 'std::complex' class as its element
		double dummyReal, dummyImag;
		for(int iev=0; iev<nev; ++iev) {
			infile.read(reinterpret_cast<char*>(&dummyReal), sizeof(double));
			infile.read(reinterpret_cast<char*>(&dummyImag), sizeof(double));
			Evals[iev] = complex<double>(dummyReal, dummyImag);
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
			for(int iEqn=0; iEqn<NcontrolEqns_file; ++iEqn)
				printf("%.8f, ", lambda_file[iEqn]);
			printf("total number of CVs=%d \n", cvora_file[mpi_size_file]);
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
				printf("    [%2d] %g + %gi\n", int(iev), Evals[iev].real(), Evals[iev].imag());
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

	// Check the range of the direct vectors and the adjoint vectors
	double myMinMagDirectEvec[5+nScal];  	for(int i=0; i<5+nScal; ++i) myMinMagDirectEvec[i] =  ABSURDLY_BIG_NUMBER;
	double myMaxMagDirectEvec[5+nScal];  	for(int i=0; i<5+nScal; ++i) myMaxMagDirectEvec[i] = -ABSURDLY_BIG_NUMBER;
	double myMinMagAdjointEvec[5+nScal]; 	for(int i=0; i<5+nScal; ++i) myMinMagAdjointEvec[i] =  ABSURDLY_BIG_NUMBER;
	double myMaxMagAdjointEvec[5+nScal]; 	for(int i=0; i<5+nScal; ++i) myMaxMagAdjointEvec[i] = -ABSURDLY_BIG_NUMBER;

	double myMinRealDirectEvec[5+nScal];  	for(int i=0; i<5+nScal; ++i) myMinRealDirectEvec[i] =  ABSURDLY_BIG_NUMBER;
	double myMaxRealDirectEvec[5+nScal];  	for(int i=0; i<5+nScal; ++i) myMaxRealDirectEvec[i] = -ABSURDLY_BIG_NUMBER;
	double myMinRealAdjointEvec[5+nScal];  	for(int i=0; i<5+nScal; ++i) myMinRealAdjointEvec[i] =  ABSURDLY_BIG_NUMBER;
	double myMaxRealAdjointEvec[5+nScal];  	for(int i=0; i<5+nScal; ++i) myMaxRealAdjointEvec[i] = -ABSURDLY_BIG_NUMBER;

	double myMinImagDirectEvec[5+nScal];  	for(int i=0; i<5+nScal; ++i) myMinImagDirectEvec[i] =  ABSURDLY_BIG_NUMBER;
	double myMaxImagDirectEvec[5+nScal];  	for(int i=0; i<5+nScal; ++i) myMaxImagDirectEvec[i] = -ABSURDLY_BIG_NUMBER;
	double myMinImagAdjointEvec[5+nScal];  	for(int i=0; i<5+nScal; ++i) myMinImagAdjointEvec[i] =  ABSURDLY_BIG_NUMBER;
	double myMaxImagAdjointEvec[5+nScal];  	for(int i=0; i<5+nScal; ++i) myMaxImagAdjointEvec[i] = -ABSURDLY_BIG_NUMBER;

	for(int icv=0; icv<ncv; ++icv) {
		for(size_t ivar=0; ivar<5+nScal; ++ivar) {
			// Range of magnitude
			double mag = sqrt( pow(DirectEvecs[0][ivar][icv].real(), 2.0) + pow(DirectEvecs[0][ivar][icv].imag(), 2.0) );
			myMinMagDirectEvec[ivar] = min(myMinMagDirectEvec[ivar], mag);
			myMaxMagDirectEvec[ivar] = max(myMaxMagDirectEvec[ivar], mag);

			mag = sqrt( pow(AdjointEvecs[0][ivar][icv].real(), 2.0) + pow(AdjointEvecs[0][ivar][icv].imag(), 2.0) );
			myMinMagAdjointEvec[ivar] = min(myMinMagAdjointEvec[ivar], mag);
			myMaxMagAdjointEvec[ivar] = max(myMaxMagAdjointEvec[ivar], mag);

			// Range of real part
			myMinRealDirectEvec[ivar] = min(myMinRealDirectEvec[ivar], DirectEvecs[0][ivar][icv].real());
			myMaxRealDirectEvec[ivar] = max(myMaxRealDirectEvec[ivar], DirectEvecs[0][ivar][icv].real());

			myMinRealAdjointEvec[ivar] = min(myMinRealAdjointEvec[ivar], AdjointEvecs[0][ivar][icv].real());
			myMaxRealAdjointEvec[ivar] = max(myMaxRealAdjointEvec[ivar], AdjointEvecs[0][ivar][icv].real());

			// Range of imag part
			myMinImagDirectEvec[ivar] = min(myMinImagDirectEvec[ivar], DirectEvecs[0][ivar][icv].imag());
			myMaxImagDirectEvec[ivar] = max(myMaxImagDirectEvec[ivar], DirectEvecs[0][ivar][icv].imag());

			myMinImagAdjointEvec[ivar] = min(myMinImagAdjointEvec[ivar], AdjointEvecs[0][ivar][icv].imag());
			myMaxImagAdjointEvec[ivar] = max(myMaxImagAdjointEvec[ivar], AdjointEvecs[0][ivar][icv].imag());
		}
	}

	double minMagDirectEvec[5+nScal];
	double maxMagDirectEvec[5+nScal];
	double minMagAdjointEvec[5+nScal];
	double maxMagAdjointEvec[5+nScal];
	MPI_Allreduce(myMinMagDirectEvec,  minMagDirectEvec,  5+nScal, MPI_DOUBLE, MPI_MIN, mpi_comm);
	MPI_Allreduce(myMaxMagDirectEvec,  maxMagDirectEvec,  5+nScal, MPI_DOUBLE, MPI_MAX, mpi_comm);
	MPI_Allreduce(myMinMagAdjointEvec, minMagAdjointEvec, 5+nScal, MPI_DOUBLE, MPI_MIN, mpi_comm);
	MPI_Allreduce(myMaxMagAdjointEvec, maxMagAdjointEvec, 5+nScal, MPI_DOUBLE, MPI_MAX, mpi_comm);

	double minRealDirectEvec[5+nScal];
	double maxRealDirectEvec[5+nScal];
	double minRealAdjointEvec[5+nScal];
	double maxRealAdjointEvec[5+nScal];
	MPI_Allreduce(myMinRealDirectEvec,  minRealDirectEvec,  5+nScal, MPI_DOUBLE, MPI_MIN, mpi_comm);
	MPI_Allreduce(myMaxRealDirectEvec,  maxRealDirectEvec,  5+nScal, MPI_DOUBLE, MPI_MAX, mpi_comm);
	MPI_Allreduce(myMinRealAdjointEvec, minRealAdjointEvec, 5+nScal, MPI_DOUBLE, MPI_MIN, mpi_comm);
	MPI_Allreduce(myMaxRealAdjointEvec, maxRealAdjointEvec, 5+nScal, MPI_DOUBLE, MPI_MAX, mpi_comm);

	double minImagDirectEvec[5+nScal];
	double maxImagDirectEvec[5+nScal];
	double minImagAdjointEvec[5+nScal];
	double maxImagAdjointEvec[5+nScal];
	MPI_Allreduce(myMinImagDirectEvec,  minImagDirectEvec,  5+nScal, MPI_DOUBLE, MPI_MIN, mpi_comm);
	MPI_Allreduce(myMaxImagDirectEvec,  maxImagDirectEvec,  5+nScal, MPI_DOUBLE, MPI_MAX, mpi_comm);
	MPI_Allreduce(myMinImagAdjointEvec, minImagAdjointEvec, 5+nScal, MPI_DOUBLE, MPI_MIN, mpi_comm);
	MPI_Allreduce(myMaxImagAdjointEvec, maxImagAdjointEvec, 5+nScal, MPI_DOUBLE, MPI_MAX, mpi_comm);

	if(mpi_rank == 0) {
		cout<<"  RANGE OF 1ST DIRECT-VECTORS :"<<endl;
		cout<<"    MAGNITUDE = ";
		for(int ivar=0; ivar<5+nScal-1; ++ivar)
			printf("%.3e~%.3e | ", minMagDirectEvec[ivar], maxMagDirectEvec[ivar]);
		printf("%.3e~%.3e\n", minMagDirectEvec[5+nScal-1], maxMagDirectEvec[5+nScal-1]);
		cout<<"    REAL = ";
		for(int ivar=0; ivar<5+nScal-1; ++ivar)
			printf("%.3e~%.3e | ", minRealDirectEvec[ivar], maxRealDirectEvec[ivar]);
		printf("%.3e~%.3e\n", minRealDirectEvec[5+nScal-1], maxRealDirectEvec[5+nScal-1]);
		cout<<"    IMAG = ";
		for(int ivar=0; ivar<5+nScal-1; ++ivar)
			printf("%.3e~%.3e | ", minImagDirectEvec[ivar], maxImagDirectEvec[ivar]);
		printf("%.3e~%.3e\n", minImagDirectEvec[5+nScal-1], maxImagDirectEvec[5+nScal-1]);

		cout<<"  RANGE OF 1ST ADJOINT-VECTORS :"<<endl;
		cout<<"    MAGNITUDE = ";
		for(int ivar=0; ivar<5+nScal-1; ++ivar)
			printf("%.3e~%.3e | ", minMagAdjointEvec[ivar], maxMagAdjointEvec[ivar]);
		printf("%.3e~%.3e\n", minMagAdjointEvec[5+nScal-1], maxMagAdjointEvec[5+nScal-1]);
		cout<<"    REAL = ";
		for(int ivar=0; ivar<5+nScal-1; ++ivar)
			printf("%.3e~%.3e | ", minRealAdjointEvec[ivar], maxRealAdjointEvec[ivar]);
		printf("%.3e~%.3e\n", minRealAdjointEvec[5+nScal-1], maxRealAdjointEvec[5+nScal-1]);
		cout<<"    IMAG = ";
		for(int ivar=0; ivar<5+nScal-1; ++ivar)
			printf("%.3e~%.3e | ", minImagAdjointEvec[ivar], maxImagAdjointEvec[ivar]);
		printf("%.3e~%.3e\n", minImagAdjointEvec[5+nScal-1], maxImagAdjointEvec[5+nScal-1]);

		cout<<endl;
	}

	MPI_Barrier(mpi_comm);
}

/*
 * Method: readPsALCdumpedDataSerial
 * ---------------------------------
 * Read field data from file
 * Original code: IkeWithModels::readPsALCdumpedDataSerial()
 */
void VanGoghWithModels::readPsALCdumpedDataSerial(const char filename[], double* qVec, const int NcontrolEqns) {
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

		if(mpi_rank==0) {
			printf("\nReading \'%s\': step=%d, lambda1=", filename, step0);
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
				printf("%.8f, ", lambda1[iEqn]);
			printf("ds=%.8f, total number of CVs=%d, mpi_size of the file=%d\n", dsTry, cvora_file[mpi_size_file], mpi_size_file);
			cout<<"  Grid tolerance = "<<gridTol<<endl;
			if(VANGOGH_DEBUG_LEVEL > 1) {
				cout<<"  XCV_BOX = "<<endl;
				for(int impi=0; impi<mpi_size_file; ++impi)
					printf("    %3d: min=(%.3f, %.3f, %.3f)   max=(%.3f, %.3f, %.3f)\n", impi,
							xcvMinArray_file[impi*3],xcvMinArray_file[impi*3+1],xcvMinArray_file[impi*3+2],
							xcvMaxArray_file[impi*3],xcvMaxArray_file[impi*3+1],xcvMaxArray_file[impi*3+2]);
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
						for(int i=0; i<m; ++i)
							qVec[m*icv+i] = qVec_file[m*icv_file+i];
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
 * Method: writeJOEDataParallel
 * ----------------------------
 * Write JOE data on a binary file.
 * Original function = writeJOEDataParallel() in JoeUgpWithCvCompFlow.h
 */
void VanGoghWithModels::writeJOEDataParallel(char filename[]) {
	string funcID = "VanGoghWithModels::writeJOEDataParallel";

	// Number of variables
	int nScal = scalarTranspEqVector.size();
	int m = nScal+5; // number of variables (rho, rhou, rhov, rhow, rhoE, and scalars)

	if(mpi_rank==0)
		cout<<funcID<<"(): write data on "<<filename<<endl;

	// xcvMin and xcvMax
	double xcvMin[3] = {2.2e22, 2.2e22, 2.2e22}, 	xcvMax[3] = {-2.2e22, -2.2e22, -2.2e22};
	double *xcvMinArray = new double [mpi_size*3];
	double *xcvMaxArray = new double [mpi_size*3];
	for(int icv=0; icv<ncv; ++icv) {
		for(int i=0; i<3; ++i) {
			xcvMin[i] = min(xcvMin[i], x_cv[icv][i]);
			xcvMax[i] = max(xcvMax[i], x_cv[icv][i]);
		}
	}
	MPI_Gather(xcvMin, 3, MPI_DOUBLE, xcvMinArray, 3, MPI_DOUBLE, 0, mpi_comm);
	MPI_Gather(xcvMax, 3, MPI_DOUBLE, xcvMaxArray, 3, MPI_DOUBLE, 0, mpi_comm);
	MPI_Barrier(mpi_comm);

	/***********
	 ** Write the data
	 ***********/
	double wtime0, wtimeF;
	if(mpi_rank==0)
		wtime0 = MPI_Wtime();

	// 1. Header
	//      Structure of the header part:
	//              1. mpi_size (int)                  2. cvora (int*(mpi_size+1))
	//              3. xMinArray (double*3*mpi_size)   4. xMaxArray (double*3*mpi_size)
	//              5. nScal (int)
	if(mpi_rank==0) {
		// Open the file
		ofstream ofile;
		ofile.open(filename, ios_base::out | ios_base::trunc | ios_base::binary);

		// Write on the file
		int dummyInt;
		double dummyDouble;

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

		dummyInt=nScal; 			ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));

		// Close the file
		ofile.close();
	}
	MPI_Barrier(mpi_comm);
	int initDisp = sizeof(int) * (3 + mpi_size) + sizeof(double) * (2*mpi_size*3);

	// 2. Body
	//       Structure of the body part:
	//            For each mpi,
	//              1. x_cv (double*3*ncv)
	//              2. data (double*(5+nScal)*ncv)
	MPI_Status status;
	MPI_File fh;
	MPI_Offset displacement;

	// Open the file
	if (MPI_File_open(mpi_comm, filename, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh)!=0) {
			// Note: MPI_MODE_WRONLY = write only
			//       MPI_MODE_CREATE = create the file if it does not exist.
		cerr<<"ERROR "<<funcID<<"(): Cannot open "<<filename<<endl;
		throw(-1);
	}

	// Write the CV coordinates first
	int myOffsetNcv;
	MPI_Scan(&ncv, &myOffsetNcv, 1, MPI_INT, MPI_SUM, mpi_comm);
	myOffsetNcv -= ncv;
	displacement = MPI_Offset(initDisp + sizeof(double)*myOffsetNcv*(3+m));
	if(MPI_File_set_view(fh, displacement, MPI_INT, MPI_INT, "native", MPI_INFO_NULL)!=0) {
		// Note: native     = Data in this representation are stored in a file exactly as it is in memory
		//       Internal   = Data in this representation can be used for I/O operations in a homogeneous or heterogeneous environment
		//       External32 = This data representation states that read and write operations convert all data from and to the "external32" representation
		cerr<<"ERROR! "<<funcID<<"(): Cannot set the first MPI_File_set_view -- offset="<<displacement<<endl;
		throw(-1);
	}
	double* bufferDouble = new double [ncv*3];
	for(int icv=0; icv<ncv; ++icv) {
		int indexStart = icv*3;
		for(int i=0; i<3; ++i)
			bufferDouble[indexStart+i] = x_cv[icv][i];
	}
	MPI_File_write(fh, bufferDouble, ncv*3, MPI_DOUBLE, &status);
	delete [] bufferDouble; 	bufferDouble = NULL;

	// Write the flow data
	displacement += MPI_Offset(3*ncv*sizeof(double));
	bufferDouble = new double [ncv*m];
	for(int icv=0; icv<ncv; ++icv) {
		int indexStart = icv*m;
		bufferDouble[indexStart] = rho[icv];
		for(int i=0; i<3; ++i)
			bufferDouble[indexStart+1+i] = rhou[icv][i];
		bufferDouble[indexStart+4] = rhoE[icv];
		for(int iScal=0; iScal<nScal; ++iScal) {
			double *phi = scalarTranspEqVector[iScal].phi;
			bufferDouble[indexStart+5+iScal] = phi[icv];
		}
	}
	MPI_File_write(fh, bufferDouble, ncv*m, MPI_DOUBLE, &status);
	delete [] bufferDouble; 	bufferDouble = NULL;

	// 2-4. Close the file
	MPI_File_close(&fh);
	MPI_Barrier(mpi_comm);

	/*
	 * 3. Write the foot
	 */
	if(mpi_rank == mpi_size-1) {
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
	if(mpi_rank==0)
		cout<<"> The flow data is written on "<<filename<<" (RUN TIME = "<<wtimeF-wtime0<<" [sec]) \n"<<endl;

	/***********
	 ** Free the memory
	 ***********/
	delete [] xcvMinArray;
	delete [] xcvMaxArray;
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
 * ------------------------
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
		perturbParams.disturbSmoothXmin = getParam("SMOOTHING")->getDouble("SMOOTH_XMIN");
		perturbParams.disturbSmoothXmax = getParam("SMOOTHING")->getDouble("SMOOTH_XMAX");
		perturbParams.disturbSmoothEdgeSize = getParam("SMOOTHING")->getDouble("SMOOTH_EDGE_SIZE");
	}
	
	string boolStringOM = getStringParam("CALC_OPTIMAL_METRIC", "TRUE");
	std::transform(boolStringOM.begin(), boolStringOM.end(), boolStringOM.begin(), ::tolower);
	if(boolStringOM.compare("yes")==0 || boolStringOM.compare("y")==0 || boolStringOM.compare("true")==0 || boolStringOM.compare("t")==0) {
		perturbParams.calcOptimalMetric = true;
	}
}

/*
 * Method: reinitialHook
 * ---------------------
 * Similar role to initialHook(): Reinitialize the flow field from the previous simulation.
 */
void VanGoghWithModels::reinitialHook() {
	assert(rho_init!=NULL && rhou_init!=NULL && rhoE_init!=NULL);

	for (int icv=0; icv<ncv; ++icv) {
		rho[icv] = rho_init[icv];
		for(int i=0; i<3; ++i)
			rhou[icv][i] = rhou_init[icv][i];
		rhoE[icv] = rhoE_init[icv];
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

	double myMinPress = ABSURDLY_BIG_NUMBER;
	double myMinTemp  = ABSURDLY_BIG_NUMBER;
	for(int icv=0; icv<ncv; ++icv) {
		myMinPress = max(min(myMinPress, UgpWithCvCompFlow::press[icv]), 1.0e-6);
		myMinPress = max(min(myMinTemp,  UgpWithCvCompFlow::temp[icv]),  1.0e-6);
	}
	double minPress;
	double minTemp;
	MPI_Allreduce(&myMinPress, &minPress, 1, MPI_DOUBLE, MPI_MIN, mpi_comm);
	MPI_Allreduce(&myMinTemp,  &minTemp,  1, MPI_DOUBLE, MPI_MIN, mpi_comm);

	// ---------------
	// Perturb rho
	// ---------------
	double disturbRatio = perturbParams.disturbMag;
	bool applyClipping = true;
	perturbScalar(rho, array_perturb, disturbRatio, "rho", applyClipping);  // Note: updateCvDataG1G2() was called in this method

	// write the perturbation for Tecplot output
	for(int icv=0; icv<ncv; ++icv)
		rho_perturb[icv] = array_perturb[icv];
	updateCvDataG1G2(rho_perturb, REPLACE_DATA);

	// ---------------
	// Perturb rhou-X
	// ---------------
	disturbRatio = perturbParams.disturbMag;
	int coord = 0;
	applyClipping = false;
	perturbVector(rhou, coord, array_perturb, disturbRatio, "rhou-X", applyClipping);  // Note: updateCvDataG1G2() was called in this method

	// write the perturbation for Tecplot output
	for(int icv=0; icv<ncv; ++icv)
		rhou_perturb[icv][coord] = array_perturb[icv];

	// ---------------
	// Perturb rhou-Y
	// ---------------
	coord = 1;
	disturbRatio = perturbParams.disturbMag;
	applyClipping = false;
	perturbVector(rhou, coord, array_perturb, disturbRatio, "rhou-Y", applyClipping);  // Note: updateCvDataG1G2() was called in this method

	// write the perturbation for Tecplot output
	for(int icv=0; icv<ncv; ++icv)
		rhou_perturb[icv][coord] = array_perturb[icv];
	updateCvDataG1G2(rhou_perturb, REPLACE_ROTATE_DATA); // Update all the three coordinates!!

	// ---------------
	// Perturb rhoE
	// ---------------
	disturbRatio = perturbParams.disturbMag;
	applyClipping = true;
	perturbScalar(rhoE, array_perturb, disturbRatio, "rhoE", applyClipping);  // Note: updateCvDataG1G2() was called in this method

	// write the perturbation for Tecplot output
	for(int icv=0; icv<ncv; ++icv)
		rhoE_perturb[icv] = array_perturb[icv];
	updateCvDataG1G2(rhoE_perturb, REPLACE_DATA);

	// Check if pressure or temperature becomes negative even though rhoE has been already clipped.
	double *kineArray = NULL;
	int kine_Index = getScalarTransportIndex("kine");
	if(kine_Index>=0)
		kineArray = scalarTranspEqVector[kine_Index].phi;

	double MIN_PRESS = max(0.01 * minPress, 1.0e-6);
	double MIN_TEMP  = max(0.01 * minTemp,  1.0e-6);

	int myCountClip = 0;
	double myMaxAddEnergy = 0.0;
	for(int icv=0; icv<ncv; ++icv) {
		double addEnergy = 0.0;

		// Check pressure
		double kinecv = 0.0;
		if(kine_Index>=0)
			kinecv = kineArray[icv];
		double pressCheck = calcPress(UgpWithCvCompFlow::gamma[icv], rhoE[icv], rhou[icv], rho[icv], kinecv);
		if(pressCheck < 0.0)
			addEnergy = fabs(pressCheck + MIN_PRESS) / (UgpWithCvCompFlow::gamma[icv] - 1.0);

		// Check energy
		double tempCheck = calcTemp(pressCheck, rho[icv], UgpWithCvCompFlow::RoM[icv]);
		if(tempCheck < 0.0)
			addEnergy = max(addEnergy,  rho[icv] * fabs(tempCheck + MIN_TEMP) * UgpWithCvCompFlow::RoM[icv] / (UgpWithCvCompFlow::gamma[icv] - 1.0));

		// Add energy
		if(addEnergy != 0.0) {
			rhoE[icv] += addEnergy;

			++myCountClip;
			myMaxAddEnergy = max(myMaxAddEnergy, addEnergy);
		}
	}
	updateCvDataG1G2(rhoE, REPLACE_DATA);

	int countClip = 0;
	double maxAddEnergy = 0.0;
	MPI_Allreduce(&myCountClip,    &countClip,     1, MPI_INT,    MPI_SUM, mpi_comm);
	MPI_Allreduce(&myMaxAddEnergy, &maxAddEnergy,  1, MPI_DOUBLE, MPI_MAX, mpi_comm);
	if(countClip > 0 && mpi_rank==0)
		cout<<"    rhoE is clipped again due to negative pressure or temperature: count CV = "<<countClip<<", max added energy = "<<maxAddEnergy<<endl;


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
 *           array_perturb[icv] = coeff * RANDOM_VARIABLE
 */
double VanGoghWithModels::perturbScalar(double* scalarArray, double* array_perturb, const double disturbRatio, const char varName[], const bool applyClipping) {
	double rmsVal = calcRMS(scalarArray, false); // get RMS
	double coeff = rmsVal * disturbRatio; // calculate the coefficient based on the RMS

	perturbScalar(scalarArray, array_perturb, coeff, varName, applyClipping);

	return coeff;
}

/*
 * Method: perturbScalarArray
 * --------------------------
 * Perturb the initial field with random variables and apply filter to the perturbations
 * Note: The array "array_perturb" is used only in this method, but sometimes the user may want to take a look at it.
 *       Thus, it is given as an argument.
 *
 * Arguments: coeff = interval size for UNIFORM, std for GAUSSIAN
 *
 * Return: array_perturb[icv] = coeff * RANDOM_VARIABLE
 */
void VanGoghWithModels::perturbScalarArray(double* scalarArray, double* array_perturb, const double coeff, const char varName[], const bool applyClipping) {
	assert(scalarArray != NULL && array_perturb != NULL);

	// Smoothing: Some parameters for the hat-shaped function
	//   f(x) = 1.0 - 0.5*(exp(-a*(x-xMin)) + exp(a*(x-xMax))), where a = ln(100)/edgeSize
	// Note: At xMin & xMax, f(x) = 0.0. At the "edges", f(x) = 0.99. At the center of the domain, f(x) ~= 1.0
	//       The smoothing occurs only if useSmoothing is TRUE.
	double a = 0.0;
	if(perturbParams.useSmoothing) {
		a = log(100.0)/perturbParams.disturbSmoothEdgeSize;
		assert(perturbParams.disturbSmoothXmin < perturbParams.disturbSmoothXmax);
		assert(a>0);
	}

	// -----------------------
	// Generate a random array
	// -----------------------
	if(perturbParams.randDistribFunc == UNIFORM_DISTRIB) {
		randArrayUniform(array_perturb, ncv); // note: randArrayUniform() generates random variables in [0,1]
		for(int icv=0; icv<ncv; ++icv)
			array_perturb[icv] -= 0.5;        // shifting
	} else
		randArrayNormal(array_perturb, ncv);  // Note: default = normal distribution
	updateCvDataG1G2(array_perturb, REPLACE_DATA);

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

		applyDiffFilterG1G2(tempFiltered, array_perturb, perturbParams.filterAbsTol, perturbParams.filterMaxIter);
		for(int icv=0; icv<ncv; ++icv)
			array_perturb[icv] = tempFiltered[icv];
//		updateCvDataG1G2(array_perturb, REPLACE_DATA);

		delete [] tempFiltered;
	}

	// ------------------------------
	// Update the designated variable (with clipping if necessary)
	// ------------------------------
	int myCountClip = 0;
	for(int icv=0; icv<ncv; ++icv) {
		array_perturb[icv] *= coeff;
		if(perturbParams.useSmoothing) {
			double x = x_cv[icv][0];
			if(x<perturbParams.disturbSmoothXmin || x>perturbParams.disturbSmoothXmax)
				array_perturb[icv] = 0.0;
			else
				array_perturb[icv] *= 1.0 - 0.5*(exp(-a*(x-perturbParams.disturbSmoothXmin)) + exp(a*(x-perturbParams.disturbSmoothXmax)));
		}
		if(applyClipping && scalarArray[icv]+array_perturb[icv] < perturbParams.disturbClip*scalarArray[icv]) {
			scalarArray[icv] = perturbParams.disturbClip*scalarArray[icv]; // clipping
			myCountClip++;
		} else {
			scalarArray[icv] += array_perturb[icv];
		}
	}
	updateCvDataG1G2(scalarArray,  REPLACE_DATA);

	int countClip; 	MPI_Allreduce(&myCountClip, &countClip, 1, MPI_INT, MPI_SUM, mpi_comm);
	if(mpi_rank == 0) {
		if(applyClipping) {
			if(perturbParams.randDistribFunc == UNIFORM_DISTRIB)
				printf("    Perturb the %s field with UNIFORM:  min disturb = %.3e, max disturb = %.3e, # clip = %d \n", varName, -0.5*coeff, 0.5*coeff, countClip);
			else
				printf("    Perturb the %s field with GAUSSIAN:  std = %.3e, # clip = %d \n", varName, coeff, countClip);
		} else {
			if(perturbParams.randDistribFunc == UNIFORM_DISTRIB)
				printf("    Perturb the %s field with UNIFORM:  min disturb = %.3e, max disturb = %.3e\n", varName, -0.5*coeff, 0.5*coeff);
			else
				printf("    Perturb the %s field with GAUSSIAN:  std = %.3e\n", varName, coeff);
		}
	}

	// Error check: negative value
	if(applyClipping) {
		int myCountNegative = 0;
		double myMostNegative = ABSURDLY_BIG_NUMBER;

		for(int icv=0; icv<ncv; ++icv) {
			if(scalarArray[icv] < 0.0) {
				++myCountNegative;
				myMostNegative = min(myMostNegative, scalarArray[icv]);
			}
		}

		int countNegative;
		double mostNegative;
		MPI_Allreduce(&myCountNegative, &countNegative, 1, MPI_INT,    MPI_SUM, mpi_comm);
		MPI_Allreduce(&myMostNegative,  &mostNegative,  1, MPI_DOUBLE, MPI_MIN, mpi_comm);

		if(countNegative > 0 && mpi_rank==0)
			printf("    WARNING! %s still has negative value(s) even after clipping: count=%d, most negative value=%g\n", varName, countNegative, mostNegative);
	}
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
 *           array_perturb[icv] = coeff * RANDOM_VARIABLE
 */
double VanGoghWithModels::perturbVector(double (*vectorArray)[3], const int coord, double* array_perturb, const double disturbRatio, const char varName[], const bool applyClipping) {
	double rmsVal = calcRMS3D(vectorArray, coord, false); // get RMS
	double coeff = rmsVal * disturbRatio; // calculate the coefficient based on the RMS

	perturbVectorArray(vectorArray, coord, array_perturb, coeff, varName, applyClipping);

	return coeff;
}

/*
 * Method: perturbVectorArray
 * --------------------------
 * Perturb the initial field with random variables and apply filter to the perturbations
 * Note: The array "array_perturb" is used only in this method, but sometimes the user may want to take a look at it.
 *       Thus, it is given as an argument.
 *
 * Arguments: coeff = interval size for UNIFORM, std for GAUSSIAN
 *
 * Return: array_perturb[icv] = coeff * RANDOM_VARIABLE
 */
void VanGoghWithModels::perturbVectorArray(double (*vectorArray)[3], const int coord, double* array_perturb, const double coeff, const char varName[], const bool applyClipping) {
	assert(vectorArray != NULL && array_perturb != NULL);
	assert(coord>=0 && coord<3);

	// Smoothing: Some parameters for the hat-shaped function
	//   f(x) = 1.0 - 0.5*(exp(-a*(x-xMin)) + exp(a*(x-xMax))), where a = ln(100)/edgeSize
	// Note: At xMin & xMax, f(x) = 0.0. At the "edges", f(x) = 0.99. At the center of the domain, f(x) ~= 1.0
	//       The smoothing occurs only if useSmoothing is TRUE.
	double a = 0.0;
	if(perturbParams.useSmoothing) {
		a = log(100.0)/perturbParams.disturbSmoothEdgeSize;
		assert(perturbParams.disturbSmoothXmin < perturbParams.disturbSmoothXmax);
		assert(a>0);
	}

	// -----------------------
	// Generate a random array
	// -----------------------
	if(perturbParams.randDistribFunc == UNIFORM_DISTRIB) {
		randArrayUniform(array_perturb, ncv); // note: randArrayUniform() generates random variables in [0,1]
		for(int icv=0; icv<ncv; ++icv)
			array_perturb[icv] -= 0.5;        // shifting
	} else
		randArrayNormal(array_perturb, ncv);  // Note: default = normal distribution
	updateCvDataG1G2(array_perturb, REPLACE_DATA);

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

		applyDiffFilterG1G2(tempFiltered, array_perturb, perturbParams.filterAbsTol, perturbParams.filterMaxIter);
		for(int icv=0; icv<ncv; ++icv)
			array_perturb[icv] = tempFiltered[icv];
//		updateCvDataG1G2(array_perturb, REPLACE_DATA);

		delete [] tempFiltered;
	}

	// ------------------------------
	// Update the designated variable (with clipping if necessary)
	// ------------------------------
	int myCountClip = 0;
	for(int icv=0; icv<ncv; ++icv) {
		array_perturb[icv] *= coeff;
		if(perturbParams.useSmoothing) {
			double x = x_cv[icv][0];
			if(x<perturbParams.disturbSmoothXmin || x>perturbParams.disturbSmoothXmax)
				array_perturb[icv] = 0.0;
			else
				array_perturb[icv] *= 1.0 - 0.5*(exp(-a*(x-perturbParams.disturbSmoothXmin)) + exp(a*(x-perturbParams.disturbSmoothXmax)));
		}

		if(applyClipping && vectorArray[icv][coord]+array_perturb[icv] < perturbParams.disturbClip*vectorArray[icv][coord]) {
			vectorArray[icv][coord] = perturbParams.disturbClip * vectorArray[icv][coord]; // clipping
			myCountClip++;
		} else {
			vectorArray[icv][coord] += array_perturb[icv];
		}
	}
	updateCvDataG1G2(vectorArray, REPLACE_ROTATE_DATA); // Update all the three coordinates!!

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
}

/*
 * Method: specify_filter_width
 * ----------------------------
 * Original method = DiffFilter::specify_filter_width()
 *     delta^2
 * p = -------, where delta = filter width, if p is constant everywhere
 *      40.0
 */
double VanGoghWithModels::specify_filter_width(const int ifa) {
	assert(ifa < nfa_b2gg);  // The size of the wdFace array is equal to nfa_b2gg.
	assert(wdFace[ifa] > 0.0);
	
	if(wdFace[ifa] >= perturbParams.filterBoundaryLength)
		return pow(perturbParams.maxFilterWidth, 2.0) / 40.0;
	
	double ratio = wdFace[ifa] / perturbParams.filterBoundaryLength;
	double filterWidth = perturbParams.minFilterWidth + ratio*(perturbParams.maxFilterWidth - perturbParams.minFilterWidth);
	assert(filterWidth > 0.0 && filterWidth <= perturbParams.maxFilterWidth);
	
	return pow(filterWidth, 2.0) / 40.0;
}

/****************************
 * OPTIMAL METRIC
 ****************************/

/*
 * Method: calcOptMetrics
 * ----------------------
 * Calculate the optimal metrics using adjoint. 
 * Formulation: Let phi=flow field, phi0=flow field on the unstable branch, y=least stable eigenvector, a=adjoint vector,
 *                       (phi-phi0)*a
 *              metric = ------------
 *                           y*a
 * Return: metricNumer = (phi-phi0)*a
 *         metricDenom = y*a
 */
void VanGoghWithModels::calcOptMetrics(double* metricNumer, double* metricDenom) {
	assert(metricNumer != NULL && metricDenom != NULL);
	
	// define variables
	int nScal = scalarTranspEqVector.size();
	int m = 5+nScal; // number of variables

	// initialize the solution vectors
	for(int i=0; i<m; ++i) {
		metricNumer[i] = 0.0;
		metricDenom[i] = 0.0;
	}

	// calculate numerator(=(phi-phi0)*a) and denominator(=y*a)
	calcOptMetricsNS(metricNumer, metricDenom, nScal);
	calcOptMetricsScalarRansTurbModel(metricNumer, metricDenom, nScal);
	calcOptMetricsScalarRansCombModel(metricNumer, metricDenom, nScal);
}

/*
 * Method: calcOptMetricsNS
 * ------------------------
 * Calculate the NS part of the optimal metrics using adjoint. 
 * Formulation: Let phi=flow field, phi0=flow field on the unstable branch, y=least stable eigenvector, a=adjoint vector,
 *                       (phi-phi0)*a
 *              metric = ------------
 *                           y*a
 * Return: metricNumer = (phi-phi0)*a
 *         metricDenom = y*a
 */
void VanGoghWithModels::calcOptMetricsNS(double* metricNumer, double* metricDenom, const int nScal) {
	assert(rho_unstable != NULL && rhou_unstable != NULL && rhoE_unstable != NULL);
	assert(!AdjointEvecs.empty());
	assert(!DirectEvecs.empty());
	
	// define variables
	double myTemp[5];

	// calculate numerator = (phi-phi0)*a
	for(int i=0; i<5; ++i)
		myTemp[i] = 0.0;
	for(int icv=0; icv<ncv; ++icv) {
		myTemp[0] += (rho[icv]-rho_unstable[icv]) * AdjointEvecs[0][0][icv].real();
		myTemp[1] += (rhou[icv][0]-rhou_unstable[icv][0]) * AdjointEvecs[0][1][icv].real();
		myTemp[2] += (rhou[icv][1]-rhou_unstable[icv][1]) * AdjointEvecs[0][2][icv].real();
		myTemp[3] += (rhou[icv][2]-rhou_unstable[icv][2]) * AdjointEvecs[0][3][icv].real();
		myTemp[4] += (rhoE[icv]-rhoE_unstable[icv]) * AdjointEvecs[0][4][icv].real();
	}
	double totTemp[5]; 	for(int i=0; i<5; ++i) totTemp[i] = 0.0;
	MPI_Allreduce(myTemp, totTemp, 5, MPI_DOUBLE, MPI_SUM, mpi_comm);
	for(int i=0; i<5; ++i)
		metricNumer[i] = totTemp[i];

	// calculate denominator = y*a
	for(int i=0; i<5; ++i)
		myTemp[i] = 0.0;
	for(int icv=0; icv<ncv; ++icv) {
		myTemp[0] += DirectEvecs[0][0][icv].real() * AdjointEvecs[0][0][icv].real();
		myTemp[1] += DirectEvecs[0][1][icv].real() * AdjointEvecs[0][1][icv].real();
		myTemp[2] += DirectEvecs[0][2][icv].real() * AdjointEvecs[0][2][icv].real();
		myTemp[3] += DirectEvecs[0][3][icv].real() * AdjointEvecs[0][3][icv].real();
		myTemp[4] += DirectEvecs[0][4][icv].real() * AdjointEvecs[0][4][icv].real();
	}
	for(int i=0; i<5; ++i) totTemp[i] = 0.0;
	MPI_Allreduce(myTemp, totTemp, 5, MPI_DOUBLE, MPI_SUM, mpi_comm);
	for(int i=0; i<5; ++i)
		metricDenom[i] = totTemp[i];
}

/*
 * Method: writeOptimalMeticOnFile
 * -------------------------------
 * Write the optimal metrics (for rho,rhou,etc.) on file
 */
void VanGoghWithModels::writeOptimalMeticOnFile(const int itest, char filename[], bool rewrite, const double* metricNumer, const double* metricDenom, const int nScal) {
	assert(metricNumer!=NULL && metricDenom!=NULL);

	if(mpi_rank==0) {
		FILE *fp;
		if(rewrite) {
			fp = fopen(filename, "w");
			fprintf(fp, "ITEST,    METRIC_RHO, METRIC_RHOU-X, METRIC_RHOU-Y, METRIC_RHOU-Z,   METRIC_RHOE");
			for(int i=0; i<nScal; ++i)
				fprintf(fp, ", METRIC_SCAL-%d", i);
			fprintf(fp, ",   METRIC_FLOW,  TOT\n");
		} else
			fp = fopen(filename, "a");

		fprintf(fp, "%5d", itest);
		for(int i=0; i<5+nScal; ++i) {
			if(i==3) // Note: rhou-Z can be zero
				fprintf(fp, ", %13.6e", metricNumer[i]/(metricDenom[i]+MACHINE_EPS));
			else
				fprintf(fp, ", %13.6e", metricNumer[i]/metricDenom[i]);
		}
		
		double metricNumerSum = 0.0;
		double metricDenomSum = 0.0;

		for(int i=0; i<5; ++i) {
			metricNumerSum += metricNumer[i];
			metricDenomSum += metricDenom[i];
		}
		fprintf(fp, ", %13.6e", metricNumerSum/metricDenomSum); // Write flow metric

		for(int iScal=0; iScal<nScal; ++iScal) {
			metricNumerSum += metricNumer[5+iScal];
			metricDenomSum += metricDenom[5+iScal];
		}
		fprintf(fp, ", %13.6e", metricNumerSum/metricDenomSum); // Write tot metric
		
		fprintf(fp, "\n");

		fclose(fp);
	}

	MPI_Barrier(mpi_comm);
}

/****************************
 * VANGOGH OVERLOADED METHODS
 ****************************/
/*
 * Method: writeResidualOnFile
 * ---------------------------
 *
 */
void VanGoghWithModels::writeResidualOnFile(const int itest, char filename[], bool rewrite) {
	int nScal = scalarTranspEqVector.size();

	if(mpi_rank==0) {
		FILE* fp;
		if(rewrite) {
			fp = fopen(filename, "w");
			if (fp != NULL) {
				fprintf(fp, "  STEP,  rho,          rhou-X,       rhou-Y,       rhou-Z,       rhoE");
				for (int iScal = 0; iScal < nScal; iScal++)
					fprintf(fp, ", %12s", scalarTranspEqVector[iScal].getName());
				fprintf(fp, "\n");
			} else {
				cout<<"ERROR VanGoghWithModels::writeResidualOnFile(): cannot open "<<filename<<endl;
				assert(false);
			}
		} else {
			if (fp != NULL) {
				fp = fopen(filename, "a");
			} else {
				cout<<"ERROR VanGoghWithModels::writeResidualOnFile(): cannot open "<<filename<<endl;
				assert(false);
			}
		}

		if(JoeWithModels::Residual == NULL) {
			cout<<"WARNING "<<classID<<"::writeResidualOnFile(): Residual is NULL"<<endl;
		} else {
			fprintf(fp, "%6d, %12.4e, %12.4e, %12.4e, %12.4e, %12.4e", 
					step, JoeWithModels::Residual[0], JoeWithModels::Residual[1], JoeWithModels::Residual[2], JoeWithModels::Residual[3], JoeWithModels::Residual[4]);
			for (int iScal = 0; iScal < nScal; iScal++) {
				fprintf(fp, ", %12.4e", JoeWithModels::Residual[5+iScal]);
			}
			fprintf(fp, "\n");
		}

		fclose(fp);
	}

	MPI_Barrier(mpi_comm);
}

/*
 * Method: writeFlowHistory
 * --------------------------
 * Write flow history (the optimal metric, other QoI's) on a file
 * This is called by temporalHook() -- i.e. This can be called at each step in a test case.
 * If the end-users want to add more QoI's, they should write their own writeFlowHistory().
 */
void VanGoghWithModels::writeFlowHistory(const char filename[], const double* metricNumer, const double* metricDenom, const bool rewrite) {
	if(perturbParams.calcOptimalMetric) {
		assert(metricNumer != NULL && metricDenom != NULL);
	}

	int nScal = scalarTranspEqVector.size();

	// Calculate the optimal metric from the given numeriators and denomiators
	double OptMetricFlow = 0.0, OptMetricTot = 0.0;
	if(perturbParams.calcOptimalMetric) {
		double metricNumerSum = 0.0;
		double metricDenomSum = 0.0;

		for(int i=0; i<5; ++i) {
			metricNumerSum += metricNumer[i];
			metricDenomSum += metricDenom[i];
		}
		OptMetricFlow = metricNumerSum / metricDenomSum;

		for(int iScal=0; iScal<nScal; ++iScal) {
			metricNumerSum += metricNumer[5+iScal];
			metricDenomSum += metricDenom[5+iScal];
		}
		OptMetricTot = metricNumerSum / metricDenomSum;
	}

	if(mpi_rank==0) {
		FILE *fp;
		if(rewrite) {
			fp = fopen(filename, "w");
			fprintf(fp, "STEP,    METRIC_FLOW,     METRIC_TOT\n");
		} else
			fp = fopen(filename, "a");

		fprintf(fp, "%5d", step);
		fprintf(fp, ", %13.6e, %13.6e",
				OptMetricFlow, OptMetricTot);

		fprintf(fp, "\n");

		fclose(fp);
	}
}
