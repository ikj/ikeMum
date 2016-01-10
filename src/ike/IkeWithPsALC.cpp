/*
 * IkeWithPsALC.cpp
 *
 *  Created on: Nov 10, 2012
 *      Author: ikj
 */

#include "IkeWithPsALC.h"

/*
 * Note: CATCH and RE-THROW
 * ------------------------
 * main() {
 *    ...
 *    try {
 *       runPsALC() {
 *          try {
 *             initFirstTwoPts() {
 *                getSteadySolnByNewton() {
 *                   try {
 *                      calcJacobian1DAD() {
 *                         try {
 *                            ...
 *                         }
 *                         catch {
 *                            throw
 *                         }
 *                      }
 *                   }
 *                   catch {
 *                      throw
 *                   }
 *                }
 *             }
 *          }
 *          catch
 *
 *          try {
 *             initFirstTwoPts() {
 *                getSteadySolnByNewton() {
 *                   try {
 *                      calcJacobian1DAD() {
 *                         try {
 *                            ...
 *                         }
 *                         catch {
 *                            throw
 *                         }
 *                      }
 *                   }
 *                   catch {
 *                      throw
 *                   }
 *                }
 *             }
 *          }
 *          catch
 *       }
 *    }
 *    catch
 * }
 */

/*
 * Method: init
 * ------------
 *
 */
void IkeWithPsALC_AD::init() {
	nScal = -1000;

	// Minimum and Maximum of x_cv
	for(int i=0; i<3; ++i) {
		xcvMin[i] = ABSURDLY_BIG_NUMBER;
		xcvMax[i] = -ABSURDLY_BIG_NUMBER;
	}

#ifdef USE_VOLUME_WEIGHTED_INNER_PROD
	// Minimum and Maximum of cv_volume
	totMaxCvVol = -ABSURDLY_BIG_NUMBER;
	totMinCvVol = ABSURDLY_BIG_NUMBER;
#endif

	// PetscSolvers
	petscSolver2 = NULL;
#ifdef USE_SLEPC_WITH_PETSC
	slepcSolver2 = NULL;
#endif

	// Pseudo-arclength continuation method
	NcontrolEqns = 0;
	lambda = NULL;
	lambda_AD= NULL;

	// RHS 1D Array
	rhs1DArray = NULL;

	// RHS weight
	weightRhs = NULL;
	weightRhsMethod = NO_RHSWEIGHT;

	// RHS arrays (for JOE RHS calculations and Tecplot output)
	RHSrho   = NULL;
	RHSrhou  = NULL;
	RHSrhoE  = NULL;

	registerScalar(RHSrho,  "RHSRHO",  CV_DATA);
	registerVector(RHSrhou, "RHSRHOU", CV_DATA);
	registerScalar(RHSrhoE, "RHSRHOE", CV_DATA);

	RHSrhoScal = NULL;

	// Scalar RHS arrays (for Tecplot output)
	RHSkine  = NULL;
	RHSomega = NULL;
	RHSsa    = NULL;
    RHSZMean = NULL;
    RHSZVar  = NULL;
    RHSCMean = NULL;

	// Sum of the source terms due to Barrier functions
	myBarrierMassSourceSumJOE   = 0.0;
	myBarrierEnergySourceSumJOE = 0.0;
	myBarrierMassSourceSum1D_AD   = 0.0;
	myBarrierEnergySourceSum1D_AD = 0.0;

	// QoI
	rewriteQoI = true;

	// Debugging level
	debugLevel = 0; // default is 0: it will be updated in IkeWithPsALC_AD::run()
}

/*
 * Method: clear
 * -------------
 *
 */
void IkeWithPsALC_AD::clear() {
	if(debugLevel>0 && mpi_rank==0)
		cout<<"IkeWithPsALC_AD::clear(): Free allocated memory..."<<endl;

	// Jacobian matrix
	if (!jacMatrix.empty()) {
		jacMatrix.clear();
	}
	if (!jacMatrixSTL.empty()) {
		jacMatrixSTL.clear();
	}

	// PetscSolvers
	if (petscSolver2 != NULL) {
		delete petscSolver2;
		petscSolver2 = NULL;
	}

	// Pseudo-arclength continuation method
	if(NcontrolEqns > 0) {
		if(lambda != NULL) {
			delete [] lambda; 		lambda = NULL;
		}
		if(lambda_AD != NULL) {
			delete [] lambda_AD; 	lambda_AD = NULL;
		}
	}

	// 1D Arrays
	if(rhs1DArray != NULL) {
		delete [] rhs1DArray;     rhs1DArray = NULL;
	}

	// RHS weight
	if(weightRhs != NULL) {
		delete [] weightRhs;  	weightRhs = NULL;
	}

	// Variables from UgpWithCvCompFlow: These variables are not deleted in the UgpWithCvCompFlow class (This bug will be reported)
	if(nbocv_v_global != NULL) {
		delete [] nbocv_v_global; nbocv_v_global = NULL;
	}
}

/*
 * Method: run
 * -----------
 *
 */
void IkeWithPsALC_AD::run() {
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

	// Get xcvMin and xcvMax
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
			cout<<"CAUTION!! cvora is NULL -- Make sure that you will initialize it"<<endl;
	}

	if (nbocv_v_global == NULL) {
		nbocv_v_global = new int[ncv_gg];
		for (int icv = 0; icv < ncv; icv++)
			nbocv_v_global[icv] = cvora[mpi_rank] + icv;
		updateCvDataG1G2(nbocv_v_global, REPLACE_DATA);

		if(mpi_rank==0)
			cout<<"nbocv_v_global is allocated and initialized"<<endl<<endl;
	}

	// Build the 2-layerd CSR structure (nbocv2_i & nbocv2_v)
	if(mpi_rank==0){
		cout<<"======================="<<endl;
		cout<<"BUILD THE 2-CSR STRUCT."<<endl;
		cout<<"======================="<<endl;
	}
	bool addFakeCells = false;
	build2layerCSRstruct(addFakeCells);
	MPI_Barrier(mpi_comm);
	if(mpi_rank==0)
//		cout<<"nbocv2_v_global, nbocv2_i, and nbocv2_v are successfully constructed"<<endl<<endl;
		cout<<"nbocv2_i, and nbocv2_v are successfully constructed"<<endl<<endl;

	nScal = scalarTranspEqVector.size();

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

	// +++++++++++++++++++++++++++++++++++++
	// Launch the pseudo-arclength algorithm
	// +++++++++++++++++++++++++++++++++++++
	runPsALC();
}

/*
 * Method: runPsALC
 * ----------------
 * Argument: howToCalcJac = How to calculate the Jacobian matrix: ROW_1D, ORDINARY_2D -- default = ORDINARY_2D
 */
void IkeWithPsALC_AD::runPsALC() {
	initWtime = MPI_Wtime();

	if(mpi_rank==0){
		cout<<"======================="<<endl;
		cout<<"PSEUDO-ARCLENGTH SOLVER"<<endl;
		cout<<"ikj AT stanford DOT edu"<<endl;
		cout<<"======================="<<endl;
	}
	string correctorMethod = getStringParam("CORRECTOR_METHOD", "NEWTON");
	string booleanString = getStringParam("AFTER_TIMESTEPPING_LAUNCH_NT", "FALSE");  // This will be used for TIME_STEPPING
	std::transform(booleanString.begin(), booleanString.end(), booleanString.begin(), ::tolower);
	bool AfterTSlaunchNT = (booleanString.compare("true")==0 || booleanString.compare("yes")==0);

	int nVars = 5+nScal;

	if (!checkParam("WEIGHT_LAMBDA_IN_ARCLENGTH")) {
		ParamMap::add("WEIGHT_LAMBDA_IN_ARCLENGTH=1.0"); // add default values
		if (mpi_rank == 0)
			cout << "WARNING: added keyword \"WEIGHT_LAMBDA_IN_ARCLENGTH=1.0\" to parameter map!"<<endl;
	}
	double weightLambda = getDoubleParam("WEIGHT_LAMBDA_IN_ARCLENGTH", "1.0");
	if(mpi_rank==0 && debugLevel>0) cout<<"> WEIGHT_LAMBDA_IN_ARCLENGTH = "<<weightLambda<<endl;

	// Arclength control algorithm
	bool controlArclength = false;
	double initResidNorm_opt;
	int iterNewton_opt;
	if (!checkParam("ARCLENGTH_CONTROL")) {
		ParamMap::add("ARCLENGTH_CONTROL  USE=NO  INIT_RESID_OPT=1.0  ITER_NEWTON_OPT=10  MIN_ARCLENGTH=0.0  MAX_ARCLENGTH=1.0e22"); // add default values
		if (mpi_rank == 0)
			cout<< "WARNING: added keyword \"ARCLENGTH_CONTROL  USE=NO  INIT_RESID_OPT=1.0  ITER_NEWTON_OPT=10  MIN_ARCLENGTH=0.0  MAX_ARCLENGTH=1.0e22\""<< " to parameter map!" << endl;
	}
	string tempString = getParam("ARCLENGTH_CONTROL")->getString("USE");
	if(mpi_rank==0 && debugLevel>0) cout<<"> Use ARCLENGTH_CONTROL = "<<tempString;
	std::transform(tempString.begin(), tempString.end(), tempString.begin(), ::toupper);
	if(tempString.compare("YES")==0 || tempString.compare("Y")==0) {
		controlArclength  = true;
		initResidNorm_opt = getParam("ARCLENGTH_CONTROL")->getDouble("INIT_RESID_OPT");
		iterNewton_opt    = getParam("ARCLENGTH_CONTROL")->getInt("ITER_NEWTON_OPT");
		if(mpi_rank==0 && debugLevel>0) cout<<"  INIT_RESID_OPT="<<initResidNorm_opt<<"  ITER_NEWTON_OPT="<<iterNewton_opt;
	}
	if(mpi_rank==0 && debugLevel>0) cout<<endl;
	double minArclength = getParam("ARCLENGTH_CONTROL")->getDouble("MIN_ARCLENGTH");
	double maxArclength = getParam("ARCLENGTH_CONTROL")->getDouble("MAX_ARCLENGTH");

	// Parameters related to Newton's method
	if(correctorMethod.compare("NEWTON") == 0 || (correctorMethod.compare("TIME_STEPPING")==0 && AfterTSlaunchNT)) {
		// Interval to check the convergence of PETSc KSP
		if (!checkParam("PETSC_CONVERG_MONITOR_INTERVAL")) {
			ParamMap::add("PETSC_CONVERG_MONITOR_INTERVAL=0"); // add default values
			if (mpi_rank == 0)
				cout << "WARNING: added keyword \"PETSC_CONVERG_MONITOR_INTERVAL=0\" to parameter map!"<<endl;
		}
		monitorConvergInterval = getIntParam("PETSC_CONVERG_MONITOR_INTERVAL", "0");

		// The restart value for GMRes
		gmresRestart = getIntParam("GMRES_RESTART", "100");
	}

	/*****
	 ** Allocation of memory
	 *****/
	//---------------------
	// Arrays and variables: PsALC is basically a predictor-corrector method.
	//                       In the predictor step, two previous solutions, (lambda0,q0) and (lambda1,q1) are to be used to get Qguess.
	//                       In the corrector step, Qguess is corrected by using the Newton's method
	// Note: The arrays will be cleaned by calling the freePsALCMemory() function:
	//       (q0, q1, Qguess, lambda0, lambda1, lambdaInit0, lambdaInit1, rhs1DArray, q_tangent, lambda_tangent)

	NcontrolEqns = getNumControlVar(); 	assert(NcontrolEqns >= 0);
	double *lambda0=NULL, *lambda1=NULL;
	double *lambdaInit0=NULL, *lambdaInit1=NULL;
	if(NcontrolEqns > 0) {
		lambda  = new double [NcontrolEqns]; // lambda is the input parameter to the system (e.g. fuel equivalence ratio)
							 				 // lambda is defined in the UgpWithCvCompFlow_AD class (for details, see IkeUgpWithCvCompFlow.h)
											 // lambda_AD will be allocated in the IkeWithPsALC_AD::initialize_adoubles() method
		lambda0 = new double [NcontrolEqns];
		lambda1 = new double [NcontrolEqns];
		lambdaInit0 = new double [NcontrolEqns];
		lambdaInit1 = new double [NcontrolEqns];
	} else {
		if(mpi_rank==0) cout<<"WARNING in IkeWithPsALC_AD::runPsALC(): NcontrolEqns = 0! The PsALC program might crush at some point."<<endl;
		lambda  = NULL;
		lambda0 = NULL;
		lambda1 = NULL;
		lambdaInit0 = NULL;
		lambdaInit1 = NULL;
	}

	double arclength              = -ABSURDLY_BIG_NUMBER;
	double normSqRatioInArclength = -ABSURDLY_BIG_NUMBER;

	double *q0 = NULL; // q0 is the 1D solution array at the (n-1)th step
	if(mpi_rank == mpi_size-1) 
		q0 = new double[ncv*nVars+NcontrolEqns];
	else 
		q0 = new double[ncv*nVars];
	double *q1 = NULL; // q1 is the 1D solution array at the (n)th step
	if(mpi_rank == mpi_size-1) 
		q1 = new double[ncv*nVars+NcontrolEqns];
	else
		q1 = new double[ncv*nVars];
	
	double *Qguess = NULL; // Qguess is the 1D solution array at the (n+1)th step
	if(mpi_rank == mpi_size-1)
		Qguess = new double[ncv*nVars+NcontrolEqns]; // note: Qguess = [q_guess; lambda_guess], where q_quess is the state vector
	else
		Qguess = new double[ncv*nVars];

	// Initialize arrays: use a very small number for a debugging purpose
	for(int i=0; i<ncv*nVars; ++i) {
		q0[i]     = -ABSURDLY_BIG_NUMBER;
		q1[i]     = -ABSURDLY_BIG_NUMBER;
		Qguess[i] = -ABSURDLY_BIG_NUMBER;
	}
	if(mpi_rank == mpi_size-1) {
		for(int iParam=0; iParam<NcontrolEqns; ++iParam) {
			q0[ncv*nVars+iParam]     = -ABSURDLY_BIG_NUMBER;
			q1[ncv*nVars+iParam]     = -ABSURDLY_BIG_NUMBER;
			Qguess[ncv*nVars+iParam] = -ABSURDLY_BIG_NUMBER;
		}
	}

	// Memory related to Newton's method
	if(correctorMethod.compare("NEWTON") == 0 || (correctorMethod.compare("TIME_STEPPING")==0 && AfterTSlaunchNT)) {
		// RHS of the N-S system: "rhs1DArray" is defined as a member variable since it can be very useful for debugging if it can be seen outside of this method
		assert(rhs1DArray==NULL);
		if(getStringParam("HOW_TO_CALC_JAC")=="ORDINARY_2D") {
			if(mpi_rank == mpi_size-1) {
				rhs1DArray = new double[ncv_gg*nVars+NcontrolEqns]; // Since the ordinary method was developed for Adjoint studies, it is also considering the ghost cells (ncv_gg). For more details, contact with Dr. Duraisamy
				for(int iParam=0; iParam<NcontrolEqns; ++iParam)
					rhs1DArray[ncv_gg*nVars+iParam] = -ABSURDLY_BIG_NUMBER; // Initialize with a very small number for a debugging purpose
			} else
				rhs1DArray = new double[ncv_gg*nVars]; // Since the ordinary method was developed for Adjoint studies, it is also considering the ghost cells. For more details, contact with Dr. Duraisamy
		} else if (getStringParam("HOW_TO_CALC_JAC")=="ROW_1D") {
			if(mpi_rank == mpi_size-1) {
				rhs1DArray = new double[ncv*nVars+NcontrolEqns];
				for(int iParam=0; iParam<NcontrolEqns; ++iParam)
					rhs1DArray[ncv*nVars+iParam] = -ABSURDLY_BIG_NUMBER; // Initialize with a very small number for a debugging purpose
			} else
				rhs1DArray = new double[ncv*nVars];
		} else {
			if(mpi_rank==0)
				cout<<"ERROR in runPsALC(): HOW_TO_CALC_JAC="<<getStringParam("HOW_TO_CALC_JAC")<<" cannot be supported"<<endl;
			return;
		}
	}

	// Tangential vectors
	double **q_tangent     = NULL;
	double *lambda_tangent = NULL;
	if(NcontrolEqns > 0) {
		getMem2D(&q_tangent, 0, NcontrolEqns-1, 0, (ncv*nVars)-1, "q_tangent"); // The size of q_tangent is the same for both ORDINARY_2D and ROW_1D: In the matrix used for Newton's method, you don't need to have ncv_gg even for ORDINARY_2D
		for(int iParam=0; iParam<NcontrolEqns; ++iParam)
			for(int i=0; i<ncv*nVars; ++i)
				q_tangent[iParam][i] = 0.0;

		lambda_tangent = new double [NcontrolEqns];
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
			lambda_tangent[iEqn] = 0.0;
	}

	// RHS weight
	weightRhs = new double [ncv*(5+nScal)];

	// For RHS calculations and tecplot output of residuals
//	assert(RHSrho==NULL && RHSrhou==NULL && RHSrhoE==NULL);
//	registerScalar(RHSrho,  "RHSRHO",  CV_DATA);
//	registerVector(RHSrhou, "RHSRHOU", CV_DATA);
//	registerScalar(RHSrhoE, "RHSRHOE", CV_DATA);
	assert(RHSrhoScal == NULL);
	if (nScal > 0) getMem2D(&RHSrhoScal,  0, nScal-1, 0, ncv_g-1, "rhsScal");

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

	//---------------------
	// QoI
	string QoIfilename = "QoI.txt";

	//---------------------
	// Status file
	if(correctorMethod.compare("NEWTON") == 0 || (correctorMethod.compare("TIME_STEPPING")==0 && AfterTSlaunchNT)) {
		if(mpi_rank==0) {
			FILE *fp = fopen(JAC1D_STATUS_FILENAME, "w");
			fprintf(fp, "STEP   NEWTON_ITER         ICV   WTIME[sec]\n");
			fclose(fp);
		}
	}

	//---------------------
	// step
	startingStep = getIntParam("STARTING_STEP", "0");
	step = startingStep;

	bool done = false;
	if (nsteps-step < 2) {
		if(mpi_rank==0)
			cout<<"IkeWithPsALC_AD::runPsALC() -- code will NOT produce any result since NSTEPS(="<<nsteps-step<<") < 2"<<endl;
		done = true; // Note: We need at least two initial points (say, step=0 and step=1) to calculate the the tangential vector.
	}

	/*****
	 ** Main body
	 *****/
	//---------------------
	// wall time
	double wtime;
	MPI_Barrier( mpi_comm);
	if (mpi_rank == 0)
		wtime = MPI_Wtime();

	//---------------------
	// get (lambda0, q0) and (lambda1, q1)
	try {
		// Get q0, q1, lambda0, lambda1, and arclength
		initFirstTwoPts(q0, q1, rhs1DArray, lambdaInit0, lambdaInit1, NcontrolEqns, arclength, &normSqRatioInArclength, weightLambda, QoIfilename);
			// Note1: After calling this function, "step" becomes step+2
			// Note2: initFirstTwoPts() also write a tecplot file and a restart file
			// Note3: "arclength" is also calculated in this method since arclength is stored in the binary file for the 2nd bifurcation point

		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
			lambda0[iEqn] = lambdaInit0[iEqn];
			lambda1[iEqn] = lambdaInit1[iEqn];
		}
	}
	catch (int e) {
		freePsALCMemory(&q0, &q1, &Qguess, &lambda0, &lambda1, &lambdaInit0, &lambdaInit1, &lambda_tangent, nVars);
		if(q_tangent != NULL) {
			freeMem2D(q_tangent, 0, NcontrolEqns-1, 0, (ncv*nVars)-1); 	q_tangent = NULL;
		}
		throw(e);
	}
	catch (...) {
		freePsALCMemory(&q0, &q1, &Qguess, &lambda0, &lambda1, &lambdaInit0, &lambdaInit1, &lambda_tangent, nVars);
		if(q_tangent != NULL) {
			freeMem2D(q_tangent, 0, NcontrolEqns-1, 0, (ncv*nVars)-1); 	q_tangent = NULL;
		}
		throw;
	}

	// Reduce the arclength size if required
	double factorReduceTangential = getDoubleParam("REDUCE_ARCLENTH_THIRD", "1.0"); // If you know that the pseudo-arclength continuation method cannot converge at point 3,
	                                                                                // you may want to restart the continuation from 2 and 1 but reduce the arclength.
	assert(factorReduceTangential>0.0);
	if(mpi_rank==0) {
		cout<<endl;
		if(factorReduceTangential<1.0-MACHINE_EPS)
			cout<<">> ARCLENGTH IS REDUCED BY FACTOR OF "<<factorReduceTangential<<" ("<<arclength<<" --> "<<arclength*factorReduceTangential<<")"<<endl;
		else if(factorReduceTangential>1.0+MACHINE_EPS)
			cout<<">> ARCLENGTH IS INCREASED BY FACTOR OF "<<factorReduceTangential<<" ("<<arclength<<" --> "<<arclength*factorReduceTangential<<")"<<endl;
	}
	arclength *= factorReduceTangential;

	// Show the arclength on the screen
	if(mpi_rank==0) {
		// Calculate the flow part and lambda part from the inner product to help the user to decide WEIGHT_FOR_LAMBDA
		double lambdaProdUnweight = 0.0;
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
			lambdaProdUnweight += (lambda1[iEqn] - lambda0[iEqn])*(lambda1[iEqn] - lambda0[iEqn]); // note: lambda_tangent has not been assigned and is still zero

		double qProd = sqrt(arclength * arclength - weightLambda * lambdaProdUnweight);
		lambdaProdUnweight = sqrt(lambdaProdUnweight);

		printf("\n>> ARCLENGTH TO BE USED DURING THE CONTINUATION = %.6e\n", arclength);
		cout<<"   (ARCLENGTH = "<<arclength<<" = sqrt( (FLOW_PART="<<qProd<<")^2 + WEIGHT_FOR_LAMBDA * (LAMBDA_PART="<<lambdaProdUnweight<<")^2 ),"<<endl;
		cout<<"    WEIGHT FOR LAMBDA = "<<weightLambda<<","<<endl;
		cout<<"    RATIO BETWEEN THE NORM OF THE FLOW VARIABLES AND THE LAMBDAS (BEFORE WEIGHT_FOR_LAMBDA) = "<< sqrt(normSqRatioInArclength) <<")"<<endl;
	}

	// Set the R.H.S. weighting for the tangential equation
	if (!checkParam("AREA_LAMBDA")) {
		ParamMap::add("AREA_LAMBDA  METHOD=INVERSE_ARCLENGTH  COEFF=100.0  CUTOFF_VAL=1.0"); // add default values
		if (mpi_rank == 0)
			cout<< "WARNING: added keyword \"AREA_LAMBDA  METHOD=INVERSE_ARCLENGTH  COEFF=100.0  CUTOFF_VAL=1.0\""<< " to parameter map!" << endl;
	}
	string areaLambda_method = getParam("AREA_LAMBDA")->getString("METHOD");
	double areaLambda_coeff  = getParam("AREA_LAMBDA")->getDouble("COEFF");

	if(areaLambda_method.compare("FIXED_VALUE") == 0)
		AreaLambda = areaLambda_coeff;
	else {
		double areaLambda_cutoff = getParam("AREA_LAMBDA")->getDouble("CUTOFF_VAL");

		AreaLambda = min(areaLambda_cutoff, areaLambda_coeff/arclength);
	}

	if(mpi_rank == 0)
		cout<<endl
		    <<">> AreaLambda: Method = "<<areaLambda_method<<", Value = "<<AreaLambda<<endl;

	//---------------------
	// Pseudo-arclength continuation
	if (nsteps-step < 0) {
		if(mpi_rank==0)
			cout<<endl<<"IkeWithPsALC_AD::runPsALC() -- code will not run bifurcation since NSTEPS(="<<nsteps<<") < STEP(="<<step<<")+1"<<endl;
		done = true;
	}

	string tempFor3rdPt = getStringParam("INIT_THIRD_FROM_Q1", "NO");
	std::transform(tempFor3rdPt.begin(), tempFor3rdPt.end(), tempFor3rdPt.begin(), ::tolower);
	if(tempFor3rdPt.compare("yes")==0)
		initThirdFromQ1 = true;
	else
		initThirdFromQ1 = false;

	PsALCinitialHOOK_debug(); // This function is only for debugging: the user can have a opportunity to look into the variables

	double wtime0, wtimeF;
	while (!done) {
		// 1. one step pseudo-arclength continuation (updating Qguess) with a backtracking algorithm
		if(mpi_rank == 0 && debugLevel>0)
			cout<<endl
				<<"Starting PsALC step = "<<step<<endl;

		// 1.1. calculate "Qguess" by running single step pseudo-arclength continuation
		if(mpi_rank==0)
			wtime0 = MPI_Wtime();

		try {
			bool writeJacMat = true;
			if(correctorMethod.compare("TIME_STEPPING") == 0)
				writeJacMat = false;

			singlePsALCstep(Qguess, q0, q1, rhs1DArray, q_tangent, lambda_tangent, arclength, lambda0, lambda1, weightLambda, writeJacMat, factorReduceTangential);
				// Note: Even though factorReduceTangential != 1.0, [q_tangent; lambda_tangent] must be a unit vector.
				//       Thus, [q_tangent; lambda_tangent] is calculated by 1/ds * factorReduceTangential * ([q1; lambda1] - [q0; lambda0])
				//      in singlePsALCstep(), where ds is the arc-length that was already reduce by factorReduceTangential
				//      (i.e. ds = arclength*factorReduceTangential).
		}
		catch (int e) {
			freePsALCMemory(&q0, &q1, &Qguess, &lambda0, &lambda1, &lambdaInit0, &lambdaInit1, &lambda_tangent, nVars);
			if(q_tangent != NULL) {
				freeMem2D(q_tangent, 0, NcontrolEqns-1, 0, (ncv*nVars)-1); 	q_tangent = NULL;
			}
			throw(PSALC_ERROR_CODE);
		}
		catch (...) {
			freePsALCMemory(&q0, &q1, &Qguess, &lambda0, &lambda1, &lambdaInit0, &lambdaInit1, &lambda_tangent, nVars);
			if(q_tangent != NULL) {
				freeMem2D(q_tangent, 0, NcontrolEqns-1, 0, (ncv*nVars)-1); 	q_tangent = NULL;
			}
			throw(PSALC_ERROR_CODE);
		}

		// 1.2. update q and lambda
		for(int i=0; i<ncv*(nScal+5); ++i) {
			q0[i] = q1[i];
			q1[i] = Qguess[i];
		}
		if(mpi_rank == mpi_size-1) {
			for(int iParam=0; iParam<NcontrolEqns; ++iParam) {
				q0[ncv*nVars+iParam] = q1[ncv*nVars+iParam];
				q1[ncv*nVars+iParam] = Qguess[ncv*nVars+iParam];
			}
		}

		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
			lambda0[iEqn] = lambda1[iEqn];
		MPI_Barrier(mpi_comm);
		if(mpi_rank == mpi_size-1) {
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {

				lambda1[iEqn] = Qguess[ncv*(nScal+5+iEqn)];
			}
		}
		MPI_Bcast(lambda1, NcontrolEqns, MPI_DOUBLE, mpi_size-1, mpi_comm );

		// 1.3. recalculate the arc-length
		double normSqRatioInArclength; // Actually, this variable will not be used for the current version
		double dsNew = calcArclength(&normSqRatioInArclength, q0, q1, lambda0, lambda1, NcontrolEqns, nScal, weightLambda);
			// Note: You should not pass q_tangent & lambda_tangent to the method in this case, because they are old (not updated).

		// 1.4. report the result on the screen
		if(mpi_rank==0) {
			wtimeF = MPI_Wtime();
			printf("> Summary of PSALC STEP=%3d: LAMBDA=%13.6e  RUNTIME[sec]=%.2f  ARCLENGTH=%.6e\n", step, lambda1[0], wtimeF-wtime0, dsNew);
		}

		// 1.5. warning error if arclength is not consistent or too small
		double arclengthError = fabs((arclength-dsNew)/arclength);
//		if(arclengthError > ARCLENGTH_ERROR_EPS) {
//			if(mpi_rank==0)
//				printf("ERROR! IkeWithPsALC_AD::runPsALC(): New ds(=%.8e) is not equal to original arclength(=%.8e) -- error=%.6f% \n", dsNew, arclength, arclengthError*100);
//			freePsALCMemory(&q0, &q1, &Qguess, &lambda0, &lambda1, &lambdaInit0, &lambdaInit1, &lambda_tangent, nVars);
//			if(q_tangent != NULL) {
//				freeMem2D(q_tangent, 0, NcontrolEqns-1, 0, (ncv*nVars)-1); 	q_tangent = NULL;
//			}
//			throw(PSALC_ERROR_CODE);
//		}
		if(arclengthError > ARCLENGTH_WARNING_EPS) {
			if(mpi_rank==0)
				printf("WARNING! IkeWithPsALC_AD::runPsALC(): New ds(=%.6e) is not equal to original ds(=%.6e) -- difference=%.4f\n", dsNew, arclength, arclengthError*100);
		}
		if(dsNew < MACHINE_EPS) {
			if(mpi_rank==0)
				cerr<<"ERROR in IkeWithPsALC_AD::runPsALC(): New ds(="<<dsNew<<") is too small"<<endl;
			freePsALCMemory(&q0, &q1, &Qguess, &lambda0, &lambda1, &lambdaInit0, &lambdaInit1, &lambda_tangent, nVars);
			if(q_tangent != NULL) {
				freeMem2D(q_tangent, 0, NcontrolEqns-1, 0, (ncv*nVars)-1); 	q_tangent = NULL;
			}
			throw(PSALC_ERROR_CODE);
		}

		// 1.6. Update the arclength
		arclength = dsNew;

		// 2. write a tecPlot file
		writeData(step);

		// 3. write restart
		if (write_restart > 0) {
			if(step%write_restart==0) {
				writeRestart(step); // Write the current flow field as a restart file
				char filename[40];
				sprintf(filename, "Q1_PT%05d.bin", step);
//				writePsALCdumpedDataSerial(filename, step, arclength, lambda0, lambda1, NcontrolEqns, q1);
				writePsALCdumpedDataParallel(filename, step, arclength, lambda0, lambda1, NcontrolEqns, q1);
			} else if(write_restart>1 && (step+1)%write_restart==0) {
				writeRestart(step);
				char filename[40];
				sprintf(filename, "Q1_PT%05d.bin", step);
//				writePsALCdumpedDataSerial(filename, step, arclength, lambda0, lambda1, NcontrolEqns, q1);
				writePsALCdumpedDataParallel(filename, step, arclength, lambda0, lambda1, NcontrolEqns, q1);
			}
		}

		// 4. calculate QoI and dump it on a file
		rewriteQoI = false;
		writeQoIOnFile(step, QoIfilename, NcontrolEqns, rewriteQoI);

		// 5. launch temporalHook
		temporalHook();

		// 6. check if the continuation should stop
		step++;
		if (step > nsteps)
			done = true;

		// 7. clear the PETSc object
		if(petscSolver2 != NULL) {
			delete petscSolver2; 	petscSolver2 = NULL;
		}

		// 8. Update the arclength
		if(controlArclength) {
			if(step-1-startingStep==2 && initThirdFromQ1) { // If the 3rd point was initialized by a binary file in singlePsALCstep(),
				                                            // you cannot use current residual and iterNewton to decide the next arc-length.
				double startingResidual = getDoubleParam("INIT_THIRD_RESID", "0.0");
				if(mpi_rank==0) cout<<"Note: Since the initial guess for the 3rd point was a binary file, use the user-provided INIT_THIRD_RESID for arclength-control"<<endl;

				factorReduceTangential = arclengthControl(startingResidual, initResidNorm_opt, iterNewton, iterNewton_opt);
						// Note: iterNewton was already added by INIT_THIRD_NEWTON_ITERS in getSteadySolnByNewton(), so you don't have to update it here.
			} else
				factorReduceTangential = arclengthControl(initResidNorm,    initResidNorm_opt, iterNewton, iterNewton_opt);
						// Note: "initResidNorm" and "iterNewton" are Member variables

			if(factorReduceTangential < minArclength/arclength) {
				if(mpi_rank==0) cout<<"WARNING in runPsALC(): too small arclength after applying arclength-control (new arclength="<<factorReduceTangential*arclength<<" < MIN_ARCLENGTH="<<minArclength<<") -- CLIP TO "<<minArclength<<endl;
				factorReduceTangential = minArclength/arclength;
			} else if(factorReduceTangential > maxArclength/arclength) {
				if(mpi_rank==0) cout<<"WARNING in runPsALC(): too big arclength after applying arclength-control (new arclength="<<factorReduceTangential*arclength<<" < MAX_ARCLENGTH="<<maxArclength<<") -- CLIP TO "<<maxArclength<<endl;
				factorReduceTangential = maxArclength/arclength;
			}

			if(mpi_rank==0)
				printf("Change arclength by %.6f percents: %.4e --> %.4e\n", factorReduceTangential*100, arclength, arclength*factorReduceTangential);

			arclength *= factorReduceTangential;
		} else
			factorReduceTangential = 1.0;

		MPI_Barrier(mpi_comm);
	}

	/****/
	/* show the total runtime */
	/****/
	MPI_Barrier(mpi_comm);
	if (mpi_rank == 0) {
		double wtime0 = wtime;
		wtime = MPI_Wtime();
		cout << " > runtime for iterations[sec]: " << wtime - wtime0 << endl;
	}

	finalHook();

	/*****
	 ** Finalize
	 *****/
	MPI_Barrier(mpi_comm);
	if (mpi_rank == 0) {
		double wtime0 = wtime;
		wtime = MPI_Wtime();
		cout << " > runtime for iterations[sec]: " << wtime - wtime0 << endl;
	}

	finalHook();

	freePsALCMemory(&q0, &q1, &Qguess, &lambda0, &lambda1, &lambdaInit0, &lambdaInit1, &lambda_tangent, nVars);
	if(q_tangent != NULL) {
		freeMem2D(q_tangent, 0, NcontrolEqns-1, 0, (ncv*nVars)-1); 	q_tangent = NULL;
	}
}

/*
 * Method: freePsALCMemory
 * -----------------------
 * Clear the arrays used in the runPsALC class.
 * This method is called not only at the end of the runPsALC class but also whenever an exception occurs.
 */
void IkeWithPsALC_AD::freePsALCMemory(double **p_q0, double **p_q1, double **p_Qguess, double **p_lambda0, double **p_lambda1,
		double **p_lambdaInit0, double **p_lambdaInit1, double **p_lambda_tangent, const int nVars) {
	if(debugLevel>0 && mpi_rank==0)
		cout<<"IkeWithPsALC_AD::freePsALCMemory(): Free memory allocated for the variables in runPsALC()..."<<endl;

	if((*p_q0) != NULL) {
		delete [] (*p_q0); 			(*p_q0) = NULL;
	}
	if((*p_q1) != NULL) {
		delete [] (*p_q1); 			(*p_q1) = NULL;
	}
	if((*p_Qguess) != NULL) {
		delete [] (*p_Qguess); 		(*p_Qguess) = NULL;
	}
	if((*p_lambda0) != NULL) {
		delete [] (*p_lambda0);	 	(*p_lambda0) = NULL;
	}
	if((*p_lambda1) != NULL) {
		delete [] (*p_lambda1); 	(*p_lambda1) = NULL;
	}
	if((*p_lambdaInit0) != NULL) {
		delete [] (*p_lambdaInit0); 	(*p_lambdaInit0) = NULL;
	}
	if((*p_lambdaInit1) != NULL) {
		delete [] (*p_lambdaInit1); 	(*p_lambdaInit1) = NULL;
	}

	if(rhs1DArray != NULL) {
		delete [] rhs1DArray; 	rhs1DArray = NULL;
	}

	if (RHSrhoScal != NULL && nScal > 0)
		freeMem2D(RHSrhoScal, 0, nScal-1, 0, ncv-1);

	if((*p_lambda_tangent) != NULL) {
		delete [] (*p_lambda_tangent); 	(*p_lambda_tangent) = NULL;
	}

	if(!jacMatrix.empty())
		jacMatrix.clear();
	if(!jacMatrixSTL.empty())
		jacMatrixSTL.clear();

	MPI_Barrier(mpi_comm);
}

/*
 * Method: initFirstTwoPts
 * -----------------------
 * Obtain the first point and second point in the bifurcation curve.
 * This method has the option to call HOOK functions for both of the two points.
 *
 * Return:
 *   By reference:
 *     q0
 *     q1
 *     rhs (depends on the initialization method speficied in Ike.in)
 *     lambdaInit0
 *     lambdaInit1
 *     arclength
 *     normSqRatioInArclength = || q1-q0 ||_{W}^2 / || lambdaInit1-lambdaInit0 ||_{W}^2
 */
void IkeWithPsALC_AD::initFirstTwoPts(double* q0, double* q1, double *rhs, double *lambdaInit0, double *lambdaInit1,
		const int NcontrolEqns, double &arclength, double *normSqRatioInArclength, const double weightLambda, string &QoIfilename) {
	assert(q0!=NULL && q1!=NULL);
	int m = nScal+5; // number of variables (rho, rhou, rhov, rhow, rhoE, and scalars)
	int kine_Index = getScalarTransportIndex("kine");
	if(kine_Index>=0 && UgpWithCvCompFlow::kine==NULL)
		UgpWithCvCompFlow::kine = scalarTranspEqVector[kine_Index].phi;
	initialHookScalarRansTurbModel(); // Calculate wallDist, connect the member variables (e.g. kine, omega, grad_kine, grad_omega) in the turbulence model to arrays
	initialHookScalarRansCombModel(); // Update gamma and RoM

	initHookPsALC();

	// ------------
	// First point
	// ------------
	int howToInitFirst = getIntParam("HOW_TO_INIT_FIRST", "0");
	if(mpi_rank==0) cout<<"> You specify HOW_TO_INIT_FIRST = "<<howToInitFirst<<": ";
	string filename;
	char paramName[40];

	double lambdaInit0temp[NcontrolEqns]; // Used only if howToInitFirst == 2

	switch(howToInitFirst) {
	case 1 : // from a hook function
		if(mpi_rank==0) cout<<" Initialize q0 from initHookPsALC_1st()"<<endl;
		initHookPsALC_1st();
		update1DVecFromPrimVars(q0, nScal);

		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
			sprintf(paramName, "CONTROL_PARAM%d", iEqn);
			lambdaInit0[iEqn] = getParam("LAMBDA_INITIAL_0") -> getDouble(paramName);
		}
		break;
	case 2 : // from a previous PsALC simulation by IKE
		if(mpi_rank==0) cout<<" Initialize q0 from a binary file dumped by IKE";
		filename = getStringParam("Q0_FILENAME", "Q1_PT00000.bin");
		if(mpi_rank==0) cout<<" - filename = "<<filename;

		readPsALCdumpedDataSerial(filename, q0, step, lambdaInit0temp, lambdaInit0, NcontrolEqns, arclength);

		if(mpi_rank==0) {
			printf("  lambda0=");
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
				printf("%.3e,", lambdaInit0temp[iEqn]);
			printf("  lambda1=");
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
				printf("%.3e,", lambdaInit0[iEqn]);
			printf("  arclength=%.3e \n", arclength);
		}

		updatePrimVarsFrom1DVec(q0, nScal);
		break;
	case 3 : // from a bindary file dumped by normal JOE
		if(mpi_rank==0) cout<<" Initialize q0 from a binary file dumped by JOE";
		filename = getStringParam("Q0_FILENAME", "finalData.bin");
		if(mpi_rank==0) cout<<" - filename ="<<filename<<endl;

		readJOEdumpedDataSerial(filename, q0);

		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
			sprintf(paramName, "CONTROL_PARAM%d", iEqn);
			lambdaInit0[iEqn] = getParam("LAMBDA_INITIAL_0") -> getDouble(paramName);
		}
		if(mpi_rank==0) {
			printf("  lambda0=");
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
				printf("%.3e ", lambdaInit0[iEqn]);
			cout<<endl;
		}

		updatePrimVarsFrom1DVec(q0, nScal);
		break;
	default : // from a restart file
		if(mpi_rank==0) cout<<" Initialize q0 from a restart file";
		assert(checkScalarFlag("RHO")); // the restart file should exist

		update1DVecFromPrimVars(q0, nScal);

		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
			sprintf(paramName, "CONTROL_PARAM%d", iEqn);
			lambdaInit0[iEqn] = getParam("LAMBDA_INITIAL_0") -> getDouble(paramName);
		}
		if(mpi_rank==0) {
			printf(" - lambda0=");
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
				printf("%.3e ", lambdaInit0[iEqn]);
			cout<<endl;
		}
		break;
	}
	// Update ghost cells in the 2nd-layer
    updateCvDataG1G2(rho, REPLACE_DATA);
    updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
    updateCvDataG1G2(rhoE, REPLACE_DATA);
    for (int iScal = 0; iScal < nScal; iScal++)
      updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);

	calcStateVariables();     // Compute velocity, pressure, temperature, enthalpy and speed of sound at cell centers
//	calcMaterialProperties(); // Compute viscosity and thermal diffusivity at cell faces
	setBC();
	for(int ifa=nfa; ifa<nfa_b2; ++ifa)
		setBC1D(ifa, UgpWithCvCompFlow::rho, UgpWithCvCompFlow::rhou, UgpWithCvCompFlow::rhoE); // Note: ncvggf's cannot be not updated, so the user should take care of it
	calcRansTurbViscMuet();

	for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
		lambda[iEqn] = lambdaInit0[iEqn]; // This is required since lambda will be written on QoI.txt
		if(mpi_rank==mpi_size-1)
			q0[ncv*(5+nScal)+iEqn] = lambda[iEqn];
	}

	temporalHook();
	MPI_Barrier(mpi_comm);

	rewriteQoI = true;
	writeQoIOnFile(step, QoIfilename, NcontrolEqns, rewriteQoI);

//	initHookPsALC();

//char filenameT[40];
//sprintf(filenameT, "coordInfo%05d.bin\0", step);
//write1DArrayWithXcvSerial(filenameT, step, q0);
//MPI_Barrier(mpi_comm);

//	// Solve the system by Backward Euler: just for fun
//	{
//		int maxStepBEuler   = 0;
//
//		double absTolBEuler = 100.0*((double) cvora[mpi_size])*MACHINE_EPS;
//		double relTolBEuler = 1.0e-7;
//		int maxIterLS_BEuler;
//		double zeroAbsLS_BEuler, zeroRelLS_BEuler;
//		double CflBEuler;
//		int startIncCFLBEuler, intervalIncCFLBEuler;
//		double incCFLBEuler, maxCFLBEuler;
//		getParamsBEuler(maxIterLS_BEuler, zeroAbsLS_BEuler, zeroRelLS_BEuler,
//				CflBEuler, startIncCFLBEuler, intervalIncCFLBEuler, incCFLBEuler, maxCFLBEuler);
//
//		double (*dq_BEuler)[5]   = new double[ncv_g][5];
//		double (*rhs_BEuler)[5]  = new double[ncv][5];
//		double (*A_BEuler)[5][5] = new double[nbocv_s][5][5];
//		double **dScal_BEuler   = NULL;  if (nScal > 0) getMem2D(&dScal_BEuler,   0, nScal-1, 0, ncv_g-1, "dScal_BEuler");
//		double **rhsScal_BEuler = NULL;  if (nScal > 0) getMem2D(&rhsScal_BEuler, 0, nScal-1, 0, ncv-1,   "rhsScal_BEuler");
//		double ***AScal_BEuler  = NULL;  if (nScal > 0) getMem3D(&AScal_BEuler,   0, nScal-1, 0, 5, 0, nbocv_s-1, "AScal_BEuler");
//		double *Residual_BEuler = new double[5+nScal];
//		for(int i=0; i<5+nScal; ++i)
//			Residual_BEuler[i] = 0.0;
//
//		if(debugLevel>0) {
//			if(mpi_rank==0) {
//				cout<<endl
//					<<"  >> RUNNING BACKWARD_EULER for "<<maxStepBEuler<<" STEPS: absTol="<<absTolBEuler<<", relTol="<<relTolBEuler<<endl
//					<<      "     CFL_START="<<CflBEuler<<", CFL_MAX="<<maxCFLBEuler<<endl;
//			}
//		}
//		MPI_Barrier(mpi_comm);
//
//		double initTotalResidual_BEuler = 0.0;
//		bool isBEulerConverged = getSteadySolnByBackwardEuler(dq_BEuler, dScal_BEuler, rhs_BEuler, rhsScal_BEuler, A_BEuler, AScal_BEuler,
//						initTotalResidual_BEuler, Residual_BEuler,
//						maxStepBEuler, absTolBEuler, relTolBEuler, maxIterLS_BEuler, zeroAbsLS_BEuler, zeroRelLS_BEuler,
//						CflBEuler, startIncCFLBEuler, intervalIncCFLBEuler, incCFLBEuler, maxCFLBEuler, 1);
//
//		if(maxStepBEuler>0) {
//			if(mpi_rank==0)
//				cout<<"     BACKWARD_EULER TEST WAS DONE"<<endl
//					<<endl;
//			MPI_Barrier(mpi_comm);
//
//			if(debugLevel>0) {
//				int whichNorm = 1;
//				double totResidual = 0.0;
//				for(int i=0; i<5+nScal; ++i)
//					totResidual += Residual_BEuler[i];
//				if(mpi_rank==0) {
//					cout<<"  >> FINAL RESIDUAL FROM BACKWARD_EULER ("<<whichNorm<<"-norm): ";
//					if(isBEulerConverged)
//						printf("CONVERGED: InitResidual=%.5e --> FinalResidual=%.5e \n", initTotalResidual_BEuler, totResidual);
//					else
//						printf("NOT-CONVERGED: InitResidual=%.5e --> FinalResidual=%.5e \n", initTotalResidual_BEuler, totResidual);
//				}
//			}
//		}
//
//		delete [] dq_BEuler;
//		delete [] rhs_BEuler;
//		delete [] A_BEuler;
//		if (nScal > 0) freeMem2D(dScal_BEuler,   0, nScal-1, 0, ncv_g-1);
//		if (nScal > 0) freeMem2D(rhsScal_BEuler, 0, nScal-1, 0, ncv-1);
//		if (nScal > 0) freeMem3D(AScal_BEuler,   0, nScal-1, 0, 5, 0, nbocv_s-1);
//		delete [] Residual_BEuler;
//		if(petscSolver != NULL) {
//			delete petscSolver;
//			petscSolver = NULL;
//		}
//	}

	// Write a summary of the Bifurcation code on a file
	if(mpi_rank==0) {
		FILE *fp = fopen(BIFUR_SUMMARY_FILENAME, "w");
		fprintf(fp, "STEP,  ARC_LENGTH,  RESID_INIT, RESID_FINAL, TOT_ITER, DIFF_LAMBDA, LAMBDA \n");
		fprintf(fp, "%4d,            ,            ,            ,         ,            , %.15f\n", step, lambdaInit0[0]);
		fclose(fp);
	}

	// Write the initial field on a (tecplot) file
	if (initial_flowfield_output == "YES") {
		writeRestart(step); // Restart file
		writeData(step);    // Tecplot file
	}
	++step;

	// ------------
	// Second point
	// ------------
	arclength = -ABSURDLY_BIG_NUMBER;
	int howToInitSecond = getIntParam("HOW_TO_INIT_SECOND", "0");
	if(mpi_rank==0)
		cout<<endl
			<<"> You specify HOW_TO_INIT_SECOND = "<<howToInitSecond<<": ";
	char paramName2nd[20];
	double *lambda0Temp2nd=NULL, *lambda1Temp2nd=NULL;
	double errorAbs2nd = 0.0;
	double errorArclength = 0.0;

	// Update lambdaInit1, lambda
	switch(howToInitSecond) {
		case 1 : // from a hook function
			if(mpi_rank==0) cout<<" Initialize q1 from initHookPsALC_2nd()"<<endl;

        	for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
        		sprintf(paramName2nd, "CONTROL_PARAM%d", iEqn);
        		lambdaInit1[iEqn] = getParam("LAMBDA_INITIAL_1") -> getDouble(paramName2nd);
        		lambda[iEqn]      = lambdaInit1[iEqn]; // This is required since lambda will be written on QoI.txt

    			if(mpi_rank==mpi_size-1)
    				q1[ncv*(5+nScal)+iEqn] = lambda[iEqn];
        	}

            initHookPsALC_2nd();
			update1DVecFromPrimVars(q1, nScal);

			arclength = calcArclength(normSqRatioInArclength, q0, q1, lambdaInit0, lambdaInit1, NcontrolEqns, nScal, weightLambda);

			break;
		case 2 : // from a previous PsALC simulation by IKE
			if(mpi_rank==0) cout<<" Initialize q1 from a binary file dumped by IKE";
			filename = getStringParam("Q1_FILENAME", "Q1_PT00001.bin");
			if(mpi_rank==0) cout<<" - filename = "<<filename<<endl;

			lambda0Temp2nd = new double [NcontrolEqns];
			lambda1Temp2nd = new double [NcontrolEqns];
			double arclengthTemp2nd;

			readPsALCdumpedDataSerial(filename, q1, step, lambda0Temp2nd, lambda1Temp2nd, NcontrolEqns, arclengthTemp2nd);

            for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
				sprintf(paramName2nd, "CONTROL_PARAM%d", iEqn);
				lambdaInit1[iEqn] = lambda1Temp2nd[iEqn];
				lambda[iEqn]      = lambdaInit1[iEqn];

    			if(mpi_rank==mpi_size-1)
    				q1[ncv*(5+nScal)+iEqn] = lambda[iEqn];
			}

			if(mpi_rank==0) {
				for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
					printf("  lambda%d=%.3e", iEqn, lambda[iEqn]);
				printf("  arclength=%.4e \n", arclengthTemp2nd);
			}

			// Check possible errors
			if(mpi_rank==0) cout<<"  >> Checking the compatibility of the binary file to the simulation...";
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
				errorAbs2nd += fabs((lambda1Temp2nd[iEqn] - lambdaInit1[iEqn])/lambdaInit1[iEqn]);
			if(errorAbs2nd > LAMBDA_COMPATIBILITY_EPS) {
				if(mpi_rank==0) {
					printf("\n  Error! The binary file does not match with the Ike.in file \n");
					printf("       Binary file:");
					for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
						printf(" lambda%d=%.5e ", iEqn, lambdaInit1[iEqn]);
					printf("\n       Input file : ");
					for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
						printf(" lambda%d=%.5e ", iEqn, lambda1Temp2nd[iEqn]);
					cout<<endl;
				}
				throw(PSALC_ERROR_CODE);
			}

			// Calculate the arclength using q0, q1, lambdaInit0, and lambdaInit1
			arclength = calcArclength(normSqRatioInArclength, q0, q1, lambdaInit0, lambdaInit1, NcontrolEqns, nScal, weightLambda);

			errorArclength = fabs((arclength-arclengthTemp2nd)/arclength);
			if(errorArclength > sqrt(cvora[mpi_size])*MACHINE_EPS) {
				if(mpi_rank==0)
					printf("\n     WARNING! arclength - binary file=%.5e, calculation=%.5e, error=%.4e  \n", arclengthTemp2nd, arclength, errorArclength);
			} else {
				if(mpi_rank==0) cout<<" OK"<<endl;
			}

			// update the flow field
			updatePrimVarsFrom1DVec(q1, nScal);

			// free memory
			delete [] lambda0Temp2nd;
			delete [] lambda1Temp2nd;

			break;
		case 3 : // from a bindary file dumped by normal JOE
			if(mpi_rank==0) cout<<" Initialize q1 from a binary file dumped by JOE";
			filename = getStringParam("Q1_FILENAME", "finalData.bin");
			if(mpi_rank==0) cout<<" - filename = "<<filename<<endl;

			readJOEdumpedDataSerial(filename, q1);
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
				sprintf(paramName2nd, "CONTROL_PARAM%d", iEqn);
				lambdaInit1[iEqn] = getParam("LAMBDA_INITIAL_1") -> getDouble(paramName2nd);
                lambda[iEqn]      = lambdaInit1[iEqn];

    			if(mpi_rank==mpi_size-1)
    				q1[ncv*(5+nScal)+iEqn] = lambda[iEqn];
            }
			if(mpi_rank==0) {
				printf("  lambda1=");
				for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
					printf("%.3e ", lambda[iEqn]);
				cout<<endl;
			}

			arclength = calcArclength(normSqRatioInArclength, q0, q1, lambdaInit0, lambdaInit1, NcontrolEqns, nScal, weightLambda);

			updatePrimVarsFrom1DVec(q1, nScal);

			break;
		case 4 : // from a restart file
			if(mpi_rank==0) cout<<" Initialize q1 from a restart file";
			assert(checkScalarFlag("RHO")); // the restart file should exist
			update1DVecFromPrimVars(q1, nScal);

			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
				sprintf(paramName2nd, "CONTROL_PARAM%d", iEqn);
				lambdaInit1[iEqn] = getParam("LAMBDA_INITIAL_1") -> getDouble(paramName2nd);
                lambda[iEqn]      = lambdaInit1[iEqn];
			}
			if(mpi_rank==0) {
				printf(" - lambda1=");
				for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
					printf("%.3e ", lambdaInit1[iEqn]);
				cout<<endl;
			}

			arclength = calcArclength(normSqRatioInArclength, q0, q1, lambdaInit0, lambdaInit1, NcontrolEqns, nScal, weightLambda);

			break;
		default : // calculate q1 by using the Newton's method
			if(mpi_rank==0) cout<<" Calculate q1 by using the Newton method";
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
				sprintf(paramName2nd, "CONTROL_PARAM%d", iEqn);
				lambdaInit1[iEqn] = getParam("LAMBDA_INITIAL_1") -> getDouble(paramName2nd);
				lambda[iEqn] = lambdaInit1[iEqn];
			}
			if(mpi_rank==0) {
				printf(" - lambda1=");
				for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
					printf("%.3e ", lambdaInit1[iEqn]);
				cout<<endl;
			}

			if (!checkParam("NEWTON_PARAMETERS_2NDPT")) {
				ParamMap::add("NEWTON_PARAMETERS_2NDPT  MAX_ITER=100  ABS_RESID=1.0e-14  REL_RESID=1.0e-8"); // add default values
				if (mpi_rank == 0)
					cout<< "WARNING: added keyword \"NEWTON_PARAMETERS_2NDPT  MAX_ITER=100  ABS_RESID=1.0e-14  REL_RESID=1.0e-8\""<< " to parameter map!" << endl;
			}
			int maxIterNewton   = getParam("NEWTON_PARAMETERS_2NDPT")->getInt("MAX_ITER");
			double absTolNewton = getParam("NEWTON_PARAMETERS_2NDPT")->getDouble("ABS_RESID");
			double relTolNewton = getParam("NEWTON_PARAMETERS_2NDPT")->getDouble("REL_RESID");

			// Initialize q1 with an initial guess
			int initGuessNewton = getIntParam("INIT_GUESS_NEWTON", "0"); // Note: 0 = Previous rho, rhou, etc.
			                                                             //       1 = binary file by JOE (filename is defined by Q1_FILENAME)
			                                                             //       2 = binary file by IKE (filename is defined by Q1_FILENAME)
			if(initGuessNewton > 0) {
				filename = getStringParam("Q1_FILENAME", "Q1_PT00001.bin");

				if(initGuessNewton == 1) {
					if(mpi_rank==0) cout<<"   Initial guess is provided from a binary file dumped by JOE - filename = "<<filename<<endl;

					readJOEdumpedDataSerial(filename, q1);  // Note: writeJOEDataParallel() is defined in Joe
				} else if(initGuessNewton == 2) {
					if(mpi_rank==0) cout<<"   Initial guess is provided from a binary file dumped by IKE - filename = "<<filename<<endl;

					lambda0Temp2nd = new double [NcontrolEqns];
					lambda1Temp2nd = new double [NcontrolEqns];
					double arclengthTemp2nd;

					readPsALCdumpedDataSerial(filename, q1, step, lambda0Temp2nd, lambda1Temp2nd, NcontrolEqns, arclengthTemp2nd);

					if(mpi_rank==0) {
						for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
							printf("  lambda%d=%.3e", iEqn, lambda1Temp2nd[iEqn]);
						printf("  arclength=%.4e \n", arclengthTemp2nd);
					}

					// Check possible errors
					errorAbs2nd = 0.0;
					for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
						errorAbs2nd += fabs((lambda1Temp2nd[iEqn] - lambdaInit1[iEqn])/lambdaInit1[iEqn]);
					if(errorAbs2nd > LAMBDA_COMPATIBILITY_EPS) {
						if(mpi_rank==0) {
							printf("\n  CAUTION! The binary file does not match with the Ike.in file \n");
							printf("       Binary file:");
							for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
								printf(" lambda%d=%.5e ", iEqn, lambda1Temp2nd[iEqn]);
							printf("\n       Input file : ");
							for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
								printf(" lambda%d=%.5e ", iEqn, lambdaInit1[iEqn]);
							cout<<endl;
						}
					}

					// Free memory
					delete [] lambda0Temp2nd;
					delete [] lambda1Temp2nd;
				}

				// update the flow field
				updatePrimVarsFrom1DVec(q1, nScal);

				// Update ghost cells in the 2nd-layer
			    updateCvDataG1G2(rho, REPLACE_DATA);
			    updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
			    updateCvDataG1G2(rhoE, REPLACE_DATA);
			    for (int iScal = 0; iScal < nScal; iScal++)
			      updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);

				calcStateVariables();
//				calcMaterialProperties();
				setBC();
				for(int ifa=nfa; ifa<nfa_b2; ++ifa)
					setBC1D(ifa, UgpWithCvCompFlow::rho, UgpWithCvCompFlow::rhou, UgpWithCvCompFlow::rhoE); // Note: ncvggf's cannot be not updated, so the user should take care of it
				calcRansTurbViscMuet();
				MPI_Barrier(mpi_comm);
			} else {
				for(int icv=0; icv<ncv; ++icv) {
					int tempInt = icv*m;
					q1[tempInt] = rho[icv];
					for(int i=0; i<3; ++i)
						q1[tempInt+1+i] = rhou[icv][i];
					q1[tempInt+4] = rhoE[icv];

					if(nScal > 0) {
						for(int iScal=0; iScal<nScal; ++iScal) {
							double *scalArray = scalarTranspEqVector[iScal].phi;
							q1[tempInt+5+iScal] = scalArray[icv];
						}
					}
				}
			}
			if(mpi_rank==mpi_size-1)
				for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
					q1[ncv*(5+nScal)+iEqn] = lambda[iEqn];
			MPI_Barrier(mpi_comm);

			// Update q1 by Newton's method
			assert(rhs!=NULL);
			getSteadySolnByNewton(q1, rhs, maxIterNewton, absTolNewton, relTolNewton);

			delete petscSolver2; 	petscSolver2 = NULL;

			arclength = calcArclength(normSqRatioInArclength, q0, q1, lambdaInit0, lambdaInit1, NcontrolEqns, nScal, weightLambda);

			break;
	}
	// Write a summary of the Bifurcation code on a file
	if(howToInitSecond!=0 && mpi_rank==0) { // If Newton's method, it will be handled in the getSteadySolnByNewton() method
		FILE *fp = fopen(BIFUR_SUMMARY_FILENAME, "a");
		fprintf(fp, "%4d,            ,            ,            ,         , %11.4e, %.15f\n", step, lambdaInit1[0]-lambdaInit0[0], lambdaInit1[0]);
		fclose(fp);
	}

	assert(arclength > 0.0);

	temporalHook();
	rewriteQoI = false;
	writeQoIOnFile(step, QoIfilename, NcontrolEqns, rewriteQoI);
	string second_field_data = getStringParam("SECOND_FIELD_DATA", "NO");
	if(second_field_data == "YES" || second_field_data == "yes" || second_field_data == "Y" || second_field_data == "y")
		writeData(step);

	if (write_restart > 0) {
		string second_field_restart = getStringParam("SECOND_FIELD_RESTART", "NO");
		if(second_field_restart == "YES" || second_field_restart == "yes" || second_field_restart == "Y" || second_field_restart == "y")
			writeRestart(step); // Write the current flow field as a restart file

		string second_field_binary  = getStringParam("SECOND_FIELD_BINARY",  "NO");
		if(second_field_binary  == "YES" || second_field_binary  == "yes" || second_field_binary  == "Y" || second_field_binary  == "y") {
			char filename[40];
			sprintf(filename, "Q1_PT%05d.bin", step);
			string filenameQ1JOE = getStringParam("Q1_FILENAME", ""); // The filename used in the case that howToInitSecond == 3
			if(filenameQ1JOE.compare(filename)==0)
				sprintf(filename, "Q1_PT%05d.IKE.bin", step);
			writePsALCdumpedDataParallel(filename, step, arclength, lambdaInit0, lambdaInit1, NcontrolEqns, q1);
		}
	}

	++step;

	MPI_Barrier(mpi_comm);
}

/*
 * Method: singlePsALCstep()
 * -------------------------
 * Perform single continuation on the solution (or bifurcation) curve
 * Arguments: qVec = At the initial step of this function, the guessed solution is stored in this array.
 *                   It will be updated during the Newton process, and eventually stores the solution of this system.
 *                   At mpi_rank==mpi_size-1, the size of qVec must be [ncv*nVars+NcontrolEqns]
 *                   Otherwise, the size must be [ncv*nVars]
 *            q0 = The solution vector at two steps before.
 *                 At mpi_rank==mpi_size-1, the size of q0 must be [ncv*nVars+NcontrolEqns].
 *                 However, don't try to use the last NcontrolEqns elements of the array: the code might be modified 
 *            q1 = The solution vector at the previous step
 *                 At mpi_rank==mpi_size-1, the size of q1 must be [ncv*nVars+NcontrolEqns].
 *                 However, don't try to use the last NcontrolEqns elements of the array: the code might be modified 
 *            rhs =
 *            q_tangent =
 *            lambda_tangent =
 *            ds = arc-length: this is re-calculated in this method (arcLength) for error-checking, but ds is used in the calculations
 *            lambda0 =
 *            lambda1 =
 *            weightLambda =
 *            writeJacOnFile =
 *            factorReduceArclength = The user can reduce the length of the tangential vectors by setting this argument less than 1.0
 *                                    (The arclength itself must have been reduced in the function calling this method before calling this.
 *                                     However, the tangential vectors calculated in this method can be only reduced here)
 */
void IkeWithPsALC_AD::singlePsALCstep(double *qVec, const double *q0, const double *q1, double *rhs, double **q_tangent, double *lambda_tangent,
		const double arcLength, const double* lambda0, const double* lambda1, const double weightLambda, const int writeJacOnFile,
		const double factorReduceTangential /* = 1.0 */) {
	//const double lambda0, const double lambda1, double *dRHSdlambda, const double weightLambda) {
	/*
	 * Note: We want to solve two equations for the pseudo-arclength continuation:
	 *       1. The RANS system with a parameter lambda (lambda can be heat addition, etc.)
	 *          i.e., RHS(q, lambda) = 0   where q is the solution vector (which contains rho, rhou, rhoE, and scalars).
	 *
	 *       2. The tangential property of the solution curve:
	 *            We have the solution curve on the lambda - q plane and want to march a small length (ds) from the initial solution on the curve.
	 *            Let c = [q_t; lambda_t] is the tangential vector of the curve on the current solution.
	 *
	 *            Then, the tangential property is
	 *               c^{T} * [dq; dlambda] = ds
	 *
	 *       One way to solve the solution which satisfies above two equations is the Newton's method.
	 *       The standard Newton's method is the following:
	 *         If you want to find the solution of "f(x) = 0" iteratively,
	 *           then guess the solution (say, x_{0}) and iterate x_{n+1} = x_{n} - J^{-1}*f(x_{n}),  where J is the Jacobian matrix of f with respect to x
	 *                                                    until norm(f(x_{n+1})) becomes below some tolerance,
	 *
	 *       Let Q = [q; lambda]. You want to find Q which satisfies above two equations.
	 *       After some algebra, the you will end up with the Newton's iteration, Q_{n+1} = Q_{n} - L^{-1}*[RHS(q, lambda);  c^{T}*[dq; dlambda]-ds],
	 *       Where L = [ J       d(RHS)/d(lambda)
	 *                  q_t^{T}      lambda_t    ],      J = the Jacobian matrix of RHS w.r.t. q
	 *
	 * Note: residual = absTolNewton + relTolNewton * norm(rhs_initial)
	 */
	static bool firstCall = true;

	if(debugLevel>1 && mpi_rank==0) cout<<"IkeWithPsALC_AD::singlePsALCstep()"<<endl;

	assert(lambda0!=NULL && lambda1!=NULL);
	assert(jacMatrix.empty() && jacMatrixSTL.empty());
	assert(NcontrolEqns > 0);

	int nVars = nScal+5; // number of variables (rho, rhou, rhov, rhow, rhoE, and scalars)

	/***********************************************/
	/* Tangential conditions                       *
	 *   Q_tangent = (Q_{n} - Q_{n-1}) / arclength *
	 *                                             *
	 *         where  Q = [q                       *
	 *                     lambda ]                *
	 ***********************************************/
	// Calculate the tangent of the bifurcation curve: Note: (q_tangent, lambda_tangent) must be an unit vector, i.e., ||(q_tangent,lambda_tangent)||_W = 1.0
	for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
		for(int i=0; i<ncv*nVars; ++i)
			q_tangent[iEqn][i] = factorReduceTangential*(q1[i] - q0[i]);
	        // Since ds was already reduced by factorReduceTangential, this will make [q_tangent; lambda_tangent] a unit vector
	        // after calling calcUnitVecWithWeightedNorm() several lines below:
	        //   factor*(q1-q0)   factor*(q1-q0)    q1-q0
	        //   -------------- = -------------- =  ------
	        //       ds_new       factor*ds_old     ds_old

	if(mpi_rank==mpi_size-1)
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
			lambda_tangent[iEqn] = factorReduceTangential*(lambda1[iEqn] - lambda0[iEqn]);

	// Normalize q_tangent and lambda_tangent !!!!!!!!
	calcUnitVecWithWeightedNorm(q_tangent[0], lambda_tangent, q_tangent[0], lambda_tangent, arcLength, nVars, NcontrolEqns);

	MPI_Bcast(lambda_tangent, NcontrolEqns, MPI_DOUBLE, mpi_size-1, mpi_comm);
		// Note: lambda_tangent will be used only at mpi_rank==mpi_size-1, but broadcasting it just in case

	if(debugLevel>0 && firstCall) {
		showOverviewOnTangentVectors(q_tangent[0], lambda_tangent, nVars, NcontrolEqns);
	}

	if(debugLevel>1) {
		double totalLength = vecDotVecWithWeight(q_tangent[0], q_tangent[0], lambda_tangent, lambda_tangent, weightLambda, nVars, NcontrolEqns);
		if(fabs(totalLength-1.0) > 1.0e3 * sqrt(cvora[mpi_size]*7.0) * MACHINE_EPS) {
			if(mpi_rank==0)
				cout<<"WARNING!! singlePsALCstep(): the tangent vector (=[q_tangent, lambda_tangent]^T) is not an unit vector -- length = "<<totalLength<<endl;
		}
	}

	/***********************************************/
	/* Predictor Step                              *
	 *   Q = Q_{n} + arclength*Q_tangent           *
	 *                                             *
	 *         where  Q = [q                       *
	 *                     lambda ]                *
	 ***********************************************/
	// initial guess of the solution ( note: qVec = [q_guess; lambda_guess] )
	if(NcontrolEqns>1) {
		if(mpi_rank==0)
			cerr<<"ERROR in IkeWithPsALC_AD::singlePsALCstep(): N_PSALC_CONTROL_PARAMS="<<NcontrolEqns<<" > 1 -- This it not supported by now"<<endl;
		throw(PSALC_ERROR_CODE);
	}

	if(step-startingStep==2 && initThirdFromQ1) { // From a binary file
		if(mpi_rank==0) cout<<"> Initial guess from a binary file dumped by IKE";
		string filename = getStringParam("Q2_FILENAME", "Q1_PT00002.dumped.00100.bin");
		if(mpi_rank==0) cout<<" - filename = "<<filename<<endl;

		double *lambda0Temp3rd = new double [NcontrolEqns];
		double *lambda1Temp3rd = new double [NcontrolEqns];
		double arclengthTemp3rd;
		int stepTemp;

		readPsALCdumpedDataSerial(filename, qVec, stepTemp, lambda0Temp3rd, lambda1Temp3rd, NcontrolEqns, arclengthTemp3rd);

		if(mpi_rank == mpi_size-1) {
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)       // Even though only NcontrolEqns==1 is supported, the code has a general form for future uses
				qVec[ncv*nVars+iEqn] = lambda1Temp3rd[iEqn];
		}

		delete [] lambda0Temp3rd;
		delete [] lambda1Temp3rd;

		MPI_Barrier(mpi_comm);
	} else { // From extrapolation
		if(mpi_rank == mpi_size-1) {
			for(int i=0; i<ncv*nVars; ++i)
				qVec[i] = q1[i] + arcLength*q_tangent[0][i]; // TO DO: Only NcontrolEqns==1 is supported by now: i.e. q_tangent[0]
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)       // Even though only NcontrolEqns==1 is supported, the code has a general form for future uses
				qVec[ncv*nVars+iEqn] = lambda1[iEqn] + arcLength*lambda_tangent[iEqn];
		} else {
			for(int i=0; i<ncv*nVars; ++i)
				qVec[i] = q1[i] + arcLength*q_tangent[0][i]; // TO DO: Only NcontrolEqns==1 is supported by now: i.e. q_tangent[0]
		}
	}

	// update the solution
	for(int icv=0; icv<ncv; ++icv) {
		rho[icv] = qVec[icv*nVars];
		for(int i=0; i<3; ++i)
			rhou[icv][i] = qVec[icv*nVars+1+i];
		rhoE[icv] = qVec[icv*nVars+4];
		for(int iScal=0; iScal<nScal; ++iScal) {
			double *phi = scalarTranspEqVector[iScal].phi;
			phi[icv] = qVec[icv*nVars+5+iScal];
		}
	}
	updateCvDataG1G2(rho, REPLACE_DATA);
	updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
	updateCvDataG1G2(rhoE, REPLACE_DATA);
	for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
		updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);

	calcStateVariables();     // Compute velocity, pressure, temperature, enthalpy and speed of sound at cell centers
//	calcMaterialProperties(); // Compute viscosity and thermal diffusivity at cell faces
	setBC();
	for(int ifa=nfa; ifa<nfa_b2; ++ifa)
		setBC1D(ifa, UgpWithCvCompFlow::rho, UgpWithCvCompFlow::rhou, UgpWithCvCompFlow::rhoE); // Note: ncvggf's cannot be not updated, so the user should take care of it
	calcRansTurbViscMuet();

	if(mpi_rank == mpi_size-1) {
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
			lambda[iEqn] = qVec[ncv*nVars+iEqn];
	}
	MPI_Bcast(lambda, NcontrolEqns, MPI_DOUBLE, mpi_size-1, mpi_comm);

	if(debugLevel > 1 && mpi_rank == 0) {
		printf("> IkeWithPsALC_AD::singlePsALCstep(): NcontrolEqns = %d \n", NcontrolEqns);
	}
	if(debugLevel > 0 && mpi_rank == 0) {
		printf("> IkeWithPsALC_AD::singlePsALCstep(): initially guessed lambda = ");
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
			printf("%.5e  ", lambda[iEqn]);
		printf("\n");
	}

	/***********************************************/
	/* Corrector step                              *
	 ***********************************************/
	string correctorMethod = getStringParam("CORRECTOR_METHOD", "NEWTON");

	/*-----------------------------------------------*/
	/* Corrector step (Newton iteration)             *
	 *   Solve Q_{n+1} = Q_{n} - L^{-1}*RESID        *
	 *                                               *
	 *   Where L = [      J          dF/dlambda      *
	 *               -q_tangent'  -lambda_tangent ]  *
	 *                                               *
	 *         RESID = [ RHS(q, lambda)              *
	 *                   Nres           ]            *
	 *         c.f.) Nres = ds-c^{T}*[dq; dlambda]   *
	 *                                               *
	 *         Q = [q                                *
	 *              lambda ]                         *
	 *-----------------------------------------------*/
	if(correctorMethod.compare("NEWTON") == 0) {
		assert(rhs!=NULL);

		if (!checkParam("NEWTON_PARAMETERS_PSALC")) {
			ParamMap::add("NEWTON_PARAMETERS_PSALC  MAX_ITER=100  ABS_RESID=1.0e-14  REL_RESID=1.0e-8"); // add default values
			if (mpi_rank == 0)
				cout<< "WARNING: added keyword \"NEWTON_PARAMETERS_PSALC  MAX_ITER=100  ABS_RESID=1.0e-14  REL_RESID=1.0e-8\""<< " to parameter map!" << endl;
		}
		int maxIterNewton   = getParam("NEWTON_PARAMETERS_PSALC")->getInt("MAX_ITER");
		double absTolNewton = getParam("NEWTON_PARAMETERS_PSALC")->getDouble("ABS_RESID");
		double relTolNewton = getParam("NEWTON_PARAMETERS_PSALC")->getDouble("REL_RESID");

		getSteadySolnByNewton(qVec, rhs, maxIterNewton, absTolNewton, relTolNewton,
				writeJacOnFile, NcontrolEqns, q1, q_tangent, lambda_tangent, lambda0, lambda1, weightLambda, arcLength);
		// Note: getSteadySolnByNewton() can access to the member variable NcontrolEqns, so you don't have to pass it to the method
	}

	/*----------------------------------------------------------------------------------------------------------------------*/
	/* Corrector step (Time stepping with implicit Euler)                                                                   *
	 *   Solve dQ/dt = f(Q)                                                                                                 *
	 *   The tangential equation is                                                                                         *
	 *           dLambda                                                                                                    *
	 *         V ------- = -( q_tangent^T*(u-u_guess) + lambda_tangent*(lambda - lambda_tangent) ) * sign(lambda_tangent),  *
	 *              dt                                                                                                      *
	 *         where V = 1 / weightTangentCond / |lambda_tangent|.                                                          *
	 *----------------------------------------------------------------------------------------------------------------------*/
	else if (correctorMethod.compare("TIME_STEPPING") == 0) {
		// Allocate for initial guess: they will be used for the tangential condition (they are used to implement the penalty method in the flow RHS)
		double* q_guess = new double [ncv*nVars];
		for(int i=0; i<ncv*nVars; ++i)
			q_guess[i] = qVec[i];
		double* lambda_guess = NULL;
		if(NcontrolEqns > 0) {
			lambda_guess = new double [NcontrolEqns];
			for(int iScal=0; iScal<NcontrolEqns; ++iScal)
				lambda_guess[iScal] = lambda[iScal];
		}

		// Solve the corrector step using time-stepping
		if (!checkParam("TIMESTEPPING_PARAMETERS_PSALC")) {
			ParamMap::add("TIMESTEPPING_PARAMETERS_PSALC  MAX_STEP=100  ABS_RESID=1.0e-14  REL_RESID=1.0e-8"); // add default values
			if (mpi_rank == 0)
				cout<< "WARNING: added keyword \"TIMESTEPPING_PARAMETERS_PSALC  MAX_STEP=100  ABS_RESID=1.0e-14  REL_RESID=1.0e-8\""<< " to parameter map!" << endl;
		}
		int maxTimeStep = getParam("TIMESTEPPING_PARAMETERS_PSALC")->getInt("MAX_STEP");
		double absTolTS = getParam("TIMESTEPPING_PARAMETERS_PSALC")->getDouble("ABS_RESID");
		double relTolTS = getParam("TIMESTEPPING_PARAMETERS_PSALC")->getDouble("REL_RESID");

		double startCFL = getDoubleParam("CFL", "1.0");

		getSteadySolnByTimeStepping(qVec, maxTimeStep, absTolTS, relTolTS, startCFL,
				q_guess, q_tangent, NcontrolEqns, lambda1, lambda_guess, lambda_tangent,
				weightLambda, arcLength);  // Note: lambda1 is required only to save a binary file.

		// Clear memory
		delete [] q_guess;
		if(lambda_guess != NULL) 	delete [] lambda_guess;

		// Boost convergence with Newton's method
		string booleanString = getStringParam("AFTER_TIMESTEPPING_LAUNCH_NT", "FALSE");
		std::transform(booleanString.begin(), booleanString.end(), booleanString.begin(), ::tolower);
		if(booleanString.compare("true")==0 || booleanString.compare("yes")==0) {
			assert(rhs!=NULL);

			int maxIterNewton   = getParam("NEWTON_PARAMETERS_PSALC")->getInt("MAX_ITER");
			double absTolNewton = getParam("NEWTON_PARAMETERS_PSALC")->getDouble("ABS_RESID");
			double relTolNewton = getParam("NEWTON_PARAMETERS_PSALC")->getDouble("REL_RESID");

			getSteadySolnByNewton(qVec, rhs, maxIterNewton, absTolNewton, relTolNewton,
					writeJacOnFile, NcontrolEqns, q1, q_tangent, lambda_tangent, lambda0, lambda1, weightLambda, arcLength);
			// Note: getSteadySolnByNewton() can access to the member variable NcontrolEqns, so you don't have to pass it to the method
		}
	}

	else {
		if(mpi_rank == 0)
			cout << "ERROR in IkeWithPsALC_AD::singlePsALCstep(): CORRECTOR_METHOD = " << correctorMethod << " is not supported" << endl;
		throw(PSALC_ERROR_CODE);
	}

	/***********************************************/
	/* Clean up                                    */
	/***********************************************/
	// Update normal Joe variables (e.g. press, vel, etc. & BC) for tecplot output
//	updateJoeWithModelsProperties();

	// Free memory
	if(!jacMatrix.empty())
		jacMatrix.clear();
	if(!jacMatrixSTL.empty())
		jacMatrixSTL.clear();

	MPI_Barrier(mpi_comm);

	firstCall = false;
}

/*
 * Method: arclengthControl
 * ------------------------
 * Calculate the ratio for the arclength to be reduced or increased
 * Note: Increasing the arclength is more conservative than reducing it.
 */
double IkeWithPsALC_AD::arclengthControl(const double residNormInit, const double residNormInit_opt, const int niterNewton, const int niterNewton_opt) {
	double upperBound = 2.0;
	double lowerBound = 0.5;

	double exponent0 = 0.138;
	double exponent1 = 0.5;

	// Ratio based on the residual norm of the initial guess
	double ratio0 = pow(residNormInit_opt / residNormInit, exponent0);
	// Ratio based on the number of Newton-iterations (Arie de Niet, Dept. of Mathemantics, University of Groningen, Mater's thesis, 2002)
	double ratio1 = pow(niterNewton_opt / niterNewton,     exponent1);

	double newRatio = 1.0;
	if(ratio0 >= 1.0)
		newRatio = min(max(max(ratio0, ratio1), 1.02), upperBound); // At least, 2% increase
	else
		newRatio = max(min(min(ratio0, ratio1), 0.98), lowerBound); // At least, 2% decrease

	// Force to reduce the arclength if too many Newton iterations were used!
	if(ratio1 < 1.0)
		newRatio = max(min(min(ratio0, ratio1), 0.98), lowerBound);

	return newRatio;
}

/*
 * Method: getSteadySolnByNewton
 * -----------------------------
 * Newton-Raphson method with barriers and backtracking
 * Note: tolerance = absTolNewton + relTolNewton * norm(rhs_initial)   -- C.T.Kelley, SIAM 2003, Chap.1.5
 *
 * In the pseudo-arclength method, additional system control parameters must be introduced (e.g. heat release rate).
 * If NcontrolEqns > 0, Newton-Raphson will also take care of it.
 * Caution: if NcontrolEqns > 0, the size of rhs at mpi_rank==mpi_size-1 should ncv*(5+nScal)+NcontrolEqns or ncv_gg*(5+nScal)+NcontrolEqns
 */
void IkeWithPsALC_AD::getSteadySolnByNewton(double *q, double* rhs, const int maxIterNewton, const double absTolNewton, const double relTolNewton,
		const bool writeJacOnFile /* = false */,
		const int NcontrolEqns /* =0 */, const double *q1 /* =NULL */, double** q_tangent /* =NULL */,
		double *lambda_tangent /* =NULL */, const double *lambda0 /* =NULL */, const double *lambda1 /* =NULL */, const double weightLambda /* =0.0 */,
        const double arcLength /* =0 */) {
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
	int maxIterLS = newtonParam.maxIterLS;    // Since PETSc GMRes tresholds can be changed, initialize the value here
	double zeroAbsLS = newtonParam.zeroAbsLS; // Since PETSc GMRes tresholds can be changed, initialize the value here
	double zeroRelLS = newtonParam.zeroRelLS; // Since PETSc GMRes tresholds can be changed, initialize the value here

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
					cout<<"ERROR in IkeWithPsALC_AD::getSteadySolnByNewton(): compatibility check by compatibilityCheck() FAILs"<<endl;
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

	// Scaling factor of the RHS: Currently, this will affect only the "1D-style"
//	if(weightRhsMethod == REF_VALUES_AND_DT_OVER_VOL) {
//		double dtOverVol_min, dtOverVol_max;
//		double cfl_target = getDoubleParam("RHS_SCALING_CFL_TARGET", "1.0");
//		calcDtOverVol(dtOverVol_min, dtOverVol_max, cfl_target);
//
//		if(debugLevel>0 && mpi_rank==0)
//			printf("> RHS will be also scaled to match CFL_TARGET=%.2f: Max=%.4e, Min=%.4e\n", cfl_target, dtOverVol_max, dtOverVol_min);
//	}

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
	if(NcontrolEqns>0) {
		calcNres(Nres, q, rhs, m, NcontrolEqns, q1, q_tangent, lambda_tangent, lambda1, weightLambda, arcLength);
		MPI_Barrier(mpi_comm);

		if(mpi_rank==0) {
			for(int iParam=0; iParam<NcontrolEqns; ++iParam) {
				if(fabs(Nres[iParam]/arcLength) > 1.0e4*double(NcontrolEqns*cvora[mpi_size])*MACHINE_EPS) { // Note: Nres must be very close to zero
					cout<<"WARNING in getSteadySolnByNewton(): iterNewton=="<<iterNewton<<", The tangential condition["<<iParam<<"] shows a big residual."<<endl
					    <<"                                    Tangential-residual = "<<Nres[iParam]<<" (target arclength="<<arcLength<<")"<<endl;
				}
				if(isNaN(Nres[iParam])) {
					cout<<"ERROR in getSteadySolnByNewton(): iterNewton=="<<iterNewton<<", The residual of the tangential condition becomes NaN for iParam="<<iParam<<endl;
					throw(PSALC_ERROR_CODE);
				}
			}
		}

		if(mpi_rank == mpi_size-1) {
			for(int iParam=0; iParam<NcontrolEqns; ++iParam)
				rhs[ncv*m+iParam] = Nres[iParam];
		}
	}

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
	if(NcontrolEqns > 0) { // If NcontrolEqns==0, q1 was not passed (i.e. q1 is NULL if NcontrolEqns==0)
        double mySumPhiSq = 0.0, sumPhiSq;

		for(int i=0; i<ncv*m; ++i)
				mySumPhiSq += pow(q[i]-q1[i], 2.0);

		if(mpi_rank==mpi_size-1 && NcontrolEqns>0)
			for(int i=0; i<NcontrolEqns; ++i)
				mySumPhiSq += weightLambda * pow(q[ncv*(5+nScal)+i]-q1[ncv*(5+nScal)+i], 2.0);

		MPI_Allreduce(&mySumPhiSq, &sumPhiSq, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

		deltaQnorm = sqrt(sumPhiSq);
	}

	if(newtonParam.trustRegionSize > 8.0*deltaQnorm && NcontrolEqns > 0) {
		if(mpi_rank==0)
			cout<<"WARNING in getSteadySolnByNewton(): TRUST_REGION_SIZE="<<newtonParam.trustRegionSize<<" is larger than 8 * DELTA_Q of the init. guess="<<deltaQnorm<<endl;
		newtonParam.trustRegionSize = 8.0*deltaQnorm;
		if(mpi_rank==0)
			cout<<"                                    Reduce TRUST_REGION_SIZE to "<<newtonParam.trustRegionSize<<endl
			    <<endl;
	}

	// Show the result on the screen
	if(isNaN(residNormTot)) {
		if(mpi_rank==0) cerr<<"getSteadySolnByNewton(): The residual of the initial field is NaN"<<endl;
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
				    <<"WARNING in getSteadySolnByNewton(): The initial residual is greater than INIT_THIRD_RESID="<<startingResidual<<endl
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

	int iterNewton_stab_newton = iterNewton; // This is used in the stabilized Newton's method to increase "alpha" in the diagonal.

	while(!done) {
		/******
		 ** Check iterations and report error
		 ******/
		if(iterNewton > maxIterNewton+startingNewtonIter) { // If the Newton procedure fails to converge (Note: startingNewtonIter=0 unless initial guess for the thrid point)
			if(mpi_rank==0) {
				cout<<"ERROR in getSteadySolnByNewton(): Newton iter has NOT been converged until the "<<iterNewton-1<<"th iter"<<": residual="<<residNormTot;
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
//			if (iterNewton >= newtonParam.startDecLS) {
				int rampCount = int(abs(iterNewton - newtonParam.startDecLS) / newtonParam.intervalDecLS);

				maxIterLS = min(newtonParam.maxIterLS + newtonParam.incIterLS * rampCount, newtonParam.maxFinalIterLS);
				zeroAbsLS = max(newtonParam.zeroAbsLS * pow(newtonParam.decZeroLS, double(rampCount)), newtonParam.minZeroAbsLS);
				zeroRelLS = max(newtonParam.zeroRelLS * pow(newtonParam.decZeroLS, double(rampCount)), newtonParam.minZeroRelLS);

				petscSolver2->setTresholds(zeroAbsLS, zeroRelLS, maxIterLS);

				if(debugLevel>1 && mpi_rank==0) {
					printf("RAMP_LSNT_PARAMETERS: MAX_ITER=%d  MIN_ABS_RESID=%.2e  MIN_REL_RESID=%.2e\n", maxIterLS, zeroAbsLS, zeroRelLS);
				}
			}

		/******
		 ** Get q_tangent_inJacMat and lambda_tangent_inJacMat that will be actually used in the Jacobian matrix
		 ******/
		if(NcontrolEqns > 0) {
			assert(q_tangent_inJacMat != NULL && q_tangent != NULL && lambda_tangent_inJacMat != NULL && lambda_tangent != NULL);

			if(howToCalcJac == ROW_1D) {
				for(int iParam=0; iParam<NcontrolEqns; ++iParam)
					for(int i=0; i<ncv*m; ++i)
//						q_tangent_inJacMat[iParam][i] = -q_tangent[iParam][i];
						q_tangent_inJacMat[iParam][i] = -q_tangent[iParam][i] * AreaLambda;
				for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
//					lambda_tangent_inJacMat[iEqn] = -lambda_tangent[iEqn];
					lambda_tangent_inJacMat[iEqn] = -lambda_tangent[iEqn] * AreaLambda;
			} else if(howToCalcJac == ORDINARY_2D) { // Note: ORDINARY_2D is still using old formulation unlike ROW_1D
				for(int iParam=0; iParam<NcontrolEqns; ++iParam)
					for(int i=0; i<ncv*m; ++i)
						q_tangent_inJacMat[iParam][i] = q_tangent[iParam][i];
				for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
					lambda_tangent_inJacMat[iEqn] = lambda_tangent[iEqn];
			} else {
				throw(PSALC_ERROR_CODE);
			}
		}

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

		if(howToCalcJac == ROW_1D) {
			if(firstCall && iterNewton==startingNewtonIter+1 && debugLevel>0 && mpi_rank==0)
				wtime0 = MPI_Wtime();

			// Note: If you want to print out the matrix stored in PETSc on the screen, set the ShowPetscMatrixMatlab variable as true
			//       ShowPetscMatrixMatlab is defined as a static variable in the PetscSolver2.h file
			//       e.g. if(iterNewton==1 && NcontrolEqns>0) ShowPetscMatrixMatlab = true;
			bool useOldSolnAsInitGuess = false; // If true, use the previous solution vector as the initial guess.
			                                    // Otherwise, apply PC(pre-conditioner) to get the inital guess (the Knoll trick) -- For Bifurcation, this works better.

			solveLinSysNSCoupled2<MatComprsedSTL>(phi, jacMatrixSTL, rhs, nIter, absResid, kspMonitorHistory,
					useOldSolnAsInitGuess, zeroAbsLS, zeroRelLS, maxIterLS, nScal, monitorConvergInterval,
					NcontrolEqns, ncv_gg, q_tangent_inJacMat, lambda_tangent_inJacMat, step, iterNewton);
			// Note that you should pass nScal to this function instead of m(=5+nScal)

			if(firstCall && iterNewton==startingNewtonIter+1 && debugLevel>0 && mpi_rank==0)
				wtimeLS = MPI_Wtime();
			if(firstCall && iterNewton==startingNewtonIter+1 && debugLevel>0 && mpi_rank==0)
				cout<<"   > getSteadySolnByNewton(): Runtime for the LS solver   [sec] = "<<wtimeLS - wtime0<<endl;
		} else if(howToCalcJac == ORDINARY_2D) {
			// Note: If you want to print out the matrix stored in PETSc on the screen, set the ShowPetscMatrixMatlab variable as true
			//       ShowPetscMatrixMatlab is defined as a static variable in the PetscSolver2.h file
			//       e.g. if(iterNewton==1 && NcontrolEqns>0) ShowPetscMatrixMatlab = true;
			bool useOldSolnAsInitGuess = false; // If true, use the previous solution vector as the initial guess.
			                                    // Otherwise, apply PC(pre-conditioner) to get the inital guess (the Knoll trick) -- For Bifurcation, this works better.
			int saveConvergInterval = int(max(monitorConvergInterval/10, 1.01));

			solveLinSysNSCoupled2<MatComprsed>(phi, jacMatrix, rhs, nIter, absResid, kspMonitorHistory,
					useOldSolnAsInitGuess, zeroAbsLS, zeroRelLS, maxIterLS, nScal, monitorConvergInterval,
					NcontrolEqns, ncv_gg, q_tangent_inJacMat, lambda_tangent_inJacMat, step, iterNewton);
			// Note that you should pass nScal to this function instead of m(=5+nScal)
		} else {
			throw(PSALC_ERROR_CODE);
		}

		if(nIter == maxIterLS  &&  absResid > zeroAbsLS) {
			if(debugLevel>0 || (newtonParam.incIterLS!=0 && newtonParam.incIterLS!=1.0)) {
				// Show the warning message on the screen
				double trueAbsResid = petscSolver2->calcTrueResidualNorm(); // Since "absResid" contains the last approximate & left-preconditioned residual norm,
				                                                            // we want to calculate the true residual here.

				if(mpi_rank==0)
					cout<<"           >> WARNING in getSteadySolnByNewton(): Linear solver reaches MAX ITER (iterNewton="<<iterNewton<<"): nIter="<<nIter<<", absResid="<<trueAbsResid<<endl; // absResid<<endl;

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
        		cout<<"WARNING in getSteadySolnByNewton(): After solving the linear system,"<<endl
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
                printf("           >> WARNING! getSteadySolnByNewton(): Checking Negative rho, press, or kine - Negative occurs at total %d CVs / %d FAs\n", negativeValCount_CV, negativeValCount_FA);
                printf("                                                                                        Reduce RELAXATION to %.3e\n", relaxation);
            }
            
            updateFlow_primVars(q, phi, -relaxation, nScal); // note: You should not update q here! (because of backtracking)
                                                             //       Only rho, rhou, rhoE, and scalars should be updated here
            
            if(NcontrolEqns > 0)
				updateLambda(q, phi, -relaxation, NcontrolEqns, nScal);
        }
        MPI_Barrier(mpi_comm);

        /******
         ** Check if Delta(Lambda) is too large
         ******/
        double relaxBeforeDlambda = relaxation;

        // If delta_lambda is not that big, then the backtracking algorithm can take care of it.
        // However, if delta_lambda is extremely large, we need to do something here.
		if(NcontrolEqns>0) {
			if(mpi_rank==mpi_size-1) {
                for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
                    double delta_lambda      = relaxBeforeDlambda * phi[ncv*(5+nScal)+iEqn];
                    double delta_lambda_prev = lambda1[iEqn] - lambda0[iEqn];
                    if(fabs(delta_lambda) > 7.77*fabs(delta_lambda_prev)) {
                        cout<<"WARNING in getSteadySolnByNewton(): delta_lambda["<<iEqn<<"]="<< delta_lambda
                        	<<" is greater than 777% of lambda["<<iEqn<<"]-lambda1["<<iEqn<<"]="<<delta_lambda_prev<<")"<<endl;

                        if(fabs(delta_lambda_prev)>1.0e-10) {
                        	cout<<"                                    reduce relaxation from "<<relaxation;
                            relaxation = min(max(RELAXATION_EPS, 7.77/fabs(delta_lambda/delta_lambda_prev)*relaxation), relaxBeforeDlambda);
                            cout<<" to "<<relaxation<<endl;
                        }
                    }
                }
			}
		}
        MPI_Bcast(&relaxation, 1, MPI_DOUBLE, mpi_size-1, mpi_comm);

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
//				if(fabs(Nres[iParam]/arcLength) > 1.0e4*double(NcontrolEqns*cvora[mpi_size])*MACHINE_EPS) // Note: Nres must be very close to zero
//					if(mpi_rank==0) {
//						cout<<"WARNING in getSteadySolnByNewton(): iterNewton=="<<iterNewton<<", The tangential condition["<<iParam<<"] was NOT satisfied."<<endl
//						    <<"                                    Tangetial-residual = "<<Nres[iParam]<<" (target arclength="<<arcLength<<")"<<endl;
//					}

				if(isNaN(Nres[iParam])) {
					if(mpi_rank==0)
						cerr<<"ERROR in getSteadySolnByNewton(): iterNewton=="<<iterNewton<<", The residual of the tangential condition becomes NaN for iParam="<<iParam<<endl;
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
				if(mpi_rank==0) cout<<"ERROR! getSteadySolnByNewton(): After negative vars = "<<countNegative<<endl;
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
		if(NcontrolEqns>0) {
			calcNres(Nres, q, rhs, m, NcontrolEqns, q1, q_tangent, lambda_tangent, lambda1, weightLambda, arcLength);

			for(int iParam=0; iParam<NcontrolEqns; ++iParam)
				if(fabs(Nres[iParam]/arcLength) > 1.0e-10 * sqrt(double(NcontrolEqns*cvora[mpi_size]))) {
					if(mpi_rank==0) {
						double arcLengthNew = arcLength - Nres[iParam];
						double lambdaPart = weightLambda * lambda_tangent[iParam] * (lambda[iParam] - lambda1[iParam]);
						cout<<"WARNING in getSteadySolnByNewton(): The tangential condition is NOT satisfied: Nres["<<iParam<<"] = "<<Nres[iParam]<<" (target arclength="<<arcLength<<")"<<endl
							<<"                                    (lambdaPart + flowPart = "<<lambdaPart<<" + "<<arcLengthNew-lambdaPart<<" = "<<arcLengthNew<<")"<<endl;
					}
				}

			if(mpi_rank == mpi_size-1) {
				for(int iParam=0; iParam<NcontrolEqns; ++iParam)
					rhs[ncv*m+iParam] = Nres[iParam];
			}
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

////IKJ
//for(int iParam=0; iParam<NcontrolEqns; ++iParam)
//for(int icv=0; icv<ncv; ++icv) {
//	int indexStart = icv*m;
//	residField[icv]  = fabs(q_tangent[iParam][indexStart]);
//	residField2[icv] = fabs(q_tangent[iParam][indexStart+1]);
//	residField3[icv] = fabs(q_tangent[iParam][indexStart+2]);
//	residField4[icv] = fabs(q_tangent[iParam][indexStart+4]);
//	residField5[icv] = fabs(q_tangent[iParam][indexStart+5]);
//	residField6[icv] = fabs(q_tangent[iParam][indexStart+6]);
//}
//updateCvDataG1G2(residField,  REPLACE_DATA);
//updateCvDataG1G2(residField2, REPLACE_DATA);
//updateCvDataG1G2(residField3, REPLACE_DATA);
//updateCvDataG1G2(residField4, REPLACE_DATA);
//updateCvDataG1G2(residField5, REPLACE_DATA);
//updateCvDataG1G2(residField6, REPLACE_DATA);

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
				printf("     %15.8e\n", lambda[0]);
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
			// Note: For now, GENERAL_NEWTON_HU and GENERAL_NEWTON_HUSEO are not affected by this ramping
			// Check ramping
			int startIncAlpha    = newtonParam.stabilizationAlphaRamp_startInc;
			int intervalIncAlpha = newtonParam.stabilizationAlphaRamp_interval;
			double incAlpha    = newtonParam.stabilizationAlphaRamp_incAlpha;
			double targetAlpha = newtonParam.stabilizationAlphaRamp_targetAlpha;

//			if ((iterNewton > startIncAlpha) && (iterNewton%intervalIncAlpha == 0)) {
			if ((iterNewton_stab_newton > startIncAlpha) && (iterNewton_stab_newton%intervalIncAlpha == 0)) {
				if(incAlpha>1.0 && alpha_stab_newton<targetAlpha) {
					if(relaxation > incAlpha-1.0) {
//						int exponent = int(abs(iterNewton - startIncAlpha) / std::max<int>(intervalIncAlpha, 1));
						int exponent = int(abs(iterNewton_stab_newton - startIncAlpha) / std::max<int>(intervalIncAlpha, 1));
						alpha_stab_newton = min(newtonParam.stabilizationAlpha * pow(incAlpha, double(exponent)), targetAlpha);
						++iterNewton_stab_newton;
					} else {
						if(mpi_rank==0) printf("           >> Note: Doesn't increase alpha_stab_newton because relaxation(=%g) <= FACTOR_INC-1.0(=%g)\n", relaxation, incAlpha-1.0);
					}
				}
				if(incAlpha<1.0 && alpha_stab_newton>targetAlpha) {
					if(relaxation > incAlpha) {
//						int exponent = int(abs(iterNewton - startIncAlpha) / std::max<int>(intervalIncAlpha, 1));
						int exponent = int(abs(iterNewton_stab_newton - startIncAlpha) / std::max<int>(intervalIncAlpha, 1));
						alpha_stab_newton = max(newtonParam.stabilizationAlpha * pow(incAlpha, double(exponent)), targetAlpha);
						++iterNewton_stab_newton;
					} else {
						if(mpi_rank==0) printf("           >> Note: Doesn't decrease alpha_stab_newton because relaxation(=%g) <= FACTOR_INC(=%g)\n", relaxation, incAlpha);
					}
				}
			} else {
				++iterNewton_stab_newton;
			}
		}
////IKJ
//writeData(step, iterNewton);
	}

	if(debugLevel>0 && mpi_rank == 0)
		printf("Newton solver converged after %5dth outer iter: residual = %12.5e \n", iterNewton-1, residNormTot);

	// finalHookNewton()
	finalHookNewton();

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
					cout<<"ERROR! IkeWithPsALC_AD::getSteadySolnByNewton(): Both jacMatrixSTL and jacMatrix are empty!"<<endl;
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
 * Method: freeMemGetSteadySolnByNewton
 * ---------------------------------------
 * Clear memory allocated for the getSteadySolnByNewton() function
 * Note: The arrays will NOT be re-initialized by NULL
 */
void IkeWithPsALC_AD::freeMemGetSteadySolnByNewton(double *phi, double* residNorm, double* residNormOld, double *Nres,
		MatComprsed& jacMatrix, MatComprsedSTL& jacMatrixSTL) {
	delete [] phi;
	delete [] residNorm;
	delete [] residNormOld;
	delete [] Nres;

	if(!jacMatrix.empty())
		jacMatrix.clear();
	if(!jacMatrixSTL.empty())
		jacMatrixSTL.clear();

	MPI_Barrier(mpi_comm);
}

void IkeWithPsALC_AD::freeMemGetSteadySolnByNewton(double *phi, double* residNorm, double* residNormOld, double *Nres,
		double *turbScalWeightVec, double** weighted_q_tangent, double* weighted_lambda_tangent, const int nVars, const int NcontrolEqns,
		MatComprsed& jacMatrix, MatComprsedSTL& jacMatrixSTL) {
	delete [] phi;
	delete [] residNorm;
	delete [] residNormOld;
	delete [] Nres;

	if(turbScalWeightVec != NULL)
		delete [] turbScalWeightVec;
	if(weighted_lambda_tangent != NULL)
		delete [] weighted_lambda_tangent;
	if(weighted_q_tangent      != NULL)
		freeMem2D(weighted_q_tangent, 0, NcontrolEqns-1, 0, (ncv*nVars)-1);

	if(!jacMatrix.empty())
		jacMatrix.clear();
	if(!jacMatrixSTL.empty())
		jacMatrixSTL.clear();

	MPI_Barrier(mpi_comm);
}

/*
 * Method: getSteadySolnByTimeStepping
 * -----------------------------------
 * Time stepping method (for example, implicit Euler) with a penalty term.
 * Barrier and trust-region are also used to stabilize the calculation.
 * Note: tolerance = absTolTimeStepping + relTolTimeStepping * norm(rhs_initial)   -- C.T.Kelley, SIAM 2003, Chap.1.5
 *
 * In the pseudo-arclength method, additional system control parameters must be introduced (e.g. heat release rate).
 * If NcontrolEqns > 0, this function will also take care of it.
 * Caution: if NcontrolEqns > 0, you need one more equation for the system control paramter.
 * Note that lambda1 is required only to save a binary file.
 */
void IkeWithPsALC_AD::getSteadySolnByTimeStepping(double *q, const int maxTimeStep, const double absTolTS, const double relTolTS, const double startCFL,
		const double *q_guess, double** q_tangent, const int NcontrolEqns /* =0 */,
		const double *lambda1 /* =NULL */, const double *lambda_guess /* =NULL */, double *lambda_tangent /* =NULL */,
		const double weightLambda /* =0.0 */, const double arcLength /* =0 */) {
	static bool firstCall = true;

	assert(q!=NULL && q_tangent!=NULL && q_guess != NULL);
	if(NcontrolEqns==0)
		assert(lambda_tangent==NULL);
	else
		assert(lambda_tangent!=NULL && arcLength!=0.0);

	assert(weightLambda == 1.0 || weightLambda == 0.0);  // Note: The formulation in time-stepping assumes that weightLambda (WEIGHT_LAMBDA_IN_ARCLENGTH) is always 1.0

	int m = nScal+5; // number of variables (rho, rhou, rhov, rhow, rhoE, and scalars)

	int countNegative = 0; // number of negative rho or press. If this is greater than 0, backtracking becomes active

	/********************/
	/*  Initialization  *
	 ********************/
	// Weight in the tangential condition
	// Note: The tangential equation is
	//           dLambda
	//         V ------- = -( q_tangent^T*(u-u_guess) + lambda_tangent*(lambda-lambda_guess) ) * sign(lambda_tangent),
	//              dt
	//         where V = 1 / weightTangentCond / |lambda_tangent|.
	//       With this setting (and weightTangentCond = 1), RHS is the same as the tangential condition of the original
	//       pseudo-arclength continuation.
	//       Note that in the code q_tangent and lambda_tangent have been already normalized before calling this method,
	//       i.e. q_tangent = q_tangent / arcLength,  lambda_tangent = lambda_tangent / arcLength.
	double vol_lambda = 1.0;
	if(NcontrolEqns > 0) {
		double weightTangentCond = getDoubleParam("WEIGHT_TANGENT_COND", "1.0"); // Small weight means that dLambda is small: Optimal value should be 1.0
		if(firstCall && mpi_rank==0)
			cout << "> WEIGHT_TANGENT_COND = " << weightTangentCond << endl;

//		double myVolSum = 0.0;
//		for(int icv=0; icv<ncv; ++icv)
//			myVolSum += cv_volume[icv];
//		double totVolSum;
//		MPI_Allreduce(&myVolSum, &totVolSum, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

		vol_lambda *= 1.0 / weightTangentCond / fabs(lambda_tangent[0]);  // Note: Just take the magnitude of the lambda_tangent of the first control Eqn.!!
		if(firstCall && mpi_rank==0)
			cout << "> VOL_LAMBDA = " << vol_lambda << endl;
	}

	// Set up MAX_DLAMBDA_MAG
    double MAX_DLAMBDA_MAG = 0.0;
    if(NcontrolEqns > 0)
    	MAX_DLAMBDA_MAG = 0.05 * arcLength * fabs(lambda_tangent[0]);  // Note: 5 percents of dLambda (Note: lambda_tangent has been already normalized)

	// Set iterTS to zero
	int iterTS = 0;

	int startingTSstep = 0;
	if(initThirdFromQ1) {
		if(step - startingStep == 2) {
			startingTSstep = getIntParam("INIT_THIRD_TS_STEP", "0");
			if(firstCall && mpi_rank==0)
				cout << "> INIT_THIRD_TS_STEP = " << startingTSstep << endl;
			iterTS += startingTSstep;

			MPI_Barrier(mpi_comm);
		}
	}

	/******************************************/
	/*  Time stepping:                        *
	 ******************************************/
	double *myResidual = new double[5+nScal];
	assert(Residual == NULL);
	Residual = new double[5+nScal];

	// RHS vectors are member variables (thus you don't need to define them) since they are also used in other methods.
	// On the other hand, the Jacobian matrix must be defined and assigned here.
	assert(RHSrho  != NULL);
	assert(RHSrhou != NULL);
	assert(RHSrhoE != NULL);
	if(nScal > 0)
		assert(RHSrhoScal != NULL);

	double (*A)[5][5] = new double[nbocv_s][5][5];
	double (*dq)[5]   = new double[ncv_g][5];
	double (*rhs)[5]  = new double[ncv][5];

	double ***AScal  = NULL;  if (nScal > 0) getMem3D(&AScal, 0, nScal-1, 0, 5, 0, nbocv_s-1, "AScal");
	double **dScal   = NULL;  if (nScal > 0) getMem2D(&dScal, 0, nScal-1, 0, ncv_g-1, "dScal");

	double *RHSlambda = NULL;
	double **ALambda  = NULL;
	double *dLambda   = NULL;
	if(NcontrolEqns > 0) {
		RHSlambda = new double[NcontrolEqns];
		getMem2D(&ALambda, 0, NcontrolEqns-1, 0, NcontrolEqns-1, "ALambda");
		dLambda = new double[NcontrolEqns];
	}

	//------------------------------------
	// some parameters that are only relevant for backward Euler
	//------------------------------------
	double underRelax = getDoubleParam("UNDER_RELAXATION", "0.5");

	if (!checkParam("LINEAR_SOLVER_NS_TRESHOLDS")) {
		ParamMap::add("LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=10  ABS_RESID=1.0e-8  REL_RESID=1.0e-4");    // add default values
		if (firstCall && mpi_rank == 0)
			cout << "WARNING: added keyword \"LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=10  ABS_RESID=1.0e-8  REL_RESID=1.0e-4\"" <<
			" to parameter map!" << endl;
	}
	int maxIterLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getInt("MAX_ITER");
	double zeroAbsLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("ABS_RESID");
	double zeroRelLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("REL_RESID");

	if (!checkParam("CFL_RAMP")) {
		ParamMap::add("CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0");    // add default: no increase of CFL number!
		if (firstCall && mpi_rank == 0)
			cout << "WARNING: added \"CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0\" to parameter map" << endl;
	}
	int startIncCFL = getParam("CFL_RAMP")->getInt("AFTER_ITER");
	int intervalIncCFL = getParam("CFL_RAMP")->getInt("INTERVAL_ITER");
	double incCFL = getParam("CFL_RAMP")->getDouble("FACTOR_CFL");
	double maxCFL = getParam("CFL_RAMP")->getDouble("MAX_CFL");

	//------------------------------------
	// Run the iterations for few more steps even after it converges
	//------------------------------------
	int moreTSstepsBelowConv_moreSteps;
	double moreTSstepsBelowConv_resid, moreTSstepsBelowConv_deltaQ;
	if (!checkParam("MORE_TS_STEPS_BELOW_CONV")) {
		ParamMap::add("MORE_TS_STEPS_BELOW_CONV  MORE_STEPS=0  RESID=1.0e-15  DELTA_Q=1.0e-15"); // add default values
		if (firstCall && mpi_rank == 0)
			cout<< "WARNING: added keyword \"MORE_TS_STEPS_BELOW_CONV  MORE_STEPS=0  RESID=1.0e-15  DELTA_Q=1.0e-15\""<< " to parameter map!" << endl;
	}
	moreTSstepsBelowConv_moreSteps = getParam("MORE_TS_STEPS_BELOW_CONV")->getInt("MORE_STEPS");
	moreTSstepsBelowConv_resid     = getParam("MORE_TS_STEPS_BELOW_CONV")->getDouble("RESID");
	moreTSstepsBelowConv_deltaQ    = getParam("MORE_TS_STEPS_BELOW_CONV")->getDouble("DELTA_Q");

	if(firstCall && mpi_rank==0)
		cout<<"> MORE_TS_STEPS_BELOW_CONV  MORE_STEPS="<<moreTSstepsBelowConv_moreSteps<<"  RESID="<<moreTSstepsBelowConv_resid<<"  DELTA_Q="<<moreTSstepsBelowConv_deltaQ<<endl;

	// -------------------------------------------------------------------------------------------
	// Note: The flow field (rho, rhou, rhoE, and scalars) should have been updated
	//       before calling this method
	// -------------------------------------------------------------------------------------------

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

	if ((maxTimeStep >= 0) && (iterTS >= maxTimeStep))  done = 1;

	if (initial_flowfield_output == "YES")
		writeData(step, maxTimeStep);

	// provide total runtime
	double wtime, wtime0;
	MPI_Barrier(mpi_comm);
	if (mpi_rank == 0)
		wtime = MPI_Wtime();

	cfl = startCFL;

	double initTotResidual = 0.0;
	double TotResidual;

	if(firstCall && mpi_rank == 0)
		cout<<"> Staring Time-stepping with TIME STEP MODE = " << timeStepMode << endl;

	int TSdumpQ1Interval = getIntParam("TS_DUMP_Q1_TECPLOT_INTERVAL", "0");

	while (done != 1) {
		iterTS++;
		if ((iterTS >= startIncCFL) && (iterTS%intervalIncCFL == 0) && (cfl < maxCFL))      cfl *= incCFL;
		double dt_min = calcDt(cfl);

		// find dt_max, too
		double dt_maxCPU = 0.0;
		for(int icv=0; icv<ncv; ++icv)
			dt_maxCPU = max(dt_maxCPU, local_dt[icv]);
		double dt_max;
		MPI_Allreduce(&dt_maxCPU, &dt_max, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);

		// ---------------------------------------------------------------------------------
		// Set matrices and vectors to zero
		// ---------------------------------------------------------------------------------
		for (int noc = 0; noc < nbocv_s; noc++)                             // set A, dq to zero! rhs is set zero in calcRHS
			for (int i = 0; i < 5; i++)
				for (int j = 0; j < 5; j++)
					A[noc][i][j] = 0.0;
		for (int icv = 0; icv < ncv_g; icv++)
			for (int i = 0; i < 5; i++)
				dq[icv][i] = 0.0;

		for (int iScal = 0; iScal < nScal; iScal++) {                        // set AScal, dScal to zero! rhs is set zero in calcRHS
			for (int i = 0; i <= 5; i++)
				for (int noc = 0; noc < nbocv_s; noc++)
					AScal[iScal][i][noc] = 0.0;
			for (int icv = 0; icv < ncv_g; icv++)
				dScal[iScal][icv] = 0.0;
		}

		for (int iEqn = 0; iEqn < NcontrolEqns; ++iEqn) {
			ALambda[iEqn][iEqn] = 0.0;
			dLambda[iEqn] = 0.0;
		}

		// ---------------------------------------------------------------------------------
		// calculate RHS for both NSE and scalars
		// ---------------------------------------------------------------------------------
		countNegative = calcRhs(RHSrho, RHSrhou, RHSrhoE, RHSrhoScal, A, AScal, true);
		if(countNegative > 0 && mpi_rank == 0)
			printf(">> WARNING! getSteadySolnByTimeStepping(): Negative rho, press, or kine occurs at total %d places\n", countNegative);

		bool useBarrier = true;
		if(useBarrier) {
			barrierSourceNS(RHSrho, RHSrhou, RHSrhoE, RHSrhoScal, A, AScal, true);     // Add barrier functions
			barrierSourceTurbScalars(RHSrhoScal, nScal, iterTS, ABSURDLY_BIG_NUMBER);  // Add barrier functions
		}

//		// Add penalty terms to RHS: This will make the solution not to go too far from the initial guess
//		double pentaltyStrength = 1.0; // This is a very important constant but hasn't been tested yet.
//
//		if(NcontrolEqns > 0) {
//			assert(q_tangent != NULL);
//			for(int iEqn = 0; iEqn < NcontrolEqns; ++iEqn) {
//				double tangentDotUvecFromUguess = calcUnweightedTangentVecDotFlowVecMinusGuess(iEqn, q_tangent, q_guess,
//						lambda_tangent, lambda, lambda_guess,
//						m, NcontrolEqns); // Note: tangentDotUvecFromUguess = [q_tangent; lambda_tangent]^T [flowVec - q_guess; lambda - lambda_guess]
//				                          //       is calculated here to avoid passing too many arguments to sourcePenalty().
//
//				// TO DO: sourcePenalty() has not been implemented yet!! Currently this method is empty!!
//				sourcePenalty(iEqn, pentaltyStrength, tangentDotUvecFromUguess, q_tangent, lambda_tangent, RHSrho, RHSrhou, RHSrhoE, RHSrhoScal, A, AScal, true);
//			}
//		}

		// ---------------------------------------------------------------------------------
		// calculate RHS for lambda
		// ---------------------------------------------------------------------------------
		if(NcontrolEqns > 0)
			calcRhsLambda(RHSlambda, ALambda, true, q_tangent, lambda_tangent, q_guess, lambda_guess, vol_lambda, NcontrolEqns, arcLength);

		// ---------------------------------------------------------------------------------
		// solve linear system for the NSE
		// ---------------------------------------------------------------------------------
		for (int icv=0; icv<ncv; ++icv) {                                // prepare rhs and A
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

		for (int icv=0; icv<ncv; icv++) {                                // update solution: Q_new = Q_old + Delta_Q
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

		for (int iScal = 0; iScal < nScal; iScal++) {                              // prepare rhs and A
			string scalname = scalarTranspEqVector[iScal].getName();
			for (int icv = 0; icv < ncv; ++icv) {
				//RHSrhoScal[iScal][icv] *= underRelax;
				RHSrhoScal[iScal][icv] *= scalarTranspEqVector[iScal].relax;

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
					RHSrhoScal[iScal][icv] = RHSrhoScal[iScal][icv]
					                      - AScal[iScal][0][noc] * dq[nbocv_v[noc]][0]
					                      - AScal[iScal][1][noc] * dq[nbocv_v[noc]][1]
					                      - AScal[iScal][2][noc] * dq[nbocv_v[noc]][2]
					                      - AScal[iScal][3][noc] * dq[nbocv_v[noc]][3]
					                      - AScal[iScal][4][noc] * dq[nbocv_v[noc]][4];
			}

			solveLinSysScalar(dScal[iScal], AScal[iScal][5], RHSrhoScal[iScal],
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
		// solve linear system for lambdas: every cpu calculate dLambda because it makes code simpler
		// -------------------------------------------------------------------------------------------
		// Note: The tangential equation is
		//           dLambda
		//         V ------- = -( q_tangent^T*(u-u_guess) + lambda_tangent*(lambda-lambda_guess) ) * sign(lambda_tangent),
		//              dt
		//         where V = 1 / weightTangentCond / |lambda_tangent|.
		//       With this setting (and weightTangentCond = 1), RHS is the same as the tangential condition of the original
		//       pseudo-arclength continuation.
		//       Note that in the code q_tangent and lambda_tangent have been already normalized before calling this method,
		//       i.e. q_tangent = q_tangent / arcLength,  lambda_tangent = lambda_tangent / arcLength.
//		double relaxLambda = 1.0; // This variable is defined here since it will be used for dq and dScal

		if(NcontrolEqns > 0) {
			for (int iEqn = 0; iEqn < NcontrolEqns; ++iEqn)
				RHSlambda[iEqn] *= underRelax;

			for (int iEqn = 0; iEqn < NcontrolEqns; ++iEqn) {
				double tmp = vol_lambda / dt_max;  // Take dt_max for lambda
				ALambda[iEqn][iEqn] += tmp;        // diagonal part ( vol/dt + A )
			}

			// move the other implicit terms to the RHS as you did for scalars: q_tangent*dq*sign(lambda_tangent)
			for (int iEqn = 0; iEqn < NcontrolEqns; ++iEqn) {
				double sign_lambda_tangent = (lambda_tangent[iEqn]>=0.0 ? 1.0 : -1.0);

				double flowImplicitTerms = -calcUnweightedTangentVecDotJOEdq(iEqn, q_tangent, lambda_tangent, dq, dScal, m) * sign_lambda_tangent;
				RHSlambda[iEqn] += flowImplicitTerms;
			}

			// Solve Linear system for lambda
			for (int iEqn = 0; iEqn < NcontrolEqns; ++iEqn)
				dLambda[iEqn] = RHSlambda[iEqn] / ALambda[iEqn][iEqn];   // Assuming that ALambda is a diagonal matrix

			// update lambdas and clip
			for (int iEqn = 0; iEqn < NcontrolEqns; ++iEqn) {
				if(fabs(dLambda[iEqn]) > MAX_DLAMBDA_MAG) {
					if(mpi_rank == mpi_size - 1)
						cout << "WARNING: dLambda[" << iEqn << "] = " << dLambda[iEqn] << " is too big: clip to ";
					dLambda[iEqn] = min(max(dLambda[iEqn], -MAX_DLAMBDA_MAG), MAX_DLAMBDA_MAG);  // note: -MAX_DLAMBDA_MAG <= dLambda <= MAX_DLAMBDA_MAG
					if(mpi_rank == mpi_size - 1)
						cout << dLambda[iEqn] << endl;
				}

				lambda[iEqn] += dLambda[iEqn];
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
		for (int i = 0; i < 5+nScal; i++) {
			myResidual[i] = 0.0;
			Residual[i] = 0.0;
		}

		for (int icv = 0; icv < ncv; icv++) {
			myResidual[0] += fabs(RHSrho[icv]);
			for (int i=0; i<3; i++)
				myResidual[i+1] += fabs(RHSrhou[icv][i]);
			myResidual[4] += fabs(RHSrhoE[icv]);
		}
		for (int iScal = 0; iScal < nScal; iScal++)
			for (int icv = 0; icv < ncv; icv++)
				myResidual[5+iScal] += fabs(RHSrhoScal[iScal][icv]/underRelax);

		MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);

		TotResidual = 0.0;
		for(int i=0; i<5+nScal; ++i)
			TotResidual += Residual[i];
		for (int iEqn = 0; iEqn < NcontrolEqns; ++iEqn)
			TotResidual += fabs(RHSlambda[iEqn]);

		// The following part is adapted from JoeWithModels::showResidue()
		if (iterTS%check_interval == 0 || iterTS == startingTSstep+1) {
			if ((mpi_rank == 0) && (iterTS%(check_interval*10) == 0))
				cout << "\ndone step: "<< iterTS << ", cfl: " << cfl << ", min. dt: " << dt_min << " time: " << time << endl;

			// residual label at every 10 output steps
			if ((mpi_rank == 0) && ((iterTS%(check_interval*10) == 0) || (iterTS == startingTSstep+1))) {
				printf("               rho         rhou-X      rhou-Y      rhou-Z      rhoE      ");
				for (int iScal = 0; iScal < nScal; iScal++)
					printf("%12s", scalarTranspEqVector[iScal].getName());
				if(NcontrolEqns > 0)
					printf("   RHSlambda      lambda");
				cout << endl;
			}

			// residual value at each output step
			if (mpi_rank == 0) {
				printf("RESID: %6d%12.4e%12.4e%12.4e%12.4e%12.4e", iterTS, Residual[0], Residual[1], Residual[2], Residual[3], Residual[4]);
				for (int iScal = 0; iScal < nScal; iScal++)
					printf("%12.4e", Residual[5+iScal]);
				if(NcontrolEqns > 0)
					printf("%12.4e%12.4e", fabs(RHSlambda[0]), lambda[0]);
				cout << endl;
			}
		}

		// IKJ: For output display (Metric to see the convergence of the simulation)
		if(log10_resid_rhoE != NULL)
			for (int icv = 0; icv < ncv; icv++)
				log10_resid_rhoE[icv] = log10(fabs(RHSrhoE[icv]) + MACHINE_EPS);
		if(log10_resid_scalar0 != NULL)
			for (int iScal = 0; iScal < nScal; iScal++)
				for (int icv = 0; icv < ncv; icv++)
					log10_resid_scalar0[icv] = log10(fabs(RHSrhoScal[iScal][icv]/underRelax) + MACHINE_EPS);
		updateCvDataG1G2(residField,       REPLACE_DATA);
		updateCvDataG1G2(log10_resid_rhoE, REPLACE_DATA);
		if(nScal>0)
			updateCvDataG1G2(log10_resid_scalar0, REPLACE_DATA);

		temporalHook();
//		dumpProbes(iterTS, 0.0);
//		writeData(step, iterTS);

		if(iterTS > startingTSstep+1 && TSdumpQ1Interval > 0) { // If the user wants to dump the data at every "newtonDumpQ1Interval" Newton step
			if(iterTS%TSdumpQ1Interval == 0) {
				if(NcontrolEqns>0) {
					char filename[50];
					sprintf(filename, "Q1_PT%05d.dumped.TS%05d.bin", step, iterTS);
					if(mpi_rank==0)
						cout<<"   > Dump the current Q vector on "<<filename<<": lambda[0]="<<std::setprecision(15)<<lambda[0]<<endl;
					writePsALCdumpedDataParallel(filename, step, arcLength, lambda1, lambda, NcontrolEqns, q);
				}

				writeData(step, iterTS);
			}
		}


		// ---------------------------------------------------------------------------------
		// save convergence history on a file
		// ---------------------------------------------------------------------------------
		if(mpi_rank==0) {
			if(iterTS == startingTSstep + 1) {
				FILE *fp;
				if(firstCall)
					fp = fopen(TIMESTEP_STATUS_FILENAME, "w");
				else
					fp = fopen(TIMESTEP_STATUS_FILENAME, "a");

				fprintf(fp, "CONTINUATION STEP = %d\n", step);
				fprintf(fp, "  ITER,  total     ,  rho       ,  rhou-X    ,  rhou-Y    ,  rhou-Z    ,  rhoE      ");
				for (int iScal = 0; iScal < nScal; iScal++)
					fprintf(fp, ",%12s", scalarTranspEqVector[iScal].getName());
				if(NcontrolEqns > 0)
					fprintf(fp, ",%12s,%12s", "tangential", "lambda");
				fprintf(fp, "\n");
				fclose(fp);
			}

			FILE *fp = fopen(TIMESTEP_STATUS_FILENAME, "a");
			fprintf(fp, "%6d,%12.4e,%12.4e,%12.4e,%12.4e,%12.4e,%12.4e", iterTS, TotResidual, Residual[0], Residual[1], Residual[2], Residual[3], Residual[4]);
			for (int iScal = 0; iScal < nScal; iScal++)
				fprintf(fp, ",%12.4e", Residual[5+iScal]);
			if(NcontrolEqns > 0)
				fprintf(fp, ",%12.4e,%12.4e", fabs(RHSlambda[0]), lambda[0]);
			fprintf(fp, "\n");
			fclose(fp);
		}

		// ---------------------------------------------------------------------------------
		// check tolerance
		// ---------------------------------------------------------------------------------
		double tolerance = absTolTS + initTotResidual * relTolTS;
		if ((iterTS >= maxTimeStep) || (TotResidual <= tolerance))   done = 1;

		// Update at the first iteration
		if(initTotResidual == 0.0)
			initTotResidual = TotResidual;
	}

	double tolerance = absTolTS + initTotResidual * relTolTS;
	if(TotResidual > tolerance) {
		if(mpi_rank == 0)
			cout << "TIME STEPPING has not been converged: maxTimeStep = " << maxTimeStep << ", total residual = " << TotResidual << " > tolerance = " << tolerance << endl;

		string booleanString = getStringParam("AFTER_TIMESTEPPING_LAUNCH_NT", "FALSE");
		std::transform(booleanString.begin(), booleanString.end(), booleanString.begin(), ::tolower);
		bool AfterTSlaunchNT = (booleanString.compare("true")==0 || booleanString.compare("yes")==0);
		if(!AfterTSlaunchNT)
			throw(PSALC_ERROR_CODE);
	}

	MPI_Barrier(mpi_comm);
	if (mpi_rank == 0) {
		double wtime0 = wtime;
		wtime = MPI_Wtime();
		cout << " > runtime for iterations[s]: " << wtime - wtime0 << endl;
	}

//	// ---------------------------------------------------------------------------------
//	// output
//	// ---------------------------------------------------------------------------------
//
//	temporalHook();
//	finalHook();
//
//	writeRestart(step);
//	writeData(step, iterTS);

	// ---------------------------------------------------------------------------------
	// Update q
	// ---------------------------------------------------------------------------------
	for (int icv = 0; icv < ncv; ++icv) {
		int indexStart = icv*m;
		q[indexStart] = rho[icv];
		for(int i=0; i<3; ++i)
			q[indexStart+1+i] = rhou[icv][i];
		q[indexStart+4] = rhoE[icv];

		for (int iScal = 0; iScal < nScal; iScal++) {
			double *phi = scalarTranspEqVector[iScal].phi;
			q[indexStart+5+iScal] = phi[icv];
		}
	}
	if(NcontrolEqns > 0 && mpi_rank == mpi_size-1) {
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
			q[ncv*m+iEqn] = lambda[iEqn];
	}

	// ---------------------------------------------------------------------------------
	// delete memory
	// ---------------------------------------------------------------------------------

	delete [] A;
	delete [] rhs;

	delete [] dq;

	delete [] myResidual;
	delete [] Residual; 	Residual = NULL;

	if (nScal > 0) freeMem3D(AScal, 0, nScal-1, 0, 5, 0, nbocv_s-1);
	if (nScal > 0) freeMem2D(dScal, 0, nScal-1, 0, ncv_g-1);

	delete [] RHSlambda;
	freeMem2D(ALambda, 0, NcontrolEqns-1, 0, NcontrolEqns-1);
	delete [] dLambda;

    if (petscSolver != NULL)        	{ delete petscSolver;         petscSolver = NULL; }
    if (petscSolverScalars != NULL) 	{ delete petscSolverScalars;  petscSolverScalars = NULL; }

	firstCall = false;

	MPI_Barrier(mpi_comm);
}

/*
 * Method: sourcePenalty
 * ---------------------
 * Add penalty to prevent diverging solution
 */
void IkeWithPsALC_AD::sourcePenalty(const int iEqn, const double pentaltyStrength, const double tangentDotUvecFromUguess,
		double** q_tangent, const double *lambda_tangent,
		double *RHSrho, double (*RHSrhou)[3], double *RHSrhoE, double **RHSrhoScal, double (*A)[5][5], double ***AScal, const bool flagImplicit) {
	int nScal = scalarTranspEqVector.size();
	for(int icv=0; icv<ncv; ++icv) {
		int indexStart = icv*(5+nScal);
		RHSrho[icv] -= pentaltyStrength * tangentDotUvecFromUguess * q_tangent[iEqn][indexStart];
		for(int i=0; i<3; ++i)
			RHSrhou[icv][i] -= pentaltyStrength * tangentDotUvecFromUguess * q_tangent[iEqn][indexStart+1+i];
		RHSrhoE[icv] -= pentaltyStrength * tangentDotUvecFromUguess * q_tangent[iEqn][indexStart+4];
		for(int iScal = 0; iScal <nScal; ++iScal)
			RHSrhoScal[iScal][icv] -= pentaltyStrength * tangentDotUvecFromUguess * q_tangent[iEqn][indexStart+5+iScal];
	}

	if(flagImplicit) { // Actually, it must be implicit...
		for(int icv=0; icv<ncv; ++icv) {
			int indexStart = icv*(5+nScal);
			int noc00 = nbocv_i[icv];

			for(int i=0; i<5; ++i)
				A[noc00][i][i] += pentaltyStrength * pow(q_tangent[iEqn][indexStart+i], 2.0);

			for (int iScal = 0; iScal < nScal; iScal++)
				AScal[iScal][5][noc00] += pentaltyStrength * pow(q_tangent[iEqn][indexStart+5+iScal], 2.0);
					// Note: The size of AScal for implicit Euler = [nScal][6][nbocv_s]
		}
	}
}

/*
 * Method: calcRhsLambda
 * ---------------------
 * Calculate RHS for lambda evolution equation
 */
void IkeWithPsALC_AD::calcRhsLambda(double* RHSlambda, double** ALambda, const bool flagImplicit,
		double** q_tangent, const double* lambda_tangent, const double* q_initGuess, const double* lambda_initGuess,
		const double vol_lambda, const int NcontrolEqns, const double arcLength) {
	assert(RHSlambda!=NULL && q_tangent!=NULL && lambda_tangent!=NULL && q_initGuess!=NULL && lambda_initGuess!=NULL);
	assert(NcontrolEqns != 0);

	// Set RHS zero
	for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
		RHSlambda[iEqn] = 0.0;

	// Solve RHS
	// Note: The tangential equation is
	//           dLambda
	//         V ------- = -( q_tangent^T*(u-u_guess) + lambda_tangent*(lambda - lambda_guess) ) * sign(lambda_tangent),
	//              dt
	//         where V = 1 / weightTangentCond / |lambda_tangent|.
	//       Note that in the code q_tangent and lambda_tangent have been already normalized before calling this method,
	//       i.e. q_tangent = q_tangent / arcLength,  lambda_tangent = lambda_tangent / arcLength.
	for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
		double tangentDotUvecFromUguess = calcUnweightedTangentVecDotFlowVecMinusGuess(iEqn, q_tangent, q_initGuess,
				lambda_tangent, lambda, lambda_initGuess,
				5+nScal, NcontrolEqns);
//		RHSlambda[iEqn] -= arcLength * tangentDotUvecFromUguess * ((lambda_tangent[iEqn]>=0.0) ? 1.0 : -1.0);
		RHSlambda[iEqn] -= tangentDotUvecFromUguess * ((lambda_tangent[iEqn]>=0.0) ? 1.0 : -1.0);

		if(flagImplicit) {
			assert(ALambda != NULL);
//			ALambda[iEqn][iEqn] += arcLength * lambda_tangent[iEqn] * (lambda_tangent[iEqn]>=0.0 ? 1.0 : -1.0); // Note: In Joe, A is actually minus of Jacobian
			ALambda[iEqn][iEqn] += lambda_tangent[iEqn] * (lambda_tangent[iEqn]>=0.0 ? 1.0 : -1.0); // Note: In Joe, A is actually minus of Jacobian
		}
	}
}

/*
 * Method: getSteadySolnByBackwardEuler
 * ------------------------------------
 * Original code = runBackwardEuler() in JoeWithModels.cpp
 * Note: this method cannot update lambda
 *
 * Return:
 *   By value    : isConverged
 *   By reference: dq, dScal, rhs, rhsScal, A, Ascal, initTotResidual, Residual
 */
bool IkeWithPsALC_AD::getSteadySolnByBackwardEuler(double (*dq)[5], double **dScal, double (*rhs)[5], double **rhsScal,
		double (*A)[5][5], double ***AScal, double &initTotResidual, double *Residual,
		const int maxStepEuler, const double absTolEuler, const double relTolEuler,
		const int maxIterLS, const double zeroAbsLS, const double zeroRelLS,
		double &cfl, const int startIncCFL, const int intervalIncCFL, const double incCFL, const double maxCFL,
		const int check_interval) {
	int nScal = scalarTranspEqVector.size();

	assert(dq!=NULL && rhs!=NULL && Residual!=NULL);
	if(nScal>0)
		assert(AScal!=NULL && dScal!=NULL && rhsScal!=NULL);
	else
		assert(AScal==NULL && dScal==NULL && rhsScal==NULL);

	// Note: RHSrho, RHSrhou, RHSrhoE (and rhsScal) are used in calcRhs()
	//       rhs, and rhsScal are used when solving the linear system
	assert(RHSrho != NULL && RHSrhou != NULL && RHSrhoE != NULL);

	bool isConverged = false;

	double *myResidual = new double[5+nScal];

	//------------------------------------
	// some parameters that are only relevant for backward euler
	//------------------------------------
	double underRelax = getDoubleParam("UNDER_RELAXATION", "1.0");

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
	for(int ifa=nfa; ifa<nfa_b2; ++ifa)
		setBC1D(ifa, UgpWithCvCompFlow::rho, UgpWithCvCompFlow::rhou, UgpWithCvCompFlow::rhoE); // Note: ncvggf's cannot be not updated, so the user should take care of it

	// -------------------------------------------------------------------------------------------
	// update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
	// -------------------------------------------------------------------------------------------
	calcRansTurbViscMuet();

	// -------------------------------------------------------------------------------------------
	//
	//   Loop over time steps
	//
	// -------------------------------------------------------------------------------------------

	int step = 0;   // set to zero, even if read from restart file, as for steady state time step doesn't matter

	int done = 0;
	if (maxStepEuler == 0)
		done = 1;

	double totResidual;

//	if (initial_flowfield_output == "YES")
//		writeData(0);

	// provide total runtime
	double wtime, wtime0;
	MPI_Barrier(mpi_comm);
	if (mpi_rank == 0)
		wtime = MPI_Wtime();

	while (done != 1) {
		step++;
		if ((step >= startIncCFL) && (step%intervalIncCFL == 0) && (cfl < maxCFL))
			cfl *= incCFL;
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

		// For tecplot output
		if(nScal>0)
			convertRHSrhoScalToSeparatedRhs(rhsScal);

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

				double tmp = cv_volume[icv]/(local_dt[icv]);

				if(iScal == getScalarTransportIndex("f")) {
					// do nothing
				}else{
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
		for(int ifa=nfa; ifa<nfa_b2; ++ifa)
			setBC1D(ifa, UgpWithCvCompFlow::rho, UgpWithCvCompFlow::rhou, UgpWithCvCompFlow::rhoE); // Note: ncvggf's cannot be not updated, so the user should take care of it

		// -------------------------------------------------------------------------------------------
		// update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
		// -------------------------------------------------------------------------------------------
		calcRansTurbViscMuet();


		// =========================================================================================
		// calculate and show residual
		// =========================================================================================
		for (int i = 0; i < 5+nScal; i++) {
			myResidual[i] = 0.0;
			Residual[i] = 0.0;
		}

		for (int icv = 0; icv < ncv; icv++) {
			myResidual[0] += fabs(RHSrho[icv]);
			for (int i=0; i<3; i++)
				myResidual[i+1] += fabs(RHSrhou[icv][i]);
			myResidual[4] += fabs(RHSrhoE[icv]);
		}
		for (int iScal = 0; iScal < nScal; iScal++)
			for (int icv = 0; icv < ncv; icv++)
				myResidual[5+iScal] += fabs(rhsScal[iScal][icv]/underRelax);

		MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);

		totResidual = 0.0;
		for(int i=0; i<5+nScal; ++i)
			totResidual = Residual[i];

		if(step==1)
			initTotResidual = totResidual;

		if(debugLevel > 0) {
			if (step%check_interval == 0 || step==1) {
				if ((mpi_rank == 0) && (step%(check_interval*10) == 0))
					cout << "\n          done step: "<< step << ", cfl: " << cfl << ", min. dt: " << dt_min << " time: " << time << endl;

				showResidue2(Residual, step, check_interval);
			}
		}

		temporalHook();
//		dumpProbes(step, 0.0);

		if (step >= maxStepEuler)
			done = 1;
		if (totResidual <= absTolEuler+relTolEuler*initTotResidual) {
			done = 1;
			isConverged = true;
		}
	}

	MPI_Barrier(mpi_comm);
	if (mpi_rank == 0) {
		double wtime0 = wtime;
		wtime = MPI_Wtime();
		cout << "          > runtime for iterations[s]: " << wtime - wtime0 << endl;
	}

	// ---------------------------------------------------------------------------------
	// output
	// ---------------------------------------------------------------------------------
	temporalHook();
	finalHook();

	// ---------------------------------------------------------------------------------
	// delete memory
	// ---------------------------------------------------------------------------------
	delete [] myResidual;

	return isConverged;
}

/*
 * Method: calcJacobian1DAD
 * ------------------------
 * calculate Rhs and get the Jacobian matrix using AD using the so-called "1D style"
 * Return: If backtracking line search is required (due to negative density, etc.), return true.
 * Original code of Rhs part = JoeWithModels_AD::calcResidualDerivative(double***, double***, int)
 */
int IkeWithPsALC_AD::calcJacobian1DAD(MatComprsedSTL &jacMatrixSTL, double *rhsSingleArray, const int debugLevel, const int NcontrolEqns) {
	assert(lambda_AD == NULL);

	//==============================
	// Initialization step-1 (BASIC): required for both RHS calculation and Jacobian calculation
	//==============================
	int nScal = scalarTranspEqVector.size();
	int nVars = 5+nScal;

	int myCountReducedOrder = 0;
	int CountReducedOrder;

	// booleans for first-call
	static bool firstCall = true;
	static bool firstCallScalarTurb = true, firstCallScalarComb = true;

#ifndef USE_MEM_SAVING_ADVAR_1D_
	// allocate memory for flow variables: e.g. rho_AD, rhou_AD, vel, mul_fa, ...
	if(NcontrolEqns==0)
		UgpWithCvCompFlow_AD::initialize_adoubles();
	else {
		initialize_adoubles(); // The pseudo-arclength continuation method requires additional AD variables for the control eqns
		assert(lambda_AD != NULL);
		assert(lambda != NULL);
	}
#endif

	// variables for checking wall-clock time
	double myWtime0, myWtime1;
	myWTimeJacCalc.wallTimeRHSinit  = 0.0;
	myWTimeJacCalc.wallTimeRHScalc  = 0.0;
	myWTimeJacCalc.wallTimeRHSfinal = 0.0;
	myWTimeJacCalc.wallTimeJac      = 0.0;
	myWTimeJacCalc.wallTimeCleanup  = 0.0;

	//==============================
	// Initialization step-2 (JACOBI): required for only Jacobian calculation
	//==============================
	// allocate memory for indev and initialize it with the flow variables
#ifndef USE_MEM_SAVING_ADVAR_1D_
	double *indvec = new double [ncv_gg*nVars+NcontrolEqns];
	for (int icvTemp = 0; icvTemp < ncv_gg; icvTemp++) {
		int icv_i = nVars*icvTemp;
		indvec[icv_i]   = rho[icvTemp];
		indvec[icv_i+1] = rhou[icvTemp][0];
		indvec[icv_i+2] = rhou[icvTemp][1];
		indvec[icv_i+3] = rhou[icvTemp][2];
		indvec[icv_i+4] = rhoE[icvTemp];
		for(int iScal=0; iScal<nScal; iScal++)
			indvec[icv_i+5+iScal] = scalarTranspEqVector[iScal].phi[icvTemp];
	}
	for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
		indvec[nVars*ncv_gg+iEqn] = lambda[iEqn];
#endif

	// allocate memory for the Jacobian matrix
	jacMatrixSTL.setMatSize(ncv_gg*nVars+NcontrolEqns, ncv*nVars);
	jacMatrixSTL.initMat();

	//==============================
	// Looping over the all the ICVs
	//==============================
	myBarrierMassSourceSum1D_AD   = 0.0; // Note: The barrier source term will be added to this variable in barrierSourceNS1D_AD()
	myBarrierEnergySourceSum1D_AD = 0.0; // Note: The barrier source term will be added to this variable in barrierSourceNS1D_AD()

	warningsVectorIkeUgp1D.clear();
	warningsVectorIkeModels1D.clear(); // Note: If the simulation is about to break down for some reasons,
                                       //      there will be a massive amount of warning/error messages.
                                       //       Thus, this "vector<std::string> warningsVectorIkeModels1D" object
	                                   //      defined in IkeWithModels.h stores the (only) warning messages so that
	                                   //      the user can manipulate the amount of messages that should be
	                                   //      printed out on the screen.

	myNonDiagAdd_count     = 0;   // Collect the statistics of artificially added diagonal entries
	myNonDiagAdd_rhsAbsSum = 0.0;

	for(int icv=0; icv<ncv; ++icv) {
		int tag = mpi_rank;

		// +++++++++++++++++++++++++++++++++
		// Calculate RHS
		// +++++++++++++++++++++++++++++++++
		calcJacobian1DAD_calcRhs(rhsSingleArray, myCountReducedOrder, myWTimeJacCalc, icv, tag, nScal, debugLevel, NcontrolEqns,
				firstCall, firstCallScalarTurb, firstCallScalarComb);

		// +++++++++++++++++++++++++++++++++
		// Calculate Jacobian
		// +++++++++++++++++++++++++++++++++
		myWtime0 = MPI_Wtime();

		/* indvec for the MEM_SAVING mode */
#ifdef USE_MEM_SAVING_ADVAR_1D_
		int n_nbocv2gg = nbocv2_eachIcv.size();
		assert(n_nbocv2gg > 0);
		double *indvec = new double [n_nbocv2gg*nVars+NcontrolEqns];
		for (size_t index=0; index<n_nbocv2gg; ++index) {
			int icv = nbocv2_eachIcv[index];
			int arrayStartIndex = nVars*index;
			indvec[arrayStartIndex] = rho[icv];
			for(int i=0; i<3; ++i)
				indvec[arrayStartIndex+1+i] = rhou[icv][i];
			indvec[arrayStartIndex+4] = rhoE[icv];
			for(int iScal=0; iScal<nScal; ++iScal)
				indvec[arrayStartIndex+5+iScal] = scalarTranspEqVector[iScal].phi[icv];
		}
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
			indvec[n_nbocv2gg*nVars+iEqn] = lambda[iEqn];
#endif

		/* coordinate format for Jacobian */
		unsigned int *rind   = NULL;        /* row indices    */
		unsigned int *cind   = NULL;        /* column indices */
		double       *values = NULL;        /* values         */
		int          nnz     = 0;

		int options[8]; /* option for Jacobian computation (ADOL-C) */
		options[0] = 0; /* sparsity pattern computation by index domains (default) */
		options[1] = 0; /*                                     safe mode (default) */
		options[2] = 0; /*                          not required if options[0] = 0 */
		options[3] = 1; /*     way of compression: 0 = column compression, 1 = row */

#ifdef USE_MEM_SAVING_ADVAR_1D_
		sparse_jac(tag, nVars, nbocv2_eachIcv.size()*nVars+NcontrolEqns, 0, indvec, &nnz, &rind, &cind, &values, options); // To DO
#else
		sparse_jac(tag, nVars, ncv_gg*nVars+NcontrolEqns,                0, indvec, &nnz, &rind, &cind, &values, options);
#endif
				/* Arguments:
				 *   tag = tape identification
				 *   dims1 = number of dependent variables
				 *   dims2 = number of independent variables
				 *   repeat = indicate repeated call at same argument
				 *            repeat=0 : The functions are called at a point with a new sparsity structure
				 *            repeat=1 : Re-usage of the sparsity pattern from the previous call.
				 *   indvec[dims2] = independent vector
				 *   nnz = number of nonzeros
				 *   rind[nnz] = row index
				 *   cind[nnz] = column index
				 *   values[nnz] = non-zero values
				 *   option[4] = array of control parameters
				 * Note: The sparse_jac() function calls ColPack::JacobianRecovery1D::RecoverD2Row_CoordinateFormat_unmanaged()
				 *       And rind, cind, values are allocated by the following commands in the function (see JacobianRecovery1D.cpp):
				 *           (*ip2_RowIndex) = (unsigned int*) malloc(numOfNonZeros * sizeof(unsigned int));
				 *           (*ip2_ColumnIndex) = (unsigned int*) malloc(numOfNonZeros * sizeof(unsigned int));
				 *           (*dp2_JacobianValue) = (double*) malloc(numOfNonZeros * sizeof(double));
				 *       However, the output is NOT managed by ColPack.
				 *       Thus, the user should free the output manually using free() (NOT delete) function when it is no longer needed.
				 */

		if(icv==0 && mpi_rank==0 && (debugLevel>1 || (debugLevel>0 && firstCall))) {
			int icv_i = nVars*icv;

			cout<<"================================================"<<endl;
			cout<<"Jacobian stats (mpi_rank=0, icv=0)"<<endl;
			cout<<"================================================"<<endl;
			cout<<"(Rhs values =";
			for(int ii=0; ii<5+nScal; ++ii) printf("  %.4e", rhsSingleArray[icv_i+ii]);
			cout<<")"<<endl;
			cout<<"Number of independents = "<<ncv_gg*nVars+NcontrolEqns<<endl;
			cout<<"Number of dependents   = "<<nVars<<endl;
			cout<<"Number of non-zeros    = "<<nnz<<endl;
			cout<<"================================================"<<endl<<endl;
		}

		// +++++++++++++++++++++++++++++++++
		// Copy the calculated Jacobian to the data container
		// +++++++++++++++++++++++++++++++++
		int firstRindReal = icv*nVars;
		jacMatrixSTL.blockCopy(nnz, firstRindReal, rind, cind, values);

		jacMatrixSTL.convert_blockCindices(nbocv2_eachIcv, jacMatrixSTL.get_nnz()-nnz, nnz, nVars, nVars, NcontrolEqns, ncv_gg);

		myWtime1 = MPI_Wtime();
		myWTimeJacCalc.wallTimeJac += myWtime1-myWtime0;
		myWtime0 = myWtime1;

		// +++++++++++++++++++++++++++++++++
		// Free memory
		// +++++++++++++++++++++++++++++++++
		removeTape(tag, ADOLC_REMOVE_COMPLETELY);  // Note: "int removeTape(short tapeID, short type)" in tape_handling.cpp
		                                           //       If type==ADOLC_REMOVE_COMPLETELY, remove the tape instead of free it.
		                                           //       Returns 0 if the removing process is completed okay.

#ifdef USE_MEM_SAVING_ADVAR_1D_
		destroy_adoubles(NcontrolEqns);

		finalHookScalarRansTurbModel1D_AD();  // Reset internal pointers (e.g. kine, omega, etc) in a turb model
		finalHookScalarRansCombModel1D_AD();  // Reset internal pointers (e.g. ZMean, Zvar, CMean, etc) in a turb model

		delete [] indvec;

		nbocv2_eachIcv.clear();
		nbocv2ff_eachIcv.clear();
		fa2_eachIcv.clear();
#endif

		free(rind); 	//rind = NULL;
		free(cind); 	//cind = NULL;
		free(values); 	//values = NULL;

		if(mpi_rank==0 && icv%50==0) {
			FILE *fp = fopen(JAC1D_STATUS_FILENAME, "a");
			fprintf(fp, "%4d  %12d  %10d   %10.3f\n", step, iterNewton, icv, MPI_Wtime()-initWtime);
			fclose(fp);
		}

		myWtime1 = MPI_Wtime();
		myWTimeJacCalc.wallTimeCleanup += myWtime1-myWtime0;

		if(firstCall && debugLevel>0 && mpi_rank==0) {
			if(icv==100-1) {
				cout<<"     > Runtime for the first   100 ICVs [sec] = "<<
						myWTimeJacCalc.wallTimeRHSinit + myWTimeJacCalc.wallTimeRHScalc + myWTimeJacCalc.wallTimeRHSfinal +myWTimeJacCalc.wallTimeJac + myWTimeJacCalc.wallTimeCleanup <<endl;
			} else if(icv==1000-1) {
				cout<<"     > Runtime for the first  1000 ICVs [sec] = "<<
						myWTimeJacCalc.wallTimeRHSinit + myWTimeJacCalc.wallTimeRHScalc + myWTimeJacCalc.wallTimeRHSfinal +myWTimeJacCalc.wallTimeJac + myWTimeJacCalc.wallTimeCleanup <<endl;
			} else if(icv==10000-1) {
				cout<<"     > Runtime for the first 10000 ICVs [sec] = "<<
						myWTimeJacCalc.wallTimeRHSinit + myWTimeJacCalc.wallTimeRHScalc + myWTimeJacCalc.wallTimeRHSfinal +myWTimeJacCalc.wallTimeJac + myWTimeJacCalc.wallTimeCleanup <<endl;
			} else if(icv==100000-1) {
				cout<<"     > Runtime for the first 100000 ICVs [sec] = "<<
						myWTimeJacCalc.wallTimeRHSinit + myWTimeJacCalc.wallTimeRHScalc + myWTimeJacCalc.wallTimeRHSfinal +myWTimeJacCalc.wallTimeJac + myWTimeJacCalc.wallTimeCleanup <<endl;
			} else if(icv==ncv-1) {
				cout<<"     > Runtime for ncv="<<ncv<<" ICVs [sec] = "<<
						myWTimeJacCalc.wallTimeRHSinit + myWTimeJacCalc.wallTimeRHScalc + myWTimeJacCalc.wallTimeRHSfinal +myWTimeJacCalc.wallTimeJac + myWTimeJacCalc.wallTimeCleanup <<endl;
			}
		}
	}

	jacMatrixSTL.finalizeMat(); // Set jacMatrixSTL.ncols_eachRow_i for empty rows on the bottom of the matrix

	// Show some statistics to the user
//	double totBarrierMassSourceSum1D_AD;
//	double totBarrierEnergySourceSum1D_AD;
//	MPI_Allreduce(&myBarrierMassSourceSum1D_AD,   &totBarrierMassSourceSum1D_AD,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
//	MPI_Allreduce(&myBarrierEnergySourceSum1D_AD, &totBarrierEnergySourceSum1D_AD, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
//
//	bool useBarrierNS = checkBarrierParamNS(iterNewton);
//	if(debugLevel>0 && mpi_rank==0) {
//		if(useBarrierNS)
//			printf("              >> Barrier: Total mass source = %.4e, Total energy source = %.4e \n", totBarrierMassSourceSum1D_AD, totBarrierEnergySourceSum1D_AD);
//		else {
//			if(debugLevel>1)
//				printf("              >> Barrier is not used for NS \n");
//		}
//	}

	probeInCalcJacobian1DAD(jacMatrixSTL, rhsSingleArray, debugLevel, NcontrolEqns);

	// free memory allocated for flow variables
#ifndef USE_MEM_SAVING_ADVAR_1D_
	// free indvec
	delete [] indvec;
#endif

#ifdef USE_MEM_SAVING_ADVAR_1D_
	delete [] scalarTranspEqVector_AD;
	scalarTranspEqVector_AD = NULL;
#else
	if(NcontrolEqns==0)
		UgpWithCvCompFlow_AD::destroy_adoubles();
	else
		destroy_adoubles();
#endif

	if(firstCall && debugLevel>0 && mpi_rank==0) {
		cout<<"   > calcJacobian1DAD(): Runtime for the RHS calculations [sec] = "<<myWTimeJacCalc.wallTimeRHSinit+myWTimeJacCalc.wallTimeRHScalc+myWTimeJacCalc.wallTimeRHSfinal<<endl;
		printf("                           Runtime Details [sec]\n");
		printf("                           1. Initialize AD     =%9.4f\n", myWTimeJacCalc.wallTimeRHSinit);
		printf("                           2. Calc RHS          =%9.4f\n", myWTimeJacCalc.wallTimeRHScalc);
		printf("                           3. Finalize AD       =%9.4f\n", myWTimeJacCalc.wallTimeRHSfinal);
		cout<<"   > calcJacobian1DAD(): Runtime for the Jac evaluations  [sec] = "<<myWTimeJacCalc.wallTimeJac<<endl;
		cout<<"   > calcJacobian1DAD(): Runtime for the clean-up process [sec] = "<<myWTimeJacCalc.wallTimeCleanup<<endl;
	}

	//
	firstCall = false;
	firstCallScalarTurb = false;
	firstCallScalarComb = false;

	MPI_Allreduce(&myCountReducedOrder, &CountReducedOrder, 1, MPI_INT, MPI_SUM, mpi_comm);

	return CountReducedOrder;
}

///*
// * Method: calcJacobian1DAD_new
// * ----------------------------
// * Calculate Rhs and get the Jacobian matrix using AD using the so-called "1D style"
// * For each CV, the adouble memory is allocated only for the neighboring CVs
// *
// * Return: If backtracking line search is required (due to negative density, etc.), return true.
// * Original code of Rhs part = JoeWithModels_AD::calcResidualDerivative(double***, double***, int)
// */
//bool IkeWithPsALC_AD::calcJacobian1DAD_new(MatComprsedSTL &jacMatrixSTL, double *rhsSingleArray, bool printStats, const int NcontrolEqns) {
//	assert(lambda_AD == NULL);
//
//	int nScal = scalarTranspEqVector.size();
//	int nVars = 5+nScal;
//
//	int myCountReducedOrder = 0;
//	int CountReducedOrder;
//	static bool firstCall = true;
//
//	double wtime0, wtimeRHS, wtimeJac;
//	double runtimeRHS = 0.0, runtimeJac = 0.0;
//
//	// allocate memory for flow variables: e.g. rho_AD, rhou_AD, vel, mul_fa, ...
//	if(NcontrolEqns==0)
//		UgpWithCvCompFlow_AD::initialize_adoubles();
//	else {
//		initialize_adoubles(); // The pseudo-arclength continuation method requires additional AD variables for the control eqns
//		assert(lambda_AD != NULL);
//	}
//
////	// initialize the flow field for the base variables
////	initialHook();
//
//	// allocate memory for indev and initialize it with the flow variables
//	double *indvec = new double [ncv_gg*nVars+NcontrolEqns];
//	for (int icvTemp = 0; icvTemp < ncv_gg; icvTemp++) {
//		int icv_i = nVars*icvTemp;
//		indvec[icv_i]   = rho[icvTemp];
//		indvec[icv_i+1] = rhou[icvTemp][0];
//		indvec[icv_i+2] = rhou[icvTemp][1];
//		indvec[icv_i+3] = rhou[icvTemp][2];
//		indvec[icv_i+4] = rhoE[icvTemp];
//		for(int iScal=0; iScal<nScal; iScal++)
//			indvec[icv_i+5+iScal] = scalarTranspEqVector[iScal].phi[icvTemp];
//	}
//	for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
//		indvec[nVars*ncv_gg+iEqn] = lambda[iEqn];
//
//	// booleans for first-call
//	bool firstCallScalarTurb = true, firstCallScalarComb = true;
//
//	// allocate memory for the Jacobian matrix
//	jacMatrixSTL.setMatSize(ncv_gg*nVars+NcontrolEqns, ncv*nVars);
//	jacMatrixSTL.initMat();
//
//if(mpi_rank==0) cout<<"IkeWithPsALC_AD::calcJacobian1DAD_new(): checking wTime() -- also in calcResidual1D_AD()"<<endl;
//double myWtimeAllocADvars = 0.0;
//double myWtimeInitADvars = 0.0;
//double myWtimeTraceOn = 0.0;
//double myWtimeAssgnIndep = 0.0;
//double myWtimeInitHooks = 0.0;
//double myWtimeFinalADvars = 0.0;
//double myWtimeTraceOff = 0.0;
//// The following variables are defined in IkeWithModels.h and calculated in calcResidual1D_AD()
//myWtimeFindFaces = 0.0;
//myWtimeStateVars = 0.0;
//myWtimeMaterPropert = 0.0;
//myWtimeSetBC = 0.0;
//myWtimeCalcGrad = 0.0;
//myWtimeSetupZoneArray = 0.0;
//myWtimeCalcTurbVisc = 0.0;
//myWtimeCalcRhs = 0.0;
//
//	for(int icv=0; icv<ncv; ++icv) {
//		if(firstCall && debugLevel>0 && mpi_rank==0) {
//			wtime0 = MPI_Wtime();
//		}
//double myWtime0 = MPI_Wtime();
//double myWtime1;
//
//		// +++++++++++++++++++++++++++++++++
//		// Allocate memory for Flow residual
//		// +++++++++++++++++++++++++++++++++
//		REALA rhs_rho_AD;
//		REALA rhs_rhou_AD[3];
//		REALA rhs_rhoE_AD;
//		REALAS *rhs_rhoScal_AD = new REALA [nScal];
//myWtime1 = MPI_Wtime();
//myWtimeAllocADvars += myWtime1-myWtime0;
//myWtime0 = myWtime1;
//
//		rhs_rho_AD = 0.0;
//		for(int i=0; i<3; ++i)
//			rhs_rhou_AD[i] = 0.0;
//		rhs_rhoE_AD = 0.0;
//		for(int i=0; i<nScal; ++i)
//			rhs_rhoScal_AD[i] = 0.0;
//myWtime1 = MPI_Wtime();
//myWtimeAssgnIndep += myWtime1-myWtime0;
//myWtime0 = myWtime1;
//
//		// +++++++++++++++++++++++++++++++++
//		// Start the calculation of Residual wrt flow
//		// +++++++++++++++++++++++++++++++++
//		int tag = mpi_rank;
//		try {
//			trace_on(tag, 1);
//myWtime1 = MPI_Wtime();
//myWtimeTraceOn += myWtime1-myWtime0;
//myWtime0 = myWtime1;
//
//	//		// initialize the flow field for the base variables
//	//		initialHook();
//
//			// Independent variables
//			for (int icvTemp = 0; icvTemp < ncv_gg; icvTemp++) {
//				rho_AD[icvTemp]      <<= rho[icvTemp] ;
//				rhou_AD[icvTemp][0]  <<= rhou[icvTemp][0] ;
//				rhou_AD[icvTemp][1]  <<= rhou[icvTemp][1] ;
//				rhou_AD[icvTemp][2]  <<= rhou[icvTemp][2] ;
//				rhoE_AD[icvTemp]     <<= rhoE[icvTemp] ;
//				for(int iScal=0; iScal<nScal; iScal++)
//					scalarTranspEqVector_AD[iScal].phi[icvTemp] <<= scalarTranspEqVector[iScal].phi[icvTemp];
//			}
//			if(NcontrolEqns > 0)
//				for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
//					lambda_AD[iEqn] <<= lambda[iEqn];
//myWtime1 = MPI_Wtime();
//myWtimeInitADvars += myWtime1-myWtime0;
//myWtime0 = myWtime1;
//
//			initialHook_AD(); // Note: The method is for adoubles not for doubles
//							  //       Flow field must be initialized with initialHook() ... not with this function!!
//
//			initialHookScalarRansTurbModel1D_AD(firstCallScalarTurb); // This method must be called in order to get the adouble pointers, kine, grad_kine, kine_diff, omega, grad_omega, omega_diff, and omega
//																	  // This method also update wallDist in the ghost cells
//																	  // Also, "firstCallScalarTurb" is passed as an reference: it will be updated as "false"
//			initialHookScalarRansCombModel1D_AD(0, firstCallScalarComb); // This method must be called in order to update gamma and RoM
//
//myWtime1 = MPI_Wtime();
//myWtimeInitHooks += myWtime1-myWtime0;
//myWtime0 = myWtime1;
//	//		initializeAdjoint();
//
//			// Check possible errors: scalars = NULL
//			checkScalarsMemAD();
//
//			// Calculate Rhs
//			double (*A)[5][5] = NULL, ***AScal=NULL;
//			int flagImplicit = false;
//			myCountReducedOrder += calcResidual1D_AD(icv, rhs_rho_AD, rhs_rhou_AD, rhs_rhoE_AD, rhs_rhoScal_AD, rho_AD, rhou_AD, rhoE_AD, A, AScal, flagImplicit);
//				/* Note: calcResidual1D_AD() successively calls the following methods:
//				 *       1. calcStateVariables1D_AD()                                     -- for the 2-layer neighboring cells
//				 *       2. calcMaterialProperties1D_AD()                                 -- for the 1-layer neighboring faces
//				 *       3. setBC1D_AD()                                                  -- for the 1-layer neighboring faces
//				 *       4. if(mu_ref>0.0 || sndOrder==true): calcCv2Grad1D_AD(grad_u)    -- for the 1-layer neighboring cells
//				 *       5. if(mu_ref>0.0)                  : calcRansTurbViscMuet1D_AD() -- for the 1-layer neighboring faces (virtual function!)
//				 *       6. calcRhs1D_AD()
//				 */
//			barrierSourceNS1D_AD(icv, rhs_rho_AD, rhs_rhou_AD, rhs_rhoE_AD);
//			barrierSourceTurbScalars1D_AD(icv, rhs_rhoScal_AD, nScal, iterNewton, residNormTotOld);
//#ifdef USE_DT_OVER_VOL_SCALING // Note: Since local_dtOverVol was initialized by 1.0, you don't actually need this preprocessor
//			rhs_rho_AD *= local_dtOverVol[icv];
//			for(int i=0; i<3; ++i)
//				rhs_rhou_AD[i] *= local_dtOverVol[icv];
//			rhs_rhoE_AD *= local_dtOverVol[icv];
//			for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++) {
//				rhs_rhoScal_AD[iScal] *= local_dtOverVol[icv];
//			}
//#endif
//
//myWtime1 = MPI_Wtime();
//myWtime0 = myWtime1;
//
//			//Dependent variables
//			int icv_i = nVars*icv;
//			rhs_rho_AD     >>= rhsSingleArray[icv_i] ;
//			rhs_rhou_AD[0] >>= rhsSingleArray[icv_i+1] ;
//			rhs_rhou_AD[1] >>= rhsSingleArray[icv_i+2] ;
//			rhs_rhou_AD[2] >>= rhsSingleArray[icv_i+3] ;
//			rhs_rhoE_AD    >>= rhsSingleArray[icv_i+4] ;
//			for(int iScal=0; iScal<nScal; iScal++)
//				rhs_rhoScal_AD[iScal] >>= rhsSingleArray[icv_i+5+iScal];
//myWtime1 = MPI_Wtime();
//myWtimeFinalADvars += myWtime1-myWtime0;
//myWtime0 = myWtime1;
//
//			trace_off();
//myWtime1 = MPI_Wtime();
//myWtimeTraceOff += myWtime1-myWtime0;
//myWtime0 = myWtime1;
//
//			if(icv==0 && printStats && mpi_rank==0) {
//				cout<<"======================================================"<<endl;
//				cout<<"Tape stats for Residual evaluation (mpi_rank=0, icv=0)"<<endl;
//				cout<<"======================================================"<<endl;
//				print_tapestats(tag);
//			}
//		}
//		catch (int e) { // Catch and Re-throw
//			trace_off();
//			delete [] indvec;
//			delete [] rhs_rhoScal_AD;
//			if(NcontrolEqns==0)
//				UgpWithCvCompFlow_AD::destroy_adoubles();
//			else
//				destroy_adoubles();
//
//			throw (e);
//		}
//		catch (...) { // Catch and Re-throw
//			trace_off();
//			delete [] indvec;
//			delete [] rhs_rhoScal_AD;
//			if(NcontrolEqns==0)
//				UgpWithCvCompFlow_AD::destroy_adoubles();
//			else
//				destroy_adoubles();
//
//			throw;
//		}
//
//		if(firstCall && debugLevel>0 && mpi_rank==0) {
//			wtimeRHS = MPI_Wtime();
//			runtimeRHS += wtimeRHS - wtime0;
//		}
//
//		// +++++++++++++++++++++++++++++++++
//		// Calculate Jacobian
//		// +++++++++++++++++++++++++++++++++
//		/* coordinate format for Jacobian */
//		unsigned int *rind   = NULL;        /* row indices    */
//		unsigned int *cind   = NULL;        /* column indices */
//		double       *values = NULL;        /* values         */
//		int          nnz     = 0;
//
//		int options[8]; /* option for Jacobian computation (ADOL-C) */
//		options[0] = 0; /* sparsity pattern computation by index domains (default) */
//		options[1] = 0; /*                                     safe mode (default) */
//		options[2] = 0; /*                          not required if options[0] = 0 */
//		options[3] = 1; /*     way of compression: 0 = column compression, 1 = row */
//
//		sparse_jac(tag, nVars, ncv_gg*nVars+NcontrolEqns, 0, indvec, &nnz, &rind, &cind, &values, options);
//				/* Arguments:
//				 *   tag = tape identification
//				 *   dims1 = number of dependent variables
//				 *   dims2 = number of independent variables
//				 *   repeat = indicate repeated call at same argument
//				 *            repeat=0 : The functions are called at a point with a new sparsity structure
//				 *            repeat=1 : Re-usage of the sparsity pattern from the previous call.
//				 *   indvec[dims2] = independent vector
//				 *   nnz = number of nonzeros
//				 *   rind[nnz] = row index
//				 *   cind[nnz] = column index
//				 *   values[nnz] = non-zero values
//				 *   option[4] = array of control parameters
//				 * Note: The sparse_jac() function calls ColPack::JacobianRecovery1D::RecoverD2Row_CoordinateFormat_unmanaged()
//				 *       And rind, cind, values are allocated by the following commands in the function (see JacobianRecovery1D.cpp):
//				 *           (*ip2_RowIndex) = (unsigned int*) malloc(numOfNonZeros * sizeof(unsigned int));
//				 *           (*ip2_ColumnIndex) = (unsigned int*) malloc(numOfNonZeros * sizeof(unsigned int));
//				 *           (*dp2_JacobianValue) = (double*) malloc(numOfNonZeros * sizeof(double));
//				 *       However, the output is NOT managed by ColPack.
//				 *       Thus, the user should free the output manually using free() (NOT delete) function when it is no longer needed.
//				 */
//		if(icv==0 && printStats && mpi_rank==0) {
//			int icv_i = nVars*icv;
//
//			cout<<"================================================"<<endl;
//			cout<<"Jacobian stats (mpi_rank=0, icv=0)"<<endl;
//			cout<<"================================================"<<endl;
//			cout<<"(Rhs values =";
//			for(int ii=0; ii<5+nScal; ++ii) printf("  %.4e", rhsSingleArray[icv_i+ii]);
//			cout<<")"<<endl;
//			cout<<"Number of independents = "<<ncv_gg*nVars<<endl;
//			cout<<"Number of dependents   = "<<nVars<<endl;
//			cout<<"Number of non-zeros    = "<<nnz<<endl;
//			cout<<"================================================"<<endl<<endl;
//		}
//		int firstRindReal = icv*nVars;
//		jacMatrixSTL.blockCopy(nnz, firstRindReal, rind, cind, values);
//
//		if(firstCall && debugLevel>0 && mpi_rank==0) {
//			wtimeJac = MPI_Wtime();
//			runtimeJac += wtimeJac - wtimeRHS;
//		}
//
//		// +++++++++++++++++++++++++++++++++
//		// Free memory
//		// +++++++++++++++++++++++++++++++++
//		free(rind); 	//rind = NULL;
//		free(cind); 	//cind = NULL;
//		free(values); 	//values = NULL;
//
//		delete [] rhs_rhoScal_AD;
//
//		if(mpi_rank==0 && icv%50==0) {
//			FILE *fp = fopen(JAC1D_STATUS_FILENAME, "a");
//			fprintf(fp, "%4d  %12d  %10d   %10.3f\n", step, iterNewton, icv, MPI_Wtime()-initWtime);
//			fclose(fp);
//		}
//	}
//	jacMatrixSTL.finalizeMat();
//
//	// free indvec
//	delete [] indvec;
//
//	// free memory allocated for flow variables
//	if(NcontrolEqns==0)
//		UgpWithCvCompFlow_AD::destroy_adoubles();
//	else
//		destroy_adoubles();
//
//	if(firstCall && debugLevel>0 && mpi_rank==0) {
//		cout<<"   > calcJacobian1DAD(): Runtime for the RHS calculations [sec] = "<<runtimeRHS<<endl;
//		printf("                           Runtime Details [sec]\n");
//		printf("                           1. Allocate AD vars         =%9.4f\n", myWtimeAllocADvars);
//		printf("                           2. Initialize AD vars       =%9.4f\n", myWtimeInitADvars);
//		printf("                           3. Trace-on                 =%9.4f\n", myWtimeTraceOn);
//		printf("                           4. Assign indep vars        =%9.4f\n", myWtimeAssgnIndep);
//		printf("                           5. InitialHook functions    =%9.4f\n", myWtimeInitHooks);
//		printf("                           -------------------------------------\n");
//		printf("                           6. Find neighbor faces      =%9.4f\n", myWtimeFindFaces);
//		printf("                           7. Calc state vars          =%9.4f\n", myWtimeStateVars);
//		printf("                           8. Calc material properties =%9.4f\n", myWtimeMaterPropert);
//		printf("                           9. Set BC                   =%9.4f\n", myWtimeSetBC);
//		printf("                           10. Calc gradients          =%9.4f\n", myWtimeCalcGrad);
//		printf("                           11. Set up zone arrays      =%9.4f\n", myWtimeSetupZoneArray);
//		printf("                           12. Calc turb visc          =%9.4f\n", myWtimeCalcTurbVisc);
//		printf("                           13. Calc RHS                =%9.4f\n", myWtimeCalcRhs);
//		printf("                           -------------------------------------\n");
//		printf("                           14. Finalize AD vars        =%9.4f\n", myWtimeFinalADvars);
//		printf("                           15. Trace-off               =%9.4f\n\n", myWtimeTraceOff);
//
//		cout<<"   > calcJacobian1DAD(): Runtime for the Jac evaluations  [sec] = "<<runtimeJac<<endl;
//	}
//
//	//
//	firstCall = false;
//
//	MPI_Allreduce(&myCountReducedOrder, &CountReducedOrder, 1, MPI_INT, MPI_SUM, mpi_comm);
//	return CountReducedOrder>0;
//}

/*
 * Method: calcJacobianAD
 * ----------------------
 * calculate Rhs and get the Jacobian matrix using AD
 * Original code of Rhs part = JoeWithModels_AD::calcResidualDerivative(double***, double***, int)
 * Arguments: jacMatrix      = The Jacobian matrix (NOTE: the size of the matrix is (ncv*nVars)*(ncv_gg*nVars), where nVars=5+nScal)
 *            rhsSingleArray = R.H.S. of the N-S system with scalar equations (NOTE: The size is ncv_gg*(5+nScal) because using Dr. Duraisamy's adjoint code)
 *            printStats     = If yes, show the Jacobian statistics at the first cpu core on the screen (for a debugging purpose)
 * Return: If backtracking line search is required (due to negative density, etc.), return true.
 *
 * Note: For some reason, calcResidual_AD() calculates for ncv_gg CVs, which finally yields ncv_gg rows in Jacobian instead of ncv rows
 *       (For details of calcResidual_AD(), contact Dr. Duraisamy)
 */
int IkeWithPsALC_AD::calcJacobianAD(MatComprsed &jacMatrix, double *rhsSingleArray, const int debugLevel, const int NcontrolEqns)
{
	assert(lambda_AD == NULL);

	int nScal = scalarTranspEqVector.size();
	int nVars = 5+nScal;

	int countReducedCV;

	static bool firstCall = true;

	// +++++++++++++++++++++++++++++++++
	// Allocate memory for Flow residual
	// +++++++++++++++++++++++++++++++++
	REALA *rhs_rho_AD       = new REALA[ncv_gg];
	REALA (*rhs_rhou_AD)[3] = new REALA[ncv_gg][3];
	REALA *rhs_rhoE_AD      = new REALA[ncv_gg];

	REALAS **rhs_rhoScal_AD = NULL;
	if (nScal > 0)
		getMem2D(&rhs_rhoScal_AD, 0, nScal-1, 0, ncv_gg-1, "rhs_rhoScal_AD");

	double *indvec = new double [nVars*ncv_gg+NcontrolEqns];

	// +++++++++++++++++++++++++++++++++
	// Start the calculation of Residual wrt flow
	// +++++++++++++++++++++++++++++++++
	int tag = mpi_rank;
	trace_on(tag, 0); // Note: trace_on(tag, keep) -- If keep==1, the numerical values of all active variables
                      //       are recorded in a buffered temporary file before they will be overwritten.
                      //       (preparing the scene for an immediately following reverse mode differentiation)

	// allocate memory
	// allocate memory for flow variables: e.g. rho_AD, rhou_AD, vel, mul_fa, ...
	if(NcontrolEqns==0)
		UgpWithCvCompFlow_AD::initialize_adoubles();
	else {
		initialize_adoubles(); // The pseudo-arclength continuation method requires additional AD variables for the control eqns
		assert(lambda_AD != NULL);
	}
//	// initialize the flow field for the base variables
//	initialHook();

	// write the tecplot file
	if (firstCall && initial_flowfield_output == "YES")
		writeData(0);

	// Independent variables
	for (int icv = 0; icv < ncv_gg; icv++) {
		rho_AD[icv]      <<= rho[icv] ;
		rhou_AD[icv][0]  <<= rhou[icv][0] ;
		rhou_AD[icv][1]  <<= rhou[icv][1] ;
		rhou_AD[icv][2]  <<= rhou[icv][2] ;
		rhoE_AD[icv]     <<= rhoE[icv] ;
		for(int iScal=0; iScal<nScal; iScal++)
			scalarTranspEqVector_AD[iScal].phi[icv] <<= scalarTranspEqVector[iScal].phi[icv];
	}
	for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
		lambda_AD[iEqn] <<= lambda[iEqn];

	initialHook_AD(); // Note: The flow field should be initialized with initialHook()
	                  // To DO: in the original code initialHook_AD() isn't called here ( see JoeWithModels_AD::calcResidualDerivative() )
	initialHookScalarRansTurbModel_AD();
	initialHookScalarRansCombModel_AD(0);

	// Check possible errors: scalars = NULL
	checkScalarsMemAD();

	// Calculate Rhs
	double (*A)[5][5] = NULL, ***AScal=NULL;
	int flagImplicit = false;
	countReducedCV = calcResidual_AD(rhs_rho_AD, rhs_rhou_AD, rhs_rhoE_AD, rhs_rhoScal_AD, rho_AD, rhou_AD, rhoE_AD, A, AScal, flagImplicit);
          // Note: calcaResidual_AD() successively calls the following methods:
          //       1. calcStateVariables_AD(),   2. calcMaterialProperties_AD(),  
          //       3. setBC_AD(),
          //       if(mu_ref>0.0 || sndOrder==true): 4. calcCv2Grad_AD(),
          //       if(mu_ref>0.0):                   5. calcRansTurbViscMuet_AD(),
          //       6. calcRhs_AD()
	barrierSourceNS_AD(rhs_rho_AD, rhs_rhou_AD, rhs_rhoE_AD);
	barrierSourceTurbScalars_AD(rhs_rhoScal_AD, nScal, iterNewton, residNormTotOld);

	//Dependent variables
	for (int icv = 0; icv < ncv_gg; icv++) {
		int icv_i = nVars*icv;
		rhs_rho_AD[icv]      >>= rhsSingleArray[icv_i] ;
		rhs_rhou_AD[icv][0]  >>= rhsSingleArray[icv_i+1] ;
		rhs_rhou_AD[icv][1]  >>= rhsSingleArray[icv_i+2] ;
		rhs_rhou_AD[icv][2]  >>= rhsSingleArray[icv_i+3] ;
		rhs_rhoE_AD[icv]     >>= rhsSingleArray[icv_i+4] ;
		for(int iScal=0; iScal<nScal; iScal++) {
			rhs_rhoScal_AD[iScal][icv] >>= rhsSingleArray[icv_i+5+iScal];
		}
	}
	trace_off(); // Note: trace_off(file) -- If the argument "file" is omitted, it defaults to 0,
                 //                          so that the tape array is written onto an external file
                 //                          only if the length of any of the buffers exceeds BUFSIZE

	if(mpi_rank==0 && (debugLevel>1 || (debugLevel>0 && firstCall))) {
		cout<<"=================================="<<endl;
		cout<<"Tape stats for Residual evaluation"<<endl;
		cout<<"=================================="<<endl;
		print_tapestats(tag);
	}
	int memoryDeficit = countMemoryDeficit_fromTapeStats(tag);
	if(memoryDeficit > 0)
		cout<<"WARNING! Too much memory is requires for AD at mpi_rank=="<<mpi_rank<<": memory deficit = "<<memoryDeficit<<endl;

	// +++++++++++++++++++++++++++++++++
	// Calculate Jacobian
	// +++++++++++++++++++++++++++++++++
	if(mpi_rank==0 && (debugLevel>1 || (debugLevel>0 && firstCall)))
		printf("\nSparse Jacobians:\n");

	for (int icv = 0; icv < ncv_gg; icv++) {
		int icv_i = nVars*icv;
		indvec[icv_i]   = rho[icv] ;
		indvec[icv_i+1] = rhou[icv][0];
		indvec[icv_i+2] = rhou[icv][1];
		indvec[icv_i+3] = rhou[icv][2];
		indvec[icv_i+4] = rhoE[icv];
		for(int iScal=0; iScal<nScal; iScal++)
			indvec[icv_i+5+iScal] = scalarTranspEqVector[iScal].phi[icv];
	}
	for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
		indvec[nVars*ncv_gg+iEqn] = lambda[iEqn];

	/* coordinate format for Jacobian */
	unsigned int *rind  = NULL;        /* row indices    */
	unsigned int *cind  = NULL;        /* column indices */
	double       *values = NULL;       /* values         */
	int nnz;

	int options[8]; /* option for Jacobian computation (ADOL-C) */
	options[0] = 0; /* sparsity pattern computation by index domains (default) */
	options[1] = 0; /*                                     safe mode (default) */
	options[2] = 0; /*                          not required if options[0] = 0 */
	options[3] = 1; /*     way of compression: 0 = column compression, 1 = row */

	sparse_jac(tag, ncv_gg*nVars, ncv_gg*nVars+NcontrolEqns, 0, indvec, &nnz, &rind, &cind, &values, options);
			/* Arguments:
			 *   tag = tape identification
			 *   dims1 = number of dependent variables
			 *   dims2 = number of independent variables
			 *   repeat = indicate repeated call at same argument
			 *            repeat=0 : The functions are called at a point with a new sparsity structure
			 *            repeat=1 : Re-usage of the sparsity pattern from the previous call.
			 *   indvec[dims2] = independent vector
			 *   nnz = number of nonzeros
			 *   rind[nnz] = row index
			 *   cind[nnz] = column index
			 *   values[nnz] = non-zero values
			 *   option[4] = array of control parameters
			 * Note: The sparse_jac() function calls ColPack::JacobianRecovery1D::RecoverD2Row_CoordinateFormat_unmanaged()
			 *       And rind, cind, values are allocated by the following commands in the function (see JacobianRecovery1D.cpp):
			 *           (*ip2_RowIndex) = (unsigned int*) malloc(numOfNonZeros * sizeof(unsigned int));
			 *           (*ip2_ColumnIndex) = (unsigned int*) malloc(numOfNonZeros * sizeof(unsigned int));
			 *           (*dp2_JacobianValue) = (double*) malloc(numOfNonZeros * sizeof(double));
			 *       However, the output is NOT managed by ColPack.
			 *       Thus, the user should free the output manually using free() (NOT delete) function when it is no longer needed.
			 */

	jacMatrix.setMatSize(ncv_gg*nVars+NcontrolEqns, ncv*nVars); // NOTE: the number of rows is ncv*nVars
	jacMatrix.matCopy(nnz, rind, cind, values);    // NOTE: jacMatrix stores only 0 ~ (ncv*nVars) rows

	// +++++++++++++++++++++++++++++++++
	// Free memory
	// +++++++++++++++++++++++++++++++++
	free(rind);
	free(cind);
	free(values);
	removeTape(tag, ADOLC_REMOVE_COMPLETELY);  // Note: "int removeTape(short tapeID, short type)" in tape_handling.cpp
	                                           //       If type==ADOLC_REMOVE_COMPLETELY, remove the tape instead of free it.
	                                           //       Returns 0 if the removing process is completed okay.
	destroy_adoubles();

//delete [] rhs_rhoScal_AD ;
	freeMem2D(rhs_rhoScal_AD, 0, nScal-1, 0, ncv_gg-1);
	delete [] rhs_rhoE_AD;
	delete [] rhs_rhou_AD;
	delete [] rhs_rho_AD;

	delete [] indvec;

	//
	firstCall = false;

	return countReducedCV;
}

/*
 * Method: calcJacobian1DAD_calcRhs
 * --------------------------------
 * Calculate RHS for Jacobian calculation
 * Return:
 *   rhsSingleArray
 *   myCountReducedOrder
 *   wall-clock times for RHS calculation
 * Note: "firstCallScalarTurb" will be modified in this method (when initialHookScalarRansTurbModel1D_AD() is called)
 */
void IkeWithPsALC_AD::calcJacobian1DAD_calcRhs(double *rhsSingleArray, int &myCountReducedOrder, wallTimeJacCalc &myWTimeJacCalc,
		const int icv, const int tag, const int nScal, const int debugLevel, const int NcontrolEqns,
		const bool firstCall, bool &firstCallScalarTurb, bool &firstCallScalarComb) {
	double myWtime0 = MPI_Wtime();
	double myWtime1;

	int nVars = 5+nScal;

	// +++++++++++++++++++++++++++++++++
	// Allocate memory for Flow residual
	// +++++++++++++++++++++++++++++++++
#ifdef USE_MEM_SAVING_ADVAR_1D_
	assert(nbocv2_eachIcv.empty()); 	find2layerCSRstruct_eachIcv(icv, nbocv2_eachIcv, false);
	assert(nbocv2ff_eachIcv.empty()); 	find2layerCSRstruct_eachIcv(icv, nbocv2ff_eachIcv, true);
	assert(fa2_eachIcv.empty()); 		find2LayerFaceIndices(icv, fa2_eachIcv);

	initialize_adoubles(icv, NcontrolEqns);
	// The pseudo-arclength continuation method requires additional AD variables for the control eqns
	if(NcontrolEqns>0)
		assert(lambda_AD != NULL);
#endif

#ifdef USE_MEM_SAVING_ADVAR_1D_
	REALA *rhs_rho_AD  = new REALA;
	REALA *rhs_rhou_AD = new REALA[3];
	REALA *rhs_rhoE_AD = new REALA;
#else
	REALA rhs_rho_AD;
	REALA rhs_rhou_AD[3];
	REALA rhs_rhoE_AD;
#endif
	REALAS *rhs_rhoScal_AD = new REALA [nScal];

#ifdef USE_MEM_SAVING_ADVAR_1D_
	(*rhs_rho_AD) = 0.0;
	for(int i=0; i<3; ++i)
		rhs_rhou_AD[i] = 0.0;
	(*rhs_rhoE_AD) = 0.0;
	for(int i=0; i<nScal; ++i)
		rhs_rhoScal_AD[i] = 0.0;
#else
	rhs_rho_AD = 0.0;
	for(int i=0; i<3; ++i)
		rhs_rhou_AD[i] = 0.0;
	rhs_rhoE_AD = 0.0;
	for(int i=0; i<nScal; ++i)
		rhs_rhoScal_AD[i] = 0.0;
#endif

	// +++++++++++++++++++++++++++++++++
	// Start the residual calculation
	// +++++++++++++++++++++++++++++++++
	try {
		trace_on(tag, 0); // Note: trace_on(tag, keep) -- If keep==1, the numerical values of all active variables
		                  //       are recorded in a buffered temporary file before they will be overwritten.
		                  //       (preparing the scene for an immediately following reverse mode differentiation)

		// Independent variables
#ifdef USE_MEM_SAVING_ADVAR_1D_
		for (size_t i=0; i<nbocv2_eachIcv.size(); ++i) {
			int icvTemp = nbocv2_eachIcv[i];
			rho_AD[icvTemp]      <<= rho[icvTemp] ;
			rhou_AD[icvTemp][0]  <<= rhou[icvTemp][0] ;
			rhou_AD[icvTemp][1]  <<= rhou[icvTemp][1] ;
			rhou_AD[icvTemp][2]  <<= rhou[icvTemp][2] ;
			rhoE_AD[icvTemp]     <<= rhoE[icvTemp] ;
			for(int iScal=0; iScal<nScal; iScal++)
				scalarTranspEqVector_AD[iScal].phi[icvTemp] <<= scalarTranspEqVector[iScal].phi[icvTemp];
		}
#else
		for (int icvTemp = 0; icvTemp < ncv_gg; icvTemp++) {
			rho_AD[icvTemp]      <<= rho[icvTemp] ;
			rhou_AD[icvTemp][0]  <<= rhou[icvTemp][0] ;
			rhou_AD[icvTemp][1]  <<= rhou[icvTemp][1] ;
			rhou_AD[icvTemp][2]  <<= rhou[icvTemp][2] ;
			rhoE_AD[icvTemp]     <<= rhoE[icvTemp] ;
			for(int iScal=0; iScal<nScal; iScal++)
				scalarTranspEqVector_AD[iScal].phi[icvTemp] <<= scalarTranspEqVector[iScal].phi[icvTemp];
		}
#endif
		if(NcontrolEqns > 0) {
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
				lambda_AD[iEqn] <<= lambda[iEqn];
		}

		// Check wall-clock time
		myWtime1 = MPI_Wtime();
		myWTimeJacCalc.wallTimeRHSinit += myWtime1-myWtime0;
		myWtime0 = myWtime1;

		initialHook_AD(); // Note: The method is for adoubles not for doubles
						  //       Flow field must be initialized with initialHook() ... never with this function!!

#ifdef USE_MEM_SAVING_ADVAR_1D_
		initialHookScalarRansTurbModel1D_AD(nbocv2ff_eachIcv, firstCallScalarTurb); // This method must be called in order to get the adouble pointers, kine, grad_kine, kine_diff, omega, grad_omega, omega_diff, and omega

		initialHookScalarRansCombModel1D_AD(0, nbocv2ff_eachIcv, firstCallScalarComb); // This method must be called in order to update gamma and RoM
#else
		initialHookScalarRansTurbModel1D_AD(firstCallScalarTurb); // This method must be called in order to get the adouble pointers, kine, grad_kine, kine_diff, omega, grad_omega, omega_diff, and omega
																  // This method also update wallDist in the ghost cells
																  // Also, "firstCallScalarTurb" is passed as an reference: it will be updated as "false"
		initialHookScalarRansCombModel1D_AD(0, firstCallScalarComb); // This method must be called in order to update gamma and RoM
#endif
//			initializeAdjoint();

		// Check possible errors: scalars = NULL
		checkScalarsMemAD();

		// Check wall-clock time
		myWtime1 = MPI_Wtime();
		myWTimeJacCalc.wallTimeRHSinit = myWtime1-myWtime0;
		myWtime0 = myWtime1;

		// Calculate Rhs
		double (*A)[5][5] = NULL, ***AScal=NULL;
		int flagImplicit = false;
#ifdef USE_MEM_SAVING_ADVAR_1D_
		myCountReducedOrder += calcResidual1D_AD(icv, *rhs_rho_AD, rhs_rhou_AD, *rhs_rhoE_AD, rhs_rhoScal_AD, rho_AD, rhou_AD, rhoE_AD, A, AScal, flagImplicit);
			/* Note: calcResidual1D_AD() successively calls the following methods:
			 *       1. calcStateVariables1D_AD()                                     -- for the 2-layer neighboring cells
			 *       2. calcMaterialProperties1D_AD()                                 -- for the 1-layer neighboring faces
			 *       3. setBC1D_AD()                                                  -- for the 1-layer neighboring faces
			 *       4. if(mu_ref>0.0 || sndOrder==true): calcCv2Grad1D_AD(grad_u)    -- for the 1-layer neighboring cells
			 *       5. if(mu_ref>0.0)                  : calcRansTurbViscMuet1D_AD() -- for the 1-layer neighboring faces (virtual function!)
			 *       6. calcRhs1D_AD()
			 */
		barrierSourceNS1D_AD(icv, *rhs_rho_AD, rhs_rhou_AD, *rhs_rhoE_AD);
		barrierSourceTurbScalars1D_AD(icv, rhs_rhoScal_AD, nScal, iterNewton, residNormTotOld);
#else
		myCountReducedOrder += calcResidual1D_AD(icv, rhs_rho_AD, rhs_rhou_AD, rhs_rhoE_AD, rhs_rhoScal_AD, rho_AD, rhou_AD, rhoE_AD, A, AScal, flagImplicit);
			/* Note: calcResidual1D_AD() successively calls the following methods:
			 *       1. calcStateVariables1D_AD()                                     -- for the 2-layer neighboring cells
			 *       2. calcMaterialProperties1D_AD()                                 -- for the 1-layer neighboring faces
			 *       3. setBC1D_AD()                                                  -- for the 1-layer neighboring faces
			 *       4. if(mu_ref>0.0 || sndOrder==true): calcCv2Grad1D_AD(grad_u)    -- for the 1-layer neighboring cells
			 *       5. if(mu_ref>0.0)                  : calcRansTurbViscMuet1D_AD() -- for the 1-layer neighboring faces (virtual function!)
			 *       6. calcRhs1D_AD()
			 */
		barrierSourceNS1D_AD(icv, rhs_rho_AD, rhs_rhou_AD, rhs_rhoE_AD);
		barrierSourceTurbScalars1D_AD(icv, rhs_rhoScal_AD, nScal, iterNewton, residNormTotOld);
#endif
#ifdef USE_DT_OVER_VOL_SCALING // Note: Since local_dtOverVol was initialized by 1.0, you don't actually need this preprocessor
		rhs_rho_AD *= local_dtOverVol[icv];
		for(int i=0; i<3; ++i)
			rhs_rhou_AD[i] *= local_dtOverVol[icv];
		rhs_rhoE_AD *= local_dtOverVol[icv];
		for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++) {
			rhs_rhoScal_AD[iScal] *= local_dtOverVol[icv];
		}
#endif

		// Calculate weight for RHS
		assert(weightRhs != NULL);
		int startingIndex = icv*nVars;
		for(int i=0; i<nVars; ++i)
			weightRhs[startingIndex+i] = 1.0;

		calcWeightRhs(icv, weightRhsMethod, nScal);
		calcWeightRhsHook(icv, nScal);

#ifdef USE_MEM_SAVING_ADVAR_1D_
		*rhs_rho_AD *= weightRhs[startingIndex];

		for(int i=0; i<3; ++i)
			rhs_rhou_AD[i] *= weightRhs[startingIndex+1+i];

		*rhs_rhoE_AD *= weightRhs[startingIndex+4];

		for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
			rhs_rhoScal_AD[iScal] *= weightRhs[startingIndex+5+iScal];
#else
		rhs_rho_AD *= weightRhs[startingIndex];

		for(int i=0; i<3; ++i)
			rhs_rhou_AD[i] *= weightRhs[startingIndex+1+i];

		rhs_rhoE_AD *= weightRhs[startingIndex+4];

		for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
			rhs_rhoScal_AD[iScal] *= weightRhs[startingIndex+5+iScal];
#endif

		// Check wall-clock time
		myWtime1 = MPI_Wtime();
		myWTimeJacCalc.wallTimeRHScalc += myWtime1-myWtime0;
		myWtime0 = myWtime1;

		//Dependent variables
		int icv_i = nVars*icv;
#ifdef USE_MEM_SAVING_ADVAR_1D_
		(*rhs_rho_AD)  >>= rhsSingleArray[icv_i] ;
		rhs_rhou_AD[0] >>= rhsSingleArray[icv_i+1] ;
		rhs_rhou_AD[1] >>= rhsSingleArray[icv_i+2] ;
		rhs_rhou_AD[2] >>= rhsSingleArray[icv_i+3] ;
		(*rhs_rhoE_AD) >>= rhsSingleArray[icv_i+4] ;
		for(int iScal=0; iScal<nScal; iScal++)
			rhs_rhoScal_AD[iScal] >>= rhsSingleArray[icv_i+5+iScal];
#else
		rhs_rho_AD     >>= rhsSingleArray[icv_i] ;
		rhs_rhou_AD[0] >>= rhsSingleArray[icv_i+1] ;
		rhs_rhou_AD[1] >>= rhsSingleArray[icv_i+2] ;
		rhs_rhou_AD[2] >>= rhsSingleArray[icv_i+3] ;
		rhs_rhoE_AD    >>= rhsSingleArray[icv_i+4] ;
		for(int iScal=0; iScal<nScal; iScal++)
			rhs_rhoScal_AD[iScal] >>= rhsSingleArray[icv_i+5+iScal];
#endif
		trace_off(); // Note: trace_off(file) -- If the argument "file" is omitted, it defaults to 0,
		             //                          so that the tape array is written onto an external file
		             //                          only if the length of any of the buffers exceeds BUFSIZE

		if(icv==0 && mpi_rank==0 && (debugLevel>1 || (debugLevel>0 && firstCall))) {
			cout<<"======================================================"<<endl;
			cout<<"Tape stats for Residual evaluation (mpi_rank=0, icv=0)"<<endl;
			cout<<"======================================================"<<endl;
			print_tapestats(tag);
		}
		int memoryDeficit = countMemoryDeficit_fromTapeStats(tag);
		if(memoryDeficit > 0)
			cout<<"WARNING! Too much memory is requires for AD at mpi_rank=="<<mpi_rank<<", icv="<<icv<<": memory deficit = "<<memoryDeficit<<endl;
	}
	catch (int e) { // Catch and Re-throw
		trace_off();
		delete [] rhs_rhoScal_AD;
		if(NcontrolEqns==0)
			UgpWithCvCompFlow_AD::destroy_adoubles();
		else
			destroy_adoubles();

		throw (e);
	}
	catch (...) { // Catch and Re-throw
		trace_off();
		delete [] rhs_rhoScal_AD;
		if(NcontrolEqns==0)
			UgpWithCvCompFlow_AD::destroy_adoubles();
		else
			destroy_adoubles();

		throw;
	}

#ifdef USE_MEM_SAVING_ADVAR_1D_
	delete rhs_rho_AD;
	delete [] rhs_rhou_AD;
	delete rhs_rhoE_AD;
#endif
	delete [] rhs_rhoScal_AD;

	// Check wall-clock time
	myWtime1 = MPI_Wtime();
	myWTimeJacCalc.wallTimeRHSfinal += myWtime1-myWtime0;
}


/*
 * Method: probeInCalcJacobian1DAD
 * -------------------------------
 * This method will be called in the end of calcJacobian1DAD() -- just before the memory-releasing step.
 * If the user wants to probe a quantity (usually something related to the RHS calculation) at every Newton-Raphson iteration,
 * she/he may want to call this method.
 */
void IkeWithPsALC_AD::probeInCalcJacobian1DAD(MatComprsedSTL &jacMatrixSTL, double *rhsSingleArray, const int debugLevel, const int NcontrolEqns) {
	/* empty */
}

/*
 * Method: solveLinSysNSCoupled2
 * -----------------------------
 * Solve coupled linear system for Navier-Stokes equations and scalars together
 * original code: solveCoupledLinSysNSCoupled() in UgpWithCvCompFlow.h
 *
 * Note: If useOldSolnForInitGuess is TRUE, phi will be used as the initial guess.
 */
template <class MatT>
void IkeWithPsALC_AD::solveLinSysNSCoupled2(double *phi, MatT &A, double *rhs,
		bool useOldSolnForInitGuess, const double zeroAbs, const double zeroRel, const int maxIter, int nScal, const int monitorConvergInterval,
		int NcontrolParams, int ncv_gg, double **vecC, const double *d, const int step, const int newtonIter) {
	switch (linearSolverNS) { // Note: BCGSTAB, PETSC_GMRES, and PETSC_LU_MUMPS are defined in UgpWithCv2.h, whereas linearSolverNS is set in UgpWithCvCompFlow.h
	case PETSC_GMRES:
		// on the first time, initialize the petsc solver...
		if (petscSolver2 == NULL) {
			if (nbocv_v_global == NULL) {
				nbocv_v_global = new int[ncv_gg];
				for (int icv = 0; icv < ncv; icv++)
					nbocv_v_global[icv] = cvora[mpi_rank] + icv;
				updateCvDataG1G2(nbocv_v_global, REPLACE_DATA);
			}

			string pcType; // "BJACOBI", "ASM", "ILU", "AMG", "NONE"
			bool pcReuse;
			int pcLevels;  // Note: this is required only for the ILU pre-conditioner. Otherwise, any number doesn't matter.
			getPCcontext(pcType, pcReuse, pcLevels);

			petscSolver2 = new PetscSolver2(cvora, nbocv2_i, nbocv2_v, 5 + nScal, pcType, pcLevels, pcReuse, useOldSolnForInitGuess, NcontrolParams, monitorConvergInterval);
			petscSolver2->setTresholds(zeroAbs, zeroRel, maxIter);

			petscSolver2->setGmresRestart(gmresRestart);
		}
		if(useOldSolnForInitGuess)
			petscSolver2->solveGMRES<MatT>(A, phi, rhs, phi,  cvora, nbocv_v_global, 5 + nScal, NcontrolParams, ncv_gg, vecC, d, step, newtonIter);
		else
			petscSolver2->solveGMRES<MatT>(A, phi, rhs, NULL, cvora, nbocv_v_global, 5 + nScal, NcontrolParams, ncv_gg, vecC, d, step, newtonIter);

		break;

#if defined(PETSC_HAVE_MUMPS) // This flag, PETSC_HAVE_MUMPS, is defined if MUMPS is installed while configuring PETSc.
                              // Note: For superLU -- PETSC_HAVE_SUPERLU or PETSC_HAVE_SUPERLU_DIST.
	case PETSC_LU_MUMPS:
		// on the first time, initialize the petsc solver...
		if (petscSolver2 == NULL) {
			if (nbocv_v_global == NULL) {
				nbocv_v_global = new int[ncv_gg];
				for (int icv = 0; icv < ncv; icv++)
					nbocv_v_global[icv] = cvora[mpi_rank] + icv;
				updateCvDataG1G2(nbocv_v_global, REPLACE_DATA);
			}

			int printLevel = 1; // MUMPS print level (0-4): MUMPS default = 2
			bool pcReuseDummy = false;
			petscSolver2 = new PetscSolver2(cvora, nbocv2_i, nbocv2_v, 5 + nScal, "LU", printLevel, pcReuseDummy, false, NcontrolParams);
//			petscSolver2->setTresholds(zeroAbs, zeroRel, maxIter);
		}
		petscSolver2->solveLUMUMPS<MatT>(A, phi, rhs, cvora, nbocv_v_global, 5 + nScal, NcontrolParams, ncv_gg, vecC, d, step, newtonIter);

		break;
#endif
	default:
		if (mpi_rank == 0)
			cerr << "Error: unrecognized solver: " << linearSolverNS
					<< endl;
		throw(-1);
		break;
	}
}

template <class MatT>
void IkeWithPsALC_AD::solveLinSysNSCoupled2(double *phi, MatT &A, double *rhs, int &nIter, double &absResid,
		bool useOldSolnForInitGuess, const double zeroAbs, const double zeroRel, const int maxIter, int nScal, const int monitorConvergInterval,
		int NcontrolParams, int ncv_gg, double **vecC, const double *d, const int step, const int newtonIter) {
	switch (linearSolverNS) { // Note: BCGSTAB, PETSC_GMRES, and PETSC_LU_MUMPS are defined in UgpWithCv2.h, whereas linearSolverNS is set in UgpWithCvCompFlow.h
	case PETSC_GMRES:
		// on the first time, instantiate the petsc solver...
		if (petscSolver2 == NULL) {
			if (nbocv_v_global == NULL) {
				nbocv_v_global = new int[ncv_gg];
				for (int icv = 0; icv < ncv; icv++)
					nbocv_v_global[icv] = cvora[mpi_rank] + icv;
				updateCvDataG1G2(nbocv_v_global, REPLACE_DATA);
			}

			string pcType; // "BJACOBI", "ASM", "ILU", "AMG", "NONE"
			bool pcReuse;
			int pcLevels;  // Note: this is required only for the ILU pre-conditioner. Otherwise, any number doesn't matter.
			getPCcontext(pcType, pcReuse, pcLevels);

	#ifdef USE_SLEPC_WITH_PETSC
			petscSolver2 = new SlepcSolver2(cvora, nbocv2_i, nbocv2_v, 5 + nScal, pcType, pcLevels, pcReuse, useOldSolnForInitGuess, NcontrolParams, monitorConvergInterval);
	#else
			petscSolver2 = new PetscSolver2(cvora, nbocv2_i, nbocv2_v, 5 + nScal, pcType, pcLevels, pcReuse, useOldSolnForInitGuess, NcontrolParams, monitorConvergInterval);
	#endif

			petscSolver2->setTresholds(zeroAbs, zeroRel, maxIter);

			int gmresRestart = getIntParam("GMRES_RESTART", "100");
			petscSolver2->setGmresRestart(gmresRestart);
		}
		if(useOldSolnForInitGuess)
			petscSolver2->solveGMRES<MatT>(nIter, absResid, A, phi, rhs, phi,  cvora, nbocv_v_global, 5 + nScal, NcontrolParams, ncv_gg, vecC, d, step, newtonIter);
		else
			petscSolver2->solveGMRES<MatT>(nIter, absResid, A, phi, rhs, NULL, cvora, nbocv_v_global, 5 + nScal, NcontrolParams, ncv_gg, vecC, d, step, newtonIter);

		break;

#if defined(PETSC_HAVE_MUMPS) // This flag, PETSC_HAVE_MUMPS, is defined if MUMPS is installed while configuring PETSc.
                              // Note: For superLU -- PETSC_HAVE_SUPERLU or PETSC_HAVE_SUPERLU_DIST.
	case PETSC_LU_MUMPS:
		// on the first time, initialize the petsc solver...
		if (petscSolver2 == NULL) {
			if (nbocv_v_global == NULL) {
				nbocv_v_global = new int[ncv_gg];
				for (int icv = 0; icv < ncv; icv++)
					nbocv_v_global[icv] = cvora[mpi_rank] + icv;
				updateCvDataG1G2(nbocv_v_global, REPLACE_DATA);
			}

			int printLevel = 1; // MUMPS print level (0-4): MUMPS default = 2
			bool pcReuseDummy = false;
			petscSolver2 = new PetscSolver2(cvora, nbocv2_i, nbocv2_v, 5 + nScal, "LU", printLevel, pcReuseDummy, false, NcontrolParams);
//			petscSolver2->setTresholds(zeroAbs, zeroRel, maxIter);
		}
		petscSolver2->solveLUMUMPS<MatT>(A, phi, rhs, cvora, nbocv_v_global, 5 + nScal, NcontrolParams, ncv_gg, vecC, d, step, newtonIter);

		break;
#endif

	default:
		if (mpi_rank == 0)
			cerr << "Error: unrecognized solver: " << linearSolverNS
					<< endl;
		throw(-1);
		break;
	}
}

template <class MatT>
void IkeWithPsALC_AD::solveLinSysNSCoupled2(double *phi, MatT &A, double *rhs, int &nIter, double &absResid, vector<pair<int, double> >& kspMonitorHistory,
		bool useOldSolnForInitGuess, const double zeroAbs, const double zeroRel, const int maxIter, int nScal, const int monitorConvergInterval,
		int NcontrolParams, int ncv_gg, double **vecC, const double *d, const int step, const int newtonIter) {
	switch (linearSolverNS) { // Note: BCGSTAB, PETSC_GMRES, and PETSC_LU_MUMPS are defined in UgpWithCv2.h, whereas linearSolverNS is set in UgpWithCvCompFlow.h
	case PETSC_GMRES:
		// on the first time, instantiate the petsc solver...
		if (petscSolver2 == NULL) {
			if (nbocv_v_global == NULL) {
				nbocv_v_global = new int[ncv_gg];
				for (int icv = 0; icv < ncv; icv++)
					nbocv_v_global[icv] = cvora[mpi_rank] + icv;
				updateCvDataG1G2(nbocv_v_global, REPLACE_DATA);
			}

			string pcType; // "BJACOBI", "ASM", "ILU", "AMG", "NONE"
			bool pcReuse;
			int pcLevels;  // Note: this is required only for the ILU pre-conditioner. Otherwise, any number doesn't matter.
			getPCcontext(pcType, pcReuse, pcLevels);

//	#ifdef USE_SLEPC_WITH_PETSC
//			petscSolver2 = new SlepcSolver2(cvora, nbocv2_i, nbocv2_v, 5 + nScal, pcType, pcLevels, pcReuse, useOldSolnForInitGuess, NcontrolParams, monitorConvergInterval);
//	#else
			petscSolver2 = new PetscSolver2(cvora, nbocv2_i, nbocv2_v, 5 + nScal, pcType, pcLevels, pcReuse, useOldSolnForInitGuess, NcontrolParams, monitorConvergInterval);
//	#endif

			petscSolver2->setTresholds(zeroAbs, zeroRel, maxIter);

			int gmresRestart = getIntParam("GMRES_RESTART", "100");
			petscSolver2->setGmresRestart(gmresRestart);
		}

		if(useOldSolnForInitGuess)
			petscSolver2->solveGMRES<MatT>(kspMonitorHistory, nIter, absResid, A, phi, rhs, phi,  cvora, nbocv_v_global, 5 + nScal, NcontrolParams, ncv_gg, vecC, d, step, newtonIter);
		else
			petscSolver2->solveGMRES<MatT>(kspMonitorHistory, nIter, absResid, A, phi, rhs, NULL, cvora, nbocv_v_global, 5 + nScal, NcontrolParams, ncv_gg, vecC, d, step, newtonIter);

		break;

#if defined(PETSC_HAVE_MUMPS) // This flag, PETSC_HAVE_MUMPS, is defined if MUMPS is installed while configuring PETSc.
                              // Note: For superLU -- PETSC_HAVE_SUPERLU or PETSC_HAVE_SUPERLU_DIST.
	case PETSC_LU_MUMPS:
		// on the first time, initialize the petsc solver...
		if (petscSolver2 == NULL) {
			if (nbocv_v_global == NULL) {
				nbocv_v_global = new int[ncv_gg];
				for (int icv = 0; icv < ncv; icv++)
					nbocv_v_global[icv] = cvora[mpi_rank] + icv;
				updateCvDataG1G2(nbocv_v_global, REPLACE_DATA);
			}

			int printLevel = 1; // MUMPS print level (0-4): MUMPS default = 2
			bool pcReuseDummy = false;
			petscSolver2 = new PetscSolver2(cvora, nbocv2_i, nbocv2_v, 5 + nScal, "LU", printLevel, pcReuseDummy, false, NcontrolParams);
//			petscSolver2->setTresholds(zeroAbs, zeroRel, maxIter);
		}
		petscSolver2->solveLUMUMPS<MatT>(A, phi, rhs, cvora, nbocv_v_global, 5 + nScal, NcontrolParams, ncv_gg, vecC, d, step, newtonIter);

		break;
#endif

	default:
		if (mpi_rank == 0)
			cerr << "Error: unrecognized solver: " << linearSolverNS
					<< endl;
		throw(-1);
		break;
	}
}

/*
 * Method: solveBasicMatOperationsPetsc
 * ------------------------------------
 * Solve for basic matrix operations (such as multiplication between a matrix and a vector) using Petsc
 *
 * whichOperation == 0: A*b=x
 */
template <class MatT>
void IkeWithPsALC_AD::solveBasicMatOperationsPetsc(double *xVec, MatT &A, double *bVec, double *cVec, const int whichOperation,
		int NcontrolParams /* = 0 */, int ncv_gg /* = 0 */, double **vecC /* = NULL */, const double *d /* = NULL */) {
	switch (whichOperation) {
	case 0: // A*b = x
		assert(cVec == NULL); // The cVec will not be used for this operation
		assert(xVec != NULL && bVec != NULL);

		// For the first time, initialize the petsc solver...
		if (petscSolver2 == NULL) {
			if (nbocv_v_global == NULL) {
				nbocv_v_global = new int[ncv_gg];
				for (int icv = 0; icv < ncv; icv++)
					nbocv_v_global[icv] = cvora[mpi_rank] + icv;
				updateCvDataG1G2(nbocv_v_global, REPLACE_DATA);
			}

			// Even though the basic matrix operations do not need PC information, you must pass it to the PetscSolver2 (or SlepcSolver2) object
			// for KSP operations since the PC information will not be updated if you once construct the object now.
			string pcType; // "BJACOBI", "ASM", "ILU", "AMG", "NONE"
			bool pcReuse;
			int pcLevels;  // Note: this is required only for the ILU pre-conditioner. Otherwise, any number doesn't matter.
			getPCcontext(pcType, pcReuse, pcLevels);

//	#ifdef USE_SLEPC_WITH_PETSC
//			petscSolver2 = new SlepcSolver2(cvora, nbocv2_i, nbocv2_v, 5 + nScal, pcType, pcLevels, pcReuse, false, NcontrolParams, monitorConvergInterval);
//	#else
			petscSolver2 = new PetscSolver2(cvora, nbocv2_i, nbocv2_v, 5 + nScal, pcType, pcLevels, pcReuse, false, NcontrolParams, monitorConvergInterval);
//	#endif
		}

		// Solve for the linear operation
		petscSolver2->matMultiplyVec(A, xVec, bVec, cvora, nbocv_v_global, 5 + nScal, NcontrolParams, ncv_gg, vecC, d);
		break;

	default:
		if (mpi_rank == 0)
			cerr << "Error: unrecognized solver: " << whichOperation << endl;
		throw(-1);
		break;
	}
}

/*
 * Method: getPCcontext
 * --------------------
 * Get the preconditioner contexts from the input file
 */
void IkeWithPsALC_AD::getPCcontext(string& pcType, bool &pcReuse, int& pcLevels) {
	if (!checkParam("PRECONDITIONER")) {
		ParamMap::add("PRECONDITIONER  METHOD=ASM  PC_REUSE=NO"); // add default values
		if (mpi_rank == 0)
			cout << "WARNING: added keyword \"PRECONDITIONER  METHOD=ASM  PC_REUSE=NO\" to parameter map!"<<endl;
	}

	pcType = getParam("PRECONDITIONER")->getString("METHOD"); // "BJACOBI", "ASM", "ILU", "AMG", "NONE"

	string tempString = getParam("PRECONDITIONER")->getString("PC_REUSE");
	pcReuse = false; // initially no!
	if(tempString=="YES" || tempString=="yes" || tempString=="Y" || tempString=="y")
		pcReuse = true;

	pcLevels = -1; // Note: this is required only for the ILU pre-conditioner. Otherwise, any number doesn't matter.
	if(pcType == "ILU")
		pcLevels = getIntParam("ILU_PC_LEVEL", "4");
}

#ifdef USE_SLEPC_WITH_PETSC
/*
 * Method: solveEigenProblem
 * -------------------------
 * Perform eigen-decomposition of the given matrix and get eigenvalues and eigenvectors
 *
 * Return:
 *   By reference:
 *     evalsReal = real part of the eigenvalues ( size: nev )
 *     evalsImag = imag part of the eigenvalues ( size: nev )
 *     evecsReal = real part of the eigenvectors ( size: nev * (ncv*(5+nScal)) )
 *     evecsImag = imag part of the eigenvectors ( size: nev * (ncv*(5+nScal)) )
 *     relError  = relative error (to the eigenvalue or to the matrix norms) for each eigenpair
 *     numIter   = total number of iterations for the iterative eigen solver
 *   By value:
 *     number of converged eigen-pairs
 *
 * Arguements:
 *   A     = matrix
 *   cvora =
 *   cv_gl =
 *   nScal =
 *   nev = the number of eigenvalues to compute
 *   ncv = the number of column vectors to be used by the solution algorithm
 *         It is recommended (depending on the method) to use ncv >= 2*nev or more.
 *   mpd = maximum projected dimension (for a more advanced usage -- Chaprter 2.6.4 in the SLEPc manual)
 *         In the case of a large number of nev, the computation costs can be reduced by setting mpd << nev.
 *   tol      = stopping criterion for the iterative eigen-solver, 2-norm (Default value = 1.0e-7)
 *   max_iter = maximum number of iteration for the iterative eigen-solver
 *   step       =
 *   newtonIter =
 */
template <class MatT>
int IkeWithPsALC_AD::solveEigenProblem(double *evalsReal, double *evalsImag, double **evecsReal, double **evecsImag, double *relError, int &numIter,
		MatT &A, int *cvora, const int nScal,
		const int nev, const int ncv, const int mpd, const double tol, const int max_iter,
		const int step /*=-1*/, const int newtonIter /*=-1*/) {
	int nconv = 0;

	SlepcSolver2 *slepcSolver = NULL; // Note: You must define a variable outside a switch command

	switch (linearSolverNS) {
	case PETSC_GMRES:
#ifdef USE_SLEPC_WITH_PETSC
		if (nbocv_v_global == NULL) {
			nbocv_v_global = new int[ncv_gg];
			for (int icv = 0; icv < ncv; icv++)
				nbocv_v_global[icv] = cvora[mpi_rank] + icv;
			updateCvDataG1G2(nbocv_v_global, REPLACE_DATA);
		}

		// Initiate the SLEPc solver...
//		if(petscSolver2 == NULL)
//			petscSolver2 = new SlepcSolver2(cvora, nbocv2_i, nbocv2_v, 5 + nScal);
//		else
		slepcSolver = new SlepcSolver2(cvora, nbocv2_i, nbocv2_v, 5 + nScal);

		// Solve the linear system
//		if(slepcSolver == NULL) // Note that petscSolver2 is seldom NULL since it should be allocated during the Newton-Raphson iterations
//			nconv = petscSolver2->solveEigenProblemSlepc<MatT>(evalsReal, evalsImag, evecsReal, evecsImag, relError, numIter,
//					                                      A, cvora, nbocv_v_global, nScal, ncv_gg, nev, ncv, mpd, tol, max_iter);
//		else
		nconv = slepcSolver ->solveEigenProblemSlepc<MatT>(evalsReal, evalsImag, evecsReal, evecsImag, relError, numIter,
					                                      A, cvora, nbocv_v_global, nScal, ncv_gg, nev, ncv, mpd, tol, max_iter);
			// Note: If ncv>nev, nconv can be larger than nev

		if(nconv == 0)
			if(mpi_rank==0)
				cout<<"WARNING in IkeWithPsALC_AD::solveEigenProblem(): No converged eigen-pair was found"<<endl;

		// Free the SLEPc solver
//		if(petscSolver2 != NULL) {
//			delete petscSolver2; 	petscSolver2 = NULL;
//		} else {
		delete slepcSolver; 	slepcSolver = NULL;
//		}
#else
		if (mpi_rank == 0)
			cerr << "ERROR in IkeWithPsALC_AD::solveEigenProblem(): USE_SLEPC_WITH_PETSC was not defined"<< endl;
		throw(-1);
#endif
		break;

	default:
		if (mpi_rank == 0)
			cerr << "ERROR in IkeWithPsALC_AD::solveEigenProblem(): unrecognized solver: " << linearSolverNS << endl;
		throw(-1);
		break;
	}

	return std::min<int>(nconv, nev);
}
#endif

///*
// * Method: initializeAdjoint
// * -------------------------
// * Original code: initializeAdjoint() in UgpWithCvCompFlowAD.h
// */
//void IkeWithPsALC_AD::initializeAdjoint(){
//	int nScal = scalarTranspEqVector.size();
//
//	scalarTranspEqVector_psi = new psi_AD[nScal];
//	for(int i=0; i<nScal ; i++)
//		scalarTranspEqVector_psi[i].phi      = new double[ncv_ggff];
//}

/*
 * Method: updatePrimVarsFrom1DVec
 * -------------------------------
 * Update the flow field primary variables (rho, rhou, rhoE, scals) from a 1D vector (q0 or q1)
 * Either q0 or q1 has been gotten from a binary file
 * Note: You must update nScal before calling this method
 */
void IkeWithPsALC_AD::updatePrimVarsFrom1DVec(const double* qVec, const int nScal) {
	int m = 5+nScal;
	for(int icv=0; icv<ncv; ++icv) {
		rho[icv] = qVec[icv*m];
		rhou[icv][0] = qVec[icv*m+1];
		rhou[icv][1] = qVec[icv*m+2];
		rhou[icv][2] = qVec[icv*m+3];
		rhoE[icv] = qVec[icv*m+4];

		if(nScal > 0) {
			for(int iScal=0; iScal<nScal; ++iScal) {
				double *phi = scalarTranspEqVector[iScal].phi;
				phi[icv] = qVec[icv*m+5+iScal];
			}
		}
	}
}

/*
 * Method: update1DVecFromPrimVars
 * -------------------------------
 * Update a 1D vector (q0 or q1) from the flow field primary variables (rho, rhou, rhoE, scals)
 * The primary variables can be obtained from a Hook function or a restart file, etc.
 * Note: You must update nScal before calling this method
 */
void IkeWithPsALC_AD::update1DVecFromPrimVars(double* qVec, const int nScal) {
	int m = 5+nScal;
	for(int icv=0; icv<ncv; ++icv) {
		qVec[icv*m] = rho[icv];
		qVec[icv*m+1] = rhou[icv][0];
		qVec[icv*m+2] = rhou[icv][1];
		qVec[icv*m+3] = rhou[icv][2];
		qVec[icv*m+4] = rhoE[icv];

		if(nScal > 0) {
			for(int iScal=0; iScal<nScal; ++iScal) {
				double *phi = scalarTranspEqVector[iScal].phi;
				qVec[icv*m+5+iScal] = phi[icv];
			}
		}
	}
}

/*
 * Method: updateFlow_primVars
 * ---------------------------
 * Update the flow field primary variables (Q won't be changed).
 * Arguments: Q = flow vector
 *            delQ = increment in Q
 *            relaxation = relaxation due to Backtracking. Usually sign = 1.0 (be careful in the case of PsALC: sign must be -1.0)
 *            nScal = number of scalars (e.g. If you use the k-omega turbulence model, iScal might be 2)
 * Note: you must call updateLambda(double*, double*, double, double, int) after calling this method
 */
void IkeWithPsALC_AD::updateFlow_primVars(const double* Q, const double* delQ, const double relaxation, const int nScal) {
	int m = 5 + nScal;

	/* update the flow field */
	for(int icv=0; icv<ncv; ++icv) {
		int icvTemp = icv*m;
		rho[icv] = Q[icvTemp] + relaxation*delQ[icvTemp];
		for(int i=0; i<3; ++i)
			rhou[icv][i] = Q[icvTemp+1+i] + relaxation*delQ[icvTemp+1+i];
		rhoE[icv] = Q[icvTemp+4] + relaxation*delQ[icvTemp+4];
		for(int iScal=0; iScal<nScal; ++iScal) {
			double *phi = scalarTranspEqVector[iScal].phi;
			phi[icv] = Q[icvTemp+5+iScal] + relaxation*delQ[icvTemp+5+iScal];
		}
	}

	/* update ghost cells */
	updateCvDataG1G2(rho, REPLACE_DATA);
	updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
	updateCvDataG1G2(rhoE, REPLACE_DATA);
	for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
		updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);

	calcStateVariables();     // Compute velocity, pressure, temperature, enthalpy and speed of sound at cell centers
//	calcMaterialProperties(); // Compute viscosity and thermal diffusivity at cell faces
}

/*
 * Method: updateLambda
 * --------------------
 * Update the flow field primary variables (Q won't be changed).
 * Arguments: Q = flow vector
 *            delQ = increment in Q
 *            relaxation = relaxation due to Backtracking. Usually sign = 1.0 (be careful in the case of PsALC: sign must be -1.0)
 *            NcontrolParams = number of control parameters
 *            nScal = number of scalars (e.g. If you use the k-omega turbulence model, nScal should be 2)
 */
void IkeWithPsALC_AD::updateLambda(const double* Q, const double* delQ, const double relaxation, const int NcontrolParams, const int nScal) {
	if(mpi_rank==mpi_size-1) {
		int Nflow = ncv*(5+nScal);
		for(int iEqn=0; iEqn<NcontrolParams; ++iEqn) {
			double lambdaBeforeUpdate = Q[Nflow+iEqn];
			if(fabs(relaxation*delQ[Nflow+iEqn]) > 0.1*fabs(lambdaBeforeUpdate))
				cout<<"           WARNING in updateLambda(): delta_lambda["<<iEqn<<"] (="<<relaxation*delQ[Nflow+iEqn]<<") is greater than 10% of lambda["<<iEqn<<"] (="<<lambdaBeforeUpdate<<")"<<endl;
			else if(fabs(relaxation*delQ[Nflow+iEqn]) < 1.0e-10)
				cout<<"           WARNING in updateLambda(): delta_lambda["<<iEqn<<"] (="<<relaxation*delQ[Nflow+iEqn]<<") is too small"<<endl;
			lambda[iEqn] = Q[Nflow+iEqn] + relaxation*delQ[Nflow+iEqn];
		}
	}
	MPI_Bcast(lambda, NcontrolParams, MPI_DOUBLE, mpi_size-1, mpi_comm);
	MPI_Barrier(mpi_comm);
}

/*
 * Method: updateFlow_1Dvec
 * ------------------------
 * Update the flow field (Q = Q + sign*delQ). The primary variables (rho, rhou, rhoE, scalars) won't be changed
 * Arguments: Q = flow vector
 *            delQ = increment in Q
 *            relaxation = magnitude of the increment. Usually sign = 1.0 (be careful in the case of PsALC: sign must be -1.0)
 *            N = size of the 1D flow vectors (Q and delQ)
 */
void IkeWithPsALC_AD::updateFlow_1Dvec(double* Q, const double* delQ, const double relaxation, const int N) {
	/* increase Q by delQ */
	for(int i=0; i<N; ++i)
		Q[i] = Q[i] + relaxation*delQ[i]; // be careful when you call this method while doing PsALC: "sign" should be negative (i.e. -1.0)
}

/*
 * Method: backtrackForNegativeVals
 * --------------------------------
 * Update relaxation to prevent negative density, pressure, or kine at both CVs and FAs
 */
double IkeWithPsALC_AD::backtrackForNegativeVals(int &negativeValCount_CV, int &negativeValCount_FA, const double clipParameter, const double safeParameter, 
                                        const double relaxation, const int kine_index, const int nScal, const double* qArray, const double* delQ, const int iterMax) {
    double newRelaxation = relaxation;

    int m = 5 + nScal;

    /* Parameters */
    double rhoRef   = RefFlowParams.rho_ref;
    double pressRef = RefFlowParams.press_ref;
    double kineRef  = RefFlowParams.kine_ref;

    double rhoClipBound   = clipParameter*rhoRef;
    double pressClipBound = clipParameter*pressRef;
    double kineClipBound  = clipParameter*kineRef;

    /* Check negative values at CVs and calculate new relaxation */
    int myNegativeCount_CV = 0;
    double myRelax = relaxation;
    
    for(int icv=0; icv<ncv; ++icv) {
        int icvTemp = icv*m;
   
        bool negativeFound = false;

        if(rho[icv] < rhoClipBound) {
            double oldRho = qArray[icvTemp];
        
            double tempRelax = safeParameter*(oldRho-rhoClipBound)/fabs(delQ[icvTemp]);
            myRelax = min(myRelax, tempRelax);
            ++myNegativeCount_CV;

            negativeFound = true;
        }
    
        // Note: You cannot simply use UgpWithCvCompFlow::press instead of "newKinecv" below
        //       because UgpWithCvCompFlow::calcStateVariables() gives up updating "press" if press becomes negative
        double newKinecv  = 0.0;
//        if(UgpWithCvCompFlow::kine != NULL)
//            newKinecv = UgpWithCvCompFlow::kine[icv];
        if(kine_index > -1)
            newKinecv = UgpWithCvCompFlow::scalarTranspEqVector[kine_index].phi[icv];

        double newRhouMag = vecDotVec3d(UgpWithCvCompFlow::rhou[icv], UgpWithCvCompFlow::rhou[icv]);
        double newPress = (UgpWithCvCompFlow::rhoE[icv] - 0.5*newRhouMag/UgpWithCvCompFlow::rho[icv] - UgpWithCvCompFlow::rho[icv]*newKinecv) * (UgpWithCvCompFlow::gamma[icv] - 1.0);
        
        double oldKinecv = 0.0;
        if(kine_index>-1)
            oldKinecv = qArray[icvTemp+5+kine_index];
           
		if(newPress < pressClipBound) {
			double oldRho = qArray[icvTemp];
			double oldRhouMag  = 0.0;
			for(int i=0; i<3; ++i)
				oldRhouMag += pow(qArray[icvTemp+1+i], 2.0);
            double oldPress = (qArray[icvTemp+4] - 0.5*oldRhouMag/oldRho - qArray[icvTemp]*oldKinecv) * (UgpWithCvCompFlow::gamma[icv]-1.0);
                    
            double newRelax = safeParameter*(oldPress-pressClipBound)/fabs(newPress-oldPress);
            myRelax = min(myRelax, newRelax);

            if(!negativeFound) 
                ++myNegativeCount_CV;
            negativeFound = true;
		}
    
        if(kine_index>-1) {
            if(newKinecv < kineClipBound) {
                double newRelax = safeParameter*(oldKinecv-kineClipBound)/fabs(delQ[icvTemp+5+kine_index]);
                myRelax = min(myRelax, newRelax);

                if(!negativeFound)
                    ++myNegativeCount_CV;
                negativeFound = true;
            }
        }

        if(nScal>0 && turbModel != NONE) {
        	int negativeScalar = checkNegativeTurbScalarCv_JOE(myRelax, myRelax, icv, clipParameter, safeParameter, qArray, delQ); // Check negative scalar EXCEPT kine
        	if(negativeScalar > 0) {
        		if(!negativeFound)
        			++myNegativeCount_CV;
        		negativeFound = true;
        	}
        }
    }
    
    negativeValCount_CV = 0;
    MPI_Allreduce(&myNegativeCount_CV, &negativeValCount_CV, 1, MPI_INT, MPI_SUM, mpi_comm);
    if(negativeValCount_CV>0) {
        MPI_Allreduce(&myRelax, &newRelaxation, 1, MPI_DOUBLE, MPI_MIN, mpi_comm);
        newRelaxation = min(relaxation, max(newRelaxation, clipParameter)); // To prevent 0.0: Note that clipParameter is usually a very small number
    }
    
    /* Check negative values at FAs and iteratively calculate local new relaxation */
    int myNegativeCount_FA = 0;
    negativeValCount_FA = 0;
    if (sndOrder == true) { // Note: Unless sndOrder, face values are NOT re-constructed by using gradient of the variable BUT adopted cell-center values
        for(int iIter=0; iIter<iterMax; ++iIter) {
            myNegativeCount_FA  = 0; // Reset to zero
            negativeValCount_FA = 0;
        
            // Update the flow field based on the relaxation (Note: lambda will not be updated here)
            updateFlow_primVars(qArray, delQ, -newRelaxation, nScal); // note: You should not update q here! (because of backtracking)
                                                                      //       Only rho, rhou, rhoE, and scalars should be updated here

            // Calculate state variables
            calcStateVariables(false);
        
            // Calculate gradients
            calcCv2Grad(UgpWithCvCompFlow::grad_rho,  UgpWithCvCompFlow::rho,   limiterNavierS, UgpWithCvCompFlow::rho,   epsilonSDWLS);
#ifdef temp_reconstruction
            calcCv2Grad(UgpWithCvCompFlow::grad_temp, UgpWithCvCompFlow::temp,  limiterNavierS, UgpWithCvCompFlow::temp,  epsilonSDWLS);
#else
            calcCv2Grad(UgpWithCvCompFlow::grad_p,    UgpWithCvCompFlow::press, limiterNavierS, UgpWithCvCompFlow::press, epsilonSDWLS);
#endif

            for (int iScal = 0; iScal < nScal; iScal++) {
                if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE") {
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
                } else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD") {
                    double *phi = scalarTranspEqVector[iScal].phi;
                    double (*grad_phi)[3] = scalarTranspEqVector[iScal].grad_phi;
                    calcCv2Grad(grad_phi, phi, limiterNavierS, phi, epsilonSDWLS);
                }
            }

            // Calculate the values at the faces
            // ===============================================================================================
            // cycle through internal faces
            // ===============================================================================================
            for (int ifa = nfa_b; ifa < nfa; ifa++) {
                int icv0 = cvofa[ifa][0];
                int icv1 = cvofa[ifa][1];
                assert( icv0 >= 0 );
                assert( icv1 >= 0 );

                bool negativeFound = false;
                
                // face unit normal and area...
                double nVec[3] = {0.0, 0.0, 0.0};
                double area    = normVec3d(nVec, fa_normal[ifa]);
                
                // reconstruction at the faces
                double rho0 = UgpWithCvCompFlow::rho[icv0];
                double p0   = UgpWithCvCompFlow::press[icv0];
                double T0   = UgpWithCvCompFlow::temp[icv0];
                
                double rho1 = UgpWithCvCompFlow::rho[icv1];
                double p1   = UgpWithCvCompFlow::press[icv1];
                double T1   = UgpWithCvCompFlow::temp[icv1];
                
                double r0[3] = {0.0, 0.0, 0.0}, r1[3] = {0.0, 0.0, 0.0};
                vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);
                vecMinVec3d(r1, x_fa[ifa], x_cv[icv1]);
                
                // ----------------------------------------
                // left side
                // ----------------------------------------
                rho0 += vecDotVec3d(r0, UgpWithCvCompFlow::grad_rho[icv0]);
                if(rho0 <= 0.0) {
                    ++myNegativeCount_FA;
                    negativeFound = true;
                }
#ifdef temp_reconstruction
                T0 += vecDotVec3d(r0, UgpWithCvCompFlow::grad_temp[icv0]);
                if (T0 <= 0.0) {
                    if(!negativeFound)
                        ++myNegativeCount_FA;
                    negativeFound = true;
                }
#else
                p0 += vecDotVec3d(r0, UgpWithCvCompFlow::grad_p[icv0]);
                if (p0 <= 0.0) {
                    if(!negativeFound)
                        ++myNegativeCount_FA;
                    negativeFound = true;
                }
#endif

                if(kine_index > -1) {
                    int iScal = kine_index;

                    double (*grad_phi)[3] = UgpWithCvCompFlow::scalarTranspEqVector[iScal].grad_phi;

                    double scalar0 = UgpWithCvCompFlow::scalarTranspEqVector[iScal].phi[icv0];
                    double rho0    = UgpWithCvCompFlow::rho[icv0];
                    if (strcmp(scalarTranspEqVector[iScal].reconstruction.c_str(), "CONSERVATIVE")==0)
                        scalar0 = (rho0 * scalar0 + vecDotVec3d(r0, grad_phi[icv0])) / rho0;
                    else if (strcmp(scalarTranspEqVector[iScal].reconstruction.c_str(), "STANDARD")==0)
                        scalar0 += vecDotVec3d(r0, grad_phi[icv0]);

                    if(scalar0 <= 0.0) {
                        if(!negativeFound)
                            ++myNegativeCount_FA;
                        negativeFound = true;
                    }
                }

                // ----------------------------------------
                // right side
                // ----------------------------------------
                rho1 += vecDotVec3d(r1, UgpWithCvCompFlow::grad_rho[icv1]);
                if(rho1 <= 0.0) {
                    if(!negativeFound)
                        ++myNegativeCount_FA;
                    negativeFound = true;
                }
#ifdef temp_reconstruction
                T1 += vecDotVec3d(r1, UgpWithCvCompFlow::grad_temp[icv1]);
                if (T1 <= 0.0) {
                    if(!negativeFound)
                        ++myNegativeCount_FA;
                    negativeFound = true;
                }
#else
                p1 += vecDotVec3d(r1, UgpWithCvCompFlow::grad_p[icv1]);
                if (p1 <= 0.0) {
                    if(!negativeFound)
                        ++myNegativeCount_FA;
                    negativeFound = true;
                }
#endif

                if(kine_index > -1) {
                    int iScal = kine_index;

                    double (*grad_phi)[3] = UgpWithCvCompFlow::scalarTranspEqVector[iScal].grad_phi;

                    double scalar1 = UgpWithCvCompFlow::scalarTranspEqVector[iScal].phi[icv1];
                    double rho1    = UgpWithCvCompFlow::rho[icv1];
                    if (strcmp(scalarTranspEqVector[iScal].reconstruction.c_str(), "CONSERVATIVE")==0)
                        scalar1 = (rho1 * scalar1 + vecDotVec3d(r1, grad_phi[icv1])) / rho1;
                    else if (strcmp(scalarTranspEqVector[iScal].reconstruction.c_str(), "STANDARD")==0)
                        scalar1 += vecDotVec3d(r1, grad_phi[icv1]);

                    if(scalar1 <= 0.0) {
                        if(!negativeFound)
                            ++myNegativeCount_FA;
                        negativeFound = true;
                    }
                }

                // ----------------------------------------
                // turbulent scalar EXCEPT kine at both faces
                // ----------------------------------------
                if(nScal>0 && turbModel != NONE) {
					int negativeTurbScalar = checkNegativeTurbScalarFa_JOE(ifa, icv0, icv1, true);
					if(negativeTurbScalar > 0) {
						if(!negativeFound)
							++myNegativeCount_FA;
						negativeFound  = true;
					}
                }
            }
            
            // ===============================================================================================
            // cycle through boundary faces
            // ===============================================================================================
            for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
                if (zone->getKind() == FA_ZONE_BOUNDARY) {
                    Param *param;
                    
                    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ++ifa) {
                        int icv0 = cvofa[ifa][0];
                        int icv1 = cvofa[ifa][1];
                        assert( icv0 >= 0 );

                        bool negativeFound = false;
                        
                        // face unit normal and area...
                        double nVec[3] = {0.0, 0.0, 0.0};
                        double area    = normVec3d(nVec, fa_normal[ifa]);
                        
                        // reconstruction at the faces
                        double rho0 = UgpWithCvCompFlow::rho[icv0];
                        double p0   = UgpWithCvCompFlow::press[icv0];
                        double T0   = UgpWithCvCompFlow::temp[icv0];
                        
                        double r0[3] = {0.0, 0.0, 0.0};
                        vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);
                        
                        // ----------------------------------------
                        // left side
                        // ----------------------------------------
                        rho0 += vecDotVec3d(r0, UgpWithCvCompFlow::grad_rho[icv0]);
                        if(rho0 <= 0.0) {
                            ++myNegativeCount_FA;
                            negativeFound = true;
                        }
#ifdef temp_reconstruction
                        T0 += vecDotVec3d(r0, UgpWithCvCompFlow::grad_temp[icv0]);
                        if (T0 <= 0.0) {
                            if(!negativeFound)
                                ++myNegativeCount_FA;
                            negativeFound = true;
                        }
#else
                        p0 += vecDotVec3d(r0, UgpWithCvCompFlow::grad_p[icv0]);
                        if (p0 <= 0.0) {
                            if(!negativeFound)
                                ++myNegativeCount_FA;
                            negativeFound = true;
                        }
#endif

                        if(kine_index > -1) {
                            int iScal = kine_index;

                            double (*grad_phi)[3] = UgpWithCvCompFlow::scalarTranspEqVector[iScal].grad_phi;

                            double scalar0 = UgpWithCvCompFlow::scalarTranspEqVector[iScal].phi[icv0];
                            double rho0    = UgpWithCvCompFlow::rho[icv0];
                            if (strcmp(scalarTranspEqVector[iScal].reconstruction.c_str(), "CONSERVATIVE")==0)
                                scalar0 = (rho0 * scalar0 + vecDotVec3d(r0, grad_phi[icv0])) / rho0;
                            else if (strcmp(scalarTranspEqVector[iScal].reconstruction.c_str(), "STANDARD")==0)
                                scalar0 += vecDotVec3d(r0, grad_phi[icv0]);

                            if(scalar0 <= 0.0) {
                                if(!negativeFound)
                                    ++myNegativeCount_FA;
                                negativeFound = true;
                            }
                        }

                        // ----------------------------------------
                        // turbulent scalar EXCEPT kine at both faces
                        // ----------------------------------------
                        if(nScal>0 && turbModel != NONE) {
							int negativeTurbScalar = checkNegativeTurbScalarFa_JOE(ifa, icv0, icv1, false);
							if(negativeTurbScalar > 0) {
								if(!negativeFound)
									++myNegativeCount_FA;
								negativeFound = true;
							}
                        }
                    }
                }
            }
            
            MPI_Allreduce(&myNegativeCount_FA, &negativeValCount_FA, 1, MPI_INT, MPI_SUM, mpi_comm);
            if(negativeValCount_FA>0)
                newRelaxation *= 0.5;
            
            if(negativeValCount_FA==0 || newRelaxation<safeParameter) // To prevent 0.0
                break;
        }

        MPI_Barrier(mpi_comm);
    }

    /* Return */
    return newRelaxation;
}

/*
 * Method: backtrackWithJOE_calcRelaxAndRHS
 * ----------------------------------------
 * Launch the backtracking algorithm: calculate "relax" and recalculate RHS and Resid
 * This method is called by the Newton-Raphson method when
 *   1. residual increases
 *   2. density, pressure ( or kine) become negative
 * For details of the algorithm in this method, see C.T.Kelley, Iterative Methods for Linear and Nonlinear Equations, SIAM 1995, Chap.8.3.1, "Three-point parabolic model", pp.143-144
 *
 * Note:
 *   1. The primary variables (i.e., rho, rhou, rhoE, scalars) will be updated, but the Q array won't.
 *   2. Use JoeWithModels::calcRhs()
 *
 * Arguments: rhs = the residual array
 *            qArray = flow vector
 *            delQ = increment in Q
 *            notConverged = true if backtracking fails to drop residual
 *            maxBacktrackingIter = maximum number of backtracking
 *            relaxationOld =
 *            relaxationLowerBound =
 *            relaxationUpperBound =
 *            ResidVecOld = the residual array (previous step)
 *            ResidVecNew = the residual array (next step) -- Return by argument
 *            whichNorm
 *            backtrackBarrierCoeff = barrier coefficient which is used to check if the residual increased or not
 *                                    (i.e. (1.0-backtrackBarrierCoeff*relaxation)*residNormTotOld)
 *            ----- BELOW is FOR TANGENTIAL CONDITIONS -----
 *            NcontrolEqns   = Number of control parameters
 *            q_tangent      = Tangential vector of q
 *            lambda_tangent = Tangential vector of lambda
 *            weightLambda   = the weight for lambdas in the weighted 2-norm
 *            qArrayOld      = Flow vector at the previous time step
 *            NresOld        = Residual of the tangential equations at the previous time step
 *            lambdaOld      = lambda at the previous time step
 *            arcLength      = pesudo-arclength
 * Return: By arguments - The following 2 variables will be updated
 *                        1. rhs
 *                        2. ResidNew
 *                        3. notConverged
 *         By return - final relaxation in step
 */
double IkeWithPsALC_AD::backtrackWithJOE_calcRelaxAndRHS(double* rhs, const double* qArray, const double* delQ, bool &notConverged,
		const int maxBacktrackingIter, const double relaxationOld, const double relaxationLowerBound, const double relaxationUpperBound,
		const double* ResidVecOld, double* ResidVecNew, const int whichNorm, const double backtrackBarrierCoeff,
		const int NcontrolEqns, double** q_tangent, const double* lambda_tangent, const double weightLambda,
		const double* qArrayOld, const double* NresOld, const double* lambdaOld, const double arcLength) {
	// Store the history of backtracking
	vector<pair<double, double> > backtrackHistory; // This vector stores the history of <relax, residNormTot> from relax1 to the final backtracking
	                                                // Note: If negative values (e.g. rho) occurs, then ABSURDLY_BIG_NUMBER will be stored instead of residNormTot.

	// ----------------------------
	int nVars = 5 + nScal;
	int N = ncv*nVars; // size of the flow vector (Q and delQ)

	double relax0 = 0.0;
	double relax1 = relaxationOld; // Usually, 1.0. But if other relaxation algorithm (e.g. trust-region) was active, it can be less than 1.0.
	double relax2 = relaxationOld * relaxationUpperBound;

	double *ResidVecTemp = new double [nVars];
	double *NresTemp = NULL;
	if (NcontrolEqns>0)
		NresTemp = new double [NcontrolEqns];

	// ----------------------------
	double residNormFlow0 = calcSumResidual(ResidVecOld, whichNorm);
	double residNormTot0;
	if (NcontrolEqns > 0)
		residNormTot0 = updateTotResidualWithNres(residNormFlow0, NresOld, NcontrolEqns, whichNorm);
	else
		residNormTot0 = residNormFlow0;

	double residNormFlow1 = calcSumResidual(ResidVecNew, whichNorm);
	double residNormTot1;
#ifdef USE_TOT_NORM_WITH_LAMBDA
	if (NcontrolEqns > 0) {
		bool tangentialSatisfied = calcNresForJOE(NresTemp, nVars, NcontrolEqns, qArrayOld, q_tangent, lambda_tangent, lambda, lambdaOld, weightLambda, arcLength);

		residNormTot1 = updateTotResidualWithNres(residNormFlow1, NresTemp, NcontrolEqns, whichNorm);
	} else
		residNormTot1 = residNormFlow1;
#else
	residNormTot1 = residNormFlow1;
#endif
	backtrackHistory.push_back(std::make_pair(relax1, residNormTot1));

	double residNormFlow2 = ABSURDLY_BIG_NUMBER; // Initialize with a very large value for debugging
	double residNormTot2  = ABSURDLY_BIG_NUMBER; // Initialize with a very large value for debugging

	double residNormFlowOld = residNormFlow0;
	double residNormTotOld  = residNormTot0;
	double residNormFlowInit = residNormFlow1; // The initial residual of next(?) step before launching the backtracking algorithm
	double residNormTotInit  = residNormTot1;  // The initial residual of next(?) step before launching the backtracking algorithm

	// ----------------------------
	/*
	 * Initial Backtracking
	 *   In order to use the three-point parabolic model, you should launch the first backtracking anyway
	 */
	// Update the flow field and lambda
	updateFlow_primVars(qArray, delQ, -relax2, nScal);

#ifdef UPDATE_LAMBDA_IN_BACKTRACK
	if (NcontrolEqns > 0)
		updateLambda(qArray, delQ, -relax2, NcontrolEqns, nScal);
#endif

	// Calculate RHS and residual
	bool useBarrier = true;

	int countNegative = calcRhsWithBarrier(rhs, useBarrier);
	calcResidualsFrom1Drhs(ResidVecTemp, rhs, whichNorm); // Residual in terms of each flow variable (size = nVars)
	residNormFlow2 = calcSumResidual(ResidVecTemp, whichNorm);

#ifdef USE_TOT_NORM_WITH_LAMBDA
	if (NcontrolEqns > 0) {
		bool tangentialSatisfied = calcNresForJOE(NresTemp, nVars, NcontrolEqns, qArrayOld, q_tangent, lambda_tangent, lambda, lambdaOld, weightLambda, arcLength);

		residNormTot2 = updateTotResidualWithNres(residNormFlow2, NresTemp, NcontrolEqns, whichNorm);
	} else
		residNormTot2 = residNormFlow2;
#else
	residNormTot2 = residNormFlow2;
#endif
	if(countNegative > 0)
		backtrackHistory.push_back(std::make_pair(relax2, ABSURDLY_BIG_NUMBER));
	else
		backtrackHistory.push_back(std::make_pair(relax2, residNormTot2));

	// show residual on screen
	if(debugLevel > 1 && mpi_rank == 0)
		printf("                  >> backtrackWithJOE_calcRelaxAndRHS(): MaxIter=%d, relaxLowBound=%.3e, relaxUpperBound=%.3e \n", maxBacktrackingIter, relaxationLowerBound, relaxationUpperBound);
	if(debugLevel > 1 && mpi_rank == 0) {
#ifdef USE_TOT_NORM_WITH_LAMBDA
		printf("                    >> Initial backtracking (RELAX=%.4e) -- RESID_FLOW=%.8e, RESID_TOT=%.8e", relax2, residNormFlow2, residNormTot2);
#else
		printf("                    >> Initial backtracking (RELAX=%.4e) -- RESID_TOT=%.8e", relax2, residNormTot2);
#endif
		if(lambda != NULL)
			printf(", LAMBDA=%.7e", lambda[0]);
		printf("\n");
	}

	/*
	 * Three-point parabolic backtracking algorithm
	 */
	int iter = 2; // Note that you already run one backtracking before jumping into the following while loop

	bool notConverged_check1 = checkBacktrackFromResidInc(residNormTot1, residNormTot2, relax2, backtrackBarrierCoeff);

	bool notConverged_check2 = checkBacktrackFromResidInc(residNormTotOld, residNormTot2, relax2, backtrackBarrierCoeff);

	// Make sure that the final result is less than the last step
	notConverged = (notConverged_check1 || notConverged_check2);
	while(countNegative>0 || notConverged) {
		// If relation becomes too small, break the iterative procedure
		if(relax2 < RELAXATION_EPS) {
			// Note: You should not make "notConverged" true since it will be returned by reference!

			if(debugLevel > 0 && mpi_rank == 0)
				printf("                  >> backtrackWithJOE_calcRelaxAndRHS(): RELAXATION=%.3e becomes too small at iter = %d\n", relax2, iter-1);

			relax2 = min(RELAXATION_EPS, relaxationUpperBound);

			// Update the flow field
			updateFlow_primVars(qArray, delQ, -relax2, nScal);

#ifdef UPDATE_LAMBDA_IN_BACKTRACK
			if (NcontrolEqns > 0)
				updateLambda(qArray, delQ, -relax2, NcontrolEqns, nScal);
#endif

			// calculate the norm of the residual
			MPI_Barrier(mpi_comm);
			countNegative = calcRhsWithBarrier(rhs, useBarrier);
			calcResidualsFrom1Drhs(ResidVecTemp, rhs, whichNorm);
			residNormFlow2 = calcSumResidual(ResidVecTemp, whichNorm);
#ifdef USE_TOT_NORM_WITH_LAMBDA
			if (NcontrolEqns > 0) {
				bool tangentialSatisfied = calcNresForJOE(NresTemp, nVars, NcontrolEqns, qArrayOld, q_tangent, lambda_tangent, lambda, lambdaOld, weightLambda, arcLength);

				residNormTot2 = updateTotResidualWithNres(residNormFlow2, NresTemp, NcontrolEqns, whichNorm);
			} else
				residNormTot2 = residNormFlow2;
#else
			residNormTot2 = residNormFlow2;
#endif
			if(countNegative > 0)
				backtrackHistory.push_back(std::make_pair(relax2, ABSURDLY_BIG_NUMBER));
			else
				backtrackHistory.push_back(std::make_pair(relax2, residNormTot2));

			// Break the while() routine
			if(debugLevel > 0 && mpi_rank == 0)
				printf("                                                         Set RELAXATION = %.3e (RESID_TOT=%.8e) and Break\n", relax2, residNormTot2);
			break;
		}

		/*
		 * check maximum iter
		 */
		if(iter > maxBacktrackingIter) {
			if(countNegative>0) {
				if(mpi_rank==0) {
					cerr<<"ERROR in backtrackWithJOE_calcRelaxAndRHS()! Backtracking can never go to the feasible region within iter="<<iter-1<<endl;
					throw(PSALC_ERROR_CODE);
				}
			} else {
				if(mpi_rank==0)
					cout<<"                  backtrackWithJOE_calcRelaxAndRHS(): WARNING!! Too many iterations for backtracking: iter="<<iter-1<<". Break."<<endl;
				break;
			}
		}

		/*
		 * Calculate "relax": For details of the the three-point parabolic model, see C.T.Kelley
		 * Short descriptions:
		 *  First, construct a 2nd order polynomial p(relax) in terms of "relax" using the points of (relax0, resid0), (relax1, resid1), and (relax2, resid2).
		 *  (Note that in the initial iteration of the while() roop, relax0 = 0.0(previous step), relax1 = 1.0(no backtracking), relax2=relaxation (user specified argument)
		 *  Then, get "relax" for the smallest p(relax)
		 *  Finally, apply a safeguarding step.
		 */
		// get "relax" from the 2nd order polynomial
		double relaxNew;
		double pDD = 2*(relax2*(residNormTot1-residNormTot0)-relax1*(residNormTot2-residNormTot0)) / (relax1*relax2*(relax1-relax2));
		if(pDD <= 0.0) {
			relaxNew = relaxationUpperBound*relax2;
		} else {
			double pD = (relax1*relax1*(residNormTot2-residNormTot0) - relax2*relax2*(residNormTot1-residNormTot0)) / (relax1*relax2*(relax1-relax2));
			relaxNew = -pD/pDD;
		}
		// apply the safeguarding step
		relaxNew = max(relaxNew, relaxationLowerBound*relax2);
		relaxNew = min(relaxNew, relaxationUpperBound*relax2);

		/*
		 * Update relax0, relax1, relax2, resid0, and resid1
		 */
		relax0 = relax1;
		relax1 = relax2;
		relax2 = relaxNew;
		residNormFlow0 = residNormFlow1;
		residNormTot0  = residNormTot1;
		residNormFlow1 = residNormFlow2;
		residNormTot1  = residNormTot2;

		/*
		 * Update the flow field and lambda
		 */
		updateFlow_primVars(qArray, delQ, -relax2, nScal);

#ifdef UPDATE_LAMBDA_IN_BACKTRACK
		if (NcontrolEqns > 0)
			updateLambda(qArray, delQ, -relax2, NcontrolEqns, nScal);
#endif

		/*
		 * calculate the norm of the residual
		 */
		MPI_Barrier(mpi_comm);
		countNegative = calcRhsWithBarrier(rhs, useBarrier);
		calcResidualsFrom1Drhs(ResidVecTemp, rhs, whichNorm);
		residNormFlow2 = calcSumResidual(ResidVecTemp, whichNorm);
#ifdef USE_TOT_NORM_WITH_LAMBDA
		if (NcontrolEqns > 0) {
			bool tangentialSatisfied = calcNresForJOE(NresTemp, nVars, NcontrolEqns, qArrayOld, q_tangent, lambda_tangent, lambda, lambdaOld, weightLambda, arcLength);

			residNormTot2 = updateTotResidualWithNres(residNormFlow2, NresTemp, NcontrolEqns, whichNorm);
		} else
			residNormTot2 = residNormFlow2;
#else
		residNormTot2 = residNormFlow2;
#endif
		if(countNegative > 0)
			backtrackHistory.push_back(std::make_pair(relax2, ABSURDLY_BIG_NUMBER));
		else
			backtrackHistory.push_back(std::make_pair(relax2, residNormTot2));

		/*
		 * show residual on screen
		 */
		if(debugLevel > 1 && mpi_rank == 0) {
			if(NcontrolEqns>0) {
				switch(iter) {
#ifdef USE_TOT_NORM_WITH_LAMBDA
				case 2:
					printf("                    >> %5dnd backtracking (RELAX=%.4e) -- RESID_FLOW=%.8e, RESID_TOT=%.8e", iter,relax2,residNormFlow2,residNormTot2);
					break;
				case 3:
					printf("                    >> %5drd backtracking (RELAX=%.4e) -- RESID_FLOW=%.8e, RESID_TOT=%.8e", iter,relax2,residNormFlow2,residNormTot2);
					break;
				default:
					printf("                    >> %5dth backtracking (RELAX=%.4e) -- RESID_FLOW=%.8e, RESID_TOT=%.8e", iter,relax2,residNormFlow2,residNormTot2);
					break;
#else
				case 2:
					printf("                    >> %5dnd backtracking (RELAX=%.4e) -- RESID_TOT=%.8e", iter,relax2,residNormTot2);
					break;
				case 3:
					printf("                    >> %5drd backtracking (RELAX=%.4e) -- RESID_TOT=%.8e", iter,relax2,residNormTot2);
					break;
				default:
					printf("                    >> %5dth backtracking (RELAX=%.4e) -- RESID_TOT=%.8e", iter,relax2,residNormTot2);
					break;
#endif
				}
				printf(", LAMBDA=%.7e \n", lambda[0]);
			} else {
				switch(iter) {
				case 2:
					printf("                    >> %5dnd backtracking (RELAX=%.4e) -- RESID_TOT=%.8e \n", iter,relax2,residNormTot2);
					break;
				case 3:
					printf("                    >> %5drd backtracking (RELAX=%.4e) -- RESID_TOT=%.8e \n", iter,relax2,residNormTot2);
					break;
				default:
					printf("                    >> %5dth backtracking (RELAX=%.4e) -- RESID_TOT=%.8e \n", iter,relax2,residNormTot2);
					break;
				}
			}
		}

		/*
		 * increase iter
		 */
		++iter;

		/*
		 * Check if converged:
		 *   Criterion: 1. Residual decreases compared to the previous backtracking step
		 *              2. Residual is smaller than the initial residual (Important !!)
		 */
		bool notConverged_check1 = checkBacktrackFromResidInc(residNormTot1, residNormTot2, relax2, backtrackBarrierCoeff);

		bool notConverged_check2 = checkBacktrackFromResidInc(residNormTotOld, residNormTot2, relax2, backtrackBarrierCoeff);

		notConverged = (notConverged_check1 || notConverged_check2);
	}

	/*
	 * Search for the minimum residual in "backtrackHistory" if backtracking fails
	 * Note: backtrackHistory[i].second for i>=1 was set as ABSURDLY_BIG_NUMBER instead of residNormTot2 if countNegative>0
	 */
	if(countNegative>0 || notConverged) {
		double minResid     = (countNegative>0) ? ABSURDLY_BIG_NUMBER : residNormTot2;
		int    minResidIter = backtrackHistory.size()-1;
		for(size_t i=0; i<backtrackHistory.size(); ++i) {
			if(backtrackHistory[i].second < minResid) {
				relax2        = backtrackHistory[i].first;
				residNormTot2 = backtrackHistory[i].second;

				minResid     = backtrackHistory[i].second;
				minResidIter = i;
			}
		}

		// Update the flow field using the "optimal" relax2
		updateFlow_primVars(qArray, delQ, -relax2, nScal);

#ifdef UPDATE_LAMBDA_IN_BACKTRACK
		if (NcontrolEqns > 0)
			updateLambda(qArray, delQ, -relax2, NcontrolEqns, nScal);
#endif

		// calculate the norm of the residual at the "optimal" relax2: You should calculate "rhs" anyways...
		MPI_Barrier(mpi_comm);
		int countNegativeTemp = calcRhsWithBarrier(rhs, useBarrier);
		calcResidualsFrom1Drhs(ResidVecTemp, rhs, whichNorm);
		residNormFlow2 = calcSumResidual(ResidVecTemp, whichNorm);
#ifdef USE_TOT_NORM_WITH_LAMBDA
		if (NcontrolEqns > 0) {
			bool tangentialSatisfied = calcNresForJOE(NresTemp, nVars, NcontrolEqns, qArrayOld, q_tangent, lambda_tangent, lambda, lambdaOld, weightLambda, arcLength);

			residNormTot2 = updateTotResidualWithNres(residNormFlow2, NresTemp, NcontrolEqns, whichNorm);
		} else
			residNormTot2 = residNormFlow2;
#else
		residNormTot2 = residNormFlow2;
#endif

		// Show the result on the screen
		if(debugLevel > 0 && mpi_rank == 0) {
			double relax2Max = relaxationOld*relaxationUpperBound;
			double threshold1 = calcBacktrackThreshold(residNormTot1,   relax2Max, backtrackBarrierCoeff);
			double threshold2 = calcBacktrackThreshold(residNormTotOld, relax2Max, backtrackBarrierCoeff);
			printf("                  >> backtrackWithJOE_calcRelaxAndRHS(): Backtracking result: iter = %d  RELAX=%.3e  RESID_TOT=%.5e --> %.5e (> THRESHOLD=%.4e)\n",
					minResidIter, relax2, residNormTotInit, residNormTot2, min(threshold1, threshold2));
		}
	}

	/*
	 * Show the backtracking result on the screen if backtracking is successful
	 */
	if(countNegative==0 && !notConverged)
		if(debugLevel > 0 && mpi_rank == 0) {

			printf("                  >> backtrackWithJOE_calcRelaxAndRHS(): Backtracking successful! iter = %d  RELAX=%.3e  RESID_TOT=%.5e --> %.5e\n",
					iter-1, relax2, residNormTotInit, residNormTot2);
		}

	/* update residNorm */
	for(int i=0; i<5+nScal; ++i)
		ResidVecNew[i] = ResidVecTemp[i];

	/* clear memory */
	delete [] ResidVecTemp;
	if (NcontrolEqns>0)
		delete [] NresTemp;

	/* return */
	return relax2;
}

/*
 * Method: backtrackWithJOEcoupled_calcRelaxAndRHS
 * -----------------------------------------------
 * Launch the backtracking algorithm: calculate "relax" and recalculate RHS and Resid
 * This method is called by the Newton-Raphson method when
 *   1. residual increases
 *   2. density, pressure ( or kine) become negative
 *
 * For details of the algorithm in this method, see C.T.Kelley, Iterative Methods for Linear and Nonlinear Equations, SIAM 1995, Chap.8.3.1, "Three-point parabolic model", pp.143-144
 *
 * Note:
 *   The primary variables (i.e., rho, rhou, rhoE, scalars) will be updated, but the Q array won't.
 *   This is the vanilla version (without the tangential condition due to lambda)
 *   Use JoeWithModels::calcRhsCoupled()
 *
 * Arguments: rhs = the residual array
//	 *            qm1 = flow vector (previous step)
 *            qArray = flow vector
 *            delQ = increment in Q
 *            maxBacktrackingIter = maximum number of backtracking
 *            relaxationLowerBound =
 *            relaxationUpperBound =
 *            ResidOld = the residual array (previous step)
 *            ResidNew = the residual array (next step)
 *            whichNorm
 *            backtrackBarrierCoeff = barrier coefficient which is used to check if the residual increased or not
 *                                    (i.e. (1.0-backtrackBarrierCoeff*relaxation)*residNormTotOld)
 * Return: By arguments - The following 2 variables will be updated
 *                        1. rhs
 *                        2. ResidNew
 *         By return - final relaxation in step
 */
double IkeWithPsALC_AD::backtrackWithJOEcoupled_calcRelaxAndRHS(double* rhs, const double* qArray, const double* delQ,
		const int maxBacktrackingIter, const double relaxationLowerBound, const double relaxationUpperBound,
		const double* ResidOld, double* ResidNew, const int whichNorm, const double backtrackBarrierCoeff) {
	assert(relaxationUpperBound > relaxationLowerBound);

	int nVars = 5 + nScal;
	int N = ncv*nVars; // size of the flow vector (Q and delQ)

	double relax0 = 0.0, relax1 = 1.0, relax2 = relaxationUpperBound;
	double *ResidTemp = new double [nVars];
	double residNormTot0 = calcSumResidual(ResidOld, whichNorm);
	double residNormTot1 = calcSumResidual(ResidNew, whichNorm);
	double residNormTot2 = ABSURDLY_BIG_NUMBER;

	double **rhs2D; 	getMem2D(&rhs2D, 0, ncv-1, 0, 5+nScal-1, "IkeWithPsALC_AD::backtrackWithJOEcoupled_calcRelaxAndRHS->rhs2D", true);

	/*
	 * Initial Backtracking
	 *   In order to use the three-point parabolic model, you should launch the first backtracking anyway
	 */
	// Update flow field
	updateFlow_primVars(qArray, delQ, -relax2, nScal);

	// Calculate RHS and residual
	double ***Adummy = NULL;

	calcStateVariables(); 		// update state properties: velocity, pressure, temperature, enthalpy, gamma and R
	calcMaterialProperties(); 	// update material properties: laminar viscosity and heat conductivity
	setBC(); 					// set BC's for NS and scalars
	for(int ifa=nfa; ifa<nfa_b2; ++ifa)
		setBC1D(ifa, UgpWithCvCompFlow::rho, UgpWithCvCompFlow::rhou, UgpWithCvCompFlow::rhoE); // Note: ncvggf's cannot be not updated, so the user should take care of it
	calcRansTurbViscMuet(); 	// update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
	int countNegative = calcRhsCoupled(rhs2D, Adummy, nScal, false);
	convert2DrhsTo1Drhs(rhs, rhs2D, nScal);

	barrierSourceNS(rhs);      // Add barrier functions
	barrierSourceTurbScalars(rhs, nScal, iterNewton, residNormTotOld); // Add barrier functions

	calcResidualsFrom1Drhs(ResidTemp, rhs, whichNorm);
	residNormTot2 = calcSumResidual(ResidTemp, whichNorm);

	// show residual on screen
	if(debugLevel > 1 && mpi_rank == 0)
		printf("                  >> backtrackWithJOEcoupled_calcRelaxAndRHS(): MaxIter=%d, relaxLowBound=%.3e, relaxUpperBound=%.3e \n", maxBacktrackingIter, relaxationLowerBound, relaxationUpperBound);
	if(debugLevel > 1 && mpi_rank ==0)
		printf("                    >> Initial backtracking (RELAX=%.4e) -- RESID_TOT=%.8e \n", relax2, residNormTot2);

	/*
	 * Three-point parabolic backtracking algorithm
	 */
	int iter = 2; // note that you already run one backtracking before jumping into the following while loop

	bool needBacktracking = checkBacktrackFromResidInc(residNormTot1, residNormTot2, relax2, backtrackBarrierCoeff);

	while(countNegative>0 || needBacktracking) {
		/*
		 * check maximum iter
		 */
		if(iter > maxBacktrackingIter) {
			if(countNegative>0) {
				if(mpi_rank==0) {
					cerr<<"ERROR in backtrackWithJOEcoupled_calcRelaxAndRHS()! Backtracking can never go to the feasible region within iter="<<iter-1<<endl;
					throw(PSALC_ERROR_CODE);
				}
			} else {
				if(mpi_rank==0)
					cout<<"                  backtrackWithJOEcoupled_calcRelaxAndRHS(): WARNING!! Too many iterations for backtracking: iter="<<iter-1<<endl;
				break;
			}
		}

		/*
		 * calculate "relax": For details of the the three-point parabolic model, see C.T.Kelley
		 * Short descriptions:
		 *  First, construct a 2nd order polynomial p(relax) in terms of "relax" using the points of (relax0, resid0), (relax1, resid1), and (relax2, resid2).
		 *  (Note that in the initial iteration of the while() roop, relax0 = 0.0(previous step), relax1 = 1.0(no backtracking), relax2=relaxation (user specified argument)
		 *  Then, get "relax" for the smallest p(relax)
		 *  Finally, apply a safeguarding step.
		 */
		// get "relax" from the 2nd order polynomial
		double relaxNew;
		double pDD = 2*(relax2*(residNormTot1-residNormTot0)-relax1*(residNormTot2-residNormTot0)) / (relax1*relax2*(relax1-relax2));
		if(pDD <= 0.0) {
			relaxNew = relaxationUpperBound*relax2;
		} else {
			double pD = (relax1*relax1*(residNormTot2-residNormTot0) - relax2*relax2*(residNormTot1-residNormTot0)) / (relax1*relax2*(relax1-relax2));
			relaxNew = -pD/pDD;
		}
		// apply the safeguarding step
		if(relaxNew < relaxationLowerBound*relax2)
			relaxNew = relaxationLowerBound*relax2;
		else if (relaxNew > relaxationUpperBound*relax2)
			relaxNew = relaxationUpperBound*relax2;

		/*
		 * Update relax0, relax1, relax2, resid0, and resid1
		 */
		relax0 = relax1;
		relax1 = relax2;
		relax2 = relaxNew;
		residNormTot0 = residNormTot1;
		residNormTot1 = residNormTot2;

		/*
		 * Update the flow field
		 */
		updateFlow_primVars(qArray, delQ, -relax2, nScal);

		/*
		 * calculate the norm of the residual
		 */
		calcStateVariables(); 		// update state properties: velocity, pressure, temperature, enthalpy, gamma and R
		calcMaterialProperties(); 	// update material properties: laminar viscosity and heat conductivity
		setBC(); 					// set BC's for NS and scalars
		for(int ifa=nfa; ifa<nfa_b2; ++ifa)
			setBC1D(ifa, UgpWithCvCompFlow::rho, UgpWithCvCompFlow::rhou, UgpWithCvCompFlow::rhoE); // Note: ncvggf's cannot be not updated, so the user should take care of it
		calcRansTurbViscMuet(); 	// update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars
		countNegative = calcRhsCoupled(rhs2D, Adummy, nScal, false);
		convert2DrhsTo1Drhs(rhs, rhs2D, nScal);

		barrierSourceNS(rhs);      // Add barrier functions
		barrierSourceTurbScalars(rhs, nScal, iterNewton, residNormTotOld); // Add barrier functions

		calcResidualsFrom1Drhs(ResidTemp, rhs, whichNorm);
		residNormTot2 = calcSumResidual(ResidTemp, whichNorm);

		/*
		 * show residual on screen
		 */
		if(debugLevel > 1 && mpi_rank == 0) {
			switch(iter) {
			case 2:
				printf("                    >> %5dnd backtracking (RELAX=%.4e) -- RESID_TOT=%.8e\n", iter,relax2,residNormTot2);
				break;
			case 3:
				printf("                    >> %5drd backtracking (RELAX=%.4e) -- RESID_TOT=%.8e\n", iter,relax2,residNormTot2);
				break;
			default:
				printf("                    >> %5dth backtracking (RELAX=%.4e) -- RESID_TOT=%.8e\n", iter,relax2,residNormTot2);
				break;
			}
		}

		/*
		 * increase iter
		 */
		++iter;
	}
	if(debugLevel > 0 && mpi_rank == 0) {
		if(iter<=maxBacktrackingIter)
			printf("                  >> backtrackWithJOEcoupled_calcRelaxAndRHS(): Success! iter=%d  RELAX=%.4e  RESID_TOT=%.8e\n",iter-1,relax2,residNormTot2);
		else
			printf("                  >> backtrackWithJOEcoupled_calcRelaxAndRHS(): Unsatisfactory! iter=%d  RELAX=%.4e  RESID_TOT=%.8e\n",iter-1,relax2,residNormTot2);
	}

	/* update residNorm */
	for(int i=0; i<5+nScal; ++i)
		ResidNew[i] = ResidTemp[i];

	/* clear memory */
	delete [] ResidTemp;
	freeMem2D(rhs2D, 0, ncv-1, 0, 5+nScal-1); 	rhs2D = NULL;

	/* return */
	return relax2;
}

/*
 * Method: calcRhsWithBarrier
 * --------------------------
 * Calculate the single "rhs" array
 * Barrier functions are also called at the end of this method
 *
 * Return: countNegative = the number of negative-valued cells found during the RHS calculation
 */
int IkeWithPsALC_AD::calcRhsWithBarrier(double* rhs, const bool useBarrier) {
	/* Allocate memory */
	assert(RHSrho!=NULL && RHSrhou!=NULL && RHSrhoE!=NULL); // Note that RHSrho, RHSrhou, RHSrhoE, and RHSrhoScal were already allocated for tecplot output
	assert(RHSrhoScal!= NULL);
	double (*A)[5][5] = NULL;
	double ***AScal   = NULL;

	calcStateVariables(); 		// update state properties: velocity, pressure, temperature, enthalpy, gamma and R
	calcMaterialProperties(); 	// update material properties: laminar viscosity and heat conductivity
	setBC(); 					// set BC's for NS and scalars
	for(int ifa=nfa; ifa<nfa_b2; ++ifa)
		setBC1D(ifa, UgpWithCvCompFlow::rho, UgpWithCvCompFlow::rhou, UgpWithCvCompFlow::rhoE); // Note: ncvggf's cannot be not updated, so the user should take care of it
	calcRansTurbViscMuet(); 	// update turbulent properties: mut and initialize strain rate, vort mag, ... for turb scalars

	int countNegative = calcRhs(RHSrho, RHSrhou, RHSrhoE, RHSrhoScal, A, AScal, false);
	convertSeparatedRhsTo1Drhs(rhs, RHSrho, RHSrhou, RHSrhoE, RHSrhoScal, nScal);

	if(useBarrier) {
		barrierSourceNS(rhs);      // Add barrier functions
		barrierSourceTurbScalars(rhs, nScal, iterNewton, residNormTotOld); // Add barrier functions
	}

	for(int icv=0; icv<ncv; ++icv) {
		calcWeightRhs(icv, weightRhsMethod, nScal);
		calcWeightRhsHook(icv, nScal);
	}
	for(int icv=0; icv<ncv; ++icv) {
		int indexStart = icv*(5+nScal);

		rhs[indexStart] *= weightRhs[indexStart];

		for(int i=0; i<3; ++i)
			rhs[indexStart+1+i] *= weightRhs[indexStart+1+i];

		rhs[indexStart+4] *= weightRhs[indexStart+4];

		for (int iScal = 0; iScal < nScal; iScal++)
			rhs[indexStart+5+iScal] *= weightRhs[indexStart+5+iScal];
	}

	// Just for tecplot output: Since rhs is updated by barriers, RHSrho and etc. should be updated
	convert1DrhsToSeparatedRhs(rhs);
	if(nScal>0)
		convertRHSrhoScalToSeparatedRhs(RHSrhoScal);

//	double totBarrierMassSourceSumJOE;
//	double totBarrierEnergySourceSumJOE;
//	MPI_Allreduce(&myBarrierMassSourceSumJOE,   &totBarrierMassSourceSumJOE,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
//	MPI_Allreduce(&myBarrierEnergySourceSumJOE, &totBarrierEnergySourceSumJOE, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
//	if(debugLevel>0 && mpi_rank==0) {
//		printf("              >> Barrier: Total mass source = %.4e, Total energy source = %.4e \n", totBarrierMassSourceSumJOE, totBarrierEnergySourceSumJOE);
//	}

	/* Return */
	return countNegative;
}

/*
 * Method: calcBacktrackThreshold
 * ------------------------------
 * Calculate the thereshold for the backtracking method (the new residual norm should be smaller than this threshold).
 * This method will be called by checkBacktrackFromResidInc()
 */
double IkeWithPsALC_AD::calcBacktrackThreshold(const double residNormTotOld, const double relaxation, const double backtrackBarrierCoeff) {
	return (1.0-backtrackBarrierCoeff*relaxation)*residNormTotOld;
}

/*
 * Method: checkBacktrackFromResidInc
 * ----------------------------------
 * Check if there is any increment in the norm of the residual
 * For details, see C.T.Kelley, Iterative Methods for Linear and Nonlinear Equations, SIAM 1995, Chap.8.1, p.137
 *
 * Return: true if there is increment (you need backtracking in this case)
 *         false if the residual drops (you don't need backtracking in this case)
 */
bool IkeWithPsALC_AD::checkBacktrackFromResidInc(const double residNormTotOld, const double residNormTot, const double relaxation, const double backtrackBarrierCoeff) {
	assert(backtrackBarrierCoeff*relaxation>0.0 && backtrackBarrierCoeff*relaxation<1.0);

	double threshold = calcBacktrackThreshold(residNormTotOld, relaxation, backtrackBarrierCoeff);

	return !(residNormTot < threshold);
}

/*
 * Method: calcArclength
 * ---------------------
 * Calculate the arc-length
 * Original formulation: ds = q_tangent*(q1-q0)+weightLambda*lambda_tangent*(lambda1-lambda0),
 *                       where [q_tangent, lambda_tangent]^T is an unit vector
 * If [q_tangent, lambda_tangent]^T is not given, it can be easily calculated: q_tangent is equal to (q1-q0)/ds, and lambda_tangent is equal to (lambda1-lambda0)/ds.
 *
 * Return:
 *   By value     = arclength
 *   By reference = normSqRatioInArclength ( || q1-q0 || / || lambda1-lambda0 || )
 */
double IkeWithPsALC_AD::calcArclength(double* normSqRatioInArclength,
		const double* q0, const double* q1, const double* lambda0, const double* lambda1,
		const int NcontrolEqns, const int nScal, const double weightForLambda,
		double** q_tangent /*=NULL*/, const double* lambda_tangent /*=NULL*/) {
	if(q_tangent==NULL)
		assert(lambda_tangent==NULL);
	if(lambda_tangent==NULL)
		assert(q_tangent==NULL);
	assert(NcontrolEqns == 1); // Note: If NcontrolEqns!=1, you must modify the following lines !!
	                           //       (See calcNres(): it is not compatible with this method if NcontrolEqns!=1)

	int nVars = 5 + nScal;

	double arclength;
	if(q_tangent==NULL) {
		arclength = calcWeighted2NormForVecMinusVec(normSqRatioInArclength, q1, q0, lambda1, lambda0, weightForLambda, nVars, NcontrolEqns);
			// Note: arclength = the weighted 2-norm of [qVec1-qVec0; lambdaVec1-lambdaVec0]
	} else {
		// Allocate memory
		double* qVecDiff      = new double [ncv*nVars];
		double* lambdaVecDiff = NULL;
		if(mpi_rank == mpi_size-1)
			lambdaVecDiff = new double [NcontrolEqns];

		// Calculate the arclength
		for(int i=0; i<ncv*nVars; ++i)
			qVecDiff[i]	 = q1[i] - q0[i];
		if(mpi_rank == mpi_size-1)
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
				lambdaVecDiff[iEqn] = lambda1[iEqn] - lambda0[iEqn];

		arclength = vecDotVecWithWeight(normSqRatioInArclength, qVecDiff, q_tangent[0], lambdaVecDiff, lambda_tangent, weightForLambda, nVars, NcontrolEqns);
			// Note: arclength = the weighted inner product between [qVec1-qVec0; lambdaVec1-lambdaVec0] and [q_tangent[0]; lambda_tangent]

		// Free memory
		delete [] qVecDiff;
		if(mpi_rank == mpi_size-1)
			delete [] lambdaVecDiff;
	}

	return arclength;
}

double IkeWithPsALC_AD::calcArclength(const double* q0, const double* q1, const double* lambda0, const double* lambda1,
		const int NcontrolEqns, const int nScal, const double weightForLambda,
		double** q_tangent /*=NULL*/, const double* lambda_tangent /*=NULL*/) {
	return calcArclength(NULL, q0, q1, lambda0, lambda1, NcontrolEqns, nScal, weightForLambda, q_tangent, lambda_tangent);
}

/****************************
 * OVERLOADED FUNCTIONS
 ****************************/
/*
 * Method: UgpWithCvCompFlow_AD_init
 * ---------------------------------
 * Initialize some UgpWithCvCompFlow_AD variables
 */
void IkeWithPsALC_AD::UgpWithCvCompFlow_AD_init() {
	lambda = NULL;
	lambda_AD = NULL;
}

/*
 * Method: initialize_adoubles
 * ---------------------------
 * Original function = initialize_adoubles() in UgpWithCvCompFlowAD.h
 * System control parameters are added
 */
void IkeWithPsALC_AD::initialize_adoubles() {
	UgpWithCvCompFlow_AD::initialize_adoubles();

	if(NcontrolEqns>0)
		lambda_AD = new adouble [NcontrolEqns];
	else {
		if(mpi_rank==0 && debugLevel > 0)
			cout<<"WARNING IkeWithPsALC_AD::initialize_adoubles(): NcontrolEqns=0, so why are you calling this function?"<<endl;
	}
}

/*
 * Method: destroy_adoubles
 * ------------------------
 * Original function = destroy_adoubles() in UgpWithCvCompFlowAD.h
 * System control parameters are added
 */
void IkeWithPsALC_AD::destroy_adoubles() {
	UgpWithCvCompFlow_AD::destroy_adoubles();

	delete [] lambda_AD; 	lambda_AD = NULL;
}

#ifdef USE_MEM_SAVING_ADVAR_1D_
/*
 * Method: initialize_adoubles
 * ---------------------------
 * Original function = initialize_adoubles() in UgpWithCvCompFlowAD.h
 * System control parameters are added
 * Memory saving version with the ADvar classes
 */
void IkeWithPsALC_AD::initialize_adoubles(const int icvCenter, const int NcontrolEqns) {
	// If nbocv2_eachIcv, nbocv2ff_eachIcv, or fa2_eachIcv is empty, it must be an error.
	// However, we will perform the allocation here just to make this method run anyway.
	if(nbocv2_eachIcv.empty())
		find2layerCSRstruct_eachIcv(icvCenter, nbocv2_eachIcv, false);
	if(nbocv2ff_eachIcv.empty())
		find2layerCSRstruct_eachIcv(icvCenter, nbocv2ff_eachIcv, true);
	if(fa2_eachIcv.empty())
		find2LayerFaceIndices(icvCenter, fa2_eachIcv);

	rho_AD.allocate(nbocv2ff_eachIcv, icvCenter);
	rhou_AD.allocate(nbocv2ff_eachIcv, icvCenter);
	rhoE_AD.allocate(nbocv2ff_eachIcv, icvCenter);
	// Initialize as zero
	for(ADscalar<adouble>::iterator iter=rho_AD.begin(); iter!=rho_AD.end(); ++iter)
		*iter = 0.0;
	for(ADvector<adouble>::iterator iter=rhou_AD.begin(); iter!=rhou_AD.end(); ++iter) {
		for(int i=0; i<3; ++i)
			(*iter)[i] = 0.0;
	}
	for(ADscalar<adouble>::iterator iter=rhoE_AD.begin(); iter!=rhoE_AD.end(); ++iter)
		*iter = 0.0;

	vel.allocate(nbocv2ff_eachIcv, icvCenter);
	press.allocate(nbocv2ff_eachIcv, icvCenter);
	temp.allocate(nbocv2ff_eachIcv, icvCenter);
	enthalpy.allocate(nbocv2ff_eachIcv, icvCenter);
	// Initialize as zero
	for(ADvector<adouble>::iterator iter=vel.begin(); iter!=vel.end(); ++iter) {
		for(int i=0; i<3; ++i)
			(*iter)[i] = 0.0;
	}
	for(ADscalar<adouble>::iterator iter=press.begin(); iter!=press.end(); ++iter)
		*iter = 0.0;
	for(ADscalar<adouble>::iterator iter=temp.begin(); iter!=temp.end(); ++iter)
		*iter = 0.0;
	for(ADscalar<adouble>::iterator iter=enthalpy.begin(); iter!=enthalpy.end(); ++iter)
		*iter = 0.0;

	gamma.allocate(nbocv2ff_eachIcv, icvCenter);
	RoM.allocate(nbocv2ff_eachIcv, icvCenter);
	sos.allocate(nbocv2ff_eachIcv, icvCenter);
	// Initialize as zero
	for(ADscalar<adouble>::iterator iter=gamma.begin(); iter!=gamma.end(); ++iter)
		*iter = 0.0;
	for(ADscalar<adouble>::iterator iter=RoM.begin(); iter!=RoM.end(); ++iter)
		*iter = 0.0;
	for(ADscalar<adouble>::iterator iter=sos.begin(); iter!=sos.end(); ++iter)
		*iter = 0.0;

	strMag.allocate(nbocv2_eachIcv, icvCenter);
	vortMag.allocate(nbocv2_eachIcv, icvCenter);
	diverg.allocate(nbocv2_eachIcv, icvCenter);
	// Initialize as zero
	for(ADscalar<adouble>::iterator iter=strMag.begin(); iter!=strMag.end(); ++iter)
		*iter = 0.0;
	for(ADscalar<adouble>::iterator iter=vortMag.begin(); iter!=vortMag.end(); ++iter)
		*iter = 0.0;
	for(ADscalar<adouble>::iterator iter=diverg.begin(); iter!=diverg.end(); ++iter)
		*iter = 0.0;

	blendFuncF1_AD.allocate(nbocv2_eachIcv, icvCenter);
	blendFuncF2_AD.allocate(nbocv2_eachIcv, icvCenter);
	crossDiff_AD.allocate(nbocv2_eachIcv, icvCenter);
	// Initialize as zero
	for(ADscalar<adouble>::iterator iter=blendFuncF1_AD.begin(); iter!=blendFuncF1_AD.end(); ++iter)
		*iter = 0.0;
	for(ADscalar<adouble>::iterator iter=blendFuncF2_AD.begin(); iter!=blendFuncF2_AD.end(); ++iter)
		*iter = 0.0;
	for(ADscalar<adouble>::iterator iter=crossDiff_AD.begin(); iter!=crossDiff_AD.end(); ++iter)
		*iter = 0.0;

	CmeanSource_AD.allocate(nbocv2_eachIcv, icvCenter);
	muLam_AD.allocate(nbocv2_eachIcv, icvCenter);
	LambdaOverCp_AD.allocate(nbocv2_eachIcv, icvCenter);
	chi_AD.allocate(nbocv2_eachIcv, icvCenter);
	// Initialize as zero
	for(ADscalar<adouble>::iterator iter=CmeanSource_AD.begin(); iter!=CmeanSource_AD.end(); ++iter)
		*iter = 0.0;
	for(ADscalar<adouble>::iterator iter=muLam_AD.begin(); iter!=muLam_AD.end(); ++iter)
		*iter = 0.0;
	for(ADscalar<adouble>::iterator iter=LambdaOverCp_AD.begin(); iter!=LambdaOverCp_AD.end(); ++iter)
		*iter = 0.0;
	for(ADscalar<adouble>::iterator iter=chi_AD.begin(); iter!=chi_AD.end(); ++iter)
		*iter = 0.0;

	assert(kine == NULL);

#ifdef USE_ARTIF_VISC_WITH_MEM_SAVING
	artifViscMag_AD.allocate(nbocv2ff_eachIcv, icvCenter);
	// Initialize as zero
	for(ADscalar<adouble>::iterator iter=artifViscMag_AD.begin(); iter!=artifViscMag_AD.end(); ++iter)
		*iter = 0.0;
#endif

	// ----------------------------------------------------------------------------------------
	// init memory for face-based data
	// ----------------------------------------------------------------------------------------
	mul_fa.allocate(fa2_eachIcv, icvCenter);
	lamOcp_fa.allocate(fa2_eachIcv, icvCenter);
	mut_fa.allocate(fa2_eachIcv, icvCenter);
	// Initialize as zero
	for(ADscalar<adouble>::iterator iter=mul_fa.begin(); iter!=mul_fa.end(); ++iter)
		*iter = 0.0;
	for(ADscalar<adouble>::iterator iter=lamOcp_fa.begin(); iter!=lamOcp_fa.end(); ++iter)
		*iter = 0.0;
	for(ADscalar<adouble>::iterator iter=mut_fa.begin(); iter!=mut_fa.end(); ++iter)
		*iter = 0.0;

	// ----------------------------------------------------------------------------------------
	// gradients
	// ----------------------------------------------------------------------------------------
	if (sndOrder == true) {
		grad_rho.allocate(nbocv2_eachIcv, icvCenter);
		// Initialize as zero -- Note: In some parts of the RHS calculations, there is an assumption that gradients have been initialized as zero
		for(ADvector<adouble>::iterator iter=grad_rho.begin(); iter!=grad_rho.end(); ++iter) {
			for(int i=0; i<3; ++i)
				(*iter)[i] = 0.0;
		}
#ifdef temp_reconstruction
		grad_temp.allocate(nbocv2_eachIcv, icvCenter);
		// Initialize as zero -- Note: In some parts of the RHS calculations, there is an assumption that gradients have been initialized as zero
		for(ADvector<adouble>::iterator iter=grad_temp.begin(); iter!=grad_temp.end(); ++iter) {
			for(int i=0; i<3; ++i)
				(*iter)[i] = 0.0;
		}
#else
		grad_p.allocate(nbocv2_eachIcv, icvCenter);
		// Initialize as zero -- Note: In some parts of the RHS calculations, there is an assumption that gradients have been initialized as zero
		for(ADvector<adouble>::iterator iter=grad_p.begin(); iter!=grad_p.end(); ++iter) {
			for(int i=0; i<3; ++i)
				(*iter)[i] = 0.0;
		}
#endif
	}

	grad_u.allocate(nbocv2_eachIcv, icvCenter);         // allocate always
	// Initialize as zero -- Note: In some parts of the RHS calculations, there is an assumption that gradients have been initialized as zero
	for(ADtensor<adouble>::iterator iter=grad_u.begin(); iter!=grad_u.end(); ++iter) {
		for(int i=0; i<3; ++i)
			for(int j=0; j<3; ++j)
				(*iter)[i][j] = 0.0;
	}

	if ((sndOrder == true) || (mu_ref > 0.0)) {
		grad_enthalpy.allocate(nbocv2_eachIcv, icvCenter);
		// Initialize as zero -- Note: In some parts of the RHS calculations, there is an assumption that gradients have been initialized as zero
		for(ADvector<adouble>::iterator iter=grad_enthalpy.begin(); iter!=grad_enthalpy.end(); ++iter) {
			for(int i=0; i<3; ++i)
				(*iter)[i] = 0.0;
		}
	}

	int nScal = scalarTranspEqVector.size();

	if(nScal>0) {
		if(scalarTranspEqVector_AD == NULL)
			scalarTranspEqVector_AD = new ScalarTranspEq_AD[nScal];
	} else
		scalarTranspEqVector_AD = NULL;

	for(int i=0; i<nScal ; i++){
		scalarTranspEqVector_AD[i].phi.allocate(nbocv2ff_eachIcv, icvCenter);
		scalarTranspEqVector_AD[i].diff.allocate(fa2_eachIcv, icvCenter);
		scalarTranspEqVector_AD[i].rhophi.allocate(nbocv2ff_eachIcv, icvCenter);
		scalarTranspEqVector_AD[i].grad_phi.allocate(nbocv2_eachIcv, icvCenter);
		scalarTranspEqVector_AD[i].grad_rhophi.allocate(nbocv2_eachIcv, icvCenter);
		// Initialize as zero -- Note: In some parts of the RHS calculations, there is an assumption that gradients have been initialized as zero
		// Initialize as zero
		for(ADscalar<adouble>::iterator iter=scalarTranspEqVector_AD[i].phi.begin(); iter!=scalarTranspEqVector_AD[i].phi.end(); ++iter)
			*iter = 0.0;
		for(ADscalar<adouble>::iterator iter=scalarTranspEqVector_AD[i].diff.begin(); iter!=scalarTranspEqVector_AD[i].diff.end(); ++iter)
			*iter = 0.0;
		for(ADscalar<adouble>::iterator iter=scalarTranspEqVector_AD[i].rhophi.begin(); iter!=scalarTranspEqVector_AD[i].rhophi.end(); ++iter)
			*iter = 0.0;
		for(ADvector<adouble>::iterator iter=scalarTranspEqVector_AD[i].grad_phi.begin(); iter!=scalarTranspEqVector_AD[i].grad_phi.end(); ++iter) {
			for(int i=0; i<3; ++i)
				(*iter)[i] = 0.0;
		}
		for(ADvector<adouble>::iterator iter=scalarTranspEqVector_AD[i].grad_rhophi.begin(); iter!=scalarTranspEqVector_AD[i].grad_rhophi.end(); ++iter) {
			for(int i=0; i<3; ++i)
				(*iter)[i] = 0.0;
		}
//		assert(scalarTranspEqVector_AD[i].dpress_dphi.empty());
		strcpy(scalarTranspEqVector_AD[i].name, scalarTranspEqVector[i].getName());
	}

	set_adouble_names();

	if(NcontrolEqns>0) {
		if(lambda_AD == NULL)
			lambda_AD = new adouble [NcontrolEqns];
	}
}

/*
 * Method: destroy_adoubles
 * ------------------------
 * Original function = destroy_adoubles() in UgpWithCvCompFlowAD.h
 * System control parameters are added
 */
void IkeWithPsALC_AD::destroy_adoubles(const int NcontrolEqns) {
	if(scalarTranspEqVector_AD != NULL) {
		for(int i=0; i<nScal ; i++){
			scalarTranspEqVector_AD[i].phi.clear();
			scalarTranspEqVector_AD[i].diff.clear();
			scalarTranspEqVector_AD[i].rhophi.clear();
			scalarTranspEqVector_AD[i].grad_phi.clear();
			scalarTranspEqVector_AD[i].grad_rhophi.clear();
		}
	}

	if ((sndOrder == true) || (mu_ref > 0.0))
		grad_enthalpy.clear();

	grad_u.clear();

	if (sndOrder == true) {
#ifdef temp_reconstruction
		grad_temp.clear();
#else
		grad_p.clear();
#endif
		grad_rho.clear();
	}

	mut_fa.clear();
	lamOcp_fa.clear();
	mul_fa.clear();

	chi_AD.clear();
	LambdaOverCp_AD.clear();
	muLam_AD.clear();
	CmeanSource_AD.clear();

//	(*kine).clear(); 	kine = NULL;

	crossDiff_AD.clear();
	blendFuncF2_AD.clear();
	blendFuncF1_AD.clear();

	diverg.clear();
	vortMag.clear();
	strMag.clear();

	sos.clear();
	RoM.clear();
	gamma.clear();

	enthalpy.clear();
	temp.clear();
	press.clear();
	vel.clear();

	rhoE_AD.clear();
	rhou_AD.clear();
	rho_AD.clear();

#ifdef USE_ARTIF_VISC_WITH_MEM_SAVING
	artifViscMag_AD.clear();
#endif

	if(NcontrolEqns > 0) {
		delete [] lambda_AD; 	lambda_AD = NULL;
	}
}
#endif

/*
 * Method: print_tapestats
 * -----------------------
 * Original code: print_tapestats() in the JoeWithModels_AD class.
 */
void IkeWithPsALC_AD::print_tapestats(int tag) {
	int count[100];
	tapestats(tag, count);
	printf("Number of independents         = %6d\n", count[0]);
	printf("Number of dependents           = %6d\n", count[1]);
	printf("Max number of live active vars = %6d\n", count[2]);
	printf("Size of value stack            = %6d\n", count[3]);
	printf("Buffer size                    = %6d\n", count[4]);
	printf("Total number of operations     = %6d\n", count[5]);
	cout<<"======================================================"<<endl;
}

/*
 * Method: countMemoryDeficit_fromTapeStats
 * ----------------------------------------
 * Return the size of the memory deficit = the number of live active variables - the buffer size
 *                                         (Positive = memory is not enough).
 */
int IkeWithPsALC_AD::countMemoryDeficit_fromTapeStats(int tag) {
	int count[100];  // count[0] = number of independents,              count[1] = number of dependents,
	// count[2] = maximum number of linve active vars, count[3] = size of value stack,
	// count[4] = buffer size,                         count[5] = total number of operations
	// count[6-10] = other internal info
	tapestats(tag, count);

	return count[2] - count[4];
}

/****************************
 * HOOK FUNCTIONS
 ****************************/
/*
 * Method: initHookPsALC
 * ---------------------
 *
 */
void IkeWithPsALC_AD::initHookPsALC() {
	/* empty */
}

/*
 * Method: initHookPsALC_1st
 * -------------------------
 *
 */
void IkeWithPsALC_AD::initHookPsALC_1st() {
	/* empty */
}

/*
 * Method: initHookPsALC_2nd
 * -------------------------
 *
 */
void IkeWithPsALC_AD::initHookPsALC_2nd() {
	/* empty */
}

/*
 * Method: PsALCinitialHOOK_debug
 * ------------------------------
 * This function is for debugging purpose only
 */
void IkeWithPsALC_AD::PsALCinitialHOOK_debug() {
	/* empty */
}

/*
 * Method: initialHookNewton
 * -------------------------
 *
 */
void IkeWithPsALC_AD::initialHookNewton() {
	/* empty */
}

/*
 * Method: initialHookNewton_firstCall
 * -----------------------------------
 *
 */
void IkeWithPsALC_AD::initialHookNewton_firstCall() {
	/* empty */
}

/*
 * Method: temporalHookNewton
 * --------------------------
 * The user has an chance to apply an ad-hoc modification to the solution during the Newton iterations.
 * The user only needs to flow variables (rho, rhou, etc.), qVec, and lambda, but also provides the access to delQ and relaxation just in case.
 *
 * You must return TRUE if the solution is modified so that the Newton solver can update flow variables in the CPU boundaries and care other stuff.
 */
bool IkeWithPsALC_AD::temporalHookNewton(double *qVec, double *lambda, double *delQ, const double relaxation) {
	return false;
}

/*
 * Method: finalHookNewton
 * -----------------------
 *
 */
void IkeWithPsALC_AD::finalHookNewton() {
	/* empty */
}

/*
 * Method: checkBarrierParamNS
 * ---------------------------
 * Read parameters for N-S barrier functions and just return if the barrier functions will be used or not
 * Note: Wrapper function for readBarrierParamNS()
 */
bool IkeWithPsALC_AD::checkBarrierParamNS(const int iterNewton) {
	string dummyFunctionName;
	int    dummyMaxIter;
	double dummyThreshold_resid;
	double dummyCoeffRho;
	double dummyCoeffPress;

	return readBarrierParamNS(dummyFunctionName, dummyMaxIter, dummyThreshold_resid, dummyCoeffRho, dummyCoeffPress, iterNewton, residNormTotOld);
}

/*
 * Method: checkBarrierParamTurbScalars
 * ------------------------------------
 * Read parameters for turbuent scalar barrier functions and just return if the barrier functions will be used or not
 * Note: Wrapper function for readBarrierParamScalars()
 */
bool IkeWithPsALC_AD::checkBarrierParamTurbScalars(const int iterNewton) {
	string dummyFunctionName;
	int    dummyMaxIter;
	double dummyThreshold_resid;
    vector<double> dummyCoeffScalars;

	return readBarrierParamTurbScalars(dummyFunctionName, dummyMaxIter, dummyThreshold_resid, dummyCoeffScalars, iterNewton, residNormTotOld);
}

/*
 * Method: readBarrierParamNS
 * --------------------------
 * Read parameters for N-S barrier functions
 * Note: barrierSourceNS(), barrierSourceNS_AD(), and barrierSourceNS1D_AD() are calling this method
 * Return: true if the barreir functions will be used
 */
bool IkeWithPsALC_AD::readBarrierParamNS(string& barrierFunctionName, int& maxIter, double& threshold_resid, double& coeffRho, double& coeffPress,
		const int iterNewton, const double resid_norm_tot /* = ABSURDLY_BIG_NUMBER = 2.22e22 */ ) {
	static int iterNewtonSaved = -1; // readBarrierParamNS() can be called many times during single iterNewton (for example, called for ncv-times if the "1D-style" is used).
	                                 // Thus, all the messages will be also shown many times.
	                                 // This variable will be used to print-out messages just once for a iterNewton.
	static int iterNewtonBarrierStart = -1;

	if(iterNewton == 0) {
		iterNewtonSaved        = -1; // Reset to -1
		iterNewtonBarrierStart = -1; // Reset to -1
	}

	if (!checkParam("BARRIER_SOURCE_NS")) {
		ParamMap::add("BARRIER_SOURCE_NS  FUNCTION=NO_METHOD  MAX_ITER=0  THRESHOLD_RESID=1.0e-6  COEFF_RHO=1.0  COEFF_PRESS=1.0"); // add default values
		if (mpi_rank == 0)
			cout<< "WARNING: added keyword \"BARRIER_SOURCE_NS  METHOD=NO_METHOD  MAX_ITER=0  THRESHOLD_RESID=1.0e-6  COEFF_RHO=1.0  COEFF_PRESS=1.0\""<< " to parameter map!" << endl;
	}

	barrierFunctionName = getParam("BARRIER_SOURCE_NS")->getString("FUNCTION"); // If functionName=="NO_METHOD", skip the barrier

	if(strcmp(barrierFunctionName.c_str(),"NO_METHOD") != 0) {
		maxIter         = getParam("BARRIER_SOURCE_NS")->getInt("MAX_ITER"); // If maxIter==0, also skip the barrier
		threshold_resid = getParam("BARRIER_SOURCE_NS")->getDouble("THRESHOLD_RESID"); // If resid_norm_tot < threshold_resid, also skip the barrier
		coeffRho        = getParam("BARRIER_SOURCE_NS")->getDouble("COEFF_RHO");
		coeffPress      = getParam("BARRIER_SOURCE_NS")->getDouble("COEFF_PRESS");

		if(maxIter>0 && iterNewton<=maxIter && resid_norm_tot>=threshold_resid) {
			if (checkParam("BARRIER_SOURCE_NS_RAMP")) {  // Note: "BARRIER_SOURCE_NS_RAMP  AFTER_NEWTON_ITER=2  FACTOR_RESID=0.1  MIN_COEFF_RHO=1.0e-12  MIN_COEFF_PRESS=1.0e-12"
				int    afterNewtonIter = getParam("BARRIER_SOURCE_NS_RAMP")->getInt("AFTER_NEWTON_ITER");
				double threshResid     = getParam("BARRIER_SOURCE_NS_RAMP")->getDouble("THRESH_RESID");
				double factorResid     = getParam("BARRIER_SOURCE_NS_RAMP")->getDouble("FACTOR_RESID");
				double minCoeffRho     = getParam("BARRIER_SOURCE_NS_RAMP")->getDouble("MIN_COEFF_RHO");
				double minCoeffPress   = getParam("BARRIER_SOURCE_NS_RAMP")->getDouble("MIN_COEFF_PRESS");

				if((factorResid<=0.0 || factorResid>1.0) && mpi_rank==0)
					if(iterNewtonSaved != iterNewton)
						cout<<"WARNING in IkeWithPsALC_AD::readBarrierParamNS(): BARRIER_SOURCE_NS_RAMP-->FACTOR_RESID is "<<factorResid<<endl;

				if(iterNewton>afterNewtonIter || resid_norm_tot<threshResid) {
					if(iterNewtonBarrierStart == -1)
						iterNewtonBarrierStart = iterNewton-1;

					if(iterNewtonBarrierStart >= iterNewton) {
						if(mpi_rank==0)
							cout<<"WARNING in readBarrierParamNS(): iterNewtonBarrierStart="<<iterNewtonBarrierStart<<" is larger than iterNewton="<<iterNewton<<"."<<endl
							    <<"                                 This can cause NaN"<<endl;
					}

					coeffRho   = max(coeffRho  *pow(factorResid, (double) iterNewton-iterNewtonBarrierStart),  minCoeffRho  );
					coeffPress = max(coeffPress*pow(factorResid, (double) iterNewton-iterNewtonBarrierStart),  minCoeffPress);

					if(coeffRho>1.0 || coeffPress>1.0) {
						if(mpi_rank==0)
							cout<<"WARNING in readBarrierParamNS(): The barrier coefficients become too large after ramping..."<<endl
							    <<"                                 Ramp Parameters - AFTER_NEWTON_ITER="<<afterNewtonIter<<", THRESH_RESID="<<threshResid<<endl
							    <<"                                                   FACTOR_RESID="<<factorResid<<", MIN_COEFF_RHO="<<minCoeffRho<<", MIN_COEFF_PRESS="<<minCoeffPress<<endl;
					}

					if(debugLevel>0 && mpi_rank==0)
						if(iterNewtonSaved != iterNewton)
                            if(coeffRho>minCoeffRho || coeffPress>minCoeffPress)
    							printf("                         Reduce NS BARRIER COEFF: COEFF_RHO=%.5e, COEFF_PRESS=%.5e\n", coeffRho, coeffPress);
				}
			}

			iterNewtonSaved = iterNewton;
			return true;
		} else {
			iterNewtonSaved = iterNewton;
			return false;
		}
	}

	iterNewtonSaved = iterNewton;
	return false;
}

// Note: readBarrierParamTurbScalars() is in IkeUgpWithCvCompFlow.h

/*
 * Method: barrierSourceNS
 * -----------------------
 * Add barrier functions to the RHS of the N-S equations
 */
void IkeWithPsALC_AD::barrierSourceNS(double* rhs) {
	string functionName;    // If functionName=="NO_METHOD",            skip the barrier
	int    maxIter;         // If maxIter==0,                           skip the barrier
	double threshold_resid; // If resid_norm_tot_old < threshold_resid, skip the barrier
	double coeffRho;
	double coeffPress;
	bool useBarrier = readBarrierParamNS(functionName, maxIter, threshold_resid, coeffRho, coeffPress, iterNewton, residNormTotOld);
	// Note: iterNewton and residNormTotOld are member variables of the IkeWithPsALC_AD class

	if(useBarrier) {
		int nVars = 5+nScal;

		double *kineArray = NULL;
		int kine_Index = getScalarTransportIndex("kine");
		if(kine_Index>=0)
			kineArray = scalarTranspEqVector[kine_Index].phi;

		myBarrierMassSourceSumJOE   = 0.0;
		myBarrierEnergySourceSumJOE = 0.0;

		if (strcmp(functionName.c_str(), "LOG") == 0) {
			for (int icv=0; icv<ncv; ++icv) {
				// source term in the mass equation
				double massSource = max(-coeffRho*log(fabs(rho[icv]/rho_ref)), 0.0)*cv_volume[icv];
				rhs[icv*nVars] += massSource;

				// source term in the energy equation
				double kinecv = 0.0;
				if(kine_Index>=0)
					kinecv = kineArray[icv];
				double tempPress = calcPress(UgpWithCvCompFlow::gamma[icv], rhoE[icv], rhou[icv], rho[icv], kinecv);
				double energySource = max(-coeffPress*log(fabs(tempPress/p_ref)), 0.0)*cv_volume[icv];
				rhs[icv*nVars+4] += energySource;

				// statistics
				myBarrierMassSourceSumJOE   += massSource;
				myBarrierEnergySourceSumJOE += energySource;
			}
		} else if (strcmp(functionName.c_str(), "RECIPROCAL") == 0) {
			for (int icv=0; icv<ncv; ++icv) {
				// source term in the mass equation
				double massSource = coeffRho/(fabs(rho[icv]/rho_ref))*cv_volume[icv];
				rhs[icv*nVars] += massSource;

				// source term in the energy equation
				double kinecv = 0.0;
				if(kine_Index>=0)
					kinecv = kineArray[icv];
				double tempPress = calcPress(UgpWithCvCompFlow::gamma[icv], rhoE[icv], rhou[icv], rho[icv], kinecv);
				double energySource = coeffPress/(fabs(tempPress/p_ref))*cv_volume[icv];
				rhs[icv*nVars+4] += energySource;

				// statistics
				myBarrierMassSourceSumJOE   += massSource;
				myBarrierEnergySourceSumJOE += energySource;
			}
		} else {
			if(debugLevel>1 && mpi_rank==0)
				cout<<"barrierSourceNS() is not active: No barrier method"<<endl;
		}
	} else {
		if(debugLevel>1 && mpi_rank==0)
			cout<<"IkeWithPsALC_AD::barrierSourceNS(): Barrier function = NO_METHOD"<<endl;
	}
}

/*
 * Method: barrierSourceNS
 * -----------------------
 * Add barrier functions to the RHS of the N-S equations for JOE-type RHS vectors (with Jacobian matrix update)
 */
void IkeWithPsALC_AD::barrierSourceNS(double *RHSrho, double (*RHSrhou)[3], double *RHSrhoE, double **rhsScal, double (*A)[5][5], double ***AScal, int flagImplicit) {
	string functionName;    // If functionName=="NO_METHOD",            skip the barrier
	int    maxIter;         // If maxIter==0,                           skip the barrier
	double threshold_resid; // If resid_norm_tot_old < threshold_resid, skip the barrier
	double coeffRho;
	double coeffPress;
	bool useBarrier = readBarrierParamNS(functionName, maxIter, threshold_resid, coeffRho, coeffPress, iterNewton, residNormTotOld);
	// Note: iterNewton and residNormTotOld are member variables of the IkeWithPsALC_AD class

	if(useBarrier) {
		int nVars = 5+nScal;

		double *kineArray = NULL;
		int kine_Index = getScalarTransportIndex("kine");
		if(kine_Index>=0)
			kineArray = scalarTranspEqVector[kine_Index].phi;

		myBarrierMassSourceSumJOE   = 0.0;
		myBarrierEnergySourceSumJOE = 0.0;

		if (strcmp(functionName.c_str(), "LOG") == 0) {
			for (int icv=0; icv<ncv; ++icv) {
				// source term in the mass equation
				double massSource = max(-coeffRho*log(fabs(rho[icv]/rho_ref)), 0.0)*cv_volume[icv];
				RHSrho[icv] += massSource;

				// Jacobian matrix for the mass conservation equation
				if(-coeffRho*log(fabs(rho[icv]/rho_ref)) > 0.0)
					A[icv][0][0] += -coeffRho / rho[icv];

				// source term in the energy equation
				double kinecv = 0.0;
				if(kine_Index>=0)
					kinecv = kineArray[icv];
				double tempPress = calcPress(UgpWithCvCompFlow::gamma[icv], rhoE[icv], rhou[icv], rho[icv], kinecv);
				double energySource = max(-coeffPress*log(fabs(tempPress/p_ref)), 0.0)*cv_volume[icv];
				RHSrhoE[icv] += energySource;

				// Jacobian matrix for the energy equation
				if(-coeffPress*log(fabs(tempPress/p_ref)) > 0.0) {
					// Note: Even though pressure can be calculated in a different way than calcPress(),
					//       here I assume calcPress() for simplicity.
					double coeffOverPressTimesGammaMinusOne = coeffPress/UgpWithCvCompFlow::press[icv] * (UgpWithCvCompFlow::gamma[icv]-1.0);
					double rhouSq = 0.0;
					for(int i=0; i<3; ++i)
						rhouSq += rhou[icv][i]*rhou[icv][i];

					A[icv][4][0] += -coeffOverPressTimesGammaMinusOne * (0.5*rhouSq/rho[icv]/rho[icv] - kinecv) *cv_volume[icv];
					for(int i=1; i<4; ++i)
						A[icv][4][i] += coeffOverPressTimesGammaMinusOne * (rhou[icv][i]/rho[icv]) *cv_volume[icv];
					A[icv][4][4] += -coeffOverPressTimesGammaMinusOne *cv_volume[icv];
				}

				// statistics
				myBarrierMassSourceSumJOE   += massSource;
				myBarrierEnergySourceSumJOE += energySource;
			}
		} else if (strcmp(functionName.c_str(), "RECIPROCAL") == 0) {
			for (int icv=0; icv<ncv; ++icv) {
				// source term in the mass equation
				double massSource = coeffRho/(fabs(rho[icv]/rho_ref))*cv_volume[icv];
				RHSrho[icv] += massSource;

				// Jacobian matrix for the mass conservation equation
				double signRho = (rho[icv]>=0.0) ? 1.0 : -1.0;
				A[icv][0][0] += -signRho*coeffRho / pow(rho[icv], 2.0);

				// source term in the energy equation
				double kinecv = 0.0;
				if(kine_Index>=0)
					kinecv = kineArray[icv];
				double tempPress = calcPress(UgpWithCvCompFlow::gamma[icv], rhoE[icv], rhou[icv], rho[icv], kinecv);
				double energySource = coeffPress/(fabs(tempPress/p_ref))*cv_volume[icv];
				RHSrhoE[icv] += energySource;

				// Jacobian matrix for the energy equation
				double signPress = (tempPress>=0.0) ? 1.0 : -1.0;
				// Note: Even though pressure can be calculated in a different way than calcPress(),
				//       here I assume calcPress() for simplicity.
				double signedCoeffOverPressSqTimesGammaMinusOne = signPress * coeffPress/pow(UgpWithCvCompFlow::press[icv], 2.0) * (UgpWithCvCompFlow::gamma[icv]-1.0);
				double rhouSq = 0.0;
				for(int i=0; i<3; ++i)
					rhouSq += rhou[icv][i]*rhou[icv][i];

				A[icv][4][0] += -signedCoeffOverPressSqTimesGammaMinusOne * (0.5*rhouSq/rho[icv]/rho[icv] - kinecv) *cv_volume[icv];
				for(int i=1; i<4; ++i)
					A[icv][4][i] += signedCoeffOverPressSqTimesGammaMinusOne * (rhou[icv][i]/rho[icv]) *cv_volume[icv];
				A[icv][4][4] += -signedCoeffOverPressSqTimesGammaMinusOne *cv_volume[icv];

				// statistics
				myBarrierMassSourceSumJOE   += massSource;
				myBarrierEnergySourceSumJOE += energySource;
			}
		} else {
			if(debugLevel>1 && mpi_rank==0)
				cout<<"barrierSourceNS() is not active: No barrier method"<<endl;
		}
	} else {
		if(debugLevel>1 && mpi_rank==0)
			cout<<"IkeWithPsALC_AD::barrierSourceNS(): Barrier function = NO_METHOD"<<endl;
	}
}

// Note: barrierSourceTurbScalars() is in IkeUgpWithCvCompFlow.h

/*
 * Method: barrierSourceNS_AD
 * --------------------------
 * Add barrier functions to the RHS of the N-S equations
 */
void IkeWithPsALC_AD::barrierSourceNS_AD(REALA *rhs_rho_AD, REALA (*rhs_rhou_AD)[3], REALA *rhs_rhoE_AD) {
	string functionName; // If functionName=="NO_METHOD", skip the barrier
	int    maxIter;      // If maxIter==0,                skip the barrier
	double threshold_resid; // If resid_norm_tot_old < threshold_resid, skip the barrier
	double coeffRho;
	double coeffPress;
	bool useBarrier = readBarrierParamNS(functionName, maxIter, threshold_resid, coeffRho, coeffPress, iterNewton, residNormTotOld);
	// Note: iterNewton and residNormTotOld are member variables of the IkeWithPsALC_AD class

	if(useBarrier) {
		int nVars = 5+nScal;

		if (strcmp(functionName.c_str(), "LOG") == 0) {
			for (int icv=0; icv<ncv; ++icv) {
				// source term in the mass equation
				adouble massSource = max(-coeffRho*log(fabs(rho_AD[icv]/rho_ref)), 0.0)*cv_volume[icv];
				rhs_rho_AD[icv] += massSource;

				// source term in the energy equation
				adouble kinecv_AD = 0.0;
#ifdef USE_MEM_SAVING_ADVAR_1D_
				if (kine != NULL)
					kinecv_AD = (*kine)[icv]; // Note: here, "kine" is a pointer to a class
#else
				if (kine != NULL)
					kinecv_AD = kine[icv]; // Note: here, "kine" is an adouble array
#endif
				adouble tempPress = calcPress(gamma[icv], rhoE_AD[icv], rhou_AD[icv], rho_AD[icv], kinecv_AD);
				adouble energySource = max(-coeffPress*log(fabs(tempPress/p_ref)), 0.0)*cv_volume[icv];
				rhs_rhoE_AD[icv] += energySource;
			}
		} else if (strcmp(functionName.c_str(), "RECIPROCAL") == 0) {
			for (int icv=0; icv<ncv; ++icv) {
				// source term in the mass equation
				adouble massSource = coeffRho/(fabs(rho_AD[icv]/rho_ref))*cv_volume[icv];
				rhs_rho_AD[icv*nVars] += massSource;

				// source term in the energy equation
				adouble kinecv_AD = 0.0;
#ifdef USE_MEM_SAVING_ADVAR_1D_
				if (kine != NULL)
					kinecv_AD = (*kine)[icv]; // Note: here, "kine" is an adouble array
#else
				if (kine != NULL)
					kinecv_AD = kine[icv]; // Note: here, "kine" is an adouble array
#endif
				adouble tempPress = calcPress(gamma[icv], rhoE_AD[icv], rhou_AD[icv], rho_AD[icv], kinecv_AD);
				adouble energySource = coeffPress/(fabs(tempPress/p_ref))*cv_volume[icv];
				rhs_rhoE_AD[icv] += energySource;
			}
		} else {
			if(debugLevel>1 && mpi_rank==0)
				cout<<"barrierSourceNS_AD() is empty: No barrier method"<<endl;
		}
	} else {
		if(debugLevel>1 && mpi_rank==0)
			cout<<"IkeWithPsALC_AD::barrierSourceNS_AD(): Barrier function = NO_METHOD"<<endl;
	}
}

// Note: barrierSourceTurbScalars_AD() is in IkeUgpWithCvCompFlow.h

/*
 * Method: barrierSourceNS1D_AD
 * ----------------------------
 * Add barrier functions to the RHS of the N-S equations
 */
void IkeWithPsALC_AD::barrierSourceNS1D_AD(const int icvCenter, REALA &rhs_rho_AD, REALA rhs_rhou_AD[3], REALA &rhs_rhoE_AD) {
	string functionName; // If functionName=="NO_METHOD", skip the barrier
	int    maxIter;      // If maxIter==0,                skip the barrier
	double threshold_resid; // If resid_norm_tot_old < threshold_resid, skip the barrier
	double coeffRho;
	double coeffPress;
	bool useBarrier = readBarrierParamNS(functionName, maxIter, threshold_resid, coeffRho, coeffPress, iterNewton, residNormTotOld);
	// Note: iterNewton and residNormTotOld are member variables of the IkeWithPsALC_AD class

	if(useBarrier) {
		int nVars = 5+nScal;

		if (strcmp(functionName.c_str(), "LOG") == 0) {
			// source term in the mass equation
			adouble massSource = max(-coeffRho*log(fabs(rho_AD[icvCenter]/rho_ref)), 0.0)*cv_volume[icvCenter];
			rhs_rho_AD += massSource;
			myBarrierMassSourceSum1D_AD += massSource.value(); // Just to check how much mass is added to the system
			// Note: You must reinitialize this variable to zero before calling this method

			// source term in the energy equation
			adouble kinecv_AD = 0.0;
#ifdef USE_MEM_SAVING_ADVAR_1D_
			if (kine != NULL)
				kinecv_AD = (*kine)[icvCenter]; // Note: here, "kine" is an adouble array
#else
			if (kine != NULL)
				kinecv_AD = kine[icvCenter]; // Note: here, "kine" is an adouble array
#endif
			adouble tempPress = calcPress(gamma[icvCenter], rhoE_AD[icvCenter], rhou_AD[icvCenter], rho_AD[icvCenter], kinecv_AD);
			adouble energySource = max(-coeffPress*log(fabs(tempPress/p_ref)), 0.0)*cv_volume[icvCenter];
			rhs_rhoE_AD += energySource;
			myBarrierEnergySourceSum1D_AD += energySource.value(); // Just to check how much energy is added to the system
			// Note: You must reinitialize this variable to zero before calling this method
		} else if (strcmp(functionName.c_str(), "RECIPROCAL") == 0) {
			// source term in the mass equation
			adouble massSource = coeffRho/(fabs(rho_AD[icvCenter]/rho_ref))*cv_volume[icvCenter];
			rhs_rho_AD += massSource;
			myBarrierMassSourceSum1D_AD += massSource.value(); // Just to check how much mass is added to the system
			// Note: You must reinitialize this variable to zero before calling this method

			// source term in the energy equation
			adouble kinecv_AD = 0.0;
#ifdef USE_MEM_SAVING_ADVAR_1D_
			if (kine != NULL)
				kinecv_AD = (*kine)[icvCenter]; // Note: here, "kine" is an adouble array
#else
			if (kine != NULL)
				kinecv_AD = kine[icvCenter]; // Note: here, "kine" is an adouble array
#endif
			adouble tempPress = calcPress(gamma[icvCenter], rhoE_AD[icvCenter], rhou_AD[icvCenter], rho_AD[icvCenter], kinecv_AD);
			adouble energySource = coeffPress/(fabs(tempPress/p_ref))*cv_volume[icvCenter];
			rhs_rhoE_AD += energySource;
			myBarrierEnergySourceSum1D_AD += energySource.value(); // Just to check how much energy is added to the system
			// Note: You must reinitialize this variable to zero before calling this method
		} else {
			if(debugLevel>1 && mpi_rank==0)
				cout<<"barrierSourceNS1D_AD() is empty: No barrier method"<<endl;
		}
	} else {
		if(debugLevel>1 && mpi_rank==0 && icvCenter==0)
			cout<<"IkeWithPsALC_AD::barrierSourceNS1D_AD(): Barrier function = NO_METHOD"<<endl;
	}
}

// Note: barrierSourceTurbScalars1D_AD() is in IkeUgpWithCvCompFlow.h


/*
 * Method: calcWeightRhs
 * ---------------------
 * calculate "weightRhs" (weight for RHS) at the given icv
 */
void IkeWithPsALC_AD::calcWeightRhs(const int icvCenter, WEIGHT_RHS_METHOD weightRhsMethod, const int nScal) {
	int startingIndex = icvCenter * (5+nScal);
	static bool firstCall = true;

	if(weightRhsMethod == NO_RHSWEIGHT) {
		for(int i=0; i<5+nScal; ++i)
			weightRhs[startingIndex+i] = 1.0;
	} else if(weightRhsMethod == REF_VALUES) {
		weightRhs[startingIndex] = max( 1.0/(ADDITIONAL_SCALING_VALUE*RefFlowParams.rho_ref),  WEIGHT_RHS_MIN );
		for(int i=0; i<3; ++i)
//			weightRhs[startingIndex+1+i] = max( 1.0/(ADDITIONAL_SCALING_VALUE*max(RefFlowParams.rhou_ref[i], 0.02*RefFlowParams.rhouMag_ref)),  WEIGHT_RHS_MIN );
			weightRhs[startingIndex+1+i] = max( 1.0/(ADDITIONAL_SCALING_VALUE*RefFlowParams.rhouMag_ref),  WEIGHT_RHS_MIN );
		weightRhs[startingIndex+4] = max( 1.0/(ADDITIONAL_SCALING_VALUE*RefFlowParams.rhoE_ref),  WEIGHT_RHS_MIN );
		for(int iScal=0; iScal<nScal; ++iScal)
			weightRhs[startingIndex+5+iScal] = max( 1.0/(ADDITIONAL_SCALING_VALUE*RefFlowParams.scalar_ref[iScal]),  WEIGHT_RHS_MIN);

		if(mpi_rank==0 && firstCall) {
			cout<<"           > IkeWithPsALC_AD::calcWeightRhs(): REF_VALUES weighting =";
			for(int i=0; i<5+nScal; ++i) cout<<" "<<std::setprecision(3)<<weightRhs[startingIndex+i]<<" ";
			cout<<endl;
		}
	} else if(weightRhsMethod == LOCAL_VALUES) {
		weightRhs[startingIndex] = max( 1.0/(ADDITIONAL_SCALING_VALUE*fabs(rho[icvCenter])),  WEIGHT_RHS_MIN );
		double rhouMag = sqrt(vecDotVec3d(rhou[icvCenter], rhou[icvCenter]));
		for(int i=0; i<3; ++i)
//			weightRhs[startingIndex+1+i] = max( 1.0/(ADDITIONAL_SCALING_VALUE*max(0.02, fabs(rhou[icvCenter][i])/rhouMag)*rhouMag),  WEIGHT_RHS_MIN );
			weightRhs[startingIndex+1+i] = max( 1.0/(ADDITIONAL_SCALING_VALUE*rhouMag),  WEIGHT_RHS_MIN );
		weightRhs[startingIndex+4] = max( 1.0/(ADDITIONAL_SCALING_VALUE*fabs(rhoE[icvCenter])),  WEIGHT_RHS_MIN );
		for(int iScal=0; iScal<nScal; ++iScal)
			weightRhs[startingIndex+5+iScal] = max( 1.0/(ADDITIONAL_SCALING_VALUE*max(fabs(scalarTranspEqVector[iScal].phi[icvCenter]), 0.02*RefFlowParams.scalar_ref[iScal])),  WEIGHT_RHS_MIN);

		if(mpi_rank==0 && firstCall) {
			cout<<"           > IkeWithPsALC_AD::calcWeightRhs(): LOCAL_VALUES weighting (at icv="<<icvCenter<<") =";
			for(int i=0; i<5+nScal; ++i) cout<<" "<<weightRhs[startingIndex+i]<<" ";
			cout<<endl;
		}
	} else {
		if(mpi_rank==0) {
			cerr<<"ERROR IkeWithPsALC_AD::calcWeightRhs() - weightRhsMethod = "<<weightRhsMethod<<" cannot be supported"<<endl;
			cerr<<"      Supported methods: "<<NO_RHSWEIGHT<<" = NO_RHSWEIGHT"<<endl
			    <<"      Supported methods: "<<REF_VALUES<<" = REF_VALUES"<<endl
			    <<"      Supported methods: "<<LOCAL_VALUES<<" = LOCAL_VALUES"<<endl;
		}
		throw(PSALC_ERROR_CODE);
	}

	firstCall = false;
}

/*
 * Method: calcWeightRhsHook
 * -------------------------
 * calculate "weightRhs" (weight for RHS) at the given icv
 * This method is called right after calling calcWeightRhs() so that the user can modify weightRhs
 */
void IkeWithPsALC_AD::calcWeightRhsHook(const int icvCenter, const int nVars) { /* empty */ }

/****************************
 * UTILITY FUNCTIONS
 ****************************/
/*
 * Method: getNewtonParam
 * ----------------------
 * Read parameters for Newton's method from the input file
 */
void IkeWithPsALC_AD::getNewtonParam() {
	static bool firstCall = true;

	// Linear solver setting (PETSc)
	//		 *   1. Linear solver thresholds
	//		 *   	maxIterLS = maximum number of Newton (outer) iterations
	//		 *   	zeroAbsLS = absolute residual of Newton (outer) iterations
	//		 *   	zeroRelLS = relative residual of Newton (outer) iterations
	//		 *   2. Ramp
	//		 * 		startDecLS     = After "startDecLS" iteration, ramp starts
	//		 * 		intervalDecLS  = At every "intervalDecLS" iteration, ramp occurs
	//		 * 		incIterLS      = When ramps occurs, "maxIterLS" increases by "incIterLS"
	//		 * 		maxFinalIterLS = After "maxFinalIterLS" iteration, ramp no longer occurs
	//		 * 		decZeroLS      = When ramps occurs, "zeroAbsLS" and "zeroRelLS" drops by "decZeroLS"
	//		 * 		minZeroAbsLS   =
	//		 * 		minZeroRelLS   =
	if (!checkParam("LINEAR_SOLVER_NEWTON_TRESHOLDS")) {
		ParamMap::add("LINEAR_SOLVER_NEWTON_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-3"); // add default values
		if (mpi_rank == 0)
			cout<< "WARNING: added keyword \"LINEAR_SOLVER_NEWTON_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-3\""<< " to parameter map!" << endl;
	}
	newtonParam.maxIterLS = getParam("LINEAR_SOLVER_NEWTON_TRESHOLDS")->getInt("MAX_ITER");
	newtonParam.zeroAbsLS = getParam("LINEAR_SOLVER_NEWTON_TRESHOLDS")->getDouble("ABS_RESID");
	newtonParam.zeroRelLS = getParam("LINEAR_SOLVER_NEWTON_TRESHOLDS")->getDouble("REL_RESID");
	if(firstCall && mpi_rank==0)
		printf("> LINEAR_SOLVER_NEWTON_TRESHOLDS  MAX_ITER=%d  ABS_RESID=%.2e  REL_RESID=%.2e\n", newtonParam.maxIterLS, newtonParam.zeroAbsLS, newtonParam.zeroRelLS);

	// Ramp maximum iterations and residual
	if (!checkParam("LSNT_RAMP")) {
		ParamMap::add("LSNT_RAMP  AFTER_NEWTON_ITER=10  INTERVAL_NEWTON_ITER=1  FACTOR_ITER=0  MAX_ITER=30  FACTOR_RESID=1.0  MIN_ABS_RESID=1.0e-8  MIN_REL_RESID=1.0e-3"); // add default values
		if (mpi_rank == 0)
			cout<< "WARNING: added keyword \"LSNT_RAMP  AFTER_NEWTON_ITER=10  INTERVAL_NEWTON_ITER=1  FACTOR_ITER=0  MAX_ITER=30  FACTOR_RESID=1.0  MIN_ABS_RESID=1.0e-8  MIN_REL_RESID=1.0e-3\""<< " to parameter map!" << endl;
	}

	newtonParam.startDecLS = getParam("LSNT_RAMP")->getInt("AFTER_NEWTON_ITER");
	newtonParam.intervalDecLS = getParam("LSNT_RAMP")->getInt("INTERVAL_NEWTON_ITER");
	newtonParam.incIterLS = getParam("LSNT_RAMP")->getInt("FACTOR_ITER");
	newtonParam.maxFinalIterLS = getParam("LSNT_RAMP")->getInt("MAX_ITER");
	newtonParam.decZeroLS = getParam("LSNT_RAMP")->getDouble("FACTOR_RESID");
	newtonParam.minZeroAbsLS = getParam("LSNT_RAMP")->getDouble("MIN_ABS_RESID");
	newtonParam.minZeroRelLS = getParam("LSNT_RAMP")->getDouble("MIN_REL_RESID");
	if(firstCall && mpi_rank==0)
		printf("> LSNT_RAMP  AFTER_NEWTON_ITER=%d  INTERVAL_NEWTON_ITER=%d  FACTOR_ITER=%d  MAX_ITER=%d  FACTOR_RESID=%.2e  MIN_ABS_RESID=%.2e  MIN_REL_RESID=%.2e\n",
				newtonParam.startDecLS, newtonParam.intervalDecLS, newtonParam.incIterLS, newtonParam.maxIterLS,
				newtonParam.decZeroLS, newtonParam.minZeroAbsLS, newtonParam.minZeroRelLS);
	if(newtonParam.incIterLS<0 && mpi_rank==0)
		printf("  WARNING: FACTOR_ITER is a negative number = %d \n", newtonParam.incIterLS);
	if((newtonParam.decZeroLS>1.0 || newtonParam.decZeroLS<=0.0) && mpi_rank==0)
		printf("  WARNING: FACTOR_RESID is out of range (should be 0.0<FACTOR_RESID<=1.0) = %.2e \n", newtonParam.decZeroLS);


	// Relaxation reducing based on negative values
	if (!checkParam("RELAXATION_CLIP_FOR_NEGATIVE_VALS")) {
		ParamMap::add("RELAXATION_CLIP_FOR_NEGATIVE_VALS  CLIP_PARAM=1.0e-9  SAFE_PARAM=0.8"); // add default values
		if (mpi_rank == 0)
			cout << "WARNING: added keyword \"RELAXATION_CLIP_NEGATIVE_VALS  CLIP_PARAM=1.0e-9  SAFE_PARAM=0.8\" to parameter map!"<<endl;
	}

	newtonParam.clipParameter  = getParam("RELAXATION_CLIP_FOR_NEGATIVE_VALS")->getDouble("CLIP_PARAM");
	newtonParam.safeParameter  = getParam("RELAXATION_CLIP_FOR_NEGATIVE_VALS")->getDouble("SAFE_PARAM");

	if(firstCall && mpi_rank==0) {
		cout<<"> RELAXATION_CLIP_FOR_NEGATIVE_VALS (if some scalars become negative during NT, reduce relax): CLIP_PARAM = "
		    <<newtonParam.clipParameter<<"  SAFE_PARAM="<<newtonParam.safeParameter<<endl;
	}

	// Backtracking
	//		 *   backtrackMaxIter 			= "BACKTRACKING_PARAMETERS"->"MAX_ITER"
	//		 *   backtrackRelax_LowerBound 	= "BACKTRACKING_PARAMETERS"->"RELAX_LOWER_BOUND"
	//		 *   backtrackRelax_UpperBound 	= "BACKTRACKING_PARAMETERS"->"RELAX_UPPER_BOUND"
	//		 *   backtrackBarrierCoeff 		= "BACKTRACKING_PARAMETERS"->"BARRIER_COEFF"
	//		 *
	//		 *   skipBT_firstITer  = "BACKTRACKING_SKIP"->"FIRST_ITER"
	//		 *   skipBT_freq       = "BACKTRACKING_SKIP"->"FREQUENCY"
	if (!checkParam("BACKTRACKING_PARAMETERS")) {
		ParamMap::add("BACKTRACKING_PARAMETERS  MAX_ITER=10  RELAX_LOWER_BOUND=0.1  RELAX_UPPER_BOUND=0.5  BARRIER_COEFF=0.01"); // add default values
		if (mpi_rank == 0)
			cout << "WARNING: added keyword \"BACKTRACKING_PARAMETERS  MAX_ITER=10  RELAX_LOWER_BOUND=0.1  RELAX_UPPER_BOUND=0.5  BARRIER_COEFF=0.01\" to parameter map!"<<endl;
	}
	newtonParam.backtrackMaxIter = getParam("BACKTRACKING_PARAMETERS")->getInt("MAX_ITER");
	newtonParam.backtrackRelax_LowerBound = getParam("BACKTRACKING_PARAMETERS")->getDouble("RELAX_LOWER_BOUND");
	newtonParam.backtrackRelax_UpperBound = getParam("BACKTRACKING_PARAMETERS")->getDouble("RELAX_UPPER_BOUND");
	newtonParam.backtrackBarrierCoeff = getParam("BACKTRACKING_PARAMETERS")->getDouble("BARRIER_COEFF");
	if(firstCall && mpi_rank==0)
		printf("> BACKTRACKING_PARAMETERS  MAX_ITER=%d  RELAX_LOWER_BOUND=%.2e  RELAX_UPPER_BOUND=%.2e  BARRIER_COEFF=%.2e\n",
				newtonParam.backtrackMaxIter, newtonParam.backtrackRelax_LowerBound, newtonParam.backtrackRelax_UpperBound, newtonParam.backtrackBarrierCoeff);

	if (!checkParam("BACKTRACKING_SKIP")) {
		ParamMap::add("BACKTRACKING_SKIP  FIRST_ITER=FALSE  FREQUENCY=0"); // add default values
		if (mpi_rank == 0)
			cout << "WARNING: added keyword \"BACKTRACKING_SKIP  FIRST_ITER=FALSE  FREQUENCY=0\" to parameter map!"<<endl;
	}
	string skipBacktrack_firstIter = getParam("BACKTRACKING_SKIP")->getString("FIRST_ITER");
	if(strcmp(skipBacktrack_firstIter.c_str(), "TRUE")==0 || strcmp(skipBacktrack_firstIter.c_str(), "T")==0)
		newtonParam.skipBT_firstITer = true;
	else
		newtonParam.skipBT_firstITer = false;
	newtonParam.skipBT_freq = getParam("BACKTRACKING_SKIP")->getInt("FREQUENCY");
	if(firstCall && mpi_rank==0)
		printf("> BACKTRACKING_SKIP  FIRST_ITER=%d  FREQUENCY=%d\n", newtonParam.skipBT_firstITer, newtonParam.skipBT_freq);

	//  Stabilized Newton's method for the problem with a singular Jacobian
	//		 *   "STABILIZATION_FOR_SINGULAR_JAC"-->"METHOD_TYPE"
	//		 *   "STABILIZATION_FOR_SINGULAR_JAC"-->"ALPHA"
	//		 *   "STABILIZATION_FOR_SINGULAR_JAC"-->"ALPHA_EPS"
	if (!checkParam("STABILIZATION_FOR_SINGULAR_JAC")) {
		ParamMap::add("STABILIZATION_FOR_SINGULAR_JAC  METHOD_TYPE=NO_STAB_NEWTON  ALPHA=0.0  ALPHA_EPS=0.0"); // add default values
		if (mpi_rank == 0)
			cout << "WARNING: added keyword \"STABILIZATION_FOR_SINGULAR_JAC  METHOD_TYPE=NO_STAB_NEWTON  ALPHA=0.0  ALPHA_EPS=0.0\" to parameter map!"<<endl;
	}

	string methodTypeString = getParam("STABILIZATION_FOR_SINGULAR_JAC")->getString("METHOD_TYPE");
		// Note: Possible options are: NO_STAB_NEWTON, CONST_DIAG, CFL_BASED_DIAG, GENERAL_NEWTON_HU, GENERAL_NEWTON_HUSEO
	std::transform(methodTypeString.begin(), methodTypeString.end(), methodTypeString.begin(), ::toupper);
	if(methodTypeString.compare("CONST_DIAG") == 0)
		newtonParam.stabilizationMethodType = CONST_DIAG;
	else if(methodTypeString.compare("CFL_BASED_DIAG") == 0)
		newtonParam.stabilizationMethodType = CFL_BASED_DIAG;
	else if(methodTypeString.compare("GENERAL_NEWTON_HU") == 0)
		newtonParam.stabilizationMethodType = GENERAL_NEWTON_HU;
	else if(methodTypeString.compare("GENERAL_NEWTON_HUSEO") == 0)
		newtonParam.stabilizationMethodType = GENERAL_NEWTON_HUSEO;
	else
		newtonParam.stabilizationMethodType = NO_STAB_NEWTON;

	if(newtonParam.stabilizationMethodType != NO_STAB_NEWTON)
		newtonParam.stabilizationAlpha    = getParam("STABILIZATION_FOR_SINGULAR_JAC")->getDouble("ALPHA"); // Model coefficient
	if(newtonParam.stabilizationMethodType == GENERAL_NEWTON_HU || newtonParam.stabilizationMethodType == GENERAL_NEWTON_HUSEO)
		newtonParam.stabilizationAlphaEps = getParam("STABILIZATION_FOR_SINGULAR_JAC")->getDouble("ALPHA_EPS"); // Minimum coefficient

	if(firstCall && mpi_rank==0) {
		cout<<"> STABILIZATION_FOR_SINGULAR_JAC : METHOD_TYPE=";
		switch(newtonParam.stabilizationMethodType) {
			case NO_STAB_NEWTON:       cout<<"NO_STAB_NEWTON";       break;
			case CONST_DIAG:           cout<<"CONST_DIAG";           break;
			case CFL_BASED_DIAG:       cout<<"CFL_BASED_DIAG";       break;
			case GENERAL_NEWTON_HU:    cout<<"GENERAL_NEWTON_HU:";   break;
			case GENERAL_NEWTON_HUSEO: cout<<"GENERAL_NEWTON_HUSEO"; break;
			default: cerr<<"ERROR"; throw(PSALC_ERROR_CODE);
		}

		if(newtonParam.stabilizationMethodType != NO_STAB_NEWTON)
			cout<<"  ALPHA="<<newtonParam.stabilizationAlpha;
		if(newtonParam.stabilizationMethodType == GENERAL_NEWTON_HU || newtonParam.stabilizationMethodType == GENERAL_NEWTON_HUSEO)
			cout<<" (ALPHA_EPS = "<<newtonParam.stabilizationAlphaEps<<")";
		cout<<endl;
	}

	// Ramping alpha for STABILIZATION_FOR_SINGULAR_JAC
	if(newtonParam.stabilizationMethodType != NO_STAB_NEWTON) {
		// Check ramping
		if (!checkParam("RAMP_STAB_NEWTON_ALPHA")) {
			ParamMap::add("RAMP_STAB_NEWTON_ALPHA AFTER_ITER=10  INTERVAL_ITER=10  FACTOR_INC=1.0  TARGET_ALPHA=1.0");    // add default: no increase of CFL number!
			if (mpi_rank == 0)
				cout << "WARNING: added \"RAMP_STAB_NEWTON_ALPHA AFTER_ITER=10  INTERVAL_ITER=10  FACTOR_INC=1.0  TARGET_ALPHA=1.0\" to parameter map" << endl;
		}

		newtonParam.stabilizationAlphaRamp_startInc = getParam("RAMP_STAB_NEWTON_ALPHA")->getInt("AFTER_ITER");
		newtonParam.stabilizationAlphaRamp_interval = getParam("RAMP_STAB_NEWTON_ALPHA")->getInt("INTERVAL_ITER");
		newtonParam.stabilizationAlphaRamp_incAlpha    = getParam("RAMP_STAB_NEWTON_ALPHA")->getDouble("FACTOR_INC");
		newtonParam.stabilizationAlphaRamp_targetAlpha = getParam("RAMP_STAB_NEWTON_ALPHA")->getDouble("TARGET_ALPHA");

		if(firstCall && mpi_rank==0) {
			printf("  Ramping of stabilization ALPHA: AFTER_ITER=%d  INTERVAL_ITER=%d  FACTOR_INC=%g  TARGET_ALPHA=%g\n",
					newtonParam.stabilizationAlphaRamp_startInc, newtonParam.stabilizationAlphaRamp_interval,
					newtonParam.stabilizationAlphaRamp_incAlpha, newtonParam.stabilizationAlphaRamp_targetAlpha);
		}
	}

	// Trust-region size
	newtonParam.trustRegionSize = getDoubleParam("TRUST_REGION_SIZE", "1.0e6");

	if(firstCall && mpi_rank==0)
		cout<<"> TRUST_REGION_SIZE="<<newtonParam.trustRegionSize<<endl;

	firstCall = false;
}

/*
 * Method: getParamMoreNewtonSteps
 * -------------------------------
 * Obtain the parameters to run Newton's method for more steps after convergence
 *   "MORE_NTSTEPS_BELOW_CONV"-->"MORE_STEPS"
 *   "MORE_NTSTEPS_BELOW_CONV"-->"RESID"
 *   "MORE_NTSTEPS_BELOW_CONV"-->"RESID"
 */
void IkeWithPsALC_AD::getParamMoreNewtonSteps(int &moreNTstepsBelowConv_moreSteps, double &moreNTstepsBelowConv_resid, double &moreNTstepsBelowConv_deltaQ,
		const bool firstCall) {
	if (!checkParam("MORE_NTSTEPS_BELOW_CONV")) {
		ParamMap::add("MORE_NTSTEPS_BELOW_CONV  MORE_STEPS=0  RESID=1.0e-15  DELTA_Q=1.0e-15"); // add default values
		if (mpi_rank == 0)
			cout<< "WARNING: added keyword \"MORE_NTSTEPS_BELOW_CONV  MORE_STEPS=0  RESID=1.0e-15  DELTA_Q=1.0e-15\""<< " to parameter map!" << endl;
	}
	moreNTstepsBelowConv_moreSteps = getParam("MORE_NTSTEPS_BELOW_CONV")->getInt("MORE_STEPS");
	moreNTstepsBelowConv_resid     = getParam("MORE_NTSTEPS_BELOW_CONV")->getDouble("RESID");
	moreNTstepsBelowConv_deltaQ    = getParam("MORE_NTSTEPS_BELOW_CONV")->getDouble("DELTA_Q");

	if(firstCall && mpi_rank==0)
		cout<<"> MORE_NTSTEPS_BELOW_CONV  MORE_STEPS="<<moreNTstepsBelowConv_moreSteps<<"  RESID="<<moreNTstepsBelowConv_resid<<"  DELTA_Q = "<<moreNTstepsBelowConv_deltaQ<<endl;
}

/*
 * Method: getParamsModifiedNewtonMethod
 * -------------------------------------
 * Obtain the Newton method to be used.
 * Currently supported methods = BASIC, SHAMANSKII, MODIFIED_SHAMANSKII
 */
void IkeWithPsALC_AD::getParamsModifiedNewtonMethod(MODIFIED_NEWTON& modifiedNewtonMethod, double& modifNewtonRelResid, int& modifNewtonFreq) {
	modifNewtonRelResid = 1.0e10;
	modifNewtonFreq = 0;

	string tempoString = getStringParam("MODIFIED_NEWTON_NS", "BASIC");
	if(tempoString == "BASIC") {
		modifiedNewtonMethod = BASIC;
		if(mpi_rank==0)
			cout<<"> MODIFIED_NEWTON_NS = BASIC "<<endl;
	} else if (tempoString == "SHAMANSKII") {
		modifiedNewtonMethod = SHAMANSKII;
		if (!checkParam("MODIFIED_NEWTON_NS_PARAMETERS")) {
			ParamMap::add("MODIFIED_NEWTON_NS_PARAMETERS  REL_RESID=0.1  FREQ=3"); // add default values
			if (mpi_rank == 0)
				cout << "WARNING: added keyword \"MODIFIED_NEWTON_NS_PARAMETERS  REL_RESID=0.1  FREQ=3\" to parameter map!"<<endl;
		}
		modifNewtonRelResid = getParam("MODIFIED_NEWTON_NS_PARAMETERS")->getDouble("REL_RESID");
		modifNewtonFreq     = getParam("MODIFIED_NEWTON_NS_PARAMETERS")->getInt("FREQ");
		if(mpi_rank==0)
			printf("> MODIFIED_NEWTON_NS = SHAMANSKII, REL_RESID=%.3e  FREQ=%d\n", modifNewtonRelResid, modifNewtonFreq);
	} else if (tempoString == "MODIFIED_SHAMANSKII") {
		modifiedNewtonMethod = MODIFIED_SHAMANSKII;
		if (!checkParam("MODIFIED_NEWTON_NS_PARAMETERS")) {
			ParamMap::add("MODIFIED_NEWTON_NS_PARAMETERS  REL_RESID=0.1  FREQ=3"); // add default values
			if (mpi_rank == 0)
				cout << "WARNING: added keyword \"MODIFIED_NEWTON_NS_PARAMETERS  REL_RESID=0.1  FREQ=3\" to parameter map!"<<endl;
		}
		modifNewtonRelResid = getParam("MODIFIED_NEWTON_NS_PARAMETERS")->getDouble("REL_RESID");
		modifNewtonFreq     = getParam("MODIFIED_NEWTON_NS_PARAMETERS")->getInt("FREQ");
		if(mpi_rank==0)
			printf("> MODIFIED_NEWTON_NS = MODIFIED_SHAMANSKII, REL_RESID=%.3e  FREQ=%d\n", modifNewtonRelResid, modifNewtonFreq);
	} else {
		if(mpi_rank == 0) {
			cout<<"getModifiedNewtonMethod(): Unrecognizable MODIFIED_NEWTON_NS - "<<tempoString<<endl;
			cout<<"                           BASIC method will be used by default"<<endl;
		}
		modifiedNewtonMethod = BASIC;
	}
}


/*
 * Method: getParamsBEuler
 * -----------------------
 * Obtain some parameters for the backward-Euler time-integration method
 */
void IkeWithPsALC_AD::getParamsBEuler(int &maxIterLS_BEuler, double &zeroAbsLS_BEuler, double &zeroRelLS_BEuler,
		double &CflBEuler, int &startIncCFLBEuler, int &intervalIncCFLBEuler, double &incCFLBEuler, double &maxCFLBEuler) {
	if (!checkParam("LINEAR_SOLVER_BACKWARDEULER_TRESHOLDS")) {
		ParamMap::add("LINEAR_SOLVER_BACKWARDEULER_TRESHOLDS  MAX_ITER=10  ABS_RESID=1.0e-10  REL_RESID=1.0e-4");    // add default values
		if (mpi_rank == 0)
			cout << "WARNING: added keyword \"LINEAR_SOLVER_BACKWARDEULER_TRESHOLDS  MAX_ITER=10  ABS_RESID=1.0e-10  REL_RESID=1.0e-4\"" <<
			" to parameter map!" << endl;
	}
	maxIterLS_BEuler = getParam("LINEAR_SOLVER_BACKWARDEULER_TRESHOLDS")->getInt("MAX_ITER");
	zeroAbsLS_BEuler = getParam("LINEAR_SOLVER_BACKWARDEULER_TRESHOLDS")->getDouble("ABS_RESID");
	zeroRelLS_BEuler = getParam("LINEAR_SOLVER_BACKWARDEULER_TRESHOLDS")->getDouble("REL_RESID");

	CflBEuler = getDoubleParam("BACKWARDEULER_CFL", "10.0");

	if (!checkParam("BACKWARDEULER_CFL_RAMP")) {
		ParamMap::add("BACKWARDEULER_CFL_RAMP AFTER_ITER=10  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0");    // add default: no increase of CFL number!
		if (mpi_rank == 0)
			cout << "WARNING: added \"BACKWARDEULER_CFL_RAMP AFTER_ITER=10  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0\" to parameter map" << endl;
	}
	startIncCFLBEuler    = getParam("BACKWARDEULER_CFL_RAMP")->getInt("AFTER_ITER");
	intervalIncCFLBEuler = getParam("BACKWARDEULER_CFL_RAMP")->getInt("INTERVAL_ITER");
	incCFLBEuler         = getParam("BACKWARDEULER_CFL_RAMP")->getDouble("FACTOR_CFL");
	maxCFLBEuler         = getParam("BACKWARDEULER_CFL_RAMP")->getDouble("MAX_CFL");
}

/*
 * Method: getHowToCalcJac
 * -----------------------
 * Obtain the method to calculate the Jacobian matrix.
 * Currently supported methods = ROW_1D, ORDINARY_2D
 */
HOW_TO_CALC_JAC IkeWithPsALC_AD::getHowToCalcJac() {
	static bool firstCall = true;

	if(firstCall) {
		if (!checkParam("HOW_TO_CALC_JAC")) {
			ParamMap::add("HOW_TO_CALC_JAC=ORDINARY_2D"); // add default values
			if (mpi_rank == 0)
				cout << "WARNING: added keyword \"HOW_TO_CALC_JAC=ORDINARY_2D\" to parameter map!"<<endl;
		}
	}

	HOW_TO_CALC_JAC howToCalcJac;
	string tempString = getStringParam("HOW_TO_CALC_JAC", "ORDINARY_2D");
	if(tempString=="ORDINARY_2D") {
		howToCalcJac = ORDINARY_2D;
		if(firstCall && mpi_rank==0) {
			cout<<"> HOW_TO_CALC_JAC = ORDINARY_2D"<<endl;
			if(int(cvora[mpi_rank+1]/mpi_size)>2000 && mpi_rank==0) {
				cout<<"    Averaged ncv per core is "<<int(cvora[mpi_rank+1]/mpi_size)<<":"<<endl;
				cout<<"        it is highly recommended to use HOW_TO_CALC_JAC = ROW_1D"<<endl;
			}
		}
	} else if (tempString=="ROW_1D") {
		howToCalcJac = ROW_1D;
		if(firstCall && mpi_rank==0)
			cout<<"> HOW_TO_CALC_JAC = ROW_1D  (i.e. memory saving mode)"<<endl;
	} else {
		if(mpi_rank==0)
			cout<<"> unsupported type: HOW_TO_CALC_JAC = "<<tempString<<endl;
		throw(PSALC_ERROR_CODE);
	}

	firstCall = false;

	return howToCalcJac;
}

/*
 * Method: vecDotVecWithoutWeight
 * ------------------------------
 * Calculate the inner product between [q0, lambda0] and [q1, lambda1] without weight,
 * i.e., calculate [q0, lambda0]^{T} * [q1, lambda1]
 *
 * Note: 1. The dimension of q must be ncv*nVars and the dimension of lambda must be NcontrolEqns.
 *       2. This is a MPI call, i.e., all the cores must have the q vectors,
 *         and only the last core (mpi_rank==mpi_size-1) must have the lambda vectors.
 *
 * Return:
 *   By value     = the inner product
 *   By reference = innerProductRatio, i.e., (q0^{T} * q1) / (lambda0^{T} * lambda1)
 */
double IkeWithPsALC_AD::vecDotVecWithoutWeight (double* innerProductRatio,
		const double* qVec0, const double* qVec1, const double* lambdaVec0, const double* lambdaVec1,
		const int Nvars, const int NcontrolEqns) {
	assert(Nvars >= 5);
	if(mpi_rank == mpi_size-1)
		assert(lambdaVec0 != NULL && lambdaVec1 != NULL);

	/* Calculate the inner product between the two q vectors at each core */
	double myQvecCalc = 0.0;
	for(int icv=0; icv<ncv; ++icv) {
		int indexStart = icv*Nvars;

		for(int ivar=0; ivar<Nvars; ++ivar) {
			int index = indexStart + ivar;
			myQvecCalc += qVec0[index] * qVec1[index];
		}
	}

	/* Calculate the inner product between the two lambda vectors (only mpi_rank==mpi_size-1 has it!) */
	double lambdaCalc = 0.0;
	if(mpi_rank == mpi_size -1) {
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
			lambdaCalc += lambdaVec0[iEqn] * lambdaVec1[iEqn];
		}
	}

	/* Sum up the inner product of the q vectors from all the cpu cores, and broadcast the inner product of the lambda vectors */
	double totQvecCalc;
	MPI_Allreduce(&myQvecCalc, &totQvecCalc, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

	MPI_Bcast(&lambdaCalc, 1, MPI_DOUBLE, mpi_size-1, mpi_comm);

	/* Return */
	if(innerProductRatio != NULL) {
		if(fabs(lambdaCalc) < MACHINE_EPS) {
			if(mpi_rank==0)
				cout<<"WARNING in IkeWithPsALC_AD::vecDotVecWithoutWeight(): You asked to calculate innerProductRatio,"<<endl
					<<"                                                      i.e., (q0^{T}*q1)/(lambda0^{T}*lambda1),"<<endl
					<<"                                                      but the demoninator (inner product between the lambda vectors)"<<endl
					<<"                                                      is almost zero -- value = "<<lambdaCalc<<endl;
		}
		*innerProductRatio = totQvecCalc / lambdaCalc;
	}

	return totQvecCalc + lambdaCalc;
}

double IkeWithPsALC_AD::vecDotVecWithoutWeight (const double* qVec0, const double* qVec1, const double* lambdaVec0, const double* lambdaVec1,
		const int Nvars, const int NcontrolEqns) {
	return vecDotVecWithoutWeight (NULL, qVec0, qVec1, lambdaVec0, lambdaVec1, Nvars, NcontrolEqns);
}

/*
 * Method: getFlowWeightsForInnerProduct
 * -------------------------------------
 * Get weights for the inner product between flow vectors
 */
void IkeWithPsALC_AD::getFlowWeightsForInnerProduct(double* weightVec, const double defaultConst, const int icv, const int nVars) {
	assert(weightVec!=NULL && nVars>0);

	int nScal = nVars-5;

	// Initialize
	for(int i=0; i<nVars; ++i) {
#ifdef USE_VOLUME_WEIGHTED_INNER_PROD
//		weightVec[i] = defaultConst/cv_volume[icv];
		weightVec[i] = defaultConst/sqrt(max(cv_volume[icv],MIN_CV_VOL_CUTOFF));
#else
		weightVec[i] = defaultConst;
#endif
	}

	// Get weights for each variable
#ifdef USE_REF_WEIGHTED_INNER_PROD
	weightVec[0]       *= 1.0/fabs(RefFlowParams.rho_ref);
	for(int i=1; i<=3; ++i)
		weightVec[1+i] *= 1.0/RefFlowParams.rhouMag_ref;
	weightVec[4]       *= 1.0/fabs(RefFlowParams.rhoE_ref);
	if(nScal>0) {
		for(size_t iScal=0; iScal<RefFlowParams.nScal; ++iScal)
			weightVec[5+iScal] *= 1.0 / fabs(RefFlowParams.scalar_ref[iScal]);
	}
#endif
}

/*
 * Method: vecDotVecWithWeight
 * ---------------------------
 * Calculate the inner product between [q0, lambda0] and [q1, lambda1] with weights for lambda,
 * i.e., calculate [q0, lambda0]^T W [q1, lambda1],  where W = [ I 0
 *                                                               0 w ]
 * Note: 1. The dimension of q must be ncv*nVars and the dimension of lambda must be NcontrolEqns.
 *       2. This is a MPI call, i.e., all the cores must have the q vectors,
 *         and only the last core (mpi_rank==mpi_size-1) must have the lambda vectors.
 *
 *       3. If you modify the weights here, you must commit the same change for q_tangent and lambda_tangent in getSteadySolnByNewton()
 *
 * Return:
 *   By value     = the weighted inner product
 *   By reference = innerProductRatio, i.e., (q0^{T} * q1) / (lambda0^{T} * lambda1)  (note that (lambda0^{T} * lambda1) is NOT a weighted norm)
 */
double IkeWithPsALC_AD::vecDotVecWithWeight (double* innerProductRatio,
		const double* qVec0, const double* qVec1, const double* lambdaVec0, const double* lambdaVec1, const double weightForLambda,
		const int Nvars, const int NcontrolEqns) {
	assert(weightForLambda > 0.0);
	assert(Nvars >= 5);
	if(NcontrolEqns > 0 && mpi_rank == mpi_size-1)
		assert(lambdaVec0 != NULL && lambdaVec1 != NULL);

	static bool firstCall = true;

	/* Allocate memory */
	double* weightVec = new double [Nvars];

	/* Calculate the inner product between the two q vectors at each core */
	double myQvecCalc = 0.0;
	for(int icv=0; icv<ncv; ++icv) {
		int indexStart = icv*Nvars;

#ifdef USE_VOLUME_WEIGHTED_INNER_PROD
		getFlowWeightsForInnerProduct(weightVec, sqrt(max(totMinCvVol,MIN_CV_VOL_CUTOFF)), icv, Nvars);
#else
		getFlowWeightsForInnerProduct(weightVec, 1.0,                                      icv, Nvars);
#endif

		for(int ivar=0; ivar<Nvars; ++ivar) {
			int index = indexStart + ivar;
			myQvecCalc += weightVec[ivar] * qVec0[index] * qVec1[index];
		}
	}

	/* Calculate the inner product (not-weighted by weightLambda) between the two lambda vectors (only mpi_rank==mpi_size-1 has it!) */
	double lambdaCalc = 0.0;
	if(mpi_rank == mpi_size -1) {
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
			lambdaCalc += lambdaVec0[iEqn] * lambdaVec1[iEqn];
		}
	}

	/* Sum up the inner product of the q vectors from all the cpu cores, and broadcast the inner product of the lambda vectors */
	double totQvecCalc;
	MPI_Allreduce(&myQvecCalc, &totQvecCalc, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

	MPI_Bcast(&lambdaCalc, 1, MPI_DOUBLE, mpi_size-1, mpi_comm);

	/* Clear memory if necessary */
	delete [] weightVec;

	/* Return */
	if(innerProductRatio != NULL) {
		if(fabs(lambdaCalc) < MACHINE_EPS) {
			if(mpi_rank==0)
				cout<<"WARNING in IkeWithPsALC_AD::vecDotVecWithWeight(): You asked to calculate innerProductRatio,"<<endl
					<<"                                                   i.e., (q0^{T}*W_q*q1)/(lambda0^{T}*W_lambda*lambda1),"<<endl
					<<"                                                   but the demoninator (inner product between the lambda vectors)"<<endl
					<<"                                                   is almost zero -- value = "<<lambdaCalc<<endl;
		}
		*innerProductRatio = totQvecCalc / lambdaCalc;
	}

	firstCall = false;

	return totQvecCalc + weightForLambda*lambdaCalc;
}

double IkeWithPsALC_AD::vecDotVecWithWeight (const double* qVec0, const double* qVec1, const double* lambdaVec0, const double* lambdaVec1, const double weightForLambda,
		const int Nvars, const int NcontrolEqns) {
	return vecDotVecWithWeight (NULL, qVec0, qVec1, lambdaVec0, lambdaVec1, weightForLambda, Nvars, NcontrolEqns);
}

/*
 * Method: calcUnweighted2NormForVecMinusVec
 * -----------------------------------------
 * Calculate unweighted 2-norm for a vector which is given in the form of [q1-q0; lambda1-lambda0],
 * i.e., calculate sqrt( [q1-q0; lambda1-lambda0]^T * [q1-q0; lambda1-lambda0] )
 *
 * Note: 1. The dimension of q must be ncv*nVars and the dimension of lambda must be NcontrolEqns.
 *       2. This is a MPI call, i.e., all the cores must have the q vectors,
 *         and only the last core (mpi_rank==mpi_size-1) must have the lambda vectors.
 *
 * Return:
 *   By value     = the unweighted norm
 *   By reference = normSqRatio, i.e., ||q1-q0||^2 / ||lambda1-lambda0||^2
 */
double IkeWithPsALC_AD::calcUnweighted2NormForVecMinusVec(double *normSqRatio,
		const double* qVec1, const double* qVec0, const double* lambdaVec1, const double* lambdaVec0,
		const int Nvars, const int NcontrolEqns) {
	assert(Nvars >= 5);

	// Allocate memory
	double* qVecDiff      = new double [ncv*Nvars];
	double* lambdaVecDiff = NULL;
	if(mpi_rank == mpi_size-1)
		lambdaVecDiff = new double [NcontrolEqns];

	// Calculate the squre of the weighted 2-norm
	for(int i=0; i<ncv*Nvars; ++i)
		qVecDiff[i] = qVec1[i] - qVec0[i];
	if(mpi_rank == mpi_size-1)
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
			lambdaVecDiff[iEqn] = lambdaVec1[iEqn] - lambdaVec0[iEqn];

	double normSq = vecDotVecWithoutWeight(normSqRatio, qVecDiff, qVecDiff, lambdaVecDiff, lambdaVecDiff, Nvars, NcontrolEqns);

	// Free memory
	delete [] qVecDiff;
	if(mpi_rank == mpi_size-1)
		delete [] lambdaVecDiff;

	// Return
	return sqrt(normSq);
}

double IkeWithPsALC_AD::calcUnweighted2NormForVecMinusVec(const double* qVec1, const double* qVec0, const double* lambdaVec1, const double* lambdaVec0,
		const int Nvars, const int NcontrolEqns) {
	return calcUnweighted2NormForVecMinusVec(NULL, qVec1, qVec0, lambdaVec1, lambdaVec0, Nvars, NcontrolEqns);
}

/*
 * Method: calcWeighted2NormForVecMinusVec
 * ---------------------------------------
 * Calculate weighted 2-norm for a vector which is given in the form of [q1-q0; lambda1-lambda0],
 * i.e., calculate sqrt( [q1-q0; lambda1-lambda0]^T W [q1-q0; lambda1-lambda0] ),  where W = [ I 0
 *                                                                                             0 w ]
 * Note: 1. The dimension of q must be ncv*nVars and the dimension of lambda must be NcontrolEqns.
 *       2. This is a MPI call, i.e., all the cores must have the q vectors,
 *         and only the last core (mpi_rank==mpi_size-1) must have the lambda vectors.
 *
 * Return:
 *   By value     = the weighted norm
 *   By reference = normSqRatio, i.e., ||q1-q0||^2 / ||lambda1-lambda0||^2  (note that ||lamda1-lambda0|| is NOT a weighted norm)
 */
double IkeWithPsALC_AD::calcWeighted2NormForVecMinusVec(double *normSqRatio,
		const double* qVec1, const double* qVec0, const double* lambdaVec1, const double* lambdaVec0, const double weightForLambda,
		const int Nvars, const int NcontrolEqns) {
	assert(weightForLambda > 0.0);
	assert(Nvars >= 5);

	// Allocate memory
	double* qVecDiff      = new double [ncv*Nvars];
	double* lambdaVecDiff = NULL;
	if(mpi_rank == mpi_size-1 && NcontrolEqns > 0)
		lambdaVecDiff = new double [NcontrolEqns];

	// Calculate the squre of the weighted 2-norm
	for(int i=0; i<ncv*Nvars; ++i)
		qVecDiff[i] = qVec1[i] - qVec0[i];
	if(mpi_rank == mpi_size-1)
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
			lambdaVecDiff[iEqn] = lambdaVec1[iEqn] - lambdaVec0[iEqn];

	double normSq = vecDotVecWithWeight(normSqRatio, qVecDiff, qVecDiff, lambdaVecDiff, lambdaVecDiff, weightForLambda, Nvars, NcontrolEqns);

	// Free memory
	delete [] qVecDiff;
	if(mpi_rank == mpi_size-1 && NcontrolEqns > 0)
		delete [] lambdaVecDiff;

	// Return
	return sqrt(normSq);
}

double IkeWithPsALC_AD::calcWeighted2NormForVecMinusVec(const double* qVec1, const double* qVec0, const double* lambdaVec1, const double* lambdaVec0, const double weightForLambda,
		const int Nvars, const int NcontrolEqns) {
	return calcWeighted2NormForVecMinusVec(NULL, qVec1, qVec0, lambdaVec1, lambdaVec0, weightForLambda, Nvars, NcontrolEqns);
}

/*
 * Method: calcUnweightedTangentVecDotFlowVecMinusGuess
 * ----------------------------------------------------
 * Calculate unweighted vector multiplication of [q_tangent; lambda_tangent]^T * [flowVec - q_guess; lambda - lambda_guess],
 * where flowVec = [rho; rhou; rhoE; scalars].
 *
 * This will be called mainly if time stepping is used instead of Newton's method
 *
 * Note: 1. The dimension of qTangent and qGuess must be ncv*nVars and the dimension of lambda and lambdaGuess can be either zero or NcontrolEqns.
 *       2. This is a MPI call, i.e., all the cores must have the q vectors,
 *          and only the last core (mpi_rank==mpi_size-1) can have the lambda vectors.
 *       3. If NcontrolEqns == 0, innerProductRatio will be set as -1.
 *
 * Return:
 *   By value     = the weighted norm
 *   By reference = normSqRatio
 */
double IkeWithPsALC_AD::calcUnweightedTangentVecDotFlowVecMinusGuess(const int iEqn, double** q_tangent, const double* q_guess,
		const double* lambda_tangent, const double* lambda, const double* lambda_guess,
		const int Nvars, const int NcontrolEqns) {
	assert(q_tangent != NULL && q_guess != NULL);
	if(NcontrolEqns != 0 && mpi_rank == mpi_size-1)
		assert(lambda_tangent != NULL && lambda != NULL && lambda_guess != NULL);
	assert(Nvars >= 5);

	// Allocate memory for qVecDiff and lambdaVecDiff
	double* qVecDiff      = new double [ncv*Nvars];
	double* lambdaVecDiff = NULL;
	if(mpi_rank == mpi_size-1 && NcontrolEqns > 0)
		lambdaVecDiff = new double [NcontrolEqns];

	// Calculate qVecDiff and lambdaVecDiff
	for(int icv=0; icv<ncv; ++icv) {
		int index = icv*Nvars;

		// Primary variables
		qVecDiff[index] = rho[icv] - q_guess[index];
		for (int i=0; i<3; ++i)
			qVecDiff[index+1+i] = rhou[icv][i] - q_guess[index+1+i];
		qVecDiff[index+4] = rhoE[icv] - q_guess[index+4];
		// Scalars
		for (int iScal=0; iScal < Nvars-5; ++iScal) {
			double *phi = scalarTranspEqVector[iScal].phi;
			qVecDiff[index+5+iScal] = phi[icv] - q_guess[index+5+iScal];
		}
	}
	if(mpi_rank == mpi_size-1 && NcontrolEqns > 0)
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
			lambdaVecDiff[iEqn] = lambda[iEqn] - lambda_guess[iEqn];

	// Calculate the square of the unweighted 2-norm
	double vecDotVec;
	if(NcontrolEqns == 0) {
		double myQvecCalc = 0.0;
		for(int i=0; i<ncv*Nvars; ++i)
			myQvecCalc += q_tangent[iEqn][i] * qVecDiff[i];
		MPI_Allreduce(&myQvecCalc, &vecDotVec, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
	} else {
		double myQvecCalc = 0.0;
		for(int i=0; i<ncv*Nvars; ++i)
			myQvecCalc += q_tangent[iEqn][i] * qVecDiff[i];
		if(mpi_rank == mpi_size -1) {
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
				myQvecCalc += lambda_tangent[iEqn] * lambdaVecDiff[iEqn];
			}
		}
		MPI_Allreduce(&myQvecCalc, &vecDotVec, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
	}

	// Free memory
	delete [] qVecDiff;
	if(mpi_rank == mpi_size-1 && NcontrolEqns > 0)
		delete [] lambdaVecDiff;

	// Return
	return vecDotVec;
}

/*
 * Method: calcUnweightedTangentVecDotJOEdq
 * ----------------------------------------------------
 * Calculate unweighted vector multiplication of [q_tangent; lambda_tangent]^T [dq; dLambda],
 * where dq is JOE vector
 */
double IkeWithPsALC_AD::calcUnweightedTangentVecDotJOEdq(const int iEqn, double** q_tangent, double* lambda_tangent, double (*dq)[5], double** dScal, const int Nvars) {
	double myVecDotVec = 0.0;

	for(int icv=0; icv<ncv; ++icv) {
		int indexStart = icv*Nvars;
		for(int i=0; i<5; ++i)
			myVecDotVec += q_tangent[iEqn][indexStart+i] * dq[icv][i];
		for(int iScal=0; iScal<nScal; ++iScal)
			myVecDotVec += q_tangent[iEqn][indexStart+5+iScal] * dScal[iScal][icv];
	}

	double vecDotVec;
	MPI_Allreduce(&myVecDotVec, &vecDotVec, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

	return vecDotVec;
}

/*
 * Method: calcUnitVecWithWeightedNorm
 * -----------------------------------
 * Calculate the normalized form of a given vector [q; lambda],
 * i.e., calculate [u_q; u_lambda] = 1/(|| [q; lambda] ||_W) * [q; lambda],  where || [q; lambda] ||_W is the weighted 2-norm
 *
 * Note: 1. The dimension of q must be ncv*nVars and the dimension of lambda must be NcontrolEqns.
 *       2. This is a MPI call, i.e., all the cores must have the q vector,
 *         and only the last core (mpi_rank==mpi_size-1) must have the lambda vector.
 *       3. This method does not calculate the weighted 2-norm (i.e., || [q; lambda] ||_W).
 *          It must be given by the user.
 * Return:
 *   By reference = u_q and u_lambda
 */
void IkeWithPsALC_AD::calcUnitVecWithWeightedNorm(double* u_q, double* u_lambda, const double* qVec, const double* lambdaVec, const double weighted2norm,
		const int Nvars, const int NcontrolEqns) {
	assert(Nvars >= 5);
	assert(u_q != NULL && qVec != NULL);
	if(mpi_rank == mpi_size-1)
		assert(u_lambda != NULL && lambdaVec != NULL);

	for(int i=0; i<ncv*Nvars; ++i)
		u_q[i] = qVec[i] / weighted2norm;

	if(mpi_rank==mpi_size-1)
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
			u_lambda[iEqn] = lambdaVec[iEqn] / weighted2norm;
}

/*
 * Method: showOverviewOnTangentVectors
 * ------------------------------------
 * Here we give the overview on the tangential vectors
 */
void IkeWithPsALC_AD::showOverviewOnTangentVectors(const double* u_q, const double* u_lambda, const int Nvars, const int NcontrolEqns) {
	// BINning
	int TOT_BIN_SIZE = 16;
	int my_bin_count[TOT_BIN_SIZE];
	for (int ibin=0; ibin<TOT_BIN_SIZE; ++ibin)
		my_bin_count[ibin] = 0;

	for(int i=0; i<ncv*Nvars; ++i) {
		double val = fabs(u_q[i]);
		for (int ibin=0; ibin<TOT_BIN_SIZE-1; ++ibin) {
			if(val<=pow(10.0, -double(ibin)) && val>pow(10.0, -double(ibin+1))) {
				++my_bin_count[ibin];
				break;
			}
		}
		if(val<=pow(10.0, -double(TOT_BIN_SIZE-1)))
			++my_bin_count[TOT_BIN_SIZE-1];
	}

	int bin_count[TOT_BIN_SIZE];
	MPI_Allreduce(my_bin_count, bin_count, TOT_BIN_SIZE, MPI_INT, MPI_SUM, mpi_comm);

	// Min and Max for each variables
	double my_min_uq[Nvars], my_max_uq[Nvars];
	for(int i=0; i<Nvars; ++i) {
		my_min_uq[i] = ABSURDLY_BIG_NUMBER;
		my_max_uq[i] = 0;
	}

	for(int icv=0; icv<ncv; ++icv) {
		int index = icv*Nvars;
		for(int i=0; i<Nvars; ++i) {
			double val = fabs(u_q[index]);
			my_min_uq[i] = min(my_min_uq[i], val);
			my_max_uq[i] = max(my_max_uq[i], val);
			++index;
		}
	}

	double min_uq[Nvars], max_uq[Nvars];
	MPI_Allreduce(my_min_uq, min_uq, Nvars, MPI_DOUBLE, MPI_MIN, mpi_comm);
	MPI_Allreduce(my_max_uq, max_uq, Nvars, MPI_DOUBLE, MPI_MAX, mpi_comm);

	// Show the information on the screen
	if(mpi_rank==mpi_size-1) {
		cout<<endl
				<<"IkeWithPsALC_AD::showOverviewOnTangentVectors(): Information on q_tangent >>"<<endl;
		for (int ibin=0; ibin<TOT_BIN_SIZE-1; ++ibin)
			printf("     BIN%d (10^%d ~ 10^%d): count = %d\n", ibin, -(ibin+1), -ibin, bin_count[ibin]);
		printf("     BIN%d (      ~ 10^%d): count = %d\n", TOT_BIN_SIZE-1, -(TOT_BIN_SIZE-1), bin_count[TOT_BIN_SIZE-1]);
		cout<<"     tot # = "<<cvora[mpi_size]*Nvars<<endl;
		cout<<"     ------------------"<<endl;
		for(int i=0; i<Nvars; ++i)
			printf("     VAR%d: %e ~ %e\n", i, min_uq[i], max_uq[i]);

		if(NcontrolEqns > 0)
			cout<<"---------------------"<<endl
			<<"IkeWithPsALC_AD::showOverviewOnTangentVectors(): Information on lambda_tangent >>"<<endl
			<<"     |lambda_tangent[0]| = "<<u_lambda[0]<<endl
			<<endl;
	}

	MPI_Barrier(mpi_comm);
}

/*
 * Method: calcUnitVec
 * -------------------
 * Calculate the normalized form of a given vector [q; lambda],
 * i.e., calculate [u_q; u_lambda] = 1/(|| [q; lambda] ||_W) * [q; lambda],  where || [q; lambda] ||_W is the weighted 2-norm
 *
 * Note: 1. The dimension of q must be ncv*nVars and the dimension of lambda must be NcontrolEqns.
 *       2. This is a MPI call, i.e., all the cores must have the q vector,
 *         and only the last core (mpi_rank==mpi_size-1) must have the lambda vector.
 * Return:
 *   By value     = the weighted 2-norm
 *   By reference = normSqRatio, i.e., ||q1-q0||^2 / ||lambda1-lambda0||^2  (note that ||lamda1-lambda0|| is NOT a weighted norm)
 *                  u_q and u_lambda
 */
double IkeWithPsALC_AD::calcUnitVec(double *normSqRatio,
		double* u_q, double* u_lambda, const double* qVec, const double* lambdaVec, const double weightForLambda,
		const int Nvars, const int NcontrolEqns) {
	assert(weightForLambda > 0.0);
	assert(Nvars >= 5);
	assert(u_q != NULL && qVec != NULL);
	if(mpi_rank == mpi_size-1)
		assert(u_lambda != NULL && lambdaVec != NULL);

	// Calculate the weighted 2-norm
	double weighted2normSq = vecDotVecWithWeight(normSqRatio, qVec, qVec, lambdaVec, lambdaVec, weightForLambda, Nvars, NcontrolEqns);
	double weighted2norm   = sqrt(weighted2normSq);

	// Get the normalized vectors and return
	calcUnitVecWithWeightedNorm(u_q, u_lambda, qVec, lambdaVec, weighted2norm, Nvars, NcontrolEqns);

	return weighted2norm;
}

double IkeWithPsALC_AD::calcUnitVec(double* u_q, double* u_lambda, const double* qVec, const double* lambdaVec, const double weightForLambda,
		const int Nvars, const int NcontrolEqns) {
	return calcUnitVec(NULL, u_q, u_lambda, qVec, lambdaVec, weightForLambda, Nvars, NcontrolEqns);
}

/*
 * Method: calcResidualsFrom1Drhs
 * --------------------------------
 * Calculate residual
 * Supported norm: inf-norm, one-norm, two-norm
 *                   whichNorm==0 -> infinity-norm
 *                   whichNorm==1 -> one-norm
 *                   whichNorm==2 -> two-norm
 * Note: "Residual" should have been defined as " double *Residual = new double[5+nScal]; "
 */
void IkeWithPsALC_AD::calcResidualsFrom1Drhs(double *Residual, double *rhs1Darray, const int whichNorm) {
	int nVars = 5+nScal;
	double *myResidual = new double[nVars];

	switch (whichNorm) {
	case 0: // note: infinity-norm
		for (int i = 0; i < 5+nScal; i++) {
			myResidual[i] = -ABSURDLY_BIG_NUMBER;
			Residual[i]   = -ABSURDLY_BIG_NUMBER;
		}

		for (int icv = 0; icv < ncv; icv++) {
			int icv_i = nVars*icv;

			for(int i=0; i<5+nScal; ++i)
				if(fabs(rhs1Darray[icv_i+i]) > myResidual[i])
					myResidual[i] = fabs(rhs1Darray[icv_i+i]);
		}
		MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_MAX, mpi_comm);

		break;
	case 1:
		for (int i = 0; i < 5+nScal; i++) {
			myResidual[i] = 0.0;
			Residual[i]   = 0.0;
		}

		for (int icv = 0; icv < ncv; icv++) {
			int icv_i = nVars*icv;

			for(int i=0; i<5+nScal; ++i)
				myResidual[i] += fabs(rhs1Darray[icv_i+i]);
		}
		MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);

		break;
	case 2:
		for (int i = 0; i < 5+nScal; i++) {
			myResidual[i] = 0.0;
			Residual[i]   = 0.0;
		}

		for (int icv = 0; icv < ncv; icv++) {
			int icv_i = nVars*icv;

			for(int i=0; i<5+nScal; ++i)
				myResidual[i] += fabs(rhs1Darray[icv_i+i])*fabs(rhs1Darray[icv_i+i]);
		}
		MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);
		for (int i = 0; i < 5+nScal; i++) {
			Residual[i] = sqrt(Residual[i]);
		}

		break;
	default:
		if(mpi_rank==0)
			cout<<"ERROR in IkeWithPsALC_AD::calcResidualsFrom1Drhs(): unsupported norm = "<<whichNorm<<endl;
		throw(PSALC_ERROR_CODE);
		break;
	}

	delete [] myResidual;
}

/*
 * Method: calcSumResidual
 * -----------------------
 * Calculate total residual from the residual vector
 * Supported norm: inf-norm, one-norm, two-norm
 *                   whichNorm==0 -> infinity-norm
 *                   whichNorm==1 -> one-norm
 *                   whichNorm==2 -> two-norm
 * Note: "Residual" should have been defined as " double *Residual = new double[5+nScal]; "
 */
double IkeWithPsALC_AD::calcSumResidual(const double *Residual, const int whichNorm) {
	assert(Residual != NULL);

	int nVars = 5+nScal;

	double residTot;
	switch(whichNorm) {
	case 0:
		residTot = -ABSURDLY_BIG_NUMBER;
		for(int i=0; i<nVars; ++i) {
			if(residTot < Residual[i])
				residTot = Residual[i];
		}
		break;
	case 1:
		residTot = 0.0;
		for(int i=0; i<nVars; ++i)
			residTot += Residual[i];
		break;
	case 2:
		residTot = 0.0;
		for(int i=0; i<nVars; ++i)
			residTot += Residual[i]*Residual[i];
		residTot = sqrt(residTot);
		break;
	default:
		if(mpi_rank==0)
			cout<<"ERROR in IkeWithPsALC_AD::calcSumResidual(): unsupported norm = "<<whichNorm<<endl;
		throw(PSALC_ERROR_CODE);
		break;
	}

	return residTot;
}

/*
 * Method: calcNres
 * ----------------
 * Calculate Nres == ds - { q_tangent*(q_guess-q1) + lambda_tangent*(lambda_guess-lambda1) }
 * (Nres is the last element of the rhs vector at mpi_rank=mpi_size-1)
 *
 * Return: True if the tangential condition is satisfied
 *
 * Note: Even though this part of the code is designed to handle multiple parameters (or lambdas), the actual formulation here is only for one parameters
 */
bool IkeWithPsALC_AD::calcNres(double* Nres, const double *q, double* rhs, const int Nvars, const int NcontrolEqns, const double *q1, double** q_tangent,
		const double* lambda_tangent, const double *lambda1, const double weightLambda, const double arcLength) {
	assert(NcontrolEqns>0);
	assert(Nres != NULL);
	assert(q_tangent!=NULL && lambda_tangent!=NULL);

	// Allocate memory
	double *lambdaNew = NULL;
	if(mpi_rank == mpi_size -1)
		lambdaNew = new double [NcontrolEqns];

	// Calculate the new arclength based on updated q
	if(mpi_rank == mpi_size-1) {
		for(int iParam=0; iParam<NcontrolEqns; ++iParam)
			lambdaNew[iParam] = q[ncv*Nvars+iParam];
	}

	double arcLengthNew = calcArclength(q1, q, lambda1, lambdaNew, NcontrolEqns, nScal, weightLambda, q_tangent, lambda_tangent);

	// Error check
	if(isNaN(arcLengthNew)) {
		if(mpi_rank==0) cerr<<"Error in IkeWithPsALC_AD::calcNres(): The new arclength calculated for the tangential eqns is NaN"<<endl;
		delete [] lambdaNew;
		throw(PSALC_ERROR_CODE);
	}

	// Calculate Nres = (old arclength) - (new arclength)
	double NresVal = arcLength - arcLengthNew;

//IKJ
NresVal *= AreaLambda;

	for(int iParam=0; iParam<NcontrolEqns; ++iParam)
		Nres[iParam] = NresVal;

	// Free memory
	if(lambdaNew != NULL)
		delete [] lambdaNew;

	// Return true if the tangential condition is satisfied, i.e., check if (new arclength) - (old arclength) == 0
	if(NresVal > 1.0e3*double(NcontrolEqns*cvora[mpi_size])*MACHINE_EPS) {
		return false;
	}

	return true;
}

/*
 * Method: calcNresForJOE
 * ----------------------
 * Calculate Nres == ds - {q_tangent*(q_guess-q1) + lambda_tangent*(lambda_guess-lambda1)}
 * (Nres is the last element of the rhs vector at mpi_rank=mpi_size-1)
 *
 * The calcNres() method uses the q arrays (q_guess and q1). Thus, it cannot be used if q is not updated (e.g. in the backtracking algorithm).
 * However, this method uses flow variables (rho, rhou, rhoE, ...) so that it can be called during the backtracking algorithm.
 * (Instead of updating q, you must update the flow variables before calling this method)
 *
 * Return: True if the tangential condition is satisfied
 *
 * Note: Even though this part of the code is designed to handle multiple parameters (or lambdas), the actual formulation here is only for one parameters
 */
bool IkeWithPsALC_AD::calcNresForJOE(double* Nres, const int Nvars, const int NcontrolEqns, const double *q1, double** q_tangent,
		const double* lambda_tangent, const double *lambdaNew, const double *lambda1, const double weightLambda, const double arcLength) {
	assert(NcontrolEqns>0);
	assert(Nres != NULL);
	assert(q_tangent!=NULL && lambda_tangent!=NULL);

	// Allocate memory
	double *qVec = new double [ncv*Nvars];

	// Calculate the new arclength based on updated flow variables
	for (int icv=0; icv<ncv; ++icv) {
		int dummy = icv*(5+nScal);
		qVec[dummy+0] = rho[icv];
		for(int i=0; i<3; ++i)
			qVec[dummy+1+i] = rhou[icv][i];
		qVec[dummy+4] = rhoE[icv];
		for(int iScal=0; iScal<nScal; ++iScal) {
			qVec[dummy+5+iScal] = UgpWithCvCompFlow::scalarTranspEqVector[iScal].phi[icv];
		}
	}

	double arcLengthNew = calcArclength(q1, qVec, lambda1, lambdaNew, NcontrolEqns, nScal, weightLambda, q_tangent, lambda_tangent);

	// Error check
	if(isNaN(arcLengthNew)) {
		if(mpi_rank==0) cerr<<"Error in IkeWithPsALC_AD::calcNresForJOE(): The new arclength calculated for the tangential eqns is NaN"<<endl;
		delete [] qVec;
		throw(PSALC_ERROR_CODE);
	}

	// Calculate Nres = (new arclength) - (old arclength)
	double NresVal = arcLength - arcLengthNew;

//IKJ
NresVal *= AreaLambda;

	for(int iParam=0; iParam<NcontrolEqns; ++iParam)
		Nres[iParam] = NresVal;

	// Free memory
	delete [] qVec;

	// Return true if the tangential condition is satisfied, i.e., check if (new arclength) - (old arclength) == 0
	if(NresVal > 1.0e3*double(NcontrolEqns*cvora[mpi_size])*MACHINE_EPS) {
		return false;
	}

	return true;
}

/*
 * Method: updateTotResidualWithNres
 * ---------------------------------
 * Update total residual norm with the Nres array
 */
double IkeWithPsALC_AD::updateTotResidualWithNres(double residNormTot, const double* Nres, const int NcontrolEqns, const int whichNorm) {
	assert(NcontrolEqns>0 && Nres!=NULL);

	// update residNorm to include Nres
	double residTotNew;

	switch(whichNorm) {
	case 0:
		residTotNew = residNormTot;
		for(int i=0; i<NcontrolEqns; ++i) {
			if(residTotNew < Nres[i])
				residTotNew = Nres[i];
		}
		break;
	case 1:
		residTotNew = residNormTot;
		for(int i=0; i<NcontrolEqns; ++i)
			residTotNew += fabs(Nres[i]);
		break;
	case 2:
		residTotNew = pow(residNormTot, 2.0);
		for(int i=0; i<NcontrolEqns; ++i)
			residTotNew += pow(Nres[i], 2.0);
		residTotNew = sqrt(residTotNew);
		break;
	default:
		if(mpi_rank==0)
			cout<<"ERROR in IkeWithPsALC_AD::UpdateTotResidualWithNres(): unsupported norm = "<<whichNorm<<endl;
		throw(PSALC_ERROR_CODE);
		break;
	}

	return residTotNew;
}

/*
 * Method: write1DArrayWithXcvSerial
 * ---------------------------------
 * Write 1D data array (size = ncv*(5+nScal))
 * Format:
 *   Header:
 *       step (int),
 *       mpi_size (int), cvora (int[mpi_size+1]), xcvMinArray(double[mpi_size*3]), xcvMaxArray(double[mpi_size*3])
 *   Body: For each mpi_rank,
 *       x_cv
 *       qVec
 *   Foot:
 *       EOF_ERROR_CHECK_CODE
 */
void IkeWithPsALC_AD::write1DArrayWithXcvSerial(const string &filename, const int step, const double *qVec) {
	char *filenameArray = new char[filename.size() + 1]; // Since MPI_File_open() doesn't accept "const char *", you cannot use std::string.c_str()
	std::copy(filename.begin(), filename.end(), filenameArray);
	filenameArray[filename.size()] = '\0';
	write1DArrayWithXcvSerial(filenameArray, step, qVec);

	delete [] filenameArray;
}

void IkeWithPsALC_AD::write1DArrayWithXcvSerial(const char filename[], const int step, const double *qVec) {
	// Number of variables
	int nScal = scalarTranspEqVector.size();
	int m = nScal+5; // number of variables (rho, rhou, rhov, rhow, rhoE, and scalars)

	// xcvMin and xcvMax
	double *xcvMinArray = new double [mpi_size*3];
	double *xcvMaxArray = new double [mpi_size*3];
	MPI_Gather(xcvMin, 3, MPI_DOUBLE, xcvMinArray, 3, MPI_DOUBLE, 0, mpi_comm);
	MPI_Gather(xcvMax, 3, MPI_DOUBLE, xcvMaxArray, 3, MPI_DOUBLE, 0, mpi_comm);

	if(mpi_rank != 0) { MPI_Status status; int dummy; MPI_Recv(&dummy, 1, MPI_INT, mpi_rank-1, 829, mpi_comm, &status); }

	// 1. Open the file
	ofstream ofile;
	if(mpi_rank == 0)
		ofile.open(filename, ios_base::out | ios_base::trunc | ios_base::binary);
	else
		ofile.open(filename, ios_base::out | ios_base::app | ios_base::binary);

	// 2. Write data on the file
	int dummyInt;
	double dummyDouble;
	// 2-1. Header
	if(mpi_rank==0) {
		dummyInt=step;				ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int)); // Note: if you simply use the "<<" operator, it will cause some problem due to bugs in a g++ library on linux
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
	}
	// 2-2. Body
	for(int icv=0; icv<ncv; ++icv) {
		for(int i=0; i<3; ++i) {
			dummyDouble=x_cv[icv][i]; 	ofile.write(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
		}
	}
	for(int icv=0; icv<ncv; ++icv) {
		for(int i=0; i<m; ++i) {
			dummyDouble=qVec[m*icv+i]; 	ofile.write(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
		}
	}

	// 2-3. Foot
	if(mpi_rank==mpi_size-1) {
		dummyInt = EOF_ERROR_CHECK_CODE; 	ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
	}

	// 3. Close the file
	ofile.close();

	if(mpi_rank < mpi_size-1) { int dummy; MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 829, mpi_comm); }

	delete [] xcvMinArray;
	delete [] xcvMaxArray;

	MPI_Barrier(mpi_comm);
}

/*
 * Method: writePsALCdumpedDataSerial
 * ----------------------------------
 * Write data for an IKE simulation
 * Format:
 *   Header:
 *       step (int), NcontrolEqns (int), lambda0 (double[NcontrolEqns]), lambda1 (double[NcontrolEqns]), dsTry (double),
 *       mpi_size (int), cvora (int[mpi_size+1]), xcvMinArray(double[mpi_size*3]), xcvMaxArray(double[mpi_size*3])
 *   Body: For each mpi_rank,
 *       x_cv
 *       qVec
 *   Foot:
 *       EOF_ERROR_CHECK_CODE
 */
void IkeWithPsALC_AD::writePsALCdumpedDataSerial(const string &filename, const int step, const double dsTry, const double *lambda0, const double *lambda1, const int NcontrolEqns, const double *qVec) {
	char *filenameArray = new char[filename.size() + 1]; // Since MPI_File_open() doesn't accept "const char *", you cannot use std::string.c_str()
	std::copy(filename.begin(), filename.end(), filenameArray);
	filenameArray[filename.size()] = '\0';
	writePsALCdumpedDataSerial(filenameArray, step, dsTry, lambda0, lambda1, NcontrolEqns, qVec);

	delete [] filenameArray;
}

void IkeWithPsALC_AD::writePsALCdumpedDataSerial(const char filename[], const int step, const double dsTry, const double *lambda0, const double *lambda1, const int NcontrolEqns, const double *qVec) {
	// Number of variables
	int nScal = scalarTranspEqVector.size();
	int m = nScal+5; // number of variables (rho, rhou, rhov, rhow, rhoE, and scalars)

	// xcvMin and xcvMax
	double *xcvMinArray = new double [mpi_size*3];
	double *xcvMaxArray = new double [mpi_size*3];
	MPI_Gather(xcvMin, 3, MPI_DOUBLE, xcvMinArray, 3, MPI_DOUBLE, 0, mpi_comm);
	MPI_Gather(xcvMax, 3, MPI_DOUBLE, xcvMaxArray, 3, MPI_DOUBLE, 0, mpi_comm);

	if(mpi_rank != 0) { MPI_Status status; int dummy; MPI_Recv(&dummy, 1, MPI_INT, mpi_rank-1, 829, mpi_comm, &status); }

	// 1. Open the file
	ofstream ofile;
	if(mpi_rank == 0)
		ofile.open(filename, ios_base::out | ios_base::trunc | ios_base::binary);
	else
		ofile.open(filename, ios_base::out | ios_base::app | ios_base::binary);

	// 2. Write data on the file
	int dummyInt;
	double dummyDouble;
	// 2-1. Header: Structure of the header part
	//              1. step (int)                      2. NcontrolEqns (int)
	//              3. lambda0 (double*NcontrolEqns)   4. lambda1 (double*NcontrolEqns)
	//              5. ds (double)
	//              6. mpi_size (int)                  7. cvora (double*(mpi_size+1)
	//              8. xMinArray (double*3*mpi_size)   9. xMaxArray (double*3*mpi_size)
	if(mpi_rank==0) {
		dummyInt=step;				ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int)); // Note: if you simply use the "<<" operator, it will cause some problem due to bugs in a g++ library on linux
		dummyInt=NcontrolEqns;		ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
			dummyDouble=lambda0[iEqn];		ofile.write(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
		}
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
			dummyDouble=lambda1[iEqn];		ofile.write(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
		}
		dummyDouble=dsTry;			ofile.write(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
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
	}
	// 2-2. Body: Structure of the body part
	//            For each mpi,
	//              1. x_cv (double*3*ncv)
	//              2. qVec (double*m*ncv)
	for(int icv=0; icv<ncv; ++icv) {
		for(int i=0; i<3; ++i) {
			dummyDouble=x_cv[icv][i]; 	ofile.write(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
		}
	}
	for(int icv=0; icv<ncv; ++icv) {
		for(int i=0; i<m; ++i) {
			dummyDouble=qVec[m*icv+i]; 	ofile.write(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
		}
	}

	// 2-3. Foot
	if(mpi_rank==mpi_size-1) {
		dummyInt = EOF_ERROR_CHECK_CODE; 	ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
	}

	// 3. Close the file
	ofile.close();

	if(mpi_rank < mpi_size-1) { int dummy; MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 829, mpi_comm); }

	delete [] xcvMinArray;
	delete [] xcvMaxArray;

	MPI_Barrier(mpi_comm);
}

/*
 * Method: writePsALCdumpedDataParallel
 * ------------------------------------
 * Write data for an IKE simulation
 * Format:
 *   Header:
 *       step (int), NcontrolEqns (int), lambda0 (double[NcontrolEqns]), lambda1 (double[NcontrolEqns]), dsTry (double),
 *       mpi_size (int), cvora (int[mpi_size+1]), xcvMinArray(double[mpi_size*3]), xcvMaxArray(double[mpi_size*3])
 *   Body: For each mpi_rank,
 *       x_cv
 *       qVec
 *   Foot:
 *       EOF_ERROR_CHECK_CODE
 */
void IkeWithPsALC_AD::writePsALCdumpedDataParallel(const string &filename, const int step, const double dsTry, const double *lambda0, const double *lambda1, const int NcontrolEqns, double *qVec) {
	char *filenameArray = new char[filename.size() + 1]; // Since MPI_File_open() doesn't accept "const char *", you cannot use std::string.c_str()
	std::copy(filename.begin(), filename.end(), filenameArray);
	filenameArray[filename.size()] = '\0';

	writePsALCdumpedDataParallel(filenameArray, step, dsTry, lambda0, lambda1, NcontrolEqns, qVec);

	delete [] filenameArray;
}

void IkeWithPsALC_AD::writePsALCdumpedDataParallel(char filename[], const int step, const double dsTry, const double *lambda0, const double *lambda1, const int NcontrolEqns, double *qVec) {
	// Number of variables
	int nScal = scalarTranspEqVector.size();
	int m = nScal+5; // number of variables (rho, rhou, rhov, rhow, rhoE, and scalars)

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
	//              3. lambda0 (double*NcontrolEqns)   4. lambda1 (double*NcontrolEqns)
	//              5. ds (double)
	//              6. mpi_size (int)                  7. cvora (int*(mpi_size+1))
	//              8. xMinArray (double*3*mpi_size)   9. xMaxArray (double*3*mpi_size)
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
			dummyDouble=lambda0[iEqn];		ofile.write(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
		}
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
			dummyDouble=lambda1[iEqn];		ofile.write(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
		}
		dummyDouble=dsTry;			ofile.write(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
		dummyInt=mpi_size; 			ofile.write(reinterpret_cast<char*>(&dummyInt),    sizeof(int));
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
	int initDisp = sizeof(int)*(3 + mpi_size+1) + sizeof(double)*(NcontrolEqns*2 + 1 + 2*mpi_size*3);

	// 2. Body
	//       Structure of the body part:
	//            For each mpi,
	//              1. x_cv (double*3*ncv)
	//              2. qVec (double*m*ncv)
	MPI_Status status;
	MPI_File fh;
	MPI_Offset displacement;

	// Open the file
	if (MPI_File_open(mpi_comm, filename, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh)!=0) {
			// Note: MPI_MODE_WRONLY = write only
			//       MPI_MODE_CREATE = create the file if it does not exist.
		stringstream ss;
		ss<<"ERROR! IkeWithPsALC_AD::writePsALCdumpedDataParallel(): Cannot open "<<filename<<endl;
		errorMessage = ss.str();
		showMessageParallel(errorMessage);

		delete [] xcvMinArray; 	delete [] xcvMaxArray;
		throw(PSALC_ERROR_CODE);
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
		stringstream ss;
		ss<<"ERROR! IkeWithPsALC_AD::writePsALCdumpedDataParallel(): Cannot set the first MPI_File_set_view -- offset="<<displacement<<endl;
		errorMessage = ss.str();
		showMessageParallel(errorMessage);

		delete [] xcvMinArray; 	delete [] xcvMaxArray;
		throw(PSALC_ERROR_CODE);
	}
	double* bufferDouble = new double [ncv*3];
	for(int icv=0; icv<ncv; ++icv) {
		int indexStart = icv*3;
		for(int i=0; i<3; ++i)
			bufferDouble[indexStart+i] = x_cv[icv][i];
	}
	MPI_File_write(fh, bufferDouble, ncv*3, MPI_DOUBLE, &status);
	delete [] bufferDouble; 	bufferDouble = NULL;

	// Write the qVec vector
	displacement += MPI_Offset(3*ncv*sizeof(double));
	MPI_File_write(fh, qVec, ncv*m, MPI_DOUBLE, &status);

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
		cout<<"> The Q vector is written on "<<filename<<" (RUN TIME = "<<wtimeF-wtime0<<" [sec]) \n"<<endl;

	/***********
	 ** Free the memory
	 ***********/
	delete [] xcvMinArray;
	delete [] xcvMaxArray;
}

/*
 * Method: readPsALCdumpedDataSerial
 * ---------------------------------
 * Read data from previous IKE simulation
 */
void IkeWithPsALC_AD::readPsALCdumpedDataSerial(const string &filename, double* qVec, int& step0, double* lambda0, double* lambda1, const int NcontrolEqns, double& dsTry) {
	char *filenameArray = new char[filename.size() + 1]; // Since MPI_File_open() doesn't accept "const char *", you cannot use std::string.c_str()
	std::copy(filename.begin(), filename.end(), filenameArray);
	filenameArray[filename.size()] = '\0';
	readPsALCdumpedDataSerial(filenameArray, qVec, step0, lambda0, lambda1, NcontrolEqns, dsTry);

	delete [] filenameArray;
}

void IkeWithPsALC_AD::readPsALCdumpedDataSerial(const char filename[], double* qVec, int& step0, double* lambda0, double* lambda1, const int NcontrolEqns, double& dsTry) {
	// Check some possible errors
	for(int i=0; i<3; ++i)
		if(xcvMax[i] < xcvMin[i]) {
			cerr<<"ERROR in IkeWithPsALC_AD::readPsALCdumpedDataSerial(): mpi_rank == "<<mpi_rank<<" : xcvMax["<<i<<"] < xcvMin["<<i<<"]"<<endl;
			throw(PSALC_ERROR_CODE);
		}

	// number of variables
	int nScal = scalarTranspEqVector.size();
	int m = nScal+5; // number of variables (rho, rhou, rhov, rhow, rhoE, and scalars)

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
				printf("ERROR in IkeWithPsALC_AD::readPsALCdumpedData: NcontrolEqns from file=%d is not equal to given NcontrolEqns=%d\n", NcontrolEqns_file, NcontrolEqns);
			delete [] countFound;
			throw(PSALC_ERROR_CODE);
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
				cout<<"ERROR! IkeWithPsALC_AD::readPsALCdumpedDataSerial(): file does not end with EOF_ERROR_CHECK_CODE="<<EOF_ERROR_CHECK_CODE<<endl;
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

	MPI_Barrier(mpi_comm);
}

/*
 * Method: readJOEdumpedDataSerial
 * -------------------------------
 * Read data from previous JOE simulation
 */
void IkeWithPsALC_AD::readJOEdumpedDataSerial(const string &filename, double* q) {
	// number of variables
	int nScal = scalarTranspEqVector.size();
	int m = nScal+5; // number of variables (rho, rhou, rhov, rhow, rhoE, and scalars)

	// variables related to the Mesh
	const double gridTol = GRID_TOL;

	int *countFound = new int [ncv];
	for(int icv=0; icv<ncv; ++icv)
		countFound[icv] = 0;

	// read the file
	if(mpi_rank==0)
		cout<<">> readJOEdumpedDataSerial(): ";

	int mpi_size_file;
	int *cvora_file = NULL;
	double *xcvMinArray_file = NULL;
	double *xcvMaxArray_file = NULL;
	int nScal_file;

	ifstream infile(filename.c_str(), ios::in | ios::binary);
	if(infile.is_open()) {
		// read some parameters from the file
		infile.read(reinterpret_cast<char*>(&mpi_size_file), sizeof(int)); // Note: if you simply use the ">>" operator, it will cause some problem due to bugs in a g++ library on linux

		cvora_file = new int [mpi_size_file+1];
		infile.read(reinterpret_cast<char*>(cvora_file), sizeof(int)*(mpi_size_file+1));

		xcvMinArray_file = new double [mpi_size_file*3];
		infile.read(reinterpret_cast<char*>(xcvMinArray_file), sizeof(double)*(mpi_size_file*3));
		xcvMaxArray_file = new double [mpi_size_file*3];
		infile.read(reinterpret_cast<char*>(xcvMaxArray_file), sizeof(double)*(mpi_size_file*3));

		infile.read(reinterpret_cast<char*>(&nScal_file), sizeof(int));

		if(mpi_rank==0) {
			printf("Reading \'%s\': total number of CVs=%d with # of fields=%d and nScal=%d ...\n", filename.c_str(), cvora_file[mpi_size_file], mpi_size_file, nScal_file);
		}
		assert(nScal_file == nScal);
		assert(cvora_file[mpi_size_file] == cvora[mpi_size]);

		// variables
		int ncv_file;
		double (*x_cv_file)[3] = NULL;
		double *q_file = NULL;

		// read each data
		for(int impi_file=0; impi_file<mpi_size_file; ++impi_file) {
			// read data dumped by a CPU
			ncv_file = cvora_file[impi_file+1]-cvora_file[impi_file];
			assert(ncv_file>0 && ncv_file<=cvora[mpi_size]);
			if(mpi_rank==0)
				cout<<"  Data field "<<impi_file<<": ncv = "<<ncv_file<<endl;

			assert(x_cv_file == NULL);	x_cv_file = new double [ncv_file][3];
			assert(q_file == NULL);	    q_file    = new double [ncv_file*m];

			for(int icv=0; icv<ncv_file; ++icv) {
				for(int i=0; i<3; ++i)
					infile.read(reinterpret_cast<char*>(&x_cv_file[icv][i]), sizeof(double));
			}
			infile.read(reinterpret_cast<char*>(q_file), sizeof(double)*(ncv_file*m));

			// find the matached data
			for(int icv=0; icv<ncv; ++icv) {
				int icv_file;
				int count = findMatchedIndex(icv_file, icv, ncv_file, x_cv_file, gridTol);
				if(count>0) {
					for(int i=0; i<m; ++i)
						q[m*icv+i] = q_file[m*icv_file+i];
					countFound[icv] += count;
				}
			}

			delete [] x_cv_file;	x_cv_file = NULL;
			delete [] q_file;		q_file = NULL;
		}

		// Read the foot
		int dummyInt;
		infile.read(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		if(dummyInt != EOF_ERROR_CHECK_CODE) {
			if(mpi_rank==0)
				cout<<"ERROR! IkeWithPsALC_AD::readJOEdumpedDataSerial(): file does not end with EOF_ERROR_CHECK_CODE="<<EOF_ERROR_CHECK_CODE<<endl;
			assert(false);
		}

		// check if each CV finds a matched CV
		for (int icv=0; icv<ncv; ++icv) {
			if(countFound[icv]==0) {
				printf("  Cannot find a matched point for icv=%d (%.6f,%.6f,%.6f) at mpi_rank=%d \n", icv, x_cv[icv][0],x_cv[icv][1],x_cv[icv][2], mpi_rank);
				throw(-1);
			} else if(countFound[icv]>1) {
				printf("  Too many matched points(%d) for icv=%d (%.6f,%.6f,%.6f) at mpi_rank=%d \n", countFound[icv], icv, x_cv[icv][0],x_cv[icv][1],x_cv[icv][2], mpi_rank);
				throw(-1);
			}
		}
	} else {
		if(mpi_rank==0)
			printf("Cannot read \'%s\' \n", filename.c_str());
		throw(-1);
	}

	// close the file
	infile.close();

	delete [] countFound;
	delete [] cvora_file;
	delete [] xcvMinArray_file;
	delete [] xcvMaxArray_file;

	MPI_Barrier(mpi_comm);
}

/*
 * Method: findMatchedIndex
 * ------------------------
 * This method is used to find a matched index for the readPsALCdumpedData() method and the readJOEdumpedData() method
 */
int IkeWithPsALC_AD::findMatchedIndex(int& foundIndex, const int icv, const int ncv_file, const double (*x_cv_file)[3], const double gridTol) {
	int count = 0;
	foundIndex = -1;
	for(int index =0; index<ncv_file; ++index) {
		if( fabs(x_cv_file[index][0]-x_cv[icv][0])<gridTol && fabs(x_cv_file[index][1]-x_cv[icv][1])<gridTol && fabs(x_cv_file[index][2]-x_cv[icv][2])<gridTol ) {
			foundIndex = index;
			count++;
		}
	}
	return count;
}

/*
 * Method: writeMatrixMfile
 * ------------------------
 * write the matrix on files as a matlab format (mfile format)
 * Sometimes, the size of the matrix is too large to be saved on a single file.
 */
template <class MatT>
void IkeWithPsALC_AD::writeMatrixMfile(string &filename, MatT &jacMatrix, const int nScal, const int *cvora, const int *cv_gl) {
	assert(NcontrolEqns > 0); // If you delete this line, you will have some compile errors. I don't have enough time to fix this right now.

	int m = 5 + nScal;

	// get nCols, nRows, nnz
	int nCols = jacMatrix.get_nCols(); // Each CPU might have different number of columns
	int nRows = jacMatrix.get_nRows();
	int nnz = jacMatrix.get_nnz();

	// open the file
	FILE *fp;
	stringstream ss;
	string newFilename;
	ss<<filename<<"_"<<mpi_rank<<".m";
	newFilename = ss.str();

	fp = fopen(newFilename.c_str(), "w");

	// write total number of global Rows, global Columns, and global number of non-zeros at the beginning of the file
	fprintf(fp, "function [nRows, nCols, nnz, A] = %s_%d() \n", filename.c_str(),mpi_rank);
	fprintf(fp, "nRows = %d; \n",nRows);
	fprintf(fp, "nCols = %d; \n",nCols);
	fprintf(fp, "nnz = %d; \n",nnz);

	fprintf(fp, "A = [ \n");

	// write row index, column index, and values as a MATLAB sparse matrix format
	for(int i=0; i<nnz; ++i) {
		int globalCind = jacMatrix.get_global_cind(i, mpi_rank, m, cvora, cv_gl, NcontrolEqns, ncv_gg);

		int globalRind = jacMatrix.get_global_rind(i, mpi_rank, m, cvora);

		fprintf(fp, "%d\t%d\t%.8e \n", globalRind+1, globalCind+1, jacMatrix.get_values(i)); // MATLAB index starts from 1 instead of 0
	}

	fprintf(fp, "]; \n");
	fprintf(fp, "end \n\n");

	// close the file
	fclose(fp);

	MPI_Barrier(mpi_comm);
}

/*
 * Method: writeMatrixBinarySerial
 * --------------------------------
 * Write the matrix on a binary file in a serial fashion (i.e. One CPU writes its data while the other CPUs wait)
 * File format:
 *   Header:
 *     nRows(int), nCols(int), nnz(int), NcontrolEqns(int)
 *     mpi_size, localNnzArray(int [])
 *   Matrix:
 *     global_rind(int) vector
 *     global_cind(int) vector
 *     values(double) vector
 *   Foot:
 *     EOF_ERROR_CHECK_CODE
 */
template <class MatT>
void IkeWithPsALC_AD::writeMatrixBinarySerial(string &filename, MatT &jacMatrix, const int nScal, const int *cvora, const int *cv_gl) {
	assert(NcontrolEqns > 0); // If you delete this line, you will have some compile errors. I don't have enough time to fix this right now.

	int m = 5 + nScal;

	/*
	 * Get nRows, nCols, nnz
	 */
	int myNRows = jacMatrix.get_nRows();
//	unsigned long int myNnz = (unsigned long int) jacMatrix.get_nnz();
	int myNnz = jacMatrix.get_nnz(); // Note: Since nnz can be a gigantic number, it is better to use "unsigned int" instead of "int".
	                                 //       However, "unsigned int" gives an weird value after MPI_Allreduce and contaminates the following data for some reasons.

	int nRows;
	MPI_Allreduce(&myNRows, &nRows, 1, MPI_INT, MPI_SUM, mpi_comm);
//	unsigned long int nnz;
	int nnz;
//	MPI_Allreduce(&myNnz, &nnz, 1, MPI_UNSIGNED, MPI_SUM, mpi_comm);
	MPI_Allreduce(&myNnz, &nnz, 1, MPI_INT, MPI_SUM, mpi_comm);

	int *localNnzArray = NULL;
	if(mpi_rank == 0)
		localNnzArray = new int [mpi_size];
	MPI_Gather (&myNnz, 1, MPI_INT, localNnzArray, 1, MPI_INT, 0, mpi_comm);

	/*
	 * Write the data
	 */
	ofstream ofile;
	int dummyInt;
//	unsigned long int dummyUint;
	double dummyDouble;

	// Write the row indices first (Write the header if mpi_rank==0)
	if(mpi_rank != 0) { MPI_Status status; int dummy; MPI_Recv(&dummy, 1, MPI_INT, mpi_rank-1, 430, mpi_comm, &status); }
	if(mpi_rank == 0) {
		ofile.open(filename.c_str(), ios_base::out | ios_base::trunc | ios_base::binary);

		// Write the header
		dummyInt = nRows; 			ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int)); // Note: if you simply use the "<<" operator, it will cause some problem due to bugs in a g++ library on linux
		dummyInt = cvora[mpi_size]*m+NcontrolEqns; 			ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
//		dummyUint = nnz; 			ofile.write(reinterpret_cast<char*>(&dummyUint), sizeof(unsigned long int));
		dummyInt = nnz; 			ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		dummyInt = NcontrolEqns; 	ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		dummyInt = mpi_size; 								ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		for(int i=0; i<mpi_size; ++i)
			dummyInt = localNnzArray[i]; 					ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
	} else
		ofile.open(filename.c_str(), ios_base::out | ios_base::app | ios_base::binary);

	for(int i=0; i<(int)myNnz; ++i) {
		dummyInt = jacMatrix.get_global_rind(i, mpi_rank, m, cvora);
		ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
	}

	ofile.close();
	if(mpi_rank < mpi_size-1) { int dummy; MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 430, mpi_comm); }
	MPI_Barrier(mpi_comm);

	// Write the column indices
	if(mpi_rank != 0) { MPI_Status status; int dummy; MPI_Recv(&dummy, 1, MPI_INT, mpi_rank-1, 431, mpi_comm, &status); }
	ofile.open(filename.c_str(), ios_base::out | ios_base::app | ios_base::binary);

	for(int i=0; i<(int)myNnz; ++i) {
		dummyInt = jacMatrix.get_global_cind(i, mpi_rank, m, cvora, cv_gl, NcontrolEqns, ncv_gg);
		ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
	}

	ofile.close();
	if(mpi_rank < mpi_size-1) { int dummy; MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 431, mpi_comm); }
	MPI_Barrier(mpi_comm);

	// Write the values
	if(mpi_rank != 0) { MPI_Status status; int dummy; MPI_Recv(&dummy, 1, MPI_INT, mpi_rank-1, 432, mpi_comm, &status); }
	ofile.open(filename.c_str(), ios_base::out | ios_base::app | ios_base::binary);

	for(int i=0; i<(int)myNnz; ++i) {
		dummyDouble = jacMatrix.get_values(i);
		ofile.write(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
	}

	if(mpi_rank==mpi_size-1) {
		dummyInt = EOF_ERROR_CHECK_CODE; // The last element is a pre-defined integer value --> For the Error check
		ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
	}

	ofile.close();
	if(mpi_rank < mpi_size-1) { int dummy; MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 432, mpi_comm); }
	MPI_Barrier(mpi_comm);
}

/*
 * Method: writeMatrixBinaryParallel
 * ---------------------------------
 * Write the matrix on a binary file using MPI-2 parallel I/O
 * File format:
 *   Header:
 *     nRows(int), nCols(int), nnz(int), NcontrolEqns(int)
 *     mpi_size, localNnzArray(int [])
 *   Matrix:
 *     global_rind(int) vector
 *     global_cind(int) vector
 *     values(double) vector
 *   Foot:
 *     EOF_ERROR_CHECK_CODE
 */
template <class MatT>
void IkeWithPsALC_AD::writeMatrixBinaryParallel(string &filename, MatT &jacMatrix, const int nScal, const int *cvora, const int *cv_gl) {
	assert(NcontrolEqns > 0); // If you delete this line, you will have some compile errors. I don't have enough time to fix this right now.

	int m = 5 + nScal;

	char *filenameArray = new char[filename.size() + 1]; // Since MPI_File_open() doesn't accept "const char *", you cannot use std::string.c_str()
	std::copy(filename.begin(), filename.end(), filenameArray);
	filenameArray[filename.size()] = '\0';

	string errorMessage;

	/***********
	 ** Initialization: Get nRows, nCols, nnz, myNnzArray
	 ***********/
	int myNRows = jacMatrix.get_nRows();
	int myNnz = jacMatrix.get_nnz();

	int nRows;
	MPI_Allreduce(&myNRows, &nRows, 1, MPI_INT, MPI_SUM, mpi_comm);
	int nnz;
	MPI_Allreduce(&myNnz, &nnz, 1, MPI_INT, MPI_SUM, mpi_comm);

	int *localNnzArray = NULL;
	if(mpi_rank == 0)
		localNnzArray = new int [mpi_size];
	MPI_Gather(&myNnz, 1, MPI_INT, localNnzArray, 1, MPI_INT, 0, mpi_comm);

	/***********
	 ** Write the data
	 ***********/
	double wtime0, wtimeF;
	if(mpi_rank==0)
		wtime0 = MPI_Wtime();

	/*
	 * 1. Write the header part
	 */
	ofstream ofile;
	int dummyInt;

	if(mpi_rank == 0) {
		ofile.open(filenameArray, ios_base::out | ios_base::trunc | ios_base::binary);

		// Write the header (number of global rows, number of global columns, global nnz, NcontrolEqns)
		dummyInt = nRows; 									ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int)); // Note: if you simply use the "<<" operator, it will cause some problem due to bugs in a g++ library on linux
		dummyInt = cvora[mpi_size]*m+NcontrolEqns; 			ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		dummyInt = nnz; 									ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		dummyInt = NcontrolEqns; 							ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		dummyInt = mpi_size; 								ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		for(int i=0; i<mpi_size; ++i) {
			dummyInt = localNnzArray[i]; 					ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		}
		ofile.close();
	}
	MPI_Barrier(mpi_comm);
	int initDisp = sizeof(int)*(5 + mpi_size);

	/*
	 * 2. Write the body (matrix)
	 */
	MPI_Status status;
	MPI_File fh;
	MPI_Offset displacement;

	// Open the file
	if (MPI_File_open(mpi_comm, filenameArray, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh)!=0) {
			// Note: MPI_MODE_WRONLY = write only
			//       MPI_MODE_CREATE = create the file if it does not exist.
		stringstream ss;
		ss<<"ERROR! IkeWithPsALC_AD::writeMatrixBinaryParallel(): Cannot open "<<filenameArray;
		errorMessage = ss.str();
		showMessageParallel(errorMessage);

		delete [] filenameArray;
		throw(PSALC_ERROR_CODE);
	}

	// 2-1. Write the row indices first
	int myOffsetNnz;
	MPI_Scan(&myNnz, &myOffsetNnz, 1, MPI_INT, MPI_SUM, mpi_comm);
	myOffsetNnz -= myNnz;
	displacement = MPI_Offset(initDisp + sizeof(int)*myOffsetNnz);
	if(MPI_File_set_view(fh, displacement, MPI_INT, MPI_INT, "native", MPI_INFO_NULL)!=0) {
		// Note: native     = Data in this representation are stored in a file exactly as it is in memory
		//       Internal   = Data in this representation can be used for I/O operations in a homogeneous or heterogeneous environment
		//       External32 = This data representation states that read and write operations convert all data from and to the "external32" representation
		stringstream ss;
		ss<<"ERROR! IkeWithPsALC_AD::writeMatrixBinaryParallel(): Cannot open "<<filenameArray;
		errorMessage = ss.str();
		showMessageParallel(errorMessage);
		delete [] filenameArray;
		if(localNnzArray != NULL) { delete [] localNnzArray; 	localNnzArray = NULL; }
		throw(PSALC_ERROR_CODE);
	}

	int *bufferInt = new int [myNnz];
	for(int i=0; i<myNnz; ++i)
		bufferInt[i] = jacMatrix.get_global_rind(i, mpi_rank, m, cvora);
	MPI_File_write(fh, bufferInt, myNnz, MPI_INT, &status);
	delete [] bufferInt; bufferInt = NULL;

	// 2-2. Write the column indices
	displacement += MPI_Offset(nnz*sizeof(int));
	if(MPI_File_set_view(fh, displacement, MPI_INT, MPI_INT, "native", MPI_INFO_NULL)!=0) {
		stringstream ss;
		ss<<"ERROR! IkeWithPsALC_AD::writeMatrixBinaryParallel(): Cannot set the second MPI_File_set_view -- offset="<<displacement<<endl;
		errorMessage = ss.str();
		showMessageParallel(errorMessage);

		delete [] filenameArray;
		if(localNnzArray != NULL) { delete [] localNnzArray; 	localNnzArray = NULL; }
		throw(PSALC_ERROR_CODE);
	}

	bufferInt = new int [myNnz];
	for(int i=0; i<myNnz; ++i)
		bufferInt[i] = jacMatrix.get_global_cind(i, mpi_rank, m, cvora, cv_gl, NcontrolEqns, ncv_gg);
	MPI_File_write(fh, bufferInt, myNnz, MPI_INT, &status);
	delete [] bufferInt; bufferInt = NULL;

	// 2-3. Write the values
	displacement += MPI_Offset((nnz-myOffsetNnz)*sizeof(int) + myOffsetNnz*sizeof(double));
	if(MPI_File_set_view(fh, displacement, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL)!=0) {
		stringstream ss;
		ss<<"ERROR! IkeWithPsALC_AD::writeMatrixBinaryParallel(): Cannot set the third MPI_File_set_view -- offset="<<displacement<<endl;
		errorMessage = ss.str();
		showMessageParallel(errorMessage);

		delete [] filenameArray;
		if(localNnzArray != NULL) { delete [] localNnzArray; 	localNnzArray = NULL; }
		throw(PSALC_ERROR_CODE);
	}

	double *bufferDouble = new double [myNnz];
	for(int i=0; i<myNnz; ++i)
		bufferDouble[i] = jacMatrix.get_values(i);
	MPI_File_write(fh, bufferDouble, myNnz, MPI_DOUBLE, &status);
	delete [] bufferDouble; bufferDouble = NULL;

	// 2-4. Close the file
	MPI_File_close(&fh);
	MPI_Barrier(mpi_comm);

	/*
	 * 3. Write the foot
	 */
	if(mpi_rank == mpi_size-1) {
		ofile.open(filenameArray, ios_base::out | ios_base::app | ios_base::binary);
		dummyInt = EOF_ERROR_CHECK_CODE; 	ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		ofile.close();
	}
	MPI_Barrier(mpi_comm);

	if(mpi_rank==0)
		wtimeF = MPI_Wtime();

	/***********
	 ** Show the summary on the screen if the debugging level is high
	 ***********/
	if(debugLevel>0 && mpi_rank==0)
		printf("     > The Jacobian matrix is written on %s (RUN TIME = %.3f [sec]): nRows=%d, nCols=%d, nnz=%d \n",
				filenameArray, wtimeF-wtime0, nRows, cvora[mpi_size]*m+NcontrolEqns, nnz);

	/***********
	 ** Free the memory
	 ***********/
	delete [] filenameArray;
	if(localNnzArray != NULL)
		delete [] localNnzArray;
}

/*
 * Method: readMatrixBinary
 * ------------------------
 * Read a matrix from a binary file and store it in a matrix container
 */
template <class MatT>
void IkeWithPsALC_AD::readMatrixBinary(string &filename, MatT &jacMatrix, const int nScal, const int *cvora, const int *cv_gl) {
	int nVars = 5 + nScal;

	/***********
	 ** Read the data
	 ***********/
	int myNnz_fromFile;
	vector<unsigned int> rind;
	vector<unsigned int> cind;
	vector<double> values;

	readMatrixBinaryParallel(myNnz_fromFile, rind, cind, values, filename, nVars, cvora, cv_gl);

	/***********
	 ** Store the matrix in the data container
	 ***********/
	assert(jacMatrix.empty());

	if(debugLevel>0 && mpi_rank==0)
		cout<<"    > Storing the matrix in the data container ... ";

	jacMatrix.setMatSize(ncv_gg*nVars+NcontrolEqns, ncv*nVars);

	for(int i=0; i<myNnz_fromFile; ++i) {
		rind[i] = jacMatrix.get_local_rind(rind[i], mpi_rank, nVars, cvora);
		cind[i] = jacMatrix.get_local_cind(cind[i], mpi_rank, nVars, cvora, cv_gl, NcontrolEqns, ncv, ncv_gg);
	}

	jacMatrix.matCopy(myNnz_fromFile, rind, cind, values);

	if(debugLevel>0 && mpi_rank==0)
		cout<<"DONE!"<<endl;

	/***********
	 ** Clear the memory
	 ***********/
	rind.clear();
	cind.clear();
	values.clear();

	MPI_Barrier(mpi_comm);
}

/*
 * Method: readMatrixBinaryParallel
 * --------------------------------
 * Read a matrix from a binary file using MPI-2 parallel I/O
 * Output: myNnz_fromFile = local number of non-zeros (NOT the total nnz -- it is the nnz for a core)
 *         rind           = GLOBAL row indices array
 *         cind           = GLOBAL column indices array
 *         values         = values array
 */
void IkeWithPsALC_AD::readMatrixBinaryParallel(int &myNnz_fromFile, vector<unsigned int> &rind, vector<unsigned int> &cind, vector<double> &values,
		string &filename, const int nVars, const int *cvora, const int *cv_gl) {
	/***********
	 ** Read the data
	 ***********/
	char *filenameArray = new char[filename.size() + 1]; // Since MPI_File_open() doesn't accept "const char *", you cannot use std::string.c_str()
	std::copy(filename.begin(), filename.end(), filenameArray);
	filenameArray[filename.size()] = '\0';

	/*
	 * 1. Read the header part
	 */
	int nRows_fromFile = -1;
	int nCols_fromFile = -1;
	int nnz_fromFile = -1;
	int NcontrolEqns_fromFile = -1;
	int mpi_size_fromFile = -1;
	int *localNnzArray_fromFile = NULL;

	ifstream infile(filename.c_str(), ios::in | ios::binary);
	if(infile.is_open()) {
		infile.read(reinterpret_cast<char*>(&nRows_fromFile), sizeof(int)); // Note: if you simply use the ">>" operator, it will cause some problem due to bugs in a g++ library on linux
		infile.read(reinterpret_cast<char*>(&nCols_fromFile), sizeof(int));
		infile.read(reinterpret_cast<char*>(&nnz_fromFile),   sizeof(int));
		infile.read(reinterpret_cast<char*>(&NcontrolEqns_fromFile), sizeof(int));
		if(NcontrolEqns_fromFile != NcontrolEqns) {
			if(mpi_rank==0)
				printf("ERROR in IkeWithPsALC_AD::readMatrixBinaryParallel: NcontrolEqns from file=%d is not equal to given NcontrolEqns=%d\n", NcontrolEqns_fromFile, NcontrolEqns);
			delete [] filenameArray;
			throw(PSALC_ERROR_CODE);
		}

		infile.read(reinterpret_cast<char*>(&mpi_size_fromFile), sizeof(int));
		assert(mpi_size_fromFile>0);
		if(mpi_size_fromFile != mpi_size) {
			if(mpi_rank==0)
				cerr<<"ERROR! IkeWithPsALC_AD::readMatrixBinaryParallel: mpi_size_fromFile(="<<mpi_size_fromFile<<") is not equal to mpi_size(="<<mpi_size<<")"<<endl;
			delete [] filenameArray;
			throw(PSALC_ERROR_CODE);
		}
		localNnzArray_fromFile = new int [mpi_size_fromFile];
		infile.read(reinterpret_cast<char*>(localNnzArray_fromFile), sizeof(int)*mpi_size_fromFile);
		if(debugLevel>0 && mpi_rank==0)
			cout<<"    > Reading the Jacobian matrix from "<<filename<<": nRows="<<nRows_fromFile<<", nCols="<<nCols_fromFile<<", nnz="<<nnz_fromFile<<endl;
		if(debugLevel>0 && mpi_rank==0) {
			cout<<"      > local nnz =";
			for(int i=0; i<mpi_size_fromFile; ++i)
				cout<<" "<<localNnzArray_fromFile[i]<<" ";
			cout<<endl;
		}

		assert(nRows_fromFile>0 && nCols_fromFile>0 && nnz_fromFile>0);
		assert(localNnzArray_fromFile[mpi_rank] > 0);
	} else {
		if(mpi_rank==0)
			printf("Cannot read \'%s\' \n", filename.c_str());
		assert(false);
	}
	infile.close();
	MPI_Barrier(mpi_comm);

	myNnz_fromFile = localNnzArray_fromFile[mpi_rank];
	int initDisp = sizeof(int)*(5 + mpi_size_fromFile);

	/*
	 * 2. Read the body (matrix)
	 */
	MPI_Status status;
	MPI_File fh;
	MPI_Offset displacement;

	// Open the file
	if (MPI_File_open(mpi_comm, filenameArray, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh)!=0) {
			// Note: MPI_MODE_RDONLY = read only
		cerr<<"ERROR! IkeWithPsALC_AD::readMatrixBinary(): Cannot open "<<filenameArray<<endl;
		delete [] filenameArray;
		delete [] localNnzArray_fromFile;
		throw(PSALC_ERROR_CODE);
	}

	// 2-1. Read the row indices
	int myOffsetNnz = 0;
	for(int i=0; i<mpi_rank; ++i)
		myOffsetNnz += localNnzArray_fromFile[i];
	displacement = MPI_Offset(initDisp + sizeof(int)*myOffsetNnz);
	if(MPI_File_set_view(fh, displacement, MPI_INT, MPI_INT, "native", MPI_INFO_NULL)!=0) {
		// Note: native     = Data in this representation are stored in a file exactly as it is in memory
		//       Internal   = Data in this representation can be used for I/O operations in a homogeneous or heterogeneous environment
		//       External32 = This data representation states that read and write operations convert all data from and to the "external32" representation
		cerr<<"ERROR! IkeWithPsALC_AD::readMatrixBinary(): Cannot set the first MPI_File_set_view -- offset="<<displacement<<endl;
		delete [] filenameArray;
		delete [] localNnzArray_fromFile;
		throw(PSALC_ERROR_CODE);
	}

	unsigned int *bufferRind = new unsigned int [myNnz_fromFile];
	MPI_File_read(fh, bufferRind, myNnz_fromFile, MPI_INT, &status);

	// 2-2. Read the column indices
	displacement += MPI_Offset(nnz_fromFile*sizeof(int));
	if(MPI_File_set_view(fh, displacement, MPI_INT, MPI_INT, "native", MPI_INFO_NULL)!=0) {
		cerr<<"ERROR! IkeWithPsALC_AD::readMatrixBinary(): Cannot set the first MPI_File_set_view -- offset="<<displacement<<endl;
		delete [] filenameArray;
		delete [] localNnzArray_fromFile;
		delete [] bufferRind;
		throw(PSALC_ERROR_CODE);
	}

	unsigned int *bufferCind = new unsigned int [myNnz_fromFile];
	MPI_File_read(fh, bufferCind, myNnz_fromFile, MPI_INT, &status);

	// 2-3. Read the values
	displacement += MPI_Offset((nnz_fromFile-myOffsetNnz)*sizeof(int) + myOffsetNnz*sizeof(double));
	if(MPI_File_set_view(fh, displacement, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL)!=0) {
		cerr<<"ERROR! IkeWithPsALC_AD::readMatrixBinary(): Cannot set the first MPI_File_set_view -- offset="<<displacement<<endl;
		delete [] filenameArray;
		delete [] localNnzArray_fromFile;
		delete [] bufferRind;
		delete [] bufferCind;
		throw(PSALC_ERROR_CODE);
	}

	double *bufferValues = new double [myNnz_fromFile];
	MPI_File_read(fh, bufferValues, myNnz_fromFile, MPI_DOUBLE, &status);

	/*
	 * 3. Read the foot
	 */
	displacement = (MPI_Offset) (initDisp + sizeof(int)*nnz_fromFile*2 + + sizeof(double)*nnz_fromFile);
	if(MPI_File_set_view(fh, displacement, MPI_INT, MPI_INT, "native", MPI_INFO_NULL)!=0) {
		cerr<<"ERROR! IkeWithPsALC_AD::readMatrixBinary(): Cannot set the first MPI_File_set_view -- offset="<<displacement<<endl;
		delete [] filenameArray;
		delete [] localNnzArray_fromFile;
		delete [] bufferRind;
		delete [] bufferCind;
		delete [] bufferValues;
		throw(PSALC_ERROR_CODE);
	}

	int dummyInt;
	MPI_File_read(fh, &dummyInt, 1, MPI_INT, &status);
	if(dummyInt != EOF_ERROR_CHECK_CODE) {
		cerr<<"ERROR! IkeWithPsALC_AD::readMatrixBinary(): The foot of the file (="<<dummyInt<<") is not equal to "<<EOF_ERROR_CHECK_CODE<<endl;
		delete [] filenameArray;
		delete [] localNnzArray_fromFile;
		delete [] bufferRind;
		delete [] bufferCind;
		delete [] bufferValues;
		throw(PSALC_ERROR_CODE);
	}

	/*
	 * Close the file and clear some arrays
	 */
	MPI_File_close(&fh);

	rind.resize(myNnz_fromFile);
	cind.resize(myNnz_fromFile);
	values.resize(myNnz_fromFile);
	for(int i=0; i<myNnz_fromFile; ++i) {
		rind[i] = bufferRind[i];
		cind[i] = bufferCind[i];
		values[i] = bufferValues[i];
	}

	delete [] filenameArray;
	delete [] localNnzArray_fromFile;
	delete [] bufferRind;
	delete [] bufferCind;
	delete [] bufferValues;

	MPI_Barrier(mpi_comm);
}

/*
 * Method: writeCompatibilityEveryPtsOnTxtFile
 * -------------------------------------------
 * Compare RHS from IKE and JOE, and write the difference if the difference is greater than threshold.
 * This is useful when Tecplot is not avilable or when you find the ICVs where the differnce is huge.
 * Filename is always ERROR_IKE_JOE.csv
 * The format of the file is "mpi_rank, icv, x_cv[0], x_cv[1], x_cv[2], var(int), diffABS(double)"
 * The first 2*nVars lines are always the min & max coordinate of the given computational domain
 * (This information is useful when you make a plot with GnuPlot), and its diffABS is always zero.
 */
void IkeWithPsALC_AD::writeCompatibilityEveryPtsOnTxtFile(const double threshold, const int nVars, double* rhsForAD, double* rhsForNormal) {
	// Get the min & max (x,y,z)
	double my_x_min[3] = { ABSURDLY_BIG_NUMBER,  ABSURDLY_BIG_NUMBER,  ABSURDLY_BIG_NUMBER};
	double my_x_max[3] = {-ABSURDLY_BIG_NUMBER, -ABSURDLY_BIG_NUMBER, -ABSURDLY_BIG_NUMBER};
	for(int icv=0; icv<ncv; ++icv) {
		my_x_min[0] = min(x_cv[icv][0], my_x_min[0]);
		my_x_min[1] = min(x_cv[icv][1], my_x_min[1]);
		my_x_min[2] = min(x_cv[icv][2], my_x_min[2]);
		my_x_max[0] = max(x_cv[icv][0], my_x_max[0]);
		my_x_max[1] = max(x_cv[icv][1], my_x_max[1]);
		my_x_max[2] = max(x_cv[icv][2], my_x_max[2]);
	}
	double tot_x_min[3], tot_x_max[3];
	MPI_Allreduce(my_x_min, tot_x_min, 3, MPI_DOUBLE, MPI_MIN, mpi_comm);
	MPI_Allreduce(my_x_max, tot_x_max, 3, MPI_DOUBLE, MPI_MAX, mpi_comm);

	// Write the difference on file
	if(mpi_rank != 0) { MPI_Status status; int dummy; MPI_Recv(&dummy, 1, MPI_INT, mpi_rank-1, 103, mpi_comm, &status); }

	// 1. Open the file
	string filename = "ERROR_IKE_JOE.csv";
	ofstream ofile;
	if(mpi_rank == 0)
		ofile.open(filename.c_str(), ios_base::out | ios_base::trunc);
	else
		ofile.open(filename.c_str(), ios_base::out | ios_base::app);

	// 2. Write data on the file
	if(mpi_rank==0) {
		printf("  Write difference on %s for thereshold>=%e\n", filename.c_str(), threshold);

		ofile<<"mpi_rank, icv, x_cv[0], x_cv[1], x_cv[2], var, diffABS"<<endl;
		for(int i=0; i<nVars; ++i) {
			ofile<<-1<<", "<<-1<<", "<<tot_x_min[0]<<", "<<tot_x_min[1]<<", "<<tot_x_min[2]<<", "<<i<<", 0.0"<<endl;
			ofile<<-1<<", "<<-1<<", "<<tot_x_max[0]<<", "<<tot_x_max[1]<<", "<<tot_x_max[2]<<", "<<i<<", 0.0"<<endl;
		}
		// dummyInt=step;				ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int)); // Note: if you simply use the "<<" operator, it will cause some problem due to bugs in a g++ library on linux
	}
	for(int icv=0; icv<ncv; ++icv) {
		int tempInt = nVars*icv;
		for(int i=0; i<nVars; ++i) {
			double tempDiff = rhsForAD[tempInt+i] - rhsForNormal[tempInt+i];
			if(fabs(tempDiff) >= threshold)
				ofile<<mpi_rank<<", "<<icv<<", "<<x_cv[icv][0]<<", "<<x_cv[icv][1]<<", "<<x_cv[icv][2]<<", "<<i<<", "<<fabs(tempDiff)<<endl;
		}
	}
	// 3. Close the file
	ofile.close();

	if(mpi_rank < mpi_size-1) { int dummy; MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 103, mpi_comm); }
	MPI_Barrier(mpi_comm);
}

/*
 * Method: writeQoIOnFile
 * ----------------------
 * By default, averaged density, maximum pressure, and subsonic portion will be stored
 */
void IkeWithPsALC_AD::writeQoIOnFile(const int step, string filename, const int NcontrolEqns, bool rewrite) {
	double avgRho = calcAvgRho();
	double maxPress = calcMaxPress();
	double subsonicVolPortion = calcSubsonicPortion();
	if(mpi_rank==0) {
		FILE *fp;
		if(rewrite)
			fp = fopen(filename.c_str(), "w");
		else
			fp = fopen(filename.c_str(), "a");

		fprintf(fp, "%d,\t", step);
		for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
			fprintf(fp, "%18.11e,\t",lambda[iEqn]);
		fprintf(fp, "%15.8e,\t%15.8e,\t%15.8e\n", avgRho, maxPress, subsonicVolPortion);

		fclose(fp);
	}

	MPI_Barrier(mpi_comm);
}

/*
 * Method: calcAvgRho
 * ------------------
 */
double IkeWithPsALC_AD::calcAvgRho() {
	double mySumOfRho = 0.0;
	double sumOfRho;
	double myVolume = 0.0;
	double Volume;

	for(int icv=0; icv<ncv; ++icv) {
//		if(x_cv[icv][0]>=xMin && x_cv[icv][0]<=xMax && x_cv[icv][1]>=yMin && x_cv[icv][1]<=yMax && x_cv[icv][2]>=zMin && x_cv[icv][2]<=zMax) {
			mySumOfRho += rho[icv]*cv_volume[icv];
			myVolume += cv_volume[icv];
//		}
	}

	MPI_Allreduce(&mySumOfRho, &sumOfRho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
	MPI_Allreduce(&myVolume, &Volume, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

	return sumOfRho/Volume;
}

/*
 * Method: calcMaxPress
 * --------------------
 */
double IkeWithPsALC_AD::calcMaxPress() {
	double* joePress = UgpWithCvCompFlow::press;
	assert(joePress != NULL);

	double myMaxPress = -ABSURDLY_BIG_NUMBER;
	double maxPress;

	for(int icv=0; icv<ncv; ++icv) {
//		if(x_cv[icv][0]>=xMin && x_cv[icv][0]<=xMax && x_cv[icv][1]>=yMin && x_cv[icv][1]<=yMax && x_cv[icv][2]>=zMin && x_cv[icv][2]<=zMax) {
			if(myMaxPress < joePress[icv])
				myMaxPress = joePress[icv];
//		}
	}

	MPI_Allreduce(&myMaxPress, &maxPress, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);

	return maxPress;
}

/*
 * Method: calcSubsonicPortion
 * ---------------------------
 */
double IkeWithPsALC_AD::calcSubsonicPortion() {
	double* joeSos = UgpWithCvCompFlow::sos;
	assert(joeSos != NULL);

	double mySubsonicVol = 0.0;
	double SubsonicVol;
	double myVolume = 0.0;
	double Volume;

	for(int icv=0; icv<ncv; ++icv) {
//		if(x_cv[icv][0]>=xMin && x_cv[icv][0]<=xMax && x_cv[icv][1]>=yMin && x_cv[icv][1]<=yMax && x_cv[icv][2]>=zMin && x_cv[icv][2]<=zMax) {
//			double velMag = vel[icv][0].value()*vel[icv][0].value() + vel[icv][1].value()*vel[icv][1].value() + vel[icv][2].value()*vel[icv][2].value();
			double velMag = (vecDotVec3d(rhou[icv], rhou[icv]))/rho[icv];
			if(velMag/joeSos[icv] < 1.0)
				SubsonicVol += cv_volume[icv];

			myVolume += cv_volume[icv];
//		}
	}

	MPI_Allreduce(&mySubsonicVol, &SubsonicVol, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
	MPI_Allreduce(&myVolume, &Volume, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

	return SubsonicVol/Volume;
}

/*
 * Method: convertSeparatedRhsTo1Drhs
 * ----------------------------------
 * Convert separated rhs (used in calcRhs()) to 1D rhs
 */
void IkeWithPsALC_AD::convertSeparatedRhsTo1Drhs(double* rhs1D, double *RHSrho, double (*RHSrhou)[3], double *RHSrhoE, double **RHSrhoScal, const int nScal) {
	for (int icv=0; icv<ncv; ++icv) {
		int dummy = icv*(5+nScal);
		rhs1D[dummy+0] = RHSrho[icv];
		for(int i=0; i<3; ++i)
			rhs1D[dummy+1+i] = RHSrhou[icv][i];
		rhs1D[dummy+4] = RHSrhoE[icv];
		for(int iScal=0; iScal<nScal; ++iScal)
			rhs1D[dummy+5+iScal] = RHSrhoScal[iScal][icv];
	}
}

/*
 * Method: convert1DrhsToSeparatedRhs
 * ----------------------------------
 * Update RHSrho, RHSrhou, RHSrhoE, (and scalars) from rhs1Darray.
 * These residual arrays will be used for tecplot output
 *
 * Note: If nScal==1, assume SA model is used
 *       If nScal==2, assume KW model is used
 */
void IkeWithPsALC_AD::convert1DrhsToSeparatedRhs(const double *rhs1Darray) {
	assert(RHSrho!=NULL && RHSrhou!=NULL && RHSrhoE!=NULL && rhs1Darray!=NULL);

	int nVars = 5+nScal;

	for(int icv=0; icv<ncv; ++icv) {
		int icvStart = icv*nVars;
		RHSrho[icv] = rhs1Darray[icvStart];
		for(int i=0; i<3; ++i)
			RHSrhou[icv][i] = rhs1Darray[icvStart+1+i];
		RHSrhoE[icv] = rhs1Darray[icvStart+4];

		if(nScal==1) {
			assert(RHSsa != NULL);

			RHSsa[icv] = rhs1Darray[icvStart+5];
		} else if (nScal==2) {
			assert(RHSkine != NULL && RHSomega != NULL);

			RHSkine[icv]  = rhs1Darray[icvStart+5];
			RHSomega[icv] = rhs1Darray[icvStart+6];
		}
	}
}

/*
 * Method: convertRHSrhoScalToSeparatedRhs
 * ---------------------------------------
 * Only for scalars!! : Update RHSsa or (RHSkine, RHSomega) from RHSrhoScal.
 * These residual arrays will be used for tecplot output
 *
 * Note: If nScal==1, assume SA model is used
 *       If nScal==2, assume KW model is used
 */
void IkeWithPsALC_AD::convertRHSrhoScalToSeparatedRhs(double **RHSrhoScal) {
	assert(RHSrhoScal!=NULL);

	for(int icv=0; icv<ncv; ++icv) {
		if(nScal==1) {
			assert(RHSsa != NULL);

			RHSsa[icv] = RHSrhoScal[0][icv];
		} else if (nScal==2) {
			assert(RHSkine != NULL && RHSomega != NULL);

			RHSkine[icv]  = RHSrhoScal[0][icv];
			RHSomega[icv] = RHSrhoScal[1][icv];
		}
	}
}

/*
 * Method: convert2DrhsTo1Drhs
 * ---------------------------
 * Convert 2D rhs (used in calcRhsCoupled()) to 1D rhs
 */
void IkeWithPsALC_AD::convert2DrhsTo1Drhs(double* rhs1D, double** rhs2D, const int nScal) {
	for(int icv=0; icv<ncv; ++icv) {
		int dummy = icv*(5+nScal);
		for (int i=0; i<5+nScal; i++)
			rhs1D[dummy+i] = rhs2D[icv][i];
	}
}

/*
 * Method: showResidue2
 * --------------------
 * Show residuals of the flow variables on the screen
 * Original code = showResidue() in JoeWithModels.cpp
 */
void IkeWithPsALC_AD::showResidue2(double *rhsResid, const bool showLabel) {
	int nScal = scalarTranspEqVector.size();

	// residual label
	if(mpi_rank==0 && showLabel) {
		printf("         rho          rhou-X       rhou-Y       rhou-Z       rhoE      ");
		for (int iScal = 0; iScal < nScal; iScal++)
			printf("%12s", scalarTranspEqVector[iScal].getName());
		cout << endl;
	}

	// residual value at each output step
	if (mpi_rank == 0) {
		printf("RESID: %12.4e %12.4e %12.4e %12.4e %12.4e", rhsResid[0], rhsResid[1], rhsResid[2], rhsResid[3], rhsResid[4]);
		for (int iScal = 0; iScal < nScal; iScal++)
			printf("%12.4e", rhsResid[5+iScal]);
		cout << endl;
	}
}

void IkeWithPsALC_AD::showResidue2(double *rhsResid, const int step, const int check_interval) {
	int nScal = scalarTranspEqVector.size();

	// residual label at every 10 output steps
	if (mpi_rank == 0)
		if( step%(check_interval*10) == 0 || step == 1 ) {
			printf("                          rho          rhou-X       rhou-Y       rhou-Z       rhoE      ");
			for (int iScal = 0; iScal < nScal; iScal++)
				printf("%12s", scalarTranspEqVector[iScal].getName());
			cout << endl;
		}

	// residual value at each output step
	if (mpi_rank == 0) {
		printf("          RESID: %6d %12.4e %12.4e %12.4e %12.4e %12.4e", step, rhsResid[0], rhsResid[1], rhsResid[2], rhsResid[3], rhsResid[4]);
		for (int iScal = 0; iScal < nScal; iScal++)
			printf("%12.4e", rhsResid[5+iScal]);
		cout << endl;
	}
}

/****************************
 * FOR DEBUGGING
 ****************************/
/*
 * Method: compatibilityCheck
 * --------------------------
 * The Psdueo-arclength continuation code uses two different types of RHS calculations: ADOL-C (for Jacobian calculation) and normal JOE (for backtracking)
 * This method checks if the two different calculation can give the same results or not.
 * Note: The RHS calculations should be same as those in calcJacobian1DAD()/calcJacobianAD() and backtrackWithJOE_calcRelaxAndRHS()
 *
 * Return: false if the two calculations are different
 */
bool IkeWithPsALC_AD::compatibilityCheck() {
	int nVars = 5+nScal;

	double* rhsForAD; // This will be allocated later since the size depends on "howToCalcJac"
	double* rhsForNormal = new double [ncv*nVars];
	for(int i=0; i<ncv*nVars; ++i)
		rhsForNormal[i] = 0.0;

	double *ResidAD     = new double [nVars];
	double *ResidNormal = new double [nVars];

	int countReducedCV_AD, countNegative_Normal;

	int debugLevel = getDebugLevel();

	if(mpi_rank==0) {
		cout<<endl<<"-----------------------------------------------------------------------"<<endl;
		cout<<">> Checking compatibility: RHS calculations by AD-IKE and by normal-JOE"<<endl;
		cout<<"-----------------------------------------------------------------------"<<endl;
	}

	/* ~~~~~~~~~~~~~~~~~~~~~~~~ */
	/* 1. Calculation by ADOL-C */
	/* ~~~~~~~~~~~~~~~~~~~~~~~~ */
	if(mpi_rank==0) cout<<"   1. RHS calculation by ADOL-C"<<endl;

	HOW_TO_CALC_JAC howToCalcJac = getHowToCalcJac();
	if(howToCalcJac == ROW_1D) {
		// UNIQUE-VARIABLE
		//----------------
		rhsForAD = new double[(5+nScal)*ncv];
		for(int i=0; i<(5+nScal)*ncv; ++i)
			rhsForAD[i] = 0.0;

		// THE FOLLOWING SHOULD BE THE SAME AS calcJacobian1DAD
		//-----------------------------------------------------
		assert(lambda_AD == NULL);

		int myCountReducedOrder = 0;
		int CountReducedOrder;

		// booleans for first-call
		bool firstCall = true;
		bool firstCallScalarTurb = true, firstCallScalarComb = true;

		// Wall-time
		double wtime0, wtime100, wtime1000, wtime10000;
		myWTimeJacCalc.wallTimeRHSinit  = 0.0;
		myWTimeJacCalc.wallTimeRHScalc  = 0.0;
		myWTimeJacCalc.wallTimeRHSfinal = 0.0;
		myWTimeJacCalc.wallTimeJac      = 0.0;
		myWTimeJacCalc.wallTimeCleanup  = 0.0;

		// Looping over the all the ICVs
		myBarrierMassSourceSum1D_AD   = 0.0; // Note: The barrier source term will be added to this variable in barrierSourceNS1D_AD()
		myBarrierEnergySourceSum1D_AD = 0.0; // Note: The barrier source term will be added to this variable in barrierSourceNS1D_AD()

		for(int icv=0; icv<ncv; ++icv) {
			if(debugLevel>0 && mpi_rank==0 && ncv>=100) {
				if(icv==0)
					wtime0 = MPI_Wtime();
				else if(icv==100) {
					wtime100 = MPI_Wtime();
					cout<<"     > Runtime for the first 100 RHS evaluations [sec] = "<<wtime100 - wtime0<<endl;
				}
				else if(icv==1000) {
					wtime1000 = MPI_Wtime();
					cout<<"     > Runtime for the first 1000 RHS evaluations [sec] = "<<wtime1000 - wtime0<<endl;
				}
				else if(icv==10000) {
					wtime10000 = MPI_Wtime();
					cout<<"     > Runtime for the first 10000 RHS evaluations [sec] = "<<wtime10000 - wtime0<<endl;
				}
			}

			//
			int tag = mpi_rank;
			calcJacobian1DAD_calcRhs(rhsForAD, myCountReducedOrder, myWTimeJacCalc, icv, tag, nScal, debugLevel, NcontrolEqns,
					firstCall, firstCallScalarTurb, firstCallScalarComb);

			//
			destroy_adoubles(NcontrolEqns);

			finalHookScalarRansTurbModel1D_AD();  // Reset internal pointers (e.g. kine, omega, etc) in a turb model
			finalHookScalarRansCombModel1D_AD();  // Reset internal pointers (e.g. ZMean, Zvar, CMean, etc) in a turb model

			nbocv2_eachIcv.clear();
			nbocv2ff_eachIcv.clear();
			fa2_eachIcv.clear();
		}

		// Show some statistics to the user
//		double totBarrierMassSourceSum1D_AD;
//		double totBarrierEnergySourceSum1D_AD;
//		MPI_Allreduce(&myBarrierMassSourceSum1D_AD,   &totBarrierMassSourceSum1D_AD,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
//		MPI_Allreduce(&myBarrierEnergySourceSum1D_AD, &totBarrierEnergySourceSum1D_AD, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
//
//		if(debugLevel>0 && mpi_rank==0) {
//			printf("              >> Barrier: Total mass source = %.4e, Total energy source = %.4e \n", totBarrierMassSourceSum1D_AD, totBarrierEnergySourceSum1D_AD);
//		}

		probeInCalcJacobian1DAD(jacMatrixSTL, rhsForAD, debugLevel, NcontrolEqns);
		    // Call probeInCalcJacobian1DAD() for the case that the user want to print out some statistics (on the screen)

		// free memory allocated for flow variables
		delete [] scalarTranspEqVector_AD;
		scalarTranspEqVector_AD = NULL;

		MPI_Allreduce(&myCountReducedOrder, &countReducedCV_AD, 1, MPI_INT, MPI_SUM, mpi_comm);
	} else if(howToCalcJac == ORDINARY_2D) {
		rhsForAD = new double[(5+nScal)*ncv_gg];

		// +++++++++++++++++++++++++++++++++
		// Allocate memory for Flow residual
		// +++++++++++++++++++++++++++++++++
		REALA *rhs_rho_AD       = new REALA[ncv_gg];
		REALA (*rhs_rhou_AD)[3] = new REALA[ncv_gg][3];
		REALA *rhs_rhoE_AD      = new REALA[ncv_gg];

		REALAS **rhs_rhoScal_AD = NULL;
		if (nScal > 0)
			getMem2D(&rhs_rhoScal_AD, 0, nScal-1, 0, ncv_gg-1, "rhs_rhoScal_AD");

		// +++++++++++++++++++++++++++++++++
		// Start the calculation of Residual wrt flow
		// +++++++++++++++++++++++++++++++++
		int tag = mpi_rank;
		trace_on(tag, 1); // Note: trace_on(tag, keep) -- If keep==1, the numerical values of all active variables
                          //       are recorded in a buffered temporary file before they will be overwritten.
                          //       (preparing the scene for an immediately following reverse mode differentiation)

		// allocate memory
		initialize_adoubles();
		assert(lambda_AD != NULL);

		// Independent variables
		for (int icv = 0; icv < ncv_gg; icv++) {
			rho_AD[icv]      <<= rho[icv] ;
			rhou_AD[icv][0]  <<= rhou[icv][0] ;
			rhou_AD[icv][1]  <<= rhou[icv][1] ;
			rhou_AD[icv][2]  <<= rhou[icv][2] ;
			rhoE_AD[icv]     <<= rhoE[icv] ;
			for(int iScal=0; iScal<nScal; iScal++)
				scalarTranspEqVector_AD[iScal].phi[icv] <<= scalarTranspEqVector[iScal].phi[icv];
		}
		if(NcontrolEqns > 0)
			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
				lambda_AD[iEqn] <<= lambda[iEqn];

		initialHook_AD(); // Note: The flow field should be initialized with initialHook()
		                  // To DO: in the original code initialHook_AD() isn't called here ( see JoeWithModels_AD::calcResidualDerivative() )
		initialHookScalarRansTurbModel_AD();
		initialHookScalarRansCombModel_AD(0);

		// Check possible errors: scalars = NULL
		checkScalarsMemAD();

		// Calculate Rhs
		double (*A)[5][5] = NULL, ***AScal=NULL;
		int flagImplicit = false;
		countReducedCV_AD = calcResidual_AD(rhs_rho_AD, rhs_rhou_AD, rhs_rhoE_AD, rhs_rhoScal_AD, rho_AD, rhou_AD, rhoE_AD, A, AScal, flagImplicit);
	          // Note: calcaResidual_AD() successively calls the following methods:
	          //       1. calcStateVariables_AD(),   2. calcMaterialProperties_AD(),
	          //       3. setBC_AD(),
	          //       if(mu_ref>0.0 || sndOrder==true): 4. calcCv2Grad_AD(),
	          //       if(mu_ref>0.0):                   5. calcRansTurbViscMuet_AD(),
	          //       6. calcRhs_AD()
		barrierSourceNS_AD(rhs_rho_AD, rhs_rhou_AD, rhs_rhoE_AD);
		barrierSourceTurbScalars_AD(rhs_rhoScal_AD, nScal, iterNewton, residNormTotOld);

		//Dependent variables
		for (int icv = 0; icv < ncv_gg; icv++) {
			int icv_i = nVars*icv;
			rhs_rho_AD[icv]      >>= rhsForAD[icv_i] ;
			rhs_rhou_AD[icv][0]  >>= rhsForAD[icv_i+1] ;
			rhs_rhou_AD[icv][1]  >>= rhsForAD[icv_i+2] ;
			rhs_rhou_AD[icv][2]  >>= rhsForAD[icv_i+3] ;
			rhs_rhoE_AD[icv]     >>= rhsForAD[icv_i+4] ;
			for(int iScal=0; iScal<nScal; iScal++) {
				rhs_rhoScal_AD[iScal][icv] >>= rhsForAD[icv_i+5+iScal];
			}
		}
		trace_off(); // Note: trace_off(file) -- If the argument "file" is omitted, it defaults to 0,
                     //                          so that the tape array is written onto an external file
                     //                          only if the length of any of the buffers exceeds BUFSIZE

		checkNanOrInfRhs_AD(rhs_rho_AD, rhs_rhou_AD, rhs_rhoE_AD, rhs_rhoScal_AD);

		if(debugLevel>0 && mpi_rank==0) {
			cout<<"=================================="<<endl;
			cout<<"Tape stats for Residual evaluation"<<endl;
			cout<<"=================================="<<endl;
			print_tapestats(tag);
		}
		int memoryDeficit = countMemoryDeficit_fromTapeStats(tag);
		if(memoryDeficit > 0)
			cout<<"WARNING! Too much memory is requires for AD at mpi_rank=="<<mpi_rank<<": memory deficit = "<<memoryDeficit<<endl;

		// +++++++++++++++++++++++++++++++++
		// Free memory
		// +++++++++++++++++++++++++++++++++
		removeTape(tag, ADOLC_REMOVE_COMPLETELY);  // Note: "int removeTape(short tapeID, short type)" in tape_handling.cpp
		                                           //       If type==ADOLC_REMOVE_COMPLETELY, remove the tape instead of free it.
		                                           //       Returns 0 if the removing process is completed okay.
		destroy_adoubles();

		freeMem2D(rhs_rhoScal_AD, 0, nScal-1, 0, ncv_gg-1);
		delete [] rhs_rhoE_AD;
		delete [] rhs_rhou_AD;
		delete [] rhs_rho_AD;
	} else
		throw(PSALC_ERROR_CODE);

	if(debugLevel>0) {
		calcResidualsFrom1Drhs(ResidAD, rhsForAD, 1);
		showResidue2(ResidAD, true);
	}
	if(countReducedCV_AD>0 && mpi_rank==0)
		cout<<"     countReducedCV_AD(="<<countReducedCV_AD<<") > 0, but proceeds anyway"<<endl;
	MPI_Barrier(mpi_comm);

	/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
	/* 2. Calculation by normal JOE */
	/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
	// To Do: calcRhsCoupled() is not working properly for some reasons
	if(mpi_rank==0) cout<<endl<<"   2. RHS calculation by normal JOE (calcRhs)"<<endl;

	bool useBarrier = true;
	countNegative_Normal = calcRhsWithBarrier(rhsForNormal, useBarrier);
	temporalHook(); // Call temporalHook() for the case that the user want to print out some statistics (on the screen)

	if(debugLevel>0) {
		calcResidualsFrom1Drhs(ResidNormal, rhsForNormal, 1);
		showResidue2(ResidNormal, true);
	}
	if(countNegative_Normal>0 && mpi_rank==0)
		cout<<"     countNegative_Normal(="<<countNegative_Normal<<") > 0, but proceeds anyway"<<endl;

	/* ~~~~~~~~~~~~~~~~~~~~~~~~ */
	/* 3. Compare the two calcs */
	/* ~~~~~~~~~~~~~~~~~~~~~~~~ */
//	double myDiffSq = 0.0, DiffSq;
	double myDiffAbs =  0.0, DiffAbs;
	double myDiffMax = -1.0, DiffMax;
	for(int i=0; i<ncv*nVars; ++i) {
//		myDiffSq += pow(rhsForAD[i]-rhsForNormal[i], 2.0);
		myDiffAbs += fabs(rhsForAD[i]-rhsForNormal[i]);
		myDiffMax = max(myDiffMax, fabs(rhsForAD[i]-rhsForNormal[i]));
	}
//	MPI_Allreduce(&myDiffSq,  &DiffSq,  1, MPI_DOUBLE, MPI_SUM, mpi_comm);
	MPI_Allreduce(&myDiffAbs, &DiffAbs, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
	MPI_Allreduce(&myDiffMax, &DiffMax, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);

	if(mpi_rank==0) {
		cout<<endl;
		printf("   *****************************************************************\n");
		printf("   Comparison Result: \n");
		printf("     MAXIMUM DIFFERENCE = %.4e \n", DiffMax);
		printf("     SUM OF DIFFERENCE  = %.4e \n", DiffAbs);
		printf("     RESIDUAL SUMMARY (1-norm): \n");
		printf("        VAR  RESID_AD               RESID_NORMAL          RESID_DIFF \n");
		for(int i=0; i<nVars; ++i)
			printf("        %3d  %20.13e  %20.13e  %11.4e \n", i, ResidAD[i], ResidNormal[i], ResidAD[i]-ResidNormal[i]);
		printf("   *****************************************************************\n\n");
	}

	bool negligibleError = true;
	double tolCompatibility = getDoubleParam("TOL_COMPATIBILITY", "1.0e4");
//	if(isNaN(DiffSq) || (sqrt(DiffSq) > tolCompatibility*sqrt((double)cvora[mpi_size])*MACHINE_EPS)) {
	if(isNaN(DiffAbs) || (DiffAbs > tolCompatibility*((double) cvora[mpi_size])*MACHINE_EPS)) {
		if(mpi_rank==0) {
//			printf("  WARNING! RHS vectors are not the same: error 2-norm = %.3e (> THRESHOLD=%.3e)\n", sqrt(DiffSq), tolCompatibility*sqrt((double)cvora[mpi_size])*MACHINE_EPS);
			printf("  WARNING! RHS vectors are not the same: error 1-norm = %.3e (> THRESHOLD=%.3e)\n", DiffAbs, tolCompatibility*((double) cvora[mpi_size])*MACHINE_EPS);
			printf("           * POSSIBLE ERROR SOURCES: Hook functions, CPU boarder near a domain boundary, etc. \n");
		}
		if(debugLevel>1) { // If debugLevel is high, show the CV information
			if(mpi_rank != 0) { int dummy_Recv; 	MPI_Status status; 	MPI_Recv(&dummy_Recv, 1, MPI_INT, mpi_rank-1, 313, mpi_comm, &status); }
			if(mpi_rank == 0) {
				cout<<"  ERROR DETAILS -- show the CVs whose errors are not negligible:"<<endl;
				cout<<"  mpi_rank  icv  x          y          z          ";
				for(int i=0; i<nVars; ++i)
					printf("var%d        ",i);
				cout<<endl;
			}
			for(int icv=0; icv<ncv; ++icv) {
				double sumAbsErrors = 0.0;
				int tempInt = nVars*icv;
				for(int i=0; i<nVars-1; ++i)
					sumAbsErrors += fabs(rhsForAD[tempInt+i]-rhsForNormal[tempInt+i]);
				if(isNaN(sumAbsErrors) || (sumAbsErrors > 1000.0*(double)nVars*MACHINE_EPS)) {
					printf("  %3d   %5d   %10.3e %10.3e %10.3e", mpi_rank, icv, x_cv[icv][0], x_cv[icv][1], x_cv[icv][2]);
					for(int i=0; i<nVars; ++i) 
						printf(" %11.4e", rhsForAD[tempInt+i]-rhsForNormal[tempInt+i]);
					cout<<endl;
				}
			}
			if(mpi_rank < mpi_size-1) { int dummy_Send = 1; 	MPI_Send(&dummy_Send, 1, MPI_INT, mpi_rank+1, 313, mpi_comm); }

			MPI_Barrier(mpi_comm);
		}
		negligibleError = false;
	} else {
		if(mpi_rank==0) {
//			printf("  Error 2-norm = %.3e -- OK \n", sqrt(DiffSq));
			printf("  Error 1-norm = %.3e -- OK (<Threshold=%.3e)\n", DiffAbs, tolCompatibility*((double) cvora[mpi_size])*MACHINE_EPS);
			cout<<"-----------------------------------------------------------------------"<<endl;
		}
		negligibleError = true;
	}

	if(countReducedCV_AD != countNegative_Normal) {
		if(mpi_rank==0) printf("  WARNING! countReducedCV_AD(=%d) and countNegative_Normal(=%d) are not the same \n", countReducedCV_AD, countNegative_Normal);
		negligibleError = false;
	}

    if(!negligibleError) {
    	if(debugLevel>1) {
    		// Write the difference on a txt file only for the CVs where the difference for each variable is greater than thresholdWriteTxt
    		double thresholdWriteTxt = 1.0e-7;
    		writeCompatibilityEveryPtsOnTxtFile(thresholdWriteTxt, nVars, rhsForAD, rhsForNormal);
    	}

    	// Write the difference on a output (Tecplot) file
        assert(RHSrho != NULL && RHSrhou != NULL && RHSrhoE != NULL);
        for(int icv=0; icv<ncv; ++icv) {
			int tempInt = nVars*icv;

            RHSrho[icv]  = fabs(rhsForAD[tempInt]   - rhsForNormal[tempInt]);
            for(int i=0; i<3; ++i)
                RHSrhou[icv][i] = fabs(rhsForAD[tempInt+1+i]  - rhsForNormal[tempInt+1+i]);

            RHSrhoE[icv] = fabs(rhsForAD[tempInt+4] - rhsForNormal[tempInt+4]);

            if(nScal==1) {
                RHSsa[icv]  = fabs(rhsForAD[tempInt+5] - rhsForNormal[tempInt+5]);
            } else if (nScal==2) {
                RHSkine[icv]  = fabs(rhsForAD[tempInt+5] - rhsForNormal[tempInt+5]);
                RHSomega[icv] = fabs(rhsForAD[tempInt+6] - rhsForNormal[tempInt+6]);
            } else if (nScal==3) {
                assert(RHSZMean != NULL && RHSZVar != NULL && RHSCMean != NULL);
                RHSZMean[icv] = fabs(rhsForAD[tempInt+5] - rhsForNormal[tempInt+5]);
                RHSZVar[icv]  = fabs(rhsForAD[tempInt+6] - rhsForNormal[tempInt+6]);
                RHSCMean[icv] = fabs(rhsForAD[tempInt+7] - rhsForNormal[tempInt+7]);
            } else if (nScal==5) {
                RHSkine[icv]  = fabs(rhsForAD[tempInt+5] - rhsForNormal[tempInt+5]);
                RHSomega[icv] = fabs(rhsForAD[tempInt+6] - rhsForNormal[tempInt+6]);
                RHSZMean[icv] = fabs(rhsForAD[tempInt+7] - rhsForNormal[tempInt+7]);
                RHSZVar[icv]  = fabs(rhsForAD[tempInt+8] - rhsForNormal[tempInt+8]);
                RHSCMean[icv] = fabs(rhsForAD[tempInt+9] - rhsForNormal[tempInt+9]);
            }
        }

        if(mpi_rank==0) {
            printf("  SAVING the N-S difference between AD and JOE in RHSRHO, RHSRHOU, RHSRHOE\n");
            if (nScal==1) 
                printf("  SAVING the SA difference between AD and JOE in RHSSA\n");
            else if (nScal==2) 
                printf("  SAVING the KOm difference between AD and JOE in RHSKINE, RHSOMEGA\n");
            else if (nScal==3) 
                printf("  SAVING the FPVA difference between AD and JOE in RHSZMEAN, RHSZVAR, RHSCMEAN\n");
            else if (nScal==5) 
                printf("  SAVING the KOm and FPVA difference between AD and JOE in RHSKINE, RHSOMEGA, RHSZMEAN, RHSZVAR, RHSCMEAN\n");
        }

        MPI_Barrier(mpi_comm);
        
        writeData(0);
        MPI_Barrier(mpi_comm);
    }

	/* ~~~~~~~~~~~~ */
	/* Clear memory */
	/* ~~~~~~~~~~~~ */
	delete [] rhsForAD;
	delete [] rhsForNormal;

	delete [] ResidAD;
	delete [] ResidNormal;

	if(mpi_rank==0)
		cout<<endl;

	MPI_Barrier(mpi_comm);

	return negligibleError;
}

/*
 * Method: checkNanOrInfRhs1D_AD
 * -----------------------------
 *
 */
bool IkeWithPsALC_AD::checkNanOrInfRhs1D_AD(const int icv, adouble &rhs_rho_AD, adouble *rhs_rhou_AD, adouble &rhs_rhoE_AD, adouble *rhs_rhoScal_AD) {
	bool showComment;
	bool nonDetected = false;
//	bool infDetected = false;

	if(debugLevel>0) {
		// Check for NaN
		showComment = true;
		nonDetected = false;
		if(isNaN(rhs_rho_AD.value()) || isInfinite(rhs_rho_AD.value())) {
			if(showComment) {
				printf("!CAUTION! NaN or Inf occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
				showComment = false;
			}
			cout<<" RHS_RHO="<<rhs_rho_AD.value()<<" ";
			nonDetected = true;
		}
		if(isNaN(rhs_rhou_AD[0].value()) || isInfinite(rhs_rhou_AD[0].value())) {
			if(showComment) {
				printf("!CAUTION! NaN or Inf occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
				showComment = false;
			}
			cout<<" RHS_RHOU[0]="<<rhs_rhou_AD[0].value()<<" ";
			nonDetected = true;
		}
		if(isNaN(rhs_rhou_AD[1].value())) {
			if(showComment) {
				printf("!CAUTION! NaN or Inf occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
				showComment = false;
			}
			cout<<" RHS_RHOU[1]="<<rhs_rhou_AD[1].value()<<" ";
			nonDetected = true;
		}
		if(isNaN(rhs_rhou_AD[2].value())) {
			if(showComment) {
				printf("!CAUTION! NaN or Inf occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
				showComment = false;
			}
			cout<<" RHS_RHOU[2]="<<rhs_rhou_AD[2].value()<<" ";
			nonDetected = true;
		}
		if(isNaN(rhs_rhoE_AD.value())) {
			if(showComment) {
				printf("!CAUTION! NaN or Inf occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
				showComment = false;
			}
			cout<<" RHS_RHOE="<<rhs_rhoE_AD.value()<<" ";
			nonDetected = true;
		}
		for(int iScal=0; iScal<nScal; ++iScal) {
			if(isNaN(rhs_rhoScal_AD[iScal].value())) {
				if(showComment) {
					printf("!CAUTION! NaN or Inf occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
					showComment = false;
				}
				cout<<" RHS_RHOscal["<<iScal<<"]="<<rhs_rhoScal_AD[iScal].value()<<" ";
				nonDetected = true;
			}
		}
		if(nonDetected)
			cout<<endl;

//		// Check for Infinite
//		showComment = true;
//		infDetected = false;
//		if(isInfinite(rhs_rho_AD.value())) {
//			if(showComment) {
//				printf("!CAUTION! Inf occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
//				showComment = false;
//			}
//			cout<<" RHS_RHO ";
//			infDetected = true;
//		}
//		if(isNaN(rhs_rhou_AD[0].value())) {
//			if(showComment) {
//				printf("!CAUTION! Inf occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
//				showComment = false;
//			}
//			cout<<" RHS_RHOU[0] ";
//			infDetected = true;
//		}
//		if(isInfinite(rhs_rhou_AD[1].value())) {
//			if(showComment) {
//				printf("!CAUTION! Inf occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
//				showComment = false;
//			}
//			cout<<" RHS_RHOU[1] ";
//			infDetected = true;
//		}
//		if(isInfinite(rhs_rhou_AD[2].value())) {
//			if(showComment) {
//				printf("!CAUTION! Inf occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
//				showComment = false;
//			}
//			cout<<" RHS_RHOU[2] ";
//			infDetected = true;
//		}
//		if(isInfinite(rhs_rhoE_AD.value())) {
//			if(showComment) {
//				printf("!CAUTION! Inf occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
//				showComment = false;
//			}
//			cout<<" RHS_RHOE ";
//			infDetected = true;
//		}
//		for(int iScal=0; iScal<nScal; ++iScal) {
//			if(isInfinite(rhs_rhoScal_AD[iScal].value())) {
//				if(showComment) {
//					printf("!CAUTION! Inf occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
//					showComment = false;
//				}
//				cout<<" RHS_RHOscal["<<iScal<<"] ";
//				infDetected = true;
//			}
//		}
//		if(infDetected)
//			cout<<endl;
	}

	if(nonDetected) // || infDetected)
		return true;
	return false;
}

/*
 * Method: checkNanOrInfRhs_AD
 * -----------------------------
 * To DO: make this as checkNanOrInfRhs1D_AD()
 */
bool IkeWithPsALC_AD::checkNanOrInfRhs_AD(REALA *rhs_rho_AD, REALA (*rhs_rhou_AD)[3], REALA *rhs_rhoE_AD, REALA **rhs_rhoScal_AD) {
	int myCountNaN = 0, countNaN;
	int myCountInf = 0, countInf;
	bool foundError = false;
	if(debugLevel>0) {
		// Check for NaN
		for(int icv=0; icv<ncv; ++icv) {
			bool showCommentForIcv = true;
			bool nonDetectedForIcv = false;
			if(isNaN(rhs_rho_AD[icv].value())) {
				if(showCommentForIcv) {
					printf("!CAUTION! NaN occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
					showCommentForIcv = false;
				}
				cout<<" RHS_RHO ";
				nonDetectedForIcv = true;
				++myCountNaN;
			}
			if(isNaN(rhs_rhou_AD[icv][0].value())) {
				if(showCommentForIcv) {
					printf("!CAUTION! NaN occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
					showCommentForIcv = false;
				}
				cout<<" RHS_RHOU[0] ";
				nonDetectedForIcv = true;
				++myCountNaN;
			}
			if(isNaN(rhs_rhou_AD[icv][1].value())) {
				if(showCommentForIcv) {
					printf("!CAUTION! NaN occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
					showCommentForIcv = false;
				}
				cout<<" RHS_RHOU[1] ";
				nonDetectedForIcv = true;
				++myCountNaN;
			}
			if(isNaN(rhs_rhou_AD[icv][2].value())) {
				if(showCommentForIcv) {
					printf("!CAUTION! NaN occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
					showCommentForIcv = false;
				}
				cout<<" RHS_RHOU[2] ";
				nonDetectedForIcv = true;
				++myCountNaN;
			}
			if(isNaN(rhs_rhoE_AD[icv].value())) {
				if(showCommentForIcv) {
					printf("!CAUTION! NaN occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
					showCommentForIcv = false;
				}
				cout<<" RHS_RHOE ";
				nonDetectedForIcv = true;
				++myCountNaN;
			}
			for(int iScal=0; iScal<nScal; ++iScal) {
				if(isNaN(rhs_rhoScal_AD[iScal][icv].value())) {
					if(showCommentForIcv) {
						printf("!CAUTION! NaN occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
						showCommentForIcv = false;
					}
					cout<<" RHS_RHOscal["<<iScal<<"] ";
					nonDetectedForIcv = true;
					++myCountNaN;
				}
			}
			if(nonDetectedForIcv)
				cout<<endl;
		}

		// Check for Infinite
		for(int icv=0; icv<ncv; ++icv) {
			bool showCommentForIcv = true;
			bool infDetectedForIcv = false;
			if(isInfinite(rhs_rho_AD[icv].value())) {
				if(showCommentForIcv) {
					printf("!CAUTION! Inf occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
					showCommentForIcv = false;
				}
				cout<<" RHS_RHO ";
				infDetectedForIcv = true;
				++myCountInf;
			}
			if(isNaN(rhs_rhou_AD[icv][0].value())) {
				if(showCommentForIcv) {
					printf("!CAUTION! Inf occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
					showCommentForIcv = false;
				}
				cout<<" RHS_RHOU[0] ";
				infDetectedForIcv = true;
				++myCountInf;
			}
			if(isInfinite(rhs_rhou_AD[icv][1].value())) {
				if(showCommentForIcv) {
					printf("!CAUTION! Inf occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
					showCommentForIcv = false;
				}
				cout<<" RHS_RHOU[1] ";
				infDetectedForIcv = true;
				++myCountInf;
			}
			if(isInfinite(rhs_rhou_AD[icv][2].value())) {
				if(showCommentForIcv) {
					printf("!CAUTION! Inf occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
					showCommentForIcv = false;
				}
				cout<<" RHS_RHOU[2] ";
				infDetectedForIcv = true;
				++myCountInf;
			}
			if(isInfinite(rhs_rhoE_AD[icv].value())) {
				if(showCommentForIcv) {
					printf("!CAUTION! Inf occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
					showCommentForIcv = false;
				}
				cout<<" RHS_RHOE ";
				infDetectedForIcv = true;
				++myCountInf;
			}
			for(int iScal=0; iScal<nScal; ++iScal) {
				if(isInfinite(rhs_rhoScal_AD[iScal][icv].value())) {
					if(showCommentForIcv) {
						printf("!CAUTION! Inf occurs at icv=%d(%.3f,%.3f,%.3f): ",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2]);
						showCommentForIcv = false;
					}
					cout<<" RHS_RHOscal["<<iScal<<"] ";
					infDetectedForIcv = true;
					++myCountInf;
				}
			}
			if(infDetectedForIcv)
				cout<<endl;
		}

		MPI_Allreduce(&myCountNaN, &countNaN, 1, MPI_INT, MPI_SUM, mpi_comm);
		MPI_Allreduce(&myCountInf, &countInf, 1, MPI_INT, MPI_SUM, mpi_comm);

		if(countNaN>0 || countInf>0)
			foundError = true;
	}
	return foundError;
}

/*
 * Method: showBoundaryInfoOnScreen
 * --------------------------------
 *
 */
void IkeWithPsALC_AD::showBoundaryInfoOnScreen(string& boundaryName) {
	if(mpi_rank == 0) {
		cout<<"IkeWithPsALC_AD::showBoundaryInfo(): boundary information for "<<boundaryName<<endl;
		cout<<"  mpi_rank   ifa  x          y          z"<<endl;
	}

	if(mpi_rank != 0) {
		int dummy_Recv;
		MPI_Status status;
		MPI_Recv(&dummy_Recv, 1, MPI_INT, mpi_rank-1, 1201, mpi_comm, &status);
	}

	{
		for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
			if (zone->getKind() == FA_ZONE_BOUNDARY) {
				if( zone->getName() == boundaryName ) {
					for (int index = 0; index < zone->faVec.size(); ++index) {
						int ifa = zone->faVec[index];
						int icv1 = cvofa[ifa][1];
						printf("  %7d  %5d  %9.3e  %9.3e  %9.3e\n", mpi_rank, ifa, x_fa[ifa][0], x_fa[ifa][1], x_fa[ifa][2]);
					}
				}
			}
		}
	}

	if(mpi_rank < mpi_size-1) {
		int dummy_Send = 1;
		MPI_Send(&dummy_Send, 1, MPI_INT, mpi_rank+1, 1201, mpi_comm);
	}

	MPI_Barrier(mpi_comm);
}
void IkeWithPsALC_AD::showBoundaryInfoOnScreen(char boundaryName[]) {
	string bName(boundaryName);
	showBoundaryInfoOnScreen(bName);
}

/*
 * Method: showFlowfieldOnScreen
 * -----------------------------
 * Display data on the screen
 * You must provide gammaRef and RoMRef for the case that gamma and RoM were not established yet.
 */
void IkeWithPsALC_AD::showFlowfieldOnScreen(double gammaRef, double RoMRef, bool showGhosts, bool showFakes) {
	int myUseRoM = 1, countUseRoM;
	if(RoMRef!=UgpWithCvCompFlow::RoM[0])
		myUseRoM = 0;
	MPI_Allreduce(&myUseRoM, &countUseRoM, 1, MPI_INT, MPI_SUM, mpi_comm);

	if(mpi_rank != 0) {
		int dummy_Recv;
		MPI_Status status;
		MPI_Recv(&dummy_Recv, 1, MPI_INT, mpi_rank-1, 1111, mpi_comm, &status);
	}


	if(mpi_rank == 0) {
		cout<<endl<<"IkeWithPsALC_AD::showFlowfieldOnScreen(): CV data information "<<endl;
		cout<<"mpi_rank  icv  x          y          z          rho          x-vel        y-vel        z-vel        rhoE         press        temp"<<endl;

		if(countUseRoM!=mpi_size)
			printf("  WARNING! RoM[icv=0]=%.3e is not equal to RoMRef=P_REF/(RHO_REF*T_REF)=%.3e:  use RoMRef and gammaRef here\n", UgpWithCvCompFlow::RoM[0], RoMRef);
	}
	{
		double *kineArray = NULL;
		int kine_Index = getScalarTransportIndex("kine");
		if(kine_Index>=0)
			kineArray = scalarTranspEqVector[kine_Index].phi;
		for(int icv=0; icv<ncv; ++icv) {
			double kinecv = 0.0;
			if(kine_Index>=0)
				kinecv = kineArray[icv];
			double press, temp;
			if(countUseRoM!=mpi_size) {
				press = calcPress(gammaRef, rhoE[icv], rhou[icv], rho[icv], kinecv);
				temp = calcTemp(press, rho[icv], RoMRef);
			} else {
				press = calcPress(UgpWithCvCompFlow::gamma[icv], rhoE[icv], rhou[icv], rho[icv], kinecv);
				temp = calcTemp(press, rho[icv], UgpWithCvCompFlow::RoM[icv]);
			}

			printf("%2d    %5d   %10.3e %10.3e %10.3e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e",
					mpi_rank, icv, x_cv[icv][0], x_cv[icv][1], x_cv[icv][2],
					rho[icv], rhou[icv][0]/rho[icv], rhou[icv][1]/rho[icv], rhou[icv][2]/rho[icv], rhoE[icv], press, temp);
			for(int iScal=0; iScal<nScal; ++iScal)
				printf(" %12.5e", scalarTranspEqVector[iScal].phi[icv]);
			printf("\n");
		}

		if(showGhosts) {
			for(int icv=ncv; icv<ncv_gg; ++icv) {
				if(ncv<=icv && icv<ncv_g)
					cout<<"G1     ";
				else if(ncv_g<=icv && icv<ncv_gg)
					cout<<"G2     ";

				printf("%5d   %10.3e %10.3e %10.3e %12.5e %12.5e %12.5e %12.5e %12.5e                          ",
						icv, x_cv[icv][0], x_cv[icv][1], x_cv[icv][2],
						rho[icv], rhou[icv][0]/rho[icv], rhou[icv][1]/rho[icv], rhou[icv][2]/rho[icv], rhoE[icv]);
				for(int iScal=0; iScal<nScal; ++iScal)
					printf(" %12.5e", scalarTranspEqVector[iScal].phi[icv]);
				printf("\n");
			}
		}

		if(showFakes) {
			for(int icv=ncv_gg; icv<ncv_ggff; ++icv) {
				if(ncv_gg<=icv && icv<ncv_ggf)
					cout<<"F1     ";
				else if(ncv_ggf<=icv && icv<ncv_ggff)
					cout<<"F2     ";

				printf("%5d   %10.3e %10.3e %10.3e %12.5e %12.5e %12.5e %12.5e %12.5e                          ",
						icv, x_cv[icv][0], x_cv[icv][1], x_cv[icv][2],
						rho[icv], rhou[icv][0]/rho[icv], rhou[icv][1]/rho[icv], rhou[icv][2]/rho[icv], rhoE[icv]);
				for(int iScal=0; iScal<nScal; ++iScal)
					printf(" %12.5e", scalarTranspEqVector[iScal].phi[icv]);
				printf("\n");
			}
		}
	}
	if(mpi_rank==mpi_size-1){
		cout<<endl;
	}
	if(mpi_rank < mpi_size-1) {
		int dummy_Send = 1;
		MPI_Send(&dummy_Send, 1, MPI_INT, mpi_rank+1, 1111, mpi_comm);
	}
	MPI_Barrier(mpi_comm);
}

/*
 * Method: show1DArrayOnScreen
 * -----------------------------
 * Display data on the screen with the x_cv information
 * Argument: Q = 1D array (size = ncv*nVars)
 *
 * Note: This is a MPI call. Thus, all the cpu cores print out its data on the screen.
 */
void IkeWithPsALC_AD::show1DArrayOnScreen(const double* Q, const int nVars) {
	if(mpi_rank != 0) {
		int dummy_Recv;
		MPI_Status status;
		MPI_Recv(&dummy_Recv, 1, MPI_INT, mpi_rank-1, 1111, mpi_comm, &status);
	}

	if(mpi_rank == 0) {
		cout<<endl<<"IkeWithPsALC_AD::show1DArrayOnScreen() "<<endl;
		cout<<"mpi_rank  icv  x          y          z          ";
		for(int i=0; i<nVars; ++i)
			printf("var%d         ", i);
		cout<<endl;
	}
	{
		for(int icv=0; icv<ncv; ++icv) {
			printf("%2d    %5d   %10.3e %10.3e %10.3e", mpi_rank, icv, x_cv[icv][0], x_cv[icv][1], x_cv[icv][2]);
			int tempInt = icv*nVars;
			for(int i=0; i<nVars; ++i)
				printf(" %12.5e", Q[tempInt+i]);
			printf("\n");
		}
	}
	if(mpi_rank==mpi_size-1){
		cout<<endl;
	}
	if(mpi_rank < mpi_size-1) {
		int dummy_Send = 1;
		MPI_Send(&dummy_Send, 1, MPI_INT, mpi_rank+1, 1111, mpi_comm);
	}

	MPI_Barrier(mpi_comm);
	for(int i=0; i<10000; ++i) {}
	MPI_Barrier(mpi_comm);
}

/*
 * Method: show1DArrayOnScreen
 * ---------------------------
 * Display data on the screen with the x_cv information
 * (For each entry in Q, if icv<ncv, show the x_cv information. Otherwise, just show the value of the entry)
 * Argument: Q = 1D array (size = arrayLength)
 *
 * Note: This is a MPI call. Thus, all the cpu cores print out its data on the screen.
 */
void IkeWithPsALC_AD::show1DArrayOnScreen(const double* Q, const int nVars, int arrayLength, const bool showIcv/*=true*/) {
	int totArrayLength = 0;
	MPI_Allreduce(&arrayLength, &totArrayLength, 1, MPI_INT, MPI_SUM, mpi_comm);
	MPI_Barrier(mpi_comm);

	if(mpi_rank != 0) {
		int dummy_Recv;
		MPI_Status status;
		MPI_Recv(&dummy_Recv, 1, MPI_INT, mpi_rank-1, 1111, mpi_comm, &status);
	}

	if(mpi_rank == 0) {
		cout<<"IkeWithPsALC_AD::show1DArrayOnScreen() "<<endl;
		if(showIcv)
			cout<<"mpi_rank  icv  x          y          z          value"<<endl;
	}
	{
		for(int i=0; i<arrayLength; ++i) {
			int icv = int(i/nVars);
			assert(icv<ncv+NcontrolEqns && icv>=0);

			if(showIcv) {
				if(i%nVars == 0) {
					if(icv<ncv)
						printf("%2d    %5d   %10.3e %10.3e %10.3e %12.5e\n", mpi_rank, icv, x_cv[icv][0], x_cv[icv][1], x_cv[icv][2], Q[i]);
					else
						printf("%2d    DUNNO                               %12.5e\n", mpi_rank, Q[i]);
				} else
					printf("                                                   %12.5e\n", Q[i]);
			} else {
				printf("%12.5e\n", Q[i]);
			}
		}
	}

	if(mpi_rank==mpi_size-1){
		cout<<">> TOTAL ARRAY LENGTH="<<totArrayLength<<endl<<endl;
	}
	if(mpi_rank < mpi_size-1) {
		int dummy_Send = 1;
		MPI_Send(&dummy_Send, 1, MPI_INT, mpi_rank+1, 1111, mpi_comm);
	}

	MPI_Barrier(mpi_comm);
	for(int i=0; i<10000; ++i) {}
	MPI_Barrier(mpi_comm);
}
