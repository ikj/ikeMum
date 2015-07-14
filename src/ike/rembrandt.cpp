#include "JoeWithModels.h"
#include "JOE/ADJOINT_FILES/JoeWithModelsAD.h"
#include "adolc.h"

//#include "IkeTurbModel_KOM.h"
//#include "turbModels/TurbModel_KOM.h"
//#include "ADJOINT_FILES/TurbModel_KOM_AD.h"

#include "IkeWithPsALC.h"

#include <sstream>
#include <complex>

#include <stdexcept>
#include <iomanip>


//#define NVAL 7
#define VERY_SMALL_DOUBLE 1.0e-10

#define REMBRANDT_ERROR_CODE -3

// ###########################################################################################
// ------                                                                               ------
// ------                     Dump a Tecplot file of the global-modes                   ------
// ------                                                                               ------
// ###########################################################################################
class Rembrandt : public IkeWithPsALC_AD { //, public IkeRansTurbKOm_AD {
public:
	Rembrandt(char *name) : IkeWithPsALC_AD(name), JoeWithModels(name), JoeWithModels_AD(name), UgpWithCvCompFlow(name) {
		if(mpi_rank==0)
			cout<<"Rembrandt()"<<endl;
		init();
	}
	~Rembrandt() {
		if(mpi_rank==0)
			cout<<"~Rembrandt()"<<endl;
		clear();
	}

	void run() {
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

		// Get xcvMin and xcvMax
		for(int i=0; i<3; ++i) {
			xcvMin[i] = ABSURDLY_BIG_NUMBER;
			xcvMax[i] = -ABSURDLY_BIG_NUMBER;
		}
		for(int icv=0; icv<ncv; ++icv) {
			for(int i=0; i<3; ++i) {
				xcvMin[i] = min(xcvMin[i], x_cv[icv][i]);
				xcvMax[i] = max(xcvMax[i], x_cv[icv][i]);
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
			cout<<"nbocv2_v_global, nbocv2_i, and nbocv2_v are successfully constructed"<<endl<<endl;

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

		// Run the eigendecomposition
		if(mpi_rank==0){
			cout<<"======================="<<endl;
			cout<<"EIGEN-DECOMPOSITION"<<endl;
			cout<<"======================="<<endl;
		}

		NcontrolEqns = getNumControlVar(); 	assert(NcontrolEqns >= 0);

		initialHook();

		runEigenAnalysis();

		finalHook();
	}

	/*
	 * Method: runEigenAnalysis
	 */
	void runEigenAnalysis() {
//		// ===============
//		// Read the matrix
//		// ===============
//		string filename = getStringParam("JACOBIAN_MATRIX_FILENAME", "jacMatrix0.bin");
//		readMatrixBinary<MatComprsedSTL>(filename, jacMatrixSTL, nScal, cvora, nbocv_v_global);

		// =============================
		// Obtain the Jacobian matrix
		// =============================
		NcontrolEqns = getNumControlVar(); 	assert(NcontrolEqns >= 0);
		int nVars = 5+nScal;

//		double *qVec   = new double[ncv*nVars];
//		double *rhsVec = new double[ncv*nVars];
//
//		string filename = getStringParam("Q1_FILENAME", "Q1_PT00001.bin");
//		if(mpi_rank==0) cout<<">> Initialize from = "<<filename<<endl;
//
//		double *lambda0Temp = NULL;
//		if(NcontrolEqns > 0) {
//			lambda      = new double [NcontrolEqns];
//			lambda0Temp = new double [NcontrolEqns];
//		}
//		double arclengthTemp;
//
//		// Read the flow field from a binary file
//		readPsALCdumpedDataSerial(filename, qVec, step, lambda0Temp, lambda, NcontrolEqns, arclengthTemp);
//		if(mpi_rank==0) {
//			for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn)
//				printf("  lambda%d=%.3e", iEqn, lambda[iEqn]);
//			printf("  arclength=%.4e \n", arclengthTemp);
//		}
//
//		for(int icv=0; icv<ncv; ++icv) {
//			int tempInt = icv*nVars;
//			rho[icv] = qVec[tempInt];
//			for(int i=0; i<3; ++i)
//				rhou[icv][i] = qVec[tempInt+1+i];
//			rhoE[icv] = qVec[tempInt+4];
//
//			if(nScal > 0) {
//				for(int iScal=0; iScal<nScal; ++iScal) {
//					double *scalArray = scalarTranspEqVector[iScal].phi;
//					scalArray[icv] = qVec[tempInt+5+iScal];
//				}
//			}
//		}
//		MPI_Barrier(mpi_comm);
//
//		// Calculate the Jacobian matrix
//		if(mpi_rank==0) cout<<">> Calculate Jacobian: "<<endl;
//		calcJacobian1DAD(jacMatrixSTL, rhsVec, debugLevel, 0); // Note: arguments = matrix, rhs, debugLevel, number of AD system parameters
//
//		delete [] lambda0Temp;
//		delete [] qVec;
//		delete [] rhsVec;

		string jacFilename = getStringParam("JAC_MAT_FILENAME", "jacMatrix0.bin");
		if(mpi_rank==0) cout<<">> Matrix from = "<<jacFilename<<endl;

//		writeMatrixBinaryParallel<MatComprsedSTL>(jacFilename, jacMatrixSTL, nScal, cvora, nbocv_v_global);
		readMatrixBinary<MatComprsedSTL>(jacFilename, jacMatrixSTL, nScal, cvora, nbocv_v_global);

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
		int nconv = calcGlobalModes(nev);

		// ===============
		// Post-processing
		// ===============
		if(nconv>0) {
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

			for(int icv=0; icv<ncv; ++icv) {
				int indexTemp = icv*nVars;

				// Tecplot output for the least-stable direct & adjoint modes
				FirstDirectMode_rho[icv]   = directEvecsReal[directLeastStableIndex][indexTemp];
				for(int i=0; i<3; ++i) {
					FirstDirectMode_rhou[icv][i] = directEvecsReal[directLeastStableIndex][indexTemp+1+i];
					FirstDirectMode_vel[icv][i]  = directEvecsReal[directLeastStableIndex][indexTemp+1+i]/FirstDirectMode_rho[icv];
				}
				FirstDirectMode_rhoE[icv]  = directEvecsReal[directLeastStableIndex][indexTemp+4];
//				FirstDirectMode_kine[icv]  = directEvecsReal[directLeastStableIndex][indexTemp+5];
//				FirstDirectMode_omega[icv] = directEvecsReal[directLeastStableIndex][indexTemp+6];

				FirstAdjointMode_rho[icv]   = adjointEvecsReal[adjointLeastStableIndex][indexTemp];
				for(int i=0; i<3; ++i) {
					FirstAdjointMode_rhou[icv][i] = adjointEvecsReal[adjointLeastStableIndex][indexTemp+1+i];
					FirstAdjointMode_vel[icv][i]  = adjointEvecsReal[adjointLeastStableIndex][indexTemp+1+i]/FirstAdjointMode_rho[icv];
				}
				FirstAdjointMode_rhoE[icv]  = adjointEvecsReal[adjointLeastStableIndex][indexTemp+4];
//				FirstAdjointMode_kine[icv]  = adjointEvecsReal[adjointLeastStableIndex][indexTemp+5];
//				FirstAdjointMode_omega[icv] = adjointEvecsReal[adjointLeastStableIndex][indexTemp+6];

				// Structural sensitivity
				rho_sens[icv] 	= fabs(FirstDirectMode_rho[icv] * FirstAdjointMode_rho[icv]);
				for(int i=1; i<3; ++i) {
					rhou_sens[icv][i] = fabs(FirstDirectMode_rhou[icv][i] * FirstAdjointMode_rhou[icv][i]);
					vel_sens[icv][i]  = fabs(FirstDirectMode_vel[icv][i] * FirstAdjointMode_vel[icv][i]);
				}
				rhoE_sens[icv]	= fabs(FirstDirectMode_rhoE[icv] * FirstAdjointMode_rhoE[icv]);
//				kine_sens[icv] 	= fabs(FirstDirectMode_kine[icv] * FirstAdjointMode_kine[icv]);
//				omega_sens[icv]	= fabs(FirstDirectMode_omega[icv] * FirstAdjointMode_omega[icv]);
			}

			updateCvDataG1G2(FirstDirectMode_rho,  REPLACE_DATA);
			updateCvDataG1G2(FirstDirectMode_rhou, REPLACE_ROTATE_DATA);
			updateCvDataG1G2(FirstDirectMode_vel,  REPLACE_ROTATE_DATA);
			updateCvDataG1G2(FirstDirectMode_rhoE, REPLACE_DATA);
//			updateCvDataG1G2(FirstDirectMode_kine, REPLACE_DATA);
//			updateCvDataG1G2(FirstDirectMode_omega,REPLACE_DATA);

			updateCvDataG1G2(FirstAdjointMode_rho,  REPLACE_DATA);
			updateCvDataG1G2(FirstAdjointMode_rhou, REPLACE_ROTATE_DATA);
			updateCvDataG1G2(FirstAdjointMode_vel,  REPLACE_ROTATE_DATA);
			updateCvDataG1G2(FirstAdjointMode_rhoE, REPLACE_DATA);
//			updateCvDataG1G2(FirstAdjointMode_kine, REPLACE_DATA);
//			updateCvDataG1G2(FirstAdjointMode_omega,REPLACE_DATA);

			updateCvDataG1G2(rho_sens,  REPLACE_DATA);
			updateCvDataG1G2(rhou_sens, REPLACE_ROTATE_DATA);
			updateCvDataG1G2(vel_sens,  REPLACE_ROTATE_DATA);
			updateCvDataG1G2(rhoE_sens, REPLACE_DATA);
//			updateCvDataG1G2(kine_sens, REPLACE_DATA);
//			updateCvDataG1G2(omega_sens,REPLACE_DATA);
		}
	}

	/*
	 * Method: initialHook()
	 */
	void initialHook() {
		if(mpi_rank == 0)
			cout<<"Rembrandt::InitialHook() "<<endl;

//		// ==================
//		// Heat release model
//		// ==================
//		fSt   = 0.028; //Stoichiometric fuel/air ratio (0.28 for H2/air)
//		hFuel = 1.20e8; //Fuel heating value (120MJ/kg for H2)
//		mAir  = 0.053*20; //Mass flow rate of injected air
//
//		kComb = 0.95; //Free parameter: Fraction of completed combustion (DLR simulation result: 95%)
//		dComb = getDoubleParam("SHAPE_HEAT_RELEASE",0.75); //Free parameter: Shape of heat release (best result from Vince's experiment - he tried dComb=0.5,0.75,1.0,1.25,1.5)
//
//		xCombStart = 0.408; //Combustion ignition starting position
//
//		xThroat  = 0.64983;
//		xMax     = 0.773; //0.775856;
//		yBottom  = 0.1195;
//		yTop     = 0.1293;
//		yMax     = 0.155571; //0.156088;
//
//		lengthComb = xMax - xCombStart; //Combustor length
//		width      = 0.2; //Combustor width
//
//		if(mpi_rank==0) {
//			cout<<endl<<">> SCRAMJET HEAT RELEASE MODEL: FROM C.Doolan & R.Boyce, AIAA 2008-2603"<<endl;
//			if(NcontrolEqns > 0) {
//				printf("**    lambda (phi) = ");
//				for(int i=0; i<NcontrolEqns; ++i)
//					printf("%.3e ", lambda[i]);
//				cout<<endl;
//			} else {
//				printf("**    lambda (phi) = EMPTY \n");
//			}
//			cout<<"**    SHAPE_HEAT_RELEASE = "<<dComb<<endl;
//			cout<<"**    xCombStart         = "<<xCombStart<<endl<<endl;
//		}
	}

//	/*
//	 * Method: boundaryHookScalarRansTurb
//	 * ----------------------------------
//	 *
//	 */
//	void boundaryHookScalarRansTurb(double *phi, FaZone *zone, const string &name) {
//		if (zone->getKind() == FA_ZONE_BOUNDARY) {
//			Param *param;
//			if (getParam(param, zone->getName())) {
//				if ((param->getString() == "WALL") && (name == "omega")) {
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//						int icv0 = cvofa[ifa][0];
//						int icv1 = cvofa[ifa][1];
//						double muLamCV = calcMuLam(UgpWithCvCompFlow::temp[icv0]);
//						double omegaNew = 6.0*muLamCV/(rho[icv0]*beta0*wallDist[icv0]*wallDist[icv0]);
//
//						phi[icv1] = omegaNew;
//					}
//				}
//			}
//		}
//	}
//
//	/*
//	 * Method: boundaryHookScalarRansTurb1D
//	 * ------------------------------------
//	 * Original code = boundaryHookScalarRansTurb_AD() in IkeUgpWithCvCompFlow.h and boundaryHookScalarRansTurb() in UgpWithCvCompFlowAD.h
//	 * In some situation, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
//	 */
//	void boundaryHookScalarRansTurb1D(const int ifa, double *phi, FaZone *zone, const string &scalName) {
//		if (zone->getKind() == FA_ZONE_BOUNDARY) {
//			Param *param;
//			if (getParam(param, zone->getName())) {
//				if ((param->getString() == "WALL") && (scalName == "omega")) {
//					int icv0 = cvofa[ifa][0];
//					int icv1 = cvofa[ifa][1];
//					double muLamCV = calcMuLam(UgpWithCvCompFlow::temp[icv0]);
//					double omegaNew = 6.0*muLamCV/(rho[icv0]*beta0*wallDist[icv0]*wallDist[icv0]);
//					phi[icv1] = omegaNew;
//				}
//			}
//		}
//	}
//
//	/*
//	 * Method: boundaryHookScalarRansTurb1D_AD
//	 * ---------------------------------------
//	 * Acutally this function doesn't need to be overloaded: Currently just a debugging purpose
//	 */
//#ifdef USE_MEM_SAVING_ADVAR
//	virtual void boundaryHookScalarRansTurb1D_AD(const int ifa, ADscalar<adouble> &phi, FaZone *zone, const string &name) {
//		static bool firstCall = true;
//		int debugLevel = getDebugLevel();
//
//		if (zone->getKind() == FA_ZONE_BOUNDARY) {
//			Param *param;
//			if (getParam(param, zone->getName()))
//				if ((param->getString() == "WALL") && (name == "omega")) {
//					if(debugLevel>0) {
//						if(mpi_rank==0 && firstCall)
//							cout<<"Calling PsALChyShot2D::boundaryHookScalarRansTurb1D_AD() for boundary "<<zone->getName()<<endl;
//						firstCall = false;
//					}
//
//					int icv0 = cvofa[ifa][0];
//					int icv1 = cvofa[ifa][1];
//					adouble omegaNew_AD = 6.0*calcMuLam_AD(icv0) / (rho_AD[icv0]*beta0*wallDist[icv0]*wallDist[icv0]);
//					phi[icv1] = omegaNew_AD;
//				}
//		}
//	}
//#endif
//
//	/*
//	 * Method: boundaryHookScalarRansTurb1D_AD
//	 * ---------------------------------------
//	 * This method is for old IKE calculations: it will NOT be called if USE_MEM_SAVING_ADVAR is used
//	 */
//	virtual void boundaryHookScalarRansTurb1D_AD(const int ifa, adouble *phi, FaZone *zone, const string &scalName)  {
//		if (zone->getKind() == FA_ZONE_BOUNDARY) {
//			Param *param;
//			if (getParam(param, zone->getName())) {
//				if (strcmp(param->getString().c_str(),"WALL")==0 && strcmp(scalName.c_str(),"omega")==0) {
//					int icv0 = cvofa[ifa][0];
//					int icv1 = cvofa[ifa][1];
//					adouble muLamCV_AD = calcMuLam_AD(icv0);
//					adouble omegaNew_AD = 6.0*muLamCV_AD/(rho_AD[icv0]*beta0*wallDist[icv0]*wallDist[icv0]);
//					phi[icv1] = omegaNew_AD;
//				}
//			}
//		}
//	}

protected:
	/*
	 * Method: calcGlobalModes
	 * -----------------------
	 * Calculate both direct and adjoint global-modes
	 */
	int calcGlobalModes(const int nevSlepc) {
		// Get the slepc parameters
		if (!checkParam("SLEPC_EIGEN_PARAM")) {
			ParamMap::add("SLEPC_EIGEN_PARAM  NCV=15  MPD=15  TOL=1.0e-7  MAX_ITER=1000"); // add default values
			if (mpi_rank == 0)
				cout<< "WARNING: added keyword \"SLEPC_EIGEN_PARAM  NCV=15  MPD=15  TOL=1.0e-7  MAX_ITER=1000\""<< " to parameter map!" << endl;
		}
		int ncvSlepc      = getParam("SLEPC_EIGEN_PARAM")->getInt("NCV");
		int mpdSlepc      = getParam("SLEPC_EIGEN_PARAM")->getInt("MPD");
		double tolSlepc   = getParam("SLEPC_EIGEN_PARAM")->getDouble("TOL");
		int max_iterSlepc = getParam("SLEPC_EIGEN_PARAM")->getInt("MAX_ITER");

		if(mpi_rank==0)
			cout<<"Rembrandt::calcGlobalModes(): nev="<<nevSlepc<<", ncv="<<ncvSlepc<<", mpd="<<mpdSlepc<<", tol="<<tolSlepc<<", max_iter="<<max_iterSlepc<<endl;

		// Allocate variables
		double* relError = new double [nevSlepc];
		int numIter;

		// Solve for direct global-modes
		int nconvDirect = petscSolver2->solveEigenProblemSlepc<MatComprsedSTL>(directEvalsReal, directEvalsImag, directEvecsReal, directEvecsImag, relError, numIter,
											jacMatrixSTL, cvora, nbocv_v_global, nScal, ncv_gg,
											nevSlepc, ncvSlepc, mpdSlepc, tolSlepc, max_iterSlepc);

		if(mpi_rank==0) {
			cout<<endl
				<<">> "<<nconvDirect<<" converged direct eigen-pairs after "<<numIter<<"iterations"<<endl;
			cout<<">> ERROR for each eigen-pairs:"<<endl;
			for(int i=0; i<nevSlepc; ++i)
				printf("    %.2d  %.5e \n", i, relError[i]);
			cout<<endl;
		}

		// Solve for adjoint global-modes
		petscSolver2->transpose();
		int nconvAdjoint = petscSolver2->solveEigenProblemSlepc<MatComprsedSTL>(directEvalsReal, directEvalsImag, directEvecsReal, directEvecsImag, relError, numIter,
											jacMatrixSTL, cvora, nbocv_v_global, nScal, ncv_gg,
											nevSlepc, ncvSlepc, mpdSlepc, tolSlepc, max_iterSlepc);

		if(mpi_rank==0) {
			cout<<endl
				<<">> "<<nconvAdjoint<<" converged adjoint eigen-pairs after "<<numIter<<"iterations"<<endl;
			cout<<">> ERROR for each eigen-pairs:"<<endl;
			for(int i=0; i<nevSlepc; ++i)
				printf("    %.2d  %.5e \n", i, relError[i]);
			cout<<endl;
		}

		// Free memory
		delete [] relError;

		// Return
		return int( min(nconvDirect, nconvAdjoint) );
	}

	void init() {
		// Global-modes
		directEvalsReal = NULL;
		directEvalsImag = NULL;
		directEvecsReal = NULL;
		directEvecsImag = NULL;

		adjointEvalsReal = NULL;
		adjointEvalsImag = NULL;
		adjointEvecsReal = NULL;
		adjointEvecsImag = NULL;

		// Structural sensitivity
		rho_sens = NULL;	registerScalar(rho_sens, "SENSITIVITY_RHO", CV_DATA);
		rhou_sens = NULL;	registerVector(rhou_sens, "SENSITIVITY_RHOU", CV_DATA);
		vel_sens = NULL;	registerVector(vel_sens, "SENSITIVITY_VEL", CV_DATA);
		rhoE_sens = NULL;	registerScalar(rhoE_sens, "SENSITIVITY_RHOE", CV_DATA);
//		kine_sens = NULL;	registerScalar(kine_sens, "SENSITIVITY_kine", CV_DATA);
//		omega_sens = NULL;	registerScalar(omega_sens, "SENSITIVITY_omega", CV_DATA);

		// Tecplot output
		FirstDirectMode_rho = NULL; 	registerScalar(FirstDirectMode_rho,  "FirstDirectMode_RHO", CV_DATA);
		FirstDirectMode_rhou = NULL;	registerVector(FirstDirectMode_rhou, "FirstDirectMode_RHOU", CV_DATA);
		FirstDirectMode_vel = NULL; 	registerVector(FirstDirectMode_vel,  "FirstDirectMode_VEL", CV_DATA);
		FirstDirectMode_rhoE = NULL; 	registerScalar(FirstDirectMode_rhoE, "FirstDirectMode_RHOE", CV_DATA);
//		FirstDirectMode_kine = NULL; 	registerScalar(FirstDirectMode_kine, "FirstDirectMode_KINE", CV_DATA);
//		FirstDirectMode_omega = NULL; 	registerScalar(FirstDirectMode_omega,"FirstDirectMode_OMEGA", CV_DATA);

		FirstAdjointMode_rho = NULL; 	registerScalar(FirstAdjointMode_rho,  "FirstAdjointMode_RHO", CV_DATA);
		FirstAdjointMode_rhou = NULL; 	registerVector(FirstAdjointMode_rhou, "FirstAdjointMode_RHOU", CV_DATA);
		FirstAdjointMode_vel = NULL; 	registerVector(FirstAdjointMode_vel,  "FirstAdjointMode_VEL", CV_DATA);
		FirstAdjointMode_rhoE = NULL; 	registerScalar(FirstAdjointMode_rhoE, "FirstAdjointMode_RHOE", CV_DATA);
//		FirstAdjointMode_kine = NULL; 	registerScalar(FirstAdjointMode_kine, "FirstAdjointMode_KINE", CV_DATA);
//		FirstAdjointMode_omega = NULL; 	registerScalar(FirstAdjointMode_omega,"FirstAdjointMode_OMEGA", CV_DATA);
	}

	void clear() {
		// Global-modes
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
	}

protected:
//	// Heat Release Model: Based on C.Doolan & R.Boyce, AIAA 2008-2603
//	double fSt; //Stoichiometric fuel/air ratio (0.28 for H2/air)
//	double hFuel; //Fuel heating value (H2)
//	double mAir; //Mass flow rate of injected air
//
//	double lengthComb; //Combustor length
//	double kComb; //Free parameter: Fraction of completed combustion (BLR exp.: 95%)
//	double dComb; //Free parameter: Shape of heat release
//
//	double xCombStart; //Combustion ignition starting position
//	double xThroat, xMax;
//	double yBottom, yTop, yMax;
//	double width;

	// Number of global modes that the user wants to calculate
	int nev;

	// Global-modes
	double*  directEvalsReal;
	double*  directEvalsImag;
	double** directEvecsReal;
	double** directEvecsImag;

	double*  adjointEvalsReal;
	double*  adjointEvalsImag;
	double** adjointEvecsReal;
	double** adjointEvecsImag;

	// Structural sensitivity
	double *rho_sens;
	double (*rhou_sens)[3];
	double (*vel_sens)[3];
	double *rhoE_sens;
//	double *kine_sens;
//	double *omega_sens;

	// Tecplot output
	double *FirstDirectMode_rho;
	double (*FirstDirectMode_rhou)[3];
	double (*FirstDirectMode_vel)[3];
	double *FirstDirectMode_rhoE;
//	double *FirstDirectMode_kine;
//	double *FirstDirectMode_omega;

	double *FirstAdjointMode_rho;
	double (*FirstAdjointMode_rhou)[3];
	double (*FirstAdjointMode_vel)[3];
	double *FirstAdjointMode_rhoE;
//	double *FirstAdjointMode_kine;
//	double *FirstAdjointMode_omega;
};

// ###########################################################################################
// ------                                                                               ------
// ------                                      MAIN                                     ------
// ------                                                                               ------
// ###########################################################################################

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);

	initMpiStuff();

	// the run number specifies the class which is going to be instantiated
	// set run to default value of 0 => instantiation of JoeWithModels
	int run = 0;
	int ptOnSCurve = 0;

	// set input name to default "Joe.in"
	char inputFileName[50];
	sprintf(inputFileName, "Rembrandt.in");

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

		// declare pointer to JoeWithModels
		JoeWithModels *joe;

		switch (run) {
			case 0:
				joe = new Rembrandt(inputFileName);      break;
			default:
				if (mpi_rank == 0)
					cerr << "ERROR: run number not available!" << endl;
				throw(-1);
		}
    
		// provide total runtime
		double wtime, wtime0;
		MPI_Barrier(mpi_comm);
		if (mpi_rank == 0)
			wtime = MPI_Wtime();

		// run joe
		joe->run();
    
		MPI_Barrier(mpi_comm);
		if (mpi_rank == 0) {
			double wtime0 = wtime;
			wtime = MPI_Wtime();
			cout << " > total runtime [s]: " << wtime - wtime0 << endl;
		}

		// delete joe (make sure memory is deallocated in destructors
		delete joe;
	}

	catch (int e) {
		cerr << "Exception: " << e << endl;
		MPI_Finalize();
		return(-1);
	}
	catch (const invalid_argument& problem) {
		if(mpi_rank==0)
			cout << problem.what() <<endl;
		MPI_Finalize();
		return(-1);
	}
	catch (const runtime_error& problem) {
		if(mpi_rank==0)
			cout << problem.what() <<endl;
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

