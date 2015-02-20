/*
 * standaloneNewton.cpp
 *
 *  Created on: Feb 9, 2015
 *      Author: ikj
 */

#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include <assert.h>
#include <algorithm>
#include <math.h>
#include <vector>

#include "MpiStuff.h"
using namespace MpiStuff;
#include "myMem.h"

using namespace std;

#include "MatComprsed.h"
#include "PetscSolver2.h"

//enum MPI_Datatype {MPI_INT, MPI_DOUBLE};
//enum MPI_Op {MPI_MIN, MPI_MAX, MPI_SUM};

#define BCGSTAB              301
#define PETSC_GMRES          302

#ifndef WITH_PETSC
#define WITH_PETSC
#endif

enum HOW_TO_CALC_JAC {ROW_1D, ORDINARY_2D};
enum MODIFIED_NEWTON {BASIC, SHAMANSKII, MODIFIED_SHAMANSKII};
enum WEIGHT_RHS_METHOD {NO_RHSWEIGHT, REF_VALUES, LOCAL_VALUES, REF_VALUES_AND_DT_OVER_VOL};

#define JAC1D_STATUS_FILENAME  "Jac1D_status.dat"
#define BIFUR_SUMMARY_FILENAME "Bifur_summary.csv"
#define NEWTON_STATUS_FILENAME "Newton_converg.csv"
#define LS_CONVERGENCE_FILENAME "LinearSystem_converg.csv"

#ifndef ABSURDLY_BIG_NUMBER
	#define ABSURDLY_BIG_NUMBER 2.22e22
#endif

#define PSALC_ERROR_CODE -2
#define IKE_BACKTRACKING_CODE -22

#define USE_TOT_NORM_WITH_LAMBDA

#ifndef MACHINE_EPS
#define MACHINE_EPS 2.0e-16 // Note: machine epsilon for double precision = 1.11e-16
#endif

#define RELAXATION_EPS 7.777e-6

class PsalcStandAlone {
public:
	// Variables that should be easily accessed by PsALC
	int ncv, ncv_gg, ncv_ggff;
//	int mpi_rank;
//	int mpi_size;
//	int mpi_comm;

	int nScal;
	int debugLevel;
	int iterNewton;

	int *cvora;

	int *nbocv_i;
//	int nbocv_s;
	int *nbocv_v;

	int *nbocv_i_global;
	int *nbocv_v_global;

	int *nbocv2_i;
	vector<int> nbocv2_v;

	double *rho;
	double (*rhou)[3];
	double *rhoE;

	int linearSolverNS;         ///< linear solver for Navier-Stokes equations

	MatComprsed    jacMatrix;       // This will be used if the Jacobian matrix is calculated as whole
	MatComprsedSTL jacMatrixSTL;

	vector<int> nbocv2_eachIcv;

	PetscSolver2 *petscSolver2;

	int step;

	double* lambda;

	double initResidNorm;
	double residNormTot;
	double residNormTotOld;

	int monitorConvergInterval;
	int gmresRestart;

	int myNonDiagAdd_count;     // Statistics: the number of modified diagonal

public:
	PsalcStandAlone() {
		ncv      = 1;
		ncv_gg   = 1;
		ncv_ggff = 1;

//		mpi_rank = 0;
//		mpi_size = 1;
//		mpi_comm = 0;

		nScal = 0;
		debugLevel = 1;
		iterNewton = 0;

		cvora = new int [2];
		cvora[0] = 0; 	cvora[1] = ncv;

		nbocv_i = new int [2];
		nbocv_i[0] = 0; nbocv_i[1] = 1;
		nbocv_v = new int [1];
		nbocv_v[0] = 0;

		nbocv_i_global = new int [2];
		nbocv_i_global[0] = 0; nbocv_i_global[1] = 1;
		nbocv_v_global = new int [1];
		nbocv_v_global[0] = 0;

		nbocv2_i = new int [2];
		nbocv2_i[0] = 0; nbocv2_i[1] = 1;
		nbocv2_v.resize(1);
		nbocv2_v[0] = 0;

		rho  = new double [ncv];
		rhou = new double [ncv][3];
		rhoE = new double [ncv];

		linearSolverNS = PETSC_GMRES;

		petscSolver2 = NULL;

		step = 0;

		lambda = new double [1];
		lambda[0] = 0.0;

		initResidNorm   = 0.0;
		residNormTot    = 0.0;
		residNormTotOld = 0.0;

		monitorConvergInterval = 10;
		gmresRestart           = 50;

		myNonDiagAdd_count = 0.0;
	}
	~PsalcStandAlone() {
		delete [] cvora;

		delete [] nbocv_i;
		delete [] nbocv_v;

		delete [] nbocv_i_global;
		delete [] nbocv_v_global;

		delete [] nbocv2_i;
		nbocv2_v.clear();

		delete [] rho;
		delete [] rhou;
		delete [] rhoE;

		delete [] lambda;
	}

public:
	// Variables that can come from Ike.in
	struct IkeDotIn {
		// HOW_TO_CALC_JAC
		HOW_TO_CALC_JAC howToCalcJac;

		// LINEAR_SOLVER_NEWTON_TRESHOLDS
		int maxIterLS;
		double zeroAbsLS;
		double zeroRelLS;
		// LSNT_RAMP
		int startDecLS;
		int intervalDecLS;
		int incIterLS;
		int maxFinalIterLS;
		double decZeroLS;
		double minZeroAbsLS;
		double minZeroRelLS;

		// PRECONDITIONER
		string pcType;
		bool pcReuse;
		// ILU_PC_LEVEL
		int pcLevels;

		// STABILIZATION_FOR_SINGULAR_JAC
		int    stabilizationMethodType;
		double stabilizationAlpha;
		double stabilizationAlphaEps;
		int    stabilizationStartingIter;

		// RELAXATION_CLIP_FOR_NEGATIVE_VALS
		double clipParameter;
		double safeParameter;

		// BACKTRACKING_PARAMETERS
		int backtrackMaxIter;
		double backtrackRelax_LowerBound;
		double backtrackRelax_UpperBound;
		double backtrackBarrierCoeff;
		// BACKTRACKING_SKIP
		bool skipBT_firstITer;
		int skipBT_freq;

		// A_PRIORI_WEIGHT_REGULARIZ
		bool useAprioriWeightFunction;
		double digitAprioriWeight;
		double x1AprioriWeight;
		double x2AprioriWeight;

		// WEIGHT_TANGENT_COND
		double weightTangentCond;

		// TRUST_REGION_SIZE
		double trustRegionSize;

		// MORE_NTSTEPS_BELOW_CONV
		int moreNTstepsBelowConv_moreSteps;
		double moreNTstepsBelowConv_resid;
		double moreNTstepsBelowConv_deltaQ;

		// UNDER_RELAXATION
		double relaxation;

		IkeDotIn() {
			// HOW_TO_CALC_JAC
			howToCalcJac = ROW_1D;

			// LINEAR_SOLVER_NEWTON_TRESHOLDS
			maxIterLS = 10;
			zeroAbsLS = 1.0e-6;
			zeroRelLS = 1.0e-6;
			// LSNT_RAMP
			startDecLS = 10;
			intervalDecLS = 1;
			incIterLS = 0;
			maxFinalIterLS = 30;
			decZeroLS = 1.0;
			minZeroAbsLS = 1.0e-8;
			minZeroRelLS = 1.0e-3;

			pcType = "ILU";
			pcReuse = false;
			pcLevels = 2;

			// STABILIZATION_FOR_SINGULAR_JAC
			stabilizationMethodType   = 4;
			stabilizationAlpha        = 0.0;
			stabilizationAlphaEps     = 0.0;
			stabilizationStartingIter = 0;

			// RELAXATION_CLIP_FOR_NEGATIVE_VALS CLIP_PARAM=1.0e-9  SAFE_PARAM=0.8
			clipParameter = 1.0e-9;
			safeParameter = 0.8;

			// BACKTRACKING_PARAMETERS  MAX_ITER=10  RELAX_LOWER_BOUND=0.1  RELAX_UPPER_BOUND=0.5  BARRIER_COEFF=0.01
			backtrackMaxIter          = 3;
			backtrackRelax_LowerBound = 0.3;
			backtrackRelax_UpperBound = 0.8;
			backtrackBarrierCoeff     = 0.01;
			// BACKTRACKING_SKIP  FIRST_ITER=FALSE  FREQUENCY=0
			skipBT_firstITer = false;
			skipBT_freq      = 0;

			// A_PRIORI_WEIGHT_REGULARIZ
			useAprioriWeightFunction = false;
			digitAprioriWeight = 3.0;
			x1AprioriWeight    = 0.35;
			x2AprioriWeight    = 0.4;

			// WEIGHT_TANGENT_COND
			weightTangentCond = 1.0; // Small weight means that the tangential condition is less important

			// TRUST_REGION_SIZE
			trustRegionSize = 1.0e6;

			// MORE_NTSTEPS_BELOW_CONV
			moreNTstepsBelowConv_moreSteps = 0;
			moreNTstepsBelowConv_resid     = 1.0e-15;
			moreNTstepsBelowConv_deltaQ    = 1.0e-15;

			// UNDER_RELAXATION
			relaxation = 1.0;
		}
	};

	IkeDotIn ikeDotIn;

public:
//	// Fake MPI methods
//	int MPI_Barrier(int mpi_comm) {
//		return 1;
//	}
//
//	int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int mpi_comm) {
//		if(datatype == MPI_INT) {
//			int dummy      = *(reinterpret_cast<int*>(&sendbuf));
//			int *recvDummy = reinterpret_cast<int*>(recvbuf);
//			for(int i=0; i<count; ++i)
//				recvbuf[i] = recvDummy[i];
//		} else if(datatype == MPI_DOUBLE) {
//			double dummy      = *(reinterpret_cast<double*>(&sendbuf));
//			double *recvDummy = reinterpret_cast<double*>(recvbuf);
//			for(int i=0; i<count; ++i)
//				recvbuf[i] = recvDummy[i];
//		}
//
//		return 1;
//	}
//
//	int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, int comm) {
//		return 1;
//	}

	int getScalarTransportIndex(string varName) {
		return -1;
	}

	int find2layerCSRstruct_eachIcv(const int icvCenter, vector<int> &nbocv2_v, bool addFakeCells) {
		nbocv2_v.push_back(0);

		return nbocv2_v.size();
	}

	void getPCcontext(string& pcType, bool &pcReuse, int& pcLevels) {
		pcType   = ikeDotIn.pcType;
		pcReuse  = ikeDotIn.pcReuse;
		pcLevels = ikeDotIn.pcLevels;
	}

	template <class MatT>
	void solveLinSysNSCoupled2(double *phi, MatT &A, double *rhs, int &nIter, double &absResid, vector<pair<int, double> >& kspMonitorHistory,
			bool useOldSolnForInitGuess, const double zeroAbs, const double zeroRel, const int maxIter, int nScal, const int monitorConvergInterval,
			int NcontrolParams, int ncv_gg, double **vecC, const double *d, const int step, const int newtonIter) {
		switch (linearSolverNS) {
		case PETSC_GMRES:
			// on the first time, instantiate the petsc solver...
			if (petscSolver2 == NULL) {
				if (nbocv_v_global == NULL) {
					nbocv_v_global = new int[ncv_gg];
					for (int icv = 0; icv < ncv; icv++)
						nbocv_v_global[icv] = cvora[mpi_rank] + icv;
//					updateCvDataG1G2(nbocv_v_global, REPLACE_DATA);
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

				petscSolver2->setGmresRestart(gmresRestart);
			}

			if(useOldSolnForInitGuess)
				petscSolver2->solveGMRES<MatT>(kspMonitorHistory, nIter, absResid, A, phi, rhs, phi,  cvora, nbocv_v_global, 5 + nScal, NcontrolParams, ncv_gg, vecC, d, step, newtonIter);
			else
				petscSolver2->solveGMRES<MatT>(kspMonitorHistory, nIter, absResid, A, phi, rhs, NULL, cvora, nbocv_v_global, 5 + nScal, NcontrolParams, ncv_gg, vecC, d, step, newtonIter);

			break;

		default:
			if (mpi_rank == 0)
				cerr << "Error: unrecognized solver: " << linearSolverNS
						<< endl;
			throw(-1);
			break;
		}
	}

	// Fake methods that make this class same as IkeWithPsALC
	HOW_TO_CALC_JAC getHowToCalcJac() {
		return ikeDotIn.howToCalcJac;
	}

	void getParamsLinSolverNewton(int &maxIterLS, double &zeroAbsLS, double &zeroRelLS,
			int &startDecLS, int &intervalDecLS, int &incIterLS, int &maxFinalIterLS, double &decZeroLS, double &minZeroAbsLS, double &minZeroRelLS,
			const bool firstCall, const int debugLevel) {
		maxIterLS      = ikeDotIn.maxIterLS;
		zeroAbsLS      = ikeDotIn.zeroAbsLS;
		zeroRelLS      = ikeDotIn.zeroRelLS;
		startDecLS     = ikeDotIn.startDecLS;
		intervalDecLS  = ikeDotIn.intervalDecLS;
		incIterLS      = ikeDotIn.incIterLS;
		maxFinalIterLS = ikeDotIn.maxFinalIterLS;
		decZeroLS      = ikeDotIn.decZeroLS;
		minZeroAbsLS   = ikeDotIn.minZeroAbsLS;
		minZeroRelLS   = ikeDotIn.minZeroRelLS;
	}

	void getParamsReducingRelax(double &clipParameter, double &safeParameter, bool showOnScreen) {
		clipParameter = ikeDotIn.clipParameter;
		safeParameter = ikeDotIn.safeParameter;
	}

	void getParamsBacktrackNewton(int &backtrackMaxIter, double &backtrackRelax_LowerBound, double &backtrackRelax_UpperBound, double &backtrackBarrierCoeff,
			bool &skipBT_firstITer, int &skipBT_freq,
			const bool firstCall, const int debugLevel) {
		backtrackMaxIter          = ikeDotIn.backtrackMaxIter;
		backtrackRelax_LowerBound = ikeDotIn.backtrackRelax_LowerBound;
		backtrackRelax_UpperBound = ikeDotIn.backtrackRelax_UpperBound;
		backtrackBarrierCoeff     = ikeDotIn.backtrackBarrierCoeff;
		skipBT_firstITer          = ikeDotIn.skipBT_firstITer;
		skipBT_freq               = ikeDotIn.skipBT_freq;
	}

	void getParamStabilizedNewton(int &stabilizationMethodType, double &stabilizationAlpha, double &stabilizationAlphaEps, int &stabilizationStartingIter,
			const bool firstCall, const int debugLevel) {
		stabilizationMethodType   = ikeDotIn.stabilizationMethodType;
		stabilizationAlpha        = ikeDotIn.stabilizationAlpha;
		stabilizationAlphaEps     = ikeDotIn.stabilizationAlphaEps;
		stabilizationStartingIter = ikeDotIn.stabilizationStartingIter;
	}

	void getParamMoreNewtonSteps(int &moreNTstepsBelowConv_moreSteps, double &moreNTstepsBelowConv_resid, double &moreNTstepsBelowConv_deltaQ,
			const bool firstCall) {
		moreNTstepsBelowConv_moreSteps = ikeDotIn.moreNTstepsBelowConv_moreSteps;
		moreNTstepsBelowConv_resid     = ikeDotIn.moreNTstepsBelowConv_resid;
		moreNTstepsBelowConv_deltaQ    = ikeDotIn.moreNTstepsBelowConv_deltaQ;
	}

	void initialHookNewton() {
		int icv = 0;

		rho[icv]     = 1.1;
		rhou[icv][0] = 1.2;
		rhou[icv][1] = 0.9;
		rhou[icv][2] = 0.0;
		rhoE[icv]    = 3.8;
	}
	void initialHookNewton_firstCall() {}

	void getFlowWeightsForInnerProduct(double* weightVec, const double defaultConst, const int icv, const int nVars) {
		assert(weightVec!=NULL && nVars>0);

		int nScal = nVars-5;

		// Initialize
		for(int i=0; i<nVars; ++i)
			weightVec[i] = defaultConst;
	}

	int calcRhsWithBarrier(double* rhs, const bool useBarrier) {
		int icv = 0;
		// The solution of the following example is [1 1 1 0 4]
		rhs[icv+0] = -2*rho[icv] + rhou[icv][0] + rhou[icv][1] + rhou[icv][2];
		rhs[icv+1] = 5*pow(rhou[icv][0],2) - pow(rhou[icv][0],22)/pow(rhou[icv][1],20) + pow(rhou[icv][2],2) - rhoE[icv]/rho[icv];
		rhs[icv+2] = rhou[icv][0]*rhou[icv][1] + pow(rhou[icv][1],2) + pow(rhou[icv][2],2) - rhoE[icv]/rho[icv] + rho[icv]*exp(pow(rhou[icv][0]/rhou[icv][1],2)-1.0) + 1.0;
		rhs[icv+3] = rhou[icv][2];
		rhs[icv+4] = rho[icv]*rhoE[icv]*rhou[icv][0] + rho[icv]*rhoE[icv]*rhou[icv][1] + pow(rhou[icv][0]-1.0,3) - pow(rhou[icv][1]+1,3) + pow(rhou[icv][2],3) + pow(rhou[icv][0]-rhou[icv][1],33)/pow(rho[icv],30);

		int countNegative;
		for(int icv=0; icv<ncv; ++icv) {
			if(rho[icv] < 0.0)
				countNegative++;
			if(rhoE[icv] < 0.0)
				countNegative++;
		}

		return countNegative;
	}

	int calcJacobian1DAD(MatComprsedSTL &jacMatrixSTL, double *rhsSingleArray, const int debugLevel, const int NcontrolEqns) {
		int icv = 0;

		// The solution of the following example is [1 1 1 0 4]
		rhsSingleArray[icv+0] = -2*rho[icv] + rhou[icv][0] + rhou[icv][1] + rhou[icv][2];
		rhsSingleArray[icv+1] = 5*pow(rhou[icv][0],2) - pow(rhou[icv][0],22)/pow(rhou[icv][1],20) + pow(rhou[icv][2],2) - rhoE[icv]/rho[icv];
		rhsSingleArray[icv+2] = rhou[icv][0]*rhou[icv][1] + pow(rhou[icv][1],2) + pow(rhou[icv][2],2) - rhoE[icv]/rho[icv] + rho[icv]*exp(pow(rhou[icv][0]/rhou[icv][1],2)-1.0) + 1.0;
		rhsSingleArray[icv+3] = rhou[icv][2];
		rhsSingleArray[icv+4] = rho[icv]*rhoE[icv]*rhou[icv][0] + rho[icv]*rhoE[icv]*rhou[icv][1] + pow(rhou[icv][0]-1.0,3) - pow(rhou[icv][1]+1,3) + pow(rhou[icv][2],3) + pow(rhou[icv][0]-rhou[icv][1],33)/pow(rho[icv],30);

		// Jacobian
		int nVars = 5+nScal;

		assert(nbocv2_eachIcv.empty());

		nbocv2_eachIcv.push_back(icv);

		jacMatrixSTL.setMatSize(ncv_gg*nVars+NcontrolEqns, ncv*nVars);
		jacMatrixSTL.initMat();

		int nnz = 20;
		unsigned int* rind   = new unsigned int [nnz];
		unsigned int* cind   = new unsigned int [nnz];
		double* values = new double [nnz];

		double r  = rho[icv];
		double m1 = rhou[icv][0],  m2 = rhou[icv][1],  m3 = rhou[icv][2];
		double E  = rhoE[icv];

		int row = 0;
		int s = 0; // starting index
		{
			for (int col=0; col<4; ++col) { // First row has only 4 non-zero entries
				rind[s+col] = row;
				cind[s+col] = col;
			}
			values[s+0] = -2;
			values[s+1] =  1;
			values[s+2] =  1;
			values[s+3] =  1;
		}
		{
			row = 1;
			s += 4;
			for (int col=0; col<5; ++col) { // First row has only 4 non-zero entries
				rind[s+col] = row;
				cind[s+col] = col;
			}
			values[s+0] = E/pow(r,2);
			values[s+1] = 10*m1 - 22*pow(m1,21)/pow(m2,20);
			values[s+2] = 20*pow(m1,22)/pow(m2,21);
			values[s+3] = 2*m3;
			values[s+4] = -1/r;
		}
		{
			row = 2;
			s += 5;
			for (int col=0; col<5; ++col) { // First row has only 4 non-zero entries
				rind[s+col] = row;
				cind[s+col] = col;
			}
			values[s+0] = E/pow(r,2) + exp(pow(m1/m2,2)-1);
			values[s+1] = m2 + 2*r*m1/m2*exp(pow(m1/m2,2)-1);
			values[s+2] = m1 + 2*m2 - 2*r*pow(m1,2)/pow(m2,3)*exp(pow(m1/m2,2)-1);
			values[s+3] = 2*m3;
			values[s+4] = -1/r;
		}
		{
			row = 3;
			s += 5;

			rind[s]   = row;
			cind[s]   = 3;
			values[s] = 1;
		}
		{
			row = 4;
			s += 1;
			for (int col=0; col<5; ++col) { // First row has only 4 non-zero entries
				rind[s+col] = row;
				cind[s+col] = col;
			}
			values[s+0] = E*m1 + E*m2 - 30*pow(m1-m2,33)/pow(r,31);
			values[s+1] = r*E + 3*pow(m1-1,2) + 33*pow(m1-m2,32)/pow(r,30);
			values[s+2] = r*E - 3*pow(m2+1,2) - 33*pow(m1-m2,32)/pow(r,30);
			values[s+3] = 3*pow(m3,2);
			values[s+4] = r*m1+r*m2;
		}



		int firstRindReal = icv*nVars;
		jacMatrixSTL.blockCopy(nnz, firstRindReal, rind, cind, values);

		jacMatrixSTL.convert_blockCindices(nbocv2_eachIcv, jacMatrixSTL.get_nnz()-nnz, nnz, nVars, nVars, NcontrolEqns, ncv_gg);
		jacMatrixSTL.finalizeMat();

		delete [] rind;
		delete [] cind;
		delete [] values;

		nbocv2_eachIcv.clear();

		return 0;
	}

	bool calcNres(double* Nres, const double *q, double* rhs, const int Nvars, const int NcontrolEqns, const double *q1, double** q_tangent,
			const double* lambda_tangent, const double *lambda1, const double weightLambda, const double arcLength) {
		return true;
	}

	bool calcNresForJOE(double* Nres, const int Nvars, const int NcontrolEqns, const double *q1, double** q_tangent,
			const double* lambda_tangent, const double *lambdaNew, const double *lambda1, const double weightLambda, const double arcLength) {
		return true;
	}

	void freeMemGetSteadySolnByNewton(double *phi, double* residNorm, double* residNormOld, double *Nres,
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

	void calcResidualsFrom1Drhs(double *Residual, double *rhs1Darray, const int whichNorm) {
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
				cout<<"ERROR in calcResidualsFrom1Drhs(): unsupported norm = "<<whichNorm<<endl;
			throw(PSALC_ERROR_CODE);
			break;
		}

		delete [] myResidual;
	}

	double calcTotResidual(const double *Residual, const int whichNorm) {
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
				cout<<"ERROR in calcTotResidual(): unsupported norm = "<<whichNorm<<endl;
			throw(PSALC_ERROR_CODE);
			break;
		}

		return residTot;
	}

	double updateTotResidualWithNres(double residNormTot, const double* Nres, const int NcontrolEqns, const int whichNorm) {
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

	void updateFlow_primVars(const double* Q, const double* delQ, const double relaxation, const int nScal) {
		int m = 5 + nScal;

		/* update the flow field */
		for(int icv=0; icv<ncv; ++icv) {
			int icvTemp = icv*m;
			rho[icv] = Q[icvTemp] + relaxation*delQ[icvTemp];
			for(int i=0; i<3; ++i)
				rhou[icv][i] = Q[icvTemp+1+i] + relaxation*delQ[icvTemp+1+i];
			rhoE[icv] = Q[icvTemp+4] + relaxation*delQ[icvTemp+4];
//			for(int iScal=0; iScal<nScal; ++iScal) {
//				double *phi = scalarTranspEqVector[iScal].phi;
//				phi[icv] = Q[icvTemp+5+iScal] + relaxation*delQ[icvTemp+5+iScal];
//			}
		}
	}

	void updateLambda(const double* Q, const double* delQ, const double relaxation, const int NcontrolParams, const int nScal) {
		if(mpi_rank==mpi_size-1) {
			int Nflow = ncv*(5+nScal);
			for(int i=0; i<NcontrolParams; ++i) {
				if(fabs(relaxation*delQ[Nflow+i]) > 0.1*fabs(lambda[i]))
					cout<<"           WARNING in updateLambda(): delta_lambda["<<i<<"] (="<<relaxation*delQ[Nflow+i]<<") is greater than 10% of lambda["<<i<<"] (="<<lambda[i]<<")"<<endl;
				else if(fabs(relaxation*delQ[Nflow+i]) < 1.0e-10)
					cout<<"           WARNING in updateLambda(): delta_lambda["<<i<<"] (="<<relaxation*delQ[Nflow+i]<<") is too small"<<endl;
				lambda[i] = Q[Nflow+i] + relaxation*delQ[Nflow+i];
			}
		}
		MPI_Bcast(lambda, NcontrolParams, MPI_DOUBLE, mpi_size-1, mpi_comm);
	}

	void updateFlow_1Dvec(double* Q, const double* delQ, const double relaxation, const int N) {
		/* increase Q by delQ */
		for(int i=0; i<N; ++i)
			Q[i] = Q[i] + relaxation*delQ[i]; // be careful when you call this method while doing PsALC: "sign" should be negative (i.e. -1.0)
	}

	bool checkBacktrackFromResidInc(const double residNormTotOld, const double residNormTot, const double relaxation, const double backtrackBarrierCoeff) {
		assert(backtrackBarrierCoeff*relaxation>0.0 && backtrackBarrierCoeff*relaxation<1.0);

		double threshold = calcBacktrackThreshold(residNormTotOld, relaxation, backtrackBarrierCoeff);

		return !(residNormTot < threshold);
	}

	double calcBacktrackThreshold(const double residNormTotOld, const double relaxation, const double backtrackBarrierCoeff) {
		return (1.0-backtrackBarrierCoeff*relaxation)*residNormTotOld;
	}

	double backtrackForNegativeVals(int &negativeValCount_CV, int &negativeValCount_FA, const double clipParameter, const double safeParameter,
	                                const double relaxation, const int kine_index, const int nScal, const double* qArray, const double* delQ, const int iterMax) {
		double newRelaxation = relaxation;

		int m = 5 + nScal;

		double rhoClipBound  = clipParameter * 1.0;
		double rhoEClipBound = clipParameter * 4.0;

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

	        if(rhoE[icv] < rhoEClipBound) {
	            double oldRhoE = qArray[icvTemp+4];

	            double tempRelax = safeParameter*(oldRhoE-rhoEClipBound)/fabs(delQ[icvTemp+4]);
	            myRelax = min(myRelax, tempRelax);
	            ++myNegativeCount_CV;

	            negativeFound = true;
	        }
	    }

	    negativeValCount_CV = 0;
	    MPI_Allreduce(&myNegativeCount_CV, &negativeValCount_CV, 1, MPI_INT, MPI_SUM, mpi_comm);
assert(myNegativeCount_CV == negativeValCount_CV);
	    if(negativeValCount_CV>0) {
	        MPI_Allreduce(&myRelax, &newRelaxation, 1, MPI_DOUBLE, MPI_MIN, mpi_comm);
assert(myRelax == myRelax);
	        newRelaxation = min(relaxation, max(newRelaxation, clipParameter)); // To prevent 0.0: Note that clipParameter is usually a very small number
	    }

	    negativeValCount_FA = 0;

	    /* Return */
	    return newRelaxation;
	}

	double backtrackWithJOE_calcRelaxAndRHS(double* rhs, const double* qArray, const double* delQ, bool &notConverged,
			const int maxBacktrackingIter, const double relaxationOld, const double relaxationLowerBound, const double relaxationUpperBound,
			const double* ResidVecOld, double* ResidVecNew, const int whichNorm, const double backtrackBarrierCoeff,
			const int NcontrolEqns, double** q_tangent, const double* lambda_tangent, const double weightLambda, const double WeightTangentCond,
			const double* qArrayOld, const double* NresOld, const double* lambdaOld, const double arcLength) {
		// Store the history of backtracking
		vector<pair<double, double> > backtrackHistory; // This vector stores the history of <relax, residNormTot> from relax1 to the final backtracking
		                                                // Note: If negative values (e.g. rho) occurs, then ABSURDLY_BIG_NUMBER will be stored instead of residNormTot.

		// ----------------------------
		int nVars = 5 + nScal;
		int N = ncv*nVars; // size of the flow vector (Q and delQ)

		double relax0 = 0.0;
		double relax1 = relaxationOld; // Usually, 1.0. But if other relaxation algorithm (e.g. trust-region) was active, it can be less than 1.0.
		double relax2 = relaxationOld*relaxationUpperBound;

		double *ResidVecTemp = new double [nVars];
		double *NresTemp = NULL;
		if (NcontrolEqns>0)
			NresTemp = new double [NcontrolEqns];

		// ----------------------------
		double residNormFlow0 = calcTotResidual(ResidVecOld, whichNorm);
		double residNormTot0;
		if (NcontrolEqns > 0)
			residNormTot0 = updateTotResidualWithNres(residNormFlow0, NresOld, NcontrolEqns, whichNorm);
		else
			residNormTot0 = residNormFlow0;

		double residNormFlow1 = calcTotResidual(ResidVecNew, whichNorm);
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

		double residNormFlow2 = ABSURDLY_BIG_NUMBER; // Initialization with a very large value for debugging
		double residNormTot2  = ABSURDLY_BIG_NUMBER; // Initialization with a very large value for debugging

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
		calcResidualsFrom1Drhs(ResidVecTemp, rhs, whichNorm); // Residual in terms of each flow variable (size = nvars)
		residNormFlow2 = calcTotResidual(ResidVecTemp, whichNorm);

	#ifdef USE_TOT_NORM_WITH_LAMBDA
		if (NcontrolEqns > 0) {
			bool tangentialSatisfied = calcNresForJOE(NresTemp, nVars, NcontrolEqns, qArrayOld, q_tangent, lambda_tangent, lambda, lambdaOld, weightLambda, arcLength);
			for(int iParam=0; iParam<NcontrolEqns; ++iParam)
				NresTemp[iParam] *= WeightTangentCond;

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

		bool notConverged_check1;
//		if(useWatchdogBacktrack)
//			notConverged_check1 = checkBacktrackFromResidIncUsingWatchdog(residNormTot1, residNormTot2, relax2, gammaWatchdog, NcontrolEqns, 5+nScal, rhs, delQ);
//		else
			notConverged_check1 = checkBacktrackFromResidInc             (residNormTot1, residNormTot2, relax2, backtrackBarrierCoeff);

		bool notConverged_check2;
//		if(useWatchdogBacktrack)
//			notConverged_check2 = checkBacktrackFromResidIncUsingWatchdog(residNormTotOld, residNormTot2, relax2, gammaWatchdog, NcontrolEqns, 5+nScal, rhs, delQ);
//		else
			notConverged_check2 = checkBacktrackFromResidInc             (residNormTotOld, residNormTot2, relax2, backtrackBarrierCoeff);

		// Make sure that the final result is less than last step
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
				residNormFlow2 = calcTotResidual(ResidVecTemp, whichNorm);
	#ifdef USE_TOT_NORM_WITH_LAMBDA
				if (NcontrolEqns > 0) {
					bool tangentialSatisfied = calcNresForJOE(NresTemp, nVars, NcontrolEqns, qArrayOld, q_tangent, lambda_tangent, lambda, lambdaOld, weightLambda, arcLength);
					for(int iParam=0; iParam<NcontrolEqns; ++iParam)
						NresTemp[iParam] *= WeightTangentCond;

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
			residNormFlow2 = calcTotResidual(ResidVecTemp, whichNorm);
	#ifdef USE_TOT_NORM_WITH_LAMBDA
			if (NcontrolEqns > 0) {
				bool tangentialSatisfied = calcNresForJOE(NresTemp, nVars, NcontrolEqns, qArrayOld, q_tangent, lambda_tangent, lambda, lambdaOld, weightLambda, arcLength);
				for(int iParam=0; iParam<NcontrolEqns; ++iParam)
					NresTemp[iParam] *= WeightTangentCond;

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
			 *   Creterion: 1. Residual decreases compared to the previous backtracking step
			 *              2. Residual is smaller than the initial residual (Important !!)
			 */
			bool notConverged_check1;
//			if(useWatchdogBacktrack)
//				notConverged_check1 = checkBacktrackFromResidIncUsingWatchdog(residNormTot1, residNormTot2, relax2, gammaWatchdog, NcontrolEqns, 5+nScal, rhs, delQ);
//			else
				notConverged_check1 = checkBacktrackFromResidInc             (residNormTot1, residNormTot2, relax2, backtrackBarrierCoeff);

			bool notConverged_check2; // Make sure that the final result is less than last step
//			if(useWatchdogBacktrack)
//				notConverged_check2 = checkBacktrackFromResidIncUsingWatchdog(residNormTotOld, residNormTot2, relax2, gammaWatchdog, NcontrolEqns, 5+nScal, rhs, delQ);
//			else
				notConverged_check2 = checkBacktrackFromResidInc             (residNormTotOld, residNormTot2, relax2, backtrackBarrierCoeff);

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
			residNormFlow2 = calcTotResidual(ResidVecTemp, whichNorm);
	#ifdef USE_TOT_NORM_WITH_LAMBDA
			if (NcontrolEqns > 0) {
				bool tangentialSatisfied = calcNresForJOE(NresTemp, nVars, NcontrolEqns, qArrayOld, q_tangent, lambda_tangent, lambda, lambdaOld, weightLambda, arcLength);
				for(int iParam=0; iParam<NcontrolEqns; ++iParam)
					NresTemp[iParam] *= WeightTangentCond;

				residNormTot2 = updateTotResidualWithNres(residNormFlow2, NresTemp, NcontrolEqns, whichNorm);
			} else
				residNormTot2 = residNormFlow2;
	#else
			residNormTot2 = residNormFlow2;
	#endif

			// Show the result on the screen
			if(debugLevel > 0 && mpi_rank == 0) {
				double threshold;
//				if(useWatchdogBacktrack)
//					threshold = calcBacktrackThresholdUsingWatchdog(residNormTotOld, relax2, gammaWatchdog, NcontrolEqns, 5+nScal, rhs, delQ);
//				else
					threshold = calcBacktrackThreshold(residNormTotOld, relax2, backtrackBarrierCoeff);

				printf("                  >> backtrackWithJOE_calcRelaxAndRHS(): Backtracking result: iter = %d  RELAX=%.3e  RESID_TOT=%.5e --> %.5e (> THRESHOLD=%.4e)\n",
						minResidIter, relax2, residNormTotInit, residNormTot2, threshold);
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

public:
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
	void getSteadySolnByNewton(double *q, double* rhs, const int maxIterNewton, const double absTolNewton, const double relTolNewton,
			const bool writeJacOnFile = false,
			const int NcontrolEqns = 0, const double *q1 = NULL, double** q_tangent = NULL,
			double *lambda_tangent = NULL, const double *lambda0 = NULL, const double *lambda1 = NULL, const double weightLambda = 0.0,
			const double arcLength = 0) {
		static bool firstCall = true;

		assert(q!=NULL && rhs!=NULL);
		if(NcontrolEqns==0)
			assert(q1==NULL && q_tangent==NULL && lambda_tangent==NULL && lambda1==NULL);
		else
			assert(q1!=NULL && q_tangent!=NULL && lambda_tangent!=NULL && lambda1!=NULL && arcLength!=0.0);

		int m = nScal+5; // number of variables (rho, rhou, rhov, rhow, rhoE, and scalars)

		int countNegative = 0; // number of negative rho or press. If this is greater than 0, the backtracking line search becomes active

		/********************/
		/*  Initialization  *
		 ********************/
		/* Read parameters from input file */
		// How to calculate the Jacobian matrix (ROW_1D, ORDINARY_2D): For a better maintenance of the code, call a separated function
		HOW_TO_CALC_JAC howToCalcJac = getHowToCalcJac();

		// Linear solver setting (PETSc)
		int maxIterLS; 		// "LINEAR_SOLVER_NEWTON_TRESHOLDS"-->"MAX_ITER" in the Ike.in file
		double zeroAbsLS; 	// "LINEAR_SOLVER_NEWTON_TRESHOLDS"-->"ABS_RESID"
		double zeroRelLS; 	// "LINEAR_SOLVER_NEWTON_TRESHOLDS"-->"REL_RESID"
		int startDecLS; 		// "LSNT_RAMP"-->"AFTER_NEWTON_ITER"
		int intervalDecLS; 		// "LSNT_RAMP"-->"INTERVAL_NEWTON_ITER"
		int incIterLS; 			// "LSNT_RAMP"-->"FACTOR_ITER"
		int maxFinalIterLS; 	// "LSNT_RAMP"-->"MAX_ITER"
		double decZeroLS; 		// "LSNT_RAMP"-->"FACTOR_RESID"
		double minZeroAbsLS; 	// "LSNT_RAMP"-->"MIN_ABS_RESID"
		double minZeroRelLS; 	// "LSNT_RAMP"-->"MIN_REL_RESID"
		getParamsLinSolverNewton(maxIterLS, zeroAbsLS, zeroRelLS,
				startDecLS, intervalDecLS, incIterLS, maxFinalIterLS, decZeroLS, minZeroAbsLS, minZeroRelLS,
				firstCall, debugLevel);

		// Modified Newton method
		MODIFIED_NEWTON modifiedNewtonMethod = BASIC;
		double modifNewtonRelResid = 0.0;
		int    modifNewtonFreq     = 0;

		// RHS weight
		WEIGHT_RHS_METHOD weightRhsMethod = REF_VALUES;

		// Relaxation reducing based on negative values
		double clipParameter; // "RELAXATION_CLIP_THRESHOLD"
		double safeParameter; // "RELAXATION_CLIP_SAFETY"
		getParamsReducingRelax(clipParameter, safeParameter, firstCall);

		// Backtracking
		int backtrackMaxIter; 				// "BACKTRACKING_PARAMETERS"-->"MAX_ITER"
		double backtrackRelax_LowerBound; 	// "BACKTRACKING_PARAMETERS"-->"RELAX_LOWER_BOUND"
		double backtrackRelax_UpperBound; 	// "BACKTRACKING_PARAMETERS"-->"RELAX_UPPER_BOUND"
		double backtrackBarrierCoeff; 		// "BACKTRACKING_PARAMETERS"-->"BARRIER_COEFF"
		bool skipBT_firstITer;   // "BACKTRACKING_SKIP"-->"FIRST_ITER"
		int  skipBT_freq;        // "BACKTRACKING_SKIP"-->"FREQUENCY"
		getParamsBacktrackNewton(backtrackMaxIter, backtrackRelax_LowerBound, backtrackRelax_UpperBound, backtrackBarrierCoeff, skipBT_firstITer, skipBT_freq, firstCall, debugLevel);

		//  Stabilized Newton's method for the problem with a singular Jacobian
		int    stabilizationMethodType;   // "STABILIZATION_FOR_SINGULAR_JAC"-->"MATRIX_TYPE"
		double stabilizationAlpha;        // "STABILIZATION_FOR_SINGULAR_JAC"-->"ALPHA"
		double stabilizationAlphaEps;     // "STABILIZATION_FOR_SINGULAR_JAC"-->"ALPHA_EPS"
		int    stabilizationStartingIter; // "STABILIZATION_FOR_SINGULAR_JAC"-->"STARTING_ITER"
		getParamStabilizedNewton(stabilizationMethodType, stabilizationAlpha, stabilizationAlphaEps, stabilizationStartingIter, firstCall, debugLevel);

		// A-priori regularization
		//      Function: If x<=x_2, w = 10^(-digit) + (1.0-10^(-digit)) * exp(-beta*(x-x2)^2)
		//                Otherwise, w = 1.0
		//      Note: 1. The function is uni-directional Gaussian.
		//            2. Maximum function value is 1.0, and Minimum function value is 10^(-digit)
		//      If beta = digit*ln(10) / (x_2-x_1)^2, w(x_1) = 2*10^(-digit)
		//                                            w(x_2) = 1.0
		bool useAprioriWeightFunction = ikeDotIn.useAprioriWeightFunction;
		double digitAprioriWeight     = ikeDotIn.digitAprioriWeight;
		double x1AprioriWeight        = ikeDotIn.x1AprioriWeight;
		double x2AprioriWeight        = ikeDotIn.x2AprioriWeight;
		double betaAprioriWeight = digitAprioriWeight*2.302585093/pow(x2AprioriWeight - x1AprioriWeight, 2.0); // Note: ln(10) = 2.302585093

		// Weight in the tangential condition
		double weightTangentCond = ikeDotIn.weightTangentCond; // Small weight means that the tangential condition is less important

		// Trust-region size
		double trustRegionSize = ikeDotIn.trustRegionSize;

		// Run the interations for few more steps even after they converges
		int moreNTstepsBelowConv_moreSteps;
		double moreNTstepsBelowConv_resid, moreNTstepsBelowConv_deltaQ;
		getParamMoreNewtonSteps(moreNTstepsBelowConv_moreSteps, moreNTstepsBelowConv_resid, moreNTstepsBelowConv_deltaQ, firstCall);

		/* Initialize q */
		//	for(int icv=0; icv<ncv; ++icv) {
		//		int tempInt = icv*m;
		//		q[tempInt] = rho[icv];
		//		for(int i=0; i<3; ++i)
		//			q[tempInt+1+i] = rhou[icv][i];
		//		q[tempInt+4] = rhoE[icv];
		//
		//		if(nScal > 0) {
		//			for(int iScal=0; iScal<nScal; ++iScal) {
		//				double *scalArray = scalarTranspEqVector[iScal].phi;
		//				q[tempInt+5+iScal] = scalArray[icv];
		//			}
		//		}
		//	}
		//	if(NcontrolEqns>0 && mpi_rank==mpi_size-1) {
		//		for(int iParam=0; iParam<NcontrolEqns; ++iParam)
		//			q[ncv*m+iParam] = lambda[iParam];
		//	}
		//	MPI_Barrier(mpi_comm);

		/* Initialize rhs */
		for(int i=0; i<ncv*m; ++i)
			rhs[i] = 0.0;
		if(mpi_rank==mpi_size-1)
			for(int iParam=0; iParam<NcontrolEqns; ++iParam)
				rhs[ncv*(5+nScal)+iParam] = 0.0;

		/* Assign an array for the SAS-form stabilized Newton's method */
		double *sasScaling   = NULL; // Scaling is required for the diagonal of the Jacobian
		double *sasScaledPHI = NULL;
		if(stabilizationMethodType == 2 || stabilizationMethodType == 5) { // Only for the 2nd or 5th method
			if(mpi_rank != mpi_size-1)
				sasScaling = new double [ncv*(5+nScal)];
			else
				sasScaling = new double [ncv*(5+nScal) + NcontrolEqns];
		}
		if(stabilizationMethodType == 5) { // Only for the 5th method
			if(mpi_rank != mpi_size-1)
				sasScaledPHI = new double [ncv*(5+nScal)];
			else
				sasScaledPHI = new double [ncv*(5+nScal) + NcontrolEqns];
		}

		// Set iterNewton to zero (Some messages will be printed in some methods if iterNewton == 0. Also the barrier works correctly for the first time only if iterNewton is reset)
		iterNewton = 0;

		int startingNewtonIter = 0;
//		if(initThirdFromQ1) {
//			if(step-startingStep==2) {
//				startingNewtonIter = getIntParam("INIT_THIRD_NEWTON_ITERS", "0");
//				iterNewton += startingNewtonIter;
//
//				MPI_Barrier(mpi_comm);
//
//				if(firstCall && mpi_rank==0)
//					cout<<"> For the THIRD point, iterNewton starts from "<<startingNewtonIter<<endl;
//			}
//		}

		// InitialHookNewton() & initialHookNewton_firstCall()
		initialHookNewton();
		if(firstCall)
			initialHookNewton_firstCall();

		/*********************************************************/
		/*  Check if RHS's by ADOLC and normal JOE are the same  *
		 *********************************************************/
//		if(firstCall) {
//			string temp = getStringParam("SKIP_COMPATIBILITY_CHECK", "NO");
//			if(temp=="YES" || temp=="yes" || temp=="Y" || temp=="y") {
//				if(mpi_rank==0)
//					cout<<endl
//					<<"> CAUTION!! SKIP_COMPATIBILITY_CHECK = "<<temp<<": skip the compatibility-check between IKE and JOE"<<endl;
//			} else {
//				bool negligibleError = true;
//				if(firstCall)
//					negligibleError = compatibilityCheck();
//				if(!negligibleError) {
//					if(mpi_rank==0)
//						cout<<"ERROR in IkeWithPsALC_AD::getSteadySolnByNewton(): compatibility check by compatibilityCheck() FAILs"<<endl;
//					throw(PSALC_ERROR_CODE);
//				}
//			}
//		}

		/******************************************/
		/*  Newton iteration:                     *
		 *    Solve Q_{n+1} = Q_{n} - L^{-1}*RES  *
		 ******************************************/
		// Settings for the modified Newton's method
		bool useExactJac = true, useExactJacTempo = true;
		int countModifiedNewton = 0;
		int myNnzJac = 0, nnzJac = 0;

		// Allocate and initialize phi
		double *phi; // note: this array will be used while solving the linear system
		if(mpi_rank==mpi_size-1)
			phi = new double [ncv*m + NcontrolEqns];
		else
			phi = new double [ncv*m];

		if(q1 != NULL) {
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

		double* residNormVec    = new double [5+nScal];
		residNormTot            = ABSURDLY_BIG_NUMBER;
		double* residNormVecOld = new double [5+nScal];
		residNormTotOld         = ABSURDLY_BIG_NUMBER;

		// The residuals of the tangential conditions (For the details, google "pseudo-arclength continuation method")
		double* Nres = NULL; // Note: 1. Nres = q_tangent*(q_guess-q1) + lambda_tangent*(lambda_guess-lambda1) - ds
		//       2. Nres must be zero if the tangential condition is satisfied.
		if(NcontrolEqns>0) {
			Nres = new double [NcontrolEqns];
			for(int iParam=0; iParam<NcontrolEqns>0; ++iParam)
				Nres[iParam] = -ABSURDLY_BIG_NUMBER; // Initialize with a very small number for a debugging purpose
		}

		// Weighted tangential vector for q:
		//   1. Two options -- USE_VOLUME_WEIGHTED_INNER_PROD and USE_REF_WEIGHTED_INNER_PROD defined in IkeWithPsALC.h
		//                     If neither USE_VOLUME_WEIGHTED_INNER_PROD nor USE_REF_WEIGHTED_INNER_PROD is defined, no weight (i.e., weight == 1.0).
		//   2. The weighted_q_tangent vector will be used only for the construction of the Jacobian matrix of the pseudo-arclength cont. method
		//
		//   Note that this part must be the same as what can be found in vecDotVecWithWeight() -- Actually, getFlowWeightsForInnerProduct() called by vecDotVecWithWeight()
		//   since the arclength will be calculated by calling vecDotVecWithWeight()
		double **weighted_q_tangent = NULL;

		double rhoWeight  = 1.0;
		double rhouWeight = 1.0;
		double rhoEWeight = 1.0;
		double *turbScalWeightVec = NULL;

		if(q_tangent != NULL) {
			getMem2D(&weighted_q_tangent, 0, NcontrolEqns-1, 0, (ncv*m)-1, "weighted_q_tangent");
			// The size of weighted_q_tangent is same as that of q_tangent

			/* Allocate memory */
			double* weightVec = new double [m];

			/* Calculate weighted q_tangent */
			for(int iParam=0; iParam<NcontrolEqns; ++iParam)
				for(int icv=0; icv<ncv; ++icv) {
					int indexStart = icv*m;

#ifdef USE_VOLUME_WEIGHTED_INNER_PROD
					getFlowWeightsForInnerProduct(weightVec, sqrt(max(totMinCvVol,MIN_CV_VOL_CUTOFF)), icv, m);
#else
					getFlowWeightsForInnerProduct(weightVec, 1.0,                                      icv, m);
#endif

					for(int ivar=0; ivar<m; ++ivar) {
						int index = indexStart + ivar;
						weighted_q_tangent[iParam][index] = weightTangentCond * weightVec[ivar] * q_tangent[iParam][index];
						// If neither USE_VOLUME_WEIGHTED_INNER_PROD nor USE_REF_WEIGHTED_INNER_PROD is defined, weighted_q_tangent becomes same as q_tangent
					}
				}

			/* Free memory */
			delete [] weightVec;
		}

		// Weighted tangential vector for lambda:  weighted_lambda_tangent = weightLambda * lambda_tangent
		//                                         This vector will be used only for the construction of the Jacobian matrix of the pseudo-arclength cont. method
		double* weighted_lambda_tangent = NULL;
		if(lambda_tangent != NULL) {
			weighted_lambda_tangent = new double [NcontrolEqns];
			for(int iParam=0; iParam<NcontrolEqns>0; ++iParam)
				weighted_lambda_tangent[iParam] = weightTangentCond * weightLambda * lambda_tangent[iParam];
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
			}
//			else if(howToCalcJac == ORDINARY_2D) {
//				countNegative = calcJacobianAD(jacMatrix, rhs, debugLevel, NcontrolEqns);
//				myNnzJac      = jacMatrix.get_nnz();
//			}
			else { // ERROR
				freeMemGetSteadySolnByNewton(phi, residNormVec, residNormVecOld, Nres,
						turbScalWeightVec, weighted_q_tangent, weighted_lambda_tangent, m, NcontrolEqns, jacMatrix, jacMatrixSTL);
				assert(true);
			}
			MPI_Allreduce(&myNnzJac, &nnzJac, 1, MPI_INT, MPI_SUM, mpi_comm);
			nnzJac = myNnzJac;

			calcResidualsFrom1Drhs(residNormVec, rhs, whichNorm);
			residNormTot = calcTotResidual(residNormVec, whichNorm);
			// Note: The actual total residual is the sum of residual from the flow variables plus and the residual from the tangential condition.
			//       "residNormTot" will be updated soon after calculating Nres (the residual from the tangential condition).
			//       However, Nres should not change residNormTot much since Nres must be very close to zero
		}
		catch(int e) { // Catch and re-throw
			freeMemGetSteadySolnByNewton(phi, residNormVec, residNormVecOld, Nres,
					turbScalWeightVec, weighted_q_tangent, weighted_lambda_tangent, m, NcontrolEqns, jacMatrix, jacMatrixSTL);
			throw(e);
		}
		catch(...) { // Catch and re-throw
			freeMemGetSteadySolnByNewton(phi, residNormVec, residNormVecOld, Nres,
					turbScalWeightVec, weighted_q_tangent, weighted_lambda_tangent, m, NcontrolEqns, jacMatrix, jacMatrixSTL);
			throw;
		}

		// Fill up the last element of rhs (Nres == q_tangent*(q_guess-q1) + weightLambda*lambda_tangent*(lambda_guess-lambda1) - ds) for mpi_rank == mpi_size-1
		// Note: 1. Even though this part of the code is designed to handle multiple parameters (or lambdas), the actual formulation here is only for one parameters.
		//       2. You don't have to pass weighted_q_tangent instead of q_tangent since it will be weighted in getFlowWeightsForInnerProduct() called by vecDotVecWithWeight().
		if(NcontrolEqns>0) {
			calcNres(Nres, q, rhs, m, NcontrolEqns, q1, q_tangent, lambda_tangent, lambda1, weightLambda, arcLength);
			for(int iParam=0; iParam<NcontrolEqns; ++iParam)
				Nres[iParam] *= weightTangentCond;
			MPI_Barrier(mpi_comm);

			if(mpi_rank==0) {
				for(int iParam=0; iParam<NcontrolEqns; ++iParam) {
					if(fabs(Nres[iParam]/arcLength) > 1.0e4*double(NcontrolEqns*cvora[mpi_size])*MACHINE_EPS) { // Note: Nres must be very close to zero
						cout<<"WARNING in getSteadySolnByNewton(): iterNewton=="<<iterNewton<<", The tangential condition["<<iParam<<"] was NOT satisfied."<<endl
								<<"                                    Tangetial-residual = "<<Nres[iParam]<<" (target arclength="<<arcLength<<")"<<endl;
					}
					if(isnan(Nres[iParam])) {
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

		// Show the total residual calculated with the flow variables
		if(debugLevel>0 && mpi_rank==0)
			printf("           >> Residual of the flow RHS = %.5e", residNormTot);

		// Update residNormTot with Nres
		if(NcontrolEqns>0) {
			residNormTot = updateTotResidualWithNres(residNormTot, Nres, NcontrolEqns, whichNorm);
			if(debugLevel>0 && mpi_rank==0)
				printf(",  Total residual (includ. tangential) = %.5e\n", residNormTot);
		} else {
			if(debugLevel>0 && mpi_rank==0)
				printf("\n");
		}

		// Update residNormOld and residNormTotOld (It will be used for backtracking)
		for(int i=0; i<5+nScal; ++i)
			residNormVecOld[i] = residNormVec[i];
		residNormTotOld = residNormTot;

		// Calculate deltaQ from the previous point (q1) to the initial guess (q)
		double deltaQnorm = 0.0;
		if(NcontrolEqns > 0) { // If NcontrolEqns==0, q1 was not passed (i.e. q1 is NULL if NcontrolEqns==0)
			double mySumPhiSq= 0.0, sumPhiSq;

			double *weightVec = new double [5+nScal];
			for(int icv=0; icv<ncv; ++icv) {
				getFlowWeightsForInnerProduct(weightVec, 1.0, icv, 5+nScal);  // 1.0 --> without volume scaling
				int startIndex = icv*(5+nScal);
				for(int i=0; i<5+nScal; ++i)
					mySumPhiSq += pow(weightVec[i] * (q[i]-q1[i]), 2.0);
			}
			delete [] weightVec;
			if(mpi_rank==mpi_size-1 && NcontrolEqns>0)
				for(int i=0; i<NcontrolEqns; ++i)
					mySumPhiSq += weightLambda * pow(q[ncv*(5+nScal)+i]-q1[ncv*(5+nScal)+i], 2.0);

			MPI_Allreduce(&mySumPhiSq, &sumPhiSq, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

			deltaQnorm = sqrt(sumPhiSq);
		}

		if(trustRegionSize > 8.0*deltaQnorm && NcontrolEqns > 0) {
			if(mpi_rank==0)
				cout<<"WARNING in getSteadySolnByNewton(): TRUST_REGION_SIZE="<<trustRegionSize<<" is larger than 8 * DELTA_Q of the init. guess="<<deltaQnorm<<endl;
			trustRegionSize = 8.0*deltaQnorm;
			if(mpi_rank==0)
				cout<<"                                    Reduce TRUST_REGION_SIZE to "<<trustRegionSize<<endl
				<<endl;
		}

		// Show the result on the screen
		if(isnan(residNormTot)) {
			if(mpi_rank==0) cerr<<"getSteadySolnByNewton(): The residual of the initial field is NaN"<<endl;
			if(debugLevel>1)
				for(int i=0; i<ncv*(5+nScal); ++i)
					if(isnan(rhs[i]))
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
//		if(mpi_rank==0) {
//			FILE *fp;
//			if(firstCall) {
//				fp = fopen(NEWTON_STATUS_FILENAME, "w");
//				fprintf(fp, "STEP,  ITER,      JAC_NNZ,    RESID_INIT,  RELAXATION,   RESID_FINAL,      NRES,        LAMBDA \n");
//			} else
//				fp = fopen(NEWTON_STATUS_FILENAME, "a");
//
//			fprintf(fp, "%4d, %5d,             ,              ,            , %13.6e,          , %13.6e\n", step, 0, residNormTot, lambda[0]);
//			fclose(fp);
//		}

		// Gives a warning message if the residual is too high even though the initial guess comes from a binary
//		if(step - startingStep==2 && initThirdFromQ1) {
//			double startingResidual = getDoubleParam("INIT_THIRD_RESID", "0.0");
//			if(residNormTot > startingResidual)
//				if(mpi_rank==0)
//					cout<<endl
//					<<"WARNING in getSteadySolnByNewton(): The initial residual is greater than INIT_THIRD_RESID="<<startingResidual<<endl
//					<<"                                    Possible reasons - 1. Check if REDUCE_ARCLENTH_THIRD is correct"<<endl
//					<<"                                                       2. Check if LAMBDA_INITIAL_0 and LAMBDA_INITIAL_1 are correct"<<endl
//					<<endl;
//		}

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
		double tolNewton = absTolNewton + relTolNewton*initResidNorm;

		MPI_Barrier(mpi_comm);

		/*
		 * Iterations
		 */
		++iterNewton;
		int succeessiveBacktrack = 0; // Count the number of successive backtracking

		bool done = false;
		int  countConvergedSteps = 0; // In the case that the user wants to run the iterations even after the simulations converges.
//		int  newtonDumpQ1Interval = getIntParam("NEWTON_DUMP_Q1_INTERVAL", "100");
//
//		double relaxation = getDoubleParam("UNDER_RELAXATION", "1.0");
		double relaxation = ikeDotIn.relaxation;

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
//					writePsALCdumpedDataParallel(filename, step, arcLength, lambda1, lambda, NcontrolEqns, q);
				}

//				writeData(step, iterNewton-1);

				freeMemGetSteadySolnByNewton(phi, residNormVec, residNormVecOld, Nres,
						turbScalWeightVec, weighted_q_tangent, weighted_lambda_tangent, m, NcontrolEqns, jacMatrix, jacMatrixSTL);
				throw (PSALC_ERROR_CODE);
			} else {
//				if(iterNewton>startingNewtonIter+1 && newtonDumpQ1Interval>0) { // If the user wants to dump the data at every "newtonDumpQ1Interval" Newton step
//					if((iterNewton-1)%newtonDumpQ1Interval==0) {
//						if(NcontrolEqns>0) {
//							char filename[50];
//							sprintf(filename, "Q1_PT%05d.dumped.%05d.bin", step, iterNewton-1);
//							if(mpi_rank==0)
//								cout<<"   > Dump the current Q vector on "<<filename<<": lambda[0]="<<std::setprecision(15)<<lambda[0]<<endl;
//							writePsALCdumpedDataParallel(filename, step, arcLength, lambda1, lambda, NcontrolEqns, q);
//						}
//
//						writeData(step, iterNewton-1);
//					}
//				}
			}

			/******
			 ** Ramp maximum iterations and residual (For iterative solver)
			 ******/
			if(intervalDecLS > 0)
				if ((iterNewton > startingNewtonIter+1) && (iterNewton >= startDecLS)) {
					if((iterNewton-startDecLS)%intervalDecLS == 0) {
						if(maxIterLS < maxFinalIterLS)
							maxIterLS += incIterLS;
						if(zeroAbsLS > minZeroAbsLS+MACHINE_EPS || zeroRelLS > minZeroRelLS+MACHINE_EPS) {
							zeroAbsLS *= decZeroLS;
							zeroRelLS *= decZeroLS;
						}
						petscSolver2->setTresholds(zeroAbsLS, zeroRelLS, maxIterLS);

						if(debugLevel>1 && mpi_rank==0) {
							printf("RAMP_LSNT_PARAMETERS: MAX_ITER=%d  MIN_ABS_RESID=%.2e  MIN_REL_RESID=%.2e\n", maxIterLS, zeroAbsLS, zeroRelLS);
						}
					}
				}

			/******
			 ** Stabilized Newton's method for ill-conditioned (or almost singular) Jacobians: This will modify the Jacobian matrix
			 ******/
			if(stabilizationMethodType > 0) {
				MPI_Barrier(mpi_comm);

				if(stabilizationMethodType == 1) { // Just apply the diagonal matrix: Often unstable!! (SAS form)
					// Calculate the actual alpha
					double addedValue = max( min(residNormTot, stabilizationAlpha), stabilizationAlphaEps );

					// Matrix
					for(int icv=0; icv<ncv; ++icv) {
						for(int ivar=0; ivar<m; ++ivar) {
							int row = icv*m + ivar;
							int diagIndex = jacMatrixSTL.get_diag_index(row);

							double oldValue = jacMatrixSTL.get_values(diagIndex);

							jacMatrixSTL.set_values(diagIndex, oldValue - addedValue); // Note that the Jacobian matrix is negative (semi-) definite

							++myNonDiagAdd_count; // Just for Statistics
						}
					}
					if(mpi_rank==mpi_size-1) {
						for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
							weighted_lambda_tangent[iEqn] -= addedValue;
						}
					}

					MPI_Barrier(mpi_comm);
				} else if(stabilizationMethodType == 2) { // Apply a scaled diagonal matrix - scaled with local flow variables (Second SAS form)
//					// Calculate the actual alpha
//					double addedValue = max( min(residNormTot, stabilizationAlpha), stabilizationAlphaEps );
//
//					// Calculate the scaling factor
//					for(int icv=0; icv<ncv; ++icv) {
//						int row = icv*m;
//
//						sasScaling[row]   = ADDITIONAL_SCALING_VALUE / (fabs(rho[icv]) + MACHINE_EPS);
//						++row;
//
//						double rhouMag = sqrt(vecDotVec3d(rhou[icv], rhou[icv]) + MACHINE_EPS);
//						for(int i=0; i<3; ++i) {
//							sasScaling[row]   = ADDITIONAL_SCALING_VALUE / (max(0.001, fabs(rhou[icv][i])/rhouMag)*rhouMag);
//							++row;
//						}
//
//						sasScaling[row]   = ADDITIONAL_SCALING_VALUE / (fabs(rhoE[icv]) + MACHINE_EPS);
//						++row;
//
//						for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++) {
//							sasScaling[row]   = ADDITIONAL_SCALING_VALUE / (max(fabs(scalarTranspEqVector[iScal].phi[icv]), 0.001*RefFlowParams.scalar_ref[iScal] + MACHINE_EPS));
//							++row;
//						}
//					}
//					if(mpi_rank==mpi_size-1)
//						for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
//							int row = ncv*m+iEqn;
//							sasScaling[row]   = ADDITIONAL_SCALING_VALUE / q[row];
//						}
//					MPI_Barrier(mpi_comm);
//
//					// Matrix
//					for(int icv=0; icv<ncv; ++icv) {
//						for(int ivar=0; ivar<m; ++ivar) {
//							int row = icv*m + ivar;
//							int diagIndex = jacMatrixSTL.get_diag_index(row);
//
//							double oldValue = jacMatrixSTL.get_values(diagIndex);
//							jacMatrixSTL.set_values(diagIndex, oldValue - addedValue*sasScaling[row]); // Note that the Jacobian matrix is negative (semi-) definite
//
//							++myNonDiagAdd_count; // Just for Statistics
//						}
//					}
//					if(mpi_rank==mpi_size-1) {
//						for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
//							weighted_lambda_tangent[iEqn] -= addedValue*sasScaling[ncv*m+iEqn];
//						}
//					}
//
//					MPI_Barrier(mpi_comm);
				} else if(stabilizationMethodType == 3) { // Apply a scaled diagonal matrix - scaled with rhs: See X.Wu, Applied Mathematics and Computation, vol.189, 2007.
					// Basically finding x that minimize "0.5 * || diag(exp(w*(x-x*))) * f(x) ||_2^2"
					//           or finding the root of "diag(exp(w*(x-x*))) * f(x) = 0"
					// Calculate the actual alpha
					double alpha = max( min(residNormTot, stabilizationAlpha), stabilizationAlphaEps);

					// Matrix
					if(iterNewton>=stabilizationStartingIter) {
						for(int icv=0; icv<ncv; ++icv) {
							for(int ivar=0; ivar<m; ++ivar) {
								int row = icv*m + ivar;
								int diagIndex = jacMatrixSTL.get_diag_index(row);

								double oldValue = jacMatrixSTL.get_values(diagIndex);
								double addVal   = - max(MACHINE_EPS, alpha*fabs(rhs[row])) - stabilizationAlphaEps;
								jacMatrixSTL.set_values(diagIndex, oldValue + addVal); // Note that the Jacobian matrix is negative (semi-) definite

								++myNonDiagAdd_count; // Just for Statistics
							}
						}
						if(mpi_rank==mpi_size-1) {
							for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
								weighted_lambda_tangent[iEqn] -= max(alpha*fabs(rhs[ncv*m+iEqn]), MACHINE_EPS) + stabilizationAlphaEps;

								++myNonDiagAdd_count; // Just for Statistics
							}
						}
					}
				} else if(stabilizationMethodType == 4) { // Apply both scaled diagonal matrix and scaled rhs: See J.L.Hueso, J. Comp. and Applied Math., vol.224, 2009.
					// More generalized version than X.Wu, 2007.
					// Note: 1. We will use the same weight in the A-priori regularization for the exponent (m in the formula)
					//       2. If this version is active, you must deactivate the A-priori regularization in the next part of Newton iteration.
					// Basically finding x that minimize "0.5 * || diag(exp(w*(x-x*))) * f(x)^(1/m) ||_2^2"
					//           or finding the root of "diag(exp(w*m*(x-x*))) * diag(m) * f(x) = 0"
					// Then the Newton algorithm: dx = - ( diag(w*m*f) + A )^(-1) * diag(m) * f

					// Calculate alpha (same as the SAS method , i.e. stabilizationMethodType == 1)
//					double alpha = max( min(residNormTot, stabilizationAlpha), stabilizationAlphaEps); // We will add (subtract) a small number to the diagonal
					double alpha = 1.0/10.0;

					// The coefficients to calculate the exponent (m in the formulation)
					double exponent[5] = {1, 1, 1, 1, 1};
					double weight[5]   = {-1.0e-5, -1.0e-5, -1.0e-5, -1.0e-5, -1.0e-5};

					// If a region should not be modfied based on a priori knowledge (e.g. upstream of scramjet inlet), apply a weight function that is close to zero:
					//     Function: If x<=x_2, w = 10^(-digit) + (1.0-10^(-digit)) * exp(-beta*(x-x2)^2)
					//               Otherwise, w = 1.0
					//     i.e., The function is uni-directional Gaussian
//					double minWeight = pow(10.0, -digitAprioriWeight); // If digitAprioriWeight==3, minWeight=0.001

					// Matrix
					if(iterNewton>=stabilizationStartingIter) {
						for(int icv=0; icv<ncv; ++icv) {
//							double xCoord = x_cv[icv][0];

//							if(xCoord <= x2AprioriWeight) {
//								exponent = minWeight + max( (1-minWeight)*exp(-betaAprioriWeight*pow(xCoord-x2AprioriWeight, 2.0)), MACHINE_EPS );
//								exponent = max(MACHINE_EPS, min(exponent, 1.0)); // Make sure that 0 < exponent <= 1.0
//							}

							for(int ivar=0; ivar<m; ++ivar) {
								int row = icv*m + ivar;
								int diagIndex = jacMatrixSTL.get_diag_index(row);

								double oldValue = jacMatrixSTL.get_values(diagIndex);
								double newValue = oldValue + weight[ivar]*exponent[ivar]*rhs[row];
//								double newValueSign = (newValue>0.0) ? 1.0 : -1.0;

//								jacMatrixSTL.set_values(diagIndex, newValue + newValueSign*alpha);
								jacMatrixSTL.set_values(diagIndex, newValue - alpha);

								++myNonDiagAdd_count; // Just for Statistics
							}
						}
//						if(mpi_rank==mpi_size-1) {
//							for(int iEqn=0; iEqn<NcontrolEqns; ++iEqn) {
//								weighted_lambda_tangent[iEqn] -= max(alpha*fabs(rhs[ncv*m+iEqn]), MACHINE_EPS) + stabilizationAlphaEps;
//
//								++myNonDiagAdd_count; // Just for Statistics
//							}
//						}
					}

					// Rhs
					if(iterNewton>=stabilizationStartingIter) {
						for(int icv=0; icv<ncv; ++icv) {
							int index = icv*m;
							for(int ivar=0; ivar<m; ++ivar) {
								rhs[index+ivar] *= exponent[ivar];

								if(fabs(rhs[index+ivar]) <= MACHINE_EPS)
									rhs[index+ivar] = 0.0;
							}
						}
					}
				}

				MPI_Barrier(mpi_comm);
			}

			/******
			 ** solve the linear (NEWTON) system :
			 ** 1. Normal Newton:    q_{n+1} = q_{n} - J^{-1}*RHS
			 ** 2. Pseudo-arclength: Q_{n+1} = Q_{n} - L^{-1}*RESID
			 **                          Where L = [          J                 dF/dlambda
			 **                                      weighted_q_tangent'  weighted_lambda_tangent ]
			 **                                RESID = [ RHS(q, lambda)
			 **                                          Nres           ]
			 **                          c.f.) Nres = [q_tangent', lambda_tangent']*W*[dq; dlambda]-ds,  where W is the (diagonal) weight matrix
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
				//			solveLinSysNSCoupled2<MatComprsedSTL>(phi, jacMatrixSTL, rhs, nIter, absResid,
				//					zeroAbsLS, zeroRelLS, maxIterLS, nScal, monitorConvergInterval,
				//					NcontrolEqns, ncv_gg, q_tangent, weighted_lambda_tangent, step, iterNewton);
				//			solveLinSysNSCoupled2<MatComprsedSTL>(phi, jacMatrixSTL, rhs, nIter, absResid,
				//					zeroAbsLS, zeroRelLS, maxIterLS, nScal, monitorConvergInterval,
				//					NcontrolEqns, ncv_gg, weighted_q_tangent, weighted_lambda_tangent, step, iterNewton);
				solveLinSysNSCoupled2<MatComprsedSTL>(phi, jacMatrixSTL, rhs, nIter, absResid, kspMonitorHistory,
						useOldSolnAsInitGuess, zeroAbsLS, zeroRelLS, maxIterLS, nScal, monitorConvergInterval,
						NcontrolEqns, ncv_gg, weighted_q_tangent, weighted_lambda_tangent, step, iterNewton);
				// Note that you should pass nScal to this function instead of m(=5+nScal)

				if(firstCall && iterNewton==startingNewtonIter+1 && debugLevel>0 && mpi_rank==0)
					wtimeLS = MPI_Wtime();
				if(firstCall && iterNewton==startingNewtonIter+1 && debugLevel>0 && mpi_rank==0)
					cout<<"   > getSteadySolnByNewton(): Runtime for the LS solver   [sec] = "<<wtimeLS - wtime0<<endl;
			}
//			else if(howToCalcJac == ORDINARY_2D) {
//				// Note: If you want to print out the matrix stored in PETSc on the screen, set the ShowPetscMatrixMatlab variable as true
//				//       ShowPetscMatrixMatlab is defined as a static variable in the PetscSolver2.h file
//				//       e.g. if(iterNewton==1 && NcontrolEqns>0) ShowPetscMatrixMatlab = true;
//				bool useOldSolnAsInitGuess = false; // If true, use the previous solution vector as the initial guess.
//				// Otherwise, apply PC(pre-conditioner) to get the inital guess (the Knoll trick) -- For Bifurcation, this works better.
//				int saveConvergInterval = int(max(monitorConvergInterval/10, 1.01));
//				//			solveLinSysNSCoupled2<MatComprsed>(phi, jacMatrix, rhs, nIter, absResid,
//				//					zeroAbsLS, zeroRelLS, maxIterLS, nScal, monitorConvergInterval,
//				//					NcontrolEqns, ncv_gg, q_tangent, weighted_lambda_tangent, step, iterNewton);
//				solveLinSysNSCoupled2<MatComprsed>(phi, jacMatrix, rhs, nIter, absResid, kspMonitorHistory,
//						useOldSolnAsInitGuess, zeroAbsLS, zeroRelLS, maxIterLS, nScal, monitorConvergInterval,
//						NcontrolEqns, ncv_gg, weighted_q_tangent, weighted_lambda_tangent, step, iterNewton);
//				// Note that you should pass nScal to this function instead of m(=5+nScal)
//			}
			else {
				throw(PSALC_ERROR_CODE);
			}

			if(nIter==maxIterLS && absResid>zeroAbsLS) {
				if(debugLevel>0 || (incIterLS==0 && incIterLS==1.0)) { // TO DO
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

					//				// Calculate the Gershgorin disk
					//				if(howToCalcJac == ROW_1D) {
					//					vector<std::pair<double, double> > gershgorinDisk;
					//					jacMatrixSTL.calcGershgorinDisks(gershgorinDisk);
					//
					////					cout << "    Gershgorin disk (center, radius):";
					////					for(int i=0; i<jacMatrixSTL.get_nRows(); ++i)
					////						cout<<" ("<<gershgorinDisk[i].first<<","<<gershgorinDisk[i].second<<") ";
					////					cout << endl;
					//				}
					//
					//				MPI_Barrier(mpi_comm);
				}
			}

{
int icv = 0;
int m = 5+nScal;
cout<<"            iter = "<<iterNewton<<": phi = "<<phi[icv*m]<<"  "<<phi[icv*m+1]<<"  "<<phi[icv*m+2]<<"  "<<phi[icv*m+3]<<"  "<<phi[icv*m+4]<<endl;
}

			// (Re)Initialize relaxation
			relaxation = ikeDotIn.relaxation;
			MPI_Barrier(mpi_comm);

			/******
			 ** A-priori regularization
			 ******/
//
//			// If a region should not be modfied based on a priori knowledge (e.g. upstream of scramjet inlet), apply a weight function that is close to zero:
//			// Function: If x<=x_2, w = 10^(-digit) + (1.0-10^(-digit)) * exp(-beta*(x-x2)^2)
//			//           Otherwise, w = 1.0
//			// i.e., The function is uni-directional Gaussian
//
//			if(stabilizationMethodType != 4) {
//				if(useAprioriWeightFunction) {
//					double minWeight = pow(10.0, -digitAprioriWeight);
//
//					for(int icv=0; icv<ncv; ++icv) {
//						double xCoord = x_cv[icv][0];
//
//						if(xCoord <= x2AprioriWeight) {
//							double weight = minWeight + max( (1-minWeight)*exp(-betaAprioriWeight*pow(xCoord-x2AprioriWeight, 2.0)), MACHINE_EPS );
//							weight = max(0.0, min(weight, 1.0)); // Make sure that weight never exceed 1.0
//
//							int index = icv*m;
//							for(int i=0; i<m; ++i) {
//								phi[index+i] *= weight;
//								if(fabs(phi[index+i]) <= MACHINE_EPS)
//									phi[index+i] = 0.0;
//							}
//						}
//					}
//					MPI_Barrier(mpi_comm);
//				}
//			}

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

			if(deltaQnorm >= trustRegionSize) {
				if(mpi_rank==0)
					cout<<"WARNING in getSteadySolnByNewton(): After solving the linear system,"<<endl
					<<"                                    deltaQ="<<deltaQnorm<<" is greater than trustRegionSize="<<trustRegionSize<<endl
					<<"                                    Reduce relaxation from "<<relaxation;
				relaxation *= trustRegionSize / deltaQnorm;
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
			relaxation = backtrackForNegativeVals(negativeValCount_CV, negativeValCount_FA, clipParameter, safeParameter, relaxation, kine_index, nScal, q, phi, maxIterFAcheck);

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
						double delta_lambda      = phi[ncv*(5+nScal)+iEqn];
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

			calcResidualsFrom1Drhs(residNormVec, rhs, whichNorm);
			residNormTot = calcTotResidual(residNormVec, whichNorm);

#ifdef USE_TOT_NORM_WITH_LAMBDA
			if(NcontrolEqns>0) {
				// Note: 1. Even though this part of the code is designed to handle multiple parameters (or lambdas), the actual formulation here is only for one parameters.
				//       2. You don't have to pass weighted_q_tangent instead of q_tangent since it will be weighted in getFlowWeightsForInnerProduct() called by vecDotVecWithWeight().
				//			calcNres(Nres, q, rhs, m, NcontrolEqns, q1, q_tangent, lambda_tangent, lambda1, weightLambda, arcLength);
				calcNresForJOE(Nres, m, NcontrolEqns, q1, q_tangent, lambda_tangent, lambda, lambda1, weightLambda, arcLength);
				// Use calcNresForJOE() instead of calcNres() since q has not been updated yet!

				for(int iParam=0; iParam<NcontrolEqns; ++iParam) {
					if(fabs(Nres[iParam]/arcLength) > 1.0e4*double(NcontrolEqns*cvora[mpi_size])*MACHINE_EPS) // Note: Nres must be very close to zero
						if(mpi_rank==0) {
							cout<<"WARNING in getSteadySolnByNewton(): iterNewton=="<<iterNewton<<", The tangential condition["<<iParam<<"] was NOT satisfied."<<endl
									<<"                                    Tangetial-residual = "<<Nres[iParam]<<" (target arclength="<<arcLength<<")"<<endl;
						}

					if(isnan(Nres[iParam])) {
						if(mpi_rank==0)
							cerr<<"ERROR in getSteadySolnByNewton(): iterNewton=="<<iterNewton<<", The residual of the tangential condition becomes NaN for iParam="<<iParam<<endl;
						throw(PSALC_ERROR_CODE);
					}
				}

				for(int iParam=0; iParam<NcontrolEqns; ++iParam)
					Nres[iParam] *= weightTangentCond;
				MPI_Barrier(mpi_comm);
			}

{
int icv = 0;
int m = 5+nScal;
cout<<"            iter = "<<iterNewton<<": rho="<<rho[icv]<<"  rhou="<<rhou[icv][0]<<"  "<<rhou[icv][1]<<"  "<<rhou[icv][2]<<"  rhoE="<<rhoE[icv]<<endl;
cout<<"                       rhs="<<rhs[icv*m]<<"  "<<rhs[icv*m+1]<<"  "<<rhs[icv*m+2]<<"  "<<rhs[icv*m+3]<<"  "<<rhs[icv*m+4]<<endl;
}

			if(debugLevel>0 && mpi_rank==0)
				printf("           >> Residual of the flow RHS = %.5e", residNormTot);

			if(NcontrolEqns>0) {
				residNormTot = updateTotResidualWithNres(residNormTot, Nres, NcontrolEqns, whichNorm);

				if(debugLevel>0 && mpi_rank==0)
					printf(",  Total residual (includ. tangential) = %.5e\n", residNormTot);
			} else {
				if(debugLevel>0 && mpi_rank==0)
					printf("\n");
			}
#else
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
			if(backtrackMaxIter > 0) {
				bool needBacktracking;
//				if(useWatchdogBacktrack)
//					needBacktracking = checkBacktrackFromResidIncUsingWatchdog(residNormTotOld, residNormTot, min(relaxation, backtrackRelax_UpperBound), gammaWatchdog, NcontrolEqns, 5+nScal, rhs, phi);
//				else
					needBacktracking = checkBacktrackFromResidInc             (residNormTotOld, residNormTot, min(relaxation, backtrackRelax_UpperBound), backtrackBarrierCoeff);

				if(countNegative>0 || needBacktracking) {
					++succeessiveBacktrack;

					// Check if backtracking is required: If the conditions are satisfied, backtracking will be omitted
					if(skipBT_freq>0)
						forgiveBacktrack = ( (skipBT_firstITer && iterNewton==startingNewtonIter+1) || succeessiveBacktrack%skipBT_freq==0 );
					else
						forgiveBacktrack = ( skipBT_firstITer && iterNewton==startingNewtonIter+1 );

					// Check if the user want to skip backtraking for first few times
					if(forgiveBacktrack && countNegative==0) { // Note: If countNegative>0, backtracking will be active anyway
						relaxation = 0.2 * ( backtrackRelax_LowerBound + min(relaxation, backtrackRelax_UpperBound) );

						updateFlow_primVars(q, phi, -relaxation, nScal);
						if(NcontrolEqns>0)
							updateLambda(q, phi, -relaxation, NcontrolEqns, nScal);

						if(debugLevel > 0 && mpi_rank == 0)
							cout<<"                  Forgive BACKTRACKING to escape from the current point: Instead, set relaxation="<<relaxation<<" (lambda="<<lambda[0]<<")"<<endl;
					} else {
						MPI_Barrier(mpi_comm);

						// Show the reason to use the backtracking algorithm on the screen if the debug level is high
						if(debugLevel > 0) {
							if(mpi_rank==0) {
								cout<<"                  Calling BACKTRACKING ";
								if(debugLevel > 1) {
									cout<<"due to";
									if(countNegative>0)
										printf("  negative vars=%d", countNegative);
									if(needBacktracking)
										printf("  residual=%.8e", residNormTot);
								}
								cout<<endl;
							}
						}

						// Launch the backtracking algorithm -------
						// First, calculate the relaxation size by calling backtrackWithJOE_calcRelaxAndRHS()
						relaxation = backtrackWithJOE_calcRelaxAndRHS(rhs, q, phi, btNotConverged,
								backtrackMaxIter, relaxation, backtrackRelax_LowerBound, backtrackRelax_UpperBound,
								residNormVecOld, residNormVec, whichNorm, backtrackBarrierCoeff,
								NcontrolEqns, q_tangent, lambda_tangent, weightLambda, weightTangentCond, q1, Nres, lambda1, arcLength);
						// Since the reduction algorithm based on delta_lambda is heuristics, use "relaxBeforeDlambda" instead of "relaxation"

						calcResidualsFrom1Drhs(residNormVec, rhs, whichNorm);
						residNormTot = calcTotResidual(residNormVec, whichNorm);

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
							for(int iParam=0; iParam<NcontrolEqns; ++iParam)
								Nres[iParam] *= weightTangentCond;
							MPI_Barrier(mpi_comm);
							residNormTot = updateTotResidualWithNres(residNormTot, Nres, NcontrolEqns, whichNorm);
						}
#endif
						//					// If backtracking fails, launch simple reducing algorithm
						//					if(residNormTot > residNormTotFresh) { // If the result got worse than before backtracking
						//						relaxation = relaxBeforeBacktrack;
						//
						//						updateFlow_primVars(q, phi, -relaxation, nScal);
						//						if(NcontrolEqns>0) {
						//#ifdef UPDATE_LAMBDA_IN_BACKTRACK
						//							updateLambda(q, phi, -relaxation, NcontrolEqns, nScal);
						//#endif
						//						}
						//
						//#ifdef USE_TOT_NORM_WITH_LAMBDA
						//						if(NcontrolEqns>0) {
						//							// Finally, update the total residual using the tangential condition
						//							calcNresForJOE(Nres, m, NcontrolEqns, q1, q_tangent, lambda_tangent, lambda, lambda1, weightLambda, arcLength);
						//									// Use calcNresForJOE() instead of calcNres() since q has not been updated yet!
						//							for(int iParam=0; iParam<NcontrolEqns; ++iParam)
						//								Nres[iParam] *= weightTangentCond;
						//							MPI_Barrier(mpi_comm);
						//							residNormTot = updateTotResidualWithNres(residNormTot, Nres, NcontrolEqns, whichNorm);
						//						}
						//#endif
						//					}
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
//			bool solnModified = temporalHookNewton(q, lambda, phi, relaxation);
//			if(solnModified) {
//				updateCvDataG1G2(rho, REPLACE_DATA);
//				updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
//				updateCvDataG1G2(rhoE, REPLACE_DATA);
//
//				for (int iScal = 0; iScal < nScal; iScal++)
//					updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);
//
//				calcStateVariables();
//			}

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
					}
//					else if(howToCalcJac == ORDINARY_2D) {
//						jacMatrix.clear();
//						countNegative = calcJacobianAD(jacMatrix, rhs, debugLevel, NcontrolEqns);
//						myNnzJac = jacMatrix.get_nnz();
//					}
					else { // ERROR
						freeMemGetSteadySolnByNewton(phi, residNormVec, residNormVecOld, Nres,
								turbScalWeightVec, weighted_q_tangent, weighted_lambda_tangent, m, NcontrolEqns, jacMatrix, jacMatrixSTL);
						throw(PSALC_ERROR_CODE);
					}
					MPI_Allreduce(&myNnzJac, &nnzJac, 1, MPI_INT, MPI_SUM, mpi_comm);
				}
				catch(int e) { // Catch and re-throw
					freeMemGetSteadySolnByNewton(phi, residNormVec, residNormVecOld, Nres,
							turbScalWeightVec, weighted_q_tangent, weighted_lambda_tangent, m, NcontrolEqns, jacMatrix, jacMatrixSTL);
					throw(e);
				}
				catch(...) { // Catch and re-throw
					freeMemGetSteadySolnByNewton(phi, residNormVec, residNormVecOld, Nres,
							turbScalWeightVec, weighted_q_tangent, weighted_lambda_tangent, m, NcontrolEqns, jacMatrix, jacMatrixSTL);
					throw;
				}

				if(countNegative>0) {
					freeMemGetSteadySolnByNewton(phi, residNormVec, residNormVecOld, Nres,
							turbScalWeightVec, weighted_q_tangent, weighted_lambda_tangent, m, NcontrolEqns, jacMatrix, jacMatrixSTL);
					if(mpi_rank==0) cout<<"ERROR! getSteadySolnByNewton(): After negative vars = "<<countNegative<<endl;
					throw(PSALC_ERROR_CODE);
				}
			} else { // Calculate rhs, residNorm, and residNormTot by using normal JOE
				bool useBarrier = true;
				int countNegative = calcRhsWithBarrier(rhs, useBarrier);
				calcResidualsFrom1Drhs(residNormVec, rhs, whichNorm);
				residNormTot = calcTotResidual(residNormVec, whichNorm);
			}

			calcResidualsFrom1Drhs(residNormVec, rhs, whichNorm);

			/******
			 ** fill up the last element of rhs (Nres == q_tangent*(q_guess-q1) + weightLambda*lambda_tangent*(lambda_guess-lambda1) - ds) for mpi_rank == mpi_size-1
			 ** Note: Even though this part of the code is designed to handle multiple parameters (or lambdas), the actual formulation here is implemented for one parameter now.
			 ******/
			if(NcontrolEqns>0) {
				calcNres(Nres, q, rhs, m, NcontrolEqns, q1, q_tangent, lambda_tangent, lambda1, weightLambda, arcLength);

				if(mpi_rank==0) {
					for(int iParam=0; iParam<NcontrolEqns; ++iParam)
						if(fabs(Nres[iParam]/arcLength) > 1.0e4*double(NcontrolEqns*cvora[mpi_size])*MACHINE_EPS)
							cout<<"WARNING in getSteadySolnByNewton(): iterNewton=="<<iterNewton<<", The tangential condition["<<iParam<<"] was NOT satisfied."<<endl
							<<"                                    Tangetial-residual = "<<Nres[iParam]<<" (target arclength="<<arcLength<<")"<<endl;
				}

				for(int iParam=0; iParam<NcontrolEqns; ++iParam)
					Nres[iParam] *= weightTangentCond;

				if(mpi_rank == mpi_size-1) {
					for(int iParam=0; iParam<NcontrolEqns; ++iParam)
						rhs[ncv*m+iParam] = Nres[iParam];
				}
			}

			if(debugLevel>0 && mpi_rank==0)
				printf("           >> Residual of the flow field RHS = %.5e", residNormTot);

			if(NcontrolEqns>0) {
				residNormTot = updateTotResidualWithNres(residNormTot, Nres, NcontrolEqns, whichNorm);

				if(debugLevel>0 && mpi_rank==0)
					printf(",  Total residual including the tangential condition = %.5e\n", residNormTot);
			} else {
				if(debugLevel>0 && mpi_rank==0)
					printf("\n");
			}

			/******
			 ** Update residNormVecOld and residNormTotOld:
			 ** Note: "residNormTotOld" is used only to compare residuals from RHS calculation without the tangential condition,
			 **       thus it should be updated before the tangential condition.
			 ******/
			if(residNormTot < residNormTotOld) {
				for(int i=0; i<5+nScal; ++i)
					residNormVecOld[i] = residNormVec[i];
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
{
int icv = 0;
int m = 5+nScal;
cout<<"            iter = "<<iterNewton<<": rho="<<rho[icv]<<"  rhou="<<rhou[icv][0]<<"  "<<rhou[icv][1]<<"  "<<rhou[icv][2]<<"  rhoE="<<rhoE[icv]<<endl;
}
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
					printf("     %15.8e\n\n",lambda[0]);
				} else {
					printf("\n\n");
				}
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
						if(mpi_rank==0 && countConvergedSteps==1)
							cout<<"Newton-solver converged but keep running the simulation since MORE_STEPS_BELOW_CONV="<<moreNTstepsBelowConv_moreSteps<<": DELTA_Q="<<moreNTstepsBelowConv_deltaQ<<endl;

						if(deltaQnorm <= moreNTstepsBelowConv_deltaQ)
							done = true;
					}
				}

				if(residNormTot <= moreNTstepsBelowConv_resid)
					done=true;

				if(iterNewton > startingNewtonIter+maxIterNewton)
					done = true;
			}
			////IKJ
			//writeData(step, iterNewton);
		}

		if(debugLevel>0 && mpi_rank == 0)
			printf("Newton solver converged after %5dth outer iter: residual = %12.5e \n", iterNewton-1, residNormTot);

//		// Write the Jacobian matrix on a file for post-processing(eigen-decomposition)
//		if(writeJacOnFile) {
//			//		if(modifiedNewtonMethod == NO_MODIFIED_METHOD) {
//			//			assert(jacMatrix.empty());
//			//			calcJacobianADOLC(tagNum, jacMatrix);
//			//		} else {
//			//			if(!useExactJacTempo) { // if Jacobian is NOT exact
//			//				jacMatrix.clear();
//			//				calcJacobianADOLC(tagNum, jacMatrix);
//			//			}
//			//		}
//
//			if(stabilizationMethodType != 0 || !useExactJacTempo) { // i.e., If the Jacobian matrix was modified in order not to have a singula matrix
//				//       or an approximate Jacobian is used, ...
//				jacMatrix.clear();
//
//				if(howToCalcJac == ROW_1D) {
//					jacMatrixSTL.clear();
//					countNegative = calcJacobian1DAD(jacMatrixSTL, rhs, debugLevel, NcontrolEqns);
//				} else if(howToCalcJac == ORDINARY_2D) {
//					jacMatrix.clear();
//					countNegative = calcJacobianAD(jacMatrix, rhs, debugLevel, NcontrolEqns);
//				}
//			}
//
//			stringstream ss;
//			ss<<"jacMatrix"<<step<<".bin";
//			string filename = ss.str();
//
//			if(jacMatrixSTL.empty()) {
//				if(jacMatrix.empty()) { // This is an error!! Both jacMatrixSTL and jacMatrix are empty!
//					if(mpi_rank==0) {
//						cout<<"ERROR! IkeWithPsALC_AD::getSteadySolnByNewton(): Both jacMatrixSTL and jacMatrix are empty!"<<endl;
//						cout<<"                                                 Cannot write the Jacobian matrix on a file"<<endl;
//					}
//					throw(PSALC_ERROR_CODE);
//				}
//				//			writeMatrixMfile<MatComprsed>(filename, jacMatrix, nScal, cvora, nbocv_v_global);
//				writeMatrixBinaryParallel<MatComprsed>(filename, jacMatrix, nScal, cvora, nbocv_v_global);
//			} else {
//				//			writeMatrixMfile<MatComprsedSTL>(filename, jacMatrixSTL, nScal, cvora, nbocv_v_global);
//				writeMatrixBinaryParallel<MatComprsedSTL>(filename, jacMatrixSTL, nScal, cvora, nbocv_v_global);
//			}
//		}

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
		freeMemGetSteadySolnByNewton(phi, residNormVec, residNormVecOld, Nres,
				turbScalWeightVec, weighted_q_tangent, weighted_lambda_tangent, m, NcontrolEqns, jacMatrix, jacMatrixSTL);

		//IKJ
		if(sasScaling != NULL) {
			delete [] sasScaling; 		sasScaling = NULL;
		}
		if(sasScaledPHI != NULL) {
			delete [] sasScaledPHI; 	sasScaledPHI = NULL;
		}

		firstCall = false;

		MPI_Barrier(mpi_comm);
	}
};

/***********
 * MAIN
 ***********/
int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	initMpiStuff();

	if(mpi_rank == 0)
		cout<<endl
		    <<"Welcome to standaloneNewton.cpp ... "<<endl
		    <<"MPI initialization check ... "<<"mpi_size = "<<mpi_size<<endl
		    <<endl;

	PsalcStandAlone psalcStandAlone;

	int ncv = psalcStandAlone.ncv;
	int nVars = 5 + psalcStandAlone.nScal;

	double *q1  = new double[ncv*nVars];
	double *rhs = new double[ncv*nVars];
	int maxIterNewton   = 20; //psalcStandAlone.ikeDotIn.maxIterLS;
	double absTolNewton = 1.0e-10; //psalcStandAlone.ikeDotIn.zeroAbsLS;
	double relTolNewton = 1.0e-12; //psalcStandAlone.ikeDotIn.zeroRelLS;

	psalcStandAlone.initialHookNewton();
	for(int icv=0; icv<ncv; ++icv) {
		q1[icv*nVars]   = psalcStandAlone.rho[icv];
		q1[icv*nVars+1] = psalcStandAlone.rhou[icv][0];
		q1[icv*nVars+2] = psalcStandAlone.rhou[icv][1];
		q1[icv*nVars+3] = psalcStandAlone.rhou[icv][2];
		q1[icv*nVars+4] = psalcStandAlone.rhoE[icv];
	}
	try {
		psalcStandAlone.getSteadySolnByNewton(q1, rhs, maxIterNewton, absTolNewton, relTolNewton);
	}
	catch (int e) {
		cout<<"catched: "<<e<<endl;
	}
	catch (...) {

	}

	delete [] q1;
	delete [] rhs;

	MPI_Barrier(MPI_COMM_WORLD);
}
