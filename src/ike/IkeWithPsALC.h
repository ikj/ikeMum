/*
 * IkeWithPsALC.h
 *
 *  Created on: Nov 10, 2012
 *      Author: ikj
 */

#ifndef IKEWITHPSALC_H_
#define IKEWITHPSALC_H_

#include <memory> // Use shared_pointer

#include "ADJOINT_FILES/JoeWithModelsAD.h"
#include <adolc_sparse.h>
#include <set>
#include "MatComprsed.h"
#include "PetscSolver2.h"

#include "IkeUtilsAD.h"
#include "IkeWithModels.h"

#include "ADvar.h"

#ifndef WITH_PETSC
#define WITH_PETSC
#endif


/* Internal Buffer sizes */
// <adolc/usrparms.h>
//#define OBUFSIZE 524288  // default=65536
//#define LBUFSIZE 524288  // default=65536
//#define VBUFSIZE 524288  // default=65536
///* Internal buffer size for temporary Taylor store */
//#define TBUFSIZE 524288  // default=65536
///* Maximal number of Taylor stacks */
//#define TBUFNUM  256     // default=32

enum HOW_TO_CALC_JAC {ROW_1D, ORDINARY_2D};
enum MODIFIED_NEWTON {BASIC, SHAMANSKII, MODIFIED_SHAMANSKII};

#ifndef ABSURDLY_BIG_NUMBER
	#define ABSURDLY_BIG_NUMBER 2.22e22
#endif

#define PSALC_ERROR_CODE -2
#define IKE_BACKTRACKING_CODE -22

#ifndef MACHINE_EPS
#define MACHINE_EPS 2.0e-16 // Note: machine epsilon for double precision = 1.11e-16
#endif

#define JAC1D_STATUS_FILENAME  "Jac1D_status.dat"
#define BIFUR_SUMMARY_FILENAME "Bifur_summary.csv"
#define NEWTON_STATUS_FILENAME "Newton_converg.csv"
#define LS_CONVERGENCE_FILENAME "LinearSystem_converg.csv"

#ifndef EOF_ERROR_CHECK_CODE
#define EOF_ERROR_CHECK_CODE 503
#endif

#define GRID_TOL 1.0e-9

#define ARCLENGTH_WARNING_EPS 1.0e-5 // 0.001% error
#define ARCLENGTH_ERROR_EPS 1.0e-3   // 0.1%   error

#define LAMBDA_COMPATIBILITY_EPS 1.0e-4 // 0.01% error

#define USE_LOCAL_VALUE_SCALING // The Jacobian of mass conservation has zero diagonal entry (in the analytic expression)
                                // and small off-diagonal entries, which generates almost zero eigenvalue according to the Gershgorin disk theorem.
                                // However, the Jacobian of energy conservation has a huge diagonal entry that can produce a very large eigenvalue.
                                // Thus, if the equations are not scaled properly, they can produce a huge condition number.
                                // Note: USE_DT_OVER_VOL_SCALING can be defined in IkeWithModels.h
//#define USE_REF_VALUE_SCALING // Similar to USE_LOCAL_VALUE_SCALING, but use REFERENCE values instead of local rho, rhou, rhoE, etc.
#define ADDITIONAL_SCALING_VALUE 1.0 // If the scaling is applied, the flow residual drops and sometimes it becomes too small compared to the tangential equation.
                                      // Then, the total residual (flow residual + tangential residual) is domiated by the tangential residual. Thus, we need to scale it up.

#define USE_TOT_NORM_WITH_LAMBDA // Sometimes, it seems to be natural not to use the tangential condition when calculating the total norm.
                                 // If you DON'T define this variable, then the tangential condition will not be considered in some parts of
                                 // the getSteadySolnByNewton() method and the backtrackWithJOE_calcRelaxAndRHS() methods.

#define USE_VOLUME_WEIGHTED_INNER_PROD // Use weighted inner product with CV volumes,
                                       // i.e., || [q; lambda] ||_{W} = [q^{T}, lambda^{T}] * [ cv_volume      0         * [  q
                                       //                                                           0     weightLambda ]    lambda]
#ifdef USE_VOLUME_WEIGHTED_INNER_PROD
#define MIN_CV_VOL_CUTOFF 1.0e-12
#endif

#define USE_REF_WEIGHTED_INNER_PROD // Similar idea as using USE_VOLUME_WEIGHTED_INNER_PROD,
                                    // but each variables are normalized by its reference value (e.g. rho / RHO_REF),
                                    // i.e., || [q; lambda] ||_{W} = [q^{T}, lambda^{T}] * [ 1/REF_VAL      0         * [  q
                                    //                                                           0     weightLambda ]    lambda]

#define UPDATE_LAMBDA_IN_BACKTRACK

#define RELAXATION_EPS 7.777e-6


struct wallTimeJacCalc {
	double wallTimeRHSinit;
	double wallTimeRHScalc;
	double wallTimeRHSfinal;
	double wallTimeJac;
	double wallTimeCleanup;
};


/*
 * Bug report:
 *   In JoeWithModelsAD.cpp
 *     1. cannot detect "kine_index" in the calcViscousFluxNS_AD() method:
 *        Original: if (scalarTranspEqVector[iScal].getName() == "kine")
 *        Modified: if (strcmp(scalarTranspEqVector[iScal].getName(),"kine")==0)
 *
 *     2. "kine_fa" is wrong in the calcViscousFluxNS_AD() method:
 *        Original: if (kine_index > -1)
 *			        {
 *	                   kine0 = scalarTranspEqVector_AD[kine_index].phi[icv0];
 *                     kine1 = scalarTranspEqVector_AD[kine_index].phi[icv1];
 *                     REALQS kine_fa = w1*kine0 + w0*kine1;
 *                  }
 *        Modified: if (kine_index > -1)
 *			        {
 *	                   kine0 = scalarTranspEqVector_AD[kine_index].phi[icv0];
 *                     kine1 = scalarTranspEqVector_AD[kine_index].phi[icv1];
 *                     kine_fa = w1*kine0 + w0*kine1;
 *                  }
 */

/*
 * Modifications in JOE and ADKOINT_FILES:
 *   UgpWithCvCompFlowAD.h:
 *     1. insert the following two lines at the beginning of UgpWithCvCompFlowAD.h
 *       enum BOUNDARY_TYPE {WRONG, INTERNAL, HOOK, CBC, CBC_SUBSONIC_INLET, CBC_SUBSONIC_OUTLET, SYMMETRY, NEUMANN, WALL};
 *       enum BOUNDARY_TYPE_SCALAR {INTERNAL_SCALAR, HOOK_SCALAR, DIRICHLET_SCALAR, FLUX_SCALAR, OTHER_SCALAR};
 *     2. insert the following line in the UgpWithCvCompFlowAD class as public
 *       #include "../../IkeUgpWithCvCompFlow.h"
 *     3. make initialize_adoubles() and destroy_adoubles() virtual
 *     4. add UgpWithCvCompFlow_AD_init() in UgpWithCvCompFlow_AD::init() so that we can initialize some variables later
 *     		void init() {
 *     			....
 *     			UgpWithCvCompFlow_AD_init();
 *     		}
 *     5. add UgpWithCvCompFlow_AD_clear() in the destructor of UgpWithCvCompFlow_AD
 *     		virtual ~UgpWithCvCompFlow_AD() {
 *     			UgpWithCvCompFlow_AD_clear();
 *     		}
 *   JoeWithModels class:
 *     1. make run() virtual
 *
 *   UgpWithCvCompFlowAD.h:
 *     1. Comment out the following line in the initialize_adoubles() method
 *         if(mpi_rank==0) cout<<"Allocated Scalar "<<scalarTranspEqVector_AD[i].name<<endl
 *
 *   And there are more...
 */

/* Note: If you want to use the SLEPc package, define USE_SLEPC_WITH_PETSC in PetscSolver2.h
 *       (SLEPc is an extension of PETSc which can be used to solve large scale sparse eigenvalue problems on parallel computers.
 *       It can be used for either standard or generalized eigen-problems, a partial SVD of a large, sparse, rectangular matrix,
 *       and quadratic eigenvalue problems: http://www.grycap.upv.es/slepc )
 */

/* Note: Numbering convention in Joe (the version of two-layer-ghosts)
 *   1. Faces
 *     0 ~ nfa_b-1 : boundary faces
 *     nfa_b ~ nfa_bpi-1 :
 *     nfa_bpi ~ nfa-1 :
 *     nfa ~ nfa_b2g : additional boundary faces between g1:f2 cvs
 *     nfa_b2g ~ nfa_b2g-1: faces between g1:g1 cvs
 *     nfa_b2g-1 ~ nfa_b2gg-1: faces between g1:g2 cvs
 *   2. CVs
 *     0 ~ ncv-1: internal cells
 *     ncv ~ ncv_g-1 : includes level-1 ghosts
 *     ncv_g ~ ncv_gg-1 : includes level-2 ghosts (face-based nbrs of level-1 ghosts)
 *     ncv_gg ~ ncv_ggf-1 : includes fakes owned by nfa_b faces: i.e. ncv_ggf-ncv_gg == nfa_b
 *     ncv_ggf ~ ncv_ggff-1 : includes fakes that are nbrs of level-1 ghosts
 */

/* Note: Multiple inheritance occurs here. The problem of ambiguous hierarchy composition is avoided by using virtual inheritance
 *       (In the diagram below, doubled lines are virtual public, single lines are public)
 *
 *                        UgpWithCvCompFlow
 *                  //           ||           \\
 *     JoeWithModels    UgpWithCvCompFlow_AD    RansTurbKOm
 *              \\       //              \\       /?
 *         JoeWithModels_AD              RansTurbKOm_AD
 *                     ||                  |
 *         IkeWithModels_AD              IkeRansTurbKOm_AD
 *                      |                  |
 *         IkeWithPsALC_AD                 |
 *                       \                /
 *                      endUserClass(ike.cpp)
 */

/*
 * Class: IkeWithPsALC_AD
 */
class IkeWithPsALC_AD: public IkeWithModels_AD //, virtual public JoeWithModels, virtual public UgpWithCvCompFlow_AD
{
public:
	/*
	 * constructor, pass ParamMap
	 */
	IkeWithPsALC_AD(ParamMap &p) : IkeWithModels_AD(p), JoeWithModels_AD(p), JoeWithModels(p), UgpWithCvCompFlow_AD(p) {
		if(mpi_rank==0)
			cout<<"IkeWithPsALC_AD()"<<endl;
		this->init();
	}
	/*
	 * constructor, pass name of joe's input file
	 */
	IkeWithPsALC_AD(char *name) : IkeWithModels_AD(name), JoeWithModels_AD(name), JoeWithModels(name), UgpWithCvCompFlow_AD(name) {
		if(mpi_rank==0)
			cout<<"IkeWithPsALC_AD()"<<endl;
		this->init();
	}

	/*
	 * destructor
	 */
	virtual ~IkeWithPsALC_AD() {
		if(mpi_rank==0)
			cout<<"~IkeWithPsALC_AD()"<<endl;
		this->clear();
	}

private:
	/*
	 * Method: init
	 * ------------
	 *
	 */
	void init();

	/*
	 * Method: clear
	 * -------------
	 *
	 */
	void clear();

public:
	/****************************
	 * PSEUDO-ARCLENGTH FUNCTIONS
	 ****************************/
	/*
	 * Method: run
	 * -----------
	 *
	 */
	virtual void run();

	/*
	 * Method: runPsALC
	 * ----------------
	 * Argument: howToCalcJac = How to calculate the Jacobian matrix: ROW_1D, ORDINARY_2D -- default = ORDINARY_2D
	 */
	void runPsALC();

	/*
	 * Method: freePsALCMemory
	 * -----------------------
	 * Clear the arrays used in the runPsALC class.
	 * This method is called not only at the end of the runPsALC class but also whenever an exception occurs.
	 */
	void freePsALCMemory(double **p_q0, double **p_q1, double **p_Qguess, double **p_lambda0, double **p_lambda1,
			double **p_lambdaInit0, double **p_lambdaInit1, double **p_lambda_tangent, const int nVars);

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
	 *     normSqRatioInArclength = || q1-q0 ||^2 / || lambdaInit1-lambdaInit0 ||^2
	 */
	void initFirstTwoPts(double* q0, double* q1, double *rhs, double *lambdaInit0, double *lambdaInit1,
			const int NcontrolEqns, double &arclength, double *normSqRatioInArclength, const double weightLambda, string &QoIfilename);

	/*
	 * Method: singlePsALCstep()
	 * -----------------------------------
	 * Perform one continuation on the solution curve
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
	 *            arcLength = arc-length
	 *            lambda0 =
	 *            lambda1 =
	 *            weightLambda =
	 *            writeJacOnFile =
	 *            factorReduceArclength = The user can reduce the length of the tangential vectors by setting this argument less than 1.0
     *                                    (The arclength itself must have been reduced in the function calling this method before calling this.
     *                                     However, the tangential vectors calculated in this method can be only reduced here)
	 */
	void singlePsALCstep(double *qVec, const double *q0, const double *q1, double *rhs, double **q_tangent, double *lambda_tangent,
			const double arcLength, const double* lambda0, const double* lambda1, const double weightLambda, const int writeJacOnFile,
			const double factorReduceTangential = 1.0);

	/*
	 * Method: arclengthControl
	 * ------------------------
	 * Calculate the ratio for the arclength to be reduced or increased
	 * Note: Increasing the arclength is more conservative than reducing it.
	 */
	double arclengthControl(const double residNormInit, const double residNormInit_opt, const int niterNewton, const int niterNewton_opt);

	/*
	 * Method: getSteadySolnByNewton
	 * -----------------------------
	 *
	 * Note: tolerance = absTolNewton + relTolNewton * norm(rhs_initial)   -- C.T.Kelley, SIAM 2003, Chap.1.5
	 */
	void getSteadySolnByNewton(double* q, double* rhs, const int maxIterNewton, const double absTolNewton, const double relTolNewton,
			const bool writeJacOnFile=false,
			const int NcontrolEqns=0, const double *q1=NULL, double** q_tangent=NULL,
			double* lambda_tangent=NULL, const double *lambda0=NULL, const double *lambda1=NULL, const double weightLambda=0.0, const double arclength=0.0);

	/*
	 * Method: freeMemGetSteadySolnByNewton
	 * ---------------------------------------
	 * Clear memory allocated for the getSteadySolnByNewton() function
	 * Note: The arrays will NOT be re-initialized by NULL
	 */
	void freeMemGetSteadySolnByNewton(double *phi, double* residNorm, double* residNormOld, double *Nres,
			MatComprsed& jacMatrix, MatComprsedSTL& jacMatrixSTL);

	void freeMemGetSteadySolnByNewton(double *phi, double* residNorm, double* residNormOld, double *Nres,
			double *turbScalWeightVec, double** weighted_q_tangent, double* weighted_lambda_tangent, const int nVars, const int NcontrolEqns,
			MatComprsed& jacMatrix, MatComprsedSTL& jacMatrixSTL);

	/*
	 * Method: getSteadySolnByBackwardEuler
	 * ------------------------------------
	 * Original code = runBackwardEuler() in JoeWithModels.cpp
	 * Note: this method cannot update lambda
	 *
	 * Return:
	 *   By value    : isConverged
	 *   By reference: dq, dScal, rhs, rhsScal, A, Ascal, Residual
	 */
	bool getSteadySolnByBackwardEuler(double (*dq)[5], double **dScal, double (*rhs)[5], double **rhsScal,
			double (*A)[5][5], double ***AScal, double &initTotResidual, double *Residual,
			const int maxStepEuler, const double absTolEuler, const double relTolEuler,
			const int maxIterLS, const double zeroAbsLS, const double zeroRelLS,
			double &cfl, const int startIncCFL, const int intervalIncCFL, const double incCFL, const double maxCFL,
			const int check_interval = 1);

	/*
	 * Method: calcJacobian1DAD
	 * ------------------------
	 * calculate Rhs and get the Jacobian matrix using AD using the so-called "1D style"
	 * Return: If backtracking line search is required (due to negative density, etc.), return true.
	 * Original code of Rhs part = JoeWithModels_AD::calcResidualDerivative(double***, double***, int)
	 */
	bool calcJacobian1DAD(MatComprsedSTL &jacMatrixSTL, double *rhsSingleArray, const int debugLevel, const int NcontrolEqns=0);

//	/*
//	 * Method: calcJacobian1DAD_new
//	 * ----------------------------
//	 * calculate Rhs and get the Jacobian matrix using AD using the so-called "1D style"
//	 * Return: If backtracking line search is required (due to negative density, etc.), return true.
//	 * Original code of Rhs part = JoeWithModels_AD::calcResidualDerivative(double***, double***, int)
//	 */
//	bool calcJacobian1DAD_new(MatComprsedSTL &jacMatrixSTL, double *rhsSingleArray, bool printStats, const int NcontrolEqns);

	/*
	 * Method: calcJacobianAD
	 * ----------------------
	 * calculate Rhs and get the Jacobian matrix using AD
	 * Arguments: jacMatrix      = The Jacobian matrix (NOTE: the size of the matrix is (ncv*nVars)*(ncv_gg*nVars), where nVars=5+nScal)
	 *            rhsSingleArray = R.H.S. of the N-S system with scalar equations (NOTE: The size is ncv_gg*(5+nScal) because using Dr. Duraisamy's adjoint code)
	 *            printStats     = If yes, show the Jacobian statistics at the first cpu core on the screen (for a debugging purpose)
	 * Return: If backtracking line search is required (due to negative density, etc.), return true.
	 * Original code of Rhs part = JoeWithModels_AD::calcResidualDerivative(double***, double***, int)
	 */
	bool calcJacobianAD(MatComprsed &jacMatrix, double *rhsSingleArray, const int debugLevel, const int NcontrolEqns=0);

	/*
	 * Method: calcJacobian1DAD_calcRhs
	 * --------------------------------
	 * Calculate RHS for Jacobian calculation
	 * Return:
	 *   rhsSingleArray
	 *   myCountReducedOrder
	 *   wall-clock times for RHS calculation
	 */
	void calcJacobian1DAD_calcRhs(double *rhsSingleArray, int &myCountReducedOrder, wallTimeJacCalc &myWTimeJacCalc,
			const int icv, const int tag, const int nScal, const int debugLevel, const int NcontrolEqns,
			const bool firstCall, bool &firstCallScalarTurb);

	/*
	 * Method: probeInCalcJacobian1DAD
	 * -------------------------------
	 * This method will be called in the end of calcJacobian1DAD() -- just before the memory-releasing step.
	 * If the user wants to probe a quantity (usually something related to the RHS calculation) at every Newton-Raphson iteration,
	 * she/he may want to call this method.
	 */
	virtual void probeInCalcJacobian1DAD(MatComprsedSTL &jacMatrixSTL, double *rhsSingleArray, const int debugLevel, const int NcontrolEqns);

	/*
	 * Method: solveLinSysNSCoupled2
	 * -----------------------------
	 * solve coupled linear system for Navier-Stokes equations and scalars together
	 * original code: solveCoupledLinSysNSCoupled() in UgpWithCvCompFlow.h
	 *
	 * Matrix = [ ---A---
	 *            vecC  d ]
	 *
	 * Note: If useOldSolnForInitGuess is TRUE, phi will be used as the initial guess.
	 */
	template <class T>
	void solveLinSysNSCoupled2(double *phi, T &A, double *rhs,
			bool useOldSolnForInitGuess, const double zeroAbs, const double zeroRel, const int maxIter, int nScal, const int monitorConvergInterval,
			int NcontrolParams = 0, int ncv_gg = 0, double **vecC = NULL, const double *d = NULL, const int step=-1, const int newtonIter=-1);

	template <class T>
	void solveLinSysNSCoupled2(double *phi, T &A, double *rhs, int &nIter, double &absResid,
			bool useOldSolnForInitGuess, const double zeroAbs, const double zeroRel, const int maxIter, int nScal, const int monitorConvergInterval,
			int NcontrolParams = 0, int ncv_gg = 0, double **vecC = NULL, const double *d = NULL, const int step=-1, const int newtonIter=-1);

	template <class T>
	void solveLinSysNSCoupled2(double *phi, T &A, double *rhs, int &nIter, double &absResid, vector<pair<int, double> >& kspMonitorHistory,
			bool useOldSolnForInitGuess, const double zeroAbs, const double zeroRel, const int maxIter, int nScal, const int monitorConvergInterval,
			int NcontrolParams, int ncv_gg, double **vecC, const double *d, const int step, const int newtonIter);

	/*
	 * Method: solveBasicMatOperationsPetsc
	 * ------------------------------------
	 * Solve for basic matrix operations (such as multiplication between a matrix and a vector) using Petsc
	 *
	 * whichOperation == 0: A*b=x
	 */
	template <class MatT>
	void solveBasicMatOperationsPetsc(double *xVec, MatT &A, double *bVec, double *cVec, const int whichOperation,
			int NcontrolParams = 0, int ncv_gg = 0, double **vecC = NULL, const double *d = NULL);

	/*
	 * Method: getPCcontext
	 * --------------------
	 * Get the preconditioner contexts from the input file
	 */
	void getPCcontext(string& pcType, bool &pcReuse, int& pcLevels);

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
	 *     relError  = relative error (to the eigenvalue or to the matrix norms) for each eigen-pair
	 *     numIter   = total number of iterations for the iterative eigen-solver
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
	template <class T>
	int solveEigenProblem(double *evalsReal, double *evalsImag, double **evecsReal, double **evecsImag, double *relError, int &numIter,
			T &A, int *cvora, const int nScal,
			const int nev, const int ncv, const int mpd, const double tol, const int max_iter,
			const int step=-1, const int newtonIter=-1);
#endif

	/*
	 * Method: updatePrimVarsFrom1DVec
	 * -------------------------------
	 * Update the flow field primary variables (rho, rhou, rhoE, scals) from a 1D vector (q0 or q1)
	 * Either q0 or q1 has been gotten from a binary file
	 * Note: You must update nScal before calling this method
	 */
	void updatePrimVarsFrom1DVec(const double* qVec, const int nScal);

	/*
	 * Method: update1DVecFromPrimVars
	 * -------------------------------
	 * Update a 1D vector (q0 or q1) from the flow field primary variables (rho, rhou, rhoE, scals)
	 * The primary variables can be obtained from a Hook function or a restart file, etc.
	 * Note: You must update nScal before calling this method
	 */
	void update1DVecFromPrimVars(double* qVec, const int nScal);

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
	virtual void updateFlow_primVars(const double* Q, const double* delQ, const double relaxation, const int nScal);

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
	virtual void updateLambda(const double* Q, const double* delQ, const double relaxation, const int NcontrolParams, const int nScal);
	
	/*
	 * Method: updateFlow_1Dvec
	 * ------------------------
	 * Update the flow field (Q = Q + sign*delQ). The primary variables (rho, rhou, rhoE, scalars) won't be changed
	 * Arguments: Q = flow vector
	 *            delQ = increment in Q
	 *            relaxation = magnitude of the increment. Usually sign = 1.0 (be careful in the case of PsALC: sign must be -1.0)
	 *            N = size of the 1D flow vectors (Q and delQ)
	 */
	virtual void updateFlow_1Dvec(double* Q, const double* delQ, const double relaxation, const int N);

    /*
     * Method: backtrackForNegativeVals
     * --------------------------------
     * Update relaxation to prevent negative density, pressure, or kine at both CVs and FAs
     */
    double backtrackForNegativeVals(int &negativeValCount_CV, int &negativeValCount_FA, const double clipParameter, const double safeParameter, 
                            const double relaxation, const int kine_index, const int nScal, const double* qArray, const double* delQ, const int iterMax = 1);
    
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
	 *            ResidOld = the residual array (previous step)
	 *            ResidNew = the residual array (next step)
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
	double backtrackWithJOE_calcRelaxAndRHS(double* rhs, const double* qArray, const double* delQ, bool &notConverged,
			const int maxBacktrackingIter, const double relaxationOld, const double relaxationLowerBound, const double relaxationUpperBound,
			const double* ResidOld, double* ResidNew, const int whichNorm, const double backtrackBarrierCoeff,
			const int NcontrolEqns, double** q_tangent, const double* lambda_tangent, const double weightLambda, const double WeightTangentCond,
			const double* qArrayOld, const double* NresOld, const double* lambdaOld, const double arcLength);

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
	double backtrackWithJOEcoupled_calcRelaxAndRHS(double* rhs, const double* qArray, const double* delQ,
			const int maxBacktrackingIter, const double relaxationLowerBound, const double relaxationUpperBound,
			const double* ResidOld, double* ResidNew, const int whichNorm, const double backtrackBarrierCoeff);

	/*
	 * Method: calcRhsWithBarrier
	 * --------------------------
	 * Calculate the single "rhs" array
	 * Barrier functions are also called at the end of this method
	 *
	 * Return: countNegative = the number of negative-valued cells found during the RHS calculation
	 */
	int calcRhsWithBarrier(double* rhs, const bool useBarrier);

	/*
	 * Method: calcBacktrackThreshold
	 * ------------------------------
	 * Calculate the thereshold for the backtracking method (the new residual norm should be smaller than this threshold).
	 * This method will be called by checkBacktrackFromResidInc()
	 */
	double calcBacktrackThreshold(const double residNormTotOld, const double relaxation, const double backtrackBarrierCoeff);

	/*
	 * Method: checkBacktrackFromResidInc
	 * ----------------------------------
	 * Check if there is any increment in the norm of the residual
	 * For details, see C.T.Kelley, Iterative Methods for Linear and Nonlinear Equations, SIAM 1995, Chap.8.1, p.137
	 *
	 * Return: true if there is increment (you need backtracking in this case)
	 *         false if the residual drops (you don't need backtracking in this case)
	 */
	bool checkBacktrackFromResidInc(const double residNormTotOld, const double residNormTot, const double relaxation, const double backtrackBarrierCoeff);

//IKJ
	/*
	 * Method: calcBacktrackThresholdUsingWatchdog
	 * -------------------------------------------
	 * Calculate the thereshold for the backtracking method using the watchdog technique (the new residual norm should be smaller than this threshold).
	 * This method will be called by checkBacktrackFromResidIncUsingWatchdog()
	 */
	double calcBacktrackThresholdUsingWatchdog(const double residNormTotOld, const double relaxation,
			const double gammaWatchdog, const int NcontrolEqns, const int nVars, const double *rhs, const double *phi);

//IKJ
	/*
	 * Method: checkBacktrackFromResidIncUsingWatchdog
	 * -----------------------------------------------
	 * Check if there is any increment in the norm of the residual using the "watchdog" technique:
	 *   Let d_k = the Newton increment (or search direction vector) and h(x) is the objective function (e.g. ||f(x)||_1),
	 *   Then h(x_k + relaxation*d_k) <=   max (h_j(x_k)) + gammaWatchdog*relaxation*g(x_k)^T*d_k,
	 *                                   0<=j<=m
	 *   	where g(x_k) is the gradient vector of h(x):
	 *   		e.g. If 1-norm is used, g(x) = D^T*sign(f(x))
	 *               If 2-norm is used, g(x) = D^T*f(x)
	 *               However, the make the code simpler, g(x) is set to be f(x) here.
	 *
	 * For details, see L.Grippo et al., SIAM J. Numer. Anal., vol.123, 707-716, 1986.
	 *               or R.M.Chamberlain et al., Math. Prog. Study, Vol.16, 1-17, 1982.
	 *
	 * Return: true if there is increment (you need backtracking in this case)
	 *         false if the residual drops (you don't need backtracking in this case)
	 */
	bool checkBacktrackFromResidIncUsingWatchdog(const double residNormTotOld, const double residNormTot, const double relaxation,
			const double gammaWatchdog, const int NcontrolEqns, const int nVars, const double *rhs, const double *phi);

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
	 *   By reference = normSqRatioInArclength ( || q1-q0 ||^2 / || lambda1-lambda0 ||^2 )
	 */
	double calcArclength(double* normSqRatioInArclength,
			const double* q0, const double* q1, const double* lambda0, const double* lambda1,
			const int NcontrolEqns, const int nScal, const double weightForLambda,
			double** q_tangent = NULL, const double* lambda_tangent = NULL);

	/*
	 * Method: calcArclength
	 * ---------------------
	 * Calculate the arc-length
	 * Original formulation: ds = q_tangent*(q1-q0)+weightLambda*lambda_tangent*(lambda1-lambda0),
	 *                       where [q_tangent, lambda_tangent]^T is an unit vector
	 * If [q_tangent, lambda_tangent]^T is not given, it can be easily calculated: q_tangent is equal to (q1-q0)/ds, and lambda_tangent is equal to (lambda1-lambda0)/ds.
	 */
	double calcArclength(const double* q0, const double* q1, const double* lambda0, const double* lambda1,
			const int NcontrolEqns, const int nScal, const double weightForLambda,
			double** q_tangent=NULL, const double* lambda_tangent=NULL);

	/****************************
	 * OVERLOADED FUNCTIONS
	 ****************************/
    /*
     * Method: UgpWithCvCompFlow_AD_init
     * ---------------------------------
     * Initialize some UgpWithCvCompFlow_AD variables
     */
    void UgpWithCvCompFlow_AD_init();

    /*
	 * Method: UgpWithCvCompFlow_AD_clear
	 * ---------------------------------
	 * Clear some UgpWithCvCompFlow_AD variables
	 */
    void UgpWithCvCompFlow_AD_clear() {/*empty*/}

	/*
	 * Method: initialize_adoubles
	 * ---------------------------
	 * Original function = initialize_adoubles() in UgpWithCvCompFlowAD.h
	 * System control parameters are added
	 */
	void initialize_adoubles();

	/*
	 * Method: destroy_adoubles
	 * ------------------------
	 * Original function = destroy_adoubles() in UgpWithCvCompFlowAD.h
	 * System control parameters are added
	 */
	void destroy_adoubles();

#ifdef USE_MEM_SAVING_ADVAR_1D_
	/*
	 * Method: initialize_adoubles
	 * ---------------------------
	 * Original function = initialize_adoubles() in UgpWithCvCompFlowAD.h
	 * System control parameters are added
	 * Memory saving version with the ADvar classes
	 */
	void initialize_adoubles(const int icvCenter, const int NcontrolEqns);

	/*
	 * Method: destroy_adoubles
	 * ------------------------
	 * Original function = destroy_adoubles() in UgpWithCvCompFlowAD.h
	 * System control parameters are added
	 */
	void destroy_adoubles(const int NcontrolEqns);
#endif

	/*
	 * Method: print_tapestats
	 * -----------------------
	 * Original code: print_tapestats() in the JoeWithModels_AD class
	 */
	void print_tapestats(int tag);

	/****************************
	 * HOOK FUNCTIONS
	 ****************************/
	/*
	 * Method: initHookPsALC
	 * ---------------------
	 *
	 */
	virtual void initHookPsALC();

	/*
	 * Method: initHookPsALC_1st
	 * -------------------------
	 *
	 */
	virtual void initHookPsALC_1st();

	/*
	 * Method: initHookPsALC_2nd
	 * -------------------------
	 *
	 */
	virtual void initHookPsALC_2nd();

	/*
	 * Method: PsALCinitialHOOK_debug
	 * ------------------------------
	 *
	 */
	virtual void PsALCinitialHOOK_debug();

	/*
	 * Method: initialHookNewton
	 * -------------------------
	 *
	 */
	virtual void initialHookNewton();

	/*
	 * Method: initialHookNewton_firstCall
	 * -----------------------------------
	 *
	 */
	virtual void initialHookNewton_firstCall();

	/*
	 * Method: temporalHookNewton
	 * --------------------------
	 * The user has an chance to apply an ad-hoc modification to the solution during the Newton iterations.
	 * The user only needs to flow variables (rho, rhou, etc.), qVec, and lambda, but also provides the access to delQ and relaxation just in case.
	 *
	 * You must return TRUE if the solution is modified so that the Newton solver can update flow variables in the CPU boundaries and care other stuff.
	 */
	virtual bool temporalHookNewton(double *qVec, double *lambda, double *delQ, const double relaxation);

	/*
	 * Method: checkBarrierParamNS
	 * ---------------------------
	 * Read parameters for N-S barrier functions and just return if the barrier functions will be used or not
	 * Note: Wrapper function for readBarrierParamNS()
	 */
	bool checkBarrierParamNS(const int iterNewton);

    /*
     * Method: checkBarrierParamTurbScalars
     * ------------------------------------
     * Read parameters for turbuent scalar barrier functions and just return if the barrier functions will be used or not
     * Note: Wrapper function for readBarrierParamScalars()
     */
    bool checkBarrierParamTurbScalars(const int iterNewton);

	/*
	 * Method: readBarrierParamNS
	 * --------------------------
	 * Read parameters for N-S barrier functions
	 * Note: barrierSourceNS(), barrierSourceNS_AD(), and barrierSourceNS1D_AD() are calling this method
	 * Return: true if the barreir functions will be used
	 */
	bool readBarrierParamNS(string& barrierFunctionName, int& maxIter, double& threshold_resid, double& coeffRho, double& coeffPress,
			const int iterNewton, const double resid_norm_tot = ABSURDLY_BIG_NUMBER);

    // Note: readBarrierParamTurbScalars() is in IkeUgpWithCvCompFlow.h

	/*
	 * Method: barrierSourceNS
	 * -----------------------
	 * Add barrier functions to the RHS of the N-S equations
	 */
	virtual void barrierSourceNS(double* rhs);

    // Note: barrierSourceTurbScalars() is in IkeUgpWithCvCompFlow.h

	/*
	 * Method: barrierSourceNS_AD
	 * --------------------------
	 * Add barrier functions to the RHS of the N-S equations
	 */
	virtual void barrierSourceNS_AD(REALA *rhs_rho_AD, REALA (*rhs_rhou_AD)[3], REALA *rhs_rhoE_AD);

    // Note: barrierSourceTurbScalars_AD() is in IkeUgpWithCvCompFlow.h

	/*
	 * Method: barrierSourceNS1D_AD
	 * ----------------------------
	 * Add barrier functions to the RHS of the N-S equations
	 */
	virtual void barrierSourceNS1D_AD(const int icvCenter, REALA &rhs_rho_AD, REALA rhs_rhou_AD[3], REALA &rhs_rhoE_AD);

    // Note: barrierSourceTurbScalars1D_AD() is in IkeUgpWithCvCompFlow.h

	/****************************
	 * UTILITY FUNCTIONS
	 ****************************/
	/*
	 * Method: getParamsLinSolverNewton
	 * --------------------------------
	 * Obtain the linear solver parameters for Newton's method
	 * Parameters:
	 *   1. Linear solver thresholds
	 *   	maxIterLS = maximum number of Newton (outer) iterations
	 *   	zeroAbsLS = absolute residual of Newton (outer) iterations
	 *   	zeroRelLS = relative residual of Newton (outer) iterations
	 *   2. Ramp
	 * 		startDecLS     = After "startDecLS" iteration, ramp starts
	 * 		intervalDecLS  = At every "intervalDecLS" iteration, ramp occurs
	 * 		incIterLS      = When ramps occurs, "maxIterLS" increases by "incIterLS"
	 * 		maxFinalIterLS = After "maxFinalIterLS" iteration, ramp no longer occurs
	 * 		decZeroLS      = When ramps occurs, "zeroAbsLS" and "zeroRelLS" drops by "decZeroLS"
	 * 		minZeroAbsLS   =
	 * 		minZeroRelLS   =
	 */
	void getParamsLinSolverNewton(int &maxIterLS, double &zeroAbsLS, double &zeroRelLS,
			int &startDecLS, int &intervalDecLS, int &incIterLS, int &maxFinalIterLS, double &decZeroLS, double &minZeroAbsLS, double &minZeroRelLS,
			const bool firstCall, const int debugLevel);

	/*
	 * Method: getParamsBacktrackNewton
	 * --------------------------------
	 * Obtain the backtracking parameters for the Newton method
	 * Parameters:
	 *   backtrackMaxIter 			= "BACKTRACKING_PARAMETERS"->"MAX_ITER"
	 *   backtrackRelax_LowerBound 	= "BACKTRACKING_PARAMETERS"->"RELAX_LOWER_BOUND"
	 *   backtrackRelax_UpperBound 	= "BACKTRACKING_PARAMETERS"->"RELAX_UPPER_BOUND"
	 *   backtrackBarrierCoeff 		= "BACKTRACKING_PARAMETERS"->"BARRIER_COEFF"
	 *
	 *   skipBT_firstITer  = "BACKTRACKING_SKIP"->"FIRST_ITER"
	 *   skipBT_freq       = "BACKTRACKING_SKIP"->"FREQUENCY"
	 */
	void getParamsBacktrackNewton(int &backtrackMaxIter, double &backtrackRelax_LowerBound, double &backtrackRelax_UpperBound, double &backtrackBarrierCoeff,
			bool &skipBT_firstITer, int &skipBT_freq,
			const bool firstCall, const int debugLevel);

	/*
	 * Method: getParamStabilizedNewton
	 * --------------------------------
	 * Obtain the parameters for the stabilized Newton's method for the problems with singular Jacobians
	 *   "STABILIZATION_FOR_SINGULAR_JAC"-->"MATRIX_TYPE"
	 *   "STABILIZATION_FOR_SINGULAR_JAC"-->"ALPHA"
	 *   "STABILIZATION_FOR_SINGULAR_JAC"-->"ALPHA_EPS"
	 *   "STABILIZATION_FOR_SINGULAR_JAC"-->"STARTING_ITER"
	 */
	void getParamStabilizedNewton(int &stabilizationMethodType, double &stabilizationAlpha, double &stabilizationAlphaEps, int &stabilizationStartingIter,
			const bool firstCall, const int debugLevel);

	/*
	 * Method: getParamsModifiedNewtonMethod
	 * -------------------------------------
	 * Obtain the Newton method to be used.
	 * Currently supported methods = BASIC, SHAMANSKII, MODIFIED_SHAMANSKII
	 */
	void getParamsModifiedNewtonMethod(MODIFIED_NEWTON& modifiedNewtonMethod, double& modifNewtonRelResid, int& modifNewtonFreq);

	/*
	 * Method: getParamsReducingRelax
	 * ------------------------------
	 * Get the parameters for the relaxation reducing algorithm:
	 *   If positive values (rho, press, kine) become negative during Newton's method, reduce "relaxtion"
	 */
	void getParamsReducingRelax(double &clipParameter, double &safeParameter, bool showOnScreen);

	/*
	 * Method: getParamsBEuler
	 * -----------------------
	 * Obtain some parameters for the backward-Euler time-integration method
	 */
	void getParamsBEuler(int &maxIterLS_BEuler, double &zeroAbsLS_BEuler, double &zeroRelLS_BEuler,
			double &CflBEuler, int &startIncCFLBEuler, int &intervalIncCFLBEuler, double &incCFLBEuler, double &maxCFLBEuler);

	/*
	 * Method: getHowToCalcJac
	 * -----------------------
	 * Obtain the method to calculate the Jacobian matrix.
	 * Currently supported methods = ROW_1D, ORDINARY_2D
	 */
	HOW_TO_CALC_JAC getHowToCalcJac();

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
	double vecDotVecWithoutWeight (double* innerProductRatio,
			const double* qVec0, const double* qVec1, const double* lambdaVec0, const double* lambdaVec1,
			const int Nvars, const int NcontrolEqns);

	double vecDotVecWithoutWeight (const double* qVec0, const double* qVec1, const double* lambdaVec0, const double* lambdaVec1,
			const int Nvars, const int NcontrolEqns);

	/*
	 * Method: getFlowWeightsForInnerProduct
	 * -------------------------------------
	 * Get weights for the inner product between flow vectors
	 */
	void getFlowWeightsForInnerProduct(double* weightVec, const double defaultConst, const int icv, const int nVars);

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
	double vecDotVecWithWeight(double* innerProductRatio,
			const double* qVec0, const double* qVec1, const double* lambdaVec0, const double* lambdaVec1, const double weightForLambda,
			const int Nvars, const int NcontrolEqns);

	double vecDotVecWithWeight(const double* qVec0, const double* qVec1, const double* lambdaVec0, const double* lambdaVec1, const double weightForLambda,
			const int Nvars, const int NcontrolEqns);

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
	double calcUnweighted2NormForVecMinusVec(double *normSqRatio,
			const double* qVec1, const double* qVec0, const double* lambdaVec1, const double* lambdaVec0,
			const int Nvars, const int NcontrolEqns);

	double calcUnweighted2NormForVecMinusVec(const double* qVec1, const double* qVec0, const double* lambdaVec1, const double* lambdaVec0,
			const int Nvars, const int NcontrolEqns);

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
	double calcWeighted2NormForVecMinusVec(double *normSqRatio,
			const double* qVec1, const double* qVec0, const double* lambdaVec1, const double* lambdaVec0, const double weightForLambda,
			const int Nvars, const int NcontrolEqns);

	double calcWeighted2NormForVecMinusVec(const double* qVec1, const double* qVec0, const double* lambdaVec1, const double* lambdaVec0, const double weightForLambda,
			const int Nvars, const int NcontrolEqns);

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
	void calcUnitVecWithWeightedNorm(double* u_q, double* u_lambda, const double* qVec, const double* lambdaVec, const double weighted2norm,
			const int Nvars, const int NcontrolEqns);

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
	double calcUnitVec(double *normSqRatio,
			double* u_q, double* u_lambda, const double* qVec, const double* lambdaVec, const double weightForLambda,
			const int Nvars, const int NcontrolEqns);

	double calcUnitVec(double* u_q, double* u_lambda, const double* qVec, const double* lambdaVec, const double weightForLambda,
			const int Nvars, const int NcontrolEqns);

	/*
	 * Method: calcResidualsFrom1Drhs
	 * ------------------------------
	 * Calculate residual
	 * Supported norm: inf-norm, one-norm, two-norm
	 *                   whichNorm==0 -> infinity-norm
	 *                   whichNorm==1 -> one-norm
	 *                   whichNorm==2 -> two-norm
	 * Note: "Residual" should have been defined as " double *Residual = new double[5+nScal]; "
	 */
	virtual void calcResidualsFrom1Drhs(double *Residual, double *rhs1Darray, const int whichNorm);

	/*
	 * Method: calcTotResidual
	 * -----------------------
	 * Calculate total residual from the residual vector
	 * Supported norm: inf-norm, one-norm, two-norm
	 *                   whichNorm==0 -> infinity-norm
	 *                   whichNorm==1 -> one-norm
	 *                   whichNorm==2 -> two-norm
	 * Note: "Residual" should have been defined as " double *Residual = new double[5+nScal]; "
	 */
	double calcTotResidual(const double *Residual, const int whichNorm);

	/*
	 * Method: calcNres
	 * ----------------
	 * Calculate Nres == q_tangent*(q_guess-q1) + lambda_tangent*(lambda_guess-lambda1) - ds)
	 * (Nres is the last element of the rhs vector at mpi_rank=mpi_size-1)
	 *
	 * Return: True if the tangential condition is satisfied
	 *
	 * Note: Even though this part of the code is designed to handle multiple parameters (or lambdas), the actual formulation here is only for one parameters
	 */
	bool calcNres(double* Nres, const double *q, double* rhs, const int Nvars, const int NcontrolEqns, const double *q1, double** q_tangent,
			const double* lambda_tangent, const double *lambda1, const double weightLambda, const double arcLength);

	/*
	 * Method: calcNresForJOE
	 * ----------------------
	 * Calculate Nres == q_tangent*(q_guess-q1) + lambda_tangent*(lambda_guess-lambda1) - ds)
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
	bool calcNresForJOE(double* Nres, const int Nvars, const int NcontrolEqns, const double *q1, double** q_tangent,
			const double* lambda_tangent, const double *lambdaNew, const double *lambda1, const double weightLambda, const double arcLength);

	/*
	 * Method: UpdateTotResidualWithNres
	 * ---------------------------------
	 * Update total residual norm with the Nres array
	 */
	double updateTotResidualWithNres(double residNormTot, const double* Nres, const int NcontrolEqns, const int whichNorm);

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
	void write1DArrayWithXcvSerial(const string &filename, const int step, const double *qVec);
	void write1DArrayWithXcvSerial(const char filename[], const int step, const double *qVec);

	/*
	 * Method: writePsALCdumpedDataSerial
	 * ----------------------------------
	 * Write data for an IKE simulation
	 */
	void writePsALCdumpedDataSerial(const string &filename, const int step, const double dsTry, const double *lambda0, const double *lambda1, const int NcontrolEqns, const double *qVec);
	void writePsALCdumpedDataSerial(const char filename[], const int step, const double dsTry, const double* lambda0, const double* lambda1, const int NcontrolEqns, const double *qVec);

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
	void writePsALCdumpedDataParallel(const string &filename, const int step, const double dsTry, const double *lambda0, const double *lambda1, const int NcontrolEqns, double *qVec);
	void writePsALCdumpedDataParallel(char filename[], const int step, const double dsTry, const double *lambda0, const double *lambda1, const int NcontrolEqns, double *qVec);

	/*
	 * Method: readPsALCdumpedDataSerial
	 * ---------------------------------
	 * Read data from previous IKE simulation
	 */
	void readPsALCdumpedDataSerial(const string &filename, double* qVec, int& step0, double *lambda0, double *lambda1, const int NcontrolEqns, double& dsTry);
	void readPsALCdumpedDataSerial(const char filename[], double* qVec, int& step0, double* lambda0, double* lambda1, const int NcontrolEqns, double& dsTry);

	/*
	 * Method: readJOEdumpedData
	 * -------------------------
	 * Read data from previous JOE simulation
	 */
	void readJOEdumpedData(const string &filename, double* q);

	/*
	 * Method: findMatchedIndex
	 * ------------------------
	 * This method is used to find a matched index for the readPsALCdumpedData() method
	 */
	int findMatchedIndex(int& foundIndex, const int icv, const int ncv_file, const double (*x_cv_file)[3], const double gridTol);

	/*
	 * Method: writeMatrixMfile
	 * ------------------------
	 * write the matrix on files as a matlab format (mfile format)
	 * Sometimes, the size of the matrix is too large to be saved on a single file.
	 */
	template <class MatT>
	void writeMatrixMfile(string &filename, MatT &jacMatrix, const int nScal, const int *cvora, const int *nbocv_v_global);

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
	void writeMatrixBinarySerial(string &filename, MatT &jacMatrix, const int nScal, const int *cvora, const int *cv_gl);

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
	void writeMatrixBinaryParallel(string &filename, MatT &jacMatrix, const int nScal, const int *cvora, const int *cv_gl);

	/*
	 * Method: readMatrixBinary
	 * ------------------------
	 * Read a matrix from a binary file and store it in a matrix container
	 */
	template <class MatT>
	void readMatrixBinary(string &filename, MatT &jacMatrix, const int nScal, const int *cvora, const int *cv_gl);

	/*
	 * Method: readMatrixBinaryParallel
	 * --------------------------------
	 * Read a matrix from a binary file using MPI-2 parallel I/O
	 * Output: myNnz_fromFile = local number of non-zeros (NOT the total nnz -- it is the nnz for a core)
	 *         rind           = GLOBAL row indices array
	 *         cind           = GLOBAL column indices array
	 *         values         = values array
	 */
	void readMatrixBinaryParallel(int &myNnz_fromFile, vector<unsigned int> &rind, vector<unsigned int> &cind, vector<double> &values,
			string &filename, const int nVars, const int *cvora, const int *cv_gl);

	/*
	 * Method: writeQoIOnFile
	 * ----------------------
	 * By default, averaged density, maximum pressure, and subsonic portion will be stored
	 */
	virtual void writeQoIOnFile(const int step, string filename, const int NcontrolEqns, bool rewrite);

	/*
	 * Method: calcAvgRho
	 * ------------------
	 */
	virtual double calcAvgRho();

	/*
	 * Method: calcMaxPress
	 * --------------------
	 */
	virtual double calcMaxPress();

	/*
	 * Method: calcSubsonicPortion
	 * ---------------------------
	 */
	virtual double calcSubsonicPortion();

	/*
	 * Method: convertSeparatedRhsTo1Drhs
	 * ----------------------------------
	 * Convert separated rhs (used in calcRhs()) to 1D rhs
	 */
	void convertSeparatedRhsTo1Drhs(double* rhs1D, double *RHSrho, double (*RHSrhou)[3], double *RHSrhoE, double **RHSrhoScal, const int nScal);

	/*
	 * Method: convert1DrhsToSeparatedRhs
	 * ----------------------------------
	 * Update RHSrho, RHSrhou, RHSrhoE, (and scalars) from rhs1Darray.
	 * These residual arrays will be used for tecplot output
	 *
	 * Note: If nScal==1, assume SA model is used
	 *       If nScal==2, assume KW model is used
	 */
	void convert1DrhsToSeparatedRhs(const double *rhs1Darray);

	/*
	 * Method: convertRHSrhoScalToSeparatedRhs
	 * ---------------------------------------
	 * Only for scalars!! : Update RHSsa or (RHSkine, RHSomega) from RHSrhoScal.
	 * These residual arrays will be used for tecplot output
	 *
	 * Note: If nScal==1, assume SA model is used
	 *       If nScal==2, assume KW model is used
	 */
	void convertRHSrhoScalToSeparatedRhs(double **RHSrhoScal);

	/*
	 * Method: convert2DrhsTo1Drhs
	 * ---------------------------
	 * Convert 2D rhs (used in calcRhsCoupled()) to 1D rhs
	 */
	void convert2DrhsTo1Drhs(double* rhs1D, double** rhs2D, const int nScal);

	/*
	 * Method: showResidue2
	 * --------------------
	 * Show residuals of the flow variables on the screen
	 * Original code = showResidue() in JoeWithModels.cpp
	 */
	void showResidue2(double *rhsResid, const bool showLabel);
	void showResidue2(double *rhsResid, const int step, const int check_interval);

	/****************************
	 * FOR DEBUGGING
	 ****************************/
	/*
	 * Method: compatibilityCheck
	 * --------------------------
	 * The Psdueo-arclength continuation code uses two different types of RHS calculations: ADOL-C (for Jacobian calculation) and normal JOE (for backtracking)
	 * This method checks if the two different calculation can give the same results or not.
	 * Note: The RHS calculations should be same as those in calcJacobian1DAD()/calcJacobianAD() and backtrackWithJOE()
	 *
	 * Return: false if the two calculations are different
	 */
	bool compatibilityCheck();

	/*
	 * Method: checkNanOrInfRhs1D_AD
	 * -----------------------------
	 *
	 */
	bool checkNanOrInfRhs1D_AD(const int icv, adouble &rhs_rho_AD, adouble *rhs_rhou_AD, adouble &rhs_rhoE_AD, adouble *rhs_rhoScal_AD);

	/*
	 * Method: checkNanOrInfRhs_AD
	 * ---------------------------
	 *
	 */
	bool checkNanOrInfRhs_AD(REALA *rhs_rho_AD, REALA (*rhs_rhou_AD)[3], REALA *rhs_rhoE_AD, REALA **rhs_rhoScal_AD);

	/*
	 * Method: showBoundaryInfoOnScreen
	 * --------------------------------
	 *
	 */
	void showBoundaryInfoOnScreen(string& boundaryName);
	void showBoundaryInfoOnScreen(char boundaryName[]);

	/*
	 * Method: showFlowfieldOnScreen
	 * -----------------------------
	 * Display data on the screen
	 * You must provide gammaRef and RoMRef for the case that gamma and RoM were not established yet.
	 */
	void showFlowfieldOnScreen(double gammaRef, double RoMRef, bool showGhosts = false, bool showFakes = false);

	/*
	 * Method: show1DArrayOnScreen
	 * -----------------------------
	 * Display data on the screen with the x_cv information
	 * Argument: Q = 1D array (size = ncv*nVars, where nVars = 5+nScal)
	 *
	 * Note: This is a MPI call. Thus, all the cpu cores print out its data on the screen.
	 */
	void show1DArrayOnScreen(const double* Q, const int nVars);

	/*
	 * Method: show1DArrayOnScreen
	 * ---------------------------
	 * Display data on the screen with the x_cv information
	 * (For each entry in Q, if icv<ncv, show the x_cv information. Otherwise, just show the value of the entry)
	 * Argument: Q = 1D array (size = arrayLength)
	 *
	 * Note: This is a MPI call. Thus, all the cpu cores print out its data on the screen.
	 */
	void show1DArrayOnScreen(const double* Q, const int nVars, int arrayLength, const bool showIcv=true);

protected:
	int nScal; // This will be initialized in IkeWithPsALC_AD::run()

	// Minimum and Maximum of x_cv among icv=0~ncv (This will be used when the data is written on or read from files)
	double xcvMin[3], xcvMax[3]; // This values will be updated in IkeWithPsALC_AD::run()

#ifdef USE_VOLUME_WEIGHTED_INNER_PROD
	// Minimum and Maximum of cv_volume
	double totMaxCvVol, totMinCvVol;
#endif

	// Jacobian matrix
	MatComprsed jacMatrix;       // This will be used if the Jacobian matrix is calculated as whole
	MatComprsedSTL jacMatrixSTL; // This will be used if the Jacobian matrix is calculated in the so-called "1D style"

	// PetscSolver
#ifdef USE_SLEPC_WITH_PETSC
	SlepcSolver2 *petscSolver2;
#else
	PetscSolver2 *petscSolver2;
#endif

	// Pseudo-arclength continuation method
	int NcontrolEqns; // NcontrolEqns will be given in the runPsALC() method

	// RHS 1D array
	double *rhs1DArray;

	// RHS arrays (for Tecplot output)
	double *RHSrho;
	double (*RHSrhou)[3];
	double *RHSrhoE;
	double *RHSkine;
	double *RHSomega;
	double *RHSsa;

	// Heat: lambda and lambda_AD were moved to the UgpWithCvCompFlow_AD class so that other functions can access to them
	//       See the IkeUgpWithCvCompFlow.h file for the details

	// Total residual (flow residual + tangential condition)
	double residNormTotOld; // Total residual at the previous Newton step
	double residNormTot;    // Total residual at the current Newton step

	// Sum of the source terms due to Barrier functions
	double myBarrierMassSourceSumJOE;
	double myBarrierEnergySourceSumJOE;
	double myBarrierMassSourceSum1D_AD;
	double myBarrierEnergySourceSum1D_AD;

	// QoI
	bool rewriteQoI;

	// Starting step and initialization method for the 3rd point
	int startingStep;     // STARTING_STEP
	bool initThirdFromQ1; // INIT_THIRD_FROM_Q1

	// Debugging level
	int debugLevel;

	// Wall time
	double initWtime;
	wallTimeJacCalc myWTimeJacCalc;

	//
	double initResidNorm; // residual of the initial guess of the Newton's method (Declare this as a member variable so that any functions can get the access to it)
	int iterNewton; // iteration number of the Newton's method (Declare this as a member variable so that any functions can get the access to it)

	int    myNonDiagAdd_count;     // Statistics: the number of modified diagonal
	double myNonDiagAdd_rhsAbsSum; // Statistics: the increase of RHS due to this treatment

	// Interval to check the convergence of PETSc KSP
	int monitorConvergInterval; // This will be initialized in IkeWithPsALC_AD::runPsALC()

	// GMRes restart
	int gmresRestart;
};

#endif /* IKEWITHPSALC_H_ */
