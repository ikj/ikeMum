/*
 * rembrandtWithModels.h
 *
 *  Created on: Jul 21, 2015
 *      Author: ikj
 */

#ifndef REMBRANDTWITHMODELS_H_
#define REMBRANDTWITHMODELS_H_

#include "JoeWithModels.h"
#include "JOE/ADJOINT_FILES/JoeWithModelsAD.h"
#include "adolc.h"

#include "IkeTurbModel_KOM.h"
#include "turbModels/TurbModel_KOM.h"
#include "ADJOINT_FILES/TurbModel_KOM_AD.h"

#include "IkeWithPsALC.h"

#include <sstream>
#include <complex>

#include <stdexcept>
#include <iomanip>

enum HOW_TO_GET_JACOBIAN {FROM_MAT_BINARY_FILE, CALC_USING_ADOLC};

#define VERY_SMALL_DOUBLE 1.0e-10

#define REMBRANDT_ERROR_CODE -3

#define A_LARGE_TOT_FLOW_RESID 1.0e-4

#define THRESHOLD_EVAL_DIFF 1.0e-8

#define COEFF_BOOST_SENS 1.0e6  // Sensitivity becomes very tiny (because it is a multiplication between two small numbers).
                                // Thus, scale the sensitivity so that it is much larger than the machine zero.

#define FILENAME_DIRECT_EVALS "SLEPC_DIRECT_EVALS.csv"

/*
 * Inline function: calcComplexMag
 * -------------------------------
 * Calculate the magnitude of a complex number
 */
inline double calcComplexMag(const double RealPart, const double ImagPart) {
	return sqrt(RealPart*RealPart + ImagPart*ImagPart);
}

// ###########################################################################################
// ------                                                                               ------
// ------                                    REMBRANDT                                  ------
// ------     Perform eigen-decomposition with SLEPSC and Export the data in TECPLOT    ------
// ------                                                                               ------
// ###########################################################################################
class RembrandtWithModels : public IkeWithPsALC_AD {
protected:
	// Class name
	string classID;

	// Compatibility check
	bool compCheckOn;  // We want to allow lambda_AD only during the compatibility check

	// Number of global modes that the user wants to calculate
	int nev;

	// -------------------------
	// Global-modes calculations
	// -------------------------
	// DIRECT eigen-values and eigen-vectors
	double*  directEvalsReal;
	double*  directEvalsImag;
	double** directEvecsReal;
	double** directEvecsImag;

	// ADJOINT eigen-values and eigen-vectors
	double*  adjointEvalsReal;
	double*  adjointEvalsImag;
	double** adjointEvecsReal;
	double** adjointEvecsImag;

	// -------------------------
	// Tecplot output
	// -------------------------
	// The least-stable DIRECT modes (Since the least-stable modes are real, we don't store their imaginary parts)
	double *FirstDirectMode_real_rho;
	double (*FirstDirectMode_real_rhou)[3];
	double (*FirstDirectMode_real_vel)[3];
	double *FirstDirectMode_real_rhoE;
	// NOTE: scalars should be registered in the child class

	// The least-stable ADJOINT modes (Since the least-stable modes are real, we don't store their imaginary parts)
	double *FirstAdjointMode_real_rho;
	double (*FirstAdjointMode_real_rhou)[3];
	double (*FirstAdjointMode_real_vel)[3];
	double *FirstAdjointMode_real_rhoE;
	// NOTE: scalars should be registered in the child class

	// Structural sensitivity
	double *rho_sens;
	double (*rhou_sens)[3];
	double (*vel_sens)[3];
	double *rhoE_sens;
	// NOTE: scalars should be registered in the child class



public:
	RembrandtWithModels(char *name) : IkeWithPsALC_AD(name), JoeWithModels(name), JoeWithModels_AD(name), UgpWithCvCompFlow(name) {
		classID = "RembrandtWithModels";
		if(mpi_rank==0)
			cout<<classID<<"()"<<endl;

		compCheckOn = false;  // We want to allow lambda_AD only during the compatibility check

		init();
	}
	virtual ~RembrandtWithModels() {
		if(mpi_rank==0)
			cout<<"~"<<classID<<"()"<<endl;
		clear();
	}

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

	/*
	 * Method: run
	 * -----------
	 *
	 */
	void run();

	/*
	 * Method: obtainJacobianMat
	 * -------------------------
	 * Obtain the Jacobian matrix: There are two possible options --
	 *   1. Calculate the Jacobian matrix using ADOL-C
	 *   2. Read the Jacobian matrix from a file (still in development)
	 */
	void obtainJacobianMat();

	/*
	 * Method: runEigenAnalysis
	 * ------------------------
	 * Calculate both direct and adjoint global-modes and perform the required post-processing.
	 * Return = the number of converged eigenpairs (min in Direct and Adjont)
	 */
	int runEigenAnalysis();

	/*
	 * Method: postEigenDecompScalarRansTurbModel
	 * ------------------------------------------
	 *
	 */
	virtual void postEigenDecompScalarRansTurbModel(double** directEvecsReal, double** adjointEvecsReal,
			const int directLeastStableIndex, const int adjointLeastStableIndex) { /* empty */ };

	/*
	 * Method: postEigenDecompScalarRansCombModel
	 * ------------------------------------------
	 *
	 */
	virtual void postEigenDecompScalarRansCombModel(double** directEvecsReal, double** adjointEvecsReal,
			const int directLeastStableIndex, const int adjointLeastStableIndex) { /* empty */ };

	/*
	 * Method: showDataRangeEigenAnalysis
	 * ----------------------------------
	 * Show data range of the least-stable direct-mode, the least-stable ajoint-mode, and the sensitivity
	 */
	void showDataRangeEigenAnalysis();

	/*
	 * Method: readMatrixBinaryRembrandt
	 * ---------------------------------
	 * Original function = IkeWithPsALC_AD::readMatrixBinary()
	 * For some reason, "undefined reference to" error occurs by the linker. Thus, define the same function in Rembrandt.
	 */
	template <class MatT>
	void readMatrixBinaryRembrandt(string &Q1Filename, string &JacFilename, MatT &jacMatrix,
			const int nScal, const int *cvora, const int *cv_gl);

	/*
	 * Method: readMatchedPtsFromPsALCdumpedDataSerial
	 * -----------------------------------------------
	 * Read data from previous IKE simulation
	 */
	void readMatchedPtsFromPsALCdumpedDataSerial(const string &filename,
			int* mpi_rank_inFile, int* icv_inFile, set<int>& Set_mpi_rank_inFile, const int NcontrolEqns);

	void readMatchedPtsFromPsALCdumpedDataSerial(const char filename[],
			int* mpi_rank_inFile, int* icv_inFile, set<int>& Set_mpi_rank_inFile, const int NcontrolEqns);

//	/*
//	 * Method: readMatrixBinaryForRembrandtParallel
//	 * --------------------------------------------
//	 * Read a matrix from a binary file using MPI-2 parallel I/O
//	 * Output: myNnz_fromFile = local number of non-zeros (NOT the total nnz -- it is the nnz for a core)
//	 *         rind           = GLOBAL row indices array
//	 *         cind           = GLOBAL column indices array
//	 *         values         = values array
//	 */
//	void readMatrixBinaryForRembrandtParallel(int &nnz_local,
//			vector<unsigned int> &rind_local_unsorted, vector<unsigned int> &cind_local_unsorted, vector<double> &values_unsorted,
//			string &JacFilename, const int* mpi_rank_inFile, const int* icv_inFile, set<int>& Set_mpi_rank_inFile,
//			const int nVars, const int *cvora_inFile, const int *cv_gl_inFile);

	/*
	 * Method: calcMatrixFromField
	 * ---------------------------
	 * Wrapper function for runEigenAnalysis
	 */
	void calcMatrixFromField();

	/*
	 * Method: calcRhsWithJacMatrixRembrandt
	 * -------------------------------------
	 * Calculate the Jacobian matrix as IKE
	 */
	template <class MatT>
	void calcRhsWithJacMatrixRembrandt(MatT &jacMatrix, double* rhsVec, const int nScal);

	/*
	 * Method: calcRhsWithJacMatrixRembrandt
	 * -------------------------------------
	 * Original code = IkeWithPsALC_AD::getSteadySolnByNewton
	 * Since IkeWithPsALC_AD::getSteadySolnByNewton() deletes jacMatrix in the routine,
	 * this function must be called in order to get the jacMatrix.
	 */
	void calcSteadySolnByNewtonRembrandt(double* q, double* rhs,
			const int maxIterNewton, const double absTolNewton, const double relTolNewton, const int nScal);

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
	 *       eigen-vectors   (variable names = directEvecsReal, directEvalsImag)
	 *       adjoint vectors (variable names = adjointEvalsReal, adjointEvalsImag)
	 *   Foot:
	 *       EOF_ERROR_CHECK_CODE
	 */
	void writeEigenPairsParallel(const string &filename, const int step, const int nev);

	void writeEigenPairsParallel(char filename[], const int step, const int nev);

protected:
	/*
	 * Method: calcGlobalModes
	 * -----------------------
	 * Calculate both direct and adjoint global-modes
	 */
	int calcGlobalModes(const int nevSlepc);

	/*
	 * Method: getSlepcEPSParams
	 * -------------------------
	 * Get SLEPc EPS (eigen problem solver parameters) from the input file
	 *
	 * Arguments:
	 *   nevSlepc = This is required only for printing output
	 */
	void getSlepcEPSParams(SlepcEPSParams &slepcEPSParams, const int nevSlepc);

	/*
	 * Method: getEigenOfInterestFromString
	 * ------------------------------------
	 *
	 */
	EPSWhich getEigenOfInterestFromString(string& tempString);
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

	/*
	 * Method: getStringFromEigenOfInterest
	 * ------------------------------------
	 *
	 */
	string getStringFromEigenOfInterest(const EPSWhich eigenOfInterest);

	/*
	 * Method: getEPSsolverTypeFromString
	 * ----------------------------------
	 *
	 */
	EPSType getEPSsolverTypeFromString(string &tempString);
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
	void getSlepcSTparams(SlepcSTparams &slepcSTparams, const SlepcEPSParams &slepcEPSParams);

	/*
	 * Method: getEvecRange
	 * --------------------
	 * Get the ranges of eigen-vectors
	 * Return = minVal_Real[], maxVal_Real[], minVal_Imag[], maxVal_Imag[]
	 */
	void getEvecsRange(double* minVal_Real, double* maxVal_Real, double* minVal_Imag, double* maxVal_Imag,
			double** EvecsReal, double** EvecsImag, const int nevSlepc);

	/*
	 * Method: checkTwoEvecsDifferenceRange
	 * ------------------------------------
	 * Get the ranges of the difference between two eigen-vectors
	 * Return = 1-norm of real and imaginary parts
	 */
	void checkTwoEvecsDifferenceRange(double* diff1Norm_real, double* diff1Norm_imag, double* diffInfNorm_real, double* diffInfNorm_imag,
			double** EvecsReal_1st, double** EvecsImag_1st, double** EvecsReal_2nd, double** EvecsImag_2nd,
			const int nevSlepc);
};

#endif
