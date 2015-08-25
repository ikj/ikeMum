/*
 * vanGoghWithModels.h
 *
 *  Created on: Aug, 2015
 *      Author: ikj
 */

#ifndef VANGOGHWITHMODELS_H_
#define VANGOGHWITHMODELS_H_

#include "JoeWithModels.h"

#include <cstdlib>
#include <ctime>
#include <math.h>

#include <sstream>
#include <complex>

#include <stdexcept>
#include <iomanip>

#include "vanGoghUtils.h"
#include "DiffFilter.h"

enum RAND_DISTRIB_FUNC {NORMAL_DISTRIB, UNIFORM_DISTRIB};

#ifndef PI
#define PI 3.14159265359
#endif

#ifndef ABSURDLY_BIG_NUMBER
#define ABSURDLY_BIG_NUMBER 2.2e22
#endif

#ifndef MACHINE_EPS
#define MACHINE_EPS 2.0e-16 // Note: machine epsilon for double precision = 1.11e-16
#endif

#ifndef EOF_ERROR_CHECK_CODE
#define EOF_ERROR_CHECK_CODE 503   // This value must be the same as IkeWithPsALC.h
#endif

#define VANGOGH_ERROR_CODE -4

#define VANGOGH_DEBUG_LEVEL 1

// TO DO: MATLAB version of file read!!!!
#define NVARS_EIGEN 7
#define NEIGENS 30

/*
 * Note:
 * You must make JoeWithModels::run() virtual
 * You must make the "Residual" array as a member variable of the JoeWithModels class
  and comment out the "Residual" array in runForwardEuler(), runRK(), runBackwardEuler(), runBackwardEulerCoupled(), and runBDF2();
 */

struct PerturbParams {
	PerturbParams() {
		useSmoothing = false;
		useFiltering = false;
	}

	RAND_DISTRIB_FUNC randDistribFunc;
	double disturbMag;
	double disturbClip;

	bool useSmoothing;
	double disturbSmoothXmin;
	double disturbSmoothXmax;
	double disturbSmoothEdgeSize;

	bool useFiltering;
    double minFilterWidth;
	double maxFilterWidth;
	double filterBoundaryLength;

	RAND_DISTRIB_FUNC RandDistribFunc;
};

// ###########################################################################################
// ------                                                                               ------
// ------                       Perturbation tests called VanGogh                       ------
// ------                                                                               ------
// ###########################################################################################
class VanGoghWithModels : virtual public JoeWithModels, public DiffFilter {
protected:
	int NcontrolEqns;
	double *lambda;
	
	string classID;

	//**************
	// VAN-GOGH
	//**************
	// test number
	int itest;
	int ntests;
	
	// xcvMin & xcvMax (will be used while reading Q1.bin file from IKE or EigenPairs.bin file)
	double xcvMin[3];
	double xcvMax[3];
	
	// initial field (from a restart file)
	double *rho_init;
	double (*rhou_init)[3];
	double *rhoE_init;
	
	// data on the unstable(hyperbolic) branch (from Q0_PT000**.bin)
	double *rho_unstable;
	double (*rhou_unstable)[3];
	double *rhoE_unstable;

	// -------------------------
	// Perturbations
	// -------------------------
	// perturbations
	double *array_perturb;

	// filter
	PerturbParams perturbParams;
	double *wdFace;

	// -------------------------
	// Eigenvecs: Data container
	// -------------------------
	vector<vector<vector<std::complex<double> > > >  DirectEvecs;
	vector<vector<vector<std::complex<double> > > >  AdjointEvecs;

	// -------------------------
	// Eigenvecs: Tecplot output
	// -------------------------
	// least-stable direct global modes
	double *rho_Dct_1stReal;
	double (*rhou_Dct_1stReal)[3];
	double *rhoE_Dct_1stReal;

	// first adjoint global modes
	double *rho_Adj_1stReal;
	double (*rhou_Adj_1stReal)[3];
	double *rhoE_Adj_1stReal;
	
public:
	/*
	 * constructor
	 */
	VanGoghWithModels(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name) {
		classID = "VanGoghWithModels";
		if(mpi_rank==0)
			cout<<classID<<"()"<<endl;
		
		//****************
		// CONTROL_PARAMS
		//****************
		NcontrolEqns = getIntParam("N_PSALC_CONTROL_PARAMS", "0");
		if(NcontrolEqns > 0) {
			lambda = new double [NcontrolEqns]; // lambda is the input parameter to the system (e.g. fuel equivalence ratio)
			char paramName[20];
			for(int i=0; i<NcontrolEqns; ++i) {
				sprintf(paramName, "CONTROL_PARAM%d", i);
				lambda[i] = getParam("LAMBDA_INITIAL_0") -> getDouble(paramName); // Note: lambda is same as phi in the model
			}

			if(mpi_rank == 0)
				cout<<endl
				    <<"> LAMBDA0 = "<<lambda[0]<<endl
				    <<endl;
		} else {
			lambda = NULL;
		}
		
		//****************
		// VAN_GOGH
		//****************
		// test number
		itest = 0;
		ntests = 0;

		for(int i=0; i<3; ++i) {
			xcvMin[i] =  ABSURDLY_BIG_NUMBER;
			xcvMax[i] = -ABSURDLY_BIG_NUMBER;
		}

		this->init();
	}
	
	/*
	 * destructor
	 */
	virtual ~VanGoghWithModels() {
		if(mpi_rank==0)
			cout<<"~"<<classID<<"()"<<endl;
		this->clear();
	}
	
protected:
	/*
	 * Method: init
	 * ------------
	 * Initialize VanGoghWithModels member variables
	 */
	void init() ;

	/*
	 * Method: clear()
	 * ---------------
	 *
	 */
	void clear() ;
	
public:
	/**********************
	 * OVERLOADED METHODS
	 **********************/
	/*
	 * Method: run
	 * -----------
	 * Original code = run() in JoeWithModels.h
	 */
	virtual void run();
	
	/****************************
	 * UTILITY METHODS
	 ****************************/
	/*
	 * Method: calcRMS
	 * -------------------
	 * Return the root mean square of a field variable (Option: 1. simple  2. volume averaged)
	 */
	double calcRMS(const double *vec, bool volumeAvg);
	
	/*
	 * Method: calcRMS3D
	 * -----------------
	 * Return the root mean square of a field variable (Option: 1. simple  2. volume averaged)
	 */
	double calcRMS3D(const double (*vec)[3], const int dim, bool volumeAvg);
	
	/*
	 * Method: findMinValue
	 * --------------------
	 * Find the minimun value of an array and return it
	 */
	double findMinValue(double* array, const int arrSize);

	/*
	 * Method: findMaxalue
	 * --------------------
	 * Find the maximum value of an array and return it
	 */
	double findMaxValue(double* array, const int arrSize);

	/*
	 * Method: calcWallDistanceFace
	 * ----------------------------
	 * calculate wall distance: minimum distance of cv to wall face by looping over all cv's and faces
	 */
	void calcWallDistanceFace(double *wd);
	
	/****************************
	 * VANGOGH SPECIFIC METHODS
	 ****************************/
	/*
	 * Method: findMatchedIndex
	 * ------------------------
	 *
	 */
	int findMatchedIndex(int& foundIndex, const int icv, const int nCVs, const double (*x_cv_eigen)[3], const double epsil);

	/*
	 * Method: readLinSystemMATLAB
	 * ---------------------
	 * Read eigenmodes and adjoint modes from "dumpedEigen_pt*****.bin"
	 */
	void readLinSystemMATLAB(const int ptOnSCurve);

	/*
	 * Method: readEigenPairsRembrandtSerial
	 * -------------------------------------
	 * Read data dumped by RembrandtWithModels::writeEigenPairsParallel()
	 */
	void readEigenPairsRembrandtSerial(string &filename, int &ptOnSCurve_file, double* lambda_file, const int NcontrolEqns);

	void readEigenPairsRembrandtSerial(const char filename[], int &ptOnSCurve_file, double* lambda_file, const int NcontrolEqns);

	/*
	 * Method: readPsALCdumpedDataSerial
	 * --------------------------
	 * Read field data from file
	 * Original code: IkeWithModels::readPsALCdumpedDataSerial()
	 */
	void readPsALCdumpedDataSerial(const char filename[], double* qVec, const int NcontrolEqns);

	/*
	 * Method: updateEigenVecTecplotNS
	 * --------------------------
	 * Update NS variabls (rho_Dct_1stReal, rho_Adj_1stReal, etc.) from DirectEvecs and AdjointEvecs.
	 * This is just for Tecplot output
	 */
	void updateEigenVecTecplotNS();

	/*
	 * Method: updateEigenVecTecplotScalars
	 * ------------------------------------
	 * Update scalar variabls (e.g. kine_Dct_1stReal, kine_Adj_1stReal, etc.) from DirectEvecs and AdjointEvecs.
	 * This is just for Tecplot output
	 */
	virtual void updateEigenVecTecplotScalars() {}
	
	/*
	 * Method: updateUnstableVecNS
	 * ---------------------------
	 * Update NS variabls (rho_unstable, etc.) from a qVec that has been generated by reading a IKE binary file.
	 * This is NOT ONLY for Tecplot output BUT ALSO linear analysis
	 */
	void updateUnstableVecNS(double* qVecTemp, const int nVarsl);

	/*
	 * Method: updateUnstableVecScalars
	 * --------------------------------
	 * Update scalar variabls (kine_unstable, etc.) from a qVec that has been generated by reading a IKE binary file.
	 * This is NOT ONLY for Tecplot output BUT ALSO linear analysis
	 */
	virtual void updateUnstableVecScalars(double* qVecTemp, const int nVars) {}

	/*
	 * Method: storeInitRestartNS
	 * --------------------------
	 * Store initial NS data in arrays (e.g. rho_init, etc.)
	 */
	void storeInitRestartNS();

	/*
	 * Method: storeInitRestartScalars
	 * -------------------------------
	 * Store initial scalar data in arrays (e.g. kine_init, etc.)
	 */
	virtual void storeInitRestartScalars() {}
	
	/****************************
	 * FUNCTIONS FOR PERTURBATION
	 ****************************/
	/*
	 * Method: getPerturbParams
	 * ----------------------
	 * Get the parameters for perturbations
	 */
	void getPerturbParams(PerturbParams &perturbParams);

	/*
	 * Method: perturbFieldNS
	 * ----------------------
	 * Perturb the initial field with random variables and apply filter to the perturbations
	 */
	void perturbFieldNS();

	/*
	 * Method: perturbFieldScalars
	 * ---------------------------
	 * Perturb the initial field with random variables and apply filter to the perturbations
	 */
	virtual void perturbFieldScalars() {};

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
	double perturbScalar(double* scalarArray, double* array_perturb, const char varName[], const bool applyClipping);

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
	double perturbVector(double (*vectorArray)[3], const int coord, double* array_perturb, const char varName[], const bool applyClipping);

	/****************************
	 * OPTIMAL METRIC
	 ****************************/
	/*
	 * Method: calcMetricsAdj
	 * ----------------------
	 * Calculate metrics using adjoint: Let phi=flow field, phi0=flow field on the unstable branch, y=least stable eigenvector, a=adjoint vector,
	 *                                           (phi-phi0)*a
	 *                                  metric = ------------
	 *                                                y*a
	 * Return: metricNumer = (phi-phi0)*a
	 *         metricDenom = y*a
	 */
	void calcMetricsAdj(double* metricNumer, double* metricDenom,
			const double *rho_unstable, const double (*rhou_unstable)[3], const double *rhoE_unstable, const double *kine_unstable, const double *omega_unstable,
			const double *rho_1stMode, const double (*rhou_1stMode)[3], const double *rhoE_1stMode, const double *kine_1stMode, const double *omega_1stMode,
			const double *rho_adj, const double (*rhou_adj)[3], const double *rhoE_adj, const double *kine_adj, const double *omega_adj);

	/****************************
	 * HOOK FUNCTIONS
	 ****************************/
	/*
	 * Method: writeQoIOnFile
	 * ----------------------
	 * Averaged density and Averaged Mach number
	 */
	virtual void writeQoIOnFile(const int itest, char filename[], bool rewrite) {}

	/*
	 * Method: writeResidualOnFile
	 * ---------------------------
	 * Averaged density and Averaged Mach number
	 */
	virtual void writeResidualOnFile(const int itest, char filename[], bool rewrite) {}
};

#endif
