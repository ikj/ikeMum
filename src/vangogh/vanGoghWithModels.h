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
 * Check if the time-advancement scheme that you want to use has the option to be stopped with "resid_energ_th"
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

class Complex1Darray
{
private:
	std::complex<double>* array1D;
	int arraySize;

public:
	Complex1Darray() {
		array1D = NULL;
		arraySize = 0;
	}
	~Complex1Darray() {
		clear();
	}

	/*
	 * Copy Constructors and Assignment Operators
	 */
	Complex1Darray(const Complex1Darray& other) {
		copyOther(other);
	}
	Complex1Darray& operator = (const Complex1Darray& other) {
		if(this != &other) {
			clear();
			// Note: When we cover inheritance, there's one more step here.
			copyOther(other);
		}
		return *this;
	}

	/*
	 * Method: copyOther
	 * -----------------
	 *
	 */
	void copyOther(const Complex1Darray& other) {
		array1D = new std::complex<double> [other.size()];
		copy(other.begin(), other.end(), begin());
	}

	/*
	 * Method: allocate
	 * ----------------
	 *
	 */
	void allocate(const int n) {
		assert(array1D == NULL && arraySize == 0);
		arraySize = n;
		array1D = new std::complex<double> [arraySize];
	}

	/*
	 * Method: clear
	 * -------------
	 *
	 */
	void clear() {
		if(array1D != NULL) {
			assert(arraySize != 0);
			delete [] array1D; 	array1D = NULL;
		} else {
			assert(arraySize == 0);
		}
		arraySize = 0;
	}

	/*
	 * Method: empty
	 * -------------
	 *
	 */
	bool empty() const {
		if(array1D==NULL) {
			assert(arraySize == 0);
			return true;
		}
		return false;
	}

	/*
	 * Method: getAt
	 * -------------
	 *
	 */
	complex<double>& getAt(const int index) {
		if(index<0)
			cerr<<"Complex1Darray::getAt(): ERROR -- index="<<index<<" is smaller than 0"<<endl;
		if(index>=arraySize)
			cerr<<"Complex1Darray::getAt(): ERROR -- index="<<index<<" is greater than size()="<<arraySize<<endl;
		assert(index>=0 && index<arraySize);

		return array1D[index];
	}
	const complex<double>& getAt(const int index) const {
		if(index<0)
			cerr<<"Complex1Darray::getAt(): ERROR -- index="<<index<<" is smaller than 0"<<endl;
		if(index>=arraySize)
			cerr<<"Complex1Darray::getAt(): ERROR -- index="<<index<<" is greater than size()="<<arraySize<<endl;
		assert(index>=0 && index<arraySize);

		return array1D[index];
	}

	/*
	 * Method: size
	 * ------------
	 *
	 */
	int size() const {
		return arraySize;
	}

	/***********************************************************
	 * Element Selection Operator: []
	 * ------------------------------
	 * Two overloaded operators for the mutable and immutable cases
	 ***********************************************************/
	complex<double>& operator[] (const int i) {
		return getAt(i);
	}

	const complex<double>& operator[] (const int i) const {
		return getAt(i);
	}

	/***********
	 * Iterators
	 * ---------
	 * Sample code from http://www.dreamincode.net/forums/topic/58468-making-your-own-iterators/
	 * For example, if you want to iterate over the elements of "ADscalar<elemType> arr" and print them out:
	 *   for(ADscalar<elemType>::iterator iter=arr.begin(); iter!=arr.end(); ++iter) {
	 *   	printf("  %.2e", *iter); // Note: for the "adouble" elemType: printf("  %.2e", iter->value());
	 *   }
	 ***********/
	class iterator
	{
	public:
		iterator(complex<double>* ptr) : ptr_(ptr) { }
		iterator operator ++() { iterator i = *this; ptr_++; return i; }
		iterator operator ++(int junk) { ptr_++; return *this; }
		complex<double>& operator *() { return *ptr_; }
		complex<double>* operator ->() { return ptr_; }
		bool operator ==(const iterator& rhs) { return ptr_ == rhs.ptr_; }
		bool operator !=(const iterator& rhs) { return ptr_ != rhs.ptr_; }
	private:
		complex<double>* ptr_;
	};

	class const_iterator
	{
	public:
		const_iterator(complex<double>* ptr) : ptr_(ptr) { }
		const_iterator operator ++() { const_iterator i = *this; ptr_++; return i; }
		const_iterator operator ++(int junk) { ptr_++; return *this; }
		const complex<double>& operator *() { return *ptr_; }
		const complex<double>* operator ->() { return ptr_; }
		bool operator ==(const const_iterator& rhs) { return ptr_ == rhs.ptr_; }
		bool operator !=(const const_iterator& rhs) { return ptr_ != rhs.ptr_; }
	private:
		complex<double>* ptr_;
	};

	friend class iterator;
	friend class const_iterator;

	iterator begin() {
		return iterator(array1D);
	}

	iterator end() {
		return iterator(array1D + this->size());
	}

	const_iterator begin() const {
		return const_iterator(array1D);
	}

	const_iterator end() const {
		return const_iterator(array1D + this->size());
	}

	/*
	 * Method: copy
	 * ------------
	 *
	 */
	iterator copy(const_iterator first, const_iterator last, iterator result) {
		while(first!=last) {
			*result = *first;
			++result;
			++first;
		}
		return result;
	}
};

//class ComplexTensors {
//
//};

/*
 * Method: calcPress
 * -----------------
 *
 */
inline double calcPress(double gamma, double rhoE, double *rhou, double rho, double kinecv) {
	double rhouSq = 0.0;
	for(int i=0; i<3; ++i)
		rhouSq += rhou[i]*rhou[i];
	return (gamma - 1.0) * (rhoE - 0.5 * rhouSq/rho - rho * kinecv);
}

/*
 * Method: calcTemp
 * ----------------
 *
 */
inline double calcTemp(double press, double rho, double RoM) {
	return press / (rho * RoM);
}

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
	// NOTE: scalars should be registered in the child class
	
	// data on the unstable(hyperbolic) branch (from Q0_PT000**.bin)
	double *rho_unstable;
	double (*rhou_unstable)[3];
	double *rhoE_unstable;
	// NOTE: scalars should be registered in the child class

	// -------------------------
	// Perturbations
	// -------------------------
	// perturbations
	double *array_perturb;

	double *rho_perturb;
	double (*rhou_perturb)[3];
	double *rhoE_perturb;
	// NOTE: scalars should be registered in the child class

	// filter
	PerturbParams perturbParams;
	double *wdFace;

	// -------------------------
	// Eigenvals: Data container
	// -------------------------
	Complex1Darray Evals;   // Note: A data container that has the 'std::complex' class as its element

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
	// NOTE: scalars should be registered in the child class

	// first adjoint global modes
	double *rho_Adj_1stReal;
	double (*rhou_Adj_1stReal)[3];
	double *rhoE_Adj_1stReal;
	// NOTE: scalars should be registered in the child class
	
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
	virtual void updateEigenVecTecplotScalars() { /* empty */ }
	
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
	virtual void updateUnstableVecScalars(double* qVecTemp, const int nVars) { /* empty */ }

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
	virtual void storeInitRestartScalars() { /* empty */ }
	
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
	 * Method: reinitialHook
	 * ---------------------
	 * Similar role to initialHook(): Reinitialize the flow field from the previous simulation.
	 */
	virtual void reinitialHook();

	/*
	 * Method: reinitialHookScalarRansTurbModel
	 * ----------------------------------------
	 * Similar role to initialHook(): Reinitialize the flow field from the previous simulation.
	 */
	virtual void reinitialHookScalarRansTurbModel() { /* empty */ }

	/*
	 * Method: reinitialHookScalarRansCombModel
	 * ----------------------------------------
	 * Similar role to initialHook(): Reinitialize the flow field from the previous simulation.
	 */
	virtual void reinitialHookScalarRansCombModel() { /* empty */ }

	/*
	 * Method: perturbFieldNS
	 * ----------------------
	 * Perturb the initial field with random variables and apply filter to the perturbations
	 */
	virtual void perturbFieldNS();

	/*
	 * Method: perturbFieldScalars
	 * ---------------------------
	 * Perturb the initial field with random variables and apply filter to the perturbations
	 */
	virtual void perturbFieldScalars() { /* empty */ }

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
	 *           array_perturb[icv] = coeff * RANDOM_VARIABLE
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
	 * VANGOGH OVERLOADED METHODS
	 ****************************/
	/*
	 * Method: writeQoIOnFile
	 * ----------------------
	 * Averaged density and Averaged Mach number
	 */
	virtual void writeQoIOnFile(const int itest, char filename[], bool rewrite) { /* empty */ }

	/*
	 * Method: writeResidualOnFile
	 * ---------------------------
	 *
	 */
	virtual void writeResidualOnFile(const int itest, char filename[], bool rewrite);
};

#endif
