/*
 * rubensWithModels.h
 *
 *  Created on: Sep, 2015
 *      Author: ikj
 */

#ifndef RUBENSWITHMODELS_H_
#define RUBENSWITHMODELS_H_

#include "JoeWithModels.h"
#include "myMem.h"

#include <math.h>
#include <sstream>
#include <iomanip>

#ifndef MPI_WTIME_IS_GLOBAL
#define MPI_WTIME_IS_GLOBAL // MPI_Wtime() is synchronized across all processes in MPI_COMM_WORLD.
#endif

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

#define MAX_Q1_FILENAME_FULL_LENGTH 200

#define RUBENS_ERROR_CODE -5

#define RUBENS_DEBUG_LEVEL 1

//// TO DO: MATLAB version of file read!!!!
//#define NVARS_EIGEN 7
//#define NEIGENS 30

/*
 * Note:
 * You must make JoeWithModels::run() virtual
 */

// ###########################################################################################
// ------                                                                               ------
// ------                 Rubens continuously reads binary Q1 files                     ------
// ------                                                                               ------
// ###########################################################################################
class RubensWithModels : virtual public JoeWithModels {
protected:
	int NcontrolEqns;
	double *lambda;
	
	string classID;

	//**************
	// Rubens
	//**************
	// point number
	int ipt;
	int init_pt;
	int final_pt;

	// directory of the Q1 files
	string PathQ1dir;

	// Q1 filename
	char filenameQ1[MAX_Q1_FILENAME_FULL_LENGTH];

	// q vector (the data from the Q1 file is stored in this array)
	double* qVec;

	// RHS arrays (for rhs calculation and Tecplot output)
	double *RHSrho;
	double (*RHSrhou)[3];
	double *RHSrhoE;
	double **RHSrhoScal;

	// xcvMin & xcvMax (will be used while reading Q1.bin file from IKE file)
	double xcvMin[3];
	double xcvMax[3];
	
	// QoI filename
	char filenameQoI[20];

public:
	/*
	 * constructor
	 */
	RubensWithModels(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name) {
		classID = "RubensWithModels";
		if(mpi_rank==0)
			cout<<classID<<"()"<<endl;
		
		//****************
		// CONTROL_PARAMS
		//****************
		NcontrolEqns = getIntParam("N_PSALC_CONTROL_PARAMS", "0");
		if(NcontrolEqns > 0) {
			lambda = new double [NcontrolEqns]; // lambda is the input parameter to the system (e.g. fuel equivalence ratio)
		} else {
			lambda = NULL;
		}
		
		//****************
		// RUBENS
		//****************
		// point number
		ipt = -1;
		init_pt  = -1;
		final_pt = -1;

		// xcvMin & xcvMax (will be used while reading Q1.bin file from IKE file)
		for(int i=0; i<3; ++i) {
			xcvMin[i] =  ABSURDLY_BIG_NUMBER;
			xcvMax[i] = -ABSURDLY_BIG_NUMBER;
		}

		// QoI filename
		sprintf(filenameQoI, "QoI.txt");

		this->init();
	}
	
	/*
	 * destructor
	 */
	virtual ~RubensWithModels() {
		if(mpi_rank==0)
			cout<<"~"<<classID<<"()"<<endl;
		this->clear();
	}
	
protected:
	/*
	 * Method: init
	 * ------------
	 * Initialize RubensWithModels member variables
	 */
	void init();

	/*
	 * Method: clear()
	 * ---------------
	 *
	 */
	void clear();
	
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

	/**********************
	 * BUG-FIX in JoeWithModels
	 **********************/
	/*
	 * Method: setBC1D
	 * ---------------
	 * Original code = setBC1D in IkeWithModels_AD.cpp and setBC() in JoeWithModels.cpp
	 * In some situations, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
	 *
	 * CAUTION: For now, this method is empty since calculating correct RHS is not important in Rubens!!
	 */
	void setBC1D(const int ifa, double *rho, double (*rhou)[3], double *rhoE);

	/****************************
	 * RUBENS SPECIFIC METHODS
	 ****************************/
	/*
	 * Method: findMatchedIndex
	 * ------------------------
	 *
	 */
	int findMatchedIndex(int& foundIndex, const int icv, const int nCVs, const double (*x_cv_eigen)[3], const double epsil);

	/*
	 * Method: readPsALCdumpedDataSerial
	 * --------------------------
	 * Read field data from file
	 * Original code: IkeWithModels::readPsALCdumpedDataSerial()
	 */
	void readPsALCdumpedDataSerial(const char filename[], int &stepFromFile, double* lambda1FromFile, double* qVec, const int NcontrolEqns);

	/****************************
	 * RUBENS OVERLOADED METHODS
	 ****************************/
	/*
	 * Method: writeQoIOnFile
	 * ----------------------
	 * 
	 */
	virtual void writeQoIOnFile(const int itest, char filename[], bool rewrite) { /* empty */ }
};

#endif
