/*
 * DiffFilter.h
 *
 *  Created on: Jan 23, 2013
 *      Author: ikj
 */

#ifndef DIFFFILTER_H_
#define DIFFFILTER_H_

//#include "UgpWithCv2.h"
//#include "Param.h"
//#include "myMem.h"
#include "UgpWithCvCompFlow.h"

/*
 * Class: DiffFilter
 * -----------------
 * Applying differential filter (Germano: Phy.Fluid. 1986, letter) to a scalar (or a vector)
 * Note: The original code was developed by Sanjeeb Bose (stbose AT stanford DOT edu), 
 *       but this class is taken from Jeff O'Brien(jobrien AT stanford DOT edu)'s version.
 */
class DiffFilter : virtual public UgpWithCvCompFlow // virtual public UgpWithCv2Op, virtual public ParamMap
{
protected:
	double *cvDf; 
	
public:
	/*
	 * Constructor: standard constructor needed for Models that virtual inherits from UgpWithCvCompFlow
	 */
	DiffFilter() {
		if ( mpi_rank == 0 ) 
			cout << "DiffFilter() " << endl ; 
		init();  
	}

	/*
	 * Destructor
	 */
	virtual ~DiffFilter() {
		if ( mpi_rank == 0 ) 
			cout << "~DiffFilter() " << endl ; 
		clear();
	}
	
private:
	/*
	 * Method: init
	 * ------------
	 * 
	 */
	void init() {
		cvDf = NULL;
	}
	
	/*
	 * Method: clear
	 * -------------
	 * 
	 */
	void clear() {
		if(cvDf != NULL) {
			delete [] cvDf;
			cvDf = NULL;
		}
	}
	
public:
	/*
	 * Method: filterConsturcted
	 * -------------------------
	 *
	 */
	bool filterConsturcted();

	/*
	 * Method: appplyDiffFilterG1G2
	 * ----------------------------
	 *
	 */
	void applyDiffFilterG1G2(double *phif, const double *phi, const double absTol=1.0e-6, const int maxIter=1000);
	/*
	 * Method: appplyDiffFilterG1G2
	 * ----------------------------
	 *
	 */
	void applyDiffFilterG1G2(double (*phif)[3], const double (*phi)[3], const double absTol=1.0e-6, const int maxIter=1000);

	/*
	 * Method: specify_filter_width
	 * ----------------------------
	 *     delta^2
	 * p = -------, where delta = filter width, if p is constant everywhere
	 *      40.0
	 */
	virtual double specify_filter_width(const int ifa);

	/*
	 * Method: buildCvDifferentialFilter
	 * ---------------------------------
	 * 
	 */
	void buildCvDifferentialFilter();

protected:
	/*
	 * Method: solveCvScalarCgG1G2
	 * ---------------------------
	 * conjugate gradient
	 * Original code = solveCvScalarCg in UgpWithCv2.h (original version)
	 */
	int solveCvScalarCgG1G2(double * phi, const double *Ap, const double *rhs, const int mode, const double zero, const int maxiter);
};

#endif
