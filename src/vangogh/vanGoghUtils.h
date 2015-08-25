/*
 * vanGoghUtils.h
 *
 *  Created on: Jan 9, 2013
 *      Author: ikj
 */

#ifndef VANGOGHUTILS_H_
#define VANGOGHUTILS_H_

#include <cstdlib>
#include <ctime>
#include <math.h>
#ifndef PI
	#define PI 3.14159265358979323846
#endif

/*
 * Function: randArrayUniform
 * --------------------------
 * Generate an uniformly distributed random array whose values are in the interval of [0,1]
 */
inline void randArrayUniform(double* array, const int size) {
	// seed random variables
	srand((unsigned)time(NULL));

	// generate random variables
	for(int index=0; index<size; index++){
		array[index] = (double) rand()/RAND_MAX;
	}
}

/*
 * Function: randArrayUniformPositive
 * ----------------------------------
 * Generate an uniformly distributed random array whose values are in the interval of (0,1]
 */
inline void randArrayUniformPositive(double* array, const int size) {
	// seed random variables
	srand((unsigned)time(NULL));

	// generate random variables
	for(int index=0; index<size; index++){
		int dummy = rand();
		while(dummy == 0) {
			dummy = rand();
		}
		array[index] = (double) dummy/RAND_MAX;
	}
}

/*
 * Function: randArrayNormal
 * -------------------------
 * Generate a normally distributed (Gaussian) random array by using the standard Box–Muller transform
 */
inline void randArrayNormal(double* array, const int size) {
	// Get two uniformly distributed random arrays in the interval of (0,1]
	double* unif1 = new double [size];
	double* unif2 = new double [size];
	randArrayUniformPositive(unif1, size);
	randArrayUniformPositive(unif2, size);

	// Generate a normally distributed array by using the standard Box–Muller transform
	for(int index=0; index<size; index++){
		double R = sqrt(-2.0*log(unif1[index]));
		double theta = 2*PI*unif2[index];

		array[index] = R*cos(theta); // Box–Muller transform
	}

	delete [] unif1;
	delete [] unif2;
}

/*
 * Function: randArrayNormalBounded
 * --------------------------------
 * Generate a normally distributed (Gaussian) random array by using the standard Box–Muller transform
 */
inline void randArrayNormalBounded(double* array, const int size, const double minVal, const double maxVal) {
	// seed random variables
	srand((unsigned)time(NULL));

	// Generate a normally distributed array by using the standard Box–Muller transform
	for(int index=0; index<size; index++){
		bool notBounded = true;
		double dummyDouble;
		while(notBounded) {
			// Get two uniformly distributed random variables in the interval of (0,1]
			double unif1, unif2;

			int dummyInt = rand();
			while(dummyInt == 0) {
				dummyInt = rand();
			}
			unif1 = (double) dummyInt/RAND_MAX;

			dummyInt = rand();
			while(dummyInt == 0) {
				dummyInt = rand();
			}
			unif2 = (double) dummyInt/RAND_MAX;

			// Box–Muller transform
			double R = sqrt(-2.0*log(unif1));
			double theta = 2*PI*unif2;

			dummyDouble = R*cos(theta);

			// Check if the value is in the range
			if(dummyDouble>=minVal && dummyDouble<=maxVal)
				notBounded = false;
		}
		array[index] = dummyDouble;
	}
}

#endif /* VANGOGHUTILS_H_ */
