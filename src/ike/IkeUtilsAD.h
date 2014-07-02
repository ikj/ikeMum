/*
 * IkeUtilsAD.h
 *
 *  Created on: Nov 3, 2011
 *      Author: ikj
 */

#ifndef IKEUTILSAD_H_
#define IKEUTILSAD_H_

#include "JoeWithModels.h" // this header file is required for MPI (e.g. writeMatrixBinary())
#include "adolc.h"
#include "MatComprsed.h"
#include <sstream>
#include <string>

/*
 * Function: isNaN
 * ---------------
 * Check if the given double variable is NaN
 */
inline bool isNaN(const double num) {
	return (num) != (num);
}

/*
 * Function: isInfinite
 * ---------------
 * Check if the given double variable is inf
 */
inline bool isInfinite(const double num) {
	return !((num+1)!=num && num==num);
}

/*
 * Function: isBigEndian()
 * -----------------------
 * Determine which type is the machine (Big endian or Little endian)
 */
bool isBigEndian();

/*
 * Function: convertEndian()
 * -------------------------
 * Convert between big-endian and little-endian values
 * There are two possible remedies found from the web:
 *   http://stackoverflow.com/questions/1001307/detecting-endianness-programmatically-in-a-c-program
 *   http://stackoverflow.com/questions/3823921/convert-big-endian-to-little-endian-when-reading-from-a-binary-file
 * Here, I adopt the second method which uses the reverse() function.
 */
template <class T>
void convertEndian(T *objp);

#ifndef MY_MEM_H
#include "../common/myMem.h"
///*
// * Function: getMem1D
// * ------------------
// * This function is actually same as the getMem1D() function defined in myMem.h
// */
//template <class VarType>
//inline void getMem1D(VarType **ptr, int i1, int i2, char *name, bool setZero = false) {
//    VarType *ar = new VarType[i2-i1+1];
//
//    if(ar == NULL)
//        throw(-1);
//
//    if(setZero)
//        memset(ar, 0, sizeof(VarType)*(i2-i1+1));
//
//    ar -= i1;
//    *ptr = ar;
//}
//
///*
// * Function: getMem2D
// * ------------------
// * This function is actually same as the getMem2D() function defined in myMem.h
// */
//template <class VarType>
//inline void getMem2D(VarType ***ptr, int i1, int i2, int j1, int j2, char *name, bool setZero = false) {
//    VarType **ar = new VarType*[i2-i1+1];
//
//    if(ar == NULL)
//    	throw(-1);
//
//    ar -= i1;
//
//    for( int i=i1; i<=i2; i++ )
//        getMem1D(&ar[i],j1,j2,name,setZero);
//
//    *ptr = ar;
//}
//
///*
// * Function: getMem3D
// * ------------------
// * This function is actually same as the getMem3D() function defined in myMem.h
// */
//template <class VarType>
//inline void getMem3D(VarType ****ptr, int i1, int i2, int j1, int j2, int k1, int k2, char *name, bool setZero = false) {
//	VarType ***ar = new VarType**[i2-i1+1];
//
//	if(ar == NULL)
//	throw(-1);
//
//	ar -= i1;
//
//	for( int i=i1; i<=i2; i++ )
//		getMem2D(&ar[i],j1,j2,k1,k2,name,setZero);
//
//	*ptr = ar;
//}
//
///*
// * Function: freeMem2D
// * -------------------
// * This function is actually same as the freeMem2D() function defined in myMem.h
// */
//template<class VarType>
//inline void freeMem2D(VarType **ptr, int i1, int i2, int j1, int j2) {
//	for( int i=i1; i<=i2; i++)
//		delete [] (ptr[i] + j1);
//
//	delete [] (ptr+i1);
//}
//
///*
// * Function: freeMem3D
// * -------------------
// * This function is actually same as the freeMem3D() function defined in myMem.h
// */
//template<class VarType>
//inline void freeMem3D(VarType ***ptr, int i1, int i2, int j1, int j2, int k1, int k2) {
//	for( int i=i1; i<=i2; i++) {
//		for( int j=j1; j<=j2; j++)
//			delete [] (ptr[i][j] + k1);
//
//		delete [] (ptr[i] + j1);
//	}
//	delete [] (ptr+i1);
//}
#endif

/*
 * Method: max
 * -----------
 * Compare adouble and double variable
 */
inline adouble max(adouble a, double b) {
	adouble c;
	c = (a.value()>b) ? a : b;
	return c;
}
inline adouble max(double a, adouble b) {
	adouble c;
	c = (a>b.value()) ? a : b;
	return c;
}
inline adouble max(adouble a, adouble b) {
	adouble c;
	c = (a.value()>b.value()) ? a : b;
	return c;
}

/*
 * Method: min
 * -----------
 * Compare adouble and double variable
 */
inline adouble min(adouble a, double b) {
	adouble c;
	c = (a.value()<b) ? a : b;
	return c;
}
inline adouble min(double a, adouble b) {
	adouble c;
	c = (a<b.value()) ? a : b;
	return c;
}
inline adouble min(adouble a, adouble b) {
	adouble c;
	c = (a.value()<b.value()) ? a : b;
	return c;
}

/*
 * Method: vecDotVec3dADOLC
 * ------------------------
 * Dot product of two vectors which have 3 elements
 * Original code: vecDotVec3d() in tc_vec3d.h
 */
inline adouble vecDotVec3dADOLC(const adouble *v1, const adouble *v2) {
	adouble dotproduct = 0.0;

	for (int i = 0; i < 3; i++)
		dotproduct += v1[i]*v2[i];

	return(dotproduct);
}

/*
 * Method: vecDotVec3dADOLC
 * ------------------------
 * Dot product of two vectors which have 3 elements
 * Original code: vecDotVec3d() in tc_vec3d.h
 *
 * Just to avoid the problem of "const"
 */
inline adouble vecDotVec3dADOLC(adouble *v1, adouble *v2, bool dummy) {
	adouble dotproduct = 0.0;

	for (int i = 0; i < 3; i++)
		dotproduct += v1[i]*v2[i];

	return(dotproduct);
}

/*
 * Method: vecDotVec3dADOLC
 * ------------------------
 * Dot product of two vectors which have 3 elements
 * Original code: vecDotVec3d() in tc_vec3d.h
 */
inline adouble vecDotVec3dADOLC(const double *v1, const adouble *v2) {
	adouble dotproduct = 0.0;

	for (int i = 0; i < 3; i++)
		dotproduct += v1[i]*v2[i];

	return(dotproduct);
}

/*
 * Method: vecDotVec3dADOLC
 * ------------------------
 * Dot product of two vectors which have 3 elements
 * Original code: vecDotVec3d() in tc_vec3d.h
 */
inline adouble vecDotVec3dADOLC(const adouble *v1, const double *v2) {
	adouble dotproduct = 0.0;

	for (int i = 0; i < 3; i++)
		dotproduct += v1[i]*v2[i];

	return(dotproduct);
}

/*
 * Method: vecDotVec3dADOLC
 * ------------------------
 * Dot product of two vectors which have 3 elements
 * Original code: vecDotVec3d() in tc_vec3d.h
 *
 * Just to avoid the problem of "const"
 */
inline adouble vecDotVec3dADOLC(adouble *v1, double *v2, bool dummy) {
	adouble dotproduct = 0.0;

	for (int i = 0; i < 3; i++)
		dotproduct += v1[i]*v2[i];

	return(dotproduct);
}

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
 * Method: calcPress
 * -----------------
 *
 */
inline adouble calcPress(adouble gamma, adouble rhoE, adouble *rhou, adouble rho, adouble kinecv) {
	adouble rhouSq = 0.0;
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

/*
 * Method: writeMatrixForMatlab
 * --------------------------------
 * write the matrix on files as a matlab format
 * Sometimes, the size of the matrix is too large to be saved on a single file.
 */
void writeMatrixForMatlab(string &filename, MatComprsed &jacMatrix, bool mfileFormat, const int mpi_rank, const int nScal, const int *cvora, const int *nbocv_v_global);
void writeMatrixForMatlab(string &filename, MatComprsed &jacMatrix, bool mfileFormat, const int mpi_rank, const int nScal, const int *cvora, const int *nbocv_v_global, const int NcontrolEqns, const int ncv_gg);


/*
 * Method: writeMatrixForMatlab
 * --------------------------------
 * write the matrix on files as a matlab format
 * Sometimes, the size of the matrix is too large to be saved on a single file.
 */
void writeMatrixForMatlab(string &filename, MatComprsedSTL &jacMatrixSTL, bool mfileFormat, const int mpi_rank, const int nScal, const int *cvora, const int *nbocv_v_global);
void writeMatrixForMatlab(string &filename, MatComprsedSTL &jacMatrixSTL, bool mfileFormat, const int mpi_rank, const int nScal, const int *cvora, const int *nbocv_v_global, const int NcontrolEqns, const int ncv_gg);

/*
 * Method: writeMatrixBinary
 * --------------------------------
 * write the matrix on binary files
 * Sometimes, the size of the matrix is too large to be saved on a single file.
 */
void writeMatrixBinary(string &filename, MatComprsed &jacMatrix, const int mpi_rank, const int nScal, const int *cvora, const int *nbocv_v_global);
void writeMatrixBinary(string &filename, MatComprsed &jacMatrix, const int mpi_rank, const int nScal, const int *cvora, const int *nbocv_v_global, const int NcontrolEqns, const int ncv_gg);

/*
 * Method: writeMatrixBinary
 * --------------------------------
 * write the matrix on binary files
 * Sometimes, the size of the matrix is too large to be saved on a single file.
 */
void writeMatrixBinary(string &filename, MatComprsedSTL &jacMatrixSTL, const int mpi_rank, const int nScal, const int *cvora, const int *nbocv_v_global);
void writeMatrixBinary(string &filename, MatComprsedSTL &jacMatrixSTL, const int mpi_rank, const int nScal, const int *cvora, const int *nbocv_v_global, const int NcontrolEqns, const int ncv_gg);

/*
 * Method: writeCVcoordForMatlab
 * -----------------------------
 * write the coordinates of CVs on a file (as a matlab format)
 */
void writeCVcoordForMatlab(string &filename, bool mfileFormat, const int ncv, double (*x_cv)[3]);

/*
 * Method: writeCVcoordBinary
 * -----------------------------
 * write the coordinates of CVs on a binary file
 */
void writeCVcoordBinary(string &filename, const int ncv, int ncvGlobal, double (*x_cv)[3]);

#endif /* IKEUTILSAD_H_ */
