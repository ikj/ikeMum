/*
 * IkeUtilsAD.cpp
 *
 *  Created on: Dec 18, 2012
 *      Author: ikj
 */

#include "IkeUtilsAD.h"

/*
 * Function: isBigEndian()
 * -----------------------
 * Determine which type is the machine (Big endian or Little endian)
 */
bool isBigEndian() {
	int num = 1;
	if(*(char *)&num == 1)
		return false; // little endian
	else
		return true; // big endian
}

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
void convertEndian(T *objp) {
	unsigned char *memp = reinterpret_cast<unsigned char*>(objp);
	std::reverse(memp, memp + sizeof(T));
}


/*
 * Method: writeMatrixForMatlab
 * --------------------------------
 * write the matrix on files as a matlab format
 * Sometimes, the size of the matrix is too large to be saved on a single file.
 */
void writeMatrixForMatlab(string &filename, MatComprsed &jacMatrix, bool mfileFormat, const int mpi_rank, const int nScal, const int *cvora, const int *cv_gl) {
	writeMatrixForMatlab(filename, jacMatrix, mfileFormat, mpi_rank, nScal, cvora, cv_gl, 0, 0);
}
void writeMatrixForMatlab(string &filename, MatComprsed &jacMatrix, bool mfileFormat, const int mpi_rank, const int nScal, const int *cvora, const int *cv_gl, const int NcontrolEqns, const int ncv_gg) {
	assert(cvora!=NULL && cv_gl!=NULL);
	int m = 5 + nScal;

	// get nCols, nRows, nnz
	int nCols = jacMatrix.get_nCols(); // Each CPU might have different number of columns
	int nRows = jacMatrix.get_nRows();
	int nnz = jacMatrix.get_nnz();

	// open the file
	FILE *fp;
	stringstream ss;
	string newFilename;
	if(mfileFormat)
		ss<<filename<<"_"<<mpi_rank<<".m";
	else
		ss<<filename<<"_"<<mpi_rank<<".txt";
	newFilename = ss.str();

	fp = fopen(newFilename.c_str(), "w");

	// write total number of global Rows, global Columns, and global number of non-zeros at the beginning of the file
	if(mfileFormat) {
		fprintf(fp, "function [nRows, nCols, nnz, A] = %s_%d() \n", filename.c_str(),mpi_rank);
		fprintf(fp, "nRows = %d; \n",nRows);
		fprintf(fp, "nCols = %d; \n",nCols);
		fprintf(fp, "nnz = %d; \n",nnz);

		fprintf(fp, "A = [ \n");
	} else {
		fprintf(fp, "%d; \n",nRows);
		fprintf(fp, "%d; \n",nCols);
		fprintf(fp, "%d; \n",nnz);
	}

	// write row index, column index, and values as a MATLAB sparse matrix format
	for(int i=0; i<nnz; ++i) {
		fprintf(fp, "%d\t%d\t%.8e \n", jacMatrix.get_global_rind(i, mpi_rank, m, cvora)+1, jacMatrix.get_global_cind(i, mpi_rank, m, cvora, cv_gl, NcontrolEqns, ncv_gg)+1, jacMatrix.get_values(i)); // MATLAB index starts from 1 instead of 0
	}

	if(mfileFormat) {
		fprintf(fp, "]; \n");
		fprintf(fp, "end \n\n");
	}

	// close the file
	fclose(fp);
}

/*
 * Method: writeMatrixForMatlab
 * --------------------------------
 * write the matrix on files as a matlab format
 * Sometimes, the size of the matrix is too large to be saved on a single file.
 */
void writeMatrixForMatlab(string &filename, MatComprsedSTL &jacMatrixSTL, bool mfileFormat, const int mpi_rank, const int nScal, const int *cvora, const int *cv_gl) {
	writeMatrixForMatlab(filename, jacMatrixSTL, mfileFormat, mpi_rank, nScal, cvora, cv_gl, 0, 0);
}
void writeMatrixForMatlab(string &filename, MatComprsedSTL &jacMatrixSTL, bool mfileFormat, const int mpi_rank, const int nScal, const int *cvora, const int *cv_gl, const int NcontrolEqns, const int ncv_gg) {
	assert(cvora!=NULL && cv_gl!=NULL);
	int m = 5 + nScal;

	// get nCols, nRows, nnz
	int nCols = jacMatrixSTL.get_nCols(); // Each CPU might have different number of columns
	int nRows = jacMatrixSTL.get_nRows();
	int nnz = jacMatrixSTL.get_nnz();

	// open the file
	FILE *fp;
	stringstream ss;
	string newFilename;
	if(mfileFormat)
		ss<<filename<<"_"<<mpi_rank<<".m";
	else
		ss<<filename<<"_"<<mpi_rank<<".txt";
	newFilename = ss.str();

	fp = fopen(newFilename.c_str(), "w");

	// write total number of global Rows, global Columns, and global number of non-zeros at the beginning of the file
	if(mfileFormat) {
		fprintf(fp, "function [nRows, nCols, nnz, A] = %s_%d() \n", filename.c_str(),mpi_rank);
		fprintf(fp, "nRows = %d; \n",nRows);
		fprintf(fp, "nCols = %d; \n",nCols);
		fprintf(fp, "nnz = %d; \n",nnz);

		fprintf(fp, "A = [ \n");
	} else {
		fprintf(fp, "%d; \n",nRows);
		fprintf(fp, "%d; \n",nCols);
		fprintf(fp, "%d; \n",nnz);
	}

	// write row index, column index, and values as a MATLAB sparse matrix format
	for(int i=0; i<nnz; ++i) {
		fprintf(fp, "%d\t%d\t%.8e \n", jacMatrixSTL.get_global_rind(i, mpi_rank, m, cvora)+1, jacMatrixSTL.get_global_cind(i, mpi_rank, m, cvora, cv_gl, NcontrolEqns, ncv_gg)+1, jacMatrixSTL.get_values(i)); // MATLAB index starts from 1 instead of 0
	}

	if(mfileFormat) {
		fprintf(fp, "]; \n");
		fprintf(fp, "end \n\n");
	}

	// close the file
	fclose(fp);
}

/*
 * Method: writeMatrixBinary
 * --------------------------------
 * write the matrix on binary files
 * Sometimes, the size of the matrix is too large to be saved on a single file.
 */
void writeMatrixBinary(string &filename, MatComprsed &jacMatrix, const int mpi_rank, const int nScal, const int *cvora, const int *cv_gl) {
	writeMatrixBinary(filename, jacMatrix, mpi_rank, nScal, cvora, cv_gl, 0, 0);
}
void writeMatrixBinary(string &filename, MatComprsed &jacMatrix, const int mpi_rank, const int nScal, const int *cvora, const int *cv_gl, const int NcontrolEqns, const int ncv_gg) {
	assert(cvora!=NULL && cv_gl!=NULL);
	int m = 5 + nScal;
	int ncv = cvora[mpi_rank+1] - cvora[mpi_rank];

	// open the file
	stringstream ss;
	ss<<filename<<"_"<<mpi_rank<<".bin";
	filename = ss.str();
	ofstream outfile(filename.c_str(),ofstream::binary);

	// write nRows, nCols, and nnz
	int nCols = jacMatrix.get_nCols(); // Each CPU might have different number of columns
	int nRows = jacMatrix.get_nRows();
	int nnz = jacMatrix.get_nnz();

	outfile.write(reinterpret_cast<char*>(&nRows),sizeof(int));
	outfile.write(reinterpret_cast<char*>(&nCols),sizeof(int));
	outfile.write(reinterpret_cast<char*>(&nnz),sizeof(int));

	// write row index, column index, and values as a MATLAB sparse matrix format
	for(int i=0; i<nnz; ++i) {
		int localRind = jacMatrix.get_rind(i);
		if(localRind < ncv) { // For some reason, jacMatrix can have more than ncv rows (actually ncv_gg rows)
			int globalRind = jacMatrix.get_global_rind(i, mpi_rank, m, cvora)+1;
			int globalCind = jacMatrix.get_global_cind(i, mpi_rank, m, cvora, cv_gl, NcontrolEqns, ncv_gg)+1;
			double value =  jacMatrix.get_values(i);
			outfile.write(reinterpret_cast<char*>(&globalRind),sizeof(int));
			outfile.write(reinterpret_cast<char*>(&globalCind),sizeof(int));
			outfile.write(reinterpret_cast<char*>(&value),sizeof(double));
		}
	}

	// close the file
	outfile.close();
}

/*
 * Method: writeMatrixBinary
 * --------------------------------
 * write the matrix on binary files
 * Sometimes, the size of the matrix is too large to be saved on a single file.
 */
void writeMatrixBinary(string &filename, MatComprsedSTL &jacMatrixSTL, const int mpi_rank, const int nScal, const int *cvora, const int *cv_gl) {
	writeMatrixBinary(filename, jacMatrixSTL, mpi_rank, nScal, cvora, cv_gl, 0, 0);
}
void writeMatrixBinary(string &filename, MatComprsedSTL &jacMatrixSTL, const int mpi_rank, const int nScal, const int *cvora, const int *cv_gl, const int NcontrolEqns, const int ncv_gg) {
	assert(cvora!=NULL && cv_gl!=NULL);
	int m = 5 + nScal;
	int ncv = cvora[mpi_rank+1] - cvora[mpi_rank];

	// open the file
	stringstream ss;
	ss<<filename<<"_"<<mpi_rank<<".bin";
	filename = ss.str();
	ofstream outfile(filename.c_str(),ofstream::binary);

	// write nRows, nCols, and nnz
	int nCols = jacMatrixSTL.get_nCols(); // Each CPU might have different number of columns
	int nRows = jacMatrixSTL.get_nRows();
	int nnz = jacMatrixSTL.get_nnz();

	outfile.write(reinterpret_cast<char*>(&nRows),sizeof(int));
	outfile.write(reinterpret_cast<char*>(&nCols),sizeof(int));
	outfile.write(reinterpret_cast<char*>(&nnz),sizeof(int));

	// write row index, column index, and values as a MATLAB sparse matrix format
	for(int i=0; i<nnz; ++i) {
		int localRind = jacMatrixSTL.get_rind(i);
		if(localRind < ncv) { // For some reason, jacMatrix can have more than ncv rows (actually ncv_gg rows)
			int globalRind = jacMatrixSTL.get_global_rind(i, mpi_rank, m, cvora)+1;
			int globalCind = jacMatrixSTL.get_global_cind(i, mpi_rank, m, cvora, cv_gl, NcontrolEqns, ncv_gg)+1;
			double value =  jacMatrixSTL.get_values(i);
			outfile.write(reinterpret_cast<char*>(&globalRind),sizeof(int));
			outfile.write(reinterpret_cast<char*>(&globalCind),sizeof(int));
			outfile.write(reinterpret_cast<char*>(&value),sizeof(double));
		}
	}

	// close the file
	outfile.close();
}

/*
 * Method: writeCVcoordForMatlab
 * -----------------------------
 * write the coordinates of CVs on a file (as a matlab format)
 */
void writeCVcoordForMatlab(string &filename, bool mfileFormat, const int ncv, double (*x_cv)[3]) {
	// get global numbering of CV from the previous CPUs
	int ncvSum_previousCore = 0;
	if(mpi_rank != 0) {
		MPI_Status status;
		MPI_Recv(&ncvSum_previousCore, 1, MPI_INT, mpi_rank-1, 820, mpi_comm, &status);
	}

	// open the file
	FILE *fp;
	stringstream ss;
	string newFilename;
	if(mfileFormat)
		ss<<filename<<".m";
	else
		ss<<filename<<".txt";

	newFilename = ss.str();
	if(mpi_rank==0)
		fp = fopen(newFilename.c_str(), "w");
	else
		fp = fopen(newFilename.c_str(), "a");

	if(mfileFormat)
		if(mpi_rank==0) {
			fprintf(fp, "function [CVs] = %s() \n", filename.c_str());
			fprintf(fp, "CVs = [ \n");
		}

	// write the coordinates
	for(int icv=0; icv<ncv; ++icv) {
		fprintf(fp, "%d\t%.8e\t%.8e\t%.8e \n", ncvSum_previousCore+icv, x_cv[icv][0], x_cv[icv][1], x_cv[icv][2]);
	}

	if(mfileFormat)
		if(mpi_rank==mpi_size-1) {
			fprintf(fp, "]; \n");
			fprintf(fp, "end \n\n");
		}

	// close the file
	fclose(fp);

	if(mpi_rank < mpi_size-1) {
		int myGlobalNcvSum = ncvSum_previousCore + ncv;
		MPI_Send(&myGlobalNcvSum, 1, MPI_INT, mpi_rank+1, 820, mpi_comm);
	}
}

/*
 * Method: writeCVcoordBinary
 * -----------------------------
 * write the coordinates of CVs on a binary file
 */
void writeCVcoordBinary(string &filename, const int ncv, int ncvGlobal, double (*x_cv)[3]) {
	// get global numbering of CV from the previous CPUs
	int ncvSum_previousCore = 0;
	if(mpi_rank != 0) {
		MPI_Status status;
		MPI_Recv(&ncvSum_previousCore, 1, MPI_INT, mpi_rank-1, 1111, mpi_comm, &status);
	}

	// open the file
	stringstream ss;
	ss<<filename<<".bin";
	filename = ss.str();
	ofstream outfile;
	if(mpi_rank==0)
		outfile.open(filename.c_str(), ofstream::trunc | ofstream::binary);
	else
		outfile.open(filename.c_str(), ofstream::app | ofstream::binary);

	// write the coordinates
	if(mpi_rank==0)
		outfile.write(reinterpret_cast<char*>(&ncvGlobal), sizeof(int));

	for(int icv=0; icv<ncv; ++icv) {
		int icvGlobal = ncvSum_previousCore+icv;
		outfile.write(reinterpret_cast<char*>(&icvGlobal),sizeof(int));
		outfile.write(reinterpret_cast<char*>(&x_cv[icv][0]),sizeof(double)*3);
	}

	// close the file
	outfile.close();

	if(mpi_rank < mpi_size-1) {
		int myGlobalNcvSum = ncvSum_previousCore + ncv; // actually, cvora has this information
		MPI_Send(&myGlobalNcvSum, 1, MPI_INT, mpi_rank+1, 1111, mpi_comm);
	}
}


