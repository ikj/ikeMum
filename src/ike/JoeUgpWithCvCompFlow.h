/*
 * JoeUgpWithCvCompFlow.h
 *
 *  Created on: Sep 3, 2014
 *      Author: ikj
 */

// This file will be included in JOE/ADJOINT_FILES/UgpWithCvCompFlowAD.h

#ifndef JOEUGPWITHCVCOMPFLOW_H_
#define JOEUGPWITHCVCOMPFLOW_H_

#ifndef EOF_ERROR_CHECK_CODE
#define EOF_ERROR_CHECK_CODE 503
#endif

#ifndef ABSURDLY_BIG_NUMBER
#define ABSURDLY_BIG_NUMBER 2.22e22
#endif

// =========================
//  UTILITY FUNCTIONS
// =========================

/*
 * Method: showBarrieredMessage
 * ----------------------------
 * One CPU core shows its message (usually, error or warning message).
 * MPI_Barrier() function will be called both at the beginning and in the end.
 *
 * Note: 1. All the CPU cores must call this method (Otherwise, communication error)
 *       2. If the string is empty, NO message will be shown
 */
void showBarrieredMessage(const char message[], const int mpi_rank_to_show) {
	MPI_Barrier(mpi_comm);
	for(int i=0; i<100000; ++i) {} // Pause for a while
	if(mpi_rank == mpi_rank_to_show)
		cout<<message<<endl;
	MPI_Barrier(mpi_comm);
}
void showBarrieredMessage(string& message, const int mpi_rank_to_show) {
	MPI_Barrier(mpi_comm);
	for(int i=0; i<100000; ++i) {} // Pause for a while
	if(mpi_rank == mpi_rank_to_show)
		cout<<message<<endl;
	MPI_Barrier(mpi_comm);
}

void writeJOEDataParallel(char filename[]) {
	string funcID = "UgpWithCvCompFlow::writeJOEDataParallel";

	// Number of variables
	int nScal = scalarTranspEqVector.size();
	int m = nScal+5; // number of variables (rho, rhou, rhov, rhow, rhoE, and scalars)

	if(mpi_rank==0)
		cout<<funcID<<"(): write data on "<<filename<<endl;

	// xcvMin and xcvMax
	double xcvMin[3] = {2.2e22, 2.2e22, 2.2e22}, 	xcvMax[3] = {-2.2e22, -2.2e22, -2.2e22};
	double *xcvMinArray = new double [mpi_size*3];
	double *xcvMaxArray = new double [mpi_size*3];
	for(int icv=0; icv<ncv; ++icv) {
		for(int i=0; i<3; ++i) {
			xcvMin[i] = min(xcvMin[i], x_cv[icv][i]);
			xcvMax[i] = max(xcvMax[i], x_cv[icv][i]);
		}
	}
	MPI_Gather(xcvMin, 3, MPI_DOUBLE, xcvMinArray, 3, MPI_DOUBLE, 0, mpi_comm);
	MPI_Gather(xcvMax, 3, MPI_DOUBLE, xcvMaxArray, 3, MPI_DOUBLE, 0, mpi_comm);
	MPI_Barrier(mpi_comm);

	/***********
	 ** Write the data
	 ***********/
	double wtime0, wtimeF;
	if(mpi_rank==0)
		wtime0 = MPI_Wtime();

	// 1. Header
	//      Structure of the header part:
	//              1. mpi_size (int)                  2. cvora (int*(mpi_size+1))
	//              3. xMinArray (double*3*mpi_size)   4. xMaxArray (double*3*mpi_size)
	//              5. nScal (int)
	if(mpi_rank==0) {
		// Open the file
		ofstream ofile;
		ofile.open(filename, ios_base::out | ios_base::trunc | ios_base::binary);

		// Write on the file
		int dummyInt;
		double dummyDouble;

		dummyInt=mpi_size; 			ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));

		for(int irank=0; irank<mpi_size+1; ++irank) {
			dummyInt=cvora[irank]; 		ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		}

		for(int irank=0; irank<mpi_size; ++irank) {
			for(int i=0; i<3; ++i) {
				dummyDouble=xcvMinArray[irank*3+i]; 	ofile.write(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
			}
		}
		for(int irank=0; irank<mpi_size; ++irank) {
			for(int i=0; i<3; ++i) {
				dummyDouble=xcvMaxArray[irank*3+i]; 	ofile.write(reinterpret_cast<char*>(&dummyDouble), sizeof(double));
			}
		}

		dummyInt=nScal; 			ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));

		// Close the file
		ofile.close();
	}
	MPI_Barrier(mpi_comm);
	int initDisp = sizeof(int) * (3 + mpi_size) + sizeof(double) * (2*mpi_size*3);

	// 2. Body
	//       Structure of the body part:
	//            For each mpi,
	//              1. x_cv (double*3*ncv)
	//              2. data (double*(5+nScal)*ncv)
	MPI_Status status;
	MPI_File fh;
	MPI_Offset displacement;

	// Open the file
	if (MPI_File_open(mpi_comm, filename, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh)!=0) {
			// Note: MPI_MODE_WRONLY = write only
			//       MPI_MODE_CREATE = create the file if it does not exist.
		cerr<<"ERROR "<<funcID<<"(): Cannot open "<<filename<<endl;
		throw(-1);
	}

	// Write the CV coordinates first
	int myOffsetNcv;
	MPI_Scan(&ncv, &myOffsetNcv, 1, MPI_INT, MPI_SUM, mpi_comm);
	myOffsetNcv -= ncv;
	displacement = MPI_Offset(initDisp + sizeof(double)*myOffsetNcv*(3+m));
	if(MPI_File_set_view(fh, displacement, MPI_INT, MPI_INT, "native", MPI_INFO_NULL)!=0) {
		// Note: native     = Data in this representation are stored in a file exactly as it is in memory
		//       Internal   = Data in this representation can be used for I/O operations in a homogeneous or heterogeneous environment
		//       External32 = This data representation states that read and write operations convert all data from and to the "external32" representation
		cerr<<"ERROR! "<<funcID<<"(): Cannot set the first MPI_File_set_view -- offset="<<displacement<<endl;
		throw(-1);
	}
	double* bufferDouble = new double [ncv*3];
	for(int icv=0; icv<ncv; ++icv) {
		int indexStart = icv*3;
		for(int i=0; i<3; ++i)
			bufferDouble[indexStart+i] = x_cv[icv][i];
	}
	MPI_File_write(fh, bufferDouble, ncv*3, MPI_DOUBLE, &status);
	delete [] bufferDouble; 	bufferDouble = NULL;

	// Write the flow data
	displacement += MPI_Offset(3*ncv*sizeof(double));
	bufferDouble = new double [ncv*m];
	for(int icv=0; icv<ncv; ++icv) {
		int indexStart = icv*m;
		bufferDouble[indexStart] = rho[icv];
		for(int i=0; i<3; ++i)
			bufferDouble[indexStart+1+i] = rhou[icv][i];
		bufferDouble[indexStart+4] = rhoE[icv];
		for(int iScal=0; iScal<nScal; ++iScal) {
			double *phi = scalarTranspEqVector[iScal].phi;
			bufferDouble[indexStart+5+iScal] = phi[icv];
		}
	}
	MPI_File_write(fh, bufferDouble, ncv*m, MPI_DOUBLE, &status);
	delete [] bufferDouble; 	bufferDouble = NULL;

	// 2-4. Close the file
	MPI_File_close(&fh);
	MPI_Barrier(mpi_comm);

	/*
	 * 3. Write the foot
	 */
	if(mpi_rank == mpi_size-1) {
		ofstream ofile;
		ofile.open(filename, ios_base::out | ios_base::app | ios_base::binary);
		int dummyInt = EOF_ERROR_CHECK_CODE; 	ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		ofile.close();
	}
	MPI_Barrier(mpi_comm);

	if(mpi_rank==0)
		wtimeF = MPI_Wtime();

	/***********
	 ** Show the summary on the screen if the debugging level is high
	 ***********/
	if(mpi_rank==0)
		cout<<"> The flow data is written on "<<filename<<" (RUN TIME = "<<wtimeF-wtime0<<" [sec]) \n"<<endl;

	/***********
	 ** Free the memory
	 ***********/
	delete [] xcvMinArray;
	delete [] xcvMaxArray;
}

/*
 * Method: writeJoeUncoupledNSMatrixClearDiagBinaryParallel
 * --------------------------------------------------------
 * Similar to writeJoeUncoupledNSMatrixBinaryParallel.
 * The difference is
 *   Extracting  tmp=cv_volume/local_dt  from the diagonal of A
 * Arguments:
 *   cvora          = number of ncv's at each mpi_rank
 *   nbocv_i        = indices information of nbocv_v. This can be also used to find the number of neighbors of a given icv
 *   nbocv_v_gloval = global neighboring icv information of a given icv. To access this array, nbocv_i is required.
 */
void writeJoeUncoupledNSMatrixClearDiagBinaryParallel(const char filename[], double (*A)[5][5], const int *cvora,
		const int ncv, const int *nbocv_i, const int *nbocv_v_golbal,
		const double *cv_volume, const double *local_dt) {
	string filenameString(filename);
	writeJoeUncoupledNSMatrixClearDiagBinaryParallel(filenameString, A, cvora, ncv, nbocv_i, nbocv_v_golbal,
			cv_volume, local_dt);
}

void writeJoeUncoupledNSMatrixClearDiagBinaryParallel(string &filename, double (*A)[5][5], const int *cvora,
		const int ncv, const int *nbocv_i, const int *nbocv_v_golbal,
		const double *cv_volume, const double *local_dt) {
	string funcID = "JoeUgpWithCvCompFlow.h::writeJoeUncoupledNSMatrixClearDiagBinaryParallel()";

	if(A == NULL) {
		if(mpi_rank==0)
			cerr<<"ERROR "<<funcID<<": The given matrix A is NULL"<<endl;
		throw(-1);
	}
	if(cvora == NULL || ncv == 0) {
		if(mpi_rank==0)
			cerr<<"ERROR "<<funcID<<": The given CV information (cvora, ncv) is NULL"<<endl;
		throw(-1);
	}
	if(nbocv_i == NULL || nbocv_v_golbal == NULL) {
		if(mpi_rank==0)
			cerr<<"ERROR "<<funcID<<": The given neighbor information (nbocv_i, nbocv_v_golbal) is NULL"<<endl;
		throw(-1);
	}

	char *filenameArray = new char[filename.size() + 1]; // Since MPI_File_open() doesn't accept "const char *", you cannot use std::string.c_str()
	std::copy(filename.begin(), filename.end(), filenameArray);
	filenameArray[filename.size()] = '\0';

	/***********
	 ** Initialization: Get nRows, nCols, nnz, localNnzArray
	 ***********/
	int nRowsPerBlock = 5;
	int nColsPerBlock = 5;
	int nbocv_s = nbocv_i[ncv]; // Total number of non-zero blocks in a mpi_rank

	int myNRows = ncv     * nRowsPerBlock;
	int myNnz   = nbocv_s * nRowsPerBlock*nColsPerBlock;

	int nRows;
	MPI_Allreduce(&myNRows, &nRows, 1, MPI_INT, MPI_SUM, mpi_comm);
	int nnz;
	MPI_Allreduce(&myNnz,   &nnz,   1, MPI_INT, MPI_SUM, mpi_comm);

	int nCols = cvora[mpi_size]*nColsPerBlock;

	int *localNnzArray = NULL;
	if(mpi_rank == 0)
		localNnzArray = new int [mpi_size];
	MPI_Gather(&myNnz, 1, MPI_INT, localNnzArray, 1, MPI_INT, 0, mpi_comm);

	/***********
	 ** Write the data
	 ***********/
	double wtime0, wtimeF;
	if(mpi_rank==0)
		wtime0 = MPI_Wtime();

	/*
	 * 1. Write the header part
	 */
	ofstream ofile;
	int dummyInt;

	if(mpi_rank == 0) {
		ofile.open(filenameArray, ios_base::out | ios_base::trunc | ios_base::binary);

		// Write the header (number of global rows, number of global columns, global nnz, local nnz array)
		dummyInt = nRows; 					ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int)); // Note: if you simply use the "<<" operator, it will cause some problem due to bugs in a g++ library on linux
		dummyInt = nCols; 					ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		dummyInt = nnz; 					ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		dummyInt = mpi_size; 				ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		for(int i=0; i<mpi_size; ++i) {
			dummyInt = localNnzArray[i]; 	ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		}
		ofile.close();
	}
	MPI_Barrier(mpi_comm);

	int initDisp = sizeof(int)*(4 + mpi_size);  // Initial displacement due to the header

	/*
	 * 2. Write the body (matrix)
	 */
	MPI_Status status;
	MPI_File   fh;
	MPI_Offset displacement;

	// Open the file
	if (MPI_File_open(mpi_comm, filenameArray, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh)!=0) {
			// Note: MPI_MODE_WRONLY = write only
			//       MPI_MODE_CREATE = create the file if it does not exist.
		stringstream ss;
		ss<<"ERROR! "<<funcID<<": Cannot open "<<filenameArray;
		string errorMessage = ss.str();

		if(mpi_rank==0) cerr<<errorMessage;

		delete [] filenameArray;
		if(localNnzArray != NULL) { delete [] localNnzArray; 	localNnzArray = NULL; }
		throw(-1);
	}

	// 2-1. Write the row indices first
	int myOffsetNnz;
	MPI_Scan(&myNnz, &myOffsetNnz, 1, MPI_INT, MPI_SUM, mpi_comm);
	myOffsetNnz -= myNnz;
	displacement = MPI_Offset(initDisp + sizeof(int)*myOffsetNnz);
	if(MPI_File_set_view(fh, displacement, MPI_INT, MPI_INT, "native", MPI_INFO_NULL)!=0) {
		// Note: native     = Data in this representation are stored in a file exactly as it is in memory
		//       Internal   = Data in this representation can be used for I/O operations in a homogeneous or heterogeneous environment
		//       External32 = This data representation states that read and write operations convert all data from and to the "external32" representation
		stringstream ss;
		ss<<"ERROR! "<<funcID<<": Cannot open "<<filenameArray;
		string errorMessage = ss.str();

		if(mpi_rank==0) cerr<<errorMessage;

		delete [] filenameArray;
		if(localNnzArray != NULL) { delete [] localNnzArray; 	localNnzArray = NULL; }
		throw(-1);
	}

	int *bufferInt = new int [myNnz];
	int count = 0;
	for(int icv=0; icv<ncv; ++icv) {
	    int noc_f = nbocv_i[icv];
	    int noc_l = nbocv_i[icv+1]-1;

	    for(int noc = noc_f; noc <= noc_l; ++noc) {
			for(int i=0; i<nRowsPerBlock; ++i) {
				int currGlobalRow = (cvora[mpi_rank] + icv)*nRowsPerBlock + i;
				for(int j=0; j<nColsPerBlock; ++j) {
					bufferInt[count] = currGlobalRow;
					++count;
				}
			}
	    }
	}
	assert(count == myNnz);
	MPI_File_write(fh, bufferInt, myNnz, MPI_INT, &status);
	delete [] bufferInt; 	bufferInt = NULL;

	// 2-2. Write the column indices
	displacement += MPI_Offset(nnz*sizeof(int));
	if(MPI_File_set_view(fh, displacement, MPI_INT, MPI_INT, "native", MPI_INFO_NULL)!=0) {
		stringstream ss;
		ss<<"ERROR! "<<funcID<<": Cannot set the second MPI_File_set_view -- offset="<<displacement<<endl;
		string errorMessage = ss.str();

		cerr<<errorMessage;

		delete [] filenameArray;
		if(localNnzArray != NULL) { delete [] localNnzArray; 	localNnzArray = NULL; }
		throw(-1);
	}

	bufferInt = new int [myNnz];
	count = 0;
	for(int icv=0; icv<ncv; ++icv) {
	    int noc_f = nbocv_i[icv];
	    int noc_l = nbocv_i[icv+1]-1;

	    for(int noc = noc_f; noc <= noc_l; ++noc) {
	    	int icvNbr = nbocv_v[noc];
			for(int i=0; i<nRowsPerBlock; ++i) {
				for(int j=0; j<nColsPerBlock; ++j) {
					int currGlobalCol = nbocv_v_golbal[icvNbr]*nColsPerBlock + j;
					assert(currGlobalCol >= 0 && currGlobalCol <= nCols);
					bufferInt[count] = currGlobalCol;
					++count;
				}
			}
	    }
	}
	MPI_File_write(fh, bufferInt, myNnz, MPI_INT, &status);
	delete [] bufferInt; 	bufferInt = NULL;

	// 2-3. Write the values
	displacement += MPI_Offset((nnz-myOffsetNnz)*sizeof(int) + myOffsetNnz*sizeof(double));
	if(MPI_File_set_view(fh, displacement, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL)!=0) {
		stringstream ss;
		ss<<"ERROR! "<<funcID<<": Cannot set the third MPI_File_set_view -- offset="<<displacement<<endl;
		string errorMessage = ss.str();
		cerr<<errorMessage;

		delete [] filenameArray;
		if(localNnzArray != NULL) { delete [] localNnzArray; 	localNnzArray = NULL; }
		throw(-1);
	}

	double *bufferDouble = new double [myNnz];
	count = 0;
	int diagCount = 0;
	for(int icv=0; icv<ncv; ++icv) {
	    int noc_f = nbocv_i[icv];
	    int noc_l = nbocv_i[icv+1]-1;

	    for(int noc = noc_f; noc <= noc_l; ++noc) {
			for(int i=0; i<nRowsPerBlock; ++i) {
				for(int j=0; j<nColsPerBlock; ++j) {
					bufferDouble[count] = A[noc][i][j];
					++count;

					if(noc == noc_f) {
						assert(icv == nbocv_v[noc]);
						int currGlobalRow = (cvora[mpi_rank] + icv)*nRowsPerBlock + i;
						int currGlobalCol = nbocv_v_golbal[icv]*nColsPerBlock + j;

						if(currGlobalRow == currGlobalCol) {
							double tmp = cv_volume[icv]/(local_dt[icv]);
							bufferDouble[count] -= tmp;

							++diagCount;
						}
					}
				}
			}
	    }
	}
	assert(count == myNnz);
	assert(diagCount == myNRows);
	MPI_File_write(fh, bufferDouble, myNnz, MPI_DOUBLE, &status);
	delete [] bufferDouble; 	bufferDouble = NULL;

	// 2-4. Close the file
	MPI_File_close(&fh);
	MPI_Barrier(mpi_comm);

	/*
	 * 3. Write the foot
	 */
	if(mpi_rank == mpi_size-1) {
		ofile.open(filenameArray, ios_base::out | ios_base::app | ios_base::binary);
		dummyInt = EOF_ERROR_CHECK_CODE; 	ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		ofile.close();
	}
	MPI_Barrier(mpi_comm);

	if(mpi_rank==0)
		wtimeF = MPI_Wtime();

	/***********
	 ** Show the summary on the screen if the debugging level is high
	 ***********/
	if(mpi_rank==0)
		printf("> The Jacobian matrix is written on %s (RUN TIME = %.3f [sec]): nRows=%d, nCols=%d, nnz=%d \n",
				filenameArray, wtimeF-wtime0, nRows, nCols, nnz);

	/***********
	 ** Free the memory
	 ***********/
	delete [] filenameArray;
	if(localNnzArray != NULL)
		delete [] localNnzArray;
}

/*
 * Method: writeJoeUncoupledMatrixBinaryParallel
 * ---------------------------------------------
 * Arguments:
 *   cvora          = number of ncv's at each mpi_rank
 *   nbocv_i        = indices information of nbocv_v. This can be also used to find the number of neighbors of a given icv
 *   nbocv_v_gloval = global neighboring icv information of a given icv. To access this array, nbocv_i is required.
 */
void writeJoeUncoupledNSMatrixBinaryParallel(const char filename[], double (*A)[5][5], const int *cvora,
		const int ncv, const int *nbocv_i, const int *nbocv_v_golbal) {
	string filenameString(filename);
	writeJoeUncoupledNSMatrixBinaryParallel(filenameString, A, cvora, ncv, nbocv_i, nbocv_v_golbal);
}

void writeJoeUncoupledNSMatrixBinaryParallel(string &filename, double (*A)[5][5], const int *cvora,
		const int ncv, const int *nbocv_i, const int *nbocv_v_golbal) {
	string funcID = "JoeUgpWithCvCompFlow.h::writeJoeUncoupledNSMatrixBinaryParallel()";

	if(A == NULL) {
		if(mpi_rank==0)
			cerr<<"ERROR "<<funcID<<": The given matrix A is NULL"<<endl;
		throw(-1);
	}
	if(cvora == NULL || ncv == 0) {
		if(mpi_rank==0)
			cerr<<"ERROR "<<funcID<<": The given CV information (cvora, ncv) is NULL"<<endl;
		throw(-1);
	}
	if(nbocv_i == NULL || nbocv_v_golbal == NULL) {
		if(mpi_rank==0)
			cerr<<"ERROR "<<funcID<<": The given neighbor information (nbocv_i, nbocv_v_golbal) is NULL"<<endl;
		throw(-1);
	}

	char *filenameArray = new char[filename.size() + 1]; // Since MPI_File_open() doesn't accept "const char *", you cannot use std::string.c_str()
	std::copy(filename.begin(), filename.end(), filenameArray);
	filenameArray[filename.size()] = '\0';

	/***********
	 ** Initialization: Get nRows, nCols, nnz, localNnzArray
	 ***********/
	int nRowsPerBlock = 5;
	int nColsPerBlock = 5;
	int nbocv_s = nbocv_i[ncv]; // Total number of non-zero blocks in a mpi_rank

	int myNRows = ncv     * nRowsPerBlock;
	int myNnz   = nbocv_s * nRowsPerBlock*nColsPerBlock;

	int nRows;
	MPI_Allreduce(&myNRows, &nRows, 1, MPI_INT, MPI_SUM, mpi_comm);
	int nnz;
	MPI_Allreduce(&myNnz,   &nnz,   1, MPI_INT, MPI_SUM, mpi_comm);

	int nCols = cvora[mpi_size]*nColsPerBlock;

	int *localNnzArray = NULL;
	if(mpi_rank == 0)
		localNnzArray = new int [mpi_size];
	MPI_Gather(&myNnz, 1, MPI_INT, localNnzArray, 1, MPI_INT, 0, mpi_comm);

	/***********
	 ** Write the data
	 ***********/
	double wtime0, wtimeF;
	if(mpi_rank==0)
		wtime0 = MPI_Wtime();

	/*
	 * 1. Write the header part
	 */
	ofstream ofile;
	int dummyInt;

	if(mpi_rank == 0) {
		ofile.open(filenameArray, ios_base::out | ios_base::trunc | ios_base::binary);

		// Write the header (number of global rows, number of global columns, global nnz, local nnz array)
		dummyInt = nRows; 					ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int)); // Note: if you simply use the "<<" operator, it will cause some problem due to bugs in a g++ library on linux
		dummyInt = nCols; 					ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		dummyInt = nnz; 					ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		dummyInt = mpi_size; 				ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		for(int i=0; i<mpi_size; ++i) {
			dummyInt = localNnzArray[i]; 	ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		}
		ofile.close();
	}
	MPI_Barrier(mpi_comm);

	int initDisp = sizeof(int)*(4 + mpi_size);  // Initial displacement due to the header

	/*
	 * 2. Write the body (matrix)
	 */
	MPI_Status status;
	MPI_File   fh;
	MPI_Offset displacement;

	// Open the file
	if (MPI_File_open(mpi_comm, filenameArray, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh)!=0) {
			// Note: MPI_MODE_WRONLY = write only
			//       MPI_MODE_CREATE = create the file if it does not exist.
		stringstream ss;
		ss<<"ERROR! "<<funcID<<": Cannot open "<<filenameArray;
		string errorMessage = ss.str();

		if(mpi_rank==0) cerr<<errorMessage;

		delete [] filenameArray;
		if(localNnzArray != NULL) { delete [] localNnzArray; 	localNnzArray = NULL; }
		throw(-1);
	}

	// 2-1. Write the row indices first
	int myOffsetNnz;
	MPI_Scan(&myNnz, &myOffsetNnz, 1, MPI_INT, MPI_SUM, mpi_comm);
	myOffsetNnz -= myNnz;
	displacement = MPI_Offset(initDisp + sizeof(int)*myOffsetNnz);
	if(MPI_File_set_view(fh, displacement, MPI_INT, MPI_INT, "native", MPI_INFO_NULL)!=0) {
		// Note: native     = Data in this representation are stored in a file exactly as it is in memory
		//       Internal   = Data in this representation can be used for I/O operations in a homogeneous or heterogeneous environment
		//       External32 = This data representation states that read and write operations convert all data from and to the "external32" representation
		stringstream ss;
		ss<<"ERROR! "<<funcID<<": Cannot open "<<filenameArray;
		string errorMessage = ss.str();

		if(mpi_rank==0) cerr<<errorMessage;

		delete [] filenameArray;
		if(localNnzArray != NULL) { delete [] localNnzArray; 	localNnzArray = NULL; }
		throw(-1);
	}

	int *bufferInt = new int [myNnz];
	int count = 0;
	for(int icv=0; icv<ncv; ++icv) {
	    int noc_f = nbocv_i[icv];
	    int noc_l = nbocv_i[icv+1]-1;

	    for(int noc = noc_f; noc <= noc_l; ++noc) {
			for(int i=0; i<nRowsPerBlock; ++i) {
				int currGlobalRow = (cvora[mpi_rank] + icv)*nRowsPerBlock + i;
				for(int j=0; j<nColsPerBlock; ++j) {
					bufferInt[count] = currGlobalRow;
					++count;
				}
			}
	    }
	}
	assert(count == myNnz);
	MPI_File_write(fh, bufferInt, myNnz, MPI_INT, &status);
	delete [] bufferInt; 	bufferInt = NULL;

	// 2-2. Write the column indices
	displacement += MPI_Offset(nnz*sizeof(int));
	if(MPI_File_set_view(fh, displacement, MPI_INT, MPI_INT, "native", MPI_INFO_NULL)!=0) {
		stringstream ss;
		ss<<"ERROR! "<<funcID<<": Cannot set the second MPI_File_set_view -- offset="<<displacement<<endl;
		string errorMessage = ss.str();

		cerr<<errorMessage;

		delete [] filenameArray;
		if(localNnzArray != NULL) { delete [] localNnzArray; 	localNnzArray = NULL; }
		throw(-1);
	}

	bufferInt = new int [myNnz];
	count = 0;
	for(int icv=0; icv<ncv; ++icv) {
	    int noc_f = nbocv_i[icv];
	    int noc_l = nbocv_i[icv+1]-1;

	    for(int noc = noc_f; noc <= noc_l; ++noc) {
	    	int icvNbr = nbocv_v[noc];
			for(int i=0; i<nRowsPerBlock; ++i) {
				for(int j=0; j<nColsPerBlock; ++j) {
					int currGlobalCol = nbocv_v_golbal[icvNbr]*nColsPerBlock + j;
					assert(currGlobalCol >= 0 && currGlobalCol <= nCols);
					bufferInt[count] = currGlobalCol;
					++count;
				}
			}
	    }
	}
	MPI_File_write(fh, bufferInt, myNnz, MPI_INT, &status);
	delete [] bufferInt; 	bufferInt = NULL;

	// 2-3. Write the values
	displacement += MPI_Offset((nnz-myOffsetNnz)*sizeof(int) + myOffsetNnz*sizeof(double));
	if(MPI_File_set_view(fh, displacement, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL)!=0) {
		stringstream ss;
		ss<<"ERROR! "<<funcID<<": Cannot set the third MPI_File_set_view -- offset="<<displacement<<endl;
		string errorMessage = ss.str();
		cerr<<errorMessage;

		delete [] filenameArray;
		if(localNnzArray != NULL) { delete [] localNnzArray; 	localNnzArray = NULL; }
		throw(-1);
	}

	double *bufferDouble = new double [myNnz];
	count = 0;
	for(int icv=0; icv<ncv; ++icv) {
	    int noc_f = nbocv_i[icv];
	    int noc_l = nbocv_i[icv+1]-1;

	    for(int noc = noc_f; noc <= noc_l; ++noc) {
			for(int i=0; i<nRowsPerBlock; ++i) {
				for(int j=0; j<nColsPerBlock; ++j) {
					bufferDouble[count] = A[noc][i][j];
					++count;
				}
			}
	    }
	}
	assert(count == myNnz);
	MPI_File_write(fh, bufferDouble, myNnz, MPI_DOUBLE, &status);
	delete [] bufferDouble; 	bufferDouble = NULL;

	// 2-4. Close the file
	MPI_File_close(&fh);
	MPI_Barrier(mpi_comm);

	/*
	 * 3. Write the foot
	 */
	if(mpi_rank == mpi_size-1) {
		ofile.open(filenameArray, ios_base::out | ios_base::app | ios_base::binary);
		dummyInt = EOF_ERROR_CHECK_CODE; 	ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		ofile.close();
	}
	MPI_Barrier(mpi_comm);

	if(mpi_rank==0)
		wtimeF = MPI_Wtime();

	/***********
	 ** Show the summary on the screen if the debugging level is high
	 ***********/
	if(mpi_rank==0)
		printf("> The Jacobian matrix is written on %s (RUN TIME = %.3f [sec]): nRows=%d, nCols=%d, nnz=%d \n",
				filenameArray, wtimeF-wtime0, nRows, nCols, nnz);

	/***********
	 ** Free the memory
	 ***********/
	delete [] filenameArray;
	if(localNnzArray != NULL)
		delete [] localNnzArray;
}

/*
 * Method: writeJoeUncoupledNSMatrixWeightedDiagBinaryParallel
 * -----------------------------------------------------------
 * Similar to writeJoeUncoupledNSMatrixBinaryParallel().
 * The difference is each row is multiplied by 1 / (vol/dt) for some eigen-analysis and stability analysis.
 * Arguments:
 *   cvora          = number of ncv's at each mpi_rank
 *   nbocv_i        = indices information of nbocv_v. This can be also used to find the number of neighbors of a given icv
 *   nbocv_v_gloval = global neighboring icv information of a given icv. To access this array, nbocv_i is required.
 */
void writeJoeUncoupledNSMatrixWeightedDiagBinaryParallel(const char filename[], double (*A)[5][5], const int *cvora,
		const int ncv, const int *nbocv_i, const int *nbocv_v_golbal,
		const double *cv_volume, const double *local_dt) {
	string filenameString(filename);
	writeJoeUncoupledNSMatrixWeightedDiagBinaryParallel(filenameString, A, cvora, ncv, nbocv_i, nbocv_v_golbal, cv_volume, local_dt);
}

void writeJoeUncoupledNSMatrixWeightedDiagBinaryParallel(string &filename, double (*A)[5][5], const int *cvora,
		const int ncv, const int *nbocv_i, const int *nbocv_v_golbal,
		const double *cv_volume, const double *local_dt) {
	string funcID = "JoeUgpWithCvCompFlow.h::writeJoeUncoupledNSMatrixWeightedDiagBinaryParallel()";

	if(A == NULL) {
		if(mpi_rank==0)
			cerr<<"ERROR "<<funcID<<": The given matrix A is NULL"<<endl;
		throw(-1);
	}
	if(cvora == NULL || ncv == 0) {
		if(mpi_rank==0)
			cerr<<"ERROR "<<funcID<<": The given CV information (cvora, ncv) is NULL"<<endl;
		throw(-1);
	}
	if(nbocv_i == NULL || nbocv_v_golbal == NULL) {
		if(mpi_rank==0)
			cerr<<"ERROR "<<funcID<<": The given neighbor information (nbocv_i, nbocv_v_golbal) is NULL"<<endl;
		throw(-1);
	}

	char *filenameArray = new char[filename.size() + 1]; // Since MPI_File_open() doesn't accept "const char *", you cannot use std::string.c_str()
	std::copy(filename.begin(), filename.end(), filenameArray);
	filenameArray[filename.size()] = '\0';

	/***********
	 ** Initialization: Get nRows, nCols, nnz, localNnzArray
	 ***********/
	int nRowsPerBlock = 5;
	int nColsPerBlock = 5;
	int nbocv_s = nbocv_i[ncv]; // Total number of non-zero blocks in a mpi_rank

	int myNRows = ncv     * nRowsPerBlock;
	int myNnz   = nbocv_s * nRowsPerBlock*nColsPerBlock;

	int nRows;
	MPI_Allreduce(&myNRows, &nRows, 1, MPI_INT, MPI_SUM, mpi_comm);
	int nnz;
	MPI_Allreduce(&myNnz,   &nnz,   1, MPI_INT, MPI_SUM, mpi_comm);

	int nCols = cvora[mpi_size]*nColsPerBlock;

	int *localNnzArray = NULL;
	if(mpi_rank == 0)
		localNnzArray = new int [mpi_size];
	MPI_Gather(&myNnz, 1, MPI_INT, localNnzArray, 1, MPI_INT, 0, mpi_comm);

	/***********
	 ** Write the data
	 ***********/
	double wtime0, wtimeF;
	if(mpi_rank==0)
		wtime0 = MPI_Wtime();

	/*
	 * 1. Write the header part
	 */
	ofstream ofile;
	int dummyInt;

	if(mpi_rank == 0) {
		ofile.open(filenameArray, ios_base::out | ios_base::trunc | ios_base::binary);

		// Write the header (number of global rows, number of global columns, global nnz, local nnz array)
		dummyInt = nRows; 					ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int)); // Note: if you simply use the "<<" operator, it will cause some problem due to bugs in a g++ library on linux
		dummyInt = nCols; 					ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		dummyInt = nnz; 					ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		dummyInt = mpi_size; 				ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		for(int i=0; i<mpi_size; ++i) {
			dummyInt = localNnzArray[i]; 	ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		}
		ofile.close();
	}
	MPI_Barrier(mpi_comm);

	int initDisp = sizeof(int)*(4 + mpi_size);  // Initial displacement due to the header

	/*
	 * 2. Write the body (matrix)
	 */
	MPI_Status status;
	MPI_File   fh;
	MPI_Offset displacement;

	// Open the file
	if (MPI_File_open(mpi_comm, filenameArray, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh)!=0) {
			// Note: MPI_MODE_WRONLY = write only
			//       MPI_MODE_CREATE = create the file if it does not exist.
		stringstream ss;
		ss<<"ERROR! "<<funcID<<": Cannot open "<<filenameArray;
		string errorMessage = ss.str();

		if(mpi_rank==0) cerr<<errorMessage;

		delete [] filenameArray;
		if(localNnzArray != NULL) { delete [] localNnzArray; 	localNnzArray = NULL; }
		throw(-1);
	}

	// 2-1. Write the row indices first
	int myOffsetNnz;
	MPI_Scan(&myNnz, &myOffsetNnz, 1, MPI_INT, MPI_SUM, mpi_comm);
	myOffsetNnz -= myNnz;
	displacement = MPI_Offset(initDisp + sizeof(int)*myOffsetNnz);
	if(MPI_File_set_view(fh, displacement, MPI_INT, MPI_INT, "native", MPI_INFO_NULL)!=0) {
		// Note: native     = Data in this representation are stored in a file exactly as it is in memory
		//       Internal   = Data in this representation can be used for I/O operations in a homogeneous or heterogeneous environment
		//       External32 = This data representation states that read and write operations convert all data from and to the "external32" representation
		stringstream ss;
		ss<<"ERROR! "<<funcID<<": Cannot open "<<filenameArray;
		string errorMessage = ss.str();

		if(mpi_rank==0) cerr<<errorMessage;

		delete [] filenameArray;
		if(localNnzArray != NULL) { delete [] localNnzArray; 	localNnzArray = NULL; }
		throw(-1);
	}

	int *bufferInt = new int [myNnz];
	int count = 0;
	for(int icv=0; icv<ncv; ++icv) {
	    int noc_f = nbocv_i[icv];
	    int noc_l = nbocv_i[icv+1]-1;

	    for(int noc = noc_f; noc <= noc_l; ++noc) {
			for(int i=0; i<nRowsPerBlock; ++i) {
				int currGlobalRow = (cvora[mpi_rank] + icv)*nRowsPerBlock + i;
				for(int j=0; j<nColsPerBlock; ++j) {
					bufferInt[count] = currGlobalRow;
					++count;
				}
			}
	    }
	}
	assert(count == myNnz);
	MPI_File_write(fh, bufferInt, myNnz, MPI_INT, &status);
	delete [] bufferInt; 	bufferInt = NULL;

	// 2-2. Write the column indices
	displacement += MPI_Offset(nnz*sizeof(int));
	if(MPI_File_set_view(fh, displacement, MPI_INT, MPI_INT, "native", MPI_INFO_NULL)!=0) {
		stringstream ss;
		ss<<"ERROR! "<<funcID<<": Cannot set the second MPI_File_set_view -- offset="<<displacement<<endl;
		string errorMessage = ss.str();

		cerr<<errorMessage;

		delete [] filenameArray;
		if(localNnzArray != NULL) { delete [] localNnzArray; 	localNnzArray = NULL; }
		throw(-1);
	}

	bufferInt = new int [myNnz];
	count = 0;
	for(int icv=0; icv<ncv; ++icv) {
	    int noc_f = nbocv_i[icv];
	    int noc_l = nbocv_i[icv+1]-1;

	    for(int noc = noc_f; noc <= noc_l; ++noc) {
	    	int icvNbr = nbocv_v[noc];
			for(int i=0; i<nRowsPerBlock; ++i) {
				for(int j=0; j<nColsPerBlock; ++j) {
					int currGlobalCol = nbocv_v_golbal[icvNbr]*nColsPerBlock + j;
					assert(currGlobalCol >= 0 && currGlobalCol <= nCols);
					bufferInt[count] = currGlobalCol;
					++count;
				}
			}
	    }
	}
	MPI_File_write(fh, bufferInt, myNnz, MPI_INT, &status);
	delete [] bufferInt; 	bufferInt = NULL;

	// 2-3. Write the values
	displacement += MPI_Offset((nnz-myOffsetNnz)*sizeof(int) + myOffsetNnz*sizeof(double));
	if(MPI_File_set_view(fh, displacement, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL)!=0) {
		stringstream ss;
		ss<<"ERROR! "<<funcID<<": Cannot set the third MPI_File_set_view -- offset="<<displacement<<endl;
		string errorMessage = ss.str();
		cerr<<errorMessage;

		delete [] filenameArray;
		if(localNnzArray != NULL) { delete [] localNnzArray; 	localNnzArray = NULL; }
		throw(-1);
	}

	double *bufferDouble = new double [myNnz];
	count = 0;
	for(int icv=0; icv<ncv; ++icv) {
	    int noc_f = nbocv_i[icv];
	    int noc_l = nbocv_i[icv+1]-1;

	    double tmp = cv_volume[icv] / local_dt[icv];

	    for(int noc = noc_f; noc <= noc_l; ++noc) {
			for(int i=0; i<nRowsPerBlock; ++i) {
				for(int j=0; j<nColsPerBlock; ++j) {
					bufferDouble[count] = A[noc][i][j] / tmp;
					++count;
				}
			}
	    }
	}
	assert(count == myNnz);
	MPI_File_write(fh, bufferDouble, myNnz, MPI_DOUBLE, &status);
	delete [] bufferDouble; 	bufferDouble = NULL;

	// 2-4. Close the file
	MPI_File_close(&fh);
	MPI_Barrier(mpi_comm);

	/*
	 * 3. Write the foot
	 */
	if(mpi_rank == mpi_size-1) {
		ofile.open(filenameArray, ios_base::out | ios_base::app | ios_base::binary);
		dummyInt = EOF_ERROR_CHECK_CODE; 	ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		ofile.close();
	}
	MPI_Barrier(mpi_comm);

	if(mpi_rank==0)
		wtimeF = MPI_Wtime();

	/***********
	 ** Show the summary on the screen if the debugging level is high
	 ***********/
	if(mpi_rank==0)
		printf("> The Jacobian matrix is written on %s (RUN TIME = %.3f [sec]): nRows=%d, nCols=%d, nnz=%d \n",
				filenameArray, wtimeF-wtime0, nRows, nCols, nnz);

	/***********
	 ** Free the memory
	 ***********/
	delete [] filenameArray;
	if(localNnzArray != NULL)
		delete [] localNnzArray;
}

/*
 * Method: writeJoeUncoupledNSrhsBinaryParallel
 * --------------------------------------------
 * Arguments: cvora = number of ncv's at each mpi_rank
 */
void writeJoeUncoupledNSrhsBinaryParallel(const char filename[], double *RHSrho, double (*RHSrhou)[3], double *RHSrhoE,
		const int *cvora, const int ncv) {
	string filenameString(filename);
	writeJoeUncoupledNSrhsBinaryParallel(filenameString, RHSrho, RHSrhou, RHSrhoE, cvora, ncv);
}

void writeJoeUncoupledNSrhsBinaryParallel(string &filename, double *RHSrho, double (*RHSrhou)[3], double *RHSrhoE,
		const int *cvora, const int ncv) {
	string funcID = "JoeUgpWithCvCompFlow.h::writeJoeUncoupledNSrhsBinaryParallel()";

	if(RHSrho == NULL || RHSrhou == NULL || RHSrhoE == NULL) {
		if(mpi_rank==0)
			cerr<<"ERROR "<<funcID<<": The given rhs vectors (RHSrho, RHSrhou, RHSrhoE) is NULL"<<endl;
		throw(-1);
	}
	if(cvora == NULL || ncv == 0) {
		if(mpi_rank==0)
			cerr<<"ERROR "<<funcID<<": The given CV information (cvora, ncv) is NULL"<<endl;
		throw(-1);
	}

	char *filenameArray = new char[filename.size() + 1]; // Since MPI_File_open() doesn't accept "const char *", you cannot use std::string.c_str()
	std::copy(filename.begin(), filename.end(), filenameArray);
	filenameArray[filename.size()] = '\0';

	/***********
	 ** Initialization: Get NrowsTotal, localNrowsArray
	 ***********/
	int nRowsPerBlock = 5;
	int myNrows    = ncv * nRowsPerBlock;
	int NrowsTotal = cvora[mpi_size] * nRowsPerBlock;

	int *localNrowsArray = NULL;
	if(mpi_rank == 0)
		localNrowsArray = new int [mpi_size];
	MPI_Gather(&myNrows, 1, MPI_INT, localNrowsArray, 1, MPI_INT, 0, mpi_comm);

	/***********
	 ** Write the data
	 ***********/
	double wtime0, wtimeF;
	if(mpi_rank==0)
		wtime0 = MPI_Wtime();

	/*
	 * 1. Write the header part
	 */
	ofstream ofile;
	int dummyInt;

	if(mpi_rank == 0) {
		ofile.open(filenameArray, ios_base::out | ios_base::trunc | ios_base::binary);

		// Write the header (number of global rows)
		dummyInt = NrowsTotal; 				ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int)); // Note: if you simply use the "<<" operator, it will cause some problem due to bugs in a g++ library on linux
		dummyInt = nRowsPerBlock; 			ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		dummyInt = mpi_size; 				ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		for(int i=0; i<mpi_size; ++i) {
			dummyInt = localNrowsArray[i]; 	ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		}
		ofile.close();
	}
	MPI_Barrier(mpi_comm);

	int initDisp = sizeof(int)*(3 + mpi_size);  // Initial displacement due to the header

	/*
	 * 2. Write the body (vector)
	 */
	MPI_Status status;
	MPI_File   fh;
	MPI_Offset displacement;

	// Open the file
	if (MPI_File_open(mpi_comm, filenameArray, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh) != 0) {
			// Note: MPI_MODE_WRONLY = write only
			//       MPI_MODE_CREATE = create the file if it does not exist.
		if(mpi_rank==0) cerr<<"ERROR! "<<funcID<<": Cannot open "<<filenameArray<<endl;

		delete [] filenameArray;
		if(localNrowsArray != NULL) { delete [] localNrowsArray; 	localNrowsArray = NULL; }
		throw(-1);
	}

	// 2-1. Write the row indices first
	int myOffsetNnz;
	MPI_Scan(&myNrows, &myOffsetNnz, 1, MPI_INT, MPI_SUM, mpi_comm);
	myOffsetNnz -= myNrows;
	displacement = MPI_Offset(initDisp + sizeof(double)*myOffsetNnz);

	if(MPI_File_set_view(fh, displacement, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL)!=0) {
		stringstream ss;
		ss<<"ERROR! "<<funcID<<": Cannot set the third MPI_File_set_view -- offset="<<displacement<<endl;
		string errorMessage = ss.str();
		cerr<<errorMessage;

		delete [] filenameArray;
		if(localNrowsArray != NULL) { delete [] localNrowsArray; 	localNrowsArray = NULL; }
		throw(-1);
	}

	double *bufferDouble = new double [myNrows];
	int count = 0;
	for(int icv=0; icv<ncv; ++icv) {
		bufferDouble[count] = RHSrho[icv];
		++count;

		for(int i=0; i<3; ++i) {
			bufferDouble[count] = RHSrhou[icv][i];
			++count;
		}

		bufferDouble[count] = RHSrhoE[icv];
		++count;
	}
	assert(count == myNrows);
	MPI_File_write(fh, bufferDouble, myNrows, MPI_DOUBLE, &status);
	delete [] bufferDouble; 	bufferDouble = NULL;

	// 2-4. Close the file
	MPI_File_close(&fh);
	MPI_Barrier(mpi_comm);

	/*
	 * 3. Write the foot
	 */
	if(mpi_rank == mpi_size-1) {
		ofile.open(filenameArray, ios_base::out | ios_base::app | ios_base::binary);
		dummyInt = EOF_ERROR_CHECK_CODE; 	ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		ofile.close();
	}
	MPI_Barrier(mpi_comm);

	if(mpi_rank==0)
		wtimeF = MPI_Wtime();

	/***********
	 ** Show the summary on the screen if the debugging level is high
	 ***********/
	if(mpi_rank==0)
		printf("> The given RHS vector is written on %s (RUN TIME = %.3f [sec]): nRows=%d\n",
				filenameArray, wtimeF-wtime0, NrowsTotal);

	/***********
	 ** Free the memory
	 ***********/
	delete [] filenameArray;
	if(localNrowsArray != NULL)
		delete [] localNrowsArray;
}

/*
 * Method: writeJoe2DvectorBinaryParallel
 * --------------------------------------
 * Arguments: rhs   = 2D vector: You can get it by calling  getMem2D(&rhs, 0, ncv-1, 0, nRowsPerBlock-1, "Name", true);
 *                               and free it by calling     freeMem2D(rhs, 0, ncv-1, 0, nRowsPerBlock-1); 	rhs = NULL;
 *            cvora = number of ncv's at each mpi_rank
 */
void writeJoe2DvectorBinaryParallel(const char filename[], double **rhs, const int nRowsPerBlock,
		const int *cvora, const int ncv) {
	string filenameString(filename);
	writeJoe2DvectorBinaryParallel(filenameString, rhs, nRowsPerBlock, cvora, ncv);
}

void writeJoe2DvectorBinaryParallel(string &filename, double **rhs, const int nRowsPerBlock,
		const int *cvora, const int ncv) {
	string funcID = "JoeUgpWithCvCompFlow.h::writeJoe2DrhsBinaryParallel()";

	if(rhs == NULL) {
		if(mpi_rank==0)
			cerr<<"ERROR "<<funcID<<": The given rhs vectors is NULL"<<endl;
		throw(-1);
	}
	if(nRowsPerBlock == 0) {
		if(mpi_rank==0)
			cerr<<"ERROR "<<funcID<<": The given nRowsPerBlock is ZERO"<<endl;
		throw(-1);
	}
	if(cvora == NULL || ncv == 0) {
		if(mpi_rank==0)
			cerr<<"ERROR "<<funcID<<": The given CV information (cvora, ncv) is NULL"<<endl;
		throw(-1);
	}

	char *filenameArray = new char[filename.size() + 1]; // Since MPI_File_open() doesn't accept "const char *", you cannot use std::string.c_str()
	std::copy(filename.begin(), filename.end(), filenameArray);
	filenameArray[filename.size()] = '\0';

	/***********
	 ** Initialization: Get NrowsTotal, localNrowsArray
	 ***********/
	int myNrows    = ncv * nRowsPerBlock;
	int NrowsTotal = cvora[mpi_size] * nRowsPerBlock;

	int *localNrowsArray = NULL;
	if(mpi_rank == 0)
		localNrowsArray = new int [mpi_size];
	MPI_Gather(&myNrows, 1, MPI_INT, localNrowsArray, 1, MPI_INT, 0, mpi_comm);

	/***********
	 ** Write the data
	 ***********/
	double wtime0, wtimeF;
	if(mpi_rank==0)
		wtime0 = MPI_Wtime();

	double *myMinValues = new double [nRowsPerBlock];
	double *myMaxValues = new double [nRowsPerBlock];
	for(int i=0; i<nRowsPerBlock; ++i) {
		myMinValues[i] =  ABSURDLY_BIG_NUMBER;
		myMaxValues[i] = -ABSURDLY_BIG_NUMBER;
	}

	/*
	 * 1. Write the header part
	 */
	ofstream ofile;
	int dummyInt;

	if(mpi_rank == 0) {
		ofile.open(filenameArray, ios_base::out | ios_base::trunc | ios_base::binary);

		// Write the header (number of global rows)
		dummyInt = NrowsTotal; 				ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int)); // Note: if you simply use the "<<" operator, it will cause some problem due to bugs in a g++ library on linux
		dummyInt = nRowsPerBlock; 			ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		dummyInt = mpi_size; 				ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		for(int i=0; i<mpi_size; ++i) {
			dummyInt = localNrowsArray[i]; 	ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		}
		ofile.close();
	}
	MPI_Barrier(mpi_comm);

	int initDisp = sizeof(int)*(3 + mpi_size);  // Initial displacement due to the header

	/*
	 * 2. Write the body (vector)
	 */
	MPI_Status status;
	MPI_File   fh;
	MPI_Offset displacement;

	// Open the file
	if (MPI_File_open(mpi_comm, filenameArray, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh) != 0) {
			// Note: MPI_MODE_WRONLY = write only
			//       MPI_MODE_CREATE = create the file if it does not exist.
		if(mpi_rank==0) cerr<<"ERROR! "<<funcID<<": Cannot open "<<filenameArray<<endl;

		delete [] filenameArray;
		if(localNrowsArray != NULL) { delete [] localNrowsArray; 	localNrowsArray = NULL; }
		delete [] myMinValues; 	delete [] myMaxValues;
		throw(-1);
	}

	// 2-1. Write the row indices first
	int myOffsetNnz;
	MPI_Scan(&myNrows, &myOffsetNnz, 1, MPI_INT, MPI_SUM, mpi_comm);
	myOffsetNnz -= myNrows;
	displacement = MPI_Offset(initDisp + sizeof(double)*myOffsetNnz);

	if(MPI_File_set_view(fh, displacement, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL)!=0) {
		stringstream ss;
		ss<<"ERROR! "<<funcID<<": Cannot set the third MPI_File_set_view -- offset="<<displacement<<endl;
		string errorMessage = ss.str();
		cerr<<errorMessage;

		delete [] filenameArray;
		if(localNrowsArray != NULL) { delete [] localNrowsArray; 	localNrowsArray = NULL; }
		delete [] myMinValues; 	delete [] myMaxValues;
		throw(-1);
	}

	double *bufferDouble = new double [myNrows];
	int count = 0;
	for(int icv=0; icv<ncv; ++icv) {
		for(int i=0; i<nRowsPerBlock; ++i) {
			bufferDouble[count] = rhs[icv][i];
			++count;

			myMinValues[i] = min(myMinValues[i], rhs[icv][i]);
			myMaxValues[i] = max(myMaxValues[i], rhs[icv][i]);
		}
	}
	assert(count == myNrows);
	MPI_File_write(fh, bufferDouble, myNrows, MPI_DOUBLE, &status);
	delete [] bufferDouble; 	bufferDouble = NULL;

	// 2-4. Close the file
	MPI_File_close(&fh);
	MPI_Barrier(mpi_comm);

	/*
	 * 3. Write the foot
	 */
	if(mpi_rank == mpi_size-1) {
		ofile.open(filenameArray, ios_base::out | ios_base::app | ios_base::binary);
		dummyInt = EOF_ERROR_CHECK_CODE; 	ofile.write(reinterpret_cast<char*>(&dummyInt), sizeof(int));
		ofile.close();
	}
	MPI_Barrier(mpi_comm);

	if(mpi_rank==0)
		wtimeF = MPI_Wtime();

	/***********
	 ** Show the summary on the screen if the debugging level is high
	 ***********/
	double *minValues = new double [nRowsPerBlock];
	double *maxValues = new double [nRowsPerBlock];
	MPI_Allreduce(myMinValues, minValues, nRowsPerBlock, MPI_DOUBLE, MPI_MIN, mpi_comm);
	MPI_Allreduce(myMaxValues, maxValues, nRowsPerBlock, MPI_DOUBLE, MPI_MAX, mpi_comm);

	if(mpi_rank==0) {
		printf("> The given vector is written on %s (RUN TIME = %.3f [sec]): nRows=%d, nRowsPerBlock=%d\n",
				filenameArray, wtimeF-wtime0, NrowsTotal, nRowsPerBlock);
		printf("  Data range:");
		for(int i=0; i<nRowsPerBlock; ++i)
			cout<<" ["<<minValues[i]<<","<<maxValues[i]<<"]";
		cout<<endl;
	}

	/***********
	 ** Free the memory
	 ***********/
	delete [] filenameArray;
	if(localNrowsArray != NULL)
		delete [] localNrowsArray;
	delete [] myMinValues; 	delete [] myMaxValues;
	delete [] minValues; 	delete [] maxValues;
}

#endif /* JOEUGPWITHCVCOMPFLOW_H_ */
