/*
 * MatComprsed.h
 *
 *  Created on: Oct 26, 2011
 *      Author: ikj
 */

#ifndef MATCOMPRSED_H_
#define MATCOMPRSED_H_

#include <assert.h>
#include <vector>
#include <utility>      // std::pair
//#include <iostream>

#define SLOW_BUT_ERROR_CHECK
#define MATCOMPRESED_NOT_FOUND -2
//#define USE_MPI_
#define USE_GERSHGORIN  // Use the Gershroin circle theorem to estimate eigenvalues

static bool ShowMatcompresedDetails = true;

/*
 * class: MatComprsed
 * ------------------
 * Matrix in a compressed form (row form)
 */
class MatComprsed
{
public:
	MatComprsed() {
		nnz = 0;
		rind = NULL;
		cind = NULL;
		values = NULL;

		ncols_eachRow_i = NULL;
		diag_index = NULL;

		nRows = 0;
		nCols = 0;
	}
	~MatComprsed() {
		clear();
	}

	bool empty() {
		if (nnz==0 && rind==NULL && cind==NULL && values==NULL && ncols_eachRow_i==NULL && diag_index==NULL && nRows==0 && nCols==0)
			return true;
		else
			return false;
	}

	int get_nnz() {
		return nnz;
	}
	int get_nCols() {
		return nCols;
	}
	int get_nRows() {
		return nRows;
	}
	unsigned int get_rind(const int index) {
		return rind[index];
	}
	unsigned int get_cind(const int index) {
		return cind[index];
	}
	double get_values(const int index) {
		return values[index];
	}
	int get_ncols_eachRow_i(const int index) {
		return ncols_eachRow_i[index];
	}
	int get_diag_index(const int index) {
		assert(index>=0 && index<get_nRows());
		return diag_index[index];
	}

	int get_global_rind(const int index, const int mpi_rank, const int nVars, const int *cvora) {
		int rindLocal = get_rind(index);
		return cvora[mpi_rank]*nVars + rindLocal;
	}
	int get_global_cind(const int index, const int mpi_rank, const int nVars, const int *cvora, const int *cv_gl, const int NcontrolEqns=0, const int ncv_gg=0) {
		int cindLocal = get_cind(index);
		int globalCind = MATCOMPRESED_NOT_FOUND; // This variable will be updated and returned

		// first, get icvGlobal(global icv)
		int icvLocal = (int) (cindLocal/nVars); // Note: icvLocal can be greater than ncv_gg if NcontrolEqns>0

		if(icvLocal<0) { 	// THIS IS AN ERROR!
			printf("ERROR in MatComprsed::get_global_cind(): mpi_rank=%d -- icvLocal=%d < 0 \n", mpi_rank, icvLocal);
			assert(icvLocal >= 0);
		} else if(icvLocal<ncv_gg) {
			int rest = cindLocal%nVars;
			int icvGlobal = cv_gl[icvLocal];
#ifdef SLOW_BUT_ERROR_CHECK
			if(icvGlobal<0 || icvGlobal>cvora[mpi_size]) {
				printf("ERROR in MatComprsed::get_global_cind(): mpi_rank=%d -- icvLocal=%d, icvGlobal=%d, cvora[mpi_size]=%d \n", mpi_rank, icvLocal, icvGlobal, cvora[mpi_size]);
				assert(icvGlobal>=0);
				assert(icvGlobal<=cvora[mpi_size]);
			}
#endif
			// then, calculate globalCind using icvGlobal and rest
			globalCind = icvGlobal*nVars + rest;
		} else if (icvLocal<ncv_gg+NcontrolEqns) {
			int iEqn = cindLocal - nVars*ncv_gg;
#ifdef SLOW_BUT_ERROR_CHECK
			assert(NcontrolEqns>0);
			if(iEqn<0 || iEqn>=NcontrolEqns) {
				printf("ERROR in MatComprsed::get_global_cind(): mpi_rank=%d -- iEqn=%d (!), NcontrolEqns=%d, cindLocal=%d, nVars*ncv_gg=%d \n", mpi_rank, iEqn, NcontrolEqns, cindLocal, nVars*ncv_gg);
				assert(iEqn>=0);
				assert(iEqn<NcontrolEqns);
			}
#endif
			// then, calculate globalCind using iEqn, get_nCols(), & NcontrolEqns
			globalCind = cvora[mpi_size]*nVars + iEqn;
		} else { 			// THIS IS AN ERROR!
			printf("ERROR in MatComprsed::get_global_cind(): mpi_rank=%d -- icvLocal=%d >= ncv_gg(=%d)+NcontrolEqns(=%d) \n", mpi_rank, icvLocal, ncv_gg, NcontrolEqns);
			assert(icvLocal<ncv_gg+NcontrolEqns);
		}

		return globalCind;
	}

	int get_local_rind(const int globalIndex, const int mpi_rank, const int nVars, const int *cvora) {
		return globalIndex - cvora[mpi_rank]*nVars;
	}
	int get_local_cind(const int globalIndex, const int mpi_rank, const int nVars, const int *cvora, const int *cv_gl, const int NcontrolEqns=0, const int ncv=0, const int ncv_gg=0) {
		int localIndex = MATCOMPRESED_NOT_FOUND; // This variable will be updated and returned

		int icvGlobal = (int) globalIndex/nVars;
		int rest = globalIndex%nVars;

		// First, get icvLocal(local icv)
		int icvLocal = icvGlobal - cvora[mpi_rank]; // For interior CVs, 0<=icvLocal<ncv.

		// There are three possibility: 1. interior CV   2. control Eqn   3. ghost CV
		if(icvLocal>=0 && icvLocal<ncv) {
			localIndex = icvLocal*nVars + rest;
		} else if(globalIndex >= cvora[mpi_size]*nVars) {
			assert(NcontrolEqns > 0);

			int iEqn = globalIndex - cvora[mpi_size]*nVars;
			assert(iEqn < NcontrolEqns);

			localIndex = ncv_gg*nVars + iEqn;
		} else {
			int foundICVindex = MATCOMPRESED_NOT_FOUND;
			for(int icv=ncv; icv<ncv_gg; ++icv) {
				if(cv_gl[icv] == icvGlobal) {
					foundICVindex = icv;
					break;
				}
			}

			if(foundICVindex == MATCOMPRESED_NOT_FOUND) {
				cerr<<"ERROR in MatComprsedSTL::get_local_cind(): cannot find local index"<<endl;
				assert(false);
			} else {
				localIndex = foundICVindex*nVars + rest;
			}
		}

		return localIndex;
	}

	void allocateMem(const int nnz) {
		assert(this->nnz == 0);
		this->nnz = nnz;
		assert(rind==NULL && cind==NULL && values==NULL);
		rind = new unsigned int [nnz];
		cind = new unsigned int [nnz];
		values = new double [nnz];

		assert(nRows>0 && nCols>0);
		assert(ncols_eachRow_i==NULL);
		ncols_eachRow_i = new int [nRows+1];
		assert(diag_index==NULL);
		diag_index = new int [nRows];
	}

	void set_rind(const int index, const unsigned int ind) {
		assert(index < nnz && ind < get_nRows());
		rind[index] = ind;
	}
	void set_cind(const int index, const unsigned int ind) {
		assert(index < nnz && ind < get_nCols());
		cind[index] = ind;
	}
	void convert_cindices(const vector<int> &myNbocv2, const int nVars) { // This method only for the MEM_SAVING_ADVAR case
		for(int i=0; i<nRows; ++i)
			diag_index[i] = MATCOMPRESED_NOT_FOUND;

		for(int i=0; i<get_nnz(); ++i) {
			int fake_icv = int(get_cind(i)/nVars);
			int rest = get_cind(i)%nVars;

			unsigned int tempCind = myNbocv2[fake_icv]*nVars+rest;
			set_cind(i, tempCind);

			if(get_rind(i) == tempCind) {
				set_diag_index(tempCind, i);
			}
		}
	}

	void set_values(const int index, const double val) {
		assert(index < nnz);
		values[index] = val;
	}
	void set_ncols_eachRow_i(const int index, const int ind) {
		assert(index < nnz);
		ncols_eachRow_i[index] = ind;
	}
	void set_diag_index(const int index, const int ind) {
		assert(index < nnz);
		diag_index[index] = ind;
	}

	void setMatSize(const int nCols, const int nRows) {
		assert( empty() );

		this->nCols = nCols;
		this->nRows = nRows;
	}

	void matCopy(const int nnz, const unsigned int *rind, const unsigned int *cind, const double *values) {
		// check if setMatSize() was already called
		assert(nRows>0 && nCols>0);

		// allocate memory
		allocateMem(nnz);

		// copy the matrix and set ncols_eachRow_i and diag_index
		/*   ncols_eachRow_i[row] indicates the accumulated number of the first non-zero element of the row in rind, cind, and values
		 *   e.g. A = [0 0 0 0
		 *             0 1 1 0
		 *             1 1 0 1
		 *             0 0 0 0]
		 *        ncols_eachRow_i[5] = {0, 0, 2, 5, 5}
		 *   Thus, the number of non-zeros in a row: ncols_eachRow_i[row+1] - ncols_eachRow_i[row]
		 */
		/*
		 *   diag_index[row] indicates the index of the diagonal element if the row in rind, cind, and values
		 *   e.g. A = [9 9 0 0
		 *             0 9 9 0
		 *             0 0 0 0
		 *             0 0 0 9]
		 *        rind[5] = {0, 0, 1, 1, 3}
		 *        cind[5] = {0, 1, 1, 2, 3}
		 *        values[5] = {9, 9, 9, 9, 9}
		 *
		 *        In the case that MATCOMPRESED_NOT_FOUND is defined as -1, diag_index[4] = {0, 2, -1, 4}
		 */
		for(int i=0; i<nRows+1; ++i)
			ncols_eachRow_i[i] = 0;
		for(int i=0; i<nRows; ++i)
			diag_index[i] = MATCOMPRESED_NOT_FOUND;

		int currRow = 0;
		int count = 0;
		for(int i=0; i<nnz; ++i) {
			if(rind[i]>=nRows)
				break;	// If the number of rows of the Jacobian matrix is greater than nRows, just stop feeding the matrix

			this->rind[i] = rind[i];
			this->cind[i] = cind[i];
			this->values[i] = values[i];

			if(currRow < rind[i]) {
				ncols_eachRow_i[currRow+1] = ncols_eachRow_i[currRow] + count;
				while(++currRow < rind[i]) {
					ncols_eachRow_i[currRow+1] = ncols_eachRow_i[currRow];
				}
				count = 0;
			}
			++count;

			if(currRow == cind[i]) {
				diag_index[currRow] = i;
			}
		}
		if(currRow == nRows-1) {
			ncols_eachRow_i[nRows] = ncols_eachRow_i[nRows-1] + count;
		} else { // check here!!!!!!
			ncols_eachRow_i[currRow+1] = ncols_eachRow_i[currRow] + count;
			while(currRow < nRows-1) {
				ncols_eachRow_i[currRow+1] = ncols_eachRow_i[currRow];
				++currRow;
			}
			ncols_eachRow_i[nRows] = ncols_eachRow_i[nRows-1];
		}
	}
	void matCopy(const int nnz, const vector<unsigned int> &rind, const vector<unsigned int> &cind, const vector<double> &values) {
		// check if setMatSize() was already called
		assert(nRows>0 && nCols>0);

		// allocate memory
		allocateMem(nnz);

		// copy the matrix and set ncols_eachRow_i and diag_index
		for(int i=0; i<nRows+1; ++i)
			ncols_eachRow_i[i] = 0;
		for(int i=0; i<nRows; ++i)
			diag_index[i] = MATCOMPRESED_NOT_FOUND;

		int currRow = 0;
		int count = 0;
		for(int i=0; i<nnz; ++i) {
			if(rind[i]>=nRows)
				break;	// If the number of rows of the Jacobian matrix is greater than nRows, just stop feeding the matrix

			this->rind[i] = rind[i];
			this->cind[i] = cind[i];
			this->values[i] = values[i];

			if(currRow < rind[i]) {
				ncols_eachRow_i[currRow+1] = ncols_eachRow_i[currRow] + count;
				while(++currRow < rind[i]) {
					ncols_eachRow_i[currRow+1] = ncols_eachRow_i[currRow];
				}
				count = 0;
			}
			++count;

			if(currRow == cind[i]) {
				diag_index[currRow] = i;
			}
		}
		if(currRow == nRows-1) {
			ncols_eachRow_i[nRows] = ncols_eachRow_i[nRows-1] + count;
		} else { // check here!!!!!!
			ncols_eachRow_i[currRow+1] = ncols_eachRow_i[currRow] + count;
			while(currRow < nRows-1) {
				ncols_eachRow_i[currRow+1] = ncols_eachRow_i[currRow];
				++currRow;
			}
			ncols_eachRow_i[nRows] = ncols_eachRow_i[nRows-1];
		}
	}

#ifdef USE_MPI_
	void matCopyUGP(const int ncv, const int nbocv_s, const int* nbocv_i, const int* nbocv_v, const int nScal, double (*A)[5][5], double ***AScal) {
		// check if setMatSize() was already called
		assert(nRows>0 && nCols>0);

		// check if matrices are empty
		assert(A != NULL);
		if(nScal>0) {
			assert(AScal != NULL);
			if(mpi_rank==0)
				cout<<"MatComprsed::matCopyUGP() - Scalar cases have NOT been tested!"<<endl;
		} else {
			assert(AScal == NULL);
		}

		// allocate memory
		int m = 5+nScal; //
		int nnzUGP = nbocv_s*m*m;
		allocateMem(nnzUGP);

		// copy the matrix and set ncols_eachRow_i and diag_index
		/*   To see how ncols_eachRow_i and diag_index are constructed, refer to the matCopy() method */
		for(int i=0; i<nRows+1; ++i)
			ncols_eachRow_i[i] = 0;
		for(int i=0; i<nRows; ++i)
			diag_index[i] = MATCOMPRESED_NOT_FOUND;

		int count = 0;
		for (int icv = 0; icv < ncv; ++icv) {
			int noc0 = nbocv_i[icv]; // the index of nbocv_v for the first block in the row
			int noc1 = nbocv_i[icv+1]-1; // the index of nbocv_v for the last block in the row
										/* For example, if we have the following 4 CVs such as
										 *     2 -- 3
										 *     |    |
										 *     0 -- 1 ,
										 * Then, nbocv_v = { (0,1,2),(1,0,3),(2,0,3),(3, 1, 2) }
										 *                    0 1 2   3 4 5   6 7 8   9 10 11
										 *       nbocv_i = { 0,3,6,9,? }
										 */

			assert(noc1 > noc0); // note : this algorithm is based on the assumption that no row is empty (in other words, noc1 == noc0 is not allowed)

			// save rind, cind, and values for diagonal blocks of A
			int startingRowIndex = icv*m;
			int startingColIndex = startingRowIndex;
			matCopyUGP_copyEachIcv(count, noc0, nScal, startingRowIndex, startingColIndex, A, AScal, true);

			// save rind, cind, and values for off-diagonal blocks of A
			for (int noc=noc0+1; noc<=noc1; ++noc) {
				int currIcv = nbocv_v[noc];
				int startingColIndexOffDiag = currIcv*m;
				matCopyUGP_copyEachIcv(count, noc, nScal, startingRowIndex, startingColIndexOffDiag, A, AScal, false);
			}

			// get ncols_eachRow_i
			for (int i=0; i<m; ++i) {
				int currRow = startingRowIndex + i;
				ncols_eachRow_i[currRow+1] = ncols_eachRow_i[currRow] + (noc1-noc0+1)*m;
			}
		}

		assert(count == nnzUGP);
	}

	void matCopyUGP_copyEachIcv(int &count, const int noc, const int nScal, const int startingRowIndex, const int startingColIndex, double (*A)[5][5], double ***AScal, const bool isDiagonalBlock) {
		/*   Note: 1. A and AScal in the JOE code are actually the negative of Jacobian
		 *         2. The size of A is always [nbocv_s][5][5], which means that there is no dependency of A on scalars
		 *            The size of AScal is always [nScal][6][nbocv_s], which means that there is no dependency of AScal on other scalars
		 *         Thus, the Jacobian of an icv should looks like this (let's assume that you have 2 scalars):
		 *           Jac = [  -A00      -A01      -A02      -A03      -A04        0         0
		 *                    -A10      -A11      -A12      -A13      -A14        0         0
		 *                    -A20      -A21      -A22      -A23      -A24        0         0
		 *                    -A30      -A31      -A32      -A33      -A34        0         0
		 *                    -A40      -A41      -A42      -A43      -A44        0         0
		 *                  -AScal0_0 -AScal0_1 -AScal0_2 -AScal0_3 -AScal0_4 -AScal0_5     0
		 *                  -AScal1_0 -AScal1_1 -AScal1_2 -AScal1_3 -AScal1_4     0     -AScal1_5 ]
		 */
		// Jacobian from A
		for (int i = 0; i < 5; i++) { // Jacobian in terms of rho, rhou, rhoE
			int currRow = startingRowIndex + i;
			for (int j = 0; j < 5; j++) {
				this->rind[count]   = currRow;
				this->cind[count]   = startingColIndex + j;
				this->values[count] = -A[noc][i][j];

				// get diag_index
				if(isDiagonalBlock) {
					if(i==j)
						this->diag_index[currRow] = count;
				}

				++count;
			}
			for (int iScal=0; iScal<nScal; ++iScal) { // Jacobian in terms of scalars
				this->rind[count]   = currRow;
				this->cind[count]   = startingColIndex + 5 + iScal;
				this->values[count] = 0.0;
				++count;
			}
		}
		// Jacobian from AScal
		for (int iScal=0; iScal<nScal; ++iScal) {
			int currRow = startingRowIndex + 5 + iScal;
			for (int j=0; j<5; j++) { // Jacobian in terms of rho, rhou, rhoE
				this->rind[count]   = currRow;
				this->cind[count]   = startingColIndex + j;
				this->values[count] = -AScal[iScal][j][noc];
				++count;
			}
			for (int j=0; j<nScal; j++) { // Jacobian in terms of scalars
				this->rind[count] = currRow;
				this->cind[count] = startingColIndex + 5 + j;
				if(j==iScal) {
					this->values[count] = -AScal[iScal][5][noc];

					// get diag_index
					if(isDiagonalBlock) {
						this->diag_index[currRow] = count;
					}
				} else {
					this->values[count] = 0.0;
				}
				++count;
			}
		}
	}

	void matCopyUGPcoupled(const int ncv, const int nbocv_s, const int* nbocv_i, const int* nbocv_v, const int nPriVar, const int nScal, double*** A) {
		// check if setMatSize() was already called
		assert(nRows>0 && nCols>0);

		// allocate memory
		int m = nPriVar+nScal; // note: In most of the cases, nPriVal == 5
		int nnzUGP = nbocv_s*m*m;
		allocateMem(nnzUGP);

		// copy the matrix and set ncols_eachRow_i and diag_index
		/*   To see how ncols_eachRow_i and diag_index are constructed, refer to the matCopy() method */
		for(int i=0; i<nRows+1; ++i)
			ncols_eachRow_i[i] = 0;
		for(int i=0; i<nRows; ++i)
			diag_index[i] = MATCOMPRESED_NOT_FOUND;

		int count = 0;
		for (int icv = 0; icv < ncv; ++icv) {
			int noc0 = nbocv_i[icv]; // the index of nbocv_v for the first block in the row
			int noc1 = nbocv_i[icv+1]-1; // the index of nbocv_v for the last block in the row
										/* For example, if we have the following 4 CVs such as
										 *     2 -- 3
										 *     |    |
										 *     0 -- 1 ,
										 * Then, nbocv_v = { (0,1,2),(1,0,3),(2,0,3),(3,1,2) }
										 *       nbocv_i = { 0,3,6,9 }
										 */

			assert(noc1 > noc0); // note : this algorithm is based on the assumption that no row is empty

			// save rind, cind, and values for diagonal blocks of A
			int startingRowIndex = icv*m;
			int startingColIndex = startingRowIndex;
				for (int i = 0; i < m; i++)
					for (int j = 0; j < m; j++) {
						int currRow = startingRowIndex + i;
						this->rind[count]   = currRow;
						this->cind[count]   = startingColIndex + j;
						this->values[count] = A[noc0][i][j];

						// get diag_index
						if(i==j)
							diag_index[currRow] = count;

						++count;
					}

			// save rind, cind, and values for off-diagonal blocks of A
			for (int noc=noc0+1; noc<=noc1; ++noc) {
				int currIcv = nbocv_v[noc];
				int startingColIndex = currIcv*m;
				for (int i = 0; i < m; i++)
					for (int j = 0; j < m; j++) {
						this->rind[count]   = startingRowIndex + i;
						this->cind[count]   = startingColIndex + j;
						this->values[count] = A[noc][i][j];
						++count;
					}
			}

			// get ncols_eachRow_i
			for (int i=0; i<m; ++i) {
				int currRow = startingRowIndex + i;
				ncols_eachRow_i[currRow+1] = ncols_eachRow_i[currRow] + (noc1-noc0+1)*m;
			}
		}

		assert(count == nnzUGP);
	}
#endif

	void clear() {
		if(rind != NULL) {
			delete [] rind;
			rind = NULL;
		}
		if(cind != NULL) {
			delete [] cind;
			cind = NULL;
		}
		if(values != NULL) {
			delete [] values;
			values = NULL;
		}

		if(ncols_eachRow_i != NULL) {
			delete [] ncols_eachRow_i;
			ncols_eachRow_i = NULL;
		}
		if(diag_index != NULL) {
			delete [] diag_index;
			diag_index = NULL;
		}

		nnz = 0;

		nCols = 0;
		nRows = 0;
	}

	void showMatrix() {
		cout << endl << "***** "<<endl;
		cout << "Size of the Matrix: "<<endl;
		cout << " nnz = "<<get_nnz()<<endl;
		cout << " nRows = "<<get_nRows()<<endl;
		cout << " nCols = "<<get_nCols()<<endl;
		cout <<endl; 

		cout << "Pattern of the Matrix: "<<endl;
		int index = 0;
		for (int row = 0; row < get_nRows(); ++row) {
			cout<<" row "<<row<<": ";
			for (int col = 0; col < get_nCols(); ++col) {
				if (row == (int) get_rind(index) && col == (int) get_cind(index)) {
					cout << " 1";
					++index;
				} else {
					cout << " 0";
				}
			}
			cout << endl;
		}
		if (index != get_nnz())
			cout << "Warning! Total number of non-zero elements shown on the screen (=" << index << ") is not same as nnz(="<< get_nnz() << ")" << endl
		         << "         The last two non-zero element shown on the screen: index=" << index-1 << ", rind=" << get_rind(index-1) << ", cind=" << get_cind(index-1) << endl
		         << "                                                            index=" << index   << ", rind=" << get_rind(index)   << ", cind=" << get_cind(index)   << endl;


		cout << endl;
		cout << "Values of the Matrix: "<<endl;

		{
			int index = 0;
			for (int row = 0; row < get_nRows(); ++row) {
				cout<<" row "<<row<<": ";
				for (int col = 0; col < get_nCols(); ++col) {
					if (row == (int) get_rind(index) && col == (int) get_cind(index)) {
						cout << " "<<get_values(index);
						++index;
					} else {	
						cout << " 0.0";
					}
				}
				cout << endl;
			}
		}

		cout << endl;
		cout << "Structure of the Matrix: "<<endl;
		cout << " ncols_eachRow_i: ";
		for(int i=0; i<get_nRows()+1; ++i)
			cout<<get_ncols_eachRow_i(i)<<"  ";
		cout << endl;
		cout << " diag_index: ";
		for(int i=0; i<get_nRows(); ++i)
			cout<<get_diag_index(i)<<"  ";
		cout << endl;

		cout << endl;
		cout << "***** " << endl <<endl;
	}


	void showAddress() {
		if(rind == NULL) {
			cout<<"**************************"<<endl;
			cout<<"MatComprsed::showAddress()"<<endl;
			cout<<" &rind="<<rind<<endl;
			cout<<" &cind="<<cind<<endl;
			cout<<" &values="<<values<<endl;
			cout<<" &ncols_eachRow_i="<<ncols_eachRow_i<<endl;
			cout<<" &diag_index="<<diag_index<<endl;
			cout<<"**************************"<<endl;
		} else {
			cout<<"MatComprsed::showAddress()"<<endl;
			cout<<" &rind="<<rind<<": first 5 elements=";
			for(int i=0; i<5; ++i) cout<<rind[i]<<"  ";
			cout<<endl;
			cout<<" &cind="<<cind<<": first 5 elements=";
			for(int i=0; i<5; ++i) cout<<cind[i]<<"  ";
			cout<<endl;
			cout<<" &values="<<values<<": first 5 elements=";
			for(int i=0; i<5; ++i) cout<<values[i]<<"  ";
			cout<<endl;
			cout<<" &ncols_eachRow_i="<<ncols_eachRow_i<<": first 5 elements=";
			for(int i=0; i<5; ++i) cout<<ncols_eachRow_i[i]<<"  ";
			cout<<endl;
			cout<<" &diag_index="<<diag_index<<": first 5 elements=";
			for(int i=0; i<5; ++i) cout<<diag_index[i]<<"  ";
			cout<<endl;
		}
	}

protected:
	int nnz;               /* number of nonzeros */
	unsigned int *rind;    /* row indices        */
	unsigned int *cind;    /* column indices     */
	double       *values;  /* values             */

	int *ncols_eachRow_i;
	int *diag_index;

	int nCols, nRows;
};

/*
 * class: MatComprsedSTL
 * ---------------------
 * Matrix in a compressed form (row form) using STL containers
 */
class MatComprsedSTL
{
public:
	MatComprsedSTL() {
		nnz = 0;
		nRows = 0; // Note: this is local (NOT global) matrix size
		nCols = 0; // Note: this is local (NOT global) matrix size

		currMaxRow = 0;

		showDetails = false;
	}
	~MatComprsedSTL() {
		clear();
	}

	bool empty() {
		if (nnz==0 && rind.empty() && cind.empty() && values.empty() && ncols_eachRow_i.empty() && diag_index.empty() && nRows==0 && nCols==0)
			return true;
		else
			return false;
	}

	int get_nnz() {
		return nnz;
	}
	int get_nCols() {
		return nCols;
	}
	int get_nRows() {
		return nRows;
	}
	int get_currMaxRow() {
		return currMaxRow;
	}
	unsigned int get_rind(const size_t index) {
//		assert(index>=0 && index < get_nnz());
		return rind[index];
	}
	unsigned int get_cind(const size_t index) {
//		assert(index>=0 && index < get_nnz());
		return cind[index];
	}
	double get_values(const size_t index) {
//		assert(index>=0 && index < get_nnz());
		return values[index];
	}
	int get_ncols_eachRow_i(const size_t index) {
		return ncols_eachRow_i[index];
	}
	int get_diag_index(const size_t index) {
		assert(index<get_nRows());
		return diag_index[index];
	}
	int get_global_rind(const int index, const int mpi_rank, const int nVars, const int *cvora) {
		int rindLocal = get_rind(index);
		return cvora[mpi_rank]*nVars + rindLocal;
	}
	int get_global_cind(const int index, const int mpi_rank, const int nVars, const int *cvora, const int *cv_gl, const int NcontrolEqns=0, const int ncv_gg=0) {
		int cindLocal = get_cind(index);

		int globalCind = MATCOMPRESED_NOT_FOUND; // This variable will be updated and returned

		// first, get icvGlobal(global icv)
		int icvLocal = (int) (cindLocal/nVars); // icvLocal can be greater than ncv_gg if NcontrolEqns>0

		if(icvLocal<0) { 	// THIS IS AN ERROR!
			printf("ERROR in MatComprsedSTL::get_global_cind(): mpi_rank=%d -- icvLocal=%d < 0 \n", mpi_rank, icvLocal);
			assert(icvLocal >= 0);
		} else if(icvLocal<ncv_gg) {
			int rest = cindLocal%nVars;
			int icvGlobal = cv_gl[icvLocal];
#ifdef SLOW_BUT_ERROR_CHECK
			if(icvGlobal<0 || icvGlobal>cvora[mpi_size]) {
				printf("ERROR in MatComprsedSTL::get_global_cind(): mpi_rank=%d -- icvLocal=%d, icvGlobal=%d (!), cvora[mpi_size]=%d \n", mpi_rank, icvLocal, icvGlobal, cvora[mpi_size]);
				assert(icvGlobal>=0);
				assert(icvGlobal<=cvora[mpi_size]);
			}
#endif
			// then, calculate globalCind using icvGlobal and rest
			globalCind = icvGlobal*nVars + rest;
		} else if (icvLocal<ncv_gg+int(NcontrolEqns/nVars)+1) {
			int iEqn = cindLocal - nVars*ncv_gg;
#ifdef SLOW_BUT_ERROR_CHECK
			assert(NcontrolEqns>0);
			if(iEqn<0 || iEqn>=NcontrolEqns) {
				printf("ERROR in MatComprsedSTL::get_global_cind(): mpi_rank=%d -- iEqn=%d (!), NcontrolEqns=%d, cindLocal=%d, nVars*ncv_gg=%d \n", mpi_rank, iEqn, NcontrolEqns, cindLocal, nVars*ncv_gg);
				assert(iEqn>=0);
				assert(iEqn<NcontrolEqns);
			}
#endif
			// then, calculate globalCind using iEqn, get_nCols(), & NcontrolEqns
			globalCind = cvora[mpi_size]*nVars + iEqn;
		} else { 			// THIS IS AN ERROR!
			printf("ERROR in MatComprsedSTL::get_global_cind(): mpi_rank=%d -- icvLocal=%d >= ncv_gg(=%d)+ControlEqnsBlock(=%d) \n", mpi_rank, icvLocal, ncv_gg, int(NcontrolEqns/nVars)+1);
			assert(icvLocal<ncv_gg+int(NcontrolEqns/nVars)+1);
		}

		return globalCind;
	}

	int get_local_rind(const int globalIndex, const int mpi_rank, const int nVars, const int *cvora) {
		return globalIndex - cvora[mpi_rank]*nVars;
	}
	int get_local_cind(const int globalIndex, const int mpi_rank, const int nVars, const int *cvora, const int *cv_gl, const int NcontrolEqns=0, const int ncv=0, const int ncv_gg=0) {
		int localIndex = MATCOMPRESED_NOT_FOUND; // This variable will be updated and returned

		int icvGlobal = (int) globalIndex/nVars;
		int rest = globalIndex%nVars;

		// First, get icvLocal(local icv)
		int icvLocal = icvGlobal - cvora[mpi_rank]; // For interior CVs, 0<=icvLocal<ncv.

		// There are three possibility: 1. interior CV   2. control Eqn   3. ghost CV
		if(icvLocal>=0 && icvLocal<ncv) {
			localIndex = icvLocal*nVars + rest;
		} else if(globalIndex >= cvora[mpi_size]*nVars) {
			assert(NcontrolEqns > 0);

			int iEqn = globalIndex - cvora[mpi_size]*nVars;
			assert(iEqn < NcontrolEqns);

			localIndex = ncv_gg*nVars + iEqn;
		} else {
			int foundICVindex = MATCOMPRESED_NOT_FOUND;
			for(int icv=ncv; icv<ncv_gg; ++icv) {
				if(cv_gl[icv] == icvGlobal) {
					foundICVindex = icv;
					break;
				}
			}

			if(foundICVindex == MATCOMPRESED_NOT_FOUND) {
				cerr<<"ERROR in MatComprsedSTL::get_local_cind(): cannot find the local column index for the global index="<<globalIndex<<" (mpi_rank="<<mpi_rank<<")"<<endl;
				cerr<<"                                                   (global icv="<<icvGlobal<<", possible range in this core="
						<<cvora[mpi_rank]<<"~"<<cvora[mpi_rank]+ncv<<"(or "<<cvora[mpi_rank]+ncv_gg<<" for the ghost cells)"<<endl;
				assert(false);
			} else {
				localIndex = foundICVindex*nVars + rest;
			}
		}

		return localIndex;
	}

	void set_nnz(const int nnz) {
		this->nnz = nnz;
	}
	void setMatSize(const int nCols, const int nRows) {
		assert( empty() );

		this->nCols = nCols;
		this->nRows = nRows;
	}
	void set_currMaxRow(const int currMaxRow) {
		this->currMaxRow = currMaxRow;
	}
	void resize_ncolsEachRowI_diagIndex() {
		assert(nRows>0 && nCols>0);

		assert(ncols_eachRow_i.empty());
		ncols_eachRow_i.resize(nRows+1);

		assert(diag_index.empty());
		resize_diagIndex();
	}
	void resize_diagIndex() {
		diag_index.resize(nRows, MATCOMPRESED_NOT_FOUND);
	}

	void set_rind(const size_t index, const unsigned int ind) {
		assert(ind < get_nRows());
		rind[index] = ind;
	}
	void set_cind(const size_t index, const unsigned int ind) {
		assert(ind < get_nCols());
		cind[index] = ind;
	}
	void convert_blockCindices(const vector<int> &myNbocv2, const int indexStart, const int nnzBlock, const int nRowsBlock, const int nVars, const int NcontrolEqns=0, const int ncv_gg=0) { // This method only for the MEM_SAVING_ADVAR case
		/* For example,
		 *   Full allocation:            ncv=5, nVars=2, NcontrolEqns=1.
		 *   MEM_SAVING_ADVAR algorithm: myNbocv2={3,1,4} is allocated
		 *
		 *   If a row of the fully-allocated Jacobian matrix is [0.0  0.0  1.0  1.0  0.0  0.0  3.0  3.0  4.0  4.0  9.0]   (including zeros for a demostration purpose),
		 *                                                        0    0    1    1    2    2    3    3    4    4          <-- icv
		 *   the MEM_SAVING_ADVAR algorithm gives [3.0  3.0  1.0  1.0  4.0  4.0  9.0] after calling blockCopy()
		 * Thus, you must convert the column indices so that you can call other routines which were designed to be compatible with the fully-allocated case.
		 *
		 * Note that this routine automatically reset the column indices only in the very last "nRowsBlock" rows
		 */

		assert(indexStart+nnzBlock <= get_nnz()); // Check: You must call this method after calling blockCopy()
		assert(get_rind(indexStart)+nRowsBlock <= get_nRows());
#ifdef SLOW_BUT_ERROR_CHECK
		if(nRowsBlock != nVars)
			cout<<"WARNING in MatComprsedSTL::convert_blockCindices(): nRowsBlock(="<<nRowsBlock<<") is not equal to nVars(="<<nVars<<")"<<endl;
#endif

		// First, reset diag_index for the block
		int rindStart = get_rind(indexStart);
		for(size_t iRow=rindStart; iRow<rindStart+nRowsBlock; ++iRow)
			set_diag_index(iRow, MATCOMPRESED_NOT_FOUND);

		// Reset cind and diag_index
		for(size_t i=indexStart; i<indexStart+nnzBlock; ++i) {
			int fake_cind = get_cind(i); // Note: In the above example, fake_cind = {0 1 2 3 4 5 6}
			int fake_icv  = int(fake_cind/nVars); // Get the fake icv due to the local assignment
			                                      // In the above example, fake_icv = {0 0 1 1 2 2 3}
			int n_myNbocv2 = myNbocv2.size();

			unsigned int tempCind;
			if(fake_icv < n_myNbocv2) {
				int rest = fake_cind%nVars;
				tempCind = myNbocv2[fake_icv]*nVars+rest; // myNbocv2[fake_icv] is the real icv number
				                                          // In the above example, myNbocv2[fake_icv] = {3 3 1 1 4 4}
				set_cind(i, tempCind);
			} else {
				if(NcontrolEqns==0) {
					cout<<"ERROR in MatComprsedSTL::convert_blockCindices(): try to access to an index(="<<fake_icv<<") > size of the vector(="<<n_myNbocv2<<")"<<endl;
					assert(false);
				}
				int iEqn = fake_cind - nVars*n_myNbocv2; // Since NcontrolEqns>=Nvars, you should not use "rest" instead of "iEqn"
				                                         // Note: Actually, you must add 1 to iEqn (i.e. iEqn = fake_cind - nVars*n_myNbocv2 + 1),
				                                         //       but in that case you will also subtract 1 from tempCind anyway.
				assert(iEqn>=0);
				tempCind = ncv_gg*nVars + iEqn;
				set_cind(i, tempCind);
			}

			if(get_rind(i) == tempCind) {
				set_diag_index(tempCind, i);
			}
		}

		// Sort in the order of ascending cind
		for(size_t iRow=rindStart; iRow<rindStart+nRowsBlock; ++iRow)
			sortCind_ascendingOrder(iRow);
	}
	void set_values(const size_t index, const double val) {
		values[index] = val;
	}
	void set_ncols_eachRow_i(const size_t index, const int ind) {
		ncols_eachRow_i[index] = ind;
	}
	void set_diag_index(const size_t index, const int ind) {
		diag_index[index] = ind;
	}

	void set_showDetails(const bool input) {
		showDetails = input;
	}

	void add_rind(const unsigned int ind) {
		assert(ind < get_nRows());
		rind.push_back(ind);
	}
	void add_cind(const unsigned int ind) {
		assert(ind < get_nCols());
		cind.push_back(ind);
	}
	void add_values(const double val) {
		values.push_back(val);
	}

	void initMat() {
		resize_ncolsEachRowI_diagIndex();
	}
	void initMat(const int nnz) {
		set_nnz(nnz);

		assert(rind.empty() && cind.empty() && values.empty());
		rind.resize(nnz);
		cind.resize(nnz);
		values.resize(nnz);

		resize_ncolsEachRowI_diagIndex();
	}

	void matCopy(const int nnz, const unsigned int *rind, const unsigned int *cind, const double *values) {
		// check if setMatSize() was already called
		assert(nRows>0 && nCols>0);

		// allocate memory
		initMat(nnz);

		// copy the matrix and set ncols_eachRow_i and diag_index
		/*   ncols_eachRow_i[row] indicates the accumulated number of the first non-zero element of the row in rind, cind, and values
		 *   e.g. A = [0 0 0 0
		 *             0 1 1 0
		 *             1 1 0 1
		 *             0 0 0 0]
		 *        ncols_eachRow_i[5] = {0, 0, 2, 5, 5}
		 *   Thus, the number of non-zeros in a row: ncols_eachRow_i[row+1] - ncols_eachRow_i[row]
		 */
		/*
		 *   diag_index[row] indicates the index of the diagonal element if the row in rind, cind, and values
		 *   e.g. A = [9 9 0 0
		 *             0 9 9 0
		 *             0 0 0 0
		 *             0 0 0 9]
		 *        rind[5] = {0, 0, 1, 1, 3}
		 *        cind[5] = {0, 1, 1, 2, 3}
		 *        values[5] = {9, 9, 9, 9, 9}
		 *
		 *        In the case that MATCOMPRESED_NOT_FOUND is defined as -1, diag_index[4] = {0, 2, -1, 4}
		 */
		int currRow = 0;
		int count = 0;
		for(int i=0; i<nnz; ++i) {
			set_rind(i, rind[i]);
			set_cind(i, cind[i]);
			set_values(i, values[i]);

			if(currRow < rind[i]) {
				set_ncols_eachRow_i(currRow+1, ncols_eachRow_i[currRow] + count);
				while(++currRow < rind[i]) {
					set_ncols_eachRow_i(currRow+1, ncols_eachRow_i[currRow]);
				}
				count = 0;
			}
			++count;

			if(currRow == cind[i]) {
				set_diag_index(currRow, i);
			}
		}
		if(currRow == nRows-1) {
			set_ncols_eachRow_i(nRows, get_ncols_eachRow_i(nRows-1) + count);
		} else { // check here!!!!!!
			set_ncols_eachRow_i(currRow+1, get_ncols_eachRow_i(currRow) + count);
			set_ncols_eachRow_i(nRows, get_ncols_eachRow_i(nRows-1));
		}
	}
	void matCopy(const int nnz, const vector<unsigned int> &rind, const vector<unsigned int> &cind, const vector<double> &values) {
		// check if setMatSize() was already called
		assert(nRows>0 && nCols>0);

		// allocate memory
		initMat(nnz);

		// copy the matrix and set ncols_eachRow_i and diag_index
		int currRow = 0;
		int count = 0;
		for(int i=0; i<nnz; ++i) {
			set_rind(i, rind[i]);
			set_cind(i, cind[i]);
			set_values(i, values[i]);

			if(currRow < rind[i]) {
				set_ncols_eachRow_i(currRow+1, ncols_eachRow_i[currRow] + count);
				while(++currRow < rind[i]) {
					set_ncols_eachRow_i(currRow+1, ncols_eachRow_i[currRow]);
				}
				count = 0;
			}
			++count;

			if(currRow == cind[i]) {
				set_diag_index(currRow, i);
			}
		}
		if(currRow == nRows-1) {
			set_ncols_eachRow_i(nRows, get_ncols_eachRow_i(nRows-1) + count);
		} else { // check here!!!!!!
			set_ncols_eachRow_i(currRow+1, get_ncols_eachRow_i(currRow) + count);
			set_ncols_eachRow_i(nRows, get_ncols_eachRow_i(nRows-1));
		}
	}

	void rowCopy(const int nnzRow, const unsigned int currRow, const unsigned int *cind, const double *values) {
		// check if resize_ncolsEachRowI_diagIndex() was already called
		assert(ncols_eachRow_i.size()>0 && diag_index.size()>0);

		// allocate memory
		this->cind.resize(nnz+nnzRow);
		this->rind.resize(nnz+nnzRow);
		this->values.resize(nnz+nnzRow);

		// copy the row and set ncols_eachRow_i and diag_index
		/*   ncols_eachRow_i[row] indicates the accumulated number of the first non-zero element of the row in rind, cind, and values
		 *   e.g. A = [0 0 0 0
		 *             0 1 1 0
		 *             1 1 0 1
		 *             0 0 0 0]
		 *        ncols_eachRow_i[5] = {0, 0, 2, 5, 5}
		 *   Thus, the number of non-zeros in a row: ncols_eachRow_i[row+1] - ncols_eachRow_i[row]
		 */
		/*
		 *   diag_index[row] indicates the index of the diagonal element if the row in rind, cind, and values
		 *   e.g. A = [9 9 0 0
		 *             0 9 9 0
		 *             0 0 0 0
		 *             0 0 0 9]
		 *        rind[5] = {0, 0, 1, 1, 3}
		 *        cind[5] = {0, 1, 1, 2, 3}
		 *        values[5] = {9, 9, 9, 9, 9}
		 *
		 *        In the case that MATCOMPRESED_NOT_FOUND is defined as -1, diag_index[4] = {0, 2, -1, 4}
		 */
		while(++currMaxRow < currRow) {
			set_ncols_eachRow_i(currMaxRow+1, get_ncols_eachRow_i(currMaxRow));
		}

		int count = 0;
		for(int i=0; i<nnzRow; ++i) {
			size_t currIndex = get_nnz()+i;
			set_rind(currIndex, currRow);
			set_cind(currIndex, cind[i]);
			set_values(currIndex, values[i]);
			++count;

			if(currRow == cind[i])
				set_diag_index(currRow, currIndex);
		}

		set_ncols_eachRow_i(currRow+1, get_ncols_eachRow_i(currRow) + count);

//		// Sort in the order of ascending cind
//		sortCind_ascendingOrder(currRow);

		currMaxRow=currRow;

		set_nnz(get_nnz()+nnzRow);
	}

	void blockCopy(const int nnzBlock, const int firstRindReal, const unsigned int *rindBlock, const unsigned int *cind, const double *values) {
		// check if nnzBlock, rindBlock, cind, and values are correct
		if(nnzBlock==0 || rindBlock==NULL || cind==NULL || values==NULL) {
			cout<<"Error in MatComprsedSTL::blockCopy() at mpi="<<mpi_rank<<" => passed block is empty!"<<endl;
			cout<<"                                                 nnzBlock="<<nnzBlock<<"  &rind="<<rindBlock<<"  &cind="<<cind<<endl;
			throw(-1);
		}

		int currRowReal = firstRindReal+rindBlock[0];

		// check if resize_ncolsEachRowI_diagIndex() was already called properly
		if(ncols_eachRow_i.size()<currRowReal || diag_index.size()<currRowReal) {
			printf("  ncols_eachRow_i.size()=%d, currRowReal=%d, diag_index.size()=%d, currRowReal=%d \n", (int)ncols_eachRow_i.size(), currRowReal, (int)diag_index.size(), currRowReal);
		}
		assert(ncols_eachRow_i.size()>currRowReal && diag_index.size()>currRowReal);

		// check possible mistakes
		assert(currMaxRow <= currRowReal);

		// allocate memory
		this->cind.resize(nnz+nnzBlock);
		this->rind.resize(nnz+nnzBlock);
		this->values.resize(nnz+nnzBlock);

		// copy the matrix and set ncols_eachRow_i and diag_index
		/*   ncols_eachRow_i[row] indicates the accumulated number of the first non-zero element of the row in rind, cind, and values
		 *   e.g. A = [0 0 0 0
		 *             0 1 1 0
		 *             1 1 0 1
		 *             0 0 0 0]
		 *        ncols_eachRow_i[5] = {0, 0, 2, 5, 5}
		 *   Thus, the number of non-zeros in a row: ncols_eachRow_i[row+1] - ncols_eachRow_i[row]
		 */
		/*
		 *   diag_index[row] indicates the index of the diagonal element if the row in rind, cind, and values
		 *   e.g. A = [9 9 0 0
		 *             0 9 9 0
		 *             0 0 0 0
		 *             0 0 0 9]
		 *        rind[5] = {0, 0, 1, 1, 3}
		 *        cind[5] = {0, 1, 1, 2, 3}
		 *        values[5] = {9, 9, 9, 9, 9}
		 *
		 *        In the case that MATCOMPRESED_NOT_FOUND is defined as -1, diag_index[4] = {0, 2, -1, 4}
		 */
		if(currMaxRow < currRowReal) {
			while(++currMaxRow < currRowReal) {
				set_ncols_eachRow_i(currMaxRow+1, get_ncols_eachRow_i(currMaxRow));
			} //Thus, currMaxRow == firstRindReal+rindBlock[0] after this while() loop
		}

		int count = 0;
		for(int i=0; i<nnzBlock; ++i) {
			int currIndexReal = get_nnz()+i;
			currRowReal = firstRindReal+rindBlock[i];
			set_rind(currIndexReal, currRowReal);
			set_cind(currIndexReal, cind[i]);
			set_values(currIndexReal, values[i]);

			if(currMaxRow < currRowReal) { // i.e. new row
				set_ncols_eachRow_i(currMaxRow+1, ncols_eachRow_i[currMaxRow] + count); // set the number of elements in the previous row

				// Take care of the case that a row is empty
				while(++currMaxRow < currRowReal) {
					set_ncols_eachRow_i(currMaxRow+1, ncols_eachRow_i[currMaxRow]);
				} //Thus, currMaxRow == firstRindReal+rindBlock[i] after this while() loop

				count = 0;
			}
			++count;

			if(currRowReal == cind[i]) {
				set_diag_index(currRowReal, currIndexReal);
			}
		}
		set_ncols_eachRow_i(currMaxRow+1, get_ncols_eachRow_i(currMaxRow) + count);
		currMaxRow = currRowReal; // actually, this line is irrelevant

		set_nnz(get_nnz()+nnzBlock);
	}

	void sortCind_ascendingOrder(const int rindGiven) {
		int ncolsLocal = get_ncols_eachRow_i(rindGiven+1) - get_ncols_eachRow_i(rindGiven);
		assert(ncolsLocal >= 0);
		if(ncolsLocal == 0) {
			cout<<"WARNING in sortCind_ascendingOrder(): The number of columns in the "<<rindGiven<<"th row is ZERO"<<endl;
			return;
		}

		int startingIndex = get_ncols_eachRow_i(rindGiven);
		int endingIndex   = get_ncols_eachRow_i(rindGiven+1); // Actual end-index + 1
		assert(ncolsLocal>=0);

		// Reset diag_index
		set_diag_index(rindGiven, MATCOMPRESED_NOT_FOUND);

		// Find the sorting order
		//   For example, cind = [..., 32, 71, 12, 45, ...]
		//                             10  11  12  13        <-- index
		//                              0   1   2   3        <-- local index
		vector<pair<size_t, unsigned int> > cindPairLocal(ncolsLocal);  // (index, cind)

		for(size_t i=0; i<ncolsLocal; ++i) {
			cindPairLocal[i] = std::make_pair(i, get_cind(startingIndex + i)); // For example, cindPairLocal = [(0,32), (1,71), (2,12), (3,45)]
		}

		std::sort(cindPairLocal.begin(), cindPairLocal.end(), comparisonObject); // For example, cindPairLocal becomes [(2,12), (0,32), (3,45), (1,71)]

		// Reorder
		vector<unsigned int> cindLocal(  cind.begin()+startingIndex,   cind.begin()+endingIndex); // For example, cindLocal = [32, 71, 12, 45]
		vector<unsigned int> rindLocal(  rind.begin()+startingIndex,   rind.begin()+endingIndex);
		vector<double>       valuesLocal(values.begin()+startingIndex, values.begin()+endingIndex);

		size_t indexNew = startingIndex; // For example, 10
		for(size_t i=0; i<ncolsLocal; ++i) {
			size_t iSaved = cindPairLocal[i].first; // For example, [2, 0, 3, 1]

			set_cind(indexNew, cindLocal[iSaved]); // For example, set_cind([10, 11, 12, 13], [12, 32, 45, 71]);
			assert(rindLocal[iSaved] == rindGiven);
			set_values(indexNew, valuesLocal[iSaved]);

			// Reset diag_index
			if(rindGiven == cindLocal[iSaved]) {
				set_diag_index(rindGiven, indexNew);
			}

			++indexNew; // For example, [10, 11, 12, 13]
		}
	}

	void finalizeMat() {
		if(currMaxRow+1 < get_nRows()) {
			while(++currMaxRow < get_nRows()) {
				set_ncols_eachRow_i(currMaxRow+1, ncols_eachRow_i[currMaxRow]);
			}
		}
	}

#ifdef USE_GERSHGORIN
	/* Estimate eigenvalues using the Gershgorin circle theorem
	 * Note: Be careful about the end indices --
	 *       e.g. If you want to calculate the Gershgorin disks for the whole matrix,
	 *            you should call calcGershgorinDisk(gershgorinDisk, 0, nRows-1, 0, nnz-1);
	 */
	void calcGershgorinDisks(vector<std::pair<double, double> >& gershgorinDisk,
			const int startRow, const int endRow, const int startIndex, const int endIndex) {
		assert(startRow>=0 && endRow>=startRow && startIndex>=0 && endIndex>=startIndex);
		assert(endRow<get_nRows() && endIndex<get_nnz());

		// Resize the Gershgorin disk vector
		assert(gershgorinDisk.empty());
		gershgorinDisk.resize(endRow-startRow+1, std::make_pair(0.0, 0.0));

		// Reinitialize with (0.0, 0.0)
		for(size_t iRow=startRow; iRow<=endRow; ++iRow) {
			size_t vecIndex = iRow-startRow;
			gershgorinDisk[vecIndex] = std::make_pair(0.0, 0.0); // Initialize all the rows with (0.0, 0.0)
		}

		// Iterative over i = startIndex : endIndex
		double diskCenter = 0.0;
		double diskRadius = 0.0;
		int currRow       = startRow; // Start from the first element in the current block
		for(size_t i=startIndex; i<=endIndex; ++i) {
			int cindTemp = get_cind(i);
			int rindTemp = get_rind(i);
			if(currRow != rindTemp) { // If a new row comes out,
				size_t vecIndex = currRow-startRow;
				gershgorinDisk[vecIndex] = std::make_pair(diskCenter, diskRadius); // Update the disk with the values that have been calculated so far

				currRow = rindTemp;

				diskCenter = 0.0;  // Initialize the center of the disk as zero whenever a new row comes out
				diskRadius = 0.0;  // Initialize the radius of the disk as zero whenever a new row comes out
			}

			if(cindTemp == rindTemp) {
				diskCenter = get_values(i); // Center of the disk is the value of the diagonal entry
			} else {
				diskRadius += fabs(get_values(i)); // Radius of the disk is the sum of the absolute values of the non-diagonal entries
			}
		}
		size_t vecIndex = currRow-startRow;
		gershgorinDisk[vecIndex] = std::make_pair(diskCenter, diskRadius); // Update for the last row
	}
	void calcGershgorinDisks(vector<std::pair<double, double> >& gershgorinDisk) {
		calcGershgorinDisks(gershgorinDisk, 0, nRows-1, 0, nnz-1);
	}
#endif

	void clear() {
		rind.clear();
		cind.clear();
		values.clear();

		ncols_eachRow_i.clear();
		diag_index.clear();

		set_nnz(0);

		nCols = 0;
		nRows = 0;

		set_currMaxRow(0);
	}

	void showMatrix() {
		cout << endl << "***** "<<endl;
		cout << "Size of the Matrix: "<<endl;
		cout << " nnz = "<<get_nnz()<<endl;
		cout << " nRows = "<<get_nRows()<<endl;
		cout << " nCols = "<<get_nCols()<<endl;
		cout <<endl;

		cout << "Pattern of the Matrix: "<<endl;
		int index = 0;
		for (int row = 0; row < get_nRows(); ++row) {
			cout<<" row "<<row<<": ";
			for (int col = 0; col < get_nCols(); ++col) {
				if (row == (int) get_rind(index) && col == (int) get_cind(index)) {
					cout << " 1";
					++index;
				} else {
					cout << " 0";
				}
			}
			cout << endl;
		}
		if (index != get_nnz())
			cout << "Warning! Total number of non-zero elements shown on the screen (=" << index << ") is not same as nnz(="<< get_nnz() << ")" << endl
			     << "         The last two non-zero element shown on the screen: index=" << index-1 << ", rind=" << get_rind(index-1) << ", cind=" << get_cind(index-1) << endl
			     << "                                                            index=" << index   << ", rind=" << get_rind(index)   << ", cind=" << get_cind(index)   << endl;

		cout << endl;
		cout << "Values of the Matrix: "<<endl;

		{
			int index = 0;
			for (int row = 0; row < get_nRows(); ++row) {
				cout<<" row "<<row<<": ";
				for (int col = 0; col < get_nCols(); ++col) {
					if (row == (int) get_rind(index) && col == (int) get_cind(index)) {
						cout << " "<<get_values(index);
						++index;
					} else {
						cout << " 0.0";
					}
				}
				cout << endl;
			}
		}

		cout << endl;
		cout << "Structure of the Matrix: "<<endl;
		cout << " ncols_eachRow_i: ";
		for(int i=0; i<get_nRows()+1; ++i)
			cout<<get_ncols_eachRow_i(i)<<"  ";
		cout << endl;
		cout << " diag_index: ";
		for(int i=0; i<get_nRows(); ++i)
			cout<<get_diag_index(i)<<"  ";
		cout << endl;

#ifdef USE_GERSHGORIN
		vector<std::pair<double, double> > gershgorinDisk;
		calcGershgorinDisks(gershgorinDisk);
		cout << " Gershgorin disk (center, radius):";
		for(int i=0; i<get_nRows(); ++i)
			cout<<" ("<<gershgorinDisk[i].first<<","<<gershgorinDisk[i].second<<") ";
		cout << endl;
#endif
		cout << endl;
		cout << "***** " << endl <<endl;
	}

protected:
	int nnz;               		  /* number of nonzeros */
	vector<unsigned int> rind;    /* row indices        */
	vector<unsigned int> cind;    /* column indices     */
	vector<double>       values;  /* values             */

	vector<int> ncols_eachRow_i;
	vector<int> diag_index;

	int nCols, nRows;  // Note: this is local (NOT global) matrix size
	int currMaxRow;

	struct ComparisonStruct {
		bool operator() (std::pair<size_t, unsigned int> one, std::pair<size_t, unsigned int> two) {
			return (one.second < two.second);
		}
	};
	ComparisonStruct comparisonObject;

	bool showDetails;
};

#endif /* MATCOMPRSED_H_ */
