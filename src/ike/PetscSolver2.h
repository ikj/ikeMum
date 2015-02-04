/*
 * PetscSolver2.h
 *
 *  Created on: Oct 24, 2011
 *      Author: ikj
 */

#ifndef PETSCSOLVER2_H_
#define PETSCSOLVER2_H_

#include "Logging.h"
using namespace logging;
#include "petscksp.h"
#include "petscmat.h"

#include "MatComprsed.h"

#ifndef WITH_PETSC
#define WITH_PETSC
#endif

#ifndef MACHINE_EPS
#define MACHINE_EPS 2.0e-16 // Note: machine epsilon for double precision = 1.11e-16
#endif

#define PETSC_DEBUG_LEVEL 1

#define CALC_SINGULAR_VALUES

static bool ShowPetscMatrixMatlab = false; // If this is true, the matrix will be shown on the screen in the MATLAB format
                                           // right after the matrix is constructed in setLinSysForPetsc()
static bool ShowPetscRhsMatlab = false; // If this is true, the RHS vector will be shown on the screen
                                        // right after the vector is constructed in setLinSysForPetsc()
static bool ShowPetscXMatlab = false; // If this is true, the solution vector will be shown on the screen
		                              // right after the vector is constructed in setLinSysForPetsc()

//#define USE_SLEPC_WITH_PETSC
#ifdef USE_SLEPC_WITH_PETSC
	//#include "slepcsys.h"
	#include <slepceps.h>
	#include <slepcsvd.h>
#endif

/**********************
 * GLOBAL OBJECTS AND FUNCS
 **********************/

/*
 * Struct: PetscTol
 */
struct PetscTol {
	PetscTol() {
		convergMonitorInterv = -1;
		kspMonitorHistory.clear();
	}
	~PetscTol() { }

	double rtol, abstol, dtol;
	int maxits;
	/*  rtol   = the relative convergence tolerance (relative decrease in the residual norm)
	 *  abstol = the absolute convergence tolerance (absolute size of the residual norm)
	 *  dtol   = the divergence tolerance (amount residual can increase before KSPDefaultConverged() concludes that the method is diverging)
	 *  maxits = maximum number of iterations to use
	 */

	vector<pair<int,double> > kspMonitorHistory;
	int convergMonitorInterv;
	double initResid;

	int psalcStep, newtonIter;
};


/*
 * Function: MyKSPMonitor
 * ----------------------
 * This is a user-defined routine for monitoring the KSP iterative solvers: return a vector that contains the history of residual
 *
 * Input Parameters:
 *  ksp   - iterative context
 *  n     - iteration number
 *  rnorm - 2-norm (preconditioned) residual value (may be estimated)
 *  dummy - optional user-defined monitor context
 *
 * Return:
 *  kspMonitorStat - pair of (iteration number, residual)
 */
inline PetscErrorCode MyKSPMonitor(KSP ksp, PetscInt nIter, PetscReal rNorm, void *dummy) {
	PetscTol *KspTol = static_cast<PetscTol*>(dummy);

	Vec resid;
	PetscReal trueNorm;

	int monitorInterv = KspTol->convergMonitorInterv;
	if(monitorInterv>0) {
		double absResid = KspTol->abstol;
		double relResid = KspTol->rtol;
		int maxIter = KspTol->maxits;

		if(nIter==0) {
			KspTol->initResid = rNorm;
			KspTol->kspMonitorHistory.clear();
		}

		if(nIter==0 || nIter%monitorInterv==0 || nIter==maxIter) {
			if(PETSC_DEBUG_LEVEL < 1)
				KspTol->kspMonitorHistory.push_back(std::make_pair(nIter, rNorm));
			else {
				// Calculate the true norm of the (left-preconditioned) residual
				KSPBuildResidual(ksp, NULL, NULL, &resid); // This can be expensive for some iterative solvers such as GMRes
				VecNorm(resid, NORM_2, &trueNorm);
				VecDestroy(&resid);

				KspTol->kspMonitorHistory.push_back(std::make_pair(nIter, trueNorm));
			}
		}
	}

	return 0;
}

/**********************
 * PETSCSOLVER2 OBJECT
 **********************/

#ifdef WITH_PETSC
/*
 * Class: PetscSolver2
 * -------------------
 * Two-layer-ghost version of PetscSolver
 * Original code: PetscSolver in PetscSolver.h
 */
class PetscSolver2
{
protected:
	KSP ksp;        // linear solver context
	PC pc;          // preconditioner context

	Vec x_, b_;     // solution, residual vector
	Mat A_;         // implicit operator matrix

	bool hasInitialGuess;

	PetscViewer viewer;
	//VecScatter scatterContext;

	PetscTol KspTols;

public:
	/*
	 * Constructor
	 */
	PetscSolver2(int *cvora, int *nbocv2_i, vector<int> &nbocv2_v, int m, bool hasInitialGuess=false, int NcontrolParams=0, int monitorConvergInterval=0) {
		KspTols.convergMonitorInterv = monitorConvergInterval;

	#if PETSC_DEBUG_LEVEL > 1
		if(mpi_rank==0) cout<<"PetscSolver2()"<<endl;
		if(mpi_rank==0 && monitorConvergInterval>0) cout<<"PETSc convergence history will be saved in a vector<pair<int, double>>"<<endl;
	#endif
		initPetscSolver(cvora, nbocv2_i, nbocv2_v, m, "ASM", false, hasInitialGuess, NcontrolParams); // By default, use the Additive Schwarz pre-conditioner and do not use PCreuse
	}

	PetscSolver2(int *cvora, int *nbocv2_i, vector<int> &nbocv2_v, int m, string pcType, int pcLevels, bool pcReuse, bool hasInitialGuess=false, int NcontrolParams=0, int monitorConvergInterval=0) {
		KspTols.convergMonitorInterv = monitorConvergInterval;

	#if PETSC_DEBUG_LEVEL > 1
		if(mpi_rank==0) cout<<"PetscSolver2()"<<endl;
		if(mpi_rank==0 && monitorConvergInterval>0) cout<<"PETSc convergence history will be saved in a vector<pair<int, double>>"<<endl;
	#endif
		initPetscSolver(cvora, nbocv2_i, nbocv2_v, m, pcType, pcReuse, hasInitialGuess, NcontrolParams, pcLevels); // Note: pcLevels is required only for the ILU pre-conditioner
	}

	/*
	 * Destructor
	 */
	~PetscSolver2() {
#if PETSC_DEBUG_LEVEL > 1
		if(mpi_rank==0)
			cout<<"~PetscSolver2()"<<endl;
#endif
		finalizePetscSolver();
	}

	/*
	 * Method: setTresholds
	 * --------------------
	 * set tolerances for the KSP linear solver (Why the name is Tresholds (not Thresholds)? I don't know.. )
	 */
	void setTresholds(double zeroAbs, double zeroRel, int maxIter) {
		KspTols.abstol = zeroAbs;
		KspTols.rtol = zeroRel;
		KspTols.dtol = 1.0e10;
		KspTols.maxits = maxIter;

		KSPSetTolerances(ksp, zeroRel, zeroAbs, 1.e10, maxIter);
		/*
		 * PetscErrorCode KSPSetTolerances(KSP ksp,PetscReal rtol,PetscReal abstol,PetscReal dtol,PetscInt maxits)
		 *   ksp    = the Krylov subspace context
		 *   rtol   = the relative convergence tolerance (relative decrease in the residual norm)
		 *   abstol = the absolute convergence tolerance (absolute size of the residual norm)
		 *   dtol   = the divergence tolerance (amount residual can increase before KSPDefaultConverged() concludes that the method is diverging)
		 *   maxits = maximum number of iterations to use
		 */
	}

	/*
	 * Method: setGmresRestart()
	 * -------------------------
	 * Set number of iterations at which GMRes restarts.
	 * Note: The larger the value of m for GMRes(m), the fewer iterations are required for convergence.
 	 *       Moreover, a large enough m can reduce the impediment to superlinear convergence and may be required to avoid stalling.
	 *       However, if m is too large, the goal of restarting as a means of reducing computational and storage costs is negated.
	 *       It was also recently shown that a smaller m can actually result in fewer iterations for some problems (See M. Elermann et al., JCAM, 2000).
	 *       PETSC default = 30
	 */
	void setGmresRestart(const int gmresRestart) {
		KSPGMRESSetRestart(ksp, gmresRestart);
	}

	/*
	 * Method: setLinSysForPetsc
	 * -------------------------
	 * Set the linear system: basically upload the matrix to the memory as the PETSC multi-core form
	 * Called by the solveGMRES_solver() method
	 * note: 1. Matrix L = [     A
	 *                      vecC'   d ]
	 *
	 *       2. If NcontrolParams==0, actually ncv_gg doesn't matter
	 *       3. Matrix structures (MatComprsed and MatComprsed) are defined in MatComprsed.h
	 */
	template <class MatT>
	void setLinSysForPetsc(MatT &A, const double *rhs, const double *xInit, const int *cvora, const int *cv_gl, const int m,
			const int NcontrolParams = 0, const int ncv_gg = 0, double **vecC = NULL, const double *d = NULL) {
// Note: due to ghost cells, A.get_nCols() is neigher cvora[mpi_size]*m+NcontrolParams nor (cvora[mpi_rank+1]-cvora[mpi_rank])*m+NcontrolParams).
//		if(A.get_nCols() != (cvora[mpi_rank+1]-cvora[mpi_rank])*m+NcontrolParams) {
//			if(mpi_rank==0)
//				cerr<<"PetscSolver2::setLinSysForPetsc(): A.get_nCols()="<<A.get_nCols()<<", cvora[mpi_rank]="<<cvora[mpi_rank]<<", m="<<m<<", NcontrolParams="<<NcontrolParams<<endl;
//			assert(A.get_nCols() != (cvora[mpi_rank+1]-cvora[mpi_rank])*m+NcontrolParams);
//		}
		if(NcontrolParams==0)
			assert(vecC==NULL && d==NULL);
		else
			assert(vecC!=NULL && d!=NULL);

		int nno = (int) (A.get_nRows()/m); // local number of CVs
		assert( nno == cvora[mpi_rank+1] - cvora[mpi_rank] );

		// assemble rhs
		if(rhs != NULL) { // In the case of not solving a linear system (e.g. Eigen-decomposition), rhs is not required and thus NULL
			int rowGlobal;
			for (int rowLocal=0; rowLocal<A.get_nRows(); ++rowLocal) {
				rowGlobal = cvora[mpi_rank]*m + rowLocal;
				VecSetValues(b_, 1, &rowGlobal, &rhs[rowLocal], INSERT_VALUES);
			}
			if(mpi_rank == mpi_size-1) {
				rowGlobal = cvora[mpi_size]*m;
				for(int iParam=0; iParam<NcontrolParams; ++iParam)
					VecSetValues(b_, 1, &rowGlobal, &rhs[A.get_nRows()+iParam], INSERT_VALUES); // add one more only for mpi_rank == mpi_size-1
			}

			VecAssemblyBegin(b_);
			VecAssemblyEnd(b_);

			if(ShowPetscRhsMatlab)
				VecView(b_,PETSC_VIEWER_STDOUT_WORLD);
		}

		// assemble xInit
		if(hasInitialGuess) {
			if(xInit == NULL) {
				if(mpi_rank==0)
					cerr<<"ERROR PetscSolver2::setLinSysForPetsc(): Initial guess is not provided by the user."<<endl
						<<"                                         If you dont want to use an initial guess,"<<endl
						<<"                                         you must give hasInitialGuess=FALSE when calling PetscSolver2()"<<endl
						<<"                                         and make phiInit=NULL when calling solveGMRES()"<<endl;
				assert(xInit != NULL);
			}
		} else {
			if(xInit != NULL)
				cout<<"WARNING PetscSolver2::setLinSysForPetsc(): xInit is NOT null even though hasInitialGuess==FALSE."<<endl
				    <<"                                           Did you make this by intention?"<<endl;
		}

		if(xInit != NULL) { // In the case that the initial guess for x_ cannot be estimated, x_ is not required and thus NULL
			int rowGlobal;
			for (int rowLocal=0; rowLocal<A.get_nRows(); ++rowLocal) {
				rowGlobal = cvora[mpi_rank]*m + rowLocal;
				if(fabs(xInit[rowLocal]) < 100*MACHINE_EPS)
					VecSetValue(x_, rowGlobal, 0.0, INSERT_VALUES);
				else
					VecSetValues(x_, 1, &rowGlobal, &xInit[rowLocal], INSERT_VALUES);
			}
			if(mpi_rank == mpi_size-1) {
				rowGlobal = cvora[mpi_size]*m;
				for(int iParam=0; iParam<NcontrolParams; ++iParam) {
					if(fabs(xInit[A.get_nRows()+iParam]) < 100*MACHINE_EPS)
						VecSetValue(x_, rowGlobal, 0.0, INSERT_VALUES);
					else
						VecSetValues(x_, 1, &rowGlobal, &xInit[A.get_nRows()+iParam], INSERT_VALUES); // add one more only for mpi_rank == mpi_size-1
				}
			}

			VecAssemblyBegin(x_);
			VecAssemblyEnd(x_);

			if(ShowPetscXMatlab)
				VecView(x_, PETSC_VIEWER_STDOUT_WORLD);
		}

		// Build matrix
		// 1. Matrix A part
		int myEmptyRows = 0;
		int emptyRows;
		for (int rowLocal=0; rowLocal<A.get_nRows(); ++rowLocal) {
			int ncols_eachRow = A.get_ncols_eachRow_i(rowLocal+1) - A.get_ncols_eachRow_i(rowLocal); // note: number of columns in each row of matrix A

			if(ncols_eachRow > 0) {
				// Arrays
				int *cindGlobalArr_eachRow = new int [ncols_eachRow]; // global column index array
				double *valuesArr_eachRow = new double [ncols_eachRow]; // values array

				// Update global column index array
				int rowGlobal = cvora[mpi_rank]*m + rowLocal; // global row index

				int index_f = A.get_ncols_eachRow_i(rowLocal);
				int index_l = A.get_ncols_eachRow_i(rowLocal+1)-1;

				int arrIndex = 0; // array index for cindGlobalArr_eachRow or valuesArr_eachRow
				for (int i=index_f; i<=index_l; ++i) {
	#if PETSC_DEBUG_LEVEL > 0
					assert(arrIndex<ncols_eachRow);
					assert(rowGlobal==A.get_global_rind(i, mpi_rank, m, cvora));
	#endif
					cindGlobalArr_eachRow[arrIndex] = A.get_global_cind(i, mpi_rank, m, cvora, cv_gl, NcontrolParams, ncv_gg);

					// increase arrIndex
					++arrIndex;
				}

				// Update values array
				arrIndex = 0;
				for (int i=index_f; i<=index_l; ++i) {
					valuesArr_eachRow[arrIndex] = A.get_values(i);
					++arrIndex;
				}

				// call MatSetValues
				MatSetValues(A_, 1, &rowGlobal, ncols_eachRow, cindGlobalArr_eachRow, valuesArr_eachRow, INSERT_VALUES);
				/*
				 * MatSetValues(Mat A,int m,const int idxm[],int n,const int idxn[],const PetscScalar values[], INSERT VALUES or ADD VALUES);
				 *   MatSetValues() uses the standard C convention, where the row and column matrix indices begin with zero.
				 *   The array values is logically two-dimensional and is given in row major order, meaning that the value to be put in row idxm[i] and column idxn[j] is located in values[i*n+j].
				 */

				delete [] cindGlobalArr_eachRow; 	cindGlobalArr_eachRow = NULL;
				delete [] valuesArr_eachRow; 		valuesArr_eachRow = NULL;
			} else {
				++myEmptyRows;
			}
		}

		// 2. vecC and d part
		{
			int sizeOfA = cvora[mpi_size]*m; // number of rows in A

			// vecC: vecC is distributed among CPUs (CAUTION: this can be problematic!!)
			// send vecC
			for(int iParam=0; iParam<NcontrolParams; ++iParam) {
				if(mpi_rank < mpi_size-1) {
					MPI_Request request;
					MPI_Status status;
					MPI_Isend(vecC[iParam], nno*m, MPI_DOUBLE, mpi_size-1, 320+iParam*mpi_size + mpi_rank, mpi_comm, &request);
					MPI_Wait( &request, &status );
				} else {
					int *cind = new int [sizeOfA+NcontrolParams];
					double *values = new double [sizeOfA+NcontrolParams];

					// receive vecC
					for(int icpu=0; icpu<mpi_size-1; ++icpu) {
						int startIndex = cvora[icpu]*m;
						int count = (cvora[icpu+1]-cvora[icpu])*m;
						MPI_Request request;
						MPI_Status status;
						MPI_Irecv(values+startIndex, count, MPI_DOUBLE, icpu, 320+iParam*mpi_size + icpu, mpi_comm, &request);
						MPI_Wait( &request, &status );
					}

					// fill up the rest vecC
					for(int i=0; i<nno*m; ++i)
						values[cvora[mpi_size-1]*m+i] = vecC[iParam][i];

					// d
					for(int iAdd=0; iAdd<NcontrolParams; ++iAdd) {
						if(iAdd==iParam)
							values[sizeOfA+iAdd] = d[iParam];
						else
							values[sizeOfA+iAdd] = 0.0;
					}

					// make cind
					for(int i=0; i<sizeOfA+NcontrolParams; ++i)
						cind[i] = i;

					MatSetValues(A_, 1, &sizeOfA, sizeOfA+1, cind, values, INSERT_VALUES);

					delete [] cind;
					delete [] values;
				}
			}
		}

		MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);

		// Show matrix A_ on screen !
		if(ShowPetscMatrixMatlab) { // ShowPetscMatrixMatlab is defined as a static variable at the beginning of this file and set as "false" by default
			PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);
			MatView(A_,PETSC_VIEWER_STDOUT_WORLD);
		}

		MPI_Allreduce(&myEmptyRows, &emptyRows, 1, MPI_INT, MPI_SUM, mpi_comm);
		if(emptyRows>0 && mpi_rank == 0)
			cout<<"CAUTION in PetscSolvers2::setLinSysForPetsc():  -- total "<<emptyRows<<" rows are empty"<<endl;
		MPI_Barrier(mpi_comm);
	}

	/*
	 * Method: writePetscLinSysOnFile
	 * ----------------------------
	 * Write the linear system stored in the Petsc object on a file (binary or ascII) in a Matlab-readable format
	 */
	 void writePetscLinSysOnFile (const char filename[], bool binary) {
		 PetscViewerCreate(mpi_comm, &viewer);
		 if(binary) {
			 PetscViewerSetType(viewer, PETSCVIEWERBINARY); 
			 			// possible methods (Petsc ver.3.2.7): PETSCVIEWERSOCKET, PETSCVIEWERASCII, PETSCVIEWERBINARY, 
			 			//                                     PETSCVIEWERSTRING, PETSCVIEWERDRAW, PETSCVIEWERVU, 
						//                                     PETSCVIEWERMATHEMATICA, PETSCVIEWERNETCDF, PETSCVIEWERHDF5, 
			 			//                                     PETSCVIEWERMATLAB, PETSCVIEWERAMS 

			 PetscViewerBinaryOpen(mpi_comm, filename, FILE_MODE_WRITE, &viewer);	
		 } else {
			 PetscViewerSetType(viewer, PETSCVIEWERASCII); 
		 	 PetscViewerASCIIOpen(mpi_comm, filename, &viewer);
		 }

		 MatView(A_,viewer);

		 PetscViewerDestroy(&viewer);
	 }

	/*
	 * Method: solveGMRES
	 * ------------------
	 * interface to call solveGMRES_solver()
	 */
	 template <class MatT>
	 void solveGMRES(MatT &A, double *phi, const double *rhs, double *phiInit, const int *cvora, const int *cv_gl, const int m,
			 const int NcontrolParams = 0, const int ncv_gg = 0, double **vecC = NULL, const double *d = NULL, const int step=-1, const int newtonIter=-1) {
		 KspTols.psalcStep = step;
		 KspTols.newtonIter = newtonIter;

		 int nIter;
		 double absResid;
		 solveGMRES_solver(nIter, absResid, A, phi, rhs, phiInit, cvora, cv_gl, m, NcontrolParams, ncv_gg, vecC, d);
	 }

	 template <class MatT>
	 void solveGMRES(int &nIter, double &absResid,
			 MatT &A, double *phi, const double *rhs, double *phiInit, const int *cvora, const int *cv_gl, const int m,
			 const int NcontrolParams = 0, const int ncv_gg = 0, double **vecC = NULL, const double *d = NULL, const int step=-1, const int newtonIter=-1) {
		 KspTols.psalcStep = step;
		 KspTols.newtonIter = newtonIter;

		 solveGMRES_solver(nIter, absResid, A, phi, rhs, phiInit, cvora, cv_gl, m, NcontrolParams, ncv_gg, vecC, d);
	 }

	 template <class MatT>
	 void solveGMRES(vector<pair<int, double> >& kspMonitorHistory,
			 int &nIter, double &absResid, MatT &A, double *phi, const double *rhs, double *phiInit, const int *cvora, const int *cv_gl, const int m,
			 const int NcontrolParams = 0, const int ncv_gg = 0, double **vecC = NULL, const double *d = NULL, const int step=-1, const int newtonIter=-1) {
		 solveGMRES(nIter, absResid, A, phi, rhs, phiInit, cvora, cv_gl, m, NcontrolParams, ncv_gg, vecC, d, step, newtonIter);
		 kspMonitorHistory = KspTols.kspMonitorHistory;
	 }

	 /*
	  * Method: matMultiplyVec
	  * ----------------------
	  * Solve for A*b = x (simple multiplication between a matrix A and a vector b)
	  */
	 template <class MatT>
	 void matMultiplyVec(MatT &A, double *xVec, const double *bVec, const int *cvora, const int *cv_gl, const int m,
			 const int NcontrolParams = 0, const int ncv_gg = 0, double **vecC = NULL, const double *d = NULL) {
		 bool oldHasInitialGuess = hasInitialGuess;
		 hasInitialGuess = false;

		 setLinSysForPetsc(A, bVec, NULL, cvora, cv_gl, m, NcontrolParams, ncv_gg, vecC, d);

		 hasInitialGuess = oldHasInitialGuess;

		 MatMult(A_, b_, x_);
		 // VecView(x_,PETSC_VIEWER_STDOUT_WORLD);

		 VecAssemblyBegin(x_); // ? Do we need this?
		 VecAssemblyEnd(x_);

		 int nno = cvora[mpi_rank+1] - cvora[mpi_rank];
		 for (int i = 0; i < nno*m; i++) {
			 int row = cvora[mpi_rank]*m + i;
			 VecGetValues(x_, 1, &row, &xVec[i]);   // inefficient, check how to change in future
		 }
		 if(NcontrolParams>0 && mpi_rank == mpi_size-1) {
			 for(int iParam=0; iParam<NcontrolParams; ++iParam) {
				 int row = cvora[mpi_size]*m + iParam;
				 VecGetValues(x_, 1, &row, &xVec[nno*m+iParam]);
			 }
		 }
	 }

	 /*
	  * Method: calcTrueResidualNorm
	  * ----------------------------
	  * Calculate the norm of the true residual: || A*x - b ||_2
	  */
	 double calcTrueResidualNorm() {
		 Vec Ax_;
		 VecDuplicate(b_, &Ax_);
		 VecSet(Ax_, 0.0);

		 // Calculate A*x
		 MatMult(A_, x_, Ax_);

		 // Calculate A*x-b
		 VecAXPY(Ax_, -1.0, b_); // VecAXPY(Vec y, PetscScalar a, Vec x): y = y + a*x

		 // Calculate the 2-norm of A*x-b
		 double trueNorm;
		 VecNorm(Ax_, NORM_2, &trueNorm);

		 // Clear and Return
		 VecDestroy(&Ax_);

		 return trueNorm;
	 }

protected:
	/*
	 * Method: solveGMRES_solver
	 * -------------------------
	 *
	 */
	template <class MatT>
	void solveGMRES_solver(int &nIter, double &absResid, MatT &A, double *phi, const double *rhs, double *phiInit, const int *cvora, const int *cv_gl, const int m,
			const int NcontrolParams, const int ncv_gg, double **vecC, const double *d) {
		setLinSysForPetsc(A, rhs, phiInit, cvora, cv_gl, m, NcontrolParams, ncv_gg, vecC, d);
//MatView(A_,PETSC_VIEWER_STDOUT_WORLD); // Note: ASCII matrix output not allowed for matrices with more than 1024 rows
//if(mpi_rank==0)
//	cout<<endl
//	    <<"Prinitng b_:"<<endl;
//VecView(b_,PETSC_VIEWER_STDOUT_WORLD);
//if(mpi_rank==0)
//	cout<<endl
//	    <<"Prinitng x_:"<<endl;
//VecView(x_, PETSC_VIEWER_STDOUT_WORLD);

//		KSPSetUp(ksp); // KSPSetUP() sets up the internal data structures for the later use of an iterative solver.
//		               // However, you don't have to call it since it will be called in KSPSolve() anyways.

		//KSPSetOperators(ksp, A_, A_, DIFFERENT_NONZERO_PATTERN);
		KSPSolve(ksp, b_, x_);
		// VecView(x_,PETSC_VIEWER_STDOUT_WORLD);

		KSPGetIterationNumber(ksp, &nIter);
		KSPGetResidualNorm(ksp, &absResid); // Gets the last (approximate preconditioned) residual norm that has been computed.

		int maxiter;
		double rel, abs, div;
		KSPGetTolerances(ksp, &rel, &abs, &div, &maxiter);
		lout(INFO_LO) << "tresholds: rel/abs/div/iter:\t" << rel << " " << abs << " " << div << " " << maxiter;
		lout(INFO_LO) << "\tPETSC: iter/absResidual:\t" << nIter << "\t" << absResid << endl;

		//VecView(x_,PETSC_VIEWER_STDOUT_WORLD);

		VecAssemblyBegin(x_);
		VecAssemblyEnd(x_);
		int nno = cvora[mpi_rank+1] - cvora[mpi_rank];
		for (int i = 0; i < nno*m; i++) {
			int row = cvora[mpi_rank]*m + i;
			VecGetValues(x_, 1, &row, &phi[i]);   // inefficient, check how to change in future
		}
		if(NcontrolParams>0 && mpi_rank == mpi_size-1) {
			for(int iParam=0; iParam<NcontrolParams; ++iParam) {
				int row = cvora[mpi_size]*m + iParam;
				VecGetValues(x_, 1, &row, &phi[nno*m+iParam]);
			}
		}

		// Get the extreme singular values
#ifdef CALC_SINGULAR_VALUES
		if(PETSC_DEBUG_LEVEL <= 1) { // If PETSC_DEBUG_LEVEL>1, KSPMonitorSingularValue() is called instead.
			                         // KSPMonitorSingularValue() prints the two norm of the true residual and estimation of the extreme singular values of the preconditioned problem at each iteration.
			double sigmaMax, sigmaMin;
			KSPComputeExtremeSingularValues(ksp, &sigmaMax, &sigmaMin);
					// KSPComputeExtremeSingularValues() computes the extreme singular values for the preconditioned operator. Called after or during KSPSolve().
					// Note: 1. One must call KSPSetComputeSingularValues() before calling KSPSetUp()
					//       2. Estimates of the smallest singular value may be very inaccurate, especially if the Krylov method has not converged.
					//          The largest singular value is usually accurate to within a few percent if the method has converged.
					//       3. You may want to disable restarts if using KSPGMRES, otherwise this estimate will only be using those iterations after the last restart.
			if(mpi_rank==0)
				cout<<"           >> PetscSolver2::solveGMRES_solver(): Maximum singular value="<<sigmaMax<<", Minimum singular value="<<sigmaMin<<endl;
		}
#endif
	}

public:
	/*
	 * Method: setPCreuse
	 * ------------------
	 *
	 */
	void setPCreuse(const bool pcReuse) {
		/* set KSP options */
		if(pcReuse) {
			KSPSetOperators(ksp, A_, A_, SAME_PRECONDITIONER);
#if PETSC_DEBUG_LEVEL > 0
			if(mpi_rank==0)
				cout<<"= PC operation: SAME_PRECONDITIONER"<<endl;
#endif
		} else {
			KSPSetOperators(ksp, A_, A_, DIFFERENT_NONZERO_PATTERN);
		                                  /* Note: KSPSetOperators(KSP ksp, Mat Amat, Mat Pmat, MatStructure flag)
		                                   *  	     ksp -the KSP context
		                                   *         Amat - the matrix associated with the linear system
		                                   *         Pmat - the matrix to be used in constructing the preconditioner, usually the same as Amat.
		                                   *         flag - flag indicating information about the preconditioner matrix structure during successive linear solves.
		                                   *                This flag is ignored the first time a linear system is solved, and thus is irrelevant when solving just one linear system.
		                                   *                >> POSSIBLE OPTIONS
		                                   *                   SAME_PRECONDITIONER - Pmat is identical during successive linear solves.
		                                   *                   SAME_NONZERO_PATTERN - Pmat has the same nonzero structure during successive linear solves.
		                                   *                   DIFFERENT_NONZERO_PATTERN - Pmat does not have the same nonzero structure.
		                                   */
#if PETSC_DEBUG_LEVEL > 0
			if(mpi_rank==0)
				cout<<"= PC operation: DIFFERENT_NONZERO_PATTERN"<<endl;
#endif
		}
	}

protected:
	void initPetscSolver(const int *cvora, const int *nbocv2_i, const vector<int> &nbocv2_v, const int m, const string pcType, const bool pcReuse, const bool hasInitialGuess, const int NcontrolParams,
			const int pcLevels=0) {
		int argc = 0; char **argv;
		static char help[] = "nothing!";
		PetscInitialize(&argc, &argv, (char *)0, help);

		this->hasInitialGuess = hasInitialGuess;

		/* Create linear solver context (a Krylov sub-space methods) */
		KSPCreate(mpi_comm, &ksp);

		/* allcoate mem for the vectors and the matrix */
		initMatVec(cvora, nbocv2_i, nbocv2_v, m, NcontrolParams);

		/* set preconditioner */
		// step 1. extracting the preconditioner context from the KSP context
		KSPGetPC(ksp, &pc);
		// step 1-2. set a flag so that the extreme singular values will be calculated via a Arnoldi process as the linear system is solved by GMRes
#ifdef CALC_SINGULAR_VALUES
		KSPSetComputeSingularValues(ksp, PETSC_TRUE); // Note: Only CG (Lanczos process) or GMRes(Arnoldi) supports this option
#endif
		// step 2. set the method
		PCType pcMethod;
		if(pcType=="BJACOBI")
			pcMethod = PCBJACOBI; // Block Jacobi
		else if(pcType=="ASM")
			pcMethod = PCASM; // Additive Schwarz
		else if(pcType=="ILU") {
			if(mpi_size==1)
				pcMethod = PCILU;   // Incomplete LU -- The ILU preconditioner in Petsc only works for serial -- For parallel, call an external package
			else
				pcMethod = PCHYPRE; // Incomplete LU -- For parallel, use PCHYPRE with the PCHYPRESetType(pc, "euclid") command
		} else if(pcType=="AMG")  // Algebraic Multi-grid -- Use PCHYPRE with the PCHYPRESetType(pc, "boomeramg") command
			                      // Note that the interfaces for this PC is still in development (Feb, 2014)
			pcMethod = PCHYPRE;
		else {
			if(mpi_rank==0)
				cout<<"WARNING in PetscSolver2::initPetscSolver(): No Pre-conditioner will be used"<<endl;
			pcMethod = PCNONE;
		}

#if PETSC_DEBUG_LEVEL > 0
		if(mpi_rank==0) {
			cout<<"============================================="<<endl;
			cout<<"= PC Type: "<<pcMethod<<endl;
		}
#endif
		PCSetType(pc, pcMethod);

		// step 3. set options for HYPRE-ILU or HYPRE-AMG if HYPRE is required
		if(pcMethod == PCHYPRE) {
			// Note: hypre error codes
			//          HYPRE_ERROR_GENERIC         1    generic error
			//          HYPRE_ERROR_MEMORY          2    unable to allocate memory
			//          HYPRE_ERROR_ARG             4    argument error
			//                                           bits 4-8 are reserved for the index of the argument error
			//          HYPRE_ERROR_CONV          256    method did not converge as expected

			if(pcType=="ILU") {
//				PCFactorSetShiftType(pc, MAT_SHIFT_POSITIVE_DEFINITE); // avoid zero-pivot
				PCHYPRESetType(pc, "euclid"); // set hypre : pilut, parasails, boomeramg, euclid
											  //             -- "euclid" is the ILU
											  //             -- "boomeramg" (boomerAMG) is the AMG
				                              //             -- "pilut" is Y. Saads's dual-threshold incomplete factorization algorithm (The algorithm produces an approximate factorization LU)
				                              //             -- "parasails" is a sparse approximate inverse (SPAI) preconditioner (using a priori sparsity patterns and least-squares (Frobenius norm) minimization)

				stringstream ss;
				ss<<pcLevels;
				string buffer = ss.str();
				PetscOptionsSetValue("-pc_hypre_euclid_levels", buffer.c_str()); // Number of levels of fill ILU(k)
//				PetscOptionsSetValue("-pc_hypre_euclid_bj", "0");        // Don't use block Jacobi ILU(k) instead of parallel ILU
				PetscOptionsSetValue("-pc_hypre_euclid_bj", PETSC_NULL); // Use block Jacobi ILU(k) instead of parallel ILU
#if PETSC_DEBUG_LEVEL > 1
				PetscOptionsSetValue("-pc_hypre_euclid_print_statistics", 2);
#else
				PetscOptionsSetValue("-pc_hypre_euclid_print_statistics", "NO");
#endif

#if PETSC_DEBUG_LEVEL > 0
				showHYPREoptions_Euclid();
#endif
			} else if(pcType=="AMG") { // Algebraic multigrid - using an external package "hypre"
				                       // Note: Still testing this PC -- Currently gives "Inf/NaN" errors
				PCHYPRESetType(pc, "boomeramg"); // set hypre : pilut, parasails, boomeramg, euclid
											     //             -- euclid is the ILU
											     //             -- boomeramg (boomerAMG) is the AMG

//				PetscOptionsSetValue("-pc_hypre_boomeramg_cycle_type", "V"); // Cycle type: "", "V", "W"

				PetscOptionsSetValue("-pc_hypre_boomeramg_max_levels", "3"); // At least 2
				PetscOptionsSetValue("-pc_hypre_boomeramg_max_iter", "5"); // The number of iterations (V-cycles) per each boomerAMG call.
				                                                            // "-ksp_max_it" and "-ksp_rtol STILL" determine the total number of iterations and tolerance for the Krylov solver:
				                                                            // For example, if -pc_hypre_boomeramg_max_iter is 2 and -ksp_max_it is 10, then AT MOST twenty V-cycles of boomeramg will be called.
				PetscOptionsSetValue("-pc_hypre_boomeramg_rtol", "0.0"); // Convergence tolerance for each boomerAMG call.
				                                                         // If -pc_hypre_boomeramg_max_iter is set to 2,
				                                                         // -pc_hypre_boomeramg_rtol should be set to 0.0 - the default - to strictly use a fixed number of iterations per hypre call.
				PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold", "0.25"); // Important parameter that will affect both coarsening and interpolation.
				                                                                      // Threshold for being strongly connected.
				                                                                      // Default = 0.25 (a good choice for 2-dimensional problems and the low complexity coarsening algorithms)
				                                                                      // For 3-dimensional problems, a better choice appears to be 0.5, if one uses the default coarsening algorithm or CLJP
//				PetscOptionsSetValue("-pc_hypre_boomeramg_truncfactor", "0.0"); // Truncation factor for interpolation (0=no truncation)

//				PetscOptionsSetValue("-pc_hypre_boomeramg_agg_nl",        "0"); // Number of levels of aggressive coarsening
//				PetscOptionsSetValue("-pc_hypre_boomeramg_agg_num_paths", "1"); // Number of paths for aggressive coarsening
//				                                                                // At least 1

//				PetscOptionsSetValue("-pc_hypre_boomeramg_max_row_sum", "0.9"); // Maximum row sum

//				PetscOptionsSetValue("-pc_hypre_boomeramg_grid_sweeps_down",   "1"); // Number of sweeps for the down cycles
//				PetscOptionsSetValue("-pc_hypre_boomeramg_grid_sweeps_up",     "1"); // Number of sweeps for the up cycles
//				PetscOptionsSetValue("-pc_hypre_boomeramg_grid_sweeps_coarse", "1"); // Number of sweeps for the coarse cycles

				PetscOptionsSetValue("-pc_hypre_boomeramg_relax_type_all", "SOR/Jacobi"); // Relax type for the up and down cycles: You can separately gives options for up and down: -pc_hypre_boomeramg_relax_type_up  pc_hypre_boomeramg_relax_type_down
				                                                                // Options: "Jacobi", "sequential-Gauss-Seidel", "seqboundary-Gauss-Seidel", "SOR/Jacobi", "backward-SOR/Jacobi",
				                                                                //          "" ([5] hybrid chaotic Gauss-Seidel -- works only with OpenMP), "symmetric-SOR/Jacobi",
				                                                                //          "" ([7]), "l1scaled-SOR/Jacobi", "Gaussian-elimination",
				                                                                //          "" ([10]), "" ([11]), "" ([12]), "" ([13), "" ([14]),
				                                                                //          "CG" (non-stationary), "Chebyshev", "FCF-Jacobi", "l1scaled-Jacobi"
				                                                                // Default = 9
				                                                                // Note that the option -pc_hypre_boomeramg_relax_type_all defaults to symmetric relaxation
				                                                                // (symmetric-SOR/Jacobi), which is required for Krylov solvers like CG that expect symmetry.
				                                                                // Otherwise, you may want to use -pc_hypre_boomeramg_relax_type_all SOR/Jacobi.
				PetscOptionsSetValue("-pc_hypre_boomeramg_relax_type_coarse", "SOR/Jacobi"); // Relax type on coarse grid
				PetscOptionsSetValue("-pc_hypre_boomeramg_relax_weight_all",       "0"); // Relaxation weight for all levels (0 = hypre estimates, -k = determined with k CG steps)
//				PetscOptionsSetValue("-pc_hypre_boomeramg_outer_relax_weight_all", "0"); // Outer relaxation weight for all levels (-k = determined with k CG steps)
//				PetscOptionsSetValue("-pc_hypre_boomeramg_no_CF",                  "?"); // Do not use CF-relaxation

//				PetscOptionsSetValue("-pc_hypre_boomeramg_measure_type", "local"); // "local", "global"

				PetscOptionsSetValue("-pc_hypre_boomeramg_coarsen_type", "modifiedRuge-Stueben"); // Coarsening method
				                                                                  // Options: "CLJP", "Ruge-Stueben", "", "modifiedRuge-Stueben", "", "", "Falgout", "", "PMIS", "", "HMIS"
				                                                                  //   CLJP   : (Cleary-Luby-Jones-Plassman) Entirely parallel algorithm. Large complexities particularly on the interior of CPU domains (requires more memory)
				                                                                  //   RS     : Has some issues at CPU boundaries (not entirely parallel).
				                                                                  //   Falgout: Use RS on the interior and use CLJP near CPU boundaries
				                                                                  //   BC-RS  : First treat CPU boundaries using CLJP and then fill out the interior using RS (not efficient)
                                                                                  // Default: Falgout (a combination of CLJP and the classical RS)
				                                                                  // For structured (and symmetric) 2D, Falgout is the best choice. For unstructured 3D, CLJP is also a good choice.

				PetscOptionsSetValue("-pc_hypre_boomeramg_interp_type", "classical"); // Interpolation type: "classical", "", "", "direct", "multipass", "multipass-wts", "ext+i",
				                                                                      //                     "ext+i-cc", "standard", "standard-wts", "", "", "FF", "FF1"};


//				PetscOptionsSetValue("-pc_hypre_boomeramg_print_statistics", "2"); // Print statistics
//				PetscOptionsSetValue("-pc_hypre_boomeramg_print_debug", "2"); // Print debug information

				// Note: If you wish to use BoomerAMG WITHOUT a Krylov method use -ksp_type richardson NOT -ksp_type preonly
				//       and use -ksp_max_it to control the number of V-cycles.

				showHYPREoptions_bommerAMG();
			}
		}

		// step 4. set reuse
		setPCreuse(pcReuse); // Note: There are three possible options: SAME_PRECONDITIONER, SAME_NONZERO_PATTERN, and DIFFERENT_NONZERO_PATTERN
		                     //       If "pcResue" is true, SAME_PRECONDITIONER. Otherwise, DIFFERENT_NONZERO_PATTERN

		// step 5. set monitor
		if(KspTols.convergMonitorInterv>0)
			KSPMonitorSet(ksp, MyKSPMonitor, &KspTols, 0);

#ifdef CALC_SINGULAR_VALUES
		if(PETSC_DEBUG_LEVEL > 1) {
			KSPMonitorSet(ksp, KSPMonitorSingularValue, NULL, NULL); // KSPMonitorSet() sets an ADDITIONAL function to be called at every iteration to monitor the residual/error etc.
			                                                         // KSPMonitorSingularValue() prints the two norm of the true residual and estimation of the extreme singular values
			                                                         // of the preconditioned problem at each iteration.
		}
#endif

		MPI_Barrier(mpi_comm);

		/* set initial-guess options and GMRES options */
		if(hasInitialGuess)
			KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
		else
			KSPSetInitialGuessKnoll(ksp, PETSC_TRUE); // Tells the iterative solver to use PCApply(pc,b,..) to compute the initial guess (The Knoll trick)

		KSPSetType(ksp, KSPGMRES);

		setGmresRestart(100);

		KSPGMRESSetOrthogonalization(ksp, KSPGMRESClassicalGramSchmidtOrthogonalization); // use classical (unmodified) Gram-Schmidt to orthogonalize against the Krylov space (fast) (the default)
//		KSPGMRESSetOrthogonalization(ksp, KSPGMRESModifiedGramSchmidtOrthogonalization);  // use modified Gram-Schmidt in the orthogonalization (more stable, but slower)
		KSPSetFromOptions(ksp); // set runtime options e.g. -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol> : These options will override those specified above

#if PETSC_DEBUG_LEVEL > 0
		if(mpi_rank==0)
			cout<<"============================================="<<endl<<endl;
#endif

#if PETSC_DEBUG_LEVEL > 1
		PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);
		KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);
#endif
	}

	/*
	 * Method: initMatVec
	 * ------------------
	 * Upload the matrix and the vectors to the memory. This method will be called by initPetscSolver
	 */
	void initMatVec(const int *cvora, const int *nbocv2_i, const vector<int> &nbocv2_v, const int m, const int NcontrolParams) {
		int nno = cvora[mpi_rank+1] - cvora[mpi_rank]; // number of CV in the current CPU

		// allcoate mem for vec (b and x)
		if(mpi_rank == mpi_size-1)
			VecCreateMPI(mpi_comm, nno*m+NcontrolParams, cvora[mpi_size]*m+NcontrolParams, &b_);
		else
			VecCreateMPI(mpi_comm, nno*m,                cvora[mpi_size]*m+NcontrolParams, &b_);
		VecSetFromOptions(b_);
		VecDuplicate(b_, &x_);
		VecSet(b_, 0.0);
		VecSet(x_, 0.0);

		// allocate memory for the Matrix
		/*
		 * Note: For exapmle, A = [ S1  O_12 O_13
		 *                         O_21  S2  O_23
		 *                         O_31 O_32  S3 ]
		 *                    where S is square sub-matrix on the diagonal, O's are off-diagonal sub-matrices,
		 *                    "Square matrices" are {S1, S2, S3} and "Rectangular matrices" are {[O_12 O_13], [O_21 O_23], [O_31 O_32]}.
		 *                    ( each CPU has one diagonal block, i.e. the number of rows of S1 is equal to ncv*m of CPU1 )
		 *       The user need not preallocate exactly the correct amount of space; as long as a sufficiently close estimate
		 *     is given, the high efficiency for matrix assembly will remain
		 *
		 *       For more details, see the PetSc user manual.
		 */
		int *diag_sizes; // number of non-zeros in the each row of the square matrix
		int *offd_sizes; // number of non-zeros in the each row of the rectangular matrix
		if(mpi_rank == mpi_size-1) { // The number of rows in the last block is greater than the others
			diag_sizes = new int[nno*m+NcontrolParams];
			offd_sizes = new int[nno*m+NcontrolParams];
		} else {
			diag_sizes = new int[nno*m];
			offd_sizes = new int[nno*m];
		}

		// Give the block-wise non-zero information to PetSc
		for (int ino=0; ino<nno; ino++) { // looping over all the CVs in each CPU
			// First, count the (maximum possible) number of non-zeros in the first row of each block
			if(mpi_rank == mpi_size-1) {
				diag_sizes[ino*m] = NcontrolParams; // The last column belongs to the diagonal block
				offd_sizes[ino*m] = 0;
			} else {
				diag_sizes[ino*m] = 0;
				offd_sizes[ino*m] = NcontrolParams; // The last column belongs to the off-diagonal block
			}

			int non_f = nbocv2_i[ino];
			int non_l = nbocv2_i[ino+1]-1;

			for (int non=non_f; non<=non_l; non++) {
				int ino_nbr = cvora[mpi_rank] + nbocv2_v[non];

				if ((ino_nbr >= cvora[mpi_rank]) && (ino_nbr < cvora[mpi_rank+1]))
					diag_sizes[ino*m] += m;
				else
					offd_sizes[ino*m] += m;
			}

			// Second, copy the information of the first row to the other rows
			for (int i=1; i<m; i++) {
				diag_sizes[ino*m+i] = diag_sizes[ino*m];
				offd_sizes[ino*m+i] = offd_sizes[ino*m];
			}
		}
		if(mpi_rank == mpi_size-1) { // The last block has one more row
			for(int iParam=0; iParam<NcontrolParams; ++iParam) {
				diag_sizes[nno*m+iParam] = nno*m + NcontrolParams;
				offd_sizes[nno*m+iParam] = (cvora[mpi_size] - nno)*m;
			}
		}

		if(mpi_rank == mpi_size-1) {
			MatCreateMPIAIJ(mpi_comm, nno*m+NcontrolParams, nno*m+NcontrolParams, cvora[mpi_size]*m+NcontrolParams, cvora[mpi_size]*m+NcontrolParams, 0, &diag_sizes[0], 0, &offd_sizes[0], &A_);
			/*
			 * MatCreateMPIAIJ(MPI Comm comm,int m,int n,int M,int N,int d_nz,int *d_nnz,int o_nz,int *o_nnz,Mat *A);
			 * 	m: the number of local rows (determined by PETSc if m is PETSC_DECIDE.
			 *                               If PETSC DECIDE is not used for the arguments m and n,
			 *                               then the user must ensure that they are chosen to be compatible with the vectors)
			 *     This value should be the same as the local size used in creating the y vector for the matrix-vector product y = Ax
			 *  n: the number of local columns (determined by PETSc if n is PETSC_DECIDE.
			 *                                  If PETSC DECIDE is not used for the arguments m and n,
			 *                                  then the user must ensure that they are chosen to be compatible with the vectors)
			 *     This value should be the same as the local size used in creating the x vector for the matrix-vector product y = Ax
			 *  M: the number of global rows
			 *  N: the number of global columns
			 *  d_nz: number of non-zeros per row in the DIAGONAL sub-matrix (same value is used for all local rows)
			 *  d_nnz: array containing the number of non-zeros in the various rows of the DIAGONAL sub-matrix (possibly different for each row)
			 *  o_nz: number of non-zeros per row in the OFF-DIAGONAL sub-matrix (same value is used for all local rows)
			 *  o_nnz: array containing the number of non-zeros in the various rows of the OFF-DIAGONAL sub-matrix (possibly different for each row)
			 */
		} else {
			MatCreateMPIAIJ(mpi_comm, nno*m, nno*m, cvora[mpi_size]*m+NcontrolParams, cvora[mpi_size]*m+NcontrolParams, 0, &diag_sizes[0], 0, &offd_sizes[0], &A_);
		}

		delete [] diag_sizes; 	diag_sizes = NULL;
		delete [] offd_sizes; 	offd_sizes = NULL;
	}

	/*
	 * Method: finalizePetscSolver
	 * ---------------------------
	 * finalize the PetscSolver object
	 */
	PetscErrorCode finalizePetscSolver() {
		PetscErrorCode ierr;

		ierr = VecDestroy(&x_);  CHKERRQ(ierr);
		ierr = VecDestroy(&b_);  CHKERRQ(ierr);
		ierr = MatDestroy(&A_);  CHKERRQ(ierr);
		ierr = KSPDestroy(&ksp); CHKERRQ(ierr); // note: you should NOT destroy PC if you call KSPDestroy

		PetscFinalize();
				/* This routine
				 * - finalizes the PETSc libraries as well as MPI
				 * - provides summary and diagnostic information if certain runtime options are chosen (e.g., -log_summary).
				 */

		return 0;
	}

	/*
	 * Method: showHYPREoptions_Euclid()
	 * ---------------------------------
	 * Print some HYPRE PC options on the screen : level, Block-Jacobi, statistics
	 */
	PetscErrorCode showHYPREoptions_Euclid() {
		int levels;
		PetscBool bjilu, printStats;
		PetscBool flag;

		if(mpi_rank==0)
			printf("= HYPRE PCILU-Euclid Options: \n");

		PetscOptionsGetInt(PETSC_NULL, "-pc_hypre_euclid_levels", &levels, &flag);
		if(mpi_rank==0)
			if (flag)
				printf("=   Number of levels = %d \n", levels);
			else
				printf("=   Cannot find Number of levels \n");

		PetscOptionsBool("-pc_hypre_euclid_bj", PETSC_NULL, 0, PETSC_TRUE, &bjilu, 0);
		if(mpi_rank==0)
			if (bjilu)
				printf("=   Use block-Jacobi ILU(k) instead of parallel ILU = TRUE \n");
			else
				printf("=   Use block-Jacobi ILU(k) instead of parallel ILU = FALSE \n");

		PetscOptionsBool("-pc_hypre_euclid_print_statistics", PETSC_NULL, 0, PETSC_FALSE, &printStats, 0);
		if(mpi_rank==0)
			if (printStats)
				printf("=   Print statistics = TRUE \n");
			else
				printf("=   Print statistics = FALSE \n");

		MPI_Barrier(mpi_comm);
		return 0;
	}

	/*
	 * Method: showHYPREoptions_bommerAMG()
	 * ---------------------------------
	 * Print some HYPRE PC options on the screen : coarsen_type, max_iter, max_levels, strong_threshold
	 */
	PetscErrorCode showHYPREoptions_bommerAMG() {
		int hypreMaxIter, hypreMaxLevels;
		PetscReal alpha;
		char hypreCoarsenName[10];
		PetscBool flag;

		if(mpi_rank==0)
			printf("= HYPRE PCAMG-boomerAMG Options: \n");

		PetscOptionsGetString(NULL, "-pc_hypre_boomeramg_coarsen_type", hypreCoarsenName, sizeof(hypreCoarsenName), &flag);
		if(mpi_rank==0)
			if (flag)
				printf("=   coarsening method = %s \n", hypreCoarsenName);
			else
				printf("=   coarsening method: Cannot find coarsening method from Petsc options \n");
		PetscOptionsGetInt(PETSC_NULL, "-pc_hypre_boomeramg_max_iter", &hypreMaxIter, &flag);
		if(mpi_rank==0)
			if (flag)
				printf("=   max_iter = %d \n", hypreMaxIter);
			else
				printf("=   Cannot find max_iter from Petsc options \n");
		PetscOptionsGetInt(PETSC_NULL, "-pc_hypre_boomeramg_max_levels", &hypreMaxLevels, &flag);
		if(mpi_rank==0)
			if (flag)
				printf("=   max_levels = %d \n", hypreMaxLevels);
			else
				printf("=   Cannot find max_levels from Petsc options \n");
		PetscOptionsGetReal(PETSC_NULL, "-pc_hypre_boomeramg_strong_threshold", &alpha, &flag);
		if(mpi_rank==0)
			if (flag)
				printf("=   strong_threshold = %.3f \n", alpha);
			else
				printf("=   strong_threshold: Cannot find strong_threshold from Petsc options \n");

		return 0;
	}
};

#else
/*
 * Class: PetscSolver2
 * -------------------
 * Two-layer-ghost version of PetscSolver
 * Original code: PetscSolver in PetscSolver.h
 */
class PetscSolver2
{
protected:
	KSP ksp;        // linear solver context
	Vec x_, b_;     // solution, residual vector
	Mat A_;         // implicit operator matrix
	//VecScatter scatterContext;

public:
	/*
	 * Constructor
	 */
	PetscSolver2(int *cvora, int *nbocv2_i, vector<int> &nbocv2_v, int m, bool hasInitialGuess=false, int NcontrolParams=0, int monitorConvergInterval=0) {
		if(mpi_rank == 0)
			cout<<"PetscSolver2 constructor: do nothing since LIN_SOLVER_PETSC_ is not defined"<<endl;
		initPetscSolver(cvora, nbocv2_i, nbocv2_v, m, "ASM", false, hasInitialGuess, NcontrolParams);
	}
	/*
	 * Destructor
	 */
	~PetscSolver2() {
		if(mpi_rank == 0)
			cout<<"PetscSolver2 destructor: do nothing since LIN_SOLVER_PETSC_ is not defined"<<endl;
	}

	/*
	 * Method: setTresholds
	 * --------------------
	 * set tolerances for the KSP linear solver (Why the name is Tresholds (not Thresholds)? I don't know.. )
	 */
	void setTresholds(double zeroAbs, double zeroRel, int maxIter) {
		/* do nothing */
	}

	/*
	 * Method: setLinSysForPetsc
	 * -------------------------
	 * set the linear system: basically upload the matrix to the memory as the PETSC multi-core form
	 */
	void setLinSysForPetsc(MatT &A, const double *rhs, const double *xInit, const int *cvora, const int *cv_gl, const int m,
			const int NcontrolParams = 0, const int ncv_gg = 0, double **vecC = NULL, const double *d = NULL) {
		/* do nothing */
	}

	/*
	 * Method: writePetscLinSysOnFile
	 * ----------------------------
	 * Write the linear system stored in the Petsc object on a file (binary or ascII) in a Matlab-readable format
	 */
	 void writePetscLinSysOnFile (const char filename[], bool binary) {
		 /* do nothing */
	 }

	/*
	 * Method: solveGMRES
	 * ------------------
	 *
	 */
	void solveGMRES(MatT &A, double *phi, const double *rhs, double *phiInit, const int *cvora, const int *cv_gl, const int m,
	 			 const int NcontrolParams = 0, const int ncv_gg = 0, double **vecC = NULL, const double *d = NULL, const int step=-1, const int newtonIter=-1) {
		/* do nothing */
	}
	void solveGMRES(int &nIter, double &absResid,
	 			 MatT &A, double *phi, const double *rhs, double *phiInit, const int *cvora, const int *cv_gl, const int m,
	 			 const int NcontrolParams = 0, const int ncv_gg = 0, double **vecC = NULL, const double *d = NULL, const int step=-1, const int newtonIter=-1) {
		/* do nothing */
	}
	 void solveGMRES(vector<pair<int, double> >& kspMonitorHistory,
				 int &nIter, double &absResid, MatT &A, double *phi, const double *rhs, double *phiInit, const int *cvora, const int *cv_gl, const int m,
				 const int NcontrolParams = 0, const int ncv_gg = 0, double **vecC = NULL, const double *d = NULL, const int step=-1, const int newtonIter=-1) {
		/* do nothing */
	}

	/*
	 * Method: setPCreuse
	 * ------------------
	 *
	 */
	void setPCreuse(const bool pcReuse) {
		/* do nothing */
	}

protected:
	/*
	 * Method: initPetscSolver
	 * -----------------------
	 * initialize the PetscSolver object: basically, allocate memory for the object
	 */
	void initPetscSolver(const int *cvora, const int *nbocv2_i, const vector<int> &nbocv2_v, const int m, const string pcType, const bool pcReuse, const bool hasInitialGuess, const int NcontrolParams, const int pcLevels=0) {
		/* do nothing */
	}
};
#endif

#ifdef USE_SLEPC_WITH_PETSC
/*
 * Class: SlepcSolver2
 * -------------------
 *
 */
class SlepcSolver2: public PetscSolver2
{
public:
	/*
	 * Constructor
	 */
	SlepcSolver2(int *cvora, int *nbocv2_i, vector<int> &nbocv2_v, int m, bool hasInitialGuess=false, int NcontrolParams=0, int monitorConvergInterval=0) : PetscSolver2(cvora, nbocv2_i, nbocv2_v, m, hasInitialGuess, NcontrolParams, monitorConvergInterval) {
#if PETSC_DEBUG_LEVEL > 1
		if(mpi_rank==0) cout<<"SlepcSolver2()"<<endl;
#endif
		initSlepcSolver();
	}
	SlepcSolver2(int *cvora, int *nbocv2_i, vector<int> &nbocv2_v, int m, string &pcType, int pcLevels, bool pcReuse, bool hasInitialGuess=false, int NcontrolParams=0, int monitorConvergInterval=0) : PetscSolver2(cvora, nbocv2_i, nbocv2_v, m, pcType, pcLevels, pcReuse, hasInitialGuess, NcontrolParams, monitorConvergInterval){
#if PETSC_DEBUG_LEVEL > 1
		if(mpi_rank==0) cout<<"SlepcSolver2()"<<endl;
#endif
		initSlepcSolver();
	}

	/*
	 * Destructor
	 */
	~SlepcSolver2() {
#if PETSC_DEBUG_LEVEL > 1
		if(mpi_rank==0) cout<<"~SlepcSolver2()"<<endl;
#endif
		finalizeSlepcSolver();
	}

protected:
	/* SVD */
	SVD svd_; /* singular value solver context */

	/* EV */
	PetscBool eps_empty;
	EPS eps_; /* eigensolver context */
	PetscScalar kr, ki; /* eigenvalue,  k */
	Vec         xr, xi; /* eigenvector, x */

	/* Matrix information */
	PetscBool alreadyHasMatrix;

public:
	/* ==================================
	 * Basic matrix operations
	 * ================================== */

	/*
	 * Method: transpose
	 * -----------------
	 * Take the transpose of the matrix
	 */
	void transpose() {
		MatTranspose(A_, MAT_REUSE_MATRIX, &A_);
	}

	/* ==================================
	 * Singular Value Decomposition (SVD)
	 * ================================== */

	/*
	 * Method: getLargestSmallestSV
	 * ----------------------------
	 * Get the largest and smallest singular values
	 */
	bool getLargestSmallestSV(double& sigmaLargest, double&sigmaSmallest) {
		PetscInt nconv1, nconv2;

		// Create the singular value solver and set the solution method
		initSVDsolver();

		// First request a singular value from one end of the spectrum (largest)
		SVDSetWhichSingularTriplets(svd_, SVD_LARGEST);
		SVDSolve(svd_);
		SVDGetConverged(svd_, &nconv1); // Get number of converged singular values
		if (nconv1 > 0) {
			SVDGetSingularTriplet(svd_, 0, &sigmaLargest, PETSC_NULL, PETSC_NULL);
		}

		// First request a singular value from the other end of the spectrum (smallest)
		SVDSetWhichSingularTriplets(svd_, SVD_SMALLEST);
		SVDSolve(svd_);
		SVDGetConverged(svd_, &nconv2); // Get number of converged singular values
		if (nconv2 > 0) {
			SVDGetSingularTriplet(svd_, 0, &sigmaSmallest, PETSC_NULL, PETSC_NULL);
		}

		// Free work space
		SVDDestroy(&svd_);

		if(nconv1>0 && nconv2>0)
			return true;
		return false;
	}

	/* ==================================
	 * Eigen-decomposition
	 * ================================== */

	/*
	 * Method: setJacobianMatrixForSlepc
	 * ---------------------------------
	 * Set the linear system: basically upload the matrix to the memory in the form of the PETSC multi-core.
	 * Called by the solveEigenProblem() method
	 * Note: 1. Matrix structures (MatComprsed and MatComprsed) are defined in MatComprsed.h
	 */
	template <class MatT>
	void setJacobianMatrixForSlepc(MatT &A, const int *cvora, const int *cv_gl, const int nVars, const int ncv_gg) {
		int nno = (int) (A.get_nRows()/nVars); // local number of CVs
		assert( nno == cvora[mpi_rank+1] - cvora[mpi_rank] );

		// Check matrix size (Given matrix A and Petsc matrix A_)
		bool isOnlyJacMatrix = (A.get_nCols() == ncv_gg*nVars); // Check if the matrix is exactly a Jacobian matrix
		                                                        // (Note that the size of a Pseudo-arclength matrix is bigger than a Jacobian matrix)
		if(!isOnlyJacMatrix && mpi_rank==0)
			cout<<"WARNING in SlepcSolvers2::setJacobianMatrixForSlepc(): Given matrix is bigger than Jacobian"<<endl;

//		PetscInt m_, n_;
//		MatGetSize(A_, &m_, &n_);
//		assert(m_ == n_);
//		if(n_ != cvora[mpi_size]*nVars) {
//			PetscPrintf(mpi_comm, "SlepcSolvers2::setJacobianMatrixForSlepc(): Petsc matrix size : %d x %d\n", m_, n_);
//			PetscPrintf(mpi_comm, "                                            Target matrix size: %d x %d\n", cvora[mpi_size]*nVars, cvora[mpi_size]*nVars);
//			assert(n_ == cvora[mpi_size]*nVars);
//		}

		// Build matrix
		int myEmptyRows = 0;
		int emptyRows;
		int myNonJacElements = 0;
		int nonJacElements;

		for (int rowLocal=0; rowLocal<A.get_nRows(); ++rowLocal) {
			int ncols_eachRow = A.get_ncols_eachRow_i(rowLocal+1) - A.get_ncols_eachRow_i(rowLocal); // note: number of columns in each row of matrix A

			if(ncols_eachRow > 0) {
				// Arrays
				int *cindGlobalArr_eachRow = new int [ncols_eachRow]; // global column index array
				double *valuesArr_eachRow  = new double [ncols_eachRow]; // values array

				// Update global column index array
				int rowGlobal = cvora[mpi_rank]*nVars + rowLocal; // global row index

				int index_f = A.get_ncols_eachRow_i(rowLocal);
				int index_l = A.get_ncols_eachRow_i(rowLocal+1)-1;

				int arrIndex = 0; // array index for cindGlobalArr_eachRow or valuesArr_eachRow
				int ncols_eachRow_JacOnly = 0;

				for (int i=index_f; i<=index_l; ++i) {
#if PETSC_DEBUG_LEVEL > 0
					assert(arrIndex<ncols_eachRow);
					assert(rowGlobal==A.get_global_rind(i, mpi_rank, nVars, cvora));
#endif
					if(!isOnlyJacMatrix) {
						int localCind = A.get_cind(i);
						if(localCind < ncv_gg*nVars) {
							cindGlobalArr_eachRow[arrIndex] = A.get_global_cind(i, mpi_rank, nVars, cvora, cv_gl, 0, ncv_gg);
							assert(cindGlobalArr_eachRow[arrIndex] < cvora[mpi_size]*nVars);

							++ncols_eachRow_JacOnly;
						} else {
							cindGlobalArr_eachRow[arrIndex] = MATCOMPRESED_NOT_FOUND;

							++myNonJacElements;
						}
					} else
						cindGlobalArr_eachRow[arrIndex] = A.get_global_cind(i, mpi_rank, nVars, cvora, cv_gl, 0, ncv_gg);

					// increase arrIndex
					++arrIndex;
				}

				// Update values array
				arrIndex = 0;
				for (int i=index_f; i<=index_l; ++i) {
					if(!isOnlyJacMatrix) {
						int localCind = A.get_cind(i);
						if(localCind < ncv_gg*nVars)
							valuesArr_eachRow[arrIndex] = A.get_values(i);
						else
							valuesArr_eachRow[arrIndex] = 0.0;
					} else
						valuesArr_eachRow[arrIndex] = A.get_values(i);
					++arrIndex;
				}

				// call MatSetValues
				if(!isOnlyJacMatrix)
					MatSetValues(A_, 1, &rowGlobal, ncols_eachRow_JacOnly, cindGlobalArr_eachRow, valuesArr_eachRow, INSERT_VALUES);
				else
					MatSetValues(A_, 1, &rowGlobal, ncols_eachRow,         cindGlobalArr_eachRow, valuesArr_eachRow, INSERT_VALUES);
					/*
					 * MatSetValues(Mat A,int m,const int idxm[],int n,const int idxn[],const PetscScalar values[], INSERT VALUES or ADD VALUES);
					 *   MatSetValues() uses the standard C convention, where the row and column matrix indices begin with zero.
					 *   The array values is logically two-dimensional and is given in row major order, meaning that the value to be put in row idxm[i] and column idxn[j] is located in values[i*n+j].
					 */

				delete [] cindGlobalArr_eachRow; 	cindGlobalArr_eachRow = NULL;
				delete [] valuesArr_eachRow; 		valuesArr_eachRow = NULL;
			} else {
				++myEmptyRows;
			}
		}

		MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(A_,   MAT_FINAL_ASSEMBLY);

		// Show matrix A_ on screen !
		if(ShowPetscMatrixMatlab) { // ShowPetscMatrixMatlab is defined as a static variable at the beginning of this file and set as "false" by default
			PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);
			MatView(A_, PETSC_VIEWER_STDOUT_WORLD);
		}

		MPI_Allreduce(&myEmptyRows, &emptyRows, 1, MPI_INT, MPI_SUM, mpi_comm);
		if(emptyRows>0 && mpi_rank == 0)
			cout<<"CAUTION in SlepcSolvers2::setJacobianMatrixForSlepc(): total "<<emptyRows<<" rows are empty"<<endl;

		MPI_Allreduce(&myNonJacElements, &nonJacElements, 1, MPI_INT, MPI_SUM, mpi_comm);
		if(nonJacElements>0)
			if(nonJacElements != cvora[mpi_size]*nVars && mpi_rank == 0)
				cout<<"CAUTION in SlepcSolvers2::setJacobianMatrixForSlepc(): Total "<<nonJacElements<<" elements from the given matrix will be excluded"<<endl;

		MPI_Barrier(mpi_comm);

		// Set the length of eigenvectors so that they becomes compatible with the matrix, i.e. with the same parallel layout
		MatGetVecs(A_, PETSC_NULL, &xr);
		MatGetVecs(A_, PETSC_NULL, &xi);

		// Update alreadyHasMatrix
		alreadyHasMatrix = PETSC_TRUE;
	}

	/*
	 * Method: solveEigenProblemSlepc
	 * ------------------------------
	 * Solve the eigen-problem to get (complex) eigenvalues and eigenvectors
	 *
	 * Return:
	 *   By reference:
	 *     evalsReal = real part of the eigenvalues
	 *     evalsImag = imag part of the eigenvalues
	 *     evecsReal = real part of the eigenvectors
	 *     evecsImag = imag part of the eigenvectors
	 *     relError  = relative error to the eigenvalue or to the matrix norms for each eigenpair
	 *     numIter   = total number of iterations for the iterative eigen solver
	 *   By value:
	 *     number of converged eigen-pairs
	 *
	 * Arguements:
	 *   A      = matrix
	 *   cvora  =
	 *   cv_gl  =
	 *   nScal  =
	 *   ncv_gg =
	 *   nev = the number of eigenvalues to compute
	 *   ncv = the number of column vectors to be used by the solution algorithm
	 *         It is recommended (depending on the method) to use ncv >= 2*nev or more.
	 *   mpd = maximum projected dimension (for a more advanced usage -- Chaprter 2.6.4 in the SLEPc manual)
	 *         In the case of a large number of nev, the computation costs can be reduced by setting mpd << nev.
	 *   tol      = stopping criterion for the iterative eigen-solver, 2-norm (Default value = 1.0e-7)
	 *   max_iter = maximum number of iteration for the iterative eigen-solver
	 */
	template <class MatT>
	int solveEigenProblemSlepc(double *evalsReal, double *evalsImag, double **evecsReal, double **evecsImag, double *relError, int &numIter,
			MatT A, const int *cvora, const int *cv_gl, const int nScal, const int ncv_gg,
			const PetscInt nev, const PetscInt ncv, const PetscInt mpd, const double tol, const int max_iter) {
		// Set-up the matrix
		if(!alreadyHasMatrix)
			setJacobianMatrixForSlepc(A, cvora, cv_gl, 5+nScal, ncv_gg);
		else
			if(mpi_rank==0) cout<<"WARNING in SlepcSolver2::solveEigenProblemSlepc(): Matrix already exists. Use the pre-existing one"<<endl;

		// Initialize the eigen-problem context
		setEigenSolver(nev, ncv, mpd, tol, max_iter);

		// Solve the eigen-problem and get eigenvalues and eigenvectors
		return solveAndGetEivenpairs(evalsReal, evalsImag, evecsReal, evecsImag, relError, numIter, cvora, 5+nScal, nev);
	}

protected:
	/* ==================================
	 * Singular Value Decomposition (SVD)
	 * ================================== */

	/*
	 * Method: initSVDsolver
	 * ---------------------
	 * Create the singular value solver and set the solution method
	 */
	void initSVDsolver() {
		// Create singular value context
		SVDCreate(mpi_comm, &svd_);

		// Set operator
		SVDSetOperator(svd_, A_);

		// Set solver parameters at runtime
		SVDSetFromOptions(svd_);
		SVDSetDimensions(svd_, 1, PETSC_IGNORE, PETSC_IGNORE);
	}

	/* ==================================
	 * Eigen-decomposition
	 * ================================== */

	/*
	 * Method: initEigenSolver
	 * -----------------------
	 * Create the eigen-decomposition solver
	 */
	void initEigenSolver() {
		// Create the eigenpair solver context
		EPSCreate(mpi_comm, &eps_);
		eps_empty = PETSC_FALSE;

		// Set operator
		EPSSetOperators(eps_, A_, PETSC_NULL); // Note: In the case of a generalized problem (A*x = lambda*B*x), it would be necessary
		                                       //       also to provide matrix B as the third argument
		// Set problem type
		EPSSetProblemType(eps_, EPS_NHEP); // A standard problem (A*x = lambda*x) is Hermitian if matrix A is Hermitian.
		                                   // A generalized problem (A*x = lambda*B*x) is Hermitian if both matrix A and B are Hermitian
		                                   // and positive (semi-)definite ( If B is not Hermitian, the problem cannot be considered Hermitian
		                                   // but symmetry can still be exploited to some extent in some solvers -- problem type EPS_GHIEP ).
		                                   // By default, slepc assumes that the problem is non-Hermitian, but symmetric solvers may be more accurate.
		        // Problem types: EPS_HEP     Hermitian
		        //                EPS_NHEP    Non-Hermitian
		        //                EPS_GHEP    Generalized Hermitian
		        //                EPS_GHIEP   Generalized Hermitian indefinite
		        //                EPS_GNHEP   Generalized Non-Hermitian
		        //                EPS_PGNHEP  GNHEP with positive (semi-)definite B
//		EPSIsGeneralized(EPS eps,PetscBool *gen);
//		EPSIsHermitian(EPS eps,PetscBool *her);
//		EPSIsPositive(EPS eps,PetscBool *pos);

		EPSSetFromOptions(eps_);
	}

	/*
	 * Method: setEigenSolver
	 * ----------------------
	 *
	 */
	void setEigenSolver(const PetscInt nev, const PetscInt ncv, const PetscInt mpd, const double tol, const int max_iter) {
		assert(nev>=0 && ncv>=nev);

		// Initialize the eps context
		if(eps_empty)
			initEigenSolver();

		// Set the number of eigenvalues to compute
		EPSSetDimensions(eps_, nev, ncv, mpd); // nev = the number of eigenvalues to compute
		                                       //       The default is to compute only one.
		                                       // ncv = the number of column vectors to be used by the solution algorithm
		                                       //       It is recommended (depending on the method) to use ncv >= 2*nev or more.
		                                       // mpd = maximum projected dimension (for a more advanced usage -- Chaprter 2.6.4 in the SLEPc manual)
		                                       //       In the case of a large number of nev, the computation costs can be reduced by setting mpd << nev.

//		EPSSetLeftVectorsWanted(eps_, PETSC_TRUE); // To compute also left eigenvectors (For ver3.3, support for left eigenvectors is limited and will be extended in future versions of SLEPc)

		// Eigenvalues of interest
		EPSSetWhichEigenpairs(eps_, EPS_LARGEST_REAL);
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

//		PetscScalar tau;         // Target value (tau) -- Since the target (tau) is defined as a PetscScalar, complex values of tau are allowed only
//		                         //                       in the case of complex scalar builds of the SLEPc library
//		EPSSetTarget(eps_, tau); // If "EPS_TARGET_MAGNITUDE", "EPS_TARGET_REAL", or "EPS_TARGET_IMAGINARY" is used, set the target-value tau.

//		PetscScalar a;
//		PetscScalar b;
//		EPSSetInterval(eps_, a, b); // If "EPS_ALL" is used, set the interval [a, b]

		// Select the eigensolver
		EPSSetType(eps_, EPSKRYLOVSCHUR); // Default = EPSKRYLOVSCHUR
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

		// Control convergence
		EPSSetConvergenceTest(eps_, EPS_CONV_ABS); // Default is EPS_CONV_EIG. However, for computing eigenvalues close to the origin,
		                                           // this default criterion will give very poor accuracy.
		      // Available criterion: Absolute                :  EPS_CONV_ABS    ||r||
		      //                      Relative to eigenvalue  :  EPS_CONV_EIG    ||r|| / |lambda|
		      //                      Relative to matrix norms:  EPS_CONV_NORM   ||r|| / (||A|| + |lambda| ||B||)

		EPSSetTolerances(eps_, tol, max_iter); // Default value of tol = 1.0e-7
	}

	/*
	 * Method: solveAndGetEivenpairs
	 * -----------------------------
	 * Get some eigenvalue-eigenvector pairs (spectrum = largest Re(eigenvalues))
	 *
	 * Return:
	 *   By reference:
	 *     evalsReal = real part of the eigenvalues
	 *     evalsImag = imag part of the eigenvalues
	 *     evecsReal = real part of the eigenvectors
	 *     evecsImag = imag part of the eigenvectors
	 *     relError  = relative error to the eigenvalue or to the matrix norms for each eigenpair
	 *     numIter   = total number of iterations for the iterative eigen solver
	 *   By value:
	 *     number of converged eigen-pairs
	 */
	int solveAndGetEivenpairs(double *evalsReal, double *evalsImag, double **evecsReal, double **evecsImag, double *relError, int &numIter,
			const int *cvora, const int nVars, const PetscInt nev) {
		PetscErrorCode ierr;
#if PETSC_DEBUG_LEVEL > 1
		PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);
#endif

		// Initialize eigenvalues, eigenvectors, errors, and number of iterations
		assert(evalsReal!=NULL && evalsImag!=NULL && evecsReal!=NULL && evecsImag!=NULL && relError!=NULL);
		assert(cvora!=NULL);

		for(int i=0; i<nev; ++i)
			evalsReal[i] = 0.0;
		for(int i=0; i<nev; ++i)
			evalsImag[i] = 0.0;

		for(int i=0; i<nev; ++i)
			for(int j=0; j<(cvora[mpi_rank+1]-cvora[mpi_rank])*nVars; ++j)
				evecsReal[i][j] = 0.0;
		for(int i=0; i<nev; ++i)
			for(int j=0; j<(cvora[mpi_rank+1]-cvora[mpi_rank])*nVars; ++j)
				evecsImag[i][j] = 0.0;

		for(int i=0; i<nev; ++i)
			relError[i] = 0.0;

		numIter = 0;

		// If the number of eigenvalues to compute is zero, return 0
		if(nev==0)
			return 0;

		// Solve the eigen-problem
		ierr = EPSSolve(eps_);
		CHKERRQ(ierr);

		// Number of iterations required by the eigensolver
		EPSGetIterationNumber(eps_, &numIter);
#if PETSC_DEBUG_LEVEL > 1
		const EPSType type;
		PetscInt nev_, maxit;
		PetscReal tol;
		EPSGetType(eps_, &type);                                PetscPrintf(mpi_comm," Solution method                   : %s\n",type);
		EPSGetDimensions(eps_, &nev_, PETSC_NULL, PETSC_NULL);  PetscPrintf(mpi_comm," Number of requested eigenvalues   : %D\n", nev);
		EPSGetTolerances(eps_, &tol, &maxit);                   PetscPrintf(mpi_comm," Stopping condition                : tol=%.4G, maxit=%D\n", tol, maxit);
		EPSGetIterationNumber(eps_, &numIter);                  PetscPrintf(mpi_comm," Number of iterations of the method: %d\n", numIter);
		PetscPrintf(mpi_comm," Number of iterations of the method: %d\n", numIter);

		// Display solution
		ierr = EPSPrintSolution(eps_, PETSC_NULL); 	CHKERRQ(ierr);
#endif

		// Get the number of converged pairs (Note that the number of computed solutions depends on the convergence and,
		//                                    therefore, it may be different from the number of solutions requested by the user)
		PetscInt nconv;
		EPSGetConverged(eps_, &nconv);
#if PETSC_DEBUG_LEVEL > 1
		PetscPrintf(mpi_comm," Number of converged eigenpairs: %D\n\n", nconv);
#endif

		// Get the computed solution
		for (int iEigen=0; iEigen < std::min<int>(nev, nconv); ++iEigen) { // Note: If nev<ncv, nconv can be larger than nev
			EPSGetEigenpair(eps_, iEigen, &kr, &ki, xr, xi); // The eigenvectors are normalized so that they have a unit 2-norm,
			                                                 // except for problem type EPS_GHEP in which case returned eigenvectors
			                                                 // have a unit B-norm.
			    // Note:
			    //   Real SLEPc    -- kr and ki contain the real and imaginary parts of the eigenvalue, respectively,
			    //                    and xr and xi contain the associated eigenvector.
			    //   Complex SLEPc -- kr contains the (complex) eigenvalue and xr contains the corresponding (complex) eigenvector.
			    //                    In this case, ki and xi are not used (set to all zeros).

			ierr = EPSComputeRelativeError(eps_, iEigen, &relError[iEigen]); // error vector = A*x - lambda*x  or  A*x - lambda*B*x
			                                                                // Relative error = the norm of the error vector relative to the eigenvalue or to the matrix norms
			CHKERRQ(ierr);

			// Move the results to the solution variables to be returned
#if defined (PETSC_USE_COMPLEX)
			evalsReal[iEigen] = PetscRealPart(kr);
			evalsImag[iEigen] = PetscImaginaryPart(kr);
#else
			evalsReal[iEigen] = kr;
			evalsImag[iEigen] = ki;
#endif

#if defined(PETSC_USE_COMPLEX)
			for(int localIndex = 0; localIndex < (cvora[mpi_rank+1]-cvora[mpi_rank])*nVars; ++localIndex) {
				PetscScalar temp;
				int globalIndex = cvora[mpi_rank]*nVars + localIndex;
				VecGetValues(xr, 1, &globalIndex, &temp); // inefficient, need to be improved

				evecsReal[iEigen][localIndex] = PetscRealPart(temp);
				evecsImag[iEigen][localIndex] = PetscImaginaryPart(temp);
			}
#else
			for(int localIndex=0; localIndex<(cvora[mpi_rank+1]-cvora[mpi_rank])*nVars; ++localIndex) {
				int globalIndex = cvora[mpi_rank]*nVars + localIndex;
				VecGetValues(xr, 1, &globalIndex, &evecsReal[iEigen][localIndex]); // inefficient, need to be improved
			}
			for(int localIndex=0; localIndex<(cvora[mpi_rank+1]-cvora[mpi_rank])*nVars; ++localIndex) {
				int globalIndex = cvora[mpi_rank]*nVars + localIndex;
				VecGetValues(xi, 1, &globalIndex, &evecsImag[iEigen][localIndex]); // inefficient, need to be improved
			}
#endif
		}

		return nconv;
	}

private:
	void initSlepcSolver() {
		int argc = 0; char **argv;
		static char slepcHelp[] = "nothing!";
		SlepcInitialize(&argc, &argv, (char*)0, slepcHelp);

		eps_empty = PETSC_TRUE;

		alreadyHasMatrix = PETSC_FALSE;
	}

	PetscErrorCode finalizeSlepcSolver() {
		PetscErrorCode ierr;

//		PetscInt vecSize;
//		ierr = VecGetSize(xr, &vecSize); 	CHKERRQ(ierr);
//		if(vecSize != 0){
//			ierr = VecDestroy(&xr);  CHKERRQ(ierr);
//		}
//		ierr = VecGetSize(xi, &vecSize); 	CHKERRQ(ierr);
//		if(vecSize != 0){
//			ierr = VecDestroy(&xi);  CHKERRQ(ierr);
//		}

		if(!eps_empty) {
			ierr = EPSDestroy(&eps_); CHKERRQ(ierr);
			eps_empty = PETSC_TRUE;
		}

		ierr = SlepcFinalize();  CHKERRQ(ierr);

		return 0;
	}
};
#endif

#endif /* PETSCSOLVER2_H_ */
