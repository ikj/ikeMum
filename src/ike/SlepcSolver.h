/*
 * SlepcSolver.h
 *
 *  Created on: Nov 26, 2013
 *      Author: ikj
 */

#ifndef SLEPCSOLVER_H_
#define SLEPCSOLVER_H_

#include "Logging.h"
using namespace logging;
#include "petscksp.h"

#include "MatComprsed.h"

#define SLEPC_DEBUG_LEVEL 1

#ifndef WITH_PETSC
#define WITH_PETSC
#endif

#define SLEPC_DEBUG_LEVEL 1

#ifdef USE_SLEPC_WITH_PETSC
	//#include "slepcsys.h"
	#include <slepceps.h>
	#include <slepcsvd.h>
#endif

/*
 * Class: SlepcSolver
 * ------------------
 *
 */
class SlepcSolver
{
public:
	/*
	 * Constructor
	 */
	SlepcSolver(int *cvora, int *nbocv2_i, vector<int> &nbocv2_v, int m, int NcontrolParams=0, int monitorConvergInterval=0) {
#if SLEPC_DEBUG_LEVEL > 1
		if(mpi_rank==0) cout<<"SlepcSolver()"<<endl;
#endif
		initSlepcSolver();
	}

	/*
	 * Destructor
	 */
	~SlepcSolver() {
#if PETSC_DEBUG_LEVEL > 1
		if(mpi_rank==0) cout<<"~SlepcSolver()"<<endl;
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

public:
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
	 * Method: setLinSysForPetsc
	 * -------------------------
	 * Set the linear system: basically upload the matrix to the memory in the form of the PETSC multi-core.
	 * Called by the solveEigenProblem() method
	 * Note: 1. Matrix structures (MatComprsed and MatComprsed) are defined in MatComprsed.h
	 */
	template <class MatT>
	void setJacobianMatrixForSlepc(MatT &A, const int *cvora, const int *cv_gl, const int nVars, const int ncv_gg) {


		int nno = (int) (A.get_nRows()/nVars); // local number of CVs
		assert( nno == cvora[mpi_rank+1] - cvora[mpi_rank] );

		bool isOnlyJacMatrix = (A.get_nCols() == ncv_gg*nVars); // Check if the matrix is exactly a Jacobian matrix
		                                                        // (Note that the size of a Pseudo-arclength matrix is bigger than a Jacobian matrix)
		if(!isOnlyJacMatrix && mpi_rank==0)
			cout<<"WARNING in SlepcSolvers2::setJacobianMatrixForSlepc(): Given matrix is bigger than Jacobian"<<endl;

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
		MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);

		// Show matrix A_ on screen !
		if(ShowPetscMatrixMatlab) { // ShowPetscMatrixMatlab is defined as a static variable at the beginning of this file and set as "false" by default
			PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_MATLAB);
			MatView(A_,PETSC_VIEWER_STDOUT_WORLD);
		}

		MPI_Allreduce(&myEmptyRows, &emptyRows, 1, MPI_INT, MPI_SUM, mpi_comm);
		if(emptyRows>0 && mpi_rank == 0)
			cout<<"CAUTION in SlepcSolvers2::setJacobianMatrixForSlepc(): total "<<emptyRows<<" rows are empty"<<endl;

		MPI_Allreduce(&myNonJacElements, &nonJacElements, 1, MPI_INT, MPI_SUM, mpi_comm);
		if(nonJacElements>0)
			if(nonJacElements != cvora[mpi_size]*nVars && mpi_rank == 0)
				cout<<"CAUTION in SlepcSolvers2::setJacobianMatrixForSlepc(): total "<<nonJacElements<<" elements from the given Jacobian matrix will be omitted"<<endl;

		MPI_Barrier(mpi_comm);

ShowPetscMatrixMatlab = false;
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
		setJacobianMatrixForSlepc(A, cvora, cv_gl, 5+nScal, ncv_gg);

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
		PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);

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

		const EPSType type;
		PetscInt nev_, maxit;
		PetscReal tol;
		EPSGetType(eps_, &type);                                PetscPrintf(mpi_comm," Solution method                   : %s\n",type);
		EPSGetDimensions(eps_, &nev_, PETSC_NULL, PETSC_NULL);  PetscPrintf(mpi_comm," Number of requested eigenvalues   : %D\n", nev);
		EPSGetTolerances(eps_, &tol, &maxit);                   PetscPrintf(mpi_comm," Stopping condition                : tol=%.4G, maxit=%D\n", tol, maxit);
		EPSGetIterationNumber(eps_, &numIter);                  PetscPrintf(mpi_comm," Number of iterations of the method: %d\n", numIter);

		// Number of iterations required by the eigensolver
		ierr = EPSGetIterationNumber(eps_, &numIter);
		CHKERRQ(ierr);

		// Display solution
		ierr = EPSPrintSolution(eps_, PETSC_NULL);
		CHKERRQ(ierr);

		// Get the number of converged pairs (Note that the number of computed solutions depends on the convergence and,
		//                                    therefore, it may be different from the number of solutions requested by the user)
		PetscInt nconv;
		EPSGetConverged(eps_, &nconv);
		PetscPrintf(mpi_comm," Number of converged eigen-pairs: %D\n\n", nconv);

if(mpi_rank==0)
	cout<<"solveAndGetEivenpairs()-1: "<<endl;

//		// Get the computed solution
//		for (int iEigen=0; iEigen<nconv; iEigen++) {
//			PetscInt iTemp = iEigen;
//			EPSGetEigenpair(eps_, iTemp, &kr, &ki, xr, xi); // The eigenvectors are normalized so that they have a unit 2-norm,
//			                                                 // except for problem type EPS_GHEP in which case returned eigenvectors
//			                                                 // have a unit B-norm.
//			    // Note:
//			    //   Real SLEPc    -- kr and ki contain the real and imaginary parts of the eigenvalue, respectively,
//			    //                    and xr and xi contain the associated eigenvector.
//			    //   Complex SLEPc -- kr contains the (complex) eigenvalue and xr contains the corresponding (complex) eigenvector.
//			    //                    In this case, ki and xi are not used (set to all zeros).
//
//			PetscReal error;
//			ierr = EPSComputeRelativeError(eps_, iTemp, &error); // error vector = A*x - lambda*x  or  A*x - lambda*B*x
//			                                                      // Relative error = the norm of the error vector relative to the eigenvalue or to the matrix norms
//			CHKERRQ(ierr);
//			relError[iEigen] = error;
//
//			// Move the results to the solution variables to be returned
//#if defined (PETSC_USE_COMPLEX)
//			PetscReal re = PetscRealPart(kr);
//			PetscReal im = PetscImaginaryPart(kr);
//			PetscPrintf(mpi_comm," %9F%+9Fj %12G\n", re, im, error);
//
//			evalsReal[iEigen] = PetscRealPart(kr);
//			evalsImag[iEigen] = PetscImaginaryPart(kr);
//#else
//			PetscPrintf(mpi_comm," %9F%+9Fj %12G\n", kr, ki, error);
//			VecView(xr, viewer);
//
//			evalsReal[iEigen] = kr;
//			evalsImag[iEigen] = ki;
//#endif
//
////#if defined(PETSC_USE_COMPLEX)
////			for(int localIndex = 0; localIndex < (cvora[mpi_rank+1]-cvora[mpi_rank])*nVars; ++localIndex) {
////				PetscScalar temp;
////				int globalIndex = cvora[mpi_rank]*nVars + localIndex;
////				VecGetValues(xr, 1, &globalIndex, &temp); // inefficient, need to be improved
////
////				evecsReal[iEigen][localIndex] = PetscRealPart(temp);
////				evecsImag[iEigen][localIndex] = PetscImaginaryPart(temp);
////			}
////#else
////			for(int localIndex=0; localIndex<(cvora[mpi_rank+1]-cvora[mpi_rank])*nVars; ++localIndex) {
////				int globalIndex = cvora[mpi_rank]*nVars + localIndex;
////				VecGetValues(xr, 1, &globalIndex, &evecsReal[iEigen][localIndex]); // inefficient, need to be improved
////			}
////			for(int localIndex=0; localIndex<(cvora[mpi_rank+1]-cvora[mpi_rank])*nVars; ++localIndex) {
////				int globalIndex = cvora[mpi_rank]*nVars + localIndex;
////				VecGetValues(xi, 1, &globalIndex, &evecsImag[iEigen][localIndex]); // inefficient, need to be improved
////			}
////#endif
//		}
if(mpi_rank==0)
	cout<<"solveAndGetEivenpairs()-2"<<endl;

		return nconv;
	}

private:
	void initSlepcSolver() {
		int argc = 0; char **argv;
		static char slepcHelp[] = "nothing!";
		SlepcInitialize(&argc, &argv, (char*)0, slepcHelp);

		eps_empty = PETSC_TRUE;
	}

	PetscErrorCode finalizeSlepcSolver() {
		PetscErrorCode ierr;

		ierr = VecDestroy(&xr);  CHKERRQ(ierr);
		ierr = VecDestroy(&xi);  CHKERRQ(ierr);
		if(!eps_empty) {
			ierr = EPSDestroy(&eps_); CHKERRQ(ierr);
			eps_empty = PETSC_TRUE;
		}
		ierr = SlepcFinalize();  CHKERRQ(ierr);

		return 0;
	}
};

#endif /* SLEPCSOLVER_H_ */
