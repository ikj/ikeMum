#ifndef PETSC_SOLVER_H
#define PETSC_SOLVER_H


#ifdef WITH_PETSC

#include "Logging.h"
using namespace logging;


#include "petscksp.h"

static bool ShowPetscMatrixOnScreenMatlab = false; // If this is true, the matrix will be shown on the screen in the MATLAB format
                                                   // right after the matrix is constructed in setLinSysForPetsc()

class PetscSolver
{
private:
  KSP ksp;        // linear solver context
  Vec x_, b_;     // solution, residual vector
  Mat A_;         // implicit operator matrix
  VecScatter scatterContext;

public:

  PetscSolver(int *noora, int *nbono_i, int *nbono_v, int m)
  {
    initPetscSolver(noora, nbono_i, nbono_v, m);
  }

  void setTresholds(double zeroAbs, double zeroRel, int maxIter)
  {
    KSPSetTolerances(ksp, zeroRel, zeroAbs, 1.e10, maxIter);
  }

  void setLinSysForPetsc( double *Ap, double *rhs,
                          int *noora, int *nbocv_i,int *nbocv_v, int *nbono_v_gl)
  {
    int nno = noora[mpi_rank+1] - noora[mpi_rank];

    // assemble rhs
    for (int ino=0; ino<nno; ++ino)
    {
      int row = noora[mpi_rank] + ino;
      VecSetValues(b_, 1, &row, &rhs[ino], INSERT_VALUES);
    }

    VecAssemblyBegin(b_);
    VecAssemblyEnd(b_);

    static int row_nnz_max = 0;
    if (row_nnz_max == 0)
    {
      for (int ino = 0; ino < nno; ino++)
      {
        int non_f = nbocv_i[ino];
        int non_l = nbocv_i[ino+1]-1;
        row_nnz_max = max(row_nnz_max, non_l-non_f+1);
      }
    }

    static int * cols = new int[row_nnz_max];
    static double *values = new double[row_nnz_max];
    int row;

    for (int ino=0; ino<nno; ++ino)
    {
      int non_f = nbocv_i[ino];
      int non_l = nbocv_i[ino+1]-1;
      int nnz = non_l-non_f+1;

      row = noora[mpi_rank] + ino;

      for (int non=non_f; non<=non_l; ++non)
      {
        int ino_nbr = nbocv_v[non];
        cols[non-non_f] = nbono_v_gl[ino_nbr];
      }

      for (int non=non_f; non<=non_l; ++non)
        values[non-non_f] = Ap[non];

      MatSetValues(A_, 1, &row, nnz, cols, values, INSERT_VALUES);
    }

    MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);

	// Show matrix A_ on screen !
	if(ShowPetscMatrixOnScreenMatlab) { // ShowPetscMatrixOnScreenMatlab is defined as a static variable at the beginning of this file and set as "false" by default
		PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);
		MatView(A_,PETSC_VIEWER_STDOUT_WORLD);
	}
  }

  void setLinSysForPetsc( double (*Ap)[5][5], double (*rhs)[5],
                          int *noora, int * nbocv_i,int * nbocv_v,
                          int *nbono_v_gl, int m)
  {
    int nno = noora[mpi_rank+1] - noora[mpi_rank];

    // assemble rhs
    for (int ino=0; ino<nno; ++ino)
    for (int i=0; i<m; ++i)
    {
      int row = noora[mpi_rank]*m + ino*m + i;
      VecSetValues(b_, 1, &row, &rhs[ino][i], INSERT_VALUES);
    }

    VecAssemblyBegin(b_);
    VecAssemblyEnd(b_);

    static int row_nnz_max = 0;
    if (row_nnz_max == 0)
    {
      for (int ino = 0; ino < nno; ino++)
      {
        int non_f = nbocv_i[ino];
        int non_l = nbocv_i[ino+1]-1;
        row_nnz_max = max(row_nnz_max, non_l-non_f+1);
      }
    }
    row_nnz_max *= m;

    static int * cols = new int[row_nnz_max];
    static int * rows = new int[m];
    static double *values = new double[row_nnz_max*m]; // times m for 5 rows

    for (int ino=0; ino<nno; ++ino)
    {
      int non_f = nbocv_i[ino];
      int non_l = nbocv_i[ino+1]-1;
      int nnz = (non_l-non_f+1)*m;

      for (int i=0; i<m; ++i)
        rows[i] = noora[mpi_rank]*m + ino*m + i;

      for (int non=non_f; non<=non_l; ++non)
      {
        int ino_nbr = nbocv_v[non];
        for (int i=0; i<m; ++i)
          cols[(non-non_f)*m+i] = nbono_v_gl[ino_nbr]*m + i;
      }

      for (int i=0; i<m; ++i)
      for (int non=non_f; non<=non_l; ++non)
      for (int j=0; j<m; ++j)
        values[nnz*i + (non-non_f)*m + j] = Ap[non][i][j];

      MatSetValues(A_, m, rows, nnz, cols, values, INSERT_VALUES);
    }

    MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);

	// Show matrix A_ on screen !
	if(ShowPetscMatrixOnScreenMatlab) { // ShowPetscMatrixOnScreenMatlab is defined as a static variable at the beginning of this file and set as "false" by default
		PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);
		MatView(A_,PETSC_VIEWER_STDOUT_WORLD);
	}
  }

   void setLinSysForPetscCoupled( double ***Ap, double **rhs,
                          int *noora, int * nbocv_i, int * nbocv_v,
                          int *nbono_v_gl, int m)
  {
    int nno = noora[mpi_rank+1] - noora[mpi_rank];

    // assemble rhs
    for (int ino=0; ino<nno; ++ino)
    for (int i=0; i<m; ++i)
    {
      int row = noora[mpi_rank]*m + ino*m + i;
      VecSetValues(b_, 1, &row, &rhs[ino][i], INSERT_VALUES);
    }

    VecAssemblyBegin(b_);
    VecAssemblyEnd(b_);

    static int row_nnz_max = 0;
    if (row_nnz_max == 0)
    {
      for (int ino = 0; ino < nno; ino++)
      {
        int non_f = nbocv_i[ino];
        int non_l = nbocv_i[ino+1]-1;
        row_nnz_max = max(row_nnz_max, non_l-non_f+1);
      }
    }
    row_nnz_max *= m;

    static int * cols = new int[row_nnz_max];
    static int * rows = new int[m];
    static double *values = new double[row_nnz_max*m]; // times m for 5 rows

    for (int ino=0; ino<nno; ++ino)
    {
      int non_f = nbocv_i[ino];
      int non_l = nbocv_i[ino+1]-1;
      int nnz = (non_l-non_f+1)*m;

      for (int i=0; i<m; ++i)
        rows[i] = noora[mpi_rank]*m + ino*m + i;

      for (int non=non_f; non<=non_l; ++non)
      {
        int ino_nbr = nbocv_v[non];
        for (int i=0; i<m; ++i)
          cols[(non-non_f)*m+i] = nbono_v_gl[ino_nbr]*m + i;
      }

      for (int i=0; i<m; ++i)
      for (int non=non_f; non<=non_l; ++non)
      for (int j=0; j<m; ++j)
        values[nnz*i + (non-non_f)*m + j] = Ap[non][i][j];

      MatSetValues(A_, m, rows, nnz, cols, values, INSERT_VALUES);
    }

    MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
  }



  void solveGMRES(double *Ap, double *phi, double *rhs,
      int *noora, int *nbono_i, int *nbono_v, int *nbono_v_gl)
  {
    setLinSysForPetsc(Ap, rhs, noora, nbono_i, nbono_v, nbono_v_gl);

    /*
    PC pc; // preconditioner context
    KSPGetPC(ksp,&pc);
    //PCSetType(pc, PCNONE);
    //PCSetType(pc, PCICC);
    //PCSetType(pc, PCLU);
    //PCSetType(pc, PCASM);
    //PCSetType(pc, PCJACOBI);
    PCSetType(pc, PCML);
    KSPSetOperators(ksp, A_, A_, SAME_NONZERO_PATTERN);
    KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);
    KSPSetType(ksp, KSPGMRES);
    //KSPSetType(ksp, KSPRICHARDSON);
    KSPGMRESSetRestart(ksp, 100);
    KSPSetFromOptions(ksp);
    */
    
     
    /*
    PC pc; // preconditioner context
    KSPGetPC(ksp,&pc);
    PCSetType(pc, PCMG);
    PetscOptionsSetValue("-pc_ml_maxNlevels", "3");
    PetscOptionsSetValue("-pc_ml_PrintLevel", "1");
    KSPSetOperators(ksp, A_, A_, SAME_NONZERO_PATTERN);
    KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);
    KSPSetType(ksp, KSPGMRES);
    //KSPSetType(ksp, KSPRICHARDSON);
    KSPGMRESSetRestart(ksp, 100);
    KSPSetFromOptions(ksp);
    */

    KSPSolve(ksp, b_, x_);

    int nIter;
    double relTol, absTol, dTol;

    KSPGetIterationNumber(ksp,&nIter);
    KSPGetResidualNorm(ksp,&absTol);

    int maxiter;
    double rel, abs, div;
    KSPGetTolerances(ksp, &rel, &abs, &div, &maxiter);
    lout(INFO_LO) << "tresholds 1: rel/abs/div/iter:\t" << rel << " " << abs << " " << div << " " << maxiter;
    double bNorm = 1.0;
    VecNorm(b_, NORM_2, &bNorm);
    lout(INFO_LO) << "\tPETSC: rel/abs/iter:\t" << absTol/bNorm << " " << absTol << " " << nIter << endl;
    // Note: absTol is (pre-conditioned) estimated norm, whereas bNorm is an exact norm
    //KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);


    int nno = noora[mpi_rank+1] - noora[mpi_rank];
    for (int ino = 0; ino < nno; ino++)
    {
      int row = noora[mpi_rank] + ino;
      VecGetValues(x_, 1, &row, &phi[ino]);   // inefficient, check how to change in future
    }

    KSPSetOperators(ksp, A_, A_, SAME_NONZERO_PATTERN);
  }


  void solveGMRES(double (*Ap)[5][5], double (*phi)[5], double (*rhs)[5],
      int *noora, int *nbono_i, int *nbono_v, int *nbono_v_gl, int m)
  {
    setLinSysForPetsc(Ap, rhs, noora, nbono_i, nbono_v, nbono_v_gl, m);

    /*
    PC pc; // preconditioner context
    KSPGetPC(ksp,&pc);
    //PCSetType(pc, PCILU);
    //PCSetType(pc, PCICC);
    PCSetType(pc, PCLU);
    //PCSetType(pc, PCASM);
    //PCSetType(pc, PCJACOBI);
    KSPSetOperators(ksp, A_, A_, SAME_NONZERO_PATTERN);
    KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);
    KSPSetType(ksp, KSPGMRES);
    //KSPSetType(ksp, KSPRICHARDSON);
    KSPGMRESSetRestart(ksp, 100);
    KSPSetFromOptions(ksp);
    */    
    

    KSPSolve(ksp, b_, x_);

    int nIter;
    double absTol;

    KSPGetIterationNumber(ksp,&nIter);
    KSPGetResidualNorm(ksp,&absTol);

    int maxiter;
    double rel, abs, div;
    KSPGetTolerances(ksp, &rel, &abs, &div, &maxiter);
    lout(INFO_LO) << "tresholds 2: rel/abs/div/iter:\t" << rel << " " << abs << " " << div << " " << maxiter;
    double bNorm = 1.0;
    VecNorm(b_, NORM_2, &bNorm);
    lout(INFO_LO) << "\tPETSC: rel/abs/iter:\t" << absTol/bNorm << " " << absTol << " " << nIter << endl;
    // Note: absTol is (pre-conditioned) estimated norm, whereas bNorm is an exact norm
    //KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);

    int nno = noora[mpi_rank+1] - noora[mpi_rank];
    for (int ino = 0; ino < nno; ino++)
    for (int i=0; i<m; ++i)
    {
      int row = noora[mpi_rank]*m + ino*m + i;
      VecGetValues(x_, 1, &row, &phi[ino][i]);
    }

    KSPSetOperators(ksp, A_, A_, SAME_NONZERO_PATTERN);
  }

void solveGMRESCoupled(double ***Ap, double **phi, double **rhs,
      int *noora, int *nbono_i, int *nbono_v, int *nbono_v_gl, int m)
  {

    setLinSysForPetscCoupled(Ap, rhs, noora, nbono_i, nbono_v, nbono_v_gl, m);

    
    /*
    PC pc; // preconditioner context
    KSPGetPC(ksp,&pc);
    //PCSetType(pc, PCICC);
    PCSetType(pc, PCLU);
    //PCSetType(pc, PCASM);
    //PCSetType(pc, PCJACOBI);
    KSPSetOperators(ksp, A_, A_, SAME_NONZERO_PATTERN);
    KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);
    KSPSetType(ksp, KSPGMRES);
    //KSPSetType(ksp, KSPRICHARDSON);
    KSPGMRESSetRestart(ksp, 100);
    KSPSetFromOptions(ksp);
    */
    

    KSPSolve(ksp, b_, x_);
    
    int nIter;
    double absTol;

    KSPGetIterationNumber(ksp,&nIter);
    KSPGetResidualNorm(ksp,&absTol);

    int maxiter;
    double rel, abs, div; 
    KSPGetTolerances(ksp, &rel, &abs, &div, &maxiter);
    lout(INFO_LO) << "tresholds 3: rel/abs/div/iter:\t" << rel << " " << abs << " " << div << " " << maxiter;
    lout(INFO_LO) << "\tPETSC: iter/absResidual:\t" << nIter << "\t" << absTol << endl;

    
    
    int nno = noora[mpi_rank+1] - noora[mpi_rank];
    for (int ino = 0; ino < nno; ino++)
    for (int i=0; i<m; ++i)
    {
      int row = noora[mpi_rank]*m + ino*m + i;
      VecGetValues(x_, 1, &row, &phi[ino][i]);
    }

    KSPSetOperators(ksp, A_, A_, SAME_NONZERO_PATTERN); 
  }


private:

  void initPetscSolver( int *noora, int *nbono_i, int *nbono_v, int m)
  {
    int argc = 0; char **argv;
    static char help[] = "nothing!";
    PetscInitialize(&argc,&argv,(char *)0,help);

    PC pc; // preconditioner context
    PetscErrorCode ierr;

    // Create linear solver context
    KSPCreate(mpi_comm, &ksp);

    int nno = noora[mpi_rank+1] - noora[mpi_rank];

    // allcoate mem for vec
    VecCreateMPI(mpi_comm, nno*m, noora[mpi_size]*m, &b_);
    VecSetFromOptions(b_);
    VecDuplicate(b_, &x_);
    VecSet(b_, 0.0);
    VecSet(x_, 0.0);

    // allcoate mem for implicit operator
    int * diag_sizes = new int[nno*m];
    int * offd_sizes = new int[nno*m];

    for (int ino=0; ino<nno; ino++) //nnoS
    {
     diag_sizes[ino*m] = offd_sizes[ino*m] = 0;

     int non_f = nbono_i[ino];
     int non_l = nbono_i[ino+1]-1;

     for (int non=non_f; non<=non_l; non++)
     {
       int ino_nbr = noora[mpi_rank] + nbono_v[non];

       if ((ino_nbr >= noora[mpi_rank]) && (ino_nbr < noora[mpi_rank+1]))
         diag_sizes[ino*m] += m;
       else
         offd_sizes[ino*m] += m;
     }

     for (int i=1; i<m; i++)
     {
       diag_sizes[ino*m+i] = diag_sizes[ino*m];
       offd_sizes[ino*m+i] = offd_sizes[ino*m];
     }
    }

    MatCreateMPIAIJ(mpi_comm,
       nno*m, nno*m, noora[mpi_size]*m, noora[mpi_size]*m,
       0, &diag_sizes[0], 0, &offd_sizes[0], &A_);

    delete [] diag_sizes;
    delete [] offd_sizes;

   
    // Case 1
    
     /*
    KSPSetOperators(ksp, A_, A_, SAME_NONZERO_PATTERN);
    KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);
    //PetscOptionsSetValue("-ksp_initial_guess_nonzero", "true");
    KSPSetType(ksp, KSPGMRES);
    KSPGMRESSetRestart(ksp, 200);
    PetscOptionsSetValue("-pc_type","bjacobi");
    PetscOptionsSetValue("-sub_pc_type","lu");
    PetscOptionsSetValue("-sub_ksp_rtol","1.e-15");
    PetscOptionsSetValue("-ksp_monitor_true_residual","resid");
    PetscOptionsSetValue("-ksp_monitor_singular_value","singval");
    KSPSetFromOptions(ksp);
    */
    
    
    

    /*
    KSPSetOperators(ksp, A_, A_, SAME_NONZERO_PATTERN);
    KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);
    PetscOptionsSetValue("-ksp_type","gmres");
    KSPGMRESSetRestart(ksp, 100);
    //PetscOptionsSetValue("-pc_type","ilu");
    //PetscOptionsSetValue("-pc_factor_levels","8");
    PetscOptionsSetValue("-pc_type","lu");
    //PetscOptionsSetValue("-pc_factor_mat_solver_package","superlu_dist");
    PetscOptionsSetValue("-ksp_monitor_true_residual","resid");
    PetscOptionsSetValue("-ksp_monitor_singular_value","singval");
    KSPSetFromOptions(ksp);
    */
   
    
    
    

   /*
    KSPSetOperators(ksp, A_, A_, SAME_NONZERO_PATTERN);
    KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);
    KSPGetPC(ksp,&pc);
    PCSetType(pc, PCHYPRE);
    PetscOptionsSetValue("-pc_hypre_type", "euclid");
    PetscOptionsSetValue("-pc_hypre_euclid_levels", "2");
    PetscOptionsSetValue("-pc_hypre_euclid_print_statistics","true");
    KSPSetType(ksp, KSPFGMRES);
    KSPGMRESSetRestart(ksp, 100);
    PetscOptionsSetValue("-ksp_monitor_true_residual","resid");
    PetscOptionsSetValue("-ksp_monitor_singular_value","singval");
    KSPSetFromOptions(ksp);
    */
   
   

    // Case 0 
    
    
    KSPSetOperators(ksp, A_, A_, SAME_NONZERO_PATTERN);
    KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);
    KSPSetType(ksp, KSPGMRES);
    KSPGMRESSetRestart(ksp, 100);
    //PetscOptionsSetValue("-ksp_monitor_true_residual","1");
    KSPSetFromOptions(ksp);
   
    
    
    
    
    
   
    
    
   
    

  }
};

#else

// if the user has not built the hypre libraries, provide an empty
// HypreSolver class that errors is instantiated. This means you
// can still explicitly include a hypre solver as a POINTER in your
// solver, but do not instantiate it...

class PetscSolver {

public:

  PetscSolver(int *noora, int *nbono_i, int *nbono_v, int m)
  {
    std::cerr << "Error: to use petsc you need to compile with -D WITH_PETSC." << std::endl;
    throw(-1);
  }

  // only need to provide blank interface to public methods...

  void setTresholds(double zeroAbs, double zeroRel, int maxIter) {}

  void solveGMRES(double *Ap, double *phi, double *rhs,
        int *noora, int *nbono_i, int *nbono_v, int *nbono_v_gl) {}

  void solveGMRES(double (*Ap)[5][5], double (*phi)[5], double (*rhs)[5],
        int *noora, int *nbono_i, int *nbono_v, int *nbono_v_gl, int m) {}

    void solveGMRESCoupled(double ***Ap, double **phi, double **rhs,
        int *noora, int *nbono_i, int *nbono_v, int *nbono_v_gl, int m) {}



};

#endif


#endif
