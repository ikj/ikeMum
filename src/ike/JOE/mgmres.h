void atx_cr ( int n, int nz_num, double a[], int ia[], int ja[], double x[], 
  double w[] );
void atx_st ( int n, int nz_num, double a[], unsigned int ia[], unsigned int ja[], double x[], 
  double w[] );
void ax_cr ( int n, int nz_num, double a[], int ia[], int ja[], double x[], 
  double w[] );
void ax_st ( int n, int nz_num, double a[], unsigned int ia[], unsigned int ja[], double x[], 
  double w[] );
void diagonal_pointer_cr ( int n, int nz_num,  int ia[],  int ja[], int ua[] );
void lus_cr ( int n, int nz_num,  int ia[],  int ja[], double l[], int ua[], 
  double r[], double z[] );
void mgmres_st ( int n, int nz_num, unsigned int ia[], unsigned int ja[], double a[], double x[],
  double rhs[], int itr_max, int mr, double tol_abs, double tol_rel );
void mult_givens ( double c, double s, int k, double g[] );
void pmgmres_ilu_cr ( int n, int nz_num, int ia[], int ja[], double a[], 
  double x[], double rhs[], int itr_max, int mr, double tol_abs, 
  double tol_rel );
void ilu_cr ( int n, int nz_num, int ia[], int ja[], double a[], int ua[],
  double l[] );
double r8vec_dot ( unsigned int n, double a1[], double a2[] );
double r8vec_dot ( int n, double a1[], double a2[] );
double *r8vec_uniform_01 ( int n, int *seed );
void rearrange_cr ( int n, int nz_num, int ia[], int ja[], double a[] );
void timestamp ( void );

void sort(unsigned int *col_idx, double *a, int start, int end);
void coo2csr_in(int n, int nz, double *a, unsigned int *i_idx, unsigned int *j_idx);
void csr2csc(int n, int m, int nz, double *a, unsigned int *col_idx, unsigned int *row_start,
             double *csc_a, int *row_idx, int *col_start);
void sort(int *col_idx, double *a, int start, int end);
void coo2csr_in(int n, int nz, double *a, int *i_idx, int *j_idx);


