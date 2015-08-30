#include "DiffFilter.h"
#define DONT_SHOW_CG_DETAILS

/*
 * Method: filterConsturcted
 * -------------------------
 *
 */
bool DiffFilter::filterConsturcted() {
	if(cvDf==NULL)
		return false;
	return true;
}


/*
 * Method: appplyDiffFilterG1G2
 * ----------------------------
 *
 */
void DiffFilter::applyDiffFilterG1G2(double *phif, const double *phi, const double absTol, const int maxIter) {
	assert(cvDf != NULL);
	assert(phif != NULL && phi != NULL);

	double *rhs = new double [ncv];

	FOR_ICV {
		rhs[icv] = cv_volume[icv] * phi[icv];
		phif[icv] = phi[icv]; // initial guess
	}
	updateCvDataG1G2(phif,REPLACE_DATA);

	solveCvScalarCgG1G2(phif, cvDf, rhs, ABSOLUTE_RESIDUAL, absTol, maxIter);
	delete [] rhs;
}

/*
 * Method: appplyDiffFilterG1G2
 * ----------------------------
 *
 */
void DiffFilter::applyDiffFilterG1G2(double (*phif)[3], const double (*phi)[3], const double absTol, const int maxIter) {
	assert(cvDf != NULL);
	assert(phif != NULL && phi != NULL);

	double *rhs = new double[ncv];
	double *soln1D = new double[ncv_ggff];
	FOR_I3 {
		FOR_ICV {
			rhs[icv] = cv_volume[icv] * phi[icv][i];
			soln1D[icv] = phi[icv][i]; // initial guess
		}
		updateCvDataG1G2(soln1D, REPLACE_DATA);

		solveCvScalarCgG1G2(soln1D, cvDf, rhs, ABSOLUTE_RESIDUAL, absTol, maxIter);

		FOR_ICV {
			phif[icv][i] = soln1D[icv];
		}
	}

	delete [] rhs;
	delete [] soln1D;
}

/*
 * Method: specify_filter_width
 * ----------------------------
 *     delta^2
 * p = -------, where delta = filter width, if p is constant everywhere
 *      40.0
 */
double DiffFilter::specify_filter_width(const int ifa) {
	return 1.0;
}

/*
 * Method: buildCvDifferentialFilter
 * ---------------------------------
 *
 */
void DiffFilter::buildCvDifferentialFilter() {
	if (mpi_rank == 0)
		cout << "DiffFilter::buildDifferentialFilter()"<<endl;

	assert( cvDf == NULL );
	cvDf = new double[nbocv_s]; // note: nbocv_s == size of nbocv_v

	for (int i = 0 ; i < nbocv_s ; ++i)
		cvDf[i] = 0.0;

	// now add the scaled laplacian ..
	FOR_IFA_NONB {
		int icv0 = cvofa[ifa][0];
		int icv1 = cvofa[ifa][1];

		assert((icv0 >= 0) && (icv0 < ncv));
		assert((icv1 >= 0) && (icv1 < ncv_g));

		int noc00, noc01, noc11, noc10;
		getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);
				/* Note: noc00 -> diagonal element of the icv0 row
				 *       noc11 -> diagonal element of the icv1 row
				 *       noc01 -> flux contribution of icv1 to the icv0 balance
				 *       noc10 -> flux contribution of icv0 to the icv1 balance
				 */
		double unit_s[3];
		vecMinVec3d( unit_s, x_cv[icv1], x_cv[icv0]) ;
		double smag = normVec3d(unit_s) ;

		// determine the function p(x) -- this is related to the filter width  px = Delta_f^2/40 (approx) ..
		double px = specify_filter_width(ifa);

		double sdotn = 0.0 ;
		FOR_I3 {
			sdotn += unit_s[i] * fa_normal[ifa][i];
		}

		// simple gradient at the face for now ... there is a more
		// involved procedure based on the way that the viscous fluxes are calculated
		cvDf[noc00] -= px* sdotn / smag;
		cvDf[noc01] += px* sdotn/ smag ;
		if ( icv1 < ncv ) {
			cvDf[noc11] -= px* sdotn/ smag ;
			cvDf[noc10] += px* (sdotn/smag);
		} 

		// add further terms here
	}

	for (int i=0 ; i < nbocv_s ; ++i)
		cvDf[i] *= -1.0;   // flip the sign...

	FOR_ICV {
		int noc_f = nbocv_i[icv];
		cvDf[noc_f] += cv_volume[icv];
	}

	MPI_Barrier(mpi_comm);
}

/*
 * Method: solveCvScalarCgG1G2
 * ---------------------------
 * conjugate gradient
 * Original code = solveCvScalarCg in UgpWithCv2.h (original version)
 */
int DiffFilter::solveCvScalarCgG1G2(double *phi, const double *Ap, const double *rhs, const int mode, const double zero, const int maxiter) {
	// we need the following work arrays...
	double *res = new double[ncv];
	double *v   = new double[ncv];
	double *p   = new double[ncv_ggff];

	// initialize...
	for (int icv = 0; icv < ncv_ggff; icv++)
		p[icv] = 0.0;

	double rho = 1.0;

	// calculate the residual in rhs format...
	for (int icv = 0; icv < ncv; icv++) {
		int noc_f = nbocv_i[icv];
		int noc_l = nbocv_i[icv + 1] - 1;
		res[icv] = Ap[noc_f] * phi[icv]; // diagonal
		// cv neighbors...
		for (int noc = noc_f + 1; noc <= noc_l; noc++)
			res[icv] += Ap[noc] * phi[nbocv_v[noc]];
		res[icv] = rhs[icv] - res[icv];
	}

	double res0_max;
	if (mode == RELATIVE_RESIDUAL) {
		double my_res0_max = 0.0;

		for (int icv = 0; icv < ncv; icv++)
			my_res0_max = max(my_res0_max, fabs(res[icv] / Ap[nbocv_i[icv]]));

		MPI_Reduce(&my_res0_max, &res0_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);

#ifndef DONT_SHOW_CG_DETAILS
		if(mpi_rank==0)
			printf("cg initial   res_max(l-inf norm)= %.5e \n", res0_max);
#endif
	}

	int iter = 0;
	int done = 0;
	while (done == 0) {
		iter++;

		// diagonal precon...
		for (int icv = 0; icv < ncv; icv++)
			v[icv] = res[icv] / Ap[nbocv_i[icv]];

		double rho_prev = rho;

		double my_rho = 0.0;
		for (int icv = 0; icv < ncv; icv++)
			my_rho += res[icv] * v[icv];

		MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

		if (fabs(rho_prev) < 1.0E-20)
			rho_prev = 1.0E-20;
		double beta = rho / rho_prev;

		for (int icv = 0; icv < ncv; icv++)
			p[icv] = v[icv] + beta * p[icv];
		updateCvDataG1G2(p, REPLACE_DATA);

		// v = [Ap]{p}...
		for (int icv = 0; icv < ncv; icv++) {
			int noc_f = nbocv_i[icv];
			int noc_l = nbocv_i[icv + 1] - 1;
			v[icv] = Ap[noc_f] * p[icv]; // diagonal
			for (int noc = noc_f + 1; noc <= noc_l; noc++)
				v[icv] += Ap[noc] * p[nbocv_v[noc]];
		}

		double my_gamma = 0.0;
		for (int icv = 0; icv < ncv; icv++)
			my_gamma += p[icv] * v[icv];
		double gamma;
		MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

		if (fabs(gamma) < 1.0E-20)
			gamma = 1.0E-20;

		double alpha = rho / gamma;
		for (int icv = 0; icv < ncv_gg; icv++)
			phi[icv] += alpha * p[icv];

		// check if we are done...
		if (iter % 5 == 0) {

			double my_res_max = 0.0;

			// recompute the residual...
			updateCvDataG1G2(phi, REPLACE_DATA);
			for (int icv = 0; icv < ncv; icv++) {
				int noc_f = nbocv_i[icv];
				int noc_l = nbocv_i[icv + 1] - 1;
				res[icv] = Ap[noc_f] * phi[icv]; // diagonal
				for (int noc = noc_f + 1; noc <= noc_l; noc++)
					res[icv] += Ap[noc] * phi[nbocv_v[noc]];
				// above is LHS. residual is then...
				res[icv] = rhs[icv] - res[icv];
				my_res_max = max(my_res_max, fabs(res[icv] / Ap[nbocv_i[icv]]));
			}

			double res_max;
			MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
			if (mpi_rank == 0) {
#ifndef DONT_SHOW_CG_DETAILS
				printf("cg iter=%3d, res_max(l-inf norm)= %.5e \n", iter, res_max);
#endif

				if ((mode == ABSOLUTE_RESIDUAL) && (res_max <= zero)) {
					done = 1;
				} else if ((mode == RELATIVE_RESIDUAL) && (res_max / (res0_max + 1.0E-20) <= zero)) {
					done = 1;
				} else if (iter > maxiter) {
					cout << "Warning: solveCvScalarCg did not converge after " << maxiter << " iters, res_max: " << res_max
							<< endl;
					done = 2;
				}
			}
			MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

		} else {
			// on the other iterations, use this approximation...
			for (int icv = 0; icv < ncv; icv++)
				res[icv] -= alpha * v[icv];
		}
	}

	//delete[]
	delete[] res;
	delete[] v;
	delete[] p;

	// let the calling routine know if we were successful...
	return (done == 1);
}
