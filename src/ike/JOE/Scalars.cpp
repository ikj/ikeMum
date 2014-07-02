#include "UgpWithCvCompFlow.h"



void UgpWithCvCompFlow::setScalarBC(FaZone *zone)
{
  for (ScalarTranspEqIterator scal = scalarTranspEqVector.begin(); scal < scalarTranspEqVector.end(); scal++)
  {
    string scalName(scal->getName());
    double phiBCval = 0.0, phiBCflux = 0.0;

    double *phi = scal->phi;
    double (*grad_phi)[3] = scal->grad_phi;

    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      // .............................................................................................
      // HOOK BOUNDARY CONDITION
      // .............................................................................................
      if (scalarZoneIsHook(zone->getName(), scalName))
      {
        // user defined boundary hook
        boundaryHookScalarRansTurb(phi, &(*zone), scalName);
        boundaryHookScalarRansComb(phi, &(*zone), scalName);
      }
      // .............................................................................................
      // DIRICHLET BOUNDARY CONDITION
      // .............................................................................................
      else if (scalarZoneIsDirichlet(phiBCval, zone->getName(), scalName))
      {
      //for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      //{
        for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
          int icv0 = cvofa[ifa][0];
          int icv1 = cvofa[ifa][1];
          if (vecDotVec3d(vel[icv1], fa_normal[ifa]) > 1.0e-8 && scal->convTerm)
          {
            double r[3] = {0.0, 0.0, 0.0};
            vecMinVec3d(r, x_fa[ifa], x_cv[icv0]);
            phi[icv1] = phi[icv0];// + vecDotVec3d(r, grad_phi[icv0]);
          }
          else
            phi[icv1] = phiBCval;
        }
      }
      // .............................................................................................
      // FLUX BOUNDARY CONDITION: NOT IMPLEMENTED
      // .............................................................................................
      else if (scalarZoneIsFlux(phiBCflux, zone->getName(), scalName))
      {
        cerr << "scalarZoneIsFlux not implemented yet!" << endl;
        throw(-1);
      }
      // .............................................................................................
      // OTHER BOUNDARY CONDITIONS: NO FLUX
      // .............................................................................................
      else // apply zero flux if no boundary condition is specified
      {
        //for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
        //{
        for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
          int icv0 = cvofa[ifa][0];
          int icv1 = cvofa[ifa][1];

          double r[3] = {0.0, 0.0, 0.0};
          vecMinVec3d(r, x_fa[ifa], x_cv[icv0]);
 	  phi[icv1] = phi[icv0];
        }
      }
    }
  }
}


void UgpWithCvCompFlow::calcViscousFluxScalar(double *rhs_rhoScal, double *AScal, ScalarTranspEq &transpScal, int flagImplicit) {
	if (!transpScal.diffTerm)
		return;

	string scalName(transpScal.getName());
	double phiBCval = 0.0, phiBCflux = 0.0;

	double *phi = transpScal.phi;
	double (*grad_phi)[3] = transpScal.grad_phi;
	double *diff = transpScal.diff;

	// =============================================================================================
	// compute gradients, with boundary values
	// =============================================================================================
	calcCv2Grad(grad_phi, phi, limiterNavierS, phi, epsilonSDWLS);

	// =============================================================================================
	// compute user-defined diffusivity
	// =============================================================================================
	diffusivityHookScalarRansTurb(scalName);
	diffusivityHookScalarRansComb(scalName);

	// =============================================================================================
	// cycle trough internal faces and compute viscous flux and implicit matrix
	// =============================================================================================
	for (int ifa = nfa_b; ifa < nfa; ifa++)
	{
		int icv0 = cvofa[ifa][0];
		int icv1 = cvofa[ifa][1];
		int noc00, noc01, noc11, noc10;
		if (flagImplicit)
			getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

		double n[3] = {0.0, 0.0, 0.0}, s[3] = {0.0, 0.0, 0.0};
		double nmag = normVec3d(n, fa_normal[ifa]);
		vecMinVec3d(s, x_cv[icv1], x_cv[icv0]);
		double smag = normVec3d(s);
		double alpha = vecDotVec3d(n, s);
		assert((alpha > 0.0) && (alpha < 1.0+1.0E-12));   // alpha should now contain s dot n...

		double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
		vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
		vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);

		double viscFlux = diff[ifa] * nmag * (alpha * (phi[icv1] - phi[icv0]) / smag
				+ 0.5 * (grad_phi[icv0][0] + grad_phi[icv1][0]) * (n[0] - alpha * s[0])
				+ 0.5 * (grad_phi[icv0][1] + grad_phi[icv1][1]) * (n[1] - alpha * s[1])
				+ 0.5 * (grad_phi[icv0][2] + grad_phi[icv1][2]) * (n[2] - alpha * s[2]));

		rhs_rhoScal[icv0] += viscFlux;

		if (icv1 < ncv)
			rhs_rhoScal[icv1] -= viscFlux;

		if (flagImplicit)
		{
			AScal[noc00] +=   diff[ifa] * nmag * alpha / smag / rho[icv0];
			AScal[noc01] += - diff[ifa] * nmag * alpha / smag / rho[icv1];

			if (icv1 < ncv)
			{
				AScal[noc11] +=   diff[ifa] * nmag * alpha / smag / rho[icv1];
				AScal[noc10] += - diff[ifa] * nmag * alpha / smag / rho[icv0];
			}
		}
	}

	// =============================================================================================
	// cycle trough boundary faces and compute viscous flux and implicit matrix
	// =============================================================================================
	for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
	{
		if (zone->getKind() == FA_ZONE_BOUNDARY)
		{
			// .............................................................................................
			// HOOK BOUNDARY CONDITION
			// .............................................................................................
			if (scalarZoneIsHook(zone->getName(), scalName))
			{
				if ((mpi_rank == 0) && firstCall_JOEcalcViscousFluxScalar)
					cout << scalName << ": " << zone->getName() << " is HOOK" << endl;

				for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
				{
					int icv0 = cvofa[ifa][0];
					int icv1 = cvofa[ifa][1];
					int noc00 = nbocv_i[icv0];

					double n[3] = {0.0, 0.0, 0.0}, s_half[3] = {0.0, 0.0, 0.0};
					double nmag = normVec3d(n, fa_normal[ifa]);
					vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
					double smag_half = fabs(vecDotVec3d(s_half, n));   //normVec3d(s_half);
					//double alpha = vecDotVec3d(n, s_half);           // HACK alfa should be always 1 at the boundary (provide such a grid!)

					double viscFlux = diff[ifa]*nmag*(phi[icv1]-phi[icv0])/smag_half;
					rhs_rhoScal[icv0] += viscFlux;

					if (flagImplicit)
						AScal[noc00] += diff[ifa]*nmag/smag_half/rho[icv0];
				}
			}
			// .............................................................................................
			// DIRICHLET BOUNDARY CONDITION
			// .............................................................................................
			else if (scalarZoneIsDirichlet(phiBCval, zone->getName(), scalName))
			{
				if ((mpi_rank == 0) && firstCall_JOEcalcViscousFluxScalar)
					cout << scalName << ": " << zone->getName() << " is Dirichlet: VALUE: " << phiBCval << endl;

				for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
				{
					int icv0 = cvofa[ifa][0];
					int icv1 = cvofa[ifa][1];
					int noc00 = nbocv_i[icv0];

					double n[3] = {0.0, 0.0, 0.0}, s_half[3] = {0.0, 0.0, 0.0};
					double nmag = normVec3d(n, fa_normal[ifa]);
					vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
					double smag_half = fabs(vecDotVec3d(s_half, n));   //normVec3d(s_half);
					//double alpha = vecDotVec3d(n, s_half);           // HACK alfa should be always 1 at the boundary (provide such a grid!)

					double viscFlux = diff[ifa]*nmag*(phi[icv1]-phi[icv0])/smag_half; // TODO: phi_bfa or phiBCval??
					rhs_rhoScal[icv0] += viscFlux;

					if (flagImplicit)
						AScal[noc00] += diff[ifa]*nmag/smag_half/rho[icv0];
				}
			}
			// .............................................................................................
			// FLUX BOUNDARY CONDITION: NOT IMPLEMENTED YET
			// .............................................................................................
			else if (scalarZoneIsFlux(phiBCflux, zone->getName(), scalName))
			{
				cerr << "scalarZoneIsFlux not implemented yet!" << endl;
				throw(-1);
			}
			// .............................................................................................
			// OTHER BOUNDARY CONDITIONS: NO VISCOUS FLUX
			// .............................................................................................
			else
			{
				if ((mpi_rank == 0) && firstCall_JOEcalcViscousFluxScalar)
					cout << scalName << ": " << zone->getName() << " is Symmetry" << endl;
			}
		}
	}
}



