#include "UgpWithCvCompFlowAD.h"

#ifdef USE_MEM_SAVING_ADVAR
void UgpWithCvCompFlow_AD::setScalarBC_AD(FaZone *zone)
{
  int nScal = scalarTranspEqVector.size();

  for (int scal = 0; scal < nScal; scal++)
  {
    string scalName = scalarTranspEqVector[scal].getName();
    double phiBCval = 0.0, phiBCflux = 0.0;

    ADscalar<adouble> *phi = &(scalarTranspEqVector_AD[scal].phi);

    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      // .............................................................................................
      // HOOK BOUNDARY CONDITION
      // .............................................................................................
      if (scalarZoneIsHook(zone->getName(), scalName))
      {
        // user defined boundary hook
        boundaryHookScalarRansTurb_AD((*phi), &(*zone), scalName);
        boundaryHookScalarRansComb_AD((*phi), &(*zone), scalName);

      }
      // .............................................................................................
      // DIRICHLET BOUNDARY CONDITION
      // .............................................................................................
      else if (scalarZoneIsDirichlet(phiBCval, zone->getName(), scalName))
      {
    	  for (int index = 0; index < zone->faVec.size(); ++index) {
    		  int ifa = zone->faVec[index];
    		  int icv0 = cvofa[ifa][0];
    		  int icv1 = cvofa[ifa][1];

    		  if (vecDotVec3d_AD(vel[icv1], fa_normal[ifa]) > 1.0e-8)
    			  (*phi)[icv1] = (*phi)[icv0];
    		  else
    			  (*phi)[icv1] = phiBCval;
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
    	  for (int index = 0; index < zone->faVec.size(); ++index) {
    		  int ifa = zone->faVec[index];
    		  int icv0 = cvofa[ifa][0];
    		  int icv1 = cvofa[ifa][1];

    		  (*phi)[icv1] = (*phi)[icv0];
    	  }
      }
    }
  }
}
#else
void UgpWithCvCompFlow_AD::setScalarBC_AD(FaZone *zone)
{
  int nScal = scalarTranspEqVector.size();

  for (int scal = 0; scal < nScal; scal++)
  {
    string scalName = scalarTranspEqVector[scal].getName();
    double phiBCval = 0.0, phiBCflux = 0.0;

    adouble *phi = scalarTranspEqVector_AD[scal].phi;

    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      // .............................................................................................
      // HOOK BOUNDARY CONDITION
      // .............................................................................................
      if (scalarZoneIsHook(zone->getName(), scalName))
      {
        // user defined boundary hook
        boundaryHookScalarRansTurb_AD(phi, &(*zone), scalName);
        boundaryHookScalarRansComb_AD(phi, &(*zone), scalName);

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

          if (vecDotVec3d_AD(vel[icv1], fa_normal[ifa]) > 1.0e-8)
        	  phi[icv1] = phi[icv0];
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

          phi[icv1] = phi[icv0];
        }
      }
    }
  }
}
#endif

#ifdef USE_MEM_SAVING_ADVAR
void UgpWithCvCompFlow_AD::calcViscousFluxScalar_new_AD(adouble *rhs_rhoScal, double *AScal, ScalarTranspEq &transpScal, ScalarTranspEq_AD &transpScal_AD, int flagImplicit)
{
  if (!transpScal.diffTerm)
    return;

  string scalName(transpScal.getName());
  double phiBCval = 0.0, phiBCflux = 0.0;

  ADscalar<adouble> (*phi) = &(transpScal_AD.phi);
  ADvector<adouble> (*grad_phi) = &(transpScal_AD.grad_phi);
  ADscalar<adouble> (*diff) = &(transpScal_AD.diff);

  // =============================================================================================
  // compute gradients, with boundary values
  // =============================================================================================
  calcCv2Grad_AD((*grad_phi), (*phi), limiterNavierS, (*phi), epsilonSDWLS);

  // =============================================================================================
  // compute user-defined diffusivity
  // =============================================================================================
  diffusivityHookScalarRansTurb_AD(scalName);
  diffusivityHookScalarRansComb_AD(scalName);

  // =============================================================================================
  // cycle trough internal faces and compute viscous flux and implicit matrix
  // =============================================================================================
  //for (int ifa = nfa_b; ifa < nfa; ifa++)
  for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) // Ensure we sweep through internal faces of ghosts too
  {
   if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg))
   {

    int icv0 = cvofa[ifa][0];
    int icv1 = cvofa[ifa][1];
    int noc00, noc01, noc11, noc10;

     if (flagImplicit && icv0 < ncv && icv1 < ncv_g)
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

    adouble viscFlux = (*diff)[ifa] * nmag * (alpha * ((*phi)[icv1] - (*phi)[icv0]) / smag
                    + 0.5 * ((*grad_phi)[icv0][0] + (*grad_phi)[icv1][0]) * (n[0] - alpha * s[0])
                    + 0.5 * ((*grad_phi)[icv0][1] + (*grad_phi)[icv1][1]) * (n[1] - alpha * s[1])
                    + 0.5 * ((*grad_phi)[icv0][2] + (*grad_phi)[icv1][2]) * (n[2] - alpha * s[2]));

    rhs_rhoScal[icv0] += viscFlux;

    if (icv1 < ncv_gg)
      rhs_rhoScal[icv1] -= viscFlux;

    if (flagImplicit && icv0 < ncv && icv1 < ncv_g)
    {
      AScal[noc00] +=   (*diff)[ifa].value() * nmag * alpha / smag / rho_AD[icv0].value();
      AScal[noc01] += - (*diff)[ifa].value() * nmag * alpha / smag / rho_AD[icv0].value();

      if (icv1 < ncv)
      {
        AScal[noc11] +=   (*diff)[ifa].value() * nmag * alpha / smag / rho_AD[icv1].value();
        AScal[noc10] += - (*diff)[ifa].value() * nmag * alpha / smag / rho_AD[icv1].value();
      }
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
        if ((step == 1) && (mpi_rank == 0))
          cout << scalName << ": " << zone->getName() << " is HOOK" << endl;

        for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
          int icv0 = cvofa[ifa][0];
          int icv1 = cvofa[ifa][1];
          int noc00 = nbocv_i[icv0];

          double n[3] = {0.0, 0.0, 0.0}, s_half[3] = {0.0, 0.0, 0.0};
          double nmag = normVec3d(n, fa_normal[ifa]);
          vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
          double smag_half = fabs(vecDotVec3d(s_half, n));   //normVec3d(s_half);
          //double alpha = vecDotVec3d(n, s_half);           // HACK alfa should be always 1 at the boundary (provide such a grid!)

          adouble viscFlux = (*diff)[ifa] * nmag * ((*phi)[icv1] - (*phi)[icv0]) / smag_half;
          rhs_rhoScal[icv0] += viscFlux;

          if (flagImplicit && icv0<ncv)
            AScal[noc00] += (*diff)[ifa].value() * nmag / smag_half / rho_AD[icv0].value();
        }
      }
      // .............................................................................................
      // DIRICHLET BOUNDARY CONDITION
      // .............................................................................................
      else if (scalarZoneIsDirichlet(phiBCval, zone->getName(), scalName))
      {
        if ((step == 1) && (mpi_rank == 0))
          cout << scalName << ": " << zone->getName() << " is Dirichlet: VALUE: " << phiBCval << endl;

        for (int index = 0; index < zone->faVec.size(); ++index) {
        int ifa = zone->faVec[index];
          int icv0 = cvofa[ifa][0];
          int icv1 = cvofa[ifa][1];
          int noc00 = nbocv_i[icv0];

          double n[3] = {0.0, 0.0, 0.0}, s_half[3] = {0.0, 0.0, 0.0};
          double nmag = normVec3d(n, fa_normal[ifa]);
          vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
          double smag_half = fabs(vecDotVec3d(s_half, n));   //normVec3d(s_half);
          //double alpha = vecDotVec3d(n, s_half);           // HACK alfa should be always 1 at the boundary (provide such a grid!)

          adouble viscFlux = (*diff)[ifa] * nmag * ((*phi)[icv1] - (*phi)[icv0]) / smag_half;                // TODO: phi_bfa or phiBCval??
          rhs_rhoScal[icv0] += viscFlux;

          if (flagImplicit && icv0<ncv)
            AScal[noc00] += (*diff)[ifa].value() * nmag / smag_half / rho_AD[icv0].value();
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
        if ((step == 1) && (mpi_rank == 0))
          cout << scalName << ": " << zone->getName() << " is Symmetry" << endl;
      }
    }
  }
}
#else
void UgpWithCvCompFlow_AD::calcViscousFluxScalar_new_AD(adouble *rhs_rhoScal, double *AScal, ScalarTranspEq &transpScal, ScalarTranspEq_AD &transpScal_AD, int flagImplicit)
{
  if (!transpScal.diffTerm)
    return;

  string scalName(transpScal.getName());
  double phiBCval = 0.0, phiBCflux = 0.0;

  adouble *phi = transpScal_AD.phi;
  adouble (*grad_phi)[3] = transpScal_AD.grad_phi;
  adouble *diff = transpScal_AD.diff;

  // =============================================================================================
  // compute gradients, with boundary values
  // =============================================================================================
  calcCv2Grad_AD(grad_phi, phi, limiterNavierS, phi, epsilonSDWLS);

  // =============================================================================================
  // compute user-defined diffusivity
  // =============================================================================================
  diffusivityHookScalarRansTurb_AD(scalName);
  diffusivityHookScalarRansComb_AD(scalName);

  // =============================================================================================
  // cycle trough internal faces and compute viscous flux and implicit matrix
  // =============================================================================================
  //for (int ifa = nfa_b; ifa < nfa; ifa++)
  for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) // Ensure we sweep through internal faces of ghosts too
  {
   if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg))
   {

    int icv0 = cvofa[ifa][0];
    int icv1 = cvofa[ifa][1];
    int noc00, noc01, noc11, noc10;

     if (flagImplicit && icv0 < ncv && icv1 < ncv_g)
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

    adouble viscFlux = diff[ifa] * nmag * (alpha * (phi[icv1] - phi[icv0]) / smag
                    + 0.5 * (grad_phi[icv0][0] + grad_phi[icv1][0]) * (n[0] - alpha * s[0])
                    + 0.5 * (grad_phi[icv0][1] + grad_phi[icv1][1]) * (n[1] - alpha * s[1])
                    + 0.5 * (grad_phi[icv0][2] + grad_phi[icv1][2]) * (n[2] - alpha * s[2]));

    rhs_rhoScal[icv0] += viscFlux;

    if (icv1 < ncv_gg)
      rhs_rhoScal[icv1] -= viscFlux;

    if (flagImplicit && icv0 < ncv && icv1 < ncv_g)
    {
      AScal[noc00] +=   diff[ifa].value() * nmag * alpha / smag / rho_AD[icv0].value();
      AScal[noc01] += - diff[ifa].value() * nmag * alpha / smag / rho_AD[icv0].value();

      if (icv1 < ncv)
      {
        AScal[noc11] +=   diff[ifa].value() * nmag * alpha / smag / rho_AD[icv1].value();
        AScal[noc10] += - diff[ifa].value() * nmag * alpha / smag / rho_AD[icv1].value();
      }
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
        if ((step == 1) && (mpi_rank == 0))
          cout << scalName << ": " << zone->getName() << " is HOOK" << endl;

        for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
          int icv0 = cvofa[ifa][0];
          int icv1 = cvofa[ifa][1];
          int noc00 = nbocv_i[icv0];

          double n[3] = {0.0, 0.0, 0.0}, s_half[3] = {0.0, 0.0, 0.0};
          double nmag = normVec3d(n, fa_normal[ifa]);
          vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
          double smag_half = fabs(vecDotVec3d(s_half, n));   //normVec3d(s_half);
          //double alpha = vecDotVec3d(n, s_half);           // HACK alfa should be always 1 at the boundary (provide such a grid!)

          adouble viscFlux = diff[ifa] * nmag * (phi[icv1] - phi[icv0]) / smag_half;
          rhs_rhoScal[icv0] += viscFlux;

          if (flagImplicit && icv0<ncv)
            AScal[noc00] += diff[ifa].value() * nmag / smag_half / rho_AD[icv0].value();
        }
      }
      // .............................................................................................
      // DIRICHLET BOUNDARY CONDITION
      // .............................................................................................
      else if (scalarZoneIsDirichlet(phiBCval, zone->getName(), scalName))
      {
        if ((step == 1) && (mpi_rank == 0))
          cout << scalName << ": " << zone->getName() << " is Dirichlet: VALUE: " << phiBCval << endl;

        for (int index = 0; index < zone->faVec.size(); ++index) {
        int ifa = zone->faVec[index];
          int icv0 = cvofa[ifa][0];
          int icv1 = cvofa[ifa][1];
          int noc00 = nbocv_i[icv0];

          double n[3] = {0.0, 0.0, 0.0}, s_half[3] = {0.0, 0.0, 0.0};
          double nmag = normVec3d(n, fa_normal[ifa]);
          vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
          double smag_half = fabs(vecDotVec3d(s_half, n));   //normVec3d(s_half);
          //double alpha = vecDotVec3d(n, s_half);           // HACK alfa should be always 1 at the boundary (provide such a grid!)

          adouble viscFlux = diff[ifa] * nmag * (phi[icv1] - phi[icv0]) / smag_half;                // TODO: phi_bfa or phiBCval??
          rhs_rhoScal[icv0] += viscFlux;

          if (flagImplicit && icv0<ncv)
            AScal[noc00] += diff[ifa].value() * nmag / smag_half / rho_AD[icv0].value();
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
        if ((step == 1) && (mpi_rank == 0))
          cout << scalName << ": " << zone->getName() << " is Symmetry" << endl;
      }
    }
  }

}
#endif


