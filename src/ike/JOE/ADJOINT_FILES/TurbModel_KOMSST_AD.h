#ifndef RANSTURBMODEL_KOMSST_AD_H
#define RANSTURBMODEL_KOMSST_AD_H

#include "turbModels/TurbModel_KOMSST.h"



/**
 * wilcox k-omega model as in Turbulence Modeling for CFD, Third Edition, Wilcox 2006
 */

class RansTurbKOmSST_AD : virtual public UgpWithCvCompFlow_AD, public RansTurbKOmSST
{
public:   // constructors
  RansTurbKOmSST_AD()
  {
    if (mpi_rank == 0)
      cout << "RansTurbKOmSST_AD()" << endl;
    
    psi_kine    = NULL;     registerScalar(psi_kine, "PSI_KINE", CV_DATA);
    psi_omega   = NULL;     registerScalar(psi_omega, "PSI_OMEGA", CV_DATA);


  }

  virtual ~RansTurbKOmSST_AD() {}

public:   // member vars

  adouble (*grad_kine)[3], *kine_diff;
  adouble *omega, (*grad_omega)[3], *omega_diff;
  adouble *limiterFunc_AD;



public:   // member functions

  
  virtual void initialize_turb_adjoint()
  {

     if (!checkDataFlag(psi_kine)) {
        if(mpi_rank==0) cout<<"INITIALIZING ADJOINT TKE TO ZERO"<<endl;
      for (int icv=0; icv<ncv_ggff; icv++)
        psi_kine[icv] = 0.;
    }

    if (!checkDataFlag(psi_omega)) {
        if(mpi_rank==0) cout<<"INITIALIZING ADJOINT OMEGA TO ZERO"<<endl;
      for (int icv=0; icv<ncv_ggff; icv++)
        psi_omega[icv] = 0.;
    }

        updateCvDataG1G2(psi_kine,  REPLACE_DATA);
        updateCvDataG1G2(psi_omega,  REPLACE_DATA);


       int nScal = scalarTranspEqVector.size();

         for(int i=0; i<nScal ; i++){
                if(strcmp(scalarTranspEqVector_AD[i].name,"kine")==0){
                  for (int icv = 0; icv < ncv_gg; icv++) {
                        scalarTranspEqVector_psi[i].phi[icv]=psi_kine[icv];
                  }
                }
                if(strcmp(scalarTranspEqVector_AD[i].name,"omega")==0){
                  for (int icv = 0; icv < ncv_gg; icv++) {
                        scalarTranspEqVector_psi[i].phi[icv]=psi_omega[icv];
                  }
                }
         }
  }

   virtual void copy_turb_adjoint()
  {   
    int nScal = scalarTranspEqVector.size();

         for(int i=0; i<nScal ; i++){
                if(strcmp(scalarTranspEqVector_AD[i].name,"kine")==0){
                  for (int icv = 0; icv < ncv_gg; icv++)
                        psi_kine[icv]=scalarTranspEqVector_psi[i].phi[icv];
                }
                if(strcmp(scalarTranspEqVector_AD[i].name,"omega")==0){
                  for (int icv = 0; icv < ncv_gg; icv++)
                        psi_omega[icv]=scalarTranspEqVector_psi[i].phi[icv];
                }
         }

  }

  virtual void initialHookScalarRansTurbModel_AD()
  {
 // connect pointers 

    int nScal = scalarTranspEqVector.size();

         for(int i=0; i<nScal ; i++){
                if(strcmp(scalarTranspEqVector_AD[i].name,"kine")==0){
                        kine      = scalarTranspEqVector_AD[i].phi ;
                        grad_kine = scalarTranspEqVector_AD[i].grad_phi ;
                        kine_diff = scalarTranspEqVector_AD[i].diff ;
                if(mpi_rank==0) cout<<"Connected Scalar Pointer "<<scalarTranspEqVector_AD[i].name<<endl;
                }
                if(strcmp(scalarTranspEqVector_AD[i].name,"omega")==0){
                        omega      = scalarTranspEqVector_AD[i].phi ;
                        grad_omega = scalarTranspEqVector_AD[i].grad_phi ;
                        omega_diff = scalarTranspEqVector_AD[i].diff ;
                if(mpi_rank==0) cout<<"Connected Scalar Pointer "<<scalarTranspEqVector_AD[i].name<<endl;
                }
         }

                  updateCvDataG1G2(wallDist, REPLACE_DATA);

    if (SST_limiterShearStress == 1)      limiterFunc_AD = strMag;   // shear stress limiter = strain rate magnitude
    else                                  limiterFunc_AD = vortMag;  // shear stress limiter = vorticity magnitude


  }


   virtual void calcMenterBlendingFunctions_AD(adouble *rho, adouble (*rhou)[3])
  {

    calcCv2Grad_AD(grad_kine, kine,   limiterNavierS, kine, epsilonSDWLS);
    calcCv2Grad_AD(grad_omega, omega, limiterNavierS, omega, epsilonSDWLS);

    for (int icv=0; icv<ncv_g; icv++)
    {
      double d = wallDist[icv];
      adouble mue = InterpolateAtCellCenterFromFaceValues(mul_fa, icv);

      crossDiff_AD[icv] = max(2.0*rho[icv]*sigma_om2/omega[icv]*vecDotVec3d_AD(grad_kine[icv], grad_omega[icv]),  1.0e-10);
//      crossDiff_AD[icv] = min(crossDiff_AD[icv], 10000.0);

      adouble gamma1 = 500.0*mue/(pow(d, 2.0)*rho[icv]*omega[icv]);
      adouble gamma2 = 4.0*sigma_om2*rho[icv]*kine[icv]/(d*d*crossDiff_AD[icv]);
      adouble gamma3 = sqrt(kine[icv])/(betaStar*omega[icv]*d);

      adouble gamma = min(max(gamma1, gamma3), gamma2);
      blendFuncF1_AD[icv] = tanh(gamma*gamma*gamma*gamma);
      gamma = max(2.0*gamma3, gamma1);
      blendFuncF2_AD[icv] = tanh(gamma*gamma);
    }

       double (*buf) = new double[ncv_ggff];

    FOR_ICV
        buf[icv] = crossDiff_AD[icv].value();
        updateCvDataG1G2(buf, REPLACE_DATA);
    FOR_ICV_G2_ONLY
        crossDiff_AD[icv] = buf[icv];

    FOR_ICV
        buf[icv] = blendFuncF1_AD[icv].value();
        updateCvDataG1G2(buf, REPLACE_DATA);
    FOR_ICV_G2_ONLY
        blendFuncF1_AD[icv] = buf[icv];

    FOR_ICV
        buf[icv] = blendFuncF2_AD[icv].value();
        updateCvDataG1G2(buf, REPLACE_DATA);
    FOR_ICV_G2_ONLY
        blendFuncF2_AD[icv] = buf[icv];

   delete [] buf;



    
  }



  virtual void calcRansTurbViscMuet_AD(adouble *rho, adouble (*rhou)[3])
  {
    // update strain rate tensor 
    calcStrainRateAndDivergence_AD();
    calcVorticity_AD();

    calcMenterBlendingFunctions_AD(rho, rhou);

    // internal faces
       //for (int ifa = nfa_b; ifa < nfa; ifa++)
  for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) // Ensure we sweep through internal faces of ghosts too
  {
   if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg))
   {

      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];

      double dx0[3], dx1[3];
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      double w0 = sqrt(vecDotVec3d(dx0, dx0));
      double w1 = sqrt(vecDotVec3d(dx1, dx1));

      adouble rho_fa = (w1*rho[icv0] + w0*rho[icv1])/(w0+w1);
      adouble kine_fa = (w1*kine[icv0] + w0*kine[icv1])/(w0+w1);
      adouble om_fa = (w1*omega[icv0] + w0*omega[icv1])/(w0+w1);
      adouble limitFunc_fa = (w1*limiterFunc_AD[icv0] + w0*limiterFunc_AD[icv1])/(w0+w1);
      adouble blendFuncF2_AD_fa = (w1*blendFuncF2_AD[icv0] + w0*blendFuncF2_AD[icv1])/(w0+w1);

      adouble zeta = min(1.0/om_fa, a1/(limitFunc_fa*blendFuncF2_AD_fa+1.e-20));

      mut_fa[ifa] = min(max(rho_fa*kine_fa*zeta, 0.0), 1.0);
    }
   }

    // boundary faces
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      if (zoneIsWall(zone->getName()))
	for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
          mut_fa[ifa] = 0.0;                                     // set mut zero at walls
	}
      else
	for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
          int icv0 = cvofa[ifa][0];
          int icv1 = cvofa[ifa][1];
          
          adouble zeta = min(1.0/omega[icv1], a1/(limiterFunc_AD[icv0]*blendFuncF2_AD[icv0]+1.e-20));
          mut_fa[ifa] = min(max(rho[icv1]*kine[icv1]*zeta, 0.0), 1.0); 
        }
    }


/*
        // TBD_AD move this elsewhere later
	delete [] blendFuncF1_AD;
	delete [] blendFuncF2_AD;
	delete [] crossDiff_AD;
*/

  }



  virtual void diffusivityHookScalarRansTurb_AD(const string &name)
  {

    if ((name == "kine") || (name == "omega"))
    {
      //ScalarTranspEq *eq = getScalarTransportData(name);
      double sigma_1,sigma_2;

      if (name == "kine")  { sigma_1 = sigma_k1;  sigma_2 = sigma_k2; }
      if (name == "omega") { sigma_1 = sigma_om1; sigma_2 = sigma_om2;}

      // internal faces
	  //for (int ifa = nfa_b; ifa < nfa; ifa++)
  for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) // Ensure we sweep through internal faces of ghosts too
  {
   if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg))
   {

        int icv0 = cvofa[ifa][0];
        int icv1 = cvofa[ifa][1];

        double dx0[3], dx1[3];
        vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
        vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
        double w0 = sqrt(vecDotVec3d(dx0, dx0));
        double w1 = sqrt(vecDotVec3d(dx1, dx1));

        adouble f1_fa = (w1*blendFuncF1_AD[icv0] + w0*blendFuncF1_AD[icv1])/(w0+w1);
        adouble coeff = sigma_1*f1_fa + (1. - f1_fa)*sigma_2;


        if(name=="kine")   kine_diff[ifa] = mul_fa[ifa] + coeff*mut_fa[ifa];
        if(name=="omega") omega_diff[ifa] = mul_fa[ifa] + coeff*mut_fa[ifa];
      }
   }
      // boundary faces
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))
        {
	for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
            int icv0 = cvofa[ifa][0];
            if(name=="kine")   kine_diff[ifa] = mul_fa[ifa];
            if(name=="omega")  omega_diff[ifa] = mul_fa[ifa];
          }
        }
        else
        {
	for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
            int icv0 = cvofa[ifa][0];
            adouble coeff = sigma_1*blendFuncF1_AD[icv0] + (1. - blendFuncF1_AD[icv0])*sigma_2;

            if(name=="kine")   kine_diff[ifa] = mul_fa[ifa] + coeff*mut_fa[ifa]; 
            if(name=="omega") omega_diff[ifa] = mul_fa[ifa] + coeff*mut_fa[ifa]; 
          }
        }
      }
    }
  }

     virtual adouble calcTurbProd_AD(int icv, int limiter)
  {
//    adouble ome = max(omega[icv], 1.0*strMag[icv]/sqrt(betaStar));
    adouble ome = omega[icv];

    adouble zeta = min(1.0/ome, a1/(limiterFunc_AD[icv]*blendFuncF2_AD[icv]+1.e-20));
    adouble mut = min(max(rho_AD[icv]*kine[icv]*zeta, 0.0), 1.0e5);
    adouble Pk = mut*strMag[icv]*strMag[icv] - 2./3.*rho_AD[icv]*kine[icv]*diverg[icv];

    if (limiter == 1)      return min(Pk, 10.0*betaStar*rho_AD[icv]*kine[icv]*omega[icv]);
    else                   return Pk;
  }



    virtual void sourceHookScalarRansTurb_new_AD(adouble *rhs, double *A, const string &name, int flagImplicit)
  {
    if (name == "kine")
      for (int icv=0; icv<ncv_g; icv++)
      {
        adouble src = calcTurbProd_AD(icv, SST_limitPK) - betaStar*rho_AD[icv]*omega[icv]*kine[icv];
        rhs[icv] += src*cv_volume[icv];

        if (flagImplicit && icv<ncv)
        {
          int noc00 = nbocv_i[icv];
          adouble dsrcdphi = - betaStar*omega[icv];
          A[noc00] -= dsrcdphi.value()*cv_volume[icv];
        }
      }

    if (name == "omega")
    {
      for (int icv=0; icv<ncv_g; icv++)
      {
        adouble F1 = blendFuncF1_AD[icv];
        adouble alfa = F1*alfa_1 + (1.0 - F1)*alfa_2;
        adouble beta = F1*beta_1 + (1.0 - F1)*beta_2;

        adouble zeta = min(1.0/omega[icv], a1/(limiterFunc_AD[icv]*blendFuncF2_AD[icv]+1.e-20));
        adouble mut = min(max(rho_AD[icv]*kine[icv]*zeta, 1.e-08), 1.0e5);

        adouble src = alfa*rho_AD[icv]/mut*calcTurbProd_AD(icv,SST_limitPK) + (1.0-F1)*crossDiff_AD[icv] - beta*rho_AD[icv]*omega[icv]*omega[icv];
        rhs[icv] += src*cv_volume[icv];

        if (flagImplicit && icv<ncv)
        {
          int noc00 = nbocv_i[icv];
          adouble dsrcdphi =  - 2.0*beta*omega[icv];
          A[noc00] -= dsrcdphi.value()*cv_volume[icv];
        }
      }
    }
  }



    virtual void sourceHookRansTurbCoupled_AD(adouble **rhs, double ***A, int flagImplicit)
  {
    int kine_Index = getScalarTransportIndex("kine");
    int omega_Index = getScalarTransportIndex("omega");

      for (int icv=0; icv<ncv_g; icv++)
      {
        adouble src = calcTurbProd_AD(icv, SST_limitPK) - betaStar*rho_AD[icv]*omega[icv]*kine[icv];
        rhs[icv][5+kine_Index] += src*cv_volume[icv];

        if (flagImplicit && icv<ncv)
        {
          int noc00 = nbocv_i[icv];
          adouble dsrcdphi = - betaStar*omega[icv];
          A[noc00][5+kine_Index][5+kine_Index] -= dsrcdphi.value()*cv_volume[icv];
        }
      }

      for (int icv=0; icv<ncv_g; icv++)
      {
        adouble F1 = blendFuncF1_AD[icv];
        adouble alfa = F1*alfa_1 + (1.0 - F1)*alfa_2;
        adouble beta = F1*beta_1 + (1.0 - F1)*beta_2;

        adouble zeta = min(1.0/omega[icv], a1/(limiterFunc_AD[icv]*blendFuncF2_AD[icv]+1.e-20));
        adouble mut = min(max(rho_AD[icv]*kine[icv]*zeta, 1.e-08), 1.0e5);

        adouble src = alfa*rho_AD[icv]/mut*calcTurbProd_AD(icv,SST_limitPK) + (1.0-F1)*crossDiff_AD[icv] - beta*rho_AD[icv]*omega[icv]*omega[icv];
        rhs[icv][5+omega_Index] += src*cv_volume[icv];

        if (flagImplicit && icv<ncv)
        {
          int noc00 = nbocv_i[icv];
          adouble dsrcdphi =  - 2.0*beta*omega[icv];
          A[noc00][5+omega_Index][5+omega_Index] -= dsrcdphi.value()*cv_volume[icv];
        }
      }

  }


  virtual void boundaryHookScalarRansTurb_AD(adouble *phi_fa, FaZone *zone, const string &name)
  {
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      Param *param;
      if (getParam(param, zone->getName()))
        if ((param->getString() == "WALL") && (name == "omega"))
        {
          //for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          //{
          for (int index = 0; index < zone->faVec.size(); ++index) {
            int ifa = zone->faVec[index];
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];

            double nVec[3], s_half[3];
            normVec3d(nVec, fa_normal[ifa]);
            vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
            double wallDist = fabs(vecDotVec3d(s_half, nVec));
            phi_fa[icv1] = 60.0*calcMuLam_AD(icv0)/(rho_AD[icv0]*beta_1*wallDist*wallDist);

          }
        }
    }
  }

};



#endif
