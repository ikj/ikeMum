#ifndef RANSTURBMODEL_KOM_AD_H
#define RANSTURBMODEL_KOM_AD_H

#include "../turbModels/TurbModel_KOM.h"



/**
 * wilcox k-omega model as in Turbulence Modeling for CFD, Third Edition, Wilcox 2006
 */

class RansTurbKOm_AD : virtual public UgpWithCvCompFlow_AD, public RansTurbKOm
{
public:   // constructors
  RansTurbKOm_AD()
  {
    if (mpi_rank == 0)
      cout << "RansTurbKOm_AD()" << endl;
    
    psi_kine    = NULL;     registerScalar(psi_kine,  "PSI_KINE",  CV_DATA);
    psi_omega   = NULL;     registerScalar(psi_omega, "PSI_OMEGA", CV_DATA);



  }

  virtual ~RansTurbKOm_AD() {}

public:   // member vars

  adouble  (*grad_kine)[3], *kine_diff;
  adouble *omega, (*grad_omega)[3], *omega_diff;
  
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
       int nScal = scalarTranspEqVector.size();

        updateCvDataG1G2(psi_kine,  REPLACE_DATA);
        updateCvDataG1G2(psi_omega, REPLACE_DATA);

         for(int i=0; i<nScal ; i++){
                if(strcmp(scalarTranspEqVector_AD[i].name,"kine")==0){
		if(mpi_rank==0) cout<<"Initializing  KINE"<<endl;
                  for (int icv = 0; icv < ncv_gg; icv++) {
                        scalarTranspEqVector_psi[i].phi[icv]=psi_kine[icv];
                  }
                }
                if(strcmp(scalarTranspEqVector_AD[i].name,"omega")==0){
		if(mpi_rank==0) cout<<"Initializing  OMEGA"<<endl;
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
  }


  virtual void calcRansTurbViscMuet_AD(adouble *rho, adouble (*rhou)[3])
  {
    // update strain rate tensor 
    calcStrainRateAndDivergence_AD();
    calcVorticity_AD();

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
      double ws = w0 + w1;
      w0 /= ws;
      w1 /= ws;

      adouble rho_fa    = w1*rho[icv0] + w0*rho[icv1];
      adouble kine_fa   = w1*kine[icv0] + w0*kine[icv1];
      adouble om_fa     = w1*omega[icv0] + w0*omega[icv1];
      adouble strMag_fa = w1*strMag[icv0] + w0*strMag[icv1];

      if (KOM_RealizableConstraint == 1)
      {
        adouble omega_tilde = max(om_fa, cLim*strMag_fa/sqrt(betaStar));
        adouble mut_f = min(rho_fa*kine_fa/omega_tilde, 100.0);
        mut_fa[ifa] = mut_f.value();
      }
      else if (KOM_RealizableConstraint == 2)
      {
        adouble TS = min(1.0/om_fa, 0.6/(sqrt(6.0)*strMag_fa));
        adouble mut_f = rho_fa*kine_fa*TS;
        mut_fa[ifa] = mut_f.value();
      }
      else
      {
        adouble mut_f = min(rho_fa*kine_fa/om_fa, 100.0);
        mut_fa[ifa] = mut_f.value();
      }
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
          
          if (KOM_RealizableConstraint == 1)
          {
            adouble omega_tilde = max(omega[icv1], cLim*strMag[icv0]/sqrt(betaStar));
            adouble mut_f = min(rho[icv1]*kine[icv1]/omega_tilde, 100.0);    // zero order extrapolation for others
            mut_fa[ifa] = mut_f.value();
          }
          else if (KOM_RealizableConstraint == 2)
          {
            adouble TS = min(1.0/omega[icv1], 0.6/(sqrt(6.0)*strMag[icv0]));
            adouble mut_f = rho[icv1]*kine[icv1]*TS;
            mut_fa[ifa] = mut_f.value();
          }
          else {
            adouble mut_f = min(rho[icv1]*kine[icv1]/omega[icv1], 100.0);    // zero order extrapolation for others
            mut_fa[ifa] = mut_f.value();
          }
        }
    }

  }


  virtual void diffusivityHookScalarRansTurb_AD(const string &name)
  {

    if ((name == "kine") || (name == "omega"))
    {
      double sigma;
      if (name == "kine")       sigma = sigmaStar;
      if (name == "omega")      sigma = sigmaOmega;

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

        adouble rho_fa = (w1*rho_AD[icv0] + w0*rho_AD[icv1])/(w0+w1);
        adouble kine_fa = (w1*kine[icv0] + w0*kine[icv1])/(w0+w1);
        adouble om_fa = (w1*omega[icv0] + w0*omega[icv1])/(w0+w1);

        if(name=="kine")   kine_diff[ifa] = mul_fa[ifa] + sigma*rho_fa*kine_fa/om_fa;
        if(name=="omega") omega_diff[ifa] = mul_fa[ifa] + sigma*rho_fa*kine_fa/om_fa;
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
            int icv1 = cvofa[ifa][1];
            if(name=="kine") kine_diff[ifa]   = mul_fa[ifa] + sigma*rho_AD[icv1]*kine[icv1]/omega[icv1]; 
            if(name=="omega") omega_diff[ifa] = mul_fa[ifa] + sigma*rho_AD[icv1]*kine[icv1]/omega[icv1]; 
          }
        }
      }
    }
  }

	
  
 

  virtual adouble calcTurbProd_AD(int icv)
  {
    if (KOM_RealizableConstraint == 1)
    {
      adouble omega_tilde = max(omega[icv], cLim*strMag[icv]/sqrt(betaStar));
      adouble mu_t = rho_AD[icv]*kine[icv]/omega_tilde;
      return max(mu_t*strMag[icv]*strMag[icv] - 2./3.*rho_AD[icv]*kine[icv]*diverg[icv], 0.0);
    }
    else if (KOM_RealizableConstraint == 2)
    {
      adouble TS = min(1.0/(omega[icv]), 0.6/(sqrt(6.0)*strMag[icv]));
      adouble mu_t = rho_AD[icv]*kine[icv]*TS;
      return max(mu_t*strMag[icv]*strMag[icv] - 2./3.*rho_AD[icv]*kine[icv]*diverg[icv], 0.0);
    }
    else
    {
      adouble mu_t = rho_AD[icv]*kine[icv]/omega[icv];
      return max(mu_t*strMag[icv]*strMag[icv] - 2./3.*rho_AD[icv]*kine[icv]*diverg[icv], 0.0);
    }
  }


  virtual void sourceHookScalarRansTurb_new_AD(adouble *rhs, double *A, const string &name, int flagImplicit)
  {
    if (name == "kine")
    for (int icv=0; icv<ncv_g; icv++)
    {
      adouble src = calcTurbProd_AD(icv) - betaStar*rho_AD[icv]*omega[icv]*kine[icv];
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
      adouble OM[3][3], STR_hat[3][3];

      for (int icv=0; icv<ncv_g; icv++)
      {

        adouble TS = 1.0;
        if (KOM_RealizableConstraint == 2)
          TS = min(1.0/(omega[icv]), 0.6/(sqrt(6.0)*strMag[icv]));

        adouble sigmad;
        adouble crossDiff = vecDotVec3d_AD(grad_kine[icv], grad_omega[icv]);

        if (crossDiff <= 0.0)   sigmad = 0.0;
        else                    sigmad = sigmad0;

        adouble divU = grad_u[icv][0][0] + grad_u[icv][1][1] + grad_u[icv][2][2];

        for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
        {
          OM[i][j] = 0.5*(grad_u[icv][i][j] - grad_u[icv][j][i]);

          // strange: STR_hat subtracts 0.5 divergence instead 1/3, no bug!!!
          if (i==j)   STR_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i] - divU);
          else        STR_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]);
        }

        adouble chiOm = 0.0;
        for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
        for (int k=0; k<3; k++)
          chiOm += OM[i][j]*OM[j][k]*STR_hat[k][i];
        chiOm = fabs(chiOm/(betaStar*omega[icv]*betaStar*omega[icv]*betaStar*omega[icv]));

        adouble fbeta = (1.0 + 85.0*chiOm)/(1.0 + 100.0*chiOm);
        adouble beta = beta0*fbeta;

        adouble src =  TS*alfa*omega[icv]/kine[icv]*calcTurbProd_AD(icv)
                    - TS*beta*rho_AD[icv]*omega[icv]*omega[icv]
                    + sigmad*rho_AD[icv]/omega[icv]*crossDiff;

        rhs[icv] += src*cv_volume[icv];

        if (flagImplicit && icv<ncv)
        {
          int noc00 = nbocv_i[icv];
          adouble dsrcdphi = -2.0*TS*beta*omega[icv];
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
      adouble src = calcTurbProd_AD(icv) - betaStar*rho_AD[icv]*omega[icv]*kine[icv];
      rhs[icv][5+kine_Index] += src*cv_volume[icv];

      if (flagImplicit && icv<ncv)
      {
        int noc00 = nbocv_i[icv];
        adouble dsrcdphi = - betaStar*omega[icv];
        A[noc00][5+kine_Index][5+kine_Index] -= dsrcdphi.value()*cv_volume[icv];
      }
    }

      adouble OM[3][3], STR_hat[3][3];

      for (int icv=0; icv<ncv_g; icv++)
      {

        adouble TS = 1.0;
        if (KOM_RealizableConstraint == 2)
          TS = min(1.0/(omega[icv]), 0.6/(sqrt(6.0)*strMag[icv]));

        adouble sigmad;
        adouble crossDiff = vecDotVec3d_AD(grad_kine[icv], grad_omega[icv]);

        if (crossDiff <= 0.0)   sigmad = 0.0;
        else                    sigmad = sigmad0;

        adouble divU = grad_u[icv][0][0] + grad_u[icv][1][1] + grad_u[icv][2][2];

        for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
        {
          OM[i][j] = 0.5*(grad_u[icv][i][j] - grad_u[icv][j][i]);

          // strange: STR_hat subtracts 0.5 divergence instead 1/3, no bug!!!
          if (i==j)   STR_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i] - divU);
          else        STR_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]);
        }

        adouble chiOm = 0.0;
        for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
        for (int k=0; k<3; k++)
          chiOm += OM[i][j]*OM[j][k]*STR_hat[k][i];
        chiOm = fabs(chiOm/(betaStar*omega[icv]*betaStar*omega[icv]*betaStar*omega[icv]));

        adouble fbeta = (1.0 + 85.0*chiOm)/(1.0 + 100.0*chiOm);
        adouble beta = beta0*fbeta;

        adouble src =  TS*alfa*omega[icv]/kine[icv]*calcTurbProd_AD(icv)
                    - TS*beta*rho_AD[icv]*omega[icv]*omega[icv]
                    + sigmad*rho_AD[icv]/omega[icv]*crossDiff;

        rhs[icv][5+omega_Index] += src*cv_volume[icv];


        if (flagImplicit && icv<ncv)
        {
          int noc00 = nbocv_i[icv];
          adouble dsrcdphi = -2.0*TS*beta*omega[icv];
          A[noc00][5+omega_Index][5+omega_Index] -= dsrcdphi.value()*cv_volume[icv];
        }
      }


  }


  virtual void boundaryHookScalarRansTurb_AD(adouble *phi, FaZone *zone, const string &name)
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
            phi[icv1] = 6.0*calcMuLam_AD(icv0)/(rho_AD[icv0]*beta0*wallDist[icv0]*wallDist[icv0]);
          }
        }
    }
  }

};



#endif
