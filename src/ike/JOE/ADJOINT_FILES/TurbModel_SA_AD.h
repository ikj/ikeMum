#ifndef RANSTURBMODEL_SA_AD_H
#define RANSTURBMODEL_SA_AD_H

#include "../turbModels/TurbModel_SA.h"




class RansTurbSA_AD : virtual public UgpWithCvCompFlow_AD, public RansTurbSA
{
public:   // constructors
  RansTurbSA_AD()
  {
    if (mpi_rank == 0)
      cout << "RansTurbSA_AD()" << endl;

    psi_nuSA   = NULL;     registerScalar(psi_nuSA, "PSI_NUSA", CV_DATA);


  }

  virtual ~RansTurbSA_AD() {}

public:   // member vars

   adouble *nuSA;
   adouble *nuSA_diff;
   adouble (*grad_nuSA)[3];

public:   // member functions

  virtual void initialize_turb_adjoint() 
  {
  if (!checkDataFlag(psi_nuSA)) {
        if(mpi_rank==0) cout<<"INITIALIZING ADJOINT TKE TO ZERO"<<endl;
      for (int icv=0; icv<ncv_ggff; icv++)
        psi_nuSA[icv] = 0.;
    }

    updateCvDataG1G2(psi_nuSA,  REPLACE_DATA);

    int nScal = scalarTranspEqVector.size();

         for(int i=0; i<nScal ; i++){
		if(strcmp(scalarTranspEqVector_AD[i].name,"nuSA")==0){
		  for (int icv = 0; icv < ncv_gg; icv++) {
			scalarTranspEqVector_psi[i].phi[icv]=psi_nuSA[icv];
		  }
                }
         }

  }

  virtual void copy_turb_adjoint() 
  {
    int nScal = scalarTranspEqVector.size();

         for(int i=0; i<nScal ; i++){
		if(strcmp(scalarTranspEqVector_AD[i].name,"nuSA")==0){
		  for (int icv = 0; icv < ncv_gg; icv++)
			psi_nuSA[icv]=scalarTranspEqVector_psi[i].phi[icv];
                }
         }

  }

  virtual void initialHookScalarRansTurbModel_AD()
  {
 // connect pointers 

    int nScal = scalarTranspEqVector.size();

         for(int i=0; i<nScal ; i++){
		if(strcmp(scalarTranspEqVector_AD[i].name,"nuSA")==0){
			nuSA     = scalarTranspEqVector_AD[i].phi ;
			nuSA_diff = scalarTranspEqVector_AD[i].diff ;
			grad_nuSA = scalarTranspEqVector_AD[i].grad_phi ;
	 	if(mpi_rank==0) cout<<"Connected Scalar Pointer "<<scalarTranspEqVector_AD[i].name<<endl;
                }
         }

        updateCvDataG1G2(wallDist, REPLACE_DATA);
  }

  virtual void calcRansTurbViscMuet_AD(adouble *rho, adouble (*rhou)[3])
  {


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

      adouble rho_fa  = (w1*rho[icv0]   + w0*rho[icv1])/(w0+w1);
      adouble nuSA_fa = (w1*nuSA[icv0]  + w0*nuSA[icv1])/(w0+w1);

      adouble muSA_fa = rho_fa*nuSA_fa;
      adouble chi_3   = pow(muSA_fa/mul_fa[ifa], 3.0);
      adouble fv1     = chi_3/(chi_3 + cv1_3);
      
      mut_fa[ifa]  = max(muSA_fa*fv1, 0.0);              
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
          int icv1 = cvofa[ifa][1];

          adouble muSA  = rho[icv1]*nuSA[icv1];
          //adouble mul_icv0 = InterpolateAtCellCenterFromFaceValues(mul_fa, icv0);
          adouble chi_3 = pow(muSA/mul_fa[ifa], 3.0);
          adouble fv1   = chi_3/(chi_3 + cv1_3);
          mut_fa[ifa] = max(muSA*fv1, 0.0);    // zero order extrapolation for others
        }
    }

  }




  virtual void diffusivityHookScalarRansTurb_AD(const string &name)
  {
    if (name == "nuSA")
    {
      //ScalarTranspEq *eq = getScalarTransportData(name);
      double invSigma = 1.0/sigma;
      
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

        adouble rho_fa  = (w1*rho_AD[icv0] + w0*rho_AD[icv1])/(w0+w1);
        adouble nuSA_fa = (w1*nuSA[icv0] + w0*nuSA[icv1])/(w0+w1);

        nuSA_diff[ifa] = invSigma*(mul_fa[ifa] + rho_fa*nuSA_fa);
      }
  }

      
      // boundary faces
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))
   	for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
            int icv0 = cvofa[ifa][0];
            nuSA_diff[ifa] = invSigma*mul_fa[ifa];
          }
        else
   	for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
            int icv1 = cvofa[ifa][1];
            nuSA_diff[ifa] = invSigma*(mul_fa[ifa] + rho_AD[icv1]*nuSA[icv1]);
          }
      }
    }
  }

  virtual void sourceHookScalarRansTurb_new_AD(adouble *rhs, double *A, const string &name, int flagImplicit)
  {

    if (name == "nuSA")
    {
    // update strain rate tensor  and vorticity
    calcStrainRate_AD();
    calcVorticity_AD();

      // calculate gradient of nuSA
      //calcCvScalarGrad_AD(grad_nuSA, nuSA, NULL, gradreconstruction, limiterNavierS, nuSA, epsilonSDWLS);

      for (int icv=0; icv<ncv_g; icv++)
      {
        // model variables
        adouble mue        = InterpolateAtCellCenterFromFaceValues(mul_fa, icv);
        adouble Om         = max(vortMag[icv], 1.e-8);   // vorticity
        adouble Strainrate = max(strMag[icv], 1.e-8);    // strainrate
        double d = wallDist[icv];

        // model functions
        adouble chi   = nuSA[icv]*rho_AD[icv]/mue;
        adouble chi_3 = chi*chi*chi;

        double inv_KarmanConst2_d2 = 1.0/(KarmanConst_2*d*d);

        adouble fv1 = chi_3/(chi_3+cv1_3);
        adouble fv2 = 1.0 - chi/(1.0+chi*fv1);

        double BlendTurbulentProduction = 0.0;
        adouble S = max(Om + min(0.0, Strainrate- Om)*BlendTurbulentProduction, 0.0);
        adouble Shat = S + nuSA[icv]*fv2*inv_KarmanConst2_d2;

        adouble inv_Shat =  1.0/max(Shat, 1.0e-8);

        /*  if(fabs(Shat) < 1.0e-13)
          inv_Shat = 1.0/(1.0e-13*sign(Shat));
        else
          inv_Shat = 1.0/Shat;  */

        adouble r   = min(nuSA[icv]*inv_Shat*inv_KarmanConst2_d2, 10.0);
        adouble g   = r + cw2*(pow(r, 6.0)-r);
        adouble g_6 = pow(g, 6.0);
        adouble fw  = g*pow((1.0+cw3_6)/(g_6+cw3_6), 1.0/6.0);

        // model sources
        adouble H1 =  cb1*Shat*nuSA[icv];
        adouble H2 = -cw1*fw*nuSA[icv]*nuSA[icv]/pow(d, 2.0);
        adouble H3 =  cb2/sigma*vecDotVec3d_AD(grad_nuSA[icv], grad_nuSA[icv]);

        rhs[icv] += rho_AD[icv]*cv_volume[icv]*(H1 + H2 + H3);

        if (flagImplicit && icv<ncv)
        {
          double d_H1 =  0.0;//cb1*Shat;
          double d_H2 = -2.0*cw1*fw.value()*nuSA[icv].value()/pow(d, 2.0);

          int noc00 = nbocv_i[icv];
          A[noc00] -= (d_H1+d_H2)*cv_volume[icv];
        }
      }

    }
  }

  virtual void sourceHookRansTurbCoupled_AD(adouble **rhs, double ***A, int flagImplicit)
  {
    int nuSA_Index = getScalarTransportIndex("nuSA");
    
    // update strain rate tensor  and vorticity
    calcStrainRate_AD();
    calcVorticity_AD();

      // calculate gradient of nuSA
      //calcCvScalarGrad_AD(grad_nuSA, nuSA, NULL, gradreconstruction, limiterNavierS, nuSA, epsilonSDWLS);

      for (int icv=0; icv<ncv_g; icv++)
      {
        // model variables
        adouble mue        = InterpolateAtCellCenterFromFaceValues(mul_fa, icv);
        adouble Om         = max(vortMag[icv], 1.e-8);   // vorticity
        adouble Strainrate = max(strMag[icv], 1.e-8);    // strainrate
        double d = wallDist[icv];

        // model functions
        adouble chi   = nuSA[icv]*rho_AD[icv]/mue;
        adouble chi_3 = chi*chi*chi;

        double inv_KarmanConst2_d2 = 1.0/(KarmanConst_2*d*d);

        adouble fv1 = chi_3/(chi_3+cv1_3);
        adouble fv2 = 1.0 - chi/(1.0+chi*fv1);

        double BlendTurbulentProduction = 0.0;
        adouble S = max(Om + min(0.0, Strainrate- Om)*BlendTurbulentProduction, 0.0);
        adouble Shat = S + nuSA[icv]*fv2*inv_KarmanConst2_d2;

        adouble inv_Shat =  1.0/max(Shat, 1.0e-8);

        /*  if(fabs(Shat) < 1.0e-13)
          inv_Shat = 1.0/(1.0e-13*sign(Shat));
        else
          inv_Shat = 1.0/Shat;  */

        adouble r   = min(nuSA[icv]*inv_Shat*inv_KarmanConst2_d2, 10.0);
        adouble g   = r + cw2*(pow(r, 6.0)-r);
        adouble g_6 = pow(g, 6.0);
        adouble fw  = g*pow((1.0+cw3_6)/(g_6+cw3_6), 1.0/6.0);

        // model sources
        adouble H1 =  cb1*Shat*nuSA[icv];
        adouble H2 = -cw1*fw*nuSA[icv]*nuSA[icv]/pow(d, 2.0);
        adouble H3 =  cb2/sigma*vecDotVec3d_AD(grad_nuSA[icv], grad_nuSA[icv]);

        rhs[icv][5+nuSA_Index] += rho_AD[icv]*cv_volume[icv]*(H1 + H2 + H3);

        if (flagImplicit && icv<ncv)
        {
          double d_H1 =  0.0;//cb1*Shat;
          double d_H2 = -2.0*cw1*fw.value()*nuSA[icv].value()/pow(d, 2.0);

          int noc00 = nbocv_i[icv];
          A[noc00][5+nuSA_Index][5+nuSA_Index] -= (d_H1+d_H2)*cv_volume[icv];
        }
      }

  }





};


#endif
