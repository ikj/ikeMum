#ifndef RANSTURBMODEL_KOMSST_H
#define RANSTURBMODEL_KOMSST_H

#include "UgpWithCvCompFlow.h"



class RansTurbKOmSST : virtual public UgpWithCvCompFlow
{
public:   // constructors

  RansTurbKOmSST()
  {
    if (mpi_rank == 0)
      cout << "RansTurbKOmSST()" << endl;

    turbModel = KOMSST;

    sigma_k1  = getDoubleParam("sigma_k1",   "0.85");
    sigma_k2  = getDoubleParam("sigma_k2",   "1.0");
    sigma_om1 = getDoubleParam("sigma_om1",  "0.5");
    sigma_om2 = getDoubleParam("sigma_om2",  "0.856");
    beta_1    = getDoubleParam("beta_1",     "0.075");
    beta_2    = getDoubleParam("beta_2",     "0.0828");
    betaStar  = getDoubleParam("betaStar",   "0.09");
    a1        = getDoubleParam("a1", "0.31");
    alfa_1 = beta_1/betaStar - sigma_om1*pow(0.41, 2.0)/sqrt(betaStar);
    alfa_2 = beta_2/betaStar - sigma_om2*pow(0.41, 2.0)/sqrt(betaStar);

    SST_limitPK = getIntParam("SST_LIMIT_PK", "1");  // limiter for turbulent production: 0 ... no limiter
                                                     //                                   1 ... with Pk=min(Pk, 10*betastar*rho*kine*omega);

    SST_limiterShearStress = getIntParam("SST_LIMITER_SHEARSTRESS", "0");  // limit shear stress with: 0 ... vort magnitude,
                                                                           //                          1 ... strain rate magnitude

    ScalarTranspEq *eq;
    eq = registerScalarTransport("kine", CV_DATA);
    eq->relax = getDoubleParam("RELAX_kine", "0.6");
    eq->phiZero = getDoubleParam("ZERO_kine", "1.0e-8");
    eq->phiZeroRel = getDoubleParam("ZERO_REL_kine", "1.0e-2");
    eq->phiMaxiter = 1000;
    eq->lowerBound = 1.0e-8;
    eq->upperBound = 1.0e10;
    eq->reconstruction = getStringParam("SCALAR_TURB_RECONSTRUCTION", "STANDARD");

    eq = registerScalarTransport("omega", CV_DATA);
    eq->relax = getDoubleParam("RELAX_omega", "0.4");
    eq->phiZero = getDoubleParam("ZERO_omega", "1.0e-8");
    eq->phiZeroRel = getDoubleParam("ZERO_REL_omega", "1.0e-2");
    eq->phiMaxiter = 1000;
    eq->lowerBound = 1.0e-4;
    eq->upperBound = 1.0e15;
    eq->reconstruction = getStringParam("SCALAR_TURB_RECONSTRUCTION", "STANDARD");

    strMag      = NULL;       registerScalar(strMag, "strMag", CV_DATA);
    vortMag     = NULL;       registerScalar(vortMag, "vortMag", CV_DATA);
    diverg      = NULL;       registerScalar(diverg, "diverg", CV_DATA);
    muT         = NULL;       registerScalar(muT, "muT", CV_DATA);
    crossDiff   = NULL;       registerScalar(crossDiff, "crossDiff", CV_DATA);
    blendFuncF1 = NULL;       registerScalar(blendFuncF1, "blendFuncF1", CV_DATA);
    blendFuncF2 = NULL;       registerScalar(blendFuncF2, "blendFuncF2", CV_DATA);
    wallDist    = NULL;       registerScalar(wallDist, "wallDist",  CV_DATA);
  }

  virtual ~RansTurbKOmSST() {}

public:   // member vars

  double *omega;                                        ///< turbulent scalars, introduced to have access to variables, results in to more readable code
  double (*grad_kine)[3], (*grad_omega)[3];
  double *muT;                                          ///< turbulent viscosity at cell center for output
  double *wallDist;                                     ///< wall distance
  double *crossDiff, *blendFuncF1, *blendFuncF2;
  
  // limiter for shear stress
  int SST_limiterShearStress;
  double *limiterFunc;

  double sigma_k1, sigma_k2, sigma_om1, sigma_om2, alfa_1, alfa_2, beta_1, beta_2, betaStar, a1;

  // limiter for Pk
  int SST_limitPK;

public:   // member functions

  virtual void initialHookScalarRansTurbModel()
  {
    if (mpi_rank == 0) 
      cout << "initialHook KOM SST model" << endl;

    if (!checkParam("DO_NOT_CALC_WALLDIST"))
    {
      if (!checkDataFlag(wallDist))
      {
        for (int icv=0; icv<ncv; icv++)
          wallDist[icv] = 0.0;

        calcWallDistance(NULL, wallDist);
      }
    }
    else
    {
      for (int icv=0; icv<ncv; icv++)
        wallDist[icv] = 1.0e10;
    }
      updateCvDataG1G2(wallDist, REPLACE_DATA);  // VERY IMPORTANT TO DO THIS AS THIS IS REQUIRED TO SET BOUNDARY INFO

    if (SST_limiterShearStress == 1)      limiterFunc = strMag;   // shear stress limiter = strain rate magnitude
    else                                  limiterFunc = vortMag;  // shear stress limiter = vorticity magnitude
      
    // connect pointers
    ScalarTranspEq *eq;
    eq = getScalarTransportData("kine");      kine = eq->phi;    grad_kine = eq->grad_phi;
    eq = getScalarTransportData("omega");     omega = eq->phi;   grad_omega = eq->grad_phi;
    

    // set initial condition if parameter given in Joe.in

    double kineInit, omegaInit;
    if (checkParam("INITIAL_CONDITION_TURB"))
    {
      kineInit = getParam("INITIAL_CONDITION_TURB")->getDouble(1);
      omegaInit = getParam("INITIAL_CONDITION_TURB")->getDouble(2);
    }
    else
    {
      cout << "Could not find the parameter INITIAL_CONDITION_TURB to set the initial field "<< endl;
      throw(-1);
    }

    if (!checkScalarFlag("kine"))
      for (int icv=0; icv<ncv_ggff; icv++)
        kine[icv] = kineInit;

    if (!checkScalarFlag("omega"))
      for (int icv=0; icv<ncv_ggff; icv++)
        omega[icv] = omegaInit;


  }

  virtual void calcMenterBlendingFunctions()
  {

    calcCv2Grad(grad_kine,  kine,  limiterNavierS, kine,  epsilonSDWLS);
    calcCv2Grad(grad_omega, omega, limiterNavierS, omega, epsilonSDWLS);

    for (int icv=0; icv<ncv; icv++)
    {
      double d = wallDist[icv];
      double mue = InterpolateAtCellCenterFromFaceValues(mul_fa, icv);

      crossDiff[icv] = max(2.0*rho[icv]*sigma_om2/omega[icv]*vecDotVec3d(grad_kine[icv], grad_omega[icv]),  1.0e-10);
//      crossDiff[icv] = min(crossDiff[icv], 10000.0);

      double gamma1 = 500.0*mue/(pow(d, 2.0)*rho[icv]*omega[icv]);
      double gamma2 = 4.0*sigma_om2*rho[icv]*kine[icv]/(d*d*crossDiff[icv]);
      double gamma3 = sqrt(kine[icv])/(betaStar*omega[icv]*d);

      double gamma = min(max(gamma1, gamma3), gamma2);
      blendFuncF1[icv] = tanh(pow(gamma,4.0));
      gamma = max(2.0*gamma3, gamma1);
      blendFuncF2[icv] = tanh(pow(gamma,2.0));
    }

    updateCvDataG1G2(crossDiff, REPLACE_DATA);
    updateCvDataG1G2(blendFuncF1, REPLACE_DATA);
    updateCvDataG1G2(blendFuncF2, REPLACE_DATA);
  }
  
  virtual void calcRansTurbViscMuet()
  {
    calcGradVel();
    
    // update strain rate tensor 
    calcStrainRateAndDivergence();
    calcVorticity();    
    
    calcMenterBlendingFunctions();

    
    // internal faces
    for (int ifa=nfa_b; ifa<nfa; ifa++)
    {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];

      double dx0[3], dx1[3];
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      double w0 = sqrt(vecDotVec3d(dx0, dx0));
      double w1 = sqrt(vecDotVec3d(dx1, dx1));

      double rho_fa = (w1*rho[icv0] + w0*rho[icv1])/(w0+w1);
      double kine_fa = (w1*kine[icv0] + w0*kine[icv1])/(w0+w1);
      double om_fa = (w1*omega[icv0] + w0*omega[icv1])/(w0+w1);
      double limitFunc_fa = (w1*limiterFunc[icv0] + w0*limiterFunc[icv1])/(w0+w1);
      double blendFuncF2_fa = (w1*blendFuncF2[icv0] + w0*blendFuncF2[icv1])/(w0+w1);

      double zeta = min(1.0/om_fa, a1/(limitFunc_fa*blendFuncF2_fa));

      mut_fa[ifa] = min(max(rho_fa*kine_fa*zeta, 0.0), 1.0);
    }

    // boundary faces
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      if (zoneIsWall(zone->getName()))
        for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          mut_fa[ifa] = 0.0;                                     // set mut zero at walls
      else
        for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
        {
          int icv0 = cvofa[ifa][0];
          int icv1 = cvofa[ifa][1];

          double zeta = min(1.0/omega[icv1], a1/(limiterFunc[icv0]*blendFuncF2[icv0]));
          mut_fa[ifa] = min(max(rho[icv1]*kine[icv1]*zeta, 0.0), 1.0);
        }
    }
    
    // just for output
    for (int icv=0; icv<ncv; icv++)
      muT[icv] = InterpolateAtCellCenterFromFaceValues(mut_fa, icv);
  }

  virtual void diffusivityHookScalarRansTurb(const string &name)
  {
    ScalarTranspEq *eq;

    if (name == "kine")
    {
      eq = getScalarTransportData(name);
      
      // internal faces
      for (int ifa=nfa_b; ifa<nfa; ifa++)
      {
        int icv0 = cvofa[ifa][0];
        int icv1 = cvofa[ifa][1];

        double dx0[3], dx1[3];
        vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
        vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
        double w0 = sqrt(vecDotVec3d(dx0, dx0));
        double w1 = sqrt(vecDotVec3d(dx1, dx1));

        double f1_fa = (w1*blendFuncF1[icv0] + w0*blendFuncF1[icv1])/(w0+w1);

        double coeff = sigma_k1*f1_fa + (1. - f1_fa)*sigma_k2;
        eq->diff[ifa] = mul_fa[ifa] + coeff*mut_fa[ifa];
      }

      // boundary faces
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            eq->diff[ifa] = mul_fa[ifa];
          }
        else
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            double coeff = sigma_k1*blendFuncF1[icv0] + (1. - blendFuncF1[icv0])*sigma_k2;
            eq->diff[ifa] = mul_fa[ifa] + coeff*mut_fa[ifa];
          }
      }
    }

    if (name == "omega")
    {
      eq = getScalarTransportData(name);

      // internal faces
      for (int ifa=nfa_b; ifa<nfa; ifa++)
      {
        int icv0 = cvofa[ifa][0];
        int icv1 = cvofa[ifa][1];

        double dx0[3], dx1[3];
        vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
        vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
        double w0 = sqrt(vecDotVec3d(dx0, dx0));
        double w1 = sqrt(vecDotVec3d(dx1, dx1));

        double f1_fa = (w1*blendFuncF1[icv0] + w0*blendFuncF1[icv1])/(w0+w1);

        double coeff = sigma_om1*f1_fa + (1. - f1_fa)*sigma_om2;
        eq->diff[ifa] = mul_fa[ifa] + coeff*mut_fa[ifa];
      }

      // boundary faces
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            eq->diff[ifa] = mul_fa[ifa];
          }
        else
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            double coeff = sigma_om1*blendFuncF1[icv0] + (1. - blendFuncF1[icv0])*sigma_om2;
            eq->diff[ifa] = mul_fa[ifa] + coeff*mut_fa[ifa];
          }
      }
    }
  }



    virtual double calcTurbProd(int icv, int limiter)
  {
//    double ome = max(omega[icv], 1.0*strMag[icv]/sqrt(betaStar));
    double ome = omega[icv];
    
    double zeta = min(1.0/ome, a1/(limiterFunc[icv]*blendFuncF2[icv]));
    double mut = min(max(rho[icv]*kine[icv]*zeta, 0.0), 1.0e5);
    double Pk = mut*strMag[icv]*strMag[icv] - 2./3.*rho[icv]*kine[icv]*diverg[icv];
    
    if (limiter == 1)      return min(Pk, 10.0*betaStar*rho[icv]*kine[icv]*omega[icv]);
    else                   return Pk;
  } 
  


virtual void sourceHookScalarRansTurb(double *rhs, double *A, const string &name, int flagImplicit)
  {
    if (name == "kine")
      for (int icv=0; icv<ncv; icv++)
      {
        double src = calcTurbProd(icv, SST_limitPK) - betaStar*rho[icv]*omega[icv]*kine[icv];
        rhs[icv] += src*cv_volume[icv];
        
        if (flagImplicit)
        {
          int noc00 = nbocv_i[icv];
          double dsrcdphi = - betaStar*omega[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }


    if (name == "omega")    // gamma1*str*str - beta1*rho*om*om
    {
      for (int icv=0; icv<ncv; icv++)
      {
        double F1 = blendFuncF1[icv];
        double alfa = F1*alfa_1 + (1.0 - F1)*alfa_2;
        double beta = F1*beta_1 + (1.0 - F1)*beta_2;

        double zeta = min(1.0/omega[icv], a1/(limiterFunc[icv]*blendFuncF2[icv]));
        double mut  = min(max(rho[icv]*kine[icv]*zeta, 1.0e-8), 1.0e5);

        double src = alfa*rho[icv]/mut*calcTurbProd(icv, SST_limitPK) + (1.0-F1)*crossDiff[icv] - beta*rho[icv]*omega[icv]*omega[icv];
        rhs[icv] += src*cv_volume[icv];

        if (flagImplicit)
        {
          int noc00 = nbocv_i[icv];
          double dsrcdphi =  - 2.0*beta*omega[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }
    }





  }


    virtual void sourceHookRansTurbCoupled(double **rhs, double ***A, int nScal, int flagImplicit)
  {
    int kine_Index = getScalarTransportIndex("kine");
    int omega_Index = getScalarTransportIndex("omega");

      for (int icv=0; icv<ncv; icv++)
      {
        double src = calcTurbProd(icv, SST_limitPK) - betaStar*rho[icv]*omega[icv]*kine[icv];
        rhs[icv][5+kine_Index] += src*cv_volume[icv];
        
        if (flagImplicit)
        {
          int noc00 = nbocv_i[icv];
          double dsrcdphi = - betaStar*omega[icv];
          A[noc00][5+kine_Index][5+kine_Index] -= dsrcdphi*cv_volume[icv];
        }
      }


      for (int icv=0; icv<ncv; icv++)
      {
        double F1 = blendFuncF1[icv];
        double alfa = F1*alfa_1 + (1.0 - F1)*alfa_2;
        double beta = F1*beta_1 + (1.0 - F1)*beta_2;

        double zeta = min(1.0/omega[icv], a1/(limiterFunc[icv]*blendFuncF2[icv]));
        double mut  = min(max(rho[icv]*kine[icv]*zeta, 1.0e-8), 1.0e5);

        double src = alfa*rho[icv]/mut*calcTurbProd(icv, SST_limitPK) + (1.0-F1)*crossDiff[icv] - beta*rho[icv]*omega[icv]*omega[icv];
        rhs[icv][5+omega_Index] += src*cv_volume[icv];

        if (flagImplicit)
        {
          int noc00 = nbocv_i[icv];
          double dsrcdphi =  - 2.0*beta*omega[icv];
          A[noc00][5+omega_Index][5+omega_Index] -= dsrcdphi*cv_volume[icv];
        }
      }



  }






  virtual void boundaryHookScalarRansTurb(double *phi_fa, FaZone *zone, const string &name)
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
            int icv  = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];
//            double muLamCV = calcMuLam(temp[icv]);
            double muLamCV = calcMuLam(icv);
            phi_fa[icv1] = 60.0*muLamCV/(rho[icv]*beta_1*wallDist[icv]*wallDist[icv]);
          }
        }
    }
  }
};



#endif


