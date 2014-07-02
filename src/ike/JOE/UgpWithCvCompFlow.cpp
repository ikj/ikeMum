#include <math.h>
#include "UgpWithCvCompFlow.h"

double UgpWithCvCompFlow::calcDt(double cfl_target)
{
  static int calcTimeStep = 1;    // if timeStepMode == "DT" set calcTimeStep = 0, to avoid resetting dt
  static double dt;

  if (calcTimeStep == 0)
    return dt;                    // return dt, which has been set in timeStepMode == "DT"


  if (timeStepMode == "DT")
  {
    dt = getDoubleParam("DT");    // set dt
    calcTimeStep = 0;             // set calcTimeStep = 0 to avoid resetting local_dt[icv] = dt

    for (int icv = 0; icv < ncv; icv++)
      local_dt[icv] = dt;
  }


  if ((timeStepMode == "CFL") || (timeStepMode == "CFL_MIN"))
  {
    double dt_minCPU = 1.0e20;

   string TimeStep = getStringParam("TIME_STEP","MAX_BASED");

   if (TimeStep == "MAX_BASED") 
    for (int icv = 0; icv < ncv; icv++)
    {
      double dt_cv = 1.0e20;
      double lambdaMax = 0.0;
      
      double c = sqrt(gamma[icv]*press[icv]/rho[icv]);

      int foc_f = faocv_i[icv];
      int foc_l = faocv_i[icv+1]-1;

      for (int foc=foc_f; foc<=foc_l; foc++)
      {
        int ifa = faocv_v[foc];

        double nVec[3];
        double area = normVec3d(nVec, fa_normal[ifa]);

        double Uk = vecDotVec3d(rhou[icv], nVec)/rho[icv];
        double lambda = (fabs(Uk) + c)*area;
	double lamdiff = 0.0;
        if (mu_ref > 0.0) 
             lamdiff = 4.0 * (mul_fa[ifa] + mut_fa[ifa]) / rho[icv] * area * area / cv_volume[icv];
        lambda = max(lambda, lamdiff);
        lambdaMax = max(lambdaMax, lambda);
      }

      dt_cv = cfl_target*cv_volume[icv]/lambdaMax;

      dt_minCPU = min(dt_minCPU, dt_cv);

      local_dt[icv] = dt_cv;
    }


      if (TimeStep == "SUM_BASED") 
                  for (int icv = 0; icv < ncv; icv++) {
                          double dt_cv, dt_conv, dt_visc;
                          double c = sqrt(gamma[icv]*press[icv]/rho[icv]);

                          int foc_f = faocv_i[icv];
                          int foc_l = faocv_i[icv+1]-1;

                          double lambda_conv = 0.0;
                          double lambda_diff = 0.0;
                          for (int foc=foc_f; foc<=foc_l; foc++)
                          {
                                  int ifa = faocv_v[foc];
                                  double nVec[3];
                                  double area = normVec3d(nVec, fa_normal[ifa]);
                                  double Uk = vecDotVec3d(rhou[icv], nVec)/rho[icv];
                                  double volume = 0.0;
                                  int icv0 = cvofa[ifa][0]; int icv1 = cvofa[ifa][1];
                                  if (( icv0 >= 0 && icv0<ncv_gg) && ( icv1 >= 0 && icv1 <ncv_gg )) volume = 0.5*(cv_volume[icv0] + cv_volume[icv1]);
                                  else volume = 0.5*cv_volume[icv];

                                  lambda_conv += (fabs(Uk) + c)*area;
                                  if (mu_ref > 0.0) {
                                          double KinVisc = (4.0/3.0) * (mul_fa[ifa] + mut_fa[ifa]) / rho[icv];
                                          lambda_diff +=  area * area * KinVisc / volume;
                                  }
                                  else lambda_diff = 0.0;
			
                          }

                          dt_conv = cfl_target*cv_volume[icv] / lambda_conv;
                          if (mu_ref > 0.0) dt_visc = 0.25*cfl_target*cv_volume[icv] / lambda_diff;
                          else dt_visc = 1.0e20;

                          dt_cv = min(dt_conv, dt_visc);
//                        dt_cv = cfl_target*cv_volume[icv]*(1.0/(lambda_conv*lambda_diff))/((1.0/lambda_conv)+(1.0/lambda_diff));
                          dt_minCPU = min(dt_minCPU, dt_cv);
                          local_dt[icv] = dt_cv;
                  }
          
            updateCvDataG1G2(local_dt, REPLACE_DATA);


    	    MPI_Allreduce(&dt_minCPU, &dt, 1, MPI_DOUBLE, MPI_MIN, mpi_comm);
    
    if (timeStepMode == "CFL_MIN")
      for (int icv = 0; icv < ncv_gg; icv++)
        local_dt[icv] = dt;
  }
  

  return dt;
}

// IKJ
/*
 * Method: calcDtForCFLone
 * -----------------------
 * Original code = UgpWithCvCompFlow::calcDt()
 * Calculate Dt corresponding to CFL=1.0 using the "MAX_BASED" formulation
 * This will not affect the "local_dt" variables, so it is safe to call this function anytime
 */
double UgpWithCvCompFlow::calcDtForCFLone() {
	calcStateVariables();     // Compute velocity, pressure, temperature, enthalpy and speed of sound at cell centers
	calcMaterialProperties(); // Compute viscosity and thermal diffusivity at cell faces

	double dt_minCPU = 2.2e22;

	double cfl_target = 1.0;
	for (int icv = 0; icv < ncv; icv++) { // "MAX_BASED" calculation of CFL
		double dt_cv = 1.0e20;
		double lambdaMax = 0.0;

		double c = sqrt(gamma[icv]*press[icv]/rho[icv]);

		int foc_f = faocv_i[icv];
		int foc_l = faocv_i[icv+1]-1;

		for (int foc=foc_f; foc<=foc_l; foc++) {
			int ifa = faocv_v[foc];

			double nVec[3];
			double area = normVec3d(nVec, fa_normal[ifa]);

			double Uk = vecDotVec3d(rhou[icv], nVec)/rho[icv];
			double lambda = (fabs(Uk) + c)*area;
			double lamdiff = 0.0;
			if (mu_ref > 0.0)
				lamdiff = 4.0 * (mul_fa[ifa] + mut_fa[ifa]) / rho[icv] * area * area / cv_volume[icv];
			lambda = max(lambda, lamdiff);
			lambdaMax = max(lambdaMax, lambda);
		}

		dt_cv = cfl_target*cv_volume[icv]/lambdaMax;

		dt_minCPU = min(dt_minCPU, dt_cv);
	}

	double dt_min;
	MPI_Allreduce(&dt_minCPU, &dt_min, 1, MPI_DOUBLE, MPI_MIN, mpi_comm);

	return dt_min;
}


/**
 *
 *  explicit euler flux HLLC
 *
 */
int UgpWithCvCompFlow::calcEulerFlux_HLLC(double &Frho, double *Frhou, double &FrhoE, double *FrhoScal,
        const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *ScalL, const double kL,
        const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *ScalR, const double kR,
        const double area, const double *nVec, const int nScal, const double surfVeloc)
{

  double surfv = 0.0;
  //double surfv = -0.1*nVec[0];

  double unL  = vecDotVec3d(uL, nVec) - surfv;
  double unLr = unL - surfv;
  double uLuL = vecDotVec3d(uL, uL);
  double cL   = sqrt(gammaL*pL/rhoL);
  double hL   = gammaL/(gammaL-1.0)*pL/rhoL + 0.5*uLuL + kL;
//  double hL   = h0 + 0.5*uLuL + kL;
  double eL   = hL*rhoL-pL;

  double unR  = vecDotVec3d(uR, nVec) - surfv;
  double unRr = unR - surfv;
  double uRuR = vecDotVec3d(uR, uR);
  double cR   = sqrt(gammaR*pR/rhoR);
  double hR   = gammaR/(gammaR-1.0)*pR/rhoR + 0.5*uRuR + kR;
//  double hR   = h1 + 0.5*uRuR + kR;
  double eR   = hR*rhoR-pR;


  // Roe's aveaging
  double Rrho = sqrt(rhoR/rhoL);
  double tmp = 1.0/(1.0+Rrho);
  double velRoe[3];
  for (int i=0; i<3; i++)
    velRoe[i] = tmp*(uL[i] + uR[i]*Rrho);
  double uRoe  = vecDotVec3d(velRoe, nVec);
//  double hRoe = tmp*(hL + hR*Rrho);
//  double cRoe  = sqrt((gammaL-1.0)*(hRoe- 0.5*vecDotVec3d(velRoe, velRoe)));
  double gamPdivRho = tmp*((gammaL*pL/rhoL+0.5*(gammaL-1.0)*uLuL) + (gammaR*pR/rhoR+0.5*(gammaR-1.0)*uRuR)*Rrho);
  double cRoe  = sqrt(gamPdivRho - ((gammaL+gammaR)*0.5-1.0)*0.5*vecDotVec3d(velRoe, velRoe));

  // speed of sound at L and R
  double sL = min(uRoe-cRoe, unLr-cL);
  double sR = max(uRoe+cRoe, unRr+cR);

  // speed of contact surface
  double sM = (pL-pR-rhoL*unL*(sL-unL)+rhoR*unR*(sR-unR))/(rhoR*(sR-unR)-rhoL*(sL-unL));

  // pressure at right and left (pR=pL) side of contact surface
  double pStar = rhoR*(unR-sR)*(unR-sM)+pR;

  if (sM >= 0.0)
  {
    if (sL > 0.0)
    {
      Frho = rhoL*unLr;
      for (int i=0; i<3; i++)
        Frhou[i] = rhoL*uL[i]*unLr + pL*nVec[i];
      FrhoE = eL*unLr + pL*unL;
      for (int iScal=0; iScal<nScal; iScal++)
        FrhoScal[iScal] = Frho * ScalL[iScal];
    }
    else
    {
      double invSLmSs = 1.0/(sL-sM);
      double sLmuL = sL-unL;
      double rhoSL = rhoL*sLmuL*invSLmSs;
      double rhouSL[3];
      for (int i=0; i<3; i++)
        rhouSL[i] = (rhoL*uL[i]*sLmuL+(pStar-pL)*nVec[i])*invSLmSs;
      double eSL = (sLmuL*eL-pL*unL+pStar*sM)*invSLmSs;

      Frho = rhoSL*sM;
      for (int i=0; i<3; i++)
        Frhou[i] = rhouSL[i]*sM + pStar*nVec[i];
      FrhoE = (eSL+pStar)*sM;
      for (int iScal=0; iScal<nScal; iScal++)
        FrhoScal[iScal] = Frho * ScalL[iScal];
    }
  }
  else
  {
    if (sR >= 0.0)
    {
      double invSRmSs = 1.0/(sR-sM);
      double sRmuR = sR-unR;
      double rhoSR = rhoR*sRmuR*invSRmSs;
      double rhouSR[3];
      for (int i=0; i<3; i++)
        rhouSR[i] = (rhoR*uR[i]*sRmuR+(pStar-pR)*nVec[i])*invSRmSs;
      double eSR = (sRmuR*eR-pR*unR+pStar*sM)*invSRmSs;

      Frho = rhoSR*sM;
      for (int i=0; i<3; i++)
        Frhou[i] = rhouSR[i]*sM + pStar*nVec[i];
      FrhoE = (eSR+pStar)*sM;
      for (int iScal=0; iScal<nScal; iScal++)
        FrhoScal[iScal] = Frho * ScalR[iScal];
    }
    else
    {
      Frho = rhoR*unRr;
      for (int i=0; i<3; i++)
        Frhou[i] = rhoR*uR[i]*unRr + pR*nVec[i];
      FrhoE = eR*unRr + pR*unR;
      for (int iScal=0; iScal<nScal; iScal++)
        FrhoScal[iScal] = Frho * ScalR[iScal];
    }
  }

  Frho *= area;
  for (int i=0; i<3; i++)
    Frhou[i] *= area;
  FrhoE *= area;
  for (int iScal=0; iScal<nScal; iScal++)
    FrhoScal[iScal] *= area;

  return 0;
}

/**
 *
 *  implicit euler flux HLLC matrix
 *
 */
int UgpWithCvCompFlow::calcEulerFluxMatrices_HLLC(double (*A_L)[5], double (*A_R)[5], double (*A_L_Scal)[6], double (*A_R_Scal)[6],
        const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *scalL, const double kL,
        const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *scalR, const double kR,
        const double area, const double *nVec, const int nScal, const double surfVeloc)
{  
  double unL  = vecDotVec3d(uL, nVec);
  double uLuL = vecDotVec3d(uL, uL);
  double cL   = sqrt(gammaL*pL/rhoL);
  double hL   = gammaL/(gammaL-1.0)*pL/rhoL + 0.5*uLuL + kL;
//  double hL   = h0 + 0.5*uLuL + kL;
  double eL   = hL*rhoL-pL;

  double unR  = vecDotVec3d(uR, nVec);
  double uRuR = vecDotVec3d(uR, uR);
  double cR   = sqrt(gammaR*pR/rhoR);
  double hR   = gammaR/(gammaR-1.0)*pR/rhoR + 0.5*uRuR + kR;
//  double hR   = h1 + 0.5*uRuR + kR;
  double eR   = hR*rhoR-pR;


  // Roe's aveaging
  double Rrho = sqrt(rhoR/rhoL);
  double tmp = 1.0/(1.0+Rrho);
  double velRoe[3];
  for (int i=0; i<3; i++)
    velRoe[i] = tmp*(uL[i] + uR[i]*Rrho);
  double uRoe  = vecDotVec3d(velRoe, nVec);
  double hRoe = tmp*(hL + hR*Rrho);
  
//  double cRoe  = sqrt((gammaL-1.0)*(hRoe- 0.5*vecDotVec3d(velRoe, velRoe)));
  double gamPdivRho = tmp*((gammaL*pL/rhoL+0.5*(gammaL-1.0)*uLuL) + (gammaR*pR/rhoR+0.5*(gammaR-1.0)*uRuR)*Rrho);
  double cRoe  = sqrt(gamPdivRho - ((gammaL+gammaR)*0.5-1.0)*0.5*vecDotVec3d(velRoe, velRoe));

  // speed of sound at L and R
  double sL = min(uRoe-cRoe, unL-cL);
  double sR = max(uRoe+cRoe, unR+cR);

  // speed of contact surface
  double sM = (pL-pR-rhoL*unL*(sL-unL)+rhoR*unR*(sR-unR))/(rhoR*(sR-unR)-rhoL*(sL-unL));

  // pressure at right and left (pR=pL) side of contact surface
  double pStar = rhoR*(unR-sR)*(unR-sM)+pR;

  if (sM >= 0.0)
  {
    if (sL > 0.0)
    {
      if (A_L != NULL)
      {
        double nVecArea[3];
        for (int i=0; i<3; i++)        nVecArea[i] = nVec[i]*area;
        calcJacobianA(A_L, uL, pL, rhoL, nVecArea, 0.5*(gammaL+gammaR), 0.0);

        for (int iScal = 0; iScal < nScal; iScal++)
        {
          // F(rhophi_L) = F(rho_L) * phi_L = F(rho_L) * rhophi_L / rho_L = rho_L * unL * area * rhophi_L / rho_L 
          
          // dF(rhophi_L) / drho_L    = dF(rho_L) / drho_L  * phi_L - F(rho_L) * phi_L / rho_L = A_L[0,0] * phi_L - F(rho_L) * phi_L / rho_L
          // dF(rhophi_L) / drhou_L   = dF(rho_L) / drhou_L * phi_L                            = A_L[0,1] * phi_L
          // dF(rhophi_L) / drhoe_L   = dF(rho_L) / drhoe_L * phi_L                            = A_L[0,4] * phi_L
          // dF(rhophi_L) / drhophi_L =  F(rho_L) / rho_L                                      = F(rho_L) / rho_L
          
          double FrhostarOverrho = unL * area;
          A_L_Scal[iScal][0] = (A_L[0][0] - FrhostarOverrho) * scalL[iScal];    // d(rho phi un)/d(rho)
          A_L_Scal[iScal][1] =  A_L[0][1] * scalL[iScal];                       // d(rho phi un)/d(rho u)
          A_L_Scal[iScal][2] =  A_L[0][2] * scalL[iScal];                       // d(rho phi un)/d(rho v)
          A_L_Scal[iScal][3] =  A_L[0][3] * scalL[iScal];                       // d(rho phi un)/d(rho w)
          A_L_Scal[iScal][4] =  A_L[0][4] * scalL[iScal];                       // d(rho phi un)/d(rho e)
          A_L_Scal[iScal][5] =  FrhostarOverrho;                                // d(rho phi un)/d(rho phi)
        }
      }

      if (A_R != NULL)
      {
        for (int i=0; i<5; i++)
          for (int j=0; j<5; j++)
            A_R[i][j] = 0.0;

        for (int iScal = 0; iScal < nScal; iScal++)
          // fluxes are not dependent on the right values
          for (int i = 0; i <= 5; i++)
            A_R_Scal[iScal][i] = 0.0;
      }
    }
    else
    {
      double invSLmSs = 1.0/(sL-sM);
      double sLmuL = sL-unL;
      double rhoSL = rhoL*sLmuL*invSLmSs;
      double rhouSL[3];
      for (int i=0; i<3; i++)
        rhouSL[i] = (rhoL*uL[i]*sLmuL+(pStar-pL)*nVec[i])*invSLmSs;
      double eSL = (sLmuL*eL-pL*unL+pStar*sM)*invSLmSs;
      double gammaLM1 = (gammaL-1.0);
      double gammaRM1 = (gammaR-1.0);
      double invrhotld = 1.0/(rhoR*(sR-unR)-rhoL*(sL-unL));

      double dSMdUL[5], dSMdUR[5];
      double dpsdUL[5], dpsdUR[5];
      
      if (A_L != NULL)
      {
        dSMdUL[0] = -unL*unL + uLuL*gammaLM1/2.0 + sM*sL;
        dSMdUL[1] =  nVec[0]*(2.0*unL-sL-sM) - gammaLM1*uL[0];
        dSMdUL[2] =  nVec[1]*(2.0*unL-sL-sM) - gammaLM1*uL[1];
        dSMdUL[3] =  nVec[2]*(2.0*unL-sL-sM) - gammaLM1*uL[2];
        dSMdUL[4] =  gammaLM1;
        
        for (int i=0; i<5; i++)
        {
          dSMdUL[i] *= invrhotld;
          dpsdUL[i] = rhoR*(sR-unR)*dSMdUL[i];
        }
      }

      if (A_R != NULL)
      {
        dSMdUR[0] =  unR*unR - uRuR*gammaRM1/2.0 - sM*sR;
        dSMdUR[1] = -nVec[0]*(2.0*unR-sR-sM) + gammaRM1*uR[0];
        dSMdUR[2] = -nVec[1]*(2.0*unR-sR-sM) + gammaRM1*uR[1];
        dSMdUR[3] = -nVec[2]*(2.0*unR-sR-sM) + gammaRM1*uR[2];
        dSMdUR[4] = -gammaRM1;
        
        for (int i=0; i<5; i++)
        {
          dSMdUR[i] *= invrhotld;
          dpsdUR[i] = rhoL*(sL-unL)*dSMdUR[i];
        }
      }

      calcSubSonicJacobeanHLLC(A_L, A_R,
          rhoL, uL, pL, eL, unL, uLuL, sL,
          rhoSL, rhouSL, eSL, dSMdUL, 
          dSMdUR, dpsdUL, dpsdUR, sM, pStar, 0.5*(gammaL+gammaR), nVec);

      if (A_L != NULL)
      {
        for (int i=0; i<5; i++)  A_L[0][i] =  A_L[0][i]*sM + dSMdUL[i]*rhoSL;
        for (int i=0; i<5; i++)  A_L[1][i] =  A_L[1][i]*sM + dSMdUL[i]*rhouSL[0] + dpsdUL[i]*nVec[0];
        for (int i=0; i<5; i++)  A_L[2][i] =  A_L[2][i]*sM + dSMdUL[i]*rhouSL[1] + dpsdUL[i]*nVec[1];
        for (int i=0; i<5; i++)  A_L[3][i] =  A_L[3][i]*sM + dSMdUL[i]*rhouSL[2] + dpsdUL[i]*nVec[2];
        for (int i=0; i<5; i++)  A_L[4][i] = (A_L[4][i]+dpsdUL[i])*sM + (eSL+pStar)*dSMdUL[i];
        
        for (int i=0; i<5; i++)
          for (int j=0; j<5; j++)
            A_L[i][j] *= area;

        for (int iScal = 0; iScal < nScal; iScal++)
        {
          // F(rhophi_L^star) = F(rho_L^star) * phi_L = F(rho_L^star) * rhophi_L / rho_L = invSLmSM * (SL - q_L) * SM * area * rho_L * rhophi_L / rho_L
          
          // dF(rhophi_L^star) / drho_L    = dF(rho_L^star) / drho_L  * phi_L - F(rho_L^star) * phi_L / rho_L = A_L[0,0] * phi_L - F(rho_L^star) * phi_L / rho_L
          // dF(rhophi_L^star) / drhou_L   = dF(rho_L^star) / drhou_L * phi_L                                 = A_L[0,1] * phi_L
          // dF(rhophi_L^star) / drhoe_L   = dF(rho_L^star) / drhoe_L * phi_L                                 = A_L[0,4] * phi_L
          // dF(rhophi_L^star) / drhophi_L =  F(rho_L^star) / rho_L                                           = F(rho_L^star) / rho_L
          
          double FrhostarOverrho = invSLmSs * sLmuL * sM * area;
          A_L_Scal[iScal][0] = (A_L[0][0] - FrhostarOverrho) * scalL[iScal];     // d(rho phi un)/d(rho)
          A_L_Scal[iScal][1] =  A_L[0][1] * scalL[iScal];                        // d(rho phi un)/d(rho u)
          A_L_Scal[iScal][2] =  A_L[0][2] * scalL[iScal];                        // d(rho phi un)/d(rho v)
          A_L_Scal[iScal][3] =  A_L[0][3] * scalL[iScal];                        // d(rho phi un)/d(rho w)
          A_L_Scal[iScal][4] =  A_L[0][4] * scalL[iScal];                        // d(rho phi un)/d(rho e)
          A_L_Scal[iScal][5] =  FrhostarOverrho;                                 // d(rho phi un)/d(rho phi)
        }
      }

      if (A_R != NULL)
      {
        for (int i=0; i<5; i++)  A_R[0][i] =  A_R[0][i]*sM + dSMdUR[i]*rhoSL;
        for (int i=0; i<5; i++)  A_R[1][i] =  A_R[1][i]*sM + dSMdUR[i]*rhouSL[0] + dpsdUR[i]*nVec[0];
        for (int i=0; i<5; i++)  A_R[2][i] =  A_R[2][i]*sM + dSMdUR[i]*rhouSL[1] + dpsdUR[i]*nVec[1];
        for (int i=0; i<5; i++)  A_R[3][i] =  A_R[3][i]*sM + dSMdUR[i]*rhouSL[2] + dpsdUR[i]*nVec[2];
        for (int i=0; i<5; i++)  A_R[4][i] = (A_R[4][i]+dpsdUR[i])*sM + (eSL+pStar)*dSMdUR[i];
        
        for (int i=0; i<5; i++)
          for (int j=0; j<5; j++)
            A_R[i][j] *= area;

        for (int iScal = 0; iScal < nScal; iScal++)
        {
          // F(rhophi_L^star) = F(rho_L^star) * phi_L = F(rho_L^star) * rhophi_L / rho_L = invSLmSM * (SL - q_L) * SM * area * rho_L * rhophi_L / rho_L
          
          // dF(rhophi_L^star) / drho_R    = dF(rho_L^star) / drho_R  * phi_L = A_R[0,0] * phi_L
          // dF(rhophi_L^star) / drhou_R   = dF(rho_L^star) / drhou_R * phi_L = A_R[0,1] * phi_L
          // dF(rhophi_L^star) / drhoe_R   = dF(rho_L^star) / drhoe_R * phi_L = A_R[0,4] * phi_L
          // dF(rhophi_L^star) / drhophi_R                                    = 0
          
          A_R_Scal[iScal][0] = A_R[0][0] * scalL[iScal];    // d(rho phi un)/d(rho)
          A_R_Scal[iScal][1] = A_R[0][1] * scalL[iScal];    // d(rho phi un)/d(rho u)
          A_R_Scal[iScal][2] = A_R[0][2] * scalL[iScal];    // d(rho phi un)/d(rho v)
          A_R_Scal[iScal][3] = A_R[0][3] * scalL[iScal];    // d(rho phi un)/d(rho w)
          A_R_Scal[iScal][4] = A_R[0][4] * scalL[iScal];    // d(rho phi un)/d(rho e)
          A_R_Scal[iScal][5] = 0.0;                         // d(rho phi un)/d(rho phi)
        }
      }
    }
  }
  else
  {
    if (sR >= 0.0)
    {
      double invSRmSs = 1.0/(sR-sM);
      double sRmuR = sR-unR;
      double rhoSR = rhoR*sRmuR*invSRmSs;
      double rhouSR[3];
      for (int i=0; i<3; i++)
        rhouSR[i] = (rhoR*uR[i]*sRmuR+(pStar-pR)*nVec[i])*invSRmSs;
      double eSR = (sRmuR*eR-pR*unR+pStar*sM)*invSRmSs;      
      double gammaLM1 = (gammaL-1.0);
      double gammaRM1 = (gammaR-1.0);
      double invrhotld = 1.0/(rhoR*(sR-unR)-rhoL*(sL-unL));

      double dSMdUL[5], dSMdUR[5];
      double dpsdUL[5], dpsdUR[5];
      
      if (A_L != NULL)
      {
        dSMdUL[0] = -unL*unL + uLuL*gammaLM1/2.0 + sM*sL;
        dSMdUL[1] =  nVec[0]*(2.0*unL-sL-sM) - gammaLM1*uL[0];
        dSMdUL[2] =  nVec[1]*(2.0*unL-sL-sM) - gammaLM1*uL[1];
        dSMdUL[3] =  nVec[2]*(2.0*unL-sL-sM) - gammaLM1*uL[2];
        dSMdUL[4] =  gammaLM1;

        for (int i=0; i<5; i++)
        {
          dSMdUL[i] *= invrhotld;
          dpsdUL[i] = rhoR*(sR-unR)*dSMdUL[i];
        }
      }

      if (A_R != NULL)
      {
        dSMdUR[0] =  unR*unR - uRuR*gammaRM1/2.0 - sM*sR;
        dSMdUR[1] = -nVec[0]*(2.0*unR-sR-sM) + gammaRM1*uR[0];
        dSMdUR[2] = -nVec[1]*(2.0*unR-sR-sM) + gammaRM1*uR[1];
        dSMdUR[3] = -nVec[2]*(2.0*unR-sR-sM) + gammaRM1*uR[2];
        dSMdUR[4] = -gammaRM1;

        for (int i=0; i<5; i++)
        {
          dSMdUR[i] *= invrhotld;
          dpsdUR[i] = rhoL*(sL-unL)*dSMdUR[i];
        }
      }

      calcSubSonicJacobeanHLLC(A_R, A_L,
          rhoR, uR, pR, eR, unR, uRuR, sR,
          rhoSR, rhouSR, eSR,
          dSMdUR, dSMdUL, dpsdUR, dpsdUL, sM, pStar, 0.5*(gammaL+gammaR), nVec);

      if (A_L != NULL)
      {
        for (int i=0; i<5; i++)  A_L[0][i] =  A_L[0][i]*sM + dSMdUL[i]*rhoSR;
        for (int i=0; i<5; i++)  A_L[1][i] =  A_L[1][i]*sM + dSMdUL[i]*rhouSR[0] + dpsdUL[i]*nVec[0];
        for (int i=0; i<5; i++)  A_L[2][i] =  A_L[2][i]*sM + dSMdUL[i]*rhouSR[1] + dpsdUL[i]*nVec[1];
        for (int i=0; i<5; i++)  A_L[3][i] =  A_L[3][i]*sM + dSMdUL[i]*rhouSR[2] + dpsdUL[i]*nVec[2];
        for (int i=0; i<5; i++)  A_L[4][i] = (A_L[4][i]+dpsdUL[i])*sM + (eSR+pStar)*dSMdUL[i];
        
        for (int i=0; i<5; i++)
          for (int j=0; j<5; j++)
            A_L[i][j] *= area;

        for (int iScal = 0; iScal < nScal; iScal++)
        {
          // F(rhophi_R^star) = F(rho_R^star) * phi_R = F(rho_R^star) * rhophi_R / rho_R = invSRmSM * (SR - q_R) * SM * area * rho_R * rhophi_R / rho_R
          
          // dF(rhophi_R^star) / drho_L    = dF(rho_R^star) / drho_L  * phi_R = A_L[0,0] * phi_R
          // dF(rhophi_R^star) / drhou_L   = dF(rho_R^star) / drhou_L * phi_R = A_L[0,1] * phi_R
          // dF(rhophi_R^star) / drhoe_L   = dF(rho_R^star) / drhoe_L * phi_R = A_L[0,4] * phi_R
          // dF(rhophi_R^star) / drhophi_L                                    = 0
          
          A_L_Scal[iScal][0] = A_L[0][0] * scalR[iScal];    // d(rho phi un)/d(rho)
          A_L_Scal[iScal][1] = A_L[0][1] * scalR[iScal];    // d(rho phi un)/d(rho u)
          A_L_Scal[iScal][2] = A_L[0][2] * scalR[iScal];    // d(rho phi un)/d(rho v)
          A_L_Scal[iScal][3] = A_L[0][3] * scalR[iScal];    // d(rho phi un)/d(rho w)
          A_L_Scal[iScal][4] = A_L[0][4] * scalR[iScal];    // d(rho phi un)/d(rho e)
          A_L_Scal[iScal][5] = 0.0;                         // d(rho phi un)/d(rho phi)
        }
      }
      
      if (A_R != NULL)
      {
        for (int i=0; i<5; i++)  A_R[0][i] =  A_R[0][i]*sM + dSMdUR[i]*rhoSR;
        for (int i=0; i<5; i++)  A_R[1][i] =  A_R[1][i]*sM + dSMdUR[i]*rhouSR[0] + dpsdUR[i]*nVec[0];
        for (int i=0; i<5; i++)  A_R[2][i] =  A_R[2][i]*sM + dSMdUR[i]*rhouSR[1] + dpsdUR[i]*nVec[1];
        for (int i=0; i<5; i++)  A_R[3][i] =  A_R[3][i]*sM + dSMdUR[i]*rhouSR[2] + dpsdUR[i]*nVec[2];
        for (int i=0; i<5; i++)  A_R[4][i] = (A_R[4][i]+dpsdUR[i])*sM + (eSR+pStar)*dSMdUR[i];

        for (int i=0; i<5; i++)
          for (int j=0; j<5; j++)
            A_R[i][j] *= area;

        for (int iScal = 0; iScal < nScal; iScal++)
        {
          // F(rhophi_R^star) = F(rho_R^star) * phi_R = F(rho_R^star) * rhophi_R / rho_R = invSRmSM * (SR - q_R) * SM * area * rho_R * rhophi_R / rho_R
          
          // dF(rhophi_R^star) / drho_R    = dF(rho_R^star) / drho_R  * phi_R - F(rho_R^star) * phi_R / rho_R = A_R[0,0] * phi_R - F(rho_R^star) * phi_R / rho_R
          // dF(rhophi_R^star) / drhou_R   = dF(rho_R^star) / drhou_R * phi_R                                 = A_R[0,1] * phi_R
          // dF(rhophi_R^star) / drhoe_R   = dF(rho_R^star) / drhoe_R * phi_R                                 = A_R[0,4] * phi_R
          // dF(rhophi_R^star) / drhophi_R =  F(rho_R^star) / rho_R                                           = F(rho_R^star) / rho_R
          
          double FrhostarOverrho = invSRmSs * sRmuR * sM * area;
          A_R_Scal[iScal][0] = (A_R[0][0] - FrhostarOverrho) * scalR[iScal];     // d(rho phi un)/d(rho)
          A_R_Scal[iScal][1] =  A_R[0][1] * scalR[iScal];                        // d(rho phi un)/d(rho u)
          A_R_Scal[iScal][2] =  A_R[0][2] * scalR[iScal];                        // d(rho phi un)/d(rho v)
          A_R_Scal[iScal][3] =  A_R[0][3] * scalR[iScal];                        // d(rho phi un)/d(rho w)
          A_R_Scal[iScal][4] =  A_R[0][4] * scalR[iScal];                        // d(rho phi un)/d(rho e)
          A_R_Scal[iScal][5] =  FrhostarOverrho;                                 // d(rho phi un)/d(rho phi)
        }
      }
    }
    else
    {
      if (A_R != NULL)
      {
        double nVecArea[3];
        for (int i=0; i<3; i++)        nVecArea[i] = nVec[i]*area;
        calcJacobianA(A_R, uR, pR, rhoR, nVecArea, 0.5*(gammaL+gammaR), 0.0);

        for (int iScal = 0; iScal < nScal; iScal++)
        {
          // F(rhophi_R) = F(rho_R) * phi_R = F(rho_R) * rhophi_R / rho_R = rho_R * unR * area * rhophi_R / rho_R 
          
          // dF(rhophi_R) / drho_R    = dF(rho_R) / drho_R  * phi_R - F(rho_R) * phi_R / rho_R = A_R[0,0] * phi_R - F(rho_R) * phi_R / rho_R
          // dF(rhophi_R) / drhou_R   = dF(rho_R) / drhou_R * phi_R                            = A_R[0,1] * phi_R
          // dF(rhophi_R) / drhoe_R   = dF(rho_R) / drhoe_R * phi_R                            = A_R[0,4] * phi_R
          // dF(rhophi_R) / drhophi_R =  F(rho_R) / rho_R                                      = F(rho_R) / rho_R
          
          double FrhostarOverrho = unR * area;
          A_R_Scal[iScal][0] = (A_R[0][0] - FrhostarOverrho) * scalR[iScal];    // d(rho phi un)/d(rho)
          A_R_Scal[iScal][1] =  A_R[0][1] * scalR[iScal];                       // d(rho phi un)/d(rho u)
          A_R_Scal[iScal][2] =  A_R[0][2] * scalR[iScal];                       // d(rho phi un)/d(rho v)
          A_R_Scal[iScal][3] =  A_R[0][3] * scalR[iScal];                       // d(rho phi un)/d(rho w)
          A_R_Scal[iScal][4] =  A_R[0][4] * scalR[iScal];                       // d(rho phi un)/d(rho e)
          A_R_Scal[iScal][5] =  FrhostarOverrho;                                // d(rho phi un)/d(rho phi)
        }
      }

      if (A_L != NULL)
      {
        for (int i=0; i<5; i++)
          for (int j=0; j<5; j++)
            A_L[i][j] = 0.0;

        for (int iScal = 0; iScal < nScal; iScal++)
          // Fluxes are not dependent on the left values
          for (int i = 0; i <= 5; i++)
            A_L_Scal[iScal][i] = 0.0;
      }
    }
  }

  return 0;
}




/**
 *
 *  implicit euler flux HLLC matrix for subsonic part
 *
 */
void UgpWithCvCompFlow::calcSubSonicJacobeanHLLC(double (*AL)[5], double (*AR)[5],
    double rhoL, const double *uL, double pL, double eL, double qL, double psiL, double SL,
    double rhoSL, double *rhouSL, double eSL,
    double *dSMdUL, double *dSMdUR, double *dpsdUL, double *dpsdUR, double SM, double pS,
    double gamma, const double *nV) // nV is not normalized
{
  double gammaMinus1 = (gamma-1.0);
  double omL = 1.0/(SL-SM);

  if (AL != NULL)
  {
    AL[0][0] =  SL    + rhoSL*dSMdUL[0];
    AL[0][1] = -nV[0] + rhoSL*dSMdUL[1];
    AL[0][2] = -nV[1] + rhoSL*dSMdUL[2];
    AL[0][3] = -nV[2] + rhoSL*dSMdUL[3];
    AL[0][4] =        + rhoSL*dSMdUL[4];

    AL[1][0] =    qL*uL[0]       - nV[0]*psiL*gammaMinus1/2.0   + nV[0]*dpsdUL[0] + rhouSL[0]*dSMdUL[0];
    AL[1][1] =  SL - qL          + nV[0]*(gamma-2.0)*uL[0]      + nV[0]*dpsdUL[1] + rhouSL[0]*dSMdUL[1];
    AL[1][2] =     - uL[0]*nV[1] + nV[0]*gammaMinus1*uL[1]      + nV[0]*dpsdUL[2] + rhouSL[0]*dSMdUL[2];
    AL[1][3] =     - uL[0]*nV[2] + nV[0]*gammaMinus1*uL[2]      + nV[0]*dpsdUL[3] + rhouSL[0]*dSMdUL[3];
    AL[1][4] = -gammaMinus1*nV[0]                               + nV[0]*dpsdUL[4] + rhouSL[0]*dSMdUL[4];

    AL[2][0] =    qL*uL[1]       - nV[1]*psiL*gammaMinus1/2.0   + nV[1]*dpsdUL[0] + rhouSL[1]*dSMdUL[0];
    AL[2][1] =     - uL[1]*nV[0] + nV[1]*gammaMinus1*uL[0]      + nV[1]*dpsdUL[1] + rhouSL[1]*dSMdUL[1];
    AL[2][2] =  SL - qL          + nV[1]*(gamma-2.0)*uL[1]      + nV[1]*dpsdUL[2] + rhouSL[1]*dSMdUL[2];
    AL[2][3] =     - uL[1]*nV[2] + nV[1]*gammaMinus1*uL[2]      + nV[1]*dpsdUL[3] + rhouSL[1]*dSMdUL[3];
    AL[2][4] = -gammaMinus1*nV[1]                               + nV[1]*dpsdUL[4] + rhouSL[1]*dSMdUL[4];

    AL[3][0] =    qL*uL[2]       - nV[2]*psiL*gammaMinus1/2.0   + nV[2]*dpsdUL[0] + rhouSL[2]*dSMdUL[0];
    AL[3][1] =     - uL[2]*nV[0] + nV[2]*gammaMinus1*uL[0]      + nV[2]*dpsdUL[1] + rhouSL[2]*dSMdUL[1];
    AL[3][2] =     - uL[2]*nV[1] + nV[2]*gammaMinus1*uL[1]      + nV[2]*dpsdUL[2] + rhouSL[2]*dSMdUL[2];
    AL[3][3] =  SL - qL          + nV[2]*(gamma-2.0)*uL[2]      + nV[2]*dpsdUL[3] + rhouSL[2]*dSMdUL[3];
    AL[3][4] = -gammaMinus1*nV[2]                               + nV[2]*dpsdUL[4] + rhouSL[2]*dSMdUL[4];

    AL[4][0] =      qL*(eL+pL)/rhoL - qL*psiL*(gamma-1.0)/2.0   + SM*dpsdUL[0] + (pS+eSL)*dSMdUL[0];
    AL[4][1] = - nV[0]*(eL+pL)/rhoL + gammaMinus1*uL[0]*qL      + SM*dpsdUL[1] + (pS+eSL)*dSMdUL[1];
    AL[4][2] = - nV[1]*(eL+pL)/rhoL + gammaMinus1*uL[1]*qL      + SM*dpsdUL[2] + (pS+eSL)*dSMdUL[2];
    AL[4][3] = - nV[2]*(eL+pL)/rhoL + gammaMinus1*uL[2]*qL      + SM*dpsdUL[3] + (pS+eSL)*dSMdUL[3];
    AL[4][4] =   SL-qL*gamma                                    + SM*dpsdUL[4] + (pS+eSL)*dSMdUL[4];

    for (int i=0; i<5; i++)
      for (int j=0; j<5; j++)
        AL[i][j] *= omL;
  }

  if (AR != NULL)
  {
    for (int i=0; i<5; i++)    AR[0][i] = omL*rhoSL*dSMdUR[i];
    for (int i=0; i<5; i++)    AR[1][i] = omL*(nV[0]*dpsdUR[i]+rhouSL[0]*dSMdUR[i]);
    for (int i=0; i<5; i++)    AR[2][i] = omL*(nV[1]*dpsdUR[i]+rhouSL[1]*dSMdUR[i]);
    for (int i=0; i<5; i++)    AR[3][i] = omL*(nV[2]*dpsdUR[i]+rhouSL[2]*dSMdUR[i]);
    for (int i=0; i<5; i++)    AR[4][i] = omL*(dpsdUR[i]*SM+(pS+eSL)*dSMdUR[i]);
  }
}

/**
 *
 * jacobean
 *
 */
void UgpWithCvCompFlow::calcJacobianA(double (*A)[5], const double *vel, double pp, double rrho, const double *nV, double gamma, double surfVeloc) // nV is not normalized
{
  double kapm1 = (gamma - 1.0);

  double nVel[3];
  nVel[0] = vel[0]*nV[0];
  nVel[1] = vel[1]*nV[1];
  nVel[2] = vel[2]*nV[2];
  double U_k = nVel[0]+nVel[1]+nVel[2];
  double vSquHlf = 0.5*vecDotVec3d(vel, vel);
  double c = sqrt(gamma*pp/rrho);
  double inv_kap_m1 = 1.0/kapm1;

  // A Matrix
  A[0][0] =-surfVeloc;
  A[0][1] = nV[0];
  A[0][2] = nV[1];
  A[0][3] = nV[2];
  A[0][4] = 0.0;

  A[1][0] = -vel[0]*(nVel[1]+nVel[2])+nV[0]*(kapm1*vSquHlf-vel[0]*vel[0]);
  A[1][1] = (2.-gamma)*nVel[0]+U_k-surfVeloc;
  A[1][2] = vel[0]*nV[1]-kapm1*vel[1]*nV[0];
  A[1][3] = vel[0]*nV[2]-kapm1*vel[2]*nV[0];
  A[1][4] = kapm1*nV[0];

  A[2][0] = -vel[1]*(nVel[0]+nVel[2])+nV[1]*(kapm1*vSquHlf-vel[1]*vel[1]);
  A[2][1] = -kapm1*vel[0]*nV[1]+ vel[1]*nV[0];
  A[2][2] = (2.-gamma)*nVel[1]+U_k-surfVeloc;
  A[2][3] = vel[1]*nV[2]-kapm1*vel[2]*nV[1];
  A[2][4] = kapm1*nV[1];

  A[3][0] = -vel[2]*(nVel[0]+nVel[1])+nV[2]*(kapm1*vSquHlf-vel[2]*vel[2]);
  A[3][1] = -kapm1*vel[0]*nV[2]+vel[2]*nV[0];
  A[3][2] = -kapm1*vel[1]*nV[2]+vel[2]*nV[1];
  A[3][3] = (2.-gamma)*nVel[2]+U_k-surfVeloc;
  A[3][4] = kapm1*nV[2];

  A[4][0] = U_k*((gamma-2.)*vSquHlf-c*c*inv_kap_m1);
  A[4][1] = c*c*inv_kap_m1*nV[0]-kapm1*vel[0]*(nVel[1]+nVel[2])-(kapm1*vel[0]*vel[0]-vSquHlf)*nV[0];
  A[4][2] = c*c*inv_kap_m1*nV[1]-kapm1*vel[1]*(nVel[0]+nVel[2])-(kapm1*vel[1]*vel[1]-vSquHlf)*nV[1];
  A[4][3] = c*c*inv_kap_m1*nV[2]-kapm1*vel[2]*(nVel[0]+nVel[1])-(kapm1*vel[2]*vel[2]-vSquHlf)*nV[2];
  A[4][4] = gamma*U_k-surfVeloc;
}

//--------------------------------------------------------------------------//
//                                                                          //
//                            dF/dU x dU/dQ                                 //
//           multiplications are only performed on non-zero entries         //
//                                                                          //
//--------------------------------------------------------------------------//
#define dF_dQ_10(A,fU,UQ)  {A[1][0] = fU[1][1]*UQ[1][0] + fU[1][2]*UQ[2][0] + fU[1][3]*UQ[3][0];                       } //df(u_imp)/d(rho)
#define dF_dQ_20(A,fU,UQ)  {A[2][0] = fU[2][1]*UQ[1][0] + fU[2][2]*UQ[2][0] + fU[2][3]*UQ[3][0];                       } //df(v_imp)/d(rho)
#define dF_dQ_30(A,fU,UQ)  {A[3][0] = fU[3][1]*UQ[1][0] + fU[3][2]*UQ[2][0] + fU[3][3]*UQ[3][0];                       } //df(w_imp)/d(rho)
#define dF_dQ_40(A,fU,UQ)  {A[4][0] = fU[4][1]*UQ[1][0] + fU[4][2]*UQ[2][0] + fU[4][3]*UQ[3][0] + fU[4][4]*UQ[4][0];   } //df(energ)/d(rho)

#define dF_dQ_11(A,fU,UQ)  {A[1][1] = fU[1][1]*UQ[1][1];                                                               } //df(u_imp)/d(rho*u)
#define dF_dQ_21(A,fU,UQ)  {A[2][1] = fU[2][1]*UQ[1][1];                                                               } //df(v_imp)/d(rho*u)
#define dF_dQ_31(A,fU,UQ)  {A[3][1] = fU[3][1]*UQ[1][1];                                                               } //df(w_imp)/d(rho*u)
#define dF_dQ_41(A,fU,UQ)  {A[4][1] = fU[4][1]*UQ[1][1]                                         + fU[4][4]*UQ[4][1];   } //df(energ)/d(rho*u)

#define dF_dQ_12(A,fU,UQ)  {A[1][2] =                   fU[1][2]*UQ[2][2];                                             } //df(u_imp)/d(rho*v)
#define dF_dQ_22(A,fU,UQ)  {A[2][2] =                   fU[2][2]*UQ[2][2];                                             } //df(v_imp)/d(rho*v)
#define dF_dQ_32(A,fU,UQ)  {A[3][2] =                   fU[3][2]*UQ[2][2];                                             } //df(w_imp)/d(rho*v)
#define dF_dQ_42(A,fU,UQ)  {A[4][2] =                   fU[4][2]*UQ[2][2]                       + fU[4][4]*UQ[4][2];   } //df(energ)/d(rho*v)

#define dF_dQ_13(A,fU,UQ)  {A[1][3] =                                       fU[1][3]*UQ[3][3];                         } //df(u_imp)/d(rho*w)
#define dF_dQ_23(A,fU,UQ)  {A[2][3] =                                       fU[2][3]*UQ[3][3];                         } //df(v_imp)/d(rho*w)
#define dF_dQ_33(A,fU,UQ)  {A[3][3] =                                       fU[3][3]*UQ[3][3];                         } //df(w_imp)/d(rho*w)
#define dF_dQ_43(A,fU,UQ)  {A[4][3] =                                       fU[4][3]*UQ[3][3]   + fU[4][4]*UQ[4][3];   } //df(energ)/d(rho*w)

#define dF_dQ_44(A,fU,UQ)  {A[4][4] = fU[4][1]*UQ[1][4] + fU[4][2]*UQ[2][4] + fU[4][3]*UQ[3][4] + fU[4][4]*UQ[4][4];   } //df(energ)/d(energ)

#define dF_dQ_0(A,fU,UQ) { dF_dQ_10(A,fU,UQ)  dF_dQ_20(A,fU,UQ) dF_dQ_30(A,fU,UQ) dF_dQ_40(A,fU,UQ) }
#define dF_dQ_1(A,fU,UQ) { dF_dQ_11(A,fU,UQ)  dF_dQ_21(A,fU,UQ) dF_dQ_31(A,fU,UQ) dF_dQ_41(A,fU,UQ) }
#define dF_dQ_2(A,fU,UQ) { dF_dQ_12(A,fU,UQ)  dF_dQ_22(A,fU,UQ) dF_dQ_32(A,fU,UQ) dF_dQ_42(A,fU,UQ) }
#define dF_dQ_3(A,fU,UQ) { dF_dQ_13(A,fU,UQ)  dF_dQ_23(A,fU,UQ) dF_dQ_33(A,fU,UQ) dF_dQ_43(A,fU,UQ) }
#define dF_dQ_4(A,fU,UQ) {                                                        dF_dQ_44(A,fU,UQ) }

#define dFdU_times_dUdQ(A,fU,UQ)   { dF_dQ_0(A,fU,UQ) dF_dQ_1(A,fU,UQ) dF_dQ_2(A,fU,UQ) dF_dQ_3(A,fU,UQ) dF_dQ_4(A,fU,UQ) }



/**
 *   calc viscous flux...
 *
 * ------------------------------------------------------------------------------------------
 *
 *   explicit part
 *
 *   viscous flux in RHS form is....
 *      Frhou_i = area*mu*( dui/dxj + duj/dxi - 2/3*deltaij*duk/dxk )*nj
 *              = area*( tauij_nj )
 *
 *   energy flux
 *      FrhoE = area*( k*dT/dj*nj + ui*tauij_nj )
 *
 * ------------------------------------------------------------------------------------------
 *
 *   implicit part
 *
 *   calculates viscous Jacobeans
 *   first calculate dFlux/dU, with U=[rho,u,v,w,T],
 *   then calculate transformation matrix dU/dQ, with Q=[rho,rhou,rhov,rhow,e]
 *   as given in Pulliam, AIAA-85-0360.
 *   matrix multiplication dFdU_times_dUdQ gives Av = dFlux/dU x dU/dQ
 *
 *   Jacobians for implicit part are obtained with simplified viscous fluxes - terms with alfa
 *   are omitted. coefficients are taken from the least square gradient
 */
void UgpWithCvCompFlow::addViscFlux(double *Frhou, double &FrhoE, double (*A0)[5], double (*A1)[5],
    const double rho0, const double *u0, const double (&grad_u0)[3][3], const double h0, const double *grad_h0, const double T0, const double R0, const double gam0, const double kine0,
    const double rho1, const double *u1, const double (&grad_u1)[3][3], const double h1, const double *grad_h1, const double T1, const double R1, const double gam1, const double kine1,
    const double mul, const double mut, const double lambdaOverCp, const double kine_fa, const double *u_fa, 
    const double area, const double *nVec, const double smag, const double *sVec)
{
  double alpha = vecDotVec3d(nVec, sVec);
  assert((alpha > 0.0) && (alpha < 1.000000001));

  double grad_u_f[3][3];
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      grad_u_f[i][j] = 0.5*(grad_u0[i][j] + grad_u1[i][j]);


  double fCorr[3] = { nVec[0] - alpha*sVec[0], 
                      nVec[1] - alpha*sVec[1],
                      nVec[2] - alpha*sVec[2]};

  /*double fCorr[3] = { 0.0, 0.0, 0.0};
  alpha = 1.0;*/

  
  // ========================================================================
  // momentum equation
  // ========================================================================

  double muTotalMomentum = mul + mut;
  double tauij_nj[3];

  for (int i = 0; i < 3; i++)
  {
    tauij_nj[i] = muTotalMomentum * (alpha * (u1[i] - u0[i]) / smag
                      + grad_u_f[i][0] * fCorr[0] + grad_u_f[0][i] * nVec[0]
                      + grad_u_f[i][1] * fCorr[1] + grad_u_f[1][i] * nVec[1]
                      + grad_u_f[i][2] * fCorr[2] + grad_u_f[2][i] * nVec[2]);
  }

  // viscosity times trace of strain rate tensor times 2/3...
  double tmp = 2.0 / 3.0 * muTotalMomentum * (grad_u_f[0][0] + grad_u_f[1][1] + grad_u_f[2][2]);
  if (turbModel > NONE)
    tmp += 1.0 / 3.0 * (rho0 + rho1) * kine_fa;  // and 2/3*rho*kine if turb model is on

  tauij_nj[0] -= tmp * nVec[0];
  tauij_nj[1] -= tmp * nVec[1];
  tauij_nj[2] -= tmp * nVec[2];

  // subtract from momentum flux (require LHS form - see convective term above)...
  Frhou[0] = -area * tauij_nj[0];
  Frhou[1] = -area * tauij_nj[1];
  Frhou[2] = -area * tauij_nj[2];


  // ========================================================================
  // energy equation
  // ========================================================================

  double keff = lambdaOverCp + mut / PrTurb;
  double enth = keff * (alpha * (h1 - h0) / smag
                + 0.5 * ((grad_h0[0] + grad_h1[0]) * fCorr[0]
                        +(grad_h0[1] + grad_h1[1]) * fCorr[1]
                        +(grad_h0[2] + grad_h1[2]) * fCorr[2]));

  // model for triple correlation
  double psi = 0.0;
  if (turbModel > NONE)
    psi = muTotalMomentum * (alpha * (kine1 - kine0) / smag);

  FrhoE = -area * (enth + tauij_nj[0] * u_fa[0] + tauij_nj[1] * u_fa[1] + tauij_nj[2] * u_fa[2] + psi);

  // ========================================================================
  // implicit part
  // ========================================================================
  double dfdU[5][5] = {0., 0., 0., 0., 0.,  0., 0., 0., 0., 0.,  0., 0., 0., 0., 0.,  0., 0., 0., 0., 0.,  0., 0., 0., 0., 0.};
  double dUdQ[5][5] = {0., 0., 0., 0., 0.,  0., 0., 0., 0., 0.,  0., 0., 0., 0., 0.,  0., 0., 0., 0., 0.,  0., 0., 0., 0., 0.};

  // ========================================================================
  // calculate A0
  // ========================================================================
  if (A0 != NULL)
  {
    // d(fx)/duvw
    dfdU[1][1] = - muTotalMomentum * alpha / smag * area;
    dfdU[1][2] = 0.0;
    dfdU[1][3] = 0.0;
    // d(fy)/duvw
    dfdU[2][1] = 0.0;
    dfdU[2][2] = - muTotalMomentum * alpha / smag * area;
    dfdU[2][3] = 0.0;
    // d(fz)/duvw
    dfdU[3][1] = 0.0;
    dfdU[3][2] = 0.0;
    dfdU[3][3] = - muTotalMomentum * alpha / smag * area;
  
    // d(fE)/duvwTemp
    dfdU[4][1] = area * tauij_nj[0] * 0.5 + u_fa[0] * dfdU[1][1];
    dfdU[4][2] = area * tauij_nj[1] * 0.5 + u_fa[1] * dfdU[2][2];
    dfdU[4][3] = area * tauij_nj[2] * 0.5 + u_fa[2] * dfdU[3][3];
  
    dfdU[4][4] = -area * keff * alpha / smag;
  
    // define dU/dQ on the left side, with: U=[rho, u_i, h], Q=[rho, rhou_i, rhoE]
    double invRho = 1.0 / rho0;

    dUdQ[1][0] = -u0[0] * invRho;
    dUdQ[2][0] = -u0[1] * invRho;
    dUdQ[3][0] = -u0[2] * invRho;
    dUdQ[4][0] =  gam0 * invRho * (- h0 - kine0 + R0*T0 + 0.5 * vecDotVec3d(u0,u0));

    dUdQ[1][1] =  invRho;
    dUdQ[4][1] = -gam0 * u0[0] * invRho;

    dUdQ[2][2] =  invRho;
    dUdQ[4][2] = -gam0 * u0[1] * invRho;

    dUdQ[3][3] =  invRho;
    dUdQ[4][3] = -gam0 * u0[2] * invRho;

    dUdQ[4][4] =  gam0 * invRho;
    
    dFdU_times_dUdQ(A0, dfdU, dUdQ);
  }

  // ========================================================================
  // calculate A1
  // ========================================================================
  if (A1 != NULL)
  {
    dfdU[1][1] = muTotalMomentum * alpha / smag * area;
    dfdU[1][2] = 0.0;
    dfdU[1][3] = 0.0;
  
    dfdU[2][1] = 0.0;
    dfdU[2][2] = muTotalMomentum * alpha / smag * area;
    dfdU[2][3] = 0.0;
  
    dfdU[3][1] = 0.0;
    dfdU[3][2] = 0.0;
    dfdU[3][3] = muTotalMomentum * alpha / smag * area;

    dfdU[4][1] = area * tauij_nj[0] * 0.5 + u_fa[0] * dfdU[1][1];
    dfdU[4][2] = area * tauij_nj[1] * 0.5 + u_fa[1] * dfdU[2][2];
    dfdU[4][3] = area * tauij_nj[2] * 0.5 + u_fa[2] * dfdU[3][3];
  
    dfdU[4][4] = area * keff * alpha / smag;
  
    // define dU/dQ on the right side, with: U=[rho, u_i, h], Q=[rho, rhou_i, rhoE]
    double invRho = 1.0 / rho1;

    dUdQ[1][0] = -u1[0] * invRho;
    dUdQ[2][0] = -u1[1] * invRho;
    dUdQ[3][0] = -u1[2] * invRho;
    dUdQ[4][0] =  gam1 * invRho * (- h1 - kine1 + R1 * T1 + 0.5 * vecDotVec3d(u1,u1));

    dUdQ[1][1] =  invRho;
    dUdQ[4][1] = -gam1 * u1[0] * invRho;

    dUdQ[2][2] =  invRho;
    dUdQ[4][2] = -gam1 * u1[1] * invRho;

    dUdQ[3][3] =  invRho;
    dUdQ[4][3] = -gam1 * u1[2] * invRho;

    dUdQ[4][4] =  gam1 * invRho;

    dFdU_times_dUdQ(A1, dfdU, dUdQ);
  }
}


//
//#undef dU_dQ_0
//#undef dU_dQ_1
//#undef dU_dQ_2
//#undef dU_dQ_3
//#undef dU_dQ_4
//#undef dU_dQ
//
//#undef dF_dQ_10
//#undef dF_dQ_20
//#undef dF_dQ_30
//#undef dF_dQ_40
//
//#undef dF_dQ_11
//#undef dF_dQ_21
//#undef dF_dQ_31
//#undef dF_dQ_41
//
//#undef dF_dQ_12
//#undef dF_dQ_22
//#undef dF_dQ_32
//#undef dF_dQ_42
//
//#undef dF_dQ_13
//#undef dF_dQ_23
//#undef dF_dQ_33
//#undef dF_dQ_43
//#undef dF_dQ_44
//
//#undef dF_dQ_0
//#undef dF_dQ_1
//#undef dF_dQ_2
//#undef dF_dQ_3
//#undef dF_dQ_4
//
//#undef dFdU_times_dUdQ


void UgpWithCvCompFlow::calcEulerFluxCoupled_HLLC(double *EulerFlux,
         const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *ScalL, const double kL,
         const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *ScalR, const double kR,
         const double area, const double *nVec, const int nScal, const double *ConvTerm, const double surfVeloc, const string BC_treatment)
{

  if (BC_treatment == "ALL_TERMS")
  {
    double unL   = vecDotVec3d(uL, nVec);
    double uLuL  = vecDotVec3d(uL, uL);
    double cL    = sqrt(gammaL * pL / rhoL);
    double rhoeL = (h0 + 0.5 * uLuL + kL) * rhoL - pL;
  
    double unR   = vecDotVec3d(uR, nVec);
    double uRuR  = vecDotVec3d(uR, uR);
    double cR    = sqrt(gammaR * pR / rhoR);
    double rhoeR = (h1 + 0.5 * uRuR + kR) * rhoR - pR;
  
  
    // Roe's aveaging
    double Rrho = sqrt(rhoR / rhoL);
    double tmp = 1.0 / (1.0 + Rrho);
    double velRoe[3];
    for (int i = 0; i < 3; i++)
      velRoe[i] = tmp * (uL[i] + uR[i] * Rrho);
    double uRoe  = vecDotVec3d(velRoe, nVec);
  
    //    double cRoe  = sqrt((gammaL-1.0)*(hRoe- 0.5*vecDotVec3d(velRoe, velRoe)));
    double gamPdivRho = tmp * ((gammaL * pL / rhoL + 0.5 * (gammaL - 1.0) * uLuL) + (gammaR * pR / rhoR + 0.5 * (gammaR - 1.0) * uRuR) * Rrho);
    double cRoe  = sqrt(gamPdivRho - ((gammaL + gammaR) * 0.5 - 1.0) * 0.5 * vecDotVec3d(velRoe, velRoe));
  
    // speed of sound at L and R
    double sL = min(uRoe - cRoe, unL - cL);
    double sR = max(uRoe + cRoe, unR + cR);
  
    // speed of contact surface
    double sM = (pL - pR - rhoL * unL * (sL - unL) + rhoR * unR * (sR - unR)) / (rhoR * (sR - unR) - rhoL * (sL - unL));
  
    // pressure at right and left (pR=pL) side of contact surface
    double pStar = rhoR * (unR - sR) * (unR - sM) + pR;
  
    if (sM >= 0.0)
    {
      if (sL > 0.0)
      {
        EulerFlux[0] = rhoL * unL;
        for (int i = 0; i < 3; i++)
          EulerFlux[1+i] = rhoL * uL[i] * unL + pL * nVec[i];
        EulerFlux[4] = (rhoeL + pL) * unL;
        for (int iScal = 0; iScal < nScal; iScal++)
          EulerFlux[5+iScal] = ConvTerm[iScal] * EulerFlux[0] * ScalL[iScal];                      // Flux = 0 if ConvTerm = 0 (no advection)
      }
      else
      {
        double invSLmSs = 1.0 / (sL - sM);
        double sLmuL = sL - unL;
        double rhoSL = rhoL * sLmuL * invSLmSs;
        double rhouSL[3];
        for (int i = 0; i < 3; i++)
          rhouSL[i] = (rhoL * uL[i] * sLmuL + (pStar - pL) * nVec[i]) * invSLmSs;
        double rhoeSL = (sLmuL * rhoeL - pL * unL + pStar * sM) * invSLmSs;
  
        EulerFlux[0] = rhoSL * sM;
        for (int i = 0; i < 3; i++)
          EulerFlux[1+i] = rhouSL[i] * sM + pStar * nVec[i];
        EulerFlux[4] = (rhoeSL + pStar) * sM;
        for (int iScal = 0; iScal < nScal; iScal++)
          EulerFlux[5+iScal] = ConvTerm[iScal] * EulerFlux[0] * ScalL[iScal];
      }
    }
    else
    {
      if (sR >= 0.0)
      {
        double invSRmSs = 1.0 / (sR - sM);
        double sRmuR = sR - unR;
        double rhoSR = rhoR * sRmuR * invSRmSs;
        double rhouSR[3];
        for (int i = 0; i < 3; i++)
          rhouSR[i] = (rhoR * uR[i] * sRmuR + (pStar - pR) * nVec[i]) * invSRmSs;
        double rhoeSR = (sRmuR * rhoeR - pR * unR + pStar * sM) * invSRmSs;
  
        EulerFlux[0] = rhoSR * sM;
        for (int i = 0; i < 3; i++)
          EulerFlux[1+i] = rhouSR[i] * sM + pStar * nVec[i];
        EulerFlux[4] = (rhoeSR + pStar) * sM;
        for (int iScal = 0; iScal < nScal; iScal++)
          EulerFlux[5+iScal] = ConvTerm[iScal] * EulerFlux[0] * ScalR[iScal];
      }
      else
      {
        EulerFlux[0] = rhoR * unR;
        for (int i = 0; i < 3; i++)
          EulerFlux[1+i] = rhoR * uR[i] * unR + pR * nVec[i];
        EulerFlux[4] = (rhoeR + pR) * unR;
        for (int iScal = 0; iScal < nScal; iScal++)
          EulerFlux[5+iScal] = ConvTerm[iScal] * EulerFlux[0] * ScalR[iScal];
      }
    }
  }
  
  // Only pressure term present in Euler Flux (e.g., symmetry or wall boundary condition); assumes that phiL=phiR -> SM=0.
  else if (BC_treatment == "ONLY_PRESSURE")
  {
    EulerFlux[0] = 0.0;
    for (int i = 1; i < 4; i++)
      EulerFlux[i] = pL * nVec[i-1];
    for (int i = 4; i < 5+nScal; i++)
      EulerFlux[i] = 0.0;
  }
  
  else
  {
    cerr << "### BC treatment flag '" << BC_treatment << "' not recognized for the Euler flux! ###" << endl;
    throw(-1);    
  }

  for (int i = 0; i < 5+nScal; i++)
    EulerFlux[i] *= area;
}


void UgpWithCvCompFlow::calcEulerFluxMatricesCoupled_HLLC(double **A_L, double **A_R,
         const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *scalL, const double *dpress_dscalL, const double kL,
         const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *scalR, const double *dpress_dscalR, const double kR,
         const double area, const double *nVec, const int nScal, const double *ConvTerm, const double surfVeloc, const string BC_treatment)
{
  
  // A_L = dF / dQ_L = (dF / dU_L) * (dU_L / dQ_L)
  // A_R = dF / dQ_R = (dF / dU_R) * (dU_R / dQ_R)
  // Q = (rho, rho*u, rho*E, rho*Z) and U = (rho, rho*u, rho*E, Z)
  // Only the scalars are different between Q and U


  if (BC_treatment == "ALL_TERMS")
  {
    double *dSMdUL = new double[5+nScal]; 
    double *dSMdUR = new double[5+nScal];
    
    // First compute dF / dU_L and dF / dU_R and save them in A_L and A_R
    double unL   = vecDotVec3d(uL, nVec);
    double uLuL  = vecDotVec3d(uL, uL);
    double cL    = sqrt(gammaL * pL / rhoL);
    double hL    = h0 + kL;
    double rhoeL = (hL + 0.5 * uLuL) * rhoL - pL;
  
    double unR   = vecDotVec3d(uR, nVec);
    double uRuR  = vecDotVec3d(uR, uR);
    double cR    = sqrt(gammaR * pR / rhoR);
    double hR    = h1 + kR;
    double rhoeR = (hR + 0.5 * uRuR) * rhoR - pR;
  
  
    // Roe's aveaging
    double Rrho = sqrt(rhoR / rhoL);
    double tmp = 1.0 / (1.0 + Rrho);
    double velRoe[3];
    for (int i = 0; i < 3; i++)
      velRoe[i] = tmp * (uL[i] + uR[i] * Rrho);
    double uRoe  = vecDotVec3d(velRoe, nVec);
  
    double gamPdivRho = tmp * ((gammaL * pL / rhoL + 0.5 * (gammaL - 1.0) * uLuL) + (gammaR * pR / rhoR + 0.5 * (gammaR - 1.0) * uRuR) * Rrho);
    double cRoe  = sqrt(gamPdivRho - ((gammaL + gammaR) * 0.5 - 1.0) * 0.5 * vecDotVec3d(velRoe, velRoe));
  
    // speed of sound at L and R
    double sL = min(uRoe - cRoe, unL - cL);
    double sR = max(uRoe + cRoe, unR + cR);
  
    // speed of contact surface
    double invrhotld = 1.0 / (rhoR * (sR - unR) - rhoL * (sL - unL));
    double sM = (pL - pR - rhoL * unL * (sL - unL) + rhoR * unR * (sR - unR)) / (rhoR * (sR - unR) - rhoL * (sL - unL));
  
    // pressure at right and left (pR=pL) side of contact surface
    double pStar = rhoR * (unR - sR) * (unR - sM) + pR;
  
    // normal vector scaled by area
    double nVecArea[3];
    for (int i = 0; i < 3; i++)        
      nVecArea[i] = nVec[i] * area;
    
    if (sM >= 0.0)
    {
      if (sL > 0.0)
      {
        if (A_L != NULL)
          calcJacobianCoupled_HLLC(A_L, uL, pL, rhoL, hL, scalL, dpress_dscalL, nVecArea, gammaL, 0.0, nScal);
        
        if (A_R != NULL)
          for (int i = 0; i < 5+nScal; i++)
            for (int j = 0; j < 5+nScal; j++)
              A_R[i][j] = 0.0;
      }
      else
      {
        double OmegaL = 1.0 / (sL - sM);
        double sLmuL = sL - unL;
        double rhoSL = rhoL * sLmuL * OmegaL;
        double rhouSL[3];
        for (int i = 0; i < 3; i++)
          rhouSL[i] = (rhoL * uL[i] * sLmuL + (pStar - pL) * nVec[i]) * OmegaL;
        double rhoeSL = (sLmuL * rhoeL - pL * unL + pStar * sM) * OmegaL;
        
        if (A_L != NULL)
        {
          calcdSMdU_HLLC(dSMdUL, uL, pL, rhoL, hL, dpress_dscalL, nVec, gammaL, sM, sL, invrhotld, 1.0, nScal);
          calcJacobianCoupled_HLLC(A_L, uL, pL, rhoL, hL, scalL, dpress_dscalL, nVecArea, gammaL, 0.0, nScal);
          double rhoSmqR = rhoR * (sR - unR);
          calcSubSonicJacobianCoupled_HLLC(A_L, dSMdUL, rhoSL, rhouSL, rhoeSL, scalL, pStar, OmegaL, sM, rhoSmqR, nVec, area, nScal);
          for (int i = 0; i < 5; i++)
            A_L[i][i] += area * OmegaL * sM * sL;
          for (int iScal = 0; iScal < nScal; iScal++)
          {
            A_L[5+iScal][0]       += area * OmegaL * sM * sL * scalL[iScal];
            A_L[5+iScal][5+iScal] += area * OmegaL * sM * sL * rhoL;
          }
        }
  
        if (A_R != NULL)
        {
          calcdSMdU_HLLC(dSMdUR, uR, pR, rhoR, hR, dpress_dscalR, nVec, gammaR, sM, sR, invrhotld, -1.0, nScal);
          for (int i = 0; i < 5+nScal; i++)
            for (int j = 0; j < 5+nScal; j++)
              A_R[i][j] = 0.0;
          double rhoSmqL = rhoL * (sL - unL);
          calcSubSonicJacobianCoupled_HLLC(A_R, dSMdUR, rhoSL, rhouSL, rhoeSL, scalL, pStar, OmegaL, sM, rhoSmqL, nVec, area, nScal);
        }
      }
    }
    else
    {
      if (sR >= 0.0)
      {
        double OmegaR = 1.0 / (sR - sM);
        double sRmuR = sR - unR;
        double rhoSR = rhoR * sRmuR * OmegaR;
        double rhouSR[3];
        for (int i = 0; i < 3; i++)
          rhouSR[i] = (rhoR * uR[i] * sRmuR + (pStar - pR) * nVec[i]) * OmegaR;
        double rhoeSR = (sRmuR * rhoeR - pR * unR + pStar * sM) * OmegaR;
        
        if (A_R != NULL)
        {
          calcdSMdU_HLLC(dSMdUR, uR, pR, rhoR, hR, dpress_dscalR, nVec, gammaR, sM, sR, invrhotld, -1.0, nScal);
          calcJacobianCoupled_HLLC(A_R, uR, pR, rhoR, hR, scalR, dpress_dscalR, nVecArea, gammaR, 0.0, nScal);
          double rhoSmqL = rhoL * (sL - unL);
          calcSubSonicJacobianCoupled_HLLC(A_R, dSMdUR, rhoSR, rhouSR, rhoeSR, scalR, pStar, OmegaR, sM, rhoSmqL, nVec, area, nScal);
          for (int i = 0; i < 5; i++)
            A_R[i][i] += area * OmegaR * sM * sR;
          for (int iScal = 0; iScal < nScal; iScal++)
          {
            A_R[5+iScal][0]       += area * OmegaR * sM * sR * scalR[iScal];
            A_R[5+iScal][5+iScal] += area * OmegaR * sM * sR * rhoR;
          }
        }
  
        if (A_L != NULL)
        {
          calcdSMdU_HLLC(dSMdUL, uL, pL, rhoL, hL, dpress_dscalL, nVec, gammaL, sM, sL, invrhotld, 1.0, nScal);
          for (int i = 0; i < 5+nScal; i++)
            for (int j = 0; j < 5+nScal; j++)
              A_L[i][j] = 0.0;
          double rhoSmqR = rhoR * (sR - unR);
          calcSubSonicJacobianCoupled_HLLC(A_L, dSMdUL, rhoSR, rhouSR, rhoeSR, scalR, pStar, OmegaR, sM, rhoSmqR, nVec, area, nScal);
        }
         
      }
      else
      {
        if (A_R != NULL)
          calcJacobianCoupled_HLLC(A_R, uR, pR, rhoR, hR, scalR, dpress_dscalR, nVecArea, gammaR, 0.0, nScal);
        
        if (A_L != NULL)
          for (int i = 0; i < 5+nScal; i++)
            for (int j = 0; j < 5+nScal; j++)
              A_L[i][j] = 0.0;
      }
    }
    
    if (dSMdUL != NULL)   {delete [] dSMdUL;   dSMdUL = NULL;}
    if (dSMdUR != NULL)   {delete [] dSMdUR;   dSMdUR = NULL;}
  }
  
  // Only pressure term present in Euler Flux (e.g., symmetry or wall boundary condition); assumes that phiL=phiR -> SM=0.
  // A_R = NULL
  else if (BC_treatment == "ONLY_PRESSURE")
  {
    if (A_L != NULL)
    {
      double uLuL  = vecDotVec3d(uL, uL);
      double hL    = h0 + kL;
  
      // No Euler flux for density
      for (int j = 0; j < 5+nScal; j++)
        A_L[0][j] = 0.0;
  
      double dpdrho = 0.5 * (gammaL - 1.0) * uLuL - (gammaL - 1.0) * hL + gammaL * pL / rhoL;
      
      A_L[1][0] = dpdrho *                   nVec[0] * area;
      A_L[1][1] = - (gammaL - 1.0) * uL[0] * nVec[0] * area;
      A_L[1][2] = - (gammaL - 1.0) * uL[1] * nVec[0] * area;
      A_L[1][3] = - (gammaL - 1.0) * uL[2] * nVec[0] * area;
      A_L[1][4] =   (gammaL - 1.0) *         nVec[0] * area;
      for (int iScal = 0; iScal < nScal; iScal++)
        A_L[1][5+iScal] = dpress_dscalL[iScal] * nVec[0] * area;
      
      A_L[2][0] = dpdrho *                   nVec[1] * area;
      A_L[2][1] = - (gammaL - 1.0) * uL[0] * nVec[1] * area;
      A_L[2][2] = - (gammaL - 1.0) * uL[1] * nVec[1] * area;
      A_L[2][3] = - (gammaL - 1.0) * uL[2] * nVec[1] * area;
      A_L[2][4] =   (gammaL - 1.0) *         nVec[1] * area;
      for (int iScal = 0; iScal < nScal; iScal++)
        A_L[2][5+iScal] = dpress_dscalL[iScal] * nVec[1] * area;
      
      A_L[3][0] = dpdrho *                   nVec[2] * area;
      A_L[3][1] = - (gammaL - 1.0) * uL[0] * nVec[2] * area;
      A_L[3][2] = - (gammaL - 1.0) * uL[1] * nVec[2] * area;
      A_L[3][3] = - (gammaL - 1.0) * uL[2] * nVec[2] * area;
      A_L[3][4] =   (gammaL - 1.0) *         nVec[2] * area;
      for (int iScal = 0; iScal < nScal; iScal++)
        A_L[3][5+iScal] = dpress_dscalL[iScal] * nVec[2] * area;
          
      // No Euler flux for energy
      for (int i = 4; i < 5+nScal; i++)
        for (int j = 0; j < 5+nScal; j++)
          A_L[i][j] = 0.0;
    }
    
    if (A_R != NULL)
      for (int i = 0; i < 5+nScal; i++)
        for (int j = 0; j < 5+nScal; j++)
          A_R[i][j] = 0.0;
  }
  
  // Finally, multiply dF / dU_L and dF / dU_R with transformation matrix dU_L / dQ_L and dU_R / dQ_R
  // Also multiply by 0 if no convective term
  // Only the scalars are concerned!
  if (A_L != NULL)
  {
    for (int iScal = 0; iScal < nScal; iScal++)
      for (int j = 0; j < 5+nScal; j++)
         A_L[5+iScal][j] *= ConvTerm[iScal];

    for (int j = 0; j < 5+nScal; j++)
      for (int iScal = 0; iScal < nScal; iScal++)
      {
        A_L[j][0] -= scalL[iScal] / rhoL * A_L[j][5+iScal];
        A_L[j][5+iScal] /= rhoL;
      }
  }
  
  if (A_R != NULL)
  {
    for (int iScal = 0; iScal < nScal; iScal++)
      for (int j = 0; j < 5+nScal; j++)
        A_R[5+iScal][j] *= ConvTerm[iScal];

    for (int j = 0; j < 5+nScal; j++)
      for (int iScal = 0; iScal < nScal; iScal++)
      {
        A_R[j][0] -= scalR[iScal] / rhoR * A_R[j][5+iScal];
        A_R[j][5+iScal] /= rhoR;
      }
  }
}

/**
 * supersonic jacobian for Euler flux, coupled NSE and scalars 
 */  
void UgpWithCvCompFlow::calcJacobianCoupled_HLLC(double **A, const double *vel, const double pp, const double rrho, const double hh, const double *scal, const double *dpress_dscal, 
                                                 const double *nV, const double gamma, const double surfVeloc, const int nScal) // nV is not normalized
{
  double kapm1 = (gamma - 1.0);
  double nVel[3];
  nVel[0] = vel[0] * nV[0];
  nVel[1] = vel[1] * nV[1];
  nVel[2] = vel[2] * nV[2];
  double U_k = nVel[0] + nVel[1] + nVel[2];
  double vSquHlf = 0.5 * vecDotVec3d(vel, vel);

  A[0][0] = 0.0;
  A[0][1] = nV[0];
  A[0][2] = nV[1];
  A[0][3] = nV[2];
  A[0][4] = 0.0;
  for (int iScal = 0; iScal < nScal; iScal++)
    A[0][5+iScal] = 0.0;

  A[1][0] = - vel[0] * U_k + nV[0] * (kapm1 * vSquHlf + gamma * pp / rrho - kapm1 * hh);   //Changed here for new energy definition
  A[1][1] = (2.0 - gamma) * nVel[0] + U_k;
  A[1][2] = vel[0] * nV[1] - kapm1 * vel[1] * nV[0];
  A[1][3] = vel[0] * nV[2] - kapm1 * vel[2] * nV[0];
  A[1][4] = kapm1 * nV[0];
  for (int iScal = 0; iScal < nScal; iScal++)
    A[1][5+iScal] = dpress_dscal[iScal] * nV[0];

  A[2][0] = - vel[1] * U_k + nV[1] * (kapm1 * vSquHlf + gamma * pp / rrho - kapm1 * hh);   //Changed here for new energy definition
  A[2][1] = vel[1] * nV[0] - kapm1 * vel[0] * nV[1];
  A[2][2] = (2.0 - gamma) * nVel[1] + U_k;
  A[2][3] = vel[1] * nV[2] - kapm1 * vel[2] * nV[1];
  A[2][4] = kapm1 * nV[1];
  for (int iScal = 0; iScal < nScal; iScal++)
    A[2][5+iScal] = dpress_dscal[iScal] * nV[1];

  A[3][0] = - vel[2] * U_k + nV[2] * (kapm1 * vSquHlf + gamma * pp / rrho - kapm1 * hh);   //Changed here for new energy definition
  A[3][1] = vel[2] * nV[0] - kapm1 * vel[0] * nV[2];
  A[3][2] = vel[2] * nV[1] - kapm1 * vel[1] * nV[2];
  A[3][3] = (2.0 - gamma) * nVel[2] + U_k;
  A[3][4] = kapm1 * nV[2];
  for (int iScal = 0; iScal < nScal; iScal++)
    A[3][5+iScal] = dpress_dscal[iScal] * nV[2];

  A[4][0] = U_k * ((gamma - 2.0) * vSquHlf - gamma * (hh - pp / rrho));                    //Changed here for new energy definition
  A[4][1] = (hh + vSquHlf) * nV[0] - kapm1 * vel[0] * U_k;                                 //Changed here for new energy definition
  A[4][2] = (hh + vSquHlf) * nV[1] - kapm1 * vel[1] * U_k;                                 //Changed here for new energy definition
  A[4][3] = (hh + vSquHlf) * nV[2] - kapm1 * vel[2] * U_k;                                 //Changed here for new energy definition
  A[4][4] = gamma * U_k;
  for (int iScal = 0; iScal < nScal; iScal++)
    A[4][5+iScal] = dpress_dscal[iScal] * U_k;
  
  for (int iScal = 0; iScal < nScal; iScal++)
  {
    A[5+iScal][0] = 0.0;
    A[5+iScal][1] = scal[iScal] * nV[0];
    A[5+iScal][2] = scal[iScal] * nV[1];
    A[5+iScal][3] = scal[iScal] * nV[2];
    A[5+iScal][4] = 0.0;
    for (int jScal = 0; jScal < nScal; jScal++)
      A[5+iScal][5+jScal] = 0.0;
    A[5+iScal][5+iScal] = rrho * U_k;
  }
}

/**
 *  derivative of contact surface speed
 */  
void UgpWithCvCompFlow::calcdSMdU_HLLC(double *dSMdU, const double *vel, const double pp, const double rrho, const double hh, const double *dpress_dscal,
                                       const double *nV, const double gamma, const double SM, const double S, const double invertrho, const double Factor, const int nScal)            // nV is not normalized
{
  double kapm1 = (gamma - 1.0);
  double nVel[3];
  nVel[0] = vel[0] * nV[0];
  nVel[1] = vel[1] * nV[1];
  nVel[2] = vel[2] * nV[2];
  double U_k = nVel[0] + nVel[1] + nVel[2];
  double vSquHlf = 0.5 * vecDotVec3d(vel, vel);
  
  dSMdU[0] = Factor * invertrho * (-U_k*U_k + kapm1*vSquHlf + gamma*pp/rrho - kapm1*hh + SM*S);
  dSMdU[1] = Factor * invertrho * ((2*U_k - S - SM)*nV[0] - kapm1*vel[0]);
  dSMdU[2] = Factor * invertrho * ((2*U_k - S - SM)*nV[1] - kapm1*vel[1]);
  dSMdU[3] = Factor * invertrho * ((2*U_k - S - SM)*nV[2] - kapm1*vel[2]);
  dSMdU[4] = Factor * invertrho * kapm1;
  for (int iScal = 0; iScal < nScal; iScal++)
    dSMdU[5+iScal] = Factor * invertrho * dpress_dscal[iScal];
}

/**
 *  subsonic jacobian for Euler flux, coupled NSE and scalars
 */
void UgpWithCvCompFlow::calcSubSonicJacobianCoupled_HLLC(double **A, const double *dSMdU, const double rhoS, const double *rhouS, const double rhoeS, const double *scal,
                                                         const double pS, const double Omega, const double sM, const double rhoSmq, const double *nV, const double area, const int nScal)
{
  double OmSM = Omega * sM;
  double OmSMp1 = OmSM + 1.0;
  
  for (int i = 0; i < 5+nScal; i++)
  {
    A[0][i] = -OmSM * A[0][i] + area * OmSMp1 * dSMdU[i] * rhoS;
    A[1][i] = -OmSM * A[1][i] + area * OmSMp1 * dSMdU[i] * (rhouS[0] + rhoSmq * nV[0]);
    A[2][i] = -OmSM * A[2][i] + area * OmSMp1 * dSMdU[i] * (rhouS[1] + rhoSmq * nV[1]);
    A[3][i] = -OmSM * A[3][i] + area * OmSMp1 * dSMdU[i] * (rhouS[2] + rhoSmq * nV[2]);
    A[4][i] = -OmSM * A[4][i] + area * OmSMp1 * dSMdU[i] * (rhoeS + pS + rhoSmq * sM);
    for (int iScal = 0; iScal < nScal; iScal++)
      A[5+iScal][i] = -OmSM * A[5+iScal][i] + area * OmSMp1 * dSMdU[i] * rhoS * scal[iScal];
  }
}

/**
 *  viscous flux and jacobian, coupled NSE and scalars
 */
void UgpWithCvCompFlow::calcViscousFluxCoupled(double *ViscousFlux, double **A0, double **A1,
         const double rho0, const double *u0, const double (&grad_u0)[3][3], const double h0, const double *grad_h0, const double T0, const double R0, const double gam0, const double *Scal0, const double (*gradScal0)[3], const double *dpress_dscal0, const double kine0,
         const double rho1, const double *u1, const double (&grad_u1)[3][3], const double h1, const double *grad_h1, const double T1, const double R1, const double gam1, const double *Scal1, const double (*gradScal1)[3], const double *dpress_dscal1, const double kine1,
         const double mul, const double mut, const double lambdaOverCp, const double kine_fa, const double *u_fa, const double *diff, const double *DiffTerm,
         const double area, const double *nVec, const double smag, const double *sVec, const double alpha, const int nScal)
{
  double grad_u_f[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      grad_u_f[i][j] = 0.5 * (grad_u0[i][j] + grad_u1[i][j]);


  double fCorr[3] = {nVec[0] - alpha * sVec[0], 
                     nVec[1] - alpha * sVec[1],
                     nVec[2] - alpha * sVec[2]};

  
  // ************************************************************************
  // FLUXES
  // ************************************************************************
  
  // ========================================================================
  // momentum equation
  // ========================================================================

  double muTotalMomentum = mul + mut;
  double tauij_nj[3];

  // no flux for density
  ViscousFlux[0] = 0.0;
  
  for (int i = 0; i < 3; i++)
    tauij_nj[i] = muTotalMomentum * (alpha * (u1[i] - u0[i]) / smag
                      + grad_u_f[i][0] * fCorr[0] + grad_u_f[0][i] * nVec[0]
                      + grad_u_f[i][1] * fCorr[1] + grad_u_f[1][i] * nVec[1]
                      + grad_u_f[i][2] * fCorr[2] + grad_u_f[2][i] * nVec[2]);

  // viscosity times trace of strain rate tensor times 2/3... 
  double tmp = 2.0 / 3.0 * muTotalMomentum * (grad_u_f[0][0] + grad_u_f[1][1] + grad_u_f[2][2]); 
  if (turbModel > NONE)
    tmp += 1.0 / 3.0 * (rho0 + rho1) * kine_fa;  // and 2/3*rho*kine if turb model is on

  tauij_nj[0] -= tmp * nVec[0];
  tauij_nj[1] -= tmp * nVec[1];
  tauij_nj[2] -= tmp * nVec[2];

  // subtract from momentum flux (require LHS form - see convective term above)...
  ViscousFlux[1] = - area * tauij_nj[0];
  ViscousFlux[2] = - area * tauij_nj[1];
  ViscousFlux[3] = - area * tauij_nj[2];


  // ========================================================================
  // energy equation
  // ========================================================================

  double keff = lambdaOverCp + mut / PrTurb;
  double enth = keff * (alpha * (h1 - h0) / smag
                + 0.5 * ((grad_h0[0] + grad_h1[0]) * fCorr[0]
                        +(grad_h0[1] + grad_h1[1]) * fCorr[1]
                        +(grad_h0[2] + grad_h1[2]) * fCorr[2]));

  // model for triple correlation
  double psi = 0.0;
  if (turbModel > NONE)
    psi = muTotalMomentum * (alpha * (kine1 - kine0) / smag);


  ViscousFlux[4] = - area * (enth + tauij_nj[0] * u_fa[0] + tauij_nj[1] * u_fa[1] + tauij_nj[2] * u_fa[2] + psi);

  
  // ========================================================================
  // scalars equation
  // ========================================================================
  
  for (int iScal = 0; iScal < nScal; iScal++)
    ViscousFlux[5+iScal] = - area * DiffTerm[iScal] * diff[iScal] * (alpha * (Scal1[iScal] - Scal0[iScal]) / smag
                           + 0.5 * ((gradScal0[iScal][0] + gradScal1[iScal][0]) * fCorr[0]
                           +        (gradScal0[iScal][1] + gradScal1[iScal][1]) * fCorr[1]
                           +        (gradScal0[iScal][2] + gradScal1[iScal][2]) * fCorr[2]));

  
  
  // ************************************************************************
  // JACOBI MATRIX         dF/dU  with U = (rho, u_i, h, Z_j)
  // TRANSFORMATION MATRIX dU/dQ  with Q = (rho, rho u_i, rho_e_t, rho Z_j)
  // IMPLICIT MATRIX       dF/dQ = dF/dU * dU/dQ
  // ************************************************************************

  // To save space and computing time, only non-zero entries
  
  double dFdU_1_1, dFdU_2_2, dFdU_3_3, dFdU_4_1, dFdU_4_2, dFdU_4_3, dFdU_4_4;
  double *dFdU_iScal_iScal = new double[nScal];
  
  double dUdQ_1_0, dUdQ_1_1, dUdQ_2_0, dUdQ_2_2, dUdQ_3_0, dUdQ_3_3;
  double dUdQ_4_0, dUdQ_4_1, dUdQ_4_2, dUdQ_4_3, dUdQ_4_4;
  double *dUdQ_4_iScal     = new double[nScal];
  double *dUdQ_iScal_0     = new double[nScal];
  double *dUdQ_iScal_iScal = new double[nScal];
  
  if (A0 != NULL)
  {
    // dFdU
    // ****
    double area_alpha_smag = area * alpha / smag;
    
    dFdU_1_1 = + muTotalMomentum * area_alpha_smag;
    dFdU_2_2 = + muTotalMomentum * area_alpha_smag;
    dFdU_3_3 = + muTotalMomentum * area_alpha_smag;
    
    dFdU_4_1 = - area * tauij_nj[0] * 0.5 + u_fa[0] * dFdU_1_1;
    dFdU_4_2 = - area * tauij_nj[1] * 0.5 + u_fa[1] * dFdU_2_2;
    dFdU_4_3 = - area * tauij_nj[2] * 0.5 + u_fa[2] * dFdU_3_3;
    dFdU_4_4 = + keff * area_alpha_smag;
    
    for (int iScal = 0; iScal < nScal; iScal++)
      dFdU_iScal_iScal[iScal] = DiffTerm[iScal] * diff[iScal] * area_alpha_smag;
    
    // dUdQ
    // ****
    double invRho = 1.0 / rho0;

    dUdQ_1_0 = - invRho * u0[0];
    dUdQ_1_1 =   invRho;
    dUdQ_2_0 = - invRho * u0[1];
    dUdQ_2_2 =   invRho;
    dUdQ_3_0 = - invRho * u0[2];
    dUdQ_3_3 =   invRho;
    
    dUdQ_4_0 =   invRho * gam0 * (- h0 - kine0 + R0 * T0 + 0.5 * vecDotVec3d(u0,u0));
    dUdQ_4_1 = - invRho * u0[0] * gam0;
    dUdQ_4_2 = - invRho * u0[1] * gam0;
    dUdQ_4_3 = - invRho * u0[2] * gam0;
    dUdQ_4_4 =   invRho * gam0;
    
    for (int iScal = 0; iScal < nScal; iScal++)
    {
      dUdQ_4_iScal[iScal]     =   invRho * invRho * dpress_dscal0[iScal];
      dUdQ_iScal_0[iScal]     = - invRho * Scal0[iScal];
      dUdQ_iScal_iScal[iScal] =   invRho;
    }
    
    // A = dFdU * dUdQ
    // ***************
    A0[1][0] = dFdU_1_1 * dUdQ_1_0;
    A0[1][1] = dFdU_1_1 * dUdQ_1_1;
    A0[2][0] = dFdU_2_2 * dUdQ_2_0;
    A0[2][2] = dFdU_2_2 * dUdQ_2_2;
    A0[3][0] = dFdU_3_3 * dUdQ_3_0;
    A0[3][3] = dFdU_3_3 * dUdQ_3_3;
    
    A0[4][0] = dFdU_4_1 * dUdQ_1_0 + dFdU_4_2 * dUdQ_2_0 + dFdU_4_3 * dUdQ_3_0 + dFdU_4_4 * dUdQ_4_0;
    A0[4][1] = dFdU_4_1 * dUdQ_1_1 + dFdU_4_4 * dUdQ_4_1;
    A0[4][2] = dFdU_4_2 * dUdQ_2_2 + dFdU_4_4 * dUdQ_4_2;
    A0[4][3] = dFdU_4_3 * dUdQ_3_3 + dFdU_4_4 * dUdQ_4_3;
    A0[4][4] = dFdU_4_4 * dUdQ_4_4;
    
    for (int iScal = 0; iScal < nScal; iScal++)
    {
      A0[4][5+iScal]       = dFdU_4_4 * dUdQ_4_iScal[iScal];
      A0[5+iScal][0]       = dFdU_iScal_iScal[iScal] * dUdQ_iScal_0[iScal];
      A0[5+iScal][5+iScal] = dFdU_iScal_iScal[iScal] * dUdQ_iScal_iScal[iScal];
    }
  }
  
  if (A1 != NULL)
  {
    // dFdU
    // ****
    double area_alpha_smag = area * alpha / smag;
    
    dFdU_1_1 = - muTotalMomentum * area_alpha_smag;
    dFdU_2_2 = - muTotalMomentum * area_alpha_smag;
    dFdU_3_3 = - muTotalMomentum * area_alpha_smag;
    
    dFdU_4_1 = - area * tauij_nj[0] * 0.5 + u_fa[0] * dFdU_1_1;
    dFdU_4_2 = - area * tauij_nj[1] * 0.5 + u_fa[1] * dFdU_2_2;
    dFdU_4_3 = - area * tauij_nj[2] * 0.5 + u_fa[2] * dFdU_3_3;
    dFdU_4_4 = - keff * area_alpha_smag;
    
    for (int iScal = 0; iScal < nScal; iScal++)
      dFdU_iScal_iScal[iScal] = - DiffTerm[iScal] * diff[iScal] * area_alpha_smag;
    
    // dUdQ
    // ****
    double invRho = 1.0 / rho1;

    dUdQ_1_0 = - invRho * u1[0];
    dUdQ_1_1 =   invRho;
    dUdQ_2_0 = - invRho * u1[1];
    dUdQ_2_2 =   invRho;
    dUdQ_3_0 = - invRho * u1[2];
    dUdQ_3_3 =   invRho;
    
    dUdQ_4_0 =   invRho * gam1 * (- h1 - kine1 + R1 * T1 + 0.5 * vecDotVec3d(u1,u1));
    dUdQ_4_1 = - invRho * u1[0] * gam1;
    dUdQ_4_2 = - invRho * u1[1] * gam1;
    dUdQ_4_3 = - invRho * u1[2] * gam1;
    dUdQ_4_4 =   invRho * gam1;
    
    for (int iScal = 0; iScal < nScal; iScal++)
    {
      dUdQ_4_iScal[iScal]     =   invRho * invRho * dpress_dscal1[iScal];
      dUdQ_iScal_0[iScal]     = - invRho * Scal1[iScal];
      dUdQ_iScal_iScal[iScal] =   invRho;
    }
    
    // A = dFdU * dUdQ
    // ***************
    A1[1][0] = dFdU_1_1 * dUdQ_1_0;
    A1[1][1] = dFdU_1_1 * dUdQ_1_1;
    A1[2][0] = dFdU_2_2 * dUdQ_2_0;
    A1[2][2] = dFdU_2_2 * dUdQ_2_2;
    A1[3][0] = dFdU_3_3 * dUdQ_3_0;
    A1[3][3] = dFdU_3_3 * dUdQ_3_3;
    
    A1[4][0] = dFdU_4_1 * dUdQ_1_0 + dFdU_4_2 * dUdQ_2_0 + dFdU_4_3 * dUdQ_3_0 + dFdU_4_4 * dUdQ_4_0;
    A1[4][1] = dFdU_4_1 * dUdQ_1_1 + dFdU_4_4 * dUdQ_4_1;
    A1[4][2] = dFdU_4_2 * dUdQ_2_2 + dFdU_4_4 * dUdQ_4_2;
    A1[4][3] = dFdU_4_3 * dUdQ_3_3 + dFdU_4_4 * dUdQ_4_3;
    A1[4][4] = dFdU_4_4 * dUdQ_4_4;
    
    for (int iScal = 0; iScal < nScal; iScal++)
    {
      A1[4][5+iScal]       = dFdU_4_4 * dUdQ_4_iScal[iScal];
      A1[5+iScal][0]       = dFdU_iScal_iScal[iScal] * dUdQ_iScal_0[iScal];
      A1[5+iScal][5+iScal] = dFdU_iScal_iScal[iScal] * dUdQ_iScal_iScal[iScal];
    }
  }

  delete [] dFdU_iScal_iScal;
  delete [] dUdQ_4_iScal;
  delete [] dUdQ_iScal_0;
  delete [] dUdQ_iScal_iScal;
}
