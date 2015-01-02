#include <math.h>
#include "UgpWithCvCompFlowAD.h"

/**
 *
 *  explicit euler flux HLLC
 *
 */
int UgpWithCvCompFlow_AD::calcEulerFlux_HLLC_orig_AD(REALQ &Frho, REALQ *Frhou, REALQ &FrhoE, REALQS *FrhoScal,
        REALQ rhoL, REALQ *uL, REALQ pL, REALQ TL, REALQ h0, REALQ RL, REALQ gammaL, REALQS *ScalL, REALQS kL,
        REALQ rhoR, REALQ *uR, REALQ pR, REALQ TR, REALQ h1, REALQ RR, REALQ gammaR, REALQS *ScalR, REALQS kR,
         double area,  double *nVec, int nScal,  double surfVeloc)
{
  REALQ unL  = vecDotVec3d_AD(uL, nVec);
  REALQ uLuL = vecDotVec3d_AD(uL, uL);
  REALQ cL   = sqrt(gammaL*pL/rhoL);
  REALQ hL   = gammaL/(gammaL-1.0)*pL/rhoL + 0.5*uLuL + kL;
  REALQ eL   = hL*rhoL-pL;

  REALQ unR  = vecDotVec3d_AD(uR, nVec);
  REALQ uRuR = vecDotVec3d_AD(uR, uR);
  REALQ cR   = sqrt(gammaR*pR/rhoR);
  REALQ hR   = gammaR/(gammaR-1.0)*pR/rhoR + 0.5*uRuR + kR;
  REALQ eR   = hR*rhoR-pR;


  // Roe's aveaging
  REALQ Rrho = sqrt(rhoR/rhoL);
  REALQ tmp = 1.0/(1.0+Rrho);
  REALQ velRoe[3];
  for (int i=0; i<3; i++)
    velRoe[i] = tmp*(uL[i] + uR[i]*Rrho);
  REALQ uRoe  = vecDotVec3d_AD(velRoe, nVec);
  REALQ hRoe = tmp*(hL + hR*Rrho);

//  REALQ cRoe  = sqrt((gammaL-1.0)*(hRoe- 0.5*vecDotVec3d_AD(velRoe, velRoe)));
  REALQ gamPdivRho = tmp*((gammaL*pL/rhoL+0.5*(gammaL-1.0)*uLuL) + (gammaR*pR/rhoR+0.5*(gammaR-1.0)*uRuR)*Rrho);
  REALQ cRoe  = sqrt(gamPdivRho - ((gammaL+gammaR)*0.5-1.0)*0.5*vecDotVec3d_AD(velRoe, velRoe));

  // speed of sound at L and R
  REALQ sL = min(uRoe-cRoe, unL-cL);
  REALQ sR = max(uRoe+cRoe, unR+cR);

  // speed of contact surface
  REALQ sM = (pL-pR-rhoL*unL*(sL-unL)+rhoR*unR*(sR-unR))/(rhoR*(sR-unR)-rhoL*(sL-unL));

  // pressure at right and left (pR=pL) side of contact surface
  REALQ pStar = rhoR*(unR-sR)*(unR-sM)+pR;

  if (sM >= 0.0)
  {
    if (sL > 0.0)
    {
      Frho = rhoL*unL;
      for (int i=0; i<3; i++)
        Frhou[i] = rhoL*uL[i]*unL + pL*nVec[i];
      FrhoE = (eL+pL)*unL;
      for (int iScal=0; iScal<nScal; iScal++)
        FrhoScal[iScal] = Frho * ScalL[iScal];            //CHANGE FOR SCALAR AD
    }
    else
    {
      REALQ invSLmSs = 1.0/(sL-sM);
      REALQ sLmuL = sL-unL;
      REALQ rhoSL = rhoL*sLmuL*invSLmSs;
      REALQ rhouSL[3];
      for (int i=0; i<3; i++)
        rhouSL[i] = (rhoL*uL[i]*sLmuL+(pStar-pL)*nVec[i])*invSLmSs;
      REALQ eSL = (sLmuL*eL-pL*unL+pStar*sM)*invSLmSs;

      Frho = rhoSL*sM;
      for (int i=0; i<3; i++)
        Frhou[i] = rhouSL[i]*sM + pStar*nVec[i];
      FrhoE = (eSL+pStar)*sM;
      for (int iScal=0; iScal<nScal; iScal++)
        FrhoScal[iScal] = Frho * ScalL[iScal];		//CHANGE FOR SCALAR AD
    }
  }
  else
  {
    if (sR >= 0.0)
    {
      REALQ invSRmSs = 1.0/(sR-sM);
      REALQ sRmuR = sR-unR;
      REALQ rhoSR = rhoR*sRmuR*invSRmSs;
      REALQ rhouSR[3];
      for (int i=0; i<3; i++)
        rhouSR[i] = (rhoR*uR[i]*sRmuR+(pStar-pR)*nVec[i])*invSRmSs;
      REALQ eSR = (sRmuR*eR-pR*unR+pStar*sM)*invSRmSs;

      Frho = rhoSR*sM;
      for (int i=0; i<3; i++)
        Frhou[i] = rhouSR[i]*sM + pStar*nVec[i];
      FrhoE = (eSR+pStar)*sM;
      for (int iScal=0; iScal<nScal; iScal++)
        FrhoScal[iScal] = Frho * ScalR[iScal];			// CHANGE FOR SCALAR AD
    }
    else
    {
      Frho = rhoR*unR;
      for (int i=0; i<3; i++)
        Frhou[i] = rhoR*uR[i]*unR + pR*nVec[i];
      FrhoE = (eR+pR)*unR;
      for (int iScal=0; iScal<nScal; iScal++)
        FrhoScal[iScal] = Frho * ScalR[iScal];			// CHANGE FOR SCALAR AD
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
   *  explicit euler flux HLLC
   *
   */
  int UgpWithCvCompFlow_AD::calcEulerFlux_HLLC_AD(REALQ &Frho, REALQ *Frhou, REALQ &FrhoE, REALQS *FrhoScal,
       REALQ rhoL,  REALQ *uL,  REALQ pL,  REALQ TL,  REALQ h0,  REALQ RL,  REALQ gammaL,  REALQS *ScalL,  REALQS kineL,
       REALQ rhoR,  REALQ *uR,  REALQ pR,  REALQ TR,  REALQ h1,  REALQ RR,  REALQ gammaR,  REALQS *ScalR,  REALQS kineR,
       double area,  double *nVec,  int nScal,  double surfVeloc)
  {
	  REALQ unL   = vecDotVec3d_AD(uL, nVec);
	  REALQ uLuL  = vecDotVec3d_AD(uL, uL);
	  REALQ cL    = sqrt(gammaL*pL/rhoL);
	  REALQ rhoeL = (h0 + 0.5*uLuL + kineL)*rhoL-pL;

	  REALQ unR   = vecDotVec3d_AD(uR, nVec);
	  REALQ uRuR  = vecDotVec3d_AD(uR, uR);
	  REALQ cR    = sqrt(gammaR*pR/rhoR);
	  REALQ rhoeR = (h1 + 0.5*uRuR + kineR)*rhoR-pR;


	  // Roe's aveaging
	  REALQ Rrho = sqrt(rhoR/rhoL);
	  REALQ tmp = 1.0/(1.0+Rrho);
	  REALQ velRoe[3];
	  for (int i=0; i<3; i++)
		  velRoe[i] = tmp*(uL[i] + uR[i]*Rrho);
	  REALQ uRoe  = vecDotVec3d_AD(velRoe, nVec);

	  //    double cRoe  = sqrt((gammaL-1.0)*(hRoe- 0.5*vecDotVec3d(velRoe, velRoe)));
	  REALQ gamPdivRho = tmp*((gammaL*pL/rhoL+0.5*(gammaL-1.0)*uLuL) + (gammaR*pR/rhoR+0.5*(gammaR-1.0)*uRuR)*Rrho);
	  REALQ cRoe  = sqrt(gamPdivRho - ((gammaL+gammaR)*0.5-1.0)*0.5*vecDotVec3d_AD(velRoe, velRoe));

	  // speed of sound at L and R
	  REALQ sL = min(uRoe-cRoe, unL-cL);
	  REALQ sR = max(uRoe+cRoe, unR+cR);

	  // speed of contact surface
	  REALQ sM = (pL-pR-rhoL*unL*(sL-unL)+rhoR*unR*(sR-unR))/(rhoR*(sR-unR)-rhoL*(sL-unL));

	  // pressure at right and left (pR=pL) side of contact surface
	  REALQ pStar = rhoR*(unR-sR)*(unR-sM)+pR;

#ifdef USE_CONDASSIGN
	  if (sM.value() >= 0.0) {
		  if (sL.value() > 0.0) {
			  Frho = rhoL*unL;
			  for (int i=0; i<3; i++)
				  Frhou[i] = rhoL*uL[i]*unL + pL*nVec[i];
			  FrhoE = (rhoeL+pL)*unL;
			  for (int iScal=0; iScal<nScal; iScal++)
				  FrhoScal[iScal] = Frho * ScalL[iScal];
		  } else {
			  REALQ invSLmSs = 1.0/(sL-sM);
			  REALQ sLmuL = sL-unL;
			  REALQ rhoSL = rhoL*sLmuL*invSLmSs;
			  REALQ rhouSL[3];
			  for (int i=0; i<3; i++)
				  rhouSL[i] = (rhoL*uL[i]*sLmuL+(pStar-pL)*nVec[i])*invSLmSs;
			  REALQ rhoeSL = (sLmuL*rhoeL-pL*unL+pStar*sM)*invSLmSs;

			  Frho = rhoSL*sM;
			  for (int i=0; i<3; i++)
				  Frhou[i] = rhouSL[i]*sM + pStar*nVec[i];
			  FrhoE = (rhoeSL+pStar)*sM;
			  for (int iScal=0; iScal<nScal; iScal++)
				  FrhoScal[iScal] = Frho * ScalL[iScal];
		  }
	  } else {
		  if (sR.value() >= 0.0) {
			  REALQ invSRmSs = 1.0/(sR-sM);
			  REALQ sRmuR = sR-unR;
			  REALQ rhoSR = rhoR*sRmuR*invSRmSs;
			  REALQ rhouSR[3];
			  for (int i=0; i<3; i++)
				  rhouSR[i] = (rhoR*uR[i]*sRmuR+(pStar-pR)*nVec[i])*invSRmSs;
			  REALQ rhoeSR = (sRmuR*rhoeR-pR*unR+pStar*sM)*invSRmSs;

			  Frho = rhoSR*sM;
			  for (int i=0; i<3; i++)
				  Frhou[i] = rhouSR[i]*sM + pStar*nVec[i];
			  FrhoE = (rhoeSR+pStar)*sM;
			  for (int iScal=0; iScal<nScal; iScal++)
				  FrhoScal[iScal] = Frho * ScalR[iScal];
		  } else {
			  Frho = rhoR*unR;
			  for (int i=0; i<3; i++)
				  Frhou[i] = rhoR*uR[i]*unR + pR*nVec[i];
			  FrhoE = (rhoeR+pR)*unR;
			  for (int iScal=0; iScal<nScal; iScal++)
				  FrhoScal[iScal] = Frho * ScalR[iScal];
		  }
	  }
//CASE1
//	  if (sM.value() >= 0.0) {
//		  REALQ invSLmSs = 1.0/(sL-sM);
//		  REALQ sLmuL = sL-unL;
//		  REALQ rhoSL = rhoL*sLmuL*invSLmSs;
//		  REALQ rhouSL[3];
//		  for (int i=0; i<3; i++)
//			  rhouSL[i] = (rhoL*uL[i]*sLmuL+(pStar-pL)*nVec[i])*invSLmSs;
//		  REALQ rhoeSL = (sLmuL*rhoeL-pL*unL+pStar*sM)*invSLmSs;
//
//		  condassign(Frho, sL, rhoL*unL, rhoSL*sM);
//		  for(int i=0; i<3; ++i)
//			  condassign(Frhou[i], sL, rhoL*uL[i]*unL + pL*nVec[i], rhouSL[i]*sM + pStar*nVec[i]);
//		  condassign(FrhoE, sL, (rhoeL+pL)*unL, (rhoeSL+pStar)*sM);
//		  for (int iScal=0; iScal<nScal; iScal++)
//			  FrhoScal[iScal] = Frho * ScalL[iScal];
//	  } else {
//		  REALQ invSRmSs = 1.0/(sR-sM);
//		  REALQ sRmuR = sR-unR;
//		  REALQ rhoSR = rhoR*sRmuR*invSRmSs;
//		  REALQ rhouSR[3];
//		  for (int i=0; i<3; i++)
//			  rhouSR[i] = (rhoR*uR[i]*sRmuR+(pStar-pR)*nVec[i])*invSRmSs;
//		  REALQ rhoeSR = (sRmuR*rhoeR-pR*unR+pStar*sM)*invSRmSs;
//
//		  condassign(Frho, sR+MACHINE_EPS, rhoSR*sM, rhoR*unR);
//		  for(int i=0; i<3; ++i)
//			  condassign(Frhou[i], sR+MACHINE_EPS, rhouSR[i]*sM + pStar*nVec[i], rhoR*uR[i]*unR + pR*nVec[i]);
//		  condassign(FrhoE, sR+MACHINE_EPS, (rhoeSR+pStar)*sM, (rhoeR+pR)*unR);
//		  for(int iScal=0; iScal<nScal; ++iScal)
//			  FrhoScal[iScal] = Frho * ScalR[iScal];
//	  }

#else
	  if (sM >= 0.0)
	  {
		  if (sL > 0.0)
		  {
			  Frho = rhoL*unL;
			  for (int i=0; i<3; i++)
				  Frhou[i] = rhoL*uL[i]*unL + pL*nVec[i];
			  FrhoE = (rhoeL+pL)*unL;
			  for (int iScal=0; iScal<nScal; iScal++)
				  FrhoScal[iScal] = Frho * ScalL[iScal];
		  }
		  else
		  {
			  REALQ invSLmSs = 1.0/(sL-sM);
			  REALQ sLmuL = sL-unL;
			  REALQ rhoSL = rhoL*sLmuL*invSLmSs;
			  REALQ rhouSL[3];
			  for (int i=0; i<3; i++)
				  rhouSL[i] = (rhoL*uL[i]*sLmuL+(pStar-pL)*nVec[i])*invSLmSs;
			  REALQ rhoeSL = (sLmuL*rhoeL-pL*unL+pStar*sM)*invSLmSs;

			  Frho = rhoSL*sM;
			  for (int i=0; i<3; i++)
				  Frhou[i] = rhouSL[i]*sM + pStar*nVec[i];
			  FrhoE = (rhoeSL+pStar)*sM;
			  for (int iScal=0; iScal<nScal; iScal++)
				  FrhoScal[iScal] = Frho * ScalL[iScal];
		  }
	  }
	  else
	  {
		  if (sR >= 0.0)
		  {
			  REALQ invSRmSs = 1.0/(sR-sM);
			  REALQ sRmuR = sR-unR;
			  REALQ rhoSR = rhoR*sRmuR*invSRmSs;
			  REALQ rhouSR[3];
			  for (int i=0; i<3; i++)
				  rhouSR[i] = (rhoR*uR[i]*sRmuR+(pStar-pR)*nVec[i])*invSRmSs;
			  REALQ rhoeSR = (sRmuR*rhoeR-pR*unR+pStar*sM)*invSRmSs;

			  Frho = rhoSR*sM;
			  for (int i=0; i<3; i++)
				  Frhou[i] = rhouSR[i]*sM + pStar*nVec[i];
			  FrhoE = (rhoeSR+pStar)*sM;
			  for (int iScal=0; iScal<nScal; iScal++)
				  FrhoScal[iScal] = Frho * ScalR[iScal];
		  }
		  else
		  {
			  Frho = rhoR*unR;
			  for (int i=0; i<3; i++)
				  Frhou[i] = rhoR*uR[i]*unR + pR*nVec[i];
			  FrhoE = (rhoeR+pR)*unR;
			  for (int iScal=0; iScal<nScal; iScal++)
				  FrhoScal[iScal] = Frho * ScalR[iScal];
		  }
	  }
#endif

	  Frho *= area;
	  for (int i=0; i<3; i++)
		  Frhou[i] *= area;
	  FrhoE *= area;
	  for (int iScal=0; iScal<nScal; iScal++)
		  FrhoScal[iScal] *= area;

	  return 0;
  }








int UgpWithCvCompFlow_AD::calcEulerFluxMatrices_HLLC_AD(double (*A_L)[5], double (*A_R)[5], double (*A_L_Scal)[6], double (*A_R_Scal)[6],
        REALQ rhoL_AD, REALQ *uL_AD, REALQ pL_AD, REALQ TL_AD, REALQ h0_AD, REALQ RL_AD, REALQ gammaL_AD, REALQS *scalL_AD, REALQS kL_AD,
        REALQ rhoR_AD, REALQ *uR_AD, REALQ pR_AD, REALQ TR_AD, REALQ h1_AD, REALQ RR_AD, REALQ gammaR_AD, REALQS *scalR_AD, REALQS kR_AD,
        double area, double *nVec, int nScal, double surfVeloc){

        double rhoL = rhoL_AD.value();
        double rhoR = rhoR_AD.value();
        double pL   = pL_AD.value();
        double pR   = pR_AD.value();
        double TL   = TL_AD.value();
        double TR   = TR_AD.value();
        double h0   = h0_AD.value();
        double h1   = h1_AD.value();
        double RL   = RL_AD.value();
        double RR   = RR_AD.value();
        double gammaL = gammaL_AD.value();
        double gammaR = gammaR_AD.value();

        double *scalL = new double[nScal];
        double *scalR = new double[nScal];

        double uL[3], uR[3];
        for(int i=0;i<3;i++){
	    uL[i] = uL_AD[i].value();
	    uR[i] = uR_AD[i].value();
        }

      for (int iScal=0; iScal<nScal; iScal++) {
	    scalL[iScal]=scalL_AD[iScal].value();
	    scalR[iScal]=scalR_AD[iScal].value();
      }
	    double kL          =kL_AD.value();
	    double kR          =kR_AD.value();

       calcEulerFluxMatrices_HLLC(A_L, A_R, A_L_Scal, A_R_Scal,
                rhoL, uL, pL, TL, h0, RL, gammaL, scalL, kL,
                rhoR, uR, pR, TR, h1, RR, gammaR, scalR, kR,
                area, nVec, nScal, surfVeloc);

       delete [] scalL;
       delete [] scalR;

return 0;
}




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
 */
#ifdef USE_ARTIF_VISC_WITH_MEM_SAVING
void UgpWithCvCompFlow_AD::addViscFlux_AD(REALQ *Frhou, REALQ &FrhoE, double (*A0)[5], double (*A1)[5],
     REALQ rho0,  REALQ *u0,  REALQ (&grad_u0)[3][3],  REALQ h0,  REALQ *grad_h0,  REALQ T0,  REALQ R0,  REALQ gam0,  REALQS kine0,
     REALQ rho1,  REALQ *u1,  REALQ (&grad_u1)[3][3],  REALQ h1,  REALQ *grad_h1,  REALQ T1,  REALQ R1,  REALQ gam1,  REALQS kine1,
     REALQ mul,  REALQS mut,  REALQ lambdaOverCp,  REALQS kine_fa,  REALQ *u_fa,
     double area,  double *nVec,  double smag,  double *sVec, REALQ artifBulkViscosity)
#else
void UgpWithCvCompFlow_AD::addViscFlux_AD(REALQ *Frhou, REALQ &FrhoE, double (*A0)[5], double (*A1)[5],
     REALQ rho0,  REALQ *u0,  REALQ (&grad_u0)[3][3],  REALQ h0,  REALQ *grad_h0,  REALQ T0,  REALQ R0,  REALQ gam0,  REALQS kine0,
     REALQ rho1,  REALQ *u1,  REALQ (&grad_u1)[3][3],  REALQ h1,  REALQ *grad_h1,  REALQ T1,  REALQ R1,  REALQ gam1,  REALQS kine1,
     REALQ mul,  REALQS mut,  REALQ lambdaOverCp,  REALQS kine_fa,  REALQ *u_fa, 
     double area,  double *nVec,  double smag,  double *sVec)
#endif
{
	REALX alpha = vecDotVec3d(nVec, sVec);
	assert((alpha > 0.0) && (alpha < 1.000000001));

	REALQ grad_u_f[3][3];
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			grad_u_f[i][j] = 0.5*(grad_u0[i][j] + grad_u1[i][j]);


	REALX fCorr[3] = { nVec[0] - alpha*sVec[0],
			nVec[1] - alpha*sVec[1],
			nVec[2] - alpha*sVec[2]};

	/*REALX fCorr[3] = { 0.0, 0.0, 0.0};
  alpha = 1.0;*/


	// ========================================================================
	// momentum equation
	// ========================================================================

	REALQ muTotalMomentum = mul + mut;
	REALQ tauij_nj[3];

	for (int i = 0; i < 3; i++)
	{
		tauij_nj[i] = muTotalMomentum * (alpha * (u1[i] - u0[i]) / smag
				+ grad_u_f[i][0] * fCorr[0] + grad_u_f[0][i] * nVec[0]
				+ grad_u_f[i][1] * fCorr[1] + grad_u_f[1][i] * nVec[1]
				+ grad_u_f[i][2] * fCorr[2] + grad_u_f[2][i] * nVec[2]);
	}

	// viscosity times trace of strain rate tensor times 2/3...
#ifdef USE_ARTIF_VISC_WITH_MEM_SAVING
	REALQ tmp = 0.0;
	if(turnOnArtifVisc && artifVisc_bulkViscOnly)
		tmp = (2.0/3.0*muTotalMomentum + artifBulkViscosity) * (grad_u_f[0][0] + grad_u_f[1][1] + grad_u_f[2][2]);
	else
		tmp = (2.0/3.0*muTotalMomentum) * (grad_u_f[0][0] + grad_u_f[1][1] + grad_u_f[2][2]);
#else
	REALQ tmp = 2.0 / 3.0 * muTotalMomentum * (grad_u_f[0][0] + grad_u_f[1][1] + grad_u_f[2][2]);
#endif

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

	REALQ keff = lambdaOverCp + mut / PrTurb;
	REALQ enth = keff * (alpha * (h1 - h0) / smag
			+ 0.5 * ((grad_h0[0] + grad_h1[0]) * fCorr[0]
			                                           +(grad_h0[1] + grad_h1[1]) * fCorr[1]
			                                                                              +(grad_h0[2] + grad_h1[2]) * fCorr[2]));

	// model for triple correlation
	REALQ psi = 0.0;
	if (turbModel > NONE)
		psi = muTotalMomentum * (alpha * (kine1 - kine0) / smag);

	FrhoE = -area * (enth + tauij_nj[0] * u_fa[0] + tauij_nj[1] * u_fa[1] + tauij_nj[2] * u_fa[2] +psi);

	// TBD THIS INVOLVES SOME EXTRA WORK. SHOULD OPTIMIZE.

	addViscFluxJacobians(A0, A1,
			rho0, u0, grad_u0, h0, grad_h0, T0, R0, gam0, kine0,
			rho1, u1, grad_u1, h1, grad_h1, T1, R1, gam1, kine1,
			mul,  mut, lambdaOverCp, kine_fa, u_fa,
			area, nVec, smag, sVec);
}

//--------------------------------------------------------------------------//
//                                                                          //
//                            dF/dU x dU/dQ                                 //
//           multiplications are only performed on non-zero entries         //
//                                                                          //
//--------------------------------------------------------------------------//
#define dF_dQ_10(A,fU,UQ)  {A[1][0] = fU[1][1]*UQ[1][0] + fU[1][2]*UQ[2][0] + fU[1][3]*UQ[3][0];                       } //df(u_im p)/d(rho)
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

 //   implicit part of viscous fluxes

void UgpWithCvCompFlow_AD::addViscFluxJacobians(double (*A0)[5], double (*A1)[5],
    REALQ rho0_AD, REALQ *u0_AD, REALQ (&grad_u0_AD)[3][3], REALQ h0_AD, REALQ *grad_h0_AD, REALQ T0_AD, REALQ R0_AD, REALQ gam0_AD, REALQS kine0_AD,
    REALQ rho1_AD, REALQ *u1_AD, REALQ (&grad_u1_AD)[3][3], REALQ h1_AD, REALQ *grad_h1_AD, REALQ T1_AD, REALQ R1_AD, REALQ gam1_AD, REALQS kine1_AD,
    REALQ mul_AD,  REALQS mut_AD, REALQ lambdaOverCp_AD, REALQS kine_fa_AD, REALQ *u_fa_AD, 
    double area, double *nVec, double smag, double *sVec)
{

  double rho0  = rho0_AD.value() ;
  double rho1  = rho1_AD.value() ;
  double h0    = h0_AD.value() ;
  double h1    = h1_AD.value() ;
  double T0    = T0_AD.value() ;
  double T1    = T1_AD.value() ;
  double R0    = R0_AD.value() ;
  double R1    = R1_AD.value() ;
  double kine0 = kine0_AD.value() ;       //CHANGE FOR SCALAR AD
  double kine1 = kine1_AD.value() ;	  //CHANGE FOR SCALAR AD
  double gam0  = gam0_AD.value() ;
  double gam1  = gam1_AD.value() ;

  double mul     = mul_AD.value() ;
  double mut     = mut_AD.value() ;     // CHANGE FOR SCALAR AD
  double kine_fa = kine_fa_AD.value() ; // CHANGE FOR SCALAR AD
  double lambdaOverCp   = lambdaOverCp_AD.value() ;

  double u0[3],grad_u0[3][3], grad_h0[3], u1[3], grad_u1[3][3], grad_h1[3], u_fa[3] ;


  for (int i=0; i<3; i++){
      u0[i] = u0_AD[i].value();
      u1[i] = u1_AD[i].value();
       for (int j=0; j<3; j++){
         grad_u0[i][j] = grad_u0_AD[i][j].value();
         grad_u1[i][j] = grad_u1_AD[i][j].value();
       }
      grad_h0[i] = grad_h0_AD[i].value();
      grad_h1[i] = grad_h1_AD[i].value();
      u_fa[i]    = u_fa_AD[i].value();
  }


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

  // ========================================================================
  // energy equation
  // ========================================================================

  double keff = lambdaOverCp + mut / PrTurb;

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


void UgpWithCvCompFlow_AD::calcEulerFluxCoupled_HLLC_AD(adouble *EulerFlux,
         adouble rhoL, adouble *uL, adouble pL, adouble TL, adouble h0, adouble RL, adouble gammaL, adouble *ScalL, adouble kL,
         adouble rhoR, adouble *uR, adouble pR, adouble TR, adouble h1, adouble RR, adouble gammaR, adouble *ScalR, adouble kR,
         double area, double *nVec, const int nScal, double *ConvTerm, adouble surfVeloc, const string BC_treatment)
{       
        
  if (BC_treatment == "ALL_TERMS")
  {     
    adouble unL   = vecDotVec3d_AD(uL, nVec);
    adouble uLuL  = vecDotVec3d_AD(uL, uL);
    adouble cL    = sqrt(gammaL * pL / rhoL);
    adouble rhoeL = (h0 + 0.5 * uLuL + kL) * rhoL - pL;
        
    adouble unR   = vecDotVec3d_AD(uR, nVec);
    adouble uRuR  = vecDotVec3d_AD(uR, uR);
    adouble cR    = sqrt(gammaR * pR / rhoR);
    adouble rhoeR = (h1 + 0.5 * uRuR + kR) * rhoR - pR;


    // Roe's aveaging
    adouble Rrho = sqrt(rhoR / rhoL);
    adouble tmp = 1.0 / (1.0 + Rrho);
    adouble velRoe[3];
    for (int i = 0; i < 3; i++)
      velRoe[i] = tmp * (uL[i] + uR[i] * Rrho);
    adouble uRoe  = vecDotVec3d_AD(velRoe, nVec);

    //    double cRoe  = sqrt((gammaL-1.0)*(hRoe- 0.5*vecDotVec3d_AD(velRoe, velRoe)));
    adouble gamPdivRho = tmp * ((gammaL * pL / rhoL + 0.5 * (gammaL - 1.0) * uLuL) + (gammaR * pR / rhoR + 0.5 * (gammaR - 1.0) * uRuR) * Rrho);
    adouble cRoe  = sqrt(gamPdivRho - ((gammaL + gammaR) * 0.5 - 1.0) * 0.5 * vecDotVec3d_AD(velRoe, velRoe));

    // speed of sound at L and R
    adouble sL = min(uRoe - cRoe, unL - cL);
    adouble sR = max(uRoe + cRoe, unR + cR);

    // speed of contact surface
    adouble sM = (pL - pR - rhoL * unL * (sL - unL) + rhoR * unR * (sR - unR)) / (rhoR * (sR - unR) - rhoL * (sL - unL));

  // pressure at right and left (pR=pL) side of contact surface
    adouble pStar = rhoR * (unR - sR) * (unR - sM) + pR;
    
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
        adouble invSLmSs = 1.0 / (sL - sM);
        adouble sLmuL = sL - unL;
        adouble rhoSL = rhoL * sLmuL * invSLmSs;
        adouble rhouSL[3];
        for (int i = 0; i < 3; i++)
          rhouSL[i] = (rhoL * uL[i] * sLmuL + (pStar - pL) * nVec[i]) * invSLmSs;
        adouble rhoeSL = (sLmuL * rhoeL - pL * unL + pStar * sM) * invSLmSs;

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
        adouble invSRmSs = 1.0 / (sR - sM);
        adouble sRmuR = sR - unR;
        adouble rhoSR = rhoR * sRmuR * invSRmSs;
        adouble rhouSR[3];
        for (int i = 0; i < 3; i++)
          rhouSR[i] = (rhoR * uR[i] * sRmuR + (pStar - pR) * nVec[i]) * invSRmSs;
        adouble rhoeSR = (sRmuR * rhoeR - pR * unR + pStar * sM) * invSRmSs;

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



int UgpWithCvCompFlow_AD::calcEulerFluxMatricesCoupled_HLLC_AD(double **A_L, double **A_R,
        REALQ rhoL_AD, REALQ *uL_AD, REALQ pL_AD, REALQ TL_AD, REALQ h0_AD, REALQ RL_AD, REALQ gammaL_AD, REALQS *scalL_AD, REALQ *dpress_dscalL_AD, REALQS kL_AD,
        REALQ rhoR_AD, REALQ *uR_AD, REALQ pR_AD, REALQ TR_AD, REALQ h1_AD, REALQ RR_AD, REALQ gammaR_AD, REALQS *scalR_AD, REALQ *dpress_dscalR_AD, REALQS kR_AD,
        double area, double *nVec, int nScal, double *ConvTerm, double surfVeloc, const string BC_treatment){

        double rhoL = rhoL_AD.value();
        double rhoR = rhoR_AD.value();
        double pL   = pL_AD.value();
        double pR   = pR_AD.value();
        double TL   = TL_AD.value();
        double TR   = TR_AD.value();
        double h0   = h0_AD.value();
        double h1   = h1_AD.value();
        double RL   = RL_AD.value();
        double RR   = RR_AD.value();
        double gammaL = gammaL_AD.value();
        double gammaR = gammaR_AD.value();

        double *scalL = new double[nScal];
        double *scalR = new double[nScal];

        double *dpress_dscalL = new double[nScal];
        double *dpress_dscalR = new double[nScal];

        double uL[3], uR[3];
        for(int i=0;i<3;i++){
            uL[i] = uL_AD[i].value();
            uR[i] = uR_AD[i].value();
        }

      for (int iScal=0; iScal<nScal; iScal++) {
            scalL[iScal]=scalL_AD[iScal].value();
            scalR[iScal]=scalR_AD[iScal].value();
            dpress_dscalL[iScal]=dpress_dscalL_AD[iScal].value();
            dpress_dscalR[iScal]=dpress_dscalR_AD[iScal].value();
      }
            double kL          =kL_AD.value();
            double kR          =kR_AD.value();

       calcEulerFluxMatricesCoupled_HLLC(A_L, A_R,
                rhoL, uL, pL, TL, h0, RL, gammaL, scalL, dpress_dscalL, kL,
                rhoR, uR, pR, TR, h1, RR, gammaR, scalR, dpress_dscalR, kR,
                area, nVec, nScal, ConvTerm, surfVeloc, BC_treatment);

       delete [] scalL;
       delete [] scalR;
       delete [] dpress_dscalL;
       delete [] dpress_dscalR;

return 0;
}


/**
 *  viscous flux and jacobian, coupled NSE and scalars
 */
void UgpWithCvCompFlow_AD::calcViscousFluxCoupled_AD(REALQ *ViscousFlux, double **A0, double **A1,
//void UgpWithCvCompFlow_AD::calcViscousFluxCoupled_AD(REALQ *ViscousFlux, double (*A0)[5], double (*A1)[5],
         REALQ rho0, REALQ *u0, REALQ (&grad_u0)[3][3], REALQ h0, REALQ *grad_h0, REALQ T0, REALQ R0, REALQ gam0, REALQS *Scal0, REALQS (*gradScal0)[3], REALQ *dpress_dscal0, REALQS kine0,
         REALQ rho1, REALQ *u1, REALQ (&grad_u1)[3][3], REALQ h1, REALQ *grad_h1, REALQ T1, REALQ R1, REALQ gam1, REALQS *Scal1, REALQS (*gradScal1)[3], REALQ *dpress_dscal1, REALQS kine1,
         REALQ  mul, REALQ mut, REALQ lambdaOverCp, REALQS kine_fa, REALQ *u_fa, REALQS *diff, double *DiffTerm,
         double area, double *nVec, double smag, double *sVec, double alpha, const int nScal)
{
  REALQ grad_u_f[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      grad_u_f[i][j] = 0.5 * (grad_u0[i][j] + grad_u1[i][j]);


  REALQ fCorr[3] = {nVec[0] - alpha * sVec[0],
                     nVec[1] - alpha * sVec[1],
                     nVec[2] - alpha * sVec[2]};


  // ************************************************************************
  // FLUXES
  // ************************************************************************

  // ========================================================================
  // momentum equation
  // ========================================================================

  REALQ muTotalMomentum = mul + mut;
  REALQ tauij_nj[3];

  // no flux for density
  ViscousFlux[0] = 0.0;

 for (int i = 0; i < 3; i++)
    tauij_nj[i] = muTotalMomentum * (alpha * (u1[i] - u0[i]) / smag
                      + grad_u_f[i][0] * fCorr[0] + grad_u_f[0][i] * nVec[0]
                      + grad_u_f[i][1] * fCorr[1] + grad_u_f[1][i] * nVec[1]
                      + grad_u_f[i][2] * fCorr[2] + grad_u_f[2][i] * nVec[2]);
                      
  // viscosity times trace of strain rate tensor times 2/3... 
  REALQ tmp = 2.0 / 3.0 * muTotalMomentum * (grad_u_f[0][0] + grad_u_f[1][1] + grad_u_f[2][2]);
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

  REALQS keff = lambdaOverCp + mut / PrTurb;
  REALQ  enth = keff * (alpha * (h1 - h0) / smag
                + 0.5 * ((grad_h0[0] + grad_h1[0]) * fCorr[0]
                        +(grad_h0[1] + grad_h1[1]) * fCorr[1]
                        +(grad_h0[2] + grad_h1[2]) * fCorr[2]));

  // model for triple correlation
  REALQ psi = 0.0;
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

  // TBD THIS INVOLVES SOME EXTRA WORK. SHOULD OPTIMIZE.

    addViscFluxJacobiansCoupled(A0, A1,
    rho0, u0, grad_u0, h0, grad_h0, T0, R0, gam0, Scal0, dpress_dscal0, kine0,
    rho1, u1, grad_u1, h1, grad_h1, T1, R1, gam1, Scal1, dpress_dscal1, kine1,
    mul,  mut, lambdaOverCp, kine_fa, u_fa, diff, DiffTerm, 
    area, nVec, smag, sVec, alpha, nScal);

}

void UgpWithCvCompFlow_AD::addViscFluxJacobiansCoupled(double **A0, double **A1,
    REALQ rho0_AD, REALQ *u0_AD, REALQ (&grad_u0_AD)[3][3], REALQ h0_AD, REALQ *grad_h0_AD, REALQ T0_AD, REALQ R0_AD, REALQ gam0_AD, REALQS *Scal0_AD, REALQS *dpress_dscal0_AD, REALQS kine0_AD,
    REALQ rho1_AD, REALQ *u1_AD, REALQ (&grad_u1_AD)[3][3], REALQ h1_AD, REALQ *grad_h1_AD, REALQ T1_AD, REALQ R1_AD, REALQ gam1_AD, REALQS *Scal1_AD, REALQS *dpress_dscal1_AD, REALQS kine1_AD,
    REALQ mul_AD,  REALQS mut_AD, REALQ lambdaOverCp_AD, REALQS kine_fa_AD, REALQ *u_fa_AD, REALQS *diff_AD, double *DiffTerm,
    double area, double *nVec, double smag, double *sVec, double alpha, const int nScal)
{

  double rho0  = rho0_AD.value() ;
  double rho1  = rho1_AD.value() ;
  double h0    = h0_AD.value() ;
  double h1    = h1_AD.value() ;
  double T0    = T0_AD.value() ;
  double T1    = T1_AD.value() ;
  double R0    = R0_AD.value() ;
  double R1    = R1_AD.value() ;
  double kine0 = kine0_AD.value() ;       //CHANGE FOR SCALAR AD
  double kine1 = kine1_AD.value() ;	  //CHANGE FOR SCALAR AD
  double gam0  = gam0_AD.value() ;
  double gam1  = gam1_AD.value() ;

  double mul     = mul_AD.value() ;
  double mut     = mut_AD.value() ;     // CHANGE FOR SCALAR AD
  double kine_fa = kine_fa_AD.value() ; // CHANGE FOR SCALAR AD
  double lambdaOverCp   = lambdaOverCp_AD.value() ;

  double u0[3],grad_u0[3][3], grad_h0[3], u1[3], grad_u1[3][3], grad_h1[3], u_fa[3] ;
  double *Scal0 = new double[nScal];
  double *Scal1 = new double[nScal];
  double *dpress_dscal0 = new double[nScal];
  double *dpress_dscal1 = new double[nScal];
  double *diff = new double[nScal];

  for (int i=0; i<3; i++){
      u0[i] = u0_AD[i].value();
      u1[i] = u1_AD[i].value();
       for (int j=0; j<3; j++){
         grad_u0[i][j] = grad_u0_AD[i][j].value();
         grad_u1[i][j] = grad_u1_AD[i][j].value();
       }
      grad_h0[i] = grad_h0_AD[i].value();
      grad_h1[i] = grad_h1_AD[i].value();
      u_fa[i]    = u_fa_AD[i].value();
  }

   for (int iScal=0; iScal<nScal; iScal++) {
            Scal0[iScal]=Scal0_AD[iScal].value();
            Scal1[iScal]=Scal1_AD[iScal].value();
            dpress_dscal0[iScal]=dpress_dscal0_AD[iScal].value();
            dpress_dscal1[iScal]=dpress_dscal1_AD[iScal].value();
            diff[iScal]=diff_AD[iScal].value();
      }


  double grad_u_f[3][3];
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      grad_u_f[i][j] = 0.5*(grad_u0[i][j] + grad_u1[i][j]);


  double fCorr[3] = { nVec[0] - alpha*sVec[0], 
                      nVec[1] - alpha*sVec[1],
                      nVec[2] - alpha*sVec[2]};

  // ========================================================================
  // momentum equation
  // ========================================================================

  double muTotalMomentum = mul + mut;
  double tauij_nj[3];

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

  // ========================================================================
  // energy equation
  // ========================================================================

  double keff = lambdaOverCp + mut / PrTurb;


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


  





#undef dU_dQ_0
#undef dU_dQ_1
#undef dU_dQ_2
#undef dU_dQ_3
#undef dU_dQ_4
#undef dU_dQ

#undef dF_dQ_10
#undef dF_dQ_20
#undef dF_dQ_30
#undef dF_dQ_40

#undef dF_dQ_11
#undef dF_dQ_21
#undef dF_dQ_31
#undef dF_dQ_41

#undef dF_dQ_12
#undef dF_dQ_22
#undef dF_dQ_32
#undef dF_dQ_42

#undef dF_dQ_13
#undef dF_dQ_23
#undef dF_dQ_33
#undef dF_dQ_43
#undef dF_dQ_44

#undef dF_dQ_0
#undef dF_dQ_1
#undef dF_dQ_2
#undef dF_dQ_3
#undef dF_dQ_4

#undef dFdU_times_dUdQ
