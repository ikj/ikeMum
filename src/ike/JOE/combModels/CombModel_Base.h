/*! \brief Combustion base class.
 *
 * \author Vincent Terrapon 
 * \date May 2009
 * \version 2.0
 */

#ifndef RANSCOMBMODEL_BASE_H
#define RANSCOMBMODEL_BASE_H

#include "UgpWithCvCompFlow.h"
#include "Combustion.h"

// #############################################################################################
// ------                                                                                 ------
// ----                                  RansCombBase                                       ----
// ------                                                                                 ------
// #############################################################################################

/*! \brief Base class containing the combustion specific HLLC fluxes and jacobians. 
 *
 *  Base combustion class with combustion specific HLLC flux and corresponding Jacobians calculation.
 */
class RansCombBase: virtual public UgpWithCvCompFlow {
public:
  // member variables
  
public:
  // constructors

  RansCombBase()
  {
    if (mpi_rank == 0)
      cout << "RansCombBase()" << endl;
  }

  virtual ~RansCombBase()
  {
  }

public:
  // member functions 

  virtual void initialHookScalarRansCombModel() {/*empty*/} 
  virtual void sourceHookRansComb(double *rhsRho, double(*rhsRhou)[3], double *rhsRhoE, double(*A)[5][5]) {/*empty*/}
  virtual void sourceHookScalarRansComb(double *rhs, double *A, const string &name,int flagImplicit) {/*empty*/}
  virtual void sourceHookRansCombCoupled(double **rhs, double ***A, int nScal, int flagImplicit) {/*empty*/}
  virtual void boundaryHookScalarRansComb(double *phi_ph, FaZone *zone, const string &name) {/*empty*/}  
  virtual void UserDefinedScalarClipping(const string &name) {/*empty*/}

  
  

  /**
   *
   *  explicit euler flux HLLC
   *
   */
  virtual int calcEulerFlux_HLLC(double &Frho, double *Frhou, double &FrhoE, double *FrhoScal,
      const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *ScalL, const double kineL,
      const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *ScalR, const double kineR,
      const double area, const double *nVec, const int nScal, const double surfVeloc)
  {
    double unL   = vecDotVec3d(uL, nVec);
    double uLuL  = vecDotVec3d(uL, uL);
    double cL    = sqrt(gammaL*pL/rhoL);
    double rhoeL = (h0 + 0.5*uLuL + kineL)*rhoL-pL;

    double unR   = vecDotVec3d(uR, nVec);
    double uRuR  = vecDotVec3d(uR, uR);
    double cR    = sqrt(gammaR*pR/rhoR);
    double rhoeR = (h1 + 0.5*uRuR + kineR)*rhoR-pR;


    // Roe's aveaging
    double Rrho = sqrt(rhoR/rhoL);
    double tmp = 1.0/(1.0+Rrho);
    double velRoe[3];
    for (int i=0; i<3; i++)
      velRoe[i] = tmp*(uL[i] + uR[i]*Rrho);
    double uRoe  = vecDotVec3d(velRoe, nVec);

    //    double cRoe  = sqrt((gammaL-1.0)*(hRoe- 0.5*vecDotVec3d(velRoe, velRoe)));
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
        Frho = rhoL*unL;
        for (int i=0; i<3; i++)
          Frhou[i] = rhoL*uL[i]*unL + pL*nVec[i];
        FrhoE = (rhoeL+pL)*unL;
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
        double rhoeSL = (sLmuL*rhoeL-pL*unL+pStar*sM)*invSLmSs;

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
        double invSRmSs = 1.0/(sR-sM);
        double sRmuR = sR-unR;
        double rhoSR = rhoR*sRmuR*invSRmSs;
        double rhouSR[3];
        for (int i=0; i<3; i++)
          rhouSR[i] = (rhoR*uR[i]*sRmuR+(pStar-pR)*nVec[i])*invSRmSs;
        double rhoeSR = (sRmuR*rhoeR-pR*unR+pStar*sM)*invSRmSs;

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
   * HLLC Euler flux and Jacobian routines used by all combustion models
   *
   */
  virtual int calcEulerFluxMatrices_HLLC(double (*A_L)[5], double (*A_R)[5], double (*A_L_Scal)[6], double (*A_R_Scal)[6],
          const double rhoL, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *scalL, const double kineL,
          const double rhoR, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *scalR, const double kineR,
          const double area, const double *nVec, const int nScal, const double surfVeloc)
  {
    double unL   = vecDotVec3d(uL, nVec);
    double uLuL  = vecDotVec3d(uL, uL);
    double cL    = sqrt(gammaL*pL/rhoL);
    double hL    = h0 + kineL;
    double rhoeL = (hL + 0.5*uLuL)*rhoL-pL;

    double unR   = vecDotVec3d(uR, nVec);
    double uRuR  = vecDotVec3d(uR, uR);
    double cR    = sqrt(gammaR*pR/rhoR);
    double hR    = h1 + kineR;
    double rhoeR = (hR + 0.5*uRuR)*rhoR-pR;


    // Roe's aveaging
    double Rrho = sqrt(rhoR/rhoL);
    double tmp = 1.0/(1.0+Rrho);
    double velRoe[3];
    for (int i=0; i<3; i++)
      velRoe[i] = tmp*(uL[i] + uR[i]*Rrho);
    double uRoe  = vecDotVec3d(velRoe, nVec);

    double gamPdivRho = tmp*((gammaL*pL/rhoL+0.5*(gammaL-1.0)*uLuL) + (gammaR*pR/rhoR+0.5*(gammaR-1.0)*uRuR)*Rrho);
    double cRoe  = sqrt(gamPdivRho - ((gammaL+gammaR)*0.5-1.0)*0.5*vecDotVec3d(velRoe, velRoe));

    // speed of sound at L and R
    double sL = min(uRoe-cRoe, unL-cL);
    double sR = max(uRoe+cRoe, unR+cR);

    // speed of contact surface
    double invrhotld = 1.0 / (rhoR*(sR-unR)-rhoL*(sL-unL));
    double sM = (pL-pR-rhoL*unL*(sL-unL)+rhoR*unR*(sR-unR)) * invrhotld;

    // pressure at right and left (pR=pL) side of contact surface
    double pStar = rhoR*(unR-sR)*(unR-sM)+pR;

    // normal vector scaled by area
    double nVecArea[3];
    for (int i=0; i<3; i++)        
      nVecArea[i] = nVec[i]*area;
    
    if (sM >= 0.0)
    {
      if (sL > 0.0)
      {
        if (A_L != NULL)
        {
          calcJacobianA(A_L, uL, pL, rhoL, hL, nVecArea, gammaL, 0.0);

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
        double OmegaL = 1.0/(sL-sM);
        double sLmuL = sL-unL;
        double rhoSL = rhoL*sLmuL*OmegaL;
        double rhouSL[3];
        for (int i=0; i<3; i++)
          rhouSL[i] = (rhoL*uL[i]*sLmuL+(pStar-pL)*nVec[i])*OmegaL;
        double rhoeSL = (sLmuL*rhoeL-pL*unL+pStar*sM)*OmegaL;

        double dSMdUL[5], dSMdUR[5];
        
        if (A_L != NULL)
          {
            calcdSMdU(dSMdUL, uL, pL, rhoL, hL, nVec, gammaL, sM, sL, invrhotld, 1.0);
            calcJacobianA(A_L, uL, pL, rhoL, hL, nVecArea, gammaL, 0.0);
            double rhoSmqR = rhoR * (sR - unR);
            calcSubSonicJacobeanHLLC_COMB(A_L, dSMdUL, rhoSL, rhouSL, rhoeSL, pStar, OmegaL, sM, rhoSmqR, nVec, area);
            for (int i=0; i<5; i++)
              A_L[i][i] += area*OmegaL*sM*sL;

            for (int iScal = 0; iScal < nScal; iScal++)
            {
              // F(rhophi_L^star) = F(rho_L^star) * phi_L = F(rho_L^star) * rhophi_L / rho_L = OmegaL * (SL - q_L) * SM * area * rho_L * rhophi_L / rho_L
              
              // dF(rhophi_L^star) / drho_L    = dF(rho_L^star) / drho_L  * phi_L - F(rho_L^star) * phi_L / rho_L = A_L[0,0] * phi_L - F(rho_L^star) * phi_L / rho_L
              // dF(rhophi_L^star) / drhou_L   = dF(rho_L^star) / drhou_L * phi_L                                 = A_L[0,1] * phi_L
              // dF(rhophi_L^star) / drhoe_L   = dF(rho_L^star) / drhoe_L * phi_L                                 = A_L[0,4] * phi_L
              // dF(rhophi_L^star) / drhophi_L =  F(rho_L^star) / rho_L                                           = F(rho_L^star) / rho_L
              
              double FrhostarOverrho = OmegaL * sLmuL * sM * area;
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
            calcdSMdU(dSMdUR, uR, pR, rhoR, hR, nVec, gammaR, sM, sR, invrhotld, -1.0);
            for (int i=0; i<5; i++)
              for (int j=0; j<5; j++)
                A_R[i][j] = 0.0;
            double rhoSmqL = rhoL * (sL - unL);
            calcSubSonicJacobeanHLLC_COMB(A_R, dSMdUR, rhoSL, rhouSL, rhoeSL, pStar, OmegaL, sM, rhoSmqL, nVec, area);

            for (int iScal = 0; iScal < nScal; iScal++)
            {
              // F(rhophi_L^star) = F(rho_L^star) * phi_L = F(rho_L^star) * rhophi_L / rho_L = OmegaL * (SL - q_L) * SM * area * rho_L * rhophi_L / rho_L
              
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
        double OmegaR = 1.0/(sR-sM);
        double sRmuR = sR-unR;
        double rhoSR = rhoR*sRmuR*OmegaR;
        double rhouSR[3];
        for (int i=0; i<3; i++)
          rhouSR[i] = (rhoR*uR[i]*sRmuR+(pStar-pR)*nVec[i])*OmegaR;
        double rhoeSR = (sRmuR*rhoeR-pR*unR+pStar*sM)*OmegaR;

        double dSMdUL[5], dSMdUR[5];
        
        if (A_R != NULL)
          {
            calcdSMdU(dSMdUR, uR, pR, rhoR, hR, nVec, gammaR, sM, sR, invrhotld, -1.0);
            calcJacobianA(A_R, uR, pR, rhoR, hR, nVecArea, gammaR, 0.0);
            double rhoSmqL = rhoL * (sL - unL);
            calcSubSonicJacobeanHLLC_COMB(A_R, dSMdUR, rhoSR, rhouSR, rhoeSR, pStar, OmegaR, sM, rhoSmqL, nVec, area);
            for (int i=0; i<5; i++)
              A_R[i][i] += area*OmegaR*sM*sR;

            for (int iScal = 0; iScal < nScal; iScal++)
            {
              // F(rhophi_R^star) = F(rho_R^star) * phi_R = F(rho_R^star) * rhophi_R / rho_R = OmegaR * (SR - q_R) * SM * area * rho_R * rhophi_R / rho_R
              
              // dF(rhophi_R^star) / drho_R    = dF(rho_R^star) / drho_R  * phi_R - F(rho_R^star) * phi_R / rho_R = A_R[0,0] * phi_R - F(rho_R^star) * phi_R / rho_R
              // dF(rhophi_R^star) / drhou_R   = dF(rho_R^star) / drhou_R * phi_R                                 = A_R[0,1] * phi_R
              // dF(rhophi_R^star) / drhoe_R   = dF(rho_R^star) / drhoe_R * phi_R                                 = A_R[0,4] * phi_R
              // dF(rhophi_R^star) / drhophi_R =  F(rho_R^star) / rho_R                                           = F(rho_R^star) / rho_R
              
              double FrhostarOverrho = OmegaR * sRmuR * sM * area;
              A_R_Scal[iScal][0] = (A_R[0][0] - FrhostarOverrho) * scalR[iScal];     // d(rho phi un)/d(rho)
              A_R_Scal[iScal][1] =  A_R[0][1] * scalR[iScal];                        // d(rho phi un)/d(rho u)
              A_R_Scal[iScal][2] =  A_R[0][2] * scalR[iScal];                        // d(rho phi un)/d(rho v)
              A_R_Scal[iScal][3] =  A_R[0][3] * scalR[iScal];                        // d(rho phi un)/d(rho w)
              A_R_Scal[iScal][4] =  A_R[0][4] * scalR[iScal];                        // d(rho phi un)/d(rho e)
              A_R_Scal[iScal][5] =  FrhostarOverrho;                                 // d(rho phi un)/d(rho phi)
            }
        }

        if (A_L != NULL)
          {
            calcdSMdU(dSMdUL, uL, pL, rhoL, hL, nVec, gammaL, sM, sL, invrhotld, 1.0);
            for (int i=0; i<5; i++)
              for (int j=0; j<5; j++)
                A_L[i][j] = 0.0;
            double rhoSmqR = rhoR * (sR - unR);
            calcSubSonicJacobeanHLLC_COMB(A_L, dSMdUL, rhoSR, rhouSR, rhoeSR, pStar, OmegaR, sM, rhoSmqR, nVec, area);

            for (int iScal = 0; iScal < nScal; iScal++)
            {
              // F(rhophi_R^star) = F(rho_R^star) * phi_R = F(rho_R^star) * rhophi_R / rho_R = OmegaR * (SR - q_R) * SM * area * rho_R * rhophi_R / rho_R
              
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
         
      }
      else
      {
        if (A_R != NULL)
        {
          calcJacobianA(A_R, uR, pR, rhoR, hR, nVecArea, gammaR, 0.0);

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
   *  implizit euler flux HLLC2
   *
   */
  void calcSubSonicJacobeanHLLC_COMB(double (*A)[5], double *dSMdU, double rhoS, double *rhouS, double rhoeS, 
                                                   double pS, double Omega, double sM, double rhoSmq, const double *nV, double area)
  {
    double OmSM = Omega * sM;
    double OmSMp1 = OmSM + 1.0;
    double temp;
    
    for (int i=0; i<5; i++)
      {
        A[0][i] = -OmSM * A[0][i] + area * OmSMp1 * dSMdU[i] * rhoS;
        A[1][i] = -OmSM * A[1][i] + area * OmSMp1 * dSMdU[i] * (rhouS[0] + rhoSmq * nV[0]);
        A[2][i] = -OmSM * A[2][i] + area * OmSMp1 * dSMdU[i] * (rhouS[1] + rhoSmq * nV[1]);
        A[3][i] = -OmSM * A[3][i] + area * OmSMp1 * dSMdU[i] * (rhouS[2] + rhoSmq * nV[2]);
        A[4][i] = -OmSM * A[4][i] + area * OmSMp1 * dSMdU[i] * (rhoeS + pS + rhoSmq * sM);
     }
  }


  /**
   *
   * jacobean 
   *
   */  
  void calcJacobianA(double (*A)[5], const double *vel, double pp, double rrho, double hh, const double *nV, double gamma, double surfVeloc) // nV is not normalized
  {
    double kapm1 = (gamma - 1.0);
    double nVel[3];
    nVel[0] = vel[0]*nV[0];
    nVel[1] = vel[1]*nV[1];
    nVel[2] = vel[2]*nV[2];
    double U_k = nVel[0]+nVel[1]+nVel[2];
    double vSquHlf = 0.5*vecDotVec3d(vel, vel);
  //  double c = sqrt(gamma*pp/rrho);
  //  double inv_kap_m1 = 1.0/kapm1;

    // A Matrix
  //  A[0][0] =-surfVeloc;
    A[0][0] = 0.0;
    A[0][1] = nV[0];
    A[0][2] = nV[1];
    A[0][3] = nV[2];
    A[0][4] = 0.0;

    A[1][0] = -vel[0]*U_k+nV[0]*(kapm1*vSquHlf+gamma*pp/rrho-kapm1*hh);   //Changed here for new energy definition
  //  A[1][1] = (2.0-gamma)*nVel[0]+U_k-surfVeloc;
    A[1][1] = (2.0-gamma)*nVel[0]+U_k;
    A[1][2] = vel[0]*nV[1]-kapm1*vel[1]*nV[0];
    A[1][3] = vel[0]*nV[2]-kapm1*vel[2]*nV[0];
    A[1][4] = kapm1*nV[0];

    A[2][0] = -vel[1]*U_k+nV[1]*(kapm1*vSquHlf+gamma*pp/rrho-kapm1*hh);   //Changed here for new energy definition
    A[2][1] = vel[1]*nV[0]-kapm1*vel[0]*nV[1];
  //  A[2][2] = (2.-gamma)*nVel[1]+U_k-surfVeloc;
    A[2][2] = (2.-gamma)*nVel[1]+U_k;
    A[2][3] = vel[1]*nV[2]-kapm1*vel[2]*nV[1];
    A[2][4] = kapm1*nV[1];

    A[3][0] = -vel[2]*U_k+nV[2]*(kapm1*vSquHlf+gamma*pp/rrho-kapm1*hh);   //Changed here for new energy definition
    A[3][1] = vel[2]*nV[0]-kapm1*vel[0]*nV[2];
    A[3][2] = vel[2]*nV[1]-kapm1*vel[1]*nV[2];
  //  A[3][3] = (2.-gamma)*nVel[2]+U_k-surfVeloc;
    A[3][3] = (2.-gamma)*nVel[2]+U_k;
    A[3][4] = kapm1*nV[2];

    A[4][0] = U_k*((gamma-2.)*vSquHlf-gamma*(hh-pp/rrho));                //Changed here for new energy definition
    A[4][1] = (hh+vSquHlf)*nV[0]-kapm1*vel[0]*U_k;                        //Changed here for new energy definition
    A[4][2] = (hh+vSquHlf)*nV[1]-kapm1*vel[1]*U_k;                        //Changed here for new energy definition
    A[4][3] = (hh+vSquHlf)*nV[2]-kapm1*vel[2]*U_k;                        //Changed here for new energy definition
  //  A[4][4] = gamma*U_k-surfVeloc;
    A[4][4] = gamma*U_k;
  }


  /**
   *
   * jacobean 
   *
   */  
  void calcdSMdU(double *dSMdU, const double *vel, double pp, double rrho, double hh, 
                                    const double *nV, double gamma, double SM, double S, double invertrho, double Factor) // nV is not normalized
  {
    double kapm1 = (gamma - 1.0);
    double nVel[3];
    nVel[0] = vel[0]*nV[0];
    nVel[1] = vel[1]*nV[1];
    nVel[2] = vel[2]*nV[2];
    double U_k = nVel[0]+nVel[1]+nVel[2];
    double vSquHlf = 0.5*vecDotVec3d(vel, vel);
    
    dSMdU[0] = Factor*invertrho * (-U_k*U_k + kapm1*vSquHlf + gamma*pp/rrho - kapm1*hh + SM*S);
    dSMdU[1] = Factor*invertrho * ((2*U_k - S - SM)*nV[0] - kapm1*vel[0]);
    dSMdU[2] = Factor*invertrho * ((2*U_k - S - SM)*nV[1] - kapm1*vel[1]);
    dSMdU[3] = Factor*invertrho * ((2*U_k - S - SM)*nV[2] - kapm1*vel[2]);
    dSMdU[4] = Factor*invertrho * kapm1;
    
    return;
  }
};

#endif
