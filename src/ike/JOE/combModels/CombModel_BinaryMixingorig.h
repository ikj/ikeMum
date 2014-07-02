/*
 * CombModel_BinaryMixing.h
 *
 *  Created on: May 6, 2009
 *      Author: Vincent Terrapon
 *      Update: Jan 11, 2010
 */

#ifndef COMBMODEL_BINARYMIXING_H_
#define COMBMODEL_BINARYMIXING_H_

// Analytical deviation from linear mixing rules for single analytical flamelet
#define ROM_DEF R  = zm * RoM_1 + (1.0 - zm) * RoM_2 - dRoM_0 * ( (1.0 - exp(-coeff_1 * zm))                             - pow(zm, coeff_4) * (1.0 - exp(-coeff_1)) );
#define CP_DEF  cp = zm * Cp_1  + (1.0 - zm) * Cp_2  - dcp_0  * ( (1.0 - exp(-coeff_1 * zm)) + coeff_2 * zm * (1.0 - zm) - pow(zm, coeff_3) * (1.0 - exp(-coeff_1)) );

#include "UgpWithCvCompFlow.h"
#include "CombModel_Base.h"


// #############################################################################################
// ------                                                                                 ------
// ----                         RansCombBinaryMixing                                        ----
// ------                                                                                 ------
// #############################################################################################

/*! \brief Binary Mixing model based on mixture fraction. 
 *
 *  The species mass fractions Species::Yk are functions of the mean scalar mixture fraction \p Zmean
 *  only through specific mixing rules. It can represent pure mixing (linear mixing rule) or combustion
 *  (non-linear rule). Viscosity and thermal diffusion use simple temperature dependence rule. 
 *  Specific heat coefficients are not a direct function of the temperature but of the mixture fraction.
 *  Only 1 flamelet can be represented.
 */
//class RansCombBinaryMixing : virtual public UgpWithCvCompFlow {
class RansCombBinaryMixing: public RansCombBase  {
public:
  // member variables

  double  *ZMean;                       ///< Mixture fraction Z.
//  double  *rhoDkZ;                      ///< Diffusion coefficient of mixture fraction Z.
  double  Schmidt_turb_Z;               ///< Turbulent Schmidt number for mixture fraction Z.
  double  Gamma_1;                      ///< Gamma of gas 1.
  double  Gamma_2;                      ///< Gamma of gas 2.
  double  MolMass_1;                    ///< Molecular mass of gas 1 (Z=1).
  double  MolMass_2;                    ///< Molecular mass of gas 2 (Z=0).
  double  RoM_1;                        ///< Gas constant of gas 1 (Z=1).
  double  RoM_2;                        ///< Gas constant of gas 2 (Z=0).
  double  Cp_1;                         ///< Constant specific heat capacity of gas 1 (Z=1).
  double  Cp_2;                         ///< Constant specific heat capacity of gas 2 (Z=0).
  double  rho_1;                        ///< Density of gas 1 (Z=1).
  double  rho_2;                        ///< Density of gas 2 (Z=0).
  int     ZMean_Index;                  ///< Index in scalar vector of ZMean.
  
  double  coeff_1;                      ///< Coefficient to describe analytically non-linear mixing rule
  double  coeff_2;                      ///< Coefficient to describe analytically non-linear mixing rule
  double  coeff_3;                      ///< Coefficient to describe analytically non-linear mixing rule
  double  coeff_4;                      ///< Coefficient to describe analytically non-linear mixing rule
  double  dcp_0;                        ///< Coefficient to describe amplitude of the cp deviation in non-linear mixing rule
  double  dRoM_0;                       ///< Coefficient to describe amplitude of the RoM deviation in non-linear mixing rule

public:
  // constructors

  RansCombBinaryMixing()
  {
    if (mpi_rank == 0)
      cout << "RansCombBinaryMixing()" << endl;

    // Overall constants
    Schmidt_turb_Z = getDoubleParam("SCTURB_ZMean", "1.0");
    if (mpi_rank == 0) 
      cout << "Turbulent Schmidt number for Z Sc_t=" << Schmidt_turb_Z << endl;
  
    // Scalar transport equation for ZMean
    ScalarTranspEq *eq;
    eq = registerScalarTransport("ZMean", CV_DATA);
    eq->diffTerm = getIntParam("SCAL_DIFF", "1");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-4;
    eq->phiMaxiter = 100;
    eq->lowerBound = -10.0;
    eq->upperBound = 10.0;
    eq->reconstruction = getStringParam("SCALAR_RECONSTRUCTION", "CONSERVATIVE");
  }

  /*! \brief Unload myMixing and clean memory. */
  virtual ~RansCombBinaryMixing()
  {
    if (mpi_rank == 0)
      cout << endl << "***** Binary Mixing finalized *****" << endl << endl;
  }

public:
  // member functions 

  /*! \brief Read thermo input file and load mixing boundary conditions. */
  virtual void initialHookScalarRansCombModel()
  {
    MolMass_1 = getDoubleParam("MOLMASS_1", "4.0026");
    MolMass_2 = getDoubleParam("MOLMASS_2", "2.016");
    
    Gamma_1 = getDoubleParam("GAMMA_1", "1.66");
    Gamma_2 = getDoubleParam("GAMMA_2", "1.41");
    
    RoM_1 = 8314.4 / MolMass_1;
    RoM_2 = 8314.4 / MolMass_2;
    
    Cp_1 = Gamma_1 * RoM_1 / (Gamma_1 - 1.0);
    Cp_2 = Gamma_2 * RoM_2 / (Gamma_2 - 1.0);   
    
    coeff_1 = getDoubleParam("NON_LIN_MIX_COEFF_1", "50.0");
    coeff_2 = getDoubleParam("NON_LIN_MIX_COEFF_2", "10.0");
    coeff_3 = getDoubleParam("NON_LIN_MIX_COEFF_3", "10.0");
    coeff_4 = getDoubleParam("NON_LIN_MIX_COEFF_4", "10.0");
    
    // Default values correspond to linear mixing rule
    dcp_0  = getDoubleParam("NON_LIN_MIX_D_CP", "0.0");
    dRoM_0 = getDoubleParam("NON_LIN_MIX_D_R", "0.0");
  
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");           ZMean = eq->phi;           ZMean_Index = getScalarTransportIndex("ZMean");

    if (mpi_rank == 0)
      cout << endl << "***** Binary Mixing initialized *****" << endl << endl;
  }

  /*! \brief Compute for each cell the mixture properties as function of \p Zmean.
   *
   *  Compute composition Species::Yk, temperature, viscosity Mixture::mul, heat coefficient Mixture::lambda,
   *  enthalpy Mixture::H, heat capacity ratio Mixture::gama, pressure and gas constant for the mixture.
   */
  
  virtual void calcStateVariables()
  {
    double zm, R, E, cp;
    double kinecv = 0.0;

    for (int icv = 0; icv < ncv_gg; icv++)
    {
      for (int i=0; i<3; i++)
        vel[icv][i] = rhou[icv][i]/rho[icv];

      if (kine != NULL)
        kinecv = kine[icv];
      
      // Compute mixture gas constant and specific heat coefficient
      zm = ZMean[icv];
      ROM_DEF
      RoM[icv] = R;
      CP_DEF

      E = rhoE[icv] / rho[icv] - 0.5 * vecDotVec3d(vel[icv], vel[icv]) - kinecv;

      // Compute mixture temperature from energy and species mass fractions
      temp[icv] = E / (cp - RoM[icv]);

      // Compute mixture enthalpy, pressure, gamma and speed of sound
      enthalpy[icv] = cp * temp[icv];
      press[icv] = rho[icv] * RoM[icv] * temp[icv];
      gamma[icv] = cp / (cp - RoM[icv]);
      sos[icv] = sqrt(gamma[icv]*press[icv]/rho[icv]);
   }

    
    // viscosity and thermal diffusivity computed in UgpWithCvCompFlow::calcMaterialProperties() at cell faces
  }
  
  // \brief Compute for a given temperature, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_T(double &p, double &h, double &R, double &gam, double &rho, double &T, double *Scal, int nScal)
  {
    double cp;
    double zm = Scal[ZMean_Index];
    ROM_DEF
    CP_DEF
    gam = cp / (cp - R);
    p = rho * R * T;
    h = cp * T;
  }

  // \brief Compute for a given pressure, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_p(double &T, double &h, double &R, double &gam, double &rho, double &p, double *Scal, int nScal)
  {
    double cp;
    double zm = Scal[ZMean_Index];
    ROM_DEF
    CP_DEF
    gam = cp / (cp - R);
    T = p / (rho * R);
    h = cp * T;
  }
  
 
  // \brief Compute for a given temperature the properties of the mixture at a face.
  virtual void ComputeBCProperties_T(FaZone *zone)
  {
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");         double *zM_fa = eq->phi;
    
      for (int index = 0; index < zone->faVec.size(); ++index) {
      int ifa = zone->faVec[index];
      int icv1 = cvofa[ifa][1];
      ComputeProperties_T(enthalpy[icv1],RoM[icv1],gamma[icv1],temp[icv1],zM_fa[icv1],0.0,0.0);
      }
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
      int icv1 = cvofa[ifa][1];
      ComputeProperties_vis(mul_fa[ifa],lamOcp_fa[ifa],temp[icv1]);
      }
  }
 
  
  // \brief Compute for a given enthalpy the properties of the mixture at a face.
  virtual void ComputeBCProperties_H(FaZone *zone)
  {
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");         double *zM_fa = eq->phi;
    
      for (int index = 0; index < zone->faVec.size(); ++index) {
      int ifa = zone->faVec[index];
      int icv1 = cvofa[ifa][1];
      ComputeProperties_H(temp[icv1],RoM[icv1],gamma[icv1],enthalpy[icv1],zM_fa[icv1],0.0,0.0);
      }
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
      int icv1 = cvofa[ifa][1];
      ComputeProperties_vis(mul_fa[ifa],lamOcp_fa[ifa],temp[icv1]);
      }
  }
  
  /*! \brief Compute for a given temperature the properties of the mixture.
   *
   *  The enthalpy (chemical + sensible) is a function of temperature and mixture composition. 
   *  \param[out] H Enthalpy (chemical + sensible) in [J/kg].
   *  \param[out] RoverM Gas constant R/M of the mixture in [J/(kg K)].
   *  \param[out] gam Heat capacity ratio in [-].
   *  \param[out] mul Laminar viscosity in [kg/(m s)].
   *  \param[out] lam_o_cp Heat conductivity divided by heat capacity in [kg/(m s)].
   *  \param[in]  T Temperature in [K].
   *  \param[in]  zm Mean mixture fraction.
   *  \param[in]  zv Variance of mixture fraction.
   *  \param[in]  chi Scalar dissipation rate.
   */
 
  void ComputeProperties_T(double &h, double &R, double &gam, double T, double zm, double zv, double chi)
  {
    double cp;
    ROM_DEF
    CP_DEF
    h = cp * T;
    gam = cp / (cp - R);
  }

  void ComputeProperties_vis(double &mu, double &lamOcp, double T)
  {
    if (mu_ref > 0.0)
    {
      mu = calcMuLam(T);
      lamOcp = mu / Pr;
    }
  }

  /*! \brief Compute for a given enthalpy the properties of the mixture.
   *
   *  The enthalpy (chemical + sensible) is a function of temperature and mixture composition. 
   *  \param[out] T Temperature in [K].
   *  \param[out] RoverM Gas constant R/M of the mixture in [J/(kg K)].
   *  \param[out] gam Heat capacity ratio in [-].
   *  \param[out] mul Laminar viscosity in [kg/(m s)].
   *  \param[out] lam_o_cp Heat conductivity divided by heat capacity in [kg/(m s)].
   *  \param[in]  H Enthalpy (chemical + sensible) in [J/kg].
   *  \param[in]  zm Mean mixture fraction.
   *  \param[in]  zv Variance of mixture fraction.
   *  \param[in]  chi Scalar dissipation rate.
   */
 
  void ComputeProperties_H(double &T, double &R, double &gam, double h, double zm, double zv, double chi)
  {
    double cp;
    ROM_DEF
    CP_DEF
    T = h / cp;
    gam = cp / (cp - R);
  }
 
 
  virtual void diffusivityHookScalarRansComb(const string &name)
  {
    ScalarTranspEq *eq;

    if (name == "ZMean")
    {
      eq = getScalarTransportData(name);
 
      // internal faces
      for (int ifa = nfa_b; ifa < nfa; ifa++)
        eq->diff[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;    

      // boundary faces
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
            eq->diff[ifa] = lamOcp_fa[ifa];
        else
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
           eq->diff[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;
      }
    }
  }
  


};


#endif  /* COMBMODEL_BINARYMIXING_H_ */
