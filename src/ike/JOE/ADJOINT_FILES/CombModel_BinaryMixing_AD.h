#ifndef COMBMODEL_BINARYMIXING_AD_H
#define COMBMODEL_BINARYMIXING_AD_H

#include "../combModels/CombModel_Base.h"
#include "../combModels/CombModel_BinaryMixing.h"


class RansCombBinaryMixing_AD : virtual public UgpWithCvCompFlow_AD, public RansCombBinaryMixing {

public:
  // constructors

  RansCombBinaryMixing_AD()
  {
    if (mpi_rank == 0)
      cout << "RansCombBinaryMixing_AD()" << endl;

    psi_ZMean = NULL; 	registerScalar(psi_ZMean,"PSI_ZMEAN", CV_DATA); ZMean_Index = getScalarTransportIndex("ZMean");
  }

  /*! \brief Unload myMixing and clean memory. */
  virtual ~RansCombBinaryMixing_AD() {}

public:   // member vars

  adouble *ZMean_diff, *zM_fa, *ZMean;

public:
  // member functions 

   virtual void initialize_comb_adjoint()
  {
    if (!checkDataFlag(psi_ZMean)) {
        if(mpi_rank==0) cout<<"INITIALIZING ADJOINT ZMEAN TO ZERO"<<endl;
      for (int icv=0; icv<ncv_gg; icv++)
        psi_ZMean[icv] = 0.;
    }

        updateCvDataG1G2(psi_ZMean,  REPLACE_DATA);

       int nScal = scalarTranspEqVector.size();

         for(int i=0; i<nScal ; i++){
                if(strcmp(scalarTranspEqVector_AD[i].name,"ZMean")==0){
                  for (int icv = 0; icv < ncv_gg; icv++) {
                        scalarTranspEqVector_psi[i].phi[icv]=psi_ZMean[icv];
                  }
                }
         }
  }

   virtual void copy_comb_adjoint()
  {
    int nScal = scalarTranspEqVector.size();

         for(int i=0; i<nScal ; i++){
                if(strcmp(scalarTranspEqVector_AD[i].name,"ZMean")==0){
                  for (int icv = 0; icv < ncv_gg; icv++)
                        psi_ZMean[icv]=scalarTranspEqVector_psi[i].phi[icv];
                }
         }
  }


  virtual void initialHookScalarRansCombModel_AD(int loadtable)
  {
     // connect pointers 

    int nScal = scalarTranspEqVector.size();

         for(int i=0; i<nScal ; i++){
                if(strcmp(scalarTranspEqVector_AD[i].name,"ZMean")==0){
                        ZMean      = scalarTranspEqVector_AD[i].phi ;
                        ZMean_diff = scalarTranspEqVector_AD[i].diff ;
                if(mpi_rank==0) cout<<"Connected Scalar Pointer "<<scalarTranspEqVector_AD[i].name<<endl;
                }
         }

  }

  /*! \brief Compute for each cell the mixture properties as function of \p Zmean.
   *
   *  Compute composition Species::Yk, temperature, viscosity Mixture::mul, heat coefficient Mixture::lambda,
   *  enthalpy Mixture::H, heat capacity ratio Mixture::gama, pressure and gas constant for the mixture.
   */
  virtual void calcStateVariables_AD(adouble *rho, adouble (*rhou)[3], adouble *rhoE)
  {
    adouble zm, R, E, cp;
    adouble kinecv = 0.0;

    for (int icv = 0; icv < ncv_gg; icv++)
    {
      for (int i=0; i<3; i++)
        vel[icv][i] = rhou[icv][i]/rho[icv];

      if (kine != NULL)
        kinecv = kine[icv];
      
      // Compute mixture gas constant and specific heat coefficient
      zm = ZMean[icv];
      R  = zm * RoM_1 + (1.0 - zm) * RoM_2 - dRoM_0 * ( (1.0 - exp(-coeff_1 * zm))                             - pow(zm, coeff_4) * (1.0 - exp(-coeff_1)) );
      RoM[icv] = R;
      cp = zm * Cp_1  + (1.0 - zm) * Cp_2  - dcp_0  * ( (1.0 - exp(-coeff_1 * zm)) + coeff_2 * zm * (1.0 - zm) - pow(zm, coeff_3) * (1.0 - exp(-coeff_1)) );

      E = rhoE[icv] / rho[icv] - 0.5 * vecDotVec3d_AD(vel[icv], vel[icv]) - kinecv;

      // Compute mixture temperature from energy and species mass fractions
      temp[icv] = E / (cp - RoM[icv]);

      // Compute mixture enthalpy, pressure, gamma and speed of sound
      enthalpy[icv] = cp * temp[icv];
      press[icv] = rho[icv] * RoM[icv] * temp[icv];
      gamma[icv] = cp / (cp - RoM[icv]);
      sos[icv] = sqrt(gamma[icv]*press[icv]/rho[icv]);
   }

  }
  
  // \brief Compute for a given temperature, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_T_AD(adouble &p, adouble &h, adouble &R, adouble &gam, adouble &rho, adouble &T, adouble *Scal, int nScal)
  {
    adouble cp;
    adouble zm = Scal[ZMean_Index];
    R  = zm * RoM_1 + (1.0 - zm) * RoM_2 - dRoM_0 * ( (1.0 - exp(-coeff_1 * zm))                             - pow(zm, coeff_4) * (1.0 - exp(-coeff_1)) );
    cp = zm * Cp_1  + (1.0 - zm) * Cp_2  - dcp_0  * ( (1.0 - exp(-coeff_1 * zm)) + coeff_2 * zm * (1.0 - zm) - pow(zm, coeff_3) * (1.0 - exp(-coeff_1)) );
    gam = cp / (cp - R);
    p = rho * R * T;
    h = cp * T;
  }
 
  // \brief Compute for a given pressure, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_p_AD(adouble &T, adouble &h, adouble &R, adouble &gam, adouble &rho, adouble &p, adouble *Scal, int nScal)
  {
    adouble cp;
    adouble zm = Scal[ZMean_Index];
    R  = zm * RoM_1 + (1.0 - zm) * RoM_2 - dRoM_0 * ( (1.0 - exp(-coeff_1 * zm))                             - pow(zm, coeff_4) * (1.0 - exp(-coeff_1)) );
    cp = zm * Cp_1  + (1.0 - zm) * Cp_2  - dcp_0  * ( (1.0 - exp(-coeff_1 * zm)) + coeff_2 * zm * (1.0 - zm) - pow(zm, coeff_3) * (1.0 - exp(-coeff_1)) );
    gam = cp / (cp - R);
    T = p / (rho * R);
    h = cp * T;
  }
  
  // \brief Compute for a given temperature the properties of the mixture at a face.
  virtual void ComputeBCProperties_T_AD(FaZone *zone)
  {
    for (int index = 0; index < zone->faVec.size(); ++index) {
      int ifa = zone->faVec[index];
      int icv1=cvofa[ifa][1];
      ComputeProperties_T_AD(enthalpy[icv1],RoM[icv1],gamma[icv1],temp[icv1],ZMean[icv1],0.0,0.0);
    }
    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
      int icv1 = cvofa[ifa][1];
      ComputeProperties_vis_AD(mul_fa[ifa],lamOcp_fa[ifa],temp[icv1]);
    }

  }
  
  // \brief Compute for a given enthalpy the properties of the mixture at a face.
  virtual void ComputeBCProperties_H_AD(FaZone *zone)
  {
     for (int index = 0; index < zone->faVec.size(); ++index) {
      int ifa = zone->faVec[index];
      int icv1=cvofa[ifa][1];
      ComputeProperties_H_AD(temp[icv1],RoM[icv1],gamma[icv1],enthalpy[icv1],ZMean[icv1],0.0,0.0);
     }
     for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
      int icv1 = cvofa[ifa][1];
      ComputeProperties_vis_AD(mul_fa[ifa],lamOcp_fa[ifa],temp[icv1]);
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
  void ComputeProperties_T_AD(adouble &h, adouble &R, adouble &gam, adouble T, adouble zm, adouble zv, adouble chi)
  {
    adouble cp;
    R  = zm * RoM_1 + (1.0 - zm) * RoM_2 - dRoM_0 * ( (1.0 - exp(-coeff_1 * zm))                             - pow(zm, coeff_4) * (1.0 - exp(-coeff_1)) );
    cp = zm * Cp_1  + (1.0 - zm) * Cp_2  - dcp_0  * ( (1.0 - exp(-coeff_1 * zm)) + coeff_2 * zm * (1.0 - zm) - pow(zm, coeff_3) * (1.0 - exp(-coeff_1)) );
    h = cp * T;
    gam = cp / (cp - R);
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
  void ComputeProperties_H_AD(adouble &T, adouble &R, adouble &gam, adouble h, adouble zm, adouble zv, adouble chi)
  {
    adouble cp;
    R  = zm * RoM_1 + (1.0 - zm) * RoM_2 - dRoM_0 * ( (1.0 - exp(-coeff_1 * zm))                             - pow(zm, coeff_4) * (1.0 - exp(-coeff_1)) );
    cp = zm * Cp_1  + (1.0 - zm) * Cp_2  - dcp_0  * ( (1.0 - exp(-coeff_1 * zm)) + coeff_2 * zm * (1.0 - zm) - pow(zm, coeff_3) * (1.0 - exp(-coeff_1)) );
    T = h / cp;
    gam = cp / (cp - R);
  }

   void ComputeProperties_vis_AD(adouble &mu, adouble &lamOcp, adouble T)
  {
    if (mu_ref > 0.0)
    {
      mu = calcMuLam_AD(T);
      lamOcp = mu / Pr;
    }
  }


  virtual void diffusivityHookScalarRansComb_AD(const string &name)
  {

    if (name == "ZMean")
    {
 
      // internal faces
      //for (int ifa = nfa_b; ifa < nfa; ifa++)
	 for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) // Ensure we sweep through internal faces of ghosts too
   		if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg))
        		ZMean_diff[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;    


      // boundary faces
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
            ZMean_diff[ifa] = lamOcp_fa[ifa];
        else
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
           ZMean_diff[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;
      }
    }
  }

};


#endif  /* COMBMODEL_BINARYMIXING_H_ */
