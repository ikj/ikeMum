#ifndef COMBMODEL_FPVA_COEFF_AD_H
#define COMBMODEL_FPVA_COEFF_AD_H

#include "../combModels/CombModel_Base.h"
#include "../combModels/CombModel_FPVA_Coeff.h"



template <class Chemtable>
class RansCombFPVA_Coeff_AD : virtual public UgpWithCvCompFlow_AD, public RansCombFPVA_Coeff<ChemtableAdaptiveLinear> {
//class RansCombFPVA_Coeff_AD : virtual public UgpWithCvCompFlow_AD, public RansCombFPVA_Coeff<ChemtableCartesianLinear> {

public:   // member vars
  Chemtable          myChemTable;
  adouble *ZMean, *ZVar, *CMean;
  adouble *ZMean_diff, *ZVar_diff, *CMean_diff;
  adouble (*grad_ZMean)[3], (*grad_ZVar)[3], (*grad_CMean)[3];

public:
  // constructors

  RansCombFPVA_Coeff_AD()
  {
    if (mpi_rank == 0)
      cout << "RansCombFPVA_Coeff_AD()" << endl;

    psi_ZMean = NULL; 	registerScalar(psi_ZMean,"PSI_ZMEAN", CV_DATA); ZMean_Index = getScalarTransportIndex("ZMean");
    psi_ZVar  = NULL; 	registerScalar(psi_ZVar, "PSI_ZVAR",  CV_DATA); ZVar_Index  = getScalarTransportIndex("ZVar");
    psi_CMean = NULL; 	registerScalar(psi_CMean,"PSI_CMEAN", CV_DATA); CMean_Index = getScalarTransportIndex("CMean");


  }

  /*! \brief Unload myMixing and clean memory. */
  virtual ~RansCombFPVA_Coeff_AD() {}



public:
  // member functions 

   virtual void initialize_comb_adjoint()
  {
      if (!checkDataFlag(psi_ZMean)) {
        if(mpi_rank==0) cout<<"INITIALIZING ADJOINT ZMEAN TO ZERO"<<endl;
      for (int icv=0; icv<ncv_gg; icv++)
        psi_ZMean[icv] = 0.;
    }
      if (!checkDataFlag(psi_ZVar)) {
        if(mpi_rank==0) cout<<"INITIALIZING ADJOINT ZVAR TO ZERO"<<endl;
      for (int icv=0; icv<ncv_gg; icv++)
        psi_ZVar[icv] = 0.;
    }
      if (!checkDataFlag(psi_CMean)) {
        if(mpi_rank==0) cout<<"INITIALIZING ADJOINT CMEAN TO ZERO"<<endl;
      for (int icv=0; icv<ncv_gg; icv++)
        psi_CMean[icv] = 0.;
    }
       updateCvDataG1G2(psi_ZMean,  REPLACE_DATA);
       updateCvDataG1G2(psi_ZVar,   REPLACE_DATA);
       updateCvDataG1G2(psi_CMean,  REPLACE_DATA);

       int nScal = scalarTranspEqVector.size();

         for(int i=0; i<nScal ; i++){
                if(strcmp(scalarTranspEqVector_AD[i].name,"ZMean")==0){
                if(mpi_rank==0) cout<<"Initializing  ZMean"<<endl;
                  for (int icv = 0; icv < ncv_gg; icv++) {
                        scalarTranspEqVector_psi[i].phi[icv]=psi_ZMean[icv];
                  }
                }
                if(strcmp(scalarTranspEqVector_AD[i].name,"ZVar")==0){
                if(mpi_rank==0) cout<<"Initializing  ZVar"<<endl;
                  for (int icv = 0; icv < ncv_gg; icv++) {
                        scalarTranspEqVector_psi[i].phi[icv]=psi_ZVar[icv];
                  }
                }
                if(strcmp(scalarTranspEqVector_AD[i].name,"CMean")==0){
                if(mpi_rank==0) cout<<"Initializing  CMean"<<endl;
                  for (int icv = 0; icv < ncv_gg; icv++) {
                        scalarTranspEqVector_psi[i].phi[icv]=psi_CMean[icv];
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
                if(strcmp(scalarTranspEqVector_AD[i].name,"ZVar")==0){
                  for (int icv = 0; icv < ncv_gg; icv++)
                        psi_ZVar[icv]=scalarTranspEqVector_psi[i].phi[icv];
                }
                if(strcmp(scalarTranspEqVector_AD[i].name,"CMean")==0){
                  for (int icv = 0; icv < ncv_gg; icv++)
                        psi_CMean[icv]=scalarTranspEqVector_psi[i].phi[icv];
                }
         }
  }


  virtual void initialHookScalarRansCombModel_AD(int loadtable)
  {
        if(loadtable==1) {
        myChemTable.Load(getStringParam("CHEMTABLE_FILE"));

        Preference = getDoubleParam("REFERENCE_PRESSURE", "0.0");
    	if (Preference == 0.0)
      	Preference = myChemTable.GetReferencePressure();
        }

     // connect pointers 

    int nScal = scalarTranspEqVector.size();

         for(int i=0; i<nScal ; i++){
                if(strcmp(scalarTranspEqVector_AD[i].name,"ZMean")==0){
                        ZMean      = scalarTranspEqVector_AD[i].phi ;
                        grad_ZMean = scalarTranspEqVector_AD[i].grad_phi ;
                        ZMean_diff = scalarTranspEqVector_AD[i].diff ;
                if(mpi_rank==0) cout<<"Connected Scalar Pointer "<<scalarTranspEqVector_AD[i].name<<endl;
                }
                if(strcmp(scalarTranspEqVector_AD[i].name,"ZVar")==0){
                        ZVar      = scalarTranspEqVector_AD[i].phi ;
                        grad_ZVar = scalarTranspEqVector_AD[i].grad_phi ;
                        ZVar_diff = scalarTranspEqVector_AD[i].diff ;
                if(mpi_rank==0) cout<<"Connected Scalar Pointer "<<scalarTranspEqVector_AD[i].name<<endl;
                }
                if(strcmp(scalarTranspEqVector_AD[i].name,"CMean")==0){
                        CMean      = scalarTranspEqVector_AD[i].phi ;
                        grad_CMean = scalarTranspEqVector_AD[i].grad_phi ;
                        CMean_diff = scalarTranspEqVector_AD[i].diff ;
                if(mpi_rank==0) cout<<"Connected Scalar Pointer "<<scalarTranspEqVector_AD[i].name<<endl;
                }
         }

  }


  virtual void calcStateVariables_AD(adouble *rho, adouble (*rhou)[3], adouble *rhoE)
  {
    adouble E, cp;
    adouble kinecv = 0.0;
	  
    ComputeScalarDissipation_AD(chi_AD);
	  
    for (int icv = 0; icv < ncv_gg; icv++)
    {

      // Read and interpolate coefficients from chemistry table and save source term for progress variable
      myChemTable.LookupCoeff_AD(RoM[icv], T0, E0, GAMMA0, AGAMMA, MU0, AMU, LOC0, ALOC, CmeanSource_AD[icv], ZMean[icv], ZVar[icv], CMean[icv]);
      
      for (int i=0; i<3; i++)
        vel[icv][i] = rhou[icv][i]/rho[icv];

      if (kine != NULL)
        kinecv = kine[icv];
      
      E = rhoE[icv] / rho[icv] - 0.5 * vecDotVec3d_AD(vel[icv], vel[icv]) - kinecv;

      // Compute mixture temperature from energy and species mass fractions
      if (AGAMMA == 0.0)
        temp[icv] = T0 + (GAMMA0 - 1.0) / RoM[icv] * (E - E0);
      else
        temp[icv] = T0 + (GAMMA0 - 1.0) / AGAMMA * (exp(AGAMMA * (E - E0) / RoM[icv]) - 1.0);
      
      // Temperature clipping
      temp[icv] = max(temp[icv], Tminimum);

      // Compute mixture enthalpy and pressure
      enthalpy[icv] = E + RoM[icv] * temp[icv];
      press[icv] = rho[icv] * RoM[icv] * temp[icv];
          
      // Compute mixture viscosity, heat conductivity, gamma, speed of sound
      muLam_AD[icv] = MU0 * pow(temp[icv] / T0, 0.7);
      LambdaOverCp_AD[icv] = LOC0 * pow(temp[icv] / T0, 0.62);
      gamma[icv] = GAMMA0 + AGAMMA * (temp[icv] - T0);
      sos[icv] = sqrt(gamma[icv] * press[icv] / rho[icv]);

      if ((isnan(press[icv].value()) || (press[icv]<= 0.0) ))
      {
        cout << "WARNING! :" << endl;
        cout << "x = " << x_cv[icv][0] << " / " << x_cv[icv][1]  << " / " << x_cv[icv][2] << endl;
        cout << "press = " << press[icv] << "  rho = " << rho[icv] << "  temp = " << temp[icv] << "  RoM = " << RoM[icv] << "  E = " << E <<
        "  rhoE = " << rhoE[icv] << "  1/2u^2 = " << 0.5 * vecDotVec3d_AD(vel[icv], vel[icv]) << "  Z = " << ZMean[icv] << endl;
        throw(-1);
      }
       
     

      
      
#ifdef PRESSURE_SCALING_FPVA_COEFF
      // rescale source term
      CmeanSource_AD[icv] = CmeanSource_AD[icv]*pow( press[icv]/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
#endif      

    }
	  
	  
  }
  
#ifdef MIXING_VISCOSITY  
  /*! \brief Compute for each face the material properties.
   *
   *  Compute laminar viscosity mul_fa and heat conductivity lamOcp_fa at the faces.
   *  
   */
  virtual void calcMaterialProperties_AD(adouble *rho, adouble (*rhou)[3], adouble *rhoE)
  {
    if (mu_ref > 0.0)
    {
      // internal faces
    for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) // Ensure we sweep through internal faces of ghosts too
        {
        if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg))
        {

        int icv0 = cvofa[ifa][0];
        int icv1 = cvofa[ifa][1];
  
        double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
        vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
        vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
        double w0 = sqrt(vecDotVec3d(dx0, dx0));
        double w1 = sqrt(vecDotVec3d(dx1, dx1));
 
        mul_fa[ifa] = (w1*muLam_AD[icv0] + w0*muLam_AD[icv1])/(w0+w1);
        lamOcp_fa[ifa] = (w1*LambdaOverCp_AD[icv0] + w0*LambdaOverCp_AD[icv1])/(w0+w1);
      }    
      
      // boundary faces computed in setBC
    }
  }
 }
#endif
  
  
  // \brief Compute for a given temperature, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_T_AD(adouble &p, adouble &h, adouble &R, adouble &gam, adouble &rho, adouble &T, adouble *Scal, int nScal)
  {
    adouble zm = Scal[ZMean_Index];
    adouble zv = Scal[ZVar_Index];
    adouble cm = Scal[CMean_Index];
    
    myChemTable.LookupSelectedCoeff_AD(R, T0, E0, GAMMA0, AGAMMA, zm, zv, cm);
    if (AGAMMA == 0.0)
      h = E0 + R / (GAMMA0 - 1.0) * (T - T0) + R * T;
    else
      h = E0 + R / AGAMMA * log(1.0 + AGAMMA * (T - T0) / (GAMMA0 - 1.0)) + R * T;
    gam = GAMMA0 + AGAMMA * (T - T0);
    p = rho * R * T;
  }
 
  // \brief Compute for a given pressure, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_p_AD(adouble &T, adouble &h, adouble &R, adouble &gam, adouble &rho, adouble &p, adouble *Scal, int nScal)
  {
    adouble zm = Scal[ZMean_Index];
    adouble zv = Scal[ZVar_Index];
    adouble cm = Scal[CMean_Index];
    
    myChemTable.LookupSelectedCoeff_AD(R, T0, E0, GAMMA0, AGAMMA, zm, zv, cm);
    T = p / (rho * R);
    if (AGAMMA == 0.0)
      h = E0 + R / (GAMMA0 - 1.0) * (T - T0) + R * T;
    else    
      h = E0 + R / AGAMMA * log(1.0 + AGAMMA * (T - T0) / (GAMMA0 - 1.0)) + R * T;
    gam = GAMMA0 + AGAMMA * (T - T0);
  }
  
  virtual adouble calcMuLam_AD(int icv)
  {
    return muLam_AD[icv];
  }
  
  virtual adouble calcMuLam_AD(adouble temp)
  {
    cout << "### RansCombFPVA::calcMuLam(double temp) does not work! ####" << endl;
    throw(-1);
  }
  
  // \brief Compute for a given temperature the properties of the mixture at a face.
  virtual void ComputeBCProperties_T_AD(FaZone *zone)
  {

    for (int index = 0; index < zone->faVec.size(); ++index) {
      int ifa = zone->faVec[index];
      int icv1 = cvofa[ifa][1];
      ComputeProperties_T_AD(enthalpy[icv1],RoM[icv1],gamma[icv1],temp[icv1],ZMean[icv1],ZVar[icv1],CMean[icv1]);
    }
   if (mu_ref > 0.0)
    {
     // for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
     for (int index = 0; index < zone->faVec.size(); ++index) {
      int ifa = zone->faVec[index];
      int icv1 = cvofa[ifa][1];
      ComputeProperties_vis(mul_fa[ifa],lamOcp_fa[ifa],temp[icv1],ZMean[icv1],ZVar[icv1],CMean[icv1]);
      }
    }
  }
  
  // \brief Compute for a given enthalpy the properties of the mixture at a face.
  virtual void ComputeBCProperties_H_AD(FaZone *zone)
  {    
    for (int index = 0; index < zone->faVec.size(); ++index) {
      int ifa = zone->faVec[index];
      int icv1 = cvofa[ifa][1];
      ComputeProperties_H_AD(temp[icv1],RoM[icv1],gamma[icv1],enthalpy[icv1],ZMean[icv1],ZVar[icv1],CMean[icv1]);
    }
   if (mu_ref > 0.0)
    {
     // for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
     for (int index = 0; index < zone->faVec.size(); ++index) {
      int ifa = zone->faVec[index];
      int icv1 = cvofa[ifa][1];
      ComputeProperties_vis(mul_fa[ifa],lamOcp_fa[ifa],temp[icv1],ZMean[icv1],ZVar[icv1],CMean[icv1]);
      }
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
  void ComputeProperties_T_AD(adouble &h, adouble &RoM, adouble &gam, adouble T, adouble zm, adouble zv, adouble cm)
  {
    adouble dummy;

    myChemTable.LookupCoeff_AD(RoM, T0, E0, GAMMA0, AGAMMA, MU0, AMU, LOC0, ALOC, dummy, zm, zv, cm);
    if (AGAMMA == 0.0)
      h = E0 + RoM / (GAMMA0 - 1.0) * (T - T0) + RoM * T;
    else    
      h = E0 + RoM / AGAMMA * log(1.0 + AGAMMA * (T - T0) / (GAMMA0 - 1.0)) + RoM * T;
    gam = GAMMA0 + AGAMMA * (T - T0);
    //mu = MU0 * pow(T / T0, 0.7);
    //lamOcp = LOC0 * pow(T / T0, 0.62);
  }

  void ComputeProperties_H_AD(adouble &T, adouble &RoM, adouble &gam, adouble h, adouble zm, adouble zv, adouble cm)
  {
    adouble dummy;
	
    myChemTable.LookupCoeff_AD(RoM, T0, E0, GAMMA0, AGAMMA, MU0, AMU, LOC0, ALOC, dummy, zm, zv, cm);
    T = (h - E0) / RoM;
    T = SolveNewtonTemperature_H_AD(RoM, T0, E0, GAMMA0, AGAMMA, h, T);
    T = max(T, Tminimum);
    gam = GAMMA0 + AGAMMA * (T - T0);
    // mu = MU0 * pow(T / T0, 0.7);
    // lamOcp = LOC0 * pow(T / T0, 0.62);

  }

    adouble SolveNewtonTemperature_H_AD(adouble R, adouble t0, adouble e0, adouble gamma0, adouble agamma, adouble h, adouble Tguess)
  {
    int    itermax = 50, iter;
    double errormax = 1.0e-7;
    adouble Tresidual, hh, cpp, error;

    iter = 0;
    do
    {
      ++iter;
      if (agamma == 0.0)
        hh = e0 + R / (gamma0 - 1.0) * (Tguess - t0) + R * Tguess;
      else
        hh = e0 + R / agamma * log(1.0 + agamma * (Tguess - t0) / (gamma0 - 1.0)) + R * Tguess;
      cpp = R * (1.0 + 1.0 / (gamma0 + agamma * (Tguess - t0) - 1.0));

      Tresidual = Tguess - (hh - h) / cpp;

      error = Tguess - Tresidual;
      if (fabs(error) <= errormax)
        return Tresidual;

      Tguess = Tresidual;

    }
    while (iter < itermax);

    /* if the iteration has not converged, exit with error */
    cerr << "### Computation of temperature in Newton iteration has not converged: Tguess=" << Tguess + error
         << ", Tresidual=" << Tresidual << ", enthalpy=" << h << ", #iterations=" << iter << " ###" << endl;
    throw(-1);
  }


  void ComputeProperties_vis(adouble &mu, adouble &lamOcp, adouble T, adouble zm, adouble zv, adouble cm)
  {
    adouble dummy;
    adouble h,RoM,gam;

    myChemTable.LookupCoeff_AD(RoM, T0, E0, GAMMA0, AGAMMA, MU0, AMU, LOC0, ALOC, dummy, zm, zv, cm);
    if (AGAMMA == 0.0)
      h = E0 + RoM / (GAMMA0 - 1.0) * (T - T0) + RoM * T;
    else    
      h = E0 + RoM / AGAMMA * log(1.0 + AGAMMA * (T - T0) / (GAMMA0 - 1.0)) + RoM * T;
    gam = GAMMA0 + AGAMMA * (T - T0);
    mu = MU0 * pow(T / T0, 0.7);
    lamOcp = LOC0 * pow(T / T0, 0.62);
  }


  void ComputeScalarDissipation_AD(adouble *chi_AD)
  {
    if ((turbModel == KOM) || (turbModel == KOMSST))
    {
      int om_index = getScalarTransportIndex("omega");
      for (int icv = 0; icv < ncv_gg; icv++)
        chi_AD[icv] = Cchi * Cmu * scalarTranspEqVector_AD[om_index].phi[icv];
    }
    
    else if (turbModel == KEPS)
    {
      int epsi_index = getScalarTransportIndex("eps");
      for (int icv = 0; icv < ncv_gg; icv++)
        chi_AD[icv] = Cchi * scalarTranspEqVector_AD[epsi_index].phi[icv] / kine[icv];
    }
    
    else if (turbModel == SA)
    {
      calcVorticity();
      for (int icv = 0; icv < ncv_gg; icv++)
        chi_AD[icv] = Cchi * sqrt(Cmu) * vortMag[icv];
    }
    
    else
    {
      for (int icv = 0; icv < ncv_gg; icv++)
        chi_AD[icv] = 0.0;
    }
    
  }

  virtual void diffusivityHookScalarRansComb_AD(const string &name)
  {

    if (name == "ZMean")
    {
 
      // internal faces
        for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) // Ensure we sweep through internal faces of ghosts too
                if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg))
        ZMean_diff[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;    

      // boundary faces
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))
        for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
            ZMean_diff[ifa] = lamOcp_fa[ifa];
         }
        else
        for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
            ZMean_diff[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;
         }
      }
    }
    if (name == "ZVar")
    {
 
      // internal faces
        for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) // Ensure we sweep through internal faces of ghosts too
                if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg))
        ZVar_diff[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;    

      // boundary faces
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))
        for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
            ZVar_diff[ifa] = lamOcp_fa[ifa];
        }
        else
        for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
            ZVar_diff[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;
        }
      }
    }
    if (name == "CMean")
    {
 
      // internal faces
        for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) // Ensure we sweep through internal faces of ghosts too
                if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg))
        CMean_diff[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;    

      // boundary faces
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))
        for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
            CMean_diff[ifa] = lamOcp_fa[ifa];
        }
        else
        for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
            CMean_diff[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;
        }
      }
    }
  }
  
	virtual void sourceHookScalarRansComb_new_AD(adouble *rhs, double *A, const string &name, int flagImplicit)
	{
		if (name == "ZVar")
		{
			for (int icv = 0; icv < ncv_g; icv++)
			{
				adouble src = 2.0*InterpolateAtCellCenterFromFaceValues(mut_fa, icv)/Schmidt_turb_Z2 * vecDotVec3d_AD(grad_ZMean[icv], grad_ZMean[icv]) - chi_AD[icv] * rho_AD[icv] * ZVar[icv];
				rhs[icv] += src * cv_volume[icv];
			}
			
			if (flagImplicit)
			{
				for (int icv = 0; icv < ncv; icv++)
				{
					int noc00 = nbocv_i[icv];
					double dsrcdphi = - chi_AD[icv].value();
					A[noc00] -= dsrcdphi*cv_volume[icv];
				}
			}
		}
		
		if (name == "CMean")
		{
			for (int icv = 0; icv < ncv_g; icv++)
			{
				rhs[icv] += rho_AD[icv] * CmeanSource_AD[icv] * cv_volume[icv];
			}
			
			if (flagImplicit) {
				// Derivative with respect to CMean
				for (int icv = 0; icv < ncv; icv++) {
					adouble E = enthalpy[icv] - RoM[icv] * temp[icv];
					double E1 = E.value();
      					double ZM = ZMean[icv].value();
      					double ZV = ZMean[icv].value();
					double delta_CM = 1.0e-5;
					double CMp = CMean[icv].value() + delta_CM;
					double CMm = CMean[icv].value() - delta_CM;
					double pressp, pressm, CmeanSourcep, CmeanSourcem;
					pressure_scalSource_AD(pressp, CmeanSourcep, rho[icv], E1, ZM, ZV, CMp);
					pressure_scalSource_AD(pressm, CmeanSourcem, rho[icv], E1, ZM, ZV, CMm);
#ifdef PRESSURE_SCALING_FPVA_COEFF
					CmeanSourcep *= pow( pressp/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
					CmeanSourcem *= pow( pressm/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
#endif      
					int noc00 = nbocv_i[icv];
					A[noc00] -= (CmeanSourcep - CmeanSourcem) / (2.0 * delta_CM) * cv_volume[icv] ;
				}
			}
		}
	}
	
	
	
  virtual void sourceHookRansCombCoupled_AD(adouble **rhs, double ***A, int flagImplicit)
  {

    int ZMean_Coupled_Index = getScalarTransportIndex("ZMean");
    if (ZMean_Coupled_Index == -1)
    {
      cerr << "### Error ZMean_Coupled_Index is " << ZMean_Coupled_Index << " !###" << endl;
      throw(-1);
    }
    int ZVar_Coupled_Index = getScalarTransportIndex("ZVar");
    if (ZVar_Coupled_Index == -1)
    {
      cerr << "### Error ZVar_Coupled_Index is " << ZVar_Coupled_Index << " !###" << endl;
      throw(-1);
    }
    int CMean_Coupled_Index = getScalarTransportIndex("CMean");
    if (CMean_Coupled_Index == -1)
    {
      cerr << "### Error CMean_Coupled_Index is " << CMean_Coupled_Index << " !###" << endl;
      throw(-1);
    }
    
    for (int icv = 0; icv < ncv_g; icv++)
    {
      rhs[icv][5+ZVar_Coupled_Index ] += (2.0 * InterpolateAtCellCenterFromFaceValues(mut_fa, icv) / Schmidt_turb_Z2 * vecDotVec3d_AD(grad_ZMean[icv], grad_ZMean[icv]) - chi_AD[icv] * rho_AD[icv] * ZVar[icv]) * cv_volume[icv];
      rhs[icv][5+CMean_Coupled_Index] += rho_AD[icv] * CmeanSource_AD[icv] * cv_volume[icv];

      if (flagImplicit && icv<ncv)
      {


// ********************* Derivative computation ***********************

      double dCmeanSource_dZM, dCmeanSource_dZV, dCmeanSource_dCM,E;
      double delta_ZM = 1.0e-5;
      double delta_ZV = 1.0e-8;
      double delta_CM = 1.0e-5;
      double pressp, pressm, CmeanSourcep, CmeanSourcem;
      double ZV=ZVar[icv].value();
      double ZM=ZMean[icv].value();
      double CM=CMean[icv].value();

        E = enthalpy[icv].value() - RoM[icv].value() * temp[icv].value();


        // Derivative with respect to ZMean
        double ZMp = ZMean[icv].value() + delta_ZM;
        double ZMm = ZMean[icv].value() - delta_ZM;

        pressure_scalSource_AD(pressp, CmeanSourcep, rho[icv], E, ZMp, ZV, CM);
        pressure_scalSource_AD(pressm, CmeanSourcem, rho[icv], E, ZMm, ZV, CM);
#ifdef PRESSURE_SCALING_FPVA_COEFF
        CmeanSourcep *= pow( pressp/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
        CmeanSourcem *= pow( pressm/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
#endif
        dCmeanSource_dZM = (CmeanSourcep - CmeanSourcem) / (2.0 * delta_ZM);

        // Derivative with respect to ZVar
        double ZVp = ZVar[icv].value() + delta_ZV;
        double ZVm = ZVar[icv].value() - delta_ZV;
        pressure_scalSource_AD(pressp, CmeanSourcep, rho[icv], E, ZM, ZVp, CM);
        pressure_scalSource_AD(pressm, CmeanSourcem, rho[icv], E, ZM, ZVm, CM);
#ifdef PRESSURE_SCALING_FPVA_COEFF
        CmeanSourcep *= pow( pressp/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
        CmeanSourcem *= pow( pressm/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
#endif
        dCmeanSource_dZV = (CmeanSourcep - CmeanSourcem) / (2.0 * delta_ZV);

        // Derivative with respect to CMean
        double CMp = CM + delta_CM;
        double CMm = CM - delta_CM;
        pressure_scalSource_AD(pressp, CmeanSourcep, rho[icv], E, ZM, ZV, CMp);
        pressure_scalSource_AD(pressm, CmeanSourcem, rho[icv], E, ZM, ZV, CMm);
#ifdef PRESSURE_SCALING_FPVA_COEFF
        CmeanSourcep *= pow( pressp/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
        CmeanSourcem *= pow( pressm/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
#endif
        dCmeanSource_dCM = (CmeanSourcep - CmeanSourcem) / (2.0 * delta_CM);


// ********************* Derivative computation ***********************


        int noc00 = nbocv_i[icv];
        A[noc00][5+ZVar_Coupled_Index ][5+ZVar_Coupled_Index ] -= - chi_AD[icv].value() * cv_volume[icv];
        A[noc00][5+ZMean_Coupled_Index][5+CMean_Coupled_Index] -= dCmeanSource_dZM * cv_volume[icv];
        A[noc00][5+ZVar_Coupled_Index][5+CMean_Coupled_Index ] -= dCmeanSource_dZV * cv_volume[icv];
        A[noc00][5+CMean_Coupled_Index][5+CMean_Coupled_Index] -= dCmeanSource_dCM * cv_volume[icv];
        A[noc00][0][5+CMean_Coupled_Index] -= (CmeanSource_AD[icv].value() - ZMean[icv].value() * dCmeanSource_dZM - ZVar[icv].value() * dCmeanSource_dZV - CMean[icv].value() * dCmeanSource_dCM) * cv_volume[icv];
     }
    }
  }
  

  void pressure_scalSource_AD(double &pp, double &CMSource, double &rrho, double &E, double &ZM, double &ZV, double &CM)
  {
    double r0, T,t0,e0,gamma0,agamma;
    myChemTable.LookupSelectedCoeff(r0, t0, e0, gamma0, agamma, ZM, ZV, CM);
    CMSource = myChemTable.Lookup(ZM, ZV, CM, "SRC_PROG");
    if (agamma == 0.0)
      T = t0 + (gamma0 - 1.0) / r0 * (E - e0);
    else
      T = t0 + (gamma0 - 1.0) / agamma * (exp(agamma * (E - e0) / r0) - 1.0);

    // Temperature clipping
    T = max(T, Tminimum);
    pp = rrho * r0 * T;
  }



private:
  adouble T0;
  adouble E0;
  adouble GAMMA0;
  adouble AGAMMA;
  adouble MU0;
  adouble AMU;
  adouble LOC0;
  adouble ALOC;

};


#endif  /* COMBMODEL_FPVA_COEFF_H_AD */
