/*
 * CombModel_Powers.h
 *
 *  Created on: Dec 22, 2010
 *      Author: Karthik Duraisamy
 *      Update: Dec 22, 2010
 */

#ifndef COMBMODEL_POWERS_H_
#define COMBMODEL_POWERS_H_

#include "UgpWithCvCompFlow.h"
#include "CombModel_Base.h"


// #############################################################################################
// ------                                                                                 ------
// ----                         RansCombPowers                                        ----
// ------                                                                                 ------
// #############################################################################################

/*! \brief Combustion model based on mixture fraction, with simple Arrhenius chemistry
 *
 */
class RansCombPowers: public RansCombBase  {
public:
  // member variables

  double  *ZMean;                      ///< Mixture fraction Z.
  int     ZMean_Index;                 ///< Index in scalar vector of ZMean.
  double  q_release;                   ///< Heat release.
  double  alpha;                       ///< alpha.
  

public:
  // constructors

  RansCombPowers()
  {
    if (mpi_rank == 0)
      cout << "RansCombPowers()" << endl;

    q_release = getDoubleParam("q_release", "3.e5");
    alpha     = getDoubleParam("alpha", "1000");
  
    // Scalar transport equation for ZMean
    ScalarTranspEq *eq;
    eq = registerScalarTransport("ZMean", CV_DATA);
    eq->diffTerm = getIntParam("SCAL_DIFF", "0");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-4;
    eq->phiMaxiter = 100;
    eq->lowerBound = -10.0;
    eq->upperBound = 10.0;
    eq->reconstruction = getStringParam("SCALAR_RECONSTRUCTION", "CONSERVATIVE");
  }

  /*! \brief Unload myMixing and clean memory. */
  virtual ~RansCombPowers()
  {
    if (mpi_rank == 0)
      cout << endl << "***** Powers Mixing finalized *****" << endl << endl;
  }

public:
  // member functions 

  /*! \brief Read thermo input file and load mixing boundary conditions. */
  virtual void initialHookScalarRansCombModel()
  {
    
  
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");           ZMean = eq->phi;           ZMean_Index = getScalarTransportIndex("ZMean");
    
    if (!checkScalarFlag("ZMean"))
     for (int icv = 0; icv < ncv_gg; icv++)
	ZMean[icv]=0.0;

     for (int icv = 0; icv < ncv_gg; icv++) {
      RoM[icv] = R_gas;
      gamma[icv] = GAMMA;
    }


    if (mpi_rank == 0)
      cout << endl << "***** Powers Mixing initialized *****" << endl << endl;
  }

  /*! \brief Compute for each cell the mixture properties as function of \p Zmean.
   *
   *  Compute composition Species::Yk, temperature, viscosity Mixture::mul, heat coefficient Mixture::lambda,
   *  enthalpy Mixture::H, heat capacity ratio Mixture::gama, pressure and gas constant for the mixture.
   */
  virtual void calcStateVariables()
  {
    double zm, R, E, cp, gam, h, T;
    double kinecv = 0.0;

    for (int icv = 0; icv < ncv_gg; icv++)
    {
      for (int i=0; i<3; i++)
        vel[icv][i] = rhou[icv][i]/rho[icv];

      if (kine != NULL)
        kinecv = kine[icv];
      
       double pr = (gamma[icv]-1.0)*(rhoE[icv]+ZMean[icv]*q_release 
        - 0.5 * vecDotVec3d(rhou[icv], rhou[icv]) / rho[icv] - rho[icv] * kinecv);
      if (pr <= 0.0)
        cout << "State var fail negative pressure at xcv: " << x_cv[icv][0] << ", " << x_cv[icv][1] << ", " << x_cv[icv][2] << " "<<rhou[icv][0]<<" "<<rhou[icv][1]<<" "<<rhou[icv][2]<<" "<<ZMean[icv]<<endl;

      press[icv] = pr;
      temp[icv] = pr / (rho[icv] * RoM[icv]);
      enthalpy[icv] = (rhoE[icv]+pr)/rho[icv];
      sos[icv] = sqrt(gamma[icv] * press[icv] / rho[icv]);

   }

    // viscosity and thermal diffusivity computed in UgpWithCvCompFlow::calcMaterialProperties() at cell faces
  }
  
  // \brief Compute for a given temperature, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_T(double &p, double &h, double &R, double &gam, double &rho, double &T, double *Scal, int nScal, int dum1, int dum2)
  {
    R = R_gas;
    gam = GAMMA;
    p = rho * R * T;
    double zm = Scal[ZMean_Index];
    h = gam * R / (gam - 1.0) * T - zm*q_release;
  }
  
  // \brief Compute for a given pressure, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_p(double &T, double &h, double &R, double &gam, double &rho, double &p, double *Scal, int nScal, int dum1, int dum2)
  { 
    R = R_gas;
    gam = GAMMA; 
    T = p / (rho * R); 
    double zm = Scal[ZMean_Index];
    h = gam * R / (gam - 1.0) * T - zm*q_release;
  }

  
  // \brief Compute for a given temperature the properties of the mixture at a face.
  virtual void ComputeBCProperties_T(FaZone *zone)
  {
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");         double *zM_fa = eq->phi;

    double dum;
    for (int index = 0; index < zone->faVec.size(); ++index) {
      int ifa = zone->faVec[index];
      int icv1 = cvofa[ifa][1];
      ComputeProperties_T(enthalpy[icv1],RoM[icv1],gamma[icv1],dum,dum,temp[icv1],zM_fa[icv1],0.0,0.0);
    }
  }
  
  // \brief Compute for a given enthalpy the properties of the mixture at a face.
  virtual void ComputeBCProperties_H(FaZone *zone)
  {
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");         double *zM_fa = eq->phi;
    
    double dum;
    for (int index = 0; index < zone->faVec.size(); ++index) {
      int ifa = zone->faVec[index];
      int icv1 = cvofa[ifa][1];
      ComputeProperties_H(temp[icv1],RoM[icv1],gamma[icv1],dum,dum,enthalpy[icv1],zM_fa[icv1],0.0,0.0);
    }
  }

    /*! \brief Compute for a given temperature the properties of the mixture.
   *
   *  The enthalpy (chemical + sensible) is a function of temperature and mixture composition. 
   *  \param[out] H Enthalpy (chemical + sensible) in [J/kg].
   *  \param[out] RoverM Gas constant R/M of the mixture in [J/(kg K)].
   *  \param[out] gam Heat capacity ratio in [-].
   *  \param[in]  T Temperature in [K].
   *  \param[in]  zm Mean mixture fraction.
   *  \param[in]  zv Variance of mixture fraction.
   *  \param[in]  chi Scalar dissipation rate.
   */

  
     void ComputeProperties_T(double &h, double &R, double &gam, double &mu, double &lamOcp, double T, double zm, double zv, double chi)
  {
    R = R_gas;
    gam = GAMMA;
    h = gam * R / (gam - 1.0) * T - zm*q_release;
  }
    
  /*! \brief Compute for a given enthalpy the properties of the mixture.
   *
   *  The enthalpy (chemical + sensible) is a function of temperature and mixture composition. 
   *  \param[out] T Temperature in [K].
   *  \param[out] RoverM Gas constant R/M of the mixture in [J/(kg K)].
   *  \param[out] gam Heat capacity ratio in [-].
   *  \param[in]  H Enthalpy (chemical + sensible) in [J/kg].
   *  \param[in]  zm Mean mixture fraction.
   *  \param[in]  zv Variance of mixture fraction.
   *  \param[in]  chi Scalar dissipation rate.
   */
  void ComputeProperties_H(double &T, double &R, double &gam, double &mu, double &lamOcp, double h, double zm, double zv, double chi)
  {
    R = R_gas;
    gam = GAMMA;
    T = (h+zm*q_release)*(gam-1.)/(gam*R);
  }

  virtual void sourceHookScalarRansComb_new(double *rhs, double *A, const string &name, int flagImplicit)
  {

	/* Model of new2 */
	double a         = 1.173958365311995e+03;
	double a2        = a*a;
        double To        = 423.4159779614325; 
        double To2       = To*To; 
        double To3       = To2*To; 
        double To4       = To3*To; 
        double base_fact = exp(a/To); 
	double delta_rho = 0.3/16;
	double start_rho = 2.0;
	double delta_temp= 30./16;
	double start_temp= 400.;

      for (int icv = 0; icv < ncv; icv++)
      {
        if(temp[icv]>310.0) 
        {

          //double fact = 1.0 						//orig;
          double fact = exp(a/temp[icv]) ;				//new2;
	  //double fact1 = base_fact*(1.+(temp[icv]-To)*(-a/To2));
	  //double fact2 = fact1+base_fact*0.5*pow((temp[icv]-To),2)*(a2/To4+2*a/To3);
	  /*
          double fact = exp(a/temp[icv])+40./rho[icv]*exp(-pow((rho[icv]-3.1)/0.25,2)-pow((temp[icv]-450.)/20.,2));
	  */


	  // Interpolate
	  /*
	  int index_temp = (temp[icv]-start_temp)/delta_temp;
	  int index_rho = (rho[icv]-start_rho)/delta_rho;
	  double t0 = start_temp + index_temp*delta_temp;
	  double t1 = t0+delta_temp;
	  double tf = (temp[icv]-t0)/delta_temp;
	  double r0 = start_rho + index_rho*delta_rho;
	  double r1 = r0+delta_rho;
	  double rf = (rho[icv]-r0)/delta_rho;
	  double fact_00 = exp(a/t0)+40./r0*exp(-pow((r0-3.1)/0.25,2)-pow((t0-450.)/20.,2));
	  double fact_01 = exp(a/t0)+40./r1*exp(-pow((r1-3.1)/0.25,2)-pow((t0-450.)/20.,2));
	  double fact_10 = exp(a/t1)+40./r0*exp(-pow((r0-3.1)/0.25,2)-pow((t1-450.)/20.,2));
	  double fact_11 = exp(a/t1)+40./r1*exp(-pow((r1-3.1)/0.25,2)-pow((t1-450.)/20.,2));
	  double fact = (1.-tf)*(1.-rf)*fact_00 + (1.-tf)*rf*fact_01 + tf*(1.-rf)*fact_10 + tf*rf*fact_11 ;
	  */

          double src = alpha*(1.-ZMean[icv])*rho[icv]*fact*0.0625; 
           rhs[icv] += src * cv_volume[icv];
        }
      }

      if (flagImplicit)
      {
        for (int icv = 0; icv < ncv; icv++)
        {
         if(temp[icv]>310.0) 
         {
          int noc00 = nbocv_i[icv];
          double dsrcdphi = - alpha;
          A[noc00] -= dsrcdphi * cv_volume[icv];
         }
        }
      }

  }

  virtual void sourceHookRansCombCoupled(double **rhs, double ***A, int nScal, int flagImplicit)
  {
      for (int icv = 0; icv < ncv; icv++)
      {
        if(temp[icv]>0.0) 
        {
        double src = alpha*(1.-ZMean[icv])*rho[icv]*exp(-120.0/temp[icv]);
        rhs[icv][5] += src * cv_volume[icv];
        }
      }

      if (flagImplicit)
      {
        for (int icv = 0; icv < ncv; icv++)
        {
         if(temp[icv]>0.0) 
         {
          int noc00 = nbocv_i[icv];
          double dsrcdphi = - alpha*exp(-120.0/temp[icv]);
          A[noc00][5][5] -= dsrcdphi * cv_volume[icv];
         }
        }
      }

  }

};


#endif  /* COMBMODEL_POWERS_H_ */
