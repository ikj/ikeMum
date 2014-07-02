#ifndef COMBMODEL_BINARYMIXING_AD_H
#define COMBMODEL_BINARYMIXING_AD_H

#include "../combModels/CombModel_Base.h"
#include "../combModels/CombModel_Powers.h"


class RansCombPowers_AD : virtual public UgpWithCvCompFlow_AD, public RansCombPowers {

public:
  // constructors

  RansCombPowers_AD()
  {
    if (mpi_rank == 0)
      cout << "RansCombPowers_AD()" << endl;

    psi_ZMean = NULL; 	registerScalar(psi_ZMean,"PSI_ZMEAN", CV_DATA); ZMean_Index = getScalarTransportIndex("ZMean");
 

  }

  /*! \brief Unload myMixing and clean memory. */
  virtual ~RansCombPowers_AD() {}

public:   // member vars

  adouble *ZMean;

public:
  // member functions 

   virtual void initialize_comb_adjoint()
  {
     if (!checkDataFlag(psi_ZMean)) {
        if(mpi_rank==0) cout<<"INITIALIZING ADJOINT ZMEAN TO ZERO"<<endl;
      for (int icv=0; icv<ncv_ggff; icv++)
        psi_ZMean[icv] = 0.;
    }
 
       updateCvDataG1G2(psi_ZMean,  REPLACE_DATA);

       int nScal = scalarTranspEqVector.size();

         for(int i=0; i<nScal ; i++){
                if(strcmp(scalarTranspEqVector_AD[i].name,"ZMean")==0){
                  for (int icv = 0; icv < ncv_gg; icv++) {
                        scalarTranspEqVector_psi[i].phi[icv]=psi_ZMean[icv];
                  }
                cout<<"Initialized Mixing adjoint "<<scalarTranspEqVector_AD[i].name<<endl;
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
                cout<<"Connected Scalar Pointer to Mixing adjoint "<<scalarTranspEqVector_AD[i].name<<endl;
                }
         }

     for (int icv = 0; icv < ncv_gg; icv++) {
      RoM[icv] = R_gas;
      gamma[icv] = GAMMA;
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
      
     adouble pr = (gamma[icv]-1.0)*(rhoE[icv]+ZMean[icv]*q_release
        - 0.5 * vecDotVec3d_AD(rhou[icv], rhou[icv]) / rho[icv] - rho[icv] * kinecv);
      if (pr <= 0.0)
        cout << "State var fail negative pressure at xcv: " << x_cv[icv][0] << ", " << x_cv[icv][1] << ", " << x_cv[icv][2]<<endl;


      press[icv] = pr;
      temp[icv] = pr / (rho[icv] * RoM[icv]);
      enthalpy[icv] = (rhoE[icv]+pr)/rho[icv];
      sos[icv] = sqrt(gamma[icv] * press[icv] / rho[icv]);
   }

  }
  
  // \brief Compute for a given temperature, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_T_AD(adouble &p, adouble &h, adouble &R, adouble &gam, adouble &rho, adouble &T, adouble *Scal, int nScal, int dum1, int dum2)
  {
    R = R_gas;
    gam = GAMMA;
    p = rho * R * T;
    adouble zm = Scal[ZMean_Index];
    h = gam * R / (gam - 1.0) * T - zm*q_release;

  }
 
  // \brief Compute for a given pressure, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_p_AD(adouble &T, adouble &h, adouble &R, adouble &gam, adouble &rho, adouble &p, adouble *Scal, int nScal, int dum1, int dum2)
  {
    R = R_gas;
    gam = GAMMA;
    T = p / (rho * R);
    adouble zm = Scal[ZMean_Index];
    h = gam * R / (gam - 1.0) * T - zm*q_release;

  }
  
  // \brief Compute for a given temperature the properties of the mixture at a face.
  virtual void ComputeBCProperties_T_AD(FaZone *zone)
  {
    for (int index = 0; index < zone->faVec.size(); ++index) {
      int ifa = zone->faVec[index];
      int icv1=cvofa[ifa][1];
      ComputeProperties_T_AD(enthalpy[icv1],RoM[icv1],gamma[icv1],temp[icv1],ZMean[icv1],0.0,0.0);
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
  }
  
  void ComputeProperties_T_AD(adouble &h, adouble &R, adouble &gam, adouble T, adouble zm, adouble zv, adouble chi)
  {
    R = R_gas;
    gam = GAMMA;
    h = gam * R / (gam - 1.0) * T - zm*q_release;
  }

  void ComputeProperties_H_AD(adouble &T, adouble &R, adouble &gam, adouble h, adouble zm, adouble zv, adouble chi)
  {
    R = R_gas;
    gam = GAMMA;
    T = (h+zm*q_release)*(gam-1.)/(gam*R);
  }


  virtual void sourceHookScalarRansComb_new_AD(adouble *rhs, double *A, const string &name, int flagImplicit)
  {


	  /* Model of new2 */
        double a         = 1.173958365311995e+03;
        double a2        = a*a;
        double To        = 423.4159779614325;
        double To2       = To*To;
        double To3       = To2*To;
        double To4       = To3*To;
        double base_fact = exp(a/To);



      for (int icv = 0; icv < ncv_g; icv++)
      {
        if(temp[icv]>310.0)
        {
        //adouble src = alpha*(1.-ZMean[icv])*rho_AD[icv]; //orig
        //adouble src = alpha*(1.-ZMean[icv])*rho_AD[icv]*exp(-293.4895913279988/temp[icv])*2.; //new
        //adouble src = alpha*(1.-ZMean[icv])*rho_AD[icv]*(0.5+0.1*ZMean[icv]*ZMean[icv])*2.; //new1
        adouble src = alpha*(1.-ZMean[icv])*rho_AD[icv]*exp(1.173958365311995e+03/temp[icv])*0.0625; //new2

	  //double fact  = exp(1.173958365311995e+03/temp[icv].value());
          //double fact1 = base_fact*(1.+(temp[icv].value()-To)*(-a/To2));
          //double fact2 = fact1+base_fact*0.5*pow((temp[icv].value()-To),2)*(a2/To4+2*a/To3);

          //adouble src = alpha*(1.-ZMean[icv])*rho_AD[icv]*fact2*0.0625;


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
          double dsrcdphi = - alpha;//*exp(-293.4895913279988/temp[icv])*2.;
          A[noc00] -= dsrcdphi * cv_volume[icv];
         }
        }
      }

  }


};


#endif  /* COMBMODEL_BINARYMIXING_H_ */
