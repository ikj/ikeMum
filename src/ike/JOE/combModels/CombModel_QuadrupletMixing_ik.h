/*
 * CombModel_QuadrupletMixing_ik.h
 *
 *  Created on: Jun 19, 2014
 *      Author: ikj
 */

#ifndef COMBMODEL_QUADRUPLETMIXING_IK_H_
#define COMBMODEL_QUADRUPLETMIXING_IK_H_

#define YO2_DEF    yO2 = YO2_PURE * (PHI*yfuel/YFUEL_PURE - (PHI+1.0)*zm + 1.0 );
#define ROM_DEF    R   = yfuel * RoM_1 + (1.0 - yfuel - yO2 - yprod) * RoM_2 + yO2 * RoM_3 + yprod * RoM_4;
#define CP_GAM_DEF cp  = yfuel * Cp_1  + (1.0 - yfuel - yO2 - yprod) * Cp_2 + yO2 * Cp_3 + yprod * Cp_4;    gam = cp / (cp - R);
#define H_DEF      h   = gam * R * (T - Treference) / (gam - 1.0);
#define T_DEF      T   = (gam - 1.0) * h / (gam * R) + Treference;
#define T_E_DEF    T   = (gam - 1.0) * E / R + gam * Treference;


#include "UgpWithCvCompFlow.h"
#include "CombModel_Base.h"

// #############################################################################################
// ------                                                                                 ------
// ----                         RansCombQuadrupletMixing                                   ----
// ------                                                                                 ------
// #############################################################################################

/*! \brief Quadruple Mixing model based on mixture fraction.
 *
 *  Quadruple species: Fuel, Oxygen, Product, and other molecules (N2, Ar, etc.: regarded as one species)
 *  The species mass fractions of Fuel and Oxygen are functions of the mean scalar mixture fraction (Zmean)
 *  and the Fuel mass fraction through specific mixing rules.
 *            s*Yf - Yo + Yo,0     1            Yf     Yo
 *    ZMean = ---------------- = ----- * ( phi*---- - ---- + 1 )  (Note that ZMean=0 if air only, and Zmean=1 if fuel only)
 *             s*Yf,0 - Yo,0     phi+1         Yf,0   Yo,0
 *
 *    We have three transport equations: one for rho*Yf, one for rho*ZMean, and the other for rho*Yp
 *    1) Amount of fuel (rho*Yf) can be obtained from the transport equation of Yf.
 *    2) Then amount of oxygen (rho*Yo) can be calculated using Yf and Zmean:
 *                            Yf
 *         Yo = Yo,0 * ( phi*---- - (phi+1)*z + 1 )
 *                           Yf,0
 *    3) Amount of product (rho*Yp) can be obtained from the transport equation of Yp.
 *    4) Finally, Yother = 1 - (Yf + Yo + Yp)
 *
 *  Index: 1st = Fuel (H2)
 *         2nd = Other molecules (N2, Ar, etc.)
 *         3rd = Oxygen (O2)
 *         4th = Water vapor (H2O)
 *
 *  It can represent pure mixing, and viscosity and thermal diffusion use simple temperature dependence
 *  rule. Specific heat coefficients are not a direct function of the temperature but of the mixture fraction.
 *
 *  Fuel reaction rate is not included in this class: You may calculate it with a simple Arrhenius chemistry or
 *                                                    Eddy dissipation Concept (EDC).
 */

class RansCombQuadrupletMixing: public RansCombBase
{
public:
	// member variables
	double *ZMean;                       ///< Mixture fraction Z.
	double *YFuel;                       ///< Mass fraction of fuel (gas 1)
	double *YProd;                       ///< Mass fraction of water vapor (gas 4)

	double YO2_PURE;                     ///< Mass fraction of oxygen in the pure air stream (e.g. Oxygen mass fraction in standard air is 0.2315)
	double YFUEL_PURE;                   ///< Mass fraction of fuel in the pure fuel stream
	double PHI;
//
//	double  stoi_ratio;                   ///< Stoichiometric ratio

	double Sc, ScTurb;                   ///< Schmidt number and turbulent Schmidt number

	double Gamma_1;                      ///< Gamma of gas 1 (Z=1): Fuel
	double Gamma_2;                      ///< Gamma of gas 2      : Others
	double Gamma_3;                      ///< Gamma of gas 3 (Z=0): Oxygen
	double Gamma_4;                      ///< Gamma of gas 4      : Product
	double MolMass_1;                    ///< Molecular mass of gas 1 (Z=1): Fuel
	double MolMass_2;                    ///< Molecular mass of gas 2      : Others
	double MolMass_3;                    ///< Molecular mass of gas 3 (Z=0): Oxygen
	double MolMass_4;                    ///< Molecular mass of gas 4      : Product
	double RoM_1;                        ///< Gas constant of gas 1 (Z=1): Fuel
	double RoM_2;                        ///< Gas constant of gas 2      : Others
	double RoM_3;                        ///< Gas constant of gas 3 (Z=0): Oxygen
	double RoM_4;                        ///< Gas constant of gas 4      : Product
	double Cp_1;                         ///< Constant specific heat capacity of gas 1 (Z=1): Fuel
	double Cp_2;                         ///< Constant specific heat capacity of gas 2      : Others
	double Cp_3;                         ///< Constant specific heat capacity of gas 3 (Z=0): Oxygen
	double Cp_4;                         ///< Constant specific heat capacity of gas 4      : Product
	double rho_1;                        ///< Density of gas 1 (Z=1): Fuel
	double rho_2;                        ///< Density of gas 2      : Others
	double rho_3;                        ///< Density of gas 3 (Z=0): Oxygen
	double rho_4;                        ///< Density of gas 4      : Product

	double muref_1;                      ///< gas 1 (Z=1)
	double muref_2;                      ///< gas 2
	double muref_3;                      ///< gas 3 (Z=0)
	double muref_4;                      ///< gas 4
	double SL_Tref_1;                    ///< gas 1 (Z=1)
	double SL_Tref_2;                    ///< gas 2
	double SL_Tref_3;                    ///< gas 3 (Z=0)
	double SL_Tref_4;                    ///< gas 4
	double SL_Sref_1;                    ///< gas 1 (Z=1)
	double SL_Sref_2;                    ///< gas 2
	double SL_Sref_3;                    ///< gas 3 (Z=0)
	double SL_Sref_4;                    ///< gas 4

	int ZMean_Index;                  ///< Index in scalar vector of ZMean.
	int YFuel_Index;                  ///< Index in scalar vector of YFuel.
	int YProd_Index;                  ///< Index in scalar vector of YProd.

	double Treference;                   ///< Reference temperature for which enthalpy is zero ( H(Treference) = 0 ).
public:
	// constructors

	RansCombQuadrupletMixing() {
		if (mpi_rank == 0)
			cout << "RansCombQuadrupletMixing()" << endl;

		// Scalar transport equation for ZMean
		ScalarTranspEq *eq;
		eq = registerScalarTransport("ZMean", CV_DATA);
		eq->diffTerm = getIntParam("SCAL_DIFF_ZMean", "1");
		eq->phiZero = 1.0e-8;
		eq->phiZeroRel = 1.0e-4;
		eq->phiMaxiter = 100;
		eq->lowerBound = -10.0;
		eq->upperBound = 10.0;
//		eq->lowerBound = 0.0;
//		eq->upperBound = 1.0;
		eq->reconstruction = getStringParam("SCALAR_RECONSTRUCTION", "CONSERVATIVE");

		// Scalar transport equation for YFuel
		ScalarTranspEq *eq2;
		eq2 = registerScalarTransport("YFuel", CV_DATA);
		eq2->diffTerm = getIntParam("SCAL_DIFF_YFuel", "1");
		eq2->phiZero = 1.0e-8;
		eq2->phiZeroRel = 1.0e-4;
		eq2->phiMaxiter = 100;
		eq2->lowerBound = -10.0;
		eq2->upperBound = 10.0;
//		eq2->lowerBound = 0.0;
//		eq2->upperBound = 1.0;
		eq2->reconstruction = getStringParam("SCALAR_RECONSTRUCTION", "CONSERVATIVE");

		// Scalar transport equation for YProd
		ScalarTranspEq *eq3;
		eq3 = registerScalarTransport("YProd", CV_DATA);
		eq3->diffTerm = getIntParam("SCAL_DIFF_YFuel", "1");
		eq3->phiZero = 1.0e-8;
		eq3->phiZeroRel = 1.0e-4;
		eq3->phiMaxiter = 100;
		eq3->lowerBound = -10.0;
		eq3->upperBound = 10.0;
//		eq3->lowerBound = 0.0;
//		eq3->upperBound = 1.0;
		eq3->reconstruction = getStringParam("SCALAR_RECONSTRUCTION", "CONSERVATIVE");


		// Overall constants
		YO2_PURE   = getDoubleParam("YO2_PURE", "0.232");
		YFUEL_PURE = getDoubleParam("YFUEL_PURE", "1.0");

		if(getIntParam("N_PSALC_CONTROL_PARAMS", "0") > 0)
			PHI = getParam("LAMBDA_INITIAL_0") -> getDouble("CONTROL_PARAM0");
		else
			PHI = 0.0;

		Sc     = getDoubleParam("SC",     "1.0");
		ScTurb = getDoubleParam("SCTURB", "1.0");

		if (mpi_rank == 0)
			cout << "** YO2_PURE  =" << YO2_PURE << endl
			     << "** YFUEL_PURE=" << YFUEL_PURE << endl
			     << "** PHI       =" << PHI << endl
			     << "** Schmidt number(nu/D) = "<<Sc<<", Turbulent Schmidt number = " << ScTurb << endl;
	}

	/*! \brief Unload myMixing and clean memory. */
	virtual ~RansCombQuadrupletMixing() {
		if (mpi_rank == 0)
			cout << endl << "~RansCombQuadrupletMixing()" << endl << endl;

		if (scalarTranspEqVector[ZMean_Index].dpress_dphi != NULL) {
			delete [] scalarTranspEqVector[ZMean_Index].dpress_dphi;
			scalarTranspEqVector[ZMean_Index].dpress_dphi = NULL;
		}

		if (scalarTranspEqVector[YFuel_Index].dpress_dphi != NULL) {
			delete [] scalarTranspEqVector[YFuel_Index].dpress_dphi;
			scalarTranspEqVector[YFuel_Index].dpress_dphi = NULL;
		}

		if (scalarTranspEqVector[YProd_Index].dpress_dphi != NULL) {
			delete [] scalarTranspEqVector[YProd_Index].dpress_dphi;
			scalarTranspEqVector[YProd_Index].dpress_dphi = NULL;
		}
	}

public:
	// member functions

	/*! \brief Read thermo input file and load mixing boundary conditions. */
	virtual void initialHookScalarRansCombModel() {
		if(mpi_rank==0)
			cout<<"RansCombQuadrupletMixing::initialHookScalarRansCombModel()"<<endl;

		MolMass_1 = getDoubleParam("MOLMASS_1", "2.016");  // H2
		MolMass_2 = getDoubleParam("MOLMASS_2", "28.17");  // Air except O2 (Note: air is 28.97)
		MolMass_3 = getDoubleParam("MOLMASS_3", "32.00");  // O2
		MolMass_4 = getDoubleParam("MOLMASS_4", "18.02");  // H2O

		Gamma_1 = getDoubleParam("GAMMA_1", "1.4");
		Gamma_2 = getDoubleParam("GAMMA_2", "1.4");
		Gamma_3 = getDoubleParam("GAMMA_3", "1.4");
		Gamma_4 = getDoubleParam("GAMMA_4", "1.33"); // Specific-heat ratio for H2O: Table A-4 of White appendix

		RoM_1 = 8314.4 / MolMass_1;
		RoM_2 = 8314.4 / MolMass_2;
		RoM_3 = 8314.4 / MolMass_3;
		RoM_4 = 8314.4 / MolMass_4;

		Cp_1 = Gamma_1 * RoM_1 / (Gamma_1 - 1.0);
		Cp_2 = Gamma_2 * RoM_2 / (Gamma_2 - 1.0);
		Cp_3 = Gamma_3 * RoM_3 / (Gamma_3 - 1.0);
		Cp_4 = Gamma_4 * RoM_4 / (Gamma_4 - 1.0);

		Treference = getDoubleParam("T_REFERENCE", "0.0");

		ScalarTranspEq *eq;
		eq = getScalarTransportData("ZMean");               ZMean = eq->phi;               ZMean_Index = getScalarTransportIndex("ZMean");
		eq = getScalarTransportData("YFuel");               YFuel = eq->phi;               YFuel_Index = getScalarTransportIndex("YFuel");
		eq = getScalarTransportData("YProd");               YProd = eq->phi;               YProd_Index = getScalarTransportIndex("YProd");

		// For coupled solution
		if (getStringParam("TIME_INTEGRATION") == "BACKWARD_EULER_COUPLED") {
			scalarTranspEqVector[ZMean_Index].dpress_dphi = new double[ncv_g];
			scalarTranspEqVector[YFuel_Index].dpress_dphi = new double[ncv_g];
			scalarTranspEqVector[YProd_Index].dpress_dphi = new double[ncv_g];
		}

		muref_1 = getDoubleParam("MUREF_1", "0.0");
		muref_2 = getDoubleParam("MUREF_2", "0.0");
		muref_3 = getDoubleParam("MUREF_3", "0.0");
		muref_4 = getDoubleParam("MUREF_4", "0.0");
		SL_Tref_1 = getDoubleParam("SL_TREF_1", "1.0");
		SL_Tref_2 = getDoubleParam("SL_TREF_2", "1.0");
		SL_Tref_3 = getDoubleParam("SL_TREF_3", "1.0");
		SL_Tref_4 = getDoubleParam("SL_TREF_4", "1.0");
		SL_Sref_1 = getDoubleParam("SL_SREF_1", "1.0");
		SL_Sref_2 = getDoubleParam("SL_SREF_2", "1.0");
		SL_Sref_3 = getDoubleParam("SL_SREF_3", "1.0");
		SL_Sref_4 = getDoubleParam("SL_SREF_4", "1.0");

		if (mpi_rank == 0) {
			cout << muref_1 << "  MUREF_1"<< endl;
			cout << muref_2 << "  MUREF_2"<< endl;
			cout << muref_3 << "  MUREF_3"<< endl;
			cout << muref_4 << "  MUREF_4"<< endl;
			cout << SL_Tref_1 << "  SL_TREF_1"<< endl;
			cout << SL_Tref_2 << "  SL_TREF_2"<< endl;
			cout << SL_Tref_3 << "  SL_TREF_3"<< endl;
			cout << SL_Tref_4 << "  SL_TREF_4"<< endl;
			cout << SL_Sref_1 << "  SL_SREF_1"<< endl;
			cout << SL_Sref_2 << "  SL_SREF_2"<< endl;
			cout << SL_Sref_3 << "  SL_SREF_3"<< endl;
			cout << SL_Sref_4 << "  SL_SREF_4"<< endl;
		}
		if (mpi_rank == 0)
			cout << endl << "***** Quadruple Mixing initialized *****" << endl << endl;
	}

	/*! \brief Compute for each cell the mixture properties as function of \p Zmean.
	 *
	 *  Compute composition Species::Yk, temperature, viscosity Mixture::mul, heat coefficient Mixture::lambda,
	 *  enthalpy Mixture::H, heat capacity ratio Mixture::gama, pressure and gas constant for the mixture.
	 */
	virtual void calcStateVariables() {
		double yfuel, yprod, yO2, zm, R, E, cp, gam, h, T;
		double kinecv = 0.0;

		for (int icv = 0; icv < ncv_gg; icv++) {
			for (int i=0; i<3; i++)
				vel[icv][i] = rhou[icv][i]/rho[icv];

			if (kine != NULL)
				kinecv = kine[icv];

			// Compute mixture gas constant and specific heat coefficient
			yfuel = YFuel[icv];
			yprod = YProd[icv];
			zm = ZMean[icv];
			YO2_DEF

			ROM_DEF
			RoM[icv] = R;
			CP_GAM_DEF
			gamma[icv] = gam;

			E = rhoE[icv] / rho[icv] - 0.5 * vecDotVec3d(vel[icv], vel[icv]) - kinecv;

			// Compute mixture temperature from energy and species mass fractions
			T_E_DEF
			temp[icv] = T;

			// Compute mixture enthalpy, pressure, gamma and speed of sound
			H_DEF
			enthalpy[icv] = h;
			press[icv] = rho[icv] * RoM[icv] * temp[icv];
			sos[icv] = sqrt(gamma[icv]*press[icv]/rho[icv]);
		}

		// viscosity and thermal diffusivity computed in UgpWithCvCompFlow::calcMaterialProperties() at cell faces
	}

	// \brief Compute derivative of pressure with respect to scalar for coupled solution (Jacobi matrix)
	virtual void pressureDerivativeHookScalarRansComb() {
		cout<<"ERROR in RansCombQuadrupletMixing::pressureDerivativeHookScalarRansComb(): This method has not been implemented!!"<<endl;
		throw(-1);
//		double dcp_dscal = Cp_1 - Cp_2;
//		double dR_dscal = RoM_1 - RoM_2;
//
//		for (int icv = 0; icv < ncv; icv++) {
//			double E = enthalpy[icv] - RoM[icv] * temp[icv];
//			double dgam_dscal = (gamma[icv] - 1.0) / RoM[icv] * ((1.0 - gamma[icv]) * dcp_dscal + gamma[icv] * dR_dscal);
//			double dRT_dscal = (E + RoM[icv] * Treference) * dgam_dscal + gamma[icv] * Treference * dR_dscal;
//			scalarTranspEqVector[ZMean_Index].dpress_dphi[icv] = rho[icv] * dRT_dscal;
//		}
	}

	// \brief Compute for a given temperature, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
	virtual void calcThermoProp_T(double &p, double &h, double &R, double &gam, double &rho, double &T, double *Scal, int nScal) {
		double cp, yO2;
		double yfuel = Scal[YFuel_Index];
		double zm    = Scal[ZMean_Index];
		double yprod = Scal[YProd_Index];
		YO2_DEF

		ROM_DEF
		CP_GAM_DEF
		H_DEF
		p = rho * R * T;
	}

	// \brief Compute for a given pressure, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
	virtual void calcThermoProp_p(double &T, double &h, double &R, double &gam, double &rho, double &p, double *Scal, int nScal) {
		double cp, yO2;
		double yfuel = Scal[YFuel_Index];
		double zm    = Scal[ZMean_Index];
		double yprod = Scal[YProd_Index];
		YO2_DEF

		ROM_DEF
		CP_GAM_DEF
		T = p / (rho * R);
		H_DEF
	}

	// \brief Compute for a given temperature the properties of the mixture at a face.
	virtual void ComputeBCProperties_T(FaZone *zone) {
		ScalarTranspEq *eq;
		eq = getScalarTransportData("YFuel");         double *yfuel_fa = eq->phi;
		eq = getScalarTransportData("ZMean");         double *zM_fa = eq->phi;
		eq = getScalarTransportData("YProd");         double *yprod_fa = eq->phi;

		for (int index = 0; index < zone->faVec.size(); ++index) {
			int ifa = zone->faVec[index];
			int icv1 = cvofa[ifa][1];
			ComputeProperties_T(enthalpy[icv1], RoM[icv1], gamma[icv1], temp[icv1], yfuel_fa[icv1], yprod_fa[icv1], zM_fa[icv1], 0.0, 0.0);
		}
		for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
			int icv1 = cvofa[ifa][1];
			ComputeProperties_vis(mul_fa[ifa], lamOcp_fa[ifa], temp[icv1], yfuel_fa[icv1], yprod_fa[icv1], zM_fa[icv1]);
		}
	}

	// \brief Compute for a given enthalpy the properties of the mixture at a face.
	virtual void ComputeBCProperties_H(FaZone *zone) {
		ScalarTranspEq *eq;
		eq = getScalarTransportData("YFuel");         double *yfuel_fa = eq->phi;
		eq = getScalarTransportData("ZMean");         double *zM_fa    = eq->phi;
		eq = getScalarTransportData("YProd");         double *yprod_fa = eq->phi;

		for (int index = 0; index < zone->faVec.size(); ++index) {
			int ifa = zone->faVec[index];
			int icv1 = cvofa[ifa][1];
			ComputeProperties_H(temp[icv1], RoM[icv1], gamma[icv1], enthalpy[icv1], yfuel_fa[icv1], yprod_fa[icv1], zM_fa[icv1], 0.0, 0.0);
		}
		for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
			int icv1 = cvofa[ifa][1];
			ComputeProperties_vis(mul_fa[ifa], lamOcp_fa[ifa], temp[icv1], yfuel_fa[icv1], yprod_fa[icv1], zM_fa[icv1]);
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
	void ComputeProperties_T(double &h, double &R, double &gam, double T, double yfuel, double yprod, double zm, double zv, double chi) {
		double cp, yO2;
		YO2_DEF

		ROM_DEF
		CP_GAM_DEF
		H_DEF
	}


	void ComputeProperties_vis(double &mu, double &lamOcp, double T, double yfuel, double yprod, double zm) {
		if (mu_ref > 0.0) {
			double yO2;
			YO2_DEF

			//mu = calcMuLam(T);
			double mu_1 = muref_1*pow((T/SL_Tref_1),1.5)*((SL_Tref_1+SL_Sref_1)/(T+SL_Sref_1));
			double mu_2 = muref_2*pow((T/SL_Tref_2),1.5)*((SL_Tref_2+SL_Sref_2)/(T+SL_Sref_2));
			double mu_3 = muref_3*pow((T/SL_Tref_3),1.5)*((SL_Tref_3+SL_Sref_3)/(T+SL_Sref_3));
			double mu_4 = muref_4*pow((T/SL_Tref_4),1.5)*((SL_Tref_4+SL_Sref_4)/(T+SL_Sref_4));

			double MolMass_tot = 1.0 / (yfuel/MolMass_1 + (1.0-yfuel-yO2-yprod)/MolMass_2 + yO2/MolMass_3 + yprod/MolMass_4);
			double x1 = MolMass_tot / MolMass_1 * yfuel;   // Mole fraction of the fuel
			double x2 = MolMass_tot / MolMass_2 * (1.0 - yfuel - yO2 - yprod);
			double x3 = MolMass_tot / MolMass_1 * yO2;   // Mole fraction of oxygen
			double x4 = 1.0 - x1 - x2 - x3; // Mole fraction of product (H2O)

			assert(x1>=0.0 && x2>=0.0 && x3>=0.0 && x4>=0.0 && x4<=1.0);

			mu = (mu_1*x1*sqrt(MolMass_1) + mu_2*x2*sqrt(MolMass_2) + mu_3*x3*sqrt(MolMass_3) + mu_4*x4*sqrt(MolMass_4)) / (x1*sqrt(MolMass_1) + x2*sqrt(MolMass_2) + x3*sqrt(MolMass_3) + x4*sqrt(MolMass_4));

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
	 *  \param[in]  yfuel. Mass fraction of the fuel.
	 *  \param[in]  zv Variance of mixture fraction. !!!
	 *  \param[in]  chi Scalar dissipation rate. !!!
	 */
	void ComputeProperties_H(double &T, double &R, double &gam, double h, double yfuel, double yprod, double zm, double zv, double chi) {
		double yO2;
		YO2_DEF

		double cp;
		ROM_DEF
		CP_GAM_DEF
		T_DEF
	}


	virtual void diffusivityHookScalarRansComb(const string &name) {
		ScalarTranspEq *eq;

		if (name.compare("ZMean")==0 || name.compare("YFuel")==0 || name.compare("YProd")==0) {
			eq = getScalarTransportData(name);

			// internal faces
			for (int ifa = nfa_b; ifa < nfa; ifa++)
				eq->diff[ifa] = mul_fa[ifa]/Sc + mut_fa[ifa]/ScTurb; // Note: Schmidt number = nu/D

			// boundary faces
			for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
				if (zone->getKind() == FA_ZONE_BOUNDARY) {
					if (zoneIsWall(zone->getName()))
						for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
							eq->diff[ifa] = mul_fa[ifa]/Sc;
					else
						for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
							eq->diff[ifa] = mul_fa[ifa]/Sc + mut_fa[ifa]/ScTurb;
				}
		}
	}

	/*
	 * Method: calcMuLam
	 * -----------------
	 * Originally defined in UgpWithCvCompFlow::calcMuLam(double)
	 * Code adopted from calcMuLam() in CombModel_BinaryMixing_new.h
	 */
	double calcMuLam(int icv) {
		double T = temp[icv];
		double yfuel = YFuel[icv];
		double zm    = ZMean[icv];
		double yprod = YProd[icv];
		double yO2;
		YO2_DEF

		double mu_1 = muref_1*pow((T/SL_Tref_1),1.5)*((SL_Tref_1+SL_Sref_1)/(T+SL_Sref_1));
		double mu_2 = muref_2*pow((T/SL_Tref_2),1.5)*((SL_Tref_2+SL_Sref_2)/(T+SL_Sref_2));
		double mu_3 = muref_3*pow((T/SL_Tref_3),1.5)*((SL_Tref_3+SL_Sref_3)/(T+SL_Sref_3));
		double mu_4 = muref_4*pow((T/SL_Tref_4),1.5)*((SL_Tref_4+SL_Sref_4)/(T+SL_Sref_4));

		double mu = (mu_1*yfuel*pow(MolMass_1,0.5)+mu_2*(1-yfuel-yO2-yprod)*pow(MolMass_2,0.5)+mu_3*yfuel*pow(MolMass_3,0.5)+mu_4*yfuel*pow(MolMass_4,0.5)) / (yfuel*pow(MolMass_1,0.5)+(1-yfuel-yO2-yprod)*pow(MolMass_2,0.5)+yO2*pow(MolMass_3,0.5)+yprod*pow(MolMass_4,0.5));
		return mu;
	}
};

#endif /* COMBMODEL_QUADRUPLETMIXING_IK_H_ */
