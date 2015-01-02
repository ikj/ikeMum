/*
 * IkeUgpWithCvCompFlow.h
 *
 *  Created on: Jan 15, 2013
 *      Author: ikj
 */

// This file will be included in JOE/ADJOINT_FILES/UgpWithCvCompFlowAD.h

//#ifndef IKEUGPWITHCVCOMPFLOW_H_
//#define IKEUGPWITHCVCOMPFLOW_H_

#ifdef USE_MEM_SAVING_ADVAR
	// =========================
	//  REFERENCE VALUES
	// =========================
	// Note: Actually, the following reference values were already defined in the UgpWithCvCompFlow class:
	//                     viscMode, mu_ref, SL_Tref, SL_Tref,
	//                     T_ref, p_ref, rho_ref,
	//                     GAMMA, R_gas, Pr, PrTurb.
	struct REF_FLOW_PARAMS {
		REF_FLOW_PARAMS(): nScal(0), rho_ref(1.0), press_ref(1.0), temp_ref(1.0), gamma_ref(1.0), Rgas_ref(1.0) {
			for(int i=0; i<3; ++i)
				u_ref[i]    = 1.0;
			for(int i=0; i<3; ++i)
				rhou_ref[i] = 1.0;
			kine_ref = 0.0;
		}
		REF_FLOW_PARAMS(int nScal): rho_ref(1.0), press_ref(1.0), temp_ref(1.0), gamma_ref(1.0), Rgas_ref(1.0) {
			for(int i=0; i<3; ++i)
				u_ref[i]    = 1.0;
			for(int i=0; i<3; ++i)
				rhou_ref[i] = 1.0;
			kine_ref = 0.0;

			scalar_ref.resize(nScal, 1.0);
		}
		~REF_FLOW_PARAMS() {
			clear();
		}

		void clear() {
			nScal = 0;
			scalar_ref.clear();
		}

		void calcRhouRefFromU() {
			for(int i=0; i<3; ++i)
				rhou_ref[i] = rho_ref*u_ref[i];
		}

		void calcVelMagRef() {
			double umagRefSq = 0.0;
			for(int i=0; i<3; ++i)
				umagRefSq += u_ref[i] * u_ref[i];
			uMag_ref = sqrt(umagRefSq);

			if(uMag_ref <= 1.0e-10)
				cout<<"WARNING in REF_FLOW_PARAMS::calcVelMagRef(): uMag_ref is near ZERO  "<<uMag_ref<<endl;

			rhouMag_ref = rho_ref * uMag_ref;
		}

		void calcRhoERef() {
			rhoE_ref = fabs(1/(gamma_ref-1.0)*press_ref) + fabs(0.5*rho_ref*uMag_ref) + fabs(rho_ref*kine_ref);
		}

		void showOnScreen() {
			cout<<endl
				<<"--------------------------------------------------"<<endl
				<<"FLOW REFERENCE VALUES IN THE REF_FLOW_PARAMS STRUCT:"<<endl
				<<"  RHO_REF   = "<<rho_ref<<endl
				<<"  U_REF     = "<<u_ref[0]<<"  "<<u_ref[1]<<"  "<<u_ref[2]<<endl
				<<"  P_REF     = "<<press_ref<<endl
				<<"  T_REF     = "<<temp_ref<<endl
				<<"  RGAS_REF  = "<<Rgas_ref<<endl
				<<"  GAMMA_REF = "<<gamma_ref<<endl
				<<"  KINE_REF  = "<<kine_ref<<endl
				<<"  NSCAL     = "<<nScal<<endl;
			if(nScal>0) {
				cout<<"  SCAL_REF  =";
				for(size_t iScal=0; iScal<scalar_ref.size(); ++iScal)
					cout<<" "<<scalar_ref[iScal]<<" ";
				cout<<endl;
			}
			cout<<"--------------------------------------------------"<<endl
				<<endl;
		}

		int nScal;

		double rho_ref;

		double u_ref[3];
		double uMag_ref;
		double rhou_ref[3];
		double rhouMag_ref;
		double rhoE_ref;

		double press_ref;
		double temp_ref;

		double gamma_ref;
		double Rgas_ref;

		double kine_ref;
		vector<double> scalar_ref;
	};

	REF_FLOW_PARAMS RefFlowParams;

	// =========================
	//  MEMBER VARIABLES
	// =========================
	ADscalar<adouble> rho_AD;       ///< density derivative
	ADvector<adouble> rhou_AD;  ///< momentum derivative
	ADscalar<adouble> rhoE_AD;      ///< energy derivative

	ADvector<adouble> vel;   ///< velocity
	ADscalar<adouble> press;      ///< pressure e
	ADscalar<adouble> temp;       ///< temperature
	ADscalar<adouble> enthalpy;   ///< chemical + sensible enthalpy (no kinetic energy!) derivative

	ADscalar<adouble> gamma;      ///< gamma derivative
	ADscalar<adouble> RoM;        ///< Gas constant derivative
	ADscalar<adouble> sos;        ///< speed of sound derivative

	ADscalar<adouble> mul_fa;           ///< laminar viscosity at cell faces
	ADscalar<adouble> lamOcp_fa;        ///< heat conductivity at cell faces

	ADscalar<REALQS> mut_fa;     ///< turbulent viscosity derivative
	ADscalar<REALQS>* kine;      ///< turbulent kinetic energy derivative

	ADscalar<adouble> strMag;
	ADscalar<adouble> vortMag;
	ADscalar<adouble> diverg;  ///< magnitude strain rate, divergence derivatives
	ADscalar<adouble> blendFuncF1_AD;
	ADscalar<adouble> blendFuncF2_AD;
	ADscalar<adouble> crossDiff_AD; ///<SST model functions

	ADscalar<adouble> CmeanSource_AD;
	ADscalar<adouble> muLam_AD;
	ADscalar<adouble> LambdaOverCp_AD;
	ADscalar<adouble> chi_AD; ///<FPVA model functions

	// ----------------------------------------------------------------------------------------
	// boundary faces, allocate memory, TODO: should be incorporated in to the associated arrays ???
	// ----------------------------------------------------------------------------------------

	ADvector<adouble> grad_rho;      ///< density  gradient der
	ADtensor<adouble> grad_u;     ///< velocity gradient tensor der
	ADvector<adouble> grad_p;        ///< pressure gradient der
	ADvector<adouble> grad_temp;     ///< temperature gradient der
	ADvector<adouble> grad_enthalpy; ///< enthalpy gradient der

	// Two-layered CSR struct: It used to be in IkeWithModels but was moved to here (Dec 2014)
	int *nbocv2_i;
	vector<int> nbocv2_v;

	// For the memory saving ADvars ("1D-style")
	vector<int> nbocv2_eachIcv;   // The 2-layered neighbors of an icv
	vector<int> nbocv2ff_eachIcv; // The 2-layered neighbors of an icv including the fake cells
	vector<int> fa2_eachIcv;      // The 2-layered faces of an icv

#ifdef USE_ARTIF_VISC_WITH_MEM_SAVING
	// --------------------
	// Artificial viscosity
	// --------------------
	// Artificial viscosity for shock-capturing: cell-centered value
	ADscalar<adouble> artifViscMag_AD;

	// Statistics for artificial viscosity: These values are updated in the IkeWithModel class
	double myArtifViscMagMin_AD;
	double myArtifViscMagMax_AD;
	int myArtifViscSavedStep;
#endif

	class ScalarTranspEq_AD {
	public:
		ScalarTranspEq_AD() {
		}
		~ScalarTranspEq_AD() {
			phi.clear();
			diff.clear();
			grad_phi.clear();
			rhophi.clear();
			grad_rhophi.clear();
			dpress_dphi.clear();
		}

		ADscalar<adouble> phi;
		ADscalar<adouble> diff;
		ADvector<adouble> grad_phi;
		ADscalar<adouble> rhophi;
		ADvector<adouble> grad_rhophi;
		ADscalar<adouble> dpress_dphi;
		char name[DATA_NAME_LEN];
	};

	/*
	 * If the memory saving mode (so-called the 1D-style) is not used:
	 *     If control parameters are used (i.e. NcontrolEqns != 0),
	 *     the IkeWithPsALC_AD class will call IkeWithPsALC_AD::initialize_adoubles().
	 *     That method first assigns an adouble variable for the control parameter
	 *     and call this method.
	 * Otherwise (1D-style):
	 *     The IkeWithPsALC_AD class will call IkeWithPsALC_AD::initialize_adoubles(const int icvCenter, const int NcontrolEqns).
	 */
	virtual void initialize_adoubles() {
		rho_AD.allocate(ncv_ggff);
		rhou_AD.allocate(ncv_ggff);
		rhoE_AD.allocate(ncv_ggff);

		vel.allocate(ncv_ggff);
		press.allocate(ncv_ggff);
		temp.allocate(ncv_ggff);
		enthalpy.allocate(ncv_ggff);

		gamma.allocate(ncv_ggff);
		RoM.allocate(ncv_ggff);
		sos.allocate(ncv_ggff);

		strMag.allocate(ncv_gg);
		vortMag.allocate(ncv_gg);
		diverg.allocate(ncv_gg);

		//TBD_AD MOVE TO KOMSST
		blendFuncF1_AD.allocate(ncv_gg);
		blendFuncF2_AD.allocate(ncv_gg);
		crossDiff_AD.allocate(ncv_gg);
		//TBD_AD MOVE TO FPVA
		CmeanSource_AD.allocate(ncv_gg);
		muLam_AD.allocate(ncv_gg);
		LambdaOverCp_AD.allocate(ncv_gg);
		chi_AD.allocate(ncv_gg);

//		kine.clear();
		kine = NULL;

		// ----------------------------------------------------------------------------------------
		// init memory for face-based data
		// ----------------------------------------------------------------------------------------
		mul_fa.allocate(nfa_b2gg);
		lamOcp_fa.allocate(nfa_b2gg);
		mut_fa.allocate(nfa_b2gg);

		// ----------------------------------------------------------------------------------------
		// gradients
		// ----------------------------------------------------------------------------------------
		if (sndOrder == true) {
			grad_rho.allocate(ncv_gg);
#ifdef temp_reconstruction
			grad_temp.allocate(ncv_gg);
#else
			grad_p.allocate(ncv_gg);
#endif
		}

		grad_u.allocate(ncv_gg);         // allocate always

		if ((sndOrder == true) || (mu_ref > 0.0))
			grad_enthalpy.allocate(ncv_gg);

		int nScal = scalarTranspEqVector.size();

		if(nScal>0)
			scalarTranspEqVector_AD = new ScalarTranspEq_AD[nScal];
		else
			scalarTranspEqVector_AD = NULL;

		for(int i=0; i<nScal ; i++){
			scalarTranspEqVector_AD[i].phi.allocate(ncv_ggff);
			scalarTranspEqVector_AD[i].diff.allocate(nfa_b2gg);
			scalarTranspEqVector_AD[i].rhophi.allocate(ncv_ggff);
			scalarTranspEqVector_AD[i].grad_phi.allocate(ncv_gg);
			scalarTranspEqVector_AD[i].grad_rhophi.allocate(ncv_gg);
//			scalarTranspEqVector_AD[i].dpress_dphi.clear();
			strcpy(scalarTranspEqVector_AD[i].name, scalarTranspEqVector[i].getName());
		}

		set_adouble_names();
	}

	virtual void destroy_adoubles() {
		if(scalarTranspEqVector_AD != NULL)
			delete [] scalarTranspEqVector_AD;

		if ((sndOrder == true) || (mu_ref > 0.0))
			grad_enthalpy.clear();

		grad_u.clear();

		if (sndOrder == true) {
#ifdef temp_reconstruction
			grad_temp.clear();
#else
			grad_p.clear();
#endif
			grad_rho.clear();
		}

		mut_fa.clear();
		lamOcp_fa.clear();
		mul_fa.clear();

		chi_AD.clear();
		LambdaOverCp_AD.clear();
		muLam_AD.clear();
		CmeanSource_AD.clear();

		crossDiff_AD.clear();
		blendFuncF2_AD.clear();
		blendFuncF1_AD.clear();

		diverg.clear();
		vortMag.clear();
		strMag.clear();

		sos.clear();
		RoM.clear();
		gamma.clear();

		enthalpy.clear();
		temp.clear();
		press.clear();
		vel.clear();

		rhoE_AD.clear();
		rhou_AD.clear();
		rho_AD.clear();
	}

	void set_adouble_names() {
		rho_AD.setName("rho_AD");       ///< density derivative
		rhou_AD.setName("rhou_AD");  ///< momentum derivative
		rhoE_AD.setName("rhoE_AD");      ///< energy derivative

		vel.setName("vel");   ///< velocity
		press.setName("press");      ///< pressure e
		temp.setName("temp");       ///< temperature
		enthalpy.setName("enthalpy");   ///< chemical + sensible enthalpy (no kinetic energy!) derivative

		gamma.setName("gamma");      ///< gamma derivative
		RoM.setName("RoM");        ///< Gas constant derivative
		sos.setName("sos");        ///< speed of sound derivative

		mul_fa.setName("mul_fa");           ///< laminar viscosity at cell faces
		lamOcp_fa.setName("lamOcp_fa");        ///< heat conductivity at cell faces

		mut_fa.setName("mut_fa");     ///< turbulent viscosity derivative
		if(kine != NULL)
			(*kine).setName("kine");      ///< turbulent kinetic energy derivative

		strMag.setName("strMag");
		vortMag.setName("vortMag");
		diverg.setName("diverg");  ///< magnitude strain rate, divergence derivatives
		blendFuncF1_AD.setName("blendFuncF1_AD");
		blendFuncF2_AD.setName("blendFuncF2_AD");
		crossDiff_AD.setName("crossDiff_AD"); ///<SST model functions

		CmeanSource_AD.setName("CmeanSource_AD");
		muLam_AD.setName("muLam_AD");
		LambdaOverCp_AD.setName("LambdaOverCp_AD");
		chi_AD.setName("chi_AD"); ///<FPVA model functions

#ifdef USE_ARTIF_VISC_WITH_MEM_SAVING
		artifViscMag_AD.setName("artifViscMag_AD");
#endif

		grad_rho.setName("grad_rho");      ///< density  gradient der
		grad_u.setName("grad_u");     ///< velocity gradient tensor der
		grad_p.setName("grad_p");        ///< pressure gradient der
		grad_temp.setName("grad_temp");     ///< temperature gradient der
		grad_enthalpy.setName("grad_enthalpy"); ///< enthalpy gradient der

		int nScal = scalarTranspEqVector.size();
		for(int i=0; i<nScal ; i++){
			scalarTranspEqVector_AD[i].phi.setName("phi");
			scalarTranspEqVector_AD[i].diff.setName("diff");
			scalarTranspEqVector_AD[i].grad_phi.setName("grad_phi");
			scalarTranspEqVector_AD[i].rhophi.setName("rhophi");
			scalarTranspEqVector_AD[i].grad_rhophi.setName("grad_rhophi");
			scalarTranspEqVector_AD[i].dpress_dphi.setName("dpress_dphi");
		}
	}
#endif

	// =========================
	//  ORIGINAL FUNCTIONS
	// =========================
#ifdef USE_MEM_SAVING_ADVAR
	/*
	 * Method: calcCv2Grad_AD
	 * ----------------------
	 *
	 */
	void calcCv2Grad_AD(ADvector<adouble> &gradPhi, ADscalar<adouble> &phi, const int limiter, ADscalar<adouble>& refValue, const double epsilon) {
		calcCv2ScalGrad_AD(gradPhi, phi);

		switch (limiter) {
		case BARTH_JESPERSEN_MOD:
			limitCv2GradBJ_AD(gradPhi, phi, refValue, epsilon);
			break;

		default:
			break;
		}

		// Make sure we have gradients in ncv_g:ncv_gg-1

		double (*buf) = new double[ncv_ggff];
		FOR_I3 {
			FOR_ICV
				buf[icv] = gradPhi[icv][i].value();

			updateCvDataG1G2(buf, REPLACE_DATA);

			FOR_ICV_G2_ONLY
				gradPhi[icv][i] = buf[icv];
		}

		delete [] buf;
	}

	void calcCv2Grad_AD(ADtensor<adouble> &gradPhi, ADvector<adouble> &phi, const int limiter, ADscalar<adouble> &refValue, const double epsilon) {
		calcCv2VecGrad_AD(gradPhi, phi);

		switch (limiter) {
		case BARTH_JESPERSEN_MOD:
			limitCv2GradBJ_AD(gradPhi, phi, refValue, epsilon);
			break;

		default:
			break;
		}

		// Make sure we have gradients in ncv_g:ncv_gg-1

		double (*buf) = new double[ncv_ggff];
		FOR_I3 {
			FOR_J3 {
				FOR_ICV
					buf[icv] = gradPhi[icv][i][j].value();

				updateCvDataG1G2(buf, REPLACE_DATA);

				FOR_ICV_G2_ONLY
					gradPhi[icv][i][j] = buf[icv];
			}
		}

		delete [] buf;
	}


	/*
	 * Method: calcCv2ScalGrad_AD
	 * --------------------------
	 * populates the gradient in cv range [0:ncv_g-1]. Note that this requires
	 * phi to be defined over all cv's [0:ncv_ggff-1]...
	 */
	void calcCv2ScalGrad_AD(ADvector<adouble> &dphidxi, ADscalar<adouble> &phi) {
		// populates the gradient in cv range [0:ncv_g-1]. Note that this requires
		// phi to be defined over all cv's [0:ncv_ggff-1]...

		FOR_ICV_G {
			int noc_f = nbocv_all_i[icv];
			FOR_I3 dphidxi[icv][i] = nbocv_all_grad_coeff[noc_f][i]*phi[icv];

			int noc_l = nbocv_all_i[icv+1]-1;
			for (int noc = noc_f+1; noc <= noc_l; ++noc) {
				int icv_nbr = nbocv_all_v[noc];
				FOR_I3 dphidxi[icv][i] += nbocv_all_grad_coeff[noc][i]*phi[icv_nbr];
			}
		}
	}

	void calcCv2VecGrad_AD(ADtensor<adouble> &duidxj, ADvector<adouble> &u) {
		FOR_ICV_G {
			int noc_f = nbocv_all_i[icv];
			FOR_I3
				FOR_J3
					duidxj[icv][i][j] = nbocv_all_grad_coeff[noc_f][j]*u[icv][i];

			int noc_l = nbocv_all_i[icv+1]-1;
			for (int noc = noc_f+1; noc<=noc_l; ++noc)	{
				int icv_nbr = nbocv_all_v[noc];
				FOR_I3
					FOR_J3
						duidxj[icv][i][j] += nbocv_all_grad_coeff[noc][j]*u[icv_nbr][i];
			}
		}
	}

	/*
	 * Method: limitCv2GradBJ_AD
	 * -------------------------
	 * apply smooth barth jespersen limiters
	 */
	void limitCv2GradBJ_AD(ADvector<adouble> &grad_p, ADscalar<adouble> &phi, ADscalar<adouble> &refValue, const double epsilonSDWLS) {
		for (int icv = 0; icv < ncv_g; icv++) {
			adouble phiMax = phi[icv];
			adouble phiMin = phiMax;

			int noc_f = nbocv_all_i[icv];
			int noc_l = nbocv_all_i[icv + 1] - 1;
			// skip the diagonal in this loop...
			for (int noc = noc_f + 1; noc <= noc_l; noc++) {
				int icv_nbr = nbocv_all_v[noc];
				phiMax = max(phiMax, phi[icv_nbr]);
				phiMin = min(phiMin, phi[icv_nbr]);
			}

			adouble alfa = 1.0;

			if ((phiMax-phiMin) > 1.0e-12) {
				int foc_f = faocv_i[icv];
				int foc_l = faocv_i[icv + 1] - 1;
				for (int foc = foc_f; foc <= foc_l; foc++) {
					int ifa = faocv_v[foc];

					//if (ifa >= nfa_b)
					{
						double r0[3] = {0.0, 0.0, 0.0};
						vecMinVec3d(r0, x_fa[ifa], x_cv[icv]);
						adouble phifa = phi[icv] + vecDotVec3d_AD(r0, grad_p.getAt(icv));
#ifdef BJ
						if ((phifa - phi[icv]) > 0.0)               alfa = min(alfa, (phiMax-phi[icv])/(phifa-phi[icv]));
						else if ((phifa - phi[icv]) < 0.0)          alfa = min(alfa, (phiMin-phi[icv])/(phifa-phi[icv]));
#else
						adouble dp;
						adouble eps = max(epsilonSDWLS*fabs(refValue[icv]), SDWLSClippingScalar);
						adouble eps2 = eps*eps;
						adouble dm = phifa - phi[icv];

						if   ((phifa - phi[icv]) >= 0.0)    dp = phiMax - phi[icv];
						else                                dp = phiMin - phi[icv];

						alfa = min(alfa, (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2));
#endif
					}
				}
			}
			else alfa = 0.0;

			for (int i=0; i<3; i++)
				grad_p[icv][i] *= alfa;
		}
	}

	void limitCv2GradBJ_AD(ADtensor<adouble> &grad_p, ADvector<adouble> &phi, ADscalar<adouble> &refValue, const double epsilonSDWLS) {
		for (int icv = 0; icv < ncv_g; icv++) {
			for (int i=0; i<3; i++) {
				adouble phiMax = phi[icv][i];
				adouble phiMin = phiMax;

				int noc_f = nbocv_all_i[icv];
				int noc_l = nbocv_all_i[icv + 1] - 1;
				// skip the diagonal in this loop...
				for (int noc = noc_f + 1; noc <= noc_l; noc++) {
					int icv_nbr = nbocv_all_v[noc];
					phiMax = max(phiMax, phi[icv_nbr][i]);
					phiMin = min(phiMin, phi[icv_nbr][i]);
				}

				adouble alfa = 1.0;

				if ((phiMax-phiMin) > 1.0e-8) {
					int foc_f = faocv_i[icv];
					int foc_l = faocv_i[icv + 1] - 1;
					for (int foc = foc_f; foc <= foc_l; foc++) {
						int ifa = faocv_v[foc];

						//if (ifa >= nfa_b)
						{
							double r0[3] = {0.0, 0.0, 0.0};
							vecMinVec3d(r0, x_fa[ifa], x_cv[icv]);
							adouble phifa = phi[icv][i] + vecDotVec3d_AD(r0, grad_p.getAt(icv,i));

#ifdef BJ
							if ((phifa - phi[icv][i]) > 0.0)               alfa = min(alfa, (phiMax-phi[icv][i])/(phifa-phi[icv][i]));
							else if ((phifa - phi[icv][i]) < 0.0)          alfa = min(alfa, (phiMin-phi[icv][i])/(phifa-phi[icv][i]));
#else
							adouble dp;
							adouble eps  = epsilonSDWLS*refValue[icv];
							adouble eps2 = eps*eps;
							adouble dm = phifa - phi[icv][i];

							if ((phifa - phi[icv][i]) >= 0.0)   dp = phiMax - phi[icv][i];
							else                                dp = phiMin - phi[icv][i];

							alfa = min(alfa, (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2));
#endif
						}
					}
				}
				else alfa = 0.0;

				for (int j=0; j<3; j++)
					grad_p[icv][i][j] *= alfa;
			}
		}
	}

	/*
	 * Method: calcStateVariables_AD
	 * -----------------------------
	 *
	 */
	virtual void calcStateVariables_AD(ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<REALQ> &rhoE) {
		for (int icv = 0; icv < ncv_gg; icv++) {
#ifdef USE_CONDASSIGN
			if (rho[icv].value() <= 0.0)
				cout << "Negative density at xcv: " << x_cv[icv][0] << ", " << x_cv[icv][1] << ", " << x_cv[icv][2] << endl;
#else
			if (rho[icv] <= 0.0)
				cout << "Negative density at xcv: " << x_cv[icv][0] << ", " << x_cv[icv][1] << ", " << x_cv[icv][2] << endl;
#endif

			for (int i=0; i<3; i++)
				vel[icv][i] = rhou[icv][i]/rho[icv];

			REALQS kinecv = 0.0;

			if (kine != NULL)
				kinecv = (*kine)[icv];

			REALQ pr = (gamma[icv]-1.0)*(rhoE[icv]
			           - 0.5*(rhou[icv][0]*rhou[icv][0]+rhou[icv][1]*rhou[icv][1]+rhou[icv][2]*rhou[icv][2])/rho[icv]
			           - rho[icv]*kinecv);
			if (pr <= 0.0) {
				cout << "Calc RANS negative pressure at xcv: " << x_cv[icv][0] << ", " << x_cv[icv][1] << ", " << x_cv[icv][2] << gamma[icv]<<endl;
			}
			else
				press[icv] = pr;

			temp[icv] = press[icv]/(rho[icv]*RoM[icv]);
			enthalpy[icv] = gamma[icv]*RoM[icv]/(gamma[icv]-1.0)*temp[icv];
			sos[icv] = sqrt(gamma[icv]*press[icv]/rho[icv]);
		}
	}

	/*
	 * Method: calcMaterialProperties_AD
	 * ---------------------------------
	 *
	 */
	virtual void calcMaterialProperties_AD(ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<REALQ> &rhoE) {
		if (mu_ref > 0.0) {
			if (viscMode == "SUTHERLAND") {
				// internal faces
				//for (int ifa = nfa_b; ifa < nfa; ifa++)
				for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) // Ensure we sweep through internal faces of ghosts too
				{
					if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg)) {
						int icv0 = cvofa[ifa][0];
						int icv1 = cvofa[ifa][1];

						REALX dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
						vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
						vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
						REALX w0 = sqrt(vecDotVec3d(dx0, dx0));
						REALX w1 = sqrt(vecDotVec3d(dx1, dx1));

						REALQ temperature  = (w1*temp[icv0] + w0*temp[icv1])/(w0+w1);
						mul_fa[ifa] = mu_ref*pow(temperature/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(temperature + SL_Sref);
						lamOcp_fa[ifa] = mul_fa[ifa] / Pr;
					}
				}
				// boundary faces computed next in setBC
			} else if (viscMode == "POWERLAW") {
				// internal faces
				//for (int ifa = nfa_b; ifa < nfa; ifa++)
				for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) // Ensure we sweep through internal faces of ghosts too
				{
					if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg)) {
						int icv0 = cvofa[ifa][0];
						int icv1 = cvofa[ifa][1];

						REALX dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
						vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
						vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
						REALX w0 = sqrt(vecDotVec3d(dx0, dx0));
						REALX w1 = sqrt(vecDotVec3d(dx1, dx1));

						REALQ temperature  = (w1*temp[icv0] + w0*temp[icv1])/(w0+w1);
						mul_fa[ifa] = mu_ref*pow(temperature/T_ref, mu_power_law);
						lamOcp_fa[ifa] = mul_fa[ifa] / Pr;
					}
				}
				// boundary faces computed next in setBC
			} else {
				cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
				throw(-1);
			}
		}
	}

	virtual void calcRansTurbViscMuet_AD(ADscalar<adouble> &rho, ADvector<adouble> &rhou) {
		static int flag = 1;

		// provide zero mut for laminar calculations
		if (flag == 1) {
			flag = 0;
			for (int ifa=0; ifa<nfa; ifa++)
				mut_fa[ifa] = 0.0;
		}
	}
#endif

	// =========================
	//  1D FUNCTIONS
	// =========================
	/*
	 * Method: calcCv2Grad1D_AD
	 * ------------------------
	 * Gradient routine
	 * Original code: calcCv2Grad_AD() in UgpWithCvCompFlowAD.h
	 */
	void calcCv2Grad1D_AD(const int icvCenter, adouble (*gradPhi)[3], adouble *phi, const int limiter, adouble *refValue, const double epsilon) {
		calcCv2ScalGrad1D_AD(icvCenter, gradPhi, phi);

		switch (limiter) {
			case BARTH_JESPERSEN_MOD:
				limitCv2GradBJ1D_AD(icvCenter, gradPhi, phi, refValue, epsilon);
				break;

			default:
				break;
		}
//TEMP__


	//	// Make sure we have gradients in ncv_g:ncv_gg-1
	//	double (*buf) = new double[ncv_ggff];
	//	FOR_I3 {
	//		FOR_ICV
	//		buf[icv] = gradPhi[icv][i].value();
	//
	//		updateCvDataG1G2(buf, REPLACE_DATA);
	//
	//		FOR_ICV_G2_ONLY
	//		gradPhi[icv][i] = buf[icv];
	//	}
	//
	//	delete [] buf;
	}

#ifdef USE_MEM_SAVING_ADVAR
	void calcCv2Grad1D_AD(const int icvCenter, ADvector<adouble>& gradPhi, ADscalar<adouble>& phi, const int limiter, ADscalar<adouble>& refValue, const double epsilon) {
		calcCv2ScalGrad1D_AD(icvCenter, gradPhi, phi);

		switch (limiter) {
			case BARTH_JESPERSEN_MOD:
				limitCv2GradBJ1D_AD(icvCenter, gradPhi, phi, refValue, epsilon);
				break;

			default:
				break;
		}
	}

	void calcCv2Grad1D_AD(const int icvCenter, ADvector<adouble>& gradPhi, adouble *phi, const int limiter, adouble *refValue, const double epsilon) {
		if(mpi_rank==0)
			cout<<"POSSIBLE ERROR in IkeUgpWithCvCompFlow.h: calcCv2Grad1D_AD(int, ADvector, adouble, ...): You shouldn't call this function! Do nothing!"<<endl;
	}
#endif

	/*
	 * Method: calcCv2Grad1D_AD
	 * ------------------------
	 * Gradient routine
	 * Original code: calcCv2Grad_AD() in UgpWithCvCompFlowAD.h
	 */
	void calcCv2Grad1D_AD(const int icvCenter, adouble (*gradPhi)[3][3], adouble (*phi)[3], const int limiter, adouble *refValue, const double epsilon) {
		calcCv2VecGrad1D_AD(icvCenter, gradPhi, phi);

		switch (limiter) {
			case BARTH_JESPERSEN_MOD:
				limitCv2GradBJ1D_AD(icvCenter, gradPhi, phi, refValue, epsilon);
				break;

			default:
				break;
		}

	//	// Make sure we have gradients in ncv_g:ncv_gg-1
	//	double (*buf) = new double[ncv_ggff];
	//	FOR_I3 {
	//		FOR_J3 {
	//			FOR_ICV
	//			buf[icv] = gradPhi[icv][i][j].value();
	//
	//			updateCvDataG1G2(buf, REPLACE_DATA);
	//
	//			FOR_ICV_G2_ONLY
	//			gradPhi[icv][i][j] = buf[icv];
	//		}
	//	}
	//
	//	delete [] buf;
	}

#ifdef USE_MEM_SAVING_ADVAR
	void calcCv2Grad1D_AD(const int icvCenter, ADtensor<adouble>& gradPhi, ADvector<adouble>& phi, const int limiter, ADscalar<adouble>& refValue, const double epsilon) {
		calcCv2VecGrad1D_AD(icvCenter, gradPhi, phi);

		switch (limiter) {
			case BARTH_JESPERSEN_MOD:
				limitCv2GradBJ1D_AD(icvCenter, gradPhi, phi, refValue, epsilon);
				break;

			default:
				break;
		}
	}
	void calcCv2Grad1D_AD(const int icvCenter, ADtensor<adouble>& gradPhi, adouble (*phi)[3], const int limiter, adouble *refValue, const double epsilon) {
		if(mpi_rank==0)
			cout<<"POSSIBLE ERROR in IkeUgpWithCvCompFlow.h: calcCv2Grad1D_AD(int, ADtensor, adouble, ...): You shouldn't call this function! Do nothing!"<<endl;
	}
#endif

	/*
	 * Method: calcCv2ScalGrad1D_AD
	 * ----------------------------
	 * Calculate gradients for the 1-layer neighboring CVs of a center CV
	 * Original code: calcCv2ScalGrad_AD() in UgpWithCvCompFlowAD.h
	 */
	void calcCv2ScalGrad1D_AD(const int icvCenter, adouble (*dphidxi)[3], adouble * phi) {
		int noc_f_center = nbocv_all_i[icvCenter];
		int noc_l_center = nbocv_all_i[icvCenter+1]-1;

		for (int noc_center = noc_f_center; noc_center <= noc_l_center; ++noc_center) {
			int icv = nbocv_all_v[noc_center]; // Neighboring CVs

			// Calculate gradients
			if(icv<ncv_g) {
				int noc_f = nbocv_all_i[icv];
				FOR_I3 dphidxi[icv][i] = nbocv_all_grad_coeff[noc_f][i]*phi[icv];

				int noc_l = nbocv_all_i[icv+1]-1;
				for (int noc = noc_f+1; noc <= noc_l; ++noc) {
					int icv_nbr = nbocv_all_v[noc];
					FOR_I3 dphidxi[icv][i] += nbocv_all_grad_coeff[noc][i]*phi[icv_nbr];

				}
			}
		}
	}

#ifdef USE_MEM_SAVING_ADVAR
	void calcCv2ScalGrad1D_AD(const int icvCenter, ADvector<adouble>& dphidxi, ADscalar<adouble>& phi) {
		int noc_f_center = nbocv_all_i[icvCenter];
		int noc_l_center = nbocv_all_i[icvCenter+1]-1;

		for (int noc_center = noc_f_center; noc_center <= noc_l_center; ++noc_center) {
			int icv = nbocv_all_v[noc_center]; // Neighboring CVs

			// Calculate gradients
			if(icv<ncv_g) {
				int noc_f = nbocv_all_i[icv];
				FOR_I3
					dphidxi[icv][i] = nbocv_all_grad_coeff[noc_f][i]*phi[icv];

				int noc_l = nbocv_all_i[icv+1]-1;
				for (int noc = noc_f+1; noc <= noc_l; ++noc) {
					int icv_nbr = nbocv_all_v[noc];
					FOR_I3
						dphidxi[icv][i] += nbocv_all_grad_coeff[noc][i]*phi[icv_nbr];
				}
			}
		}
	}
#endif

	/*
	 * Method: calcCv2VecGrad1D_AD
	 * ---------------------------
	 * Calculate gradients for the 1-layer neighboring CVs of a center CV
	 * Original code: calcCv2VecGrad_AD() in UgpWithCvCompFlowAD.h
	 */
	void calcCv2VecGrad1D_AD(const int icvCenter, adouble(*duidxj)[3][3], adouble(*u)[3]) {
		int noc_f_center = nbocv_all_i[icvCenter];
		int noc_l_center = nbocv_all_i[icvCenter+1]-1;
		for (int noc_center = noc_f_center; noc_center <= noc_l_center; ++noc_center) {
			int icv = nbocv_all_v[noc_center]; // Neighboring CVs

			// Calculate gradients
			if(icv<ncv_g) {
				int noc_f = nbocv_all_i[icv];

				FOR_I3
					FOR_J3
						duidxj[icv][i][j] = nbocv_all_grad_coeff[noc_f][j]*u[icv][i];

				int noc_l = nbocv_all_i[icv+1]-1;
				for (int noc = noc_f+1; noc<=noc_l; ++noc) {
					int icv_nbr = nbocv_all_v[noc];
					FOR_I3
						FOR_J3
							duidxj[icv][i][j] += nbocv_all_grad_coeff[noc][j]*u[icv_nbr][i];
				}
			}
		}
	}

#ifdef USE_MEM_SAVING_ADVAR
	void calcCv2VecGrad1D_AD(const int icvCenter, ADtensor<adouble>& duidxj, ADvector<adouble>& u) {
		int noc_f_center = nbocv_all_i[icvCenter];
		int noc_l_center = nbocv_all_i[icvCenter+1]-1;
		for (int noc_center = noc_f_center; noc_center <= noc_l_center; ++noc_center) {
			int icv = nbocv_all_v[noc_center]; // Neighboring CVs

			// Calculate gradients
			if(icv<ncv_g) {
				int noc_f = nbocv_all_i[icv];

				FOR_I3
					FOR_J3
						duidxj[icv][i][j] = nbocv_all_grad_coeff[noc_f][j]*u[icv][i];

				int noc_l = nbocv_all_i[icv+1]-1;
				for (int noc = noc_f+1; noc<=noc_l; ++noc) {
					int icv_nbr = nbocv_all_v[noc];
					FOR_I3
						FOR_J3
							duidxj[icv][i][j] += nbocv_all_grad_coeff[noc][j]*u[icv_nbr][i];
				}
			}
		}
	}
#endif

	/*
	 * Method: limitCv2GradBJ1D_AD
	 * ---------------------------
	 * Apply smooth Barth-Jespersen limiters
	 * original code: limitCv2GradBJ_AD() in UgpWithCvCompFlowAD.h
	 */
	void limitCv2GradBJ1D_AD(const int icvCenter, adouble (*grad_p)[3], adouble *phi, adouble *refValue, const double epsilonSDWLS) {
		int noc_f_center = nbocv_all_i[icvCenter];
		int noc_l_center = nbocv_all_i[icvCenter+1]-1;
		for (int noc_center = noc_f_center; noc_center <= noc_l_center; ++noc_center) {
			int icv = nbocv_all_v[noc_center];
			if(icv<ncv_g) {
				adouble phiMax = phi[icv];
				adouble phiMin = phiMax;

				int noc_f = nbocv_all_i[icv];
				int noc_l = nbocv_all_i[icv + 1] - 1;
				// skip the diagonal in this loop...
				for (int noc = noc_f + 1; noc <= noc_l; noc++) {
					int icv_nbr = nbocv_all_v[noc];
					phiMax = max(phiMax, phi[icv_nbr]);
					phiMin = min(phiMin, phi[icv_nbr]);
				}

				adouble alfa = 1.0;

				if ((phiMax-phiMin) > 1.0e-12) {
					int foc_f = faocv_i[icv];
					int foc_l = faocv_i[icv + 1] - 1;
					for (int foc = foc_f; foc <= foc_l; foc++) {
						int ifa = faocv_v[foc];

						//if (ifa >= nfa_b)
						{
							double r0[3] = {0.0, 0.0, 0.0};
							vecMinVec3d(r0, x_fa[ifa], x_cv[icv]);
							adouble phifa = phi[icv] + vecDotVec3d_AD(r0, grad_p[icv]);
		#ifdef BJ
							if ((phifa - phi[icv]) > 0.0)               alfa = min(alfa, (phiMax-phi[icv])/(phifa-phi[icv]));
							else if ((phifa - phi[icv]) < 0.0)          alfa = min(alfa, (phiMin-phi[icv])/(phifa-phi[icv]));
		#else
							adouble dp;
							adouble eps = max(epsilonSDWLS*fabs(refValue[icv]), SDWLSClippingScalar);
							adouble eps2 = eps*eps;
							adouble dm = phifa - phi[icv];

							if   ((phifa - phi[icv]) >= 0.0)    dp = phiMax - phi[icv];
							else                                dp = phiMin - phi[icv];

							alfa = min(alfa, (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2));
		#endif
						}
					}
				}
				else alfa = 0.0;

				for (int i=0; i<3; i++)
					grad_p[icv][i] *= alfa;
			}
		}
	}
#ifdef USE_MEM_SAVING_ADVAR
	void limitCv2GradBJ1D_AD(const int icvCenter, ADvector<adouble>& grad_p, ADscalar<adouble>& phi, ADscalar<adouble>& refValue, const double epsilonSDWLS) {
		int noc_f_center = nbocv_all_i[icvCenter];
		int noc_l_center = nbocv_all_i[icvCenter+1]-1;
		for (int noc_center = noc_f_center; noc_center <= noc_l_center; ++noc_center) {
			int icv = nbocv_all_v[noc_center];
			if(icv<ncv_g) {
				adouble phiMax = phi[icv];
				adouble phiMin = phiMax;

				int noc_f = nbocv_all_i[icv];
				int noc_l = nbocv_all_i[icv + 1] - 1;
				// skip the diagonal in this loop...
				for (int noc = noc_f + 1; noc <= noc_l; noc++) {
					int icv_nbr = nbocv_all_v[noc];
					phiMax = max(phiMax, phi[icv_nbr]);
					phiMin = min(phiMin, phi[icv_nbr]);
				}

#ifdef USE_CONDASSIGN // Note: condassign(a,b,c,d) <==> a = (b>0) ? c : d;
				adouble alfa = 1.0;
				if ((phiMax-phiMin) > 1.0e-12) {
					int foc_f = faocv_i[icv];
					int foc_l = faocv_i[icv + 1] - 1;
					for (int foc = foc_f; foc <= foc_l; foc++) {
						int ifa = faocv_v[foc];

						//if (ifa >= nfa_b)
						{
							double r0[3] = {0.0, 0.0, 0.0};
							vecMinVec3d(r0, x_fa[ifa], x_cv[icv]);
							adouble phifa = phi[icv] + vecDotVec3d_AD(r0, grad_p[icv]);
		#ifdef BJ
							if ((phifa - phi[icv]) > 0.0)               alfa = min(alfa, (phiMax-phi[icv])/(phifa-phi[icv]));
							else if ((phifa - phi[icv]) < 0.0)          alfa = min(alfa, (phiMin-phi[icv])/(phifa-phi[icv]));
		#else
							adouble dp;
							adouble eps = max(epsilonSDWLS*fabs(refValue[icv]), SDWLSClippingScalar);
							adouble eps2 = eps*eps;
							adouble dm = phifa - phi[icv];

							if   ((phifa - phi[icv]) >= 0.0)    dp = phiMax - phi[icv];
							else                                dp = phiMin - phi[icv];

							alfa = min(alfa, (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2));
		#endif
						}
					}
				}
				else alfa = 0.0;
//CASE1
//				adouble activeLimVal = 1.0;
//
//				int foc_f = faocv_i[icv];
//				int foc_l = faocv_i[icv + 1] - 1;
//				for (int foc = foc_f; foc <= foc_l; foc++) {
//					int ifa = faocv_v[foc];
//
//					double r0[3] = {0.0, 0.0, 0.0};
//					vecMinVec3d(r0, x_fa[ifa], x_cv[icv]);
//					adouble phifa = phi[icv] + vecDotVec3d_AD(r0, grad_p[icv]);
//		#ifdef BJ
//					if ((phifa - phi[icv]) > 0.0)               alfa = min(alfa, (phiMax-phi[icv])/(phifa-phi[icv]));
//					else if ((phifa - phi[icv]) < 0.0)          alfa = min(alfa, (phiMin-phi[icv])/(phifa-phi[icv]));
//		#else
//					adouble dp;
//					adouble eps = fmax(epsilonSDWLS*fabs(refValue[icv]), SDWLSClippingScalar);
//					adouble eps2 = eps*eps;
//					adouble dm = phifa - phi[icv];
//
//					condassign(dp, phifa - phi[icv] + MACHINE_EPS, phiMax - phi[icv], phiMin - phi[icv]);
//
//					activeLimVal = fmin(activeLimVal, (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2));
//		#endif
//				}
//
//				adouble alfa;
//				condassign(alfa, phiMax-phiMin-1.0e-12, activeLimVal, 0.0);
#else
				adouble alfa = 1.0;
				if ((phiMax-phiMin) > 1.0e-12) {
					int foc_f = faocv_i[icv];
					int foc_l = faocv_i[icv + 1] - 1;
					for (int foc = foc_f; foc <= foc_l; foc++) {
						int ifa = faocv_v[foc];

						//if (ifa >= nfa_b)
						{
							double r0[3] = {0.0, 0.0, 0.0};
							vecMinVec3d(r0, x_fa[ifa], x_cv[icv]);
							adouble phifa = phi[icv] + vecDotVec3d_AD(r0, grad_p[icv]);
		#ifdef BJ
							if ((phifa - phi[icv]) > 0.0)               alfa = min(alfa, (phiMax-phi[icv])/(phifa-phi[icv]));
							else if ((phifa - phi[icv]) < 0.0)          alfa = min(alfa, (phiMin-phi[icv])/(phifa-phi[icv]));
		#else
							adouble dp;
							adouble eps = max(epsilonSDWLS*fabs(refValue[icv]), SDWLSClippingScalar);
							adouble eps2 = eps*eps;
							adouble dm = phifa - phi[icv];

							if   ((phifa - phi[icv]) >= 0.0)    dp = phiMax - phi[icv];
							else                                dp = phiMin - phi[icv];

							alfa = min(alfa, (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2));
		#endif
						}
					}
				}
				else alfa = 0.0;
#endif

				for (int i=0; i<3; i++)
					grad_p[icv][i] *= alfa;
			}
		}
	}
#endif

	/*
	 * Method: limitCv2GradBJ1D_AD
	 * ---------------------------
	 * Apply smooth Barth-Jespersen limiters
	 * original code: limitCv2GradBJ_AD() in UgpWithCvCompFlowAD.h
	 */
	void limitCv2GradBJ1D_AD(const int icvCenter, adouble (*grad_p)[3][3], adouble (*phi)[3], adouble *refValue, const double epsilonSDWLS)
	{
		int noc_f_center = nbocv_all_i[icvCenter];
		int noc_l_center = nbocv_all_i[icvCenter+1]-1;
		for (int noc_center = noc_f_center; noc_center <= noc_l_center; ++noc_center) {
			int icv = nbocv_all_v[noc_center];
			if(icv<ncv_g) {
				for (int i=0; i<3; i++) {
					adouble phiMax = phi[icv][i];
					adouble phiMin = phiMax;

					int noc_f = nbocv_all_i[icv];
					int noc_l = nbocv_all_i[icv + 1] - 1;
					// skip the diagonal in this loop...
					for (int noc = noc_f + 1; noc <= noc_l; noc++) {
						int icv_nbr = nbocv_all_v[noc];
						phiMax = max(phiMax, phi[icv_nbr][i]);
						phiMin = min(phiMin, phi[icv_nbr][i]);
					}

					adouble alfa = 1.0;

					if ((phiMax-phiMin) > 1.0e-8) {
						int foc_f = faocv_i[icv];
						int foc_l = faocv_i[icv + 1] - 1;
						for (int foc = foc_f; foc <= foc_l; foc++) {
							int ifa = faocv_v[foc];

							//if (ifa >= nfa_b)
							{
								double r0[3] = {0.0, 0.0, 0.0};
								vecMinVec3d(r0, x_fa[ifa], x_cv[icv]);
								adouble phifa = phi[icv][i] + vecDotVec3d_AD(r0, grad_p[icv][i]);

		#ifdef BJ
								if ((phifa - phi[icv][i]) > 0.0)               alfa = min(alfa, (phiMax-phi[icv][i])/(phifa-phi[icv][i]));
								else if ((phifa - phi[icv][i]) < 0.0)          alfa = min(alfa, (phiMin-phi[icv][i])/(phifa-phi[icv][i]));
		#else
								adouble dp;
								adouble eps  = epsilonSDWLS*refValue[icv];
								adouble eps2 = eps*eps;
								adouble dm = phifa - phi[icv][i];

								if ((phifa - phi[icv][i]) >= 0.0)   dp = phiMax - phi[icv][i];
								else                                dp = phiMin - phi[icv][i];

								alfa = min(alfa, (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2));
		#endif
							}
						}
					}
					else alfa = 0.0;

					for (int j=0; j<3; j++)
						grad_p[icv][i][j] *= alfa;
				}
			}
		}
	}
#ifdef USE_MEM_SAVING_ADVAR
	void limitCv2GradBJ1D_AD(const int icvCenter, ADtensor<adouble>& grad_p, ADvector<adouble>& phi, ADscalar<adouble>& refValue, const double epsilonSDWLS) {
		int noc_f_center = nbocv_all_i[icvCenter];
		int noc_l_center = nbocv_all_i[icvCenter+1]-1;
		for (int noc_center = noc_f_center; noc_center <= noc_l_center; ++noc_center) {
			int icv = nbocv_all_v[noc_center];
			if(icv<ncv_g) {
				for (int i=0; i<3; i++) {
					adouble phiMax = phi[icv][i];
					adouble phiMin = phiMax;

					int noc_f = nbocv_all_i[icv];
					int noc_l = nbocv_all_i[icv + 1] - 1;
					// skip the diagonal in this loop...
					for (int noc = noc_f + 1; noc <= noc_l; noc++) {
						int icv_nbr = nbocv_all_v[noc];
#ifdef USE_CONDASSIGN
						phiMax = fmax(phiMax, phi[icv_nbr][i]);
						phiMin = fmin(phiMin, phi[icv_nbr][i]);
#else
						phiMax = max(phiMax, phi[icv_nbr][i]);
						phiMin = min(phiMin, phi[icv_nbr][i]);
#endif
					}

#ifdef USE_CONDASSIGN
					adouble alfa = 1.0;
					if ((phiMax-phiMin) > 1.0e-8) {
						int foc_f = faocv_i[icv];
						int foc_l = faocv_i[icv + 1] - 1;
						for (int foc = foc_f; foc <= foc_l; foc++) {
							int ifa = faocv_v[foc];

							//if (ifa >= nfa_b)
							{
								double r0[3] = {0.0, 0.0, 0.0};
								vecMinVec3d(r0, x_fa[ifa], x_cv[icv]);
								adouble phifa = phi[icv][i] + vecDotVec3d_AD(r0, grad_p[icv][i]);
		#ifdef BJ
								if ((phifa - phi[icv][i]) > 0.0)               alfa = min(alfa, (phiMax-phi[icv][i])/(phifa-phi[icv][i]));
								else if ((phifa - phi[icv][i]) < 0.0)          alfa = min(alfa, (phiMin-phi[icv][i])/(phifa-phi[icv][i]));
		#else
								adouble dp;
								adouble eps  = epsilonSDWLS*refValue[icv];
								adouble eps2 = eps*eps;
								adouble dm = phifa - phi[icv][i];

								if ((phifa - phi[icv][i]) >= 0.0)   dp = phiMax - phi[icv][i];
								else                                dp = phiMin - phi[icv][i];

								alfa = fmin(alfa, (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2));
		#endif
							}
						}
					} else {
						alfa = 0.0;
					}
//CASE1
//					adouble activeLimVal = 1.0;
//					int foc_f = faocv_i[icv];
//					int foc_l = faocv_i[icv + 1] - 1;
//					for (int foc = foc_f; foc <= foc_l; foc++) {
//						int ifa = faocv_v[foc];
//
//						double r0[3] = {0.0, 0.0, 0.0};
//						vecMinVec3d(r0, x_fa[ifa], x_cv[icv]);
//						adouble phifa = phi[icv][i] + vecDotVec3d_AD(r0, grad_p[icv][i]);
//	#ifdef BJ
//						if ((phifa - phi[icv][i]) > 0.0)               alfa = min(alfa, (phiMax-phi[icv][i])/(phifa-phi[icv][i]));
//						else if ((phifa - phi[icv][i]) < 0.0)          alfa = min(alfa, (phiMin-phi[icv][i])/(phifa-phi[icv][i]));
//	#else
//						adouble dp;
//						adouble eps  = epsilonSDWLS*refValue[icv];
//						adouble eps2 = eps*eps;
//						adouble dm = phifa - phi[icv][i];
//
//						condassign(dp, phifa - phi[icv][i] - MACHINE_EPS, phiMax - phi[icv][i], phiMin - phi[icv][i]); // Note: condassign(a,b,c,d) <==> a = (b>0) ? c : d;
//
//						activeLimVal = fmin(activeLimVal, (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2));
//	#endif
//					}
//					adouble alfa;
//					condassign(alfa, phiMax-phiMin-1.0e-8, activeLimVal, 0.0);
#else
					adouble alfa = 1.0;
					if ((phiMax-phiMin) > 1.0e-8) {
						int foc_f = faocv_i[icv];
						int foc_l = faocv_i[icv + 1] - 1;
						for (int foc = foc_f; foc <= foc_l; foc++) {
							int ifa = faocv_v[foc];

							//if (ifa >= nfa_b)
							{
								double r0[3] = {0.0, 0.0, 0.0};
								vecMinVec3d(r0, x_fa[ifa], x_cv[icv]);
								adouble phifa = phi[icv][i] + vecDotVec3d_AD(r0, grad_p[icv][i]);
		#ifdef BJ
								if ((phifa - phi[icv][i]) > 0.0)               alfa = min(alfa, (phiMax-phi[icv][i])/(phifa-phi[icv][i]));
								else if ((phifa - phi[icv][i]) < 0.0)          alfa = min(alfa, (phiMin-phi[icv][i])/(phifa-phi[icv][i]));
		#else
								adouble dp;
								adouble eps  = epsilonSDWLS*refValue[icv];
								adouble eps2 = eps*eps;
								adouble dm = phifa - phi[icv][i];

								if ((phifa - phi[icv][i]) >= 0.0)   dp = phiMax - phi[icv][i];
								else                                dp = phiMin - phi[icv][i];

#ifdef USE_CONDASSIGN
								alfa = fmin(alfa, (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2));
#else
								alfa = min(alfa, (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2));
#endif
		#endif
							}
						}
					} else {
						alfa = 0.0;
					}
#endif

					for (int j=0; j<3; j++)
						grad_p[icv][i][j] *= alfa;
				}
			}
		}
	}
#endif

    /*
     * Method: checkNegativeScalarCv_JOE
     * ---------------------------------
     * Check if negative omega occurs at the given icv.
     * Note that kine should be checked in the IkeWithPsAlc method
     *
     * Return: by value     = number of negative scalar found
     *                                                   oldOmega - newOmega
     *         by reference = newRelax = safeParameter * -------------------
     *                                                        deltaOmega
     *         Note: If no negative scalar found, newRelax is same as oldRelax
     */
    virtual int checkNegativeTurbScalarCv_JOE(double& newRelax, const double oldRelax, const int icv, const double clipParameter, const double safeParameter, const double* qArray, const double* delQ) {
        if(mpi_rank==0) cout<<"WARNING checkNegativeTurbScalarCv_JOE(): empty method"<<endl;
        newRelax = oldRelax;
        return 0;
    }

    /*
     * Method: checkNegativeTurbScalarFa_JOE
     * -------------------------------------
     * Return: number of negative values found at the face 
     */
    virtual int checkNegativeTurbScalarFa_JOE(const int ifa, const int icv0, const int icv1, const bool bothFaces) {
        if(mpi_rank==0) cout<<"WARNING checkNegativeTurbScalarFa_JOE(): empty method"<<endl;
        return 0;
    }

    // BARRIER FOR TURBULENT SCALARS
    /*
     * Method: readBarrierParamTurbScalars
     * -----------------------------------
     * Read parameters for turbulent scalar barrier functions
     * Note: barrierSourceTurbScalars(), barrierSourceTurbScalars_AD(), and barrierSourceTurbScalars1D_AD() are calling this method
     * Return: true if the barreir functions will be used
     */
    virtual bool readBarrierParamTurbScalars(string& barrierFunctionName, int& maxIter, double& threshold_resid, vector<double>& coeffScalars,
            const int iterNewton, const double resid_norm_tot /* = ABSURDLY_BIG_NUMBER = 2.22e22 */ ) {
        return false;
    }

    /*
     * Method: barrierSourceTurbScalars
     * ----------------------------
     * Add barrier functions to the RHS of the scalars equations
     */
    virtual void barrierSourceTurbScalars(double* rhs, const int nScal, const int iterNewton, const double residNormTotOld) {
    //	if(debugLevel>0 && mpi_rank==0) cout<<"barrierSourceTurbScalars() is currently empty!"<<endl;
    }

    /*
     * Method: barrierSourceTurbScalars_AD
     * -----------------------------------
     * Add barrier functions to the RHS of the scalars equations
     */
    virtual void barrierSourceTurbScalars_AD(REALAS **rhs_rhoScal_AD, const int nScal, const int iterNewton, const double residNormTotOld) {

    //	if(debugLevel>0 && mpi_rank==0) cout<<"barrierSourceTurbScalars_AD() is currently empty!"<<endl;
    }

    /*
     * Method: barrierSourceTurbScalars1D_AD
     * ---------------------------------
     * Add barrier functions to the RHS of the scalars equations
     */
    virtual void barrierSourceTurbScalars1D_AD(const int icvCenter, REALAS* rhs_rhoScal_AD, const int nScal, const int iterNewton, const double residNormTotOld) {

    //	if(debugLevel>0 && mpi_rank==0 && icvCenter==0) cout<<"barrierSourceTurbScalars1D_AD() is currently empty!"<<endl;
    }

	// =========================
	//  HOOK FUNCTIONS
	// =========================
#ifdef USE_MEM_SAVING_ADVAR
	/*
	 * Method: boundaryHookScalarRansTurb_AD
	 * -------------------------------------
	 *
	 */
	virtual void boundaryHookScalarRansTurb_AD(ADscalar<adouble>& phi_ph, FaZone *zone, const string &name)  {
		/* empty */
	}
	/*
	 * Method: sourceHookScalarRansTurb_new_AD
	 * ---------------------------------------
	 *
	 */
	virtual void sourceHookScalarRansTurb_new_AD(ADscalar<adouble>& rhs, double *A, const string &name, int flagImplicit) {
		/* empty */
	}

	/*
	 * Method: sourceHookScalarRansTurb_new_AD
	 * ---------------------------------------
	 *
	 */
	virtual void boundaryHookScalarRansComb_AD(ADscalar<adouble>& phi_ph, FaZone *zone, const string &name)  {
		/* empty */
	}
	/*
	 * Method: sourceHookScalarRansComb_new_AD
	 * ---------------------------------------
	 *
	 */
	virtual void sourceHookScalarRansComb_new_AD(ADscalar<adouble>& rhs, double *A, const string &name, int flagImplicit)  {
		/* empty */
	}
#endif

	/*
	 * Method: initialHookScalarRansCombModel1D_AD
	 * -------------------------------------------
	 * This method is doing the exact same thing as initialHookScalarRansCombModel_AD()
	 * but has extra consideration that makes this method run only at the first call
	 * Original code = initialHookScalarRansCombModel_AD in UgpWithCvCompFlowAD.h
	 */
	virtual void initialHookScalarRansCombModel1D_AD(int loadtable, bool &firstCall) {
		if(firstCall) {
			for (int icv = 0; icv < ncv_ggff; icv++) {
				RoM[icv] = R_gas;
				gamma[icv] = GAMMA;
			}
			firstCall = false;
		}
	}

	/*
	 * Method: initialHookScalarRansTurbModel1D_AD()
	 * ---------------------------------------------
	 * This method is doing the exact same thing as initialHookScalarRansCombModel_AD()
	 * but has extra consideration that makes this method run only at the first call
	 * Original code = initialHookScalarRansTurbModel_AD() in UgpWithCvCompFlowAD.h
	 */
	virtual void initialHookScalarRansTurbModel1D_AD(bool &firstCall) {
		if(firstCall)
			firstCall = false;
	}

#ifdef USE_MEM_SAVING_ADVAR
	/*
	 * Method: initialHookScalarRansCombModel1D_AD
	 * -------------------------------------------
	 * This method is doing the exact same thing as initialHookScalarRansCombModel_AD()
	 * but the implementation that works only for the ICVs stored in nbocv2ff
	 * Original code = initialHookScalarRansCombModel_AD in UgpWithCvCompFlowAD.h
	 */
	virtual void initialHookScalarRansCombModel1D_AD(int loadtable, vector<int> &nbocv2ff, bool &firstCall) {
		for(size_t i=0; i<nbocv2ff.size(); ++i) {
			int icv = nbocv2ff[i];
			RoM[icv] = R_gas;
			gamma[icv] = GAMMA;
		}

		if(firstCall)
			firstCall = false;
	}

	/*
	 * Method: initialHookScalarRansTurbModel1D_AD()
	 * ---------------------------------------------
	 * This method is doing the exact same thing as initialHookScalarRansCombModel_AD()
	 * but the implementation that works only for the ICVs stored in nbocv2ff
	 * Original code = initialHookScalarRansTurbModel_AD() in UgpWithCvCompFlowAD.h
	 */
	virtual void initialHookScalarRansTurbModel1D_AD(vector<int> &nbocv2ff, bool &firstCall) {
		/* empty */
	}

	/*
	 * Method: finalHookScalarRansTurbModel1D_AD()
	 * ---------------------------------------------
	 *
	 */
	virtual void finalHookScalarRansTurbModel1D_AD() {
		/* empty */
	}
#endif

	/*
	 * Method: boundaryHookScalarRansTurb1D_AD
	 * ---------------------------------------
	 * Original code = boundaryHookScalarRansTurb_AD() in UgpWithCvCompFlowAD.h
	 */
    virtual void boundaryHookScalarRansTurb1D_AD(const int ifa, adouble *phi_ph, FaZone *zone, const string &scalName)  {/*empty*/}

	/*
	 * Method: boundaryHookScalarRansComb1D_AD
	 * ---------------------------------------
	 * Original code = boundaryHookScalarRansComb_AD in UgpWithCvCompFlowAD.h
	 */
    virtual void boundaryHookScalarRansComb1D_AD(const int ifa, adouble *phi_ph, FaZone *zone, const string &scalName)  {/*empty*/}

	/*
	 * Method: boundaryHookScalarRansTurb1D
	 * ------------------------------------
	 * Original code = boundaryHookScalarRansTurb_AD() in IkeUgpWithCvCompFlow.h and boundaryHookScalarRansTurb() in UgpWithCvCompFlowAD.h
	 * In some situation, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
	 */
    virtual void boundaryHookScalarRansTurb1D(const int ifa, double *phi_ph, FaZone *zone, const string &scalName)  {/*empty*/}

	/*
	 * Method: boundaryHookScalarRansComb1D
	 * ------------------------------------
	 * Original code = boundaryHookScalarRansComb1D_AD in IkeUgpWithCvCompFlow.h and boundaryHookScalarRansComb() in UgpWithCvCompFlowAD.h
	 * In some situation, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
	 */
    virtual void boundaryHookScalarRansComb1D(const int ifa, double *phi_ph, FaZone *zone, const string &scalName)  {/*empty*/}

	/*
	 * Method: diffusivityHookScalarRansTurb1D_AD
	 * ------------------------------------------
	 * Original code: diffusivityHookScalarRansTurb_AD in UgpWithCvCompFlowAD.h
	 */
    virtual void diffusivityHookScalarRansTurb1D_AD(const int ifa, FaZone* zone, const string &name)  {/*empty*/}

	/*
	 * Method: diffusivityHookScalarRansComb1D_AD
	 * ------------------------------------------
	 * Original code: diffusivityHookScalarRansComb_AD in UgpWithCvCompFlowAD.h
	 */
    virtual void diffusivityHookScalarRansComb1D_AD(const int ifa, FaZone* zone, const string &name)  {/*empty*/}

	/*
	 * Method: sourceHookScalarRansTurb_new1D_AD
	 * ------------------------------------------
	 * Original code: sourceHookScalarRansTurb_new_AD in UgpWithCvCompFlowAD.h
	 */
    virtual void sourceHookScalarRansTurb_new1D_AD(const int icvCenter, adouble &rhs, double *A, const string &name, int flagImplicit) { /*empty*/ }

	/*
	 * Method: sourceHookScalarRansComb_new1D_AD
	 * ------------------------------------------
	 * Original code: sourceHookScalarRansComb_new_AD in UgpWithCvCompFlowAD.h
	 */
    virtual void sourceHookScalarRansComb_new1D_AD(const int icvCenter, adouble &rhs, double *A, const string &name, int flagImplicit)  {/*empty*/}

    /*
     * Method: sourceHookRansTurbCoupled1D_AD
     * --------------------------------------
     * Original code: sourceHookRansTurbCoupled_AD in UgpWithCvCompFlowAD.h
     */
    virtual void sourceHookRansTurbCoupled1D_AD(const int icvCenter, adouble *rhs, double ***A, int flagImplicit) {/*empty*/}

    /*
     * Method: sourceHookRansCombCoupled1D_AD
     * --------------------------------------
     * Original code: sourceHookRansCombCoupled_AD in UgpWithCvCompFlowAD.h
     */
    virtual void sourceHookRansCombCoupled1D_AD(const int icvCenter, adouble *rhs, double ***A, int flagImplicit) {/*empty*/}

#ifdef USE_MEM_SAVING_ADVAR
	/*
	 * Method: boundaryHookScalarRansTurb1D_AD
	 * ---------------------------------------
	 * Original code = boundaryHookScalarRansTurb_AD() in UgpWithCvCompFlowAD.h
	 */
    virtual void boundaryHookScalarRansTurb1D_AD(const int ifa, ADscalar<adouble>& phi_ph, FaZone *zone, const string &scalName)  {
		/* empty */
    }

	/*
	 * Method: boundaryHookScalarRansComb1D_AD
	 * ---------------------------------------
	 * Original code = boundaryHookScalarRansComb_AD in UgpWithCvCompFlowAD.h
	 */
    virtual void boundaryHookScalarRansComb1D_AD(const int ifa, ADscalar<adouble>& phi_ph, FaZone *zone, const string &scalName)  {
		/* empty */
    }

	/*
	 * Method: boundaryHookScalarRansTurb1D
	 * ------------------------------------
	 * Original code = boundaryHookScalarRansTurb_AD() in IkeUgpWithCvCompFlow.h and boundaryHookScalarRansTurb() in UgpWithCvCompFlowAD.h
	 * In some situation, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
	 */
    virtual void boundaryHookScalarRansTurb1D(const int ifa, ADscalar<adouble>& phi_ph, FaZone *zone, const string &scalName)  {
		/* empty */
    }

	/*
	 * Method: boundaryHookScalarRansComb1D
	 * ------------------------------------
	 * Original code = boundaryHookScalarRansComb1D_AD in IkeUgpWithCvCompFlow.h and boundaryHookScalarRansComb() in UgpWithCvCompFlowAD.h
	 * In some situation, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
	 */
    virtual void boundaryHookScalarRansComb1D(const int ifa, ADscalar<adouble>& phi_ph, FaZone *zone, const string &scalName)  {
		/* empty */
    }

    /*
     * Method: sourceHookRansTurbCoupled1D_AD
     * --------------------------------------
     * Original code: sourceHookRansTurbCoupled_AD in UgpWithCvCompFlowAD.h
     */
    virtual void sourceHookRansTurbCoupled1D_AD(const int icv, ADscalar<adouble>& rhs, double ***A, int flagImplicit) {
		/* empty */
    }
#endif

    // =========================
    //  OTHER VIRTUAL FUNCTIONS
    // =========================
    /*
     * Method: UgpWithCvCompFlow_AD_init
     * ---------------------------------
     * Initialize some UgpWithCvCompFlow_AD variables
     */
    virtual void UgpWithCvCompFlow_AD_init() {/*empty*/}

    /*
	 * Method: UgpWithCvCompFlow_AD_clear
	 * ---------------------------------
	 * Clear some UgpWithCvCompFlow_AD variables
	 */
    virtual void UgpWithCvCompFlow_AD_clear() {/*empty*/}

#ifdef USE_MEM_SAVING_ADVAR
    /*
     * Method: calcRansTurbViscMuet1D_AD
     * ---------------------------------
     * Calculate turbulence viscosity at some faces
     * By default, this virtual method initialize mut_fa as zero. Thus, it must be overloaded later in a turbulence-model class
     * Original code = calcRansTurbViscMuet_AD() in UgpWithCvCompFlowAD.h
     */
    virtual void calcRansTurbViscMuet1D_AD(const int icvCenter, FaZone** zoneArrayBoundary, vector<int>& faInternal, vector<int>& faBoundary, ADscalar<REALQ> &rho_AD, ADvector<REALQ> &rhou_AD) {
    	// provide zero mut for laminar calculations
    	for (size_t index=0; index<faInternal.size(); ++index) {
    		int ifa = faInternal[index];
    		mut_fa[ifa] = 0.0;
    	}
    	for (size_t index=0; index<faBoundary.size(); ++index) {
    		int ifa = faBoundary[index];
    		mut_fa[ifa] = 0.0;
    	}
    }
#endif
    /*
     * Method: calcRansTurbViscMuet1D_AD
     * ---------------------------------
     * Calculate turbulence viscosity at some faces
     * By default, this virtual method initialize mut_fa as zero. Thus, it must be overloaded later in a turbulence-model class
     * Original code = calcRansTurbViscMuet_AD() in UgpWithCvCompFlowAD.h
     */
    virtual void calcRansTurbViscMuet1D_AD(const int icvCenter, FaZone** zoneArrayBoundary, vector<int>& faInternal, vector<int>& faBoundary, REALQ *rho, REALQ (*rhou)[3]) {
    	// provide zero mut for laminar calculations
    	for (size_t index=0; index<faInternal.size(); ++index) {
    		int ifa = faInternal[index];
    		mut_fa[ifa] = 0.0;
    	}
    	for (size_t index=0; index<faBoundary.size(); ++index) {
    		int ifa = faBoundary[index];
    		mut_fa[ifa] = 0.0;
    	}
    }

    /*
     * Method: calcGradVel1D_AD
     * ------------------------
     * Original code = calcGradVel in UgpWithCvCompFlow.h
     * (Note: The current version of UgpWithCvCompFlowAD.h doesn't have this function)
     */
    virtual void calcGradVel1D_AD(const int icvCenter) {
    	calcCv2Grad1D_AD(icvCenter, grad_u, vel, limiterNavierS, sos, epsilonSDWLS);
    }

    /*
     * Method: InterpolateAtCellCenterFromFaceValues_AD
     * ------------------------------------------------
     * Original code = InterpolateAtCellCenterFromFaceValues() in UgpWithCvCompFlow.h
     *                 InterpolateAtCellCenterFromFaceValues() in UgpWithCvCompFlowAD.h
     * (Note: The current version of UgpWithCvCompFlowAD.h doesn't support ADscalar)
     */
    adouble InterpolateAtCellCenterFromFaceValues_AD(ADscalar<adouble> &phi, int icv) {
    	adouble phiC1 = 0.0;

    	int foc_f = faocv_i[icv];
    	int foc_l = faocv_i[icv + 1] - 1;
    	for (int foc=foc_f; foc<=foc_l; foc++) {
    		int ifa = faocv_v[foc];
    		phiC1 += phi[ifa];
    	}
    	adouble phiC = phiC1/(double)(foc_l-foc_f+1);

    	return phiC;
    }

    /*
     * Method: calcStrainRateAndDivergence1D_AD
     * ----------------------------------------
     * Original code = calcStrainRateAndDivergence_AD() in UgpWithCvCompFlowAD.h
     */
    virtual void calcStrainRateAndDivergence1D_AD(const int icvCenter) {
    	// Prevent calling this method twice for an icvCenter: Note that this algorithm works only for the 1D-style
    	static int myFlag = -1;
    	if(myFlag == icvCenter && ncv > 1)
    		return;  // If this method has been already called once for the current icvCenter, just return.
    	myFlag = icvCenter;


#ifdef USE_MEM_SAVING_ADVAR
    	if ((strMag.empty()) || (diverg.empty())) {
    		cout << "Error in calcStrainRateAndDivergence1D_AD() => Either strMag or diverg has not been registered!" << endl;
    		throw(-1);
    	}
#else
    	if (strMag==NULL || diverg==NULL) {
    		cout << "Error in calcStrainRateAndDivergence1D_AD() => Either strMag or diverg has not been registered!" << endl;
    		throw(-1);
    	}
#endif

    	int noc_f = nbocv_i[icvCenter];
    	int noc_l = nbocv_i[icvCenter+1] - 1;
    	for (int noc = noc_f; noc <= noc_l; noc++) {
    		int icv = nbocv_v[noc];

			diverg[icv] = grad_u[icv][0][0] + grad_u[icv][1][1] + grad_u[icv][2][2];

			strMag[icv] = 0.0;
			for (int i=0; i<3; i++)
				for (int j=0; j<3; j++) {
					if (i == j)  strMag[icv] += pow(0.5*(grad_u[icv][i][j] + grad_u[icv][j][i])-1./3.*diverg[icv], 2.0);
					else         strMag[icv] += pow(0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]), 2.0);
				}

			strMag[icv] = sqrt(2.0*strMag[icv]);
    	}
    }

    /*
     * Method:calcVorticity_AD
     * ----------------------------------------
     * Original code = calcVorticity_AD() in UgpWithCvCompFlowAD.h
     */
    virtual void calcVorticity1D_AD(const int icvCenter) {
#ifdef USE_MEM_SAVING_ADVAR
    	if (vortMag.empty()) {
    		cout << "Error in calcVorticity1D_AD() => vortMag has not been registered!" << endl;
    		throw(-1);
    	}
#else
    	if (vortMag == NULL) {
    		cout << "Error in calcVorticity1D_AD() => vortMag has not been registered!" << endl;
    		throw(-1);
    	}
#endif

    	int noc_f = nbocv_i[icvCenter];
    	int noc_l = nbocv_i[icvCenter+1] - 1;
    	for (int noc = noc_f; noc <= noc_l; noc++) {
    		int icv = nbocv_v[noc];

    		vortMag[icv] = 0.0;
    		for (int i=0; i<3; i++)
    			for (int j=0; j<3; j++)
    				vortMag[icv] += pow(0.5*(grad_u[icv][i][j] - grad_u[icv][j][i]), 2.0);

    		vortMag[icv] = sqrt(2.0*vortMag[icv]);
    	}
    }

    // =========================
    //  UTILITY FUNCTIONS
    // =========================
    /*
     * Method: getNumControlVar
     * ------------------------
     * Get the number of control variables (e.g. heat release rate)
     */
    int getNumControlVar() {
    	return getIntParam("N_PSALC_CONTROL_PARAMS", "0");
    }

    /*
     * Method: getBoundaryType
     * -----------------------
     * Get the boundary condition for N-S
     * If this method cannot find a proper BC, it will report an warning and return "WRONG"
     */
    inline BOUNDARY_TYPE getBoundaryType(const int ifa, Param *param) {
    	BOUNDARY_TYPE btype = WRONG;

    	if((ifa>=nfa_b && ifa<nfa) || (ifa>=nfa_b2 && ifa<nfa_b2gg)) {
    		btype = INTERNAL;
    	} else {
    		if (param->getString() == "HOOK")
    			btype = HOOK;
    		else if (param->getString() == "CBC")
    			btype = CBC;
    		else if (param->getString() == "CBC_SUBSONIC_INLET")
    			btype = CBC_SUBSONIC_INLET;
    		else if (param->getString() == "CBC_SUBSONIC_OUTLET")
    			btype = CBC_SUBSONIC_OUTLET;
    		else if (param->getString() == "SYMMETRY")
    			btype = SYMMETRY;
    		else if (param->getString() == "NEUMANN")
    			btype = NEUMANN;
    		else if (param->getString() == "WALL")
    			btype = WALL;
    		else {
    			cout<<"IkeWithModels_AD::getBoundaryType() => WARNING! Not a proper bc for ifa=" << ifa << " bc: " << param->getString() << endl;
    		}
    	}

    	return btype;
    }

    /*
     * Method: getBoundaryType
     * -----------------------
     * Get the boundary condition for N-S
     */
    inline BOUNDARY_TYPE getBoundaryType(const int ifa, FaZone* zone) {
    	Param *param;
    	getParam(param, zone->getName());

    	return getBoundaryType(ifa, param);
    }

    /*
     * Method: getBoundaryTypeScalar
     * -----------------------------
     * Get the scalar boundary condition for scalars
     * Original code = scalarZoneIsHook(), scalarZoneIsDirichlet(), and scalarZoneIsFlux() in UgpWithCvCompFlow.h
     * If this method cannot find a proper BC, it will return "OTHER_SCALAR".
     * Recommendation: It is highly recommended to use scalarZoneIsHook(), scalarZoneIsDirichlet(), and scalarZoneIsFlux() instead of this method.
     *                 They can also give the boundary value, whereas this method can provide boundary type only.
     * Note: By default, this method will return "OTHER_SCALAR" if the face is not in among "HOOK_SCALAR", "DIRICHLET_SCALAR", and "FLUX_SCALAR"
     */
    inline BOUNDARY_TYPE_SCALAR getBoundaryTypeScalar(const int ifa, FaZone* zone, const string& scalName) {
    	BOUNDARY_TYPE_SCALAR btypeScalar;

    	if((ifa>=nfa_b && ifa<nfa) || (ifa>=nfa_b2 && ifa<nfa_b2gg)) {
    		btypeScalar = INTERNAL_SCALAR;
    	} else {
    		string faName(zone->getName());
    		Param *param;
    		if(getParam(param, faName+"."+scalName)) {
    			string name = param->getString(1);
    			if((name == "HOOK") || (name == "Hook")|| (name == "hook"))
    				btypeScalar = HOOK_SCALAR;
    			else
    				btypeScalar = DIRICHLET_SCALAR;
    		} else {
    			if (getParam(param, faName+"."+scalName+".flux"))
    				btypeScalar = FLUX_SCALAR;
    			else
    				btypeScalar = OTHER_SCALAR;
    		}
    	}

    	return btypeScalar;
    }

	/*
	 * Method: checkScalarsMemAD
	 * -------------------------
	 * This is an original code of IkeWithPsALC
	 */
	virtual void checkScalarsMemAD() { /*empty*/ }

	/*
	 * Method: getDebugLevel
	 * ---------------------
	 * Obtain the debugging level.
	 * The higher is the level, shows more details on the screen
	 */
	int getDebugLevel() {
		static bool firstCall = true;

		if(firstCall) {
			if (!checkParam("DEBUG_LEVEL")) {
				ParamMap::add("DEBUG_LEVEL=0"); // add default values
				if (mpi_rank == 0)
					cout << "WARNING: added keyword \"DEBUG_LEVEL=0\" to parameter map!"<<endl;
			}
		}

		int debugLevel = getIntParam("DEBUG_LEVEL", "0");
		if(firstCall && mpi_rank==0)
			cout<<"> DEBUG_LEVEL = "<<debugLevel<<endl;

		firstCall = false;
		return debugLevel;
	}

	// =========================
	//  MEMBER VARIABLES
	// =========================

	// --------------------
	// Heat
	// --------------------
	double *lambda; // lambda will be allocated in the runPsALC() method
	REALQ *lambda_AD;

//#endif /* IKEUGPWITHCVCOMPFLOW_H_ */
