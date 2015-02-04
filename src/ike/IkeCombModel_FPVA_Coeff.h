#ifndef IKECOMBMODEL_FPVA_COEFF_AD_H
#define IKECOMBMODEL_FPVA_COEFF_AD_H

//#include "JOE/ADJOINT_FILES/UgpWithCvCompFlowAD.h"
//#include "JOE/ADJOINT_FILES/CombModel_FPVA_Coeff_AD.h"

#include "JOE/combModels/CombModel_Base.h"
#include "JOE/combModels/CombModel_FPVA_Coeff.h"

#ifdef USE_MEM_SAVING_ADVAR
	#define USE_MEM_SAVING_ADVAR_1D_ // Note: Actually, USE_MEM_SAVING_ADVAR is defined in UgpWithCvCompFlowAD.h (thus, it is readable here)
	                                 //       and USE_MEM_SAVING_ADVAR_1D_ is defined in IkeWithModels.h (thus, it is not readable here).
	                                 //       You should make USE_MEM_SAVING_ADVAR_1D_ available here to use the 1D memory saving technique.
	#include "ADvar.h"
#endif

#define IKECOMBMODEL_FPVA_COEFF_ERROR_CODE 1208

//#include "IkeWithPsALC.h"

template <class Chemtable>
//class IkeRansCombFPVA_Coeff_AD : virtual public UgpWithCvCompFlow_AD, public RansCombFPVA_Coeff<ChemtableAdaptiveLinear> {
class IkeRansCombFPVA_Coeff_AD : virtual public UgpWithCvCompFlow_AD, public RansCombFPVA_Coeff<ChemtableCartesianLinear> {
public:   // member vars
	string classID;

//	Chemtable myChemTable;
	bool chemtableAlreadyLoaded;

#ifdef USE_MEM_SAVING_ADVAR_1D_
	ADscalar<adouble>* ZMean;
	ADscalar<adouble>* ZVar;
	ADscalar<adouble>* CMean;

	ADscalar<adouble>* ZMean_diff;
	ADscalar<adouble>* ZVar_diff;
	ADscalar<adouble>* CMean_diff;

	ADvector<adouble>* grad_ZMean;
	ADvector<adouble>* grad_ZVar;
	ADvector<adouble>* grad_CMean;
#else
//	adouble *ZMean, *ZVar, *CMean;
//	adouble *ZMean_diff, *ZVar_diff, *CMean_diff;
//	adouble (*grad_ZMean)[3], (*grad_ZVar)[3], (*grad_CMean)[3];
#endif

public:
	// constructors
	IkeRansCombFPVA_Coeff_AD() {
		classID = "IkeRansCombFPVA_Coeff_AD";
		if (mpi_rank == 0)
			cout << classID<<"()" << endl;

		chemtableAlreadyLoaded = false;

#ifdef USE_MEM_SAVING_ADVAR_1D_
		ZMean = NULL;
		ZVar  = NULL;
		CMean = NULL;

		ZMean_diff = NULL;
		ZVar_diff  = NULL;
		CMean_diff = NULL;

		grad_ZMean = NULL;
		grad_ZVar  = NULL;
		grad_CMean = NULL;
#endif

		ZMean_Index = getScalarTransportIndex("ZMean");
		ZVar_Index  = getScalarTransportIndex("ZVar");
		CMean_Index = getScalarTransportIndex("CMean");
	}

	/*! \brief Unload myMixing and clean memory. */
	virtual ~IkeRansCombFPVA_Coeff_AD() {
		if (mpi_rank == 0)
			cout << "~"<<classID<<"()" << endl;
	}

public:
	// ===================================================
	// METHODS WHICH COME FROM RansCombFPVA_Coeff_AD
	// ===================================================
	/*
	 * Method: initialHookScalarRansCombModel
	 * --------------------------------------
	 * Wrapper function for RansCombFPVA_Coeff::initialHookScalarRansCombModel().
	 * You should call this function at the beginning of IkeWithPsALC to initialize the chemical table
	 */
	virtual void initialHookScalarRansCombModel() {
		if(!chemtableAlreadyLoaded)
			RansCombFPVA_Coeff<ChemtableCartesianLinear>::initialHookScalarRansCombModel();
		chemtableAlreadyLoaded = true;
	}

	/*
	 * Method: initialHookScalarRansCombModel1D_AD()
	 * ---------------------------------------------
	 * This method is doing the exact same thing as initialHookScalarRansCombModel_AD()
	 * but the implementation that works only for the ICVs stored in nbocv2ff.
	 * Original code = initialHookScalarRansCombModel_AD in RansCombFPVA_Coeff_AD.h
	 * Note: Actually the input argument "loadtable" is not used here.
	 */
	virtual void initialHookScalarRansCombModel1D_AD(int loadtable, vector<int> &nbocv2ff, bool &firstCallScalarComb) {
//		if(firstCall && loadtable==1) {
//			myChemTable.Load(getStringParam("CHEMTABLE_FILE"));
//
//			Preference = getDoubleParam("REFERENCE_PRESSURE", "0.0");
//			if (Preference == 0.0)
//				Preference = myChemTable.GetReferencePressure();
//		}
		if(!chemtableAlreadyLoaded) {
			if(mpi_rank==0)
				cerr<<"ERROR " << classID<<"::initialHookScalarRansCombModel1D_AD(): chemTable has not been loaded by RansCombFPVA_Coeff::initialHookScalarRansCombModel()."<<endl
				    <<"Note that "<<classID<<"::initialHookScalarRansCombModel() should call RansCombFPVA_Coeff::initialHookScalarRansCombModel()."<<endl;
			throw(IKECOMBMODEL_FPVA_COEFF_ERROR_CODE);
		}

		// connect pointers
		int nScal = scalarTranspEqVector.size();

		for(int i=0; i<nScal ; i++) {
			if(strcmp(scalarTranspEqVector_AD[i].name,"ZMean")==0){
                if(ZMean != NULL) {
                    cout<<"ERROR "<<classID<<"::initialHookScalarRansCombModel1D_AD(): ZMean is not NULL, icv="<<nbocv2ff[0]<<", mpi_rank="<<mpi_rank<<endl;
			        throw(IKECOMBMODEL_FPVA_COEFF_ERROR_CODE);
                }
				assert(scalarTranspEqVector_AD[i].phi.size() == nbocv2ff.size());

				ZMean      = &(scalarTranspEqVector_AD[i].phi) ;
				grad_ZMean = &(scalarTranspEqVector_AD[i].grad_phi) ;
				ZMean_diff = &(scalarTranspEqVector_AD[i].diff) ;
				if(mpi_rank==0 && firstCallScalarComb)
					cout<<"Connected Scalar Pointer "<<scalarTranspEqVector_AD[i].name<<endl;
			}
			if(strcmp(scalarTranspEqVector_AD[i].name,"ZVar")==0){
				assert(ZVar == NULL);
				assert(scalarTranspEqVector_AD[i].phi.size() == nbocv2ff.size());

				ZVar      = &(scalarTranspEqVector_AD[i].phi) ;
				grad_ZVar = &(scalarTranspEqVector_AD[i].grad_phi) ;
				ZVar_diff = &(scalarTranspEqVector_AD[i].diff) ;
				if(mpi_rank==0 && firstCallScalarComb)
					cout<<"Connected Scalar Pointer "<<scalarTranspEqVector_AD[i].name<<endl;
			}
			if(strcmp(scalarTranspEqVector_AD[i].name,"CMean")==0){
				assert(CMean == NULL);
				assert(scalarTranspEqVector_AD[i].phi.size() == nbocv2ff.size());

				CMean      = &(scalarTranspEqVector_AD[i].phi) ;
				grad_CMean = &(scalarTranspEqVector_AD[i].grad_phi) ;
				CMean_diff = &(scalarTranspEqVector_AD[i].diff) ;
				if(mpi_rank==0 && firstCallScalarComb)
					cout<<"Connected Scalar Pointer "<<scalarTranspEqVector_AD[i].name<<endl;
			}
		}

		// I don't know if commenting this out will cause ERROR or not...
		for(size_t i=0; i<nbocv2ff.size(); ++i) {
			int icv = nbocv2ff[i];
			RoM[icv] = UgpWithCvCompFlow::RoM[icv];
			gamma[icv] = UgpWithCvCompFlow::gamma[icv];
		}
		if(firstCallScalarComb) {
			for(size_t i=0; i<nbocv2ff.size(); ++i) {
				int icv = nbocv2ff[i];
				if(RoM[icv].value() <= 0.0 || isnan(RoM[icv].value())) {
					cerr<<"ERROR "<<classID<<"::initialHookScalarRansCombModel1D_AD(): RoM["<<icv<<"]="<<RoM[icv]<<endl;
					throw(IKECOMBMODEL_FPVA_COEFF_ERROR_CODE);
				}

				if(gamma[icv].value() <= 0.0 || isnan(gamma[icv].value())) {
					cerr<<"ERROR "<<classID<<"::initialHookScalarRansCombModel1D_AD(): gamma["<<icv<<"]="<<gamma[icv]<<endl;
					throw(IKECOMBMODEL_FPVA_COEFF_ERROR_CODE);
				}
			}
		}

		firstCallScalarComb = false;
	}
/*
	virtual void initialHookScalarRansCombModel1D_AD(int loadtable, vector<int> &nbocv2ff, bool &firstCallScalarComb) {
//		if(firstCall && loadtable==1) {
//			myChemTable.Load(getStringParam("CHEMTABLE_FILE"));
//
//			Preference = getDoubleParam("REFERENCE_PRESSURE", "0.0");
//			if (Preference == 0.0)
//				Preference = myChemTable.GetReferencePressure();
//		}
		if(!chemtableAlreadyLoaded) {
			if(mpi_rank==0)
				cerr<<"ERROR " << classID<<"::initialHookScalarRansCombModel1D_AD(): chemTable has not been loaded by RansCombFPVA_Coeff::initialHookScalarRansCombModel()."<<endl
				    <<"Note that "<<classID<<"::initialHookScalarRansCombModel() should call RansCombFPVA_Coeff::initialHookScalarRansCombModel()."<<endl;
			throw(IKECOMBMODEL_FPVA_COEFF_ERROR_CODE);
		}

		// connect pointers
		int nScal = scalarTranspEqVector.size();

		for(int i=0; i<nScal ; i++) {
#ifdef USE_MEM_SAVING_ADVAR
			if(strcmp(scalarTranspEqVector_AD[i].name,"ZMean")==0){
                if(ZMean != NULL) cout<<"ERROR "<<classID<<"::initialHookScalarRansCombModel1D_AD(): ZMean is not NULL, icv="<<nbocv2ff[0]<<", mpi_rank="<<mpi_rank<<endl;
				ZMean      = &(scalarTranspEqVector_AD[i].phi) ;
				grad_ZMean = &(scalarTranspEqVector_AD[i].grad_phi) ;
				ZMean_diff = &(scalarTranspEqVector_AD[i].diff) ;
				if(mpi_rank==0 && firstCallScalarComb)
					cout<<"Connected Scalar Pointer "<<scalarTranspEqVector_AD[i].name<<endl;
			}
			if(strcmp(scalarTranspEqVector_AD[i].name,"ZVar")==0){
				assert(ZVar == NULL);
				ZVar      = &(scalarTranspEqVector_AD[i].phi) ;
				grad_ZVar = &(scalarTranspEqVector_AD[i].grad_phi) ;
				ZVar_diff = &(scalarTranspEqVector_AD[i].diff) ;
				if(mpi_rank==0 && firstCallScalarComb)
					cout<<"Connected Scalar Pointer "<<scalarTranspEqVector_AD[i].name<<endl;
			}
			if(strcmp(scalarTranspEqVector_AD[i].name,"CMean")==0){
				assert(CMean == NULL);
				CMean      = &(scalarTranspEqVector_AD[i].phi) ;
				grad_CMean = &(scalarTranspEqVector_AD[i].grad_phi) ;
				CMean_diff = &(scalarTranspEqVector_AD[i].diff) ;
				if(mpi_rank==0 && firstCallScalarComb)
					cout<<"Connected Scalar Pointer "<<scalarTranspEqVector_AD[i].name<<endl;
			}
#else
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
#endif
		}

		// I don't know if commenting this out will cause ERROR or not...
		for(size_t i=0; i<nbocv2ff.size(); ++i) {
			int icv = nbocv2ff[i];
			RoM[icv] = UgpWithCvCompFlow::RoM[icv];
			gamma[icv] = UgpWithCvCompFlow::gamma[icv];
		}
		if(firstCallScalarComb) {
			for(size_t i=0; i<nbocv2ff.size(); ++i) {
				int icv = nbocv2ff[i];
				if(RoM[icv].value() <= 0.0 || isnan(RoM[icv].value())) {
					cerr<<"ERROR "<<classID<<"::initialHookScalarRansCombModel1D_AD(): RoM["<<icv<<"]="<<RoM[icv]<<endl;
					throw(IKECOMBMODEL_FPVA_COEFF_ERROR_CODE);
				}

				if(gamma[icv].value() <= 0.0 || isnan(gamma[icv].value())) {
					cerr<<"ERROR "<<classID<<"::initialHookScalarRansCombModel1D_AD(): gamma["<<icv<<"]="<<gamma[icv]<<endl;
					throw(IKECOMBMODEL_FPVA_COEFF_ERROR_CODE);
				}
			}
		}

		firstCallScalarComb = false;
	}
*/

	/*
	 * Method: finalHookScalarRansCombModel1D_AD()
	 * -------------------------------------------
	 *
	 */
	virtual void finalHookScalarRansCombModel1D_AD() {
		// debug level
		int debugLevel = getDebugLevel();

		if(!(*ZMean).empty())
			(*ZMean).clear();
		ZMean = NULL;
		if(!(*grad_ZMean).empty())
			(*grad_ZMean).clear();
		grad_ZMean = NULL;
		if(!(*ZMean_diff).empty())
			(*ZMean_diff).clear();
		ZMean_diff = NULL;
		if (debugLevel > 1 && mpi_rank == 0)
			cout<<classID<<"::finalHookScalarRansCombModel1D_AD(): Reset ZMean-related pointers as NULL "<<endl;

		if(!(*ZVar).empty())
			(*ZVar).clear();
		ZVar = NULL;
		if(!(*grad_ZVar).empty())
			(*grad_ZVar).clear();
		grad_ZVar = NULL;
		if(!(*ZVar_diff).empty())
			(*ZVar_diff).clear();
		ZVar_diff = NULL;
		if (debugLevel > 1 && mpi_rank == 0)
			cout<<classID<<"::finalHookScalarRansCombModel1D_AD(): Reset ZVar-related pointers as NULL "<<endl;

		if(!(*CMean).empty())
			(*CMean).clear();
		CMean = NULL;
		if(!(*grad_CMean).empty())
			(*grad_CMean).clear();
		grad_CMean = NULL;
		if(!(*CMean_diff).empty())
			(*CMean_diff).clear();
		CMean_diff = NULL;
		if (debugLevel > 1 && mpi_rank == 0)
			cout<<classID<<"::finalHookScalarRansCombModel1D_AD(): Reset CMean-related pointers as NULL "<<endl;
	}

	/*
	 * Method: initialHookScalarRansCombModel_AD()
	 * -------------------------------------------
	 * Note: Actually the input argument "loadtable" is not used here.
	 */
	virtual void initialHookScalarRansCombModel_AD(int loadtable) {
		if(!chemtableAlreadyLoaded) {
			if(mpi_rank==0)
				cerr<<"ERROR " << classID<<"::initialHookScalarRansCombModel1D_AD(): chemTable has not been loaded by RansCombFPVA_Coeff::initialHookScalarRansCombModel()."<<endl
				    <<"Note that "<<classID<<"::initialHookScalarRansCombModel() should call RansCombFPVA_Coeff::initialHookScalarRansCombModel()."<<endl;
			throw(IKECOMBMODEL_FPVA_COEFF_ERROR_CODE);
		}

#ifndef USE_MEM_SAVING_ADVAR
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
#endif
		if(mpi_rank==0)
			cout<<"WARNING "<<classID<<"::initialHookScalarRansCombModel_AD(): This method is not yet ready to use!"<<endl;
	}

	/*
	 * Method: calcStateVariables1D_AD()
	 * -------------------------------
	 * Calculate state variables (vel, press, temp, enthalpy, sos) for the neighbors of icvCenter
	 * Original code = calcStateVariables_AD() in UgpWithCvCompFlowAD.h
	 *                 calcStateVariables1D_AD() in IkeWithModels.cpp
	 */
	virtual void calcStateVariables1D_AD(const int icvCenter, ADscalar<adouble> &rho, ADvector<adouble> &rhou, ADscalar<adouble> &rhoE) {
		// Check if the 2-layer CSR structure has already been developed
		assert(nbocv2_i != NULL && !nbocv2_v.empty());

		int debugLevel = getDebugLevel();

		adouble E, cp;
		adouble kinecv = 0.0;

		// Note: IkeUgpWithCvCompFlow.h has CmeanSource_AD, chi_AD
		ComputeScalarDissipation1D_AD(icvCenter, chi_AD);

		// Calculate state variables only for the 2-layer neighbors of icvCenter
		int noc_f = nbocv2_i[icvCenter];
		int noc_l = nbocv2_i[icvCenter+1] - 1;

		for (int noc = noc_f; noc <= noc_l; noc++) {
			int icv_nbr = nbocv2_v[noc];

			// Read and interpolate coefficients from chemistry table and save source term for progress variable
			myChemTable.LookupCoeff_AD(RoM[icv_nbr], T0, E0, GAMMA0, AGAMMA, MU0, AMU, LOC0, ALOC, CmeanSource_AD[icv_nbr], (*ZMean)[icv_nbr], (*ZVar)[icv_nbr], (*CMean)[icv_nbr]);

			for (int i=0; i<3; i++)
				vel[icv_nbr][i] = rhou[icv_nbr][i] / rho[icv_nbr];

			if (kine != NULL)
				kinecv = (*kine)[icv_nbr];

			E = rhoE[icv_nbr] / rho[icv_nbr] - 0.5 * vecDotVec3d_AD(vel[icv_nbr], vel[icv_nbr]) - kinecv;

			// Compute mixture temperature from energy and species mass fractions
			if (AGAMMA == 0.0)
				temp[icv_nbr] = T0 + (GAMMA0 - 1.0) / RoM[icv_nbr] * (E - E0);
			else
				temp[icv_nbr] = T0 + (GAMMA0 - 1.0) / AGAMMA * (exp(AGAMMA * (E - E0) / RoM[icv_nbr]) - 1.0);

			// Temperature clipping
			temp[icv_nbr] = max(temp[icv_nbr], Tminimum);

			// Compute mixture enthalpy and pressure
			enthalpy[icv_nbr] = E + RoM[icv_nbr] * temp[icv_nbr];
			press[icv_nbr] = rho[icv_nbr] * RoM[icv_nbr] * temp[icv_nbr];

			// Compute mixture viscosity, heat conductivity, gamma, speed of sound
			muLam_AD[icv_nbr] = MU0 * pow(temp[icv_nbr] / T0, 0.7);
			LambdaOverCp_AD[icv_nbr] = LOC0 * pow(temp[icv_nbr] / T0, 0.62);
			gamma[icv_nbr] = GAMMA0 + AGAMMA * (temp[icv_nbr] - T0);
			sos[icv_nbr] = sqrt(gamma[icv_nbr] * press[icv_nbr] / rho[icv_nbr]);

			if ((isnan(press[icv_nbr].value()) || (press[icv_nbr]<= 0.0) ))
			{
				cout << "WARNING! " << classID << "::calcStateVariables1D_AD(): Negative pressure at xcv:"
				     << x_cv[icv_nbr][0] << ", " << x_cv[icv_nbr][1]  << ", " << x_cv[icv_nbr][2] << endl;
				if(debugLevel > 0) {
					stringstream ss;
					ss<<std::scientific;
					ss<<std::setprecision(3);
					ss<<"                                    rhoE="<<rhoE[icv_nbr]
					  <<", rho="<<rho[icv_nbr]
					  <<", rhou=("<<rhou[icv_nbr][0]<<","<<rhou[icv_nbr][1]<<","<<rhou[icv_nbr][2]<<")"
					  <<", gamma="<<gamma[icv_nbr]
					  <<", kinecv="<<kinecv<<endl;
				}
				throw(-1);
			}
#ifdef PRESSURE_SCALING_FPVA_COEFF
			// rescale source term
			CmeanSource_AD[icv_nbr] = CmeanSource_AD[icv_nbr]*pow( press[icv_nbr]/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
#endif
		}

		if(debugLevel > 0) {
			if (gamma[icvCenter].value()<MACHINE_EPS) {
				cout<<"ERROR! "<<classID<<"::calcStateVariables1D_AD(): gamma["<<icvCenter<<"] is zero or negative = "<<gamma[icvCenter]<<endl;
				throw(IKECOMBMODEL_FPVA_COEFF_ERROR_CODE);
			}
			if (RoM[icvCenter].value()<MACHINE_EPS) {
				cout<<"ERROR! "<<classID<<"::calcStateVariables1D_AD(): RoM["<<icvCenter<<"] is zero or negative = "<<RoM[icvCenter]<<endl;
				throw(IKECOMBMODEL_FPVA_COEFF_ERROR_CODE);
			}
		}
	}

	/*
	 * Method: calcStateVariables_AD()
	 * -------------------------------
	 *
	 */
	virtual void calcStateVariables_AD(adouble *rho, adouble (*rhou)[3], adouble *rhoE) {
		adouble E, cp;
		adouble kinecv = 0.0;

		ComputeScalarDissipation_AD(chi_AD);

		for (int icv = 0; icv < ncv_gg; icv++)
		{

			// Read and interpolate coefficients from chemistry table and save source term for progress variable
			myChemTable.LookupCoeff_AD(RoM[icv], T0, E0, GAMMA0, AGAMMA, MU0, AMU, LOC0, ALOC, CmeanSource_AD[icv], (*ZMean)[icv], (*ZVar)[icv], (*CMean)[icv]);

			for (int i=0; i<3; i++)
				vel[icv][i] = rhou[icv][i]/rho[icv];

			if (kine != NULL)
				kinecv = (*kine)[icv];

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
						"  rhoE = " << rhoE[icv] << "  1/2u^2 = " << 0.5 * vecDotVec3d_AD(vel[icv], vel[icv]) << "  Z = " << (*ZMean)[icv] << endl;
				throw(-1);
			}

#ifdef PRESSURE_SCALING_FPVA_COEFF
			// rescale source term
			CmeanSource_AD[icv] = CmeanSource_AD[icv]*pow( press[icv]/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
#endif      
		}
	}

#ifdef MIXING_VISCOSITY  
	/*
	 * Method: calcMaterialProperties1D_AD
	 * -----------------------------------
	 * Compute laminar viscosity mul_fa and heat conductivity lamOcp_fa at the faces.
	 * Original code: 1. calcMaterialProperties_AD() in CombModel_FPVA_Coeff_AD.h
	 *                2. calcMaterialProperties1D_AD() in IkeWithModel.cpp
	 */
	virtual void calcMaterialProperties1D_AD(const int icvCenter, vector<int>& faInternal, vector<int>& faBoundary, ADscalar<adouble> &rho, ADvector<adouble> &rhou, ADscalar<adouble> &rhoE) {
		if (mu_ref > 0.0) {
			// internal faces
			for (size_t index = 0; index < faInternal.size(); index++) {
				int ifa = faInternal[index];
				if((ifa>=nfa_b && ifa<nfa) || (ifa >= nfa_b2 && ifa < nfa_b2gg)) {
					 /* Note: 0 ~ nfa_b-1           : boundary faces
					  *       nfa_b ~ nfa_bpi-1     :
					  *       nfa_bpi ~ nfa-1       :
					  *       nfa ~ nfa_b2          : additional boundary faces between g1:f2 cvs
					  *       nfa_b2 ~ nfa_b2g-1    : faces between g1:g1 cvs
					  *       nfa_b2g-1 ~ nfa_b2gg-1: faces between g1:g2 cvs
					  */
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

	virtual void calcMaterialProperties_AD(adouble *rho, adouble (*rhou)[3], adouble *rhoE) {
		if (mu_ref > 0.0) {
			// internal faces
			for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) { // Ensure we sweep through internal faces of ghosts too
				if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg)) {
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

	/*
	 * Method: calcThermoProp_T_AD
	 * ---------------------------
	 * brief Compute for a given temperature, density and scalars
	 */
	virtual void calcThermoProp_T_AD(adouble &p, adouble &h, adouble &R, adouble &gam, adouble &rho, adouble &T, adouble *Scal, int nScal) {
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

	/*
	 * Method: calcThermoProp_p_AD
	 * ---------------------------
	 * brief Compute for a given pressure, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
	 */
	virtual void calcThermoProp_p_AD(adouble &T, adouble &h, adouble &R, adouble &gam, adouble &rho, adouble &p, adouble *Scal, int nScal) {
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

	/*
	 * Method: calcMuLam_AD
	 * --------------------
	 *
	 */
	virtual adouble calcMuLam_AD(int icv) {
		return muLam_AD[icv];
	}

	virtual adouble calcMuLam_AD(adouble temp) {
		cout << "ERROR " << classID<<"::calcMuLam(double temp) does not work! ####" << endl;
		throw(-1);
	}

	/*
	 * Method: ComputeBCProperties1D_T_AD
	 * --------------------------------
	 * brief Compute for a given temperature the properties of the mixture at a face.
	 * Original code = 1. ComputeBCProperties_T_AD in CombModel_FPVA_Coeff_AD.h
	 *                 2. ComputeBCProperties1D_T_AD in IkeWithModels.cpp
	 */
	virtual void ComputeBCProperties1D_T_AD(const int ifa) {
		int icv1 = cvofa[ifa][1];
		ComputeProperties_T_AD(enthalpy[icv1], RoM[icv1], gamma[icv1], temp[icv1], (*ZMean)[icv1], (*ZVar)[icv1], (*CMean)[icv1]);

		if (mu_ref > 0.0) {
			ComputeProperties_vis_AD(mul_fa[ifa], lamOcp_fa[ifa], temp[icv1], (*ZMean)[icv1], (*ZVar)[icv1], (*CMean)[icv1]);
		}
	}

	/*
	 * Method: ComputeBCProperties_T_AD
	 * --------------------------------
	 *
	 */
	virtual void ComputeBCProperties_T_AD(FaZone *zone) {
		for (int index = 0; index < zone->faVec.size(); ++index) {
			int ifa = zone->faVec[index];
			int icv1 = cvofa[ifa][1];
			ComputeProperties_T_AD(enthalpy[icv1], RoM[icv1], gamma[icv1], temp[icv1], (*ZMean)[icv1], (*ZVar)[icv1], (*CMean)[icv1]);
		}
		if (mu_ref > 0.0)
		{
			// for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
			for (int index = 0; index < zone->faVec.size(); ++index) {
				int ifa = zone->faVec[index];
				int icv1 = cvofa[ifa][1];
				ComputeProperties_vis_AD(mul_fa[ifa], lamOcp_fa[ifa], temp[icv1], (*ZMean)[icv1], (*ZVar)[icv1], (*CMean)[icv1]);
			}
		}
	}

	/*
	 * Method: ComputeBCProperties1D_H_AD
	 * ----------------------------------
	 * brief Compute for a given enthalpy the properties of the mixture at a face.
	 * Original code = 1. ComputeBCProperties_H_AD in CombModel_FPVA_Coeff_AD.h
	 *                 2. ComputeBCProperties1D_H_AD in IkeWithModels.cpp
	 */
	 virtual void ComputeBCProperties1D_H_AD(const int ifa) {
		 int icv1 = cvofa[ifa][1];

		 if (mu_ref > 0.0) {
			 ComputeProperties_vis_AD(mul_fa[ifa], lamOcp_fa[ifa], temp[icv1], (*ZMean)[icv1], (*ZVar)[icv1], (*CMean)[icv1]);
		 }
	 }

	virtual void ComputeBCProperties_H_AD(FaZone *zone) {
		for (int index = 0; index < zone->faVec.size(); ++index) {
			int ifa = zone->faVec[index];
			int icv1 = cvofa[ifa][1];
			ComputeProperties_H_AD(temp[icv1], RoM[icv1], gamma[icv1], enthalpy[icv1], (*ZMean)[icv1], (*ZVar)[icv1], (*CMean)[icv1]);
		}
		if (mu_ref > 0.0) {
			// for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++) {
			for (int index = 0; index < zone->faVec.size(); ++index) {
				int ifa = zone->faVec[index];
				int icv1 = cvofa[ifa][1];
				ComputeProperties_vis_AD(mul_fa[ifa], lamOcp_fa[ifa], temp[icv1], (*ZMean)[icv1], (*ZVar)[icv1], (*CMean)[icv1]);
			}
		}
	}

	/*
	 * Method: ComputeProperties_T_AD
	 * ------------------------------
	 * brief Compute for a given temperature the properties of the mixture.
	 * The enthalpy (chemical + sensible) is a function of temperature and mixture composition.
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
	void ComputeProperties_T_AD(adouble &h, adouble &RoM, adouble &gam, adouble T, adouble zm, adouble zv, adouble cm) {
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

	/*
	 * Method: ComputeProperties_H_AD
	 * ------------------------------
	 *
	 */
	void ComputeProperties_H_AD(adouble &T, adouble &RoM, adouble &gam, adouble h, adouble zm, adouble zv, adouble cm) {
		adouble dummy;

		myChemTable.LookupCoeff_AD(RoM, T0, E0, GAMMA0, AGAMMA, MU0, AMU, LOC0, ALOC, dummy, zm, zv, cm);
		T = (h - E0) / RoM;
		// Note: Solve T by using a Newton's method !!!
		T = SolveNewtonTemperature_H_AD(RoM, T0, E0, GAMMA0, AGAMMA, h, T);
		T = max(T, Tminimum);
		gam = GAMMA0 + AGAMMA * (T - T0);
		// mu = MU0 * pow(T / T0, 0.7);
		// lamOcp = LOC0 * pow(T / T0, 0.62);
	}

	/*
	 * Method: SolveNewtonTemperature_H_AD
	 * -----------------------------------
	 * Solve for Temperature by using Newton's method
	 */
	adouble SolveNewtonTemperature_H_AD(adouble R_given, adouble t0_given, adouble e0_given, adouble gamma0_given, adouble agamma_given, adouble h_given, adouble Tguess_given) {
		adouble Tresidual_AD;

		// The original code from Karthik uses only AD during the following Newton-iteration.
		// However, if AD stores all the operations during the Newton-interation, the data size can be gigantic.
		// Thus, here we will just use normal double variables until the Newton iteration coverges, then convert the final operations to AD.
		double R  = R_given.value();
		double t0 = t0_given.value();
		double e0 = e0_given.value();
		double gamma0 = gamma0_given.value();
		double agamma = agamma_given.value();
		double h  = h_given.value();
		double Tguess = Tguess_given.value();

		int    itermax = 50, iter;
		double errormax = 1.0e-7;
		double Tresidual, hh, cpp, error;

		iter = 0;
		while (iter < itermax) {
			++iter;
			if (agamma == 0.0)
				hh = e0 + R / (gamma0 - 1.0) * (Tguess - t0) + R * Tguess;
			else
				hh = e0 + R / agamma * log(1.0 + agamma * (Tguess - t0) / (gamma0 - 1.0)) + R * Tguess;
			cpp = R * (1.0 + 1.0 / (gamma0 + agamma * (Tguess - t0) - 1.0));

			Tresidual = Tguess - (hh - h) / cpp;

			error = Tguess - Tresidual;
			if (fabs(error) <= errormax) {
				adouble hh_AD, cpp_AD;
				if (agamma == 0.0)
					hh_AD = e0_given + R_given / (gamma0_given - 1.0) * (Tguess - t0_given) + R_given * Tguess;
				else
					hh_AD = e0_given + R_given / agamma_given * log(1.0 + agamma_given * (Tguess - t0_given) / (gamma0_given - 1.0)) + R_given * Tguess;
				cpp_AD = R_given * (1.0 + 1.0 / (gamma0_given + agamma_given * (Tguess - t0_given) - 1.0));

				Tresidual_AD = Tguess - (hh_AD - h_given) / cpp_AD;
				return Tresidual_AD;
			}

			Tguess = Tresidual;
		}

		/* if the iteration has not converged, exit with error */
		cerr << "ERROR " << classID<<"::SolveNewtonTemperature_H_AD(): Computation of temperature in Newton iteration has not converged: Tguess=" << Tguess + error
				<< ", Tresidual=" << Tresidual << ", enthalpy=" << h << ", #iterations=" << iter << " ###" << endl;
		throw(-1);
	}
//	adouble SolveNewtonTemperature_H_AD(adouble R, adouble t0, adouble e0, adouble gamma0, adouble agamma, adouble h, adouble Tguess) {
//		static bool firstCall = true;
//		if(firstCall && mpi_rank==0)
//			cout<<"WARNING! "<<classID<<"::SolveNewtonTemperature_H_AD() is running an iterative method!"<<endl;
//		firstCall = false;
//
//		int    itermax = 50, iter;
//		double errormax = 1.0e-7;
//		adouble Tresidual, hh, cpp, error;
//
//		iter = 0;
//		do
//		{
//			++iter;
//			if (agamma == 0.0)
//				hh = e0 + R / (gamma0 - 1.0) * (Tguess - t0) + R * Tguess;
//			else
//				hh = e0 + R / agamma * log(1.0 + agamma * (Tguess - t0) / (gamma0 - 1.0)) + R * Tguess;
//			cpp = R * (1.0 + 1.0 / (gamma0 + agamma * (Tguess - t0) - 1.0));
//
//			Tresidual = Tguess - (hh - h) / cpp;
//
//			error = Tguess - Tresidual;
//			if (fabs(error) <= errormax)
//				return Tresidual;
//
//			Tguess = Tresidual;
//
//		}
//		while (iter < itermax);
//
//		/* if the iteration has not converged, exit with error */
//		cerr << "ERROR " << classID<<"::SolveNewtonTemperature_H_AD(): Computation of temperature in Newton iteration has not converged: Tguess=" << Tguess + error
//				<< ", Tresidual=" << Tresidual << ", enthalpy=" << h << ", #iterations=" << iter << " ###" << endl;
//		throw(-1);
//	}

	/*
	 * Method: ComputeProperties_vis_AD
	 * --------------------------------
	 * Exactly same as ComputeProperties_vis() in CombModel_FPVA_Coeff_AD.h,
	 * but the name of the method is changed to become more clear about using AD.
	 */
	void ComputeProperties_vis_AD(adouble &mu, adouble &lamOcp, adouble T, adouble zm, adouble zv, adouble cm) {
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

	/*
	 * Method: ComputeScalarDissipation1D_AD
	 * -------------------------------------
	 *
	 */
	void ComputeScalarDissipation1D_AD(const int icvCenter, ADscalar<adouble> &chi_AD) {
		if ((turbModel == KOM) || (turbModel == KOMSST)) {
			int om_index = getScalarTransportIndex("omega");
			chi_AD[icvCenter] = Cchi * Cmu * scalarTranspEqVector_AD[om_index].phi[icvCenter];
		} else if (turbModel == KEPS) {
			int epsi_index = getScalarTransportIndex("eps");
			chi_AD[icvCenter] = Cchi * scalarTranspEqVector_AD[epsi_index].phi[icvCenter] / (*kine)[icvCenter];
		} else if (turbModel == SA) {
			calcVorticity1D_AD(icvCenter);
			chi_AD[icvCenter] = Cchi * sqrt(Cmu) * vortMag[icvCenter];
		} else {
			chi_AD[icvCenter] = 0.0;
		}
	}

	/*
	 * Method: ComputeScalarDissipation_AD
	 * -----------------------------------
	 *
	 */
#ifdef USE_MEM_SAVING_ADVAR
	void ComputeScalarDissipation_AD(ADscalar<adouble> &chi_AD) {
		if ((turbModel == KOM) || (turbModel == KOMSST)) {
			int om_index = getScalarTransportIndex("omega");
			for (int icv = 0; icv < ncv_gg; icv++)
				chi_AD[icv] = Cchi * Cmu * scalarTranspEqVector_AD[om_index].phi[icv];
		} else if (turbModel == KEPS) {
			int epsi_index = getScalarTransportIndex("eps");
			for (int icv = 0; icv < ncv_gg; icv++)
				chi_AD[icv] = Cchi * scalarTranspEqVector_AD[epsi_index].phi[icv] / (*kine)[icv];
		} else if (turbModel == SA) {
			calcVorticity();
			for (int icv = 0; icv < ncv_gg; icv++)
				chi_AD[icv] = Cchi * sqrt(Cmu) * vortMag[icv];
		} else {
			for (int icv = 0; icv < ncv_gg; icv++)
				chi_AD[icv] = 0.0;
		}
	}
#else
//	void ComputeScalarDissipation_AD(adouble *chi_AD) {
//		if ((turbModel == KOM) || (turbModel == KOMSST)) {
//			int om_index = getScalarTransportIndex("omega");
//			for (int icv = 0; icv < ncv_gg; icv++)
//				chi_AD[icv] = Cchi * Cmu * scalarTranspEqVector_AD[om_index].phi[icv];
//		} else if (turbModel == KEPS) {
//			int epsi_index = getScalarTransportIndex("eps");
//			for (int icv = 0; icv < ncv_gg; icv++)
//				chi_AD[icv] = Cchi * scalarTranspEqVector_AD[epsi_index].phi[icv] / kine[icv];
//		} else if (turbModel == SA) {
//			calcVorticity();
//			for (int icv = 0; icv < ncv_gg; icv++)
//				chi_AD[icv] = Cchi * sqrt(Cmu) * vortMag[icv];
//		} else {
//			for (int icv = 0; icv < ncv_gg; icv++)
//				chi_AD[icv] = 0.0;
//		}
//	}
#endif

	/*
	 * Method: diffusivityHookScalarRansComb1D_AD
	 * ------------------------------------------
	 * Original code = 1. diffusivityHookScalarRansComb_AD() in CombModel_FPVA_Coeff_AD.h
	 *                 2. diffusivityHookScalarRansComb1D_AD() in IkeUgpWithCvCompFlow.h
	 */
	virtual void diffusivityHookScalarRansComb1D_AD(const int ifa, FaZone* zone, const string &name) {
		/* Note: 0 ~ nfa_b-1           : boundary faces
		 *       nfa_b ~ nfa_bpi-1     :
		 *       nfa_bpi ~ nfa-1       :
		 *       nfa ~ nfa_b2          : additional boundary faces between g1:f2 cvs
		 *       nfa_b2 ~ nfa_b2g-1    : faces between g1:g1 cvs
		 *       nfa_b2g-1 ~ nfa_b2gg-1: faces between g1:g2 cvs
		 */

		if (name == "ZMean") {
			// internal faces
			if((ifa>=nfa_b && ifa<nfa) || (ifa>=nfa_b2 && ifa<nfa_b2gg)) {
				(*ZMean_diff)[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;
			}
			// boundary faces
			else {
				if (zone->getKind() == FA_ZONE_BOUNDARY) {
					if (getBoundaryType(ifa, zone) == WALL) {
						(*ZMean_diff)[ifa] = lamOcp_fa[ifa];
					} else {
						(*ZMean_diff)[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;
					}
				}
			}
		}
		if (name == "ZVar") {
			// internal faces
			if((ifa>=nfa_b && ifa<nfa) || (ifa>=nfa_b2 && ifa<nfa_b2gg)) {
				(*ZVar_diff)[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;
			}
			// boundary faces
			else {
				if (zone->getKind() == FA_ZONE_BOUNDARY) {
					if (getBoundaryType(ifa, zone) == WALL) {
						(*ZVar_diff)[ifa] = lamOcp_fa[ifa];
					} else {
						(*ZVar_diff)[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;
					}
				}
			}
		}
		if (name == "CMean") {
			// internal faces
			if((ifa>=nfa_b && ifa<nfa) || (ifa>=nfa_b2 && ifa<nfa_b2gg)) {
				(*CMean_diff)[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;
			}
			// boundary faces
			else {
				if (zone->getKind() == FA_ZONE_BOUNDARY) {
					if (getBoundaryType(ifa, zone) == WALL) {
						(*CMean_diff)[ifa] = lamOcp_fa[ifa];
					} else {
						(*CMean_diff)[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;
					}
				}
			}
		}
	}

	/*
	 * Method: diffusivityHookScalarRansComb_AD
	 * ----------------------------------------
	 *
	 */
	virtual void diffusivityHookScalarRansComb_AD(const string &name) {
		if (name == "ZMean") {
			// internal faces
			for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) // Ensure we sweep through internal faces of ghosts too
				if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg))
					(*ZMean_diff)[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;

			// boundary faces
			for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
				if (zone->getKind() == FA_ZONE_BOUNDARY) {
					if (zoneIsWall(zone->getName()))
						for (int index = 0; index < zone->faVec.size(); ++index) {
							int ifa = zone->faVec[index];
							(*ZMean_diff)[ifa] = lamOcp_fa[ifa];
						}
					else
						for (int index = 0; index < zone->faVec.size(); ++index) {
							int ifa = zone->faVec[index];
							(*ZMean_diff)[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;
						}
				}
		}
		if (name == "ZVar") {
			// internal faces
			for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) // Ensure we sweep through internal faces of ghosts too
				if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg))
					(*ZVar_diff)[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;

			// boundary faces
			for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
				if (zone->getKind() == FA_ZONE_BOUNDARY) {
					if (zoneIsWall(zone->getName()))
						for (int index = 0; index < zone->faVec.size(); ++index) {
							int ifa = zone->faVec[index];
							(*ZVar_diff)[ifa] = lamOcp_fa[ifa];
						}
					else
						for (int index = 0; index < zone->faVec.size(); ++index) {
							int ifa = zone->faVec[index];
							(*ZVar_diff)[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;
						}
				}
		}
		if (name == "CMean") {
			// internal faces
			for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) // Ensure we sweep through internal faces of ghosts too
				if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg))
					(*CMean_diff)[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;

			// boundary faces
			for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
				if (zone->getKind() == FA_ZONE_BOUNDARY) {
					if (zoneIsWall(zone->getName()))
						for (int index = 0; index < zone->faVec.size(); ++index) {
							int ifa = zone->faVec[index];
							(*CMean_diff)[ifa] = lamOcp_fa[ifa];
						}
					else
						for (int index = 0; index < zone->faVec.size(); ++index) {
							int ifa = zone->faVec[index];
							(*CMean_diff)[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;
						}
				}
		}
	}

	/*
	 * Method: sourceHookScalarRansComb_new1D_AD
	 * -----------------------------------------
	 * Original code = 1. sourceHookScalarRansComb_new_AD() in CombModel_FPVA_Coeff_AD.h
	 *                 2. sourceHookScalarRansTurb_new1D_AD in IkeUgpWithCvCompFlow.h
	 */
	virtual void sourceHookScalarRansComb_new1D_AD(const int icv, adouble &rhs, double *A, const string &name, int flagImplicit) {
		if (name == "ZVar") {
			adouble src = 2.0*InterpolateAtCellCenterFromFaceValues_AD(mut_fa, icv)/Schmidt_turb_Z2 * vecDotVec3d_AD((*grad_ZMean)[icv], (*grad_ZMean)[icv]) - chi_AD[icv] * rho_AD[icv] * (*ZVar)[icv];
			rhs += src * cv_volume[icv];

			if (flagImplicit) {
				int noc00 = nbocv_i[icv];
				double dsrcdphi = - chi_AD[icv].value();
				A[noc00] -= dsrcdphi*cv_volume[icv];
			}
		}

		if (name == "CMean") {
			rhs += rho_AD[icv] * CmeanSource_AD[icv] * cv_volume[icv];

			if (flagImplicit) {
				// Derivative with respect to CMean
				adouble E = enthalpy[icv] - RoM[icv] * temp[icv];
				double E1 = E.value();
				double ZM = (*ZMean)[icv].value();
				double ZV = (*ZMean)[icv].value();
				double delta_CM = 1.0e-5;
				double CMp = (*CMean)[icv].value() + delta_CM;
				double CMm = (*CMean)[icv].value() - delta_CM;
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

	/*
	 * Method: sourceHookScalarRansComb_new_AD
	 * ---------------------------------------
	 *
	 */
	virtual void sourceHookScalarRansComb_new_AD(adouble *rhs, double *A, const string &name, int flagImplicit) {
		if (name == "ZVar") {
			for (int icv = 0; icv < ncv_g; icv++) {
				adouble src = 2.0*InterpolateAtCellCenterFromFaceValues_AD(mut_fa, icv)/Schmidt_turb_Z2 * vecDotVec3d_AD((*grad_ZMean)[icv], (*grad_ZMean)[icv]) - chi_AD[icv] * rho_AD[icv] * (*ZVar)[icv];
				rhs[icv] += src * cv_volume[icv];
			}

			if (flagImplicit) {
				for (int icv = 0; icv < ncv; icv++) {
					int noc00 = nbocv_i[icv];
					double dsrcdphi = - chi_AD[icv].value();
					A[noc00] -= dsrcdphi*cv_volume[icv];
				}
			}
		}

		if (name == "CMean") {
			for (int icv = 0; icv < ncv_g; icv++) {
				rhs[icv] += rho_AD[icv] * CmeanSource_AD[icv] * cv_volume[icv];
			}

			if (flagImplicit) {
				// Derivative with respect to CMean
				for (int icv = 0; icv < ncv; icv++) {
					adouble E = enthalpy[icv] - RoM[icv] * temp[icv];
					double E1 = E.value();
					double ZM = (*ZMean)[icv].value();
					double ZV = (*ZMean)[icv].value();
					double delta_CM = 1.0e-5;
					double CMp = (*CMean)[icv].value() + delta_CM;
					double CMm = (*CMean)[icv].value() - delta_CM;
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

	/*
	 * Method: sourceHookRansCombCoupled1D_AD
	 * -----------------------------------------
	 * Original code = 1. sourceHookRansCombCoupled_AD() in CombModel_FPVA_Coeff_AD.h
	 *                 2. sourceHookRansCombCoupled1D_AD in IkeUgpWithCvCompFlow.h
	 */
	virtual void sourceHookRansCombCoupled1D_AD(const int icv, adouble *rhs, double ***A, int flagImplicit) {
		string funcID = "sourceHookRansCombCoupled1D_AD";

		int ZMean_Coupled_Index = getScalarTransportIndex("ZMean");
		if (ZMean_Coupled_Index == -1) {
			cerr << "ERROR "<<classID<<"::"<<funcID<<"(): ZMean_Coupled_Index is " << ZMean_Coupled_Index << " !###" << endl;
			throw(-1);
		}
		int ZVar_Coupled_Index = getScalarTransportIndex("ZVar");
		if (ZVar_Coupled_Index == -1) {
			cerr << "ERROR "<<classID<<"::"<<funcID<<"(): ZVar_Coupled_Index is " << ZVar_Coupled_Index << " !###" << endl;
			throw(-1);
		}
		int CMean_Coupled_Index = getScalarTransportIndex("CMean");
		if (CMean_Coupled_Index == -1) {
			cerr << "ERROR "<<classID<<"::"<<funcID<<"(): CMean_Coupled_Index is " << CMean_Coupled_Index << " !###" << endl;
			throw(-1);
		}

		rhs[5+ZVar_Coupled_Index ] += (2.0 * InterpolateAtCellCenterFromFaceValues_AD(mut_fa, icv) / Schmidt_turb_Z2 * vecDotVec3d_AD((*grad_ZMean)[icv], (*grad_ZMean)[icv]) - chi_AD[icv] * rho_AD[icv] * (*ZVar)[icv]) * cv_volume[icv];
		rhs[5+CMean_Coupled_Index] += rho_AD[icv] * CmeanSource_AD[icv] * cv_volume[icv];

		if (flagImplicit && icv<ncv) {
			// ********************* Derivative computation ***********************

			double dCmeanSource_dZM, dCmeanSource_dZV, dCmeanSource_dCM,E;
			double delta_ZM = 1.0e-5;
			double delta_ZV = 1.0e-8;
			double delta_CM = 1.0e-5;
			double pressp, pressm, CmeanSourcep, CmeanSourcem;
			double ZV = (*ZVar)[icv].value();
			double ZM = (*ZMean)[icv].value();
			double CM = (*CMean)[icv].value();

			E = enthalpy[icv].value() - RoM[icv].value() * temp[icv].value();

			// Derivative with respect to ZMean
			double ZMp = (*ZMean)[icv].value() + delta_ZM;
			double ZMm = (*ZMean)[icv].value() - delta_ZM;

			pressure_scalSource_AD(pressp, CmeanSourcep, rho[icv], E, ZMp, ZV, CM);
			pressure_scalSource_AD(pressm, CmeanSourcem, rho[icv], E, ZMm, ZV, CM);
#ifdef PRESSURE_SCALING_FPVA_COEFF
			CmeanSourcep *= pow( pressp/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
			CmeanSourcem *= pow( pressm/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
#endif
			dCmeanSource_dZM = (CmeanSourcep - CmeanSourcem) / (2.0 * delta_ZM);

			// Derivative with respect to ZVar
			double ZVp = (*ZVar)[icv].value() + delta_ZV;
			double ZVm = (*ZVar)[icv].value() - delta_ZV;
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
			A[noc00][0][5+CMean_Coupled_Index] -= (CmeanSource_AD[icv].value() - (*ZMean)[icv].value() * dCmeanSource_dZM - (*ZVar)[icv].value() * dCmeanSource_dZV - (*CMean)[icv].value() * dCmeanSource_dCM) * cv_volume[icv];
		}
	}

	/*
	 * Method: sourceHookRansCombCoupled_AD
	 * ------------------------------------
	 *
	 */
	virtual void sourceHookRansCombCoupled_AD(adouble **rhs, double ***A, int flagImplicit) {
		string funcID = "sourceHookRansCombCoupled_AD";

		int ZMean_Coupled_Index = getScalarTransportIndex("ZMean");
		if (ZMean_Coupled_Index == -1) {
			cerr << "ERROR "<<classID<<"::"<<funcID<<"(): ZMean_Coupled_Index is " << ZMean_Coupled_Index << " !###" << endl;
			throw(-1);
		}
		int ZVar_Coupled_Index = getScalarTransportIndex("ZVar");
		if (ZVar_Coupled_Index == -1) {
			cerr << "ERROR "<<classID<<"::"<<funcID<<"(): ZVar_Coupled_Index is " << ZVar_Coupled_Index << " !###" << endl;
			throw(-1);
		}
		int CMean_Coupled_Index = getScalarTransportIndex("CMean");
		if (CMean_Coupled_Index == -1) {
			cerr << "ERROR "<<classID<<"::"<<funcID<<"(): CMean_Coupled_Index is " << CMean_Coupled_Index << " !###" << endl;
			throw(-1);
		}

		for (int icv = 0; icv < ncv_g; icv++) {
			rhs[icv][5+ZVar_Coupled_Index ] += (2.0 * InterpolateAtCellCenterFromFaceValues_AD(mut_fa, icv) / Schmidt_turb_Z2 * vecDotVec3d_AD((*grad_ZMean)[icv], (*grad_ZMean)[icv]) - chi_AD[icv] * rho_AD[icv] * (*ZVar)[icv]) * cv_volume[icv];
			rhs[icv][5+CMean_Coupled_Index] += rho_AD[icv] * CmeanSource_AD[icv] * cv_volume[icv];

			if (flagImplicit && icv<ncv) {
				// ********************* Derivative computation ***********************

				double dCmeanSource_dZM, dCmeanSource_dZV, dCmeanSource_dCM,E;
				double delta_ZM = 1.0e-5;
				double delta_ZV = 1.0e-8;
				double delta_CM = 1.0e-5;
				double pressp, pressm, CmeanSourcep, CmeanSourcem;
				double ZV = (*ZVar)[icv].value();
				double ZM = (*ZMean)[icv].value();
				double CM = (*CMean)[icv].value();

				E = enthalpy[icv].value() - RoM[icv].value() * temp[icv].value();

				// Derivative with respect to ZMean
				double ZMp = (*ZMean)[icv].value() + delta_ZM;
				double ZMm = (*ZMean)[icv].value() - delta_ZM;

				pressure_scalSource_AD(pressp, CmeanSourcep, rho[icv], E, ZMp, ZV, CM);
				pressure_scalSource_AD(pressm, CmeanSourcem, rho[icv], E, ZMm, ZV, CM);
#ifdef PRESSURE_SCALING_FPVA_COEFF
				CmeanSourcep *= pow( pressp/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
				CmeanSourcem *= pow( pressm/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
#endif
				dCmeanSource_dZM = (CmeanSourcep - CmeanSourcem) / (2.0 * delta_ZM);

				// Derivative with respect to ZVar
				double ZVp = (*ZVar)[icv].value() + delta_ZV;
				double ZVm = (*ZVar)[icv].value() - delta_ZV;
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
				A[noc00][0][5+CMean_Coupled_Index] -= (CmeanSource_AD[icv].value() - (*ZMean)[icv].value() * dCmeanSource_dZM - (*ZVar)[icv].value() * dCmeanSource_dZV - (*CMean)[icv].value() * dCmeanSource_dCM) * cv_volume[icv];
			}
		}
	}

	/*
	 * Method: pressure_scalSource_AD
	 * ------------------------------
	 *
	 */
	void pressure_scalSource_AD(double &pp, double &CMSource, double &rrho, double &E, double &ZM, double &ZV, double &CM) {
		double r0, T, t0, e0, gamma0, agamma;
		myChemTable.LookupSelectedCoeff(r0, t0, e0, gamma0, agamma, ZM, ZV, CM);
		CMSource = myChemTable.Lookup(ZM, ZV, CM, "SRC_PROG");
		if (agamma == 0.0)
			T = t0 + (gamma0 - 1.0) / r0 * (E - e0);
		else
			T = t0 + (gamma0 - 1.0) / agamma * (exp(agamma * (E - e0) / r0) - 1.0);

		// Temperature clipping
		T  = max(T, Tminimum);
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


#endif  /* IKECOMBMODEL_FPVA_COEFF_AD_H */
