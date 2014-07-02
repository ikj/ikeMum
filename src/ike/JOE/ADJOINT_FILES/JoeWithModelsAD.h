#ifndef JOEWITHMODELSAD_H
#define JOEWITHMODELSAD_H

#include "MiscUtils.h"
#include "tc_vec3d.h"

#include "UgpWithCvCompFlowAD.h"
#include "JoeWithModels.h"
#include "adolc.h"


class JoeWithModels_AD: virtual public JoeWithModels, virtual public UgpWithCvCompFlow_AD
{
private:

	int icv_debug;

	void setDebugCv() {

		icv_debug = -1;
		double x_debug[3] = {46.2617,0.5,0.376642};

		for (int icv = 0; icv < ncv_ggff; ++icv) {
			double d2 = 0.0;
			FOR_I3 {
				double dx = x_cv[icv][i] - x_debug[i];
				d2 += dx*dx;
			}
			if (sqrt(d2) < 1.0E-4) {
				cout << "found it on rank: " << mpi_rank << " " << icv<< " "<<x_cv[icv][0] << " " << x_cv[icv][1] << " " << x_cv[icv][2];
				if (icv < ncv)
					cout << " internal cv" << endl;
				else if (icv < ncv_g)
					cout << " ghost level 1 cv" << endl;
				else if (icv < ncv_gg)
					cout << " ghost level 2 cv" << endl;
				else if (icv < ncv_ggf)
					cout << " fake level 1 cv" << endl;
				else
					cout << " fake level 2 cv" << endl;
				icv_debug = icv;
			}
		}
		MPI_Pause("did we find one?");

	}

	void dumpDebugState(char * message) {

		int nScal =  scalarTranspEqVector.size();
		if (icv_debug >= 0) {
			cout << "On processor "<<mpi_rank<<" "<<"debug state: " << rho_AD[icv_debug] << " " << rhou_AD[icv_debug][0] << " " << rhou_AD[icv_debug][1] << " " << rhoE_AD[icv_debug] << endl;
			for(int iScal=0; iScal<nScal; iScal++)
				cout<<"Scalar "<<iScal<<" "<<scalarTranspEqVector_AD[iScal].phi[icv_debug]<<endl;
		}
		MPI_Pause(message);
	}


public:
	/**
	 * constructor, pass ParamMap
	 */
	JoeWithModels_AD(ParamMap &p) : JoeWithModels(p), UgpWithCvCompFlow_AD(p) {
		if (mpi_rank == 0)
			cout<<"JoeWithModels_AD()"<<endl;
		init();
	}

	/**
	 * constructor, pass name of joe's input file
	 */
	JoeWithModels_AD(char *name) : JoeWithModels(name), UgpWithCvCompFlow_AD(name) {
		if (mpi_rank == 0)
			cout<<"JoeWithModels_AD()"<<endl;
		init();
	}

	virtual ~JoeWithModels_AD() {
		if (mpi_rank == 0)
			cout<<"~JoeWithModels_AD()"<<endl;
	}

	virtual void init() {
		if (mpi_rank == 0)
			cout << "JoeWithModels_AD()::init()"<< endl;
	}

public:

	void runAdjoint();

	void calcFunctionalDerivative();

	virtual void initialHook_AD() {/* empty */};

#ifdef USE_MEM_SAVING_ADVAR
	virtual void calcFunctional_AD(ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<REALQ> &rhoE) {/* empty */};
#endif
	virtual void calcFunctional_AD(REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE) {/* empty */};

	void calcResidualDerivative(double (*A)[5][5], double ***AScal, int flagImplicit);

#ifdef USE_MEM_SAVING_ADVAR
	int calcResidual_AD(REALA *rhs_rho_AD, REALA (*rhs_rhou_AD)[3],REALA *rhs_rhoE_AD, REALAS **rhs_rhoScal_AD,
			ADscalar<REALQ> &rho_AD, ADvector<REALQ> &rhou_AD, ADscalar<REALQ> &rhoE_AD, double (*A)[5][5], double ***AScal, int flagImplicit);
#endif
	int calcResidual_AD(REALA *rhs_rho_AD, REALA (*rhs_rhou_AD)[3], REALA *rhs_rhoE_AD, REALAS **rhs_rhoScal_AD,
			REALQ *rho_AD, REALQ (*rhou_AD)[3], REALQ *rhoE_AD, double (*A)[5][5], double ***AScal, int flagImplicit);

#ifdef USE_MEM_SAVING_ADVAR
	virtual void boundaryHook_AD(ADscalar<REALQ> &T_fa, ADvector<REALQ> &vel_fa, ADscalar<REALQ> &p_fa, FaZone *zone) {
		if(mpi_rank==0)
			cout<<"WARNING in JoeWithModels_AD::boundaryHook_AD(): empty function"<<endl;
	}
	virtual void sourceHook_AD(ADscalar<REALA> &rhs_rho, ADvector<REALA> &rhs_rhou, ADscalar<REALA> &rhs_rhoE, double (*A)[5][5]) {
		if(mpi_rank==0)
			cout<<"WARNING in JoeWithModels_AD::sourceHook_AD(): empty function"<<endl;
	}
	virtual void sourceHookRansTurb_AD(ADscalar<REALA> &rhs_rho, ADvector<REALA> &rhs_rhou, ADscalar<REALA> &rhs_rhoE, double (*A)[5][5]) {
		if(mpi_rank==0)
			cout<<"WARNING in JoeWithModels_AD::sourceHookRansTurb_AD(): empty function"<<endl;
	}
	virtual void sourceHookRansComb_AD(ADscalar<REALA> &rhs_rho, ADvector<REALA> &rhs_rhou, ADscalar<REALA> &rhs_rhoE, double (*A)[5][5]) {
		if(mpi_rank==0)
			cout<<"WARNING in JoeWithModels_AD::sourceHookRansComb_AD(): empty function"<<endl;
	}
#endif
	virtual void boundaryHook_AD(REALQ *T_fa, REALQ (*vel_fa)[3], REALQ *p_fa, FaZone *zone) {/* empty */}
	virtual void sourceHook_AD(REALA *rhs_rho, REALA (*rhs_rhou)[3], REALA *rhs_rhoE, double (*A)[5][5]) {/* empty */}
	virtual void sourceHookRansTurb_AD(REALA *rhs_rho, REALA (*rhs_rhou)[3], REALA *rhs_rhoE, double (*A)[5][5]) {/* empty */}
	virtual void sourceHookRansComb_AD(REALA *rhs_rho, REALA (*rhs_rhou)[3], REALA *rhs_rhoE, double (*A)[5][5]) {/* empty */}
	//virtual void sourceHookCoupled_AD(REALA **rhs,          double ***A,  int flagImplicit) {/* empty */}
	//virtual void sourceHookRansTurbCoupled_AD(REALA **rhs,          double ***A,  int flagImplicit) {/* empty */}
	//virtual void sourceHookRansCombCoupled_AD(REALA **rhs,          double ***A,  int flagImplicit) {/* empty */}

#ifdef USE_MEM_SAVING_ADVAR
	int calcRhs_AD(REALA *rhs_rho_AD, REALA (*rhs_rhou_AD)[3], REALA *rhs_rhoE_AD, REALAS **rhs_rhoScal_AD,
			ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<REALQ> &rhoE, double (*A)[5][5], double ***AScal, int flagImplicit);

	int calcEulerFlux_AD(REALA *rhs_rho, REALA (*rhs_rhou)[3], REALA *rhs_rhoE, REALAS **rhs_rhoScal, ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<REALQ> &rhoE, double (*A)[5][5], double ***AScal, int flagImplicit);

	void calcViscousFluxNS_AD(REALA *rhs_rho, REALA (*rhs_rhou)[3], REALA *rhs_rhoE, ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<REALQ> &rhoE, double (*A)[5][5], int flagImplicit);

	void setBC_AD(ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<adouble> &rhoE);
#endif
	int calcRhs_AD(REALA *rhs_rho_AD, REALA (*rhs_rhou_AD)[3], REALA *rhs_rhoE_AD, REALAS **rhs_rhoScal_AD,
				REALQ *rho_AD, REALQ (*rhou_AD)[3], REALQ *rhoE_AD, double (*A)[5][5], double ***AScal, int flagImplicit);

	int calcEulerFlux_AD(REALA *rhs_rho, REALA (*rhs_rhou)[3], REALA *rhs_rhoE, REALAS **rhs_rhoScal, REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE, double (*A)[5][5], double ***AScal, int flagImplicit);

	void calcViscousFluxNS_AD(REALA *rhs_rho, REALA (*rhs_rhou)[3], REALA *rhs_rhoE, REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE, double (*A)[5][5], int flagImplicit);

	void setBC_AD(REALQ *rho_AD, REALQ (*rhou_AD)[3], REALQ *rhoE_AD);

	void solveADSystemExplicitEuler();

	void solveADSystemExplicitRK();

	void solveADSystemImplicitEuler();

	void solveADSystem_Sparse();
	void solveADSystem_Sparse1();


	virtual void calcFunctionalGradient(double *adj_vars) {/* empty */};

	virtual void print_tapestats(int tag);

	void showResidue(double *rhsResid, int step);

	void solveADSystemImplicitEuler_Coupled();

	void calcResidualDerivative_Coupled(double ***A, int flagImplicit);

#ifdef USE_MEM_SAVING_ADVAR
	void calcResidual_AD_Coupled(REALA **rhs_AD, ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<REALQ> &rhoE, double ***A, int flagImplicit);

	void calcRhs_AD_Coupled(REALA **rhs, ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<REALQ> &rhoE, double ***A, int flagImplicit);

	void calcFluxCoupled_AD(REALA **rhs, ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<REALQ> &rhoE, double ***A, int nScal, int flagImplicit);
#endif
	void calcResidual_AD_Coupled(REALA **rhs_AD, REALQ *rho_AD, REALQ (*rhou_AD)[3], REALQ *rhoE_AD, double ***A, int flagImplicit);

	void calcRhs_AD_Coupled(REALA **rhs, REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE, double ***A, int flagImplicit);

	void calcFluxCoupled_AD(REALA **rhs, REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE, double ***A, int nScal, int flagImplicit);

	virtual void finalHook_AD() {/*empty */}

};


#endif



