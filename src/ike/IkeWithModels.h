/*
 * IkeWithModels.h
 *
 *  Created on: Jan 4, 2013
 *      Author: ikj
 */

#ifndef IKEWITHMODELS_H_
#define IKEWITHMODELS_H_

#include "JOE/ADJOINT_FILES/JoeWithModelsAD.h"
#include <adolc_sparse.h>
#include <set>
//#include "MatComprsed.h"
//#include "PetscSolver2.h"

#include "IkeUtilsAD.h"

#ifndef WITH_PETSC
#define WITH_PETSC
#endif

#ifndef ABSURDLY_BIG_NUMBER
#define ABSURDLY_BIG_NUMBER 2.22e22
#endif

#ifndef MACHINE_EPS
#define MACHINE_EPS 2.0e-16 // Note: machine epsilon for double precision = 1.11e-16
#endif

#define IKEWITHMODELS_ERROR_CODE 0

//#define USE_DT_OVER_VOL_SCALING
#define DT_OVER_VOL_MAX 1.0e20

#ifdef USE_MEM_SAVING_ADVAR
#define USE_MEM_SAVING_ADVAR_1D_ // Note: Actually, USE_MEM_SAVING_ADVAR is defined in UgpWithCvCompFlowAD.h
#ifdef USE_MEM_SAVING_ADVAR_1D_
	#include "ADvar.h"
#endif
#endif

/*
 * Bug report:
 *   In JoeWithModelsAD.cpp (Don't forget to fix JoeWithModels.cpp in the same way)
 *     1. cannot detect "kine_index" in the calcViscousFluxNS_AD() method:
 *        Original: if (scalarTranspEqVector[iScal].getName() == "kine")
 *        Modified: if (strcmp(scalarTranspEqVector[iScal].getName(),"kine")==0)
 *
 *     2. "kine_fa" is wrong in the calcViscousFluxNS_AD() method:
 *        Original: if (kine_index > -1)
 *			        {
 *	                   kine0 = scalarTranspEqVector_AD[kine_index].phi[icv0];
 *                     kine1 = scalarTranspEqVector_AD[kine_index].phi[icv1];
 *                     REALQS kine_fa = w1*kine0 + w0*kine1;
 *                  }
 *        Modified: if (kine_index > -1)
 *			        {
 *	                   kine0 = scalarTranspEqVector_AD[kine_index].phi[icv0];
 *                     kine1 = scalarTranspEqVector_AD[kine_index].phi[icv1];
 *                     kine_fa = w1*kine0 + w0*kine1;
 *                  }
 *   In UgpWithCvCompFlowAD.h
 *     1. The member variables of the ScalarTranspEq_AD class are not deleted when the class is destroyed:
 *        Make the class delete its member variables in its destructor if they are not equal to NULL
 *        (Do the same thing for the ScalarTranspEq class in UgpWithCv2.h
 *         -- However, the destructor of the UgpWithCvCompFlow class in UgpWithCvCompFlow.h tries to delete
 *            the "diff" array of the ScalarTranspEq_AD class.
 *            Thus, comment out that part from the UgpWithCvCompFlow class)
 *   In UgpWithTools.h
 *     1. The member variable "geom" of the WriteData class is not freed: free it in the destructor if it is not equal to NULL
 *        (??? -> This makes an segmentation fault error)
 *        Modified: virtual ~WriteData() {
 *                     if(geom != NULL) {
 *                        delete geom; geom = NULL;
 *                     }
 *                  }
 *   In UgpWithCvCompFlow.h
 *     1. The size of the member variables mul_fa, lamOcp_fa, and mut_fa must be nfa_b2gg.
 *        Thus, in the initializeFromRestartFile() function,
 *          mul_fa    = new double[nfa_b2gg];
 *          lamOcp_fa = new double[nfa_b2gg];
 *          mut_fa    = new double[nfa_b2gg];
 */

/*
 * Modifications in JOE and ADKOINT_FILES:
 *   UgpWithCvCompFlowAD.h:
 *     1. insert the following two lines at the beginning of UgpWithCvCompFlowAD.h
 *       enum BOUNDARY_TYPE {WRONG, INTERNAL, HOOK, CBC, CBC_SUBSONIC_INLET, CBC_SUBSONIC_OUTLET, SYMMETRY, NEUMANN, WALL};
 *       enum BOUNDARY_TYPE_SCALAR {INTERNAL_SCALAR, HOOK_SCALAR, DIRICHLET_SCALAR, FLUX_SCALAR, OTHER_SCALAR};
 *     2. insert the following line in the UgpWithCvCompFlowAD class as public
 *       #include "../../IkeUgpWithCvCompFlow.h"
 *
 *   JoeWithModelsAD.h and JoeWithModelsAD.cpp:
 *   (This is for the backtracking algorithm in the pseudo-arclength code)
 *     1. make the return type of calcRhs_AD(), calcEulerFlux_AD(), and calcResidual_AD() as integer:
 *       int calcRhs_AD(REALA *rhs_rho_AD, ...)
 *       int calcEulerFlux_AD(REALA *rhs_rho, ...)
 *       int calcResidual_AD(REALA *rhs_rho_AD, ...)
 *     2. make the two function returns the number of reduced order
 *       2-1. In the calcEulerFlux_AD() function,
 *         int JoeWithModels_AD::calcEulerFlux_AD(REALA *rhs_rho, ...) {
 *             ....
 *             Return CountReducedOrder;
 *         }
 *       2-2. In the calcRhs_AD() function,
 *         int JoeWithModels_AD::calcRhs_AD(REALA *rhs_rho, ...) {
 *             ...
 *             int CountReducedOrder = calcEulerFlux_AD(rhs_rho, ...);
 *             ...
 *             Return CountReducedOrder;
 *         }
 *       3-3. In the calcResidual_AD() function,
 *         int JoeWithModels_AD::calcResidual_AD(REALA *rhs_rho_AD, ...) {
 *             ...
 *             int CountReducedOrder = calcRhs_AD(rhs_rho_AD, ...)
 *             ...
 *             Return CountReducedOrder;
 *         }
 *
 *   JoeWithModels.h and JoeWithModels.cpp:
 *   (This is for the backtracking algorithm in the pseudo-arclength code):
 *     1. make the return type of calcRhs(),
 *       int calcRhs(double *rhs_rho, ...)
 *       int calcEulerFlux(double *rhs_rho, ...)
 *     2. make the two function returns the number of reduced order
 *       4.1. In the calcRhs() function,
 *         int JoeWithModels::calcRhs(double *rhs_rho, ...) {
 *             ....
 *             int CountReducedOrder = calcEulerFlux(rhs_rho, rhs_rhou, rhs_rhoE, rhs_rhoScal, A, AScal, flagImplicit);
 *             ....
 *             return CountReducedOrder;
 *         }
 *       4.2. In the calcEulerFlux() function,
 *         int JoeWithModels::calcEulerFlux(double *rhs_rho, ...) {
 *             ...
 *             return CountReducedOrder;
 *         }
 *
 *     3. make the return type of calcRhsCoupled(), calcFluxCoupled(), and calcRhsCoupled() as integer:
 *       int calcFluxCoupled(double **rhs, ...)
 *       int calcRhsCoupled(double **rhs, ...)
 *     4. make the two function returns the number of reduced order
 *       4-1. In the calcRhsCoupled() function,
 *         int JoeWithModels::calcRhsCoupled(double **rhs, ...) {
 *             ....
 *             int CountReducedOrder = calcFluxCoupled(rhs, A, nScal, flagImplicit);
 *             ....
 *             return CountReducedOrder
 *         }
 *       4-2. In the calcFluxCoupled() function,
 *         int JoeWithModel::calcFluxCoupled(double **rhs, ...) {
 *             ....
 *             Return CountReducedOrder;
 *         }
 *
 *   And there are more...
 */

/* Note: Numbering convention in Joe (the version of two-layer-ghosts)
 *   1. Faces
 *     0 ~ nfa_b-1 : boundary faces
 *     nfa_b ~ nfa_bpi-1 :
 *     nfa_bpi ~ nfa-1 :
 *     nfa ~ nfa_b2 : additional boundary faces between g1:f2 cvs
 *     nfa_b2 ~ nfa_b2g-1: faces between g1:g1 cvs
 *     nfa_b2g-1 ~ nfa_b2gg-1: faces between g1:g2 cvs
 *   2. CVs
 *     0 ~ ncv-1: internal cells
 *     ncv ~ ncv_g-1 : includes level-1 ghosts
 *     ncv_g ~ ncv_gg-1 : includes level-2 ghosts (face-based nbrs of level-1 ghosts)
 *     ncv_gg ~ ncv_ggf-1 : includes fakes owned by nfa_b faces: i.e. ncv_ggf-ncv_gg == nfa_b
 *     ncv_ggf ~ ncv_ggff-1 : includes fakes that are nbrs of level-1 ghosts
 */

/*
 * Class: IkeWithModels_AD
 */
class IkeWithModels_AD: virtual public JoeWithModels_AD //, virtual public JoeWithModels, virtual public UgpWithCvCompFlow_AD
{
public:
	/*
	 * constructor, pass ParamMap
	 */
	IkeWithModels_AD(ParamMap &p) : JoeWithModels_AD(p), JoeWithModels(p), UgpWithCvCompFlow_AD(p) {
		if(mpi_rank==0)
			cout<<"IkeWithModels_AD()"<<endl;
		this->init();
	}
	/*
	 * constructor, pass name of joe's input file
	 */
	IkeWithModels_AD(char *name) : JoeWithModels_AD(name), JoeWithModels(name), UgpWithCvCompFlow_AD(name) {
		if(mpi_rank==0)
			cout<<"IkeWithModels_AD()"<<endl;
		this->init();
	}

	/*
	 * destructor
	 */
	virtual ~IkeWithModels_AD() {
		if(mpi_rank==0)
			cout<<"~IkeWithModels_AD()"<<endl;
		this->clear();
	}

public:
	// Two-layerd CSR struct
	int *nbocv2_i;
	vector<int> nbocv2_v;

	// Global icv index
	int *nbocv2_v_global;

	// Boundary condition
//	BOUNDARY_TYPE* btofa;
//	BOUNDARY_TYPE_SCALAR** btofaScalar;
//	double** bvalScalar;
	FaZone** zofa;

	vector<std::string> warningsVectorIkeModels1D; // If the simulation breaks down for some reasons,
			                                       // it will dump out a massive amount of warning/error messages.
			                                       // Thus, this object stores the (only) warning messages so that
			                                       // the user can manipulate the amount of messages that should be
			                                       // printed out on the screen.

#ifdef USE_DT_OVER_VOL_SCALING
	// Local dt/dV: u^{n+1}-u{n} = RHS * (dt/dV)
	double *local_dtOverVol;
#endif

private:
	/*
	 * Method: init
	 * ------------
	 * Initialize member variables
	 */
	void init();

	/*
	 * Method: clear
	 * -------------
	 * Clear member variables
	 */
	void clear();

public:
#ifdef USE_MEM_SAVING_ADVAR_1D_
	/*
	 * Method: boundaryHook1D_AD
	 * -------------------------
	 * Original code = boundaryHook_AD() in JoeWithModelsAD.h
	 */
	virtual void boundaryHook1D_AD(const int ifa, ADscalar<REALQ> &T_fa, ADvector<REALQ> &vel_fa, ADscalar<REALQ> &p_fa, FaZone *zone) {/* empty */}

	/*
	 * Method: boundaryHook1D_AD
	 * -------------------------
	 * Original code = boundaryHook() in JoeWithModels.h
	 * In some situation, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
	 */
	virtual void boundaryHook1D(const int ifa, ADscalar<REALQ> &T_fa, ADvector<REALQ> &vel_fa, ADscalar<REALQ> &p_fa, FaZone *zone) {/* empty */}
#endif

	/*
	 * Method: boundaryHook1D_AD
	 * -------------------------
	 * Original code = boundaryHook_AD() in JoeWithModelsAD.h
	 */
	virtual void boundaryHook1D_AD(const int ifa, REALQ *T_fa, REALQ (*vel_fa)[3], REALQ *p_fa, FaZone *zone) {/* empty */}

	/*
	 * Method: boundaryHook1D_AD
	 * -------------------------
	 * Original code = boundaryHook() in JoeWithModels.h
	 * In some situation, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
	 */
	virtual void boundaryHook1D(const int ifa, double *T_fa, double (*vel_fa)[3], double *p_fa, FaZone *zone) {/* empty */}

	/*
	 * Method: sourceHook1D_AD
	 * -----------------------
	 * Original code = sourceHook_AD() in JoeWithModelsAD.h
	 */
	virtual void sourceHook1D_AD(const int icvCenter, REALA &rhs_rho, REALA rhs_rhou[3], REALA &rhs_rhoE, double (*A)[5][5]) {/* empty */}

	/*
	 * Method: sourceHookRansTurb1D_AD
	 * -------------------------------
	 * Original code = sourceHookRansTurb_AD() in JoeWithModelsAD.h
	 */
	virtual void sourceHookRansTurb1D_AD(const int icvCenter, REALA &rhs_rho, REALA rhs_rhou[3], REALA &rhs_rhoE, double (*A)[5][5]) {/* empty */}

	/*
	 * Method: sourceHookRansComb1D_AD
	 * -------------------------------
	 * Original code = sourceHookRansComb_AD() in JoeWithModelsAD.h
	 */
	virtual void sourceHookRansComb1D_AD(const int icvCenter, REALA &rhs_rho, REALA rhs_rhou[3], REALA &rhs_rhoE, double (*A)[5][5]) {/* empty */}

protected:
	/*
	 * Method: getReferenceParams
	 * --------------------------
	 * Get reference parameters from the input file
	 * Note: If the reference value for a variable (except uRefVec and turbScalRefVec) cannot be found,
	 *       1.0 will be stored (The value 1.0 is chosen for normalization).
	 */
	void getReferenceParams(REF_FLOW_PARAMS& RefFlowParams);

	/*
	 * Method: build2layerCSRstruct
	 * ----------------------------
	 * Build a bigger CSR structure which contains two-layered neighbors of each CV:
	 *   The results are nbocv2_i = same as nbocv_i but 2-layered
	 *                   nbocv2_v = same as nbocv_v but 2-layered
	 */
	void build2layerCSRstruct(bool addFakeCells);

	/*
	 * Method: find2layerCSRstruct_eachIcv
	 * -----------------------------------
	 * Find all the neighbors in the two-layer of a given CV
	 */
	int find2layerCSRstruct_eachIcv(int icvCenter, vector<int> &nbocv2_v, bool addFakeCells);

//	/*
//	 * Method: setBtofa
//	 * ----------------
//	 * Set up the boundary type of faces (btofa, btofaScalar, bvalScalar)
//	 * The variables are first initialized as internal faces. Then, by iterating zones, a proper BC is applied to each face.
//	 * i.e. set up btofa as one of BOUNDARY_TYPEs
//	 *      set up btofaScalar as one of BOUNDARY_TYPE_SCALARs
//	 *      store given scalar boundary value in bvalScalar
//	 */
//	void setBtofa();

	/*
	 * Method: setZofa
	 * ---------------
	 * Set up the zone type of faces (zofa): zofa is an array which stores the pointers to a FaZone at each face.
	 */
	void setZofa();

	/*
	 * Method: find1LayerFaceIndices
	 * -----------------------------
	 * Find all the faces in the one-layer neighbors.
	 * For example, 4 faces from the total 50 faces will be found for the following 2D structure,
	 *            -----
	 *            |   |
	 *        -------------
	 *        |   |   |   |
	 *    ---------------------
	 *    |   |   | x |   |   |
	 *    ---------------------
	 *        |   |   |   |
	 *        -------------
	 *            |   |
	 *            -----
	 * Return: total number of faces
	 */
	int find1LayerFaceIndices(const int icvCenter, vector<int>& faInternal, vector<int>& faBoundary);
	int find1LayerFaceIndices(const int icvCenter, vector<int>& fa2);

	/*
	 * Method: find2LayerFaceIndices
	 * -----------------------------
	 * Find all the faces in the two-layer neighbors.
	 * For example, 16 faces from the total 50 faces will be found for the following 2D structure,
	 *            -----
	 *            |   |
	 *        -------------
	 *        |   |   |   |
	 *    ---------------------
	 *    |   |   | x |   |   |
	 *    ---------------------
	 *        |   |   |   |
	 *        -------------
	 *            |   |
	 *            -----
	 * Return: total number of faces
	 */
	int find2LayerFaceIndices(const int icvCenter, vector<int>& faInternal, vector<int>& faBoundary);
	int find2LayerFaceIndices(const int icvCenter, vector<int>& fa2);

#ifdef USE_DT_OVER_VOL_SCALING
	/*
	 * Method: calcDtOverVol
	 * ---------------------
	 * Calculate dt/cv_vol: u^{n+1}-u{n} = RHS * (dt/cv_vol)
	 * "local_dtOverVol" will be used to scale RHS
	 * Original code: UgpWithCvCompFlow::calcDt()
	 */
	void calcDtOverVol(double &dtOverVol_min, double &dtOverVol_max, const double cfl_target);
#endif

#ifdef USE_MEM_SAVING_ADVAR_1D_
	/*
	 * Method: calcResidual1D_AD
	 * -------------------------
	 * Original code = calcResidual1D_AD() in JoeWithModelsAD.cpp
	 * Note: gamma and RoM must have been updated before calling this method
	 *       calcResidual1D_AD() successively calls the following methods:
	 *         1. calcStateVariables1D_AD()                                     -- for the 2-layer neighboring cells
	 *         2. calcMaterialProperties1D_AD()                                 -- for the 1-layer neighboring faces
	 *         3. setBC1D_AD()                                                  -- for the 1-layer neighboring faces
	 *         4. if(mu_ref>0.0 || sndOrder==true): calcCv2Grad1D_AD(grad_u)    -- for the 1-layer neighboring cells
	 *         5. if(mu_ref>0.0)                  : calcRansTurbViscMuet1D_AD() -- for the 1-layer neighboring faces (virtual function!)
	 *         6. calcRhs1D_AD()
	 */
	int calcResidual1D_AD(const int icvCenter,
			REALA &rhs_rho_AD, REALA rhs_rhou_AD[3],REALA &rhs_rhoE_AD, REALAS *rhs_rhoScal_AD,
			ADscalar<REALQ> &rho_AD, ADvector<REALQ> &rhou_AD, ADscalar<REALQ> &rhoE_AD, double (*A)[5][5], double ***AScal, int flagImplicit);
	/*
	 * Method: calcStateVariables1D_AD
	 * -------------------------------
	 * Calculate state variables (vel, press, temp, enthalpy, sos) for the neighbors of icvCenter
	 * Original code = calcStateVariables_AD() in UgpWithCvCompFlowAD.h
	 */
	void calcStateVariables1D_AD(const int icvCenter, ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<REALQ> &rhoE);

	/*
	 * Method: calcMaterialProperties1D_AD
	 * -----------------------------------
	 * Calculate material properties (mul_fa, lamOcp_fa) only for the 1-layer neighboring faces of icvCenter
	 * Original code: calcMaterialProperties1D_AD in UgpWithCvCompFlowAD.h
	 */
	void calcMaterialProperties1D_AD(const int icvCenter, vector<int>& faInternal, vector<int>& faBoundary, ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<REALQ> &rhoE);

	/*
	 * Method: setBC1D_AD
	 * ------------------
	 * Original code = setBC_AD in JoeWithModelsAD.cpp
	 */
	void setBC1D_AD(const int icvCenter, vector<int>& faInternal, vector<int>& faBoundary, ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<REALQ> &rhoE);
#endif

	/*
	 * Method: calcResidual1D_AD
	 * -------------------------
	 * Original code = calcResidual1D_AD() in JoeWithModelsAD.cpp
	 * Note: gamma and RoM must have been updated before calling this method
	 *       calcResidual1D_AD() successively calls the following methods:
	 *         1. calcStateVariables1D_AD()                                     -- for the 2-layer neighboring cells
	 *         2. calcMaterialProperties1D_AD()                                 -- for the 1-layer neighboring faces
	 *         3. setBC1D_AD()                                                  -- for the 1-layer neighboring faces
	 *         4. if(mu_ref>0.0 || sndOrder==true): calcCv2Grad1D_AD(grad_u)    -- for the 1-layer neighboring cells
	 *         5. if(mu_ref>0.0)                  : calcRansTurbViscMuet1D_AD() -- for the 1-layer neighboring faces (virtual function!)
	 *         6. calcRhs1D_AD()
	 */
	int calcResidual1D_AD(const int icvCenter,
			REALA &rhs_rho_AD, REALA rhs_rhou_AD[3], REALA &rhs_rhoE_AD, REALAS *rhs_rhoScal_AD,
			REALQ *rho_AD, REALQ (*rhou_AD)[3], REALQ *rhoE_AD, double (*A)[5][5], double ***AScal, int flagImplicit);

	/*
	 * Method: calcStateVariables1D_AD
	 * -------------------------------
	 * Calculate state variables (vel, press, temp, enthalpy, sos) only for the 2-layer neighbors of icvCenter
	 * Original code = calcStateVariables_AD() in UgpWithCvCompFlowAD.h
	 */
	void calcStateVariables1D_AD(const int icvCenter, REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE);

	/*
	 * Method: calcMaterialProperties1D_AD
	 * -----------------------------------
	 * Calculate material properties (mul_fa, lamOcp_fa) only for the 1-layer neighboring faces of icvCenter
	 * Original code: calcMaterialProperties1D_AD in UgpWithCvCompFlowAD.h
	 */
	void calcMaterialProperties1D_AD(const int icvCenter, vector<int>& faInternal, vector<int>& faBoundary, REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE);

	/*
	 * Method: setBC1D_AD
	 * ------------------
	 * Original code = setBC_AD in JoeWithModelsAD.cpp
	 */
	void setBC1D_AD(const int icvCenter, vector<int>& faInternal, vector<int>& faBoundary, REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE);

	/*
	 * Method: setBC1D
	 * ---------------
	 * Original code = setBC1D_AD in IkeWithModelsAD.cpp and setBC() in JoeWithModels.cpp
	 * In some situation, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
	 */
	void setBC1D(const int ifa, double *rho, double (*rhou)[3], double *rhoE);

	/*
	 * Method: setScalarBC1D_AD
	 * ------------------------
	 * Original code = setScalarBC_AD() in ScalarsAD.cpp
	 */
	void setScalarBC1D_AD(const int ifa);

	/*
	 * Method: setScalarBC1D
	 * ------------------------
	 * Original code = setScalarBC1D_AD() in IkeWithModels.cpp and setBC() in Scalars.cpp
	 * In some situation, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
	 */
	void setScalarBC1D(const int ifa);

	/*
	 * Method: ComputeBCProperties1D_T_AD
	 * ----------------------------------
	 * brief Compute for a given temperature the properties of the mixture at a face.
	 * Original code: ComputeBCProperties_T_AD() in UgpWithCvCompFlowAD.h
	 */
	void ComputeBCProperties1D_T_AD(const int ifa);

	/*
	 * Method: ComputeBCProperties1D_T
	 * -------------------------------
	 * brief Compute for a given temperature the properties of the mixture at a face.
	 * Original code: ComputeBCProperties1D_T_AD() in IkeWithModel.cpp and ComputeBCProperties_T() UgpWithCvCompFlow.h
	 * In some situation, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
	 */
	void ComputeBCProperties1D_T(const int ifa);

	/*
	 * Method: ComputeBCProperties1D_H_AD
	 * ----------------------------------
	 * brief Compute for a given enthalpy the properties of the mixture at a face.
	 * Original code: ComputeBCProperties1D_H_AD() in UgpWithCvCompFlowAD.h
	 */
	void ComputeBCProperties1D_H_AD(const int ifa);

	/*
	 * Method: ComputeBCProperties1D_H
	 * -------------------------------
	 * brief Compute for a given enthalpy the properties of the mixture at a face.
	 * Original code: ComputeBCProperties1D_H_AD() in UgpWithCvCompFlowAD.h
	 */
	void ComputeBCProperties1D_H(const int ifa);

#ifdef USE_MEM_SAVING_ADVAR_1D_
	/*
	 * Method: calcRhs1D_AD
	 * --------------------
	 * Calculate the Rhs of the N-S equations in the so-called "1D-style"
	 * Original code = calcRhs_AD() in JoeWithModelsAD.cpp
	 */
	int calcRhs1D_AD(int icvCenter, REALA &rhs_rho, REALA rhs_rhou[3], REALA &rhs_rhoE, REALAS *rhs_rhoScal,
			ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<REALQ> &rhoE,
			double (*A)[5][5], double ***AScal, int flagImplicit);

	/*
	 * Method: calcEulerFlux1D_AD
	 * --------------------------
	 * Calculate Euler flux in the so-called "1D style"
	 * Return: the number of times switched back to first order at faces
	 * Original code = calcEulerFlux_AD() in JoeWithModelsAD.cpp
	 */
	int calcEulerFlux1D_AD(const int icvCenter, REALA &rhs_rho, REALA rhs_rhou[3], REALA &rhs_rhoE, REALAS *rhs_rhoScal,
			ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<REALQ> &rhoE,
			double (*A)[5][5], double ***AScal, int flagImplicit);

	/*
	 * Method: calcViscousFluxNS1D_AD
	 * ------------------------------
	 * Calculate viscous flux in the so-called "1D style"
	 * Original code = calcViscousFluxNS_AD() in JoeWithModelsAD.cpp
	 */
	void calcViscousFluxNS1D_AD(const int icvCenter, REALA &rhs_rho, REALA rhs_rhou[3], REALA &rhs_rhoE,
			ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<REALQ> &rhoE,
			double (*A)[5][5], int flagImplicit);

//	/*
//	 * Method: calcViscousFluxScalar_new1D_AD
//	 * --------------------------------------
//	 * Original code: calcViscousFluxScalar_new_AD in ScalarsAD.cpp
//	 */
//	void calcViscousFluxScalar_new1D_AD(const int icvCenter, const int iScal, adouble &rhs_rhoScal, double *AScal, ScalarTranspEq &transpScal, ScalarTranspEq_AD &transpScal_AD, int flagImplicit);
#endif

	/*
	 * Method: calcRhs1D_AD
	 * --------------------
	 * Calculate the Rhs of the N-S equations in the so-called "1D-style"
	 * Original code = calcRhs_AD() in JoeWithModelsAD.cpp
	 */
	int calcRhs1D_AD(int icvCenter, REALA &rhs_rho, REALA rhs_rhou[3], REALA &rhs_rhoE, REALAS *rhs_rhoScal,
			REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE,
			double (*A)[5][5], double ***AScal, int flagImplicit);

	/*
	 * Method: calcEulerFlux1D_AD
	 * --------------------------
	 * Calculate Euler flux in the so-called "1D style"
	 * Return: the number of times switched back to first order at faces
	 * Original code = calcEulerFlux_AD() in JoeWithModelsAD.cpp
	 */
	int calcEulerFlux1D_AD(const int icvCenter, REALA &rhs_rho, REALA rhs_rhou[3], REALA &rhs_rhoE, REALAS *rhs_rhoScal,
			REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE,
			double (*A)[5][5], double ***AScal, int flagImplicit);

	/*
	 * Method: calcViscousFluxNS1D_AD
	 * ------------------------------
	 * Calculate viscous flux in the so-called "1D style"
	 * Original code = calcViscousFluxNS_AD() in JoeWithModelsAD.cpp
	 */
	void calcViscousFluxNS1D_AD(const int icvCenter, REALA &rhs_rho, REALA rhs_rhou[3], REALA &rhs_rhoE,
			REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE,
			double (*A)[5][5], int flagImplicit);

	/*
	 * Method: freeMemForCalcEulerFlux1D_AD
	 * ------------------------------------
	 * Delete the memory allocated in the calcEulerFlux1D_AD() function
	 * (The arrays will NOT be re-initialized by NULL)
	 */
	void freeMemForCalcEulerFlux1D_AD(double (*Apl)[5], double (*Ami)[5], double (*AplScal)[6], double (*AmiScal)[6],
		REALQS *FrhoScal, REALQS *Scalar0, REALQS *Scalar1, REALQS *ScalCV0, REALQS *ScalCV1, double *ScalConvTerm, const int nScal);

	/*
	 * Method: calcViscousFluxScalar_new1D_AD
	 * --------------------------------------
	 * Original code: calcViscousFluxScalar_new_AD in ScalarsAD.cpp
	 */
	void calcViscousFluxScalar_new1D_AD(const int icvCenter, const int iScal, adouble &rhs_rhoScal, double *AScal, ScalarTranspEq &transpScal, ScalarTranspEq_AD &transpScal_AD, int flagImplicit);

	/*
	 * Method: showMessageParallel
	 * ---------------------------
	 * Each CPU core shows its message (usually, error or warning message) in a mpi_rank order
	 *
	 * Note: 1. All the CPU cores must call this method (Otherwise, communication error)
	 *       2. If the string is empty, NO message will be shown
	 */
	void showMessageParallel(string& message);
	void showMessageParallel(vector<string>& messageVector, const size_t maxMessageNum, char* messageVectorName = NULL);
};

#endif /* IKEWITHMODELS_H_ */
