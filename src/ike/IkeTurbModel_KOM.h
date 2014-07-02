/*
 * IkeTurbModel_KOM.h
 *
 *  Created on: Jan 15, 2013
 *      Author: ikj
 */

#ifndef IKETURBMODEL_KOM_H_
#define IKETURBMODEL_KOM_H_

#include "JOE/ADJOINT_FILES/UgpWithCvCompFlowAD.h"
#include "JOE/ADJOINT_FILES/TurbModel_KOM_AD.h"

#include "IkeWithPsALC.h"

/**
 * wilcox k-omega model as in Turbulence Modeling for CFD, Third Edition, Wilcox 2006
 */

class IkeRansTurbKOm_AD: public RansTurbKOm_AD //, public UgpWithCvCompFlow_AD //, public RansTurbKOm
{
public:
	/*
	 * constructors
	 */
	IkeRansTurbKOm_AD() {
		if (mpi_rank == 0)
			cout << "IkeRansTurbKOm_AD()" << endl;

#ifdef USE_MEM_SAVING_ADVAR_1D_
		kine = NULL;
		grad_kine = NULL;
		kine_diff = NULL;

		omega = NULL;
		grad_omega = NULL;
		omega_diff = NULL;
#endif
	}

	/*
	 * destructors
	 */
	virtual ~IkeRansTurbKOm_AD() {
	}

public:
	// ============================================
	// METHODS WHICH DON'T COME FROM RansTurbKOm_AD
	// ============================================
	/*
	 * Method: checkScalarsMemAD
	 * -------------------------
	 * This method will be called during the Newton iterations in the pseudo-arclength method just before calculating the RHS
	 * (For the details, see IkeWithPsALC_AD::calcJacobian1DAD(), etc.
	 */
	virtual void checkScalarsMemAD() {
		if (kine == NULL || omega == NULL) {
			cout<< "ERROR in IkeRansTurbKOm_AD::checkScalarsMemAD(): kine or omega is NULL"<< endl;
			assert(false);
		}
	}

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
        newRelax = oldRelax;

        int omega_index = getScalarTransportIndex("omega");
        double omegaRef = 1.0e-14;
        Param *pmy;
        if (getParam(pmy, "INITIAL_CONDITION_TURB"))
            omegaRef = pmy->getDouble(2); 

        double omegaClipBound = clipParameter*omegaRef; // If omega is less than this value, consider it as negative

        double newOmegacv = UgpWithCvCompFlow::scalarTranspEqVector[omega_index].phi[icv]; // omega after solving the linear system in the Newton method
        double oldOmegacv = qArray[icv+5+omega_index];  // omega before solving the linear system in the Newton method

        if(newOmegacv < omegaClipBound) {
            double tempRelax = safeParameter*(oldOmegacv-omegaClipBound)/fabs(delQ[icv+5+omega_index]);
            newRelax = min(oldRelax, tempRelax);

            return 1;
        }

        return 0;
    }

    /*
     * Method: checkNegativeTurbScalarFa_JOE
     * -------------------------------------
     * Return: number of negative values found at the face 
     */
    virtual int checkNegativeTurbScalarFa_JOE(const int ifa, const int icv0, const int icv1, const bool bothFaces) {
        assert( icv0 >= 0 );
        if(bothFaces)
            assert( icv1 >= 0 );

		int kine_index = getScalarTransportIndex("kine");
		int omega_index = getScalarTransportIndex("omega");

        int negativeFound = 0;

        // face unit normal and area...
        double nVec[3] = {0.0, 0.0, 0.0};
        double area    = normVec3d(nVec, fa_normal[ifa]);
        double r0[3] = {0.0, 0.0, 0.0}, r1[3] = {0.0, 0.0, 0.0};
        vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);
        vecMinVec3d(r1, x_fa[ifa], x_cv[icv1]);

        // check omega at each face
        int iScal = omega_index;
        // ----------------------------------------
        // left side: check omega
        // ----------------------------------------
        double (*grad_phi)[3] = UgpWithCvCompFlow::scalarTranspEqVector[iScal].grad_phi;

        double scalar0 = UgpWithCvCompFlow::scalarTranspEqVector[iScal].phi[icv0];
        double rho0    = UgpWithCvCompFlow::rho[icv0];

        if (strcmp(scalarTranspEqVector[iScal].reconstruction.c_str(), "CONSERVATIVE")==0)
            scalar0 = (rho0 * scalar0 + vecDotVec3d(r0, grad_phi[icv0])) / rho0;
        else if (strcmp(scalarTranspEqVector[iScal].reconstruction.c_str(), "STANDARD")==0)
            scalar0 += vecDotVec3d(r0, grad_phi[icv0]);

        if(scalar0 <= 0.0) 
            ++negativeFound;

        // ----------------------------------------
        // right side: check omega if bothFaces==true
        // ----------------------------------------
        if(bothFaces) {
            double scalar1 = UgpWithCvCompFlow::scalarTranspEqVector[iScal].phi[icv1];
            double rho1    = UgpWithCvCompFlow::rho[icv1];

            if (strcmp(scalarTranspEqVector[iScal].reconstruction.c_str(), "CONSERVATIVE")==0)
                scalar1 = (rho1 * scalar1 + vecDotVec3d(r1, grad_phi[icv1])) / rho1;
            else if (strcmp(scalarTranspEqVector[iScal].reconstruction.c_str(), "STANDARD")==0)
                scalar1 += vecDotVec3d(r1, grad_phi[icv1]);

            if(scalar1 <= 0.0)
                ++negativeFound;
        }

        // ----------------------------------------
        // At the face: check mut_fa
        // ----------------------------------------
        double mut_fa = 0.0;

        double kine0 = UgpWithCvCompFlow::scalarTranspEqVector[kine_index].phi[icv0]; 
        double kine1 = UgpWithCvCompFlow::scalarTranspEqVector[kine_index].phi[icv1]; 
        double omega0 = UgpWithCvCompFlow::scalarTranspEqVector[omega_index].phi[icv0]; 
        double omega1 = UgpWithCvCompFlow::scalarTranspEqVector[omega_index].phi[icv1]; 

        double rho1 = UgpWithCvCompFlow::rho[icv1];
        double strMag0 = UgpWithCvCompFlow::strMag[icv0];
        double strMag1 = UgpWithCvCompFlow::strMag[icv1];

 		// internal faces
		if(bothFaces) {
			double dx0[3], dx1[3];
			vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
			vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
			double w0 = sqrt(vecDotVec3d(dx0, dx0));
			double w1 = sqrt(vecDotVec3d(dx1, dx1));
			double ws = w0 + w1;
			w0 /= ws;
			w1 /= ws;

			double rho_fa  = w1 * rho0 + w0 * rho1;
			double kine_fa = w1 * kine0 + w0 * kine1;
			double om_fa   = w1 * omega0 + w0 * omega1;
			double strMag_fa = w1 * strMag0 + w0 * strMag1;

			if (KOM_RealizableConstraint == 1) {
				double omega_tilde = max(om_fa, cLim * strMag_fa / sqrt(betaStar));
				mut_fa = min(rho_fa * kine_fa / omega_tilde, 100.0);
			} else if (KOM_RealizableConstraint == 2) {
				double TS = min(1.0 / om_fa, 0.6 / (sqrt(6.0) * strMag_fa));
				mut_fa = rho_fa * kine_fa * TS;
			} else {
				mut_fa = min(rho_fa * kine_fa / om_fa, 100.0);
			}

            if(mut_fa <= 0.0)
                ++negativeFound;
		} else { // boundary faces: Check only for the ghost boundaries
                 //                 Note that mut_fa is zero at WALL, and it is usually not a problem for the other fake boundaries (e.g. SYMMETRY)
            if(icv1>=ncv_g && icv1<ncv_ggf) {
                if (KOM_RealizableConstraint == 1) {
                    double omega_tilde = max(omega1, cLim * strMag0 / sqrt(betaStar));
                    mut_fa = min(rho1 * kine1 / omega_tilde, 100.0); // zero order extrapolation for others
                } else if (KOM_RealizableConstraint == 2) {
                    double TS = min(1.0 / omega1, 0.6 / (sqrt(6.0) * strMag0));
                    mut_fa = rho1 * kine1 * TS;
                } else
                    mut_fa = min(rho1 * kine1 / omega1, 100.0); // zero order extrapolation for others
            }
		}

        return negativeFound;
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
        static int iterNewtonSaved = -1; // readBarrierParamTurbScalars() can be called many times during single iterNewton (for example, called for ncv-times if the "1D-style" is used).
                                         // Thus, all the messages will be also shown many times.
                                         // This variable will be used to print-out messages just once for a iterNewton.
        static int iterNewtonBarrierStart = -1;

        int debugLevel = getDebugLevel();

        coeffScalars.resize(2);

        if(iterNewton == 0) {
            iterNewtonSaved        = -1; // Reset to -1
            iterNewtonBarrierStart = -1; // Reset to -1
        }

        if (!checkParam("BARRIER_SOURCE_TURB_SCALARS")) {
            ParamMap::add("BARRIER_SOURCE_TURB_SCALARS  FUNCTION=NO_METHOD  MAX_ITER=0  THRESHOLD_RESID=1.0e-6  COEFF_KINE=1.0  COEFF_OMEGA=1.0"); // add default values
            if (mpi_rank == 0)
                cout<< "WARNING: added keyword \"BARRIER_SOURCE_TURB_SCALARS  METHOD=NO_METHOD  MAX_ITER=0  THRESHOLD_RESID=1.0e-6  COEFF_KINE=1.0  COEFF_OMEGA=1.0\""<< " to parameter map!" << endl;
        }

        barrierFunctionName = getParam("BARRIER_SOURCE_TURB_SCALARS")->getString("FUNCTION"); // If functionName=="NO_METHOD", skip the barrier

        if(strcmp(barrierFunctionName.c_str(),"NO_METHOD") != 0) {
            maxIter         = getParam("BARRIER_SOURCE_TURB_SCALARS")->getInt("MAX_ITER"); // If maxIter==0, also skip the barrier
            threshold_resid = getParam("BARRIER_SOURCE_TURB_SCALARS")->getDouble("THRESHOLD_RESID"); // If resid_norm_tot < threshold_resid, also skip the barrier
            coeffScalars[0] = getParam("BARRIER_SOURCE_TURB_SCALARS")->getDouble("COEFF_KINE");
            coeffScalars[1] = getParam("BARRIER_SOURCE_TURB_SCALARS")->getDouble("COEFF_OMEGA");

            if(maxIter>0 && iterNewton<=maxIter && resid_norm_tot>=threshold_resid) {
                if (checkParam("BARRIER_SOURCE_TURB_SCALARS_RAMP")) {  // Note: "BARRIER_SOURCE_NS_RAMP  AFTER_NEWTON_ITER=2  FACTOR_RESID=0.1  MIN_COEFF_RHO=1.0e-12  MIN_COEFF_PRESS=1.0e-12"
                    int    afterNewtonIter = getParam("BARRIER_SOURCE_TURB_SCALARS_RAMP")->getInt("AFTER_NEWTON_ITER");
                    double threshResid     = getParam("BARRIER_SOURCE_TURB_SCALARS_RAMP")->getDouble("THRESH_RESID");
                    double factorResid     = getParam("BARRIER_SOURCE_TURB_SCALARS_RAMP")->getDouble("FACTOR_RESID");
                    double minCoeffKine    = getParam("BARRIER_SOURCE_TURB_SCALARS_RAMP")->getDouble("MIN_COEFF_KINE");
                    double minCoeffOmega   = getParam("BARRIER_SOURCE_TURB_SCALARS_RAMP")->getDouble("MIN_COEFF_OMEGA");

                    if((factorResid<=0.0 || factorResid>1.0) && mpi_rank==0)
                        if(iterNewtonSaved != iterNewton)
                            cout<<"WARNING in IkeWithPsALC_AD::readBarrierParamTurbScalars(): BARRIER_SOURCE_NS_RAMP-->FACTOR_RESID is "<<factorResid<<endl;

                    if(iterNewton>afterNewtonIter || resid_norm_tot<threshResid) {
                        if(iterNewtonBarrierStart == -1)
                            iterNewtonBarrierStart = iterNewton-1;

                        coeffScalars[0] = max(coeffScalars[0]*pow(factorResid, (double) iterNewton-iterNewtonBarrierStart), minCoeffKine );
                        coeffScalars[1] = max(coeffScalars[1]*pow(factorResid, (double) iterNewton-iterNewtonBarrierStart), minCoeffOmega);

                        if(debugLevel>0 && mpi_rank==0)
                            if(iterNewtonSaved != iterNewton)
                                if(coeffScalars[0]>minCoeffKine || coeffScalars[1]>minCoeffOmega)
                                    printf("                         Reduce TURB SCALAR BARRIER COEFF: COEFF_KINE=%.5e, COEFF_OMEGA=%.5e\n", coeffScalars[0], coeffScalars[1]);
                    }
                }

                iterNewtonSaved = iterNewton;
                return true;
            } else {
                iterNewtonSaved = iterNewton;
                return false;
            }
        }

        iterNewtonSaved = iterNewton;
        return false;
    }

    /*
     * Method: barrierSourceTurbScalars
     * --------------------------------
     * Add barrier functions to the RHS of the scalars equations
     */
    void barrierSourceTurbScalars(double* rhs, const int nScal, const int iterNewton, const double residNormTotOld) {
        int debugLevel = getDebugLevel();

		int kine_index = getScalarTransportIndex("kine");       assert(kine_index>-1);
		int omega_index = getScalarTransportIndex("omega");     assert(omega_index>-1);

        double kineRef = 1.0e-14;
        Param *pmy;
        if (getParam(pmy, "INITIAL_CONDITION_TURB"))
            kineRef = pmy->getDouble(1);
        double omegaRef = 1.0e-14;
        if (getParam(pmy, "INITIAL_CONDITION_TURB"))
            omegaRef = pmy->getDouble(2);


        string functionName;    // If functionName=="NO_METHOD",            skip the barrier
        int    maxIter;         // If maxIter==0,                           skip the barrier
        double threshold_resid; // If resid_norm_tot_old < threshold_resid, skip the barrier
        vector<double> coeffScalars;
        bool useBarrier = readBarrierParamTurbScalars(functionName, maxIter, threshold_resid, coeffScalars, iterNewton, residNormTotOld);
        // Note: iterNewton and residNormTotOld are member variables of the IkeWithPsALC_AD class

        if(useBarrier) {
            int nVars = 5+nScal; 

            double *kineArray  = NULL;
            double *omegaArray = NULL;
            if(kine_index>=0)
                kineArray  = scalarTranspEqVector[kine_index].phi;
            if(omega_index>=0)
                omegaArray = scalarTranspEqVector[omega_index].phi;

//            myBarrierKineSourceSumJOE  = 0.0;
//            myBarrierOmegaSourceSumJOE = 0.0;

            if (strcmp(functionName.c_str(), "LOG") == 0) {
                for (int icv=0; icv<ncv; ++icv) {
                    // source term in the kine equation
                    double kineSource = max(-coeffScalars[0]*log(fabs(kineArray[icv]/kineRef)), 0.0)*cv_volume[icv];
                    rhs[icv*nVars+5+kine_index]  += kineSource;

                    // source term in the omega equation
                    double omegaSource = max(-coeffScalars[1]*log(fabs(omegaArray[icv]/omegaRef)), 0.0)*cv_volume[icv];
                    rhs[icv*nVars+5+omega_index] += omegaSource;

                    // statistics
                }
            } else if (strcmp(functionName.c_str(), "RECIPROCAL") == 0) {
                for (int icv=0; icv<ncv; ++icv) {
                    // source term in the kine equation
                    double kineSource = coeffScalars[0]/(fabs(kineArray[icv]/kineRef))*cv_volume[icv];
                    rhs[icv*nVars+5+kine_index]  += kineSource;
                    // source term in the omega equation
                    double omegaSource = coeffScalars[1]/(fabs(omegaArray[icv]/omegaRef))*cv_volume[icv];
                    rhs[icv*nVars+5+omega_index] += omegaSource;

                    // statistics
                }
            } else {
                if(debugLevel>1 && mpi_rank==0)
                    cout<<"barrierSourceTurbScalars() is not active: No barrier method"<<endl;
            }
        } else {
            if(debugLevel>1 && mpi_rank==0)
                cout<<"IkeWithPsALC_AD::barrierSourceTurbScalars(): Barrier function = NO_METHOD"<<endl;
        }
    }

    /*
     * Method: barrierSourceTurbScalars_AD
     * -----------------------------------
     * Add barrier functions to the RHS of the scalars equations
     */
    void barrierSourceTurbScalars_AD(REALAS **rhs_rhoScal_AD, const int nScal, const int iterNewton, const double residNormTotOld) {
        int debugLevel = getDebugLevel();

		int kine_index = getScalarTransportIndex("kine");       assert(kine_index>-1);
		int omega_index = getScalarTransportIndex("omega");     assert(omega_index>-1);

        double kineRef = 1.0e-14;
        Param *pmy;
        if (getParam(pmy, "INITIAL_CONDITION_TURB"))
            kineRef = pmy->getDouble(1);
        double omegaRef = 1.0e-14;
        if (getParam(pmy, "INITIAL_CONDITION_TURB"))
            omegaRef = pmy->getDouble(2);

        string functionName;    // If functionName=="NO_METHOD",            skip the barrier
        int    maxIter;         // If maxIter==0,                           skip the barrier
        double threshold_resid; // If resid_norm_tot_old < threshold_resid, skip the barrier
        vector<double> coeffScalars;
        bool useBarrier = readBarrierParamTurbScalars(functionName, maxIter, threshold_resid, coeffScalars, iterNewton, residNormTotOld);
        // Note: iterNewton and residNormTotOld are member variables of the IkeWithPsALC_AD class

        if(useBarrier) {
            int nVars = 5+nScal; 

//            myBarrierKineSourceSumJOE  = 0.0;
//            myBarrierOmegaSourceSumJOE = 0.0;

            if (strcmp(functionName.c_str(), "LOG") == 0) {
                for (int icv=0; icv<ncv; ++icv) {
                    adouble kinecv_AD  = 0.0;
                    adouble omegacv_AD = 0.0;
        #ifdef USE_MEM_SAVING_ADVAR_1D_
                    kinecv_AD  = (*kine)[icv]; // Note: here, "kine" is an adouble array
                    omegacv_AD = (*omega)[icv]; 
        #else
                    kinecv_AD  = kine[icv]; // Note: here, "kine" is an adouble array
                    omegacv_AD = omega[icv];
        #endif

                    // source term in the kine equation
                    adouble kineSource = max(-coeffScalars[0]*log(fabs(kinecv_AD/kineRef)), 0.0)*cv_volume[icv];
                    rhs_rhoScal_AD[kine_index][icv]  += kineSource;

                    // source term in the omega equation
                    adouble omegaSource = max(-coeffScalars[1]*log(fabs(omegacv_AD/omegaRef)), 0.0)*cv_volume[icv];
                    rhs_rhoScal_AD[omega_index][icv] += omegaSource;

                    // statistics
                }
            } else if (strcmp(functionName.c_str(), "RECIPROCAL") == 0) {
                for (int icv=0; icv<ncv; ++icv) {
                    adouble kinecv_AD  = 0.0;
                    adouble omegacv_AD = 0.0;
        #ifdef USE_MEM_SAVING_ADVAR_1D_
                    kinecv_AD  = (*kine)[icv]; // Note: here, "kine" is an adouble array
                    omegacv_AD = (*omega)[icv]; 
        #else
                    kinecv_AD  = kine[icv]; // Note: here, "kine" is an adouble array
                    omegacv_AD = omega[icv];
        #endif

                    // source term in the kine equation
                    adouble kineSource = coeffScalars[0]/(fabs(kinecv_AD/kineRef))*cv_volume[icv];
                    rhs_rhoScal_AD[kine_index][icv]  += kineSource;
                    // source term in the omega equation
                    adouble omegaSource = coeffScalars[1]/(fabs(omegacv_AD/omegaRef))*cv_volume[icv];
                    rhs_rhoScal_AD[omega_index][icv] += omegaSource;

                    // statistics
                }
            } else {
                if(debugLevel>1 && mpi_rank==0)
                    cout<<"barrierSourceTurbScalars_AD() is not active: No barrier method"<<endl;
            }
        } else {
            if(debugLevel>1 && mpi_rank==0)
                cout<<"IkeWithPsALC_AD::barrierSourceTurbScalars_AD(): Barrier function = NO_METHOD"<<endl;
        }
    }

    /*
     * Method: barrierSourceTurbScalars1D_AD
     * -------------------------------------
     * Add barrier functions to the RHS of the scalars equations
     */
    void barrierSourceTurbScalars1D_AD(const int icvCenter, REALAS* rhs_rhoScal_AD, const int nScal, const int iterNewton, const double residNormTotOld) {
        int debugLevel = getDebugLevel();

		int kine_index = getScalarTransportIndex("kine");       assert(kine_index>-1);
		int omega_index = getScalarTransportIndex("omega");     assert(omega_index>-1);

        double kineRef = 1.0e-14;
        Param *pmy;
        if (getParam(pmy, "INITIAL_CONDITION_TURB"))
            kineRef = pmy->getDouble(1);
        double omegaRef = 1.0e-14;
        if (getParam(pmy, "INITIAL_CONDITION_TURB"))
            omegaRef = pmy->getDouble(2);

        string functionName;    // If functionName=="NO_METHOD",            skip the barrier
        int    maxIter;         // If maxIter==0,                           skip the barrier
        double threshold_resid; // If resid_norm_tot_old < threshold_resid, skip the barrier
        vector<double> coeffScalars;
        bool useBarrier = readBarrierParamTurbScalars(functionName, maxIter, threshold_resid, coeffScalars, iterNewton, residNormTotOld);
        // Note: iterNewton and residNormTotOld are member variables of the IkeWithPsALC_AD class

        if(useBarrier) {
            int nVars = 5+nScal;

//            myBarrierKineSourceSumJOE  = 0.0;
//            myBarrierOmegaSourceSumJOE = 0.0;

            adouble kinecv_AD  = 0.0;
            adouble omegacv_AD = 0.0;
#ifdef USE_MEM_SAVING_ADVAR_1D_
            kinecv_AD  = (*kine)[icvCenter]; // Note: here, "kine" is an adouble array
            omegacv_AD = (*omega)[icvCenter]; 
#else
            kinecv_AD  = kine[icvCenter]; // Note: here, "kine" is an adouble array
            omegacv_AD = omega[icvCenter];
#endif

            if (strcmp(functionName.c_str(), "LOG") == 0) {
                // source term in the kine equation
                adouble kineSource = max(-coeffScalars[0]*log(fabs(kinecv_AD/kineRef)), 0.0)*cv_volume[icvCenter];
                rhs_rhoScal_AD[kine_index]  += kineSource;

                // source term in the omega equation
                adouble omegaSource = max(-coeffScalars[1]*log(fabs(omegacv_AD/omegaRef)), 0.0)*cv_volume[icvCenter];
                rhs_rhoScal_AD[omega_index] += omegaSource;

                // statistics
            } else if (strcmp(functionName.c_str(), "RECIPROCAL") == 0) {
                for (int icv=0; icv<ncv; ++icv) {
                    // source term in the kine equation
                    adouble kineSource = coeffScalars[0]/(fabs(kinecv_AD/kineRef))*cv_volume[icv];
                    rhs_rhoScal_AD[kine_index]  += kineSource;
                    // source term in the omega equation
                    adouble omegaSource = coeffScalars[1]/(fabs(omegacv_AD/omegaRef))*cv_volume[icv];
                    rhs_rhoScal_AD[omega_index] += omegaSource;

                    // statistics
                }
            } else {
                if(debugLevel>1 && mpi_rank==0)
                    cout<<"barrierSourceTurbScalars1D_AD() is not active: No barrier method"<<endl;
            }
        } else {
            if(debugLevel>1 && mpi_rank==0)
                cout<<"IkeWithPsALC_AD::barrierSourceTurbScalars1D_AD(): Barrier function = NO_METHOD"<<endl;
        }
    }

	// ============================================
	// METHODS WHICH COME FROM RansTurbKOm_AD
	// ============================================
#ifdef USE_MEM_SAVING_ADVAR_1D_
	/*
	 * Method: initialHookScalarRansTurbModel1D_AD()
	 * ---------------------------------------------
	 * This method is doing the exact same thing as initialHookScalarRansCombModel_AD()
	 * but the implementation that works only for the ICVs stored in nbocv2ff
	 * Original code = 1. initialHookScalarRansCombModel_AD in UgpWithCvCompFlowAD.h
	 *                 2. initialHookScalarRansCombModel1D_AD in IkeTurbModel_KOM.h
	 */
	virtual void initialHookScalarRansTurbModel1D_AD(vector<int> &nbocv2ff, bool &firstCall) {
		// debug level
		int debugLevel = getDebugLevel();

		// connect pointers
		int nScal = scalarTranspEqVector.size();

		for (int i = 0; i < nScal; i++) {
			if (strcmp(scalarTranspEqVector_AD[i].name, "kine") == 0) {
				assert(kine == NULL);
				assert(scalarTranspEqVector_AD[i].phi.size() == nbocv2ff.size());

				kine = &(scalarTranspEqVector_AD[i].phi);
				grad_kine = &(scalarTranspEqVector_AD[i].grad_phi);
				kine_diff = &(scalarTranspEqVector_AD[i].diff);
				if (debugLevel > 1 && mpi_rank == 0)
					cout<< "IkeRansTurbKOm_AD::initialHookScalarRansTurbModel1D_AD(): Connected Scalar Pointer "
					<< scalarTranspEqVector_AD[i].name << endl;
			}
			if (strcmp(scalarTranspEqVector_AD[i].name, "omega") == 0) {
				assert(omega == NULL);
				assert(scalarTranspEqVector_AD[i].phi.size() == nbocv2ff.size());

				omega = &(scalarTranspEqVector_AD[i].phi);
				grad_omega = &(scalarTranspEqVector_AD[i].grad_phi);
				omega_diff = &(scalarTranspEqVector_AD[i].diff);
				if (debugLevel > 1 && mpi_rank == 0)
					cout<< "IkeRansTurbKOm_AD::initialHookScalarRansTurbModel1D_AD(): Connected Scalar Pointer "
					<< scalarTranspEqVector_AD[i].name << endl;
			}
		}

		// update wall distance which is required for the k-omega model
		if(firstCall)
			updateCvDataG1G2(wallDist, REPLACE_DATA); // VERY IMPORTANT TO DO THIS AS THIS IS REQUIRED TO SET BOUNDARY INFO

		// update "initialCall"
		firstCall = false;
	}

	/*
	 * Method: finalHookScalarRansTurbModel1D_AD()
	 * ---------------------------------------------
	 *
	 */
	virtual void finalHookScalarRansTurbModel1D_AD() {
		// debug level
		int debugLevel = getDebugLevel();

		if(!(*kine).empty())
			(*kine).clear();
		kine = NULL;
		if(!(*grad_kine).empty())
			(*grad_kine).clear();
		grad_kine = NULL;
		if(!(*kine_diff).empty())
			(*kine_diff).clear();
		kine_diff = NULL;
		if (debugLevel > 1 && mpi_rank == 0)
			cout<<"IkeRansTurbKOm_AD::finalHookScalarRansTurbModel1D_AD(): Reset kine-related pointers as NULL "<<endl;

		if(!(*omega).empty())
			(*omega).clear();
		omega = NULL;
		if(!(*grad_omega).empty())
			(*grad_omega).clear();
		grad_omega = NULL;
		if(!(*omega_diff).empty())
			(*omega_diff).clear();
		omega_diff = NULL;
		if (debugLevel > 1 && mpi_rank == 0)
			cout<<"IkeRansTurbKOm_AD::finalHookScalarRansTurbModel1D_AD(): Reset omega-related pointers as NULL "<<endl;
	}
#endif

	/*
	 * Method: initialHookScalarRansTurbModel1D_AD()
	 * ---------------------------------------------
	 * This method is doing the exact same thing as initialHookScalarRansCombModel_AD()
	 * but has extra consideration that makes this method run only at the first call
	 * Original code = 1. initialHookScalarRansCombModel_AD in UgpWithCvCompFlowAD.h
	 *                 2. initialHookScalarRansCombModel1D_AD in IkeTurbModel_KOM.h
	 */
	virtual void initialHookScalarRansTurbModel1D_AD(bool &firstCall) {
		if (firstCall) {
			// debug level
			int debugLevel = getDebugLevel();

			// connect pointers
			int nScal = scalarTranspEqVector.size();

			for (int i = 0; i < nScal; i++) {
				if (strcmp(scalarTranspEqVector_AD[i].name, "kine") == 0) {
#ifdef USE_MEM_SAVING_ADVAR
					assert(kine == NULL);
					kine = &(scalarTranspEqVector_AD[i].phi);
					grad_kine = &(scalarTranspEqVector_AD[i].grad_phi);
					kine_diff = &(scalarTranspEqVector_AD[i].diff);
#else
					assert(kine == NULL);
					kine = scalarTranspEqVector_AD[i].phi;
					grad_kine = scalarTranspEqVector_AD[i].grad_phi;
					kine_diff = scalarTranspEqVector_AD[i].diff;
#endif
					if (debugLevel > 1 && mpi_rank == 0)
						cout<< "IkeRansTurbKOm_AD::initialHookScalarRansTurbModel1D_AD(): Connected Scalar Pointer "
							<< scalarTranspEqVector_AD[i].name << endl;
				}
				if (strcmp(scalarTranspEqVector_AD[i].name, "omega") == 0) {
#ifdef USE_MEM_SAVING_ADVAR
					omega = &(scalarTranspEqVector_AD[i].phi);
					grad_omega = &(scalarTranspEqVector_AD[i].grad_phi);
					omega_diff = &(scalarTranspEqVector_AD[i].diff);
#else
//					assert(omega == NULL);
					omega = scalarTranspEqVector_AD[i].phi;
					grad_omega = scalarTranspEqVector_AD[i].grad_phi;
					omega_diff = scalarTranspEqVector_AD[i].diff;
#endif
					if (debugLevel > 1 && mpi_rank == 0)
						cout<< "IkeRansTurbKOm_AD::initialHookScalarRansTurbModel1D_AD(): Connected Scalar Pointer "
							<< scalarTranspEqVector_AD[i].name << endl;
				}
			}

			// update wall distance which is required for the k-omega model
			updateCvDataG1G2(wallDist, REPLACE_DATA); // VERY IMPORTANT TO DO THIS AS THIS IS REQUIRED TO SET BOUNDARY INFO

			// update "initialCall"
			firstCall = false;
		}
	}

	/*
	 * Method: initialHookScalarRansTurbModel_AD()
	 * -------------------------------------------
	 *
	 */
	virtual void initialHookScalarRansTurbModel_AD() {
		// bebug level
		int debugLevel = getDebugLevel();

		// connect pointers
		int nScal = scalarTranspEqVector.size();

		for (int i = 0; i < nScal; i++) {
			if (strcmp(scalarTranspEqVector_AD[i].name, "kine") == 0) {
#ifdef USE_MEM_SAVING_ADVAR
				kine = &(scalarTranspEqVector_AD[i].phi);
				grad_kine = &(scalarTranspEqVector_AD[i].grad_phi);
				kine_diff = &(scalarTranspEqVector_AD[i].diff);
#else
				kine = scalarTranspEqVector_AD[i].phi;
				grad_kine = scalarTranspEqVector_AD[i].grad_phi;
				kine_diff = scalarTranspEqVector_AD[i].diff;
#endif
				if (debugLevel > 1 && mpi_rank == 0)
					cout
							<< "IkeRansTurbKOm_AD::initialHookScalarRansTurbModel_AD(): Connected Scalar Pointer "
							<< scalarTranspEqVector_AD[i].name << endl;
			}
			if (strcmp(scalarTranspEqVector_AD[i].name, "omega") == 0) {
#ifdef USE_MEM_SAVING_ADVAR
				omega = &(scalarTranspEqVector_AD[i].phi);
				grad_omega = &(scalarTranspEqVector_AD[i].grad_phi);
				omega_diff = &(scalarTranspEqVector_AD[i].diff);
#else
				omega = scalarTranspEqVector_AD[i].phi;
				grad_omega = scalarTranspEqVector_AD[i].grad_phi;
				omega_diff = scalarTranspEqVector_AD[i].diff;
#endif
				if (debugLevel > 1 && mpi_rank == 0)
					cout
							<< "IkeRansTurbKOm_AD::initialHookScalarRansTurbModel_AD(): Connected Scalar Pointer "
							<< scalarTranspEqVector_AD[i].name << endl;
			}
		}

        // Caculate wall distance
		updateCvDataG1G2(wallDist, REPLACE_DATA); // VERY IMPORTANT TO DO THIS AS THIS IS REQUIRED TO SET BOUNDARY INFO
	}

#ifdef USE_MEM_SAVING_ADVAR
	virtual void calcRansTurbViscMuet1D_AD(const int icvCenter, FaZone** zoneArrayBoundary, vector<int>& faInternal,
			vector<int>& faBoundary, ADscalar<REALQ> &rho, ADvector<REALQ> &rhou) {
		if (kine == NULL || omega == NULL) {
			cout<< "ERROR in IkeRansTurbKOm_AD::calcRansTurbViscMuet1D_AD(): kine or omega is NULL"<< endl;
			throw(IKEWITHMODELS_ERROR_CODE);
		}

		// update strain rate tensor
		calcCv2Grad1D_AD(icvCenter, grad_u, vel, limiterNavierS, sos, epsilonSDWLS);

		calcStrainRateAndDivergence1D_AD(icvCenter);

		calcVorticity1D_AD(icvCenter);

		// internal faces
		for (size_t index = 0; index < faInternal.size(); ++index) {
			int ifa = faInternal[index];
			int icv0 = cvofa[ifa][0];
			int icv1 = cvofa[ifa][1];

			double dx0[3], dx1[3];
			vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
			vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
			double w0 = sqrt(vecDotVec3d(dx0, dx0));
			double w1 = sqrt(vecDotVec3d(dx1, dx1));
			double ws = w0 + w1;
			w0 /= ws;
			w1 /= ws;

			adouble rho_fa = w1 * rho[icv0] + w0 * rho[icv1];
			adouble kine_fa = w1 * (*kine)[icv0] + w0 * (*kine)[icv1];
			adouble om_fa = w1 * (*omega)[icv0] + w0 * (*omega)[icv1];
			adouble strMag_fa = w1 * strMag[icv0] + w0 * strMag[icv1];

#ifdef USE_CONDASSIGN
			if (KOM_RealizableConstraint == 1) {
				adouble omega_tilde = fmax(om_fa, cLim * strMag_fa / sqrt(betaStar));
				mut_fa[ifa] = fmin(rho_fa * kine_fa / omega_tilde, 100.0);
			} else if (KOM_RealizableConstraint == 2) {
				adouble TS = fmin(1.0 / om_fa, 0.6 / (sqrt(6.0) * strMag_fa));
				mut_fa[ifa] = rho_fa * kine_fa * TS;
			} else {
				mut_fa[ifa] = fmin(rho_fa * kine_fa / om_fa, 100.0);
			}
#else
			if (KOM_RealizableConstraint == 1) {
				adouble omega_tilde = max(om_fa, cLim * strMag_fa / sqrt(betaStar));
				mut_fa[ifa] = min(rho_fa * kine_fa / omega_tilde, 100.0);
			} else if (KOM_RealizableConstraint == 2) {
				adouble TS = min(1.0 / om_fa, 0.6 / (sqrt(6.0) * strMag_fa));
				mut_fa[ifa] = rho_fa * kine_fa * TS;
			} else {
				mut_fa[ifa] = min(rho_fa * kine_fa / om_fa, 100.0);
			}
#endif
		}

		// boundary faces
		for (size_t index = 0; index < faBoundary.size(); ++index) {
			if (zoneArrayBoundary == NULL) {
				cout << "Error in calcRansTurbViscMuet1D_AD => zoneArrayBoundary is NULL" << endl;
				throw(-1);
			}
			int ifa = faBoundary[index];
			FaZone* zone = zoneArrayBoundary[index];

			if (getBoundaryType(ifa, zone) == WALL) {
				mut_fa[ifa] = 0.0; // set mut zero at walls
			} else {
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];

#ifdef USE_CONDASSIGN
				if (KOM_RealizableConstraint == 1) {
					adouble omega_tilde = fmax((*omega)[icv1], cLim * strMag[icv0] / sqrt(betaStar));
					mut_fa[ifa] = fmin(rho[icv1] * (*kine)[icv1] / omega_tilde, 100.0); // zero order extrapolation for others
				} else if (KOM_RealizableConstraint == 2) {
					adouble TS = fmin(1.0 / (*omega)[icv1], 0.6 / (sqrt(6.0) * strMag[icv0]));
					mut_fa[ifa] = rho[icv1] * (*kine)[icv1] * TS;
				} else {
					mut_fa[ifa] = fmin(rho[icv1] * (*kine)[icv1] / (*omega)[icv1], 100.0); // zero order extrapolation for others
				}
#else
				if (KOM_RealizableConstraint == 1) {
					adouble omega_tilde = max((*omega)[icv1], cLim * strMag[icv0] / sqrt(betaStar));
					mut_fa[ifa] = min(rho[icv1] * (*kine)[icv1] / omega_tilde, 100.0); // zero order extrapolation for others
				} else if (KOM_RealizableConstraint == 2) {
					adouble TS = min(1.0 / (*omega)[icv1], 0.6 / (sqrt(6.0) * strMag[icv0]));
					mut_fa[ifa] = rho[icv1] * (*kine)[icv1] * TS;
				} else {
					mut_fa[ifa] = min(rho[icv1] * (*kine)[icv1] / (*omega)[icv1], 100.0); // zero order extrapolation for others
				}
#endif
			}
		}
	}
#endif

#ifdef USE_MEM_SAVING_ADVAR
	virtual void calcRansTurbViscMuet1D_AD(const int icvCenter, FaZone** zoneArrayBoundary, vector<int>& faInternal,
			vector<int>& faBoundary, adouble *rho, adouble(*rhou)[3]) {
		if (kine == NULL || omega == NULL) {
			cout<< "ERROR in IkeRansTurbKOm_AD::calcRansTurbViscMuet1D_AD(): kine or omega is NULL"<< endl;
			throw(IKEWITHMODELS_ERROR_CODE);
		}
		// update strain rate tensor
		calcCv2Grad1D_AD(icvCenter, grad_u, vel, limiterNavierS, sos, epsilonSDWLS);
		calcStrainRateAndDivergence1D_AD(icvCenter);

		calcVorticity1D_AD(icvCenter);

		// internal faces
		for (size_t index = 0; index < faInternal.size(); ++index) {
			int ifa = faInternal[index];
			int icv0 = cvofa[ifa][0];
			int icv1 = cvofa[ifa][1];

			double dx0[3], dx1[3];
			vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
			vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
			double w0 = sqrt(vecDotVec3d(dx0, dx0));
			double w1 = sqrt(vecDotVec3d(dx1, dx1));
			double ws = w0 + w1;
			w0 /= ws;
			w1 /= ws;

			adouble rho_fa = w1 * rho[icv0] + w0 * rho[icv1];
			adouble kine_fa = w1 * (*kine)[icv0] + w0 * (*kine)[icv1];
			adouble om_fa = w1 * (*omega)[icv0] + w0 * (*omega)[icv1];
			adouble strMag_fa = w1 * strMag[icv0] + w0 * strMag[icv1];

			if (KOM_RealizableConstraint == 1) {
				adouble omega_tilde = max(om_fa, cLim * strMag_fa / sqrt(betaStar));
				mut_fa[ifa] = min(rho_fa * kine_fa / omega_tilde, 100.0);
			} else if (KOM_RealizableConstraint == 2) {
				adouble TS = min(1.0 / om_fa, 0.6 / (sqrt(6.0) * strMag_fa));
				mut_fa[ifa] = rho_fa * kine_fa * TS;
			} else {
				mut_fa[ifa] = min(rho_fa * kine_fa / om_fa, 100.0);
			}
		}

		// boundary faces
		for (size_t index = 0; index < faBoundary.size(); ++index) {
			if (zoneArrayBoundary == NULL) {
				cout << "Error in calcRansTurbViscMuet1D_AD => zoneArrayBoundary is NULL" << endl;
				throw(-1);
			}
			int ifa = faBoundary[index];
			FaZone* zone = zoneArrayBoundary[index];

			if (getBoundaryType(ifa, zone) == WALL)
				mut_fa[ifa] = 0.0; // set mut zero at walls
			else {
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];

				if (KOM_RealizableConstraint == 1) {
					adouble omega_tilde = max((*omega)[icv1], cLim * strMag[icv0] / sqrt(betaStar));
					mut_fa[ifa] = min(rho[icv1] * (*kine)[icv1] / omega_tilde, 100.0); // zero order extrapolation for others
				} else if (KOM_RealizableConstraint == 2) {
					adouble TS = min(1.0 / (*omega)[icv1], 0.6 / (sqrt(6.0) * strMag[icv0]));
					mut_fa[ifa] = rho[icv1] * (*kine)[icv1] * TS;
				} else
					mut_fa[ifa] = min(rho[icv1] * (*kine)[icv1] / (*omega)[icv1], 100.0); // zero order extrapolation for others
			}
		}
	}
#else
	virtual void calcRansTurbViscMuet1D_AD(const int icvCenter, FaZone** zoneArrayBoundary, vector<int>& faInternal,
			vector<int>& faBoundary, adouble *rho, adouble(*rhou)[3]) {
		if (kine == NULL || omega == NULL) {
			cout<< "ERROR in IkeRansTurbKOm_AD::calcRansTurbViscMuet1D_AD(): kine or omega is NULL"<< endl;
			throw(IKEWITHMODELS_ERROR_CODE);
		}
		// update strain rate tensor
		calcCv2Grad1D_AD(icvCenter, grad_u, vel, limiterNavierS, sos, epsilonSDWLS);
		calcStrainRateAndDivergence1D_AD(icvCenter);

		calcVorticity1D_AD(icvCenter);

		// internal faces
		for (size_t index = 0; index < faInternal.size(); ++index) {
			int ifa = faInternal[index];
			int icv0 = cvofa[ifa][0];
			int icv1 = cvofa[ifa][1];

			double dx0[3], dx1[3];
			vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
			vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
			double w0 = sqrt(vecDotVec3d(dx0, dx0));
			double w1 = sqrt(vecDotVec3d(dx1, dx1));
			double ws = w0 + w1;
			w0 /= ws;
			w1 /= ws;

			adouble rho_fa = w1 * rho[icv0] + w0 * rho[icv1];
			adouble kine_fa = w1 * kine[icv0] + w0 * kine[icv1];
			adouble om_fa = w1 * omega[icv0] + w0 * omega[icv1];
			adouble strMag_fa = w1 * strMag[icv0] + w0 * strMag[icv1];

			if (KOM_RealizableConstraint == 1) {
				adouble omega_tilde = max(om_fa, cLim * strMag_fa / sqrt(betaStar));
				mut_fa[ifa] = min(rho_fa * kine_fa / omega_tilde, 100.0);
			} else if (KOM_RealizableConstraint == 2) {
				adouble TS = min(1.0 / om_fa, 0.6 / (sqrt(6.0) * strMag_fa));
				mut_fa[ifa] = rho_fa * kine_fa * TS;
			} else {
				mut_fa[ifa] = min(rho_fa * kine_fa / om_fa, 100.0);
			}
		}

		// boundary faces
		for (size_t index = 0; index < faBoundary.size(); ++index) {
			if (zoneArrayBoundary == NULL) {
				cout << "Error in calcRansTurbViscMuet1D_AD => zoneArrayBoundary is NULL" << endl;
				throw(-1);
			}
			int ifa = faBoundary[index];
			FaZone* zone = zoneArrayBoundary[index];

			if (getBoundaryType(ifa, zone) == WALL)
				mut_fa[ifa] = 0.0; // set mut zero at walls
			else {
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];

				if (KOM_RealizableConstraint == 1) {
					adouble omega_tilde = max(omega[icv1], cLim * strMag[icv0] / sqrt(betaStar));
					mut_fa[ifa] = min(rho[icv1] * kine[icv1] / omega_tilde, 100.0); // zero order extrapolation for others
				} else if (KOM_RealizableConstraint == 2) {
					adouble TS = min(1.0 / omega[icv1], 0.6 / (sqrt(6.0) * strMag[icv0]));
					mut_fa[ifa] = rho[icv1] * kine[icv1] * TS;
				} else
					mut_fa[ifa] = min(rho[icv1] * kine[icv1] / omega[icv1], 100.0); // zero order extrapolation for others
			}
		}
	}
#endif

	/*
	 * Method: diffusivityHookScalarRansTurb1D_AD
	 * ------------------------------------------
	 *
	 */
#ifdef USE_MEM_SAVING_ADVAR
	virtual void diffusivityHookScalarRansTurb1D_AD(const int ifa, FaZone* zone, const string &name) {
		if ((name == "kine") || (name == "omega")) {
			double sigma;
			if (name == "kine")
				sigma = sigmaStar;
			if (name == "omega")
				sigma = sigmaOmega;

			// internal faces
			if ((ifa >= nfa_b && ifa < nfa) || (ifa >= nfa_b2 && ifa < nfa_b2gg)) {
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];

				double dx0[3], dx1[3];
				vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
				vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
				double w0 = sqrt(vecDotVec3d(dx0, dx0));
				double w1 = sqrt(vecDotVec3d(dx1, dx1));

				adouble rho_fa = (w1 * rho_AD[icv0] + w0 * rho_AD[icv1]) / (w0 + w1);
				adouble kine_fa = (w1 * (*kine)[icv0] + w0 * (*kine)[icv1]) / (w0 + w1);
				adouble om_fa = (w1 * (*omega)[icv0] + w0 * (*omega)[icv1]) / (w0 + w1);

				if (name == "kine")
					(*kine_diff)[ifa] = mul_fa[ifa] + sigma * rho_fa * kine_fa / om_fa;
				if (name == "omega")
					(*omega_diff)[ifa] = mul_fa[ifa] + sigma * rho_fa * kine_fa / om_fa;
			}
			// boundary faces
			else {
				if (zone->getKind() == FA_ZONE_BOUNDARY) {
					if (getBoundaryType(ifa, zone) == WALL) {
						int icv0 = cvofa[ifa][0];
						if (name == "kine")
							(*kine_diff)[ifa] = mul_fa[ifa];
						if (name == "omega")
							(*omega_diff)[ifa] = mul_fa[ifa];
					} else {
						int icv1 = cvofa[ifa][1];
						if (name == "kine")
							(*kine_diff)[ifa] = mul_fa[ifa] + sigma * rho_AD[icv1] * (*kine)[icv1] / (*omega)[icv1];
						if (name == "omega")
							(*omega_diff)[ifa] = mul_fa[ifa] + sigma * rho_AD[icv1] * (*kine)[icv1] / (*omega)[icv1];
					}
				}
			}
		}
	}
#else
	virtual void diffusivityHookScalarRansTurb1D_AD(const int ifa, FaZone* zone, const string &name) {
		if ((name == "kine") || (name == "omega")) {
			double sigma;
			if (name == "kine")
				sigma = sigmaStar;
			if (name == "omega")
				sigma = sigmaOmega;

			// internal faces
			if ((ifa >= nfa_b && ifa < nfa) || (ifa >= nfa_b2 && ifa < nfa_b2gg)) {
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];

				double dx0[3], dx1[3];
				vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
				vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
				double w0 = sqrt(vecDotVec3d(dx0, dx0));
				double w1 = sqrt(vecDotVec3d(dx1, dx1));

				adouble rho_fa = (w1 * rho_AD[icv0] + w0 * rho_AD[icv1]) / (w0 + w1);
				adouble kine_fa = (w1 * kine[icv0] + w0 * kine[icv1]) / (w0 + w1);
				adouble om_fa = (w1 * omega[icv0] + w0 * omega[icv1]) / (w0 + w1);

				if (name == "kine")
					kine_diff[ifa] = mul_fa[ifa] + sigma * rho_fa * kine_fa / om_fa;
				if (name == "omega")
					omega_diff[ifa] = mul_fa[ifa] + sigma * rho_fa * kine_fa / om_fa;
			}
			// boundary faces
			else {
				if (zone->getKind() == FA_ZONE_BOUNDARY) {
					if (getBoundaryType(ifa, zone) == WALL) {
						int icv0 = cvofa[ifa][0];
						if (name == "kine")
							kine_diff[ifa] = mul_fa[ifa];
						if (name == "omega")
							omega_diff[ifa] = mul_fa[ifa];
					} else {
						int icv1 = cvofa[ifa][1];
						if (name == "kine")
							kine_diff[ifa] = mul_fa[ifa] + sigma * rho_AD[icv1] * kine[icv1] / omega[icv1];
						if (name == "omega")
							omega_diff[ifa] = mul_fa[ifa] + sigma * rho_AD[icv1] * kine[icv1] / omega[icv1];
					}
				}
			}
		}
	}
#endif

	/*
	 * Method: calcTurbProd_AD
	 * -----------------------
	 *
	 */
#ifdef USE_MEM_SAVING_ADVAR
	virtual adouble calcTurbProd_AD(int icv) {
		if (KOM_RealizableConstraint == 1) {
			adouble omega_tilde = fmax((*omega)[icv],cLim * strMag[icv] / sqrt(betaStar));
			adouble mu_t = rho_AD[icv]*(*kine)[icv]/omega_tilde;
			return fmax(mu_t*strMag[icv]*strMag[icv] - 2.0/3.0*rho_AD[icv]*(*kine)[icv]*diverg[icv], 0.0);
		} else if (KOM_RealizableConstraint == 2) {
			adouble TS = fmin(1.0/((*omega)[icv]), 0.6/(sqrt(6.0)*strMag[icv]));
			adouble mu_t = rho_AD[icv]*(*kine)[icv]*TS;
			return fmax(mu_t*strMag[icv]*strMag[icv] - 2.0/3.0*rho_AD[icv]*(*kine)[icv]*diverg[icv], 0.0);
		} else {
			adouble mu_t = rho_AD[icv]*(*kine)[icv]/(*omega)[icv];
			return fmax(mu_t*strMag[icv]*strMag[icv] - 2.0/3.0*rho_AD[icv]*(*kine)[icv]*diverg[icv], 0.0);
		}
	}
#else
	virtual adouble calcTurbProd_AD(int icv) {
		if (KOM_RealizableConstraint == 1) {
			adouble omega_tilde = max(omega[icv],cLim * strMag[icv] / sqrt(betaStar));
			adouble mu_t = rho_AD[icv]*kine[icv]/omega_tilde;
			return max(mu_t*strMag[icv]*strMag[icv] - 2.0/3.0*rho_AD[icv]*kine[icv]*diverg[icv], 0.0);
		} else if (KOM_RealizableConstraint == 2) {
			adouble TS = min(1.0/(omega[icv]), 0.6/(sqrt(6.0)*strMag[icv]));
			adouble mu_t = rho_AD[icv]*kine[icv]*TS;
			return max(mu_t*strMag[icv]*strMag[icv] - 2.0/3.0*rho_AD[icv]*kine[icv]*diverg[icv], 0.0);
		} else {
			adouble mu_t = rho_AD[icv]*kine[icv]/omega[icv];
			return max(mu_t*strMag[icv]*strMag[icv] - 2.0/3.0*rho_AD[icv]*kine[icv]*diverg[icv], 0.0);
		}
	}
#endif

	/*
	 * Method: sourceHookScalarRansTurb_new_AD
	 * ---------------------------------------
	 *
	 */
#ifdef USE_MEM_SAVING_ADVAR
	virtual void sourceHookScalarRansTurb_new1D_AD(const int icv, adouble& rhs, double *A, const string &name, int flagImplicit) {
		if (name == "kine") {
			adouble src = calcTurbProd_AD(icv) - betaStar * rho_AD[icv] * (*omega)[icv] * (*kine)[icv];
			rhs += src * cv_volume[icv];

			if (flagImplicit && icv < ncv) {
				int noc00 = nbocv_i[icv];
				adouble dsrcdphi = -betaStar * (*omega)[icv];
				A[noc00] -= dsrcdphi.value() * cv_volume[icv];
			}
		}

		if (name == "omega") {
			adouble OM[3][3], STR_hat[3][3];

			adouble TS = 1.0;
			if (KOM_RealizableConstraint == 2)
				TS = fmin(1.0 / ((*omega)[icv]), 0.6 / (sqrt(6.0) * strMag[icv]));

			adouble sigmad;
			adouble crossDiff = vecDotVec3d_AD((*grad_kine)[icv], (*grad_omega)[icv]);

#ifdef USE_CONDASSIGN
			condassign(sigmad, crossDiff, sigmad0, 0.0);
#else
			if (crossDiff <= 0.0)
				sigmad = 0.0;
			else
				sigmad = sigmad0;
#endif

			adouble divU = grad_u[icv][0][0] + grad_u[icv][1][1] + grad_u[icv][2][2];

			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++) {
					OM[i][j] = 0.5 * (grad_u[icv][i][j] - grad_u[icv][j][i]);

					// strange: STR_hat subtracts 0.5 divergence instead 1/3, no bug!!!
					if (i == j)
						STR_hat[i][j] = 0.5 * (grad_u[icv][i][j] + grad_u[icv][j][i] - divU);
					else
						STR_hat[i][j] = 0.5 * (grad_u[icv][i][j] + grad_u[icv][j][i]);
				}

			adouble chiOm = 0.0;
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					for (int k = 0; k < 3; k++)
						chiOm += OM[i][j] * OM[j][k] * STR_hat[k][i];
			chiOm = fabs(chiOm / (betaStar * (*omega)[icv] * betaStar * (*omega)[icv] * betaStar * (*omega)[icv]));

			adouble fbeta = (1.0 + 85.0 * chiOm) / (1.0 + 100.0 * chiOm);
			adouble beta = beta0 * fbeta;

			adouble src = TS*alfa*(*omega)[icv]/(*kine)[icv]*calcTurbProd_AD(icv) - TS*beta*rho_AD[icv]*(*omega)[icv]*(*omega)[icv]
			             + sigmad*rho_AD[icv]/(*omega)[icv]*crossDiff;
			rhs += src * cv_volume[icv];

			if (flagImplicit && icv < ncv) {
				int noc00 = nbocv_i[icv];
				adouble dsrcdphi = -2.0 * TS * beta * (*omega)[icv];
				A[noc00] -= dsrcdphi.value() * cv_volume[icv];
			}
		}
	}
#else
	virtual void sourceHookScalarRansTurb_new1D_AD(const int icv, adouble& rhs, double *A, const string &name, int flagImplicit) {
		if (name == "kine") {
			adouble src = calcTurbProd_AD(icv)
					- betaStar * rho_AD[icv] * omega[icv] * kine[icv];
			rhs += src * cv_volume[icv];

			if (flagImplicit && icv < ncv) {
				int noc00 = nbocv_i[icv];
				adouble dsrcdphi = -betaStar * omega[icv];
				A[noc00] -= dsrcdphi.value() * cv_volume[icv];
			}
		}

		if (name == "omega") {
			adouble OM[3][3], STR_hat[3][3];

			adouble TS = 1.0;
			if (KOM_RealizableConstraint == 2)
				TS = min(1.0 / (omega[icv]), 0.6 / (sqrt(6.0) * strMag[icv]));

			adouble sigmad;
			adouble crossDiff = vecDotVec3d_AD(grad_kine[icv], grad_omega[icv]);

			if (crossDiff <= 0.0)
				sigmad = 0.0;
			else
				sigmad = sigmad0;

			adouble divU = grad_u[icv][0][0] + grad_u[icv][1][1]
					+ grad_u[icv][2][2];

			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++) {
					OM[i][j] = 0.5 * (grad_u[icv][i][j] - grad_u[icv][j][i]);

					// strange: STR_hat subtracts 0.5 divergence instead 1/3, no bug!!!
					if (i == j)
						STR_hat[i][j] =
								0.5
										* (grad_u[icv][i][j] + grad_u[icv][j][i]
												- divU);
					else
						STR_hat[i][j] = 0.5
								* (grad_u[icv][i][j] + grad_u[icv][j][i]);
				}

			adouble chiOm = 0.0;
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					for (int k = 0; k < 3; k++)
						chiOm += OM[i][j] * OM[j][k] * STR_hat[k][i];
			chiOm = fabs(
					chiOm
							/ (betaStar * omega[icv] * betaStar * omega[icv]
									* betaStar * omega[icv]));

			adouble fbeta = (1.0 + 85.0 * chiOm) / (1.0 + 100.0 * chiOm);
			adouble beta = beta0 * fbeta;

			adouble src = TS * alfa * omega[icv] / kine[icv]
					* calcTurbProd_AD(icv)
					- TS * beta * rho_AD[icv] * omega[icv] * omega[icv]
					+ sigmad * rho_AD[icv] / omega[icv] * crossDiff;

			rhs += src * cv_volume[icv];

			if (flagImplicit && icv < ncv) {
				int noc00 = nbocv_i[icv];
				adouble dsrcdphi = -2.0 * TS * beta * omega[icv];
				A[noc00] -= dsrcdphi.value() * cv_volume[icv];
			}
		}
	}
#endif

	/*
	 * Method: sourceHookRansTurbCoupled_AD
	 * ------------------------------------
	 *
	 */
#ifdef USE_MEM_SAVING_ADVAR
	virtual void sourceHookRansTurbCoupled1D_AD(const int icv, adouble *rhs, double ***A, int flagImplicit) {
		int kine_Index = getScalarTransportIndex("kine");
		int omega_Index = getScalarTransportIndex("omega");

		//
		{
			adouble src = calcTurbProd_AD(icv) - betaStar * rho_AD[icv] * (*omega)[icv] * (*kine)[icv];
			rhs[5 + kine_Index] += src * cv_volume[icv];

			if (flagImplicit && icv < ncv) {
				int noc00 = nbocv_i[icv];
				adouble dsrcdphi = -betaStar * (*omega)[icv];
				A[noc00][5 + kine_Index][5 + kine_Index] -= dsrcdphi.value() * cv_volume[icv];
			}
		}

		//
		{
			adouble OM[3][3], STR_hat[3][3];

			adouble TS = 1.0;
			if (KOM_RealizableConstraint == 2)
				TS = min(1.0 / ((*omega)[icv]), 0.6 / (sqrt(6.0) * strMag[icv]));

			adouble sigmad;
			adouble crossDiff = vecDotVec3d_AD((*grad_kine)[icv], (*grad_omega)[icv]);

			if (crossDiff <= 0.0)
				sigmad = 0.0;
			else
				sigmad = sigmad0;

			adouble divU = grad_u[icv][0][0] + grad_u[icv][1][1] + grad_u[icv][2][2];

			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++) {
					OM[i][j] = 0.5 * (grad_u[icv][i][j] - grad_u[icv][j][i]);

					// strange: STR_hat subtracts 0.5 divergence instead 1/3, no bug!!!
					if (i == j)
						STR_hat[i][j] = 0.5 * (grad_u[icv][i][j] + grad_u[icv][j][i] - divU);
					else
						STR_hat[i][j] = 0.5 * (grad_u[icv][i][j] + grad_u[icv][j][i]);
				}

			adouble chiOm = 0.0;
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					for (int k = 0; k < 3; k++)
						chiOm += OM[i][j] * OM[j][k] * STR_hat[k][i];
			chiOm = fabs(chiOm / (betaStar * (*omega)[icv] * betaStar * (*omega)[icv] * betaStar * (*omega)[icv]));

			adouble fbeta = (1.0 + 85.0 * chiOm) / (1.0 + 100.0 * chiOm);
			adouble beta = beta0 * fbeta;

			adouble src = TS * alfa * (*omega)[icv] / (*kine)[icv]
					* calcTurbProd_AD(icv)
					- TS * beta * rho_AD[icv] * (*omega)[icv] * (*omega)[icv]
					+ sigmad * rho_AD[icv] / (*omega)[icv] * crossDiff;

			rhs[5 + omega_Index] += src * cv_volume[icv];

			if (flagImplicit && icv < ncv) {
				int noc00 = nbocv_i[icv];
				adouble dsrcdphi = -2.0 * TS * beta * (*omega)[icv];
				A[noc00][5+omega_Index][5+omega_Index] -= dsrcdphi.value() * cv_volume[icv];
			}
		}
	}
#else
	virtual void sourceHookRansTurbCoupled1D_AD(const int icv, adouble *rhs, double ***A, int flagImplicit) {
		int kine_Index = getScalarTransportIndex("kine");
		int omega_Index = getScalarTransportIndex("omega");

		//
		{
			adouble src = calcTurbProd_AD(icv)
					- betaStar * rho_AD[icv] * omega[icv] * kine[icv];
			rhs[5 + kine_Index] += src * cv_volume[icv];

			if (flagImplicit && icv < ncv) {
				int noc00 = nbocv_i[icv];
				adouble dsrcdphi = -betaStar * omega[icv];
				A[noc00][5 + kine_Index][5 + kine_Index] -= dsrcdphi.value()
						* cv_volume[icv];
			}
		}

		//
		{
			adouble OM[3][3], STR_hat[3][3];

			adouble TS = 1.0;
			if (KOM_RealizableConstraint == 2)
				TS = min(1.0 / (omega[icv]), 0.6 / (sqrt(6.0) * strMag[icv]));

			adouble sigmad;
			adouble crossDiff = vecDotVec3d_AD(grad_kine[icv], grad_omega[icv]);

			if (crossDiff <= 0.0)
				sigmad = 0.0;
			else
				sigmad = sigmad0;

			adouble divU = grad_u[icv][0][0] + grad_u[icv][1][1]
					+ grad_u[icv][2][2];

			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++) {
					OM[i][j] = 0.5 * (grad_u[icv][i][j] - grad_u[icv][j][i]);

					// strange: STR_hat subtracts 0.5 divergence instead 1/3, no bug!!!
					if (i == j)
						STR_hat[i][j] =
								0.5
										* (grad_u[icv][i][j] + grad_u[icv][j][i]
												- divU);
					else
						STR_hat[i][j] = 0.5
								* (grad_u[icv][i][j] + grad_u[icv][j][i]);
				}

			adouble chiOm = 0.0;
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					for (int k = 0; k < 3; k++)
						chiOm += OM[i][j] * OM[j][k] * STR_hat[k][i];
			chiOm = fabs(
					chiOm
							/ (betaStar * omega[icv] * betaStar * omega[icv]
									* betaStar * omega[icv]));

			adouble fbeta = (1.0 + 85.0 * chiOm) / (1.0 + 100.0 * chiOm);
			adouble beta = beta0 * fbeta;

			adouble src = TS * alfa * omega[icv] / kine[icv]
					* calcTurbProd_AD(icv)
					- TS * beta * rho_AD[icv] * omega[icv] * omega[icv]
					+ sigmad * rho_AD[icv] / omega[icv] * crossDiff;

			rhs[5 + omega_Index] += src * cv_volume[icv];

			if (flagImplicit && icv < ncv) {
				int noc00 = nbocv_i[icv];
				adouble dsrcdphi = -2.0 * TS * beta * omega[icv];
				A[noc00][5+omega_Index][5+omega_Index] -= dsrcdphi.value() * cv_volume[icv];
			}
		}
	}
#endif

	/*
	 * Method: boundaryHookScalarRansTurb1D
	 * ------------------------------------
	 * Original code = boundaryHookScalarRansTurb_AD() in IkeUgpWithCvCompFlow.h and boundaryHookScalarRansTurb() in UgpWithCvCompFlowAD.h
	 * In some situation, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
	 */
	void boundaryHookScalarRansTurb1D(const int ifa, double *phi, FaZone *zone, const string &scalName) {
		static bool firstCall = true;
		int debugLevel = getDebugLevel();

		if (zone->getKind() == FA_ZONE_BOUNDARY) {
			Param *param;
			if (getParam(param, zone->getName())) {
				if ((param->getString() == "WALL") && (scalName == "omega")) {
					if(debugLevel>0) {
						if(mpi_rank==0 && firstCall)
							cout<<"Calling boundaryHookScalarRansTurb1D() for boundary "<<zone->getName()<<endl;
						firstCall = false;
					}

					int icv0 = cvofa[ifa][0];
					int icv1 = cvofa[ifa][1];
					double muLamCV = calcMuLam(UgpWithCvCompFlow::temp[icv0]);
					phi[icv1] = 6.0*muLamCV/(rho[icv0]*beta0*wallDist[icv0]*wallDist[icv0]);
				}
			}
		}
	}

	/*
	 * Method: boundaryHookScalarRansTurb1D_AD
	 * ---------------------------------------
	 *
	 */
#ifdef USE_MEM_SAVING_ADVAR
	virtual void boundaryHookScalarRansTurb1D_AD(const int ifa, ADscalar<adouble> &phi, FaZone *zone, const string &name) {
		static bool firstCall = true;
		int debugLevel = getDebugLevel();

		if (zone->getKind() == FA_ZONE_BOUNDARY) {
			Param *param;
			if (getParam(param, zone->getName()))
				if ((param->getString() == "WALL") && (name == "omega")) {
					if(debugLevel>0) {
						if(mpi_rank==0 && firstCall)
							cout<<"Calling boundaryHookScalarRansTurb1D_AD() for boundary "<<zone->getName()<<endl;
						firstCall = false;
					}

					int icv0 = cvofa[ifa][0];
					int icv1 = cvofa[ifa][1];
					phi[icv1] = 6.0*calcMuLam_AD(icv0) / (rho_AD[icv0]*beta0*wallDist[icv0]*wallDist[icv0]);
				}
		}
	}
#endif
	virtual void boundaryHookScalarRansTurb1D_AD(const int ifa, adouble *phi, FaZone *zone, const string &name) {
		static bool firstCall = true;
		int debugLevel = getDebugLevel();

		if (zone->getKind() == FA_ZONE_BOUNDARY) {
			Param *param;
			if (getParam(param, zone->getName()))
				if ((param->getString() == "WALL") && (name == "omega")) {
					if(debugLevel>0) {
						if(mpi_rank==0 && firstCall)
							cout<<"Calling boundaryHookScalarRansTurb1D_AD() for boundary "<<zone->getName()<<endl;
						firstCall = false;
					}

					int icv0 = cvofa[ifa][0];
					int icv1 = cvofa[ifa][1];
					phi[icv1] = 6.0*calcMuLam_AD(icv0) / (rho_AD[icv0]*beta0*wallDist[icv0]*wallDist[icv0]);
				}
		}
	}
};

#endif /* IKETURBMODEL_KOM_H_ */
