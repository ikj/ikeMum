#include "JoeWithModels.h"
#include "JOE/ADJOINT_FILES/JoeWithModelsAD.h"
#include "adolc.h"

#include "IkeTurbModel_KOM.h"
#include "turbModels/TurbModel_KOM.h"
#include "ADJOINT_FILES/TurbModel_KOM_AD.h"
//#include "turbModels/TurbModel_KOMSST.h"
//#include "ADJOINT_FILES/TurbModel_KOMSST_AD.h"

#include "IkeWithPsALC.h"

#include <cstdlib>
#include <ctime>
#include <math.h>

#include <iomanip>

#ifndef PI
#define PI 3.14159265359
#endif

// ###########################################################################################
// ------                                                                               ------
// ------                      JOE with k-omega called JOEhyShot2D                      ------
// ------                                                                               ------
// ###########################################################################################
#define NVAL 7
#define VERY_SMALL_DOUBLE 1.0e-10

class JOEhyShot2DkOm : public JoeWithModels, public RansTurbKOm
{
protected:
//	double gammaRef, RoMRef;
//	double heatMaxX;

	int NcontrolEqns;
	double *lambda;

	// Heat Release Model: Based on C.Doolan & R.Boyce, AIAA 2008-2603
	double fSt; //Stoichiometric fuel/air ratio (0.28 for H2/air)
	double hFuel; //Fuel heating value (H2)
	double mAir; //Mass flow rate of injected air

	double lengthComb; //Combustor length
	double kComb; //Free parameter: Fraction of completed combustion (BLR exp.: 95%)
	double dComb; //Free parameter: Shape of heat release

	double xCombStart; //Combustion ignition starting position
	double xThroat, xMax;
	double yBottom, yTop, yMax;
	double width;

	// Statistics
	double mySumOmegaWall_JOE;
	double myMaxOmegaWall_JOE;
	double mySumHeatRelease_JOE;

public:
	JOEhyShot2DkOm(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name) {
		if (mpi_rank == 0)
			cout << "JOEhyShot2DkOm()" << endl;

//		/* Flow parameter - gammaRef and RoMRef */
//		gammaRef = getDoubleParam("GAMMA", "1.4");
//		double pressRef = getDoubleParam("P_REF", "1.4");
//		double tempRef = getDoubleParam("T_REF", "1.0");
//		double rhoRef = getDoubleParam("RHO_REF", "1.0");
//		RoMRef = pressRef/rhoRef/tempRef;

		/* Heat release */
		NcontrolEqns = getIntParam("N_PSALC_CONTROL_PARAMS", "0");
		if(NcontrolEqns > 0) {
			lambda = new double [NcontrolEqns]; // lambda is the input parameter to the system (e.g. fuel equivalence ratio)
			char paramName[20];
			for(int i=0; i<NcontrolEqns; ++i) {
				sprintf(paramName, "CONTROL_PARAM%d", i);
				lambda[i] = getParam("LAMBDA_INITIAL_0") -> getDouble(paramName); // Note: lambda is same as phi in the model
			}
		} else {
			lambda = NULL;
		}

		fSt   = 0.028; //Stoichiometric fuel/air ratio (0.28 for H2/air)
		hFuel = 1.20e8; //Fuel heating value (120MJ/kg for H2)
		mAir  = 0.053*20; //Mass flow rate of injected air

		kComb = 0.95; //Free parameter: Fraction of completed combustion (DLR simulation result: 95%)
		dComb = getDoubleParam("SHAPE_HEAT_RELEASE",0.75); //Free parameter: Shape of heat release (best result from Vince's experiment - he tried dComb=0.5,0.75,1.0,1.25,1.5)

		xCombStart = 0.408; //Combustion ignition starting position

		xThroat  = 0.64983;
		xMax     = 0.773; //0.775856;
		yBottom  = 0.1195;
		yTop     = 0.1293;
		yMax     = 0.155571; //0.156088;

		lengthComb = xMax - xCombStart; //Combustor length
		width      = 0.2; //Combustor width

		// Statistics
		mySumOmegaWall_JOE   = 0.0;
		mySumHeatRelease_JOE = 0.0;
		myMaxOmegaWall_JOE   = -1.0;

		if(mpi_rank==0) {
			cout<<endl<<"**** SCRAMJET HEAT RELEASE MODEL: FROM C.Doolan & R.Boyce, AIAA 2008-2603"<<endl;
			if(NcontrolEqns > 0) {
				printf("**    lambda (phi) = ");
				for(int i=0; i<NcontrolEqns; ++i)
					printf("%.3e ", lambda[i]);
				cout<<endl;
			} else {
				printf("**    lambda (phi) = EMPTY \n");
			}
			cout<<"**    SHAPE_HEAT_RELEASE = "<<dComb<<endl;
			cout<<"**    xCombStart         = "<<xCombStart<<endl<<endl;
		}
	}

	/*
	 * destructor
	 */
	virtual ~JOEhyShot2DkOm() {
		if (mpi_rank == 0)
			cout << "~JOEhyShot2DkOm()" << endl;
	}

	/**************************
	 *---  HOOK FUNCTIONS  ---*
	 **************************/
	/*
	 * Method: initHookPsALC_1st
	 * -------------------------
	 *
	 */
	virtual void initialHook() {
		if(mpi_rank == 0)
			cout<<"initialHook(): ";
		if(!checkScalarFlag("RHO")) {
			JoeWithModels::initialHook();

			if(mpi_rank == 0)
				cout<<"  initialize from reference values "<<endl;
			double u_init[3] ={2313.0, 0.0, 0.0};

			for (int icv=0; icv<ncv; ++icv) {
				rho[icv] = rho_ref;
				rhou[icv][0] = rho_ref*u_init[0];
				rhou[icv][1] = rho_ref*u_init[1];
				rhou[icv][2] = rho_ref*u_init[2];
				rhoE[icv] = p_ref/(gamma[icv]-1.0)+ 0.5*rho_ref*vecDotVec3d(u_init, u_init);
			}

			double *kine = getR1("kine");
			double *omega = getR1("omega");

			for(int icv=0; icv<ncv; ++icv) {
				kine[icv] = 32.0;
				omega[icv] = 1.0e4;
			}

			updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
			updateCvDataG1G2(rho,  REPLACE_DATA);
			updateCvDataG1G2(rhoE, REPLACE_DATA);
			updateCvDataG1G2(kine, REPLACE_DATA);
			updateCvDataG1G2(omega, REPLACE_DATA);

			if(mpi_rank == 0)
				cout<<" initialize by constant profiles "<<endl;
		} else {
			if(mpi_rank == 0)
				cout<<" initialize from a restart file "<<endl;
		}
	}

	/*
	 * Method: boundaryHook
	 * --------------------
	 *
	 */
	virtual void boundaryHook(double *temp, double (*vel)[3], double *press, FaZone *zone) {
		/* empty */
	}

	/*
	 * Method: boundaryHook1D
	 * -------------------------
	 *
	 */
	virtual void boundaryHook1D(const int ifa, double *T_fa, double (*vel_fa)[3], double *p_fa, FaZone *zone) {
		/* empty */
	}

	/*
	 * Method: boundaryHookScalarRansTurb
	 * ----------------------------------
	 *
	 */
	void boundaryHookScalarRansTurb(double *phi, FaZone *zone, const string &name) {
		if (zone->getKind() == FA_ZONE_BOUNDARY) {
			Param *param;
			if (getParam(param, zone->getName())) {
				if ((param->getString() == "WALL") && (name == "omega")) {
					for (int index = 0; index < zone->faVec.size(); ++index) {
						int ifa = zone->faVec[index];
						int icv0 = cvofa[ifa][0];
						int icv1 = cvofa[ifa][1];
						double muLamCV = calcMuLam(UgpWithCvCompFlow::temp[icv0]);
						double omegaNew = 6.0*muLamCV/(rho[icv0]*beta0*wallDist[icv0]*wallDist[icv0]);

						phi[icv1] = omegaNew;
						mySumOmegaWall_JOE += omegaNew;
						myMaxOmegaWall_JOE = max(myMaxOmegaWall_JOE, fabs(omegaNew));
					}
				}
			}
		}
	}

	/*
	 * Method: boundaryHookScalarRansTurb1D
	 * ------------------------------------
	 * Original code = boundaryHookScalarRansTurb_AD() in IkeUgpWithCvCompFlow.h and boundaryHookScalarRansTurb() in UgpWithCvCompFlowAD.h
	 * In some situation, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
	 */
	void boundaryHookScalarRansTurb1D(const int ifa, double *phi, FaZone *zone, const string &scalName) {
		if (zone->getKind() == FA_ZONE_BOUNDARY) {
			Param *param;
			if (getParam(param, zone->getName())) {
				if ((param->getString() == "WALL") && (scalName == "omega")) {
					int icv0 = cvofa[ifa][0];
					int icv1 = cvofa[ifa][1];
					double muLamCV = calcMuLam(UgpWithCvCompFlow::temp[icv0]);
					double omegaNew = 6.0*muLamCV/(rho[icv0]*beta0*wallDist[icv0]*wallDist[icv0]);

					phi[icv1] = omegaNew;
					mySumOmegaWall_JOE += omegaNew;
					myMaxOmegaWall_JOE = max(myMaxOmegaWall_JOE, fabs(omegaNew));
				}
			}
		}
	}

	/*
	 * Method: sourceHook
	 * ------------------
	 *
	 */
	void sourceHook(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5]) {
		double constant = lambda[0]*fSt*mAir*hFuel;

		for (int icv=0; icv<ncv; ++icv) {
			double xCoord = x_cv[icv][0];
			double yCoord = x_cv[icv][1];

			double area = 0.0;
			if(xCoord < xThroat) {
				if(xCoord > xCombStart)
					area = (yTop-yBottom)*width;
			} else {
				area = ((yTop-yBottom) + (yMax-yTop)/(xMax-xThroat)*(xCoord-xThroat))*width;
			}

			if (xCoord >= xCombStart && yCoord > yBottom && yCoord < yBottom+area/width) {
				double xN = (xCoord-xCombStart)/lengthComb;
				double func = dComb*kComb/lengthComb*pow(kComb*xN,dComb-1)*exp(-pow(kComb*xN,dComb));
				rhs_rhoE[icv] += constant/area*func*cv_volume[icv];

				mySumHeatRelease_JOE += constant/area*func*cv_volume[icv];
			}
		}
	}

	/*
	 * Method: temporalHook
	 * --------------------
	 *
	 */
	void temporalHook() {
		if(step%100 == 0) {
			double totSumOmegaWall_JOE;
			MPI_Allreduce(&mySumOmegaWall_JOE, &totSumOmegaWall_JOE, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
			double totMaxOmegaWall_JOE;
			MPI_Allreduce(&myMaxOmegaWall_JOE, &totMaxOmegaWall_JOE, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);

			double totSumHeatRelease_JOE;
			MPI_Allreduce(&mySumHeatRelease_JOE, &totSumHeatRelease_JOE, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

			if(mpi_rank==0) {
				cout<<"JOEhyShot2DkOm::temporalHook()"<<endl;
				printf("    TOTAL OMEGA FROM WALL-FUNC (JOE): %.5e\n", totSumOmegaWall_JOE);
				printf("    MAX   OMEGA FROM WALL-FUNC (JOE): %.5e\n", totMaxOmegaWall_JOE);
				printf("    TOTAL HEAT-RELEASE FROM WALL-FUNC (JOE): %.5e\n", totSumHeatRelease_JOE);
			}

			mySumOmegaWall_JOE = 0.0;
			myMaxOmegaWall_JOE = -1.0;
			mySumHeatRelease_JOE = 0.0;
		}
	}

	/*
	 * Method: finalHook
	 * -----------------
	 *
	 */
	void finalHook() {
	}
};


// ###########################################################################################
// ------                                                                               ------
// ------                  2D Nozzle with PsALC called PsALChyShot2D                    ------
// ------                                                                               ------
// ###########################################################################################
/* Note: Multiple inheritance occurs here. The problem of ambiguous hierarchy composition is avoided by using virtual inheritance
 *       (In the diagram below, doubled lines are virtual public, single lines are public)
 *
 *                        UgpWithCvCompFlow
 *                  //           ||           \\
 *     JoeWithModels    UgpWithCvCompFlow_AD    RansTurbKOm
 *              \\       //              \\       /?
 *        JoeWithModels_AD               RansTurbKOm_AD
 *                     ||                  |
 *        IkeWithModels_AD               IkeRansTurbKOm_AD
 *                      |                  |
 *        IkeWithPsALC_AD                  |
 *                       \                /
 *                      endUserClass(ike.cpp)
 */
class PsALChyShot2D : public IkeWithPsALC_AD, public IkeRansTurbKOm_AD
{
protected:
//	double gammaRef, RoMRef;
//	double heatMaxX;

	// Heat Release Model: Based on C.Doolan & R.Boyce, AIAA 2008-2603
	double fSt; //Stoichiometric fuel/air ratio (0.28 for H2/air)
	double hFuel; //Fuel heating value (H2)
	double mAir; //Mass flow rate of injected air

	double lengthComb; //Combustor length
	double kComb; //Free parameter: Fraction of completed combustion (BLR exp.: 95%)
	double dComb; //Free parameter: Shape of heat release

	double xCombStart; //Combustion ignition starting position
	double xThroat, xMax;
	double yBottom, yTop, yMax;
	double width;

	// Statistics
	int    myCountOmegaWall_AD;
	double mySumOmegaWall_JOE;
	double mySumOmegaWall_AD;
	double mySumHeatRelease_JOE;
	double mySumHeatRelease_AD;

public:
	PsALChyShot2D(char *name) : IkeWithPsALC_AD(name), JoeWithModels(name), JoeWithModels_AD(name), UgpWithCvCompFlow(name) {
		if (mpi_rank == 0)
			cout << "PsALChyShot2D()" << endl;

//		/* Flow parameter - gammaRef and RoMRef */
//		gammaRef = getDoubleParam("GAMMA", "1.4");
//		double pressRef = getDoubleParam("P_REF", "1.4");
//		double tempRef = getDoubleParam("T_REF", "1.0");
//		double rhoRef = getDoubleParam("RHO_REF", "1.0");
//		RoMRef = pressRef/rhoRef/tempRef;

		fSt   = 0.028; //Stoichiometric fuel/air ratio (0.28 for H2/air)
		hFuel = 1.20e8; //Fuel heating value (120MJ/kg for H2)
		mAir  = 0.053*20; //Mass flow rate of injected air

		kComb = 0.95; //Free parameter: Fraction of completed combustion (DLR simulation result: 95%)
		dComb = getDoubleParam("SHAPE_HEAT_RELEASE",0.75); //Free parameter: Shape of heat release (best result from Vince's experiment - he tried dComb=0.5,0.75,1.0,1.25,1.5)

		xCombStart = 0.408; //Combustion ignition starting position

		xThroat  = 0.64983;
		xMax     = 0.773; //0.775856;
		yBottom  = 0.1195;
		yTop     = 0.1293;
		yMax     = 0.155571; //0.156088;

		lengthComb = xMax - xCombStart; //Combustor length
		width      = 0.2; //Combustor width

		if(mpi_rank==0) {
			cout<<endl<<"**** SCRAMJET HEAT RELEASE MODEL: FROM C.Doolan & R.Boyce, AIAA 2008-2603"<<endl;
			if(NcontrolEqns > 0) {
				printf("**    lambda (phi) = ");
				for(int i=0; i<NcontrolEqns; ++i)
					printf("%.3e ", lambda[i]);
				cout<<endl;
			} else {
				printf("**    lambda (phi) = EMPTY \n");
			}
			cout<<"**    SHAPE_HEAT_RELEASE = "<<dComb<<endl;
			cout<<"**    xCombStart         = "<<xCombStart<<endl<<endl;
		}

		// Statistics
		myCountOmegaWall_AD = 0;
		mySumOmegaWall_JOE   = 0.0;
		mySumOmegaWall_AD    = 0.0;
		mySumHeatRelease_JOE = 0.0;
		mySumHeatRelease_AD  = 0.0;
	}

	/*
	 * destructor
	 */
	virtual ~PsALChyShot2D() {
		if (mpi_rank == 0)
			cout << "~PsALChyShot2D()" << endl;
	}

	/**************************
	 *---  HOOK FUNCTIONS  ---*
	 **************************/
	/*
	 * Method: initHookPsALC_1st
	 * -------------------------
	 *
	 */
	virtual void initHookPsALC_1st() {
		if(mpi_rank == 0)
			cout<<"initHookPsALC_1st(): ";
		if(!checkScalarFlag("RHO")) {
			cerr<<"ERROR!!! Cannot detect the initial field"<<endl;
			throw(-328);
		} else {
			if(mpi_rank == 0)
				cout<<"  initialize from a restart file "<<endl;
		}
	}

	/*
	 * Method: boundaryHook
	 * --------------------
	 *
	 */
	virtual void boundaryHook(double *temp, double (*vel)[3], double *press, FaZone *zone) {
		/* empty */
	}

	/*
	 * Method: boundaryHook1D
	 * -------------------------
	 *
	 */
	virtual void boundaryHook1D(const int ifa, double *T_fa, double (*vel_fa)[3], double *p_fa, FaZone *zone) {
		/* empty */
	}

//	/*
//	 * Method: boundaryHook_AD
//	 * -----------------------
//	 *
//	 */
//	virtual void boundaryHook_AD(REALQ *T_fa, REALQ (*vel_fa)[3], REALQ *p_fa, FaZone *zone) {
//		/* empty */
//	}

//	/*
//	 * Method: boundaryHook1D_AD
//	 * -------------------------
//	 *
//	 */
//	virtual void boundaryHook1D_AD(const int ifa, REALQ *T_fa, REALQ (*vel_fa)[3], REALQ *p_fa, FaZone *zone) {
//		/* empty */
//	}

	/*
	 * Method: boundaryHookScalarRansTurb
	 * ----------------------------------
	 *
	 */
	void boundaryHookScalarRansTurb(double *phi, FaZone *zone, const string &name) {
		if (zone->getKind() == FA_ZONE_BOUNDARY) {
			Param *param;
			if (getParam(param, zone->getName())) {
				if ((param->getString() == "WALL") && (name == "omega")) {
					for (int index = 0; index < zone->faVec.size(); ++index) {
						int ifa = zone->faVec[index];
						int icv0 = cvofa[ifa][0];
						int icv1 = cvofa[ifa][1];
						double muLamCV = calcMuLam(UgpWithCvCompFlow::temp[icv0]);
						double omegaNew = 6.0*muLamCV/(rho[icv0]*beta0*wallDist[icv0]*wallDist[icv0]);

						phi[icv1] = omegaNew;
						mySumOmegaWall_JOE += omegaNew;
					}
				}
			}
		}
	}

	/*
	 * Method: boundaryHookScalarRansTurb1D
	 * ------------------------------------
	 * Original code = boundaryHookScalarRansTurb_AD() in IkeUgpWithCvCompFlow.h and boundaryHookScalarRansTurb() in UgpWithCvCompFlowAD.h
	 * In some situation, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
	 */
	void boundaryHookScalarRansTurb1D(const int ifa, double *phi, FaZone *zone, const string &scalName) {
		if (zone->getKind() == FA_ZONE_BOUNDARY) {
			Param *param;
			if (getParam(param, zone->getName())) {
				if ((param->getString() == "WALL") && (scalName == "omega")) {
					int icv0 = cvofa[ifa][0];
					int icv1 = cvofa[ifa][1];
					double muLamCV = calcMuLam(UgpWithCvCompFlow::temp[icv0]);
					double omegaNew = 6.0*muLamCV/(rho[icv0]*beta0*wallDist[icv0]*wallDist[icv0]);
					phi[icv1] = omegaNew;
					mySumOmegaWall_JOE += omegaNew;
					++myCountOmegaWall_AD;
				}
			}
		}
	}

//	/*
//	 * Method: boundaryHookScalarRansTurb_AD
//	 * -------------------------------------
//	 *
//	 */
//	void boundaryHookScalarRansTurb_AD(adouble *phi, FaZone *zone, const string &scalName) {
//			if (zone->getKind() == FA_ZONE_BOUNDARY) {
//			Param *param;
//			if (getParam(param, zone->getName())) {
//				if (strcmp(param->getString().c_str(),"WALL")==0 && strcmp(scalName.c_str(),"omega")==0) {
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//						int icv0 = cvofa[ifa][0];
//						int icv1 = cvofa[ifa][1];
//						adouble muLamCV_AD = calcMuLam_AD(temp[icv0]);
//						adouble omegaNew_AD = 6.0*muLamCV_AD/(rho_AD[icv0]*beta0*wallDist[icv0]*wallDist[icv0]);
//						phi[icv1] = omegaNew_AD;
//						mySumOmegaWall_AD += omegaNew_AD.value();
//						++myCountOmegaWall_AD;
//					}
//				}
//			}
//		}
//	}

	/*
	 * Method: boundaryHookScalarRansTurb1D_AD
	 * ---------------------------------------
	 * Acutally this function doesn't need to be overloaded: Currently just a debugging purpose
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
							cout<<"Calling PsALChyShot2D::boundaryHookScalarRansTurb1D_AD() for boundary "<<zone->getName()<<endl;
						firstCall = false;
					}

					int icv0 = cvofa[ifa][0];
					int icv1 = cvofa[ifa][1];
					adouble omegaNew_AD = 6.0*calcMuLam_AD(icv0) / (rho_AD[icv0]*beta0*wallDist[icv0]*wallDist[icv0]);
					phi[icv1] = omegaNew_AD;

					mySumOmegaWall_AD += omegaNew_AD.value();
				}
		}
	}
#endif

	/*
	 * Method: boundaryHookScalarRansTurb1D_AD
	 * ---------------------------------------
	 * This method is for old IKE calculations: it will NOT be called if USE_MEM_SAVING_ADVAR is used
	 */
	virtual void boundaryHookScalarRansTurb1D_AD(const int ifa, adouble *phi, FaZone *zone, const string &scalName)  {
		if (zone->getKind() == FA_ZONE_BOUNDARY) {
			Param *param;
			if (getParam(param, zone->getName())) {
				if (strcmp(param->getString().c_str(),"WALL")==0 && strcmp(scalName.c_str(),"omega")==0) {
					int icv0 = cvofa[ifa][0];
					int icv1 = cvofa[ifa][1];
					adouble muLamCV_AD = calcMuLam_AD(icv0);
					adouble omegaNew_AD = 6.0*muLamCV_AD/(rho_AD[icv0]*beta0*wallDist[icv0]*wallDist[icv0]);
					phi[icv1] = omegaNew_AD;

					mySumOmegaWall_AD += omegaNew_AD.value();
				}
			}
		}
	}

	/*
	 * Method: sourceHook
	 * ------------------
	 *
	 */
	void sourceHook(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5]) {
		if(lambda==NULL) {
			/* nothing by now */
		} else {
			double constant = lambda[0]*fSt*mAir*hFuel;

			for (int icv=0; icv<ncv; ++icv) {
				double xCoord = x_cv[icv][0];
				double yCoord = x_cv[icv][1];

				double area = 0.0;
				if(xCoord < xThroat) {
					if(xCoord > xCombStart)
						area = (yTop-yBottom)*width;
				} else {
					area = ((yTop-yBottom) + (yMax-yTop)/(xMax-xThroat)*(xCoord-xThroat))*width;
				}

				if (xCoord >= xCombStart && yCoord > yBottom && yCoord < yBottom+area/width) {
					double xN = (xCoord-xCombStart)/lengthComb;
					double func = dComb*kComb/lengthComb*pow(kComb*xN,dComb-1)*exp(-pow(kComb*xN,dComb));
					rhs_rhoE[icv] += constant/area*func*cv_volume[icv];
					mySumHeatRelease_JOE += constant/area*func*cv_volume[icv];
				}
			}
		}
	}

//	/*
//	 * Method: sourceHook_AD
//	 * ---------------------
//	 *
//	 */
//	void sourceHook_AD(REALA *rhs_rho_AD, REALA (*rhs_rhou_AD)[3], REALA *rhs_rhoE_AD, double (*A)[5][5]) {
//		if(lambda_AD==NULL) { // If the Jacobian doesn't need to produce the dRHSdLambda terms (i.e. Newton method for the second point)
//			                  // (Note that NcontrolEqns is always non-zero; so never use something like "if(NcontrolEqns==0)")
//			double constant = lambda[0]*fSt*mAir*hFuel;
//
//			for (int icv=0; icv<ncv; ++icv) {
//				double xCoord = x_cv[icv][0];
//				double yCoord = x_cv[icv][1];
//
//				double area = 0.0;
//				if(xCoord < xThroat) {
//					if(xCoord > xCombStart)
//						area = (yTop-yBottom)*width;
//				} else {
//					area = ((yTop-yBottom) + (yMax-yTop)/(xMax-xThroat)*(xCoord-xThroat))*width;
//				}
//
//				if (xCoord >= xCombStart && yCoord > yBottom && yCoord < yBottom+area/width) {
//					double xN = (xCoord-xCombStart)/lengthComb;
//					double func = dComb*kComb/lengthComb*pow(kComb*xN,dComb-1)*exp(-pow(kComb*xN,dComb));
//					rhs_rhoE_AD[icv] += constant/area*func*cv_volume[icv];
//
//					mySumHeatRelease_AD += constant/area*func*cv_volume[icv];
//				}
//			}
//		} else {
//			adouble constant = lambda_AD[0]*fSt*mAir*hFuel;
//
//			for (int icv=0; icv<ncv; ++icv) {
//				double xCoord = x_cv[icv][0];
//				double yCoord = x_cv[icv][1];
//
//				double area = 0.0;
//				if(xCoord < xThroat) {
//					if(xCoord > xCombStart)
//						area = (yTop-yBottom)*width;
//				} else {
//					area = ((yTop-yBottom) + (yMax-yTop)/(xMax-xThroat)*(xCoord-xThroat))*width;
//				}
//
//				if (xCoord >= xCombStart && yCoord > yBottom && yCoord < yBottom+area/width) {
//					double xN = (xCoord-xCombStart)/lengthComb;
//					double func = dComb*kComb/lengthComb*pow(kComb*xN,dComb-1)*exp(-pow(kComb*xN,dComb));
//					rhs_rhoE_AD[icv] += constant/area*func*cv_volume[icv];
//					mySumHeatRelease_AD += constant.value()/area*func*cv_volume[icv];
//				}
//			}
//		}
//	}

	/*
	 * Method: sourceHook1D_AD
	 * -----------------------
	 * Original code = sourceHook_AD() in JoeWithModelsAD.h
	 */
	virtual void sourceHook1D_AD(const int icvCenter, REALA &rhs_rho_AD, REALA rhs_rhou_AD[3], REALA &rhs_rhoE_AD, double (*A)[5][5]) {
		if(lambda_AD==NULL) { // If the Jacobian doesn't need to produce the dRHSdLambda terms (i.e. Newton method for the second point)
		                      // (Note that NcontrolEqns is always non-zero; so never use something like "if(NcontrolEqns==0)")
			assert(lambda!=NULL);
			double constant = lambda[0]*fSt*mAir*hFuel;

			double xCoord = x_cv[icvCenter][0];
			double yCoord = x_cv[icvCenter][1];

			double area = 0.0;
			if(xCoord < xThroat) {
				if(xCoord > xCombStart)
					area = (yTop-yBottom)*width;
			} else {
				area = ((yTop-yBottom) + (yMax-yTop)/(xMax-xThroat)*(xCoord-xThroat))*width;
			}

			if (xCoord >= xCombStart && yCoord > yBottom && yCoord < yBottom+area/width) {
				double xN = (xCoord-xCombStart)/lengthComb;
				double func = dComb*kComb/lengthComb*pow(kComb*xN,dComb-1)*exp(-pow(kComb*xN,dComb));
				rhs_rhoE_AD += constant/area*func*cv_volume[icvCenter];
				mySumHeatRelease_AD += constant/area*func*cv_volume[icvCenter];
			}
		} else {
			adouble constant = lambda_AD[0]*fSt*mAir*hFuel;

			double xCoord = x_cv[icvCenter][0];
			double yCoord = x_cv[icvCenter][1];

			double area = 0.0;
			if(xCoord < xThroat) {
				if(xCoord > xCombStart)
					area = (yTop-yBottom)*width;
			} else {
				area = ((yTop-yBottom) + (yMax-yTop)/(xMax-xThroat)*(xCoord-xThroat))*width;
			}

			if (xCoord >= xCombStart && yCoord > yBottom && yCoord < yBottom+area/width) {
				double xN = (xCoord-xCombStart)/lengthComb;
				double func = dComb*kComb/lengthComb*pow(kComb*xN,dComb-1)*exp(-pow(kComb*xN,dComb));
				rhs_rhoE_AD += constant/area*func*cv_volume[icvCenter];
				mySumHeatRelease_AD += constant.value()/area*func*cv_volume[icvCenter];
			}
		}
	}

	/*
	 * Method: temporalHook
	 * --------------------
	 * Basically, this method is not called in AD
	 */
	virtual void temporalHook() {
		double totSumOmegaWall_JOE, totSumOmegaWall_AD;
		MPI_Allreduce(&mySumOmegaWall_JOE, &totSumOmegaWall_JOE, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(&mySumOmegaWall_AD, &totSumOmegaWall_AD, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

		double totSumHeatRelease_JOE, totSumHeatRelease_AD;
		MPI_Allreduce(&mySumHeatRelease_JOE, &totSumHeatRelease_JOE, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(&mySumHeatRelease_AD,  &totSumHeatRelease_AD,  1, MPI_DOUBLE, MPI_SUM, mpi_comm);

		if(mpi_rank==0) {
			cout<<"PsALChyShot2D::temporalHook() --";
			if(fabs(totSumOmegaWall_JOE)>MACHINE_EPS)
				printf("  TOTAL OMEGA FROM WALL-FUNC (JOE)= %.5e", totSumOmegaWall_JOE);
			if(fabs(totSumOmegaWall_AD)>MACHINE_EPS)
				printf("  TOTAL OMEGA FROM WALL-FUNC (AD) = %.5e", totSumOmegaWall_AD);

			if(fabs(totSumHeatRelease_JOE)>MACHINE_EPS)
				printf("   TOTAL HEAT-RELEASE (JOE)= %.5e", totSumHeatRelease_JOE);
			if(fabs(totSumHeatRelease_AD)>MACHINE_EPS)
				printf("   TOTAL HEAT-RELEASE (AD) = %.5e", totSumHeatRelease_AD);
			cout<<endl;
		}

		mySumOmegaWall_JOE = 0.0;
		mySumOmegaWall_AD  = 0.0;
		mySumHeatRelease_JOE = 0.0;
		mySumHeatRelease_AD  = 0.0;
	}

	/*
	 * Method: probeInCalcJacobian1DAD
	 * -------------------------------
	 * This method will be called in the end of calcJacobian1DAD() -- just before the memory-releasing step.
	 * If the user wants to probe a quantity (usually something related to the RHS calculation) at every Newton-Raphson iteration,
	 * she/he may want to call this method.
	 */
	void probeInCalcJacobian1DAD(MatComprsedSTL &jacMatrixSTL, double *rhsSingleArray, const int debugLevel, const int NcontrolEqns) {
		int totCountOmegaWall_AD;
		MPI_Allreduce(&myCountOmegaWall_AD, &totCountOmegaWall_AD, 1, MPI_INT, MPI_SUM, mpi_comm);

		double totSumOmegaWall_AD;
		MPI_Allreduce(&mySumOmegaWall_AD, &totSumOmegaWall_AD, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

		double totSumHeatRelease_AD;
		MPI_Allreduce(&mySumHeatRelease_AD,  &totSumHeatRelease_AD,  1, MPI_DOUBLE, MPI_SUM, mpi_comm);

		if(mpi_rank==0) {
			cout<<"PsALChyShot2D::probeInCalcJacobian1DAD()"<<endl;
			printf("    TOTAL OMEGA FROM WALL-FUNC (AD) : %.5e  from %d CVs\n", totSumOmegaWall_AD, totCountOmegaWall_AD);
			printf("    TOTAL HEAT-RELEASE(AD      )    : %.5e\n", totSumHeatRelease_AD);
		}

		myCountOmegaWall_AD = 0;
		mySumOmegaWall_AD   = 0.0;
		mySumHeatRelease_AD = 0.0;
	}

	/*
	 * Method: finalHook
	 * -----------------
	 *
	 */
	void finalHook() {
	}
};

// ###########################################################################################
// ------                                                                               ------
// ------                                      MAIN                                     ------
// ------                                                                               ------
// ###########################################################################################

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);

	initMpiStuff();

	// the run number specifies the class which is going to be instantiated
	int run = 0;

	// set input name to default "Ike.in"
	char inputFileName[50];
	sprintf(inputFileName, "Ike.in");

	for (int i=1; i<argc; i++) {
		string str(argv[i]);
		if (from_string<int>(run, str, std::dec)) {
			if (mpi_rank == 0)
				cout << "You have specified run number = " << run << endl;
		}
		else
			strcpy(inputFileName, argv[i]);
	}

	if (mpi_rank == 0) {
		cout << "SPECIFIED INPUT NAME = " << inputFileName << endl;
		cout << "SPECIFIED RUN = " << run << endl;
	}


	try {
		// provide total runtime
		double wtime, wtime0;
		MPI_Barrier(mpi_comm);
		if (mpi_rank == 0)
			wtime = MPI_Wtime();

		// declare pointer to JoeWithModels_AD
		JoeWithModels *ike;

		switch (run) {
			case 1:  ike = new JOEhyShot2DkOm(inputFileName); 	break;
			case 11: ike = new PsALChyShot2D(inputFileName); break;
			default:
				if (mpi_rank == 0)
					cerr << "ERROR: run number not available!" << endl;
				throw(-1);
		}

		// run ike
		ike->run();

		MPI_Barrier(mpi_comm);
		if (mpi_rank == 0) {
			double wtime0 = wtime;
			wtime = MPI_Wtime();
			cout << " > total runtime [s]: " << wtime - wtime0 << endl;
		}

		// delete ike (make sure memory is deallocated in destructors
		delete ike;
	}
	catch (int e) {
		cerr << "Exception: " << e << endl;
		MPI_Finalize();
		return(-1);
	}
	catch(...) {
		cerr << "unhandled Exception.\n" << endl;
		MPI_Finalize();
		return(-1);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return (0);
}

