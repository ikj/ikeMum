/*
 * vanGogh.h
 *
 *  Created on: Jan 10, 2013
 *      Author: ikj
 */

#ifndef VANGOGH_H_
#define VANGOGH_H_

#include "JoeWithModels.h"
//#include "JOE/ADJOINT_FILES/JoeWithModelsAD.h"
//#include "adolc.h"

//#include "IkeTurbModel_KOM.h"
//#include "turbModels/TurbModel_KOM.h"
//#include "ADJOINT_FILES/TurbModel_KOM_AD.h"
//#include "turbModels/TurbModel_KOMSST.h"
//#include "ADJOINT_FILES/TurbModel_KOMSST_AD.h"

//#include "IkeWithPsALC.h"

#include <cstdlib>
#include <ctime>
#include <math.h>

#include <sstream>
#include <complex>

#include <stdexcept>
#include <iomanip>

#include "vanGoghUtils.h"
#include "DiffFilter.h"

enum RAND_DISTRIB_FUNC {NORMAL_DISTRIB, UNIFORM_DISTRIB};

#define EOF_ERROR_CHECK_CODE 503

#ifndef PI
#define PI 3.14159265359
#endif

#define NVARS_EIGEN 7
#define NEIGENS 30

/*
 * Note:
 * You must make JoeWithModels::run() virtual
 * You must add the writeData2() method to the UgpWithTools class in UgpWithTools.h:
	  void writeData2(const int caseNum, const int step) {
		for (list<WriteData>::iterator wd = writeDataList.begin(); wd!=writeDataList.end(); wd++)
		  wd->write2(caseNum, step, this);
	  }
 * You must make add the following method to the WriteData class in UgpWithTools.h:
	  void write2(const int caseNum, const int step, UgpWithTools1 * ugp) {
		if (step%interval==0) {
		  // set the data flags...
		  ugp->clearDataFlags();
		  for (list<string>::iterator var = varList.begin(); var!=varList.end(); var++) {
			if (var->compare("STATS")==0) ugp->setStatsFlags();
			else ugp->setDataFlag(var->c_str());
		  }

		  // if a geom is specified, then we flag a bunch of cvs and write...
		  if (geom!=NULL) {
			for (int icv = 0; icv<ugp->getNcv(); icv++)
			  ugp->cv_flag[icv] = 0;
			geom->flagCvs(ugp);

			char filename[32];
			sprintf(filename, "%s.%case04d.%06d.plt", name.c_str(), caseNum, step);
			ugp->writeFlaggedCvsTecplot(filename);
		  }
		}
	  }
 * You must make the "Residual" array as a member variable of the JoeWithModels class
  and comment out the "Residual" array in runForwardEuler(), runRK(), runBackwardEuler(), runBackwardEulerCoupled(), and runBDF2();
 */

// ###########################################################################################
// ------                                                                               ------
// ------                       Perturbation tests called VanGogh                       ------
// ------                                                                               ------
// ###########################################################################################
class VanGogh : public JoeWithModels, public DiffFilter {
protected:
	int NcontrolEqns;
	double *lambda;

	//**************
	// HEAT_RELEASE
	//**************
	
	// Heat Release Model: Based on C.Doolan & R.Boyce, AIAA 2008-2603
	double fSt; //Stoichiometric fuel/air ratio (0.28 for H2/air)
	double hFuel; //Fuel heating value (H2)
	double mAir; //Mass flow rate of injected air

	double lengthComb; //Combustor length
	double kComb; //Free parameter: Fraction of completed combustion (BLR exp.: 95%)
	double dComb; //Free parameter: Shape of heat release

	double xCombStart; //Combustion ignition starting position
	double xMin;
	double xMax;
	double xInlet;
	double xThroat;
	double yInlet;

	// Heat release
	double *heatReleaseFunc;

	// Statistics
	double mySumHeatRelease_JOE;
	
	//**************
	// VAN-GOGH
	//**************
	
	// xcvMin & xcvMax (will be used while reading Q1_PT*****.bin)
	double xcvMin[3];
	double xcvMax[3];
	
	// initial field (from a restart file)
	double *rho_init;
	double (*rhou_init)[3];
	double *rhoE_init;
//	double *kine_init;
//	double *omega_init;
	
	// data on the unstable branch (from Q0_PT000**.bin)
	double *rho_unstable;
	double (*rhou_unstable)[3];
	double *rhoE_unstable;
//	double *kine_unstable;
//	double *omega_unstable;

	// least-stable global modes
	double *rho_1stMode;
	double (*rhou_1stMode)[3];
	double *rhoE_1stMode;
//	double *kine_1stMode;
//	double *omega_1stMode;

	// first adjoint global modes
	double *rho_adj;
	double (*rhou_adj)[3];
	double *rhoE_adj;
//	double *kine_adj;
//	double *omega_adj;

	// past residual
	double residRhoEPast;
	double residRhoEChange;

	// perturbations
	double *array_perturb;

	// filter
    double minFilterWidth;
	double maxFilterWidth;
	double filterBoundaryLength;

//	 // wall distance at the faces (this is required for the filter)
//	double *wdFace;

	// test number 
	int itest;
	int ntests;

	// wall time
	double wallTime;
	
public:
	/*
	 * constructor
	 */
	VanGogh(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name) {
		if(mpi_rank==0)
			cout<<"VanGogh()"<<endl;
		
		//**************
		// HEAT_RELEASE
		//**************
		
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
		hFuel = 1.20e8 * (1/GAMMA)/R_gas * 2.0e-6; //Fuel heating value (120MJ/kg for H2)
		mAir  = 0.02156; //Mass flow rate of injected air

		dComb = getDoubleParam("SHAPE_HEAT_RELEASE",0.75); //Free parameter: Shape of heat release (best result from Vince's experiment - he tried dComb=0.5,0.75,1.0,1.25,1.5)
		kComb = -(1/dComb)*log(1-0.95); //Free parameter: Fraction of completed combustion (DLR simulation result: 95%)

		xCombStart = 0.069; //Combustion ignition starting position

		xMin    = -0.3;
		xMax    = 0.423; //0.775856;
		xInlet  = 0.0;
		xThroat = 0.3;
		yInlet  = 0.0098;

		lengthComb = 0.368; //Combustor length

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
			cout<<"**    Stoichiometic ratio (fSt)                = "<<fSt<<endl;
			cout<<"**    Fuel heating value of H2 (hFuel)         = "<<hFuel<<endl;
			cout<<"**    Mass flow rate of injected air (mAir)    = "<<mAir<<endl;
			cout<<"**    Fraction of completed combustion (kComb) = "<<kComb<<endl;
			cout<<"**    Shape of heat release (dComb)            = "<<dComb<<endl;
			cout<<"**    Ignition starting point (xCombStart)     = "<<xCombStart<<endl;
			cout<<"**    xMin                          = "<<xMin<<endl;
			cout<<"**    xMax                          = "<<xMax<<endl;
			cout<<"**    xInlet                        = "<<xInlet<<endl;
			cout<<"**    xThroat                       = "<<xThroat<<endl;
			cout<<"**    yInlet                        = "<<yInlet<<endl;
			cout<<"**    Combustor length (lengthComb) = "<<lengthComb<<endl<<endl;
		}

		// Heat release
		heatReleaseFunc = NULL;
		registerScalar(heatReleaseFunc,  "HEAT_RELEASE_FUNC",  CV_DATA);

		// Statistics
		mySumHeatRelease_JOE = 0.0;
		
		//**************
		// VAN_GOGH
		//**************
		
		// initial field (from a restart file)
		rho_init = NULL; 	registerScalar(rho_init, "INIT_RHO", CV_DATA);
		rhou_init = NULL; 	registerVector(rhou_init, "INIT_RHOU", CV_DATA);
		rhoE_init = NULL; 	registerScalar(rhoE_init, "INIT_RHOE", CV_DATA);
//		kine_init = NULL; 	registerScalar(kine_init, "INIT_kine", CV_DATA);
//		omega_init = NULL; 	registerScalar(omega_init, "INIT_omega", CV_DATA);

		// data on the unstable branch (from Q0_PT000**.bin)
		rho_unstable = NULL; 	registerScalar(rho_unstable, "UNSTABLE_RHO", CV_DATA);
		rhou_unstable = NULL;	registerVector(rhou_unstable, "UNSTABLE_RHOU", CV_DATA);
		rhoE_unstable = NULL; 	registerScalar(rhoE_unstable, "UNSTABLE_RHOE", CV_DATA);
//		kine_unstable = NULL; 	registerScalar(kine_unstable, "UNSTABLE_kine", CV_DATA);
//		omega_unstable = NULL; 	registerScalar(omega_unstable, "UNSTABLE_omega", CV_DATA);

		// least-stable global modes
		rho_1stMode = NULL; 	registerScalar(rho_1stMode, "FirstMODE_RHO", CV_DATA);
		rhou_1stMode = NULL;	registerVector(rhou_1stMode, "FirstMODE_RHOU", CV_DATA);
		rhoE_1stMode = NULL;	registerScalar(rhoE_1stMode, "FirstMODE_RHOE", CV_DATA);
//		kine_1stMode = NULL;	registerScalar(kine_1stMode, "FirstMODE_kine", CV_DATA);
//		omega_1stMode = NULL;	registerScalar(omega_1stMode, "FirstMODE_omega", CV_DATA);

		// first adjoint global modes
		rho_adj = NULL;		registerScalar(rho_adj, "ADJOINT_RHO", CV_DATA);
		rhou_adj = NULL;	registerVector(rhou_adj, "ADJOINT_RHOU", CV_DATA);
		rhoE_adj = NULL;	registerScalar(rhoE_adj, "ADJOINT_RHOE", CV_DATA);
//		kine_adj = NULL;	registerScalar(kine_adj, "ADJOINT_kine", CV_DATA);
//		omega_adj = NULL;	registerScalar(omega_adj, "ADJOINT_omega", CV_DATA);
		
		// past residual
		residRhoEPast = 0.0;
		residRhoEChange = 0.0;

		// perturbations
		array_perturb = NULL;

		// filter
//		wdFace = NULL;

		// test number
		itest = 0;
		ntests = 0;

		// wall time
		wallTime = 0.0;
	}
	
	/*
	 * destructor
	 */
	virtual ~VanGogh() {
		if(mpi_rank==0)
			cout<<"~VanGogh()"<<endl;
		this->clear();
	}
	
protected:
	/*
	 * Method: clear()
	 * ---------------
	 *
	 */
	void clear() {
		if(lambda != NULL)
			delete [] lambda;	lambda = NULL;
	}
	
	/**********************
	 * OVERLOADED METHODS
	 **********************/
	/*
	 * Method: run
	 * -----------
	 * Original code = run() in JoeWithModels.h
	 */
	virtual void run();
	
	/**************************
	 *---  HOOK FUNCTIONS  ---*
	 **************************/
	/*
	 * Method: initalHook
	 * ------------------
	 *
	 */
	void initialHook() {
		static bool firstCall = true;

		//**************
		// HEAT_RELEASE
		//**************
		double constant1 = fSt*mAir;
		double constant2 = kComb*dComb/lengthComb;

		double dx = (xMax-xMin) / cvora[mpi_size];

		int myCount = 0;
		double myHeatRelease = 0.0;
		for (int icv=0; icv<ncv; ++icv) {
			double xCoord = x_cv[icv][0];

			if (xCoord >= xCombStart) {
				double xN = (xCoord-xCombStart)/lengthComb;
				double func = max( 0.0, pow(kComb*xN,dComb-1)*exp(-pow(kComb*xN,dComb)) );

				heatReleaseFunc[icv] = constant1*constant2*func * dx;

				++myCount;
				myHeatRelease += heatReleaseFunc[icv];
			} else
				heatReleaseFunc[icv] = 0.0;
		}
		updateCvDataG1G2(heatReleaseFunc, REPLACE_DATA);

		int totCount;
		MPI_Allreduce(&myCount, &totCount, 1, MPI_INT, MPI_SUM, mpi_comm);
		double totHeatRelease;
		MPI_Allreduce(&myHeatRelease, &totHeatRelease, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

		if(mpi_rank==0)
			cout<<endl
			    <<"Heat will be released in total "<<totCount<<" CVs: Total heat-release per unit (LAMBDA * hFuel) = "<<totHeatRelease<<endl
			    <<endl;

		//**************
		// VAN_GOGH
		//**************

		if(mpi_rank == 0)
			cout<<"VanGoghKOm::InitialHook() "<<endl;
		assert(checkScalarFlag("RHO")); // The restart file should exist!!

		// -------------------------------------------------
		// update wall time (it will be used in finalHook())
		// -------------------------------------------------
		MPI_Barrier(mpi_comm);
		if (mpi_rank == 0)
			wallTime = MPI_Wtime();

		// -------------------------------------------------
		// random perturbation
		// -------------------------------------------------
		randomPerturbField();
		
		//**************
		// Record 
		//**************
		writeData(0); // Save each case on a different file
		
		char filenameQoISummary[30];
		sprintf(filenameQoISummary, "QoIsummary.csv");
		writeQoIsummaryOnFile(filenameQoISummary, 0, true, firstCall);

		char filenameQoI[30];
		sprintf(filenameQoI, "QoI.case%04d.csv", itest);
		writeQoIOnFile(0, filenameQoI, 0, true); 
		
		char filenameProfile[40];
		sprintf(filenameProfile, "profile.case%04d.step%05d.csv", itest, 0);
		write1Dprofile(filenameProfile);

		char filenameResidual[30];
		sprintf(filenameResidual, "residual.case%04d.csv", itest);
		if(mpi_rank==0) {
			FILE *fp = fopen(filenameResidual, "w");

			fclose(fp);
		}

		firstCall = false;
	}

	//	virtual void boundaryHook(double *temp, double (*vel)[3], double *press, FaZone *zone) {}
	//	virtual void boundaryHook1D(const int ifa, double *T_fa, double (*vel_fa)[3], double *p_fa, FaZone *zone) {}

	/*
	 * Method: sourceHook
	 * ------------------
	 *
	 */
	void sourceHook(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5]) {
		if(lambda != NULL) {
			for (int icv=0; icv<ncv; ++icv) {
				if (x_cv[icv][0] >= xCombStart) {
					rhs_rhoE[icv] += lambda[0]*hFuel*heatReleaseFunc[icv];

					mySumHeatRelease_JOE += lambda[0]*hFuel*heatReleaseFunc[icv];
				}
			}
		}
	}

	/*
	 * Method: temporalHook
	 * --------------------
	 *
	 */
	void temporalHook() {
		// Sum heat release
		double totSumHeatRelease_JOE;
		MPI_Allreduce(&mySumHeatRelease_JOE, &totSumHeatRelease_JOE, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

		if(mpi_rank==0 && step%(check_interval*100)==0) {
			if(fabs(totSumHeatRelease_JOE)>1.0e-10)
				printf("   TOTAL HEAT-RELEASE (JOE)= %.5e\n", totSumHeatRelease_JOE);
		}

		mySumHeatRelease_JOE = 0.0;
				
		//**************
		// Record 
		//**************
		int saveQoI = getIntParam("INTERVAL_SAVE_QOI", "100");
		if(step%saveQoI == 0) {
			char filenameQoI[30];
			sprintf(filenameQoI, "QoI.case%04d.csv", itest);
			writeQoIOnFile(step, filenameQoI, 0, false); 
		}
		
		int saveProfile = getIntParam("INTERVAL_SAVE_PROFILE", "100");
		if(step%saveProfile == 0) {
			char filenameProfile[40];
			sprintf(filenameProfile, "profile.case%04d.step%05d.csv", itest, step);
			write1Dprofile(filenameProfile);
		}
		
		if(step%check_interval==0) {
			char filenameResidual[30];
			sprintf(filenameResidual, "residual.case%04d.csv", itest);
			if(mpi_rank==0) {
				FILE *fp = fopen(filenameResidual, "a");
	
				fprintf(fp,"%d, ", step);
				
				int nScal = scalarTranspEqVector.size();
				for(int i=0; i<5+nScal-1; ++i)
					fprintf(fp, "%11.4e, ", Residual[i]);
				fprintf(fp,"%11.4e\n", Residual[5+nScal-1]);
				
				fclose(fp);
			}
		}
	}

	/*
	 * Method: finalHook
	 * -----------------
	 *
	 */
	void finalHook() {
		char filenameQoISummary[30];
		sprintf(filenameQoISummary, "QoIsummary.csv");
		writeQoIsummaryOnFile(filenameQoISummary, 0, false, false);
	}

	/*
	 * Method: calcPTM
	 * ---------------
	 *
	 */
	void calcPTM(double &p, double &T, double &M, const double gamma, const double Rgas, const double rho, const double rhou, const double rhoE) {
		p = (gamma-1) * (rhoE - 0.5*pow(rhou, 2.0)/rho);
		T = p/rho/Rgas;
		M = rhou/rho / sqrt(fabs(gamma*p/rho));
	}
	
	/*
	 * Method: writeQoIOnFile
	 * ----------------------
	 * Averaged density and Averaged Mach number
	 */
	void writeQoIOnFile(const int step, char filename[], const int NcontrolEqns, bool rewrite) {
		double avgRho;
		double avgMa;
		double subsonicVolPortion;
		double pressMax;
		double pressMin;
		double totVolume;

		double myAvgRho      = 0.0;
		double myAvgMa       = 0.0;
		double mySubsonicVol = 0.0;
		double myPressMax = -1.0e22;
		double myPressMin =  1.0e22;
		double myVolume = 0.0;

		for(int icv=0; icv<ncv; ++icv) {
			if(x_cv[icv][0] >= xInlet && x_cv[icv][0] <= xThroat) {
				myAvgRho += UgpWithCvCompFlow::rho[icv] * cv_volume[icv];

				double localMa = fabs(UgpWithCvCompFlow::rhou[icv][0] / UgpWithCvCompFlow::rho[icv]) / UgpWithCvCompFlow::sos[icv];
				myAvgMa  += localMa*cv_volume[icv];

				myPressMax = max(myPressMax, UgpWithCvCompFlow::press[icv]);
				myPressMin = min(myPressMin, UgpWithCvCompFlow::press[icv]);

				if(localMa < 1.0)
					mySubsonicVol += cv_volume[icv];
				myVolume += cv_volume[icv];
			}
		}

		MPI_Allreduce(&myAvgRho,      &avgRho,             1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(&myAvgMa,       &avgMa,              1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(&myPressMax,    &pressMax,           1, MPI_DOUBLE, MPI_MAX, mpi_comm);
		MPI_Allreduce(&myPressMin,    &pressMin,           1, MPI_DOUBLE, MPI_MIN, mpi_comm);
		MPI_Allreduce(&mySubsonicVol, &subsonicVolPortion, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(&myVolume,      &totVolume,          1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		avgRho /= totVolume;
		avgMa  /= totVolume;
		subsonicVolPortion /= totVolume;

		if(mpi_rank==0) {
			FILE *fp;
			if(rewrite) {
				fp = fopen(filename, "w");
				fprintf(fp, "STEP, AVG_RHO, AVG_MA, SUBSONIC_POR, MAX_PRESS, MIN_PRESS\n", step);
			} else
				fp = fopen(filename, "a");

			fprintf(fp, "%d, ", step);
			fprintf(fp, "%15.8e, %15.8e, %15.8e, %15.8e, %15.8e\n", avgRho, avgMa, subsonicVolPortion, pressMax, pressMin);

			fclose(fp);
		}

		MPI_Barrier(mpi_comm);
	}
	
	/*
	 * Method: writeQoIsummaryOnFile
	 * -----------------------------
	 * Write the initial and final QoI on a file
	 * This method should be called twice: initialHook() and finalHook()
	 */
	void writeQoIsummaryOnFile(char filename[], const int NcontrolEqns, bool initialCalc, bool rewrite) {
		double avgRho;
		double avgMa;
		double subsonicVolPortion;
		double pressMax;
		double pressMin;
		double totVolume;

		double myAvgRho      = 0.0;
		double myAvgMa       = 0.0;
		double mySubsonicVol = 0.0;
		double myPressMax = -1.0e22;
		double myPressMin =  1.0e22;
		double myVolume = 0.0;

		for(int icv=0; icv<ncv; ++icv) {
			if(x_cv[icv][0] >= xInlet && x_cv[icv][0] <= xThroat) {
				myAvgRho += UgpWithCvCompFlow::rho[icv] * cv_volume[icv];

				double localMa = fabs(UgpWithCvCompFlow::rhou[icv][0] / UgpWithCvCompFlow::rho[icv]) / UgpWithCvCompFlow::sos[icv];
				myAvgMa  += localMa*cv_volume[icv];

				myPressMax = max(myPressMax, UgpWithCvCompFlow::press[icv]);
				myPressMin = min(myPressMin, UgpWithCvCompFlow::press[icv]);

				if(localMa < 1.0)
					mySubsonicVol += cv_volume[icv];
				myVolume += cv_volume[icv];
			}
		}

		MPI_Allreduce(&myAvgRho,      &avgRho,             1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(&myAvgMa,       &avgMa,              1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(&myPressMax,    &pressMax,           1, MPI_DOUBLE, MPI_MAX, mpi_comm);
		MPI_Allreduce(&myPressMin,    &pressMin,           1, MPI_DOUBLE, MPI_MIN, mpi_comm);
		MPI_Allreduce(&mySubsonicVol, &subsonicVolPortion, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(&myVolume,      &totVolume,          1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		avgRho /= totVolume;
		avgMa  /= totVolume;
		subsonicVolPortion /= totVolume;

		if(mpi_rank==0) {
			FILE *fp;
			if(rewrite) {
				fp = fopen(filename, "w");
				fprintf(fp, "TEST_NUM, INIT_AVG_RHO, INIT_AVG_MA, INIT_SUBSONIC, INIT_MAX_PRESS, INIT_MIN_PRESS, FIN_AVG_RHO, FIN_AVG_MA, FIN_SUBSONIC, FIN_MAX_PRESS, FIN_MIN_PRESS\n", step);
			} else
				fp = fopen(filename, "a");

			if(initialCalc) {
				fprintf(fp, "%d, ", itest);
				fprintf(fp, "%15.8e, %15.8e, %15.8e, %15.8e, %15.8e,", avgRho, avgMa, subsonicVolPortion, pressMax, pressMin);
			} else {
				fprintf(fp, "%15.8e, %15.8e, %15.8e, %15.8e, %15.8e\n", avgRho, avgMa, subsonicVolPortion, pressMax, pressMin);
			}

			fclose(fp);
		}

		MPI_Barrier(mpi_comm);
	}

	/*
	 * Method: write1Dprofile
	 * ----------------------
	 *
	 */
	void write1Dprofile(const char filename[]) {
		// open the file
		MPI_Status status;
		if(mpi_rank != 0) { MPI_Status status; int dummy; MPI_Recv(&dummy, 1, MPI_INT, mpi_rank-1, 120, mpi_comm, &status); }

		FILE *fp;
		if(mpi_rank==0) {
			fp = fopen(filename, "w");
			fprintf(fp, "X,           PRESSURE,    TEMPERATURE, MACH_NUMBER \n", ncv);
		} else
			fp = fopen(filename, "a");

		// write data
		for(int icv=0; icv<ncv; ++icv) {
			double p, T, M;
			calcPTM(p, T, M, GAMMA, R_gas, rho[icv], rhou[icv][0], rhoE[icv]);
			fprintf(fp, "%.5e, %.5e, %.5e, %.5e\n", x_cv[icv][0], p, T, M);
		}

		// close the file
		fclose(fp);

		if(mpi_rank < mpi_size-1) { int dummy; MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 120, mpi_comm); }
		MPI_Barrier(mpi_comm);
	}
	
	/****************************
	 * UTILITY METHODS
	 ****************************/
	/*
	 * Method: calcRMS
	 * -------------------
	 * Return the root mean square of a field variable (Option: 1. simple  2. volume averaged)
	 */
	double calcRMS(const double *vec, bool volumeAvg);
	
	/*
	 * Method: calcRMS3D
	 * -----------------
	 * Return the root mean square of a field variable (Option: 1. simple  2. volume averaged)
	 */
	double calcRMS3D(const double (*vec)[3], const int dim, bool volumeAvg);
	
	/*
	 * Method: findMinValue
	 * --------------------
	 * Find the minimun value of an array and return it
	 */
	double findMinValue(double* array, const int arrSize);

	/*
	 * Method: findMaxalue
	 * --------------------
	 * Find the maximum value of an array and return it
	 */
	double findMaxValue(double* array, const int arrSize);

	/*
	 * Method: calcWallDistanceFace
	 * ----------------------------
	 * calculate wall distance: minimum distance of cv to wall face by looping over all cv's and faces
	 */
	void calcWallDistanceFace(double *wd);
	
	/****************************
	 * VANGOGH SPECIFIC METHODS
	 ****************************/
	/*
	 * Method: findMatchedIndex
	 * ------------------------
	 *
	 */
	int findMatchedIndex(int& foundIndex, const int icv, const int nCVs, const double (*x_cv_eigen)[3], const double epsil);
	
	/*
	 * Method: readEigenBinaryFileHeader()
	 * ------------------------------
	 * Read the header part of a given file.
	 */
	void readEigenBinaryFileHeader(const string &filename, int &nCVs, int &dimEigen, int &nEigens);
	
	/*
	 * Method: readEigenBinaryFileBody()
	 * ------------------------
	 * Read the body part of a given file.
	 */
	void readEigenBinaryFileBody(const string &filename, double (*x_cv_eigen)[3], complex<double> (*EVecDirect)[NEIGENS], complex<double> (*EVecAdjoint)[NEIGENS]);
	
	/*
	 * Method: readLinSystem
	 * ---------------------
	 * Read eigenmodes and adjoint modes from "dumpedEigen_pt*****.bin"
	 */
	void readLinSystem(const int ptOnSCurve,
			double *rho_1stMode, double (*rhou_1stMode)[3], double *rhoE_1stMode, double *kine_1stMode, double *omega_1stMode,
			double *rho_adj, double (*rhou_adj)[3], double *rhoE_adj, double *kine_adj, double *omega_adj);

	/*
	 * Method: readPsALCdumpedDataSerial
	 * --------------------------
	 * Read field data from file
	 * Original code: readQ0() in joeWithPsALcont.h
	 */
	void readPsALCdumpedDataSerial(const char filename[], double* rho_unstable, double (*rhou_unstable)[3], double* rhoE_unstable);

	
	/*
	 * Method: storeInitRestart
	 * ------------------------
	 * Store initial data in arrays and calculate the norms of them for future uses
	 */
	void storeInitRestart();
	
	/*
	 * Method: reinitializeNS
	 * ----------------------
	 *
	 */
	void reinitializeNS(const double* rho_init, const double (*rhou_init)[3], const double* rhoE_init);
	
	/*
	 * Method: reinitializeScalars
	 * ---------------------------
	 *
	 */
	void reinitializeScalars(const double* kine_init, const double* omega_init);
	
	/*
	 * Method: perturbFieldScalars
	 * ---------------------------
	 * Perturb the initial field with random variables and apply filter to the perturbations
	 */
	void perturbFieldScalars(const RAND_DISTRIB_FUNC randDistribFunc, const bool useSmoothing, const bool useFiltering);

	/*
	 * Method: perturbFieldNS
	 * ----------------------
	 * Perturb the initial field with random variables and apply filter to the perturbations
	 */
	void perturbFieldNS(const RAND_DISTRIB_FUNC randDistribFunc, const bool useSmoothing, const bool useFiltering);

	/*
	 * Method: calcMetricsAdj
	 * ----------------------
	 * Calculate metrics using adjoint: Let phi=flow field, phi0=flow field on the unstable branch, y=least stable eigenvector, a=adjoint vector,
	 *                                           (phi-phi0)*a
	 *                                  metric = ------------
	 *                                                y*a
	 * Return: metricNumer = (phi-phi0)*a
	 *         metricDenom = y*a
	 */
	void calcMetricsAdj(double* metricNumer, double* metricDenom,
			const double *rho_unstable, const double (*rhou_unstable)[3], const double *rhoE_unstable, const double *kine_unstable, const double *omega_unstable,
			const double *rho_1stMode, const double (*rhou_1stMode)[3], const double *rhoE_1stMode, const double *kine_1stMode, const double *omega_1stMode,
			const double *rho_adj, const double (*rhou_adj)[3], const double *rhoE_adj, const double *kine_adj, const double *omega_adj);

	/****************************
	 * FUNCTIONS FOR HOOK FUNCS
	 ****************************/
	/*
	 * Method: randomPerturbField
	 * --------------------------
	 * Reinitialize the flow field with the restart file, perturb the field, and finally filter it.
	 */
	void randomPerturbField() {
		// reinitialize the flow field
		string reinitField = getStringParam("REINIT_FIELD", "YES");
		if(reinitField=="YES" || reinitField=="yes"  || reinitField=="Y" || reinitField=="y") {
			reinitializeNS(rho_init, rhou_init, rhoE_init);
//			reinitializeScalars(kine_init, omega_init);
		}

		// Perturb the field with filtering
		string randomDistribFunc = getStringParam("RANDOM_DISTRIB_FUNC", "NORMAL");
		
		bool useSmoothing = false;
		string smoothing = getStringParam("USE_SMOOTHING", "NO");
		if(smoothing=="YES" || smoothing=="yes" || smoothing=="Y" || smoothing=="y")
			useSmoothing = true;

//		string filtering = getStringParam("APLLY_FILTER", "NO");
		bool useFiltering = false;
//		if(filtering=="YES" || filtering=="yes" || filtering=="Y" || filtering=="y") {
//			useFiltering = true;
//
//			// calculate wall distance at faces
//			if(mpi_rank==0)
//				cout<<"allocate memory for wdFace"<<endl;
//			wdFace = new double [nfa_b2gg];
//			calcWallDistanceFace(wdFace);
//
//			// get filter size
//            minFilterWidth = getDoubleParam("MIN_FILTER_WIDTH", 0.01);
//			maxFilterWidth = getDoubleParam("MAX_FILTER_WIDTH", 1.0);
//			filterBoundaryLength = getDoubleParam("FILTER_BOUNDARY_LENGTH", 0.0);
//		}

		RAND_DISTRIB_FUNC randDistribFunc;
		if(randomDistribFunc == "UNIFORM")
			randDistribFunc = UNIFORM_DISTRIB;
		else
			randDistribFunc = NORMAL_DISTRIB;

//		perturbFieldScalars(randDistribFunc, useSmoothing, useFiltering);
		perturbFieldNS(randDistribFunc, useSmoothing, useFiltering);

		// calculate state variables: press, etc.
		calcStateVariables();
	}
};

#endif
