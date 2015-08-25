#include "vanGoghWithModels.h"
#include "TurbModel_KOM.h"

#ifndef PI
#define PI 3.14159265359
#endif

#define NVAL 8

/*
 * Struct: HyShotGeometry
 * HyShot mesh geometry
 */
struct HyShotGeometry {
	HyShotGeometry() {
		initialized = false;
	}

	bool initialized;

	double xMin;
	double xMax;

	double xInlet;
	double xThroat;

	double combustorLength;

	double yLowerwall;
	double yUpperwall;

	double width;

	double slopeNozzle;

	void showDataOnScreen() {
		cout << "**** HYSHOT GEOMETRY" << endl;
		cout << "**    X: xMin    = " << xMin << endl;
		cout << "**       xMax    = " << xMax << endl;
		cout << "**       xInlet  = " << xInlet << endl;
		cout << "**       xThroat = " << xThroat << endl;
		cout << "**    combustorLength = " << combustorLength << endl;
		cout << "**    Y: yLowerwall = "<<yLowerwall << endl;
		cout << "**       yUpperwall = "<<yUpperwall << endl;
		cout << "**    Width of combustor (width) = " << width << endl;
		cout << "**    Slope of downstream nozzle (slopeNozzle) = " << slopeNozzle << endl;
		cout << "**" << endl;
	}
};

/*
 * Struct: HyShotHeatReleaseModel
 * 1D heat release model
 */
struct HyShotHeatReleaseModel {
	HyShotHeatReleaseModel() {
		initialized = false;
	}

	bool initialized;

	double fSt; //Stoichiometric fuel/air ratio (0.28 for H2/air)
	double hFuel; //Fuel heating value (H2)
	double mAir; //Mass flow rate of injected air

	double lengthComb; //Combustor length
//	double kComb; //Free parameter: Fraction of completed combustion (BLR exp.: 95%)
//	double dComb; //Free parameter: Shape of heat release

	double xCombStart; //Combustion ignition starting position
	double xCombEnd;

	void showDataOnScreen() {
		cout<<"**** SCRAMJET HEAT RELEASE MODEL: FROM C.Doolan & R.Boyce, AIAA 2008-2603"<<endl;
		cout<<"**    Stoichiometic ratio (fSt)                = "<<fSt<<endl;
		cout<<"**    Fuel heating value of H2 (hFuel)         = "<<hFuel<<endl;
		cout<<"**    Mass flow rate of injected air (mAir)    = "<<mAir<<endl;
//		cout<<"**    Fraction of completed combustion (kComb) = "<<kComb<<endl;
//		cout<<"**    Shape of heat release (dComb)            = "<<dComb<<endl;
		cout<<"**    Ignition starting point (xCombStart) = "<<xCombStart<<endl;
		cout<<"**    Combustion end point (xCombEnd)      = "<<xCombEnd<<endl;
		cout<<"**    Combustion zone length (lengthComb)  = "<<lengthComb<<endl;
	}
};

// ###########################################################################################
// ------                                                                               ------
// ------                       Perturbation tests called VanGogh                       ------
// ------                                                                               ------
// ###########################################################################################
class VanGoghWithKOm : public VanGoghWithModels, public RansTurbKOm {
protected:
	// Class name
	string classID;

	// initial field (from a restart file)
	double *kine_init;
	double *omega_init;

	// data on the unstable branch (from Q0_PT000**.bin)
	double *kine_unstable;
	double *omega_unstable;

	// -------------------------
	// Eigenvecs: Tecplot output
	// -------------------------
	// least-stable direct global modes
	double *kine_Dct_1stReal;
	double *omega_Dct_1stReal;

	// first adjoint global modes
	double *kine_Adj_1stReal;
	double *omega_Adj_1stReal;

	// -------------------------
	// HyShot related parameters
	// -------------------------
	// HyShot Geometry
	HyShotGeometry hyShotGeometry;

	// Heat Release Model: Based on C.Doolan & R.Boyce, AIAA 2008-2603
	HyShotHeatReleaseModel hyShotHeatReleaseModel;
	double *heatReleaseFunc;

	// Inlet profile
	int npos;
	vector<vector<double> > inletProfile;

public:
	/*
	 * constructor
	 */
	VanGoghWithKOm(char *name) : VanGoghWithModels(name), JoeWithModels(name), UgpWithCvCompFlow(name) {
		classID = "VanGoghWithKOm";
		if(mpi_rank==0)
			cout<<classID<<"()"<<endl;

		init();

		// -------------------------
		// HyShot related parameters
		// -------------------------
		// HyShot geometry
		setHyShotGeometry(hyShotGeometry);

		// HyShot heat-release model
		double shapeHeatRelease = getDoubleParam("SHAPE_HEAT_RELEASE", 0.75); //Free parameter: Shape of heat release (best result from Vince's experiment - he tried dComb=0.5,0.75,1.0,1.25,1.5)
		setHyShotHeatReleaseModel(hyShotHeatReleaseModel, hyShotGeometry, shapeHeatRelease);

		showHyShotParamsOnScreen(); // Show both geometry and heat-release model

		heatReleaseFunc = NULL; 	registerScalar(heatReleaseFunc, "HEAT_RELEASE_FUNC", CV_DATA);

		// -------------------
		// Read inlet profiles
		// -------------------
		npos = 0;
		string inletFilename = "inletData.txt";
		readInletProfile(inletFilename, inletProfile);

		MPI_Barrier(mpi_comm);
	}

	/*
	 * destructor
	 */
	virtual ~VanGoghWithKOm() {
		if(mpi_rank==0)
			cout<<"~"<<classID<<"()"<<endl;
		clear();

		// -------------------------
		// HyShot related parameters
		// -------------------------
		npos = 0;
		inletProfile.clear();
	}

protected:
	/*
	 * Method: init
	 * ------------
	 *
	 */
	void init() {
		// initial field (from a restart file)
		kine_init  = NULL; 	registerScalar(kine_init,  "INIT_KINE",  CV_DATA);
		omega_init = NULL; 	registerScalar(omega_init, "INIT_OMEGA", CV_DATA);

		// data on the unstable branch (from Q0_PT000**.bin)
		kine_unstable = NULL; 	registerScalar(kine_unstable,  "UNSTABLE_KINE",  CV_DATA);
		omega_unstable = NULL; 	registerScalar(omega_unstable, "UNSTABLE_OMEGA", CV_DATA);

		// -------------------------
		// Eigenvecs: Tecplot output
		// -------------------------
		// least-stable global modes
		kine_Dct_1stReal  = NULL; 	registerScalar(kine_Dct_1stReal,  "DCT1ST_KINE",  CV_DATA);
		omega_Dct_1stReal = NULL;	registerScalar(omega_Dct_1stReal, "DCT1ST_OMEGA", CV_DATA);

		// first adjoint global modes
		kine_Adj_1stReal  = NULL;	registerScalar(kine_Adj_1stReal,  "ADJ1ST_KINE",  CV_DATA);
		omega_Adj_1stReal = NULL;	registerScalar(omega_Adj_1stReal, "ADJ1ST_OMEGA", CV_DATA);
	}

	/*
	 * Method: clear
	 * -------------
	 *
	 */
	void clear() {}


	// -------------------------
	// HyShot related parameters
	// -------------------------
	void setHyShotGeometry(HyShotGeometry &hyShotGeometry) {
		hyShotGeometry.initialized = true;

		hyShotGeometry.xMin = 0.345717;
		hyShotGeometry.xInlet = 0.350028;
		hyShotGeometry.combustorLength = 0.035;
		hyShotGeometry.xThroat = hyShotGeometry.xInlet + hyShotGeometry.combustorLength;
		hyShotGeometry.xMax = hyShotGeometry.xThroat + 0.123;
		hyShotGeometry.yLowerwall = 0.11952;
		hyShotGeometry.yUpperwall = 0.12932;
		hyShotGeometry.width = 0.01; //Combustor width
		hyShotGeometry.slopeNozzle = tan(12*PI/180);
	}

	void setHyShotHeatReleaseModel(HyShotHeatReleaseModel &heatReleaseParams, HyShotGeometry &GeomParams, const double shapeHeatRelease) {
		assert(GeomParams.initialized == true);
		heatReleaseParams.initialized = true;

		heatReleaseParams.fSt   = 0.028; //Stoichiometric fuel/air ratio (0.28 for H2/air)
//		heatReleaseParams.hFuel = 1.20e8; //Fuel heating value (120MJ/kg for H2)
		heatReleaseParams.hFuel = 1.20e-1; //Fuel heating value (120MJ/kg for H2)
		heatReleaseParams.mAir  = 0.053*7.5; //Mass flow rate of injected air (The viscous 2D hyshot case has actual width of 0.01: Thus, it has 0.053*20)

//		heatReleaseParams.kComb = 0.95; //Free parameter: Fraction of completed combustion (DLR simulation result: 95%)
//		heatReleaseParams.dComb = shapeHeatRelease; //Free parameter: Shape of heat release (best result from Vince's experiment - he tried dComb=0.5,0.75,1.0,1.25,1.5)

		heatReleaseParams.xCombStart = GeomParams.xInlet + 0.002851; //Combustion ignition starting position
		heatReleaseParams.xCombEnd   = GeomParams.xThroat;

		heatReleaseParams.lengthComb = heatReleaseParams.xCombEnd - heatReleaseParams.xCombStart; //Combustor length
	}

	void showHyShotParamsOnScreen() {
		if(mpi_rank==0) {
			cout<<endl;
			hyShotGeometry.showDataOnScreen();
			cout<<endl;

			cout<<endl;
			hyShotHeatReleaseModel.showDataOnScreen();
			cout<<endl;
		}
	}

	void initHeatReleaseFunc() {
		int myCountHR = 0;
		double mySumHeatRelease = 0.0;

		// Note: In SourceHook(), heatSource = lambda[0] * hyShotHeatReleaseModel.hFuel * heatReleaseFunc[icv]
		double constant = hyShotHeatReleaseModel.fSt * hyShotHeatReleaseModel.mAir;
		for (int icv=0; icv<ncv; ++icv) {
			double xCoord = x_cv[icv][0];
			double yCoord = x_cv[icv][1];

			if (xCoord >= hyShotHeatReleaseModel.xCombStart && xCoord <= hyShotHeatReleaseModel.xCombEnd && yCoord >= hyShotGeometry.yLowerwall && yCoord <= hyShotGeometry.yUpperwall) {
				double area = (hyShotGeometry.yUpperwall - hyShotGeometry.yLowerwall) * hyShotGeometry.width;

//				double xN = (xCoord-hyShotHeatReleaseModel.xCombStart) / hyShotHeatReleaseModel.lengthComb;
//
//				double constFunc = hyShotHeatReleaseModel.dComb * hyShotHeatReleaseModel.kComb / hyShotHeatReleaseModel.lengthComb;
//				double func = constFunc * pow(hyShotHeatReleaseModel.kComb*xN, hyShotHeatReleaseModel.dComb-1) * exp(-pow(hyShotHeatReleaseModel.kComb*xN, hyShotHeatReleaseModel.dComb));
				double func = 1.0 / hyShotHeatReleaseModel.lengthComb;

				heatReleaseFunc[icv] = constant/area*func*cv_volume[icv];

				++myCountHR;
				mySumHeatRelease += heatReleaseFunc[icv];
			} else
				heatReleaseFunc[icv] = 0.0;
		}

		int totCountHR;
		double totSumHeatRelease
		MPI_Allreduce(&myCountHR,        &totCountHR,        1, MPI_INT,    MPI_SUM, mpi_comm);
		MPI_Allreduce(&mySumHeatRelease, &totSumHeatRelease, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

		if(mpi_rank==0)
			cout << "Heat will be released in total " << totCountHR << " CVs: Total heat-release per unit (LAMBDA * hFuel) = " << totSumHeatRelease << endl
			     << endl;
	}

	void readInletProfile(string &filename, vector<vector<double> > &inletProfile) {
		// read inlet profile
		// X Y TEMP PRESS VEL-X VEL-Y KINE OMEGA
		// 0 1 2    3     4     5     6    7
		FILE *fp;
		if ((fp=fopen(filename.c_str(), "rt")) == NULL) {
			if (mpi_rank == 0)
				cerr << "could not open " << filename << endl;
			throw(-1);
		}

		fscanf(fp, "n=%d\n", &npos);
		inletProfile.resize(npos, vector<double>(NVAL));

		double dummy;
		for (int i=0; i<npos; i++)
			for (int v=0; v<NVAL; v++) {
				fscanf(fp, "%lf", &dummy);
				inletProfile[i][v] = dummy;
			}

		fclose(fp);

		// Show it on the screen
		if (mpi_rank == 0) {
			printf("n=%d\td=%d\n", npos, NVAL);
			for (int i=0; i<npos; i++) {
				for (int v=0; v<NVAL; v++)
					printf("%.6e\t", inletProfile[i][v]);
				printf("\n");
			}
		}
	}

	double getValuesLinearInterp(double ycoord, int ival) {
		int pos=1;
		while ((inletProfile[pos][1] < ycoord) && (pos<npos-1))     pos++;

		double f = (ycoord-inletProfile[pos-1][1]) / (inletProfile[pos][1]-inletProfile[pos-1][1]);
		f = max(0.0, min(1.0, f));
		return ( inletProfile[pos-1][ival] + f*(inletProfile[pos][ival]-inletProfile[pos-1][ival]) ) ;
	}

public:
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

		if(firstCall) {
			if(mpi_rank==0)
				cout<<endl
				<<classID<<"::initialHook()"<<endl;

			// Heat release
			// Check if hyShotGeometry & hyShotHeatReleaseModel have data
			assert(hyShotGeometry.initialized         == true);
			assert(hyShotHeatReleaseModel.initialized == true);

			// Initialization of the shape function of the heat release
			bool isHRFregistered = false;
			for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
				if (strcmp("HEAT_RELEASE_FUNC", data->getName())==0) {
					isHRFregistered = true;
					break;
				}
			}
			if(!isHRFregistered) {
				if(mpi_rank==0) cout<<"ERROR! HEAT_RELEASE_FUNC has not been registered!"<<endl;
				throw(VANGOGH_ERROR_CODE);
			}

			initHeatReleaseFunc();
		}

		firstCall = false;
	}


	/*
	 * Method: boundaryHook
	 * --------------------
	 *
	 */
	void boundaryHook(double *temp, double (*vel)[3], double *press, FaZone *zone) {
		string zoneName = zone->getNameString();
		if(zoneName.compare("inlet")==0) {
			for (int index = 0; index < zone->faVec.size(); ++index) {
				int ifa = zone->faVec[index];
				int icv1 = cvofa[ifa][1];

				temp[icv1]   = getValuesLinearInterp(x_fa[ifa][1], 2);
				vel[icv1][0] = getValuesLinearInterp(x_fa[ifa][1], 4);
				vel[icv1][1] = getValuesLinearInterp(x_fa[ifa][1], 5);
				vel[icv1][2] = 0.0;
				press[icv1]  = getValuesLinearInterp(x_fa[ifa][1], 3);
			}
		}
	}
	/*
	 * Method: boundaryHook1D
	 * ----------------------
	 *
	 */
	void boundaryHook1D(const int ifa, double *T_fa, double (*vel_fa)[3], double *p_fa, FaZone *zone) {
		string zoneName = zone->getNameString();
		if(zoneName.compare("inlet")==0) {
			int icv1 = cvofa[ifa][1];

			T_fa[icv1]      = getValuesLinearInterp(x_fa[ifa][1], 2);
			vel_fa[icv1][0] = getValuesLinearInterp(x_fa[ifa][1], 4);
			vel_fa[icv1][1] = getValuesLinearInterp(x_fa[ifa][1], 5);
			vel_fa[icv1][2] = 0.0;
			p_fa[icv1]      = getValuesLinearInterp(x_fa[ifa][1], 3);
		}
	}

	/*
	 * Method: boundaryHookScalarRansTurb
	 * ----------------------------------
	 *
	 */
	void boundaryHookScalarRansTurb(double *phi, FaZone *zone, const string &scalName) {
		double beta0 = 0.0708;
		double *wallDist = getScalar("wallDist");
		if (zone->getKind() == FA_ZONE_BOUNDARY) {
			Param *param;
			string zoneName = zone->getNameString();

			if(zoneName.compare("inlet")==0) {
				if(scalName.compare("kine") == 0) {
					for (int index = 0; index < zone->faVec.size(); ++index) {
						int ifa = zone->faVec[index];
						int icv1 = cvofa[ifa][1];

						phi[icv1] = getValuesLinearInterp(x_fa[ifa][1], 6);
					}
				} else if(scalName.compare("omega") == 0) {
					for (int index = 0; index < zone->faVec.size(); ++index) {
						int ifa = zone->faVec[index];
						int icv1 = cvofa[ifa][1];

						phi[icv1] = getValuesLinearInterp(x_fa[ifa][1], 7);
					}
				}
			} else if (getParam(param, zoneName)) {
				if ((param->getString()=="WALL" || param->getString()=="HOOK") && (scalName == "omega")) {
					for (int index = 0; index < zone->faVec.size(); ++index) {
						int ifa = zone->faVec[index];
						int icv0 = cvofa[ifa][0];
						int icv1 = cvofa[ifa][1];

//						double muLamCV = calcMuLam(UgpWithCvCompFlow::temp[icv0]);
						double muLamCV = calcMuLam(icv0); // calcMuLam(temp) generates error with FPVA combustion
						phi[icv1] = 6.0*muLamCV/(rho[icv0]*beta0*wallDist[icv0]*wallDist[icv0]);
					}
				}
			}
		}
	}
	/*
	 * Method: boundaryHookScalarRansTurb1D
	 * ------------------------------------
	 *
	 */
	void boundaryHookScalarRansTurb1D(const int ifa, double *phi, FaZone *zone, const string &scalName) {
		double beta0 = 0.0708;
		double *wallDist = getScalar("wallDist");
		if (zone->getKind() == FA_ZONE_BOUNDARY) {
			Param *param;
			string zoneName = zone->getNameString();

			if(zoneName.compare("inlet")==0) {
				if(scalName.compare("kine") == 0) {
					int icv1 = cvofa[ifa][1];

					phi[icv1] = getValuesLinearInterp(x_fa[ifa][1], 6);
				} else if(scalName.compare("omega") == 0) {
					int icv1 = cvofa[ifa][1];

					phi[icv1] = getValuesLinearInterp(x_fa[ifa][1], 7);
				}
			} else if (getParam(param, zone->getName())) {
				if ((param->getString()=="WALL" || param->getString()=="HOOK") && (scalName == "omega")) {
					int icv0 = cvofa[ifa][0];
					int icv1 = cvofa[ifa][1];

					//					double muLamCV = calcMuLam(UgpWithCvCompFlow::temp[icv0]);
					double muLamCV = calcMuLam(icv0); // calcMuLam(temp) generates error with FPVA combustion
					phi[icv1] = 6.0*muLamCV/(rho[icv0]*beta0*wallDist[icv0]*wallDist[icv0]);
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
		static bool firstCall = true;

		if(lambda==NULL) {
			if(mpi_rank == 0)
				cerr<<"ERROR! lambda was not registered!"<<endl;
			throw(VANGOGH_ERROR_CODE);
		} else {
			if(firstCall && mpi_rank == 0)
				cout<<"> sourceHook(): lambda[0] = "<<lambda[0]<<endl;

			for (int icv=0; icv<ncv; ++icv) {
				if (x_cv[icv][0] >= hyShotHeatReleaseModel.xCombStart && x_cv[icv][0] <= hyShotHeatReleaseModel.xCombEnd) {
					rhs_rhoE[icv] += lambda[0] * hyShotHeatReleaseModel.hFuel * heatReleaseFunc[icv];
				}
			}
		}

		firstCall = false;
	}

	/*
	 * Method: temporalHook
	 * --------------------
	 *
	 */
	void temporalHook() {
//		// Sum heat release
//		double totSumHeatRelease_JOE;
//		MPI_Allreduce(&mySumHeatRelease_JOE, &totSumHeatRelease_JOE, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
//
//		if(mpi_rank==0 && step%(check_interval*100)==0) {
//			if(fabs(totSumHeatRelease_JOE)>1.0e-10)
//				printf("   TOTAL HEAT-RELEASE (JOE)= %.5e\n", totSumHeatRelease_JOE);
//		}
//
//		mySumHeatRelease_JOE = 0.0;
//
//		//**************
//		// Record
//		//**************
//		int saveQoI = getIntParam("INTERVAL_SAVE_QOI", "100");
//		if(step%saveQoI == 0) {
//			char filenameQoI[30];
//			sprintf(filenameQoI, "QoI.case%04d.csv", itest);
//			writeQoIOnFile(step, filenameQoI, 0, false);
//		}
//
//		int saveProfile = getIntParam("INTERVAL_SAVE_PROFILE", "100");
//		if(step%saveProfile == 0) {
//			char filenameProfile[40];
//			sprintf(filenameProfile, "profile.case%04d.step%05d.csv", itest, step);
//			write1Dprofile(filenameProfile);
//		}
//
//		if(step%check_interval==0) {
//			char filenameResidual[30];
//			sprintf(filenameResidual, "residual.case%04d.csv", itest);
//			if(mpi_rank==0) {
//				FILE *fp = fopen(filenameResidual, "a");
//
//				fprintf(fp,"%d, ", step);
//
//				int nScal = scalarTranspEqVector.size();
//				for(int i=0; i<5+nScal-1; ++i)
//					fprintf(fp, "%11.4e, ", Residual[i]);
//				fprintf(fp,"%11.4e\n", Residual[5+nScal-1]);
//
//				fclose(fp);
//			}
//		}
	}

	///*
	// * Method: perturbFieldScalars
	// * ---------------------------
	// * Perturb the initial field with random variables and apply filter to the perturbations
	// */
	//void VanGoghKOm::perturbFieldScalars(const RAND_DISTRIB_FUNC randDistribFunc, const bool useSmoothing, const bool useFiltering) {
	//	if(mpi_rank==0)
	//		cout<<"VanGoghKOm::perturbFieldScalars()"<<endl;
	//
	//	double *kine = getR1("kine");
	//	double *omega = getR1("omega");
	//
	//	// Random variable array and filter array
	//	assert(array_perturb == NULL);
	//	array_perturb = new double [ncv_ggff];
	//
	//	double *tempFiltered = NULL;
	//	if(useFiltering)
	//		tempFiltered = new double [ncv_ggff];
	//
	//	// Some parameters for the disturbances
	//	double disturbMag = getDoubleParam("DISTURB_MAGNITUDE", "0.05");
	//	double disturbClip = getDoubleParam("DISTURB_CLIP", "0.001");
	//
	//	// Some parameters for the hat-shaped smoothing function: f(x) = 1.0 - 0.5*(exp(-a*(x-xMin)) + exp(a*(x-xMax))), where a = ln(100)/edgeSize
	//	// Note: At xMin & xMax, f(x) = 0.0. At the "edges", f(x) = 0.99. At the center of the domain, f(x) ~= 1.0
	//	//       The smoothing occurs only if useSmoothing is TRUE.
	//	double disturbSmoothXmin = getDoubleParam("DISTURB_SMOOTH_XMIN", "-1.0e10");
	//	double disturbSmoothXmax = getDoubleParam("DISTURB_SMOOTH_XMAX", "1.0e10");
	//	double disturbSmoothEdgeSize = getDoubleParam("DISTURB_SMOOTH_EDGE_SIZE", "1.0e-5");
	//	double a = log(100.0)/disturbSmoothEdgeSize;
	//	assert(disturbSmoothXmin<disturbSmoothXmax);
	//	assert(a>0);
	//
	//	// -------------------
	//	// Generate kine_perturb
	//	// -------------------
	//	if(randDistribFunc == UNIFORM_DISTRIB) {
	//		randArrayUniform(array_perturb, ncv); // note: randArrayUniform() generates random variables in [0,1]
	//		for(int icv=0; icv<ncv; ++icv)
	//			array_perturb[icv] -= 0.5;        // shifting
	//	} else
	//		randArrayNormal(array_perturb, ncv);  // default = normal distribution
	//
	//	updateCvDataG1G2(array_perturb, REPLACE_DATA);
	//
	//	// -------------------
	//	// Filter kine_perturb
	//	// -------------------
	//	// note: parameters for the filtering such as maxFilterWidth and filterBoundaryLength were already set
	//	if(useFiltering) {
	//		// Check if filter has been constructed.
	//		if(!filterConsturcted())
	//			buildCvDifferentialFilter();
	//		assert(maxFilterWidth > 0.0);
	//
	//		applyDiffFilterG1G2(tempFiltered, array_perturb);
	//		for(int icv=0; icv<ncv; ++icv)
	//			array_perturb[icv] = tempFiltered[icv];
	//		updateCvDataG1G2(array_perturb, REPLACE_DATA);
	//	}
	//
	//	// -------------------
	//	// Update kine
	//	// -------------------
	//	double rmsVal = calcRMS(kine, false); // get RMS
	//	double coeff = rmsVal*disturbMag; // calculate the coefficient based on the RMS
	//	int myCountClip = 0;
	//	for(int icv=0; icv<ncv; ++icv) {
	//		double disturb = coeff*array_perturb[icv];
	//		if(useSmoothing) {
	//			double x = x_cv[icv][0];
	//			if(x<disturbSmoothXmin || x>disturbSmoothXmax)
	//				disturb = 0.0;
	//			else
	//				disturb *= 1.0 - 0.5*(exp(-a*(x-disturbSmoothXmin)) + exp(a*(x-disturbSmoothXmax)));
	//		}
	//		if(kine[icv]+disturb < disturbClip*kine[icv]) {
	//			kine[icv] = disturbClip*kine[icv];
	//			++myCountClip;
	//		} else
	//			kine[icv] += disturb;
	//	}
	//	updateCvDataG1G2(kine,  REPLACE_DATA);
	//	int countClip; 	MPI_Allreduce(&myCountClip, &countClip, 1, MPI_INT, MPI_SUM, mpi_comm);
	//	if(mpi_rank == 0) {
	//		if(randDistribFunc == UNIFORM_DISTRIB)
	//			printf("    Perturb the kine field with Uniform:  min disturb = %.3e, max disturb = %.3e, # clip = %d \n", -0.5*coeff, 0.5*coeff, countClip);
	//		else
	//			printf("    Perturb the kine field with Gaussian:  std = %.3e, # clip = %d \n", coeff, countClip);
	//	}
	//
	//	// -------------------
	//	// Generate omega_perturb
	//	// -------------------
	//	if(randDistribFunc == UNIFORM_DISTRIB) {
	//		randArrayUniform(array_perturb, ncv); // note: randArrayUniform() generates random variables in [0,1]
	//		for(int icv=0; icv<ncv; ++icv)
	//			array_perturb[icv] -= 0.5;        // shifting
	//	} else
	//		randArrayNormal(array_perturb, ncv);  // default = normal distribution
	//	updateCvDataG1G2(array_perturb, REPLACE_DATA);
	//
	//	// -------------------
	//	// Filter omega_perturb
	//	// -------------------
	//	if(useFiltering) {
	//		applyDiffFilterG1G2(tempFiltered, array_perturb);
	//		for(int icv=0; icv<ncv; ++icv)
	//			array_perturb[icv] = tempFiltered[icv];
	//		updateCvDataG1G2(array_perturb, REPLACE_DATA);
	//	}
	//
	//	// -------------------
	//	// Update omega
	//	// -------------------
	//	rmsVal = calcRMS(omega, false); // get RMS
	//	coeff = rmsVal*disturbMag; // calculate the coefficient based on the RMS
	//
	//	myCountClip = 0;
	//	for(int icv=0; icv<ncv; ++icv) {
	//		double disturb = coeff*array_perturb[icv];
	//		if(useSmoothing) {
	//			double x = x_cv[icv][0];
	//			if(x<disturbSmoothXmin || x>disturbSmoothXmax)
	//				disturb = 0.0;
	//			else
	//				disturb *= 1.0 - 0.5*(exp(-a*(x-disturbSmoothXmin)) + exp(a*(x-disturbSmoothXmax)));
	//		}
	//		if(omega[icv]+disturb < disturbClip*omega[icv]) {
	//			omega[icv] = disturbClip*omega[icv];
	//			myCountClip++;
	//		} else
	//			omega[icv] += disturb;
	//	}
	//	updateCvDataG1G2(omega,  REPLACE_DATA);
	//	MPI_Allreduce(&myCountClip, &countClip, 1, MPI_INT, MPI_SUM, mpi_comm);
	//	if(mpi_rank == 0) {
	//		if(randDistribFunc == UNIFORM_DISTRIB)
	//			printf("    Perturb the omega field with Uniform: min disturb = %.3e, max disturb = %.3e, # clip = %d \n", -0.5*coeff, 0.5*coeff, countClip);
	//		else
	//			printf("    Perturb the omega field with Gaussian: std = %.3e, # clip = %d \n", coeff, countClip);
	//	}
	//
	//	// -------------
	//	// Clear memory
	//	// -------------
	//	delete [] array_perturb; array_perturb = NULL;
	//	if(useFiltering)
	//		delete [] tempFiltered;
	//
	//	if(mpi_rank==0)
	//		cout<<endl;
	//	MPI_Barrier(mpi_comm);
	//}


	/*
	 * Method: calcGradRhoIfNecessary
	 * -------------------------------
	 * Calculate gradRho if it is not done by the solver (first order in space)
	 * Check if gradRho is a NULL pointer. If NULL, assign memory for it
	 */
	void calcGradRhoIfNecessary() {
		static bool firstCall = true;
		static bool calcGradRho = false;

		if(firstCall) {
			if(sndOrder == false)
				assert(UgpWithCvCompFlow::grad_rho == NULL); // In this case, grad_rho should not have been assigned.

			if(UgpWithCvCompFlow::grad_rho == NULL) {
				if(mpi_rank == 0)
					cout<<classID<<"::calcGradRhoIfNecessary(): grad_rho has not been assigned. Allocate grad_rho."<<endl;

				UgpWithCvCompFlow::grad_rho = new double[ncv_g][3]; // This will be cleared in the destructor of UgpCvCompFlow
				calcGradRho = true;
			}
		}

		if(calcGradRho)
			calcCv2Grad(UgpWithCvCompFlow::grad_rho, rho, limiterNavierS, rho, epsilonSDWLS);

		firstCall = false;
	}

	/*
	 * Method: showQoI
	 * ---------------
	 *
	 */
	void showQoI() {
		string funcID = classID + "::showQoI()";

		if(UgpWithCvCompFlow::diverg==NULL ) {
			cerr<<funcID<<": diverg has not been assigned"<<endl;
			throw(-1);
		}

		calcGradRhoIfNecessary();

		//=================
		// Calculate QoI's
		//=================
		// 1 - Shock impingement point on the upperwall
		string faZoneName = "upperwall";

		double shockLocUpperwall;
		double shockStrengthUpperwall;

		// Find local shock location mainly based on the magnitude of grad_rho
		double myShockLocUpperwall      =  ABSURDLY_BIG_NUMBER;
		double myShockStrengthUpperwall = -ABSURDLY_BIG_NUMBER;
		for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
			if (zone->getKind() == FA_ZONE_BOUNDARY) {
				if (faZoneName.compare(zone->getName()) == 0) {
					for (int index = 0; index < zone->faVec.size(); ++index) {
						int ifa = zone->faVec[index];
						int icv0 = cvofa[ifa][0];

						double grad_rho_mag = sqrt(vecDotVec3d(UgpWithCvCompFlow::grad_rho[icv0], UgpWithCvCompFlow::grad_rho[icv0]));
						if(grad_rho_mag >= 100.0 && UgpWithCvCompFlow::diverg[icv0] <= 0.0) {   // This value of 100 works for inviscid Hyshot
							if(myShockStrengthUpperwall < grad_rho_mag) {
								myShockStrengthUpperwall = grad_rho_mag;
								myShockLocUpperwall      = x_fa[ifa][0];
							}
						}
					}
				}
			}
		}
		// Find global shock location
		MPI_Allreduce(&myShockStrengthUpperwall, &shockStrengthUpperwall, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
		if(myShockStrengthUpperwall < shockStrengthUpperwall)
			myShockLocUpperwall = ABSURDLY_BIG_NUMBER;
		MPI_Allreduce(&myShockLocUpperwall, &shockLocUpperwall, 1, MPI_DOUBLE, MPI_MIN, mpi_comm);

		// 2 - Volume-averaged density, Averaged density gradient
		double avgRho,  avgRhoGradMag;
		double totVolume;

		double myAvgRho = 0.0, myAvgRhoGradMag = 0.0;
		double myVolume = 0.0;

		for(int icv=0; icv<ncv; ++icv) {
			if(x_cv[icv][0] >= hyShotGeometry.xMin && x_cv[icv][0] <= hyShotGeometry.xMax) {
				myAvgRho += UgpWithCvCompFlow::rho[icv] * cv_volume[icv];

				if(UgpWithCvCompFlow::grad_rho != NULL)
					myAvgRhoGradMag += sqrt(vecDotVec3d(UgpWithCvCompFlow::grad_rho[icv], UgpWithCvCompFlow::grad_rho[icv])) * cv_volume[icv];

				myVolume += cv_volume[icv];
			}
		}

		MPI_Allreduce(&myAvgRho,        &avgRho,        1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(&myAvgRhoGradMag, &avgRhoGradMag, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

		MPI_Allreduce(&myVolume, &totVolume, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

		avgRho        /= totVolume;
		avgRhoGradMag /= totVolume;

		// 3 - Wall-pressures at the throat
		double myUpperwallPress = 0.0;
		double myLowerwallPress = 0.0;

		double dxNearThroat = 0.75 * 0.0002;
		for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
			if (zone->getKind() == FA_ZONE_BOUNDARY) {
				if (faZoneName.compare("upperwall") == 0) {
					for (int index = 0; index < zone->faVec.size(); ++index) {
						int ifa = zone->faVec[index];
						int icv0 = cvofa[ifa][0];

						if(x_fa[ifa][0] < hyShotGeometry.xThroat && fabs(x_fa[ifa][0] - hyShotGeometry.xThroat) < dxNearThroat)
							myUpperwallPress = UgpWithCvCompFlow::press[icv0];
					}
				}
				if (faZoneName.compare("lowerwall") == 0) {
					for (int index = 0; index < zone->faVec.size(); ++index) {
						int ifa = zone->faVec[index];
						int icv0 = cvofa[ifa][0];

						if(x_fa[ifa][0] < hyShotGeometry.xThroat && fabs(x_fa[ifa][0] - hyShotGeometry.xThroat) < dxNearThroat)
							myLowerwallPress = UgpWithCvCompFlow::press[icv0];
					}
				}
			}
		}

		double upperwallPress;
		double lowerwallPress;
		MPI_Allreduce(&myUpperwallPress, &upperwallPress, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
		MPI_Allreduce(&myLowerwallPress, &lowerwallPress, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);

		//=================
		// Show QoI's on the screen
		//=================
		if(mpi_rank==0) {
			printf("\n");
			printf("> QoI: \n");
			printf(" STEP,            LAMBDA,  SHOCK_LOC_UP,       AVG_RHO,  AVG_RHO_GRAD, THR_PRESS_UPWALL, THR_PRESS_DWNWALL\n");
			printf("%5d, %17.10e, %13.6e, %13.6e, %13.6e, %13.6e, %13.6e\n", step, lambda[0], shockLocUpperwall, avgRho, avgRhoGradMag, upperwallPress, lowerwallPress);
			printf("\n");
		}

		MPI_Barrier(mpi_comm);
	}


	/*
	 * Method: writeQoIOnFile
	 * ----------------------
	 * Averaged density and Averaged Mach number
	 */
	void writeQoIOnFile(const int itest, char filename[], bool rewrite) {
//		double avgRho;
//		double avgMa;
//		double subsonicVolPortion;
//		double pressMax;
//		double pressMin;
//		double totVolume;
//
//		double myAvgRho      = 0.0;
//		double myAvgMa       = 0.0;
//		double mySubsonicVol = 0.0;
//		double myPressMax = -1.0e22;
//		double myPressMin =  1.0e22;
//		double myVolume = 0.0;
//
//		for(int icv=0; icv<ncv; ++icv) {
//			if(x_cv[icv][0] >= xInlet && x_cv[icv][0] <= xThroat) {
//				myAvgRho += UgpWithCvCompFlow::rho[icv] * cv_volume[icv];
//
//				double localMa = fabs(UgpWithCvCompFlow::rhou[icv][0] / UgpWithCvCompFlow::rho[icv]) / UgpWithCvCompFlow::sos[icv];
//				myAvgMa  += localMa*cv_volume[icv];
//
//				myPressMax = max(myPressMax, UgpWithCvCompFlow::press[icv]);
//				myPressMin = min(myPressMin, UgpWithCvCompFlow::press[icv]);
//
//				if(localMa < 1.0)
//					mySubsonicVol += cv_volume[icv];
//				myVolume += cv_volume[icv];
//			}
//		}
//
//		MPI_Allreduce(&myAvgRho,      &avgRho,             1, MPI_DOUBLE, MPI_SUM, mpi_comm);
//		MPI_Allreduce(&myAvgMa,       &avgMa,              1, MPI_DOUBLE, MPI_SUM, mpi_comm);
//		MPI_Allreduce(&myPressMax,    &pressMax,           1, MPI_DOUBLE, MPI_MAX, mpi_comm);
//		MPI_Allreduce(&myPressMin,    &pressMin,           1, MPI_DOUBLE, MPI_MIN, mpi_comm);
//		MPI_Allreduce(&mySubsonicVol, &subsonicVolPortion, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
//		MPI_Allreduce(&myVolume,      &totVolume,          1, MPI_DOUBLE, MPI_SUM, mpi_comm);
//		avgRho /= totVolume;
//		avgMa  /= totVolume;
//		subsonicVolPortion /= totVolume;
//
//		if(mpi_rank==0) {
//			FILE *fp;
//			if(rewrite) {
//				fp = fopen(filename, "w");
//				fprintf(fp, "STEP, AVG_RHO, AVG_MA, SUBSONIC_POR, MAX_PRESS, MIN_PRESS\n");
//			} else
//				fp = fopen(filename, "a");
//
//			fprintf(fp, "%d, ", step);
//			fprintf(fp, "%15.8e, %15.8e, %15.8e, %15.8e, %15.8e\n", avgRho, avgMa, subsonicVolPortion, pressMax, pressMin);
//
//			fclose(fp);
//		}
//
//		MPI_Barrier(mpi_comm);
	}

	/*
	 * Method: writeResidualOnFile
	 * ---------------------------
	 * Averaged density and Averaged Mach number
	 */
	void writeResidualOnFile(const int itest, char filename[], bool rewrite) {

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
	sprintf(inputFileName, "VanGogh.in");

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
		JoeWithModels *vanGogh;

		switch (run) {
			case 1:  vanGogh = new VanGoghWithKOm(inputFileName); 	break;
			default:
				if (mpi_rank == 0)
					cerr << "WARNING: run number not available!" << endl;
				vanGogh = new VanGoghWithKOm(inputFileName);
		}

		// run vanGogh
		vanGogh->run();

		MPI_Barrier(mpi_comm);
		if (mpi_rank == 0) {
			double wtime0 = wtime;
			wtime = MPI_Wtime();
			cout << " > total runtime [s]: " << wtime - wtime0 << endl;
		}

		// delete vanGogh (make sure memory is deallocated in destructors
		delete vanGogh;
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

