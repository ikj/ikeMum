#include "rembrandtWithModels.h"

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
	double kComb; //Free parameter: Fraction of completed combustion (BLR exp.: 95%)
	double dComb; //Free parameter: Shape of heat release

	double xCombStart; //Combustion ignition starting position
	double xCombEnd;

	void showDataOnScreen() {
		cout<<"**** SCRAMJET HEAT RELEASE MODEL: FROM C.Doolan & R.Boyce, AIAA 2008-2603"<<endl;
		cout<<"**    Stoichiometic ratio (fSt)                = "<<fSt<<endl;
		cout<<"**    Fuel heating value of H2 (hFuel)         = "<<hFuel<<endl;
		cout<<"**    Mass flow rate of injected air (mAir)    = "<<mAir<<endl;
		cout<<"**    Fraction of completed combustion (kComb) = "<<kComb<<endl;
		cout<<"**    Shape of heat release (dComb)            = "<<dComb<<endl;
		cout<<"**    Ignition starting point (xCombStart) = "<<xCombStart<<endl;
		cout<<"**    Combustion end point (xCombEnd)      = "<<xCombEnd<<endl;
		cout<<"**    Combustion zone length (lengthComb)  = "<<lengthComb<<endl;
	}
};

// ###########################################################################################
// ------                                                                               ------
// ------      RembrandtWithKOm: Eigen-decomposition with the Wilcox k-omega model      ------
// ------                                                                               ------
// ###########################################################################################

class RembrandtWithKOm : public RembrandtWithModels, public IkeRansTurbKOm_AD {
protected:
	// Class name
	string classID;

	// -------------------------
	// Tecplot output
	// -------------------------
	// The least-stable DIRECT modes (Since the least-stable modes are real, we don't store their imaginary parts)
	double* FirstDirectMode_real_kine;
	double* FirstDirectMode_real_omega;

	// The least-stable ADJOINT modes (Since the least-stable modes are real, we don't store their imaginary parts)
	double* FirstAdjointMode_real_kine;
	double* FirstAdjointMode_real_omega;

	// Structural sensitivity
	double* kine_sens;
	double* omega_sens;

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
	RembrandtWithKOm(char *name) : RembrandtWithModels(name), JoeWithModels(name), JoeWithModels_AD(name), UgpWithCvCompFlow(name) {
		classID = "RembrandtWithKOm";
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

		// ----------------------------
		// Allocate memory for "lambda"
		// ----------------------------
		NcontrolEqns = getIntParam("N_PSALC_CONTROL_PARAMS", "0");
		if(NcontrolEqns > 0) {
			assert(lambda == NULL);

			lambda = new double [NcontrolEqns]; // lambda is the input parameter to the system (e.g. fuel equivalence ratio)
			char paramName[20];
			for(int i=0; i<NcontrolEqns; ++i) {
				sprintf(paramName, "CONTROL_PARAM%d", i);
				lambda[i] = getParam("LAMBDA_INITIAL_0") -> getDouble(paramName); // Note: lambda is same as phi in the model
			}

			if(mpi_rank == 0)
				cout<<endl
				    <<"> LAMBDA0 = "<<lambda[0]<<endl
				    <<endl;
		} else {
			lambda = NULL;
		}

		// -------------------
		// Read inlet profiles
		// -------------------
		npos = 0;
		string inletFilename = "inletData.txt";
		readInletProfile(inletFilename, inletProfile);

		MPI_Barrier(mpi_comm);
	}
	virtual ~RembrandtWithKOm() {
		if(mpi_rank==0)
			cout<<"~"<<classID<<"()"<<endl;
		clear();

		// -------------------------
		// HyShot related parameters
		// -------------------------
		npos = 0;
		inletProfile.clear();
	}

	void init() {
		// -------------------------
		// Tecplot output
		// -------------------------
		FirstDirectMode_real_kine  = NULL; 	registerScalar(FirstDirectMode_real_kine, "FirstDirectMode_REAL_KINE", CV_DATA);
		FirstDirectMode_real_omega = NULL; 	registerScalar(FirstDirectMode_real_omega,"FirstDirectMode_REAL_OMEGA", CV_DATA);

		FirstAdjointMode_real_kine  = NULL; 	registerScalar(FirstAdjointMode_real_kine, "FirstAdjointMode_REAL_KINE", CV_DATA);
		FirstAdjointMode_real_omega = NULL; 	registerScalar(FirstAdjointMode_real_omega,"FirstAdjointMode_REAL_OMEGA", CV_DATA);

		kine_sens  = NULL;	registerScalar(kine_sens,  "SENSITIVITY_kine", CV_DATA);
		omega_sens = NULL;	registerScalar(omega_sens, "SENSITIVITY_omega", CV_DATA);
	}

	void clear() {}

	void postEigenDecompScalarRansTurbModel(double** directEvecsReal, double** adjointEvecsReal, const int directLeastStableIndex, const int adjointLeastStableIndex) {
		assert(directLeastStableIndex  >= 0);
		assert(adjointLeastStableIndex >= 0);

		int nVars = 5 + nScal;

		for(int icv=0; icv<ncv; ++icv) {
			int indexTemp = icv*nVars;

			// -------------------------------------------------
			// Tecplot output for the least-stable direct modes
			// -------------------------------------------------
			FirstDirectMode_real_kine[icv]  = directEvecsReal[directLeastStableIndex][indexTemp+5];
			FirstDirectMode_real_omega[icv] = directEvecsReal[directLeastStableIndex][indexTemp+6];

			// -------------------------------------------------
			// Tecplot output for the least-stable adjoint modes
			// -------------------------------------------------
			FirstAdjointMode_real_kine[icv]  = adjointEvecsReal[adjointLeastStableIndex][indexTemp+5];
			FirstAdjointMode_real_omega[icv] = adjointEvecsReal[adjointLeastStableIndex][indexTemp+6];

			// -------------------------------------------------
			// Structural sensitivity
			// -------------------------------------------------
			// Note: COEFF_BOOST_SENS is defined in rembrandtWithModels.h
			kine_sens[icv] 	= fabs(COEFF_BOOST_SENS * FirstDirectMode_real_kine[icv] * FirstAdjointMode_real_kine[icv]);
			omega_sens[icv]	= fabs(COEFF_BOOST_SENS * FirstDirectMode_real_omega[icv] * FirstAdjointMode_real_omega[icv]);
		}

		updateCvDataG1G2(FirstDirectMode_real_kine, REPLACE_DATA);
		updateCvDataG1G2(FirstDirectMode_real_omega,REPLACE_DATA);

		updateCvDataG1G2(FirstAdjointMode_real_kine, REPLACE_DATA);
		updateCvDataG1G2(FirstAdjointMode_real_omega,REPLACE_DATA);

		updateCvDataG1G2(kine_sens, REPLACE_DATA);
		updateCvDataG1G2(omega_sens,REPLACE_DATA);
	}

	// -------------------------
	// HyShot related parameters
	// -------------------------
	void setHyShotGeometry(HyShotGeometry &hyShotGeometry) {
		hyShotGeometry.initialized = true;

		hyShotGeometry.xMin = 0.345717;
		hyShotGeometry.xInlet = 0.350028;
		hyShotGeometry.combustorLength = 0.090; // This is important because "combustorLength" defines "xThroat"
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

		heatReleaseParams.kComb = 0.95; //Free parameter: Fraction of completed combustion (DLR simulation result: 95%)
		heatReleaseParams.dComb = shapeHeatRelease; //Free parameter: Shape of heat release (best result from Vince's experiment - he tried dComb=0.5,0.75,1.0,1.25,1.5)

		heatReleaseParams.xCombStart = GeomParams.xInlet + 0.057972; //Combustion ignition starting position
		heatReleaseParams.xCombEnd   = GeomParams.xThroat;

		heatReleaseParams.lengthComb = 0.365; // The length of combustion for the cae of duct length = 300mm
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

			if (xCoord >= hyShotHeatReleaseModel.xCombStart && yCoord >= hyShotGeometry.yLowerwall) {
				if ((xCoord >= hyShotGeometry.xThroat) || (xCoord < hyShotGeometry.xThroat && yCoord <= hyShotGeometry.yUpperwall)) {
					double area = 1.0;
					if(xCoord <= hyShotGeometry.xThroat)
						area = (hyShotGeometry.yUpperwall - hyShotGeometry.yLowerwall) * hyShotGeometry.width;
					else {
						double height = (hyShotGeometry.yUpperwall - hyShotGeometry.yLowerwall)
									+ hyShotGeometry.slopeNozzle * (xCoord - hyShotGeometry.xThroat);
						area = height * hyShotGeometry.width;
					}

					double xN = (xCoord-hyShotHeatReleaseModel.xCombStart) / hyShotHeatReleaseModel.lengthComb;

					double constFunc = hyShotHeatReleaseModel.dComb * hyShotHeatReleaseModel.kComb / hyShotHeatReleaseModel.lengthComb;
					double func = constFunc * pow(hyShotHeatReleaseModel.kComb*xN, hyShotHeatReleaseModel.dComb-1) * exp(-pow(hyShotHeatReleaseModel.kComb*xN, hyShotHeatReleaseModel.dComb));

					heatReleaseFunc[icv] = constant / area * func;

					++myCountHR;
					mySumHeatRelease += heatReleaseFunc[icv] * cv_volume[icv];
				}
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

	/*
	 * Method: initialHook
	 * -------------------
	 *
	 */
	void initialHook() {
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
			throw(REMBRANDT_ERROR_CODE);
		}

		initHeatReleaseFunc();
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
	 * Method: boundaryHook1D_AD
	 * -------------------------
	 * Original code = boundaryHook_AD() in JoeWithModelsAD.h
	 */
	void boundaryHook1D_AD(const int ifa, ADscalar<REALQ> &T_fa, ADvector<REALQ> &vel_fa, ADscalar<REALQ> &p_fa, FaZone *zone) {
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
	 * Method: boundaryHookScalarRansTurb1D_AD
	 * ---------------------------------------
	 *
	 */
	void boundaryHookScalarRansTurb1D_AD(const int ifa, ADscalar<adouble> &phi, FaZone *zone, const string &scalName) {
		int debugLevel = getDebugLevel();

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
					phi[icv1] = 6.0*calcMuLam_AD(icv0) / (rho_AD[icv0]*beta0*wallDist[icv0]*wallDist[icv0]);
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
			throw(REMBRANDT_ERROR_CODE);
		} else {
			if(firstCall && mpi_rank == 0)
				cout<<"> sourceHook(): lambda[0] = "<<lambda[0]<<endl;

			for (int icv=0; icv<ncv; ++icv) {
				if (x_cv[icv][0] >= hyShotHeatReleaseModel.xCombStart) {
					rhs_rhoE[icv] += lambda[0] * hyShotHeatReleaseModel.hFuel * heatReleaseFunc[icv] * cv_volume[icv];
				}
			}
		}

		firstCall = false;
	}

	/*
	 * Method: sourceHook1D_AD
	 * -----------------------
	 * Original code = sourceHook_AD() in JoeWithModelsAD.h
	 */
	void sourceHook1D_AD(const int icvCenter, REALA &rhs_rho_AD, REALA rhs_rhou_AD[3], REALA &rhs_rhoE_AD, double (*A)[5][5]) {
        static bool firstCall = true;

		if(lambda_AD==NULL) { // If the Jacobian doesn't need to produce the dRHSdLambda terms (i.e. Newton method for the second point)
		                      // (Note that NcontrolEqns is always non-zero; so never use something like "if(NcontrolEqns==0)")
			assert(lambda!=NULL);
			if(mpi_rank == 0 && firstCall)
				cout<<endl
				    <<">> sourceHook1D_AD(): lambda_AD[] is NULL. Use lambda[]: lambda[0] = "<<lambda[0]<<endl
				    <<endl;

			if (x_cv[icvCenter][0] >= hyShotHeatReleaseModel.xCombStart) {
				rhs_rhoE_AD += lambda[0] * hyShotHeatReleaseModel.hFuel * heatReleaseFunc[icvCenter] * cv_volume[icvCenter];;
			}
		} else {
			if(compCheckOn) { // We want to allow lambda_AD only during the compatibility check
				if(mpi_rank == 0 && firstCall)
					cout<<endl
					    <<">> sourceHook1D_AD(): lambda_AD[] is not NULL for compatibility check"<<endl
					    <<endl;

				if (x_cv[icvCenter][0] >= hyShotHeatReleaseModel.xCombStart) {
					rhs_rhoE_AD += lambda_AD[0] * hyShotHeatReleaseModel.hFuel * heatReleaseFunc[icvCenter] * cv_volume[icvCenter];;
				}
			} else {
				if(mpi_rank == 0)
					cerr<<"ERROR in sourceHook1D_AD(): lambda_AD[] is alive"<<endl;
				throw(REMBRANDT_ERROR_CODE);
			}
		}

        firstCall = false;
	}

	/*
	 * Method: finalHook
	 * -----------------
	 *
	 */
	void finalHook() {
		showQoI();
	}

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
	// set run to default value of 0 => instantiation of JoeWithModels
	int run = 0;
	int ptOnSCurve = 0;

	// set input name to default "Joe.in"
	char inputFileName[50];
	sprintf(inputFileName, "Rembrandt.in");

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

		// declare pointer to JoeWithModels
		JoeWithModels *rembrandt;

		switch (run) {
			case 0:
				rembrandt = new RembrandtWithModels(inputFileName);   break;
			case 1:
				rembrandt = new RembrandtWithKOm(inputFileName);      break;
			default:
				if (mpi_rank == 0)
					cerr << "ERROR: run number not available!" << endl;
				throw(-1);
		}
    
		// provide total runtime
		double wtime, wtime0;
		MPI_Barrier(mpi_comm);
		if (mpi_rank == 0)
			wtime = MPI_Wtime();

		// run joe
		rembrandt->run();
    
		MPI_Barrier(mpi_comm);
		if (mpi_rank == 0) {
			double wtime0 = wtime;
			wtime = MPI_Wtime();
			cout << " > total runtime [s]: " << wtime - wtime0 << endl;
		}

		// delete joe (make sure memory is deallocated in destructors
		delete rembrandt;
	}

	catch (int e) {
		cerr << "Exception: " << e << endl;
		MPI_Finalize();
		return(-1);
	}
	catch (const invalid_argument& problem) {
		if(mpi_rank==0)
			cout << problem.what() <<endl;
		MPI_Finalize();
		return(-1);
	}
	catch (const runtime_error& problem) {
		if(mpi_rank==0)
			cout << problem.what() <<endl;
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

