#include "yangWithModels.h"
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
		cout << "**    Y: yLowerwall = " << yLowerwall << endl;
		cout << "**       yUpperwall = " << yUpperwall << endl;
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
// ------                   Natural parameter continuation called yang                  ------
// ------                                                                               ------
// ###########################################################################################
class YangWithKOm : public YangWithModels, public RansTurbKOm {
protected:
	// Class name
	string classID;

	// -------------------------
	// HyShot related parameters
	// -------------------------
	// HyShot Geometry
	HyShotGeometry hyShotGeometry;

	// Heat Release Model: Based on C.Doolan & R.Boyce, AIAA 2008-2603
	HyShotHeatReleaseModel hyShotHeatReleaseModel;
	double *heatReleaseFunc;

	// Statistics
	double sumHeatRelease;

	// Inlet profile
	int npos;
	vector<vector<double> > inletProfile;

public:
	/*
	 * constructor
	 */
	YangWithKOm(char *name) : YangWithModels(name), JoeWithModels(name), UgpWithCvCompFlow(name) {
		classID = "YangWithKOm";
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
	virtual ~YangWithKOm() {
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
	void init() { /* empty */ }

	/*
	 * Method: clear
	 * -------------
	 *
	 */
	void clear() { /* empty */ }

	// -------------------------
	// HyShot related parameters
	// -------------------------
	void setHyShotGeometry(HyShotGeometry &hyShotGeometry) {
		hyShotGeometry.initialized = true;

		hyShotGeometry.xMin = 0.345717;
		hyShotGeometry.xInlet = 0.350028;
		hyShotGeometry.combustorLength = 0.3; // This is important because "combustorLength" defines "xThroat"
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
		heatReleaseParams.hFuel = 1.20e-1;
		heatReleaseParams.mAir  = 0.053*7.5; //Mass flow rate of injected air (The viscous 2D hyshot case has actual width of 0.01: Thus, it has 0.053*20)

		heatReleaseParams.kComb = 0.95; //Free parameter: Fraction of completed combustion (DLR simulation result: 95%)
		heatReleaseParams.dComb = shapeHeatRelease; //Free parameter: Shape of heat release (best result from Vince's experiment - he tried dComb=0.5,0.75,1.0,1.25,1.5)

		heatReleaseParams.xCombStart = GeomParams.xInlet + 0.057972; //Combustion ignition starting position = 0.408
		heatReleaseParams.xCombEnd   = GeomParams.xMax;

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
				throw(YANG_ERROR_CODE);
			}

			initHeatReleaseFunc();
		}

		// mAir
		if(firstCall) {
			double mAirInlet, mAirOutlet;
			double rhouAirInlet, rhouAirOutlet; // Actually these are momentum fluxes.
			getFluxesInletOutlet(mAirInlet, mAirOutlet, rhouAirInlet, rhouAirOutlet);
			if(mpi_rank==0) {
				cout<<"> mAir in the simulation:  Engine inlet = "<<mAirInlet<<",  Engine outlet = "<<mAirOutlet<<endl;
				double mAirAvg = 0.5*(mAirInlet + mAirOutlet);
				if( fabs((mAirAvg - hyShotHeatReleaseModel.mAir) / hyShotHeatReleaseModel.mAir) > 0.01)
					cout<<"WARNING! mAir in the simulation is not same as mAir in the heat-release model(="<<hyShotHeatReleaseModel.mAir<<")"<<endl;
				cout<<endl;
			}
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
		double mySumHeatRelease;

		if(lambda != NULL) {
			for (int icv=0; icv<ncv; ++icv) {
				if (x_cv[icv][0] >= hyShotHeatReleaseModel.xCombStart) {
					rhs_rhoE[icv] += lambda[0]*hyShotHeatReleaseModel.hFuel*heatReleaseFunc[icv] * cv_volume[icv];

					mySumHeatRelease += lambda[0]*hyShotHeatReleaseModel.hFuel*heatReleaseFunc[icv] * cv_volume[icv];
				}
			}
		}

		MPI_Allreduce(&mySumHeatRelease, &sumHeatRelease, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
	}


	/*
	 * Method: writeQoIOnFile
	 * ----------------------
	 * Averaged density and Averaged Mach number
	 */
	void writeQoIOnFile(const int itest, char filename[], bool rewrite) {
        static bool firstCall = true;

		if(UgpWithCvCompFlow::diverg==NULL ) {
			cerr<<classID<<"::writeQoIOnFile(): diverg has not been assigned"<<endl;
			throw(-1);
		}

		if(UgpWithCvCompFlow::grad_rho == NULL && mpi_rank==0)
			cout<<"WARNING "<<classID<<"::writeQoIOnFile(): grad_rho is NULL"<<endl;

		calcGradRhoIfNecessary();

		//=================
		// Calculate QoI's
		//=================
		// 1 - Shock impingement point on the upperwall
		string faZoneName = "upperwall";

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
		double shockLocUpperwall;
		double shockStrengthUpperwall;
		MPI_Allreduce(&myShockStrengthUpperwall, &shockStrengthUpperwall, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
		if(myShockStrengthUpperwall < shockStrengthUpperwall)
			myShockLocUpperwall = ABSURDLY_BIG_NUMBER;
		MPI_Allreduce(&myShockLocUpperwall, &shockLocUpperwall, 1, MPI_DOUBLE, MPI_MIN, mpi_comm);

		// 2 - Volume-averaged density, Volume-averaged density gradient, Subsonic volume portion
		double avgRho,  avgRhoGradMag,  subsonicPortion;
		getEngineIntegQuantities(avgRho, avgRhoGradMag, subsonicPortion);

		// 3 - Wall-pressures at the throat
		double upperwallPress, lowerwallPress;
		getThroatWallPress(upperwallPress, lowerwallPress);

		// 4 - Thrust of the engine
		double thrust_X, thrust_Y;
		get2DengineThrustQuantities(thrust_X, thrust_Y);

		//=================
		// Write QoI's on the file
		//=================
		if(mpi_rank==0) {
			FILE *fp;
			if(rewrite) {
				fp = fopen(filename, "w");
				fprintf(fp, "ITEST,               PHI,  SHOCK_LOC_UP,       AVG_RHO,  AVG_RHO_GRAD,  THR_P_UPWALL, THR_P_DWNWALL,  SUBSONIC_POR,      THRUST_X,      THRUST_Y,           PHI_ADJ,  THRUST_X_ADJ,  THRUST_Y_ADJ\n");
			} else
				fp = fopen(filename, "a");

			double LAMBDA_ADJ_COEFF = 7.939916e-9;
			double THRUST_ADJ_COEFF = 7.5;

			fprintf(fp, "%5d, %17.10e, %13.6e, %13.6e, %13.6e, %13.6e, %13.6e, %13.6e, %13.6e, %13.6e, %17.10e, %13.6e, %13.6e\n",
					itest, lambda[0], shockLocUpperwall, avgRho, avgRhoGradMag, upperwallPress, lowerwallPress,
					subsonicPortion, thrust_X, thrust_Y,
					lambda[0]*LAMBDA_ADJ_COEFF, thrust_X*THRUST_ADJ_COEFF, thrust_Y*THRUST_ADJ_COEFF);

			fclose(fp);
		}

		MPI_Barrier(mpi_comm);

        //=======================================================
        // Write the local subsonic ratio distribution on a file
        //=======================================================
        vector<double> xDistInCombustor;
        vector<double> localSubsonicRatio;
        const bool writeCountOnFile = firstCall;
        getLocalSubsonicRatio(xDistInCombustor, localSubsonicRatio, writeCountOnFile);
        assert(localSubsonicRatio.size() == localSubsonicRatio.size());

        if(mpi_rank ==0) {
        	//
        	for(size_t i=0; i<localSubsonicRatio.size(); ++i) {
        		if(localSubsonicRatio[i] >= 0.99) {
        			cout<<endl
        			    <<classID<<"::writeQoIOnFile(): Full subsonic location occurs!"<<endl
        			    <<endl;
        			break;
        		}
        	}

        	//
            char filename[70];
            sprintf(filename, "./DATASET_SUBSONIC_RATIO_DIST/SUBSONIC_RATIO_DIST%03d.csv", itest);
            FILE* fp = fopen(filename, "w");

            fprintf(fp, "X_LOC, SUBSONIC_RATIO\n");
            for(size_t i=0; i<localSubsonicRatio.size(); ++i)
                fprintf(fp, "%e, %e\n", xDistInCombustor[i], localSubsonicRatio[i]);

            fclose(fp);
        }

        //=======================================================
        // Write the wall pressure distribution on a file
        //=======================================================
        vector<double> lowerwallPressDist;
        vector<double> upperwallPressDist;
        getWallPressDist(xDistInCombustor, lowerwallPressDist, "lowerwall");
        getWallPressDist(xDistInCombustor, upperwallPressDist, "upperwall");

        assert(xDistInCombustor.size() == lowerwallPressDist.size());
        assert(xDistInCombustor.size() == upperwallPressDist.size());

        if(mpi_rank ==0) {
        	{
        		char filename[70];
        		sprintf(filename, "./DATASET_WALLPRESS_DIST_DWNWALL/PRESSDIST_DWNWALL%03d.csv", itest);
        		FILE* fp = fopen(filename, "w");

        		fprintf(fp, "X_LOC, PRESS\n");
        		for(size_t i=0; i<lowerwallPressDist.size(); ++i)
        			fprintf(fp, "%e, %e\n", xDistInCombustor[i], lowerwallPressDist[i]);

        		fclose(fp);
        	}

        	{
        		char filename[70];
        		sprintf(filename, "./DATASET_WALLPRESS_DIST_UPWALL/PRESSDIST_UPWALL%03d.csv", itest);
        		FILE* fp = fopen(filename, "w");

        		fprintf(fp, "X_LOC, PRESS\n");
        		for(size_t i=0; i<upperwallPressDist.size(); ++i)
        			fprintf(fp, "%e, %e\n", xDistInCombustor[i], upperwallPressDist[i]);

        		fclose(fp);
        	}
        }

		MPI_Barrier(mpi_comm);

        firstCall= false;
	}

	/****************************
	 *---  UTILITY FUNCTION  ---*
	 ****************************/
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
	 * Method: getFluxesInletOutlet
	 * ----------------------------
	 * Calculate mAir and rhouAir
	 */
	void getFluxesInletOutlet(double &mAirInlet, double &mAirOutlet, double &rhouAirInlet, double &rhouAirOutlet) {
		int INLET_OUTLET_FACE_COUNTS = 130;

		// Inlet
		double dxNearInletMin = 0.75 * 0.000028;
		double dxNearInletMax = 0.75 * 0.000040;
		double xInletCvMax = hyShotGeometry.xInlet + dxNearInletMax;

		double mySum_mAirInlet    = 0.0;
		double mySum_rhouAirInlet = 0.0;
		int myCountInletCv = 0;

		for(int icv=0; icv<ncv; ++icv) {
			bool possiblyInletCv = false;
			if(x_cv[icv][0] >= hyShotGeometry.xInlet && x_cv[icv][0] <= xInletCvMax) {
				if(x_cv[icv][1] > hyShotGeometry.yLowerwall && x_cv[icv][1] < hyShotGeometry.yUpperwall)
					possiblyInletCv = true;
			}

			if(possiblyInletCv) {
				// Check if the CV has a face attached to the inlet surface
				int nofa_f = faocv_i[icv];
				int nofa_l = faocv_i[icv+1]-1;
				for (int foc = nofa_f; foc <= nofa_l; foc++) {
					int ifa = faocv_v[foc];

					if(fabs(x_fa[ifa][0]-hyShotGeometry.xInlet) <= dxNearInletMin) {
						// face unit normal and area...
						double nVec[3] = {0.0, 0.0, 0.0};
						double area    = normVec3d(nVec, fa_normal[ifa]);
						assert(area > 0.0);

						if (nVec[0] > 0.966 || nVec[0] < -0.966) {  // Note: angle < 15 degrees
							mySum_mAirInlet    += area * rhou[icv][0];
							mySum_rhouAirInlet += area * rhou[icv][0]*fabs(rhou[icv][0]/rho[icv]); // momentum flux = Area * rho * u^2 * direct
							++myCountInletCv;
							break;
						}
					}
				}
			}
		}

		// Outlet
		double mySum_mAirOutlet    = 0.0;
		double mySum_rhouAirOutlet = 0.0;
		int myCountOutletFa = 0;

		for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
			if (zone->getKind() == FA_ZONE_BOUNDARY) {
				string tempFaZoneName = zone->getNameString();

				if (tempFaZoneName.compare("outlet") == 0) {
					for (int index = 0; index < zone->faVec.size(); ++index) {
						int ifa = zone->faVec[index];
						int icv0 = cvofa[ifa][0];

						if(icv0 < ncv) {
							if(x_fa[ifa][0] > hyShotGeometry.xThroat) {
								// face unit normal and area...
								double nVec[3] = {0.0, 0.0, 0.0};
								double area    = normVec3d(nVec, fa_normal[ifa]);
								assert(area > 0.0);

								if (nVec[0] > 0.866 || nVec[0] < -0.866) {  // Note: angle < 30 degrees
									mySum_mAirOutlet    += area * rhou[icv0][0];
									mySum_rhouAirOutlet += area * rhou[icv0][0]*fabs(rhou[icv0][0]/rho[icv0]); // momentum flux = Area * rho * u^2 * direct
									++myCountOutletFa;
								}
							}
						}
					}
				}
			}
		}

		// MPI sum
		int totCountInletCv, totCountOutletFa;
		MPI_Allreduce(&mySum_mAirInlet,  &mAirInlet,  1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(&mySum_mAirOutlet, &mAirOutlet, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(&mySum_rhouAirInlet,  &rhouAirInlet,  1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(&mySum_rhouAirOutlet, &rhouAirOutlet, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(&myCountInletCv,  &totCountInletCv,  1, MPI_INT, MPI_SUM, mpi_comm);
		MPI_Allreduce(&myCountOutletFa, &totCountOutletFa, 1, MPI_INT, MPI_SUM, mpi_comm);

		if(totCountInletCv != INLET_OUTLET_FACE_COUNTS && mpi_rank == 0)
			cout<<"WARNING "<<classID<<"::getFluxesInletOutlet(): totCountInletCv(="<<totCountInletCv<<") is not equal to "<<INLET_OUTLET_FACE_COUNTS<<endl;
		if(totCountOutletFa != INLET_OUTLET_FACE_COUNTS && mpi_rank == 0)
			cout<<"WARNING "<<classID<<"::getFluxesInletOutlet(): totCountOutletFa(="<<totCountOutletFa<<") is not equal to "<<INLET_OUTLET_FACE_COUNTS<<endl;
	}

	/*
	 * Method: getEngineIntegQuantities
	 * --------------------------------
	 * Calculate QoI's integral quantities (calculated in the engine): avgRho, avgRhoGradMag, subsonicPortion
	 */
	void getEngineIntegQuantities(double &avgRho, double &avgRhoGradMag, double &subsonicPortion) {
		assert(UgpWithCvCompFlow::grad_rho != NULL);

		double myAvgRho = 0.0, myAvgRhoGradMag = 0.0, mySubsonicPortion = 0.0;
		double myVolume = 0.0;

		for(int icv=0; icv<ncv; ++icv) {
			bool inTheEngine = false;
			if(x_cv[icv][0] >= hyShotGeometry.xInlet) {
				if(x_cv[icv][0] > hyShotGeometry.xThroat) {
					if(x_cv[icv][0] <= hyShotGeometry.xMax)
						inTheEngine = true;
				} else {
					if(x_cv[icv][1] >= hyShotGeometry.yLowerwall && x_cv[icv][1] <= hyShotGeometry.yUpperwall)
						inTheEngine = true;
				}
			}

			if(inTheEngine) {
				myAvgRho += UgpWithCvCompFlow::rho[icv] * cv_volume[icv];

				if(UgpWithCvCompFlow::grad_rho != NULL)
					myAvgRhoGradMag += sqrt(vecDotVec3d(UgpWithCvCompFlow::grad_rho[icv], UgpWithCvCompFlow::grad_rho[icv])) * cv_volume[icv];

				double uMag = sqrt(vecDotVec3d(UgpWithCvCompFlow::vel[icv], UgpWithCvCompFlow::vel[icv]));
				if(uMag / UgpWithCvCompFlow::sos[icv] < 1.0)
					mySubsonicPortion += cv_volume[icv];

				myVolume += cv_volume[icv];
			}
		}

		double totVolume;
		MPI_Allreduce(&myAvgRho,          &avgRho,          1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(&myAvgRhoGradMag,   &avgRhoGradMag,   1, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(&mySubsonicPortion, &subsonicPortion, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

		MPI_Allreduce(&myVolume, &totVolume, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

		avgRho          /= totVolume;
		avgRhoGradMag   /= totVolume;
		subsonicPortion /= totVolume;
	}

	/*
	 * Method: getThroatWallPress
	 * --------------------------
	 * Obtain the throat pressure on the lower wall and the upper wall
	 */
	void getThroatWallPress(double &upperwallPress, double &lowerwallPress) {
		double myUpperwallPress = 0.0;
		double myLowerwallPress = 0.0;
		int myUpperwallFoundCount = 0;
		int myLowerwallFoundCount = 0;

		double dxNearThroat = 0.75 * 0.0002;
		for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
			if (zone->getKind() == FA_ZONE_BOUNDARY) {
				string tempFaZoneName = zone->getNameString();

				if (tempFaZoneName.compare("upperwall") == 0) {
					for (int index = 0; index < zone->faVec.size(); ++index) {
						int ifa = zone->faVec[index];
						int icv0 = cvofa[ifa][0];

						if(icv0 < ncv) {
							if(x_fa[ifa][0] < hyShotGeometry.xThroat && fabs(x_fa[ifa][0] - hyShotGeometry.xThroat) < dxNearThroat) {
								myUpperwallPress = UgpWithCvCompFlow::press[icv0];
								++myUpperwallFoundCount;
							}
						}
					}
				}
				if (tempFaZoneName.compare("lowerwall") == 0) {
					for (int index = 0; index < zone->faVec.size(); ++index) {
						int ifa = zone->faVec[index];
						int icv0 = cvofa[ifa][0];

						if(icv0 < ncv) {
							if(x_fa[ifa][0] < hyShotGeometry.xThroat && fabs(x_fa[ifa][0] - hyShotGeometry.xThroat) < dxNearThroat) {
								myLowerwallPress = UgpWithCvCompFlow::press[icv0];
								++myLowerwallFoundCount;
							}
						}
					}
				}
			}
		}

		int totUpperwallFoundCount;
		int totLowerwallFoundCount;
		MPI_Allreduce(&myUpperwallPress, &upperwallPress, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
		MPI_Allreduce(&myLowerwallPress, &lowerwallPress, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
		MPI_Allreduce(&myUpperwallFoundCount, &totUpperwallFoundCount, 1, MPI_INT, MPI_SUM, mpi_comm);
		MPI_Allreduce(&myLowerwallFoundCount, &totLowerwallFoundCount, 1, MPI_INT, MPI_SUM, mpi_comm);

		if(totUpperwallFoundCount != 1 && mpi_rank == 0)
			cout<<"WARNING "<<classID<<"::getThroatWallPress(): "<<totUpperwallFoundCount<<" face(s) is(are) found on the upperwall"<<endl;
		if(totLowerwallFoundCount != 1 && mpi_rank == 0)
			cout<<"WARNING "<<classID<<"::getThroatWallPress(): "<<totLowerwallFoundCount<<" face(s) is(are) found on the lowerwall"<<endl;
	}

	/*
	 * Method: get2DengineThrustQuantities
	 * -----------------------------------
	 * Calculate the quantities related to the engine thrust
	 * Note: the given arguments must have only 2 elements (2D: X-Y)
	 */
	void get2DengineThrustQuantities(double &thrust_X, double &thrust_Y) {
		int INLET_OUTLET_FACE_COUNTS = 130;

		double dxNearInletMin = 0.75 * 0.000028;
		double dxNearInletMax = 0.75 * 0.000040;
		double xInletCvMax = hyShotGeometry.xInlet + dxNearInletMax;

		double mySumRhouInlet[2]  = {0.0, 0.0};
		double mySumPressInlet[2] = {0.0, 0.0};
		int myCountInletCv = 0;

		for(int icv=0; icv<ncv; ++icv) {
			bool possiblyInletCv = false;
			if(x_cv[icv][0] >= hyShotGeometry.xInlet && x_cv[icv][0] <= xInletCvMax) {
				if(x_cv[icv][1] > hyShotGeometry.yLowerwall && x_cv[icv][1] < hyShotGeometry.yUpperwall)
					possiblyInletCv = true;
			}

			if(possiblyInletCv) {
				// Check if the CV has a face attached to the inlet surface
				int nofa_f = faocv_i[icv];
				int nofa_l = faocv_i[icv+1]-1;
				for (int foc = nofa_f; foc <= nofa_l; foc++) {
					int ifa = faocv_v[foc];

					if(fabs(x_fa[ifa][0]-hyShotGeometry.xInlet) <= dxNearInletMin) {
						// face unit normal and area...
						double nVec[3] = {0.0, 0.0, 0.0};
						double area    = normVec3d(nVec, fa_normal[ifa]);
						assert(area > 0.0);

						if (nVec[0] > 0.966 || nVec[0] < -0.966) {  // Note: angle < 15 degrees
							for(int i=0; i<2; ++i)
								mySumRhouInlet[i] += area*rhou[icv][i]*fabs(rhou[icv][i]/rho[icv]);  // momentum flux = Area * rho * u^2 * direct
							for(int i=0; i<2; ++i)
								mySumPressInlet[i] += area*fabs(nVec[i])*UgpWithCvCompFlow::press[icv];
							++myCountInletCv;

							break;
						}
					}
				}
			}
		}

		double mySumRhouOutlet[2]  = {0.0, 0.0};
		double mySumPressOutlet[2] = {0.0, 0.0};
		int myCountOutletFa = 0;

		for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
			if (zone->getKind() == FA_ZONE_BOUNDARY) {
				string tempFaZoneName = zone->getNameString();

				if (tempFaZoneName.compare("outlet") == 0) {
					for (int index = 0; index < zone->faVec.size(); ++index) {
						int ifa = zone->faVec[index];
						int icv0 = cvofa[ifa][0];

						if(icv0 < ncv) {
							if(x_fa[ifa][0] > hyShotGeometry.xThroat) {
								// face unit normal and area...
								double nVec[3] = {0.0, 0.0, 0.0};
								double area    = normVec3d(nVec, fa_normal[ifa]);
								assert(area > 0.0);

								if (nVec[0] > 0.866 || nVec[0] < -0.866) {  // Note: angle < 30 degrees
									for(int i=0; i<2; ++i)
										mySumRhouOutlet[i] += area*rhou[icv0][i]*fabs(rhou[icv0][i]/rho[icv0]);  // momentum flux = Area * rho * u^2 * direct
									for(int i=0; i<2; ++i)
										mySumPressOutlet[i] += area*fabs(nVec[i])*UgpWithCvCompFlow::press[icv0]; // Here, it is assumed that the faces on the outlet face to the same direction.
									++myCountOutletFa;
								}
							}
						}
					}
				}
			}
		}

		double totSumRhouInlet[2];  	double totSumRhouOutlet[2];
		double totSumPressInlet[2]; 	double totSumPressOutlet[2];
		int totCountInletCv, totCountOutletFa;
		MPI_Allreduce(mySumRhouInlet,  totSumRhouInlet,  2, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(mySumRhouOutlet, totSumRhouOutlet, 2, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(mySumPressInlet,  totSumPressInlet,  2, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(mySumPressOutlet, totSumPressOutlet, 2, MPI_DOUBLE, MPI_SUM, mpi_comm);
		MPI_Allreduce(&myCountInletCv,  &totCountInletCv,  1, MPI_INT, MPI_SUM, mpi_comm);
		MPI_Allreduce(&myCountOutletFa, &totCountOutletFa, 1, MPI_INT, MPI_SUM, mpi_comm);

		if(totCountInletCv != INLET_OUTLET_FACE_COUNTS && mpi_rank == 0)
			cout<<"WARNING "<<classID<<"::get2DengineThrustQuantities(): totCountInletCv(="<<totCountInletCv<<") is not equal to "<<INLET_OUTLET_FACE_COUNTS<<endl;
		if(totCountOutletFa != INLET_OUTLET_FACE_COUNTS && mpi_rank == 0)
			cout<<"WARNING "<<classID<<"::get2DengineThrustQuantities(): totCountOutletFa(="<<totCountOutletFa<<") is not equal to "<<INLET_OUTLET_FACE_COUNTS<<endl;

		// Calculate thrust
		thrust_X = totSumRhouOutlet[0]-totSumRhouInlet[0] + totSumPressOutlet[0]-totSumPressInlet[0];
		thrust_Y = totSumRhouOutlet[1]-totSumRhouInlet[1] - totSumPressOutlet[1]+totSumPressInlet[1];
	}

    /*
     * Method: getLocalSubsonicRatio
     * -----------------------------
     * Get the ratio of subsonic volume at each x location
     */
    void getLocalSubsonicRatio(vector<double>& xDistInCombustor, vector<double>& localSubsonicRatio, const bool writeCountOnFile) {
        string SUB_RATIO_CV_COUNT_FILE_NAME = "SUB_RATIO_CV_COUNT.csv";
        double dxCell = 0.0002;
        double xMin = hyShotGeometry.xInlet;
        double xMax = hyShotGeometry.xThroat;

        int Nbins = (int) ceil ((xMax - xMin) / dxCell);

        xDistInCombustor.resize(Nbins, 0.0);
        localSubsonicRatio.resize(Nbins, 0.0);

        for(size_t i=0; i<Nbins; ++i)
            xDistInCombustor[i] = xMin + 0.5*dxCell + dxCell*i;

        //
        double* myLocalSubsonicVolume = new double [Nbins];
        double* myLocalVolume = new double [Nbins];
        int* myCount = new int [Nbins];

        for (int i=0; i<Nbins; ++i) {
            myLocalSubsonicVolume[i] = 0.0;
            myLocalVolume[i] = 0.0;
            myCount[i] = 0;
        }

        //
        for(int icv=0; icv<ncv; ++icv) {
            bool inTheEngine = false;
            if(x_cv[icv][0] >= xMin && x_cv[icv][0] <= hyShotGeometry.xThroat) {
                if(x_cv[icv][1] >= hyShotGeometry.yLowerwall && x_cv[icv][1] <= hyShotGeometry.yUpperwall)
                    inTheEngine = true;
            } else if (x_cv[icv][0] > hyShotGeometry.xThroat && x_cv[icv][0] < xMax) {
                inTheEngine = true;
            }

            if(inTheEngine) {
                int binIndex = int( (x_cv[icv][0] - xMin) / dxCell );
                assert(binIndex >= 0 && binIndex < Nbins);

                double M = sqrt(vecDotVec3d(UgpWithCvCompFlow::vel[icv], UgpWithCvCompFlow::vel[icv])) / UgpWithCvCompFlow::sos[icv];
                if(M < 1.0) {
                    myLocalSubsonicVolume[binIndex] += cv_volume[icv];
                }
                myLocalVolume[binIndex] += cv_volume[icv];
                ++myCount[binIndex];
            }
        }

        //
        double* localSubsonicVolume = new double [Nbins];
        double* localVolume = new double [Nbins];
        int* count = new int [Nbins];
        MPI_Allreduce(myLocalSubsonicVolume, localSubsonicVolume, Nbins, MPI_DOUBLE, MPI_SUM, mpi_comm);
        MPI_Allreduce(myLocalVolume,         localVolume,         Nbins, MPI_DOUBLE, MPI_SUM, mpi_comm);
        MPI_Allreduce(myCount,               count,               Nbins, MPI_INT,    MPI_SUM, mpi_comm);

        //
        for(int i=0; i<Nbins; ++i)
            localSubsonicRatio[i] = localSubsonicVolume[i] / localVolume[i];

        //
        if(writeCountOnFile && mpi_rank==0) {
            cout<<endl
                <<classID<<"::getLocalSubsonicRatio(): write the number of CVs in each bins on "<<SUB_RATIO_CV_COUNT_FILE_NAME<<" -- TOT # OF BINS = "<<Nbins<<endl
                <<endl;

            FILE* fp = fopen(SUB_RATIO_CV_COUNT_FILE_NAME.c_str(), "w");
            fprintf(fp, "INDEX, X_LOC, CV_COUNT\n");
            for(int i=0; i<Nbins; ++i)
                fprintf(fp, "%d, %e, %d\n", i, xDistInCombustor[i], count[i]);

            fclose(fp);
        }

        if(mpi_rank == 0) {
            for(int i=0; i<Nbins; ++i)
                if(count[i] == 0)
                    cout<<"WARNING "<<classID<<"::getLocalSubsonicRatio(): 0 CV found in the "<<i<<"th bin: x_bin="<<xDistInCombustor[i]<<endl;
        }

        delete [] myLocalSubsonicVolume;
        delete [] myLocalVolume;
        delete [] myCount;

        delete [] localSubsonicVolume;
        delete [] localVolume;
        delete [] count;

        MPI_Barrier(mpi_comm);
    }

    /*
      * Method: getWallPressDist
      * ------------------------
      * Get the pressure distribution along a wall
      */
     void getWallPressDist(vector<double>& xDistInCombustor, vector<double>& wallPressDist, const char wallName[]) {
         double dxCell = 0.0002;
         double dyTol  = 1.0e-7;
         double xMin = hyShotGeometry.xInlet;
         double xMax = hyShotGeometry.xThroat;

         int Nbins = (int) ceil ((xMax - xMin) / dxCell);

         xDistInCombustor.resize(Nbins, 0.0);
         wallPressDist.resize(Nbins, 0.0);

         for(size_t i=0; i<Nbins; ++i)
             xDistInCombustor[i] = xMin + 0.5*dxCell + dxCell*i;

         //
         double* myLocalPressAreaSum = new double [Nbins];
         double* myLocalAreaSum      = new double [Nbins];
         int* myCount = new int [Nbins];

         for (int i=0; i<Nbins; ++i) {
         	myLocalPressAreaSum[i] = 0.0;
         	myLocalAreaSum[i]      = 0.0;
             myCount[i] = 0;
         }

         //
 		for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
 			if (zone->getKind() == FA_ZONE_BOUNDARY) {
 				string tempFaZoneName = zone->getNameString();

 				if (tempFaZoneName.compare(wallName) == 0) {
 					for (int index = 0; index < zone->faVec.size(); ++index) {
 						int ifa = zone->faVec[index];
 						int icv0 = cvofa[ifa][0];

 						if(icv0 < ncv) {
 				            bool inTheEngine = false;
 				            if(x_fa[ifa][0] >= xMin && x_fa[ifa][0] <= hyShotGeometry.xThroat) {
 				                if(x_fa[ifa][1] >= hyShotGeometry.yLowerwall-dyTol && x_fa[ifa][1] <= hyShotGeometry.yUpperwall+dyTol)
 				                    inTheEngine = true;
 				            } else if (x_fa[ifa][0] > hyShotGeometry.xThroat && x_fa[ifa][0] < xMax) {
 				                inTheEngine = true;
 				            }

 				            if(inTheEngine) {
 				                int binIndex = int( (x_fa[ifa][0] - xMin) / dxCell );
 				                if(binIndex < 0 || binIndex >= Nbins) {
 				                	cerr<<"ERROR "<<classID<<"::getWallPressDist(): binIndex = "<<binIndex<<" is not in the range of [0 ~ Nbins="<<Nbins<<"] @ x="<<x_fa[ifa][0]<<endl
 				                	    <<"      Note: xMin="<<xMin<<", xMax="<<xMax<<endl;
 				                	throw(YANG_ERROR_CODE);
 				                }

 								// face unit normal and area...
 								double nVec[3] = {0.0, 0.0, 0.0};
 								double area    = normVec3d(nVec, fa_normal[ifa]);
 								assert(area > 0.0);

 								//
 				                myLocalPressAreaSum[binIndex] += UgpWithCvCompFlow::press[icv0] * area;
 				                myLocalAreaSum[binIndex]      += area;
 				                ++myCount[binIndex];
 				            }
 						}
 					}
 				}
 			}
 		}

         //
         double* localPressAreaSum = new double [Nbins];
         double* localAreaSum      = new double [Nbins];
         int* count = new int [Nbins];
         MPI_Allreduce(myLocalPressAreaSum, localPressAreaSum, Nbins, MPI_DOUBLE, MPI_SUM, mpi_comm);
         MPI_Allreduce(myLocalAreaSum,      localAreaSum,      Nbins, MPI_DOUBLE, MPI_SUM, mpi_comm);
         MPI_Allreduce(myCount,             count,             Nbins, MPI_INT,    MPI_SUM, mpi_comm);

         //
         for(int i=0; i<Nbins; ++i)
         	wallPressDist[i] = localPressAreaSum[i] / localAreaSum[i];

         //
         if(mpi_rank == 0) {
             for(int i=0; i<Nbins; ++i)
                 if(count[i] == 0)
                     cout<<"WARNING "<<classID<<"::getWallPressDist(): 0 CV found in the "<<i<<"th bin: x_bin="<<xDistInCombustor[i]<<endl;
         }

         delete [] myLocalPressAreaSum;
         delete [] myLocalAreaSum;
         delete [] myCount;

         delete [] localPressAreaSum;
         delete [] localAreaSum;
         delete [] count;

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
	int run = 0;

	// set input name to default "Ike.in"
	char inputFileName[50];
	sprintf(inputFileName, "Yang.in");

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
		JoeWithModels *Yang;

		switch (run) {
			case 1:  Yang = new YangWithKOm(inputFileName); 	break;
			default:
				if (mpi_rank == 0)
					cerr << "WARNING: run number not available!" << endl;
				Yang = new YangWithKOm(inputFileName);
		}

		// run vanGogh
		Yang->run();

		MPI_Barrier(mpi_comm);
		if (mpi_rank == 0) {
			double wtime0 = wtime;
			wtime = MPI_Wtime();
			cout << " > total runtime [s]: " << wtime - wtime0 << endl;
		}

		// delete vanGogh (make sure memory is deallocated in destructors
		delete Yang;
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

