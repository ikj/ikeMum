/*
 * IkeWithModels.cpp
 *
 *  Created on: Jan 4, 2013
 *      Author: ikj
 */

#include "IkeWithModels.h"

/*
 * Method: init
 * ------------
 * Initialize member variables
 */
void IkeWithModels_AD::init() {
	// CSR struct
	nbocv2_i = NULL;

	// global icv index
	nbocv2_v_global = NULL;

	// BC information
//		btofa = NULL;
//		btofaScalar = NULL;
//		bvalScalar = NULL;
	zofa = NULL;

#ifdef USE_ARTIF_VISC_WITH_MEM_SAVING
	// These values are defined in the IkeUgpWithCvCompFlow.h file
	myArtifViscMagMin_AD =  ABSURDLY_BIG_NUMBER;
	myArtifViscMagMax_AD = -ABSURDLY_BIG_NUMBER;
	myArtifViscSavedStep = -1;
#endif

#ifdef USE_DT_OVER_VOL_SCALING
	// Local dt/dV
	local_dtOverVol = NULL;	registerScalar(local_dtOverVol, "local_dtOverVol", CV_DATA);
	for(int icv=0; icv<ncv; ++icv)
		local_dtOverVol[icv] = 1.0;
	updateCvDataG1G2(local_dtOverVol, REPLACE_DATA);
#endif
}

/*
 * Method: clear
 * -------------
 * Clear member variables
 */
void IkeWithModels_AD::clear() {
	// CSR struct
	if (nbocv2_i != NULL) {
		delete[] nbocv2_i;
		nbocv2_i = NULL;
	}
	if (!nbocv2_v.empty()) {
		nbocv2_v.clear();
	}

	if (nbocv2_v_global != NULL)
		delete [] nbocv2_v_global; nbocv2_v_global = NULL;

	// BC information
	//		if (btofa != NULL) {
	//			delete [] btofa;
	//			btofa = NULL;
	//		}
	//		if (btofaScalar != NULL) {
	//			int nScal = scalarTranspEqVector.size();
	//			for (int scal=0; scal<nScal; ++scal)
	//				delete [] btofaScalar[scal];
	//
	//			delete [] btofaScalar;
	//			btofaScalar = NULL;
	//		}
	//		if (bvalScalar != NULL) {
	//			int nScal = scalarTranspEqVector.size();
	//			for (int scal=0; scal<nScal; ++scal)
	//				delete [] bvalScalar[scal];
	//
	//			delete [] bvalScalar;
	//			bvalScalar = NULL;
	//		}
	if (zofa != NULL) {
		delete [] zofa;
		zofa = NULL;
	}
}


/*
 * Method: getReferenceParams
 * --------------------------
 * Get reference parameters from the input file
 * Note: If the reference value for a variable (except uRefVec and turbScalRefVec) cannot be found,
 *       1.0 will be returned (The value 1.0 is chosen for normalization).
 */
void IkeWithModels_AD::getReferenceParams(REF_FLOW_PARAMS& RefFlowParams) {
	RefFlowParams.rho_ref   = fabs(UgpWithCvCompFlow::rho_ref);
	RefFlowParams.press_ref = fabs(UgpWithCvCompFlow::p_ref);
	RefFlowParams.temp_ref  = fabs(UgpWithCvCompFlow::T_ref);
	RefFlowParams.gamma_ref = fabs(UgpWithCvCompFlow::GAMMA);
	RefFlowParams.Rgas_ref  = fabs(UgpWithCvCompFlow::R_gas);

	Param *pmy;
	if (getParam(pmy, "U_INITIAL")) {
		for(int i=0; i<3; ++i)
			RefFlowParams.u_ref[i] = fabs(pmy->getDouble(i+1));
	}

	RefFlowParams.calcRhouRefFromU();
	RefFlowParams.calcVelMagRef();
	RefFlowParams.calcRhoERef();

	if(scalarTranspEqVector.size() > 0) {
		if(RefFlowParams.nScal == 0) {
			RefFlowParams.nScal = scalarTranspEqVector.size();
			RefFlowParams.scalar_ref.resize(RefFlowParams.nScal, 1.0);

			if (getParam(pmy, "INITIAL_CONDITION_TURB")) {
				for(size_t iScal=0; iScal<scalarTranspEqVector.size(); ++iScal)
					RefFlowParams.scalar_ref[iScal] = fabs(pmy->getDouble(iScal+1));

				RefFlowParams.kine_ref = RefFlowParams.scalar_ref[0];
			} else {
				if(mpi_rank==0)
					cout<<"WARNING in IkeWithModels_AD::getReferenceParams(): Reference values for turbulence scalars were not found -- Just use 1.0 as reference values"<<endl;
			}
		}
	}
}

/*
 * Method: build2layerCSRstruct
 * ----------------------------
 * Build a bigger CSR structure which contains two-layered neighbors of each CV
 *   The results are nbocv2_v_global = same as nbocv2_v_global but 2-layered
 *                   nbocv2_i = same as nbocv_i but 2-layered
 *                   nbocv2_v = same as nbocv_v but 2-layered
 */
void IkeWithModels_AD::build2layerCSRstruct(bool addFakeCells) {
	// global icv index (original code: solveCoupledLinSysNS() in UgpWithCvCompFlow.h
	assert( nbocv2_v_global == NULL );
	nbocv2_v_global = new int[ncv_gg];
	for (int icv = 0; icv < ncv; icv++)
		nbocv2_v_global[icv] = cvora[mpi_rank] + icv;
	updateCvDataG1G2(nbocv2_v_global, REPLACE_DATA);

	// prep the CSR struct...
	assert( nbocv2_i == NULL );
	nbocv2_i = new int[ncv + 1];
	for (int i = 0; i < ncv + 1; ++i)
		nbocv2_i[i] = 0;

	assert( nbocv2_v.empty() );

	// build the struct
	for (int icv = 0; icv < ncv; ++icv) {
		int nNbrs = find2layerCSRstruct_eachIcv(icv, nbocv2_v, addFakeCells);
		nbocv2_i[icv + 1] = nbocv2_i[icv] + nNbrs;
	}
}

/*
 * Method: find2layerCSRstruct_eachIcv
 * -----------------------------------
 * Find all the neighbors in the two-layer of a given CV
 */
int IkeWithModels_AD::find2layerCSRstruct_eachIcv(int icvCenter, vector<int> &nbocv2_v, bool addFakeCells) {
	set<int> visited;

	// first, add icvCenter itself
	nbocv2_v.push_back(icvCenter);
	visited.insert(icvCenter);

	// second, add first-layer neighbors
	int noc_f = nbocv_all_i[icvCenter];
	int noc_l = nbocv_all_i[icvCenter + 1] - 1;
	for (int noc = noc_f + 1; noc <= noc_l; noc++) {
		int icvNbr = nbocv_all_v[noc];
		if (addFakeCells) {
			nbocv2_v.push_back(icvNbr);
			visited.insert(icvNbr);
		} else {
			if (icvNbr < ncv_gg) {
				nbocv2_v.push_back(icvNbr);
				visited.insert(icvNbr);
			}
		}
	}

	// finally, add second-layer neighbors
	for (int noc = noc_f + 1; noc <= noc_l; noc++) {
		int icvNbr = nbocv_all_v[noc];
//		if (addFakeCells) {
//			int noc_f2 = nbocv_all_i[icvNbr];
//			int noc_l2 = nbocv_all_i[icvNbr + 1] - 1;
//			for (int noc2 = noc_f2 + 1; noc2 <= noc_l2; noc2++) {
//				int icvNbr2 = nbocv_all_v[noc2];
//				if (visited.count(icvNbr2) == 0) {
//					nbocv2_v.push_back(icvNbr2);
//					visited.insert(icvNbr2);
//				}
//			}
//		} else {
//			if (icvNbr < ncv_gg) {
//				int noc_f2 = nbocv_all_i[icvNbr];
//				int noc_l2 = nbocv_all_i[icvNbr + 1] - 1;
//				for (int noc2 = noc_f2 + 1; noc2 <= noc_l2; noc2++) {
//					int icvNbr2 = nbocv_all_v[noc2];
//					if (icvNbr2 < ncv_gg && visited.count(icvNbr2) == 0) {
//						nbocv2_v.push_back(icvNbr2);
//						visited.insert(icvNbr2);
//					}
//				}
//			}
//		}
		if (icvNbr < ncv_gg) {
			int noc_f2 = nbocv_all_i[icvNbr];
			int noc_l2 = nbocv_all_i[icvNbr + 1] - 1;
			for (int noc2 = noc_f2 + 1; noc2 <= noc_l2; noc2++) {
				int icvNbr2 = nbocv_all_v[noc2];
				if (visited.count(icvNbr2) == 0) {
					if (addFakeCells) {
						nbocv2_v.push_back(icvNbr2);
						visited.insert(icvNbr2);
					} else {
						if (icvNbr2 < ncv_gg) {
							nbocv2_v.push_back(icvNbr2);
							visited.insert(icvNbr2);
						}
					}
				}
			}
		}
	}

	return visited.size();
}

///*
// * Method: setBtofa
// * ----------------
// * Set the boundary type of faces (btofa, btofaScalar, bvalScalar)
// * The variables are first initialized as internal faces. Then, by iterating zones, a proper BC is applied to each face.
// */
//void IkeWithModels_AD::setBtofa() {
//	int nScal = scalarTranspEqVector.size();
//
//	assert(btofa == NULL);
//	assert(btofaScalar == NULL);
//	assert(bvalScalar == NULL);
//
//	// allocate memory
//	btofa = new BOUNDARY_TYPE [nfa_b2gg];
//	btofaScalar = new BOUNDARY_TYPE_SCALAR* [nScal];
//	for(int scal=0; scal<nScal; ++scal)
//		btofaScalar[scal] = new BOUNDARY_TYPE_SCALAR [nfa_b2gg];
//	bvalScalar = new double* [nScal];
//	for(int scal=0; scal<nScal; ++scal)
//		bvalScalar[scal] = new double [nfa_b2gg];
//
//	// initialize as an internal face
//	for (int ifa=0; ifa<nfa_b2gg; ++ifa)
//		btofa[ifa] = INTERNAL;
//	for (int scal=0; scal<nScal; ++scal)
//		for (int ifa=0; ifa<nfa_b2gg; ++ifa)
//			btofaScalar[scal][ifa] = INTERNAL_SCALAR;
//	for (int scal=0; scal<nScal; ++scal)
//		for (int ifa=0; ifa<nfa_b2gg; ++ifa)
//			bvalScalar[scal][ifa] = 0.0;
//
//	// update boundary type for boundary faces
//	// -- btofa
//	for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
//		if (zone->getKind() == FA_ZONE_BOUNDARY) {
//			Param *param;
//
//			if (getParam(param, zone->getName())) {
//				if (param->getString() == "HOOK") {
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//						btofa[ifa] = HOOK;
//					}
//				} else if (param->getString() == "CBC") {
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//						btofa[ifa] = CBC;
//					}
//				} else if (param->getString() == "CBC_SUBSONIC_INLET") {
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//						btofa[ifa] = CBC_SUBSONIC_INLET;
//					}
//				} else if (param->getString() == "CBC_SUBSONIC_OUTLET") {
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//						btofa[ifa] = CBC_SUBSONIC_OUTLET;
//					}
//				} else if (param->getString() == "SYMMETRY") {
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//						btofa[ifa] = SYMMETRY;
//					}
//				} else if (param->getString() == "NEUMANN") {
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//						btofa[ifa] = NEUMANN;
//					}
//				} else if (param->getString() == "WALL") {
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//						btofa[ifa] = WALL;
//					}
//				} else {
//					if (mpi_rank == 0)
//						cerr << "Error: cannot recognize bc for zone: " << zone->getName() << " bc: " << param->getString() << endl;
//				}
//			}
//		}
//
//	// -- btofaScalar
//	for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
//		for (int scal = 0; scal < nScal; scal++) {
//			string scalName = scalarTranspEqVector[scal].getName();
//			double dummy;
//
//			if (zone->getKind() == FA_ZONE_BOUNDARY) {
//				// HOOK BOUNDARY CONDITION
//				if (scalarZoneIsHook(zone->getName(), scalName)) {
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//						btofaScalar[scal][ifa] = HOOK_SCALAR;
//					}
//				}
//				// DIRICHLET BOUNDARY CONDITION
//				else if (scalarZoneIsDirichlet(dummy, zone->getName(), scalName)) {
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//						btofaScalar[scal][ifa] = DIRICHLET_SCALAR;
//						bvalScalar[scal][ifa] = dummy;
//					}
//				}
//				// FLUX BOUNDARY CONDITION (ACTUALLY, NOT IMPLEMENTED YET IN THE ORIGINAL JOE CODE)
//				else if (scalarZoneIsFlux(dummy, zone->getName(), scalName))
//				{
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//						btofaScalar[scal][ifa] = FLUX_SCALAR;
//						bvalScalar[scal][ifa] = dummy;
//					}
//				}
//				// OTHER BOUNDARY CONDITIONS: NO FLUX
//				else // apply zero flux if no boundary condition is specified
//				{
//					for (int index = 0; index < zone->faVec.size(); ++index) {
//						int ifa = zone->faVec[index];
//						btofaScalar[scal][ifa] = OTHER_SCALAR;
//					}
//				}
//			}
//		}
//}

/*
 * Method: setZofa
 * ---------------
 * Set up the zone type of faces (zofa): zofa is an array which stores the pointers to a FaZone at each face.
 */
void IkeWithModels_AD::setZofa() {
	int nScal = scalarTranspEqVector.size();

	// allocate memory
	assert(zofa == NULL);
	zofa = new FaZone* [nfa_b2gg];

	// initialize the array as NULL
	for (int ifa=0; ifa<nfa_b2gg; ++ifa)
		zofa[ifa] = NULL;

	// find the pointers
	for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
		if (zone->getKind() == FA_ZONE_BOUNDARY) {
			for (int index = 0; index < zone->faVec.size(); ++index) {
				int ifa = zone->faVec[index];
				zofa[ifa] = &(*zone);
			}
		}
	}
}

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
int IkeWithModels_AD::find1LayerFaceIndices(const int icvCenter, vector<int>& faInternal, vector<int>& faBoundary) {
	// Check if the 2-layer CSR structure has already been developed
	assert(faocv_i != NULL && faocv_v != NULL); // these are defined in Ugp

	// Clear faInternal and faBoundary
	faInternal.clear();		faBoundary.clear();

	// Find faces
	set<int> visited;

	int nofa_f = faocv_i[icvCenter];
	int nofa_l = faocv_i[icvCenter+1]-1;

	for (int nofa = nofa_f; nofa <= nofa_l; nofa++) {
		int ifa = faocv_v[nofa];

		if (visited.count(ifa) == 0) {
			if((ifa>=0 && ifa<nfa_b) || (ifa>=nfa && ifa<nfa_b2))
				faBoundary.push_back(ifa);
			else if((ifa>=nfa_b && ifa<nfa) || (ifa>=nfa_b2 && ifa<nfa_b2gg))
				faInternal.push_back(ifa);
			else {
				cerr<<"ERROR in IkeWithModels_AD::find1LayerFaceIndices() => cannot recognize the given ifa="<<ifa<<endl;
				printf("  (Note: nfa_b=%d, nfa_bpi=%d, nfa=%d, nfa_b2=%d, nfa_b2g=%d, nfa_b2gg=%d) \n",nfa_b, nfa_bpi, nfa, nfa_b2, nfa_b2g, nfa_b2gg);
				throw(-1);
			}
				 /* Note: 0 ~ nfa_b-1           : boundary faces
				  *       nfa_b ~ nfa_bpi-1     :
				  *       nfa_bpi ~ nfa-1       :
				  *       nfa ~ nfa_b2         : additional boundary faces between g1:f2 cvs
				  *       nfa_b2 ~ nfa_b2g-1   : faces between g1:g1 cvs
				  *       nfa_b2g-1 ~ nfa_b2gg-1: faces between g1:g2 cvs
				  */
			visited.insert(ifa);
		}
	}

	if(visited.size() == 0)
		printf("Warning in IkeWithModels_AD::find1LayerFaceIndices() => cannot find any face at icv=%d",icvCenter);

	return visited.size();
}
int IkeWithModels_AD::find1LayerFaceIndices(const int icvCenter, vector<int>& fa2) {
	// Check if the 2-layer CSR structure has already been developed
	assert(faocv_i != NULL && faocv_v != NULL); // these are defined in Ugp

	// Clear faInternal and faBoundary
	fa2.clear();

	// Find faces
	set<int> visited;

	int nofa_f = faocv_i[icvCenter];
	int nofa_l = faocv_i[icvCenter+1]-1;

	for (int nofa = nofa_f; nofa <= nofa_l; nofa++) {
		int ifa = faocv_v[nofa];

		if (visited.count(ifa) == 0) {
			if(ifa>=0 && nfa_b2gg)
				fa2.push_back(ifa);
			else {
				cerr<<"ERROR in IkeWithModels_AD::find1LayerFaceIndices() => cannot recognize the given ifa="<<ifa<<endl;
				printf("  (Note: nfa_b=%d, nfa_bpi=%d, nfa=%d, nfa_b2=%d, nfa_b2g=%d, nfa_b2gg=%d) \n",nfa_b, nfa_bpi, nfa, nfa_b2, nfa_b2g, nfa_b2gg);
				throw(-1);
			}
				 /* Note: 0 ~ nfa_b-1           : boundary faces
				  *       nfa_b ~ nfa_bpi-1     :
				  *       nfa_bpi ~ nfa-1       :
				  *       nfa ~ nfa_b2         : additional boundary faces between g1:f2 cvs
				  *       nfa_b2 ~ nfa_b2g-1   : faces between g1:g1 cvs
				  *       nfa_b2g-1 ~ nfa_b2gg-1: faces between g1:g2 cvs
				  */
			visited.insert(ifa);
		}
	}

	if(visited.size() == 0)
		printf("Warning in IkeWithModels_AD::find1LayerFaceIndices() => cannot find any face at icv=%d",icvCenter);

	return visited.size();
}

/*
 * Method: find2LayerFaceIndices
 * -----------------------------
 * Find all the faces in the two-layer neighbors.
 * For example, 16 faces from the total 50 faces will be found for the following 2D structure
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
int IkeWithModels_AD::find2LayerFaceIndices(const int icvCenter, vector<int>& faInternal, vector<int>& faBoundary) {
	// Check if the 2-layer CSR structure has already been developed
	assert(nbocv_all_i!= NULL && nbocv_all_v!= NULL); // these are defined in Ugp
	assert(faocv_i != NULL && faocv_v != NULL); // these are defined in Ugp

	// Clear faInternal and faBoundary
	faInternal.clear();		faBoundary.clear();

	// Find faces
	set<int> visited;

	int noc_f = nbocv_all_i[icvCenter];
	int noc_l = nbocv_all_i[icvCenter+1] - 1;
	for (int noc = noc_f; noc <= noc_l; noc++) {
		int icv_nbr = nbocv_all_v[noc];

		if(icv_nbr < ncv_g) {
			int nofa_f = faocv_i[icv_nbr];
			int nofa_l = faocv_i[icv_nbr+1]-1;

			for (int nofa = nofa_f; nofa <= nofa_l; nofa++) {
				int ifa = faocv_v[nofa];

				if (visited.count(ifa) == 0) {
					if((ifa>=0 && ifa<nfa_b) || (ifa>=nfa && ifa<nfa_b2))
						faBoundary.push_back(ifa);
					else if((ifa>=nfa_b && ifa<nfa) || (ifa>=nfa_b2 && ifa<nfa_b2gg))
						faInternal.push_back(ifa);
					else {
						cerr<<"ERROR in IkeWithModels_AD::find2LayerFaceIndices() => cannot recognize the given ifa="<<ifa<<" (face of icv="<<icv_nbr<<")"<<endl;
						printf("  (Note: nfa_b=%d, nfa_bpi=%d, nfa=%d, nfa_b2=%d, nfa_b2g=%d, nfa_b2gg=%d) \n",nfa_b, nfa_bpi, nfa, nfa_b2, nfa_b2g, nfa_b2gg);
						throw(-1);
					}
						 /* Note: 0 ~ nfa_b-1           : boundary faces
						  *       nfa_b ~ nfa_bpi-1     :
						  *       nfa_bpi ~ nfa-1       :
						  *       nfa ~ nfa_b2         : additional boundary faces between g1:f2 cvs
						  *       nfa_b2 ~ nfa_b2g-1   : faces between g1:g1 cvs
						  *       nfa_b2g-1 ~ nfa_b2gg-1: faces between g1:g2 cvs
						  */
					visited.insert(ifa);
				}
			}
		}
	}

	if(visited.size() == 0)
		printf("Warning in IkeWithModels_AD::find2LayerFaceIndices() => cannot find any face at icv=%d",icvCenter);

	return visited.size();
}

int IkeWithModels_AD::find2LayerFaceIndices(const int icvCenter, vector<int>& fa2) {
	// Check if the 2-layer CSR structure has already been developed
	assert(nbocv_all_i!= NULL && nbocv_all_v!= NULL); // these are defined in Ugp
	assert(faocv_i != NULL && faocv_v != NULL); // these are defined in Ugp

	// Clear faInternal and faBoundary
	fa2.clear();

	// Find faces
	set<int> visited;

	int noc_f = nbocv_all_i[icvCenter];
	int noc_l = nbocv_all_i[icvCenter+1] - 1;
	for (int noc = noc_f; noc <= noc_l; noc++) {
		int icv_nbr = nbocv_all_v[noc];

		if(icv_nbr < ncv_g) {
			int nofa_f = faocv_i[icv_nbr];
			int nofa_l = faocv_i[icv_nbr+1]-1;

			for (int nofa = nofa_f; nofa <= nofa_l; nofa++) {
				int ifa = faocv_v[nofa];

				if (visited.count(ifa) == 0) {
					if(ifa>=0 && ifa<nfa_b2gg)
						fa2.push_back(ifa);
					else {
						cerr<<"ERROR in IkeWithModels_AD::find2LayerFaceIndices() => cannot recognize the given ifa="<<ifa<<" (face of icv="<<icv_nbr<<")"<<endl;
						printf("  (Note: nfa_b=%d, nfa_bpi=%d, nfa=%d, nfa_b2=%d, nfa_b2g=%d, nfa_b2gg=%d) \n",nfa_b, nfa_bpi, nfa, nfa_b2, nfa_b2g, nfa_b2gg);
						throw(-1);
					}
						 /* Note: 0 ~ nfa_b-1           : boundary faces
						  *       nfa_b ~ nfa_bpi-1     :
						  *       nfa_bpi ~ nfa-1       :
						  *       nfa ~ nfa_b2         : additional boundary faces between g1:f2 cvs
						  *       nfa_b2 ~ nfa_b2g-1   : faces between g1:g1 cvs
						  *       nfa_b2g-1 ~ nfa_b2gg-1: faces between g1:g2 cvs
						  */
					visited.insert(ifa);
				}
			}
		}
	}

	if(visited.size() == 0)
		printf("Warning in IkeWithModels_AD::find2LayerFaceIndices() => cannot find any face at icv=%d",icvCenter);

	return visited.size();
}

#ifdef USE_DT_OVER_VOL_SCALING
/*
 * Method: calcDtOverVol
 * ---------------------
 * Calculate dt/cv_vol: u^{n+1}-u{n} = RHS * (dt/cv_vol)
 * "local_dtOverVol" will be used to scale RHS
 * Original code: UgpWithCvCompFlow::calcDt()
 */
void IkeWithModels_AD::calcDtOverVol(double &dtOverVol_min, double &dtOverVol_max, const double cfl_target) {
	double *press = UgpWithCvCompFlow::press;
	double *gamma = UgpWithCvCompFlow::gamma;
	double *mul_fa = UgpWithCvCompFlow::mul_fa;
	double *mut_fa = UgpWithCvCompFlow::mut_fa;

	double myDtOverVol_min = ABSURDLY_BIG_NUMBER;
	double myDtOverVol_max = -ABSURDLY_BIG_NUMBER;

	string scalingMethod = getStringParam("RHS_SCALING","CFL_MAX_BASED");

	if (scalingMethod == "CFL_MAX_BASED") { // The scaling is large when 1. small cell (~ 1/area), 2. small acoustic velocity (~ 1/(|U|+c))
		for (int icv = 0; icv < ncv; icv++) {
			double dt_cv = DT_OVER_VOL_MAX;
			double lambdaMax = 0.0;

			double c = sqrt(gamma[icv]*press[icv]/rho[icv]);

			int foc_f = faocv_i[icv];
			int foc_l = faocv_i[icv+1]-1;

			for (int foc=foc_f; foc<=foc_l; foc++) {
				int ifa = faocv_v[foc];

				double nVec[3];
				double area = normVec3d(nVec, fa_normal[ifa]);

				double Uk = vecDotVec3d(rhou[icv], nVec)/rho[icv];
				double lambda = (fabs(Uk) + c)*area;
				double lamdiff = 0.0;
				if (mu_ref > 0.0)
					lamdiff = 4.0 * (mul_fa[ifa] + mut_fa[ifa]) / rho[icv] * area * area / cv_volume[icv];
				lambda = max(lambda, lamdiff);
				lambdaMax = max(lambdaMax, lambda);
			}

			dt_cv = min(DT_OVER_VOL_MAX, cfl_target/lambdaMax); // Note: Original formulation is multiplied by cv_volume[icv]
			local_dtOverVol[icv] = dt_cv;

			myDtOverVol_min = min(myDtOverVol_min, dt_cv);
			myDtOverVol_max = max(myDtOverVol_max, dt_cv);
		}
	} else if (scalingMethod == "CFL_SUM_BASED") { // The scaling is large when 1. small cell (~ 1/area), 2. small acoustic velocity (~ 1/(|U|+c))
		for (int icv = 0; icv < ncv; icv++) {
			double dt_cv, dt_conv, dt_visc;
			double c = sqrt(gamma[icv]*press[icv]/rho[icv]);

			int foc_f = faocv_i[icv];
			int foc_l = faocv_i[icv+1]-1;

			double lambda_conv = 0.0;
			double lambda_diff = 0.0;
			for (int foc=foc_f; foc<=foc_l; foc++) {
				int ifa = faocv_v[foc];
				double nVec[3];
				double area = normVec3d(nVec, fa_normal[ifa]);
				double Uk = vecDotVec3d(rhou[icv], nVec)/rho[icv];
				double volume = 0.0;
				int icv0 = cvofa[ifa][0]; int icv1 = cvofa[ifa][1];
				if (( icv0 >= 0 && icv0<ncv_gg) && ( icv1 >= 0 && icv1 <ncv_gg )) volume = 0.5*(cv_volume[icv0] + cv_volume[icv1]);
				else volume = 0.5*cv_volume[icv];

				lambda_conv += (fabs(Uk) + c)*area;
				if (mu_ref > 0.0) {
					double KinVisc = (4.0/3.0) * (mul_fa[ifa] + mut_fa[ifa]) / rho[icv];
					lambda_diff +=  area * area * KinVisc / volume;
				}
				else
					lambda_diff = 0.0;
			}

			dt_conv = cfl_target / lambda_conv; // Note: Original formulation is multiplied by cv_volume[icv]
			if (mu_ref > 0.0) dt_visc = 0.25*cfl_target / lambda_diff; // Note: Original formulation is multiplied by cv_volume[icv]
			else dt_visc = DT_OVER_VOL_MAX;

			dt_cv = min(DT_OVER_VOL_MAX, min(dt_conv, dt_visc));

			local_dtOverVol[icv] = dt_cv;

			myDtOverVol_min = min(dtOverVol_min, dt_cv);
			myDtOverVol_max = max(dtOverVol_max, dt_cv);
		}
	}

	updateCvDataG1G2(local_dtOverVol, REPLACE_DATA);

	MPI_Allreduce(&myDtOverVol_min, &dtOverVol_min, 1, MPI_DOUBLE, MPI_MIN, mpi_comm);
	MPI_Allreduce(&myDtOverVol_max, &dtOverVol_max, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
}
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
int IkeWithModels_AD::calcResidual1D_AD(const int icvCenter,
		REALA &rhs_rho_AD, REALA rhs_rhou_AD[3],REALA &rhs_rhoE_AD, REALAS *rhs_rhoScal_AD,
		ADscalar<REALQ> &rho_AD, ADvector<REALQ> &rhou_AD, ADscalar<REALQ> &rhoE_AD, double (*A)[5][5], double ***AScal, int flagImplicit) {
	//--------------
	// find faces
	//--------------
	vector<int> faInternal;
	vector<int> faBoundary;
	find1LayerFaceIndices(icvCenter, faInternal, faBoundary);

	vector<int> fa2LayerInternal;
	vector<int> fa2LayerBoundary;
	find2LayerFaceIndices(icvCenter, fa2LayerInternal, fa2LayerBoundary);

	//--------------
	// calculate rhs
	//--------------
	calcStateVariables1D_AD(icvCenter, rho_AD, rhou_AD, rhoE_AD);
	calcMaterialProperties1D_AD(icvCenter, faInternal, faBoundary, rho_AD, rhou_AD, rhoE_AD);
//	setBC1D_AD(icvCenter, fa2LayerBoundary, fa2LayerBoundary, rho_AD, rhou_AD, rhoE_AD);
	setBC1D_AD(icvCenter, fa2LayerInternal, fa2LayerBoundary, rho_AD, rhou_AD, rhoE_AD);

	if(mu_ref>0.0 || sndOrder==true)
		calcCv2Grad1D_AD(icvCenter, grad_u, vel, limiterNavierS, sos, epsilonSDWLS);

//	if(mu_ref>0.0) {
		FaZone** zoneArrayBoundary = NULL;
		if(!faBoundary.empty()) {
			zoneArrayBoundary = new FaZone* [faBoundary.size()];
			for(size_t index=0; index<faBoundary.size(); ++index) {
				int ifa = faBoundary[index];
				zoneArrayBoundary[index] = zofa[ifa];
			}
		}
		calcRansTurbViscMuet1D_AD(icvCenter, zoneArrayBoundary, faInternal, faBoundary, rho_AD, rhou_AD);
		if(zoneArrayBoundary != NULL) {
			delete [] zoneArrayBoundary;
			zoneArrayBoundary = NULL;
		}
//	}

	int CountReducedOrder = calcRhs1D_AD(icvCenter, rhs_rho_AD, rhs_rhou_AD, rhs_rhoE_AD, rhs_rhoScal_AD, rho_AD, rhou_AD, rhoE_AD, A, AScal, flagImplicit);

	return CountReducedOrder;
}

#ifdef USE_ARTIF_VISC_WITH_MEM_SAVING
/*
 * Method: calcArtifVisc1D_AD
 * --------------------------
 * Original code = UgpWithCvCompFlow::calcArtifVisc()
 * Calculating artificial viscosity -- This method calculates the artifical viscosity at the cell center.
 *                                     If user DON'T want to use bulk-visc only, also interpolate it on the faces.
 *                                     (For the bulk-visc only mode, interpolation occurs later)
 * It takes the parameters only related to the artificial viscosity,
 * but it can access to the all the member variables of UgpWithCvCompFlows and UgpWithCvCompFlows_AD.
 */
void IkeWithModels_AD::calcArtifVisc1D_AD(const int icvCenter, ADscalar<REALQ> &artifViscMag_AD,
		const bool artifVisc_bulkViscOnly, const bool artifVisc_shockOnly,
		const double artifVisc_smoothstepThresh, const string &artifVisc_type, const double artifVisc_coeff) {
	static bool firstCall = true;

	assert(!diverg.empty() && !strMag.empty());
//  Note: Unlike UgpWithCvCompFlow::calcArtifVisc(), we will not use the grad_rho based grid-size calculation.
//	assert(!grad_rho.empty());
	assert(!artifViscMag_AD.empty());

	calcStrainRateAndDivergence1D_AD(icvCenter);

//	if(!sndOrder && !grad_rho.empty())
//		calcCv2Grad1D_AD(grad_rho, rho_AD, limiterNavierS, rho_AD, epsilonSDWLS);

	// Calculate artifViscMag_AD in the CV cells
	int noc_f = nbocv_i[icvCenter];
	int noc_l = nbocv_i[icvCenter+1] - 1;
	for (int noc = noc_f; noc <= noc_l; noc++) {
		int icv = nbocv_v[noc];
		// If SHOCK_ONLY, clip for expansion waves and weak shocks
		adouble clipFunc = 1.0;
		if(artifVisc_shockOnly) {
			if(diverg[icv].value() > 0.0) // expansion waves
				clipFunc = 0.0;
			else {
				if(diverg[icv].value() > artifVisc_smoothstepThresh) {
					adouble divergNormalized = -diverg[icv] / artifVisc_smoothstepThresh;
					clipFunc = 3.0*pow(divergNormalized, 2.0) - 2.0*pow(divergNormalized, 3.0);
					// A smooth-step function whose domain is between 0.0 and 1.0
					// (For details, search Wikipedia with "smoothstep")
				}
			}
		}
		if(clipFunc.value()<0.0 || clipFunc.value()>1.0)
			cout<<"WARNING IkeWithModels_AD::calcArtifVisc1D_AD(): The clip function is not in the correct range (0~1) = "<<clipFunc.value()<<endl;

		// Calculate grid-size.
		double delta = 0.0;

//		if(!grad_rho.empty()) {
//			double grad_rho_ugp[3];
//			for(int i=0; i<3; ++i)
//				grad_rho_ugp[i] = grad_rho[icvCenter][i].value();
//			double unitGradRho[3];
//			double gradRhoMag = normVec3d(unitGradRho, grad_rho_ugp);
//			if(gradRhoMag < 1.0e-12)  // Sometimes unitGradRho becomes NaN when gradRhoMag==0.0
//				for(int i=0; i<3; ++i)
//					unitGradRho[i] = 0.0;
//
//			int foc_f = faocv_i[icvCenter];
//			int foc_l = faocv_i[icvCenter+1]-1;
//			for (int foc=foc_f; foc<=foc_l; foc++) {
//				int ifa = faocv_v[foc];
//
//				double nVec[3];
//				double area = normVec3d(nVec, fa_normal[ifa]);
//				double dx = cv_volume[icv] / area; // Approximated grid-spacing in the face-normal direction: works only for a hex mesh
//
//				double dxProjected = fabs(dx * vecDotVec3d(unitGradRho, nVec));
//
//				if( dxProjected<0.0 || isnan(dxProjected) ) {
//					cerr<<"UgpWithCvCompFlow::calcMaterialProperties(): wrong dxProjected is calculated = "<<dxProjected<<endl;
//					throw(-1);
//				}
//
//				//					  delta = max(delta, dxProjected);  // MAX_BASED
//				delta += dxProjected;  // SUM_BASED
//			}
//
//			delta /= double(foc_l - foc_f + 1);
//		} else { // If grad_rho is not defined, then you must find another (but inaccurate) method
			double unitRhou[3];
			double rhouMag = normVec3d(unitRhou, rhou[icv]);
			assert(!isnan(rhouMag));
			int countFaces = 0;

			int foc_f = faocv_i[icv];
			int foc_l = faocv_i[icv+1]-1;
			for (int foc=foc_f; foc<=foc_l; foc++) {
				int ifa = faocv_v[foc];

				double nVec[3];
				double area = normVec3d(nVec, fa_normal[ifa]);

				double flowdirectInnerprod = fabs(vecDotVec3d(unitRhou, nVec));

				if(fabs(rhouMag) < 1.0e-8 || flowdirectInnerprod > 2.0e-3) { // Note: cos(0.1 degrees) = 0.0017
					double dx = cv_volume[icv] / area; // Approximated grid-spacing in the face-normal direction: works only for a hex mesh
					delta += dx;  // Note: MAX_BASE cannot work since it will only take the span-wise grid-spacing in most of the 2D calculations
					++countFaces;
				}
			}

			delta /= double(countFaces);
//		}

		// Calculate CV-based artificial viscosity
		if(artifVisc_type.compare("STRAIN_BASE") == 0) {
			artifViscMag_AD[icv] = artifVisc_coeff * rho_AD[icv] * pow(delta, 2.0) * fabs(strMag[icv]) * clipFunc;
			artifVisc_mag[icv]   = artifViscMag_AD[icv].value();
		} else if(artifVisc_type.compare("DIVERG_BASE") == 0) {
			artifViscMag_AD[icv] = artifVisc_coeff * rho_AD[icv] * pow(delta, 2.0) * fabs(diverg[icv]) * clipFunc;
			artifVisc_mag[icv]   = artifViscMag_AD[icv].value();
		} else {
			if(mpi_rank==0) cerr<<"IkeWithModels_AD::calcArtifVisc1D_AD(): artifVisc_type is not supported = "<<artifVisc_type<<endl;
			throw(-1);
		}

		if(step != myArtifViscSavedStep) {
			myArtifViscSavedStep = step;
			myArtifViscMagMin_AD =  ABSURDLY_BIG_NUMBER;
			myArtifViscMagMax_AD = -ABSURDLY_BIG_NUMBER;
		}
		myArtifViscMagMin_AD = min(myArtifViscMagMin_AD, artifViscMag_AD[icv].value());
		myArtifViscMagMax_AD = max(myArtifViscMagMax_AD, artifViscMag_AD[icv].value());

		// Check possible errors
		if( artifViscMag_AD[icv].value() < 0.0 || isnan(artifViscMag_AD[icv].value()) ) {
			cerr<<"IkeWithModels_AD::calcArtifVisc1D_AD(): artifVisc_mag becomes negative = "<<artifViscMag_AD[icv].value()<<endl;
			throw(-1);
		}
	}

	// Interpolate artifVisc_mag on the faces
	if(!artifVisc_bulkViscOnly) { // If the user DOES NOT want to add the artificial viscosity only on the bulk viscosity
		int nofa_f = faocv_i[icvCenter];
		int nofa_l = faocv_i[icvCenter+1]-1;
		for (int foc = nofa_f; foc <= nofa_l; foc++) {
			int ifa = faocv_v[foc];

			// internal faces
			if((ifa >= nfa_b && ifa<nfa) || (ifa >= nfa_b2 && ifa < nfa_b2gg)) {
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];

				double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
				vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
				vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
				double w0 = sqrt(vecDotVec3d(dx0, dx0));
				double w1 = sqrt(vecDotVec3d(dx1, dx1));

				if(mu_ref == 0)
					mul_fa[ifa] = 0.0;  // Make sure that the baseline viscosity is equal to zero for inviscid calculations

				mul_fa[ifa] += (w1*artifViscMag_AD[icv0] + w0*artifViscMag_AD[icv1]) / (w0+w1);
			}
		}

		// boundary faces computed next in setBC
	}

	firstCall = false;
}
#endif

/*
 * Method: setBC1D_AD
 * ------------------
 * Original code = setBC_AD in JoeWithModelsAD.cpp
 */
void IkeWithModels_AD::setBC1D_AD(const int icvCenter, vector<int>& faInternal, vector<int>& faBoundary, ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<REALQ> &rhoE) {
	static int first = 1;
	int bc_err = 0;

	if(zofa==NULL)
		setZofa();

	// Show B.C. at each boundary on screen
	if(first && mpi_rank==0) {
		for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
			if (zone->getKind() == FA_ZONE_BOUNDARY) {
				Param *param;

				if (getParam(param, zone->getName())) {
					if (param->getString() == "HOOK") {
						cout << "Applying HOOK                to zone: "<< zone->getName() << endl;
					} else if (param->getString() == "CBC") {
						double u_bc[3], T_bc, p_bc;

						for (int i=0; i<3; i++)
							u_bc[i] = param->getDouble(i+2);
						T_bc = param->getDouble(5);
						p_bc = param->getDouble(6);

						cout << "Applying CBC                 to zone: "<< zone->getName() <<"\t u_bc: "<< u_bc[0]<< " "<< u_bc[1]<< " "
								<< u_bc[2] <<" T_bc: " << T_bc << " p_bc: " << p_bc << endl;
					} else if (param->getString() == "CBC_SUBSONIC_INLET") {
						double angleU[3], Ttot, htot, ptot;

						for (int i=0; i<3; i++)
							angleU[i] = param->getDouble(i+2);
						Ttot = param->getDouble(5);
						ptot = param->getDouble(6);

						cout << "Applying CBC_SUBSONIC        to zone: "<< zone->getName() <<"\t angleU: "<<angleU[0]<<" "<<angleU[1]<<" "<<angleU[2]
						                                                                                                                           << " Ttot: " << Ttot << " Ptot: "<< ptot << endl;
					} else if (param->getString() == "CBC_SUBSONIC_OUTLET") {
						double p_bc = param->getDouble(2);

						cout << "Applying CBC_SUBSONIC_OUTLET to zone: "<< zone->getName() << "\t pOut: "<< p_bc << endl;
					} else if (param->getString() == "SYMMETRY") {
						cout << "Applying SYMMETRY            to zone: "<< zone->getName() << endl;
					} else if (param->getString() == "NEUMANN") {
						cout << "Applying NEUMANN             to zone: "<< zone->getName() << endl;
					} else if (param->getString() == "WALL") {
						int i=0;
						double T_bc = 0.0;
						if ((i = param->findString("TEMP")) != 0) {
							T_bc = param->getDouble(i+1);
						}

						if (T_bc > 0.0)    cout << "Applying WALL isothermal     to zone: "<< zone->getName() << "\t Temp_BC: " << T_bc << endl;
						else               cout << "Applying WALL adiabatic      to zone: "<< zone->getName() << endl;
					}
				} else {
					cerr << "Error: no bc set for: "<< zone->getName() << endl;
					bc_err = 1;
				}
			}
		}
	}

	// Set up B.C. at each boundary
	for(size_t index=0; index<faBoundary.size(); ++index) {
		int ifa = faBoundary[index];

		FaZone* zone = zofa[ifa];
		Param* param;
		if(getParam(param, zone->getName())) {
			if (getBoundaryType(ifa, param) == HOOK) {
				boundaryHook1D_AD(ifa, temp, vel, press, zone);

				setScalarBC1D_AD(ifa);
				ComputeBCProperties1D_T_AD(ifa);
			} else if (getBoundaryType(ifa, param) == CBC) {
				double u_bc[3], T_bc, p_bc;

				for (int i=0; i<3; i++)
					u_bc[i] = param->getDouble(i+2);
				T_bc = param->getDouble(5);
				p_bc = param->getDouble(6);

				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];

				double nVec[3];
				double area = normVec3d(nVec, fa_normal[ifa]);

				if (vecDotVec3d(u_bc, nVec) > 0.0) {  // outlet
					REALQ velMagN = vecDotVec3d_AD(vel[icv0], nVec);
					REALQ mach = fabs(velMagN)/sos[icv0];

					temp[icv1] = temp[icv0];
					for (int i=0; i<3; i++)
						vel[icv1][i] = vel[icv0][i];

#ifdef USE_CONDASSIGN
					if (mach >= 1.0)  press[icv1] = press[icv0];
					else              press[icv1] = p_bc;
//CASE1
//					condassign(press[icv1], mach+MACHINE_EPS, press[icv0], p_bc);
#else
					if (mach >= 1.0)  press[icv1] = press[icv0];
					else              press[icv1] = p_bc;
#endif
				} else {     // inlet
					temp[icv1] = T_bc;
					for (int i=0; i<3; i++)
						vel[icv1][i] = u_bc[i];
					press[icv1] = p_bc;
				}

				setScalarBC1D_AD(ifa);
				ComputeBCProperties1D_T_AD(ifa);
			} else if (getBoundaryType(ifa, param) == CBC_SUBSONIC_INLET) {
				double angleU[3], Ttot, htot, ptot;

				for (int i=0; i<3; i++)
					angleU[i] = param->getDouble(i+2);
				Ttot = param->getDouble(5);
				ptot = param->getDouble(6);

				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];
				//assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

				{
					REALQ u[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
					REALQ wPow2 = vecDotVec3d_AD(u, u);            // velocity squared
					REALQ velMag = sqrt(wPow2);                   // apply angle to extrapolated velocity
					for (int i=0; i<3; i++)
						vel[icv1][i] = angleU[i]*velMag;
					temp[icv1] = Ttot;                              // save temporary total temperature (to compute then total enthalpy)
				}

				setScalarBC1D_AD(ifa);
				ComputeBCProperties1D_T_AD(ifa);                  // total enthalpy from total temperature ?

				{
					REALQ wPow2 = vecDotVec3d_AD(vel[icv1], vel[icv1]);         // velocity squared
					enthalpy[icv1] -= 0.5* wPow2;                               // static enthalpy
				}

				ComputeBCProperties1D_H_AD(ifa);                  // static temperature and thermo properties from static enthalpy ?

				// Assumes isentropic relations to determine static pressure (= constant cp)
				// At first approximation ok, but could be improved; should for now be considered in defining p_bc
				press[icv1] = ptot*pow(temp[icv1]/Ttot, gamma[icv1]/(gamma[icv1]-1.0));
			} else if (getBoundaryType(ifa, param) == CBC_SUBSONIC_OUTLET) {
				double p_bc = param->getDouble(2);

				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];
				//assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

				// Assumes that the temperature is constant across outlet boundary (extrapolation) while p_bc is prescribed
				press[icv1] = p_bc;
				for (int i=0; i<3; i++)
					vel[icv1][i] = vel[icv0][i];

				temp[icv1] = temp[icv0];

				setScalarBC1D_AD(ifa);
				ComputeBCProperties1D_T_AD(ifa);
			} else if (getBoundaryType(ifa, param) == SYMMETRY) {
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];
				//assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

				REALX nVec[3];
				REALX area = normVec3d(nVec, fa_normal[ifa]);

				// flip u, APPROXIMATION ---> take velocity at the cell center
				REALQ u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
				REALQ un = vecDotVec3d_AD(nVec, u0);
				for (int i = 0; i < 3; i++)
					vel[icv1][i] = u0[i] - 1.0*un*nVec[i];
				//assert(fabs(vecDotVec3d(vel_bfa[ifa], nVec)) < 1.0e-10);

				temp[icv1]  = temp[icv0];
				press[icv1] = press[icv0];

				setScalarBC1D_AD(ifa);
				ComputeBCProperties1D_T_AD(ifa);
			} else if (getBoundaryType(ifa, param) == NEUMANN) {
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];
				//assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

				for (int i = 0; i < 3; i++)
					vel[icv1][i] = vel[icv0][i];

				temp[icv1] = temp[icv0];
				press[icv1] = press[icv0];

				setScalarBC1D_AD(ifa);
				ComputeBCProperties1D_T_AD(ifa);
			} else if (getBoundaryType(ifa, param) == WALL) {
				int i=0;
				double T_bc = 0.0;
				if ((i = param->findString("TEMP")) != 0)
					T_bc = param->getDouble(i+1);

				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];
				//assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

				for (int i = 0; i < 3; i++)
					vel[icv1][i] = 0.0;

				if (T_bc > 0.0)   temp[icv1] = T_bc;           // wall temperature
				else              temp[icv1] = temp[icv0];     // adiabatic wall

				press[icv1] = press[icv0];                      // Assumes zero pressure gradient at the wall

				setScalarBC1D_AD(ifa);
				ComputeBCProperties1D_T_AD(ifa);
			} else {
				if (mpi_rank == 0)
					cerr << "Error: cannot recognize bc for zone: " << zone->getName() << " bc: " << param->getString() << endl;
				bc_err = 1;
			}
		}

		// update density at boundary using EOS
		if((ifa>= 0 && ifa<nfa_b) || (ifa>=nfa && ifa<nfa_b2)) {
			int icv1 = cvofa[ifa][1];
			rho[icv1] = press[icv1]/(RoM[icv1]*temp[icv1]);
			sos[icv1] = sqrt(gamma[icv1]*press[icv1]/RoM[icv1]);
		}
	}
	if (bc_err != 0)
		throw(IKEWITHMODELS_ERROR_CODE);
	first = 0;
}
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
int IkeWithModels_AD::calcResidual1D_AD(const int icvCenter,
		REALA &rhs_rho_AD, REALA rhs_rhou_AD[3],REALA &rhs_rhoE_AD, REALAS *rhs_rhoScal_AD,
		REALQ *rho_AD, REALQ (*rhou_AD)[3], REALQ *rhoE_AD, double (*A)[5][5], double ***AScal, int flagImplicit) {
	//--------------
	// find faces
	//--------------
	vector<int> faInternal;
	vector<int> faBoundary;
	find1LayerFaceIndices(icvCenter, faInternal, faBoundary);

	vector<int> fa2LayerInternal;
	vector<int> fa2LayerBoundary;
	find2LayerFaceIndices(icvCenter, fa2LayerInternal, fa2LayerBoundary);

	//--------------
	// calculate rhs
	//--------------
	calcStateVariables1D_AD(icvCenter, rho_AD, rhou_AD, rhoE_AD);
	calcMaterialProperties1D_AD(icvCenter, faInternal, faBoundary, rho_AD, rhou_AD, rhoE_AD);
//	setBC1D_AD(icvCenter, fa2LayerBoundary, fa2LayerBoundary, rho_AD, rhou_AD, rhoE_AD);
	setBC1D_AD(icvCenter, fa2LayerInternal, fa2LayerBoundary, rho_AD, rhou_AD, rhoE_AD);
	if(mu_ref>0.0 || sndOrder==true)
		calcCv2Grad1D_AD(icvCenter, grad_u, vel, limiterNavierS, sos, epsilonSDWLS);
//	if(mu_ref>0.0) {
		FaZone** zoneArrayBoundary = NULL;
		if(!faBoundary.empty()) {
			zoneArrayBoundary = new FaZone* [faBoundary.size()];
			for(size_t index=0; index<faBoundary.size(); ++index) {
				int ifa = faBoundary[index];
				zoneArrayBoundary[index] = zofa[ifa];
			}
		}
		calcRansTurbViscMuet1D_AD(icvCenter, zoneArrayBoundary, faInternal, faBoundary, rho_AD, rhou_AD);
		if(zoneArrayBoundary != NULL) {
			delete [] zoneArrayBoundary;
			zoneArrayBoundary = NULL;
		}
//	}

	int CountReducedOrder = calcRhs1D_AD(icvCenter, rhs_rho_AD, rhs_rhou_AD, rhs_rhoE_AD, rhs_rhoScal_AD, rho_AD, rhou_AD, rhoE_AD, A, AScal, flagImplicit);

// IKJ
if(mpi_rank==0 && icvCenter<2) {
	int nScal = scalarTranspEqVector.size();
    ADscalar<adouble>* ZMean = NULL;
    ADscalar<adouble>* ZVar  = NULL;
    ADscalar<adouble>* CMean = NULL;

    for(int iScal=0; iScal<nScal; ++iScal) {
        if(iScal == getScalarTransportIndex("ZMean"))
    		ZMean = &(scalarTranspEqVector_AD[iScal].phi) ;
        if(iScal == getScalarTransportIndex("ZVar"))
		    ZVar  = &(scalarTranspEqVector_AD[iScal].phi) ;
        if(iScal == getScalarTransportIndex("CMean"))
		    CMean = &(scalarTranspEqVector_AD[iScal].phi) ;
    }
    assert(ZMean != NULL && ZVar != NULL && CMean != NULL);

    cout<<"IkeWithModels_AD::calcResidual1D_AD(): check for icv="<<icvCenter<<endl
        <<"  rho="<<rho[icvCenter]<<", rhou="<<rhou[icvCenter][0]<<", rhov="<<rhou[icvCenter][1]<<", rhoE="<<rhoE[icvCenter]<<", ZMean="<<(*ZMean)[icvCenter]<<", Zvar="<<(*ZVar)[icvCenter]<<", CMean="<<(*CMean)[icvCenter]<<endl
        <<"  press="<<press[icvCenter]<<", temp="<<temp[icvCenter]<<endl
        <<"  rhs_rhou_AD[1]="<<rhs_rhou_AD[1]<<endl;
}


	return CountReducedOrder;
}


/*
 * Method: setBC1D_AD
 * ------------------
 * Original code = setBC_AD in JoeWithModelsAD.cpp
 */
void IkeWithModels_AD::setBC1D_AD(const int icvCenter, vector<int>& faInternal, vector<int>& faBoundary, REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE) {
	static int first = 1;
	int bc_err = 0;

//	if(btofa==NULL)
//		setBtofa();
	if(zofa==NULL)
		setZofa();

	// Show B.C. at each boundary on screen
	if(first && mpi_rank==0) {
		for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
			if (zone->getKind() == FA_ZONE_BOUNDARY) {
				Param *param;

				if (getParam(param, zone->getName())) {
					if (param->getString() == "HOOK") {
						cout << "Applying HOOK                to zone: "<< zone->getName() << endl;
					} else if (param->getString() == "CBC") {
						double u_bc[3], T_bc, p_bc;

						for (int i=0; i<3; i++)
							u_bc[i] = param->getDouble(i+2);
						T_bc = param->getDouble(5);
						p_bc = param->getDouble(6);

						cout << "Applying CBC                 to zone: "<< zone->getName() <<"\t u_bc: "<< u_bc[0]<< " "<< u_bc[1]<< " "
								<< u_bc[2] <<" T_bc: " << T_bc << " p_bc: " << p_bc << endl;
					} else if (param->getString() == "CBC_SUBSONIC_INLET") {
						double angleU[3], Ttot, htot, ptot;

						for (int i=0; i<3; i++)
							angleU[i] = param->getDouble(i+2);
						Ttot = param->getDouble(5);
						ptot = param->getDouble(6);

						cout << "Applying CBC_SUBSONIC        to zone: "<< zone->getName() <<"\t angleU: "<<angleU[0]<<" "<<angleU[1]<<" "<<angleU[2]
						                                                                                                                           << " Ttot: " << Ttot << " Ptot: "<< ptot << endl;
					} else if (param->getString() == "CBC_SUBSONIC_OUTLET") {
						double p_bc = param->getDouble(2);

						cout << "Applying CBC_SUBSONIC_OUTLET to zone: "<< zone->getName() << "\t pOut: "<< p_bc << endl;
					} else if (param->getString() == "SYMMETRY") {
						cout << "Applying SYMMETRY            to zone: "<< zone->getName() << endl;
					} else if (param->getString() == "NEUMANN") {
						cout << "Applying NEUMANN             to zone: "<< zone->getName() << endl;
					} else if (param->getString() == "WALL") {
						int i=0;
						double T_bc = 0.0;
						if ((i = param->findString("TEMP")) != 0) {
							T_bc = param->getDouble(i+1);
						}

						if (T_bc > 0.0)    cout << "Applying WALL isothermal     to zone: "<< zone->getName() << "\t Temp_BC: " << T_bc << endl;
						else               cout << "Applying WALL adiabatic      to zone: "<< zone->getName() << endl;
					}
				} else {
					cerr << "Error: no bc set for: "<< zone->getName() << endl;
					bc_err = 1;
				}
			}
		}
	}

	// Set up B.C. at each boundary
	for(size_t index=0; index<faBoundary.size(); ++index) {
		int ifa = faBoundary[index];

		FaZone* zone = zofa[ifa];
		Param* param;
		if(getParam(param, zone->getName())) {
			if (getBoundaryType(ifa, param) == HOOK) {
				boundaryHook1D_AD(ifa, temp, vel, press, zone);

				setScalarBC1D_AD(ifa);
				ComputeBCProperties1D_T_AD(ifa);
			} else if (getBoundaryType(ifa, param) == CBC) {
				double u_bc[3], T_bc, p_bc;

				for (int i=0; i<3; i++)
					u_bc[i] = param->getDouble(i+2);
				T_bc = param->getDouble(5);
				p_bc = param->getDouble(6);

				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];

				double nVec[3];
				double area = normVec3d(nVec, fa_normal[ifa]);

				if (vecDotVec3d(u_bc, nVec) > 0.0) {  // outlet
					REALQ velMagN = vecDotVec3d_AD(vel[icv0], nVec);
					REALQ mach = fabs(velMagN)/sos[icv0];

					temp[icv1] = temp[icv0];
					for (int i=0; i<3; i++)
						vel[icv1][i] = vel[icv0][i];

					if (mach >= 1.0)  press[icv1] = press[icv0];
					else              press[icv1] = p_bc;
				} else {     // inlet
					temp[icv1] = T_bc;
					for (int i=0; i<3; i++)
						vel[icv1][i] = u_bc[i];
					press[icv1] = p_bc;
				}

				setScalarBC1D_AD(ifa);
				ComputeBCProperties1D_T_AD(ifa);
			} else if (getBoundaryType(ifa, param) == CBC_SUBSONIC_INLET) {
				double angleU[3], Ttot, htot, ptot;

				for (int i=0; i<3; i++)
					angleU[i] = param->getDouble(i+2);
				Ttot = param->getDouble(5);
				ptot = param->getDouble(6);

				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];
				//assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

				{
					REALQ u[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
					REALQ wPow2 = vecDotVec3d_AD(u, u);            // velocity squared
					REALQ velMag = sqrt(wPow2);                   // apply angle to extrapolated velocity
					for (int i=0; i<3; i++)
						vel[icv1][i] = angleU[i]*velMag;
					temp[icv1] = Ttot;                              // save temporary total temperature (to compute then total enthalpy)
				}

				setScalarBC1D_AD(ifa);
				ComputeBCProperties1D_T_AD(ifa);                  // total enthalpy from total temperature ?

				{
					REALQ wPow2 = vecDotVec3d_AD(vel[icv1], vel[icv1]);         // velocity squared
					enthalpy[icv1] -= 0.5* wPow2;                               // static enthalpy
				}

				ComputeBCProperties1D_H_AD(ifa);                  // static temperature and thermo properties from static enthalpy ?

				// Assumes isentropic relations to determine static pressure (= constant cp)
				// At first approximation ok, but could be improved; should for now be considered in defining p_bc
				press[icv1] = ptot*pow(temp[icv1]/Ttot, gamma[icv1]/(gamma[icv1]-1.0));
			} else if (getBoundaryType(ifa, param) == CBC_SUBSONIC_OUTLET) {
				double p_bc = param->getDouble(2);

				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];
				//assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

				// Assumes that the temperature is constant across outlet boundary (extrapolation) while p_bc is prescribed
				press[icv1] = p_bc;
				for (int i=0; i<3; i++)
					vel[icv1][i] = vel[icv0][i];

				temp[icv1] = temp[icv0];

				setScalarBC1D_AD(ifa);
				ComputeBCProperties1D_T_AD(ifa);
			} else if (getBoundaryType(ifa, param) == SYMMETRY) {
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];
				//assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

				REALX nVec[3];
				REALX area = normVec3d(nVec, fa_normal[ifa]);

				// flip u, APPROXIMATION ---> take velocity at the cell center
				REALQ u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
				REALQ un = vecDotVec3d_AD(nVec, u0);
				for (int i = 0; i < 3; i++)
					vel[icv1][i] = u0[i] - 1.0*un*nVec[i];
				//assert(fabs(vecDotVec3d(vel_bfa[ifa], nVec)) < 1.0e-10);

				temp[icv1]  = temp[icv0];
				press[icv1] = press[icv0];

				setScalarBC1D_AD(ifa);
				ComputeBCProperties1D_T_AD(ifa);
			} else if (getBoundaryType(ifa, param) == NEUMANN) {
//				if ((first)&&(mpi_rank == 0))
//					cout << "Applying NEUMANN             to zone: "<< zone->getName() << endl;

				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];
				//assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

				for (int i = 0; i < 3; i++)
					vel[icv1][i] = vel[icv0][i];

				temp[icv1] = temp[icv0];
				press[icv1] = press[icv0];

				setScalarBC1D_AD(ifa);
				ComputeBCProperties1D_T_AD(ifa);
			} else if (getBoundaryType(ifa, param) == WALL) {
				int i=0;
				double T_bc = 0.0;
				if ((i = param->findString("TEMP")) != 0)
					T_bc = param->getDouble(i+1);

				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];
				//assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

				for (int i = 0; i < 3; i++)
					vel[icv1][i] = 0.0;

				if (T_bc > 0.0)   temp[icv1] = T_bc;           // wall temperature
				else              temp[icv1] = temp[icv0];     // adiabatic wall

				press[icv1] = press[icv0];                      // Assumes zero pressure gradient at the wall

				setScalarBC1D_AD(ifa);
				ComputeBCProperties1D_T_AD(ifa);
			} else {
				if (mpi_rank == 0)
					cerr << "Error: cannot recognize bc for zone: " << zone->getName() << " bc: " << param->getString() << endl;
				bc_err = 1;
			}
		}

		// update density at boundary using EOS
		if((ifa>= 0 && ifa<nfa_b) || (ifa>=nfa && ifa<nfa_b2)) {
			int icv1 = cvofa[ifa][1];
			rho[icv1] = press[icv1]/(RoM[icv1]*temp[icv1]);
			sos[icv1] = sqrt(gamma[icv1]*press[icv1]/RoM[icv1]);
		}
	}
	if (bc_err != 0)
		throw(IKEWITHMODELS_ERROR_CODE);
	first = 0;
}

/*
 * Method: setBC1D
 * ---------------
 * Original code = setBC1D_AD in IkeWithModelsAD.cpp and setBC() in JoeWithModels.cpp
 * In some situations, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
 */
void IkeWithModels_AD::setBC1D(const int ifa, double *rho, double (*rhou)[3], double *rhoE) {
	int bc_err = 0;

	if(zofa==NULL)
		setZofa();

	FaZone* zone = zofa[ifa];
	Param* param;

	double *temp = UgpWithCvCompFlow::temp;
	double (*velUgp)[3] =  UgpWithCvCompFlow::vel;  // Note: if the name is just "vel", the compiler is confused with the adouble version vel
	double *press = UgpWithCvCompFlow::press;
	double *sos = UgpWithCvCompFlow::sos;
	double *enthalpy = UgpWithCvCompFlow::enthalpy;
	double *gamma = UgpWithCvCompFlow::gamma;
	double *RoM = UgpWithCvCompFlow::RoM;

	if(getParam(param, zone->getName())) {
		if (getBoundaryType(ifa, param) == HOOK) {
			boundaryHook1D(ifa, temp, velUgp, press, zone);

			setScalarBC1D(ifa);
			ComputeBCProperties1D_T(ifa);
		} else if (getBoundaryType(ifa, param) == CBC) {
			double u_bc[3], T_bc, p_bc;

			for (int i=0; i<3; i++)
				u_bc[i] = param->getDouble(i+2);
			T_bc = param->getDouble(5);
			p_bc = param->getDouble(6);

			int icv0 = cvofa[ifa][0];
			int icv1 = cvofa[ifa][1];

			double nVec[3];
			double area = normVec3d(nVec, fa_normal[ifa]);

			if (vecDotVec3d(u_bc, nVec) > 0.0) {  // outlet
				double velMagN = vecDotVec3d(velUgp[icv0], nVec);
				double mach = fabs(velMagN)/sos[icv0];

				temp[icv1] = temp[icv0];
				for (int i=0; i<3; i++)
					velUgp[icv1][i] = velUgp[icv0][i];

				if (mach >= 1.0)  press[icv1] = press[icv0];
				else              press[icv1] = p_bc;
			} else {     // inlet
				temp[icv1] = T_bc;
				for (int i=0; i<3; i++)
					velUgp[icv1][i] = u_bc[i];
				press[icv1] = p_bc;
			}

			setScalarBC1D(ifa);
			ComputeBCProperties1D_T(ifa);
		} else if (getBoundaryType(ifa, param) == CBC_SUBSONIC_INLET) {
			double angleU[3], Ttot, htot, ptot;

			for (int i=0; i<3; i++)
				angleU[i] = param->getDouble(i+2);
			Ttot = param->getDouble(5);
			ptot = param->getDouble(6);

			int icv0 = cvofa[ifa][0];
			int icv1 = cvofa[ifa][1];
			//assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

			{
				double u[3] = {velUgp[icv0][0], velUgp[icv0][1], velUgp[icv0][2]};
				double wPow2 = vecDotVec3d_AD(u, u);            // velocity squared
				double velMag = sqrt(wPow2);                   // apply angle to extrapolated velocity
				for (int i=0; i<3; i++)
					velUgp[icv1][i] = angleU[i]*velMag;
				temp[icv1] = Ttot;                              // save temporary total temperature (to compute then total enthalpy)
			}

			setScalarBC1D(ifa);
			ComputeBCProperties1D_T(ifa);                  // total enthalpy from total temperature ?

			{
				double wPow2 = vecDotVec3d(velUgp[icv1], velUgp[icv1]);         // velocity squared
				enthalpy[icv1] -= 0.5* wPow2;                               // static enthalpy
			}

			ComputeBCProperties1D_H(ifa);                  // static temperature and thermo properties from static enthalpy ?

			// Assumes isentropic relations to determine static pressure (= constant cp)
			// At first approximation ok, but could be improved; should for now be considered in defining p_bc
			press[icv1] = ptot*pow(temp[icv1]/Ttot, gamma[icv1]/(gamma[icv1]-1.0));
		} else if (getBoundaryType(ifa, param) == CBC_SUBSONIC_OUTLET) {
			double p_bc = param->getDouble(2);

			int icv0 = cvofa[ifa][0];
			int icv1 = cvofa[ifa][1];
			//assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

			// Assumes that the temperature is constant across outlet boundary (extrapolation) while p_bc is prescribed
			press[icv1] = p_bc;
			for (int i=0; i<3; i++)
				velUgp[icv1][i] = velUgp[icv0][i];

			temp[icv1] = temp[icv0];

			setScalarBC1D(ifa);
			ComputeBCProperties1D_T(ifa);
		} else if (getBoundaryType(ifa, param) == SYMMETRY) {
			int icv0 = cvofa[ifa][0];
			int icv1 = cvofa[ifa][1];
			//assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

			double nVec[3];
			double area = normVec3d(nVec, fa_normal[ifa]);

			// flip u, APPROXIMATION ---> take velocity at the cell center
			double u0[3] = {velUgp[icv0][0], velUgp[icv0][1], velUgp[icv0][2]};
			double un = vecDotVec3d_AD(nVec, u0);
			for (int i = 0; i < 3; i++)
				velUgp[icv1][i] = u0[i] - 1.0*un*nVec[i];
			//assert(fabs(vecDotVec3d(vel_bfa[ifa], nVec)) < 1.0e-10);

			temp[icv1]  = temp[icv0];
			press[icv1] = press[icv0];

			setScalarBC1D(ifa);
			ComputeBCProperties1D_T(ifa);
		} else if (getBoundaryType(ifa, param) == NEUMANN) {
			//				if ((first)&&(mpi_rank == 0))
			//					cout << "Applying NEUMANN             to zone: "<< zone->getName() << endl;

			int icv0 = cvofa[ifa][0];
			int icv1 = cvofa[ifa][1];
			//assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

			for (int i = 0; i < 3; i++)
				velUgp[icv1][i] = velUgp[icv0][i];

			temp[icv1] = temp[icv0];
			press[icv1] = press[icv0];

			setScalarBC1D(ifa);
			ComputeBCProperties1D_T(ifa);
		} else if (getBoundaryType(ifa, param) == WALL) {
			int i=0;
			double T_bc = 0.0;
			if ((i = param->findString("TEMP")) != 0)
				T_bc = param->getDouble(i+1);

			int icv0 = cvofa[ifa][0];
			int icv1 = cvofa[ifa][1];
			//assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

			for (int i = 0; i < 3; i++)
				velUgp[icv1][i] = 0.0;

			if (T_bc > 0.0)   temp[icv1] = T_bc;           // wall temperature
			else              temp[icv1] = temp[icv0];     // adiabatic wall

			press[icv1] = press[icv0];                      // Assumes zero pressure gradient at the wall

			setScalarBC1D(ifa);
			ComputeBCProperties1D_T(ifa);
		} else {
			if (mpi_rank == 0)
				cerr << "Error: cannot recognize bc for zone: " << zone->getName() << " bc: " << param->getString() << endl;
			bc_err = 1;
		}
	}

	// update density at boundary using EOS
	if((ifa>= 0 && ifa<nfa_b) || (ifa>=nfa && ifa<nfa_b2)) {
		int icv1 = cvofa[ifa][1];
		rho[icv1] = press[icv1]/(RoM[icv1]*temp[icv1]);
		sos[icv1] = sqrt(gamma[icv1]*press[icv1]/RoM[icv1]);
	}

	if (bc_err != 0)
		throw(IKEWITHMODELS_ERROR_CODE);
}

/*
 * Method: setScalarBC1D_AD
 * ------------------------
 * Original code = setScalarBC_AD() in ScalarsAD.cpp
 */
void IkeWithModels_AD::setScalarBC1D_AD(const int ifa) {
	int nScal = scalarTranspEqVector.size();

	for (int scal = 0; scal < nScal; scal++) {
		string scalName = scalarTranspEqVector[scal].getName();
		double phiBCval = 0.0, phiBCflux = 0.0;

#ifdef USE_MEM_SAVING_ADVAR // note: USE_MEM_SAVING_ADVAR is defined in UgpWithCvCompFlow.h
		ADscalar<adouble> *phi = &(scalarTranspEqVector_AD[scal].phi);
#else
		REALQ *phi = scalarTranspEqVector_AD[scal].phi;
#endif

		FaZone* zone = zofa[ifa];
		if (zone->getKind() == FA_ZONE_BOUNDARY) {
			// .............................................................................................
			// HOOK BOUNDARY CONDITION
			// .............................................................................................
			if (scalarZoneIsHook(zone->getName(), scalName)) {
				// user defined boundary hook
#ifdef USE_MEM_SAVING_ADVAR
				boundaryHookScalarRansTurb1D_AD(ifa, *phi, zone, scalName);
				boundaryHookScalarRansComb1D_AD(ifa, *phi, zone, scalName);
#else
				boundaryHookScalarRansTurb1D_AD(ifa, phi, zone, scalName);
				boundaryHookScalarRansComb1D_AD(ifa, phi, zone, scalName);
#endif
			}
			// .............................................................................................
			// DIRICHLET BOUNDARY CONDITION
			// .............................................................................................
			else if (scalarZoneIsDirichlet(phiBCval, zone->getName(), scalName)) {
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];
#ifdef USE_MEM_SAVING_ADVAR
	#ifdef USE_CONDASSIGN
				condassign((*phi)[icv1], vecDotVec3d_AD(vel[icv1], fa_normal[ifa])-1.0e-8, (*phi)[icv0], phiBCval);
	#else
				if (vecDotVec3d_AD(vel[icv1], fa_normal[ifa]) > 1.0e-8)
					(*phi)[icv1] = (*phi)[icv0];
				else
					(*phi)[icv1] = phiBCval;
	#endif
#else
				if (vecDotVec3d_AD(vel[icv1], fa_normal[ifa]) > 1.0e-8)
					phi[icv1] = phi[icv0];
				else
					phi[icv1] = phiBCval;
#endif
			}
			// .............................................................................................
			// FLUX BOUNDARY CONDITION: NOT IMPLEMENTED
			// .............................................................................................
			else if (scalarZoneIsFlux(phiBCflux, zone->getName(), scalName)) {
				cerr << "scalarZoneIsFlux not implemented yet!" << endl;
				throw(-1);
			}
			// .............................................................................................
			// OTHER BOUNDARY CONDITIONS: NO FLUX
			// .............................................................................................
			else { // OTHER BOUNDARY CONDITIONS: NO FLUX
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];

#ifdef USE_MEM_SAVING_ADVAR
				(*phi)[icv1] = (*phi)[icv0];
#else
				phi[icv1] = phi[icv0];
#endif
			}
		}
	}
}

/*
 * Method: setScalarBC1D
 * ------------------------
 * Original code = setScalarBC1D_AD() in IkeWithModels.cpp and setBC() in Scalars.cpp
 * In some situation, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
 */
void IkeWithModels_AD::setScalarBC1D(const int ifa) {
	int nScal = scalarTranspEqVector.size();

	double (*velUgp)[3] =  UgpWithCvCompFlow::vel; // Note: if the name is just "vel", the compiler is confused with the adouble version vel

	for (int scal = 0; scal < nScal; scal++) {
		string scalName = scalarTranspEqVector[scal].getName();
		double phiBCval = 0.0, phiBCflux = 0.0;

		double *phi = scalarTranspEqVector[scal].phi;

		FaZone* zone = zofa[ifa];
		if (zone->getKind() == FA_ZONE_BOUNDARY) {
			// .............................................................................................
			// HOOK BOUNDARY CONDITION
			// .............................................................................................
			if (scalarZoneIsHook(zone->getName(), scalName)) {
				// user defined boundary hook
				boundaryHookScalarRansTurb1D(ifa, phi, zone, scalName);
				boundaryHookScalarRansComb1D(ifa, phi, zone, scalName);
			}
			// .............................................................................................
			// DIRICHLET BOUNDARY CONDITION
			// .............................................................................................
			else if (scalarZoneIsDirichlet(phiBCval, zone->getName(), scalName)) {
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];
				if (vecDotVec3d(velUgp[icv1], fa_normal[ifa]) > 1.0e-8)
					phi[icv1] = phi[icv0];
				else
					phi[icv1] = phiBCval;
			}
			// .............................................................................................
			// FLUX BOUNDARY CONDITION: NOT IMPLEMENTED
			// .............................................................................................
			else if (scalarZoneIsFlux(phiBCflux, zone->getName(), scalName)) {
				cerr << "scalarZoneIsFlux not implemented yet!" << endl;
				throw(-1);
			}
			// .............................................................................................
			// OTHER BOUNDARY CONDITIONS: NO FLUX
			// .............................................................................................
			else { // OTHER BOUNDARY CONDITIONS: NO FLUX
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];

				phi[icv1] = phi[icv0];
			}
		}
	}
}

/*
 * Method: ComputeBCProperties1D_T_AD
 * ----------------------------------
 * brief Compute for a given temperature the properties of the mixture at a face.
 * Original code: ComputeBCProperties_T_AD() in UgpWithCvCompFlowAD.h
 */
void IkeWithModels_AD::ComputeBCProperties1D_T_AD(const int ifa) {
	int icv1=cvofa[ifa][1];
	gamma[icv1] = GAMMA;
	RoM[icv1] = R_gas;
	enthalpy[icv1] = GAMMA * R_gas / (GAMMA - 1.0) * temp[icv1];

	if (mu_ref > 0.0) {
		if (viscMode == "SUTHERLAND") {
			int icv1=cvofa[ifa][1];
			mul_fa[ifa] = mu_ref*pow(temp[icv1]/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(temp[icv1] + SL_Sref);
		} else if (viscMode == "POWERLAW") {
			int icv1=cvofa[ifa][1];
			mul_fa[ifa] = mu_ref*pow(temp[icv1]/T_ref, mu_power_law);
		} else {
			cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
			throw(-1);
		}

		lamOcp_fa[ifa] = mul_fa[ifa] / Pr;
	}
}

/*
 * Method: ComputeBCProperties1D_T
 * -------------------------------
 * brief Compute for a given temperature the properties of the mixture at a face.
 * Original code: ComputeBCProperties1D_T_AD() in IkeWithModel.cpp and ComputeBCProperties_T() UgpWithCvCompFlow.h
 * In some situation, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
 */
void IkeWithModels_AD::ComputeBCProperties1D_T(const int ifa) {
	double *gamma = UgpWithCvCompFlow::gamma;
	double *RoM = UgpWithCvCompFlow::RoM;
	double *enthalpy = UgpWithCvCompFlow::enthalpy;
	double *temp = UgpWithCvCompFlow::temp;
	double *mul_fa = UgpWithCvCompFlow::mul_fa;
	double *lamOcp_fa = UgpWithCvCompFlow::lamOcp_fa;

	int icv1=cvofa[ifa][1];
	gamma[icv1] = GAMMA;
	RoM[icv1] = R_gas;
	enthalpy[icv1] = GAMMA * R_gas / (GAMMA - 1.0) * temp[icv1];

	if (mu_ref > 0.0) {
		if (viscMode == "SUTHERLAND") {
			int icv1=cvofa[ifa][1];
			mul_fa[ifa] = mu_ref*pow(temp[icv1]/SL_Tref, 1.5)*(SL_Tref + SL_Sref)/(temp[icv1] + SL_Sref);
		} else if (viscMode == "POWERLAW") {
			int icv1=cvofa[ifa][1];
			mul_fa[ifa] = mu_ref*pow(temp[icv1]/T_ref, mu_power_law);
		} else {
			cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
			throw(-1);
		}

		lamOcp_fa[ifa] = mul_fa[ifa] / Pr;
	}
}

/*
 * Method: ComputeBCProperties1D_H_AD
 * ----------------------------------
 * brief Compute for a given enthalpy the properties of the mixture at a face.
 * Original code: ComputeBCProperties1D_H_AD() in IkeWithModels.cpp and ComputeBCProperties1D_H() in UgpWithCvCompFlow.h
 * In some situation, setBC() in JoeWithModels.cpp cannot update ggf faces. Thus, the user needs to manually update them.
 */
void IkeWithModels_AD::ComputeBCProperties1D_H_AD(const int ifa) {
	int icv1 = cvofa[ifa][1];
	gamma[icv1] = GAMMA;
	RoM[icv1] = R_gas;
	temp[icv1] = enthalpy[icv1]*(GAMMA-1.0)/(GAMMA*R_gas);

	if (mu_ref > 0.0) {
		if (viscMode == "SUTHERLAND") {
			int icv1 = cvofa[ifa][1];
			mul_fa[ifa] = mu_ref*pow(temp[icv1]/SL_Tref, 1.5)*(SL_Tref+SL_Sref)/(temp[icv1]+SL_Sref);
		} else if (viscMode == "POWERLAW") {
			int icv1 = cvofa[ifa][1];
			mul_fa[ifa] = mu_ref*pow(temp[icv1]/T_ref, mu_power_law);
		} else {
			cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
			throw(-1);
		}

		lamOcp_fa[ifa] = mul_fa[ifa] / Pr;
	}
}

/*
 * Method: ComputeBCProperties1D_H
 * -------------------------------
 * brief Compute for a given enthalpy the properties of the mixture at a face.
 * Original code: ComputeBCProperties1D_H_AD() in UgpWithCvCompFlowAD.h
 */
void IkeWithModels_AD::ComputeBCProperties1D_H(const int ifa) {
	double *gamma = UgpWithCvCompFlow::gamma;
	double *RoM = UgpWithCvCompFlow::RoM;
	double *enthalpy = UgpWithCvCompFlow::enthalpy;
	double *temp = UgpWithCvCompFlow::temp;
	double *mul_fa = UgpWithCvCompFlow::mul_fa;
	double *lamOcp_fa = UgpWithCvCompFlow::lamOcp_fa;

	int icv1 = cvofa[ifa][1];
	gamma[icv1] = GAMMA;
	RoM[icv1] = R_gas;
	temp[icv1] = enthalpy[icv1]*(GAMMA-1.0)/(GAMMA*R_gas);

	if (mu_ref > 0.0) {
		if (viscMode == "SUTHERLAND") {
			int icv1 = cvofa[ifa][1];
			mul_fa[ifa] = mu_ref*pow(temp[icv1]/SL_Tref, 1.5)*(SL_Tref+SL_Sref)/(temp[icv1]+SL_Sref);
		} else if (viscMode == "POWERLAW") {
			int icv1 = cvofa[ifa][1];
			mul_fa[ifa] = mu_ref*pow(temp[icv1]/T_ref, mu_power_law);
		} else {
			cerr << "viscosity mode not recognized, current options are \"MU_MODE = SUTHERLAND\" and \"MU_MODE = POWERLAW\"" << endl;
			throw(-1);
		}

		lamOcp_fa[ifa] = mul_fa[ifa] / Pr;
	}
}

#ifdef USE_MEM_SAVING_ADVAR_1D_
/*
 * Method: calcRhs1D_AD
 * --------------------
 * Calculate the Rhs of the N-S equations in the so-called "1D-style"
 * Original code = calcRhs_AD() in JoeWithModelsAD.cpp
 */
int IkeWithModels_AD::calcRhs1D_AD(int icvCenter, REALA &rhs_rho, REALA rhs_rhou[3], REALA &rhs_rhoE, REALAS *rhs_rhoScal,
		ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<REALQ> &rhoE,
		double (*A)[5][5], double ***AScal, int flagImplicit) {
	static bool firstCall = true;

	// set RHS to zero
	rhs_rho = 0.0;
	for (int i = 0; i < 3; i++)
		rhs_rhou[i] = 0.0;
	rhs_rhoE = 0.0;

	// set scalars RHS to zero
	for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
		for (int icv = 0; icv < ncv_gg; icv++)
			rhs_rhoScal[iScal] = 0.0;

	// =======================================================================================
	// NAVIER-STOKES
	// =======================================================================================
	// compute Euler Flux for NS and scalars
	int myCountReducedOrder = calcEulerFlux1D_AD(icvCenter, rhs_rho, rhs_rhou, rhs_rhoE, rhs_rhoScal, rho, rhou, rhoE, A, AScal, flagImplicit);

// IKJ
if(mpi_rank==0 && icvCenter==0) {
	int icv = icvCenter;

	cout<<"IkeWithModels_AD::calcRhs1D_AD(): rhs from calcEulerFlux1D_AD()"<<endl;
	cout<<"    rhs = "<<rhs_rho;
	for(int i=0; i<3; ++i) cout<<", "<<rhs_rhou[i];
	cout<<", "<<rhs_rhoE;
	for(int i=0; i<3; ++i) cout<<", "<<rhs_rhoScal[i];
	cout<<endl
		<<endl;
}

	// compute viscous Flux for NS
#ifdef USE_ARTIF_VISC_WITH_MEM_SAVING
	if(UgpWithCvCompFlow::turnOnArtifVisc) {
		calcArtifVisc1D_AD(icvCenter, artifViscMag_AD,
				artifVisc_bulkViscOnly, artifVisc_shockOnly,
				artifVisc_smoothstepThresh, artifVisc_type, artifVisc_coeff);
	}
	if (mu_ref > 0.0 || UgpWithCvCompFlow::turnOnArtifVisc)
		calcViscousFluxNS1D_AD(icvCenter, rhs_rho, rhs_rhou, rhs_rhoE, rho, rhou, rhoE, A, flagImplicit);
#else
	if (mu_ref > 0.0)
		calcViscousFluxNS1D_AD(icvCenter, rhs_rho, rhs_rhou, rhs_rhoE, rho, rhou, rhoE, A, flagImplicit);
#endif

	// add source terms to RHS of Navier-Stokes equations
	sourceHook1D_AD(icvCenter, rhs_rho, rhs_rhou, rhs_rhoE, A);
	sourceHookRansTurb1D_AD(icvCenter, rhs_rho, rhs_rhou, rhs_rhoE, A);
	sourceHookRansComb1D_AD(icvCenter, rhs_rho, rhs_rhou, rhs_rhoE, A);

	// =======================================================================================
	// SCALARS
	// =======================================================================================
	for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++) {
		// compute viscous Flux for scalars and add source terms to RHS of scalar equations
		if (!scalarTranspEqVector[iScal].diffTerm) {
			/* Do nothing */
		} else {
			/*** From here ... */
			string scalName(scalarTranspEqVector[iScal].getName());
			double phiBCvalDummy = 0.0, phiBCfluxDummy = 0.0;

			for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
				if (zone->getKind() == FA_ZONE_BOUNDARY) {
					if (scalarZoneIsHook(zone->getName(), scalName)) {
						if(firstCall && mpi_rank == 0)
							cout << scalName << ": " << zone->getName() << " is HOOK" << endl;
					} else if (scalarZoneIsDirichlet(phiBCvalDummy, zone->getName(), scalName)) {
						if(firstCall && mpi_rank == 0)
							cout << scalName << ": " << zone->getName() << " is Dirichlet: VALUE: " << phiBCvalDummy << endl;
					} else if (scalarZoneIsFlux(phiBCfluxDummy, zone->getName(), scalName)) {
						cerr << "scalarZoneIsFlux not implemented yet!" << endl;
						throw(-1);
					} else {
						if(firstCall && mpi_rank == 0)
							cout << scalName << ": " << zone->getName() << " is Symmetry" << endl;
					}
				}
			}
			/*** ... to here, this part comes from calcViscousFluxScalar_new_AD() in ScalarsAD.cpp */
			if (AScal == NULL) {
				if (mu_ref > 0.0)
					calcViscousFluxScalar_new1D_AD(icvCenter, iScal, rhs_rhoScal[iScal], NULL, scalarTranspEqVector[iScal], scalarTranspEqVector_AD[iScal], flagImplicit);
			} else {
				if (mu_ref > 0.0)
					calcViscousFluxScalar_new1D_AD(icvCenter, iScal, rhs_rhoScal[iScal], AScal[iScal][5], scalarTranspEqVector[iScal], scalarTranspEqVector_AD[iScal], flagImplicit);
			}
		}

		if (AScal == NULL) {
			sourceHookScalarRansTurb_new1D_AD(icvCenter, rhs_rhoScal[iScal], NULL, scalarTranspEqVector[iScal].getName(), flagImplicit);
			sourceHookScalarRansComb_new1D_AD(icvCenter, rhs_rhoScal[iScal], NULL, scalarTranspEqVector[iScal].getName(), flagImplicit);
		} else {
			sourceHookScalarRansTurb_new1D_AD(icvCenter, rhs_rhoScal[iScal], AScal[iScal][5], scalarTranspEqVector[iScal].getName(), flagImplicit);
			sourceHookScalarRansComb_new1D_AD(icvCenter, rhs_rhoScal[iScal], AScal[iScal][5], scalarTranspEqVector[iScal].getName(), flagImplicit);
		}
	}
	firstCall = false;

	return myCountReducedOrder;
}

/*
 * Method: calcEulerFlux1D_AD
 * --------------------------
 * Calculate Euler flux in the so-called "1D style"
 * Return: the number of times switched back to first order at faces
 * Original code = calcEulerFlux_AD() in JoeWithModelsAD.cpp
 */
int IkeWithModels_AD::calcEulerFlux1D_AD(const int icvCenter, REALA &rhs_rho, REALA rhs_rhou[3], REALA &rhs_rhoE, REALAS *rhs_rhoScal,
		ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<REALQ> &rhoE,
		double (*A)[5][5], double ***AScal, int flagImplicit) {
	static bool firstCall = true;

	int nScal = scalarTranspEqVector.size();

	double (*Apl)[5] = NULL;
	double (*Ami)[5] = NULL;
	if (flagImplicit) {
		Apl = new double[5][5];
		Ami = new double[5][5];

		for (int i = 0; i < 5; i++)
			for (int j = 0; j < 5; j++)
				Apl[i][j] = Ami[i][j] = 0.0;
	}

	// Implicit matrix for scalars: the definition is not easy, so ask an expert (e.g., R. Pecnik)
	double (*AplScal)[6] = NULL;
	double (*AmiScal)[6] = NULL;
	if (flagImplicit) {
		if (nScal > 0) AplScal = new double[nScal][6];
		if (nScal > 0) AmiScal = new double[nScal][6];

		for (int iScal = 0; iScal < nScal; iScal++)
			for (int i = 0; i <= 5; i++)
				AplScal[iScal][i] = AmiScal[iScal][i] = 0.0;
	}

	REALQ Frho, Frhou[3], FrhoE;
	REALQS *FrhoScal     = NULL;         if (nScal > 0) FrhoScal     = new REALQS[nScal];
	REALQS *Scalar0      = NULL;         if (nScal > 0) Scalar0      = new REALQS[nScal];            // cell face if second order
	REALQS *Scalar1      = NULL;         if (nScal > 0) Scalar1      = new REALQS[nScal];
	REALQS *ScalCV0      = NULL;         if (nScal > 0) ScalCV0      = new REALQS[nScal];            // cell face if second order
	REALQS *ScalCV1      = NULL;         if (nScal > 0) ScalCV1      = new REALQS[nScal];
	double *ScalConvTerm = NULL;         if (nScal > 0) ScalConvTerm = new double[nScal];            // 0 if convective term no considered, otherwise 1

	// count how many cells switched back to first order due to extrapolated negative density or pressure/temperature at the faces
	int myCountReducedOrder = 0;

	// save the index of kine if defined and save convTerm for speed
	int kine_Index = getScalarTransportIndex("kine");
	for (int iScal = 0; iScal < nScal; iScal++)
		ScalConvTerm[iScal] = double(scalarTranspEqVector[iScal].convTerm);

	// =============================================================================================
	// compute gradients, with boundary values
	// =============================================================================================
	if (sndOrder == true) {
		calcCv2Grad1D_AD(icvCenter, grad_rho, rho, limiterNavierS, rho, epsilonSDWLS);
#ifdef temp_reconstruction
		calcCv2Grad1D_AD(icvCenter, grad_temp, temp, limiterNavierS, temp, epsilonSDWLS);
#else
		calcCv2Grad1D_AD(icvCenter, grad_p, press, limiterNavierS, press, epsilonSDWLS);
#endif

		// Check if grad_rho & grad_p don't give NaN (just do this only for firstCall which also means icvCenter==0)
		if(firstCall) {
			int noc_f_center = nbocv_all_i[icvCenter];
			int noc_l_center = nbocv_all_i[icvCenter+1]-1;
			for (int noc_center = noc_f_center; noc_center <= noc_l_center; ++noc_center) {
				int icv = nbocv_all_v[noc_center]; // Neighboring CVs

				if(icv<ncv_g) {
					if(isNaN(grad_rho[icv][0].value()) || isNaN(grad_rho[icv][1].value()) || isNaN(grad_rho[icv][2].value())) {
						printf("Error in IkeWithModels_AD::calcEulerFlux1D_AD(): NaN detected in grad_rho[%d]=(%.3e,%.3e,%.3e) @x_cv=(%.2e,%.2e,%.2e) - mpi_rank=%d\n",
								icv, grad_rho[icv][0].value(),grad_rho[icv][1].value(),grad_rho[icv][2].value(), x_cv[icv][0],x_cv[icv][1],x_cv[icv][2], mpi_rank);
						printf("\n");
						freeMemForCalcEulerFlux1D_AD(Apl, Ami, AplScal, AmiScal, FrhoScal, Scalar0, Scalar1, ScalCV0, ScalCV1, ScalConvTerm, nScal);
						throw(IKEWITHMODELS_ERROR_CODE);
					}

					if(isNaN(grad_p[icv][0].value()) || isNaN(grad_p[icv][1].value()) || isNaN(grad_p[icv][2].value())) {
						printf("Error in IkeWithModels_AD::calcEulerFlux1D_AD(): NaN detected in grad_p[%d]=(%.3e,%.3e,%.3e) @x_cv=(%.2e,%.2e,%.2e) - mpi_rank=%d\n",
								icv, grad_p[icv][0].value(),grad_p[icv][1].value(),grad_p[icv][2].value(), x_cv[icv][0],x_cv[icv][1],x_cv[icv][2], mpi_rank);
						printf("\n");
						freeMemForCalcEulerFlux1D_AD(Apl, Ami, AplScal, AmiScal, FrhoScal, Scalar0, Scalar1, ScalCV0, ScalCV1, ScalConvTerm, nScal);
						throw(IKEWITHMODELS_ERROR_CODE);
					}
				}
			}
		}

		for (int iScal = 0; iScal < nScal; iScal++) {
			if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE") {
				// For the scalar reconstruction, Grad(rho*Phi) is required
				// Gradients of rho*Phi are saved in grad_phi temporarily
				// Gradients of rho*Phi are also limited like rho with alpha_rho
				// Boundary face values rho_bfa*Phi_fa are saved in Phi_fa temporarily
				ADscalar<REALQS> rhoPhi("rhoPhi");
				rhoPhi.allocate(scalarTranspEqVector_AD[iScal].phi.getMyNbocv2(), scalarTranspEqVector_AD[iScal].phi.getMyIcv());

				// Compute rho*Phi
				int noc_f_center = nbocv_all_i[icvCenter];
				int noc_l_center = nbocv_all_i[icvCenter+1]-1;
				for (int noc_center = noc_f_center; noc_center <= noc_l_center; ++noc_center) {
					int icv = nbocv_all_v[noc_center]; // Neighboring CVs
					rhoPhi[icv] = rho[icv] * (scalarTranspEqVector_AD[iScal].phi)[icv];
				}

				// Compute gradients of rho*Phi and limit based on rho*Phi
				calcCv2Grad1D_AD(icvCenter, scalarTranspEqVector_AD[iScal].grad_phi, rhoPhi, limiterNavierS, rhoPhi, epsilonSDWLS);

				rhoPhi.clear();
			} else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD") {
				calcCv2Grad1D_AD(icvCenter, scalarTranspEqVector_AD[iScal].grad_phi, scalarTranspEqVector_AD[iScal].phi, limiterNavierS, scalarTranspEqVector_AD[iScal].phi, epsilonSDWLS);
			} else {
				if(mpi_rank==0)
					cerr << "### IkeWithModels::calcEulerFlux1D_AD() => Wrong reconstruction type for scalars! ###" << endl;
				freeMemForCalcEulerFlux1D_AD(Apl, Ami, AplScal, AmiScal, FrhoScal, Scalar0, Scalar1, ScalCV0, ScalCV1, ScalConvTerm, nScal);
				throw(-1);
			}
		}
	}

	// Calculate flux
	int nofa_f = faocv_i[icvCenter];
	int nofa_l = faocv_i[icvCenter+1]-1;
	for (int foc = nofa_f; foc <= nofa_l; foc++) {
		int ifa = faocv_v[foc];
		// ===============================================================================================
		// cycle through internal faces, assembling flux to both sides
		// ===============================================================================================
		if((ifa >= nfa_b && ifa<nfa) || (ifa >= nfa_b2 && ifa < nfa_b2gg)) {
			int icv0 = cvofa[ifa][0];
			int icv1 = cvofa[ifa][1];
			assert( icv0 >= 0 && icv0 < ncv_gg);
			assert( icv1 >= 0 && icv1 < ncv_gg);

			int noc00, noc01, noc11, noc10;
			if (flagImplicit && icv0 < ncv && icv1 < ncv_g)
				getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

			// face unit normal and area...
			double nVec[3] = {0.0, 0.0, 0.0};
			double area = normVec3d(nVec, fa_normal[ifa]);

			// .............................................................................................
			// reconstruction of variables at faces: rho, u, T or P, scalars
			// .............................................................................................
			REALQ rho0 = rho[icv0];
			REALQ u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
			REALQ p0 = press[icv0];
			REALQ T0 = temp[icv0];
			REALQ h0 = enthalpy[icv0];
			REALQ gam0 = gamma[icv0];
			REALQ R0 = RoM[icv0];
			REALQS kineCV0 = 0.0;          // cell center
			REALQS kineFA0 = 0.0;          // cell face if second order

			REALQ rho1 = rho[icv1];
			REALQ u1[3] = {vel[icv1][0], vel[icv1][1], vel[icv1][2]};
			REALQ p1 = press[icv1];
			REALQ T1 = temp[icv1];
			REALQ h1 = enthalpy[icv1];
			REALQ gam1 = gamma[icv1];
			REALQ R1 = RoM[icv1];
			REALQS kineCV1 = 0.0;
			REALQS kineFA1 = 0.0;

			for (int iScal = 0; iScal < nScal; iScal++) {
				ScalCV0[iScal] = Scalar0[iScal] = scalarTranspEqVector_AD[iScal].phi[icv0];
				ScalCV1[iScal] = Scalar1[iScal] = scalarTranspEqVector_AD[iScal].phi[icv1];
			}

			if (sndOrder == true) {
				double r0[3] = {0.0, 0.0, 0.0}, r1[3] = {0.0, 0.0, 0.0};
				vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);
				vecMinVec3d(r1, x_fa[ifa], x_cv[icv1]);

				// ----------------------------------------
				// left side
				// ----------------------------------------
				rho0 += vecDotVec3d_AD(r0, grad_rho[icv0]);

#ifdef temp_reconstruction
				T0 += vecDotVec3d_AD(r0, grad_temp[icv0]);
	#ifdef USE_CONDASSIGN
				if ((T0 <= 0.0) || (rho0 <= 0.0)) {
					T0 = temp[icv0];
					rho0 = rho[icv0];
					myCountReducedOrder++;
				}
//CASE1
//				condassign(T0, -T0*rho0+MACHINE_EPS, temp[icv0]);
//				condassign(rho0, -T0*rho0+MACHINE_EPS, rho[icv0]);
//
//				if ((T0.value() <= 0.0) || (rho0.value() <= 0.0))
//					myCountReducedOrder++;
	#else
				if ((T0 <= 0.0) || (rho0 <= 0.0)) {
					T0 = temp[icv0];
					rho0 = rho[icv0];
					myCountReducedOrder++;
				}
	#endif
#else
				p0 += vecDotVec3d_AD(r0, grad_p[icv0]);
	#ifdef USE_CONDASSIGN
				if ((p0 <= 0.0) || (rho0 <= 0.0)) {
					p0 = press[icv0];
					rho0 = rho[icv0];
					myCountReducedOrder++;
				}
//CASE1
//				condassign(p0, -p0*rho0+MACHINE_EPS, press[icv0]);
//				condassign(rho0, -p0*rho0+MACHINE_EPS, rho[icv0]);
//				if ((p0.value() <= 0.0) || (rho0.value() <= 0.0))
//					myCountReducedOrder++;
	#else
				if ((p0 <= 0.0) || (rho0 <= 0.0)) {
					p0 = press[icv0];
					rho0 = rho[icv0];
					myCountReducedOrder++;
				}
	#endif
#endif
				else {
					for (int i = 0; i < 3; i++)
						u0[i] += vecDotVec3d_AD(r0, grad_u[icv0][i]);
					for (int iScal = 0; iScal < nScal; iScal++) {
						if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
							Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d_AD(r0, scalarTranspEqVector_AD[iScal].grad_phi[icv0])) / rho0;
						else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
							Scalar0[iScal] += vecDotVec3d_AD(r0, scalarTranspEqVector_AD[iScal].grad_phi[icv0]);
					}
				}

				// ----------------------------------------
				// right side
				// ----------------------------------------
				rho1 += vecDotVec3d_AD(r1, grad_rho[icv1]);
#ifdef temp_reconstruction
				T1 += vecDotVec3d_AD(r1, grad_temp[icv1]);
	#ifdef USE_CONDASSIGN
				if ((T1 <= 0.0) || (rho1 <= 0.0)) {
					T1 = temp[icv1];
					rho1 = rho[icv1];
					myCountReducedOrder++;
				}
//CASE1
//				condassign(T1, -T1*rho1+MACHINE_EPS, temp[icv1]);
//				condassign(rho1, -T1*rho1+MACHINE_EPS, rho[icv1]);
//				if ((T1.value() <= 0.0) || (rho1.value() <= 0.0))
//					myCountReducedOrder++;
	#else
				if ((T1 <= 0.0) || (rho1 <= 0.0)) {
					T1 = temp[icv1];
					rho1 = rho[icv1];
					myCountReducedOrder++;
				}
	#endif
#else
				p1 += vecDotVec3d_AD(r1, grad_p[icv1]);
	#ifdef USE_CONDASSIGN
				if ((p1 <= 0.0) || (rho1 <= 0.0)) {
					p1 = press[icv1];
					rho1 = rho[icv1];
					myCountReducedOrder++;
				}
//CASE1
//				condassign(p1, -p1*rho1+MACHINE_EPS, press[icv1]);
//				condassign(rho1, -p1*rho1+MACHINE_EPS, rho[icv1]);
//				if ((p1.value() <= 0.0) || (rho1.value() <= 0.0))
//					myCountReducedOrder++;
	#else
				if ((p1 <= 0.0) || (rho1 <= 0.0)) {
					p1 = press[icv1];
					rho1 = rho[icv1];
					myCountReducedOrder++;
				}
	#endif
#endif
				else {
					for (int i = 0; i < 3; i++)
						u1[i] += vecDotVec3d_AD(r1, grad_u[icv1][i]);
					for (int iScal = 0; iScal < nScal; iScal++) {
						if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
							Scalar1[iScal] = (rho[icv1] * Scalar1[iScal] + vecDotVec3d_AD(r1, scalarTranspEqVector_AD[iScal].grad_phi[icv1])) / rho1;    // CHANGE FOR SCALAR AD

						else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
							Scalar1[iScal] += vecDotVec3d_AD(r1, scalarTranspEqVector_AD[iScal].grad_phi[icv1]);
					}
				}

				// .............................................................................................
				// calculation of other variables at faces: p/T, h, R, gam
				// .............................................................................................
#ifdef temp_reconstruction
				calcThermoProp_T_AD(p0, h0, R0, gam0, rho0, T0, Scalar0, nScal);
				calcThermoProp_T_AD(p1, h1, R1, gam1, rho1, T1, Scalar1, nScal);
#else
				calcThermoProp_p_AD(T0, h0, R0, gam0, rho0, p0, Scalar0, nScal);
				calcThermoProp_p_AD(T1, h1, R1, gam1, rho1, p1, Scalar1, nScal);
#endif
			}


			if (kine_Index > -1) {  // save kine if defined
				kineCV0 = ScalCV0[kine_Index];         // cell center left for implicit side (Jacobi is computed with cell center)
				kineFA0 = Scalar0[kine_Index];         // cell face left, if second order
				kineCV1 = ScalCV1[kine_Index];         // cell center left for implicit side (Jacobi is computed with cell center)
				kineFA1 = Scalar1[kine_Index];         // cell face right,
			}

			// .............................................................................................
			// calculation of Euler Flux explicit using HLLC
			// .............................................................................................
			calcEulerFlux_HLLC_AD(Frho, Frhou, FrhoE, FrhoScal,
					rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
					rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
					area, nVec, nScal, 0.0);

			if(icvCenter == icv0) {
				rhs_rho -= Frho;
				for (int i = 0; i < 3; i++)
					rhs_rhou[i] -= Frhou[i];
				rhs_rhoE -= FrhoE;
				for (int iScal = 0; iScal < nScal; iScal++)
					rhs_rhoScal[iScal] -= ScalConvTerm[iScal] * FrhoScal[iScal];
			} else if(icvCenter == icv1) {
				rhs_rho += Frho;
				for (int i = 0; i < 3; i++)
					rhs_rhou[i] += Frhou[i];
				rhs_rhoE += FrhoE;
				for (int iScal = 0; iScal < nScal; iScal++)
					rhs_rhoScal[iScal] += ScalConvTerm[iScal] * FrhoScal[iScal];
			} else {
				printf("Error in IkeWithModels_AD::calcEulerFlux1D_AD() => icvCenter(%d) is not equal to icv0(%d) nor icv1(%d) \n",icvCenter, icv0, icv1);
				freeMemForCalcEulerFlux1D_AD(Apl, Ami, AplScal, AmiScal, FrhoScal, Scalar0, Scalar1, ScalCV0, ScalCV1, ScalConvTerm, nScal);
				throw(-1);
			}

			if(isNaN(rhs_rho.value()) || isNaN(rhs_rhou[0].value()) || isNaN(rhs_rhou[1].value()) || isNaN(rhs_rhou[2].value()) || isNaN(rhs_rhoE.value())) {
				printf("Error in IkeWithModels_AD::calcEulerFlux1D_AD: NaN detected in flux at internal ifa=%d(%.2e,%.2e,%.2e)\n", ifa, x_fa[ifa][0], x_fa[ifa][1], x_fa[ifa][2]);
				printf("      Details: mpi_rank=%d, icvCenter=%d, Frho=%.3e, Frhou=(%.3e,%.3e,%.3e), FrhoE=%3e\n", mpi_rank, icvCenter, Frho.value(), Frhou[0].value(),Frhou[1].value(),Frhou[2].value(), FrhoE.value());
				printf("               icv0=%d(%.2e,%.2e,%.2e); icv1=%d(%.2e,%.2e,%.2e) \n", icv0, x_cv[icv0][0],x_cv[icv0][1],x_cv[icv0][2], icv1, x_cv[icv1][0],x_cv[icv1][1],x_cv[icv1][2]);
				printf("               rho0=%.3e, u0=(%.3e, %.3e, %.3e), p0=%.3e, T0=%.3e, h0=%.3e, R0=%.3e, gam0=%.3e, kineFA0=%.3e\n", rho0.value(), u0[0].value(), u0[1].value(), u0[2].value(), p0.value(), T0.value(), h0.value(), R0.value(), gam0.value(), kineFA0.value());
				printf("               rho1=%.3e, u1=(%.3e, %.3e, %.3e), p1=%.3e, T1=%.3e, h1=%.3e, R1=%.3e, gam1=%.3e, kineFA1=%.3e\n", rho1.value(), u1[0].value(), u1[1].value(), u1[2].value(), p1.value(), T1.value(), h1.value(), R1.value(), gam1.value(), kineFA1.value());
				printf("\n");
				freeMemForCalcEulerFlux1D_AD(Apl, Ami, AplScal, AmiScal, FrhoScal, Scalar0, Scalar1, ScalCV0, ScalCV1, ScalConvTerm, nScal);
				throw(IKEWITHMODELS_ERROR_CODE);
			}


//IKJ
if(mpi_rank==0 && icv0==0 && ifa==51535) {
	printf("=> IkeWithModels_AD::calcEulerFlux1D_AD: INTERNAL ifa=%d  icvCenter=%d, mpi_rank=%d\n", ifa, icvCenter, mpi_rank);
	printf("               Frho=%.6e, Frhou=(%.4e,%.4e,%.4e), FrhoE=%.4e\n", Frho.value(), Frhou[0].value(),Frhou[1].value(),Frhou[2].value(), FrhoE.value());
	printf("               rho0=%.6e, u0=(%.6e, %e, %e), p0=%.6e, T0=%.6e, h0=%e, R0=%e, gam0=%e, kineFA0=%e\n", rho0.value(), u0[0].value(), u0[1].value(), u0[2].value(), p0.value(), T0.value(), h0.value(), R0.value(), gam0.value(), kineFA0.value());
	cout<<"               grad_rho[icv0]="<<grad_rho[icv0][0]<<", "<<grad_rho[icv0][1]<<", "<<grad_rho[icv0][2]<<endl;
//	printf("               rho1=%.4e, u1=(%.4e, %.4e, %.4e), p1=%.4e, T1=%.4e, h1=%.4e, R1=%.4e, gam1=%.4e, kineFA1=%.3e\n", rho1.value(), u1[0].value(), u1[1].value(), u1[2].value(), p1.value(), T1.value(), h1.value(), R1.value(), gam1.value(), kineFA1.value());	printf("\n");
}

			// .............................................................................................
			// calculate implicit matrix using HLLC
			// .............................................................................................
			if (flagImplicit && icv0 < ncv && icv1 < ncv_g) {
				calcEulerFluxMatrices_HLLC_AD(Apl, Ami, AplScal, AmiScal,
						rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, kineCV0,
						rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, kineCV1,
						area, nVec, nScal, 0.0);

				if (icv0 < ncv) { // if icv0 is internal...
					for (int i = 0; i < 5; i++)
						for (int j = 0; j < 5; j++) {
							A[noc00][i][j] += Apl[j][i];
							A[noc01][i][j] -= Apl[j][i];
						}
				}

				if (icv1 < ncv) { // if icv1 is internal...
					for (int i = 0; i < 5; i++)
						for (int j = 0; j < 5; j++) {
							A[noc11][i][j] -= Ami[j][i];
							A[noc10][i][j] += Ami[j][i];
						}
				}

				for (int iScal = 0; iScal < nScal; iScal++)
					for (int i = 0; i <= 5; i++) {
						if (icv0 < ncv) {
							AScal[iScal][i][noc00] += ScalConvTerm[iScal] * AplScal[iScal][i];
							AScal[iScal][i][noc01] -= ScalConvTerm[iScal] * AplScal[iScal][i];
						}

						if (icv1 < ncv) {
							AScal[iScal][i][noc11] -= ScalConvTerm[iScal] * AmiScal[iScal][i];
							AScal[iScal][i][noc10] += ScalConvTerm[iScal] * AmiScal[iScal][i];
						}
					}
			}
		}
		// ===============================================================================================
		// cycle through boundary faces, assembling flux
		// ===============================================================================================
		else {
			FaZone* zone = zofa[ifa];
			if (zone->getKind() == FA_ZONE_BOUNDARY) {
				Param *param;

				if (getParam(param, zone->getName())) {
					// .............................................................................................
					// SYMMETRY BOUNDARY CONDITION OR WALL BOUNDARY, ATTENTION: ONLY FIRST ORDER!!!
					// .............................................................................................
					if(getBoundaryType(ifa, param) == SYMMETRY || getBoundaryType(ifa, param) == WALL) {
						int icv0 = cvofa[ifa][0];
						int icv1 = cvofa[ifa][1];

						double nVec[3] = {0.0, 0.0, 0.0};
						double area = normVec3d(nVec, fa_normal[ifa]);

						for (int iScal = 0; iScal < nScal; iScal++)
							Scalar0[iScal] = scalarTranspEqVector_AD[iScal].phi[icv1];

						REALQS kineFA = 0.0;
						if (kine_Index > -1)
							kineFA = scalarTranspEqVector_AD[kine_Index].phi[icv1];

						calcEulerFlux_HLLC_AD(Frho, Frhou, FrhoE, FrhoScal,
								rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
								rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
								area, nVec, nScal, 0.0);

						rhs_rho -= Frho;
						for (int i = 0; i < 3; i++)
							rhs_rhou[i] -= Frhou[i];
						rhs_rhoE -= FrhoE;
						for (int iScal = 0; iScal < nScal; iScal++)
							rhs_rhoScal[iScal] -= ScalConvTerm[iScal] * FrhoScal[iScal];

						if(isNaN(rhs_rho.value()) || isNaN(rhs_rhou[0].value()) || isNaN(rhs_rhou[1].value()) || isNaN(rhs_rhou[2].value()) || isNaN(rhs_rhoE.value())) {
							printf("Error in IkeWithModels_AD::calcEulerFlux1D_AD: NaN detected in flux at boundary (SYMM or WALL) ifa=%d(%.2e,%.2e,%.2e)\n", ifa, x_fa[ifa][0], x_fa[ifa][1], x_fa[ifa][2]);
							printf("      Details: mpi_rank=%d, icvCenter=%d, Frho=%.3e, Frhou=(%.3e,%.3e,%.3e), FrhoE=%.3e\n", mpi_rank, icvCenter, Frho.value(), Frhou[0].value(),Frhou[1].value(),Frhou[2].value(), FrhoE.value());
							printf("               icv1=%d(%.2e,%.2e,%.2e) \n", icv1, x_cv[icv1][0],x_cv[icv1][1],x_cv[icv1][2]);
							printf("               rho=%.3e, vel=(%.3e, %.3e, %.3e), press=%.3e, temp=%.3e, enthalpy=%.3e, RoM=%.3e, gamma=%.3e, kineFA=%.3e\n", rho[icv1].value(), vel[icv1][0].value(),vel[icv1][1].value(),vel[icv1][2].value(), press[icv1].value(), temp[icv1].value(), enthalpy[icv1].value(), RoM[icv1].value(), gamma[icv1].value(), kineFA.value());
							printf("\n");
							freeMemForCalcEulerFlux1D_AD(Apl, Ami, AplScal, AmiScal, FrhoScal, Scalar0, Scalar1, ScalCV0, ScalCV1, ScalConvTerm, nScal);
							throw(IKEWITHMODELS_ERROR_CODE);
						}

//IKJ
if(mpi_rank==0 && icv0==0 && ifa==227) {
    printf("=> IkeWithModels_AD::calcEulerFlux1D_AD: SYMM or WALL ifa=%d  icvCenter=%d, mpi_rank=%d\n", ifa, icvCenter, mpi_rank);
    printf("               Frho=%.6e, Frhou=(%.4e,%.4e,%.4e), FrhoE=%.4e\n", Frho.value(), Frhou[0].value(),Frhou[1].value(),Frhou[2].value(), FrhoE.value());
    printf("               rho=%.6e, vel=(%.6e, %e, %e), press=%.6e, temp=%.6e, enthalpy=%e, RoM=%e, gamma=%e, kineFA=%e\n", rho[icv1].value(), vel[icv1][0].value(),vel[icv1][1].value(),vel[icv1][2].value(), press[icv1].value(), temp[icv1].value(), enthalpy[icv1].value(), RoM[icv1].value(), gamma[icv1].value(), kineFA.value());
    printf("\n");
}



						if (flagImplicit && icv0 < ncv) {
							calcEulerFluxMatrices_HLLC_AD(Apl, NULL, AplScal, AmiScal,
									rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
									rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
									area, nVec, nScal, 0.0);

							int noc00 = nbocv_i[icv0]; // icv0's diagonal
							for (int i = 0; i < 5; i++)
								for (int j = 0; j < 5; j++)
									A[noc00][i][j] += Apl[j][i];

							for (int iScal = 0; iScal < nScal; iScal++)
								for (int i = 0; i <= 5; i++)
									AScal[iScal][i][noc00] += ScalConvTerm[iScal] * AplScal[iScal][i];
						}
					}
					// .............................................................................................
					// OTHER BOUNDARY CONDITIONS (HOOK, DIRICHLET(CBC), NEUMANN, ...)
					// .............................................................................................
					else {
						int icv0 = cvofa[ifa][0];
						int icv1 = cvofa[ifa][1];

						double nVec[3] = {0.0, 0.0, 0.0};
						double area = normVec3d(nVec, fa_normal[ifa]);

						REALQ rho0 = rho[icv0];
						REALQ u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
						REALQ p0 = press[icv0];
						REALQ T0 = temp[icv0];
						REALQ h0 = enthalpy[icv0];
						REALQ gam0 = gamma[icv0];
						REALQ R0 = RoM[icv0];
						REALQS kineCV0 = 0.0;           // cell center
						REALQS kineFA0 = 0.0;           // cell face

						REALQS kineFA1 = 0.0;

						for (int iScal = 0; iScal < nScal; iScal++) {
							ScalCV0[iScal] = Scalar0[iScal] = scalarTranspEqVector_AD[iScal].phi[icv0];
							Scalar1[iScal] = scalarTranspEqVector_AD[iScal].phi[icv1];
						}

						if (sndOrder == true) {
							double r0[3] = {0.0, 0.0, 0.0};
							vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);

							// left side
							rho0 += vecDotVec3d_AD(r0, grad_rho[icv0]);
#ifdef temp_reconstruction
							T0 += vecDotVec3d_AD(r0, grad_temp[icv0]);
	#ifdef USE_CONDASSIGN
							if ((T0.value() <= 0.0) || (rho0.value() <= 0.0)) {
								T0 = temp[icv0];
								rho0 = rho[icv0];
								myCountReducedOrder++;
							}
	#else
							if ((T0 <= 0.0) || (rho0 <= 0.0)) {
								T0 = temp[icv0];
								rho0 = rho[icv0];
								myCountReducedOrder++;
							}
	#endif
#else
							p0 += vecDotVec3d_AD(r0, grad_p[icv0]);
	#ifdef USE_CONDASSIGN
							if ((p0.value() <= 0.0) || (rho0.value() <= 0.0)) {
								p0 = press[icv0];
								rho0 = rho[icv0];
								myCountReducedOrder++;
							}
	#else
							if ((p0 <= 0.0) || (rho0 <= 0.0)) {
								p0 = press[icv0];
								rho0 = rho[icv0];
								myCountReducedOrder++;
							}
	#endif
#endif
							else {
								for (int i = 0; i < 3; i++)
									u0[i] += vecDotVec3d_AD(r0, grad_u[icv0][i]);
								for (int iScal = 0; iScal < nScal; iScal++) {
									if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
										Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d_AD(r0, scalarTranspEqVector_AD[iScal].grad_phi[icv0])) / rho0;    // CHANGE FOR SCALAR AD
									else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
										Scalar0[iScal] += vecDotVec3d_AD(r0, scalarTranspEqVector_AD[iScal].grad_phi[icv0]);
								}
							}

							// calculation of other variables at faces: p/T, h, R, gam
#ifdef temp_reconstruction
							calcThermoProp_T_AD(p0, h0, R0, gam0, rho0, T0, Scalar0, nScal);
#else
							calcThermoProp_p_AD(T0, h0, R0, gam0, rho0, p0, Scalar0, nScal);
#endif
						}

						if (kine_Index > -1) {  // save kine if defined
							kineCV0 = ScalCV0[kine_Index];          // cell center
							kineFA0 = Scalar0[kine_Index];                                 // cell face
							kineFA1 = Scalar1[kine_Index];
						}

						calcEulerFlux_HLLC_AD(Frho, Frhou, FrhoE, FrhoScal,
								rho0,      u0,        p0,          T0,         h0,             R0,        gam0,        Scalar0, kineFA0,
								rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar1, kineFA1,
								area, nVec, nScal, 0.0);

						rhs_rho -= Frho;
						for (int i = 0; i < 3; i++)
							rhs_rhou[i] -= Frhou[i];
						rhs_rhoE -= FrhoE;
						for (int iScal = 0; iScal < nScal; iScal++)
							rhs_rhoScal[iScal] -= ScalConvTerm[iScal] * FrhoScal[iScal];

						if(isNaN(rhs_rho.value()) || isNaN(rhs_rhou[0].value()) || isNaN(rhs_rhou[1].value()) || isNaN(rhs_rhou[2].value()) || isNaN(rhs_rhoE.value())) {
							printf("Error in IkeWithModels_AD::calcEulerFlux1D_AD: NaN detected in flux at boundary (HOOK, DIRI(CBC), NEUM) ifa=%d(%.2e,%.2e,%.2e)\n", ifa, x_fa[ifa][0], x_fa[ifa][1], x_fa[ifa][2]);
							printf("      Details: mpi_rank=%d, icvCenter=%d, Frho=%.3e, Frhou=(%.3e,%.3e,%.3e), FrhoE=%.3e\n", mpi_rank, icvCenter, Frho.value(), Frhou[0].value(),Frhou[1].value(),Frhou[2].value(), FrhoE.value());
							printf("               icv1=%d(%.2e,%.2e,%.2e) \n", icv1, x_cv[icv1][0],x_cv[icv1][1],x_cv[icv1][2]);
							printf("               rho0=%.3e, u0=(%.3e, %.3e, %.3e), p0=%.3e, T0=%.3e, h0=%.3e, R0=%.3e, gam0=%.3e, kineFA0=%.3e\n", rho0.value(), u0[0].value(), u0[1].value(), u0[2].value(), p0.value(), T0.value(), h0.value(), R0.value(), gam0.value(), kineFA0.value());
							printf("               rho=%.3e, vel=(%.3e, %.3e, %.3e), press=%.3e, temp=%.3e, enthalpy=%.3e, RoM=%.3e, gamma=%.3e, kineFA=%.3e\n", rho[icv1].value(), vel[icv1][0].value(),vel[icv1][1].value(),vel[icv1][2].value(), press[icv1].value(), temp[icv1].value(), enthalpy[icv1].value(), RoM[icv1].value(), gamma[icv1].value(), kineFA1.value());
							printf("\n");
							freeMemForCalcEulerFlux1D_AD(Apl, Ami, AplScal, AmiScal, FrhoScal, Scalar0, Scalar1, ScalCV0, ScalCV1, ScalConvTerm, nScal);
							throw(IKEWITHMODELS_ERROR_CODE);
						}

if(mpi_rank==0 && icv0==0) {
    printf("=> IkeWithModels_AD::calcEulerFlux1D_AD: flux at boundary (HOOK, DIRI(CBC), NEUM) ifa=%d(%.2e,%.2e,%.2e)\n", ifa, x_fa[ifa][0], x_fa[ifa][1], x_fa[ifa][2]);
    printf("      Details: mpi_rank=%d, icvCenter=%d, Frho=%.3e, Frhou=(%.3e,%.3e,%.3e), FrhoE=%3e\n", mpi_rank, icvCenter, Frho.value(), Frhou[0].value(),Frhou[1].value(),Frhou[2].value(), FrhoE.value());
    printf("               icv1=%d(%.2e,%.2e,%.2e) \n", icv1, x_cv[icv1][0],x_cv[icv1][1],x_cv[icv1][2]);
    printf("               rho0=%.3e, u0=(%.3e, %.3e, %.3e), p0=%.3e, T0=%.3e, h0=%.3e, R0=%.3e, gam0=%.3e, kineFA0=%.3e\n", rho0.value(), u0[0].value(), u0[1].value(), u0[2].value(), p0.value(), T0.value(), h0.value(), R0.value(), gam0.value(), kineFA0.value());
    printf("               rho=%.3e, vel=(%.3e, %.3e, %.3e), press=%.3e, temp=%.3e, enthalpy=%.3e, RoM=%.3e, gamma=%.3e, kineFA=%.3e\n", rho[icv1].value(), vel[icv1][0].value(),vel[icv1][1].value(),vel[icv1][2].value(), press[icv1].value(), temp[icv1].value(), enthalpy[icv1].value(), RoM[icv1].value(), gamma[icv1].value(), kineFA1.value());
    printf("\n");
}



						if (flagImplicit && icv0 < ncv) {
							calcEulerFluxMatrices_HLLC_AD(Apl, NULL, AplScal, AmiScal,
									rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, kineCV0,
									rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar1, kineFA1,
									area, nVec, nScal, 0.0);

							int noc00 = nbocv_i[icv0]; // icv0's diagonal
							for (int i = 0; i < 5; i++)
								for (int j = 0; j < 5; j++)
									A[noc00][i][j] += Apl[j][i];

							for (int iScal = 0; iScal < nScal; iScal++)
								for (int i = 0; i <= 5; i++)
									AScal[iScal][i][noc00] += ScalConvTerm[iScal] * AplScal[iScal][i];
						}
					}
				}
			}
		}
	}

	freeMemForCalcEulerFlux1D_AD(Apl, Ami, AplScal, AmiScal, FrhoScal, Scalar0, Scalar1, ScalCV0, ScalCV1, ScalConvTerm, nScal);

	firstCall = false;
	return myCountReducedOrder;
}

/*
 * Method: calcViscousFluxNS1D_AD
 * ------------------------------
 * Calculate viscous flux in the so-called "1D style"
 * Original code =calcViscousFluxNS_AD() in JoeWithModels.cpp
 */
void IkeWithModels_AD::calcViscousFluxNS1D_AD(const int icvCenter, REALA &rhs_rho, REALA rhs_rhou[3], REALA &rhs_rhoE,
		ADscalar<REALQ> &rho, ADvector<REALQ> &rhou, ADscalar<REALQ> &rhoE,
		double (*A)[5][5], int flagImplicit) {
	double (*A0)[5];
	double (*A1)[5];

	if (flagImplicit) {
		A0 = new double[5][5];
		A1 = new double[5][5];

		for (int i = 0; i < 5; i++)
			for (int j = 0; j < 5; j++)
				A0[i][j] = A1[i][j] = 0.0;
	} else {
		A0 = A1 = NULL;
	}

	REALQ Frhou[3] = {0.0, 0.0, 0.0}, FrhoE = 0.0;

	// save the index of kine if defined
	int kine_index = -1;
	for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
//		if (scalarTranspEqVector[iScal].getName() == "kine")
		if (strcmp(scalarTranspEqVector[iScal].getName(),"kine")==0)
			kine_index = iScal;

	// ====================================================================
	//        compute gradients, with boundary values
	// ====================================================================
	calcCv2Grad1D_AD(icvCenter, grad_enthalpy, enthalpy, limiterNavierS, enthalpy, epsilonSDWLS);

	// Calculate flux
	int nofa_f = faocv_i[icvCenter];
	int nofa_l = faocv_i[icvCenter+1]-1;
	for (int foc = nofa_f; foc <= nofa_l; foc++) {
		int ifa = faocv_v[foc];

		// ====================================================================
		// cycle through internal faces, assembling flux and matrix
		// ====================================================================
		if((ifa >= nfa_b && ifa<nfa) || (ifa >= nfa_b2 && ifa < nfa_b2gg)) {
			int icv0 = cvofa[ifa][0];
			int icv1 = cvofa[ifa][1];

			int noc00, noc01, noc11, noc10;
			if (flagImplicit && icv0 < ncv && icv1 < ncv_g)
				getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

			// face unit normal and area...
			double nVec[3] = {0.0, 0.0, 0.0};
			double area = normVec3d(nVec, fa_normal[ifa]);
			double sVec[3] = {0.0, 0.0, 0.0};
			vecMinVec3d(sVec, x_cv[icv1], x_cv[icv0]);
			double smag = normVec3d(sVec);

			double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
			vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
			vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
			double w0 = sqrt(vecDotVec3d(dx0, dx0));
			double w1 = sqrt(vecDotVec3d(dx1, dx1));
			double ws = w0 + w1;
			w0 /= ws;
			w1 /= ws;
			REALQ uAux_fa[3] = { w1*vel[icv0][0]+ w0*vel[icv1][0],
					w1*vel[icv0][1]+ w0*vel[icv1][1],
					w1*vel[icv0][2]+ w0*vel[icv1][2]};

			// kinetic energy, if defined
			REALQS kine0 = 0.0;
			REALQS kine1 = 0.0;
			REALQS kine_fa = 0.0;
			if (kine_index > -1) {
				kine0 = scalarTranspEqVector_AD[kine_index].phi[icv0];
				kine1 = scalarTranspEqVector_AD[kine_index].phi[icv1];
				kine_fa = w1*kine0 + w0*kine1;
			}

			// calculate viscous flux
#ifdef USE_ARTIF_VISC
		// If the user wants to have artificial bulk viscosity, pass the information to the addViscFlux() method
		// (If the user wants to use artificial viscosity (not bulk only), it was already reflected in mul_fa)
			REALQS artifBulkViscosity_AD = 0.0;
			if(turnOnArtifVisc && artifVisc_bulkViscOnly)
				artifBulkViscosity_AD = w1*artifViscMag_AD[icv0] + w0*artifViscMag_AD[icv1];

			addViscFlux_AD(Frhou, FrhoE, A0, A1,
					rho[icv0], vel[icv0], grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0], gamma[icv0], kine0,
					rho[icv1], vel[icv1], grad_u[icv1], enthalpy[icv1], grad_enthalpy[icv1], temp[icv1], RoM[icv1], gamma[icv1], kine1,
					mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kine_fa, uAux_fa,
					area, nVec, smag, sVec, artifBulkViscosity_AD);
#else
			addViscFlux_AD(Frhou, FrhoE, A0, A1,
					rho[icv0], vel[icv0], grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0], gamma[icv0], kine0,
					rho[icv1], vel[icv1], grad_u[icv1], enthalpy[icv1], grad_enthalpy[icv1], temp[icv1], RoM[icv1], gamma[icv1], kine1,
					mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kine_fa, uAux_fa,
					area, nVec, smag, sVec);
#endif

			if (flagImplicit && icv0 < ncv && icv1 < ncv_g) {
				for (int i=0; i<5; i++)
					for (int j=0; j<5; j++) {
						A[noc00][i][j] -= A0[j][i];
						A[noc01][i][j] += A0[j][i];
					}

				if (icv1 < ncv) { // if icv1 is internal...
					for (int i=0; i<5; i++)
						for (int j=0; j<5; j++) {
							A[noc11][i][j] += A1[j][i];
							A[noc10][i][j] -= A1[j][i];
						}
				}
			}

			if(icvCenter == icv0) {
				for (int i = 0; i < 3; i++)
					rhs_rhou[i] -= Frhou[i];
				rhs_rhoE -= FrhoE;
			} else if (icv1 == icv1) {
				for (int i = 0; i < 3; i++)
					rhs_rhou[i] += Frhou[i];
				rhs_rhoE += FrhoE;
			}
		}
		// ====================================================================
		// cycle through boundary faces, assembling flux and matrix
		// ====================================================================
		else {
			FaZone* zone = zofa[ifa];
			// .............................................................................................
			// SYMMETRY BOUNDARY CONDITION
			// .............................................................................................
			if(getBoundaryType(ifa, zone) == SYMMETRY) {
				if (kine_index > -1) {
					int icv0 = cvofa[ifa][0];
					int icv1 = cvofa[ifa][1];
					assert( icv0 >= 0 );
					double nVec[3] = {0.0, 0.0, 0.0};
					double area = normVec3d(nVec, fa_normal[ifa]);

					adouble kine_fa  = scalarTranspEqVector_AD[kine_index].phi[icv1];

					adouble tmp = 1.0/3.0*(rho[icv0] + rho[icv1])*kine_fa;

					for (int i = 0; i < 3; i++)
						rhs_rhou[i] -= tmp*fa_normal[ifa][i];
				}
			}
			// .............................................................................................
			// WALL BOUNDARY CONDITION
			// .............................................................................................
			else if(getBoundaryType(ifa, zone) == WALL) {
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];

				REALX nVec[3] = {0.0, 0.0, 0.0};
				REALX area = normVec3d(nVec, fa_normal[ifa]);
				REALX sVec[3] = {0.0, 0.0, 0.0};
				vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
				REALX smag = fabs(vecDotVec3d(sVec, nVec));  // project sVec to wall face normal

				REALQS kine0 = 0.0;
				REALQS kine1 = 0.0;
				REALQS kine_fa = 0.0;
				if (kine_index > -1)
					kine0 = scalarTranspEqVector_AD[kine_index].phi[icv0];

				// calculate viscous flux
				addViscFlux_AD(Frhou, FrhoE, A0, NULL,
						rho[icv0],    vel[icv0],    grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0],    gamma[icv0],  kine0,
						rho[icv1],    vel[icv1],    grad_u[icv0], enthalpy[icv1], grad_enthalpy[icv0], temp[icv1], RoM[icv1],    gamma[icv1],  kine1,
						mul_fa[ifa], 0.0, lamOcp_fa[ifa], kine_fa, vel[icv1],
						area, nVec, smag, nVec);  /* <- use nVec here instead of sVec, to avoid inaccurate correction*/

				if (flagImplicit && icv0<ncv) {
					int noc00 = nbocv_i[icv0]; // icv0's diagonal

					for (int i=0; i<5; i++)
						for (int j=0; j<5; j++)
							A[noc00][i][j] -= A0[j][i];
				}

				if(icvCenter == icv0) {
					for (int i = 0; i < 3; i++)
						rhs_rhou[i] -= Frhou[i];
					rhs_rhoE -= FrhoE;
				}
			}
			// .............................................................................................
			// OTHER BOUNDARY CONDITIONS
			// .............................................................................................
			else {
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];

				REALX nVec[3] = {0.0, 0.0, 0.0};
				REALX area = normVec3d(nVec, fa_normal[ifa]);
				REALX sVec[3] = {0.0, 0.0, 0.0};
				vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
				REALX smag = normVec3d(sVec);

				REALQS kine0 = 0.0;
				REALQS kine1 = 0.0;
				REALQS kine_fa = 0.0;
				if (kine_index > -1) {
					kine0 = scalarTranspEqVector_AD[kine_index].phi[icv0];
					kine1 = scalarTranspEqVector_AD[kine_index].phi[icv1];
					kine_fa = kine1;
				}

				// calculate viscous flux
				addViscFlux_AD(Frhou, FrhoE, A0, NULL,
						rho[icv0],    vel[icv0],    grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0],    gamma[icv0],  kine0,
						rho[icv1],    vel[icv1],    grad_u[icv0], enthalpy[icv1], grad_enthalpy[icv0], temp[icv1], RoM[icv1],    gamma[icv1],  kine1,
						mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kine_fa, vel[icv1],
						area, nVec, smag, sVec);

				if (flagImplicit && icv0<ncv) {
					int noc00 = nbocv_i[icv0]; // icv0's diagonal

					for (int i=0; i<5; i++)
						for (int j=0; j<5; j++)
							A[noc00][i][j] -= A0[j][i];
				}

				if(icvCenter == icv0) {
					for (int i = 0; i < 3; i++)
						rhs_rhou[i] -= Frhou[i];
					rhs_rhoE -= FrhoE;
				}
			}
		}
	}

	if (A0  != NULL)  delete [] A0;
	if (A1  != NULL)  delete [] A1;
}

///*
// * Method: calcViscousFluxScalar_new1D_AD
// * --------------------------------------
// * Original code: calcViscousFluxScalar_new_AD in ScalarsAD.cpp
// */
//void IkeWithModels_AD::calcViscousFluxScalar_new1D_AD(const int icvCenter, const int iScal, adouble &rhs_rhoScal, double *AScal, ScalarTranspEq &transpScal, ScalarTranspEq_AD &transpScal_AD, int flagImplicit) {
//	if (!transpScal.diffTerm) { // Note: transScal_AD doesn't have "diffTerm"
//		return;
//	}
//
//	string scalName(transpScal.getName());
//	double phiBCval = 0.0, phiBCflux = 0.0;
//
//	adouble *phi = transpScal_AD.phi;
//	adouble (*grad_phi)[3] = transpScal_AD.grad_phi;
//	adouble *diff = transpScal_AD.diff;
//
//	// =============================================================================================
//	// compute gradients, with boundary values
//	// =============================================================================================
//	calcCv2Grad1D_AD(icvCenter, grad_phi, phi, limiterNavierS, phi, epsilonSDWLS);
//
//	// =============================================================================================
//	// compute user-defined diffusivity
//	// =============================================================================================
//	vector<int> faInternal;
//	vector<int> faBoundary;
//	//find2LayerFaceIndices(icvCenter, faInternal, faBoundary);
//	find1LayerFaceIndices(icvCenter, faInternal, faBoundary);
//	// internal faces
//	for(size_t index=0; index<faInternal.size(); ++index) {
//		int ifa = faInternal[index];
//		FaZone* zone = zofa[ifa];
//
//		diffusivityHookScalarRansTurb1D_AD(ifa, zone, scalName);
//		diffusivityHookScalarRansComb1D_AD(ifa, zone, scalName);
//	}
//	for(size_t index=0; index<faBoundary.size(); ++index) {
//		int ifa = faBoundary[index];
//		FaZone* zone = zofa[ifa];
//		diffusivityHookScalarRansTurb1D_AD(ifa, zone, scalName);
//		diffusivityHookScalarRansComb1D_AD(ifa, zone, scalName);
//	}
//
//	// Calculate flux
//	int nofa_f = faocv_i[icvCenter];
//	int nofa_l = faocv_i[icvCenter+1]-1;
//	for (int nofa = nofa_f; nofa <= nofa_l; nofa++) {
//		int ifa = faocv_v[nofa];
//
//		// =============================================================================================
//		// cycle trough internal faces and compute viscous flux and implicit matrix
//		// =============================================================================================
//		if((ifa >= nfa_b && ifa<nfa) || (ifa >= nfa_b2 && ifa < nfa_b2gg)) {
//			int icv0 = cvofa[ifa][0];
//			int icv1 = cvofa[ifa][1];
//			int noc00, noc01, noc11, noc10;
//
//			if (flagImplicit && icv0 < ncv && icv1 < ncv_g)
//				getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);
//
//			double n[3] = {0.0, 0.0, 0.0}, s[3] = {0.0, 0.0, 0.0};
//			double nmag = normVec3d(n, fa_normal[ifa]);
//			vecMinVec3d(s, x_cv[icv1], x_cv[icv0]);
//			double smag = normVec3d(s);
//			double alpha = vecDotVec3d(n, s);
//			assert((alpha > 0.0) && (alpha < 1.0+1.0E-12));   // alpha should now contain s dot n...
//
//			double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
//			vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
//			vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
//
//			adouble viscFlux = diff[ifa] * nmag * (alpha * (phi[icv1] - phi[icv0]) / smag
//					+ 0.5 * (grad_phi[icv0][0] + grad_phi[icv1][0]) * (n[0] - alpha * s[0])
//					+ 0.5 * (grad_phi[icv0][1] + grad_phi[icv1][1]) * (n[1] - alpha * s[1])
//					+ 0.5 * (grad_phi[icv0][2] + grad_phi[icv1][2]) * (n[2] - alpha * s[2]));
//
//			if(icvCenter == icv0)
//				rhs_rhoScal += viscFlux;
//			else if(icvCenter == icv1)
//				rhs_rhoScal -= viscFlux;
//			else {
//				printf("ERROR in IkeWithModels_AD::calcViscousFluxScalar_new1D_AD() => icvCenter=%d is not equal to icv0=%d nor icv1=%d \n", icvCenter, icv0, icv1);
//			}
//
//			if (flagImplicit && icv0 < ncv && icv1 < ncv_g) {
//				AScal[noc00] +=   diff[ifa].value() * nmag * alpha / smag / rho_AD[icv0].value();
//				AScal[noc01] += - diff[ifa].value() * nmag * alpha / smag / rho_AD[icv0].value();
//
//				if (icv1 < ncv) {
//					AScal[noc11] +=   diff[ifa].value() * nmag * alpha / smag / rho_AD[icv1].value();
//					AScal[noc10] += - diff[ifa].value() * nmag * alpha / smag / rho_AD[icv1].value();
//				}
//			}
//		}
//		// =============================================================================================
//		// cycle trough boundary faces and compute viscous flux and implicit matrix
//		// =============================================================================================
//		else {
//			FaZone* zone = zofa[ifa];
//			// .............................................................................................
//			// HOOK BOUNDARY CONDITION
//			// .............................................................................................
//			if(getBoundaryTypeScalar(ifa, zone, scalName) == HOOK_SCALAR)  {
//				int icv0 = cvofa[ifa][0];
//				int icv1 = cvofa[ifa][1];
//				int noc00 = nbocv_i[icv0];
//
//				double n[3] = {0.0, 0.0, 0.0}, s_half[3] = {0.0, 0.0, 0.0};
//				double nmag = normVec3d(n, fa_normal[ifa]);
//				vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
//				double smag_half = fabs(vecDotVec3d(s_half, n));   //normVec3d(s_half);
//				//double alpha = vecDotVec3d(n, s_half);           // HACK alfa should be always 1 at the boundary (provide such a grid!)
//
//				adouble viscFlux = diff[ifa] * nmag * (phi[icv1] - phi[icv0]) / smag_half;
//
//				if(icvCenter == icv0)
//					rhs_rhoScal += viscFlux;
//
//				if (flagImplicit && icv0<ncv)
//					AScal[noc00] += diff[ifa].value() * nmag / smag_half / rho_AD[icv0].value();
//			}
//			// .............................................................................................
//			// DIRICHLET BOUNDARY CONDITION
//			// .............................................................................................
//			else if (getBoundaryTypeScalar(ifa, zone, scalName) == DIRICHLET_SCALAR){
//			    Param *p;
//			    string fName = zone->getName();
//			    if (getParam(p, fName + "." + scalName))
//			    	phiBCval = p->getDouble(1);
//			    else {
//			    	cout<<"Error in IkeWithModels_AD::calcViscousFluxScalar_new1D_AD() => Cannot find Param for "<<fName+"."+scalName<<" at ifa="<<ifa<<endl;
//			    	throw(-1);
//			    }
//
//				int icv0 = cvofa[ifa][0];
//				int icv1 = cvofa[ifa][1];
//				int noc00 = nbocv_i[icv0];
//
//				double n[3] = {0.0, 0.0, 0.0}, s_half[3] = {0.0, 0.0, 0.0};
//				double nmag = normVec3d(n, fa_normal[ifa]);
//				vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
//				double smag_half = fabs(vecDotVec3d(s_half, n));   //normVec3d(s_half);
//				//double alpha = vecDotVec3d(n, s_half);           // HACK alfa should be always 1 at the boundary (provide such a grid!)
//
//				adouble viscFlux = diff[ifa] * nmag * (phi[icv1] - phi[icv0]) / smag_half;                // TODO: phi_bfa or phiBCval??
//
//				if(icvCenter == icv0)
//					rhs_rhoScal += viscFlux;
//
//				if (flagImplicit && icv0<ncv)
//					AScal[noc00] += diff[ifa].value() * nmag / smag_half / rho_AD[icv0].value();
//			}
//			// .............................................................................................
//			// FLUX BOUNDARY CONDITION: NOT IMPLEMENTED YET
//			// .............................................................................................
//			else if (getBoundaryTypeScalar(ifa, zone, scalName) == FLUX_SCALAR) {
////				Param *p;
////				string fName = zone->getName();
////				if (getParam(p, fName + "." + scalName))
////					phiBCflux = p->getDouble(1);
////				else {
////					cout<<"Error in IkeWithModels_AD::calcViscousFluxScalar_new1D_AD() => Cannot find Param for "<<fName+"."+scalName+".flux"<<" at ifa="<<ifa<<endl;
////					throw(-1);
////				}
//				cerr << "scalarZoneIsFlux not implemented yet!" << endl;
//				throw(-1);
//			}
//			// .............................................................................................
//			// OTHER BOUNDARY CONDITIONS: NO VISCOUS FLUX
//			// .............................................................................................
//			else {
//				// note: SYMMETRY
//			}
//		}
//	}
//}
#endif

/*
 * Method: calcRhs1D_AD
 * --------------------
 * Calculate the Rhs of the N-S equations in the so-called "1D-style"
 * Original code = calcRhs_AD() in JoeWithModelsAD.cpp
 */
int IkeWithModels_AD::calcRhs1D_AD(int icvCenter, REALA &rhs_rho, REALA rhs_rhou[3], REALA &rhs_rhoE, REALAS *rhs_rhoScal, REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE, double (*A)[5][5], double ***AScal, int flagImplicit) {
	static bool firstCall = true;

	// set RHS to zero
	rhs_rho = 0.0;
	for (int i = 0; i < 3; i++)
		rhs_rhou[i] = 0.0;
	rhs_rhoE = 0.0;

	// set scalars RHS to zero
	for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
		for (int icv = 0; icv < ncv_gg; icv++)
			rhs_rhoScal[iScal] = 0.0;

	// =======================================================================================
	// NAVIER-STOKES
	// =======================================================================================
	// compute Euler Flux for NS and scalars
	int myCountReducedOrder = calcEulerFlux1D_AD(icvCenter, rhs_rho, rhs_rhou, rhs_rhoE, rhs_rhoScal, rho, rhou, rhoE, A, AScal, flagImplicit);

	// compute viscous Flux for NS
#ifdef USE_ARTIF_VISC_WITH_MEM_SAVING
	if (mu_ref > 0.0 || turnOnArtifVisc)
		calcViscousFluxNS1D_AD(icvCenter, rhs_rho, rhs_rhou, rhs_rhoE, rho, rhou, rhoE, A, flagImplicit);
#else
	if (mu_ref > 0.0)
		calcViscousFluxNS1D_AD(icvCenter, rhs_rho, rhs_rhou, rhs_rhoE, rho, rhou, rhoE, A, flagImplicit);
#endif

	// add source terms to RHS of Navier-Stokes equations
	sourceHook1D_AD(icvCenter, rhs_rho, rhs_rhou, rhs_rhoE, A);
	sourceHookRansTurb1D_AD(icvCenter, rhs_rho, rhs_rhou, rhs_rhoE, A);
	sourceHookRansComb1D_AD(icvCenter, rhs_rho, rhs_rhou, rhs_rhoE, A);

	// =======================================================================================
	// SCALARS
	// =======================================================================================
	for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++) {
		// compute viscous Flux for scalars and add source terms to RHS of scalar equations
		if (!scalarTranspEqVector[iScal].diffTerm) {
			/* Do nothing */
		} else {
			/*** From here ... */
			string scalName(scalarTranspEqVector[iScal].getName());
			double phiBCvalDummy = 0.0, phiBCfluxDummy = 0.0;

			for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
				if (zone->getKind() == FA_ZONE_BOUNDARY) {
					if (scalarZoneIsHook(zone->getName(), scalName)) {
						if(firstCall && mpi_rank == 0)
							cout << scalName << ": " << zone->getName() << " is HOOK" << endl;
					} else if (scalarZoneIsDirichlet(phiBCvalDummy, zone->getName(), scalName)) {
						if(firstCall && mpi_rank == 0)
							cout << scalName << ": " << zone->getName() << " is Dirichlet: VALUE: " << phiBCvalDummy << endl;
					} else if (scalarZoneIsFlux(phiBCfluxDummy, zone->getName(), scalName)) {
						cerr << "scalarZoneIsFlux not implemented yet!" << endl;
						throw(-1);
					} else {
						if(firstCall && mpi_rank == 0)
							cout << scalName << ": " << zone->getName() << " is Symmetry" << endl;
					}
				}
			}
			/*** ... to here, this part comes from calcViscousFluxScalar_new_AD() in ScalarsAD.cpp */

			if (AScal == NULL) {
				if (mu_ref > 0.0)
					calcViscousFluxScalar_new1D_AD(icvCenter, iScal, rhs_rhoScal[iScal], NULL, scalarTranspEqVector[iScal], scalarTranspEqVector_AD[iScal], flagImplicit);
			} else {
				if (mu_ref > 0.0)
					calcViscousFluxScalar_new1D_AD(icvCenter, iScal, rhs_rhoScal[iScal], AScal[iScal][5], scalarTranspEqVector[iScal], scalarTranspEqVector_AD[iScal], flagImplicit);
			}
		}

		if (AScal == NULL) {
			sourceHookScalarRansTurb_new1D_AD(icvCenter, rhs_rhoScal[iScal], NULL, scalarTranspEqVector[iScal].getName(), flagImplicit);
			sourceHookScalarRansComb_new1D_AD(icvCenter, rhs_rhoScal[iScal], NULL, scalarTranspEqVector[iScal].getName(), flagImplicit);
		} else {
			sourceHookScalarRansTurb_new1D_AD(icvCenter, rhs_rhoScal[iScal], AScal[iScal][5], scalarTranspEqVector[iScal].getName(), flagImplicit);
			sourceHookScalarRansComb_new1D_AD(icvCenter, rhs_rhoScal[iScal], AScal[iScal][5], scalarTranspEqVector[iScal].getName(), flagImplicit);
		}
	}

	firstCall = false;

	return myCountReducedOrder;
}

/*
 * Method: calcEulerFlux1D_AD
 * --------------------------
 * Calculate Euler flux in the so-called "1D style"
 * Return: the number of times switched back to first order at faces
 * Original code = calcEulerFlux_AD() in JoeWithModelsAD.cpp
 */
int IkeWithModels_AD::calcEulerFlux1D_AD(const int icvCenter, REALA &rhs_rho, REALA rhs_rhou[3], REALA &rhs_rhoE, REALAS *rhs_rhoScal,
		REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE, double (*A)[5][5], double ***AScal, int flagImplicit) {
	int nScal = scalarTranspEqVector.size();

	double (*Apl)[5] = NULL;
	double (*Ami)[5] = NULL;
	if (flagImplicit) {
		Apl = new double[5][5];
		Ami = new double[5][5];

		for (int i = 0; i < 5; i++)
			for (int j = 0; j < 5; j++)
				Apl[i][j] = Ami[i][j] = 0.0;
	}

	// Implicit matrix for scalars: the definition is not easy, so ask an expert (e.g., R. Pecnik)
	double (*AplScal)[6] = NULL;
	double (*AmiScal)[6] = NULL;
	if (flagImplicit){
		if (nScal > 0) AplScal = new double[nScal][6];
		if (nScal > 0) AmiScal = new double[nScal][6];

		for (int iScal = 0; iScal < nScal; iScal++)
			for (int i = 0; i <= 5; i++)
				AplScal[iScal][i] = AmiScal[iScal][i] = 0.0;
	}

	REALQ Frho, Frhou[3], FrhoE;
	REALQS *FrhoScal     = NULL;         if (nScal > 0) FrhoScal     = new REALQS[nScal];
	REALQS *Scalar0      = NULL;         if (nScal > 0) Scalar0      = new REALQS[nScal];            // cell face if second order
	REALQS *Scalar1      = NULL;         if (nScal > 0) Scalar1      = new REALQS[nScal];
	REALQS *ScalCV0      = NULL;         if (nScal > 0) ScalCV0      = new REALQS[nScal];            // cell face if second order
	REALQS *ScalCV1      = NULL;         if (nScal > 0) ScalCV1      = new REALQS[nScal];
	double *ScalConvTerm = NULL;         if (nScal > 0) ScalConvTerm = new double[nScal];            // 0 if convective term no considered, otherwise 1

	// count how many cells switched back to first order due to extrapolated negative density or pressure/temperature at the faces
	int myCountReducedOrder = 0;

	// save the index of kine if defined and save convTerm for speed
	int kine_Index = getScalarTransportIndex("kine");
	for (int iScal = 0; iScal < nScal; iScal++)
		ScalConvTerm[iScal] = double(scalarTranspEqVector[iScal].convTerm);

	// =============================================================================================
	// compute gradients, with boundary values
	// =============================================================================================
	if (sndOrder == true) {
		calcCv2Grad1D_AD(icvCenter, grad_rho, rho, limiterNavierS, rho, epsilonSDWLS);
#ifdef temp_reconstruction
		calcCv2Grad1D_AD(icvCenter, grad_temp, temp, limiterNavierS, temp, epsilonSDWLS);
#else
		calcCv2Grad1D_AD(icvCenter, grad_p, press, limiterNavierS, press, epsilonSDWLS);
#endif

		for (int iScal = 0; iScal < nScal; iScal++) {
			if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE") {
				// For the scalar reconstruction, Grad(rho*Phi) is required
				// Gradients of rho*Phi are saved in grad_phi temporarily
				// Gradients of rho*Phi are also limited like rho with alpha_rho
				// Boundary face values rho_bfa*Phi_fa are saved in Phi_fa temporarily
				REALQS *rhoPhi = new REALQS[ncv_ggff];
#ifndef USE_MEM_SAVING_ADVAR
				REALQS *phi = scalarTranspEqVector_AD[iScal].phi;
				REALQS (*grad_phi)[3] = scalarTranspEqVector_AD[iScal].grad_phi;
#endif

				// Compute rho*Phi
				int noc_f_center = nbocv_all_i[icvCenter];
				int noc_l_center = nbocv_all_i[icvCenter+1]-1;
				for (int noc_center = noc_f_center; noc_center <= noc_l_center; ++noc_center) {
					int icv = nbocv_all_v[noc_center]; // Neighboring CVs
#ifdef USE_MEM_SAVING_ADVAR
					rhoPhi[icv] = rho[icv] * scalarTranspEqVector_AD[iScal].phi[icv];
#else
					rhoPhi[icv] = rho[icv] * phi[icv];
#endif
				}

				// Compute gradients of rho*Phi and limit based on rho*Phi
#ifdef USE_MEM_SAVING_ADVAR
				calcCv2Grad1D_AD(icvCenter, scalarTranspEqVector_AD[iScal].grad_phi, rhoPhi, limiterNavierS, rhoPhi, epsilonSDWLS);
#else
				calcCv2Grad1D_AD(icvCenter, grad_phi, rhoPhi, limiterNavierS, rhoPhi, epsilonSDWLS);
#endif

				delete [] rhoPhi;
			}
			else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD") {
#ifdef USE_MEM_SAVING_ADVAR
				calcCv2Grad1D_AD(icvCenter, scalarTranspEqVector_AD[iScal].grad_phi, scalarTranspEqVector_AD[iScal].phi, limiterNavierS, scalarTranspEqVector_AD[iScal].phi, epsilonSDWLS);
#else
				REALQS *phi = scalarTranspEqVector_AD[iScal].phi;
				REALQS (*grad_phi)[3] = scalarTranspEqVector_AD[iScal].grad_phi;
				calcCv2Grad1D_AD(icvCenter, grad_phi, phi, limiterNavierS, phi, epsilonSDWLS);
#endif
			} else {
				if(mpi_rank==0)
					cerr << "### IkeWithModels::calcEulerFlux1D_AD() => Wrong reconstruction type for scalars! ###" << endl;
				freeMemForCalcEulerFlux1D_AD(Apl, Ami, AplScal, AmiScal, FrhoScal, Scalar0, Scalar1, ScalCV0, ScalCV1, ScalConvTerm, nScal);
				throw(-1);
			}
		}
	}

	// Calculate flux
	int nofa_f = faocv_i[icvCenter];
	int nofa_l = faocv_i[icvCenter+1]-1;
	for (int foc = nofa_f; foc <= nofa_l; foc++) {
		int ifa = faocv_v[foc];
		// ===============================================================================================
		// cycle through internal faces, assembling flux to both sides
		// ===============================================================================================
		if((ifa >= nfa_b && ifa<nfa) || (ifa >= nfa_b2 && ifa < nfa_b2gg)) {
			int icv0 = cvofa[ifa][0];
			int icv1 = cvofa[ifa][1];
			assert( icv0 >= 0 && icv0 < ncv_gg);
			assert( icv1 >= 0 && icv1 < ncv_gg);

			int noc00, noc01, noc11, noc10;
			if (flagImplicit && icv0 < ncv && icv1 < ncv_g)
				getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

			// face unit normal and area...
			double nVec[3] = {0.0, 0.0, 0.0};
			double area = normVec3d(nVec, fa_normal[ifa]);

			// .............................................................................................
			// reconstruction of variables at faces: rho, u, T or P, scalars
			// .............................................................................................
			REALQ rho0 = rho[icv0];
			REALQ u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
			REALQ p0 = press[icv0];
			REALQ T0 = temp[icv0];
			REALQ h0 = enthalpy[icv0];
			REALQ gam0 = gamma[icv0];
			REALQ R0 = RoM[icv0];
			REALQS kineCV0 = 0.0;          // cell center
			REALQS kineFA0 = 0.0;          // cell face if second order

			REALQ rho1 = rho[icv1];
			REALQ u1[3] = {vel[icv1][0], vel[icv1][1], vel[icv1][2]};
			REALQ p1 = press[icv1];
			REALQ T1 = temp[icv1];
			REALQ h1 = enthalpy[icv1];
			REALQ gam1 = gamma[icv1];
			REALQ R1 = RoM[icv1];
			REALQS kineCV1 = 0.0;
			REALQS kineFA1 = 0.0;

			for (int iScal = 0; iScal < nScal; iScal++) {
#ifdef USE_MEM_SAVING_ADVAR
				ScalCV0[iScal] = Scalar0[iScal] = scalarTranspEqVector_AD[iScal].phi[icv0];
				ScalCV1[iScal] = Scalar1[iScal] = scalarTranspEqVector_AD[iScal].phi[icv1];
#else
				REALQS *phi = scalarTranspEqVector_AD[iScal].phi;
				ScalCV0[iScal] = Scalar0[iScal] = phi[icv0];
				ScalCV1[iScal] = Scalar1[iScal] = phi[icv1];
#endif
			}

			if (sndOrder == true) {
				double r0[3] = {0.0, 0.0, 0.0}, r1[3] = {0.0, 0.0, 0.0};
				vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);
				vecMinVec3d(r1, x_fa[ifa], x_cv[icv1]);

				// ----------------------------------------
				// left side
				// ----------------------------------------
				rho0 += vecDotVec3d_AD(r0, grad_rho[icv0]);
#ifdef temp_reconstruction
				T0 += vecDotVec3d_AD(r0, grad_temp[icv0]);
				if ((T0 <= 0.0) || (rho0 <= 0.0)) {
					T0 = temp[icv0];
					rho0 = rho[icv0];
					myCountReducedOrder++;
				}
#else
				p0 += vecDotVec3d_AD(r0, grad_p[icv0]);
				if ((p0 <= 0.0) || (rho0 <= 0.0)) {
					p0 = press[icv0];
					rho0 = rho[icv0];
					myCountReducedOrder++;
				}
#endif
				else {
					for (int i = 0; i < 3; i++)
						u0[i] += vecDotVec3d_AD(r0, grad_u[icv0][i]);
					for (int iScal = 0; iScal < nScal; iScal++) {
#ifdef USE_MEM_SAVING_ADVAR
						if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
							Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d_AD(r0, scalarTranspEqVector_AD[iScal].grad_phi[icv0])) / rho0;
						else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
							Scalar0[iScal] += vecDotVec3d_AD(r0, scalarTranspEqVector_AD[iScal].grad_phi[icv0]);
#else
						REALQS (*grad_phi)[3] = scalarTranspEqVector_AD[iScal].grad_phi;

						if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
							Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d_AD(r0, grad_phi[icv0])) / rho0;
						else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
							Scalar0[iScal] += vecDotVec3d_AD(r0, grad_phi[icv0]);
#endif
					}
				}

				// ----------------------------------------
				// right side
				// ----------------------------------------
				rho1 += vecDotVec3d_AD(r1, grad_rho[icv1]);
#ifdef temp_reconstruction
				T1 += vecDotVec3d_AD(r1, grad_temp[icv1]);
				if ((T1 <= 0.0) || (rho1 <= 0.0)) {
					T1 = temp[icv1];
					rho1 = rho[icv1];
					myCountReducedOrder++;
				}
#else
				p1 += vecDotVec3d_AD(r1, grad_p[icv1]);
				if ((p1 <= 0.0) || (rho1 <= 0.0)) {
					p1 = press[icv1];
					rho1 = rho[icv1];
					myCountReducedOrder++;
				}
#endif
				else {
					for (int i = 0; i < 3; i++)
						u1[i] += vecDotVec3d_AD(r1, grad_u[icv1][i]);
					for (int iScal = 0; iScal < nScal; iScal++) {
#ifdef USE_MEM_SAVING_ADVAR
						if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
							Scalar1[iScal] = (rho[icv1] * Scalar1[iScal] + vecDotVec3d_AD(r1, scalarTranspEqVector_AD[iScal].grad_phi[icv1])) / rho1;    // CHANGE FOR SCALAR AD

						else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
							Scalar1[iScal] += vecDotVec3d_AD(r1, scalarTranspEqVector_AD[iScal].grad_phi[icv1]);
#else
						REALQS (*grad_phi)[3] = scalarTranspEqVector_AD[iScal].grad_phi;

						if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
							Scalar1[iScal] = (rho[icv1] * Scalar1[iScal] + vecDotVec3d_AD(r1, grad_phi[icv1])) / rho1;    // CHANGE FOR SCALAR AD

						else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
							Scalar1[iScal] += vecDotVec3d_AD(r1, grad_phi[icv1]);
#endif
					}
				}

				// .............................................................................................
				// calculation of other variables at faces: p/T, h, R, gam
				// .............................................................................................
#ifdef temp_reconstruction
				calcThermoProp_T_AD(p0, h0, R0, gam0, rho0, T0, Scalar0, nScal);
				calcThermoProp_T_AD(p1, h1, R1, gam1, rho1, T1, Scalar1, nScal);
#else
				calcThermoProp_p_AD(T0, h0, R0, gam0, rho0, p0, Scalar0, nScal);
				calcThermoProp_p_AD(T1, h1, R1, gam1, rho1, p1, Scalar1, nScal);
#endif
			}


			if (kine_Index > -1) {  // save kine if defined
				kineCV0 = ScalCV0[kine_Index];         // cell center left for implicit side (Jacobi is computed with cell center)
				kineFA0 = Scalar0[kine_Index];         // cell face left, if second order
				kineCV1 = ScalCV1[kine_Index];         // cell center left for implicit side (Jacobi is computed with cell center)
				kineFA1 = Scalar1[kine_Index];         // cell face right,
			}

			// .............................................................................................
			// calculation of Euler Flux explicit using HLLC
			// .............................................................................................
			calcEulerFlux_HLLC_AD(Frho, Frhou, FrhoE, FrhoScal,
					rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
					rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
					area, nVec, nScal, 0.0);

			if(icvCenter == icv0) {
				rhs_rho -= Frho;
				for (int i = 0; i < 3; i++)
					rhs_rhou[i] -= Frhou[i];
				rhs_rhoE -= FrhoE;
				for (int iScal = 0; iScal < nScal; iScal++)
					rhs_rhoScal[iScal] -= ScalConvTerm[iScal] * FrhoScal[iScal];
			} else if(icvCenter == icv1) {
				rhs_rho += Frho;
				for (int i = 0; i < 3; i++)
					rhs_rhou[i] += Frhou[i];
				rhs_rhoE += FrhoE;
				for (int iScal = 0; iScal < nScal; iScal++)
					rhs_rhoScal[iScal] += ScalConvTerm[iScal] * FrhoScal[iScal];
			} else {
				printf("Error in IkeWithModels_AD::calcEulerFlux1D_AD() => icvCenter(%d) is not equal to icv0(%d) nor icv1(%d) \n",icvCenter, icv0, icv1);
				freeMemForCalcEulerFlux1D_AD(Apl, Ami, AplScal, AmiScal, FrhoScal, Scalar0, Scalar1, ScalCV0, ScalCV1, ScalConvTerm, nScal);
				throw(-1);
			}

			if(isNaN(rhs_rho.value()) || isNaN(rhs_rhou[0].value()) || isNaN(rhs_rhou[1].value()) || isNaN(rhs_rhou[2].value()) || isNaN(rhs_rhoE.value())) {
				printf("Error in IkeWithModels_AD::calcEulerFlux1D_AD: NaN detected in flux at internal ifa=%d(%.2e,%.2e,%.2e)\n", ifa, x_cv[ifa][0], x_cv[ifa][1], x_cv[ifa][2]);
				printf("      Details: mpi_rank=%d, icvCenter=%d, Frho=%.3e, Frhou=(%.3e,%.3e,%.3e), FrhoE=%3e\n", mpi_rank, icvCenter, Frho.value(), Frhou[0].value(),Frhou[1].value(),Frhou[2].value(), FrhoE.value());
				printf("               rho0=%.3e, u0=(%.3e, %.3e, , %.3e), p0=%.3e, T0=%.3e, h0=%.3e, R0=%.3e, gam0=%.3e, kineFA0=%.3e\n", rho0.value(), u0[0].value(), u0[1].value(), u0[2].value(), p0.value(), T0.value(), h0.value(), R0.value(), gam0.value(), kineFA0.value());
				printf("               rho1=%.3e, u1=(%.3e, %.3e, , %.3e), p1=%.3e, T1=%.3e, h1=%.3e, R1=%.3e, gam1=%.3e, kineFA1=%.3e\n", rho1.value(), u1[0].value(), u1[1].value(), u1[2].value(), p1.value(), T1.value(), h1.value(), R1.value(), gam1.value(), kineFA1.value());
				freeMemForCalcEulerFlux1D_AD(Apl, Ami, AplScal, AmiScal, FrhoScal, Scalar0, Scalar1, ScalCV0, ScalCV1, ScalConvTerm, nScal);
				throw(IKEWITHMODELS_ERROR_CODE);
			}

			// .............................................................................................
			// calculate implicit matrix using HLLC
			// .............................................................................................
			if (flagImplicit && icv0 < ncv && icv1 < ncv_g) {
				calcEulerFluxMatrices_HLLC_AD(Apl, Ami, AplScal, AmiScal,
						rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, kineCV0,
						rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, kineCV1,
						area, nVec, nScal, 0.0);

				if (icv0 < ncv) { // if icv0 is internal...
					for (int i = 0; i < 5; i++)
						for (int j = 0; j < 5; j++) {
							A[noc00][i][j] += Apl[j][i];
							A[noc01][i][j] -= Apl[j][i];
						}
				}

				if (icv1 < ncv) { // if icv1 is internal...
					for (int i = 0; i < 5; i++)
						for (int j = 0; j < 5; j++) {
							A[noc11][i][j] -= Ami[j][i];
							A[noc10][i][j] += Ami[j][i];
						}
				}

				for (int iScal = 0; iScal < nScal; iScal++)
					for (int i = 0; i <= 5; i++) {
						if (icv0 < ncv) {
							AScal[iScal][i][noc00] += ScalConvTerm[iScal] * AplScal[iScal][i];
							AScal[iScal][i][noc01] -= ScalConvTerm[iScal] * AplScal[iScal][i];
						}

						if (icv1 < ncv) {
							AScal[iScal][i][noc11] -= ScalConvTerm[iScal] * AmiScal[iScal][i];
							AScal[iScal][i][noc10] += ScalConvTerm[iScal] * AmiScal[iScal][i];
						}
					}
			}
		}
		// ===============================================================================================
		// cycle through boundary faces, assembling flux
		// ===============================================================================================
		else {
			FaZone* zone = zofa[ifa];
			if (zone->getKind() == FA_ZONE_BOUNDARY) {
				Param *param;

				if (getParam(param, zone->getName())) {
					// .............................................................................................
					// SYMMETRY BOUNDARY CONDITION OR WALL BOUNDARY, ATTENTION: ONLY FIRST ORDER!!!
					// .............................................................................................
					if(getBoundaryType(ifa, param) == SYMMETRY || getBoundaryType(ifa, param) == WALL) {
						int icv0 = cvofa[ifa][0];
						int icv1 = cvofa[ifa][1];

						double nVec[3] = {0.0, 0.0, 0.0};
						double area = normVec3d(nVec, fa_normal[ifa]);

						for (int iScal = 0; iScal < nScal; iScal++)
							Scalar0[iScal] = scalarTranspEqVector_AD[iScal].phi[icv1];

						REALQS kineFA = 0.0;
						if (kine_Index > -1)
							kineFA = scalarTranspEqVector_AD[kine_Index].phi[icv1];

						calcEulerFlux_HLLC_AD(Frho, Frhou, FrhoE, FrhoScal,
								rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
								rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
								area, nVec, nScal, 0.0);

						rhs_rho -= Frho;
						for (int i = 0; i < 3; i++)
							rhs_rhou[i] -= Frhou[i];
						rhs_rhoE -= FrhoE;
						for (int iScal = 0; iScal < nScal; iScal++)
							rhs_rhoScal[iScal] -= ScalConvTerm[iScal] * FrhoScal[iScal];

						if(isNaN(rhs_rho.value()) || isNaN(rhs_rhou[0].value()) || isNaN(rhs_rhou[1].value()) || isNaN(rhs_rhou[2].value()) || isNaN(rhs_rhoE.value())) {
							printf("Error in IkeWithModels_AD::calcEulerFlux1D_AD: NaN detected in flux at boundary (SYMM or WALL) ifa=%d(%.2e,%.2e,%.2e)\n", ifa, x_cv[ifa][0], x_cv[ifa][1], x_cv[ifa][2]);
							printf("      Details: mpi_rank=%d, icvCenter=%d, Frho=%.3e, Frhou=(%.3e,%.3e,%.3e), FrhoE=%3e\n", mpi_rank, icvCenter, Frho.value(), Frhou[0].value(),Frhou[1].value(),Frhou[2].value(), FrhoE.value());
							printf("               rho=%.3e, vel=(%.3e, %.3e, %.3e), press=%.3e, temp=%.3e, enthalpy=%.3e, RoM=%.3e, gamma=%.3e, kineFA=%.3e\n", rho[icv1].value(), vel[icv1][0].value(),vel[icv1][1].value(),vel[icv1][2].value(), press[icv1].value(), temp[icv1].value(), enthalpy[icv1].value(), RoM[icv1].value(), gamma[icv1].value(), kineFA.value());
							freeMemForCalcEulerFlux1D_AD(Apl, Ami, AplScal, AmiScal, FrhoScal, Scalar0, Scalar1, ScalCV0, ScalCV1, ScalConvTerm, nScal);
							throw(IKEWITHMODELS_ERROR_CODE);
						}

						if (flagImplicit && icv0 < ncv) {
							calcEulerFluxMatrices_HLLC_AD(Apl, NULL, AplScal, AmiScal,
									rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
									rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
									area, nVec, nScal, 0.0);

							int noc00 = nbocv_i[icv0]; // icv0's diagonal
							for (int i = 0; i < 5; i++)
								for (int j = 0; j < 5; j++)
									A[noc00][i][j] += Apl[j][i];

							for (int iScal = 0; iScal < nScal; iScal++)
								for (int i = 0; i <= 5; i++)
									AScal[iScal][i][noc00] += ScalConvTerm[iScal] * AplScal[iScal][i];
						}
					}
					// .............................................................................................
					// OTHER BOUNDARY CONDITIONS (HOOK, DIRICHLET(CBC), NEUMANN, ...)
					// .............................................................................................
					else {
						int icv0 = cvofa[ifa][0];
						int icv1 = cvofa[ifa][1];

						double nVec[3] = {0.0, 0.0, 0.0};
						double area = normVec3d(nVec, fa_normal[ifa]);

						REALQ rho0 = rho[icv0];
						REALQ u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
						REALQ p0 = press[icv0];
						REALQ T0 = temp[icv0];
						REALQ h0 = enthalpy[icv0];
						REALQ gam0 = gamma[icv0];
						REALQ R0 = RoM[icv0];
						REALQS kineCV0 = 0.0;           // cell center
						REALQS kineFA0 = 0.0;           // cell face

						REALQS kineFA1 = 0.0;

						for (int iScal = 0; iScal < nScal; iScal++) {
							ScalCV0[iScal] = Scalar0[iScal] = scalarTranspEqVector_AD[iScal].phi[icv0];
							Scalar1[iScal] = scalarTranspEqVector_AD[iScal].phi[icv1];
						}

						if (sndOrder == true) {
							double r0[3] = {0.0, 0.0, 0.0};
							vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);

							// left side
							rho0 += vecDotVec3d_AD(r0, grad_rho[icv0]);
#ifdef temp_reconstruction
							T0 += vecDotVec3d_AD(r0, grad_temp[icv0]);
							if ((T0 <= 0.0) || (rho0 <= 0.0)) {
								T0 = temp[icv0];
								rho0 = rho[icv0];
								myCountReducedOrder++;
							}
#else
							p0 += vecDotVec3d_AD(r0, grad_p[icv0]);
							if ((p0 <= 0.0) || (rho0 <= 0.0)) {
								p0 = press[icv0];
								rho0 = rho[icv0];
								myCountReducedOrder++;
							}
#endif
							else {
								for (int i = 0; i < 3; i++)
									u0[i] += vecDotVec3d_AD(r0, grad_u[icv0][i]);
								for (int iScal = 0; iScal < nScal; iScal++) {
#ifdef USE_MEM_SAVING_ADVAR
									if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
										Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d_AD(r0, scalarTranspEqVector_AD[iScal].grad_phi[icv0])) / rho0;    // CHANGE FOR SCALAR AD
									else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
										Scalar0[iScal] += vecDotVec3d_AD(r0, scalarTranspEqVector_AD[iScal].grad_phi[icv0]);
#else
									REALQS (*grad_phi)[3] = scalarTranspEqVector_AD[iScal].grad_phi;
									if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
										Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d_AD(r0, grad_phi[icv0])) / rho0;    // CHANGE FOR SCALAR AD
									else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
										Scalar0[iScal] += vecDotVec3d_AD(r0, grad_phi[icv0]);
#endif
								}
							}

							// calculation of other variables at faces: p/T, h, R, gam
#ifdef temp_reconstruction
							calcThermoProp_T_AD(p0, h0, R0, gam0, rho0, T0, Scalar0, nScal);
#else
							calcThermoProp_p_AD(T0, h0, R0, gam0, rho0, p0, Scalar0, nScal);
#endif
						}

						if (kine_Index > -1) {  // save kine if defined
							kineCV0 = ScalCV0[kine_Index];          // cell center
							kineFA0 = Scalar0[kine_Index];                                 // cell face
							kineFA1 = Scalar1[kine_Index];
						}

						calcEulerFlux_HLLC_AD(Frho, Frhou, FrhoE, FrhoScal,
								rho0,      u0,        p0,          T0,         h0,             R0,        gam0,        Scalar0, kineFA0,
								rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar1, kineFA1,
								area, nVec, nScal, 0.0);

						rhs_rho -= Frho;
						for (int i = 0; i < 3; i++)
							rhs_rhou[i] -= Frhou[i];
						rhs_rhoE -= FrhoE;
						for (int iScal = 0; iScal < nScal; iScal++)
							rhs_rhoScal[iScal] -= ScalConvTerm[iScal] * FrhoScal[iScal];

						if(isNaN(rhs_rho.value()) || isNaN(rhs_rhou[0].value()) || isNaN(rhs_rhou[1].value()) || isNaN(rhs_rhou[2].value()) || isNaN(rhs_rhoE.value())) {
							printf("Error in IkeWithModels_AD::calcEulerFlux1D_AD: NaN detected in flux at boundary (HOOK, DIRI(CBC), NEUM) ifa=%d(%.2e,%.2e,%.2e)\n", ifa, x_cv[ifa][0], x_cv[ifa][1], x_cv[ifa][2]);
							printf("      Details: mpi_rank=%d, icvCenter=%d, Frho=%.3e, Frhou=(%.3e,%.3e,%.3e), FrhoE=%3e\n", mpi_rank, icvCenter, Frho.value(), Frhou[0].value(),Frhou[1].value(),Frhou[2].value(), FrhoE.value());
							printf("               rho0=%.3e, u0=(%.3e, %.3e, , %.3e), p0=%.3e, T0=%.3e, h0=%.3e, R0=%.3e, gam0=%.3e, kineFA0=%.3e\n", rho0.value(), u0[0].value(), u0[1].value(), u0[2].value(), p0.value(), T0.value(), h0.value(), R0.value(), gam0.value(), kineFA0.value());
							printf("               rho=%.3e, vel=(%.3e, %.3e, %.3e), press=%.3e, temp=%.3e, enthalpy=%.3e, RoM=%.3e, gamma=%.3e, kineFA=%.3e\n", rho[icv1].value(), vel[icv1][0].value(),vel[icv1][1].value(),vel[icv1][2].value(), press[icv1].value(), temp[icv1].value(), enthalpy[icv1].value(), RoM[icv1].value(), gamma[icv1].value(), kineFA1.value());
							freeMemForCalcEulerFlux1D_AD(Apl, Ami, AplScal, AmiScal, FrhoScal, Scalar0, Scalar1, ScalCV0, ScalCV1, ScalConvTerm, nScal);
							throw(IKEWITHMODELS_ERROR_CODE);
						}

						if (flagImplicit && icv0 < ncv) {
							calcEulerFluxMatrices_HLLC_AD(Apl, NULL, AplScal, AmiScal,
									rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, kineCV0,
									rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar1, kineFA1,
									area, nVec, nScal, 0.0);

							int noc00 = nbocv_i[icv0]; // icv0's diagonal
							for (int i = 0; i < 5; i++)
								for (int j = 0; j < 5; j++)
									A[noc00][i][j] += Apl[j][i];

							for (int iScal = 0; iScal < nScal; iScal++)
								for (int i = 0; i <= 5; i++)
									AScal[iScal][i][noc00] += ScalConvTerm[iScal] * AplScal[iScal][i];
						}
					}
				}
			}
		}
	}

	freeMemForCalcEulerFlux1D_AD(Apl, Ami, AplScal, AmiScal, FrhoScal, Scalar0, Scalar1, ScalCV0, ScalCV1, ScalConvTerm, nScal);
	return myCountReducedOrder;
}

/*
 * Method: freeMemForCalcEulerFlux1D_AD
 * ------------------------------------
 * Delete the memory allocated in the calcEulerFlux1D_AD() function
 * (The arrays will NOT be re-initialized by NULL)
 */
void IkeWithModels_AD::freeMemForCalcEulerFlux1D_AD(double (*Apl)[5], double (*Ami)[5], double (*AplScal)[6], double (*AmiScal)[6],
	REALQS *FrhoScal, REALQS *Scalar0, REALQS *Scalar1, REALQS *ScalCV0, REALQS *ScalCV1, double *ScalConvTerm, const int nScal) {
	if (Apl != NULL)  delete [] Apl;
	if (Ami != NULL)  delete [] Ami;

	if (AplScal != NULL) delete [] AplScal;
	if (AmiScal != NULL) delete [] AmiScal;

	if (nScal > 0) delete [] FrhoScal;
	if (nScal > 0) delete [] Scalar0;
	if (nScal > 0) delete [] Scalar1;
	if (nScal > 0) delete [] ScalCV0;
	if (nScal > 0) delete [] ScalCV1;
	if (nScal > 0) delete [] ScalConvTerm;
}

/*
 * Method: calcViscousFluxNS1D_AD
 * ------------------------------
 * Calculate viscous flux in the so-called "1D style"
 * Original code =calcViscousFluxNS_AD() in JoeWithModels.cpp
 */
void IkeWithModels_AD::calcViscousFluxNS1D_AD(const int icvCenter, REALA &rhs_rho, REALA rhs_rhou[3], REALA &rhs_rhoE, REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE, double (*A)[5][5], int flagImplicit) {
	double (*A0)[5];
	double (*A1)[5];

	if (flagImplicit) {
		A0 = new double[5][5];
		A1 = new double[5][5];

		for (int i = 0; i < 5; i++)
			for (int j = 0; j < 5; j++)
				A0[i][j] = A1[i][j] = 0.0;
	} else {
		A0 = A1 = NULL;
	}

	REALQ Frhou[3] = {0.0, 0.0, 0.0}, FrhoE = 0.0;

	// save the index of kine if defined
	int kine_index = -1;
	for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
//		if (scalarTranspEqVector[iScal].getName() == "kine")
		if (strcmp(scalarTranspEqVector[iScal].getName(),"kine")==0)
			kine_index = iScal;

	// ====================================================================
	//        compute gradients, with boundary values
	// ====================================================================
	calcCv2Grad1D_AD(icvCenter, grad_enthalpy, enthalpy, limiterNavierS, enthalpy, epsilonSDWLS);

	// Calculate flux
	int nofa_f = faocv_i[icvCenter];
	int nofa_l = faocv_i[icvCenter+1]-1;
	for (int foc = nofa_f; foc <= nofa_l; foc++) {
		int ifa = faocv_v[foc];

		// ====================================================================
		// cycle through internal faces, assembling flux and matrix
		// ====================================================================
		if((ifa >= nfa_b && ifa<nfa) || (ifa >= nfa_b2 && ifa < nfa_b2gg)) {
			int icv0 = cvofa[ifa][0];
			int icv1 = cvofa[ifa][1];

			int noc00, noc01, noc11, noc10;
			if (flagImplicit && icv0 < ncv && icv1 < ncv_g)
				getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

			// face unit normal and area...
			double nVec[3] = {0.0, 0.0, 0.0};
			double area = normVec3d(nVec, fa_normal[ifa]);
			double sVec[3] = {0.0, 0.0, 0.0};
			vecMinVec3d(sVec, x_cv[icv1], x_cv[icv0]);
			double smag = normVec3d(sVec);

			double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
			vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
			vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
			double w0 = sqrt(vecDotVec3d(dx0, dx0));
			double w1 = sqrt(vecDotVec3d(dx1, dx1));
			double ws = w0 + w1;
			w0 /= ws;
			w1 /= ws;
			REALQ uAux_fa[3] = { w1*vel[icv0][0]+ w0*vel[icv1][0],
					w1*vel[icv0][1]+ w0*vel[icv1][1],
					w1*vel[icv0][2]+ w0*vel[icv1][2]};

			// kinetic energy, if defined
			REALQS kine0 = 0.0;
			REALQS kine1 = 0.0;
			REALQS kine_fa = 0.0;
			if (kine_index > -1) {
				kine0 = scalarTranspEqVector_AD[kine_index].phi[icv0];
				kine1 = scalarTranspEqVector_AD[kine_index].phi[icv1];
				kine_fa = w1*kine0 + w0*kine1;
			}

			// calculate viscous flux
#ifdef USE_ARTIF_VISC_WITH_MEM_SAVING
			// If the user wants to have artificial bulk viscosity, pass the information to the addViscFlux() method
			// (If the user wants to use artificial viscosity (not bulk only), it was already reflected in mul_fa)
			REALQS artifBulkViscosity_AD = 0.0;
			if(turnOnArtifVisc && artifVisc_bulkViscOnly)
				artifBulkViscosity_AD = w1*artifViscMag_AD[icv0] + w0*artifViscMag_AD[icv1];

			addViscFlux_AD(Frhou, FrhoE, A0, A1,
					rho[icv0], vel[icv0], grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0], gamma[icv0], kine0,
					rho[icv1], vel[icv1], grad_u[icv1], enthalpy[icv1], grad_enthalpy[icv1], temp[icv1], RoM[icv1], gamma[icv1], kine1,
					mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kine_fa, uAux_fa,
					area, nVec, smag, sVec, artifBulkViscosity_AD);
#else
			addViscFlux_AD(Frhou, FrhoE, A0, A1,
					rho[icv0], vel[icv0], grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0], gamma[icv0], kine0,
					rho[icv1], vel[icv1], grad_u[icv1], enthalpy[icv1], grad_enthalpy[icv1], temp[icv1], RoM[icv1], gamma[icv1], kine1,
					mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kine_fa, uAux_fa,
					area, nVec, smag, sVec);
#endif

			if (flagImplicit && icv0 < ncv && icv1 < ncv_g) {
				for (int i=0; i<5; i++)
					for (int j=0; j<5; j++) {
						A[noc00][i][j] -= A0[j][i];
						A[noc01][i][j] += A0[j][i];
					}

				if (icv1 < ncv) { // if icv1 is internal...
					for (int i=0; i<5; i++)
						for (int j=0; j<5; j++) {
							A[noc11][i][j] += A1[j][i];
							A[noc10][i][j] -= A1[j][i];
						}
				}
			}

			if(icvCenter == icv0) {
				for (int i = 0; i < 3; i++)
					rhs_rhou[i] -= Frhou[i];
				rhs_rhoE -= FrhoE;
			} else if (icv1 == icv1) {
				for (int i = 0; i < 3; i++)
					rhs_rhou[i] += Frhou[i];
				rhs_rhoE += FrhoE;
			}
		}
		// ====================================================================
		// cycle through boundary faces, assembling flux and matrix
		// ====================================================================
		else {
			FaZone* zone = zofa[ifa];
			// .............................................................................................
			// SYMMETRY BOUNDARY CONDITION
			// .............................................................................................
			if(getBoundaryType(ifa, zone) == SYMMETRY) {
				if (kine_index > -1) {
					int icv0 = cvofa[ifa][0];
					int icv1 = cvofa[ifa][1];
					assert( icv0 >= 0 );
					double nVec[3] = {0.0, 0.0, 0.0};
					double area = normVec3d(nVec, fa_normal[ifa]);

#ifdef USE_MEM_SAVING_ADVAR
					adouble kine_fa  = scalarTranspEqVector_AD[kine_index].phi[icv1];
#else
					adouble *phi = scalarTranspEqVector_AD[kine_index].phi;
					adouble kine_fa  = phi[icv1];
#endif

					adouble tmp = 1.0/3.0*(rho[icv0] + rho[icv1])*kine_fa;

					for (int i = 0; i < 3; i++)
						rhs_rhou[i] -= tmp*fa_normal[ifa][i];
				}
			}
			// .............................................................................................
			// WALL BOUNDARY CONDITION
			// .............................................................................................
			else if(getBoundaryType(ifa, zone) == WALL) {
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];

				REALX nVec[3] = {0.0, 0.0, 0.0};
				REALX area = normVec3d(nVec, fa_normal[ifa]);
				REALX sVec[3] = {0.0, 0.0, 0.0};
				vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
				REALX smag = fabs(vecDotVec3d(sVec, nVec));  // project sVec to wall face normal

				REALQS kine0 = 0.0;
				REALQS kine1 = 0.0;
				REALQS kine_fa = 0.0;
				if (kine_index > -1)
					kine0 = scalarTranspEqVector_AD[kine_index].phi[icv0];

				// calculate viscous flux
				addViscFlux_AD(Frhou, FrhoE, A0, NULL,
						rho[icv0],    vel[icv0],    grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0],    gamma[icv0],  kine0,
						rho[icv1],    vel[icv1],    grad_u[icv0], enthalpy[icv1], grad_enthalpy[icv0], temp[icv1], RoM[icv1],    gamma[icv1],  kine1,
						mul_fa[ifa], 0.0, lamOcp_fa[ifa], kine_fa, vel[icv1],
						area, nVec, smag, nVec);  /* <- use nVec here instead of sVec, to avoid inaccurate correction*/

				if (flagImplicit && icv0<ncv) {
					int noc00 = nbocv_i[icv0]; // icv0's diagonal

					for (int i=0; i<5; i++)
						for (int j=0; j<5; j++)
							A[noc00][i][j] -= A0[j][i];
				}

				if(icvCenter == icv0) {
					for (int i = 0; i < 3; i++)
						rhs_rhou[i] -= Frhou[i];
					rhs_rhoE -= FrhoE;
				}
			}
			// .............................................................................................
			// OTHER BOUNDARY CONDITIONS
			// .............................................................................................
			else {
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];

				REALX nVec[3] = {0.0, 0.0, 0.0};
				REALX area = normVec3d(nVec, fa_normal[ifa]);
				REALX sVec[3] = {0.0, 0.0, 0.0};
				vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
				REALX smag = normVec3d(sVec);

				REALQS kine0 = 0.0;
				REALQS kine1 = 0.0;
				REALQS kine_fa = 0.0;
				if (kine_index > -1) {
#ifdef USE_MEM_SAVING_ADVAR
					kine_fa = scalarTranspEqVector_AD[kine_index].phi[icv1];
#else
					kine0 = scalarTranspEqVector_AD[kine_index].phi[icv0];
					kine1 = scalarTranspEqVector_AD[kine_index].phi[icv1];
					kine_fa = kine1;
#endif
				}

				// calculate viscous flux
				addViscFlux_AD(Frhou, FrhoE, A0, NULL,
						rho[icv0],    vel[icv0],    grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0],    gamma[icv0],  kine0,
						rho[icv1],    vel[icv1],    grad_u[icv0], enthalpy[icv1], grad_enthalpy[icv0], temp[icv1], RoM[icv1],    gamma[icv1],  kine1,
						mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kine_fa, vel[icv1],
						area, nVec, smag, sVec);

				if (flagImplicit && icv0<ncv) {
					int noc00 = nbocv_i[icv0]; // icv0's diagonal

					for (int i=0; i<5; i++)
						for (int j=0; j<5; j++)
							A[noc00][i][j] -= A0[j][i];
				}

				if(icvCenter == icv0) {
					for (int i = 0; i < 3; i++)
						rhs_rhou[i] -= Frhou[i];
					rhs_rhoE -= FrhoE;
				}
			}
		}
	}

	if (A0  != NULL)  delete [] A0;
	if (A1  != NULL)  delete [] A1;
}

/*
 * Method: calcViscousFluxScalar_new1D_AD
 * --------------------------------------
 * Original code: calcViscousFluxScalar_new_AD in ScalarsAD.cpp
 */
void IkeWithModels_AD::calcViscousFluxScalar_new1D_AD(const int icvCenter, const int iScal, adouble &rhs_rhoScal, double *AScal, ScalarTranspEq &transpScal, ScalarTranspEq_AD &transpScal_AD, int flagImplicit) {
	if (!transpScal.diffTerm) { // Note: transScal_AD doesn't have "diffTerm"
		return;
	}

	string scalName(transpScal.getName());
	double phiBCval = 0.0, phiBCflux = 0.0;

#ifdef USE_MEM_SAVING_ADVAR
	// =============================================================================================
	// compute gradients, with boundary values
	// =============================================================================================
	calcCv2Grad1D_AD(icvCenter, transpScal_AD.grad_phi, transpScal_AD.phi, limiterNavierS, transpScal_AD.phi, epsilonSDWLS);

	// =============================================================================================
	// compute user-defined diffusivity
	// =============================================================================================
	vector<int> faInternal;
	vector<int> faBoundary;
	//find2LayerFaceIndices(icvCenter, faInternal, faBoundary);
	find1LayerFaceIndices(icvCenter, faInternal, faBoundary);
	// internal faces
	for(size_t index=0; index<faInternal.size(); ++index) {
		int ifa = faInternal[index];
		FaZone* zone = zofa[ifa];

		diffusivityHookScalarRansTurb1D_AD(ifa, zone, scalName);
		diffusivityHookScalarRansComb1D_AD(ifa, zone, scalName);
	}
	for(size_t index=0; index<faBoundary.size(); ++index) {
		int ifa = faBoundary[index];
		FaZone* zone = zofa[ifa];
		diffusivityHookScalarRansTurb1D_AD(ifa, zone, scalName);
		diffusivityHookScalarRansComb1D_AD(ifa, zone, scalName);
	}

	// Calculate flux
	int nofa_f = faocv_i[icvCenter];
	int nofa_l = faocv_i[icvCenter+1]-1;
	for (int nofa = nofa_f; nofa <= nofa_l; nofa++) {
		int ifa = faocv_v[nofa];

		// =============================================================================================
		// cycle trough internal faces and compute viscous flux and implicit matrix
		// =============================================================================================
		if((ifa >= nfa_b && ifa<nfa) || (ifa >= nfa_b2 && ifa < nfa_b2gg)) {
			int icv0 = cvofa[ifa][0];
			int icv1 = cvofa[ifa][1];
			int noc00, noc01, noc11, noc10;

			if (flagImplicit && icv0 < ncv && icv1 < ncv_g)
				getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

			double n[3] = {0.0, 0.0, 0.0}, s[3] = {0.0, 0.0, 0.0};
			double nmag = normVec3d(n, fa_normal[ifa]);
			vecMinVec3d(s, x_cv[icv1], x_cv[icv0]);
			double smag = normVec3d(s);
			double alpha = vecDotVec3d(n, s);
			assert((alpha > 0.0) && (alpha < 1.0+1.0E-12));   // alpha should now contain s dot n...

			double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
			vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
			vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);

			adouble viscFlux = transpScal_AD.diff[ifa] * nmag * (alpha * (transpScal_AD.phi[icv1] - transpScal_AD.phi[icv0]) / smag
					+ 0.5 * (transpScal_AD.grad_phi[icv0][0] + transpScal_AD.grad_phi[icv1][0]) * (n[0] - alpha * s[0])
					+ 0.5 * (transpScal_AD.grad_phi[icv0][1] + transpScal_AD.grad_phi[icv1][1]) * (n[1] - alpha * s[1])
					+ 0.5 * (transpScal_AD.grad_phi[icv0][2] + transpScal_AD.grad_phi[icv1][2]) * (n[2] - alpha * s[2]));

			if(icvCenter == icv0)
				rhs_rhoScal += viscFlux;
			else if(icvCenter == icv1)
				rhs_rhoScal -= viscFlux;
			else {
				printf("ERROR in IkeWithModels_AD::calcViscousFluxScalar_new1D_AD() => icvCenter=%d is not equal to icv0=%d nor icv1=%d \n", icvCenter, icv0, icv1);
			}

			if (flagImplicit && icv0 < ncv && icv1 < ncv_g) {
				AScal[noc00] +=   transpScal_AD.diff[ifa].value() * nmag * alpha / smag / rho_AD[icv0].value();
				AScal[noc01] += - transpScal_AD.diff[ifa].value() * nmag * alpha / smag / rho_AD[icv0].value();

				if (icv1 < ncv) {
					AScal[noc11] +=   transpScal_AD.diff[ifa].value() * nmag * alpha / smag / rho_AD[icv1].value();
					AScal[noc10] += - transpScal_AD.diff[ifa].value() * nmag * alpha / smag / rho_AD[icv1].value();
				}
			}
		}
		// =============================================================================================
		// cycle trough boundary faces and compute viscous flux and implicit matrix
		// =============================================================================================
		else {
			FaZone* zone = zofa[ifa];
			// .............................................................................................
			// HOOK BOUNDARY CONDITION
			// .............................................................................................
			if(getBoundaryTypeScalar(ifa, zone, scalName) == HOOK_SCALAR)  {
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];
				int noc00 = nbocv_i[icv0];

				double n[3] = {0.0, 0.0, 0.0}, s_half[3] = {0.0, 0.0, 0.0};
				double nmag = normVec3d(n, fa_normal[ifa]);
				vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
				double smag_half = fabs(vecDotVec3d(s_half, n));   //normVec3d(s_half);
				//double alpha = vecDotVec3d(n, s_half);           // HACK alfa should be always 1 at the boundary (provide such a grid!)

				adouble viscFlux = transpScal_AD.diff[ifa] * nmag * (transpScal_AD.phi[icv1] - transpScal_AD.phi[icv0]) / smag_half;

				if(icvCenter == icv0)
					rhs_rhoScal += viscFlux;

				if (flagImplicit && icv0<ncv)
					AScal[noc00] += transpScal_AD.diff[ifa].value() * nmag / smag_half / rho_AD[icv0].value();
			}
			// .............................................................................................
			// DIRICHLET BOUNDARY CONDITION
			// .............................................................................................
			else if (getBoundaryTypeScalar(ifa, zone, scalName) == DIRICHLET_SCALAR){
			    Param *p;
			    string fName = zone->getName();
			    if (getParam(p, fName + "." + scalName))
			    	phiBCval = p->getDouble(1);
			    else {
			    	cout<<"Error in IkeWithModels_AD::calcViscousFluxScalar_new1D_AD() => Cannot find Param for "<<fName+"."+scalName<<" at ifa="<<ifa<<endl;
			    	throw(-1);
			    }

				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];
				int noc00 = nbocv_i[icv0];

				double n[3] = {0.0, 0.0, 0.0}, s_half[3] = {0.0, 0.0, 0.0};
				double nmag = normVec3d(n, fa_normal[ifa]);
				vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
				double smag_half = fabs(vecDotVec3d(s_half, n));   //normVec3d(s_half);
				//double alpha = vecDotVec3d(n, s_half);           // HACK alfa should be always 1 at the boundary (provide such a grid!)

				adouble viscFlux = transpScal_AD.diff[ifa] * nmag * (transpScal_AD.phi[icv1] - transpScal_AD.phi[icv0]) / smag_half;                // TODO: phi_bfa or phiBCval??

				if(icvCenter == icv0)
					rhs_rhoScal += viscFlux;

				if (flagImplicit && icv0<ncv)
					AScal[noc00] += transpScal_AD.diff[ifa].value() * nmag / smag_half / rho_AD[icv0].value();
			}
			// .............................................................................................
			// FLUX BOUNDARY CONDITION: NOT IMPLEMENTED YET
			// .............................................................................................
			else if (getBoundaryTypeScalar(ifa, zone, scalName) == FLUX_SCALAR) {
//				Param *p;
//				string fName = zone->getName();
//				if (getParam(p, fName + "." + scalName))
//					phiBCflux = p->getDouble(1);
//				else {
//					cout<<"Error in IkeWithModels_AD::calcViscousFluxScalar_new1D_AD() => Cannot find Param for "<<fName+"."+scalName+".flux"<<" at ifa="<<ifa<<endl;
//					throw(-1);
//				}
				cerr << "scalarZoneIsFlux not implemented yet!" << endl;
				throw(-1);
			}
			// .............................................................................................
			// OTHER BOUNDARY CONDITIONS: NO VISCOUS FLUX
			// .............................................................................................
			else {
				// note: SYMMETRY
			}
		}
	}
#else
	adouble *phi = transpScal_AD.phi;
	adouble (*grad_phi)[3] = transpScal_AD.grad_phi;
	adouble *diff = transpScal_AD.diff;

	// =============================================================================================
	// compute gradients, with boundary values
	// =============================================================================================
	calcCv2Grad1D_AD(icvCenter, grad_phi, phi, limiterNavierS, phi, epsilonSDWLS);

	// =============================================================================================
	// compute user-defined diffusivity
	// =============================================================================================
	vector<int> faInternal;
	vector<int> faBoundary;
	//find2LayerFaceIndices(icvCenter, faInternal, faBoundary);
	find1LayerFaceIndices(icvCenter, faInternal, faBoundary);
	// internal faces
	for(size_t index=0; index<faInternal.size(); ++index) {
		int ifa = faInternal[index];
		FaZone* zone = zofa[ifa];

		diffusivityHookScalarRansTurb1D_AD(ifa, zone, scalName);
		diffusivityHookScalarRansComb1D_AD(ifa, zone, scalName);
	}
	for(size_t index=0; index<faBoundary.size(); ++index) {
		int ifa = faBoundary[index];
		FaZone* zone = zofa[ifa];
		diffusivityHookScalarRansTurb1D_AD(ifa, zone, scalName);
		diffusivityHookScalarRansComb1D_AD(ifa, zone, scalName);
	}

	// Calculate flux
	int nofa_f = faocv_i[icvCenter];
	int nofa_l = faocv_i[icvCenter+1]-1;
	for (int nofa = nofa_f; nofa <= nofa_l; nofa++) {
		int ifa = faocv_v[nofa];

		// =============================================================================================
		// cycle trough internal faces and compute viscous flux and implicit matrix
		// =============================================================================================
		if((ifa >= nfa_b && ifa<nfa) || (ifa >= nfa_b2 && ifa < nfa_b2gg)) {
			int icv0 = cvofa[ifa][0];
			int icv1 = cvofa[ifa][1];
			int noc00, noc01, noc11, noc10;

			if (flagImplicit && icv0 < ncv && icv1 < ncv_g)
				getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

			double n[3] = {0.0, 0.0, 0.0}, s[3] = {0.0, 0.0, 0.0};
			double nmag = normVec3d(n, fa_normal[ifa]);
			vecMinVec3d(s, x_cv[icv1], x_cv[icv0]);
			double smag = normVec3d(s);
			double alpha = vecDotVec3d(n, s);
			assert((alpha > 0.0) && (alpha < 1.0+1.0E-12));   // alpha should now contain s dot n...

			double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
			vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
			vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);

			adouble viscFlux = diff[ifa] * nmag * (alpha * (phi[icv1] - phi[icv0]) / smag
					+ 0.5 * (grad_phi[icv0][0] + grad_phi[icv1][0]) * (n[0] - alpha * s[0])
					+ 0.5 * (grad_phi[icv0][1] + grad_phi[icv1][1]) * (n[1] - alpha * s[1])
					+ 0.5 * (grad_phi[icv0][2] + grad_phi[icv1][2]) * (n[2] - alpha * s[2]));

			if(icvCenter == icv0)
				rhs_rhoScal += viscFlux;
			else if(icvCenter == icv1)
				rhs_rhoScal -= viscFlux;
			else {
				printf("ERROR in IkeWithModels_AD::calcViscousFluxScalar_new1D_AD() => icvCenter=%d is not equal to icv0=%d nor icv1=%d \n", icvCenter, icv0, icv1);
			}

			if (flagImplicit && icv0 < ncv && icv1 < ncv_g) {
				AScal[noc00] +=   diff[ifa].value() * nmag * alpha / smag / rho_AD[icv0].value();
				AScal[noc01] += - diff[ifa].value() * nmag * alpha / smag / rho_AD[icv0].value();

				if (icv1 < ncv) {
					AScal[noc11] +=   diff[ifa].value() * nmag * alpha / smag / rho_AD[icv1].value();
					AScal[noc10] += - diff[ifa].value() * nmag * alpha / smag / rho_AD[icv1].value();
				}
			}
		}
		// =============================================================================================
		// cycle trough boundary faces and compute viscous flux and implicit matrix
		// =============================================================================================
		else {
			FaZone* zone = zofa[ifa];
			// .............................................................................................
			// HOOK BOUNDARY CONDITION
			// .............................................................................................
			if(getBoundaryTypeScalar(ifa, zone, scalName) == HOOK_SCALAR)  {
				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];
				int noc00 = nbocv_i[icv0];

				double n[3] = {0.0, 0.0, 0.0}, s_half[3] = {0.0, 0.0, 0.0};
				double nmag = normVec3d(n, fa_normal[ifa]);
				vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
				double smag_half = fabs(vecDotVec3d(s_half, n));   //normVec3d(s_half);
				//double alpha = vecDotVec3d(n, s_half);           // HACK alfa should be always 1 at the boundary (provide such a grid!)

				adouble viscFlux = diff[ifa] * nmag * (phi[icv1] - phi[icv0]) / smag_half;

				if(icvCenter == icv0)
					rhs_rhoScal += viscFlux;

				if (flagImplicit && icv0<ncv)
					AScal[noc00] += diff[ifa].value() * nmag / smag_half / rho_AD[icv0].value();
			}
			// .............................................................................................
			// DIRICHLET BOUNDARY CONDITION
			// .............................................................................................
			else if (getBoundaryTypeScalar(ifa, zone, scalName) == DIRICHLET_SCALAR){
			    Param *p;
			    string fName = zone->getName();
			    if (getParam(p, fName + "." + scalName))
			    	phiBCval = p->getDouble(1);
			    else {
			    	cout<<"Error in IkeWithModels_AD::calcViscousFluxScalar_new1D_AD() => Cannot find Param for "<<fName+"."+scalName<<" at ifa="<<ifa<<endl;
			    	throw(-1);
			    }

				int icv0 = cvofa[ifa][0];
				int icv1 = cvofa[ifa][1];
				int noc00 = nbocv_i[icv0];

				double n[3] = {0.0, 0.0, 0.0}, s_half[3] = {0.0, 0.0, 0.0};
				double nmag = normVec3d(n, fa_normal[ifa]);
				vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
				double smag_half = fabs(vecDotVec3d(s_half, n));   //normVec3d(s_half);
				//double alpha = vecDotVec3d(n, s_half);           // HACK alfa should be always 1 at the boundary (provide such a grid!)

				adouble viscFlux = diff[ifa] * nmag * (phi[icv1] - phi[icv0]) / smag_half;                // TODO: phi_bfa or phiBCval??

				if(icvCenter == icv0)
					rhs_rhoScal += viscFlux;

				if (flagImplicit && icv0<ncv)
					AScal[noc00] += diff[ifa].value() * nmag / smag_half / rho_AD[icv0].value();
			}
			// .............................................................................................
			// FLUX BOUNDARY CONDITION: NOT IMPLEMENTED YET
			// .............................................................................................
			else if (getBoundaryTypeScalar(ifa, zone, scalName) == FLUX_SCALAR) {
//				Param *p;
//				string fName = zone->getName();
//				if (getParam(p, fName + "." + scalName))
//					phiBCflux = p->getDouble(1);
//				else {
//					cout<<"Error in IkeWithModels_AD::calcViscousFluxScalar_new1D_AD() => Cannot find Param for "<<fName+"."+scalName+".flux"<<" at ifa="<<ifa<<endl;
//					throw(-1);
//				}
				cerr << "scalarZoneIsFlux not implemented yet!" << endl;
				throw(-1);
			}
			// .............................................................................................
			// OTHER BOUNDARY CONDITIONS: NO VISCOUS FLUX
			// .............................................................................................
			else {
				// note: SYMMETRY
			}
		}
	}
#endif
}

/*
 * Method: showMessageParallel
 * ---------------------------
 * Each CPU core shows its message (usually, error or warning message) in a mpi_rank order
 *
 * Note: 1. All the CPU cores must call this method (Otherwise, communication error)
 *       2. If the string is empty, NO message will be shown
 */
void IkeWithModels_AD::showMessageParallel(string& message) {
	if(mpi_rank != 0) { MPI_Status status; int dummy; MPI_Recv(&dummy, 1, MPI_INT, mpi_rank-1, 1115, mpi_comm, &status); }

	if(!message.empty())
		cout<<message<<" (mpi_rank="<<mpi_rank<<")"<<endl;
	for(int i=0; i<10000; ++i) {} // Pause for a while

	if(mpi_rank < mpi_size-1) { int dummy; MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 1115, mpi_comm); }

	MPI_Barrier(mpi_comm);
}

void IkeWithModels_AD::showMessageParallel(vector<string>& messageVector, const size_t maxMessageNum, char* messageVectorName /*= NULL*/) {
	int mySize = messageVector.size();
	int totSize;
	MPI_Allreduce(&mySize, &totSize, 1, MPI_INT, MPI_SUM, mpi_comm);

	if(mpi_rank != 0) { MPI_Status status; int dummy; MPI_Recv(&dummy, 1, MPI_INT, mpi_rank-1, 1118, mpi_comm, &status); }

	if(!messageVector.empty()) {
		if(totSize > maxMessageNum) {
    		if(mpi_rank == 0) {
    			cout<<"IkeWithModels_AD::showMessageParallel(): TOO MANY MESSAGES ("<< messageVector.size() <<")";
    			if(messageVectorName != NULL)
    				cout<<" FOR "<<messageVectorName;
    			cout<<" at mpi_rank="<<mpi_rank<<endl
    					<<"                                         SHOW ONLY FIRST 100 MESSGAES: "<<endl;
    		}

			for(size_t i=0; i<maxMessageNum; ++i)
				cout<<messageVector[i];
		} else {
			for(size_t i=0; i<messageVector.size(); ++i)
				cout<<messageVector[i];
		}
	}
	for(int i=0; i<10000; ++i) {} // Pause for a while

	if(mpi_rank < mpi_size-1) { int dummy; MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 1118, mpi_comm); }

	MPI_Barrier(mpi_comm);
}
