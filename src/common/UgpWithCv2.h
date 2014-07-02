#ifndef UGPWITHCV2_H
#define UGPWITHCV2_H

#include "UgpWithTools.h"
#include "tc_vec3d.h"

#define FOR_ICV_G       for (int icv = 0; icv < ncv_g; ++icv)
#define FOR_ICV_GG      for (int icv = 0; icv < ncv_gg; ++icv)
#define FOR_ICV_GGF     for (int icv = 0; icv < ncv_ggf; ++icv)
#define FOR_ICV_GGFF    for (int icv = 0; icv < ncv_ggff; ++icv)
#define FOR_ICV_ALL     for (int icv = 0; icv < ncv_ggff; ++icv)
#define FOR_ICV_G1_ONLY for (int icv = ncv; icv < ncv_g; ++icv)
#define FOR_ICV_G2_ONLY for (int icv = ncv_g; icv < ncv_gg; ++icv)
#define FOR_ICV_F1_ONLY for (int icv = ncv_gg; icv < ncv_ggf; ++icv)
#define FOR_ICV_F2_ONLY for (int icv = ncv_ggf; icv < ncv_ggff; ++icv)
#define FOR_IFA_B2_ONLY for (int ifa = nfa; ifa < nfa_b2; ++ifa)
#define FOR_IFA_GG_ONLY for (int ifa = nfa_b2; ifa < nfa_b2gg; ++ifa)
#define FOR_IFA_ALL     for (int ifa = 0; ifa < nfa_b2gg; ++ifa)

#define NO_GRADIENTS 100
#define LSG_GRADIENTS 101
#define MGG_GRADIENTS 102
#define GG_GRADIENTS  103  ///< Green-Gauss with arithmetic mean


class UgpWithCv2: public UgpWithTools {
  
 private:
  
  list<Prcomm> cvPrcommListG1;
  list<Prcomm> cvPrcommListG1G2;
  list<Prcomm> cvPrcommListG1G2F2;
  
 protected:  // member variables
  
  // ghost and fake cv ranges...
  int ncv_g;    /// includes level-1 ghosts
  int ncv_gg;   /// includes level-2 ghosts (face-based nbrs of level-1 ghosts) 
  int ncv_ggf;  /// includes fakes owned by nfa_b faces: i.e. ncv_ggf-ncv_gg == nfa_b 
  int ncv_ggff; /// includes fakes that are nbrs of level-1 ghosts

  // new boundary faces associated with additional fakes of level-0 ghost cvs...
  // int nfa_b are the first boundary faces already defined 
  // int nfa are all face associated with the ncv cvs
  int nfa_b2;   // additional boundary faces are between g1:f2 cvs
  int nfa_b2g;  // faces between g1:g1 cvs
  int nfa_b2gg; // faces between g1:g2 cvs
    
  // cv geometry - eventually gets defined over full range of cvs [0:ncv_ggff-1]...
  double *cv_volume;          ///< cell volume
  double (*x_cv)[3];          ///< cell centroid coordinate
  
  // fa geometry: defined over [0:nfa-1]...
  double (*fa_normal)[3];     ///< face normal area weighted
  double (*x_fa)[3];          ///< face center coordinate
  
  // neighbor-of-cv CSR structure: includes self (diagonal) as
  // first element, but does NOT include any fake cvs. defined only
  // over locally owned cvs, i.e. [0:ncv-1]...
  int * nbocv_i;
  int nbocv_s;
  int * nbocv_v;

  // FULL neighbor-of-cv CSR structure: includes self (diagonal) as
  // first element and is also defined for level-0 ghosts, where nbrs will include
  // level-1 ghosts and potentially additional fakes...
  int * nbocv_all_i;
  int * nbocv_all_v;
  int nbocv_all_s;
  
  // gradient operator - defined in all local cvs AND
  // level-0 ghosts...
  double (*nbocv_all_grad_coeff)[3];
  
  // filter operator - defined in all local cvs AND
  // level-0 ghosts...
  double *nbocv_all_filter_coeff;
  
 public:   // constructors/destructors

  UgpWithCv2() {
    
    if (mpi_rank == 0)
      cout << "UgpWithCv2()" << endl;

    ncv_g = 0;
    ncv_gg = 0;
    ncv_ggf = 0;
    ncv_ggff = 0;

    nfa_b2 = 0;
    nfa_b2g = 0;
    nfa_b2gg = 0;

    cv_volume = NULL;
    x_cv = NULL;
    
    fa_normal = NULL;
    x_fa = NULL;

    nbocv_i = nbocv_v = NULL;
    nbocv_s = 0;

    nbocv_all_i = nbocv_all_v = NULL;
    nbocv_all_s = 0;
    
    nbocv_all_grad_coeff = NULL;
    
    nbocv_all_filter_coeff = NULL;
    
  }
  
  virtual ~UgpWithCv2() {

    if (cv_volume != NULL) delete[] cv_volume;
    if (x_cv != NULL)      delete[] x_cv;
    
    if (fa_normal != NULL) delete[] fa_normal;
    if (x_fa != NULL)      delete[] x_fa;
    
    if (nbocv_i != NULL)   delete[] nbocv_i;
    if (nbocv_v != NULL)   delete[] nbocv_v;
    
    if (nbocv_all_i != NULL)   delete[] nbocv_all_i;
    if (nbocv_all_v != NULL)   delete[] nbocv_all_v;
    
    if (nbocv_all_grad_coeff != NULL) delete[] nbocv_all_grad_coeff;
    
    if (nbocv_all_filter_coeff != NULL) delete[] nbocv_all_filter_coeff;
    
  }
  
 public:   // member functions
  
  void init(const int grad_flag = NO_GRADIENTS) {
    
    // call UgpWithTools's init routine...
    UgpWithTools::init();
    
    // add 2 levels of ghost and fakes, build
    // cv-based paired communicators, etc...
    addGhostAndFakeCvs();

    // add new boundary faces...
    addB2AndGhostFaces();
    
    // allocate and compute the geometry...
    calcGeometry();
    
    // if this is a uniform hex grid, you can do this check...
    uniformHexCheck();
    
    // build nbocv_i/v CSR structure that contains local and ghosts
    // only...
    buildNbocv();

    // gradient and filter operators...
    switch (grad_flag) {
    case GG_GRADIENTS:
          calcGradCoeffGG();
          break;
    case NO_GRADIENTS:
      // default is no gradient operators...
      if (mpi_rank == 0)
	cout << " > no gradient operators built" << endl;
      break;
    case MGG_GRADIENTS:
      // modified Green-Gauss gradients...
      calcGradAndFilterCoeffMGG();
      checkCvGrad();
      checkCv2Grad();
      checkCvFilter();
      checkCv2Filter();
      break;
    default:
      if (mpi_rank == 0)
	cout << "Error: unrecognized grad_flag: " << grad_flag << endl;
      throw(-1);
    }

  }

 public:
  /*
   * Function: isNaN0
   * ----------------
   * Check if the given double variable is NaN
   */
  inline bool isNaN0(const double num) {
  	return (num) != (num);
  }

  
  void updateFaDataB2G1G2(double (*v)[3],const int action) {
    
    // do not use this routine inside a solver. Recommended just
    // for preprocessing (e.g. exchanging geometry)...
    
    int my_nfoc_max = 0;
    FOR_ICV my_nfoc_max = max( my_nfoc_max, faocv_i[icv+1]-faocv_i[icv] );
    int nfoc_max;
    MPI_Allreduce(&my_nfoc_max,&nfoc_max,1,MPI_INT,MPI_MAX,mpi_comm);
    
    double (*buf)[3] = new double[ncv_ggff][3]; 
    FOR_ICV_ALL FOR_I3 buf[icv][i] = 0.0;
    
    // report the max diff as well...
    double my_d2_max = 0.0;

    // for deciding when the diff can be calculated...
    FOR_IFA_ALL fa_flag[ifa] = 0;
    
    // packing/unpacking is different depending on whether the quantity is an absolute quantity (e.g. x) or
    // relative signed quantity (eg normal)...
    
    if (action == REPLACE_TRANSLATE_DATA) {
     
      // absolute quantity...
      
      for (int foc_inc = 0; foc_inc < nfoc_max; ++foc_inc) {

	// pack...
	// here we should be able to cycle through just a subset of these internal
	// cvs that actually correspond to ghost cvs somewhere (and for that matter
	// nfoc_max may be too large) but to simplfy this routine, just do all internal
	// cvs...
	FOR_ICV {
	  int foc = faocv_i[icv] + foc_inc;
	  if (foc < faocv_i[icv+1]) {
	    int ifa = faocv_v[foc];
	    assert( (ifa >= 0)&&(ifa < nfa) );
	    FOR_I3 buf[icv][i] = v[ifa][i];
	  }
	}
	// also cycle through the boundary faces and populate F1 cvs with the 
	// boundary face normals. Note that these are outward pointing by default...
	FOR_IFA_B {
	  // because this is done every iteration, don't check, just set...
	  fa_flag[ifa] = 1;
	  int icv1 = cvofa[ifa][1];
	  assert( (icv1 >= ncv_gg)&&(icv1 < ncv_ggf) );
	  FOR_I3 buf[icv1][i] = v[ifa][i];
	}

	// exchange...
	updateCvDataG1G2F2(buf,REPLACE_TRANSLATE_DATA);

	// F2 cvs...
	FOR_IFA_B2_ONLY {
	  // because this is done every iteration, don't check, just set...
	  fa_flag[ifa] = 1;
	  int icv1 = cvofa[ifa][1];
	  assert( (icv1 >= ncv_ggf)&&(icv1 < ncv_ggff) );
	  FOR_I3 v[ifa][i] = buf[icv1][i];
	}
	
	// unpack...
	FOR_ICV_G1_ONLY {
	  int foc = faocv_i[icv] + foc_inc;
	  if (foc < faocv_i[icv+1]) {
	    int ifa = faocv_v[foc];
	    if (ifa < nfa) {
	      // this is an old face, so its value should be already present...
	      assert( fa_flag[ifa] == 0 );
	      fa_flag[ifa] = 1;
	      double this_d2 = 0.0;
	      FOR_I3 {
		double dx = v[ifa][i] - buf[icv][i];
		this_d2 += dx*dx;
	      }
	      my_d2_max = max(my_d2_max,this_d2);
	    }
	    else if (ifa < nfa_b2) {
	      // b2 boundary face... 
	      assert( cvofa[ifa][0] == icv );
	      assert( fa_flag[ifa] == 1 );
	      double this_d2 = 0.0;
	      FOR_I3 {
		double dx = v[ifa][i] - buf[icv][i];
		this_d2 += dx*dx;
	      }
	      my_d2_max = max(my_d2_max,this_d2);
	    }
	    else if (ifa < nfa_b2g) {
	      // face between G1 and another G1 cv. these ones will be exchanged
	      // twice, and the order is not guaranteed to be smaller icv first, because this
	      // routine is done one neighbor at a time, so use the fa_flag to indicate 
	      // is a value is there already...
	      if (fa_flag[ifa] == 0) {
		fa_flag[ifa] = 1;
		FOR_I3 v[ifa][i] = buf[icv][i];
	      }
	      else {
		double this_d2 = 0.0;
		FOR_I3 {
		  double dx = v[ifa][i] - buf[icv][i];
		  this_d2 += dx*dx;
		}
		my_d2_max = max(my_d2_max,this_d2);
	      }
	    }
	    else {
	      // face between G1 and G2 cv - should be the only one...
	      assert( ifa < nfa_b2gg );
	      assert( fa_flag[ifa] == 0 );
	      fa_flag[ifa] = 1;
	      assert( cvofa[ifa][0] == icv );
	      FOR_I3 v[ifa][i] = buf[icv][i];
	    }
	  }
	}

      }
    }
    else if (action == REPLACE_ROTATE_DATA) {
      
      // vector quantity the sign in the face direction...

      for (int foc_inc = 0; foc_inc < nfoc_max; ++foc_inc) {

	// pack...
	// here we should be able to cycle through just a subset of these internal
	// cvs that actually correspond to ghost cvs somewhere (and for that matter
	// nfoc_max may be too large) but to simplify this routine, just do all internal
	// cvs...
	FOR_ICV {
	  int foc = faocv_i[icv] + foc_inc;
	  if (foc < faocv_i[icv+1]) {
	    int ifa = faocv_v[foc];
	    assert( (ifa >= 0)&&(ifa < nfa) );
	    // back an outward-pointing normal wrt icv...
	    if (cvofa[ifa][0] == icv) {
	      // face is already outard-pointing...
	      FOR_I3 buf[icv][i] = v[ifa][i];
	    }
	    else {
	      // face is inward-pointing wrt icv, so flip...
	      assert(cvofa[ifa][1] == icv);
	      FOR_I3 buf[icv][i] = -v[ifa][i];
	    }
	  }
	}
	// also cycle through the boundary faces and populate F1 cvs with the 
	// boundary face normals. Note that these are outward pointing by default...
	FOR_IFA_B {
	  // because this is done every iteration, don't check, just set...
	  fa_flag[ifa] = 1;
	  int icv1 = cvofa[ifa][1];
	  assert( (icv1 >= ncv_gg)&&(icv1 < ncv_ggf) );
	  FOR_I3 buf[icv1][i] = v[ifa][i];
	}

	// exchange...
	updateCvDataG1G2F2(buf,REPLACE_ROTATE_DATA);
	
	// F2 cvs...
	FOR_IFA_B2_ONLY {
	  // because this is done every iteration, don't check, just set...
	  fa_flag[ifa] = 1;
	  int icv1 = cvofa[ifa][1];
	  assert( (icv1 >= ncv_ggf)&&(icv1 < ncv_ggff) );
	  FOR_I3 v[ifa][i] = buf[icv1][i];
	}

	// unpack...
	FOR_ICV_G1_ONLY {
	  int foc = faocv_i[icv] + foc_inc;
	  if (foc < faocv_i[icv+1]) {
	    int ifa = faocv_v[foc];
	    if (ifa < nfa) {
	      // this is an old face, so its value should be already present...
	      assert( fa_flag[ifa] == 0 );
	      fa_flag[ifa] = 1;
	      // also, it should be inward-pointing with respect to icv (this is
	      // the G1 cv), the compare the flipped version of the normal...
	      assert( cvofa[ifa][1] == icv );
	      double this_d2 = 0.0;
	      FOR_I3 {
		double dx = v[ifa][i] + buf[icv][i]; // changed to + to flip normal
		this_d2 += dx*dx;
	      }
	      my_d2_max = max(my_d2_max,this_d2);
	    }
	    else if (ifa < nfa_b2) {
	      // b2 boundary face... 
	      assert( cvofa[ifa][0] == icv );
	      assert( fa_flag[ifa] == 1 );
	      double this_d2 = 0.0;
	      FOR_I3 {
		double dx = v[ifa][i] - buf[icv][i];
		this_d2 += dx*dx;
	      }
	      my_d2_max = max(my_d2_max,this_d2);
	    }
	    else if (ifa < nfa_b2g) {
	      // face between G1 and another G1 cv. these ones will be exchanged
	      // twice, and the order is not guaranteed to be smaller icv first, because this
	      // routine is done one neighbor at a time, so use the fa_flag to indicate 
	      // is a value is there already...
	      if (fa_flag[ifa] == 0) {
		fa_flag[ifa] = 1;
		if (cvofa[ifa][0] == icv) {
		  FOR_I3 v[ifa][i] = buf[icv][i];
		}
		else {
		  assert(cvofa[ifa][1] == icv);
		  FOR_I3 v[ifa][i] = -buf[icv][i];
		}
	      }
	      else {
		// this one has been already set, so just check...
		if (cvofa[ifa][0] == icv) {
		  double this_d2 = 0.0;
		  FOR_I3 {
		    double dx = v[ifa][i] - buf[icv][i];
		    this_d2 += dx*dx;
		  }
		  my_d2_max = max(my_d2_max,this_d2);
		}
		else {
		  assert(cvofa[ifa][1] == icv);
		  double this_d2 = 0.0;
		  FOR_I3 {
		    double dx = v[ifa][i] + buf[icv][i];
		    this_d2 += dx*dx;
		  }
		  my_d2_max = max(my_d2_max,this_d2);
		}
	      }
	    }
	    else {
	      // face between G1 and G2 cv - should be the only one...
	      assert( ifa < nfa_b2gg );
	      assert( fa_flag[ifa] == 0 );
	      fa_flag[ifa] = 1;
	      assert( cvofa[ifa][0] == icv );
	      FOR_I3 v[ifa][i] = buf[icv][i];
	    }
	  }
	}

      }

    }
    else {
      
      if (mpi_rank == 0)
	cout << "Error: unsupported action: " << action << endl;
      throw(-1);
      
    }
   
    delete[] buf;
    
    // report d2...
    double d2_max;
    MPI_Reduce(&my_d2_max,&d2_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    if (mpi_rank == 0)
      cout << " > updateFaDataB2GG: max error (should be zero): " << sqrt(d2_max) << endl;
    
    // all faces should be set...
    for (int ifa = nfa; ifa < nfa_b2gg; ++ifa)
      assert( fa_flag[ifa] == 1 );
    
  }
  
 private:

  // ---------------------------
  // initialization routines
  // ---------------------------
  
  void addGhostAndFakeCvs();

  void addB2AndGhostFaces();

  void reorderPrcomm(list<Prcomm>& prcommList,const int * order,const int n) {
    
    // here we assume that any new order does not impact the range dilineation
    // that mat be present in this prcomm... 
    
    for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm != prcommList.end(); prcomm++) {
      // pack...
      for (int i = 0; i < prcomm->packIndexVec.size(); i++) {
	int index_old = prcomm->packIndexVec[i];
	assert( (index_old >= 0)&&(index_old < n) );
	prcomm->packIndexVec[i] = order[index_old];
      }
      // unpack...
      for (int i = 0; i < prcomm->unpackIndexVec.size(); i++) {
	int index_old = prcomm->unpackIndexVec[i];
	assert( (index_old >= 0)&&(index_old < n) );
	prcomm->unpackIndexVec[i] = order[index_old];
      }
    }
    
  }

  void calcGeometry();
  
  void uniformHexCheck();

  void buildNbocv();

  void calcGradCoeffGG();

  void calcGradAndFilterCoeffMGG();

  int getNextNodeOfFaceCCW(const int ifa,const int ino,const int icv) {

    int nof_f = noofa_i[ifa];
    int nof_l = noofa_i[ifa+1]-1;
    for (int nof = nof_f; nof <= nof_l; ++nof) {
      int this_ino = noofa_v[nof];
      if (this_ino == ino) {
	if (cvofa[ifa][0] == icv) {
	  if (nof == nof_l)
	    return( noofa_v[nof_f] );
	  else
	    return( noofa_v[nof+1] );
	}
	else {
	  assert( cvofa[ifa][1] == icv );
	  if (nof == nof_f)
	    return( noofa_v[nof_l] );
	  else
	    return( noofa_v[nof-1] );
	}
      }
    }

    cerr << "Error: could not find node " << ino << " in face " << ifa << endl;
    throw(-1);
    
  }

  int getPrevNodeOfFaceCCW(const int ifa,const int ino,const int icv) {

    int nof_f = noofa_i[ifa];
    int nof_l = noofa_i[ifa+1]-1;
    for (int nof = nof_f; nof <= nof_l; ++nof) {
      int this_ino = noofa_v[nof];
      if (this_ino == ino) {
	if (cvofa[ifa][0] == icv) {
	  if (nof == nof_f)
	    return( noofa_v[nof_l] );
	  else
	    return( noofa_v[nof-1] );
	}
	else {
	  assert( cvofa[ifa][1] == icv );
	  if (nof == nof_l)
	    return( noofa_v[nof_f] );
	  else
	    return( noofa_v[nof+1] );
	}
      }
    }

    cerr << "Error: could not find node " << ino << " in face " << ifa << endl;
    throw(-1);
    
  }
  
  void checkCvGrad();
  
  void checkCv2Grad();
  
  void checkCvFilter();
  
  void checkCv2Filter();
  
 public:

  void calcCvGrad(double (*dphidxi)[3],const double * phi) {

    // populates the gradient in cv range [0:ncv-1]... i.e. there is no setting
    // or reduction to ghost cvs, and there is nothing set in fake data (if present)... 
    
    FOR_ICV {
      int noc_f = nbocv_all_i[icv];
      FOR_I3 dphidxi[icv][i] = nbocv_all_grad_coeff[noc_f][i]*phi[icv];
      int noc_l = nbocv_all_i[icv+1]-1;
      for (int noc = noc_f+1; noc <= noc_l; ++noc) {
	int icv_nbr = nbocv_all_v[noc];
	FOR_I3 dphidxi[icv][i] += nbocv_all_grad_coeff[noc][i]*phi[icv_nbr];
      }
    }

  }

  void calcCv2ScalGrad(double (*dphidxi)[3],const double * phi)
  {
	  // populates the gradient in cv range [0:ncv_g-1]. Note that this requires
	  // phi to be defined over all cv's [0:ncv_ggff-1]...
	  FOR_ICV_G {
		  int noc_f = nbocv_all_i[icv];
		  FOR_I3 dphidxi[icv][i] = nbocv_all_grad_coeff[noc_f][i]*phi[icv];
		  int noc_l = nbocv_all_i[icv+1]-1;
		  for (int noc = noc_f+1; noc <= noc_l; ++noc) {
			  int icv_nbr = nbocv_all_v[noc];
			  FOR_I3 dphidxi[icv][i] += nbocv_all_grad_coeff[noc][i]*phi[icv_nbr];
		  }

		  if(isNaN0(dphidxi[icv][0]) || isNaN0(dphidxi[icv][1]) || isNaN0(dphidxi[icv][2])) {
			  printf("ERROR in UgpWithCv2::calcCv2ScalGrad(): NaN detected in dphidxi at icv=%d(%.2e,%.2e,%.2e) \n", icv, x_cv[icv][0], x_cv[icv][1], x_cv[icv][2]);
			  printf("         Details: mpi_rank=%d;  ncv=%d, ncv_g=%d, ncv_gg=%d, ncv_ggf=%d, ncv_ggff=%d\n", mpi_rank, ncv, ncv_g, ncv_gg, ncv_ggf, ncv_ggff);
			  cout<< "                  ADDRESS of phi="<<phi<<", ADDRESS of dphidxi="<<dphidxi<<endl;
			  printf("                  icv=%d (%.4e,%.4e,%.4e): phi=%.3e \n",icv,x_cv[icv][0],x_cv[icv][1],x_cv[icv][2],phi[icv]);
			  int noc_l = nbocv_all_i[icv+1]-1;
			  for (int noc = noc_f+1; noc <= noc_l; ++noc) {
				  int icv_nbr = nbocv_all_v[noc];
				  printf("                  icv_nbr=%d (%.4e,%.4e,%.4e): phi=%.3e \n",icv_nbr,x_cv[icv_nbr][0],x_cv[icv_nbr][1],x_cv[icv_nbr][2],phi[icv_nbr]);
			  }
			  throw(-11);
		  }
	  }
  }

  /**
   * populates the gradient in cv range [0:ncv_g-1]. Note that this requires
   * phi to be defined over all cv's [0:ncv_ggff-1]...
   */
  void calcCv2VecGrad(double(*duidxj)[3][3], const double(*u)[3])
  {
    FOR_ICV_G
    {
      int noc_f = nbocv_all_i[icv];
      FOR_I3
        FOR_J3
          duidxj[icv][i][j] = nbocv_all_grad_coeff[noc_f][j]*u[icv][i];

      int noc_l = nbocv_all_i[icv+1]-1;
      for (int noc = noc_f+1; noc<=noc_l; ++noc) 
      { 
        int icv_nbr = nbocv_all_v[noc];
        FOR_I3
          FOR_J3
            duidxj[icv][i][j] += nbocv_all_grad_coeff[noc][j]*u[icv_nbr][i];
      }
    }
  }

  
  void filterCvData(double * phi_hat,const double * phi) {
    
    /*
    // uncomment following and run on 1 processor 
    // to inspect filter coefficients...
    FOR_ICV {
    int noc_f = nbocv_all_i[icv];
    cout << "diag: " << nbocv_all_filter_coeff[noc_f] << endl;
    int noc_l = nbocv_all_i[icv+1]-1;
    for (int noc = noc_f+1; noc <= noc_l; ++noc) {
    cout << "filter coeff: " << nbocv_all_filter_coeff[noc] << endl;
    }
    getchar();
    }
    */
    
    FOR_ICV {
      int noc_f = nbocv_all_i[icv];
      phi_hat[icv] = nbocv_all_filter_coeff[noc_f]*phi[icv];
      int noc_l = nbocv_all_i[icv+1]-1;
      for (int noc = noc_f+1; noc <= noc_l; ++noc) {
	int icv_nbr = nbocv_all_v[noc];
	phi_hat[icv] += nbocv_all_filter_coeff[noc]*phi[icv_nbr];
      }
    }

  }

  void filterCv2Data(double * phi_hat,const double * phi) {

    // include filter in ghost level-0 data - note that this requires
    // phi to be defined over all cv's [0:ncv_ggff-1]...
    
    FOR_ICV_G {
      int noc_f = nbocv_all_i[icv];
      phi_hat[icv] = nbocv_all_filter_coeff[noc_f]*phi[icv];
      int noc_l = nbocv_all_i[icv+1]-1;
      for (int noc = noc_f+1; noc <= noc_l; ++noc) {
	int icv_nbr = nbocv_all_v[noc];
	phi_hat[icv] += nbocv_all_filter_coeff[noc]*phi[icv_nbr];
      }
    }

  }

  // ---------------------------
  // ghost/fake data updates
  // ---------------------------

  void updateCvDataG1(double * d, const int action) {
    updateR1(d, action, cvPrcommListG1);
  }
  
  void updateCvDataG1(int * d, const int action) {
    updateI1(d, action, cvPrcommListG1);
  }

  void updateCvDataG1(double(*d)[3], const int action) {
    updateR2(d, action, cvPrcommListG1);
  }

  void updateCvDataG1(double (*d)[3][3],const int action) {
    updateR3(d,action, cvPrcommListG1);
  }

  void updateCvSymmetricR3G1(double (*diag)[3],double (*offd)[3],const int action) {
    updateSymmetricR3(diag,offd,action,cvPrcommListG1);
  }
  
  void updateCvSymmetricR3G1(double (*d)[6],const int action) {
    // here the 6 tensor elements are stored as 
    // xx,xy,xz,yy,yz,zz
    updateSymmetricR3(d,action,cvPrcommListG1);
  }
  
  // these routines update operators into the G1 cvs...

  void updateCvOperatorG1(double * A);
  void updateCvOperatorG1(double (*A)[3]);

  // ---------------------------
  // full ghost data updates
  // ---------------------------
  
  void updateCvDataG1G2(double * d, const int action) {
    updateR1(d, action, cvPrcommListG1G2);
  }
  
  void updateCvDataG1G2(int * d, const int action) {
    updateI1(d, action, cvPrcommListG1G2);
  }
  
  void updateCvDataG1G2(double(*d)[3], const int action) {
    updateR2(d, action, cvPrcommListG1G2);
  }
  
  void updateCvDataG1G2(double (*d)[3][3],const int action) {
    updateR3(d,action, cvPrcommListG1G2);
  }
  
  void updateCvSymmetricR3G1G2(double (*d)[6],const int action) {
    // here the 6 tensor elements are stored as 
    // xx,xy,xz,yy,yz,zz
    updateSymmetricR3(d,action,cvPrcommListG1G2);
  }
  
  // ---------------------------
  // full ghost and fake data updates
  // ---------------------------
  
  void updateCvDataG1G2F2(double * d, const int action) {
    updateR1(d, action, cvPrcommListG1G2F2);
  }
  
  void updateCvDataG1G2F2(int * d, const int action) {
    updateI1(d, action, cvPrcommListG1G2F2);
  }
  
  void updateCvDataG1G2F2(double(*d)[3], const int action) {
    updateR2(d, action, cvPrcommListG1G2F2);
  }
  
  void updateCvDataG1G2F2(double (*d)[3][3],const int action) {
    updateR3(d,action, cvPrcommListG1G2F2);
  }
  
  void updateCvSymmetricR3G1G2F2(double (*d)[6],const int action) {
    // here the 6 tensor elements are stored as 
    // xx,xy,xz,yy,yz,zz
    updateSymmetricR3(d,action,cvPrcommListG1G2F2);
  }

  void updateCvDataByName(char *name, const int action) {
    double *scal = getScalar(name);
    updateCvDataG1G2(scal, action);
  }
    
};



//######################################################################################
//######################################################################################
//######################################################################################
// put into another file
//######################################################################################
//######################################################################################
//######################################################################################

#define NOLIMITER            201
#define BARTH_JESPERSEN_MOD  202

#define BCGSTAB              301
#define PETSC_GMRES          302

#define SDWLSClippingScalar   1.0e-8

#ifndef sign
#define sign(a)       ( ((a) > 0) ? (1.0) : ((a) < 0.0 ? -1.0 : 0.0))
#endif



class ScalarTranspEq: public Data
{
public:
  double *phi;              ///< scalar field
  double *diff;             ///< scalar diffusion coefficient D ... d/dxi(D dphi/dxi)
  double (*grad_phi)[3];    ///< scalar gradients

  double *rhophi;                         ///< conserved scalar field (rho * phi)
  double (*grad_rhophi)[3];               ///< conserved scalar gradients (rho * phi)

  double *dpress_dphi;                    ///< derivative of pressure with respect to scalar, used for coupled solver

  double turbSchmidtNumber;
  bool convTerm, diffTerm;

  int phiMaxiter;
  double relax;
  double phiZero, phiZeroRel;
  double resMax, resAve;
  double lowerBound, upperBound;
  string reconstruction;
  string coupling;             ///< determine if scalar is solved coupled or uncoupled from NSE ("COUPLED" or "UNCOUPLED")


  ScalarTranspEq(const char * name, int datatype) : Data(name, datatype)
  {
    phi = NULL;
    diff = NULL;
    grad_phi = NULL;
    rhophi = NULL;
    grad_rhophi = NULL;

   dpress_dphi = NULL;                    ///< derivative of pressure with respect to scalar, used for coupled solver


    turbSchmidtNumber = 1.0;

    convTerm = true;
    diffTerm = true;

    phiMaxiter = 100;
    relax = 1.0;
    phiZero = 1.0e-6;
    phiZeroRel = 1.0e-3;

    lowerBound = 0.0;
    upperBound = 1.0e10;

    reconstruction = "STANDARD";
    coupling = "UNCOUPLED";
  }

  ~ScalarTranspEq() {
	  if(phi != NULL) {
		  delete [] phi; 			phi = NULL;
	  }
	  if(diff != NULL) { // Note: diff is freed somewhere else
		  delete [] diff; 			diff = NULL;
	  }
	  if(grad_phi != NULL) {
		  delete [] grad_phi; 		grad_phi = NULL;
	  }
	  if(rhophi != NULL) {
		  delete [] rhophi; 		rhophi = NULL;
	  }
	  if(grad_rhophi != NULL) {
		  delete [] grad_rhophi; 	grad_rhophi = NULL;
	  }
	  if(dpress_dphi != NULL) {
		  delete [] dpress_dphi; 	dpress_dphi = NULL;
	  }
  }
};


class UgpWithCv2Op: public UgpWithCv2 {

public:
  // data registration for transport variables ...
  vector<ScalarTranspEq> scalarTranspEqVector;                        ///< vector of scalar transport equations
  typedef vector<ScalarTranspEq>::iterator ScalarTranspEqIterator;    ///< iterator for vector of scalar transport equations


  /**
   *    noc00 -> diagonal element of the icv0 row
   *    noc11 -> diagonal element of the icv1 row
   *    noc01 -> flux contribution of icv1 to the icv0 balance
   *    noc10 -> flux contribution of icv0 to the icv1 balance
   */
  inline void getImplDependencyIndex(int &noc00, int &noc01, int &noc11, int &noc10, int icv0, int icv1) {
    noc00 = nbocv_i[icv0];
    noc01 = noc00;
    while (nbocv_v[++noc01] != icv1);

    if (icv1 < ncv) {
      noc11 = nbocv_i[icv1];
      noc10 = noc11;
      while (nbocv_v[++noc10] != icv0);
    }
  }


  void calcCv2Grad(double (*gradPhi)[3], const double *phi,
        const int limiter, const double *refValue, const double epsilon)
  {
    // calculate gradient in ghost level-0 as well, hence no update is necessary!
    //calcCv2GradLimited(gradPhi, phi, refValue, epsilon);
    calcCv2ScalGrad(gradPhi, phi);

    switch (limiter)
    {
    case BARTH_JESPERSEN_MOD:
      limitCv2GradBJ(gradPhi, phi, refValue, epsilon);
      break;

    default:
      break;
    }
  }

  void calcCv2Grad(double (*gradPhi)[3][3], const double (*phi)[3],
          const int limiter, const double *refValue, const double epsilon)
  {
    // calculate gradient in ghost level-0 as well, hence no update is necessary!
    //calcCv2GradLimited(gradPhi, phi, refValue, epsilon);
    calcCv2VecGrad(gradPhi, phi);

    switch (limiter)
    {
    case BARTH_JESPERSEN_MOD:
      limitCv2GradBJ(gradPhi, phi, refValue, epsilon);
      break;

    default:
     break;
    }
  }




  /**
   *   populates the gradient in cv range [0..ncv-1]... i.e. there is no reduction
   *   to ghost cvs, and there is nothing set in fake data (if present)...
   *
   *   the gradient requires a scaling, related to the square of the
   *   undivided difference at which a nbr's influence should be reduced
   */
  void calcCv2GradLimited(double(*dphidxi)[3], const double * phi, const double * eps, const double factor)
  {
    FOR_ICV_G
    {
      // build the local system
      // [A]{grad}={B}
      double A[3][3];
      FOR_I3
        FOR_J3
          A[i][j] = 0.0;
      FOR_I3
        A[i][i] = 1.0;

      double B[3] = {0.0, 0.0, 0.0};

      // build the system: internal part first...
      int noc_f = nbocv_all_i[icv];
      int noc_l = nbocv_all_i[icv+1]-1;
      for (int noc = noc_f+1; noc<=noc_l; ++noc) {
        int icv_nbr = nbocv_all_v[noc];
        double dphi = phi[icv_nbr]-phi[icv];
        // this weight is always less than 1 and gets smaller for larger dphi
        double eps2 = pow(eps[icv]*factor, 2.0);
        double weight = eps2/(eps2+dphi*dphi);
        FOR_I3
          B[i] += nbocv_all_grad_coeff[noc][i]*weight*dphi;
        FOR_J3 {
          double dx = x_cv[icv_nbr][j]-x_cv[icv][j];
          FOR_I3
            A[i][j] -= nbocv_all_grad_coeff[noc][i]*(1.0-weight)*dx;
        }
      }

      double inv_denom = 1.0/(A[0][0]*A[1][1]*A[2][2]+A[0][1]*A[1][2]*A[2][0]
                             +A[1][0]*A[2][1]*A[0][2]-A[0][0]*A[1][2]*A[2][1]
                             -A[1][1]*A[2][0]*A[0][2]-A[2][2]*A[0][1]*A[1][0]);

      dphidxi[icv][0] = inv_denom*(B[0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
                                  +B[1]*(A[2][1]*A[0][2]-A[2][2]*A[0][1])
                                  +B[2]*(A[0][1]*A[1][2]-A[0][2]*A[1][1]));

      dphidxi[icv][1] = inv_denom*(B[0]*(A[1][2]*A[2][0]-A[1][0]*A[2][2])
                                  +B[1]*(A[2][2]*A[0][0]-A[2][0]*A[0][2])
                                  +B[2]*(A[0][2]*A[1][0]-A[0][0]*A[1][2]));

      dphidxi[icv][2] = inv_denom*(B[0]*(A[1][0]*A[2][1]-A[1][1]*A[2][0])
                                  +B[1]*(A[2][0]*A[0][1]-A[2][1]*A[0][0])
                                  +B[2]*(A[0][0]*A[1][1]-A[0][1]*A[1][0]));
    }
  }


  /**
   *   populates the gradient in cv range [0..ncv-1]... i.e. there is no reduction
   *   to ghost cvs, and there is nothing set in fake data (if present)...
   *
   *   the gradient requires a scaling, related to the square of the
   *   undivided difference at which a nbr's influence should be reduced
   *   double eps2 = 1.0E-2; // now a pssed argument to this routine...
   */
  void calcCv2GradLimited(double(*dudxi)[3][3], const double(*u)[3], const double * eps, const double factor)
  {
    // we can skip the diagonal here -- it is not set...
    FOR_ICV_G
    {
      // build the local system
      // [A]{grad}={B}
      double A[3][3];
      FOR_I3
        FOR_J3
          A[i][j] = 0.0;
      FOR_I3
        A[i][i] = 1.0;

      double B[3][3];
      FOR_I3
        FOR_J3
          B[i][j] = 0.0;

      // build the system...
      int noc_f = nbocv_all_i[icv];
      int noc_l = nbocv_all_i[icv+1]-1;
      for (int noc = noc_f+1; noc<=noc_l; ++noc) {
        int icv_nbr = nbocv_all_v[noc];
        double du[3];
        FOR_I3
          du[i] = u[icv_nbr][i]-u[icv][i];
        // this weight is always less than 1 and gets smaller for larger du. Note also that we
        // weight all velocity coeff's the same based on comparison of the sum of the squares
        // of all components...
        double eps2 = pow(eps[icv]*factor, 2.0);
        double weight = eps2/(eps2+du[0]*du[0]+du[1]*du[1]+du[2]*du[2]);
        FOR_I3
          FOR_J3
            B[i][j] += nbocv_all_grad_coeff[noc][j]*weight*du[i];
        FOR_J3 {
          double dx = x_cv[icv_nbr][j]-x_cv[icv][j];
          FOR_I3
            A[i][j] -= nbocv_all_grad_coeff[noc][i]*(1.0-weight)*dx;
        }
      }

      double inv_denom = 1.0/(A[0][0]*A[1][1]*A[2][2]+A[0][1]*A[1][2]*A[2][0]
                             +A[1][0]*A[2][1]*A[0][2]-A[0][0]*A[1][2]*A[2][1]
                             -A[1][1]*A[2][0]*A[0][2]-A[2][2]*A[0][1]*A[1][0]);

      FOR_I3
      {
        dudxi[icv][i][0] = inv_denom*(B[i][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
                                     +B[i][1]*(A[2][1]*A[0][2]-A[2][2]*A[0][1])
                                     +B[i][2]*(A[0][1]*A[1][2]-A[0][2]*A[1][1]));

        dudxi[icv][i][1] = inv_denom*(B[i][0]*(A[1][2]*A[2][0]-A[1][0]*A[2][2])
                                     +B[i][1]*(A[2][2]*A[0][0]-A[2][0]*A[0][2])
                                     +B[i][2]*(A[0][2]*A[1][0]-A[0][0]*A[1][2]));

        dudxi[icv][i][2] = inv_denom*(B[i][0]*(A[1][0]*A[2][1]-A[1][1]*A[2][0])
                                     +B[i][1]*(A[2][0]*A[0][1]-A[2][1]*A[0][0])
                                     +B[i][2]*(A[0][0]*A[1][1]-A[0][1]*A[1][0]));
      }
    }
  }


  
  /**
   * apply smooth barth jespersen limiters
   */
  void limitCv2GradBJ(double (*grad_p)[3], const double *phi, const double *refValue, const double epsilonSDWLS)
  {
    for (int icv = 0; icv < ncv_g; icv++)
    {
      double phiMax = phi[icv];
      double phiMin = phiMax;

      int noc_f = nbocv_all_i[icv];
      int noc_l = nbocv_all_i[icv + 1] - 1;
      // skip the diagonal in this loop...
      for (int noc = noc_f + 1; noc <= noc_l; noc++)
      {
        int icv_nbr = nbocv_all_v[noc];
        phiMax = max(phiMax, phi[icv_nbr]);
        phiMin = min(phiMin, phi[icv_nbr]);
      }

      double alfa = 1.0;

      if ((phiMax-phiMin) > 1.0e-12)
      {
        int foc_f = faocv_i[icv];
        int foc_l = faocv_i[icv + 1] - 1;
        for (int foc = foc_f; foc <= foc_l; foc++)
        {
          int ifa = faocv_v[foc];

          //if (ifa >= nfa_b)
          {
            double r0[3] = {0.0, 0.0, 0.0};
            vecMinVec3d(r0, x_fa[ifa], x_cv[icv]);
            double phifa = phi[icv] + vecDotVec3d(r0, grad_p[icv]);
#ifdef BJ
            if ((phifa - phi[icv]) > 0.0)               alfa = min(alfa, (phiMax-phi[icv])/(phifa-phi[icv]));
            else if ((phifa - phi[icv]) < 0.0)          alfa = min(alfa, (phiMin-phi[icv])/(phifa-phi[icv]));
#else
            double dp;
            double eps2 = max(epsilonSDWLS*fabs(refValue[icv]), SDWLSClippingScalar);
            eps2 = pow(eps2, 2.0);
            double dm = phifa - phi[icv];

            if   ((phifa - phi[icv]) >= 0.0)    dp = phiMax - phi[icv];
            else                                dp = phiMin - phi[icv];

            alfa = min(alfa, (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2));
#endif
          }
        }
      }
      else alfa = 0.0;

      for (int i=0; i<3; i++)
        grad_p[icv][i] *= alfa;
    }
  }

  /**
   * apply smooth barth jespersen limiters
   */
  void limitCv2GradBJ(double (*grad_p)[3][3], const double (*phi)[3], const double *refValue, const double epsilonSDWLS)
  {
    for (int icv = 0; icv < ncv_g; icv++)
    {
      for (int i=0; i<3; i++)
      {
        double phiMax = phi[icv][i];
        double phiMin = phiMax;

        int noc_f = nbocv_all_i[icv];
        int noc_l = nbocv_all_i[icv + 1] - 1;
        // skip the diagonal in this loop...
        for (int noc = noc_f + 1; noc <= noc_l; noc++)
        {
          int icv_nbr = nbocv_all_v[noc];
          phiMax = max(phiMax, phi[icv_nbr][i]);
          phiMin = min(phiMin, phi[icv_nbr][i]);
        }

        double alfa = 1.0;

        if ((phiMax-phiMin) > 1.0e-8)
        {
          int foc_f = faocv_i[icv];
          int foc_l = faocv_i[icv + 1] - 1;
          for (int foc = foc_f; foc <= foc_l; foc++)
          {
            int ifa = faocv_v[foc];

            //if (ifa >= nfa_b)
            {
              double r0[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(r0, x_fa[ifa], x_cv[icv]);
              double phifa = phi[icv][i] + vecDotVec3d(r0, grad_p[icv][i]);

#ifdef BJ
              if ((phifa - phi[icv][i]) > 0.0)               alfa = min(alfa, (phiMax-phi[icv][i])/(phifa-phi[icv][i]));
              else if ((phifa - phi[icv][i]) < 0.0)          alfa = min(alfa, (phiMin-phi[icv][i])/(phifa-phi[icv][i]));
#else
              double dp;
              double eps2 = pow(epsilonSDWLS*refValue[icv], 2.0);
              double dm = phifa - phi[icv][i];

              if ((phifa - phi[icv][i]) >= 0.0)   dp = phiMax - phi[icv][i];
              else                                dp = phiMin - phi[icv][i];

              alfa = min(alfa, (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2));
#endif
            }
          }
        }
        else alfa = 0.0;

        for (int j=0; j<3; j++)
          grad_p[icv][i][j] *= alfa;
      }
    }
  }

  //
  // transport scalar registration
  //
  ScalarTranspEq * registerScalarTransport(const char * name, int datatype)
  {
    // check that the name does not conflict with any other registered data...
    for (ScalarTranspEqIterator data = scalarTranspEqVector.begin(); data < scalarTranspEqVector.end(); data++)
    {
      if (strcmp(name, data->getName()) == 0)
      {
        if (mpi_rank == 0)
          cerr << "Error: transport scalar already registered with name: " << name << endl;
        throw(-1);
      }
    }

    //-----------------------------------------------------------------------------------------------------------------
    // push it into vector,
    // ATTENTION: Reallocations of a stl vector invalidate all previously obtained iterators, references and pointers.
    // see: http://www.cplusplus.com/reference/stl/vector/push_back/
    scalarTranspEqVector.push_back(ScalarTranspEq(name, datatype));

    //-----------------------------------------------------------------------------------------------------------------
    // therefore, we need to re-connect the doubleScalarListPointer (phi) to the pointer sitting in ScalarList (**ptr)
    for (int i = 0; i < scalarTranspEqVector.size()-1; i++)
      for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data != doubleScalarList.end(); data++)
        if (strcmp(scalarTranspEqVector[i].getName(), data->getName()) == 0)
          data->ptr = &scalarTranspEqVector[i].phi; //data->connectPointers(scalarTranspEqVector[i].phi);


    ScalarTranspEq *transScal = &(scalarTranspEqVector.back());

    registerScalar(transScal->phi, name, datatype);

    return transScal;
  }

  ScalarTranspEq * getScalarTransportData(const string &name)
  {
    return getScalarTransportData(name.c_str());
  }

  ScalarTranspEq * getScalarTransportData(const char * name)
  {
    for (ScalarTranspEqIterator data = scalarTranspEqVector.begin(); data < scalarTranspEqVector.end(); data++)
    {
      if (strcmp(name, data->getName()) == 0)
        return (&(*data));
    }
    return (NULL);
  }

  int getScalarTransportIndex(const string &name)
  {
    return getScalarTransportIndex(name.c_str());
  }

  int getScalarTransportIndex(const char * name)
  {
    for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
    {
      if (strcmp(name, scalarTranspEqVector[iScal].getName()) == 0)
        return iScal;
    }
    return (-1);
  }

    int getScalarTransportCoupledIndex(const string &name)
  { 
    return getScalarTransportIndex(name.c_str());
  }     
    
  int getScalarTransportCoupledIndex(const char * name)
  {
    int iScalCoupled = 0;
    for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
    {
      if (scalarTranspEqVector[iScal].coupling == "COUPLED")
      {
        if (strcmp(name, scalarTranspEqVector[iScal].getName()) == 0)
          return iScalCoupled;
        iScalCoupled++;
      } 
    }
    return (-1);
  }    

  /**
   * scalar transport equation -> calc residual
   */
  void calcResidual(double &res_ave, double &res_max, double *phi, double *A, double *rhs)
  {
    double this_res;
    res_max = -1.0;
    res_ave =  0.0;
    for (int icv = 0; icv < ncv; icv++)
    {
      this_res = rhs[icv];

      int noc_f = nbocv_i[icv]; // diagonal
      int noc_l = nbocv_i[icv+1]-1;

      this_res -= A[noc_f]*phi[icv];

      for (int noc = noc_f+1; noc <= noc_l; noc++)
      {
        int icv_nbr = nbocv_v[noc];
        this_res -= A[noc]*phi[icv_nbr];
      }

      res_ave += fabs(this_res);
      res_max = max(res_max, this_res);
    }
  }



  //##############################################################################
  //
  //      home cooked math functions used for BCGSTABR5,
  //      ludeco ... LU decomposition
  //      lusolv ... LU back substitution
  //
  //##############################################################################
  void ludeco(double(*A)[5], int order) {
    for (int jc = 1; jc < order; jc++)
      A[0][jc] /= A[0][0];

    int jrjc = 0;

    for (;;) {
      jrjc++;
      int jrjcm1 = jrjc - 1;
      int jrjcp1 = jrjc + 1;
      for (int jr = jrjc; jr < order; jr++) {
        double sum = A[jr][jrjc];
        for (int jm = 0; jm <= jrjcm1; jm++)
          sum -= A[jr][jm] * A[jm][jrjc];
        A[jr][jrjc] = sum;
      }
      if (jrjc == (order - 1))
        return;

      for (int jc = jrjcp1; jc < order; jc++) {
        double sum = A[jrjc][jc];

        for (int jm = 0; jm <= jrjcm1; jm++)
          sum -= A[jrjc][jm] * A[jm][jc];

        A[jrjc][jc] = sum / A[jrjc][jrjc];
      }
    }
  }

  void lusolv(double(*A)[5], double *c, int order) {
    //              ...First l(inv)*b...
    c[0] = c[0] / A[0][0];
    for (int jr = 1; jr < order; jr++) {
      int jrm1 = jr - 1;
      double sum = c[jr];
      for (int jm = 0; jm <= jrm1; jm++)
        sum -= A[jr][jm] * c[jm];
      c[jr] = sum / A[jr][jr];
    }

    //             ...Next u(inv) of l(inv)*b...
    for (int jrjr = 1; jrjr < order; jrjr++) {
      int jr = (order - 1) - jrjr;
      int jrp1 = jr + 1;
      double sum = c[jr];
      for (int jmjm = jrp1; jmjm < order; jmjm++) {
        int jm = (order - 1) - jmjm + jrp1;
        sum -= A[jr][jm] * c[jm];
      }
      c[jr] = sum;
    }
  }

  /**
   * calculates A[icv][5][5]*phi[icv][5] = res[icv][5]
   * used for solveCvVectorR5Bcgstab()
   */
  void matTimesVecOverCVs(double(*res)[5], double(*Ap)[5][5], double(*phi)[5])
  {
    for (int icv = 0; icv < ncv; icv++)
    {
      for (int i = 0; i < 5; i++)
        res[icv][i] = 0.0;

      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;
      for (int noc = noc_f; noc <= noc_l; noc++)      // cv diagonal + neighbors...
      {
        const int icv_nbr = nbocv_v[noc];

        for (int i = 0; i < 5; i++)
          for (int j = 0; j < 5; j++)
            res[icv][i] += Ap[noc][i][j] * phi[icv_nbr][j];
      }
    }
  }

  /**
   *  update vector phi[icv][5], where [1][2][3] is the velocity vector that needs REPLACE_ROTATE_DATA
   *  used for solveCvVectorR5Bcgstab()
   */
  void UpdateCvDataStateVec(double(*phi)[5]) {
    static double *scal = new double[ncv_g];
    static double (*vec)[3] = new double[ncv_g][3];

    for (int icv = 0; icv < ncv_g; icv++)
      scal[icv] = phi[icv][0];
    updateCvDataG1(scal, REPLACE_DATA);
    for (int icv = 0; icv < ncv_g; icv++)
      phi[icv][0] = scal[icv];

    for (int icv = 0; icv < ncv_g; icv++)
      scal[icv] = phi[icv][4];
    updateCvDataG1(scal, REPLACE_DATA);
    for (int icv = 0; icv < ncv_g; icv++)
      phi[icv][4] = scal[icv];

    for (int icv = 0; icv < ncv_g; icv++)
      for (int i = 0; i < 3; i++)
        vec[icv][i] = phi[icv][i + 1];
    updateCvDataG1(vec, REPLACE_ROTATE_DATA);
    for (int icv = 0; icv < ncv_g; icv++)
      for (int i = 0; i < 3; i++)
        phi[icv][i + 1] = vec[icv][i];
  }

  /**
   * bcgstab
   */
  int solveCvScalarBcgstab(double * phi, double * Ap, double * rhs,
      const double zero, const double zeroRel, const int maxiter, char *scalName);


  /**
   * bcgstab for coupled 5x5 matrices. uses diagonal pre-conditioning using ludeco and lusolv (lower-upper deco and solve)
   */
  int solveCvVectorR5Bcgstab(double(*phi)[5], double(*Ap)[5][5], double(*rhs)[5],
      const double zeroAbs, const double zeroRel, const int maxiter, char *scalarName);


};





#endif

