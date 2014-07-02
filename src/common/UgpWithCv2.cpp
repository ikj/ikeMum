#include "UgpWithCv2.h"
#include <math.h>



void UgpWithCv2::addGhostAndFakeCvs() {

  if (mpi_rank == 0)
    cout << " > addGhostAndFakeCvs()" << endl;
  
  // this routine requires ordered faces.
  
  // we will use the standard cvora to distinguish fake cvs on each processor...
  assert( cvora[mpi_rank+1] - cvora[mpi_rank] == ncv );
  
  // put a global index into cv_flag that INCLUDES the fake cvs on each processor...
  int * cvora_fake = NULL;
  buildXora(cvora_fake,ncv+nfa_b);
  FOR_ICV cv_flag[icv] = icv + cvora_fake[mpi_rank];
    
  // now copy this global index into fa_flag along inter-processor boundaries...
  // 1. boundary faces get the global fake index that assumes no intermediate
  // ghost cvs for now (this gets corrected below -- at present we do not know
  // how many ghost cvs there will be)...
  for (int ifa = 0; ifa < nfa_b; ++ifa) {
    fa_flag[ifa] = ifa + ncv + cvora_fake[mpi_rank];
    assert( (cvofa[ifa][0] >= 0)&&(cvofa[ifa][0] < ncv) );
    assert( cvofa[ifa][1] == -1 );
  }
  // 2. inter-processor/periodic faces...
  for (int ifa = nfa_b; ifa < nfa_bpi; ++ifa) {
    // put our local cv in -- this gets flipped with the updateFaData below...
    fa_flag[ifa] = cv_flag[cvofa[ifa][0]];
    assert( (cvofa[ifa][0] >= 0)&&(cvofa[ifa][0] < ncv) );
    assert( cvofa[ifa][1] == -1 );
  }
  // 3. internal faces...
  for (int ifa = nfa_bpi; ifa < nfa; ++ifa) {
    fa_flag[ifa] = -1;
    assert( (cvofa[ifa][0] >= 0)&&(cvofa[ifa][0] < ncv) );
    assert( (cvofa[ifa][1] >= 0)&&(cvofa[ifa][1] < ncv) );
  }
  // exchange global icv0 indices between inter-processor/periodic face pairs...
  updateFaData(fa_flag,REPLACE_DATA);
    
  // check faces are not overwritten...
  for (int ifa = 0; ifa < nfa_bpi; ++ifa)   assert(fa_flag[ifa] >= 0);
  for (int ifa = nfa_bpi; ifa < nfa; ++ifa) assert(fa_flag[ifa] == -1);
    
  // we are also going to store the periodic transform (if any) associated
  // with faces...
  int * fa_bits = new int[nfa_bpi];
  for (int ifa = 0; ifa < nfa_bpi; ++ifa) fa_bits[ifa] = -1;
    
  // start by thinking about the pack-side of the cvPrcommList. On the pack side, the
  // ghost cvs are owned. 

  // pack stuff...
  for (int iter = 0; iter < 2; ++iter) {
      
    // cycle through the ranges of each face's prcomm...
    for (list<Prcomm>::iterator fa_prcomm = facePrcommList.begin(); 
	 fa_prcomm != facePrcommList.end(); ++fa_prcomm) {
	
      // nbr rank stuff...
      int nbr_rank = fa_prcomm->getNbrRank();
      fa_prcomm->npack_v = 0;
	
      // cycle through the ranges... each range is associated with a unique transformation...
      for (list<Range>::iterator fa_range = fa_prcomm->packRangeList.begin(); 
	   fa_range != fa_prcomm->packRangeList.end(); ++fa_range) {
	  
	// the periodicity associated with this face is...
	int periodic_bits = fa_range->getBits(); // may be zero (INTERNAL) or a single specific bit indicating a FULL or PERIODIC range 
	assert( periodic_bits >= 0 );
	int inverse_periodic_bits = getInversePeriodicBits(periodic_bits);
	  
	// the faces in this range are...
	int fa_index_f = fa_range->getIndexFirst();
	int fa_index_l = fa_range->getIndexLast();
	for (int fa_index = fa_index_f; fa_index <= fa_index_l; ++fa_index) {

	  // this face...
	  int ifa = fa_prcomm->packIndexVec[fa_index];
	  assert( (ifa >= nfa_b)&&(ifa < nfa_bpi) );
	  
	  // grab the internal cv...
	  int icv0 = cvofa[ifa][0];
	  assert( (icv0 >= 0)&&(icv0 < ncv) );

	  if (iter == 0) {
	    // on the first iter, put any transformation into the fa_bits... 
	    assert( fa_bits[ifa] == -1 );
	    // put the local transform - i.e. periodic bits. This gets switched to
	    // inverse periodic bits by doing a updateFace on fa_bits after iter==1
	    fa_bits[ifa] = periodic_bits; 
	    // also count how many cvs we need to send...
	    fa_prcomm->npack_v += 1 + 2*(1+faocv_i[icv0+1]-faocv_i[icv0]); // count plus 2 ints per nbr...
	  }
	  else {
	    // now the fa_bits should contain the inverse periodic bits - verify...
	    assert( fa_bits[ifa] == inverse_periodic_bits ); 
	    // on the second iter, we actually pack...
	    fa_prcomm->packBufferInt[fa_prcomm->npack_v++] = 1+faocv_i[icv0+1]-faocv_i[icv0];
	    // the first thing we send is icv0, the diagonal...
	    fa_prcomm->packBufferInt[fa_prcomm->npack_v++] = cv_flag[icv0];
	    fa_prcomm->packBufferInt[fa_prcomm->npack_v++] = 0; // no additional transformation
	    // now cycle through icv0's faces and send remaining cv's - it is important to send them
	    // in the same order as the faces because this is the basis for the nbocv_all CSR structure
	    // that we want to reporoduce in ghost cvs...
	    int foc_f = faocv_i[icv0];
	    int foc_l = faocv_i[icv0+1]-1;
	    for (int foc = foc_f; foc <= foc_l; ++foc) {
	      int ifa_nbr = faocv_v[foc];
	      if ( ifa_nbr >= nfa_bpi ) {
		// ------------------------------------------------------------------
		// this is an internal face, so send the global cv number and no
		// periodic transform...
		// ------------------------------------------------------------------
		assert( fa_flag[ifa_nbr] == -1 );
		int icv_nbr;
		if (cvofa[ifa_nbr][0] == icv0) {
		  icv_nbr = cvofa[ifa_nbr][1];
		}
		else {
		  assert( cvofa[ifa_nbr][1] == icv0 );
		  icv_nbr = cvofa[ifa_nbr][0];
		}
		assert( (icv_nbr >= 0)&&(icv_nbr < ncv) );
		fa_prcomm->packBufferInt[fa_prcomm->npack_v++] = cv_flag[icv_nbr];
		fa_prcomm->packBufferInt[fa_prcomm->npack_v++] = 0;
	      }
	      else if ( ifa_nbr >= nfa_b ) {
		// ------------------------------------------------------------------
		// this is a ghost face - its global index is in fa_flag... 
		// ------------------------------------------------------------------
		assert( fa_flag[ifa_nbr] >= 0 );
		assert( cvofa[ifa_nbr][0] == icv0 );
		fa_prcomm->packBufferInt[fa_prcomm->npack_v++] = fa_flag[ifa_nbr];
		// and it may have a transform associated with it...
		assert( fa_bits[ifa_nbr] >= 0 );
		fa_prcomm->packBufferInt[fa_prcomm->npack_v++] = fa_bits[ifa_nbr];
	      }
	      else {
		// ------------------------------------------------------------------
		// this is a fake face. it should have no associated transform...
		// ------------------------------------------------------------------
		assert( fa_flag[ifa_nbr] >= 0 );
		assert( fa_bits[ifa_nbr] == -1 );
		fa_prcomm->packBufferInt[fa_prcomm->npack_v++] = fa_flag[ifa_nbr];
		fa_prcomm->packBufferInt[fa_prcomm->npack_v++] = 0;
	      }
	    }
	  }
	  
	}
	
      }
      
      // on the first iter, resize the packBufferInt....
      if ( iter == 0 )
	fa_prcomm->ensurePackBufferIntSize(fa_prcomm->npack_v);
      
    }
      
    // this is used to test the periodicity...
    if ( iter == 0 )
      updateFaData(fa_bits,REPLACE_DATA);
    
  }

  delete[] fa_bits;
    
  exchangePrcommBufferIntV( facePrcommList );
    
  // ================================================
  // fill a local vector that contains the following:
  // global cv number
  // global periodicity (in pack form)
  // local cv number
  // ================================================
    
  vector<IntQuad> quadVec;
  
  for (int iter = 0; iter < 2; ++iter) {
    
    // the index in the quadVec...
    int ii = 0;
      
    // cycle through the ranges of each face's prcomm...
    for (list<Prcomm>::iterator fa_prcomm = facePrcommList.begin(); 
	 fa_prcomm != facePrcommList.end(); ++fa_prcomm) {
	
      // nbr rank stuff...
      int nbr_rank = fa_prcomm->getNbrRank();
      int nunpack_v_check = 0;
	
      // cycle through the ranges... each range is associated with a unique transformation...
      for (list<Range>::iterator fa_range = fa_prcomm->unpackRangeList.begin(); 
	   fa_range != fa_prcomm->unpackRangeList.end(); ++fa_range) {
	  
	// the periodicity associated with this face is...
	int periodic_bits = fa_range->getBits(); // may be zero (INTERNAL) or a single specific bit indicating a FULL or PERIODIC range 
	int inverse_periodic_bits = getInversePeriodicBits(periodic_bits);
	
	// the faces in this range are...
	int fa_index_f = fa_range->getIndexFirst();
	int fa_index_l = fa_range->getIndexLast();
	for (int fa_index = fa_index_f; fa_index <= fa_index_l; ++fa_index) {
	  // this face...
	  int ifa = fa_prcomm->unpackIndexVec[fa_index];
	  assert( (ifa >= nfa_b)&&(ifa < nfa_bpi) );
	  // the number of nbs in the unpack buffer...
	  int nnoc = fa_prcomm->unpackBufferInt[nunpack_v_check++];
	  //assert( nnoc == 7 ); // only on hex grids
	  // the first should be the ghost cv...
	  if (iter == 0) {
	    ii += nnoc;
	    nunpack_v_check += 2*nnoc;
	  }
	  else {
	    for (int i = 0; i < nnoc; ++i) {
	      quadVec[ii].i1 = fa_prcomm->unpackBufferInt[nunpack_v_check++]; // the global icv
	      quadVec[ii].i2 = fa_prcomm->unpackBufferInt[nunpack_v_check++]; // any incremental periodicity
	      assert( quadVec[ii].i2 >= 0 );
	      if (inverse_periodic_bits != 0) {
		if  ( quadVec[ii].i2 == periodic_bits )
		  // this must be a local cv - kill the periodicity...
		  quadVec[ii].i2 = 0;
		else if (quadVec[ii].i2 == inverse_periodic_bits)
		  // For the case of 1-cv thick periodicity, it is possible to have a double 
		  // application of the same transform - use the bit at position 30 (rather than 31, the sign bit)
		  // to indicate a double transform (we can do this because we never have more than
		  // two transforms total, because the nbrs are never more than 2 cvs away)...
		  quadVec[ii].i2 |= (1<<30);
		else
		  quadVec[ii].i2 |= inverse_periodic_bits;
	      }
	      // if this is the FIRST cv, then it is in the first layer of ghosts...
	      if (i == 0)
		quadVec[ii].i3 = 0; // put a zero here so it is first in the sorted list
	      else
		quadVec[ii].i3 = 1;
	      // and store the current index for re-sorting...
	      quadVec[ii].i4 = ii;
	      ++ii;
	    }
	  }
	}
	
      }
	
      assert( nunpack_v_check == fa_prcomm->nunpack_v );
	
    }
    
    if ( iter == 0 )
      quadVec.resize(ii);
    
  }
    
  // recall...
  // i1: global icv
  // i2: periodicity
  // i3: 0 in first layer ghosts, 1 otherwise
  // i4: original index of full vector
  
  std::sort(quadVec.begin(),quadVec.end(),IntQuadCompare123());
  
  // now loop on the sorted quadVec and add a new cv for every unique global index AND transformation
  // do two loops to put the local index into i3, count first and then set...
  
  for (int iter = 0; iter <= 1; ++iter) {
    int n_g = 0;
    int n_gg = 0;
    // n_f is the number of level 1 fakes - we already know this is exactly nfa_b
    int n_ff = 0; // level 2 fakes
    for (int ii = 0; ii < quadVec.size(); ++ii) {
      if ( (ii > 0)&&(quadVec[ii].i1 == quadVec[ii-1].i1)&&(quadVec[ii].i2 == quadVec[ii-1].i2) ) {
	// a simple duplicate...
	if (iter == 1)
	  quadVec[ii].i3 = quadVec[ii-1].i3;
      }
      else {
	// this is a cv with unique data...
	// get the rank based on the global index...
	assert( quadVec[ii].i1 < cvora_fake[mpi_size] );
	int rank = 0; while (cvora_fake[rank+1] <= quadVec[ii].i1) ++rank;
	bool is_fake = quadVec[ii].i1-cvora_fake[rank] >= cvora[rank+1]-cvora[rank];
	// and classify...
	if ( (quadVec[ii].i2 == 0)&&(rank == mpi_rank) ) {
	  // this is a local cv - it should never be a fake...
	  assert( !is_fake );
	  if (iter == 1)
	    quadVec[ii].i3 = quadVec[ii].i1 - cvora_fake[mpi_rank];
	}
	else if (quadVec[ii].i3 == 0) {
	  // should be a level-0 ghost.
	  // this should also never be fake...
	  assert( !is_fake );
	  if (iter == 0)
	    n_g += 1;
	  else
	    quadVec[ii].i3 = ncv_g++;
	}
	else if (is_fake) {
	  // this is a fake cv in the second group (i.e. non-local)...
	  if (iter == 0)
	    n_ff += 1;
	  else
	    quadVec[ii].i3 = ncv_ggff++;
	}
	else {
	  if (iter == 0)
	    n_gg += 1;
	  else
	    quadVec[ii].i3 = ncv_gg++;
	}
      }
    }
    if (iter == 0) {
      // we now have all sub-counts, so 
      ncv_g = ncv;
      ncv_gg = ncv+n_g;
      ncv_ggf = ncv+n_g+n_gg+nfa_b; // not incremented above
      ncv_ggff = ncv_ggf; // local "fakes" are first on the list 
    }
  }
  
  /*
  // we now want to build the paired communicator structures. 
  cout << "ncv     : " << ncv << endl;
  cout << "ncv_g   : " << ncv_g << " " << ncv_g-ncv << endl;
  cout << "ncv_gg  : " << ncv_gg << " " << ncv_gg-ncv_g << endl;
  cout << "ncv_ggf : " << ncv_ggf << " " << ncv_ggf-ncv_gg << " nfa_b: " << nfa_b << endl;
  cout << "ncv_ggff: " << ncv_ggff << " " << ncv_ggff-ncv_ggf << endl;
  */

  // ==================================================================
  // build parallel connectivity...
  // at his point the new ghost cvs are in quadVec, and it is sorted
  // i1: global icv - i.e. who owns the cv
  // i2: incremental periodicity
  // i3: local icv
  // i4: original order index
  // ==================================================================
  
  {
    
    // to build this structure, we need to know the total ghost count on each 
    // processor. This way we can correctly set the local index on the nbr_rank...
    
    int ng = ncv_gg-ncv; // total number of ghosts (both levels)
    int * ngora = new int[mpi_size];
    MPI_Allgather(&ng,1,MPI_INT,ngora,1,MPI_INT,mpi_comm);
    
    // here we reuse the structure used to do the nodes, so we can make use
    // of common periodicity-setting routines, etc. Put this stuff in its own
    // scope to clear the daono_v...
    
    vector<Daono> daono_ggff(ncv_ggff-ncv-nfa_b);
    vector<Daono> daono_gg(ncv_gg-ncv);
    vector<Daono> daono_g(ncv_g-ncv);
    
    // note - they are sorted in rank-order (global icv index) so we can do
    // rank identification quite efficiently. Also, daono_ggff is automagically 
    // sorted (I think - note that I sort anyways, but should return after checking)...
    
    int rank_nbr = 0;
    int iggff = 0;
    int igg = 0;
    int ig = 0;
    for (int ii = 0; ii < quadVec.size(); ++ii) {
      if ( ((ii == 0)||(quadVec[ii].i3 != quadVec[ii-1].i3)) && (quadVec[ii].i3 >= ncv) ) {
	
	// this is a new unique cv that has to come from somewhere...
	// from the global i1 we can get the rank...
	while ( cvora_fake[rank_nbr+1] <= quadVec[ii].i1 ) {
	  ++rank_nbr;
	  assert( rank_nbr < mpi_size );
	}
	// set daono_ggff data structure...
	daono_ggff[iggff].rank_nbr = rank_nbr;
	daono_ggff[iggff].ino      = quadVec[ii].i3; // local index
	
	int icv_nbr = quadVec[ii].i1 - cvora_fake[rank_nbr];
	if (icv_nbr >= cvora[rank_nbr+1]-cvora[rank_nbr]) {
	  // this is a level-0 fake on the nbr, it should be 
	  // a level-1 fake here...
	  assert( (quadVec[ii].i3 >= ncv_ggf)&&(quadVec[ii].i3 < ncv_ggff) );
	  // increment the nbr by the total number of ghosts...
	  icv_nbr += ngora[rank_nbr];
	}
	daono_ggff[iggff].ino_nbr  = icv_nbr;
	daono_ggff[iggff].bits     = quadVec[ii].i2; // any periodicity
	++iggff;

	// ===============================================
	// for daono_gg, only add level-1 and level-2 ghost data...
	// ===============================================
	
	if ( quadVec[ii].i3 < ncv_gg ) {
	  daono_gg[igg].rank_nbr = rank_nbr;
	  daono_gg[igg].ino      = quadVec[ii].i3; // local index
	  assert( icv_nbr < cvora[rank_nbr+1]-cvora[rank_nbr] );
	  daono_gg[igg].ino_nbr  = icv_nbr;
	  daono_gg[igg].bits     = quadVec[ii].i2; // any periodicity
	  ++igg;
	}
	
	// ===============================================
	// for daono_g, only add level-1 ghost data...
	// ===============================================
	
	if ( quadVec[ii].i3 < ncv_g ) {
	  daono_g[ig].rank_nbr = rank_nbr;
	  daono_g[ig].ino      = quadVec[ii].i3; // local index
	  assert( icv_nbr < cvora[rank_nbr+1]-cvora[rank_nbr] );
	  daono_g[ig].ino_nbr  = icv_nbr;
	  daono_g[ig].bits     = quadVec[ii].i2; // any periodicity
	  ++ig;
	}
	
      }
    }
    
    assert( iggff == ncv_ggff-ncv-nfa_b );
    assert( igg == ncv_gg-ncv );
    assert( ig == ncv_g-ncv );
    
    // cleanup...
    delete[] ngora;
    
    // I think this should already be sorted...
    if (daono_ggff.size() > 0) 
      std::sort(daono_ggff.begin(),daono_ggff.end(),DaonoCompare());
    // even if our local daono_ggff is zero, we should call this collectively...
    buildPrcommsAndRanges(daono_ggff,cvPrcommListG1G2F2);
      
    // I think this should already be sorted...
    if (daono_gg.size() > 0) 
      std::sort(daono_gg.begin(),daono_gg.end(),DaonoCompare());
    // even if our local daono_gg is zero, we should call this collectively...
    buildPrcommsAndRanges(daono_gg,cvPrcommListG1G2);
    
    // I think this should already be sorted...
    if (daono_g.size() > 0) 
      std::sort(daono_g.begin(),daono_g.end(),DaonoCompare());
    // even if our local daono_g is zero, we should call this collectively...
    buildPrcommsAndRanges(daono_g,cvPrcommListG1);
      
  }
    
  delete[] cvora_fake;
    
  // ============================================================
  // now build the nbocv_all_i struct.
  // Unlike nbocv_i, this one includes ALL possible cv nbrs (including
  // fake, ghosts layers 0 and 1) and is defines over cvs 0:ncv_g-1,
  // so you can build a ghost gradient locally if you need to.
  // ============================================================
  
  // prep the CSR struct...
  assert( nbocv_all_i == NULL );
  nbocv_all_i = new int[ncv_g+1];
  for (int i = 0; i < ncv_g; ++i)
    nbocv_all_i[i+1] = 0;
  assert( nbocv_all_v == NULL );
  
  // also set the fake cvs to the correct value in cvofa[ifa][1] of boundary faces
  FOR_IFA_B {
    // should be -1...
    assert( cvofa[ifa][1] == -1 );
    // set to the local face index...
    cvofa[ifa][1] = ncv_gg + ifa;
  }
  
  // now that we have i3 set to the LOCAL cv index, return to original order (stored in i4)...
  // note that we could do this faster (i.e. linear), but who cares, and build up 
  // the ghost connectivity structures...
  std::sort(quadVec.begin(),quadVec.end(),IntQuadCompare4());
  
  // resize cv_flag - it is used here -- note we discard all previous values (i.e. the global indexing)...
  resize(cv_flag,ncv_ggff);
  
  // resize registered data, retaining values in the current [0:ncv-1]...
  resizeRegisteredData(CV_DATA,ncv,ncv_ggff);
  
  // do it in 2 iters - count the first, build the second...
  for (int iter = 0; iter < 2; ++iter) {
    
    // -----------------------------
    // local part first...
    // -----------------------------
    
    FOR_ICV {
      // diagonal first...
      if (iter == 0) {
	nbocv_all_i[icv+1] += 1;
      }
      else {
	assert( nbocv_all_v[nbocv_all_i[icv]] == -1 );
	nbocv_all_v[nbocv_all_i[icv]] = icv;
	nbocv_all_i[icv] += 1;
      }
      // get nbrs using faces...
      const int foc_f = faocv_i[icv];
      const int foc_l = faocv_i[icv+1]-1;
      for (int foc = foc_f; foc <= foc_l; ++foc) {
	const int ifa = faocv_v[foc];
	if (iter == 0) {
	  nbocv_all_i[icv+1] += 1;
	}
	else {
	  int icv_nbr;
	  if (cvofa[ifa][0] == icv) {
	    icv_nbr = cvofa[ifa][1];
	  }
	  else {
	    assert(cvofa[ifa][1] == icv);
	    icv_nbr = cvofa[ifa][0];
	  }
	  assert( (icv_nbr >= 0)&&(icv_nbr < ncv_ggff) );
	  assert( nbocv_all_v[nbocv_all_i[icv]] == -1 );
	  nbocv_all_v[nbocv_all_i[icv]] = icv_nbr;
	  nbocv_all_i[icv] += 1;
	}
      }
    }
    
    // -----------------------------
    // ghost part
    // -----------------------------
    
    // the index in the quadVec...
    int ii = 0;
    
    // zero the ghost cv_flag values. They are used to prevent double counting of ghost cvs, 
    // which can exist twice in the unpack data...
    FOR_ICV_GGFF cv_flag[icv] = 2;
    for (int icv = ncv; icv < ncv_g; ++icv)
      cv_flag[icv] = 0;
    
    // cycle through the ranges of each face's prcomm...
    for (list<Prcomm>::iterator fa_prcomm = facePrcommList.begin(); 
	 fa_prcomm != facePrcommList.end(); ++fa_prcomm) {
      
      // nbr rank stuff...
      int nbr_rank = fa_prcomm->getNbrRank();
      int nunpack_v_check = 0;
      
      // cycle through the ranges... each range is associated with a unique transformation...
      for (list<Range>::iterator fa_range = fa_prcomm->unpackRangeList.begin(); 
	   fa_range != fa_prcomm->unpackRangeList.end(); ++fa_range) {
	// the faces in this range are...
	int fa_index_f = fa_range->getIndexFirst();
	int fa_index_l = fa_range->getIndexLast();
	for (int fa_index = fa_index_f; fa_index <= fa_index_l; ++fa_index) {
	  // this face...
	  int ifa = fa_prcomm->unpackIndexVec[fa_index];
	  assert( (ifa >= nfa_b)&&(ifa < nfa_bpi) );
	  // the number of ints in the unpack buffer...
	  int nnoc = fa_prcomm->unpackBufferInt[nunpack_v_check++];
	  // the first should be the ghost cv...
	  if (iter == 0) {
	    int icv = quadVec[ii].i3;
	    assert( (icv >= ncv)&&(icv < ncv_g) );
	    // on the first iter, set the cvofa for this face...
	    assert( cvofa[ifa][1] == -1 );
	    cvofa[ifa][1] = icv;
	    // and count for nbocv_all_i/v...
	    if ( cv_flag[icv] == 0 ) {
	      cv_flag[icv] = 1;
	      assert( nbocv_all_i[icv+1] == 0 );
	      nbocv_all_i[icv+1] = nnoc; // already includes diagonal
	    }
	    else {
	      assert( cv_flag[icv] == 1 );
	      assert( nbocv_all_i[icv+1] == nnoc );
	    }
	    ii += nnoc;
	    nunpack_v_check += 2*nnoc;
	  }
	  else {
	    int icv = quadVec[ii++].i3;
	    assert( (icv >= ncv)&&(icv < ncv_g) );
	    if ( cv_flag[icv] == 0 ) {
	      cv_flag[icv] = 1;
	      // add the diagonal...
	      assert( nbocv_all_v[nbocv_all_i[icv]] == -1 );
	      nbocv_all_v[nbocv_all_i[icv]] = icv;
	      nbocv_all_i[icv] += 1;
	      // the remaining nbrs...
	      for (int i = 1; i < nnoc; ++i) {
		int icv_nbr = quadVec[ii++].i3;
		assert( nbocv_all_v[nbocv_all_i[icv]] == -1 );
		nbocv_all_v[nbocv_all_i[icv]] = icv_nbr;
		nbocv_all_i[icv] += 1;
	      }
	    }
	    else {
	      assert( cv_flag[icv] == 1 );
	      ii += nnoc-1; // ii already incremented by 1 above
	    }
	    nunpack_v_check += 2*nnoc;
	  }
	}
      }
      assert( nunpack_v_check == fa_prcomm->nunpack_v );
    }
    
    if (iter == 0) {
      nbocv_all_i[0] = 0;
      FOR_ICV_G nbocv_all_i[icv+1] += nbocv_all_i[icv];
      nbocv_all_s = nbocv_all_i[ncv_g];
      nbocv_all_v = new int[nbocv_all_s];
      for (int noc = 0; noc < nbocv_all_s; ++noc) nbocv_all_v[noc] = -1;
    } 
    else {
      // return nbocv_all_i...
      for (int icv = ncv_g-1; icv > 0; --icv)
	nbocv_all_i[icv] = nbocv_all_i[icv-1];
      nbocv_all_i[0] = 0;
    }
    
  }

  // ============================================
  // use cv_flag to do a little checking...
  // ============================================

  // 1. check standard cv update
  FOR_ICV_GGFF cv_flag[icv] = 0;
  FOR_ICV cv_flag[icv] = 1;
  updateCvDataG1(cv_flag,ADD_DATA);
  FOR_ICV assert(cv_flag[icv] == 1);
  FOR_ICV_G1_ONLY assert(cv_flag[icv] == 1);
  FOR_ICV_G2_ONLY assert(cv_flag[icv] == 0);
  FOR_ICV_F1_ONLY assert(cv_flag[icv] == 0);
  FOR_ICV_F2_ONLY assert(cv_flag[icv] == 0);
  FOR_IFA_B {
    assert( cv_flag[cvofa[ifa][0]] == 1 );
    assert( cv_flag[cvofa[ifa][1]] == 0 );
  }
  FOR_IFA_NONB {
    assert( cv_flag[cvofa[ifa][0]] == 1 );
    assert( cv_flag[cvofa[ifa][1]] == 1 );
  }
  
  // 2. check G1G2 cv update
  FOR_ICV_GGFF cv_flag[icv] = 0;
  FOR_ICV cv_flag[icv] = 1;
  updateCvDataG1G2(cv_flag,ADD_DATA);
  FOR_ICV assert(cv_flag[icv] == 1);
  FOR_ICV_G1_ONLY assert(cv_flag[icv] == 1);
  FOR_ICV_G2_ONLY assert(cv_flag[icv] == 1);
  FOR_ICV_F1_ONLY assert(cv_flag[icv] == 0);
  FOR_ICV_F2_ONLY assert(cv_flag[icv] == 0);
  FOR_IFA_B {
    assert( cv_flag[cvofa[ifa][0]] == 1 );
    assert( cv_flag[cvofa[ifa][1]] == 0 );
  }
  FOR_IFA_NONB {
    assert( cv_flag[cvofa[ifa][0]] == 1 );
    assert( cv_flag[cvofa[ifa][1]] == 1 );
  }

  // 3. check cv update with G1G2F2 
  FOR_ICV_GGFF cv_flag[icv] = 0;
  FOR_ICV cv_flag[icv] = 1;
  FOR_ICV_F1_ONLY cv_flag[icv] = 1;
  updateCvDataG1G2F2(cv_flag,ADD_DATA);
  FOR_ICV_GGFF assert( cv_flag[icv] == 1 );
  
} // addGhostAndFakeCvs

void UgpWithCv2::addB2AndGhostFaces() {

  // ===========================================================================================
  // given the stuff above, now we add new faces to the G1 cvs. This requires reordering the
  // F2 cvs from their present order as set (somewhat arbitrarily) in the above routine because
  // these cvs must be contiguous by boundary zone so we can reference their corresponding
  // faces using the zone->ifa_i_f:zone->ifa_i_l indices.
  // ===========================================================================================
  
  if (mpi_rank == 0)
    cout << " > addB2AndGhostFaces()" << endl;
  
  // loop on local boundary faces and set the boundary zone index...
  FOR_IFA fa_flag[ifa] = -1;
  int n_bzones = 0;
  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
    if (zone->getKind() == FA_ZONE_BOUNDARY) {
      zone->flag = n_bzones++;
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ++ifa) {
	assert( fa_flag[ifa] == -1 ); // should not yet be touched
	fa_flag[ifa] = zone->flag;
	int icv1 = cvofa[ifa][1];
	assert( (icv1 >= ncv_gg)&&(icv1 < ncv_ggf) );
      }
    }
    else {
      zone->flag = -1;
    }
  }
  // check that all boundary faces were touched...
  FOR_IFA_B assert(fa_flag[ifa] >= 0);
  // and no other faces were touched...
  FOR_IFA_NONB assert(fa_flag[ifa] == -1);

  // use the cv_flag to exchange the boundary zone indicies...
  FOR_ICV_GGFF cv_flag[icv] = -1;
  FOR_IFA_B {
    // icv1 is the level-0 fake... 
    int icv1 = cvofa[ifa][1];
    assert( (icv1 >= ncv_gg)&&(icv1 < ncv_ggf) );
    assert( cv_flag[icv1] == -1 ); // should be untouched
    cv_flag[icv1] = fa_flag[ifa];
  }
  // all level-0 fakes should be set...
  FOR_ICV_F1_ONLY assert(cv_flag[icv] >= 0);  
  // make sure level-1 fakes are still untouched...
  FOR_ICV_F2_ONLY assert(cv_flag[icv] == -1);
  updateCvDataG1G2F2(cv_flag,REPLACE_DATA);
  FOR_ICV_F2_ONLY assert(cv_flag[icv] >= 0);
  
  // the list of cvs is not guaranteed monotonic. XXXXX We could have done this by 
  // passing this fa_zone stuff during the previous routine's sort, but do it here
  // with a partial reordering instead...
  
  int bzone_count[n_bzones+1];
  for (int i = 0; i < n_bzones; ++i)
    bzone_count[i+1] = 0;
  for (int iter = 0; iter < 2; ++iter) {
    FOR_ICV_F2_ONLY {
      int izone = cv_flag[icv];
      assert( (izone >= 0)&&(izone < n_bzones) );
      if (iter == 0) {
	bzone_count[izone+1] += 1;
      }
      else {
	cv_flag[icv] = bzone_count[izone];
	bzone_count[izone] += 1;
      }
    }
    if (iter == 0) {
      // turn into displacement...
      bzone_count[0] = ncv_ggf;
      for (int i = 0; i < n_bzones; ++i)
	bzone_count[i+1] += bzone_count[i];
    }
    else {
      // return to displacement...
      for (int i = n_bzones-1; i > 0; --i)
	bzone_count[i] = bzone_count[i-1];
      bzone_count[0] = ncv_ggf;
    }
  }
  // cv_flag now contains a new order for the cvs in the last range. The zone offsets are
  // also in bzone_count.
  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
    if (zone->getKind() == FA_ZONE_BOUNDARY) {
      // put the indices of this additional set of faces in ifa_i_f and ifa_i_l...
      // here we note that these faces are the first faces after the current nfa and have
      // 1-to-1 correspondence with the g2 cvs, so their offset can be computed easily...
      zone->ifa_i_f = bzone_count[zone->flag]-ncv_ggf+nfa; // turn into a b2 face count
      zone->ifa_i_l = bzone_count[zone->flag+1]-1-ncv_ggf+nfa; // same
      //cout << "rank: " << mpi_rank << " zone: " << zone->getName() << " face range: " << 
      //  zone->ifa_i_f << ":" << zone->ifa_i_l << " nfa: " <<  zone->ifa_i_l-zone->ifa_i_f+1 << endl;
    }
  }
  
  // order for all other cv's remains unchanged...
  FOR_ICV_GGF cv_flag[icv] = icv;

  // check order in cv_flag is single-valued...
  {
    int * cv_check = new int[ncv_ggff];
    FOR_ICV_GGFF cv_check[icv] = 0;
    FOR_ICV_GGFF {
      int icv_new = cv_flag[icv];
      cv_check[icv_new] += 1;
    }
    FOR_ICV_GGFF assert( cv_check[icv] == 1 );
    delete[] cv_check;
  }
  
  // reorder the cvPrcommListG1G2F2...
  reorderPrcomm(cvPrcommListG1G2F2,cv_flag,ncv_ggff);
  
  // and the nbocv_all_i/v struct...
  assert( nbocv_all_s == nbocv_all_i[ncv_g] );
  for (int noc = 0; noc < nbocv_all_s; ++noc) {
    int icv_old = nbocv_all_v[noc];
    assert( (icv_old >= 0)&&(icv_old < ncv_ggff) );
    nbocv_all_v[noc] = cv_flag[icv_old];
  }
  // there is no need to reorder the csr struct, as it only spans the
  // internal and G1 cvs, so its order is unchanged. It is only
  // F2 cvs that need the order change...
  //reorder_csr(nbocv_all_i,nbocv_all_v,cv_flag,ncv_g);

  // check: make sure the first element is still the diagonal...
  FOR_ICV_G assert( nbocv_all_v[nbocv_all_i[icv]] == icv );
  
  // ----------------------------------------------
  // did the cv ordering work out?...
  // fa_flag still has the bzone index in it...
  // ----------------------------------------------
  
  // use the cv_flag to exchange the boundary zone indicies...
  FOR_ICV_GGFF cv_flag[icv] = -1;
  FOR_IFA_B {
    // icv1 is the level-0 fake... 
    int icv1 = cvofa[ifa][1];
    assert( (icv1 >= ncv_gg)&&(icv1 < ncv_ggf) );
    assert( cv_flag[icv1] == -1 ); // should be untouched
    assert( fa_flag[ifa] >= 0 ); // should stil have zone index
    cv_flag[icv1] = fa_flag[ifa];
  }
  // all level-0 fakes should be set...
  FOR_ICV_F1_ONLY assert(cv_flag[icv] >= 0);  
  // make sure level-1 fakes are still untouched...
  FOR_ICV_F2_ONLY assert(cv_flag[icv] == -1);
  updateCvDataG1G2F2(cv_flag,REPLACE_DATA);
  FOR_ICV_F2_ONLY assert(cv_flag[icv] >= 0);

  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
    if (zone->getKind() == FA_ZONE_BOUNDARY) {
      int icv_f = zone->ifa_i_f-nfa+ncv_ggf;
      int icv_l = zone->ifa_i_l-nfa+ncv_ggf;
      for (int icv = icv_f; icv <= icv_l; ++icv) {
	assert( (icv >= ncv_ggf)&&(icv < ncv_ggff) );
	assert( cv_flag[icv] == zone->flag );
      }
    }
  }

  // =============================================================
  // now we need to extend the faocv_i/v struct to include
  // ALL faces surrounding G1 cvs...
  // =============================================================
  
  int nfa_g1g1 = 0;
  int nfa_g1g2 = 0;
  int nfa_g1f2 = 0;
  assert(faocv_s == faocv_i[ncv]);
  int faocv_new_s = faocv_s;
  FOR_ICV_G1_ONLY {
    int noc_f = nbocv_all_i[icv];
    int noc_l = nbocv_all_i[icv+1]-1;
    faocv_new_s += noc_l-noc_f; // do not include diagonal...
    for (int noc = noc_f+1; noc <= noc_l; ++noc) {
      int icv_nbr = nbocv_all_v[noc];
      assert( icv_nbr != icv );
      if (icv_nbr < ncv) {
	// this face is already present...
      }
      else if (icv_nbr < ncv_g) {
	// this is a g1:g1 face. It will eventually be
	// double counted, so only add if icv < icv_nbr...
	if (icv < icv_nbr) nfa_g1g1 += 1;
      }
      else if (icv_nbr < ncv_gg) {
	// this is a g1:g2 face...
	nfa_g1g2 += 1;
      }
      else {
	// must be an f2 nbr...
	assert( (icv_nbr >= ncv_ggf)&&(icv_nbr < ncv_ggff) );
	nfa_g1f2 += 1;
      }
    }
  }

  // the number of g1f2 faces should be the same as the number of g2 cvs...
  assert( nfa_g1f2 == ncv_ggff-ncv_ggf );

  // now extend the faces...
  nfa_b2 = nfa + nfa_g1f2;
  nfa_b2g = nfa_b2 + nfa_g1g1;
  nfa_b2gg = nfa_b2g + nfa_g1g2;

  // and resize face-based stuff...
  
  resize(fa_flag,nfa_b2gg);
  resize(faocv_i,ncv+1,ncv_g+1);
  resize(cvofa,nfa,nfa_b2gg);
  for (int ifa = nfa; ifa < nfa_b2gg; ++ifa)
    cvofa[ifa][0] = cvofa[ifa][1] = -1;
  
  resize(faocv_v,faocv_s,faocv_new_s);
  for (int foc = faocv_s; foc < faocv_new_s; ++foc)
    faocv_v[foc] = -1;

  // check if faocv_i/v and nbocv_all_i/v correspond on
  // internal cvs...
  FOR_ICV {
    int noc_f = nbocv_all_i[icv];
    int noc_l = nbocv_all_i[icv+1]-1;
    int foc_f = faocv_i[icv];
    int foc_l = faocv_i[icv+1]-1;
    assert( (noc_l-noc_f) == (foc_l-foc_f+1) );
    for (int noc = noc_f+1; noc <= noc_l; ++noc) {
      int icv_nbr = nbocv_all_v[noc];
      int foc = foc_f + noc - (noc_f+1);
      int ifa = faocv_v[foc];
      assert( ((cvofa[ifa][0] == icv)&&(cvofa[ifa][1] == icv_nbr)) ||
	      ((cvofa[ifa][1] == icv)&&(cvofa[ifa][0] == icv_nbr)) );
    }
  }
  
  // populate faocv_i/v and cvofa. For boundary B2 faces, we don't want the
  // order of the ghosh cvs determining the order of these. They must be synchronized
  // with F2 cvs...

  int ifa_g1g1 = nfa_b2;
  int ifa_g1g2 = nfa_b2g;
  FOR_ICV_G1_ONLY {
    int noc_f = nbocv_all_i[icv];
    int noc_l = nbocv_all_i[icv+1]-1;
    int foc_f = faocv_i[icv];
    faocv_i[icv+1] = foc_f + noc_l-noc_f;
    for (int noc = noc_f+1; noc <= noc_l; ++noc) {
      int icv_nbr = nbocv_all_v[noc];
      int foc = foc_f + noc - (noc_f+1);
      assert( faocv_v[foc] == -1 );
      if (icv_nbr < ncv) {
	// this face is already present...
	int noc_nbr = nbocv_all_i[icv_nbr]+1; // skip diag
	while (nbocv_all_v[noc_nbr] != icv) ++noc_nbr;
	int foc_nbr = faocv_i[icv_nbr] + noc_nbr - (nbocv_all_i[icv_nbr]+1); 
	int ifa = faocv_v[foc_nbr];
	assert( (ifa >= 0)&&(ifa < nfa) );
	assert( (cvofa[ifa][1] == icv)&&(cvofa[ifa][0] == icv_nbr) ); // should be outward-pointing.
	faocv_v[foc] = ifa;
      }
      else if (icv_nbr < ncv_g) {
	// this is a g1:g1 face. It will be double counted, so to
	// avoid this, the first icv adds...
	if (icv < icv_nbr) {
	  // has not been set...
	  faocv_v[foc] = ifa_g1g1;
	  cvofa[ifa_g1g1][0] = icv;
	  cvofa[ifa_g1g1][1] = icv_nbr;
	  ifa_g1g1++;
	}
	else {
	  int noc_nbr = nbocv_all_i[icv_nbr]+1; // skip diag
	  while (nbocv_all_v[noc_nbr] != icv) ++noc_nbr;
	  int foc_nbr = faocv_i[icv_nbr] + noc_nbr - (nbocv_all_i[icv_nbr]+1); 
	  int ifa = faocv_v[foc_nbr];
	  assert( (cvofa[ifa][1] == icv)&&(cvofa[ifa][0] == icv_nbr) );
	  faocv_v[foc] = ifa;
	}
      }
      else if (icv_nbr < ncv_gg) {
	// this is a g1:g2 face...
	faocv_v[foc] = ifa_g1g2;
	cvofa[ifa_g1g2][0] = icv;
	cvofa[ifa_g1g2][1] = icv_nbr;
	ifa_g1g2++;
      }
      else {
	// must be an f2 nbr...
	assert( (icv_nbr >= ncv_ggf)&&(icv_nbr < ncv_ggff) );
	// as mentioned above, here we DO NOT us the G1 cv ordering to
	// arbitrarily determine the B2 face ordering, because
	// the B2 face must be related to the F2 cv (icv_nbr in this case) 
	// by a constant offset...
	int ifa = icv_nbr-ncv_ggf+nfa;
	assert( (ifa >= nfa)&&(ifa < nfa_b2) );
	faocv_v[foc] = ifa;
	assert( cvofa[ifa][0] == -1 );
	cvofa[ifa][0] = icv;
	assert( cvofa[ifa][1] == -1 );
	cvofa[ifa][1] = icv_nbr;
      }
    }
  }

  // check if faocv_i/v and nbocv_all_i/v correspond on
  // G1 cvs...
  FOR_ICV_G1_ONLY {
    int noc_f = nbocv_all_i[icv];
    int noc_l = nbocv_all_i[icv+1]-1;
    int foc_f = faocv_i[icv];
    int foc_l = faocv_i[icv+1]-1;
    assert( (noc_l-noc_f) == (foc_l-foc_f+1) );
    for (int noc = noc_f+1; noc <= noc_l; ++noc) {
      int icv_nbr = nbocv_all_v[noc];
      int foc = foc_f + noc - (noc_f+1);
      int ifa = faocv_v[foc];
      assert( ((cvofa[ifa][0] == icv)&&(cvofa[ifa][1] == icv_nbr)) ||
	      ((cvofa[ifa][1] == icv)&&(cvofa[ifa][0] == icv_nbr)) );
    }
  }

  // as one last check, check that the face boundary zones match...

  FOR_IFA_ALL fa_flag[ifa] = -1;
  FOR_ICV_ALL cv_flag[icv] = -1;
  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
    if (zone->getKind() == FA_ZONE_BOUNDARY) {
      // the range of B faces is...
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ++ifa) {
	assert( fa_flag[ifa] == -1 );
	fa_flag[ifa] = zone->flag;
	// also set the associated F1 cv to the same...
	int icv1 = cvofa[ifa][1];
	// G1 cvs are related to B faces by a simple offset...
	assert( ifa == icv1 - ncv_gg );
	assert( (icv1 >= ncv_gg)&&(icv1 < ncv_ggf) );
	assert( cv_flag[icv1] == -1 );
	cv_flag[icv1] = zone->flag;
      }
      // the range of B2 faces is...
      for (int ifa = zone->ifa_i_f; ifa <= zone->ifa_i_l; ++ifa) {
	assert( (ifa >= nfa)&&(ifa < nfa_b2) );
	assert( fa_flag[ifa] == -1 );
	fa_flag[ifa] = zone->flag;
	int icv1 = cvofa[ifa][1];
	// G2 cvs are related to B2 faces by a simple offset...
	assert( ifa-nfa == icv1-ncv_ggf );
	assert( cv_flag[icv1] == -1 );
      }
    }
  }

  // check face...
  FOR_IFA_B assert(fa_flag[ifa] >= 0);
  FOR_IFA_NONB assert( fa_flag[ifa] == -1 );
  FOR_IFA_B2_ONLY assert( fa_flag[ifa] >= 0 );
  FOR_IFA_GG_ONLY assert( fa_flag[ifa] == -1 );
  // check cvs...
  FOR_ICV_F1_ONLY assert( cv_flag[icv] >= 0 );
  FOR_ICV_F2_ONLY assert( cv_flag[icv] == -1 );
  // exchange cv2's...
  updateCvDataG1G2F2(cv_flag,REPLACE_DATA);
  // only F2's should have changed...
  FOR_ICV_F2_ONLY assert( cv_flag[icv] >= 0 );
  // the zones of B2 faces should match the F2 zones just exchanged...
  FOR_IFA_B2_ONLY {
    assert( fa_flag[ifa] >= 0 );
    int icv1 = cvofa[ifa][1];
    assert( cv_flag[icv1] >= 0 );
    assert( fa_flag[ifa] == cv_flag[icv1] );
    cv_flag[icv1] = -1; // uniqueness check
  }
  FOR_ICV_F2_ONLY assert( cv_flag[icv] == -1 );

  // everything looks great. Finally, verify that the updateFaData will 
  // access the faces in the proper order. Start by building a global cv index,
  // including fake cvs (which should be globally ordered at the end)...

  int ncv_global;
  MPI_Allreduce(&ncv,&ncv_global,1,MPI_INT,MPI_SUM,mpi_comm);
  assert( ncv_ggf-ncv_gg == nfa_b ); // the number of F1 cvs should match the face b count 
  int count[2] = { ncv, ncv_ggf-ncv_g };
  int offset[2];
  MPI_Scan(count,offset,2,MPI_INT,MPI_SUM,mpi_comm);
  // the active cv offset...
  offset[0] -= count[0];
  // the fake cv offset...
  offset[1] += ncv_global-count[1];
  // now populate cv_flag with the global index and exchange...
  FOR_ICV cv_flag[icv] = icv + offset[0];
  FOR_ICV_G1_ONLY cv_flag[icv] = -1;
  FOR_ICV_G2_ONLY cv_flag[icv] = -1;
  FOR_ICV_F1_ONLY cv_flag[icv] = icv + offset[1];
  FOR_ICV_F2_ONLY cv_flag[icv] = -1;
  updateCvDataG1G2F2(cv_flag,REPLACE_DATA);
  FOR_ICV_ALL assert( cv_flag[icv] >= 0 );

  // the idea behind the updateFaData exchange is to loop on face neighbors of cvs
  // and buffer the face data, then exchange it. To do this in a synchronized way, we need to 
  // know the maximum number of faces any one cv has...
  int my_nfoc_max = 0;
  FOR_ICV my_nfoc_max = max( my_nfoc_max, faocv_i[icv+1]-faocv_i[icv] );
  int nfoc_max;
  MPI_Allreduce(&my_nfoc_max,&nfoc_max,1,MPI_INT,MPI_MAX,mpi_comm);

  int * ibuf = new int[ncv_ggff]; // XXXXX once updateCvDataG1G2F2 does not have fake update, this can be ncv_g in size
  FOR_ICV_ALL ibuf[icv] = -1;
  for (int foc_inc = 0; foc_inc < nfoc_max; ++foc_inc) {
    // here we should be able to cycle through just a subset of these internal
    // cvs that actually correspond to ghost cvs somewhere (and for that matter
    // nfoc_max may be too large) but to simplfy this routine, just do all internal
    // cvs...
    FOR_ICV {
      int foc = faocv_i[icv] + foc_inc;
      if (foc < faocv_i[icv+1]) {
	int ifa = faocv_v[foc];
	if (cvofa[ifa][0] == icv) {
	  ibuf[icv] = cv_flag[cvofa[ifa][1]];
	}
	else {
	  assert( cvofa[ifa][1] == icv );
	  ibuf[icv] = cv_flag[cvofa[ifa][0]];
	}
      }
    }
    updateCvDataG1G2F2(ibuf,REPLACE_DATA);
    FOR_ICV_G1_ONLY {
      int foc = faocv_i[icv] + foc_inc;
      if (foc < faocv_i[icv+1]) {
	int ifa = faocv_v[foc];
	if (cvofa[ifa][0] == icv) {
	  assert( ibuf[icv] == cv_flag[cvofa[ifa][1]] );
	}
	else {
	  assert( cvofa[ifa][1] == icv );
	  assert( ibuf[icv] == cv_flag[cvofa[ifa][0]] );
	}
      }
    }
  }
  delete[] ibuf;

  // as a final step to simplify the application of bcs, put these faces into 
  // one list (vector actually) described by faVec...
  
  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++) {
    if (zone->getKind() == FA_ZONE_BOUNDARY) {
      assert( zone->faVec.size() == 0 );
      zone->faVec.resize(zone->ifa_l-zone->ifa_f+1+zone->ifa_i_l-zone->ifa_i_f+1);
      int index = 0;
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ++ifa)
	zone->faVec[index++] = ifa;
      for (int ifa = zone->ifa_i_f; ifa <= zone->ifa_i_l; ++ifa)
	zone->faVec[index++] = ifa;
      assert( index == zone->faVec.size() );
    }
  }


} // void UgpWithCv2::addBoundaryFaces2()

void UgpWithCv2::calcGeometry() {
  
  if (mpi_rank == 0)
    cout << " > calcGeometry()" << endl;
  
  // ======================================
  // face normal and centroid...
  // normal has area magnitude
  // ======================================

  assert(x_fa == NULL);
  x_fa = new double[nfa_b2gg][3];

  assert(fa_normal == NULL);
  fa_normal = new double[nfa_b2gg][3];

  FOR_IFA {
  
    FOR_I3 x_fa[ifa][i] = 0.0;
    for (int nof = noofa_i[ifa]; nof < noofa_i[ifa+1]; ++nof) {
      int ino = noofa_v[nof];
      FOR_I3 x_fa[ifa][i] += x_no[ino][i];
    }
    FOR_I3 x_fa[ifa][i] /= (double) (noofa_i[ifa+1] - noofa_i[ifa]);
    
    // we can compute the face normal directly with this approx centroid...
    
    FOR_I3 fa_normal[ifa][i] = 0.0;
    int ino2 = noofa_v[noofa_i[ifa+1] - 1]; // last node in the list
    for (int nof = noofa_i[ifa]; nof < noofa_i[ifa+1]; ++nof) {
      int ino1 = ino2;
      ino2 = noofa_v[nof];
      double v1[3];
      double v2[3];
      FOR_I3 {
        v1[i] = x_no[ino1][i] - x_fa[ifa][i];
        v2[i] = x_no[ino2][i] - x_fa[ifa][i];
      }
      fa_normal[ifa][0] += 0.5 * (v1[1] * v2[2] - v1[2] * v2[1]);
      fa_normal[ifa][1] += 0.5 * (v1[2] * v2[0] - v1[0] * v2[2]);
      fa_normal[ifa][2] += 0.5 * (v1[0] * v2[1] - v1[1] * v2[0]);
    }

    double fa_area = sqrt( fa_normal[ifa][0]*fa_normal[ifa][0] + 
			   fa_normal[ifa][1]*fa_normal[ifa][1] + 
			   fa_normal[ifa][2]*fa_normal[ifa][2] );
    double unit_normal[3];
    FOR_I3 unit_normal[i] = fa_normal[ifa][i]/fa_area;
    
    // find true centroid...
    
    bool done = false;
    int iter = 0;
    while (!done) {
      double dx_fa[3] = { 0.0, 0.0, 0.0 };
      double x_fa_tmp[3] = { 0.0, 0.0, 0.0 };
      double area_sum = 0.0;
      int ino2 = noofa_v[noofa_i[ifa+1]-1]; // last node in the list
      for (int nof = noofa_i[ifa]; nof < noofa_i[ifa+1]; ++nof) {
	int ino1 = ino2;
	ino2 = noofa_v[nof];
	double v1[3];
	double v2[3];
	FOR_I3 {
	  v1[i] = x_no[ino1][i] - x_fa[ifa][i];
	  v2[i] = x_no[ino2][i] - x_fa[ifa][i];
	}
	double this_normal[3];
	this_normal[0] = (v1[1] * v2[2] - v1[2] * v2[1]);
	this_normal[1] = (v1[2] * v2[0] - v1[0] * v2[2]);
	this_normal[2] = (v1[0] * v2[1] - v1[1] * v2[0]);
	// I think this is how fluent does it. This guarantees positive 
	// weights, but can converge quite slowly for highly non-coplanar 
	// faces...
	/*
	  double this_area = sqrt( this_normal[0]*this_normal[0] + 
	  this_normal[1]*this_normal[1] + 
	  this_normal[2]*this_normal[2] );
	*/
	// this is the one I prefer...
	double this_area = 
	  this_normal[0]*unit_normal[0] + 
	  this_normal[1]*unit_normal[1] + 
	  this_normal[2]*unit_normal[2];
	FOR_I3 x_fa_tmp[i] += this_area*( x_no[ino1][i] + x_no[ino2][i] + x_fa[ifa][i] );
	FOR_I3 dx_fa[i] += this_area*( v1[i] + v2[i] ); 
	area_sum += this_area;
      }
      double d2 = 0.0;
      FOR_I3 {
	x_fa_tmp[i] /= 3.0*area_sum;
	dx_fa[i] /= 2.0*area_sum;
	x_fa[ifa][i] += dx_fa[i];
	d2 += dx_fa[i]*dx_fa[i];
      }
      // just set to the x_fa_tmp value...
      FOR_I3 x_fa[ifa][i] = x_fa_tmp[i];
      if (++iter > 200) {
	cerr << "Error: face centroid did not converge: d2/area_sum: " << d2/area_sum << endl;
	throw(-1);
      }
      else if (iter > 150) {
	cout << " > face centroid not converging: " << ifa << " " << iter << " " << d2/area_sum << endl;
      }
      done = (d2/area_sum < 1.0E-20);
    }

  }

  // and extend face centroids and normals...

  // now we need to get the centroids and normals over to G1 faces...
  // note that the updateFaData routine is a hack built on updateCvDataG1G2F2, so I
  // would not recommend using it during anything by grid preprocessing...
  
  updateFaDataB2G1G2(x_fa,REPLACE_TRANSLATE_DATA);
  updateFaDataB2G1G2(fa_normal,REPLACE_ROTATE_DATA);
  
  // ======================================
  // x_cv, cv_volume in ghost as well...
  // ======================================
  
  assert(x_cv == NULL);
  x_cv = new double[ncv_ggff][3]; // x_cv available in ALL cvs
  
  assert(cv_volume == NULL);
  cv_volume = new double[ncv_g]; // Note: volume only defined over ghost level-0
  
  FOR_ICV {
    
    // compute an approximate center as the mean of the face x_fa's...
    double x_cv_approx[3];
    FOR_I3 x_cv_approx[i] = 0.0;
    for (int foc = faocv_i[icv]; foc < faocv_i[icv + 1]; ++foc) {
      int ifa = faocv_v[foc];
      FOR_I3 x_cv_approx[i] += x_fa[ifa][i];
    }
    // divide by the number of faces...
    FOR_I3 x_cv_approx[i] /= (double) (faocv_i[icv + 1] - faocv_i[icv]);
    
    // add tet volumes...
    cv_volume[icv] = 0.0;
    FOR_I3 x_cv[icv][i] = 0.0;
    // loop on the faces of this cv...
    for (int foc = faocv_i[icv]; foc < faocv_i[icv + 1]; ++foc) {
      int ifa = faocv_v[foc];
      double v1[3];
      FOR_I3 v1[i] = x_fa[ifa][i] - x_cv_approx[i];
      // check if the face is inward or outward wrt this cv...
      if (cvofa[ifa][0] == icv) {
        // face is outward - loop through edges in forward direction...
        int ino2 = noofa_v[noofa_i[ifa + 1] - 1];
        for (int nof = noofa_i[ifa]; nof < noofa_i[ifa + 1]; ++nof) {
          int ino1 = ino2;
          ino2 = noofa_v[nof];
          double v2[3];
          double v3[3];
          FOR_I3 {
            v2[i] = x_no[ino1][i] - x_cv_approx[i];
            v3[i] = x_no[ino2][i] - x_cv_approx[i];
          }
          // 2 nodes, the face, and the approx cv form a tet...
          double this_volume = 
	    v1[0] * (v2[1] * v3[2] - v2[2] * v3[1]) +
	    v1[1] * (v2[2] * v3[0] - v2[0] * v3[2]) +
	    v1[2] * (v2[0] * v3[1] - v2[1] * v3[0]);
          assert(this_volume > 0.0); // check on the grid ordering/right-handedness
          cv_volume[icv] += this_volume;
          FOR_I3 x_cv[icv][i] += this_volume * (x_cv_approx[i] + x_fa[ifa][i] + x_no[ino1][i] + x_no[ino2][i]);
        }
      }
      else {
        assert(cvofa[ifa][1] == icv);
        // face is outward, loop through edges in backward direction...
        int ino2 = noofa_v[noofa_i[ifa]];
        for (int nof = noofa_i[ifa + 1] - 1; nof >= noofa_i[ifa]; nof--) {
          int ino1 = ino2;
          ino2 = noofa_v[nof];
          double v2[3];
          double v3[3];
          FOR_I3 {
            v2[i] = x_no[ino1][i] - x_cv_approx[i];
            v3[i] = x_no[ino2][i] - x_cv_approx[i];
          }
          double this_volume = 
	    v1[0] * (v2[1] * v3[2] - v2[2] * v3[1]) +
	    v1[1] * (v2[2] * v3[0] - v2[0] * v3[2]) +
	    v1[2] * (v2[0] * v3[1] - v2[1] * v3[0]);
          assert(this_volume > 0.0);
          cv_volume[icv] += this_volume;
          FOR_I3 x_cv[icv][i] += this_volume * (x_cv_approx[i] + x_fa[ifa][i] + x_no[ino1][i] + x_no[ino2][i]);
        }
      }
    }
    // normalize both...
    FOR_I3 x_cv[icv][i] /= 4.0 * cv_volume[icv];
    cv_volume[icv] /= 6.0;
  }

  // update volumes accross processors...
  updateCvDataG1(cv_volume, REPLACE_DATA); // volume in just G1 ghosts?

  // before update x_cv, we must define it in the fake level-0 cvs. here we put the 
  // fake cv collocated with the boundary face...
  FOR_IFA_B {
    int icv1 = cvofa[ifa][1];
    assert( icv1 == ncv_gg+ifa );
    FOR_I3 x_cv[icv1][i] = x_fa[ifa][i];
  }

  // now update into all GGFF cvs...
  updateCvDataG1G2F2(x_cv, REPLACE_TRANSLATE_DATA); // "TRANSLATE" is for any periodicity

  // ===============================
  // geometry checks...
  // ===============================

  double my_volume_sum = 0.0;
  FOR_ICV my_volume_sum += cv_volume[icv];
  double volume_sum;
  MPI_Reduce(&my_volume_sum, &volume_sum, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);
  
  // also check s dot n...
  double my_buf[2];
  my_buf[0] = 1.0;
  my_buf[1] = 0.0;
  FOR_IFA_B {
    // the first faces are boundary faces, and should have icv1 >= ncv_gg...
    assert( (cvofa[ifa][1] >= ncv_gg)&&(cvofa[ifa][1] < ncv_ggf) );
    // infact, these cvs are contiguous like the faces...
    assert( cvofa[ifa][1] == ncv_gg+ifa );
  }
  FOR_IFA_NONB {
    // the remaining faces should have valid cvofa, with a value
    // in the ghost range for the periodic cvs, and in the
    // owned cv range for totally internal cvs...
    assert(cvofa[ifa][1] >= 0);
    if (ifa < nfa_bpi) {
      assert((cvofa[ifa][1] >= ncv) && (cvofa[ifa][1] < ncv_g));
    } 
    else {
      assert(cvofa[ifa][1] < ncv);
    }
    int icv0 = cvofa[ifa][0];
    int icv1 = cvofa[ifa][1];
    double s[3];
    double smag = 0.0;
    double nmag = 0.0;
    FOR_I3 {
      s[i] = x_cv[icv1][i] - x_cv[icv0][i];
      smag += s[i] * s[i];
      nmag += fa_normal[ifa][i] * fa_normal[ifa][i];
    }
    smag = sqrt(smag);
    nmag = sqrt(nmag);
    //cout << "ifa, smag, nmag: " << ifa << " " << smag << " " << nmag << endl;
    double n[3];
    FOR_I3 {
      s[i] /= smag;
      n[i] = fa_normal[ifa][i] / nmag;
    }
    double dp = s[0] * n[0] + s[1] * n[1] + s[2] * n[2];
    assert(dp > 0.0);
    my_buf[0] = min(my_buf[0], dp);
    my_buf[1] = min(my_buf[1], -dp);
  }
  
  double buf[2];
  MPI_Reduce(my_buf, buf, 2, MPI_DOUBLE, MPI_MIN, 0, mpi_comm);
  
  if (mpi_rank == 0)
    cout << " > calcGeometry(), cv volume: " << volume_sum << ", s dot n min/max: " << buf[0] << " " << -buf[1] << endl;
  
  // also compute the min cv distance...
  double my_min_delta = 1.0E+20;
  FOR_IFA {
    const int icv0 = cvofa[ifa][0];
    assert( (icv0 >= 0)&&(icv0 < ncv) );
    const int icv1 = cvofa[ifa][1];
    assert( icv1 >= 0 );
    double d2 = 0.0;
    FOR_I3 {
      double dx = x_cv[icv1][i] - x_cv[icv0][i];
      d2 += dx*dx;
    }
    my_min_delta = min( my_min_delta, sqrt(d2) );
  }
  double min_delta;
  MPI_Reduce(&my_min_delta,&min_delta, 1, MPI_DOUBLE, MPI_MIN, 0, mpi_comm);
  if (mpi_rank == 0)
    cout << " > calcGeometry(), min distance between cvs: " << min_delta << endl;

  // check gcl...

  double (*gcl)[3] = new double[ncv_g][3]; // level-1 ghosts only
  FOR_ICV_G {
    FOR_I3 gcl[icv][i] = 0.0;
    int foc_f = faocv_i[icv]; // available in ghost-1 also 
    int foc_l = faocv_i[icv+1]-1; 
    for (int foc = foc_f; foc <= foc_l; ++foc) {
      int ifa = faocv_v[foc];
      if (cvofa[ifa][0] == icv) { 
	// face normal is outward-pointing, so add to GCL..
	FOR_I3 gcl[icv][i] += fa_normal[ifa][i];
      }
      else {
	assert( cvofa[ifa][1] == icv );
	// face normal is inward-pointing, so subtract...
	FOR_I3 gcl[icv][i] -= fa_normal[ifa][i];
      }
    }
  }
  dumpVectorRange(gcl,ncv_g,"GCL (should be zero)");
  delete[] gcl;
  
  // also check 1 - vol_ratio...
  
  double (*vol_ratio) = new double[ncv_g];
  FOR_ICV_G {
    vol_ratio[icv] = 0.0;
    int foc_f = faocv_i[icv];
    int foc_l = faocv_i[icv+1]-1;
    for (int foc = foc_f; foc <= foc_l; ++foc) {
      int ifa = faocv_v[foc];
      if (cvofa[ifa][0] == icv) {
	// face normal is outward-pointing, so add to GCL..
	FOR_I3  vol_ratio[icv] += (x_fa[ifa][i]-x_cv[icv][i])*fa_normal[ifa][i];
      }
      else {
	assert( cvofa[ifa][1] == icv );
	// face normal is inward-pointing, so subtract...
	FOR_I3  vol_ratio[icv] -= (x_fa[ifa][i]-x_cv[icv][i])*fa_normal[ifa][i];
      }
    }
    vol_ratio[icv] = vol_ratio[icv]/3.0/cv_volume[icv] - 1.0;
  }
  
  dumpScalarRange(vol_ratio,ncv_g,"VOL_RATIO-1 (should be 0 for planar-faced meshes)");
  delete[] vol_ratio;

} // calcGeometry

void UgpWithCv2::uniformHexCheck() {

  // ==================================================================
  // if this is uniform Hex grid with or without boundaries, we 
  // can do this test. It simply takes the average of all nbr cvs and
  // compares to the actual cv for all cvs including level-0 ghosts. The
  // fact that fake cvs have zero volume (i.e. are collocated with the 
  // boundary faces) is accounted for by doubling their associated dx... 
  // ==================================================================

  if (mpi_rank == 0)
    cout << " > uniformHexCheck()...";

  double my_dx_max[3] = { 0.0, 0.0, 0.0 };
  
  FOR_ICV_G {
    int noc_f = nbocv_all_i[icv];
    int noc_l = nbocv_all_i[icv+1]-1;
    //assert( noc_l-noc_f+1 == 7 ); // if pure hex
    assert( nbocv_all_v[noc_f] == icv ); // first is diagonal
    double dx[3] = { 0.0, 0.0, 0.0 };
    for (int noc = noc_f+1; noc <= noc_l; ++noc) {
      int icv_nbr = nbocv_all_v[noc];
      assert( (icv_nbr >= 0)&&(icv_nbr < ncv_ggff) );
      if ( icv_nbr < ncv_gg ) {
	// this is an internal or ghost nbr, so add full delta to dx...
	FOR_I3 dx[i] += x_cv[icv_nbr][i]-x_cv[icv][i];
      }
      else {
	// this is a fake nbr, so double delta and add...
	FOR_I3 dx[i] += 2.0*(x_cv[icv_nbr][i]-x_cv[icv][i]);
      }
    }
    FOR_I3 my_dx_max[i] = max( fabs(dx[i]), my_dx_max[i] );
  }

  double dx_max[3];
  MPI_Reduce(my_dx_max, dx_max, 3, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
  if (mpi_rank == 0)
    cout << " (should be zero on Cartesian hex meshes): " << 
      dx_max[0] << " " << dx_max[1] << " " << dx_max[2] << endl;
  
}

void UgpWithCv2::buildNbocv() {

  //
  // nbocv_i/v (neighbor-of-cv) is the CSR (compact-storage-row) format
  // matrix that holds the immediate cv neighbours. The first element
  // of any row is the diagonal element, i.e.
  //
  // double * A = new double[nbocv_s];
  // for (int icv = 0; icv < ncv; ++icv) {
  //   int noc_f = nbocv_i[icv];
  //   assert( nbocv_v[noc_f] == icv );
  //   A[noc_f] = ... // the diagonal
  //   ...
  //
  // requirements: face data must be ordered.
  //
  // Note that this structure does is defined over locally-owned
  // cv data only, and does NOT contain any reference to ghost data.
  //
  
  if (mpi_rank == 0)
    cout << " > buildNbocv()" << endl;
  
  assert(nbocv_i == NULL);
  assert(nbocv_v == NULL);
  
  nbocv_i = new int[ncv+1];
  FOR_ICV nbocv_i[icv+1] = 0; // used for counts in first iteration below...
  
  for (int icv = 0; icv < ncv; ++icv)          cv_flag[icv] = -1; // locally owned
  for (int icv = ncv; icv < ncv_g; ++icv)      cv_flag[icv] = -2; // ghost level-0
  for (int icv = ncv_g; icv < ncv_ggff; ++icv) cv_flag[icv] = -3; // ghost level-1 and all fake
  
  for (int iter = 0; iter < 2; ++iter) {
    FOR_ICV {
      // diagonal first...
      assert( cv_flag[icv] == -1 );
      cv_flag[icv] = 1;
      if (iter == 0) {
	nbocv_i[icv+1] += 1;
      }
      else {
	nbocv_v[nbocv_i[icv]] = icv;
	nbocv_i[icv] += 1;
      }
      // get nbrs using faces...
      const int foc_f = faocv_i[icv];
      const int foc_l = faocv_i[icv+1]-1;
      for (int foc = foc_f; foc <= foc_l; ++foc) {
	const int ifa = faocv_v[foc];
	int icv_nbr;
	if (cvofa[ifa][0] == icv) {
	  icv_nbr = cvofa[ifa][1];
	}
	else {
	  assert(cvofa[ifa][1] == icv);
	  icv_nbr = cvofa[ifa][0];
	}
	if (ifa < nfa_b) {
	  // this is a boundary face, and should have a fake cv. It
	  // is NOT included in the nbocv CSR structure...
	  assert(cv_flag[icv_nbr] == -3);
	  cv_flag[icv_nbr] = 3; // switch to positive for checking
	}
	else if (ifa < nfa_bpi) {
	  // this is a ghost face. It may already be included in the
	  // list because two faces can point to the same cv now:
	  // (addresses non-coplanar face problem)
	  if (cv_flag[icv_nbr] == -2) {
	    cv_flag[icv_nbr] = 2; // flip
	    if (iter == 0) {
	      nbocv_i[icv+1] += 1;
	    }
	    else {
	      nbocv_v[nbocv_i[icv]] = icv_nbr;
	      nbocv_i[icv] += 1;
	    }
	  }
	  else {
	    assert(cv_flag[icv_nbr] == 2);
	  }
	}
	else {
	  // this is an internal face. It may already be included in the
	  // list for the same reason as above...
	  if (cv_flag[icv_nbr] == -1) {
	    cv_flag[icv_nbr] = 1; // flip
	    if (iter == 0) {
	      nbocv_i[icv+1] += 1;
	    }
	    else {
	      nbocv_v[nbocv_i[icv]] = icv_nbr;
	      nbocv_i[icv] += 1;
	    }
	  }
	  else {
	    assert(cv_flag[icv_nbr] == 1);
	  }
	}
      }
      // clear flags...
      for (int foc = foc_f; foc <= foc_l; ++foc) {
	const int ifa = faocv_v[foc];
	if (cv_flag[cvofa[ifa][0]] > 0)
	  cv_flag[cvofa[ifa][0]] = -cv_flag[cvofa[ifa][0]];
	if (cv_flag[cvofa[ifa][1]] > 0)
	  cv_flag[cvofa[ifa][1]] = -cv_flag[cvofa[ifa][1]];
      }
    }
    if (iter == 0) {
      nbocv_i[0] = 0;
      FOR_ICV nbocv_i[icv+1] += nbocv_i[icv];
      nbocv_s = nbocv_i[ncv];
      nbocv_v = new int[nbocv_s];
    } 
    else {
      for (int icv = ncv-1; icv > 0; --icv)
        nbocv_i[icv] = nbocv_i[icv-1];
      nbocv_i[0] = 0;
    }
  }
  
}

void UgpWithCv2::calcGradCoeffGG()
{
  if (mpi_rank==0) cout<<" > calcGradCoeffGG()"<<endl;

  // a Green-Gauss gradient reconstruction for cv's
  assert(nbocv_all_grad_coeff==NULL);
  nbocv_all_grad_coeff = new double[nbocv_all_s][3];

  for (int nb = 0; nb < nbocv_all_s; ++nb)
    for (int i = 0; i < 3; ++i)
      nbocv_all_grad_coeff[nb][i] = 0.0;


  for (int ifa = nfa_b; ifa < nfa; ++ifa)
  {
    // these faces have both icv0 and icv1 local...
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];

    int noc00, noc01, noc11, noc10;
    noc00 = nbocv_all_i[icv0];
    noc01 = noc00;
    while (nbocv_all_v[++noc01] != icv1);

    for (int i = 0; i < 3; ++i)
    {
      nbocv_all_grad_coeff[noc00][i] += 0.5*fa_normal[ifa][i]/cv_volume[icv0];
      nbocv_all_grad_coeff[noc01][i] += 0.5*fa_normal[ifa][i]/cv_volume[icv0];
    }

    if (icv1 < ncv)
    {
      noc11 = nbocv_all_i[icv1];
      noc10 = noc11;
      while (nbocv_all_v[++noc10] != icv0);

      for (int i = 0; i < 3; ++i)
      {
        nbocv_all_grad_coeff[noc11][i] -= 0.5*fa_normal[ifa][i]/cv_volume[icv1];
        nbocv_all_grad_coeff[noc10][i] -= 0.5*fa_normal[ifa][i]/cv_volume[icv1];
      }
    }
  }

  for (int ifa = 0; ifa < nfa_b; ++ifa)
  {
    const int icv0 = cvofa[ifa][0];
    const int icv1 = cvofa[ifa][1];
    const int noc00 = nbocv_all_i[icv0];
    int noc01 = noc00;
    while (nbocv_all_v[++noc01] != icv1);

    for (int i = 0; i < 3; ++i)
      nbocv_all_grad_coeff[noc01][i] += fa_normal[ifa][i]/cv_volume[icv0];
  }

  // update the operators into the level-0 ghost data...
  updateCvOperatorG1(nbocv_all_grad_coeff);
}


void UgpWithCv2::calcGradAndFilterCoeffMGG() {
  
  // ======================================================================================
  // MGG: modified green gauss gradients
  //
  // Accurate gradient for cv-based scheme involving only von Neumann (face-based) neighbors
  // and fake cv values (with specified centroids - they do not have to be collocated
  // with the boundary face centroid, but this makes the most sense to me).
  //
  // The same principle can be applied to develop a linear-preserving cv-based filter
  // also involving only von Neumann (face-based) neighbors.
  //
  // F. Ham
  // April 2010
  // ======================================================================================

  if (mpi_rank == 0)
    cout << " > calcGradAndFilterCoeffMGG()" << endl;
  
  // a Green-Gauss gradient reconstruction for cv's
  assert( nbocv_all_grad_coeff == NULL );
  nbocv_all_grad_coeff = new double[nbocv_all_s][3];
  
  assert( nbocv_all_filter_coeff == NULL );
  nbocv_all_filter_coeff = new double[nbocv_all_s];
  
  // the key to these gradient and filter operators is to rebuild nodal quantities
  // using the nearest cv values by linear interpolation/extrapolation (in some cases).
  // The gradient and filter operators are then computed based on manipulating 
  // the nodal values, and finally written in terms of the cv neighbors.
  
  int faono_i_tmp[NNOC_MAX+1];
  const int FAONO_S_MAX = 4*NNOC_MAX; // assume 4 faces per node
  int faono_v_tmp[FAONO_S_MAX];
  double faono_weight_tmp[FAONO_S_MAX]; 
    
  FOR_INO no_flag[ino] = -1;
  FOR_IFA fa_flag[ifa] = -1;
  FOR_ICV {
      
    int noc_f = noocv_i[icv];
    int noc_l = noocv_i[icv+1]-1;
    for (int noc = noc_f; noc <= noc_l; ++noc) {
      int ino = noocv_v[noc];
      assert( no_flag[ino] == -1 );
      no_flag[ino] = noc-noc_f;
    }
    int nnoc = noc_l-noc_f+1;
    assert( nnoc <= NNOC_MAX );

    // at the end of this routine, the grad coeff's for each nbr
    // cv will be stored in faocv_grad. This can accomodate the 
    // fact that fake cvs are not cv neighbors...

    int foc_f = faocv_i[icv];
    int foc_l = faocv_i[icv+1]-1;
    for (int foc = foc_f; foc <= foc_l; ++foc) {
      int ifa = faocv_v[foc];
      assert( fa_flag[ifa] == -1 );
      fa_flag[ifa] = foc-foc_f;
    }
    int nfoc = foc_l-foc_f+1;
    assert( nfoc <= NFOC_MAX );

    // build local faono for this cv...

    for (int ii = 0; ii < nnoc; ++ii)
      faono_i_tmp[ii+1] = 0;
    int faono_s_tmp;
      
    for (int iter = 0; iter <= 1; ++iter) {
      for (int foc = foc_f; foc <= foc_l; ++foc) {
	int ifa = faocv_v[foc];	
	int nof_f = noofa_i[ifa];
	int nof_l = noofa_i[ifa+1]-1;
	for (int nof = nof_f; nof <= nof_l; ++nof) {
	  int ino = noofa_v[nof];
	  int ii = no_flag[ino];
	  assert( ii >= 0 );
	  if (iter == 0) {
	    faono_i_tmp[ii+1] += 1;
	  }
	  else {
	    faono_v_tmp[faono_i_tmp[ii]] = ifa;
	    faono_i_tmp[ii] += 1;
	  }
	}
      }
      if (iter == 0) {
	faono_i_tmp[0] = 0;
	for (int ii = 0; ii < nnoc; ++ii)
	  faono_i_tmp[ii+1] += faono_i_tmp[ii];
	faono_s_tmp = faono_i_tmp[nnoc];
	assert( faono_s_tmp <= FAONO_S_MAX );
      }
      else {
	for (int ii = nnoc; ii > 0; --ii)
	  faono_i_tmp[ii] = faono_i_tmp[ii-1];
	faono_i_tmp[0] = 0;
      }
    }

    // for nodes with just 2 nbrs, additional nbrs associated with edge-shared nodes must be 
    // added...
    int vshift = 0;
    for (int noc = noc_f; noc <= noc_l; ++noc) {
      int ii = noc-noc_f;
      if (faono_i_tmp[ii+1]-faono_i_tmp[ii] == 2) {
	int fon_f = faono_i_tmp[ii];
	int ifa0 = faono_v_tmp[fon_f  ];
	int ifa1 = faono_v_tmp[fon_f+1];
	int ino = noocv_v[noc];
	int ino0 = getPrevNodeOfFaceCCW(ifa0,ino,icv);
	assert( getNextNodeOfFaceCCW(ifa1,ino,icv) == ino0 );
	int ii0 = no_flag[ino0];
	vshift += faono_i_tmp[ii0+1]-faono_i_tmp[ii0]-2;
	int ino1 = getNextNodeOfFaceCCW(ifa0,ino,icv);
	assert( getPrevNodeOfFaceCCW(ifa1,ino,icv) == ino1 );
	int ii1 = no_flag[ino1];
	vshift += faono_i_tmp[ii1+1]-faono_i_tmp[ii1]-2;
      }
    }
    if (vshift > 0) {
      faono_s_tmp += vshift;
      assert( faono_s_tmp <= FAONO_S_MAX );
      // cycle through in reverse...
      for (int noc = noc_l; noc >= noc_f; --noc) {
	int ii = noc-noc_f;
	int fon_f = faono_i_tmp[ii];
	int fon_l = faono_i_tmp[ii+1]-1;
	faono_i_tmp[ii+1] += vshift;
	if (fon_l-fon_f+1 > 2) {
	  // for any node that already has 3 or more nbrs, just shift up...
	  for (int fon = fon_l; fon >= fon_f; --fon)
	    faono_v_tmp[fon+vshift] = faono_v_tmp[fon];
	}
	else {
	  assert(fon_l-fon_f+1 == 2);
	  int ifa0 = faono_v_tmp[fon_f  ];
	  int ifa1 = faono_v_tmp[fon_f+1];
	  int ino = noocv_v[noc];
	  int ino0 = getPrevNodeOfFaceCCW(ifa0,ino,icv);
	  assert( getNextNodeOfFaceCCW(ifa1,ino,icv) == ino0 );
	  int ii0 = no_flag[ino0];
	  int fon0_f = faono_i_tmp[ii0];
	  int fon0_l = faono_i_tmp[ii0+1]-1;
	  for (int fon0 = fon0_l; fon0 >= fon0_f; --fon0) {
	    int ifa_nbr = faono_v_tmp[fon0];
	    if ( (ifa_nbr != ifa0)&&(ifa_nbr != ifa1) ) {
	      faono_v_tmp[fon_l+vshift] = ifa_nbr;
	      --vshift;
	    }
	  }
	  int ino1 = getNextNodeOfFaceCCW(ifa0,ino,icv);
	  assert( getPrevNodeOfFaceCCW(ifa1,ino,icv) == ino1 );
	  int ii1 = no_flag[ino1];
	  int fon1_f = faono_i_tmp[ii1];
	  int fon1_l = faono_i_tmp[ii1+1]-1;
	  for (int fon1 = fon1_l; fon1 >= fon1_f; --fon1) {
	    int ifa_nbr = faono_v_tmp[fon1];
	    if ( (ifa_nbr != ifa0)&&(ifa_nbr != ifa1) ) {
	      faono_v_tmp[fon_l+vshift] = ifa_nbr;
	      --vshift;
	    }
	  }
	  // and add original pair of faces...
	  faono_v_tmp[fon_f+1+vshift] = ifa1;
	  faono_v_tmp[fon_f+vshift]   = ifa0;
	}
      }
      assert( vshift == 0 );
    }
      
    // now sort faces and build the geomtric interpolation factors 
    // that determine each node's value as a function of surrounding cv's...
      
    for (int ii = 0; ii < faono_s_tmp; ++ii)
      faono_weight_tmp[ii] = -1000.0;
      
    for (int noc = noc_f; noc <= noc_l; ++noc) {
      int ino = noocv_v[noc];
      int ii = noc-noc_f;
      int fon_f = faono_i_tmp[ii];
      int fon_l = faono_i_tmp[ii+1]-1;

      switch (fon_l-fon_f+1) {

      case 3:
	{
	  int ifa0 = faono_v_tmp[fon_f  ];
	  int ifa1 = faono_v_tmp[fon_f+1];
	  int ifa2 = faono_v_tmp[fon_f+2];
	  int ino01 = getNextNodeOfFaceCCW(ifa0,ino,icv);
	  if ( getPrevNodeOfFaceCCW(ifa1,ino,icv) != ino01 ) {
	    // flip ifa1 and ifa2...
	    faono_v_tmp[fon_f+1] = ifa2;
	    faono_v_tmp[fon_f+2] = ifa1;
	    ifa1 = faono_v_tmp[fon_f+1];
	    ifa2 = faono_v_tmp[fon_f+2];
	    assert( getPrevNodeOfFaceCCW(ifa1,ino,icv) == ino01 );
	  }
	  // now should all be in order...
	  int ino12 = getNextNodeOfFaceCCW(ifa1,ino,icv);
	  assert( ino12 == getPrevNodeOfFaceCCW(ifa2,ino,icv) );
	  int ino20 = getNextNodeOfFaceCCW(ifa2,ino,icv);
	  assert( ino20 == getPrevNodeOfFaceCCW(ifa0,ino,icv) );
	  // and now the interpolation values...
	  int icv0;
	  if (cvofa[ifa0][0] == icv)
	    icv0 = cvofa[ifa0][1];
	  else {
	    assert( cvofa[ifa0][1] == icv );
	    icv0 = cvofa[ifa0][0];
	  }
	  int icv1;
	  if (cvofa[ifa1][0] == icv)
	    icv1 = cvofa[ifa1][1];
	  else {
	    assert( cvofa[ifa1][1] == icv );
	    icv1 = cvofa[ifa1][0];
	  }
	  int icv2;
	  if (cvofa[ifa2][0] == icv)
	    icv2 = cvofa[ifa2][1];
	  else {
	    assert( cvofa[ifa2][1] == icv );
	    icv2 = cvofa[ifa2][0];
	  }

	  double dx0[3]; FOR_I3 dx0[i] = x_cv[icv0][i] - x_cv[icv][i];
	  double dx1[3]; FOR_I3 dx1[i] = x_cv[icv1][i] - x_cv[icv][i];
	  double dx2[3]; FOR_I3 dx2[i] = x_cv[icv2][i] - x_cv[icv][i];
	  double dx[3]; FOR_I3 dx[i] = x_no[ino][i] - x_cv[icv][i];
	  double w0 = 
	    dx[0]*(dx1[1]*dx2[2]-dx1[2]*dx2[1]) +
	    dx[1]*(dx1[2]*dx2[0]-dx1[0]*dx2[2]) +
	    dx[2]*(dx1[0]*dx2[1]-dx1[1]*dx2[0]);
	  double w1 = 
	    dx[0]*(dx2[1]*dx0[2]-dx2[2]*dx0[1]) +
	    dx[1]*(dx2[2]*dx0[0]-dx2[0]*dx0[2]) +
	    dx[2]*(dx2[0]*dx0[1]-dx2[1]*dx0[0]);
	  double w2 = 
	    dx[0]*(dx0[1]*dx1[2]-dx0[2]*dx1[1]) +
	    dx[1]*(dx0[2]*dx1[0]-dx0[0]*dx1[2]) +
	    dx[2]*(dx0[0]*dx1[1]-dx0[1]*dx1[0]);
	  double w = 
	    dx0[0]*(dx1[1]*dx2[2]-dx1[2]*dx2[1]) +
	    dx0[1]*(dx1[2]*dx2[0]-dx1[0]*dx2[2]) +
	    dx0[2]*(dx1[0]*dx2[1]-dx1[1]*dx2[0]);
	  assert( faono_weight_tmp[fon_f  ] == -1000.0 );
	  assert( faono_weight_tmp[fon_f+1] == -1000.0 );
	  assert( faono_weight_tmp[fon_f+2] == -1000.0 );
	  faono_weight_tmp[fon_f  ] = w0/w;
	  faono_weight_tmp[fon_f+1] = w1/w;
	  faono_weight_tmp[fon_f+2] = w2/w;
	}
	break;

      case 4:
	{
	  // here we are going to use a inverse-distance weighted least squares,
	  // so ordering doesn't matter...
	  int ifa0 = faono_v_tmp[fon_f  ];
	  int ifa1 = faono_v_tmp[fon_f+1];
	  int ifa2 = faono_v_tmp[fon_f+2];
	  int ifa3 = faono_v_tmp[fon_f+3];

	  // and the interpolation values...
	  int icv0;
	  if (cvofa[ifa0][0] == icv)
	    icv0 = cvofa[ifa0][1];
	  else {
	    assert( cvofa[ifa0][1] == icv );
	    icv0 = cvofa[ifa0][0];
	  }
	  int icv1;
	  if (cvofa[ifa1][0] == icv)
	    icv1 = cvofa[ifa1][1];
	  else {
	    assert( cvofa[ifa1][1] == icv );
	    icv1 = cvofa[ifa1][0];
	  }
	  int icv2;
	  if (cvofa[ifa2][0] == icv)
	    icv2 = cvofa[ifa2][1];
	  else {
	    assert( cvofa[ifa2][1] == icv );
	    icv2 = cvofa[ifa2][0];
	  }
	  int icv3;
	  if (cvofa[ifa3][0] == icv)
	    icv3 = cvofa[ifa3][1];
	  else {
	    assert( cvofa[ifa3][1] == icv );
	    icv3 = cvofa[ifa3][0];
	  }
	  // we need to get weights for the nodal value in terms of the
	  // surrounding cv values. Weight these by inverse distance from the node...
	  double w0 = 0.0;
	  FOR_I3 {
	    double dx = x_cv[icv0][i]-x_no[ino][i];
	    w0 += dx*dx;
	  }
	  w0 = 1.0/sqrt(w0);
	  double w1 = 0.0;
	  FOR_I3 {
	    double dx = x_cv[icv1][i]-x_no[ino][i];
	    w1 += dx*dx;
	  }
	  w1 = 1.0/sqrt(w1);
	  double w2 = 0.0;
	  FOR_I3 {
	    double dx = x_cv[icv2][i]-x_no[ino][i];
	    w2 += dx*dx;
	  }
	  w2 = 1.0/sqrt(w2);
	  double w3 = 0.0;
	  FOR_I3 {
	    double dx = x_cv[icv3][i]-x_no[ino][i];
	    w3 += dx*dx;
	  }
	  w3 = 1.0/sqrt(w3);
	  // dx's are relative to the icv, where the polynomial is constrained...
	  double dx0[3]; FOR_I3 dx0[i] = x_cv[icv0][i] - x_cv[icv][i];
	  double dx1[3]; FOR_I3 dx1[i] = x_cv[icv1][i] - x_cv[icv][i];
	  double dx2[3]; FOR_I3 dx2[i] = x_cv[icv2][i] - x_cv[icv][i];
	  double dx3[3]; FOR_I3 dx3[i] = x_cv[icv3][i] - x_cv[icv][i];
	  // sums...
	  double swxx = w0*dx0[0]*dx0[0] + w1*dx1[0]*dx1[0] + w2*dx2[0]*dx2[0] + w3*dx3[0]*dx3[0];
	  double swxy = w0*dx0[0]*dx0[1] + w1*dx1[0]*dx1[1] + w2*dx2[0]*dx2[1] + w3*dx3[0]*dx3[1];
	  double swzx = w0*dx0[0]*dx0[2] + w1*dx1[0]*dx1[2] + w2*dx2[0]*dx2[2] + w3*dx3[0]*dx3[2];
	  double swyy = w0*dx0[1]*dx0[1] + w1*dx1[1]*dx1[1] + w2*dx2[1]*dx2[1] + w3*dx3[1]*dx3[1];
	  double swyz = w0*dx0[1]*dx0[2] + w1*dx1[1]*dx1[2] + w2*dx2[1]*dx2[2] + w3*dx3[1]*dx3[2];
	  double swzz = w0*dx0[2]*dx0[2] + w1*dx1[2]*dx1[2] + w2*dx2[2]*dx2[2] + w3*dx3[2]*dx3[2];
	  // denom...
	  double denom = 
	    swxy*swxy*swzz + 
	    swzx*swzx*swyy + 
	    swyz*swyz*swxx - 
	    swxx*swyy*swzz - 
	    2.0*swxy*swyz*swzx;
	  double dx[3]; FOR_I3 dx[i] = x_no[ino][i] - x_cv[icv][i];
	  double cx = ( dx[0]*(swyz*swyz - swyy*swzz) + dx[1]*(swxy*swzz - swzx*swyz) + dx[2]*(swzx*swyy - swxy*swyz) )/denom;
	  double cy = ( dx[1]*(swzx*swzx - swzz*swxx) + dx[2]*(swyz*swxx - swxy*swzx) + dx[0]*(swxy*swzz - swyz*swzx) )/denom;
	  double cz = ( dx[2]*(swxy*swxy - swxx*swyy) + dx[0]*(swzx*swyy - swyz*swxy) + dx[1]*(swyz*swxx - swzx*swxy) )/denom;
	  w0 *= cx*dx0[0] + cy*dx0[1] + cz*dx0[2];
	  w1 *= cx*dx1[0] + cy*dx1[1] + cz*dx1[2];
	  w2 *= cx*dx2[0] + cy*dx2[1] + cz*dx2[2];
	  w3 *= cx*dx3[0] + cy*dx3[1] + cz*dx3[2];
	  assert( faono_weight_tmp[fon_f  ] == -1000.0 );
	  assert( faono_weight_tmp[fon_f+1] == -1000.0 );
	  assert( faono_weight_tmp[fon_f+2] == -1000.0 );
	  assert( faono_weight_tmp[fon_f+3] == -1000.0 );
	  faono_weight_tmp[fon_f  ] = w0;
	  faono_weight_tmp[fon_f+1] = w1;
	  faono_weight_tmp[fon_f+2] = w2;
	  faono_weight_tmp[fon_f+3] = w3;
	}
	break;

      case 5:
	{
	  // here we are going to use a inverse-distance weighted least squares,
	  // so ordering doesn't matter. XXXX In the future, we should use a version
	  // of this reconstruction that constrains the reconstructions associated
	  // with shared faces, but for now, just use inverse distance and hope for the
	  // best...
	  int ifa0 = faono_v_tmp[fon_f  ];
	  int ifa1 = faono_v_tmp[fon_f+1];
	  int ifa2 = faono_v_tmp[fon_f+2];
	  int ifa3 = faono_v_tmp[fon_f+3];
	  int ifa4 = faono_v_tmp[fon_f+4];

	  // and the interpolation values...
	  int icv0;
	  if (cvofa[ifa0][0] == icv)
	    icv0 = cvofa[ifa0][1];
	  else {
	    assert( cvofa[ifa0][1] == icv );
	    icv0 = cvofa[ifa0][0];
	  }
	  int icv1;
	  if (cvofa[ifa1][0] == icv)
	    icv1 = cvofa[ifa1][1];
	  else {
	    assert( cvofa[ifa1][1] == icv );
	    icv1 = cvofa[ifa1][0];
	  }
	  int icv2;
	  if (cvofa[ifa2][0] == icv)
	    icv2 = cvofa[ifa2][1];
	  else {
	    assert( cvofa[ifa2][1] == icv );
	    icv2 = cvofa[ifa2][0];
	  }
	  int icv3;
	  if (cvofa[ifa3][0] == icv)
	    icv3 = cvofa[ifa3][1];
	  else {
	    assert( cvofa[ifa3][1] == icv );
	    icv3 = cvofa[ifa3][0];
	  }
	  int icv4;
	  if (cvofa[ifa4][0] == icv)
	    icv4 = cvofa[ifa4][1];
	  else {
	    assert( cvofa[ifa4][1] == icv );
	    icv4 = cvofa[ifa4][0];
	  }
	  // we need to get weights for the nodal value in terms of the
	  // surrounding cv values. Weight these by inverse distance from the node...
	  double w0 = 0.0;
	  FOR_I3 {
	    double dx = x_cv[icv0][i]-x_no[ino][i];
	    w0 += dx*dx;
	  }
	  w0 = 1.0/sqrt(w0);
	  double w1 = 0.0;
	  FOR_I3 {
	    double dx = x_cv[icv1][i]-x_no[ino][i];
	    w1 += dx*dx;
	  }
	  w1 = 1.0/sqrt(w1);
	  double w2 = 0.0;
	  FOR_I3 {
	    double dx = x_cv[icv2][i]-x_no[ino][i];
	    w2 += dx*dx;
	  }
	  w2 = 1.0/sqrt(w2);
	  double w3 = 0.0;
	  FOR_I3 {
	    double dx = x_cv[icv3][i]-x_no[ino][i];
	    w3 += dx*dx;
	  }
	  w3 = 1.0/sqrt(w3);
	  double w4 = 0.0;
	  FOR_I3 {
	    double dx = x_cv[icv4][i]-x_no[ino][i];
	    w4 += dx*dx;
	  }
	  w4 = 1.0/sqrt(w4);
	  // dx's are relative to the icv, where the polynomial is constrained...
	  double dx0[3]; FOR_I3 dx0[i] = x_cv[icv0][i] - x_cv[icv][i];
	  double dx1[3]; FOR_I3 dx1[i] = x_cv[icv1][i] - x_cv[icv][i];
	  double dx2[3]; FOR_I3 dx2[i] = x_cv[icv2][i] - x_cv[icv][i];
	  double dx3[3]; FOR_I3 dx3[i] = x_cv[icv3][i] - x_cv[icv][i];
	  double dx4[3]; FOR_I3 dx4[i] = x_cv[icv4][i] - x_cv[icv][i];
	  // sums...
	  double swxx = w0*dx0[0]*dx0[0] + w1*dx1[0]*dx1[0] + w2*dx2[0]*dx2[0] + w3*dx3[0]*dx3[0] + w4*dx4[0]*dx4[0];
	  double swxy = w0*dx0[0]*dx0[1] + w1*dx1[0]*dx1[1] + w2*dx2[0]*dx2[1] + w3*dx3[0]*dx3[1] + w4*dx4[0]*dx4[1];
	  double swzx = w0*dx0[0]*dx0[2] + w1*dx1[0]*dx1[2] + w2*dx2[0]*dx2[2] + w3*dx3[0]*dx3[2] + w4*dx4[0]*dx4[2];
	  double swyy = w0*dx0[1]*dx0[1] + w1*dx1[1]*dx1[1] + w2*dx2[1]*dx2[1] + w3*dx3[1]*dx3[1] + w4*dx4[1]*dx4[1];
	  double swyz = w0*dx0[1]*dx0[2] + w1*dx1[1]*dx1[2] + w2*dx2[1]*dx2[2] + w3*dx3[1]*dx3[2] + w4*dx4[1]*dx4[2];
	  double swzz = w0*dx0[2]*dx0[2] + w1*dx1[2]*dx1[2] + w2*dx2[2]*dx2[2] + w3*dx3[2]*dx3[2] + w4*dx4[2]*dx4[2];
	  // denom...
	  double denom = 
	    swxy*swxy*swzz + 
	    swzx*swzx*swyy + 
	    swyz*swyz*swxx - 
	    swxx*swyy*swzz - 
	    2.0*swxy*swyz*swzx;
	  double dx[3]; FOR_I3 dx[i] = x_no[ino][i] - x_cv[icv][i];
	  double cx = ( dx[0]*(swyz*swyz - swyy*swzz) + dx[1]*(swxy*swzz - swzx*swyz) + dx[2]*(swzx*swyy - swxy*swyz) )/denom;
	  double cy = ( dx[1]*(swzx*swzx - swzz*swxx) + dx[2]*(swyz*swxx - swxy*swzx) + dx[0]*(swxy*swzz - swyz*swzx) )/denom;
	  double cz = ( dx[2]*(swxy*swxy - swxx*swyy) + dx[0]*(swzx*swyy - swyz*swxy) + dx[1]*(swyz*swxx - swzx*swxy) )/denom;
	  w0 *= cx*dx0[0] + cy*dx0[1] + cz*dx0[2];
	  w1 *= cx*dx1[0] + cy*dx1[1] + cz*dx1[2];
	  w2 *= cx*dx2[0] + cy*dx2[1] + cz*dx2[2];
	  w3 *= cx*dx3[0] + cy*dx3[1] + cz*dx3[2];
	  w4 *= cx*dx4[0] + cy*dx4[1] + cz*dx4[2];
	  assert( faono_weight_tmp[fon_f  ] == -1000.0 );
	  assert( faono_weight_tmp[fon_f+1] == -1000.0 );
	  assert( faono_weight_tmp[fon_f+2] == -1000.0 );
	  assert( faono_weight_tmp[fon_f+3] == -1000.0 );
	  assert( faono_weight_tmp[fon_f+4] == -1000.0 );
	  faono_weight_tmp[fon_f  ] = w0;
	  faono_weight_tmp[fon_f+1] = w1;
	  faono_weight_tmp[fon_f+2] = w2;
	  faono_weight_tmp[fon_f+3] = w3;
	  faono_weight_tmp[fon_f+4] = w4;
	}
	break;

      default:
	cerr << "Error: strange number of faces per node: " << fon_l-fon_f+1 << endl;
	throw(-1);
      }
    }

    // ====================================================================
    // we now have the interpolation weights associated with each node...
    // now loop through the faces and build the gradient...
    // ====================================================================
      
    double faocv_grad[NFOC_MAX][3];
    for (int ii = 0; ii < nfoc; ++ii)
      FOR_I3 faocv_grad[ii][i] = 0.0;
    double cv_vol = 0.0;
          
    double faocv_filter[NFOC_MAX];
    for (int ii = 0; ii < nfoc; ++ii)
      faocv_filter[ii] = 0.0;

    // the filter requires noocv_weight, the weight of the cv centroid wrt the
    // cv's nodes. Since the cv centroid is already calculated, we can compute this
    // using the following formula...
    //
    // 
    // (sum_{sub-tet}(mag))*x_cv = sum_{sub-tet}(mag*(x_cv+x_fa+x_no0+x_no1)/4);
    // 
    // and solve for x_cv.
    
    double noocv_weight[NNOC_MAX];
    for (int ii = 0; ii < nnoc; ++ii)
      noocv_weight[ii] = 0.0;
    double noocv_weight_sum = 0.0;

    for (int foc = foc_f; foc <= foc_l; ++foc) {
	
      int ifa = faocv_v[foc];
      double fa_area = sqrt( fa_normal[ifa][0]*fa_normal[ifa][0] + 
			     fa_normal[ifa][1]*fa_normal[ifa][1] + 
			     fa_normal[ifa][2]*fa_normal[ifa][2] );
      double unit_normal[3];
      FOR_I3 unit_normal[i] = fa_normal[ifa][i]/fa_area;
	
      // loop thru face nodes in CCW order...
      int nof_f = noofa_i[ifa];
      int nof_begin,nof_end,nof_inc,nnof;
      if (cvofa[ifa][0] == icv) {
	nof_begin = noofa_i[ifa];
	nof_end = noofa_i[ifa+1];
	nof_inc = 1;
	nnof = nof_end - nof_begin;
      }
      else {
	assert(cvofa[ifa][1] == icv);
	nof_begin = noofa_i[ifa+1]-1;
	nof_end = noofa_i[ifa]-1;
	nof_inc = -1;
	nnof = nof_begin - nof_end;
	FOR_I3 unit_normal[i] = -unit_normal[i];
      }
      
      // compute the set of of weights required to interpolate to the
      // face centroid. Since the centroid is already calculated, just use the
      // following formula....
      //
      // (sum_{sub-tri}(mag))*x_fa = sum_{sub-tri}(mag*(x_fa+x_no0+x_no1)/3);
      // 
      // and solve for x_fa.

      double noofa_weight[NNOF_MAX];
      for (int ii = 0; ii < nnof; ++ii)
	noofa_weight[ii] = 0.0;
      double noofa_weight_sum = 0.0;
	
      int nof1 = nof_end-nof_inc; 
      int ino1 = noofa_v[nof1];
      for (int nof = nof_begin; nof != nof_end; nof += nof_inc) {
	int nof0 = nof1;
	int ino0 = ino1;
	nof1 = nof;
	ino1 = noofa_v[nof1];
	double dx0[3]; FOR_I3 dx0[i] = x_no[ino0][i] - x_fa[ifa][i];
	double dx1[3]; FOR_I3 dx1[i] = x_no[ino1][i] - x_fa[ifa][i];
	double this_normal[3];
	this_normal[0] = dx0[1]*dx1[2] - dx0[2]*dx1[1];
	this_normal[1] = dx0[2]*dx1[0] - dx0[0]*dx1[2];
	this_normal[2] = dx0[0]*dx1[1] - dx0[1]*dx1[0];
	// note this is the fluent definition (I think) that always produces
	// positive weights, although the convergence is potentially very slow
	/*
	  double mag = sqrt( this_normal[0]*this_normal[0] + 
	  this_normal[1]*this_normal[1] + 
	  this_normal[2]*this_normal[2] );
	*/
	// dotting on the uniquely defined unit_normal (does not depend on the 
	// face center/centroid) results in faster convergence, but some of these 
	// weights could be negative if the face-handedness is right on the edge. The
	// only important thing is that this reconstruction is the same as that used in 
	// the calcGeometry routine...
	double mag = 
	  this_normal[0]*unit_normal[0] +
	  this_normal[1]*unit_normal[1] +
	  this_normal[2]*unit_normal[2];
	// and add sub-area to sub-tri corners...
	noofa_weight[nof0-nof_f] += mag;
	noofa_weight[nof1-nof_f] += mag;
	noofa_weight_sum += mag*2.0;
      }
      // normalize and check...
      double x_fa_check[3] = { 0.0, 0.0, 0.0 };
      for (int nof = nof_begin; nof != nof_end; nof += nof_inc) {
	noofa_weight[nof-nof_f] /= noofa_weight_sum;
	int ino = noofa_v[nof];
	FOR_I3 x_fa_check[i] += noofa_weight[nof-nof_f]*x_no[ino][i];
      }
      double d2 = 0.0;
      FOR_I3 {
	double dx = x_fa_check[i] - x_fa[ifa][i];
	d2 += dx*dx;
      }
      //cout << "dist: " << sqrt(d2) << endl; 
      assert( sqrt(d2) < 1.0E-10*sqrt(noofa_weight_sum) ); // noofa_weight_sum is proportional to the face area
      
      // -----------------------------
      // now loop on sub-tets...
      // -----------------------------
      ino1 = noofa_v[nof_end-nof_inc];
      int noc1 = no_flag[ino1]; // already has noc_f subtracted
      for (int nof = nof_begin; nof != nof_end; nof += nof_inc) {
	int ino0 = ino1;
	int noc0 = noc1;
	// this is the subtet from ino0-ino1-ifa-icv...
	ino1 = noofa_v[nof];
	noc1 = no_flag[ino1]; // already has noc_f subtracted
	double dx0[3]; FOR_I3 dx0[i] = x_no[ino0][i] - x_fa[ifa][i];
	double dx1[3]; FOR_I3 dx1[i] = x_no[ino1][i] - x_fa[ifa][i];
	double this_normal[3];
	this_normal[0] = dx0[1]*dx1[2] - dx0[2]*dx1[1];
	this_normal[1] = dx0[2]*dx1[0] - dx0[0]*dx1[2];
	this_normal[2] = dx0[0]*dx1[1] - dx0[1]*dx1[0];
	double this_vol = 
	  (x_fa[ifa][0]-x_cv[icv][0])*this_normal[0] + 
	  (x_fa[ifa][1]-x_cv[icv][1])*this_normal[1] + 
	  (x_fa[ifa][2]-x_cv[icv][2])*this_normal[2];
	assert( this_vol > 0.0 );
	// and add sub-volume to sub-tet corners...
	noocv_weight[noc0] += this_vol;
	noocv_weight[noc1] += this_vol;
	for (int nof2 = nof_begin; nof2 != nof_end; nof2 += nof_inc) {
	  int ino = noofa_v[nof2];
	  int noc = no_flag[ino]; // already has noc_f subtracted
	  noocv_weight[noc] += this_vol*noofa_weight[nof2-nof_f];
	}
	noocv_weight_sum += this_vol*3.0;
      }

    }
    // normalize and check...
    double x_cv_check[3] = { 0.0, 0.0, 0.0 };
    for (int noc = noc_f; noc <= noc_l; ++noc) {
      noocv_weight[noc-noc_f] /= noocv_weight_sum;
      int ino = noocv_v[noc];
      FOR_I3 x_cv_check[i] += noocv_weight[noc-noc_f]*x_no[ino][i];
    }
    double d2 = 0.0;
    FOR_I3 {
      double dx = x_cv_check[i] - x_cv[icv][i];
      d2 += dx*dx;
    }
    //cout << "dist: " << sqrt(d2) << endl; 
    assert( sqrt(d2) < 1.0E-10*pow(noocv_weight_sum,1.0/3.0) ); // noocv_weight_sum is proportional to the cv volume

    // produces a ratio wrt cv_volume of exactly 18...
    //cout << "done icv: " << icv << " noocv_weight_sum/cv_volume: " << noocv_weight_sum/cv_volume[icv] << endl;
    //getchar();
    
    // ======================================================================
    // now build operators...
    // ======================================================================

    for (int foc = foc_f; foc <= foc_l; ++foc) {
	
      int ifa = faocv_v[foc];
      double fa_area = sqrt( fa_normal[ifa][0]*fa_normal[ifa][0] + 
			     fa_normal[ifa][1]*fa_normal[ifa][1] + 
			     fa_normal[ifa][2]*fa_normal[ifa][2] );
      double unit_normal[3];
      FOR_I3 unit_normal[i] = fa_normal[ifa][i]/fa_area;
	
      // loop thru face nodes in CCW order...
      int nof_f = noofa_i[ifa];
      int nof_begin,nof_end,nof_inc,nnof;
      if (cvofa[ifa][0] == icv) {
	nof_begin = noofa_i[ifa];
	nof_end = noofa_i[ifa+1];
	nof_inc = 1;
	nnof = nof_end - nof_begin;
      }
      else {
	assert(cvofa[ifa][1] == icv);
	nof_begin = noofa_i[ifa+1]-1;
	nof_end = noofa_i[ifa]-1;
	nof_inc = -1;
	nnof = nof_begin - nof_end;
	FOR_I3 unit_normal[i] = -unit_normal[i];
      }
      
      // re-compute the set of of weights required to interpolate to the
      // face centroid. Since the centroid is already calculated, just use the
      // following formula....
      //
      // (sum_{sub-tri}(mag))*x_fa = sum_{sub-tri}(mag*(x_fa+x_no0+x_no1)/3);
      // 
      // and solve for x_fa.

      double noofa_weight[NNOF_MAX];
      for (int ii = 0; ii < nnof; ++ii)
	noofa_weight[ii] = 0.0;
      double noofa_weight_sum = 0.0;
	
      int nof1 = nof_end-nof_inc; 
      int ino1 = noofa_v[nof1];
      for (int nof = nof_begin; nof != nof_end; nof += nof_inc) {
	int nof0 = nof1;
	int ino0 = ino1;
	nof1 = nof;
	ino1 = noofa_v[nof1];
	double dx0[3]; FOR_I3 dx0[i] = x_no[ino0][i] - x_fa[ifa][i];
	double dx1[3]; FOR_I3 dx1[i] = x_no[ino1][i] - x_fa[ifa][i];
	double this_normal[3];
	this_normal[0] = dx0[1]*dx1[2] - dx0[2]*dx1[1];
	this_normal[1] = dx0[2]*dx1[0] - dx0[0]*dx1[2];
	this_normal[2] = dx0[0]*dx1[1] - dx0[1]*dx1[0];
	// note this is the fluent definition (I think) that always produces
	// positive weights, although the convergence is potentially very slow
	/*
	  double mag = sqrt( this_normal[0]*this_normal[0] + 
	  this_normal[1]*this_normal[1] + 
	  this_normal[2]*this_normal[2] );
	*/
	// dotting on the uniquely defined unit_normal (does not depend on the 
	// face center/centroid) results in faster convergence, but some of these 
	// weights could be negative if the face-handedness is right on the edge. The
	// only important thing is that this reconstruction is the same as that used in 
	// the calcGeometry routine...
	double mag = 
	  this_normal[0]*unit_normal[0] +
	  this_normal[1]*unit_normal[1] +
	  this_normal[2]*unit_normal[2];
	noofa_weight[nof0-nof_f] += mag;
	noofa_weight[nof1-nof_f] += mag;
	noofa_weight_sum += mag*2.0;
      }
      // normalize and check...
      double x_fa_check[3] = { 0.0, 0.0, 0.0 };
      for (int nof = nof_begin; nof != nof_end; nof += nof_inc) {
	noofa_weight[nof-nof_f] /= noofa_weight_sum;
	int ino = noofa_v[nof];
	FOR_I3 x_fa_check[i] += noofa_weight[nof-nof_f]*x_no[ino][i];
      }
      double d2 = 0.0;
      FOR_I3 {
	double dx = x_fa_check[i] - x_fa[ifa][i];
	d2 += dx*dx;
      }
      //cout << "dist: " << sqrt(d2) << endl; 
      assert( sqrt(d2) < 1.0E-10*sqrt(noofa_weight_sum) ); // noofa_weight_sum is proportional to the face area
	
      // now assemble gradient...
      nof1 = nof_end-nof_inc; 
      ino1 = noofa_v[nof1];
      for (int nof = nof_begin; nof != nof_end; nof += nof_inc) {
	int nof0 = nof1;
	int ino0 = ino1;
	nof1 = nof;
	ino1 = noofa_v[nof1];
	double dx0[3]; FOR_I3 dx0[i] = x_no[ino0][i] - x_fa[ifa][i];
	double dx1[3]; FOR_I3 dx1[i] = x_no[ino1][i] - x_fa[ifa][i];
	double this_normal[3];
	this_normal[0] = dx0[1]*dx1[2] - dx0[2]*dx1[1];
	this_normal[1] = dx0[2]*dx1[0] - dx0[0]*dx1[2];
	this_normal[2] = dx0[0]*dx1[1] - dx0[1]*dx1[0];
	// this is twice the outward pointing normal between x_fa, x_no0 and x_no1...
	// ino0...
	int fon_f = faono_i_tmp[no_flag[ino0]];
	int fon_l = faono_i_tmp[no_flag[ino0]+1]-1;
	for (int fon = fon_f; fon <= fon_l; ++fon) {
	  int ifa_nbr = faono_v_tmp[fon];
	  double this_weight = faono_weight_tmp[fon];
	  FOR_I3 faocv_grad[fa_flag[ifa_nbr]][i] += this_weight*this_normal[i];
	}
	// ino1...
	fon_f = faono_i_tmp[no_flag[ino1]];
	fon_l = faono_i_tmp[no_flag[ino1]+1]-1;
	for (int fon = fon_f; fon <= fon_l; ++fon) {
	  int ifa_nbr = faono_v_tmp[fon];
	  double this_weight = faono_weight_tmp[fon];
	  FOR_I3 faocv_grad[fa_flag[ifa_nbr]][i] += this_weight*this_normal[i];
	}
	// face...
	for (int nof2 = nof_begin; nof2 != nof_end; nof2 += nof_inc) {
	  int ino2 = noofa_v[nof2];
	  fon_f = faono_i_tmp[no_flag[ino2]];
	  fon_l = faono_i_tmp[no_flag[ino2]+1]-1;
	  for (int fon = fon_f; fon <= fon_l; ++fon) {
	    int ifa_nbr = faono_v_tmp[fon];
	    double this_weight = faono_weight_tmp[fon]*noofa_weight[nof2-nof_f];
	    FOR_I3 faocv_grad[fa_flag[ifa_nbr]][i] += this_weight*this_normal[i];
	  }
	}
	
	// to get the sub-tet volume for normalization, dot with the vector from the 
	// cell centroid to the face centroid
	double this_vol = 
	  this_normal[0]*(x_fa[ifa][0]-x_cv[icv][0]) +
	  this_normal[1]*(x_fa[ifa][1]-x_cv[icv][1]) +
	  this_normal[2]*(x_fa[ifa][2]-x_cv[icv][2]);
	assert( this_vol > 0.0 );
	cv_vol += this_vol;
	
	// for the filtering operator, add 0.25*(ino0+ino1+fa+cv)*this_vol...
	// ino0. Note factor of 0.25 is handled during normalization below...
	fon_f = faono_i_tmp[no_flag[ino0]];
	fon_l = faono_i_tmp[no_flag[ino0]+1]-1;
	for (int fon = fon_f; fon <= fon_l; ++fon) {
	  int ifa_nbr = faono_v_tmp[fon];
	  double this_weight = faono_weight_tmp[fon];
	  faocv_filter[fa_flag[ifa_nbr]] += this_weight*this_vol;
	}
	// ino1...
	fon_f = faono_i_tmp[no_flag[ino1]];
	fon_l = faono_i_tmp[no_flag[ino1]+1]-1;
	for (int fon = fon_f; fon <= fon_l; ++fon) {
	  int ifa_nbr = faono_v_tmp[fon];
	  double this_weight = faono_weight_tmp[fon];
	  faocv_filter[fa_flag[ifa_nbr]] += this_weight*this_vol;
	}
	// face...
	for (int nof2 = nof_begin; nof2 != nof_end; nof2 += nof_inc) {
	  int ino2 = noofa_v[nof2];
	  fon_f = faono_i_tmp[no_flag[ino2]];
	  fon_l = faono_i_tmp[no_flag[ino2]+1]-1;
	  for (int fon = fon_f; fon <= fon_l; ++fon) {
	    int ifa_nbr = faono_v_tmp[fon];
	    double this_weight = faono_weight_tmp[fon]*noofa_weight[nof2-nof_f];
	    faocv_filter[fa_flag[ifa_nbr]] += this_weight*this_vol;
	  }
	}
	// cv...
	// Could comment this out to try skipping cv. On a Cartesian uniform 3D grid, this still
	// leads to a negative diagonal, so instead modify the operators for diagonal 
	// dominance below.
	for (int noc2 = noc_f; noc2 <= noc_l; ++noc2) {
	  int ino2 = noocv_v[noc2];
	  fon_f = faono_i_tmp[no_flag[ino2]];
	  fon_l = faono_i_tmp[no_flag[ino2]+1]-1;
	  for (int fon = fon_f; fon <= fon_l; ++fon) {
	    int ifa_nbr = faono_v_tmp[fon];
	    double this_weight = faono_weight_tmp[fon]*noocv_weight[noc2-noc_f];
	    faocv_filter[fa_flag[ifa_nbr]] += this_weight*this_vol;
	  }
	}
      }
    }
    
    // check on cv_vol...
    assert( fabs(cv_vol-6.0*cv_volume[icv]) < 1.0E-10*cv_volume[icv] );

    // for the filter coefficient, we need to make it weakly diagonally dominant - i.e.
    // we need diag = sum( fabs(coeff) ). To do this we take some fraction, eps, of the
    // existing coeff, and add in (1-eps) of [1, 0, 0, 0...] where 1 is on the diagonal.
    
    // compute eps, the multiplier applied to the filter coeff nbrs...
    double sum_abs_filter_coeff = 0.0;
    double sum_filter_coeff = 0.0;
    int nboc_f = nbocv_all_i[icv]; // diagonal...
    for (int foc = foc_f; foc <= foc_l; ++foc) {
      double filter_coeff = 0.25*faocv_filter[foc-foc_f]/cv_vol;
      sum_filter_coeff += filter_coeff;
      sum_abs_filter_coeff += fabs(filter_coeff);
    }
    double eps = 1.0/(sum_filter_coeff+sum_abs_filter_coeff);

    // now build...
    double diag[3] = { 0.0, 0.0, 0.0 };
    double filter_diag = 1.0;
    for (int foc = foc_f; foc <= foc_l; ++foc) {
      int ifa = faocv_v[foc];
      FOR_I3 diag[i] -= faocv_grad[foc-foc_f][i]/cv_vol;
      filter_diag -= eps*0.25*faocv_filter[foc-foc_f]/cv_vol;
      // the cv-nbr in the all structure should have the same 
      // offset (plus 1 because of the diag) as the faocv_i/v structure...
      int nboc = nboc_f+foc-foc_f+1;
      int icv_nbr = nbocv_all_v[nboc];
      assert( ((cvofa[ifa][0] == icv)&&(cvofa[ifa][1] == icv_nbr)) ||
	      ((cvofa[ifa][1] == icv)&&(cvofa[ifa][0] == icv_nbr)) );
      FOR_I3 nbocv_all_grad_coeff[nboc][i] = faocv_grad[foc-foc_f][i]/cv_vol;
      nbocv_all_filter_coeff[nboc] = eps*0.25*faocv_filter[foc-foc_f]/cv_vol;
    }
    FOR_I3 nbocv_all_grad_coeff[nboc_f][i] = diag[i];
    nbocv_all_filter_coeff[nboc_f] = filter_diag;
    
    // cleanup...
    
    for (int noc = noc_f; noc <= noc_l; ++noc) {
      int ino = noocv_v[noc];
      no_flag[ino] = -1;
    }
    for (int foc = foc_f; foc <= foc_l; ++foc) {
      int ifa = faocv_v[foc];
      fa_flag[ifa] = -1;
    }
    
  }

  // finally update the operators into the level-0 ghost data...
  updateCvOperatorG1(nbocv_all_grad_coeff);
  updateCvOperatorG1(nbocv_all_filter_coeff);
  
}

void UgpWithCv2::updateCvOperatorG1(double * A) {
  
  // pack the operator...
  
  for (list<Prcomm>::iterator prcomm = cvPrcommListG1.begin();
       prcomm != cvPrcommListG1.end(); ++prcomm) {
    
    // count how many doubles we need to pack and send...
    prcomm->npack_v = 0;
    int npack = prcomm->packIndexVec.size();
    for (int i = 0; i < npack; ++i) {
      int icv = prcomm->packIndexVec[i];
      // should be locally owned...
      assert( (icv >= 0)&&(icv < ncv) );
      prcomm->npack_v += nbocv_all_i[icv+1]-nbocv_all_i[icv];
    }

    // ensure size of packBufferDouble...
    prcomm->ensurePackBufferDoubleSize(prcomm->npack_v);

    int ii = 0;
    for (int i = 0; i < npack; ++i) {
      int icv = prcomm->packIndexVec[i];
      int noc_f = nbocv_all_i[icv];
      int noc_l = nbocv_all_i[icv+1]-1;
      for (int noc = noc_f; noc <= noc_l; ++noc) {
	prcomm->packBufferDouble[ii] = A[noc];
	++ii;
      }
    }
    assert( ii == prcomm->npack_v );
    
  }

  // exchange...

  exchangePrcommBufferDoubleV(cvPrcommListG1);

  // unpack...

  for (list<Prcomm>::iterator prcomm = cvPrcommListG1.begin();
       prcomm != cvPrcommListG1.end(); ++prcomm) {
    
    int ii = 0;
    int nunpack = prcomm->unpackIndexVec.size();
    for (int i = 0; i < nunpack; ++i) {
      int icv = prcomm->unpackIndexVec[i];
      // should be level-0 ghost...
      assert( (icv >= ncv)&&(icv < ncv_g) );
      int noc_f = nbocv_all_i[icv];
      int noc_l = nbocv_all_i[icv+1]-1;
      for (int noc = noc_f; noc <= noc_l; ++noc) {
	A[noc] = prcomm->unpackBufferDouble[ii];
	++ii;
      }
    }
    assert( ii == prcomm->nunpack_v );
    
  }
  
}

void UgpWithCv2::updateCvOperatorG1(double (*A)[3]) {
  
  // pack the operator...
  
  for (list<Prcomm>::iterator prcomm = cvPrcommListG1.begin();
       prcomm != cvPrcommListG1.end(); ++prcomm) {
    
    // count how many doubles we need to pack and send...
    prcomm->npack_v = 0;
    int npack = prcomm->packIndexVec.size();
    for (int i = 0; i < npack; ++i) {
      int icv = prcomm->packIndexVec[i];
      // should be locally owned...
      assert( (icv >= 0)&&(icv < ncv) );
      prcomm->npack_v += nbocv_all_i[icv+1]-nbocv_all_i[icv];
    }

    // ensure size of packBufferDouble (accounting for vector nature of the operator)...
    prcomm->ensurePackBufferDoubleSize(3*prcomm->npack_v);
    
    // loop on Ranges to allow operator transformation for
    // periodic FaZones...
    int ii = 0;
    for (list<Range>::iterator range = prcomm->packRangeList.begin();
	 range != prcomm->packRangeList.end(); ++range) {
      int ii_f = ii;
      for (int i = range->getIndexFirst(); i <= range->getIndexLast(); ++i) {
	int icv = prcomm->packIndexVec[i];
	int noc_f = nbocv_all_i[icv];
	int noc_l = nbocv_all_i[icv+1]-1;
	for (int noc = noc_f; noc <= noc_l; ++noc) {
	  FOR_J3 prcomm->packBufferDouble[3*ii+j] = A[noc][j];
	  ++ii;
	}
      }
      // rangeRotateVector is a virtual function that must be defined in 
      // any instantiated version of Gp...
      rangeRotateVector(prcomm->packBufferDouble+3*ii_f, // start
			ii-ii_f,                       // size
			range->getFlag(),range->dxyz);   // info about the transform
    }
    
    assert( ii == prcomm->npack_v );
    
  }
  
  // exchange...
  
  exchangePrcommBufferDoubleV(cvPrcommListG1,3);

  // unpack...

  for (list<Prcomm>::iterator prcomm = cvPrcommListG1.begin();
       prcomm != cvPrcommListG1.end(); ++prcomm) {
    
    int ii = 0;
    int nunpack = prcomm->unpackIndexVec.size();
    for (int i = 0; i < nunpack; ++i) {
      int icv = prcomm->unpackIndexVec[i];
      // should be level-0 ghost...
      assert( (icv >= ncv)&&(icv < ncv_g) );
      int noc_f = nbocv_all_i[icv];
      int noc_l = nbocv_all_i[icv+1]-1;
      for (int noc = noc_f; noc <= noc_l; ++noc) {
	FOR_J3 A[noc][j] = prcomm->unpackBufferDouble[3*ii+j];
	++ii;
      }
    }
    
    assert( ii == prcomm->nunpack_v );
    
  }
  
}

void UgpWithCv2::checkCvGrad() {
  
  if (mpi_rank == 0)
    cout << " > checkCvGrad()...";
  
  const double phi0 = 1.01;
  const double phix = 2.012;
  const double phiy = 3.0123;
  const double phiz = 4.01234;
  
  double * phi = new double[ncv_ggff];
  FOR_ICV_GGFF phi[icv] = phi0 + phix*x_cv[icv][0] + phiy*x_cv[icv][1] + phiz*x_cv[icv][2];
  
  double (*grad_phi)[3] = new double[ncv][3];
  calcCvGrad(grad_phi,phi);
    
  double my_max_error = 0.0;
  FOR_ICV {
    my_max_error = max( my_max_error, fabs(grad_phi[icv][0]-phix) ); 
    my_max_error = max( my_max_error, fabs(grad_phi[icv][1]-phiy) ); 
    my_max_error = max( my_max_error, fabs(grad_phi[icv][2]-phiz) ); 
  }
  
  double max_error;
  MPI_Reduce(&my_max_error,&max_error,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
  if (mpi_rank == 0)
    cout << "max_error (should be small): " << max_error <<endl;
  
  delete[] phi;
  delete[] grad_phi;

}

void UgpWithCv2::checkCv2Grad() {
  
  if (mpi_rank == 0)
    cout << " > checkCv2Grad()...";
  
  const double phi0 = 1.01;
  const double phix = 2.012;
  const double phiy = 3.0123;
  const double phiz = 4.01234;
  
  double * phi = new double[ncv_ggff];
  FOR_ICV_GGFF phi[icv] = phi0 + phix*x_cv[icv][0] + phiy*x_cv[icv][1] + phiz*x_cv[icv][2];
  
  double (*grad_phi)[3] = new double[ncv_g][3];
  calcCv2ScalGrad(grad_phi,phi);
    
  double my_max_error = 0.0;
  FOR_ICV_G {
    my_max_error = max( my_max_error, fabs(grad_phi[icv][0]-phix) ); 
    my_max_error = max( my_max_error, fabs(grad_phi[icv][1]-phiy) ); 
    my_max_error = max( my_max_error, fabs(grad_phi[icv][2]-phiz) ); 
  }
  
  double max_error;
  MPI_Reduce(&my_max_error,&max_error,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
  if (mpi_rank == 0)
    cout << "max_error (should be small): " << max_error <<endl;
  
  delete[] phi;
  delete[] grad_phi;

}

void UgpWithCv2::checkCvFilter() {
  
  if (mpi_rank == 0)
    cout << " > checkCvFilter()...";
  
  const double phi0 = 1.01;
  const double phix = 2.012;
  const double phiy = 3.0123;
  const double phiz = 4.01234;
  
  double * phi = new double[ncv_ggff];
  FOR_ICV_GGFF phi[icv] = phi0 + phix*x_cv[icv][0] + phiy*x_cv[icv][1] + phiz*x_cv[icv][2];
  
  double * phi_hat = new double[ncv];
  filterCvData(phi_hat,phi);
  
  double my_max_error = 0.0;
  FOR_ICV my_max_error = max( my_max_error, fabs(phi_hat[icv]-phi[icv]) ); 
    
  double max_error;
  MPI_Reduce(&my_max_error,&max_error,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
  if (mpi_rank == 0)
    cout << "max_error (should be small): " << max_error <<endl;
  
  delete[] phi;
  delete[] phi_hat;
  
}
  
void UgpWithCv2::checkCv2Filter() {
  
  if (mpi_rank == 0)
    cout << " > checkCv2Filter()...";
  
  const double phi0 = 1.01;
  const double phix = 2.012;
  const double phiy = 3.0123;
  const double phiz = 4.01234;
  
  double * phi = new double[ncv_ggff];
  FOR_ICV_GGFF phi[icv] = phi0 + phix*x_cv[icv][0] + phiy*x_cv[icv][1] + phiz*x_cv[icv][2];
  
  double * phi_hat = new double[ncv_g];
  filterCv2Data(phi_hat,phi);
  
  double my_max_error = 0.0;
  FOR_ICV_G my_max_error = max( my_max_error, fabs(phi_hat[icv]-phi[icv]) ); 
  
  double max_error;
  MPI_Reduce(&my_max_error,&max_error,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
  if (mpi_rank == 0)
    cout << "max_error (should be small): " << max_error <<endl;
  
  delete[] phi;
  delete[] phi_hat;
  
}




//######################################################################################
//######################################################################################
//######################################################################################
// put into another file
//######################################################################################
//######################################################################################
//######################################################################################

#include "Logging.h"
using namespace logging;



int UgpWithCv2Op::solveCvScalarBcgstab(double * phi, double * Ap, double * rhs,
    const double zero, const double zeroRel, const int maxiter, char *scalarName) {

  // we need the following work arrays...
  double *res = new double[ncv];
  double *res0 = new double[ncv];
  double *p = new double[ncv];
  double *v = new double[ncv];
  double *t = new double[ncv];
  double *s = new double[ncv];
  double *phat = new double[ncv_g];
  double *shat = new double[ncv_g];

  // initialize...
  for (int icv = 0; icv < ncv; icv++) {
    p[icv] = 0.0;
    v[icv] = 0.0;
  }
  for (int icv = 0; icv < ncv_g; icv++) {
    phat[icv] = 0.0;
    shat[icv] = 0.0;
  }

  double rho = 1.0;
  double alpha = 1.0;
  double omega = 1.0;

  // calculate the residual in rhs format...
  for (int icv = 0; icv < ncv; icv++)
  {
    int noc_f = nbocv_i[icv];
    int noc_l = nbocv_i[icv + 1] - 1;
    res[icv] = Ap[noc_f] * phi[icv]; // diagonal
    // cv neighbors...
    for (int noc = noc_f + 1; noc <= noc_l; noc++)
      res[icv] += Ap[noc] * phi[nbocv_v[noc]];

    res[icv] = rhs[icv] - res[icv];
    res0[icv] = res[icv];
  }

  double res0_max;
  double my_res0_max = 0.0;
  for (int icv = 0; icv < ncv; icv++)
    my_res0_max = max(my_res0_max, fabs(res[icv] / Ap[nbocv_i[icv]]));
  MPI_Reduce(&my_res0_max, &res0_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);



  int iter = 0;
  int done = 0;

  // check if we are already done
  if (mpi_rank == 0)
    if (res0_max <= zero)
    {
      cout << scalarName << ": " << iter << endl;
      done = 1;
    }
  MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

  while (done == 0)
  {
    iter++;

    double rho_prime = -omega * rho;

    double my_rho = 0.0;
    for (int icv = 0; icv < ncv; icv++)
      my_rho += res[icv] * res0[icv];
    MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (fabs(rho_prime) < 1.0E-20)
      rho_prime = 1.0E-20;
    double beta = alpha * rho / rho_prime;

    for (int icv = 0; icv < ncv; icv++)
      p[icv] = res[icv] - beta * (p[icv] - omega * v[icv]);

    // diagonal precon...
    for (int icv = 0; icv < ncv; icv++)
      phat[icv] = p[icv] / Ap[nbocv_i[icv]];
    updateCvDataG1(phat, REPLACE_DATA);

    // v = [A]{phat}
    for (int icv = 0; icv < ncv; icv++) {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;
      v[icv] = Ap[noc_f] * phat[icv]; // diagonal
      for (int noc = noc_f + 1; noc <= noc_l; noc++)
        v[icv] += Ap[noc] * phat[nbocv_v[noc]];
    }

    double my_gamma = 0.0;
    for (int icv = 0; icv < ncv; icv++)
      my_gamma += v[icv] * res0[icv];
    double gamma;
    MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (fabs(gamma) < 1.0E-20)
      gamma = 1.0E-20;

    alpha = rho / gamma;

    for (int icv = 0; icv < ncv; icv++)
      s[icv] = res[icv] - alpha * v[icv];

    // diagonal precon...
    for (int icv = 0; icv < ncv; icv++)
      shat[icv] = s[icv] / Ap[nbocv_i[icv]];
    updateCvDataG1(shat, REPLACE_DATA);

    // t = [A] shat...
    for (int icv = 0; icv < ncv; icv++)
    {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;
      t[icv] = Ap[noc_f] * shat[icv]; // diagonal
      for (int noc = noc_f + 1; noc <= noc_l; noc++)
        t[icv] += Ap[noc] * shat[nbocv_v[noc]];
    }

    double my_buf[2];
    my_buf[0] = my_buf[1] = 0.0;
    for (int icv = 0; icv < ncv; icv++)
    {
      my_buf[0] += s[icv] * t[icv];
      my_buf[1] += t[icv] * t[icv];
    }
    double buf[2];
    MPI_Allreduce(my_buf, buf, 2, MPI_DOUBLE, MPI_SUM, mpi_comm);

    omega = buf[0] / (buf[1] + 1.0E-20);

    // update phi...
    for (int icv = 0; icv < ncv; icv++)
      phi[icv] += alpha * phat[icv] + omega * shat[icv];
    updateCvDataG1(phi, REPLACE_DATA);

    double my_res_max = 0.0;

    // recompute the residual...
    for (int icv = 0; icv < ncv; icv++)
    {
      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;
      res[icv] = Ap[noc_f] * phi[icv]; // diagonal
      for (int noc = noc_f + 1; noc <= noc_l; noc++)
        res[icv] += Ap[noc] * phi[nbocv_v[noc]];
      // above is LHS. residual is then...
      res[icv] = rhs[icv] - res[icv];
      my_res_max = max(my_res_max, fabs(res[icv] / Ap[nbocv_i[icv]]));
    }

    // check if we are done...
    if (iter % 1 == 0)
    {
      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);

      if (mpi_rank == 0)
      {

        lout(DEBUG_HI) << scalarName << ":\t" << "bcgstab iter, res_max: " << iter << " " << res_max / (res0_max + 1.0E-12) << endl;

        if ((res_max <= zero) || (res_max/(res0_max+1.0E-12) <= zeroRel))
        {
          lout(INFO_LO) << scalarName << "\t iter/resAbs/resRel/tresholds = "
              << iter << "\t" << res_max << "\t" << res_max/(res0_max+1.0E-12) << "\t" << zero << "\t" << zeroRel << endl;
          done = 1;
        }

        if (iter > maxiter)
        {
          cout << "\nWarning: " << scalarName << " solveCvScalarBcgstab did not converge after " << maxiter << " iters, res_max: " << res_max << endl;
          done = 2;
        }
      }

      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
    }
  }


  delete[] res;
  delete[] res0;
  delete[] p;
  delete[] v;
  delete[] t;
  delete[] s;
  delete[] phat;
  delete[] shat;

  // let the calling routine know if we were successful...
  return (done == 1);

}




int UgpWithCv2Op::solveCvVectorR5Bcgstab(double(*phi)[5], double(*Ap)[5][5], double(*rhs)[5],
                                      const double zeroAbs, const double zeroRel, const int maxiter,
                                      char *scalarName) {

  // we need the following work arrays...
  static double (*res)[5] = new double[ncv][5];
  static double (*res0)[5] = new double[ncv][5];
  static double (*p)[5] = new double[ncv][5];
  static double (*v)[5] = new double[ncv][5];
  static double (*t)[5] = new double[ncv][5];
  static double (*s)[5] = new double[ncv][5];
  static double (*phat)[5] = new double[ncv_g][5];
  static double (*shat)[5] = new double[ncv_g][5];

  static double (*LUDEC)[5][5] = new double[ncv][5][5]; // LU decomposed diagonal for speed up
  double b[5];

  // LU decompose diagonal and store
  for (int icv = 0; icv < ncv; icv++) {
    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 5; j++)
        LUDEC[icv][i][j] = Ap[nbocv_i[icv]][i][j];

    ludeco(LUDEC[icv], 5);
  }

  // initialize...
  for (int icv = 0; icv < ncv; icv++)
    for (int i = 0; i < 5; i++) {
      p[icv][i] = 0.0;
      v[icv][i] = 0.0;
    }

  for (int icv = 0; icv < ncv_g; icv++)
    for (int i = 0; i < 5; i++) {
      phat[icv][i] = 0.0;
      shat[icv][i] = 0.0;
    }

  double rho = 1.0;
  double alpha = 1.0;
  double omega = 1.0;

  // calculate the residual in rhs format ...
  matTimesVecOverCVs(res, Ap, phi);

  for (int icv = 0; icv < ncv; icv++)
  for (int i = 0; i < 5; i++)
  {
    res[icv][i] = (rhs[icv][i] - res[icv][i]);
    res0[icv][i] = res[icv][i];
  }

  // compute first residual
  double res0_max, my_res0_max = 0.0;

  for (int icv = 0; icv < ncv; icv++)
  {
    double tmp[5];

    for (int i = 0; i < 5; i++)
      tmp[i] = res[icv][i];

    lusolv(LUDEC[icv], tmp, 5);

    for (int i = 0; i < 5; i++)
      my_res0_max = max(my_res0_max, fabs(tmp[i]));
  }

  MPI_Reduce(&my_res0_max, &res0_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);

  // start iteration

  int iter = 0;
  int done = 0;
  double res_max;

  while (done == 0) {

    iter++;

    double rho_prime = -omega*rho;

    double my_rho = 0.0;
    for (int icv = 0; icv < ncv; icv++)
      for (int i = 0; i < 5; i++)
        my_rho += res[icv][i] * res0[icv][i];

    MPI_Allreduce(&my_rho, &rho, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

    if (fabs(rho_prime) < 1.0E-20)
      rho_prime = 1.0E-20;

    double beta = alpha * rho / rho_prime;

    for (int icv = 0; icv < ncv; icv++)
      for (int i = 0; i < 5; i++)
        p[icv][i] = res[icv][i] - beta * (p[icv][i] - omega * v[icv][i]);

    // block-diagonal precon...
    // solve the system
    for (int icv = 0; icv < ncv; icv++) {
      for (int i = 0; i < 5; i++)
        phat[icv][i] = p[icv][i];
      lusolv(LUDEC[icv], phat[icv], 5);
    }

    UpdateCvDataStateVec(phat);

    // v = [A]{phat}
    matTimesVecOverCVs(v, Ap, phat);

    double my_gamma = 0.0;

    for (int icv = 0; icv < ncv; icv++)
      for (int i = 0; i < 5; i++)
        my_gamma += v[icv][i] * res0[icv][i];

    double gamma;
    MPI_Allreduce(&my_gamma, &gamma, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

//    gamma = my_gamma;

    if (fabs(gamma) < 1.0E-20)
      gamma = 1.0E-20;

    alpha = rho/gamma;

    for (int icv = 0; icv < ncv; icv++)
      for (int i = 0; i < 5; i++)
        s[icv][i] = res[icv][i] - alpha * v[icv][i];

    // diagonal precon...
    for (int icv = 0; icv < ncv; icv++) {
      for (int i = 0; i < 5; i++)
        shat[icv][i] = s[icv][i];
      lusolv(LUDEC[icv], shat[icv], 5);
    }

    UpdateCvDataStateVec(shat);

    // t = [A] shat...
    matTimesVecOverCVs(t, Ap, shat);

    double my_buf[2];
    my_buf[0] = my_buf[1] = 0.0;
    for (int icv = 0; icv < ncv; icv++)
      for (int i = 0; i < 5; i++) {
        my_buf[0] += s[icv][i] * t[icv][i];
        my_buf[1] += t[icv][i] * t[icv][i];
      }

    double buf[2];
    MPI_Allreduce(my_buf, buf, 2, MPI_DOUBLE, MPI_SUM, mpi_comm);

    omega = buf[0] / (buf[1] + 1.0E-20);

    // update phi...
    for (int icv = 0; icv < ncv; icv++)
      for (int i = 0; i < 5; i++)
        phi[icv][i] += alpha * phat[icv][i] + omega * shat[icv][i];

    UpdateCvDataStateVec(phi);


    double my_res_max = 0.0;

    // recompute the residual...
    matTimesVecOverCVs(res, Ap, phi);

    for (int icv = 0; icv < ncv; icv++)           // above is LHS. residual is then...
    for (int i = 0; i < 5; i++)
      res[icv][i] = rhs[icv][i] - res[icv][i];

    // compute normalized res
    for (int icv = 0; icv < ncv; icv++)
    {
      double tmp[5];

      for (int i = 0; i < 5; i++)
        tmp[i] = res[icv][i];

      lusolv(LUDEC[icv], tmp, 5);

      for (int i = 0; i < 5; i++)
        my_res_max = max(my_res_max, fabs(tmp[i]));
    }

    MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);

    // check residual
    if (mpi_rank == 0)
    {
//      cout << "iter: " << iter << " " << res_max << " " << res_max / (res0_max + 1.0E-12) << endl;

      if ((res_max <= zeroAbs) || (res_max / (res0_max + 1.0E-12) <= zeroRel))
        done = 1;
      else if (iter > maxiter)
      {
/*        cout << "Warning: " << scalarName
             << " solveCvScalarBcgstab did not converge after " << maxiter
             << " iters, res_max: " << res_max << endl;*/

        done = 2;
      }
    }
    MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
  }

  lout(INFO_LO) << "\tbcgstab: iter: " << iter << " res_max: " << res_max << endl;

  /*    delete[] res;
   delete[] res0;
   delete[] p;
   delete[] v;
   delete[] t;
   delete[] s;
   delete[] phat;
   delete[] shat;

   delete [] LU;*/

  // let the calling routine know if we were successful...
  return (done == 1);
}





