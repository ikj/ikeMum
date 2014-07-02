#include "JoeWithModels.h"
#include "ADJOINT_FILES/JoeWithModelsAD.h"
#include "adolc.h"
#include "turbModels/TurbModel_KOMSST.h"
#include "ADJOINT_FILES/TurbModel_KOMSST_AD.h"



// ###########################################################################################
// ------                                                                               ------
// ------                    Joe for mixing layer (binary mixing)                       ------
// ------                                                                               ------
// ###########################################################################################
class MyJoe_flow : public JoeWithModels, public RansTurbKOmSST
{
protected:
  // identifiers for inlet profiles
#define NO_INIT_FILE "NO_INIT_FILE"
#define NO_BC_FILE "NO_BC_FILE"
#define NO_TRANSFORM_MESH_FILE "NO_TRANSFORM_MESH_FILE"

// memory for initial conditions
#define NVALUES 4
  double (*initVal)[NVALUES];

// memory for boundary conditions (inlet)
#define NVALUES_BC 6
  double (*bcVal)[NVALUES_BC];


public:
   int ownID;
   double *Parent;

public:

  MyJoe_flow(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name), ownID(0)
  {
    if (mpi_rank == 0)      cout << "MyJoe_flow()" << endl;
     Parent = NULL;
     registerScalar(Parent,  "Parent",  CV_DATA);
    init();
  }


  virtual ~MyJoe_flow()  {}



    void initialHook() {JoeWithModels::initialHook();

FILE *fp;
    string fname;
    int nPos;

    // initialization

    fname = getStringParam("INIT_FILE", NO_INIT_FILE);
    if (fname == NO_INIT_FILE) {
      cout << "no initial file specified. skip initialization " << endl;
    } else { // DO INITIALIZE FROM FILE

      if ((fp=fopen(fname.c_str(),"r")) == NULL) {
        cout << "could not open file " << fname << endl;
        throw(-1);
      }

      fscanf(fp, "n=%d", &nPos);
      cout << "read " << nPos << " points" << endl;

      initVal = new double[nPos][NVALUES];

      // file has values
      // x, rhou_x, rhoE, rho   
      for (int i=0; i<nPos; i++)
        for (int v=0; v<NVALUES; v++)
          fscanf(fp, "%lf", &initVal[i][v]);

      fclose(fp);

      // apply initial conditions

      for (int icv = 0; icv < ncv_ggff; icv++) {
        int pos = 0;
        while(initVal[pos][0] <= x_cv[icv][0]) pos++;
        double f = (x_cv[icv][0]-initVal[pos-1][0])/(initVal[pos][0]-initVal[pos-1][0]);

        double rhou_x1 = initVal[pos-1][1];
        double rhou_x2 = initVal[pos  ][1];
        double my_rhou_x = rhou_x1 + f*(rhou_x2 - rhou_x1);
        rhou[icv][0] = my_rhou_x;

        double rhoe1 = initVal[pos-1][2];
        double rhoe2 = initVal[pos  ][2];
        double my_rhoe = rhoe1 + f*(rhoe2 - rhoe1);
        rhoE[icv] = my_rhoe;

        double rho1 = initVal[pos-1][3];
        double rho2 = initVal[pos  ][3];
        double my_rho  = rho1 + f*(rho2 - rho1);
        rho[icv] = my_rho;
      }
    }
    // read inlet boundary condition

 

 fname = getStringParam("BC_FILE", NO_BC_FILE);
    if (fname == NO_BC_FILE) {
       cout << "no bc file specified." << endl;
    } else { // DO READ BC FILE

      if ((fp=fopen(fname.c_str(),"r")) == NULL) {
        cout << "could not open file " << fname << endl;
        throw(-1);
      }

      fscanf(fp, "%d", &nPos);
      cout << "read " << nPos << "points" << endl;

      bcVal = new double[nPos][NVALUES_BC];

      // file has values
      // y, vel_x, vel_y, vel_z, press, temp
      for (int i=0; i<nPos; i++)
        for (int v=0; v<NVALUES_BC; v++)
          fscanf(fp, "%lf", &bcVal[i][v]);

      fclose(fp);
    }
}

  

  virtual void initialHookScalarRansTurbModel()
  {
	  RansTurbKOmSST::initialHookScalarRansTurbModel();
  }
  
  void temporalHook()
  {
    if (step%10 == 0)
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
        if (zone->getKind() == FA_ZONE_BOUNDARY)
	  {
	    Param *param;
	    if (getParam(param, zone->getName()))
	      if (param->getString() == "WALL")
		writeWallValues(zone->getName());
	  } 
  }
  
  void finalHook()
  {
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
        {
          Param *param;
          if (getParam(param, zone->getName()))
            if (param->getString() == "WALL")
              writeWallValues(zone->getName());
        }
  }


// inlet Hook
  virtual void boundaryHook(double *temp, double (*vel)[3], double *press, FaZone *zone)
  {
    if (zone->getNameString() == "inlet")   // if more HOOK boundaries are defined
    {
	 for (int index = 0; index < zone->faVec.size(); ++index) {
              int ifa = zone->faVec[index];
              int icv1 = cvofa[ifa][1];

        int pos=0;
        while(bcVal[pos][0] < x_fa[ifa][1]) pos++;
        double f = (x_fa[ifa][1]-bcVal[pos-1][0])/(bcVal[pos][0]-bcVal[pos-1][0]);
        
        double u1 = bcVal[pos-1][1];
        double u2 = bcVal[pos  ][1];
        double u = u1 + f*(u2-u1);
        vel[icv1][0] = u;
        
        double v1 = bcVal[pos-1][2];
        double v2 = bcVal[pos  ][2];
        double v = v1 + f*(v2-v1);
        vel[icv1][1] = v;
        
        double w1 = bcVal[pos-1][3];
        double w2 = bcVal[pos  ][3];
        double w = w1 + f*(w2-w1);
        vel[icv1][2] = w;
        
        double p1 = bcVal[pos-1][4];
        double p2 = bcVal[pos  ][4];
        double p = p1 + f*(p2-p1);
        press[icv1] = p;
        
        double t1 = bcVal[pos-1][5];
        double t2 = bcVal[pos  ][5];
        double t = t1 + f*(t2-t1);
        temp[icv1] = t; 
      } 
    } 
  } 



  
  /*
   * write wall values to file in tecplot format
   */
  virtual void writeWallValues(string name)
  {
    // count the walls for each rank
    int my_wall_faces = 0;
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
        if (zone->getNameString() == name)
          for (int ifa = zone->ifa_f; ifa<=zone->ifa_l; ifa++)
            my_wall_faces++;

    int tot_wall_faces = 0;
    MPI_Allreduce(&my_wall_faces, &tot_wall_faces, 1, MPI_INT, MPI_SUM, mpi_comm);

    // write to file in tecplot format
    FILE *fp;
    char fname[200];
    sprintf(fname, "%s.dat", name.c_str());
    if ( mpi_rank == 0 )
    {
      if ( (fp=fopen(fname,"wt"))==NULL )
      {
        cerr << "Error: cannot open file " << fname << endl;
        throw(-1);
      }
      fprintf(fp, "VARIABLES =\"X\" \"Y\" \"Z\" \"rho\" \"press\" \"temp\" \"tau_wall\" \"qdot\" \"yplus\" \"muLam\"\n");
      fprintf(fp, "Zone T =\"wall\" I = %d , F = point\n", tot_wall_faces);
    }
    else
    {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
      if ( (fp=fopen(fname,"a"))==NULL )
      {
        cerr << "Error: cannot open file " << fname << endl;
        throw(-1);
      }
    }

    for (list<FaZone>::iterator faZone = faZoneList.begin(); faZone != faZoneList.end(); faZone++)  // loop over all boundary zones
      if (faZone->getKind() == FA_ZONE_BOUNDARY)
        if (faZone->getNameString() == name)                            // identify zone = name
        {
          for (int ifa = faZone->ifa_f; ifa <= faZone->ifa_l; ifa++)    // loop over boundary faces
          {
            int icv = cvofa[ifa][0];                                    // get cv of boundary face

            double n[3], s_half[3], vel[3], velTang[3];

            normVec3d(n, fa_normal[ifa]);                               // get area weighted face normal (outward pointing)
            vecMinVec3d(s_half, x_fa[ifa], x_cv[icv]);

            vel[0] = rhou[icv][0]/rho[icv];
            vel[1] = rhou[icv][1]/rho[icv];
            vel[2] = rhou[icv][2]/rho[icv];
            double un = vecDotVec3d(n, vel);

            velTang[0] = vel[0] - un*n[0];
            velTang[1] = vel[1] - un*n[1];
            velTang[2] = vel[2] - un*n[2];

            double velMag = sqrt(vecDotVec3d(velTang, velTang));
            double walld = fabs(vecDotVec3d(s_half, n));

            double wallTemp = 300.0;
            double mulam = mul_fa[ifa];

            double tau = mulam*velMag/walld*sign(vel[0]);
            double utau = sqrt(fabs(tau)/rho[icv]);
            double yplus = walld*utau/(mulam/rho[icv]);

            double kTotal = R_gas*gamma[icv]/(gamma[icv] - 1.0)*mulam/Pr;
            double qDot = kTotal*(temp[icv]-wallTemp)/walld;

            fprintf(fp, "%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\n",
                x_fa[ifa][0], x_fa[ifa][1], x_fa[ifa][2], rho[icv], press[icv], temp[icv], tau, qDot, yplus, mul_fa[ifa]);
          }
        }

    fclose(fp);

    if ( mpi_rank < mpi_size-1 )
    {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }
    MPI_Barrier(mpi_comm);
  }






};



class MyJoe_adj : public JoeWithModels_AD, public RansTurbKOmSST_AD
{ 
protected:
// identifiers for inlet profiles
  int nPos, nval;

#define NVALUES_BC 6
  double (*bcVal)[NVALUES_BC];

public:
   int ownID;
   double *Parent;

public:

  MyJoe_adj(char *name) : JoeWithModels(name), JoeWithModels_AD(name), UgpWithCvCompFlow(name), ownID(0) {

    if (mpi_rank == 0)      cout << "MyJoe_adj()" << endl;
    Parent = NULL;
    registerScalar(Parent,  "Parent",  CV_DATA);

    init();

    FILE *fp;
    string fname;

    fname = getStringParam("BC_FILE", NO_BC_FILE);
    if (fname == NO_BC_FILE) {
       cout << "no bc file specified." << endl;
    } else { // DO READ BC FILE

      if ((fp=fopen(fname.c_str(),"r")) == NULL) {
        cout << "could not open file " << fname << endl;
        throw(-1);
      }

      fscanf(fp, "%d", &nPos);
      cout << "read " << nPos << "points" << endl;

      bcVal = new double[nPos][NVALUES_BC];

      // file has values
      // y, vel_x, vel_y, vel_z, press, temp
      for (int i=0; i<nPos; i++)
        for (int v=0; v<NVALUES_BC; v++)
          fscanf(fp, "%lf", &bcVal[i][v]);

      fclose(fp);
    }


  }


  virtual ~MyJoe_adj() {}

    


//*************** User defined objective functional ************
virtual void calcFunctional_AD(adouble *rho, adouble (*rhou)[3], adouble *rhoE)
//**************************************************************
{
      int icv, totcv=0;
      functional=0.;

        //double xIn  =  0.2;
        //double xOut =  0.25;
	//A
        //double yIn  =  0.0036;
        //double yOut =  0.0037;
	//B
        //double yIn  =  0.0025;
        //double yOut =  0.0045;

        double xIn  = getDoubleParam("XIN",  "0.2");
        double xOut = getDoubleParam("XOUT", "0.25");
        double yIn  = getDoubleParam("YIN",  "0.0036");
        double yOut = getDoubleParam("YOUT", "0.0037");

        double Denominator = 0.;
        double Numerator   = 0.;
        adouble LocalFuncN  = 0.;
        adouble LocalFuncD  = 0.;

	adouble *muTdvdy;
	muTdvdy   = new adouble [ncv_gg];


	int kine_index  = getScalarTransportIndex("kine");
	int omega_index = getScalarTransportIndex("omega");

          for (int icv = 0; icv < ncv_g; icv++)
          {
            adouble vely0 = rhou[icv][1]/rho[icv];

              double volume = cv_volume[icv];

              adouble kine_cv    = scalarTranspEqVector_AD[kine_index].phi[icv];
              adouble omega_cv   = scalarTranspEqVector_AD[omega_index].phi[icv];

	      adouble muT = rho[icv]*kine_cv/omega_cv;

		double sum_weights = 0.;
		int num_weights = 0;
		adouble dvdy = 0.;

		int noc_f = nbocv_all_i[icv];
		int noc_l = nbocv_all_i[icv+1]-1;

		for (int noc = noc_f+1; noc <= noc_l; ++noc) {

		int icv_nbr = nbocv_all_v[noc];
		double ds[3]={0.,0.,0.};
		vecMinVec3d(ds, x_cv[icv_nbr], x_cv[icv]);

		  if(abs(ds[1])>1.e-6) {
            		adouble vely = rhou[icv_nbr][1]/rho[icv_nbr];
   			dvdy += (vely-vely0)/ds[1]*(ds[1]*ds[1]);
   			sum_weights += ds[1]*ds[1];
   			num_weights++;
 		  }

		}

    		dvdy /= sum_weights; 

		muTdvdy[icv] = muT*dvdy; 
		//muTdvdy[icv] = press[icv]; 

       }

          for (int icv = 0; icv < ncv; icv++)
          {

		if(x_cv[icv][0] > xIn && x_cv[icv][0] < xOut && x_cv[icv][1] > yIn && x_cv[icv][1] < yOut) {

            	adouble q0 = muTdvdy[icv];

		LocalFuncN    += pow(q0,8)*cv_volume[icv]*x_cv[icv][0];
		LocalFuncD    += pow(q0,8)*cv_volume[icv];
                Denominator   += pow(q0.value(),8)*cv_volume[icv];
                Numerator     += pow(q0.value(),8)*cv_volume[icv]*x_cv[icv][0];

                 totcv++;

	 	//cout<<icv<<" "<<x_cv[icv][0]<<" "<<x_cv[icv][1]<<" "<<q0.value()<<" "<<dqdx.value()<<" "<<num_weights<<endl; 
		}
          }


	         double functional_local,functional_total,Denominator_total,Numerator_total;
	         int    totcv_total;
	         MPI_Allreduce(&Denominator,           &Denominator_total,  1, MPI_DOUBLE, MPI_SUM, mpi_comm);
	         MPI_Allreduce(&Numerator,             &Numerator_total,    1, MPI_DOUBLE, MPI_SUM, mpi_comm);

            	 functional = LocalFuncN/Denominator_total - LocalFuncD*Numerator_total/Denominator_total/Denominator_total;

            	 functional_local = LocalFuncN.value()/Denominator_total;

	         MPI_Allreduce(&functional_local, &functional_total,1, MPI_DOUBLE, MPI_SUM, mpi_comm);
	         MPI_Allreduce(&totcv,            &totcv_total,     1, MPI_INTEGER,MPI_SUM, mpi_comm);

         
		if(mpi_rank==0) {
               cout<<"Total CVs= "<<totcv_total<<endl;
	       printf("Denominator  = %.16e\n", Denominator_total);               
	       printf("Functional_total  = %.16e\n", functional_total);               
               }
}


/*

//*** You will have to write your own routine for other gradients**
virtual void calcFunctionalGradient(double *adj_vars)
//**************************************************************
// dJ/d U_fs =(V)^T . dR/d U_fs
{
      int nScal = scalarTranspEqVector.size();
      int nAdjvars = 5+nScal;

      int icv, totcv_local=0;
      double functional_gradient_local =0.;

      double functional_gradient_global =0.;
      int totcv_global = 0;

          for (int icv= 0; icv <ncv; icv++)
          {
        adouble TS = 1.0;
        if (KOM_RealizableConstraint == 2)
          TS = min(1.0/(omega[icv]), 0.6/(sqrt(6.0)*strMag[icv]));
    
        adouble src =  TS*omega[icv]/kine[icv]*calcTurbProd_AD(icv);
    
        double delta_gradient  = src.value()*cv_volume[icv]*adj_vars[icv*nAdjvars+6];
      
        functional_gradient_local += delta_gradient;
             
        totcv_local++;
        
           }
        
               MPI_Allreduce(&functional_gradient_local, &functional_gradient_global, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
               MPI_Allreduce(&totcv_local, &totcv_global, 1, MPI_INT, MPI_SUM, mpi_comm);
                                
               if(mpi_rank==0) {
               cout<<"Total CVs= "<<totcv_global<<endl;
               printf("Functional Gradient  = %.16e\n", functional_gradient_global);
               }                
                                

}               
*/


//*** You will have to write your own routine for other gradients**
virtual void calcFunctionalGradient(double *adj_vars)
//**************************************************************
// dJ/d U_fs =(V)^T . dR/d U_fs
{
      int nScal = scalarTranspEqVector.size();
      int nAdjvars = 5+nScal;

      int icv, totcv_local=0;
      double functional_gradient_local =0.;

      double functional_gradient_global =0.;
      int totcv_global = 0;

          for (int icv= 0; icv <ncv; icv++)
          {
      adouble F1 = blendFuncF1_AD[icv];
      adouble dalfa = F1;

      adouble zeta = min(1.0/omega[icv], a1/(limiterFunc_AD[icv]*blendFuncF2_AD[icv]));
      adouble mut = min(max(rho_AD[icv]*kine[icv]*zeta, 1.0e-8), 1.0e5);

      adouble src = dalfa*rho_AD[icv]/mut*calcTurbProd_AD(icv,SST_limitPK);

        double delta_gradient  = src.value()*cv_volume[icv]*adj_vars[icv*nAdjvars+6];

        functional_gradient_local += delta_gradient;

        totcv_local++;

           }

               MPI_Allreduce(&functional_gradient_local, &functional_gradient_global, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
               MPI_Allreduce(&totcv_local, &totcv_global, 1, MPI_INT, MPI_SUM, mpi_comm);

               if(mpi_rank==0) {
               cout<<"Total CVs= "<<totcv_global<<endl;
               printf("Functional Gradient  = %.16e\n", functional_gradient_global);
               }

}

  virtual void boundaryHook_AD(adouble *temp, adouble (*vel)[3], adouble *press, FaZone *zone)
  {
    if (zone->getNameString() == "inlet")   // if more HOOK boundaries are defined
    {
	 for (int index = 0; index < zone->faVec.size(); ++index) {
              int ifa = zone->faVec[index];
              int icv1 = cvofa[ifa][1];

        int pos=0;
        while(bcVal[pos][0] < x_fa[ifa][1]) pos++;
        double f = (x_fa[ifa][1]-bcVal[pos-1][0])/(bcVal[pos][0]-bcVal[pos-1][0]);
        
        double u1 = bcVal[pos-1][1];
        double u2 = bcVal[pos  ][1];
        double u = u1 + f*(u2-u1);
        vel[icv1][0] = u;
        
        double v1 = bcVal[pos-1][2];
        double v2 = bcVal[pos  ][2];
        double v = v1 + f*(v2-v1);
        vel[icv1][1] = v;
        
        double w1 = bcVal[pos-1][3];
        double w2 = bcVal[pos  ][3];
        double w = w1 + f*(w2-w1);
        vel[icv1][2] = w;
        
        double p1 = bcVal[pos-1][4];
        double p2 = bcVal[pos  ][4];
        double p = p1 + f*(p2-p1);
        press[icv1] = p;
        
        double t1 = bcVal[pos-1][5];
        double t2 = bcVal[pos  ][5];
        double t = t1 + f*(t2-t1);
        temp[icv1] = t; 
      } 
    } 
  } 


/*

//*** Routine to evaluate derivative wrt to Freestream velocity****
//*** You will have to write your own routine for other gradients**
virtual void calcFunctionalGradient(double *adj_vars)
//**************************************************************
// dJ/d U_fs =(V)^T . dR/d U_fs
// In this case (supersonic inflow), dR/d U_fs = d (flux_freestream)/d U_fs
{
      int nScal = scalarTranspEqVector.size();
      int nAdjvars = 5+nScal;
      
      int icv, totcv_local=0;
      double functional_gradient_local =0.;
      
      double functional_gradient_global =0.;
      int totcv_global = 0;
      
          for (int icv= 0; icv <ncv; icv++)
          {
          
        adouble mue       = InterpolateAtCellCenterFromFaceValues(mul_fa, icv);
        double Om         = max(vortMag[icv].value(), 1.e-8);   // vorticity
        double Strainrate = max(strMag[icv].value(), 1.e-8);    // strainrate
        double chi   = nuSA[icv].value()*rho[icv]/mue.value();
        double chi_3 = pow(chi, 3.0);
        double d = wallDist[icv];

        double inv_KarmanConst2_d2 = 1.0/(KarmanConst_2*d*d);

        double fv1 = chi_3/(chi_3 + cv1_3);
        double fv2 = 1.0 - chi/(1.0 + chi*fv1);

        double BlendTurbulentProduction = 0.0;
        double S = max(Om + min(0.0, Strainrate - Om)*BlendTurbulentProduction, 0.0);
        double Shat = S + nuSA[icv].value()*fv2*inv_KarmanConst2_d2;

        double inv_Shat = 1.0/max(Shat, 1.0e-8);

	     // model sources
        double dH1bydcb1 = rho[icv]*Shat*nuSA[icv].value();

        double delta_gradient  = cv_volume[icv]*dH1bydcb1*adj_vars[icv*nAdjvars+5];

        functional_gradient_local += delta_gradient;

        totcv_local++;

           }

               MPI_Allreduce(&functional_gradient_local, &functional_gradient_global, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
               MPI_Allreduce(&totcv_local, &totcv_global, 1, MPI_INT, MPI_SUM, mpi_comm);

               if(mpi_rank==0) {
               cout<<"Total CVs= "<<totcv_global<<endl;
               printf("Functional Gradient  = %.16e\n", functional_gradient_global);
               }



}

*/



  
};



int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  
  initMpiStuff();

  int run = 0;

  // set input name to default "Joe.in"
  char inputFileName[50];
  sprintf(inputFileName, "Joe.in");

  for (int i=1; i<argc; i++)
  {
    string str(argv[i]);
    if (from_string<int>(run, str, std::dec))
    {
      if (mpi_rank == 0)
        cout << "You have specified run number = " << run << endl;
    }
    else
      strcpy(inputFileName, argv[i]);
  }

  if (mpi_rank == 0)
  {
    cout << "SPECIFIED INPUT NAME = " << inputFileName << endl;
    cout << "SPECIFIED RUN = " << run << endl;
  }


  try {

    
    // provide total runtime 
    double wtime, wtime0;
    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0)
      wtime = MPI_Wtime();

        if(run==0) {
    	  JoeWithModels *joe;
    	  joe = new MyJoe_flow(inputFileName);    
          joe->run();
          delete joe;
        }
        else {
    	  JoeWithModels_AD *joe;
    	  joe = new MyJoe_adj(inputFileName);    
          joe->runAdjoint();
    	  delete joe;
        }

    

    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0) 
    {
      double wtime0 = wtime;
      wtime = MPI_Wtime();
      cout << " > total runtime [s]: " << wtime - wtime0 << endl;
    }

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
