/**
 * \brief Chemtable class for chemistry table: load table and lookup value based on 3 coordinates.
 * 
 *  This version of the chemistry table is based on a 3D structured Cartesian tabulation,
 *  i.e., the 3 dimensions are simple vectors
 * 
 * \author Vincent Terrapon 
 * \date August 2009
 * \version 1.1
 */

#ifndef CHEMISTRYTABLECARTESIAN_H
#define CHEMISTRYTABLECARTESIAN_H

#include "CombustionGeneralDefinitions.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////

/*! \brief Class containing the chemistry table and associated functions.
 * 
 *  Read and load 3D table from chemistry table input file.
 *  Lookup values based on 3 dimensions by interpolating structured Cartesian table
 *  using simple linear interpolation 
 *  
 *  Class assumes a structured Cartesian table: x1=x1(i), x2=x2(j), x3=x3(k).
 */
//template <class Chemtable>
class ChemtableCartesian
{
private:
  int     n1, n2, n3;                   ///< Dimensions of the table (corresponding to Zmean, Zvar and chi/progress variable)
  int     nvar;                         ///< Number of variables stored in table
  char    CombustionModel[kStrMed];     ///< Combustion model used: Steady Flamelet or FPVA
  char    (*VarNames)[kStrMed];         ///< Names of the stored variables
  double  *x1, *x2, *x3;                ///< Axis coordinates of the table (Zmean, Zvar, chi/progress variable)
  int     ***iMask;                     ///< Indicating cells of table with no physical meaning: iMask(k,j,i) where k->x3, j->x2, i->x1
  bool    is_iMask;                     ///< True if iMask (blanked cells) should be loaded
  double  ****Data;                     ///< Actual table data: Data(l,k,j,i) where l->variable, k->x3, j->x2, i->x1
  double  Coord1, Coord2, Coord3;       ///< Coordinates used for interpolation.
  InterpolationIndex interp;            ///< Interpolation indices and weights (set by SetCoordinates once to avoid recomputing it for each species).
  
public:
  map<string, int> myTableVar;                   ///< Map containing the species linking the name (string) to table first index (int)
  map<string, int>::iterator itTp;               ///< Iterators to loop over the species
  pair<map<string, int>::iterator, bool> retT;   ///< Return value for map insertion to check if species is already in map
  
  
  ChemtableCartesian()
  {
    if (mpi_rank == 0)
      cout << "ChemtableCartesian()" << endl;
    
    x1 = NULL;
    x2 = NULL;
    x3 = NULL;
    iMask = NULL;
    Data = NULL;
  }
  
  /* Accessors */
  /*************/
  
  int GetChemtableDimension1() {return n1;}
  int GetChemtableDimension2() {return n2;}
  int GetChemtableDimension3() {return n3;}
  int GetChemtableDimension4() {return nvar;}

  double GetChemtableCoordinate1(int i1, int i2, int i3) {return x1[i1];}
  double GetChemtableCoordinate2(int i1, int i2, int i3) {return x2[i2];}
  double GetChemtableCoordinate3(int i1, int i2, int i3) {return x3[i3];}

  double GetChemtableValue(int i1, int i2, int i3, int ivar) {return Data[ivar][i3][i2][i1];}


  /***********************************************************************************************************/

  /*! \brief Open, read and load chemistry table.
   *
   *  Read input file to determine table file name and some parameters, then read table from file and allocate 
   *  large array to contain table. Table is a 4D array with each species contiguous in memory.
   *  Also output information on table for verification.
   *  
   *  Table might use blanking if created with \p CreateChemtable tool from NGA. Blanking are dummy values
   *  if table is created with tool \p CreateChemTable in the tool directory.
   *  
   *  \param myInputPtr A pointer of type ParamMap pointing on the general input file.
   */
  void Load(ParamMap *myInputPtr)
  {
    FILE *inFp = 0;
    int i, j, k, l, mi1, mi2, mi3, ma1, ma2, ma3;
    string key, buffer_str;
    char filename[kStrLong], buffer_m[kStrMed];
    char *p = 0, dummy;
    bool is_iMask;
    size_t dum;
    double minval, maxval;

    /* Read input file for chemistry table information */
    buffer_str = myInputPtr->getStringParam("CHEMTABLE_FILE");
    dum = buffer_str.copy(&filename[0], kStrLong);
    dum = buffer_str.size();
    if (dum < kStrLong)
      strcpy(&filename[dum], "\0");
    is_iMask = (bool) myInputPtr->getIntParam("LOAD_MASK");

    /* Read chemistry table file */
    if (mpi_rank == 0)
    {
      cout << endl << "Chemtable is " << filename << endl;
    }
    if (!(inFp = fopen(filename, "rb")))
    {
      cerr << "### Could not open input file " << filename << " ###\n";
      throw(-1);
    }
    /* Dimensions */
    dum = fread(&n1, sizeof(int), 1, inFp);
    if (dum != 1)
    {
      cerr << "### Error reading file: n1 ###" << endl;
      throw(-1);
    }
    dum = fread(&n2, sizeof(int), 1, inFp);
    if (dum != 1)
    {
      cerr << "### Error reading file: n2 ###" << endl;
      throw(-1);
    }
    dum = fread(&n3, sizeof(int), 1, inFp);
    if (dum != 1)
    {
      cerr << "### Error reading file: n3 ###" << endl;
      throw(-1);
    }
    dum = fread(&nvar, sizeof(int), 1, inFp);
    if (dum != 1)
    {
      cerr << "### Error reading file: nvar ###" << endl;
      throw(-1);
    }

    /* Coordinates */
    getMem1D(&x1, 0, n1 - 1, "Chemtable::Load x1", true);
    getMem1D(&x2, 0, n2 - 1, "Chemtable::Load x2", true);
    getMem1D(&x3, 0, n3 - 1, "Chemtable::Load x3", true);
    dum = fread(x1, sizeof(double), n1, inFp);
    if (dum != n1)
    {
      cerr << "### Error reading file: x1 ###" << endl;
      throw(-1);
    }
    dum = fread(x2, sizeof(double), n2, inFp);
    if (dum != n2)
    {
      cerr << "### Error reading file: x2 ###" << endl;
      throw(-1);
    }
    dum = fread(x3, sizeof(double), n3, inFp);
    if (dum != n3)
    {
      cerr << "### Error reading file: x3 ###" << endl;
      throw(-1);
    }

    /* Masks */
    getMem3D(&iMask, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1, "Chemtable::Load iMask", false);
    for (k = 0; k < n3; k++)
    {
      for (j = 0; j < n2; j++)
      {
        for (i = 0; i < n1; i++)
        {
          dum = fread(&iMask[k][j][i], sizeof(int), 1, inFp);
          if (dum != 1)
          {
            cerr << "### Error reading file: mask (" << i << j << k
            << ") ###" << endl;
            throw(-1);
          }
        }
      }
    }
    if (!(is_iMask))
    {
      freeMem3D(iMask, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1);
      iMask = NULL;
    } /* iMask not used, save memory */

    /* Other input data */
    for (i = 0; i < kStrMed - 1; i++)
    {
      strcpy(&CombustionModel[i], " ");
    }
    dum = fread(&CombustionModel[0], sizeof(char), kStrMed, inFp);
    if (dum != kStrMed)
    {
      cerr << "### Error reading file: Combustion model ###" << endl;
      throw(-1);
    }
    strcpy(&CombustionModel[kStrMed - 1], "\0");

    VarNames = new char[nvar][kStrMed];
    for (l = 0; l < nvar; l++)
    {
      strcpy(&VarNames[l][0], "\0");
      dum = fread(&buffer_m, sizeof(char), kStrMed, inFp);
      if (dum != kStrMed)
      {
        cerr << "### Error reading file: VarNames no. " << l << " ###"
        << endl;
        throw(-1);
      }
      p = strchr(buffer_m, ' ');
      if (p != 0)
        strcpy(p, "\0");
      strcpy(&VarNames[l][0], buffer_m);
      /* Create map to link variable names and index */
      key.assign(buffer_m);
      retT = myTableVar.insert(pair<string, int> (string(key), l));
      if (retT.second == false)
      {
        cout << "### Variable " << key
        << " already in myTableVar in loop ####" << endl;
      }
    }

    /* Table data */
    getMem4D(&Data, 0, nvar - 1, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1,
        "Chemtable::Load Data", false);
    for (l = 0; l < nvar; l++)
    {
      for (k = 0; k < n3; k++)
      {
        for (j = 0; j < n2; j++)
        {
          for (i = 0; i < n1; i++)
          {
            dum = fread(&Data[l][k][j][i], sizeof(double), 1, inFp);
            if (dum != 1)
            {
              cerr << "### Error reading file: data (" << i << j << k
              << l << ") ### " << endl;
              throw(-1);
            }
          }
        }
      }
    }
    fclose(inFp);

    /* Output useful information on table */
    if (mpi_rank == 0)
    {
      cout << "Chemistry table is structured Cartesian" << endl;
      if (is_iMask)
      {
        cout << "Chemistry table loaded with masks" << endl;
      }
      else
      {
        cout << "Chemistry table loaded without masks" << endl;
      }
      cout << "The combustion model is: " << CombustionModel << endl;
      cout << "Chemistry table has size " << n1 << " x " << n2 << " x " << n3
      << endl;
      cout << "x1 = Zmean:  " << "\t min=" << MinVal(x1, n1) << "  at i="
      << MinIndex(x1, n1) << "\t\t max=" << MaxVal(x1, n1) << "  at i="
      << MaxIndex(x1, n1) << endl;
      cout << "x2 = Zvar:  " << "\t min=" << MinVal(x2, n2) << "  at j="
      << MinIndex(x2, n2) << "\t\t max=" << MaxVal(x2, n2) << "  at j="
      << MaxIndex(x2, n2) << endl;
      strcpy(buffer_m, &VarNames[nvar - 1][0]);
      cout << "x3 = " << buffer_m << ":  " << "\t min=" << MinVal(x3, n3)
      << "  at k=" << MinIndex(x3, n3) << "\t\t max=" << MaxVal(x3, n3)
      << "  at k=" << MaxIndex(x3, n3) << endl;
      cout << "and contains " << nvar << " variables: " << endl;
      for (l = 0; l < nvar; l++)
      {
        strcpy(buffer_m, &VarNames[l][0]);
        GetArrayInfo(Data[l], n3, n2, n1, minval, maxval, mi3, mi2, mi1,
            ma3, ma2, ma1);
        cout.width(3);
        cout << l << ". ";
        cout.width(10);
        cout << buffer_m;
        cout << "\t min=" << minval << "\t at i=" << mi1 << "  j=" << mi2
        << "  k=" << mi3;
        cout << "\t max=" << maxval << "\t at i=" << ma1 << "  j=" << ma2
        << "  k=" << ma3 << endl;
      }
    }
    return;
  }

  /***********************************************************************************************************/

  /*! \brief Clean table and deallocate memory */
  void Unload()
  {
    if (x1 != NULL)       {freeMem1D(x1, 0, n1 - 1);                                         x1 = NULL;}
    if (x2 != NULL)       {freeMem1D(x2, 0, n2 - 1);                                         x2 = NULL;}
    if (x3 != NULL)       {freeMem1D(x3, 0, n3 - 1);                                         x3 = NULL;}
    if (iMask != NULL)    {freeMem3D(iMask, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1);                iMask = NULL;}
    if (VarNames != NULL) {delete[] VarNames;                                                VarNames = NULL;}
    if (Data != NULL)     {freeMem4D(Data, 0, nvar - 1, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1);    Data = NULL;}

    return;
  }

  /***********************************************************************************************************/

  /*! \brief Compute index and interpolation weights based on 3 input parameters.
   *
   *  Computes lower indices of cell containing the 3 input parameters (usually: Zmean, Zvar and chi/progress variable)
   *  and interpolation weights to compute linear interpolation. If coordinates (input) are out of table bounds,
   *  the coordinates are projected to boundary of table: if x < x(0), then x=x(0).
   *  Data is saved in variable #interp.
   *  \param[in]  A1  Value along first dimension (Zmean).
   *  \param[in]  A2  Value along second dimension (Zvar).
   *  \param[in]  A3  Value along third dimension (chi/progress variable).
   */
  void SetCoordinates(double A1, double A2, double A3)
  {
    int k, j, i;

    /* Determine index corresponding to direct lower value than input parameter */
    if (A1 <= x1[0])
    {
      interp.ii = 0;
      interp.wi = 1.0;
    }
    else
    {
      if (A1 >= x1[n1 - 1])
      {
        interp.ii = n1 - 2;
        interp.wi = 0.0;
      }
      else
      {
        for (i = 0; i < n1 - 1; i++)
        {
          if (A1 < x1[i + 1])
          {
            interp.ii = i;
            break;
          }
        }
        interp.wi = (x1[interp.ii + 1] - A1) / (x1[interp.ii + 1] - x1[interp.ii]);
      }
    }

    if (A2 <= x2[0])
    {
      interp.jj = 0;
      interp.wj = 1.0;
    }
    else
    {
      if (A2 >= x2[n2 - 1])
      {
        interp.jj = n2 - 2;
        interp.wj = 0.0;
      }
      else
      {
        for (j = 0; j < n2 - 1; j++)
        {
          if (A2 < x2[j + 1])
          {
            interp.jj = j;
            break;
          }
        }
        interp.wj = (x2[interp.jj + 1] - A2) / (x2[interp.jj + 1] - x2[interp.jj]);
      }
    }

    if (A3 <= x3[0])
    {
      interp.kk = 0;
      interp.wk = 1.0;
    }
    else
    {
      if (A3 >= x3[n3 - 1])
      {
        interp.kk = n3 - 2;
        interp.wk = 0.0;
      }
      else
      {
        for (k = 0; k < n3 - 1; k++)
        {
          if (A3 < x3[k + 1])
          {
            interp.kk = k;
            break;
          }
        }
        interp.wk = (x3[interp.kk + 1] - A3) / (x3[interp.kk + 1] - x3[interp.kk]);
      }
    }

    return;
  }

  /***********************************************************************************************************/

  /*! \brief Interpolate chemistry table based on indices, weights and name of variable.
   *
   *  Interpolation is linear. Any variable contained in table can be looked up. Indices and weights must be
   *  first computed using Chemtable::SetCoordinates . Returns the interpolated value of the variable.
   *  \param[in] tag          Name of variable to lookup.
   *  \return    Interpolated value of variable \p tag.
   */
  double Lookup(string tag)
  {
    int ll;
    double val;
    InterpolationIndex *p;

    /* Determine the index of the variable to lookup */
    itTp = myTableVar.find(tag);
    if (itTp == myTableVar.end())
    {
      cerr << "### Unknown variable in chemistry table (Lookup): " << tag << " ###" << endl;
      throw(-1);
    }
    ll = itTp->second;

    /* Interpolate value based on index and weights for given variable */
    p = &interp;
    val =          p->wk  * (          p->wj  * (          p->wi  * Data[ll][p->kk][p->jj][p->ii] 
                                                  + (1.0 - p->wi) * Data[ll][p->kk][p->jj][p->ii + 1]) 
                              + (1.0 - p->wj) * (          p->wi  * Data[ll][p->kk][p->jj + 1][p->ii] 
                                                  + (1.0 - p->wi) * Data[ll][p->kk][p->jj + 1][p->ii + 1]))
          + (1.0 - p->wk) * (          p->wj  * (          p->wi  * Data[ll][p->kk + 1][p->jj][p->ii] 
                                                  + (1.0 - p->wi) * Data[ll][p->kk + 1][p->jj][p->ii + 1]) 
                              + (1.0 - p->wj) * (          p->wi  * Data[ll][p->kk + 1][p->jj + 1][p->ii] 
                                                  + (1.0 - p->wi) * Data[ll][p->kk + 1][p->jj + 1][p->ii + 1]));
    return val;
  }

  /***********************************************************************************************************/

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
