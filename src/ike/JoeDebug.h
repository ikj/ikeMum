/*
 * JoeDebug.h
 *
 *  Created on: Jun 13, 2014
 *      Author: ikj
 */

#ifndef JOEDEBUG_H_
#define JOEDEBUG_H_

#include <iomanip> // std::setprecision

// ###########################################################################################
// ------                                    JoeDebug                                   ------
// ------                 This class is developed to help killing bugs                  ------
// ###########################################################################################
class JoeWithUtils : virtual public UgpWithCvCompFlow
{
protected:

public:
	// constructor
	JoeWithUtils() {
		if (mpi_rank == 0)
			cout << "JoeWithUtils()" << endl;

		init();
	}

	// Destructor
	~JoeWithUtils() {
		if (mpi_rank == 0)
			cout << "~JoeWithUtils()" << endl;

		clear();
	}

	void init() {
		MPI_Barrier(mpi_comm);
	}

	void clear() {
	}


	int calcNumNbrFAs(const int icv) {
		int nofa_f    = faocv_i[icv];
		int nofa_l_p1 = faocv_i[icv+1];
//		for (int nofa = nofa_f; nofa < nofa_l_p1; nofa++) {
//			int ifa = faocv_v[nofa];
//		}

		return max(1, nofa_l_p1 - nofa_f);
	}
};

// ###########################################################################################
// ------                                    JoeDebug                                   ------
// ------                 This class is developed to help killing bugs                  ------
// ###########################################################################################
class JoeDebug : public JoeWithUtils
{
protected:

public:
	// constructor
	JoeDebug() {
		if (mpi_rank == 0)
			cout << "JoeDebug()" << endl;

		init();
	}

	// Destructor
	~JoeDebug() {
		if (mpi_rank == 0)
			cout << "~JoeDebug()" << endl;

		clear();
	}

	void init() {
		MPI_Barrier(mpi_comm);
	}

	void clear() {
	}

	/*
	 * Method: showAddressRegisteredData
	 * ---------------------------------
	 * Show the list of registered data for all the CPU cores
	 */
	void showAddressRegisteredData() {
		MPI_Barrier(mpi_comm);
		if(mpi_rank != 0) { MPI_Status status; int dummy; MPI_Recv(&dummy, 1, MPI_INT, mpi_rank-1, 412, mpi_comm, &status); }

		if(mpi_rank==0) cout<<endl<<"JoeDebug::showAddressRegisteredData() "<<endl;
		cout<<">> mpi_rank="<<mpi_rank<<endl;
		myShowAddressRegisteredData();
		cout<<endl;
		if(mpi_rank==mpi_size-1) cout<<endl;

		if(mpi_rank < mpi_size-1) { int dummy; MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 412, mpi_comm); }
		MPI_Barrier(mpi_comm);
	}
	/*
	 * Method: myShowAddressRegisteredData
	 * -----------------------------------
	 * Show the list of registered data on a single CPU core
	 * Details of the DATA type: Ugp.h
	 */
	void myShowAddressRegisteredData() {
		cout<<"  > INT_VALUE_LIST:  ==========="<<endl;
		int nameLength = 0;
		for (list<IntValue>::iterator data = intValueList.begin(); data!=intValueList.end(); data++)
			nameLength = max(nameLength, (int)strlen(data->getName()));
		for (list<IntValue>::iterator data = intValueList.begin(); data!=intValueList.end(); data++) {
			cout<<"    NAME="<<std::setw(nameLength)<<data->getName()<<": ADDRESS="<<data->ptr<<" (DATA_TYPE=";
			switch(data->getDatatype()) {
			case VALUE_DATA: cout<<"VALUE_DATA"; break;
			case CV_DATA: cout<<"CV_DATA"; break;
			case FA_DATA: cout<<"FA_DATA"; break;
			case NO_DATA: cout<<"NO_DATA"; break;
			default: cout<<"UNKNOWN"; break;
			}
			cout<<", FLAG="<<data->getFlag()<<")"<<endl;
		}

		cout<<"  > INT_SCALAR_LIST:  ==========="<<endl;
		nameLength = 0;
		for (list<IntScalar>::iterator data = intScalarList.begin(); data!=intScalarList.end(); data++)
			nameLength = max(nameLength, (int)strlen(data->getName()));
		for (list<IntScalar>::iterator data = intScalarList.begin(); data!=intScalarList.end(); data++) {
			cout<<"    NAME="<<std::setw(nameLength)<<data->getName()<<": ADDRESS="<<*(data->ptr)<<" (DATA_TYPE=";
			switch(data->getDatatype()) {
			case VALUE_DATA: cout<<"VALUE_DATA"; break;
			case CV_DATA: cout<<"CV_DATA"; break;
			case FA_DATA: cout<<"FA_DATA"; break;
			case NO_DATA: cout<<"NO_DATA"; break;
			default: cout<<"UNKNOWN"; break;
			}
			cout<<", FLAG="<<data->getFlag()<<")"<<endl;
		}

		cout<<"  > DOUBLE_VALUE_LIST:  ==========="<<endl;
		nameLength = 0;
		for (list<DoubleValue>::iterator data = doubleValueList.begin(); data!=doubleValueList.end(); data++)
			nameLength = max(nameLength, (int)strlen(data->getName()));
		for (list<DoubleValue>::iterator data = doubleValueList.begin(); data!=doubleValueList.end(); data++) {
			cout<<"    NAME="<<std::setw(nameLength)<<data->getName()<<": ADDRESS="<<data->ptr<<" (DATA_TYPE=";
			switch(data->getDatatype()) {
			case VALUE_DATA: cout<<"VALUE_DATA"; break;
			case CV_DATA: cout<<"CV_DATA"; break;
			case FA_DATA: cout<<"FA_DATA"; break;
			case NO_DATA: cout<<"NO_DATA"; break;
			default: cout<<"UNKNOWN"; break;
			}
			cout<<", FLAG="<<data->getFlag()<<")"<<endl;
		}

		cout<<"  > DOUBLE_SCALAR_LIST:  ==========="<<endl;
		nameLength = 0;
		for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++)
			nameLength = max(nameLength, (int)strlen(data->getName()));
		for (list<DoubleScalar>::iterator data = doubleScalarList.begin(); data!=doubleScalarList.end(); data++) {
			cout<<"    NAME="<<std::setw(nameLength)<<data->getName()<<": ADDRESS="<<*(data->ptr)<<" (DATA_TYPE=";
			switch(data->getDatatype()) {
			case VALUE_DATA: cout<<"VALUE_DATA"; break;
			case CV_DATA: cout<<"CV_DATA"; break;
			case FA_DATA: cout<<"FA_DATA"; break;
			case NO_DATA: cout<<"NO_DATA"; break;
			default: cout<<"UNKNOWN"; break;
			}
			cout<<", FLAG="<<data->getFlag()<<")"<<endl;
		}

		cout<<"  > DOUBLE_VECTOR_LIST:  ==========="<<endl;
		nameLength = 0;
		for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++)
			nameLength = max(nameLength, (int)strlen(data->getName()));
		for (list<DoubleVector>::iterator data = doubleVectorList.begin(); data!=doubleVectorList.end(); data++) {
			cout<<"    NAME="<<std::setw(nameLength)<<data->getName()<<": ADDRESS="<<*(data->ptr)<<" (DATA_TYPE=";
			switch(data->getDatatype()) {
			case VALUE_DATA: cout<<"VALUE_DATA"; break;
			case CV_DATA: cout<<"CV_DATA"; break;
			case FA_DATA: cout<<"FA_DATA"; break;
			case NO_DATA: cout<<"NO_DATA"; break;
			default: cout<<"UNKNOWN"; break;
			}
			cout<<", FLAG="<<data->getFlag()<<")"<<endl;
		}

		cout<<"  > DOUBLE_TENSOR_LIST:  ==========="<<endl;
		nameLength = 0;
		for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data!=doubleTensorList.end(); data++)
			nameLength = max(nameLength, (int)strlen(data->getName()));
		for (list<DoubleTensor>::iterator data = doubleTensorList.begin(); data!=doubleTensorList.end(); data++) {
			cout<<"    NAME="<<std::setw(nameLength)<<data->getName()<<": ADDRESS="<<*(data->ptr)<<" (DATA_TYPE=";
			switch(data->getDatatype()) {
			case VALUE_DATA: cout<<"VALUE_DATA"; break;
			case CV_DATA: cout<<"CV_DATA"; break;
			case FA_DATA: cout<<"FA_DATA"; break;
			case NO_DATA: cout<<"NO_DATA"; break;
			default: cout<<"UNKNOWN"; break;
			}
			cout<<", FLAG="<<data->getFlag()<<")"<<endl;
		}
	}


	/*
	 * Method: showAddressUgpWithCvCompFlowArray
	 * -----------------------------------------
	 * Show the list of the member array in the UgpWithCvCompFlow class
	 */
	void showAddressUgpWithCvCompFlowArray() {
		MPI_Barrier(mpi_comm);
		if(mpi_rank != 0) { MPI_Status status; int dummy; MPI_Recv(&dummy, 1, MPI_INT, mpi_rank-1, 413, mpi_comm, &status); }

		if(mpi_rank==0) cout<<endl<<"JoeDebug::showAddressUgpWithCvCompFlowArray() "<<endl;
		cout<<">> mpi_rank="<<mpi_rank<<endl;
		myShowAddressUgpWithCvCompFlowArray();
		cout<<endl;
		if(mpi_rank==mpi_size-1) cout<<endl;

		if(mpi_rank < mpi_size-1) { int dummy; MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 413, mpi_comm); }
		MPI_Barrier(mpi_comm);
	}
	/*
	 * Method: myShowAddressUgpWithCvCompFlowArray
	 * -------------------------------------------
	 * Show the list of the member array in the UgpWithCvCompFlow class
	 * Details of the DATA type: UgpWithCvCompFlow.h
	 */
	void myShowAddressUgpWithCvCompFlowArray() {
		cout<<"  > SCALAR ARRAY:  ==========="<<endl;
		cout<<"    NAME=mul_fa    : ADDRESS="<<mul_fa<<endl;
		cout<<"    NAME=lamOcp_fa : ADDRESS="<<lamOcp_fa<<endl;
		cout<<"    NAME=mut_fa    : ADDRESS="<<mut_fa<<endl;
		cout<<"    NAME=kine      : ADDRESS="<<kine<<endl;
		cout<<"    NAME=strMag    : ADDRESS="<<strMag<<endl;
		cout<<"    NAME=vortMag   : ADDRESS="<<vortMag<<endl;
		cout<<"    NAME=diverg    : ADDRESS="<<diverg<<endl;

		cout<<"  > VECTOR ARRAY:  ==========="<<endl;
		cout<<"    NAME=grad_rho      : ADDRESS="<<grad_rho<<endl;
		cout<<"    NAME=grad_p        : ADDRESS="<<grad_p<<endl;
		cout<<"    NAME=grad_temp     : ADDRESS="<<grad_temp<<endl;
		cout<<"    NAME=grad_enthalpy : ADDRESS="<<grad_enthalpy<<endl;

		cout<<"  > TENSOR ARRAY:  ==========="<<endl;
		cout<<"    NAME=grad_u : ADDRESS="<<grad_u<<endl;
	}


	/*
	 * Method: showAddressScalarTransport
	 * ----------------------------------
	 * Show the list of scalar transport equation
	 */
	void showAddressScalarTransport() {
		MPI_Barrier(mpi_comm);
		if(mpi_rank != 0) { MPI_Status status; int dummy; MPI_Recv(&dummy, 1, MPI_INT, mpi_rank-1, 414, mpi_comm, &status); }

		if(mpi_rank==0) cout<<endl<<"JoeDebug::showAddressScalarTransport() "<<endl;
		cout<<">> mpi_rank="<<mpi_rank<<endl;
		myShowAddressScalarTransport();
		cout<<endl;
		if(mpi_rank==mpi_size-1) cout<<endl;

		if(mpi_rank < mpi_size-1) { int dummy; MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 414, mpi_comm); }
		MPI_Barrier(mpi_comm);
	}
	/*
	 * Method: myShowAddressScalarTransport
	 * ------------------------------------
	 * Show the list of scalar transport equation
	 * Details of the DATA type: UgpWithCv2.h
	 */
	void myShowAddressScalarTransport() {
		for (int i = 0; i < scalarTranspEqVector.size(); i++) {
			ScalarTranspEq *transScal = &(scalarTranspEqVector[i]);
			cout<<"  > SCALAR NAME = "<<transScal->getName()<<" ==========="<<endl;
			cout<<"    DATA_NAME=phi         : ADDRESS="<<transScal->phi        <<" (SCALAR ARRAY)"<<endl;
			cout<<"    DATA_NAME=diff        : ADDRESS="<<transScal->diff       <<" (SCALAR ARRAY)"<<endl;
			cout<<"    DATA_NAME=rhophi      : ADDRESS="<<transScal->rhophi     <<" (SCALAR ARRAY)"<<endl;
			cout<<"    DATA_NAME=dpress_dphi : ADDRESS="<<transScal->dpress_dphi<<" (SCALAR ARRAY)"<<endl;
			cout<<"    DATA_NAME=grad_phi    : ADDRESS="<<transScal->grad_phi   <<" (VECTOR ARRAY)"<<endl;
			cout<<"    DATA_NAME=grad_rhophi : ADDRESS="<<transScal->grad_rhophi<<" (VECTOR ARRAY)"<<endl;
		}
	}


	/*
	 * Method: calcResidualsFrom1Drhs
	 * ------------------------------
	 * Original code = IkeWithPsALC::calcResidualsFrom1Drhs()
	 *
	 * Calculate residual
	 * Supported norm: inf-norm, one-norm, two-norm
	 *                   whichNorm==0 -> infinity-norm
	 *                   whichNorm==1 -> one-norm
	 *                   whichNorm==2 -> two-norm
	 * Note: "Residual" should have been defined as " double *Residual = new double[5+nScal]; "
	 */
	void calcResidualsFrom1Drhs(double *Residual, double *rhs1Darray, const int whichNorm) {
		int nScal = scalarTranspEqVector.size();
		int nVars = 5+nScal;
		double *myResidual = new double[nVars];

		switch (whichNorm) {
		case 0: // note: infinity-norm
			for (int i = 0; i < 5+nScal; i++) {
				myResidual[i] = -ABSURDLY_BIG_NUMBER;
				Residual[i]   = -ABSURDLY_BIG_NUMBER;
			}

			for (int icv = 0; icv < ncv; icv++) {
				if(x_cv[icv][0] > 0.0) { // Note: Consider only the combustor
					int icv_i = nVars*icv;

					for(int i=0; i<5+nScal; ++i)
						if(fabs(rhs1Darray[icv_i+i]) > myResidual[i])
							myResidual[i] = fabs(rhs1Darray[icv_i+i]);
				}
			}
			MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_MAX, mpi_comm);

			break;
		case 1:
			for (int i = 0; i < 5+nScal; i++) {
				myResidual[i] = 0.0;
				Residual[i]   = 0.0;
			}

			for (int icv = 0; icv < ncv; icv++) {
				if(x_cv[icv][0] > 0.0) { // Note: Consider only the combustor
					int icv_i = nVars*icv;

					for(int i=0; i<5+nScal; ++i)
						myResidual[i] += fabs(rhs1Darray[icv_i+i]);
				}
			}
			MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);

			break;
		case 2:
			for (int i = 0; i < 5+nScal; i++) {
				myResidual[i] = 0.0;
				Residual[i]   = 0.0;
			}

			for (int icv = 0; icv < ncv; icv++) {
				if(x_cv[icv][0] > 0.0) { // Note: Consider only the combustor
					int icv_i = nVars*icv;

					for(int i=0; i<5+nScal; ++i)
						myResidual[i] += fabs(rhs1Darray[icv_i+i])*fabs(rhs1Darray[icv_i+i]);
				}
			}
			MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);
			for (int i = 0; i < 5+nScal; i++) {
				Residual[i] = sqrt(Residual[i]);
			}

			break;
		default:
			if(mpi_rank==0)
				cout<<"ERROR in IkeWithPsALC_AD::calcResidualsFrom1Drhs(): unsupported norm = "<<whichNorm<<endl;
			throw(PSALC_ERROR_CODE);
			break;
		}

		delete [] myResidual;
	}
};


#endif /* JOEDEBUG_H_ */
