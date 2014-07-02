/*
 * ADvar.h
 *
 *  Created on: May 22, 2013
 *      Author: ikj
 */

#ifndef ADVAR_H_
#define ADVAR_H_

#include <iostream>
#include <vector>
#include <string>
#include <assert.h>
#include "adolc.h"

#define ADVAR_ERROR_CODE 2147483647 // Note: Maximum number for (signed) long int = 2147483647
#define ADVAR_DEBUG_LEVEL 1

//#define USE_PROXY_OBJECT_

/*
 * Class: ADscalar
 * ---------------
 * A class to store some AD scalar variables (e.g. rho) for the "1D-style" memory-saving code.
 * The algorithm of this class came from the Stanford CS-106L course reader chapter 10:
 * http://www.stanford.edu/class/cs106l/course-reader/Ch10_OperatorOverloading.pdf
 */
template <typename T>
class ADscalar {
private:
	bool isMemSavingMode;
	string varName;
	T* arrayAD;

	// For the memory saving mode (the "1D-style")
	int myIcv;
	vector<int> myNbocv2;

	// For the regular mode (normal adjoint)
	int myNcv; // e.g. ncv_ggff, ncv_gg, etc.

public:
	/*
	 * Constructor and Destructor
	 */
	ADscalar(char variableName[]="NOT_DEFINED"){
		varName = variableName;
		arrayAD = NULL;

		myIcv = ADVAR_ERROR_CODE;
		myNcv = 0;
	}
	~ADscalar() {
		clear();
	}

	/*
	 * Copy Constructors and Assignment Operators
	 */
	ADscalar(const ADscalar& other) {
		copyOther(other);
	}
	ADscalar<T>& operator = (const ADscalar& other) {
		if(this != &other) {
			clear();
			// Note: When we cover inheritance, there's one more step here.
			copyOther(other);
		}
		return *this;
	}

	/*
	 * Method: copyOther
	 * -----------------
	 *
	 */
	void copyOther(const ADscalar& other) {
		isMemSavingMode = other.isMemSavingMode;
		varName = other.varName;

		myIcv = other.myIcv;
		myNbocv2 = other.myNbocv2;

		myNcv = other.myNcv;

		arrayAD = new T[other.size()];
		copy(other.begin(), other.end(), begin());
	}

	/*
	 * Method: allocate
	 * ----------------
	 * There are two versions: 1. For the memory saving mode (the "1D-style")
	 *                         2. For the normal mode (normal adjoint)
	 */
	void allocate(vector<int>& localNbocv2, const int icvCenter = ADVAR_ERROR_CODE) {
		isMemSavingMode = true;

		myNbocv2 = localNbocv2;
		myIcv = icvCenter;

		assert(arrayAD == NULL);
		arrayAD = new T [localNbocv2.size()];
	}
	void allocate(const int ncv) {
		isMemSavingMode = false;

		myNcv = ncv;

		assert(arrayAD == NULL);
		arrayAD = new T [myNcv];
	}

	/*
	 * Method: clear
	 * -------------
	 *
	 */
	void clear() {
		if(!empty()) {
//			varName.clear();
			if(arrayAD!=NULL) {
				delete [] arrayAD;
				arrayAD=NULL;
			}

			myNbocv2.clear();
			myNcv = 0;
		}

		myIcv = ADVAR_ERROR_CODE;
	}

	/*
	 * Method: empty
	 * -------------
	 *
	 */
	bool empty() const {
		if(arrayAD==NULL) {
			if(isMemSavingMode) {
				if(!myNbocv2.empty())
					cout<<"WARNING! ADscalar::empty(): the array is empty, but (vector<int> myNbocv2) is not empty"<<endl;
			} else {
				if(myNcv!=0)
					cout<<"WARNING! ADscalar::empty(): the array is empty, but (myNcv) is not zero"<<endl;
			}

			return true;
		}
		return false;
	}

	/*
	 * Method: getAt
	 * -------------
	 *
	 */
	T& getAt(const int icv) {
		if(isMemSavingMode) {
			size_t index = convertIcvToIndex(icv);
			return arrayAD[index];
		}

		assert(icv<myNcv);
		return arrayAD[icv];
	}
	const T& getAt(const int icv) const {
		if(isMemSavingMode) {
			size_t index = convertIcvToIndex(icv);
			return arrayAD[index];
		}

		assert(icv<myNcv);
		return arrayAD[icv];
	}

	/*
	 * Method: getMyIcv
	 * ----------------
	 *
	 */
	int getMyIcv() const {
		if(myIcv == ADVAR_ERROR_CODE)
			cout<<"WARNING in ADscalar::getMyIcv(): myIcv has not been defined"<<endl;
		return myIcv;
	}

	/*
	 * Method: size()
	 * --------------
	 *
	 */
	int size() const {
		if(isMemSavingMode)
			return myNbocv2.size();

		return myNcv;
	}

	/*
	 * Method: getName
	 * ---------------
	 *
	 */
	string& getName() {
		return varName;
	}
	const string& getName() const {
		return varName;
	}

	/*
	 * Method: setName
	 * ---------------
	 *
	 */
	void setName(char variableName[]){
		varName.clear();
		varName = variableName;
	}

	/*
	 * Method: getMyNbocv2
	 * -------------------
	 *
	 */
	vector<int>& getMyNbocv2() {
		return myNbocv2;
	}
	const vector<int>& getMyNbocv2() const {
		return myNbocv2;
	}

	/***********************************************************
	 * Element Selection Operator: []
	 * ------------------------------
	 * Two overloaded operators for the mutable and immutable cases
	 ***********************************************************/
	T& operator[] (const int icv) {
		return getAt(icv);
	}

	const T& operator[] (const int icv) const {
		return getAt(icv);
	}

	/***********
	 * Iterators
	 * ---------
	 * Sample code from http://www.dreamincode.net/forums/topic/58468-making-your-own-iterators/
	 * For example, if you want to iterate over the elements of "ADscalar<elemType> arr" and print them out:
	 *   for(ADscalar<elemType>::iterator iter=arr.begin(); iter!=arr.end(); ++iter) {
	 *   	printf("  %.2e", *iter); // Note: for the "adouble" elemType: printf("  %.2e", iter->value());
	 *   }
	 ***********/
	class iterator
	{
	public:
		iterator(T* ptr) : ptr_(ptr) { }
		iterator operator ++() { iterator i = *this; ptr_++; return i; }
		iterator operator ++(int junk) { ptr_++; return *this; }
		T& operator *() { return *ptr_; }
		T* operator ->() { return ptr_; }
		bool operator ==(const iterator& rhs) { return ptr_ == rhs.ptr_; }
		bool operator !=(const iterator& rhs) { return ptr_ != rhs.ptr_; }
	private:
		T* ptr_;
	};

	class const_iterator
	{
	public:
		const_iterator(T* ptr) : ptr_(ptr) { }
		const_iterator operator ++() { const_iterator i = *this; ptr_++; return i; }
		const_iterator operator ++(int junk) { ptr_++; return *this; }
		const T& operator *() { return *ptr_; }
		const T* operator ->() { return ptr_; }
		bool operator ==(const const_iterator& rhs) { return ptr_ == rhs.ptr_; }
		bool operator !=(const const_iterator& rhs) { return ptr_ != rhs.ptr_; }
	private:
		T* ptr_;
	};

	friend class iterator;
	friend class const_iterator;

	iterator begin() {
		return iterator(arrayAD);
	}

	iterator end() {
		return iterator(arrayAD + this->size());
	}

	const_iterator begin() const {
		return const_iterator(arrayAD);
	}

	const_iterator end() const {
		return const_iterator(arrayAD + this->size());
	}

	/*
	 * Method: copy
	 * ------------
	 *
	 */
	iterator copy(const_iterator first, const_iterator last, iterator result) {
		while(first!=last) {
			*result = *first;
			++result;
			++first;
		}
		return result;
	}

private:
	/*
	 * Method: convertIcvToIndex
	 * -------------------------
	 *
	 */
	size_t convertIcvToIndex(const int icv) const {
		size_t index = ADVAR_ERROR_CODE;
		bool found = false;

		for(size_t i=0; i<myNbocv2.size(); ++i) {
			if(myNbocv2[i] == icv) {
				index = i;
				found = true;
				break;
			}
		}

		if(!found) {
			cerr<<"ERROR in ADscalar::convertIcvToIndex(): Cannot find the given icv = "<<icv<<", neighgor of icvCenter="<<getMyIcv()<<", varName="<<getName()<<endl;
#if ADVAR_DEBUG_LEVEL>0
			cerr<<"                                        Members of myNbocv2 = ";
			for(size_t i=0; i<myNbocv2.size(); ++i)
				cerr<<myNbocv2[i]<<"  ";
			cerr<<endl;
#endif
			assert(false);
		}

		return index;
	}
};

/*
 * Class: ADvector
 * ---------------
 * A class to store some AD vector variables (e.g. rhou) for the "1D-style" memory-saving code.
 *
 * You can use something like rhou[10][1]
 *
 * The algorithm of this class came from Stanford CS-106L course reader chapter 10:
 * http://www.stanford.edu/class/cs106l/course-reader/Ch10_OperatorOverloading.pdf
 */
template <typename T>
class ADvector {
private:
	bool isMemSavingMode;
	string varName;
	T (*arrayAD)[3];

	// For the memory saving mode (the "1D-style")
	int myIcv;
	vector<int> myNbocv2;

	// For the regular mode (normal adjoint)
	int myNcv; // e.g. ncv_ggff, ncv_gg, etc.

public:
	/*
	 * Constructor and Destructor
	 */
	ADvector(char variableName[]="NOT_DEFINED"){
		varName = variableName;
		arrayAD = NULL;

		myIcv = ADVAR_ERROR_CODE;
		myNcv = 0;
	}
	~ADvector() {
		clear();
	}

	/*
	 * Copy Constructors and Assignment Operators
	 */
	ADvector(const ADvector& other) {
		copyOther(other);
	}
	ADvector<T>& operator = (const ADvector& other) {
		if(this != &other) {
			clear();
			// Note: When we cover inheritance, there's one more step here.
			copyOther(other);
		}
		return *this;
	}

	/*
	 * Method: copyOther
	 * -----------------
	 *
	 */
	void copyOther(const ADvector& other) {
		isMemSavingMode = other.isMemSavingMode;
		varName = other.varName;

		myIcv = other.myIcv;
		myNbocv2 = other.myNbocv2;

		myNcv = other.myNcv;

		arrayAD = new T[other.size()][3];
		copy(other.begin(), other.end(), begin());
	}

	/*
	 * Method: allocate
	 * ----------------
	 * There are two versions: 1. For the memory saving mode (the "1D-style")
	 *                         2. For the normal mode (normal adjoint)
	 */
	void allocate(vector<int>& localNbocv2, const int icvCenter = ADVAR_ERROR_CODE) {
		isMemSavingMode = true;

		myNbocv2 = localNbocv2;
		myIcv = icvCenter;

		assert(arrayAD == NULL);
		arrayAD = new T [localNbocv2.size()][3];
	}
	void allocate(const int ncv) {
		isMemSavingMode = false;

		myNcv = ncv;

		assert(arrayAD == NULL);
		arrayAD = new T [myNcv][3];
	}

	/*
	 * Method: clear
	 * -------------
	 *
	 */
	void clear() {
		if(!empty()) {
//			varName.clear();
			if(arrayAD!=NULL) {
				delete [] arrayAD;
				arrayAD=NULL;
			}

			myNbocv2.clear();
			myNcv = 0;
		}

		myIcv = ADVAR_ERROR_CODE;
	}

	/*
	 * Method: empty
	 * -------------
	 *
	 */
	bool empty() const {
		if(arrayAD==NULL) {
			if(isMemSavingMode) {
				if(!myNbocv2.empty())
					cout<<"WARNING! ADvector::empty(): the array is empty, but (vector<int> myNbocv2) is not empty"<<endl;
			} else {
				if(myNcv!=0)
					cout<<"WARNING! ADvector::empty(): the array is empty, but (myNcv) is not zero"<<endl;
			}

			return true;
		}
		return false;
	}

	/*
	 * Method: getAt
	 * -------------
	 *
	 */
	T& getAt(const int icv, const int dim) {
		if(isMemSavingMode) {
			size_t index = convertIcvToIndex(icv);
			assert(dim>=0 && dim<3);
			return arrayAD[index][dim];
		}

		assert(icv<myNcv && dim>=0 && dim<3);
		return arrayAD[icv][dim];
	}
	const T& getAt(const int icv, const int dim) const {
		if(isMemSavingMode) {
			size_t index = convertIcvToIndex(icv);
			assert(dim>=0 && dim<3);
			return arrayAD[index][dim];
		}

		assert(icv<myNcv && dim>=0 && dim<3);
		return arrayAD[icv][dim];
	}
//	T* getAt(const int icv) {
//		if(isMemSavingMode) {
//			size_t index = convertIcvToIndex(icv);
//			return arrayAD[index];
//		}
//
//		assert(icv<myNcv);
//		return arrayAD[icv];
//	}
//	const T* getAt(const int icv) const {
//		if(isMemSavingMode) {
//			size_t index = convertIcvToIndex(icv);
//			return arrayAD[index];
//		}
//
//		assert(icv<myNcv);
//		return arrayAD[icv];
//	}
	T (&getAt(const int icv))[3] {
		if(isMemSavingMode) {
			size_t index = convertIcvToIndex(icv);
			return arrayAD[index];
		}

		assert(icv<myNcv);
		return arrayAD[icv];
	}
	const T (&getAt(const int icv) const)[3] {
		if(isMemSavingMode) {
			size_t index = convertIcvToIndex(icv);
			return arrayAD[index];
		}

		assert(icv<myNcv);
		return arrayAD[icv];
	}

	/*
	 * Method: getMyIcv
	 * ----------------
	 *
	 */
	int getMyIcv() const {
		if(myIcv == ADVAR_ERROR_CODE)
			cout<<"WARNING in ADvector::getMyIcv: myIcv has not been defined"<<endl;
		return myIcv;
	}

	/*
	 * Method: size()
	 * --------------
	 *
	 */
	int size() const {
		if(isMemSavingMode)
			return myNbocv2.size();

		return myNcv;
	}

	/*
	 * Method: getName
	 * ---------------
	 *
	 */
	string& getName() {
		return varName;
	}
	const string& getName() const {
		return varName;
	}

	/*
	 * Method: setName
	 * ---------------
	 *
	 */
	void setName(char variableName[]){
		varName.clear();
		varName = variableName;
	}

	/*
	 * Method: getMyNbocv2
	 * -------------------
	 *
	 */
	vector<int>& getMyNbocv2() {
		return myNbocv2;
	}
	const vector<int>& getMyNbocv2() const {
		return myNbocv2;
	}

	/*****************************
	 * Element Selection Operator
	 *****************************/
	// Note: In general, we need the so-called "proxy object" for a multi-dimensional container.
	//       In the case of 3D CFD, we don't need to use it due to the special data structure
	//       (i.e., something like "double (*rhou)[3]").
	//       However, I leave the proxy-object implementation for the future reference.
#ifdef USE_PROXY_OBJECT_
	/*
	 * We need the so-called "proxy object" for this vector case:
	 * If we write
	 *   value = ADvector[10][2],
	 * The following sequences of actions occurs:
	 *   1. ADvector.operator[] is invoked with the parameter 10.
	 *   2. ADvector.operator[] creates a MutableReference storing the icv 10 and a means for
	 *     communicating back with the ADvector object.
	 *   3. ADvector.operator[] returns this MutableReference.
	 *   4. The returned MutableReference then has its operator[] function called with parameter 2.
	 *   5. The returned MutableReference then calls back to the ADvector object and asks for
	 *     the element at position (10, 2).
	 */
	/*
	 * Class: MutableReference
	 * -----------------------
	 * For the mutable case
	 */
	class MutableReference {
	public:
		friend class ADvector;
		friend class MutableProxy;
		T& operator[] (const int dim) {
			return owner->getAt(icv, dim);
		}

	private:
		MutableReference(ADvector* owner, int icv): owner(owner), icv(icv) {};
		ADvector* const owner;
		const int icv;
	};

	/*
	 * Operator: []
	 * ------------
	 * Element selection operator for the mutable case
	 */
	MutableReference operator [] (int icv) {
		return MutableReference(this, icv);
	}

	/*
	 * Class: ImmutableReference
	 * -------------------------
	 * For the immutable case
	 */
	class ImmutableReference {
	public:
		friend class ADvector;
		const T& operator[] (const int dim) const {
			return owner->getAt(icv, dim);
		}
	private:
		ImmutableReference(const ADvector* owner, int icv): owner(owner), icv(icv) {};
		const ADvector* const owner;
		const int icv;
	};

	/*
	 * Operator: []
	 * ------------
	 * Element selection operator for the immutable case
	 */
	ImmutableReference operator[] (int icv) const {
		return ImmutableReference(this, icv);
	}
#else
	/*
	 * Operator: []
	 * ------------
	 * Element selection operator for the mutable case
	 */
	T (&operator [] (int icv))[3] {
		return getAt(icv);
	}
	/*
	 * Operator: []
	 * ------------
	 * Element selection operator for the immutable case
	 */
	T (&operator[] (int icv) const)[3] {
		return getAt(icv);
	}
#endif

	/***********
	 * Iterators
	 * ---------
	 * For example, if you want to iterate over the elements of "ADvector<elemType> vec" and print them out:
	 *   for(ADvector<elemType>::iterator iter=vec.begin(); iter!=vec.end(); ++iter) {
	 *   	for(int i=0; i<3; ++i)
	 *   		printf("  %.2e", (*iter)[i]); // Note: for the "adouble" elemType: printf("  %.2e", (*iter)[i].value());
	 *   	printf("\n");
	 *   }
	 ***********/
	class iterator
	{
	public:
		iterator(T *ptr) : ptr_(ptr) { }
		iterator operator ++() { iterator i = *this; ptr_+=3; return i; }
		iterator operator ++(int junk) { ptr_+=3; return *this; }
		T* operator *() { return ptr_; }
		bool operator ==(const iterator& rhs) { return ptr_ == rhs.ptr_; }
		bool operator !=(const iterator& rhs) { return ptr_ != rhs.ptr_; }
	private:
		T *ptr_;
	};

	class const_iterator
	{
	public:
		const_iterator(T *ptr) : ptr_(ptr) { }
		const_iterator operator ++() { const_iterator i = *this; ptr_+=3; return i; }
		const_iterator operator ++(int junk) { ptr_+=3; return *this; }
		const T* operator *() const { return ptr_; }
		bool operator ==(const const_iterator& rhs) { return ptr_ == rhs.ptr_; }
		bool operator !=(const const_iterator& rhs) { return ptr_ != rhs.ptr_; }
	private:
		T *ptr_;
	};

	friend class iterator;
	friend class const_iterator;

	iterator begin() {
		return iterator(arrayAD[0]);
	}

	iterator end() {
		return iterator(arrayAD[0] + (this->size())*3);
	}

	const_iterator begin() const {
		return const_iterator(arrayAD[0]);
	}

	const_iterator end() const {
		return const_iterator(arrayAD[0] + (this->size())*3);
	}

	/*
	 * Method: copy
	 * ------------
	 *
	 */
	iterator copy(const_iterator first, const_iterator last, iterator result) {
		while(first!=last) {
			for(int i=0; i<3; ++i)
				(*result)[i] = (*first)[i];

			++result;
			++first;
		}
		return result;
	}

private:
	/*
	 * Method: convertIcvToIndex
	 * -------------------------
	 *
	 */
	size_t convertIcvToIndex(const int icv) const {
		size_t index = ADVAR_ERROR_CODE;
		bool found = false;

		for(size_t i=0; i<myNbocv2.size(); ++i) {
			if(myNbocv2[i] == icv) {
				index = i;
				found = true;
				break;
			}
		}

		if(!found) {
			cerr<<"ERROR in ADvector::convertIcvToIndex(): Cannot find the given icv = "<<icv<<", neighgor of icvCenter="<<getMyIcv()<<", varName="<<getName()<<endl;
#if ADVAR_DEBUG_LEVEL>0
			cerr<<"                                        Members of myNbocv2 = ";
			for(size_t i=0; i<myNbocv2.size(); ++i)
				cerr<<myNbocv2[i]<<"  ";
			cerr<<endl;
#endif
			assert(false);
		}

		return index;
	}
};

/*
 * Class: ADtensor
 * ---------------
 * A class to store a AD tensor variable (e.g. grad_u) for the "1D-style" memory-saving code.
 *
 * You can use something like grad_u[111][2][2]
 *
 * The algorithm of this class came from Stanford CS-106L course reader chapter 10:
 * http://www.stanford.edu/class/cs106l/course-reader/Ch10_OperatorOverloading.pdf
 */
template <typename T>
class ADtensor {
private:
	bool isMemSavingMode;
	string varName;
	T (*arrayAD)[3][3];

	// For the memory saving mode (the "1D-style")
	int myIcv;
	vector<int> myNbocv2;

	// For the regular mode (normal adjoint)
	int myNcv; // e.g. ncv_ggff, ncv_gg, etc.

public:
	/*
	 * Constructor and Destructor
	 */
	ADtensor(char variableName[]="NOT_DEFINED"){
		varName = variableName;
		arrayAD = NULL;

		myIcv = ADVAR_ERROR_CODE;
		myNcv = 0;
	}
	~ADtensor() {
		clear();
	}

	/*
	 * Copy Constructors and Assignment Operators
	 */
	ADtensor(const ADtensor& other) {
		copyOther(other);
	}
	ADtensor<T>& operator = (const ADtensor& other) {
		if(this != &other) {
			clear();
			// Note: When we cover inheritance, there's one more step here.
			copyOther(other);
		}
		return *this;
	}

	/*
	 * Method: copyOther
	 * -----------------
	 *
	 */
	void copyOther(const ADtensor& other) {
		isMemSavingMode = other.isMemSavingMode;
		varName = other.varName;

		myIcv = other.myIcv;
		myNbocv2 = other.myNbocv2;

		myNcv = other.myNcv;

		arrayAD = new T[other.size()][3][3];
		copy(other.begin(), other.end(), begin());
	}

	/*
	 * Method: allocate
	 * ----------------
	 * There are two versions: 1. For the memory saving mode (the "1D-style")
	 *                         2. For the normal mode (normal adjoint)
	 */
	void allocate(vector<int>& localNbocv2, const int icvCenter = ADVAR_ERROR_CODE) {
		isMemSavingMode = true;

		myNbocv2 = localNbocv2;
		myIcv = icvCenter;

		assert(arrayAD == NULL);
		arrayAD = new T [localNbocv2.size()][3][3];
	}
	void allocate(const int ncv) {
		isMemSavingMode = false;

		myNcv = ncv;

		assert(arrayAD == NULL);
		arrayAD = new T [myNcv][3][3];
	}

	/*
	 * Method: clear
	 * -------------
	 *
	 */
	void clear() {
		if(!empty()) {
//			varName.clear();
			if(arrayAD!=NULL) {
				delete [] arrayAD;
				arrayAD=NULL;
			}

			myNbocv2.clear();
			myNcv = 0;
		}

		myIcv = ADVAR_ERROR_CODE;
	}

	/*
	 * Method: empty
	 * -------------
	 *
	 */
	bool empty() const {
		if(arrayAD==NULL) {
			if(isMemSavingMode) {
				if(!myNbocv2.empty())
					cout<<"WARNING! ADtensor::empty(): the array is empty, but (vector<int> myNbocv2) is not empty"<<endl;
			} else {
				if(myNcv!=0)
					cout<<"WARNING! ADtensor::empty(): the array is empty, but (myNcv) is not zero"<<endl;
			}

			return true;
		}
		return false;
	}

	/*
	 * Method: getAt
	 * -------------
	 * For the mutable case and the immutable case
	 */
	T& getAt(const int icv, const int dim1, const int dim2) {
		if(isMemSavingMode) {
			size_t index = convertIcvToIndex(icv);
			assert(dim1>=0 && dim1<3 && dim2>=0 && dim2<3);
			return arrayAD[index][dim1][dim2];
		}

		assert(icv<myNcv && dim1>=0 && dim1<3 && dim2>=0 && dim2<3);
		return arrayAD[icv][dim1][dim2];
	}
	const T& getAt(const int icv, const int dim1, const int dim2) const {
		if(isMemSavingMode) {
			size_t index = convertIcvToIndex(icv);
			assert(dim1>=0 && dim1<3 && dim2>=0 && dim2<3);
			return arrayAD[index][dim1][dim2];
		}

		assert(icv<myNcv && dim1>=0 && dim1<3 && dim2>=0 && dim2<3);
		return arrayAD[icv][dim1][dim2];
	}
//
//	T* getAt(const int icv, const int dim1){
//		if(isMemSavingMode) {
//			size_t index = convertIcvToIndex(icv);
//			assert(dim1>=0 && dim1<3);
//			return arrayAD[index][dim1];
//		}
//
//		assert(icv<myNcv && dim1>=0 && dim1<3);
//		return arrayAD[icv][dim1];
//	}
//	const T* getAt(const int icv, const int dim1) const {
//		if(isMemSavingMode) {
//			size_t index = convertIcvToIndex(icv);
//			assert(dim1>=0 && dim1<3);
//			return arrayAD[index][dim1];
//		}
//
//		assert(icv<myNcv && dim1>=0 && dim1<3);
//		return arrayAD[icv][dim1];
//	}
//
//	T (*getAt(const int icv))[3] {
//		if(isMemSavingMode) {
//			size_t index = convertIcvToIndex(icv);
//			return arrayAD[index];
//		}
//
//		assert(icv<myNcv);
//		return arrayAD[icv];
//	}
//	const T (*getAt(const int icv) const)[3] {
//		if(isMemSavingMode) {
//			size_t index = convertIcvToIndex(icv);
//			return arrayAD[index];
//		}
//
//		assert(icv<myNcv);
//		return arrayAD[icv];
//	}
//
	T (&getAt(const int icv, const int dim1))[3] {
		if(isMemSavingMode) {
			size_t index = convertIcvToIndex(icv);
			assert(dim1>=0 && dim1<3);
			return arrayAD[index][dim1];
		}

		assert(icv<myNcv && dim1>=0 && dim1<3);
		return arrayAD[icv][dim1];
	}
	const T (&getAt(const int icv, const int dim1) const)[3] {
		if(isMemSavingMode) {
			size_t index = convertIcvToIndex(icv);
			assert(dim1>=0 && dim1<3);
			return arrayAD[index][dim1];
		}

		assert(icv<myNcv && dim1>=0 && dim1<3);
		return arrayAD[icv][dim1];
	}

	T (&getAt(const int icv))[3][3] {
		if(isMemSavingMode) {
			size_t index = convertIcvToIndex(icv);
			return arrayAD[index];
		}

		assert(icv<myNcv);
		return arrayAD[icv];
	}
	const T (&getAt(const int icv) const)[3][3] {
		if(isMemSavingMode) {
			size_t index = convertIcvToIndex(icv);
			return arrayAD[index];
		}

		assert(icv<myNcv);
		return arrayAD[icv];
	}

	/*
	 * Method: getMyIcv
	 * ----------------
	 *
	 */
	int getMyIcv() const {
		if(myIcv == ADVAR_ERROR_CODE)
			cout<<"WARNING in ADtensor::getMyIcv(): myIcv has not been defined"<<endl;
		return myIcv;
	}

	/*
	 * Method: size()
	 * --------------
	 *
	 */
	int size() const {
		if(isMemSavingMode)
			return myNbocv2.size();

		return myNcv;
	}

	/*
	 * Method: getName
	 * ---------------
	 *
	 */
	string& getName() {
		return varName;
	}
	const string& getName() const {
		return varName;
	}

	/*
	 * Method: setName
	 * ---------------
	 *
	 */
	void setName(char variableName[]){
		varName.clear();
		varName = variableName;
	}

	/*
	 * Method: getMyNbocv2
	 * -------------------
	 *
	 */
	vector<int>& getMyNbocv2() {
		return myNbocv2;
	}
	const vector<int>& getMyNbocv2() const {
		return myNbocv2;
	}

	/*****************************
	 * Element Selection Operator
	 *****************************/
	// Note: In general, we need the so-called "proxy object" for a multi-dimensional container.
	//       In the case of 3D CFD, we don't need to use it due to the special data structure
	//       (i.e., something like "double (*rhou)[3]").
	//       However, I leave the proxy-object implementation for the future reference.
	//       For the details of the algorithm, see the ADvector class.
#ifdef USE_PROXY_OBJECT_
	/*
	 * Class: MutableReferenceMother
	 * -----------------------------
	 * For the mutable case
	 */
	class MutableReferenceMother {
	public:
		friend class ADtensor;

		/*
		 * Class: MutableReferenceChild
		 * ----------------------------
		 *
		 */
		class MutableReferenceChild {
		public:
			friend class ADtensor;
			friend class MutableReferenceMother;
			T& operator[] (const int dim2) {
				return owner->getAt(icv, dim1, dim2);
			}
		private:
			MutableReferenceChild(ADtensor* owner, int icv, int dim1): owner(owner), icv(icv), dim1(dim1) {};
			ADtensor* const owner;
			const int icv;
			const int dim1;
		};

		MutableReferenceChild operator[] (const int dim1) {
			return MutableReferenceChild(owner, icv, dim1);
		}
	private:
		MutableReferenceMother(ADtensor* owner, int icv): owner(owner), icv(icv) {};
		ADtensor* const owner;
		const int icv;
	};

	/*
	 * Operator: []
	 * ------------
	 * Element selection operator for the mutable case
	 */
	MutableReferenceMother operator[] (int icv) {
		return MutableReferenceMother(this, icv);
	}

	/*
	 * Class: ImmutableReferenceMother
	 * -------------------------------
	 * For the "right values"
	 */
	class ImmutableReferenceMother {
	public:
		friend class ADtensor;

		/*
		 * Class: ImmutableReferenceChild
		 * ------------------------------
		 * For the "right values"
		 */
		class ImmutableReferenceChild {
		public:
			friend class ADtensor;
			friend class ImmutableReferenceMother;
			const T& operator[] (const int dim2) const {
				return owner->getAt(icv, dim1, dim2);
			}
		private:
			ImmutableReferenceChild(ADtensor* owner, int icv, int dim1): owner(owner), icv(icv), dim1(dim1) {};
			ADtensor* const owner;
			const int icv;
			const int dim1;
		};

		const ImmutableReferenceChild operator[] (const int dim1) const {
			return ImmutableReferenceChild(owner, icv, dim1);
		}

	private:
		ImmutableReferenceMother(const ADtensor* owner, int icv): owner(owner), icv(icv) {};
		const ADtensor* const owner;
		const int icv;
	};

	/*
	 * Operator: []
	 * ------------
	 * Element selection operator for the immutable case
	 */
	ImmutableReferenceMother operator[] (int icv) const {
		return ImmutableReferenceMother(this, icv);
	}
#else
	/*
	 * Operator: []
	 * ------------
	 * Element selection operator for the mutable case
	 */
	T (&operator[] (int icv))[3][3] {
		return getAt(icv);
	}

	/*
	 * Operator: []
	 * ------------
	 * Element selection operator for the immutable case
	 */
	T (&operator[] (int icv) const)[3][3] {
		return getAt(icv);
	}
#endif

	/***********
	 * Iterators
	 * ---------
	 * For example, if you want to iterate over the elements of "ADtensor<elemType> mat" and print them out:
	 *   for(ADtensor<elemType>::iterator iter=mat.begin(); iter!=mat.end(); ++iter) {
	 *   	for(int i=0; i<3; ++i) {
	 *   		for(int j=0; j<3; ++j)
	 *   			printf("  %.2e", (*iter)[i][j]); // Note: for the "adouble" elemType: printf("  %.2e", (*iter)[i][j].value());
	 *   		printf("\n");
	 *   	}
	 *   	printf("\n");
	 *   }
	 ***********/
//	class iterator
//	{
//	public:
//		iterator(T (*ptr)[3]) : ptr_(ptr) { }
//		iterator operator ++() { iterator i = *this; ptr_+=3; return i; }
//		iterator operator ++(int junk) { ptr_+=3; return *this; }
//		T& getAt(const int dim1, const int dim2) { assert(dim1>=0 && dim1<3 && dim2>=0 && dim2<3); return ptr_[dim1][dim2]; }
//		class MutableRef {
//		public:
//			friend class iterator;
//			T& operator [](const int dim2) { return owner->getAt(dim1, dim2); }
//		private:
//			MutableRef(iterator* owner, int dim1): owner(owner), dim1(dim1) {};
//			iterator* const owner;
//			const int dim1;
//		};
//		MutableRef operator [](int dim1) { return MutableRef(this, dim1); }
//		bool operator ==(const iterator& rhs) { return ptr_ == rhs.ptr_; }
//		bool operator !=(const iterator& rhs) { return ptr_ != rhs.ptr_; }
//	private:
//		T (*ptr_)[3];
//	};
//
//	class const_iterator
//	{
//	public:
//		const_iterator(T (*ptr)[3]) : ptr_(ptr) { }
//		const_iterator operator ++() { const_iterator i = *this; ptr_+=3; return i; }
//		const_iterator operator ++(int junk) { ptr_+=3; return *this; }
//		const T& getAt(const int dim1, const int dim2) const { assert(dim1>=0 && dim1<3 && dim2>=0 && dim2<3); return ptr_[dim1][dim2]; }
//		class ImmutableRef {
//		public:
//			friend class const_iterator;
//			const T& operator [](const int dim2) const { return owner->getAt(dim1, dim2); }
//		private:
//			ImmutableRef(const const_iterator* owner, int dim1): owner(owner), dim1(dim1) {};
//			const const_iterator* const owner;
//			const int dim1;
//		};
//		ImmutableRef operator [](int dim1) const { return ImmutableRef(this, dim1); }
//		bool operator ==(const const_iterator& rhs) { return ptr_ == rhs.ptr_; }
//		bool operator !=(const const_iterator& rhs) { return ptr_ != rhs.ptr_; }
//	private:
//		T (*ptr_)[3];
//	};
	class iterator
	{
	public:
		iterator(T (*ptr)[3]) : ptr_(ptr) { }
		iterator operator ++() { iterator i = *this; ptr_+=3; return i; }
		iterator operator ++(int junk) { ptr_+=3; return *this; }
		T (*operator *())[3] { return ptr_; }
		bool operator ==(const iterator& rhs) { return ptr_ == rhs.ptr_; }
		bool operator !=(const iterator& rhs) { return ptr_ != rhs.ptr_; }
	private:
		T (*ptr_)[3];
	};

	class const_iterator
	{
	public:
		const_iterator(T (*ptr)[3]) : ptr_(ptr) { }
		const_iterator operator ++() { const_iterator i = *this; ptr_+=3; return i; }
		const_iterator operator ++(int junk) { ptr_+=3; return *this; }
		const T (*operator *() const)[3] { return ptr_; }
		bool operator ==(const const_iterator& rhs) { return ptr_ == rhs.ptr_; }
		bool operator !=(const const_iterator& rhs) { return ptr_ != rhs.ptr_; }
	private:
		T (*ptr_)[3];
	};

	friend class iterator;
	friend class const_iterator;

	iterator begin() {
		return iterator(arrayAD[0]);
	}

	iterator end() {
		return iterator(arrayAD[0] + (this->size())*3);
	}

	const_iterator begin() const {
		return const_iterator(arrayAD[0]);
	}

	const_iterator end() const {
		return const_iterator(arrayAD[0] + (this->size())*3);
	}

	/*
	 * Method: copy
	 * ------------
	 *
	 */
	iterator copy(const_iterator first, const_iterator last, iterator result) {
		while(first!=last) {
			for(int i=0; i<3; ++i)
				for(int j=0; j<3; ++j)
					(*result)[i][j] = (*first)[i][j];

			++result;
			++first;
		}
		return result;
	}

private:
	/*
	 * Method: convertIcvToIndex
	 * -------------------------
	 *
	 */
	size_t convertIcvToIndex(const int icv) const {
		size_t index = ADVAR_ERROR_CODE;
		bool found = false;

		for(size_t i=0; i<myNbocv2.size(); ++i) {
			if(myNbocv2[i] == icv) {
				index = i;
				found = true;
				break;
			}
		}

		if(!found) {
			cerr<<"ERROR in ADtensor::convertIcvToIndex(): Cannot find the given icv = "<<icv<<", neighgor of icvCenter="<<getMyIcv()<<", varName="<<getName()<<endl;
#if ADVAR_DEBUG_LEVEL>0
			cerr<<"                                        Members of myNbocv2 = ";
			for(size_t i=0; i<myNbocv2.size(); ++i)
				cerr<<myNbocv2[i]<<"  ";
			cerr<<endl;
#endif
			assert(false);
		}

		return index;
	}
};

#endif /* ADVAR_H_ */
