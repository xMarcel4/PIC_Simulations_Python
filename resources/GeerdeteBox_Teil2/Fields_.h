#pragma once
#include <iostream>

#include "vec3.h"

template <typename T>
class Field_ {
public:
	// constructor
	Field_(int ni, int nj, int nk) : 
		ni{ ni }, nj{ nj }, nk{ nk }{
		data = new T[ni * nj * nk];	// pointer to array containing all the data on heap

		for (int index = 0; index < ni * nj * nk; index++)
			data[index] = 0;	//initialization with 0
	}

	// destructor, frees memory
	~Field_() {
		if (data == nullptr) return;	// return if unallocated
		delete[] data;
		data = nullptr;	// mark as free
	}

	// copy constructor
	Field_(const Field_& other) : Field_{ other.ni, other.nj, other.nk } {
		for (int i = 0; i < ni * nj * nk; i++)
			data[i] = other(i);
	}

	// move construtor
	Field_(Field_&& other) noexcept :
		ni{ other.ni }, nj{ other.nj }, nk{ other.nk } {
			if (data) this->~Field_();	// deallocate own data
			data = other.data;		// steal the data
			other.data = nullptr;	// invalidate
	}

	// move assigment operator
	Field_& operator=(Field_&& f) {
		if (data) ~Field_;	// deallocatore own data
		data = f.data;
		f.data = nullptr;
		return *this;
	}

	// data acces operator
	// read-only acces to data at index i,j,k
	T r_at(int i, int j, int k) const { return data[i + j * ni + k * ni * nj]; } // read only
	T operator() (int i, int j, int k) const { return data[i + j * ni + k * ni * nj]; }
	// data acces operator (read and write)
	T& w_at(int i, int j, int k) { return data[i + j * ni + k * ni * nj]; } // read and write

	// data access operators with one parameter
	T& operator[] (int i) { return data[i]; }
	T operator()(int i) const { return data[i]; }
	
	// overload the assignment operator
	Field_<T>& operator= (const T s) {
		for (int i = 0; i < ni * nj * nk; i++)
			data[i] = s; // set all values to scalar s
		return *this; // return refernce to self
	}

	// elementwise division by another field
	void operator /=(const Field_& other) {
		for (int i = 0; i < ni*nj*nk; i++) {
			if (other.data[i] != 0)
				data[i] /= other(i);
			else 
				data[i] = 0;
		}
	}
	
	Field_& operator += (const Field_& other) {
		for (int i = 0; i < ni * nj * nk; i++)
			data[i] += other(i);
		return (*this);
	}

	// compound multiplication
	Field_& operator *= (double s) {
		for (int i = 0; i < ni * nj * nk; i++)
			data[i] *= s;
		return *this;
	}

	// multiplikation operator, returns new Field set to f*s
	friend Field_<T> operator*(double s, const Field_<T>& f) {
		Field_<T> r(f);
		return std::move(r *= s);	// force move
	}

	// writes data to a file stream
	friend std::ostream& operator<<(std::ostream& out, Field_<T>& f) {
		for (int k = 0; k < f.nk; k++, out << "\n")	// new line after each "k"
			for (int j = 0; j < f.nj; j++)
				for (int i = 0; i < f.ni; i++)
					out << f.data[i + j * f.ni + k * f.ni * f.nj] << " ";
		return out;
	}

	void scatter(double3 lc, double value);
	
	// gathers field value at a logical coordinate lc
	T gather(double3 lc);

protected:
	const int ni, nj, nk; // number of nodes in x, y, z
	T* data = nullptr;	// pointer of type T
};

using Field = Field_<double>;	// field of doubles
using FieldI = Field_<int>;		// field of integers
using Field3 = Field_<double3>;	// vector field of doubles


template<typename T>
void Field_<T>::scatter(double3 lc, double value) {
	// make sure we are in domain
	if (lc[0]<0 || lc[0]>(static_cast<double>(ni) - 1) || lc[1]<0 || (lc[1]>static_cast<double>(nj) - 1) || lc[2]<0 || lc[2]>(static_cast<double>(nk) - 1)) return;

	// compute the cell index and the fractional distances
	int i = (int)lc[0];
	double di = lc[0] - i;
	int j = (int)lc[1];
	double dj = lc[1] - j;
	int k = (int)lc[2];
	double dk = lc[2] - k;

	// deposit fractional values to the 8 surrounding nodes
	data[i + j*ni + k*ni*nj] += value * (1 - di) * (1 - dj) * (1 - dk);
	data[(i+1) + j*ni + k*ni*nj] += value * (di) * (1 - dj) * (1 - dk);
	data[i + (j+1)*ni + k*ni*nj] += value * (1 - di) * (dj) * (1 - dk);
	data[(i+1) + (j+1)*ni + k*ni*nj] += value * (di) * (dj) * (1 - dk);
	data[i + j*ni + (k+1)*ni*nj] += value * (1 - di) * (1 - dj) * (dk);
	data[(i+1) + j*ni + (k+1)*ni*nj] += value * (di) * (1 - dj) * (dk);
	data[i + (j+1)*ni + (k+1)*ni*nj] += value * (1 - di) * (dj) * (dk);
	data[(i+1) + (j+1)*ni + (k+1)*ni*nj] += value * (di) * (dj) * (dk);
}

template<typename T>
T Field_<T>::gather(double3 lc) {
	// compute the cell index and the fractional distances
	int i = (int)lc[0];
	double di = lc[0] - i;
	int j = (int)lc[1];
	double dj = lc[1] - j;
	int k = (int)lc[2];
	double dk = lc[2] - k;

	// interpolate data onto particle position
	T val = r_at(i, j, k) * (1 - di) * (1 - dj) * (1 - dk) +
		r_at(i + 1, j, k) * (di) * (1 - dj) * (1 - dk) +
		r_at(i, j + 1, k) * (1 - di) * (dj) * (1 - dk) +
		r_at(i + 1, j + 1, k) * (di) * (dj) * (1 - dk) +
		r_at(i, j, k + 1) * (1 - di) * (1 - dj) * (dk)+
		r_at(i + 1, j, k + 1) * (di) * (1 - dj) * (dk)+
		r_at(i, j + 1, k + 1) * (1 - di) * (dj) * (dk)+
		r_at(i + 1, j + 1, k + 1) * (di) * (dj) * (dk);

	return val;
}