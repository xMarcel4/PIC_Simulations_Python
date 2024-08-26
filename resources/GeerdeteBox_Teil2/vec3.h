#pragma once
#include <ostream>

template <typename T>
struct vec3 {
	// constructor: assign data
	vec3(const T u, const T v, const T w) : d{ u,v,w } {}
	vec3(const T a[3]) : d{ a[0],a[1],a[2] } {}
	vec3() : d{ 0,0,0 } {}
	
	// read (and write) data
	T& operator[] (int i) { return d[i]; } // r&w
	T operator() (int i) const { return d[i]; } // r

	// operator overload for calculations
	vec3<T>& operator=(double s) {
		d[0] = s;
		d[1] = s;
		d[2] = s;
		return *this;
	}

	vec3<T>& operator+=(vec3<T> o) {
		d[0] += o[0];
		d[1] += o[1];
		d[2] += o[2];
		return *this;
	}

	vec3<T>& operator-=(vec3<T> o) {
		d[0] -= o[0];
		d[1] -= o[1];
		d[2] -= o[2];
		return *this;
	}
	
protected:
	T d[3] = {};
};

// vec3-vec3 operations
template<typename T>	// addition of two vec3s
vec3<T> operator+(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T>(a(0) + b(0), a(1) + b(1), a(2) + b(2)); }
template<typename T>	// subtraction of two vec3s
vec3<T> operator-(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T>(a(0) - b(0), a(1) - b(1), a(2) - b(2)); }
template<typename T>	// elementwise multiplication of two vec3s
vec3<T> operator*(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T>(a(0) * b(0), a(1) * b(1), a(2) * b(2)); }
template<typename T>	// elementwise davivsion of two vec3s
vec3<T> operator/(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T>(a(0) / b(0), a(1) / b(1), a(2) / b(2)); }

// vec3-scalar operations
template<typename T>	// scalar multiplication
vec3<T> operator*(const vec3<T>& a, const T& s) {
	return vec3<T>(a(0) * s, a(1) * s, a(2) * s); }
template<typename T>	// scalar multiplication 2
vec3<T> operator*(const T& s, const vec3<T>& a) {
	return vec3<T>(a(0) * s, a(1) * s, a(2) * s); }

// output
template<typename T>
std::ostream& operator<<(std::ostream& out, const vec3<T>& v) {
	out << v(0) << " " << v(1) << " " << v(2);
	return out;
}

using double3 = vec3<double>;
using int3 = vec3<int>;
