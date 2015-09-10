#ifndef COORDINATE_H__
#define COORDINATE_H__

#include <stdint.h>
#include <list>
#include <iostream>
#include <cmath>
#include <bh3binaryfile.h>

using namespace std;



template <typename T>
class Vector;

template <typename T>
class Coordinate;

template <typename T>
void write(ostream &o, const Coordinate<T> &c);

template <typename T>
void read(istream &i, Coordinate<T> &c);

template <typename T>
void write(ostream &o, const Vector<T> &c);

template <typename T>
void read(istream &i, Vector<T> &c);

template <typename T>
class Coordinate {
	private:
		T ix[3];
		int32_t dim[3];
	public:
		Coordinate() {ix[0] = 0; ix[1] = 0; ix[2] = 0; dim[0] = 0; dim[1] = 0; dim[2] = 0;}
		Coordinate(T nx, T ny, T nz, uint32_t nw, uint32_t nh, uint32_t nd);
		
		Vector<T> operator-(const Coordinate<T> &c) const;
		Coordinate<T> operator+ (const Vector<T> &v) const;
		Coordinate<T> operator- (const Vector<T> &v) const;
		Coordinate<T> & operator+= (const Vector<T> &v);
		Coordinate<T> & operator-= (const Vector<T> &v);
		bool operator== (const Coordinate<T> &c) const;
		bool operator!= (const Coordinate<T> &c) const;
		
		T x() const;
		T y() const;
		T z() const;
		T operator[] (int index) const;
		
		template <typename T2> operator Coordinate<T2>() const;
		
		static Coordinate<double> average(const list<Coordinate<T> > &vl);
		
		friend void write<T>(ostream &o, const Coordinate<T> &c);
		friend void read<T>(istream &i, Coordinate<T> &c);
};

template <typename T>
class Vector {
	private:
		T ix[3];
		int32_t dim[3];
	public:
		Vector() {for(int i = 0; i < 3; i++) {ix[i] = 0; dim[i] = 0;}}
		Vector(T nx, T ny, T nz, uint32_t nw, uint32_t nh, uint32_t nd);
		
		Vector<T> operator- (const Vector<T> &v) const;
		Vector<T> operator+ (const Vector<T> &v) const;
		Coordinate<T> operator+ (const Coordinate<T> &c) const;
		Vector<T> operator* (T r) const;
		Vector<T> operator/ (T r) const;
		Vector<T> & operator*= (T r);
		Vector<T> & operator/= (T r);
		Vector<T> & operator-= (const Vector<T> &v);
		Vector<T> & operator+= (const Vector<T> &v);
		bool operator== (const Vector<T> &v) const;
		bool operator!= (const Vector<T> &v) const;
		
		T x() const;
		T y() const;
		T z() const;
		T operator[] (int index) const;
		
		template <typename T2> operator Vector<T2>() const;
		
		T norm() const;
		friend void write<T>(ostream &o, const Vector<T> &c);
		friend void read<T>(istream &i, Vector<T> &c);
};

template <typename T>
inline Coordinate<T>::Coordinate(T nx, T ny, T nz, uint32_t nw, uint32_t nh, uint32_t nd)
{
	ix[0] = nx;
	ix[1] = ny;
	ix[2] = nz;
	dim[0] = nw;
	dim[1] = nh;
	dim[2] = nd;
}

template <typename T> template <typename T2>
Coordinate<T>::operator Coordinate<T2>() const
{
	return Coordinate<T2>(ix[0], ix[1], ix[2], dim[0], dim[1], dim[2]);
}

template <typename T>
inline Vector<T> Coordinate<T>::operator-(const Coordinate<T> &c) const
{
	Vector<T> v(ix[0] - c.ix[0], ix[1] - c.ix[1], ix[2] - c.ix[2], dim[0], dim[1], dim[2]);
	return Vector<T>(v.x(), v.y(), v.z(), dim[0], dim[1], dim[2]);
}

template <typename T>
inline Coordinate<T> Coordinate<T>::operator+ (const Vector<T> &v) const
{
	return Coordinate<T>(ix[0] + v.x(), ix[1] + v.y(), ix[2] + v.z(), dim[0], dim[1], dim[2]);
}

template <typename T>
inline Coordinate<T> Coordinate<T>::operator- (const Vector<T> &v) const
{
	return Coordinate<T>(ix[0] - v.x(), ix[1] - v.y(), ix[2] - v.z(), dim[0], dim[1], dim[2]);
}

template <typename T>
inline Coordinate<T> & Coordinate<T>::operator+= (const Vector<T> &v)
{
	ix[0] += v.x();
	ix[1] += v.y();
	ix[2] += v.z();
	return *this;
}

template <typename T>
inline Coordinate<T> & Coordinate<T>::operator-= (const Vector<T> &v)
{
	ix[0] -= v.x();
	ix[1] -= v.y();
	ix[2] -= v.z();
	return *this;
}

template <typename T>
inline T Coordinate<T>::operator[] (int index) const
{
	if(dim[index] == 0)
	{
		cout << "Warning: Coordinate operator[" << index << "] with invalid dimension!" << endl;
		return 0;
	}
	
	T ret = ix[index];
	while(ret >= dim[index])
		ret -= dim[index];
	while(ret < 0)
		ret += dim[index];
	return ret;
}

template <typename T>
inline T Coordinate<T>::x() const
{
	return (*this)[0];
}

template <typename T>
inline T Coordinate<T>::y() const
{
	return (*this)[1];
}

template <typename T>
inline T Coordinate<T>::z() const
{
	return (*this)[2];
}

template <typename T>
inline bool Coordinate<T>::operator== (const Coordinate<T> &c) const
{
	return (x() == c.x()) && (y() == c.y()) && (z() == c.z());
}

template <typename T>
inline bool Coordinate<T>::operator!= (const Coordinate<T> &c) const
{
	return !(*this == c);
}

template <typename T>
inline ostream &operator<< (ostream &str, const Coordinate<T> &c)
{
	str << c.x() << " " << c.y() << " " << c.z();
	return str;
}

template <typename T>
inline Vector<T>::Vector(T nx, T ny, T nz, uint32_t nw, uint32_t nh, uint32_t nd)
{
	ix[0] = nx;
	ix[1] = ny;
	ix[2] = nz;
	dim[0] = nw;
	dim[1] = nh;
	dim[2] = nd;
}

template <typename T> template <typename T2>
Vector<T>::operator Vector<T2>() const
{
	return Vector<T2>(ix[0], ix[1], ix[2], dim[0], dim[1], dim[2]);
}

template <typename T>
inline Vector<T> Vector<T>::operator- (const Vector<T> &v) const
{
	return Vector<T>(ix[0] - v.ix[0], ix[1] - v.ix[1], ix[2] - v.ix[2], dim[0], dim[1], dim[2]);
}

template <typename T>
inline Vector<T> Vector<T>::operator+ (const Vector<T> &v) const
{
	return Vector<T>(ix[0] + v.ix[0], ix[1] + v.ix[1], ix[2] + v.ix[2], dim[0], dim[1], dim[2]);
}

template <typename T>
inline Coordinate<T> Vector<T>::operator+ (const Coordinate<T> &c) const
{
	return Coordinate<T>(ix[0] + c.x(), ix[1] + c.y(), ix[2] + c.z(), dim[0], dim[1], dim[2]);
}

template <typename T>
inline Vector<T> Vector<T>::operator* (T r) const
{
	return Vector<T>(r * ix[0], r * ix[1], r * ix[2], dim[0], dim[1], dim[2]);
}

template <typename T>
inline Vector<T> Vector<T>::operator/ (T r) const
{
	return Vector<T>(ix[0] / r, ix[1] / r, ix[2] / r, dim[0], dim[1], dim[2]);
}

template <typename T>
inline Vector<T> & Vector<T>::operator*= (T r)
{
	ix[0] *= r;
	ix[1] *= r;
	ix[2] *= r;
	return *this;
}

template <typename T>
inline Vector<T> & Vector<T>::operator/= (T r)
{
	ix[0] /= r;
	ix[1] /= r;
	ix[2] /= r;
	return *this;
}

template <typename T>
inline Vector<T> & Vector<T>::operator-= (const Vector<T> &v)
{
	ix[0] -= v.ix[0];
	ix[1] -= v.ix[1];
	ix[2] -= v.ix[2];
	return *this;
}

template <typename T>
inline Vector<T> & Vector<T>::operator+= (const Vector<T> &v)
{
	ix[0] += v.ix[0];
	ix[1] += v.ix[1];
	ix[2] += v.ix[2];
	return *this;
}

template <typename T>
inline T Vector<T>::operator[] (int index) const
{

	if(dim[index] == 0)
	{
		cout << "Warning: Vector operator[" << index << "] with invalid dimension!" << endl;
		return 0;
	}
	
	T ret = ix[index];
	if(dim[index] != 1)
	{

		while(ret >= dim[index]/2)
			ret -= dim[index];

		while(ret < -dim[index]/2)
			ret += dim[index];

		
	}
	else
	{
		ret = 0;
	}
	return ret;
}

template <typename T>
inline T Vector<T>::x() const
{
	return (*this)[0];
}

template <typename T>
inline T Vector<T>::y() const
{
	return (*this)[1];
}

template <typename T>
inline T Vector<T>::z() const
{
	return (*this)[2];
}

template <typename T>
inline T Vector<T>::norm() const
{
	return sqrt(ix[0]*ix[0] + ix[1]*ix[1] + ix[2]*ix[2]);
}

template <typename T>
inline bool Vector<T>::operator== (const Vector<T> &v) const
{
	return (x() == v.x()) && (y() == v.y()) && (z() == v.z());
}

template <typename T>
inline bool Vector<T>::operator!= (const Vector<T> &v) const
{
	return !(*this == v);
}

template <typename T>
Coordinate<double> Coordinate<T>::average(const list<Coordinate<T> > &cl)
{
	Coordinate<double> c = cl.front();
	Vector<double> v = Vector<double>(0,0,0,c.dim[0],c.dim[1],c.dim[2]);
	_List_const_iterator<Coordinate<T> > iter;
	for(iter = cl.begin(); iter != cl.end(); iter++)
	{
		v += *iter - c;
	}
	v /= cl.size();
	return c + v;
}

template<typename T>
inline void write(ostream &o, const Coordinate<T> &c)
{
	for(int pos = 0; pos < 3; pos++)
	{
		write(o, c.dim[pos]);
		write(o, c.ix[pos]);
	}
}

template <typename T>
inline void read(istream &i, Coordinate<T> &c)
{
	for(int pos = 0; pos < 3; pos++)
	{
		read(i, c.dim[pos]);
		read(i, c.ix[pos]);
	}
}

template<typename T>
inline void write(ostream &o, const Vector<T> &v)
{
	for(int pos = 0; pos < 3; pos++)
	{
		write(o, v.dim[pos]);
		write(o, v.ix[pos]);
	}
}

template <typename T>
inline void read(istream &i, Vector<T> &v)
{
	for(int pos = 0; pos < 3; pos++)
	{
		read(i, v.dim[pos]);
		read(i, v.ix[pos]);
	}
}

namespace std {
	template <>
	struct hash<Coordinate<int32_t>>{
		size_t operator() (const Coordinate<int32_t> &c) const{
			string tmp = to_string(c.x()) + to_string(c.y()); // + to_string(c.z()); // OBVIOUS 2D ONLY
			std::istringstream istr(tmp);
			size_t hash;
			istr >> hash;
			return hash;
		}
	};
}

#endif
