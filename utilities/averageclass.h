#ifndef AVERAGECLASS_H
#define AVERAGECLASS_H

#include "complexgrid.h"

template <class T> class AverageClass {
	private:
		T x;
		T x2;
		uint32_t n;
	public:
		AverageClass() : x(), x2()
		{
			n = 0;
		}

		AverageClass(const T &value)
		{
			x = value;
			x2 = value*value;
			n = 1;
		}
		
		AverageClass(const AverageClass<T> &av)
		{
			x = av.x;
			x2 = av.x2;
			n = av.n;
		}
		
		AverageClass<T> & operator= (const T &value)
		{
			x = value;
			x2 = value*value;
			n = 1;
			return *this;
		}
		
		AverageClass<T> & operator= (const AverageClass<T> &av)
		{
			x = av.x;
			x2 = av.x2;
			n = av.n;
			return *this;
		}
		
		AverageClass<T> & operator+= (const T &value)
		{
			if(n==1)
			{
				x += value;
				x2 = x*x;
			}
			else if(n==0)
			{
				x = value;
				x2 = x*x;
				n = 1;
			}
			else
			{
				cout << "Warning: Invalid summation in AverageClass!" << endl;
			}
			return *this;
		}
		
		AverageClass<T> & operator+= (const AverageClass<T> &v)
		{
			if((v.n==1) && (n==1))
			{
				x += v.x;
				x2 = x*x;
			}
			else if(n==0)
			{
				x = v.x;
				x2 = v.x2;
				n = v.n;
			}
			else
			{
				cout << "Warning: Invalid summation in AverageClass!" << endl;
			}
			return *this;
		}
		
		AverageClass<T> & operator/= (const T &d)
		{
			x /= d;
			x2 /= d*d;
			return *this;
		}
		
		AverageClass<T> operator+ (const T &value)
		{
			AverageClass<T> ret;
			if(n==1)
			{
				ret.x = x + value;
				ret.x2 = ret.x * ret.x;
				ret.n = 1;
			}
			else if(n==0)
			{
				ret.x = value;
				ret.x2 = value*value;
				ret.n = 1;
			}
			else
			{
				cout << "Warning: Invalid summation in AverageClass!" << endl;
			}
			return ret;
		}
		
		AverageClass<T> operator+ (const AverageClass<T> &v)
		{
			AverageClass<T> ret;
			if((n==1) && (v.n==1))
			{
				ret.x = x + v.x;
				ret.x2 = ret.x * ret.x;
				ret.n = 1;
			}
			else if(n==0)
			{
				ret.x = v.x;
				ret.x2 = v.x2;
				ret.n = v.n;
			}
			else if(v.n==0)
			{
				ret.x = x;
				ret.x2 = x2;
				ret.n = n;
			}
			else
			{
				cout << "Warning: Invalid summation in AverageClass!" << endl;
			}
			return ret;
		}
		
		AverageClass<T> & operator++ ()
		{
			if(n==1)
			{
				x += 1;
				x2 = x*x;
			}
			else if(n==0)
			{
				x = 1;
				x2 = 1;
				n = 1;
			}
			else
			{
				cout << "Warning: Invalid increment in AverageClass!" << endl;
			}
			return *this;
		}
		
		operator T()
		{
			if(n==0)
				return x;
			else
				return x/(double)n;
		}
		
		T av() const
		{
			if(n==0)
				return x;
			else
				return x/(double)n;
		}
		
		T av_unorm() const
		{
	
				return x;

		}
		
		T sd() const
		{
			if(n==0)
				return x;
			else
				return sqrt(x2/(double)n - x/(double)n*x/(double)n);
		}
		
		void average(const T &v)
		{
			if(n==0)
			{
				x = v;
				x2 = v*v;
				n = 1;
			}
			else
			{
				x += v;
				x2 += v*v;
				n++;
			}
		}
		
		void average_weight(const T &v, const double weight)
		{
			if(n==0)
			{
				x = weight*v;
				x2 = weight*v*v;
				n = weight;
			}
			else
			{
				x += weight*v;
				x2 += weight*v*v;
				n += weight;
			}
		}
		
		void average(const AverageClass<T> &v)
		{
			if(n==0)
			{
				x = v.x;
				x2 = v.x2;
				n = v.n;
			}
			else
			{
				x += v.x;
				x2 += v.x2;
				n += v.n;
			}
		}
		

		
		void reset()
		{
			x = 0;
			x2 = 0;
			n = 0;
		}
		
		friend AverageClass<T> conj(const AverageClass<T> &a) {return a;}
};

template <typename T>
inline AverageClass<complex<T> > conj(const AverageClass<complex<T> > &a)
{
	AverageClass<complex<T> > res;
	res.x = conj(a.x);
	res.x2 = conj(a.x2);
	res.n = a.n;
	return res;
}

template <class T>
ostream & operator<< (ostream &o, const AverageClass<T> & av)
{
	o << av.av() << " " << av.sd();
	return o;
}

#endif
