/* Array.h:  A high-performance multi-dimensional C++ array class
Copyright (C) 1997-2005 John C. Bowman

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

#ifndef __Array_h__
#define __Array_h__ 1

#define __ARRAY_H_VERSION__ 1.39

// Defining NDEBUG improves optimization but disables argument checking.
// Defining __NOARRAY2OPT inhibits special optimization of Array2[].

#include <iostream>
#include <sstream>
#include <climits>
#include <cstdlib>
#include <cerrno>
#include <stdexcept>
#include <math.h>
#include <float.h>
#ifdef NDEBUG
#define __check(i,n,dim,m)
#define __checkSize()
#define __checkEqual(a,b,dim,m)
#define __checkActivate(i,align) this->Activate(align)
#else
#define __check(i,n,dim,m) this->Check(i,n,dim,m)
#define __checkSize() this->CheckSize()
#define __checkEqual(a,b,dim,m) this->CheckEqual(a,b,dim,m)
#define __checkActivate(i,align) this->CheckActivate(i,align)
#ifndef __NOARRAY2OPT
#define __NOARRAY2OPT
#endif
#endif
#include <crtdbg.h>

using std::istream;
using std::ostream;
using std::ostringstream;
using std::cin;
using std::cerr;
using std::endl;
using std::flush;

namespace Array 
{

	inline ostream& _newl(ostream& s) {s << '\n'; return s;}

	inline void ArrayExit(const char *x);

#ifndef __ExternalArrayExit
	inline void ArrayExit(const char *x)
	{
		cerr << _newl << "ERROR: " << x << "." << endl;
		throw std::runtime_error(x);
	} 
#endif  
}

#ifndef HAVE_POSIX_MEMALIGN

#ifdef __GLIBC_PREREQ
#if __GLIBC_PREREQ(2,3)
#define HAVE_POSIX_MEMALIGN
#endif
#else
#ifdef _POSIX_SOURCE
#define HAVE_POSIX_MEMALIGN
#endif
#endif

#endif

#ifdef HAVE_POSIX_MEMALIGN
extern "C" int posix_memalign(void **memptr, size_t alignment, size_t size);
#else
#include <malloc.h>
// Adapted from FFTW aligned malloc/free.  Assumes that malloc is at least
// sizeof(void*)-aligned. Allocated memory must be freed with free0.
inline int posix_memalign0(void **memptr, size_t alignment, size_t size)
{
	if(alignment % sizeof (void *) != 0 || (alignment & (alignment - 1)) != 0)
		return EINVAL;
	void *p0=malloc(size+alignment);
	if(!p0) return ENOMEM;
	void *p=(void *)(((size_t) p0+alignment)&~(alignment-1));
	*((void **) p-1)=p0;
	*memptr=p;
	return 0;
}
inline void free0(void *p)
{
	if(p) free(*((void **) p-1));
}
#endif

template<class T> inline void newAlign(T *&v, size_t len, size_t align)
{
	void *mem = 0;
	const char *invalid = "Invalid alignment requested";
	const char *nomem = "Memory limits exceeded";
#ifdef HAVE_POSIX_MEMALIGN
	int rc = posix_memalign(&mem, align, len*sizeof(T));
#else  
	int rc = posix_memalign0(&mem, align, len*sizeof(T));
#endif  
	if(rc == EINVAL) 
		Array::ArrayExit(invalid);
	if(rc == ENOMEM) 
		Array::ArrayExit(nomem);
	v = (T *) mem;
	for(size_t i=0; i < len; i++) 
		new(v+i) T;
}

template<class T> inline void deleteAlign(T *v, size_t len)
{
	for(size_t i = len-1; i > 0; i--) 
		v[i].~T();
	v[0].~T();
#ifdef HAVE_POSIX_MEMALIGN
	free(v);
#else
	free0(v);
#endif  
}

namespace Array 
{
	template <class T> class array2; // forward declaration

	template<class T> class array1 
	{
		friend array2<T>;
	public:
		array1() : m_size(0), m_v(0), m_state(e_unallocated), m_refCount(0) {}
		array1(unsigned int nx0, size_t align = 0) : m_state(e_unallocated), m_refCount(0) 
		{
			Allocate(nx0, align);
		}
		array1(unsigned int nx0, T *v0) : m_state(e_unallocated), m_refCount(1) 
		{ 
			Dimension(nx0, v0); 
		}
		array1(const array1<T>& A) :
			m_v(A.m_v), 
			m_size(A.m_size), 
			m_state(A.test(e_temporary)),
			m_refCount(1)
		{}

		virtual ~array1() {	Deallocate(); }

	protected:
		T *m_v;
		unsigned int m_size;
		mutable int m_state;
		unsigned int m_refCount;
		enum alloc_state {e_unallocated=0, e_allocated=1, e_temporary=2, e_aligned=4};

	public:
		virtual unsigned int Size() const {return m_size;}
		void CheckSize() const 
		{
			if(!test(e_allocated) && m_size == 0)
				ArrayExit("Operation attempted on e_unallocated array"); 
		}
		void CheckEqual(int a, int b, unsigned int dim, unsigned int m) const 
		{
			if(a != b) 
			{
				ostringstream buf;
				buf << "Array" << dim << " index ";
				if(m) buf << m << " ";
				buf << "is incompatible in assignment (" << a << " != " << b << ")";
				ArrayExit(buf.str().c_str());
			}
		}

		int test(int flag) const {return m_state & flag;}
		void clear(int flag) const {m_state &= ~flag;}
		void orstate(int flag) const {m_state |= flag;}
		void set(int flag) const {m_state = flag;}
		void Activate(size_t align=0) 
		{
			if(align) 
			{
				newAlign(m_v,m_size,align);
				orstate(e_allocated | e_aligned);
			} 
			else 
			{
				m_v = new T[m_size];
				orstate(e_allocated);
			}
		}
		void CheckActivate(int dim, size_t align=0) 
		{
			if (test(e_allocated)) 
			{
				ostringstream buf;
				buf << "Reallocation of Array" << dim
					<< " attempted (must Deallocate first)";
				ArrayExit(buf.str().c_str());
			}
			Activate(align);
		}
		void Deallocate() 
		{
			if(test(e_allocated)) 
			{
				if(test(e_aligned)) 
					deleteAlign(m_v, m_size);
				else
				{
					delete [] m_v;
					m_v = 0;
				}
				set(e_unallocated);
			}
		}
		void Dimension(unsigned int nx0) 
		{
			m_size = nx0;
		}

		void Dimension(unsigned int nx0, T *v0) 
		{
			Dimension(nx0); 
			m_v = v0; 
			clear(e_allocated);
		}

		void Dimension(const array1<T>& A) 
		{
			Dimension(A.m_size, A.m_v); 
			set(A.test(e_temporary));
		}

		void Allocate(unsigned int nx0, size_t align=0) 
		{
			Dimension(nx0);
			__checkActivate(1,align);
		}

		void Reallocate(unsigned int nx0, size_t align=0) 
		{
			Deallocate();
			Allocate(nx0,align);
		}

		void Freeze() 
		{
			set(e_unallocated);
		}
		void Hold() 
		{
			if(test(e_allocated)) 
				set(e_temporary);
		}
		void Release() 
		{
			if(test(e_temporary)) 
				set(e_allocated);
		}
		void Purge()
		{
			if(test(e_temporary)) 
			{
				Deallocate(); 
				set(e_unallocated);
			}
		}
		void Ramp(int step = 1)
		{
			__checkSize();
			for(unsigned int i = 0; i < m_size; i += step) 
				m_v[i] = i;
		}

		T maxDiff()
		{
			T tMax = -DBL_MAX;
			for (unsigned int ii = m_size - 2;ii > 0;ii--)
			{
				T tDiff = m_v[ii+1] - m_v[ii];
				if (tDiff > tMax)
					tMax = tDiff;
			}
			if (m_size)
			{
				T tDiff = m_v[1] - m_v[0];
				if (tDiff > tMax)
					tMax = tDiff;
			}

			return tMax;
		}

		virtual void Check(int i, int n, unsigned int dim, unsigned int m,
			int o=0) const 
		{
			if(i < 0 || i >= n) 
			{
				ostringstream buf;
				buf << "Array" << dim << " index ";
				if(m) buf << m << " ";
				buf << "is out of bounds (" << i+o;
				if(n == 0) 
					buf << " index given to empty array";
				else 
				{
					if(i < 0) 
						buf << " < " << o;
					else 
						buf << " > " << n+o-1;
				}
				buf << ")";
				ArrayExit(buf.str().c_str());
			}
		}

		unsigned int Nx() const {return m_size;}

#ifdef NDEBUG
		typedef T *opt;
#else
		typedef array1<T> opt;
#endif

		T& operator [] (int ix) const {__check(ix,m_size,1,1); return m_v[ix];}
		T& operator () (int ix) const {__check(ix,m_size,1,1); return m_v[ix];}
		T* operator () () const {return m_v;}
		operator T* () const {return m_v;}

		array1<T> operator + (int i) const {return array1<T>(m_size-i,m_v+i);}

		void Load(T a) const 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) 
				m_v[i]=a;
		}
		void Load(const T *a) const 
		{
			for(unsigned int i=0; i < m_size; i++) 
				m_v[i]=a[i];
		}
		void Store(T *a) const 
		{
			for(unsigned int i=0; i < m_size; i++) 
				a[i]=m_v[i];
		}
		void Set(T *a) {m_v=a; clear(e_allocated);}
		T Min() 
		{
			if(m_size == 0)
				ArrayExit("Cannot take minimum of empty array"); 
			T min = m_v[0];
			for(unsigned int i = 1; i < m_size; i++) 
				if(m_v[i] < min) 
					min = m_v[i];
			return min;
		}
		T Max() 
		{
			if(m_size == 0)
				ArrayExit("Cannot take maximum of empty array"); 
			T max=m_v[0];
			for(unsigned int i=1; i < m_size; i++) 
				if(m_v[i] > max) 
					max=m_v[i];
			return max;
		}
		int MaxIdx() 
		{
			if(m_size == 0)
				ArrayExit("Cannot take maximum of empty array"); 
			T max = m_v[0];
			int maxidx = 0;
			for(unsigned int i=1; i < m_size; i++) 
				if(m_v[i] > max) 
				{ 
					max=m_v[i];
					maxidx = i;
				}
				return maxidx;
		}
		T Mean() 
		{
			if(m_size == 0)
				ArrayExit("Cannot take mean of empty array"); 
			T mean = 0;
			for(unsigned int i = 0; i < m_size; i++) 
				mean += m_v[i];
			return mean / m_size;
		}

		array1<T> Histogram(int nHist)
		{
			array1<T> arrHist(nHist); 
			arrHist = 0.0;
			RTYPE mmin = Min();
			RTYPE Range = Max() - mmin;
			RTYPE BinSize = Range / (nHist - 1);
			for (unsigned int ii = 0; ii < m_size; ii++))
			{
					// hist algorithm
					RTYPE val = m_v[i];
					int nBin = (int)((val - mmin) / BinSize);
					if (nBin >= nHist || nBin < 0)
						DebugBreak();
					arrHist.m_v[nBin]++;
				}
			}
			arrHist.Hold();
			return arrHist;
		}

		istream& Input (istream &s) const 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) 
				s >> m_v[i];
			return s;
		}

		array1<T>& operator = (T a) {Load(a); return *this;}
		array1<T>& operator = (const T *a) {Load(a); return *this;}
		array1<T>& operator = (const array1<T>& A) 
		{
			__checkEqual(m_size,A.Size(),1,1);
			Load(A());
//			A.Purge();
			return *this;
		}

		array1<T>& operator += (const array1<T>& A) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) 
				m_v[i] += A(i);
			return *this;
		}
		array1<T>& operator -= (const array1<T>& A) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) 
				m_v[i] -= A(i);
			return *this;
		}
		array1<T>& operator *= (const array1<T>& A) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) 
				m_v[i] *= A(i);
			return *this;
		}
		array1<T>& operator /= (const array1<T>& A) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) 
				m_v[i] /= A(i);
			return *this;
		}

		array1<T>& operator += (T a) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) 
				m_v[i] += a;
			return *this;
		}
		array1<T>& operator -= (T a) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) 
				m_v[i] -= a;
			return *this;
		}
		array1<T>& operator *= (T a) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) 
				m_v[i] *= a;
			return *this;
		}
		array1<T>& operator /= (T a) 
		{
			__checkSize();
			T ainv=1.0/a;
			for(unsigned int i=0; i < m_size; i++) 
				m_v[i] *= ainv;
			return *this;
		}

		double L1() const 
		{
			__checkSize();
			double norm=0.0;
			for(unsigned int i=0; i < m_size; i++) norm += abs(m_v[i]);
			return norm;
		}
#ifdef __ArrayExtensions
		double Abs2() const 
		{
			__checkSize();
			double norm=0.0;
			for(unsigned int i=0; i < m_size; i++) norm += abs2(m_v[i]);
			return norm;
		}
		double L2() const 
		{
			return sqrt(Abs2());
		}
		double LInfinity() const 
		{
			__checkSize();
			double norm=0.0;
			for(unsigned int i=0; i < m_size; i++) 
			{
				T a = abs(m_v[i]);
				if(a > norm) 
					norm = a;
			}
			return norm;
		}

		double LMinusInfinity() const 
		{
			__checkSize();
			double norm=DBL_MAX;
			for(unsigned int i=0; i < m_size; i++) 
			{
				T a = abs(m_v[i]);
				if(a < norm) 
					norm = a;
			}
			return norm;
		}
#endif	
	};

	template<class T> void swaparray(T& A, T& B)
	{
		T C;
		C.Dimension(A);
		A.Dimension(B);
		B.Dimension(C);
	}

	template <class T> void leftshiftarray(T& A, T& B, T& C)
	{
		T D;
		D.Dimension(A);
		A.Dimension(B);
		B.Dimension(C);
		C.Dimension(D);
	}

	template<class T> void rightshiftarray(T& A, T& B, T& C)
	{
		T D;
		D.Dimension(C);
		C.Dimension(B);
		B.Dimension(A);
		A.Dimension(D);
	}

	template<class T> ostream& operator << (ostream& s, const array1<T>& A)
	{
		T *p=A();
		for(unsigned int i=0; i < A.Nx(); i++) 
		{
			s << *(p++) << " ";
		}
		return s;
	}

	template<class T> istream& operator >> (istream& s, const array1<T>& A)
	{
		return A.Input(s);
	}

	template<class T> class array2 : public array1<T> 
	{
	protected:
		unsigned int nx;
		unsigned int ny;
	public:
		void Dimension(unsigned int nx0, unsigned int ny0) 
		{
			nx=nx0; ny=ny0;
			m_size=nx*ny;
		}
		void Dimension(unsigned int nx0, unsigned int ny0, T *v0) 
		{
			Dimension(nx0,ny0);
			m_v=v0;
			clear(e_allocated);
		}
		void Dimension(const array1<T> &A) 
		{
			ArrayExit("Operation not implemented");
		} 

		void Allocate(unsigned int nx0, unsigned int ny0, size_t align=0) 
		{
			Dimension(nx0,ny0);
			__checkActivate(2,align);
		}

		array2() : nx(0), ny(0) {}
		array2(unsigned int nx0, unsigned int ny0, size_t align=0) 
		{
			Allocate(nx0,ny0,align);
		}
		array2(unsigned int nx0, unsigned int ny0, T *v0) 
		{
			Dimension(nx0,ny0,v0);
		}

		unsigned int Nx() const {return nx;}
		unsigned int Ny() const {return ny;}

#ifndef __NOARRAY2OPT
		T *operator [] (int ix) const 
		{
			return m_v+ix*ny;
		}
#else
		array1<T> operator [] (int ix) const 
		{
			__check(ix,nx,2,1);
			return array1<T>(ny,m_v+ix*ny);
		}
#endif
		T& operator () (int ix, int iy) const 
		{
			__check(ix,nx,2,1);
			__check(iy,ny,2,2);
			return m_v[ix*ny+iy];
		}
		T& operator () (int i) const 
		{
			__check(i,m_size,2,0);
			return m_v[i];
		}
		T* operator () () const 
		{
			return m_v;
		}

		array2<T>& operator = (T a) {Load(a); return *this;}
		array2<T>& operator = (T *a) {Load(a); return *this;}
		array2<T>& operator = (const array2<T>& A) 
		{
			if (test(e_allocated))
			{
				__checkEqual(nx,A.Nx(),2,1);
				__checkEqual(ny,A.Ny(),2,2);
			}
			else
			{
				Dimension(A.Nx(), A.Ny(), A);
			}
			Load(A);
//			A.Purge();
			return *this;
		}

		array2<T>& operator += (const array2<T>& A) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) m_v[i] += A(i);
			return *this;
		}
		array2<T>& operator -= (const array2<T>& A) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) m_v[i] -= A(i);
			return *this;
		}

		array2<T>& operator *= (const array2<T>& A)
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) m_v[i] *= A(i);
			return *this;
		}

		array2<T>& operator += (T a) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) 
				m_v[i] += a;
			return *this;
		}
		array2<T>& operator -= (T a) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i ++) 
				m_v[i] -= a;
			return *this;
		}
		array2<T>& operator *= (T a) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) 
				m_v[i] *= a;
			return *this;
		}

		void Identity() 
		{
			Load((T) 0);
			__checkSize();
			unsigned int inc=ny+1;
			for(unsigned int i=0; i < m_size; i += inc) 
				m_v[i]=(T) 1;
		}

		void abs() 
		{
			__checkSize();
			for(unsigned int i = 0; i < m_size; i ++) 
				m_v[i] = (T)fabs(m_v[i]);
		}

		T Sum()
		{
			T tSum = 0;
			for (int ii = 0; ii < (int)m_size;ii++)
				tSum += m_v[ii];
			return tSum;
		}
		void Clip(unsigned int nBorder)
		{	// set left and right edges to zeros
			for (unsigned int i = 0; i < nBorder; i++)
			{
				for (unsigned int j = nBorder; j < ny; j++)
				{
					m_v[j+i*nx]	= 0;
					m_v[i*nx+(nx-i)] = 0;
				}
			}
			for (unsigned int j = 0; j < nBorder; j++)
			{ // set top and bottom rows to zeros
				for (unsigned int i = 0; i < nx; i++)
				{
					m_v[i+j*nx]	= 0;
					m_v[i+(ny-j-1)*nx] = 0;
				}
			}
		}

		array2<T> Crop(unsigned int nBorder)
		{ // get area inside border
			array2<T> ret(nx-(nBorder*2), ny-(nBorder*2));
			ret.Load((T) 0);
			for (unsigned int i = nBorder; i < nx-nBorder; i++)
			{
				for (unsigned int j = nBorder; j < ny - nBorder; j++)
				{
					ret.m_v[(j - nBorder)+(i - nBorder)*(nx - nBorder*2)] = m_v[j+i*nx];
				}
			}
			ret.Hold();
			return ret;
		}
	};

	template<class T> ostream& operator << (ostream& s, const array2<T>& A)
	{
		T *p=A();
		for(unsigned int i=0; i < A.Nx(); i++) 
		{
			for(unsigned int j=0; j < A.Ny(); j++) 
			{
				s << *(p++) << " ";
			}
			s << _newl;
		}
		s << flush;
		return s;
	}

	template<class T> istream& operator >> (istream& s, const array2<T>& A)
	{
		return A.Input(s);
	}

	template <class T> class array3 : public array1<T> 
	{
	protected:
		unsigned int nx;
		unsigned int ny;
		unsigned int nz;
		unsigned int nyz;
	public:
		void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0) 
		{
			nx=nx0; ny=ny0; nz=nz0; nyz=ny*nz;
			m_size=nx*nyz;
		}
		void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0, T *v0) 
		{
			Dimension(nx0,ny0,nz0);
			m_v=v0;
			clear(e_allocated);
		}

		void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			size_t align=0) 
		{
			Dimension(nx0,ny0,nz0);
			__checkActivate(3,align);
		}

		array3() : nx(0), ny(0), nz(0), nyz(0) {}
		array3(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			size_t align=0) 
		{
			Allocate(nx0,ny0,nz0,align);
		}
		array3(unsigned int nx0, unsigned int ny0, unsigned int nz0, T *v0) 
		{
			Dimension(nx0,ny0,nz0,v0);
		}

		unsigned int Nx() const {return nx;}
		unsigned int Ny() const {return ny;}
		unsigned int Nz() const {return nz;}

		array2<T> operator [] (int ix) const 
		{
			__check(ix,nx,3,1);
			return array2<T>(ny,nz,m_v+ix*nyz);
		}
		T& operator () (int ix, int iy, int iz) const 
		{
			__check(ix,nx,3,1);
			__check(iy,ny,3,2);
			__check(iz,nz,3,3);
			return m_v[ix*nyz+iy*nz+iz];
		}
		T& operator () (int i) const 
		{
			__check(i,m_size,3,0);
			return m_v[i];
		}
		T* operator () () const {return m_v;}

		array3<T>& operator = (T a) {Load(a); return *this;}
		array3<T>& operator = (T *a) {Load(a); return *this;}
		array3<T>& operator = (const array3<T>& A) 
		{
			__checkEqual(nx,A.Nx(),3,1);
			__checkEqual(ny,A.Ny(),3,2);
			__checkEqual(nz,A.Nz(),3,3);
			Load(A());
			A.Purge(); 
			return *this;
		}

		array3<T>& operator += (array3<T>& A) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) m_v[i] += A(i);
			return *this;
		}
		array3<T>& operator -= (array3<T>& A) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) m_v[i] -= A(i);
			return *this;
		}

		array3<T>& operator += (T a) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size;i ++) 
				m_v[i] += a;
			return *this;
		}
		array3<T>& operator -= (T a) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) 
				m_v[i] -= a;
			return *this;
		}
	};

	template<class T> ostream& operator << (ostream& s, const array3<T>& A)
	{
		T *p=A();
		for(unsigned int i=0; i < A.Nx(); i++) 
		{
			for(unsigned int j=0; j < A.Ny(); j++)
			{
				for(unsigned int k=0; k < A.Nz(); k++) 
				{
					s << *(p++) << " ";
				}
				s << _newl;
			}
			s << _newl;
		}
		s << flush;
		return s;
	}

	template<class T> istream& operator >> (istream& s, const array3<T>& A)
	{
		return A.Input(s);
	}

	template<class T> class array4 : public array1<T> 
	{
	protected:
		unsigned int nx;
		unsigned int ny;
		unsigned int nz;
		unsigned int nw;
		unsigned int nyz;
		unsigned int nzw;
		unsigned int nyzw;
	public:
		void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			unsigned int nw0) 
		{
			nx=nx0; ny=ny0; nz=nz0; nw=nw0; nzw=nz*nw; nyzw=ny*nzw;
			m_size=nx*nyzw;
		}
		void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			unsigned int nw0, T *v0) 
		{
			Dimension(nx0,ny0,nz0,nw0);
			m_v=v0;
			clear(e_allocated);
		}

		void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			unsigned int nw0, size_t align=0) 
		{
			Dimension(nx0,ny0,nz0,nw0);
			__checkActivate(4,align);
		}

		array4() : nx(0), ny(0), nz(0), nw(0), nyz(0), nzw(0), nyzw(0) {}
		array4(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			unsigned int nw0, size_t align=0) 
		{Allocate(nx0,ny0,nz0,nw0,align);}
			array4(unsigned int nx0, unsigned int ny0, unsigned int nz0,
				unsigned int nw0, T *v0) 
		{
			Dimension(nx0,ny0,nz0,nw0,v0);
		}

		unsigned int Nx() const {return nx;}
		unsigned int Ny() const {return ny;}
		unsigned int Nz() const {return ny;}
		unsigned int N4() const {return nw;}

		array3<T> operator [] (int ix) const 
		{
			__check(ix,nx,3,1);
			return array3<T>(ny,nz,nw,m_v+ix*nyzw);
		}
		T& operator () (int ix, int iy, int iz, int iw) const 
		{
			__check(ix,nx,4,1);
			__check(iy,ny,4,2);
			__check(iz,nz,4,3);
			__check(iw,nw,4,4);
			return m_v[ix*nyzw+iy*nzw+iz*nw+iw];
		}
		T& operator () (int i) const 
		{
			__check(i,m_size,4,0);
			return m_v[i];
		}
		T* operator () () const {return m_v;}

		array4<T>& operator = (T a) {Load(a); return *this;}
		array4<T>& operator = (T *a) {Load(a); return *this;}
		array4<T>& operator = (const array4<T>& A) 
		{
			__checkEqual(nx,A.Nx(),4,1);
			__checkEqual(ny,A.Ny(),4,2);
			__checkEqual(nz,A.Nz(),4,3);
			__checkEqual(nw,A.N4(),4,4);
			Load(A());
			A.Purge();
			return *this;
		}

		array4<T>& operator += (array4<T>& A) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) m_v[i] += A(i);
			return *this;
		}
		array4<T>& operator -= (array4<T>& A) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) m_v[i] -= A(i);
			return *this;
		}

		array4<T>& operator += (T a) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i ++) 
				m_v[i] += a;
			return *this;
		}
		array4<T>& operator -= (T a) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) 
				m_v[i] -= a;
			return *this;
		}
	};

	template<class T> ostream& operator << (ostream& s, const array4<T>& A)
	{
		T *p=A;
		for(unsigned int i=0; i < A.Nx(); i++) 
		{
			for(unsigned int j=0; j < A.Ny(); j++) 
			{
				for(unsigned int k=0; k < A.Nz(); k++) 
				{
					for(unsigned int l=0; l < A.N4(); l++) 
					{
						s << *(p++) << " ";
					}
					s << _newl;
				}
				s << _newl;
			}
			s << _newl;
		}	
		s << flush;
		return s;
	}

	template <class T> istream& operator >> (istream& s, const array4<T>& A)
	{
		return A.Input(s);
	}

	template <class T> class array5 : public array1<T> 
	{
	protected:
		unsigned int nx;
		unsigned int ny;
		unsigned int nz;
		unsigned int nw;
		unsigned int nv;
		unsigned int nwv;
		unsigned int nzwv;
		unsigned int nyzwv;
	public:
		void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			unsigned int nw0, unsigned int nv0) 
		{
			nx=nx0; ny=ny0; nz=nz0; nw=nw0; nv=nv0; nwv=nw*nv; nzwv=nz*nwv;
			nyzwv=ny*nzwv;
			m_size=nx*nyzwv;
		}
		void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			unsigned int nw0, unsigned int nv0, T *v0) 
		{
			Dimension(nx0,ny0,nz0,nw0,nv0);
			m_v=v0;
			clear(e_allocated);
		}

		void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			unsigned int nw0, unsigned int nv0, size_t align=0) 
		{
			Dimension(nx0,ny0,nz0,nw0,nv0);
			__checkActivate(5,align);
		}

		array5() : nx(0), ny(0), nz(0), nw(0), nv(0), nwv(0), nzwv(0), nyzwv(0) {}
		array5(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			unsigned int nw0, unsigned int nv0, size_t align=0) 
		{
			Allocate(nx0,ny0,nz0,nw0,nv0,align);
		}
		array5(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			unsigned int nw0, unsigned int nv0, T *v0) 
		{
			Dimension(nx0,ny0,nz0,nw0,nv0,nv0);
		}

		unsigned int Nx() const {return nx;}
		unsigned int Ny() const {return ny;}
		unsigned int Nz() const {return ny;}
		unsigned int N4() const {return nw;}
		unsigned int N5() const {return nv;}

		array4<T> operator [] (int ix) const 
		{
			__check(ix,nx,4,1);
			return array4<T>(ny,nz,nw,nv,m_v+ix*nyzwv);
		}
		T& operator () (int ix, int iy, int iz, int iw, int iv) const 
		{
			__check(ix,nx,5,1);
			__check(iy,ny,5,2);
			__check(iz,nz,5,3);
			__check(iw,nw,5,4);
			__check(iv,nv,5,5);
			return m_v[ix*nyzwv+iy*nzwv+iz*nwv+iw*nv+iv];
		}
		T& operator () (int i) const 
		{
			__check(i,m_size,5,0);
			return m_v[i];
		}
		T* operator () () const {return m_v;}

		array5<T>& operator = (T a) {Load(a); return *this;}
		array5<T>& operator = (T *a) {Load(a); return *this;}
		array5<T>& operator = (const array5<T>& A) 
		{
			__checkEqual(nx,A.Nx(),5,1);
			__checkEqual(ny,A.Ny(),5,2);
			__checkEqual(nz,A.Nz(),5,3);
			__checkEqual(nw,A.N4(),5,4);
			__checkEqual(nv,A.N5(),5,5);
			Load(A());
			A.Purge();
			return *this;
		}

		array5<T>& operator += (array5<T>& A) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) m_v[i] += A(i);
			return *this;
		}
		array5<T>& operator -= (array5<T>& A) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) m_v[i] -= A(i);
			return *this;
		}

		array5<T>& operator += (T a) 
		{
			__checkSize();
			for(unsigned int i=0; i < m_size; i++) 
				m_v[i] += a;
			return *this;
		}
		array5<T>& operator -= (T a) 
		{
			__checkSize();
//			unsigned int inc=nyzwv+nzwv+nwv+nv+1;
			for(unsigned int i=0; i < m_size; i++)
				m_v[i] -= a;
			return *this;
		}
	};

	template<class T> ostream& operator << (ostream& s, const array5<T>& A)
	{
		T *p=A;
		for(unsigned int i=0; i < A.Nx(); i++) 
		{
			for(unsigned int j=0; j < A.Ny(); j++) 
			{
				for(unsigned int k=0; k < A.Nz(); k++) 
				{
					for(unsigned int l=0; l < A.N4(); l++) 
					{
						for(unsigned int l=0; l < A.N5(); l++) 
						{
							s << *(p++) << " ";
						}
						s << _newl;
					}
					s << _newl;
				}
				s << _newl;
			}
			s << _newl;
		}	
		s << flush;
		return s;
	}

	template<class T> istream& operator >> (istream& s, const array5<T>& A)
	{
		return A.Input(s);
	}

#undef __check

#ifdef NDEBUG
#define __check(i,n,o,dim,m)
#else
#define __check(i,n,o,dim,m) Check(i-o,n,dim,m,o)
#endif

	template <class T> class Array1 : public array1<T> 
	{
	protected:
		T *voff; // Offset pointer to memory block
		int ox;
	public:
		void Offsets() 
		{
			voff=m_v-ox;
		}
		void Dimension(unsigned int nx0, int ox0=0) 
		{
			m_size = nx0;
			ox = ox0;
			Offsets();
		}
		void Dimension(unsigned int nx0, T *v0, int ox0=0) 
		{
			m_v = v0;
			Dimension(nx0,ox0);
			clear(e_allocated);
		}
		void Dimension(const Array1<T>& A) 
		{
			Dimension(A.m_size,A.m_v,A.ox); 
			set(A.test(e_temporary));
		}

		void Allocate(unsigned int nx0, int ox0=0, size_t align=0) 
		{
			Dimension(nx0,ox0);
			__checkActivate(1,align);
			Offsets();
		}

		void Reallocate(unsigned int nx0, int ox0=0, size_t align=0) 
		{
			Deallocate();
			Allocate(nx0,ox0,align);
		}

		Array1() : ox(0) {}
		Array1(unsigned int nx0, int ox0=0, size_t align=0) 
		{
			Allocate(nx0,ox0,align);
		}
		Array1(unsigned int nx0, T *v0, int ox0=0) 
		{
			Dimension(nx0,v0,ox0);
		}
		Array1(T *v0) 
		{
			Dimension(INT_MAX,v0);
		}

#ifdef NDEBUG
		typedef T *opt;
#else
		typedef Array1<T> opt;
#endif  

		T& operator [] (int ix) const {__check(ix,m_size,ox,1,1); return voff[ix];}
		T& operator () (int i) const  {__check(i,m_size,0,1,1); return m_v[i];}
		T* operator () () const {return voff;}
		operator T* () const {return voff;}

		Array1<T> operator + (int i) const {return Array1<T>(m_size-i,m_v+i,ox);}
		void Set(T *a) {m_v=a; Offsets(); clear(e_allocated);}

		Array1<T>& operator = (T a) {Load(a); return *this;}
		Array1<T>& operator = (const T *a) {Load(a); return *this;}
		Array1<T>& operator = (const Array1<T>& A) 
		{
			__checkEqual(m_size,A.Size(),1,1);
			__checkEqual(ox,A.Ox(),1,1);
			Load(A());
			A.Purge();
			return *this;
		}
		Array1<T>& operator = (const array1<T>& A) 
		{
			__checkEqual(m_size,A.Size(),1,1);
			__checkEqual(ox,0,1,1);
			Load(A());
			A.Purge();
			return *this;
		}

		int Ox() const {return ox;}
	};

	template <class T> class Array2 : public array2<T> 
	{
	protected:
		T *voff,*vtemp;
		int ox,oy;
	public:
		void Offsets() 
		{
			vtemp=m_v-ox*(int) this->ny;
			voff=vtemp-oy;
		}
		void Dimension(unsigned int nx0, unsigned int ny0, int ox0=0, int oy0=0) 
		{
			this->nx=nx0; this->ny=ny0;
			m_size=this->nx*this->ny;
			ox=ox0; oy=oy0;
			Offsets();
		}
		void Dimension(unsigned int nx0, unsigned int ny0, T *v0, int ox0=0,
			int oy0=0) 
		{
			m_v=v0;
			Dimension(nx0,ny0,ox0,oy0);
			clear(e_allocated);
		}

		void Allocate(unsigned int nx0, unsigned int ny0, int ox0=0, int oy0=0,
			size_t align=0) 
		{
			Dimension(nx0,ny0,ox0,oy0);
			__checkActivate(2,align);
			Offsets();
		}

		Array2() : ox(0), oy(0) {}
		Array2(unsigned int nx0, unsigned int ny0, int ox0=0, int oy0=0,
			size_t align=0) 
		{
			Allocate(nx0,ny0,ox0,oy0,align);
		}
		Array2(unsigned int nx0, unsigned int ny0, T *v0, int ox0=0, int oy0=0) 
		{
			Dimension(nx0,ny0,v0,ox0,oy0);
		}

#ifndef __NOARRAY2OPT
		T *operator [] (int ix) const 
		{
			return voff+ix*(int) this->ny;
		}
#else
		Array1<T> operator [] (int ix) const 
		{
			__check(ix,this->nx,ox,2,1);
			return Array1<T>(this->ny,vtemp+ix*(int) this->ny,oy);
		}
#endif

		T& operator () (int ix, int iy) const 
		{
			__check(ix,this->nx,ox,2,1);
			__check(iy,this->ny,oy,2,2);
			return voff[ix*(int) this->ny+iy];
		}
		T& operator () (int i) const 
		{
			__check(i,m_size,0,2,0);
			return m_v[i];
		}
		T* operator () () const {return voff;}
		void Set(T *a) {m_v=a; Offsets(); clear(e_allocated);}

		Array2<T>& operator = (T a) {Load(a); return *this;}
		Array2<T>& operator = (T *a) {Load(a); return *this;}
		Array2<T>& operator = (const Array2<T>& A) 
		{
			__checkEqual(this->nx,A.Nx(),2,1);
			__checkEqual(this->ny,A.Ny(),2,2);
			__checkEqual(ox,A.Ox(),2,1);
			__checkEqual(oy,A.Oy(),2,2);
			Load(A());
			A.Purge();
			return *this;
		}
		Array2<T>& operator = (const array2<T>& A) 
		{
			__checkEqual(this->nx,A.Nx(),2,1);
			__checkEqual(this->ny,A.Ny(),2,2);
			__checkEqual(ox,0,2,1);
			__checkEqual(oy,0,2,2);
			Load(A());
			A.Purge();
			return *this;
		}

		int Ox() const {return ox;}
		int Oy() const {return oy;}

	};

	template<class T> class Array3 : public array3<T> 
	{
	protected:
		T *voff,*vtemp;
		int ox,oy,oz;
	public:
		void Offsets() 
		{
			vtemp=m_v-ox*(int) this->nyz;
			voff=vtemp-oy*(int) this->nz-oz;
		}
		void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			int ox0=0, int oy0=0, int oz0=0) 
		{
			this->nx=nx0; this->ny=ny0; this->nz=nz0; this->nyz=this->ny*this->nz;
			m_size=this->nx*this->nyz;
			ox=ox0; oy=oy0; oz=oz0;
			Offsets();
		}
		void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			T *v0, int ox0=0, int oy0=0, int oz0=0) 
		{
			m_v=v0;
			Dimension(nx0,ny0,nz0,ox0,oy0,oz0);
			clear(e_allocated);
		}

		void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			int ox0=0, int oy0=0, int oz0=0, size_t align=0) 
		{
			Dimension(nx0,ny0,nz0,ox0,oy0,oz0);
			__checkActivate(3,align);
			Offsets();
		}

		Array3() : ox(0), oy(0), oz(0) {}
		Array3(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			int ox0=0, int oy0=0, int oz0=0, size_t align=0) 
		{
			Allocate(nx0,ny0,nz0,ox0,oy0,oz0,align);
		}
		Array3(unsigned int nx0, unsigned int ny0, unsigned int nz0, T *v0,
			int ox0=0, int oy0=0, int oz0=0) 
		{
			Dimension(nx0,ny0,nz0,v0,ox0,oy0,oz0);
		}

		Array2<T> operator [] (int ix) const 
		{
			__check(ix,this->nx,ox,3,1);
			return Array2<T>(this->ny,this->nz,vtemp+ix*(int) this->nyz,oy,oz);
		}
		T& operator () (int ix, int iy, int iz) const 
		{
			__check(ix,this->nx,ox,3,1);
			__check(iy,this->ny,oy,3,2);
			__check(iz,this->nz,oz,3,3);
			return voff[ix*(int) this->nyz+iy*(int) this->nz+iz];
		}
		T& operator () (int i) const 
		{
			__check(i,m_size,0,3,0);
			return m_v[i];
		}
		T* operator () () const {return voff;}
		void Set(T *a) {m_v=a; Offsets(); clear(e_allocated);}

		Array3<T>& operator = (T a) {Load(a); return *this;}
		Array3<T>& operator = (T *a) {Load(a); return *this;}
		Array3<T>& operator = (const Array3<T>& A) 
		{
			__checkEqual(this->nx,A.Nx(),3,1);
			__checkEqual(this->ny,A.Ny(),3,2);
			__checkEqual(this->nz,A.Nz(),3,3);
			__checkEqual(ox,A.Ox(),3,1);
			__checkEqual(oy,A.Oy(),3,2);
			__checkEqual(oz,A.Oz(),3,3);
			Load(A());
			A.Purge(); 
			return *this;
		}
		Array3<T>& operator = (const array3<T>& A) 
		{
			__checkEqual(this->nx,A.Nx(),3,1);
			__checkEqual(this->ny,A.Ny(),3,2);
			__checkEqual(this->nz,A.Nz(),3,3);
			__checkEqual(ox,0,3,1);
			__checkEqual(oy,0,3,2);
			__checkEqual(oz,0,3,3);
			Load(A());
			A.Purge(); 
			return *this;
		}

		int Ox() const {return ox;}
		int Oy() const {return oy;}
		int Oz() const {return oz;}

	};

	template <class T> class Array4 : public array4<T> 
	{
	protected:
		T *voff,*vtemp;
		int ox,oy,oz,ow;
	public:
		void Offsets() 
		{
			vtemp=m_v-ox*(int) this->nyzw;
			voff=vtemp-oy*(int) this->nzw-oz*(int) this->nw-ow;
		}
		void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			unsigned int nw0, 
			int ox0=0, int oy0=0, int oz0=0, int ow0=0) 
		{
			this->nx=nx0; this->ny=ny0; this->nz=nz0; this->nw=nw0;
			this->nzw=this->nz*this->nw; this->nyzw=this->ny*this->nzw;
			m_size=this->nx*this->nyzw;
			ox=ox0; oy=oy0; oz=oz0; ow=ow0;
			Offsets();
		}
		void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			unsigned int nw0, T *v0,
			int ox0=0, int oy0=0, int oz0=0, int ow0=0) 
		{
			m_v=v0;
			Dimension(nx0,ny0,nz0,nw0,ox0,oy0,oz0,ow0);
			clear(e_allocated);
		}

		void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			unsigned int nw0,
			int ox0=0, int oy0=0, int oz0=0, int ow0=0, size_t align=0) 
		{
			Dimension(nx0,ny0,nz0,nw0,ox0,oy0,oz0,ow0);
			__checkActivate(4,align); 
			Offsets();
		}

		Array4() : ox(0), oy(0), oz(0), ow(0) {}
		Array4(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			unsigned int nw0,
			int ox0=0, int oy0=0, int oz0=0, int ow0=0, size_t align=0) 
		{
			Allocate(nx0,ny0,nz0,nw0,ox0,oy0,oz0,ow0,align);
		}
		Array4(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			unsigned int nw0, T *v0,
			int ox0=0, int oy0=0, int oz0=0, int ow0=0) 
		{
			Dimension(nx0,ny0,nz0,nw0,v0,ox0,oy0,oz0,ow0);
		}

		Array3<T> operator [] (int ix) const 
		{
			__check(ix,this->nx,ox,3,1);
			return Array3<T>(this->ny,this->nz,this->nw,vtemp+ix*(int) this->nyzw,
				oy,oz,ow);
		}
		T& operator () (int ix, int iy, int iz, int iw) const 
		{
			__check(ix,this->nx,ox,4,1);
			__check(iy,this->ny,oy,4,2);
			__check(iz,this->nz,oz,4,3);
			__check(iw,this->nw,ow,4,4);
			return voff[ix*(int) this->nyzw+iy*(int) this->nzw+iz*(int) this->nw+iw];
		}
		T& operator () (int i) const 
		{
			__check(i,m_size,0,4,0);
			return m_v[i];
		}
		T* operator () () const {return voff;}
		void Set(T *a) {m_v=a; Offsets(); clear(e_allocated);}

		Array4<T>& operator = (T a) {Load(a); return *this;}
		Array4<T>& operator = (T *a) {Load(a); return *this;}

		Array4<T>& operator = (const Array4<T>& A) 
		{
			__checkEqual(this->nx,A.Nx(),4,1);
			__checkEqual(this->ny,A.Ny(),4,2);
			__checkEqual(this->nz,A.Nz(),4,3);
			__checkEqual(this->nw,A.N4(),4,4);
			__checkEqual(ox,A.Ox(),4,1);
			__checkEqual(oy,A.Oy(),4,2);
			__checkEqual(oz,A.Oz(),4,3);
			__checkEqual(ow,A.O4(),4,4);
			Load(A());
			A.Purge();
			return *this;
		}
		Array4<T>& operator = (const array4<T>& A) 
		{
			__checkEqual(this->nx,A.Nx(),4,1);
			__checkEqual(this->ny,A.Ny(),4,2);
			__checkEqual(this->nz,A.Nz(),4,3);
			__checkEqual(this->nw,A.N4(),4,4);
			__checkEqual(this->nx,A.Nx(),4,1);
			__checkEqual(this->ny,A.Nx(),4,2);
			__checkEqual(this->nz,A.Nx(),4,3);
			__checkEqual(this->nw,A.Nx(),4,4);
			__checkEqual(ox,0,4,1);
			__checkEqual(oy,0,4,2);
			__checkEqual(oz,0,4,3);
			__checkEqual(ow,0,4,4);
			Load(A());
			A.Purge();
			return *this;
		}

		int Ox() const {return ox;}
		int Oy() const {return oy;}
		int Oz() const {return oz;}
		int O4() const {return ow;}
	};

	template<class T> class Array5 : public array5<T> 
	{
	protected:
		T *voff,*vtemp;
		int ox,oy,oz,ow,ov;
	public:
		void Offsets() 
		{
			vtemp=m_v-ox*(int) this->nyzwv;
			voff=vtemp-oy*(int) this->nzwv-oz*(int) this->nwv-ow*(int) this->nv-ov;
		}
		void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			unsigned int nw0,  unsigned int nv0,
			int ox0=0, int oy0=0, int oz0=0, int ow0=0, int ov0=0) 
		{
			this->nx=nx0; this->ny=ny0; this->nz=nz0; this->nw=nw0; this->nv=nv0;
			this->nwv=this->nw*this->nv; this->nzwv=this->nz*this->nwv;
			this->nyzwv=this->ny*this->nzwv;
			m_size=this->nx*this->nyzwv;
			ox=ox0; oy=oy0; oz=oz0; ow=ow0; ov=ov0;
			Offsets();
		}
		void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			unsigned int nw0, unsigned int nv0, T *v0,
			int ox0=0, int oy0=0, int oz0=0, int ow0=0, int ov0=0) 
		{
			m_v=v0;
			Dimension(nx0,ny0,nz0,nw0,nv0,ox0,oy0,oz0,ow0,ov0);
			clear(e_allocated);
		}

		void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			unsigned int nw0, unsigned int nv0,
			int ox0=0, int oy0=0, int oz0=0, int ow0=0, int ov0=0,
			size_t align=0) 
		{
			Dimension(nx0,ny0,nz0,nw0,nv0,ox0,oy0,oz0,ow0,ov0);
			__checkActivate(5,align); 
			Offsets();
		}

		Array5() : ox(0), oy(0), oz(0), ow(0), ov(0) {}
		Array5(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			unsigned int nw0, unsigned int nv0, int ox0=0, int oy0=0,
			int oz0=0, int ow0=0, int ov0=0, size_t align=0) 
		{
			Allocate(nx0,ny0,nz0,nw0,nv0,ox0,oy0,oz0,ow0,ov0,align);
		}
		Array5(unsigned int nx0, unsigned int ny0, unsigned int nz0,
			unsigned int nw0, unsigned int nv0, T *v0,
			int ox0=0, int oy0=0, int oz0=0, int ow0=0, int ov0=0) 
		{
			Dimension(nx0,ny0,nz0,nw0,nv0,v0,ox0,oy0,oz0,ow0,ov0);
		}

		Array4<T> operator [] (int ix) const 
		{
			__check(ix,this->nx,ox,4,1);
			return Array4<T>(this->ny,this->nz,this->nw,this->nv,
				vtemp+ix*(int) this->nyzwv,oy,oz,ow,ov);
		}
		T& operator () (int ix, int iy, int iz, int iw, int iv) const 
		{
			__check(ix,this->nx,ox,5,1);
			__check(iy,this->ny,oy,5,2);
			__check(iz,this->nz,oz,5,3);
			__check(iw,this->nw,ow,5,4);
			__check(iv,this->nv,ov,5,5);
			return voff[ix*(int) this->nyzwv+iy*(int) this->nzwv+iz*(int) this->nwv
				+iw*(int) this->nv+iv];
		}
		T& operator () (int i) const 
		{
			__check(i,m_size,0,5,0);
			return m_v[i];
		}
		T* operator () () const {return voff;}
		void Set(T *a) {m_v=a; Offsets(); clear(e_allocated);}

		Array5<T>& operator = (T a) {Load(a); return *this;}
		Array5<T>& operator = (T *a) {Load(a); return *this;}

		Array5<T>& operator = (const Array5<T>& A) 
		{
			__checkEqual(this->nx,A.Nx(),5,1);
			__checkEqual(this->ny,A.Ny(),5,2);
			__checkEqual(this->nz,A.Nz(),5,3);
			__checkEqual(this->nw,A.N4(),5,4);
			__checkEqual(this->nv,A.N5(),5,5);
			__checkEqual(ox,A.Ox(),5,1);
			__checkEqual(oy,A.Oy(),5,2);
			__checkEqual(oz,A.Oz(),5,3);
			__checkEqual(ow,A.O4(),5,4);
			__checkEqual(ov,A.O5(),5,5);
			Load(A());
			A.Purge();
			return *this;
		}
		Array5<T>& operator = (const array5<T>& A) 
		{
			__checkEqual(this->nx,A.Nx(),5,1);
			__checkEqual(this->ny,A.Ny(),5,2);
			__checkEqual(this->nz,A.Nz(),5,3);
			__checkEqual(this->nw,A.N4(),5,4);
			__checkEqual(this->nv,A.N5(),5,5);
			__checkEqual(ox,0,5,1);
			__checkEqual(oy,0,5,2);
			__checkEqual(oz,0,5,3);
			__checkEqual(ow,0,5,4);
			__checkEqual(ov,0,5,5);
			Load(A());
			A.Purge();
			return *this;
		}
		int Ox() const {return ox;}
		int Oy() const {return oy;}
		int Oz() const {return oz;}
		int O4() const {return ow;}
		int O5() const {return ov;}
	};

	template<class T> inline bool Active(array1<T>& A)
	{
		return A.Size();
	}

	template<class T> inline bool Active(T *A)
	{
		return A;
	}

	template<class T> inline void Set(T *&A, T *v)
	{
		A=v;
	}

	template<class T> inline void Set(array1<T>& A, T *v)
	{
		A.Set(v);
	}

	template<class T> inline void Set(array1<T>& A, const array1<T>& B)
	{
		A.Set(B());
	}

	template<class T> inline void Set(Array1<T>& A, T *v)
	{
		A.Set(v);
	}

	template<class T> inline void Set(Array1<T>& A, const array1<T>& B)
	{
		A.Set(B());
	}

	template<class T> inline void Null(T *&A)
	{
		A=NULL;
	}

	template<class T> inline void Null(array1<T>& A)
	{
		A.Dimension(0);
	}

	template<class T> inline void Dimension(T *&A, unsigned int)
	{
	}

	template<class T> inline void Dimension(array1<T> &A, unsigned int n)
	{
		A.Dimension(n);
	}

	template<class T> inline void Dimension(T *&A, unsigned int, T *v)
	{
		A=v;
	}

	template<class T> inline void Dimension(array1<T>& A, unsigned int n, T *v)
	{
		A.Dimension(n,v);
	}

	template<class T> inline void Dimension(Array1<T>& A, unsigned int n, T *v)
	{
		A.Dimension(n,v);
	}

	template<class T> inline void Dimension(T *&A, T *v)
	{
		A=v;
	}

	template<class T> inline void Dimension(array1<T>& A, const array1<T>& B)
	{
		A.Dimension(B);
	}

	template<class T> inline void Dimension(Array1<T>& A, const Array1<T>& B)
	{
		A.Dimension(B);
	}

	template<class T> inline void Dimension(Array1<T>& A, const array1<T>& B)
	{
		A.Dimension(B);
	}

	template<class T> inline void Dimension(array1<T>& A, unsigned int n, const array1<T>& B)
	{
		A.Dimension(n,B);
	}

	template<class T> inline void Dimension(Array1<T>& A, unsigned int n, const array1<T>& B, int o)
	{
		A.Dimension(n,B,o);
	}

	template<class T> inline void Dimension(Array1<T>& A, unsigned int n, T *v, int o)
	{
		A.Dimension(n,v,o);
	}

	template<class T> inline void Dimension(T *&A, unsigned int, T *v, int o)
	{
		A=v-o;
	}

	template<class T> inline void Allocate(T *&A, unsigned int n, size_t align=0)
	{
		if(align) newAlign(A,n,align);
		else A=new T[n];
	}

	template<class T> inline void Allocate(array1<T>& A, unsigned int n, size_t align=0)
	{  
		A.Allocate(n,align);
	}

	template<class T> inline void Allocate(Array1<T>& A, unsigned int n, size_t align=0)
	{  
		A.Allocate(n,align);
	}

	template<class T> inline void Allocate(T *&A, unsigned int n, int o, size_t align=0)
	{
		Allocate(A,n,align);
		A -= o;
	}

	template<class T> inline void Allocate(Array1<T>& A, unsigned int n, int o, size_t align=0)
	{  
		A.Allocate(n,o,align);
	}

	template<class T> inline void Deallocate(T *A)
	{
		if(A) delete [] A;
	}

	template<class T> inline void Deallocate(array1<T>& A)
	{  
		A.Deallocate();
	}

	template<class T> inline void Deallocate(Array1<T>& A)
	{  
		A.Deallocate();
	}

	template<class T> inline void Deallocate(T *A, int o)
	{
		if(A) delete [] (A+o);
	}

	template<class T> inline void Deallocate(Array1<T>& A, int)
	{  
		A.Deallocate();
	}

	template<class T> inline void Reallocate(T *&A, unsigned int n, size_t align=0)
	{  
		if(A) delete [] A;
		Allocate(A,n,align);
	}

	template<class T> inline void Reallocate(array1<T>& A, unsigned int n)
	{  
		A.Reallocate(n);
	}

	template<class T> inline void Reallocate(Array1<T>& A, unsigned int n)
	{  
		A.Reallocate(n);
	}

	template<class T>
		inline void Reallocate(T *&A, unsigned int n, int o, size_t align=0)
	{  
		if(A) delete [] A;
		Allocate(A,n,align);
		A -= o;
	}

	template<class T>
		inline void Reallocate(Array1<T>& A, unsigned int n, int o, size_t align=0)
	{  
		A.Reallocate(n,o,align);
	}

	template<class T>
		inline void CheckReallocate(T& A, unsigned int n, unsigned int& old,
		size_t align=0)
	{
		if(n > old) {A.Reallocate(n,align); old=n;}
	}

	template<class T>
		inline void CheckReallocate(T& A, unsigned int n, int o, unsigned int& old,
		size_t align=0)
	{
		if(n > old) {A.Reallocate(n,o,align); old=n;}
	}
} // namespace Array

#undef __check
#undef __checkSize
#undef __checkActivate

#endif
