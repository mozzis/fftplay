#include "Array.h"
#include "fftw++.h"

// Compile with g++ example3.cc fftw++.cc -lfftw3

using std::cout;
using Array::array3;

int main()
{
  unsigned int nx=4, ny=5, nz=6;
  unsigned int nzp=nz/2+1;
  size_t align=sizeof(Complex);
  
  array3<double> f(nx,ny,nz,align);
  array3<Complex> g(nx,ny,nzp,align);
  
  rcfft3d Forward3(f,g);
  crfft3d Backward3(g,f);
  
  for(unsigned int i=0; i < nx; i++) 
    for(unsigned int j=0; j < ny; j++) 
      for(unsigned int k=0; k < nz; k++) 
      f(i,j,k)=i+j+k;
	
  cout << f << endl;
  
  Forward3.fft(f,g);
  
  cout << g << endl;
  
  Backward3.fftNormalized(g,f);
  
  cout << f << endl;
}
