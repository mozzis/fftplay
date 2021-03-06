#include "Array.h"
#include "fftw++.h"

// Compilation: g++ example2r.cc fftw++.cc -lfftw3

using std::cout;
using Array::array1;

int main()
{
  unsigned int nx=4, ny=5;
  unsigned int nyp=ny/2+1;
  size_t align=sizeof(Complex);
  
  array2<double> f(nx,ny,align);
  array2<Complex> g(nx,nyp,align);
  
  rcfft2d Forward(f,g);
  crfft2d Backward(g,f);
  
  for(unsigned int i=0; i < nx; i++) 
    for(unsigned int j=0; j < ny; j++) 
      f(i,j)=i+j;
	
  cout << f << endl;

  Forward.fft(f,g);
  
  cout << g << endl;
  
  Backward.fftNormalized(g,f);
  
  cout << f << endl;
}
