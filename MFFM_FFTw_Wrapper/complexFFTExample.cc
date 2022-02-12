/* Copyright 2001,2002 Matt Flax <flatmax@ieee.org>
   This file is part of the MFFM FFTw Wrapper library.

   MFFM MFFM FFTw Wrapper library is free software; you can 
   redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   MFFM FFTw Wrapper library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You have received a copy of the GNU General Public License
   along with the MFFM FFTw Wrapper library
*/
#include <fstream>
#include <iomanip>

#include "complexFFT.H"
#include <iomanip.h>

int main (){
  int count =8;
  complexFFTData fftData(count);
  complexFFT fft(&fftData);

  for (int i=0;i<count;i++){
    c_re(fftData.in[i])=(double)i;
    //    fftData.in[i].im=(double)-i;
    c_im(fftData.in[i])=(double)i+5.0;
  }

  fftw_real *temp=&c_re(fftData.in[0]);
  for (int i=0; i<count; i++)
    cout << temp[i]<<endl;
  

  // forward transform :
  fft.fwdTransform();
  // inverse transform :
  fft.invTransform();

  //  for (int i=0; i<count; i++)
  //  cout << fftData.in[i].re<<' '<<fftData.in[i].im<<endl;

  // Find the power spectrum ...
  fftData.compPowerSpec();
  for (int i=0; i<count; i++)
    cout << fftData.power_spectrum[i]<<endl;
}
