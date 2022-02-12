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
#include <string.h>
#include "real2DFFT.H"
#include <iomanip.h>
using namespace std;

int main(void){
  int x=8, y=8;
  real2DFFTData *fftData = new real2DFFTData(x,y);
  real2DFFT *fft= new real2DFFT(fftData);

  // clear the data
  fftData->clearInput();
  fftData->clearOutput();

  int temp=x/2, temp2=y/2;
  for (int j=0;j<x;j++)
    fftData->in[temp2+j*x]=10000.0;
  for (int j=0;j<y;j++)
    fftData->in[temp*y+j]=10000.0;
  //  fftData->in[temp*y+(y-1)/2]=20000.0;

  for (int i=0;i<fftData->getXSize();i++){
    for (int j=0;j<fftData->getYSize();j++)
      cout<<fftData->in[i*x+j]<<'\t';
    cout<<endl;
  }
  cout<<'\n'<<endl;
  fft->fwdTransform();
  fftData->reScale();
  fftData->compPowerSpec();
  fft->invTransform();

  for (int i=0;i<fftData->getXSize();i++){
    for (int j=0;j<fftData->getYSize();j++)
      cout<<fftData->in[i*x+j]<<'\t';
    cout<<endl;
  }
  cout<<'\n'<<endl;
  /*  for (int i=0;i<fftData.getXSize();i++){
    for (int j=0;j<fftData.getYHalfSize();j++)
      cout<<fftData.out[i][j].im<<'\t';
    cout<<endl;
    }*/
  for (int i=0;i<x;i++){
    for (int j=0;j<y/2+1;j++)
      cout<<fftData->power[i*(y/2+1)+j]<<'\t';
    cout<<endl;
  }
  delete fftData;
  delete fft;
}
