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
#include "real2DFFT.H"

#include <iostream>
using namespace std;

#include <math.h>
//#include <fstream>

real2DFFTData::
real2DFFTData(int r, int c){
  //  std::cout <<"real2DFFTData init:"<<this<<std::endl;
  x=r; y=c;
  mem=NULL;
  out=NULL;
  in=power=NULL;
  timeXSum=xSum=ySum=realXSum=imagXSum=NULL;
  std::cout<<"power size="<<x*(c/2+1)<<std::endl;

  if (!(mem=new fftw_real[x*y+x*(c/2+1)+2*x*(c/2+1)+x+(c/2+1)])){
    std::cout<<"real2DFFTData: malloc 2 fail"<<std::endl;
    memDeInit();
  } else {
    in=&mem[0];
    power=&mem[x*y];
    out=(fftw_complex*)(power+x*(c/2+1));
    xSum=(fftw_real*)out+2*x*(c/2+1);
    ySum=xSum+x;
  }

  if (!(timeXSum=new fftw_real[x])){
    std::cout<<"real2DFFTData: malloc 3 fail"<<std::endl;
    memDeInit();
  }
  if (!(realXSum=new fftw_real[x])){
    std::cout<<"real2DFFTData: malloc 3a fail"<<std::endl;
    memDeInit();
  }
  if (!(imagXSum=new fftw_real[x])){
    std::cout<<"real2DFFTData: malloc 3b fail"<<std::endl;
    memDeInit();
  }

  maxPower=minPower=totalPower=0.0;
  maxXSumIndex=maxYSumIndex=0;
}

real2DFFTData::
~real2DFFTData(){
  memDeInit();
}

void real2DFFTData::
memDeInit(void){
  if (mem) delete [] mem;
  mem=in=power=NULL;
  out=NULL;
  if (timeXSum) delete [] timeXSum;
  timeXSum=NULL;
  if (realXSum) delete [] realXSum;
  realXSum=NULL;
  if (imagXSum) delete [] imagXSum;
  imagXSum=NULL;
  std::cout<<"real2DFFTData: DeInit Out"<<std::endl;
}

void real2DFFTData::
reScale(void){
  double factor=1.0/((double)x*(double)y);
  for (int i=0;i<x*(y/2+1);i++){
    c_re(out[i])*=factor; c_im(out[i])*=factor;
  }
}

void real2DFFTData::
compPowerSpec(){
  maxPower=totalPower=0.0; minPower=9999999.9;
  double temp;
  int maxIndexX, maxIndexY;
  int minIndexX, minIndexY;
  for (int i=0; i < x; i++){  // The x dimension
    int index=i*(y/2+1);
    for (int j=0;j<y/2+1;j++){ // The y dimension
      temp=c_re(out[index+j])*c_re(out[index+j])+c_im(out[index+j])*c_im(out[index+j]);
      totalPower += (power[index+j]=temp);
      //      std::cout<<temp<<'\t'<<maxPower<<'\t'<<minPower<<std::endl;
      //std::cout<<temp<<std::endl;
      if (temp>maxPower){
	  maxPower=temp; maxIndexX=i; maxIndexY=j;
      }
      if (temp<minPower){
	minPower=temp; minIndexX=i; minIndexY=j;
      }
    }
  }
  //std::cout <<"real2DFFTData: compPowerSpec: max indexes : (x, y) "<<maxIndexX<<'\t'<<maxIndexY<<std::endl;
  //std::cout <<"real2DFFTData: compPowerSpec: min indexes : (x, y) "<<minIndexX<<'\t'<<minIndexY<<std::endl;
}

#include <fstream>
void real2DFFTData::
compLogPowerSpec(){
  int maxIndexX, maxIndexY;
  int minIndexX, minIndexY;
  totalPower=0.0;
  maxPower=-999999999.9, minPower=9999999.9;
  double temp;
  std::ofstream ofs("2Dfft.dat");
  for (int i=0; i < x; i++){  // The x dimension
    int index=i*(y/2+1);
    for (int j=0;j<y/2+1;j++){ // The y dimension
      temp=c_re(out[index+j])*c_re(out[index+j])+c_im(out[index+j])*c_im(out[index+j]);
      temp=10.0*log10(temp);
      ofs<<temp<<'\t';
      totalPower += (power[index+j]=temp);
      //      std::cout<<temp<<'\t'<<maxPower<<'\t'<<minPower<<std::endl;
      //std::cout<<temp<<std::endl;
      if (temp>maxPower){
	maxPower=temp; maxIndexX=i; maxIndexY=j;
      }
      if (temp<minPower){
	minPower=temp; minIndexX=i; minIndexY=j;
      }
    }
    ofs <<std::endl;
  }
  ofs.close();
  //std::cout <<"real2DFFTData: compPowerSpec: max indexes : (x, y) "<<maxIndexX<<'\t'<<maxIndexY<<std::endl;
  //std::cout <<"real2DFFTData: compPowerSpec: min indexes : (x, y) "<<minIndexX<<'\t'<<minIndexY<<std::endl;
}

void real2DFFTData::
findYSum(int start, int stop){
  memset(ySum, 0, (y/2+1)*sizeof(fftw_real));
  for (int i=start; i < stop; i++){  // The x dimension
    int index=i*(y/2+1);
    for (int j=0;j<y/2+1;j++){ // The y dimension
      ySum[j]+=power[index+j];
    }
  }
  for (int j=0;j<y/2+1;j++){
    ySum[j]/=x;
  }
}

void real2DFFTData::
findYMax(void){
  ySumMin=99999999999.9;
  ySumMax=-99999999999.9;
  for (int j=0;j<y/2+1;j++){
    if (ySum[j]< ySumMin) ySumMin=ySum[j];
    if (ySum[j]> ySumMax && j>2){
      ySumMax=ySum[j];
      maxYSumIndex=j;
    }
  }
}

void real2DFFTData::
timeSpecAverage(){
  memset(timeXSum, 0, x*sizeof(fftw_real));
  for (int i=0; i < x; i++){  // The x dimension
    int index=i*y;
    for (int j=0;j<y;j++){ // The y dimension
      timeXSum[i]+=in[index+j];
      //      std::cout<<in[index+j]<<std::endl;
    }
    timeXSum[i]/=y;
  }
}

void real2DFFTData::
complexSpecAverage(){
  memset(realXSum, 0, x*sizeof(fftw_real));
  memset(imagXSum, 0, x*sizeof(fftw_real));
  for (int i=0; i < x; i++){  // The x dimension
    int index=i*(y/2+1);
    for (int j=0;j<y/2+1;j++){ // The y dimension
      realXSum[i]+=c_re(out[index+j]);
      imagXSum[i]+=c_im(out[index+j]);
    }
    realXSum[i]/=y/2+1;
    imagXSum[i]/=y/2+1;
  }
}

void real2DFFTData::
powerSpecAverage(){
  memset(xSum, 0, x*sizeof(fftw_real));
  memset(ySum, 0, (y/2+1)*sizeof(fftw_real));
  xSumMin=ySumMin=99999999999.9;
  xSumMax=ySumMax=-99999999999.9;

  for (int i=0; i < x; i++){  // The x dimension
    int index=i*(y/2+1);
    for (int j=0;j<y/2+1;j++){ // The y dimension
      //      if (i<x/7 || i>x*6/7)
      //    std::cout<<power[index+j]<<'\t';
	ySum[j]+=power[index+j];
      xSum[i]+=power[index+j];
    }
    xSum[i]/=y/2+1;
    if (xSum[i]< xSumMin) xSumMin=xSum[i];
    if (xSum[i]> xSumMax){
      xSumMax=xSum[i];
      maxXSumIndex=i;
    }
  }
  for (int j=0;j<y/2+1;j++){
    ySum[j]/=x;
    if (ySum[j]< ySumMin) ySumMin=ySum[j];
    if (ySum[j]> ySumMax){
      ySumMax=ySum[j];
      maxYSumIndex=j;
    }
  }
}

real2DFFT::
real2DFFT(real2DFFTData *d) {
  //std::cout <<"realFFT init:"<<this<<std::endl;
  data=d;
  // std::cout <<data->getXSize() << '\t'<<data->getYSize()<<std::endl;
  fwdPlan = fftw_plan_dft_r2c_2d(data->getXSize(), data->getYSize(), data->in, data->out, PLANTYPE);
  invPlan = fftw_plan_dft_c2r_2d(data->getXSize(), data->getYSize(), data->out, data->in, PLANTYPE);
}

real2DFFT::
~real2DFFT(){
  //  std::cout <<"realFFT DeInit:"<<this<<std::endl;
  fftw_destroy_plan(fwdPlan);
  fftw_destroy_plan(invPlan);
  std::cout <<"realFFT DeInit done"<<std::endl;
}

void real2DFFT::
fwdTransform(){
  if (!data)
    std::cerr<<"real2DFFT::fwdTransform : data not present"<<std::endl;
  else
    fftw_execute(fwdPlan);
}

void real2DFFT::
invTransform(){
  if (!data)
    std::cerr<<"real2DFFT::invTransform : data not present"<<std::endl;
  else
    fftw_execute(invPlan);
}
