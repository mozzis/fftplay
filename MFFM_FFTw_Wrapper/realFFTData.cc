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
#include "realFFTData.H"

#include <limits>
#include <iostream>
#include <math.h>
#include <crtdbg.h>
#include <stdexcept>

realFFTData::realFFTData(int sz)
{
	deleteInOutMemory = 1;
	_RPT0(_CRT_WARN, "realFFTData init:\n");
	size=sz;
	in = out = power_spectrum = NULL;
	// powerDeriv = NULL;

	in = new fftw_real[size];
	out = new fftw_real[size];
	power_spectrum = new fftw_real[size/2+1];
	if (!in || !out || !power_spectrum)
	{
		std::cerr << "Could not allocate enough mem for a realFFT"<<std::endl;
		if (in) 
			delete [] in;
		if (out) 
			delete [] out;
		if (power_spectrum) 
			delete [] power_spectrum;
		throw std::runtime_error("Could not allocate enough mem for a realFFT in realFFTData");
	}
	totalPower = 0.0;
}

realFFTData::realFFTData(int sz, fftw_real *inp, fftw_real *outp)
{
	deleteInOutMemory = 0;
	_RPT0(_CRT_WARN, "realFFTData init:\n");
	size = sz;
	if (!inp || !outp)
	{
		std::cerr<<"realFFTData::realFFTData : input or output array doesn't exist"<<std::endl;
		throw std::runtime_error("input or output array not allocated in realFFTData");
	}
	in = inp;
	out = outp;
	power_spectrum = NULL;
	//powerDeriv = NULL;

	power_spectrum = new fftw_real[size/2+1];
	if (!power_spectrum)
	{
		std::cerr << "Could not allocate enough mem for a realFFT"<<std::endl;
		throw std::runtime_error("Out of memory in realFFTData");
	}
	totalPower = 0.0;
}

realFFTData::~realFFTData()
{
	if (power_spectrum) 
		delete [] power_spectrum; 
	power_spectrum=NULL;
	//if (powerDeriv) delete [] powerDeriv; powerDeriv=NULL;
	//  std::cout<<"realFFTData::~realFFTData"<<std::endl;
	if (deleteInOutMemory)
	{
		if (in) 
			delete [] in; 
		in=NULL;
		if (out) 
			delete [] out; 
		out=NULL;
	}
	_RPT0(_CRT_WARN, "realFFTData::~realFFTData exit\n");
}

fftw_real realFFTData::findMaxIn()
{
	fftw_real max = -DBL_MAX;
	for (int i = 0; i < getSize(); i++)
		if (in[i ] > max)
			max = in[i];
	return max;
}

void realFFTData::findMaxMinPowerBins(void)
{
	double min = DBL_MAX;
	double max= -DBL_MAX;
	for (int i=0; i<getHalfSize(); i++)
	{
		if (power_spectrum[i] > max)
			max = power_spectrum[maxPowerBin = i];
		if (power_spectrum[i] < min)
			min = power_spectrum[minPowerBin = i];
	}
	_RPT4(_CRT_WARN, "min bin %f\t%d\tmax bin %f\t%d", minPowerBin, min, maxPowerBin, max);
}


int realFFTData::limitHalfPowerSpec(double lim)
{
	double max = 0.0;
	int bin = 0;
	for (int i = 0; i < getHalfSize(); i++)
		if (power_spectrum[i] > max)
			max = power_spectrum[bin = i];
	for (int i=0; i<getHalfSize(); i++)
		power_spectrum[i] /= (max/lim);
	return bin;
}


int realFFTData::compPowerSpec()
{
	totalPower = 0.0;
	double min =  DBL_MAX;
	double max = -DBL_MAX;
	power_spectrum[maxPowerBin = 0] = out[0] * out[0];  /* DC component */

	for (int k = 1; k < (getSize() + 1)/2; ++k) // (k < N/2 rounded up)
	{  
		if ((power_spectrum[k] = out[k]*out[k] + 
			out[getSize()-k] * out[getSize()-k])>max)
		{
			max=power_spectrum[maxPowerBin = k];
		}
		if (power_spectrum[k] < min) 
			min=power_spectrum[minPowerBin = k];
		totalPower += power_spectrum[k];
	}
	if (getSize() % 2 == 0)
	{ /* N is even */
		power_spectrum[getSize()/2] = out[getSize()/2]*out[getSize()/2];  /* Nyquist freq. */
		if (power_spectrum[getSize()/2] > max)
			max = power_spectrum[maxPowerBin = getSize()/2];      
		if (power_spectrum[getSize()/2]<min)
			min = power_spectrum[minPowerBin = getSize()/2];
		totalPower += power_spectrum[getSize()/2];
	}
	return maxPowerBin;
}

int realFFTData::sqrtPowerSpec()
{
	double max = -DBL_MAX;
	for (int k = 0; k < (getSize() + 1)/2; ++k)  /* (k < N/2 rounded up) */
		if ((power_spectrum[k] = sqrt(power_spectrum[k])) > max)
			max = power_spectrum[maxPowerBin = k];
	return maxPowerBin;
}

/*
int realFFTData::powerSpecDeriv()
{
	if (!powerDeriv){ // create memory if it doesn't exist
		powerDeriv = new fftw_real[size/2+1];
		if (!powerDeriv){
			std::cerr << "Could not allocate enough mem for a powerSpectrum deriv"<<std::endl;
			exit(-1);
		}
	}

	int pos=0;
	double max = 0.0;
	powerDeriv[0] = 0.0;  // DC component
	for (int k = 1; k < (getSize()+1)/2; ++k){  // (k < N/2 rounded up)
		if (fabs(powerDeriv[k] = power_spectrum[k]-power_spectrum[k-1]) > max){
			max = fabs(powerDeriv[k] = power_spectrum[k]-power_spectrum[k-1]);
			pos = k;
		}
	}
	if (getSize() % 2 == 0) // N is even
		if (fabs(powerDeriv[getSize()/2] = power_spectrum[getSize()/2]-power_spectrum[getSize()/2-1]) > max){
			max = fabs(powerDeriv[getSize()/2] = power_spectrum[getSize()/2]-power_spectrum[getSize()/2-1]);
			pos = getSize()/2;
		}
		return pos;
}
*/

void  realFFTData::zeroFFTData(void)
{
	//cout<<"here"<<std::endl;
	for (int i=0;i < getSize();i++)
		out[i] = 0.0;
}
