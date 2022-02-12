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
#include "complex2DFFT.H"

#include <fstream>
#include <iostream>
#include <math.h>
#include <crtdbg.h>

using namespace std;

complex2DFFTData::complex2DFFTData(int r, int c, fftw_complex *pMem /* = 0 */)
{
	x = r; y = c;
	out = NULL;
	in = NULL;
	power_spectrum =  NULL;

	in  = (fftw_complex *)fftw_malloc(x*y*sizeof(fftw_complex));
	out = (fftw_complex *)fftw_malloc(x*y*sizeof(fftw_complex));
	power_spectrum = (fftw_real *)fftw_malloc(x*y*sizeof(fftw_real));
	_ASSERTE(in != 0 && out != 0 && power_spectrum != 0);
	maxPower = minPower = totalPower = 0;
}

complex2DFFTData::~complex2DFFTData()
{
	memDeInit();
}

void complex2DFFTData::memDeInit(void)
{
	if (in)
		fftw_free(in);
	if (out)
		fftw_free(out);
	if (power_spectrum)
		fftw_free(power_spectrum);
	in = NULL;
	out = NULL;
	power_spectrum = NULL;
	_RPT0(_CRT_WARN, "complex2DFFTData: DeInit Out\n");
}

void complex2DFFTData::reScale(void)
{
	double factor=1.0/((double)x*(double)y);
	for (int i = 0;i < x * y;i++)
	{
		c_re(out[i]) *=  factor; 
		c_im(out[i]) *=  factor;
	}
}

void complex2DFFTData::reScaleIn(void)
{
	double factor=1.0/((double)x*(double)y);
	for (int i = 0;i < x * y;i++)
	{
		c_re(in[i]) *=  factor; 
		c_im(in[i]) *=  factor;
	}
}

void complex2DFFTData::compPowerSpec()
{  
	maxPower = totalPower = 0;
	minPower = DBL_MAX;
	double temp;
	int maxIndexX, maxIndexY;
	int minIndexX, minIndexY;
	for (int i = 0; i < x; i++)// The x dimension
	{
		int index = i * y;
		for (int j = 0;j < y;j++)// The y dimension
		{ 
			temp = c_re(out[index+j]) * c_re(out[index+j]) + 
				c_im(out[index+j]) * c_im(out[index+j]);
			totalPower += (power_spectrum[index+j] = temp);
			_RPT3(_CRT_WARN, "%f\t%f\t%f\n",temp,maxPower,minPower);
			if (temp > maxPower)
				maxPower = temp; maxIndexX = i; maxIndexY = j;
			if (temp < minPower)
				minPower = temp; minIndexX = i; minIndexY = j;
		}
	}
}


void complex2DFFTData::compLogPowerSpec()
{  
	int maxIndexX, maxIndexY;
	int minIndexX, minIndexY;
	totalPower=0.0;
	maxPower = -DBL_MAX, minPower = DBL_MAX;
	double temp;
	std::ofstream ofs("2Dfft.dat");
	for (int i = 0; i < x; i++)  // The x dimension
	{
		int index = i * y;
		for (int j = 0;j < y;j++) // The y dimension
		{
			temp = c_re(out[index+j]) * c_re(out[index+j]) + 
				c_im(out[index+j]) * c_im(out[index+j]);
			temp=10.0*log10(temp);
			ofs<<temp<<'\t';
			totalPower += (power_spectrum[index+j]=temp);
			_RPT3(_CRT_WARN, "%f\t%f\t%f\n",temp, maxPower, minPower);
			if (temp > maxPower)
				maxPower = temp; maxIndexX = i; maxIndexY = j;
			if (temp < minPower)
				minPower = temp; minIndexX = i; minIndexY = j;
		}
		ofs <<std::endl;
	}
	ofs.close();
}

complex2DFFT::complex2DFFT(complex2DFFTData *d) 
{
	data = d;
	_RPT2(_CRT_WARN, "complexFFT init: X %d Y %d\n", data->getXSize(), data->getYSize());
	fwdPlan = 
		fftw_plan_dft_2d(data->getXSize(), data->getYSize(), data->in, data->out, FFTW_FORWARD, PLANTYPE);
	invPlan = 
		fftw_plan_dft_2d(data->getXSize(), data->getYSize(), data->out, data->in, FFTW_BACKWARD, PLANTYPE);
}

complex2DFFT::~complex2DFFT()
{
	_RPT0(_CRT_WARN, "complexFFT DeInit\n");
	fftw_destroy_plan(fwdPlan);
	fftw_destroy_plan(invPlan);
	_RPT0(_CRT_WARN, "complexFFT DeInit done\n");
}

void complex2DFFT::fwdTransform()
{
	if (!data)
		std::cerr<<"complex2DFFT::fwdTransform : data not present"<<std::endl;
	else
		fftw_execute(fwdPlan);
}

void complex2DFFT::invTransform()
{  
	if (!data)
		std::cerr<<"complex2DFFT::invTransform : data not present"<<std::endl;
	else
		fftw_execute(invPlan);
}
