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
#ifndef REAL2DFFT_H_
#define REAL2DFFT_H_

#include <string.h>

#include <fftw3.h>
#ifndef fftw_real
#define fftw_real double
#endif
#define c_re(c) ((c)[0])
#define c_im(c) ((c)[1])

#include <iomanip>
using namespace std;

#ifndef PLANTYPE
#define PLANTYPE FFTW_ESTIMATE
#endif

/// class real2DFFTData controls and manipulates real 2D fft data
class real2DFFTData {
  /// x=row y=column
  int x, y;
  /// The total memory used by this class
  fftw_real *mem;
  /// Free the memory
  void memDeInit(void);
public:
  /// The input data and power spectrum
  fftw_real *in, *power;
  /// The output data
  fftw_complex *out;
  /// Arrays which sum across rows (x) and columns (y)
  fftw_real *xSum, *ySum;
  /// A sum across the input time signal
  fftw_real *timeXSum;
  /// Power spectral sums across rows (x) and columns (y)
  fftw_real *realXSum, *imagXSum;

  /// The total power in the power spectrum, the maximum and minimum powers too
  double totalPower, maxPower, minPower;
  /// The minimum/maximum row (x) and column (y) sums
  double xSumMin, xSumMax, ySumMin, ySumMax;
  /// Row (x) and Column (y) max sum indexes
  int maxXSumIndex, maxYSumIndex;

  /// Constructor with all memory to be allocated internally
  real2DFFTData(int r, int c, fftw_real * = 0);
  /// Deconstructor
  ~real2DFFTData();

  /// The row count
  int getXSize(){return x;}
  /// The column count
  int getYSize(){return y;}
  /// The half row count
  int getXHalfSize(){ if (!(x%2)) return x/2; else return x/2+1;}
  /// The half column count
  int getYHalfSize(){ if (!(y%2)) return y/2; else return y/2+1;}

  /// Scales the output down by the number of elements
  void reScale(void);
  /// This function computes the power spectrum and updates the totalPower, maxPower and minPower
  void compPowerSpec(); // Find the power spectrum
  /// Finds 10*log10(power spectrum) and updates the totalPower, maxPower and minPower
  void compLogPowerSpec(); // Find the log power spectrum

  /// Updates timeXSum
  void timeSpecAverage();
  /// Updates realXSum and imagXSum
  void complexSpecAverage();
  /// Finds the power Spectrum averages and 
  /// updates the xSumMin, xSumMax, ySumMin, ySumMax, xSum, ySum, maxXSumIndex, maxYSumIndex
  void powerSpecAverage();
  /// Finds the y-sum between columns start and stop
  void findYSum(int start, int stop);
  /// Finds the y-max for the ySum array, updates ySumMin, ySumMax, maxYSumIndex
  void findYMax(void);

  /// Zeros the in array
  void clearInput(void){memset(in, 0, x*2*(y/2+1)*sizeof(fftw_real));}
  /// Zeros the out awway
  void clearOutput(void){memset(out, 0, x*(y/2+1)*sizeof(fftw_complex));}
};

///class real2DFFT controls fftw plans and executes fwd/inv transforms
class real2DFFT {
  /// The forward and inverse plans
  fftw_plan fwdPlan, invPlan;
protected:
  /// The pointer to the relevant data
  real2DFFTData *data;
public:
  /// fft init ... data pointed to by 'd'
  real2DFFT(real2DFFTData *d);
  /// fft deconstructor
  ~real2DFFT();

  /// Forward transform the data (in to out)
  void fwdTransform(); // Forward 2D fft
  /// Inverse transform the data (out to in)
  void invTransform(); // Inverse 2D fft
};
/** \example real2DFFTExample.cc
 * This is an example of how to use the class.
 */
#endif // REAL2DFFT_H_
