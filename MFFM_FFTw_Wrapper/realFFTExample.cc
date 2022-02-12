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

#include <mffm/realFFT.H>

#define INPUTFILE "sine.1000Hz.dat"
#define OUTPUTFILE "powerSpectrum.txt"
#define OUTPUTFILE1 "in.txt"
#define OUTPUTFILE2 "out.txt"

int main (void){
  ifstream input(INPUTFILE);
  ofstream output(OUTPUTFILE);
  ofstream output1(OUTPUTFILE1);
  ofstream output2(OUTPUTFILE2);
  int count=0;
  double var;
  // Get the file size and check file exists ....
  if (!input){
    cout <<"input not opened !"<<endl;
    exit(-1);
  }
  while (input >> var)
    count++;
  //input.close();
  
  cout<<count<<" variables in file "<<INPUTFILE<<endl;

  //input.open(INPUTFILE);
  input.clear();
  input.seekg(0);
 
  realFFTData fftData(count);
  realFFT rfft(&fftData);

  // read data into data and rdata :
  for (int i=0; i<count; i++)
    input >> fftData.in[i];
  input.close();

  // forward transform :
  rfft.fwdTransform();

  // Find the power spectrum ...
  fftData.compPowerSpec();
  /*
  // inverse transform to check what happens (have to rescale too): 
  rfft.invTransform();
  */
  // output to file :
  for (int i=0; i<(count+1)/2; i++){
    output << fftData.power_spectrum[i]<<'\n';
  }
  //  cout <<(count+1)/2<<endl;
  output.close();

  for (int i=0; i<count; i++){
    output1 << fftData.in[i]<<'\n';
    output2 << fftData.out[i]<<'\n';
  }
  output1.close();
  output2.close();
  return 0;
}
