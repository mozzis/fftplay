#include "fftw++.h"

ifstream fftw::ifWisdom;
ofstream fftw::ofWisdom;
bool fftw::Wise=false;

// User settings:
unsigned int fftw::effort=FFTW_PATIENT;
const char *fftw::WisdomName="wisdom3.txt";
