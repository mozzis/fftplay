include ../Makefile.local
EXAMPLES= #complexFFTExample real2DFFTExample #rfftwExample realFFTExample envelopeExample
SRC= realFFT.cc realFFTData.cc complexFFT.cc #FFTBase.cc  envelope.cc pSCorr.cc
SRC2D= real2DFFT.cc
OBJ=$(SRC:.cc=.o)
OBJ2D=$(SRC2D:.cc=.o)

all: $(OBJ) \
	$(OBJ2D) \
	libs \
	ex

libs:
	$(MAKELIB),libfft.so.1 -o libfft.so.1 $(OBJ)
	$(MAKELIB),lib2Dfft.so.1 -o lib2Dfft.so.1 $(OBJ2D)

ex:
	$(CC) *.o -lfftw3 complexFFTExample.cc -o complexFFTExample
	$(CC)  *.o -lfftw3 realFFTExample.cc -o realFFTExample
	$(CC) *.o -lfftw3 real2DFFTExample.cc -o real2DFFTExample
#	$(CC) -lrfftw -lfftw rfftwExample.cc -o rfftwExample

clean:
	rm -f $(EXAMPLES) *.o a.out *~ audio.fft.txt *.txt lib*.so.0.1 realFFTExample complexFFTExample real2DFFTExample


