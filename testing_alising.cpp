#include "fourier_routines.h"

const int N = 1600;
const int dt = 101;

int main(void){

  fftw_complex *in_data, *transf_data;
  
  in_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  transf_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

  FFTW_r(in_data, transf_data, N, dt);
  
  return 0;
  
}
