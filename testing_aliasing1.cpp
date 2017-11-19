#include "fourier_routines.h"

double armonic_add(int n, double frec, int NMAX, int jj, double time);

int main(void){

  const double w_samp = 2.0*4.0*3.4567; //In Hz
  const double time = 5; //In sec
  const int N = floor(w_samp * time);

  const double w0 = 3.4567;
  
  fftw_complex *in_data, *transf_data, *inv_data;
  
  in_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  transf_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  
  //Fill n-tuple with data

  double y = 0;
  for(int jj = 0; jj < N; jj++){
    in_data[jj][0] = armonic_add(4,w0,N,jj,time) + cos(2.0*M_PI*2.0*(w0) * jj*time/N);
    in_data[jj][1] = 0.0;
  }
  
  //print_data(in_data,N);
  
  //Compute DFT Using FFTW v.3.3.7 (Last version)
  FFTW_r(in_data, transf_data, N);

  //Print data
  print_frecspace(transf_data,N);
  
  return 0;
  
}

double armonic_add(int n, double frec, int NMAX, int jj, double time){

  double y = 0;
  
  for(int ii = 0; ii < n+1; ii++){
    y += sin(2.0*M_PI*ii*frec * jj*time/NMAX) + cos(2.0*M_PI*ii*frec * jj*time/NMAX);
  }
  
}



