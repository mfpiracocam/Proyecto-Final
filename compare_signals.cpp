#include "fourier_routines.h"

double armonic_add(int n, double frec, int NMAX, int jj, double time);
double definition(double x, double frec);
double PeriodicArbitraryPulse(double frec, int NMAX, int jj, double time);

int main(void){

  const double w0 = 30; //In Hz
  const double T0 = 2*M_PI / w0;
  
  const double w_samp = 2.0 * 100 * w0; //In Hz. The 2.0 factor relates to w_nyquist
  const double time = 2.0 * T0; //In sec
  const int N = floor(w_samp * time);

  fftw_complex *in_data, *transf_data, *inv_data;
  
  in_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  transf_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  inv_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  
  //Fill n-tuple with data

  double y = 0;
  for(int jj = 0; jj < N; jj++){
    in_data[jj][0] = PeriodicArbitraryPulse(w0,N,jj,time);
    in_data[jj][1] = 0.0;
  }
  
  //Compute DFT Using FFTW v.3.3.7 (Last version)
  FFTW_r(in_data, transf_data, N);
  IFFTW_r(transf_data, inv_data, N);
  
  //Print data
  print_timespace(inv_data,N,time);
  
  return 0;
  
}

double armonic_add(int n, double frec, int NMAX, int jj, double time){

  double y = 0;
  
  for(int ii = 0; ii < n+1; ii++){
    y += sin(2.0*M_PI*ii*frec * jj*time/NMAX) + cos(2.0*M_PI*ii*frec * jj*time/NMAX);
  }

  return y;
  
}

double definition(double x, double frec){
  
  double t0 = 2*M_PI / frec;

  return x;

  /*
  if(x < t0/2.0)
    return x;
  if(x >= t0/2.0 && x <= t0)
    return 1-x;
  */
  
}

double PeriodicArbitraryPulse(double frec, int NMAX, int jj, double time){
  double realt = jj*time / NMAX;
  double t0 = 2*M_PI / frec;
  if(realt <= t0)
    return definition(realt,frec);
  if(realt > t0){
    return definition(realt - t0*floor(realt/t0),frec);
  }

}



