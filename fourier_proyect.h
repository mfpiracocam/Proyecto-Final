#include "fourier_routines.h"

void normalize(fftw_complex* signal, int N);
void print_frecspace(fftw_complex* transf, int N);
void print_timespace(fftw_complex* transf, int N, double dt);
void print_realFrecspace(fftw_complex* transf, int N, double sampfrec,double frecmax);

void FFTW_r(fftw_complex *in_data, fftw_complex *trans_data, int N, double dt);
void IFFTW_r(fftw_complex *in_data, fftw_complex *trans_data, int N, double dt);

double armonic_add(int nmin, int nmax, double frec, int NMAX, int jj, double time);
double definition(double x, double frec);
double PeriodicArbitraryPulse(double frec, int NMAX, int jj, double time);
double PeriodicPulse_add(int n,double frec, int NMAX, int jj, double time);
void print_RealData(const std::vector<double>& time, const std::vector<double>& func);

double definition(double x, double frec){
  double t0 = 1.0 / frec;
  return exp(-x*x);
}

double PeriodicArbitraryPulse(double frec, int NMAX, int jj, double time){
  double realt = jj*time / NMAX;
  double t0 = 1.0 / frec;
  if(realt <= t0){
    return definition(realt,frec);
  }
  if(realt > t0){
    return definition(realt - t0*floor(realt/t0),frec);
  }
}


double PeriodicPulse_add(int n, double frec, int NMAX, int jj, double time){
  double y = 0;
  for(int ii = 0; ii < n+1; ii++){
    y += PeriodicArbitraryPulse(ii*frec, NMAX, jj, time);
  }
  return y;
}
