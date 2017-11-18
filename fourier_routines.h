#include <iostream>
#include <fftw3.h>
#include <cmath>

void normalize(fftw_complex* signal, int N);
void print_data(fftw_complex* transf, int N);
void print_timespace(fftw_complex* transf, int N, double dt);

void FFTW_r(fftw_complex *in_data, fftw_complex *trans_data, int N, double dt);
void IFFTW_r(fftw_complex *in_data, fftw_complex *trans_data, int N, double dt);

void normalize(fftw_complex* transf, int N){
  
  for(int ii = 0; ii < N; ii++){

    transf[ii][0] *= 1.0/sqrt(N);
    transf[ii][1] *= 1.0/sqrt(N);
    
  }
  
}

void print_data(fftw_complex* transf, int N){

  double norm = 0;
  
  for(int ii = 0; ii < N; ii++){

    norm = sqrt(transf[ii][0]*transf[ii][0] + transf[ii][1]*transf[ii][1]);

    std::cout << ii << "\t"
	      << transf[ii][0] << "\t"
	      << transf[ii][1] << "\t"
	      << norm << " " << std::endl;

  }

  std::cout << std::endl;
    
}

void print_timespace(fftw_complex* transf, int N, double dt){

  double norm = 0;
  
  for(int ii = 0; ii < N; ii++){

    norm = sqrt(transf[ii][0]*transf[ii][0] + transf[ii][1]*transf[ii][1]);

    std::cout << ii * dt/N << "\t"
	      << transf[ii][0] << "\t"
	      << transf[ii][1] << "\t"
	      << norm << " " << std::endl;

  }

  std::cout << std::endl;
    
}

void FFTW_r(fftw_complex *in_data, fftw_complex *transf_data, int N){
  
  //Declaration of calculation plan (class defined dy fftw3)
  fftw_plan p;
  
  //Plan must be initialized before input
  p = fftw_plan_dft_1d(N, in_data, transf_data, FFTW_FORWARD, FFTW_ESTIMATE);
  
  //Computes DFT efficiently
  fftw_execute(p);
  normalize(transf_data, N);     //trasform is not normalized by default;
  
  //deallocates all
  fftw_destroy_plan(p);
  
}

void IFFTW_r(fftw_complex *in_data, fftw_complex *transf_data, int N){
  
  //Declaration of calculation plan (class defined dy fftw3)
  fftw_plan p;
  
  //Plan must be initialized before input
  p = fftw_plan_dft_1d(N, in_data, transf_data, FFTW_BACKWARD, FFTW_ESTIMATE);
  
  //Computes DFT efficiently
  fftw_execute(p);
  normalize(transf_data, N);     //trasform is not normalized by default;
  
  //deallocates all
  fftw_destroy_plan(p);
  
}
