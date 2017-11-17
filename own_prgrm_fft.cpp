#include <iostream>
#include <fftw3.h>
#include <cmath>

const int N = 160000;
const int REAL = 0;
const int IMAG = 1;

double f(double x);
void fill_data(fftw_complex* signal);
void normalize(fftw_complex* signal);
void print_data(fftw_complex* transf);

int main(void){

  //Declaration of storage arrays (class defined by fftw3)
  fftw_complex *in_data, *transf_data, *inv_trans;

  //Declaration of calculation plan (class defined dy fftw3)
  fftw_plan p, p1;
  
  //Allocation of storage arrays (assigns a size enough for data)
  in_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  transf_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  inv_trans = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

  //Plan must be initialized before input
  p = fftw_plan_dft_1d(N, in_data, transf_data, FFTW_BACKWARD, FFTW_ESTIMATE);
  p1 = fftw_plan_dft_1d(N, transf_data, inv_trans, FFTW_FORWARD, FFTW_ESTIMATE);
  //FFTW_ESTIMATE IS ENOUGH FOR A SINGLE PROCESS OF TRANSFORMING

  //Filling Process
  fill_data(in_data);
  
  //Computes DFT efficiently
  fftw_execute(p);
  normalize(transf_data);     //trasform is not normalized by default;

  //print_data(in_data);
  //print_data(transf_data);

  fftw_execute(p1);

  normalize(inv_trans);     //trasform is not normalized by default;
  print_data(inv_trans);
  
  //deallocates all
  fftw_destroy_plan(p);
  fftw_destroy_plan(p1);
  
}

double f(double x){
  return sin(400*x);
}

void fill_data(fftw_complex* signal){
  
  for(int ii = 0; ii < N; ii++){
    
    double y = ii * 0.31456;
    signal[ii][0] = f(y);
    signal[ii][1] = 0;
    
  }

}

void normalize(fftw_complex* transf){
  
  for(int ii = 0; ii < N; ii++){

    transf[ii][0] *= 1.0/sqrt(N);
    transf[ii][1] *= 1.0/sqrt(N);
    
  }
  
}

void print_data(fftw_complex* transf){

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
