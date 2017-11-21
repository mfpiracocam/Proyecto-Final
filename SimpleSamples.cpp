#include "fourier_proyect.h"
#include <fstream>

//This Program uses "Real" data obteained from sox (signal-generating app), and simulated
//sampling performed by simple c++ routines in order to illustrate aliasing and how to
//interpret DFT of a n-tuple.

void print(fftw_complex* transf, int N);

int main(void){
  
  //Read Data from input file. INSERT FILENAME BELLOW (!!!!!!!)
  std::ifstream data_file("14-output-0500-1000.dat");
  std::vector<double> time_data, func_data; //Vectors used to store data from data_file

  double t_real = 0.0;
  double f = 0.0;
  
  while(!data_file.eof()){
    data_file >> t_real;
    data_file >> f;
    time_data.push_back(t_real);
    func_data.push_back(f);
  }
  data_file.close();

  //In Hz. Represents the Sampling Rate of the sox-generated signal 
  double samp_frec = 24;
  
  const double w0 = 2*M_PI * samp_frec; //In rads per second 
  const double time = 26.0; //In sec
  const int N = floor(time*samp_frec); //Number of total data

  fftw_complex *in_data, *transf_data, *inv_data;
  
  in_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  transf_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  inv_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

  /*
  for(int ii = 0; ii < N; ii++){
    in_data[ii][0] = func_data[ii];
  }
  */
  
  //Simulated Data generated using c++ user-defined functions
  
  double x = 0;
  for(int jj = 0; jj < N; jj++){
    x = jj * time/N;
    in_data[jj][0] = PeriodicPulse_add(1,samp_frec/1000,N,jj,time); //Expected to interfere.
    in_data[jj][1] = 0.0;
  }
  
  //Compute DFT Using FFTW v.3.3.7 (Last version)
  FFTW_r(in_data, transf_data, N);
  IFFTW_r(transf_data, inv_data, N);
 
  //Print data. Print DFT spectra, in terms of wave number.
  normalize(transf_data,N);
  print_frecspace(transf_data,N);

  //print_RealData(time_data,func_data);
  
  return 0;
  
}

void print(fftw_complex* transf, int N){
  for(int ii = 0; ii < N; ii++){
    double norm = transf[ii][0]*transf[ii][0] + transf[ii][1]*transf[ii][1];
    std::cout << ii << "\t"
	      << transf[ii][0] << "\t"
	      << transf[ii][1] << "\t"
	      << norm << std::endl;
  }
}

