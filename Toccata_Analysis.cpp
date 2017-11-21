#include "fourier_proyect.h"
#include <fstream>

//This program is intended to perform fourier analysis of Toccata and Fuge in
//d-minor for different sampling frequencies. This may show the importance of
//the right choosing of the sampling rate. IMPORTANT: Data is scaled in a fac-
//tor of one thousand of the values obtained.

//.txt files from: http://cmp.fmf.uni-lj.si/. Computational Methods for Physi-
//cists. Simon Sirca and Martin Horvat.

int main(void){

  int N = 0.0; 
  double frec = 11025;
  std::vector<double> data;
  
  fftw_complex *in_data, *transf_data, *inv_data;
  
  //Read Data from input file
  std::ifstream data_file("14-chorzempa.11025.txt");

  double d = 0.0;
  while(!data_file.eof()){
    data_file >> d;
    data.push_back(d);
  }
  data_file.close();
  
  //Initialize n-tuple for DFT
  N = data.size();
  
  in_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  transf_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  inv_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

  //Fill n-tuple with data
  for(int ii = 0; ii < N; ii++){
    in_data[ii][0] = data[ii]/1000.0;
    in_data[ii][1] = 0.0;
  }

  //Compute DFT Using FFTW v.3.3.7 (Last version)
  FFTW_r(in_data, transf_data, N);

  //Print data
  normalize(transf_data,N);
  print_timecspace(transf_data,N,frec,22050);
  
  return 0;
  
}




