#include "fourier_routines.h"
#include <fstream>
#include <vector>

int main(void){

  int N = 10;
  std::vector<double> data;
  
  fftw_complex *in_data, *transf_data, *inv_data;

  //Read Data from input file
  std::ifstream data_file("14-chorzempa.882.txt");

  int ii = 0;
  double d = 0.0;
  while(!data_file.eof()){
    data_file >> d;
    data.push_back(d);
    ii++;
  }
  data_file.close();
  
  //Initialize n-tuple for DFT
  N = data.size();

  in_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  transf_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  inv_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

  //Fill n-tuple with data
  for(int ii = 0; ii < N; ii++){
    in_data[ii][0] = data[ii];
    in_data[ii][1] = 0.0;
  }
  
  //Compute DFT Using FFTW v.3.3.7 (Last version)
  FFTW_r(in_data, transf_data, N);

  //Print data
  print_data(transf_data,N);
  
  return 0;
  
}




