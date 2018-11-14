#include  <iostream>
#include <cstdio>
#include <complex>
#include <fftw3.h>
#include <chrono>


int main(int argv, char* argc[]){
	//std::cout << argv <<std::endl;
	if(argv < 2 ){
		std::cout<< "Please enter an input" <<std::endl;
		return 1;
	} 

	int length = atoi(argc[1]);
	std::cout << std::endl <<"FFTW Length: " <<length  << std::endl;

	fftw_complex *in, *out;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * length);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * length);

	//put some random data in the input array
	//std::complex<double> I(0,1);
	//for(int i=0;i<length;i++) in[i] = (fftw_complex*) rand(); 

	fftw_plan p;
	std::chrono::time_point<std::chrono::system_clock> start, end; 

	start = std::chrono::system_clock::now(); 

	p = fftw_plan_dft_1d(length,
 						reinterpret_cast<fftw_complex*>(in), 
 						reinterpret_cast<fftw_complex*>(out), 
 						FFTW_BACKWARD, 
 						FFTW_PATIENT);

	end = std::chrono::system_clock::now(); 
	std::chrono::duration<double> plan_seconds = end - start; 
	std::cout << "Patient Planning Elapsed Time: " << plan_seconds.count() << " sec" <<std::endl;

	start = std::chrono::system_clock::now(); 
	fftw_execute(p);

	end = std::chrono::system_clock::now(); 
	std::chrono::duration<double> fft_seconds = end - start; 
	std::cout << "FFT Elapsed Time: " << fft_seconds.count() << " sec" <<std::endl;
	std::cout << "Goodbye!" <<std::endl;

	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);

	return 0;
}