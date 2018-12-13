
 #include <cstdio>

 //#include <vector>
 #include <iostream>
 #include <math.h>
 #include <string>
 #include <complex>
 #include <fftw3.h>
 #include <fstream>
 #include <chrono>
 #include <vector>
 #include "freqCavitySystem.h"
//number of points in FFTW vector. Needs to be power of two.
 #define NUMPOINTS 2048
 #define PI  3.141592653589793

 int main(int argc, char* argv[]){
 	double oneTTime, fourTTime, fourTfourHTime;
 	std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
	std::chrono::duration<double> elapsed;

 	/***************************************************
 	* One 50% triangular test case
 	*
 	****************************************************/
/*	
 	fftw_plan b;
 	fftw_complex *in, *out;
 	in = fftw_alloc_complex(NUMPOINTS);
 	out = fftw_alloc_complex(NUMPOINTS);

 	b = fftw_plan_dft_1d(NUMPOINTS, 
 						in,
 						out,
 						FFTW_BACKWARD, 
 						FFTW_PATIENT);


	start = std::chrono::high_resolution_clock::now();

 	freqCavitySystem test;
 	test.addCavity(1,0.50, 0.0,0.0);
  	test.fillF(in, NUMPOINTS);
 	fftw_execute(b);
 	fftShift(out,in, NUMPOINTS);

 	end = std::chrono::high_resolution_clock::now();

 	toFile(in, NUMPOINTS, "freq_1T.txt");


	fftw_free(in);
	fftw_free(out);
	fftw_destroy_plan(b);

	elapsed = end-start;
	oneTTime = 1000000*elapsed.count();
	std::cout << "One Triangular Takes " << oneTTime << " microseconds" <<std::endl;


*/

 	/***************************************************
 	* 4T with Values from Trey's thesis test case
 	*
 	****************************************************/
/* 
 	fftw_plan(b);
 	fftw_complex *bIn, *bOut;
 	bIn  = fftw_alloc_complex(NUMPOINTS);
 	bOut = fftw_alloc_complex(NUMPOINTS);
 	b = fftw_plan_dft_1d(NUMPOINTS,
 						bIn,
 						bOut, 
 						FFTW_BACKWARD, 
 						FFTW_PATIENT);


	start = std::chrono::high_resolution_clock::now();
	double Tr[8] ={ 1,
					1,
					1,
					1};

 	double R[8] = { 0.5865,
 					0.5960,
 					0.5721,
 					0.5786};

 	double delta[8] = { -2.242,
 						-3.0307,
 						 2.4895,
 						-0.6930};
 	freqCavitySystem test;

 	for(int j=0;j<4;j++){
 		test.addCavity(Tr[j], R[j], delta[j], 0.0);
 	}
 	test.fillF(bIn, NUMPOINTS);
 	fftw_execute(b);
 	fftShift(bOut, bIn, NUMPOINTS);

 	end = std::chrono::high_resolution_clock::now();
 	toFile(bIn, NUMPOINTS, "freq_4T.txt");



	fftw_free(bIn);
	fftw_free(bOut);
	fftw_destroy_plan(b);

	elapsed = end-start;
	fourTTime = 1000000*elapsed.count();
	std::cout << "Four Triangular Takes " << fourTTime << " microseconds" <<std::endl;

*/

 	/***************************************************
 	* 4+4 with Values from Trey's thesis test case
 	*
 	****************************************************/


 	fftw_plan(b);
 	fftw_complex *bIn, *bOut;
 	bIn  = fftw_alloc_complex(NUMPOINTS);
 	bOut = fftw_alloc_complex(NUMPOINTS);
 	b = fftw_plan_dft_1d(NUMPOINTS,
 						bIn,
 						bOut, 
 						FFTW_BACKWARD, 
 						FFTW_PATIENT);

	start = std::chrono::high_resolution_clock::now();

	double Tr[8] ={ 1,
					1,
					1,
					1,
					9,
					9,
					9,
					9};

 	double R[8] = { 0.5865,
 					0.5960,
 					0.5721,
 					0.5786,
 					0.5865,
 					0.5960,
 					0.5721,
 					0.5786};

 	double delta[8] = { -2.242,
 						-3.0307,
 						 2.4895,
 						-0.6930,
 						-2.242,
 						-3.0307,
 						 2.4895,
 						-0.6930};

	freqCavitySystem test;


 	for(int j=0;j<8;j++){
 		test.addCavity(Tr[j], R[j], delta[j], 0.0);
 	}


 	test.fillF(bIn, NUMPOINTS);
 	fftw_execute(b);
 	fftShift(bOut, bIn, NUMPOINTS);

 	end = std::chrono::high_resolution_clock::now();
 	toFile(bIn, NUMPOINTS, "freq_4T4H.txt");

	elapsed = end-start;
	fourTfourHTime = 1000000*elapsed.count();
	std::cout << "4+4 Takes " << fourTfourHTime << " microseconds" <<std::endl;


 	/***************************************************
 	* 4+4 Time Testing - 100 iterations of each.
 	*
 	****************************************************/

/*
 	
 	double setupTime, fillTime, fftTime, shiftTime;
 	std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
	std::chrono::duration<double> elapsed;

 	double Tr[8] = {1.0,1.0,1.0,1.0,9.0,9.0,9.0,9.0};
 	double R[8] = {0.6855,0.6950,0.6711,0.6776,0.6855,0.6950,0.6711,0.6776};
 	//double R[8] = {0.599,0.599,0.599,0.599,0.599,0.599,0.599,0.599};
 	double delta[8] = {-2.242,-3.0307,2.4895,-0.6930,-2.242,-3.0307,2.4895,-0.6930};
 	//double delta[8] = {};

 	std::vector<freqCavitySystem> cavVec(100);
 	start = std::chrono::high_resolution_clock::now();

 	for(int i=0;i<100;i++){
 		for(int j=0;j<8;j++){
 			cavVec[i].addCavity(Tr[j], R[j] - i/1000.0, delta[j], 0.0);
 		}

 	}
 	end = std::chrono::high_resolution_clock::now();
 	elapsed = end-start;
 	setupTime = elapsed.count();
 	setupTime = 1000000*setupTime/100.0;

 	std::cout << cavVec[99].toString()<<std::endl;


 	//std::cout<<cavVec[0].toString()<<std::endl;
 	//std::cout<<cavVec[76].toString()<<std::endl;

 	fftw_complex *fIn, *fOut, *bIn, *bOut;
 	fIn  = fftw_alloc_complex(NUMPOINTS);
 	fOut = fftw_alloc_complex(NUMPOINTS);
 	bIn  = fftw_alloc_complex(NUMPOINTS);
 	bOut = fftw_alloc_complex(NUMPOINTS);

 	//allocate the fftw plan and calculate it 
 	fftw_plan f, b;
 	b = fftw_plan_dft_1d(NUMPOINTS,
 						bIn,
 						bOut, 
 						FFTW_BACKWARD, 
 						FFTW_PATIENT);
 	f = fftw_plan_dft_1d(NUMPOINTS,
 						 fIn, 
 						 fOut,
 						 FFTW_FORWARD,
 						 FFTW_ESTIMATE);



 	start = std::chrono::high_resolution_clock::now();
 	for(int i=0;i<100;i++){
 		cavVec[i].fillF(bIn, NUMPOINTS);
 	}
 	end = std::chrono::high_resolution_clock::now();
 	elapsed = end-start;
 	fillTime =elapsed.count();

 	start = std::chrono::high_resolution_clock::now();
 	for(int i=0;i<100;i++){
 		fftw_execute(b);
 	}
 	end = std::chrono::high_resolution_clock::now();
 	elapsed = end-start;
 	fftTime = elapsed.count();

 	start = std::chrono::high_resolution_clock::now();
 	for(int i=0;i<100;i++){
 		fftShift(bOut, fIn, NUMPOINTS);
 	}
 	end = std::chrono::high_resolution_clock::now();
 	elapsed = end-start;
	shiftTime=elapsed.count();
 	



 	shiftTime =1000000*shiftTime/100.0;
 	fftTime =  1000000*fftTime/100.0;
 	fillTime = 1000000*fillTime/100.0;

 	std::cout << "NUMPOINTS = " << NUMPOINTS <<std::endl;
 	std::cout << "SINELENGTH = " << SINELENGTH << std::endl;
 	std::cout<<"Setup Time = " << setupTime <<std::endl;
 	std::cout<< "Shift Time = " << shiftTime <<std::endl;
 	std::cout<< "fft Time = " << fftTime <<std::endl;
 	std::cout<< "fill Time = " << fillTime <<std::endl;
 	toFile(fIn,NUMPOINTS, "test.txt");
/*


 	std::string str1 = "fluxTest" + std::to_string(NUMPOINTS) + ".txt";
 	toFile(fIn, NUMPOINTS, str1);

 	std::cout<< NUMPOINTS << "," <<
 				SINELENGTH << "," << 
 				setupTime << "," <<
 				shiftTime << "," <<
 				fftTime << ","   <<
 				fillTime << ","  <<
 				std::endl;

 	fftw_free(fIn);
 	fftw_free(fOut);
 	fftw_free(bIn);
 	fftw_free(bOut);

/*
 	double gh, sqrtTime;
 	start = std::chrono::high_resolution_clock::now();
 	for(int i = 0; i<NUMPOINTS;i++){
 		gh = std::sqrt(i*i);
 	}
 	end = std::chrono::high_resolution_clock::now();
 	elapsed = end-start;
	sqrtTime=elapsed.count();
	sqrtTime = 1000000*sqrtTime/100.0;
	std::cout<< "sqrt takes = " << sqrtTime	<< std::endl;
 	*/


 	/*
 	double testX[NUMPOINTS] = {};
 	linspace(0, 55, NUMPOINTS,testX);
 	double the[NUMPOINTS] = {};
 	double csine[NUMPOINTS] = {};
 	for(int i=0;i<NUMPOINTS;i++){
 		the[i] = test.bSine(testX[i]);
 		csine[i] = sin(testX[i]);
 	}
 	toFile(the, NUMPOINTS, "LUTSineVals.txt");
 	toFile(testX, NUMPOINTS, "xvals.txt");
 	toFile(csine, NUMPOINTS, "sine.txt");
*/
 	return 0;
 }