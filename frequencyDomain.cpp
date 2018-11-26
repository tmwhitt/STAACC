//! Freq Domain test code
/*! \brief This code runs the frequency domain propagation of a stacker
 * 
 * \author Alex "A-drizzle"  Rainville
 * \date 11/16/2018
 * \version 1.2 fixed fftw_complex indexing and transfer function algorithm
 * \note All units mks unless noted in program
 * \par
 * loading module on flux: -bash-4.2$ module load fftw/3.3.4/gcc/4.8.5
 * \par
 * to compile on flux: -bash-4.2$ g++ -I$FFTW_INCLUDE -L$FFTW_LIB -lfftw3 -lm <filename.cpp> -o <output file name>
 * 
 * g++ frequencyDomain.cpp -o testDelete.exe -I/home/rainvila/fftw-3.3.8-lib/include -L/home/rainvila/fftw-3.3.8-lib/lib -lfftw3 -lm -std=c++11
*
 * 
 * \bug must call update on each stacker before adding another or this will not work 
 */

 #include <cstdio>

 //#include <vector>
 #include <iostream>
 #include <cmath>
 #include <string>
 #include <complex>
 #include <fftw3.h>
 #include <fstream>
 #include <chrono>
//number of points in FFTW vector. Needs to be power of two.
 #define NUMPOINTS 4096




/**
 * @brief Creates a linearly increasing array
 * @details modfies input array to have a linear increase from a to b, with a specified
 * number of steps. 
 * 
 * @param a start value
 * @param b end value
 * @param points length of array
 * @param x pointer to head of the array.
 */
void linspace(int a, int b, int points, double* x){
	//std::cout<< "called linspace" <<std::endl;
	double step = (b-a)/(points-1.0);
	double val = a;
	for(int i=0;i<=points;i++){
		//store value at that point
		*x = val;
		//increment pointer
		x++;
		//increment value
		val+=step;
	}
}


/**
 * @brief Writes a fftw_complex array to file
 * @details Writes a fftw_complex array to file name passed in. The file is is CSV format with real, imaginary parts of the number. There are no header lines so far.
 * 
 * @param a pointer to head of fftw_complex array
 * @param N length of array
 * @param filename name of file
 */
void toFile(fftw_complex* a, int N, std::string filename ){
	std::ofstream myfile (filename);
	if(myfile.is_open()){
		for(int i=0;i<N;i++){
			myfile << a[i][0] << "," << a[i][1] << "\n";;
		}

	}else std::cout << "cannot open file"<< std::endl;


	myfile.close();
}


/**
 * @brief Does a convolution with for loops
 * @details Convolves two same-length arrays and puts them in output array c. c must be of length 2N.
 * 
 * @param a Pointer to head to first array
 * @param b Pointer to head of second array
 * @param c Pointer to head of convoultion output array
 * @param N Length of arrays.
 */
void convolve(fftw_complex* a, fftw_complex *b, fftw_complex *c, int N){
	double d,e,f,g,k,l,m,n;
	for(int j=1;j<N;j++){
		double imagSumFor =0.0;
		double realSumFor =0.0;
		double imagSumBack =0.0;
		double realSumBack =0.0;
		int indy;
		for(int i=0;i<j;i++){
			indy = N-j+i;
			d = a[i][0];
			e = a[i][1];
			f = b[indy][0];
			g = b[indy][1];

			k = a[indy][0];
			l = a[indy][1];
			m = b[i][0];
			n = b[i][1];
			realSumFor += d*f-e*g;
			imagSumFor += d*g+e*f;

			realSumBack += k*m - l*n;
			imagSumBack += k*n + l*m;
		}
		c[j][0] = realSumFor;
		c[j][1] = imagSumFor;
		c[N-j][0] = realSumBack;
		c[N-j][1] = imagSumBack;
	}
	
	double realSum = 0.0;
	double imagSum = 0.0;
	for(int i=0;i<N;i++){
		d = a[i][0];
		e = a[i][1];
		f = b[i][0];
		g = b[i][1];
		realSum +=d*f - e*g;
		imagSum += d*g + e*f;
	}
	c[N/2-1][0] = realSum;
	c[N/2-1][1] = imagSum;
	
}

/**
 * @brief Re-aligns output fft array
 * @details Takes output fftw_complex array and aligns it from negative to positive frequency. Written using pointers for speed
 * 
 * @param a Input array to be aligned
 * @param b Output aligned array
 * @param Length of array
 */
void reAlign(fftw_complex *a, fftw_complex *b, int N){
	double *realIn, *imagIn, *imagOut, *realOut;
	
	realIn = &a[N/2-1][0];
	imagIn = &a[N/2-1][1];
	realOut = &b[0][0];
	imagOut = &b[0][1];

	std::cout<< realIn << std::endl;
	for(int i=0;i<N/2;i++){
		*realOut = *realIn;
		*imagOut = *imagIn;
		realOut+=2;
		imagOut+=2;
		realIn+=2;
		imagIn+=2;
	}

	std::cout <<realIn << std::endl;
	realIn = &a[0][0];
	imagIn = &a[0][1];
	realOut = &b[N/2-1][0];
	imagOut = &b[N/2-1][1];
	for(int i=N/2;i<N;i++){
		//std::cout <<i << " : "<< *realIn << "," << *imagIn <<std::endl;
		*realOut = *realIn;
		*imagOut = *imagIn;
		realOut+=2;
		imagOut+=2;
		realIn+=2;
		imagIn+=2;
	}
}

class stacker{
	private:
		double Tr; //cavity roundtrip time 
		double delta; //cavity phase
		double R; //cavity reflectivity (power)
		double alpha; //cavity absorption, in [1/cm]
		int stackerID;
		static double x[NUMPOINTS]; //linear vector from 0 to 1 
		static int number;
		static std::complex<double> F[NUMPOINTS];
		fftw_complex *toReturn;
		fftw_complex *conjToReturn;
		static bool gottenF, gottenConjF;

	public:
		//constructors
	
		//constructor with cavity absorption
		stacker(double Tr, double delta, double R, double alpha);
		//destructor 
		~stacker();

		//getters 
		int getStackerID(){return stackerID;}
		double getTr(){return Tr;}
		double getR(){return R;}
		double getAlpha(){return alpha;}
		double getDelta(){return delta;}
		//returns a pointer to the head of the overallF array
		fftw_complex* getF();
		fftw_complex* getConjF();

		//setters (do not allow for setting of stacker ID, stacker numbers)
		void setTr(double Tr){this->Tr=Tr;}
		void setR(double R){this->R=R;}
		void setAlpha(double alpha){this->alpha = alpha;}
		void setDelta(double delta){this->delta=delta;}


		//method to update the overall F(w) array.
		void updateF();

		//tostring method
		std::string toString();


};

int stacker::number = 0;
double stacker::x[] ={};
std::complex<double> stacker::F[NUMPOINTS] = {};
bool stacker::gottenF = false;
bool stacker::gottenConjF = false;

stacker::stacker(double Tr, double delta, double R, double alpha){
	this -> Tr = Tr;
	this -> delta = delta;
	this -> R = R;
	this -> alpha = alpha;
	number++;
	stackerID = number;
	//std::cout << number << " " <<stackerID << std::endl;
	if(stackerID==1) linspace(0,1,NUMPOINTS,x);
	updateF();
}



stacker::~stacker(){
	//when the last stacker is destroyed remove the fftw_complex arrays
	if(stackerID==1){
		if(gottenF==true) fftw_free(toReturn);
		else if(gottenConjF == true) fftw_free(conjToReturn);
	}
}

fftw_complex* stacker::getF(){
	if(gottenF==false){
		toReturn  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NUMPOINTS);
		for(int i=0;i<NUMPOINTS;i++){
			toReturn[i][0] = real(F[i]);
			toReturn[i][1] = imag(F[i]);
		}
	}
	return toReturn;
}

fftw_complex* stacker::getConjF(){
	if(gottenConjF==false){
		conjToReturn  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NUMPOINTS);
		for(int i=0;i<NUMPOINTS;i++){
			conjToReturn[i][0] = -1.0*imag(F[i]);
			conjToReturn[i][1] = real(F[i]);
		}
		gottenConjF = true;
	}
	return conjToReturn;
}



/**
 * @brief Calculates overall frequency transfer function
 * @details Will update the static class array F by either (1) initializing it if this is the first stacker or (2) multiplying this stacker's transfer function by the previous system transfer function.
 */
void stacker::updateF(){
	std::complex<double> r(std::sqrt(R),0.0);
	std::complex<double> I(0.0,1.0);
	std::complex<double> deltaL;
	std::complex<double> expVal;

	double *xBegin = &x[0];
	std::complex<double> *fBegin = &F[0];

	if(stackerID==1){
		//std::cout << "called first stacker, assigning values to overallF" <<std::endl;


		for(int i=0;i<NUMPOINTS;i++){
			deltaL = 2*M_PI*(Tr*(*xBegin) +delta);
			expVal = std::exp(I*deltaL);
			*fBegin = (r-expVal)/(1.0-r*expVal);
			xBegin++;
			fBegin++;
		}
	} 
	else{
		//std::cout << "called a n>1 stacker, updating overall F" <<std::endl;
		for(int i=0;i<NUMPOINTS;i++){
			deltaL = 2*M_PI*(Tr*(*xBegin) +delta);
			expVal = std::exp(I*deltaL);
			*fBegin = (r-expVal)/(1.0-r*expVal)*(*fBegin);
			fBegin++;
			xBegin++;
		}
	}
}


/**
 * @brief Returns string about stacker
 * @details Returns a string with the stacker id, number of stackers, R, Tr, delta, and alpha
 * @return See details.
 */
std::string stacker::toString(){
	std::string str1 = "Stacker " + std::to_string(stackerID) + " of " + std::to_string(number) +": Tr=" + std::to_string(Tr) + " R = "  + std::to_string(R) + " delta = " + std::to_string(delta) + " alpha = " + std::to_string(alpha);
	return str1;
}


//wanna tell me if my code is insane/bad form? it works.




 int main(int argc, char* argv[]){


 	std::chrono::time_point<std::chrono::system_clock> start, end, oStart, oEnd; 

 	oStart = std::chrono::system_clock::now();

 	//stcker takes in Tr, delta, R, alpha
 	double Tr1 = 1.0;
 	double Tr2 = 9.0;
 	double delta = 0.0;
 	double R1 = 0.58;
 	double R2 = 0.68;
 	double alpha = 0.0;
 	//construct 8 stackers (in 4+4) config

 	/*
 	stacker first (Tr1, -2.242,R1,alpha);
  	stacker second (Tr1, -3.037,R1,alpha);
 	stacker third (Tr1, 2.4895,R2,alpha);
 	stacker fourth (Tr1, -0.693,R2,alpha);
 	stacker fifth (Tr2, -2.242,R1,alpha);
 	stacker sixth (Tr2, -3.037,R1,alpha);
 	stacker seventh (Tr2, 2.4895,R2,alpha);
 	stacker eighth (Tr2, -0.693,R2,alpha);
	*/

 	start = std::chrono::system_clock::now();
 	stacker first   (Tr1, 0.0, 0.5, alpha);
 	end = std::chrono::system_clock::now();
 	std::chrono::duration<double> firstCav = end-start;

 	start  =std::chrono::system_clock::now();
  	stacker second  (Tr1, 0.0, 0.5, alpha);
  	end = std::chrono::system_clock::now();
  	std::chrono::duration<double> secondCav = end-start;

 	start  =std::chrono::system_clock::now();
 	stacker third   (Tr1, 0.0, 0.5, alpha);
 	end = std::chrono::system_clock::now();
  	std::chrono::duration<double> thirdCav = end-start;

  	start  =std::chrono::system_clock::now();
 	stacker fourth  (Tr1, 0.0, 0.5, alpha);
	end = std::chrono::system_clock::now();
  	std::chrono::duration<double> fourthCav = end-start;

 	stacker fifth   (Tr2, 0.0, 0.5, alpha);
 	stacker sixth   (Tr2, 0.0, 0.5, alpha);
 	stacker seventh (Tr2, 0.0, 0.5, alpha);
 	stacker eighth  (Tr2, 0.0, 0.5, alpha);


 	/***********************************************************************************
 	*Calculate the impulse response and time it all
 	*
 	*Put output in IR, and the properly aligned output in alignedIR
 	*
 	*********************************************************************************/
 	//create two complex fftw arrays for input and output
 	fftw_complex *IR, *alignedIR, *Fin;
 	IR = fftw_alloc_complex(NUMPOINTS);
 	alignedIR = fftw_alloc_complex(NUMPOINTS);
 	Fin = eighth.getConjF();

 	//allocate the fftw plan and calculate it 
 	fftw_plan p;
 	start = std::chrono::system_clock::now();
 	p = fftw_plan_dft_1d(NUMPOINTS,
 						Fin,
 						IR, 
 						FFTW_BACKWARD, 
 						FFTW_ESTIMATE);
 	end = std::chrono::system_clock::now();
 	std::chrono::duration<double> planTime = end-start;

 	//do the fft
 	start = std::chrono::system_clock::now();
 	fftw_execute(p);
 	end = std::chrono::system_clock::now();

 	std::chrono::duration<double> fftwTime = end-start;

 	/*

 	//optionally output information about the plan
 	fftw_print_plan(p) ;
 	std::cout<<std::endl;

	*/


 	//align the output from negative to positive frequency
 	start = std::chrono::system_clock::now();
 	reAlign(IR, alignedIR, NUMPOINTS);
 	end =  std::chrono::system_clock::now();

 	std::chrono::duration<double> alignTime = end-start;
 	oEnd = std::chrono::system_clock::now();
 	std::chrono::duration<double> overallTime = oEnd-oStart;

 	toFile(alignedIR,NUMPOINTS, "test.txt");

 	//std::cout << "Overall Time , First Cavity Setup, Second Cavity Setup, Third Cavity Setup, Fourth Cavity Setup, Plan Time, fftw Time, Align Time" << std::endl;
 	std::cout << NUMPOINTS << ","
 				<< overallTime.count() << "," 
 				<< firstCav.count() << ","  
 				<< secondCav.count() << "," 
 				<< thirdCav.count() << "," 
 				<< fourthCav.count() << "," 
 				<< planTime.count() << "," 
 				<< fftwTime.count() << "," 
 				<< alignTime.count()<< std::endl;

 	for(int i=0;i<10;i++){
 		std::cout << &IR[i] << "," << &IR[i][1] <<std::endl;
 	}


 	/*********************************************************************************
 	* Create the burst profile 
 	*
 	*Burst called burst, burst FT called burstFT
 	*
 	******************************************************************************/
/*
 	//create the stacking burst functions with 81 pulses at the center of a NUMPOINTS length array
 	fftw_complex *burst, *burstFT;
 	fftw_plan pp;
 	burst = fftw_alloc_complex(NUMPOINTS);
 	burstFT = fftw_alloc_complex(NUMPOINTS);

 	//set everything to 0
 	for(int i=0;i<NUMPOINTS;i++){
 		burst[i][0] = 0;
 		burst[i][1] = 0;
 	}

 	burst[NUMPOINTS/2-1][0] = 1.0;
/*
 	//set the middle 81 pulses to equal amplitude no phase
 	for(int i=0; i<81;i++){
 		burst[i+NUMPOINTS/2 -81][0] = 1.0;
 	}
*/
 	/*
 	pp = fftw_plan_dft_1d(NUMPOINTS,
 						burst,
 						burstFT,
 						FFTW_FORWARD,
 						FFTW_ESTIMATE);
 	fftw_execute(pp);
 	fftw_destroy_plan(pp);
 	//toFile(out,NUMPOINTS, "test.txt");


 	/******************************************************************************
 	*Multiply the burst FT by the frequency domain impulse response of the function
 	*
 	*multipled array called "multiplied"
 	*output array called "output"
 	*
 	****************************************************************************/
/*
 	fftw_plan ppp;
 	fftw_complex *multiplied, *output; 
 	multiplied = fftw_alloc_complex(NUMPOINTS);
 	output = fftw_alloc_complex(NUMPOINTS);
 	//fftw_complex *F = first.getConjF();
/*
 	for(int i=0;i<NUMPOINTS;i++){
 		multiplied[i][0] = burst[i][0] * F[i][0];
 		multiplied[i][1] = burst[i][1] * F[i][1];
 	}
*/
 	/*
 	ppp = fftw_plan_dft_1d(NUMPOINTS,
 						multiplied,
 						output,
 						FFTW_BACKWARD,
 						FFTW_ESTIMATE);

 	fftw_execute(ppp);

 	fftw_destroy_plan(ppp);
 	//toFile(output, NUMPOINTS,"test.txt");


 	//create the output convolution array
 	fftw_complex *conV;
 	conV = fftw_alloc_complex(2*NUMPOINTS);

 	//convolve 

	start = std::chrono::system_clock::now(); 
 	convolve(alignedIR, burst, conV, NUMPOINTS);
 	end = std::chrono::system_clock::now();

 	std::chrono::duration<double> convTime = end - start; 
	

 	//print to file

 	//toFile(conV, NUMPOINTS*2, "test.txt");

 	
 	std::chrono::duration<double> overallTime = oEnd-oStart;

 	std::cout<< "NUMPOINTS, Overall Time, Plan Time, FFTW time, Align Time, Convolve Time" <<std::endl;

 	std::cout << NUMPOINTS << ","<< overallTime.count() << "," <<planTime.count() << "," << fftwTime.count() <<"," << alignTime.count() << "," <<convTime.count() << std::endl;
 	*/
 	
 	fftw_free(IR);
 	fftw_free(alignedIR);
 	//fftw_free(conV);
 	//fftw_free(burst);
 	//fftw_free(multiplied);
 	//fftw_free(output);
 	
 	return 0;
 }
