//! Freq Domain test code
/*! \brief This code runs the frequency domain propagation of a stacker
 * 
 * \author Alex "A-drizzle"  Rainville
 * \date 11/27/2018
 * \version 2.0 Changed to freqCavity that represents a stacker system
 * \note All units mks unless noted in program
 * \par
 * loading module on flux: -bash-4.2$ module load fftw/3.3.4/gcc/4.8.5
 * \par
 * to compile on flux: -bash-4.2$ g++ -I$FFTW_INCLUDE -L$FFTW_LIB -lfftw3 -lm <filename.cpp> -o <output file name>
 * 
 * g++ frequencyDomain.cpp -o testDelete.exe -I/home/rainvila/fftw-3.3.8-lib/include -L/home/rainvila/fftw-3.3.8-lib/lib -lfftw3 -lm -std=c++11
*
 * 
 * \bug 
 * \TODO Implement "remove cavity" in freqCavity. 
 * \TODO Implement faster calculation of F. 
 * \TODO figure out dealing with x
 */

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
//number of points in FFTW vector. Needs to be power of two.
 #define NUMPOINTS 4096
 #define PI  3.141592653589793


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
	for(double val=a;val<=b;val+=step){
		//store value at that point
		*x = val;
		//std::cout<< val << ","<< *x << std::endl;
		//increment pointer
		x++;
	}
}

void vecMult(fftw_complex *a, fftw_complex *b, fftw_complex *c, int N){
	double *aReal, *aImag, *bReal, *bImag, *cReal, *cImag;

	aReal = &a[0][0];
	aImag = &a[0][1];
	bReal = &b[0][0];
	bImag = &b[0][1];
	cReal = &c[0][0];
	cImag = &c[0][1];

	for(int i=0;i<N;i++){
		*cReal = (*aReal)*(*bReal) - (*aImag)*(*bImag);
		*cImag = (*aImag)*(*bReal) + (*aReal)*(*bImag);	

		aReal +=2;
		aImag +=2;
		bReal +=2;
		bImag +=2;
		cReal +=2;
		cImag +=2;	
	}
}

void conj(fftw_complex *a, int N){
	double *aImag;
	for(int i=0;i<N;i++){
		*aImag = -1.0*(*aImag);
		aImag+=2;
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
 * @brief Re-aligns output fft array
 * @details Takes output fftw_complex array and aligns it from negative to positive frequency.
 * 
 *  \note  Written using pointers for speed. FFTW is 16-byte aligned and requires +=2 to go to the next element in the array. 
 * 
 * @param a Input array to be aligned
 * @param b Output aligned array
 * @param Length of array
 */
void fftShift(fftw_complex *a, fftw_complex *b, int N){
	double *realIn, *imagIn, *imagOut, *realOut;
	///swap the second half of a into the first half of b
	realIn = &a[N/2][0];
	imagIn = &a[N/2][1];
	realOut = &b[0][0];
	imagOut = &b[0][1];

	for(int i=0;i<N/2;i++){
		*realOut = *realIn;
		*imagOut = *imagIn;
		realOut+=2;
		imagOut+=2;
		realIn+=2;
		imagIn+=2;
	}

	///Swap the first half of b into the second half of b
	realIn = &a[0][0];
	imagIn = &a[0][1];
	realOut = &b[N/2][0];
	imagOut = &b[N/2][1];
	for(int i=N/2;i<N;i++){
		*realOut = *realIn;
		*imagOut = *imagIn;
		realOut+=2;
		imagOut+=2;
		realIn+=2;
		imagIn+=2;
	}
}

/**
 * @brief Frequency Domain cavity
 * @details A cavity object that describes a system of stackers. 
 * 
 */
class freqCavity{
private:
	fftw_complex *F; /*! Impulse response of the system */
	int totalNumber =0; /*! Total number of individual stackers inside */
	std::vector<int> cavity; /*! vector of ints representing each stacker*/
	std::vector<double> R; /*! vector of R for each stacker*/
	std::vector<double> Tr; /*!vector of Tr for each stacker*/
	std::vector<double> delta; /*! vector of delta for each stacker*/
	std::vector<double> alpha; /*!vector of alpha for each stacker*/
	static double x[NUMPOINTS]; /*! static x vector shared by all cavities*/
	static int totalFreqCavity; /*! integer total number of stacker systems*/
public:
	freqCavity();
	~freqCavity();
	void addStacker(double Tr, double delta, double R, double alpha);
	fftw_complex* getF(){ return F;}
	std::string toString();
	double* getX(){return x;}

};
double freqCavity::x[] = {};
int freqCavity::totalFreqCavity=0;

/**
 * @brief Constructor
 * @details Increments total number of stacker systems. If this is the first stacker system it initializes the x array.
 */
freqCavity::freqCavity(){
	totalFreqCavity++;
	if(totalFreqCavity==1) linspace(0,1,NUMPOINTS,x);
}

/**
 * @brief Destructor
 * @details Removes the instance of the freqCavity. decrements the total number of stacker systems. If this instance of freqCavity has a stacker in it, then it will free the memory allocated for the fftw_complex array. 
 * @return none
 */
freqCavity::~freqCavity(){
	totalFreqCavity--;
	if(totalNumber > 0) fftw_free(F);
}

/**
 * @brief Adds a stacker to the freqCavity instance
 * @details Adds the stacker values to the vectors. If this is the first stacker added then allocate the F array and fill it. If there is already one or more stackers, then multiply the current F array by the new impulse response. 
 * /par
 * /note The F updating loops use pointer indexing for speed. The fftw_complex array is 16-byte alinged and we use a double pointer for it so the pointer needs to be increased by two. 
 * 
 * /todo Fix the slow implementation of the exponential in the F calculation.
 * 
 * @param Tr Round-trip time
 * @param R Power reflectivity
 * @param delta Phase
 * @param alpha loss
 */
void freqCavity::addStacker(double Tr, double R,  double delta, double alpha){
	//std::cout<<"adding cavity"<<std::endl;
	this -> R.push_back(R);
	this -> Tr.push_back(Tr);
	this -> delta.push_back(delta);
	this -> alpha.push_back(alpha);
	totalNumber++;
	cavity.push_back(totalNumber);
	if(totalNumber==1) 	F = fftw_alloc_complex(NUMPOINTS);

	//now update the impulse response function

	//start by getting pointers to heads of all the arrays
	double *realBegin, *imagBegin, *xBegin;
	realBegin = &F[0][0];
	imagBegin = &F[0][1];
	xBegin = &x[0];
	double r=std::sqrt(R);
	double deltaL;
	double test;

	//if this is the first cavity then we have to initialize F
	if(totalNumber==1){
		for(int i=0;i<NUMPOINTS;i++){
			deltaL = 2.0*PI*(Tr*(*xBegin))+delta;
			*realBegin = (2.0*r-cos(deltaL)*(1.0+R))/(1.0+R-2.0*r*cos(deltaL));
			*imagBegin = -1.0*sin(deltaL)*(1.0-R)/(1.0+R-2.0*r*cos(deltaL));
			realBegin+=2;
			imagBegin+=2;
			xBegin++;
		}
	//if it is not the first cavity then we multiply the new F by the old F
	}else{
		double realVal, imagVal;
		std::complex<double> I(0,1);
		std::complex<double> temp, temp2, expVal;

		for(int i=0;i<NUMPOINTS;i++){

			temp2 = {*realBegin, *imagBegin};
			deltaL = 2.0*PI*(Tr*(*xBegin))+delta;
			expVal = std::exp(I*deltaL);
			temp = (r-expVal)/(1.0-r*expVal) * temp2;

			*realBegin = temp.real();
			*imagBegin = temp.imag();
			realBegin+=2;
			imagBegin+=2;
			xBegin++;
		}
	}
}

/**
 * @brief returns a descriptive string
 * @details returns a multi-line string that describes this freqCavity
 * @return a string
 */
std::string freqCavity::toString(){
	using namespace std;
	string str1 = "Total Number of stackers: " + to_string(totalNumber) + "\n";
	for(int i=0;i<totalNumber;i++){
		str1 +=  "stacker " + to_string(i+1) + ": Tr = " + to_string((int)Tr.at(i)) + 
				", delta = "  +to_string(delta.at(i)) + 
				", R = " + to_string(R.at(i)) +
				", alpha = " + to_string(alpha.at(i)) + "\n";
	}
	return str1;
}


//wanna tell me if my code is insane/bad form? it works.




 int main(int argc, char* argv[]){

 	//////////////////////////create cavity system ///////////////////////
 	freqCavity cav;
 	cav.addStacker(1,0.5855,0.0,0.0);
 	//cav.addStacker(1,0.5950,-3.0370,0.0);
    //cav.addStacker(1,0.5711, 2.4895,0.0);
 	//cav.addStacker(1,0.5776,-0.6930,0.0);
 	//cav.addStacker(9,0.5855,-2.2420,0.0);
 	//cav.addStacker(9,0.5950,-3.0370,0.0);
 	//cav.addStacker(9,0.5711, 2.4895,0.0);
 	//cav.addStacker(9,0.5776,-0.6930,0.0);
 	//std::cout << cav.toString()<<std::endl;

 	//////////////////////// allocate and do the fftw ///////////////////////////
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
 						FFTW_ESTIMATE);
 	f = fftw_plan_dft_1d(NUMPOINTS,
 						 fIn, 
 						 fOut,
 						 FFTW_FORWARD,
 						 FFTW_ESTIMATE);

///////////Fill with zeros except 6 time reversed and conjugated spots//////
 	fftw_complex *wind;
 	wind=fftw_alloc_complex(NUMPOINTS);
 	for(int i=0; i<NUMPOINTS;i++){
 		wind[i][0] = 0.0;
 		wind[i][1] = 0.0;
 	}
 	wind[2044][0] = 3132.41;
 	wind[2044][1] = -1.24189e-13;

 	wind[2045][0] = -1699.74;
 	wind[2045][1] = 1.30369;

 	wind[2046][0] = -1300.65;
 	wind[2046][1] = 1.99518;

 	wind[2047][0] = -995.057;
 	wind[2047][1] = 2.2896;

 	wind[2048][0] = -761.124;
 	wind[2048][1] = 2.33511;

 	wind[2049][0] = -582.083;
 	wind[2049][1] = 2.23227;

 	fftShift(wind, fIn, NUMPOINTS);

 	//now the windowed function is in fIn
 	fftw_execute(f); 
 	//now the ft of the windowed function is in fOut.
 	
 	vecMult(cav.getF(), fOut,bIn, NUMPOINTS);

 	//now bIn holds the multipled F(w) and FFT'd windowed function

 	fftw_execute(b);

 	//now bOut holds the IFFT'd multiplied function

 	fftShift(bOut, wind, NUMPOINTS);

 	//now wind holds the aligned IFFT'd multiplied function
 	toFile(wind, NUMPOINTS, "test.txt");

/*

//////////////////Generate forward fft plan, shift and do it /////////////////////

	fftShift(wind, Fin, NUMPOINTS);
 	pp = fftw_plan_dft_1d(NUMPOINTS,
 						  Fin, 
 						  ftWind,
 						  FFTW_FORWARD, 
 						  FFTW_ESTIMATE);

 	fftw_execute(pp);
//////////////////////now multiple F(w) and transformed burst//////////////////


 	vecMult(ftWind, cav.getF(), Fin, NUMPOINTS);

 	fftw_execute(p);

 	fftShift(Fout, aligned, NUMPOINTS);

 	toFile(aligned,NUMPOINTS, "test.txt");
 	/*

  ///////////////////////////convolution testing//////////////////////////////

  	fftw_complex *toConv, *conv;

  	toConv = fftw_alloc_complex(NUMPOINTS);
  	conv =   fftw_alloc_complex(NUMPOINTS*2);

  	for(int i=0;i<NUMPOINTS;i++){
  		toConv[i][0] = 0;
  		toConv[i][1] = 0;
  	}
  	double r = sqrt(0.5);
  	toConv[NUMPOINTS/2-3][0] = r;
  	toConv[NUMPOINTS/2-3][1] = 0.0;

  	toConv[NUMPOINTS/2-2][0] = 0.0;
  	toConv[NUMPOINTS/2-2][1] = -r*r;

  	toConv[NUMPOINTS/2-1][0] = 0.0;
  	toConv[NUMPOINTS/2-1][1] = -r*r*r;

  	toConv[NUMPOINTS/2][0] = 0.0;
  	toConv[NUMPOINTS/2][1] = -r*r*r*r;


  //	convolve(Fin, toConv, conv, NUMPOINTS);
  //	toFile(conv, NUMPOINTS*2, "test.txt");
  */

 	return 0;
 }
