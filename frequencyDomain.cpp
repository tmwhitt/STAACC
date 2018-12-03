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
/**
 * @brief Mutliplies two fftw_complex arrays
 * @details Multiplies one complex fftw array by the other and stores in another array. All arrays must be the same length. 
 * 
 * @param a first array to be multiplied
 * @param b second array to be multiplied
 * @param c output array
 * @param N length of the arrays
 */
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

/**
 * @brief Complex conjuages of a fftw_complex array
 * @details [long description]
 * 
 * @param a array to be conjuagted
 * @param N length of array.
 */
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
	bool gottenF = false;
public:
	freqCavity();
	~freqCavity();
	void addStacker(double Tr, double delta, double R, double alpha);
	void fillF(fftw_complex *a, int n);
	fftw_complex* getF();
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
	if(gottenF ==true) fftw_free(F);
}

/**
 * @brief Fills input array with transfer function of cavity System
 * @details Takes input array a. If the freqCavity instance already has calculated F for itself, then just copy that into a. If not, we have to calculate and fill a. 
 * 
 * @param a pointer to head of fftw_complex array.
 * @param n length of fftw_complex array a
 */
void freqCavity::fillF(fftw_complex *a, int n){
	double *aReal, *aImag, *xBegin;
	aReal = &a[0][0];
	aImag = &a[0][1];
	std::cout<<aReal<<std::endl;

	//if the instance of the class already has calculated the array then just fill a from that array.
	if(gottenF==true){
		double *bReal,*bImag;
		bReal = &F[0][0];
		bImag = &F[0][1];
		for(int i=0;i<n;i++){
			*aReal = *bReal;
			*aImag = *bImag;
			bReal +=2;
			bImag +=2;
			aReal +=2;
			aImag +=2;
		}
	//otherwise calculate and fill a
	}else{
		xBegin = &x[0];
		double r = std::sqrt(R[0]);

		//for the first stacker just initialize the array values	
		for(int i=0;i<n;i++){
			double deltaL = 2.0*PI*(Tr[0]*(*xBegin))+delta[0];
			*aReal = (2.0*r-cos(deltaL)*(1.0+R[0]))/(1.0+R[0]-2.0*r*cos(deltaL));
			*aImag = -1.0*sin(deltaL)*(1.0-R[0])/(1.0+R[0]-2.0*r*cos(deltaL));
	
			xBegin++;
			aReal+=2;
			aImag+=2;
		}

		//then for each subsequent stacker multiply the current array value by the value for this stacker
		for(int i=1; i<totalNumber;i++){
			aReal = &a[0][0];
			aImag = &a[0][1];
			xBegin = &x[0];
			double r = std::sqrt(R[i]);
			for(int j=0;j<n;j++){

				double deltaL = 2.0*PI*(Tr[i]*(*xBegin))+delta[i];
				double realTemp = (2.0*r-cos(deltaL)*(1.0+R[i]))/(1.0+R[i]-2.0*r*cos(deltaL));
				double imagTemp =  -1.0*sin(deltaL)*(1.0-R[i])/(1.0+R[i]-2.0*r*cos(deltaL));
				double aRealTemp = *aReal;
				double aImagTemp = *aImag;

				*aReal = aRealTemp*realTemp - aImagTemp*imagTemp;
				*aImag = aRealTemp*imagTemp + aImagTemp*realTemp;

				xBegin++;
				aReal+=2;
				aImag+=2;
			}
		}
	}
}


void freqCavity::addStacker(double Tr, double R,  double delta, double alpha){
	this -> R.push_back(R);
	this -> Tr.push_back(Tr);
	this -> delta.push_back(delta);
	this -> alpha.push_back(alpha);
	totalNumber++;
	cavity.push_back(totalNumber);
}

fftw_complex* freqCavity::getF(){
	if(gottenF==true) return F;
	else{
		F = fftw_alloc_complex(NUMPOINTS);
		fillF(F, NUMPOINTS);
		gottenF= true;
		return F;
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
 	cav.addStacker(1,0.5855,-2.2420,0.0);
 	cav.addStacker(1,0.5950,-3.0370,0.0);
    cav.addStacker(1,0.5711, 2.4895,0.0);
 	cav.addStacker(1,0.5776,-0.6930,0.0);
 	cav.addStacker(9,0.5855,-2.2420,0.0);
 	cav.addStacker(9,0.5950,-3.0370,0.0);
 	cav.addStacker(9,0.5711, 2.4895,0.0);
 	cav.addStacker(9,0.5776,-0.6930,0.0);
 	std::cout << cav.toString()<<std::endl;

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

//test impulse response and the new fillF() function
 	cav.fillF(bIn,NUMPOINTS);
 	fftw_execute(b);
 	fftShift(bOut, bIn, NUMPOINTS);
 	toFile(bOut, NUMPOINTS, "test.txt");



 	/*

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
*/
 //	fftw_free(wind);
 	fftw_free(fIn);
 	fftw_free(fOut);
 	//fftw_free(bIn);
 	fftw_free(bOut);
 	return 0;
 }
