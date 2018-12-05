#include <math.h>
 #include <string>
 #include <complex>
 #include <fftw3.h>
 #include <vector>
 #include "fftwComplexMath.h"

/**
 * \def NUMPOINTS
 * Number of points used for the FFTW vectors
 */

 /**
 * \def PI
 * What it looks like
 */

 /**
 * \def SINELENGTH
 * Number of points to use for the sine LUT
 */

 #define NUMPOINTS 16384 
 #define PI  3.141592653589793 
 #define SINELENGTH 4096 



/**
 * \brief Frequency Domain cavity
 * \details A class that descibes a system of stacking cavities
 * 
 */
class freqCavitySystem{
private:
	fftw_complex *F;                       /**< \brief Local copy of the impulse response vector of the system.  */
	int totalNumber =0;                    /**< \brief Total number of individual stackers inside */
	std::vector<int> cavity;               /**< \brief vector of ints representing each stacker*/
	std::vector<double> R;                 /**< \brief vector of R for each stacker*/
	std::vector<double> Tr;                /**< \brief vector of Tr for each stacker*/
	std::vector<double> delta;             /**< \brief vector of delta for each stacker*/
	std::vector<double> alpha;             /**< \brief vector of alpha for each stacker*/
	bool gottenF = false;                  /**< \brief boolean for if the stacker's F has been allocated and filled*/ 
	static int totalfreqCavitySystem;      /**< \brief integer total number of stacker systems*/
	static double x[NUMPOINTS];            /**< \brief static x vector shared by all cavities*/
	static double sineLUT[SINELENGTH];     /**< \brief static sine lookup table vector shared by all cavities */
public:
	freqCavitySystem();
	~freqCavitySystem();
	void addCavity(double Tr, double delta, double R, double alpha);
	void fillF(fftw_complex *a, int n);
	void fillFExact(fftw_complex *a, int n);
	fftw_complex* getF();
	std::string toString();
	double* getX(){return x;}
	double bSine(double x);
	double bCosine(double x);

};
double freqCavitySystem::x[] = {};
double freqCavitySystem::sineLUT[] = {};
int freqCavitySystem::totalfreqCavitySystem=0;

/**
 * @brief Default Constructor
 * @details Increments total number of stacker systems. If this is the first stacker system it initializes the x array and the sineLUT array. 
 */
freqCavitySystem::freqCavitySystem(){
	totalfreqCavitySystem++;
	if(totalfreqCavitySystem==1){
		//std::cout<<" WRITING X AND SINE TABLE" <<std::endl;
		linspace(0,1,NUMPOINTS,x);
		double val = 0;
		double step = PI/(2.0*(SINELENGTH-1.0));
		double *sineLUTStart;
		sineLUTStart = &sineLUT[0];
		for(int i=0;i<SINELENGTH;i++){
			*sineLUTStart = sin(val);
			sineLUTStart++;
			val += step; 
		}
	}
}

/**
 * @brief Destructor
 * @details Removes the instance of the freqCavitySystem. Decrements the total number of freqCavitySystems. If this instance of freqCavitySystem has allocated space for F, then it will clear it. 
 */
freqCavitySystem::~freqCavitySystem(){
	totalfreqCavitySystem--;
	if(gottenF ==true) fftw_free(F);
}

/**
 * @brief Fills input array with transfer function of cavity System using built-in sine and cosine functions. 
 * @details Fills array a with impulse response of the current stacker system. If the stacker system has initialized and filled F then this will copy that array into a. If there is no local copy of F, then this function fills a with the transfer function of the cavity system. 
 * 
 * @param a fftw_complex array.
 * @param n length of array. 
 */
void freqCavitySystem::fillFExact(fftw_complex *a, int n){
	double *aReal, *aImag, *xBegin;
	aReal = &a[0][0];
	aImag = &a[0][1];

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

/**
 * @brief Fills input array with transfer function of cavity System using sine and cosine lookup tables. 
 * @details Fills array a with impulse response of the current stacker system. If the stacker system has initialized and filled F then this will copy that array into a. If there is no local copy of F, then this function fills a with the transfer function of the cavity system. 
 * 
 * \note This utilizes the lookup tables.  
 * 
 * @param a fftw_complex array.
 * @param n length of array. 
 */
void freqCavitySystem::fillF(fftw_complex *a, int n){
	double *aReal, *aImag, *xBegin;
	aReal = &a[0][0];
	aImag = &a[0][1];

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
			*aReal = (2.0*r-bCosine(deltaL)*(1.0+R[0]))/(1.0+R[0]-2.0*r*bCosine(deltaL));
			*aImag = -1.0*bSine(deltaL)*(1.0-R[0])/(1.0+R[0]-2.0*r*bCosine(deltaL));
	
			xBegin++;
			aReal+=2;
			aImag+=2;
		}


		//then for each subsequent stacker multiply the current array value by the value for this stacker
		double *trBegin, *rBegin, *deltaBegin;
		trBegin = Tr.data();
		rBegin  = R.data();
		deltaBegin = delta.data();
		for(int i=1; i<totalNumber;i++){
			std::cout<<"NO"<<std::endl;
			aReal = &a[0][0];
			aImag = &a[0][1];
			xBegin = &x[0];

			double r = std::sqrt(R[i]);
			double deltaL,realTemp,imagTemp,aRealTemp,aImagTemp;
			for(int j=0;j<n;j++){
				deltaL = 2.0*PI*(*trBegin)*(*xBegin)+(*deltaBegin);
				realTemp = (2.0*r-cos(deltaL)*(1.0+(*rBegin))/(1.0+(*rBegin)-2.0*r*cos(deltaL)));
				imagTemp =  -1.0*bSine(deltaL)*(1.0-(*rBegin))/(1.0+(*rBegin)-2.0*r*cos(deltaL));
				aRealTemp = *aReal;
				aImagTemp = *aImag;

				*aReal = aRealTemp*realTemp - aImagTemp*imagTemp;
				*aImag = aRealTemp*imagTemp + aImagTemp*realTemp;

				xBegin++;
				aReal+=2;
				aImag+=2;
			}
			trBegin++;
			rBegin ++;
			deltaBegin++;
		}
	}
}

/**
 * @brief Adds a cavity to the cavity system. 
 * @details Increments the total number of cavities. Pushes R, Tr, delta, and alpha back on their respective vectors. Pushes back totalNumber on that vector.
 * 
 * @param Tr [description]
 * @param R [description]
 * @param delta [description]
 * @param alpha [description]
 */
void freqCavitySystem::addCavity(double Tr, double R,  double delta, double alpha){
	this -> R.push_back(R);
	this -> Tr.push_back(Tr);
	this -> delta.push_back(delta);
	this -> alpha.push_back(alpha);
	totalNumber++;
	cavity.push_back(totalNumber);
}

/**
 * @brief Returns (and possibly initializes and calculates) the transfer function of the cavity system. 
 * @details If the cavity system already has an initialized copy of F (given by gottenF==true) then the function returns that. Otherwise, it allocates F and calls fillF() to fill it. THen it sets gotten F true and returns F.
 * @return Local fftw_complex array of transfer function.
 * 
 * \warning The deconstructor will clear the memory of F if it is initialized. Clearing it from somewhere else will lead to a segfault and massive sadness. 
 */
fftw_complex* freqCavitySystem::getF(){
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
 * @details Returns a string describing the total Number of stackers in the system and the parameters of each. 
 * @return A multi-line string
 */
std::string freqCavitySystem::toString(){
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

/**
 * @brief Reutrns the sine(x) using lookup tables. 
 * @details Uses bit-manipulation to return the sine(x) from a lookup table. First folds x into 2pi, then folds into pi/2. Then, depending on where x falls, returns the proper value from the LUT. 
 * 
 * \note error is lower here than in a corresponding 4x length table of the whole sine function and I am not sure why. 
 * \bug The first if statment is used to keep from overrunning the length of the lookup table. The effect of this is to return the value at the previous index. If you happen to lookup and x that fits perfectly at the end of the array, then you will get more error than if you lookup any other. The value of sine pretty constant here so it may be negligable. 
 * 
 * @param x Value to be sine-d
 * @return sine(x)
 */
double freqCavitySystem::bSine(double x){
    int val1 = ((int) ((x/(2*PI))*SINELENGTH*4))& (4*SINELENGTH-1);
    int val2 =  (val1) & (SINELENGTH-1);
    
    if(!(val2^(SINELENGTH-1))) val2--;

    if       (val1 < SINELENGTH) return sineLUT[val2&(SINELENGTH-1)]; 
    else if(val1 < 2*SINELENGTH) return sineLUT[SINELENGTH-1-val2&(SINELENGTH-1)];
    else if(val1 < 3*SINELENGTH) return -1.0*sineLUT[val2&(SINELENGTH-1)];
    else if(val1 < 4*SINELENGTH) return -1.0*sineLUT[SINELENGTH-1-val2&(SINELENGTH-1)];
    
}

/**
 * @brief Reutrns the cos(x) using lookup tables. 
 * @details Calls bSine(PI/2-x)
 * @param x Value to be cosine-d
 * @return cos(x)
 */
double freqCavitySystem::bCosine(double x){
	return bSine(PI/2-x);
}
/*
double freqCavitySystem::bSine(double x){
    int index = ((int) ((x/(2*PI))*SINELENGTH))&(SINELENGTH-1);
    return sineLUT[index];
}
*/