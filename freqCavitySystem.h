
 #include <math.h>
 #include <string>
 #include <complex>
 #include <fftw3.h>
 #include <vector>
 #include "fftwComplexMath.h"
//number of points in FFTW vector. Needs to be power of two.
 #define NUMPOINTS 16384
 #define SINELENGTH 4096
 #define PI  3.141592653589793



/**
 * @brief Frequency Domain cavity
 * @details A cavity object that describes a system of stackers. 
 * 
 */
class freqCavitySystem{
private:
	fftw_complex *F; /*! Impulse response of the system */
	int totalNumber =0; /*! Total number of individual stackers inside */
	std::vector<int> cavity; /*! vector of ints representing each stacker*/
	std::vector<double> R; /*! vector of R for each stacker*/
	std::vector<double> Tr; /*!vector of Tr for each stacker*/
	std::vector<double> delta; /*! vector of delta for each stacker*/
	std::vector<double> alpha; /*!vector of alpha for each stacker*/
	bool gottenF = false;
	static int totalfreqCavitySystem; /*! integer total number of stacker systems*/
	static double x[NUMPOINTS]; /*! static x vector shared by all cavities*/
	static double sineLUT[SINELENGTH];
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
 * @brief Constructor
 * @details Increments total number of stacker systems. If this is the first stacker system it initializes the x array.
 */
freqCavitySystem::freqCavitySystem(){
	totalfreqCavitySystem++;
	if(totalfreqCavitySystem==1){
		//std::cout<<" WRITING X AND SINE TABLE" <<std::endl;
		linspace(0,1,NUMPOINTS,x);
		double val = 0;
		double step = 2*PI/(SINELENGTH-1.0);
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
 * @details Removes the instance of the freqCavitySystem. decrements the total number of stacker systems. If this instance of freqCavitySystem has a stacker in it, then it will free the memory allocated for the fftw_complex array. 
 * @return none
 */
freqCavitySystem::~freqCavitySystem(){
	totalfreqCavitySystem--;
	if(gottenF ==true) fftw_free(F);
}

/**
 * @brief Fills input array with transfer function of cavity System
 * @details Takes input array a. If the freqCavitySystem instance already has calculated F for itself, then just copy that into a. If not, we have to calculate and fill a. 
 * 
 * @param a pointer to head of fftw_complex array.
 * @param n length of fftw_complex array a
 */
 //with built in cosine
 
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
			*aReal = (2.0*r-cos(deltaL)*(1.0+R[0]))/(1.0+R[0]-2.0*r*sin(deltaL));
			*aImag = -1.0*cos(deltaL)*(1.0-R[0])/(1.0+R[0]-2.0*r*cos(deltaL));
	
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

//with my own cosine 
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
		double deltaL;	
		for(int i=0;i<n;i++){
			deltaL = 2.8*PI*Tr[0]*(*xBegin)+delta[0];
			*aReal = (2.0*r-bCosine(deltaL)*(1.0+R[0]))/(1.0+R[0]-2.0*r*bSine(deltaL));
			*aImag = -1.0*bCosine(deltaL)*(1.0-R[0])/(1.0+R[0]-2.0*r*bCosine(deltaL));
	
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
			aReal = &a[0][0];
			aImag = &a[0][1];
			xBegin = &x[0];

			double r = std::sqrt(R[i]);
			double deltaL,realTemp,imagTemp,aRealTemp,aImagTemp;
			for(int j=0;j<n;j++){
				deltaL = 2.0*PI*(*trBegin)*(*xBegin)+(*deltaBegin);
				realTemp = (2.0*r-bCosine(deltaL)*(1.0+(*rBegin))/(1.0+(*rBegin)-2.0*r*bCosine(deltaL)));
				imagTemp =  -1.0*bSine(deltaL)*(1.0-(*rBegin))/(1.0+(*rBegin)-2.0*r*bCosine(deltaL));
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


void freqCavitySystem::addCavity(double Tr, double R,  double delta, double alpha){
	this -> R.push_back(R);
	this -> Tr.push_back(Tr);
	this -> delta.push_back(delta);
	this -> alpha.push_back(alpha);
	totalNumber++;
	cavity.push_back(totalNumber);
}

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
 * @details returns a multi-line string that describes this freqCavitySystem
 * @return a string
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

double freqCavitySystem::bSine(double x){
    int val1 = ((int) ((x/(2*PI))*SINELENGTH*4))& (4*SINELENGTH-1);
    int val2 =  (val1) & (SINELENGTH-1);
    
    if(!(val2^(SINELENGTH-1))) val2--;

    if       (val1 < SINELENGTH) return sineLUT[val2&(SINELENGTH-1)]; 
    else if(val1 < 2*SINELENGTH) return sineLUT[SINELENGTH-1-val2&(SINELENGTH-1)];
    else if(val1 < 3*SINELENGTH) return -1.0*sineLUT[val2&(SINELENGTH-1)];
    else if(val1 < 4*SINELENGTH) return -1.0*sineLUT[SINELENGTH-1-val2&(SINELENGTH-1)];
    
}
double freqCavitySystem::bCosine(double x){
	return bSine(PI/2-x);
}