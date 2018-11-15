//! Freq Domain test code
/*! \brief This code runs the frequency domain propagation of a stacker
 * 
 * \author Alex "A-drizzle"  Rainville
 * \date 11/10/2018
 * \version 1.2 fixed fftw_complex indexing and transfer function algorithm
 * \note All units mks unless noted in program
 * \par
 * loading module on flux: -bash-4.2$ module load fftw/3.3.4/gcc/4.8.5
 * \par
 * to compile on flux: -bash-4.2$ g++ -I$FFTW_INCLUDE -L$FFTW_LIB -lfftw3 -lm <filename.cpp> -o <output file name>
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
 #define NUMPOINTS 1024



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
 * @details [long description]
 * 
 * @param a [description]
 * @param b [description]
 * @param c [description]
 * @param N [description]
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

void reAlign(fftw_complex *a, fftw_complex *b, int N){
	for(int i=0;i<N/2;i++){
		b[i][0] = a[i+N/2][0];
		b[i][1] = a[i+N/2][1];
	}
	for(int i=N/2;i<N;i++){
		b[i][0] = a[i-N/2][0];
		b[i][1] = a[i-N/2][1];
	}
}

/*
double* getFreq(int NUMPOINTS, int T){
	double w[NUMPOINTS] = {};
	for(int i=0;i<NUMPOINTS/2;i++){
		w[i] = i/T;
	}
	for(int i =NUMPOINTS/2;i<NUMPOINTS;i++){
		w[i] = 
	}
}
*/

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
		//static double tabulatedCos[NUMPOINTS];
		//static double tabulatedSin[NUMPOINTS];
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
		fftw_free(toReturn);
		fftw_free(conjToReturn);
	}
}

fftw_complex* stacker::getF(){
		toReturn  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NUMPOINTS);
		for(int i=0;i<NUMPOINTS;i++){
			toReturn[i][0] = real(F[i]);
			toReturn[i][1] = imag(F[i]);
		}
	return toReturn;
}

fftw_complex* stacker::getConjF(){
		toReturn  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NUMPOINTS);
		for(int i=0;i<NUMPOINTS;i++){
			toReturn[i][0] = -1.0*imag(F[i]);
			toReturn[i][1] = real(F[i]);
		}
	return toReturn;
}



/**
 * @brief Calculates overall transfer function
 * @details If the stacker is the first one, then memory copys F into overallF. 
 * Else it will multiply the current overallF by the current F to update. 
 */
void stacker::updateF(){
	std::complex<double> r(std::sqrt(R),0.0);
	std::complex<double> I(0.0,1.0);
	std::complex<double> deltaL;
	std::complex<double> expVal;

	if(stackerID==1){
		//std::cout << "called first stacker, assigning values to overallF" <<std::endl;

		for(int i=0;i<NUMPOINTS;i++){
			deltaL = 2*M_PI*Tr*(x[i] +delta);
			expVal = std::exp(I*deltaL);
			F[i] = (r-expVal)/(1.0-r*expVal);
		}
	} 
	else{
		//std::cout << "called a n>1 stacker, updating overall F" <<std::endl;
		for(int i=0;i<NUMPOINTS;i++){
			deltaL = 2*M_PI*Tr*(x[i] +delta);
			expVal = std::exp(I*deltaL);
			F[i] = (r-expVal)/(1.0-r*expVal)*F[i];
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




 int main(){
 	//stcker takes in Tr, delta, R, alpha
 	double Tr1 = 1.0;
 	double Tr2 = 9.0;
 	double delta = 0.0;
 	double R = 0.5;
 	double alpha = 0.0;
 	stacker first (Tr1, delta,R,alpha);
 	std::cout<< first.toString() <<std::endl;
  	stacker second (Tr1, delta,R,alpha);
 	stacker third (Tr1, delta,R,alpha);
 	stacker fourth (Tr1, delta,R,alpha);
 	stacker fifth (Tr2, delta,R,alpha);
 	stacker sixth (Tr2, delta,R,alpha);
 	stacker seventh (Tr2, delta,R,alpha);
 	stacker eighth (Tr2, delta,R,alpha);
 	//std::cout<<eighth.toString() <<std::endl;
 	//std::cout<<first.toString() <<std::endl;

 	//std::cout << first.getTr() <<std::endl;
 	//first.setTr(0.50);

 	fftw_complex *out, *alignedOut;
 	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NUMPOINTS);
 	alignedOut = fftw_alloc_complex(NUMPOINTS);
 	fftw_plan p;
 	p = fftw_plan_dft_1d(NUMPOINTS,
 						eighth.getConjF(),
 						out, 
 						FFTW_BACKWARD, 
 						FFTW_ESTIMATE);
 	fftw_execute(p);
 	//fftw_print_plan(p) ;
 
 	std::cout<<std::endl;

 	/*
 	for(int i=0;i<5;i++){
 		std::cout<< out[i] <<std::endl;
 	}
 	*/
 	reAlign(out, alignedOut, NUMPOINTS);
 	toFile(alignedOut, NUMPOINTS, "test.txt");

 	//fftw_complex

 	//create a imput fftw complex that is an impulse at 0
 //	fftw_complex *impulse;
 	//fftw_complex *ftImpulse;
 	//impulse = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NUMPOINTS);
 //	impulse = fftw_alloc_complex(NUMPOINTS);


 	//ftImpulse = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NUMPOINTS*2);


// 	for(int i=0;i<NUMPOINTS;i++){
// 		impulse[i][0] = 0;
// 		impulse[i][1] = 0;
 //	}
 //	impulse[NUMPOINTS/2-1][0] = 1;
 	//impulse[1][NUMPOINTS/2-3] = 0;
 	//std::cout<< &out[NUMPOINTS/2-3][0] <<std::endl;
 	//std::cout<< &out[NUMPOINTS/2-3][1] <<std::endl;
 	//std::cout<< &out[NUMPOINTS/2-1][0] <<std::endl;
 	//impulse[0][NUMPOINTS/2 +1] = 1;
 	//impulse[0][NUMPOINTS/2 +2] = 1;
 //	fftw_complex *convTest;
 	//convTest = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*NUMPOINTS);
 //	convolve(impulse, impulse, convTest, NUMPOINTS);
 //	toFile(convTest,NUMPOINTS, "test.txt");

/*
 	//FT the impulse response
 	fftw_plan pp;
 	pp = fftw_plan_dft_1d(NUMPOINTS,
 					impulse,
 					ftImpulse, 
 					FFTW_FORWARD, 
 					FFTW_ESTIMATE); 
 	fftw_execute(pp);
 	toFile(ftImpulse,NUMPOINTS, "test.txt");
 	std::cout<<"If we made it here its fine" <<std::endl;
 	*/
/*

//multiply the FTd impulse response by the transfer fucntion
 	fftw_complex *f;
 	f = first.getF();

 	double a,b,c,d;
 	for(int i=0;i<NUMPOINTS;i++){
 		a = f[0][i];
 		b = f[1][i];
 		c = ftImpulse[0][i];
 		d = ftImpulse[1][i];
 		ftImpulse[0][i] = a*b-c*d;
 		ftImpulse[1][i] = d*a+b*d;
 	}
//FT back

 	fftw_plan ppp;
 	ppp = fftw_plan_dft_1d(NUMPOINTS,
 					ftImpulse,
 					impulse, 
 					FFTW_BACKWARD, 
 					FFTW_ESTIMATE); 
 	fftw_execute(ppp);
 	toFile(impulse,NUMPOINTS,"test.txt");
 	std::cout<<"If we made it here its fine" <<std::endl;


 	//fftw_destroy_plan(p);
 	//fftw_free(out);

*/
 	return 0;
 	//fftw_free(convTest);
 	//fftw_free(impulse);
 //	fftw_free(ftImpulse);
 	fftw_free(out);
 	fftw_free(alignedOut);
 }
