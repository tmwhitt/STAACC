//! Freq Domain test code
/*! \brief This code runs the frequency domain propagation of a stacker
 * 
 * \author Alex "A-drizzle"  Rainville
 * \date 11/10/2018
 * \version 1.1 modified to use fftw malloc
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
 #define NUMPOINTS 16384



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
			myfile << std::to_string(a[0][i]) << "," << std::to_string(a[1][i]) << "\n";;
		}

	}else std::cout << "cannot open file"<< std::endl;


	myfile.close();
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
	fftw_free(toReturn);
}

fftw_complex* stacker::getF(){
		toReturn  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NUMPOINTS);
		for(int i=0;i<NUMPOINTS;i++){
			toReturn[0][i] = real(F[i]);
			toReturn[1][i] = imag(F[i]);
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
			deltaL = 2*M_PI*(x[i] +delta);
			expVal = std::exp(I*deltaL);
			F[i] = (r-expVal)/(1.0-r*expVal);
		}
	} 
	else{
		//std::cout << "called a n>1 stacker, updating overall F" <<std::endl;
		for(int i=0;i<NUMPOINTS;i++){
			deltaL = 2*M_PI*(x[i] +delta);
			expVal = std::exp(I*deltaL);
			F[i] = (r-expVal)/(1.0-r*expVal);
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
 	std::cout<<eighth.toString() <<std::endl;
 	std::cout<<first.toString() <<std::endl;

 	//std::cout << first.getTr() <<std::endl;
 	//first.setTr(0.50);

 	std::complex<double> out[NUMPOINTS] = {};
 	fftw_plan p;
 	p = fftw_plan_dft_1d(NUMPOINTS,
 						reinterpret_cast<fftw_complex*>(eighth.getF()), 
 						reinterpret_cast<fftw_complex*>(out), 
 						FFTW_BACKWARD, 
 						FFTW_ESTIMATE);
 	fftw_execute(p);
 	fftw_print_plan(p) ;
 	fftw_destroy_plan(p);

 	/*
 	for(int i=0;i<5;i++){
 		std::cout<< out[i] <<std::endl;
 	}
 	*/
 	//toFile(out, NUMPOINTS, "test.txt");


 	return 0;
 }
