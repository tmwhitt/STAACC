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
 * @details Takes output fftw_complex array and aligns it from negative to positive frequency. Written using pointers for speed. FFTW is 16-byte aligned and requires +=2 to go to the next element in the array. 
 * 
 * @param a Input array to be aligned
 * @param b Output aligned array
 * @param Length of array
 */
void reAlign(fftw_complex *a, fftw_complex *b, int N){
	double *realIn, *imagIn, *imagOut, *realOut;
	
	realIn = &a[N/2][0];
	imagIn = &a[N/2][1];
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
	realOut = &b[N/2][0];
	imagOut = &b[N/2][1];
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


class freqCavity{
private:
	fftw_complex *F;
	int totalNumber =0;
	std::vector<int> cavity;
	std::vector<double> R;
	std::vector<double> Tr;
	std::vector<double> delta;
	std::vector<double> alpha;
	static double x[NUMPOINTS];
public:
	freqCavity();
	~freqCavity();
	void addCavity(double Tr, double delta, double R, double alpha);
	fftw_complex* getF(){ return F;}
	std::string toString();
	double* getX(){return x;}

};
double freqCavity::x[] = {};

freqCavity::freqCavity(){

}

freqCavity::~freqCavity(){
	R.pop_back();
	Tr.pop_back();
	delta.pop_back();
	alpha.pop_back();
	cavity.pop_back();
	totalNumber--;
	if(totalNumber==0) fftw_free(F);
}

void freqCavity::addCavity(double Tr, double R,  double delta, double alpha){
	this -> R.push_back(R);
	this -> Tr.push_back(Tr);
	this -> delta.push_back(delta);
	this -> alpha.push_back(alpha);
	totalNumber++;
	cavity.push_back(totalNumber);
	if(totalNumber==1){
		//std::cout << "First Cavity, allocating F" <<std::endl;
		linspace(0,1,NUMPOINTS, x);
		F = fftw_alloc_complex(NUMPOINTS);
	} 

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

		std::cout<< "n=1" << std::endl;

		for(int i=0;i<NUMPOINTS;i++){
			deltaL = 2.0*PI*(Tr*(*xBegin))+delta;
			*realBegin = (2.0*r-cos(deltaL)*(1.0+R))/(1.0+R-2.0*r*cos(deltaL));
			*imagBegin = -sin(deltaL)*(1.0-R)/(1.0+R-2.0*r*cos(deltaL));
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

 	freqCavity cav;
 	cav.addCavity(1,0.5855,-2.2420,0.0);
 	cav.addCavity(1,0.5950,-3.0370,0.0);
 	cav.addCavity(1,0.5711, 2.4895,0.0);
 	cav.addCavity(1,0.5776,-0.6930,0.0);
 	cav.addCavity(9,0.5855,-2.2420,0.0);
 	cav.addCavity(9,0.5950,-3.0370,0.0);
 	cav.addCavity(9,0.5711, 2.4895,0.0);
 	cav.addCavity(9,0.5776,-0.6930,0.0);
 	std::cout << cav.toString()<<std::endl;

  	fftw_complex *IR, *alignedIR, *Fin;
 	IR = fftw_alloc_complex(NUMPOINTS);
 	alignedIR = fftw_alloc_complex(NUMPOINTS);
 	Fin = cav.getF();


 	//allocate the fftw plan and calculate it 
 	fftw_plan p;
 	p = fftw_plan_dft_1d(NUMPOINTS,
 						Fin,
 						IR, 
 						FFTW_BACKWARD, 
 						FFTW_ESTIMATE);

 	fftw_execute(p);
 	reAlign(IR, alignedIR, NUMPOINTS);
  	toFile(alignedIR, NUMPOINTS, "test.txt");

 	fftw_free(IR);
 	fftw_free(alignedIR);

 	return 0;
 }
