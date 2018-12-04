 #include <cstdio>
 #include <iostream>
 #include <math.h>
 #include <string>
 #include <fftw3.h>
 #include <fstream>
//number of points in FFTW vector. Needs to be power of two.



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
	for(double val=a;val<=b;val+=step){
		*x = val;
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