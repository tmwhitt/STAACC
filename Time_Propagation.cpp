#include <vector>
#include <iostream>
#include <cmath>
#include <ctime>

using namespace std;

struct cavity{

	int length;
	double refl;
	double trans;
	double phaseReal;
	double phaseImag;

	cavity(int l, double r, double p) : length(l){
		refl = sqrt(r);
		trans = sqrt(1-r);
		phaseReal = cos(p);
		phaseImag = sin(p);
	}
};

struct pulse{

	double real;
	double imag;

	pulse(double r, double i) : real(r), imag(i) {}
	pulse(double r) : real(r) {}
	pulse() : real(0), imag(0) {}
};

// Size of input pulse train must = size of output pulse train
// Size of input pulse train must also be larger than the number of cavities, cause it would be a mess of if statements otherwise
// Fill it with 0's if you have to for impulse response.
void cavityPropagation(	const vector<cavity> &cavSet, const vector<pulse> &input, 
						vector<pulse> &output){

	int i,j;
	int numCavs = cavSet.size();
	int numPulses = input.size();
	double temp;

	vector< vector<pulse> > Ecav(numCavs, vector<pulse>(numPulses + numCavs));
	vector< vector<pulse> > Eout(numCavs, vector<pulse>(numPulses + numCavs));

	// First pulse hits all the cavities, and loads them up
	Ecav[0][0].real = -input[0].imag * cavSet[0].trans;
	Ecav[0][0].imag =  input[0].real * cavSet[0].trans;
	Eout[0][0].real =  input[0].real * cavSet[0].refl;
	Eout[0][0].imag =  input[0].imag * cavSet[0].refl;

	for(j=1;j<numCavs;++j){	
		Ecav[j][0].real = -Eout[j-1][0].imag * cavSet[j].trans;
		Ecav[j][0].imag =  Eout[j-1][0].real * cavSet[j].trans;
		Eout[j][0].real =  Eout[j-1][0].real * cavSet[j].refl;
		Eout[j][0].imag =  Eout[j-1][0].imag * cavSet[j].refl;
	}

	for(i=1;i<numPulses;++i){
		if(i>=cavSet[0].length){
			Ecav[0][i].real = Ecav[0][i-cavSet[0].length].real * cavSet[0].refl - input[i].imag * cavSet[0].trans;
			Ecav[0][i].imag = Ecav[0][i-cavSet[0].length].imag * cavSet[0].refl + input[i].real * cavSet[0].trans;
			Eout[0][i].real = input[i].real * cavSet[0].refl - Ecav[0][i-cavSet[0].length].imag * cavSet[0].trans;
			Eout[0][i].imag = input[i].imag * cavSet[0].refl + Ecav[0][i-cavSet[0].length].real * cavSet[0].trans;

			temp = Ecav[0][i].real*cavSet[0].phaseReal - Ecav[0][i].imag*cavSet[0].phaseImag;
			Ecav[0][i].imag = Ecav[0][i].real*cavSet[0].phaseImag + Ecav[0][i].imag*cavSet[0].phaseReal;
			Ecav[0][i].real = temp;
		}
		else{
			Ecav[0][i].real = -input[i].imag * cavSet[0].trans;
			Ecav[0][i].imag =  input[i].real * cavSet[0].trans;
			Eout[0][i].real =  input[i].real * cavSet[0].refl;
			Eout[0][i].imag =  input[i].imag * cavSet[0].refl;
		}
		for(j=1;j<numCavs;++j){
			if(i>=cavSet[j].length){
				Ecav[j][i].real = Ecav[j][i-cavSet[j].length].real * cavSet[j].refl - Eout[j-1][i].imag * cavSet[j].trans;
				Ecav[j][i].imag = Ecav[j][i-cavSet[j].length].imag * cavSet[j].refl + Eout[j-1][i].real * cavSet[j].trans;
				Eout[j][i].real = Eout[j-1][i].real * cavSet[j].refl - Ecav[j][i-cavSet[j].length].imag * cavSet[j].trans;
				Eout[j][i].imag = Eout[j-1][i].imag * cavSet[j].refl + Ecav[j][i-cavSet[j].length].real * cavSet[j].trans;

				temp = Ecav[j][i].real*cavSet[j].phaseReal - Ecav[j][i].imag*cavSet[j].phaseImag;
				Ecav[j][i].imag = Ecav[j][i].real*cavSet[j].phaseImag + Ecav[j][i].imag*cavSet[j].phaseReal;
				Ecav[j][i].real = temp;
			}
			else{
				Ecav[j][i].real = -Eout[j-1][i].imag * cavSet[j].trans;
				Ecav[j][i].imag =  Eout[j-1][i].real * cavSet[j].trans;
				Eout[j][i].real =  Eout[j-1][i].real * cavSet[j].refl;
				Eout[j][i].imag =  Eout[j-1][i].imag * cavSet[j].refl;
			}
		}
	}

	output.swap(Eout[numCavs-1]);
	return;
}

#define PI 3.14159265359

int main(){
	clock_t start;
	double duration;
	vector<pulse> input(1000000);
	vector<pulse> output;
	double temp;

	input[0] = pulse(1);

	vector<cavity> cavSet;
	cavSet.push_back(cavity(1,0.5,0));
	cavSet.push_back(cavity(1,0.5,0));
	cavSet.push_back(cavity(1,0.5,0));
	cavSet.push_back(cavity(1,0.5,0));
	cavSet.push_back(cavity(9,0.5,0));
	cavSet.push_back(cavity(9,0.5,0));
	cavSet.push_back(cavity(9,0.5,0));
	cavSet.push_back(cavity(9,0.5,0));

	start = clock();
	cavityPropagation(cavSet,input,output);
	duration = ( clock()-start ) / (double)(CLOCKS_PER_SEC);
	cout << "Code took " << duration << " seconds." << endl;

	/*
	cout << "Output train:" << endl;
	for(int i=0;i<output.size();++i){
		temp = output[i].real*output[i].real + output[i].imag*output[i].imag;
		cout << temp << endl;
	}
	*/
	return 0;
}