#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

struct cavity{

	int length;
	double reflectivity;
	double phase;

	cavity(int l, double r, double p) : length(l), reflectivity(r), phase(p) {}
};

struct pulse{

	double real;
	double imag;

	pulse(double r, double i) : real(r), imag(i) {}
	pulse(double r) : real(r) {}
	pulse() : real(0), imag(0) {}
};

// Size of input pulse train must be <= size of output pulse train
void cavityPropagation(	const vector<cavity> &cavSet, const vector<pulse> &input, 
						vector<pulse> &output){

	int i,j;
	int numCavs = cavSet.size();
	int numPulsesIn = input.size();
	int numPulsesOut = output.size();

	vector<vector<pulse>> Ecav(numCav, vector<pulse>(numPulsesOut + numCavs));
	vector<vector<pulse>> Eout(numCav, vector<pulse>(numPulsesOut + numCavs));

	for(i=0;i<numCav;++i){
		
	}
}

int main(){
	cout << "Hello World" << endl;
	return 0;
}