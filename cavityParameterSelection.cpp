#include <vector>
#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <random>
#include "timePropagation.h"

using namespace std;

#define R_MIN 0.25
#define R_MAX 0.75
#define PI 3.141592653589793238463
#define DIV_PER_DIM 2
#define NUM_PULSES 9

void generateLandscape(	int numR, int numP, const int &numCavs, const int &numDivs, vector< vector<cavity> > &landscape, 
						const vector< vector<double> > &reflVals, const vector< vector<double> > &phaseVals, 
						const vector<int> lengths, vector<int> indexArray){
	int i;

	if(numR < numCavs){
		numR++;
		for(i=0;i<numDivs;++i){
			indexArray[numR - 1] = i;
			generateLandscape(numR,numP,numCavs,numDivs,landscape,reflVals,phaseVals,lengths,indexArray);
		}
	}
	else if(numP < numCavs){
		numP++;
		for(i=0;i<numDivs;++i){
			indexArray[numCavs + numP - 1] = i;
			generateLandscape(numR,numP,numCavs,numDivs,landscape,reflVals,phaseVals,lengths,indexArray);
		}
	}
	else{
		vector<cavity> holder;
		for(i=0;i<numCavs;++i){
			holder.push_back(cavity(lengths[i],reflVals[i][indexArray[i]],phaseVals[i][indexArray[i+numCavs]]));
		}
		landscape.push_back(holder);
	}
	return;
}

// Called from the main, this will generate the holder variables used in the recursive function
void generateLandscapeSetup(const int &numCavs, const int &numDivs, vector< vector<cavity> > &landscape, 
							const vector< vector<double> > &reflVals, const vector< vector<double> > &phaseVals, 
							const vector<int> lengths){
	int numR = 0;
	int numP = 0;
	vector<int> indexArray(numCavs * 2, -1);

	generateLandscape(numR,numP,numCavs,numDivs,landscape,reflVals,phaseVals,lengths,indexArray);

	return;
}

double MSE(const vector<double> &ideal_intensities, const vector<pulse> &calc){
	int numPoints = ideal_intensities.size();
	double error = 0;
	double temp;

	for(int i=0;i<numPoints;++i){
		temp = calc[i].real * calc[i].real + calc[i].imag * calc[i].imag;
		error += (ideal_intensities[i]-temp)*(ideal_intensities[i]-temp);
	}

	error = error / numPoints;

	return error;
}

int main (){
	int numCavs = 4;
	int numDivs = DIV_PER_DIM;
	double phaseMin = 0;
	double phaseMax = 2*PI;
	double reflMin = R_MIN;
	double reflMax = R_MAX;

	vector< vector<cavity> > landscape;
	vector< vector<double> > reflVals(numCavs, vector<double>(numDivs));
	vector< vector<double> > phaseVals(numCavs, vector<double>(numDivs));
	vector<int> lengths{1,1,1,1};

	// Create equal amplitude pulse train
	vector<double> ideal_intensities(NUM_PULSES,1/NUM_PULSES); 

	for(int i = 0; i < numCavs; ++i){
		reflVals[i][0] = 0.5;
		reflVals[i][1] = 0.6;
		phaseVals[i][0] = -0.2;
		phaseVals[i][1] = 0.9;
	}

	generateLandscapeSetup(numCavs,numDivs,landscape,reflVals,phaseVals,lengths);

	int temp = 1 << numCavs;
	
	for(int i = 0; i < temp; ++i){
		cout << "Values of cavity set " << i << endl;
		for(int j = 0; j < numCavs; ++j){
			cout << "Cavity " << j << ":" << endl;
			cout << "Reflectivity: " << landscape[i][j].refl << endl;
			cout << "Phase: " << landscape[i][j].phaseReal << ", " << landscape[i][j].phaseImag << endl;
			cout << "Length: " << landscape[i][j].length << endl << endl;
		}
		cout << endl << endl;
	}
	
	return 0;
}