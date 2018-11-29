#include <vector>
#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <random>
#include "omp.h"
#include "timePropagation.h"

using namespace std;

#define R_MIN 0.25
#define R_MAX 0.75
#define PI 3.141592653589793238463
#define DIV_PER_DIM 5
#define NUM_PULSES 9
#define ITERATIONS 20

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

double MSE(const vector<double> &ideal_intensities, vector<pulse> &calc){
	int numPoints = ideal_intensities.size();
	double error = 0;
	double temp;

	for(int i=0;i<numPoints;++i){
		temp = calc[i].getIntensity();
		error += (ideal_intensities[i]-temp)*(ideal_intensities[i]-temp);
	}

	error = error / numPoints;

	return error;
}

int main (int argc, char *argv[]){
	int i,j,k,best;
	int numCavs = 4;
	int numDivs = DIV_PER_DIM;
	int numIterations = ITERATIONS;

	ofstream outfile;
	outfile.open("testOut.txt");
	if(!outfile){cout << "Error opening file" << endl;return -1;}

	vector<double> phaseMin (numCavs,0);
	vector<double> phaseMax (numCavs,2*PI);
	vector<double> reflMin  (numCavs,R_MIN);
	vector<double> reflMax  (numCavs,R_MAX);
	double phaseStep, reflStep;
	double minError = 1000000;

	vector< vector<cavity> > landscape;
	vector< vector<double> > reflVals(numCavs, vector<double>(numDivs));
	vector< vector<double> > phaseVals(numCavs, vector<double>(numDivs));
	vector<int> lengths{1,1,1,1};

	default_random_engine generator;

	// Create equal amplitude pulse train
	vector<double> ideal_intensities(NUM_PULSES,1.0/NUM_PULSES); 
	vector<pulse> input(NUM_PULSES,0.0);
	vector<pulse> output;
	input[0].real = 1.0;

	for(k=0;k<numIterations;++k){
		// Create points to test for each cavity
		phaseStep = (phaseMax[0] - phaseMin[0]) / numDivs;
		uniform_real_distribution<double> phaseDist(0.0,phaseStep);
		reflStep = (reflMax[0] - reflMin[0]) / numDivs;
		uniform_real_distribution<double> reflDist(0.0,reflStep);
		for(i = 0; i < numCavs; ++i){
			for(j=0;j<numDivs;++j){
				reflVals[i][j] = reflMin[i] + j * reflStep + reflDist(generator);
				phaseVals[i][j] = phaseMin[i] + j * phaseStep + phaseDist(generator);
			}
		}
	
		// Generate the landscape
		generateLandscapeSetup(numCavs,numDivs,landscape,reflVals,phaseVals,lengths);
		int numPoints = landscape.size();

		// PARALLELIZE AND CALL DIFFERENT ONES HERE

		// MPI Attempt
		/*
		int rank,numProc,rowStart,rowEnd;
		MPI_Init(&argc,&argv);
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		MPI_Comm_size(MPI_COMM_WORLD,&numProc);

		rowStart = rank * (numPoints / numProc);
		if(numProc < numPoints % numProc){
			rowStart += numProc;
			rowEnd = rowStart + numPoints / numProc + 1;
		}
		else{
			rowStart += numPoints % numProc;
			rowEnd = rowStart + numPoints / numProc;
		}

		cout << "Rank=" << rank << ", start=" << rowStart << ", end=" << rowEnd << endl;
		*/

		#pragma omp parallel default(shared) firstprivate(input,output,ideal_intensities)
		{
			double minErrorLocal = 1000000;
			int bestLocal;
			double temp;
			#pragma omp for private(i)
			for(i=0;i<numPoints;++i){
				cavityPropagation(landscape[i],input,output);
				temp = MSE(ideal_intensities,output);
				if(temp < minErrorLocal){
					bestLocal = i;
					minErrorLocal = temp;
				}
			}
			#pragma omp critical
			{
				if(minErrorLocal < minError){
					minError = minErrorLocal;
					best = bestLocal;
				}
			}
		}

		// double temp;

		// for(i=0;i<numPoints;++i){
		// 	cavityPropagation(landscape[i],input,output);
		// 	temp = MSE(ideal_intensities,output);
		// 	if(temp < minError){
		// 		best = i;
		// 		minError = temp;
		// 	}
		// }

		cout << "Used " << minError << " at " << best << endl;
		cout << "Cavity 1: Reflectivity=" << landscape[best][0].getR() << endl;

		// WILL NEED TO UPDATE MIN, MAX, AND STEP
		double temp;
		for(i=0;i<numCavs;++i){
			temp = (reflMax[i] - reflMin[i])/4;
			reflMax[i] = landscape[best][i].getR() + temp;
			reflMin[i] = landscape[best][i].getR() - temp;
			temp = (phaseMax[i] - phaseMin[i])/4;
			phaseMax[i] = landscape[best][i].getR() + temp;
			phaseMin[i] = landscape[best][i].getR() - temp;
		}

		cavityPropagation(landscape[best],input,output);

		landscape.clear();
		best = -1;
		minError = 1000000;

		cout << "Cavity 1 has maxRefl=" << reflMax[0] << " and minRefl=" << reflMin[0] << endl << endl;
		outfile << "Input, Ideal, Output" << endl;
		for(i=0;i<NUM_PULSES;++i){
			outfile << input[i].getIntensity() << ", " << ideal_intensities[i] << ", " << output[i].getIntensity() << endl;
		}
	}
	
	outfile.close();
	return 0;
}