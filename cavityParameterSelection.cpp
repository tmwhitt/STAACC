#include <vector>
#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <random>
#include <algorithm>
#include "omp.h"
#include "timePropagation.h"

using namespace std;

#define R_MIN 0.25
#define R_MAX 0.75
#define PI 3.141592653589793238463
#define NUM_PULSES 81
#define NUM_CAVS 8
#define ERROR_BOUND 0.00001
#define INITIAL_POINTS 1000000

bool sortError(const pair< vector<cavity>, double> &a,
				const pair< vector<cavity>, double> &b){
	return (a.second > b.second);
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

double evaluateError(	const vector<double> &ideal_intensities, const vector<pulse> &input,
						vector<pulse> &output, vector<cavity> &cavSet){
	cavityPropagation(cavSet,input,output);
	return MSE(ideal_intensities,output);
}

int main (int argc, char *argv[]){
	int i,j,k;
	int numCavs = NUM_CAVS;

	ofstream outfile;
	outfile.open("testOut.txt");
	if(!outfile){cout << "Error opening file" << endl; return -1;}

	double phaseMin = 0;
	double phaseMax = 2*PI;
	double reflMin = R_MIN;
	double reflMax = R_MAX;
	double meanR;
	double meanP;
	vector<int> lengths{1,1,1,1,9,9,9,9};

	pair< vector<cavity>, double> cavSet(vector<cavity>(numCavs), 1.0);
	pair< vector<cavity>, double> bestCavSet(vector<cavity>(numCavs), 1.0);
	pair< vector<cavity>, double> testCavSet(vector<cavity>(numCavs), 1.0);
	pair< vector<cavity>, double> testCavSet2(vector<cavity>(numCavs), 1.0);
	vector< pair< vector<cavity>, double> > simplex(numCavs*2+1, make_pair(vector<cavity>(numCavs), 1.0));
	vector<cavity> meanCavSet(numCavs);

	default_random_engine generator;

	// Create pulse profiles for use
	vector<double> ideal_intensities(NUM_PULSES,1.0/NUM_PULSES); 
	vector<pulse> input(NUM_PULSES,0.0);
	vector<pulse> output;
	input[0].real = 1.0;

	// Set cavity lengths, they dont vary
	for(i=0;i<numCavs;++i){
		cavSet.first[i].length = lengths[i];
		bestCavSet.first[i].length = lengths[i];
		testCavSet.first[i].length = lengths[i];
		testCavSet2.first[i].length = lengths[i];
		meanCavSet[i].length = lengths[i];
	}
	bestCavSet.second = 1.0;

	uniform_real_distribution<double> reflDist(reflMin,reflMax);
	uniform_real_distribution<double> phaseDist(phaseMin,phaseMax);

	// Initial random scan to find starting point of program
	for(k=0;k<INITIAL_POINTS;++k){
		for(i=0;i<numCavs;++i){
			cavSet.first[i].setR(reflDist(generator));
			cavSet.first[i].setPhase(phaseDist(generator));
		}
		cavSet.second = evaluateError(ideal_intensities,input,output,cavSet.first);
		// cout << "Output intensities for iteration " << k << endl;
		// for(i=0;i<NUM_PULSES;++i){
		// 	cout << output[i].getIntensity() << endl;
		// }
		// cout << "Error calculated to be " << currError << endl;
		if(cavSet.second < bestCavSet.second){
			bestCavSet.swap(cavSet);
			cout << "Changed minError at " << k << endl;
		}
		// cout << endl;
	}

	cout << "Minimum error found is " << bestCavSet.second << endl;

	cavityPropagation(bestCavSet.first,input,output);
	cout << "Output intensities:" << endl;
	for(i=0;i<NUM_PULSES;++i){
		cout << output[i].getIntensity() << endl;
	}
	cout << endl;

	cout << "Ideal intensities: " << endl;
	for(i=0;i<NUM_PULSES;++i){
		cout << ideal_intensities[i] << endl;
	}

	// Set up initial simplex and sort it according to error
	double simplexSize = simplex.size();
	simplex[0].swap(bestCavSet);
	for(i=1;i<simplexSize;++i){
		simplex[i] = simplex[0];
		if(i <= NUM_CAVS){
			simplex[i].first[i-1].setR(simplex[i].first[i-1].getR() * 1.05);
		}
		else{
			simplex[i].first[i-NUM_CAVS-1].setPhase(simplex[i].first[i-NUM_CAVS-1].getPhase() * 1.05);
		}
		simplex[i].second = evaluateError(ideal_intensities, input, output, simplex[i].first);
	}
	sort(simplex.begin(), simplex.end(), sortError);

	int iteration = 0;
	// This loop will execute until the minimum error is within the bound set
	while((simplex.back().second > ERROR_BOUND) && iteration < 100000){
		// First, calculate the mean cavity set
		for(i=0;i<numCavs;++i){
			meanR = 0;
			meanP = 0;
			for(j=1;j<simplexSize;++j){
				meanR += simplex[j].first[i].getR();
				meanP += simplex[j].first[i].getPhase();
			}
			meanR = meanR / (simplexSize - 1);
			meanP = meanP / (simplexSize - 1);
			meanCavSet[i].setR(meanR);
			meanCavSet[i].setPhase(meanP);
		}

		// Then, transformations
		// First, we perform reflection, so find the reflected point
		for(i=0;i<numCavs;++i){
			testCavSet.first[i].setR(2 * meanCavSet[i].getR() - simplex[0].first[i].getR());
			testCavSet.first[i].setPhase(2 * meanCavSet[i].getPhase() - simplex[0].first[i].getPhase());
		}
		testCavSet.second = evaluateError(ideal_intensities, input, output, testCavSet.first);

		if(testCavSet.second < simplex.back().second){
			// Condition met for expansion, check it out
			for(i=0;i<numCavs;++i){
				testCavSet2.first[i].setR(2 * testCavSet.first[i].getR() - meanCavSet[i].getR());
				testCavSet2.first[i].setPhase(2 * testCavSet.first[i].getPhase() - meanCavSet[i].getPhase());
			}
			testCavSet2.second = evaluateError(ideal_intensities, input, output, testCavSet2.first);

			if(testCavSet2.second < testCavSet.second){
				simplex[0].swap(testCavSet2);
			}
			else{
				simplex[0].swap(testCavSet);
			}

		}
		else if(testCavSet.second < simplex.front().second){
			// Condition met for reflection only, simply set it
			simplex[0].swap(testCavSet);
		}
		else{
			// Condition met for contraction or shrink, do it
			for(i=0;i<numCavs;++i){
				testCavSet2.first[i].setR((simplex[0].first[i].getR() + meanCavSet[i].getR())/2.0);
				testCavSet2.first[i].setPhase((simplex[0].first[i].getPhase() + meanCavSet[i].getPhase())/2.0);
			}
			testCavSet2.second = evaluateError(ideal_intensities, input, output, testCavSet2.first);

			if(testCavSet2.second < simplex.front().second){
				simplex[0].swap(testCavSet2);
			}
			else{
				// Condition met for shrinking, unfortunately. It's gonna be a doozey
				for(i=0;i<simplexSize-1;++i){
					for(j=0;j<numCavs;++j){
						simplex[i].first[j].setR((simplex[i].first[j].getR() + simplex.back().first[j].getR())/2.0);
						simplex[i].first[j].setPhase((simplex[i].first[j].getPhase() + simplex.back().first[j].getPhase())/2.0);
					}
				}
			}

		}

		// Lastly, re-sort with the transformed node and go again
		sort(simplex.begin(), simplex.end(), sortError);
		++iteration;
	}
	cout << "Ended on iteration " << iteration << endl;
	cout << "Final error found is " << simplex.back().second << endl;

	cavityPropagation(simplex.back().first,input,output);
	cout << "Output intensities:" << endl;
	for(i=0;i<NUM_PULSES;++i){
		cout << output[i].getIntensity() << endl;
	}
	cout << endl;

	outfile.close();
	return 0;
}