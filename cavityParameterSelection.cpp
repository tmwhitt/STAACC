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
#define NUM_PULSES 9
#define NUM_CAVS 4
#define POINTS_PER_THREAD 10 		// This number is the number of points to be tested and kept per thread in the random search
#define ERROR_BOUND 0.0000001
#define ERROR_SCALE 0.000000000000000001
#define INITIAL_POINTS 10000000
#define SIMPLEX_POINTS 10000000

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
	pair< vector<cavity>, double> testCavSet(vector<cavity>(numCavs), 1.0);
	pair< vector<cavity>, double> testCavSet2(vector<cavity>(numCavs), 1.0);
	vector< pair< vector<cavity>, double> > simplex(numCavs*2+1, make_pair(vector<cavity>(numCavs), 1.0));
	vector<cavity> meanCavSet(numCavs);

	vector< pair< vector<cavity>, double> > bestCavSet;

	// Create pulse profiles for use
	vector<double> ideal_intensities(NUM_PULSES,1.0/NUM_PULSES); 
	vector<pulse> input(NUM_PULSES,0.0);
	vector<pulse> output;
	input[0].real = 1.0;

	// Set cavity lengths, they dont vary
	for(i=0;i<numCavs;++i){
		cavSet.first[i].length = lengths[i];
		testCavSet.first[i].length = lengths[i];
		testCavSet2.first[i].length = lengths[i];
		meanCavSet[i].length = lengths[i];
	}

	// Initial random scan to find starting point of program
	#pragma omp parallel default(shared) firstprivate(k,i,cavSet,input,output,ideal_intensities)
	{
		seed_seq seed1 {omp_get_thread_num()};
		default_random_engine generator (seed1);
		uniform_real_distribution<double> reflDist(reflMin,reflMax);
		uniform_real_distribution<double> phaseDist(phaseMin,phaseMax);

		pair< vector<cavity>, double> localCavSet = cavSet;
		vector< pair< vector<cavity>, double> > localBestCavSet(POINTS_PER_THREAD, make_pair(vector<cavity>(numCavs), 1.0));
		#pragma omp for
		for(k=0;k<INITIAL_POINTS;++k){
			for(i=0;i<numCavs;++i){
				localCavSet.first[i].setR(reflDist(generator));
				localCavSet.first[i].setPhase(phaseDist(generator));
			}
			localCavSet.second = evaluateError(ideal_intensities,input,output,localCavSet.first);
			// cout << "Output intensities for iteration " << k << endl;
			// for(i=0;i<NUM_PULSES;++i){
			// 	cout << output[i].getIntensity() << endl;
			// }
			// cout << "Error calculated to be " << currError << endl;
			if(localCavSet.second < localBestCavSet.front().second){
				localBestCavSet[0] = localCavSet;
				sort(localBestCavSet.begin(), localBestCavSet.end(), sortError);
				// cout << "Changed minError at " << k << " on thread " << omp_get_thread_num() << endl;
			}
			// cout << endl;
		}

		#pragma omp critical
		{
			// cout << "Minimum error found on thread " << omp_get_thread_num() << " is " << localBestCavSet.second << endl;
			for(i=0;i<POINTS_PER_THREAD;++i){
				bestCavSet.push_back(localBestCavSet[i]);
			}
		}
	}
	sort(bestCavSet.begin(), bestCavSet.end(), sortError);

	// cout << "Minimum error found in random is " << bestCavSet.back().second << endl;

	// cout << "Values used for cavities are:" << endl;
	// for(i=0;i<numCavs;++i){
	// 	cout << "R=" << bestCavSet.back().first[i].getR() << ", phi=" << bestCavSet.back().first[i].getPhase() << endl;
	// }
	// cout << endl;

	cavityPropagation(bestCavSet.back().first,input,output);
	// cout << "Output intensities:" << endl;
	// for(i=0;i<NUM_PULSES;++i){
	// 	cout << output[i].getIntensity() << endl;
	// }
	// cout << endl;

	// cout << "Ideal intensities: " << endl;
	// for(i=0;i<NUM_PULSES;++i){
	// 	cout << ideal_intensities[i] << endl;
	// }

	// Set up initial simplex and sort it according to error
	int numTests = bestCavSet.size();
	for(k=0;k<numTests;++k){
		double simplexSize = simplex.size();
		double startDeviation = 1.05;
	
		simplex[0] = bestCavSet[k];
		for(i=1;i<simplexSize;++i){
			simplex[i] = simplex[0];
			if(i <= NUM_CAVS){
				simplex[i].first[i-1].setR(simplex[i].first[i-1].getR() * startDeviation);
			}
			else{
				simplex[i].first[i-NUM_CAVS-1].setPhase(simplex[i].first[i-NUM_CAVS-1].getPhase() * startDeviation);
			}
			simplex[i].second = evaluateError(ideal_intensities, input, output, simplex[i].first);
		}
		sort(simplex.begin(), simplex.end(), sortError);
	
		double difference = 1;
		int iteration = 0;
		// This loop will execute until the minimum error is within the bound set
		while((simplex.back().second > ERROR_BOUND) && iteration < SIMPLEX_POINTS && difference > ERROR_SCALE){
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
			if(iteration % 10000 == 0){
				cout << difference << endl;
			}
			difference = simplex.front().second - simplex.back().second;
		}
		cout << "Ended on iteration " << iteration << endl;
		cout << "Final error found is " << simplex.back().second << endl << endl;
	
		// cout << "Values used for cavities are:" << endl;
		// for(i=0;i<numCavs;++i){
		// 	cout << "R=" << simplex.back().first[i].getR() << ", phi=" << simplex.back().first[i].getPhase() << endl;
		// }
	
		// cavityPropagation(simplex.back().first,input,output);
		// cout << "Output intensities:" << endl;
		// for(i=0;i<NUM_PULSES;++i){
		// 	cout << output[i].getIntensity() << endl;
		// }
		// cout << endl;
	}
	outfile.close();
	return 0;
}