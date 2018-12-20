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
#define SIMPLEX_POINTS 100000

/**
 * @brief [Custom comparator for cavity sets]
 * @details [Pairs of cavity sets & their respective Mean Squared Error are sorted according to the MSE]
 */
bool sortError(const pair< vector<cavity>, double> &a,
				const pair< vector<cavity>, double> &b){
	return (a.second > b.second);
}

/**
 * @brief [Calculates the Mean Squared Error (MSE)]
 * @details [The mean squared error of the intensities of two pulse train intensities]
 * 
 * @param ideal_intensities [A vector of doubles that represents the ideal pulse train intensities]
 * @param calc [A vector of pulses that was calculated and being evaluated]
 * @return [Returns the Mean Squared Error]
 */
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

/**
 * @brief [Helper function for calculating MSE]
 * @details [Calculates the impulse response through a cavity set, then takes the output and calculates the MSE]
 * 
 * @param ideal_intensities [A vector of doubles that represents the ideal pulse train intensities]
 * @param input [Vector of pulses that will be input to the cavities, the number of output pulses calculated is the size of this vector]
 * @param output [Vector of pulses that will be overwritten with the output train, do not store important information in this vector]
 * @param cavSet [Vector of cavities that describes the stacker system]
 * @return [Returns the Mean Squared Error]
 */
double evaluateError(	const vector<double> &ideal_intensities, const vector<pulse> &input,
						vector<pulse> &output, vector<cavity> &cavSet){
	cavityPropagation(cavSet,input,output);
	return MSE(ideal_intensities,output);
}

int main (int argc, char *argv[]){
	int i,j,k;
	int numCavs = NUM_CAVS;

	// ofstream outfile;
	// outfile.open("testOut.txt");
	// if(!outfile){cout << "Error opening file" << endl; return -1;}

	double phaseMin = 0;
	double phaseMax = 2*PI;
	double reflMin = R_MIN;
	double reflMax = R_MAX;
	double meanR;
	double meanP;
	vector<int> lengths{1,1,1,1,9,9,9,9};

	pair< vector<cavity>, double> cavSet(vector<cavity>(numCavs), 1.0);
	pair< vector<cavity>, double> reflCavSet(vector<cavity>(numCavs), 1.0);
	pair< vector<cavity>, double> expandCavSet(vector<cavity>(numCavs), 1.0);
	pair< vector<cavity>, double> contractCavSet(vector<cavity>(numCavs), 1.0);
	vector< pair< vector<cavity>, double> > simplex(numCavs*2, make_pair(vector<cavity>(numCavs), 1.0));
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
		reflCavSet.first[i].length = lengths[i];
		expandCavSet.first[i].length = lengths[i];
		contractCavSet.first[i].length = lengths[i];
		meanCavSet[i].length = lengths[i];
	}

	cavSet.first[0].setPhase(0);
	reflCavSet.first[0].setPhase(0);
	expandCavSet.first[0].setPhase(0);
	contractCavSet.first[0].setPhase(0);
	meanCavSet[0].setPhase(0);

	auto start_time = chrono::high_resolution_clock::now();

	// Initial random scan to find starting point of program
	#pragma omp parallel default(shared) firstprivate(k,i,cavSet,input,output,ideal_intensities)
	{
		#pragma omp single
			cout << "Number of threads being used: " << omp_get_num_threads() << endl;

		seed_seq seed1 {omp_get_thread_num()};
		default_random_engine generator (seed1);
		uniform_real_distribution<double> reflDist(reflMin,reflMax);
		uniform_real_distribution<double> phaseDist(phaseMin,phaseMax);

		pair< vector<cavity>, double> localCavSet = cavSet;
		vector< pair< vector<cavity>, double> > localBestCavSet(POINTS_PER_THREAD, make_pair(vector<cavity>(numCavs), 1.0));
		#pragma omp for
		for(k=0;k<INITIAL_POINTS;++k){
			localCavSet.first[0].setR(reflDist(generator));
			for(i=1;i<numCavs;++i){
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

	auto end_time = chrono::high_resolution_clock::now();
	auto duration = end_time - start_time;
	cout << "Random point testing took " << chrono::duration_cast<chrono::seconds>(duration).count() << " seconds." << endl;

	// Set up initial simplex and sort it according to error
	int numTests = bestCavSet.size();	// Size will be numThreads*POINTS_PER_THREAD
	start_time = chrono::high_resolution_clock::now();

	#pragma omp parallel for default(shared) firstprivate(i,j,simplex, ideal_intensities, input, output, meanCavSet, reflCavSet, expandCavSet, contractCavSet)
		for(k=0;k<numTests;++k){	
			bool shrinkNeeded = false;
			double simplexSize = simplex.size();
			double startDeviation = 1.05;
	
			int reflect = 0;
			int expand = 0;
			int contract = 0;
			int shrink = 0;
		
			simplex[0] = bestCavSet[k];
			for(i=1;i<simplexSize;++i){
				simplex[i] = simplex[0];
				if(i <= NUM_CAVS){
					simplex[i].first[i-1].setR(simplex[i].first[i-1].getR() * startDeviation);
				}
				else{
					simplex[i].first[i-NUM_CAVS].setPhase(simplex[i].first[i-NUM_CAVS].getPhase() * startDeviation);
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
					reflCavSet.first[i].setR(2 * meanCavSet[i].getR() - simplex[0].first[i].getR());
					reflCavSet.first[i].setPhase(2 * meanCavSet[i].getPhase() - simplex[0].first[i].getPhase());
				}
				reflCavSet.second = evaluateError(ideal_intensities, input, output, reflCavSet.first);
		
				if(reflCavSet.second < simplex.back().second){
					// Condition met for expansion, check it out
					for(i=0;i<numCavs;++i){
						expandCavSet.first[i].setR(2 * reflCavSet.first[i].getR() - meanCavSet[i].getR());
						expandCavSet.first[i].setPhase(2 * reflCavSet.first[i].getPhase() - meanCavSet[i].getPhase());
					}
					expandCavSet.second = evaluateError(ideal_intensities, input, output, expandCavSet.first);
		
					if(expandCavSet.second < reflCavSet.second){
						simplex[0].swap(expandCavSet);
						expand++;
					}
					else{
						simplex[0].swap(reflCavSet);
						reflect++;
					}
		
				}
				else if(reflCavSet.second < simplex[1].second){
					// Condition met for reflection only, simply set it
					simplex[0].swap(reflCavSet);
					reflect++;
				}
				else{
					if(reflCavSet.second < simplex[0].second){
						for(i=0;i<numCavs;++i){
							contractCavSet.first[i].setR((reflCavSet.first[i].getR() + meanCavSet[i].getR())/2.0);
							contractCavSet.first[i].setPhase((reflCavSet.first[i].getPhase() + meanCavSet[i].getPhase())/2.0);
						}
						contractCavSet.second = evaluateError(ideal_intensities, input, output, contractCavSet.first);
	
						if(contractCavSet.second <=	 reflCavSet.second){
							simplex[0].swap(contractCavSet);
							contract++;
						}
						else{
							shrinkNeeded = true;
						}
					}
					else{
						for(i=0;i<numCavs;++i){
							contractCavSet.first[i].setR((simplex[0].first[i].getR() + meanCavSet[i].getR())/2.0);
							contractCavSet.first[i].setPhase((simplex[0].first[i].getPhase() + meanCavSet[i].getPhase())/2.0);
						}
						contractCavSet.second = evaluateError(ideal_intensities, input, output, contractCavSet.first);
	
						if(contractCavSet.second <=	 simplex[0].second){
							simplex[0].swap(contractCavSet);
							contract++;
						}
						else{
							shrinkNeeded = true;
						}
					}
					if(shrinkNeeded){
						// Condition met for shrinking, unfortunately. It's gonna be a doozey
						for(i=0;i<simplexSize-1;++i){
							for(j=0;j<numCavs;++j){
								simplex[i].first[j].setR((simplex[i].first[j].getR() + simplex.back().first[j].getR())/2.0);
								simplex[i].first[j].setPhase((simplex[i].first[j].getPhase() + simplex.back().first[j].getPhase())/2.0);
							}
							simplex[i].second = evaluateError(ideal_intensities, input, output, simplex[i].first);
						}
						shrink++;
						shrinkNeeded = false;
					}
				}
		
				// Lastly, re-sort with the transformed node and go again
				sort(simplex.begin(), simplex.end(), sortError);
	
				++iteration;
				difference = simplex.front().second - simplex.back().second;
			}
			cout << "Ended on iteration " << iteration << endl;
			cout << "Final error found is " << simplex.back().second << endl;
	
			if(difference <= ERROR_SCALE){
				cout << "Terminated due to all points being almost the same" << endl;
			}
			else if(iteration >= SIMPLEX_POINTS){
				cout << "Terminated due to timeout" << endl;
			}
			else{
				cout << "Terminated due to good enough point" << endl;
			}
	
			cout << "Number of operations:" << endl;
			cout << "Reflections: " << reflect << endl;
			cout << "Expansions: " << expand << endl;
			cout << "Contractions: " << contract << endl;
			cout << "Shrinks: " << shrink << endl;
	
			cout << "Values used for cavities are:" << endl;
			for(i=0;i<numCavs;++i){
				cout << "R=" << simplex.back().first[i].getR() << ", phi=" << simplex.back().first[i].getPhase() << endl;
			}
		
			cout << endl;
		}

	end_time = chrono::high_resolution_clock::now();
	duration = end_time - start_time;
	cout << "Simplexing all points took " << chrono::duration_cast<chrono::seconds>(duration).count() << " seconds." << endl;


	// outfile.close();
	return 0;
}