#include <vector>
#include <iostream>
#include <cmath>
#include <chrono>

#include "timePropagation.h"

#define INPUT_SIZE 200

using namespace std;

int main(){
	vector<pulse> input(INPUT_SIZE);
	vector<pulse> output(INPUT_SIZE);

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

	auto start_time = chrono::high_resolution_clock::now();
	cavityPropagation(cavSet,input,output);
	auto end_time = chrono::high_resolution_clock::now();

	auto duration = end_time - start_time;
	cout << "Code took " << chrono::duration_cast<chrono::microseconds>(duration).count() << " microseconds." << endl;

	cout << "Output size = " << output.size() << endl;

	cout << "Pulse train output (real, imag):" << endl;
	for(int i=0;i<output.size();++i){
		cout << output[i].real << "," << output[i].imag << endl;
	}

	return 0;
}