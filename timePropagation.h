#include <vector>
#include <cmath>

using namespace std;

/**
 * \brief Time domain cavity class
 * \details A struct that fully describes a cavity in the time domain
 * 
 */
struct cavity{

	int length;				/**< \brief Length of cavity in pulse separations.  */
	double refl;			/**< \brief Field reflectivity of cavity beamsplitter.  */
	double trans;			/**< \brief Field transmittivity of cavity beamsplitter.  */
	double phaseReal;		/**< \brief Real portion of cavity phase.  */
	double phaseImag;		/**< \brief Imaginary portion of cavity phase.  */

	cavity() : length(0), refl(1.0), phaseReal(0), phaseImag(0) {}

	cavity(int l, double r, double p) : length(l){
		refl = sqrt(r);
		trans = sqrt(1-r);
		phaseReal = cos(p);
		phaseImag = sin(p);
	}

	void setPhase(double p){
		phaseReal = cos(p);
		phaseImag = sin(p);
	}

	double getPhase(){
		return atan2(phaseImag,phaseReal);
	}

	void setR(double r){
		refl = sqrt(r);
		trans = sqrt(1-r);
	}

	double getR(){
		return refl*refl;
	}
};

/**
 * \brief Optical pulse struct
 * \details A struct that describes an optical pulse by storing its electric field values
 * 
 */
struct pulse{

	double real;		/**< \brief Real portion of pulse electric field.  */
	double imag;		/**< \brief Imaginary portion of pulse electric field.  */

	pulse(double r, double i) : real(r), imag(i) {}
	pulse(double r) : real(r) {}
	pulse() : real(0), imag(0) {}

	double getIntensity(){
		return real*real+imag*imag;
	}

	void setPhase(double p){
		double temp = sqrt(real*real+imag*imag);
		real = cos(p)*temp;
		imag = sin(p)*temp;
	}

	double getPhase(){
		return atan2(imag,real);
	}
};

/**
 * @brief [Propagation function in time domain]
 * @details [Propagates the pulses through a cavity set by calculating the interference of electric fields at each beamsplitter]
 * 
 * @param cavSet [Vector of cavities that describes the stacker system]
 * @param input [Vector of pulses that will be input to the cavities, the number of output pulses calculated is the size of this vector]
 * @param output [Vector of pulses that will be overwritten with the output train, do not store important information in this vector]
 */
void cavityPropagation(	const vector<cavity> &cavSet, const vector<pulse> &input, 
						vector<pulse> &output){

	int i,j;
	int numCavs = cavSet.size();
	int numPulses = input.size();
	double temp;

	vector< vector<pulse> > Ecav(numCavs, vector<pulse>(numPulses));
	vector< vector<pulse> > Eout(numCavs, vector<pulse>(numPulses));

	// First pulse hits all the cavities, and loads them up

	Ecav[0][0].real = -input[0].imag * cavSet[0].trans;
	Ecav[0][0].imag =  input[0].real * cavSet[0].trans;
	Eout[0][0].real =  input[0].real * cavSet[0].refl;
	Eout[0][0].imag =  input[0].imag * cavSet[0].refl;

	temp = Ecav[0][0].real*cavSet[0].phaseReal - Ecav[0][0].imag*cavSet[0].phaseImag;
	Ecav[0][0].imag = Ecav[0][0].real*cavSet[0].phaseImag + Ecav[0][0].imag*cavSet[0].phaseReal;
	Ecav[0][0].real = temp;

	for(j=1;j<numCavs;++j){	
		Ecav[j][0].real = -Eout[j-1][0].imag * cavSet[j].trans;
		Ecav[j][0].imag =  Eout[j-1][0].real * cavSet[j].trans;
		Eout[j][0].real =  Eout[j-1][0].real * cavSet[j].refl;
		Eout[j][0].imag =  Eout[j-1][0].imag * cavSet[j].refl;

		temp = Ecav[j][0].real*cavSet[j].phaseReal - Ecav[j][0].imag*cavSet[j].phaseImag;
		Ecav[j][0].imag = Ecav[j][0].real*cavSet[j].phaseImag + Ecav[j][0].imag*cavSet[j].phaseReal;
		Ecav[j][0].real = temp;
	}

	for(i=1;i<numPulses;++i){
		if(i>=cavSet[0].length){
			Ecav[0][i].real = Ecav[0][i-cavSet[0].length].real * cavSet[0].refl - input[i].imag * cavSet[0].trans;
			Ecav[0][i].imag = Ecav[0][i-cavSet[0].length].imag * cavSet[0].refl + input[i].real * cavSet[0].trans;
			Eout[0][i].real = input[i].real * cavSet[0].refl - Ecav[0][i-cavSet[0].length].imag * cavSet[0].trans;
			Eout[0][i].imag = input[i].imag * cavSet[0].refl + Ecav[0][i-cavSet[0].length].real * cavSet[0].trans;
		}
		else{
			Ecav[0][i].real = -input[i].imag * cavSet[0].trans;
			Ecav[0][i].imag =  input[i].real * cavSet[0].trans;
			Eout[0][i].real =  input[i].real * cavSet[0].refl;
			Eout[0][i].imag =  input[i].imag * cavSet[0].refl;
		}
		temp = Ecav[0][i].real*cavSet[0].phaseReal - Ecav[0][i].imag*cavSet[0].phaseImag;
		Ecav[0][i].imag = Ecav[0][i].real*cavSet[0].phaseImag + Ecav[0][i].imag*cavSet[0].phaseReal;
		Ecav[0][i].real = temp;
		for(j=1;j<numCavs;++j){
			if(i>=cavSet[j].length){
				Ecav[j][i].real = Ecav[j][i-cavSet[j].length].real * cavSet[j].refl - Eout[j-1][i].imag * cavSet[j].trans;
				Ecav[j][i].imag = Ecav[j][i-cavSet[j].length].imag * cavSet[j].refl + Eout[j-1][i].real * cavSet[j].trans;
				Eout[j][i].real = Eout[j-1][i].real * cavSet[j].refl - Ecav[j][i-cavSet[j].length].imag * cavSet[j].trans;
				Eout[j][i].imag = Eout[j-1][i].imag * cavSet[j].refl + Ecav[j][i-cavSet[j].length].real * cavSet[j].trans;
			}
			else{
				Ecav[j][i].real = -Eout[j-1][i].imag * cavSet[j].trans;
				Ecav[j][i].imag =  Eout[j-1][i].real * cavSet[j].trans;
				Eout[j][i].real =  Eout[j-1][i].real * cavSet[j].refl;
				Eout[j][i].imag =  Eout[j-1][i].imag * cavSet[j].refl;
			}
			temp = Ecav[j][i].real*cavSet[j].phaseReal - Ecav[j][i].imag*cavSet[j].phaseImag;
			Ecav[j][i].imag = Ecav[j][i].real*cavSet[j].phaseImag + Ecav[j][i].imag*cavSet[j].phaseReal;
			Ecav[j][i].real = temp;
		}
	}
	output.swap(Eout[numCavs-1]);
	return;
}