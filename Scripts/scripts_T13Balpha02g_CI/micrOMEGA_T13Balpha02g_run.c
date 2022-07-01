// standard C headers
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
// micrOMEGAs headers
#include "../include/micromegas.h"
#include "../include/micromegas_aux.h"
#include "lib/pmodel.h"

// This minimal code writes the outputs of micrOMEGA in a file.
// The variables are stored in the "micrOMEGA_output.out" in the form of a table with the following entries:


int main(int argc, char** argv)
{  

	static const char OUTPUT_FILENAME[] = "micrOMEGA_output.out";


	int err = 0;
	// micrOMEGAs settings/initialization
	ForceUG = 1;
 	VWdecay = 0;
	VZdecay = 0;

	// use old (≤ v2.4.5) micrOMEGAs scalar form factors
	//static const bool OLD_FF = false;
	//if (OLD_FF) calcScalarQuarkFF(0.553, 18.9, 55., 243.5);
 	cleanDecayTable();


	char lspname[10], nlspname[10];

	double cut = 0.01;		// cut-off for channel output								
	int fast = 1;			/* 0 = best accuracy, 1 = "fast option" accuracy ~1% 	     */
 	double Beps = 1.E-5;  		/* Criteqrium for including co-annihilations (1 = no coann.) */


	FILE *output = fopen(OUTPUT_FILENAME, "a");
		
	err = sortOddParticles(lspname);	
	
	//printf("Dark matter candidate is %s (mass = %.2f GeV)\n", CDM1, Mcdm1);

	//printMasses(stdout,1);
	//printHiggs(stdout);
	
	// calculate relic density
	double Xf=-1;
	double relic_density = darkOmega(&Xf, fast, Beps);
	
	//printf("Xf = %.2E, Ωh² = %.2E\n\n", Xf, relic_density);

	// calculate Direct Detection cross sections
	// nucleon mass in GeV
	const double Mnuc = 0.939;
	// ℏc/GeV² in picobarn
	const double INVERSE_GEV_SQUARED_TO_PB = 3.8937966e8;
	const double SCcoeff = 4 / M_PI * INVERSE_GEV_SQUARED_TO_PB * pow(
		Mnuc * Mcdm / (Mnuc + Mcdm), 2.
	);
	double pAsi[2], pAsd[2], nAsi[2], nAsd[2];
	nucleonAmplitudes(CDM1, pAsi, pAsd, nAsi, nAsd);
	// Direct Detection Spin Independent cross section
	double DDSI = SCcoeff * pow(pAsi[0], 2);
	// Direct Detection Spin Dependent cross section
	double DDSD = 3 * SCcoeff * pow(pAsd[0], 2);


	// write data to files
	double Q = findValW("Q");
	// output file
	fprintf(output,
	// parameters
	"%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t"
	// Neutrino masses
	"%.5E\t%.5E\t%.5E\t"
	// Gauge bosons masses
	//"%.5E\t%.5E\t"
	// Higgs parameters
	"%.5E\t"
	//"%.5E\t"
	// BSM masses
	"%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t"
	// Observables
	"%.5E\t%.5E\t%.5E\t%.5E\t\n",
	// parameters
	findValW("lam4"), findValW("lam5"),
	findValW("lam111"), findValW("lam112"), 
	findValW("lam121"), findValW("lam122"),
	findValW("lam611"), findValW("lam612"),
	findValW("lam621"), findValW("lam622"),
	findValW("lam631"), findValW("lam632"),
	// Neutrino masses
	pMass("nu1"), pMass("nu2"), pMass("nu3"),
	// Gauge bosons masses
	//pMass("Z"), pMass("Wp"),
	// Higgs parameters
	pMass("h"),
	//findValW("v"),
	// BSM masses
	slhaVal("MPHI2",Q,2,1,1), slhaVal("MPHI2",Q,2,1,2), slhaVal("MPHI2",Q,2,2,1), slhaVal("MPHI2",Q,2,2,2), slhaVal("T13B",Q,1,12), slhaVal("T13B",Q,1,13), pMass("~psi"),
	pMass("~ch1"), pMass("~ch2"), pMass("~ch3"),
	// Observables
	Mcdm1, DDSI, DDSD, relic_density
	);


	fclose(output);

}
