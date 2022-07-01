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
//
// [1.mphiInput, 2.mCPsiInput, 3.mpsipsipInput, 4.lam1, 5.lam3, 6.lam4, 7.lam5, 8.lam61, 9.lam62, 10.lam63, 11.Lam,
// 12.nu1, 13.nu2, 14.nu3, 15.h, 16.~php, 17.~ph0, 18.~psi, 19.~ch1, 20.~ch2, 21.~ch3, 22.DMmass 23.relic_density, 24.Xf]


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
	
	//printMasses(stdout,1);
	
	//printf("Dark matter candidate is %s (mass = %.2f GeV)\n", CDM1, Mcdm1);
	
	double Xf=-1;
	double relic_density = darkOmega(&Xf, fast, Beps);
	
	//printf("Xf = %.2E, Ωh² = %.2E\n\n", Xf, relic_density);


	// write data to files
	double Q = findValW("Q");
	// output file
	fprintf(output,
	// parameters
	"%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t"
	"%.5E\t"
	// Neutrino masses
	"%.5E\t%.5E\t%.5E\t"
	// Gauge bosons masses
	//"%.5E\t%.5E\t"
	// Higgs parameters
	"%.5E\t"
	//"%.5E\t"
	// BSM masses
	"%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t"
	// Observables
	"%.5E\t%.5E\t%.5E\t\n",
	// parameters
	slhaVal("MINPAR", Q, 1, 11), slhaVal("MINPAR", Q, 1, 12),
	slhaVal("MINPAR", Q, 1, 13), findValW("lam1"),
	findValW("lam3"), findValW("lam4"), findValW("lam5"),
	findValW("lam61"), findValW("lam62"), findValW("lam63"),
	findValW("Lam"), 
	// Neutrino masses
	pMass("nu1"), pMass("nu2"), pMass("nu3"),
	// Gauge bosons masses
	//pMass("Z"), pMass("Wp"),
	// Higgs parameters
	pMass("h"),
	//findValW("v"),
	// BSM masses
	pMass("~php"), pMass("~ph0"), pMass("~psi"),
	pMass("~ch1"), pMass("~ch2"), pMass("~ch3"),
	// Observables
	Mcdm1, relic_density, Xf
	);


	fclose(output);

}
