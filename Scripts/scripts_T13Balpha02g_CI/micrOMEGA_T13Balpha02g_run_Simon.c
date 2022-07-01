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

/* 0 = best accuracy, 1 = "fast option" accuracy ~1% 	     */
static const int FAST = 1;
// criterion for including co-annihilations (1 = no coann.)
static const double BEPS = 1e-5;
// cut-off for channel output
static const double CUT = 0.01;



int main(int argc, char** argv)
{  

	static const char OUTPUT_FILENAME[] = "micromegas_data.out";
    static const char OUTPUT_CHANNELS[] = "micromegas_channels.out";
    static const char OUTPUT_SPHENO[] = "spheno_data.out";


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

// 	double cut = 0.01;		// cut-off for channel output								
// 	int fast = 1;			/* 0 = best accuracy, 1 = "fast option" accuracy ~1% 	     */
//  double Beps = 1.E-5;  		/* Criteqrium for including co-annihilations (1 = no coann.) */


	FILE *output = fopen(OUTPUT_FILENAME, "a");
    FILE *channels = fopen(OUTPUT_CHANNELS, "a");
    FILE *spheno = fopen(OUTPUT_SPHENO, "a");
		
	err = sortOddParticles(lspname);	
	
	//printf("Dark matter candidate is %s (mass = %.2f GeV)\n", CDM1, Mcdm1);

	//printMasses(stdout,1);
	//printHiggs(stdout);
	
	// calculate relic density
	double Xf=-1;
	double relic_density = darkOmega(&Xf, FAST, BEPS);
	
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

    
    //headers that shall be written only in the beginning
//     char input_parameters[] =
// 		"mϕ₁₁\tmϕ₁₂\tmϕ₂₁\tmϕ₂₂\tmΨ\tmψψ'\tλ₁₁₁\tλ₁₁₂\tλ₁₂₁\tλ₁₂₂\tλ₄\tλ₅\t"
// 		"λ₆₁₁\tλ₆₂₁\tλ₆₃₁\tλ₆₁₂\tλ₆₂₂\tλ₆₃₂\tλ\tv";
// 	// output file headers/comments
// 	fprintf(output,
// 		"# %s\t%s\n",
// 		// parameters
// 		input_parameters,
// 		// masses
// 		"md1\tmd2\tmd3\tmu1\tmu2\tmu3\tme1\tme2\tme3\tmν1\tmν2\tmν3\t"
// 		"mZ\tmW\tmH\tmη₁+\tmη₂+\tmη₁0\tmη₂0\tmψ\tmχ1\tmχ2\tmχ3\t"
// 		// decay widths
// 		"ΓZ\tΓW\tΓH\t"
// 		// observables
// 		"Ωh²\tXf\tσp(SI)\tσp(SD)"
// 	);
// 	fprintf(output, "%s", RUN_INFO->str);
// 	fprintf(output, "# %s\n", WORK);
// 	fprintf(output, "# ForceUG = %d; VWdecay = %d; VZdecay = %d; "
// 		"FAST = %d; BEPS = %E;\n",
// 		ForceUG, VWdecay, VZdecay, FAST, BEPS
// 	);
// 	fprintf(output, "# calcScalarQuarkFF: %d\n", OLD_FF);
// 	fprintf(channels, "# %s\n", input_parameters);
// 	fprintf(spheno, "# %s\t%s\n",
// 		input_parameters,
// 		"BR(μ→eγ)\tBR(τ→eγ)\tBR(τ→μγ)\tBR(μ→3e)\tBR(τ→3e)\tBR(τ→3μ)"
// 	);

    
    // print information to stdout
    printf("Xf = %.2E, Ωh² = %.2E\n\n", Xf, relic_density);
    printChannels(Xf, CUT, BEPS, true, stdout);
    if (relic_density > 0) {
        // write data to files
        double Q = findValW("Q");
        // output file
        fprintf(output,
            // parameters
            "%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t"
            "%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t"
            // SM masses
            "%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t"
            "%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t"
            // BSM masses
            "%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t"
            // decay widths
            "%.5E\t%.5E\t%.5E\t"
            // observables
            "%.5E\t%.5E\t%.5E\t%.5E\n",
            // parameters
            sqrt(slhaVal("MPHI2", Q, 2, 1, 1)),
            sqrt(slhaVal("MPHI2", Q, 2, 1, 2)),
            sqrt(slhaVal("MPHI2", Q, 2, 2, 1)),
            sqrt(slhaVal("MPHI2", Q, 2, 2, 2)),
            slhaVal("MINPAR", Q, 1, 12), slhaVal("MINPAR", Q, 1, 13),
            findValW("lam111"), findValW("lam112"), findValW("lam121"),
            findValW("lam122"), findValW("lam4"), findValW("lam5"),
            findValW("lam611"), findValW("lam621"), findValW("lam631"),
            findValW("lam612"), findValW("lam622"), findValW("lam632"),
            findValW("Lam"), findValW("v"),
            // masses
            pMass("d1"), pMass("d2"), pMass("d3"),
            pMass("u1"), pMass("u2"), pMass("u3"),
            pMass("e1"), pMass("e2"), pMass("e3"),
            pMass("nu1"), pMass("nu2"), pMass("nu3"),
            pMass("Z"), pMass("Wp"), pMass("h"),
            pMass("~etp1"), pMass("~etp2"), pMass("~et01"), pMass("~et02"),
            pMass("~psi"), pMass("~ch1"), pMass("~ch2"), pMass("~ch3"),
            // decay widths
            pWidth("Z", NULL), pWidth("Wp", NULL), pWidth("h", NULL),
            // observables
            relic_density, Xf,
            SCcoeff * pow(pAsi[0], 2), 3 * SCcoeff * pow(pAsd[0], 2)
        );
        // channels file
        fprintf(channels,
            // parameters
            "%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t"
            "%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\n",
            sqrt(slhaVal("MPHI2", Q, 2, 1, 1)),
            sqrt(slhaVal("MPHI2", Q, 2, 1, 2)),
            sqrt(slhaVal("MPHI2", Q, 2, 2, 1)),
            sqrt(slhaVal("MPHI2", Q, 2, 2, 2)),
            slhaVal("MINPAR", Q, 1, 12), slhaVal("MINPAR", Q, 1, 13),
            findValW("lam111"), findValW("lam112"), findValW("lam121"),
            findValW("lam122"), findValW("lam4"), findValW("lam5"),
            findValW("lam611"), findValW("lam621"), findValW("lam631"),
            findValW("lam612"), findValW("lam622"), findValW("lam632"),
            findValW("Lam"), findValW("v")
        );
        printChannels(Xf, CUT, BEPS, true, channels);
        fprintf(channels, "\n");
        // SPheno observables file
        fprintf(spheno,
            // parameters
            "%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t"
            "%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t"
            "%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\n",
            sqrt(slhaVal("MPHI2", Q, 2, 1, 1)),
            sqrt(slhaVal("MPHI2", Q, 2, 1, 2)),
            sqrt(slhaVal("MPHI2", Q, 2, 2, 1)),
            sqrt(slhaVal("MPHI2", Q, 2, 2, 2)),
            slhaVal("MINPAR", Q, 1, 12), slhaVal("MINPAR", Q, 1, 13),
            findValW("lam111"), findValW("lam112"), findValW("lam121"),
            findValW("lam122"), findValW("lam4"), findValW("lam5"),
            findValW("lam611"), findValW("lam621"), findValW("lam631"),
            findValW("lam612"), findValW("lam622"), findValW("lam632"),
            findValW("Lam"), findValW("v"),
            // observables
            // BR(μ → eγ)
            slhaVal("FlavorKitLFV", Q, 1, 701),
            // BR(τ → eγ)
            slhaVal("FlavorKitLFV", Q, 1, 702),
            // BR(τ → μγ)
            slhaVal("FlavorKitLFV", Q, 1, 703),
            // BR(μ → 3e)
            slhaVal("FlavorKitLFV", Q, 1, 901),
            // BR(τ → 3e)
            slhaVal("FlavorKitLFV", Q, 1, 902),
            // BR(τ → 3μ)
            slhaVal("FlavorKitLFV", Q, 1, 903)
        );
    }
    printf("=======================================================\n");
    fclose(output);
    fclose(channels);
    fclose(spheno);
    return 0;

}
