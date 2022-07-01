// standard C headers
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
// GLib headers
#include <gmodule.h>
// micrOMEGAs headers
#include "../include/micromegas.h"
#include "../include/micromegas_aux.h"
#include "lib/pmodel.h"

#ifndef M_PI
static const double M_PI = (3.14159265358979323846);
#endif
static const char MODEL_NAME[] = "T13Balpha02g";
static const char OUTPUT_FILENAME[] = "micromegas_data.out";
static const char CHANNELS_FILENAME[] = "micromegas_channels.out";
static const char SPHENO_FILENAME[] = "spheno_data.out";
static const char RUN_INFO_FILE[] = "parameter_scan_info.txt";
static char *TMP_DIR;
static GString *RUN_INFO;
static GString *SPHENO_OUTPUT;
static const unsigned int SLHAREAD_PATH_MAX = 60;
// use old (≤ v2.4.5) micrOMEGAs scalar form factors
static const bool OLD_FF = false;
// 0 = best accuracy, 1 = “fast option” (accuracy ≈1%)
static const int FAST = 1;
// criterion for including co-annihilations (1 = no coann.)
static const double BEPS = 1e-5;
// cut-off for channel output
static const double CUT = 0.01;

static void init_global(int argc, char **argv) {
	// turn on line buffering for stdout
	if (setvbuf(stdout, NULL, _IOLBF, 1024)) {
		printf("Could not set stdout to line buffering!\n");
		exit(10);
	}
	if (argc != 2) {
		printf("Correct usage: %s <tmpdir>\n", argv[0]);
		printf("Example: %s /tmp/micromegas-run1/\n", argv[0]);
		exit(1);
	}
	TMP_DIR = argv[1];
	// read RUN_INFO_FILE
	RUN_INFO = g_string_new(NULL);
	GString *run_info_path = g_string_new(NULL);
	g_string_printf(run_info_path, "%s/%s", TMP_DIR, RUN_INFO_FILE);
	FILE *run_info_file = fopen(run_info_path->str, "r");
	g_string_free(run_info_path, true);
	const size_t buf_len = 512;
	char buf[buf_len];
	size_t read_len;
	while ((read_len = fread(buf, sizeof(char), buf_len, run_info_file))) {
		g_string_append_len(RUN_INFO, buf, read_len);
	}
	fclose(run_info_file);
	// construct path to SPheno spectrum file
	SPHENO_OUTPUT = g_string_new(NULL);
	g_string_printf(SPHENO_OUTPUT, "%s/SPheno.spc.%s", TMP_DIR, MODEL_NAME);
	if (SPHENO_OUTPUT->len > SLHAREAD_PATH_MAX) {
		printf("Warning! Path to SPheno spectrum file is\n"
			"    %s\n"
			"    but slhaRead cannot handle paths longer than %d characters!\n",
			SPHENO_OUTPUT->str, SLHAREAD_PATH_MAX
		);
	}
}

int main(int argc, char **argv) {
	// initialization
	init_global(argc, argv);
	int err = 0;
	// micrOMEGAs settings/initialization
	ForceUG = 1;
	VWdecay = 0;
	VZdecay = 0;
	if (OLD_FF) calcScalarQuarkFF(0.553, 18.9, 55., 243.5);
	cleanDecayTable();

	// output file
	// XXX allow for variable output filenames
	FILE *output = fopen(OUTPUT_FILENAME, "w");
	FILE *channels = fopen(CHANNELS_FILENAME, "w");
	FILE *spheno = fopen(SPHENO_FILENAME, "w");
	char input_parameters[] =
		"mϕ₁₁\tmϕ₁₂\tmϕ₂₁\tmϕ₂₂\tmΨ\tmψψ'\tλ₁₁₁\tλ₁₁₂\tλ₁₂₁\tλ₁₂₂\tλ₄\tλ₅\t"
		"λ₆₁₁\tλ₆₂₁\tλ₆₃₁\tλ₆₁₂\tλ₆₂₂\tλ₆₃₂\tλ\tv";
	// output file headers/comments
	fprintf(output,
		"# %s\t%s\n",
		// parameters
		input_parameters,
		// masses
		"md1\tmd2\tmd3\tmu1\tmu2\tmu3\tme1\tme2\tme3\tmν1\tmν2\tmν3\t"
		"mZ\tmW\tmH\tmη₁+\tmη₂+\tmη₁0\tmη₂0\tmψ\tmχ1\tmχ2\tmχ3\t"
		// decay widths
		"ΓZ\tΓW\tΓH\t"
		// observables
		"Ωh²\tXf\tσp(SI)\tσp(SD)"
	);
	fprintf(output, "%s", RUN_INFO->str);
	fprintf(output, "# %s\n", WORK);
	fprintf(output, "# ForceUG = %d; VWdecay = %d; VZdecay = %d; "
		"FAST = %d; BEPS = %E;\n",
		ForceUG, VWdecay, VZdecay, FAST, BEPS
	);
	fprintf(output, "# calcScalarQuarkFF: %d\n", OLD_FF);
	fprintf(channels, "# %s\n", input_parameters);
	fprintf(spheno, "# %s\t%s\n",
		input_parameters,
		"BR(μ→eγ)\tBR(τ→eγ)\tBR(τ→μγ)\tBR(μ→3e)\tBR(τ→3e)\tBR(τ→3μ)"
	);
	// process all parameter points
	int c;
	// keep running while the calling process is still writing to stdin;
	// a character written to stdin signals that the SPheno.spc file is ready
	// to be read again
	while ((c = getchar()) != EOF) {
		// read SPheno output spectrum file
		// NOTE: slhaRead cannot handle paths longer than 60 characters!
		if ((err = slhaRead(SPHENO_OUTPUT->str, 0))) {
			printf("Could not read SPheno spectrum file %s: %d\n",
				SPHENO_OUTPUT->str, err
			);
			break;
		}
		char lopname[10];
		if ((err = sortOddParticles(lopname))) {
			// micrOMEGAs already prints a warning when an error occurs during
			// sortOddParticles
			break;
		}
		printf("Dark matter candidate is %s (mass = %.2f GeV)\n", CDM1, Mcdm1);
		// print masses to stdout
		printMasses(stdout, true);
		printHiggs(stdout);
		// calculate relic density
		double Xf = -1;
		double relic_density = darkOmega(&Xf, FAST, BEPS);
		// calculate cross sections
		// nucleon mass in GeV
		const double Mnuc = 0.939;
		// ℏc/GeV² in picobarn
		const double INVERSE_GEV_SQUARED_TO_PB = 3.8937966e8;
		const double SCcoeff = 4 / M_PI * INVERSE_GEV_SQUARED_TO_PB * pow(
			Mnuc * Mcdm / (Mnuc + Mcdm), 2.
		);
		double pAsi[2], pAsd[2], nAsi[2], nAsd[2];
		nucleonAmplitudes(CDM1, pAsi, pAsd, nAsi, nAsd);
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
	}
	fclose(output);
	fclose(channels);
	fclose(spheno);
	return 0;
}

