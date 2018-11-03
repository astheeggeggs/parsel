// Request decent POSIX version.
#define _XOPEN_SOURCE 700

#include <stdio.h>  // Needed to perform IO operations
#include <stdlib.h>
#include <stdbool.h>
#include <complex.h>
#include <string.h>
#include <strings.h> // strcasecmp.
#include <math.h>
#include <time.h>
#include <limits.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>

#include "constants.h"
#include "amino_acids.h"
#include "util.h"
#include "mcmc_moves.h"
#include "load_data.h"
#include "load_from_json.h"
#include "selection_window.h"

#define DEFAULT_CHAIN_LEN 50000
#define DEFAULT_SAMPLE_LEN 500

// Files opened by hello_parse_cmdline()
FILE *log_file = NULL;
FILE *summary_file = NULL;

int chain_length = DEFAULT_CHAIN_LEN, sample_rate = DEFAULT_SAMPLE_LEN;
int runtime_limit_hrs = 0;

bool runtime_limit_hrs_set = false;
bool number_of_states_set = false;
bool window_mcmc = false;
bool mut_mcmc = false;
bool sel_mcmc = false;
bool ind_mcmc = false;

const char *tree_path = NULL, *sequence_path = NULL, *out_dir = NULL;

static void print_usage()
{
  fprintf(stderr,
          "usage: dmodel [options] <tree_file.json> <sequence_file.fasta> <outdir>\n"
          "  options:\n"
          "    --chain_length               -n <N> Set chain length [default: %i].\n"
          "    --sample_every               -s <S> Sample every <S> steps [default: %i].\n"
          "    --runtime_hours              -t <T> Run for T hours.\n"
          "    --windows                    -w Use reversible-jump MCMC window moves [default: False].\n"
          "    --omega_exponential          -e Prior for the selection coefficients is exponential [default: False].\n"
          "    --mut_only                   -m Only perform mutation moves, for debugging [defaut: False].\n"
          "    --sel_only                   -o Only perform selection moves, for debugging [default: False].\n"
          "    --ind_only                   -i Only perform indel moves, for debugging [default: False].\n"
          "    --sample-prior               -p Sample from the prior, for debugging [default: False].\n",
          DEFAULT_CHAIN_LEN, DEFAULT_SAMPLE_LEN);
  exit(EXIT_FAILURE);
}

// Also opens output files
static void hello_parse_cmdline(int argc, char **argv)
{
  static struct option longopts[] =
  {
    {"chain_length",               required_argument, NULL, 'n'},
    {"sample_every",               required_argument, NULL, 's'},
    {"runtime_hours",              required_argument, NULL, 't'},
    {"windows",                    no_argument,       NULL, 'w'},
    {"omega_exponential",          no_argument,       NULL, 'e'},
    {"mut_only",                   no_argument,       NULL, 'm'},
    {"sel_only",                   no_argument,       NULL, 'o'},
    {"ind_only",                   no_argument,       NULL, 'i'},
    {"sample_prior",               no_argument,       NULL, 'p'},
    {NULL, 0, NULL, 0}
  };

  char shortopts[] = "n:s:t:wemoip";

  int optional;
	while ((optional = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1) {
		switch(optional) {
			case 'n': 
				chain_length = atoi(optarg);
				number_of_states_set = true;
			break;
			case 's': sample_rate = atoi(optarg); break;
			case 't':
				runtime_limit_hrs = atoi(optarg);
				runtime_limit_hrs_set = true;
			break;
			case 'w':
				window_mcmc = true;
			break;
			case 'e': omega_prior = EXPONENTIAL_PRIOR; break;
			case 'm': mut_mcmc = true; break;
			case 'o': sel_mcmc = true; break;
			case 'i': ind_mcmc = true; break;
			case 'p': sample_prior = true; break;
			default: die("Unknown option: %c", optional);
		}
	}
	printf("argc - optind %i\n", argc - optind);
	printf("optind %i\n", optind);
	if (argc - optind != 3) print_usage();
	if (runtime_limit_hrs_set == true && number_of_states_set == true)
		die("Set both a time limit and number of states");
	if (window_mcmc == true)
		printf("Using reversible-jump MCMC moves\n");

	tree_path = argv[optind];
	sequence_path = argv[optind+1];
	out_dir = argv[optind+2];

	if (!mkpath(out_dir, 0755)) die("Cannot create output dir: %s", out_dir);

	time_t rawtime;
	struct tm * timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);
	char date[100];
	sprintf(date, "%d_%d_%d_%d_%d_%d", timeinfo->tm_year+1900, 
		timeinfo->tm_mon+1, timeinfo->tm_mday, timeinfo->tm_hour,
		timeinfo->tm_min, timeinfo->tm_sec);

	char log_path[PATH_MAX+1], summary_path[PATH_MAX+1];

	sprintf(log_path, "%s/%s_log.txt", out_dir, date);
	sprintf(summary_path, "%s/%s_summary.txt", out_dir, date);
  
	while (futil_file_exists(log_path))
	{
		printf("This filename already exists, wait for a second to change the folder name.\n");
		sleep(1);
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		sprintf(date, "%d_%d_%d_%d_%d_%d\n", timeinfo->tm_year+1900, 
			timeinfo->tm_mon+1, timeinfo->tm_mday, timeinfo->tm_hour,
			timeinfo->tm_min, timeinfo->tm_sec);

		sprintf(log_path, "%s/%s_log.txt", out_dir, date); 
		sprintf(summary_path, "%s/%s_summary.txt", out_dir, date); 
	}

	printf("Writing logs to: %s\n", log_path);
	printf("Writing summary to: %s\n", summary_path);
	log_file = fopen(log_path, "w");
	summary_file = fopen(summary_path, "w");

	if (log_file == NULL)
		die("Cannot open output file: %s", log_path);

 	// Add this information to the summary file.
 	if (runtime_limit_hrs_set == true) {
 		fprintf(summary_file, "time for MCMC run (hours): %i hours\n""sample every: %i\n"
		"log file path: %s/%s_log.txt",
		runtime_limit_hrs, sample_rate, out_dir, date);
 	} else {
 		fprintf(summary_file, "chain length: %i\n""sample every: %i\n"
		"log file path: %s/%s_log.txt",
		chain_length, sample_rate, out_dir, date);
 	}
}

int main(int argc, char **argv) {

	Decimal L = 0;
	Decimal mu = 0.5;
	Decimal indel_rate = 15;
	// indel_rate = 0.007575;

	int s, n, chain;

	init_rand();
	setup_gsl();
	hello_parse_cmdline(argc, argv);

	// I can read in the sequences as previously using Isaac's code.
	char **sequences;
	int sequences_cap;

	printf("Loading sequences...\n");
	int N = load_seqs(sequence_path, &sequences, &sequences_cap);
	assert(N > 0);
	printf("Loaded %i sequences.\n", N);

	// Check all sequences have the same length.
	size_t num_bases = strlen(sequences[0]);

	if (num_bases % 3 != 0)
		die("Sequences contain partial codons [%zu mod 3 != 0].", num_bases);

	for (n = 1; n < N; n++)
		if (strlen(sequences[n]) != num_bases)
			die("Reference sequences aren't all the same length.");

	printf("kappa:"DECPRINT"\n", kappa);
	printf("phi:"DECPRINT"\n", phi);

	int num_sites = num_bases / 3;
	assert(num_sites > 0);

	printf("Loaded %i sites.\n", num_sites);
	char (*sequences_codons)[num_sites] = my_malloc(N * sizeof(char[num_sites]), __FILE__, __LINE__);

	// This is currently done stupidly - but it's a rough approximation that we can improve upon.
	char consensus[num_bases+1];
	generate_consensus(sequences, N, num_bases, consensus);

	printf("consensus sequence:\n");
	for(s=0; s < (int)num_bases; s++) printf("%c", consensus[s]);
	printf("\n");


	// map_gaps_to_consensus(sequences, N, num_bases, consensus);
 	// Want to now map anything that is not '---' or an actual codon, to consensus.
	map_partial_gaps_to_consensus(sequences, N, num_sites, consensus);
	char *gap_present = my_calloc(num_sites, sizeof(char), __FILE__, __LINE__);

	for (n = 0; n < N; n++) {
		for (s = 0; s < num_sites; s++) {
			sequences_codons[n][s] = amino_to_code(sequences[n]+s*3);
			switch(sequences_codons[n][s]) {
				case 64: gap_present[s] = 1; break;
			}
		}
	}

	for (n = 0; n < N; n++) free(sequences[n]);
	
	free(sequences);

	Decimal *node_times = my_malloc((2*N-1) * sizeof(Decimal), __FILE__, __LINE__);
	int *daughter_1 = my_malloc((N-1) * sizeof(int), __FILE__, __LINE__);
	int *daughter_2 = my_malloc((N-1) * sizeof(int), __FILE__, __LINE__);
	load_from_json(tree_path, node_times, daughter_1, daughter_2, N);

	Decimal *branch1 = my_malloc((N-1) * sizeof(Decimal), __FILE__, __LINE__);
	Decimal *branch2 = my_malloc((N-1) * sizeof(Decimal), __FILE__, __LINE__);

	for (n = 0; n < N-1; n++) {	
		// Get the branch lengths.
		branch1[n] = node_times[n+N] - node_times[daughter_1[n]-1];
		branch2[n] = node_times[n+N] - node_times[daughter_2[n]-1];
	}

	Decimal *omega = my_malloc(num_sites * sizeof(Decimal), __FILE__, __LINE__);
	Decimal *tree_site_L = my_malloc(num_sites * sizeof(Decimal), __FILE__, __LINE__);
	Decimal *tmp_tree_site_L = my_malloc(num_sites * sizeof(Decimal), __FILE__, __LINE__);

  	Decimal init_omega = 1;
  	int num_hla_types = 1;
  	int h, i;
  	Selection *omega_windows = my_malloc(num_hla_types * sizeof(Selection),
                                        __FILE__, __LINE__);
	selec_window_alloc(num_sites, omega_windows);
	Decimal log_mu_sel_prior = log(mu_sel_prior);

	if(init_num_selec_windows > num_sites) init_num_selec_windows = num_sites;

	for(h = 0; h < num_hla_types; h++) {
		selec_windows_randomise(omega_windows[h].sel, num_sites, init_num_selec_windows,
			log_mu_sel_prior, sigma_sel_prior);
	}

	// Uncomment to initialise omega coefficients at 1 (used in simulations).
	for(h = 0; h < num_hla_types; h++) {
	  for(i = 0; i < init_num_selec_windows; i++) {
	    omega_windows[h].sel->windows[i].coeff = init_omega;
	  }
	}

	// Copy window coeff into array per site.
	update_HLA_selection(omega_windows[0].sel, num_sites, omega);
  	for (s = 0; s < num_sites; s++) {
  		tree_site_L[s] = 0.0;
  		tmp_tree_site_L[s] = 0.0;
  	}

  	Decimal (*rate_out_of_codons)[num_sites]
   		= my_malloc(NUM_CODONS * sizeof(Decimal[num_sites]),__FILE__, __LINE__);
   	Decimal (*codon_to_codon)[NUM_CODONS][num_sites]
    	= my_malloc(NUM_CODONS * sizeof(Decimal[NUM_CODONS][num_sites]), __FILE__, __LINE__);
    Decimal (*codon_to_codon_tmp)[NUM_CODONS][num_sites]
    	= my_malloc(NUM_CODONS * sizeof(Decimal[NUM_CODONS][num_sites]), __FILE__, __LINE__);

    int start_site = 0;
    int end_site = num_sites;

	tree_likelihood(start_site, end_site, mu, indel_rate, num_sites, omega, gap_present, N, sequences_codons, daughter_1, daughter_2,
		codon_to_codon, rate_out_of_codons, tree_site_L, branch1, branch2, true);
	
	for (s=0; s < num_sites; s++) L += tree_site_L[s];
	printf("Initial likelihood L: %f\n", L);

	switch(omega_prior){
		case EXPONENTIAL_PRIOR:
			printf("omega prior: exponential\n"
				   "lambda_sel: "DECPRINT"\n", lambda_sel);
			break;
		case LOG_NORMAL_PRIOR:
			printf("omega prior: log normal\n"
				   "mu_sel_prior = "DECPRINT"\n"
				   "sigma_sel_prior = "DECPRINT"\n", mu_sel_prior, sigma_sel_prior);
			break;
		default:
			die("No prior detected!");

	}

	// DEV: Can also include the Gnomad probabilities in there as parameters of the model.

	// Begin MCMC.
	time_t start, end;
	Decimal diff;
	time(&start);

	bool sel_on = true, mut_on = true, ind_on = true;
	bool sel_mergesplit_on = false;
	bool sel_grow_on = false;
	bool sel_wind_coeff_on = false;

	printf("window_mcmc %i\n", window_mcmc);
	printf("mutation only %i\n", mut_mcmc);
	printf("sample prior %i\n", sample_prior);

	if (mut_mcmc) {
		printf("only mutation moves selected.\n");
		sel_on = false;
		ind_on = false;
	}

	if (sel_mcmc) {
		printf("only selection moves selected.\n");
		mut_on = false;	
		ind_on = false;
	}

	if(ind_mcmc) {
		printf("only indel moves selected.\n");
		sel_on = false;
		mut_on = false;
	}

	if (window_mcmc) {
		printf("mutation moves and window moves selected.\n");
		sel_on = false;
		sel_mergesplit_on = true;
		sel_grow_on = true;
		sel_wind_coeff_on = true;
	}

	int mut, ind, sel, sel_mergesplit, sel_grow, sel_wind_coeff;
	
	sel = (sel_on ? 20 : 0);
	mut = (mut_on ? 1 : 0 ) + sel;
	ind = (ind_on ? 1 : 0 ) + mut;
	sel_mergesplit = (sel_mergesplit_on ? 7 : 0) + ind;
	sel_grow = (sel_grow_on ? 7 : 0) + sel_mergesplit;
	sel_wind_coeff = (sel_wind_coeff_on ? 7 : 0) + sel_grow;

	int which_move = sel_wind_coeff;
	int move;
	int sample_i = 1;

	// Burn in the mutation parameter first.
	int mutation_burn = 100; // DEV: add this as an option.
	for(chain=0; chain < mutation_burn; chain++)
	{
		mutation_move(&mu, indel_rate, num_sites, omega, gap_present, N, sequences_codons, daughter_1, daughter_2,
					&codon_to_codon, rate_out_of_codons, branch1, branch2, &codon_to_codon_tmp,
					&tree_site_L, &tmp_tree_site_L, &L);
	}

	fprintf(log_file, "state likelihood mu ");
	for(s = 0; s < num_sites; s++) {
		fprintf(log_file, "omega[%i] ", s);
	}
	fprintf(log_file, "\n");

	for (chain=0; chain < chain_length; chain++, sample_i++)
	{
		int choose_move = rand_lim(which_move);
		move = choose_move < sel ? 0 :
				(choose_move < mut ? 1 :
					(choose_move < ind ? 2 : 
						(choose_move < sel_mergesplit ? 3 :
							(choose_move < sel_grow ? 4 : 5))));

		// Window action.
		int waction;

		if (sample_i == sample_rate)
		{
			// Reset the sampling number.
			sample_i = 0;
			printf("Chain: %i\n", chain);
			// Print state to log.file.
			fprintf(log_file, "%i "DECPRINT" "DECPRINT" ", chain, L, mu);
			for (s = 0; s < num_sites; s++) {
				fprintf(log_file, DECPRINT" ", omega[s]);
			}
			fprintf(log_file, "\n");
			printf("Likelihood: %f, mu: %f, indel: %f\n", L, mu, indel_rate);
		}

		switch(move)
		{	
			case 0:
				selection_move(num_sites, mu, indel_rate, omega, gap_present, N, sequences_codons,
					daughter_1, daughter_2, codon_to_codon, rate_out_of_codons, tree_site_L,
					branch1, branch2, &L);
				break;
			case 1:
				mutation_move(&mu, indel_rate, num_sites, omega, gap_present, N, sequences_codons, daughter_1, daughter_2,
					&codon_to_codon, rate_out_of_codons, branch1, branch2, &codon_to_codon_tmp,
					&tree_site_L, &tmp_tree_site_L, &L);
				break;
			case 2: 
				indel_move(mu, &indel_rate, num_sites, omega, gap_present, N, sequences_codons, daughter_1, daughter_2,
					&codon_to_codon, rate_out_of_codons, branch1, branch2, &codon_to_codon_tmp,
					&tree_site_L, &tmp_tree_site_L, &L);
				break;
			case 3:

			case 4:

			case 5:
				// Choose which window action to perform.
				waction = move <= 3 ? DO_MERGE_OR_SPLIT
									: (move <= 4 ? DO_GROW : DO_COEFF);

				window_move(num_sites, omega, gap_present, omega_windows, waction, mu, indel_rate, N, sequences_codons,
					daughter_1, daughter_2, codon_to_codon, rate_out_of_codons, tree_site_L,
					branch1, branch2, &L, tmp_tree_site_L);

				break;
		}

		if (runtime_limit_hrs_set) {
			time(&end);
			diff = difftime(end,start);
			int diff_hours = (int) (diff / (60 * 60));
			if (diff_hours >= runtime_limit_hrs) {
				break;
			} else if (chain == (chain_length - 1)) {
				chain_length *= 2;
			}
		}
	}


	printf("Proposed %i moves.\n",chain);
	// Print out the various moves accepted and rejected during the MCMC.
	printf("selection accept: %i selection reject: %i\n", sel_accept, sel_reject);
	printf("mutation accept: %i mutation reject: %i\n", mut_accept, mut_reject);
	printf("indel accept: %i indel reject: %i\n", indel_accept, indel_reject);
	// Window moves
	printf("Window moves:\n");
	printf("merge accept: %i merge reject: %i\n", merge_accept, merge_reject);
	printf("split accept: %i split reject: %i\n", split_accept, split_reject);
	printf("grow accept: %i grow reject: %i\n", grow_accept, grow_reject);
	printf("window coeff accept: %i window coeff reject: %i\n", wcoeff_accept, wcoeff_reject);

	time(&end);
	diff = difftime(end,start);
	// Print the computation time in seconds.
	printf("time taken: "DECPRINTTIME"\n", diff);

	// mallocs in this file.
	free(daughter_1);
	free(daughter_2);
	free(branch1);
	free(branch2);
	free(node_times);
	free(sequences_codons);
	free(omega);
	free(tree_site_L);
	free(tmp_tree_site_L);
	free(rate_out_of_codons);
	free(codon_to_codon);
	free(codon_to_codon_tmp);
	free(gap_present);
	// mallocs in util.c
	selec_window_dealloc(omega_windows);
  	free(omega_windows);
  	
	fclose(log_file);
	fclose(summary_file);

	return EXIT_SUCCESS;
}
