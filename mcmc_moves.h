#ifndef MCMCMOVES_H_
#define MCMCMOVES_H_

#include <stdio.h>
#include <stdbool.h>

#include "util.h"
#include "selection_window.h"

// Window move definitions.
#define DO_MERGE_OR_SPLIT 0
#define DO_MERGE 1
#define DO_SPLIT 2
#define DO_GROW 3
#define DO_COEFF 4

int sel_accept, sel_reject;
int merge_accept, merge_reject;
int split_accept, split_reject;
int grow_accept, grow_reject;
int wcoeff_accept, wcoeff_reject;
int mut_accept, mut_reject;
int indel_accept, indel_reject;

Decimal branch_length(Decimal *node_times, int from, int to);

Decimal max_vector(Decimal *x, int n);

void create_rate_and_codon_to_codon(const int num_sites, const int start_site,
	const int end_site, Decimal (*codon_to_codon)[NUM_CODONS][num_sites],
	Decimal (*rate_out_of_codons)[num_sites], const Decimal omega[num_sites]);

void tree_likelihood(const int start_site, const int end_site, Decimal mu, Decimal indel_rate, const int num_sites,
	Decimal omega[num_sites], const char gap_present[num_sites], const int N,
	char leaf_sequences[num_sites][N],int daughter_1[N-1], int daughter_2[N-1],
	Decimal (*codon_to_codon)[NUM_CODONS][num_sites], Decimal (*rate_out_of_codons)[num_sites],
	Decimal tree_site_L[num_sites], Decimal branch_1[N-1], Decimal branch_2[N-1], bool eval_rate_and_codon_to_codon);

void selection_move(int num_sites, Decimal mu, Decimal indel_rate, Decimal omega[num_sites], const char gap_present[num_sites],
	const int N, char leaf_sequences[num_sites][N], int daughter_1[N-1], int daughter_2[N-1],
	Decimal (*codon_to_codon)[NUM_CODONS][num_sites], Decimal (*rate_out_of_codons)[num_sites], Decimal tree_site_L[num_sites],
	Decimal branch1[N-1], Decimal branch2[N-1], Decimal *L);

Decimal evaluate_prior_diff_merge(int omega_prior,
								  Decimal prior_add_window,
                                  Decimal mu_sel_prior, Decimal sigma_sel_prior,
                                  Decimal old_coeff0, Decimal old_coeff1,
                                  Decimal new_coeff);

Decimal evaluate_prior_diff_split(int omega_prior,
								  Decimal prior_add_window,
                                  Decimal mu_sel_prior, Decimal sigma_sel_prior,
                                  Decimal old_coeff0,
                                  Decimal new_coeff1, Decimal new_coeff2);

Decimal evaluate_prior_diff_coeff(int omega_prior,
								  Decimal mu_sel_prior, Decimal sigma_sel_prior,
                                  Decimal old_coeff0, Decimal new_coeff);

void window_move(int num_sites,
				Decimal omega[num_sites], const char gap_present[num_sites],
				Selection *omega_windows, int window_action, Decimal mu, Decimal indel_rate, const int N,
				char leaf_sequences[num_sites][N], int daughter_1[N-1], int daughter_2[N-1],
				Decimal (*codon_to_codon)[NUM_CODONS][num_sites], Decimal (*rate_out_of_codons)[num_sites],
				Decimal tree_site_L[num_sites], Decimal branch1[N-1], Decimal branch2[N-1], Decimal *L,
				Decimal tmp_tree_site_L[num_sites]);

void mutation_move(Decimal *mu, Decimal indel_rate, int num_sites, Decimal omega[num_sites], const char gap_present[num_sites], 
	const int N, char leaf_sequences[num_sites][N], int daughter_1[N-1], int daughter_2[N-1],
	Decimal (**codon_to_codon_old)[NUM_CODONS][num_sites], Decimal (*rate_out_of_codons)[num_sites], 
	Decimal branch1[N-1], Decimal branch2[N-1], Decimal (**codon_to_codon_tmp)[NUM_CODONS][num_sites],
	Decimal **tree_site_L, Decimal **tmp_tree_site_L, Decimal *L);

void indel_move(Decimal mu, Decimal *indel_rate, int num_sites, Decimal omega[num_sites], const char gap_present[num_sites], 
	const int N, char leaf_sequences[num_sites][N], int daughter_1[N-1], int daughter_2[N-1],
	Decimal (**codon_to_codon_old)[NUM_CODONS][num_sites], Decimal (*rate_out_of_codons)[num_sites], 
	Decimal branch1[N-1], Decimal branch2[N-1], Decimal (**codon_to_codon_tmp)[NUM_CODONS][num_sites],
	Decimal **tree_site_L, Decimal **tmp_tree_site_L, Decimal *L);


#endif /* MCMCMOVES_H_ */