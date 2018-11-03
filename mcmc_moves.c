#include <stdio.h>  // Needed to perform IO operations
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>

#include "mcmc_moves.h"
#include "selection_window.h"
#include "util.h"

Decimal branch_length(Decimal *node_times, int from, int to)
{
  return node_times[to] - node_times[from];
}

Decimal max_vector(Decimal *x, int n){
  Decimal max_x = x[0];
  int i;
  for (i=1; i<n; i++){
    if (x[i] > max_x){
      max_x = x[i];
    }
  }
  return(max_x);
}

// Things to do:
// 1: Create a makefile.
// 2: Begin adding in more complex models.
// 3: Want to include some tests.
// 4: Remove constants and functions that are not required.
// 5. Decide whether we want to perform full blown MCMC or just some sort of EM or something.
// 6. Need to include the option of using distributions at the leaves of the tree - not hard.
// 7. Also need to think about indels. Not sure how to cope with them currently. Could model them?

// Think about using stuff like tensorflow for the problem....
// Why can't we just bung a ton of things into tensorflow? Need truth data though. Hmmm.

void create_rate_and_codon_to_codon(const int num_sites, const int start_site, const int end_site,
	Decimal (*codon_to_codon)[NUM_CODONS][num_sites], Decimal (*rate_out_of_codons)[num_sites],
	const Decimal omega[num_sites]) {

	int from, to, s;

	for (from = 0; from < NUM_CODONS; from++) {

		Decimal v1 = kappa * beta_S[from] + beta_V[from];
		Decimal v2 = kappa * alpha_S[from] + alpha_V[from];

		for (to = 0; to < NUM_CODONS; to++) {	
			for (s = start_site; s < end_site; s++) {
				rate_out_of_codons[from][s] = v1 + omega[s] * v2;
				if (rate_out_of_codons[from][s] < 1e-8) {
					codon_to_codon[from][to][s] =  prior_C2[to];
				} else {
				 	codon_to_codon[from][to][s]
				 		= phi * (omega[s] * (kappa * NS_TS[from][to] + NS_TV[from][to]) + 
				 			     kappa * S_TS[from][to] + S_TV[from][to]) / rate_out_of_codons[from][s] +
	                      (1-phi) * prior_C2[to];
            	}
			}
		}
	}
}

void tree_likelihood(const int start_site, const int end_site, Decimal mu, Decimal indel_rate, const int num_sites,
	Decimal omega[num_sites], const char gap_present[num_sites], const int N,
	char leaf_sequences[N][num_sites], int daughter_1[N-1], int daughter_2[N-1],
	Decimal (*codon_to_codon)[NUM_CODONS][num_sites], Decimal (*rate_out_of_codons)[num_sites],
	Decimal tree_site_L[num_sites], Decimal branch1[N-1], Decimal branch2[N-1], bool eval_rate_and_codon_to_codon) {

	// Should include some tests here.									
	if (daughter_1[0] == 0) {
		printf("Error: not passed the correct daughter vectors\n");
	}
	if (daughter_2[0] == 0) {
		printf("Error: not passed the correct daughter vectors\n");
	}

	int s, c, node, from, to;

	Decimal L_i[NUM_CODONS_AND_INDEL], L_j[NUM_CODONS_AND_INDEL];
	Decimal L_mat[2*N-1][NUM_CODONS_AND_INDEL];

	Decimal a[NUM_CODONS_AND_INDEL];
	Decimal b[NUM_CODONS_AND_INDEL];

   	Decimal P_change_branch_1[NUM_CODONS][num_sites];
   	Decimal P_change_branch_2[NUM_CODONS][num_sites];
   	Decimal P_codon_to_codon_branch_1[NUM_CODONS_AND_INDEL][NUM_CODONS_AND_INDEL];
   	Decimal P_codon_to_codon_branch_2[NUM_CODONS_AND_INDEL][NUM_CODONS_AND_INDEL];

   	// Additional variables  to consider indels.
   	Decimal P_indel_branch_1;
   	Decimal P_indel_branch_2;

   	Decimal P_no_int_indel_branch_1;
   	Decimal P_no_int_indel_branch_2;

   	Decimal prior_indel_prop_s;

   	if (eval_rate_and_codon_to_codon)
		create_rate_and_codon_to_codon(num_sites, start_site, end_site, codon_to_codon,
			rate_out_of_codons, omega);
	// printf("end site:%i\n", end_site);
	for (s = start_site; s < end_site; s++) {
		tree_site_L[s] = 0;
		prior_indel_prop_s = gap_present[s] ? prior_indel_prop : 0;
		// printf("site:%i, "DECPRINT" ", s, prior_indel_prop_s);
		// Set the likelihood at the leaves (observed information at the leaves).
		for (node = 0; node < N; node++) {
			for (c = 0; c < NUM_CODONS_AND_INDEL; c++) {
				L_mat[node][c] = 0;
				L_mat[node][(int)leaf_sequences[node][s]] = 1;
			}
		}

		// Loop over internal nodes
		for (node = 0; node < N-1; node++)
		{	
			// printf("node %i\n", node);
			// This should then be multiplied by the rows of the matrix codon_to_codon for each site.
			for (c = 0; c < NUM_CODONS_AND_INDEL; c++) {
				L_i[c] = L_mat[daughter_1[node]-1][c];
				L_j[c] = L_mat[daughter_2[node]-1][c];
			}

			// Probababilities of an indel.
			P_indel_branch_1 = 1 - exp(-omega[s] * indel_rate * branch1[node]);
			P_indel_branch_2 = 1 - exp(-omega[s] * indel_rate * branch2[node]);
			
			// Probability of no intermediate indels.
			P_no_int_indel_branch_1 = exp(-omega[s] * indel_rate * branch1[node] * prior_indel_prop_s);
			P_no_int_indel_branch_2 = exp(-omega[s] * indel_rate * branch2[node] * prior_indel_prop_s);

			// printf(DECPRINT"\n", P_no_int_indel_branch_1);
			// printf(DECPRINT"\n", P_no_int_indel_branch_2);
			// Then get the probability of going from all states to this state.
			for (from = 0; from < NUM_CODONS_AND_INDEL; from++) {

				a[from] = 0;
				b[from] = 0;

				Decimal checking_1 = 0;
				Decimal checking_2 = 0;

				if(from != 64) {

					P_change_branch_1[from][s] = 1 - exp(-rate_out_of_codons[from][s] * mu * branch1[node]);
					P_change_branch_2[from][s] = 1 - exp(-rate_out_of_codons[from][s] * mu * branch2[node]);

					for (to = 0; to < NUM_CODONS; to++) {
						P_codon_to_codon_branch_1[from][to] =
							(P_change_branch_1[from][s] * codon_to_codon[from][to][s] +
								(from == to) * (1 - P_change_branch_1[from][s])) * // Below is the additional indel contribution
								P_no_int_indel_branch_1 + 
								prior_C2[to] * (1 - prior_indel_prop_s * P_indel_branch_1 - P_no_int_indel_branch_1);
						P_codon_to_codon_branch_2[from][to] =
							(P_change_branch_2[from][s] * codon_to_codon[from][to][s] + 
								(from == to) * (1 - P_change_branch_2[from][s])) * // Below is the additional indel contribution
								P_no_int_indel_branch_2 + 
								prior_C2[to] * (1 - prior_indel_prop_s * P_indel_branch_2 - P_no_int_indel_branch_2);
						
						checking_1 += P_codon_to_codon_branch_1[from][to];
						checking_2 += P_codon_to_codon_branch_2[from][to];

						a[from] += P_codon_to_codon_branch_1[from][to] * L_i[to];
						b[from] += P_codon_to_codon_branch_2[from][to] * L_j[to];
					}
					// printf("checking 1: "DECPRINT"\n", checking_1);
					// printf("checking 2: "DECPRINT"\n", checking_2);
					// Now check to see if adding this extra indel piece gives me 1!
					P_codon_to_codon_branch_1[from][to] = prior_indel_prop_s * P_indel_branch_1;
					// checking_1 = P_codon_to_codon_branch_1[from][to];
					P_codon_to_codon_branch_2[from][to] = prior_indel_prop_s * P_indel_branch_2;
					// checking_2 = P_codon_to_codon_branch_2[from][to];
					// printf("checking 1: "DECPRINT"\n", log(P_codon_to_codon_branch_1[from][to]));
					// printf("checking 2: "DECPRINT"\n", log(P_codon_to_codon_branch_2[from][to]));
				} else {	
					// checking_1 = 0;
					// checking_2 = 0;
					for(to = 0; to < NUM_CODONS_AND_INDEL; to++) {

						if(from == to) {
							P_codon_to_codon_branch_1[from][to] = prior_indel_prop_s + (1 - prior_indel_prop_s) * (1 - P_indel_branch_1);
							P_codon_to_codon_branch_2[from][to] = prior_indel_prop_s + (1 - prior_indel_prop_s) * (1 - P_indel_branch_2);
						} else {
							P_codon_to_codon_branch_1[from][to] = prior_C2[to] * (1 - prior_indel_prop_s) * P_indel_branch_1;
							P_codon_to_codon_branch_2[from][to] = prior_C2[to] * (1 - prior_indel_prop_s) * P_indel_branch_2;
						}

						// checking_1 += P_codon_to_codon_branch_1[from][to];
						// checking_2 += P_codon_to_codon_branch_2[from][to];
						// printf("from: %i, to: %i, "DECPRINT"\n", from, to,P_codon_to_codon_branch_1[from][to]);
						a[from] += P_codon_to_codon_branch_1[from][to] * L_i[to];
						b[from] += P_codon_to_codon_branch_2[from][to] * L_j[to];
					}
					// printf("\nchecking_1:"DECPRINT"\n", checking_1);
					// printf("\nchecking_2:"DECPRINT"\n", checking_1);
				}
				L_mat[node+N][from] = a[from] * b[from];
			}
		}

		 // Lastly, need to sum over the prior on the root node.
		for (from = 0; from < NUM_CODONS; from++) {
			tree_site_L[s] += (prior_C2[from])/(1 + prior_indel_prop_s) * L_mat[2*N - 2][from];
		}
		// printf(DECPRINT"\n", tree_site_L[s]);
		// printf("hello "DECPRINT"\n", L_mat[2*N-2][from]);
		// printf("indel prior"DECPRINT"\n", prior_indel_prop_s/(1+prior_indel_prop_s));
		// printf("indel prior"DECPRINT"\n", prior_indel_prop_s);
		// printf("indel prior"DECPRINT"\n",prior_indel_prop_s/(1+prior_indel_prop_s) * L_mat[2*N - 2][from]);

		tree_site_L[s] += prior_indel_prop_s/(1+prior_indel_prop_s) * L_mat[2*N - 2][from];
		
		tree_site_L[s] = log(tree_site_L[s]);
	}
	// for(s=0; s<num_sites; s++) printf(DECPRINT" ", tree_site_L[s]);
	// printf("\n");
}

void selection_move(int num_sites, Decimal mu, Decimal indel_rate, Decimal omega[num_sites], const char gap_present[num_sites],
	const int N, char leaf_sequences[num_sites][N], int daughter_1[N-1], int daughter_2[N-1],
	Decimal (*codon_to_codon)[NUM_CODONS][num_sites], Decimal (*rate_out_of_codons)[num_sites],
	Decimal tree_site_L[num_sites], Decimal branch1[N-1], Decimal branch2[N-1], Decimal *L)             
{
	int rand_site = rand_lim(num_sites);

	// New omega.
	Decimal omega_prime = omega[rand_site] * exp(rand_lim(2)-1);
	// First store the old selection parameter.
	Decimal prev_omega = omega[rand_site];
	Decimal	prev_tree_site_L = tree_site_L[rand_site];

	omega[rand_site] = omega_prime;
	Decimal llk_ratio;
	Decimal prior_diff;

	if(sample_prior == false) {
		// Next, need to determine the new rate_out_of_codon and codon_to_codon at this site, and update the likelihood accordingly.
		// Only need evaluate the tree_likelhood at a single site - before and after.
		tree_likelihood(rand_site, rand_site+1, mu, indel_rate, num_sites, omega, gap_present, N, leaf_sequences, daughter_1, daughter_2,
			codon_to_codon, rate_out_of_codons, tree_site_L, branch1, branch2, true);

		llk_ratio = tree_site_L[rand_site] - prev_tree_site_L;
	} else {
		llk_ratio = 0;
	}

	switch(omega_prior)
	{
		case EXPONENTIAL_PRIOR:
			prior_diff = - lambda_sel * (omega_prime - prev_omega);
			break;
		case LOG_NORMAL_PRIOR:
			prior_diff = log(prev_omega / omega_prime) -
                   (pow(log(omega_prime / mu_sel_prior), 2) -
                    pow(log(prev_omega / mu_sel_prior), 2)) /
                   (2 * pow(sigma_sel_prior, 2));
           break;
        default:
        	die("No prior detected!");
	}
	
	Decimal acceptance_ratio = llk_ratio + prior_diff + log(omega_prime / prev_omega);
	Decimal alpha = exp(MIN2(0, acceptance_ratio));

	if(alpha * RAND_MAX > rand())
  	{
		// Accept the move.
		sel_accept++;
		*L = *L - prev_tree_site_L + tree_site_L[rand_site];
  	}
  	else 
  	{
  		// Reject the move.
  		sel_reject++;
  		// Reset the move.
  		omega[rand_site] = prev_omega;
  		// Reset the likelihood at the site.
  		tree_site_L[rand_site] = prev_tree_site_L;
  		// Reset codon_to_codon and rate_out_of_codons at the site.
  		create_rate_and_codon_to_codon(num_sites, rand_site, rand_site+1, codon_to_codon,
  			rate_out_of_codons, omega);
  	}
}

void mutation_move(Decimal *mu, Decimal indel_rate, int num_sites, Decimal omega[num_sites], const char gap_present[num_sites], 
	const int N, char leaf_sequences[num_sites][N], int daughter_1[N-1], int daughter_2[N-1],
	Decimal (**codon_to_codon_old)[NUM_CODONS][num_sites], Decimal (*rate_out_of_codons)[num_sites],
	Decimal branch1[N-1], Decimal branch2[N-1], Decimal (**codon_to_codon_tmp)[NUM_CODONS][num_sites],
	Decimal **tree_site_L_old, Decimal **tmp_tree_site_L, Decimal *L)
{
	int s;

	Decimal mu_prime = *mu * exp((rand_lim(1) - 0.5)/4);
	Decimal *tree_site_L_tmp = *tmp_tree_site_L;
	Decimal llk_ratio;
	Decimal L_new = 0;

	if(sample_prior == false) {
		Decimal (*tmp_codon_to_codon)[NUM_CODONS][num_sites] = *codon_to_codon_tmp;

		tree_likelihood(0, num_sites, mu_prime, indel_rate, num_sites, omega, gap_present, N, leaf_sequences,
			daughter_1, daughter_2, tmp_codon_to_codon, rate_out_of_codons, *tmp_tree_site_L, 
			branch1, branch2, true);

		
		for(s = 0; s < num_sites; s++) {
			L_new += tree_site_L_tmp[s];
		}

		llk_ratio =  L_new - *L;
	} else {
		llk_ratio = 0;
	}

	Decimal acceptance_ratio = llk_ratio - mu_prior * (mu_prime - *mu) + log(mu_prime / *mu);
	Decimal alpha = exp(MIN2((Decimal)0, acceptance_ratio));

	if(alpha * RAND_MAX > rand())
	{
		// Accept the move.
		mut_accept++;
		*mu = mu_prime;
		Decimal (*swap_codon_to_codon)[NUM_CODONS][num_sites];
		SWAP(*codon_to_codon_old, *codon_to_codon_tmp, swap_codon_to_codon);
		Decimal *swap_tree_site_L;
		SWAP(*tree_site_L_old, *tmp_tree_site_L, swap_tree_site_L);
		*L = L_new;
	}
	else
	{
		// Reject the move.
		mut_reject++;
	}
}

void indel_move(Decimal mu, Decimal *indel_rate, int num_sites, Decimal omega[num_sites], const char gap_present[num_sites], 
	const int N, char leaf_sequences[num_sites][N], int daughter_1[N-1], int daughter_2[N-1],
	Decimal (**codon_to_codon_old)[NUM_CODONS][num_sites], Decimal (*rate_out_of_codons)[num_sites],
	Decimal branch1[N-1], Decimal branch2[N-1], Decimal (**codon_to_codon_tmp)[NUM_CODONS][num_sites],
	Decimal **tree_site_L_old, Decimal **tmp_tree_site_L, Decimal *L)
{
	int s;

	// Decimal mu_prime = *mu * exp((rand_lim(1) - 0.5)/4);
	Decimal indel_rate_prime = *indel_rate * exp((rand_lim(1) - 0.5)/4);
	Decimal *tree_site_L_tmp = *tmp_tree_site_L;
	Decimal llk_ratio;
	Decimal L_new = 0;

	if(sample_prior == false) {
		Decimal (*tmp_codon_to_codon)[NUM_CODONS][num_sites] = *codon_to_codon_tmp;
		// printf("likelihood after the move\n");
		tree_likelihood(0, num_sites, mu, indel_rate_prime, num_sites, omega, gap_present, N, leaf_sequences,
			daughter_1, daughter_2, tmp_codon_to_codon, rate_out_of_codons, *tmp_tree_site_L, 
			branch1, branch2, true);

		
		for(s = 0; s < num_sites; s++) {
			L_new += tree_site_L_tmp[s];
		}

		llk_ratio =  L_new - *L;
	} else {
		llk_ratio = 0;
	}

	Decimal acceptance_ratio = llk_ratio - indel_rate_prior * (indel_rate_prime - *indel_rate) + log(indel_rate_prime / *indel_rate);
	Decimal alpha = exp(MIN2((Decimal)0, acceptance_ratio));

	if(alpha * RAND_MAX > rand())
	{
		// Accept the move.
		indel_accept++;
		*indel_rate = indel_rate_prime;
		Decimal (*swap_codon_to_codon)[NUM_CODONS][num_sites];
		SWAP(*codon_to_codon_old, *codon_to_codon_tmp, swap_codon_to_codon);
		Decimal *swap_tree_site_L;
		SWAP(*tree_site_L_old, *tmp_tree_site_L, swap_tree_site_L);
		*L = L_new;
	}
	else
	{
		// Reject the move.
		indel_reject++;
	}
}

Decimal evaluate_prior_diff_merge(int omega_prior,
                                  Decimal prior_add_window,
                                  Decimal mu_sel_prior, Decimal sigma_sel_prior,
                                  Decimal old_coeff0, Decimal old_coeff1,
                                  Decimal new_coeff)
{
	Decimal prior_diff;

	switch(omega_prior)
    {
        case LOG_NORMAL_PRIOR:
            // Log-normal prior ratio (and binomial for choosing splitting sites).
	       prior_diff = log(((1 - prior_add_window) * old_coeff0 * old_coeff1 * 
	                    sqrt(2 * M_PI) * sigma_sel_prior) / (new_coeff * prior_add_window)) 
	                  + (- pow(log(new_coeff / mu_sel_prior), 2) + 
	                    pow(log(old_coeff0 / mu_sel_prior), 2) +
	                    pow(log(old_coeff1 / mu_sel_prior), 2)) / (2 * pow(sigma_sel_prior, 2));
            break;
        case EXPONENTIAL_PRIOR:
            // Exponential prior ratio (and binomial for choosing splitting sites).
            prior_diff = log((1 - prior_add_window) / (lambda_sel * prior_add_window)) +
                        lambda_sel * ((old_coeff0 + old_coeff1) - new_coeff);
            break;
        default: die("Unknown prior!");
    }

	return(prior_diff);
}

Decimal evaluate_prior_diff_split(int omega_prior,
                                  Decimal prior_add_window,
                                  Decimal mu_sel_prior, Decimal sigma_sel_prior,
                                  Decimal old_coeff0,
                                  Decimal new_coeff0, Decimal new_coeff1)
{
	Decimal prior_diff;

    switch(omega_prior)
    {
        case LOG_NORMAL_PRIOR:
            // Log-normal prior ratio (and binomial for choosing splitting sites).
	        prior_diff = log((prior_add_window * old_coeff0) / 
	           ((1 - prior_add_window) * new_coeff0 * new_coeff1 * sqrt(2 * M_PI) * 
	             sigma_sel_prior)) 
	           + (- pow(log(new_coeff0 / mu_sel_prior), 2) - 
	              pow(log(new_coeff1 / mu_sel_prior), 2) + 
	              pow(log(old_coeff0 / mu_sel_prior), 2)) / (2 * pow(sigma_sel_prior, 2));
            break;
        case EXPONENTIAL_PRIOR:
            // Exponential prior ratio (and binomial for choosing splitting sites).
            prior_diff = log(lambda_sel * (prior_add_window) / (1 - prior_add_window)) + 
                        lambda_sel * (old_coeff0 - (new_coeff0 + new_coeff1));
            break;
        default: die("Unknown prior");
    }
  return(prior_diff);
}

Decimal evaluate_prior_diff_coeff(int omega_prior,
                                  Decimal mu_sel_prior, Decimal sigma_sel_prior,
                                  Decimal old_coeff0, Decimal new_coeff)
{
	Decimal prior_diff;
    switch(omega_prior)
    {
        case LOG_NORMAL_PRIOR:
            // Log-normal prior ratio (and binomial for choosing splitting sites).
            prior_diff = log(old_coeff0 / new_coeff) -
                       (pow(log(new_coeff / mu_sel_prior), 2) -
                        pow(log(old_coeff0 / mu_sel_prior), 2)) /
                        (2 * pow(sigma_sel_prior, 2));
            break;
        case EXPONENTIAL_PRIOR:
        	// Exponential prior ratio (and binomial for choosing splitting sites).
            prior_diff = lambda_sel * (old_coeff0 - new_coeff);
            break;
        default: die("Unknown prior!");
    }

	return(prior_diff);
}

void window_move(int num_sites,
                 Decimal omega[num_sites], const char gap_present[num_sites], 
                 Selection *omega_windows, int window_action, Decimal mu, Decimal indel_rate, const int N,
                 char leaf_sequences[num_sites][N], int daughter_1[N-1], int daughter_2[N-1],
                 Decimal (*codon_to_codon)[NUM_CODONS][num_sites], Decimal (*rate_out_of_codons)[num_sites],
                 Decimal tree_site_L[num_sites], Decimal branch1[N-1], Decimal branch2[N-1], Decimal *L,
                 // tmp variables.
                 Decimal tmp_tree_site_L[num_sites])
{
	int w, s;

	WindowSet *wset = omega_windows[0].sel;

	SelectionWindow *windows = wset->windows;
	Decimal *selection = omega;

	Decimal prob_split = 1, prob_merge = 1;
	if(window_action == DO_MERGE_OR_SPLIT)
	{
		// Choose between merge and split.
		int transition_pts = wset->num_windows - 1;
		prob_split = MIN2(1, (Decimal)(num_sites - wset->num_windows) /
		                             (Decimal)wset->num_windows *
		                              prior_add_window / (1-prior_add_window));
		status("PROB SPLIT: "DECPRINT"\n", prob_split);
		prob_merge = MIN2(1, (Decimal)transition_pts / (Decimal)(num_sites - transition_pts) *
		                             (1-prior_add_window) / prior_add_window);

		double urand = rand_lim(1);

		if(urand < prob_split / (prob_split + prob_merge))
			window_action = DO_SPLIT;
		else
			window_action = DO_MERGE;

	}

	int old_start_site;
	Decimal old_coeff0, old_coeff1;

	// Perform the move.
	switch(window_action)
	{	
		case DO_MERGE:
			status("merge move\n");
			status("NUMBER OF WINDOWS: %i\n", wset->num_windows);
			w = selec_window_merge_rnd(wset, &old_start_site, &old_coeff0, &old_coeff1);
			break;
		case DO_SPLIT:
			status("split move\n");
			status("NUMBER OF WINDOWS: %i\n", wset->num_windows);
			w = selec_window_split_rnd(wset, num_sites, &old_coeff0);
			break;
		case DO_GROW:
			status("grow move\n");
			w = selec_window_grow_rnd(wset, &old_start_site);
			if(w == -1) {
				grow_reject++;
				return;
			}
			break;
		case DO_COEFF:
			status("coeff move\n");
			w = selec_window_change_coeff_rnd(wset, &old_coeff0);
			break;
		default:
			die("No action!");
	}

	int start_site, end_site;

	if(window_action != DO_GROW) {
		start_site = windows[w].start;
		end_site = start_site + windows[w].length - 1;

		for(s = start_site; s <= end_site; s++) {
			tmp_tree_site_L[s] = tree_site_L[s];
			selection[s] = windows[w].coeff;
		}

		if(window_action == DO_SPLIT) {
			end_site += windows[w+1].length;
			for(; s <= end_site; s++) {
				tmp_tree_site_L[s] = tree_site_L[s];
				selection[s] = windows[w+1].coeff;
			}
		}
	} else {
		// The site which is moved.
		int new_start_site = windows[w+1].start;
		Decimal new_selection;

		if(old_start_site < new_start_site) {
			start_site = old_start_site;
			end_site = new_start_site -1;
			new_selection = windows[w].coeff;
		} else {
			start_site = new_start_site;
			end_site = old_start_site - 1;
			new_selection = windows[w+1].coeff;
		}

		for(s = start_site; s <= end_site; s++) {
			tmp_tree_site_L[s] = tree_site_L[s];
			selection[s] = new_selection;
		}
	}

	Decimal proposed_L = 0;
	Decimal current_L = 0;
	Decimal llk_ratio;

	if(sample_prior == false) {

		tree_likelihood(start_site, end_site+1, mu, indel_rate, num_sites, omega, gap_present, N, leaf_sequences, daughter_1, daughter_2,
			codon_to_codon, rate_out_of_codons, tree_site_L, branch1, branch2, true);

		for (s = start_site; s < end_site+1; s++) {
			proposed_L += tree_site_L[s];
			current_L += tmp_tree_site_L[s];
		} 

		llk_ratio = proposed_L - current_L;
	} else {
		llk_ratio = 0;
	}

	Decimal acceptance_ratio;
	Decimal coeff_sum, prior;
	Decimal prob_birth_in_dim_minus_one, prob_death, prob_birth_over_death;
	Decimal prob_death_in_dim_plus_one, prob_birth, prob_death_over_birth;
	Decimal prior_diff, new_coeff;
	Decimal new_coeff0, new_coeff1;

	switch(window_action)
	{
	case DO_MERGE:
	  // Note: the move has been made, so the number of windows is one less than the 
	  // 'current' number of windows. 
	  new_coeff = wset->windows[w].coeff;
	  prior = prior_add_window / (1 - prior_add_window);

	  status("PROB MERGE: "DECPRINT"\n", prob_merge);

	  prob_birth_in_dim_minus_one = ((num_sites - (Decimal)(wset->num_windows)) / 
	                                (Decimal)(wset->num_windows)) * prior;
	  prob_death = MIN2(1, 1 / prob_birth_in_dim_minus_one);

	  status("THE SAME AS PROB DEATH?: "DECPRINT"\n", prob_death);

	  prob_birth_in_dim_minus_one = MIN2(1, prob_birth_in_dim_minus_one);
	  prob_birth_over_death = prob_birth_in_dim_minus_one / prob_death;

	  prior_diff = evaluate_prior_diff_merge(omega_prior,
                                             prior_add_window,
	                                         mu_sel_prior, sigma_sel_prior,                                  
	                                         old_coeff0, old_coeff1, new_coeff);
	  
	  acceptance_ratio = llk_ratio + prior_diff +
	                     log(prob_birth_over_death *  
	                         ((wset->num_windows) * new_coeff) / 
	                         ((num_sites - (wset->num_windows)) * pow(old_coeff0 + old_coeff1, 2)));
	  break;
	case DO_SPLIT:
	  // Note: the move has been made, so the number of windows is one more than the 
	  // 'current' number of windows. 
	  new_coeff0 = wset->windows[w].coeff;
	  new_coeff1 = wset->windows[w+1].coeff;
	  status("NUMBER OF WINDOWS: %i\n", wset->num_windows);

	  coeff_sum = new_coeff0 + new_coeff1;
	  prior = (1 - prior_add_window) / prior_add_window;

	  status("PROB SPLIT: "DECPRINT"\n", prob_split);

	  prob_death_in_dim_plus_one = (Decimal)(wset->num_windows - 1) / 
	                               (Decimal)(num_sites - (wset->num_windows - 1)) * prior;
	  prob_birth = MIN2(1, 1 / prob_death_in_dim_plus_one);
	  
	  status("THE SAME AS PROB BIRTH?: "DECPRINT"\n", prob_birth);

	  prob_death_in_dim_plus_one  = MIN2(1, prob_death_in_dim_plus_one);
	  prob_death_over_birth = prob_death_in_dim_plus_one / prob_birth;

	  prior_diff = evaluate_prior_diff_split(omega_prior,
                                             prior_add_window,
	                                         mu_sel_prior, sigma_sel_prior,
	                                         old_coeff0, new_coeff0, new_coeff1);
	  
	  acceptance_ratio = llk_ratio + prior_diff +
	                     log(prob_death_over_birth *
	                         ((num_sites - (wset->num_windows - 1)) * pow(coeff_sum, 2)) / 
	                         ((wset->num_windows - 1) * old_coeff0));
	 break;
	case DO_GROW:
	  acceptance_ratio = llk_ratio;
	  break;
	case DO_COEFF:
	  new_coeff = wset->windows[w].coeff;

	  prior_diff = evaluate_prior_diff_coeff(omega_prior,
                                             mu_sel_prior, sigma_sel_prior,
                                             old_coeff0, new_coeff);

	  acceptance_ratio = llk_ratio + prior_diff +
	                     log(new_coeff / old_coeff0);
	  break;
	}

	Decimal alpha = exp(MIN2(0, acceptance_ratio));

	if(alpha * RAND_MAX > rand())
	{
		// Accept the move.
		switch(window_action) {
			case DO_MERGE: merge_accept++; break;
			case DO_SPLIT: split_accept++; break;
			case DO_GROW: 
			grow_accept++; 
			selec_window_grow_accept(wset, wset->windows[w+1].start, old_start_site, num_sites);
			break;
			case DO_COEFF: wcoeff_accept++; break;
		}
		// Accept, need to update certain entries of the likelihood over the sites.
		// Have things to fill in here.
		*L = *L + llk_ratio;
	}
	else
	{
		// Reject the move.
		switch(window_action) {
			case DO_MERGE: merge_reject++; break;
			case DO_SPLIT: split_reject++; break;
			case DO_GROW: grow_reject++; break;
			case DO_COEFF: wcoeff_reject++; break;
		}

		// Undo the window move which was applied.  
		switch(window_action) {
			case DO_MERGE:
				selec_window_merge_undo(wset, w, old_start_site, old_coeff0, old_coeff1);
				break;
			case DO_SPLIT:
				selec_window_split_undo(wset, w, old_coeff0);
				break;
			case DO_GROW:
				selec_window_grow_undo(wset, w, old_start_site);
				break;
			case DO_COEFF:
				selec_window_change_coeff_undo(wset, w, old_coeff0);
				break;
		}

  		int tmp_start_site = start_site;
		int tmp_end_site = end_site;

		start_site = windows[w].start;
		end_site = start_site + windows[w].length - 1;

		for(s = start_site; s <= end_site; s++)
			selection[s] = windows[w].coeff;

		if(window_action == DO_MERGE || window_action == DO_GROW)
		{	
			start_site = windows[w+1].start;
			end_site = start_site + windows[w+1].length - 1;
			for(s = start_site; s <= end_site; s++)
				selection[s] = windows[w+1].coeff;
		}

		create_rate_and_codon_to_codon(num_sites, tmp_start_site, tmp_end_site+1 , codon_to_codon,
			rate_out_of_codons, omega);

		for (s = tmp_start_site; s < tmp_end_site+1; s++) {
			tree_site_L[s] = tmp_tree_site_L[s];
		}
	}
}

