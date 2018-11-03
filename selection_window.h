#ifndef SELECTION_WINDOW_H_
#define SELECTION_WINDOW_H_

#include <stdio.h>
#include "constants.h"

typedef struct
{
  int start, length;
  Decimal coeff;
} SelectionWindow;

typedef struct
{
  SelectionWindow *windows;
  int num_windows;
} WindowSet;

typedef struct
{
  WindowSet *sel;
} Selection;

void setup_gsl();
void clearup_gsl();

void selec_window_alloc(int num_sites, Selection *Selections);

void selec_window_dealloc(Selection *Selections);

void selec_windows_set_lengths(WindowSet *win_set, int num_sites);
void selec_windows_sort(WindowSet *win_set);

// Returns index of first (left-hand) window.
int selec_window_split_rnd(WindowSet *win_set, int num_sites,
                           Decimal *old_coeff);

void selec_window_split_undo(WindowSet *win_set, int windex, Decimal old_coeff);

// Returns index of new bigger window.
int selec_window_merge_rnd(WindowSet *win_set, int *old_start_site,
                           Decimal *old_coeff0, Decimal *old_coeff1);

void selec_window_merge_undo(WindowSet *win_set, int windex, int old_start_site,
                             Decimal old_coeff0, Decimal old_coeff1);

// Returns index of left hand window
// or -1 if move failed.
int selec_window_grow_rnd(WindowSet *win_set, int *old_start_site);
void selec_window_grow_undo(WindowSet *win_set, int windex, int old_start_site);

void selec_window_grow_accept(WindowSet *win_set, int new_start_site, 
                              int old_start_site, int num_sites);

int selec_window_change_coeff_rnd(WindowSet *win_set, Decimal *old_coeff0);
void selec_window_change_coeff_undo(WindowSet *win_set,
                                    int windex, Decimal old_coeff);

void selec_windows_randomise(WindowSet *win_set, int num_hla_types, int num_windows,
                             Decimal log_mu_sel_prior, Decimal sigma_sel_prior);

void selec_windows_set(WindowSet *win_set, int num_sites, int num_windows, int start[num_windows], 
                       Decimal coeff[num_windows]);

void selec_window_print(WindowSet *win_set);

void selec_window_print_to_file(WindowSet *win_set, FILE *print_to_file);

void update_HLA_selection(const WindowSet *winset, int num_sites, Decimal hla_select[num_sites]);

#endif /* SELECTION_WINDOW_H_ */
