#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h> // Needed for rand().
#include <unistd.h> // Need for getpid() for getting setting rand number.
#include <limits.h> // To get PATH_MAX.
#include <assert.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "selection_window.h"
#include "util.h"
#include "constants.h"

gsl_rng *gsl_r;
double log_mu_sel_prior;

void setup_gsl()
{
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_r = gsl_rng_alloc(T);
  Decimal seed_rng = time (NULL) * getpid();
  gsl_rng_set (gsl_r, seed_rng);
  log_mu_sel_prior = log(mu_sel_prior);
}

void clearup_gsl()
{
  gsl_rng_free(gsl_r);
}

void selec_window_alloc(int num_sites,
                            Selection *Selections)
{

  WindowSet *window_sets = my_malloc(1 * sizeof(WindowSet),
                                     __FILE__,__LINE__);
  SelectionWindow *windows = my_malloc(num_sites * sizeof(SelectionWindow),
                                       __FILE__,__LINE__);

  int j;

  Selections[0].sel = &window_sets[0];

  // Set up for HLA[0].sel, hla[0].rev, hla[1].sel ...
  window_sets[0].num_windows = init_num_selec_windows;
  window_sets[0].windows = windows + num_sites * 0;

  for(j = 0; j < num_sites; j++)
    window_sets[0].windows[j].start = j;

}

void selec_window_dealloc(Selection *Selections)
{
  free(Selections[0].sel->windows);
  free(Selections[0].sel);
}

static int compare_selec_windows(const void *p1, const void *p2)
{
  const SelectionWindow *w1 = (const SelectionWindow *)p1;
  const SelectionWindow *w2 = (const SelectionWindow *)p2;
  return (w1->start - w2->start);
}

void selec_windows_sort(WindowSet *win_set)
{
  qsort(win_set->windows, win_set->num_windows, sizeof(SelectionWindow),
        compare_selec_windows);
}

static void add_rand_window(WindowSet *win_set, int num_sites)
{
  SelectionWindow tmp;
  int num_windows = win_set->num_windows;
  int index = num_windows + rand_lim(num_sites - num_windows);

  SWAP(win_set->windows[num_windows], win_set->windows[index], tmp);
}

static int bfind_window(const WindowSet *win_set,
                        const SelectionWindow *window)
{
  int mid, high = win_set->num_windows-1, low = 0;
  while(1)
  {
    mid = (low + high) / 2;
    if(win_set->windows[mid].start > window->start) high = mid - 1;
    else if(win_set->windows[mid].start + win_set->windows[mid].length
              <= window->start) low = mid + 1;
    else break;
    status("low: %i",low);
    status("high: %i",high);   
  }
  return mid;
}

// Returns index of first (left-hand) window.
int selec_window_split_rnd(WindowSet *win_set, int num_sites,
                           Decimal *old_coeff)
{
  if(win_set->num_windows == num_sites) die("This should not happen");

  status("run add_rand_window\n");
  add_rand_window(win_set, num_sites);
  status("done adding a window\n");
  
  SelectionWindow new_wind = win_set->windows[win_set->num_windows];

  // Find index to insert new window to.
  status("find window\n");
  int index = bfind_window(win_set, &new_wind);
  status("found window\n");
  status("index to put in new window: %i\n",index);

  // old_len gets split into len0, len1.
  *old_coeff = win_set->windows[index].coeff;
  int old_len = win_set->windows[index].length;
  int len0 = new_wind.start - win_set->windows[index].start;
  int len1 = old_len - len0;

  win_set->windows[index].length = len0;
  new_wind.length = len1;
  Decimal unif_rnd = rand_lim(1);

  Decimal ratio = (1 - unif_rnd) / unif_rnd;

  Decimal selec_prime_one = win_set->windows[index].coeff *
                            pow(ratio, -(Decimal)len1 / old_len);

  Decimal selec_prime_two = win_set->windows[index].coeff *
                            pow(ratio, (Decimal)len0 / old_len);

  win_set->windows[index].coeff = selec_prime_one;
  new_wind.coeff = selec_prime_two;

  size_t mem = (win_set->num_windows - index - 1) * sizeof(SelectionWindow);
  memmove(win_set->windows+index+2, win_set->windows+index+1, mem);
  win_set->windows[index+1] = new_wind;

  win_set->num_windows++;
  return index;
}

// windex is the index of the left hand window.
void selec_window_split_undo(WindowSet *win_set, int windex, Decimal old_coeff)
{
  win_set->windows[windex].length += win_set->windows[windex+1].length;
  win_set->windows[windex].coeff = old_coeff;

  // Get rid of windex+1 window.
  int old_start_site = win_set->windows[windex+1].start;
  size_t mem = (win_set->num_windows - windex - 2) * sizeof(SelectionWindow);
  memmove(win_set->windows+windex+1, win_set->windows+windex+2, mem);
  win_set->num_windows--;
  win_set->windows[win_set->num_windows].start = old_start_site;
}

// Returns index of new bigger window.
int selec_window_merge_rnd(WindowSet *win_set, int *old_start_site,
                           Decimal *old_coeff0, Decimal *old_coeff1)
{
  if(win_set->num_windows == 1) die("This should not happen");

  // Remove window at position 'index'.
  // Merge it into the window to the left.
  int index = 1 + rand_lim(win_set->num_windows-1);

  *old_start_site = win_set->windows[index].start;
  *old_coeff0 = win_set->windows[index-1].coeff;
  *old_coeff1 = win_set->windows[index].coeff;

  int len1 = win_set->windows[index-1].length;
  int len2 = win_set->windows[index].length;
  Decimal coeff1 = win_set->windows[index-1].coeff;
  Decimal coeff2 = win_set->windows[index].coeff;

  win_set->windows[index-1].coeff *= pow(coeff2 / coeff1, (Decimal)len2 / (len1+len2));
  win_set->windows[index-1].length += win_set->windows[index].length;

  size_t mem = (win_set->num_windows - index - 1) * sizeof(SelectionWindow);
  memmove(win_set->windows+index, win_set->windows+index+1, mem);
  win_set->num_windows--;
  win_set->windows[win_set->num_windows].start = *old_start_site;

  return index-1;
}

// windex is the index of the left hand window.
void selec_window_merge_undo(WindowSet *win_set, int windex, int old_start_site,
                             Decimal old_coeff0, Decimal old_coeff1)
{
  int wend = win_set->windows[windex].start + win_set->windows[windex].length;
  win_set->windows[windex].length = old_start_site - win_set->windows[windex].start;
  win_set->windows[windex].coeff = old_coeff0;

  // Move window back in.
  size_t mem = (win_set->num_windows - windex - 1) * sizeof(SelectionWindow);
  memmove(win_set->windows+windex+2, win_set->windows+windex+1, mem);
  win_set->num_windows++;
  win_set->windows[windex+1].start = old_start_site;
  win_set->windows[windex+1].coeff = old_coeff1;
  win_set->windows[windex+1].length = wend - old_start_site;
}

// Returns index of left hand window.
// or -1 if move failed.
int selec_window_grow_rnd(WindowSet *win_set, int *old_start_site)
{
  if(win_set->num_windows == 1) return -1;
  int index = 1 + rand_lim(win_set->num_windows-1);
  status("index: %i\n", index);
  int forward = (int)rand_lim(2) * 2 - 1;
  int geornd = forward * gsl_ran_geometric(gsl_r, selection_window_geom_param);
  
  if(geornd >= win_set->windows[index].length ||
     -geornd >= win_set->windows[index-1].length) return -1;

  *old_start_site = win_set->windows[index].start;

  win_set->windows[index-1].length += geornd;

  win_set->windows[index].start += geornd;

  win_set->windows[index].length += -geornd;
  
  status("\n");
  return index-1;
}

void selec_window_grow_accept(WindowSet *win_set, int new_start_site, 
                              int old_start_site, int num_sites)
{
  int i;
  for(i = win_set->num_windows; i < num_sites; i++) {
    if(win_set->windows[i].start == new_start_site) {
      win_set->windows[i].start = old_start_site;
      break;
    }
  }
}

// windex is index of left hand window.
void selec_window_grow_undo(WindowSet *win_set, int windex, int old_start_site)
{
  win_set->windows[windex].length = old_start_site - win_set->windows[windex].start;
  int wend = win_set->windows[windex+1].start + win_set->windows[windex+1].length;
  win_set->windows[windex+1].length = wend - old_start_site;
  win_set->windows[windex+1].start = old_start_site;
}

int selec_window_change_coeff_rnd(WindowSet *win_set, Decimal *old_coeff0)
{
  int w = rand_lim(win_set->num_windows);
  *old_coeff0 = win_set->windows[w].coeff;
  win_set->windows[w].coeff = *old_coeff0 * exp(rand_lim(2)-1);
  return w;
}

void selec_window_change_coeff_undo(WindowSet *win_set,
                                    int windex, Decimal old_coeff)
{
  win_set->windows[windex].coeff = old_coeff;
}

void selec_windows_set_lengths(WindowSet *win_set, int num_sites)
{
  int i;
  for(i = 0; i < win_set->num_windows-1; i++) {
    win_set->windows[i].length = win_set->windows[i+1].start -
                                 win_set->windows[i].start;
  }
  win_set->windows[win_set->num_windows-1].length
    = num_sites - win_set->windows[win_set->num_windows-1].start;
}

static void shuffle_windows(WindowSet *win_set, int num_sites, int num_windows)
{ 
  SelectionWindow tmp;
  int i, index;
  // Start from 1 so the first window is always starting from 0.
  for(i = 1; i < num_windows; i++) {
    index = i + rand_lim(num_sites-i);
    SWAP(win_set->windows[i], win_set->windows[index], tmp);
  }

  win_set->num_windows = num_windows;
  selec_windows_sort(win_set);
  selec_windows_set_lengths(win_set, num_sites);
}

void selec_windows_randomise(WindowSet *win_set, int num_sites, int num_windows,
                             Decimal log_mu_sel_prior, Decimal sigma_sel_prior)
{
  shuffle_windows(win_set, num_sites, num_windows);

  int i;
  for(i = 0; i < num_windows; i++)
  {
    win_set->windows[i].coeff = gsl_ran_lognormal(gsl_r, log_mu_sel_prior,
                                                  sigma_sel_prior);
  }
}

void selec_windows_set(WindowSet *win_set, int num_sites,
                       int num_windows, int start[num_windows], 
                       Decimal coeff[num_windows])
{
  SelectionWindow tmp;
  int i, index;
  // Start from 1 so the first window is always starting from 0.
  for(i = 1; i < num_windows; i++) {
    index = start[i];
    SWAP(win_set->windows[i], win_set->windows[index], tmp);
  }

  win_set->num_windows = num_windows;
  selec_windows_sort(win_set);

  for(i = 0; i < num_windows-1; i++) {
    win_set->windows[i].length = win_set->windows[i+1].start -
                                 win_set->windows[i].start;
  }
  win_set->windows[num_windows-1].length
    = num_sites - win_set->windows[num_windows-1].start;

  for(i = 0; i < num_windows; i++)
  {
    win_set->windows[i].coeff = coeff[i];
  }
}

void selec_window_print(WindowSet *win_set)
{
  int i;
  printf("Windows [%i]\n", win_set->num_windows);
  for(i = 0; i < win_set->num_windows; i++)
  {
    SelectionWindow *window = &win_set->windows[i];
    printf(" [%i,%i,"DECPRINT"]\n", window->start, window->length, window->coeff);
  }
}

void selec_window_print_to_file(WindowSet *win_set, FILE *print_to_file)
{
  int i;
  fprintf(print_to_file, "%i\n", win_set->num_windows);
  for(i = 0; i < win_set->num_windows - 1; i++)
  {
    SelectionWindow *window = &win_set->windows[i];
    fprintf(print_to_file, "%i,%i,"DECPRINT" | ", window->start, window->length, window->coeff);
  }
  SelectionWindow *window = &win_set->windows[i];
  fprintf(print_to_file, "%i,%i,"DECPRINT"\n", window->start, window->length, window->coeff);
  fflush(print_to_file);
}

void update_HLA_selection(const WindowSet *winset, int num_sites,
                                 Decimal hla_select[num_sites])
{
  int w, s, end;
  SelectionWindow *windows = winset->windows;
  for(w = 0; w < winset->num_windows; w++) {
    for(s = windows[w].start, end = s+windows[w].length; s < end; s++) {
      hla_select[s] = windows[w].coeff;
    }
  }
}
