#include "load_from_json.h"
#include "util.h"
#include "constants.h"
// #include "selection_window.h"

#include "cJSON.h"
#include "string_buffer.h"
#include "amino_acids.h"

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <zlib.h>

// Want to also be able to read in a newick tree and do the same...

// Definition to check existence of classes within the .json file.
#define json_assert(x,msg,...) ((x) ? (void)0 : call_die(__FILE__,__LINE__,msg,##__VA_ARGS__))

void load_from_json(const char *path, Decimal *node_times, int *daughter_1, int *daughter_2, int N) 
{
  gzFile gz = gzopen(path, "r");
  if(gz == NULL) die("Cannot open file %s.", path);

  StrBuf sbuf;
  strbuf_alloc(&sbuf, 8192); // 2^13 bytes = 8KB

  // Read in the whole file.
  while(strbuf_gzreadline(&sbuf, gz)) {}
  gzclose(gz);

  cJSON *root, *node_times_in, *daughters_1_in, *daughters_2_in;

  json_assert(sbuf.end > 0, "Empty json file.");

  root = cJSON_Parse(sbuf.b);
  node_times_in = cJSON_GetObjectItem(root, "node_time");
  json_assert(node_times_in->type == cJSON_Array, "wut");

  int node, leaf;

  daughters_1_in = cJSON_GetObjectItem(root, "daughter_1");
  daughters_2_in = cJSON_GetObjectItem(root, "daughter_2");

  cJSON *daughter_1_in = daughters_1_in->child, *daughter_2_in = daughters_2_in->child; 

  // First get the number of nodes, we'll use this for checking.
  for(node = 0; daughter_1_in!= NULL; daughter_1_in = daughter_1_in->next, node++){}

  int N_in = node+1;
  if (N_in != N) die("The number of sequences in the .fasta does not match the sequence file.");

  daughter_1_in = daughters_1_in->child;

  for(leaf = 0; (daughter_1_in!= NULL) & (daughter_2_in!=NULL); daughter_1_in = daughter_1_in->next, daughter_2_in = daughter_2_in->next, leaf++) {
    daughter_1[leaf] = daughter_1_in->valueint;
    daughter_2[leaf] = daughter_2_in->valueint;
  }

  cJSON *node_time_in = node_times_in->child;
  for(node = 0; node_time_in!= NULL; node_time_in = node_time_in->next, node++) {
    node_times[node] = node_time_in->valuedouble;
  }

  cJSON_Delete(root);
  strbuf_dealloc(&sbuf);
}
