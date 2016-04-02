#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_double.h>

#include "segments.h"
#include "genealogy.h"



int main(int argc, char *argv[]) {
 const gsl_rng_type * T;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  tree_t *sim = runsim(5, 1, 0, 0, 0);
  tree_print(sim);
  tree_destroy(sim);

  /* for (int j = 0; j < 1000; j++) { */
  /* tree_t *sim = runsim(5, 1, 0, 0, 1); */
  /* segment_t *ptr; */
  /* tree_t *sim = runsim(5, 1, 0, 0, 1); */
  /* for (int i = 0; i < sim->ninds[5]; i++) { */
  /*   /1* ind_printn(sim->tree[5][i]); *1/ */
  /*   if (sim->tree[5][i]->xsegments) { */
  /*       ptr = sim->tree[5][i]->xsegments->mum->head; */
  /*       while (ptr) { */
  /*           printf("%d\t%d\tm\t%.4f\n", j, i, ptr->end - ptr->start); */
  /*           ptr = ptr->next; */
  /*       } */
  /*       ptr = sim->tree[5][i]->xsegments->dad->head; */
  /*       while (ptr) { */
  /*           printf("%d\t%d\tp\t%.4f\n", j, i, ptr->end - ptr->start); */
  /*           ptr = ptr->next; */
  /*       } */
  /*   } */
  /* } */
  /* tree_destroy(sim); */
  /* } */

/*   segpair_t *segpair = segpair_init_x(); */
/*   recomb_t *recomb = recombine(segpair->mum, XGENLEN, rng); */
/*   uint32_t num=300; */
/*   recomb_t *ptr[num]; */
/*   ptr[0] = recomb; */
/*   for (int i = 1; i < num; i++) { */
/*     recomb = recombine(recomb->segments, XGENLEN, rng); */
/*     ptr[i] = recomb; */
/*     recomb_print(recomb);printf("\n"); */
/*   } */
/*   for (int i = 0; i < num; i++) recomb_destroy(ptr[i]); */
/*   segpair_destroy(segpair); */


  /* /1* individual_t seed = {.id = 0, .sex=female, .cid=0, .nrec=0}; *1/ */
  /* genealogy_t genealogy; */
  /* ind_print(&seed); */
  /* for (int i = 0; i < 8; i++) { */
  /*   printf("i=%d, fib(i)=%d\n", i, fib(i)); */
  /* } */
  gsl_rng_free(rng);
  return 0;
}
