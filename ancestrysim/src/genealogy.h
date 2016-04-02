#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>

#include "segments.h"

extern const double XGENLEN;
extern const double AGENLEN;
extern const double AUTOLENS[];
extern const uint32_t NAUTO;
extern gsl_rng *rng;

typedef enum _sex_t { female = 0, male = 1 } sex_t;

typedef struct {
  segments_t *mum;
  segments_t *dad;
} segpair_t;


typedef struct {
  uint32_t id;
  int32_t cid;
  uint32_t gen;
  sex_t sex;
  int8_t csex; /* unlike sex, inds can have no children, sex = -1 */
  uint32_t anrec;
  segpair_t *asegments;
  /* if X segments tracked */
  segpair_t *xsegments;
  uint32_t xnrec;
  bool is_x;
} individual_t;

typedef struct {
  uint32_t nbreaks;
  double *breaks;
  bool *states;
  segments_t *segments;
} recomb_t;

typedef struct {
  individual_t ***tree;
  /* for adding new elements, pointer to last element inserted */
  individual_t ***ptrs; 
  uint32_t *ninds;
  uint32_t ngens;
  uint32_t id;
  double xgenlen;
  double agenlen;
  bool only_x;
} tree_t;

segpair_t *segpair_init_x(void);
segpair_t *segpair_init_auto(void);
void segpair_destroy(segpair_t *segpair);
uint32_t fib(uint32_t k);
uint32_t nancestors(uint32_t k, _Bool is_x);
_Bool parent_is_x(individual_t *ind, _Bool parent_sex);
void ind_print(individual_t *ind);
void ind_print2(individual_t *ind, FILE *stream);
void ind_printn(individual_t *ind);
tree_t *tree_init(uint32_t k, _Bool only_x);
individual_t *ind_init(uint32_t id, uint32_t gen, _Bool sex, int8_t csex, int32_t cid, uint32_t anrec, segpair_t *asegs, uint32_t xnrec, segpair_t *xsegs, _Bool is_x);
void ind_destroy(individual_t *ind);
void recomb_print2(recomb_t *recomb, FILE *stream);
void recomb_print(recomb_t *recomb);
recomb_t *recombine(segments_t *segments, double genlen, gsl_rng *rng);
void recomb_destroy(recomb_t *recomb);
segpair_t *meiosis(segments_t *segs, double genlen, gsl_rng *rng);
void ancestry_recursion(individual_t *ind, uint32_t k, tree_t *tree, gsl_rng *rng);
tree_t *runsim(uint32_t k, _Bool x, _Bool autos, _Bool sex, _Bool x_only);
void tree_print(tree_t *tree);
void tree_destroy(tree_t *tree);
void rng_init(unsigned long seed);
void rng_destroy(void);
