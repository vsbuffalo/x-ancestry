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

/* for experimental purposes only, should be 0 */
#define DEPEND_SEGMENTS 0
#define PRINT_RECOMB 0 
#define PRINT_IND 0

const double XGENLEN = 1.96;
const double AGENLEN = 35.239999999999995;
const double AUTOLENS[] = {1.79, 1.6, 1.73, 1.27, 1.17, 1.31, 1.35, 1.29, 
                           1.2, 1.08, 2.78, 1.08, 0.62, 0.74, 2.63, 2.25,
                           2.13, 2.04, 1.93, 1.87, 1.7, 1.68};
const uint32_t NAUTO = 22;

gsl_rng *rng;

#define EPS 0.0000001

/* chromosome checking must come earlier! */
#define segment_contains(segment, position) (position >= segment->start && \
    position <= segment->end)

#define FOCALIND(x) (**x->tree)


segpair_t *segpair_init_x() {
  /* assumes female */
  segpair_t *segpair = malloc(sizeof(segpair_t));
  segpair->mum = segments_init2(0, XGENLEN);
  segpair->dad = segments_init2(0, XGENLEN);
  return segpair;
}

segpair_t *segpair_init_auto() {
  segpair_t *segpair = malloc(sizeof(segpair_t));
  segpair->mum = segments_init();
  segpair->dad = segments_init();
  double sum = 0;
  for (uint32_t i = 0; i < NAUTO; i++) {
    segments_pushright3(segpair->mum, sum, sum+AUTOLENS[i], (char) i+1);
    segments_pushright3(segpair->dad, sum, sum+AUTOLENS[i], (char) i+1);
    sum += AUTOLENS[i];
  }
  return segpair;
}

void segpair_destroy(segpair_t *segpair) {
  if (segpair == NULL) {
    assert(false);
    return;
  }
  segments_destroy(segpair->mum);
  segments_destroy(segpair->dad);
  free(segpair);
}

uint32_t fib(uint32_t k) {
  double phi = (1 + sqrt(5))/2;
  double psi = (1 - sqrt(5))/2;
  return (uint32_t) floor((pow(phi, k) - pow(psi, k))/sqrt(5));
}

uint32_t nancestors(uint32_t k, bool is_x) {
  if (is_x) {
    return (uint32_t) fib(k+2);
  } else {
    return (uint32_t) pow(2, k);
  }
}

bool parent_is_x(individual_t *ind, bool parent_sex) {
  if (!ind->is_x || (ind->sex && parent_sex)) /* child is not x or two fathers*/
    return false;
  /*  */
  if (!ind->sex) return true; /* every parent of female is X ancestor */
  /* only male's mothers are X ancestors */
  if (ind->sex && !parent_sex) return true; 
  return false;
}

void ind_print2(individual_t *ind, FILE *stream) {
  if (!ind) {
    fprintf(stream, "null-ind"); return;
  }
  fprintf(stream, "(id=%d, k=%d, sex=%d, cid=%d, anrec=%d, xnrec=%d, ",
          ind->id, ind->gen, ind->sex, ind->cid, ind->anrec, ind->xnrec);
  if (ind->xsegments) {
    fprintf(stream, "mum_xsegs=");
    segments_print2(ind->xsegments->mum, stream);
    fprintf(stream, ", dad_xsegs=");
    segments_print2(ind->xsegments->dad, stream);
    fprintf(stream, ", ");
  } else {
    fprintf(stream, "mum_xsegs=[], dads_xsegs=[], ");
  }

if (ind->asegments) {
    fprintf(stream, "mum_asegs=");
    segments_print2(ind->asegments->mum, stream);
    fprintf(stream, ", dad_asegs=");
    segments_print2(ind->asegments->dad, stream);
  } else {
    fprintf(stream, "mum_asegs=[], dads_asegs=[]");
  }
  fprintf(stream, ")");
}


void ind_print(individual_t *ind) {
  ind_print2(ind, stdout);
}

void ind_printn(individual_t *ind) {
  ind_print(ind); printf("\n");
}

tree_t* tree_init(uint32_t k, bool only_x) {
  /* todo check mallocs */
  tree_t *tree = malloc(sizeof(tree_t));
  tree->tree = malloc(sizeof(individual_t**)*(k+1));
  tree->id = 0;
  tree->ninds = malloc(sizeof(uint32_t)*(k+1));
  tree->ptrs = malloc(sizeof(individual_t**)*(k+1));
  tree->only_x = only_x;
  tree->ngens = k;
  for (uint32_t i = 0; i < (k+1); i++) {
    tree->ninds[i] = nancestors(i, tree->only_x);
    tree->tree[i] = malloc(tree->ninds[i]*sizeof(individual_t*));
    tree->ptrs[i] = tree->tree[i];
  }
  return tree;
}

individual_t *ind_init(uint32_t id, uint32_t gen, bool sex, int8_t csex, 
        int32_t cid, uint32_t anrec, segpair_t *asegs, uint32_t xnrec,
        segpair_t *xsegs, bool is_x) {
  individual_t *ind = malloc(sizeof(individual_t));
  ind->id = id;
  ind->gen = gen;
  ind->sex = sex;
  ind->csex = csex;
  ind->cid = cid;
  ind->anrec = anrec;
  ind->xnrec = xnrec;
  ind->is_x = is_x;
  ind->asegments = asegs;
  ind->xsegments = xsegs;
  return ind;
}

void ind_destroy(individual_t *ind) {
  if (ind->xsegments)
    segpair_destroy(ind->xsegments);
  if (ind->asegments)
    segpair_destroy(ind->asegments);
  free(ind);
}

void recomb_print2(recomb_t *recomb, FILE *stream) {
  fprintf(stream, "Recomb(nbreaks=%d, segments=", recomb->nbreaks);
  segments_print2(recomb->segments, stream);
  fprintf(stream, ", states=[");
  for (int i = 0; (size_t) i < recomb->segments->len; i++) {
    fprintf(stream, "%d", recomb->states[i]);
    if (i < (int) recomb->segments->len-1)
      fprintf(stream, ", ");
  }
  fprintf(stream, "], breaks=[");
  for (uint32_t i = 0; i < recomb->nbreaks; i++) {
    fprintf(stream, "%.2f", recomb->breaks[i]);
    if (i < recomb->nbreaks-1)
      fprintf(stream, ", ");
  }
  fprintf(stream, "])"); 
}

void recomb_print(recomb_t *recomb) {
  recomb_print2(recomb, stdout);
}

recomb_t* recombine(segments_t *segments, double genlen, gsl_rng *rng) {
  /* ALL SEGMENTS MUST BE SORTED FOR THIS ROUTINE TO WORK CORRECTLY */
  assert(segments != NULL);
  recomb_t *recomb = malloc(sizeof(recomb_t));

  if (segments->len == 0) {
    fprintf(stderr, "segments.len must be > 0\n");
    exit(1);
  }
  uint8_t curr_state = (uint8_t) gsl_ran_bernoulli(rng, 0.5);
  uint32_t nbreaks = gsl_ran_poisson(rng, genlen);
  if (nbreaks == 0) {
    recomb->segments = segments_copy(segments);
    recomb->states = calloc(sizeof(bool), segments->len);
    recomb->breaks = NULL;
    recomb->nbreaks = 0;
    return recomb;
  }

  /* number of recomb breakpoints > 0 */
  double *breaks;
  breaks = malloc(sizeof(double)*nbreaks);
  for (uint32_t i = 0; i < nbreaks; i++) {
    breaks[i] = gsl_ran_flat(rng, 0.0, genlen);
  }
  gsl_sort(breaks, 1, nbreaks);
  
  /* nbreaks > 1 */ 
  recomb->nbreaks = nbreaks;
  int bi = 0, si = 0;  /* breakpoint index, state index */
  double curr_break = breaks[bi++];
  bool need_breakpoint = false;
  /* allocate max necessary */
  bool *new_states = malloc(sizeof(bool)*(nbreaks+segments->len+1));
  segments_t *new_segs, *process_segs = segments_copy(segments);
  segment_t *seg;
  seg = segments_popleft(process_segs); 
  char last_chrom = seg->chrom;
  new_segs = segments_init();
  nbreaks--; /* for counting purposes */

  while (true) {
    if (need_breakpoint || curr_break < seg->start) {
      /* either we need a new breakpoint because we've used the last one,
         or the current breakpoint is before the current segment start. */
        curr_state = 1 - curr_state;
      if (nbreaks == 0) {
        /* no breakpoints left, so all segments are are pushed using current 
           state */
        segments_pushright(new_segs, seg);
        new_states[si++] = curr_state; 
        while (process_segs->len) {
          /* push all segments to new_segs, flipping state at each new chromosome */
          last_chrom = seg->chrom;
          seg = segments_popleft(process_segs);
          if (!DEPEND_SEGMENTS && seg->chrom != '\0' && seg->chrom != last_chrom) {
            curr_state = (uint8_t) gsl_ran_bernoulli(rng, 0.5);
          }
          segments_pushright(new_segs, seg);
          new_states[si++] = curr_state;
        }
        break;
      } else {
        /* we do have breakpoints still, so swap state and update curr_break */
        curr_break = breaks[bi++];
        /* printf("just got break: %.5f, bi=%d, nbreaks=%d\n", curr_break, bi, nbreaks); */
        nbreaks--;
        need_breakpoint = false;
      }
    } else if (segment_contains(seg, curr_break)) {
      /* our current segment contains a breakpoint, so we need to split it,
         push the left half to the new segments, and keep processing the right
         half as it may have other breakpoints in it
         */
      segments_pushright3(new_segs, seg->start, curr_break, seg->chrom); /* push left seg */
      seg->start = curr_break; /* make current segment right seg */
      /* could make this efficient through less unnecessary initiation */
      new_states[si++] = curr_state;
      need_breakpoint = true;
    } else if (seg->end < curr_break) {
      /* out current segment's end is past the breakpoint, so we need to 
         push the current segment and its state, and then try to get get 
         a new segment. If this segment is on a difference chromosome, we
         set it's state to a random Bernoulli draw (independent assortment)  */
      segments_pushright(new_segs, seg);
      new_states[si++] = curr_state;
      if (!process_segs->len) {
        /* we're out of segments -- but could still have breakpoints */
        break;
      } else {
        /* we have segments, */
        last_chrom = seg->chrom;
        seg = segments_popleft(process_segs);
        if (!DEPEND_SEGMENTS && seg->chrom != '\0' && 
            seg->chrom != last_chrom) {
            /* new chromosome, need to draw random state */
            curr_state = (uint8_t) gsl_ran_bernoulli(rng, 0.5);
        }
      }
    } else {
      /* printf("seg: "); */
      /* segment_printn(seg); */
      /* printf("break %.5f\ncontains: %d\n", curr_break, segment_contains(seg, curr_break)); */

      /* fprintf(stderr, "execution should not reach this point"); */
      exit(1);
    }
  }

  if (process_segs->len != 0) {
    fprintf(stderr, "process_segs is not empty!\n");
    exit(1);
  }

  recomb->breaks = breaks;
  /* drop off unused states */
  new_states = realloc(new_states, sizeof(bool)*new_segs->len);
  recomb->states = new_states;
  recomb->segments = new_segs;
  double tns = total_seglength(new_segs),
         tos = total_seglength(segments);
  if (fabs(tns - tos) >= EPS) {
    fprintf(stderr, "total genetic length of new segments (%.3f) != total genetic length of old segments (%.3f)", tns, tos);
    segments_printn(new_segs);
    segments_printn(segments);
    exit(1);
  }
  segments_destroy(process_segs);
  return recomb;
}

void recomb_destroy(recomb_t *recomb) {
  if (recomb == NULL) return;
  segments_destroy(recomb->segments);
  free(recomb->states);
  free(recomb->breaks);
  free(recomb); 
}

segpair_t *meiosis(segments_t *segs, double genlen, gsl_rng *rng) {
  if (segs == NULL || segs->len == 0) return NULL;
  segpair_t *segpair = malloc(sizeof(segpair_t));
  recomb_t *recomb = recombine(segs, genlen, rng);
  if (PRINT_RECOMB) {
    recomb_print2(recomb, stderr); fprintf(stderr, "\n");
  }
  segpair->mum = segments_init();
  segpair->dad = segments_init();
  segments_t *ar[2] = {segpair->mum, segpair->dad};
  
  int si = 0;
  while (recomb->segments->len) {
    /* assign segment to appropriate mum/dad haplotype based on recomb->state */
    segments_pushright(ar[recomb->states[si]], segments_popleft(recomb->segments));
    si++;
  }
  assert(recomb->segments->len == 0);
  recomb_destroy(recomb);
  /* printf("sl: %.4f, mum: %.4f, dad: %.4f\n", total_seglength(segs), total_seglength(segpair->mum), total_seglength(segpair->dad)); */
  /* check we don't lose any material */
  assert(fabs(total_seglength(segs) - total_seglength(segpair->mum) - total_seglength(segpair->dad)) < EPS);
  return segpair;
}

void ancestry_recursion(individual_t *ind, uint32_t k, tree_t *tree,
        gsl_rng *rng) {
  /* this is a generation ahead k+2, and we need k+2 gens*/
  if (tree->ngens <= (k-1)) return; 
  if (PRINT_IND) ind_printn(ind);
  individual_t *mum, *dad;
  segpair_t *ma=NULL, *da=NULL, *mx=NULL, *dx=NULL;
  bool is_mum_x, is_dad_x;
  bool only_x = tree->only_x;
  /* mother */
  is_mum_x = parent_is_x(ind, 0);
  if (!only_x || is_mum_x) {
    tree->id++; 
    if (is_mum_x) {
      if (ind->xsegments)
        mx = meiosis(ind->xsegments->mum, XGENLEN, rng);
    }
    if (ind->asegments)
      ma = meiosis(ind->asegments->mum, AGENLEN, rng);
    mum = ind_init(tree->id, k, 0, ind->sex, (int32_t) ind->id, ind->anrec+1, ma,
                   ind->xnrec+is_mum_x, mx, is_mum_x);
    *tree->ptrs[k] = mum;
    tree->ptrs[k]++;
    ancestry_recursion(mum, k+1, tree, rng);
  }

  /* father */
  is_dad_x = parent_is_x(ind, 1);
  if (!only_x || is_dad_x) {
    tree->id++; 
    if (is_dad_x) {
      if (ind->xsegments) {
        dx = malloc(sizeof(segpair_t)); 
        dx->mum = segments_copy(ind->xsegments->dad);
        dx->dad = segments_init();
      }
    }
    if (ind->asegments)
      da = meiosis(ind->asegments->dad, AGENLEN, rng);
    dad = ind_init(tree->id, k, 1, ind->sex, (int32_t) ind->id, ind->anrec+1, da, 
                   ind->xnrec, dx, is_dad_x);
    /* ind_printn(dad); */
    *tree->ptrs[k] = dad;
    tree->ptrs[k]++;
    ancestry_recursion(dad, k+1, tree, rng);
  }
}

tree_t *runsim(uint32_t k, bool x, bool autos, bool sex, bool x_only) {
  /* female ind */
  segpair_t *xsegpair=NULL, *autopair=NULL;
  tree_t *tree;
  individual_t *ind;
  if (x) {
    xsegpair = segpair_init_x();
  }
  if (autos) {
    autopair = segpair_init_auto();
  }
  tree = tree_init(k, x_only); 
  ind = ind_init(0, 0, sex, -1, -1, 0, autopair, 0, xsegpair, 1);
  /* place the seed in tree */
  tree->tree[0][0] = ind;
  /* use global RNG */
  ancestry_recursion(tree->tree[0][0], 1, tree, rng);
  return tree;
}

void tree_print(tree_t *tree) {
  const char *logical[2] = { "false", "true"};
  printf("ngens=%d, only_x=%s\n", tree->ngens, logical[tree->only_x]);
  for (uint32_t i = 0; i < (tree->ngens+1); i++) {
    for (uint32_t j = 0; j < tree->ninds[i]; j++) {
      printf("%.*s %d: ", 4*i, " ", i); ind_printn(tree->tree[i][j]);
    }
    printf("\n");
  }
}

void tree_destroy(tree_t *tree) {
  for (uint32_t i = 0; i < (tree->ngens+1); i++) {
    for (uint32_t j = 0; j < tree->ninds[i]; j++) {
     if (tree->tree[i][j] != NULL)
       ind_destroy(tree->tree[i][j]); 
     else
       assert(false);
    }
    free(tree->tree[i]);
  }
  free(tree->ptrs);
  free(tree->ninds);
  free(tree->tree);
  free(tree);
}

void rng_init(unsigned long seed) {
  const gsl_rng_type *RNGT;
  RNGT = gsl_rng_default;
  /* RNGT = gsl_rng_ranlxs0; */
  rng = gsl_rng_alloc(RNGT);
  gsl_rng_set(rng, seed);
}

void rng_destroy() {
  gsl_rng_free(rng);
}
