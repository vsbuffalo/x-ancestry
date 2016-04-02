#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_double.h>

#include "segments.h"

void valid_segment(double start, double end) {
  if (end < start || end < 0 || start < 0) {
    fprintf(stderr, "segment end must be >= segment start");
    exit(1);
  }
}

segment_t *segment_init(double start, double end, char chrom) {
  valid_segment(start, end);
  segment_t *seg = malloc(sizeof(segment_t));
  seg->start = start;
  seg->end = end;
  seg->chrom = chrom;
  seg->prev = NULL;
  seg->next = NULL;
  return seg;
}

segments_t *segments_init() {
  segments_t *segs = malloc(sizeof(segments_t));
  segs->len = 0;
  segs->head = NULL;
  segs->tail = NULL;
  return segs;
}

segments_t *segments_init2(double start, double end) {
  segments_t *segs = segments_init();
  segs->head = segment_init(start, end, '\0');
  segs->len = 1;
  segs->tail = segs->head;
  return segs;
}

segments_t *segments_init3(double start, double end, char chrom) {
  segments_t *segs = segments_init();
  segs->head = segment_init(start, end, chrom);
  segs->len = 1;
  segs->tail = segs->head;
  return segs;
}


void segments_pushright(segments_t *segs, segment_t *seg) {
  seg->prev = segs->tail;
  seg->next = NULL;
  if (segs->head == NULL) {
    segs->head = segs->tail = seg;
  } else {
    segs->tail->next = seg;
    segs->tail = seg;
  }
  segs->len++;
}

void segments_pushright2(segments_t *segs, double start, double end) {
  segment_t *seg = segment_init(start, end, '\0');
  segments_pushright(segs, seg);
}

void segments_pushright3(segments_t *segs, double start, double end, char chrom) {
  segment_t *seg = segment_init(start, end, chrom);
  segments_pushright(segs, seg);
}

void segments_pushleft(segments_t *segs, segment_t *seg) {
  seg->next = segs->head;
  seg->prev = NULL;
  if (segs->tail == NULL) {
    segs->head = segs->tail = seg;
  } else {
    segs->head->prev = seg;
    segs->head = seg;
  }
  segs->len++;
}

void segments_pushleft2(segments_t *segs, double start, double end) {
  segment_t *seg = segment_init(start, end, '\0');
  segments_pushleft(segs, seg);
}

segment_t *segments_popleft(segments_t *segs) {
  if (!segs->head) {
    if (segs->len != 0) {

      assert(false);
    }
    assert(segs->tail == NULL);
    return NULL;
  }
  segment_t *seg = segs->head;
  if (segs->len == 1)
    segs->head = segs->tail = NULL;
  else 
    segs->head = seg->next;
  seg->next = NULL;
  segs->len--;
  return seg;
}

segment_t *segments_popright(segments_t *segs) {
  if (!segs->tail) {
    assert(segs->len == 0);
    assert(segs->head == NULL);
    return NULL;
  }
  segment_t *seg = segs->tail;
  if (segs->len == 1)
    segs->head = segs->tail = NULL;
  else
    segs->tail = seg->prev;
  seg->prev = NULL;
  segs->len--;
  return seg;
}

void segment_print2(segment_t *seg, FILE *stream) {
  if (!seg) {
    fprintf(stream, "(,)");
  } else {
    if (seg->chrom == '\0')    
      fprintf(stream, "(%.2f, %.2f)", seg->start, seg->end);
    else
      fprintf(stream, "(%.2f, %.2f, %d)", seg->start, seg->end, seg->chrom);
  }
}

void segment_print(segment_t *seg) {
  segment_print2(seg, stdout);
}

void segment_printn(segment_t *seg) {
  segment_print(seg); printf("\n");
}

void segments_print2(segments_t *segs, FILE *stream) {
  if (segs == NULL) {
    fprintf(stream, "[]");
    return;
  }
  segment_t *seg = segs->head;
  fprintf(stream, "[");
  while (seg) {
    segment_print2(seg, stream);
    seg = seg->next;
    if (seg) fprintf(stream, ", ");
  }
  fprintf(stream, "]");
}

void segments_print(segments_t *segs) {
  segments_print2(segs, stdout);
}

void segments_printn(segments_t *segs) {
  segments_print(segs); printf("\n");
}

segments_t *segments_copy(segments_t *segs) {
  segments_t *segs_copy = segments_init();
  if (!segs || segs->len == 0) return segs_copy;
  segment_t *seg = segs->head;
  while (seg) {
    segments_pushright3(segs_copy, seg->start, seg->end, seg->chrom);
    seg = seg->next;
  }
  return segs_copy;
}

void segment_destroy(segment_t *seg) {
  free(seg);
}

void segments_destroy(segments_t *segs) {
  if (segs == NULL) return;
  segment_t *next, *seg = segs->head;
  int len = (int) segs->len;
  while (len) {
    next = seg->next;
    segment_destroy(seg);
    seg = next;
    len--;
  }
  free(segs);
}

double total_seglength(segments_t *segs) {
  segment_t *seg = segs->head;
  double seglength = 0;
  while (seg) {
    seglength += (seg->end - seg->start);
    seg = seg->next;
  }
  return seglength;
}

#ifdef DEBUG
int main(int argc, char *argv[]) {
  segment_t seg;
  segments_t *segs = segments_init();
  segments_t *sega = segments_init(), *segb = segments_init();
  segments_printn(segs);

  segments_pushleft2(sega, 2, 3);
  segments_pushleft2(sega, 4, 5);
  segments_printn(sega);
  segments_pushleft(segb, segments_popleft(sega));
  segments_printn(sega);
  segments_printn(segb);

  /* for (int i = 0; i < 12; i++) { */
  /*   if (i % 2) */
  /*     segments_pushright2(segs, i, i+1); */
  /*   else */
  /*     segments_pushleft2(segs, i, i+1); */
  /*   segments_printn(segs); */
  /* } */
  /* segments_t *copy = segments_copy(segs); */
  /* printf("copy: "); segments_printn(copy); */

  /* for (int i = 0; i < 12; i++) { */
  /*   if (i % 2) { */
  /*     printf("pop left:  "); */
  /*     segment_printn(segments_popleft(segs)); */
  /*   } else { */
  /*     printf("pop right:  "); */
  /*     segment_printn(segments_popright(segs)); */
  /*   } segments_printn(segs); */
  /* } */
  /* printf("copy: "); segments_printn(copy); */


  /* segments_popleft(segs); */
  /* segments_popright(segs); */
  /* segment_printn(segments_popright(segs)); */
  return 0;
}
#endif
