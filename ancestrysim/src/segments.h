#pragma once 
#include <stddef.h>
#include <stdio.h>

typedef struct _segment_t {
  double start;
  double end;
  char chrom;
  struct _segment_t *prev;
  struct _segment_t *next;
} segment_t;

typedef struct {
  size_t len;
  segment_t *head;
  segment_t *tail;
} segments_t; 


void valid_segment(double start, double end);
segment_t *segment_init(double start, double end, char chrom);
segments_t *segments_init(void);
segments_t *segments_init2(double start, double end);
segments_t *segments_init3(double start, double end, char chrom);
void segments_pushright(segments_t *segs, segment_t *seg);
void segments_pushright2(segments_t *segs, double start, double end);
void segments_pushright3(segments_t *segs, double start, double end, char chrom);
void segments_pushleft(segments_t *segs, segment_t *seg);
void segments_pushleft2(segments_t *segs, double start, double end);
segment_t *segments_popleft(segments_t *segs);
segment_t *segments_popright(segments_t *segs);
void segment_print2(segment_t *seg, FILE *stream);
void segment_print(segment_t *seg);
void segment_printn(segment_t *seg);
void segments_print2(segments_t *segs, FILE *stream);
void segments_print(segments_t *segs);
void segments_printn(segments_t *segs);
segments_t *segments_copy(segments_t *segs);
void segment_destroy(segment_t *seg);
void segments_destroy(segments_t *segs);
double total_seglength(segments_t *segs);
