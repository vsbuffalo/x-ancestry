CC=clang
LDFLAGS=-lm -lgsl
CFLAGS=-g -std=c99 -pedantic -Werror -Wall -O3  -Wconversion
SDIR = src

default: all

segments.o: $(SDIR)/segments.c $(SDIR)/segments.h
main.o: $(SDIR)/main.c $(SDIR)/genealogy.h $(SDIR)/segments.h
genealogy.o: $(SDIR)/genealogy.c $(SDIR)/genealogy.h $(SDIR)/segments.h

OBJS=segments.o genealogy.o main.o 

%.o: $(SDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

tree: $(OBJS)
		$(CC) $(CFLAGS) $? -o tree $(LDFLAGS)

pyext: $(OBJS)
		python setup.py build_ext --inplace

all: pyext

clean:
		rm -f $(OBJS) ./tree 
