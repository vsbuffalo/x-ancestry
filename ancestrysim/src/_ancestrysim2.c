#include <Python.h>

#include "segments.h"
#include "genealogy.h"

#define UNUSED(x) (void)(x)


static char simtree_doc[] = "Simulate ancestry tree of an individual.";
static char ancestrysim2_doc[] = "C routines for simulating genetic ancestry.";
static char set_seed_doc[] = "Set GSL seed.";

PyObject *segment_t2PyTuple(segment_t *seg) {
  PyObject *tuple = Py_BuildValue("dd", seg->start, seg->end);
  return tuple;
}

PyObject *segments_t2PyList(segments_t *segs) {
  PyObject *list = PyList_New((Py_ssize_t) segs->len);
  segment_t *ptr;
  ptr = segs->head;
  int i = 0;
  while (ptr) {
    PyList_SetItem(list, i, segment_t2PyTuple(ptr));
    ptr = ptr->next;
    i++;
  }
  return list;
}

PyObject *individual_t2PyTuple(individual_t *ind) {
  /* convert an individual_t to a python tuple */
  PyObject *tuple, *asegs, *xsegs, *mx, *dx, *ma, *da;

  if (ind->asegments != NULL) {
      ma = segments_t2PyList(ind->asegments->mum);
      da = segments_t2PyList(ind->asegments->dad);
      asegs = Py_BuildValue("OO", ma, da);
      Py_DECREF(ma);
      Py_DECREF(da);
  } else {
     asegs = PyTuple_New(2);
     PyTuple_SetItem(asegs, 0, PyList_New(0));
     PyTuple_SetItem(asegs, 1, PyList_New(0));
  }
  if (ind->xsegments != NULL) {
      mx = segments_t2PyList(ind->xsegments->mum);
      dx = segments_t2PyList(ind->xsegments->dad);
      xsegs = Py_BuildValue("OO", mx, dx);
      Py_DECREF(mx);
      Py_DECREF(dx);
  } else {
     xsegs = PyTuple_New(2);
     PyTuple_SetItem(xsegs, 0, PyList_New(0));
     PyTuple_SetItem(xsegs, 1, PyList_New(0));
  }
  
  tuple = Py_BuildValue("IiiIiIOOIi", ind->id, ind->sex, ind->csex, ind->gen, ind->cid, 
                                      ind->anrec, asegs, xsegs, ind->xnrec,
                                      ind->is_x);
  Py_DECREF(asegs);
  Py_DECREF(xsegs);
  return tuple;
}

PyObject *tree_t2PyTuple(tree_t *tree) {
  PyObject *treelst, *this_gen, *res;
  treelst = PyList_New(tree->ngens+1);
  for (uint32_t i = 0; i < (tree->ngens+1); i++) {
    this_gen = PyList_New(tree->ninds[i]);
    for (uint32_t j = 0; j < tree->ninds[i]; j++) {
      PyList_SetItem(this_gen, j, individual_t2PyTuple(tree->tree[i][j]));
    }
    PyList_SetItem(treelst, i, this_gen);
  }
  res = PyTuple_New(3); 
  PyTuple_SetItem(res, 0, Py_BuildValue("I", tree->ngens)); 
  PyTuple_SetItem(res, 1, Py_BuildValue("i", tree->only_x)); 
  PyTuple_SetItem(res, 2, treelst);
  return res;
}

static PyObject *ancestrysim2_set_seed(PyObject *self, PyObject *args) {
  UNUSED(self);
  int ok;
  unsigned long int seed;
  ok = PyArg_ParseTuple(args, "I", &seed);
  if (!ok) return NULL;
  rng_init(seed);
  return Py_BuildValue("");
}

static PyObject *ancestrysim2_simtree(PyObject *self, PyObject *args) {
  UNUSED(self);
  int ok;
  unsigned ngens;
  int x, autos, only_x, sex; 
  PyObject *out;
  ok = PyArg_ParseTuple(args, "Iiiii", &ngens, &x, &autos, &sex, &only_x);
  if (!ok) return NULL;
  /* const char *p = gsl_rng_default->name; */
  /* fprintf(stderr, "RNG type: %s\n", p); */
  tree_t *sim = runsim(ngens, x, autos, sex, only_x);
  /* tree_print(sim); */
  out = tree_t2PyTuple(sim);
  tree_destroy(sim);
  return out;
}

static PyMethodDef module_methods[] = {
    {"_simtree", ancestrysim2_simtree, METH_VARARGS, simtree_doc},
    {"_set_seed", ancestrysim2_set_seed, METH_VARARGS, set_seed_doc},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_ancestrysim2(void) {
  PyObject *m = Py_InitModule3("_ancestrysim2", module_methods, ancestrysim2_doc);
  if (m == NULL) return;
}
