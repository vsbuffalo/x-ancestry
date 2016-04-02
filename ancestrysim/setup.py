from distutils.core import setup, Extension

cfiles = ["src/segments.c", "src/genealogy.c", "src/_ancestrysim2.c"]
libs = ['m', 'gsl', 'clblas']

setup(
    ext_modules=[Extension("_ancestrysim2", cfiles, libraries=libs)],
)
