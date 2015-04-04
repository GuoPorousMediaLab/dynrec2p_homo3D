#ifndef SPARSE_LIB_SPARSELIB_H
#define SPARSE_LIB_SPARSELIB_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sparselib/types.h>
#include <sparselib/iterator.h>
#include <sparselib/densevector.h>
#include <sparselib/vvops.h>
#include <sparselib/givens.h>
#include <sparselib/matrixbase.h>
#include <sparselib/matrix.h>
#include <sparselib/matrix_engines/matrixengine.h>
#include <sparselib/mvops.h>

#include <sparselib/preconditioner.h>
#include <sparselib/preconditioners/diagonalpreconditioner.h>
#include <sparselib/preconditioners/incompletecholeskypreconditioner.h>
#include <sparselib/preconditioners/nopreconditioner.h>
#include <sparselib/preconditioners/additiveschwarzpreconditioner.h>

#include <sparselib/solvers/gmres_solver.h>
#include <sparselib/solvers/cg_solver.h>
#include <sparselib/solvers/schwarz_additive.h>

#endif /* end of include guard: SPARSE_LIB_SPARSELIB_H */
