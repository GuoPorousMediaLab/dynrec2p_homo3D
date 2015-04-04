#ifndef SPARSE_LIB_TYPES_H
#define SPARSE_LIB_TYPES_H

namespace sparse_lib
{

enum Symmetry {
	SYMMETRIC,
	ASYMMETRIC
};

enum MatrixType {
	GENERAL,
	UPPER,
	LOWER
};

enum MatrixStorage {
	DENSE_FULL,
	DENSE_PACKED,
	CRS,
	CDS,
	SCHWARZ_BLOCK
};

enum SolverType {
	GMRES,
	CG,
	ADDITIVE_SCHWARZ,
	MULTIPLICATIVE_SCHWARZ,
	NO_SOLVER
};

} /* sparse_lib */ 

#endif /* end of include guard: SPARSE_LIB_TYPES_H */
