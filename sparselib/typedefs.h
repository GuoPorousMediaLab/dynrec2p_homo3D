#ifndef SPARSE_LIB_TYPEDEFS_H
#define SPARSE_LIB_TYPEDEFS_H

namespace sparse_lib
{

typedef DenseVector<double> vec;
typedef DenseVector<float> s_vec;

typedef DenseVectorView<double> vec_view;
typedef DenseVectorView<float> s_vec_view;

typedef Matrix<MatrixEngine<double, GENERAL, ASYMMETRIC, DENSE_FULL> > mat;
typedef Matrix<MatrixEngine<float, GENERAL, ASYMMETRIC, DENSE_FULL> > s_mat;

typedef Matrix<MatrixEngine<double, UPPER, SYMMETRIC, DENSE_FULL> > mat_sym;
typedef Matrix<MatrixEngine<float, UPPER, SYMMETRIC, DENSE_FULL> > s_mat_sym;

typedef Matrix<MatrixEngine<double, UPPER, ASYMMETRIC, DENSE_FULL> > mat_ut;
typedef Matrix<MatrixEngine<float, UPPER, ASYMMETRIC, DENSE_FULL> > s_mat_ut;

typedef Matrix<MatrixEngine<double> > sparse;
typedef Matrix<MatrixEngine<float> > s_sparse;

typedef Matrix<MatrixEngine<double, GENERAL, ASYMMETRIC, CDS> > sparse_cds;
typedef Matrix<MatrixEngine<float, GENERAL, ASYMMETRIC, CDS> > s_spasre_cds;

typedef Matrix<MatrixEngine<double, UPPER, SYMMETRIC, CDS> > sparse_sym_cds;
typedef Matrix<MatrixEngine<float, UPPER, SYMMETRIC, CDS> > s_sparse_sym_cds;

typedef Matrix<MatrixEngine<double, LOWER, ASYMMETRIC, CDS> > sparse_lower_cds;
typedef Matrix<MatrixEngine<float, LOWER, ASYMMETRIC, CDS> > s_sparse_lower_cds;

typedef Matrix<MatrixEngine<double, GENERAL, ASYMMETRIC, SCHWARZ_BLOCK> > sparse_block;
typedef Matrix<MatrixEngine<float, GENERAL, ASYMMETRIC, SCHWARZ_BLOCK> > s_sparse_block;


} /* sparse_lib */ 

#endif /* end of include guard: SPARSE_LIB_TYPES_H */
