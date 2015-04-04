#ifndef SPARSE_LIB_MVOPS_H
#define SPARSE_LIB_MVOPS_H

#include <sparselib/matrixbase.h>
#include <sparselib/densevector.h>

namespace sparse_lib
{
	
template <typename T, typename V>
DenseVectorTemp<T> operator*(const MatrixBase<T>& lhs, const V& rhs)
{
	DenseVectorTemp<T> tmp(rhs.size(), (T)0);
	lhs.mv(rhs, tmp);
	return tmp;
}

} /* sparse_lib */ 


#endif /* end of include guard: SPARSE_LIB_MVOPS_H */
