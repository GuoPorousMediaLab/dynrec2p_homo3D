#ifndef SPARSE_LIB_SPARSEMATRIXENGINE_H
#define SPARSE_LIB_SPARSEMATRIXENGINE_H

namespace sparse_lib
{
	
template <typename T, MatrixType Tri, Symmetry Sym, MatrixStorage MatStorage>
class MatrixEngine;

template <typename T, MatrixType Tri=GENERAL, Symmetry Sym=ASYMMETRIC, MatrixStorage MatStorage=CRS>
class MatrixEngine
{
public:
	typedef T* iterator;
	typedef const T* const_iterator;
	typedef size_t size_type;
	typedef T value_type;
	
	MatrixEngine (): nRows_(0), nCols_(0), data_(0), limit_(0) { assert("Unsupported Sparse matrix engine." && false); }
	MatrixEngine (const size_type rows, const size_type cols, const size_type bw=0)
		: nRows_(rows), nCols_(cols) { assert("Unsupported Sparse matrix engine." && false); }
	~MatrixEngine () {  }
	
	inline size_type rows() const { assert("Unsupported Sparse matrix engine." && false); return nRows_; }
	inline size_type cols() const { assert("Unsupported Sparse matrix engine." && false); return nCols_; }
	
	value_type& operator()(const size_type r, const size_type c) { assert("Unsupported Sparse matrix engine." && false); return data_[c+r*nCols_]; }
	const value_type& operator()(const size_type r, const size_type c) const { assert("Unsupported Sparse matrix engine." && false); return data_[c+r*nCols_]; }

	void set(const size_type r, const size_type c, const value_type val) { assert("Unsupported Sparse matrix engine." && false); data_[c+r*nCols_] = val; }
	value_type& get(const size_type r, const size_type c) { assert("Unsupported Sparse matrix engine." && false); return data_[c+r*nCols_]; }
	const value_type& get(const size_type r, const size_type c) const { assert("Unsupported Sparse matrix engine." && false); return data_[c+r*nCols_]; }
	
	iterator data() { assert("Unsupported Sparse matrix engine." && false); return data_; }
	const_iterator data() const { assert("Unsupported Sparse matrix engine." && false); return data_; }
	
	template <typename V1, typename V2>
	void mv(const V1& x, V2& y, bool trans, T alpha, T beta) const { assert("Unsupported Sparse matrix engine." && false); }
	
	template <typename V1, typename V2>
	void solve(V1& x, const V2& rhs, bool trans) const { assert("Unsupported Sparse matrix engine." && false); }
	
	void print(std::ostream &o) const { assert("Unsupported Sparse matrix engine." && false); }
	
	void finalize() { assert("Unsupported Sparse matrix engine." && false); }

private:
	iterator data_;
	iterator limit_;
	size_type nRows_, nCols_;
};

template <typename E, typename S, typename T>
void DenseMatrixAllocate(E& TheEngine, const S size, const T val)
{
	std::allocator<T> alloc_;
	TheEngine.data_ = alloc_.allocate(size);
	TheEngine.limit_ = TheEngine.data_ + size;
	std::uninitialized_fill(TheEngine.data_, TheEngine.limit_, val);
}

template <typename E, typename S>
void DenseMatrixDeallocate(E& TheEngine, const S size)
{
	std::allocator<typename E::value_type> alloc_;
	if (TheEngine.data_) {
		typename E::iterator it = TheEngine.limit_;
		while (it != TheEngine.data_)
			alloc_.destroy(--it);

		alloc_.deallocate(TheEngine.data_, TheEngine.limit_ - TheEngine.data_);
	}
	TheEngine.data_ = TheEngine.limit_ = 0;
}

} /* sparse_lib */ 

#include <sparselib/matrix_engines/sparse_crs.tcc>
#include <sparselib/matrix_engines/sparse_cds.tcc>
#include <sparselib/matrix_engines/sparse_cds_sym_ut.tcc>
#include <sparselib/matrix_engines/sparse_cds_lower.tcc>

#include <sparselib/matrix_engines/dense_ge.tcc>
#include <sparselib/matrix_engines/dense_sy_upper.tcc>
#include <sparselib/matrix_engines/dense_tr_upper.tcc>

#include <sparselib/matrix_engines/sparse_schwarz.tcc>

#endif /* end of include guard: SPARSE_LIB_SPARSEMATRIXENGINE_H */
