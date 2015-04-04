#ifndef SPARSE_LIB_MATRIX_H
#define SPARSE_LIB_MATRIX_H

#include <iostream>
#include <sparselib/matrixbase.h>
#include <sparselib/densevector.h>

namespace sparse_lib
{

template <typename Engine>
class Matrix: public MatrixBase<typename Engine::value_type>
{
public:
	typedef typename Engine::value_type T;
	typedef T* iterator;
	typedef const T* const_iterator;
	typedef size_t size_type;
	typedef typename Engine::value_type value_type;
	typedef Engine engine_type;
	static const MatrixType MT=Engine::MT;
	static const Symmetry SYM=Engine::SYM;
	static const MatrixStorage ST=Engine::ST;
	
	Matrix (): MatrixBase<T>() { }
	Matrix (size_type r, size_type c)
		: MatrixBase<T>(r, c), theEngine_(r, c) { }
	template <typename T2>
	Matrix (size_type r, size_type c, T2 extra)
		: MatrixBase<T>(r, c), theEngine_(r, c, extra) { }
	~Matrix () { }
	
	inline size_type rows() const { return theEngine_.rows(); }
	inline size_type cols() const { return theEngine_.cols(); }
	
	inline value_type& operator()(const size_type r, const size_type c) { return theEngine_(r, c); }
	inline const value_type& operator()(const size_type r, const size_type c) const { return theEngine_(r, c); }
	
	inline void set(const size_type r, const size_type c, const value_type val) { theEngine_.set(r, c, val); }
	
	inline void mv(const DenseVector<T>& x, DenseVectorTemp<T>& y, bool trans=false, T alpha=1.0, T beta=0.0) const { assert( trans ? (theEngine_.rows() == x.size() && theEngine_.rows() == y.size()) : (theEngine_.cols() == x.size() && theEngine_.cols() == y.size()) ); theEngine_.mv(x, y, trans, alpha, beta); }
	inline void mv(const DenseVectorView<T> x, DenseVectorTemp<T>& y, bool trans=false, T alpha=1.0, T beta=0.0) const { assert(trans ? (theEngine_.rows() == x.size() && theEngine_.rows() == y.size()) : (theEngine_.cols() == x.size() && theEngine_.cols() == y.size())); theEngine_.mv(x, y, trans, alpha, beta); }
	inline void mv(const DenseVectorTemp<T>& x, DenseVectorTemp<T>& y, bool trans=false, T alpha=1.0, T beta=0.0) const { assert(trans ? (theEngine_.rows() == x.size() && theEngine_.rows() == y.size()) : (theEngine_.cols() == x.size() && theEngine_.cols() == y.size())); theEngine_.mv(x, y, trans, alpha, beta); }

	inline void mv(const DenseVector<T>& x, DenseVectorView<T> y, bool trans=false, T alpha=1.0, T beta=0.0) const { assert(trans ? (theEngine_.rows() == x.size() && theEngine_.rows() == y.size()) : (theEngine_.cols() == x.size() && theEngine_.cols() == y.size())); theEngine_.mv(x, y, trans, alpha, beta); }
	inline void mv(const DenseVectorView<T> x, DenseVectorView<T> y, bool trans=false, T alpha=1.0, T beta=0.0) const { assert(trans ? (theEngine_.rows() == x.size() && theEngine_.rows() == y.size()) : (theEngine_.cols() == x.size() && theEngine_.cols() == y.size())); theEngine_.mv(x, y, trans, alpha, beta); }
	inline void mv(const DenseVectorTemp<T>& x, DenseVectorView<T> y, bool trans=false, T alpha=1.0, T beta=0.0) const { assert(trans ? (theEngine_.rows() == x.size() && theEngine_.rows() == y.size()) : (theEngine_.cols() == x.size() && theEngine_.cols() == y.size())); theEngine_.mv(x, y, trans, alpha, beta); }

	inline void solve(DenseVector<T>& x, const DenseVector<T>& rhs, bool trans=false) const { assert( theEngine_.rows()==theEngine_.cols() && theEngine_.cols() == x.size() && theEngine_.cols() == rhs.size() ); theEngine_.solve(x, rhs, trans); }
	inline void solve(DenseVector<T>& x, const DenseVectorView<T> rhs, bool trans=false) const { assert( theEngine_.rows()==theEngine_.cols() && theEngine_.cols() == x.size() && theEngine_.cols() == rhs.size() ); theEngine_.solve(x, rhs, trans); }
	inline void solve(DenseVector<T>& x, const DenseVectorTemp<T>& rhs, bool trans=false) const { assert( theEngine_.rows()==theEngine_.cols() && theEngine_.cols() == x.size() && theEngine_.cols() == rhs.size() ); theEngine_.solve(x, rhs, trans); }

	inline void solve(DenseVectorTemp<T>& x, const DenseVector<T>& rhs, bool trans=false) const { assert( theEngine_.rows()==theEngine_.cols() && theEngine_.cols() == x.size() && theEngine_.cols() == rhs.size() ); theEngine_.solve(x, rhs, trans); }
	inline void solve(DenseVectorTemp<T>& x, const DenseVectorView<T> rhs, bool trans=false) const { assert( theEngine_.rows()==theEngine_.cols() && theEngine_.cols() == x.size() && theEngine_.cols() == rhs.size() ); theEngine_.solve(x, rhs, trans); }
	inline void solve(DenseVectorTemp<T>& x, const DenseVectorTemp<T>& rhs, bool trans=false) const { assert( theEngine_.rows()==theEngine_.cols() && theEngine_.cols() == x.size() && theEngine_.cols() == rhs.size() ); theEngine_.solve(x, rhs, trans); }

	inline void solve(DenseVectorView<T> x, const DenseVector<T>& rhs, bool trans=false) const { assert( theEngine_.rows()==theEngine_.cols() && theEngine_.cols() == x.size() && theEngine_.cols() == rhs.size() ); theEngine_.solve(x, rhs, trans); }
	inline void solve(DenseVectorView<T> x, const DenseVectorView<T> rhs, bool trans=false) const { assert( theEngine_.rows()==theEngine_.cols() && theEngine_.cols() == x.size() && theEngine_.cols() == rhs.size() ); theEngine_.solve(x, rhs, trans); }
	inline void solve(DenseVectorView<T> x, const DenseVectorTemp<T>& rhs, bool trans=false) const { assert( theEngine_.rows()==theEngine_.cols() && theEngine_.cols() == x.size() && theEngine_.cols() == rhs.size() ); theEngine_.solve(x, rhs, trans); }

	inline void finalize() { theEngine_.finalize(); }
	
	virtual void print_matlab(std::ostream &o, std::string varName, int r0 = 0, int c0 = 0) const { theEngine_.print_matlab(o, varName, r0, c0); }
	
	inline engine_type& engine() { return theEngine_; }
	inline const engine_type& engine() const { return theEngine_; }
	
	void print() const { theEngine_.print(std::cout); }
	
protected:
	void print(std::ostream& o) const { theEngine_.print(o); }

private:
	engine_type theEngine_;

};

} /* sparse_lib */ 


#endif /* end of include guard: MATRIX_H */
