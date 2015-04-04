#ifndef SPARSE_LIB_MATRIXBASE_H
#define SPARSE_LIB_MATRIXBASE_H

#include <algorithm>
#include <iostream>
#include <memory>
#include <functional>
#include <assert.h>
#include <iostream>
#include <sparselib/densevector.h>

namespace sparse_lib
{

template <typename T>
	class MatrixBase;

template<typename T>
std::ostream& operator<< (std::ostream& o, const MatrixBase<T>& M);

template <typename T>
class MatrixBase
{
public:
	typedef T* iterator;
	typedef const T* const_iterator;
	typedef size_t size_type;
	typedef T value_type;
	
	MatrixBase () { }
	MatrixBase (size_type r, size_type c) { }
	
	virtual ~MatrixBase () { }

	virtual size_type rows() const = 0;
	virtual size_type cols() const = 0;
	
	virtual value_type& operator()(const size_type, const size_type) = 0;
	virtual const value_type& operator()(const size_type, const size_type) const = 0;
	
	virtual void set(const size_type, const size_type, const value_type) = 0;
	
	virtual void mv(const DenseVector<T>& x, DenseVectorTemp<T>& y, bool trans=false, T alpha=1.0, T beta=0.0) const = 0;
	virtual void mv(const DenseVectorView<T> x, DenseVectorTemp<T>& y, bool trans=false, T alpha=1.0, T beta=0.0) const = 0;
	virtual void mv(const DenseVectorTemp<T>& x, DenseVectorTemp<T>& y, bool trans=false, T alpha=1.0, T beta=0.0) const = 0;

	virtual void mv(const DenseVector<T>& x, DenseVectorView<T> y, bool trans=false, T alpha=1.0, T beta=0.0) const = 0;
	virtual void mv(const DenseVectorView<T> x, DenseVectorView<T> y, bool trans=false, T alpha=1.0, T beta=0.0) const = 0;
	virtual void mv(const DenseVectorTemp<T>& x, DenseVectorView<T> y, bool trans=false, T alpha=1.0, T beta=0.0) const = 0;
	
	virtual void solve(DenseVectorTemp<T>& x, const DenseVector<T>& rhs, bool trans=false) const = 0;
	virtual void solve(DenseVectorTemp<T>& x, const DenseVectorView<T> rhs, bool trans=false) const = 0;
	virtual void solve(DenseVectorTemp<T>& x, const DenseVectorTemp<T>& rhs, bool trans=false) const = 0;
                       
	virtual void solve(DenseVectorView<T> x, const DenseVector<T>& rhs, bool trans=false) const = 0;
	virtual void solve(DenseVectorView<T> x, const DenseVectorView<T> rhs, bool trans=false) const = 0;
	virtual void solve(DenseVectorView<T> x, const DenseVectorTemp<T>& rhs, bool trans=false) const = 0;

	virtual void finalize() = 0;

	friend std::ostream& operator<< <>(std::ostream& o, const MatrixBase<T>& M);
	virtual void print_matlab(std::ostream &o, std::string varName, int r0 = 0, int c0 = 0) const = 0;
	
	virtual void print() const = 0;
protected:
	virtual void print(std::ostream& o) const = 0;
	
private:
	
	
	/* data */
};

template <typename T>
inline std::ostream& operator<< (std::ostream& o, const MatrixBase<T>& M)
{
	M.print(o);
	return o;
}

} /* sparse_lib */ 


#endif /* end of include guard: SPARSE_LIB_MATRIXBASE_H */
