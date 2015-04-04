#ifndef SPARSE_LIB_LINEARSOLVER_H
#define SPARSE_LIB_LINEARSOLVER_H

namespace  sparse_lib
{

template <typename T>
class LinearSolver {
public:
	LinearSolver() { }
	virtual ~LinearSolver() { }
	virtual int numIterations() const = 0;
	virtual void solve(const MatrixBase<T>* const, DenseVector<T>& x, const DenseVector<T>& b, bool setPreconditioner=true) = 0;
	virtual void solve(const MatrixBase<T>* const, DenseVector<T>& x, const DenseVectorView<T> b, bool setPreconditioner=true) = 0;
	virtual void solve(const MatrixBase<T>* const, DenseVector<T>& x, const DenseVectorTemp<T>& b, bool setPreconditioner=true) = 0;
	virtual void solve(const MatrixBase<T>* const, DenseVectorView<T> x, const DenseVector<T>& b, bool setPreconditioner=true) = 0;
	virtual void solve(const MatrixBase<T>* const, DenseVectorView<T> x, const DenseVectorView<T> b, bool setPreconditioner=true) = 0;
	virtual void solve(const MatrixBase<T>* const, DenseVectorView<T> x, const DenseVectorTemp<T>& b, bool setPreconditioner=true) = 0;

	virtual void setupPreconditioner(const MatrixBase<T>* const) = 0;
	
	virtual bool converged() const { return true; }
private:
	
};

} /*  sparse_lib */ 


#endif /* end of include guard: LINEARSOLVER_H */
