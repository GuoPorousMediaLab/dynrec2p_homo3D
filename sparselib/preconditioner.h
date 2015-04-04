#ifndef SPARSE_LIB_PRECONDITIONER_H
#define SPARSE_LIB_PRECONDITIONER_H

namespace sparse_lib
{

template <typename T>
class Preconditioner
{
public:
	Preconditioner () { }
	virtual ~Preconditioner () { }
	
	virtual void set(const MatrixBase<T>* M) = 0;
	
	virtual DenseVector<T>& apply(const DenseVector<T>& b, DenseVector<T>& x) = 0;
	virtual DenseVectorView<T>& apply(const DenseVector<T>& b, DenseVectorView<T>& x) = 0;

	virtual DenseVector<T>& apply(const DenseVectorView<T>& b, DenseVector<T>& x) = 0;
	virtual DenseVectorView<T>& apply(const DenseVectorView<T>& b, DenseVectorView<T>& x) = 0;
	
	virtual DenseVector<T>& apply(const DenseVectorTemp<T>& b, DenseVector<T>& x) = 0;
	virtual DenseVectorView<T>& apply(const DenseVectorTemp<T>& b, DenseVectorView<T>& x) = 0;
	
	virtual DenseVector<T>& applyTrans(const DenseVector<T>& b, DenseVector<T>& x) = 0;
	virtual DenseVectorView<T>& applyTrans(const DenseVector<T>& b, DenseVectorView<T>& x) = 0;

	virtual DenseVector<T>& applyTrans(const DenseVectorView<T>& b, DenseVector<T>& x) = 0;
	virtual DenseVectorView<T>& applyTrans(const DenseVectorView<T>& b, DenseVectorView<T>& x) = 0;
	
	virtual DenseVector<T>& applyTrans(const DenseVectorTemp<T>& b, DenseVector<T>& x) = 0;
	virtual DenseVectorView<T>& applyTrans(const DenseVectorTemp<T>& b, DenseVectorView<T>& x) = 0;
private:
	
};

} /* sparse_lib */ 

#endif /* end of include guard: PRECONDITIONER_H */
