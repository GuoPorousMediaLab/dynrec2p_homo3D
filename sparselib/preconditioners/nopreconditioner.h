#ifndef SPARSE_LIB_NOPRECONDITIONER_H
#define SPARSE_LIB_NOPRECONDITIONER_H

namespace sparse_lib
{

template <typename T, MatrixType Tri, Symmetry Sym, MatrixStorage MatStorage>
class NoPreconditionerBase: public Preconditioner<T>
{
public:
	typedef Matrix<MatrixEngine<T, Tri, Sym, MatStorage> > MatrixType;
	NoPreconditionerBase () {  }
	~NoPreconditionerBase () {  }
	
	void set(const MatrixBase<T>* Mat) { }
	void set(const Matrix<MatrixEngine<T, Tri, Sym, MatStorage> >& Mat) { }
	
	DenseVector<T>& apply(const DenseVector<T>& b, DenseVector<T>& x) { x=b; return x; }
	DenseVectorView<T>& apply(const DenseVector<T>& b, DenseVectorView<T>& x) { x=b; return x; }
    
	DenseVector<T>& apply(const DenseVectorView<T>& b, DenseVector<T>& x) { x=b; return x; }
	DenseVectorView<T>& apply(const DenseVectorView<T>& b, DenseVectorView<T>& x) { x=b; return x; }
	
	DenseVector<T>& apply(const DenseVectorTemp<T>& b, DenseVector<T>& x) { x=b; return x; }
	DenseVectorView<T>& apply(const DenseVectorTemp<T>& b, DenseVectorView<T>& x) { x=b; return x; }
	
	DenseVector<T>& applyTrans(const DenseVector<T>& b, DenseVector<T>& x) { x=b; return x; }
	DenseVectorView<T>& applyTrans(const DenseVector<T>& b, DenseVectorView<T>& x) { x=b; return x; }
    
	DenseVector<T>& applyTrans(const DenseVectorView<T>& b, DenseVector<T>& x) { x=b; return x; }
	DenseVectorView<T>& applyTrans(const DenseVectorView<T>& b, DenseVectorView<T>& x) { x=b; return x; }
	
	DenseVector<T>& applyTrans(const DenseVectorTemp<T>& b, DenseVector<T>& x) { x=b; return x; }
	DenseVectorView<T>& applyTrans(const DenseVectorTemp<T>& b, DenseVectorView<T>& x) { x=b; return x; }
	
private:
	/* data */
};

template <typename M>
class NoPreconditioner: public NoPreconditionerBase<typename M::value_type, M::MT, M::SYM, M::ST> {

};

} /* sparse_lib */ 


#endif /* end of include guard: SPARSE_LIB_NOPRECONDITIONER_H */
