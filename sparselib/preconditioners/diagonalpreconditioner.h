#ifndef SPARSE_LIB_DIAGONALPRECONDITIONER_H
#define SPARSE_LIB_DIAGONALPRECONDITIONER_H

namespace sparse_lib
{

template <typename T, MatrixType Tri, Symmetry Sym, MatrixStorage MatStorage>
class DiagonalPreconditionerBase: public Preconditioner<T>
{
public:
	DiagonalPreconditionerBase () { assert("Unsupported matrix type." && false); }
	~DiagonalPreconditionerBase () { }
	
	void set(const MatrixBase<T>* Mat) { } 
	
	DenseVector<T>& apply(const DenseVector<T>& b, DenseVector<T>& x) { return x; }
	DenseVectorView<T>& apply(const DenseVector<T>& b, DenseVectorView<T>& x) { return x; }
    
	DenseVector<T>& apply(const DenseVectorView<T>& b, DenseVector<T>& x) { return x; }
	DenseVectorView<T>& apply(const DenseVectorView<T>& b, DenseVectorView<T>& x) { return x; }
	
	DenseVector<T>& apply(const DenseVectorTemp<T>& b, DenseVector<T>& x) { return x; }
	DenseVectorView<T>& apply(const DenseVectorTemp<T>& b, DenseVectorView<T>& x) { return x; }
	
	DenseVector<T>& applyTrans(const DenseVector<T>& b, DenseVector<T>& x) { return x; }
	DenseVectorView<T>& applyTrans(const DenseVector<T>& b, DenseVectorView<T>& x) { return x; }
    
	DenseVector<T>& applyTrans(const DenseVectorView<T>& b, DenseVector<T>& x) { return x; }
	DenseVectorView<T>& applyTrans(const DenseVectorView<T>& b, DenseVectorView<T>& x) { return x; }
	
	DenseVector<T>& applyTrans(const DenseVectorTemp<T>& b, DenseVector<T>& x) { return x; }
	DenseVectorView<T>& applyTrans(const DenseVectorTemp<T>& b, DenseVectorView<T>& x) { return x; }
	
private:
	/* data */
};

template <typename T, MatrixType Tri, Symmetry Sym>
class DiagonalPreconditionerBase<T, Tri, Sym, CDS>: public Preconditioner<T>
{
public:
	typedef Matrix<MatrixEngine<T, Tri, Sym, CDS> > MatrixType;
	DiagonalPreconditionerBase () { }
	~DiagonalPreconditionerBase () { }
	
	void set(const MatrixBase<T>* Mat);
	void set(const Matrix<MatrixEngine<T, Tri, Sym, CDS> >& Mat);
	
	inline DenseVector<T>& apply(const DenseVector<T>& b, DenseVector<T>& x) { return apply_(b, x); }
	inline DenseVectorView<T>& apply(const DenseVector<T>& b, DenseVectorView<T>& x) { return apply_(b, x); }
	inline DenseVector<T>& apply(const DenseVectorView<T>& b, DenseVector<T>& x) { return apply_(b, x); }
	inline DenseVectorView<T>& apply(const DenseVectorView<T>& b, DenseVectorView<T>& x) { return apply_(b, x); }
	inline DenseVector<T>& apply(const DenseVectorTemp<T>& b, DenseVector<T>& x) { return apply_(b, x); }
	inline DenseVectorView<T>& apply(const DenseVectorTemp<T>& b, DenseVectorView<T>& x) { return apply_(b, x); }
	inline DenseVector<T>& applyTrans(const DenseVector<T>& b, DenseVector<T>& x) { return apply_(b, x); }
	inline DenseVectorView<T>& applyTrans(const DenseVector<T>& b, DenseVectorView<T>& x) { return apply_(b, x); }
	inline DenseVector<T>& applyTrans(const DenseVectorView<T>& b, DenseVector<T>& x) { return apply_(b, x); }
	inline DenseVectorView<T>& applyTrans(const DenseVectorView<T>& b, DenseVectorView<T>& x) { return apply_(b, x); }
	inline DenseVector<T>& applyTrans(const DenseVectorTemp<T>& b, DenseVector<T>& x) { return apply_(b, x); }
	inline DenseVectorView<T>& applyTrans(const DenseVectorTemp<T>& b, DenseVectorView<T>& x) { return apply_(b, x); }
	
private:
	DenseVector<T> invDiag_;
	template <typename V1, typename V2>
	V2& apply_(const V1& b, V2& x);
	
};

template <typename T, MatrixType Tri, Symmetry Sym>
void DiagonalPreconditionerBase<T, Tri, Sym, CDS>::set(const MatrixBase<T>* Mat)
{
	const Matrix<MatrixEngine<T, Tri, Sym, CDS> >* m = dynamic_cast<const Matrix<MatrixEngine<T, Tri, Sym, CDS> >*>(Mat);
	if(m) {
		set(*m);
	} else {
		assert("Wrong Matrix type.  Must be a CDS." && false);
	}
}

template <typename T, MatrixType Tri, Symmetry Sym>
void DiagonalPreconditionerBase<T, Tri, Sym, CDS>::set(const Matrix<MatrixEngine<T, Tri, Sym, CDS> >& m)
{
	const int* diags = m.engine().diags();
	int nDiags = m.engine().nDiags();
	const T* values = m.engine().values();
	int nRows = m.rows();
	int loc=0;
	while(loc<nRows && *diags !=0) {
		loc++; diags++;
	}
	if(loc==nRows) {
		assert("No main diagonal.");
	}
	invDiag_.resize(nRows);
	for(int i=0; i<nRows; ++i) {
		assert(values[loc*nRows+i]!=0.0);
		invDiag_(i)=T(1)/values[loc*nRows+i];
	}
}

template <typename T, MatrixType Tri, Symmetry Sym> template <typename V1, typename V2>
V2& DiagonalPreconditionerBase<T, Tri, Sym, CDS>::apply_(const V1& b, V2& x)
{
	assert(b.size() == x.size() && x.size() == invDiag_.size());
	for(int i=0; i<b.size(); ++i)
		x(i)=b(i)*invDiag_(i);
	return x;
}




template <typename T, MatrixType Tri, Symmetry Sym>
class DiagonalPreconditionerBase<T, Tri, Sym, CRS>: public Preconditioner<T>
{
public:
	typedef Matrix<MatrixEngine<T, Tri, Sym, CRS> > MatrixType;
	DiagonalPreconditionerBase () { }
	~DiagonalPreconditionerBase () { }
	
	void set(const MatrixBase<T>* Mat);
	void set(const Matrix<MatrixEngine<T, Tri, Sym, CRS> >& Mat);
	
	inline DenseVector<T>& apply(const DenseVector<T>& b, DenseVector<T>& x) { return apply_(b, x); }
	inline DenseVectorView<T>& apply(const DenseVector<T>& b, DenseVectorView<T>& x) { return apply_(b, x); }
	inline DenseVector<T>& apply(const DenseVectorView<T>& b, DenseVector<T>& x) { return apply_(b, x); }
	inline DenseVectorView<T>& apply(const DenseVectorView<T>& b, DenseVectorView<T>& x) { return apply_(b, x); }
	inline DenseVector<T>& apply(const DenseVectorTemp<T>& b, DenseVector<T>& x) { return apply_(b, x); }
	inline DenseVectorView<T>& apply(const DenseVectorTemp<T>& b, DenseVectorView<T>& x) { return apply_(b, x); }
	inline DenseVector<T>& applyTrans(const DenseVector<T>& b, DenseVector<T>& x) { return apply_(b, x); }
	inline DenseVectorView<T>& applyTrans(const DenseVector<T>& b, DenseVectorView<T>& x) { return apply_(b, x); }
	inline DenseVector<T>& applyTrans(const DenseVectorView<T>& b, DenseVector<T>& x) { return apply_(b, x); }
	inline DenseVectorView<T>& applyTrans(const DenseVectorView<T>& b, DenseVectorView<T>& x) { return apply_(b, x); }
	inline DenseVector<T>& applyTrans(const DenseVectorTemp<T>& b, DenseVector<T>& x) { return apply_(b, x); }
	inline DenseVectorView<T>& applyTrans(const DenseVectorTemp<T>& b, DenseVectorView<T>& x) { return apply_(b, x); }
	
private:
	DenseVector<T> invDiag_;
	template <typename V1, typename V2>
	V2& apply_(const V1& b, V2& x);
	
};

template <typename T, MatrixType Tri, Symmetry Sym>
void DiagonalPreconditionerBase<T, Tri, Sym, CRS>::set(const MatrixBase<T>* Mat)
{
	const Matrix<MatrixEngine<T, Tri, Sym, CRS> >* m = dynamic_cast<const Matrix<MatrixEngine<T, Tri, Sym, CRS> >*>(Mat);
	if(m) {
		set(*m);
	} else {
		assert("Wrong Matrix type.  Must be a CRS." && false);
	}
}

template <typename T, MatrixType Tri, Symmetry Sym>
void DiagonalPreconditionerBase<T, Tri, Sym, CRS>::set(const Matrix<MatrixEngine<T, Tri, Sym, CRS> >& m)
{
	//assert(m.engine().final());
	assert(m.rows()==m.cols());
	typename MatrixEngine<T, Tri, Sym, CRS>::size_type rows = m.rows(), row=0;
	typename MatrixEngine<T, Tri, Sym, CRS>::const_val_iterator vals = m.engine().values();
	typename MatrixEngine<T, Tri, Sym, CRS>::const_coord_iterator colInd = m.engine().colIndices(), col_it;
	typename MatrixEngine<T, Tri, Sym, CRS>::const_coord_iterator rowPtr = m.engine().rowPtr(), rowPtrEnd=rowPtr+rows;
	invDiag_.resize(rows);
	while(rowPtr!=rowPtrEnd) {
		col_it=std::find(colInd+*rowPtr, colInd+*(rowPtr+1), row);
		assert(col_it!=colInd+*(rowPtr+1));
		assert(vals[col_it-colInd]!=0.0);
		invDiag_(row) = T(1)/vals[col_it-colInd];
		row++;
		rowPtr++;
	}
}

template <typename T, MatrixType Tri, Symmetry Sym> template <typename V1, typename V2>
V2& DiagonalPreconditionerBase<T, Tri, Sym, CRS>::apply_(const V1& b, V2& x)
{
	if(b.size() != x.size() || x.size() != invDiag_.size())
		std::cout << b.size() << "==" << x.size() << " && " << x.size() << "==" << invDiag_.size() << std::endl;
	assert(b.size() == x.size() && x.size() == invDiag_.size());
	for(int i=0; i<b.size(); ++i)
		x(i)=b(i)*invDiag_(i);
	return x;
}

template <typename M>
class DiagonalPreconditioner: public DiagonalPreconditionerBase<typename M::value_type, M::MT, M::SYM, M::ST> {

};


} /* sparse_lib */ 



#endif /* end of include guard: DIAGONALPRECONDITIONER_H */
