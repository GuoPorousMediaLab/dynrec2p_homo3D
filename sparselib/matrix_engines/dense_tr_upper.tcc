namespace sparse_lib
{

template <typename T>
void TRMV(const MatrixEngine<T, UPPER, ASYMMETRIC, DENSE_FULL>& M, const DenseVectorView<T>& x, DenseVectorView<T>& rhs, bool trans);

template <typename T>
void TRSV(const MatrixEngine<T, UPPER, ASYMMETRIC, DENSE_FULL>& M, DenseVectorView<T>& x, const DenseVectorView<T>& rhs, bool trans);

template <typename T>
class MatrixEngine<T, UPPER, ASYMMETRIC, DENSE_FULL>
{
public:
	typedef T* iterator;
	typedef const T* const_iterator;
	typedef size_t size_type;
	typedef T value_type;
	static const MatrixType MT=UPPER;
	static const Symmetry SYM=ASYMMETRIC;
	static const MatrixStorage ST=DENSE_FULL;
	
	MatrixEngine (): nRows_(0), nCols_(0), data_(0), limit_(0) { }
	MatrixEngine (const size_type rows, const size_type cols, const value_type v=value_type())
		: nRows_(rows), nCols_(cols) { assert(nRows_==nCols_); DenseMatrixAllocate(*this, nRows_*nCols_, v); }
	~MatrixEngine () { DenseMatrixDeallocate(*this, nRows_*nCols_); }
	
	inline size_type rows() const { return nRows_; }
	inline size_type cols() const { return nCols_; }
	
	value_type& operator()(const size_type r, const size_type c) { assert(c>=r && r<nRows_ && c<nCols_); return data_[c+r*nCols_]; }
	const value_type& operator()(const size_type r, const size_type c) const { assert(c>=r && r<nRows_ && c<nCols_); return data_[c+r*nCols_]; }

	void set(const size_type r, const size_type c, const value_type val) { assert(c>=r && r<nRows_ && c<nCols_); data_[c+r*nCols_] = val; }
	value_type& get(const size_type r, const size_type c) { assert(c>=r && r<nRows_ && c<nCols_); return data_[c+r*nCols_]; }
	const value_type& get(const size_type r, const size_type c) const { assert(c>=r && r<nRows_ && c<nCols_); return data_[c+r*nCols_]; }
	
	iterator data() { return data_; }
	const_iterator data() const { return data_; }
	
	template <typename V1, typename V2>
	void mv(const V1& x, V2& y, bool trans, T alpha, T beta) const;
	
	template <typename V1, typename V2>
	void solve(V1& x, const V2& rhs, bool trans) const;
	
	void print(std::ostream &o) const;

	void finalize() { }

	friend void DenseMatrixAllocate<>(MatrixEngine&, const size_type, const value_type);
	friend void DenseMatrixDeallocate<>(MatrixEngine&, const size_type);

private:
	iterator data_;
	iterator limit_;
	size_type nRows_, nCols_;
};

template <typename T>
void MatrixEngine<T, UPPER, ASYMMETRIC, DENSE_FULL>::print(std::ostream& o) const
{
	o << "[";
	for(int i=0; i<nRows_; i++) {
		if(i>0)
			o << " [";
		else
			o << "[";
		for(int j=0; j<nCols_-1; j++) 
			o << data_[i*nCols_+j] << " ";
		o << data_[i*nCols_+nCols_-1];
		if (i<nRows_-1)
			o << "]\n";
		else
			o << "]]";
	}
}

template <typename T> template <typename V1, typename V2>
inline void MatrixEngine<T, UPPER, ASYMMETRIC, DENSE_FULL>::mv(const V1& x, V2& y, bool trans, T alpha, T beta) const
{
	const DenseVectorView<T> xv(x);
	DenseVectorView<T> yv(y);
	TRMV(*this, xv, yv, trans);
}

template <typename T> template <typename V1, typename V2>
inline void MatrixEngine<T, UPPER, ASYMMETRIC, DENSE_FULL>::solve(V1& x, const V2& rhs, bool trans) const
{
	DenseVectorView<T> xv(x);
	const DenseVectorView<T> rhsv(rhs);
	TRSV(*this, xv, rhsv, trans);
}



template <typename T>
void TRMV(const MatrixEngine<T, UPPER, ASYMMETRIC, DENSE_FULL>& M, const DenseVectorView<T>& x, DenseVectorView<T>& rhs, bool trans)
{
	const T* M_data=M.data();
	const T* x_data=x.data();
	T* rhs_data=rhs.data();
	int rows=M.rows();
	int cols=M.cols();
	
	if(!trans) {
		for(int i=0; i<rows; i++) {
			for(int j=0; j<cols; j++) {
				rhs_data[i]+=M_data[i*cols+j]*x_data[j];
			}
		}
	} else {
		for(int i=0; i<rows; i++) {
			for(int j=0; j<cols; j++) {
				rhs_data[j]+=M_data[i*cols+j]*x_data[i];
			}
		}
	}

}



#ifdef CBLAS
template <>
void TRMV(const MatrixEngine<double, UPPER, ASYMMETRIC, DENSE_FULL>& M, const DenseVectorView<double>& x, DenseVectorView<double>& rhs, bool trans)
{
	CBLAS_TRANSPOSE trans_=CblasNoTrans;
	if(trans)
		trans_=CblasTrans;
	rhs=x;
	cblas_dtrmv(CblasRowMajor, CblasUpper, trans_, CblasNonUnit, M.rows(),
	            M.data(), M.rows(), rhs.data(), 1);
}

template <>
void TRMV(const MatrixEngine<float, UPPER, ASYMMETRIC, DENSE_FULL>& M, const DenseVectorView<float>& x, DenseVectorView<float>& rhs, bool trans)
{
	CBLAS_TRANSPOSE trans_=CblasNoTrans;
	if(trans)
		trans_=CblasTrans;
	rhs=x;
	cblas_strmv(CblasRowMajor, CblasUpper, trans_, CblasNonUnit, M.rows(),
	            M.data(), M.rows(), rhs.data(), 1);
}
#endif

template <typename T>
void TRSV(const MatrixEngine<T, UPPER, ASYMMETRIC, DENSE_FULL>& M, DenseVectorView<T>& x, const DenseVectorView<T>& rhs, bool trans)
{
	assert("General Solve routine for Triangular matrices not yet implemented.  Must use type <double> or <float>." && false);
}



#ifdef CBLAS
template <>
void TRSV(const MatrixEngine<double, UPPER, ASYMMETRIC, DENSE_FULL>& M, DenseVectorView<double>& x, const DenseVectorView<double>& rhs, bool trans)
{
	CBLAS_TRANSPOSE trans_=CblasNoTrans;
	if(trans)
		trans_=CblasTrans;
	x=rhs;
	cblas_dtrsv(CblasRowMajor, CblasUpper, trans_, CblasNonUnit, M.rows(),
	            M.data(), M.rows(), x.data(), 1);
}

template <>
void TRSV(const MatrixEngine<float, UPPER, ASYMMETRIC, DENSE_FULL>& M, DenseVectorView<float>& x, const DenseVectorView<float>& rhs, bool trans)
{
	CBLAS_TRANSPOSE trans_=CblasNoTrans;
	if(trans)
		trans_=CblasTrans;
	x=rhs;
	cblas_strsv(CblasRowMajor, CblasUpper, trans_, CblasNonUnit, M.rows(),
	            M.data(), M.rows(), x.data(), 1);
}
#endif


} /* sparse_lib */ 
