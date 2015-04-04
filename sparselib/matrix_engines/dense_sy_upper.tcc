namespace sparse_lib
{

template <typename T>
void SYMV(const MatrixEngine<T, UPPER, SYMMETRIC, DENSE_FULL>& M, const DenseVectorView<T>& x, DenseVectorView<T>& rhs, T alpha, T beta);

template <typename T>
class MatrixEngine<T, UPPER, SYMMETRIC, DENSE_FULL>
{
public:
	typedef T* iterator;
	typedef const T* const_iterator;
	typedef size_t size_type;
	typedef T value_type;
	static const MatrixType MT=UPPER;
	static const Symmetry SYM=SYMMETRIC;
	static const MatrixStorage ST=DENSE_FULL;
	
	MatrixEngine (): nRows_(0), nCols_(0), data_(0), limit_(0) { }
	MatrixEngine (const size_type rows, const size_type cols, const value_type v=value_type())
		: nRows_(rows), nCols_(cols) { assert(nRows_==nCols_); DenseMatrixAllocate(*this, nRows_*nCols_, v); }
	~MatrixEngine () { DenseMatrixDeallocate(*this, nRows_*nCols_); }
	
	inline size_type rows() const { return nRows_; }
	inline size_type cols() const { return nCols_; }
	
	value_type& operator()(const size_type r, const size_type c) { assert(r<nRows_ && c<nCols_); return data_[std::max(c,r)+std::min(c,r)*nCols_]; }
	const value_type& operator()(const size_type r, const size_type c) const { assert(r<nRows_ && c<nCols_); return data_[std::max(c,r)+std::min(c,r)*nCols_]; }

	void set(const size_type r, const size_type c, const value_type val) { assert(r<nRows_ && c<nCols_); data_[std::max(c,r)+std::min(c,r)*nCols_] = val; }
	value_type& get(const size_type r, const size_type c) { assert(r<nRows_ && c<nCols_); return data_[std::max(c,r)+std::min(c,r)*nCols_]; }
	const value_type& get(const size_type r, const size_type c) const { assert(r<nRows_ && c<nCols_); return data_[std::max(c,r)+std::min(c,r)*nCols_]; }
	
	iterator data() { return data_; }
	const_iterator data() const { return data_; }
	
	template <typename V1, typename V2>
	void mv(const V1& x, V2& y, bool trans, T alpha, T beta) const;
	
	template <typename V1, typename V2>
	void solve(V1& x, const V2& rhs, bool trans) const { assert("Solve is not available for symmetric UT matrices" && false); }
	
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
void MatrixEngine<T, UPPER, SYMMETRIC, DENSE_FULL>::print(std::ostream& o) const
{
	o << "[";
	for(int i=0; i<nRows_; i++) {
		if(i>0)
			o << " [";
		else
			o << "[";
		for(int j=0; j<nCols_-1; j++) {
			if(i<=j)
				o << data_[i*nCols_+j] << " ";
			else
				o << data_[j*nCols_+i] << " ";
		}
		o << data_[i*nCols_+nCols_-1];
		if (i<nRows_-1)
			o << "]\n";
		else
			o << "]]";
	}
}

template <typename T> template <typename V1, typename V2>
inline void MatrixEngine<T, UPPER, SYMMETRIC, DENSE_FULL>::mv(const V1& x, V2& y, bool trans, T alpha, T beta) const
{
	const DenseVectorView<T> xv(x);
	DenseVectorView<T> yv(y);
	SYMV(*this, xv, yv, alpha, beta);
}

template <typename T>
void SYMV(const MatrixEngine<T, UPPER, SYMMETRIC, DENSE_FULL>& M, const DenseVectorView<T>& x, DenseVectorView<T>& rhs, T alpha, T beta)
{
	const T* M_data=M.data();
	const T* x_data=x.data();
	T* rhs_data=rhs.data();
	int rows=M.rows();
	int cols=M.cols();
	
	if(beta!=1.0)
		rhs*=beta;

	assert("Upper symv routine is not yet implemented.  Must use type <double> or <float> and cblas." && false);
}

#ifdef CBLAS
template <>
void SYMV(const MatrixEngine<double, UPPER, SYMMETRIC, DENSE_FULL>& M, const DenseVectorView<double>& x, DenseVectorView<double>& rhs, double alpha, double beta)
{
	cblas_dsymv(CblasRowMajor, CblasUpper, M.rows(), alpha, M.data(), M.rows(), x.data(), 1, beta,
	            rhs.data(), 1);
}

template <>
void SYMV(const MatrixEngine<float, UPPER, SYMMETRIC, DENSE_FULL>& M, const DenseVectorView<float>& x, DenseVectorView<float>& rhs, float alpha, float beta)
{
	cblas_ssymv(CblasRowMajor, CblasUpper, M.rows(), alpha, M.data(), M.rows(), x.data(), 1, beta,
	            rhs.data(), 1);
}
#endif


} /* sparse_lib */ 
