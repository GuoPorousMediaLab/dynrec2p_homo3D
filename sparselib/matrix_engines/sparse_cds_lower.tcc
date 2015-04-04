namespace sparse_lib
{

template <typename T, MatrixType Tri, Symmetry Sym, MatrixStorage MatStorage>
class ICCPreconditionerBase;

template <typename T>
class MatrixEngine<T, LOWER, ASYMMETRIC, CDS>
{
public:
	friend class ICCPreconditionerBase<T, UPPER, SYMMETRIC, CDS>;
	typedef T* val_iterator;
	typedef const T* const_val_iterator;
	typedef int size_type;
	typedef int* coord_iterator;
	typedef const int* const_coord_iterator;
	typedef T value_type;
	static const MatrixType MT=LOWER;
	static const Symmetry SYM=ASYMMETRIC;
	static const MatrixStorage ST=CDS;
	
	MatrixEngine (): nRows_(0), nCols_(0), zero_(T()) { allocate_(); }
	MatrixEngine (const size_type rows, const size_type cols, const size_type bw=0);
	~MatrixEngine () { release_(); }

	inline size_type rows() const { return nRows_; }
	inline size_type cols() const { return nCols_; }
	inline size_type nDiags() const { return diags_end_-diags_; }

	value_type& operator()(const size_type r, const size_type c);
	const value_type& operator()(const size_type r, const size_type c) const;

	void set(const size_type r, const size_type c, const value_type val);


	val_iterator values() { return values_; }
	const_val_iterator values() const { return values_; }

	coord_iterator diags() { return diags_; }
	const_coord_iterator diags() const { return diags_; }

	template <typename V1, typename V2>
	void mv(const V1& x, V2& y, bool trans, T alpha, T beta) const;

	template <typename V1, typename V2>
	void solve(V1& x, const V2& b, bool trans) const;
	
	template <typename V1>
	void solve(V1& x, bool trans) const;

	void print(std::ostream &o) const;
	void print_matlab(std::ostream &o, std::string varName, int r0=0, int c0=0) const { }

	void finalize() { };

private:
	val_iterator values_;
	val_iterator values_end_;
	val_iterator values_limit_;
	coord_iterator diags_;
	coord_iterator diags_end_;
	coord_iterator diags_limit_;
	std::allocator<T> values_alloc_;
	std::allocator<int> diags_alloc_;
	size_type nRows_, nCols_;
	const T zero_;
	
	void allocate_();
	void allocate_(size_type size);
	void release_();
	T& insert_diag_(int diag, size_type r, size_type c);
	void insert_diag_(int diag);
};

template <typename T>
MatrixEngine<T, LOWER, ASYMMETRIC, CDS>::MatrixEngine (const size_type rows, const size_type cols, const size_type bw)
	:nRows_(rows), nCols_(cols), zero_(T())
{
	assert(nRows_ == nCols_);
	if(bw>0)
		allocate_(bw);
	else
		allocate_();
}

template <typename T>
void MatrixEngine<T, LOWER, ASYMMETRIC, CDS>::allocate_()
{
	values_ = values_end_ = values_limit_ = 0;
	diags_ = diags_end_ = diags_limit_ = 0;
}

template <typename T>
void MatrixEngine<T, LOWER, ASYMMETRIC, CDS>::allocate_(size_type bw)
{
	values_ = new value_type[bw*nRows_];
	values_limit_ = values_ + bw*nRows_;
	values_end_ = values_;
	
	diags_ = new int[bw];
	diags_limit_ = diags_ + bw;
	diags_end_ = diags_;
}

template <typename T>
T& MatrixEngine<T, LOWER, ASYMMETRIC, CDS>::operator()(const size_type r, const size_type c)
{
	assert(r<nRows_ && c<nCols_);
	assert(r>=c);
	int diag = c-r;
	coord_iterator val_it = std::find(diags_, diags_end_, diag);
	if(val_it != diags_end_)
		return values_[(val_it-diags_)*nRows_ + r];
	else
		return insert_diag_(diag, r, c);
}

template <typename T>
const T& MatrixEngine<T, LOWER, ASYMMETRIC, CDS>::operator()(const size_type r, const size_type c) const
{
	assert(r<nRows_ && c<nCols_);
	assert(r>=c);
	int diag = c-r;
	coord_iterator val_it = std::find(diags_, diags_end_, diag);
	if(val_it != diags_end_)
		return values_[(val_it-diags_)*nRows_ + r];
	else
		return zero_;
}

template <typename T>
void MatrixEngine<T, LOWER, ASYMMETRIC, CDS>::set(const size_type r, const size_type c, const value_type val)
{
	assert(r<nRows_ && c<nCols_);
	assert(r>=c);
	int diag = c-r;
	coord_iterator val_it = std::find(diags_, diags_end_, diag);
	if(val_it != diags_end_)
		values_[(val_it-diags_)*nRows_ + r] = val;
	else
		insert_diag_(diag, r, c) = val;
}

template <typename T>
T& MatrixEngine<T, LOWER, ASYMMETRIC, CDS>::insert_diag_(int diag, size_type r, size_type c)
{
	coord_iterator it = diags_end_;
	while(it!=diags_) {
		if(*(it-1) < diag)
			break;
		--it;
	}
	
	
	if(diags_end_ != diags_limit_) {
		diags_end_++;
		std::copy_backward(it, diags_end_-1, diags_end_);
		*it=diag;
		int insert_pos = it-diags_;
		values_end_+=nRows_;
		std::copy_backward(values_+insert_pos*nRows_, values_end_-nRows_, values_end_);
		std::fill(values_+insert_pos*nRows_, values_+insert_pos*nRows_+nRows_, T());
		return values_[insert_pos*nRows_+r];
	} else {
		std::cerr << "Reallocation not completed yet!" << std::endl;
		return values_[0];
	}
}

template <typename T>
void MatrixEngine<T, LOWER, ASYMMETRIC, CDS>::insert_diag_(int diag)
{
	coord_iterator it = diags_end_;
	while(it!=diags_) {
		if(*(it-1) < diag)
			break;
		--it;
	}
	
	
	if(diags_end_ != diags_limit_) {
		diags_end_++;
		std::copy_backward(it, diags_end_-1, diags_end_);
		*it=diag;
		int insert_pos = it-diags_;
		values_end_+=nRows_;
		std::copy_backward(values_+insert_pos*nRows_, values_end_-nRows_, values_end_);
		std::fill(values_+insert_pos*nRows_, values_+insert_pos*nRows_+nRows_, T());
	} else {
		std::cerr << "Reallocation not completed yet!" << std::endl;
	}
}

template <typename T>
void MatrixEngine<T, LOWER, ASYMMETRIC, CDS>::release_()
{
	//std::cout << "releasing @ " << this << std::endl;
	//std::cout << values_ << std::endl;
	if(values_)
		delete [] values_;
	if(diags_)
		delete [] diags_;
	values_ = values_end_ = values_limit_ = 0;
	diags_ = diags_end_ = diags_limit_ = 0;
}

template <typename T>
void MatrixEngine<T, LOWER, ASYMMETRIC, CDS>::print(std::ostream &o) const
{
	o << "[" << nRows_ << ", " << nCols_ << "]: \n";
	int size = diags_end_ - diags_;
	if(size>0) {
		o << " values= \n";
		for(int i=0; i<nRows_; ++i) {
			o << "  [ ";
			for(int j=0; j<size; j++) {
				o << std::setw(8) << values_[j*nRows_ + i] << " ";
			}
			o << "]\n";
		}
		o << " diags= \n";
		for(int j=0; j<size; j++)
			o << diags_[j] << " ";
	} else {
		o << " values= \n  [ EMPTY ]\n";
		o << " diags= \n  [ EMPTY ]\n";
	}
	o << "\n";
}

template <typename T> template <typename V1, typename V2>
inline void MatrixEngine<T, LOWER, ASYMMETRIC, CDS>::mv(const V1& x, V2& y, bool trans, T alpha, T beta) const
{
	assert(trans ? (x.size()==nRows_ && y.size()==nRows_) : (x.size()==nCols_ && y.size()==nRows_) );
	int nDiags = diags_end_-diags_;
	cds_mv(trans, nRows_, nCols_, alpha, values_, nDiags, diags_, x.data(), beta, y.data());
}

template <typename T> template <typename V1, typename V2>
void MatrixEngine<T, LOWER, ASYMMETRIC, CDS>::solve(V1& x, const V2& rhs, bool trans) const
{
	assert (x.size()==nRows_ && rhs.size()==nRows_);
	assert (nRows_ == nCols_);
	int nDiags = diags_end_-diags_;
	x=rhs;
	cds_sv(trans, false, nRows_, values_, nDiags, diags_, x.data());
}

template <typename T> template <typename V1>
void MatrixEngine<T, LOWER, ASYMMETRIC, CDS>::solve(V1& x, bool trans) const
{
	assert (x.size()==nRows_);
	assert (nRows_ == nCols_);
	int nDiags = diags_end_-diags_;
	cds_sv(trans, false, nRows_, values_, nDiags, diags_, x.data());
}

} /* sparse_lib */ 
