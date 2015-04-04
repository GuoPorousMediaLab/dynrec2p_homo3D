#include <vector>
#include <utility> // for std::pair
#include <memory>
#include <algorithm>

namespace sparse_lib
{

template <typename T>
class SchwarzLinearSystemBase
{
public:
	virtual bool solve(DenseVectorView<T> &x, const DenseVectorView<T> &b) = 0;
	virtual MatrixBase<T>* getMatrix() = 0;
	virtual const MatrixBase<T>* getMatrix() const = 0;
	virtual T& operator()(const int r, const int c) = 0;
	virtual ~SchwarzLinearSystemBase() { }
private:
	
};

template <typename M, typename P, SolverType ST=NO_SOLVER>
class SchwarzLinearSystem: public SchwarzLinearSystemBase<typename M::value_type>
{
public:
	typedef typename M::value_type T;
	SchwarzLinearSystem() { std::cerr << "spasre_lib error: Schwarz linear system not defined for the solver type selected." << std::endl; assert(false); }
	
	MatrixBase<T>* getMatrix() { return &theMatrix_; }
	const MatrixBase<T>* getMatrix() const { return &theMatrix_; }
	
	T& operator()(const int r, const int c) { return theMatrix_(r,c); }
	bool solve(DenseVectorView<T> &x, const DenseVectorView<T> &b) { return false; }
private:
	M theMatrix_;
	P Preconditioner_;
	Iterator<T> it_;
};

template <typename M, typename P>
class SchwarzLinearSystem<M, P, GMRES>: public SchwarzLinearSystemBase<typename M::value_type>
{
public:
	typedef typename M::value_type T;
	SchwarzLinearSystem(int nRows, int nCols, int bandwidth, int restart, T tol, int maxIt, iter_output it_out=ITER_NONE,
						std::string iter_id="", std::ostream &it_stream=std::cout): theMatrix_(nRows,nCols,bandwidth), restart_(restart),
						it_(tol, maxIt, 0.0, it_out, it_stream) { }
	
	MatrixBase<T>* getMatrix() { return &theMatrix_; }
	const MatrixBase<T>* getMatrix() const { return &theMatrix_; }
	
	T& operator()(const int r, const int c) { return theMatrix_(r,c); }
	bool solve(DenseVectorView<T> &x, const DenseVectorView<T> &b) { P Prec; return GMRES_Solver(theMatrix_, x, b, Prec, it_, restart_, true); }
private:
	int restart_;
	M theMatrix_;
	Iterator<T> it_;
};

template <typename M, typename P>
class SchwarzLinearSystem<M, P, CG>: public SchwarzLinearSystemBase<typename M::value_type>
{
public:
	typedef typename M::value_type T;
	SchwarzLinearSystem(int nRows, int nCols, int bandwidth, T tol, int maxIt, iter_output it_out=ITER_NONE,
						std::string iter_id="", std::ostream &it_stream=std::cout): theMatrix_(nRows,nCols,bandwidth),
						it_(tol, maxIt, tol, it_out, it_stream), preconditioned_(false) { }
	
	MatrixBase<T>* getMatrix() { return &theMatrix_; }
	const MatrixBase<T>* getMatrix() const { return &theMatrix_; }
	
	T& operator()(const int r, const int c) { return theMatrix_(r,c); }
	bool solve(DenseVectorView<T> &x, const DenseVectorView<T> &b) { if(!preconditioned_) {Preconditioner_.set(&theMatrix_); preconditioned_=true; } return CG_Solver(theMatrix_, x, b, Preconditioner_, it_, false); }
private:
	M theMatrix_;
	P Preconditioner_;
	bool preconditioned_;
	Iterator<T> it_;
};


template <typename T>
class MatrixEngine<T, GENERAL, ASYMMETRIC, SCHWARZ_BLOCK>
{
public:
	typedef T* val_iterator;
	typedef const T* const_val_iterator;
	typedef int size_type;
	typedef size_type* coord_iterator;
	typedef const size_type* const_coord_iterator;
	typedef T value_type;
	typedef Matrix<MatrixEngine<T> > bgType;
	static const MatrixType MT=GENERAL;
	static const Symmetry SYM=ASYMMETRIC;
	static const MatrixStorage ST=SCHWARZ_BLOCK;
	
	MatrixEngine (): background_(1,1), mainDiagBlocks_(0), block_ptr_(0), nRows_(0), nCols_(0) { block_ptr_.push_back(0); }
	~MatrixEngine () { release_(); }
	
	inline size_type rows() const { return nRows_; }
	inline size_type cols() const { return nCols_; }
	
	inline size_type num_blocks() const { return mainDiagBlocks_.size(); }
	
	value_type& operator()(const size_type r, const size_type c);
	const value_type& operator()(const size_type r, const size_type c) const;
	void set(const size_type r, const size_type c, const value_type val);
	
	template <typename M, typename P, SolverType SolveT>
	void push_back(int nRows, int nCols, int bandwidth, T tol, int maxIt, iter_output it_out=ITER_NONE,
		std::string iter_id="", std::ostream &it_stream=std::cout)
		{
			assert(nRows == nCols);
			assert(SolveT!=GMRES);
			block_ptr_.push_back(block_ptr_[block_ptr_.size()-1]+nRows);
			nRows_+=nRows;
			nCols_+=nCols;
			mainDiagBlocks_.push_back(new SchwarzLinearSystem<M, P, SolveT>(nRows, nCols, bandwidth, tol, maxIt, it_out, iter_id, it_stream));
			background_.engine().resize(nRows_, nCols_);
		}
	
	template <typename M, typename P, SolverType SolveT>
	void push_back(int nRows, int nCols, int bandwidth, int restart, T tol, int maxIt, iter_output it_out=ITER_NONE,
		std::string iter_id="", std::ostream &it_stream=std::cout)
		{
			assert(nRows == nCols);
			assert(SolveT==GMRES);
			block_ptr_.push_back(block_ptr_[block_ptr_.size()-1]+nRows);
			nRows_+=nRows;
			nCols_+=nCols;
			mainDiagBlocks_.push_back(new SchwarzLinearSystem<M, P, SolveT>(nRows, nCols, bandwidth, restart, tol, maxIt, it_out, iter_id, it_stream));
			background_.engine().resize(nRows_, nCols_);
		}
	
	template <typename M>
	void grow_block(size_type blockId, size_type g)
	{
		M *block = dynamic_cast<M*>(get_block(blockId));
		nRows_+=g;
		nCols_+=g;
		size_type curSize = block->rows();
		assert(curSize >= 0);
		block->engine().resize(curSize+g, curSize+g);
		background_.engine().resize(nRows_, nCols_);
		for(int i=blockId+1; i<block_ptr_.size(); ++i)
			block_ptr_[i]+=g;
	}
	
	MatrixBase<T>* get_block(size_type loc);
	bgType& get_background() { return background_; }
	bgType const& get_background() const { return background_; }
	
	std::vector<SchwarzLinearSystemBase<value_type>*> const& get_blocks() const { return mainDiagBlocks_; }
	std::vector<size_type> const& get_block_ptr() const { return block_ptr_; }
	
		
	template <typename V1, typename V2>
	void mv(const V1& x, V2& y, bool trans, T alpha, T beta) const;
	
	template <typename V1, typename V2>
	void solve(V1& x, const V2& rhs, bool trans) const { assert("Solve is not available.  Only for UT and LT matrices" && false); }
	
	void print_matlab(std::ostream &o, std::string varName, int r0=0, int c0=0) const;
	void finalize();
	
	void print(std::ostream &o) const { std::cout << "SchwarzMat: [" << nRows_ << ", " << nCols_ << "]" << std::endl;
	
		for(int i=0; i< mainDiagBlocks_.size(); ++i) {
			std::cout << "   Block " << i+1 << ": [" << mainDiagBlocks_[i]->getMatrix()->rows() << ", " << mainDiagBlocks_[i]->getMatrix()->cols() << "]" << std::endl;
			std::cout << "   block_ptr[" << i << "]: " << block_ptr_[i] << std::endl;
		}
		std::cout << "   block_ptr.back(): " << block_ptr_.back() << std::endl;
	}
	
private:
	bgType background_;
	std::vector<SchwarzLinearSystemBase<value_type>*> mainDiagBlocks_;
	std::vector<size_type> block_ptr_;
	size_type nRows_, nCols_;
	void release_();
};

template <typename T>
MatrixBase<T>* MatrixEngine<T, GENERAL, ASYMMETRIC, SCHWARZ_BLOCK>::get_block(size_type loc)
{
	assert(loc<nRows_);
	std::vector<size_type>::iterator it=block_ptr_.begin();
	while(loc >= *it && it!=block_ptr_.end())
		it++;
	it--;
	return mainDiagBlocks_[it-block_ptr_.begin()]->getMatrix();
}

template <typename T>
T& MatrixEngine<T, GENERAL, ASYMMETRIC, SCHWARZ_BLOCK>::operator()(const size_type r, const size_type c)
{
	assert(r<nRows_ && c<nCols_);
	std::vector<size_type>::iterator it=block_ptr_.begin();
	while(r >= *it && it != block_ptr_.end())
		it++;
	it--;
		
	if(r < *(it+1) && r>= *it && c < *(it+1) && c>= *it)
		return mainDiagBlocks_[it-block_ptr_.begin()]->getMatrix()->operator()(r-*it, c-*it);
	else {
		return background_(r, c);
	}
}

template <typename T>
const T& MatrixEngine<T, GENERAL, ASYMMETRIC, SCHWARZ_BLOCK>::operator()(const size_type r, const size_type c) const
{
	assert(r<nRows_ && c<nCols_);
	std::vector<size_type>::const_iterator it=block_ptr_.begin();
	while(r >= *it)
		it++;
	it--;
	if(r < *(it+1) && r>= *it && c < *(it+1) && c>= *it)
		return mainDiagBlocks_[it-block_ptr_.begin()]->getMatrix()->operator()(r-*it, c-*it);
	else
		return background_(r, c);
}

template <typename T>
void MatrixEngine<T, GENERAL, ASYMMETRIC, SCHWARZ_BLOCK>::set(const size_type r, const size_type c, const T val)
{
	assert(r<nRows_ && c<nCols_);
	std::vector<size_type>::iterator it=block_ptr_.begin();
	while(r >= *it && it != block_ptr_.end())
		it++;
	it--;
	if(r < *(it+1) && r>= *it && c < *(it+1) && c>= *it)
		mainDiagBlocks_[it-block_ptr_.begin()]->getMatrix()->operator()(r-*it, c-*it) = val;
	else
		background_(r, c) = val;
}

template <typename T>
void MatrixEngine<T, GENERAL, ASYMMETRIC, SCHWARZ_BLOCK>::finalize()
{
	background_.finalize();
	for(int i=0; i<mainDiagBlocks_.size(); ++i)
		mainDiagBlocks_[i]->getMatrix()->finalize();
}

template <typename T> template <typename V1, typename V2>
void MatrixEngine<T, GENERAL, ASYMMETRIC, SCHWARZ_BLOCK>::mv(const V1& x, V2& y, bool trans, T alpha, T beta) const
{
	assert(trans ? (x.size()==nRows_ && y.size()==nRows_) : (x.size()==nCols_ && y.size()==nRows_) );
	background_.mv(x, y, trans, alpha, beta);
	int i=0, size=mainDiagBlocks_.size();
	// Cant use iterators here due to OpenMP parallelization.
	//#pragma omp parallel for default(none) private(i) shared(size, x, y, trans, alpha)
	for(i=0; i<size; ++i) {
		int start = block_ptr_[i], end = block_ptr_[i+1];
		mainDiagBlocks_[i]->getMatrix()->mv(x.range(start, end-1), y.range(start, end-1), trans, alpha, 1.0);
	}
}

template <typename T>
void MatrixEngine<T, GENERAL, ASYMMETRIC, SCHWARZ_BLOCK>::print_matlab(std::ostream &o, std::string varName, int r0, int c0) const
{
	typename std::vector<SchwarzLinearSystemBase<value_type>*>::const_iterator m_it = mainDiagBlocks_.begin();
	std::vector<size_type>::const_iterator ind_it = block_ptr_.begin();
	for(; m_it!=mainDiagBlocks_.end(); ++m_it, ++ind_it ) {
		(*m_it)->getMatrix()->print_matlab(o, varName, r0+*ind_it, c0+*ind_it);
	}
}

template <typename T>
void MatrixEngine<T, GENERAL, ASYMMETRIC, SCHWARZ_BLOCK>::release_()
{
	typename std::vector<SchwarzLinearSystemBase<value_type>*>::iterator m_it;
	for(m_it = mainDiagBlocks_.begin(); m_it!=mainDiagBlocks_.end(); ++m_it) {
		if(*m_it) {
			delete *m_it;
		}
	}
}


} /* sparse_lib */ 
