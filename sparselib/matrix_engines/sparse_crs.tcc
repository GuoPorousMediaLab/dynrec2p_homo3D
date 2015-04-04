#include <vector>
#include <utility> // for std::pair
#include <memory>
#include <algorithm>
#include <boost/timer.hpp>
#include <sparselib/sparse_blas.h>

namespace sparse_lib
{

template <typename T>
class MatrixEngine<T, GENERAL, ASYMMETRIC, CRS>;

template <typename T>
struct coord_tuple
{
	typedef typename MatrixEngine<T, GENERAL, ASYMMETRIC, CRS>::size_type size_type;
	coord_tuple(): row(0), col(0), value(T()) { }
	coord_tuple(const size_type r, const size_type c, T v=T()): row(r), col(c), value(v) { }
	size_type row, col;
	T value;
};

struct coord_tuple_cmp
{
    // return true if a < b
    //   <=>  a(1)<b(1)  or  a(1)==b(1) and a(2)<b(2)
    template <typename T>
	bool operator()(const coord_tuple<T> &a, const coord_tuple<T> &b) const;
};

template <typename T>
bool coord_tuple_cmp::operator()(const coord_tuple<T> &a, const coord_tuple<T> &b) const
{
    if (a.row<b.row) {
        return true;
    }
    if ((a.row==b.row) && (a.col<b.col)) {
        return true;
    }
    return false;
}

template <typename T>
class MatrixEngine<T, GENERAL, ASYMMETRIC, CRS>
{
public:
	typedef T* val_iterator;
	typedef const T* const_val_iterator;
	typedef int size_type;
	typedef size_type* coord_iterator;
	typedef const size_type* const_coord_iterator;
	typedef T value_type;
	static const MatrixType MT=GENERAL;
	static const Symmetry SYM=ASYMMETRIC;
	static const MatrixStorage ST=CRS;
	
	MatrixEngine (): nRows_(0), nCols_(0), finalized_(false), sorted_(false), values_(0), values_limit_(0), col_indices_(0),
	col_indices_limit_(0), row_ptrs_(0), row_ptrs_limit_(0) { }
	MatrixEngine (const size_type rows, const size_type cols, const size_type bw=0);
	~MatrixEngine () { release_(); }
	
	inline size_type rows() const { return nRows_; }
	inline size_type cols() const { return nCols_; }
	
	value_type& operator()(const size_type r, const size_type c);
	const value_type& operator()(const size_type r, const size_type c) const;

	val_iterator values() { return values_; }
	const_val_iterator values() const { return values_; }
	
	coord_iterator colIndices() { return col_indices_; }
	const_coord_iterator colIndices() const { return col_indices_; }
	
	coord_iterator rowPtr() { return row_ptrs_; }
	const_coord_iterator rowPtr() const { return row_ptrs_; }

	void set(const size_type r, const size_type c, const value_type val);
	
	template <typename V1, typename V2>
	void mv(const V1& x, V2& y, bool trans, T alpha, T beta) const;
	
	/*template <typename V1, typename V2>
	void mv(const V1& x, V2& y, size_type r_start, size_type r_end, bool trans = false, T alpha=1.0, T beta=0.0) const;*/
	
	template <typename V1, typename V2>
	void solve(V1& x, const V2& rhs, bool trans) const { assert("Solve is not available.  Only for UT and LT matrices" && false); }
	
	void print(std::ostream &o) const;
	void print_matlab(std::ostream &o, std::string varName, int r0=0, int c0=0) const { }
	
	void finalize();
	
	bool final() { return finalized_; }
	
	void resize(size_type nr, size_type nc);
	
	int nnz() const { return coord_vals_.size(); }
	
	template <typename M, typename P, SolverType SolveT>
	void push_back(int nRows, int nCols, int bandwidth, int restart, T tol, int maxIt, iter_output it_out=ITER_NONE,
		std::string iter_id="", std::ostream &it_stream=std::cout)
		{
			abort();
		}

private:
	std::vector<coord_tuple<T> > coord_vals_;
	val_iterator values_;
	val_iterator values_limit_;
	coord_iterator col_indices_;
	coord_iterator col_indices_limit_;
	coord_iterator row_ptrs_;
	coord_iterator row_ptrs_limit_;
	bool finalized_;
	bool sorted_;
	size_type nRows_, nCols_;
	coord_tuple_cmp less_;
	std::allocator<T> value_alloc_;
	std::allocator<size_type> index_alloc_;
	
	void sort_();
	void release_();
	void allocate_(const size_type size);
	void assemble_crs_();
};

template <typename T>
MatrixEngine<T, GENERAL, ASYMMETRIC, CRS>::MatrixEngine (const size_type rows, const size_type cols, const size_type bw)
	:nRows_(rows), nCols_(cols), finalized_(false), sorted_(false), values_(0), values_limit_(0), col_indices_(0),
	col_indices_limit_(0), row_ptrs_(0), row_ptrs_limit_(0)
{
	if(bw>0)
		coord_vals_.reserve(bw*nRows_);
}

template <typename T>
T& MatrixEngine<T, GENERAL, ASYMMETRIC, CRS>::operator()(const size_type r, const size_type c)
{
	assert(r<nRows_ && c<nCols_ && !finalized_);
	coord_vals_.push_back(coord_tuple<T>(r, c));
	size_type i = coord_vals_.size()-1;
	return coord_vals_[i].value;
}

template <typename T>
const T& MatrixEngine<T, GENERAL, ASYMMETRIC, CRS>::operator()(const size_type r, const size_type c) const
{
	// TODO: Fix this- it should be a search because we are doing a const ().
	// Set up to throw an assert if called right now.
	assert(r<nRows_ && c<nCols_ && !finalized_ && false);
	std::cerr << "WARNING- const operator() was called on CRS engine.  This is not yet implemented. Danger, Danger!";
	size_type i = coord_vals_.size()-1;
	return coord_vals_[i].value;
}

template <typename T>
inline void MatrixEngine<T, GENERAL, ASYMMETRIC, CRS>::set(const size_type r, const size_type c, const T val)
{
	assert(r<nRows_ && c<nCols_ && !finalized_);
	coord_vals_.push_back(coord_tuple<T>(r, c, val));
}

template <typename T>
void MatrixEngine<T, GENERAL, ASYMMETRIC, CRS>::sort_()
{
	//std::cout << "crs::sort_()" << std::flush;
	
	std::sort(coord_vals_.begin(), coord_vals_.end(), less_);
	boost::timer t;
	//typename std::vector<coord_tuple<T> >::iterator it;
	//typename std::vector<coord_tuple<T> >::iterator it_end = coord_vals_.end();
	//typename std::vector<coord_tuple<T> >::iterator it2;
	sorted_=true;
	int toDelete=0;
	t.restart();
	//std::cout << " .... sorted!" << std::flush;
	/*for(it=coord_vals_.begin(), it2=it+1; it2!=it_end; ++it, ++it2)
	{
		while(it2!=it_end && (!less_(*it, *it2)))
		{
			it->value += it2->value;
            it2->value = T(0);
            ++it2;
		}
		if(it2!=it_end)
			*(it+1) = *it2;
	}*/
	int i, I;
	for(i=0, I=1; I<coord_vals_.size(); ++i,++I) {
		while(I<coord_vals_.size() && (!less_(coord_vals_[i], coord_vals_[I]))) {
			coord_vals_[i].value += coord_vals_[I].value;
			coord_vals_[I].value = T(0);
			++I;
		}
		if(I<coord_vals_.size()) {
			coord_vals_[i+1] = coord_vals_[I];
		}
	}
	
	//std::cout << " .... added!" << std::flush;
	
	//if ((it+1!=it_end) && (it2-it-1 > 0))
	//	coord_vals_.erase(it+1, it_end);
	
	if(i<coord_vals_.size() && (I-i-1>0))
		coord_vals_.erase(coord_vals_.end()-(I-i-1), coord_vals_.end());
	//std::cout << " .... erased!" << std::endl;
	
}

template <typename T>
void MatrixEngine<T, GENERAL, ASYMMETRIC, CRS>::finalize()
{
	assert(!finalized_);
	sort_();
	
	assemble_crs_();	
	coord_vals_.clear();
	coord_vals_.reserve(0);
	finalized_=true;
}

template <typename T>
void MatrixEngine<T, GENERAL, ASYMMETRIC, CRS>::allocate_(const size_type size)
{
	
	values_= new T[size];
	values_limit_ = values_ + size;
	
	row_ptrs_= new size_type[nRows_+1];
	row_ptrs_limit_ = row_ptrs_ + nRows_ + 1;
	
	//col_indices_= row_ptrs_ + nRows_ + 1;
	col_indices_ = new size_type[size];
	col_indices_limit_ = col_indices_ + size;	

}

template <typename T>
void MatrixEngine<T, GENERAL, ASYMMETRIC, CRS>::release_()
{
	if (values_) {
		delete [] values_;
	}
	// reset pointers to indicate that the `Vec' is empty again
	values_ = values_limit_ = 0;
	
	
	if (row_ptrs_) {
		delete [] row_ptrs_;
	}
	
	if (col_indices_) {
		delete [] col_indices_;
	}
	// reset pointers to indicate that the `Vec' is empty again
	row_ptrs_ = row_ptrs_limit_ = 0;
	col_indices_ = col_indices_limit_ = 0;
	
}

template <typename T>
void MatrixEngine<T, GENERAL, ASYMMETRIC, CRS>::assemble_crs_()
{
	if(coord_vals_.size()==0) {
		values_=0;
		values_limit_=0;
		row_ptrs_ = 0;
		row_ptrs_limit_ = 0;
		col_indices_ = 0;
		col_indices_limit_ = 0;
		return;
	}
	allocate_(coord_vals_.size());
	coord_iterator colIt=col_indices_, rowIt=row_ptrs_, lastRowIt=row_ptrs_;
	val_iterator valIt=values_;
	typename std::vector<coord_tuple<T> >::iterator it=coord_vals_.begin();
	int index=0, last_row=-1, totalCount=0, lr=0;
	int i=0;
	while(it!=coord_vals_.end()) {
		if(last_row!=it->row) {
			for(i=last_row+1; i<=it->row; ++i)
				row_ptrs_[i] = totalCount;
			last_row = it->row;
		}
		*(colIt++)=it->col;
		*(valIt++)=it->value;
		totalCount++;
		row_ptrs_[it->row+1] = totalCount;
		lr = it->row+1;
		index++;
		it++;
	}
	for(i=lr+1; i<row_ptrs_limit_-row_ptrs_; ++i)
		row_ptrs_[i] = totalCount;
}

template <typename T>
void MatrixEngine<T, GENERAL, ASYMMETRIC, CRS>::print(std::ostream &o) const
{
	o << "size= [" << nRows_ << ", " << nCols_ << "], nnz~= ";
	if(!finalized_) {
		o << coord_vals_.size() << " [NOT FINALIZED]\n";
		typename std::vector<coord_tuple<T> >::const_iterator it=coord_vals_.begin();
		while(it!=coord_vals_.end()) {
			o << "  (" << it->row << ", " << it->col << ") : " << it->value << "\n";
			++it;
		}
	} else {
		o << col_indices_limit_-col_indices_ << "\n";
		int last_row_ind = 0;
		for(size_t i=0; i<row_ptrs_limit_-row_ptrs_-1; ++i) {
			int j=row_ptrs_[i];
			while(j<row_ptrs_[i+1]) {
				o << "  (" << i+1 << ", " << col_indices_[j]+1 << ") : " << values_[j] << "\n";
				++j;
			}
		}
	}
}

template <typename T> template <typename V1, typename V2>
inline void MatrixEngine<T, GENERAL, ASYMMETRIC, CRS>::mv(const V1& x, V2& y, bool trans, T alpha, T beta) const
{
	assert(finalized_);
	assert(trans ? (x.size()==nRows_ && y.size()==nRows_) : (x.size()==nCols_ && y.size()==nRows_) );
	if(values_ == 0) {
	// [TODO]: mark 2010-06-14- clean this up (implement sparse_lib::vector = scalar)
		for(int i=0; i<x.size(); ++i)
			y(i) = 0.0;
	} else {
		crs_mv(trans, nRows_, nCols_, alpha, values_, row_ptrs_, col_indices_, x.data(), beta, y.data());
	}
}

template <typename T>
void MatrixEngine<T, GENERAL, ASYMMETRIC, CRS>::resize(size_type nr, size_type nc)
{
	if(finalized_)
		std::cerr << "Finalized!" << std::endl;
	assert(!finalized_);
	if(nr >= nRows_ && nc >= nCols_) {
		nRows_=nr;
		nCols_=nc;
	} else {
		std::cerr << "Can not shrink a CRS yet!" << std::endl;
	}
}

} /* sparse_lib */
