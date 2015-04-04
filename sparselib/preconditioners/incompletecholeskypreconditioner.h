#ifndef SPARSE_LIB_INCOMPLETECHOLESKYPRECONDITIONER_H
#define SPARSE_LIB_INCOMPLETECHOLESKYPRECONDITIONER_H

#include <vector>
#include <algorithm>
#include <boost/timer.hpp>

namespace sparse_lib
{
	
template <typename T, MatrixType Tri, Symmetry Sym, MatrixStorage MatStorage>
class ICCPreconditionerBase: public Preconditioner<T>
{
public:
	ICCPreconditionerBase () { assert("Unsupported matrix type." && false); }
	~ICCPreconditionerBase () { }

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

template <typename T>
class ICCPreconditionerBase<T, UPPER, SYMMETRIC, CDS>: public Preconditioner<T>
{
public:
	typedef Matrix<MatrixEngine<T, UPPER, SYMMETRIC, CDS> > MatrixType;
	ICCPreconditionerBase () { }
	~ICCPreconditionerBase () { }

	void set(const MatrixBase<T>* Mat);
	void set(const Matrix<MatrixEngine<T, UPPER, SYMMETRIC, CDS> >& Mat);

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
	Matrix<MatrixEngine<T, LOWER, ASYMMETRIC, CDS> > CF_;
	void factor_(const Matrix<MatrixEngine<T, UPPER, SYMMETRIC, CDS> >&);
	int get_natural_index_(int r, int c);
	template <typename V1, typename V2>
	V2& apply_(const V1& b, V2& x);
	std::vector<T> F_;
	std::vector<int> D_ind_, D_row_ind_;
	int FT_start;
	int Diag_start;
};

template <typename T>
void ICCPreconditionerBase<T, UPPER, SYMMETRIC, CDS>::set(const MatrixBase<T>* Mat)
{
	const Matrix<MatrixEngine<T, UPPER, SYMMETRIC, CDS> >* m = dynamic_cast<const Matrix<MatrixEngine<T, UPPER, SYMMETRIC, CDS> >*>(Mat);
	if(m) {
		set(*m);
	} else {
		assert("Wrong Matrix type.  Must be a symmetrical upper CDS." && false);
	}
}

template <typename T>
void ICCPreconditionerBase<T, UPPER, SYMMETRIC, CDS>::set(const Matrix<MatrixEngine<T, UPPER, SYMMETRIC, CDS> >& Mat)
{
	factor_(Mat);
}

/*template <typename T>
void ICCPreconditionerBase<T, UPPER, SYMMETRIC, CDS>::factor_(const Matrix<MatrixEngine<T, UPPER, SYMMETRIC, CDS> >& Mat)
{
	int rows = Mat.rows();
	int cols = Mat.cols();
	int nDiags = Mat.engine().nDiags();
	CF_.engine().release_();
	CF_.engine().nRows_=rows;
	CF_.engine().nCols_=cols;
	CF_.engine().allocate_(nDiags);
	const int *MatDiags = Mat.engine().diags();
	const int *m_d_it = MatDiags + nDiags - 1;
	assert(nDiags > 1);
	int size = 0;
	int inset = 0;
	D_ind_.push_back(rows);
	
	while(m_d_it >= MatDiags) {
		if(*m_d_it > 0) {
			D_ind_.push_back(*m_d_it);
			size+=2*(rows-(*m_d_it));
		} else
			size+=(rows);
		inset+=*m_d_it;
		CF_.engine().insert_diag_(-*m_d_it--);
	}
	std::reverse(D_ind_.begin(), D_ind_.end());
	//std::cout << "Vector size: " << size << std::endl;
	//std::cout << "rows: " << rows << std::endl;
	//std::cout << "F_start: " << 0 << std::endl;
	//std::cout << "FT_start: " << rows*2-inset << std::endl;
	//std::cout << "Diag starg: " << rows*4-inset*2 << std::endl;
	int FT_start = rows*2-inset;
	Diag_start = rows*4-inset*2;
	
	F_.clear(); F_.resize(size);
	int * mainDiag = CF_.engine().diags_ + nDiags - 1;
	int * diagEnd = CF_.engine().diags_;
	T * values = CF_.engine().values_;
	int entries = 0;
	for(int i=0; i<rows; ++i) {
		T val = Mat(i,i);
		int * d_it = mainDiag - 1;
		while( d_it >= diagEnd && -(*d_it) <= i ) {
			val-=std::pow(values[(d_it-diagEnd)*rows+i],2);
			--d_it;
		}

		T diagVal = std::sqrt(val);
		CF_(i,i) = diagVal;
		//std::cout << "CF_(" << i << "," << i << ") = " << diagVal << std::endl;;
		F_[Diag_start+rows-(i+1)] = diagVal;
		d_it = mainDiag - 1;
		for(d_it = mainDiag-1; d_it>=diagEnd && -(*d_it) < rows-i; --d_it) {
		//while( d_it >= diagEnd && -(*d_it) < rows-i) {
			int * i_diag=d_it, * j_diag;
			val = Mat(i,-(*d_it)+i);
			for(i_diag=d_it; i_diag > diagEnd && -(*i_diag) <=i; --i_diag) {
			//while(i_diag >= diagEnd && -(*i_diag) <= i) {
				j_diag = i_diag - 1;
				if( (*(i_diag) - *(j_diag)) == -(*d_it) )
					val-=values[(i_diag-diagEnd)*rows+i]*values[(j_diag-diagEnd)*rows+i-(*d_it)];
				//--i_diag;
			}
			int r = -(*d_it)+i;
			int c = i;
			
			int ind = 0;
			for(int* d = CF_.engine().diags_; d < CF_.engine().diags_ + nDiags; ++d) {
				if(*d != 0) {
					ind+=std::max(0,(r+*d));
					if(*d < c-r && r+*d>=0)
						ind++;
				}
			}
			
			F_[ind] = val/diagVal;
			F_[Diag_start - 1 - entries] = val/diagVal;
			
			
			CF_(-(*d_it)+i, i) = val/diagVal;
			//std::cout << "CF_(" << -(*d_it)+i << "," << i << ") = " << val/diagVal << "   ";
			//std::cout << "F_[" << ind << "] = " << val/diagVal << "   ";
			//std::cout << "F_[" << Diag_start - 1 - entries << "] = " << val/diagVal << "   ";
			
			//std::cout << std::endl;
			
			//--d_it;
			++entries;
		}
	}
}*/

template <typename T>
int ICCPreconditionerBase<T, UPPER, SYMMETRIC, CDS>::get_natural_index_(int r, int c)
{
	int index = 0;
	int row_begin = 0;
	if(r > c) {
		int diag = r-c;
		//std::vector<int>::iterator sub_it = std::upper_bound(D_row_ind_.begin()+1, D_row_ind_.end(), r);
		std::vector<int>::iterator sub_it = std::lower_bound(D_row_ind_.begin()+2, D_row_ind_.end(), r);
		row_begin = 0;
		for(std::vector<int>::iterator it = D_row_ind_.begin()+1; it<sub_it; ++it) {
			row_begin += r-*it;
		}
		sub_it = std::upper_bound(D_row_ind_.begin()+2, D_row_ind_.end(), r);
		int colAdd = sub_it - std::upper_bound(D_row_ind_.begin()+2, D_row_ind_.end(), diag);
		index = row_begin + colAdd;
		//std::cout << "(" << r << ", " << c << "): "<< row_begin << ", "<< colAdd << std::endl;
		
	} else if(r < c) {
		int last = F_.size()-1;
		int rows = D_row_ind_.back() - 1;
		row_begin = FT_start;
		int diag = c-r;
		/*std::vector<int>::iterator sub_it = std::upper_bound(D_row_ind_.begin()+1, D_row_ind_.end(), rows-r);
		for(std::vector<int>::iterator it = D_row_ind_.begin()+1; it<sub_it; ++it) {
			//std::cout << "row_begin += " << (rows-r)-*it << " ";
			row_begin += (rows-r)-*it;
		}
		//std::cout << std::endl;
		int colAdd = sub_it - std::upper_bound(D_row_ind_.begin(), D_row_ind_.end(), diag);
		//std::cout << "r: " << r << "  c: " << c << "   colAdd: " << colAdd << std::endl;
		//std::cout << "(" << r << ", " << c << "): "<< row_begin << ", "<< colAdd << "    " << *sub_it << "   diag: " << diag << std::endl;
		index = row_begin + colAdd;*/
		
		std::vector<int>::iterator sub_it = std::lower_bound(D_row_ind_.begin()+2, D_row_ind_.end(), rows-r);
		for(std::vector<int>::iterator it = D_row_ind_.begin()+1; it<sub_it; ++it) {
			row_begin += (rows-r)-*it;
		}
		sub_it = std::upper_bound(D_row_ind_.begin()+2, D_row_ind_.end(), rows-r);
		int colAdd = sub_it - std::upper_bound(D_row_ind_.begin()+2, D_row_ind_.end(), diag);
		index = row_begin + colAdd;
		//std::cout << "(" << r << ", " << c << "): "<< row_begin << ", "<< colAdd << std::endl;
	} else {
		index =  FT_start - (D_row_ind_.back()) + r;
	}
	
	if(index < 0 || index >= F_.size())
		std::cout << r << ", " << c << " --> " << index << "  max is " << F_.size()-1 << "   row_begin: " << row_begin <<  std::endl;
	return index;
}

template <typename T>
void ICCPreconditionerBase<T, UPPER, SYMMETRIC, CDS>::factor_(const Matrix<MatrixEngine<T, UPPER, SYMMETRIC, CDS> >& Mat)
{
	int rows = Mat.rows();
	int cols = Mat.cols();
	assert(rows == cols);
	int nDiags = Mat.engine().nDiags();
	//CF_.engine().release_();
	//CF_.engine().nRows_=rows;
	//CF_.engine().nCols_=cols;
	//CF_.engine().allocate_(nDiags);
	const int *MatDiags = Mat.engine().diags();
	const int *m_d_it = MatDiags;
	assert(nDiags > 1);
	int size = 0;
	int inset = 0;
	
	while(m_d_it < MatDiags + nDiags) {
		if(*m_d_it > 0) {
			size+=2*(rows-(*m_d_it));
		} else
			size+=(rows);
		D_ind_.push_back(*m_d_it);
		D_row_ind_.push_back(*m_d_it);
		inset+=*m_d_it;
		//CF_.engine().insert_diag_(-*m_d_it--);
		m_d_it++;
	}
	D_row_ind_.push_back(rows);
	FT_start = rows*(nDiags)-inset;	
	
	F_.clear(); F_.resize(size);
	
	//int * mainDiag = CF_.engine().diags_ + nDiags - 1;
	//int * diagEnd = CF_.engine().diags_;
	//T * values = CF_.engine().values_;
	
	int entries = 0;
	for(int j=0; j<rows; ++j) {
		T Dval = Mat(j,j);
		
		//std::cout << "F(" << j << ", " << j << "): index = " << get_natural_index_(j, j) << std::endl;
		
		std::vector<int>::iterator k_it = D_row_ind_.begin()+1;
		
		while((j - *k_it) >= 0) {
			int k = j - *k_it;
			//std::cout << "D[" << i << "] = A[" << i << ", " << i << "] - L[" << i << ", " << k << "]*D[" << k << ", " << k << "]\n";
			Dval -= std::pow(F_[get_natural_index_(j,k)], 2) * F_[get_natural_index_(k,k)];
			++k_it;
		}
		//std::cout << "D_" << j << " = " << Dval << std::endl; 
		F_[get_natural_index_(j,j)] = Dval;
		
		std::vector<int>::iterator i_it = D_ind_.begin()+1;
		while(i_it < D_ind_.end() && *i_it+j <rows) {
			int i = *i_it + j;
			T val = Mat(i,j);
			k_it = D_row_ind_.begin()+1;
			int diag_ij = i-j;
			std::vector<int>::iterator diag_it = std::find(D_ind_.begin(), D_ind_.end(), diag_ij)+1;
			while(diag_it != D_ind_.end()) {
				int k = i - *diag_it;
				if(k<0)
					break;
				int diag_jk = j-k;
				if(std::find(D_ind_.begin(), D_ind_.end(), diag_jk) != D_ind_.end()) {
					std::cout << "Found!  " << diag_jk << "    " << diag_ij << "   " << "i: " << i << "  j: " << j << "   k: " << k << std::endl;
					val -= F_[get_natural_index_(i,k)] * F_[get_natural_index_(j,k)] * F_[get_natural_index_(k,k)];
				}
				++diag_it;
			}
			// while(k <= j-1) {
			// 	int diag1 = i-k;
			// 	int diag2 = j-k;
			// 	if(std::find(D_ind_.begin(), D_ind_.end(), diag1) != D_ind_.end() && std::find(D_ind_.begin(), D_ind_.end(), diag2) != D_ind_.end()) {
			// 		val -= F_[get_natural_index_(i,k)] * F_[get_natural_index_(j,k)]*F_[get_natural_index_(k,k)];
			// 	}
			// 	++k;
			// }
			F_[get_natural_index_(i,j)] = val/F_[get_natural_index_(j,j)];
			F_[get_natural_index_(j,i)] = val/F_[get_natural_index_(j,j)];
			//std::cout << "L[" << i << ", " << j << "] = " << val/F_[get_natural_index_(j,j)] << "\n";
			++i_it;
		}
	}
	
	for(int j=0; j<rows; ++j)
		F_[get_natural_index_(j,j)] = 1.0/F_[get_natural_index_(j,j)];
	//for(int i=0; i<F_.size(); ++i)
	//	std::cout << F_[i] << " ";
	//std::cout << std::endl;
}



template <typename T> template <typename V1, typename V2>
V2& ICCPreconditionerBase<T, UPPER, SYMMETRIC, CDS>::apply_(const V1& b, V2& x)
{
	int operations = 0;
	boost::timer t;
	DenseVector<T> y(b.size());
	//DenseVector<T> x2(b.size());
	//CF_.solve(y,b);
	//CF_.solve(x,y,true);
	
	
	int i=1;
	int ind=0;
	T* __restrict__ F = &F_[0];
	
	
	y(0)=b(0);
	int DS = D_row_ind_.size();
	double *insert = &y(1);
	for(int j=1; j<D_row_ind_.size(); ++j) {
		while(i<D_row_ind_[j]) {
			//std::cout << "y(" << i << ") = ( b(" << i << ")";
			double sum = b(i);
			for(int k=(j-1); k>=1; --k) {
				//std::cout << " - y(" << i-D_row_ind_[k] << ")*" << F[ind];
				sum -= y(i-D_row_ind_[k]) * *(F++);
			}
			y(i) = sum;
			++i;
		}
	}
	
	T* yp = &y(0);
	
	for(int j=0; j<y.size(); ++j)
		*(yp++) *= *(F++);
	
	//std::cout << y << std::endl;
	
	i = 1;
	int sz = x.size()-1;
	x(x.size()-1)=y(y.size()-1);
	for(int j=1; j<D_row_ind_.size(); ++j) {
		while(i<D_row_ind_[j]) {
			//std::cout << "x(" << (sz-i) << ") = ( y(" << (sz-i) << ")";
			double sum = y((sz-i));
			for(int k=(j-1); k>=1; --k) {
				//std::cout << " - x(" << sz-i+D_row_ind_[k] << ")*" << F[ind];
				sum -= x(sz-i+D_row_ind_[k]) * *(F++);
			}
			x((sz-i)) = sum;
			//std::cout << std::endl;
			++i;
		}
	}
	
	ind--;
	
	
	operations = 2*(F-&F_[0]) - x.size();
	//std::cout << "MFlops (tri solve) = " << (double)operations/t.elapsed()*1.0e-6 << "   operations = " << operations << std::endl;
	//std::cout << nrm2(x2) << " " << nrm2(x) << std::endl;
	return x;
}

template <typename M>
class ICCPreconditioner: public ICCPreconditionerBase<typename M::value_type, M::MT, M::SYM, M::ST> {

};


} /* sparse_lib */ 

#endif /* end of include guard: INCOMPLETECHOLESKYPRECONDITIONER_H */
