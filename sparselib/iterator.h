#ifndef SPARSELIB_ITERATOR
#define SPARSELIB_ITERATOR

#include <iostream>
#include <cmath>

namespace sparse_lib
{

enum iter_output { ITER_NONE,		// No output.
				   ITER_SUMMARY,	// Output summary at finished state.
				   ITER_ALL			// Output info at each iteration.
				 };

template <typename T>
class Iterator
{
public:
	Iterator (T rel_tol, int max_it, T abs_tol = (T)0.0, iter_output log_level=ITER_NONE, std::ostream &out=std::cout)
		: rel_tol_(rel_tol), abs_tol_(abs_tol), max_it_(max_it), log_level_(log_level), o_(out), b_norm_(1.0), r_(1.e99), it_(0) { }
	~Iterator() { }

	Iterator& operator++() { ++it_; return(*this); }
	Iterator operator++(int) { Iterator Iter(*this); ++(*this); return Iter; }

	bool hasConverged() const { return (r_/b_norm_ <= rel_tol_ || r_ <= abs_tol_); }

	bool finished() const {
		if(hasConverged()) {
			if(log_level_ == ITER_SUMMARY || log_level_ == ITER_ALL)
				o_ << id_ << " Iterator CONVERGED at step " << it_ << " with an absolute residual of " << r_ << " and a relative residual of " << r_/b_norm_ << std::endl;
			return true;
		} else if (it_ > max_it_) {
			if(log_level_ == ITER_SUMMARY || log_level_ == ITER_ALL)
				o_ << id_ << " Itereator FAILED TO CONVERGE before reaching max_it.  Absolute residual was " << r_ << " with a relative residual of " << r_/b_norm_ << std::endl;
			return true;
		} else {
			if(log_level_ == ITER_ALL)
				o_ << id_ << " Iterator at step " << it_ << " with an absolute residual of " << r_ << " and a relative residual of " << r_/b_norm_ << std::endl;
			return false;
		}
		
	}
	template <typename V>
	bool finished(const V& v) { r_ = std::abs(nrm2(v)); return finished(); }
	bool finished(T r) { r_ = r; return finished(); }
	
	int getStep() const { return it_; }
	void setStep(int ns) { it_ = ns; }
	
	T getAbsoluteResidual() const { return r_; }
	T getRelativeResidual() const { return r_/b_norm_; }
	
	T getBNorm() const { return b_norm_; }
	void setBNorm(T b_norm) { b_norm_ = b_norm; }
	
	void silence() { log_level_ = ITER_NONE; }
	iter_output getOutputLevel() { return log_level_; }
	void setOutputLevel(iter_output log_level) { log_level_ = log_level; }
	
	int getMaxIter() { return max_it_; }
	void setMaxIter(int maxIt) { max_it_ = maxIt; }
	
	std::string getId() const { return id_; }
	void setId(std::string id) { id_ = id; }
	
private:
	T rel_tol_;						// relative tolerance
	T abs_tol_;						// absolute tolerance
	int max_it_;					// maximum iterations
	iter_output log_level_;			// verbosity
	std::ostream &o_;				// Where to output to
	
	T b_norm_;						// RHS norm
	T r_;							// Last computed residual
	int it_;						// number of iterations
	std::string id_;
};

} /* sparse_lib */ 

#endif /* end of include guard: SPARSELIB_ITERATOR */
