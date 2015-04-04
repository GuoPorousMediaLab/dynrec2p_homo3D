#ifndef GIVENS_SPARSELIB
#define GIVENS_SPARSELIB

// [TODO]: mark 2010-06-21- VERY IMPORTANT: Create own implementation of givens_rotation before any kind of source code release!!!!!
namespace sparse_lib
{

template <typename T>
class givens_rotation
{
public:
	inline givens_rotation (): a_(0), b_(0), c_(0), s_(0) { }
	inline givens_rotation(const T& a_in, const T& b_in) {
		T roe;
	    if (std::abs(a_in) > std::abs(b_in))
	      roe = a_in;
	    else
	      roe = b_in;

	    T scal = std::abs(a_in) + std::abs(b_in);
	    T r, z;
	    if (scal != T(0)) {
	      T a_scl = a_in / scal;
	      T b_scl = b_in / scal;
	      r = scal * std::sqrt(a_scl * a_scl + b_scl * b_scl);
	      if (roe < T(0)) r *= -1;
	      c_ = a_in / r;
	      s_ = b_in / r;
	      z = 1;
	      if (std::abs(a_in) > std::abs(b_in))
	        z = s_;
	      else if (std::abs(b_in) >= std::abs(a_in) && c_ != T(0))
	        z = T(1) / c_;
	    } else {
	      c_ = 1; s_ = 0; r = 0; z = 0;      
	    }
	    a_ = r;
	    b_ = z;
	}

  inline void scalar_apply(T& x, T& y) {
    T tmp = c_ * x + s_ * y;
    y = c_ * y - s_ * x;
    x = tmp;
  }

private:
	T a_, b_, c_, s_;
	
};
	
} /* sparse_lib */ 

#endif /* end of include guard: GIVENS_SPARSELIB */
