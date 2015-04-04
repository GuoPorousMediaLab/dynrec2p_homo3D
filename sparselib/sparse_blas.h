#ifndef SPARSE_LIB_SPARSE_BLAS_H
#define SPARSE_LIB_SPARSE_BLAS_H

#include <boost/timer.hpp>
#ifdef HAVE_OMP
#include <omp.h>
#endif

namespace sparse_lib
{

template <typename T>
void crs_mv(bool trans, const int m, const int n, T alpha,
         const T *a, const int *ia, const int *ja,
         const T *x, T beta, T *y)
{
	const bool init = (beta==T(0));
	T sum;
	int beg, end, k;
	if (!trans) {
		#pragma omp parallel for private(sum, end, beg, k)
		for (int i=0; i<m; ++i) {
			if (init) {
				y[i] = 0;
			}
			sum = 0.0;
			beg = ia[i];
			end = ia[i+1];
			for(k=beg; k<end; ++k) {
				sum+=a[k]*x[ja[k]];
			}
			y[i] += alpha*sum;
		}
	} else {
		if (init) {
			for (int i=0; i<n; ++i) {
				y[i] = 0;
			}
		}
		for (int i=0; i<m; ++i) {
			for (int k=ia[i]; k<ia[i+1]; ++k) {
				y[ja[k]] += alpha*a[k]*x[i];
			}
		}
	}
}

template <typename T>
void cds_mv(bool trans, const int m, const int n, T alpha,
	const T *a, const int ndiag, const int *diag,
	const T *x, T beta, T *y)
{
	if(beta==T(0)) {
		for(int i=0; i<m; ++i)
			y[i]=0;
	}
	if(!trans) {
		for(int i=0; i<ndiag; i++) {
			int d = diag[i];
			int start = std::max(0, -d);
			int end = std::min(n, n-d);
			for(int loc = start; loc<end; ++loc)
				y[loc]+=a[i*m+loc]*x[loc+d];
		}
	} else {
		
	}
}

template <typename T>
void cds_sym_mv(bool trans, const int m, const int n, T alpha,
	const T * __restrict__ a, const int ndiag, const int *diag,
	const T * __restrict__ x, T beta, T * __restrict__ y)
{
	boost::timer t;
	int i=0;
	if(beta==T(0)) {
		#pragma omp parallel for private(i)
		for(int i=0; i<m; ++i)
			y[i]=0;
	}
	int operations = 0;
	i=0;
	
	/*int threads = 1;
	#ifdef HAVE_OMP
	threads = omp_get_num_threads();
	std::cout << "--------------" << std::endl;
	#endif
	#pragma omp parallel for shared(a, x, y)
	for(int k=0;k<threads;k++)  {
		int n1 = n/threads;
		int n2 = n%threads;
		int k_start = k<n2 ? (n1+1)*k : n1*k+n2;
		int k_end = k<n2 ? k_start+n1+1  : k_start+n1;
		std::cout << "from " << k_start << " to " << k_end << std::endl;
		for(int j=0;j<ndiag;j++) {
			int offset = j*n;
			int d = diag[j];
			int start = std::max(k_start, -d);
			int end = std::min(k_end, n-d);
			//T const* ao = a+offset;
			//T const* xo = x+d;
			for(int loc = start; loc<end; ++loc)
				#pragma omp atomic
				y[loc]+=a[loc+offset]*x[loc+d];
			
			if(d!=0) {
				d = -d;
				int start2 = std::max(k_start, -d);
				int end2 = std::min(k_end, n-d);
				for(int loc = start+offset, loc2=start2; loc<end+offset; ++loc, ++loc2)
					#pragma omp atomic
					y[loc2]+=a[loc]*x[loc2+d];
			}
	    }
	}*/
	
	for(i=0; i<ndiag; i++) {
		int offset = i*m;
		int d = diag[i];
		int start = std::max(0, -d);
		int end = std::min(n, n-d);
		//operations += 2*(end-start-1);
		
		T const* ao = a+offset;
		T const* xo = x+d;
		for(int loc = start; loc<end; ++loc)
			y[loc]+=ao[loc]*xo[loc];
	}
	
	for(i=0; i<ndiag; ++i) {
		int d = diag[i];
		if(d!=0) {
			int offset = i*m;
			int start = std::max(0, -d);
			int end = std::min(n, n-d);
			//operations += 2*(end-start-1);

			T const* ao = a+offset;
			T const* xo = x+d;
			
			d=-d;
			int start2 = std::max(0, -d);
			int end2 = std::min(n, n-d);
			//operations += 2*(end-start-1);
			for(int loc = start+offset, loc2=start2; loc<end+offset; ++loc, ++loc2)
				y[loc2]+=a[loc]*x[loc2+d];
		}
	}
	
	//std::cout << "Operations: " << operations << " " << m << std::endl;
	//std::cout << "MFlops (SpMv): " << (double)operations/t.elapsed() * 1.e-6 << std::endl;
}

template <typename T>
void cds_sv(bool trans, bool upper, const int m, const T *a, 
				const int ndiag, const int *diag,
				T *x)
{
	if(upper) {
		// [TODO]: mark 2009-09-13- FIX UPPER!!!!!!!
		/*const int * it_d = diag;
		const int * it_end = diag + ndiag;
		int diagCount=0;
		while(*it_d<0)
			it_d++;
		assert(*it_d == 0);
		if(!trans) {
			for(int i=m-1; i>=0; --i) {
				const int * it = it_d+1;
				while(*it > i && it!=it_end) {
					x[i]-=a[ (it-diag)*m + i ]*x[*it+i];
					++it;
				}
				x[i]/=a[(it_d-diag)*m+i];
			}
		} else {
			for(int i=0; i<m; ++i) {
				const int * it = it_d+1;
				while(*it <= i && it!=it_end) {
					x[i]-=a[(it-diag)*m + i - *it]*x[i-*it];
					++it;
				}
				x[i]/=a[(it_d-diag)*m+i];
			}
		}*/
		assert("Mark- you need to fix the upper triangular solve!" && false);
	} else {
		const int * it_d = diag;
		const int * it_end = diag + ndiag;
		int diagCount=0;
		while(*it_d<0)
			it_d++;
		assert(*it_d == 0);
		if(!trans) {
			for(int i=0; i<m; ++i) {
				const int * it = it_d-1;
				while(it >= diag && -*it <= i) {
					x[i]-=a[(it-diag)*m + i]*x[i+*it];
					//std::cout << "x[" << i << "]-=" << a[(it-diag)*m + i] << "* x[" << i+*it << "]" << std::endl;
					--it;
				}
				x[i]/=a[(it_d-diag)*m+i];
			}
		} else {
			for(int i=m-1; i>=0; --i) {
				const int * it = it_d-1;
				while(it>=diag && -*it <= m-1-i) {
					x[i]-=a[(it-diag)*m + i - *it]*x[(i-*it)];
					//std::cout << "x[" << i << "]-=" << a[(it-diag)*m + i - *it] << "* x[" << (i-*it) << "]" << std::endl;
					
					--it;
				}
				x[i]/=a[(it_d-diag)*m+i];
			}
		}
	}
}


} /* sparse_lib */ 

#endif /* end of include guard: SPARSE_BLAS_H */
