#ifndef SPARSE_LIB_VVOPS_H
#define SPARSE_LIB_VVOPS_H

#include <sparselib/densevector.h>
#include <cmath>

namespace sparse_lib
{
	
template <typename T>
DenseVectorTemp<T> operator+(const DenseVector<T>&, const DenseVector<T>&);
template <typename T>
DenseVectorTemp<T> operator+(const DenseVector<T>&, const DenseVectorView<T>&);
template <typename T>
DenseVectorTemp<T> operator+(const DenseVectorView<T>&, const DenseVector<T>&);
template <typename T>
DenseVectorTemp<T> operator+(const DenseVectorView<T>&, const DenseVectorView<T>&);
template <typename T>
DenseVectorTemp<T> operator+(const DenseVector<T>&, const DenseVectorTemp<T>&);
template <typename T>
DenseVectorTemp<T> operator+(const DenseVectorTemp<T>&, const DenseVector<T>&);
template <typename T>
DenseVectorTemp<T> operator+(const DenseVectorView<T>&, const DenseVectorTemp<T>&);
template <typename T>
DenseVectorTemp<T> operator+(const DenseVectorTemp<T>&, const DenseVectorView<T>&);
template <typename T>
DenseVectorTemp<T> operator+(const DenseVectorTemp<T>&, const DenseVectorTemp<T>&);


template <typename T>
DenseVectorTemp<T> operator-(const DenseVector<T>&, const DenseVector<T>&);
template <typename T>
DenseVectorTemp<T> operator-(const DenseVector<T>&, const DenseVectorView<T>&);
template <typename T>
DenseVectorTemp<T> operator-(const DenseVectorView<T>&, const DenseVector<T>&);
template <typename T>
DenseVectorTemp<T> operator-(const DenseVectorView<T>&, const DenseVectorView<T>&);
template <typename T>
DenseVectorTemp<T> operator-(const DenseVector<T>&, const DenseVectorTemp<T>&);
template <typename T>
DenseVectorTemp<T> operator-(const DenseVectorTemp<T>&, const DenseVector<T>&);
template <typename T>
DenseVectorTemp<T> operator-(const DenseVectorView<T>&, const DenseVectorTemp<T>&);
template <typename T>
DenseVectorTemp<T> operator-(const DenseVectorTemp<T>&, const DenseVectorView<T>&);
template <typename T>
DenseVectorTemp<T> operator-(const DenseVectorTemp<T>&, const DenseVectorTemp<T>&);


template <typename T>
DenseVectorTemp<T> operator*(const T, const DenseVector<T>&);
template <typename T>
DenseVectorTemp<T> operator*(const DenseVector<T>&, const T);
template <typename T>
DenseVectorTemp<T> operator*(const T, const DenseVectorView<T>&);
template <typename T>
DenseVectorTemp<T> operator*(const DenseVectorView<T>&, const T);
template <typename T>
DenseVectorTemp<T> operator*(const T, const DenseVectorTemp<T>&);
template <typename T>
DenseVectorTemp<T> operator*(const DenseVectorTemp<T>&, const T);


template <typename T>
DenseVectorTemp<T> operator/(const DenseVector<T>&, const T);
template <typename T>
DenseVectorTemp<T> operator/(const DenseVectorView<T>&, const T);
template <typename T>
DenseVectorTemp<T> operator/(const DenseVectorTemp<T>&, const T);

template <typename T>
T operator*(const DenseVector<T>& lhs, const DenseVector<T>& rhs);

template <typename V>
typename V::value_type nrm2(const V &x);

template <typename T>
DenseVectorTemp<T> operator+(const DenseVector<T>& lhs, const DenseVector<T>& rhs)
{
	DenseVectorTemp<T> ret(lhs);
	ret+=rhs;
	return ret;
}

template <typename T>
DenseVectorTemp<T> operator+(const DenseVector<T>& lhs, const DenseVectorView<T>& rhs)
{
	DenseVectorTemp<T> ret(lhs);
	ret+=rhs;
	return ret;
}

template <typename T>
DenseVectorTemp<T> operator+(const DenseVectorView<T>& lhs, const DenseVector<T>& rhs)
{
	DenseVectorTemp<T> ret(lhs);
	ret+=rhs;
	return ret;
}

template <typename T>
DenseVectorTemp<T> operator+(const DenseVectorView<T>& lhs, const DenseVectorView<T>& rhs)
{
	DenseVectorTemp<T> ret(lhs);
	ret+=rhs;
	return ret;
}

template <typename T>
DenseVectorTemp<T> operator+(const DenseVector<T>& lhs, const DenseVectorTemp<T>& rhs)
{
	DenseVectorTemp<T> ret(rhs);
	return ret+=lhs;
}

template <typename T>
DenseVectorTemp<T> operator+(const DenseVectorTemp<T>& lhs, const DenseVector<T>& rhs)
{
	DenseVectorTemp<T> ret(lhs);
	return ret+=rhs;
}

template <typename T>
DenseVectorTemp<T> operator+(const DenseVectorView<T>& lhs, const DenseVectorTemp<T>& rhs)
{
	DenseVectorTemp<T> ret(rhs);
	return ret+=lhs;
}

template <typename T>
DenseVectorTemp<T> operator+(const DenseVectorTemp<T>& lhs, const DenseVectorView<T>& rhs)
{
	DenseVectorTemp<T> ret(lhs);
	return ret+=rhs;
}

template <typename T>
DenseVectorTemp<T> operator+(const DenseVectorTemp<T>& lhs, const DenseVectorTemp<T>& rhs)
{
	DenseVectorTemp<T> ret(lhs);
	return ret+=rhs;
}

// - operators //

template <typename T>
DenseVectorTemp<T> operator-(const DenseVector<T>& lhs, const DenseVector<T>& rhs)
{
	DenseVectorTemp<T> ret(lhs);
	ret-=rhs;
	return ret;
}

template <typename T>
DenseVectorTemp<T> operator-(const DenseVector<T>& lhs, const DenseVectorView<T>& rhs)
{
	DenseVectorTemp<T> ret(lhs);
	ret-=rhs;
	return ret;
}

template <typename T>
DenseVectorTemp<T> operator-(const DenseVectorView<T>& lhs, const DenseVector<T>& rhs)
{
	DenseVectorTemp<T> ret(lhs);
	ret-=rhs;
	return ret;
}

template <typename T>
DenseVectorTemp<T> operator-(const DenseVectorView<T>& lhs, const DenseVectorView<T>& rhs)
{
	DenseVectorTemp<T> ret(lhs);
	ret-=rhs;
	return ret;
}

template <typename T>
DenseVectorTemp<T> operator-(const DenseVector<T>& lhs, const DenseVectorTemp<T>& rhs)
{
	DenseVectorTemp<T> ret(rhs);
	ret.negate();
	return ret+=lhs;
}

template <typename T>
DenseVectorTemp<T> operator-(const DenseVectorTemp<T>& lhs, const DenseVector<T>& rhs)
{
	DenseVectorTemp<T> ret(lhs);
	return ret-=rhs;
}

template <typename T>
DenseVectorTemp<T> operator-(const DenseVectorView<T>& lhs, const DenseVectorTemp<T>& rhs)
{
	DenseVectorTemp<T> ret(rhs);
	ret.negate();
	return ret+=lhs;
}

template <typename T>
DenseVectorTemp<T> operator-(const DenseVectorTemp<T>& lhs, const DenseVectorView<T>& rhs)
{
	DenseVectorTemp<T> ret(lhs);
	return ret-=rhs;
}

template <typename T>
DenseVectorTemp<T> operator-(const DenseVectorTemp<T>& lhs, const DenseVectorTemp<T>& rhs)
{
	DenseVectorTemp<T> ret(lhs);
	return ret-=rhs;
}


template <typename T>
DenseVectorTemp<T> operator*(const T lhs, const DenseVector<T>& rhs)
{
	DenseVectorTemp<T> ret(rhs);
	return ret*=lhs;
}

template <typename T>
DenseVectorTemp<T> operator*(const DenseVector<T>& lhs, const T rhs)
{
	DenseVectorTemp<T> ret(lhs);
	return ret*=rhs;
}

template <typename T>
DenseVectorTemp<T> operator*(const T lhs, const DenseVectorView<T>& rhs)
{
	DenseVectorTemp<T> ret(rhs);
	return ret*=lhs;
}

template <typename T>
DenseVectorTemp<T> operator*(const DenseVectorView<T>& lhs, const T rhs)
{
	DenseVectorTemp<T> ret(lhs);
	return ret*=rhs;
}

template <typename T>
DenseVectorTemp<T> operator*(const T lhs, const DenseVectorTemp<T>& rhs)
{
	DenseVectorTemp<T> ret(rhs);
	return ret*=lhs;
}

template <typename T>
DenseVectorTemp<T> operator*(const DenseVectorTemp<T>& lhs, const T rhs)
{
	DenseVectorTemp<T> ret(lhs);
	return ret*=rhs;
}

template <typename T>
T operator*(const DenseVector<T>& lhs, const DenseVector<T>& rhs)
{
	assert(lhs.size()==rhs.size());
	typename DenseVector<T>::const_iterator it_lhs=lhs.begin(), lhs_end=lhs.end();
	typename DenseVector<T>::const_iterator it_rhs=rhs.begin();
	T dp=0;
	//while(it_lhs!=lhs_end)
	//	dp+=*it_lhs++ * *it_rhs++;
	for(int i=0; i<lhs_end-it_lhs; ++i)
		dp+=it_lhs[i]*it_rhs[i];
	
	return dp;
}

template <typename T>
DenseVectorTemp<T> operator/(const DenseVector<T>& lhs, const T rhs)
{
	DenseVectorTemp<T> ret(lhs);
	return ret/=rhs;
}

template <typename T>
DenseVectorTemp<T> operator/(const DenseVectorView<T>& lhs, const T rhs)
{
	DenseVectorTemp<T> ret(lhs);
	return ret/=rhs;
}

template <typename T>
DenseVectorTemp<T> operator/(const DenseVectorTemp<T>& lhs, const T rhs)
{
	DenseVectorTemp<T> ret(lhs);
	return ret/=rhs;
}

template <typename V>
typename V::value_type nrm2(const V &x)
{
	typename V::value_type val=0;
	for(int i=0; i<x.size(); ++i)
		val+=std::pow(x(i),2);
	return std::sqrt(val);
}

} /* sparse_lib */ 


#endif /* end of include guard: VVOPS_H */
