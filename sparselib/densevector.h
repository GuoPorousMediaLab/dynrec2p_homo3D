#ifndef SPARSE_LIB_DENSEVECTOR_H
#define SPARSE_LIB_DENSEVECTOR_H

#include <algorithm>
#include <iostream>
#include <memory>
#include <functional>
#include <assert.h>

namespace sparse_lib
{

template <typename T>
class DenseVectorView;

template <typename T>
class DenseVectorTemp;

template <typename T>
class DenseVector
{

friend class DenseVectorView<T>;
friend class DenseVectorTemp<T>;
	
public:
	typedef T* iterator;
	typedef const T* const_iterator;
	typedef size_t size_type;
	typedef T value_type;
	DenseVector (): reference_ptr_(0) { create_(); }
	explicit DenseVector (size_type n, const T& fill=T()): reference_ptr_(0) { create_(n, fill); }
	DenseVector (const DenseVector& v): reference_ptr_(0) { create_(v.begin(), v.end()); }
//	DenseVector (const DenseVectorTemp<T>& v): data_(v.data_), limit_(v.limit_), reference_ptr_(v.reference_ptr_) { ++*reference_ptr_; }
//	DenseVector (const DenseVectorView<T>& v): reference_ptr_(0) { create_(v.begin(), v.end()); }
	~DenseVector () { uncreate_(); }
	
	// Assignment Operators
	DenseVector& operator=(const DenseVector&);
	DenseVector& operator=(const DenseVectorTemp<T>&);
	DenseVector& operator=(const DenseVectorView<T>&);
	DenseVector& operator=(const T v);
	
	DenseVector& operator+=(const T v);
	DenseVector& operator-=(const T v);
	DenseVector& operator*=(const T v);
	DenseVector& operator/=(const T v);
	DenseVector& operator-();
	
	template <typename V>
	DenseVector& operator+=(const V& v);
	
	template <typename V>
	DenseVector& operator-=(const V& v);
	
	template <typename V>
	DenseVector& operator*=(const V& v);
	
	// Change dimension
	void resize(size_type new_size) { grow_(new_size); }

	// Information
	size_type size() const { return limit_ - data_; }
	
	// Element Access
	T& operator()(size_type i) { return data_[i]; }
	const T& operator()(size_type i) const { return data_[i]; }
	DenseVectorView<T> range(size_type first, size_type last) { return DenseVectorView<T>(*this, first, last); }
	const DenseVectorView<T> range(size_type first, size_type last) const { return DenseVectorView<T>(*this, first, last); }
	
	// Iterators
	iterator begin() { return data_; }
	const_iterator begin() const { return data_; }
	iterator end() { return limit_; }
	const_iterator end() const { return limit_; }
	
	iterator data() { return data_; }
	const_iterator data() const { return data_; }
	
private:
	iterator data_;
	iterator limit_;
	
	std::allocator<T> alloc_;
	
	// For "GC"
	size_type *reference_ptr_;
	
	void create_();
	void create_(size_type, const T&);
	void create_(const_iterator, const_iterator);
	
	void uncreate_();
	
	void grow_(size_type);
	
};

template <typename T>
std::ostream& operator<<(std::ostream& o, const DenseVector<T>& v);

template <typename T>
class DenseVectorView
{
friend class DenseVector<T>;
friend class DenseVectorTemp<T>;

public:
	typedef typename DenseVector<T>::iterator iterator;
	typedef typename DenseVector<T>::const_iterator const_iterator;
	typedef typename DenseVector<T>::size_type size_type;
	typedef typename DenseVector<T>::value_type value_type;
	DenseVectorView() { set_(); }
	template <typename V>
	explicit DenseVectorView(const V& v, size_type first, size_type last) { set_(v, first, last); }
	DenseVectorView(const DenseVector<T>& v) { set_(v, 0, v.size()-1); }
	DenseVectorView(const DenseVectorTemp<T>& v) { set_(v, 0, v.size()-1); }
	DenseVectorView(const DenseVectorView& v): data_(v.data_), limit_(v.limit_), start_(v.start_), finish_(v.finish_) { }
	
	DenseVectorView& operator=(const DenseVectorView&);
	DenseVectorView& operator=(const DenseVector<T>&);
	DenseVectorView& operator=(const T);
	
	DenseVectorView& operator+=(const T v);
	DenseVectorView& operator-=(const T v);
	DenseVectorView& operator*=(const T v);
	DenseVectorView& operator/=(const T v);

	template <typename V>
	DenseVectorView& operator+=(const V& v);
	
	template <typename V>
	DenseVectorView& operator-=(const V& v);
	
	template <typename V>
	DenseVectorView& operator*=(const V& v);
	
	size_type size() const { return limit_ - data_; }
	
	// Element Access
	T& operator()(size_type i) { return data_[i]; }
	const T& operator()(size_type i) const { return data_[i]; }
	DenseVectorView<T> range(size_type first, size_type last) { return DenseVectorView<T>(*this, first, last); }
	const DenseVectorView<T> range(size_type first, size_type last) const { return DenseVectorView<T>(*this, first, last); }
	
	iterator data() { return data_; }
	const_iterator data() const { return data_; }
	
	// Iterators
	iterator begin() { return data_; }
	const_iterator begin() const { return data_; }
	iterator end() { return limit_; }
	const_iterator end() const { return limit_; }
	
	void print() const;
private:
	void set_() { data_ = limit_ = 0; start_ = finish_ = 0; }
	void set_(const DenseVector<T>& v, size_type first, size_type last);
	void set_(const DenseVectorTemp<T>& v, size_type first, size_type last);
	void replace_(const_iterator first, const_iterator last);
	
	iterator data_;
	iterator limit_;
	size_type start_, finish_;
};

template <typename T>
std::ostream& operator<<(std::ostream& o, const DenseVectorView<T>& v);

template <typename T>
class DenseVectorTemp
{
	
friend class DenseVector<T>;	
friend class DenseVectorView<T>;

public:
	typedef typename DenseVector<T>::iterator iterator;
	typedef typename DenseVector<T>::const_iterator const_iterator;
	typedef typename DenseVector<T>::size_type size_type;
	typedef typename DenseVector<T>::value_type value_type;
	
	DenseVectorTemp(): reference_ptr_(new size_type(1)), data_(0), limit_(0) { }
	explicit DenseVectorTemp (size_type n, const T& fill=T()): reference_ptr_(new size_type(1)) { create_(n, fill); }
	DenseVectorTemp(const DenseVector<T>& v): reference_ptr_(new size_type(1)) { create_(v.begin(), v.end()); }
	DenseVectorTemp(const DenseVectorView<T>& v): reference_ptr_(new size_type(1)) { create_(v.begin(), v.end()); }
	DenseVectorTemp(const DenseVectorTemp<T>& v): data_(v.data_), limit_(v.limit_), reference_ptr_(v.reference_ptr_) { ++(*reference_ptr_); }
	
	~DenseVectorTemp() { uncreate_(); }
	
	DenseVectorTemp& operator=(const DenseVectorTemp&);
	DenseVectorTemp& operator=(const DenseVectorView<T>&);
	DenseVectorTemp& operator=(const DenseVector<T>&);
	
	DenseVectorTemp& operator+=(const T v);
	DenseVectorTemp& operator-=(const T v);
	DenseVectorTemp& operator*=(const T v);
	DenseVectorTemp& operator/=(const T v);
	void negate();
	
	template <typename V>
	DenseVectorTemp& operator+=(const V& v);
	
	template <typename V>
	DenseVectorTemp& operator-=(const V& v);
	
	template <typename V>
	DenseVectorTemp& operator*=(const V& v);
	
	size_type size() const { return limit_ - data_; }
	
	// Element Access
	T& operator()(size_type i) { return data_[i]; }
	const T& operator()(size_type i) const { return data_[i]; }
	DenseVectorView<T> range(size_type first, size_type last) { return DenseVectorView<T>(*this, first, last); }
	const DenseVectorView<T> range(size_type first, size_type last) const { return DenseVectorView<T>(*this, first, last); }
	
	iterator data() { return data_; }
	const_iterator data() const { return data_; }
	
	// Iterators
	iterator begin() { return data_; }
	const_iterator begin() const { return data_; }
	iterator end() { return limit_; }
	const_iterator end() const { return limit_; }
private:
	iterator data_;
	iterator limit_;
	size_type *reference_ptr_;
	
	std::allocator<T> alloc_;
	
	void create_(const_iterator, const_iterator);
	void create_(size_type, const T&);
	void uncreate_();
	
};



// --------------- DenseVector ---------------//

template <typename T>
void DenseVector<T>::create_()
{
	data_ = limit_ = 0;
}

template <typename T>
void DenseVector<T>::create_(size_type n, const T& val)
{

	data_ = alloc_.allocate(n);
	limit_ = data_ + n;
	std::uninitialized_fill(data_, limit_, val);
}

template <typename T>
void DenseVector<T>::create_(const_iterator i, const_iterator j)
{
	data_ = alloc_.allocate(j - i);
	limit_ = std::uninitialized_copy(i, j, data_);
}

template <typename T>
void DenseVector<T>::uncreate_()
{
	if (data_) {
		// destroy (in reverse order) the elements that were constructed
		iterator it = limit_;
		while (it != data_)
			alloc_.destroy(--it);

		// return all the space that was allocated
		alloc_.deallocate(data_, limit_ - data_);
	}
	// reset pointers to indicate that the `Vec' is empty again
	data_ = limit_ = 0;
	if(reference_ptr_ && --(*reference_ptr_) == 0)
		delete reference_ptr_;
	reference_ptr_ = 0;

}

template <class T>
void DenseVector<T>::grow_(size_type size)
{

	// allocate new space and copy existing elements to the new space

	iterator new_data = alloc_.allocate(size);

	// [BUG?] - Is the below correct?  Or should it be limit_ - data_ - 1?
	iterator new_avail = std::uninitialized_copy(data_, data_ + std::min(size_type(limit_ - data_), size), new_data);

	// return the old space
	uncreate_();

	// reset pointers to point to the newly allocated space
	data_ = new_data;
	limit_ = data_ + size;
}

template <typename T>
DenseVector<T>& DenseVector<T>::operator=(const DenseVector& rhs)
{
	if(&rhs != this) {
		uncreate_();
		
		create_(rhs.begin(), rhs.end());
	}
	return *this;
}

template <typename T>
DenseVector<T>& DenseVector<T>::operator=(const DenseVectorTemp<T>& rhs)
{
	uncreate_();
	data_=rhs.data_;
	limit_=rhs.limit_;
	reference_ptr_=rhs.reference_ptr_;
	++*reference_ptr_;
	return *this;
}

template <typename T>
DenseVector<T>& DenseVector<T>::operator=(const DenseVectorView<T>& rhs)
{
	uncreate_();
	create_(rhs.begin(), rhs.end());
	
	return *this;
}

template <typename T>
DenseVector<T>& DenseVector<T>::operator=(const T rhs)
{

	std::fill(data_, limit_, rhs);
	return *this;
}

template <typename T>
DenseVector<T>& DenseVector<T>::operator+=(const T rhs)
{
	std::transform(begin(), end(), begin(), std::bind2nd(std::plus<T>(), rhs));
	return *this;
}

template <typename T>
DenseVector<T>& DenseVector<T>::operator-=(const T rhs)
{
	std::transform(begin(), end(), begin(), std::bind2nd(std::minus<T>(), rhs));
	return *this;
}

template <typename T>
DenseVector<T>& DenseVector<T>::operator*=(const T rhs)
{
	std::transform(begin(), end(), begin(), std::bind2nd(std::multiplies<T>(), rhs));
	return *this;
}

template <typename T>
DenseVector<T>& DenseVector<T>::operator/=(const T rhs)
{
	std::transform(begin(), end(), begin(), std::bind2nd(std::divides<T>(), rhs));
	return *this;
}

template <typename T> template <typename V>
DenseVector<T>& DenseVector<T>::operator+=(const V& rhs)
{
	assert(size()==rhs.size());
	std::transform(begin(), end(), rhs.begin(), begin(), std::plus<T>());
	return *this;
}

template <typename T> template <typename V>
DenseVector<T>& DenseVector<T>::operator-=(const V& rhs)
{
	assert(size()==rhs.size());
	std::transform(begin(), end(), rhs.begin(), begin(), std::minus<T>());
	return *this;
}

template <typename T> template <typename V>
DenseVector<T>& DenseVector<T>::operator*=(const V& rhs)
{
	assert(size()==rhs.size());
	std::transform(begin(), end(), rhs.begin(), begin(), std::multiplies<T>());
	return *this;
}

template <typename T>
std::ostream& operator<<(std::ostream& o, const DenseVector<T>& v) {
	typename DenseVector<T>::const_iterator i = v.begin(), end=v.end();
	o.precision(8);
    o.setf(std::ios::fixed);
	o << "Vector with " << v.size() << " elements:\n[  ";
	while(i!=end) {
		o << *(i++) << "  ";
	}
	std::cout << "]";
	o << std::endl;
    return o;
}

// --------------- DenseVectorView --------------- //

template <typename T>
DenseVectorView<T>& DenseVectorView<T>::operator=(const DenseVectorView<T>& rhs)
{
	assert(size() == rhs.size());
	if(&rhs != this && rhs.begin()!=begin()) {
		replace_(rhs.begin(), rhs.end());
	}
	return *this;
}

template <typename T>
DenseVectorView<T>& DenseVectorView<T>::operator=(const DenseVector<T>& rhs)
{
	assert(size() == rhs.size());
	replace_(rhs.begin(), rhs.end());
	return *this;
}

template <typename T>
DenseVectorView<T>& DenseVectorView<T>::operator=(const T rhs)
{
	std::fill(data_, limit_, rhs);
	return *this;
}


template <typename T>
void DenseVectorView<T>::set_(const DenseVector<T>& v, size_type first, size_type last)
{
	data_=v.data_+first;
	limit_=v.data_+last+1;
	start_=first;
	finish_=last;
}

template <typename T>
void DenseVectorView<T>::set_(const DenseVectorTemp<T>& v, size_type first, size_type last)
{
	data_=v.data_+first;
	limit_=v.data_+last+1;
	start_=first;
	finish_=last;
}

template <typename T>
inline void DenseVectorView<T>::replace_(const_iterator first, const_iterator last)
{
	iterator it=begin();
	std::copy(first, last, it);
}

template <typename T>
DenseVectorView<T>& DenseVectorView<T>::operator+=(const T rhs)
{
	std::transform(begin(), end(), begin(), std::bind2nd(std::plus<T>(), rhs));
	return *this;
}

template <typename T>
DenseVectorView<T>& DenseVectorView<T>::operator-=(const T rhs)
{
	std::transform(begin(), end(), begin(), std::bind2nd(std::minus<T>(), rhs));
	return *this;
}

template <typename T>
DenseVectorView<T>& DenseVectorView<T>::operator*=(const T rhs)
{
	std::transform(begin(), end(), begin(), std::bind2nd(std::multiplies<T>(), rhs));
	return *this;
}

template <typename T>
DenseVectorView<T>& DenseVectorView<T>::operator/=(const T rhs)
{
	std::transform(begin(), end(), begin(), std::bind2nd(std::divides<T>(), rhs));
	return *this;
}

template <typename T> template <typename V>
DenseVectorView<T>& DenseVectorView<T>::operator+=(const V& rhs)
{
	assert(size()==rhs.size());
	std::transform(begin(), end(), rhs.begin(), begin(), std::plus<T>());
	return *this;
}

template <typename T> template <typename V>
DenseVectorView<T>& DenseVectorView<T>::operator-=(const V& rhs)
{
	assert(size()==rhs.size());
	std::transform(begin(), end(), rhs.begin(), begin(), std::minus<T>());
	return *this;
}

template <typename T> template <typename V>
DenseVectorView<T>& DenseVectorView<T>::operator*=(const V& rhs)
{
	assert(size()==rhs.size());
	std::transform(begin(), end(), rhs.begin(), begin(), std::multiplies<T>());
	return *this;
}

template <typename T>
void DenseVectorView<T>::print() const { std::cout << *this << std::endl; }

template <typename T>
std::ostream& operator<<(std::ostream& o, const DenseVectorView<T>& v) {
	typename DenseVector<T>::const_iterator i = v.begin(), end=v.end();
	o << "[" << v.size() << "]:\n";
	while(i!=end)
		o << "  " << *(i++) << "\n";
    return o;
}


// --------------- DenseVectorTemp --------------- //

template <typename T>
DenseVectorTemp<T>& DenseVectorTemp<T>::operator=(const DenseVectorTemp& rhs)
{
	++*rhs.reference_ptr_;
	uncreate_();
	data_=rhs.data_;
	limit_=rhs.limit_;
	reference_ptr_=rhs.reference_ptr_;
	return *this;
}

template <typename T>
DenseVectorTemp<T>& DenseVectorTemp<T>::operator=(const DenseVectorView<T>& rhs)
{
	uncreate_();
	create_(rhs.begin(), rhs.end());
	reference_ptr_ = new size_type(1);
	return *this;
}

template <typename T>
DenseVectorTemp<T>& DenseVectorTemp<T>::operator=(const DenseVector<T>& rhs)
{
	uncreate_();
	create_(rhs.begin(), rhs.end());
	reference_ptr_ = new size_type(1);
	return *this;
}

template <typename T>
void DenseVectorTemp<T>::create_(size_type n, const T& val)
{

	data_ = alloc_.allocate(n);
	limit_ = data_ + n;
	std::uninitialized_fill(data_, limit_, val);
}

template <typename T>
void DenseVectorTemp<T>::create_(const_iterator i, const_iterator j)
{
	data_ = alloc_.allocate(j - i);
	limit_ = std::uninitialized_copy(i, j, data_);
}

template <typename T>
void DenseVectorTemp<T>::uncreate_()
{
	if(reference_ptr_ && --*reference_ptr_ == 0) {
		if (data_) {
			// destroy (in reverse order) the elements that were constructed
			iterator it = limit_;
			while (it != data_)
				alloc_.destroy(--it);

			// return all the space that was allocated
			alloc_.deallocate(data_, limit_ - data_);
		}
		delete reference_ptr_;
	}
		reference_ptr_ = 0;
		// reset pointers to indicate that the `Vec' is empty again
		data_ = limit_ = 0;
}

template <typename T>
DenseVectorTemp<T>& DenseVectorTemp<T>::operator+=(const T rhs)
{
	std::transform(begin(), end(), begin(), std::bind2nd(std::plus<T>(), rhs));
	return *this;
}

template <typename T>
DenseVectorTemp<T>& DenseVectorTemp<T>::operator-=(const T rhs)
{
	std::transform(begin(), end(), begin(), std::bind2nd(std::minus<T>(), rhs));
	return *this;
}

template <typename T>
DenseVectorTemp<T>& DenseVectorTemp<T>::operator*=(const T rhs)
{
	std::transform(begin(), end(), begin(), std::bind2nd(std::multiplies<T>(), rhs));
	return *this;
}

template <typename T>
DenseVectorTemp<T>& DenseVectorTemp<T>::operator/=(const T rhs)
{
	std::transform(begin(), end(), begin(), std::bind2nd(std::divides<T>(), rhs));
	return *this;
}

template <typename T> template <typename V>
DenseVectorTemp<T>& DenseVectorTemp<T>::operator+=(const V& rhs)
{
	assert(size()==rhs.size());
	std::transform(begin(), end(), rhs.begin(), begin(), std::plus<T>());
	return *this;
}

template <typename T> template <typename V>
DenseVectorTemp<T>& DenseVectorTemp<T>::operator-=(const V& rhs)
{
	assert(size()==rhs.size());
	std::transform(begin(), end(), rhs.begin(), begin(), std::minus<T>());
	return *this;
}

template <typename T> template <typename V>
DenseVectorTemp<T>& DenseVectorTemp<T>::operator*=(const V& rhs)
{
	assert(size()==rhs.size());
	std::transform(begin(), end(), rhs.begin(), begin(), std::multiplies<T>());
	return *this;
}

template <typename T>
void DenseVectorTemp<T>::negate()
{
	std::transform(begin(), end(), begin(), std::negate<T>());
}


} /* sparse_lib */ 

#endif /* end of include guard: SPARSE_LIB_DENSEVECTOR_H */
