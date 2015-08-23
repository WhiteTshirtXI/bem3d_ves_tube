/* Multi-dimensional array template
 * Note:
 *  number of dimensions <= 4 */

#ifndef MARRAY_H
#define MARRAY_H

#include <cassert>
#include <algorithm>

class IndexRange;
template<class T, int N> class MArray;
template<class T, int N> bool MArray_same_shape(const MArray<T,N> &, const MArray<T,N> &);

typedef MArray<int,1> IntArray1d;
typedef MArray<int,2> IntArray2d;
typedef MArray<int,3> IntArray3d;
typedef MArray<int,4> IntArray4d;

typedef MArray<double,1> DoubArray1d;
typedef MArray<double,2> DoubArray2d;
typedef MArray<double,3> DoubArray3d;
typedef MArray<double,4> DoubArray4d;

class IndexRange 
{
public:
    int _lbound, _ubound;

    IndexRange()  : _lbound(0), _ubound(0)
    { }

    IndexRange(int ubound) : _lbound(0), _ubound(ubound)
    {
	assert(_lbound <= _ubound);
    }

    IndexRange(int lbound, int ubound) : _lbound(lbound), _ubound(ubound)
    { 
	assert(_lbound <= _ubound);
    }

    IndexRange(const IndexRange& r) : _lbound(r._lbound), _ubound(r._ubound)
    { }

    ~IndexRange() { }

    IndexRange& operator=(const IndexRange &r) 
    {
	_lbound = r._lbound;    
	_ubound = r._ubound; 
	return *this;
    }

    int& lbound() { return _lbound; }
    const int& lbound() const { return _lbound; }

    int& ubound() { return _ubound; }
    const int& ubound() const { return _ubound; }
};


template<class T, int N> class MArray
{
public:
    int _ndim, _size, _offset;
    int _lbound[4], _ubound[4], _stride[4]; 
    T *_data, *_base;

    void setMemoryLayout();
    void allocateMemory();

    MArray();
    MArray(IndexRange);
    MArray(IndexRange, IndexRange);
    MArray(IndexRange, IndexRange, IndexRange);
    MArray(IndexRange, IndexRange, IndexRange, IndexRange); 
    MArray(const MArray<T,N> &);

    ~MArray() { if (_size > 0) delete [] _data; }

    // Memory access
    T* data() { return _data; }
    const T* data() const { return _data; }

    // Number of dimension
    int ndim() const { return _ndim; }

    // Index lower bound
    int lbound(int i) const {
	assert(i >= 0 && i < _ndim);
	return _lbound[i];
    }

    // Index upper bound 
    int ubound(int i) const 
    {
	assert(i >= 0 && i < _ndim);
	return _ubound[i];
    }

    // Total array size
    int size() const
    {
	return _size;
    }


    // Size of a dimension
    int size(int i) const
    {
	assert(i >= 0 && i < _ndim);
	return _ubound[i] - _lbound[i];
    }

    // Stride of a dimension
    int stride(int i) const
    {
	assert(i >= 0 && i < _ndim);
	return _stride[i];
    }

    //**********************************************************************
    // Access 1D array element
    T& operator()(int i0)
    {
	assert(_ndim == 1);
	assert(i0 >= _lbound[0] && i0 < _ubound[0]);
	return *(_base + i0);
    }

    const T& operator()(int i0) const
    {
	assert(_ndim == 1);
	assert(i0 >= _lbound[0] && i0 < _ubound[0]);
	return *(_base + i0);
    }


    //**********************************************************************
    // Access 2D array element
    T& operator()(int i0, int i1)
    {
	assert(_ndim == 2);
	assert(i0 >= _lbound[0] && i0 < _ubound[0]);
	assert(i1 >= _lbound[1] && i1 < _ubound[1]);
	return *(_base + i0*_stride[0] + i1);
    }

    const T& operator()(int i0, int i1) const
    {
	assert(_ndim == 2);
	assert(i0 >= _lbound[0] && i0 < _ubound[0]);
	assert(i1 >= _lbound[1] && i1 < _ubound[1]);
	return *(_base + i0*_stride[0] + i1);
    }


    //**********************************************************************
    // Access 3D array element
    T& operator()(int i0, int i1, int i2)
    {
	assert(_ndim == 3);
	assert(i0 >= _lbound[0] && i0 < _ubound[0]);
	assert(i1 >= _lbound[1] && i1 < _ubound[1]);
	assert(i2 >= _lbound[2] && i2 < _ubound[2]);

	return *(_base + i0*_stride[0] + i1*_stride[1] + i2);
    }

    const T& operator()(int i0, int i1, int i2) const
    {
	assert(_ndim == 3);
	assert(i0 >= _lbound[0] && i0 < _ubound[0]);
	assert(i1 >= _lbound[1] && i1 < _ubound[1]);
	assert(i2 >= _lbound[2] && i2 < _ubound[2]);

	return *(_base + i0*_stride[0] + i1*_stride[1] + i2);
    }


    //**********************************************************************
    // Access 4D array element
    T& operator()(int i0, int i1, int i2, int i3)
    {
	assert(_ndim == 4);
	assert(i0 >= _lbound[0] && i0 < _ubound[0]);
	assert(i1 >= _lbound[1] && i1 < _ubound[1]);
	assert(i2 >= _lbound[2] && i2 < _ubound[2]);
	assert(i3 >= _lbound[3] && i3 < _ubound[3]);
	return *(_base + i0*_stride[0] + i1*_stride[1] + i2*_stride[2] + i3);
    }
    
    const T& operator()(int i0, int i1, int i2, int i3) const
    {
	assert(_ndim == 4);
	assert(i0 >= _lbound[0] && i0 < _ubound[0]);
	assert(i1 >= _lbound[1] && i1 < _ubound[1]);
	assert(i2 >= _lbound[2] && i2 < _ubound[2]);
	assert(i3 >= _lbound[3] && i3 < _ubound[3]);
	return *(_base + i0*_stride[0] + i1*_stride[1] + i2*_stride[2] + i3);
    }

    //**********************************************************************
    // Resize
    MArray<T,N>& resize(IndexRange);
    MArray<T,N>& resize(IndexRange, IndexRange);
    MArray<T,N>& resize(IndexRange, IndexRange, IndexRange);
    MArray<T,N>& resize(IndexRange, IndexRange, IndexRange, IndexRange);

    MArray<T,N>& operator=(T);
    MArray<T,N>& operator=(const MArray<T,N>&);
    MArray<T,N>& operator+=(T);
    MArray<T,N>& operator+=(const MArray<T,N>&);
    MArray<T,N>& operator-=(T);
    MArray<T,N>& operator-=(const MArray<T,N>&);
    MArray<T,N>& operator*=(T);
};


// Compute memory layout including offset, stride and size
template<class T, int N>
void MArray<T,N>::setMemoryLayout() 
{
    if (_ndim == 0) {
        _size = 0;
	return;
    }

    _stride[_ndim-1] = 1;
    _offset = _lbound[_ndim-1];

    for (int d = _ndim-2; d >= 0; d--)  {
        _stride[d] = _stride[d+1] * (_ubound[d+1] - _lbound[d+1]);
	_offset += _stride[d]*_lbound[d];
    }

    _size = _stride[0]*(_ubound[0] - _lbound[0]);
}


// Allocate memory
template<class T, int N>
void MArray<T,N>::allocateMemory() 
{
    if (_size > 0) {
        _data = new T[_size];
        _base = _data - _offset;
    } else {
        _data = _base = NULL;
    }
}


// Default constructor
template<class T, int N>
MArray<T,N>::MArray() 
{
    assert(N > 0 && N <= 4);

    _ndim = N;
    for (int i = 0; i < _ndim; i++) {
        _ubound[i] = _lbound[i] = 0;
    }

    setMemoryLayout();
    allocateMemory();
}


// Cosntructor for 1-d array
template<class T, int N>
MArray<T,N>::MArray(IndexRange r0) 
{
    assert(N == 1);

    _ndim = N;

    _lbound[0] = r0.lbound();
    _ubound[0] = r0.ubound();

    setMemoryLayout();
    allocateMemory();
}


// Constructor for 2-d array
template<class T, int N>
MArray<T,N>::MArray(IndexRange r0, IndexRange r1) 
{
    assert(N == 2);

    _ndim = N;

    _lbound[0] = r0.lbound();
    _ubound[0] = r0.ubound();

    _lbound[1] = r1.lbound();
    _ubound[1] = r1.ubound();

    setMemoryLayout();
    allocateMemory();
}


// Constructor for 3-d array
template<class T, int N>
MArray<T,N>::MArray(IndexRange r0, IndexRange r1, IndexRange r2) 
{
    assert(N == 3);

    _ndim = N;

    _lbound[0] = r0.lbound();
    _ubound[0] = r0.ubound();

    _lbound[1] = r1.lbound();
    _ubound[1] = r1.ubound();

    _lbound[2] = r2.lbound();
    _ubound[2] = r2.ubound();

    setMemoryLayout();
    allocateMemory();
}


// Constructor for 4-d array
template<class T, int N>
MArray<T,N>::MArray(IndexRange r0, IndexRange r1, IndexRange r2, IndexRange r3) 
{
    assert(N == 4);

    _ndim = N;

    _lbound[0] = r0.lbound();
    _ubound[0] = r0.ubound();

    _lbound[1] = r1.lbound();
    _ubound[1] = r1.ubound();

    _lbound[2] = r2.lbound();
    _ubound[2] = r2.ubound();

    _lbound[3] = r3.lbound();
    _ubound[3] = r3.ubound();

    setMemoryLayout();
    allocateMemory();
}


// Copy constructor
template<class T, int N>
MArray<T,N>::MArray(const MArray<T,N> &x)
{
    _ndim = x._ndim;

    for (int l = 0; l < _ndim; ++l) {
        _lbound[l] = x._lbound[l];
	_ubound[l] = x._ubound[l];
	_stride[l] = x._stride[l];
    }
    _size = x._size;
    _offset = x._offset;

    allocateMemory();
    std::copy(x._data, x._data + x._size, _data);
}


// Resize
template<class T, int N>
MArray<T,N>& MArray<T,N>::resize(IndexRange r0)
{
    assert(_ndim == 1);

    if (_lbound[0] == r0.lbound() && _ubound[0] == r0.ubound()) return *this;

    _lbound[0] = r0.lbound();
    _ubound[0] = r0.ubound();

    if (_data) delete []_data;
    setMemoryLayout();
    allocateMemory();

    return *this;
}


template<class T, int N>
MArray<T,N>& MArray<T,N>::resize(IndexRange r0, IndexRange r1)
{
    assert(_ndim == 2);

    if (_lbound[0] == r0.lbound() && _ubound[0] == r0.ubound() &&
        _lbound[1] == r1.lbound() && _ubound[1] == r1.ubound()) return *this;

    _lbound[0] = r0.lbound();
    _ubound[0] = r0.ubound();

    _lbound[1] = r1.lbound();
    _ubound[1] = r1.ubound();

    if (_data) delete []_data;
    setMemoryLayout();
    allocateMemory();

    return *this;
}

template<class T, int N>
MArray<T,N>& MArray<T,N>::resize(IndexRange r0, IndexRange r1, IndexRange r2)
{
    assert(_ndim == 3);

    if (_lbound[0] == r0.lbound() && _ubound[0] == r0.ubound() &&
        _lbound[1] == r1.lbound() && _ubound[1] == r1.ubound() &&
        _lbound[2] == r2.lbound() && _ubound[2] == r2.ubound()) return *this;

    _lbound[0] = r0.lbound();
    _ubound[0] = r0.ubound();

    _lbound[1] = r1.lbound();
    _ubound[1] = r1.ubound();

    _lbound[2] = r2.lbound();
    _ubound[2] = r2.ubound();

    if (_data) delete []_data;
    setMemoryLayout();
    allocateMemory();

    return *this;
}


template<class T, int N>
MArray<T,N>& MArray<T,N>::resize(IndexRange r0, IndexRange r1, IndexRange r2, IndexRange r3)
{
    assert(_ndim == 4);

    if (_lbound[0] == r0.lbound() && _ubound[0] == r0.ubound() &&
        _lbound[1] == r1.lbound() && _ubound[1] == r1.ubound() &&
        _lbound[2] == r2.lbound() && _ubound[2] == r2.ubound() &&
        _lbound[3] == r3.lbound() && _ubound[3] == r3.ubound()) return *this;

    _lbound[0] = r0.lbound();
    _ubound[0] = r0.ubound();

    _lbound[1] = r1.lbound();
    _ubound[1] = r1.ubound();

    _lbound[2] = r2.lbound();
    _ubound[2] = r2.ubound();

    _lbound[3] = r3.lbound();
    _ubound[3] = r3.ubound();

    if (_data) delete []_data;
    setMemoryLayout();
    allocateMemory();

    return *this;
}


template<class T, int N>
MArray<T,N>& MArray<T,N>::operator=(T x)
{
    std::fill_n(_data, _size, x);
    return *this;
}

template<class T, int N>
MArray<T,N>& MArray<T,N>::operator=(const MArray<T,N> &x)
{
    assert( MArray_same_shape(*this, x) );

    std::copy(x._data, x._data + x._size, _data);
    return *this;
}


template<class T, int N>
MArray<T,N>& MArray<T,N>::operator+=(T x)
{
    T *p = _data;
    for (int i = 0; i < _size; i++) *(p++) += x;
    return *this;
}

template<class T, int N>
MArray<T,N>& MArray<T,N>::operator+=(const MArray<T,N> &x)
{
    assert( MArray_same_shape(*this, x) );

    T *p = _data, *q = x._data;
    for (int i = 0; i < _size; ++i) *(p++) += *(q++);
    return *this;
}

template<class T, int N>
MArray<T,N>& MArray<T,N>::operator-=(T x)
{
    T *p = _data;
    for (int i = 0; i < _size; i++) *(p++) -= x;
    return *this;
}


template<class T, int N>
MArray<T,N>& MArray<T,N>::operator-=(const MArray<T,N> &x)
{
    assert( MArray_same_shape(*this, x) );

    T *p = _data, *q = x._data;
    for (int i = 0; i < _size; ++i) *(p++) -= *(q++);
    return *this;
}


template<class T, int N>
MArray<T,N>& MArray<T,N>::operator*=(T a)
{
    T *p = _data;
    for (int i = 0; i < _size; i++) *(p++) *= a;
    return *this;
}


/* Determine if the two arrays have the same shape
 * Criterior:
 *   -- same number of dimensions 
 *   -- same size in each dimension, but not necessarily
 *      the same lower and upper bounds */
template<class T, int N>
bool MArray_same_shape(const MArray<T,N> &x, const MArray<T,N> &y)
{
    if (x._ndim != y._ndim) return false;

    for (int d = 0; d < x._ndim; d++) {
        if (x.size(d) != y.size(d)) return false;
    }

    return true;
}

#endif
