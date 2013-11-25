/**
 * @file SparseVector.cpp
 * @brief Part of sparse matrix functionality for iSAM.
 * @author Michael Kaess
 * @author Hordur Johannsson
 * @version $Id: SparseVector.cpp 4694 2011-06-09 03:16:34Z hordurj $
 *
 * Copyright (C) 2009-2013 Massachusetts Institute of Technology.
 * Michael Kaess, Hordur Johannsson, David Rosen,
 * Nicholas Carlevaris-Bianco and John. J. Leonard
 *
 * This file is part of iSAM.
 *
 * iSAM is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.
 *
 * iSAM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with iSAM.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <string>
#include <cstring> // memmove()
#include <iostream>
#include <map>

#include "isam/util.h"

#include "isam/SparseVector.h"

using namespace std;
using namespace Eigen;

namespace isam {


// will significantly influence memory usage of large, but very sparse
// matrices; also influence on execution time if chosen too small (10)
const int INITIAL_ENTRIES = 50;

void SparseVector::_copy_from(const SparseVector& vec) {
  _nnz = vec._nnz;
  _nnz_max = vec._nnz_max;

  _indices = new int[_nnz_max];
  memcpy(_indices, vec._indices, _nnz*sizeof(int));
  _values = new double[_nnz_max];
  memcpy(_values, vec._values, _nnz*sizeof(double));
}

void SparseVector::_dealloc() {
  if (_indices != NULL) {
    delete[] _indices;
    _indices = NULL;
  }
  if (_values != NULL) {
    delete[] _values;
    _values = NULL;
  }
}

int SparseVector::_search(int idx) const {
#if 0
  // naive implementation
  int n = 0;
  while (n<_nnz && (_indices[n] < idx)) {
    n++;
  }
#else
  // binary search implementation
  int n = 0;
  if (_nnz>0) {
    int n0 = 0;
    int n1 = _nnz-1;
    while (n1-n0>1) {
      n = n0+((n1-n0)>>1);
      if (_indices[n] <= idx) {
        n0 = n;
      } else {
        n1 = n;
      }
    }
    if (_indices[n1] < idx) {
      n = n1+1;
    } else if ((_indices[n1] == idx) || (_indices[n0] < idx)) {
      n = n1;
    } else {
      n = n0;
    }
  }
#endif
  return n;
}

void SparseVector::_resize(int new_nnz_max) {
  int* new_indices = new int[new_nnz_max];
  memcpy(new_indices, _indices, _nnz*sizeof(int));
  delete[] _indices;
  _indices = new_indices;
  double* new_values = new double[new_nnz_max];
  memcpy(new_values, _values, _nnz*sizeof(double));
  delete[] _values;
  _values = new_values;
  _nnz_max = new_nnz_max;
}

SparseVector::SparseVector() {
  _nnz = 0;
  _nnz_max = INITIAL_ENTRIES;

  _indices = new int[_nnz_max];
  _values = new double[_nnz_max];
}

SparseVector::SparseVector(const SparseVector& vec) {
  _copy_from(vec);
}

SparseVector::SparseVector(const SparseVector& vec, int num, int first) {
  // first have to figure out how many entries in the given range
  _nnz = 0;
  for (int i=0; i<vec._nnz; i++) {
    int idx = vec._indices[i];
    if (idx>=first && idx<first+num) {
      _nnz++;
    }
  }
  // allocate memory accordingly
  _nnz_max = _nnz;
  _indices = new int[_nnz_max];
  _values = new double[_nnz_max];
  // copy data over, renumber indices!
  int n = 0;
  for (int i=0; i<vec._nnz; i++) {
    int idx = vec._indices[i];
    if (idx>=first && idx<first+num) {
      _indices[n] = vec._indices[i]-first;
      _values[n] = vec._values[i];
      n++;
    }
  }
}

SparseVector::SparseVector(int* indices, double* values, int nnz) {
  _nnz = nnz;
  _nnz_max = nnz;

  _indices = new int[_nnz_max];
  _values = new double[_nnz_max];

  memcpy(_indices, indices, nnz*sizeof(int));
  memcpy(_values, values, nnz*sizeof(double));
}

SparseVector::~SparseVector() {
  _dealloc();
}

const SparseVector& SparseVector::operator= (const SparseVector& vec) {
  // sanity check
  if (this==&vec) {
    return *this;
  }

  // free old stuff
  _dealloc();

  // copy rhs
  _copy_from(vec);

  // return self
  return *this;
}

double SparseVector::operator()(int idx) const {
  int n = 0;
  n = _search(idx);
  if (n<_nnz && _indices[n] == idx) {
    return _values[n];
  } else {
    return 0.;
  }
}

void SparseVector::copy_raw(int* indices, double* values) const {
  memcpy(indices, _indices, _nnz*sizeof(int));
  memcpy(values, _values, _nnz*sizeof(double));
}

bool SparseVector::set(int idx, const double val) {
  VectorXd tmp(1);
  tmp << val;
  return set(idx, tmp);
}

bool SparseVector::set(int idx, const VectorXd& vals) {
  bool created_entry = false;
  int n = 0;
  int c = vals.rows();
  
  // First check if we can append
  if (_nnz > 0 && idx > _indices[_nnz-1]) {
    n = _nnz;
  } else {
    n = _search(idx);  // search closest entry
  }
  
  if (n<_nnz && _indices[n] == idx) {
    // entry exists? then overwrite!
    // BIG ASSUMPTION when writing multiple values:
    // they either all exist or they don't
    for (int i=0;i<c;i++) {
      _values[n+i] = vals(i);
    }
  } else {
    // new entries have to be created
    created_entry = true;
    // enough space left?
    if (_nnz+c > _nnz_max) {
      // automatic resizing, amortized cost O(1)
      int new_nnz_max = _nnz_max*2 + c;
      _resize(new_nnz_max);
    }
    // need to make space in the middle?
    if (n<_nnz) {
      memmove(&(_indices[n+c]), &(_indices[n]), (_nnz-n)*sizeof(int));
      memmove(&(_values[n+c]), &(_values[n]), (_nnz-n)*sizeof(double));
    }
    // insert new entry
    _nnz+=c;
    for (int i=0;i<c;i++) {
      _indices[n+i] = idx+i;
      _values[n+i] = vals(i);
    }
  }

  return created_entry;
}

void SparseVector::remove(int idx) {
  int n = 0;
  while (n<_nnz && (_indices[n] < idx)) {
    n++;
  }
  if (n<_nnz && _indices[n] == idx) {
    if (n!=_nnz-1) {
      // move all following entries
      memmove(&(_indices[n]), &(_indices[n+1]), (_nnz-n-1)*sizeof(int));
      memmove(&(_values[n]), &(_values[n+1]), (_nnz-n-1)*sizeof(double));
    }
    _indices[_nnz-1] = 0; // keep all unused entries at 0
    _nnz--;

    // potentially have to release memory; note that we require a quarter
    // to provide some threshold between reserving more and releasing it again
    int quarter = _nnz_max/4;
    if (INITIAL_ENTRIES < quarter && _nnz < quarter) {
      int new_nnz_max = _nnz_max/2;
      _resize(new_nnz_max);
    }
  }
}

int SparseVector::first() const {
  if (_nnz > 0) {
    return _indices[0];
  } else {
    return -1;
  }
}

int SparseVector::last() const {
  if (_nnz > 0) {
    return _indices[_nnz-1];
  } else {
    return -1;
  }
}

void SparseVector::add_entries(int num, int pos) {
  for (int i=0; i<_nnz; i++) {
    if (_indices[i] >= pos) {
      _indices[i] += num;
    }
  }
}

void SparseVector::print() const {
  cout << "%Vector: nnz=" << _nnz << endl;
  for (int i=0; i<_nnz; i++) {
    cout << _indices[i];
    cout << ": " << _values[i];
    cout << endl;
  }
}

}

