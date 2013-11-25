/**
 * @file SparseVector.cpp
 * @brief part of sparse matrix functionality for iSAM
 * @author Michael Kaess
 * @version $Id: SparseVector.h 3970 2011-02-16 00:16:54Z kaess $
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

#pragma once

#include <Eigen/Dense>

#include "isam/util.h"

namespace isam {

class SparseVector {
  int _nnz;
  int _nnz_max;
  int* _indices;
  double* _values;

  /**
  * Copy data from one sparse vector to a new one - private
  * @param vec Existing sparse vector to copy from.
  */
  void _copy_from(const SparseVector& vec);

  /**
  * Deallocate memory - private.
  */
  void _dealloc();

  /**
   * Search for an index - private;
   * @return Index of entry, or if not in list, index of preceding entry + 1.
   */
  int _search(int idx) const;

  /**
   * Resize the sparse vector - dangerous, only used internally - private.
   * @param new_nnz_max New size of vector.
   */
  void _resize(int new_nnz_max);

public:
  /**
   * Standard constructor.
   */
  SparseVector();

  SparseVector(int nnz_max) : _nnz(0), _nnz_max(nnz_max) {
    _indices = new int[_nnz_max];
    _values = new double[_nnz_max];
  }
  
  /**
   * Copy constructor - note that overwriting operator= is also necessary!
   * @param vec SparseVector to initialize from.
   */
  SparseVector(const SparseVector& vec);

  /**
   * Subvector copy constructor.
   * @param vec SparseVector to copy from.
   * @param num Number of entries to copy (note sparse - most will not exist!)
   * @param first Index of first entry of new vector.
   */
  SparseVector(const SparseVector& vec, int num, int first = 0);

  /**
   * Construct from raw data.
   * @param indices Array of row indices.
   * @param values Array of values.
   * @param nnz Length of vector.
   */
  SparseVector(int* indices, double* values, int nnz);

  /**
   * Destructor.
   */
  virtual ~SparseVector();

  /**
   * Assignment operator (see also copy constructor).
   * @param vec Right-hand-side vector in assignment
   * @return self.
   */
  const SparseVector& operator= (const SparseVector& vec);

  /**
   * Access value of an entry.
   * @param idx Index of entry to access.
   * @return Value of entry.
   */
  double operator()(int idx) const;

  /**
   * Copy raw data, useful for converting to Matlab format.
   * @param indices Destination for row indices.
   * @param values Destination for values.
   */
  void copy_raw(int* indices, double* values) const;

  /**
   * Create a new entry at the end of the sparse vector.
   * Used for efficient incremental creation of SparseVector in apply_givens().
   * @param idx Index of entry to add, has to be greater than last_idx().
   * @param val Value of new entry.
   */
  inline void append(int idx, const double val = 0.) {
    requireDebug(_nnz==0 || _indices[_nnz-1] < idx, "SparseVector::append: index has to be after last entry");

    if (_nnz+1 > _nnz_max) {
      // automatic resizing, amortized cost O(1)
      int new_nnz_max = _nnz_max*2;
      _resize(new_nnz_max);
    }
    // insert new entry
    _indices[_nnz] = idx;
    _values[_nnz] = val;
    _nnz++;
  }

  /**
   * Assign a value to a specific entry.
   * @param idx Index of entry to assign to.
   * @param val Value to assign.
   * @return True if new entry had to be created (needed to update row index)
   */
  bool set(int idx, const double val = 0.);
  
  
  /**
   * Assign a value to a specific entry.
   * @param idx Index of entry to assign to.
   * @param vals Vector of values to assign.
   * @param c Number of values in double array val.
   * @return True if new entries had to be created
   */
  bool set(int idx, const Eigen::VectorXd& vals);

  /**
   * Remove an entry; simply ignores non-existing entries.
   * @param idx Index of entry to remove.
   */
  void remove(int idx);

  /**
   * Find index of first non-zero entry.
   * @return Index of first non-zero entry, or -1 if vector is empty.
   */
  int first() const;

  /**
   * Find index of last non-zero entry
   * @return Index of last non-zero entry, or -1 if vector is empty.
   */
  int last() const;

  /**
   * Shift indices starting at pos by offset.
   * @param num Number of new entries.
   * @param pos First index that should be shifted.
   */
  void add_entries(int num, int pos);

  /**
   * Prints contents of vector as a list of tuples.
   */
  void print() const;

  inline int nnz() const {
    return _nnz;
  }

  friend class SparseVectorIter;
};

class SparseVectorIter {
  const SparseVector& s;
  int index;
public:
  /**
   * Iterator for SparseVector.
   * @param sv SparseVector.
   */
  inline SparseVectorIter(const SparseVector& sv) : s(sv), index(0) {}

  /**
   * Check if current element valid, ie. if we have not reached the end yet.
   * @return True if we still have a valid entry.
   */
  inline bool valid() const {
    return (index < s._nnz);
  }

  /**
   * Get current element index.
   * @return Current element index.
   */
  inline int get() const {
    requireDebug(index < s._nnz, "SparseVectorIter::get(): Index out of range.");
    int tmp = s._indices[index];
    return tmp;
  }

  /**
   * Get current element index and value.
   * @param val Current element value returned.
   * @return Current element index.
   */
  inline int get(double& val) const {
    requireDebug(index < s._nnz, "SparseVectorIter::get(): Index out of range.");
    val = s._values[index];
    return s._indices[index];
  }

  /**
   * Get current element value.
   * @return Current element value returned.
   */
  inline double get_val() const {
    requireDebug(index < s._nnz, "SparseVectorIter::get_val(): Index out of range.");
    return s._values[index];
  }

  /**
   * Go to next element.
   */
  inline void next() {
    if (index < s._nnz) {
      index++;
    }
  }
};

}

