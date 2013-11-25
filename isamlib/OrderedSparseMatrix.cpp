/**
 * @file OrderedSparseMatrix.cpp
 * @brief Adds column ordering to sparse matrix functionality for iSAM.
 * @author Michael Kaess
 * @version $Id: OrderedSparseMatrix.cpp 6377 2012-03-30 20:06:44Z kaess $
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
#include <cstring> // memset()
#include <fstream>
#include <iostream>
#include <cmath>

#include "isam/util.h"

#include "isam/OrderedSparseMatrix.h"

using namespace std;

namespace isam {

void OrderedSparseMatrix::_allocate_OrderedSparseMatrix(bool init_table) {
  _r_to_a = new int[_max_num_cols];
  _a_to_r = new int[_max_num_cols];
  if (init_table) {
    for (int col=0; col<_num_cols; col++) {
      // reset translation table
      _r_to_a[col] = col;
      _a_to_r[col] = col;
    }
  }
}

void OrderedSparseMatrix::_copy_from_OrderedSparseMatrix(const OrderedSparseMatrix& mat) {
  _allocate_OrderedSparseMatrix(false);
  memcpy(_r_to_a, mat._r_to_a, _num_cols*sizeof(int));
  memcpy(_a_to_r, mat._a_to_r, _num_cols*sizeof(int));
}

void OrderedSparseMatrix::_dealloc_OrderedSparseMatrix() {
  delete[] _r_to_a;
  _r_to_a = NULL; // recommended, as freeing a second time has unpredictable consequences
  delete[] _a_to_r;
  _a_to_r = NULL;
}

void OrderedSparseMatrix::_calc_reverse_order(int num, const int* order, int* reverse_order) const {
  for (int i=0; i<num; i++) {
    reverse_order[order[i]] = i;
  }
}

void OrderedSparseMatrix::_set_order(const int* r_to_a) {
  memcpy(_r_to_a, r_to_a, sizeof(int)*_num_cols);
  _calc_reverse_order(_num_cols, r_to_a, _a_to_r);
}

OrderedSparseMatrix::OrderedSparseMatrix(int num_rows, int num_cols) : SparseMatrix(num_rows, num_cols) {
  _allocate_OrderedSparseMatrix();
}

OrderedSparseMatrix::OrderedSparseMatrix(const OrderedSparseMatrix& mat) : SparseMatrix(mat) {
  _copy_from_OrderedSparseMatrix(mat);
}

OrderedSparseMatrix::OrderedSparseMatrix(const OrderedSparseMatrix& mat,
                                         int num_rows, int num_cols, int first_row, int first_col) :
    SparseMatrix(mat, num_rows, num_cols, first_row, first_col)
{
  _allocate_OrderedSparseMatrix(); // note: ignores original ordering
}

OrderedSparseMatrix::OrderedSparseMatrix(int num_rows, int num_cols, SparseVector_p* rows) :
    SparseMatrix(num_rows, num_cols, rows)
{
  _allocate_OrderedSparseMatrix();
}

OrderedSparseMatrix::~OrderedSparseMatrix() {
  _dealloc_OrderedSparseMatrix();
}

const OrderedSparseMatrix& OrderedSparseMatrix::operator= (const OrderedSparseMatrix& mat) {
  if (this==&mat)
    return *this;

  SparseMatrix::operator=(mat);

  _dealloc_OrderedSparseMatrix();
  _copy_from_OrderedSparseMatrix(mat);

  return *this;
}

void OrderedSparseMatrix::set_row(int row, const SparseVector& new_row) {
  // translate according to variable order:
  SparseVector reordered_row;
  for (SparseVectorIter iter(new_row); iter.valid(); iter.next()) {
    double val;
    int col = iter.get(val);
    int trans = a_to_r()[col];
    reordered_row.set(trans, val);
  }
  SparseMatrix::set_row(row, reordered_row);
}

void OrderedSparseMatrix::import_rows_ordered(int num_rows, int num_cols, SparseVector_p* rows, int* r_to_a) {
  _dealloc_OrderedSparseMatrix();
  SparseMatrix::import_rows(num_rows, num_cols, rows);
  _allocate_OrderedSparseMatrix();
  _set_order(r_to_a);
}

void OrderedSparseMatrix::append_new_cols(int num) {
  int orig_num_cols = _num_cols;
  int orig_max_num_cols = _max_num_cols;

  SparseMatrix::append_new_cols(num); // _num_cols and _max_num_cols

  if (orig_max_num_cols != _max_num_cols) { // resize arrays if needed
    int* new_r_to_a = new int[_max_num_cols];
    int* new_a_to_r = new int[_max_num_cols];
    memcpy(new_r_to_a, _r_to_a, orig_max_num_cols*sizeof(int));
    memcpy(new_a_to_r, _a_to_r, orig_max_num_cols*sizeof(int));
    delete[] _r_to_a;
    delete[] _a_to_r;
    _r_to_a = new_r_to_a;
    _a_to_r = new_a_to_r;
  }
  int next_index = orig_num_cols;
  for (int col=orig_num_cols; col<orig_num_cols+num; col++) {
    _r_to_a[col] = next_index;
    _a_to_r[next_index] = col;
    next_index++;
  }
}

const int* OrderedSparseMatrix::a_to_r() const {
  return _a_to_r;
}

const int* OrderedSparseMatrix::r_to_a() const {
  return _r_to_a;
}

}
