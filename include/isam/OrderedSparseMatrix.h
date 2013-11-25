/**
 * @file OrderedSparseMatrix.h
 * @brief Adds column ordering to sparse matrix functionality for iSAM.
 * @author Michael Kaess
 * @version $Id: OrderedSparseMatrix.h 6377 2012-03-30 20:06:44Z kaess $
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

#include "SparseVector.h"
#include "SparseMatrix.h"

#include <string>
#include <fstream>
#include <iostream>

namespace isam {

class OrderedSparseMatrix : public SparseMatrix {
  int* _r_to_a; // column variable order
  int* _a_to_r;

  /**
   * Allocate memory - private.
   * @param init_table Determines if index table is initialized with identity.
   */
  void _allocate_OrderedSparseMatrix(bool init_order = true);

  /**
   * Copy data from one sparse matrix to a new one - private.
   * @param vec Existing sparse vector to copy from.
   */
  void _copy_from_OrderedSparseMatrix(const OrderedSparseMatrix& mat);

  /**
   * Deallocate memory - private.
   */
  void _dealloc_OrderedSparseMatrix();

  /**
   * Calculate reverse ordering for O(1) access.
   * @param num Number of entries in ordering
   * @param order Input order.
   * @param reverse_order Output reverse order.
   */
  void _calc_reverse_order(int num, const int* order, int* reverse_order) const;

  /**
   * Sets the variable order, for example after batch reordering.
   * @param order Pointer to list of integer IDs as used by CSparse.
   */
  void _set_order(const int* r_to_a);

public:

  /**
   * Constructor.
   * @param num_rows Initial number of rows.
   * @param num_cols Initial number of columns.
   */
  OrderedSparseMatrix(int num_rows, int num_cols);

  /**
   * Copy constructor.
   * @param mat Matrix to copy.
   */
  OrderedSparseMatrix(const OrderedSparseMatrix& mat);

  /**
   * Submatrix copy constructor.
   * @param mat Matrix to copy.
   * @param num_rows Number of rows to copy.
   * @param num_cols Number of columns to copy.
   * @param first_row Row offset.
   * @param first_col Column offset.
   */
  OrderedSparseMatrix(const OrderedSparseMatrix& mat, int num_rows,
      int num_cols, int first_row = 0, int first_col = 0);

  OrderedSparseMatrix(int num_rows, int num_cols, SparseVector_p* rows);

  /**
  * Destructor.
  */
  virtual ~OrderedSparseMatrix();

  /**
   * Assignment operator.
   * @param mat Right-hand-side matrix in assignment
   * @return self.
   */
  const OrderedSparseMatrix& operator= (const OrderedSparseMatrix& mat);

  // overridden functions

  /**
   * Replace the given row. Also reorders the new vector according to internal ordering.
   * @param row Number of row to replace.
   * @param new_row New row vector.
   */
  void set_row(int row, const SparseVector& new_row);

  /**
   * Import externally allocated rows, and also set the ordering.
   * @param num_rows Number of rows of new matrix.
   * @param num_cols Number of columns of new matrix.
   * @param rows Array of SparseVector of length num_rows.
   * @param r_to_a Variable ordering.
   */
  void import_rows_ordered(int num_rows, int num_cols, SparseVector_p* rows, int* r_to_a);

  /**
   * Expand matrix to include new columns.
   * @param num Number of columns to add.
   */
  void append_new_cols(int num);

  // new functions

  /**
   * Return variable ordering
   * @return Array of integers that specifies for each column in A
   * which column it now corresponds to in R.
   */
  virtual const int* a_to_r() const;

  /**
   * Return inverse variable ordering
   * @return Array of integers that specifies for each column in R
   * which column it originally came from in A.
   */
  virtual const int* r_to_a() const;

};

}
