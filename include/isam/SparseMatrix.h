/**
 * @file SparseMatrix.h
 * @brief Sparse matrix functionality for iSAM
 * @author Michael Kaess
 * @version $Id: SparseMatrix.h 6376 2012-03-30 18:34:44Z kaess $
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

#include <iostream>
#include <Eigen/Dense>

#include "SparseVector.h"

namespace isam {

typedef SparseVector* SparseVector_p;

class SparseMatrix {
  int _num_rows; // current number of rows
  int _num_cols; // current number of columns
  int _max_num_rows; // allocated number of row indices
  int _max_num_cols; // allocated number of columns
  SparseVector_p* _rows;  // pointers to the actual rows

  /**
   * Allocate memory - private.
   * @param num_rows Number of active rows.
   * @param num_cols Number of active columns.
   * @param max_num_rows Number of rows to reserve memory for.
   * @param max_num_cols Number of columns to reserve memory for.
   */
  void _allocate_SparseMatrix(int num_rows, int num_cols, int max_num_rows, int max_num_cols, bool init_rows = true);

  /**
   * Copy data from one sparse matrix to a new one - private.
   * @param vec Existing sparse vector to copy from.
   */
  void _copy_from_SparseMatrix(const SparseMatrix& mat);

  /**
   * Deallocate memory - private.
   */
  void _dealloc_SparseMatrix();

public:

  /**
   * Constructor.
   * @param num_rows Initial number of rows.
   * @param num_cols Initial number of columns.
   */
  SparseMatrix(int num_rows, int num_cols);

  /**
   * Copy constructor.
   * @param mat Matrix to copy.
   */
  SparseMatrix(const SparseMatrix& mat);

  /**
   * Submatrix copy constructor
   * @param mat Matrix to copy.
   * @param num_rows Number of rows to copy.
   * @param num_cols Number of columns to copy.
   * @param first_row Row offset.
   * @param first_col Column offset.
   */
  SparseMatrix(const SparseMatrix& mat, int num_rows, int num_cols, int first_row = 0, int first_col = 0);

  SparseMatrix(int num_rows, int num_cols, SparseVector_p* rows);

  /**
   * Destructor.
   */
  virtual ~SparseMatrix();

  /**
   * Assignment operator.
   * @param mat Right-hand-side matrix in assignment
   * @return self.
   */
  const SparseMatrix& operator= (const SparseMatrix& mat);

  /**
   * Read a matrix entry.
   * @param row Row of entry.
   * @param col Column of entry.
   * @return Value of entry.
   */
  virtual double operator()(int row, int col) const;

  /**
   * Set one entry of the matrix. Non-existing entries are created.
   * @param row Row of entry.
   * @param col Column of entry.
   * @param val New value of entry.
   * @param grow_matrix Enlarge matrix if entry outside current size (default: false).
   */
  virtual void set(int row, int col, const double val, bool grow_matrix = false);

  /**
   * Append a new entry to a row. Allows efficient adding of presorted elements.
   * @param row Row of new entry.
   * @param col Column of new entry - must be after last existing one in this row.
   * @param val Value of new entry.
   */
  virtual void append_in_row(int row, int col,const double val);

  /**
   * Determine number of non-zero entries in sparse matrix.
   * @return Number of non-zero entries.
   */
  virtual int nnz() const;

  /**
   * Determine maximum number of non-zero entries in any column.
   * @return Maximum number of non-zero entries in any column.
   */
  virtual int max_nz() const;

  /**
   * Print sparse matrix as triples to stream;
   * also readable by Matlab: "load R.txt; S=spconvert(R); spy(S)"
   * @param out Output stream.
   */
  virtual void print(std::ostream& out = std::cout) const;

  /**
   * Print sparse matrix as triples to file
   * @param file_name File name to write sparse matrix to.
   */
  virtual void print(const std::string& file_name) const;

  /**
   * Print size of matrix and number of entries for debugging.
   */
  virtual void print_stats() const;

  /**
   * Prints non-zero pattern to stdout.
   */
  virtual void print_pattern() const;

  /**
   * Save eps graphics file with sparse matrix patterns.
   */
  virtual void save_pattern_eps(const std::string& file_name) const;

  /**
   * Return a sparse row.
   * @param row Number of row to return.
   * @return Sparse row vector.
   */
  virtual const SparseVector& get_row(int row) const;

  /**
   * Replace the given row.
   * @param row Number of row to replace.
   * @param new_row New row vector.
   */
  virtual void set_row(int row, const SparseVector& new_row);

  /**
   * Import externally allocated rows.
   * @param num_rows Number of rows of new matrix.
   * @param num_cols Number of columns of new matrix.
   * @param rows Array of SparseVector of length num_rows.
   */
  virtual void import_rows(int num_rows, int num_cols, SparseVector_p* rows);

  /**
   * Append new rows to matrix.
   * @param num Number of rows to add.
   */
  virtual void append_new_rows(int num);

  /**
   * Append new columns to matrix.
   * @param num Number of columns to add.
   */
  virtual void append_new_cols(int num);

  /**
   * Grow matrix to given number of rows; ignore if already at least as large.
   * @param num_rows Number of rows.
   */
  virtual void ensure_num_rows(int num_rows);

  /**
   * Grow matrix to given number of columns; ignore if already at least as large.
   * @param num_cols Number of columns.
   */
  virtual void ensure_num_cols(int num_cols);

  /**
   * Removes the last row (used in SparseSystem::add_row_givens).
   */
  virtual void remove_row();

  /**
   * Zero out an entry by applying a Givens rotation. Note that both
   * sparse rows have to be completely 0 on the left of col.
   * @param row The row from which row_bot is taken.
   * @param col The column of row_bot that should become 0.
   * @param c_givens Returns cosine of givens rotation if not NULL.
   * @param s_givens Returns sine of givens rotation if not NULL.
   */
  virtual void apply_givens(int row, int col, double* c_givens = NULL, double* s_givens = NULL);

  /**
   * Triangulate matrix by applying Givens rotations to all entries below the diagonal.
   * @return Number of Givens rotations applied (for analysis).
   */
  virtual int triangulate_with_givens();

  inline int num_rows() const {return _num_rows;}
  inline int num_cols() const {return _num_cols;}

  friend class OrderedSparseMatrix;
};

const Eigen::VectorXd operator*(const SparseMatrix& lhs, const Eigen::VectorXd& rhs);
Eigen::VectorXd mul_SparseMatrixTrans_Vector(const SparseMatrix& lhs, const Eigen::VectorXd& rhs);
SparseMatrix sparseMatrix_of_matrix(const Eigen::MatrixXd& m);
Eigen::MatrixXd matrix_of_sparseMatrix(const SparseMatrix& s);

}
