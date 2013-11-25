/**
 * @file SparseSystem.h
 * @brief Adds rhs functionality to sparse matrix for iSAM.
 * @author Michael Kaess
 * @version $Id: SparseSystem.h 4133 2011-03-22 20:40:38Z kaess $
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

#include "OrderedSparseMatrix.h"

namespace isam {

class SparseSystem : public OrderedSparseMatrix {
  Eigen::VectorXd _rhs;
public:
  SparseSystem(int num_rows, int num_cols);
  SparseSystem(const SparseSystem& mat);
  SparseSystem(const SparseSystem& mat, int num_rows, int num_cols, int first_row = 0, int first_col = 0);
  SparseSystem(int num_rows, int num_cols, SparseVector_p* rows, const Eigen::VectorXd& rhs);
  virtual ~SparseSystem();
  const SparseSystem& operator= (const SparseSystem& mat);

  const Eigen::VectorXd& rhs() const {return _rhs;}
  void set_rhs(const Eigen::VectorXd& rhs) {_rhs = rhs;}

  // overridden functions

  /**
   * Note: While rows are passed in, the rhs is required to already
   * contain the new entry - necessary because we cannot change the
   * signature of the function.
   */
  void apply_givens(int row, int col, double* c_givens = NULL, double* s_givens = NULL);

  void append_new_rows(int num);

  // new functions

  /**
   * Insert a new row
   * @param new_row The new sparse measurement row to add.
   * @param new_r New right hand side entry.
   */
  virtual void add_row(const SparseVector& new_row, double new_r);

  /**
   * Insert a new measurement row and triangulate using Givens rotations
   * @param new_row The new sparse measurement row to add.
   * @param new_r New right hand side entry.
   * @return Number of Givens rotations applied (for analysis).
   */
  virtual int add_row_givens(const SparseVector& new_row, double new_r);

  /**
   * Solve equation system by backsubstitution.
   * @return Solution for x in Rx=b'
   */
  virtual Eigen::VectorXd solve() const;

};

}
