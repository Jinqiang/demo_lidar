/**
 * @file Cholesky.h
 * @brief Cholesky batch factorization using SuiteSparse by Tim Davis.
 * @author Michael Kaess
 * @version $Id: Cholesky.h 6377 2012-03-30 20:06:44Z kaess $
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

#include "SparseSystem.h"

namespace isam {

class Cholesky {
public:
  virtual ~Cholesky() {}

  /**
   * Factorize a given system Ax=b and optionally solve.
   * @param Ab SparseSystem with measurement Jacobian A and right hand side b.
   * @param delta Optional parameter to return solution of system.
   * @param lambda Adds elements to diagonal of information matrix A'A before
   *        factorization, used for Levenberg-Marquardt algorithm.
   */
  virtual void factorize(const SparseSystem& Ab, Eigen::VectorXd* delta = NULL, double lambda = 0.) = 0;

  /**
   * Copy R into a SparseSystem data structure (expensive, so can be
   * avoided during batch factorization).
   * @param R SparseSystem that upon return will contain the R factor.
   */
  virtual void get_R(SparseSystem& R) = 0;

  /**
   * Access the variable ordering used for Cholesky factorization.
   * @return Pointer to variable ordering.
   */
  virtual int* get_order() = 0;

  static Cholesky* Create();

protected:
  Cholesky() {}
};

}
