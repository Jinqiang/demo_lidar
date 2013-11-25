/**
 * @file Properties.h
 * @brief Properties class for easy access to internal parameters.
 * @author Michael Kaess
 * @version $Id: Properties.h 6377 2012-03-30 20:06:44Z kaess $
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

namespace isam {

enum Method {GAUSS_NEWTON, LEVENBERG_MARQUARDT, DOG_LEG};

/**
 * User changeable default parameters.
 */
class Properties {
  // default copy constructor and assignment operator in use
public:

  /** print additional information if true */
  bool verbose;

  /** omit all textual output if true */
  bool quiet;

  /** Ignore any symbolic derivatives provided in factors */
  bool force_numerical_jacobian;

  /** which optimization method to use */
  Method method;

  /** Powell's Dog-Leg termination criteria: stop whenever the infinity-norm
   * of the gradient direction vector falls below this threshold. */
  double epsilon1;
  /** Powell's Dog-Leg termination criteria:  stop whenever the 2-norm of
   * the correction step falls below this threshold. */
  double epsilon2;
  /** Powell's Dog-Leg termination criteria:  stop whenever the infinity-norm
   * of the error-residual vector falls below this threshold. */
  double epsilon3;

  /** Termination criterion:  stop whenever the absolute sum-of-squares error
   * falls below this threshold. */
  double epsilon_abs;
  /** Termination criterion:  stop whenever the /difference/ in absolute
   * sum-of-squares errors between two estimates falls below this fraction
   * of the /total/ absolute sum-of-squares errors. */
  double epsilon_rel;

  /** Maximum number of iterations */
  int max_iterations;
  /** Starting value for lambda in LM */
  double lm_lambda0;
  /** Factor for multiplying (failure) or dividing (success) lambda */
  double lm_lambda_factor;

  /** Only update R matrix/solution/batch every mod_update steps */
  int mod_update;
  /** Batch solve with variable reordering and relinearization every mod_batch steps */
  int mod_batch;
  /** For incremental steps, solve by backsubstitution every mod_solve steps */
  int mod_solve;

  // default parameters
  Properties() :
    verbose(false),
    quiet(false),

    force_numerical_jacobian(false),

    method(GAUSS_NEWTON),

    epsilon1(1e-2),
    epsilon2(1e-2),
    epsilon3(1e-2),

    epsilon_abs(1e-3),
    epsilon_rel(1e-5),

    max_iterations(500),

    lm_lambda0(1e-6),
    lm_lambda_factor(10.),

    mod_update(1),
    mod_batch(100),
    mod_solve(1)
  {}
};

}
