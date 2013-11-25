/**
 * @file Optimizer.h
 * @brief Class implementing batch and incremental nonlinear equation solvers.
 * @author Michael Kaess
 * @author David Rosen
 * @version $Id: Optimizer.h 6368 2012-03-28 23:01:19Z kaess $
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

#include "Properties.h"
#include "OptimizationInterface.h"
#include "SparseSystem.h"
#include "Cholesky.h"
#include "Node.h"

namespace isam {

class Optimizer {

private:

  /**
   * Provides an interface for manipulating the linearized system in
   * the Slam class.
   */
  OptimizationInterface& function_system;

  /**
   * This data member is used to compute the thin QR decomposition of
   * the linear system during the relinearization steps.
   */
  Cholesky* _cholesky;

  /**
   * Cached gradient vector; only used with increment Powell's Dog-Leg.
   */
  Eigen::VectorXd gradient;

  /**
   * Current radius of trust region; only used with Powell's Dog-Leg.
   */
  double Delta;

  /**
   * Used to restore the previous estimate for cases in which the proposed step
   * is rejected; only used with Powell's Dog-Leg in incremental mode.
   */
  Eigen::VectorXd last_accepted_hdl;

  /**
   * This keeps a running count of the sum of squared errors at the
   * linearization point; only used with Powell's Dog-Leg in incremental mode.
   */
  double current_SSE_at_linpoint;


  void update_trust_radius(double rho, double hdl_norm);

  /**
   * Computes and returns the dog-leg step given the parameter alpha, the
   * trust region radius delta, and the steepest-descent and Gauss-Newton
   * steps. Also computes and returns the value of the denominator that
   * will be used to compute the gain ratio.
   */
  Eigen::VectorXd compute_dog_leg(double alpha, const Eigen::VectorXd& h_sd,
      const Eigen::VectorXd& h_gn, double delta,
      double& gain_ratio_denominator);

  bool powells_dog_leg_update(double epsilon1, double epsilon3,
      SparseSystem& jacobian, Eigen::VectorXd& f_x, Eigen::VectorXd& gradient);

  /**
   * Given an input vector v and an array of ints representing a permutation,
   * applies the permutation to the elements of v and returns the permuted vector
   * @param v Input vector.
   * @param permutation  An array of ints representing the permutation to be
   *                     applied to the elements of v.
   * @return  The returned permuted vector p satisfies p(permutation[i]) = v(i)
   *          i.e., p is formed by mapping the ith element of v to the
   *          permutation[i]-th element of p.
   */
  void permute_vector(const Eigen::VectorXd& v, Eigen::VectorXd& p,
      const int* permutation);

  /**
   * Helper method for computing the Gauss-Newton step h_{gn} in Gauss-Newton,
   * Levenberg-Marquardt, and Powell's dog-leg algorithms in batch mode.
   * Specifically, this function can compute the Gauss-Newton step h_gn as
   * part of the factorization of the relinearized SparseSystem computed at
   * each iteration of batch processing algorithm.
   *
   * @param jacobian  The SparseSystem representing the linearization
   * @param R
   * @param lambda
   * @return h_gn
   */
  Eigen::VectorXd compute_gauss_newton_step(const SparseSystem& jacobian,
      SparseSystem* R = NULL, double lambda = 0.);

  void gauss_newton(const Properties& prop, int* num_iterations = NULL);

  /**
   * Perform Levenberg-Marquardt
   * @param prop Properties including stopping criteria max_iterations and
   *             epsilon, and lm_... parameters
   * @param num_iterations Upon return, contains number of iterations performed.
   */
  void levenberg_marquardt(const Properties& prop, int* num_iterations = NULL);

  /**
   * Powell's dog leg algorithm, a trust region method that combines
   * Gauss-Newton and steepest descent, similar to Levenberg-Marquardt,
   * but requiring less iterations (matrix factorizations). See
   * Manolis05iccv for a comparison to LM. Implemented according to
   * Madsen, Nielson, Tingleff, "Methods for Non-Linear Least Squares
   * Problems", lecture notes, Denmark, 2004
   * (http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf).
   * @param num_iterations Contains number of iterations on return if not NULL.
   * @param delta0 Initial trust region.
   * @param max_iterations Maximum number of iterations (0 means unlimited).
   * @param epsilon1
   * @param epsilon2
   * @param epsilon3
   */
  void powells_dog_leg(int* num_iterations = NULL, double delta0 = 1.0,
      int max_iterations = 0, double epsilon1 = 1e-4, double epsilon2 = 1e-4,
      double epsilon3 = 1e-4);

public:

  Optimizer(OptimizationInterface& fs)
      : function_system(fs), Delta(1.0) {
    //Initialize the Cholesky object
    _cholesky = Cholesky::Create();
  }

  /**
   * Perform batch optimization using the method set in prop
   */
  void batch_optimize(const Properties& prop, int* num_iterations);

  /**
   * Used to augment the sparse linear system by adding new measurements.
   * Only useful in incremental mode.
   */
  void augment_sparse_linear_system(SparseSystem& W, const Properties& prop);

  /**
   * Computes the Jacobian J(x_est) of the residual error function about the
   * current estimate x_est, computes the thin QR factorization of J(x_est),
   * and then stores (R,d) as a SparseSystem, where
   *
   * Q^t f(x_est) = | d |
   *                | e |
   *
   * and || e ||^2 is the squared residual error.
   *
   *
   *
   * NOTA BENE:  (Deeply) Internally, this algorithm uses the SuiteSparse
   * library to perform the matrix decomposition shown above.  The SuiteSparse
   * library utilizes variable reordering in order to reduce fill-in and boost
   * efficiency.  Consequently, the ordering (i.e., naming) of variables in the
   * SparseSystem computed by this function differs from the ordering (naming)
   * of the variables in the Graph object contained in the Slam class.
   *
   * The mappings between the two orderings can be obtained using
   * OrderedSparseMatrix::a_to_r() and OrderedSparseMatrix::r_to_a().  These
   * functions return const int*'s that point to internal arrays encoding the
   * permutation between the naming of variables passed in to the factorization
   * routine, and the naming of variables in the factored matrices returned
   * by the relinearization.
   *
   * More precisely, if
   *
   * const int* order = function_system._R.a_to_r();
   *
   * then the variable x_i from the Graph object stored the Slam class is
   * mapped to the variable x_{order[i]} in the SparseSystem obtained after
   * linearization.
   *
   * The call
   *
   * const int* inverse_order = function_system._R.r_to_a();
   *
   * retrieves the inverse permutation.
   */
  void relinearize(const Properties& prop);

  /**
   * Updates the current estimated solution
   */
  void update_estimate(const Properties& prop);

  ~Optimizer() {
    delete _cholesky;
  }

};

}
