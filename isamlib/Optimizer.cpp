/**
 * @file Optimizer.cpp
 * @brief Class implementing batch and incremental nonlinear equation solvers.
 * @author David Rosen
 * @author Michael Kaess
 * @version $Id: Optimizer.cpp 6371 2012-03-29 22:22:23Z kaess $
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

#include <Eigen/Dense>

#include "isam/Optimizer.h"
#include "isam/OptimizationInterface.h"

using namespace std;
using namespace Eigen;

namespace isam {

/* Use Powell's Dog-Leg stopping criteria for all of the batch algorithms? */
// #define USE_PDL_STOPPING_CRITERIA

void Optimizer::permute_vector(const VectorXd& v, VectorXd& p,
    const int* permutation) {
  for (int i = 0; i < v.size(); i++) {
    p(permutation[i]) = v(i);
  }
}

VectorXd Optimizer::compute_gauss_newton_step(const SparseSystem& jacobian,
    SparseSystem* R, double lambda) {
  VectorXd delta_ordered;
  _cholesky->factorize(jacobian, &delta_ordered, lambda);
  if (R != NULL) {
    _cholesky->get_R(*R);
  }

  // delta has new ordering, need to return result with default ordering
  int nrows = delta_ordered.size();
  VectorXd delta(nrows);
  permute_vector(delta_ordered, delta, _cholesky->get_order());

  return delta;
}

VectorXd Optimizer::compute_dog_leg(double alpha, const VectorXd& h_sd,
    const VectorXd& h_gn, double delta, double& gain_ratio_denominator) {
  if (h_gn.norm() <= delta) {
    gain_ratio_denominator = current_SSE_at_linpoint;
    return h_gn;
  }

  double h_sd_norm = h_sd.norm();

  if ((alpha * h_sd_norm) >= delta) {
    gain_ratio_denominator = delta * (2 * alpha * h_sd_norm - delta)
        / (2 * alpha);
    return (delta / h_sd_norm) * h_sd;
  } else {
    // complicated case: calculate intersection of trust region with
    // line between Gauss-Newton and steepest descent solutions
    VectorXd a = alpha * h_sd;
    VectorXd b = h_gn;
    double c = a.dot(b - a);
    double b_a_norm2 = (b - a).squaredNorm();
    double a_norm2 = a.squaredNorm();
    double delta2 = delta * delta;
    double sqrt_term = sqrt(c * c + b_a_norm2 * (delta2 - a_norm2));
    double beta;
    if (c <= 0) {
      beta = (-c + sqrt_term) / b_a_norm2;
    } else {
      beta = (delta2 - a_norm2) / (c + sqrt_term);
    }

    gain_ratio_denominator = .5 * alpha * (1 - beta) * (1 - beta) * h_sd_norm
        * h_sd_norm + beta * (2 - beta) * current_SSE_at_linpoint;
    return (alpha * h_sd + beta * (h_gn - alpha * h_sd));
  }
}

void Optimizer::update_trust_radius(double rho, double hdl_norm) {
  if (rho < .25) {
    Delta /= 2.0;
  }
  if (rho > .75) {
    Delta = max(Delta, 3 * hdl_norm);
  }
}

void Optimizer::relinearize(const Properties& prop) {
  // We're going to relinearize about the current estimate.
  function_system.estimate_to_linpoint();

  // prepare factorization
  SparseSystem jac = function_system.jacobian();

  // factorization and new rhs based on new linearization point will be in _R
  VectorXd h_gn = compute_gauss_newton_step(jac, &function_system._R); // modifies _R

  if (prop.method == DOG_LEG) {
    // Compute the gradient and cache it.
    gradient = mul_SparseMatrixTrans_Vector(jac, jac.rhs());

    //Get the value of the sum-of-squared errors at the current linearization point.
    current_SSE_at_linpoint = jac.rhs().squaredNorm();

    // NB: alpha's denominator will be zero iff the gradient vector is zero
    // (since J is full-rank by hypothesis).  But the gradient is zero iff
    // we're already at the minimum, so we don't actually need to to any
    // updates; just set the estimate to be the current linearization point,
    // since we're already at the minimum.

    double alpha_denominator = (jac * gradient).squaredNorm();

    if (alpha_denominator > 0) {
      double alpha_numerator = gradient.squaredNorm();
      double alpha = alpha_numerator / alpha_denominator;

      // These values will be used to update the estimate of Delta
      double F_0, F_h;

      F_0 = current_SSE_at_linpoint;

      double rho_denominator, rho;

      VectorXd h_dl;

      do {
        // We repeat the computation of the dog-leg step, shrinking the
        // trust-region radius if necessary, until we generate a sufficiently
        // small region of trust that we accept the proposed step.

        // Compute dog-leg step.
        // NOTE: Here we use -h_gn because of the weird sign change in the exmap functions.
        h_dl = compute_dog_leg(alpha, -gradient, -h_gn, Delta, rho_denominator);

        // Update the estimate.
        // NOTE:  Here we use -h_dl because of the weird sign change in the exmap functions.
        function_system.apply_exmap(-h_dl);

        // Get the value of the sum-of-squared errors at the new estimate.
        F_h = function_system.weighted_errors(ESTIMATE).squaredNorm();

        // Compute gain ratio.
        rho = (F_0 - F_h) / (rho_denominator);

        update_trust_radius(rho, h_dl.norm());
      } while (rho < 0);

      // Cache last accepted dog-leg step
      last_accepted_hdl = h_dl;
    } else {
      function_system.linpoint_to_estimate();
      last_accepted_hdl = VectorXd::Zero(gradient.size());
    }

  } else {
    // For Gauss-Newton just apply the update directly.
    function_system.apply_exmap(h_gn);
  }
}

bool Optimizer::powells_dog_leg_update(double epsilon1, double epsilon3,
    SparseSystem& jacobian, VectorXd& f_x, VectorXd& grad) {
  jacobian = function_system.jacobian();
  f_x = function_system.weighted_errors(LINPOINT);
  grad = mul_SparseMatrixTrans_Vector(jacobian, f_x);
  return (f_x.lpNorm<Eigen::Infinity>() <= epsilon3)
      || (grad.lpNorm<Eigen::Infinity>() <= epsilon1);
}

void Optimizer::augment_sparse_linear_system(SparseSystem& W,
    const Properties& prop) {
  if (prop.method == DOG_LEG) {
    // We're using the incremental version of Powell's Dog-Leg, so we need
    // to form the updated gradient.
    const VectorXd& f_new = W.rhs();

    // Augment the running count for the sum-of-squared errors at the current
    // linearization point.
    current_SSE_at_linpoint += f_new.squaredNorm();

    // Allocate the new gradient vector
    VectorXd g_new(W.num_cols());

    // Compute W^T \cdot f_new
    VectorXd increment = mul_SparseMatrixTrans_Vector(W, f_new);

    // Set g_new = (g_old 0)^T + W^T f_new.
    g_new.head(gradient.size()) = gradient + increment.head(gradient.size());
    g_new.tail(W.num_cols() - gradient.size()) = increment.tail(
        W.num_cols() - gradient.size());

    // Cache the new gradient vector
    gradient = g_new;
  }

  // Apply Givens to QR factorize the newly augmented sparse system.
  for (int i = 0; i < W.num_rows(); i++) {
    SparseVector new_row = W.get_row(i);
    function_system._R.add_row_givens(new_row, W.rhs()(i));
  }
}

void Optimizer::update_estimate(const Properties& prop) {
  // Solve for the Gauss-Newton step.
  VectorXd h_gn_reordered = function_system._R.solve();

  // permute from R-ordering to J-ordering
  VectorXd h_gn(h_gn_reordered.size());
  permute_vector(h_gn_reordered, h_gn, function_system._R.r_to_a());

  if (prop.method == GAUSS_NEWTON) {
    function_system.apply_exmap(h_gn);
  } else { //method == DOG_LEG
    // Compute alpha.  Note that since the variable ordering of the factor
    // R differs from that of the original Jacobian J,
    // we must first rearrange the ordering of the elements in the gradient.
    VectorXd reordered_gradient(gradient.size());

    // Permute from J-ordering to R-ordering
    permute_vector(gradient, reordered_gradient, function_system._R.a_to_r());

    double alpha_denominator =
        (function_system._R * reordered_gradient).squaredNorm();

    if (alpha_denominator > 0) {
      double alpha = (gradient.squaredNorm()) / alpha_denominator;

      double rho_denominator;

      VectorXd h_dl = compute_dog_leg(alpha, -gradient, -h_gn, Delta,
          rho_denominator);
      function_system.apply_exmap(-h_dl);

      // Compute the gain ratio
      double rho = (current_SSE_at_linpoint
          - function_system.weighted_errors(ESTIMATE).squaredNorm())
          / rho_denominator;

      update_trust_radius(rho, h_dl.norm());

      if (rho < 0) {
        // The proposed update actually /increased/ the value of the
        // objective function; restore the last good estimate that we had
        VectorXd restore_step(gradient.size());
        restore_step.head(last_accepted_hdl.size()) = last_accepted_hdl;
        restore_step.tail(gradient.size() - last_accepted_hdl.size()).setZero();

        function_system.apply_exmap(-restore_step);
      } else {
        // The proposed update was accepted; cache the dog-leg step used
        // to produce it.
        last_accepted_hdl = h_dl;
      }
    }

    // NOTE:  The negatives prepended to "compute_dog_leg()" and "h_gn"
    // are due to the weird sign change in the exmap functions.
  }
}

void Optimizer::gauss_newton(const Properties& prop, int* num_iterations) {
  // Batch optimization
  int num_iter = 0;

  // Set the new linearization point to be the current estimate.
  function_system.estimate_to_linpoint();
  // Compute Jacobian about current estimate.
  SparseSystem jacobian = function_system.jacobian();

  // Get the current error residual vector
  VectorXd r = function_system.weighted_errors(LINPOINT);

#ifdef USE_PDL_STOPPING_CRITERIA
  // Compute the current gradient direction vector
  VectorXd g = mul_SparseMatrixTrans_Vector(jacobian, r);
#else
  double error = r.squaredNorm();
  double error_new;
  // We haven't computed a step yet, so this initialization is to ensure
  // that we never skip over the while loop as a result of failing
  // change-in-error check.
  double error_diff = prop.epsilon_rel * error + 1;
#endif

  // Compute Gauss-Newton step h_{gn} to get to the next estimated optimizing point.
  VectorXd delta = compute_gauss_newton_step(jacobian);

  while (
  // We ALWAYS use these criteria
  ((prop.max_iterations <= 0) || (num_iter < prop.max_iterations))
      && (delta.norm() > prop.epsilon2)

#ifdef USE_PDL_STOPPING_CRITERIA
      && (r.lpNorm<Eigen::Infinity>() > prop.epsilon3)
      && (g.lpNorm<Eigen::Infinity>() > prop.epsilon1)

#else // Custom stopping criteria for GN
      && (error > prop.epsilon_abs)
      && (fabs(error_diff) > prop.epsilon_rel * error)
#endif

  ) // end while conditional
  {
    num_iter++;

    // Apply the Gauss-Newton step h_{gn} to...
    function_system.apply_exmap(delta);
    // ...set the new linearization point to be the current estimate.
    function_system.estimate_to_linpoint();
    // Relinearize about the new current estimate.
    jacobian = function_system.jacobian();
    // Compute the error residual vector at the new estimate.
    r = function_system.weighted_errors(LINPOINT);

#ifdef USE_PDL_STOPPING_CRITERIA
    g = mul_SparseMatrixTrans_Vector(jacobian, r);
#else
    // Update the error difference in errors between the previous and
    // current estimates.
    error_new = r.squaredNorm();
    error_diff = error - error_new;
    error = error_new;  // Record the absolute error at the current estimate
#endif

    // Compute Gauss-Newton step h_{gn} to get to the next estimated
    // optimizing point.
    delta = compute_gauss_newton_step(jacobian);
    if (!prop.quiet) {
      cout << "Iteration " << num_iter << ": residual ";

#ifdef USE_PDL_STOPPING_CRITERIA
      cout << r.squaredNorm();
#else
      cout << error;
#endif
      cout << endl;
    }
  }  //end while

  if (num_iterations != NULL) {
    *num_iterations = num_iter;
  }
  _cholesky->get_R(function_system._R);
}

void Optimizer::levenberg_marquardt(const Properties& prop,
    int* num_iterations) {
  int num_iter = 0;
  double lambda = prop.lm_lambda0;
  // Using linpoint as current estimate below.
  function_system.estimate_to_linpoint();

  // Get the current Jacobian at the linearization point.
  SparseSystem jacobian = function_system.jacobian();

  // Get the error residual vector at the current linearization point.
  VectorXd r = function_system.weighted_errors(LINPOINT);

  // Record the absolute sum-of-squares error (i.e., objective function value) here.
  double error = r.squaredNorm();

#ifdef USE_PDL_STOPPING_CRITERIA
  // Compute the gradient direction vector at the current linearization point
  VectorXd g = mul_SparseMatrixTrans_Vector(jacobian, r);
#endif

  double error_diff, error_new;

  // solve at J'J + lambda*diag(J'J)
  VectorXd delta = compute_gauss_newton_step(jacobian, &function_system._R,
      lambda);

  while (
  // We ALWAYS use these stopping criteria
  ((prop.max_iterations <= 0) || (num_iter < prop.max_iterations))
      && (delta.norm() > prop.epsilon2)

#ifdef USE_PDL_STOPPING_CRITERIA
      && (r.lpNorm<Eigen::Infinity>() > prop.epsilon3)
      && (g.lpNorm<Eigen::Infinity>() > prop.epsilon1)
#else
      && (error > prop.epsilon_abs)
#endif
  )  // end while conditional
  {
    num_iter++;

    // remember the last accepted linearization point
    function_system.linpoint_to_estimate();
    // Apply the delta vector DIRECTLY TO THE LINEARIZATION POINT!
    function_system.self_exmap(delta);
    error_new = function_system.weighted_errors(LINPOINT).squaredNorm();
    error_diff = error - error_new;
    // feedback
    if (!prop.quiet) {
      cout << "LM Iteration " << num_iter << ": (lambda=" << lambda << ") ";
      if (error_diff > 0.) {
        cout << "residual: " << error_new << endl;
      } else {
        cout << "rejected" << endl;
      }
    }
    // decide if acceptable
    if (error_diff > 0.) {

#ifndef USE_PDL_STOPPING_CRITERIA
      if (error_diff < prop.epsilon_rel * error) {
        break;
      }
#endif

      // Update lambda
      lambda /= prop.lm_lambda_factor;

      // Record the error at the newly-accepted estimate.
      error = error_new;

      // Relinearize around the newly-accepted estimate.
      jacobian = function_system.jacobian();

#ifdef USE_PDL_STOPPING_CRITERIA
      r = function_system.weighted_errors(LINPOINT);
      g = mul_SparseMatrixTrans_Vector(jacobian, r);
#endif
    } else {
      // reject new estimate
      lambda *= prop.lm_lambda_factor;
      // restore previous estimate
      function_system.estimate_to_linpoint();
    }

    // Compute the step for the next iteration.
    delta = compute_gauss_newton_step(jacobian, &function_system._R, lambda);

  } // end while

  if (num_iterations != NULL) {
    *num_iterations = num_iter;
  }
  // Copy current estimate contained in linpoint.
  function_system.linpoint_to_estimate();
}

void Optimizer::powells_dog_leg(int* num_iterations, double delta0,
    int max_iterations, double epsilon1, double epsilon2, double epsilon3) {
  // Batch optimization
  int num_iter = 0;
  // current estimate is used as new linearization point
  function_system.estimate_to_linpoint();

  double delta = delta0;
  SparseSystem jacobian(1, 1);
  VectorXd f_x;
  VectorXd grad;

  bool found = powells_dog_leg_update(epsilon1, epsilon3, jacobian, f_x, grad);

  double rho_denominator;

  while ((not found) && (max_iterations == 0 || num_iter < max_iterations)) {
    num_iter++;
    cout << "PDL Iteration " << num_iter << " residual: " << f_x.squaredNorm()
        << endl;
    // compute alpha
    double alpha = grad.squaredNorm() / (jacobian * grad).squaredNorm();
    // steepest descent
    VectorXd h_sd = -grad;
    // solve Gauss Newton
    VectorXd h_gn = compute_gauss_newton_step(jacobian, &function_system._R);
    // compute dog leg h_dl
    // x0 = x: remember (and return) linearization point of R
    function_system.linpoint_to_estimate();
    VectorXd h_dl = compute_dog_leg(alpha, h_sd, h_gn, delta, rho_denominator);
    // Evaluate new solution, update estimate and trust region.
    if (h_dl.norm() <= epsilon2) {
      found = true;
    } else {
      // new estimate
      // change linearization point directly (original LP saved in estimate)
      function_system.self_exmap(h_dl);
      // calculate gain ratio rho
      VectorXd f_x_new = function_system.weighted_errors(LINPOINT);
      double rho = (f_x.squaredNorm() - f_x_new.squaredNorm())
          / (rho_denominator);
      if (rho > 0) {
        // accept new estimate
        cout << "accepted" << endl;
        f_x = f_x_new;
        found = powells_dog_leg_update(epsilon1, epsilon3, jacobian, f_x, grad);
      } else {
        // reject new estimate, overwrite with last saved one
        cout << "rejected" << endl;
        function_system.estimate_to_linpoint();
      }
      if (rho > 0.75) {
        delta = max(delta, 3.0 * h_dl.norm());
      } else if (rho < 0.25) {
        delta *= 0.5;
        found = (delta <= epsilon2);
      }
    }
  }
  if (num_iterations) {
    *num_iterations = num_iter;
  }
  // Overwrite potentially rejected linearization point with last saved one
  // (could be identical if it was accepted in the last iteration).

  function_system.swap_estimates();

}

void Optimizer::batch_optimize(const Properties& prop, int* num_iterations) {

  const double delta0 = 1.0;

  switch (prop.method) {
  case GAUSS_NEWTON:
    gauss_newton(prop, num_iterations);
    break;
  case DOG_LEG:
    powells_dog_leg(num_iterations, delta0, prop.max_iterations, prop.epsilon1,
        prop.epsilon2, prop.epsilon3); // modifies x0,R
    break;
  case LEVENBERG_MARQUARDT:
    levenberg_marquardt(prop, num_iterations);
    break;
  }

}

}
