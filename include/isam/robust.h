/**
 * @file robust.h
 * @brief Robust estimator functionality.
 * @author Michael Kaess
 * @version $Id: robust.h 6377 2012-03-30 20:06:44Z kaess $
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

#include <cmath>

#include "util.h"

namespace isam {

/*
 * robust error functions following Hartley&Zisserman book (2nd edition, pages 616-622)
 * summary:
 * - squared error is not robust
 * - Blake-Zisserman, corrupted Gaussian and Cauchy functions are
 *   non-convex, and therefore require very good initialization
 * - L1 yields a stable minimum, the median, but is non-differentiable at the origin
 * - Huber cost function is convex, combining squared near the origin, L1 further out
 * - pseudo-Huber is a good alternative to Huber, with continuous derivatives of all orders
 */

/**
 * Standard squared cost function.
 * @param d Unmodified cost/distance.
 */
inline double cost_squared(double d) {
  return d*d;
}

/**
 * Robust Blake-Zisserman cost function.
 * @param d Unmodified cost/distance.
 * @param e Approximate crossover between inliers and outliers: d^2 = -log(e).
 */
inline double cost_blake_zisserman(double d, double e) {
  return -log(exp(-d*d) + e);
}

/**
 * Robust cost function using corrupted Gaussian.
 * @param d Unmodified cost/distance.
 * @param w Ratio of standard deviations of the outliers to the inliers.
 * @param a Expected fraction of inliers.
 */
inline double cost_corrupted_gaussian(double d, double w, double a) {
  double d2 = d*d;
  double w2 = w*w;
  return -log(a*exp(-d2) + (1.-a)*exp(-d2/w2)/w);
}

/**
 * Robust Cauchy cost function.
 * @param d Unmodified cost/distance.
 * @param b Determines for which range of d the function closely
 *          approximates the squared error function.
 */
inline double cost_cauchy(double d, double b = 1.) {
  // the first term is a constant and could be omitted
  return log(M_PI/b) * log(1. + d*d/(b*b));
}

/**
 * Robust L1 cost function, the sum of absolute errors.
 * @param d Unmodified cost/distance.
 * @param b ?
 */
inline double cost_l1(double d, double b = 0.5) {
  return 2.*b*fabs(d);
}

/**
 * Robust Huber cost function.
 * @param d Unmodified cost/distance.
 * @param b Approximately the outlier threshold.
 */
inline double cost_huber(double d, double b) {
  double abs_d = fabs(d);
  if (abs_d < b) {
    return d*d;
  } else {
    return 2*b*abs_d - b*b;
  }
}

/**
 * Robust Pseudo-Huber cost function.
 * @param d Unmodified cost/distance.
 * @param b ?
 */
inline double cost_pseudo_huber(double d, double b) {
  double b2 = b*b;
  return 2*b2*(sqrt(1+d*d/b2) - 1);
}

}
