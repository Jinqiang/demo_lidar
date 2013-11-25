/**
 * @file util.h
 * @brief Basic utility functions that are independent of iSAM.
 * @author Michael Kaess
 * @version $Id: util.h 8263 2013-04-10 14:02:19Z carlevar $
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

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <string>
#include <Eigen/Dense>

namespace isam {

// some math constants
#ifndef PI
const double PI = 3.14159265358979323846;
#endif
const double TWOPI = 2.0*PI;
const double HALFPI = PI/2.0;

// values up to this constant are considered zero and removed from the matrix
const double NUMERICAL_ZERO = 1e-12;

/**
 * Return current system time in seconds.
 */
double tic();

/**
 * Remember and return system time in seconds.
 * @param id Name of time slot.
 */
double tic(std::string id);

/**
 * Return difference between current system time and t in seconds.
 * @param t0 Start time as returned by tic();
 */
double toc(double t0);

/**
 * Return difference between current system time and remembered time in
 * seconds, and add to accumulated time.
 * @param id Name of time slot.
 */
double toc(std::string id);

/**
 * Print a list of accumulated times and additional statistics
 * for each name used in tic/toc.
 */
void tictoc_print();

/**
 * Return the accumulated time.
 * @param id Name of time slot.
 */
double tictoc(std::string id);

/**
 * Return identity matrix.
 */
Eigen::MatrixXd eye(int num);

/**
 * Calculate Givens rotation so that a specific entry becomes 0.
 * @param a Diagonal entry from above.
 * @param b Entry to be zeroed out.
 * @param c Resulting cosine part.
 * @param s Resulting sine part.
 */
void givens(const double a, const double b, double& c, double& s);

/**
 * Normalize angle to be within the interval [-pi,pi].
 */
inline double standardRad(double t) {
  if (t >= 0.) {
    t = fmod(t+PI, TWOPI) - PI;
  } else {
    t = fmod(t-PI, -TWOPI) + PI;
  }
  return t;
}

inline double deg_to_rad(double d) {
  return (d/180.*PI);
}

inline double rad_to_deg(double r) {
  return (r/PI*180.);
}

/**
 * Calculate the pseudo inverse of an arbitrary matrix using the SVD
 * @input a Eigen matrix to invert
 * @input eps Numerical epsilon to determine zero singular values (defaults to std::numeric_limits::eps)
 * @return Eigen matrix of same type as a
 */
template<typename T>
T pinv(const T &a, double eps = std::numeric_limits<typename T::Scalar>::epsilon()) {

  bool m_lt_n = (a.rows() < a.cols());

  Eigen::JacobiSVD<T> svd;
  if (m_lt_n) {
      T tmp = a.transpose();
      svd = tmp.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
  } else {
      svd = a.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
  }

  typename T::Scalar tolerance = eps * std::max(a.cols(), a.rows()) * svd.singularValues().array().abs().maxCoeff();

  T result = svd.matrixV() *
             T( (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0) ).asDiagonal() *
             svd.matrixU().adjoint();

  if (m_lt_n) {
      return result.transpose();
  } else {
      return result;
  }
}

/**
 * Calculate the pseudo inverse of an positive semidefinite matrix using the eigenvalue decomposition
 * @input a Eigen matrix to invert
 * @input eps Numerical epsilon to determine zero singular values (defaults to std::numeric_limits::eps)
 * @return Eigen matrix of same type as a
 */
template<typename T>
T posdef_pinv(const T &a, double eps = std::numeric_limits<typename T::Scalar>::epsilon()) {

  Eigen::SelfAdjointEigenSolver<T> eig(a);
  if (eig.info() != Eigen::Success) {
    return T();
  }

  typename T::Scalar tolerance = eps * std::max(a.cols(), a.rows()) * eig.eigenvalues().array().abs().maxCoeff();

  T result = eig.eigenvectors() *
             T( (eig.eigenvalues().array() > tolerance).select(eig.eigenvalues().array().inverse(), 0) ).asDiagonal() *
             eig.eigenvectors().transpose();

  return result;
}


#ifdef NDEBUG
// Release mode
// remove requirements in inner loops for speed, but keep standard require functional
#define requireDebug(req,msg)
#define require(req,msg) if (!(req)) {fputs(msg, stderr);fputs("\n\n", stderr); exit(1);}
#else
// Debug mode
// cause a crash to allow backtracing
#define requireDebug(req,msg) if (!(req)) {fputs(msg, stderr);fputs("\n\n", stderr); abort();}
#define require(req,msg) requireDebug(req,msg)
#endif

}
