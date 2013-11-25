/**
 * @file numericalDiff.h
 * @brief Numerical differentiation.
 * @author Michael Kaess
 * @version $Id: numericalDiff.h 4038 2011-02-26 04:31:00Z kaess $
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

#include <vector>
#include <Eigen/Dense>

#include "isam/Node.h"

namespace isam {

// Abstract class to enforce interface for function object.
class Function {
public:
  virtual ~Function() {}
  virtual int num_measurements() const = 0;
  virtual Eigen::VectorXd evaluate() const = 0;
  virtual std::vector<Node*>& nodes() = 0;
};

/**
 * Takes a general vector valued function and returns the
 * Jacobian at the linearization point given by x0.
 * @param func Function object with evaluation function that takes and returns vectors.
 * @return Matrix containing the Jacobian of func, with
 *         dim(y) rows and dim(x) columns, where y=func(x).
 */
Eigen::MatrixXd numericalDiff(Function& func);

}
