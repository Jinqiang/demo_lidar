/**
 * @file OptimizationInterface.h
 * @brief Abstract base class for nonlinear optimizer.
 * @author Michael Kaess
 * @author David Rosen
 * @version $Id: OptimizationInterface.h 6371 2012-03-29 22:22:23Z kaess $
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
#include "Node.h"


namespace isam {

/**
 * Abstract base class providing an interface between the nonlinear system
 * to be optimized (stored in the Nodes of the Graph constructed in the SLAM)
 * and the Optimization class that actually performs the optimizations.
 */
class OptimizationInterface {

protected:

  /** Factored Jacobian about the current linearization point.*/
  SparseSystem _R;

public:
  virtual SparseSystem jacobian() = 0;
  virtual void apply_exmap(const Eigen::VectorXd& delta) = 0;
  virtual void self_exmap(const Eigen::VectorXd& delta) = 0;
  virtual void estimate_to_linpoint() = 0;
  virtual void linpoint_to_estimate() = 0;
  virtual void swap_estimates() = 0;
  virtual Eigen::VectorXd weighted_errors(Selector s = ESTIMATE) = 0;

  OptimizationInterface(): _R(1,1) {}

  virtual ~OptimizationInterface() {}

  friend class Optimizer;
};

}
