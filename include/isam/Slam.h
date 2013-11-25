/**
* @file Slam.h
* @brief SLAM implementation using iSAM
* @author Michael Kaess
* @version $Id: Slam.h 6371 2012-03-29 22:22:23Z kaess $
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

#include <string>
#include <list>
#include <Eigen/Dense>

#include "SparseSystem.h"
#include "Node.h"
#include "Factor.h"
#include "Graph.h"
#include "Properties.h"
#include "OptimizationInterface.h"
#include "Optimizer.h"
#include "Covariances.h"


namespace isam {


/**
* Return type of Slam::update() to allow future extensions without
* having to change the interface.
*/
class UpdateStats {
public:
  // current step number
  int step;

  // was batch performed?
  bool batch;

  // was the solution updated?
  bool solve;
};

/**
* The actual SLAM interface.
*/
class Slam: public Graph, OptimizationInterface {
  // Graph prohibits copy construction and assignment operator

  int _step;

  Properties _prop;

  Covariances _covariances;

public:

  //-- manipulating the graph -----------------------------

  /**
  * Default constructor.
  */
  Slam();

  /**
  * Destructor.
  */
  virtual ~Slam();

  /**
  * Returns a copy of the current properties.
  */
  Properties properties() {
    return _prop;
  }

  /**
  * Sets new properties.
  */
  void set_properties(Properties prop) {
    _prop = prop;
  }

  /**
  * Saves the graph (nodes and factors).
  * @param fname Filename with optional path to save graph to.
  */
  void save(const std::string fname) const;

  /**
  * Adds a node (variable) to the graph.
  * @param node Pointer to new node.
  */
  void add_node(Node* node);

  /**
  * Adds a factor (measurement) to the graph.
  * @param factor Pointer to new factor.
  */
  void add_factor(Factor* factor);

  /**
  * Removes a node (variable) and all adjacent factors from the graph.
  * Note that the node itself is not deallocated.
  * @param node Pointer to node.
  */
  void remove_node(Node* node);

  /**
  * Removes an factor (measurement) from the graph.
  * Note that the factor itself is not deallocated.
  * Be careful not to leave unconnected nodes behind.
  * @param factor Pointer to factor.
  */
  void remove_factor(Factor* factor);

  //-- solving the system -----------------------------

  /**
  * Update the graph by finding new solution; depending on properties
  * this might simply be a Givens update, could include a solve step,
  * or be a full batch step with reordering.
  * @return Update statistics.
  */
  virtual UpdateStats update();

  /**
  * Fully solve the system, iterating until convergence.
  * @return Number of iterations performed.
  */
  virtual int batch_optimization();

  //-- misc -----------------------------

  /**
  * Sets a cost function different from the default (quadratic).
  * @param cost_func Pointer to cost function, see util.h for a list of robust
  * cost functions. Instead of cost_squared, use NULL, which avoids calculating square roots.
  */
  void set_cost_function(cost_func_t cost_func);

  /**
  * Calculates the normalized chi-square value (weighted sum of squared
  * errors divided by degrees of freedom [# measurements - # variables])
  * for the estimate x.
  */
  double normalized_chi2();

  /**
  * Calculates the chi2 error of the last_n constraints.
  */
  double local_chi2(int last_n);

  /**
  * Weighted non-squared error vector, by default at current estimate.
  */
  Eigen::VectorXd weighted_errors(Selector s = ESTIMATE);

  /**
  * Weighted sum of squared errors, by default at the current estimate.
  */
  double chi2(Selector s = ESTIMATE);

  /**
  * Returns the current factor matrix.
  */
  virtual const SparseSystem& get_R() const;

  /**
  * Calculate the full Jacobian numerical (fast column-wise procedure).
  */
  virtual SparseSystem jacobian_numerical_columnwise();

  /**
  * Returns the last n rows of the measurement Jacobian of the SLAM system.
  * @param last_n Only return Jacobians of last n measurements (-1 returns all)
  * @return Measurement Jacobian.
  */
  virtual SparseSystem jacobian_partial(int last_n);

  /**
  * Returns the measurement Jacobian of the SLAM system.
  * @return Measurement Jacobian.
  */
  virtual SparseSystem jacobian();

  /**
   * Returns the Covariances object associated with this Slam object.
   * @return Covariances object for access to estimation covariances.
   */
  const Covariances& covariances();

  /**
  * Print statistics for debugging.
  */
  virtual void print_stats();

private:

  /**
  * Apply a delta vector to the linearization point and store result as new estimate.
  */
  void apply_exmap(const Eigen::VectorXd& x);

  /**
  * Apply a delta vector directly to the linearization point.
  */
  void self_exmap(const Eigen::VectorXd& x);

  /**
  * Set current estimate to the linearization point (needed in dog leg).
  */
  void linpoint_to_estimate();

  /**
  * Set linearization point to current estimate.
  */
  void estimate_to_linpoint();

  /**
  * Exchange linearization point and current estimate (needed in dog leg).
  */
  void swap_estimates();

  /**
  * Update the system with any newly added measurements. The measurements will be
  * appended to the existing factor matrix, and the factor is transformed into
  * triangular form again using Givens rotations.
  * Very efficient for exploration O(1), but can be more expensive otherwise >O(n).
  */
  virtual void incremental_update();

  /**
  * Resolve the system with linearization based on current estimate;
  * perform variable reordering for efficiency.
  */
  virtual void batch_optimization_step();


  // internal variable used for operations such as removing of parts of
  // the graph that currently cannot be done incrementally
  bool _require_batch;

  cost_func_t _cost_func;

  void update_starts();

protected:
  int _dim_nodes;
  int _dim_measure;
  int _num_new_measurements;
  int _num_new_rows;

  Optimizer _opt;

  friend class Covariances;

};

}
