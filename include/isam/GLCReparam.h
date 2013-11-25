/**
 * @file GLCReparam.h
 * @brief Generic Linear Constraint Reparametrization
 * @author Nicholas Carlevaris-Bianco
 * @author Michael Kaess
 * @version $Id: GLCReparam.h 8578 2013-07-01 00:28:49Z kaess $
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
#include <sstream>
#include <Eigen/Dense>

#include "Node.h"
#include "Factor.h"
#include "Anchor.h"

namespace isam {

/**
 * interface class to allow for reparametrization by the user
 * see GLC_RootShift implementation as an example
 */
class GLC_Reparam {

public:

  /**
   * set_nodes()
   * @param nodes Vector of Node pointers which support the factor
   */
  virtual void set_nodes (std::vector<Node*> nodes) = 0;
  /**
   * clone()
   * @return a new instance of the derived class
   */
  virtual GLC_Reparam* clone() = 0;
  /**
   * is_angle()
   * @return boolean vector indacating if elemetes reparameterized vector are angles
   */
  virtual Eigen::VectorXb is_angle() = 0;
  /**
   * jacobian()
   * @return jacobian of reparameterized output wrt node input
   */
  virtual Eigen::MatrixXd jacobian() = 0;
  /**
   * reparameterize()
   * @return reparameterized vector
   */
  virtual Eigen::VectorXd reparameterize (Selector s) = 0;

};

/**
 * GLC_RootShift
 * implementation of GLC_Reparam
 * performs root shift operation
 */
class GLC_RootShift : public GLC_Reparam {

private:
  std::vector<Node*> _nodes;
  Eigen::VectorXb _is_angle;
  int _dim; // input and output dim

public:

  void set_nodes (std::vector<Node*> nodes) {
    _nodes = nodes;
    _dim = 0;
    for (size_t i=0; i<_nodes.size(); i++) {
      _dim += _nodes[i]->dim();
    }
    _is_angle.resize(_dim);
    int ioff = 0;
    for (size_t i=0; i<_nodes.size(); i++) {
      _is_angle.segment(ioff, _nodes[i]->dim()) = _nodes[i]->is_angle();
      ioff += _nodes[i]->dim();
    }
  }
  
  GLC_RootShift* clone() {return new GLC_RootShift(*this);}
  Eigen::VectorXb is_angle() {return _is_angle;}
  Eigen::MatrixXd jacobian();
  Eigen::VectorXd reparameterize (Selector s);
  Eigen::VectorXd root_shift (Node* node_i, Node* node_j, Selector s);

};

} // namespace isam
