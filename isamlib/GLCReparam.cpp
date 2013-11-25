/**
 * @file GLCReparam.cpp
 * @brief Generic Linear Constraint Reparametrization
 * @author Nicholas Carlevaris-Bianco
 * @author Michael Kaess
 * @version $Id: GLCReparam.cpp 8578 2013-07-01 00:28:49Z kaess $
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

#include <vector>

#include "isam/GLCReparam.h"

#include "isam/Pose3d.h"
#include "isam/Point3d.h"
#include "isam/Pose2d.h"
#include "isam/Point2d.h"
#include "isam/slam2d.h"
#include "isam/slam3d.h"

using namespace std;
using namespace Eigen;

namespace isam {

VectorXd GLC_RootShift::reparameterize (Selector s) {

  VectorXd x;

  int np = _nodes.size();

  if (np == 1) {
    x = _nodes[0]->vector(s);
    return x;
  }

  x.resize(_dim);
  x.setZero();

  // decide on the root for the root shift
  // for now just use the first pose2d or pose3d
  int ir=0;
  for (ir=0; ir<np; ir++) {
    if (0 == strcmp(_nodes[ir]->name(), "Pose3d") ||
        0 == strcmp(_nodes[ir]->name(), "Pose2d"))
      break;
  }
  assert (ir != np); // assumes a pose3d or pose2d node is avaliable
  Node *node_i = _nodes[ir];

  int ioff = 0;
  for (int j=0; j<np; j++) {
    if (j == ir) { // root node
      x.segment(ioff, node_i->dim()) = node_i->vector(s);
      ioff += node_i->dim();
    } else {
      Node *node_j = _nodes[j];
      VectorXd x_i_j = root_shift (node_i, node_j, s);
      assert(node_j->dim() == x_i_j.size());
      x.segment(ioff, node_j->dim()) = x_i_j;
      ioff += node_j->dim();
    }
  }

  return x;
}

VectorXd GLC_RootShift::root_shift (Node* node_i, Node* node_j, Selector s) {

  VectorXd x_i_j;
  if (0 == strcmp(node_i->name(), "Pose3d")) {
    Pose3d pose3d_i = dynamic_cast<Pose3d_Node*>(node_i)->value(s);
    if (0 == strcmp(node_j->name(), "Pose3d")) {
      x_i_j = dynamic_cast<Pose3d_Node*>(node_j)->value(s).ominus(pose3d_i).vector();
    } else if (0 == strcmp(node_j->name(), "Point3d")) {
      x_i_j = pose3d_i.transform_to(dynamic_cast<Point3d_Node*>(node_j)->value(s)).vector();
    } else {
      std::cout << "Warning: root shift undefined for:" << node_i->name() << " and " << node_j->name() << std::endl;
      x_i_j = node_j->vector(s); // identity
    }
  } else if (0 == strcmp(node_i->name(), "Pose2d")) {
    Pose2d pose2d_i = dynamic_cast<Pose2d_Node*>(node_i)->value(s);
    if (0 == strcmp(node_j->name(), "Pose2d")) {
      x_i_j = dynamic_cast<Pose2d_Node*>(node_j)->value(s).ominus(pose2d_i).vector();
    } else if (0 == strcmp(node_j->name(), "Point2d")) {
      x_i_j = pose2d_i.transform_to(dynamic_cast<Point2d_Node*>(node_j)->value(s)).vector();
    } else {
      std::cout << "Warning: root shift undefined for:" << node_i->name() << " and " << node_j->name() << std::endl;
      x_i_j = node_j->vector(s); // identity
    }
  }

  return x_i_j;
}

const double epsilon = 0.0001;

// TODO replace with analytical jacobians
MatrixXd GLC_RootShift::jacobian() {

  MatrixXd J(_dim,_dim);

  int col = 0;
  // for each node...
  for (vector<Node*>::iterator it = _nodes.begin(); it!=_nodes.end(); it++) {
    Node* node = *it;
    int dim_n = node->dim();
    // for each dimension of the node...
    for (int j=0; j<dim_n; j++, col++) {
      VectorXd delta(dim_n);
      delta.setZero();
      // remember original value
      VectorXd original = node->vector0();
      // evaluate positive delta
      delta(j) = epsilon;
      //node->self_exmap(delta); // taken care of in exmap_jacobian() in glc.[h/cpp]
      node->update0(node->vector0() + delta);
      VectorXd y_plus = reparameterize(LINPOINT);
      node->update0(original);
      // evaluate negative delta
      delta(j) = -epsilon;
      //node->self_exmap(delta);  // taken care of in exmap_jacobian() in glc.[h/cpp]
      node->update0(node->vector0() + delta);
      VectorXd y_minus = reparameterize(LINPOINT);
      node->update0(original);
      // store column
      VectorXd diff = (y_plus - y_minus) / (epsilon + epsilon);
      // wrap angular difference
      for (int k=0; k<_dim; k++) {
        if (_is_angle(k))
          diff(k) = standardRad(diff(k));
      }
      J.col(col) = diff;
    }
  }

  return J;
}

} // namespace isam
