/**
 * @file numericalDiff.cpp
 * @brief Numerical differentiation.
 * @author Michael Kaess
 * @version $Id: numericalDiff.cpp 4038 2011-02-26 04:31:00Z kaess $
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

#include <iostream>

#include "isam/numericalDiff.h"

#define SYMMETRIC

const double epsilon = 0.0001;

using namespace std;
using namespace Eigen;

namespace isam {

MatrixXd numericalDiff(Function& func) {
#ifndef SYMMETRIC
  VectorXd y0 = func.evaluate();
#endif
  // number of measurement rows
  int m = func.num_measurements();
  // number of variables
  int n = 0;
  vector<Node*>& nodes = func.nodes();
  for (vector<Node*>::iterator it = nodes.begin(); it!=nodes.end(); it++) {
    n += (*it)->dim();
  }
  // result has one column per variable
  MatrixXd Jacobian(m,n);
  int col = 0;
  // for each node...
  for (vector<Node*>::iterator it = nodes.begin(); it!=nodes.end(); it++) {
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
      node->self_exmap(delta);
      VectorXd y_plus = func.evaluate();
      node->update0(original);
#ifdef SYMMETRIC
      // evaluate negative delta
      delta(j) = -epsilon;
      node->self_exmap(delta);
      VectorXd y_minus = func.evaluate();
      node->update0(original);
      // store column
      VectorXd diff = (y_plus - y_minus) / (epsilon + epsilon);
#else
      VectorXd diff = (y_plus - y0) / epsilon;
#endif
      Jacobian.col(col) = diff;
    }
  }

  return Jacobian;
}

}
