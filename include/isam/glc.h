/**
 * @file glc.h
 * @brief Generic Linear Constraint functions
 * @author Nicholas Carlevaris-Bianco
 * @author Michael Kaess
 * @version $Id: glc.h 8578 2013-07-01 00:28:49Z kaess $
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
#include "GLCReparam.h"

//#define GLC_DEBUG

namespace isam {

/**
 * Removes a node (variable) and all factors within in its elimination clique
 * replaces these factors with either a dense GLC factor or a set of sparse
 * GLC factors
 * Note that the node itself and the removed factors are not deallocated.
 *
 * We usually use GLC with USE_QUATERNIONS turned off in Rot3d.h. This
 * avoids having to account for the state transformation between euler and
 * angle axis representations. This is currently done numerically
 * (see exmap_jacobian()) and produces slightly higher KLD in the GLC results.
 *
 * @param node Pointer to node to be removed.
 * @param sparse Bool flag if new factors should be sparse approximate or dense
 * @param rp functor to reparamertize variables before linearization
 * @return vector of new factors.
 */
std::vector<Factor*> glc_remove_node(Slam& slam, Node* node, bool sparse = false,
                                     GLC_Reparam* rp = NULL);

/**
 * Find the factors that will be eliminated and replaced when node is removed
 * IMPORTANT: must be called before node is removed using glc_remove_node()
 * @param node Pointer to node to be removed.
 * @return vector of factors that will be eliminated.
 */
std::vector<Factor*> glc_elim_factors(Node* node);

#ifdef USE_QUATERNIONS
Eigen::MatrixXd exmap_jacobian (const std::vector<Node*>& nodes);
#endif

/**
 * Generic Linear Constraint
 */
class GLC_Factor : public Factor {

public:
  const Eigen::VectorXd _x;
  const Eigen::MatrixXd _G;
  GLC_Reparam* _rp;

  /**
   * Constructor.
   * @param poses Vector of Pose3d_Node pointers which support the factor
   * @param x The state vector used as the linearization point for the GLC
   * @param G The PCA transformation calculated from the target information
   * @param rp Reparametrization object, uses to implement root-shift or user defined reparameterizations
   *           to be applied before commiting to a linearization point
   */
  GLC_Factor(std::vector<Node*> nodes, const Eigen::VectorXd& x, const Eigen::MatrixXd& G, GLC_Reparam* rp = NULL)
    : Factor("GLC_Factor", G.rows(), Information(Eigen::MatrixXd::Identity(G.rows(), G.rows()))),
    _x(x), _G(G)
  {
    _nodes.resize(nodes.size());
    for (size_t i=0; i<nodes.size(); i++) {
      _nodes[i] = nodes[i];
    }
    _rp = NULL;
    if (rp != NULL) {
      _rp = rp->clone();
      _rp->set_nodes(_nodes);
    }

  }

  ~GLC_Factor () {
    if (_rp != NULL)
      delete _rp;
  }

  void initialize() {  
  }
  
  Jacobian jacobian () {

    Eigen::MatrixXd H;
    if (_rp) {
      H = sqrtinf() * _G * _rp->jacobian();
    }else {
      H = sqrtinf() * _G;
    }

#ifdef USE_QUATERNIONS
    Eigen::MatrixXd Jexmap = exmap_jacobian (_nodes);
    H = H * Jexmap;
#endif

    Eigen::VectorXd r = error(LINPOINT);
    Jacobian jac(r);
    int position = 0;
    int n_measure = dim();
    for (size_t i=0; i<_nodes.size(); i++) {
      int n_var = _nodes[i]->dim();
      Eigen::MatrixXd Hi = H.block(0, position, n_measure, n_var);
      position += n_var;
      jac.add_term(_nodes[i], Hi);
    }

#ifdef GLC_DEBUG
    // compare with numerical jacobian
    Jacobian njac = Factor::jacobian();
    Eigen::MatrixXd nH(n_measure, H.cols());
    int offset = 0;
    for (Terms::const_iterator it=njac.terms().begin(); it!=njac.terms().end(); it++) {
        int ncols = it->term().cols();
        nH.block(0, offset, n_measure, ncols) = it->term();
        offset += ncols;
    }
    if ((nH - H).array().abs().matrix().lpNorm<Eigen::Infinity>() > 1e-6) {
      std::cout << "Ja = [" << H << "];" << std::endl;
      std::cout << "Jn = [" << nH << "];" << std::endl;
    }
#endif
    return jac;
  }
  
  Eigen::VectorXd basic_error(Selector s = LINPOINT) const {
    int nn = _nodes.size();
    Eigen::VectorXd x_p;
    Eigen::VectorXb is_angle;
    
    if (_rp) {
      x_p = _rp->reparameterize(s);
      is_angle = _rp->is_angle();
    } else {
      x_p.resize(_G.cols());
      is_angle.resize(_G.cols());
      int off = 0;
      for (int i=0; i<nn; i++) {
        int d = _nodes[i]->dim();
        x_p.segment(off, d) = _nodes[i]->vector(s);
        is_angle.segment(off, d) = _nodes[i]->is_angle();
        off += _nodes[i]->dim();
      }
    }
    
    Eigen::VectorXd err_fs =  x_p - _x;
    for (int i=0; i<err_fs.size(); i++) {
      if (is_angle(i))
        err_fs(i) = standardRad(err_fs(i));
    }
    
    Eigen::VectorXd err = _G*err_fs;
    return err;
  }
  
  void write(std::ostream &out) const {
    Factor::write(out);
    out << " (";
    for (int i=0; i<_x.size(); i++)
        out << _x(i) << " ";
    out << ") {";
    for (int i=0; i<_G.rows(); i++)
        for (int j=0; j<_G.cols(); j++)
            out << _G(i,j) << " ";
    out << "}";  
  }

};

} // namespace isam
