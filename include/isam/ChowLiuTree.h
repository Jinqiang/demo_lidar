/**
 * @file ChowLiuTree.h
 * @brief ChowLiuTree for information form Gaussian distributions
 * @author Nicholas Carlevaris-Bianco
 * @author Michael Kaess
 * @version $Id: covariance.cpp 4975 2011-07-13 17:49:09Z kaess $
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
#include <vector>
#include <list>
#include <map>
#include <sstream>
#include <Eigen/Dense>

#include "Node.h"


namespace isam {


/**
 * ChowLiuTreeInfo
 * Information for Gaussian distribtuion for chow liu tree
 *
 */
class ChowLiuTreeInfo {

private:

  const Eigen::MatrixXd& _L;
  const std::vector<Node *>& _nodes;

public:

  /**
   * Constructor.
   * @param L the information matrix
   * @param nodes the nodes in the information matrix
   */
  ChowLiuTreeInfo (const Eigen::MatrixXd& L, const std::vector<Node *>& nodes)
      : _L(L), _nodes(nodes) { }

  //const std::vector<Node *>& nodes() {return _nodes;}
  int num_nodes() {return _nodes.size();}

  /**
   * Marginal distribution
   * @return Information matrix P(node[id])
   */
  Eigen::MatrixXd marginal(int id);

  /**
   * Joint distribution
   * @return Information matrix P(node[ida], node[idb])
   */
  Eigen::MatrixXd joint(int ida, int idb); // a, b

  /**
   * Conditional distribution
   * @return Information matrix P(node[ida] | node[idb])
   */
  Eigen::MatrixXd conditional(int ida, int idb);

};

/**
 * ChowLiuTreeNode
 * A node in the chowliu tree
 *
 */
class ChowLiuTreeNode {

public:
  
  int id;
  int pid;
  std::vector<int> cids;
  Eigen::MatrixXd marginal;     // marginal information
  Eigen::MatrixXd conditional;  // conditional information with parrent
  Eigen::MatrixXd joint;        // joint information with parrent
  
  bool is_root () { return (pid == -1) ? true : false; }
  bool is_leaf () { return (cids.size() == 0) ? true : false; }

};

class MI {
  
public:
  int id1;
  int id2;
  double mi;

  MI (int id1_, int id2_, double mi_)
    : id1(id1_), id2(id2_), mi(mi_){
  }

};



/**
 * ChowLiuTree
 * Chow Liu Tree class for information form gaussian distribtuions
 *
 */
class ChowLiuTree {

  ChowLiuTreeInfo _clt_info;
  std::list<MI> _edges;

  void _calc_edges();
  double _calc_mi(int ida, int idb);
  void _max_span_tree();
  void _build_tree_rec(int id, int pid);

public:

  // the resulting chow-liu tree
  std::map<int, ChowLiuTreeNode> tree;

  ChowLiuTree (const Eigen::MatrixXd &L, const std::vector<Node *>& nodes);
  
};


} // namespace isam
