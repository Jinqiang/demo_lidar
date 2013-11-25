/**
 * @file ChowLiuTree.cpp
 * @brief ChowLiuTree for information form Gaussian distributions
 * @author Michael Kaess
 * @author Nicholas Carlevaris-Bianco
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

#include "isam/ChowLiuTree.h"
#include "isam/util.h"

using namespace std;
using namespace Eigen;

namespace isam {

MatrixXd matslice (const MatrixXd A, vector<int> ii, vector<int> jj) {
  MatrixXd B(ii.size(), jj.size());
  for (size_t i=0; i<ii.size(); i++) {
    for (size_t j=0; j<jj.size(); j++) {
      B(i,j) = A(ii[i], jj[j]);
    }
  }
  return B;
}


MatrixXd ChowLiuTreeInfo::marginal(int id) {

  // get the indicies associated with the node
  vector<int> iia(0);
  vector<int> iib(0);
  int off = 0;
  for (size_t i=0; i<_nodes.size(); i++) {
    if ((int)i == id)
      for (int j=0; j<_nodes[i]->dim(); j++) iia.push_back(off+j);
    else
      for (int j=0; j<_nodes[i]->dim(); j++) iib.push_back(off+j);
    off += _nodes[i]->dim();
  }

  if (iib.size() > 0) {
      MatrixXd Laa = matslice(_L, iia, iia);
      MatrixXd Lbb = matslice(_L, iib, iib);
      MatrixXd Lab = matslice(_L, iia, iib);
      MatrixXd Lbbinv = posdef_pinv(Lbb);
      return Laa - Lab * Lbbinv * Lab.transpose();
      //return Laa - Lab * Lbb.inverse() * Lab.transpose();
  } else
      return _L;
}

MatrixXd ChowLiuTreeInfo::joint(int ida, int idb) {

  vector<int> iia(0);
  vector<int> iib(0);
  vector<int> iic(0);
  int off = 0;
  for (size_t i=0; i<_nodes.size(); i++) {
    if ((int)i == ida)
      for (int j=0; j<_nodes[i]->dim(); j++) iia.push_back(off+j);
    else if ((int)i == idb)
      for (int j=0; j<_nodes[i]->dim(); j++) iib.push_back(off+j);
    else
      for (int j=0; j<_nodes[i]->dim(); j++) iic.push_back(off+j);
    off += _nodes[i]->dim();
  }
  vector<int> iiab(0);
  iiab.insert(iiab.end(), iia.begin(), iia.end());
  iiab.insert(iiab.end(), iib.begin(), iib.end());
  
  if (iic.size() > 0) {
    MatrixXd Labab = matslice(_L, iiab, iiab);
    MatrixXd Lcc = matslice(_L, iic, iic);
    MatrixXd Labc = matslice(_L, iiab, iic);
    MatrixXd Lccinv = posdef_pinv(Lcc);
    return Labab - Labc * Lccinv * Labc.transpose();
    //return Labab - Labc * Lcc.inverse() * Labc.transpose();
  } else
      return matslice(_L, iiab, iiab);
}

MatrixXd ChowLiuTreeInfo::conditional(int ida, int idb) {

  MatrixXd Lj = joint(ida, idb);
  return Lj.block(0, 0, _nodes[ida]->dim(), _nodes[ida]->dim());

}

ChowLiuTree::ChowLiuTree (const Eigen::MatrixXd &L, const std::vector<Node *>& nodes)
  : _clt_info(L, nodes)
{
  // make sure we have at least two nodes otherwise return a trival tree
  if (nodes.size() == 1) {

    ChowLiuTreeNode node;
    node.id = 0;
    node.pid = -1;
    node.marginal = _clt_info.marginal(0);
    node.conditional = node.marginal;
    node.joint = node.marginal;
    tree[node.id] = node;

  } else {

    //calculate the parent nodes based on maximising mutual information
    _calc_edges();
    _max_span_tree();
    tree.clear();
    _build_tree_rec(_edges.front().id1, -1);

  }
}

bool mi_sort(MI &first, MI &second) {
  return first.mi > second.mi;
}

void
ChowLiuTree::_calc_edges()  {

  int nn = _clt_info.num_nodes();
  //double npairs = floor(pow((double)bb, 2.0) / 2) -  floor((double)nn/2);

  for (int i=0; i<nn; i++) {
    for (int j=(i+1); j<nn; j++) {
      MI mi_tmp (i,j, _calc_mi(i, j));
      _edges.push_back(mi_tmp);
    }
  }
  _edges.sort(mi_sort);
  
}

double ChowLiuTree::_calc_mi(int ida, int idb) {

  MatrixXd L_agb = _clt_info.conditional (ida, idb);
  MatrixXd L_a = _clt_info.marginal (ida);

  // use pdet
  //double ldL_agb = plogdet(L_agb);
  //double ldL_a = plogdet(L_a);
  //double mi = 0.5*(ldL_agb - ldL_a);

  // use normal det, must be pinned
  double ldL_agb = log ((L_agb +  MatrixXd::Identity(L_agb.rows(), L_agb.cols())).determinant());
  double ldL_a = log ((L_a +  MatrixXd::Identity(L_a.rows(), L_a.cols())).determinant());
  double mi = 0.5*(ldL_agb - ldL_a);

  return mi;
}

void ChowLiuTree::_max_span_tree() {

  // init groups: assign each id to a different group initially
  map<int, int> groups; // map <node index, group index>
  for (int i=0; i<_clt_info.num_nodes(); i++) {
      groups[i] = i;
  }
    
  int group1, group2;
  list<MI>::iterator edge = _edges.begin();
  while(edge != _edges.end()) {
    if(groups[edge->id1] != groups[edge->id2]) {
      group1 = groups[edge->id1];
      group2 = groups[edge->id2];

      // merge group2 into group 1
      map<int, int>::iterator groupIt;
      for(groupIt = groups.begin(); groupIt != groups.end(); groupIt++)
        if(groupIt->second == group2) groupIt->second = group1;

      edge++;
    } else {
      edge = _edges.erase(edge);
    }
  }
  
}

void ChowLiuTree::_build_tree_rec(int id, int pid) {

  ChowLiuTreeNode new_node;

  new_node.id = id;
  new_node.pid = pid;
  new_node.marginal = _clt_info.marginal(id);
  if (pid == -1) {
    new_node.conditional = new_node.marginal;
    new_node.joint = new_node.marginal;
  } else {
    new_node.conditional = _clt_info.conditional(id, pid);
    new_node.joint = _clt_info.joint(id, pid);
  }

  vector<int> cids; // vector of child ids
  list<MI>::iterator edge = _edges.begin();
  while(edge != _edges.end()) {
    if(edge->id1 == new_node.id) {
      cids.push_back(edge->id2);
      edge = _edges.erase(edge);
      continue;
    }
    if(edge->id2 == new_node.id) {
      cids.push_back(edge->id1);
      edge = _edges.erase(edge);
      continue;
    }
    edge++;
  }
  for(size_t i=0; i < cids.size(); i++) {
      new_node.cids.push_back(cids[i]);
      _build_tree_rec(cids[i], new_node.id);
  }

  tree[new_node.id] = new_node;
}

} //namespace isam