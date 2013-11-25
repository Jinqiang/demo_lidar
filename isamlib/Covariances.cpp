/**
 * @file Covariances.cpp
 * @brief Providing access to entries of the covariance.
 * @author Michael Kaess
 * @version $Id: Covariances.cpp 7858 2013-01-14 03:50:24Z kaess $
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
#include <list>
#include <Eigen/Dense>

#include "isam/covariance.h"
#include "isam/Slam.h"

#include "isam/Covariances.h"

using namespace std;
using namespace Eigen;

namespace isam {

Covariances::Covariances(Slam& slam) : _slam(NULL), _R(slam._R) {
  // update pointers into matrix before making a copy
  slam.update_starts();
  const list<Node*>& nodes = slam.get_nodes();
  for (list<Node*>::const_iterator it = nodes.begin(); it!=nodes.end(); it++) {
    Node* node = *it;
    int start = node->_start;
    int dim = node->dim();
    _index[node] = make_pair(start, dim);
  }  
}

int Covariances::get_start(Node* node) const {
  if (_slam) {
    return node->_start;
  } else {
    return _index.find(node)->second.first;
  }
}

int Covariances::get_dim(Node* node) const {
  if (_slam) {
    return node->_dim;
  } else {
    return _index.find(node)->second.second;
  }
}

list<MatrixXd> Covariances::marginal(const node_lists_t& node_lists) const {
  const SparseSystem& R = (_slam==NULL) ? _R : _slam->_R;
  if (_slam) {
    _slam->update_starts();
  }

  vector < vector<int> > index_lists(node_lists.size());

  if (R.num_rows()>1) { // skip if _R not calculated yet (eg. before batch step)
    int n=0;
    const int* trans = R.a_to_r();
    for (node_lists_t::const_iterator it_list = node_lists.begin();
        it_list != node_lists.end();
        it_list++, n++) {
      // assemble list of indices
      vector<int> indices;
      for (list<Node*>::const_iterator it = it_list->begin(); it!=it_list->end(); it++) {
        Node* node = *it;
        int start = get_start(node);
        int dim = get_dim(node);
        for (int i=0; i<dim; i++) {
          index_lists[n].push_back(trans[start+i]);
        }
      }
    }
    return cov_marginal(R, _cache, index_lists);
  }
  list<MatrixXd> empty_list;
  return empty_list;
}

MatrixXd Covariances::marginal(const list<Node*>& nodes) const {
  node_lists_t node_lists;
  node_lists.push_back(nodes);
  return marginal(node_lists).front();
}

list<MatrixXd> Covariances::access(const node_pair_list_t& node_pair_list) const {
  const SparseSystem& R = (_slam==NULL) ? _R : _slam->_R;
  if (_slam) {
    _slam->update_starts();
  }

  if (R.num_rows()>1) { // skip if _R not calculated yet (eg. before batch step)

    // count how many entries in requested blocks
    int num = 0;
    for (node_pair_list_t::const_iterator it = node_pair_list.begin();
         it != node_pair_list.end(); it++) {
      num += get_dim(it->first) * get_dim(it->second);
    }

    // request individual covariance entries
    vector < pair<int, int> > index_list(num);
    const int* trans = R.a_to_r();
    int n = 0;
    for (node_pair_list_t::const_iterator it = node_pair_list.begin();
         it != node_pair_list.end(); it++) {
      int start_r = get_start(it->first);
      int start_c = get_start(it->second);
      int dim_r = get_dim(it->first);
      int dim_c = get_dim(it->second);
      for (int r=0; r<dim_r; r++) {
        for (int c=0; c<dim_c; c++) {
          index_list[n] = make_pair(trans[start_r + r], trans[start_c + c]);
          n++;
        }
      }
    }
    list<double> covs = cov_marginal(R, _cache, index_list);

    // assemble into block matrices
    list<MatrixXd> result;
    list<double>::iterator it_cov = covs.begin();
    for (node_pair_list_t::const_iterator it = node_pair_list.begin();
         it != node_pair_list.end(); it++) {
      int dim_r = get_dim(it->first);
      int dim_c = get_dim(it->second);
      MatrixXd matrix(dim_r, dim_c);
      for (int r=0; r<dim_r; r++) {
        for (int c=0; c<dim_c; c++) {
          matrix(r,c) = *it_cov;
          it_cov++;
        }
      }
      result.push_back(matrix);
    }
    return result;
  }
  list<MatrixXd> empty_list;
  return empty_list;
}

}
