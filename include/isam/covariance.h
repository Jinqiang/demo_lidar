/**
 * @file covariance.h
 * @brief Recovery of marginal covariance matrix.
 * @author Michael Kaess
 * @author Nicholas Carlevaris-Bianco
 * @version $Id: covariance.h 8578 2013-07-01 00:28:49Z kaess $
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
#include <utility> // pair
#include <list>
#include <Eigen/Dense>

#define USE_TR1
#ifdef USE_TR1
#include <tr1/unordered_map>
#else
#include "boost/unordered_map.hpp"
#endif

#include "SparseMatrix.h"

namespace isam {

#ifdef USE_TR1
#include <tr1/unordered_map>
typedef std::tr1::unordered_map<int, double> umap;
#else
typedef boost::unordered_map<int, double> umap;
#endif

class CovarianceCache {
public:
  umap entries;
  // precalculated diagonal inverses
  std::vector<double> diag;
  // recovering rows is expensive, buffer results
  std::vector<SparseVector> rows;
  // avoid having to cleanup buffers each time by explicitly marking entries as valid
  std::vector<unsigned int> rows_valid;
  // avoid having to cleanup valid entries by using different indices each time
  unsigned int current_valid;
  // stats
  int num_calc;

  CovarianceCache () {
    current_valid = 1;
  }
};

typedef std::vector< std::vector<int> > index_lists_t;
typedef std::vector< std::pair<int, int> > entry_list_t;

/**
 * Takes a list of variable indices, and returns the marginal covariance matrix.
 * @param R Sparse factor matrix.
 * @param cache Covariance cache object.
 * @param index_lists List of lists of indices; a block will be recovered for each list.
 * @param debug Optional parameter to print timing information.
 * @param step Optional parameter to print statistics (default=-1, no stats printed).
 * @return List of dense marginal covariance matrices.
 */
std::list<Eigen::MatrixXd> cov_marginal(const SparseMatrix& R, CovarianceCache& cache,
                                        const index_lists_t& index_lists,
                                        bool debug=false, int step=-1);

/**
 * Takes a list of pairs of integers and returns the corresonding
 * entries of the covariance matrix.
 * @param R Sparse factor matrix.
 * @param cache Covariance cache object.
 * @param entry_lists List of pairs of integers refering to covariance matrix entries.
 * @param debug Optional parameter to print timing information.
 * @param step Optional parameter to print statistics (default=-1, no stats printed).
 * @return List of doubles corresponding to the requested covariance matrix entries.
 */
std::list<double> cov_marginal(const SparseMatrix& R, CovarianceCache& cache,
                               const entry_list_t& entry_list);

}
