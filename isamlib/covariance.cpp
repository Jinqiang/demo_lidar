/**
 * @file covariance.cpp
 * @brief Recovery of marginal covariance matrix, for details see Kaess09ras.
 * @author Michael Kaess
 * @author Nicholas Carlevaris-Bianco
 * @version $Id: covariance.cpp 8578 2013-07-01 00:28:49Z kaess $
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
#include <utility> // pair

#include "isam/covariance.h"
#include "isam/util.h"

using namespace std;
using namespace Eigen;

namespace isam {

void prepare(const SparseMatrix& R, CovarianceCache& cache) {
  // reset hash table
  cache.entries.clear();
  // speedup row access
  int n = R.num_cols();
  cache.rows.resize(n);
  cache.rows_valid.resize(n);
  cache.current_valid++; // invalidate previous entries
  // wrapped back to 0? then we do have to reset the table
  if (cache.current_valid==0) {
    for (int i=0; i<n; i++) {
      cache.rows_valid[i] = 0;
    }
    cache.current_valid = 1;
  }
  // precalculate inverses
  cache.diag.resize(n);
  for (int i=0; i<n; i++) {
    cache.diag[i] = 1. / R(i,i);
  }
  // stats
  cache.num_calc = 0;
}

const SparseVector& get_row(const SparseMatrix& R, CovarianceCache& cache, int i) {
  if (cache.rows_valid[i] != cache.current_valid) {
    // retrieve and store row
    cache.rows[i] = R.get_row(i);
    cache.rows_valid[i] = cache.current_valid;
  }
  return cache.rows[i];
}

// for recursive call
double recover(const SparseMatrix& R, CovarianceCache& cache, int n, int i, int l);

inline double sum_j(const SparseMatrix& R, CovarianceCache& cache, int n, int l, int i) {
  double sum = 0;
  for (SparseVectorIter iter(get_row(R, cache, i)); iter.valid(); iter.next()) {
    double rij;
    int j = iter.get(rij);
    if (j!=i) {
      double lj;
      if (j>l) {
        lj = recover(R, cache, n, l, j);
      } else {
        lj = recover(R, cache, n, j, l);
      }
      sum += rij * lj;
    }
  }
  return sum;
}

double recover(const SparseMatrix& R, CovarianceCache& cache, int n, int i, int l) {
  if (i>l) {int tmp=i; i=l; l=tmp;}
  int id = i*n + l; // flatten index for hash table
  umap::iterator it = cache.entries.find(id);
  double res;
  if (it == cache.entries.end()) {
    //printf("calculating entry %i,%i\n", i, l);
    // need to calculate entry
    if (i==l) {
      res = cache.diag[l] * (cache.diag[l] - sum_j(R, cache, n, l, l));
    } else {
      res = -sum_j(R, cache, n, l, i) * cache.diag[i];
    }
    // insert into hash
    pair<int, double> entry(id, res);
    cache.entries.insert(entry);
    cache.num_calc++;
  } else {
    //printf("retrieved entry %i,%i\n", i, l);
    // retrieve value from hash
    res = it->second;
  }
  return res;
}

list<MatrixXd> cov_marginal(const SparseMatrix& R, CovarianceCache& cache,
                            const index_lists_t& index_lists, bool debug, int step) {
  prepare(R, cache);
  list<MatrixXd> Cs;

  // debugging
  int requested = 0;
  double t0 = tic();

  for (unsigned int i=0; i<index_lists.size(); i++) {

    const vector<int>& indices = index_lists[i];
    unsigned int n_indices = indices.size();
    MatrixXd C(n_indices, n_indices);
    // recover entries of marginal cov
    int n = R.num_cols();
    for (int r=n_indices-1; r>=0; r--) {
      for (int c=n_indices-1; c>=r; c--) {
        C(r,c) = recover(R, cache, n, indices[r], indices[c]);
      }
    }
    for (unsigned int r=1; r<n_indices; r++) {
      for (unsigned int c=0; c<r; c++) {
        C(r,c) = C(c,r);
      }
    }
    Cs.push_back(C);
    requested += indices.size()*indices.size();
  }

  if (debug) {
    double time = toc(t0);
    // timing
    printf("cov: %i calculated for %i requested in %f seconds\n",
           cache.num_calc, requested, time);
    if (step>=0) {
      // stats for gnuplot
      printf("ggg %i %i %i %i %i %i %f ",
             step,
             R.num_cols(), // side length
             R.num_cols()*R.num_cols(), // #entries dense
             R.nnz(), // #entries sparse
             cache.num_calc, // #entries calculated
             requested, // #entries requested
             time); // #execution time
    }
  }

  return Cs;
}

list<double> cov_marginal(const SparseMatrix& R, CovarianceCache& cache,
                          const entry_list_t& entry_list) {
  prepare(R, cache);
  list<double> entries;

  int n = R.num_cols();
  for (unsigned int i=0; i<entry_list.size(); i++) {
    const pair<int, int>& index = entry_list[i];
    double entry = recover(R, cache, n, index.first, index.second);
    entries.push_back(entry);
  }

  return entries;
}

}
