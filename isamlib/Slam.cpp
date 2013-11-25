/**
 * @file Slam.cpp
 * @brief SLAM implementation using iSAM
 * @author Michael Kaess
 * @author Hordur Johannsson
 * @version $Id: Slam.cpp 7610 2012-10-25 10:21:12Z hordurj $
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

#include <iomanip>
#include <vector>
#include <map>
#include <list>

#include "isam/util.h"
#include "isam/SparseSystem.h"
#include "isam/OptimizationInterface.h"
#include "isam/covariance.h"

#include "isam/Slam.h"

using namespace std;
using namespace Eigen;

namespace isam {

// for numbering of factors and nodes
int Node::_next_id = 0;
int Factor::_next_id = 0;

struct DeleteOnReturn
{
  SparseVector** _ptr;
  explicit DeleteOnReturn(SparseVector** ptr) : _ptr(ptr) {}
  ~DeleteOnReturn() { delete [] _ptr; }
};

// for getting correct starting positions in matrix for each node,
// only needed after removing nodes
void Slam::update_starts() {
  int start = 0;
  const list<Node*>& nodes = get_nodes();
  for (list<Node*>::const_iterator it = nodes.begin(); it!=nodes.end(); it++) {
    Node* node = *it;
    node->_start = start;
    start += node->dim();
  }
}

Slam::Slam()
  : Graph(),
    _step(0), _prop(Properties()),
    _covariances(this),
    _require_batch(true), _cost_func(NULL),
    _dim_nodes(0), _dim_measure(0),
    _num_new_measurements(0), _num_new_rows(0),
    _opt(*this)
{
}

Slam::~Slam()
{
}

void Slam::save(const string fname) const {
  ofstream out(fname.c_str(), ios::out | ios::binary);
  require(out, "Slam.save: Cannot open output file.");
  write(out);
  out.close();
}

void Slam::add_node(Node* node) {
  Graph::add_node(node);
  _dim_nodes += node->dim();
}

void Slam::add_factor(Factor* factor) {
  // adds itself to factor lists of adjacent nodes; also initialized linked nodes if necessary
  factor->initialize_internal();
  // needed to change cost function
  factor->set_cost_function(&_cost_func);
  Graph::add_factor(factor);
  _num_new_measurements++;
  _num_new_rows += factor->dim();
  _dim_measure += factor->dim();
}

void Slam::remove_node(Node* node) {
  // make a copy, as the original will indirectly be modified below in remove_factor()
  list<Factor*> factors = node->factors(); 
  for (list<Factor*>::iterator factor = factors.begin(); factor!=factors.end(); factor++) {
    remove_factor(*factor);
  }
  _dim_nodes -= node->dim();
  Graph::remove_node(node);
  _require_batch = true;
}

void Slam::remove_factor(Factor* factor) {
  vector<Node*> nodes = factor->nodes();
  for (vector<Node*>::iterator node = nodes.begin(); node!=nodes.end(); node++) {
    (*node)->remove_factor(factor);
  }
  _dim_measure -= factor->dim();
  Graph::remove_factor(factor);
  _require_batch = true;
}

void Slam::incremental_update()
{
  // incremental update not possible after removing nodes or factors
  // (might change in the future)
  if (_require_batch)
  {
    batch_optimization_step();
  }
  else if (_num_new_measurements > 0)
  {
    SparseSystem jac_new = jacobian_partial(_num_new_measurements);

    _opt.augment_sparse_linear_system(jac_new, _prop);

    _num_new_measurements = 0;
    _num_new_rows = 0;
  }
}

void Slam::batch_optimization_step()
{
  _require_batch = false;
  // update linearization point x0 with current estimate x
  _num_new_measurements = 0;
  _num_new_rows = 0;

  _opt.relinearize(_prop);
}

UpdateStats Slam::update()
{
  UpdateStats stats;
  stats.batch = false;
  stats.solve = false;
  if (_step%_prop.mod_update == 0)
  {
    if (_step%_prop.mod_batch == 0)
    {
      // batch solve periodically to avoid fill-in
      if (!_prop.quiet)
      {
        cout << endl;
        cout << "step " << _step;
      }
      batch_optimization_step();
      stats.batch = true;
    }
    else
    {
      // for efficiency, incrementally update most of the time.
      if (!_prop.quiet)
      {
        cout << ".";
        fflush(stdout);
      }
      incremental_update();
      if (_step%_prop.mod_solve == 0)
      {
        stats.solve = true;

        _opt.update_estimate(_prop);
      }
    }
  }
  _step++;
  stats.step = _step;

  return stats;
}

int Slam::batch_optimization()
{
  int num_iterations = 0;

  int variables_deleted;
  int measurements_deleted;
  erase_marked(variables_deleted, measurements_deleted);
  _dim_nodes -= variables_deleted;
  _dim_measure -= measurements_deleted;

  _opt.batch_optimize(_prop, &num_iterations);
  return num_iterations;
}

void Slam::set_cost_function(cost_func_t func) {
  _cost_func = func;
}

void Slam::apply_exmap(const Eigen::VectorXd& x) {
  int pos = 0;
  for (list<Node*>::iterator node = _nodes.begin(); node != _nodes.end(); node++) {
    int dim = (*node)->dim();
    const VectorXd& xi = x.segment(pos, dim);
    (*node)->apply_exmap(xi);
    pos += dim;
  }
}

void Slam::self_exmap(const Eigen::VectorXd& x) {
  int pos = 0;
  for (list<Node*>::iterator node = _nodes.begin(); node != _nodes.end(); node++) {
    int dim = (*node)->dim();
    VectorXd xi = x.segment(pos, dim);
    (*node)->self_exmap(xi);
    pos += dim;
  }
}

void Slam::linpoint_to_estimate() {
  for (list<Node*>::iterator node = _nodes.begin(); node!=_nodes.end(); node++) {
    (*node)->linpoint_to_estimate();
  }
}

void Slam::estimate_to_linpoint() {
  for (list<Node*>::iterator node = _nodes.begin(); node!=_nodes.end(); node++) {
    (*node)->estimate_to_linpoint();
  }
}

void Slam::swap_estimates() {
  for (list<Node*>::iterator node = _nodes.begin(); node!=_nodes.end(); node++) {
    (*node)->swap_estimates();
  }
}

VectorXd Slam::weighted_errors(Selector s) {
  VectorXd werrors(_dim_measure);
  const list<Factor*>& factors = get_factors();
  int start = 0;
  for (list<Factor*>::const_iterator it = factors.begin(); it!=factors.end(); it++) {
    int dim = (*it)->dim();
    werrors.segment(start, dim) = (*it)->error(s);
    start += dim;
  }
  return werrors;
}

double Slam::chi2(Selector s) {
  return weighted_errors(s).squaredNorm();
}

double Slam::local_chi2(int last_n) {
  // avoiding two passes by allocating maximum size vector
  VectorXd werrors(_dim_measure);
  const list<Factor*>& factors = get_factors();
  int start = _dim_measure;
  int n = 0;
  for (list<Factor*>::const_reverse_iterator it = factors.rbegin();
      it!=factors.rend() && n<last_n;
      it++, n++) {
    int dim = (*it)->dim();
    start -= dim;
    werrors.segment(start, dim) = (*it)->error(ESTIMATE);
  }
  // only use actually calculated part of werrors
  return werrors.tail(_dim_measure-start).squaredNorm();
}

double Slam::normalized_chi2() {
  return chi2() / (double)(_dim_measure - _dim_nodes);
}

const SparseSystem& Slam::get_R() const {
  return _R;
}

const double epsilon = 0.0001;

SparseSystem Slam::jacobian_numerical_columnwise() {
  // label starting points of rows, and initialize sparse row vectors with
  // correct number of entries
  int num_rows = _dim_measure;
  DeleteOnReturn rows_ptr(new SparseVector* [num_rows]);
  SparseVector** rows = rows_ptr._ptr; //[num_rows];
  int pos = 0;
  vector<int> factor_offset(get_factors().size());
  for (list<Factor*>::const_iterator it = get_factors().begin();
      it!=get_factors().end();
      it++) {
    (*it)->_start = pos;
    int dimtotal = 0;
    for (vector<Node*>::const_iterator it2 = (*it)->nodes().begin();
        it2!=(*it)->nodes().end();
        it2++) {
      dimtotal += (*it2)->dim();
    }
    for (int i=0; i<(*it)->dim(); i++, pos++) {
      // do not delete, will be pulled into SparseSystem below
      rows[pos] = new SparseVector(dimtotal);
    }
  }
  // larger than needed, but avoids some book keeping and avoids many (smaller) allocations
  VectorXd y_plus(num_rows);
  VectorXd y_minus(num_rows);
  int col = 0;
  for (list<Node*>::const_iterator it = get_nodes().begin(); it!=get_nodes().end(); it++) {
    Node* node = *it;
    int dim_node = node->dim();
    VectorXd delta(dim_node);
    delta.setZero();
    VectorXd original = node->vector0();
    for (int c=0; c<dim_node; c++, col++) {
      // calculate column for +epsilon
      delta(c) = epsilon;
      node->self_exmap(delta);
      int row = 0;
      for (list<Factor*>::const_iterator it_factor = node->factors().begin();
          it_factor!=node->factors().end();
          it_factor++) {
        Factor* factor = *it_factor;
        int dim_factor = factor->dim();
        y_plus.segment(row, dim_factor) = factor->evaluate();
        row += dim_factor;
      }
      node->update0(original);
      // calculate column for -epsilon
      delta(c) = - epsilon;
      node->self_exmap(delta);
      row = 0;
      for (list<Factor*>::const_iterator it_factor = node->factors().begin();
          it_factor!=node->factors().end();
          it_factor++) {
        Factor* factor = *it_factor;
        int dim_factor = factor->dim();
        y_minus.segment(row, dim_factor) = factor->evaluate();
        row += dim_factor;
      }
      node->update0(original);
      delta(c) = 0.; // reusing delta
      // calculate derivative
      VectorXd diff = (y_plus.head(row) - y_minus.head(row)) / (epsilon + epsilon);
      // write entries into sparse Jacobian
      row = 0;
      int i = 0;
      for (list<Factor*>::const_iterator it_factor = (*it)->factors().begin();
          it_factor!=(*it)->factors().end();
          it_factor++, i++) {
        for (int r=0; r<(*it_factor)->dim(); r++, row++) {
          if (diff(row)!=0.) { // omit 0 entries
            int offset = (*it_factor)->_start;
            rows[offset+r]->append(col, diff(row)); // faster than SparseVector.set
          }
        }
      }
    }
  }
  VectorXd rhs = weighted_errors(LINPOINT);
  return SparseSystem(num_rows, _dim_nodes, rows, rhs);
}

SparseSystem Slam::jacobian() {
  if (_prop.force_numerical_jacobian) {
    // column-wise is more efficient, especially if some nodes are
    // connected to many factors
    return jacobian_numerical_columnwise();
  } else {
    // have to do row-wise if we want to use any available symbolic
    // derivatives
    return jacobian_partial(-1);
  }
}

const Covariances& Slam::covariances() {
  return _covariances;
}

SparseSystem Slam::jacobian_partial(int last_n) {
  update_starts();
  // actual assembly of Jacobian
  int num_rows = _dim_measure;
  if (last_n > 0) {
    num_rows = _num_new_rows;
  }
  DeleteOnReturn rows_ptr(new SparseVector*[num_rows]);
  SparseVector** rows = rows_ptr._ptr; //[num_rows];

  VectorXd rhs(num_rows);
  int row = 0;
  const list<Factor*>& factors = get_factors();
  list<Factor*>::const_iterator it = factors.begin();
  if (last_n != -1) {
    // skip all entries except for last_n
    for (int n = num_factors(); n>last_n; n--, it++);
  }
  for (; it!=factors.end(); it++) {
    Factor* factor = *it;
    Jacobian jac = factor->jacobian_internal(_prop.force_numerical_jacobian);
    VectorXd jac_rhs = jac.rhs();
    for (int r=0; r<jac_rhs.rows(); r++) {
      rhs(row+r) = jac_rhs(r);
      // do not delete, will be pulled into SparseSystem below
      rows[row+r] = new SparseVector(jac.dimtotal());
    }
    for (Terms::const_iterator it=jac.terms().begin(); it!=jac.terms().end(); it++) {
      int offset = it->node()->_start;
      int nr = it->term().rows();
      for (int r=0; r<nr; r++) { // 0-entries not omitted
        rows[row+r]->set(offset, it->term().row(r));
      }
    }
    row += factor->dim();
  }
  return SparseSystem(num_rows, _dim_nodes, rows, rhs);
}

void Slam::print_stats() {
  double nnz = _R.nnz();
  double max_per_col = _R.max_nz();
  double dim = _dim_nodes;
  double per_col = nnz/dim;
  double fill_in = nnz/(dim*dim);
  cout << "iSAM statistics:" << endl;
  cout << "  Normalized chi-square value: " << normalized_chi2() << endl;
  cout << "  Weighted sum of squared errors: " << chi2() << endl;
  cout << "  Number of nodes: " << num_nodes() << endl;
  cout << "  Number of factors: " << num_factors() << endl;
  cout << "  Number of variables: " << _dim_nodes << endl;
  cout << "  Number of measurements: " << _dim_measure << endl;
  cout << "  Number of non-zero entries: " << nnz << endl;
  cout << "    max per column: " << max_per_col << endl;
  cout << "    avg per column: " << per_col << endl;
  cout << "    fill in: " << fill_in << "%" << endl;
}

}
