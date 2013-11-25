/**
 * @file Factor.h
 * @brief Graph factor for iSAM.
 * @author Michael Kaess
 * @version $Id: Factor.h 7610 2012-10-25 10:21:12Z hordurj $
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
#include <string>

#include <math.h> // for sqrt
#include <Eigen/Dense>

#include "util.h"
#include "Jacobian.h"
#include "Element.h"
#include "Node.h"
#include "Noise.h"
#include "numericalDiff.h"

namespace isam {

//typedef double (*cost_func_t)(double);


// Factor of the graph of measurements between Nodes.
class Factor : public Element, Function {
  friend std::ostream& operator<<(std::ostream& output, const Factor& e) {
    e.write(output);
    return output;
  }

  cost_func_t *ptr_cost_func;

  static int _next_id;
  bool _deleted;
protected:

  const Noise _noise;

  std::vector<Node*> _nodes; // list of nodes affected by measurement

public:

  virtual Eigen::VectorXd error(Selector s = ESTIMATE) const {
    Eigen::VectorXd err = _noise.sqrtinf() * basic_error(s);
    // optional modified cost function
    if (*ptr_cost_func) {
      for (int i=0; i<err.size(); i++) {
        double val = err(i);
        err(i) = ((val>=0)?1.:(-1.)) * sqrt((*ptr_cost_func)(val));
      }
    }
    return err;
  }

  std::vector<Node*>& nodes() {return _nodes;}

  Factor(const char* name, int dim, const Noise& noise)
    : Element(name, dim), ptr_cost_func(NULL), _deleted(false), _noise(noise) {
#ifndef NDEBUG
    // all lower triagular entries below the diagonal must be 0
    for (int r=0; r<_noise.sqrtinf().rows(); r++) {
      for (int c=0; c<r; c++) {
        requireDebug(_noise.sqrtinf()(r,c)==0, "Factor::Factor: sqrtinf must be upper triangular!");
      }
    }
#endif
    _id = _next_id++;
  }

  virtual ~Factor() {}

  virtual void initialize() = 0;

  virtual void initialize_internal() {
    for (unsigned int i=0; i<_nodes.size(); i++) {
      _nodes[i]->add_factor(this);
    }
    initialize();
  }

  virtual void set_cost_function(cost_func_t* ptr) {ptr_cost_func = ptr;}

  virtual Eigen::VectorXd basic_error(Selector s = ESTIMATE) const = 0;

  virtual const Eigen::MatrixXd& sqrtinf() const {return _noise.sqrtinf();}

  Eigen::VectorXd evaluate() const {
    return error(LINPOINT);
  }

  virtual Jacobian jacobian_internal(bool force_numerical) {
    if (force_numerical) {
      // ignore any symbolic derivative provided by user
      return Factor::jacobian();
    } else {
      return jacobian();
    }
  }

  // can be replaced by symbolic derivative by user
  virtual Jacobian jacobian() {
    Eigen::MatrixXd H = numerical_jacobian();
    Eigen::VectorXd r = error(LINPOINT);
    Jacobian jac(r);
    int position = 0;
    int n_measure = dim();
    for (unsigned int i=0; i<_nodes.size(); i++) {
      int n_var = _nodes[i]->dim();
      Eigen::MatrixXd Hi = H.block(0, position, n_measure, n_var);
      position += n_var;
      jac.add_term(_nodes[i], Hi);
    }
    return jac;
  }

  int num_measurements() const {
    return dim();
  }

  void mark_deleted() { _deleted = true; }
  bool deleted() const { return _deleted; }

  virtual void write(std::ostream &out) const {
    Element::write(out);
    for (unsigned int i=0; i<_nodes.size(); i++) {
      if (_nodes[i]) {
        out << " " << _nodes[i]->unique_id();
      }
    }
  }

private:

  virtual Eigen::MatrixXd numerical_jacobian() {
    return numericalDiff(*this);
  }

};

/**
 * Convert upper triangular square root information matrix to string.
 * @param sqrtinf Upper triangular square matrix.
 */
inline std::string noise_to_string(const Noise& noise) {
  int nrows = noise.sqrtinf().rows();
  int ncols = noise.sqrtinf().cols();
  require(nrows==ncols, "slam2d::sqrtinf_to_string: matrix must be square");
  std::stringstream s;
  s << "{";
  bool first = true;
  for (int r=0; r<nrows; r++) {
    for (int c=r; c<ncols; c++) {
      if (first) {
        first = false;
      } else {
        s << ",";
      }
      s << noise.sqrtinf()(r,c);
    }
  }
  s << "}";
  return s.str();
}

// Generic template for easy instantiation of new factors
template <class T>
class FactorT : public Factor {

protected:

  const T _measure;
  
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  FactorT(const char* name, int dim, const Noise& noise, const T& measure) : Factor(name, dim, noise), _measure(measure) {}

  const T& measurement() const {return _measure;}

  void write(std::ostream &out) const {
    Factor::write(out);
    out << " " << _measure << " " << noise_to_string(_noise);
  }

};


}
