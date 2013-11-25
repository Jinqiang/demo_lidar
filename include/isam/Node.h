/**
 * @file Node.h
 * @brief Graph node for iSAM.
 * @author Michael Kaess
 * @version $Id: Node.h 8263 2013-04-10 14:02:19Z carlevar $
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

#include <list>
#include <Eigen/Dense>

#include "Element.h"

namespace Eigen {
  typedef Matrix<bool, Dynamic, 1> VectorXb;
}

namespace isam {

enum Selector {LINPOINT, ESTIMATE};

class Factor; // Factor.h not included here to avoid circular dependency

// Node of the graph also containing measurements (Factor).
class Node : public Element {
  friend std::ostream& operator<<(std::ostream& output, const Node& n) {
    n.write(output);
    return output;
  }

  static int _next_id;
  bool _deleted;

protected:

  std::list<Factor*> _factors; // list of adjacent factors

public:

  Node(const char* name, int dim) : Element(name, dim), _deleted(false) {
    _id = _next_id++;
  }

  virtual ~Node() {};

  virtual bool initialized() const = 0;

  virtual Eigen::VectorXd vector(Selector s = ESTIMATE) const = 0;

  virtual Eigen::VectorXd vector0() const = 0;
  virtual Eigen::VectorXb is_angle() const {return Eigen::VectorXb::Zero(_dim);}

  virtual void update(const Eigen::VectorXd& v) = 0;
  virtual void update0(const Eigen::VectorXd& v) = 0;

  virtual void linpoint_to_estimate() = 0;
  virtual void estimate_to_linpoint() = 0;
  virtual void swap_estimates() = 0;

  virtual void apply_exmap(const Eigen::VectorXd& v) = 0;
  virtual void self_exmap(const Eigen::VectorXd& v) = 0;

  void add_factor(Factor* e) {_factors.push_back(e);}
  void remove_factor(Factor* e) {_factors.remove(e);}

  const std::list<Factor*>& factors() {return _factors;}

  bool deleted() const { return _deleted; }

  void mark_deleted();
  void erase_marked_factors();

  virtual void write(std::ostream &out) const = 0;
};

// Generic template for easy instantiation of the multiple node types.
template <class T>
class NodeT : public Node {

 protected:
  T* _value;  // current estimate
  T* _value0; // linearization point

public:

  NodeT() : Node(T::name(), T::dim) {
    _value = NULL;
    _value0 = NULL;
  }

  NodeT(const char* name) : Node(name, T::dim) {
    _value = NULL;
    _value0 = NULL;
  }

  virtual ~NodeT() {
    delete _value;
    delete _value0;
  }

  void init(const T& t) {
    delete _value; delete _value0;
    _value = new T(t); _value0 = new T(t);
  }

  bool initialized() const {return _value != NULL;}

  T value(Selector s = ESTIMATE) const {return (s==ESTIMATE)?*_value:*_value0;}
  T value0() const {return *_value0;}

  Eigen::VectorXd vector(Selector s = ESTIMATE) const {return (s==ESTIMATE)?_value->vector():_value0->vector();}
  Eigen::VectorXd vector0() const {return _value0->vector();}
  Eigen::VectorXb is_angle() const {return _value->is_angle();}

  void update(const Eigen::VectorXd& v) {_value->set(v);}
  void update0(const Eigen::VectorXd& v) {_value0->set(v);}

  void linpoint_to_estimate() {*_value = *_value0;}
  void estimate_to_linpoint() {*_value0 = *_value;}
  void swap_estimates() {T tmp = *_value; *_value = *_value0; *_value0 = tmp;}

  void apply_exmap(const Eigen::VectorXd& v) {*_value = _value0->exmap(v);}
  void self_exmap(const Eigen::VectorXd& v) {*_value0 = _value0->exmap(v);}

  void write(std::ostream &out) const {
    out << name() << "_Node " << _id;
    if (_value != NULL) {
      out << " " << value();
    }
  }
};

}
