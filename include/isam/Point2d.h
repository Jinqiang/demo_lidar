/**
 * @file Point2d.h
 * @brief Simple 2D point class.
 * @author Michael Kaess
 * @version $Id: Point2d.h 8263 2013-04-10 14:02:19Z carlevar $
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

#include <Eigen/Dense>

namespace Eigen {
 typedef Matrix<bool, Dynamic, 1> VectorXb;
}

namespace isam {

class Point2d {
  friend std::ostream& operator<<(std::ostream& out, const Point2d& p) {
    p.write(out);
    return out;
  }

  double _x;
  double _y;

public:
  static const int dim = 2;
  static const int size = 2;
  static const char* name() {
    return "Point2d";
  }
  Point2d() : _x(0), _y(0) {}
  Point2d(double x, double y) : _x(x), _y(y) {}
  Point2d(const Eigen::Vector2d& vec) : _x(vec(0)), _y(vec(1)) {}

  double x() const {return _x;}
  double y() const {return _y;}

  void set_x(double x) {_x = x;}
  void set_y(double y) {_y = y;}

  Point2d exmap(const Eigen::Vector2d& delta) const {
    Point2d res = *this;
    res._x += delta(0);
    res._y += delta(1);
    return res;
  }

  Eigen::Vector2d vector() const {
    Eigen::Vector2d v(_x, _y);
    return v;
  }
  void set(double x, double y) {
    _x = x;
    _y = y;
  }
  void set(const Eigen::Vector2d& v) {
    _x = v(0);
    _y = v(1);
  }
  void write(std::ostream &out) const {
    out << "(" << _x << ", " << _y << ")";
  }
  
  Eigen::VectorXb is_angle() const {
    return Eigen::VectorXb::Zero(size);
  }
  
};

}
