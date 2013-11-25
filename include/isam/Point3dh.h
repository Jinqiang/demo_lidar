/**
 * @file Point3dh.h
 * @brief 3D homogeneous point class.
 * @author Michael Kaess
 * @version $Id: Point3dh.h 8263 2013-04-10 14:02:19Z carlevar $
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
#include <Eigen/Geometry>

#include "Point3d.h"

namespace isam {

class Point3dh {
  friend std::ostream& operator<<(std::ostream& out, const Point3dh& p) {
    p.write(out);
    return out;
  }

  double _x;
  double _y;
  double _z;
  double _w;

public:

  static const int dim = 3;
  static const int size = 4;
  static const char* name() {
    return "Point3dh";
  }

  Point3dh() : _x(0.), _y(0.), _z(0.), _w(1.) {} // note: 0,0,0,0 is not allowed!
  Point3dh(double x, double y, double z, double w) : _x(x), _y(y), _z(z), _w(w) {}
  Point3dh(const Eigen::Vector4d& vec) : _x(vec(0)), _y(vec(1)), _z(vec(2)), _w(vec(3)) {}
  Point3dh(const Point3d& p) : _x(p.x()), _y(p.y()), _z(p.z()), _w(1.) {}

  double x() const {return _x;}
  double y() const {return _y;}
  double z() const {return _z;}
  double w() const {return _w;}

  void x(double x) {_x = x;}
  void y(double y) {_y = y;}
  void z(double z) {_z = z;}
  void w(double w) {_w = w;}

  static Eigen::Vector4d delta3_to_homogeneous(const Eigen::Vector3d& delta) {
    double theta = delta.norm();
    double S;
    if (theta < 0.0001) { // sqrt4(machine precession)
      S = 0.5 + theta*theta/48.;
    } else {
      S = sin(0.5*theta)/theta;
    }
    double C = cos(0.5*theta);
    return Eigen::Vector4d(S*delta(0), S*delta(1), S*delta(2), C);
  }

  Point3dh exmap(const Eigen::Vector3d& delta) const {
#if 0
    // solution with Householder matrix following HZ second edition

    //    std::cout << norm() << std::endl;
    assert(norm() == 1);
    // make sure not above PI, otherwise needs to be corrected
    assert(delta.norm() < 3.14);

    Eigen::Vector4d delta4 = delta3_to_homogeneous(delta);
    // multiply with Householder matrix, without explicitly calculating the matrix
    Eigen::Vector4d e(0,0,0,1);
    Eigen::Vector4d v = vector() + (x()<0.?-1:1) * norm() * e;

    // checking if sign is correct to yield (0,0,0,1), otherwise fix
    Eigen::Vector4d test = vector() - 2*v*(v.transpose()*vector())/(v.transpose()*v);
    if (test(3)<0) {
      v = v - 2 * (x()<0.?-1:1) * norm() * e;
    }

    Eigen::Vector4d res = delta4 - 2*v*(v.transpose()*delta4)/(v.transpose()*v);
    //    std::cout << vector() << " " << res << std::endl;
    return Point3dh(res);
#else
#if 1
    // solution by mapping to unit Quaternion and back (both are unit spheres!)
    Eigen::Quaterniond delta_quat = Rot3d::delta3_to_quat(delta);
    Eigen::Quaterniond quat(_w, _x, _y, _z);
    //    quat.normalize(); // just to be safe...
    Eigen::Quaterniond res = quat * delta_quat;
    return Point3dh(res.x(), res.y(), res.z(), res.w());
#else
#if 1
    // 3DOF exmap - not recommended
    Point3dh res = *this;
    res.normalize();
    Eigen::Vector4d::Index pos;
    res.vector().cwiseAbs().maxCoeff(&pos);
    //    std::cout << *this << " -- " << pos << " -- " << res.vector() << std::endl;
    // only update 3 out of 4 entries
    int idx = 0;
    if (pos!=0) {res._x += delta(idx); idx++;}
    if (pos!=1) {res._y += delta(idx); idx++;}
    if (pos!=2) {res._z += delta(idx); idx++;}
    if (pos!=3) {res._w += delta(idx);}
    res.normalize();
    return res;
#else
    // 4DOF exmap - not feasible without constrained optimization
    Point3dh res = *this;
    res._x += delta(0);
    res._y += delta(1);
    res._z += delta(2);
    res._w += delta(3);
    return res;
#endif
#endif
#endif
  }

  Eigen::Vector4d vector() const {
    Eigen::VectorXd tmp(4);
    tmp << _x, _y, _z, _w;
    return tmp;
  }

  Eigen::VectorXb is_angle() const {
    Eigen::VectorXb isang (4);
    isang << false, false, false, false;
    return isang;
  }

  void set(double x, double y, double z, double w) {
    _x = x;
    _y = y;
    _z = z;
    _w = w;
  }

  void set(const Eigen::Vector4d& v) {
    _x = v(0);
    _y = v(1);
    _z = v(2);
    _w = v(3);
  }

  Point3d to_point3d() const {
    double w_inv = 1. / _w;
    return Point3d(_x*w_inv, _y*w_inv, _z*w_inv);
  }

  double norm() const {
    return sqrt(_x*_x + _y*_y + _z*_z + _w*_w);
  }

  Point3dh normalize() {
    double a = 1. / norm();
    _x *= a;
    _y *= a;
    _z *= a;
    _w *= a;
    return *this;
  }

  void write(std::ostream &out) const {
    out << "(" << _x << ", " << _y << ", " << _z << ", " << _w << ")";
  }
};

}
