/**
 * @file slam_stereo.h
 * @brief Provides specialized nodes and factors for stereo vision applications.
 * @author Michael Kaess
 * @version $Id: slam_stereo.h 7956 2013-01-23 16:46:03Z kaess $
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
#include <sstream>
#include <math.h>
#include <Eigen/Dense>

#include "Node.h"
#include "Factor.h"
#include "Pose3d.h"
#include "Point3dh.h"

namespace isam {

class StereoMeasurement {
  friend std::ostream& operator<<(std::ostream& out, const StereoMeasurement& t) {
    t.write(out);
    return out;
  }

public:
  double u;
  double v;
  double u2;
  bool valid;

  StereoMeasurement(double u, double v, double u2)
    : u(u), v(v), u2(u2), valid(true) {
  }
  StereoMeasurement(double u, double v, double u2, bool valid)
    : u(u), v(v), u2(u2), valid(valid) {
  }

  Eigen::Vector2d left_pixel() const {Eigen::Vector2d V(2); V << u, v; return V;}

  Eigen::Vector2d right_pixel() const {Eigen::Vector2d V(2); V << u2, v; return V;}

  double disparity() const {return u-u2;}

  Eigen::Vector3d vector() const {
    Eigen::Vector3d tmp(u, v, u2);
    return tmp;
  }

  void write(std::ostream &out) const {
    out << "(" << u << ", " << v << ", " << u2 << ")";
  }
};

class StereoCamera { // for now, camera and robot are identical
  double _f;
  Eigen::Vector2d _pp;
  double _b;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  StereoCamera() : _f(1), _pp(Eigen::Vector2d(0.5,0.5)), _b(0.1) {}
  StereoCamera(double f, const Eigen::Vector2d& pp, double b) : _f(f), _pp(pp), _b(b) {}

  inline double focalLength() const {return _f;}

  inline Eigen::Vector2d principalPoint() const {return _pp;}

  inline double baseline() const {return _b;}

  StereoMeasurement project(const Pose3d& pose, const Point3dh& Xw) const {
    Point3dh X = pose.transform_to(Xw);
    // camera system has z pointing forward, instead of x
    double x = X.y();
    double y = X.z();
    double z = X.x();
    double w = X.w();
    // left camera
    double fz = _f / z;
    double u = x * fz + _pp(0);
    double v = y * fz + _pp(1);
    // right camera
    double u2 = (x - w*_b) * fz + _pp(0);
    bool valid = ((w==0&&z>0) || (z/w) > 0.); // infront of camera?

    // Setting a point valid can make the system (hj 25/10/2012)
    // ill conditioned, i.e. if point goes behind all cameras that see it
    //valid = true;

    return StereoMeasurement(u, v, u2, valid);
  }


  Point3dh backproject(const Pose3d& pose, const StereoMeasurement& measure) const {
    double lx = (measure.u-_pp(0))*_b;
    double ly = (measure.v-_pp(1))*_b;
    double lz = _f*_b;
    double lw = measure.u - measure.u2;
    if (lw<0.) {
      std::cout << "Warning: StereoCamera.backproject called with negative disparity\n";
    }
    Point3dh X(lz, lx, ly, lw);

    return pose.transform_from(X);
  }

};

/**
 * Stereo observation of a 3D homogeneous point;
 * projective or Euclidean geometry depending on constructor used.
 */
class Stereo_Factor : public FactorT<StereoMeasurement> {
  Pose3d_Node* _pose;
  Point3d_Node* _point;
  Point3dh_Node* _point_h;
  StereoCamera* _camera;
  bool _relative;
  Pose3d_Node* _base;

public:

  // constructor for projective geometry
  Stereo_Factor(Pose3d_Node* pose, Point3dh_Node* point, StereoCamera* camera,
                         const StereoMeasurement& measure, const Noise& noise, bool relative = false)
    : FactorT<StereoMeasurement>("Stereo_Factor", 3, noise, measure),
      _pose(pose), _point(NULL), _point_h(point), _camera(camera), _relative(relative), _base(NULL) {
    // StereoCamera could also be a node later (either with 0 variables,
    // or with calibration as variables)
    _nodes.resize(2);
    _nodes[0] = pose;
    _nodes[1] = point;
    // for relative parameterization recover base pose
    if (_relative) {
      if (!point->factors().empty()) {
        _nodes.resize(3);
        _nodes[2] = point->factors().front()->nodes()[0];
        // todo: first factor might refer to a prior or other type of node...
        _base = dynamic_cast<Pose3d_Node*>(_nodes[2]);
      } else {
        _base = _pose;
      }
      point->set_base(_base);
    }
  }

  // constructor for Euclidean geometry
  // WARNING: only use for points at short range
  Stereo_Factor(Pose3d_Node* pose, Point3d_Node* point, StereoCamera* camera,
                         const StereoMeasurement& measure, const Noise& noise, bool relative = false)
    : FactorT<StereoMeasurement>("Stereo_Factor", 3, noise, measure),
      _pose(pose), _point(point), _point_h(NULL), _camera(camera), _relative(relative), _base(NULL) {
    _nodes.resize(2);
    _nodes[0] = pose;
    _nodes[1] = point;
    // for relative parameterization recover base pose
    if (_relative) {
      if (!point->factors().empty()) {
        _nodes.resize(3);
        _nodes[2] = point->factors().front()->nodes()[0];
        // todo: first factor might refer to a prior or other type of node...
        _base = dynamic_cast<Pose3d_Node*>(_nodes[2]);
      } else {
        _base = _pose;
      }
      point->set_base(_base);
    }
  }

  void initialize() {
    require(_pose->initialized(), "Stereo_Factor requires pose to be initialized");
    bool initialized = (_point_h!=NULL) ? _point_h->initialized() : _point->initialized();
    if (!initialized) {
      Point3dh predict;
      if (_relative) {
        predict = _camera->backproject(Pose3d(), _measure);
      } else {
        predict = _camera->backproject(_pose->value(), _measure);
      }
      // normalize homogeneous vector
      predict = Point3dh(predict.vector()).normalize();
      if (_point_h!=NULL) {
        _point_h->init(predict);
      } else {
        _point->init(predict.to_point3d());
      }
    }
  }

  Eigen::VectorXd basic_error(Selector s = ESTIMATE) const {
    Point3dh point = (_point_h!=NULL) ? _point_h->value(s) : _point->value(s);
    Pose3d pose = _pose->value(s);
    if (_base) {
      // pose of camera relative to base camera (which might be this one!)
      pose = pose.ominus(_base->value(s));
    }
    StereoMeasurement predicted = _camera->project(pose, point);
    if (_point_h!=NULL || predicted.valid == true) {
      return (predicted.vector() - _measure.vector());
    } else {
      //std::cout << "Warning - StereoFactor.basic_error: point behind camera, dropping term.\n";
      std::cout << "Warning - StereoFactor.basic_error: " << this
          << " #: " << _nodes[1]->factors().size() << " " << _nodes[0]->factors().size()
          << " \npoint: " << point << " ptr: " << _point
          << " pose: " << pose << " ptr: " << _pose
          << " \npredicted: " << predicted.vector().transpose()
          << " measure: " << _measure.vector().transpose() << std::endl;
      return Eigen::Vector3d::Zero();
    }
  }

};

}
