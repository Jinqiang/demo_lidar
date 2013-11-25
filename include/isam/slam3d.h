/**
 * @file slam3d.h
 * @brief Provides specialized nodes and factors for 3D SLAM
 * @author Michael Kaess
 * @version $Id: slam3d.h 7610 2012-10-25 10:21:12Z hordurj $
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
#include "Node.h"
#include "Factor.h"
#include "Pose3d.h"
#include "Point3d.h"
#include "Anchor.h"

namespace isam {

typedef NodeT<Pose3d> Pose3d_Node;

template<class T>
class Point3dT_Node : public NodeT<T>
{
  Pose3d_Node* _base;
public:
  Point3dT_Node() : _base(NULL) {}
  Point3dT_Node(const char* name) : NodeT<T>(name), _base(NULL) {}

  ~Point3dT_Node() {}

  void set_base(Pose3d_Node* base) { _base = base; }
  Pose3d_Node* base() { return _base; }
};

typedef Point3dT_Node<Point3d> Point3d_Node;
typedef Point3dT_Node<Point3dh> Point3dh_Node;

class Pose3d_Factor : public FactorT<Pose3d> {
  Pose3d_Node* _pose;

public:

  /**
   * Constructor.
   * @param pose The pose node the prior acts on.
   * @param prior The actual prior measurement.
   * @param noise The 6x6 square root information matrix (upper triangular).
   */
  Pose3d_Factor(Pose3d_Node* pose, const Pose3d& prior, const Noise& noise)
    : FactorT<Pose3d>("Pose3d_Factor", 6, noise, prior), _pose(pose) {
    _nodes.resize(1);
    _nodes[0] = pose;
  }

  void initialize() {
    if (!_pose->initialized()) {
      Pose3d predict = _measure;
      _pose->init(predict);
    }
  }

  Eigen::VectorXd basic_error(Selector s = ESTIMATE) const {
    Eigen::VectorXd err = _nodes[0]->vector(s) - _measure.vector();
    err(3) = standardRad(err(3));
    err(4) = standardRad(err(4));
    err(5) = standardRad(err(5));
    return err;
  }
};

class Pose3d_Pose3d_Factor : public FactorT<Pose3d> {
  Pose3d_Node* _pose1;
  Pose3d_Node* _pose2;

public:

  /**
   * Constructor.
   * @param pose1 The pose from which the measurement starts.
   * @param pose2 The pose to which the measurement extends.
   * @param measure The relative measurement from pose1 to pose2 (pose2 in pose1's frame).
   * @param noise The 6x6 square root information matrix (upper triangular).
   * @param anchor1 Optional anchor node for trajectory to which pose1 belongs to.
   * @param anchor2 Optional anchor node for trajectory to which pose2 belongs to.
   */
  Pose3d_Pose3d_Factor(Pose3d_Node* pose1, Pose3d_Node* pose2,
      const Pose3d& measure, const Noise& noise,
      Anchor3d_Node* anchor1 = NULL, Anchor3d_Node* anchor2 = NULL)
    : FactorT<Pose3d>("Pose3d_Pose3d_Factor", 6, noise, measure), _pose1(pose1), _pose2(pose2) {
    require((anchor1==NULL && anchor2==NULL) || (anchor1!=NULL && anchor2!=NULL),
        "slam3d: Pose3d_Pose3d_Factor requires either 0 or 2 anchor nodes");
    if (anchor1) { // offset between two relative pose graphs
      _nodes.resize(4);
      _nodes[2] = anchor1;
      _nodes[3] = anchor2;
    } else {
      _nodes.resize(2);
    }
    _nodes[0] = pose1;
    _nodes[1] = pose2;
  }

  void initialize() {
    require(_pose1->initialized() || _pose2->initialized(),
        "slam3d: Pose3d_Pose3d_Factor requires pose1 or pose2 to be initialized");

    if (!_pose1->initialized() && _pose2->initialized()) {
      // Reverse constraint 
      Pose3d p2 = _pose2->value();
      Pose3d z;
      Pose3d predict = p2.oplus(z.ominus(_measure));
      _pose1->init(predict);
    } else if (_pose1->initialized() && !_pose2->initialized()) {
      Pose3d p1 = _pose1->value();
      Pose3d predict = p1.oplus(_measure);
      _pose2->init(predict);
    }
    if (_nodes.size()==4) {
      Anchor3d_Node* anchor1 = dynamic_cast<Anchor3d_Node*>(_nodes[2]);
      Anchor3d_Node* anchor2 = dynamic_cast<Anchor3d_Node*>(_nodes[3]);
      if (!anchor1->initialized()) {
        anchor1->set_prior();
        anchor1->init(Pose3d());
      }
      if (!anchor2->initialized()) {
        Pose3d p1 = _pose1->value();
        Pose3d p2 = _pose2->value();
        Pose3d a1 = anchor1->value();
        // see notes in slam2d.h
        Pose3d zero;
        Pose3d d = a1.oplus(p1).oplus(_measure).oplus(zero.ominus(p2));
        anchor2->init(d);
      }
      if (  (anchor1->parent() != NULL && anchor2->parent() != NULL && anchor1->parent() != anchor2->parent())
              || (anchor1->parent() == NULL && anchor2->parent() == NULL)
              || (anchor1->parent() == NULL && anchor2->parent() != anchor1)
              || (anchor1->parent() != anchor2 && anchor2->parent() == NULL)) {
        Pose3d p1 = _pose1->value();
        Pose3d p2 = _pose2->value();
        Pose3d a1 = anchor1->value();
        Pose3d a2 = anchor2->value();
        // Compute the transformation from anchor1 and anchor2 frames
        //Pose2d d = (a2.oplus(p2)).ominus(a1.oplus(p1).oplus(_measure));

        Pose3d zero;
        Pose3d a2_p2 = a2.oplus(p2);
        Pose3d d = a1.oplus(p1).oplus(_measure).oplus(zero.ominus(a2_p2));

        anchor1->merge(anchor2, d);
      }
    }
  }

  Eigen::VectorXd basic_error(Selector s = ESTIMATE) const {
    const Pose3d& p1 = _pose1->value(s);
    const Pose3d& p2 = _pose2->value(s);
    Pose3d predicted;
    if (_nodes.size()==4) {
      const Pose3d& a1 = dynamic_cast<Pose3d_Node*>(_nodes[2])->value(s);
      const Pose3d& a2 = dynamic_cast<Pose3d_Node*>(_nodes[3])->value(s);
      // see notes in slam2d.h
      predicted = (a2.oplus(p2)).ominus(a1.oplus(p1));
    } else {
      predicted = p2.ominus(p1);
    }
    Eigen::VectorXd err = predicted.vector() - _measure.vector();
    err(3) = standardRad(err(3));
    err(4) = standardRad(err(4));
    err(5) = standardRad(err(5));
    return err;
  }

};

class Pose3d_Point3d_Factor : public FactorT<Point3d> {
  Pose3d_Node* _pose;
  Point3d_Node* _point;

public:

  /**
   * Constructor.
   * @param pose The pose from which the landmark is observed.
   * @param point The point or landmark that is observed
   * @param measure The relative observation of the landmark in the pose's frame.
   * @param noise The 3x3 square root information matrix (upper triangular).
   */
  Pose3d_Point3d_Factor(Pose3d_Node* pose, Point3d_Node* point,
      const Point3d& measure, const Noise& noise)
    : FactorT<Point3d>("Pose3d_Point3d_Factor", 3, noise, measure), _pose(pose), _point(point) {
    _nodes.resize(2);
    _nodes[0] = pose;
    _nodes[1] = point;
  }

  void initialize() {
    require(_pose->initialized(), "slam3d: Pose3d_Point3d_Factor requires pose to be initialized");
    if (!_point->initialized()) {
      Pose3d p = _pose->value();
      Point3d predict = p.transform_from(_measure);
      _point->init(predict);
    }
  }

  Eigen::VectorXd basic_error(Selector s = ESTIMATE) const {
    const Pose3d& po = _pose->value(s);
    const Point3d& pt = _point->value(s);
    Point3d p = po.transform_to(pt);
    Eigen::VectorXd predicted = p.vector();
    return (predicted - _measure.vector());
  }
};

}
