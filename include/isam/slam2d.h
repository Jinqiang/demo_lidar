/**
 * @file slam2d.h
 * @brief Provides specialized nodes and factors for 2D SLAM.
 * @author Michael Kaess
 * @author Hordur Johannsson
 * @version $Id: slam2d.h 6334 2012-03-22 18:53:24Z hordurj $
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
#include <Eigen/Dense>

#include "Node.h"
#include "Factor.h"
#include "Pose2d.h"
#include "Point2d.h"
#include "Anchor.h"

namespace Eigen {
typedef Matrix< double , 1 , 1> Vector1d;
}

namespace isam {

typedef NodeT<Pose2d> Pose2d_Node;
typedef NodeT<Point2d> Point2d_Node;


/**
 * Prior on Point2d.
 */
class Point2d_Factor : public FactorT<Point2d> {
  Point2d_Node* _point;

public:

  /**
   * Constructor.
   * @param point The point node the prior acts on.
   * @param prior The actual prior measurement.
   * @param noise The 2x2 square root information matrix (upper triangular).
   */
  Point2d_Factor(Point2d_Node* point, const Point2d& prior, const Noise& noise)
    : FactorT<Point2d>("Point2d_Factor", 2, noise, prior), _point(point) {
    _nodes.resize(1);
    _nodes[0] = point;
  }

  void initialize() {
    if (!_point->initialized()) {
      Point2d predict = _measure;
      _point->init(predict);
    }
  }

  Eigen::VectorXd basic_error(Selector s = LINPOINT) const {
    return (_point->vector(s) - _measure.vector());
  }

};

/**
 * Prior on Pose2d.
 */
class Pose2d_Factor : public FactorT<Pose2d> {
  Pose2d_Node* _pose;

public:

  /**
   * Constructor.
   * @param pose The pose node the prior acts on.
   * @param prior The actual prior measurement.
   * @param noise The 3x3 square root information matrix (upper triangular).
   */
  Pose2d_Factor(Pose2d_Node* pose, const Pose2d& prior, const Noise& noise)
    : FactorT<Pose2d>("Pose2d_Factor", 3, noise, prior), _pose(pose) {
    _nodes.resize(1);
    _nodes[0] = pose;
  }

  void initialize() {
    if (!_pose->initialized()) {
      Pose2d predict = _measure;
      _pose->init(predict);
    }
  }

  Eigen::VectorXd basic_error(Selector s = LINPOINT) const {
    Eigen::VectorXd err = _pose->vector(s) - _measure.vector();
    err(2) = standardRad(err(2));
    return err;
  }

  Jacobian jacobian() {
    Eigen::MatrixXd M = sqrtinf(); // derivatives are all 1 (eye)
    Eigen::VectorXd r = sqrtinf() * basic_error();
    Jacobian jac(r);
    jac.add_term(_nodes[0], M);
    return jac;
  }

};

/**
 * Odometry or loop closing constraint, from pose1 to pose2.
 */
class Pose2d_Pose2d_Factor : public FactorT<Pose2d> {
  Pose2d_Node* _pose1;
  Pose2d_Node* _pose2;
  
public:

  /**
   * Constructor.
   * @param pose1 The pose from which the measurement starts.
   * @param pose2 The pose to which the measurement extends.
   * @param measure The relative measurement from pose1 to pose2 (pose2 in pose1's frame).
   * @param noise The 3x3 square root information matrix (upper triangular).
   * @param anchor1 Optional anchor node for trajectory to which pose1 belongs to.
   * @param anchor2 Optional anchor node for trajectory to which pose2 belongs to.
   */
  Pose2d_Pose2d_Factor(Pose2d_Node* pose1, Pose2d_Node* pose2,
      const Pose2d& measure, const Noise& noise,
      Anchor2d_Node* anchor1 = NULL, Anchor2d_Node* anchor2 = NULL)
    : FactorT<Pose2d>("Pose2d_Pose2d_Factor", 3, noise, measure),
    _pose1(pose1), _pose2(pose2) {
    require((anchor1==NULL && anchor2==NULL) || (anchor1!=NULL && anchor2!=NULL),
        "slam2d: Pose2d_Pose2d_Factor requires either 0 or 2 anchor nodes");

    require((anchor1==NULL && anchor2==NULL) || (anchor1!=anchor2),
        "slam2d: Pose2d_Pose2d_Factor requires anchors to be different");

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
    Pose2d_Node* pose1 = _pose1;
    Pose2d_Node* pose2 = _pose2;
    require(pose1->initialized() || pose2->initialized(),
        "slam2d: Pose2d_Pose2d_Factor requires pose1 or pose2 to be initialized");

    if (!pose1->initialized() && pose2->initialized()) {
      // reverse constraint 
      Pose2d p2 = pose2->value();
      Pose2d z;
      Pose2d predict = p2.oplus(z.ominus(_measure));
      pose1->init(predict);
    } else if (pose1->initialized() && !pose2->initialized()) {
      Pose2d p1 = pose1->value();
      Pose2d predict = p1.oplus(_measure);
      pose2->init(predict);
    }

    if (_nodes.size()==4) {
      Anchor2d_Node* anchor1 = dynamic_cast<Anchor2d_Node*>(_nodes[2]);
      Anchor2d_Node* anchor2 = dynamic_cast<Anchor2d_Node*>(_nodes[3]);

      if (!anchor1->initialized()) {
        anchor1->set_prior();
        anchor1->init(Pose2d());
      }

      if (!anchor2->initialized()) {
        Pose2d p1 = pose1->value();
        Pose2d p2 = pose2->value();
        Pose2d a1 = anchor1->value(); 
        // note 1: ominus is tricky: a.ominus(b)=B^(-1)*A
        // note 2: unlike vectors, the order of composition matters
        //         for poses! (think in terms of matrix muliplication)
        Pose2d zero;
        Pose2d d = a1.oplus(p1).oplus(_measure).oplus(zero.ominus(p2));
        anchor2->init(d);
      }
      if (  (anchor1->parent() != NULL && anchor2->parent() != NULL && anchor1->parent() != anchor2->parent())
              || (anchor1->parent() == NULL && anchor2->parent() == NULL) 
              || (anchor1->parent() == NULL && anchor2->parent() != anchor1) 
              || (anchor1->parent() != anchor2 && anchor2->parent() == NULL)) {
        Pose2d p1 = pose1->value();
        Pose2d p2 = pose2->value();
        Pose2d a1 = anchor1->value();
        Pose2d a2 = anchor2->value();
        // Compute the transformation from anchor1 and anchor2 frames
        //Pose2d d = (a2.oplus(p2)).ominus(a1.oplus(p1).oplus(_measure));

        Pose2d zero;
        Pose2d a2_p2 = a2.oplus(p2);
        Pose2d d = a1.oplus(p1).oplus(_measure).oplus(zero.ominus(a2_p2));

        anchor1->merge(anchor2, d);
      }
    }
  }

  Eigen::VectorXd basic_error(Selector s = LINPOINT) const {
    Pose2d p1(_nodes[0]->vector(s));
    Pose2d p2(_nodes[1]->vector(s));
    Pose2d predicted;
    if (_nodes.size()==4) {
      Pose2d a1(_nodes[2]->vector(s));
      Pose2d a2(_nodes[3]->vector(s));
      // see notes above
      predicted = (a2.oplus(p2)).ominus(a1.oplus(p1));
    } else {
      predicted = p2.ominus(p1);
    }
    Eigen::VectorXd err = predicted.vector() - _measure.vector();
    err(2) = standardRad(err(2));
    return err;
  }
      
  Jacobian jacobian() {
    if (_nodes.size()==4) { // symbolic available only without anchor nodes
      return Factor::jacobian();
    } else {
      const Pose2d& p1 = _pose1->value0();
      const Pose2d& p2 = _pose2->value0();
      Pose2d p = p2.ominus(p1);
      double c = cos(p1.t());
      double s = sin(p1.t());

      double dx = p.x() - _measure.x();
      double dy = p.y() - _measure.y();
      double dt = standardRad(p.t() - _measure.t());

      Eigen::MatrixXd M1(3,3);
      M1 <<
        -c, -s,  p.y(),
        s,  -c,  -p.x(),
        0.,  0., -1.;
      M1 = sqrtinf() * M1;

      Eigen::MatrixXd M2(3,3);
      M2 <<
        c,   s,   0.,
        -s,  c,   0.,
        0.,  0.,  1.;
      M2 = sqrtinf() * M2;

      Pose2d pp(dx, dy, dt);

      Eigen::VectorXd r = sqrtinf() * pp.vector();

      Jacobian jac(r);
      jac.add_term(_pose1, M1);
      jac.add_term(_pose2, M2);

      return jac;
    }
  }
};

/**
 * Landmark observation.
 */
class Pose2d_Point2d_Factor : public FactorT<Point2d> {
  Pose2d_Node* _pose;
  Point2d_Node* _point;

public:

  /**
   * Constructor.
   * @param pose The pose from which the landmark is observed.
   * @param point The point or landmark that is observed
   * @param measure The relative observation of the landmark in the pose's frame.
   * @param noise The 2x2 square root information matrix (upper triangular).
   */
  Pose2d_Point2d_Factor(Pose2d_Node* pose, Point2d_Node* point,
      const Point2d& measure, const Noise& noise)
    : FactorT<Point2d>("Pose2d_Point2d_Factor", 2, noise, measure), _pose(pose), _point(point) {
    _nodes.resize(2);
    _nodes[0] = pose;
    _nodes[1] = point;
  }

  void initialize() {
    require(_pose->initialized(),
        "slam2d: Pose2d_Point2d_Factor requires pose to be initialized");
    if (!_point->initialized()) {
      Pose2d p = _pose->value();
      Point2d predict = p.transform_from(_measure);
      _point->init(predict);
    }
  }

  Eigen::VectorXd basic_error(Selector s = LINPOINT) const {
    Pose2d po(_nodes[0]->vector(s));
    Point2d pt(_nodes[1]->vector(s));
    Point2d p = po.transform_to(pt);
    Eigen::VectorXd predicted = p.vector();
    return (predicted - _measure.vector());
  }

  Jacobian jacobian() {
    Pose2d po = _pose->value0();
    Point2d pt = _point->value0();
    double c = cos(po.t());
    double s = sin(po.t());
    double dx = pt.x() - po.x();
    double dy = pt.y() - po.y();
    // f(x)
    double x =  c*dx + s*dy; // relative forward position of landmark point from pose
    double y = -s*dx + c*dy; // relative position to the left
    Eigen::MatrixXd M1(2,3);
    M1 <<
      -c, -s,  y,
      s,  -c, -x;
    M1 = sqrtinf() * M1;
    Eigen::MatrixXd M2(2,2);
    M2 << 
      c,   s,
      -s,  c;
    M2 = sqrtinf() * M2;
    Point2d p(x, y);
    Eigen::VectorXd r = sqrtinf() * (p.vector() - _measure.vector());
    Jacobian jac(r);
    jac.add_term(_pose, M1);
    jac.add_term(_point, M2);
    return jac;
  }

};


}
