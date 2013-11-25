/**
 * @file Rot3d.h
 * @brief 3D rotation class.
 * @author Michael Kaess
 * @version $Id: Rot3d.h 7625 2012-10-30 18:40:38Z kaess $
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

/** @class isam::Rot3d
 *
 * For conventions, see isam::Pose3d.h
 */

#define USE_QUATERNIONS

#pragma once

#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/LU>

#include <isam/util.h>

namespace isam {

class Rot3d {
  friend std::ostream& operator<<(std::ostream& out, const Rot3d& p) {
    p.write(out);
    return out;
  }

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  static Eigen::Matrix3d euler_to_wRo(double yaw, double pitch, double roll) {
    double c__ = cos(yaw);
    double _c_ = cos(pitch);
    double __c = cos(roll);
    double s__ = sin(yaw);
    double _s_ = sin(pitch);
    double __s = sin(roll);
    double cc_ = c__ * _c_;
    double cs_ = c__ * _s_;
    double sc_ = s__ * _c_;
    double ss_ = s__ * _s_;
    double c_c = c__ * __c;
    double c_s = c__ * __s;
    double s_c = s__ * __c;
    double s_s = s__ * __s;
    double _cc = _c_ * __c;
    double _cs = _c_ * __s;
    double csc = cs_ * __c;
    double css = cs_ * __s;
    double ssc = ss_ * __c;
    double sss = ss_ * __s;
    Eigen::Matrix3d wRo;
    wRo <<
      cc_  , css-s_c,  csc+s_s,
      sc_  , sss+c_c,  ssc-c_s,
     -_s_  ,     _cs,      _cc;
    return wRo;
  }

  static void wRo_to_euler(const Eigen::Matrix3d& wRo, double& yaw, double& pitch, double& roll) {
    yaw = standardRad(atan2(wRo(1,0), wRo(0,0)));
    double c = cos(yaw);
    double s = sin(yaw);
    pitch = standardRad(atan2(-wRo(2,0), wRo(0,0)*c + wRo(1,0)*s));
    roll  = standardRad(atan2(wRo(0,2)*s - wRo(1,2)*c, -wRo(0,1)*s + wRo(1,1)*c));
  }

  static Eigen::Quaterniond wRo_to_quat(const Eigen::Matrix3d& wRo) {
    return Eigen::Quaterniond(wRo);
  }

  static Eigen::Matrix3d quat_to_wRo(const Eigen::Quaterniond& quat) {
    return Eigen::Matrix3d(quat);
  }

  static Eigen::Quaterniond euler_to_quat(double yaw, double pitch, double roll) {
    double sy = sin(yaw*0.5);
    double cy = cos(yaw*0.5);
    double sp = sin(pitch*0.5);
    double cp = cos(pitch*0.5);
    double sr = sin(roll*0.5);
    double cr = cos(roll*0.5);
    double w = cr*cp*cy + sr*sp*sy;
    double x = sr*cp*cy - cr*sp*sy;
    double y = cr*sp*cy + sr*cp*sy;
    double z = cr*cp*sy - sr*sp*cy;
    return Eigen::Quaterniond(w,x,y,z);
  }

  static void quat_to_euler(Eigen::Quaterniond q, double& yaw, double& pitch, double& roll) {
    const double q0 = q.w();
    const double q1 = q.x();
    const double q2 = q.y();
    const double q3 = q.z();
//    roll = atan2(2*(q0*q1+q2*q3), 1-2*(q1*q1+q2*q2));
    roll = atan2(2.0*(q0*q1+q2*q3), q0*q0-q1*q1-q2*q2+q3*q3); // numerically more stable (thanks to Dehann for pointing this out)
    pitch = asin(2.0*(q0*q2-q3*q1));
//    yaw = atan2(2*(q0*q3+q1*q2), 1-2*(q2*q2+q3*q3));
    yaw = atan2(2.0*(q0*q3+q1*q2), q0*q0+q1*q1-q2*q2-q3*q3);
  }

  static Eigen::Quaterniond delta3_to_quat(const Eigen::Vector3d& delta) {
    double theta = delta.norm();
    double S;
    if (theta < 0.0001) { // sqrt4(machine precession)
      S = 0.5 + theta*theta/48.;
    } else {
      S = sin(0.5*theta)/theta;
    }
    double C = cos(0.5*theta);
    return Eigen::Quaterniond(C, S*delta(0), S*delta(1), S*delta(2));
  }

  static const int dim = 3;
  static const char* name() {
    return "Rot3d";
  }

#ifdef USE_QUATERNIONS

  // new quaternion (and rotation matrix) based Rot3d

  // optimization is based on quaternions, but rotations matrices and
  // Euler angles can also be set or requested; efficiency is achieved
  // by only converting once needed in combination with caching

private:

  // quat always exists, while rotation matrix and Euler angles are cached
  Eigen::Quaterniond _quat;

  // the following variables can be changed by const methods: it is
  // essential that the class still logically behaves as const !!!!
  mutable bool _wRo_cached;
  mutable Eigen::Matrix3d _wRo;
  mutable bool _ypr_cached;
  mutable double _yaw;
  mutable double _pitch;
  mutable double _roll;

  // calculate ypr if it doesn't exist
  void ensure_ypr() const {
    if (!_ypr_cached) {
      quat_to_euler(_quat, _yaw, _pitch, _roll);
      _ypr_cached = true;
    }
  }

public:

  static const int size = 4;

  Rot3d() : _quat(Eigen::Quaterniond(1,0,0,0)),
    _wRo_cached(true), _wRo(Eigen::Matrix<double, 3, 3>::Identity()),
    _ypr_cached(true), _yaw(0.), _pitch(0.), _roll(0.) {}

  Rot3d(const Eigen::Quaterniond& quat) : _quat(quat), _wRo_cached(false), _ypr_cached(false) {}

  Rot3d(double yaw, double pitch, double roll) {
    set(yaw, pitch, roll);
  }

  Rot3d(const Eigen::Matrix3d& wRo) {
    _wRo = wRo;
    _wRo_cached = true;
    _ypr_cached = false;
    _quat = wRo_to_quat(wRo);
  }

  void write(std::ostream &out) const {
    out << yaw() << " " << pitch() << " " << roll() <<
        " (quat: " << _quat.w() << " " << _quat.x() << " " << _quat.y() << " " << _quat.z() << ")";
  }

  double x() const {return _quat.x();}
  double y() const {return _quat.y();}
  double z() const {return _quat.z();}
  double w() const {return _quat.w();}

  double yaw()   const {ensure_ypr(); return _yaw;}
  double pitch() const {ensure_ypr(); return _pitch;}
  double roll()  const {ensure_ypr(); return _roll;}

  void ypr(double& yaw, double& pitch, double& roll) const {
    ensure_ypr();
    yaw   = _yaw;
    pitch = _pitch;
    roll  = _roll;
  }

  void set_yaw  (double yaw)   {ensure_ypr(); set(yaw, _pitch, _roll);}
  void set_pitch(double pitch) {ensure_ypr(); set(_yaw, pitch, _roll);}
  void set_roll (double roll)  {ensure_ypr(); set(_yaw, _pitch, roll);}

  void set(double yaw, double pitch, double roll) {
    _yaw   = yaw;
    _pitch = pitch;
    _roll  = roll;
    _ypr_cached = true;
    _wRo_cached = false;
    _quat = euler_to_quat(yaw, pitch, roll);
  }


  Rot3d exmap(const Eigen::VectorXd& delta) const {
#if 1
    // direct solution by mapping to quaternion (following Grassia98jgt)
    Rot3d rot(_quat * delta3_to_quat(delta));
    return rot;
#else
    // cumbersome and slower solution
    double theta = delta.norm();
    double a, b, c;
    if (theta>1e-10) {
      a = delta(0)/theta;
      b = delta(1)/theta;
      c = delta(2)/theta;
    } else {
      a = b = c = 0.;
    }
    double S = sin(theta);
    double C = 1-cos(theta);
    Eigen::Matrix3d m;
    m <<
      1+(-a*a-b*b)*C, a*S-b*c*C, b*S+a*c*C,
      -a*S-b*c*C, 1+(-a*a-c*c)*C, c*S-a*b*C,
      -b*S+a*c*C, -c*S-a*b*C, 1+(-b*b-c*c)*C;
    return Rot3d(wRo() * m);
#endif
  }

  Eigen::Quaterniond quaternion() const {
    return _quat;
  }

  /**
   * Generate 3x3 rotation matrix from Rot3d.
   * @return wRo
   */
  const Eigen::Matrix3d& wRo() const {
    if (!_wRo_cached) {
      _wRo = quat_to_wRo(_quat);
      _wRo_cached = true;
    }
    return _wRo;
  }

  /**
   * Return inverse rotation by transeposing.
   * @return oRw
   */
  Eigen::Matrix3d oRw() const {
    return wRo().transpose();
  }

#else

  // old Euler-angle based Rot3d

private:

  double _yaw;
  double _pitch;
  double _roll;

public:

  static const int size = 3;

  Rot3d() : _yaw(0.), _pitch(0.), _roll(0.) {}

  Rot3d(double yaw, double pitch, double roll) : _yaw(yaw), _pitch(pitch), _roll(roll) {}

  /**
   * Initialize Euler angles from 3x3 wRo rotation matrix.
   */
  Rot3d(const Eigen::MatrixXd& wRo) {
    // note that getting the sign right requires recovering both sin and cos
    // for each angle; some equations exploit the fact that sin^2+cos^2=1
    _yaw = standardRad(atan2(wRo(1,0), wRo(0,0)));
    double c = cos(_yaw);
    double s = sin(_yaw);
    _pitch = standardRad(atan2(-wRo(2,0), wRo(0,0)*c + wRo(1,0)*s));
    _roll  = standardRad(atan2(wRo(0,2)*s - wRo(1,2)*c, -wRo(0,1)*s + wRo(1,1)*c));
  }

  void write(std::ostream &out) const {
    out << yaw() << " " << pitch() << " " << roll();
  }

  double yaw()   const {return _yaw;}
  double pitch() const {return _pitch;}
  double roll()  const {return _roll;}

  void ypr(double& yaw, double& pitch, double& roll) const {
    yaw = _yaw; pitch = _pitch; roll = _roll;
  }

  void set_yaw  (double yaw)   {_yaw = yaw;}
  void set_pitch(double pitch) {_pitch = pitch;}
  void set_roll (double roll)  {_roll = roll;}

  void set(double yaw, double pitch, double roll) {
    _yaw = yaw; _pitch = pitch; _roll = roll;
  }

  Rot3d exmap(const Eigen::VectorXd& delta) const {
    Rot3d res = *this;
    res._yaw   = standardRad(res._yaw   + delta(0));
    res._pitch = standardRad(res._pitch + delta(1));
    res._roll  = standardRad(res._roll  + delta(2));
    return res;
  }

  /**
   * Generate 3x3 rotation matrix from Rot3d.
   * @return wRo
   */
  Eigen::MatrixXd wRo() const {
    Eigen::MatrixXd ret = euler_to_wRo(_yaw, _pitch, _roll);
    return ret;
  }

#endif


};

}

