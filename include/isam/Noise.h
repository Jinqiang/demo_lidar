/**
 * @file Noise.h
 * @brief Various noise models.
 * @author Michael Kaess
 * @version $Id: Noise.h 5797 2011-12-07 03:50:41Z kaess $
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
#include <Eigen/LU> 

namespace isam {

// general noise model class
class Noise {
public:
  Eigen::MatrixXd _sqrtinf;
  const Eigen::MatrixXd& sqrtinf() const {return _sqrtinf;}
};

// noise model based on square root information matrix
class SqrtInformation : public Noise {
public:
  SqrtInformation(const Eigen::MatrixXd& sqrtinf) {_sqrtinf = sqrtinf;}
};

// noise model based on information matrix
class Information : public Noise {
public:
  Information(const Eigen::MatrixXd& inf) {
    _sqrtinf = inf.llt().matrixU();
  }
};

// noise model based on covariance matrix
class Covariance : public Noise {
public:
  Covariance(const Eigen::MatrixXd& cov) {
    _sqrtinf = cov.inverse().llt().matrixU();
  }
};

}
