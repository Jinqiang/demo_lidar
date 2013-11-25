/**
 * @file util.cpp
 * @brief Basic utility functions that are independent of iSAM.
 * @author Michael Kaess
 * @version $Id: util.cpp 6335 2012-03-22 23:13:52Z kaess $
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

#include <map>
#include <sys/time.h>
#include <cmath>
#include <algorithm> // abs

#include "isam/util.h"

using namespace std;

namespace isam {

// simple class for accumulating execution timing information by name
class Timing {
  class Stats {
  public:
    double t0;
    double t;
    double t_max;
    double t_min;
    int n;
  };
  map<string, Stats> stats;
public:
  void add_t0(string id, double t0) {
    stats[id].t0 = t0;
  }
  double get_t0(string id) {
    return stats[id].t0;
  }
  void add_dt(string id, double dt) {
    Stats& s = stats[id];
    s.t += dt;
    s.n++;
    if (s.n==1 || s.t_max < dt) s.t_max = dt;
    if (s.n==1 || s.t_min > dt) s.t_min = dt;
  }
  void print() {
    map<string, Stats>::iterator it;
    for(it = stats.begin(); it!=stats.end(); it++) {
      Stats& s = it->second;
      printf("%s: %g (%i times, min: %g, max: %g)\n",
             it->first.c_str(), s.t, s.n, s.t_min, s.t_max);
    }
  }
  double time(string id) {
    Stats& s = stats[id];
    return s.t;
  }
};
Timing timing;

double tic() {
  struct timeval t;
  gettimeofday(&t, NULL);
  return ((double)t.tv_sec + ((double)t.tv_usec)/1000000.);
}

double tic(string id) {
  double t0 = tic();
  timing.add_t0(id, t0);
  return t0;
}

double toc(double t) {
  double s = tic();
  return (max(0., s-t));
}

double toc(string id) {
  double dt = toc(timing.get_t0(id));
  timing.add_dt(id, dt);
  return dt;
}

void tictoc_print() {
  timing.print();
}

double tictoc(string id) {
  return (timing.time(id));
}

Eigen::MatrixXd eye(int num) {
  return Eigen::MatrixXd::Identity(num, num);
}

void givens(const double a, const double b, double& c, double& s) {
  if (b==0) {
    c = 1.0;
    s = 0.0;
  } else if (fabs(b)>fabs(a)) {
    double t = -a/b;
    s = 1/sqrt(1+t*t);
    c = t*s;
  } else {
    double t = -b/a;
    c = 1/sqrt(1+t*t);
    s = t*c;
  }
}

}
