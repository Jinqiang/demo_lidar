/**
 * @file Element.h
 * @brief Basic functionality for nodes and factors of a graph.
 * @author Michael Kaess
 * @version $Id: Element.h 5796 2011-12-07 00:39:30Z kaess $
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
#include <ostream>

namespace isam {

class Element {
  Element(const Element& rhs); // not allowed
  const Element& operator= (const Element& rhs); // not allowed

  const char* _name;

  int _start; // needed for Slam::jacobian

protected:
  int _id;
  int _dim;

public:
  Element(const char* name, int dim) : _name(name), _dim(dim) {}
  virtual ~Element() {};

  virtual int unique_id() {return _id;}
  virtual const char* name() const {return _name;}
  inline int dim() const {return _dim;}
  inline int start() const {return _start;}

  virtual void write(std::ostream &out) const {
    out << name();
  }

  friend class Slam;
  friend class Covariances;
};

}
