/**
 * @file Node.cpp
 * @brief A single variable node
 * @author Michael Kaess
 * @author Hordur Johannsson
 * @version $Id:  $
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

#include "isam/Node.h"
#include "isam/Factor.h"

namespace isam
{
void Node::mark_deleted() {
  _deleted = true;
  for (std::list<Factor*>::iterator factor = _factors.begin(); factor != _factors.end(); ++factor)
    (*factor)->mark_deleted();
}
void Node::erase_marked_factors() {
  for (std::list<Factor*>::iterator factor = _factors.begin(); factor != _factors.end();)
    if ((*factor)->deleted()) factor = _factors.erase(factor);
    else ++factor;
}
} // namespace - isam
