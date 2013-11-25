/**
 * @file Anchor.cpp
 * @brief Implementation of the anchor node.
 * @author Hordur Johannsson
 * @author Michael Kaess
 * @version $Id: Anchor.cpp 6335 2012-03-22 23:13:52Z kaess $
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

#include <isam/Anchor.h>
#include <isam/slam2d.h>
#include <isam/slam3d.h>

namespace isam
{

Anchor2d_Node::Anchor2d_Node(Slam* slam) : Pose2d_Node("Anchor2d"), _parent(NULL), _prior(NULL), _slam(slam)
{
}

Anchor2d_Node::~Anchor2d_Node() {
  if (_prior) {
    _slam->remove_factor(_prior);
    delete _prior;
  }
}

void Anchor2d_Node::set_prior() {
  Noise noise = SqrtInformation(10000. * eye(3));
  Pose2d prior_origin(0., 0., 0.);
  _prior = new Pose2d_Factor(this, prior_origin, noise);
  _slam->add_factor(_prior);
}

void Anchor2d_Node::add_anchor(Anchor2d_Node* a) {
  if (a->_prior) {
    _slam->remove_factor(a->_prior);
    delete a->_prior;
    a->_prior = NULL;
  }
  a->_parent = this;
  _childs.push_back(a);
}

void Anchor2d_Node::merge(Anchor2d_Node* a, Pose2d old_origin) {
  if (_parent != NULL) {
    _parent->merge(a, old_origin);
  } else if (a->_parent != NULL) {
    merge(a->parent(), old_origin);
  } else {
    // this and a are root anchors.

    // Add a to the child list
    Pose2d old_anchor_pose = a->value();
    a->init(old_origin.oplus(old_anchor_pose));
    add_anchor(a);

    // Move the children of a to this
    for (size_t i = 0; i < a->_childs.size(); i++) {
      Anchor2d_Node* child = a->_childs[i];
      
      old_anchor_pose = child->value();
      child->init(old_origin.oplus(old_anchor_pose));
      add_anchor(child);
    }
    a->_childs.clear();
  }
}

Anchor3d_Node::Anchor3d_Node(Slam* slam) : Pose3d_Node("Anchor3d"), _parent(NULL), _prior(NULL), _slam(slam)
{
}

Anchor3d_Node::~Anchor3d_Node() {
  if (_prior) {
    _slam->remove_factor(_prior);
    delete _prior;
  }
}

void Anchor3d_Node::set_prior() {
  Noise noise = SqrtInformation(10000. * eye(6));
  Pose3d prior_origin;
  _prior = new Pose3d_Factor(this, prior_origin, noise);
  _slam->add_factor(_prior);
}

void Anchor3d_Node::add_anchor(Anchor3d_Node* a) {
  if (a->_prior) {
    _slam->remove_factor(a->_prior);
    delete a->_prior;
    a->_prior = NULL;
  }
  a->_parent = this;
  _childs.push_back(a);
}

void Anchor3d_Node::merge(Anchor3d_Node* a, Pose3d old_origin) {
  if (_parent != NULL) {
    _parent->merge(a, old_origin);
  } else if (a->_parent != NULL) {
    merge(a->parent(), old_origin);
  } else {
    // this and a are root anchors.

    // Add a to the child list
    Pose3d old_anchor_pose = a->value();
    a->init(old_origin.oplus(old_anchor_pose));
    add_anchor(a);

    // Move the children of a to this
    for (size_t i = 0; i < a->_childs.size(); i++) {
      Anchor3d_Node* child = a->_childs[i];

      old_anchor_pose = child->value();
      child->init(old_origin.oplus(old_anchor_pose));
      add_anchor(child);
    }
    a->_childs.clear();
  }
}

}
