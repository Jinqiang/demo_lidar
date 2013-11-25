/**
 * @file Anchor.h
 * @brief Implementation of the anchor node.
 * @author Hordur Johannsson
 * @author Michael Kaess
 * @version $Id: Anchor.h 6334 2012-03-22 18:53:24Z hordurj $
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

#include "Node.h"
#include "Factor.h"
#include "Pose2d.h"
#include "Point2d.h"
#include "Pose3d.h"
#include "Slam.h"

namespace isam {

typedef NodeT<Pose2d> Pose2d_Node;
typedef NodeT<Pose3d> Pose3d_Node;

/**
 * A anchor node for 2d poses
 */
class Anchor2d_Node : public Pose2d_Node {
public:
  Anchor2d_Node(Slam* slam);
  ~Anchor2d_Node();

  /**
   * Add a prior to the anchor. The prior will be removed
   * if this anchor is merged with another anchor's frame.
   */
  void set_prior();

  /**
   * Add a new anchor to this frame.
   */
  void add_anchor(Anchor2d_Node* a);

  /**
   * Merges the node with anchor a.
   *
   * @param a the node to merge with.
   * @param old_origin the pose of this frame in the new frame.
   *
   * Usage: b.merge(a);
   *
   * All anchors in a's frame will be merged with the b's frame.
   *
   */
  void merge(Anchor2d_Node* a, Pose2d old_origin);

  Anchor2d_Node* parent() { return _parent; }

private:
  Anchor2d_Node* _parent;
  std::vector<Anchor2d_Node*> _childs;
  Factor* _prior;
  Slam* _slam;
};

/**
 * A anchor node for 3d poses
 */
class Anchor3d_Node : public Pose3d_Node {
public:
  Anchor3d_Node(Slam* slam);
  ~Anchor3d_Node();

  /**
   * Add a prior to the anchor. The prior will be removed
   * if this anchor is merged with another anchor's frame.
   */
  void set_prior();

  /**
   * Add a new anchor to this frame.
   */
  void add_anchor(Anchor3d_Node* a);

  /**
   * Merges the node with anchor a.
   *
   * @param a the node to merge with.
   * @param old_origin the pose of this frame in the new frame.
   *
   * Usage: b.merge(a);
   *
   * All anchors in a's frame will be merged with the b's frame.
   *
   */
  void merge(Anchor3d_Node* a, Pose3d old_origin);

  Anchor3d_Node* parent() { return _parent; }

private:
  Anchor3d_Node* _parent;
  std::vector<Anchor3d_Node*> _childs;
  Factor* _prior;
  Slam* _slam;
};

}
