/**
 * @file isam.h
 * @brief Main isam header file.
 * @author Michael Kaess
 * @version $Id: isam.h 6378 2012-03-30 20:42:53Z kaess $
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

// all include files needed by a user of iSAM

#include <isam/Slam.h>
#include <isam/Anchor.h>
#include <isam/slam2d.h>
#include <isam/slam3d.h>
#include <isam/Point2d.h>
#include <isam/Pose2d.h>
#include <isam/Point3d.h>
#include <isam/Pose3d.h>
#include <isam/util.h>
#include <isam/Noise.h>
#include <isam/Properties.h>


// Doxygen documentation follows

/** @mainpage iSAM: Incremental Smoothing and Mapping (v1.6)


@section intro_what What is iSAM?

iSAM is an optimization library for sparse nonlinear problems as
encountered in simultaneous localization and mapping (SLAM). The iSAM
library provides efficient algorithms for batch and incremental
optimization, recovering the exact least-squares solution. The library
can easily be extended to new problems, and functionality for often
encountered 2D and 3D SLAM problems is already provided. The iSAM
algorithm was originally developed by <a
href="http://people.csail.mit.edu/kaess" target="_top">Michael
Kaess</a> (kaess@mit.edu) and <a
href="http://www.cc.gatech.edu/~dellaert" target="_top">Frank
Dellaert</a> (dellaert@cc.gatech.edu) at Georgia Tech.


@section intro_copyright Copyright and License

Copyright (C) 2009-2012 Massachusetts Institute of Technology.<br> <a
href="http://people.csail.mit.edu/kaess" target="_top">Michael
Kaess</a> (kaess@mit.edu), Hordur Johannsson (hordurj@mit.edu), David
Rosen (dmrosen@mit.edu) and <a href="http://cml.mit.edu/~jleonard/"
target="_top">John Leonard</a> (jleonard@mit.edu)

iSAM is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation; either version 2.1 of the License, or (at
your option) any later version.

iSAM is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public
License along with iSAM.  If not, see <http://www.gnu.org/licenses/>.


@section intro_download Download

The latest version of iSAM is available at <a
href="http://people.csail.mit.edu/kaess/isam"
target="_top">http://people.csail.mit.edu/kaess/isam</a>

The source code is available from our public subversion repository:
@verbatim
svn co https://svn.csail.mit.edu/isam
@endverbatim


@section intro_requirements Requirements

This software was tested under Ubuntu 9.04-11.04, and Mac OS X
10.5-10.6. iSAM depends on the following software:

@li CMake (version 2.8 or higher)
@li CHOLMOD (sparse matrix libraries by <a href="http://www.cise.ufl.edu/~davis/" target="_top">Tim Davis</a>, University of Florida)
@li Eigen (version 3)
@li SDL (optional, needed for 3D viewer)
@li doxygen (optional, to generate documentation)
@li graphviz (optional, to generate documentation)

To install all required packages in Ubuntu 9.10 and later, simply type:

@verbatim
sudo apt-get install cmake libsuitesparse-dev libeigen3-dev libsdl1.2-dev doxygen graphviz
@endverbatim

Note that CHOLMOD is contained in SuiteSparse. On Mac OS X, SuiteSparse
and SDL are available through MacPorts:

@verbatim
sudo port install suitesparse libsdl
@endverbatim


@section intro_installation Installation

Compile with:
@verbatim
make
@endverbatim

Directory structure:
@li @c isamlib/ contains the source code for the iSAM library
@li @c include/ contains the header files for the iSAM library
@li @c isam/ contains the source code for the iSAM executable
@li @c examples/ contains examples for using different components of the iSAM library
@li @c doc/ contains the documentation (after calling "make doc")
@li @c misc/ contains code referenced from publications
@li @c data/ contains example data files
@li @c lib/ contains the actual library "libisam.a"
@li @c bin/ contains the main executable "isam"

Usage example:
@verbatim
bin/isam -G data/sphere400.txt
@endverbatim

For more usage information:
@verbatim
bin/isam -h
@endverbatim

Install the library in your system with:
@verbatim
make install
@endverbatim

Note that make just provides a convenient wrapper for running cmake in
a separate "build" directory. Compile options can be changed with
"ccmake build". In particular, support for the 3D viewer can be
disabled by setting USE_GUI to OFF. Library and include paths can be
modified manually in case SuiteSparse/CHOLMOD was installed in a local
directory and cannot be found automatically.

@section intro_documentation Documentation

To generate detailed documentation for the source code, type:
@verbatim
make doc
@endverbatim
and open @c doc/html/index.html in your browser.


@section intro_references References (please cite when using this software)

Details of the algorithms used in this software are provided in these
publications (the latex bibliography file @c isam.bib is included for
convenience).

@li "iSAM: Incremental Smoothing and Mapping" by M. Kaess,
A. Ranganathan, and F. Dellaert, IEEE Trans. on Robotics, TRO,
vol. 24, no. 6, Dec. 2008, pp. 1365-1378, <a
href="http://www.cc.gatech.edu/~kaess/pub/Kaess08tro.pdf"
target="_top">PDF</a>

@li "Covariance Recovery from a Square Root Information Matrix for
Data Association" by M. Kaess and F. Dellaert, Journal of Robotics and
Autonomous Systems, RAS, vol. 57, Dec. 2009, pp. 1198-1210, <a
href="http://www.cc.gatech.edu/~kaess/pub/Kaess09ras.pdf"
target="_top">PDF</a>

@li "An Incremental Trust-Region Method for Robust Online Sparse
Least-Squares Estimation" by D.M. Rosen, M. Kaess and J.J. Leonard,
International Conference on Robotics and Automation (ICRA), May 2012,
<a href="http://www.cc.gatech.edu/~kaess/pub/Rosen12icra.pdf"
target="_top">PDF</a>

A full list of iSAM-related references in BibTeX format is available
here: \ref Bibliography.

Newer publications will be available from my web page at <a
href="http://people.csail.mit.edu/kaess/"
target="_top">http://people.csail.mit.edu/kaess/</a>


@section intro_ack Acknowledgments

Many thanks to Richard Roberts for his help with this software. Thanks
also to John McDonald, Ayoung Kim, Ryan Eustice, Aisha Walcott, Been
Kim and Abe Bachrach for their feedback and patience with earlier
versions.

*/
