/**
 * @file glc.cpp
 * @brief Generic Linear Constraint functions
 * @author Nicholas Carlevaris-Bianco
 * @author Michael Kaess
 * @version $Id: covariance.cpp 4975 2011-07-13 17:49:09Z kaess $
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

#include <vector>

#include "isam/glc.h"
#include "isam/util.h"
#include "isam/Rot3d.h"
#include "isam/ChowLiuTree.h"

#define GLC_EPS 1e-8
#define GLC_INCLUDE_IC_FACTORS

using namespace std;
using namespace Eigen;

namespace isam {

#ifdef USE_QUATERNIONS
// code to account for angle axis to euler state represention change when using quaternions
// this would evaluate to identity if exmap wasn't mapping between euler and angle axis
// something similar can be use to trasform covariance recovery output when using quats
MatrixXd node_exmap_jacobian (Node *n) {
  // TODO replace with analytical jacobian
  const double eps = 0.0001;
  const double eps2 = 0.0002;
  int dim_n = n->dim();
  MatrixXd J(dim_n, dim_n);
  for(int col=0; col<dim_n; col++){
    VectorXd d(dim_n);
    d.setZero();
    // remember original value
    VectorXd o = n->vector0();
    // evaluate positive delta
    d(col) = eps;
    n->self_exmap(d);
    VectorXd yp = n->vector0();
    n->update0(o);
    // evaluate negative delta
    d(col) = -eps;
    n->self_exmap(d);
    VectorXd ym = n->vector0();
    n->update0(o);
    // calculte diff
    VectorXd diff = (yp - ym) / (eps2);
    // wrap angular difference
    VectorXb is_angle = n->is_angle();
    for (int k=0; k<dim_n; k++) {
      if (is_angle(k))
        diff(k) = standardRad(diff(k));
    }
    J.col(col) = diff;
  }
  return J;
}

MatrixXd exmap_jacobian (const vector<Node*>& nodes) {
  int dim = 0;
  for (size_t i=0; i<nodes.size(); i++) {
    dim += nodes[i]->dim();;
  }
  MatrixXd J = MatrixXd::Identity(dim, dim);
  int off=0;
  for (size_t i=0; i<nodes.size(); i++) {
    int dim_n = nodes[i]->dim();
    J.block(off,off,dim_n,dim_n) = node_exmap_jacobian (nodes[i]);
    off += dim_n;
  }
  return J;
}
#endif // USE_QUATERNIONS

MatrixXd glc_cholcov (MatrixXd A, double eps = numeric_limits<float>::epsilon()) {
    
  A = (A.transpose() + A) / 2.0;
  
  SelfAdjointEigenSolver<MatrixXd> eigensolver(A);
  if (eigensolver.info() != Success) {
    cout << "ERROR: Failed to compute eigenvalue decomposition!" << endl;
    cout << "A:" << endl << A << endl;
  }
       
  VectorXd D = eigensolver.eigenvalues();
  MatrixXd U = eigensolver.eigenvectors();
  
  double tol = eps * D.size();
  
  vector<int> inds_pos;
  for (int i=0; i<D.size(); i++) {
    if (D(i) > tol) 
      inds_pos.push_back(i);
  }
  
  if (0 == inds_pos.size()) {
      return MatrixXd();
  }
  
  VectorXd Dsub_sqrt(inds_pos.size());
  MatrixXd Usub(U.rows(), inds_pos.size());
  for (size_t i=0; i<inds_pos.size(); i++) {
    Dsub_sqrt(i) = sqrt(D(inds_pos[i]));
    Usub.col(i) = U.col(inds_pos[i]);
  }
  
  MatrixXd B;
  if (inds_pos.size() == (unsigned int)D.size()) {
    B = A.llt().matrixL().transpose(); 
  } else {
    B = Dsub_sqrt.asDiagonal() * Usub.transpose();
  }
  
#ifdef GLC_DEBUG
  double recreate_test = (B.transpose()*B - A).array().abs().matrix().lpNorm<Infinity>();
  cout << "[glc]\tCholcov Recreate Test: " << recreate_test << endl;
  if (recreate_test > 10) {
    cout << "ERROR: Cholcov Recreate Test Failed!" << endl;
    cout << "A" << endl << A << endl;
  }
#endif

  return B;   
}

MatrixXd glc_get_weighted_jacobian (isam::Factor *f) {

    Jacobian jac = f->jacobian();
    // get jacobian size
    std::vector<Node*>& f_nodes = f->nodes();
    int n_measure = f->dim();
    int n_var = 0;
    for (size_t i=0; i<f_nodes.size(); i++) {
        n_var += f_nodes[i]->dim();
    }
        
    // fill jacobian matrix
    MatrixXd H (n_measure, n_var);
    int offset = 0;
    for (Terms::const_iterator it=jac.terms().begin(); it!=jac.terms().end(); it++) {
        int ncols = it->term().cols();
        H.block(0, offset, n_measure, ncols) = it->term();
        offset += ncols;
    }
    
    return H;

}


std::vector<isam::Node*> glc_elim_clique_nodes (Node *node) {
    
  vector<isam::Node*> node_vector;
  
  int id_e = node->unique_id();
  
  std::list< isam::Node* > node_list;
  const list<Factor*>& factors = node->factors();
  for (list<Factor*>::const_iterator it = factors.begin(); it!=factors.end(); it++) {
      
    std::vector<Node*>& f_nodes = (*it)->nodes();
    
    for (size_t i=0; i<f_nodes.size(); i++) {
    
      int id_i = f_nodes[i]->unique_id();
      
      if (id_e != id_i) {
          
        // try to add second node in factor
        if (node_list.end() == find (node_list.begin(), node_list.end(), f_nodes[i])) { // not already added
          node_list.push_back(f_nodes[i]);
          node_vector.push_back (f_nodes[i]);
        }       
      }   
    }        
  }
  
  return node_vector;
}

vector<Factor*>
glc_intra_clique_factors (vector<Node*> clique_nodes, Node *node) {

  vector<Factor*> ic_factors;

  // loop over each node
  for (size_t i=0; i<clique_nodes.size(); i++) {
    // get this nodes factors
    const std::list<Factor*> factors = clique_nodes[i]->factors();

    for (list<Factor*>::const_iterator it = factors.begin(); it!=factors.end(); it++) {

      // make sure factor hasnt already been added to the list
      if (ic_factors.end() != find (ic_factors.begin(), ic_factors.end(), (*it)))
          continue;

      std::vector<Node*>& f_nodes = (*it)->nodes();

      // nodes in these factors can be: the marg node, nodes in the clique, or nodes ouside the clique
      // we wish to return factors strictly between nodes in the clique, not outside nor margnode
      // is the marg node in this factor
      if (f_nodes.end() != find (f_nodes.begin(), f_nodes.end(), node))
          continue; // node found

      // strctly include in clique
      bool ic = true;
      for (size_t j=0; j<f_nodes.size() && ic; j++) {
        // if we find a factor node that is not in the clique
        if (clique_nodes.end() == find (clique_nodes.begin(), clique_nodes.end(), f_nodes[j]))
          ic = false;
      }

      if (ic)
        ic_factors.push_back(*it);

    }
  }

  return ic_factors;

}

MatrixXd glc_target_info (Node *node, vector<Node*>& clique_nodes, vector<Factor*>& ic_factors){

  vector<Node*> all_nodes (clique_nodes);
  all_nodes.push_back(node);
  
  int n_full = 0;
  for (size_t i=0; i<all_nodes.size(); i++) n_full += all_nodes[i]->dim();
  MatrixXd L (n_full, n_full); // clique nodes first then marg node at end
  L.setZero();
  
  const list<Factor*>& tmp = node->factors();
  list<Factor*> factors (tmp);
  for (size_t i=0; i<ic_factors.size(); i++)
    factors.push_back(ic_factors[i]);
  
  for (list<Factor*>::iterator it = factors.begin(); it!=factors.end(); it++) {
    
    Factor *f = *it;

    ////HACK TO THROW AWAY Z PRIORS (for comparision with compounding that also throws away z-priors)
    //if (0 == strcmp("Pose3d_z_Factor", f->name())) {
    //  continue;
    //}

    MatrixXd H = glc_get_weighted_jacobian (f);
    MatrixXd dL = H.transpose()*H;
    
    // there is probably a better/faster way to do this
    vector<Node*>& f_nodes = f->nodes();
    int ioff = 0;
    for (size_t i=0; i<f_nodes.size(); i++) {
      Node *ni = f_nodes[i];
      int joff = ioff;
      for (size_t j=i; j<f_nodes.size(); j++) {
        Node *nj = f_nodes[j];
        
        // find the position(s) of this block in the target information
        int koff = 0;
        Node *nk = NULL;
        for (size_t k=0; k<all_nodes.size(); k++) {
          if (ni == all_nodes[k]) {
            nk = all_nodes[k];
            break;
          }
          koff += all_nodes[k]->dim();
        }
        int loff = 0;
        Node *nl = NULL;
        for (size_t l=0; l<all_nodes.size(); l++) {
          if (nj == all_nodes[l]) {
            nl = all_nodes[l];
            break;
          }
          loff += all_nodes[l]->dim();
        }
        // add information to the target informaiton
        L.block(koff, loff, nk->dim(), nl->dim()) += dL.block(ioff, joff, ni->dim(), nj->dim());
        if (ni != nj) // off diagonal elements 
          L.block(loff, koff, nl->dim(), nk->dim()) += dL.block(joff, ioff, nj->dim(), ni->dim());
        
        joff += nj->dim();  
      }
      ioff += ni->dim();
    }
    
#ifdef GLC_DEBUG
    cout << "[glc]\tAdding info from " << f->name() << " factor between nodes : ";
    for (size_t n=0; n<f_nodes.size(); n++)
        cout << f_nodes[n]->unique_id() << " ";
    cout << endl;
    //cout << "L" << endl << L << endl;
#endif

  }
  
  // marginalization
  int n = n_full - node->dim();
  MatrixXd Lbb = L.bottomRightCorner(node->dim(), node->dim());
  MatrixXd L_marg = L.topLeftCorner(n,n) - (L.topRightCorner(n,node->dim()) * posdef_pinv(Lbb, GLC_EPS) * L.bottomLeftCorner(node->dim(),n));
  //cout << "L_marg" << endl << L_marg << endl;
  return L_marg;

}

Factor* glc_factor (const MatrixXd& L, const vector<Node*>& clique_nodes, GLC_Reparam* rp) {
    
    int np = clique_nodes.size();
    int n = 0;
    for (int i=0; i<np; i++) n += clique_nodes[i]->dim();

    MatrixXd Lp = L;
#ifdef USE_QUATERNIONS
    MatrixXd Jexmap = exmap_jacobian (clique_nodes);
    MatrixXd Jpinv = pinv(Jexmap, GLC_EPS);
    Lp = Jpinv.transpose() * L * Jpinv;
#endif
        
    Factor *f_add = NULL;
    if (rp == NULL) { // H with linearization on world frame xfms

      MatrixXd G = glc_cholcov(Lp, GLC_EPS);

      VectorXd x(n);
      x.setZero();
      int off = 0;
      for (int i=0; i<np; i++) {
        x.segment(off, clique_nodes[i]->dim()) = clique_nodes[i]->vector(LINPOINT);
        off += clique_nodes[i]->dim();
      }

      if (G.rows() != 0 && G.cols() != 0) {
        f_add = new GLC_Factor(clique_nodes, x, G, NULL);
      } else {
        return NULL; // can happen when removing node at end of chain
      }
    
    } else { // H with linearization on root shifted xfms

      GLC_Reparam *rptmp = rp->clone();
      rptmp->set_nodes (clique_nodes);
      VectorXd x_rp = rptmp->reparameterize(LINPOINT);
      MatrixXd F = rptmp->jacobian();
      delete rptmp;

      MatrixXd Fpinv = pinv(F, GLC_EPS);
      MatrixXd UTU = Fpinv.transpose() * Lp * Fpinv;

      MatrixXd G;
      G = glc_cholcov(UTU, GLC_EPS);

      if(G.rows() == 0) {
#ifdef GLC_DEBUG
        cout << "[GLC]\t\tComputed No Rank Factor! Skipping" << endl;
#endif
        if (Lp.rows() <= 6) {
          return NULL;
        } else {
          cout << "Lp" << endl << Lp << endl;
          exit(0);
        }
      }

#ifdef GLC_DEBUG
      MatrixXd H = G*F;
      double recreate_test = (H.transpose()*H - Lp).array().abs().matrix().lpNorm<Infinity>();
      cout << "[glc]\tTarget Recreate Test: " << recreate_test << endl;
      if (recreate_test > 10) {
        cout << "[ERROR]\tRecreate Test Failed!" << endl;
        cout << "Lp" << endl << Lp << endl;
        cout << "F" << endl << F << endl;
        exit(0);
      }
#endif

      f_add = new GLC_Factor(clique_nodes, x_rp, G, rp);
    }

#ifdef GLC_DEBUG
    if (f_add != NULL) {
      cout << "[glc]\t\tComputed GLC_Factor " << f_add->unique_id() << " between nodes ";
      for (int i=0; i<np; i++)
          cout << clique_nodes[i]->unique_id() << " ";
      cout << endl;
    }
#endif
    
    return f_add;
    
}

vector<Factor*> glc_lift_factors (const MatrixXd& L, const vector<Node*>& clique_nodes, bool sparse,
                                  GLC_Reparam* rp) {
  
  vector<Factor*> new_glc_factors;    

  if (sparse) {

    ChowLiuTree clt(L, clique_nodes);
    
    map<int, ChowLiuTreeNode>::iterator it;
    for (it=clt.tree.begin(); it!=clt.tree.end(); it++) {
      
      if (it->second.is_root()) {
        MatrixXd La = it->second.marginal;
        double mag_test = La.array().abs().matrix().lpNorm<Infinity>();
        if (mag_test > GLC_EPS) {
          vector<Node*> f_nodes;
          f_nodes.push_back(clique_nodes[it->second.id]);
          Factor * f_add = glc_factor (La, f_nodes, rp);
          if (f_add != NULL) {
            new_glc_factors.push_back (f_add);
          }
        }
      } else {
        vector<Node*> f_nodes;
        f_nodes.push_back(clique_nodes[it->second.id]);
        f_nodes.push_back(clique_nodes[it->second.pid]);
        int adim = f_nodes[0]->dim();
        int bdim = f_nodes[1]->dim();
        MatrixXd Lagb = it->second.conditional;
        MatrixXd Lab = it->second.joint;
        MatrixXd Hagb(adim, adim+bdim);
        Hagb.block(0,0,adim,adim) = MatrixXd::Identity (adim,adim);
        MatrixXd Laa = Lab.block(0,0,adim,adim);
        Hagb.block(0,adim,adim,bdim) = posdef_pinv(Laa) * Lab.block(0,adim,adim,bdim);
        MatrixXd Lt = Hagb.transpose() * Lagb * Hagb;
        Factor * f_add = glc_factor (Lt, f_nodes, rp);
        if (f_add != NULL)
            new_glc_factors.push_back (f_add);
      }
    }
    
  } else {
    Factor * f_add = glc_factor (L, clique_nodes, rp);
    if (f_add != NULL)
      new_glc_factors.push_back (f_add);
  }
      
  return new_glc_factors;
  
}

vector<Factor*> glc_elim_factors (Node* node) {

  vector<Factor*> elim_factors (node->factors().begin(), node->factors().end());
#ifdef GLC_INCLUDE_IC_FACTORS
  vector<Factor*> ic_factors;
  ic_factors = glc_intra_clique_factors (glc_elim_clique_nodes (node), node);
  elim_factors.insert(elim_factors.end(), ic_factors.begin(), ic_factors.end());
#endif
  return elim_factors;

}



vector<Factor*> glc_remove_node(Slam& slam, Node* node, bool sparse, GLC_Reparam* rp) {
  
  // get the nodes in the elimination clique, new glc factor(s) will span these nodes
  vector<Node*> clique_nodes = glc_elim_clique_nodes (node);
    
#ifdef GLC_DEBUG
  cout << "[glc]\tRemoving Node: " << node->unique_id() << endl;
  cout << "[glc]\tElimination Clique Nodes: ";
  for (size_t i=0; i<clique_nodes.size(); i++)
    cout << clique_nodes[i]->unique_id() << " ";
  cout << endl;
#endif
  
  // intra-clique factors, those in the clique that are not directly connected to the marg node
  // not strictly required, initial icra paper included them
  vector<Factor*> ic_factors;
#ifdef GLC_INCLUDE_IC_FACTORS
  ic_factors = glc_intra_clique_factors (clique_nodes, node);
#endif
  MatrixXd L = glc_target_info (node, clique_nodes, ic_factors);

  // lift glc factors
  vector<Factor*> new_glc_factors;
  new_glc_factors = glc_lift_factors (L, clique_nodes, sparse, rp);

  // remove node and delete all adjacent factors
  slam.remove_node(node);  

  // remove all ic factors
  for(size_t i=0; i<ic_factors.size(); i++) {
    slam.remove_factor(ic_factors[i]);
  }
  
  // add glc factors
  for(size_t i=0; i<new_glc_factors.size(); i++) {
    slam.add_factor(new_glc_factors[i]);
#ifdef GLC_DEBUG
    cout << "[glc]\tAdded GLC Factor: " << new_glc_factors[i]->unique_id() << endl;
#endif
  }
  

  return new_glc_factors;

}

} // namespace isam
