//===================================================================
//
// Author: Igor Chollet
//
//  This file is part of ankh.
//
//  ankh is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  ankh is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//  (see ../LICENSE.txt)
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with ankh.  If not, see <http://www.gnu.org/licenses/>
//
//====================================================================
#ifndef ANKH_GRID_HPP
#define ANKH_GRID_HPP

#include <vector>
#include <iostream>
#include <cmath>
#include "multipoles.hpp"

namespace ankh{

  inline int idx_in_grid(int i, int j, int k, int E){return i*E*E + j*E + k;}

  struct corresp{
    int true_idx;
    int part_idx;
  };

  struct part{
    double x;
    double y;
    double z;
  };

  template<typename charge_type>
  struct leaf{
    int nprt;
    std::vector<corresp> crsp;
    std::vector<part   > prts;
    std::vector<double > Qs;
    std::vector<double > Ps;
    part ctr;
  };

  template<typename charge_type>
  void allocate_preset_grid(int lvl, leaf<charge_type>*& grid,
			    double& locrad, double& cx, double& cy, double& cz){
    grid = new leaf<charge_type>[(1<<lvl)*(1<<lvl)*(1<<lvl)];
    int half_leaf_1D = (1<<(lvl-1));
    int n_leaf_1D    = (1<<lvl);
    for(int i = 0; i < n_leaf_1D; i++){
      double clx = double(1+2*(i-half_leaf_1D))*locrad + cx;
      for(int j = 0; j < n_leaf_1D; j++){
	double cly = double(1+2*(j-half_leaf_1D))*locrad + cy;
	for(int k = 0; k < n_leaf_1D; k++){
	  double clz = double(1+2*(k-half_leaf_1D))*locrad + cz;
	  int idx = idx_in_grid(i,j,k,n_leaf_1D);
	  grid[idx].nprt  = 0;
	  grid[idx].ctr.x = clx;
	  grid[idx].ctr.y = cly;
	  grid[idx].ctr.z = clz;
	}
      }
    }
  }

  template<typename charge_type>
  void compute_corresps(double& x, double& y, double& z,
			int i, leaf<charge_type>* leaves, int lvl,
			double& cx, double& cy, double& cz, double& rad){
    int E  = (1<<lvl);
    double leaf_rad = rad/double(E);
    int ix = floor((x-cx+rad)/rad*0.5*double(E));
    int iy = floor((y-cy+rad)/rad*0.5*double(E));
    int iz = floor((z-cz+rad)/rad*0.5*double(E));
    if(ix == (1<<lvl)){ix -= 1;} if(iy == (1<<lvl)){iy -= 1;} if(iz == (1<<lvl)){iz -= 1;}
    int iig = idx_in_grid(ix,iy,iz,E);
    corresp co;
    co.true_idx = i;
    co.part_idx = leaves[iig].nprt;
    leaves[iig].nprt++;
    leaves[iig].crsp.push_back(co);
    part prt;
    prt.x = x;
    prt.y = y;
    prt.z = z;
    leaves[iig].prts.push_back(prt);
  }
  
  template<typename charge_type>
  void load_charges(leaf<charge_type>* leaves, int lvl, charge_type* chrgs){
    int nleaves = (1<<lvl)*(1<<lvl)*(1<<lvl);
    int size    = n_terms<charge_type>();
    for(int ii = 0; ii < nleaves; ii++){
      leaves[ii].Qs.clear();
      leaves[ii].Qs.resize(leaves[ii].nprt*size);
      for(int jj = 0; jj < leaves[ii].nprt; jj++){
	copy_to_FLT_ptr(chrgs[leaves[ii].crsp[jj].true_idx],&(leaves[ii].Qs[jj*size]));
      }
    }
  }
  
  template<typename charge_type>
  void write_charges(leaf<charge_type>* leaves, int lvl, charge_type* chrgs){
    int nleaves = (1<<lvl)*(1<<lvl)*(1<<lvl);
    int size    = n_terms<charge_type>();
    for(int ii = 0; ii < nleaves; ii++){
      for(int jj = 0; jj < leaves[ii].nprt; jj++){
	copy_from_FLT_ptr(chrgs[leaves[ii].crsp[jj].true_idx],&(leaves[ii].Ps[jj*size]));
      }
    }
  }

} // ANKH

#endif
