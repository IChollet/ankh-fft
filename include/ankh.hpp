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
#ifndef ANKH_HPP
#define ANKH_HPP

#include "grid.hpp"
#include "short.hpp"
#include "circulant.hpp"
#include "interpolation.hpp"
#include "multipoles.hpp"
#include "wrapper.hpp"
#include "colorandtime.hpp"
#include <iostream>

#define COULOMB 0
#define EWALD 1

namespace ankh{
    
  template<typename charge_type, int KRNL>
  class matrix{
  private:
    int                 lvl;            // Number of tree level
    int                 nprt;           // Number of particles
    int                 Lu;             // Interpolation order at leaves
    int                 Li;             // Interpolation order for far images vectors
    int                 Nimages;        // Number of images in each directions
    leaf<charge_type>*  grid;           // Grid of leaves
    double              grid_leaf_rad;  // Raidus of leaves
    double              box_rad;        // Radius of the original box
    corresp*            corresps;       // Input correspondances
    double **           U;              // Interpolation matrices
    double *            diag;           // Diagonal matrix in the Fourier domain
  public :
    matrix(int Lu_, int lvl_, int Nimages_, int Li_){
      lvl = lvl_; Lu  = Lu_; Nimages = Nimages_; Li = Li_;}
    ~matrix(){}
    
    void set_kernel_params(const double& beta_){kernel::beta = beta_;}

    void prcmpt(double cx, double cy, double cz, double rad, double* parts, int N, int Lc){
      // Initialize grid -> O(N)
      double t0, t1;
      t0 = chronos::time();
      box_rad = rad;
      nprt    = N;
      grid_leaf_rad = box_rad/(1<<lvl);
      allocate_preset_grid<charge_type>(lvl,grid,grid_leaf_rad,cx,cy,cz);
      for(int i = 0; i < N; i++){
        compute_corresps<charge_type>(parts[i*3],parts[i*3+1],parts[i*3+2],
				      i,grid,lvl,cx,cy,cz,rad);
      }
      // Compute U matrix (unif interpolation) -> O(N)
      int Lu3 = Lu*Lu*Lu;
      U = new double*[(1<<lvl)*(1<<lvl)*(1<<lvl)];
      double *toUnif;
      double *Trs = new double[Lc*Lc];
      cheb_to_unif_1d_matrix(Lu, Lc, toUnif);
      getTr(Lc, Trs);
      int mult = n_terms<charge_type>();
      for(int ii = 0; ii < (1<<lvl)*(1<<lvl)*(1<<lvl); ii++){
        U[ii] = new double[Lu3*grid[ii].nprt*mult];
        for(int jj = 0; jj < grid[ii].nprt; jj++){
	  part_to_cheb_to_unif<charge_type>(
                 (grid[ii].prts[jj].x-grid[ii].ctr.x),
		 (grid[ii].prts[jj].y-grid[ii].ctr.y), 
		 (grid[ii].prts[jj].z-grid[ii].ctr.z),
		 U[ii]+jj*Lu3*mult, Lu, Lc, toUnif, Trs, grid_leaf_rad);
        }
      }
      t1 = chronos::time();
      std::cout << "Compute polynomials: " << t1-t0 << std::endl;
      // Diagonalize grid matrix in the fourier domain -> O(NlogN)
      prcmp_diag<KRNL>((1<<lvl),Lu,box_rad,lvl,diag,Nimages,Li);
    }

    double apply(charge_type* q, double* parts, int Nm){
      double t0, t1;
      t0 = chronos::time();
      // Get sorted charges -> O(N)
      load_charges<charge_type>(grid,lvl,q);
      t1 = chronos::time();
      std::cout << "                         Load charges: " << t1-t0 << std::endl;
      // Apply U to q -> O(N)
      t0 = chronos::time();
      int nleaves = (1<<lvl)*(1<<lvl)*(1<<lvl);
      int Lu3 = Lu*Lu*Lu;
      double *Uq = new double[nleaves*Lu3];
      int m_size = n_terms<charge_type>();
      for(int ii = 0; ii < nleaves; ii++){
	if(grid[ii].nprt != 0){
	  gemm(1.,U[ii],&(grid[ii].Qs)[0],0.,Uq+ii*Lu3,Lu3,grid[ii].nprt*m_size,1);}
	else{
	  for(int i = 0; i < Lu3; i++){(Uq+ii*Lu3)[i] = 0.;}}
      }
      t1 = chronos::time();
      std::cout << "                         Apply interp: " << t1-t0 << std::endl;
      // Permute entries -> O(N)
      t0 = chronos::time();
      double *tmpp = new double[nleaves*Lu3];
      for(int p0 = 0; p0 < (1<<lvl); p0++){
	for(int p1 = 0; p1 < (1 << lvl); p1++){
	  for(int p2 = 0; p2 < (1 << lvl); p2++){
	    for(int i0 = 0; i0 < Lu; i0++){
	      for(int i1 = 0; i1 < Lu; i1++){
		for(int i2 = 0; i2 < Lu; i2++){
		  tmpp[((p0*Lu+i0)*(1<<lvl)*Lu + (p1*Lu+i1))*(1<<lvl)*Lu + (p2*Lu+i2)]
		    = Uq[(p0*(1<<lvl)*(1<<lvl)+p1*(1<<lvl)+p2)*Lu3 + i0*Lu*Lu+i1*Lu+i2];
		}
	      }
	    }
	  }
	}
      }
      t1 = chronos::time();
      std::cout << "                         Pemrute entries: " << t1-t0 << std::endl;
      // Compute long range interactions -> O(NlogN)
      t0 = chronos::time();
      double far_field  = rapprox((1<<lvl),Lu,box_rad,lvl,tmpp,diag);
      t1 = chronos::time();
      std::cout << "                         Far field: " << t1-t0 << std::endl;
      // Compute short range interactions -> O(N)
      t0 = chronos::time();
      double near_field = short_range<charge_type,KRNL>(grid,lvl,Nimages,2*box_rad);
      t1 = chronos::time();
      std::cout << "                         Near field: " << t1-t0 << std::endl;
      // Ewald correction
      t0 = chronos::time();
      double ewald_rec_and_self = reciprocal_and_self_influence(parts, q, nprt, Nm, box_rad);
      t1 = chronos::time();
      std::cout << "                         Self + Rec: " << t1-t0 << std::endl;
      // Output sum of contributions -> O(1)
      return (near_field + far_field)/2. + ewald_rec_and_self;
    }

  }; // MATRIX

}//ANKH

#endif
