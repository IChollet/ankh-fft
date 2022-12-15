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
#ifndef ANKH_SHORT_HPP
#define ANKH_SHORT_HPP

#include "grid.hpp"
#include "multipoles.hpp"
#include "kernel.hpp"

namespace ankh{

  template<typename charge_type, int KRNL>
  double direct_interaction(leaf<charge_type>& trg, leaf<charge_type>& src, int Nimages, double diam){
    std::cout << "direct_interaction not implemented for the general case" << std::endl;}

  template<>
  double direct_interaction<double,0>(leaf<double>& trg, leaf<double>& src, int Nimages, double diam){
    double res = 0.;
    for(int i = 0; i < trg.nprt; i++){
      for(int j = 0; j < src.nprt; j++){
	for(int p0 = -Nimages; p0 <= Nimages; p0++){
	  for(int p1 = -Nimages; p1 <= Nimages; p1++){
	    for(int p2 = -Nimages; p2 <= Nimages; p2++){
	      double dx = trg.prts[i].x - src.prts[j].x + double(p0)*diam;
	      double dy = trg.prts[i].y - src.prts[j].y + double(p1)*diam;
	      double dz = trg.prts[i].z - src.prts[j].z + double(p2)*diam;
	      double R  = dx*dx+dy*dy+dz*dz;
	      double K  = 0.;
	      if(R > 1.e-16){K = interact<double,0>(R);}
	      res += trg.Qs[i]*src.Qs[j]*K;
	    }
	  }
	}
      }
    }
    return res;
  }

  template<>
  double direct_interaction<double,1>(leaf<double>& trg, leaf<double>& src, int Nimages, double diam){
    double res = 0.;
    for(int i = 0; i < trg.nprt; i++){
      for(int j = 0; j < src.nprt; j++){
	for(int p0 = -Nimages; p0 <= Nimages; p0++){
	  for(int p1 = -Nimages; p1 <= Nimages; p1++){
	    for(int p2 = -Nimages; p2 <= Nimages; p2++){
	      double dx = trg.prts[i].x - src.prts[j].x + double(p0)*diam;
	      double dy = trg.prts[i].y - src.prts[j].y + double(p1)*diam;
	      double dz = trg.prts[i].z - src.prts[j].z + double(p2)*diam;
	      res += energy<double,1,0>(dx,dy,dz,&(trg.Qs[i]),&(src.Qs[j]));
	    }
	  }
	}
      }
    }
    return res;
  }

  template<>
  double direct_interaction<charge_dipole_quadrupole_3<double> ,1>(leaf<charge_dipole_quadrupole_3<double> >& trg, leaf<charge_dipole_quadrupole_3<double> >& src, int Nimages, double diam){
    double res = 0.;
    for(int i = 0; i < trg.nprt; i++){
      for(int j = 0; j < src.nprt; j++){
	for(int p0 = -Nimages; p0 <= Nimages; p0++){
	  for(int p1 = -Nimages; p1 <= Nimages; p1++){
	    for(int p2 = -Nimages; p2 <= Nimages; p2++){
	      double dx = trg.prts[i].x - src.prts[j].x + double(p0)*diam;
	      double dy = trg.prts[i].y - src.prts[j].y + double(p1)*diam;
	      double dz = trg.prts[i].z - src.prts[j].z + double(p2)*diam;
	      res += energy<double,1,2>(dx,dy,dz,&(trg.Qs[10*i]),&(src.Qs[10*j]));
	    }
	  }
	}
      }
    }
    return res;
  }

  template<>
  double direct_interaction<charge_dipole_3<double>,0>(leaf<charge_dipole_3<double> >& trg, 
						       leaf<charge_dipole_3<double> >& src,
						       int Nimages, double diam){
    double res = 0.;
    for(int i = 0; i < trg.nprt; i++){
      double tq   = trg.Qs[4*i  ];
      double tmu0 = trg.Qs[4*i+1];
      double tmu1 = trg.Qs[4*i+2];
      double tmu2 = trg.Qs[4*i+3];
      for(int j = 0; j < src.nprt; j++){
	double sq   = src.Qs[4*j  ];
	double smu0 = src.Qs[4*j+1];
	double smu1 = src.Qs[4*j+2];
	double smu2 = src.Qs[4*j+3];
	for(int p0 = -Nimages; p0 <= Nimages; p0++){
	  for(int p1 = -Nimages; p1 <= Nimages; p1++){
	    for(int p2 = -Nimages; p2 <= Nimages; p2++){
	      double dx = trg.prts[i].x - src.prts[j].x + double(p0)*diam;
	      double dy = trg.prts[i].y - src.prts[j].y + double(p1)*diam;
	      double dz = trg.prts[i].z - src.prts[j].z + double(p2)*diam;
	      double R  = dx*dx+dy*dy+dz*dz;
	      double K  = 0.;
	      if(R > 1.e-16){K = interact<double,0>(R);}
	      double K3 = K*K*K;
	      double K5 = K3*K*K;
	      res += tq   * sq   * K;                 // (000,000)
	      res += tq   * smu0 * dx*K3;             // (000,100)
	      res += tq   * smu1 * dy*K3;             // (000,010)
	      res += tq   * smu2 * dz*K3;             // (000,001)
	      res -= tmu0 * sq   * dx*K3;             // (100,000)
	      res += tmu0 * smu0 * (K3-3.*dx*dx*K5);  // (100,100)
	      res -= tmu0 * smu1 * 3.*dx*dy*K5;       // (100,010)
	      res -= tmu0 * smu2 * 3.*dx*dz*K5;       // (100,001)
	      res -= tmu1 * sq   * dy*K3;             // (010,000)
	      res -= tmu1 * smu0 * 3.*dy*dx*K5;       // (010,100)
	      res += tmu1 * smu1 * (K3-3.*dy*dy*K5);  // (010,010)
	      res -= tmu1 * smu2 * 3.*dy*dz*K5;       // (010,001)
	      res -= tmu2 * sq   * dz*K3;             // (001,000)
	      res -= tmu2 * smu0 * 3.*dz*dx*K5;       // (001,100)
	      res -= tmu2 * smu1 * 3.*dz*dy*K5;       // (001,010)
	      res += tmu2 * smu2 * (K3-3.*dz*dz*K5);  // (001,001)
	    }
	  }
	}
      }
    }
    return res;
  }
  
  template<typename charge_type, int KRNL>
  double direct_interaction_0(leaf<charge_type>& trg, leaf<charge_type>& src, int u0, int u1, int u2, double diam){
    std::cout << "direct_interaction_0 not implemented for the general case" << std::endl;}

  template<>
  double direct_interaction_0<double,0>(leaf<double>& trg, leaf<double>& src, int u0, int u1, int u2, double diam){
    double res = 0.;
    for(int i = 0; i < trg.nprt; i++){
      for(int j = 0; j < src.nprt; j++){
        double dx = trg.prts[i].x - src.prts[j].x + u0*diam;
	double dy = trg.prts[i].y - src.prts[j].y + u1*diam;
	double dz = trg.prts[i].z - src.prts[j].z + u2*diam;
	double R  = dx*dx+dy*dy+dz*dz;
	double K  = 0.;
	if(R > 1.e-16){K = interact<double,0>(R);}
	res += trg.Qs[i]*src.Qs[j]*K;
      }
    }
    return res;
  }

  template<>
  double direct_interaction_0<double,1>(leaf<double>& trg, leaf<double>& src, int u0, int u1, int u2, double diam){
    double res = 0.;
    for(int i = 0; i < trg.nprt; i++){
      for(int j = 0; j < src.nprt; j++){
        double dx = trg.prts[i].x - src.prts[j].x + u0*diam;
	double dy = trg.prts[i].y - src.prts[j].y + u1*diam;
	double dz = trg.prts[i].z - src.prts[j].z + u2*diam;
	res += energy<double,1,0>(dx,dy,dz,&(trg.Qs[i]),&(src.Qs[j]));
      }
    }
    return res;
  }

  template<>
  double direct_interaction_0<charge_dipole_quadrupole_3<double>,1>(leaf<charge_dipole_quadrupole_3<double> >& trg, leaf<charge_dipole_quadrupole_3<double> >& src, int u0, int u1, int u2, double diam){
    double res = 0.;
    for(int i = 0; i < trg.nprt; i++){
      for(int j = 0; j < src.nprt; j++){
        double dx = trg.prts[i].x - src.prts[j].x + u0*diam;
	double dy = trg.prts[i].y - src.prts[j].y + u1*diam;
	double dz = trg.prts[i].z - src.prts[j].z + u2*diam;
	res += energy<double,1,2>(dx,dy,dz,&(trg.Qs[10*i]),&(src.Qs[10*j]));
      }
    }
    return res;
  }

  template<>
  double direct_interaction_0<charge_dipole_3<double>,0>(leaf<charge_dipole_3<double> >& trg, 
							 leaf<charge_dipole_3<double> >& src,
							 int u0, int u1, int u2, double diam){
    double res = 0.;
    for(int i = 0; i < trg.nprt; i++){
      double tq   = trg.Qs[4*i  ];
      double tmu0 = trg.Qs[4*i+1];
      double tmu1 = trg.Qs[4*i+2];
      double tmu2 = trg.Qs[4*i+3];
      for(int j = 0; j < src.nprt; j++){
	double sq   = src.Qs[4*j  ];
	double smu0 = src.Qs[4*j+1];
	double smu1 = src.Qs[4*j+2];
	double smu2 = src.Qs[4*j+3];
        double dx = trg.prts[i].x - src.prts[j].x + u0*diam;
	double dy = trg.prts[i].y - src.prts[j].y + u1*diam;
	double dz = trg.prts[i].z - src.prts[j].z + u2*diam;
	double R  = dx*dx+dy*dy+dz*dz;
	double K  = 0.;
	if(R > 1.e-16){K = interact<double,0>(R);}
	double K3 = K*K*K;
	double K5 = K3*K*K;
	res += tq   * sq   * K;                 // (000,000)
	res += tq   * smu0 * dx*K3;             // (000,100)
	res += tq   * smu1 * dy*K3;             // (000,010)
	res += tq   * smu2 * dz*K3;             // (000,001)
	res -= tmu0 * sq   * dx*K3;             // (100,000)
	res += tmu0 * smu0 * (K3-3.*dx*dx*K5);  // (100,100)
	res -= tmu0 * smu1 * 3.*dx*dy*K5;       // (100,010)
	res -= tmu0 * smu2 * 3.*dx*dz*K5;       // (100,001)
	res -= tmu1 * sq   * dy*K3;             // (010,000)
	res -= tmu1 * smu0 * 3.*dy*dx*K5;       // (010,100)
	res += tmu1 * smu1 * (K3-3.*dy*dy*K5);  // (010,010)
	res -= tmu1 * smu2 * 3.*dy*dz*K5;       // (010,001)
	res -= tmu2 * sq   * dz*K3;             // (001,000)
	res -= tmu2 * smu0 * 3.*dz*dx*K5;       // (001,100)
	res -= tmu2 * smu1 * 3.*dz*dy*K5;       // (001,010)
	res += tmu2 * smu2 * (K3-3.*dz*dz*K5);  // (001,001)
      }
    }
    return res;
  }

  template<typename charge_type, int KRNL>
  void n_body_direct_interaction_0(leaf<charge_type>& trg, leaf<charge_type>& src, int u0, int u1, int u2, double diam){
    std::cout << "n_body_direct_interaction_0 not implemented for the general case" << std::endl;}

  template<>
  void n_body_direct_interaction_0<double,0>(leaf<double>& trg, leaf<double>& src, int u0, int u1, int u2, double diam){
    for(int i = 0; i < trg.nprt; i++){
      for(int j = 0; j < src.nprt; j++){
        double dx = trg.prts[i].x - src.prts[j].x + u0*diam;
	double dy = trg.prts[i].y - src.prts[j].y + u1*diam;
	double dz = trg.prts[i].z - src.prts[j].z + u2*diam;
	double R  = dx*dx+dy*dy+dz*dz;
	double K  = 0.;
        if(R > 1.e-16){K = interact<double,0>(R);}
	trg.Ps[i] += src.Qs[j]*K;
      }
    }
  }

  template<>
  void n_body_direct_interaction_0<double,1>(leaf<double>& trg, leaf<double>& src, int u0, int u1, int u2, double diam){
    for(int i = 0; i < trg.nprt; i++){
      for(int j = 0; j < src.nprt; j++){
        double dx = trg.prts[i].x - src.prts[j].x + u0*diam;
	double dy = trg.prts[i].y - src.prts[j].y + u1*diam;
	double dz = trg.prts[i].z - src.prts[j].z + u2*diam;
	double R  = dx*dx+dy*dy+dz*dz;
	double K  = 0.;
        if(R > 1.e-16){K = interact<double,1>(R);}
	trg.Ps[i] += src.Qs[j]*K;
      }
    }
  }

  template<typename charge_type, int KRNL>
  double short_range(leaf<charge_type>* leaves, int lvl, int Nimages, double diam){
    int E = (1<<lvl); int nleaves = E*E*E;
    double res = 0.;
    for(int idx_trg_x = 0; idx_trg_x < E; idx_trg_x++){
      for(int idx_trg_y = 0; idx_trg_y < E; idx_trg_y++){
	for(int idx_trg_z = 0; idx_trg_z < E; idx_trg_z++){
	  int idx_trg = idx_in_grid(idx_trg_x,idx_trg_y,idx_trg_z,E);
	  for(int idx_src_x = -1; idx_src_x <= 1; idx_src_x++){
	    for(int idx_src_y = -1; idx_src_y <= 1; idx_src_y++){
	      for(int idx_src_z = -1; idx_src_z <= 1; idx_src_z++){
		int sidx = idx_trg_x+idx_src_x;
		int sidy = idx_trg_y+idx_src_y;
		int sidz = idx_trg_z+idx_src_z;
		int idx_src = idx_in_grid(sidx,sidy,sidz,E);
	        if(Nimages == 0){
		  if(sidx > -1 && sidx < E && sidy > -1 && sidy < E && sidz > -1 && sidz < E){
		    res += direct_interaction_0<charge_type,KRNL>(leaves[idx_trg],leaves[idx_src],0,0,0,diam);
		  }
		}else{
		  int u0 = (sidx < 0 ? E : 0) + (sidx >= E ? -E : 0);
		  int u1 = (sidy < 0 ? E : 0) + (sidy >= E ? -E : 0);
		  int u2 = (sidz < 0 ? E : 0) + (sidz >= E ? -E : 0);
		  idx_src = idx_in_grid(sidx+u0,sidy+u1,sidz+u2,E);
		  u0 = (u0 < 0 ? -1 : 0) + (u0 > 0 ? 1 : 0);
		  u1 = (u1 < 0 ? -1 : 0) + (u1 > 0 ? 1 : 0);
		  u2 = (u2 < 0 ? -1 : 0) + (u2 > 0 ? 1 : 0);
		  res += direct_interaction_0<charge_type,KRNL>(leaves[idx_trg],
								leaves[idx_src],
								u0,u1,u2,diam);
		}
	      }
	    }
	  }
	}
      }
    }
    return res;
  }
  
  template<typename charge_type, int KRNL>
  void n_body_short_range(leaf<charge_type>* leaves, int lvl, int Nimages, double diam){
    int E = (1<<lvl); int nleaves = E*E*E;
    double res = 0.;
    for(int idx_trg_x = 0; idx_trg_x < E; idx_trg_x++){
      for(int idx_trg_y = 0; idx_trg_y < E; idx_trg_y++){
	for(int idx_trg_z = 0; idx_trg_z < E; idx_trg_z++){
	  int idx_trg = idx_in_grid(idx_trg_x,idx_trg_y,idx_trg_z,E);
	  for(int idx_src_x = -1; idx_src_x <= 1; idx_src_x++){
	    for(int idx_src_y = -1; idx_src_y <= 1; idx_src_y++){
	      for(int idx_src_z = -1; idx_src_z <= 1; idx_src_z++){
		int sidx = idx_trg_x+idx_src_x;
		int sidy = idx_trg_y+idx_src_y;
		int sidz = idx_trg_z+idx_src_z;
		int idx_src = idx_in_grid(sidx,sidy,sidz,E);
	        if(Nimages == 0){
		  if(sidx > -1 && sidx < E && sidy > -1 && sidy < E && sidz > -1 && sidz < E){
		    n_body_direct_interaction_0<charge_type,KRNL>(leaves[idx_trg],leaves[idx_src],0,0,0,diam);
		  }
		}else{
		  int u0 = (sidx < 0 ? E : 0) + (sidx >= E ? -E : 0);
		  int u1 = (sidy < 0 ? E : 0) + (sidy >= E ? -E : 0);
		  int u2 = (sidz < 0 ? E : 0) + (sidz >= E ? -E : 0);
		  idx_src = idx_in_grid(sidx+u0,sidy+u1,sidz+u2,E);
		  u0 = (u0 < 0 ? -1 : 0) + (u0 > 0 ? 1 : 0);
		  u1 = (u1 < 0 ? -1 : 0) + (u1 > 0 ? 1 : 0);
		  u2 = (u2 < 0 ? -1 : 0) + (u2 > 0 ? 1 : 0);
		  n_body_direct_interaction_0<charge_type,KRNL>(leaves[idx_trg],
								       leaves[idx_src],
								       u0,u1,u2,diam);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  template<typename charge_type, int KRNL>
  double direct_interaction(leaf<charge_type>* leaves, int lvl, int Nimages, double diam){
    int E = (1<<lvl); int nleaves = E*E*E;
    double res = 0.;
    for(int idx_trg_x = 0; idx_trg_x < E; idx_trg_x++){
      for(int idx_trg_y = 0; idx_trg_y < E; idx_trg_y++){
	for(int idx_trg_z = 0; idx_trg_z < E; idx_trg_z++){
	  int idx_trg = idx_in_grid(idx_trg_x,idx_trg_y,idx_trg_z,E);
	  for(int idx_src_x = 0; idx_src_x < E; idx_src_x++){
	    for(int idx_src_y = 0; idx_src_y < E; idx_src_y++){
	      for(int idx_src_z = 0; idx_src_z < E; idx_src_z++){
		int idx_src = idx_in_grid(idx_src_x,idx_src_y,idx_src_z,E);
	        res += direct_interaction<charge_type,KRNL>(leaves[idx_trg],leaves[idx_src],Nimages,diam);
	      }
	    }
	  }
	}
      }
    }
    return res;
  }

  template<typename charge_type, int KRNL>
  double mutual_direct_interaction(leaf<charge_type>* leaves, int lvl, int Nimages, double diam){
    int E = (1<<lvl); int nleaves = E*E*E;
    double res = 0.;
    for(int idx_trg_x = 0; idx_trg_x < E; idx_trg_x++){
      for(int idx_trg_y = 0; idx_trg_y < E; idx_trg_y++){
	for(int idx_trg_z = 0; idx_trg_z < E; idx_trg_z++){
	  int idx_trg = idx_in_grid(idx_trg_x,idx_trg_y,idx_trg_z,E);
	  for(int idx_src_x = 0; idx_src_x < E; idx_src_x++){
	    for(int idx_src_y = 0; idx_src_y < E; idx_src_y++){
	      for(int idx_src_z = 0; idx_src_z < E; idx_src_z++){
		int idx_src = idx_in_grid(idx_src_x,idx_src_y,idx_src_z,E);
		if(idx_trg == idx_src){
		  res += direct_interaction<charge_type,KRNL>(leaves[idx_trg],leaves[idx_src],Nimages,diam);
		}else{
		  if(idx_trg < idx_src)
		    res += 2.*direct_interaction<charge_type,KRNL>(leaves[idx_trg],leaves[idx_src],Nimages,diam);
		}
	      }
	    }
	  }
	}
      }
    }
    return res;
  }
}

#endif
