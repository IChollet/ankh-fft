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
#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include <iostream>

#ifndef BLAS_PARSER_CONSTS_DOUBLE
#define BLAS_PARSER_CONSTS_DOUBLE
#define BLAS_D double
  BLAS_D D_ZERO     =  0.0;
  BLAS_D D_ONE      =  1.0;
  const char* charN = "N" ;
  const char* charT = "T" ;
#endif

  extern "C"{
    void dgemm_(const char*, const char*, unsigned*, unsigned*, unsigned*, BLAS_D*, BLAS_D*, 
		unsigned*, BLAS_D*, unsigned*, BLAS_D*, BLAS_D*, unsigned*);
  }
inline void gemm(BLAS_D u, BLAS_D* A, BLAS_D* B, BLAS_D v, BLAS_D* C, int _n, int _k, int _s){
    unsigned int n = _n;  unsigned int k = _k;  unsigned int s = _s;
    dgemm_(charN,charN,&n,&s,&k,&u,A,&n,B,&k, &v,C,&n);}
inline void gemTm(BLAS_D u, BLAS_D* A, BLAS_D* B, BLAS_D v, BLAS_D* C, int _n, int _k, int _s){
  unsigned int n = _n;  unsigned int k = _k;  unsigned int s = _s;
  dgemm_(charT,charN,&n,&s,&k,&u,A,&k,B,&k, &v,C,&n);}


namespace ankh{
  
  inline double S1D(int i, const double& x, int P){
    double res = 1.;
    double pi  = -1.+2.*((double)( i)/(double)(P-1));
    for(int ii=0; ii<i; ii++){
      double pii = -1.+2.*((double)(ii)/(double)(P-1));
      res *= (x-pii)/(pi-pii);}
    for(int ii=i+1; ii<P; ii++){
      double pii = -1.+2.*((double)(ii)/(double)(P-1));
      res *= (x-pii)/(pi-pii);}
    return res;
  }
  
  inline double unif_node(int k, int L){return -1.+double(2*k)/double(L-1);}

  inline double cheb_node(int k, int L){return std::cos((double)(2*k+1)/(double)(2*L)*M_PI);}
  inline double T(double x, int k){return cos((double)(k)*acos(x));}
  inline double U(double x, int k){if(k < 0){return 0.;}return sin(double(k+1)*acos(x))/sqrt(1. - x*x);}
  inline double dT(double x, int k){
    double res = (double)(((k%2) == 1 ? -k : 0));
    for(int j = 0; j <= floor((double)(k-1)*0.5); j++){
      res += 2.*(double)(k)*T(x,k-1-2*j);}
    return res;}
  inline double ddT(double x, int k){
    double res = 0.;
    for(int j = 0; j <= floor((double)(k-1)*0.5); j++){
      res += 2.*(double)(k)*dT(x,k-1-2*j);}
    return res;}
  inline void getTr(int L, double *Tr){
    for(int k = 0; k < L; k++){
    double r = cheb_node(k,L); for(int i = 0; i < L; i++){(Tr+k*L)[i] = T(r,i);}}}
  inline void getTs  (double u, int L, double *  Tu){
    for(int i = 0; i < L; i++){  Tu[i] =   T(u,i);}}
  inline void getdTs (double u, int L, double * dTu){
    for(int i = 0; i < L; i++){ dTu[i] =  dT(u,i);}}
  inline void getddTs(double u, int L, double *ddTu){
    for(int i = 0; i < L; i++){ddTu[i] = ddT(u,i);}}
  inline void   C1D(double *  Tx, double *Tr, int L, double *  C){
    for(int k = 0; k < L; k++){
      double res = 1.;
      double *Trk = Tr+k*L;
      for(int j = 1; j < L; j++){
	res += 2. * Tx[j] * Trk[j];}
      C[k] = res/(double)(L);}
  }

  inline void   loop_C1D(double x, int L, double *  C){
    for(int k = 0; k < L; k++){
      double res = 1.;
      for(int j = 0; j < L; j++){
	if(j!=k){res *= (x-cheb_node(j,L))/(cheb_node(k,L)-cheb_node(j,L));}
      }
      C[k] = res;
    }
  }

  inline void  dC1D(double * dTx, double *Tr, int L, double * dC){
    for(int k = 0; k < L; k++){
      double res = 0.;
      double *Trk = Tr+k*L;
      for(int j = 1; j < L; j++){
	res += dTx[j] * Trk[j];}
      dC[k] = 2.*res/(double)(L);}
  }
  inline void  loop_dC1D(double x, int L, double *  C){
    for(int k = 0; k < L; k++){
      double res = 0.;
      for(int l = 0; l < L; l++){
	if(l!=k){
	  double tmp = 1./(cheb_node(k,L)-cheb_node(l,L));
	  for(int j = 0; j < L; j++){
	    if(j!=k && j!=l){tmp *= (x-cheb_node(j,L))/(cheb_node(k,L)-cheb_node(j,L));}
	  }
	  res += tmp;
	}
      }
      C[k] = res;
    }
  }

  inline void ddC1D(double * dTx, double *Tr, int L, double * dC){
    for(int k = 0; k < L; k++){
      double res = 0.;
      double *Trk = Tr+k*L;
      for(int j = 1; j < L; j++){
	res += dTx[j] * Trk[j];}
      dC[k] = 2.*res/(double)(L);}
  }

  inline void loop_ddC1D(double x, int L, double *  C){
    for(int k = 0; k < L; k++){
      double res = 0.;
      for(int p = 0; p < L; p++){
	double ttt = 0.;
	if(p!=k){
	  for(int l = 0; l < L; l++){
	    if(l!=k && l!=p){
	      double tmp = 1./(cheb_node(k,L)-cheb_node(l,L));
	      for(int j = 0; j < L; j++){
		if(j!=k && j!=l && j!=p){tmp *= (x-cheb_node(j,L))/(cheb_node(k,L)-cheb_node(j,L));}
	      }
	      ttt += tmp;
	    }
	  }
	}
        res += ttt/(cheb_node(k,L)-cheb_node(p,L));
      }
      C[k] = res;
    }
  }
  
  inline void part_to_unif(double x, double y, double z, double* Mu, int Lu){
    double *Sx = new double[Lu];
    double *Sy = new double[Lu];
    double *Sz = new double[Lu];
    for(int i = 0; i < Lu; i++){
      Sx[i] = S1D(i,x,Lu);
      Sy[i] = S1D(i,y,Lu);
      Sz[i] = S1D(i,z,Lu);
    }
    for(int i = 0; i < Lu; i++){
      for(int j = 0; j < Lu; j++){
	for(int k = 0; k < Lu; k++){
	  Mu[i*Lu*Lu + j*Lu + k] = Sx[i] * Sy[j] * Sz[k];
	}
      }
    }
  }

  inline void cheb_to_unif_1d_matrix(int Lu, int Lc, double*& toUnif){
    toUnif = new double[Lu*Lc];
    for(int j = 0; j < Lc; j++){
      double u = cheb_node(j,Lc);
      for(int i = 0; i < Lu; i++){
	toUnif[j*Lu + i] = S1D(i,u,Lu);
      }
    }
  }

  template<typename charge_type>
  inline void part_to_cheb_to_unif(double x, double y, double z, 
				   double* Mu, int Lu, int Lc,
				   double *toUnif, double *Tr, double grid_leaf_rad){
    std::cout << "error: 'part_to_cheb_to_unif' undefined for arbitrary template" << std::endl;
    exit(1);
  }
  
  template<>
  inline void part_to_cheb_to_unif<double>(double x, double y, double z, 
					   double* Mu, int Lu, int Lc,
					   double *toUnif, double *Tr, double grid_leaf_rad){
    double nx = x / grid_leaf_rad;
    double ny = y / grid_leaf_rad;
    double nz = z / grid_leaf_rad;
    double *Txyz = new double[3*Lc];
    double *Cxyz = new double[3*Lc];
    double *Sxyz = new double[3*Lu];
    getTs(nx, Lc, Txyz     ); C1D(Txyz     , Tr, Lc, Cxyz     );
    getTs(ny, Lc, Txyz+  Lc); C1D(Txyz+  Lc, Tr, Lc, Cxyz+  Lc);
    getTs(nz, Lc, Txyz+2*Lc); C1D(Txyz+2*Lc, Tr, Lc, Cxyz+2*Lc);
    gemm(1.,toUnif,Cxyz,0.,Sxyz,Lu,Lc,3);
    double *Sx = Sxyz, *Sy = Sxyz+Lu, *Sz = Sxyz+2*Lu;
    for(int i = 0; i < Lu; i++){
      for(int j = 0; j < Lu; j++){
	for(int k = 0; k < Lu; k++){
	  Mu[i*Lu*Lu + j*Lu + k] = Sx[i] * Sy[j] * Sz[k];
	}
      }
    }
  }

  template<>
  inline void part_to_cheb_to_unif<charge_dipole_3<double> >(double x, double y, double z, 
							     double* Mu, int Lu, int Lc,
							     double *toUnif, double *Tr,
							     double grid_leaf_rad){
    double nx = x / grid_leaf_rad;
    double ny = y / grid_leaf_rad;
    double nz = z / grid_leaf_rad;
    double *  Txyz = new double[6*Lc];
    double *  Cxyz = new double[6*Lc];
    double *  Sxyz = new double[6*Lu];
    getTs (nx, Lc, Txyz     );  C1D(Txyz     , Tr, Lc, Cxyz     );
    getTs (ny, Lc, Txyz+  Lc);  C1D(Txyz+  Lc, Tr, Lc, Cxyz+  Lc);
    getTs (nz, Lc, Txyz+2*Lc);  C1D(Txyz+2*Lc, Tr, Lc, Cxyz+2*Lc);
    getdTs(nx, Lc, Txyz+3*Lc); dC1D(Txyz+3*Lc, Tr, Lc, Cxyz+3*Lc);
    getdTs(ny, Lc, Txyz+4*Lc); dC1D(Txyz+4*Lc, Tr, Lc, Cxyz+4*Lc);
    getdTs(nz, Lc, Txyz+5*Lc); dC1D(Txyz+5*Lc, Tr, Lc, Cxyz+5*Lc);
    gemm(1.,toUnif,Cxyz,0.,Sxyz,Lu,Lc,6);
    double * Sx = Sxyz     , * Sy = Sxyz+  Lu, * Sz = Sxyz+2*Lu;
    double *dSx = Sxyz+3*Lu, *dSy = Sxyz+4*Lu, *dSz = Sxyz+5*Lu;
    double *ptr = Mu;
    for(int i = 0; i < Lu; i++){
      for(int j = 0; j < Lu; j++){
	for(int k = 0; k < Lu; k++){
	  ptr[i*Lu*Lu + j*Lu + k] =  Sx[i] *  Sy[j] *  Sz[k];
	}
      }
    }
    ptr = Mu+Lu*Lu*Lu;
    for(int i = 0; i < Lu; i++){
      for(int j = 0; j < Lu; j++){
	for(int k = 0; k < Lu; k++){
	  ptr[i*Lu*Lu + j*Lu + k] = dSx[i] *  Sy[j] *  Sz[k] / grid_leaf_rad;
	}
      }
    }
    ptr = Mu+2*Lu*Lu*Lu;
    for(int i = 0; i < Lu; i++){
      for(int j = 0; j < Lu; j++){
	for(int k = 0; k < Lu; k++){
	  ptr[i*Lu*Lu + j*Lu + k] =  Sx[i] * dSy[j] *  Sz[k] / grid_leaf_rad;
	}
      }
    }
    ptr = Mu+3*Lu*Lu*Lu;
    for(int i = 0; i < Lu; i++){
      for(int j = 0; j < Lu; j++){
	for(int k = 0; k < Lu; k++){
	  ptr[i*Lu*Lu + j*Lu + k] =  Sx[i] *  Sy[j] * dSz[k] / grid_leaf_rad;
	}
      }
    }
  }

  template<>
  inline void part_to_cheb_to_unif<charge_dipole_quadrupole_3<double> >(double x, double y, double z, 
									double* Mu, int Lu, int Lc,
									double *toUnif, double *Tr,
									double grid_leaf_rad){
    double nx = x / grid_leaf_rad;
    double ny = y / grid_leaf_rad;
    double nz = z / grid_leaf_rad;
    double grid_leaf_rad_2 = grid_leaf_rad * grid_leaf_rad;
    double *  Txyz = new double[9*Lc];
    double *  Cxyz = new double[9*Lc];
    double *  Sxyz = new double[9*Lu];
    getTs  (nx, Lc, Txyz     );  C1D(Txyz     , Tr, Lc, Cxyz     );
    getTs  (ny, Lc, Txyz+  Lc);  C1D(Txyz+  Lc, Tr, Lc, Cxyz+  Lc);
    getTs  (nz, Lc, Txyz+2*Lc);  C1D(Txyz+2*Lc, Tr, Lc, Cxyz+2*Lc);
    getdTs (nx, Lc, Txyz+3*Lc); dC1D(Txyz+3*Lc, Tr, Lc, Cxyz+3*Lc);
    getdTs (ny, Lc, Txyz+4*Lc); dC1D(Txyz+4*Lc, Tr, Lc, Cxyz+4*Lc);
    getdTs (nz, Lc, Txyz+5*Lc); dC1D(Txyz+5*Lc, Tr, Lc, Cxyz+5*Lc);
    getddTs(nx, Lc, Txyz+6*Lc); ddC1D(Txyz+6*Lc, Tr, Lc, Cxyz+6*Lc);
    getddTs(ny, Lc, Txyz+7*Lc); ddC1D(Txyz+7*Lc, Tr, Lc, Cxyz+7*Lc);
    getddTs(nz, Lc, Txyz+8*Lc); ddC1D(Txyz+8*Lc, Tr, Lc, Cxyz+8*Lc);
    gemm(1.,toUnif,Cxyz,0.,Sxyz,Lu,Lc,9);
    double *  Sx = Sxyz     , *  Sy = Sxyz+  Lu, *  Sz = Sxyz+2*Lu;
    double * dSx = Sxyz+3*Lu, * dSy = Sxyz+4*Lu, * dSz = Sxyz+5*Lu;
    double *ddSx = Sxyz+6*Lu, *ddSy = Sxyz+7*Lu, *ddSz = Sxyz+8*Lu;
    double *ptr = Mu;
    for(int i = 0; i < Lu; i++){
      for(int j = 0; j < Lu; j++){
	for(int k = 0; k < Lu; k++){
	  ptr[i*Lu*Lu + j*Lu + k] =   Sx[i] *   Sy[j] *   Sz[k];
	}
      }
    }
    ptr = Mu+Lu*Lu*Lu;
    for(int i = 0; i < Lu; i++){
      for(int j = 0; j < Lu; j++){
	for(int k = 0; k < Lu; k++){
	  ptr[i*Lu*Lu + j*Lu + k] =  dSx[i] *   Sy[j] *   Sz[k] / grid_leaf_rad;
	}
      }
    }
    ptr = Mu+2*Lu*Lu*Lu;
    for(int i = 0; i < Lu; i++){
      for(int j = 0; j < Lu; j++){
	for(int k = 0; k < Lu; k++){
	  ptr[i*Lu*Lu + j*Lu + k] =   Sx[i] *  dSy[j] *   Sz[k] / grid_leaf_rad;
	}
      }
    }
    ptr = Mu+3*Lu*Lu*Lu;
    for(int i = 0; i < Lu; i++){
      for(int j = 0; j < Lu; j++){
	for(int k = 0; k < Lu; k++){
	  ptr[i*Lu*Lu + j*Lu + k] =   Sx[i] *   Sy[j] *  dSz[k] / grid_leaf_rad;
	}
      }
    }
    ptr = Mu+4*Lu*Lu*Lu;
    for(int i = 0; i < Lu; i++){
      for(int j = 0; j < Lu; j++){
	for(int k = 0; k < Lu; k++){
	  ptr[i*Lu*Lu + j*Lu + k] =  ddSx[i] *  Sy[j] *   Sz[k] / grid_leaf_rad_2;
        }
      }
    }
    ptr = Mu+5*Lu*Lu*Lu;
    for(int i = 0; i < Lu; i++){
      for(int j = 0; j < Lu; j++){
	for(int k = 0; k < Lu; k++){
	  ptr[i*Lu*Lu + j*Lu + k] =   dSx[i] * dSy[j] *   Sz[k] / grid_leaf_rad_2 * 2.;
	}
      }
    }
    ptr = Mu+6*Lu*Lu*Lu;
    for(int i = 0; i < Lu; i++){
      for(int j = 0; j < Lu; j++){
	for(int k = 0; k < Lu; k++){
	  ptr[i*Lu*Lu + j*Lu + k] =   dSx[i] *  Sy[j] *  dSz[k] / grid_leaf_rad_2 * 2.;
	}
      }
    }
    ptr = Mu+7*Lu*Lu*Lu;
    for(int i = 0; i < Lu; i++){
      for(int j = 0; j < Lu; j++){
	for(int k = 0; k < Lu; k++){
	  ptr[i*Lu*Lu + j*Lu + k] =   Sx[i] * ddSy[j] *   Sz[k] / grid_leaf_rad_2;
	}
      }
    }
    ptr = Mu+8*Lu*Lu*Lu;
    for(int i = 0; i < Lu; i++){
      for(int j = 0; j < Lu; j++){
	for(int k = 0; k < Lu; k++){
	  ptr[i*Lu*Lu + j*Lu + k] =   Sx[i] *  dSy[j] *  dSz[k] / grid_leaf_rad_2 * 2.;
	}
      }
    }
    ptr = Mu+7*Lu*Lu*Lu;
    for(int i = 0; i < Lu; i++){
      for(int j = 0; j < Lu; j++){
	for(int k = 0; k < Lu; k++){
	  ptr[i*Lu*Lu + j*Lu + k] =   Sx[i] * ddSy[j] *   Sz[k] / grid_leaf_rad_2;
	}
      }
    }
    ptr = Mu+8*Lu*Lu*Lu;
    for(int i = 0; i < Lu; i++){
      for(int j = 0; j < Lu; j++){
	for(int k = 0; k < Lu; k++){
	  ptr[i*Lu*Lu + j*Lu + k] =   Sx[i] *  dSy[j] *  dSz[k] / grid_leaf_rad_2 * 2.;
	}
      }
    }
    ptr = Mu+9*Lu*Lu*Lu;
    for(int i = 0; i < Lu; i++){
      for(int j = 0; j < Lu; j++){
	for(int k = 0; k < Lu; k++){
	  ptr[i*Lu*Lu + j*Lu + k] =   Sx[i] *   Sy[j] * ddSz[k] / grid_leaf_rad_2;
	}
      }
    }
  }

} // ANKH

#endif
