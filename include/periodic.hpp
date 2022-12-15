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
#ifndef PERIODIC_HPP
#define PERIODIC_HPP

#include "kernel.hpp"
#include <complex>

namespace ankh{
    
  inline double scaled_unif_poly(int i, const double& x, int P, double alpha){
    double res = 1.;
    double pi  = -1.+2.*((double)( i)/(double)(P-1));
    pi *= alpha;
    for(int ii=0; ii<i; ii++){
      double pii = -1.+2.*((double)(ii)/(double)(P-1));
      pii *= alpha;
      res *= (x-pii)/(pi-pii);}
    for(int ii=i+1; ii<P; ii++){
      double pii = -1.+2.*((double)(ii)/(double)(P-1));
      pii *= alpha;
      res *= (x-pii)/(pi-pii);}
    return res;
  }
  
  inline double scaled_unif_node(int k, int L, double alpha){
    return alpha*(-1.+double(2*k)/double(L-1));}

  /*
    Unscaled (w.r.t. original box size) precomputed vector
  */
  template<typename FLT, int KRNL>
  void get_periodic_interp_vector(FLT* vector, int L, int Nimages){
    for(int p0 = 0; p0 < L; p0++){
      FLT x = scaled_unif_node(p0,L,2.);
      for(int p1 = 0; p1 < L; p1++){
	FLT y = scaled_unif_node(p1,L,2.);
	for(int p2 = 0; p2 < L; p2++){
	  int idx = p0*L*L+p1*L+p2;
	  vector[idx] = 0.;
	  FLT z = scaled_unif_node(p2,L,2.);
	  for(int i = -Nimages; i <= Nimages; i++){
	    FLT xx = x + FLT(2*i);
	    for(int j = -Nimages; j <= Nimages; j++){
	      FLT yy = y + FLT(2*j);
	      for(int k = -Nimages; k <= Nimages; k++){
		FLT zz = z + FLT(2*k);
		FLT R  = xx*xx+yy*yy+zz*zz;
		if(((i*i > 1) || (j*j > 1) || (k*k > 1)) && sqrt(double(i*i+j*j+k*k)) <= Nimages){
		  vector[idx] += interact<FLT,KRNL>(R);
		}
	      }
	    }
	  }
	}
      }
    }
  }

  template<typename FLT, int KRNL>
  void sym_get_periodic_interp_vector(FLT* vector, int L, int Nimages){
    for(int p0 = 0; p0 < L; p0++){
      FLT x = scaled_unif_node(p0,L,2.);
      for(int p1 = 0; p1 <= p0; p1++){
	FLT y = scaled_unif_node(p1,L,2.);
	for(int p2 = 0; p2 <= p1; p2++){
	  int idx = p0*L*L+p1*L+p2;
	  vector[idx] = 0.;
	  FLT z = scaled_unif_node(p2,L,2.);
	  for(int i = -Nimages; i <= Nimages; i++){
	    FLT xx = x + FLT(2*i);
	    for(int j = -Nimages; j <= Nimages; j++){
	      FLT yy = y + FLT(2*j);
	      for(int k = -Nimages; k <= Nimages; k++){
		FLT zz = z + FLT(2*k);
		FLT R  = xx*xx+yy*yy+zz*zz;
	        if(((i*i > 1) || (j*j > 1) || (k*k > 1)) && sqrt(double(i*i+j*j+k*k)) <= Nimages){
		  vector[idx] += interact<FLT,KRNL>(R);
		}
	      }
	    }
	  }
	  // (1,2,3)
	  vector[     p0 *L*L +      p1 *L + (L-1-p2)] = vector[idx]; 
	  vector[     p0 *L*L + (L-1-p1)*L +      p2 ] = vector[idx]; 
	  vector[     p0 *L*L + (L-1-p1)*L + (L-1-p2)] = vector[idx];
	  vector[(L-1-p0)*L*L +      p1 *L +      p2 ] = vector[idx]; 
	  vector[(L-1-p0)*L*L +      p1 *L + (L-1-p2)] = vector[idx];
	  vector[(L-1-p0)*L*L + (L-1-p1)*L +      p2 ] = vector[idx]; 
	  vector[(L-1-p0)*L*L + (L-1-p1)*L + (L-1-p2)] = vector[idx];
	  // (1,3,2)
	  vector[     p0 *L*L +      p2 *L +      p1 ] = vector[idx];
	  vector[     p0 *L*L +      p2 *L + (L-1-p1)] = vector[idx]; 
	  vector[     p0 *L*L + (L-1-p2)*L +      p1 ] = vector[idx]; 
	  vector[     p0 *L*L + (L-1-p2)*L + (L-1-p1)] = vector[idx];
	  vector[(L-1-p0)*L*L +      p2 *L +      p1 ] = vector[idx]; 
	  vector[(L-1-p0)*L*L +      p2 *L + (L-1-p1)] = vector[idx];
	  vector[(L-1-p0)*L*L + (L-1-p2)*L +      p1 ] = vector[idx]; 
	  vector[(L-1-p0)*L*L + (L-1-p2)*L + (L-1-p1)] = vector[idx];
	  // (2,1,3)
	  vector[     p1 *L*L +      p0 *L +      p2 ] = vector[idx];
	  vector[     p1 *L*L +      p0 *L + (L-1-p2)] = vector[idx]; 
	  vector[     p1 *L*L + (L-1-p0)*L +      p2 ] = vector[idx]; 
	  vector[     p1 *L*L + (L-1-p0)*L + (L-1-p2)] = vector[idx];
	  vector[(L-1-p1)*L*L +      p0 *L +      p2 ] = vector[idx]; 
	  vector[(L-1-p1)*L*L +      p0 *L + (L-1-p2)] = vector[idx];
	  vector[(L-1-p1)*L*L + (L-1-p0)*L +      p2 ] = vector[idx]; 
	  vector[(L-1-p1)*L*L + (L-1-p0)*L + (L-1-p2)] = vector[idx];
	  // (2,3,1)
	  vector[     p1 *L*L +      p2 *L +      p0 ] = vector[idx];
	  vector[     p1 *L*L +      p2 *L + (L-1-p0)] = vector[idx]; 
	  vector[     p1 *L*L + (L-1-p2)*L +      p0 ] = vector[idx]; 
	  vector[     p1 *L*L + (L-1-p2)*L + (L-1-p0)] = vector[idx];
	  vector[(L-1-p1)*L*L +      p2 *L +      p0 ] = vector[idx]; 
	  vector[(L-1-p1)*L*L +      p2 *L + (L-1-p0)] = vector[idx];
	  vector[(L-1-p1)*L*L + (L-1-p2)*L +      p0 ] = vector[idx]; 
	  vector[(L-1-p1)*L*L + (L-1-p2)*L + (L-1-p0)] = vector[idx];
	  // (3,2,1)
	  vector[     p2 *L*L +      p1 *L +      p0 ] = vector[idx];
	  vector[     p2 *L*L +      p1 *L + (L-1-p0)] = vector[idx]; 
	  vector[     p2 *L*L + (L-1-p1)*L +      p0 ] = vector[idx]; 
	  vector[     p2 *L*L + (L-1-p1)*L + (L-1-p0)] = vector[idx];
	  vector[(L-1-p2)*L*L +      p1 *L +      p0 ] = vector[idx]; 
	  vector[(L-1-p2)*L*L +      p1 *L + (L-1-p0)] = vector[idx];
	  vector[(L-1-p2)*L*L + (L-1-p1)*L +      p0 ] = vector[idx]; 
	  vector[(L-1-p2)*L*L + (L-1-p1)*L + (L-1-p0)] = vector[idx];
	  // (3,1,2)
	  vector[     p2 *L*L +      p0 *L +      p1 ] = vector[idx];
	  vector[     p2 *L*L +      p0 *L + (L-1-p1)] = vector[idx]; 
	  vector[     p2 *L*L + (L-1-p0)*L +      p1 ] = vector[idx]; 
	  vector[     p2 *L*L + (L-1-p0)*L + (L-1-p1)] = vector[idx];
	  vector[(L-1-p2)*L*L +      p0 *L +      p1 ] = vector[idx]; 
	  vector[(L-1-p2)*L*L +      p0 *L + (L-1-p1)] = vector[idx];
	  vector[(L-1-p2)*L*L + (L-1-p0)*L +      p1 ] = vector[idx]; 
	  vector[(L-1-p2)*L*L + (L-1-p0)*L + (L-1-p1)] = vector[idx];
	}
      }
    }
  }

  template<typename FLT, int KRNL>
  void sym_get_periodic_interp_vector(FLT* vector, int L, int Nimages, FLT rad){
    FLT scale = rad;
    for(int p0 = 0; p0 < L; p0++){
      FLT x = scaled_unif_node(p0,L,2.);
      for(int p1 = 0; p1 <= p0; p1++){
	FLT y = scaled_unif_node(p1,L,2.);
	for(int p2 = 0; p2 <= p1; p2++){
	  int idx = p0*L*L+p1*L+p2;
	  vector[idx] = 0.;
	  FLT z = scaled_unif_node(p2,L,2.);
	  for(int i = -Nimages; i <= Nimages; i++){
	    FLT xx = (x + FLT(2*i))*scale;
	    for(int j = -Nimages; j <= Nimages; j++){
	      FLT yy = (y + FLT(2*j))*scale;
	      for(int k = -Nimages; k <= Nimages; k++){
		FLT zz = (z + FLT(2*k))*scale;
		FLT R  = xx*xx+yy*yy+zz*zz;
	        if((i*i > 1) || (j*j > 1) || (k*k > 1)){
		  vector[idx] += interact<FLT,KRNL>(R);
		}
	      }
	    }
	  }
	  // (1,2,3)
	  vector[     p0 *L*L +      p1 *L + (L-1-p2)] = vector[idx]; 
	  vector[     p0 *L*L + (L-1-p1)*L +      p2 ] = vector[idx]; 
	  vector[     p0 *L*L + (L-1-p1)*L + (L-1-p2)] = vector[idx];
	  vector[(L-1-p0)*L*L +      p1 *L +      p2 ] = vector[idx]; 
	  vector[(L-1-p0)*L*L +      p1 *L + (L-1-p2)] = vector[idx];
	  vector[(L-1-p0)*L*L + (L-1-p1)*L +      p2 ] = vector[idx]; 
	  vector[(L-1-p0)*L*L + (L-1-p1)*L + (L-1-p2)] = vector[idx];
	  // (1,3,2)
	  vector[     p0 *L*L +      p2 *L +      p1 ] = vector[idx];
	  vector[     p0 *L*L +      p2 *L + (L-1-p1)] = vector[idx]; 
	  vector[     p0 *L*L + (L-1-p2)*L +      p1 ] = vector[idx]; 
	  vector[     p0 *L*L + (L-1-p2)*L + (L-1-p1)] = vector[idx];
	  vector[(L-1-p0)*L*L +      p2 *L +      p1 ] = vector[idx]; 
	  vector[(L-1-p0)*L*L +      p2 *L + (L-1-p1)] = vector[idx];
	  vector[(L-1-p0)*L*L + (L-1-p2)*L +      p1 ] = vector[idx]; 
	  vector[(L-1-p0)*L*L + (L-1-p2)*L + (L-1-p1)] = vector[idx];
	  // (2,1,3)
	  vector[     p1 *L*L +      p0 *L +      p2 ] = vector[idx];
	  vector[     p1 *L*L +      p0 *L + (L-1-p2)] = vector[idx]; 
	  vector[     p1 *L*L + (L-1-p0)*L +      p2 ] = vector[idx]; 
	  vector[     p1 *L*L + (L-1-p0)*L + (L-1-p2)] = vector[idx];
	  vector[(L-1-p1)*L*L +      p0 *L +      p2 ] = vector[idx]; 
	  vector[(L-1-p1)*L*L +      p0 *L + (L-1-p2)] = vector[idx];
	  vector[(L-1-p1)*L*L + (L-1-p0)*L +      p2 ] = vector[idx]; 
	  vector[(L-1-p1)*L*L + (L-1-p0)*L + (L-1-p2)] = vector[idx];
	  // (2,3,1)
	  vector[     p1 *L*L +      p2 *L +      p0 ] = vector[idx];
	  vector[     p1 *L*L +      p2 *L + (L-1-p0)] = vector[idx]; 
	  vector[     p1 *L*L + (L-1-p2)*L +      p0 ] = vector[idx]; 
	  vector[     p1 *L*L + (L-1-p2)*L + (L-1-p0)] = vector[idx];
	  vector[(L-1-p1)*L*L +      p2 *L +      p0 ] = vector[idx]; 
	  vector[(L-1-p1)*L*L +      p2 *L + (L-1-p0)] = vector[idx];
	  vector[(L-1-p1)*L*L + (L-1-p2)*L +      p0 ] = vector[idx]; 
	  vector[(L-1-p1)*L*L + (L-1-p2)*L + (L-1-p0)] = vector[idx];
	  // (3,2,1)
	  vector[     p2 *L*L +      p1 *L +      p0 ] = vector[idx];
	  vector[     p2 *L*L +      p1 *L + (L-1-p0)] = vector[idx]; 
	  vector[     p2 *L*L + (L-1-p1)*L +      p0 ] = vector[idx]; 
	  vector[     p2 *L*L + (L-1-p1)*L + (L-1-p0)] = vector[idx];
	  vector[(L-1-p2)*L*L +      p1 *L +      p0 ] = vector[idx]; 
	  vector[(L-1-p2)*L*L +      p1 *L + (L-1-p0)] = vector[idx];
	  vector[(L-1-p2)*L*L + (L-1-p1)*L +      p0 ] = vector[idx]; 
	  vector[(L-1-p2)*L*L + (L-1-p1)*L + (L-1-p0)] = vector[idx];
	  // (3,1,2)
	  vector[     p2 *L*L +      p0 *L +      p1 ] = vector[idx];
	  vector[     p2 *L*L +      p0 *L + (L-1-p1)] = vector[idx]; 
	  vector[     p2 *L*L + (L-1-p0)*L +      p1 ] = vector[idx]; 
	  vector[     p2 *L*L + (L-1-p0)*L + (L-1-p1)] = vector[idx];
	  vector[(L-1-p2)*L*L +      p0 *L +      p1 ] = vector[idx]; 
	  vector[(L-1-p2)*L*L +      p0 *L + (L-1-p1)] = vector[idx];
	  vector[(L-1-p2)*L*L + (L-1-p0)*L +      p1 ] = vector[idx]; 
	  vector[(L-1-p2)*L*L + (L-1-p0)*L + (L-1-p1)] = vector[idx];
	}
      }
    }
  }

  inline double evaluate_periodic_influence(double x, double y, double z, 
					    double* far_influence, int L, double& rad){
    double scale = 2.*rad;
    double res = 0.;
    for(int j = 0; j < L; j++){
      double sx = scaled_unif_poly(j,x,L,scale);
      for(int k = 0; k < L; k++){
	double sy = scaled_unif_poly(k,y,L,scale);
	for(int l = 0; l < L; l++){
	  double sz = scaled_unif_poly(l,z,L,scale);
	  int idx = j*L*L+k*L+l;
	  res += sx*sy*sz*far_influence[idx];
        }
      }
    }
    return res;
  }

  template<typename charge_type>
  inline double reciprocal_and_self_influence(double *parts, charge_type* q, int N, int Nm, double& rad){
    // Reciprocal
    double L    = 2.*rad * 2.*rad * 2.*rad;
    double bpl2 = M_PI*M_PI/(kernel::beta*kernel::beta);
    double ek = 0.;
    for(int i = 0; i <= Nm; i++){
      for(int j = 0; j <= Nm; j++){
	for(int k = 0; k <= Nm; k++){
	  std::complex<double> Sm = std::complex<double>(0.,0.);
	  for(int p = 0; p < N; p++){
	    double mrj = 0.;
	    mrj += double(i) * parts[p*3  ];
	    mrj += double(j) * parts[p*3+1];
	    mrj += double(k) * parts[p*3+2];
	    Sm += exp(std::complex<double>(0.,M_PI*mrj/rad))*q[j];
	  }
	  double norn2 = double(i*i + j*j + k*k) / (4.*rad*rad);
	  if(i!=0 || j!=0 || k!=0){
	    ek += 2.*exp(-bpl2*double(norn2)) / double(norn2) * std::norm(Sm);
	  }
	}
      }
    }
    ek *= 0.5/M_PI/L;
    // Self
    double e0 = 0.;
    for(int i = 0; i < N; i++){
      e0 += q[i]*q[i];
    }
    e0 *= -kernel::beta/sqrt(M_PI);
    // Result
    return ek + e0;
  }
  
  template<typename charge_type>
  inline double self_influence(charge_type* q, int N, double& rad){
    // Self
    double e0 = 0.;
    for(int i = 0; i < N; i++){
      e0 += q[i]*q[i];
    }
    e0 *= -kernel::beta/sqrt(M_PI);
    // Result
    return e0;
  }

  template<>
  inline double reciprocal_and_self_influence<charge_dipole_quadrupole_3<double> >(double *parts, charge_dipole_quadrupole_3<double>* q, int N, int Nm, double& rad){
    // Reciprocal
    double               L    = 2.*rad * 2.*rad * 2.*rad;
    double               bpl2 = M_PI*M_PI/(kernel::beta*kernel::beta);
    double               ek   = 0.;
    std::complex<double> cste = std::complex<double>(0.,2.*M_PI);
    double               cst2 = 4.*M_PI*M_PI;
    for(int i = 0; i <= Nm; i++){
      for(int j = 0; j <= Nm; j++){
	for(int k = 0; k <= Nm; k++){
	  std::complex<double> Sm = std::complex<double>(0.,0.);
	  for(int p = 0; p < N; p++){
	    double mrj = 0.;
	    mrj += double(i) * parts[p*3  ];
	    mrj += double(j) * parts[p*3+1];
	    mrj += double(k) * parts[p*3+2];
	    std::complex<double> qtot = q[p].q;
	    qtot += cste * (q[p].mu0*double(i)+q[p].mu1*double(j)+q[p].mu2*double(k));
	    double tmp = i * (i*q[p].theta00 
			      + j*q[p].theta01
			      + k*q[p].theta02);
	    tmp  += j * (i*q[p].theta01
			 + j*q[p].theta11
			 + k*q[p].theta12);
	    tmp  += k * (i*q[p].theta02
			 + j*q[p].theta12
			 + k*q[p].theta22);
	    qtot -= tmp;
	    Sm += exp(std::complex<double>(0.,M_PI*mrj/rad))*qtot;
	  }
	  double norn2 = double(i*i + j*j + k*k) / (4.*rad*rad);
	  if(i!=0 || j!=0 || k!=0){
	    ek += 2.*exp(-bpl2*double(norn2)) / double(norn2) * std::norm(Sm);
	  }
	}
      }
    }
    ek *= 0.5/M_PI/L;
    std::cout << "Reciproque : " << ek << std::endl;
    // Self
    double e0 = 0.;
    double e1 = 0.;
    double e2 = 0.;
    for(int i = 0; i < N; i++){
      e0 +=      q[i].q       * q[i].q;
      e1 +=      q[i].mu0     * q[i].mu0;
      e1 +=      q[i].mu1     * q[i].mu1;
      e1 +=      q[i].mu2     * q[i].mu2;
      e2 +=      q[i].theta00 * q[i].theta00;
      e2 += 2. * q[i].theta01 * q[i].theta01;
      e2 += 2. * q[i].theta02 * q[i].theta02;
      e2 +=      q[i].theta11 * q[i].theta11;
      e2 += 2. * q[i].theta12 * q[i].theta12;
      e2 +=      q[i].theta22 * q[i].theta22;
    }
    e1 *= 2./3.*kernel::beta*kernel::beta;
    e2 *= 8./5.*kernel::beta*kernel::beta*kernel::beta*kernel::beta;
    std::cout << "Self :      " << kernel::beta/sqrt(M_PI)*(e0 + e1 + e2) << std::endl;
    std::cout << "Self + Rec : " << (ek - kernel::beta/sqrt(M_PI)*(e0 + e1 + e2)) << std::endl;
    // Result
    return (ek - kernel::beta/sqrt(M_PI)*(e0 + e1 + e2));
  }


} // ANKH

#endif
