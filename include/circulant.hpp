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
#ifndef CIRCULANT_HPP
#define CIRCULANT_HPP

#include <complex>
#include <fftw3.h>
#include "periodic.hpp"
#include "kernel.hpp"

namespace ankh{

  inline bool MAClist(int i, int j, int k){
    int test = 0;
    int ncmp = 3;
    ncmp -= (i ? 0 : 1);
    test += std::abs(i);
    ncmp -= (j ? 0 : 1);
    test += std::abs(j);
    ncmp -= (k ? 0 : 1);
    test += std::abs(k);
    return (test > ncmp);
  }

  inline double psi(int i, int L){return 2.*double(i)/double(L-1);}
  template<int KRNL>
  inline double G(int Lcell, int L,
		  int k0, int l0, int i0, int j0,
		  int k1, int l1, int i1, int j1,
		  int k2, int l2, int i2, int j2,
		  double& u, double& v   ){
    double R0 = psi(k0-l0,Lcell)*u + psi(i0-j0,L)*v;
    double R1 = psi(k1-l1,Lcell)*u + psi(i1-j1,L)*v;
    double R2 = psi(k2-l2,Lcell)*u + psi(i2-j2,L)*v;
    double R  = R0*R0 + R1*R1 + R2*R2;
    return (MAClist(k0-l0,k1-l1,k2-l2) ? interact<double,KRNL>(R) : 0.);
  }

  template<int KRNL>
  void generate_first_col(double *A, int L, int Li, double& rad, int lvl, int Nimages){
    double u   = double((1 << (lvl))-1) * rad / double((1 << (lvl)));
    double v   = rad / double((1 << (lvl)));
    int Lcell  = (1<<lvl);
    int Pcell  = 2*Lcell - 1;
    int P      = 2*L-1;
    int images = (Nimages > 0 ? 1 : 0);
    double *periodic_vector = new double[Li*Li*Li];
    sym_get_periodic_interp_vector<double,KRNL>(periodic_vector,Li,Nimages,rad);
    for(int i0 = 0; i0 < Pcell; i0++){
      for(int i1 = 0; i1 < P; i1++){
	int v0_0, v1_0, w0_0, w1_0;
	if(i0 < Lcell){v0_0 = i0; w0_0 = 0;}else{v0_0 = 0; w0_0 = Pcell - i0;}
	if(i1 < L    ){v1_0 = i1; w1_0 = 0;}else{v1_0 = 0; w1_0 = P     - i1;}
	for(int j0 = 0; j0 < Pcell; j0++){
	  for(int j1 = 0; j1 < P; j1++){
	    int v0_1, v1_1, w0_1, w1_1;
	    if(j0 < Lcell){v0_1 = j0; w0_1 = 0;}else{v0_1 = 0; w0_1 = Pcell - j0;}
	    if(j1 < L    ){v1_1 = j1; w1_1 = 0;}else{v1_1 = 0; w1_1 = P     - j1;}
	    for(int k0 = 0; k0 < Pcell; k0++){
	      for(int k1 = 0; k1 < P; k1++){
		int v0_2, v1_2, w0_2, w1_2;
		if(k0 < Lcell){v0_2 = k0; w0_2 = 0;}else{v0_2 = 0; w0_2 = Pcell - k0;}
		if(k1 < L    ){v1_2 = k1; w1_2 = 0;}else{v1_2 = 0; w1_2 = P     - k1;}
	        A[((i0*P+i1)*Pcell*P + (j0*P+j1))*Pcell*P + (k0*P+k1)] 
		  = evaluate_periodic_influence(psi(v0_0-w0_0,Lcell)*u + psi(v1_0-w1_0,L)*v,
						psi(v0_1-w0_1,Lcell)*u + psi(v1_1-w1_1,L)*v,
						psi(v0_2-w0_2,Lcell)*u + psi(v1_2-w1_2,L)*v,
						periodic_vector,Li,rad);
	        for(int p0 = -images; p0 <= images; p0++){
		  for(int p1 = -images; p1 <= images; p1++){
		    for(int p2 = -images; p2 <= images; p2++){  
		      A[((i0*P+i1)*Pcell*P + (j0*P+j1))*Pcell*P + (k0*P+k1)] 
			+= G<KRNL>(Lcell,L,
				   v0_0,w0_0+p0*Lcell,v1_0,w1_0,
				   v0_1,w0_1+p1*Lcell,v1_1,w1_1,
				   v0_2,w0_2+p2*Lcell,v1_2,w1_2,
				   u,v);
		    }
		  }
		}
	      }}}}}}}

  template<typename T>
  void padding(T *in, T *chi, int L0, int L1){
    int P0 = 2*L0-1;
    int P1 = 2*L1-1;
    for(int i = 0; i < P0*P1*P0*P1*P0*P1; i++){chi[i] = 0.;}
    for(int i0 = 0; i0 < L0; i0++){
      for(int i1 = 0; i1 < L1; i1++){
	for(int j0 = 0; j0 < L0; j0++){
	  for(int j1 = 0; j1 < L1; j1++){
	    for(int k0 = 0; k0 < L0; k0++){
	      for(int k1 = 0; k1 < L1; k1++){
		chi[((i0*P1+i1)*P0*P1 + j0*P1+j1)*P0*P1 + k0*P1+k1]
		  = in[((i0*L1+i1)*L0*L1 + j0*L1+j1)*L0*L1 + k0*L1+k1];
	      }}}}}}}

  template<typename T>
  void reverse_padding(T *out, T *chi, int L0, int L1){
    int P0 = 2*L0-1;
    int P1 = 2*L1-1;
    for(int i0 = 0; i0 < L0; i0++){
      for(int i1 = 0; i1 < L1; i1++){
	for(int j0 = 0; j0 < L0; j0++){
	  for(int j1 = 0; j1 < L1; j1++){
	    for(int k0 = 0; k0 < L0; k0++){
	      for(int k1 = 0; k1 < L1; k1++){
		out[((i0*L1+i1)*L0*L1 + j0*L1+j1)*L0*L1 + k0*L1+k1]
		  = chi[((i0*P1+i1)*P0*P1 + j0*P1+j1)*P0*P1 + k0*P1+k1];
	      }}}}}}}

  template<int KRNL>
  void prcmp_diag(int Lcell,int L,double& rad, int lvl, double*& diag, int Nimages, int Li){
    int Pcell  = 2*Lcell - 1;
    int P      = 2*L     - 1;
    int PP     = Pcell*P*Pcell*P*Pcell*P;
    double *indi = new double[PP];
    diag = new double[2*PP];
    fftw_plan q_forward;
    int DIM = 6; int Parray[DIM];
    Parray[0] = Pcell; Parray[1] = P; Parray[2] = Pcell; Parray[3] = P; Parray[4] = Pcell; Parray[5] = P;
    q_forward  = fftw_plan_dft_r2c(DIM,Parray,indi  ,reinterpret_cast<fftw_complex*>(diag),FFTW_ESTIMATE);
    generate_first_col<KRNL>(indi,L,Li,rad,lvl,Nimages);
    fftw_execute(q_forward);
    fftw_destroy_plan(q_forward);
  }

  double rapprox(int Lcell,int L,double& rad, int lvl, double *in, double *diag){
    double res = 0.;
    int Pcell  = 2*Lcell - 1;
    int P      = 2*L     - 1;
    int PP     = Pcell*P*Pcell*P*Pcell*P;
    double *padded = new double[PP], *out0 = new double[2*PP];
    fftw_plan p_forward, q_forward;
    int DIM = 6; int Parray[DIM];
    Parray[0] = Pcell; Parray[1] = P; Parray[2] = Pcell; Parray[3] = P; Parray[4] = Pcell; Parray[5] = P;
    p_forward  = fftw_plan_dft_r2c(DIM,Parray,padded,reinterpret_cast<fftw_complex*>(out0),FFTW_ESTIMATE);
    padding<double>(in,padded,Lcell,L);
    fftw_execute(p_forward);
    for(int i = 0; i < PP/P; i++){
      std::complex<double> qtmp(out0[2*(i*(P/2+1))],out0[2*(i*(P/2+1))+1]);
      res += norm(qtmp) * diag[2*(i*(P/2+1))] / double(PP);
      for(int j = 1; j < (P/2+1); j++){
	std::complex<double> q(out0[2*(i*(P/2+1)+j)],out0[2*(i*(P/2+1)+j)+1]);
	res += 2. * norm(q) * diag[2*(i*(P/2+1)+j)] / double(PP);
      }
    }
    fftw_destroy_plan(p_forward);
    return res;
  }

}
#endif
