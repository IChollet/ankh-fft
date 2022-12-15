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
#ifndef ANKH_KERNEL_HPP
#define ANKH_KERNEL_HPP

#include <iostream>
#include <cmath>

namespace ankh{

  /*
    Kernel definitions
  */
  namespace kernel{
    double beta;
  } // KERNEL
  
  template<typename FLT>
  class Kernel_Coulomb{
  public:
    static inline FLT interact(const FLT& R){return 1. / sqrt(R);}
  }; // Kernel_Coulomb
  
  template<typename FLT>
  class Kernel_Erfc{
  public:
    static inline FLT interact(const FLT& R){return std::erfc(kernel::beta*sqrt(R)) / sqrt(R);}
  }; // Kernel_Erfc

  template<typename FLT, int Kernel_type>
  inline FLT interact(const FLT& R){std::cout << "Abstract interaction is not defined" << std::endl; exit(1); return 0;}
  template<>
  inline double interact<double,0>(const double& R){return 1. / sqrt(R);}
  template<>
  inline double interact<double,1>(const double& R){
    double r = sqrt(R);
    return std::erfc(kernel::beta*r) / r;}

  /*
    Energy computation
       FLT  corresponds to the floating point precision
       KRNL is an identifier for the kernel type
       DER  is an identifier for the number of derivatives to handle
  */
  template<typename FLT, int KRNL, int DER>
  inline FLT energy(const FLT& dx, const FLT& dy, const FLT& dz, FLT* trg,FLT* src){
    std::cout << "Error: energy underfined for arbitrary template." << std::endl; exit(1);}
  
  // KRNL = 0
  template<> inline double energy<double,0,0>(const double& dx, const double& dy, const double& dz, double* trg, double* src){
    double R  = dx*dx+dy*dy+dz*dz;
    double K  = 0.;
    if(R > 1.e-16){K = interact<double,0>(R);}		
    return trg[0]*src[0]*K;
  }
  template<> inline double energy<double,0,1>(const double& dx, const double& dy, const double& dz, double* trg, double* src){
    double R  = dx*dx+dy*dy+dz*dz;
    double K  = 0.;
    if(R > 1.e-16){K = interact<double,0>(R);}		
    double tq   = trg[0];
    double tmu0 = trg[1];
    double tmu1 = trg[2];
    double tmu2 = trg[3];
    double sq   = src[0];
    double smu0 = src[1];
    double smu1 = src[2];
    double smu2 = src[3];
    double K3   = K*K*K;
    double res  = 0.;
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
    return res;
  }
  
  // KRNL = 1
  template<> inline double energy<double,1,0>(const double& dx, const double& dy, const double& dz, double* trg, double* src){
    double R  = dx*dx+dy*dy+dz*dz;
    double K  = 0.;
    if(R > 1.e-16){K = interact<double,1>(R);}		
    return trg[0]*src[0]*K;
  }
  
  template<> inline double energy<double,1,2>(const double& dx, const double& dy, const double& dz, double* trg, double* src){
    double R     = dx*dx+dy*dy+dz*dz;
    double K     = 0.;
    if(R > 1.e-16){K = interact<double,0>(R);}		
    double K2    = K*K;
    double ci    = trg[0];
    double dix   = trg[1];
    double diy   = trg[2];
    double diz   = trg[3];
    double qixx  = trg[4];
    double qixy  = trg[5];
    double qixz  = trg[6];
    double qiyy  = trg[7];
    double qiyz  = trg[8];
    double qizz  = trg[9];
    double ck    = src[0];
    double dkx   = src[1];
    double dky   = src[2];
    double dkz   = src[3];
    double qkxx  = src[4];
    double qkxy  = src[5];
    double qkxz  = src[6];
    double qkyy  = src[7];
    double qkyz  = src[8];
    double qkzz  = src[9];
    double dri   = - dix*dx - diy*dy - diz*dz;     // - <di ,r  >
    double drk   = - dkx*dx - dky*dy - dkz*dz;     // - <dk ,r  >
    double dik   =   dix*dkx + diy*dky + diz*dkz;  //   <di ,dk >
    double qrix  = - qixx*dx - qixy*dy - qixz*dz;  // - <qix,r  > = qrix
    double qriy  = - qixy*dx - qiyy*dy - qiyz*dz;  // - <qiy,r  > = qriy
    double qriz  = - qixz*dx - qiyz*dy - qizz*dz;  // - <qiz,r  > = qriz
    double qrkx  = - qkxx*dx - qkxy*dy - qkxz*dz;  // - <qkx,r  > = qrkx
    double qrky  = - qkxy*dx - qkyy*dy - qkyz*dz;  // - <qky,r  > = qrky
    double qrkz  = - qkxz*dx - qkyz*dy - qkzz*dz;  // - <qkz,r  > = qrkz
    double qrri  = - qrix*dx - qriy*dy - qriz*dz;  // - <qri,r  >
    double qrrk  = - qrkx*dx - qrky*dy - qrkz*dz;  // - <qrk,r  >
    double qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz; // <qri,qrk>
    double qik   = 2.0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz) + qixx*qkxx + qiyy*qkyy + qizz*qkzz;
    double diqrk = dix*qrkx + diy*qrky + diz*qrkz;
    double dkqri = dkx*qrix + dky*qriy + dkz*qriz;
    double term1 = ci*ck;
    double term2 = ck*dri - ci*drk + dik;
    double term3 = ci*qrrk + ck*qrri - dri*drk + 2.0*(dkqri-diqrk+qik);
    double term4 = dri*qrrk - drk*qrri - 4.0*qrrik;
    double term5 = qrri*qrrk;
    double bl    = 0.;
    if(R > 1.e-16){bl = interact<double,1>(R);}
    double res   = bl * term1;
    double p2a   = exp(-R*kernel::beta*kernel::beta)/(kernel::beta*sqrt(M_PI));
    double bt22  = 2.*kernel::beta*kernel::beta;
    p2a *= bt22;
    bl   = K2 * (   bl + p2a);
    res += bl * term2;
    p2a *= bt22;
    bl   = K2 * (3.*bl + p2a);
    res += bl * term3;
    p2a *= bt22;
    bl   = K2 * (5.*bl + p2a);
    res += bl * term4;
    p2a *= bt22;
    bl   = K2 * (7.*bl + p2a);
    res += bl * term5;
    return res;
  }

} // ANKH

#endif
