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
#ifndef ANKH_MULTIPOLES_HPP
#define ANKH_MULTIPOLES_HPP

#include <iostream>

namespace ankh{
  
  template<typename FLT>
  struct charge_dipole_3{
    FLT q  ;
    FLT mu0;
    FLT mu1;
    FLT mu2;
  }; // CHARGE_DIPOLE_3

  template<typename FLT>
  struct charge_dipole_quadrupole_3{
    FLT q  ;
    FLT mu0;
    FLT mu1;
    FLT mu2;
    FLT theta00;
    FLT theta01;
    FLT theta02;
    FLT theta11;
    FLT theta12;
    FLT theta22;
  }; // CHARGE_DIPOLE_QUADRUPOLE_3

  template<typename charge_type>
  inline int n_terms(){
    std::cout << "error: 'n_terms' undefined for arbitrary template" << std::endl;
    exit(1);
  }
  template<>
  inline int n_terms<double>(){return 1;}
  template<>
  inline int n_terms<charge_dipole_3<double> >(){return 4;}
  template<>
  inline int n_terms<charge_dipole_quadrupole_3<double> >(){return 10;}

  template<typename charge_type>
  inline void copy_to_FLT_ptr(charge_type& a, double* b){
    std::cout << "error: 'copy_to_FLT_ptr undefined for arbitrary template'" << std::endl;
    exit(1);
  }
  template<>
  inline void copy_to_FLT_ptr(double& a, double* b){
    (*b) = a;
  }
  template<>
  inline void copy_to_FLT_ptr(charge_dipole_3<double>& a, double* b){
    b[0] = a.q  ;
    b[1] = a.mu0;
    b[2] = a.mu1;
    b[3] = a.mu2;
  }
  template<>
  inline void copy_to_FLT_ptr(charge_dipole_quadrupole_3<double>& a, double* b){
    b[0] = a.q  ;
    b[1] = a.mu0;
    b[2] = a.mu1;
    b[3] = a.mu2;
    b[4] = a.theta00;
    b[5] = a.theta01;
    b[6] = a.theta02;
    b[7] = a.theta11;
    b[8] = a.theta12;
    b[9] = a.theta22;
  }

  template<typename charge_type>
  inline void copy_from_FLT_ptr(charge_type& a, double* b){
    std::cout << "error: 'copy_from_FLT_ptr undefined for arbitrary template'" << std::endl;
    exit(1);
  }
  template<>
  inline void copy_from_FLT_ptr(double& a, double* b){
    a = *b;
  }

} // ANKH

#endif
