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
#ifndef ANKH_WRAPPER_HPP
#define ANKH_WRAPPER_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

namespace ankh{

  void read_tinker_file_quad(const char* file_particles,
			     const char* file_dipole_quad,
			     int& Nprt, double*& prts,
			     charge_dipole_quadrupole_3<double>*& charges,
			     double& rad){
    std::string trsh;
    std::ifstream fparts;
    fparts.open(file_particles);
    std::ifstream dipfile;
    dipfile.open(file_dipole_quad);
    fparts >> Nprt;
    getline(fparts,trsh);
    getline(fparts,trsh);
    prts = new double[3*Nprt];
    charges = new charge_dipole_quadrupole_3<double>[Nprt];
    dipfile >> rad;
    rad /= 2.;
    for(int i = 0; i < Nprt; i++){
      fparts >> trsh;
      fparts >> trsh;
      for(int k = 0; k < 3; k++){
	fparts >> prts[3*i+k];
	while(prts[3*i+k] > rad){
	  prts[3*i+k] -= 2.*rad;
	}
        while(prts[3*i+k] < -rad){
	  prts[3*i+k] += 2.*rad;
	}
      }
      dipfile >> charges[i].q;      
      dipfile >> charges[i].mu0;    
      dipfile >> charges[i].mu1;    
      dipfile >> charges[i].mu2;    
      dipfile >> charges[i].theta00;
      dipfile >> charges[i].theta01;
      dipfile >> charges[i].theta02;
      dipfile >> trsh;                
      dipfile >> charges[i].theta11;
      dipfile >> charges[i].theta12;
      dipfile >> trsh;
      dipfile >> trsh;
      dipfile >> charges[i].theta22;

      getline(fparts,trsh);
    }
  }

} // ANKH

#endif
