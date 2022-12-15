#include <iomanip> 
#include "../include/ankh.hpp"

int main(int argc, char* argv[]){

  int Lu      = atoi(argv[1]); // Interpolation order at leaves
  int MD      = atoi(argv[2]); // Number of tree levels
  int Nimages = atoi(argv[3]); // Number of images in each direction
  int Li      = atoi(argv[4]); // Interpolation order for far images vectors
  int Lc      = atoi(argv[5]); // Interpolation order of the compressed grid

  int N; // Number of particles
  double *parts;
  ankh::charge_dipole_quadrupole_3<double> *q;
  double rad;double cste;
  ankh::read_tinker_file_quad("../systems/puddle.xyz","../systems/ankh_puddle.txt",N,parts,q,rad);
  std::cout << N << "\t" << rad << std::endl;
  cste = 332.063714000000;
  int Nm = 0;
  
  double t0, t1;
  ankh::matrix<ankh::charge_dipole_quadrupole_3<double>,EWALD> A(Lu,MD,Nimages,Li);
  A.set_kernel_params(0.01); // Set Ewald parameter
  std::cout << "Declaration done" << std::endl;
  t0 = chronos::time();
  A.prcmpt(0.,0.,0.,rad,parts,N,Lc);
  t1 = chronos::time();
  std::cout << "Precomputation done" << std::endl;
  std::cout << t1 - t0 << std::endl;
  t0 = chronos::time();
  double res = A.apply(q,parts,Nm);
  t1 = chronos::time();
  std::cout << "Application done" << std::endl;
  std::cout << t1 - t0 << std::endl;
  std::cout << "Result: " << std::setprecision(20) << res*cste <<std::endl;
  

  return 0;
}
