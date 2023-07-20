#ifndef MPM_MATERIAL_DAMAGE_H_
#define MPM_MATERIAL_DAMAGE_H_

#include <limits>

#include "Eigen/Dense"

#include "material.h"
class Damage : public Material<Tdim> {
  //Linear elastic degredation orthotropic
  virtual Vector6d degredation_function(Vector6d stress, double damage, ) override { 
    //Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,3,3>> es;  
    //es.compute(voigt_to_matrix(stress));
    //Eigen::Matrix<double,3,1> l = es.eigenvalues();
    //Eigen::Matrix<double,3,3> v = es.eigenvectors();
    //for(int i = 0;i < 3;++i){
    //    double esii = l[i] - reference_pressure;
    //    if(esii > 0.0){
    //        l[i] = (esii * (1.0 - damage_)) + reference_pressure;
    //    }
    //}
    ////reference_pressure = 0;
    //return matrix_to_voigt(v * l.asDiagonal() * v.transpose());
    //return stress * damage;
    if (deg_ == DegredationFunctions::Isotropic) {
	  return stress * std::max(1e-3,(1-damage));
    } else if (deg_ == DegredationFunctions::Deviatoric) {
	  Vector6d damage_stress = stress;
	  damage_stress(3) *= std::max(1e-3,(1-damage));
	  damage_stress(4) *= std::max(1e-3,(1-damage));
	  damage_stress(5) *= std::max(1e-3,(1-damage));
	  return damage_stress;
    } else {
      return stress;
    }
  };
  double critical_stress_{0};
  double damage_rate_{0};
  double local_length_{0};
  double critical_damage_{0};
  class enum DegredationFunctions {
    Isotropic,
    Deviatoric,
    Orthotropic
  };
  DegredationFunctions deg_{DegredationFunctions::Isotropic};
}
