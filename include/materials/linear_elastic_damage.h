#ifndef MPM_MATERIAL_LINEAR_ELASTIC_DAMAGE_H_
#define MPM_MATERIAL_LINEAR_ELASTIC_DAMAGE_H_

#include <limits>

#include "Eigen/Dense"

#include <string>
#include "material.h"
#include "linear_elastic.h"

namespace mpm {

//! LinearElasticDamage class
//! \brief Linear Elastic material model
//! \details LinearElasticDamage class stresses and strains
//! \tparam Tdim Dimension
template <unsigned Tdim>
class LinearElasticDamage : public LinearElastic<Tdim> {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Constructor with id
  //! \param[in] material_properties Material properties
  LinearElasticDamage(unsigned id, const Json& material_properties);

  //! Destructor
  ~LinearElasticDamage() override{};

  //! Delete copy constructor
  LinearElasticDamage(const LinearElasticDamage&) = delete;

  //! Delete assignement operator
  LinearElasticDamage& operator=(const LinearElasticDamage&) = delete;

  //Linear elastic degredation orthotropic
  virtual Vector6d degredation_function(Vector6d stress, double damage) override { 
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

  //! Compute stress
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \retval updated_stress Updated value of stress
  Vector6d compute_stress(const Vector6d& stress, const Vector6d& dstrain,
                          const ParticleBase<Tdim>* ptr,
                          mpm::dense_map* state_vars) override;

 protected:
  //! material id
  using Material<Tdim>::id_;
  //! Material properties
  using Material<Tdim>::properties_;
  //! Logger
  using Material<Tdim>::console_;
  using LinearElastic<Tdim>::de_;
  enum class DegredationFunctions {
    Isotropic,
    Deviatoric
  };
  double critical_stress_{0};
  double damage_rate_{0};
  double local_length_{0};
  double critical_damage_{0};
  DegredationFunctions deg_{DegredationFunctions::Isotropic};

};  // LinearElasticDamage class
}  // namespace mpm

#include "linear_elastic_damage.tcc"

#endif  // MPM_MATERIAL_LINEAR_ELASTIC_DAMAGE_H_
