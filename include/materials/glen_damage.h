#ifndef MPM_MATERIAL_GLEN_DAMAGE_H_
#define MPM_MATERIAL_GLEN_DAMAGE_H_

#include <limits>

#include "Eigen/Dense"

#include "material.h"
#include "glen.h"

namespace mpm {

//! LinearElasticDamage class
//! \brief Linear Elastic material model
//! \details LinearElasticDamage class stresses and strains
//! \tparam Tdim Dimension
template <unsigned Tdim>
class GlenDamage : public Glen<Tdim> {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Constructor with id
  //! \param[in] material_properties Material properties
  GlenDamage(unsigned id, const Json& material_properties);

  //! Destructor
  ~GlenDamage() override{};

  //! Delete copy constructor
  GlenDamage(const GlenDamage&) = delete;

  //! Delete assignement operator
  GlenDamage& operator=(const GlenDamage&) = delete;

  virtual Vector6d degredation_function(Vector6d stress, double damage) override { 
	  Vector6d damage_stress = stress;
	  damage_stress(3) *= std::max(1e-3,(1-damage));
	  damage_stress(4) *= std::max(1e-3,(1-damage));
	  damage_stress(5) *= std::max(1e-3,(1-damage));
	  return damage_stress;
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
  //! Elastic stiffness matrix
  using mpm::Glen<Tdim>::de_;
  //! Density
  using mpm::Glen<Tdim>::density_;
  //! Elastic youngs modulus
  using mpm::Glen<Tdim>::youngs_modulus_;
  //! Bulk modulus
  using mpm::Glen<Tdim>::bulk_modulus_;
  //! Solid poisson ratio
  using mpm::Glen<Tdim>::poisson_ratio_;
  //! Effective viscous term A
  using mpm::Glen<Tdim>::viscosity_;
  //! Viscous flow power
  using mpm::Glen<Tdim>::viscous_power_;

  double critical_stress_{0};
  double damage_rate_{0};
  double local_length_{0};
  double critical_damage_{0};

};  // GlenDamage class
}  // namespace mpm

#include "glen_damage.tcc"

#endif  // MPM_MATERIAL_LINEAR_ELASTIC_DAMAGE_H_
