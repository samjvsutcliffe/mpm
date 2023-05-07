#ifndef MPM_MATERIAL_LINEAR_ELASTIC_DAMAGE_H_
#define MPM_MATERIAL_LINEAR_ELASTIC_DAMAGE_H_

#include <limits>

#include "Eigen/Dense"

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
  double critical_stress_{0};
  double damage_rate_{0};
  double local_length_{0};

};  // LinearElasticDamage class
}  // namespace mpm

#include "linear_elastic_damage.tcc"

#endif  // MPM_MATERIAL_LINEAR_ELASTIC_DAMAGE_H_
