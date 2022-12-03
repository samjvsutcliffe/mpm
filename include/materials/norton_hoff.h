#ifndef MPM_MATERIAL_Norton_Hoff_H_
#define MPM_MATERIAL_Norton_Hoff_H_

#include <limits>

#include "Eigen/Dense"

#include "material.h"
#include "linear_elastic.h"

namespace mpm {

//! Norton_Hoff class
//! \brief Norton_Hoff fluid material model
//! \details Norton_Hoff class stresses and strains
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Norton_Hoff : public Material<Tdim> {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Constructor with id and material properties
  //! \param[in] id Material ID
  //! \param[in] material_properties Material properties
  Norton_Hoff(unsigned id, const Json& material_properties);

  //! Destructor
  ~Norton_Hoff() override{};

  //! Delete copy constructor
  Norton_Hoff(const Norton_Hoff&) = delete;

  //! Delete assignement operator
  Norton_Hoff& operator=(const Norton_Hoff&) = delete;

  //! Initialise history variables
  //! \retval state_vars State variables with history
  mpm::dense_map initialise_state_variables() override {
    mpm::dense_map state_vars;
    return state_vars;
  }

  //! State variables
  std::vector<std::string> state_variables() const override { return {}; }

  //! Compute stress
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \retval updated_stress Updated value of stress
  Vector6d compute_stress(const Vector6d& stress, const Vector6d& dstrain,
                          const ParticleBase<Tdim>* ptr,
                          mpm::dense_map* state_vars) override;
  bool compute_elastic_tensor();

 protected:
  //! material id
  using Material<Tdim>::id_;
  //! Material properties
  using Material<Tdim>::properties_;
  //! Logger
  using Material<Tdim>::console_;

 private:

  //! Elastic stiffness matrix
  Matrix6x6 de_;
  //! Density
  double density_{std::numeric_limits<double>::max()};
  //! Elastic youngs modulus
  double youngs_modulus_{std::numeric_limits<double>::max()};
  //! Solid poisson ratio
  double poisson_ratio_{std::numeric_limits<double>::max()};
  //! Effective viscous term A
  double viscosity_{std::numeric_limits<double>::max()};
  //! Viscous flow power
  double viscous_power_{std::numeric_limits<double>::max()};
};  // Norton_Hoff class
}  // namespace mpm

#include "norton_hoff.tcc"

#endif  // MPM_MATERIAL_Norton_Hoff_H_
