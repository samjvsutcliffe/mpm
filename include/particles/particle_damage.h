#ifndef MPM_PARTICLE_DAMAGE_H_
#define MPM_PARTICLE_DAMAGE_H_

#include <array>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "cell.h"
#include "logger.h"
#include "particle.h"
#include "particle_finite.h"
#include <unsupported/Eigen/MatrixFunctions>

namespace mpm {

//! ParticleDamage class
//! \brief Base class that stores the information about particles
//! \details ParticleDamage class: id_ and coordinates.
//! \tparam Tdim Dimension
template <unsigned Tdim>
class ParticleDamage : public ParticleFinite<Tdim> {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Define DOFs
  static const unsigned Tdof = (Tdim == 1) ? 1 : 3 * (Tdim - 1);

  //! Construct a particle with id and coordinates
  //! \param[in] id Particle id
  //! \param[in] coord coordinates of the particle
  ParticleDamage(Index id, const VectorDim& coord);

  //! Construct a particle with id, coordinates and status
  //! \param[in] id Particle id
  //! \param[in] coord coordinates of the particle
  //! \param[in] status Particle status (active / inactive)
  ParticleDamage(Index id, const VectorDim& coord, bool status);

  //! Destructor
  ~ParticleDamage() override{};

  //! Delete copy constructor
  ParticleDamage(const ParticleDamage<Tdim>&) = delete;

  //! Delete assignment operator
  ParticleDamage& operator=(const ParticleDamage<Tdim>&) = delete;

  //! Compute stress
  void compute_stress(float dt) noexcept override;

  //! Return damage at a particle
  virtual double damage() const override { return damage_; }

  //! Return damage at a particle
  double ybar() const { return damage_ybar_; }

  //! Return damage increment at a particle
  virtual double damage_inc() const override { return damage_inc_; }

  //! Return damage increment at a particle
  virtual double damage_inc_local() const override { return damage_inc_local_; }

  //! Initialise properties
  void initialise() override;

  //! Compute damage increment
  void compute_damage_increment(double dt, bool local) noexcept override;

  //! Apply damage increment
  void apply_damage(double dt) noexcept override;

  //! Delocalise damage
  void delocalise_damage(ParticleBase<Tdim> & pother) noexcept override;

  //! Delocalise damage
  void delocalise_damage_post() noexcept override;

  //! Type of particle
  std::string type() const override { return (Tdim == 2) ? "P2D_DAMAGE" : "P3D_DAMAGE"; }
  //using Particle<Tdim>::voigt_to_matrix;
 protected:
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
  //! Stress
  //using Particle<Tdim>::stress_;
  //using Particle<Tdim>::scalar_properties_;
  //! Material
  using ParticleBase<Tdim>::material_;
  //! Material id
  using ParticleBase<Tdim>::material_id_;
  //! State variables
  using ParticleBase<Tdim>::state_variables_;
  using ParticleFinite<Tdim>::reference_pressure;
  using ParticleFinite<Tdim>::stress_;
  using ParticleFinite<Tdim>::stress_kirchoff_;
  using ParticleFinite<Tdim>::dstrain_;
  using ParticleFinite<Tdim>::deformation_gradient_;
  using ParticleFinite<Tdim>::scalar_properties_;
  using ParticleFinite<Tdim>::vector_properties_;
  using ParticleFinite<Tdim>::tensor_properties_;

 protected:
  double damage_inc_local_{0.}; 
  //! True damage increment
  double damage_inc_{0.}; 
  //! Local damage
  double damage_{0.}; 
  //! YBar
  double damage_ybar_{0.}; 

  double acc_damage_{0};
  double acc_volume_{0};
  double local_length_{50};
  double critical_stress_{0};
  double critical_damage_{1.0};
  double damage_rate_{0};
  //! Stresses
  Eigen::Matrix<double, 6, 1> undamaged_stress_;

 public:
  //! Damage increment
  virtual double damage_local_length() const override { return local_length_; };

};  // Particle class
}  // namespace mpm

#include "particle_damage.tcc"

#endif  // MPM_PARTICLE_DAMAGE_H__
