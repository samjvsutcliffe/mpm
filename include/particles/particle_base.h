#ifndef MPM_PARTICLEBASE_H_
#define MPM_PARTICLEBASE_H_

// MPI
#ifdef USE_MPI
#include "mpi.h"
#endif

#include <array>
#include <limits>
#include <memory>
#include <vector>

#include "cell.h"
#include "data_types.h"
#include "function_base.h"
#include "hdf5_particle.h"
#include "material.h"

namespace mpm {

// Forward declaration of Material
template <unsigned Tdim>
class Material;

//! Particle phases
enum ParticlePhase : unsigned int { Solid = 0, Liquid = 1, Gas = 2 };

//! Particle type
extern std::map<std::string, int> ParticleType;
extern std::map<int, std::string> ParticleTypeName;

//! ParticleBase class
//! \brief Base class that stores the information about particleBases
//! \details ParticleBase class: id_ and coordinates.
//! \tparam Tdim Dimension
template <unsigned Tdim>
class ParticleBase {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Constructor with id and coordinates
  //! \param[in] id Particle id
  //! \param[in] coord coordinates of the particle
  ParticleBase(Index id, const VectorDim& coord);

  //! Constructor with id, coordinates and status
  //! \param[in] id Particle id
  //! \param[in] coord coordinates of the particle
  //! \param[in] status Particle status (active / inactive)
  ParticleBase(Index id, const VectorDim& coord, bool status);

  //! Destructor
  virtual ~ParticleBase(){};

  //! Delete copy constructor
  ParticleBase(const ParticleBase<Tdim>&) = delete;

  //! Delete assignement operator
  ParticleBase& operator=(const ParticleBase<Tdim>&) = delete;

  //! Initialise particle HDF5 data
  //! \param[in] particle HDF5 data of particle
  //! \retval status Status of reading HDF5 particle
  virtual bool initialise_particle(const HDF5Particle& particle) = 0;

  //! Initialise particle HDF5 data and material
  //! \param[in] particle HDF5 data of particle
  //! \param[in] material Material associated with the particle
  //! \retval status Status of reading HDF5 particle
  virtual bool initialise_particle(
      const HDF5Particle& particle,
      const std::shared_ptr<Material<Tdim>>& material) = 0;

  //! Assign material history variables
  //! \param[in] state_vars State variables
  //! \param[in] material Material associated with the particle
  //! \param[in] phase Index to indicate material phase
  //! \retval status Status of cloning HDF5 particle
  virtual bool assign_material_state_vars(
      const mpm::dense_map& state_vars,
      const std::shared_ptr<mpm::Material<Tdim>>& material,
      unsigned phase = mpm::ParticlePhase::Solid) = 0;

  //! Retrun particle data as HDF5
  //! \retval particle HDF5 data of the particle
  virtual HDF5Particle hdf5() const = 0;

  //! Return id of the particleBase
  Index id() const { return id_; }

  //! Assign coordinates
  //! \param[in] coord Assign coord as coordinates of the particleBase
  void assign_coordinates(const VectorDim& coord) { coordinates_ = coord; }

  //! Return coordinates
  //! \retval coordinates_ return coordinates of the particleBase
  VectorDim coordinates() const { return coordinates_; }

  //! Compute reference coordinates in a cell
  virtual bool compute_reference_location() = 0;

  //! Return reference location
  virtual VectorDim reference_location() const = 0;

  //! Assign cell
  virtual bool assign_cell(const std::shared_ptr<Cell<Tdim>>& cellptr) = 0;

  //! Assign cell and xi
  virtual bool assign_cell_xi(const std::shared_ptr<Cell<Tdim>>& cellptr,
                              const Eigen::Matrix<double, Tdim, 1>& xi) = 0;

  //! Assign cell id
  virtual bool assign_cell_id(Index id) = 0;

  //! Assign cell
  virtual void populate_cell_fill(const std::shared_ptr<Cell<Tdim>>& cellptr){};

  //! Return cell id
  virtual Index cell_id() const = 0;

  //! Return cell ptr status
  virtual bool cell_ptr() const = 0;

  //! Remove cell
  virtual void remove_cell() = 0;

  //! Compute shape functions
  virtual void compute_shapefn() noexcept = 0;

  //! Assign volume
  virtual bool assign_volume(double volume) = 0;

  //! Return volume
  virtual double volume() const = 0;

  //! Return size of particle in natural coordinates
  virtual VectorDim natural_size() const = 0;

  //! Compute volume of particle
  virtual void compute_volume() noexcept = 0;

  //! Update volume based on centre volumetric strain rate
  virtual void update_volume() noexcept = 0;

  //! Return mass density
  virtual double mass_density() const = 0;

  //! Compute mass of particle
  virtual void compute_mass() noexcept = 0;

  //! Map particle mass and momentum to nodes
  virtual void map_mass_momentum_to_nodes() noexcept = 0;

  //! Map multimaterial properties to nodes
  virtual void map_multimaterial_mass_momentum_to_nodes() noexcept = 0;

  //! Map multimaterial displacements to nodes
  virtual void map_multimaterial_displacements_to_nodes() noexcept = 0;

  //! Map multimaterial domain gradients to nodes
  virtual void map_multimaterial_domain_gradients_to_nodes() noexcept = 0;

  //! Assign material
  virtual bool assign_material(const std::shared_ptr<Material<Tdim>>& material,
                               unsigned phase = mpm::ParticlePhase::Solid) = 0;

  //! Return material of particle
  //! \param[in] phase Index to indicate material phase
  std::shared_ptr<Material<Tdim>> material(
      unsigned phase = mpm::ParticlePhase::Solid) const {
    return material_[phase];
  }

  //! Return material id
  //! \param[in] phase Index to indicate material phase
  unsigned material_id(unsigned phase = mpm::ParticlePhase::Solid) const {
    return material_id_[phase];
  }

  //! Return state variables
  //! \param[in] phase Index to indicate material phase
  mpm::dense_map state_variables(
      unsigned phase = mpm::ParticlePhase::Solid) const {
    return state_variables_[phase];
  }

  //! Assign status
  void assign_status(bool status) { status_ = status; }

  //! Status
  bool status() const { return status_; }

  //! Initialise properties
  virtual void initialise() = 0;

  //! Assign mass
  virtual void assign_mass(double mass) = 0;

  //! Return mass
  virtual double mass() const = 0;

  //! Return pressure
  virtual double pressure(unsigned phase = mpm::ParticlePhase::Solid) const = 0;

  //! Compute strain
  virtual void compute_strain(double dt) noexcept = 0;

  //! Compute vorticity
  virtual void compute_vorticity(double dt) noexcept = 0;

  //! Strain
  virtual Eigen::Matrix<double, 6, 1> strain() const = 0;

  //! Strain rate
  virtual Eigen::Matrix<double, 6, 1> strain_rate() const = 0;

  //! Volumetric strain of centroid
  virtual double volumetric_strain_centroid() const = 0;

  //! dvolumetric strain
  virtual double dvolumetric_strain() const = 0;

  //! Initial stress
  virtual void initial_stress(const Eigen::Matrix<double, 6, 1>&) = 0;

  //! Compute stress
  virtual void compute_stress(const float dt_) noexcept = 0;

  //! Compute damage increment
  virtual void compute_damage_increment(double dt, bool local) noexcept {};

  //! Apply damage increment
  virtual void apply_damage(double dt) noexcept {};

  //! Delocalise damage
  virtual void delocalise_damage(ParticleBase<Tdim>& pother) noexcept {};

  //! Delocalise damage
  virtual void delocalise_damage_post() noexcept {};

  //! Get damage
  virtual double damage() const { return 0; };

  //! Damage increment
  virtual double damage_inc() const { return 0; };

  //! Damage increment local
  virtual double damage_inc_local() const { return 0; };

  //! Damage increment
  virtual double damage_local_length() const { return 1; };

  //! Return stress
  virtual Eigen::Matrix<double, 6, 1> stress() const = 0;

  //! Return stress
  virtual Eigen::Matrix<double, 6, 1> stress_cauchy() const = 0;

  //! Apply a corrotational stress rate adjustment
  //! \param[in] stress increment voigt notation stress
  virtual Eigen::Matrix<double, 6, 1> objective_stress_increment(
      Eigen::Matrix<double, 6, 1> stress_inc,
      Eigen::Matrix<double, 6, 1> stress) {
    return stress;
  };

  //! Map body force
  virtual void map_body_force(const VectorDim& pgravity) noexcept = 0;

  //! Map internal force
  virtual void map_internal_force() noexcept = 0;

  //! Map particle pressure to nodes
  virtual bool map_pressure_to_nodes(
      unsigned phase = mpm::ParticlePhase::Solid) noexcept = 0;

  //! Compute pressure smoothing of the particle based on nodal pressure
  virtual bool compute_pressure_smoothing(
      unsigned phase = mpm::ParticlePhase::Solid) noexcept = 0;

  //! Assign velocity
  virtual bool assign_velocity(const VectorDim& velocity) = 0;

  //! Return velocity
  virtual VectorDim velocity() const = 0;

  //! Return displacement of the particle
  virtual VectorDim displacement() const = 0;

  //! Assign traction
  virtual bool assign_traction(unsigned direction, double traction) = 0;

  //! Return traction
  virtual VectorDim traction() const = 0;

  //! Map traction force
  virtual void map_traction_force() noexcept = 0;

  //! Compute updated position
  virtual void compute_updated_position(
      double dt, bool velocity_update = false) noexcept = 0;

  //! Return a state variable
  virtual double state_variable(
      const std::string& var,
      unsigned phase = mpm::ParticlePhase::Solid) const = 0;

  //! Return scalar data of particles
  //! \param[in] property Property string
  //! \retval data Scalar data of particle property
  virtual double scalar_data(const std::string& property) const = 0;

  //! Return vector data of particles
  //! \param[in] property Property string
  //! \retval data Vector data of particle property
  virtual VectorDim vector_data(const std::string& property) const = 0;

  //! Return tensor data of particles
  //! \param[in] property Property string
  //! \retval data Tensor data of particle property
  virtual Eigen::VectorXd tensor_data(const std::string& property) const = 0;

  //! Apply particle velocity constraints
  //! \param[in] dir Direction of particle velocity constraint
  //! \param[in] velocity Applied particle velocity constraint
  virtual void apply_particle_velocity_constraints(unsigned dir,
                                                   double velocity) = 0;

  //! Assign material id of this particle to nodes
  virtual void append_material_id_to_nodes() const = 0;

  //! Return the number of neighbour particles
  virtual unsigned nneighbours() const = 0;

  //! Assign neighbour particles
  //! \param[in] neighbours set of id of the neighbouring particles
  //! \retval insertion_status Return the successful addition of a node
  virtual void assign_neighbours(const std::vector<mpm::Index>& neighbours) = 0;

  //! Return neighbour ids
  virtual std::vector<mpm::Index> neighbours() const = 0;

  //! Type of particle
  virtual std::string type() const = 0;

  //! Serialize
  //! \retval buffer Serialized buffer data
  virtual std::vector<uint8_t> serialize() = 0;

  //! Deserialize
  //! \param[in] buffer Serialized buffer data
  //! \param[in] material Particle material pointers
  virtual void deserialize(
      const std::vector<uint8_t>& buffer,
      std::vector<std::shared_ptr<mpm::Material<Tdim>>>& materials) = 0;

  //! Reformat voigt notation symetric to matrix
  //! \param[in] voigt a symetric tensor in voigt notation
  Eigen::Matrix<double,3,3> voigt_to_matrix(Eigen::Matrix<double,6,1> voigt);
  Eigen::Matrix<double,6,1> matrix_to_voigt(Eigen::Matrix<double,3,3> mat);

  //! Reformat voigt notation vorticity to matrix
  //! \param[in] voigt an antisymetric tensor in voigt notation
  Eigen::Matrix<double,3,3> vorticity_matrix(Eigen::Matrix<double,6,1> voigt);

  //! Minus the internal force of the virtual stress field
  //! \param[in] traction Boundary traction
  //! \param[in] divergence_traction Divergence of boundary traction
  virtual void minus_virtual_stress_field(Eigen::Matrix<double, 6, 1>& traction,
                                          VectorDim& divergence_traction) {
    throw std::runtime_error(
        "Calling the base class function (minus_virtual_stress_field) "
        "in ParticleBase:: illegal operation!");
  }

 public:
     //TODO fix interface
  double reference_pressure = 0;
  bool in_damage_mesh_{false};
  //! Reference coordinates (in a cell)
  Eigen::Matrix<double, Tdim, 1> damage_position_;

 protected:
  //! particleBase id
  Index id_{std::numeric_limits<Index>::max()};
  //! coordinates
  VectorDim coordinates_;
  //! Cell id
  Index cell_id_{std::numeric_limits<Index>::max()};
  //! Status
  bool status_{true};
  //! Reference coordinates (in a cell)
  Eigen::Matrix<double, Tdim, 1> xi_;

  //! Cell
  std::shared_ptr<Cell<Tdim>> cell_;
  //! Vector of nodal pointers
  std::vector<std::shared_ptr<NodeBase<Tdim>>> nodes_;
  //! Vector of nodal pointers
  std::vector<int> nodes_local_ids_;
  //! Material
  std::vector<std::shared_ptr<Material<Tdim>>> material_;
  //! Unsigned material id
  std::vector<unsigned> material_id_;
  //! Material state history variables
  std::vector<mpm::dense_map> state_variables_;
  //! Vector of particle neighbour ids
  std::vector<mpm::Index> neighbours_;

  //Eigen::Matrix<double,6,1> matrix_to_voigt(Eigen::Matrix<double,3,3> mat) {
  //  return (Eigen::Matrix<double,6,1>() <<
  //      mat(0,0), mat(1,1),mat(2,2),
  //      mat(0,1), mat(1,2),mat(0,2)).finished();
  //};

  //Eigen::Matrix<double,3,3> voigt_to_matrix(Eigen::Matrix<double,6,1> voigt) {
  //  return (Eigen::Matrix3d() <<
  //      voigt(0), voigt(3), voigt(5),
  //      voigt(3), voigt(1), voigt(4),
  //      voigt(5), voigt(4), voigt(2)).finished();
  //};

  //Eigen::Matrix<double,3,3> vorticity_matrix(Eigen::Matrix<double,6,1> voigt) {
  //  return (Eigen::Matrix3d() <<
  //      voigt(0), voigt(3), voigt(5),
  //      -voigt(3), voigt(1), voigt(4),
  //      -voigt(5), -voigt(4), voigt(2)).finished();
  //};

};  // ParticleBase class
}  // namespace mpm

#include "particle_base.tcc"

#endif  // MPM_PARTICLEBASE_H__
