#ifndef MPM_DAMAGE_MESH_H_
#define MPM_DAMAGE_MESH_H_

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

// Eigen
#include "Eigen/Dense"
// MPI
#ifdef USE_MPI
#include "mpi.h"
#endif
// OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif
// TSL Maps
#include <tsl/robin_map.h>
// JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "cell.h"
#include "factory.h"
#include "friction_constraint.h"
#include "function_base.h"
#include "generators/injection.h"
#include "geometry.h"
#include "hdf5_particle.h"
#include "io.h"
#include "io_mesh.h"
#include "logger.h"
#include "material.h"
#include "nodal_properties.h"
#include "node.h"
#include "particle.h"
#include "particle_base.h"
#include "traction.h"
#include "vector.h"
#include "velocity_constraint.h"

namespace mpm {
//A damage node, of size h defined by the structured damage Mesh
template <unsigned Tdim>
struct DamageNode {
 public:
  SpinMutex node_mutex_;
 protected:
  std::vector<std::weak_ptr<mpm::ParticleBase>> local_list;
};

//! DamageMesh class
//! \brief Base class that stores the information about meshes
//! \details DamageMesh class which stores the particles, nodes, cells and neighbours
//! \tparam Tdim Dimension
template <unsigned Tdim>
class DamageMesh {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  // Construct a mesh with a global unique id
  //! \param[in] id Global mesh id
  //! \param[in] isoparametric DamageMesh is isoparametric
  DamageMesh(unsigned id);

  //! Default destructor
  ~DamageMesh() = default;

  //! Delete copy constructor
  DamageMesh(const DamageMesh<Tdim>&) = delete;

  //! Delete assignement operator
  DamageMesh& operator=(const DamageMesh<Tdim>&) = delete;

  int DamageMesh<1>::GetNodeRawIndex(Eigen::Matrix<int, 1, 1> index) {
    return index[1];
  }
  int DamageMesh<2>::GetNodeRawIndex(Eigen::Matrix<int, 2, 1> index) {
    return index[1] + (index[0] * mesh_size[1]);
  }
  int DamageMesh<3>::GetNodeRawIndex(Eigen::Matrix<int, 3, 1> index) {
    return index[2] + (index[1] * mesh_size[1]) + (index[0] * mesh_size[2] * mesh_size[1]);
  }
  //! Find nearest node's index
  int PositionToIndex(Eigen::Matrix<double, Tdim, 1> position) {
    Eigen::Matrix<int, Tdim, 1> index =
        (position / resolution).array().round().matrix();
    return index;
  }

  const& DamageNode GetNode(Eigen::Matrix<int, Tdim, 1> index) {
    if (index > Eigen::Zeros()) &&(index < mesh_size) {
        return raw_data[GetNodeRawIndex(index)];
    }
  };


  void PopulateMesh(Vector<ParticleBase<Tdim>>& particles) {

  };
  void AddParticleToMesh(ParticleBase<Tdim>& particle) {
    auto position = particle->position;


  };

 private:
  //! mesh id
  //! Nodal property pool
  std::shared_ptr<mpm::NodalProperties> nodal_properties_{nullptr};
  double resolution{1};
  Eigen::Matrix<int, Tdim, 1> mesh_size;
  std::vector<DamageNode> raw_data_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
  //! Maximum number of halo nodes
  unsigned nhalo_nodes_{0};
  //! Maximum number of halo nodes
  unsigned ncomms_{0};
};  // DamageMesh class
}  // namespace mpm

#include "damage_mesh.tcc"

#endif  // MPM_DAMAGE_MESH_H_
