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
 void reset(){
     local_list.clear();
 }
 void AddParticle(mpm::ParticleBase<Tdim> * p){
     node_mutex_.lock();
     local_list.emplace_back(p);
     node_mutex_.unlock();
 }
 bool RemoveParticle(mpm::ParticleBase<Tdim> * p){
   auto position = std::find(local_list.begin(), local_list.end(), p);
   if (position != local_list.end()) {
     node_mutex_.lock();
     local_list.erase(position);
     node_mutex_.unlock();
     return true;
   }
   return false;
 }
  template <typename Toper>
  inline void iterate_over_particles(Toper oper){
      for(mpm::ParticleBase<Tdim>* p : local_list){
          oper(*p);
      }
	  //oper();
    };
  //DamageNode() : node_mutex_(), local_list() {
  //}
  //~DamageNode() = default;
 protected:
  SpinMutex node_mutex_;
  std::vector<mpm::ParticleBase<Tdim>*> local_list;
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
  using IndexDim = Eigen::Matrix<int, Tdim, 1>;
 private:
  //! mesh id
  //! Nodal property pool
  std::shared_ptr<mpm::NodalProperties> nodal_properties_{nullptr};
  double resolution_{50};
  VectorDim mesh_size;
  VectorDim offset;
  std::vector<DamageNode<Tdim>> nodes_;
  //! Maximum number of halo nodes
  unsigned nhalo_nodes_{0};
  //! Maximum number of halo nodes
  unsigned ncomms_{0};

 public:
  //! LoggerWe couldn't find particle to remove!
  std::unique_ptr<spdlog::logger> console_;
  // Construct a mesh with a global unique id
  //! \param[in] id Global mesh id
  //! \param[in] isoparametric DamageMesh is isoparametric
  DamageMesh() = default;
  DamageMesh(VectorDim min, VectorDim max, double resolution) : resolution_{resolution}, nodes_(){
      offset = min;
    mesh_size = ((max - min) / resolution_).array().ceil().matrix();
      ////Default initalise nodes
      int size = mesh_size.prod();
      nodes_ = std::vector<DamageNode<Tdim>>(size);
      std::string logger = "damage_mesh";
	  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
      //nodes_.resize(size);
  };

  //! Default destructor
  ~DamageMesh() = default;

  //! Delete copy constructor
  DamageMesh(const DamageMesh<Tdim>&) = delete;

  //! Delete assignement operator
  DamageMesh& operator=(const DamageMesh<Tdim>&) = delete;

  bool InBounds(IndexDim index){
    for (int i = 0; i < Tdim; ++i)
    {
      if (index(i) < 0 || index(i) >= mesh_size(i)) {
        return false;
      }
    }
    return true;
    //return (index >= Eigen::Zeros()) && (index < mesh_size);
  }
  //! Find nearest node's index
  IndexDim PositionToIndex(Eigen::Matrix<double, Tdim, 1> position) {
    IndexDim index = (((position - offset) / resolution_)
                                            .array()
                                            .round()
                         .matrix())
                         .template cast<int>();
    return index;
  };
  //! Caculate real-world position from index
  VectorDim IndexToPosition(Eigen::Matrix<double, Tdim, 1> index) {
    return (index * resolution_) + offset;
  };
  //Calculate displaced index of a node from its index
  inline int GetNodeRawIndex(Eigen::Matrix<int, Tdim, 1> index);
  //! Get a node from an index
  DamageNode<Tdim>& GetNode(Eigen::Matrix<int, Tdim, 1> index) {
    //if (index > Eigen::Zeros()) &&(index < mesh_size) {
    //    return raw_data[GetNodeRawIndex(index)];
    //}
    return nodes_[GetNodeRawIndex(index)];
  };

  template <typename Toper>
  inline void iterate_over_nodes(Toper oper){
#pragma omp parallel for schedule(runtime)
      for (auto nitr = nodes_.begin(); nitr != nodes_.end(); ++nitr) oper(*nitr);
    };


  template <typename Toper>
  inline void iterate_over_neighbours(Toper oper,mpm::ParticleBase<Tdim> & p,double distance);

  void reset(){
      iterate_over_nodes(
          std::bind(&mpm::DamageNode<Tdim>::reset, std::placeholders::_1));
  };

  void PopulateMesh(Vector<ParticleBase<Tdim>>& particles) {
#pragma omp parallel for schedule(runtime)
      for (auto p = particles.begin(); p != p.end(); ++p)
      {
       //AddParticle(*p);
      }
  };
  void AddParticle(const std::shared_ptr<ParticleBase<Tdim>> & particle) {
    auto position = particle->coordinates();
    particle->in_damage_mesh_ = true;
    particle->damage_position_ = position;
    Eigen::Matrix<int, Tdim, 1> index = PositionToIndex(position);
    DamageNode<Tdim> & node = GetNode(index);
    node.AddParticle(particle.get());
  };
  void UpdateParticle(const std::shared_ptr<ParticleBase<Tdim>>& particle) {
    if (particle->in_damage_mesh_) {
      auto current_position = particle->coordinates();
      double diff =
          (current_position - particle->damage_position_).squaredNorm();
      if (diff > resolution_ / 4) {
        RemoveParticle(particle);
        AddParticle(particle);
      }
    } else {
      AddParticle(particle);
    }
  };
  void RemoveParticle(const std::shared_ptr<ParticleBase<Tdim>>& particle) {
	//console_->info("Remove particle from damage mesh");
    if (particle->in_damage_mesh_) {
      // auto position = particle->coordinates();
      auto position = particle->damage_position_;
      Eigen::Matrix<int, Tdim, 1> index = PositionToIndex(position);
      DamageNode<Tdim>& node = GetNode(index);
      // Happy path
      if (!node.RemoveParticle(particle.get())) {
        // Sad path - we didn't find the particle at the current pos - start
        console_->info("We couldn't find particle to remove!");
        // searching grid #pragma omp parallel for schedule(runtime)
        //for (auto nitr = nodes_.begin(); nitr != nodes_.end(); ++nitr) {
		  //if (nitr->RemoveParticle(particle.get()))
		  //{
          //  particle->in_damage_mesh_ = false;
          //  break;
		  //}
        //}
      }
	  else {
	    particle->in_damage_mesh_ = false;
	  }
    }
  };
};  // DamageMesh class
}  // namespace mpm

#include "damage_mesh.tcc"

#endif
