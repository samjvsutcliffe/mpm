//! Construct a particle with id and coordinates
template <unsigned Tdim>
mpm::ParticleDamage<Tdim>::ParticleDamage(Index id, const VectorDim& coord)
    : mpm::Particle<Tdim>(id, coord) {
  //this->initialise();
  //// Clear cell ptr
  //cell_ = nullptr;
  //// Nodes
  //nodes_.clear();
  //// Set material containers
  //this->initialise_material(1);
  // Logger
  std::string logger =
      "particle" + std::to_string(Tdim) + "d_damage::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

//! Construct a particle with id, coordinates and status
template <unsigned Tdim>
mpm::ParticleDamage<Tdim>::ParticleDamage(Index id, const VectorDim& coord, bool status)
    : mpm::Particle<Tdim>(id, coord, status) {
  //this->initialise();
  //cell_ = nullptr;
  //nodes_.clear();
  //// Set material containers
  //this->initialise_material(1);
  //! Logger
  std::string logger =
      "particle" + std::to_string(Tdim) + "d_damage::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

// Initialise particle properties
template <unsigned Tdim>
void mpm::ParticleDamage<Tdim>::initialise() {
    damage_ = 0;
    damage_local_inc_ = 0;
    damage_inc_ = 0;
}

//! Compute damage increment
template <unsigned Tdim>
void mpm::Particle<Tdim>::compute_damage_increment(double dt,bool local) noexcept {
    Eigen::Matrix<double,6,1> stress = this->stress();
    double inc = 0;
    //Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,3,3>> es;  
    //es.compute(voigt_to_matrix(stress));
    //double s1 = es.eigenvalues()[0];
    //double critical_stress =
    //    (this->material())->properties_["critical_stress"];
    //double damage_rate_ = (this->material())->properties_["damage_rate"];
    //if (s1 > critical_stress)
    //{
    //  inc += (std::max(0.0,s1) / critical_stress) * damage_rate_;
    //}

    if (local)
    {
      damage_inc_local_ = inc;
    }
    else {
      damage_inc_ = inc;
    }
}

//! Apply damage increment
template <unsigned Tdim>
void mpm::Particle<Tdim>::apply_damage(double dt) noexcept {
  damage_ += damage_inc_ * dt;
  damage_ = std::min(std::max(damage_, 0.0), 1.0);
  if (damage_ > 0.) {
    //SelfAdjointEigenSolver<Matrix3d> es;
    //es.compute(stress.adjoint());
    //auto l = es.eigenvalues();
    //auto v = es.eigenvectors();
    stress_ = stress_ * (1.0 - damage_);
  }
}
