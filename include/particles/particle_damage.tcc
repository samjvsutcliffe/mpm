//! Construct a particle with id and coordinates
template <unsigned Tdim>
mpm::ParticleDamage<Tdim>::ParticleDamage(Index id, const VectorDim& coord)
    : mpm::Particle<Tdim>(id, coord) {
  this->initialise();
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
    //mpm::Particle<Tdim>::initialise();
    damage_ = 0;
    damage_inc_local_ = 0;
    damage_inc_ = 0;
    undamaged_stress_.setZero();
    this->scalar_properties_["damage"] = [&]() { return damage(); };
    
    //for(const auto& key_value : this->scalar_properties_) {
    //    console_->info("Keys:{}",key_value.first);
    //    //std::cout << "{" << key_value.first  << "}" << std::endl;
    //}
    //scalar_properties_["damage"] = [&]() { return 0.0; };
}

//! Compute damage increment
template <unsigned Tdim>
void mpm::ParticleDamage<Tdim>::compute_damage_increment(double dt,bool local) noexcept {
    Eigen::Matrix<double,6,1> stress = this->stress();
    double inc = 0;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,3,3>> es;  
    es.compute(this->voigt_to_matrix(stress));
    double s1 = es.eigenvalues().maxCoeff();
    double critical_stress =
        (this->material())->template property<double>("critical_stress");
    double damage_rate = (this->material())->template property<double>("damage_rate");
    //double critical_stress = 1e5;
    //double damage_rate = 0.1;
    if (s1 > critical_stress*1e-3)
    {
      inc += (std::max(0.0,s1) / critical_stress) * damage_rate;
    }

    //If we are doing a local damage update, then our final damage increment is our current local one
    //If local: set inc to local damage
    //If nonlocal: store inc in local inc
    if (!local)
    {
      damage_inc_local_ = inc;
    }
    else {
      damage_inc_ = inc;
    }
}

// Compute stress
template <unsigned Tdim>
void mpm::ParticleDamage<Tdim>::compute_stress() noexcept {
  // Check if material ptr is valid
  assert(this->material() != nullptr);
  // Calculate stress
  undamaged_stress_ =
      (this->material())
          ->compute_stress(undamaged_stress_, this->dstrain_, this,
                           &this->state_variables_[mpm::ParticlePhase::Solid]);
  this->stress_ = undamaged_stress_;
}

//! Apply damage increment
template <unsigned Tdim>
void mpm::ParticleDamage<Tdim>::apply_damage(double dt) noexcept {
  Eigen::Matrix<double,6,1> & stress = this->stress_;
  damage_ += damage_inc_ * dt;
  damage_ = std::min(std::max(damage_, 0.0), 1.0);
  if (damage_ > 0.) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,3,3>> es;  
    es.compute(this->voigt_to_matrix(stress));
    Eigen::Matrix<double,3,1> l = es.eigenvalues();
    Eigen::Matrix<double,3,3> v = es.eigenvectors();
    for(int i = 0;i < 3;++i){
        if(l[i] > 0.0){
            l[i] = l[i] * (1.0 - damage_);
        }
    }
    //stress = stress * (1.0 - damage_);
    this->stress_ = this->matrix_to_voigt(v * l.asDiagonal() * v.transpose());
  }
}
//! Delocalise damage
template <unsigned Tdim>
void mpm::ParticleDamage<Tdim>::delocalise_damage(ParticleBase<Tdim> & pother) noexcept {

}
