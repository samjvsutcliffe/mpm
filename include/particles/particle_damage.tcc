//! Construct a particle with id and coordinates
template <unsigned Tdim>
mpm::ParticleDamage<Tdim>::ParticleDamage(Index id, const VectorDim& coord)
    : mpm::ParticleFinite<Tdim>(id, coord) {
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
	//local_length_ = (this->material())->template property<double>("local_length");
	//critical_stress_ = (this->material())->template property<double>("critical_stress");
	//damage_rate_ = (this->material())->template property<double>("damage_rate");
	//critical_damage_ = (this->material())->template property<double>("critical_damage");
	//local_length_ = (this->material())->template property<double>("local_length");
	//critical_stress_ = (this->material())->template property<double>("critical_stress");
	//damage_rate_ = (this->material())->template property<double>("damage_rate");
	//critical_damage_ = (this->material())->template property<double>("critical_damage");
}

//! Construct a particle with id, coordinates and status
template <unsigned Tdim>
mpm::ParticleDamage<Tdim>::ParticleDamage(Index id, const VectorDim& coord, bool status)
    : mpm::ParticleFinite<Tdim>(id, coord, status) {
  //this->initialise();
  //cell_ = nullptr;
  //nodes_.clear();
  //// Set material containers
  //this->initialise_material(1);
  //! Logger
  std::string logger =
      "particle" + std::to_string(Tdim) + "d_damage::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
	//local_length_ = (this->material())->template property<double>("local_length");
	//critical_stress_ = (this->material())->template property<double>("critical_stress");
	//damage_rate_ = (this->material())->template property<double>("damage_rate");
	//critical_damage_ = (this->material())->template property<double>("critical_damage");
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
    this->scalar_properties_["ybar"] = [&]() { return ybar(); };
    this->scalar_properties_["damage_inc"] = [&]() { return damage_inc(); };
	//local_length_ = (this->material())->template property<double>("local_length");
	//critical_stress_ = (this->material())->template property<double>("critical_stress");
	//damage_rate_ = (this->material())->template property<double>("damage_rate");
	//critical_damage_ = (this->material())->template property<double>("critical_damage");
    
    //try {
    //    local_length_ = (this->material())->template property<double>("local_length");
    //    critical_stress_ = (this->material())->template property<double>("critical_stress");
    //    damage_rate_ = (this->material())->template property<double>("damage_rate");
    //} catch (Json::exception& except) {
    //    console_->error("Material parameter not set: {} {}\n", except.what(), except.id);
    //}
    //for(const auto& key_value : this->scalar_properties_) {
    //    console_->info("Keys:{}",key_value.first);
    //    //std::cout << "{" << key_value.first  << "}" << std::endl;
    //}
    //scalar_properties_["damage"] = [&]() { return 0.0; };
}

//! Compute damage increment
template <unsigned Tdim>
void mpm::ParticleDamage<Tdim>::compute_damage_increment(double dt,bool local) noexcept {
	local_length_ = (this->material())->template property<double>("local_length");
	critical_stress_ = (this->material())->template property<double>("critical_stress");
	damage_rate_ = (this->material())->template property<double>("damage_rate");
	critical_damage_ = (this->material())->template property<double>("critical_damage");



    Eigen::Matrix<double,6,1> stress = this->stress();
    double inc = 0;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,3,3>> es;  
    es.compute(this->voigt_to_matrix(stress));
    const Eigen::Matrix<double,3,1> l = es.eigenvalues();
    const double s1 = es.eigenvalues().maxCoeff();
    const double s2 = es.eigenvalues().minCoeff();
    const double stress_crit = s1 - reference_pressure;
    //const double stress_crit = s1 - reference_pressure;
    //const double stress_crit = std::sqrt(0.5 * (std::pow(std::max(l[2],0.0) - std::max(l[1],0.0), 2) + 
    //                                            std::pow(std::max(l[1],0.0) - std::max(l[0],0.0), 2) +
    //                                            std::pow(std::max(l[0],0.0) - std::max(l[2],0.0), 2)));
    //console_->info("Stress critical {}\n", local_length_);
    //damage_ybar_ = stress_crit;
    if (stress_crit > 0)
    {
      const double integrity = 1.0d - damage_;
      if (integrity > 0) {
        inc += (stress_crit / integrity);
      }
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
    //damage_inc_ = inc;
}

// Compute stress
template <unsigned Tdim>
void mpm::ParticleDamage<Tdim>::compute_stress(float dt) noexcept {
  // Check if material ptr is valid
  assert(this->material() != nullptr);
  // Calculate stress
  stress_ = undamaged_stress_;
  mpm::ParticleFinite<Tdim>::compute_stress(dt);
  //undamaged_stress_ =
  //    (this->material())
  //        ->compute_stress(undamaged_stress_, this->dstrain_, this,
  //                         &this->state_variables_[mpm::ParticlePhase::Solid]);

  this->undamaged_stress_ = stress_ * deformation_gradient_.determinant();
  //this->stress_ = undamaged_stress_;
}

//! Apply damage increment
template <unsigned Tdim>
void mpm::ParticleDamage<Tdim>::apply_damage(double dt) noexcept {
  Eigen::Matrix<double,6,1> & stress = this->stress_;
    
  damage_ybar_ = damage_inc_;
  //Take our averaged y bar and apply damage rate to it
  damage_inc_ = std::pow(std::max(0.0,damage_inc_ - critical_stress_),3.0) * damage_rate_;
  damage_ += damage_inc_ * dt;
  damage_ = std::min(std::max(damage_, 0.0), 1.0);
  if (damage_ > critical_damage_) {
    damage_ = 1;
    damage_inc_ = 0;
  }
  if (damage_ > 0.) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,3,3>> es;  
    es.compute(this->voigt_to_matrix(stress));
    Eigen::Matrix<double,3,1> l = es.eigenvalues();
    Eigen::Matrix<double,3,3> v = es.eigenvectors();
    for(int i = 0;i < 3;++i){
        double esii = l[i] - reference_pressure;
        if(esii > 0.0){
            l[i] = (esii * (1.0 - damage_)) + reference_pressure;
        }
    }
    //reference_pressure = 0;
    this->stress_ = this->matrix_to_voigt(v * l.asDiagonal() * v.transpose());
  }
}
//! Delocalise damage
template <unsigned Tdim>
void mpm::ParticleDamage<Tdim>::delocalise_damage(ParticleBase<Tdim> & pother) noexcept {
  //Throws if particle damage isn't 
  VectorDim dist = this->coordinates() - pother.coordinates();
  if (pother.damage() < 1) {
    //console_->info("Actually delocal");
    const double distance = dist.dot(dist);
    const double weight = std::exp(((-4.0/(local_length_*local_length_)) * distance)); 
    const double weighted_volume = pother.volume() * weight;
    acc_damage_ += pother.damage_inc_local() * weighted_volume;
    acc_volume_ += weighted_volume;
    //if (pother.damage_inc_local() > 0) {
    //    console_->info("Delocalised inc: {}, dist: {}, weight: {}",pother.damage_inc_local(),std::sqrt(distance),weight);
    //}
  }
}
//! Delocalise damage post stp
template <unsigned Tdim>
void mpm::ParticleDamage<Tdim>::delocalise_damage_post() noexcept {
  //Throws if particle damage isn't 
    if (acc_volume_ > 0)
    {
        damage_inc_ = acc_damage_ / acc_volume_;
    }
    acc_volume_ = 0;
    acc_damage_ = 0;
}
