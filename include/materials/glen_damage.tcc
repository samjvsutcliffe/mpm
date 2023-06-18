//! Read material properties
template <unsigned Tdim>
mpm::GlenDamage<Tdim>::GlenDamage(unsigned id, const Json& material_properties)
    : mpm::Glen<Tdim>(id, material_properties) {
  try {
    critical_stress_ = material_properties.at("critical_stress").template get<double>();
    damage_rate_ = material_properties.at("damage_rate").template get<double>();
    local_length_ = material_properties.at("local_length").template get<double>();
    critical_damage_ = material_properties.at("critical_damage").template get<double>();
  } catch (Json::exception& except) {
    console_->error("Material parameter not set: {} {}\n", except.what(),
                    except.id);
  }
}


//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::GlenDamage<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {

  // Get strain rate
  const double dt = (*state_vars).at("dt");
  const auto& strain_rate = ptr->strain_rate();
  const auto& strain = ptr->strain();

  Vector6d log_strain_rate = dstrain * (1 / dt);
  Vector6d eng_strain_rate =
      (log_strain_rate.array() * ptr->strain().array().exp()).matrix();
  const double viscosity =
      compute_glen_viscosity_strain(eng_strain_rate, viscosity_, viscous_power_);

  Vector6d dev_stress = stress;
  const double trace_stress = (stress(0) + stress(1) + stress(2)) / 3.0;
  for (int i = 0; i < 3; ++i) {
    dev_stress(i) -= trace_stress;
  }
  // Update stress component
  Eigen::Matrix<double, 6, 1> dstress = this->degredation_function(this->de_ * dstrain,ptr->damage());
  if (viscosity > 0.0) {
	  const double relaxation_constant = youngs_modulus_ * ((*state_vars).at("dt")) / (2.0 * (1.0 - poisson_ratio_) * viscosity);
	  dstress -= dev_stress * relaxation_constant;
  }
  return stress + ptr->objective_stress_increment(dstress,stress);
}
