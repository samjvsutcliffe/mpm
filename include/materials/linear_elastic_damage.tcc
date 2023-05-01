//! Read material properties
template <unsigned Tdim>
mpm::LinearElasticDamage<Tdim>::LinearElasticDamage(unsigned id,
                                        const Json& material_properties)
    : LinearElastic<Tdim>(id, material_properties) {
  try {
    critical_stress_ = material_properties.at("critical_stress").template get<double>();
    damage_rate_ = material_properties.at("damage_rate").template get<double>();
  } catch (Json::exception& except) {
    console_->error("Material parameter not set: {} {}\n", except.what(),
                    except.id);
  }
}


//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::LinearElasticDamage<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {
  return mpm::LinearElastic<Tdim>::compute_stress(stress,dstrain,ptr,state_vars);
}
