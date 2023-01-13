//! Constructor with material properties
template <unsigned Tdim>
mpm::Glen<Tdim>::Glen(unsigned id, const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    density_ = material_properties.at("density").template get<double>();
    
    //bulk_modulus_ =
    //    material_properties.at("bulk_modulus").template get<double>();
    youngs_modulus_ =
        material_properties.at("youngs_modulus").template get<double>();
    poisson_ratio_ =
        material_properties.at("poisson_ratio").template get<double>();
    viscosity_ = 
        material_properties.at("viscosity").template get<double>();
    viscous_power_ = 
        material_properties.at("viscous_power").template get<double>();

    properties_ = material_properties;
    compute_elastic_tensor();
  } catch (Json::exception& except) {
    console_->error("Material parameter not set: {} {}\n", except.what(),
                    except.id);
  }
}

template <unsigned Tdim>
bool mpm::Glen<Tdim>::compute_elastic_tensor() {
  // Shear modulus
  const double bulk_modulus_ = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));
  const double G = youngs_modulus_ / (2.0 * (1. + poisson_ratio_));

  const double a1 = bulk_modulus_ + (4.0 / 3.0) * G;
  const double a2 = bulk_modulus_ - (2.0 / 3.0) * G;

  // clang-format off
  // compute elasticityTensor
  de_(0,0)=a1;    de_(0,1)=a2;    de_(0,2)=a2;    de_(0,3)=0;    de_(0,4)=0;    de_(0,5)=0;
  de_(1,0)=a2;    de_(1,1)=a1;    de_(1,2)=a2;    de_(1,3)=0;    de_(1,4)=0;    de_(1,5)=0;
  de_(2,0)=a2;    de_(2,1)=a2;    de_(2,2)=a1;    de_(2,3)=0;    de_(2,4)=0;    de_(2,5)=0;
  de_(3,0)= 0;    de_(3,1)= 0;    de_(3,2)= 0;    de_(3,3)=G;    de_(3,4)=0;    de_(3,5)=0;
  de_(4,0)= 0;    de_(4,1)= 0;    de_(4,2)= 0;    de_(4,3)=0;    de_(4,4)=G;    de_(4,5)=0;
  de_(5,0)= 0;    de_(5,1)= 0;    de_(5,2)= 0;    de_(5,3)=0;    de_(5,4)=0;    de_(5,5)=G;
  // clang-format on
  return true;
}

template<unsigned tdim>
double mpm::Glen<tdim>::compute_glen_viscosity(
    Eigen::Matrix<double, 6, 1> stress,
    double visc_factor,
    double visc_power) {
    const double trace_stress = (stress(0) + stress(1) + stress(2)) / 3.0;
    const Eigen::Array<double,6,1> second_invar_mult = (Eigen::Array<double, 6, 1>() << 0.5, 0.5, 0.5, 1, 1, 1).finished();
    Eigen::Matrix<double, 6, 1> dev_stress = stress;
    for (int i = 0; i < 3; ++i) {
      dev_stress(i) -= trace_stress;
    }
    double effective_stress =
        (dev_stress.dot((dev_stress.array() * second_invar_mult).matrix()));
    if (effective_stress > 0.0)
    {
      return 1.0 / (2.0 * visc_factor *
                    std::pow(effective_stress, 0.5 * (visc_power - 1.0)));
    } else {
      return -1.0;
    }
}

template<unsigned tdim>
double mpm::Glen<tdim>::compute_glen_viscosity_strain(
    Eigen::Matrix<double, 6, 1> strain,
    double visc_factor,
    double visc_power) {
    const double trace_strain = (strain(0) + strain(1) + strain(2)) / 3.0;
    const Eigen::Array<double,6,1> second_invar_mult = (Eigen::Array<double, 6, 1>() << 0.5, 0.5, 0.5, 1, 1, 1).finished();
    Eigen::Matrix<double, 6, 1> dev_strain = strain;
    for (int i = 0; i < 3; ++i) {
      dev_strain(i) -= trace_strain;
    }
    double effective_strain =
        (dev_strain.dot((dev_strain.array() * second_invar_mult).matrix()));
    return 0.5 * visc_factor *
           std::pow(effective_strain + 1e-13, 0.5 * ((1 / visc_power) - 1));
    //if (effective_stress > 0.0)
    //{
    //  return 1.0 / (2.0 * visc_factor *
    //                std::pow(effective_stress, 0.5 * (visc_power - 1.0)));
    //} else {
    //  return -1.0;
    //}
}

/*
template <>
Eigen::Matrix<double, 6, 1> mpm::Glen<2>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain, const ParticleBase<2>* ptr,
    mpm::dense_map* state_vars) {

  // Get strain rate
  const auto& strain_rate = ptr->strain_rate();

  const double volumetric_strain_rate = ptr->dvolumetric_strain();
      //strain_rate(0) + strain_rate(1) + strain_rate(2);

  const double trace_stress = (stress(0) + stress(1) + stress(2)) / 3.0;

  const double pressure_increment_ = bulk_modulus_ * volumetric_strain_rate;
  // Update stress component
  Eigen::Matrix<double, 6, 1> pstress;
  pstress(0) = stress(0) + pressure_increment_;
  pstress(1) = stress(1) + pressure_increment_;
  pstress(2) = stress(2) + pressure_increment_;
  pstress(3) = 0;
  pstress(4) = 0;
  pstress(5) = 0;

  const Eigen::Array<double,6,1> second_invar_mult =
      (Eigen::Array<double, 6, 1>() << 0.5, 0.5, 0.5, 1, 1, 1).finished();

  auto dev_strain = strain_rate;
  for (int i = 0; i < 3; ++i) {
    dev_strain(i) -= volumetric_strain_rate/3;
  }
  const double effective_strain = (dev_strain.dot((dev_strain.array() * second_invar_mult).matrix()));
  if (effective_strain > 0) {
  Vector6d viscous_stress =
      ((*state_vars).at("dt") * viscosity_ *
       std::pow(effective_strain, (0.5*((1.0/viscous_power_) - 1)))) * dev_strain;
    pstress(3) = viscous_stress(3);
  }
  //pstress(4) = 0;
  //pstress(5) = 0;
  return pstress;
}
*/
template <>
Eigen::Matrix<double, 6, 1> mpm::Glen<2>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain, const ParticleBase<2>* ptr,
    mpm::dense_map* state_vars) {

  // Get strain rate
  const double dt = (*state_vars).at("dt");
  const auto& strain_rate = ptr->strain_rate();
  const auto& strain = ptr->strain();
  const double volumetric_strain_rate = strain_rate(0) + strain_rate(1);
  // Update pressure
  //
  // Volumetric stress component
  //const double rho = youngs_modulus_ / (2.0 * (1.0 - poisson_ratio_) * viscosity_);
  //const double viscosity =
  //    compute_glen_viscosity(stress, viscosity_, viscous_power_);
  const double viscosity =
      compute_glen_viscosity_strain(dstrain * (1/dt), viscosity_, viscous_power_);

  Vector6d dev_stress = stress;
  const double trace_stress = (stress(0) + stress(1) + stress(2)) / 3.0;
  for (int i = 0; i < 3; ++i) {
    dev_stress(i) -= trace_stress;
  }
  Eigen::Matrix<double, 6, 1> elastic_inc = this->de_ * dstrain;
  const double pressure_inc = (elastic_inc(0) + elastic_inc(1) + elastic_inc(2)) / 3.0;
  for (int i = 0; i < 3; ++i) {
    elastic_inc(i) -= pressure_inc;
  }
  const double bulk_modulus_ = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));
  //const double new_pressure = trace_stress + pressure_inc;
  const double new_pressure = bulk_modulus_ * (strain(0) + strain(1) + strain(2));
  // Update stress component
  Eigen::Matrix<double, 6, 1> pstress = Eigen::Matrix<double, 6, 1>::Zero();
  pstress(0) = new_pressure;
  pstress(1) = new_pressure;
  pstress(2) = new_pressure;
  //if (viscosity > 0.0) 
  //{
    const double rho = (2.0 * (1.0 - poisson_ratio_) * viscosity) / youngs_modulus_;
    const double rho_t = dt / rho;
    const double rho_exp = std::exp(-rho_t);
    const double lambda = (1 - rho_exp) * rho / dt;
    pstress += (dev_stress.array() * rho_exp).matrix();
    pstress += (elastic_inc.array() * lambda).matrix();
  //} else {
  //  pstress += elastic_inc;
  //}

  /*
  Vector6d dev_stress = stress;
  const double trace_stress = (stress(0) + stress(1) + stress(2)) / 3.0;
  for (int i = 0; i < 3; ++i) {
    dev_stress(i) -= trace_stress;
  }
  // Update stress component
  Eigen::Matrix<double, 6, 1> pstress = stress;
  pstress += this->de_ * dstrain;
  if (viscosity > 0.0) {
	  const double relaxation_constant = youngs_modulus_ * ((*state_vars).at("dt")) / (2.0 * (1.0 - poisson_ratio_) * viscosity);
	  pstress -= dev_stress * relaxation_constant;
  }
  */
  return pstress;
}
template <>
Eigen::Matrix<double, 6, 1> mpm::Glen<3>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain, const ParticleBase<3>* ptr,
    mpm::dense_map* state_vars) {

  // Get strain rate
  const auto& strain_rate = ptr->strain_rate();

  const double volumetric_strain_rate =
      strain_rate(0) + strain_rate(1) + strain_rate(2);
  const double trace_stress = (stress(0) + stress(1) + stress(2)) / 3.0;
  const Eigen::Array<double,6,1> second_invar_mult =
      (Eigen::Array<double, 6, 1>() << 0.5, 0.5, 0.5, 1, 1, 1).finished();
  Vector6d dev_stress = stress;
  for (int i = 0; i < 3; ++i) {
    dev_stress(i) -= trace_stress;
  }
  Vector6d viscous_strain_rate =
      ((*state_vars).at("dt") * viscosity_ *
       std::pow((dev_stress.dot(
                    (dev_stress.array() * second_invar_mult).matrix())),
                (0.5*(viscous_power_ - 1)))) *
      dev_stress;

  // Update stress component
  Eigen::Matrix<double, 6, 1> pstress = stress;
  pstress += this->de_ * (dstrain - viscous_strain_rate);
  return pstress;
}
