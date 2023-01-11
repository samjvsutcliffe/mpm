//! Constructor with material properties
template <unsigned Tdim>
mpm::Maxwell<Tdim>::Maxwell(unsigned id, const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    density_ = material_properties.at("density").template get<double>();
    youngs_modulus_ =
        material_properties.at("youngs_modulus").template get<double>();
    poisson_ratio_ =
        material_properties.at("poisson_ratio").template get<double>();
    viscosity_ =
        material_properties.at("viscosity").template get<double>();

    properties_ = material_properties;
    compute_elastic_tensor();
  } catch (Json::exception& except) {
    console_->error("Material parameter not set: {} {}\n", except.what(),
                    except.id);
  }
}

template <unsigned Tdim>
bool mpm::Maxwell<Tdim>::compute_elastic_tensor() {
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

//template <>
//Eigen::Matrix<double, 6, 1> mpm::Maxwell<2>::compute_stress(
//    const Vector6d& stress, const Vector6d& dstrain, const ParticleBase<2>* ptr,
//    mpm::dense_map* state_vars) {
//
//  // Get strain rate
//  const auto& strain_rate = ptr->strain_rate();
//  const double volumetric_strain_rate = strain_rate(0) + strain_rate(1);
//  // Update pressure
//  //
//  // Volumetric stress component
//  const double relaxation_constant = youngs_modulus_ * ((*state_vars).at("dt")) /
//                                     (2.0 * (1.0 - poisson_ratio_) * viscosity_);
//  Vector6d dev_stress = stress;
//  const double trace_stress = (stress(0) + stress(1) + stress(2)) / 3.0;
//  for (int i = 0; i < 3; ++i) {
//    dev_stress(i) -= trace_stress;
//  }
//  // Update stress component
//  Eigen::Matrix<double, 6, 1> pstress = stress;
//  pstress += this->de_ * dstrain;
//  pstress -= dev_stress * relaxation_constant;
//  
//  return pstress;
//}
//! Compute stress in 2D

template <>
Eigen::Matrix<double, 6, 1> mpm::Maxwell<2>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain, const ParticleBase<2>* ptr,
    mpm::dense_map* state_vars) {

  // Get strain rate
  const auto& strain_rate = ptr->strain_rate();
  const double volumetric_strain_rate = strain_rate(0) + strain_rate(1);
  // Update pressure
  //
  // Volumetric stress component
  //const double rho = youngs_modulus_ / (2.0 * (1.0 - poisson_ratio_) * viscosity_);
  const double relaxation_constant = youngs_modulus_ * ((*state_vars).at("dt")) / (2.0 * (1.0 - poisson_ratio_) * viscosity_);
  const double rho_exp = std::exp(-relaxation_constant);
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
  const double new_pressure = trace_stress + pressure_inc;
  // Update stress component
  Eigen::Matrix<double, 6, 1> pstress = Eigen::Matrix<double, 6, 1>::Zero();
  pstress(0) = new_pressure;
  pstress(1) = new_pressure;
  pstress(2) = new_pressure;
  pstress += (dev_stress.array() * rho_exp).matrix();
  pstress += (elastic_inc.array() * ((1 - rho_exp)/relaxation_constant)).matrix();
  
  return pstress;
}
//! Compute stress in 3D
template <>
Eigen::Matrix<double, 6, 1> mpm::Maxwell<3>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain, const ParticleBase<3>* ptr,
    mpm::dense_map* state_vars) {

  // Get strain rate
  const auto& strain_rate = ptr->strain_rate();
  const double volumetric_strain_rate =
      strain_rate(0) + strain_rate(1) + strain_rate(2);
  // Update stress component
  Eigen::Matrix<double, 6, 1> pstress = stress;
  pstress += this->de_ * dstrain;
  //pstress -= stress * relaxation_constant;
  /*
  * This is the incompressible, dev only shear aproximation
  pstress(0) += volumetric_component;
  pstress(1) += volumetric_component;
  pstress(2) += volumetric_component;
  pstress(3) -= (relaxation_constant * stress(3));
  pstress(4) -= (relaxation_constant * stress(4));
  pstress(5) -= (relaxation_constant * stress(5));
  */
  return pstress;
}
