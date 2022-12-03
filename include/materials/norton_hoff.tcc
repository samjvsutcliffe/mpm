//! Constructor with material properties
template <unsigned Tdim>
mpm::Norton_Hoff<Tdim>::Norton_Hoff(unsigned id, const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    density_ = material_properties.at("density").template get<double>();
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
bool mpm::Norton_Hoff<Tdim>::compute_elastic_tensor() {
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

//! Compute stress in 3D
template <>
Eigen::Matrix<double, 6, 1> mpm::Norton_Hoff<3>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain, const ParticleBase<3>* ptr,
    mpm::dense_map* state_vars) {

  // Get strain rate
  const auto& strain_rate = ptr->strain_rate();

  const double volumetric_strain_rate =
      strain_rate(0) + strain_rate(1) + strain_rate(2);
  const double trace_stress = (stress(0) + stress(1) + stress(2)) / 3.0;
  Vector6d dev_stress = stress;
  for (int i = 0; i < 3; ++i) {
    dev_stress(i) -= trace_stress;
  }
  Vector6d viscous_strain_rate =
      (*state_vars).at("dt") * viscosity_ *
      std::pow((0.5 * dev_stress.dot(dev_stress)), (viscous_power_-1)/2) *
      dev_stress;

  // Update stress component
  Eigen::Matrix<double, 6, 1> pstress = stress;
  pstress += this->de_ * (dstrain-viscous_strain_rate);
  return pstress;
}