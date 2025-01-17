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
  this->bulk_modulus_ = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));
  const double G = youngs_modulus_ / (2.0 * (1. + poisson_ratio_));

  const double a1 = bulk_modulus_ + (4.0 / 3.0) * G;
  const double a2 = bulk_modulus_ - (2.0 / 3.0) * G;
  de_ = Matrix6x6::Zero();
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
    const Eigen::Array<double,6,1> second_invar_mult = (Eigen::Array<double, 6, 1>() << 1, 1, 1, 0.5, 0.5, 0.5).finished();
    Eigen::Matrix<double, 6, 1> dev_strain = strain;
    for (int i = 0; i < 3; ++i) {
      dev_strain(i) -= trace_strain;
    }
    double effective_strain =
        0.5 * (dev_strain.dot((dev_strain.array() * second_invar_mult).matrix()));
    return 0.5 * visc_factor *
           std::pow(effective_strain + 1e-30, 0.5 * ((1.0 / visc_power) - 1.0));
}
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::Glen<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain, const ParticleBase<Tdim>* ptr,
    mpm::dense_map* state_vars) {

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
  Eigen::Matrix<double, 6, 1> dstress = this->de_ * dstrain;
  if (viscosity > 0.0) {
	  const double relaxation_constant = youngs_modulus_ * ((*state_vars).at("dt")) / (2.0 * (1.0 - poisson_ratio_) * viscosity);
	  dstress -= dev_stress * relaxation_constant;
  }
  return stress + ptr->objective_stress_increment(dstress,stress);
}
