//! Constructor with id and coordinates
template <unsigned Tdim>
mpm::ParticleBase<Tdim>::ParticleBase(Index id, const VectorDim& coord)
    : id_{id} {
  // Check if the dimension is between 1 & 3
  static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");
  coordinates_ = coord;
  status_ = true;
}

//! Constructor with id, coordinates and status
template <unsigned Tdim>
mpm::ParticleBase<Tdim>::ParticleBase(Index id, const VectorDim& coord,
                                      bool status)
    : mpm::ParticleBase<Tdim>::ParticleBase(id, coord) {
  status_ = status;
}

template <unsigned Tdim>
Eigen::Matrix<double,6,1> mpm::ParticleBase<Tdim>::matrix_to_voigt(Eigen::Matrix<double,3,3> mat) {
  return (Eigen::Matrix<double,6,1>() <<
          mat(0,0), mat(1,1),mat(2,2),
          mat(0,1), mat(1,2),mat(0,2)
          ).finished();
}

template <unsigned Tdim>
Eigen::Matrix<double,3,3> mpm::ParticleBase<Tdim>::voigt_to_matrix(Eigen::Matrix<double,6,1> voigt) {
  return (Eigen::Matrix3d() <<
          voigt(0), voigt(3), voigt(5),
          voigt(3), voigt(1), voigt(4),
          voigt(5), voigt(4), voigt(2)).finished();
}
template <unsigned Tdim>
Eigen::Matrix<double,3,3> mpm::ParticleBase<Tdim>::vorticity_matrix(Eigen::Matrix<double,6,1> voigt) {
  return (Eigen::Matrix3d() <<
      voigt(0), voigt(3), voigt(5),
      -voigt(3), voigt(1), voigt(4),
      -voigt(5), -voigt(4), voigt(2)).finished();
};
