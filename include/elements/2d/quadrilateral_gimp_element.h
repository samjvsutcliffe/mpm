#ifndef MPM_GIMP_ELEMENT_H_
#define MPM_GIMP_ELEMENT_H_

#include "quadrilateral_element.h"

namespace mpm {

//! Quadrilateral GIMP element class derived from Quadrilateral
//! \brief Quadrilateral GIMP element
//! \details 16-noded quadrilateral GIMP element \n
//! Shape function, gradient shape function, B-matrix, indices \n
//! <pre>
//!
//!   13----------12----------11----------10
//!   |           |           |           |
//!   |           |           |           |
//!   |           |           |           |
//!   |        (-1, 1)      (1,1)         |
//!   14----------3-----------2-----------9
//!   |           |           |           |
//!   |           | particle  |           |
//!   |           | location  |           |
//!   |           |           |           |
//!   15----------0-----------1-----------8
//!   |        (-1,-1)      (1,-1)        |
//!   |           |           |           |
//!   |           |           |           |
//!   |           |           |           |
//!   4-----------5-----------6-----------7
//!
//! </pre>
//!
//!
//! \tparam Tdim Dimension
//! \tparam Tnfunctions Number of functions
template <unsigned Tdim>
class QuadrilateralGIMPElement : public QuadrilateralElement<2, 4> {

 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! constructor with number of shape functions
  QuadrilateralGIMPElement() : QuadrilateralElement<2, 4>() {
    static_assert(Tdim == 2, "Invalid dimension for a GIMP element");
    //static_assert((Tnfunctions >= 4) && (Tnfunctions <= 16),
    //              "Specified number of shape functions is not defined");

    //! Logger
    std::string logger = "quadrilateral_gimp::<" + std::to_string(Tdim) + ", " +
                         std::to_string(Tnfunctions) + ">";
    console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
  }

  //! Evaluate shape functions at given local coordinates
  //! \param[in] xi given local coordinates
  //! \param[in] particle_size Particle size
  //! \param[in] deformation_gradient Deformation gradient
  //! \retval shapefn Shape function of a given cell
  Eigen::VectorXd shapefn(const VectorDim& xi, const VectorDim& particle_size,
                          const VectorDim& deformation_gradient) const override;

  //! Evaluate local shape functions at given local coordinates
  //! \param[in] xi given local coordinates
  //! \param[in] particle_size Particle size
  //! \param[in] deformation_gradient Deformation gradient
  //! \retval shapefn Shape function of a given cell
  Eigen::VectorXd shapefn_local(
      const VectorDim& xi, const VectorDim& particle_size,
      const VectorDim& deformation_gradient) const override;

  //! Evaluate gradient of shape functions
  //! \param[in] xi given local coordinates
  //! \param[in] particle_size Particle size
  //! \param[in] deformation_gradient Deformation gradient
  //! \retval grad_shapefn Gradient of shape function of a given cell
  Eigen::MatrixXd grad_shapefn(
      const VectorDim& xi, const VectorDim& particle_size,
      const VectorDim& deformation_gradient) const override;

  //! Compute Jacobian
  //! \param[in] xi given local coordinates
  //! \param[in] nodal_coordinates Coordinates of nodes forming the cell
  //! \param[in] particle_size Particle size
  //! \param[in] deformation_gradient Deformation gradient
  //! \retval jacobian Jacobian matrix
  Eigen::Matrix<double, Tdim, Tdim> jacobian(
      const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
      const VectorDim& particle_size,
      const VectorDim& deformation_gradient) const override;

  //! Compute Jacobian local
  //! \param[in] xi given local coordinates
  //! \param[in] nodal_coordinates Coordinates of nodes forming the cell
  //! \param[in] particle_size Particle size
  //! \param[in] deformation_gradient Deformation gradient
  //! \retval jacobian Jacobian matrix
  Eigen::Matrix<double, Tdim, Tdim> jacobian_local(
      const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
      const VectorDim& particle_size,
      const VectorDim& deformation_gradient) const override;

  //! Evaluate the B matrix at given local coordinates for a real cell
  //! \param[in] xi given local coordinates
  //! \param[in] nodal_coordinates Coordinates of nodes forming the cell
  //! \param[in] particle_size Particle size
  //! \param[in] deformation_gradient Deformation gradient
  //! \retval bmatrix B matrix
  std::vector<Eigen::MatrixXd> bmatrix(
      const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
      const VectorDim& particle_size,
      const VectorDim& deformation_gradient) const override;

  //! Return the type of shape function
  mpm::ShapefnType shapefn_type() const override {
    return mpm::ShapefnType::GIMP;
  }

  //! Return number of shape functions
  unsigned nfunctions() const override { return Tnfunctions; }

  //! Return if natural coordinates can be evaluates
  bool isvalid_natural_coordinates_analytical() const override { return false; }

  //! Compute Natural coordinates of a point (analytical)
  //! \param[in] nodal_coordinates Coordinates of nodes forming the cell
  //! \param[in] point Location of the point in cell
  //! \retval xi Return the local coordinates
  VectorDim natural_coordinates_analytical(
      const VectorDim& point,
      const Eigen::MatrixXd& nodal_coordinates) const override;

  //! Compute local node matrix
  virtual void InitialiseLocalMapping(std::vector<int>& local_node_ids) override{
    //console_->info("Construct local nodes");
   Eigen::Matrix<double, 16, Tdim> local_nodes =
  (Eigen::Matrix<double, 16, Tdim>() << 
    -1, -1,
    1 ,-1 ,
    1 ,1  ,
    -1, 1 ,

    -3, -3,
    -1 ,-3,
    1 ,-3,
    3 ,-3,

    3 ,-1,
    3 ,1,
    3 ,3,

    1 ,3,
    -1 ,3,
    -3 ,3,

    -3 ,1,
    -3 ,-1


    //1 ,-1 ,
    //1 ,1  ,
    //-1, 1 ,
    //-3, -1,
    //3 ,-1 ,
    //1 ,-3 ,
    //1 ,3  ,
    //3 ,1  ,
    //-3, 1 ,
    //-1, 3 ,
    //-1, -3,
    //-3, -3,
    //3 ,-3 ,
    //3 ,3  ,
    //-3, 3
      ).finished();
    Tnfunctions = local_node_ids.size();
    //console_->info("TN functions {}",Tnfunctions);
    local_node_mapping_ =
        Eigen::MatrixXd::Zero(Tnfunctions,Tdim);
   for (int i = 0; i < Tnfunctions; ++i)
   {
      local_node_mapping_.row(i) = local_nodes.row(local_node_ids[i]);
   }
   //console_->info("Mapping complete");
  };

 private:
  //! Return natural nodal coordinates
  Eigen::MatrixXd natural_nodal_coordinates() const;

  //! Mapping matrix of local nodes
  Eigen::MatrixXd local_node_mapping_;
  //! How many shape functions exist
  int Tnfunctions = 0;

  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};

}  // namespace mpm
#include "quadrilateral_gimp_element.tcc"

#endif  // MPM_QUADRILATERAL_ELEMENT_H_
