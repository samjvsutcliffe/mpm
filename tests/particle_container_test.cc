// No-debase test
#include <functional>
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "container.h"
#include "particle.h"

//! \brief Check particle container class for 2D case
TEST_CASE("Particle container is checked for 2D case",
          "[particlecontainer][2D]") {
  // Dimension
  const unsigned Dim = 2;
  // Tolerance
  const double Tolerance = 1.E-7;

  // Particle 1
  mpm::Index id1 = 0;
  Eigen::Vector2d coords;
  coords.setZero();
  auto particle1 = std::make_shared<mpm::Particle<Dim>>(id1, coords);

  // Particle 2
  mpm::Index id2 = 1;
  auto particle2 = std::make_shared<mpm::Particle<Dim>>(id2, coords);

  // Particle container
  auto particlecontainer =
      std::make_shared<mpm::Container<mpm::Particle<Dim>>>();

  // Check insert particle
  SECTION("Check insert particle functionality") {
    // Insert particle 1 and check status
    bool status1 = particlecontainer->insert(particle1);
    REQUIRE(status1 == true);
    // Insert particle 2 and check status
    bool status2 = particlecontainer->insert(particle2);
    REQUIRE(status2 == true);
    // Check size of particle hanlder
    REQUIRE(particlecontainer->size() == 2);
    // Clear particle container
    particlecontainer->clear();
    // Check size of particle hanlder
    REQUIRE(particlecontainer->size() == 0);
  }

  // Check iterator
  SECTION("Check particle range iterator") {
    // Insert particle 1
    particlecontainer->insert(particle1);
    // Insert particle 2
    particlecontainer->insert(particle2);
    // Check size of particle hanlder
    std::size_t counter = 0;
    for (auto itr = particlecontainer->begin(); itr != particlecontainer->end();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each particle is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));
      ++counter;
    }

    // Iterate over particles and check if the number of particles is good
    REQUIRE(counter == 2);
  }

  // Check for_each
  SECTION("Check particle for_each") {
    // Insert particle 1
    particlecontainer->insert(particle1);
    // Insert particle 2
    particlecontainer->insert(particle2);
    // Check size of particle hanlder
    REQUIRE(particlecontainer->size() == 2);

    // Check coordinates before updating
    for (auto itr = particlecontainer->begin(); itr != particlecontainer->end();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each particle is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));
    }

    // Set coordinates to unity
    coords << 1., 1.;

    // Iterate through particle container to update coordinaates
    particlecontainer->for_each(  // function structure
        std::bind(static_cast<void (mpm::Particle<Dim>::*)(
                      const Eigen::Matrix<double, Dim, 1>&)>(
                      // function
                      &mpm::Particle<Dim>::coordinates),
                  // arguments
                  std::placeholders::_1, coords));

    // Check if update has gone through
    for (auto itr = particlecontainer->begin(); itr != particlecontainer->end();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each particle is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(1.).epsilon(Tolerance));
    }
  }
}

//! \brief Check particle container class for 3D case
TEST_CASE("Particle container is checked for 3D case",
          "[particlecontainer][3D]") {
  // Dimension
  const unsigned Dim = 3;
  // Tolerance
  const double Tolerance = 1.E-7;

  // Particle 1
  mpm::Index id1 = 0;
  Eigen::Vector3d coords;
  coords.setZero();
  auto particle1 = std::make_shared<mpm::Particle<Dim>>(id1, coords);

  // Particle 2
  mpm::Index id2 = 1;
  auto particle2 = std::make_shared<mpm::Particle<Dim>>(id2, coords);

  // Particle container
  auto particlecontainer =
      std::make_shared<mpm::Container<mpm::Particle<Dim>>>();

  // Check insert particle
  SECTION("Check insert particle functionality") {
    // Insert particle 1 and check status
    bool status1 = particlecontainer->insert(particle1);
    REQUIRE(status1 == true);
    // Insert particle 2 and check status
    bool status2 = particlecontainer->insert(particle2);
    REQUIRE(status2 == true);
    // Check size of particle hanlder
    REQUIRE(particlecontainer->size() == 2);
    // Clear particle container
    particlecontainer->clear();
    // Check size of particle hanlder
    REQUIRE(particlecontainer->size() == 0);
  }

  // Check iterator
  SECTION("Check particle range iterator") {
    // Insert particle 1
    particlecontainer->insert(particle1);
    // Insert particle 2
    particlecontainer->insert(particle2);
    // Check size of particle hanlder
    std::size_t counter = 0;
    for (auto itr = particlecontainer->begin(); itr != particlecontainer->end();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each particle is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));

      ++counter;
    }
    // Iterate over particles and check if the number of particles is good
    REQUIRE(counter == 2);
  }

  // Check for_each
  SECTION("Check particle for_each") {
    // Insert particle 1
    particlecontainer->insert(particle1);
    // Insert particle 2
    particlecontainer->insert(particle2);
    // Check size of particle hanlder
    REQUIRE(particlecontainer->size() == 2);

    for (auto itr = particlecontainer->begin(); itr != particlecontainer->end();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each particle is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));
    }

    // Set coordinates to unity
    coords << 1., 1., 1.;

    // Iterate through particle container to update coordinaates
    particlecontainer->for_each(
        // function structure
        std::bind(static_cast<void (mpm::Particle<Dim>::*)(
                      const Eigen::Matrix<double, Dim, 1>&)>(
                      // function
                      &mpm::Particle<Dim>::coordinates),
                  // arguments
                  std::placeholders::_1, coords));

    // Check if update has gone through
    for (auto itr = particlecontainer->begin(); itr != particlecontainer->end();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each particle is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(1.).epsilon(Tolerance));
    }
  }
}