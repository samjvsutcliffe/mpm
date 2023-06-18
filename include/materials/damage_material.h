#ifndef MPM_MATERIAL_DAMAGE_H_
#define MPM_MATERIAL_DAMAGE_H_

#include <limits>

#include "Eigen/Dense"

#include "material.h"
class Damage : public Material<Tdim> {
  double critical_stress_{0};
  double damage_rate_{0};
  double local_length_{0};
  double critical_damage_{0};
}