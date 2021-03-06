#ifndef PULSE_H
#define PULSE_H

#include <Eigen/Dense>
#include <fstream>
#include <iostream>

#include "math_utils.h"
#include "common.h"

class Pulse {
 public:
  Pulse() = default;
  Pulse(const double, const double, const double, const double, const double,
        const Eigen::Vector3d &, const Eigen::Vector3d &);

  Eigen::Vector3cd operator()(const Eigen::Vector3d &, const double, const bool, const bool) const;
 
  const Eigen::Vector3d wavevec() const { return freq/c0 * wavevector.normalized(); }

  friend std::ostream &operator<<(std::ostream &, const Pulse &);
  friend std::istream &operator>>(std::istream &, Pulse &);

 private:
  double amplitude, delay, width, freq, c0;
  Eigen::Vector3d wavevector, polarization;
};

Pulse read_pulse_config(const std::string &);

#endif
