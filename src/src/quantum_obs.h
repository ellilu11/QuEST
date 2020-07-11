#ifndef QUANTUM_OBS_H
#define QUANTUM_OBS_H

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "common.h"
// #include "quantum_dot.h"

class ObserverDot;

typedef std::vector<ObserverDot> ObsVector;

class ObserverDot {
 public:
  ObserverDot() = default;
  ObserverDot(const Eigen::Vector3d &pos);

  const Eigen::Vector3d &position() const { return pos; }

  friend class QuantumDot;
  // friend Eigen::Vector3d separation_srcobs(const QuantumDot &, const ObserverDot &);
  friend std::ostream &operator<<(std::ostream &, const ObserverDot &);
  friend std::istream &operator>>(std::istream &, ObserverDot &);

 private:
  Eigen::Vector3d pos;
};

ObsVector import_obss(const std::string &);

#endif
