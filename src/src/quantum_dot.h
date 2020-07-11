#ifndef QUANTUM_DOT_H
#define QUANTUM_DOT_H

#include <Eigen/Dense>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "common.h"
#include "quantum_obs.h"

class QuantumDot;

typedef std::vector<QuantumDot> DotVector;
typedef std::pair<DotVector::iterator, DotVector::iterator> DotRange;
typedef std::pair<DotVector::const_iterator, DotVector::const_iterator>
    const_DotRange;
typedef Eigen::Vector2cd matrix_elements;
typedef std::function<Eigen::Vector2cd(const Eigen::Vector2cd,
                                       const std::complex<double>)>
    BlochFunctionType;
enum MatrixElement { RHO_00, RHO_01 };

class QuantumDot {
 public:
  QuantumDot() = default;
  QuantumDot(const Eigen::Vector3d &pos)
      : QuantumDot(pos, 0, {0.0, 0.0}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero()){};
  QuantumDot(const Eigen::Vector3d &pos, const Eigen::Vector3d &dipr, const Eigen::Vector3d &dipi)
      : QuantumDot(pos, 0, {0.0, 0.0}, dipr, dipi){};
  QuantumDot(const Eigen::Vector3d &,
             const double,
             const std::pair<double, double> &,
             const Eigen::Vector3d &,
             const Eigen::Vector3d &);

  matrix_elements liouville_rhs(const matrix_elements &,
                                const cmplx,
                                const double, 
                                double,
                                const bool) const;

  friend class ObserverDot;

  const Eigen::Vector3d &position() const { return pos; }
  const Eigen::Vector3d &dipole() const { return dipr; }
  const Eigen::Vector3d &dipole_imag() const { return dipi; }
  void set_dipole(const Eigen::Vector3d dip) { dipr = dip; }

  friend Eigen::Vector3d separation(const QuantumDot &, const QuantumDot &);
  // const Eigen::Vector3d separation_srcobs(const QuantumDot &, const ObserverDot &);
/*  friend inline std::complex<double> dyadic_product(const QuantumDot &obs,
                                      const Eigen::Matrix3d &dyad,
                                      const QuantumDot &src)
  {
    return obs.dipr.transpose() * dyad * src.dipr + obs.dipi.transpose() * dyad * src.dipi +
     iu * (obs.dipr.transpose() * dyad * src.dipi - obs.dipi.transpose() * dyad * src.dipr );
  }*/

  friend std::ostream &operator<<(std::ostream &, const QuantumDot &);
  friend std::istream &operator>>(std::istream &, QuantumDot &);

  Eigen::Vector3d pos;
  Eigen::Vector3d dipr;
 
 private:
  double freq;
  std::pair<double, double> damping;
  Eigen::Vector3d dipi;
};

DotVector import_dots(const std::string &);
void set_dipolevec(std::shared_ptr<DotVector>, const Eigen::Vector3d dip);
std::vector<BlochFunctionType> rhs_functions(const DotVector &, const double, double, const bool);

#endif
