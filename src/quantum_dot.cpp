#include "quantum_dot.h"

using namespace std;

QuantumDot::QuantumDot(const Eigen::Vector3d &pos,
                       const double freq,
                       const std::pair<double, double> &damp,
                       const Eigen::Vector3d &dipr,
                       const Eigen::Vector3d &dipi)
    : pos(pos), freq(freq), damp(damp), dipr(dipr), dipi(dipi)
{
}

Eigen::Vector3d separation(const QuantumDot &d1, const QuantumDot &d2)
{
  return d2.pos - d1.pos;
}

int max_transit_steps_between_dots(const std::shared_ptr<DotVector> dots,
                                   const double c,
                                   const double dt)
{
  double max_separation{0};
  for(auto dot_i_it = (*dots).begin(); dot_i_it != (*dots).end(); ++dot_i_it) {
    for(auto dot_j_it = dot_i_it + 1; dot_j_it != (*dots).end(); ++dot_j_it) {
      max_separation =
          std::max(max_separation, separation(*dot_i_it, *dot_j_it).norm());
    }
  }
  return std::ceil(max_separation / (c * dt));
}

std::ostream &operator<<(std::ostream &os, const QuantumDot &qd)
{
  os << qd.pos.transpose() << " " << qd.freq << " " << qd.damp.first << " "
     << qd.damp.second << " " << qd.dipr.transpose() << std::endl; // << " " << qd.dipi.transpose();
  return os;
}

std::istream &operator>>(std::istream &is, QuantumDot &qd)
{
  is >> qd.pos[0] >> qd.pos[1] >> qd.pos[2] >> qd.freq >> qd.damp.first >>
      qd.damp.second >> qd.dipr[0] >> qd.dipr[1] >> qd.dipr[2];
        // >> qd.dipi[0] >> qd.dipi[1] >> qd.dipi[2];
  return is;
}

DotVector import_dots(const std::string &fname)
{
  std::ifstream ifs(fname);
  if(!ifs) throw std::runtime_error("Could not open " + fname);

  std::istream_iterator<QuantumDot> in_iter(ifs), eof;
  return DotVector(in_iter, eof);
}


void set_dipole_of_dots(std::shared_ptr<DotVector> dots, const Eigen::Vector3d dip)
{
  for (int dot = 0; dot < (*dots).size(); ++dot)
    (*dots)[dot].set_dipole(dip);
}

