#include "quantum_obs.h"

ObserverDot::ObserverDot(const Eigen::Vector3d &pos)
    : pos(pos)
{
}

std::ostream &operator<<(std::ostream &os, const ObserverDot &qd)
{
  os << qd.pos.transpose();
  return os;
}

std::istream &operator>>(std::istream &is, ObserverDot &qd)
{
  is >> qd.pos[0] >> qd.pos[1] >> qd.pos[2];
  return is;
}

ObsVector import_obss(const std::string &fname)
{
  std::ifstream ifs(fname);
  if(!ifs) throw std::runtime_error("Could not open " + fname);

  std::istream_iterator<ObserverDot> in_iter(ifs), eof;
  return ObsVector(in_iter, eof);
}

