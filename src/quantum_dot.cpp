#include "quantum_dot.h"

using namespace std;

QuantumDot::QuantumDot(const Eigen::Vector3d &pos,
                       const double freq,
                       const std::pair<double, double> &damping,
                       const Eigen::Vector3d &dipr,
                       const Eigen::Vector3d &dipi)
    : pos(pos), freq(freq), damping(damping), dipr(dipr), dipi(dipi)
{
}

matrix_elements QuantumDot::liouville_rhs(const matrix_elements &rho,
                                          const cmplx rabi,
                                          const int time_idx,
                                          const double laser_freq,
                                          const double damping_rr,
                                          const bool rotating) const
{
  const cmplx m0 = -iu * (rabi * std::conj(rho[1]) - std::conj(rabi) * rho[1]) -
                   0.0 * (rho[0] - 1.0) / damping.first;

  cmplx m1_temp = -iu * (rabi * (1.0 - 2.0 * rho[0])) - 0.0 * rho[1] / damping.second;

  m1_temp -=
      rotating ? iu * rho[1] * (laser_freq - freq) : iu * rho[1] * (-freq);

  const cmplx m1 = m1_temp;
 
/*  const cmplx m0 = -iu * (rabi * conj(rho[1]) - conj(rabi) * rho[1] ) -
                   (rho[0] - 1.0) / (damping.first);

  cmplx m1_temp = 
    -iu * (rabi * (1.0 - 2.0 * rho[0])) - rho[1] / (damping.second);
    
  m1_temp += rotating ? -iu * rho[1] * (laser_freq - freq) : iu * rho[1] * freq;

  const cmplx m1 = m1_temp;
*/
  return matrix_elements(m0, m1);
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
  os << qd.pos.transpose() << " " << qd.freq << " " << qd.damping.first << " "
     << qd.damping.second << " " << qd.dipr.transpose(); // << " " << qd.dipi.transpose();
  return os;
}

std::istream &operator>>(std::istream &is, QuantumDot &qd)
{
  is >> qd.pos[0] >> qd.pos[1] >> qd.pos[2] >> qd.freq >> qd.damping.first >>
      qd.damping.second >> qd.dipr[0] >> qd.dipr[1] >> qd.dipr[2];
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


void set_dipolevec(std::shared_ptr<DotVector> dots, const Eigen::Vector3d dip)
{
  for (int dot = 0; dot < (*dots).size(); ++dot)
    (*dots)[dot].set_dipole(dip);
}

std::vector<BlochFunctionType> rhs_functions(const DotVector &dots,
                                             const double laser_frequency,
                                             const double damping_rr,
                                             const bool rotating)
{
  std::vector<BlochFunctionType> funcs(dots.size());

  using std::placeholders::_1;
  using std::placeholders::_2;
  using std::placeholders::_3;
  std::transform(dots.begin(), dots.end(), funcs.begin(),
                 [laser_frequency,damping_rr,rotating](const QuantumDot &d) {
                   return std::bind(&QuantumDot::liouville_rhs, d, _1, _2, _3,
                                    laser_frequency, damping_rr, rotating);
                 });
  return funcs;
}
