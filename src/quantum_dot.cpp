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
                                          double damping_rr,
                                          const bool rotating) const
{
/*  const double T1 = 1000 * damping.first;
  const double T2 = 1000 * damping.second;

  const cmplx m0 = 2.0 * damping_rr * pow(laser_freq,3) * ( pow( real(rho[1]), 2) + pow( imag(rho[1]), 2) ) + 
                    2.0 * ( real(rho[1]) * imag(rabi) - imag(rho[1]) * real(rabi) ) + (1.0 - rho[0] ) / T1;
  
  const cmplx m1 = -iu * ( ( iu * damping_rr * pow(laser_freq,3) * rho[1] + rabi ) * ( 1.0 - 2.0 * rho[0] ) ) - rho[1] / T2;
*/

/* temporal deriv test

  const cmplx m0 = rabi;
  const cmplx m1 = 0.0;
*/

  const double time = time_idx * 5.0e-5;

  const cmplx m0 = -iu * (rabi * conj(rho[1]) - conj(rabi) * rho[1] ) -
                   0.0 * (rho[0] - 1.0) / (damping.first); // same in rot frame?
                   // - damping_rr * ( pow(laser_freq,3) / 2.0 * (pow( 1.0 - 2.0 * rho[0], 2 ) - 1.0);

  const cmplx m1 = rotating ?  
      -iu * (rabi * (1.0 - 2.0 * rho[0]) + rho[1] * (laser_freq - freq)) -
      0.0 * rho[1] / (damping.second) :
     // + damping_rr * pow(laser_freq,3) * rho[1] * ( 1.0 - 2.0 * rho[0] ) :
      -iu * (rabi * (1.0 - 2.0 * rho[0]) - rho[1] * freq) -
      0.0 * rho[1] / (damping.second);
     // + damping_rr * pow(laser_freq,3) * rho[1] * ( 1.0 - 2.0 * rho[0] ) ;

  return matrix_elements(m0, m1);
}

Eigen::Vector3d separation(const QuantumDot &d1, const QuantumDot &d2)
{
  return d2.pos - d1.pos;
}

/*const Eigen::Vector3d QuantumDot::separation_srcobs(const QuantumDot &d1, const ObserverDot &d2) 
{
  return d2.pos - d1.pos;
}*/

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
                                             double damping_rr,
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
