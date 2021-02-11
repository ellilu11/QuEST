#include "pulse.h"

Pulse::Pulse(const double amplitude, const double delay,
             const double width, const double freq, const double c0,
             const Eigen::Vector3d &wavevector,
             const Eigen::Vector3d &polarization)
    : amplitude(amplitude),
      delay(delay),
      width(width),
      freq(freq),
      c0(c0),
      wavevector(wavevector.normalized()),
      polarization(polarization.normalized())
{

}

Eigen::Vector3cd Pulse::operator()(const Eigen::Vector3d &r,
                                   const double t, const bool rotating, const bool rwa) const
{
  const double wavemag = freq/c0;
  const double arg = wavemag*wavevector.dot(r) - freq * (t - delay);
  const double width_factor = 1.0;

  Eigen::Vector3cd amp_vector = 
      amplitude/width_factor * polarization * gaussian(arg / ( width*width_factor ) );
  return amp_vector * (rotating
      ? (rwa ? 0.5 : cos(arg) * exp(iu*freq*t))
      : cos(arg) );
}


std::ostream &operator<<(std::ostream &os, const Pulse &p)
{
  os << p.amplitude << " " << p.delay << " " << p.width << " " << p.freq << " " << p.c0 << " "
     << p.wavevector.transpose() << " " << p.polarization.transpose();
  return os;
}

std::istream &operator>>(std::istream &is, Pulse &p)
{
  is >> p.amplitude >> p.delay >> p.width >> p.freq >> p.c0 >>
      p.wavevector[0] >> p.wavevector[1] >> p.wavevector[2] >> 
      p.polarization[0] >> p.polarization[1] >> p.polarization[2];
  return is;
}

Pulse read_pulse_config(const std::string &filename)
{
  std::ifstream ifs(filename);
  if(!ifs) throw std::runtime_error("Could not open " + filename);

  Pulse p;
  ifs >> p;

  return p;
}
