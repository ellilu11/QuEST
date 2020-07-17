#include "pulse.h"

Pulse::Pulse(const double amplitude, const double delay, const double width,
             const double freq, const Eigen::Vector3d &wavevector,
             const Eigen::Vector3d &polarization)
    : amplitude(amplitude),
      delay(delay),
      width(width),
      freq(freq),
      wavevector(wavevector),
      polarization(polarization.normalized())
{

}

Eigen::Vector3cd Pulse::operator()(const Eigen::Vector3d &r,
                                  const double t, const int deriv, const bool rotating) const
{
  const double arg = wavevector.dot(r) - freq * (t - delay);
  return amplitude * polarization * gaussian(arg / width); 

/*    const double arg = t - delay;
    const double arg2 = wavevector.dot(r) - freq*arg;
    const double a = pow(freq,2) / ( 2.0 * pow(width,2) );
    const double b = -wavevector.dot(r) * freq / pow(width,2);
    const double c = pow(wavevector.dot(r),2) / ( 2.0 * pow(width,2) );
    return (amplitude * polarization) * gaussian(arg, a, b, c, deriv) * 
      (rotating ? 
        // ( ( 1.0 + 0.0 * exp(2.0*iu*arg2) ) / 2.0 )
         cos(arg2) * exp( iu*freq*t )
       : cos(arg2));
*/
}


std::ostream &operator<<(std::ostream &os, const Pulse &p)
{
  os << p.amplitude << " " << p.delay << " " << p.width << " " << p.freq << " "
     << p.wavevector.transpose() << " " << p.polarization.transpose();
  return os;
}

std::istream &operator>>(std::istream &is, Pulse &p)
{
  is >> p.amplitude >> p.delay >> p.width >> p.freq >> p.wavevector[0] >>
      p.wavevector[1] >> p.wavevector[2] >> p.polarization[0] >>
      p.polarization[1] >> p.polarization[2];
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
