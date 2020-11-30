#ifndef GREEN_FUNCTION_H
#define GREEN_FUNCTION_H

#include <Eigen/Dense>
#include <array>
#include <type_traits>
#include <vector>

#include "common.h"
#include "lagrange_set.h"

namespace Propagation {
  template <class T>
  using Mat3D = Eigen::Matrix<T, 3, 3>;

  template <class T>
  class Kernel;

  template <class T>
  class Identity;

  template <class T>
  class Laplace;

  class Helmholtz;

  template <class T>
  class DelSq_Laplace;

  template <class T>
  class EFIE;

  class RotatingEFIE;

  template <class T>
  class MFIE;

  class RotatingMFIE;

  class SelfEFIE;

  class SelfRotatingEFIE;
}

template <class T>
class Propagation::Kernel {
 public:
  static_assert(std::is_same<T, double>::value || std::is_same<T, cmplx>::value,
                "Propagation kernels require numeric types");
  virtual const std::vector<Mat3D<T>> &coefficients(
      const Eigen::Vector3d &, const Interpolation::UniformLagrangeSet &) = 0;

 protected:
  std::vector<Mat3D<T>> coefs_;
};


template <class T>
class Propagation::Identity : public Propagation::Kernel<T> {
 public:
  Identity() = default;
  const std::vector<Mat3D<T>> &coefficients(
      __attribute__((unused)) const Eigen::Vector3d &dr,
      const Interpolation::UniformLagrangeSet &interp) // construct interpolation order
  {
    this->coefs_.resize(interp.order() + 1);

    for(int i = 0; i <= interp.order(); ++i) {
      this->coefs_[i] = Mat3D<T>::Identity() * interp.evaluations[0][i];
    }

    return this->coefs_;
  }
};


template <class T>
class Propagation::Laplace : public Propagation::Kernel<T> {
 public:
  explicit Laplace(const double k2 = 1) : k2(k2) {}
  const std::vector<Mat3D<T>> &coefficients(
      const Eigen::Vector3d &dr,
      const Interpolation::UniformLagrangeSet &interp)
  {
    this->coefs_.resize(interp.order() + 1);

    for(int i = 0; i <= interp.order(); ++i) {
      this->coefs_[i] =
          Mat3D<T>::Identity() * interp.evaluations[0][i] * k2 / dr.norm();
    }

    return this->coefs_;
  }

 private:
  double k2;
};

class Propagation::Helmholtz : public Propagation::Kernel<cmplx> {
 public:
  explicit Helmholtz(const double k2, const double k) : k2(k2), k(k) {}
  const std::vector<Mat3D<cmplx>> &coefficients(
      const Eigen::Vector3d &dr,
      const Interpolation::UniformLagrangeSet &interp)
  {
    this->coefs_.resize(interp.order() + 1);

    for(int i = 0; i <= interp.order(); ++i) {
      this->coefs_[i] = Mat3D<cmplx>::Identity() * interp.evaluations[0][i] *
                        k2 * std::exp(-iu * k * dr.norm()) / dr.norm();
    }

    return this->coefs_;
  }

 private:
  double k2, k;
};

template <class T>
class Propagation::DelSq_Laplace : public Propagation::Kernel<T> {
 public:
  explicit DelSq_Laplace(const double c, const double k2 = 1) : c_{c}, k2_{k2}
  {
  }
  const std::vector<Mat3D<T>> &coefficients(
      const Eigen::Vector3d &dr,
      const Interpolation::UniformLagrangeSet &interp)
  {
    this->coefs_.resize(interp.order() + 1);

    const Eigen::Matrix3d rh = dr * dr.transpose() / dr.squaredNorm();
    const double R = dr.norm();
    const double R_sq = dr.squaredNorm();
    const double R_cu = std::pow(dr.norm(), 3);

    for(int i = 0; i <= interp.order(); ++i) {
      this->coefs_[i] =
          k2_ * ((3 * interp.evaluations[0][i] / R_cu +
                  3 * interp.evaluations[1][i] / (c_ * R_sq) +
                  interp.evaluations[2][i] / (std::pow(c_, 2) * R)) *
                     rh +
                 (-1 * interp.evaluations[0][i] / R_cu -
                  interp.evaluations[1][i] / (c_ * R_sq)) *
                     Eigen::Matrix3d::Identity());
    }

    return this->coefs_;
  }

 private:
  double c_, k2_;
};

template <class T>
class Propagation::EFIE : public Propagation::Kernel<T> {
 public:
  EFIE(const double c, const double k2, const double beta, const double dist0) : c_(c), k2_(k2), beta_(beta), dist0_(dist0){};
  const std::vector<Mat3D<T>> &coefficients(
      const Eigen::Vector3d &dr,
      const Interpolation::UniformLagrangeSet &interp)
  {
    this->coefs_.resize(interp.order() + 1);

    const auto dyads(spatial_dyads(dr));
    for(int i = 0; i <= interp.order(); ++i) {
      this->coefs_[i] = Eigen::Matrix3d::Zero();

      if ( dr.norm() > 0.0 ) {
            this->coefs_[i] = -k2_ * (dyads[0] * interp.evaluations[0][i] +
                                      dyads[1] * interp.evaluations[1][i] +
                                      dyads[2] * interp.evaluations[2][i] );

      } /*else if ( dr.norm() == 0.0 ){
              this->coefs_[i] = -k2_ * 
                    (dyads[0] * interp.evaluations[0][i] +
                     dyads[3] * interp.evaluations[2][i]);
            
          this->coefs_[i] = 
                -beta_ / pow( 5.2917721e-4, 2 ) * Eigen::Matrix3d::Identity() * 
                interp.evaluations[3][i] ;
        }*/

    }

    return this->coefs_;
  }

 protected:
  double c_, k2_, beta_, dist0_;

  std::array<Eigen::Matrix3d, 4> spatial_dyads(const Eigen::Vector3d &dr) const
  {
    std::array<Eigen::Matrix3d, 4> results = {
        {identity_minus_3rsq(dr) * std::pow(c_, 2) / std::pow(dr.norm(), 3),
         identity_minus_3rsq(dr) * c_ / dr.squaredNorm(),
         identity_minus_rsq(dr) / dr.norm(),
         identity_plus_rsq(dr) / ( 2.0 * dr.norm() )}};
    return results;
  }

  static Eigen::Matrix3d rhat_dyadic(const Eigen::Vector3d &dr)
  {
    return dr * dr.transpose() / dr.squaredNorm();
  }

  static Eigen::Matrix3d identity_minus_rsq(const Eigen::Vector3d &dr)
  {
    return Eigen::Matrix3d::Identity() - rhat_dyadic(dr);
  }

  static Eigen::Matrix3d identity_minus_3rsq(const Eigen::Vector3d &dr)
  {
    return Eigen::Matrix3d::Identity() - 3 * rhat_dyadic(dr);
  }

  static Eigen::Matrix3d identity_plus_rsq(const Eigen::Vector3d &dr)
  {
    return Eigen::Matrix3d::Identity() + rhat_dyadic(dr);
  }
};

class Propagation::RotatingEFIE : public Propagation::EFIE<cmplx> {
 public:
  RotatingEFIE(const double c, const double k2, const double omega, const double beta, const double dist0)
      : EFIE<cmplx>(c, k2, beta, dist0), omega_(omega){};

  const std::vector<Eigen::Matrix3cd> &coefficients(
      const Eigen::Vector3d &dr,
      const Interpolation::UniformLagrangeSet &interp)
  {
    this->coefs_.resize(interp.order() + 1);

    const auto dyads(spatial_dyads(dr));
    for(int i = 0; i <= interp.order(); ++i) {
      this->coefs_[i] = Eigen::Matrix3cd::Zero();

	    // if ( dr.norm() > dist0_ ) {
        this->coefs_[i] = 
            -k2_ * std::exp(-iu * omega_ * dr.norm() / c_) *
              (dyads[0].cast<cmplx>() * interp.evaluations[0][i] +
               dyads[1].cast<cmplx>() * (interp.evaluations[1][i] +
                                       iu * omega_ * interp.evaluations[0][i]) +
               dyads[2].cast<cmplx>() *
                 (interp.evaluations[2][i] +
                  2.0 * iu * omega_ * interp.evaluations[1][i] -
                  std::pow(omega_, 2) * interp.evaluations[0][i]));

      /*} else if ( dr.norm() == 0.0 ){

         this->coefs_[i] = -k2_ * 
          (dyads[0].cast<cmplx>() * interp.evaluations[0][i] +
           dyads[3].cast<cmplx>() * 
              (interp.evaluations[2][i] +
               2.0 * iu * omega_ * interp.evaluations[1][i] - 
               std::pow(omega_, 2) * interp.evaluations[0][i] ) );
   
          this->coefs_[i] = 
            beta_ / pow( 5.2917721e-4, 2 ) * Eigen::Matrix3d::Identity() * 
            ( 1.0 * iu * pow(omega_,3) * interp.evaluations[0][i] +
              3.0 * pow(omega_,2) * interp.evaluations[1][i] -
              3.0 * iu * omega_ * interp.evaluations[2][i] -
              interp.evaluations[3][i] );
      } */   
    }

    return this->coefs_;
  }

 private:
  double omega_;

};

template <class T>
class Propagation::MFIE : public Propagation::Kernel<T> {
 public:
  MFIE(const double c, const double k2) : c_(c), k2_(k2){};
  const std::vector<Mat3D<T>> &coefficients(
      const Eigen::Vector3d &dr,
      const Interpolation::UniformLagrangeSet &interp)
  {
    this->coefs_.resize(interp.order() + 1);

    for(int i = 0; i <= interp.order(); ++i) {
      this->coefs_[i] = Eigen::Matrix3d::Zero();

      if ( dr.norm() > 0.0 ) {
        this->coefs_[i] = -k2_ * Eigen::Matrix3d::Identity() * 
                            ( interp.evaluations[1][i] / dr.squaredNorm() +
                              interp.evaluations[2][i] / ( c_ * dr.norm() ) );
      }
    }
    return this->coefs_;
  }

 protected:
  double c_, k2_;

};

class Propagation::RotatingMFIE : public Propagation::MFIE<cmplx> {
 public:
  RotatingMFIE(const double c, const double k2, const double omega)
       : MFIE<cmplx>(c, k2), omega_(omega){};

 const std::vector<Eigen::Matrix3cd> &coefficients(
      const Eigen::Vector3d &dr,
      const Interpolation::UniformLagrangeSet &interp)
  {
    this->coefs_.resize(interp.order() + 1);

    for(int i = 0; i <= interp.order(); ++i) {
      this->coefs_[i] = Eigen::Matrix3d::Zero();

      if ( dr.norm() > 0.0 ) {
        this->coefs_[i] = 
          -k2_ * std::exp(-iu * omega_ * dr.norm() / c_) * Eigen::Matrix3cd::Identity() *
              ( (interp.evaluations[1][i] + 
                  iu * omega_ * interp.evaluations[0][i]) / dr.squaredNorm() +
                (interp.evaluations[2][i] +
                  2.0 * iu * omega_ * interp.evaluations[1][i] -
                  std::pow(omega_, 2) * interp.evaluations[0][i]) / ( c_ * dr.norm() ) );
        // std::cout << i << " " << this->coefs_[i] << std::endl;
      }
    }
    return this->coefs_;
  }

 protected:
  double omega_;

};

class Propagation::SelfEFIE : public Propagation::Kernel<cmplx> {
 public:
  SelfEFIE(const double c, const double k2, const double beta)
      : c_(c), k2_(k2), beta_(beta){};

 const std::vector<Eigen::Matrix3cd> &coefficients(
      const Eigen::Vector3d &dr,
      const Interpolation::UniformLagrangeSet &interp)
  {
    //assert( dr.norm() == 0.0 );
    this->coefs_.resize(interp.order() + 1);

    if ( dr.norm() == 0.0){
      for(int i = 0; i <= interp.order(); ++i) {
        this->coefs_[i] = 
          -beta_ / pow( 5.2917721e-4, 2 ) * Eigen::Matrix3d::Identity() * 
          // -k2_ * 2.0 * Eigen::Matrix3d::Identity() / ( 3.0 * c_ ) * 
                interp.evaluations[3][i] ;

      }
    }
    return this->coefs_;
  }
 protected:
  double c_, k2_, beta_;
};


class Propagation::SelfRotatingEFIE : public Propagation::SelfEFIE {
 public:
  SelfRotatingEFIE(const double c, const double k2, const double omega, const double beta)
      : SelfEFIE(c, k2, beta), omega_(omega){};

  const std::vector<Eigen::Matrix3cd> &coefficients(
      const Eigen::Vector3d &dr,
      const Interpolation::UniformLagrangeSet &interp)
  {
    // assert( dr.norm() == 0.0 );
    this->coefs_.resize(interp.order() + 1);

    if ( dr.norm() == 0.0){
      for(int i = 0; i <= interp.order(); ++i) {
        this->coefs_[i] = 
          beta_ / pow( 5.2917721e-4, 2 ) * Eigen::Matrix3d::Identity() * 
          // k2_ * 2.0 * Eigen::Matrix3d::Identity() / ( 3.0 * c_ ) * 
            ( 1.0 * iu * pow(omega_,3) * interp.evaluations[0][i] +
              3.0 * pow(omega_,2) * interp.evaluations[1][i] -
              3.0 * iu * omega_ * interp.evaluations[2][i] -
              interp.evaluations[3][i] ) ;
      }
    }
  return this->coefs_;
  }

 private:
  double omega_;

};
#endif
