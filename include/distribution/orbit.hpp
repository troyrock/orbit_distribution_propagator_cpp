#pragma once
#include <array>

namespace distribution {
using Elements = std::array<double, 6>; // a, ex, ey, hx, hy, mean longitude
using Vector3 = std::array<double, 3>;
using Cartesian = std::array<double, 6>; // inertial r [m], v [m/s]
using Matrix6 = std::array<std::array<double, 6>, 6>;
inline constexpr double pi = 3.141592653589793238462643383279502884;
inline constexpr double tau = 2.0 * pi;
inline constexpr double day = 86400.0;

double wrap_angle(double angle);
double dot(const Vector3 &a, const Vector3 &b);
Vector3 cross(const Vector3 &a, const Vector3 &b);
double norm(const Vector3 &a);
std::array<Vector3, 3> rtn_basis(const Cartesian &state);
Elements from_keplerian(double a, double e, double inclination, double raan,
                        double argument_perigee, double mean_anomaly);
Cartesian to_cartesian(const Elements &elements, double mu);
Elements from_cartesian(const Cartesian &state, double mu);
void validate_elements(const Elements &elements, double mu, double minimum_perigee_m = 0.0);
// PSD factorization permits exactly deterministic axes; rejects asymmetric or
// indefinite covariances instead of silently changing the supplied distribution.
Matrix6 covariance_factor(const Matrix6 &covariance);
} // namespace distribution
