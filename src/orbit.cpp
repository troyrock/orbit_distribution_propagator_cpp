#include <algorithm>
#include <cmath>
#include <distribution/orbit.hpp>
#include <stdexcept>

namespace distribution {
double wrap_angle(double x) {
    double y = std::fmod(x, tau);
    return y < 0.0 ? y + tau : y;
}
double dot(const Vector3 &a, const Vector3 &b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector3 cross(const Vector3 &a, const Vector3 &b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double norm(const Vector3 &a) {
    return std::sqrt(dot(a, a));
}
std::array<Vector3, 3> rtn_basis(const Cartesian &s) {
    Vector3 r{s[0], s[1], s[2]}, v{s[3], s[4], s[5]};
    auto h = cross(r, v);
    const double radius = norm(r), hn = norm(h);
    if (!(radius > 0.0 && hn > 0.0))
        throw std::invalid_argument("Degenerate Cartesian orbit");
    for (int k = 0; k < 3; ++k) {
        r[k] /= radius;
        h[k] /= hn;
    }
    return {r, cross(h, r), h};
}
Elements from_keplerian(double a, double e, double i, double node, double arg, double m) {
    if (!(i >= 0.0 && i < pi - 1e-8 && e >= 0.0 && e < 1.0))
        throw std::invalid_argument("Require 0 <= eccentricity < 1 and inclination < 180 degrees");
    const double p = node + arg, h = std::tan(i / 2.0);
    return {a, e * std::cos(p), e * std::sin(p), h * std::cos(node), h * std::sin(node), p + m};
}
void validate_elements(const Elements &e, double mu, double min_perigee) {
    for (double v : e)
        if (!std::isfinite(v))
            throw std::invalid_argument("Nonfinite orbital element");
    const double ecc = std::hypot(e[1], e[2]);
    if (!(mu > 0.0 && std::isfinite(mu) && e[0] > 0.0 && ecc < 1.0))
        throw std::invalid_argument("Require finite mu > 0 and bound elliptic orbit");
    if (std::hypot(e[3], e[4]) > 1e6)
        throw std::invalid_argument("Near-retrograde equinoctial singularity is unsupported");
    if (e[0] * (1.0 - ecc) < min_perigee)
        throw std::invalid_argument("Orbit perigee below configured no-drag validity floor");
}
Cartesian to_cartesian(const Elements &e, double mu) {
    validate_elements(e, mu);
    const double a = e[0], ex = e[1], ey = e[2], q = e[3], p = e[4];
    const double lm = std::remainder(e[5], tau);
    // Solve eccentric longitude F-ex sin(F)+ey cos(F)=lambda with a
    // monotonic bracket. This also works for e=0 without choosing perigee.
    const double ecc = std::hypot(ex, ey);
    double lo = lm - ecc, hi = lm + ecc, f = lm;
    for (int it = 0; it < 80; ++it) {
        const double value = f - ex * std::sin(f) + ey * std::cos(f) - lm;
        if (std::abs(value) < 2e-15)
            break;
        if (value > 0)
            hi = f;
        else
            lo = f;
        const double next = f - value / (1.0 - ex * std::cos(f) - ey * std::sin(f));
        f = (next > lo && next < hi) ? next : 0.5 * (lo + hi);
        if (it == 79)
            throw std::runtime_error("Eccentric longitude solve did not converge");
    }
    const double b = 1.0 / (1.0 + std::sqrt(1.0 - ecc * ecc)), c = std::cos(f), s = std::sin(f);
    const double x = a * ((1 - b * ey * ey) * c + b * ex * ey * s - ex);
    const double y = a * ((1 - b * ex * ex) * s + b * ex * ey * c - ey);
    const double fd = std::sqrt(mu / a) / a / (1.0 - ex * c - ey * s);
    const double xd = a * (-(1 - b * ey * ey) * s + b * ex * ey * c) * fd;
    const double yd = a * ((1 - b * ex * ex) * c - b * ex * ey * s) * fd;
    const double den = 1 + p * p + q * q;
    const Vector3 u{(1 - p * p + q * q) / den, 2 * p * q / den, -2 * p / den};
    const Vector3 v{2 * p * q / den, (1 + p * p - q * q) / den, 2 * q / den};
    Cartesian out{};
    for (int k = 0; k < 3; ++k) {
        out[k] = x * u[k] + y * v[k];
        out[k + 3] = xd * u[k] + yd * v[k];
    }
    return out;
}
Elements from_cartesian(const Cartesian &s, double mu) {
    for (double x : s)
        if (!std::isfinite(x))
            throw std::invalid_argument("Nonfinite Cartesian state");
    if (!(mu > 0 && std::isfinite(mu)))
        throw std::invalid_argument("Invalid gravitational parameter");
    Vector3 r{s[0], s[1], s[2]}, v{s[3], s[4], s[5]}, h = cross(r, v);
    const double rn = norm(r), hn = norm(h), v2 = dot(v, v);
    if (!(rn > 0 && hn > 0))
        throw std::invalid_argument("Degenerate Cartesian orbit");
    const double a = 1.0 / (2.0 / rn - v2 / mu);
    Vector3 w = h;
    for (auto &x : w)
        x /= hn;
    if (1 + w[2] < 1e-12)
        throw std::invalid_argument("Retrograde equinoctial singularity");
    const double q = -w[1] / (1 + w[2]), p = w[0] / (1 + w[2]), den = 1 + p * p + q * q;
    const Vector3 u{(1 - p * p + q * q) / den, 2 * p * q / den, -2 * p / den};
    const Vector3 vv{2 * p * q / den, (1 + p * p - q * q) / den, 2 * q / den};
    auto ev = cross(v, h);
    for (int k = 0; k < 3; ++k)
        ev[k] = ev[k] / mu - r[k] / rn;
    const double ex = dot(ev, u), ey = dot(ev, vv), ecc = std::hypot(ex, ey);
    Elements out{a, ex, ey, q, p, 0};
    validate_elements(out, mu);
    const double beta = 1.0 / (1 + std::sqrt(1 - ecc * ecc));
    const double xx = dot(r, u) / a + ex, yy = dot(r, vv) / a + ey;
    const double m11 = 1 - beta * ey * ey, m22 = 1 - beta * ex * ex, m12 = beta * ex * ey;
    const double det = m11 * m22 - m12 * m12;
    const double cf = (m22 * xx - m12 * yy) / det, sf = (m11 * yy - m12 * xx) / det;
    const double f = std::atan2(sf, cf);
    out[5] = wrap_angle(f - ex * sf + ey * cf);
    return out;
}
Matrix6 covariance_factor(const Matrix6 &c) {
    Matrix6 correlation{}, lower{}, out{};
    std::array<double, 6> scale{};
    for (int i = 0; i < 6; ++i) {
        if (!std::isfinite(c[i][i]) || c[i][i] < 0)
            throw std::invalid_argument("Negative/nonfinite covariance diagonal");
        scale[i] = std::sqrt(c[i][i]);
    }
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) {
            if (!std::isfinite(c[i][j]))
                throw std::invalid_argument("Nonfinite covariance");
            const double tolerance = 1e-12 * std::max({std::abs(c[i][j]), std::abs(c[j][i]),
                                                       scale[i] * scale[j], 1e-300});
            if (std::abs(c[i][j] - c[j][i]) > tolerance)
                throw std::invalid_argument("Covariance must be symmetric");
            if (scale[i] == 0 || scale[j] == 0) {
                if (c[i][j] != 0)
                    throw std::invalid_argument("Zero-variance covariance row must be zero");
                correlation[i][j] = 0;
            } else
                correlation[i][j] = c[i][j] / scale[i] / scale[j];
        }
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j <= i; ++j) {
            double x = correlation[i][j];
            for (int k = 0; k < j; ++k)
                x -= lower[i][k] * lower[j][k];
            if (i == j) {
                if (x < -1e-12)
                    throw std::invalid_argument("Covariance is not positive semidefinite");
                lower[i][j] = std::sqrt(std::max(0.0, x));
            } else if (lower[j][j] > 1e-14)
                lower[i][j] = x / lower[j][j];
            else if (std::abs(x) > 1e-12)
                throw std::invalid_argument("Covariance is not positive semidefinite");
            out[i][j] = scale[i] * lower[i][j];
        }
    return out;
}
} // namespace distribution
