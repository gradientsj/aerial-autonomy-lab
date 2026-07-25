// Core value types and a portable RNG for the aerial-autonomy-lab planning core.
//
// Portability note that drives the RNG design: libstdc++, libc++ and MSVC's STL
// all implement std::uniform_real_distribution and std::normal_distribution
// differently. The standard specifies the *distribution*, not the bit-exact
// sequence, so identical seeds produce different draws on different platforms.
// Golden-trajectory tests would therefore pass on Linux and fail on Windows for
// no real reason. We generate every random number from raw mt19937_64 bits with
// our own arithmetic, which the standard *does* pin down exactly.
#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <random>

namespace aal {

using Scalar = double;

struct Vec3 {
  Scalar x{0}, y{0}, z{0};

  constexpr Vec3() = default;
  constexpr Vec3(Scalar x_, Scalar y_, Scalar z_) : x(x_), y(y_), z(z_) {}

  constexpr Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
  constexpr Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
  constexpr Vec3 operator*(Scalar s) const { return {x * s, y * s, z * s}; }
  constexpr Vec3 operator/(Scalar s) const { return {x / s, y / s, z / s}; }
  Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }

  constexpr Scalar dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
  constexpr Scalar squaredNorm() const { return x * x + y * y + z * z; }
  Scalar norm() const { return std::sqrt(squaredNorm()); }

  Vec3 normalized() const {
    const Scalar n = norm();
    return n > 0 ? *this / n : Vec3{0, 0, 0};
  }
  Scalar operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }
  Scalar& operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }
};

inline Scalar distance(const Vec3& a, const Vec3& b) { return (a - b).norm(); }

// Axis-aligned bounding box, used for the world bounds and for map extents.
struct Aabb {
  Vec3 lo, hi;

  constexpr Aabb() = default;
  constexpr Aabb(const Vec3& lo_, const Vec3& hi_) : lo(lo_), hi(hi_) {}

  Vec3 extent() const { return hi - lo; }
  Scalar volume() const {
    const Vec3 e = extent();
    return e.x * e.y * e.z;
  }
  bool contains(const Vec3& p) const {
    return p.x >= lo.x && p.x <= hi.x && p.y >= lo.y && p.y <= hi.y && p.z >= lo.z && p.z <= hi.z;
  }
};

// Deterministic, cross-platform RNG.
//
// mt19937_64 itself is bit-exact across every conforming standard library (the
// standard fixes its state transition and even mandates the 10000th output of
// the default-seeded engine), so building our own distributions on top of it
// gives reproducibility we can actually assert in tests.
class Rng {
 public:
  explicit Rng(std::uint64_t seed = 0x9E3779B97F4A7C15ULL) : eng_(seed), state_seed_(seed) {}

  std::uint64_t nextBits() { return eng_(); }

  // Uniform in [0, 1). Uses the top 53 bits, which is exactly the number of
  // bits in a double's mantissa, so every representable value in [0,1) with
  // 53-bit spacing is reachable and none is favoured by rounding.
  Scalar uniform() {
    return static_cast<Scalar>(eng_() >> 11) * (1.0 / 9007199254740992.0);  // 2^-53
  }

  Scalar uniform(Scalar lo, Scalar hi) { return lo + (hi - lo) * uniform(); }

  // Marsaglia polar method. Chosen over Box-Muller because it avoids the
  // library sin/cos, whose last-bit results are not guaranteed identical
  // across platforms. Only sqrt and log are used, and both are correctly
  // rounded under IEEE 754 in practice on the platforms we target.
  Scalar normal() {
    if (hasSpare_) {
      hasSpare_ = false;
      return spare_;
    }
    Scalar u, v, s;
    do {
      u = uniform(-1.0, 1.0);
      v = uniform(-1.0, 1.0);
      s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);
    const Scalar f = std::sqrt(-2.0 * std::log(s) / s);
    spare_ = v * f;
    hasSpare_ = true;
    return u * f;
  }

  Vec3 uniformInAabb(const Aabb& b) {
    return {uniform(b.lo.x, b.hi.x), uniform(b.lo.y, b.hi.y), uniform(b.lo.z, b.hi.z)};
  }

  // Uniform inside the unit 3-ball. Sampling a direction from three normals and
  // scaling the radius by u^(1/3) gives a uniform density; rejection sampling
  // from the cube would also work but wastes 1 - pi/6 = 48% of draws.
  Vec3 uniformInUnitBall() {
    const Vec3 dir = Vec3{normal(), normal(), normal()}.normalized();
    const Scalar r = std::cbrt(uniform());
    return dir * r;
  }

  std::uint64_t seedState() const { return state_seed_; }

 private:
  std::mt19937_64 eng_;
  Scalar spare_{0};
  bool hasSpare_{false};
  std::uint64_t state_seed_{0};
};

}  // namespace aal
