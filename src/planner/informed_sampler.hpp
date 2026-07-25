// Direct sampling of the informed set for Informed RRT* (Gammell et al., 2014).
//
// Once a solution of cost c_best exists, the only states that can possibly
// improve it are those satisfying
//
//     ||x - x_start|| + ||x - x_goal||  <=  c_best
//
// which is a prolate hyperspheroid with foci at start and goal. Sampling it by
// rejection from the bounding box is correct but degenerates badly: as c_best
// approaches c_min the acceptance rate goes to zero. This samples it directly
// by transforming a uniform unit-ball draw, so the cost per sample is constant
// no matter how tight the ellipsoid gets.
#pragma once

#include <array>
#include <cmath>

#include "core/types.hpp"

namespace aal {

class InformedSampler {
 public:
  InformedSampler(const Vec3& start, const Vec3& goal)
      : centre_((start + goal) * 0.5), c_min_((goal - start).norm()) {
    buildRotation(goal - start);
  }

  Scalar cMin() const { return c_min_; }

  // Volume of the informed set at cost c_best. In 3D,
  //   V = zeta_3 * (c_best / 2) * ( sqrt(c_best^2 - c_min^2) / 2 )^2
  // with zeta_3 = 4 pi / 3. Used to shrink the RRT* connection radius once
  // sampling is restricted, since leaving mu(X_free) in the formula would keep
  // the radius larger than it needs to be.
  Scalar volume(Scalar c_best) const {
    if (!(c_best > c_min_)) return 0.0;
    const Scalar a = c_best * 0.5;
    const Scalar b = std::sqrt(c_best * c_best - c_min_ * c_min_) * 0.5;
    return (4.0 / 3.0) * M_PI * a * b * b;
  }

  // Draws a point uniformly from the prolate hyperspheroid of cost c_best.
  // Requires c_best >= c_min; the caller checks that before switching over.
  Vec3 sample(Rng& rng, Scalar c_best) const {
    const Scalar a = c_best * 0.5;
    const Scalar transverse = std::sqrt(std::max<Scalar>(0, c_best * c_best - c_min_ * c_min_)) * 0.5;
    const Vec3 ball = rng.uniformInUnitBall();
    // Scale the unit ball into the ellipsoid, then rotate it onto the
    // start-to-goal axis and translate to the midpoint.
    const Vec3 scaled{ball.x * a, ball.y * transverse, ball.z * transverse};
    return rotate(scaled) + centre_;
  }

 private:
  // Builds a rotation whose first column is the unit start-to-goal direction.
  // The literature derives this from the SVD of a1 * e1^T, but in 3D an
  // explicit orthonormal completion is exact, cheaper and easier to verify.
  // The determinant is forced to +1 so the map is a rotation, not a reflection.
  void buildRotation(const Vec3& delta) {
    const Vec3 a1 = c_min_ > 0 ? delta / c_min_ : Vec3{1, 0, 0};
    // Pick the world axis least aligned with a1 to avoid a degenerate cross
    // product when a1 happens to lie along an axis.
    const Scalar ax = std::abs(a1.x), ay = std::abs(a1.y), az = std::abs(a1.z);
    Vec3 helper{1, 0, 0};
    if (ay <= ax && ay <= az) helper = Vec3{0, 1, 0};
    else if (az <= ax && az <= ay) helper = Vec3{0, 0, 1};

    Vec3 a2 = cross(a1, helper).normalized();
    Vec3 a3 = cross(a1, a2).normalized();
    // cross(a1, a2) already yields a right-handed triple, so det = +1.
    col_[0] = a1;
    col_[1] = a2;
    col_[2] = a3;
  }

  static Vec3 cross(const Vec3& u, const Vec3& v) {
    return {u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
  }

  Vec3 rotate(const Vec3& v) const { return col_[0] * v.x + col_[1] * v.y + col_[2] * v.z; }

  Vec3 centre_;
  Scalar c_min_;
  std::array<Vec3, 3> col_;
};

}  // namespace aal
