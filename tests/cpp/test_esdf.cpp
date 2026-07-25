// Tests for the Euclidean signed distance field and the collision checker.
//
// The central test is not "does it look right" but "does the O(n) separable
// transform agree exactly with brute force". Felzenszwalb's algorithm is exact,
// so the tolerance is floating-point noise, not a fudge factor. Chamfer or
// two-pass mask transforms would fail this test, which is the point of writing
// it: it pins down that we have the exact transform and not an approximation.

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "collision/esdf.hpp"

using namespace aal;

namespace {

Esdf makeSmallMap(int n = 12, Scalar res = 0.5) {
  Esdf f(n, n, n, res, Vec3{0, 0, 0});
  return f;
}

}  // namespace

TEST(Esdf, MatchesBruteForceExactly) {
  // Several occupancy patterns, including ones that stress the lower envelope:
  // a single voxel, a plane, and a scattered set.
  for (int pattern = 0; pattern < 3; ++pattern) {
    Esdf f = makeSmallMap();
    if (pattern == 0) {
      f.setOccupied(6, 6, 6, true);
    } else if (pattern == 1) {
      for (int j = 0; j < 12; ++j)
        for (int i = 0; i < 12; ++i) f.setOccupied(i, j, 5, true);
    } else {
      for (int k = 0; k < 12; k += 5)
        for (int j = 1; j < 12; j += 4)
          for (int i = 2; i < 12; i += 3) f.setOccupied(i, j, k, true);
    }
    f.build();
    for (int k = 0; k < 12; ++k)
      for (int j = 0; j < 12; ++j)
        for (int i = 0; i < 12; ++i) {
          EXPECT_NEAR(f.distanceAt(i, j, k), f.bruteForceDistanceAt(i, j, k), 1e-9)
              << "pattern " << pattern << " at " << i << "," << j << "," << k;
        }
  }
}

TEST(Esdf, EmptyGridDoesNotProduceNaN) {
  // Regression test. Seeding the distance transform with infinity makes the
  // parabola intersection compute INF - INF = NaN on any scanline with no
  // obstacle, which silently corrupts the whole field. A grid with no obstacles
  // at all is the extreme case.
  Esdf f = makeSmallMap(8);
  f.build();
  for (int k = 0; k < 8; ++k)
    for (int j = 0; j < 8; ++j)
      for (int i = 0; i < 8; ++i) {
        EXPECT_FALSE(std::isnan(f.distanceAt(i, j, k)));
        EXPECT_GT(f.distanceAt(i, j, k), 0.0);
      }
}

TEST(Esdf, MatchesTheAnalyticDistanceToAHalfSpace) {
  // Regression test for a half-voxel bias.
  //
  // A distance transform measures centre-to-centre distance to the nearest
  // voxel of the opposite class, but the planner needs distance to the obstacle
  // SURFACE. An occupied voxel is a cube of side res centred on its sample
  // point, so the free voxel beside it sits half a voxel from the surface and a
  // whole voxel from the centre. Reporting the latter inflates every clearance
  // by half a voxel, which at res = 0.5 m is 0.25 m the vehicle does not have.
  //
  // A half space is the one geometry where the correct answer is available in
  // closed form at every voxel, so the agreement here is exact rather than
  // approximate.
  const Scalar res = 0.5;
  const int wall_max_i = 5;  // voxels 0..5 occupied, surface at i = 5.5
  Esdf f(16, 6, 6, res, Vec3{0, 0, 0});
  for (int k = 0; k < 6; ++k)
    for (int j = 0; j < 6; ++j)
      for (int i = 0; i <= wall_max_i; ++i) f.setOccupied(i, j, k, true);
  f.build();

  for (int i = 0; i < 16; ++i) {
    const Scalar truth = (Scalar(i) - (wall_max_i + 0.5)) * res;
    EXPECT_NEAR(f.distanceAt(i, 3, 3), truth, 1e-12) << "at voxel i = " << i;
  }
}

TEST(Esdf, SignedInsideObstacles) {
  Esdf f(20, 20, 20, 0.5, Vec3{0, 0, 0});
  f.addSphere(Vec3{5, 5, 5}, 2.0);
  f.build();
  // Centre of the sphere is well inside, so the signed distance is negative and
  // roughly the radius.
  const Scalar d_centre = f.distance(Vec3{5, 5, 5});
  EXPECT_LT(d_centre, 0.0);
  EXPECT_NEAR(std::abs(d_centre), 2.0, 0.6);
  // A point far outside is positive.
  EXPECT_GT(f.distance(Vec3{0.5, 0.5, 0.5}), 0.0);
}

TEST(Esdf, IsOneLipschitz) {
  // A Euclidean distance field satisfies |d(a) - d(b)| <= |a - b|. The sphere
  // traced collision checker's soundness depends entirely on this property, so
  // it is worth asserting directly rather than assuming it.
  Esdf f(24, 24, 24, 0.4, Vec3{0, 0, 0});
  f.addSphere(Vec3{4, 4, 4}, 1.5);
  f.addBox(Aabb{Vec3{6, 1, 1}, Vec3{7, 8, 8}});
  f.build();

  Rng rng(7);
  const Aabb b = f.bounds();
  for (int t = 0; t < 20000; ++t) {
    const Vec3 a = rng.uniformInAabb(b);
    const Vec3 c = rng.uniformInAabb(b);
    const Scalar lhs = std::abs(f.distance(a) - f.distance(c));
    const Scalar rhs = (a - c).norm();
    // Trilinear interpolation of a 1-Lipschitz sampled field can exceed the
    // bound by at most a discretisation term of order the voxel diagonal.
    EXPECT_LE(lhs, rhs + 1e-9 + f.resolution() * std::sqrt(3.0));
  }
}

TEST(Esdf, GradientMatchesFiniteDifference) {
  Esdf f(24, 24, 24, 0.4, Vec3{0, 0, 0});
  f.addSphere(Vec3{4, 4, 4}, 1.5);
  f.build();

  Rng rng(11);
  const Scalar h = 1e-5;
  int checked = 0;
  for (int t = 0; t < 3000 && checked < 500; ++t) {
    const Vec3 p = rng.uniformInAabb(Aabb{Vec3{1, 1, 1}, Vec3{8, 8, 8}});
    // Skip points within h of a cell boundary: the trilinear interpolant is
    // only piecewise smooth, so its gradient is genuinely discontinuous there
    // and a central difference straddling the seam is meaningless.
    const Vec3 g = (p - f.origin()) / f.resolution();
    auto near_seam = [&](Scalar v) {
      const Scalar frac = v - std::floor(v);
      return frac < 1e-3 || frac > 1 - 1e-3;
    };
    if (near_seam(g.x) || near_seam(g.y) || near_seam(g.z)) continue;
    ++checked;

    const Vec3 grad = f.gradient(p);
    const Scalar dx = (f.distance(p + Vec3{h, 0, 0}) - f.distance(p - Vec3{h, 0, 0})) / (2 * h);
    const Scalar dy = (f.distance(p + Vec3{0, h, 0}) - f.distance(p - Vec3{0, h, 0})) / (2 * h);
    const Scalar dz = (f.distance(p + Vec3{0, 0, h}) - f.distance(p - Vec3{0, 0, h})) / (2 * h);
    EXPECT_NEAR(grad.x, dx, 1e-4);
    EXPECT_NEAR(grad.y, dy, 1e-4);
    EXPECT_NEAR(grad.z, dz, 1e-4);
  }
  EXPECT_GT(checked, 100);
}

TEST(CollisionChecker, SoundnessAgainstDenseSampling) {
  // The sphere trace takes long strides. This test verifies it never certifies
  // an edge that a brute-force dense sampling finds to be in collision.
  Esdf f(30, 30, 30, 0.4, Vec3{0, 0, 0});
  f.addSphere(Vec3{6, 6, 6}, 2.0);
  f.addCylinderZ(3.0, 8.0, 1.0, 0.0, 12.0);
  f.build();

  const Scalar radius = 0.35;
  CollisionChecker cc(f, radius);
  Rng rng(3);
  const Aabb b = f.bounds();
  int free_edges = 0, blocked_edges = 0;

  for (int t = 0; t < 4000; ++t) {
    const Vec3 a = rng.uniformInAabb(b);
    const Vec3 c = rng.uniformInAabb(b);
    const bool traced = cc.edgeFree(a, c);

    // Dense reference at a step far below the voxel size.
    const Scalar L = (c - a).norm();
    const int steps = std::max(2, static_cast<int>(L / 0.01));
    bool dense_free = true;
    for (int s = 0; s <= steps; ++s) {
      const Vec3 p = a + (c - a) * (Scalar(s) / steps);
      if (f.distance(p) < radius) {
        dense_free = false;
        break;
      }
    }
    if (traced) {
      ++free_edges;
      EXPECT_TRUE(dense_free) << "sphere trace certified an edge that is in collision";
    } else {
      ++blocked_edges;
    }
  }
  // Guard against a vacuous pass where every edge happened to be blocked.
  EXPECT_GT(free_edges, 100);
  EXPECT_GT(blocked_edges, 100);
}

TEST(CollisionChecker, StrideDoesNotDependOnEdgeDirection) {
  // edgeFree(a,b) and edgeFree(b,a) traverse the same segment with different
  // stride sequences. They must still agree, or the planner's behaviour would
  // depend on the arbitrary orientation of an edge.
  Esdf f(30, 30, 30, 0.4, Vec3{0, 0, 0});
  f.addSphere(Vec3{6, 6, 6}, 2.0);
  f.build();
  CollisionChecker cc(f, 0.3);
  Rng rng(19);
  const Aabb b = f.bounds();
  for (int t = 0; t < 3000; ++t) {
    const Vec3 a = rng.uniformInAabb(b);
    const Vec3 c = rng.uniformInAabb(b);
    EXPECT_EQ(cc.edgeFree(a, c), cc.edgeFree(c, a));
  }
}
