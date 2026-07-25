// Euclidean signed distance field over a voxel grid, and the collision checker
// built on it.
//
// Why an ESDF rather than raw occupancy: a planner that only knows "occupied /
// free" must sample an edge at a fixed resolution fine enough for the thinnest
// obstacle, which is wasteful in open space and still unsound near thin walls.
// A distance field gives, at every point, the radius of a certifiably free ball.
// That turns edge checking into sphere tracing, which takes long strides through
// open space and short ones only where it matters, and it gives the trajectory
// optimiser a smooth gradient to push against.
#pragma once

#include <cstdint>
#include <vector>

#include "core/types.hpp"

namespace aal {

class Esdf {
 public:
  // origin is the world coordinate of the *centre* of voxel (0,0,0).
  Esdf(int nx, int ny, int nz, Scalar resolution, const Vec3& origin);

  int nx() const { return nx_; }
  int ny() const { return ny_; }
  int nz() const { return nz_; }
  Scalar resolution() const { return res_; }
  const Vec3& origin() const { return origin_; }
  Aabb bounds() const;

  // Occupancy accessors, used to author a map before building the field.
  void setOccupied(int i, int j, int k, bool occupied);
  bool occupied(int i, int j, int k) const;
  void clearOccupancy();

  // Convenience authoring primitives for tests and procedural maps.
  void addSphere(const Vec3& centre, Scalar radius);
  void addBox(const Aabb& box);
  void addCylinderZ(Scalar cx, Scalar cy, Scalar radius, Scalar z_lo, Scalar z_hi);

  // Computes the signed distance field from the current occupancy.
  // Positive outside obstacles, negative inside. Exact to within the voxel
  // discretisation: the transform itself introduces no additional error beyond
  // representing obstacles as voxel centres.
  void build();

  // Distance at a world point, trilinearly interpolated. Points outside the
  // grid return the distance to the grid boundary clamped sample, which makes
  // the world bounds behave like an obstacle.
  Scalar distance(const Vec3& p) const;

  // Analytic gradient of the trilinear interpolant. Not merely a finite
  // difference of `distance`: the trilinear form is piecewise multilinear, so
  // its gradient is available in closed form and is exact within each cell.
  Vec3 gradient(const Vec3& p) const;

  // Raw voxel access to the built field, in metres.
  Scalar distanceAt(int i, int j, int k) const;

  // Fraction of the grid that is free. Used to estimate mu(X_free) for the
  // RRT* radius constant, which is otherwise unknown in an unmapped world.
  Scalar freeFraction() const;

  bool isBuilt() const { return built_; }

  // Reference implementation used only by tests: O(n_free * n_occ) brute force.
  // Slow by construction, correct by construction.
  Scalar bruteForceDistanceAt(int i, int j, int k) const;

 private:
  int idx(int i, int j, int k) const { return (k * ny_ + j) * nx_ + i; }
  bool inGrid(int i, int j, int k) const {
    return i >= 0 && i < nx_ && j >= 0 && j < ny_ && k >= 0 && k < nz_;
  }
  // Continuous grid coordinates of a world point.
  Vec3 toGrid(const Vec3& p) const;

  int nx_, ny_, nz_;
  Scalar res_;
  Vec3 origin_;
  std::vector<std::uint8_t> occ_;
  std::vector<Scalar> dist_;
  bool built_{false};
};

// Collision checking against an ESDF.
//
// `robot_radius` inflates the vehicle to a sphere. `margin` is an extra safety
// buffer. The checker is *sound* with respect to radius + margin: it never
// reports an edge free if any point on it comes within that clearance of an
// obstacle. It is conservative in the other direction only within one
// `min_step` of the boundary, which is documented and configurable.
class CollisionChecker {
 public:
  CollisionChecker(const Esdf& field, Scalar robot_radius, Scalar margin = 0.0)
      : field_(field), radius_(robot_radius), margin_(margin) {}

  Scalar clearance() const { return radius_ + margin_; }

  bool pointFree(const Vec3& p) const { return field_.distance(p) >= clearance(); }

  // Sphere-traced segment check. Returns true if every point of [a,b] has
  // clearance. Because the distance field is 1-Lipschitz, advancing by
  // d(p) - clearance can never step over a violation, which is what makes the
  // adaptive stride sound rather than merely fast.
  bool edgeFree(const Vec3& a, const Vec3& b) const;

  // Same as edgeFree but reports how far along the edge it got, in metres.
  // Used by the RRT* extend step to salvage a partial edge.
  Scalar freeArcLength(const Vec3& a, const Vec3& b) const;

  // Count of distance queries issued, for benchmarking. Mutable counter, reset
  // by the caller between runs.
  std::uint64_t queryCount() const { return queries_; }
  void resetQueryCount() const { queries_ = 0; }

 private:
  const Esdf& field_;
  Scalar radius_;
  Scalar margin_;
  mutable std::uint64_t queries_{0};
};

}  // namespace aal
