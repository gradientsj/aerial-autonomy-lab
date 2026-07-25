#include "collision/esdf.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

namespace aal {
namespace {

constexpr Scalar kInf = std::numeric_limits<Scalar>::infinity();

// Seed value for "no source here". Deliberately a large *finite* number rather
// than infinity: the parabola intersection below computes
// (f[q] + q^2) - (f[v] + v^2), and with two infinite seeds that is INF - INF,
// which is NaN. A NaN breakpoint silently corrupts the envelope and the bug
// only shows up on scanlines that contain no obstacle at all. With a finite
// sentinel the same expression collapses to the correct midpoint (q + v) / 2.
constexpr Scalar kBig = 1e18;

// Felzenszwalb and Huttenlocher, "Distance Transforms of Sampled Functions"
// (2012). Computes D(p) = min_q ( (p - q)^2 + f(q) ) for every p, in O(n).
//
// The idea: each sample q contributes a parabola rooted at (q, f(q)). The
// transform is the lower envelope of those parabolas. Scanning left to right,
// maintain the envelope as a stack of parabolas v[0..k] with breakpoints z[],
// popping any parabola the new one occludes. Each parabola is pushed and popped
// at most once, hence linear time. This is exact: no approximation, unlike the
// chamfer or two-pass masks often used for distance transforms.
void dt1d(const Scalar* f, Scalar* d, int n, int* v, Scalar* z) {
  int k = 0;
  v[0] = 0;
  z[0] = -kInf;
  z[1] = kInf;
  for (int q = 1; q < n; ++q) {
    // Intersection of the parabola from q with the one currently on top.
    Scalar s = ((f[q] + Scalar(q) * q) - (f[v[k]] + Scalar(v[k]) * v[k])) /
               (2 * Scalar(q) - 2 * Scalar(v[k]));
    // Termination is guaranteed without a k > 0 guard: z[0] is -inf and s is
    // finite, so the test fails at the bottom of the stack.
    while (s <= z[k]) {
      --k;
      s = ((f[q] + Scalar(q) * q) - (f[v[k]] + Scalar(v[k]) * v[k])) /
          (2 * Scalar(q) - 2 * Scalar(v[k]));
    }
    ++k;
    v[k] = q;
    z[k] = s;
    z[k + 1] = kInf;
  }
  k = 0;
  for (int q = 0; q < n; ++q) {
    while (z[k + 1] < Scalar(q)) ++k;
    const Scalar dq = Scalar(q) - Scalar(v[k]);
    d[q] = dq * dq + f[v[k]];
  }
}

// Separable 3D squared distance transform. Applying the 1D transform along each
// axis in turn is exact for the squared Euclidean metric because the squared
// distance separates as a sum over axes, which is precisely why the Euclidean
// case admits an O(n) algorithm while other metrics do not.
void edt3d(std::vector<Scalar>& f, int nx, int ny, int nz) {
  const int nmax = std::max({nx, ny, nz});
  std::vector<Scalar> buf_in(nmax), buf_out(nmax), z(nmax + 1);
  std::vector<int> v(nmax);
  auto at = [&](int i, int j, int k) -> Scalar& { return f[(k * ny + j) * nx + i]; };

  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) buf_in[i] = at(i, j, k);
      dt1d(buf_in.data(), buf_out.data(), nx, v.data(), z.data());
      for (int i = 0; i < nx; ++i) at(i, j, k) = buf_out[i];
    }
  }
  for (int k = 0; k < nz; ++k) {
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) buf_in[j] = at(i, j, k);
      dt1d(buf_in.data(), buf_out.data(), ny, v.data(), z.data());
      for (int j = 0; j < ny; ++j) at(i, j, k) = buf_out[j];
    }
  }
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      for (int k = 0; k < nz; ++k) buf_in[k] = at(i, j, k);
      dt1d(buf_in.data(), buf_out.data(), nz, v.data(), z.data());
      for (int k = 0; k < nz; ++k) at(i, j, k) = buf_out[k];
    }
  }
}

}  // namespace

Esdf::Esdf(int nx, int ny, int nz, Scalar resolution, const Vec3& origin)
    : nx_(nx), ny_(ny), nz_(nz), res_(resolution), origin_(origin) {
  occ_.assign(static_cast<std::size_t>(nx) * ny * nz, 0);
  dist_.assign(static_cast<std::size_t>(nx) * ny * nz, 0);
}

Aabb Esdf::bounds() const {
  const Vec3 half{res_ * 0.5, res_ * 0.5, res_ * 0.5};
  const Vec3 lo = origin_ - half;
  const Vec3 hi = origin_ + Vec3{res_ * (nx_ - 1), res_ * (ny_ - 1), res_ * (nz_ - 1)} + half;
  return {lo, hi};
}

void Esdf::setOccupied(int i, int j, int k, bool occupied) {
  if (!inGrid(i, j, k)) return;
  occ_[idx(i, j, k)] = occupied ? 1 : 0;
  built_ = false;
}

bool Esdf::occupied(int i, int j, int k) const {
  return inGrid(i, j, k) && occ_[idx(i, j, k)] != 0;
}

void Esdf::clearOccupancy() {
  std::fill(occ_.begin(), occ_.end(), 0);
  built_ = false;
}

void Esdf::addSphere(const Vec3& centre, Scalar radius) {
  for (int k = 0; k < nz_; ++k)
    for (int j = 0; j < ny_; ++j)
      for (int i = 0; i < nx_; ++i) {
        const Vec3 p = origin_ + Vec3{res_ * i, res_ * j, res_ * k};
        if ((p - centre).norm() <= radius) occ_[idx(i, j, k)] = 1;
      }
  built_ = false;
}

void Esdf::addBox(const Aabb& box) {
  for (int k = 0; k < nz_; ++k)
    for (int j = 0; j < ny_; ++j)
      for (int i = 0; i < nx_; ++i) {
        const Vec3 p = origin_ + Vec3{res_ * i, res_ * j, res_ * k};
        if (box.contains(p)) occ_[idx(i, j, k)] = 1;
      }
  built_ = false;
}

void Esdf::addCylinderZ(Scalar cx, Scalar cy, Scalar radius, Scalar z_lo, Scalar z_hi) {
  for (int k = 0; k < nz_; ++k)
    for (int j = 0; j < ny_; ++j)
      for (int i = 0; i < nx_; ++i) {
        const Vec3 p = origin_ + Vec3{res_ * i, res_ * j, res_ * k};
        const Scalar dx = p.x - cx, dy = p.y - cy;
        if (dx * dx + dy * dy <= radius * radius && p.z >= z_lo && p.z <= z_hi)
          occ_[idx(i, j, k)] = 1;
      }
  built_ = false;
}

void Esdf::build() {
  const std::size_t n = occ_.size();
  std::vector<Scalar> outside(n), inside(n);
  for (std::size_t i = 0; i < n; ++i) {
    // Seed the transform: zero cost at the source set, infinite elsewhere.
    outside[i] = occ_[i] ? 0.0 : kBig;  // distance from any voxel to an obstacle
    inside[i] = occ_[i] ? kBig : 0.0;   // distance from any voxel to free space
  }
  edt3d(outside, nx_, ny_, nz_);
  edt3d(inside, nx_, ny_, nz_);

  for (std::size_t i = 0; i < n; ++i) {
    // The transform measures centre-to-centre distance to the nearest voxel of
    // the opposite class. What the planner needs is distance to the obstacle
    // SURFACE. An obstacle voxel is a cube of side res_ centred on its sample
    // point, so the free voxel immediately beside it is 0.5 * res_ from the
    // surface but a full res_ from the centre.
    //
    // Without this correction the field is optimistic by half a voxel
    // everywhere, which at res_ = 0.5 m is 0.25 m of clearance the vehicle does
    // not actually have. Subtracting half a voxel is exact for a flat face and
    // conservative elsewhere, which is the correct direction to err.
    const Scalar half = 0.5;
    const Scalar d_out = (std::sqrt(outside[i]) - half) * res_;
    const Scalar d_in = (std::sqrt(inside[i]) - half) * res_;
    // Signed convention: positive in free space, negative inside obstacles.
    dist_[i] = occ_[i] ? -d_in : d_out;
  }
  built_ = true;
}

Scalar Esdf::distanceAt(int i, int j, int k) const {
  const int ci = std::clamp(i, 0, nx_ - 1);
  const int cj = std::clamp(j, 0, ny_ - 1);
  const int ck = std::clamp(k, 0, nz_ - 1);
  return dist_[idx(ci, cj, ck)];
}

Vec3 Esdf::toGrid(const Vec3& p) const { return (p - origin_) / res_; }

Scalar Esdf::distance(const Vec3& p) const {
  const Vec3 g = toGrid(p);
  const int i0 = static_cast<int>(std::floor(g.x));
  const int j0 = static_cast<int>(std::floor(g.y));
  const int k0 = static_cast<int>(std::floor(g.z));
  const Scalar tx = g.x - i0, ty = g.y - j0, tz = g.z - k0;

  const Scalar c000 = distanceAt(i0, j0, k0), c100 = distanceAt(i0 + 1, j0, k0);
  const Scalar c010 = distanceAt(i0, j0 + 1, k0), c110 = distanceAt(i0 + 1, j0 + 1, k0);
  const Scalar c001 = distanceAt(i0, j0, k0 + 1), c101 = distanceAt(i0 + 1, j0, k0 + 1);
  const Scalar c011 = distanceAt(i0, j0 + 1, k0 + 1), c111 = distanceAt(i0 + 1, j0 + 1, k0 + 1);

  const Scalar c00 = c000 * (1 - tx) + c100 * tx;
  const Scalar c10 = c010 * (1 - tx) + c110 * tx;
  const Scalar c01 = c001 * (1 - tx) + c101 * tx;
  const Scalar c11 = c011 * (1 - tx) + c111 * tx;
  const Scalar c0 = c00 * (1 - ty) + c10 * ty;
  const Scalar c1 = c01 * (1 - ty) + c11 * ty;
  return c0 * (1 - tz) + c1 * tz;
}

Vec3 Esdf::gradient(const Vec3& p) const {
  const Vec3 g = toGrid(p);
  const int i0 = static_cast<int>(std::floor(g.x));
  const int j0 = static_cast<int>(std::floor(g.y));
  const int k0 = static_cast<int>(std::floor(g.z));
  const Scalar tx = g.x - i0, ty = g.y - j0, tz = g.z - k0;

  const Scalar c000 = distanceAt(i0, j0, k0), c100 = distanceAt(i0 + 1, j0, k0);
  const Scalar c010 = distanceAt(i0, j0 + 1, k0), c110 = distanceAt(i0 + 1, j0 + 1, k0);
  const Scalar c001 = distanceAt(i0, j0, k0 + 1), c101 = distanceAt(i0 + 1, j0, k0 + 1);
  const Scalar c011 = distanceAt(i0, j0 + 1, k0 + 1), c111 = distanceAt(i0 + 1, j0 + 1, k0 + 1);

  // Partial derivatives of the trilinear form, in grid units, then converted to
  // metres by dividing by the resolution.
  const Scalar dx = ((c100 - c000) * (1 - ty) + (c110 - c010) * ty) * (1 - tz) +
                    ((c101 - c001) * (1 - ty) + (c111 - c011) * ty) * tz;
  const Scalar dy = ((c010 - c000) * (1 - tx) + (c110 - c100) * tx) * (1 - tz) +
                    ((c011 - c001) * (1 - tx) + (c111 - c101) * tx) * tz;
  const Scalar dz = ((c001 - c000) * (1 - tx) + (c101 - c100) * tx) * (1 - ty) +
                    ((c011 - c010) * (1 - tx) + (c111 - c110) * tx) * ty;
  return Vec3{dx, dy, dz} / res_;
}

Scalar Esdf::freeFraction() const {
  std::size_t free = 0;
  for (std::uint8_t o : occ_)
    if (!o) ++free;
  return occ_.empty() ? 0.0 : static_cast<Scalar>(free) / static_cast<Scalar>(occ_.size());
}

Scalar Esdf::bruteForceDistanceAt(int i, int j, int k) const {
  Scalar best = kInf;
  const bool self_occ = occupied(i, j, k);
  for (int c = 0; c < nz_; ++c)
    for (int b = 0; b < ny_; ++b)
      for (int a = 0; a < nx_; ++a) {
        const bool o = occ_[idx(a, b, c)] != 0;
        // Outside obstacles we seek the nearest occupied voxel; inside, the
        // nearest free voxel. Mirrors the two transforms in build().
        if (o == self_occ) continue;
        const Scalar dx = Scalar(a - i), dy = Scalar(b - j), dz = Scalar(c - k);
        best = std::min(best, dx * dx + dy * dy + dz * dz);
      }
  // Matches build()'s behaviour when a source set is empty, so the two agree
  // even in the degenerate all-free or all-occupied cases.
  if (!std::isfinite(best)) best = kBig;
  // Same half-voxel surface correction as build(), so the two agree.
  const Scalar d = (std::sqrt(best) - 0.5) * res_;
  return self_occ ? -d : d;
}

bool CollisionChecker::edgeFree(const Vec3& a, const Vec3& b) const {
  // Canonicalise the endpoint order. The sphere trace visits a different set of
  // sample points depending on which end it starts from, so without this the
  // answer could differ between edgeFree(a,b) and edgeFree(b,a) for an edge
  // that grazes the clearance boundary. The planner treats edges as undirected,
  // so the checker must be a function of the unordered pair.
  const bool swap = std::make_tuple(b.x, b.y, b.z) < std::make_tuple(a.x, a.y, a.z);
  const Vec3& p0 = swap ? b : a;
  const Vec3& p1 = swap ? a : b;
  const Scalar L = (p1 - p0).norm();
  return freeArcLength(p0, p1) >= L;
}

Scalar CollisionChecker::freeArcLength(const Vec3& a, const Vec3& b) const {
  const Scalar L = (b - a).norm();
  const Scalar clear = clearance();
  if (L == 0) {
    ++queries_;
    return field_.distance(a) >= clear ? 0.0 : -1.0;
  }
  const Vec3 dir = (b - a) / L;

  // The Lipschitz constant of the *interpolated* field, which is what we
  // actually query.
  //
  // The exact Euclidean distance function is 1-Lipschitz, and it is tempting to
  // stride by the full slack d - clearance on that basis. But we do not query
  // the exact function, we query its trilinear interpolant over the voxel grid.
  // Adjacent voxel samples differ by at most one voxel width, so each partial
  // derivative of the interpolant is bounded by 1 in magnitude, which bounds
  // the gradient norm by sqrt(3) rather than 1 along a diagonal direction.
  // Striding by the full slack therefore steps over violations, which is
  // exactly what CollisionChecker.SoundnessAgainstDenseSampling caught.
  const Scalar kLipschitz = 1.7320508075688772;  // sqrt(3)

  // Below this slack we are effectively on the clearance boundary. Reporting a
  // collision is the conservative direction, and it bounds the iteration count,
  // which a pure sphere trace does not: the stride shrinks geometrically as the
  // ray grazes a surface, so the trace can otherwise stall short of the end.
  const Scalar min_slack = 1e-4;
  constexpr int kMaxSteps = 1 << 16;

  Scalar t = 0;
  for (int step_count = 0; t <= L; ++step_count) {
    if (step_count >= kMaxSteps) return t;  // grazing, treat as blocked
    const Vec3 p = a + dir * t;
    ++queries_;
    const Scalar slack = field_.distance(p) - clear;
    if (slack < min_slack) return t;  // blocked here
    t += slack / kLipschitz;
  }
  // The last stride overshot b. Every stride was certified free over its whole
  // length by the Lipschitz bound, so b is free too.
  return L;
}

}  // namespace aal
