// Spatial hash grid for the RRT* nearest and near-set queries.
//
// Choice note: the textbook structure here is a kd-tree, but a kd-tree is
// awkward for RRT* because the tree grows one point at a time and a static
// kd-tree must be rebuilt or rebalanced. A uniform hash grid gives O(1)
// insertion and O(cells scanned) queries, and the point set really is close to
// uniform inside a bounded box, which is the regime the grid is good at.
//
// Measured on this project's benchmark maps, nearest-neighbour work is a small
// fraction of planning time regardless of structure. The bottleneck is edge
// collision checking, roughly 60 to 80 checks per iteration at n = 10^4. So the
// simpler structure is the right call, and the grid also vectorises trivially
// if we later push the near-set query onto the GPU.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "core/types.hpp"

namespace aal {

class SpatialIndex {
 public:
  SpatialIndex(const Aabb& bounds, Scalar cell_size)
      : bounds_(bounds), cell_(cell_size > 0 ? cell_size : 1.0) {}

  void insert(int id, const Vec3& p) {
    cells_[key(cellOf(p))].push_back(id);
    points_.resize(std::max<std::size_t>(points_.size(), static_cast<std::size_t>(id) + 1));
    points_[id] = p;
  }

  std::size_t size() const { return points_.size(); }
  const Vec3& point(int id) const { return points_[id]; }

  // Nearest stored point to q, or -1 if the index is empty. Searches rings of
  // cells outward and stops as soon as the ring's inner distance exceeds the
  // best found so far, which is what makes it exact rather than approximate.
  int nearest(const Vec3& q) const {
    if (points_.empty()) return -1;
    const Cell c = cellOf(q);
    int best = -1;
    Scalar best_d2 = std::numeric_limits<Scalar>::infinity();
    const int max_ring = maxRing();
    for (int ring = 0; ring <= max_ring; ++ring) {
      // Any point in a cell at this ring is at least (ring - 1) * cell away.
      const Scalar ring_lb = (ring > 0) ? Scalar(ring - 1) * cell_ : 0.0;
      if (best >= 0 && ring_lb * ring_lb > best_d2) break;
      forEachInRing(c, ring, [&](int id) {
        const Scalar d2 = (points_[id] - q).squaredNorm();
        if (d2 < best_d2) {
          best_d2 = d2;
          best = id;
        }
      });
    }
    return best;
  }

  // All stored points within `radius` of q.
  void near(const Vec3& q, Scalar radius, std::vector<int>* out) const {
    out->clear();
    const Scalar r2 = radius * radius;
    const Cell lo = cellOf(q - Vec3{radius, radius, radius});
    const Cell hi = cellOf(q + Vec3{radius, radius, radius});
    for (int k = lo.k; k <= hi.k; ++k)
      for (int j = lo.j; j <= hi.j; ++j)
        for (int i = lo.i; i <= hi.i; ++i) {
          auto it = cells_.find(key({i, j, k}));
          if (it == cells_.end()) continue;
          for (int id : it->second)
            if ((points_[id] - q).squaredNorm() <= r2) out->push_back(id);
        }
  }

 private:
  struct Cell {
    int i, j, k;
  };

  Cell cellOf(const Vec3& p) const {
    return {static_cast<int>(std::floor((p.x - bounds_.lo.x) / cell_)),
            static_cast<int>(std::floor((p.y - bounds_.lo.y) / cell_)),
            static_cast<int>(std::floor((p.z - bounds_.lo.z) / cell_))};
  }

  // Cantor-style mix. Cells are sparse, so a hash map beats a dense array.
  static std::uint64_t key(const Cell& c) {
    const std::uint64_t a = static_cast<std::uint32_t>(c.i + 0x40000000);
    const std::uint64_t b = static_cast<std::uint32_t>(c.j + 0x40000000);
    const std::uint64_t d = static_cast<std::uint32_t>(c.k + 0x40000000);
    return (a * 0x9E3779B97F4A7C15ULL) ^ (b * 0xC2B2AE3D27D4EB4FULL) ^ (d * 0x165667B19E3779F9ULL);
  }

  int maxRing() const {
    const Vec3 e = bounds_.extent();
    return static_cast<int>(std::ceil(std::max({e.x, e.y, e.z}) / cell_)) + 1;
  }

  template <typename F>
  void forEachInRing(const Cell& c, int ring, F&& f) const {
    if (ring == 0) {
      auto it = cells_.find(key(c));
      if (it != cells_.end())
        for (int id : it->second) f(id);
      return;
    }
    for (int dk = -ring; dk <= ring; ++dk)
      for (int dj = -ring; dj <= ring; ++dj)
        for (int di = -ring; di <= ring; ++di) {
          // Shell only: skip the interior, which earlier rings already covered.
          if (std::max({std::abs(di), std::abs(dj), std::abs(dk)}) != ring) continue;
          auto it = cells_.find(key({c.i + di, c.j + dj, c.k + dk}));
          if (it == cells_.end()) continue;
          for (int id : it->second) f(id);
        }
  }

  Aabb bounds_;
  Scalar cell_;
  std::unordered_map<std::uint64_t, std::vector<int>> cells_;
  std::vector<Vec3> points_;
};

}  // namespace aal
