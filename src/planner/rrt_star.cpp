#include "planner/rrt_star.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>

#include "planner/informed_sampler.hpp"
#include "planner/spatial_index.hpp"

namespace aal {
namespace {
constexpr Scalar kInf = std::numeric_limits<Scalar>::infinity();
// Volume of the unit 3-ball, zeta_3 = pi^(3/2) / Gamma(5/2) = 4 pi / 3.
constexpr Scalar kZeta3 = 4.18879020478639098462;
}  // namespace

RrtStar::RrtStar(const Esdf& field, RrtStarParams params)
    : field_(field), p_(params), checker_(field, params.robot_radius, params.margin) {}

Scalar RrtStar::gammaLowerBound(Scalar free_volume, int d) {
  const Scalar zeta_d = (d == 3) ? kZeta3
                                 : std::pow(M_PI, d / 2.0) / std::tgamma(d / 2.0 + 1.0);
  return 2.0 * std::pow(1.0 + 1.0 / d, 1.0 / d) * std::pow(free_volume / zeta_d, 1.0 / d);
}

Scalar RrtStar::connectionRadius(int n, Scalar sampling_volume) const {
  if (n < 2) return p_.step;
  const Scalar gamma = p_.gamma_scale * gammaLowerBound(sampling_volume, 3);
  const Scalar r = gamma * std::pow(std::log(Scalar(n)) / Scalar(n), 1.0 / 3.0);
  return std::min(r, p_.step);
}

Vec3 RrtStar::steer(const Vec3& from, const Vec3& to) const {
  const Vec3 d = to - from;
  const Scalar n = d.norm();
  if (n <= p_.step) return to;
  return from + d / n * p_.step;
}

void RrtStar::propagateCost(int node, Scalar delta) {
  // Iterative rather than recursive: rewiring near the root can touch a large
  // subtree and deep recursion has blown the stack in implementations of this
  // exact algorithm.
  static thread_local std::vector<int> stack;
  stack.clear();
  stack.push_back(node);
  while (!stack.empty()) {
    const int u = stack.back();
    stack.pop_back();
    for (int c : nodes_[u].children) {
      nodes_[c].cost += delta;
      stack.push_back(c);
    }
  }
}

std::vector<Vec3> RrtStar::extractPath(int goal_node) const {
  std::vector<Vec3> path;
  for (int u = goal_node; u != -1; u = nodes_[u].parent) path.push_back(nodes_[u].pos);
  std::reverse(path.begin(), path.end());
  return path;
}

PlanResult RrtStar::plan(const Vec3& start, const Vec3& goal) {
  using Clock = std::chrono::steady_clock;
  const auto t0 = Clock::now();

  PlanResult res;
  nodes_.clear();
  checker_.resetQueryCount();

  const Aabb bounds = field_.bounds();
  // mu(X_free) is not known a priori in an unmapped world. We estimate it from
  // the occupancy grid's free fraction. Underestimating it would shrink the
  // radius below the asymptotic-optimality threshold and silently break the
  // guarantee, so the bounding-box volume is kept as a floor via gamma_scale.
  const Scalar free_volume = bounds.volume() * std::max<Scalar>(field_.freeFraction(), 1e-3);

  Rng rng(p_.seed);
  InformedSampler informed(start, goal);
  SpatialIndex index(bounds, std::max(p_.step, field_.resolution()));

  if (!checker_.pointFree(start)) return res;  // start in collision, nothing to do

  nodes_.push_back(Node{start, -1, 0.0, {}});
  index.insert(0, start);

  int best_goal_node = -1;
  Scalar best_cost = kInf;
  std::vector<int> near;

  for (int it = 0; it < p_.max_iters; ++it) {
    res.iterations = it + 1;

    // --- sample ---------------------------------------------------------
    Vec3 x_rand;
    const bool have_solution = best_goal_node >= 0;
    if (rng.uniform() < p_.goal_bias) {
      x_rand = goal;
    } else if (p_.informed && have_solution && best_cost > informed.cMin()) {
      // Rejection against the world bounds keeps samples legal; the ellipsoid
      // itself is sampled directly, so this rejects only the part of the
      // ellipsoid poking outside the map.
      for (int tries = 0; tries < 16; ++tries) {
        x_rand = informed.sample(rng, best_cost);
        if (bounds.contains(x_rand)) break;
      }
    } else {
      x_rand = rng.uniformInAabb(bounds);
    }

    // --- extend ---------------------------------------------------------
    const int nearest = index.nearest(x_rand);
    if (nearest < 0) continue;
    const Vec3 x_new = steer(nodes_[nearest].pos, x_rand);
    if (!bounds.contains(x_new)) continue;
    if (!checker_.edgeFree(nodes_[nearest].pos, x_new)) continue;

    // --- near set -------------------------------------------------------
    // Once informed sampling is active the sampled measure is the ellipsoid,
    // not all of X_free, so the radius formula must use that smaller volume or
    // it stays needlessly large.
    Scalar sampling_volume = free_volume;
    if (p_.informed && have_solution) {
      const Scalar v = informed.volume(best_cost);
      if (v > 0) sampling_volume = std::min(sampling_volume, v);
    }
    const Scalar radius = connectionRadius(static_cast<int>(nodes_.size()) + 1, sampling_volume);

    int parent = nearest;
    Scalar best_new_cost = nodes_[nearest].cost + distance(nodes_[nearest].pos, x_new);

    if (p_.rewire) {
      index.near(x_new, radius, &near);
      // Choose the cheapest feasible parent among the near set.
      for (int cand : near) {
        const Scalar c = nodes_[cand].cost + distance(nodes_[cand].pos, x_new);
        if (c < best_new_cost && checker_.edgeFree(nodes_[cand].pos, x_new)) {
          best_new_cost = c;
          parent = cand;
        }
      }
    }

    const int new_id = static_cast<int>(nodes_.size());
    nodes_.push_back(Node{x_new, parent, best_new_cost, {}});
    nodes_[parent].children.push_back(new_id);
    index.insert(new_id, x_new);

    RrtEvent ev;
    if (p_.record_events) {
      ev.iter = it;
      ev.node = new_id;
      ev.parent = parent;
      ev.pos = x_new;
      ev.cost = best_new_cost;
      ev.radius = radius;
    }

    // --- rewire ---------------------------------------------------------
    if (p_.rewire) {
      for (int cand : near) {
        if (cand == parent || cand == new_id) continue;
        const Scalar through_new = best_new_cost + distance(x_new, nodes_[cand].pos);
        if (through_new < nodes_[cand].cost && checker_.edgeFree(x_new, nodes_[cand].pos)) {
          // Detach from the old parent.
          const int old_parent = nodes_[cand].parent;
          if (old_parent >= 0) {
            auto& sib = nodes_[old_parent].children;
            sib.erase(std::remove(sib.begin(), sib.end(), cand), sib.end());
          }
          const Scalar delta = through_new - nodes_[cand].cost;
          nodes_[cand].parent = new_id;
          nodes_[cand].cost = through_new;
          nodes_[new_id].children.push_back(cand);
          // Every descendant's cost shifts by the same amount, because cost is
          // additive along the path and only this one edge changed.
          propagateCost(cand, delta);
          ++res.rewires;
          if (p_.record_events) ev.rewired.emplace_back(cand, new_id);
        }
      }
    }

    // --- goal check -----------------------------------------------------
    const Scalar d_goal = distance(x_new, goal);
    if (d_goal <= p_.goal_tolerance && checker_.edgeFree(x_new, goal)) {
      const Scalar total = nodes_[new_id].cost + d_goal;
      if (total < best_cost) {
        best_cost = total;
        best_goal_node = new_id;
        res.history.push_back({it, total});
      }
    } else if (best_goal_node >= 0) {
      // Rewiring may have lowered the cost of the node the incumbent path goes
      // through, so the incumbent cost has to be re-read rather than cached.
      const Scalar refreshed = nodes_[best_goal_node].cost + distance(nodes_[best_goal_node].pos, goal);
      if (refreshed < best_cost) {
        best_cost = refreshed;
        res.history.push_back({it, refreshed});
      }
    }

    if (p_.record_events) {
      ev.best_cost = std::isfinite(best_cost) ? best_cost : 0.0;
      res.events.push_back(std::move(ev));
    }

    if (p_.stop_at_first && best_goal_node >= 0) break;
  }

  if (best_goal_node >= 0) {
    res.success = true;
    res.path = extractPath(best_goal_node);
    res.path.push_back(goal);
    res.cost = best_cost;
  }
  res.num_nodes = static_cast<int>(nodes_.size());
  res.collision_queries = checker_.queryCount();
  res.wall_ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
  return res;
}

}  // namespace aal
