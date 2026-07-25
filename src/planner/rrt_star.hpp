// RRT* in R^3, with optional informed sampling and an optional rewiring
// switch that degrades it to plain RRT for side-by-side comparison.
//
// The planner works in the 3D position space rather than the full quadrotor
// state space. That is a deliberate choice, not a simplification: the quadrotor
// is differentially flat in (x, y, z, yaw), so every state and input can be
// recovered algebraically from a sufficiently smooth position trajectory. And
// the connection radius shrinks as (log n / n)^(1/d), which at d = 9 has fallen
// by only 2.7x after 100000 samples, meaning full-state kinodynamic RRT* is
// asymptotically optimal in name only for any budget we can actually run.
#pragma once

#include <cstdint>
#include <vector>

#include "collision/esdf.hpp"
#include "core/types.hpp"

namespace aal {

struct RrtStarParams {
  int max_iters = 5000;
  Scalar step = 2.0;           // eta, the maximum steer extension in metres
  Scalar goal_bias = 0.05;     // probability of sampling the goal directly
  Scalar goal_tolerance = 1.0; // metres, radius of the goal region
  Scalar robot_radius = 0.3;
  Scalar margin = 0.0;
  bool rewire = true;          // false gives plain RRT, for the comparison plot
  bool informed = true;        // restrict sampling to the informed set
  bool stop_at_first = false;  // stop once any solution is found
  std::uint64_t seed = 1;
  // Safety factor on the asymptotic-optimality lower bound for gamma. The proof
  // requires strict inequality, so a factor slightly above 1 is required, and
  // larger values trade planning speed for faster cost convergence.
  Scalar gamma_scale = 1.1;
  bool record_events = false;  // capture a replay log for the web visualiser
};

// One recorded planner event, enough to replay tree growth in a browser.
struct RrtEvent {
  int iter = 0;
  int node = -1;
  int parent = -1;
  Vec3 pos;
  Scalar cost = 0;
  Scalar best_cost = 0;
  Scalar radius = 0;
  // Edges whose parent changed during this iteration's rewiring pass, as
  // (child, new_parent) pairs. This is the part plain RRT never produces and
  // is what the visualisation highlights.
  std::vector<std::pair<int, int>> rewired;
};

struct CostSample {
  int iter;
  Scalar cost;
};

struct PlanResult {
  bool success = false;
  std::vector<Vec3> path;
  Scalar cost = 0;
  int iterations = 0;
  int num_nodes = 0;
  int rewires = 0;
  double wall_ms = 0;
  std::uint64_t collision_queries = 0;
  // Best cost as a function of iteration. Recorded only when it improves, so
  // this is the anytime convergence curve, which is the honest way to compare
  // planners that never terminate.
  std::vector<CostSample> history;
  std::vector<RrtEvent> events;
};

class RrtStar {
 public:
  RrtStar(const Esdf& field, RrtStarParams params);

  PlanResult plan(const Vec3& start, const Vec3& goal);

  // The connection radius at tree size n, exposed for tests and for the
  // write-up. r(n) = min( gamma * (log n / n)^(1/d), eta ).
  Scalar connectionRadius(int n, Scalar sampling_volume) const;

  // The lower bound gamma must exceed for asymptotic optimality:
  //   gamma > 2 * (1 + 1/d)^(1/d) * ( mu(X_free) / zeta_d )^(1/d)
  // with zeta_d the volume of the unit d-ball. Exposed so a test can assert the
  // value actually used sits above it.
  static Scalar gammaLowerBound(Scalar free_volume, int d = 3);

 private:
  struct Node {
    Vec3 pos;
    int parent = -1;
    Scalar cost = 0;
    std::vector<int> children;
  };

  Vec3 steer(const Vec3& from, const Vec3& to) const;
  void propagateCost(int node, Scalar delta);
  std::vector<Vec3> extractPath(int goal_node) const;

  const Esdf& field_;
  RrtStarParams p_;
  CollisionChecker checker_;
  std::vector<Node> nodes_;
};

}  // namespace aal
