// Tests for RRT*.
//
// A randomised planner cannot be tested by running it once and eyeballing the
// answer. These tests fall into three groups:
//   1. Deterministic invariants that must hold on every single run.
//   2. Analytic checks of the radius formula, which is where implementations
//      most often quietly break asymptotic optimality.
//   3. Statistical claims over many seeds, stated with a margin, because any
//      single-seed comparison between RRT and RRT* proves nothing.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include "collision/esdf.hpp"
#include "planner/rrt_star.hpp"

using namespace aal;

namespace {

// An empty box map, where the optimal path is exactly the straight line and we
// can therefore compare against a known optimum rather than against ourselves.
Esdf emptyMap() {
  Esdf f(40, 40, 20, 0.5, Vec3{0, 0, 0});
  f.build();
  return f;
}

// A wall with a single window, which forces a detour and gives RRT* something
// to improve on.
Esdf wallWithWindow() {
  Esdf f(60, 40, 24, 0.5, Vec3{0, 0, 0});
  const Scalar wall_x = 15.0;
  for (int k = 0; k < 24; ++k)
    for (int j = 0; j < 40; ++j) {
      const Vec3 p = Vec3{wall_x, 0.5 * j, 0.5 * k};
      // Window centred at (y, z) = (10, 6), 3 m across.
      const bool window = std::abs(p.y - 10.0) < 1.5 && std::abs(p.z - 6.0) < 1.5;
      if (!window) {
        const int i = static_cast<int>(wall_x / 0.5);
        f.setOccupied(i, j, k, true);
      }
    }
  f.build();
  return f;
}

Scalar pathLength(const std::vector<Vec3>& path) {
  Scalar s = 0;
  for (std::size_t i = 1; i < path.size(); ++i) s += distance(path[i - 1], path[i]);
  return s;
}

}  // namespace

TEST(RrtStar, GammaExceedsTheAsymptoticOptimalityLowerBound) {
  // gamma > 2 (1 + 1/d)^(1/d) (mu_free / zeta_d)^(1/d).
  // For d = 3 and mu_free = 280000 m^3 the bound evaluates to 89.34, a number
  // worth pinning down because an off-by-one in zeta_d or a missing factor of 2
  // produces a plausible-looking radius that silently forfeits optimality.
  const Scalar bound = RrtStar::gammaLowerBound(280000.0, 3);
  EXPECT_NEAR(bound, 89.337, 0.01);

  // Reproduce it independently from the definition.
  const Scalar zeta3 = 4.0 / 3.0 * M_PI;
  const Scalar expected = 2.0 * std::pow(1.0 + 1.0 / 3.0, 1.0 / 3.0) * std::cbrt(280000.0 / zeta3);
  EXPECT_NEAR(bound, expected, 1e-9);
}

TEST(RrtStar, ConnectionRadiusShrinksLikeLogNOverN) {
  Esdf f = emptyMap();
  RrtStarParams p;
  p.step = 1e9;  // remove the eta cap so the formula itself is under test
  RrtStar planner(f, p);
  const Scalar vol = 280000.0;

  const Scalar r1e3 = planner.connectionRadius(1000, vol);
  const Scalar r1e4 = planner.connectionRadius(10000, vol);
  const Scalar r1e5 = planner.connectionRadius(100000, vol);
  EXPECT_GT(r1e3, r1e4);
  EXPECT_GT(r1e4, r1e5);

  // The ratio between decades is ( (log n1/n1) / (log n2/n2) )^(1/3), which for
  // 10^3 to 10^4 is (0.006908 / 0.000921)^(1/3) = 1.957.
  EXPECT_NEAR(r1e3 / r1e4, 1.957, 0.01);
}

TEST(RrtStar, RadiusIsCappedByStepSize) {
  Esdf f = emptyMap();
  RrtStarParams p;
  p.step = 2.0;
  RrtStar planner(f, p);
  // At small n the unbounded formula is huge, so the cap must bite.
  EXPECT_DOUBLE_EQ(planner.connectionRadius(5, 280000.0), 2.0);
}

TEST(RrtStar, FindsAPathAndTheReportedCostMatchesTheGeometry) {
  Esdf f = emptyMap();
  RrtStarParams p;
  p.max_iters = 3000;
  p.step = 2.0;
  p.seed = 1;
  RrtStar planner(f, p);

  const Vec3 start{1, 1, 1}, goal{18, 18, 8};
  const PlanResult r = planner.plan(start, goal);

  ASSERT_TRUE(r.success);
  ASSERT_GE(r.path.size(), 2u);
  EXPECT_NEAR(r.path.front().x, start.x, 1e-9);
  EXPECT_NEAR(r.path.back().x, goal.x, 1e-9);
  EXPECT_NEAR(r.path.back().y, goal.y, 1e-9);
  EXPECT_NEAR(r.path.back().z, goal.z, 1e-9);

  // The reported cost must equal the actual geometric length of the returned
  // polyline. These drift apart when rewiring updates a cost but not the parent
  // pointers, which is the classic RRT* bug and is invisible to the eye.
  EXPECT_NEAR(r.cost, pathLength(r.path), 1e-6);
}

TEST(RrtStar, ConvergesTowardTheKnownOptimumInAnEmptyMap) {
  // In an empty map the optimum is the straight line, so we can measure the
  // optimality gap against a true reference rather than a self-comparison.
  Esdf f = emptyMap();
  const Vec3 start{1, 1, 1}, goal{18, 18, 8};
  const Scalar optimal = distance(start, goal);

  Scalar gap_small = 0, gap_large = 0;
  constexpr int kSeeds = 12;
  for (int s = 0; s < kSeeds; ++s) {
    RrtStarParams p;
    p.seed = 100 + s;
    p.step = 2.0;
    p.goal_tolerance = 0.5;

    p.max_iters = 400;
    gap_small += RrtStar(f, p).plan(start, goal).cost / optimal;

    p.max_iters = 6000;
    gap_large += RrtStar(f, p).plan(start, goal).cost / optimal;
  }
  gap_small /= kSeeds;
  gap_large /= kSeeds;

  // More samples must not make the solution worse, and the converged solution
  // should be close to optimal in a map with no obstacles at all.
  EXPECT_LT(gap_large, gap_small);
  EXPECT_LT(gap_large, 1.06) << "converged cost is " << gap_large << "x optimal";
  EXPECT_GE(gap_large, 1.0 - 1e-9) << "reported a path shorter than the straight line";
}

TEST(RrtStar, ReturnedPathIsAlwaysCollisionFree) {
  Esdf f = wallWithWindow();
  const Vec3 start{2, 10, 6}, goal{25, 10, 6};

  for (int s = 0; s < 20; ++s) {
    RrtStarParams p;
    p.seed = s;
    p.max_iters = 4000;
    p.step = 1.5;
    p.robot_radius = 0.3;
    RrtStar planner(f, p);
    const PlanResult r = planner.plan(start, goal);
    if (!r.success) continue;

    CollisionChecker cc(f, p.robot_radius);
    for (std::size_t i = 1; i < r.path.size(); ++i) {
      EXPECT_TRUE(cc.edgeFree(r.path[i - 1], r.path[i]))
          << "seed " << s << " segment " << i << " passes through an obstacle";
    }
  }
}

TEST(RrtStar, IsBitwiseReproducibleForAFixedSeed) {
  Esdf f = wallWithWindow();
  const Vec3 start{2, 10, 6}, goal{25, 10, 6};
  RrtStarParams p;
  p.seed = 4242;
  p.max_iters = 2000;

  const PlanResult a = RrtStar(f, p).plan(start, goal);
  const PlanResult b = RrtStar(f, p).plan(start, goal);

  EXPECT_EQ(a.success, b.success);
  EXPECT_EQ(a.num_nodes, b.num_nodes);
  EXPECT_EQ(a.rewires, b.rewires);
  EXPECT_EQ(a.collision_queries, b.collision_queries);
  ASSERT_EQ(a.path.size(), b.path.size());
  for (std::size_t i = 0; i < a.path.size(); ++i) {
    EXPECT_DOUBLE_EQ(a.path[i].x, b.path[i].x);
    EXPECT_DOUBLE_EQ(a.path[i].y, b.path[i].y);
    EXPECT_DOUBLE_EQ(a.path[i].z, b.path[i].z);
  }
  EXPECT_DOUBLE_EQ(a.cost, b.cost);
}

TEST(RrtStar, RewiringBeatsPlainRrtOverManySeeds) {
  // The headline claim of the algorithm. Stated over 30 paired seeds rather
  // than one run, because the per-seed variance of RRT is large enough that a
  // single comparison can point either way.
  Esdf f = wallWithWindow();
  const Vec3 start{2, 10, 6}, goal{25, 10, 6};
  constexpr int kSeeds = 30;

  std::vector<Scalar> diffs;
  for (int s = 0; s < kSeeds; ++s) {
    RrtStarParams base;
    base.seed = 900 + s;
    base.max_iters = 3000;
    base.step = 1.5;
    base.informed = false;  // isolate the effect of rewiring alone

    RrtStarParams rrt = base;
    rrt.rewire = false;
    RrtStarParams star = base;
    star.rewire = true;

    const PlanResult a = RrtStar(f, rrt).plan(start, goal);
    const PlanResult b = RrtStar(f, star).plan(start, goal);
    if (a.success && b.success) diffs.push_back(a.cost - b.cost);
  }

  ASSERT_GE(diffs.size(), 20u) << "too few successes to make a claim";
  const Scalar mean = std::accumulate(diffs.begin(), diffs.end(), Scalar(0)) / diffs.size();
  Scalar var = 0;
  for (Scalar d : diffs) var += (d - mean) * (d - mean);
  var /= (diffs.size() - 1);
  const Scalar se = std::sqrt(var / diffs.size());

  // Paired t statistic against the null hypothesis that rewiring does nothing.
  const Scalar t = mean / (se > 0 ? se : 1e-12);
  EXPECT_GT(mean, 0.0) << "RRT* was not cheaper than RRT on average";
  EXPECT_GT(t, 3.0) << "improvement is not statistically significant, t = " << t;
}

TEST(RrtStar, InformedSamplingDoesNotHurtAndUsuallyHelps) {
  Esdf f = wallWithWindow();
  const Vec3 start{2, 10, 6}, goal{25, 10, 6};
  constexpr int kSeeds = 20;

  Scalar sum_plain = 0, sum_informed = 0;
  int n = 0;
  for (int s = 0; s < kSeeds; ++s) {
    RrtStarParams base;
    base.seed = 500 + s;
    base.max_iters = 4000;
    base.step = 1.5;

    RrtStarParams plain = base;
    plain.informed = false;
    RrtStarParams inf = base;
    inf.informed = true;

    const PlanResult a = RrtStar(f, plain).plan(start, goal);
    const PlanResult b = RrtStar(f, inf).plan(start, goal);
    if (a.success && b.success) {
      sum_plain += a.cost;
      sum_informed += b.cost;
      ++n;
    }
  }
  ASSERT_GT(n, 12);
  // Informed sampling only ever restricts sampling to states that could improve
  // the incumbent, so at equal iteration count it must not be worse by more
  // than noise.
  EXPECT_LE(sum_informed / n, sum_plain / n * 1.02);
}

TEST(RrtStar, ReportsFailureRatherThanACollidingPathWhenBoxedIn) {
  // Start sealed inside a solid box. There is no path, and the planner must say
  // so rather than returning something that looks like one.
  Esdf f(40, 40, 20, 0.5, Vec3{0, 0, 0});
  f.addBox(Aabb{Vec3{4, 4, 2}, Vec3{9, 9, 7}});
  // Carve a cavity, leaving a sealed shell around the start point.
  for (int k = 0; k < 20; ++k)
    for (int j = 0; j < 40; ++j)
      for (int i = 0; i < 40; ++i) {
        const Vec3 p = Vec3{0.5 * i, 0.5 * j, 0.5 * k};
        if (p.x > 5 && p.x < 8 && p.y > 5 && p.y < 8 && p.z > 3 && p.z < 6)
          f.setOccupied(i, j, k, false);
      }
  f.build();

  RrtStarParams p;
  p.max_iters = 3000;
  p.seed = 8;
  p.robot_radius = 0.2;
  RrtStar planner(f, p);
  const PlanResult r = planner.plan(Vec3{6.5, 6.5, 4.5}, Vec3{18, 18, 8});
  EXPECT_FALSE(r.success);
  EXPECT_TRUE(r.path.empty());
}

TEST(RrtStar, EventLogReplaysTheSameTreeItBuilt) {
  // The web visualisation replays this log, so it has to be faithful: every
  // event's parent must already exist, and node ids must be dense and ordered.
  Esdf f = wallWithWindow();
  RrtStarParams p;
  p.seed = 3;
  p.max_iters = 800;
  p.record_events = true;
  RrtStar planner(f, p);
  const PlanResult r = planner.plan(Vec3{2, 10, 6}, Vec3{25, 10, 6});

  ASSERT_FALSE(r.events.empty());
  int expected_id = 1;  // node 0 is the root, which is never an event
  for (const RrtEvent& e : r.events) {
    EXPECT_EQ(e.node, expected_id++);
    EXPECT_GE(e.parent, 0);
    EXPECT_LT(e.parent, e.node) << "an event references a parent that does not exist yet";
    for (const auto& [child, new_parent] : e.rewired) {
      EXPECT_LT(child, e.node);
      EXPECT_EQ(new_parent, e.node);
    }
  }
}
