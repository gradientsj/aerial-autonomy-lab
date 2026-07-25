// Benchmark harness for the planners.
//
// A randomised planner reported from a single run is not a measurement. Every
// number here is a mean over many seeds with a 95% confidence interval, and the
// RRT versus RRT* comparison is PAIRED: both planners see the same seed, so the
// per-seed difference removes the map-and-seed variance that otherwise swamps
// the effect. The output is JSON, consumed by the project write-up so that the
// published table cannot drift from what the code actually does.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <numeric>
#include <string>
#include <vector>

#include "collision/esdf.hpp"
#include "planner/rrt_star.hpp"

using namespace aal;

namespace {

struct Stats {
  double mean = 0, ci95 = 0, sd = 0;
  int n = 0;
};

Stats summarise(const std::vector<double>& v) {
  Stats s;
  s.n = static_cast<int>(v.size());
  if (s.n == 0) return s;
  s.mean = std::accumulate(v.begin(), v.end(), 0.0) / s.n;
  if (s.n < 2) return s;
  double var = 0;
  for (double x : v) var += (x - s.mean) * (x - s.mean);
  var /= (s.n - 1);
  s.sd = std::sqrt(var);
  // Normal approximation. With n >= 30 seeds the t correction is under 4%, and
  // reporting 1.96 keeps the interval comparable across rows.
  s.ci95 = 1.96 * s.sd / std::sqrt(double(s.n));
  return s;
}

struct Scenario {
  std::string name;
  Esdf field;
  Vec3 start, goal;
};

// A forest of vertical cylinders on a jittered grid. This is the canonical
// cluttered-but-solvable map for aerial planning benchmarks: obstacle density
// is tunable, and the optimal path is a non-trivial weave rather than a
// straight line.
Esdf forestMap(double density_spacing, std::uint64_t seed) {
  Esdf f(80, 80, 24, 0.5, Vec3{0, 0, 0});
  Rng rng(seed);
  const double extent = 40.0;
  for (double x = 4.0; x < extent - 2.0; x += density_spacing)
    for (double y = 2.0; y < extent - 2.0; y += density_spacing) {
      const double jx = x + rng.uniform(-0.8, 0.8);
      const double jy = y + rng.uniform(-0.8, 0.8);
      const double r = rng.uniform(0.4, 0.9);
      f.addCylinderZ(jx, jy, r, 0.0, 12.0);
    }
  f.build();
  return f;
}

// A wall with a single window. Forces one homotopy class through a narrow gap,
// which is where informed sampling should shine once a solution exists.
Esdf windowMap() {
  Esdf f(80, 60, 24, 0.5, Vec3{0, 0, 0});
  const int wall_i = static_cast<int>(20.0 / 0.5);
  for (int k = 0; k < 24; ++k)
    for (int j = 0; j < 60; ++j) {
      const double y = 0.5 * j, z = 0.5 * k;
      const bool window = std::abs(y - 15.0) < 1.6 && std::abs(z - 6.0) < 1.6;
      if (!window) f.setOccupied(wall_i, j, k, true);
    }
  f.build();
  return f;
}

Esdf emptyMap() {
  Esdf f(80, 80, 24, 0.5, Vec3{0, 0, 0});
  f.build();
  return f;
}

struct Variant {
  std::string name;
  bool rewire;
  bool informed;
};

std::string jsonStat(const std::string& key, const Stats& s) {
  char buf[256];
  std::snprintf(buf, sizeof(buf), "\"%s\": {\"mean\": %.4f, \"ci95\": %.4f, \"n\": %d}",
                key.c_str(), s.mean, s.ci95, s.n);
  return buf;
}

}  // namespace

int main(int argc, char** argv) {
  const int seeds = (argc > 1) ? std::atoi(argv[1]) : 40;
  const int iters = (argc > 2) ? std::atoi(argv[2]) : 4000;
  const std::string out_path = (argc > 3) ? argv[3] : "bench/results/planners.json";

  std::vector<Scenario> scenarios;
  scenarios.push_back({"empty", emptyMap(), Vec3{2, 2, 4}, Vec3{36, 36, 8}});
  scenarios.push_back({"forest_sparse", forestMap(6.0, 11), Vec3{2, 2, 5}, Vec3{36, 36, 5}});
  scenarios.push_back({"forest_dense", forestMap(3.5, 12), Vec3{2, 2, 5}, Vec3{36, 36, 5}});
  scenarios.push_back({"window", windowMap(), Vec3{4, 15, 6}, Vec3{34, 15, 6}});

  const std::vector<Variant> variants = {
      {"RRT", false, false},
      {"RRT*", true, false},
      {"Informed RRT*", true, true},
  };

  std::string json = "{\n  \"config\": {\"seeds\": " + std::to_string(seeds) +
                     ", \"max_iters\": " + std::to_string(iters) + "},\n  \"scenarios\": [\n";

  for (std::size_t si = 0; si < scenarios.size(); ++si) {
    Scenario& sc = scenarios[si];
    std::printf("\n=== %s (free fraction %.3f) ===\n", sc.name.c_str(), sc.field.freeFraction());
    json += "    {\"name\": \"" + sc.name + "\", \"free_fraction\": " +
            std::to_string(sc.field.freeFraction()) + ", \"variants\": [\n";

    // Paired storage so RRT* can be compared against RRT on identical seeds.
    std::vector<std::vector<double>> cost_by_variant(variants.size());

    for (std::size_t vi = 0; vi < variants.size(); ++vi) {
      const Variant& v = variants[vi];
      std::vector<double> costs, times, nodes, queries;
      int successes = 0;

      for (int s = 0; s < seeds; ++s) {
        RrtStarParams p;
        p.seed = 1000u + s;
        p.max_iters = iters;
        p.step = 2.0;
        p.goal_tolerance = 1.0;
        p.robot_radius = 0.35;
        p.rewire = v.rewire;
        p.informed = v.informed;

        RrtStar planner(sc.field, p);
        const PlanResult r = planner.plan(sc.start, sc.goal);
        if (r.success) {
          ++successes;
          costs.push_back(r.cost);
          cost_by_variant[vi].push_back(r.cost);
        } else {
          cost_by_variant[vi].push_back(std::nan(""));
        }
        times.push_back(r.wall_ms);
        nodes.push_back(r.num_nodes);
        queries.push_back(double(r.collision_queries));
      }

      const Stats c = summarise(costs), t = summarise(times), n = summarise(nodes),
                  q = summarise(queries);
      const double sr = 100.0 * successes / seeds;
      std::printf("  %-16s success %5.1f%%  cost %7.2f +- %4.2f  time %7.1f ms  nodes %6.0f\n",
                  v.name.c_str(), sr, c.mean, c.ci95, t.mean, n.mean);

      json += "      {\"planner\": \"" + v.name + "\", \"success_pct\": " + std::to_string(sr) +
              ", " + jsonStat("cost", c) + ", " + jsonStat("plan_ms", t) + ", " +
              jsonStat("nodes", n) + ", " + jsonStat("edge_queries", q) + "}";
      json += (vi + 1 < variants.size()) ? ",\n" : "\n";
    }

    // Paired comparison of RRT* against RRT, using only seeds where both found
    // a path. Reporting the paired difference rather than the difference of
    // means is what makes the claim defensible.
    std::vector<double> paired;
    for (int s = 0; s < seeds; ++s) {
      const double a = cost_by_variant[0][s], b = cost_by_variant[1][s];
      if (!std::isnan(a) && !std::isnan(b)) paired.push_back(a - b);
    }
    const Stats d = summarise(paired);
    const double t_stat = (d.n > 1 && d.sd > 0) ? d.mean / (d.sd / std::sqrt(double(d.n))) : 0.0;
    std::printf("  paired RRT - RRT*: %+.3f m (95%% CI +-%.3f, n=%d, t=%.2f)\n", d.mean, d.ci95,
                d.n, t_stat);

    json += "    ], \"paired_rrt_minus_rrtstar\": {\"mean\": " + std::to_string(d.mean) +
            ", \"ci95\": " + std::to_string(d.ci95) + ", \"n\": " + std::to_string(d.n) +
            ", \"t\": " + std::to_string(t_stat) + "}}";
    json += (si + 1 < scenarios.size()) ? ",\n" : "\n";
  }
  json += "  ]\n}\n";

  std::ofstream out(out_path);
  if (out) {
    out << json;
    std::printf("\nwrote %s\n", out_path.c_str());
  } else {
    std::printf("\ncould not open %s for writing\n", out_path.c_str());
  }
  return 0;
}
