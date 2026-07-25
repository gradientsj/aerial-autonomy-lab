# aerial-autonomy-lab

Motion planning and autonomy for a quadrotor in **unknown, GPS-denied
environments**. C++17 core, CUDA sensing, Python orchestration, and a browser
based visualisation of the algorithms.

Write-up: [stanleyjacob.dev/robotics/projects/aerial-autonomy-lab](https://stanleyjacob.dev/robotics/projects/aerial-autonomy-lab)

---

## What this is

A drone with no GPS fix has to build its own map, localise inside it, and plan
through it, all while the map is still being discovered. This repository builds
that stack from the bottom up and measures each piece rather than asserting it.

Currently implemented:

- **Exact Euclidean signed distance field.** Felzenszwalb's O(n) separable
  distance transform over a voxel grid, signed inside obstacles, with trilinear
  interpolation and a closed-form gradient. Tested for exact agreement with a
  brute-force transform, not approximate agreement.
- **Sphere-traced collision checking.** Edge checks stride by the certified free
  radius instead of sampling at fixed resolution, so open space costs almost
  nothing and only near-contact geometry costs anything.
- **RRT\*, Informed RRT\*, and plain RRT** sharing one implementation, so the
  comparison between them is controlled rather than between two codebases.
- **A paired-seed benchmark harness** reporting means with 95% confidence
  intervals and a paired t statistic.

Planned, in order: GPU batch edge checking, BIT\*, safe-flight-corridor
minimum-snap smoothing, the quadrotor dynamics and SE(3) controller, a CUDA
raycast depth sensor, an error-state Kalman filter for visual-inertial odometry,
and a learned sampling distribution for the planner.

## Two decisions worth explaining

**Plan in R³, not in the full state space.** The RRT\* connection radius shrinks
as `(log n / n)^(1/d)`. At `d = 3` and `n = 10^5` that factor is 0.049; at
`d = 9` it is 0.365, meaning the radius has shrunk by only 2.7x after a hundred
thousand samples and the near-set is effectively the entire tree. Full-state
kinodynamic RRT\* is asymptotically optimal in name only at any budget that
actually runs. The quadrotor is differentially flat in `(x, y, z, yaw)`, so
every state and input is an algebraic function of a smooth position trajectory.
Planning geometrically and recovering the dynamics through flatness is exact,
not a simplification.

**No game engine.** The development machine is a headless 2xH100 node. Hopper
has no RT cores and no display engine, and the box carries no EGL, Vulkan or
OpenGL userspace at all, so Isaac Sim, Unreal and Unity camera output are all
unavailable. This turns out not to matter: depth and range sensing are
raycasting problems rather than rasterisation problems, and a hand-written CUDA
sphere-tracer renders 1920x1080 in 0.025 ms on one H100. Visualisation lives in
the browser, where the algorithms are the thing worth showing anyway.

## Build

Requires a C++17 compiler, CMake 3.20 or newer, Eigen, and GoogleTest.

```bash
sudo apt install -y libeigen3-dev libgtest-dev libgmock-dev nlohmann-json3-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/aal_tests
./build/bench_planners 40 4000 bench/results/planners.json
```

The build pins `-ffp-contract=off`. GCC defaults to contracting `a*b+c` into a
fused multiply-add even without `-ffast-math`, which changes the low bits of
every distance computation and would break golden-trajectory comparisons
between Linux and Windows.

## Measured results

40 seeds per cell, 4000 iterations, 0.35 m robot radius, on a 40 x 40 x 12 m
world. Cost is path length in metres, plus or minus a 95% confidence interval.

| map | RRT | RRT\* | Informed RRT\* | paired RRT − RRT\* |
| --- | --- | --- | --- | --- |
| empty | 61.35 ± 1.12 | 57.52 ± 0.64 | 57.19 ± 0.83 | +3.83 ± 0.73 (t = 10.3) |
| forest, sparse | 61.63 ± 1.10 | 58.63 ± 0.83 | 57.34 ± 0.65 | +3.00 ± 0.62 (t = 9.6) |
| forest, dense | 66.58 ± 1.48 | 62.60 ± 1.06 | 63.15 ± 1.32 | +3.98 ± 0.88 (t = 8.9) |
| wall with window | 39.15 ± 1.02 | 35.22 ± 0.48 | 32.05 ± 0.32 | +3.93 ± 0.94 (t = 8.2) |

Rewiring beats plain RRT everywhere, significantly. Informed sampling is a large
win on the window map, where the informed ellipsoid collapses tightly around a
forced detour, and **no better than plain RRT\* on the dense forest**, where the
confidence intervals overlap and the point estimate is slightly worse. That
negative result is reported rather than buried: when the first solution is
already close to the straight-line lower bound, the informed set barely shrinks
and the extra sampling machinery only adds overhead.

## Testing

```bash
./build/aal_tests
```

The suite checks properties rather than outputs where it can. The distance field
is asserted equal to brute force on every voxel of several occupancy patterns.
The collision checker is asserted sound against dense sampling over thousands of
random edges. The planner is asserted to converge toward the analytically known
optimum in an empty map, to be bitwise reproducible for a fixed seed, and to
never return a path that intersects an obstacle. The claim that rewiring helps
is a paired t-test over 30 seeds, not a single run.

Two bugs this suite caught, both of which produce code that looks correct:

1. Seeding the distance transform with infinity makes the parabola intersection
   evaluate `INF - INF = NaN` on any scanline containing no obstacle.
2. The sphere trace originally strode by the full slack `d - clearance`, which
   is sound for the exact distance function because it is 1-Lipschitz. But the
   queried object is the *trilinear interpolant*, whose partial derivatives are
   each bounded by 1 and whose gradient norm therefore reaches √3 along a
   diagonal. The full stride steps over violations. Strides are now divided by
   √3, and `edgeFree` canonicalises endpoint order so it is a true function of
   the unordered pair.

## Licence

MIT.
