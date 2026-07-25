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

**Two tiers of hardware, and only one of them can draw.** The training node is
headless 2xH100. Hopper has no RT cores, no display engine and no video encoder,
and the box carries no EGL, Vulkan or OpenGL userspace at all, so Isaac Sim,
Unreal and Unity camera output cannot run there at any price. That is fine for
the part of the problem that lives there, because depth and range sensing are
raycasting problems rather than rasterisation problems, and a hand-written CUDA
sphere tracer renders 1920x1080 in 0.025 ms on one H100. Thousands of parallel
environments share one kernel launch.

Rendering is not gone from the project, it is relocated. The workstation tier is
an RTX 4090 and an RTX 5090, both of which have RT cores, display engines and
NVENC. Photorealistic imagery, which the planner does not want but the
perception front end genuinely does, is produced there. So the split is

- **H100 node**, headless: massively parallel physics, geometric sensing,
  planner benchmarks, and all network training. Deterministic and engine-free.
- **4090 / 5090 workstations**: photorealistic image generation for the
  perception stack, hardware ray tracing, PhysX cross-checks, and video encoding
  for the write-up.

CUDA therefore ships as a fatbin over sm_89 (Ada), sm_90 (Hopper) and sm_120
(consumer Blackwell), with a CPU fallback the full suite must pass under.
Note that sm_120 needs CUDA 12.8 or newer, and that PyTorch on a 5090 must be a
cu128 or later build. Earlier wheels contain no sm_120 kernels and fail at
runtime with "no kernel image is available for execution on the device".

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
| empty | 61.35 ± 1.12 | 57.52 ± 0.64 | **57.19 ± 0.83** | +3.83 ± 0.73 (t = 10.3) |
| forest, sparse | 62.73 ± 1.26 | 59.61 ± 0.93 | **58.93 ± 1.15** | +3.11 ± 0.62 (t = 9.8) |
| forest, dense | 70.23 ± 2.39 | **66.38 ± 1.44** | 66.39 ± 2.04 | +3.86 ± 1.31 (t = 5.8) |
| wall with window | 38.82 ± 0.83 | 35.19 ± 0.50 | **31.97 ± 0.32** | +3.63 ± 0.62 (t = 11.5) |

Rewiring beats plain RRT everywhere, significantly. Informed sampling is a large
win on the window map, where the informed ellipsoid collapses tightly around a
forced detour, and buys **nothing at all on the dense forest**, where the two
means differ by 0.01 m while informed sampling spends 26% more planning time
(62.4 ms against 49.5 ms). That negative result is reported rather than buried:
when the first solution is already close to the straight-line lower bound, the
informed set barely shrinks and the extra sampling machinery only adds overhead.

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

Three bugs this suite caught, all of which produce code that looks correct:

1. Seeding the distance transform with infinity makes the parabola intersection
   evaluate `INF - INF = NaN` on any scanline containing no obstacle.
2. The sphere trace originally strode by the full slack `d - clearance`, which
   is sound for the exact distance function because it is 1-Lipschitz. But the
   queried object is the *trilinear interpolant*, whose partial derivatives are
   each bounded by 1 and whose gradient norm therefore reaches √3 along a
   diagonal. The full stride steps over violations. Strides are now divided by
   √3, and `edgeFree` canonicalises endpoint order so it is a true function of
   the unordered pair.
3. The distance field was optimistic by half a voxel everywhere, because the
   transform returns centre-to-centre distance while the planner needs distance
   to the obstacle *surface*. At 0.5 m voxels that is 0.25 m of clearance the
   vehicle does not have. The exactness test in item 1's suite could not catch
   it, because the brute-force reference shared the same convention and was
   consistently wrong in the same direction. Catching it needed a test against
   a geometry with a closed-form answer, which is why there is now one against
   an analytic half space.

A note on what "reproducible" means here. Uniform draws are bit-portable across
platforms, since they use only integer operations and one multiply. Normal draws
are not, because Marsaglia's polar method needs a logarithm and neither glibc's
nor MSVC's `log` is correctly rounded. The same applies to the connection
radius. Reproducibility is bitwise within a platform and tolerance-based across
platforms, and the fixtures are written to match.

## Architecture note: two worlds

Ground truth and the agent's belief are separate, and the boundary is enforced
by the type system rather than by discipline. Most projects claiming
"GPS-denied" quietly hand the planner the true pose, and every number after that
is theatre.

- **Ground truth**, owned by the simulator. An analytic signed distance field
  over a CSG tree of primitives. One function is simultaneously the collision
  oracle, the depth-sensor model and the true clearance field, so the sensor and
  the collision checker cannot disagree, because they are the same function.
- **Belief**, owned by the agent. A rolling ego-centric voxel grid of log-odds
  occupancy plus the exact ESDF in `src/collision/`, built only from noisy depth
  projected through the *estimated* pose. It drifts and it ghosts, and the
  planner sees this and nothing else.

The planner and the mission layer accept only an `EstimatedState` carrying its
covariance. Ground truth is reachable through a debug-only accessor, and CI
fails the build if any results-generating script calls it.

## Licence

MIT.
