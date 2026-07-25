# Roadmap and decisions

Living document. Records what was decided and why, so the reasoning survives
past the commit that implemented it.

## Architecture: the Two-World Core

Ground truth and the agent's belief are separate worlds, and the boundary is
enforced by the type system rather than by discipline.

The motivation is blunt. A GPS-denied project is only credible if the plan is
computed from what the agent could actually have sensed. Most projects claiming
"GPS-denied" quietly hand the planner the true pose, and every number after that
is theatre.

**World 1, ground truth**, owned by the simulator. An exact analytic signed
distance field over a CSG tree of primitives (box, cylinder, capsule, cone,
plane), with a two-level BVH carrying triangle leaves for imported meshes. One
function, `scene_sdf(p)`, is simultaneously the collision oracle, the depth
sensor model (sphere-traced) and the true clearance field. The sensor and the
collision checker cannot disagree, because they are the same function.

Analytic primitives also settle a real memory argument against voxel ground
truth. A 480 x 480 x 80 fp16 grid is 37 MB per scene, so 1024 randomised scenes
would be 38 GB of 80 GB HBM, against roughly 100 kB per scene for primitives.
Voxels also alias: a 0.16 m trunk vanishes at 0.25 m resolution.

**World 2, belief**, owned by the agent. A rolling ego-centric voxel grid at
0.10 m holding fixed-point log-odds occupancy, plus the exact ESDF in
`src/collision/`, rebuilt each frame by the separable transform. It is built
only from noisy depth projected through the **estimated** pose, so it drifts and
ghosts and has to be replanned against.

**The boundary.** `GroundTruthView` is constructible only inside `sensors/` and
`sim/`. The planner, controller, mission layer and policy accept `EstimatedState`
(carrying its covariance) and `const BeliefMap&`. Ground truth is reachable only
through a debug-only oracle accessor, and CI fails the build if any
results-generating script calls it.

## Decisions

| decision | choice | why |
| --- | --- | --- |
| Visualisation | Browser WebGL and 2D canvas | The dev box has no graphics userspace at all, and the algorithms are what is worth showing. Doubles as the project page. |
| Estimation | Staged: drift-injected odometry, then a real error-state KF | Keeps the stack end-to-end from week one without permanently faking the defining hard part. The synthetic-error stage is labelled a sensitivity study, never a headline. |
| Planning space | R^3 geometric, dynamics via differential flatness | The radius `(log n/n)^(1/d)` has shrunk only 2.7x at d=9 after 1e5 samples. Flatness recovery is exact, not an approximation. |
| Ground-truth world | Analytic CSG SDF | One function serves as sensor, collision oracle and truth. Voxels cost 370x the memory and alias thin geometry. |
| Belief map | Voxel log-odds + exact EDT | The agent's map should be discretised and wrong in realistic ways. That is the point. |
| Generative assets | **GAN cut.** Diffusion over SDF grids for benchmark scene generation | With no rasteriser, generated appearance has no consumer, and for line-of-sight occlusion a cylinder is a cylinder. Mode collapse in a GAN is silent and catastrophic for a benchmark generator, where diffusion fails visibly. Adversarial patrolling units remain, as mission logic rather than a generative problem. |
| GPU budget | Learned RRT* sampler and the scene generator | Redirected from the cut GAN work. |
| CUDA targets | Fatbin over sm_89, sm_90, sm_120, plus a tested CPU fallback | Real target hardware is an RTX 4090 (Ada), the H100 training node (Hopper) and an RTX 5090 (consumer Blackwell). sm_120 needs CUDA 12.8+. The suite must pass under `-DAAL_DISABLE_CUDA=ON`. |
| Photorealistic rendering | On the 4090 / 5090 tier only, never in the training loop | Those cards have RT cores, display engines and NVENC. The H100 has none of the three. Pixels matter to the perception front end and not to the planner, so they are produced where they are cheap. |

## Milestones

M0 through M7 is the complete story. M8 and M9 are packaging.

| id | name | deliverable | acceptance |
| --- | --- | --- | --- |
| M0 | Vertical slice | One drone, three primitives, end to end on CPU | Reaches a goal 12 m away through an obstacle, seeded and reproducible |
| M1 | Determinism | Pinned RNG, FP-contraction policy, banned-symbol lint | `ctest -L portability` passes on GCC, Clang and MSVC |
| M2 | World and sensing | Analytic CSG SDFs, BVH, procedural scene families, depth model | Noise-free depth matches analytic ray-primitive distance to 1e-5 over 1e6 rays |
| M3 | Dynamics and control | RK4 quadrotor, mixer with thrust-priority desaturation, SE(3) controller | RK4 order fit gives p in [3.8, 4.2]; forward Euler gives p in [0.9, 1.1] |
| M4 | Estimator | 15-DoF error-state KF, FEJ toggle, chi-square gating, TUM-VI validation | Average NEES over 50 runs inside [13.52, 16.56] with FEJ on, outside with it off, as a CI assertion |
| M5 | Belief map | Fixed-point log-odds fused through the estimated pose, CUDA EDT | Matches brute force exactly on the integer field; matches analytic half-space, sphere and plane |
| M6 | Planner | Anytime RRT*, Informed, min-snap smoothing, external baselines | Every returned path collision-free against the **analytic** SDF, not merely the belief map |
| M7 | **The headline** | Estimated-versus-oracle gap over 200 held-out maps, matched seeds | The table exists, every cell has an interval, and the write-up states what the sample size cannot resolve |
| M8 | Missions | Behaviour tree, four missions, fault injection | Return-to-home lands with energy above zero on 200 of 200 seeds |
| M9 | Page and reproducibility | Visualisations, zero-console-error Playwright gate | Stable across themes and breakpoints |

## Open, needs hardware

An Orin NX 16GB (~600 USD) would close the largest structural gap in the plan.
Every performance number here is on 2xH100, and a GPS-denied drone flies a 15 W
companion computer. "Replan p99 is 47 ms on an Orin NX at 15 W" is worth more
than any throughput number on a datacenter part.

No real-flight validation in version one. TUM-VI gives real inertial data and
real feature tracks against published baselines, which is the credibility that
actually transfers. This is stated on the page rather than left for a reviewer
to notice.

## Hardware tiers

| tier | hardware | RT cores | NVENC | role |
| --- | --- | --- | --- | --- |
| Training | 2x H100 80GB (sm_90) | no | no | Parallel physics, geometric sensing, planner benchmarks, all training. Headless, no graphics userspace installed. |
| Workstation | RTX 4090 24GB (sm_89) | yes | yes | Photorealistic perception imagery, Isaac Sim or Unreal, PhysX cross-check, video encoding, interactive debugging. |
| Workstation | RTX 5090 32GB (sm_120) | yes | yes | Same, plus the newest-architecture test target. Requires CUDA 12.8+, driver 570+, and a PyTorch cu128 build. |

What this changes relative to the original plan: the perception front end was
going to synthesise observations at the disparity and feature-track level
because no machine in the project could produce pixels. That constraint is gone.
Pixel-level rendering now happens on the workstation tier, which makes learned
depth, learned feature detection and sim-to-real domain randomisation reachable
rather than aspirational.

What it does not change: the training loop stays engine-free. An engine in the
inner loop would cost determinism, throughput and portability, and would tie
every published number to a specific Omniverse or Unreal build. The engine is a
data source and a validation cross-check, not the simulator.
