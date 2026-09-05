# CONTRACT: fs-cmaes-viz-wasm

Status: **CMA schema 2 + G1 schema 9 + ARM schema 4 — source interfaces.** Stateful
packed browser boundary over the production CMA-family owner in `fs-dfo` and
the owner-composed G1 walking and KUKA household-manipulation experiments. This
crate owns validation, admission, numeric packet transport, browser-specific
resource limits, and the disclosed benchmark tasks. It owns no optimizer
recurrence or independent robot mathematics.

## Purpose and layer

Layer L6 browser transport over `fs-dfo` and the named owner-composed G1 and
household experiments; it owns packet admission and rendering data, not the
underlying optimizer or robot mathematics.

## Public types and semantics

`PackedCmaSession`, `PackedG1WalkingEvaluator`, and
`PackedManipulationEvaluator` exchange the versioned packets documented below;
the wasm32 exports are their browser boundary.

## Invariants

Accepted asks and tells contain complete populations, only one ask is pending,
and snapshots expose only state actually retained by the selected owner.

## Error model

Malformed or inadmissible packets return the documented seven-word typed
refusal without trapping; a failed tell retains its pending population for a
repair-and-retry.

## Determinism class

Deterministic within a target for the same valid packet sequence, inherited
from the seeded `fs-dfo` sampling and tie-ordering contract.

## Cancellation behavior

None: this is synchronous browser transport and adds no cancellation or
parallel execution protocol.

## Unsafe boundary

No unsafe boundary is claimed; the manifest forbids unsafe code.

## Feature flags

`g1-learned` enables the optional `fs-g1-train` transformer adapter.

## Conformance tests

The focused native, wasm32, browser-packet, replay, and UBS gates listed under
Required gates exercise the declared packet and experiment contracts.

## No-claim boundaries

This adapter adds no restart strategy, generic objective, hardware-transfer,
mesh-contact, or dense limited-memory diagnostic claim; the detailed limits
remain stated in the existing ownership text below.

## Purpose and ownership

A browser creates one optimizer session, asks for a complete candidate
population, evaluates its own black-box objective, and tells one objective per
candidate. Every accepted state transition delegates directly to
`fs_dfo::CmaOptimizer`.

The available families are:

| ID | Family | Representation | Sampling / update |
|---:|---|---|---|
| 0 | Full | dense active covariance | O(n²) / O(n³) |
| 1 | Separable | active diagonal covariance | O(n) / O(n) |
| 2 | LM-CMA | bounded direction memory | O(mn) / O(m²n) worst case |
| 3 | LM-MA | bounded direction memory | O(mn) / O(mn) |

Full CMA is capped at 256 dimensions at this browser boundary. The scalable
families are capped at 100,000 dimensions and are admitted from checked,
owner-reported storage arithmetic. The adapter rejects a conservative browser
live envelope above 16 Mi binary64 words (128 MiB): owner persistent + pending
+ update workspace, two packed candidate-population transport copies, packet
headers, and a conservative snapshot payload. The 5,040-dimensional flagship
is comfortably inside this bound.

The G1 surface composes `fs-mbd`'s source-bound mode-11 Unitree model with all
29 actuated joints and its 5,040-D policy map, `fs-ga` poses, `fs-time`
integration, `fs-contact` normal response, and `fs-tribo` friction. It returns
the owners' link poses for rendering. The binding layer adds no kinematics,
dynamics, contact law, or browser-side pose reconstruction.

The manipulation surface likewise composes the existing mathematical owners:
`fs-mbd` supplies its pinned KUKA LBR iiwa 7 R800 topology, source dimensions,
inertias, hard limits, forward kinematics, inverse dynamics, and articulated-
body forward dynamics; `fs-ga` supplies SE(3) poses and wrenches; `fs-contact`
supplies the compliant pad normal response; `fs-tribo` supplies dry-friction
capacity; and `fs-query` supplies certified convex separation over admitted
oriented-box envelopes. Its 128 coordinates are 16 uniform knots for each of seven
joint targets and one finger-separation target. The binding returns source-
ordered world poses and never reconstructs robot kinematics in JavaScript.

The household grasp is deliberately reduced and explicit. A horizontal surface
supports the object until both finite pads engage. The rollout integrates
owner-reported normal force, Coulomb stick/slip traction, the free object's
translation, its body-frame angular velocity, and the reciprocal flange wrench;
there is no Boolean latch, teleport, or rigid pose following. After release, the
object follows ballistic translation, rotation, and a one-sided support response.
Kitchen mug, living-room remote, and backyard
trowel presets change disclosed dimensions, mass, grasp width, stations, and
keep-out box while sharing the exact same physics and policy contract.

Because the source catalog intentionally omits collision meshes, the experiment
declares four equal-height compliant patches under each foot. Patch forces and
moments accumulate on the source foot links. The reference posture starts at
the analytic Hertz indentation that shares the model's static weight across
the eight patches; it does not begin above unloaded contact and rely on a
timestep-dependent impact to discover support.

The experiment's semi-implicit step projects joint velocities onto the
catalog's symmetric source speed limits. Any proposed overshoot is accumulated
in `joint_limit_integral` before projection, so the optimizer sees a graded
penalty while the next owner call remains inside its admitted state domain.
Every uncompleted fixed step adds 1,000 objective units, and a terminal guard
adds one further 1,000-unit charge. The remaining physical shaping terms are
smoothly compressed to `400 * tanh(raw / 200)`, hence lie strictly inside
±400. Survival is therefore lexicographically primary: one additional completed
step dominates every possible shaping-score difference, so early termination
cannot win by skipping future integrated costs.

The admitted task selector is explicit: balance tracks zero velocity and double
support, stepping tracks alternating support and swing clearance at zero
commanded speed, and walking adds the configured forward-speed target. Each
task applies declared weights to dimensionless speed error, stance slip,
posture, contact-schedule mismatch, swing-clearance error, lateral/heading
error, joint-limit proximity, normalized excess load, cost of transport, and
flight time. Walking additionally charges backward travel. A non-finite raw
score is a typed refusal; the bounding transform never conceals overflow.

The standing PD controller grants the residual policy 65% of each source
actuator's effort envelope. The sparse standing mean is analytically rescaled
from the earlier 32% envelope so its initial residual torques are unchanged.
This leaves enough symmetric motor authority for a learned policy to unload a
foot and flex a swing knee instead of making the standing controller an
unreported hard constraint.

Two initialization vectors are public data. `stabilizing_policy_mean` has 15
nonzero constant-bias coordinates. `walking_curriculum_mean` has exactly 105
nonzero owner-layout coordinates: those biases, 30 first-harmonic phase
weights, and 60 gravity/angular-rate feedback weights. Its full-CMA seed,
population, generation, and sigma provenance is recorded beside the constants
in `g1_walking.rs`. The live scalable-CMA stage may mutate all 5,040
coordinates; neither vector is a hidden trajectory or browser controller.

No-claims: this synchronous adapter does not add restarts, parallel execution,
cancellation, constraints, or built-in analytic landscapes. The G1 experiment
learns 15 lower-body-and-waist channels while its 14 arm joints use a disclosed
deterministic reflex; fixed head and hand shells are not source collision
meshes. It is not a full Unitree digital twin or a hardware-transfer claim. The
arm experiment uses conservative oriented-box
envelopes rather than triangle meshes and has no general impact impulse solver,
general grasp planner, deformable object, or hardware-transfer claim.
Limited-memory snapshots do not invent dense
covariance diagnostics.

## Public surface

Native Rust exposes `PackedCmaSession::{new,receipt_packet,ask_packet,tell_packet}`
and `PackedG1WalkingEvaluator::{new,receipt_packet,evaluate_packet,
evaluate_population_packet,trace_packet}`,
`PackedManipulationEvaluator::{new,receipt_packet,curriculum_policy_mean,
evaluate_packet,evaluate_population_packet,trace_packet}`, the public
curriculum-mean constructors, and `kernel_version()`.

On wasm32, wasm-bindgen exports:

| Export | Signature | Result |
|---|---|---|
| `new CmaesVizSession(config)` | packed `Float64Array` | stateful session; construction never throws for bad input |
| `session.receipt()` | no arguments | admission/current snapshot or typed refusal |
| `session.ask()` | no arguments | complete row-major candidate population or typed refusal |
| `session.tell(objectives)` | packed `Float64Array` | updated snapshot or typed refusal |
| `new G1WalkingVizEvaluator(config)` | packed `Float64Array` | reusable owner-composed walking evaluator |
| `evaluator.receipt()` | no arguments | fixed controls and exact trace layout |
| `evaluator.stabilizing_policy_mean()` | no arguments | sparse 5,040-D standing mean |
| `evaluator.walking_curriculum_mean()` | no arguments | sparse 5,040-D walking curriculum mean |
| `evaluator.evaluate(policy)` | 5,040 policy words | decomposed scalar objective receipt |
| `evaluator.evaluate_population(policies)` | up to 64 row-major policies | one objective per candidate in one boundary call |
| `evaluator.trace(policy)` | 5,040 policy words | receipt plus decimated world-from-link poses |
| `new HouseholdManipulationVizEvaluator(config)` | packed `Float64Array` | reusable owner-composed KUKA evaluator |
| `evaluator.receipt()` | no arguments | fixed controls, scene, success criteria, and exact trace layout |
| `evaluator.curriculum_policy_mean()` | no arguments | source-feasible 128-D pick/lift/transport/place mean |
| `evaluator.evaluate(policy)` | 128 policy words | decomposed scalar objective receipt |
| `evaluator.evaluate_population(policies)` | up to 64 row-major policies | one objective per candidate in one boundary call |
| `evaluator.trace(policy)` | 128 policy words | receipt plus object and source-ordered link poses |
| `cmaes_viz_kernel_version()` | no arguments | `"fs-cmaes-viz-wasm 0.6.22"` |
| `cmaes_viz_source_revision()` | no arguments | build's `FSCMAES_SOURCE_REVISION`, or `"unbound"` for development |

There is no JSON hot path and no schema-1 compatibility shim.

## Numeric packet rules

Every word is one IEEE-754 binary64 value. Count, selector, and identifier
fields must be finite exact integers. Dimension, population, memory, budget,
and packet-length fields deliberately share wasm32's portable unsigned 32-bit
domain on native and browser builds. A one-word generation is limited to
JavaScript's exact nonnegative integer domain; 64-bit seeds and Philox block
counts use low/high unsigned 32-bit words.

Input packets use:

`[magic=0x434d4132, schema=2, kind, total_words, ...]`

Output packets use:

`[magic=0x434d4132, schema=2, status, kind, total_words, ...]`

Output status is 0 for success and 1 for refusal. Unknown magic, schema, kind,
selector, non-integral field, or inconsistent packet length fails closed.

### Configuration input (kind 0)

`magic, schema, kind, total_words, family, dimension, population_or_zero,
memory_or_zero, max_evaluations, seed_low32, seed_high32, sigma, mean[n]`

Zero selects the reference default population
`4 + floor(3 ln(n))`. Zero selects the reference dimension-based memory
default for LM-CMA and LM-MA; explicit memory is invalid for Full or Separable.
The budget admits only complete populations, and at least one population must
fit. The receipt reports the exact admitted budget.

### Admission/snapshot output (kinds 1 and 4)

Words 0–30 are:

`magic, schema, status, kind, total_words, family, dimension, generation,
evaluations, sigma, population, parents, max_generations,
admitted_evaluations, stream_semantics, stream_kernel, normal_blocks_low32,
normal_blocks_high32, sampling_order, update_order, persistent_scalars,
pending_scalars, update_workspace_scalars, dense_matrix_entries,
memory_capacity, has_best, best_objective, best_generation, best_candidate,
shape_kind, shape_payload_words`

The header is followed by `mean[n]`, `best_point[n]`, and the shape payload.
Before a best point exists, `has_best=0` and the best fields/point are NaN.
Complexity order IDs are 0 linear, 1 O(mn), 2 quadratic, 3 cubic, and 4
O(m²n). LM-CMA reports ID 1 for sampling and ID 4 for its corrected temporal-
memory update; LM-MA reports ID 1 for both.

Shape payloads are representation-honest:

- kind 0 Full:
  `negative_weight_count, min_eigenvalue, max_eigenvalue, covariance_diagonal[n]`
- kind 1 Separable:
  `negative_weight_count, variances[n]`
- kind 2 LM-CMA / LM-MA:
  `stored_vectors, memory_capacity, direction_norm[stored_vectors]`

### Ask output (kind 2)

`magic, schema, status, kind, total_words, generation, evaluations_before,
dimension, population, candidates[population * dimension]`

Candidates are row-major and retain the opaque ordering owned by `fs-dfo`.
Only one ask may be outstanding.

### Tell input (kind 3)

`magic, schema, kind, total_words, generation, population,
objectives[population]`

Every objective must be finite. A malformed count, non-finite objective, wrong
generation, or owner refusal retains the pending population so the caller can
repair the packet and retry. For generation/count/objective failures, refusal
word 6 names the outstanding owner generation, never an untrusted caller value.
A successful tell consumes exactly one complete population and returns a
kind-4 snapshot.

### Refusal output

Every refusal is exactly seven words:

`magic, schema, status=1, attempted_kind, total_words=7, refusal_code,
expected_generation_or_nan`

Stable refusal codes:

| Code | Meaning |
|---:|---|
| 1 | malformed packet |
| 2 | schema mismatch |
| 3 | unknown family |
| 4 | invalid/scalable-limit dimension |
| 5 | full-CMA browser dimension limit |
| 6 | invalid population |
| 7 | invalid/inapplicable memory |
| 8 | invalid or too-small budget |
| 9 | invalid seed words |
| 10 | invalid sigma |
| 11 | non-finite mean |
| 12 | checked owner shape overflow |
| 13 | Philox counter overflow |
| 14 | dense eigensolver admission refusal |
| 15 | browser memory-envelope refusal |
| 16 | ask already pending |
| 17 | exact budget exhausted |
| 18 | tell without a pending ask |
| 19 | generation mismatch |
| 20 | objective-count mismatch |
| 21 | non-finite objective |
| 22 | opaque owner batch mismatch |
| 23 | owner numerical failure |

Nothing is silently clamped and validation failures do not trap across the
wasm boundary.

## G1 walking packets

G1 packets use magic `0x47315737` (`"G1W7"`) and schema 9. The common output
prefix is `magic, schema, status, kind, total_words`. Kinds are configuration 0,
admission 1, evaluation 2, trace 3, and population 4.

Configuration is variable length: eleven fixed words, a keep-out box
count, then eight words per box.

`magic, schema, kind=0, total_words=12+8n, step_s, duration_s,
target_forward_speed_m_per_s, gait_frequency_hz, trace_stride, task,
challenge, obstacle_count=n`, then `n` groups of
`center_xyz_m, half_extents_xyz_m, yaw_rad, body_role` (0 keep-out, 1 support).

`total_words` is self-describing and must equal the packet length. Boxes
are capped at 64, must be finite, and must have strictly positive half
extents; anything else is refused, never clamped. An empty roster leaves
the rollout identical to schema 7.

The body-vs-obstacle guard has always been implemented: every step, each
body collider sphere is tested against every declared box, the deepest
penetration is tracked, and the first penetration past the 0.01 m skin
depth terminates the rollout as `BodyObstacle` with a shaped terminal
penalty. Until schema 8 no packet could declare a box, so the guard was
unreachable from the browser and the renderer carried the whole burden of
keeping the robot out of the furniture. It no longer does: the walking
policy is scored against solid geometry, so passing through a wall costs
the optimizer its rollout. The evaluation receipt now reports the deepest
penetration the guard measured, so the browser states the number the
kernel computed instead of re-deriving contact from rendered poses.

Task IDs are 0 balance, 1 stepping, and 2 walking.
Challenge IDs are 0 flat and 1 terrain plus lateral push. The combined challenge
uses the disclosed smooth height field
`0.008 sin²(2.4 x) (1 + 0.2 sin(3 y))` metres and a 24 N peak half-sine
lateral pulse from 0.55 through 0.70 seconds. Contact indentation, normal, and
tangential slip are all evaluated in the local terrain frame; the push is an
external root wrench applied 0.42 m above the root origin.

The browser experiment admits step sizes from 1/480 through 1/30 s, durations
through 4 s, commanded forward speeds through 2 m/s, gait frequencies from
0.25 through 4 Hz, and at most 1,000 steps between trace samples. These are
transport workload limits, not claims about the robot's validated operating
envelope.

Admission appends:

`policy_dimension=5040, link_count=30, pose_words=7, trace_sample_words=213,
step_s, duration_s, target_speed, gait_frequency, trace_stride, task,
challenge, terrain_amplitude_m, terrain_wavenumber_rad_per_m, push_start_s,
push_end_s, push_peak_force_n, physical_actuators=29, learned_rows=15,
reflex_actuators=14, features_per_row=336, phase_basis_count=8,
bias_count=15, phase_count=30, feedback_count=60,
arm_swing_gate_start_s=1/3.1, arm_swing_gate_end_s=3/3.1,
curriculum_indices[105]`.

Admission has 136 words including its five-word header. The indices come
from the actual initializer: first `336*r` for rows r=0..14, then
`336*r+[1,2]`, then `336*r+[248,256,272,280]`. The learned rows follow source
joint indices 0..14; the reflex joints are source indices 15..28. Pose order
is the pelvis followed by the 29 source-actuated links. The shoulder/elbow
swing multiplier is clamped cubic smoothstep between the two fixed physical
times. Gait frequency changes the phase signal, not those gate times.

Source identity is deliberately `unbound` unless the artifact build supplies
`FSCMAES_SOURCE_REVISION`. A release consumer must compare the revision and
packet schemas with its reviewed manifest and hash the actual JavaScript and
WASM bytes before execution; a version string alone is not build provenance.

Evaluation appends:

`objective, distance_m, speed_error_integral, actuator_work_j, slip_integral,
posture_integral, joint_limit_integral, impact_integral, backward_distance_m,
lateral_error_integral, heading_error_integral,
contact_schedule_mismatch_integral, swing_clearance_error_integral,
single_support_s, double_support_s, flight_s, push_impulse_n_s,
recovery_time_s, minimum_base_height_m,
maximum_tilt_sine, maximum_abs_terrain_height_m, completed_steps,
termination_reason, maximum_body_penetration_m`.

`recovery_time_s` is elapsed time after the push until the disclosed tilt,
angular-speed, and height bands are all regained. If recovery never occurs, it
is right-censored at the available post-push horizon rather than fabricated.

Termination IDs are 0 horizon, 1 base height, 2 base tilt, 3 contact
indentation, 4 contact speed, 5 contact-model domain, and 6 joint-position
limit. A candidate-dependent contact-domain failure is a terminal evaluated
outcome with ID 5 and a finite penalty, not a transport refusal that discards
the rest of a valid CMA generation.

Trace adds `sample_count`, followed by samples containing `time_s,
left_contact, right_contact` and 30 poses in catalog link order. Every pose is
world-frame `translation_xyz, quaternion_wxyz`; the browser is a projector and
must not recompute forward kinematics.

Population success adds `population, objectives[population]`. Input is a flat
row-major sequence of 1 through 64 complete 5,040-word policies. A structural
candidate refusal still fails the whole population packet and word 6 names its
zero-based row, so a caller cannot silently replace an owner refusal with a
fabricated score. All other G1 refusals use the same seven-word envelope with a
NaN detail word.

## Household-arm packets

Arm packets use magic `0x41524d31` (`"ARM1"`) and schema 4. The common output
prefix is `magic, schema, status, kind, total_words`. Kinds are configuration 0,
admission 1, evaluation 2, trace 3, and population 4.

Configuration is variable length: twelve fixed words followed by eight words
per caller-declared keep-out box.

`magic, schema, kind=0, total_words=12+8n, step_s, duration_s, trace_stride,
task, object_mass_kg, static_mu, kinetic_mu, obstacle_count=n`, then `n`
groups of `center_xyz_m, half_extents_xyz_m, yaw_rad, body_role` (0 keep-out,
1 support).

`total_words` is self-describing and must equal the packet length; a
disagreement is a malformed packet, not a clamp. Task IDs are 0 kitchen mug,
1 living-room remote, and 2 backyard trowel. The browser experiment admits
step sizes from 1/240 through 1/45 s, durations from 3 through 6 s, and at
most 1,000 steps between retained trace samples.

`object_mass_kg`, `static_mu`, and `kinetic_mu` are overrides: **zero selects
the owner's preset value**, so a schema-3 packet with all three zero and
`n=0` reproduces the schema-2 rollout bit for bit. That parity is pinned by
`schema_three_defaults_reproduce_the_preset_owner_rollout`, which compares the
packed objective against the same rollout driven through the internal API.
Non-zero overrides admit mass in `[0.02, 20]` kg and each coefficient in
`(0.05, 2.5]`, with kinetic never exceeding static; the interface provenance
string gains a `--caller-friction` suffix whenever either coefficient is
caller-declared, so an overridden interface can never be read back as the
owner's declared dry-pad default.

Declared keep-out boxes are yawed about the world +Z axis, capped at 32, with
centres within 10 m of the origin and half extents in `(0.001, 5]` m. They are
**hard link-vs-box constraints for all seven moving link segments**, exactly
like the preset box. The manipulated object is deliberately not scored against
them: it begins and ends resting on caller-declared support geometry, so an
object-vs-box test would refuse every rollout. Any envelope violation is a
refusal (`InvalidConfig`), never a silent clamp.

Admission appends:

`policy_dimension=128, joint_count=7, policy_knots=16, link_count=8,
pose_words=7, trace_sample_words=67, step_s, duration_s, trace_stride, task,
minimum_gripper_width_m, open_gripper_width_m, placement_tolerance_m,
lift_target_m, object_mass_kg, object_dimensions_xyz_m, grasp_half_width_m,
initial_object_xyz_m, goal_object_xyz_m, support_height_m,
obstacle_center_xyz_m, obstacle_half_extents_xyz_m, static_friction_mu,
kinetic_friction_mu, declared_obstacle_count`.

This self-describing 40-word packet is the renderer's source of truth.
`object_mass_kg` and the two coefficients report the values the rollout
**actually ran with**, preset or overridden, so the renderer never has to
guess which interface produced a receipt. The object dimensions and mass are
declared benchmark estimates; robot link dimensions, poses, inertias, and
limits remain source-bound in `fs-mbd`.

Evaluation appends:

`objective, final_object_error_m, minimum_reach_error_m, maximum_lift_m,
actuator_work_j, collision_risk_integral, minimum_certified_clearance_m,
possible_collision_time_s, collision_query_iterations, control_limit_integral,
first_grasp_time_s, grasp_duration_s, peak_grip_force_n, ever_grasped,
released_after_transport, placed, completed_steps`.

Each arm step constructs conservative oriented boxes for seven moving link
segments and the object. Containing spheres provide a certified broad reject;
near pairs delegate to `fs-query` convex separation for obstacle, proximal
object, and non-adjacent self checks. `possible_collision_time_s` counts steps
where separation was not proven, while the minimum clearance and iteration
count expose the query work instead of hiding it in a generic penalty.

`placed=1` requires an established bilateral grasp, release after transport,
at least the admitted lift target, terminal position error no greater than the
admitted tolerance, no active grasp at the horizon, zero integrated collision
risk, zero possible-collision time, and a certified clearance of at least
0.045 m throughout the rollout. The disclosed curriculum mean passes that full
contract for the kitchen-mug and living-room-remote presets. The backyard-trowel
mean completes grasp, lift, transport, and release but is truthfully refused as
a placement because its current trajectory violates the collision contract; it
remains a deterministic failing curriculum case rather than a hidden browser
animation.

Trace adds `sample_count`, followed by 67-word samples containing `time_s,
gripper_width_m, grip_normal_force_n, grasped`, one object pose, and eight link
poses in source catalog order. Every pose is world-frame `translation_xyz,
quaternion_wxyz`; the browser must not recompute forward kinematics.

Population success adds `population, objectives[population]`. Input is a flat
row-major sequence of 1 through 64 complete 128-word policies. A structural
candidate refusal fails the packet and word 6 identifies the zero-based row.
All other arm refusals use the same seven-word envelope with a NaN detail word.

## Determinism and budgets

Seeded candidate generation, normal-draw semantics, and deterministic
candidate-index tie ordering come from `fs-dfo` / `fs-rand`. The adapter
adds no entropy or ordering. Replaying the same sequence of valid packets
produces word-identical ask and snapshot packets within a target.

LM-CMA follows Loshchilov's corrected July 2014 `LMCMAfixed` source, not the
known-corrupt original recurrence. Each stored rank-one projection uses the
original isotropic sample `z`, and deleting a temporal-memory predecessor
recomputes every affected `v = A^-1 p` suffix direction. The corrected source
and its note about the original stale-`v` defect are published at
<https://sites.google.com/site/lmcmaeses/>. A direct `A^-1(Az) = z` regression
and a 40-generation, 5,040-dimensional plateau regression guard both details.

`max_evaluations` is a hard input budget. Admission truncates it to
`floor(max_evaluations / population) * population`; ask refuses once that
many evaluations have been consumed. Failed tells spend nothing and preserve
the pending batch.

## Required gates

- focused native tests, including four-family dispatch, seeded replay,
  repair-and-retry ask/tell semantics, exact budget exhaustion, typed
  refusals, a real 5,040-dimensional generation for every scalable family,
  two collision-certified household placements plus the deterministic trowel
  collision refusal, and a real 128-dimensional arm generation through every
  CMA family;
- strict crate-local Clippy with warnings denied;
- locked wasm32 check;
- release `wasm-pack build --target web` to a unique external cache directory;
- browser-side packet/replay exercise against the built artifact;
- changed-file UBS scan and bundle-size report.
