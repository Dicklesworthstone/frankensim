# CONTRACT: fs-cmaes-viz-wasm

Status: **CMA schema 2 + G1 schema 2 — live surfaces.** Stateful packed browser
boundary over the production CMA-family owner in `fs-dfo` and the owner-composed
G1 walking experiment. This crate owns validation, admission, numeric packet
transport, and browser-specific resource limits. It owns no optimizer recurrence
or independent robot mathematics.

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

The G1 surface composes `fs-mbd`'s source-bound 15-DoF Unitree model and 5,040-D
policy map, `fs-ga` poses, `fs-time` integration, `fs-contact` normal response,
and `fs-tribo` friction. It returns the owners' link poses for rendering. The
binding layer adds no kinematics, dynamics, contact law, or browser-side pose
reconstruction.

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
smoothly compressed to `400 * tanh(raw / 10,000)`, hence lie strictly inside
±400. Survival is therefore lexicographically primary: one additional completed
step dominates every possible shaping-score difference, so early termination
cannot win by skipping future integrated costs.

The raw secondary score is
`-18 distance + 12 speed_error_integral + 0.008 actuator_work_j +
16 slip_integral + 30 posture_integral + 2 joint_limit_integral +
0.8 impact_integral + terminal_guard_penalty`. A non-finite raw score is a typed
refusal; the bounding transform never conceals overflow.

No-claims: this synchronous adapter does not add restarts, parallel execution,
cancellation, constraints, or built-in analytic landscapes. The G1 experiment
is a reduced lower-body-and-waist model, not a full Unitree digital twin or a
hardware-transfer claim. Limited-memory snapshots do not invent dense
covariance diagnostics.

## Public surface

Native Rust exposes `PackedCmaSession::{new,receipt_packet,ask_packet,tell_packet}`
and `PackedG1WalkingEvaluator::{new,receipt_packet,evaluate_packet,
evaluate_population_packet,trace_packet}`, plus `kernel_version()`.

On wasm32, wasm-bindgen exports:

| Export | Signature | Result |
|---|---|---|
| `new CmaesVizSession(config)` | packed `Float64Array` | stateful session; construction never throws for bad input |
| `session.receipt()` | no arguments | admission/current snapshot or typed refusal |
| `session.ask()` | no arguments | complete row-major candidate population or typed refusal |
| `session.tell(objectives)` | packed `Float64Array` | updated snapshot or typed refusal |
| `new G1WalkingVizEvaluator(config)` | packed `Float64Array` | reusable owner-composed walking evaluator |
| `evaluator.receipt()` | no arguments | fixed controls and exact trace layout |
| `evaluator.evaluate(policy)` | 5,040 policy words | decomposed scalar objective receipt |
| `evaluator.evaluate_population(policies)` | up to 64 row-major policies | one objective per candidate in one boundary call |
| `evaluator.trace(policy)` | 5,040 policy words | receipt plus decimated world-from-link poses |
| `cmaes_viz_kernel_version()` | no arguments | `"fs-cmaes-viz-wasm 0.5.7"` |

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

G1 packets use magic `0x47315732` (`"G1W2"`) and schema 2. The common output
prefix is `magic, schema, status, kind, total_words`. Kinds are configuration 0,
admission 1, evaluation 2, trace 3, and population 4.

Configuration is:

`magic, schema, kind=0, total_words=9, step_s, duration_s,
target_forward_speed_m_per_s, gait_frequency_hz, trace_stride`.

The browser experiment admits step sizes from 1/480 through 1/30 s, durations
through 4 s, commanded forward speeds through 2 m/s, gait frequencies from
0.25 through 4 Hz, and at most 1,000 steps between trace samples. These are
transport workload limits, not claims about the robot's validated operating
envelope.

Admission appends:

`policy_dimension=5040, link_count=16, pose_words=7, trace_sample_words=115,
step_s, duration_s, target_speed, gait_frequency, trace_stride`.

Evaluation appends:

`objective, distance_m, speed_error_integral, actuator_work_j, slip_integral,
posture_integral, joint_limit_integral, impact_integral, completed_steps,
termination_reason`.

Termination IDs are 0 horizon, 1 base height, 2 base tilt, 3 contact
indentation, 4 contact speed, 5 contact-model domain, and 6 joint-position
limit. A candidate-dependent contact-domain failure is a terminal evaluated
outcome with ID 5 and a finite penalty, not a transport refusal that discards
the rest of a valid CMA generation.

Trace adds `sample_count`, followed by samples containing `time_s,
left_contact, right_contact` and 16 poses in catalog link order. Every pose is
world-frame `translation_xyz, quaternion_wxyz`; the browser is a projector and
must not recompute forward kinematics.

Population success adds `population, objectives[population]`. Input is a flat
row-major sequence of 1 through 64 complete 5,040-word policies. A structural
candidate refusal still fails the whole population packet and word 6 names its
zero-based row, so a caller cannot silently replace an owner refusal with a
fabricated score. All other G1 refusals use the same seven-word envelope with a
NaN detail word.

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
  refusals, and a real 5,040-dimensional generation for every scalable family;
- strict crate-local Clippy with warnings denied;
- locked wasm32 check;
- release `wasm-pack build --target web` to a unique external cache directory;
- browser-side packet/replay exercise against the built artifact;
- changed-file UBS scan and bundle-size report.
