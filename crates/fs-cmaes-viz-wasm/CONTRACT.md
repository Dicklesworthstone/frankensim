# CONTRACT: fs-cmaes-viz-wasm

Status: **schema 2 — live surface.** Stateful packed browser boundary over the
production CMA-family owner in `fs-dfo`. This crate owns validation, admission,
numeric packet transport, and browser-specific resource limits. It owns no
optimizer recurrence.

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
| 2 | LM-CMA | bounded direction memory | O(mn) / O(mn) |
| 3 | LM-MA | bounded direction memory | O(mn) / O(mn) |

Full CMA is capped at 256 dimensions at this browser boundary. The scalable
families are capped at 100,000 dimensions and are admitted from checked,
owner-reported storage arithmetic. The adapter rejects a conservative browser
live envelope above 16 Mi binary64 words (128 MiB): owner persistent + pending
+ update workspace, two packed candidate-population transport copies, packet
headers, and a conservative snapshot payload. The 5,040-dimensional flagship
is comfortably inside this bound.

No-claims: this synchronous adapter does not add restarts, parallel objective
evaluation, cancellation, constraints, or built-in landscapes. It does not
invent dense covariance diagnostics for a limited-memory representation.

## Public surface

Native Rust exposes `PackedCmaSession::{new,receipt_packet,ask_packet,tell_packet}`
and `kernel_version()`.

On wasm32, wasm-bindgen exports:

| Export | Signature | Result |
|---|---|---|
| `new CmaesVizSession(config)` | packed `Float64Array` | stateful session; construction never throws for bad input |
| `session.receipt()` | no arguments | admission/current snapshot or typed refusal |
| `session.ask()` | no arguments | complete row-major candidate population or typed refusal |
| `session.tell(objectives)` | packed `Float64Array` | updated snapshot or typed refusal |
| `cmaes_viz_kernel_version()` | no arguments | `"fs-cmaes-viz-wasm 0.5.0"` |

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
Complexity order IDs are 0 linear, 1 O(mn), 2 quadratic, and 3 cubic.

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

## Determinism and budgets

Seeded candidate generation, normal-draw semantics, and deterministic
candidate-index tie ordering come from `fs-dfo` / `fs-rand`. The adapter
adds no entropy or ordering. Replaying the same sequence of valid packets
produces word-identical ask and snapshot packets within a target.

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
