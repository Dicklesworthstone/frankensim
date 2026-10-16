# CONTRACT: fs-truss-e2e

TrussPath — a deterministic truss iterate with an advisory, endpoint-checked
critical load path and a strictly gated positive interval certificate. Layer
L4 (ASCENT).

## Purpose and layer

Composes `fs-truss` (ground-structure LP + PDHG diagnostics), `fs-tropical`
(critical path), and `fs-evidence` (honest evidence state). Deps point downward.

## Public types and semantics

- `run_campaign(nx, ny, w, h, gap_tol, cx) -> Result<TrussReport, TrussError>` —
  optimizes a cantilever ground structure by PDHG and extracts a checked
  tropical path. The final iterates also pass through `fs-truss`'s outward
  optimum-certificate kernel before the report can carry `Color::Verified`.
  Grid dimensions must be at least 2x2. Admission bounds cubic ground-generation
  work, candidate members, sparse PDHG scalar work, certificate work and retained
  state, active tasks, and path edges. Ground rules, the support/load case, and
  the assembled sparse LP cross the immutable fallible `fs-truss` admission
  boundary before solver work.
- `analyze_load_path(...)` — shared native/WASM path analysis. It admits only
  unique in-range identities, finite positive weights, and a connected chain
  of at least two strictly support-ward bars from the indexed load node to an
  indexed support.
- `certify_load_path(...)` consumes the private `fs-truss` optimum receipt and
  its outward repaired signed/split-force boxes. It intervalizes the relative
  active threshold, separates every included/excluded member, outwardly
  multiplies each split-force sum by its admitted cost, orients every active
  member only when support-distance intervals are disjoint, removes members
  outside load-reachable/support-coreachable paths, and promotes only when one
  complete path lower bound exceeds every rival upper bound and one positive
  bottleneck lower bound exceeds every peer upper bound.
- `LoadPathCertificate` is private-by-construction. It retains the three
  collision-resistant `fs-truss` problem/input/witness identities plus the
  exact node/member/load/support/threshold identity, interval weights, active
  set, path, and bottleneck. `verifies_for` replays the bounded proof and
  requires exact receipt equality. `replay_golden` is only a 64-bit native/WASM
  drift sentinel; it is not certificate authority.
- `optimality_color_from_certificate(problem, x, y, settings, status, gap,
  eq_residual)` — the sole native/browser promotion gate. Only a structurally
  valid private certificate bound to those canonical arrays, iterates, and
  settings yields `Verified { lo, hi }`. The gate checks retained shapes and the
  certificate's operation cap before hashing caller arrays; every mismatch,
  work excess, or unavailable proof remains `Estimated`.
- `rescale_optimality_color(color, positive_divisor)` — preserves an existing
  Verified interval through outward division (used for normalized-to-physical
  yield-stress scaling), preserves weaker colors, and demotes invalid scaling.

## Invariants

- OPTIMALITY: `solver_converged` records only the declared gap and equilibrium-
  residual diagnostics. Independently, `fs-truss` proves an exactly feasible
  primal repair with a Neumann enclosure and a feasible scaled dual with outward
  slack checks. Only finite ordered endpoints bound to the exact LP, settings,
  iterates, methods, limits, and retained witness become `Verified`.
- LOAD PATH: active bars are oriented by strict distance-to-support progress.
  Reachability and co-reachability filter out disconnected components and
  interior-only chains. The max-plus witness must start at the load, end at a
  support, contain at least two bars, and be joint-continuous. A bottleneck is
  named only for a unique path with one strictly heaviest positive bar.
- A `Verified` load-path color can come only from a retained
  `LoadPathCertificate`; its interval encloses the selected material-volume
  sum. Any threshold, orientation, path, bottleneck, arithmetic, or identity
  ambiguity keeps the rounded advisory path `Estimated` with infinite
  dispersion. The color certifies the declared graph/LP arithmetic, not a
  continuum stress-flow interpretation.
- The optimizer prunes the candidate set (`num_active < num_members`).
- Deterministic (fixed ground structure + deterministic PDHG; no RNG).

## Error model

`TrussError` refuses degenerate/oversized grids, unsafe geometry scales,
non-finite or non-positive tolerances, excessive construction/solver/path work,
allocation failure, cancellation, malformed certificate state, empty member
sets, malformed path data, and incomplete load-support chains. Certificate work
or conditioning limits are a sound numerical unavailability represented by an
`Estimated` color, not an error and never a partial `Verified` result.
`LoadPathCertificateRefusal` separately names unavailable optimum evidence,
solver-identity mismatch, unseparated active membership or orientation,
incomplete paths, path/bottleneck ties, and non-finite interval arithmetic.
Caller input does not panic, and a path refusal is never converted to positive
evidence.

## Determinism class

Fully deterministic (G5).

## Cancellation behavior

Ground construction and LP assembly poll the caller's `Cx` at deterministic
bounded strides and return a structured cancellation refusal without publishing
partial state. The certificate proof also polls the same `Cx` through admission,
repair, verification, identity binding, and atomic publication. The fixed PDHG
solve remains synchronous and iteration-bounded; solver-loop cancellation is a
separate successor. Load-path certification polls admission, interval
construction, support-distance evaluation, reachability, and receipt replay;
cancelled attempts publish no `LoadPathCertificate`.

## Unsafe boundary

None; `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/truss.rs`: independent diagnostic convergence; finite positive outward
optimality and load-path bounds; exact path-receipt replay and altered-endpoint
identity refusal; interval member-product containment; deterministic native
goldens; iteration-cap truthfulness; invalid/exact-bound work admission;
index-based supports; disconnected-heavy-component filtering; and retained
near-threshold, equal-distance, direct-one-bar, and tied-bottleneck falsifiers.
A pre-cancelled context is refused before a campaign report exists.

## No-claim boundaries

The optimum interval bounds the plastic ground-structure LP optimum; it does
not certify catalog sizing, buckling, or geometric nonlinearity. The separately
verified "load path" is only a material-volume longest chain through a
certificate-separated active set, not a stress-flow, continuum load-transfer,
redundancy, fatigue, or failure-mode proof. A direct one-bar route remains
outside the declared multi-bar witness and therefore Estimated. Sizing to a
catalog and buckling checks (`fs-truss::size_and_snap`,
`rod_buckling_check`) are downstream and not exercised here. The 64-bit replay
golden detects drift but provides neither collision resistance nor authenticity;
those properties come from the retained exact identity and lower-layer BLAKE3
receipt.
