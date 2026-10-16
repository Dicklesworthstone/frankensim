# fs-regime — CONTRACT

The physics-regime and nondimensionalization kernel (plan patch Rev A):
the formal layer answering "WHICH SOLVER IS EVEN VALID for this physical
situation?" before FLUX runs. Regime checking turns silently-wrong solver
choices into structured, alternatives-ranked refusals.

Ambition tags: Pi machinery, named groups, admission gating [F per the
bead label; the implementation is exact algebra tested to [S] standards].

## Purpose and layer

Layer **L3** (FLUX support). Runtime deps: `std`, fs-qty (dimension
vectors), fs-evidence (`ModelCard`/`ValidityDomain`/`Evidence`), fs-math
(deterministic `ln`/`sqrt`). Consumers: fs-ir admission (HELM), fs-lbm's
lattice-scaling assistant, FLUX conformance, conformal hardening buckets.

## Public types and semantics

- `pi`: `pi_groups(&[Input]) -> PiBasis` — the integer nullspace of the
  6×n SI dimension matrix `[m, kg, s, K, A, mol]` by exact fraction-free
  elimination over i128.
  `PiGroup` exponents are reduced (gcd 1, leading positive); values are
  products of SI input values (hence unit-free). Exact exponents outside
  the current deterministic i32 power domain are refused before numerical
  evaluation; they are never truncated. `rank + groups.len() == n`
  (Buckingham's theorem, by construction for admitted bases).
- `groups`: `standard_groups(&[RoleInput])` — named groups from
  role-tagged inputs: Re, We, Ca, Oh, Bo, Fr, Ma, St, De, slenderness,
  damping ratio ζ, P-Delta index. Every formula runs through QtyAny
  dimension arithmetic and REFUSES if not dimensionless (wrongly-tagged
  inputs cannot silently produce a "group").
- `scaling`: `ScalingMap::recommend` (L* from Length, T* = L*/U*,
  M* = ρL*³, with explicit K*, A*, and mol* base scales), `apply`/`unapply`
  via per-dimension scale factors;
  `condition_number` — an exact 2-norm condition probe (cyclic Jacobi on
  AᵀA) for fixture-scale measurements.
- `cards`: `flux_model_cards()` — the built-in registry (Stokes/creeping,
  laminar NS, LES [F, labeled], free-surface LBM, potential flow,
  Euler-Bernoulli, Timoshenko) with validity boxes in group space;
  `admit(registry, groups, model) -> Admission { allowed, reasons,
  alternatives }` with alternatives ranked by log-decade distance to
  their validity boxes; `distance_to_validity` is the ranking metric and
  `axis_distance_to_validity` exposes its per-axis severity law.
- `output_audit`: `audit_product_output` checks every user-facing QoI at
  every named operating point against every consumed model card. Its typed
  receipt preserves exact in-domain/out-of-domain partitions, card/version,
  axis, observed value, bounds, log-distance severity, original/effective
  colors, per-partition colors, and any explicit override acknowledgement.
  One or more violations force the effective color to
  `Estimated { dispersion: infinity }`; an acknowledgement permits downstream
  policy to proceed but never restores color. Canonical JSON and deterministic
  no-claim Markdown are report/ledger handoff projections. Each exact canonical
  receipt also has a domain-separated content identity.
  `apply_output_audit_to_budget` leaves a matching fully in-domain eight-term
  budget byte-for-byte unchanged. A demoted receipt replaces its `ModelForm`
  term with explicit `Unknown` evidence naming every point/card/axis violation
  and its `fs-regime` distance, while retaining the prior model-form state and
  authority in the diagnosis. The replacement provenance is the exact receipt
  identity. If model form belonged to a covariance block, every member is made
  unknown because the finite joint representation cannot survive removal of
  one member.
  `audit_product_output_with_cards` accepts the narrower `RegimeAuditCard`
  projection (owner-defined name, version, and exact shared validity box) for
  matdb and other card classes whose authoritative schemas do not contain the
  extra ambition, discrepancy, calibration, and failure fields required by
  `fs_evidence::ModelCard`. The projection does not invent those fields;
  `audit_product_output` remains the exact evidence-card wrapper.
- `report`: `assess(&[RoleInput]) -> Evidence<RegimeReport>` — groups,
  Pi rank/count, dominant balance, valid/invalid models, recommended
  scaling, conditioning risk (decade spread of input scale factors), and
  the nearest canonical benchmark (cylinder Re=100, Stokes sphere,
  lid-driven cavity Re=1000, dam break) with expectation text and
  info/warning/far grading. `BenchmarkMatch::evidence_ref` links registered
  expectations to an executable repository battery without introducing a
  runtime dependency on the solver crate.

## Invariants

1. **Dimensionless by construction**: every Pi group's integer exponent
   combination is verified against the dimension matrix (exact integer
   arithmetic); every named group's formula is dimension-checked.
2. **Buckingham count**: group count = inputs − rank, always.
3. **Six-base atomicity**: amount of substance is an independent matrix
   row and scaling axis; mol-free inputs retain their prior rank, Pi basis,
   and values because the appended exponent is exactly zero.
4. **G3 unit-rescaling invariance**: group values and Pi bases are
   EXACTLY invariant under any coherent rescaling of the SI base units
   (tested with raw rescaled numbers, not just re-parsed spellings).
5. **Refusals are actionable**: an inadmissible model yields the violated
   bounds by name and value plus a full ranked alternative list.
6. **Reports are reproducible**: same inputs → identical report and
   identical provenance hash.
7. **Final-envelope demotion is monotone**: shrinking a card validity domain
   cannot move an operating point from the out-of-domain partition to the
   in-domain partition. Partial sweeps expose both exact partitions and are
   never averaged into a single coverage fraction.
8. **Overrides cannot launder color**: override acknowledgements are retained
   in the receipt but do not participate in effective-color computation.
9. **Out-of-domain model uncertainty is absorbing**: applying a demoted receipt
   cannot produce a finite `ModelForm` term. A QoI mismatch refuses, and prior
   model-form evidence is named rather than silently discarded.

## Error model

`RegimeError`: `NotDimensionless { context, residual }`,
`ExponentOutOfRange { context, exponent }`, `Degenerate`, `MissingRole`,
`UnknownModel`, `BadValue`. Group formation refuses
non-positive inputs (powers of signed values have no regime meaning).

## Determinism class

**D0**: exact integer elimination, deterministic `ln`/`sqrt` via
fs_math::det, BTreeMap ordering everywhere a report is rendered.

## Cancellation behavior

All bounded, small (n inputs ≤ dozens; Jacobi probe n ≲ 32). P7 by
boundedness.

## Unsafe boundary

Zero `unsafe`.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs` (JSON verdicts, suite `fs-regime/conformance`):

- **rg-001** textbook batteries (pipe flow, drag, heat convection, molar
  concentration): rank, Buckingham count, exact six-base integer
  dimensionlessness of every group.
- **rg-002** the definitive G3 test: base-unit rescaling (mm/g/min-ish
  factors plus a nontrivial mol scale applied through raw exponent
  arithmetic) leaves legacy mol-free named groups and Pi bases exactly
  invariant (< 1e−12 relative).
- **rg-003** seeded misuse: creeping-flow solver at Re = 10⁴ refused with
  the violated bound named; LES ranked as a distance-0 alternative; a
  valid LBM case admits; unknown models are structured errors.
- **rg-004** scaling improves conditioning ≥100× on a mixed-magnitude
  3-DOF fixture (before/after condition numbers ledgered as JSON).
- **rg-005** similarity: Re ≈ 98 matches the cylinder Re=100 benchmark
  with the two-width blockage-sensitivity Cd `[1.25, 1.45]` and 16D
  lift-FFT St `[0.155, 0.175]` plus 12D sensitivity targets at "info" grade,
  and points exactly to
  `crates/fs-lbm/tests/cylinder_re100.rs::lbm_109_cylinder_re100_cd_and_strouhal`;
  reports reproduce bit-identically with identical provenance.
- **rg-006** flagship fixtures vs hand calculations: spout Re/We/Ca/Oh/Bo
  (admits creeping, refuses LES), ornithoid Re/Ma/St, frame
  slenderness/ζ.

Unit tests cover the Pi machinery directly (Reynolds/pendulum/drag
recovery, mol-free basis invariance, mass isolation, nonzero-mol rank and
residual cancellation, degenerate refusals), all six scaling axes, and the
condition probe against known diagonals.

`tests/output_audit.rs` covers per-card checking, named distance-scored
violations, partial sweep partitioning, monotone domain shrinkage,
deterministic receipts, the non-restoring override law, and conservative
eight-term `ModelForm` propagation.

## No-claim boundaries

- **The registry is a v0 seed**, not a complete FLUX catalog: solvers
  register their own cards as they land; bounds here encode textbook
  regime limits (e.g. laminar pipe transition at Re ≈ 2300), not
  solver-implementation-specific calibration.
- **Dominant-balance text is a heuristic reading** of the group values
  (teaching output), not a certified asymptotic analysis.
- **The condition probe is a measurement tool** for fixture-scale
  matrices, not a linear-algebra kernel — fs-la owns real factorizations.
- **Benchmark table is small and curated**; similarity is log-distance
  over shared groups, not a learned embedding. Additions are data. An evidence
  reference identifies the intended executable gate; it does not itself assert
  that a proof-pending or ignored release lane has passed.
- **PDE-level nondimensionalization** (rewriting operators) lives with
  the solvers; this crate recommends and applies SCALES to quantities.
- **This core audit does not prove product-wide wiring.** fs-report must place
  its Markdown projection in the final no-claim section, fs-package must retain
  the receipt, and CLI/project orchestration must supply the complete consumed
  card set and final operating envelope before `f85xj.8.3` can close.
