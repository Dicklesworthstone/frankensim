# fs-dwr CONTRACT

## Purpose and layer

Layer: L3 (FLUX). Dual-weighted-residual goal-oriented adaptivity
(plan §8.6 [F], bead tfz.23): the adjoint solution weights local
residuals by their influence on the quantity being optimized, and ONE
signal drives all four refinement mechanisms — octree h-refinement,
anisotropic metric synthesis, the h-vs-p decision, and wavelet tile
thresholding. In a design-optimization system, accurate OBJECTIVE and
GRADIENT is the goal, not accurate simulation; DWR is the estimator
that knows the difference.

## Public types and semantics

- `GoalContext` / `goal_value`: volumetric goal functionals
  `J(u) = ∫ jw·u` (region averages, windowed integrals) evaluated
  with fs-cutfem's certified cut quadrature.
- `estimate` → `DwrEstimate`: primal solve on the given quadtree,
  ENRICHED adjoint solved on `Quadtree::refined_once` (the documented
  higher-resolution enrichment; patch recovery is the recorded
  alternative), and SIGNED per-cell indicators from the full discrete
  residual — interior `∫ f·w − ∇u_h·∇w` plus the Nitsche interface
  terms in fs-cutfem's exact sign convention. `eta_signed`
  approximates `J(u) − J(u_h)`; `eta_abs` is the marking mass. The
  coarse ghost-penalty residual contribution is a deliberately
  omitted O(γh)-scaled correction absorbed into the measured
  effectivity band.
- `dorfler`: fixed-energy marking, DETERMINISTIC — |indicator|
  descending, cell key ascending on ties, smallest prefix reaching
  θ·total. Bitwise-reproducible (P2).
- `adapt_loop` → `AdaptStep` rows: solve → estimate → mark → split →
  rebalance → RESTORE the uniform cut band (fs-cutfem's ghost-penalty
  precondition) at the finest cut-adjacent level; ledger-style JSON
  per step (dofs, J, η, marked).
- `synthesize_metric`: per-cell 2×2 metrics from recovered Hessians
  (two-sided second differences — a one-sided stencil vanishes by
  antisymmetry at an odd layer's inflection, a measured pathology)
  weighted by adjoint importance, anisotropy-capped 100:1, floored,
  and complexity-normalized so Σ√det(M)·|K| meets the target. The
  3D-embedded form is fs-mesh `MetricField`-compatible; unstructured
  execution through fs-mesh's remesher is the consumer wiring.
- `haar_threshold` → `ThresholdOutcome`: 2D Haar with per-coefficient
  budgets taken as the MINIMUM of the local budget over the covered
  block (conservative); DWR-weighted budgets spend accuracy where the
  adjoint says the goal cannot see.
- `h_vs_p` → `Decision`: the smoothness classifier
  `s_K = h·|H|_F/(|∇u|+δ)` routing kinks/layers to h and smooth
  regions to p. Emitting decisions only — executing local p awaits
  the high-order FEEC families.

## Invariants

1. Effectivity: signed estimate over true goal error within [0.5, 1.6]
   on known-truth MMS goals across BOTH frontends (all-embedded disk;
   strong+Nitsche strip) at two levels (dwr-001).
2. G3 monotonicity: `eta_abs` decreases under uniform refinement at
   rate ≥ 1.2 (theory 2; measured ~1.7–2.0).
3. Marking is bitwise-deterministic and the marked set is a MINIMAL
   Dörfler prefix (dwr-002).
4. Goal-oriented beats uniform on localized QoIs: strictly better
   accuracy (≤ 0.5×error) at no more DOFs, accuracy-per-DOF curves
   ledgered (dwr-003).
5. Metric synthesis: implied complexity within 5% of target; layer
   alignment ≈ 1.0; a metric-instantiated graded mesh halves the
   isotropic interpolation error at equal DOF (dwr-004).
6. Weighted thresholding: ≥5× compression with goal impact under the
   budget, and ≥2× better goal impact than unweighted thresholding at
   MATCHED compression (dwr-005).
7. h-vs-p: >90% correct routing on a kink+smooth composite (dwr-006).
8. Determinism: BTree traversal, deterministic solves and marking —
   bit-identical runs.

## Error model

fs-cutfem's `CutFemError` teaching errors propagate unchanged
(build/solve refusals). Marking on empty/zero indicators returns
empty. `synthesize_metric`/`h_vs_p` panic (structured asserts) on
non-uniform grids — the documented v1 surface, not a recoverable
state.

## Determinism class

Bit-deterministic across runs on a fixed platform (inherits
fs-cutfem's discipline; no ambient state, no threading). Cross-ISA
golden hashes not yet recorded (follow-up).

## Cancellation behavior

Bounded synchronous loops (solves are fs-cutfem's, estimator sweeps
are linear in cells, the adaptive loop has a fixed iteration count).
Chunked Cx polling belongs to the fs-exec driver (L3 discipline).

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None. The plan marks §8.6 [F]; per the crate-granular gating rule
(fs-cutfem/fs-feec precedent) the frontier surface ships as this
standalone crate.

## Conformance tests

`tests/battery.rs`: dwr-001 effectivity + monotonicity (two
frontends, known-truth goals); dwr-002 marking determinism + minimal
prefix; dwr-003 localized-QoI adaptive-vs-uniform accuracy-per-DOF;
dwr-004 metric synthesis (complexity/alignment/graded-beats-iso);
dwr-005 weighted Haar thresholding vs matched-compression unweighted;
dwr-006 h-vs-p routing. Unit tests: Haar lossless/mean roundtrips.

## No-claim boundaries

- p-enrichment EXECUTION (local high-order spaces await dcng/FEEC-p;
  this crate emits the routing decisions).
- Unstructured anisotropic remeshing execution (fs-mesh's remesh
  consumes the synthesized `MetricField`-compatible tensors; the
  graded tensor-product instantiation here is the shipped proof of
  the metric's value).
- Patch-recovery adjoint enrichment (higher-resolution solve ships;
  recovery is the cheaper documented alternative).
- The ghost-penalty residual term (O(γh) correction absorbed into
  the measured effectivity band).
- fs-opdsl-generated DWR residual terms (the DSL bead's one-source
  path; the hand path here passes the gates the generated one must).
- Error-Ledger/fs-plan budget reallocation wiring (AdaptStep rows are
  the ledger-ready shape; the composed-budget loop is Bet 12's bead).
- Graded-tree metric recovery and time-dependent (reverse-sweep)
  DWR.
