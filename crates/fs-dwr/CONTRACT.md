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
  over the built `Space` using its retained active topology and certified
  cut quadrature. Evaluation returns `Result` and refuses missing or
  non-finite nodal/goal evidence rather than substituting zero.
- `estimate` → `DwrEstimate`: primal and adjoint solves on the given
  quadtree, an ENRICHED adjoint solved on `Quadtree::refined_once` (the
  documented higher-resolution enrichment; patch recovery is the recorded
  alternative), and SIGNED per-cell indicators from the full discrete
  residual — interior `∫ f·w − ∇u_h·∇w` plus the Nitsche interface
  terms in fs-cutfem's exact sign convention. `eta_signed`
  approximates `J(u) − J(u_h)`; `eta_abs` is the marking mass. The
  coarse ghost-penalty residual contribution is a deliberately
  omitted O(γh)-scaled correction absorbed into the measured
  effectivity band. The two-level weight is evaluated pointwise as
  `z_{h/2} - z_h`; every requested corner in either active space must
  exist and be finite, so a non-nested topology never becomes a plausible
  zero coefficient.
- `estimate_elasticity_compliance` → `ElasticityDwrEstimate`: vector
  CutFEM compliance estimator with `J_h = b_h^T u_h`. Symmetry makes the
  compliance adjoint equal to the primal, so the enriched weight is the
  pointwise, fail-closed difference `w = u_{h/2} - u_h`; there is no
  second adjoint solve and no missing-node zero fallback. Signed bulk,
  symmetric-Nitsche, outer-traction, and stabilization terms are exposed
  separately. Cell indicators reconstruct `eta_signed`; `eta_abs` is only
  Dörfler marking mass; canonical coarse-face indicators reconstruct the
  reported stabilization correction. `cell_terms` exposes the same signed
  bulk, Nitsche, outer-traction, and ghost decomposition per coarse cell while
  `indicators` remains their complete marking sum.
- `ElasticityGhostMethod::CoarseConsistentEnergy`: the elasticity ghost
  term is the coarse consistent-limit correction `+g_h(u_h,u_h)` on the
  actual coarse ghost-face set, split equally between its two cells. It is
  `-g_h(u_h,w)` in the smooth-adjoint limit, where the exact adjoint has no
  normal-derivative jump and `jump(w) = -jump(u_h)`. This avoids inventing
  enriched traces on inactive halves of a coarse ghost face. It is measured
  estimator evidence, not a certified bound.
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
9. Vector residual signs match fs-cutfem assembly: bulk
   `f·w - sigma(u_h):epsilon(w)`; Nitsche
   `t(u_h)·w + t(w)·(u_h-g) - (beta*mu/h)(u_h-g)·w`; outer
   traction `t_bar·w`; coarse consistent-limit stabilization
   `+g_h(u_h,u_h)`.
   The signed cell allocation reconstructs the term sum to roundoff.
10. Vector compliance effectivity uses an independent reference solve
    finer than the estimator's enriched solve. The signed ratio
    `eta_signed / (J_ref - J_h)` must be finite, non-degenerate, and in
    the fixture's measured acceptance band; `eta_abs` is never substituted
    for the signed numerator.
11. Scalar enrichment is fail-closed with bidirectional direct coverage: every
    enriched active cell has an active direct coarse parent, every coarse
    active cell has at least one enriched active child, and every field or
    coarse/enriched adjoint corner needed by a retained quadrature rule exists and is
    finite. One active child is sufficient; four-child coverage is not
    invented as a stronger topology precondition.
12. A non-final adaptive iteration with marked cells derives its reserved
    balance/interface refinement headroom with checked arithmetic. A zero-level
    grid refuses as `InvalidFemInput`; it never underflows or silently wraps the
    refinement ceiling. An empty marked set requires no headroom.

## Error model

fs-cutfem's `CutFemError` teaching errors propagate unchanged
(build/solve refusals). The scalar estimator and `goal_value` use
`InvalidFemInput` to refuse the level-16 enrichment cap, active-space
parent/child mismatch, missing or non-finite nodal/adjoint-corner data,
non-finite goal/residual terms or totals, retained topology/rule mismatch, or a
zero-level adaptive grid that needs another marked refinement iteration.
Non-finite primal source/boundary callbacks may first surface through the
propagated scalar build/solve teaching error because they enter assembly before
the residual sweep.
The vector estimator makes the analogous refusals through its documented
elasticity-input path and additionally refuses duplicate canonical ghost
faces or failed signed allocation reconstruction. Marking on empty/zero
indicators returns empty. `synthesize_metric`/`h_vs_p` panic (structured
asserts) on non-uniform grids — the documented v1 surface, not a recoverable
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

`tests/elasticity.rs`: independent-reference vector compliance
effectivity rows across bulk/traction, embedded-Nitsche/ghost, and
plate-with-hole fixtures; exact zero/empty stabilization evidence when
`ghost_gamma = 0`; signed term/cell/face reconstruction; deterministic
replay; and falsifiable top-decile localization of the non-traction residual
at the hole boundary. That residual is bulk plus ghost in the traction-free
hole fixture; Nitsche remains in the generic decomposition but is zero there.
The separately reconstructed dead-load traction term is intentionally excluded
from that spatial claim because it correctly localizes on the loaded design-box
edge; it remains in the actual marking indicator and effectivity estimate.
Scalar estimator unit regressions cover a missing active nodal corner, a
missing enriched-field corner, non-finite nodal/goal evidence,
fine-to-coarse and coarse-to-fine
active-space mismatch, valid one-child coverage, and level-16 refusal before
`refined_once` can panic. A level-4 disk regression proves that legitimate
absent mapped fine nodes succeed through the coarse-adjoint two-level weight;
a level-15 boundary regression retains exactly one legal enrichment level.

## No-claim boundaries

- p-enrichment EXECUTION (local high-order spaces await dcng/FEEC-p;
  this crate emits the routing decisions).
- Unstructured anisotropic remeshing execution (fs-mesh's remesh
  consumes the synthesized `MetricField`-compatible tensors; the
  graded tensor-product instantiation here is the shipped proof of
  the metric's value).
- Patch-recovery adjoint enrichment (higher-resolution solve ships;
  recovery is the cheaper documented alternative).
- The scalar estimator's coarse ghost-penalty residual term (O(γh)
  correction absorbed into its measured effectivity band). The vector
  estimator instead reports the explicitly named coarse consistent-energy
  correction; it is the smooth-adjoint limit of the coarse ghost residual,
  not an exact evaluation against the non-nested enriched field.
- fs-opdsl-generated DWR residual terms (the DSL bead's one-source
  path; the hand path here passes the gates the generated one must).
- Error-Ledger/fs-plan budget reallocation wiring (AdaptStep rows are
  the ledger-ready shape; the composed-budget loop is Bet 12's bead).
- Graded-tree metric recovery and time-dependent (reverse-sweep)
  DWR.
- A certified upper/lower error bound. Every vector `eta` is Estimated;
  effectivity is fixture evidence, not a certificate.
- Physical-work semantics for nonzero embedded displacement data. In that
  case `b^T u` is only the algebraic assembled-load compliance.
- A graded vector-elasticity re-solve or repeated vector adaptive loop.
  The current vector CutFEM frontend refuses graded active trees until
  componentwise hanging constraints land; integration may drive the
  deterministic split policy but must not claim the post-split solve.
- Exact treatment of the non-nested active-space variational crime. The
  coarse consistent-energy ghost method is logged on every estimate and its
  adequacy is measured by independent-reference effectivity fixtures.
