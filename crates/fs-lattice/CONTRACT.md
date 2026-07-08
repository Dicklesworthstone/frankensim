# fs-lattice CONTRACT

## Purpose and layer

Layer: L4 (ASCENT). Lattice/infill optimization (plan §9.5 [F], bead
7tv.14), 2D SMOKE TIER and honest about it: periodic unit-cell
homogenization with physics-bound audits, graded macro-optimization
through the fitted homogenized law, and the separation-of-scales
validity flag doing real work — closed end-to-end by the
de-homogenization re-analysis (lat-007). The 3D TPMS families (gyroid
stiffness curves vs literature) are recorded successors.

## Public types and semantics

- `homogenize::UnitCell`: per-element density on an n×n grid;
  constructors for the centered-hole plate (density knob = hole
  radius) and the cross-strut cell (knob = arm half-width).
- `homogenize::Homogenizer::effective` → `EffectiveTensor` (Voigt
  3×3 + density): u = E·x + u_per split, cell stiffness from the
  fs-solid hyper2d tangent at u = 0 (exact linearization), PERIODIC
  master–slave reduction on the structured mesh's exact opposite-edge
  correspondence (corners fold to one pinned node), effective tensor
  by energy averages over the three unit-strain cell problems.
  Microstructure enters as a per-element stiffness contrast
  (void_eps; the CutFEM-exact variant is a recorded successor).
  The reference-element scatter is indexed by the element mesh's
  GLOBAL node ids (bl, br, tl, tr) — a CCW ordering was MEASURED to
  swap the top nodes and silently relax C11 from 3.5 to 1.5 while
  leaving C22 exact; lat-001's continuum-moduli gate is the tripwire.
- `graded::PropertyFit`: the property manifold — homogenization
  samples s(ρ) with piecewise-linear eval/slope, clamped to the
  sampled VALIDITY DOMAIN, plus the declared `gradation_bound` (the
  separation-of-scales model card).
- `dehomog::HoleArray` (a `CutSdf`: material = hole complement, with
  conservative per-hole interval enclosures) + `fullres_compliance`
  (CutElasticity re-analysis under the macro load convention —
  uniform right-edge traction, trapezoidal edge quadrature, so
  compliances compare like-for-like).
- `graded::graded_compliance_opt` → `GradedDesign`: cantilever
  compliance minimization over per-element cell densities through the
  fitted law; self-adjoint sensitivities, OC-style update with tight
  move limits (0.9/1.1 — wider limits measurably oscillated), volume
  bisection, and the final state re-analyzed AFTER the last update
  (the pre-update report was measured stale by one iteration). The
  scale-separation audit reports the worst adjacent density jump and
  raises `scale_separation_violated` past the fit's declared bound.

## Invariants

1. Homogeneous cells are EXACT: C_hom equals both the single-element
   probe (8.9e-15) and the continuum isotropic moduli C11 = λ+2μ,
   C12 = λ, C33 = μ, at every resolution (lat-001).
2. The hole-density sweep is monotone, SPD with square symmetry, and
   under the Voigt mixture bound at every sample; the measured dilute
   slope d(C11/C11s)/df = 2.93 at f = 0.028 is consistent with the
   classic 3f dilute-hole constant (the formal literature comparison
   row remains PENDING) (lat-002).
3. Resolution movement n = 8 → 16 stays within 16% — including the
   real density difference the element-center hole classification
   resolves (0.750 vs 0.703), printed with the gate (lat-003).
4. The cross-strut cell carries the cubic-anisotropy fingerprint:
   shear-compliant (C33/C11 = 0.082 vs solid 0.286) and 45°-soft
   (C45/C11 = 0.671 vs 1.0) (lat-004).
5. Graded beats the equal-mass uniform baseline by 26.8% on the
   smoke cantilever; volume budget respected; the aggressive fixture
   fires the separation-of-scales flag (worst jump 0.668 > 0.35) —
   REPORTED, never silent (lat-005).
6. The property fit passes through s(1) = 1 with positive slopes
   inside its validity domain (lat-006).
7. DE-HOMOGENIZATION closes the loop: the graded field realized as an
   explicit hole array (r = h·√((1−ρ)/π) per macro cell) and
   re-analyzed at full resolution with fs-cutfem CutElasticity
   (traction-free hole boundaries, matched load convention) lands
   within 4.7% of the homogenized prediction (band 15%), and the
   realized graded array still beats the realized uniform array by
   13.5% — the gradation survives realization (lat-007).
8. Empty realized hole arrays represent the solid material domain with
   finite negative SDF values/enclosures and radius range (0, 0), so
   solid designs do not leak infinities into CutFEM classification or
   diagnostics (lat-008).

## Error model

Structured panics on programmer contracts (mesh/cell mismatch,
singular reduced systems) with teaching messages; optimization
quality is reported through the design record's ledger fields, never
asserted silently.

## Determinism class

Bit-deterministic per platform: fixed assembly and iteration orders,
dense LU, no RNG anywhere.

## Cancellation behavior

Bounded synchronous loops (fixed iteration counts, bisection caps);
chunked Cx polling belongs to the fs-exec driver layer.

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None (the smoke tier ships enabled; heavier tiers will gate).

## Conformance tests

`tests/battery.rs`: lat-001 solid-cell exactness (continuum moduli +
element-probe identity); lat-002 sweep bounds + dilute ledger;
lat-003 resolution; lat-004 strut anisotropy fingerprint; lat-005
graded-vs-uniform + scale-separation flag; lat-006 property-fit
sanity; lat-007 de-homogenization prediction band + design transfer;
lat-008 empty-hole-array finite solid SDF behavior.

## No-claim boundaries

- 3D TPMS families (gyroid/Schwarz) and their literature stiffness
  curves — need 3D elasticity; recorded successor with the formal
  Hashin–Shtrikman 3D bound audit.

- CutFEM-exact cell geometry (the contrast-density approach ships;
  exact cut cells are the successor), orientation-graded anisotropic
  cells, stress-constrained objectives.
- fs-fab manufacturability audits (min feature size, powder escape)
  and FrankenNetworkx conforming lattice generation.
- fs-material registration of homogenized laws with validity-domain
  model cards (the PropertyFit carries the domain; the registry
  integration is staged with fs-material's card store).
