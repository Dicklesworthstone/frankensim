# fs-solid CONTRACT

## Purpose and layer

Layer: L3 (FLUX). The elasticity core (plan §8.2, bead tfz.13):
small-strain elasticity through finite-strain hyperelasticity, with a
MEASURED locking-free formulation, on both the body-fitted and the
CutFEM-on-SDF frontends. Constitutive laws live in fs-material; this
crate owns kinematics, weak forms, boundary operators, and the
globalized Newton loop.

## Public types and semantics

- `Mesh2` / `Patch`: structured body-fitted meshes — P1 triangles, Q1
  quads, mapped Q1 panels (`cooks_membrane`), with the four sides as
  named patches. `shapes_at`/`quad_points` expose the isoparametric
  machinery (chain rule with the transposed inverse Jacobian — the
  conformance battery pinned this exact slip on triangles and mapped
  panels).
- `LinearProblem` (plane strain/stress via `PlaneKind`, `lame`):
  small-strain elasticity with strong Dirichlet and dead-load traction
  patches; `Formulation::Standard` vs `Formulation::BBar` (Hughes'
  dilatation projection — element-average volumetric gradient,
  equivalent to condensed constant-pressure mixed elements).
- `HyperProblem` + `NewtonSettings`/`NewtonReport`: plane-strain
  finite strain; the 2D displacement gradient embeds into the 3D
  deformation gradient (F₃₃ = 1); residual `∫ P : ∇δu` and consistent
  tangent `∫ ∇δu : A : ∇Δu` come from fs-material's exact AD Piola
  stress and 9×9 tangent. Newton globalizes by Armijo backtracking on
  the potential energy (material refusals — det F ≤ 0 — halve the
  step) under linear load stepping; systems solve with MINRES
  (symmetric, indefinite-safe near buckling). `residual_and_tangent` /
  `potential_energy` are the public consistency probes.
- `CutElasticity`/`CutSolution`: vector Q1 linear elasticity on
  fs-cutfem background quadtrees — certified classification and cut
  quadrature reused verbatim; VECTOR symmetric Nitsche (full traction
  consistency `σ(N_a e_i)·n`, penalty β(λ+2μ)/h tied to certified cut
  cells); componentwise ghost penalty γ(λ+2μ)h.
- `select_formulation(RegimeIndicators)`: the documented
  element-selection policy — B-bar when ν ≥ 0.45 or slenderness ≥ 5
  (fs-regime's indicators feed this; the locking battery measures the
  boundary).
- `accept_scenario_bc`: fs-scenario integration — `Dirichlet` and
  `Traction` under `Physics::Elasticity` accepted; dimensioned-value
  resolution stays with the scenario consumer.

## Invariants

1. Patch tests exact: linear displacement fields are reproduced to
   solver tolerance on P1, Q1, B-bar, and the CutFEM frontend
   (all-embedded fixture) — sol-001.
2. MMS optimal orders per implemented family × frontend: L2 ≈ 2,
   H1 ≈ 1 for P1/Q1/CutFEM-Q1 (sol-002).
3. Objectivity through the fs-material interface: W(RF) = W(F) and
   P(RF) = R·P(F) for both cards (sol-003).
4. Consistency merge gate: the assembled tangent equals the FD
   directional derivative of the residual, and the residual is the
   exact gradient of the potential energy (rel ≤ 1e-6) — sol-004.
5. Locking is MEASURED, not assumed: on thin near-incompressible
   bending, the standard element's error degrades by ≥10× from
   ν = 0.3 to ν = 0.4999 while B-bar stays within 2× (sol-005).
6. Cook's membrane envelope: converged tip at the literature reference
   point (48, 52) inside [23.5, 24.5] (plane stress, ν = 1/3), ≤2%
   self-convergence deviation; near-incompressible B-bar
   self-converges while the standard element's shortfall is logged
   (sol-006).
7. Newton robustness: large-strain bending and a buckling-adjacent
   compressed strip converge under load stepping + line search with
   fast terminal contraction; histories and backtracks are evidence
   (sol-007).
8. Determinism: BTree/insertion-ordered assembly, deterministic
   Krylov solvers; bit-identical across runs on a platform.

## Error model

`SolidError` teaching errors: `SolveFailed` (linear gate missed),
`NewtonStalled` (carries the residual history; repair = more load
steps), `MaterialRefused` (det F ≤ 0 escaped globalization; repair =
smaller steps), `UnknownPatch`, `UnsupportedBc` (fs-scenario mapping
outside Dirichlet/Traction × Elasticity).

## Determinism class

Bit-deterministic across runs on a fixed platform (ordered assembly,
fs-solver CG/MINRES, no threading, no ambient state). Cross-ISA
golden hashes are not yet recorded (follow-up).

## Cancellation behavior

Bounded synchronous loops (assembly, capped Krylov iterations, capped
Newton with capped backtracking). Chunked Cx polling belongs to the
fs-exec driver (the L3 discipline; fs-feec/fs-cutfem precedent).

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/battery.rs`: sol-001 patch tests (body-fitted max nodal error
≤ 1e-9; CutFEM at the CG-tolerance floor ≤ 1e-7); sol-002 MMS orders
(P1 2.00/1.00, Q1 2.00/1.00, CutFEM 2.19/1.04); sol-003 objectivity
(energy rel ≤ 1e-10, stress rotation ≤ 1e-8, 8 random states × 2
cards); sol-004 tangent/energy consistency (≤ 1e-6); sol-005 locking
battery (standard degrades ≥10×, B-bar ≤2×); sol-006 Cook's envelope;
sol-007 Newton robustness (5-step large bending, 8-step compression,
terminal contraction ≥ 100×); the selection-policy probe. Unit tests:
mapped-mesh areas, selection thresholds.

## No-claim boundaries

- TDNNS-proper / weakly-imposed-symmetry FEEC stress elements: awaits
  the simplicial H(curl)/H(div) family bead (dcng); B-bar is the
  shipped, measured locking-free path with the same acceptance metric.
- 3D and shells; higher-order families (Qk/Pk, k ≥ 2); IGA frontend
  (bead tfz.9).
- Hyperelasticity on the CutFEM frontend (linear ships there; the
  finite-strain cut path composes this crate's Newton loop with
  fs-cutfem quadrature in a successor bead).
- fs-opdsl-generated residual/adjoint paths (the DSL coverage bead);
  the hand path here passes the same consistency gates the DSL output
  must pass.
- Contact, plasticity flow (fs-material's plastic module wires in a
  successor), dynamics/time integration.
- Ogden energies (fs-material's recorded no-claim propagates).
- fs-scenario dimensioned-value resolution (units plumbing stays with
  the scenario consumer; kind/physics validation ships here).
