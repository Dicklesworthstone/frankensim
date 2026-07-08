# fs-topols CONTRACT

## Purpose and layer

Layer: L4 (ASCENT). Level-set topology optimization (plan §9.5 [S/F],
bead 7tv.12): shape-gradient velocity advection with topological
derivatives for hole nucleation — genuine, mathematically justified
topology changes. The level set IS the geometry: physics is evaluated
by fs-solid's CutFEM elasticity directly on the evolving discrete SDF
with zero meshing anywhere in the loop.

## Public types and semantics

- `GridSdf`: nodal φ on a uniform `[0,1]²` lattice, bilinear between
  nodes, implementing fs-cutfem's `CutSdf`. The enclosure is EXACT up
  to an outward roundoff pad — a bilinear attains its extrema at the
  corners of any axis-aligned rectangle, so a box's range is the hull
  of per-cell clipped-corner evaluations (certified classification
  holds on the moving geometry).
- `weno::advect` + `Velocity::{Linear, Normal}`: WENO5 (Jiang–Shu)
  spatial stencils with TVD-RK3 time stepping under a CFL number;
  linear advection for the order battery, Godunov/Rouy–Tourin normal
  flow `φ_t + v_n|∇φ| = 0` for the optimizer; band-masked updates.
- `fim::redistance` → `RedistanceAudit`: fast-iterative-method eikonal
  relaxation (Godunov upwind updates, alternating-order Gauss–Seidel
  to an order-independent fixed point) with the interface FROZEN at
  values reconstructed from zero crossings; audits carry the sampled
  zero-set Hausdorff drift (units of h) and |∇φ|−1 statistics — the
  redistancing-frequency policy's inputs. `zero_crossings`/`hausdorff`
  are the audit primitives.
- `veloext::extend_velocity`: interface-normal extension by
  ascending-|φ| upwind sweeps (deterministic ordering).
- `topder::topological_derivative`: the compliance sensitivity to
  infinitesimal traction-free hole insertion (plane-strain
  Amstutz-form coefficients), NUMERICALLY GATED — the battery punches
  a real hole and checks the measured compliance change against
  `DT·πρ²` within a documented first-order band, so wrong signs or
  scales cannot ship. `nucleate`: greedy best-gain hole punching with
  spacing, box-edge margins (clamps/loads excluded), and per-event
  ledger rows (`NucleationEvent`).
- `optimize::optimize_compliance` + `OptimizeSettings`/`Cantilever` →
  `OptimizeReport`: the descent loop — CutFEM solve (traction-free Γ,
  strong box-edge clamps, box-edge dead tractions), nodal
  strain-energy densities sampled inside material, normal extension,
  fs-adjoint Sobolev H¹ smoothing, `v_n = w − ℓ`, one interface move
  per iteration, redistance + audit, augmented-Lagrangian volume
  multiplier NORMALIZED BY THE MEAN ENERGY SCALE (an O(1) multiplier
  against O(J) energies shrinks the structure to nothing at full
  speed — measured failure mode), scheduled nucleation; ledger rows
  with compliance, volume, ℓ, drift, and FNV snapshot hashes.
- `optimize::material_volume`: certified cut-quadrature area of
  `{φ < 0}`.

## Invariants

1. Enclosure containment holds on the DISCRETE moving geometry
   (bilinear corner-hull law, dense-sample battery).
2. WENO order on smooth advection: measured order > 1.6 at the finest
   pair, absolute error < 5e-4 after half a revolution (tls-001;
   design order 3 with RK3 dominates asymptotically).
3. Redistancing moves the sampled zero set < 0.2h, restores
   |∇φ| = 1 to mean deviation < 0.05, and is idempotent (< 0.05h on
   repeat) — tls-002.
4. Extension is constant along interface normals (radial fixture
   deviation < 0.15 across a 0.13-wide annulus) — tls-003.
5. Volume is conserved under rigid rotation through five
   advect+redistance cycles within 1%, with per-cycle drift audits
   < 0.25h (tls-004 — the drift POLICY evidence).
6. Sensitivities are numerically gated: topological-derivative
   prediction vs a real punched hole and shape-velocity prediction vs
   a real uniform boundary motion both land in [0.25, 4]× with the
   right SIGN (tls-005) — the bead's adjoint-vs-FD gate.
7. The cantilever descent converges to the volume target within 0.05,
   stabilizes (tail variation < 20%), fires nucleation with positive
   predicted gain producing a GENUINE interior hole (flood-fill
   verified), keeps every redistancing drift < 0.5h, and BEATS the
   trivial uniform-band design at equal volume (tls-006).
8. Determinism (P2): two descent runs produce bitwise-identical FNV
   snapshot sequences.

## Error model

fs-solid's `SolidError` propagates from physics solves. Structured
asserts (panics) guard programmer contracts: lattice/grid alignment,
nodal-array lengths, non-uniform grids where v1 requires uniform.
Audits never silently degrade — drift and gradient deviations are
returned, not clamped.

## Determinism class

Bit-deterministic across runs on a fixed platform: BTree/index-ordered
sweeps, deterministic solvers underneath, snapshot hashes asserted
bitwise in the battery. Cross-ISA goldens not yet recorded.

## Cancellation behavior

Bounded synchronous loops (fixed iteration counts, FIM sweep caps,
CFL-bounded step counts). Chunked Cx polling belongs to the fs-exec
driver (L4 consumes L3 kernels; the fs-cutfem/fs-solid discipline).

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None. §9.5's [S] core (evolution, redistancing, coupling) and the [F]
nucleation operator ship together as a standalone L4 crate per the
crate-granular gating rule.

## Conformance tests

`tests/battery.rs`: tls-001 WENO order; tls-002 FIM audits +
idempotency; tls-003 normal extension; tls-004 volume conservation +
drift policy; tls-005 numerical sensitivity gates (DT vs punched
hole, shape velocity vs FD); tls-006 the cantilever descent (volume,
stabilization, nucleation with interior-hole flood-fill proof,
determinism, beats-trivial-at-equal-volume). Unit tests: enclosure
containment, planar extension.

## No-claim boundaries

- SIMP cross-validation on a SHARED fixture: fs-topopt's SIMP is
  3D-tet-based while this crate is 2D plane strain — a
  dimension-matched pairing (2D density pipeline or 3D level set)
  closes it; recorded follow-up. The shipped quality gate is
  beats-trivial-at-equal-volume plus the sensitivity gates.
- Octree narrow bands (uniform-lattice bands ship; fs-rep-sdf tile
  storage wiring is the consumer path per the plan's shared-substrate
  design).
- 3D, multiple load cases, stress constraints, compliant mechanisms.
- fs-adjoint Hadamard boundary-form velocities (the volumetric
  energy-density form ships; the Hadamard trace form composes when
  the FEEC trace machinery lands).
- Perimeter/curvature regularization beyond Sobolev smoothing;
  velocity CFL coupling to BAND MAINTENANCE beyond the shipped
  one-move-per-iteration policy.
- The vessel flagship's lip-channel fixture (§15.3) — the cantilever
  classic ships; the flagship composes downstream.
