# CONTRACT: fs-plate

## Purpose and layer
Orthotropic thin-plate bending for instrument bodies: DKT (discrete
Kirchhoff triangle) elements with an orthotropic bending-rigidity matrix
from `fs_material::OrthotropicElastic`, membrane-prestress geometric
stiffness, offset stiffener beams, and modal analysis through
`fs_modal::slice_window` (certified intervals, inertia-certified counts).
Layer: **L3 FLUX** (deps: fs-material L3, fs-modal L1, fs-sparse L1).
Bead frankensim-fsim-plates-shells-kj3s0 (musical-acoustics program).

## Public types and semantics
- `PlateSection::{orthotropic, isotropic}` — bending rigidity D from
  engineering constants: `D11 = E1h³/12(1−ν12ν21)` etc., `D33 = G12h³/12`.
  Material axis 1 is the grain (L) axis — the matdb axis-convention
  contract. SI units (Pa, m, kg/m³); the dimensioned fs-qty front door is a
  recorded follow-up on the bead.
- `PlateMesh::{rectangle, rectangle_boundary}` — structured right-triangle
  meshes; no external mesh formats in v1 (deliberate).
- `dkt_stiffness` — the 9×9 DKT bending stiffness. DOF convention per node:
  `(w, wx, wy)` with wx = ∂w/∂x, wy = ∂w/∂y (SLOPES; Batoz's published
  θ-tables are translated internally, producing βx = −w,x, βy = −w,y whose
  sign cancels in the curvature quadratic form). Exact three-mid-side-point
  integration of the quadratic integrand.
- `assemble(mesh, section, boundary, stiffeners, opts) -> PlateModel` —
  reduced (K, M) pencil with Dirichlet elimination
  (`EdgeSupport::{SimplySupported, Clamped}`), lumped mass (translational
  ρhA/3, rotary ρh³/12·A/3 — SPD by construction, as the fs-modal pencil
  contract requires), optional uniform membrane pre-tension T entering as
  the P1 geometric stiffness on w, and `Stiffener` Hermite beams: bending
  with the parallel-axis effective rigidity `EI + EAe²` on (w, slope-along),
  torsion GJ on the cross-slope, lumped translational beam mass.
- `modes(model, window, opts)` — thin front over `fs_modal::slice_window`:
  every frequency arrives as a certified eigenvalue interval and the
  in-window count is inertia-certified.
- `PlateError` — typed refusals with stable `FS-PLATE-*` codes: bad
  section, degenerate element (with the offending element id and 2A), bad
  stiffener, forwarded modal refusals.

## Invariants
1. Element certificates (tested on an irregular triangle): stiffness
   symmetry; all three rigid-body motions store zero energy to roundoff;
   a constant-curvature field reproduces `A·κᵀDκ` EXACTLY — the patch test
   that pins the shape-function tables and the DOF convention (it caught a
   θ↔slope translation error during development).
2. Convergence to independent analytic references at two mesh densities
   with an O(h²)-or-better trend and ≤2% fine-mesh error: simply-supported
   Navier (measured 0.72% → 0.18%), clamped Leissa 35.992 (1.02% → 0.28%),
   orthotropic Navier rectangle (0.24%), continuum membrane in the
   prestress-dominated drumhead limit (0.41% → 0.10%).
3. Stiffener terms act at first-order Rayleigh scale (tested by term
   isolation — the gate that caught an asymmetric Hermite sign error which
   had inflated the fundamental 64%): weak-brace bending ≤ +1%, torsion
   ≤ +1%, lumped mass ≤ −1% on the quarter-line fixture; eccentricity
   multiplies a weak brace's rigidity visibly (measured 60.8 → 76.6 rad/s).
4. Orthotropy is load-bearing: swapping E_L/E_R moves the Navier
   fundamental by 89% on the spruce-like fixture (the wrong-grain-axis
   mutation gate).
5. All modal results flow through fs-modal: counts are inertia-certified
   and every eigenvalue carries its M⁻¹-norm residual interval.

## Error model
Structural misuse panics only inside test helpers; the public API refuses
through `PlateError`. Degenerate triangles refuse with the element id;
stiffener paths refuse on out-of-mesh nodes, zero-length segments, or
non-positive constants. Modal-stage refusals (singular window endpoints,
unresolved windows, non-SPD mass) forward the fs-modal typed error.

## Determinism class
Bit-deterministic across repeat runs on one host: sequential assembly
through `fs_sparse::Coo` canonical accumulation; modal path inherits
fs-modal's determinism (tested there). Cross-ISA goldens are not recorded.

## Cancellation behavior
None; assembly and modal calls run to completion (fs-modal budgets are the
only bounds). Joins the executor-integration seam when a consumer needs it.

## Unsafe boundary
None. Workspace `unsafe_code = "deny"`.

## Feature flags
None.

## Conformance tests
In-crate: element patch tests (symmetry, rigid-body, constant-curvature
exactness); SS Navier and clamped Leissa two-density convergence ladders
with trend assertions; orthotropic Navier + E_L/E_R swap mutation;
drumhead membrane limit; stiffener eccentricity discrimination in the
deliberately weak-brace regime (a stiff brace saturates into a line
support, where extra rigidity is invisible — measured and documented);
stiffener term isolation against Rayleigh scale; named refusals
(degenerate element, bad section, bad stiffener). All JSON-line evidence
rows; all modal counts inertia-certified through fs-modal.

## No-claim boundaries
- Bending only: NO membrane/in-plane DOFs, so no drilling stabilization is
  needed or present, and in-plane load paths (shear walls, arch action) are
  out of scope. Coupled membrane-bending (curved/arched shells, flat-facet
  shell assembly) is the recorded follow-up on the bead, triggered by the
  first curved-geometry consumer.
- MITC4 quads are NOT implemented (DKT triangles only); the bead retains
  that item. Structured rectangle meshes only in v1; no external mesh
  import.
- The stiffener is a straight Hermite beam on existing plate nodes with
  uniform cross-section: no variable cross-section (vibraphone undercuts),
  no curved braces, no beam shear deformation (Euler-Bernoulli, not
  Timoshenko in v1 — the eccentricity term dominates instrument bracing).
- The eccentricity model is the parallel-axis rigidity `EI + EAe²` acting
  through plate bending curvature; it does not model the membrane force
  the offset beam induces in the plate (requires membrane DOFs).
- Lumped mass only; consistent mass is a recorded follow-up (frequencies
  converge at the tested rates regardless).
- Uniform isotropic pre-tension only (scalar T); tensor/nonuniform
  prestress fields join the soundboard-downbearing consumer.
- No damping: pencils are (K, M); the viscoelastic bead supplies per-mode
  loss factors downstream.
- Validation is against analytic references (Navier, Leissa 35.992,
  membrane); no measured-instrument (holography) comparison yet — that
  case joins the instrument-matdb bead's consuming demo.
