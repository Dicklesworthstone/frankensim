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
  contract. Raw constructors take coherent SI (Pa, m, kg/m³); the `_qty`
  variants, `Stiffener::qty`, and `AssemblyOptions::qty` are the
  dimensioned fs-qty front doors (Length/Density/Pressure/Area/
  SecondMomentOfArea/SurfaceTension), bit-identical to the raw paths and
  pinned against a hand-computed spruce D-matrix (D11 = 27.2503 N·m,
  D22 = 2.04378, D12 = 0.715322, D33 = 1.6875 exactly). ν stays a bare
  ratio (dimensionless); the elastic law keeps fs-material's raw-Pa
  contract.
- `PlateMesh::{rectangle, rectangle_boundary, from_unstructured, structured_equivalent, boundary_nodes, boundary_edges, total_area}` —
  structured and unstructured right-triangle and general 2D triangular meshes;
  `from_unstructured` validates quality, coordinate finiteness, connectivity, and positive
  signed area (`2A > 1e-15`). `structured_equivalent` reproduces `rectangle` bit-identically.
- `PlateChart` and `PlateRegion` — mid-surface soundboard and plate chart representations
  binding geometry, thickness, material section, sub-regions, and boundary supports.
- `triangulate_soundboard` — parametric guitar/instrument soundboard mesh generator.
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
  boundary (with out-of-bounds node index and node count), bad stiffener,
  forwarded modal refusals.

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
4b. Stiffened-panel literature case: the Olson & Hazell (1977) clamped
   203×203×1.37 mm plate with one central 6.35×12.7 mm integral rib
   (E = 68.7 GPa, ν = 0.3, ρ = 2820 kg/m³; reference rows reproduced in
   Srivastava/Datta/Sheikh 2004 Table 4 and Thinh/Binh/Tu 2013 Table 1).
   The first five modes MAC-pair order-preservingly across two mesh
   densities (measured MAC ≥ 0.9997, floor 0.90) and sit within an
   authored 8% of the Olson–Hazell theory column (measured worst 1.7%,
   fine mesh: 712.2/762.9/981.7/997.6/1395.4 Hz vs
   718.1/751.4/997.4/1007.1/1419.8; experiment column recorded in the
   same JSON row, worst 5.2%). The envelope is authored against the
   model boundary: bending-only DKT + full-composite parallel-axis rib.
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

`tests/plate_conformance.rs` (bead 3ez8g.13.2 and 3ez8g.3.9; the reimplementation
contract through the public surface): pt-001 the SS isotropic plate
lands on the exact Navier frequency with mesh convergence (0.12%
fine); pt-002 the clamped square lands on Leissa 35.992 (0.28%);
pt-003 the fs-qty front door is BIT-IDENTICAL to the plain door
through assemble + modes; pt-004 stiffener term isolation (a
vanishing rib is a no-op at 3e-10, a real rib stiffens 56%); pt-005
the Olson-Hazell fundamental-region pin (709.6 vs 718.1 Hz theory,
1.2%; the full five-mode gate lives inline; corpus row
acoustic-olson-hazell-1977-mode3); pt-006 refusals by name — this
battery FOUND that an out-of-range boundary node PANICKED inside
assemble; it now refuses as `PlateError::BadBoundary`; pt-007
unstructured mesh admission (`PlateMesh::from_unstructured` /
`structured_equivalent`) reproduces rectangle meshes, boundaries, and
modal eigenvalues BIT-IDENTICALLY; pt-008 soundboard plate chart modes
on a parametric guitar-top geometry. Inline unit modules untouched.
In-crate: element patch tests (symmetry, rigid-body, constant-curvature
exactness); SS Navier and clamped Leissa two-density convergence ladders
with trend assertions; orthotropic Navier + E_L/E_R swap mutation;
drumhead membrane limit; stiffener eccentricity discrimination in the
deliberately weak-brace regime (a stiff brace saturates into a line
support, where extra rigidity is invisible — measured and documented);
stiffener term isolation against Rayleigh scale; the Olson–Hazell
stiffened-panel literature case with a MAC mode-pairing table (two mesh
densities, order-preserving pairing gate); named refusals (degenerate
element, bad section, bad boundary, bad stiffener). All JSON-line evidence rows; all
modal counts inertia-certified through fs-modal.

## No-claim boundaries
- Bending only: NO membrane/in-plane DOFs, so no drilling stabilization is
  needed or present, and in-plane load paths (shear walls, arch action) are
  out of scope. Coupled membrane-bending (curved/arched shells, flat-facet
  shell assembly) is a recorded follow-up (bead kj3s0 closed with this as
  a trigger-deferred item; tracked under 3ez8g.12.2), triggered by the first curved-geometry
  consumer.
- MITC4 quads are NOT implemented (DKT triangles only); trigger-deferred
  at the kj3s0 close (lands with the first quad-mesh consumer or a DKT
  validation gap). Arbitrary 2D planar triangulations are admitted via
  `PlateMesh::from_unstructured` (bead 3ez8g.3.9); 3D CAD surface/solid meshes
  are flattened/extracted into planar mid-surfaces before plate admission.
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
  membrane) plus the Olson–Hazell stiffened-panel literature case, whose
  experiment column (holography-measured) is recorded but gated only
  through the 8% theory-column envelope. No measured-INSTRUMENT
  (guitar-top holography) comparison yet — that case joins the
  instrument-matdb bead's consuming demo.
