# Fixed material regions in native 3-D design studies

Native `adaptive-simp` compliance studies can declare solid mounting regions and
zero-density regions that the optimizer must preserve. These are part of the
actual physical density map, not a display overlay or postprocessing repair.

```sh
cargo run --release -p fs-cli --features sdf3-study --bin frankensim -- \
  --json study examples/marquee/bracket-3d-regions.fsim regions-study.db
```

The example uses demonstration material/load values and a two-stage numerical
study, not a validated engineering specification.

## Input

Append this optional section **after the optimizer**, inside a version-1
`fsim-sdf3-study`:

```lisp
(design-regions
  :solid (((0.0 0.0 0.0) (0.5 0.5 0.5)))
  :void (((0.5 0.5 0.5) (1.0 1.0 1.0))))
```

Each box is a lower and upper coordinate vector in the study's SI frame. The
boundaries must coincide with its **initial** octree planes. For a unit box at
initial level 1, the planes are 0, 0.5 and 1; at level 2 they also include 0.25
and 0.75. Whole active cut cells are selected, not cells whose centers happen to
fall inside a partially covered region. At most 32 boxes are admitted; empty
solid/void lists are allowed. Reordered/unknown fields, zero or inverted spans,
nonaligned boundaries and overlapping solid/void interiors are refused. A box
that selects no active cells is also refused instead of being silently ignored.

Solid regions consume the material-volume allowance. If they alone exceed it,
the study refuses before equilibrium; ordinary feasibility restoration cannot
remove the prescribed material to manufacture an apparent optimization result.

## Numerical meaning

The physical map is `region_override(project(filter(raw)))`. Solid densities are
exactly one, void densities exactly zero, at every evaluated design and every
accepted continuation stage. Their local projection derivatives are zero before
the filter transpose is applied. Raw filter controls located in protected cells
may still influence neighboring free cells; raw controls are not physical density
fractions. The compliance and volume gradients differentiate this complete map.

Void means zero material in the **SIMP ersatz model**: its declared `e_min`
stiffness remains, and loads, supports and quadrature are not deleted. Use a
constructive implicit-domain difference for an actual geometric bore/cavity.
Neither this mask nor its background-cell boundary proves a manufacturing
clearance or continuum safety property.

Adaptive compliance refinement copies each protected label through the same
checked parent map used for raw density transfer. It does not reclassify a fine
grid by new cell centers. A stopped/rejected stage keeps the preceding geometry,
material labels, model and solved fields. The canonical source retains the box
policy, so the existing compliance replay-resume path re-admits the same policy.
Design exports contain the actual constrained projected densities and independent
load displacements, rather than unmasked candidate values.

Native stress-constrained minimum-volume studies currently **refuse nonempty
`design-regions`**; their newly introduced optimizer-state recovery needs a
separate constraint-preservation test before that combination is admitted.
Existing stress studies without this section are unchanged.

## Rust API and focused checks

`fs_topopt::sdf3::PhysicalRegion3::{Design, Solid, Void}` labels follow the
operator's active-cell ordering. Bind them with
`CutDensityStudy3::with_physical_regions` before numerical evaluation. Fixed-grid
OC/continuation and adaptive compliance continuation consume the same map.
Other code rebuilding a study must explicitly carry/rebind its physical policy.

```sh
cargo test --release -p fs-topopt --features cutfem-marquee --lib sdf3::regions
cargo test --release -p fs-cli --features sdf3-study --lib study::elasticity::sdf3::regions
```

The regressions require real independent CutFEM solves, coordinatewise gradient
comparisons, an actual constrained update, adaptive label inheritance, budget
refusal, cancellation and deterministic field replay. Test presence alone is not
an execution receipt.
