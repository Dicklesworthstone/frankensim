# Preserve mounting material and clearance regions

`fs-marquee-elasticity-robust --projected` accepts `--design-regions REGIONS.csv`
for NEW studies. This authors regions that must remain material or empty during
area-constrained optimization. It uses the existing fixed-node projection,
checkpoint, recovery and refinement paths; no second geometry-evolution or
elasticity solver is introduced. Stress constraints and optional stress
restoration continue to use their original complete-load-family gates.

```sh
cargo run --release -p fs-marquee --bin fs-marquee-elasticity-robust -- \
  --projected /tmp/protected-design examples/marquee/protected-loads.csv \
  4 30 0.65 6 sum 362 \
  --design-regions examples/marquee/protected-regions.csv --checkpoint
```

The supplied two-load example reserves a material attachment region and an
interior empty region. Native execution of this example is unverified. Reaching
an iteration limit is not an optimality claim; a candidate search can stall.
Use the existing stress options only with a bound appropriate to the model.

## File format and precise geometric meaning

Blank lines and `#` comments are allowed. Every other line has six fields:

```csv
# phase,x_min,y_min,x_max,y_max,phi_margin
material,0.125,0.375,0.25,0.625,0.02
void,0.5,0.4375,0.5625,0.5625,0.02
```

`phase` is exactly `material` or `void`. Bounds must be finite, strictly ordered
in each axis, and inside the normalized `[0,1]²` domain. The finite positive
margin is measured in level-set field units, NOT a distance or a stress unit.
The CLI admits at most 64 regions in a file of at most 1 MiB. The library
supports dyadic lattices with 2–256 cells per side; the executable retains its
existing stricter level limits. Malformed, empty, conflicting and oversized
inputs refuse rather than producing an unconstrained study.

Every cell with a positive-area intersection with a rectangle is covered, and
ALL of its corner values are fixed with the declared sign margin. The exact
bilinear interpolant is then negative throughout a material cell and positive
throughout a void cell. A sub-cell hole cannot disappear because no cell centre
happened to sample it. The covered rectangle may extend beyond the authored
one by less than one cell per side. The reported integer `covered_cells` bounds
are `[i_begin,j_begin,i_end,j_end]` with exclusive ends. These describe the
actual enforced support, not an invented sub-grid boundary.

Same-phase overlaps use the strongest requested margin. Opposite phases that
share any required node are incompatible at that grid resolution, including
adjacent rectangles that only share an edge. Such requests refuse; separate
regions or choose a finer input grid. Existing prescribed boundary values must
already satisfy the requested margin and are never weakened or overwritten.
Unconstrained nodes preserve their original bits. Region-controlled nodes keep
an existing stronger same-sign value or receive the minimum required change.

The union of existing fixed nodes and region nodes becomes the original
optimizer's fixed geometry. Area projection must still reach the declared
material target, and all loads must remain admissible. For example, protecting
too much material or cutting through a prescribed load can make a study
infeasible. No area target, load, stress limit or boundary condition is relaxed
to hide that conflict. Geometry authoring is not credited as optimization.

## Outputs and continuation

The existing `input-level-set.csv` remains the ORIGINAL supplied field.
`design-region-level-set.csv` is the authored field before area restoration.
`baseline-level-set.csv` is the area-restored field independently solved under
every load. `design-regions.csv` is a canonical, reusable copy of the input
records. The summary's optional `design_regions` object records coverage,
margins, node counts and the number of changed input values. The summary is
still written last; an export failure does not print success.

Existing checkpoints already retain exact prescribed nodal values. Thus
`--resume` needs no region file, cannot replace these constraints, and does not
repeat the authoring operation. `--refine` transfers fixed whole cells and their
signs using the existing bilinear fixed-node rule. Keep the initial study's
region CSV/summary for the authored rectangle descriptions: subsequent segments
retain numerical fixed-node state but do not fabricate the original input file
or repeat its authoring summary. No checkpoint format change is required.

Runs without `--design-regions` do not perform region authoring and retain their
existing output fields and numerical path. The library's controlled preparation
polls at rasterization/staging rows and before publication; interrupted work
returns no partially authored geometry and leaves inputs untouched.

## Scope and verification

This is a discrete non-design-region capability. It is NOT minimum wall-thickness
control, tool-access analysis, a machinability certificate, material-strength
validation, a physical clearance-distance bound, or a continuum stress bound.
Bilinear geometry is retained; no body-fitted mesh is created.

Seven core tests, two parser tests and four real-binary integration tests cover
whole-cell and sub-cell support, conflicting regions, fixed-boundary refusal,
order independence, cancellation, area infeasibility, actual accepted updates,
stress-enabled continuation and fine-grid baseline preservation. These are
AUTHORED, NOT EXECUTED: Rust/Cargo, rustfmt, DSR and RCH were unavailable in the
implementation environment. Independent geometry-reference and source checks
are not execution of the Rust implementation or an elasticity study.

```sh
cargo test -p fs-topols --lib design_regions
cargo test -p fs-marquee --bin fs-marquee-elasticity-robust --test design_regions
```
