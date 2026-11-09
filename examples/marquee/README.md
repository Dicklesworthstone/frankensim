# Native design studies

## Self-weight in 3-D compliance and stress studies

The `fsim-sdf3-study` compliance and stress modes accept a common fixed
acceleration such as gravity. Append this section after any `design-regions`
and `linear-solver` sections, inside the study:

```lisp
(material-acceleration
  :density-kg-m3 2700.0
  :acceleration-m-s2 (0.0 0.0 -9.81)
  :max-stencil-terms 100000)
```

The declared density is full-material reference mass per cubic metre. Actual
self-weight follows filtered/projected material occupancy: protected solids
retain full weight and protected voids have none, including on refined grids.
Each independent load case includes this common acceleration in addition to
its existing body/pressure/traction dead loads. Zero dead-load vectors are
valid for a pure self-weight study. The native compliance and stress adjoints
differentiate the changing force as well as stiffness. Stress resume checks
the rebuilt material-load identity before restoring optimizer state.
Accelerated studies retain the native material's strain-domain checks, so the
declared elastic modulus and SI loads must describe a valid small-strain state.

The acceleration stays fixed in the reference configuration; this does not
add follower loads or inertial dynamics. `max-stencil-terms` bounds each grid's
prepared sparse load map and contributes to memory admission. The compliance
OC optimizer still requires nonpositive sensitivities; self-weight can
violate that condition even at higher SIMP penalization and is refused by
its existing gradient gate. The stress optimizer has no such sign restriction. The
zero-centered `load-space-limited-simp` mode explicitly refuses this common
offset because its force-ellipsoid model does not include affine loading. Use
the `affine-load-space-limited-simp` mode below to combine self-weight with a
continuous external-load uncertainty set.

## Nearly incompressible 3-D materials

Add `:formulation bbar` as the last field of the `physics` section to select the
native B-bar volumetric-strain formulation. This admits Poisson ratios through
0.4999 and propagates through bulk, surface-loaded, embedded-support,
material-acceleration and refined-grid operators. Omitting the field preserves
the existing displacement formulation and its Poisson-ratio limit of 1/3.
The formulation is included in the study identity; it does not add nonlinear
material response or a continuum error bound.

## Robust 3-D minimum-material design

`bracket-3d-load-space.fsim` uses the existing 3-D CutFEM, density filter,
three-field projection and augmented-Lagrangian optimizer to minimize dilated
material while every simultaneous force combination `f = Fq`, `||q||₂ <= 1`,
meets the declared compliance cap on eroded material. Build with
`--features sdf3-study`, then run:

```sh
frankensim --json study examples/marquee/bracket-3d-load-space.fsim robust.db
frankensim --json report study-<receipt-hash> robust.db
frankensim --json package study-<receipt-hash> robust.db
```

Each declared load is a column of `F`; its SI force components contain its
amplitude and correlation with other columns. Weights must be exactly one.
The full cross-compliance matrix determines the worst combination. Checking
only the separately applied columns would miss correlated combinations.
`compliance-limit-j` is a physical upper bound in joules; the dimensionless
`spectral-temperature` controls the conservative smooth maximum used by the
optimizer. `threshold-offset` shifts the eroded/dilated projection thresholds;
it does not specify a geometric tolerance in metres.

The report distinguishes the last accepted optimizer point from the accepted
feasible point with least dilated material. It retains the worst combination,
cross-compliance spectrum, separate material fields, actual displacement
fields, initial directional gradient checks and original work consumption.
Accepted prefixes are durable and exportable. `--budget N` limits updates in
this invocation; an incomplete run exits 6. `study --resume` restores the
accepted optimizer state and feasible archive under the original budgets,
rebuilding geometry and verifying retained endpoint responses before continuing.
`max-optimizer-evaluations` limits the optimizer's evaluations; the initial
gradient probes additionally consume the same total Krylov and wall budgets.

These results are numerical estimates for the declared linear-elastic model.
They do not establish continuum error, stress/buckling safety, global
optimality, physical validation or a manufacturing tolerance. The example's
modulus and geometry are numerical demonstration inputs, not a calibrated part.

## Minimum material with self-weight and uncertain external loads

`bracket-3d-affine-load-space.fsim` extends the design problem to
`f = B*rho_eroded + Fq`, `||q||₂ <= 1`. Material acceleration supplies the
common self-weight, which changes as material is removed. The authored load
columns describe simultaneous uncertain external forces. They have unit
weights, while their SI components carry amplitudes and correlations.

```sh
frankensim --json study examples/marquee/bracket-3d-affine-load-space.fsim affine.db --budget 2
frankensim --json report study-<receipt-hash> affine.db
frankensim --json package study-<receipt-hash> affine.db
```

Replace the placeholder with the returned run ID. The first command limits the
run to two updates; omit `--budget 2` to use the declared update allowance.
The native affine maximizer accounts for cross terms between self-weight and
external loads when finding the worst compliance. This mode has no
`spectral-temperature` parameter. The report records the actual worst-load
coefficients, offset, numerical gap and residuals, initial gradient checks,
and separate last-accepted and least-material feasible designs.

Report export writes the retained eroded/nominal/dilated density fields and
offset, basis and worst displacements to VTU. An infeasible selected design
does not produce a field export. Incomplete runs retain their accepted results
and exit 6. Continue an accepted prefix with:

```sh
frankensim --json study --resume study-<receipt-hash> affine.db --budget 2
```

Resume restores the accepted augmented-Lagrangian state, its history and the
least-material feasible incumbent under the original source and executable.
It rebuilds geometry and re-solves the accepted and distinct incumbent endpoints
to check their physical responses before continuing. Those evaluations spend
the original optimizer, Krylov, geometry and wall allowances; earlier search
updates and the initial gradient experiment are not repeated. The new budget
adds updates without extending the original target. Historical receipts written
before checkpoint support remain exportable but cannot resume.

Stiffness and self-weight both use eroded occupancy, while the objective uses
dilated volume. Under self-weight this does not establish erosion as the worst
manufacturing realization. The example uses numerical demonstration inputs
and makes the same limited linear-elastic claims as the robust example above.

## Worst-case stiffness with a material-volume budget

`bracket-3d-affine-compliance.fsim` uses the same native self-weight, uncertain
load ball, strength, signed responses and grid-refinement workflow to minimize
worst compliance subject to a **dilated-volume fraction cap**. It declares
objective `compliance` in joules and optimizer `affine-load-space-simp`.
The ordered optimizer fields add `volume-fraction` and `compliance-reference-j`
in place of the minimum-material mode's `compliance-limit-j`.

```sh
frankensim --json study examples/marquee/bracket-3d-affine-compliance.fsim stiffness.db --budget 2
frankensim --json study --resume study-<receipt-hash> stiffness.db --budget 2
frankensim --json report study-<receipt-hash> stiffness.db
```

The compliance reference normalizes the objective and its gradients. It is
**not a compliance cap**: a design can exceed that reference while satisfying
its material-volume, strength and signed-response requirements. The native
feasible archive selects the lowest-compliance accepted design, reported as
`least-compliance-feasible`; VTU exports that design's actual physical fields.
Reports retain both physical objective values in joules and normalized values
used by the optimizer. An infeasible initial material inventory is allowed;
only designs passing all original constraints enter the feasible archive.

Resume, candidate rechecks and optional finer-grid redesign retain the same
declared objective and volume cap. Grid assessment separately gates compliance
objective change with `compliance-change-tolerance` and absolute volume-fraction
change with `inventory-change-tolerance`. The example remains a numerical
linear-elastic demonstration, without a continuum stiffness or optimality claim.

## Known service preload, self-weight and uncertain loads

`bracket-3d-affine-preload.fsim` adds a fixed mechanical preload to the
fixed-volume compliance example. Its complete force family is
`f = f_dead + B*rho_eroded + Fq`, `||q||₂ <= 1`. The optional section follows
`material-acceleration` and precedes strength, response and grid controls:

```lisp
(affine-dead-load
  :body-n-m3 (0.0 0.0 -0.1)
  :surface (pressure :pa 0.1 :x-fraction (0.5 1.0)))
```

Body force density uses N/m3; surface pressure and traction use Pa. To author
traction instead, use `(traction :pa (0.0 0.0 -0.1) :x-fraction (0.5 1.0))`.
Use `surface none` for a body-only preload. The fixed preload does not scale
with material occupancy. Self-weight still follows eroded occupancy, and the
uncertain columns still describe simultaneous deviations around the combined
common load. Omitting the section retains the original zero fixed preload.

```sh
frankensim --json study examples/marquee/bracket-3d-affine-preload.fsim preload.db --budget 2
frankensim --json study --resume study-<receipt-hash> preload.db --budget 2
frankensim --json report study-<receipt-hash> preload.db
```

Pressure acts inward on the reference implicit surface, including cavity walls
inside the authored x band; it is not a follower load. Surface bands must align
with the initial octree and cannot overlap an embedded support. Every enriched
or refitted grid reintegrates the same physical reference law. Reported offset
displacement includes both fixed preload and self-weight, with their cross terms
retained in the native worst-load calculation and every declared constraint.

## Joint compliance and strength with uncertain loads

`bracket-3d-affine-strength.fsim` adds a physical stress allowable to the
self-weight/load-ball problem. Run and continue it with the same commands:

```sh
frankensim --json study examples/marquee/bracket-3d-affine-strength.fsim strength.db --budget 2
frankensim --json report study-<receipt-hash> strength.db
frankensim --json study --resume study-<receipt-hash> strength.db --budget 2
```

The final `affine-strength` section declares `stress-limit-pa`,
`relaxation-power`, `aggregation-power`, `max-stress-point-fields` and
`joint-sharpness`, in that order. The example uses a numerical demonstration
allowable of 40 Pa. Choose the allowable from the intended material model;
the modulus does not determine it.

Strength is evaluated over the entire affine load ball at every retained bulk
quadrature point. Each point has its own adverse load coefficients, which can
differ from the worst-compliance coefficients. The conservative p-norm includes
both relaxed and physical stress envelopes. Native primal fields are shared;
stress adjoints scale with the number of load-basis fields rather than with
the number of sample points. Both compliance and strength must pass before a
design enters the least-material feasible archive.

Reports show selected and last-accepted stress aggregates, physical and
relaxed governing witnesses, numerical oracle diagnostics and the complete
joint constraint. Resume retains the strength policy and initial joint
gradient evidence under the original budgets. The point-field allowance
counts the common offset plus all uncertainty columns at every sampled point;
it is not a geometric point count. The density floor must remain above the
native stress foldback threshold after eroded projection.

For a selected feasible design, report export also writes cell-associated
`sampled_physical_stress_envelope` and `sampled_relaxed_stress_envelope` arrays
in Pa, with integer `stress_sample_count`. These are maxima over retained
quadrature-point envelopes in each background cell. They may use different
load directions within one cell; the existing `worst_displacement` field is
the compliance witness. Export reads the sealed design without another solve,
and refuses incomplete declared stress evidence.

These are fixed-grid sampled strength estimates. They do not bound the
continuum stress maximum or certify safety, fatigue, buckling or manufacturing
tolerance. The example's geometry, stiffness and allowable are demonstration
inputs, and its feasibility must be assessed from the actual computed report.

## Minimum material with displacement and support limits

`bracket-3d-affine-responses.fsim` combines compliance, sampled strength and
signed serviceability limits over the same self-weight and uncertain-load
family. It declares a free-side mean z-displacement window, a z support-force
window, and an upper bound on the y support moment about the authored origin.

```sh
frankensim --json study examples/marquee/bracket-3d-affine-responses.fsim responses.db --budget 2
frankensim --json report study-<receipt-hash> responses.db
frankensim --json study --resume study-<receipt-hash> responses.db --budget 2
```

The optional `affine-responses` section declares `sharpness`, `joint-sharpness`
and an ordered `limits` list. Each `mean-displacement`, `support-force` or
`support-moment` row has a unique name, Cartesian axis, signed `at-most` or
`at-least` sense, original bound and positive same-unit scale. Displacement
uses `bound-m`/`scale-m`, force uses `bound-n`/`scale-n`, and moment uses
`bound-n-m`/`scale-n-m`. The example demonstrates both sides of a window;
negative and zero bounds retain their ordinary signed meanings.

Mean displacement uses the fixed reference cut-volume average in the declared
initial-grid-aligned `x-fraction` band. Force and moment are native numerical
tractions on the solid's embedded support; moments use the authored
`about-m` origin. The support reaction remains material dependent during
optimization. Each inequality has its own adverse load coefficients and
remeasured combined-field witness. Reports show all original bounds and their
physical units, with selected and last-accepted values. Resume preserves the
full constraint family and original gradient evidence.

Report export includes one `response_displacement_<index>_<name>` VTU vector
field, in metres, for every original signed bound. The fields combine the
retained common offset and basis displacements with each bound's own adverse
coefficients. A displacement window consequently has two governing fields.
The physical force and moment values and their N/N*m units remain in the
report. Export checks the entire bounded expansion and requires no new solve.

The implementation shares offset/basis primal fields with compliance and
strength. Response adjoints are bounded by the smaller of the number of
limits and the number of affine fields. A response with exactly zero
uncertainty radius has no smooth norm derivative and is refused for smooth
optimization. The example's limits are demonstration inputs whose feasibility
must be read from the computed report; no continuum serviceability or physical
support-capacity certification is implied.

## Joint design with an enriched-grid physical recheck

`bracket-3d-affine-grid-assessment.fsim` adds an actual finer-grid assessment
to the complete self-weight, uncertain-load, strength and signed-serviceability
problem. The optimizer retains its initial grid. After saving its selected
design, the study refines every active cell once, transfers physical occupancy
through the native parent map, rebuilds self-weight and re-solves the whole load
family on both grids. Each displacement average uses its own grid's measured
reference cut volume over the original physical band.

```sh
frankensim --json study examples/marquee/bracket-3d-affine-grid-assessment.fsim assessed.db --budget 2
frankensim --json report study-<receipt-hash> assessed.db
frankensim --json study --resume study-<receipt-hash> assessed.db --budget 2
```

The final `affine-grid-assessment` section declares the transfer allowance and
separate compliance, response, stress and inventory change tolerances. The
existing physics `maximum-level` and `maximum-leaves` must accommodate the
one-level probe; the example raises them from 1/8 to 2/64 and uses a 1 GiB
memory admission allowance. Every solve, geometry sample and elapsed second
spends the original cumulative budgets. It adds no optimizer trial evaluation.

Read `grid_assessment.passed` together with `optimization_status`. A pass
requires both grids' original physical bounds and joint aggregate, complete
load/response/stress model stability, and stable dilated inventory. A completed
optimizer can therefore report `grid-unresolved`. Failed or interrupted probes
retain their reason and spent work, and can be retried through normal resume
without another AL step after optimization is finished. The report includes
both grids' adverse coefficients, original signed bounds, sampled stresses,
normalization volumes and frozen-load DWR/switch diagnostics. Exported design
and VTU fields remain those of the selected initial-grid design. This example
checks numerical grid dependence; it does not claim a mesh-independent optimum
or a continuum strength/serviceability certificate.

## Automatic finer-grid redesign of the joint affine problem

`bracket-3d-affine-grid-refinement.fsim` extends the assessed joint design with
one automatic native refit. After the source optimizer converges or finishes
its declared update allowance, a feasible selected design with a checked but
unresolved grid assessment starts a finer-grid candidate. The candidate inherits
raw density controls and protected regions, refilters on its own grid, rebuilds
material self-weight and preserves every original physical constraint.

```sh
frankensim --json study examples/marquee/bracket-3d-affine-grid-refinement.fsim refined.db --budget 2
frankensim --json study --resume study-<receipt-hash> refined.db --budget 2
frankensim --json report study-<receipt-hash> refined.db
```

The final `affine-grid-refinement` section declares candidate `max-updates` and
`max-optimizer-evaluations`. The invocation `--budget` caps additional accepted
updates across both source and candidate; it never renews either optimizer's
allowance. Geometry, linear-solver work and wall time remain cumulative. The
example reserves level 3 and 512 leaves for the candidate and its own subsequent
probe, with a 2 GiB memory admission allowance.

`grid_refinement.phase` distinguishes `preparing`, `candidate` and `promoted`.
Throughout preparation and candidate optimization, the selected `design` and
VTU remain exactly the original source design. `working_design` and its native
checkpoint retain the finer candidate separately. Resume restores that candidate's
accepted state, feasible archive and original gradient audit without repeating
its earlier updates. Reports keep the two grids' histories and quotas separate.

Promotion requires a completed candidate allowance and a fresh complete-family
recheck passing compliance, every signed response, sampled strength and the joint
constraint. The promoted design then receives a new finer-grid assessment;
promotion alone does not establish grid convergence. Read that assessment and
`optimization_status` alongside the phase. The CLI attempts one refit, retains
stopped prefixes for bounded continuation and does not claim cross-grid material
descent or a continuum optimum.

## Thermal cooling-hole radius study

`thermal-2d.fsim` runs normalized scalar Poisson physics through the CLI,
optimizer, ledger, report and evidence package. It solves `-laplacian(u) = 1`
on the unit square with two circular holes, with `u = 0` on every boundary.
The objective is dimensionless thermal compliance, the integral of `u` over
the material. Projected gradient steps adjust hole radii at fixed centers and
material area, using the existing bounded Armijo line search.

`bracket-2d.fsim` records the proposed elasticity study. Its physics is still
unimplemented; it is not an executable thermal example or evidence that
elasticity/topology optimization works.

With a built `frankensim` binary:

```sh
frankensim --json study examples/marquee/thermal-2d.fsim study.db
frankensim --json report study-<receipt-hash> study.db
frankensim --json package study-<receipt-hash> study.db
```

The study result supplies the full `study-<receipt-hash>` ID. Reports and
packages are exported next to the ledger; an existing file with different
bytes is refused rather than overwritten. The report contains actual accepted
compliance values, convergence SVG, DWR/algebraic estimates and final area.

For a short invocation followed by continuation:

```sh
frankensim --json study examples/marquee/thermal-2d.fsim study.db --budget 2
# Exit 6 is an explicit partial. Use the returned receipt ID:
frankensim --json study --resume study-<receipt-hash> study.db
```

Resume verifies and replays the saved prefix before continuing. It retains the
original wall and total-iteration budgets, and charges replay time as well.
The invocation cap limits new steps; it does not raise the study's grant.

Every section in the example has a concrete meaning:

| Section | Meaning |
| --- | --- |
| root/version | Study syntax version 1; executable input must preserve its complete canonical declaration. Whitespace/comments may differ. |
| metadata | Advisory scoping context; this estimated study cannot make a signoff decision. |
| versions | Declared schema. The running driver separately retains its crate version and embedded constellation-lock digest. These are provenance, not authenticated source identity. |
| seeds | Explicit root seed, retained in operations. The radius algorithm itself does not sample randomness. |
| budgets | Wall seconds, memory admission value and maximum total iterations. Memory is not a measured RSS ceiling. |
| capabilities | Explicit thermal-radius optimization, SDF geometry and CutFEM capabilities. |
| units | `normalized`; this problem does not claim a dimensional physical calibration. |
| domain | Unit square and strictly interior, disjoint circular holes. |
| physics | `thermal-poisson-2d-normalized`, with a fixed background quadtree level. |
| scenario | Unit source throughout the material, zero temperature on all boundaries. |
| objective | Minimize normalized compliance in unit `1`; joules or watts are refused. |
| constraints | Material-area equality in `(0,1)`, feasible under the declared radius bounds. |
| optimizer | Projected gradient, initial step size, radius bounds, iteration count. |

The supported execution envelope is at most mesh level 5, 32 holes and 256
iterations, with a memory admission value of at least 128 MiB. Cancellation and
wall enforcement occur between bounded iterations, not inside individual
CutFEM solves. No operating-system signal handler is installed by the CLI.

Outputs are **Estimated**. The package checker verifies structural integrity;
it does not establish a guaranteed discretization bound, physical validation,
KKT conditions or global optimality. This is fixed-center radius optimization,
not free-boundary topology or elasticity. The retained geometry consists of
exact circle/plate parameters rather than a sampled SDF grid.

The existing `scripts/e2e/marquee_01.sh` lane invokes the actual executable for
this slice. Its September 4 Linux development-build run passed eight checks:
full study, retained SVG report, checked package export, partial budget,
identical resumed trace, and wrong-unit/load/area refusals before database
creation. The [implementation evidence](../../docs/REALITY_CHECK_2026-09-01.md)
records the binary, source scope and artifacts. Release/DSR and full flagship
acceptance remain open.
