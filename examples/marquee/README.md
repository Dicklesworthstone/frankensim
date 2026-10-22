# Native design studies

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
this invocation; an incomplete run exits 6. This new fixed-grid mode does not
yet restore optimizer state across processes, so `study --resume` refuses it.
`max-optimizer-evaluations` limits the optimizer's evaluations; the initial
gradient probes additionally consume the same total Krylov and wall budgets.

These results are numerical estimates for the declared linear-elastic model.
They do not establish continuum error, stress/buckling safety, global
optimality, physical validation or a manufacturing tolerance. The example's
modulus and geometry are numerical demonstration inputs, not a calibrated part.

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
