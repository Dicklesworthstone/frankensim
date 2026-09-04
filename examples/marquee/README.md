# Thermal cooling-hole radius study

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

The existing `scripts/e2e/marquee_01.sh` lane is being wired to this executable
slice; use its recorded runtime result, not its file presence, as evidence.
