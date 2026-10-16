# Native mean-temperature bounds with finite material validity

The existing native `.fsim` `temperature-volume-mean-bound` report can now use
constant scalar material cards with a finite temperature domain. It does not
silently treat those cards as globally valid, use the nodal minimum/maximum as a
continuum proof, or change the normal conduction solve.

Build the CLI with `fs-cli/thermal-verification` and add the existing report
request alongside the original maximum objective:

```lisp
(outputs
  (qoi :name "temperature-max" :kind "scalar")
  (qoi :name "temperature-volume-mean-bound" :kind "report"))
```

Use the ordinary import/solve/report/package commands with the original mesh,
project, ledger and immutable material packs. There is no new command, project
schema, dependency, material coefficient or input tolerance.

## Admission and physical meaning

Finite-domain cards must contain the selected scalar `ConstantWithinValidity`
conductivity claim with temperature as their only validity axis. The original
native table, its narrower authored temperature span, and its property-use
receipts are retained. Flat sampled source data is not promoted to this source
claim, and an unbounded native copy cannot erase a finite source domain.

The numerical producer constructs a private mathematical constant extension,
then applies `fs-verify`'s existing outward slab comparison to that equation.
Only if the entire continuum range is inside every used material's original
domain can the regional dual and equilibrated flux majorants run. The extended
solution then also solves the original finite-domain constant law. No physical
primal solve uses the extension, and no extrapolated property receipt is made.
The supplied approximate field's algebraic error remains in the mean bound.

The native path still checks at most six deterministic comparisons. When both
ends of a proposed axis are entirely prescribed-temperature/positive-Robin
faces, it uses both ends rather than discarding one end's cooling. Unequal
coefficients and mixed prescribed/Robin end patches are supported; each end
retains its own weakest Robin coefficient. Reversing this two-ended comparison
does not spend a duplicate attempt. Partly insulated ends retain the original
one-ended comparison in the requested direction. Extension admission is shared,
and there is at most one dual/majorant, after a successful barrier. Cancellation, malformed inputs, numerical failures and resource
exhaustion propagate rather than becoming opportunities to retry a new axis.

The existing comparison requires an exterior on bounding-box planes, a fully
Dirichlet/positive-Robin cooling end, zero other Neumann flux and common uncoupled
axial conductivity. It also checks anchoring. The reusable library API preserves
constant tensors and permits tangential heterogeneity under that condition.
Finite-domain contact, radiation, feedback convection and temperature-dependent
conductivity remain unsupported. Previously supported global constant-source
and matching-contact bounds retain their original producer.

For a 0.1 m slab with k=10 W/(m K), both ends cooled with h=20 W/(m2 K)
to 300 K, and uniform source 1000 W/m3, the old one-ended comparison reaches
305.5 K. The two-ended comparison reaches 302.625 K (with outward rounding),
so an unchanged 303 K material limit can now be admitted. For equal end data
this equals the analytic maximum; general nonuniform source/cooling cases can
still have conservative bounds. A failed proof is still not evidence that the
physical solution leaves its validity domain.

The new `slab_temperature_range_two_sided` numerical producer solves the
asymmetric unit-source comparison exactly before outward interval evaluation.
If r0 and r1 are reciprocal end coefficients (zero at an all-prescribed end),
q0=L*(r1+L/(2k))/(r0+r1+L/k) is the left heat share per unit source. The barrier
maximum is r0*q0+q0^2/(2k). All original source and boundary extrema enter the
comparison; it never substitutes their averages or assumes physical symmetry.
The public one-ended verifier method is unchanged.

## Retained output

The conduction artifact's `volume_mean_bound.material_validity` identifies
`global-constant-source` or `continuum-slab-comparison`. The latter includes:

- `continuum_temperature_range_k`: outward whole-domain PDE range;
- `retained_material_domain_k`: common original table domain;
- `axis`, `reverse`, `comparison_attempts`: the successful comparison path;
- `two_sided`: whether both ends were used (`reverse` matters only otherwise).

`SlabValidatedRegionMean` also gained `two_sided`; code constructing that public
report directly must supply it. No solver input or project schema changed.

Each material row retains its card, selected claim, receipt hash and exact
receipt bytes, plus `temperature_domain_k` (null only for the global path).
The primary `enclosure_k` is still the selected volume mean, not the nodal
maximum. Normal physical fields and energy reports are unchanged. The mean
continues through the existing report and evidence-package path; unknown
material, geometric and physical uncertainty remains unknown. No complete
engineering verdict or experiment is certified by this numerical result.

## Focused verification

```sh
cargo test -p fs-verify --features certified-speculation --lib two_sided
cargo test -p fs-conduction --features thermal-verification --test bounded_material_mean
cargo test -p fs-cli --features thermal-verification --test bounded_native_mean --test native_mean_bound
cargo test -p fs-cli --features thermal-verification --lib material_range::tests
```

Four additional standalone verifier tests cover asymmetric and mixed end laws,
signed sources, partial coverage, rotation, and numerical/resource refusals.
Two additional library tests cover previously unprovable finite domains and
bounded/deduplicated direction selection. The existing actual-command tests
now require the two-ended proof for a 301 K native cube domain and retain
unchanged fields, source receipts, report/package export and real refusals.

Cargo/rustc are absent in the authoring environment; compilation and native
execution are not claimed here. Independent Python reference checks passed
400 exact-rational interval-containment comparisons and 120 polynomial/FEM
comparisons, with maximum nodal discrepancy 4.55e-11 K. These are reference
mathematics, not execution of the Rust code or physical validation.
