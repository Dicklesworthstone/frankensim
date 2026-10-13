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

The native path checks at most six deterministic axis/end directions. Extension
admission is shared, and there is at most one dual/majorant, after a successful
barrier. Cancellation, malformed inputs, numerical failures and resource
exhaustion propagate rather than becoming opportunities to retry a new axis.

The existing comparison requires an exterior on bounding-box planes, a fully
Dirichlet/positive-Robin cooling end, zero other Neumann flux and common uncoupled
axial conductivity. It also checks anchoring. The reusable library API preserves
constant tensors and permits tangential heterogeneity under that condition.
Finite-domain contact, radiation, feedback convection and temperature-dependent
conductivity remain unsupported. Previously supported global constant-source
and matching-contact bounds retain their original producer.

A failed proof does not establish physical invalidity. For example, a symmetric
slab's true maximum can be 302.625 K while this conservative one-ended barrier is
305.5 K. A 304 K material domain must therefore refuse this proof even though the
true solution is in range. A different justified comparison could admit it.
Likewise, all numerical nodes being at 300 K cannot justify a 305.4 K domain when
the independently established range reaches 305.5 K.

## Retained output

The conduction artifact's `volume_mean_bound.material_validity` identifies
`global-constant-source` or `continuum-slab-comparison`. The latter includes:

- `continuum_temperature_range_k`: outward whole-domain PDE range;
- `retained_material_domain_k`: common original table domain;
- `axis`, `reverse`, `comparison_attempts`: the successful comparison path.

Each material row retains its card, selected claim, receipt hash and exact
receipt bytes, plus `temperature_domain_k` (null only for the global path).
The primary `enclosure_k` is still the selected volume mean, not the nodal
maximum. Normal physical fields and energy reports are unchanged. The mean
continues through the existing report and evidence-package path; unknown
material, geometric and physical uncertainty remains unknown. No complete
engineering verdict or experiment is certified by this numerical result.

## Focused verification

```sh
cargo test -p fs-conduction --features thermal-verification --test bounded_material_mean
cargo test -p fs-cli --features thermal-verification --test bounded_native_mean --test native_mean_bound
cargo test -p fs-cli --features thermal-verification --lib material_range::tests
```

Five core tests, three actual-command tests and two source-admission tests cover
signed heating, reversed geometry, tensor assignments, preserved fields and
receipts, finite/global parity, report/package integration, conservative
refusals, and original budgets/cancellation. Native tests were authored but could
not run in the authoring environment because Cargo/rustc were absent. Independent
Python checks of 24 one-dimensional slab discretizations agreed with closed-form
nodal temperatures within 6.37e-12 K; those are mathematical references, not
execution of the committed Rust implementation or an experimental validation.
