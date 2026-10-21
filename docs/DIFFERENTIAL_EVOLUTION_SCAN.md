# Population search for thermal designs with flat phase responses

`fs_dfo::differential_evolution::DifferentialEvolutionState` implements
bounded **DE/rand/1/bin** (Storn and Price, 1997,
DOI `10.1023/A:1008202821328`). It fills the population from the declared
finite coordinate box, optionally retaining an exact initial point as member
zero. It does not need objective derivatives.

This is useful when liquid fraction is zero over a range of heating powers:
a local simplex may shrink entirely inside that flat region. Small local
steps then report no improvement even though larger, physically valid powers
would melt the material. Population initialization explores the whole box;
it is not a guarantee of coverage or recovery.

## Physical example

The existing `simplex_scan_fit` executable now offers either algorithm over
the **same** `ScanPowerStudy::evaluate_value` and P1 enthalpy model:

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example simplex_scan_fit -- \
  --population --initial-enthalpy 80 --seed 17

cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example simplex_scan_fit -- \
  --initial-enthalpy 80
```

The first selects the population method; the second is the local-simplex
comparison. The original invocation without options still uses Nelder-Mead
and starts at the solidus, 100 J/kg.

The explicit synthetic specimen is a unit cube, density 1 kg/m³, conductivity
2 W/(m K), and latent heat 100 J/kg. It uses six P1 tetrahedra and eight
vertices. Initial enthalpy 80 J/kg corresponds to 308 K; the solidus is 310 K.
Targets are volume-mean liquid **mass** fractions 0.1 and 0.3 at one and two
seconds, and the source position changes at one second.

Both methods invoke the original nonlinear enthalpy solver, source-power
normalization, spatial measurement, event-left sampling and physical energy
gates. The new option does not insert a lumped thermal surrogate, a penalty
for failed physics or a numerical gradient. Its two-sided derivative tape
allowance stays zero.

The default population method uses 24 members, F=0.8, CR=0.9, seed 17, and a
3,000-evaluation total allowance, including one reserved final physical replay.
Its power box is [0,50] W for each interval. Override targets, power maximum,
evaluation allowance and seed explicitly:

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example simplex_scan_fit -- \
  --population --initial-enthalpy 80 --seed 17 0.1 0.3 50 3000
```

The sensible-plus-latent mean balance for a trajectory that reaches the
plateau at all nodes gives 30 W then 20 W. The executable checks actual
spatially averaged fractions, not this analytical power pair. A fresh
value-only evaluation of the exported design must reproduce its score
bit-for-bit; unselected shutdown remains exactly zero. It prints actual
unsnapped powers, predictions, remaining loss, method, stop reason and counts.
Failure to meet the fixed objective target produces an unsuccessful exit,
even when the optimizer spent its whole budget.

The intentionally underpowered case must not claim recovery:

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example simplex_scan_fit -- \
  --population --initial-enthalpy 80 0.1 0.3 1 65
```

## Library use

```rust
use fs_dfo::differential_evolution::{
    DifferentialEvolutionConfig, DifferentialEvolutionState,
};

let mut search = DifferentialEvolutionState::new(
    cx,
    &[[0.0, 50.0], [0.0, 50.0]],
    Some(study.initial_coordinates()),
    DifferentialEvolutionConfig {
        population_size: 24,
        differential_weight: 0.8,
        crossover_probability: 0.9,
        seed: 17,
        max_evaluations: 200,
        max_population_entries: 48,
        target: Some(1e-12),
    },
)?;

let report = search.advance(cx, 200, &mut |point| {
    study.evaluate_value(cx, point).map(|r| r.physical.value)
})?;
```

The example above assumes an already constructed `Cx` and scan study.
The generic callback can be any finite, fallible simulation objective.
The optimizer does **not** infer a scan study's budget: explicitly extend
both the optimizer's allowance and the study's physical-evaluation allowance
before continuing. Budget extensions do not restart the search.

`TargetReached` means a consumed observation met the supplied target.
`EvaluationLimit` means the lifetime attempt allowance was spent.
`Paused` means the per-call allowance was spent. None is an optimality
certificate. There is deliberately no small-fitness-spread stop because an
unexplored solid phase may have exactly equal objective values.

## Determinism, bounds and interrupted work

Each trial uses three distinct donors, all different from its target.
Binomial crossover forces at least one donor coordinate, including at CR=0.
Donors always come from the previous complete generation; evaluated children
never influence their siblings. An out-of-box mutant coordinate is redrawn
uniformly within that coordinate's original interval. Equal-valued trials
replace their parent, retaining movement across plateaus; the best-observation
record keeps the earliest exact equal score.

Mutation is computed in unit-box coordinates to avoid overflow from scaled
physical differences. Finite endpoints with a nonfinite interval width refuse.
The declared population-entry cap bounds one population's shape. Two
populations, an optional retained parallel wave, initial point, best point,
and candidate buffers add storage; it is not a total byte or allocator-peak certificate.

Philox draws are keyed by the explicit seed, member ID and generation.
Integer draws have a 64-attempt rejection allowance, and counter segments
are checked before use. Refusal does not silently introduce sampling bias.
The population is published only after a full generation is scored.

Every callback attempt consumes budget, including errors and nonfinite
results. Those attempts retain their exact candidate for explicit retry.
Finite callback results are cached before checking cancellation and are not
recomputed on resume. Callbacks must keep their objective semantics fixed;
their external side effects and panics cannot be rolled back.

Execution is serial by default; the optional tile-backed parallel path is
described below. There is no durable serialization, nonlinear-constraint
adapter, KKT certificate, global optimum,
continuous-time thermal-error bound, or experimental calibration claim.

## Verification scope

The authored native regressions cover donor exclusion, forced crossover,
bounds, every early budget boundary, pause partition equivalence, atomic
generation publication, cancellation, callback retry, multimodal/flat
objectives, cold-start physical fitting, replay, and missed-target reporting.

```bash
cargo test -p fs-dfo differential_evolution
cargo test --manifest-path examples/scan_power_fit/Cargo.toml \
  --example simplex_scan_fit
```

Native compilation and these Rust tests were not executed in the authoring
environment: the required runner/toolchain commands could not start.
Independent Python/NumPy/SciPy checks are supplementary reference calculations,
not evidence that this Rust build or executable passed.

## Parallel physical evaluations

Select an existing `fs-exec::TilePool` explicitly with `--workers`:

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example simplex_scan_fit -- \
  --population --workers 8 --initial-enthalpy 80 --seed 17
```

`--workers` requires `--population` and accepts 1 through 128 requested workers.
Without it the existing serial behavior is unchanged. The example evaluates
waves of at most 24 candidates, never crossing a generation boundary, under a
64 MiB executor lease. That lease is not a total RAM cap on the physical models.
Each candidate builds an isolated trial study using its worker's own `Cx` and
arena. Only immutable specimen, phase-chart and measurement data are shared;
there is no mutex around the physical solver. Each trial is limited to one
physical evaluation, while the original DE state counts the complete search's
actual callback starts. One further serial physical replay is reserved.

The library entry point is:

```rust
let report = state.advance_parallel(
    cx, &pool, &lease, additional_evaluations, max_wave,
    &|worker_cx, point| evaluate_physics(worker_cx, point),
)?;
```

The same Philox candidate keys and ordered greedy selection are used by serial
and parallel execution. Worker completion order does not change a complete
generation. With deterministic callbacks and no target stop, changing wave
width or worker count at the same evaluation prefix gives the same population.
A target is checked **between waves**, so parallel execution may score additional
already-launched candidates after a serial run would first reach its target.
Those evaluations are counted, never refunded. This can change the reported
best candidate and evaluation count without changing DE's generation rule.

Cancellation propagates through the existing executor bridge and all workers
drain before control returns. Successful finite scores are cached before any
cancellation, callback error or executor refusal escapes. Failed or unlaunched
members retain their original coordinates; resumption retries only those gaps.
The serial `advance` API can also resume a previously parallel search. Failed
attempts stay spent. Several callback errors return the lowest failed member's
error; remaining error values are dropped, not interpreted as fitness. No parent
generation is partially replaced because a worker happened to finish early.

The new native regression targets cover serial/parallel populations, hard-budget
boundaries, cancellation with cached successful work, typed failures, serial
continuation, source export/replay and the actual cold-start thermal objective.
They have not been run in this delivery environment: Cargo/rustc, DSR and RCH
are absent. Independent Python scheduling checks are not native concurrency
verification. No wall-time speedup, cross-platform timing or total-memory claim
is made.
