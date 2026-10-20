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
populations, optional initial point, best point, and candidate buffers add
storage; it is not a total byte or allocator-peak certificate.

Philox draws are keyed by the explicit seed, member ID and generation.
Integer draws have a 64-attempt rejection allowance, and counter segments
are checked before use. Refusal does not silently introduce sampling bias.
The population is published only after a full generation is scored.

Every callback attempt consumes budget, including errors and nonfinite
results. Those attempts retain their exact candidate for explicit retry.
Finite callback results are cached before checking cancellation and are not
recomputed on resume. Callbacks must keep their objective semantics fixed;
their external side effects and panics cannot be rolled back.

Execution is serial. There is no parallel population evaluator, durable
serialization, nonlinear-constraint adapter, KKT certificate, global optimum,
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
