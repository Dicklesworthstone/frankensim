# Rare-event cooling probability learned from physical pilots

`cooling-network-uq BASE UQ --importance-pilot ROUNDS,SAMPLES,DEFENSIVE,ELITE,RELAXATION,MAX_SHIFT_NORM`
learns a Gaussian proposal from complete cooling solves and then estimates a
fixed failure event using separate, independently keyed production solves.
It calls the same cooling parser, parameter lowering, solver, completed-trajectory
selector and child watchdog as ordinary cooling UQ. There is no replacement
thermal model, polynomial prediction or clipped physical input.

```sh
cargo run -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-importance-pilot-hotspot.json \
  --importance-pilot 2,8,0.2,0.25,1,8
```

This small example spends **32 total physical calls**: two eight-call learning
pilots, then sixteen independent production calls. Its count and temperature
threshold illustrate the interface, not a resolved rare-event probability or a
qualified cooling device. The original UQ `samples` count includes every phase.

## Explicit sampling policy

The six comma-separated choices are fixed before any physical work:

- `ROUNDS`: 1..16 complete pilot rounds; no result-dependent round count.
- `SAMPLES`: at least two calls per pilot, leaving at least two production calls.
- `DEFENSIVE`: original-target mixture fraction in (0,1], with finite reciprocal.
- `ELITE`: upper-response fraction in (0,1] used by the existing weighted fit.
- `RELAXATION`: convex weight in (0,1] on each learned shift.
- `MAX_SHIFT_NORM`: Euclidean limit in (0,1e6] in independent standard-normal
  coordinates, not physical units.

The initial shift is zero. `fs_uq::PilotImportanceExecution` owns all draws,
likelihood weights, weighted proposal updates, phase keys and reductions.
The input must declare `temperature_limit_k` and explicit independent or joint
Gaussian inputs. The original numerical factorization policy, including admitted
singular laws and constant marginals, is unchanged. Uniforms and unknown dependence
refuse; Gaussian physical-domain violations fail the experiment without redrawing.
The existing command limit of 10,000 total physical calls applies to all phases.

Learning uses likelihood-weighted high-response samples before sufficient failures
are seen, then targets strict failures when enough are observed. Ties, relaxation
and norm limiting retain the library's semantics. A poor proposal can increase
variance; opposite or disconnected failure modes can defeat a single mean shift.
Defensive mixture mass preserves ideal-law support, not finite-budget discovery.

## What the result means

The strict event is **QoI > temperature_limit_k**. Equality is not a failure.
Only final production responses enter `estimated_failure_probability` and its
`failure_probability_standard_error`. Pilots never enter those statistics or
select production's allocation after seeing production results.

The `frankensim.cooling-network-uq.importance-pilot.v1` result retains the physical
QoI, input declarations, all phase seeds, final proposal, fitting policy, pilot
counts/cutoffs/effective sample counts, and separate pilot/production call counts.
`production_proposal_failures` is a count under the tilted proposal, NOT a
failure count sampled from the original target law. `mean_likelihood_weight`
and `failure_effective_samples` are diagnostics, not accuracy guarantees.
The weighted estimate is neither self-normalized nor clipped to [0,1]. No raw
proposal quantiles, target temperature distribution, complement compliance
probability or pooled-pilot probability is manufactured.

Standard errors are descriptive. Zero observed failures and zero empirical
standard error do not prove zero risk. There is no confidence interval,
statistical stopping rule, physical approval or achieved-accuracy claim.
Finite-generator, arithmetic, solver and physical-model errors are not bounded.
`Complete` and exit SUCCESS mean the fixed experiment completed successfully.

For explicit `qoi.kind=transient-sampled-peak`, one callback is the initial/
accepted-endpoint maximum of one **complete trajectory**, including fixed
repeated cycles. A partial trajectory never contributes a failure indicator.
These are numerical sampled peaks, not continuous-time peak certificates.

## Budgets and failures

One soft `wall_seconds` allowance covers setup, all pilots, production, fitting
and report preparation. The existing watchdog kills/reaps interrupted cooling
children; bounded synchronous fitting/reductions are not preemptible. Timeout
returns BUDGET and no partial probability. A real solver failure, nonfinite
result, numerical refusal or invalid physical draw returns a terminal diagnostic;
no failed observation becomes a zero, replacement sample or ordinary-MC fallback.

This CLI mode has no durable combined checkpoint. Recovery/chunk flags, other
sampling methods, mean controls, sensitivity, sequential compliance and candidate
selection refuse before input reads or output reservation. The underlying pilot
executor retains its existing in-memory continuation, but this command makes no
disk-recovery promise. Existing MC/QMC/Sobol/Shapley/PCE command paths are unchanged.

Focused native checks:

```sh
cargo test -p fs-cli --bin frankensim uq_command::execute::importance
cargo test -p fs-cli --test cooling_importance_pilot
```

The tests reconstruct a full study with actual cooling calls, cover singular
correlation and complete transient observations, and check counts, deterministic
output, invalid combinations, deadlines and failed physics. They must run natively
before claiming a Rust execution pass or any measured accuracy or speedup.

The fixed-proposal `--importance EPSILON,SHIFT...` mode remains available and
unchanged; see `IMPORTANCE_UQ.md`. The two proposal selectors are mutually
exclusive. Both use the same physical evaluator and watchdog.
