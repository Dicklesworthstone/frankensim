# Closed-enclosure radiation in coupled cooling

```bash
frankensim --json cooling-network \
  examples/cooling-network/enclosure-radiation-gap.json

frankensim --json cooling-network \
  examples/cooling-network/enclosure-radiation-pulse.json
```

Two separate solids exchange radiation across a gap while air cools both.
There is no thermal contact or welded node across that gap. Five watts enter
one body. Radiation can redistribute those watts between bodies, but cannot
remove energy from their combined closed enclosure. Only convection exports
heat to the air network. Materials, emissivities and fan data are illustrative,
not measured or experimentally validated.

The examples explicitly declare the ideal infinite-parallel-plate view-factor
limit. Their finite P1 faces represent equal-area samples of that idealization.
Their finite-geometry visibility, edge leakage and occlusion have NOT been
computed. Admitting a reciprocal matrix does not prove it describes the mesh.

## Declare the radiation model

Choose `radiation.enclosure` instead of the existing `radiation.surfaces` list
of independent fixed-temperature surroundings. Supplying both refuses.

```json
"radiation": {
  "max_iterations": 200,
  "temperature_tolerance_k": 1e-9,
  "relaxation": 0.25,
  "enclosure": {
    "surfaces": [
      {"surface":"emitter","emissivity":0.8,"source":"declared finish"},
      {"surface":"receiver","emissivity":0.6,"source":"declared finish"}
    ],
    "view_factors": [[0,1],[1,0]],
    "row_sum_tolerance": 1e-12,
    "reciprocity_relative_tolerance": 1e-12,
    "evidence": {"kind":"analytic","geometry":"ideal infinite parallel plates"}
  }
}
```

Each named patch must already be a nonempty cooling surface. Areas come from
those mesh traces, not separately entered numbers. This integration requires
convection on the radiating patches; it is not a pure-vacuum solver. Between 2
and 64 patches are admitted. Matrix rows AND columns follow the input surface
order; the complete system is canonicalized by name before use. Coefficients
must lie in [0,1], rows must close to one, and area-weighted reciprocity must
pass the existing fs-conduction admission. Both numerical tolerances must be
explicitly supplied in [0,1e-8]. No factors are clipped or renormalized to make
an invalid enclosure pass. Emissivity lies in (0,1].

`evidence` may alternatively be an externally generated QMC matrix's retained
coordinates: `{"kind":"external-qmc","seed":"73","samples":"100000",
"generator":"named-generator-and-version"}`. Seed and positive sample count
are decimal u64 strings. These labels preserve the caller's provenance; this
command neither runs that sampler nor independently verifies the named formula.
No ambient-temperature field exists in a closed enclosure, and missing external
surfaces cannot be replaced by silently discarding unclosed view-factor rows.

## Numerical and energy coupling

The existing `GrayDiffuseEnclosure` owns the radiosity solve

```
J_i - (1-epsilon_i) sum_j F_ij J_j = epsilon_i sigma T_mean_i^4
q_i = J_i - sum_j F_ij J_j.
```

Each q_i is a UNIFORM outward flux over its original patch. Both reflected
and emitted radiation participate. It is not pointwise integration of T(x)^4.
The adapter adds this flux alongside convection by applying the algebraically
equivalent Robin reference `T_air - q_i/h_air`, leaving h_air unchanged.
That shifted number is an assembly device, not a physical air temperature;
it may be negative and is never passed to the air model or reported as an inlet.
Published solid temperatures must remain finite and strictly positive.

Every iteration solves the real conduction problem, including assigned
conductivity laws and matching/nonmatching contact operators. It then recomputes
radiosity at the solved patch means. The original radiation iteration cap and
relaxation apply; the raw mean-temperature residual and EACH patch's applied
versus recomputed watt mismatch must pass. The radiosity equation residual and
closed-enclosure heat sum are separately checked. The independently accumulated
Robin boundary integral must split into convection and the applied patch flux,
and the full source/air/solid energy gates still apply.

`radiation.surfaces` in the RESULT lists emissivity, mean temperature, radiosity,
irradiation, signed outward heat and applied heat. Negative outward heat is
absorption by that solid, not an error or value to clip. `radiative_out_w` is
the near-zero closed-enclosure balance residual, NOT an external energy sink.
In a steady result, `iterations` refers to the last inner radiation solve and
`total_solid_solves` counts all inner solid evaluations across air coupling.
All iterations share the original command wall deadline. Exhaustion returns
no partial success.

## Transient heating, cooldown and repeated cycles

The pulse example uses two 30-second cycles, each with ten seconds at five watts
and twenty seconds at zero watts. It has nonlinear conductivity and different
constant heat capacities in the two bodies. The deliberately small declared
capacities (20,000 and 10,000 J/(m3 K)) expose storage effects over a short example;
they are not measured material properties. The second cycle starts from the
first cycle's accepted final field, not from a fresh 300 K state.

No additional radiation option is needed: the ordinary `transient` declaration
selects the existing backward-Euler driver. Every inner radiation and air trial
uses the SAME immutable previous physical temperature field. The old surface
means initialize radiation iteration, but emission is recomputed at the NEW
solved means. Freezing radiation at the old temperature is not this model.
The timestep acceptance equation is

```
stored-energy change = dt * (input power - air heat gain - enclosure residual).
```

Signed patch transfers cancel internally. The existing
`radiative_energy_loss_j` field therefore accumulates only the near-zero closure
residual for this CLOSED model, not the emitter's positive transfer or the sum
of absolute patch transfers. It is not a measure of heat exported by radiation.
The final `radiation` report retains both individual signed patch powers and
`temporal_scope: "final-accepted-endpoint"`. Its `iterations` and
`total_solid_solves` are null: final-field heat recomputation cannot recover
those counts. The trajectory's work counters already count every actual solid
callback, including inner radiation iterations and discarded adaptive trials.

Adaptive full steps and rejected half-step pairs never enter accepted storage
or energy history. Fixed-count repetition, ordinary periodic/controller runs,
derivative-free workload/fan sizing, and full-trajectory timestep studies all
reach the same endpoint producer; none substitutes a convection-only trajectory.
Their original timestep, cycle, trial, refinement, solver and wall limits remain
in effect. Complete trajectories retain their usual independent energy checks.

Enclosure steady/trajectory adjoints and goal-recovery spatial marking still
refuse explicitly, including direct attempts to reconstruct a derivative.
There is no frozen-radiosity gradient fallback. Existing reservoir-radiation
adjoints and other non-enclosure workflows are unchanged. General view-factor
generation, occlusion, participating media, changing geometry, continuous-time
peak certification and native .fsim lowering are not implemented by this adapter.

## Uncertain surface finish

```bash
frankensim --json cooling-network-uq \
  examples/cooling-network/enclosure-radiation-gap.json \
  examples/cooling-network/uq-enclosure-finish.json \
  --checkpoint enclosure-finish.uqcp

frankensim --json cooling-network-uq \
  examples/cooling-network/enclosure-radiation-pulse.json \
  examples/cooling-network/uq-enclosure-pulse.json \
  --checkpoint enclosure-pulse.uqcp
```

Each checkpoint destination must be new. The existing `radiation-emissivity`
target resolves against the named enclosure patch and changes the actual input
before each cooling solve. Geometry, view factors, power and other emissivities
stay fixed. For transient UQ, one draw applies throughout the entire trajectory
and its repeated cycles. The observable is the all-cycle sampled peak, not the
cooled final temperature. An interrupted trajectory contributes no partial peak
and retries the same ordinal on resume. A surroundings-temperature target
refuses for an enclosure; nonphysical samples are terminal errors, not clipped
or redrawn. Existing checkpoint and confidence policies remain in use. The
example probability distributions are illustrative, not empirical evidence.

## Focused checks and boundaries

```bash
cargo test -p fs-cli --test cooling_enclosure_radiation
cargo test -p fs-cli --test cooling_enclosure_transient
cargo test -p fs-cli --bin frankensim uq_command::model::radiation
```

The original eight enclosure regressions remain. Eight additional command tests
cover manufactured time endpoints with reversed heat flow, nonlinear repeated
storage, unrolled and matrix-axis replay, adaptive accepted-history replay,
whole-trajectory time refinement, actual workload candidates, exhausted budgets,
unsupported derivatives, and byte-identical transient-UQ checkpoint continuation.
These Rust tests have NOT been executed in the authoring environment, which
lacks Rust and network access for installing it. Compilation and actual CLI
behavior remain unverified.

Independent Python P1 calculations use direct endpoint equations with analytic
two-plate radiation and air elimination, and a separate nested radiosity/solid
iteration. Twelve transient endpoint comparisons had at most 3.65e-11 K field
difference. Two manufactured endpoints had at most 1.71e-13 K field error;
freezing radiation at the old field instead caused 0.09831 and 0.60137 K errors.

For the two-cycle nonlinear pulse example, the independent sampled peak is
315.427903 K at 40 seconds. Its 100 J input splits into 35.006204 J stored and
64.993796 J exported by air. About 15.594296 J is redistributed internally from
emitter to receiver, with zero net radiative export in the analytic reference.
Freezing old-temperature radiation changes the peak by 0.10510 K. Refining from
30 to 60 and 120 total steps gives reference peaks 315.427903, 315.564990 and
315.634980 K, illustrating why a successful algebraic solve is not a timestep
error bound. These are independent numerical references, NOT executions of
FrankenSim, runtime speedups, continuum certificates or physical validation.
