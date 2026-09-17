# Closed-enclosure radiation in coupled cooling

```bash
frankensim --json cooling-network \
  examples/cooling-network/enclosure-radiation-gap.json
```

Two separate solids exchange radiation across a gap while air cools both.
There is no thermal contact or welded node across that gap. Five watts enter
one body. Radiation can redistribute those watts between bodies, but cannot
remove energy from their combined closed enclosure. Only convection exports
heat to the air network. Materials, emissivities and fan data are illustrative,
not measured or experimentally validated.

The example explicitly declares the ideal infinite-parallel-plate view-factor
limit. Its finite P1 faces represent equal-area samples of that idealization.
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
those mesh traces, not separately entered numbers. This first integration
requires convection on the radiating patches; it is not a pure-vacuum solver.
Between 2 and 64 patches are admitted. Matrix rows AND columns follow the
input surface order; the complete system is canonicalized by name before use.
Coefficients must lie in [0,1], rows must close to one, and area-weighted
reciprocity must pass the existing fs-conduction admission. Both numerical
tolerances must be explicitly supplied in [0,1e-8]. No factors are clipped or
renormalized to make an invalid enclosure pass. Emissivity lies in (0,1].

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
`iterations` refers to the last inner radiation solve; `total_solid_solves`
counts all inner solid evaluations across air coupling. All iterations share
the original command wall deadline. Exhaustion returns no partial success.

Steady nonlinear materials, existing contact models, derivative-free fan sizing,
and uniform mesh studies use the same producer. The existing reservoir mode is
unchanged. Enclosure adjoints, goal-recovery marking and transient trajectories
currently refuse explicitly; none returns a frozen-radiosity derivative or
silently falls back to the old reservoir model. General view-factor generation,
occlusion, participating media, changing geometry, and native .fsim lowering
are not implemented by this adapter.

## Uncertain surface finish

```bash
frankensim --json cooling-network-uq \
  examples/cooling-network/enclosure-radiation-gap.json \
  examples/cooling-network/uq-enclosure-finish.json \
  --checkpoint enclosure-finish.uqcp
```

The checkpoint destination must be new. The existing `radiation-emissivity`
target resolves against the named enclosure patch and changes the actual input
before each cooling solve. Geometry, view factors, power and other emissivities
stay fixed. A reservoir-temperature target refuses for an enclosure. Nonphysical
samples are terminal model failures, not clipped or redrawn. Existing checkpoint,
exact-ordinal retry and candidate/confidence policies remain in use. The example
probability distribution is illustrative, not empirical uncertainty evidence.

## Focused checks and boundaries

```bash
cargo test -p fs-cli --test cooling_enclosure_radiation
cargo test -p fs-cli --bin frankensim uq_command::model::radiation
```

Eight new Rust regressions include six actual-command tests and two target
mutation tests. They cover independently eliminated two-plate/mixed-air balances,
black and near-mirror limits, nonlinear material influence, matrix-axis replay,
uniform refined-field replay, physical and iteration-budget refusals, actual
uncertain samples and byte-identical checkpoint/result continuation. These tests
have NOT been executed in the authoring environment, which lacks Rust. Compilation
and actual CLI behavior remain unverified.

Independent Python calculations assembled the P1 volume and boundary operators,
solved the coupled equations directly, and compared against a separate nested
mathematical implementation plus analytic two-plate heat balances. Twelve full
thermal comparisons and nine formula/null-mode cases passed. The largest
whole-field discrepancy between direct and nested methods was 2.70e-10 K.

For the example, the independent reference peak is 338.496251 K. The emitter
and receiver means are 337.883111 K and 305.958702 K, with 1.263461 W crossing
the gap and five watts exported by convection in total. Black finishes increase
internal exchange to 1.809428 W; lowering emitter emissivity to 0.2 reduces it
to 0.569927 W. A deliberately DIFFERENT model replacing mutual exchange with
two 300 K reservoirs exports 1.860575 W by radiation and is not equivalent.
These are independent numerical references, not FrankenSim executions, accuracy
certificates, measured performance or physical validation.
