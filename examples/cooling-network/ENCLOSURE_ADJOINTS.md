# Total enclosure adjoints, design and adaptive marking

```bash
frankensim --json cooling-network \
  examples/cooling-network/adjoint-enclosure-pulse.json
```

The actual cooling command now differentiates closed-enclosure radiation in
steady calculations and fixed-grid, fixed-count trajectories. This closes the
previous refusal of enclosure gradients, adjoint-guided fan/workload sizing,
and goal-recovery spatial marking. The existing forward radiosity, conduction,
contact, air transport and energy checks remain the numerical producers.
View factors are fixed declared inputs; their visibility is not computed here.

## Select and read derivatives

For a steady enclosure request, set `objective.gradient=true`. In addition to
the existing inlet, effective-h, contact and total fan-speed sensitivities,
`radiation.adjoint.surfaces` reports `dobjective_dlog_emissivity_k` and
`dobjective_demissivity_k` for each named patch. This is a derivative of the
selected temperature functional, not a derivative of net radiation heat.

For a trajectory, keep `objective.gradient=false` and select the existing
`transient.adjoint` policy:

```json
"adjoint": {
  "qoi": "sampled-peak",
  "max_checkpoint_bytes": 1048576,
  "component_power": true
}
```

Use `qoi: "final"` for the final endpoint. A single cycle reports under
`transient.adjoint`; a repeated trajectory reports under `repeated_cycles.adjoint`
and leaves the last-cycle adjoint null. Its `radiation.surfaces` rows contain
`dtemperature_dlog_emissivity_k` and `dtemperature_demissivity_k`. Emissivity
controls apply throughout the complete trajectory and every cycle. No ambient
reservoir exists in a closed enclosure, so no ambient-temperature derivative
is reported. At emissivity one, only the inward physical perturbation is
admissible; the algebra avoids division by `1-emissivity`.

Existing interval power/fan, original initial-field, common-capacity and inlet
controls use the TOTAL reflected-radiation response. Optional per-component
watts and contact-resistance controls consume the same total nodal multiplier.
An interval control changes every occurrence of that interval. Endpoints after
the selected sampled peak contribute nothing; an earlier cycle's cooldown may
still affect a later peak. Initial-state maxima give zero future controls.

## Differentiate equations, not iterations

The existing primal radiosity producer solves `M J = e sigma T_mean^4`, with
`M = I - diag(1-e) F` and net outward flux `q = (I-F) J`. The new reverse solves
`M^T z = (I-F)^T w`, including a real nonsymmetric transpose for unequal areas
and emissivities. It produces temperature and log-emissivity pullbacks from
`4 e sigma T_mean^3` and `e(sigma T_mean^4 - F J)` respectively. Flux weights
have flux units; patch-watt weights require their area factors.

The transpose uses the existing fs-solver GMRES, explicit linear iteration
limits, normalized objectives and a freshly recomputed original-equation
residual. Closed-enclosure conservation is not enforced by clipping gradients.
The coupled wall-feedback adjoint then combines this response with the actual
air-network transpose and existing IQN acceleration. A fresh unrelaxed equation
residual, not a small accelerated update, governs acceptance.

Because the boundary uses `reference_shifted = reference_air - q/h`, the h
chain rule retains the offset contribution as well as air transport. Radiation
is never treated as air heat. Conductivity K-prime and full contact integration
remain in the solid response. This is not finite differencing each parameter,
iteration-history automatic differentiation, or frozen-radiosity sensitivity.

## Reconstruction and resource limits

A steady gradient retains the exact accepted shifted Robin references and runs
one reconstruction solve with identical inputs. Its temperature bits must match
before the adjoint is used. Steady work reports distinguish forward solves,
reconstruction solves, wall iterations and radiosity Krylov iterations.

A trajectory replays the inner radiation loop at each retained accepted air
reference and immutable previous temperature field. It then linearizes the exact
final shifted rows using `C/dt + J(T_new)`. Both reconstructions must reproduce
the accepted field bits. Every replayed solid callback is counted. The existing
chronological tape propagates `C^T lambda/dt` across all cycle boundaries; no
state or adjoint is reset at a cycle boundary and no extra dt multiplies source
or contact controls.

No matrix history or extra per-frame state is added to the tape. Its existing
checkpoint-memory admission remains in force. Temporary operator/work buffers
are not claimed to fit that tape-only memory allowance. All primal and reverse
work shares the original command deadline. Linear limits apply to each solid
or radiosity solve and derivative limits to each wall-feedback solve; failures
return no partial gradient. No new per-control or per-cycle allowance is added.

## Design and spatial adaptation

The existing steady fan and fixed-trajectory fan/workload searches now consume
these total derivatives. Newton steps only propose guarded points inside an
already evaluated bracket. Every candidate still executes the full physical
model; predicted temperatures never establish feasibility. A derivative request
remains opt-in, and no universal evaluation-count or runtime speedup is claimed.

Steady goal-recovery spatial marking likewise uses the total enclosure adjoint.
Patch names, total areas, emissivities and the declared view-factor matrix remain
fixed through refinement. The global uniform-refinement confirmation remains
mandatory. Marker-only gradients stay null in user-facing sensitivity fields;
marking and reconstruction work is counted. An observed mesh agreement is not
a continuum-temperature or geometry-visibility certificate.

Adaptive-time-grid, automatic time-convergence, periodic-stopping and controller
adjoints remain explicitly unsupported. There are no view-factor or shape
controls, no unique maximum derivative at ties, no continuous-time maximum
certificate, no new native .fsim lowering and no physical-validation claim.

## Focused verification

```bash
cargo test -p fs-cli --test enclosure_radiosity_adjoint
cargo test -p fs-cli --test cooling_enclosure_radiation
cargo test -p fs-cli --test cooling_enclosure_transient
```

Thirteen new regressions comprise three radiosity tests and ten actual-command
tests. They cover nonsymmetric transpose/black limits, complete steady and
trajectory perturbations, nonmatching contact, exact accepted-field and selected
candidate replay, unrolled repeated controls, global mesh confirmation, initial
peaks and budget/adaptive-policy refusals. Existing forward and UQ regressions
remain. These Rust tests were authored but NOT executed here: no Rust toolchain
or executable CI lane was available. Cargo.lock was not regenerated after the
explicit local fs-couple/fs-solver dependency declarations; Cargo resolution,
compilation and actual command behavior remain unverified.

Independent NumPy/SciPy calculations ran 18 radiosity comparisons, 144 complete
steady-physics derivative comparisons and 64 complete-trajectory comparisons.
Their worst absolute discrepancies were 1.78e-9, 3.08e-8 and 6.01e-8 in the
respective control units. Full-FEM and wall-eliminated trajectory transposes
agreed within 1.12e-15; chronological forward tangents and reverse gradients
agreed within 1.51e-14. These are mathematical reference checks, NOT Rust runs.

For the new repeated nonlinear example the reference peak is 315.427903 K at
40 seconds. Its derivatives are -0.751476 and -1.001968 K per unit log-emissivity
change for emitter and receiver, +15.235535 K per heating-power multiplier,
-8.888973 K per common-capacity multiplier, and -0.025698 K per cooldown log-speed
change. The last quantity is nonzero because the FIRST cooldown precedes the
second-cycle peak; the single-cycle cooldown derivative is zero. Deliberately
freezing radiation feedback changes the power derivative by 1.489427 K;
resetting the adjoint at cycle boundaries changes it by 4.411346 K. Material,
fan, storage and visibility declarations in this example are illustrative.
