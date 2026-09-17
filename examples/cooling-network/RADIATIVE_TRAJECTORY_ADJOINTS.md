# Radiative trajectory adjoints

```bash
frankensim --json cooling-network \
  examples/cooling-network/adjoint-radiative-contact-pulse.json
```

Fixed-timestep radiating trajectories now support `transient.adjoint`, including
fixed repeated cycles and the existing adjoint-guided transient sizing callers.
This supersedes the earlier radiative-adjoint refusal described in
`TRANSIENT_RADIATION.md`. The radiation model itself is unchanged: each exposed
patch exchanges mean-temperature gray radiation with a fixed isothermal reservoir.
These example inputs are illustrative, not measured or validated hardware.

Keep `objective.gradient=false` (that flag requests steady gradients). Add:

```json
"adjoint": {"qoi":"sampled-peak", "max_checkpoint_bytes":1048576}
```

inside `transient`. Use `qoi=final` for the final endpoint instead. A single-cycle
result uses `transient.adjoint`; a repeated result uses `repeated_cycles.adjoint`
for the complete chronological history, leaving the last-cycle adjoint null.
The example runs two 150-second pulse/cooldown cycles, 150 accepted endpoints.

## Total feedback, not frozen radiation

The reverse endpoint uses the existing storage response `C/dt + J(T_new)`,
including material K-prime and matching contact conduction. At its accepted
combined Robin rows, the shared radiation transpose includes both the air
reference feedback and the dependence of the radiation coefficient on the
new patch mean. It retains the consistent face integrals on nonuniform patches.
The total nodal-load multiplier is then propagated through the capacity matrix
to the preceding physical field. Earlier cycles are not reset or reordered.

The method differentiates the implicit discrete equations. It does not apply AD
to the Newton, IQN or inner-radiation iteration history, and does not substitute
pointwise T(x)^4 integration for the declared mean-patch closure. Existing steady
radiative adjoints now use this same response-level transpose implementation.

## Controls and outputs

The existing interval power and fan-speed gradients, original initial-temperature
field, common capacity multiplier and inlet-temperature derivatives all include
radiative feedback. Each interval control changes every occurrence of that base
interval in a repeated schedule; it is not an independently redrawn cycle control.
An occurrence after the selected peak has no effect, while earlier occurrences
can matter through thermal history.

The adjoint additionally contains `radiation.surfaces`, with one row per patch:

* `dtemperature_dlog_emissivity_k` changes that patch's emissivity multiplicatively.
* `dtemperature_demissivity_k` is the corresponding absolute-emissivity derivative.
* `dtemperature_dambient_temperature` changes the patch's surroundings in kelvin.

Each patch control is constant throughout the whole trajectory, not a changing
surroundings schedule. Initial-state maxima have zero derivatives with respect
to all later radiation/workload/fan controls. The emissivity domain remains
(0,1]; derivatives at its upper boundary only admit inward perturbations.

Existing transient fan/workload sizing can request this sampled-peak adjoint.
Its safeguarded Newton proposals still require full candidate trajectories;
no predicted peak establishes feasibility, no iteration budget is increased,
and unusable slopes still fall back to bisection. No general speedup is claimed.

## Replay and budgets

The original accepted-field tape and its explicit memory allowance are retained.
No matrix history or extra per-frame radiation state is stored. At each reversed
endpoint, the solver replays the inner radiation boundary loop at the retained
AIR references and immutable previous physical field. It keeps the exact final
combined Robin coefficients and references, then reconstructs their storage
linearization. Both reconstructions must reproduce the accepted temperature bits
before derivatives are used. Mismatch is a refusal, not an approximate replay.

`reconstructed_solid_endpoints` counts endpoints; the added
`reconstruction_solid_solves` counts all replayed inner FEM solves plus endpoint
linearization solves. The enclosing trajectory's total work includes those
actual calls. Forward fields and accepted heat history are not advanced by
reconstruction. The original wall deadline and per-endpoint derivative/linear
budgets still apply; a failed reverse pass publishes no partial gradient.

Adaptive time grids, periodic stopping, and hysteretic-controller derivatives
remain unsupported. Ordinary forward radiating runs retain those existing
capabilities. Mesh studies and effective-h design remain excluded with radiation.
The active sampled maximum is a local selected-branch derivative, not a unique
max derivative at ties, a continuous-time maximum, a discretization certificate,
or a physical compliance guarantee. Contact resistance and material laws remain
fixed controls of this transient adjoint.

## Focused verification

```bash
cargo test -p fs-cli --test cooling_transient_radiation
cargo test -p fs-cli --test cooling_radiation adjoint
cargo test -p fs-cli --bin frankensim radiation
```

Seven new actual-command regressions compare complete perturbed trajectories,
shared repeated versus unrolled controls, unchanged forward fields and energies,
mean-objective ordering, initial maxima, bounded failure and selected workload
replay. A new unit regression rejects changed reconstruction bits. Two existing
radiation-chain-rule unit tests moved into the shared feedback module. The Rust
tests have NOT been executed in the authoring environment, which lacks Rust/RCH.

Independent NumPy P1 calculations passed 120 complete-trajectory finite-difference
comparisons across constant/nonlinear materials, contact/no-contact, single/two
cycles and final/peak objectives. The largest derivative discrepancy was
4.94e-8 in the corresponding control units. Direct and wall-eliminated transpose
solutions differed by at most 5.56e-17. Deliberately omitting radiation's mean
feedback changed a tested gradient by 1.55e-4. These are reference calculations,
not executions of the committed Rust implementation or physical validation.

For the example's two cycles, the reference sampled peak is 307.020212 K at
180 seconds. Its pulse-power-multiplier derivative is 6.999211 K; the first-patch
log-emissivity derivative is -0.241978 K. The cooldown fan's log-speed derivative
is -0.003683 K because an EARLIER cooldown influences that peak. In a single
cycle the peak occurs at 30 seconds and the later cooldown derivative is zero.
