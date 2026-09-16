# Discrete transient cooling adjoints

The `cooling-network` command can differentiate the selected final temperature
or sampled trajectory peak through its actual backward-Euler solid and mixed-air
endpoint equations. This is separate from the steady gradient fields.

```bash
frankensim --json cooling-network examples/cooling-network/adjoint-contact-pulse.json
```

Keep `objective.gradient` false: that field requests a steady adjoint. Add an
explicit policy inside the existing `transient` object instead:

```json
"adjoint": {
  "qoi": "sampled-peak",
  "max_checkpoint_bytes": 1048576
}
```

`qoi` is either `final` or `sampled-peak`. The former uses the selected spatial
objective at the final endpoint. The latter uses its largest value over the
initial state and every accepted endpoint, selecting the earliest exact time
tie. A spatial maximum keeps the existing lowest-ID active-vertex rule. The
result is the derivative of that selected branch, not a unique smooth maximum
derivative at a tie or a continuous-time peak bound. Area-mean objectives also
work, with named regions bound in network order rather than declaration order.

## Returned controls

`transient.adjoint` reports the selected value, time, state index and active
vertex, together with:

- `dtemperature_dinitial_temperatures`: one partial derivative per initial solid
  node. `dtemperature_duniform_initial_k` sums them for a uniform initial shift.
- `dtemperature_dcapacity_multiplier_k`: derivative at multiplier one when ALL
  declared element heat capacities are scaled together across the trajectory.
- `dtemperature_dinlet_temperatures`: partial derivatives per hydraulic node for
  supply-temperature changes maintained throughout the schedule; nonsupplies
  remain zero.
- Per-interval `dtemperature_dpower_multiplier_k` and
  `dtemperature_dlog_fan_speed_ratio_k`: the first scales that interval's ENTIRE
  currently declared source at multiplier one. It is not an absolute-watt or
  per-component derivative, and a zero-source interval gives zero. The second
  changes only that interval's fan speed and includes the existing single-bank
  affinity, air-capacity and supported convection response. It is null without
  a fan or when the required convection slope is unavailable.

Intervals strictly after a selected earlier peak have zero influence on that
branch. An initial-state maximum needs no reverse endpoint solves. None of the
steady `fan_speed_sensitivity`, inlet or contact-gradient fields is repurposed
as a transient result; those fields remain null for transient requests.

## Numerical implementation and resources

The endpoint equation is `C(T_new-T_old)/dt + A(T_new)T_new - b = 0`. Its solid
Jacobian is `C/dt + J(T_new)`, including the existing material K-prime and
matching-contact terms. The shared Robin-response algebra and mixed-network
implicit adjoint close the air feedback before the total nodal-load multiplier
is carried backwards as `(C/dt)^T lambda`. The old state is never advanced by
an internal coupling or Newton iteration.

The forward run is unchanged. Only accepted temperature and reference vectors
are retained. After the complete trajectory and its energy checks pass, the
reverse traversal reconstructs one endpoint at a time at its retained old
state, timestep and references. It requires exact same-profile temperature
replay before using that operator. It does not rerun a whole trajectory per
control, retain every solver matrix, or differentiate IQN/Newton iterations.

`max_checkpoint_bytes` is required and bounds retained field/reference and
frame storage, up to 512 MiB; it is not a total solver-workspace or allocator
memory guarantee. These checkpoints are in-memory derivative state, not durable
resume files. Admission precedes the first PDE solve. The original wall deadline
covers forward and reverse work, and the existing derivative/linear budgets
apply at each reverse endpoint. Failures or cancellation publish no partial
adjoint. Reconstruction solves are included in `total_solid_solves`; the result
separately reports `forward_solid_solves`, `reconstructed_solid_endpoints`,
`adjoint_sweeps` and the maximum interface-equation residual. Forward nonlinear
statistics continue to describe forward work, not reconstruction work.

This first path requires a fixed, nonrepeated schedule without nested fan/power
sizing. Combining `adjoint` with `adaptive`, `repeat`, `power_design` or
`fan_speed_design` refuses explicitly. Geometry, material laws, contact
resistance and fluid properties remain fixed. Material derivative kinks or
sampled validity endpoints refuse; no arbitrary unique slope is invented.
There is no phase change, time-grid derivative, continuous peak certificate,
physical validation or uncertainty certification claim.

## Focused checks

```bash
cargo test -p fs-conduction --test backward_euler_adjoint
cargo test -p fs-cli --test cooling_transient_adjoint
cargo test -p fs-cli --test cooling_transient_mean_adjoint
```

The ten added tests cover insulated heating, nonlinear endpoint controls,
transpose identities, history propagation, actual-binary trajectory finite
differences, exact forward-field retention, future-control causality,
initial-state maxima, region ordering, cancellation and budget refusals.
The authoring environment has no Rust toolchain: these test sources were not
executed there. Independent dense NumPy FEM checks compare 48 trajectory
control derivatives and direct versus Schur-complement adjoints; those are
numerical references, not Rust execution or runtime benchmarks.
