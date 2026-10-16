# Component-power and contact controls for cooling trajectories

```bash
frankensim --json cooling-network \
  examples/cooling-network/adjoint-component-contact-pulse.json
```

The example has chip, memory and overlapping standby-chip footprints, two
nonlinear materials, a resistive bondline, radiation and mixed airflow. Each
14-second cycle has six seconds of heating and eight seconds of cooldown.
Two fixed cycles share thermal history. Inputs are illustrative declarations,
not measured hardware or validated material data.

Inside the existing `transient.adjoint` request, select the new controls:

```json
"adjoint": {
  "qoi": "sampled-peak",
  "max_checkpoint_bytes": 1048576,
  "component_power": true,
  "contact_resistance": true
}
```

Both flags default to false. `qoi: "final"` selects the final endpoint instead.
Ordinary adjoints retain their existing controls without doing these additional
source/contact contractions. The new sections are within `transient.adjoint`
for a single cycle and `repeated_cycles.adjoint` for a complete repeated run.
They do not fill in the steady `contact_sensitivities` field on a transient.

## Independent component watts, including dormant components

`component_power_sensitivities.intervals[i].rows` names every component in the
original `solid.component_power` map. Each row reports:

- `applied_power_w`: that component's actual watts in this interval.
- `dtemperature_dpower_w_k_per_w`: the derivative with respect to ADDING watts
  to this component in every occurrence of this base interval.
- `dtemperature_dpower_multiplier_k`: multiplying only this component's actual
  interval watts; equal to the absolute derivative times `applied_power_w`.

To replay an absolute interval perturbation, express that interval using
`component_powers_w`, naming every original component and preserving the other
watts. This also works when the original interval used `power_scale`. It does
not renormalize other components to hold total power fixed. At zero watts, the
multiplier derivative is zero but the absolute derivative can be nonzero; the
admissible physical perturbation there is one-sided.

`base_rows` reports `dtemperature_dbase_power_w_k_per_w`. This changes the
original component's `solid.component_power` watts, together with its declared
system total. Only `power_scale` intervals depend on those base watts; absolute
`component_powers_w` intervals override them. Thus the base derivative sums
`power_scale * interval_absolute_derivative` over those intervals only.

The source calculation is the transpose of the actual consistent P1 volume
integration, not a lumped approximation or a ratio of current powers. The new
`StepLinearization::source_density_pullback` computes `M_source^T lambda` once
per reverse endpoint. The adapter sums those density weights over each original
footprint and divides by its retained PowerMap bound volume. Overlapping
footprints remain independent. No per-component physics solve is needed.
For a temperature objective, the intermediate density gradient has units
K/(W/m3); the normalized component gradient has units K/W.

## Persistent contact resistance

`contact_resistance_sensitivities.rows` reports each named interface's
`dtemperature_dlog_resistance_k` and `dtemperature_dresistance_w_m2`, alongside
its actual `resistance_m2_k_w`. A control changes that one resistance throughout
the complete trajectory and all cycles, with fixed geometry and overlap topology.

For `p=ln(R'')`, the discrete contact operator obeys `dK/dp=-K`. Each endpoint
therefore contributes `lambda^T K_contact T`, using the total coupled nodal-load
adjoint. Matching P1 faces reuse their full face-mass integration; nonmatching
planar traces reuse their admitted common-refinement quadrature. A product of
mean jumps is not substituted. Contact orientation changes signed heat reporting,
not the physical resistance derivative. This extends the earlier steady-only
resistance-control boundary described in the nonmatching-contact documentation.

The existing endpoint Jacobian already includes `C/dt`. Neither the source nor
contact contraction gets another dt factor. Propagating `C^T lambda/dt` backwards
carries earlier-cycle effects into later peaks. Only endpoints through the
selected peak contribute; changing a shared interval can still affect the peak
through an EARLIER occurrence even when its final occurrence is in the future.

## Same forward solve and bounded extra work

The existing primal, reconstruction, material, air/radiation, residual and energy
checks remain mandatory. Accepted temperature bits must still replay before any
endpoint derivative is used. These controls do not add a perturbed trajectory,
linear solve or coupled adjoint per component/contact, but their integrations
and reductions do consume CPU time under the original cancellation/deadline gate.

The optional accumulator arrays and vector headers are charged together with the
accepted-field tape against `max_checkpoint_bytes` before forward timesteps.
The component-by-interval grid is capped at 65,536 rows. This is not a total-process
memory cap: solver scratch and ordinary result buffers remain separate workspace.
Missing PowerMap footprints or contact declarations refuse when their control is
requested. Nonfinite derivatives and cancelled/rejected reverse work return no
partial adjoint. An initially selected maximum reports zero future-source/contact
sensitivities without inventing missing reverse solves.

Fixed timesteps, fixed cycle counts, nonlinear conductivity, matching/nonmatching
contacts and the existing surface-mean radiation model are supported. Existing
sizing candidates can retain these requested diagnostics. This does not add a
multi-component optimizer, transient contact search, geometry derivatives,
adaptive/controller derivatives, a unique maximum derivative at ties, a
continuous-time peak certificate, or hardware validation.

## Focused checks

```bash
cargo test -p fs-conduction --test backward_euler_source_controls
cargo test -p fs-cli --test cooling_transient_controls
```

Eight new Rust regressions comprise three numerical-kernel tests and five
actual-command tests. They cover consistent source integration, prescribed-row
elimination, a dormant component, complete trajectory finite differences,
repeated/unrolled control sums, split contact controls, overridden base powers,
unchanged forward fields and reverse solve counts, and memory/input refusals.
These Rust tests have NOT been executed in the authoring environment, which has
no Rust toolchain. Compilation and actual CLI behavior remain unverified.

Independent NumPy P1 calculations performed 84 symmetric complete-trajectory
comparisons and 14 admissible one-sided/base-override comparisons. Maximum errors
were 2.34e-9 and 1.40e-7 respectively, in each corresponding derivative's units.
A separately propagated forward tangent agreed with the direct transpose within
1.25e-16. Deliberately using lumped loads, adding an extra dt, or resetting history
at cycle boundaries produced clearly different component gradients.

For the example, the independent sampled peak is 302.286669 K at 20 seconds.
Heating-interval absolute sensitivities are 0.150445 K/W for chip, 0.006017 K/W
for memory, and 0.071983 K/W for the zero-watt standby footprint. The persistent
contact derivative is 0.035966 K per unit log-resistance change. The all-occurrence
cooldown chip sensitivity is nonzero because the FIRST cooldown precedes that
peak; its multiplier derivative is zero because its actual watts are zero.
These are independent mathematical references, NOT execution of FrankenSim,
measured speedups, continuum-error bounds or physical validation.
