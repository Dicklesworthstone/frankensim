# Closed transient material/electrical/thermal feedback

Enable `fs-lbm`'s opt-in `implicit-enthalpy-coupling` feature. The existing
explicit `enthalpy-coupling` entry points and their stepping are unchanged.

`step_joule_enthalpy_implicit` lives at
`fs_lbm::conjugate::enthalpy::transport::interfaces::implicit`. It advances one
stationary control-volume interval with the original material charts, current
solver, ordered electrical contacts and thermal links. Temperature, phase,
bulk resistivity, contact resistance and thermal conductance are evaluated at
one converged endpoint, rather than a source cached at the previous time.

## Physical residual and shared ownership

The unknown vector is **cell energy increments in joules**, in sorted thermal
`CellId` order. That order need not be electrical voxel order. For each cell:

```text
h_new = h_old + energy_increment / reference_mass
residual_j = reference_mass * (h_new - h_old) - dt * endpoint_net_power
```

`fs-thermal::enthalpy::backward_euler::BackwardEulerStep` exposes this residual
for other nonlinear consumers. It uses the same heat-rate routine as the
existing explicit timeline: no second contact, reservoir or convection law.
Phase comes from the original enthalpy curve; a constant Kelvin temperature
does not imply zero storage or constant resistance during melting.

The LBM adapter runs the existing `fs-solver::NewtonKrylovState`. Each residual
probe invokes the existing material/interface/current evaluation owner. Trial
Joule powers are not accumulated as physical heat. Only the accepted endpoint
power is integrated over the actual interval. Source and interface payloads
are retained without a Clone requirement. Source-backed interfaces can use
`ElectricalInterfaceResponse::try_from_source` from the existing public API.

Independent additional powers may be supplied by `external_sources`; they use
the same endpoint state. Do not duplicate a heat channel already represented
by an authored boundary or the electrical solver. Stateful observers that debit
mechanical or material histories during trials are not admitted by this API.

## Controls, refusal and event semantics

Supply an initial energy-increment guess and explicit `ImplicitJouleConfig`.
Its numerical Jacobian uses directional finite differences, with
`difference_step_j` setting the largest coordinate perturbation in joules.
These are **not analytic tangents or validated physical sensitivities**.
Choose an admitted stencil and check finite-difference refinement separately.
An unsupported stencil or material query refuses; there is no extrapolation,
source substitution, silent line-search domain repair or explicit fallback.

`input.max_current_solves` bounds every material/current epoch, including
Jacobian probes and the fresh final check. The active thermal unknown count,
Newton attempts, inner cycles and restart length have separate controls.
Even loose Newton stopping cannot bypass the independent maximum per-cell and
total energy gates. Late refusal, exhausted work or cancellation returns no
new thermal state. Original input objects remain unchanged.

Process events due at the initial time before preparing a step. Split at the
next event: an event strictly inside an interval refuses. A successful step
returns the **left limit** at its target; target-time events remain pending.
After every coupled owner accepts the interval, process those events explicitly
with the existing zero-duration `advance_to(target_s, ...)` operation. Starting
another implicit interval before doing so refuses. Drive discontinuities also
require interval boundaries; do not average an on/off pulse implicitly.

## Focused verification

From a complete checkout with pinned dependencies materialized:

```sh
cargo test -p fs-thermal
cargo test -p fs-lbm --features implicit-enthalpy-coupling --no-fail-fast \
  --test joule_enthalpy --test joule_enthalpy_transport \
  --test joule_enthalpy_composed --test joule_enthalpy_interfaces \
  --test joule_enthalpy_implicit
```

Six thermal regressions cover a heterogeneous two-body backward-Euler solution,
latent storage, cancelling local defects, nonlinear endpoint conductance,
event ordering and cancellation. Seven actual electrical integration regressions
cover analytical resistive feedback and first-order time refinement, stiff
reservoir cooling, data-only latent changes, current/contact/transport/source
consistency, rollback, independent energy acceptance and zero-current floating terminals.

These are numerical verification cases, not measurements of real specimens.
Material-source support and empirical validation remain separate obligations.

The model is fixed geometry, fixed reference mass, stationary Ohmic transport
and control-volume backward Euler. It does not simulate melt motion, magnetic
induction, mechanical work, phase kinetics or changing topology. Algebraic
convergence does not establish nonlinear uniqueness, timestep accuracy,
unconditional stability for arbitrary laws, or measured material validation.
