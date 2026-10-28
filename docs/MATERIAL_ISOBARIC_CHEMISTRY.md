# Constant-pressure reacting gas

`GasReaction::isobaric_implicit_candidate` advances a closed homogeneous gas
against a quasistatic load holding its original pressure. It requires the
existing `with_implicit_mechanism` policy and returns an immutable candidate.
Its `heat_j` is external whole-interval heat, not reaction enthalpy.

The SAME coupled extent solver now selects either the original prescribed
volume/work constraint or a constant-pressure enthalpy constraint. For every
residual and Jacobian query, the latter reconstructs

```
n1 = n0 + S xi
H1 = H0 + Q
p1 = p0
V1 = m / rho(T1, p0, n1)
xi_r = dt V1 rate_r(T1, p0, n1)
```

All channels see one material endpoint. The endpoint volume, not the initial
volume, multiplies every rate. No outer pressure iteration, additional Newton
unknown for volume, independent chemistry engine, constant-Cp approximation,
or second chemical heat source is introduced. The existing native Newton,
finite-difference Jacobian, source-backed rate scales, detailed balance,
reachability, inventory checks and complete final kinetic reports remain in
charge. Trial kinetic evaluations retain the existing receipt-free path.

The underlying material operation is
`ClosedGasParcel::replace_model_with_heat_at_constant_pressure`. It conserves
mass while retaining species formation references in total enthalpy. It admits
the new model's continuous caloric branch, uses the original safeguarded root
solver, and resolves full transport only at the accepted endpoint. It needs no
fictitious intermediate with new composition at old temperature. Identical
models use the reference-safe frozen-composition heating path. Differing
models must share a thermochemical reference convention; the chemistry owner
additionally enforces exact original species/source identities.

## Independent physical acceptance

The returned `IsobaricGasReactionStep` retains the original chemical boundary
report plus its independent enthalpy change and residual. It checks

```
Delta H - Q
Delta U + p0 (V1 - V0) - Q
Delta S - Q / T1
```

The first two use the existing absolute energy tolerance. The last retains the
original endpoint heat-entropy inequality, including backward-Euler numerical
entropy; it is not isolated chemical production. Heat reservoirs need their
own entropy accounting. No coefficient of performance or passivity theorem
is inferred. Mass remains exact; every species and channel residual is checked
against its explicit amount tolerance before any candidate is returned.

## Scope and refusal

This is first-order, homogeneous, closed ideal-gas chemistry at constant
pressure, not an open-flow reactor, flame/shock solver or inertial pressure
controller. Source-backed material/kinetic laws supply fidelity. Existing
fixed-volume scalar/subcycling selectors refuse on this entry rather than
silently changing its physical constraint. Repeated caller-selected intervals
provide temporal refinement; no arbitrary-step stability is claimed.

Composition-changing balances retain absolute formation references. This is
not a claim that all chemical differences remain resolvable under arbitrary
reference magnitudes. Unrepresentable balances, source/transport failures,
exhausted nonlinear work and cancellation return no partial physical state.

## Focused checks

`cargo test -p fs-couple --test reactive_isobaric` covers a stiff two-channel
chain, exothermic heat/work, a two-to-one molecule-count change, first-order
refinement, fixed-volume reduction, canonical ordering, budgets and cancellation.
Five material regression groups cover the lower-level composition closure.
Native compilation and tests were unavailable in the authoring environment.
Independent Python reference equations checked 512 material endpoints and 240
reactor endpoints, an analytical stiff chain, and first-order refinement
against an analytical isomer solution and DOP853. These checks do not execute
the Rust or qualify experimental chemistry.

Physical convention: Cantera 3.2, Constant Pressure Reactor, energy equation:
https://www.cantera.org/3.2/reference/reactors/constant-pressure-reactor.html
