# Complete material heat transfer in existing thermal observers

`PreparedHeatNetwork::observe(cx, inner, maximum_closed_radiation_defect_j)`
connects source-pinned bulk conduction, thermal-interface resistance and closed
radiation to the existing `DepositionObserver` protocol. This lets the native
three-dimensional thermoplastic solid, and other thermal-observer consumers,
use those same laws without a dummy circuit or a second time integrator.

```rust,ignore
let prepared = heat_network.prepare(cx, &thermal)?;
let observer = prepared.observe(cx, existing_observer, 1e-12)?;
let candidate = solid.advance(
    &thermal, &observer, target_s, &prescribed_displacements_m,
    &nodal_forces_n, mechanical_settings, thermal_settings, coupling_settings, cx,
)?;
let sample = candidate.observer().first.sample();
let internal_roundoff_j = candidate.observer().first.closed_radiation_net_j();
```

The existing `ThermoplasticSolid` supplies displacement compatibility, local
material/yield evaluation, J2 memory, physical plastic heat, reversible exchange
where admitted, and coupled acceptance. The prepared heat network supplies the
same original conductivity, series thermal-contact and radiosity computations
used by its circuit consumer. `CoupledObservers` combines source rates and
refuses two owners for the same conductance. No constitutive or heat solver is
reimplemented here.

## Physical and numerical meaning

This is the existing explicit thermal observer lane. Source powers and link
conductances are evaluated at the START of each accepted thermal substep.
The retained sample is separately refreshed at its accepted END. It is not an
implicit-radiation method and does not inherit the circuit driver's nonlinear
residual or timestep claims. The conductive CFL does not bound nonlinear
radiation stiffness: select `StepConfig::max_dt_s` accordingly and refine it.
The native solid's outer thermoplastic feedback iteration remains unchanged.

Radiative heat remains signed cell power from the original enclosure. It is not
clipped, redistributed, or added a second time. The observer integrates the
source epoch's numerical net power using the actual accepted substep duration,
retains that signed defect, and refuses an advance exceeding the caller's joule
allowance. This counter resets at `begin`, so rejected nonlinear proposals do
not accumulate physical heating. The thermal audit includes generated sources;
subtract the reported net defect once when interpreting external heat. For a
thermoplastic report, the corresponding closed-energy residual is its original
`combined_residual_j + internal_roundoff_j`; the original physical/numerical loss
accounts are otherwise unchanged. Choose the defect allowance consistently with
the consumer's total-energy tolerance. An energy gate is not an accuracy proof.

Electrical resistors, contacts or property bases refuse at adapter creation:
this observer has no electrical-storage or supply-work owner and must not silently
discard one. Use the existing circuit/enthalpy entry points for those assemblies.
The enclosing thermal timeline clones all mutable observer history. Accepted
source samples may be shared through `Arc` because they are immutable; no mutable
physical history is shared. Cancellation and concrete thermal checkpoint errors
survive the original error types.

The prepared topology must remain complete. Apply topology changes through the
existing explicit successor-network operations before resuming this observer;
it does not guess new surfaces, view factors or interface source choices.

## Focused verification

Five regressions are added under `solid_plastic_material_transport`: comparison
with manual original-owner composition in actual loaded 3D tetrahedra; data-only
emissivity/contact changes affecting temperature, free displacement and plastic
history; duplicate-conductance and cancellation rollback; preservation of rich
checkpoint errors and numerical controls; and refusal to drop an electrical
network. Existing fixtures and solvers are reused.

```sh
rch exec -- cargo test -p fs-couple --test solid_plastic_material_transport
```

These Rust tests are authored but unexecuted in the editing environment. There
is no Cargo, rustc, RCH or DSR executable available. Independent Python/SciPy
reference equations passed 27 parameter/time-step cases; they do not execute
Rust lifetimes, source admission, event sequencing, rollback, or cancellation.
The reference fixture is synthetic, not experimentally calibrated. Its parallel
plate enclosure is an algebraic limiting case, not a reconstructed tetrahedral
radiative geometry. No melting, moving contact, finite-strain deformation,
changing radiative geometry, or unconditional stability is claimed.
