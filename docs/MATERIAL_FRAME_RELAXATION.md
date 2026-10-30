# Longitudinal material memory in complete frame/fluid assemblies

`air_path::cylinder::frame::modal::relaxation` connects the original layered
material spectra to the existing complete mechanical owner. It works on the
same `MaterialFluidFrame` with or without exterior-fluid regions; it creates
neither another structural solver nor another material-history integrator.

Each selected member supplies one scalar longitudinal Young-modulus spectrum
per original concentric layer. That SAME E(t) acts on axial extension and both
bending coordinates. Area and second moment come from the existing section;
member length/orientation and reciprocal force/moment maps come from
`Frame2d::reference_strain_traces`. Native Kelvin-Voigt assembly now uses those
same geometric maps. No geometry is inferred from a material name.

The original `with_prony_bending` source resolver admits the coefficient bundles
and retains every source selector, receipt, equilibrium modulus and relaxation
time. The frame adapter uses those scalar longitudinal spectra and the original
area/moment integrals for its three strain channels. This is an explicit
longitudinal-law interpretation, not permission to treat an anisotropic
bending-only fit as an axial law. Source coefficients must be in the original
resolved section bundle and support DC through the actual sampled Nyquist band.
The source Young modulus must equal the explicit equilibrium modulus.

```rust,ignore
use fs_couple::air_path::cylinder::frame::modal::relaxation::{
    FrameLongitudinalRelaxation, FrameRelaxationBudget,
};
use fs_material::section::relaxation::LayerRelaxationProperties;

// Build/assemble the original complete frozen frame component first.
let (placement, mechanical) = source.coupled_component(
    &initial_state, sample_rate, modal_budget, coupling_budget, &gate,
)?.into_system(prepared_budget, &gate)?;
let selected = [FrameLongitudinalRelaxation {
    element: 0,
    layers: vec![LayerRelaxationProperties {
        equilibrium_modulus: "long_time_modulus".into(),
        branches: vec![("modulus_increment_1".into(), "relaxation_time_1".into())],
    }],
}];
let (material, mut mechanical) = placement.with_longitudinal_relaxation(
    mechanical, &selected, None,
    FrameRelaxationBudget {
        max_members: 16, max_source_arms: 64,
        max_memory_coordinates: 192, max_setup_terms: 1_000_000,
    }, &gate,
)?;
let report = mechanical.step_under_gate(&complete_native_force, &gate)?;
let physical = placement.physical_state(mechanical.components(), &gate)?;
let internal = material.states(&mechanical, &gate)?;
```

Each positive layer arm supplies three native Maxwell strain coordinates, with
stiffnesses delta_E*A/L, 3*delta_E*I/L^3 and delta_E*I/L^3. The equilibrium
structural springs remain exactly once. Each member's three projected columns
are computed once and reused across all its layers and arms. The full loaded
inertia congruence is applied; no per-mode damping ratio or frequency retuning
replaces the material dynamics. Native prepared/exact solvers own the complete
augmented state, instantaneous-frequency admission and work/storage/loss.

`None` explicitly means initially relaxed arms at the supplied displacement;
`Some(z)` supplies the energy-normalized histories, ordered by member, layer,
source arm, then axial/symmetric-bending/antisymmetric-bending channel. Only
exact zero modulus increments have no dynamic coordinate; their source records
remain visible. Positive arms with zero mechanical projection are retained:
they can hold and dissipate initial constrained material energy. Foreign-frame
memories remain an unchanged prefix when another placement appends its laws.

Nonzero Kelvin-Voigt material loss on a selected member is refused to avoid an
unspecified second loss mechanism. This adapter requires cold construction,
zero frozen preload, fixed reference geometry/density and fixed source
conditions. Declare a frame's selected material members together. Existing
participating relaxation, duplicate members, missing properties, altered
equilibrium springs, unsupported source bands and exceeded budgets refuse.

Transient AND harmonic queries must use the returned `CoupledModalSystem`.
The original `MaterialFluidFrame::prepare_harmonic` and source pHS remain the
immutable equilibrium/fluid image and do not acquire these additional histories.
The native harmonic convention is exp(-i*omega*t), unlike the standalone
frame harmonic adapter's exp(+i*omega*t); retain the convention at observations.
Material dissipation and stored material energy remain separate from exterior
fluid loss. `states()` reports instantaneous power, not interval heat. No heat
is automatically deposited into solid cells, and live E(t,T)/thermal-history
updates, anisotropic laws, interlayer slip, plasticity and phase change are not
implemented by this binding. Small residuals are not experimental validation.

## Focused checks

```sh
rch exec -- cargo test -p fs-solid --lib frame2d::viscosity::strain_tests
rch exec -- cargo test -p fs-couple --lib air_path::cylinder::frame::modal::relaxation
```

Three geometry and five material-coupling regression groups are authored.
They compare original canonical pHS dynamics/complex response, independent
integrated curvature, constrained initial histories, source-only substitutions,
foreign-memory retention, budgets and cancellation. Rust/Cargo/rustfmt/DSR/RCH
were absent in the editing executor; native compilation and execution are
UNVERIFIED. Independent Python/SciPy equations ran nine geometric cases, seven
single-body histories, 28 frequency points and three connected/disconnected
two-body cases. Maximum normalized state discrepancy was 6.08e-16; maximum
combined-energy discrepancy was 1.86e-17 J. These are synthetic reference
equations, not execution of the Rust source or validation of measured alloys.
