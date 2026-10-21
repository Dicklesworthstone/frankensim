# Material bending memory in reciprocal string banks

`ReciprocalStringBank::with_material_string_bending_relaxation` binds material
relaxation to the existing coupled string/body/contact state. It is a physical
constitutive selection, not a fitted per-mode damping curve or an audio effect.

For each selected original string, `StringBendingRelaxationSelection` names
one to eight pairs of properties from every layer's original resolved material
bundle. Each pair is an **additional bending modulus** in Pa and its positive
relaxation time in seconds. The original Young modulus is the equilibrium
modulus. All coefficients must be constant at the selected material state;
frequency-sampled coefficients are not silently converted into a broadband law.
An unselected string receives no new mechanism and needs no spectral properties.

Example library selection, after resolving those properties with the ordinary
material-state resolver and constructing the bank:

```rust,ignore
use fs_couple::bowed_string::runtime::reciprocal::bank::relaxation::
    StringBendingRelaxationSelection;

let bank = bank.with_material_string_bending_relaxation(
    vec![StringBendingRelaxationSelection {
        string: 0,
        properties: vec![("bending_delta_e".into(), "bending_tau".into())],
    }],
    None, // New arms start relaxed at the current displacement.
    64,   // Explicit maximum number of new material memory states.
    &gate,
)?;
```

This library API consumes source data; the existing `material_strings` file
syntax does not yet expose this selection. No instrument preset is changed.

## One energy and dynamics owner

For a mass-normalized sine mode n of a uniform span, a material layer contributes
`c_n^2 = DeltaE * I_layer * (n*pi/L)^4 / mu`. Each layer, positive arm, transverse
plane and retained mode gets an independent internal strain. The original
section owner supplies `I_layer`; the actual current transverse mass supplies
`mu`, including selected cylinder added inertia. Tension is not relaxed. The
same original `CoupledModalSystem` evolves all these states alongside bodies,
bridges and other material memories. Its Maxwell storage, restoring force and
loss quadrature feed the existing bow solve and pressure observation.

The already selected Kelvin-Voigt bending viscosity is an explicit **parallel
instantaneous dashpot**. Spectral source data must describe additional material
memory, not the same measured loss counted again. Explicit zero increments
retain their source evidence and add no states or decay. Negative increments,
missing/wrong-quantity coefficients and incompatible source domains refuse.

## Construction and physical limits

Select body memory and cylinder gas loading before string memory. The initial
history vector, when supplied, addresses only the new string channels, ordered
by selection, layer, positive arm, transverse plane and sine mode. A zero vector
on displaced material is unrelaxed. Existing body history is not reinitialized.
The source getters retain original layer sections, spectra and projected rows.

Later supported torque/topology rebuilding re-embeds the same history in the
same string basis. Appending longitudinal coordinates preserves it through the
existing common extension. A changed string mass or basis refuses stale rows;
it needs a fresh cold binding. Playback cannot be used as a material-swap point.
The prepared mechanical-only backend still refuses constitutive memory. The
ordinary exact owner's combined state/setup/energy/Nyquist limits remain binding.
Source validity covers DC through the complete coupled instantaneous stiffness
bound, and the declared clock must satisfy dt/tau <= 2 for every positive arm.

This is frozen, small-curvature bending of uniform homogeneous or perfectly
bonded circular layers. It does not model helical winding, interlayer slip,
axial/torsional hereditary loss, finite-strain viscoelasticity or heat evolution.
It does not establish measured violin, cello, bass or other instrument fidelity.

Focused native regressions:

```sh
cargo test -p fs-couple --test modal_material_composition
cargo test -p fs-couple --test material_string_relaxation
```

These Rust tests were added but not executed in the editing environment, which
has no Cargo/rustc, DSR or RCH. Independent equation/quadrature checks are not a
substitute for compilation, native renders or experimental validation.
