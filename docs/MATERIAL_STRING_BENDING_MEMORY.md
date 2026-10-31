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

The `material_strings` performance file exposes the same selection before its
`cavity` record, after any longitudinal selection:

```text
bending-memory 16 1
string-memory 0 1 bending_delta_e bending_tau
```

The first record declares the maximum new memory states and selected-string
count. Each following record names an original zero-based string, a term count,
then that many modulus/time property pairs. Every original layer of that string
supplies its own values from its existing material pack. Unrelated packs and the
body gain no new requirement. Shared packs resolve the union of explicitly
selected properties once. The existing source cache and byte limits remain.

Materials with different spectra can select each physical layer independently:

```text
bending-memory 16 1
layered-string-memory 0 2
layer-memory young_modulus 0
layer-memory young_modulus 2 coating_delta_1 coating_tau_1 coating_delta_2 coating_tau_2
```

Here the core is explicitly elastic, and the coating supplies two relaxation
arms. Layer records follow physical order: core first, then each sleeve from
inside to outside. Each names its equilibrium-modulus property, an arm count
from zero through eight, and that many modulus/time pairs. No dummy loss
properties are demanded from the elastic core. The equilibrium value must match
the already resolved Young modulus so selecting history cannot silently retune
the original elastic operator. Every selected source must supply an explicit
valid frequency band; layer count must exactly match the physical section.

Only a layer's own pack gains its selected requirements; shared packs receive
their union once. Existing `string-memory` records lower through the same shared
owner with uniform property names. Substituting original material data therefore
changes the shared dynamics and observed pressure without changing physics code.
Library callers use `LayerBendingRelaxationProperties` from
`fs_material::section::relaxation`, `LayeredStringBendingRelaxationSelection`, and
`with_material_layered_string_bending_relaxation` for the same binding.

Construction binds body memory first, then string memory, after gas and axial
coordinates have been selected. New string memories start relaxed. All original
schedules, bow interfaces, stops, pressure observation and PCM callbacks remain
in the same performer. An explicit `prepared-linear` selection retains these
histories in the prepared owner, subject to its full state, byte and work budgets.

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
The public coupled owner's `with_prepared_linear_propagation` path retains all
original mechanical and material coordinates in its augmented exponential,
including when preparation precedes material binding. Its explicit augmented
state, byte and multiply-work budgets remain binding. The `material_strings`
performer selects this owner during body construction, before original string
histories are restored, then attaches explicit layered, body and axial histories
in their existing order. Declared dimensions, source-state ceilings and setup
allowances are checked before source loading. Original bending plus torsional
string histories still share the executable's 256-state ceiling. After source
resolution, each positive history also counts in the prepared exponential,
byte and work bounds. The ordinary exact owner's combined
state/setup/energy/Nyquist limits remain binding.
Source validity covers DC through the complete coupled instantaneous stiffness
bound. Each material arm joins the exact coupled exponential, so fast decay
does not impose a midpoint `dt/tau <= 2` restriction. The declared physical
clock and the exponential's representability/work limits remain binding.

This is frozen, small-curvature bending of uniform homogeneous or perfectly
bonded circular layers. It does not model helical winding, interlayer slip,
axial/torsional hereditary loss, finite-strain viscoelasticity or heat evolution.
It does not establish measured violin, cello, bass or other instrument fidelity.

Focused native regressions:

```sh
cargo test -p fs-couple --test modal_material_composition
cargo test -p fs-couple --test material_string_relaxation
cargo test -p fs-couple --example material_strings bending_memory
```

These Rust tests were added but not executed in the editing environment, which
has no Cargo/rustc, DSR or RCH. Independent equation/quadrature checks are not a
substitute for compilation, native renders or experimental validation.
