# Source-pinned relaxation-rate temperature changes

`FrameRelaxationBinding::set_layer_temperatures` connects prescribed layer
conditions to the existing complete structural/fluid/material trajectory.
Each layer retains its original material card, selected property claims,
coordinate conventions, geometry, and every internal strain. Changing a layer's
temperature queries its original relaxation-time curves and changes that same
body's restoring-force history, including reactions onto connected bodies.
There is no second material evaluator, integrator, or accepted clock.

## Physical meaning

A native Maxwell arm stores `(sqrt(k) * b.q - z)^2 / 2`. Its relaxation time
controls the irreversible rate, not that storage. Updating only time at an
accepted boundary therefore leaves displacement, velocity, fluid momentum,
material internal strain, instantaneous stress, and mechanical energy unchanged.
It changes subsequent motion and instantaneous dissipation. Initial internal
history is never inferred again from endpoint displacement or reset to rest.

The original source must declare density, equilibrium Young's modulus and every
modulus increment as validity-wide scalar constants. Equality at two queried
temperatures is NOT enough; even a numerically flat response curve is refused
by this rate-only adapter. Temperature-dependent storage needs its own potential,
entropy, reference-state and thermal-work treatment. This path does not silently
change modulus, mass, geometry, expansion, or memory topology.

Temperature is prescribed here, not solved. A caller may supply independently
known temperatures, but this operation does not advance a caloric state, deposit
integrated relaxation heat, or close a self-heating feedback loop. The fixed
reference-geometry, perfectly bonded scalar-longitudinal law and exterior-flow
limits of `MATERIAL_FRAME_RELAXATION.md` still apply. In particular this is not
an arbitrary material swap, phase transition, or model for thermal expansion.

## Source-owned API

```rust,ignore
use fs_couple::air_path::cylinder::frame::modal::relaxation::temperature::{
    FrameRelaxationLayerTemperature,
};

// Retain the original FrameRelaxationBinding and native mechanical owner.
// Supply ALL bound layers in member declaration order, then radial layer order.
let rates = material.set_layer_temperatures(
    &mut mechanical,
    &[
        FrameRelaxationLayerTemperature {
            element: 0, layer: 0, card: &original_core_card,
            temperature_axis: "T", temperature_k: core_kelvin,
        },
        FrameRelaxationLayerTemperature {
            element: 0, layer: 1, card: &original_coating_card,
            temperature_axis: "T", temperature_k: coating_kelvin,
        },
    ],
    &gate,
)?;
let frame = mechanical.step_under_gate(&complete_native_forces, &gate)?;
```

The source evaluator pins exactly the previous claims even when another claim
might become applicable. The complete original query is preserved except for
one existing absolute-temperature coordinate. Kelvin and the source domain are
checked. Source snapshots in `material.members()` refresh even if the numerical
rates remain unchanged; the immutable source model remains reference data.
The returned rate report contains only genuine rate changes and their before/
after instantaneous powers, NOT interval heat or thermal parameter work.

All layer queries finish before native preparation. The native rate update
prepares the ordinary, prepared, or harmonic-only operator under its original
budgets before publishing. A missing/foreign/out-of-domain source, changed
storage law, stale native branch, failed solver preparation, or cancellation
leaves both accepted owners unchanged. No fallback property or partial layer
refresh is published. Foreign material branches and every fluid history remain.

The underlying reusable API is
`CoupledModalSystem::set_material_relaxation_times` with
`coefficients::rates::MaterialRelaxationTime`. It has ordered unique branch
indices and rejects edits to source-protected tensor histories. It is useful to
other constitutive producers but does not itself confer source/thermal validity.
Do not use it behind a `FrameRelaxationBinding`: that would intentionally make
the binding stale, and its next source refresh/observation must refuse.

Unchanged numerical times avoid propagator preparation. Changed times rebuild
the original native transition/loss operators; this is not an allocation-free
update or real-time qualification. Layer projections, solid spectra, fluid
quadrature, state layout, and initial conditions are not reconstructed.
Hold selected conditions through each sample. Split imposed temperature events
at sample boundaries; refine piecewise-held smooth histories independently of
structural and fluid spectral resolution. Small energy residuals are not time
accuracy bounds, source-fit evidence, or experimental validation.

## Focused qualification

```sh
rch exec -- cargo test -p fs-couple --lib render::schedule::force::coupled::coefficients::rates
rch exec -- cargo test -p fs-couple --lib air_path::cylinder::frame::modal::relaxation::temperature
```

Three native-rate and three source-temperature regression groups are authored.
They cover analytic constrained relaxation, warm continuation, prepared/ordinary/
harmonic parity, separate core/coating conditions, retained source identities,
foreign histories, late-source refusal, cancellation/retry, unchanged-rate
receipt refresh, and rejection of temperature-dependent storage.

Native Rust compilation, runtime tests, and rustfmt remain UNEXECUTED in the
editing executor: Rust/Cargo/rustfmt/RCH/DSR are unavailable, and network toolchain
hosts cannot be resolved. Source/API inspection and independent equations are
not substitutes for native qualification.

Independent NumPy/SciPy calculations ran 12 piecewise-rate cases and a constrained
internal-history case, then 12 two-layer wet/dry temperature cases. For the
layered cases, separately assembled full-momentum and normalized-coordinate
trajectories differed by at most 1.23e-16 (normalized), and DOP853 reference
states differed by at most 3.62e-13. The largest cumulative mechanical-energy
residual was 9.78e-18 J. All temperature-boundary storage jumps were zero.
In the synthetic one-coordinate wet beam, changing only the prescribed layer
conditions from 300/300 K to 320/330 K changed the final transverse displacement
from 0.984374 mm to 0.980855 mm. These are manufactured equations and source
curves, not native Rust execution or validation of any real material.
