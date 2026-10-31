# Resolved thermal and molecular fields in continuing bore playback

`fs_couple::cavity_modal::interior::thermoviscous::dg::ThermoviscousDg`
connects the original three-dimensional thermoviscous finite-element time
trajectory to a continuing `DgPortLoad`. It retains pressure, Cartesian
velocity, translational temperature, every original molecular temperature,
finite conductive-wall temperatures and structural displacement/velocity.
The bore keeps its own fields, branches, material walls, finite passages and
terminal radiation histories.

This supplies a causal composition path for spatial acoustic regions such as
small cavities, windows or junctions. It does not turn those regions into
resolved nonlinear jets, solve breath-throughflow, or qualify an instrument.
The current flute and whistle performance descriptions do not automatically
select this model; they still need suitable geometry and qualified source
connections.

## Physical connection

The field and bore have one original gas state and one exact physical clock.
The selected FE port must have unit pressure shape on an exclusive planar
surface, with the same area as the bore inlet. The mode's conjugate is actual
inward volume flow. Other modes keep their original authored tractions.
An additional area discontinuity, chimney or end correction is not inferred.
Matching area does not prove shape equivalence or justify neglecting higher
interface modes; interface placement and spatial refinement remain necessary.

With the field's mean relation `Q_FE = Q_free + beta*p` and the bore's original
relation `p = h + Z*Q_DG`, continuity gives:

```text
Q_DG = -Q_FE
p = (h - Z*Q_free) / (1 + Z*beta)
```

The original prepared pHS trajectory supplies the complete affine mean map.
There are no auxiliary probe trajectories, fitted response, added integrator,
or per-step matrix exponential. Pressure is held on the field; its actual
interval-mean flow is held on the bore. This partition requires coupling-clock
refinement independently of either spatial mesh or material model.

Both original candidates must agree on pressure and swept volume. Their
interface work cancels inside the combined account; authored external pressure
and structural force remain external supply. Each frame retains original FE
heat channels, the DG loss receipt, terminal work, complete storage and balance
residuals. An inferred DG loss remains an inferred receipt; it is not relabeled
as independently measured material heating. Rational terminal radiation retains
its original signed-work accounting.

No participant publishes until all physical and numerical checks and the final
cancellation poll succeed. Rejected steps preserve both accepted histories.

## Original material preparation

Construct `ThermoviscousInterior` from its original gas or molecular source.
Use `prepare_time` for the supported fluid/structural fields, or
`prepare_time_with_wall_diffusion` when actual finite solid walls are attached.
The latter takes an explicit element count through each original wall thickness.
It retains the source conductivity, heat capacity, backing and thickness;
it does not invent a finite wall for a half-space.

Pass that preparation, the same source, an unadvanced `DgPortLoad`, the selected
FE port index and explicit `ThermoviscousDgLimits` to `ThermoviscousDg::new`.
The constructor checks both the source FE identity and the actual prepared
calorics, wall fields and source band. Preparing the same FE without its
declared molecular reservoirs is refused. `with_initial_state` admits complete
nonzero field/wall states before the first step, including original molecular
and solid temperature perturbations.

Each direct `step` takes complete held pressure and structural-force vectors.
The connected pressure slot must be zero in the authored vector: its reaction
is solved by the coupler. Callers inspect `state`, `load` and `last_frame`
without receiving mutable access to either accepted owner.

## Pressure, microphone and WAV output

`dg::performance::ThermoviscousDgPerformance` implements the existing
`PressureRenderer` interface. Its finite `ThermoviscousDgControl` history has
one complete physical pressure/force vector per original interval. Explicit
configuration bounds the output rate, callback size and total stored controls.
Changing callback partitions does not change the control clock.

`ApertureObservation::Inlet` emits the solved interval-mean interface pressure.
`TubeBaffled` uses the existing computed DG exterior and every original
terminal flow at the declared microphone. `with_observer_gas` attaches the
shared exterior-material propagation law, including its phase and causal delay.
It neither adds a second reaction nor replaces the interior molecular dynamics
with an output filter.

The shared pressure decimator, ensemble and PCM streaming functions can consume
this renderer without retiming the physical system. A failed callback poisons
playback, preserving completed samples and the last accepted physical state.
A microphone failure after an accepted step does not undo that physical step.

## Runnable composition

From the repository root, with the pinned constellation and toolchain available:

```sh
cargo run --release -p fs-couple --example thermoviscous_bore -- output.wav
```

The example uses explicit named terrestrial-air composition and the existing
composition-limited molecular law in both the resolved region and the DG bore.
A 15 × 2 × 4 mm manufactured region with no-slip, isothermal walls receives
an authored 2 ms, one-pascal pressure pulse. A 50 mm circular bore of equal
interface area has its original viscothermal and molecular memory, computed
baffled radiation, and a microphone 0.1 m from the outlet. The shared writer
exports 20 ms at 48 kHz with a fixed one-pascal PCM scale and refuses to replace
an existing output file.

This small example exercises the public composition. It is not a real mouth,
jet, flute, tone-hole model or listening reference. Its coarse spatial mesh,
uniform interface projection, time partition and declared source/observer bands
require refinement. The input pulse is not asserted to be band-limited.

Five new native regression groups cover original-owner replay, molecular and
finite-wall evolution, actual downstream feedback, energy and rollback,
binding admission and callback-independent playback. Rust syntax was parsed
with the pinned formatter during implementation. Compilation, native numerical
tests and example execution were not run because this executor lacks the full
pinned dependency checkout and remote build runners. No audible improvement or
instrument milestone is claimed by these changes.
