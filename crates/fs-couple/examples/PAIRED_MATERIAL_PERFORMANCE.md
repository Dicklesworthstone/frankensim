# Independent material blades → physical pressure → streamed audio

The library path is
`bernoulli_aperture::paired_network::performance::PairedAperturePerformance`.
It plays the existing `PairedApertureNetwork`, not a second synthesized voice.
The scalar/symmetric `AperturePerformance` API remains unchanged.

## Assemble the physical data before playback

Reduce each original plate separately with `PlateApertureReduction`. Each may
have its own density, stiffness, regional sections, shape and material receipts.
Pass the two single-blade reductions to `PairedAperture::from_plates` with one
shared slit/contact law and the actual inlet gas/impedance. No identical-motion
constraint is applied. Both reductions must describe the same total rest gap
and slit width; these are not doubled.

Attach each desired `PlateRelaxationSpec` with
`with_blade_relaxation(blade_index, law, initial_history)`. Its regions partition
THAT blade's original triangles. Initial viscous displacements refer to its own
opening-positive displacement, never half the combined gap. Leave a blade
unattached to keep it elastic. Separate viscous damping must be zero for a
blade whose loss is supplied through material memory.

Equilibrium modulus and Poisson ratio must match the original section matrices.
Changing equilibrium material data requires rebuilding that source plate, not
changing only its additional Maxwell arms. A different relaxation spectrum can
change frequency-dependent stiffness, memory loss and sound while equilibrium
mass, stiffness and pressure area remain identical. Incompatible material bands,
unresolved time scales or excessive combined storage refuse rather than retune.

Bind that source to `PairedApertureNetwork::new` or `with_viscothermal`, using the
actual per-section `GasState` values and existing bore/load geometry.

## Schedule blowing and unequal loads

Use the ordinary `GestureSchedule` with a single `BlowingPressure` track, and the
existing `AperturePerformanceConfig`. Its sample rate must exactly match the
physical clock; its finite window must fit the source's admitted horizon.

```rust,ignore
let mut voice = PairedAperturePerformance::new(
    network,
    ApertureObservation::NetworkBaffled { node: outlet, receiver },
    blowing_schedule,
    config,
)?
.with_plate_forces(0, first_blade_force_program)?
.with_plate_forces(1, second_blade_force_program)?;

// Values are pressure in pascals, not normalized or fitted audio samples.
voice.block(&mut pressure_block)?;
```

The original `with_blade_forces(Vec<BladeForceControl>)` API remains available
for already-generalized coordinate forces, with its existing strict ordering.
Use either that representation OR physical mesh programs, not both.

Each `ApertureForceProgram` addresses nodes or triangle patches on its own blade.
Input is total transverse force in newtons on that blade. The original plate
owner derives the work-conjugate force/velocity projection from its actual mode.
Events are mechanical-sample addressed and held; same-sample assignments retain
input order, last assignment per port winning. Overlapping ports add their
physical loads. Force and pressure changes participate in the same coupled
sample. Releasing a force does not reset either blade or its material history.

`max_controls` covers pressure assignments plus BOTH force programs.
`max_compile_work` independently bounds pressure compilation and the combined
force-projection/grouping work. Neither blade receives an extra hidden budget.

## Reuse the existing output path

The player implements `PressureRenderer`: pass it directly to
`render_pressure_pcm16`, or wrap it in `DecimatedRenderer` before feeding
`Pcm16WavStream`. Use an explicit pressure-to-PCM full-scale value, shared across
material comparisons; do not normalize each result separately. The existing
streaming cancellation retains the accepted source prefix, both material
histories, force/pressure cursors, wave state and receiver delay. Callback
partitioning does not change the mechanical clock.

`Inlet` and `NetworkNode` are internal pressures. `NetworkBaffled` is the existing
one-way retarded observation of actual terminal volume flow, using that outlet's
radius and gas. It assumes an infinite rigid baffle and does not add radiation
back-loading. Retain or explicitly choose the physical acoustic terminal load.
A physical/receiver callback failure poisons the renderer; do not replay an
already accepted mechanical sample into a failed observer.

This is a supplied one-mode-per-blade, uniform-slit material model. It does not
claim measured oboe/bassoon accuracy, spatial rolling contact, automatic wet-cane
identification, live material/gas state transfer or a new paired text-file CLI.
