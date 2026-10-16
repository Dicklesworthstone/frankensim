# Physical force programs on retained plate valves

A plate valve can now be driven by independent **mechanical forces and mouth
pressure at the same time**. Pressing, releasing or pulling at a declared point
or patch changes its actual opening, distributed lay contact, material-memory
state, swept flow and returning duct pressure. No displacement, closing pressure,
rest geometry, damping or PCM gain is rewritten to imitate that actuation.

```bash
cargo run --release -p fs-couple --bin music_render -- \
  wind crates/fs-couple/examples/plate-valve-forced.performance \
  /tmp/plate-valve-forced.wav --decimate --block 37
```

The example retains the existing illustrative plate, asymmetric slit, lay,
Maxwell memory, 96 kHz clock, 1.5 kHz compact radiation boundary and receiver at
0.2 m. It adds a tip-node force and a uniform-traction patch near the tip. The
pressure phrase remains unchanged, and force release preserves ringdown. Forces
are authored Newton values, not measured lip forces. The new 1 Pa PCM full scale
is an explicit encoding choice; neither simulation pressure nor forces are
normalized. The same complete file works as `ensemble --valve` input.

## Source records

An optional block after `relaxation` and its records, before `compile_limits`:

```text
force_ports 2
force_port node 13
force_port patch 2 10 11
force_events 7
force_event 0 0 -0.0002
force_event 73 0 -0.001
force_event 512 1 -0.0008
force_event 1024 0 0
force_event 1200 1 0
force_event 4000 0 -0.0005
force_event 6000 0 0
```

`node` names one original plate node. `patch COUNT TRIANGLE...` names unique
original triangles. Its force is the **total force in N**, distributed uniformly
over their physical area, not force per triangle or a pressure in Pa. Positive
force points along positive transverse plate displacement; negative presses in
the opposite direction. Up to 32 fixed footprints are admitted. A fixed-support
site can have zero motion in the retained basis: it then does zero modal work.

Each event is `force_event SAMPLE PORT FORCE_N`. Samples are integer positions
on the mechanical clock, NOT output-WAV or pressure-control ticks. Apply before
the named sample; hold until that port's next event. Every port initially holds
zero. Events may be unsorted; simultaneous events keep input order and the last
assignment to each port wins. Different ports add in port order, even when their
footprints overlap. No intermediate same-sample assignment advances physics.
An event exactly at a callback end belongs to the next callback. Events at or
past the end of the source window refuse rather than being silently discarded.

Omitting the block preserves the pressure-only path and its arithmetic.
Nonfinite forces, absent sites, duplicate patch triangles, unsupported records,
and overflowing projected force sums refuse. The existing `max_controls` in
`compile_limits` bounds pressure assignments plus raw mechanical assignments;
the force compiler also checks its projection/summation work against
`max_compile_work` independently of the existing pressure compiler's bound.
Force magnitudes are not clamped to force convergence. Slope/contact/material,
load and receiver limits remain in force; a physically inadmissible candidate
fails rather than moving the specimen into a different model.

## Projection and work

`PlateApertureReduction::force_port` derives a signed coefficient B from the
actual retained plate shape. For a node it is that node's transverse motion per
metre of opening. For a patch it is the area-weighted mean of the existing P1
transverse displacement trace. The same B maps BOTH sides of the port:

```text
Q_opening = B F_physical
v_actuator = B v_opening
Q_opening v_opening = F_physical v_actuator
```

The equilibrium sections, supports and original physical mode determine B.
There is no guessed lever arm, pressure-area substitution or new excitation
oscillator. `DynamicAperture::step_with_force`, `ApertureTube::step_with_force`
and `ApertureNetwork::step_with_force` accept generalized opening force directly.
It enters the existing implicit momentum/pressure/contact/material equation.
`ApertureFrame::mechanical_work_j` is F times midpoint velocity times dt, kept
separate from fluid supplies and counted once by the complete balance methods.

For finite library playback, attach `ApertureForceProgram` with
`AperturePerformance::with_plate_forces`. Footprints are compiled against that
exact owned specimen, not a separately supplied reduction. Attach once before
playback; the object retains pressure and force clocks, all physical state and
receiver histories. Applied/pending force events, force ports, held generalized
force and the last accepted mechanical work remain inspectable. Compilation is
cold; scheduling allocates nothing during callbacks. Physical solver allocations
are unchanged, so this is not an allocation-free or real-time claim.

Only a successful complete sample consumes controls. Request-shape refusal and
cancellation do not consume future pressure or force events. A failed physical
callback poisons the finite pressure renderer, consistent with its existing
contract; the underlying per-sample valve/network still does not publish its
failed trial. Corrected direct solver inputs can be retried through that lower
API without losing accepted vibration or wave/material memory.

## Scope

These are prescribed point loads or uniform-traction patches on one linear
structural mode. They do not add a deformable lip/finger body, a moving lay,
frictional tissue contact, force/position servo, changing footprint or multimode
bending. They can supply or withdraw physical work and are not assumed passive.
Distributed slit/contact and all existing gas, material and radiation model
limits remain unchanged. Stepwise force programs can excite frequencies outside
a declared load/receiver band; no full-transient accuracy or measured-instrument
calibration follows from finite-band checks. Existing `plate_aperture` and
`music_render_wind` targets cover the changed mechanics and complete consumers.
