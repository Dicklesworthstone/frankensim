# Resolved thermal and molecular fields in continuing bore playback

`fs_couple::cavity_modal::interior::thermoviscous::dg::ThermoviscousDg`
connects the original three-dimensional thermoviscous finite-element time
trajectory to one or more continuing `DgPortLoad` owners. It retains pressure, Cartesian
velocity, translational temperature, every original molecular temperature,
finite conductive-wall temperatures and structural displacement/velocity.
Every bore keeps its own fields, branches, material walls, finite passages and
terminal radiation histories.

This supplies a causal composition path for spatial acoustic regions such as
small cavities, windows or junctions. It does not turn those regions into
resolved nonlinear jets, solve breath-throughflow, or qualify an instrument.
The current flute and whistle performance descriptions do not automatically
select this model; they still need suitable geometry and qualified source
connections.

## Physical connection

The field and all attached ducts have one original gas state and one exact
physical clock. Each selected FE port must have unit pressure shape on an
exclusive planar surface, with the same area as its duct inlet. The mode's
conjugate is actual inward volume flow. Unconnected modes keep their original
authored tractions.
An additional area discontinuity, chimney or end correction is not inferred.
Matching area does not prove shape equivalence or justify neglecting higher
interface modes; interface placement and spatial refinement remain necessary.

With the field's mean relation `Q_FE = Q_free + beta*p` and the bore's original
relation `p = h + Z*Q_DG`, continuity gives:

```text
Q_DG = -Q_FE
p = (h - Z*Q_free) / (1 + Z*beta)
```

For several ducts, the same complete field supplies the matrix relation
`Q_FE = Q_free + B*p`. With one original inlet impedance per duct, the
simultaneous closure is:

```text
(I + diag(Z)*B) * p = h - diag(Z)*Q_free
Q_DG = -(Q_free + B*p)
```

Every off-diagonal entry of `B` is retained. A finite cork passage or airway
therefore reacts on the main bore through the shared spatial field. The small
constant interface matrix is factored once during cold attachment using the
existing dense LU owner; each interval solves only its changing right-hand side.
No symmetry of the matrix or independence of its ports is assumed.

The original prepared pHS trajectory supplies the complete affine mean map.
There are no auxiliary probe trajectories, fitted response, added integrator,
per-step factorization or per-step matrix exponential. Pressure is held on the
field; its actual interval-mean flows are held on the original ducts. This partition requires coupling-clock
refinement independently of either spatial mesh or material model.

Every original candidate must agree on pressure and swept volume at its
interface. Each interface's work cancellation is checked independently before
the combined account is admitted; authored external pressure
and structural force remain external supply. Each frame retains original FE
heat channels, the DG loss receipt, terminal work, complete storage and balance
residuals. An inferred DG loss remains an inferred receipt; it is not relabeled
as independently measured material heating. Rational terminal radiation retains
its original signed-work accounting.

No participant publishes until all physical and numerical checks and the final
cancellation poll succeed. Rejected steps preserve the field and every duct's
accepted history, including when a later duct or the combined work check refuses.

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

Before the first step, `with_additional_duct(port, load, gate)` attaches another
unadvanced original duct at a distinct eligible FE port. This does not change
the primary `port()` or `load()`. `coupled_ports()` returns all connected FE
identities, primary first; `load_at_port()` and `additional_loads()` expose the
other accepted owners without mutable access. A frame's `additional_ducts`
contains each original output, common pressure, mean flow and interface work
residual. The frame's storage, loss, terminal-work and balance totals include
the complete field and all attached ducts.

Each direct `step` takes complete held pressure and structural-force vectors.
Every connected pressure slot must be zero in the authored vector: its reaction
is solved by the coupler. External excitation belongs in unconnected pressure
slots or original structural forces.

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

`ThermoviscousDgPerformance::new` observes the primary duct.
`new_at_port(system, port, ...)` selects another attached FE-to-duct connection
for inlet pressure or that duct's own computed-baffle microphone. All ducts
continue advancing in either case; selecting an observation does not disconnect
the others or sum unrelated terminal signals.

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
composition-limited molecular law in the resolved region and both DG ducts.
A 15 × 2 × 4 mm manufactured region connects a 20 mm closed cork-side passage
at its left face (FE port 0) and a 50 mm radiating bore at its right face
(FE port 1). Each circular duct has the same interface area as its rectangular
FE face, its own 12-pole viscothermal memory and its original molecular
reservoirs, on the same 48 kHz clock.

An independent opening on the region's top face receives an authored 2 ms,
one-pascal triangular pressure pulse at FE port 2. Controls are `[0, 0, pulse]`;
the two internal reactions come from the coupled solve. The remaining faces
are no-slip and isothermal. The existing computed baffled radiation and
microphone 0.1 m from the main outlet remain attached to the primary bore;
the closed passage's returning wave influences that microphone through the
original field and bore dynamics. The shared writer
exports 20 ms at 48 kHz with a fixed one-pascal PCM scale and refuses to replace
an existing output file. Console output includes the aggregate balance and
the closed passage's actual final interface pressure, flow and work residual.

This small example exercises the public composition. It is not a real mouth,
jet, flute, tone-hole model or listening reference. Its coarse spatial mesh,
uniform interface projection, time partition and declared source/observer bands
require refinement. The input pulse is not asserted to be band-limited.

## Resolved nonlinear flow as the source

`dg::load::ThermoviscousDgLoad` exposes an unconnected planar FE opening as a
flow-driven acoustic load. It eliminates the attached duct pressures through
the complete FE cross-port matrix and keeps every original field and duct
history. Its affine `p = H + Z*Q` response is the exact response of the admitted
mean-port partition. The effective mobility includes every internal return;
it is not just the external diagonal of the bare FE response. Other unconnected
pressures and structural forces are held at zero by this load.

`nsf_dg::thermoviscous::NsfThermoviscousDg` couples the original nonlinear
Navier–Stokes–Fourier volume to that external opening. The flow zone and FE
region must occupy distinct physical volumes. The existing characteristic
pressure solve, source material, fluid trial and interface quadrature remain
the owners of the source. The final original NSF trial, complete FE state and
all finite duct states publish together after the pressure, projection-work,
internal-work and cancellation checks. Returning waves from either duct can
therefore affect the original nonlinear source through the shared field.

Whole characteristic faces and explicit source/downstream aperture patches are
supported. Drives are original boundary pressure/characteristic controls,
original transverse reservoir states, or an explicitly attached finite material
source. `with_finite_reservoir(aperture, volume_m3, limits, gate)` binds one
existing quiescent transverse opening to the same named-mixture or NASA9
material owner. Its accepted inventory supplies the opening's pressure and
temperature; `FiniteReservoir` and `FiniteReservoirActuation` are then its
exclusive drive forms. Prescribed boundary/reservoir controls cannot overwrite
that owner.

Actual streamed mass and conservative total energy update the retained chamber.
`FiniteReservoirActuation` adds explicit volume motion and heat through the
original simultaneous midpoint law after that exchange. Positive volume rate
expands the chamber; positive heat enters its gas. This held-source, flow-first
partition requires interval refinement. It keeps no chamber momentum: incoming
resolved kinetic energy thermalizes in the original well-mixed inventory.
Pressure work already included in the NSF total-energy flux is not debited
again as acoustic work.

The finite chamber's storage, mass/energy residuals, streamed exchange and
original midpoint heat/work receipt remain in `frame.finite_reservoir`.
`finite_reservoir_inventory()` exposes the accepted named or NASA9 source
without manufacturing a NASA9 parcel for a named-mixture law. Its candidate
publishes together with the original NSF trial, complete field and all ducts.
A rejected caloric reconstruction, inventory bound, acoustic check or cancelled
trial preserves every accepted owner. A finite upstream acoustic duct remains
a distinct attachment; `UpstreamPressureRise` is still refused on this path.

This is still an acoustic projection about zero reference flow. The FE region
does not receive mean mass, species or advected enthalpy. The complete NSF mass
and total-energy accounts, streamed-mass projection difference and signed
pressure/flow covariance remain separate observable quantities. Covariance is
never supplied as invented heat to close the acoustic account.

`nsf_dg::thermoviscous::performance::NsfThermoviscousDgPerformance` implements
the existing `PressureRenderer`. It observes the actual selected attached
duct's mean inlet pressure or computed microphone while advancing the entire
source/field/duct assembly once per physical interval. `new_with_drives` binds
the matching observer clock; `new_resampled_with_drives` uses the existing
interval-mean observation clock without retiming the physics. A callback
failure preserves the accepted physical/output prefix and poisons the stream.

### Native source-to-microphone example

```sh
cargo run --release -p fs-couple --example nsf_thermoviscous_bore -- output.wav
```

The example uses an explicit named near-ambient air mixture in a manufactured
12 × 32 × 5-cell nonlinear flow region. A finite 10 mL chamber at a transverse
opening feeds an internal plenum, a two-cell-high windway, a seven-cell free
gap and a stationary staircase labium. All flow from the source plenum toward
the edge must pass through the windway. The velocity field comes from the
original nonlinear NSF solver; no prescribed jet velocity or jet card is used.

The actual NSF material law determines its timestep from the supplied cell
size. Eight original fluid substeps share each 48 kHz acoustic interval; the
coupling checks the derived clock instead of changing it. The field keeps the
same gas reference and original molecular source, and both ducts retain their
original molecular and 12-pole viscothermal states. The source chamber keeps
the original bounded named-mixture calorics used by the NSF region.

With cell size `dx`, the fluid occupies `[0,12dx] × [-32dx,0] × [0,5dx]` and
the FE region occupies `[0,6dx] × [0,2dx] × [0,5dx]`. Local integer indices
in the following table start at the low corner of the fluid volume:

| Part | Declared local geometry |
| --- | --- |
| Source opening | Low-x face, y indices 4–6 and z indices 1–3: nine original reservoir tiles |
| Internal plenum | Fluid y layers 3–7, before the windway |
| Windway | Fluid x indices 5–6 through y layers 8–12, across all five z cells |
| Free gap | Fluid y layers 13–19 |
| Labium | Solid y layers 20–29, with `6 <= x < 7 + floor((y-20)/3)` |
| Field opening | Lower half of the high-y face, x indices 0–5, adjoining FE port 2 |
| Escape opening | Upper half of the high-y face, retaining zero incoming characteristic |

The original source aperture's three-by-three tangential and two-layer inward
fluid collar is explicit. Both characteristic faces retain their required two
all-fluid layers. A solid wall at y index 2 isolates the unused low-y ambient
pocket, preventing source flow from escaping behind the windway. The two fluid
layers after the labium remain part of the original characteristic interface;
the staircase does not invent a solid edge inside that mandatory collar.

Left/right FE ports connect the 20 mm closed cork-side passage and 50 mm
radiating bore. No volume is counted in both source and acoustic owners.
The circular duct interfaces use one equal-area uniform mode and require
their own spatial comparison. The microphone observes the main bore's original
computed finite-aperture exterior 0.1 m from its outlet. The escape opening's
mass and energy remain in the NSF account; that opening is not silently assigned
a radiation impedance or added to this microphone as a point monopole.

The chamber starts at the common reference pressure and temperature. Its
volume rate ramps from zero to -1e-6 m³/s over 4 ms, remains there until 12 ms,
and returns to zero by 16 ms. A simultaneous heat supply follows the same
envelope with a 0.5 mW peak. The authored total volume reduction is 12 µL and
the total supplied heat is 6 µJ. Actual chamber pressure, temperature, mass
loss and pressure-volume work are responses of the retained material source.
Four final milliseconds continue the original dynamics with zero actuation.
The writer exports 20 ms of fixed one-pascal-scale PCM and reports clipping.
It refuses to replace an existing file.

Console output reports the accepted finite inventory, last original source
mass/total-energy exchange, midpoint heat/work receipt, source conservation
residuals, actual final windway axial velocity and prior acoustic/projection
accounts. The finite source admits at most 1% relative mass and energy movement
per interval, bounds volume movement independently at 1%, and checks original
mass and energy balances. The explicit acoustic projection allowance remains
5% of the sum of absolute exchanged works plus 1e-16 J. These declared limits
refuse inadmissible trajectories; they are not interface-transparency or
source-convergence qualifications.

The chamber, geometry, wall conditions, bulk-viscosity ratio and actuation are
declared numerical choices. The two-cell windway and staircase edge are coarse
geometry, not measured anatomy, a converged edge tone, a register-transition
result or a listening reference. Geometry, coupling/source partition and
source/observer bands require independent refinement; the actuation is not
asserted to be band-limited. The current flute and whistle descriptions are
not automatically rerouted through this example.

Native regression coverage includes original-owner replay, molecular and
finite-wall evolution, downstream and finite cork-side feedback, the full
cross-port solve, aggregate energy and all-owner rollback, binding admission
and callback-independent playback at selected ports. Rust syntax was parsed
with the pinned formatter during implementation. Compilation, native numerical
tests and example execution were not run because this executor lacks the full
pinned dependency checkout and remote build runners. No audible improvement or
instrument milestone is claimed by these changes.
