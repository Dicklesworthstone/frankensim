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
supported. Drives are original boundary pressure/characteristic controls or
original transverse reservoir states. This constructor does not silently
convert an existing finite chamber or finite upstream acoustic attachment into
prescribed data. Such compositions require their own joint attachment.

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

The new example uses an explicit named near-ambient air mixture in a
manufactured 8 × 8 × 2-cell nonlinear flow region with one interior constriction.
The actual NSF material law determines its timestep from the supplied cell
size. Eight original fluid substeps share each 48 kHz acoustic interval; the
coupling checks the derived clock instead of changing it. The field keeps the
same gas reference and original molecular source, and both ducts retain their
original molecular and 12-pole viscothermal states.

With cell size `dx`, the fluid occupies `[0,8dx] × [-8dx,0] × [0,2dx]` and
the FE region occupies `[0,8dx] × [0,2dx] × [0,2dx]`. The common opening is
the fluid's high-y face and FE port 2 at y=0. Left/right FE ports connect the
20 mm closed cork-side passage and 50 mm radiating bore. No volume is counted
in both source and acoustic owners. The circular duct interfaces still use
one equal-area uniform mode and require their own spatial comparison.

An authored pressure offset on the opposite fluid boundary rises to one pascal
over 2 ms, remains there until 6 ms and releases by 8 ms. The field and duct
pressures are solved responses to that physical boundary control. The existing
computed microphone 0.1 m from the main outlet supplies 10 ms of fixed-scale
PCM output. Console output retains original mass/energy differences and the
acoustic/projection residuals. The explicit projection allowance is 5% of the
sum of absolute exchanged works plus 1e-16 J; this exploratory budget is not
an interface-transparency qualification.

The constriction, wall conditions, bulk-viscosity ratio, geometry and controls
are declared numerical choices. They are not a measured embouchure, free jet,
flute, register-transition result or listening reference. The current flute
and whistle descriptions are not automatically rerouted through this example.

Native regression coverage includes original-owner replay, molecular and
finite-wall evolution, downstream and finite cork-side feedback, the full
cross-port solve, aggregate energy and all-owner rollback, binding admission
and callback-independent playback at selected ports. Rust syntax was parsed
with the pinned formatter during implementation. Compilation, native numerical
tests and example execution were not run because this executor lacks the full
pinned dependency checkout and remote build runners. No audible improvement or
instrument milestone is claimed by these changes.
