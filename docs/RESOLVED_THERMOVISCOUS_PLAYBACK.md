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
trial preserves every accepted owner. A separate complete upstream acoustic
assembly can coexist with this finite transverse supply.

This is still an acoustic projection about zero reference flow. The FE region
does not receive mean mass, species or advected enthalpy. The complete NSF mass
and total-energy accounts, streamed-mass projection difference and signed
pressure/flow covariance remain separate observable quantities. Covariance is
never supplied as invented heat to close the acoustic account.

### Complete acoustic feedback on both sides of the source

`with_upstream_load(load, solve, gate)` attaches a second original
`ThermoviscousDgLoad` to the opposite characteristic face, or to its explicitly
selected source patch. The upstream FE field keeps its complete cross-port
response, original molecular and thermal history, and every attached finite
DG passage. Authored nonzero initial states are retained. The two acoustic
assemblies and nonlinear source must represent distinct, non-overlapping
physical volumes; the caller supplies that geometry.

Attachment requires cold owners, matching opposite-face reference gas,
interface area and original interval duration. An upstream load with unbound
mouth coordinates is refused. Both this attachment and any downstream mouth
attachment must supply the same explicit `NsfDgRadiationSolve`; one
`joint_solve()` policy owns the common nonlinear solve and accepted warm start.
The coordinate order is downstream characteristic, independent escape-mouth
pressures, then upstream characteristic. With no mouths, the downstream
field's scalar affine response enters the same original two-face solve.
Without upstream or mouth attachments, the existing scalar path is preserved.

`step_with_upstream_pressure_rise(rise, gate)` selects an ideal series source
between the upstream field and the resolved zone. Its convention is
`p_zone - p_upstream = rise`; with `Q_out` directed out of the nonlinear zone,
the source supplies `-rise * Q_out * dt`. `Boundary` controls cannot overwrite
the solved upstream characteristic. An attached finite material supply keeps
its exclusive finite-source drive forms, so its passive upstream assembly has
zero series rise. Prescribed transverse supply can also coexist with this
passive feedback path.

Every candidate interval previews the downstream field, its reciprocal
exterior, and the upstream field against the same original NSF trial. All
interface work, complete-field balances, finite-source checks and the final
cancellation poll precede publication. Fluid, both fields, all ducts, exterior
memory and finite inventory then accept together; only those accepted controls
warm-start the next interval.

`upstream_load().system().last_frame()` retains the actual complete upstream
field and duct states. `frame.upstream` reports its selected interface, flow,
pressure, inlet work, series-source work, stored acoustic energy and signed
pressure/work/projection defects. The existing top-level acoustic storage and
balance fields remain downstream-only; upstream totals have their own receipt.
The upstream acoustic projection does not absorb mean NSF mass or enthalpy and
does not by itself establish a resolved player tract.

### Reciprocal escape mouths and their original exterior

`NsfThermoviscousDg::with_reservoir_radiation_at(aperture, input, solve, gate)`
attaches one passive transverse opening to the primary duct's existing
computed exterior. `with_reservoir_radiations` attaches several explicit
aperture identities to that same owner. Each supplied
`NsfDgReservoirRadiation` declares a common-baffle placement, rectangular
quadrature, circular-terminal quadrature and admissible ghost-pressure
bracket. The actual aperture tiles determine the physical rectangle and its
area. Every passive aperture must be distinct from the finite or prescribed
breath supply. Attachments are cold; a complete-field load with preattached,
unbound mouth coordinates is refused by the NSF constructor.

The resulting affine response has one coordinate for flow into the external
FE opening and one independent outward flow per passive mouth. The field's
complete cross-port matrix, all finite passages, and every original exterior
self/mutual term remain in this response. The shared Newton–Krylov pressure
solver evaluates the original private NSF interval with all those controls
together. It retains the declared characteristic/pressure brackets, charges
all source evaluations to the original work budget, and checks final
pressure and conservative pressure-work agreement separately at every mouth.
The solved controls warm-start the next interval only after the common
acceptance of fluid, field, ducts, exterior memory and finite supply.

This reaction has one uniform acoustic-pressure mode per rectangle. It does
not make the original prescribed-ghost boundary transparent, replace the
original source material law, or transfer mean material into the acoustic
field. `frame.reservoir_radiations` retains each mouth's original conservative
NSF exchange, solved ghost-pressure offset, actual outward flow, exterior
pressure/work and projection defects. The complete FE/DG storage account
retains its original signed circular-terminal work. Direct mouth work enters
the shared exterior separately and is not debited again as field or duct
storage. The finite supply still receives only its own aperture's original
streamed mass and total-energy exchange.

`nsf_dg::thermoviscous::performance::NsfThermoviscousDgPerformance` implements
the existing `PressureRenderer`. It observes the actual selected attached
duct's mean inlet pressure or computed microphone while advancing the entire
source/field/duct assembly once per physical interval. `new_with_drives` binds
the matching observer clock; `new_resampled_with_drives` uses the existing
interval-mean observation clock without retiming the physics. A callback
failure preserves the accepted physical/output prefix and poisons the stream.

`new_with_drives_at` and `new_resampled_with_drives_at` accept
`NsfThermoviscousDgObservationPort::Downstream(port)` or `Upstream(port)`.
The port remains the original FE identity of an attached duct. The selected
assembly supplies its actual inlet pressure, terminal flows, computed receiver
and gas propagation. Both assemblies advance once at the shared physical
clock regardless of the observation target, and both physical horizons must
cover the requested stream. Existing constructors retain their downstream
meaning. `observation_target()` exposes the assembly and port;
`observation_port()` retains the local FE port identity.

`with_reservoir_observation_at(aperture, placement)` adds that selected
opening's actual signed tile flows to the existing microphone through the
original finite-tile propagation owner. A reciprocally loaded mouth requires
the primary duct's shared exterior, the same baffle origin, zero additional
rotation and the declared circular-terminal observer quadrature. A supplied
breath aperture is not implicitly selected. The observer retains each tile's
flow and path rather than substituting bore flow or streamed mass divided by
reference density. With several selected openings their signed pressures are
summed at the same output clock. Observer gas preparation, when selected,
applies its existing shared path budget after all surfaces are attached and
aligns numerical padding without changing physical delays.
These escape surfaces belong to the downstream exterior; reservoir observation
is refused for an upstream target rather than assigning it unrelated geometry.

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
A separate manufactured escape in the free gap shares reciprocal radiation
with the bore. This extended example has not been compiled or executed here.
A second complete cavity and closed passage connect the low-y source face to
the plenum, allowing upstream acoustic feedback to participate in the same
resolved interval.

The actual NSF material law determines its timestep from the supplied cell
size. Eight original fluid substeps share each 48 kHz acoustic interval; the
coupling checks the derived clock instead of changing it. The field keeps the
same gas reference and original molecular source in both acoustic fields, and
all three ducts retain their original molecular and 12-pole viscothermal states.
The source chamber keeps the original bounded named-mixture calorics used by
the NSF region.

With cell size `dx`, the fluid occupies `[0,12dx] × [-32dx,0] × [0,5dx]` and
the downstream FE region occupies `[0,6dx] × [0,2dx] × [0,5dx]`. The upstream
FE cavity occupies `[0,12dx] × [-36dx,-32dx] × [0,5dx]`. These three interiors
do not overlap. Local integer indices in the following table start at the low
corner of the fluid volume:

| Part | Declared local geometry |
| --- | --- |
| Source opening, aperture 0 | Low-x face, y indices 4–6 and z indices 1–3: nine original reservoir tiles |
| Internal plenum | Fluid y layers 0–7, before the windway |
| Upstream field opening | Complete low-y face, adjoining upstream FE port 1 at y = -32dx |
| Windway | Fluid x indices 5–6 through y layers 8–12, across all five z cells |
| Free gap | Fluid y layers 13–19 |
| Passive escape, aperture 1 | Low-x face, y indices 16–18 and z indices 1–3: nine distinct reservoir tiles |
| Labium | Solid y layers 20–29, with `6 <= x < 7 + floor((y-20)/3)` |
| Field opening | Lower half of the high-y face, x indices 0–5, adjoining FE port 2 |
| Ambient remainder | Upper half of the high-y face, retaining zero incoming characteristic |

Both transverse apertures retain their original three-by-three tangential and
two-layer inward fluid collars. The passive escape is downstream of the
windway and upstream of the labium; it cannot bypass the windway from the
source plenum. Both characteristic faces retain their required two all-fluid
layers. The former back wall at y index 2 is removed so the upstream acoustic
field is connected to the source plenum. This face now has a solved acoustic
load, with its actual NSF mass/energy and projection account retained. The two
fluid layers after the labium remain part of the original characteristic interface;
the staircase does not invent a solid edge inside that mandatory collar.

Left/right downstream FE ports connect the 20 mm closed cork-side passage and
50 mm radiating bore. Upstream FE port 0 at y = -36dx connects a separate
40 mm closed passage directed outward from the source. Its circular interface
matches the area of the complete upstream FE cross-section. No volume is
counted in both source and acoustic owners.
The circular duct interfaces use one equal-area uniform mode and require
their own spatial comparison. Aperture 1 joins the primary bore's original
computed self/mutual exterior in an explicitly manufactured common rigid
baffle. Its grid tangential edge zero maps to `[-0.08, 0]` m, with the original
low-face reversal of the second coordinate. This placement does not claim a
free-field instrument body. The declared rectangular reaction quadrature is
16 × 16; circular outlets use eight radial rings and 32 angular points.

The microphone at `[0, 0, 0.1]` m uses those same circular quadrature counts,
the same escape origin and zero extra rotation. It receives the signed sum
of original bore-outlet pressure and all nine escape-tile pressures through
their respective original propagation paths. The finite supply aperture is
not observed as an escape. The high-y ambient remainder keeps its explicit
zero incoming characteristic and its original NSF mass/energy account; it
has no silently assigned radiation load or added microphone source.

The simultaneous three-coordinate pressure solve includes the downstream
characteristic, passive-mouth pressure and upstream characteristic. Both
characteristics use a bracket of ±0.5 m/s; the passive-mouth ghost-pressure
bracket is ±20 Pa about the source reference. The shared policy uses a
normalized finite-difference step of 1e-5, a 1e-8 Pa nonlinear absolute tolerance,
a three-vector linear restart and at most 12 Newton iterations.
The existing allowance of 64 original NSF interval evaluations also bounds
numerical Jacobian and line-search work. These are explicit unqualified
numerical settings; a trial outside a bracket or allowance refuses and does
not clip the physical solution or advance only part of the assembly.

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
accounts. The upstream receipt separately reports pressure, outward source flow,
inlet and series-source work, complete acoustic storage, balance and interface
defects. It also reports aperture 1's actual outward flow, shared-exterior
pressure/work, conservative streamed mass/total energy, pressure residual
and pressure-work covariance. The finite source admits at most 1% relative mass and energy movement
per interval, bounds volume movement independently at 1%, and checks original
mass and energy balances. The explicit acoustic projection allowance remains
5% of the sum of absolute exchanged works plus 1e-16 J. These declared limits
refuse inadmissible trajectories; they are not interface-transparency or
source-convergence qualifications.

The chamber, geometry, wall conditions, bulk-viscosity ratio and actuation are
declared numerical choices. The two-cell windway, staircase edge and upstream
cavity/passage are coarse manufactured geometry, not measured player anatomy,
a converged edge tone, a register-transition result or a listening reference.
Geometry, coupling/source partition and
source/observer bands require independent refinement; the actuation is not
asserted to be band-limited. The current flute and whistle descriptions are
not automatically rerouted through this example.

Native regression coverage includes original-owner replay, molecular and
finite-wall evolution, downstream and finite cork-side feedback, two complete
fields at whole or selected interfaces in either orientation, the full
cross-port solve including independent mouth cross terms, aggregate energy
and all-owner rollback, exclusive binding admission, signed tile/bore
microphone composition and callback-independent playback at selected ports of
either assembly. Authored upstream cases also exercise finite-supply/escape
feedback, both-field rollback and original material-microphone replay on a
separate 8 kHz to 12 kHz observation clock.
These are algebraic/replay and rollback regression cases (G0/G3/G4); they do
not establish spatial or temporal convergence. Rust syntax was parsed
with the pinned formatter during implementation. Compilation, native numerical
tests and example execution were not run because this executor lacks the full
pinned dependency checkout and remote build runners. No audible improvement or
instrument milestone is claimed by these changes.


## Native two-port resolved-field source

The schema `frankensim-resolved-thermoviscous-dg-performance-v1` now provides a
first native composition rung for the existing nonlinear NSF source and complete
three-dimensional thermoviscous field. It embeds one ordinary
`frankensim-resolved-dg-performance-v1` source after a literal
`resolved_source` separator. The ordinary source is admitted by its existing
reader first; only then can its cold NSF grid and primary DG load move into the
field composition.

The field prefix declares an SI tetrahedral mesh, an explicit wall/port
partition, a positive frequency band, an independently sourced bulk/shear
viscosity ratio, and two distinct port identities: the NSF source port and the
port receiving the original primary DG load. Every exterior triangle must be
owned by a wall or port. The source port must lie exactly on the selected NSF
face or patch, remain inside its tangential bounds, and the entire FE volume
must lie outside the nonlinear-fluid volume. The existing field/DG and NSF/field
constructors independently recheck area, gas, clock and work-conjugate matching.

This first file rung intentionally accepts only a cold direct source without
finite-reservoir, upstream-DG, reciprocal-mouth, reservoir-observer or prepared
gas-observer state. Those owners refuse the cold handoff rather than being
discarded. The embedded source can still use the existing pressure,
pressure-temperature, SGS, SI geometry, gas-mixture and aperture declarations
that do not attach one of those extra histories. Additional field ducts,
molecular field reservoirs and finite wall sources remain explicit follow-on
schema extensions.

`crates/fs-couple/examples/resolved-thermoviscous.performance` is a manufactured
integration fixture: a 3 mm long 3-D field sits immediately outside a 6 mm
resolved source volume and feeds the original 50 mm radiating DG duct. The warm
pressure history reaches the existing conservative enthalpy contact, resolved
flow, complete thermoviscous field, duct and microphone in that order. It is not
measured flute geometry, a converged jet or a listening-quality claim.
