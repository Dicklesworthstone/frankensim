# Several material plate modes in one reed source

`bernoulli_aperture::modal::ModalAperture` retains 1..32 independent structural
coordinates in one pressure/contact/flow solve. It is accepted by the existing
`LoadedAperture`, `TwoSidedAperture` and `TwoSidedPressure` paths, including
regional gas, fractional/viscothermal propagation, material bore walls, exterior
radiation and multirate pressure observation. No downstream pitch or amplitude
normalization is introduced.

The file-driven two-sided source also accepts `supply circular LENGTH_M RADIUS_M
RADIAL_CELLS`, alongside the existing planar `supply slit` declaration. A shared
finite-volume momentum owner resolves either planar strips or cylindrical
annuli from the complete inlet gas. The circular radius therefore changes real
inertial pressure, shear loss and retained flow profile, while gas density and
viscosity propagate from the same material state used at both reed faces.
Channel work enters the original coupled energy account once; physical viscous
heat remains separate from backward-Euler numerical loss. This is a fixed,
fully developed, incompressible laminar passage, with independent spatial and
temporal refinement. It does not model taper, entrance flow, moving walls or a
complete measured staple. Do not include the same passage in a duct network.

For modal displacement q and midpoint velocity v, the declared ports are:

```
gap = H + g.q
pressure generalized force = -A * delta_pressure
swept volume flow = -A.v
contact generalized force = g * contact_force
external generalized force = B * [force_0, force_1]
```

Pressure and swept volume share A, but g is a separate slit projection. A mode
whose mean slit displacement is zero can still move fluid. Signed and scaled
basis changes must transform all ports and mass/stiffness together; the physical
pressure, gap, work and energy then remain unchanged. The original Bernoulli and
contact implementations own their constitutive behavior. The scalar midpoint
condensation solves the entire source, not a series of independent reed voices.

## Original material plates

`ModalAperture::from_material_plates` accepts one or two pairs of
`(ResolvedPlateChart, ModalPlateOptions)`, the total rest opening, one coherent
inlet gas, the actual pressure/volume-flow impedance, fixed mechanical timestep,
lifetime limit, existing `PairedApertureLimits`, all initial modal states, one
explicit effective-slit `Obstacle`, and a cancellation gate.

The chart comes from `compile_plate_material_chart`, with original regional
source cards, element assignments, thickness/mass constraints and orientations.
Its optional `with_kelvin_voigt_bending_loss` selection is preserved. The modal
options reuse `PlateChartRadiation`: support/prestress, eigenwindow, requested
mode count, original nodal force footprint and explicit authored damping. When
source material damping is selected, authored damping must remain zero. The
existing diagonal weak-loss projection supplies the modal viscous coefficients.

The existing `certified_radiators_with_force_ports` performs a SINGLE eigensolve
per blade. Nodal unit-force ports extract every P1 displacement trace on that
same basis. Actual modal mass, stiffness, selected damping, signed pressure area,
mean slit motion and the blade's original force footprint are then transformed
together into metre-valued coordinates. Modes are not divided by their mean slit
displacement, which can legitimately vanish. The requested mode count must be
resolved; missing modes do not silently disappear.

Initial states and returned states are ordered by plate, then retained mode.
The first physical force input drives the first plate's original footprint;
the second drives the second plate, if present. These fixed footprint forces can
be supplied through the normal two-sided drive. The existing single-mode
`ScheduledTwoSidedPressure::with_plate_forces` attachment does not accept this
new basis; it refuses rather than applying a single-mode projection to it.
Canonical pressure phrases and held physical force inputs remain available.

`material_plates()` preserves source cards/options and the signed nodal basis.
`nodal_motion(plate, node)` reconstructs actual displacement and velocity from
the accepted modal coordinates. A density, modulus, geometry or selected
viscosity substitution is therefore a source-data change before construction,
not a change to the pressure/contact or acoustic-network code.

## Held material environments

Material source rows now accept the optional suffix
`environment COUNT AXIS VALUE ...`. For example, a source-backed plate
section can declare its own measured moisture coordinate:

```text
section blade.matpack isotropic 293.15 1000 0.0003 0 environment 1 moisture_content 0.12
```

The same suffix follows the existing fields of `solid_material`, graph
`wall`, uniform or regional DG material-wall inputs, and `boundary_plate`.
It changes the query of the original material card before any structural or
acoustic realization. It does not change pack identity, geometry, loss
selection, or the physical source/observer clocks. Different sections, walls
or keys can declare independent values even when they share one source pack.

These rows use the same `fs-material` query constructor as material strings.
Temperature remains absolute kelvin and the existing frequency coordinate
remains angular frequency in radians per second. Up to eight additional
coordinates must be explicitly typed by every declaring source claim as
static pressure, composition, moisture or concentration; conflicting,
missing, duplicate, nonfinite and out-of-domain states refuse. The source's
own normalized units and conventions are retained. Gas relative humidity
does not automatically become solid moisture, and a declaration is not a
sorption or transport model.

The original mass, elasticity, full plate viscosity, regional Maxwell memory,
solid tensor memory and compliant duct-wall owners consume that one resolved
state. Frequency-dependent material response is not frozen into a broadband
constant; the original loss law and its frequency-applicability checks remain
in force. Source coordinates are retained in the resolved material identity,
including the existing native material metadata.

Focused native regressions exercise manufactured moisture curves through the
same reed motion and microphone, characteristic and DG wall constitutive
response, and the hereditary solid-valve path. These numerical source laws
are not measured wet-cane data or an instrument-accuracy qualification.

## Full material viscosity, including cross-mode reactions

`ModalAperture::from_viscous_material_plates` accepts one or two tuples of
`(ResolvedPlateChart, ModalPlateOptions, PlateViscosity, PlateViscosityLimits)`.
Its remaining physical/state/contact arguments are the same as
`from_material_plates`. The original diagonal-loss constructor stays unchanged.

Create each `PlateViscosity` with `PlateViscosity::from_chart` and the six explicit
reduced plane-stress viscosity property names. The same shared material owner
used for instrument bodies resolves region orientation, thickness, source identity
and elastic/viscous applicability. Elastic properties do not imply viscosity;
missing data, inconsistent source identity or incompatible loss selection refuses.
Separate authored modal damping and scalar material bending loss must be absent.
An explicit zero tensor permits one blade to remain lossless.

The existing `certified_viscous_radiators_with_force_ports` projects the complete
viscous matrix on the SAME eigenbasis that supplies every physical port. For
mass-normalized C and physical coordinate conversion Q = S q, the reed receives
S^T C S. Both matrix axes are transformed, including signs and cross terms.
Matrices for independent blades occupy separate blocks; fluid/contact couples
the blades without inventing a solid-to-solid dashpot between them.

The running valve factors M + (dt/2) C + (dt/2)^2 K once through the shared
linear-algebra owner. Free motion, pressure and contact use that SAME response
inside the original nonlinear slit solve. The final force is C v_mid and the
accounted loss is dt v_mid^T C v_mid. No old-velocity force, split damping step,
modal refit or output filter substitutes for this reaction. Full pressure, force,
energy and spatial checks still precede accepted-state publication. Construction
and each linear solve reject unresolved arithmetic instead of adding a loss floor.

`material_plates()[i].viscosity()` retains the original regional tensor source;
`viscosity_matrix_n_s_m()` exposes the admitted full matrix, locally per plate or
for the complete valve. The lower-level `with_viscosity` accepts an authored
matrix but does NOT confer material identity. Explicit undamped rows are allowed;
the remaining block must be positive definite. Other nonzero singular matrices
are currently refused, not regularized. Zero full loss uses the original arithmetic.

The existing upstream/downstream networks, source-pressure scheduler, material
walls, exterior reaction and multirate microphone require no alternative runtime.
Changing viscosity data now changes that coupled mechanical trajectory. This
still requires temporal refinement: the undamped omega*dt guard is not a damping-
rate error certificate. A fixed Kelvin-Voigt tensor is not a Maxwell history,
a measured wet-cane law, or proof of real-instrument fidelity.

## Explicit approximation and failure boundaries

This is a fixed, linear plate basis with an explicitly selected diagonal or full
viscous damping model and one uniform effective slit. It is NOT the existing independent-blade hereditary or
spatial-contact image. Maxwell memory, curved shells,
distributed opposing-surface closure, nonlinear blade mechanics, live material
changes and calibrated instrument timbre are not claimed.

`maximum_trace_slope` and `maximum_slit_deviation_m` are checked on the SUM of
all modes at each candidate endpoint, per original blade. They bound the P1
pressure/displacement trace, not unexported DKT rotations or continuum error.
Both plates must supply the same physical slit width. Their local frames must
point toward increasing total opening. `maximum_angular_step` limits every
retained omega*dt; geometry, mode and projection budgets bound cold reduction.
The existing eigensolver remains synchronous, bracketed by cancellation checks.

A rejected mechanical, spatial, flow or energy trial publishes none of the
modal state. In the two-sided system neither wave network advances independently
of a rejected source. A pressure observer failure retains the existing poisoned
performer semantics, not a fabricated zero-pressure sample.

## Focused verification

```
cargo test -p fs-couple --test modal_aperture -- --test-threads=1
```

The original eleven authored tests cover original paired-mechanics correspondence, signed
basis changes, zero-gap swept-volume modes, refusal/horizon behavior, independent
wave-owner correspondence, multirate pressure, source identity, material density
and viscosity changes, spatial guards and actual material-pair microphone output.
They have NOT been executed in the implementation container: cargo/rustc are
absent. The preceding instrument CI stopped at dependency resolution, not a
passing instrument test. No native or measured-fidelity pass is inferred.

Independent SciPy references executed 240 simultaneous force/flow cases, 240
signed basis changes, 4096 source steps, 512 generalized-eigenbasis port cases
and 2048 two-sided network steps. Maximum source velocity disagreement was
1.34e-15 m/s and combined network energy residual 4.47e-23 J. These validate
reference equations, not Rust execution, the DKT implementation or real reeds.

Seven additional Rust regression sources cover full matrix mechanics, manufactured
momentum/flow/work, cross-mode transfer, signed congruence, passive decay, source
projection, data-only viscosity substitutions, 4x multirate pressure, zero-loss
parity and invalid bindings. They are part of the SAME `modal_aperture` target.
They were authored but not run in the implementation container (no Rust tools).
Independent full-matrix reference checks ran 160 simultaneous contact/flow cases,
160 signed basis cases and 1024 isolated steps. Further references ran 256
full-matrix basis projections and 3072 six-mode two-network steps, with a maximum
combined energy residual of 8.38e-24 J. These do not execute Rust, DKT, or the
pressure renderer and are not experimental validation.

### Materials on finite DG side branches

The DG flute and material-reed performances accept material walls on the actual
finite chimneys as well as the main bore. Inside the existing `dg_wall` block:

```text
dg_wall assembly 2
dg_wall_main rigid
dg_wall_chimney 7 isotropic wall.fsmcdpk 293.15 1000 0.00002 stress-free 50 2000 kelvin-voigt
```

An assembly declares exactly one main duct and distinct original finite-chimney
port IDs. Each target accepts the same rigid, isotropic, orthotropic or
`regions N` source syntax; a region's axial coordinates start at that target's
inlet and must align with its own DG element faces. Omitted chimneys retain the
existing rigid-wall model. This is a source assignment, not a new material law:
each target uses the original local membrane/DG owner with its actual radius,
regional gas and viscothermal history. All gas, radial wall and hereditary
states count toward the complete DG allowance. Material artifacts are shared
across the assembly and, in the reed path, the blade source cache.

Changing only a chimney material therefore changes its coupled impedance,
radial motion and returning pressure in the same source-to-microphone path.
The existing thin cylindrical membrane, frozen material state and linear
acoustic limits remain; this does not resolve axial shell waves or exterior
radiation from the solid wall.

## Upstream tract on the material-modal DG source

The optional `upstream_tract` block follows the downstream DG bore and its
optional material-wall declarations, before `slit_length`. It contains another
complete `dg_tube` or `dg_bore` source and optional `dg_wall` layout. For example:

```text
upstream_tract
dg_tube 0.04 0.001 closed 1 2 8 50 2000 0.01 1e-11 256
dg_wall isotropic blade.matpack 293.15 100 0.00002 stress-free 50 2000 maxwell 1 1 0.5 boundary_modulus boundary_time
```

These dimensions are an illustrative source declaration, not a measured vocal
tract. The shared DG reader also accepts a tapered geometry, axial gas field,
finite side chimneys, terminal passages, computed exterior reaction and regional
or chimney material walls. Every source retains its own explicit DG, branch,
terminal and material budgets. Both loads use the source's mechanical clock.
Tract, bore and reed material declarations share the same immutable pack cache.

The tract inlet points away from the valve. The same current through-flow leaves
the tract and enters the bore, including the original work-conjugate swept flow
of the moving material reed. The two exact mean-port impedances add inside the
existing nonlinear modal solve. Tract geometry, gas and solid-wall memory
therefore change the actual blade motion and bore sound. The original one-sided
source remains unchanged when the block is omitted.

Blowing pressure is an ideal series pressure rise with work `P_source Q dt`;
tract feedback is never counted as another prescribed source. Both acoustic
candidates and the original mechanics are validated before any accepted history
advances. Frames retain the actual pressures, signed work and terminal exchange
of each side. Inlet-pressure and baffled-outlet observations continue to refer
to the downstream bore. Upstream terminal radiation is reported in its own
exchange ledger and is not silently mixed into that microphone.

Both valve-face gas states must match completely; remote gas fields may differ
at a common equilibrium pressure. This is the existing small-signal series-tract
approximation, with frozen source composition and held mean flow per step. It
does not add lungs, transported mixture exchange, live tract shape changes or
resolved turbulent mouth flow. The original DG inferred-loss qualification and
independent spatial/time convergence requirements still apply.
