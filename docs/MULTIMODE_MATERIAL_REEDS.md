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
