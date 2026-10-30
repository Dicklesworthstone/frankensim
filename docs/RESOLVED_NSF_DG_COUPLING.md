# Resolved compressible flow and the original DG bore

`fs_couple::nsf_dg::NsfDg` connects an original material-owned
`fs_lbm::nsf3::NsfGrid3` to an original `DgPortLoad`. An initial resolved
pulse, vortex or interior obstacle field stays in that fluid owner. Incoming
pressure feedback from the evolving bore reaches the original NSF boundary and
populations during the same accepted coupled interval.

This is a reusable acoustic interface for resolved source zones. It does not add
a flute jet, a reed law, a fitted source coefficient, an output envelope or an
instrument-specific solver.

## Physical interface

The fluid zone must already have its original characteristic axis. Select either
complete face, or an explicitly declared rectangular aperture on that face;
its outward normal flow is positive into the DG inlet. The other
face receives either a caller-supplied incoming characteristic in metres per
second or a physical pressure relative to its original reference in pascals.
For a zero-velocity reference, the original boundary uses

`J− = u_normal − (p − p_reference)/(rho_reference c_reference)`.

With a prescribed external boundary, only the incoming characteristic on the coupled face is solved. For every trial,
the original NSF owner advances from the same accepted state while the original
DG owner supplies its current step-mean law `p_DG = H + Z Q`. The scalar solve
matches the mean reference-relative NSF face pressure to that law, using the
face's actual integrated continuum volume flow.

### Bore aperture and surrounding escape

`NsfGrid3::with_characteristic_patch` declares one rectangular subset of the
original face quadrature before stepping. `NsfDg::new_on_patch` couples that
subset to the original DG inlet. Its actual area must match the inlet area;
the constructor never resizes the bore or infers an aperture from it. The two
half-open bounds run along the increasing Cartesian axes other than the normal
axis. The patch retains the original face gas and two-cell fluid collar.

The remaining cells keep an independently declared pressure or incoming-wave
condition. For example, `PressureOffset(0.0)` leaves the surrounding opening at
the reference ambient pressure while bore pressure feeds back only through the
selected aperture. The original resolved field can consequently divide its
flow between the bore and the surroundings. This is one fluid solve and one
boundary reconstruction, with no empirical split fraction or added source.

The root solve uses only the patch's mean pressure and integrated volume flow.
Its acoustic-work, covariance and mass-projection checks use the same selected
quadrature. `NsfDgFrame::coupled_patch` exposes that exchange; `faces` retains
both **complete** face exchanges, already including the aperture. Never add the
patch account to the complete face account. The other complete face can still
carry the finite upstream acoustic load. Patch controls, fields, exchanges and
both duct histories publish together only after all checks succeed.

This enables a split boundary topology, not a complete flute mouth. It does not
create a labium wall at the patch edge, resolve an unmeshed jet, remove the
original first-order plane-normal characteristic approximation, or admit SGS
through open boundaries. The complement remains an ideal prescribed boundary;
escaping flow has no separate exterior-radiation observer or mutual radiation
reaction. Current pressure playback still observes the main DG bore. Geometry,
source-flow, interface and acoustic convergence remain required.

### Transverse reservoir opening

An original `NsfReservoirAperture` can also open explicit tiles in a stationary
wall on a different axis. It coexists with either the complete-face bore or a
selected bore patch. Its density, Cartesian velocity and temperature are
physical inputs; the original NSF material law supplies pressure, calorics and
transport. It is a prescribed reservoir with the same composition law as the
resolved volume, not an independently mixed atmosphere or transparent acoustic
termination. The canonical `examples/resolved-dg-side-aperture.performance`
and matching Instrument Lab starter expose this original owner through native
playback. Its declared dry-air density corresponds to 300 K and 100000 Pa;
changing other gas data does not silently change that explicit reservoir state.

Tiles must be unique fluid cells on the declared transverse face, clear of the
two-cell characteristic collars, other domain edges and solids, with the
original one-cell fluid halo. The backend retains its symmetric link ownership
at aperture/wall edges, characteristic feedback and conservative energy flux.
All these admission rules remain active in the native reader.

`NsfDgFrame::reservoir_mass_out_kg` reports the original streamed mass leaving
this opening, separately from both characteristic faces. The complete fluid
mass balance uses the two face fluxes **plus** this reservoir flux. Signed
reservoir energy entering the fluid remains `flow.open_energy_j`; it must be
included with `flow.open_boundary_energy_j`, wall heat and source work when
checking total fluid energy. Mixed-boundary unresolved internal exchange retains
the original ledger classification and is not reported as resolved viscous heat.
No reservoir flux is inserted into either acoustic duct or the bore microphone.

The signs and nonoverlapping accounts can be checked directly on each accepted
frame, up to numerical roundoff:

```text
flow_mass_change_kg + reservoir_mass_out_kg
    + faces[0].mass_out_kg + faces[1].mass_out_kg = 0
flow_energy_change_j = flow.open_energy_j + flow.open_boundary_energy_j
    + flow.wall_heat_j + flow.source_energy_j
```

The optional `coupled_patch` is already contained in its complete characteristic
face account. Native reservoir provenance names the mass and energy fields,
their SI units and their opposite sign conventions explicitly. Complementary
tests in `nsf_dg/performance/file/tests/aperture.rs` check original fluid/duct
replay, retained thermochemical sources and finite-upstream feedback together
with the opening. Numerical/interface regressions do not qualify a physical
mouth or its acoustic transparency.

### Physical pressure control

`NsfBoundaryDrive::PressureOffset(delta_p)` fixes the other face's pressure
at `p_reference + delta_p`. The original reconstruction uses its live outgoing
characteristic `J+` and the original material impedance `Z = rho_reference
c_reference`:

    J− = J+ − 2 delta_p / Z
    u_normal = u_reference_normal + J+ − delta_p / Z

This is evaluated at every original face quadrature point and every original
predictor/corrector stage. Returning acoustic waves therefore change source
flow while prescribed pressure stays fixed. Changing the gas's density,
calorics or temperature changes the physical response through the same material
owner. There is no fixed pressure-to-velocity conversion, second flow solve or
additional sound source.

`NsfGrid3::set_boundary_drives` and `preview_boundary_interval` retain the
original face temperature/entropy donor, tangential velocity, material queries,
population reconstruction, pressure force and conservative energy flux. Control
changes rebase the original force chart without changing physical density,
velocity or total energy. Rejected states and dropped/cancelled trials publish
neither a new control nor a partial field.

`NsfDg::step_with_boundary` selects this source condition while the existing
scalar interface solve still controls only the opposite DG face. `NsfDgFrame`
retains the typed external control alongside the actual source-face flow and
integrated work. A uniform held pressure has source acoustic work
`delta_p * Q_out * interval`; the original quadrature computes that work.
The existing `step(incoming_m_s)` remains the incoming-characteristic route.

This source is an **ideal prescribed-pressure boundary in the admitted local
linear characteristic domain**. It is neither stagnation pressure nor a finite
airway/reservoir's compliance, thermal capacity or mass inventory. It reflects
returning acoustic disturbances according to a fixed-pressure condition; it is
not a nonreflecting boundary. All original source, pressure-perturbation,
temperature and Mach restrictions still apply.

### Finite upstream acoustic feedback

`NsfDg::with_upstream_load` attaches a second original `DgPortLoad` to the other
complete NSF face before the first accepted interval. The upstream duct retains
its original wave, wall, molecular-memory and terminal states. Both duct inlet
flows point outward from the resolved zone. Returning upstream waves therefore
change the source-face pressure as well as its flow.

`NsfDgDrive::UpstreamPressureRise(rise_pa)` places an explicit ideal series
pressure actuator between this duct inlet and the resolved face:

    p_zone = p_upstream + rise_pa
    source_work = -rise_pa * Q_out * interval

The sign of source work is positive when this actuator supplies acoustic energy
to the combined system. The drive is a pressure **rise**, not a prescribed face
pressure. The upstream duct's instantaneous pressure remains solved. This is a
finite acoustic airway with a declared ideal actuator; it supplies neither an
anatomical lung model nor a finite reservoir's mass or thermal inventory.

`step_with_upstream_pressure_rise` solves both original incoming characteristics
against the two live current-interval DG inlet laws. A nested safeguarded scalar
solve uses one shared root-evaluation and original cell-step allowance. Every
trial starts from the same accepted NSF state. The final pair is restaged
through the original owners; both pressure residuals, both original quadrature
work accounts and cancellation are checked before the fluid and both ducts are
published together. An inadmissible bracket or exhausted shared allowance
refuses the whole interval without advancing any owner.

`NsfDgFrame::upstream` retains the upstream original step and terminal frames,
solved inlet pressure/flow, source work, and that face's pressure, covariance,
acoustic-work, port-work and mass-projection residuals. Its acoustic work residual
is `resolved_work_out + source_work - duct_inlet_work`. The existing main-face
accounts retain their original meaning. Both faces remain acoustic projections:
neither DG duct acquires NSF advected mass, entropy or enthalpy.

The complete reference gas state, including transport properties, must equal the
DG inlet gas. The actual selected area must match inlet area. The DG interval
must equal a positive integer number of original NSF steps, allowing only
floating-point roundoff in derived clocks. Both characteristic reference
velocities must be zero. These conditions prevent an unowned mean-flow interface
or independently chosen atmosphere from silently entering the acoustic path.
The original thermochemical model, caloric chart, transport support and all
original per-cell material refusals remain active.

Stationary halfway-bounceback walls may bound the transverse directions. The
two complete fluid layers at each characteristic face remain mandatory; solid
cells cannot occupy those collars. At a diagonal characteristic/wall corner,
the stationary wall takes precedence over an open link. Streaming and the
reported open mass/momentum flux use the same link selection. Finite-volume
energy, pressure work, conduction and wall heat retain their original owners.

Each candidate retains the complete NSF populations, force chart, energy,
material caches and histories. Both fluid and DG candidates remain unpublished
until root convergence, original-source admission, integrated port-work checks,
the caller's projection budget and cancellation checks all pass.

## Work and transported quantities

Every accepted frame exposes both original NSF faces, original fluid and DG
accounts, actual mass and energy changes, and original downstream terminal flows.
The NSF quadrature directly accumulates `integral((p-p_reference) u_normal) dA dt`;
it does not subtract large absolute-pressure work totals afterward.

The frame separately reports:

- **Projection covariance:** resolved acoustic work minus
  `mean(p′) mean(Q) interval`, including transverse and temporal covariance.
- **Acoustic work residual:** resolved acoustic work minus actual integrated DG
  inlet work. Its magnitude must fit the explicit absolute/relative budget.
- **Port work residual:** original DG inlet work minus its own mean-pressure and
  held-flow product.
- **Mass projection difference:** actual lattice mass flux minus
  `rho_reference mean(Q) interval`.

The acoustic duct does not transport the full NSF mass, entropy or enthalpy.
Actual NSF mass and total-energy fluxes remain visible and are not credited as
heat or forced into a false whole-system balance. The original NSF viscous-heat account and DG loss
accounts retain their existing budget-remainder qualifications.

## Pressure playback

Construct a cold `NsfDg` from the retained zone, retained load and explicit
`NsfDgLimits`. Limits bound grid cells, NSF substeps, all root-evaluation
cell-steps, iteration count, finite source duration, the incoming-characteristic
bracket and both pressure/work tolerances. Bracket endpoints must themselves be
admissible source states; an out-of-domain endpoint refuses rather than clipping
a physical amplitude.

`nsf_dg::performance::NsfDgPerformance::new` accepts that system, an
`ApertureObservation`, the matching sample rate, a finite `Vec<f64>` of incoming
characteristics on the other face, a maximum callback length and an
`Arc<CancelGate>`. One input is held per complete coupled interval. Choose
`ApertureObservation::Inlet` for the actual DG mean inlet pressure or
`TubeBaffled(receiver)` for the existing computed-radiation-bound microphone.
The microphone consumes the actual accepted terminal flows and retains its
physical/numerical delay.

`new_with_boundary` instead takes a `Vec<NsfBoundaryDrive>` so every control
value carries its physical kind and units. `new_resampled_with_boundary` uses
the same typed controls with the independent observation clock below.
`new_with_drives` and `new_resampled_with_drives` accept `Vec<NsfDgDrive>` for
either prescribed-boundary or finite-upstream control. `drive_history()` exposes
the exact typed controls; pressure offsets, series pressure rises and incoming
velocities retain distinct meanings. The control kind must match the attached
topology. The original scalar constructors select the characteristic variant and
retain their existing physical arithmetic.

Call `with_observer_gas(&propagation)` before playback to use the same physical
atmosphere law as native flute and reed observation. `GasPropagation` selects
the shared law and finite path/filter budgets; the actual DG radiation source
supplies the complete exterior gas, including viscosity, conductivity and heat
capacities. Each original surface-to-microphone path retains its own absorption,
dispersion and causal delay. `observer_gas()` exposes the original material and
path-response report. The fluid and bore reactions are unchanged, and the source's
existing cancellation gate also owns this cold preparation.

The result implements the existing `PressureRenderer`, so the original pressure
decimator, ensemble and PCM writers can consume it. Callback partitioning does
not change either physical owner or source timing. A failed physical step rolls
back both owners; a later observation failure preserves the accepted physical
prefix and poisons further playback.

### Independent physical and observation clocks

`NsfDgPerformance::new_resampled` accepts a fixed observation rate and output
sample count while retaining the original physical interval. Supply exactly
`source_intervals_for(dt, output_rate, output_samples)` incoming values, one per
**physical** interval. The original NSF step remains
`dx / sqrt(3 R T_reference)`; the DG interval is the declared integer number of
those steps. Changing gas composition or temperature therefore changes the
physical clock without changing grid geometry, source duration or playback rate.

The adapter integrates overlaps between original interval-mean pressure/terminal
flows and regular observation intervals. One physical frame can contribute to
several observation intervals, or several physical frames to one observation
interval. Only the original physical owners advance state. The shared Rayleigh
microphone consumes the projected terminal flows; exterior absorption and PCM
decimation then use their existing owners. Callback boundaries do not affect
the partitions. Matching clocks retain the original one-frame path, including
its arithmetic; derived-clock roundoff is admitted without retiming geometry.

This is a piecewise-constant temporal projection with a box-average frequency
response, **not a band-limited antialias filter**. Refine both clocks against the
intended bandwidth, and use the existing decimator for the final output where
appropriate. Receiver bandwidth must fit physical and observation Nyquist as
well as original wall/radiation bands. Physical state can lead the observation
window by less than one original interval. The source duration, physical prefix,
receiver delay basis and this approximation are explicit in native provenance.

`DgPortLoad::new_with_computed_baffled_branches_at_dt` also retains this exact
physical interval in the DG propagator and terminal realization. Its computed
radiation source uses `ceil(1 / dt)` solely as a finite source-band ceiling; actual
sampled-passivity admission uses `PI / dt`. Exact reciprocal integer clocks keep
the existing bake and arithmetic. The source ceiling never becomes a rounded
mechanical step.

## Native file and Instrument Lab

`nsf_dg::performance::file::ResolvedDgFilePerformance` assembles the original
owners from UTF-8 data. The native entry point is:

```sh
music_render resolved crates/fs-couple/examples/resolved-dg.performance output.wav --decimate --self-contained
```

For physical pressure control of the same resolved constriction, bore and
exterior microphone:

    music_render resolved crates/fs-couple/examples/resolved-dg-pressure.performance output.wav --decimate --self-contained

This second source declares a `+20 Pa`, then `−20 Pa`, then zero pressure history
relative to the authored 100 kPa source reference. The flow is solved. Changing
the source gas or pressure history is a data edit; it does not require changing
the solver or calculating a replacement characteristic velocity.

For finite upstream acoustic feedback, use:

    music_render resolved crates/fs-couple/examples/resolved-dg-upstream.performance output.wav --decimate --self-contained

This source retains the same resolved constriction and exterior receiver, adds a
closed upstream duct, and applies a `+20 Pa`, then `−20 Pa`, then zero **series
pressure rise**. Its face pressure and flow are both solved. Upstream length,
termination, wall material and admitted gas-memory declarations are source data.

The canonical example observes at 96 kHz; `--decimate` uses the shared pressure
decimator for 48 kHz PCM. Pressure remains in pascals and the authored full-scale
value controls PCM conversion; there is no peak normalization. The existing
`wind` schema dispatcher can also read the format. The explicit `resolved`
command refuses other schemas. Instrument Lab exposes this same native route,
raw-file import and attachment transport, including original solid-material
packs and gas-source associations. It does not translate resolved inputs into a
flue/reed closure. `--acceptance` remains unavailable for this source class.

The grammar below shows the fixed record order. Uppercase tokens are fields;
bracketed sections and comments explain the grammar and are **not literal input**.
Records have no ignored suffixes, blank rows or inline comments.

```text
frankensim-resolved-dg-performance-v1
clock OBSERVATION_RATE_HZ NSF_SUBSTEPS_PER_DG_INTERVAL
output OUTPUT_SAMPLES FULL_SCALE_PA
gas dry-air TEMPERATURE_K PRESSURE_PA
grid NX NY NZ DX_M PERIODIC_X PERIODIC_Y PERIODIC_Z
fluid_limits MAX_CELLS TMIN_K TMAX_K MAX_MACH MAX_TEMPERATURE_RATIO TRANSPORT_POINTS WORKERS
bulk_viscosity BULK_TO_SHEAR_RATIO SOURCE_TOKEN
wall adiabatic
characteristic AXIS low|high MAX_RELATIVE_PERTURBATION
[optional coupled_patch START0 END0 START1 END1 pressure REMAINDER_OFFSET_PA]
[optional reservoir_aperture AXIS low|high DENSITY_KG_M3 VX_M_S VY_M_S VZ_M_S TEMPERATURE_K COUNT]
[exactly COUNT reservoir_cell X Y Z records when an aperture is declared]
initial rest
solids COUNT
solid X Y Z
[repeat exactly COUNT solid records]
states COUNT
state X Y Z DENSITY_KG_M3 VX_M_S VY_M_S VZ_M_S TEMPERATURE_K
[repeat exactly COUNT state records]
coupling MAX_ROOT_EVALUATIONS MAX_CELL_STEPS JMIN_M_S JMAX_M_S PRESSURE_ABS_PA PRESSURE_REL WORK_ABS_J WORK_REL
dg_tube [existing shared DgSourceInput fields and subsidiary records]
[optional dg_molecular declaration from the shared molecular-law reader]
[optional existing dg_wall material-layout records]
[optional upstream_dg dg_tube|dg_bore followed by the shared source fields and subsidiary records]
[optional upstream_molecular declaration from the shared molecular-law reader]
[optional upstream_wall material-layout records]
observation inlet-pressure
[optional existing observer_gas record, only with a baffled observer]
boundary COUNT
hold END_OBSERVATION_SAMPLE INCOMING_CHARACTERISTIC_M_S
[repeat exactly COUNT increasing hold records]
end
```

Alternatives use the same original parsers and owners:

- `coupled_patch START0 END0 START1 END1 pressure PA` selects an explicit bore
  aperture immediately after `characteristic`; `characteristic M_S` can replace
  `pressure PA` for an incoming-wave condition on its remainder. Cell bounds are
  half-open, nonempty, and inside the original transverse grid. The selected
  aperture area must match the main DG inlet. Omitting the record preserves the
  existing complete-face path. The canonical
  `examples/resolved-dg-split-aperture.performance` and Instrument Lab starter
  exercise the native aperture, surrounding ambient condition and bore-outlet
  microphone together. This coarse manufactured fixture is not flute geometry.
- `gas thermochemical T_K P_PA mole|mass COUNT ID FRACTION ...` selects a complete
  source-backed thermochemical model. An `energy_solve ENERGY_J_KG_TOL MAXITER`
  record is then mandatory immediately after `fluid_limits`. Source packs and
  species associations come from the shared native gas options; the phase policy
  must explicitly be `FrozenGas`. A mixture snapshot alone cannot supply the
  temperature-dependent law for evolving cells.
- `wall isothermal TEMPERATURE_K` selects the original prescribed thermal wall;
  this is not a moving or finite-capacity solid. `initial uniform DENSITY_KG_M3
  VX_M_S VY_M_S VZ_M_S TEMPERATURE_K` replaces the reference-rest default. Each
  `state` record overrides one distinct fluid cell in SI units. Duplicate,
  out-of-bounds and solid-cell assignments are refused.
- `dg_bore` uses the existing shared axial geometry, regional gas and finite
  chimney records. `dg_wall` retains material-card identities and selected
  constitutive sources for the bore/chimney wall layout. The original load
  checks inlet area, complete inlet gas, material bands and work admission.
- `dg_molecular` selects the existing causal molecular-relaxation law for the
  original DG gas region, using exactly the same grammar as native reed/flue
  sources. Its species composition comes from the original thermochemical
  source record or axial field. A dry-air scalar snapshot has no species
  association and refuses this law; missing species, unsupported temperatures
  and insufficient memory allowances also refuse. This changes the DG gas
  region's memory; it does not add molecular internal-energy dynamics to the
  NSF volume or replace that volume's declared bulk-viscosity law.
- `observation baffled X_M Y_M Z_M RADIAL_RINGS ANGULAR_POINTS FMAX_HZ` binds the
  existing computed-radiation source and all its actual terminal openings. An
  optional `observer_gas` record selects the same atmospheric propagation law
  used by other native pressure renderers. The source, bore and exterior can
  have explicitly different temperatures/compositions at their admitted common
  equilibrium pressure; changing the whole atmosphere requires changing each
  explicitly authored region/exterior declaration.
  `observer_gas molecular MAX_PATHS MAX_FILTER_TAPS MAX_KERNEL_TERMS LAW`
  uses the shared material-relaxation grammar and the original `dg_exterior`
  species record. It does not borrow the resolved volume's composition or the
  upstream duct's exterior. A legacy exterior snapshot, missing species or an
  out-of-domain constitutive law refuses before preparing either physical owner.
- `boundary pressure COUNT` selects pressure offsets in Pa for each `hold`
  value. `boundary characteristic COUNT` explicitly names the original m/s
  condition. The existing bare `boundary COUNT` has the same characteristic
  meaning and arithmetic. Unknown kinds and trailing fields refuse.
- `upstream_dg dg_tube ...` or `upstream_dg dg_bore ...` attaches the second
  original load. Its subsidiary geometry, regional gas, chimney and radiation
  records use the same shared parser as the main load. `upstream_molecular` and
  `upstream_wall` retain the original gas relaxation and material-wall readers.
  Both loads must match their own complete NSF reference gas, selected area and
  actual coupled interval. An active reciprocal exterior-mouth port is refused
  at either interface. `boundary upstream-pressure COUNT` then declares series
  pressure rises in Pa. This control requires the upstream load; the prescribed
  boundary controls require its absence. The parser does not infer drive location
  or pressure meaning from a numeric history.

Cell indices are zero-based. Axes are `0`, `1`, `2`; periodic flags are `0` or
`1`, and the characteristic axis must be nonperiodic. The `low|high` field
selects the face coupled to the DG inlet. The `hold` history drives the opposite
face. Endpoints are exclusive indices on the **observation** clock and must
cover the declared window exactly. Each physical interval uses the held value
at its start, so a discontinuous control can incur up to one physical interval
of event timing error; that error refines with the original discretization.

The cold reader bounds bytes, cells, transport quadrature, workers, history,
root work and physical intervals before constructing the original solvers. It
polls the retained cancellation gate during parsing and around bounded cold
construction, then during physical solves and observation. The original NSF
constructor has no internal cancellation callback, so cancellation during that
one bounded construction is observed when it returns. No files are opened by
the library reader; native/Lab source loaders own attachment access.

The canonical file is a **synthetic numerical example**, with a coarse fixed
obstacle mask, an explicitly cited illustrative bulk-viscosity ratio and loose
projection-work tolerance. It is not a specimen description or measured flute
calibration. Its purpose is to exercise the actual fluid → bore → radiation →
material atmosphere → PCM path with data-owned geometry and source histories.

## Scope and refinement

This first composition is a uniform acoustic-mode projection with a held DG flow
and the original first-order characteristic boundary. Its plane placement,
transverse-mode content, NSF spatial/temporal resolution and projection work
must be refined for each physical source. A loose allowed work residual is not
an accuracy certificate. No full resolved flute/oboe jet qualification,
nonreflecting-boundary target, thermal-transport matching or perceptual accuracy
claim follows from these source-transaction and coupling tests.

Focused tests are in `fs-lbm/tests/nsf3_characteristic_trial.rs` and
`fs-couple/src/nsf_dg/tests.rs`. They cover exact original-step replay, partial
candidate cancellation, reference-relative work, both face signs, real downstream
feedback, gas-pressure substitution, thermochemical support, projection refusal,
cold area/clock/work admission, and physical-pressure/PCM partition identity.
Pressure-control cases additionally check local characteristic algebra on all
axes, fixed-pressure integrated work, original fluid/DG replay, returned-wave
and material-dependent flow, typed-history callback identity and transactional
domain/cancellation refusal. Native cases exercise source-only composition,
calorics and DG molecular-memory changes through the same pressure-controlled
path. Finite-upstream regressions additionally replay both original ducts and
the original fluid at both face orientations, independently check actuator work,
change upstream termination/length, verify shared-budget refusal and retry, and
exercise the canonical two-sided source through the exterior microphone.

Additional regressions live in `fs-lbm/tests/nsf3_open_boundaries.rs`,
`fs-couple/src/nsf_dg/performance/clock.rs`,
`fs-couple/src/nsf_dg/performance/file/tests.rs`,
`fs-couple/src/dg_port/radiation/tests.rs` and
`fs-couple/tests/music_render_wind/resolved.rs`. These cover finite-wall corner
selection, original interval overlap, source-backed temperature/composition
substitutions at fixed geometry, exact-step radiation, native schema refusal,
PCM partition identity and provenance clocks. These are numerical/interface
regressions; they do not establish measured instrument fidelity.
