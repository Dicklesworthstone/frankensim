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

### Localized source and finite upstream aperture

`NsfDg::new_with_source_patch` selects one original patch on the
source face, opposite the main DG bore. The prescribed boundary history then
acts only on that aperture. Its surrounding cells retain their independently
declared pressure or characteristic, while the complete downstream face feeds
the main bore. The source opening is no longer forced to fill the whole zone.

Attaching `with_upstream_load` to this explicit source-patch owner matches the
finite upstream duct to the aperture area. Both acoustic pressure equations
are solved on the same original fluid interval. Upstream flow, source work and
projection checks use only the selected source quadrature; escaping flow in
the remainder is never sent into that duct. `NsfDgFrame::source_patch` retains
this exchange separately from `coupled_patch`. Both are already included in
their complete face accounts. Use `NsfDg::new_on_patches` when the grid declares
both a source patch and a bore patch; neither acoustic load then absorbs the
other cells on its face.

The source-aperture roots test the accepted incoming controls and the preceding
conditional main-face root as initial guesses. Every trial still consumes the
same original shared budget and must meet the same pressure tolerance. Existing
complete-face and downstream-aperture owners retain their original cold-start
trial sequence.

Native `source_patch` files use the existing `music_render resolved` and
`ensemble --resolved` entry points. The prescribed-source file regression
renders its complete 32-sample history and checks the actual pressure and fluid
state across callback partitions, including the independently driven remainder.

The source-aperture roots test accepted incoming controls and the preceding
conditional main-face root as initial guesses. Prescribed controls carry no
accepted incoming wave and use the bracket midpoint. Every trial consumes the
original shared budget and must meet the original pressure tolerance. Existing
complete-face and downstream-aperture owners retain their cold-start sequence.

The local linear characteristic boundary uses the declared reference normal
velocity to select the entropy and tangential-velocity donors for pressure and
incoming-wave controls. Linearization of the convected equations gives
`ds'/dt + U_reference ds'/dn = 0` and the same convective speed for tangential
velocity. At the zero reference velocity required by `NsfDg`, those modes are
stationary and retain the adjacent fluid state. An acoustic velocity reversal
therefore does not abruptly replace the temperature and shear data supplied to
the original populations, pressure force and Fourier heat flux. The acoustic
pressure/velocity reconstruction, material law, conservative update and root
tolerances are unchanged.

The closed-upstream regression retains its original 32-sample, 2 Pa history,
geometry and numerical budgets, and requires complete playback with identical
accepted states across callback sizes. The earlier instantaneous-donor rule
had refused after 17 intervals. The longer signed `+20/-20/0 Pa` radiation
history has not been rerun with this correction; complete 2,048-sample playback
and a general backflow guarantee remain unclaimed.

Explicit `PressureTemperature` controls still use the actual flow direction
to impose the declared incoming temperature and tangential reservoir data.
Such thermal-source reversal can retain a discontinuity; this correction does
not add material transport into the acoustic DG duct. Any unresolved pressure
gap still reports `NsfDgError::Resolution` with the actual adjacent incoming
bracket and opposite nonzero residuals. A refused interval changes neither
physical owner nor its boundary history; playback retains its accepted prefix
and becomes poisoned.

This enables a split boundary topology, not a complete flute mouth. It does not
create a labium wall at the patch edge, resolve an unmeshed jet, or remove the
original first-order plane-normal characteristic approximation.
The complement remains an ideal prescribed boundary;
escaping flow has no separate exterior-radiation observer or mutual radiation
reaction. Current pressure playback still observes the main DG bore. Geometry,
source-flow, interface and acoustic convergence remain required.

### Transverse reservoir opening

An original `NsfReservoirAperture` can also open explicit tiles in a stationary
wall on a different axis. It coexists with either the complete-face bore or a
selected bore patch. Its density, Cartesian velocity and temperature can be
explicit physical inputs; the original NSF material law supplies pressure,
calorics and transport. It is a prescribed reservoir with the same composition law as the
resolved volume, not an independently mixed atmosphere or transparent acoustic
termination. The canonical `examples/resolved-dg-side-aperture.performance`
and matching Instrument Lab starter expose this original owner through native
playback. `reservoir_aperture AXIS low|high reference COUNT` binds a quiescent
opening to the original source's resolved density and temperature during native
assembly. Pressure, temperature or thermochemical composition edits then update
the opening through the same material owner, with geometry retained. The Lab's
side-opening, slit/edge and joint-radiation fixtures use this binding.

The explicit `DENSITY_KG_M3 VX_M_S VY_M_S VZ_M_S TEMPERATURE_K COUNT` form
remains available for an independently prescribed reservoir; gas edits do not
reinterpret those physical values. `reference` selects the initial state and
does not change subsequent absolute `boundary reservoir` pressure/temperature
controls. Joint radiation still requires the source and computed exterior to
share their physical gas state; it solves pressure through the existing owner.

Tiles must be unique fluid cells on the declared transverse face, clear of the
two-cell characteristic collars, other domain edges and solids, with the
original one-cell fluid halo. The backend retains its symmetric link ownership
at aperture/wall edges, characteristic feedback and conservative energy flux.
All these admission rules remain active in the native reader.

`NsfDgFrame::reservoir_mass_out_kg` reports the original streamed mass leaving
all transverse openings, separately from both characteristic faces. The complete fluid
mass balance uses the two face fluxes **plus** this reservoir flux. Signed
reservoir energy entering the fluid remains `flow.open_energy_j`; it must be
included with `flow.open_boundary_energy_j`, wall heat and source work when
checking total fluid energy. Mixed-boundary unresolved internal exchange retains
the original ledger classification and is not reported as resolved viscous heat.
The reservoir flux is not inserted into either acoustic duct. Its separate
continuum volume flow can be observed at the shared exterior microphone with
the explicit `reservoir_observation` selection described below.

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

### Independent breath and escape openings

The resolved grid accepts up to 16 independently indexed transverse openings
through `with_reservoir_apertures`. Native data declares their count before the
existing cell-list or SI-rectangle blocks:

```text
reservoir_apertures 2
reservoir_aperture 1 low reference 1
reservoir_cell 5 0 1
reservoir_aperture 1 high reference 1
reservoir_cell 5 2 1
```

Indices are zero-based in authored order. Every opening retains the original
geometry/halo checks, ghost state, streamed mass, conservative FV total energy
and acoustic pressure work. Duplicate tiles refuse; diagonal links crossing an
ownership seam remain walls. Combined tile storage is bounded by the original
grid cell budget. The legacy single-aperture block retains its original meaning.

`boundary aperture 0 COUNT` drives the selected breath opening with the same
absolute-pressure/temperature hold records as `boundary reservoir`. It can
coexist with both the bore and a finite upstream/cork-side DG load, with zero
added axial series effort. Warm breath uses the original gas energy and transport
laws; the separate escape opening can remain bound to the exterior reference.

`reservoir_observation aperture 1 baffled OX OY ANGLE_RAD` selects that
opening's original escape tiles for the microphone. Repeat this record for
every explicitly observed aperture. Each original aperture may occur once;
its independently placed tiles must be disjoint from other observed surfaces
and the bore outlets. Unselected breath openings do not acquire a receiver.

With several observations, `reservoir_feedback aperture INDEX ...` selects
exactly one observed opening for the existing live exterior reaction. The
unindexed feedback form remains valid when only one opening is observed.
The source and passive reaction indices must differ. All held breath and trial escape states enter the same original
NSF transaction with both acoustic loads. No participant publishes before the
existing pressure/work gates pass.

Typed callers use `step_with_reservoir_aperture`,
`with_reservoir_radiation_at` and `with_reservoir_observation_at`. Singular
selection methods require exactly one declared opening. Observation accessors
have indexed forms (`reservoir_observation_at`, `reservoir_receiver_at`,
`reservoir_observer_gas_at`) and `reservoir_observation_count`; their singular
compatibility views return no value when several openings are observed. Accepted
`NsfDgFrame::reservoir_exchanges` retains each owner's original interval account;
`mass_out_kg` and `total_energy_out_j` are positive outward. Their sums recover
the aggregate streamed mass and the negative of `flow.open_energy_j`.
The reaction uses the selected escape's account, never the sum with breath.

`examples/resolved-dg-separate-openings.performance` exercises these declarations
with a warm transverse source, distinct radiating escape, downstream bore,
finite closed upstream duct and common material microphone. It is a numerical
integration fixture, not measured flute geometry or a sound-quality claim.

### Explicit subgrid transport through open boundaries

The same NSF volume can carry an explicitly declared Smagorinsky, WALE or
Vreman exploration closure with characteristic faces, a selected bore patch
and a transverse reservoir aperture. Missing open neighbors use the original
reconstructed characteristic face state or prescribed reservoir ghost state,
at their actual half-cell or full-cell spacing. They are not treated as
impermeable walls. The local modeled shear transport reaches the original
collision and stress-work calculation, and modeled conductivity reaches the
same conservative open-face heat flux. Boundary-control rebinding and rejected
coupled trials preserve these fields together with the original fluid state.

Every coefficient, filter width, viscosity-ratio bound and turbulent Prandtl
number is explicit. Bounds refuse excess transport; they never clamp it.
Omitting the declaration keeps molecular transport. A zero coefficient keeps
the identical molecular trajectory and remains visible in native input
metadata as a disabled control. Complete prescribed reservoir faces retain
their separate unsupported-SGS refusal; this integration admits the localized
transverse aperture used by the resolved source.

These algebraic closures supply no fluctuating noise source, SGS kinetic-energy
equation or transition prediction. They remain exploration tools: the native
reader requires the `exploration` token and the source remains ineligible for
acceptance renders. Filter, geometry, source-flow and interface convergence are
still required; enabling a closure does not qualify a turbulent flute jet.

### Manufactured slit, free gap and solid edge

`examples/resolved-dg-slit-edge.performance` puts a solid nozzle and edge into
this same native path. Its 12 × 6 × 2.5 mm Cartesian domain uses 0.5 mm cells,
with 375 explicit stationary solid cells forming a 3 mm long, 1 mm high slit
and a staircased wedge. The slit exits at x = 4 mm; the edge starts at x = 7 mm.
A 30-tile lower downstream patch feeds the original 60 mm DG tube, with both
areas equal to 7.5 × 10⁻⁶ m². A separate 27-tile top opening retains the original
ambient reservoir. All collars and reservoir fluid halos remain explicit.

The declared 0.02 → 0.2 → 2 Pa pressure steps and release drive the actual fluid
through those solids. Four original NSF steps form each coupled interval;
96 kHz observation and optional PCM decimation retain that physical clock.
The native computed baffled load, Stokes–Kirchhoff observer and fixed 20 Pa
PCM scale are unchanged. Instrument Lab offers the same canonical file as
**Resolved flow — slit jet and solid edge**; it can also be passed directly to
`music_render resolved` with the existing command options.

This is a manufactured integration geometry. The available Bressan recorder
description does not supply its complete windway height, curvature and jet
offset, so those dimensions are not presented as measured instrument data.
The two-cell slit, voxel edge, downstream fluid collar and same-area circular
plane-mode connection require refinement and physical validation. The example
contains no jet card or synthetic noise source; it does not establish
self-sustained edge tone, a converged turbulent jet or realistic flute sound.

### Time-varying transverse breath conditions

`NsfDgDrive::Reservoir { pressure_pa, temperature_k }` and
`NsfDg::step_with_reservoir` change the held **absolute** pressure and temperature
of an attached transverse aperture. The original gas source derives density,
calorics and transport at those conditions; the aperture retains its declared
Cartesian ghost velocity. Pressure is not translated into an authored jet-speed
curve. The original population, pressure-force and conservative energy owners
determine the response.

The main bore and an optional finite upstream duct remain coupled. With both
attached, their axial interfaces solve together with zero additional series
pressure rise. This permits a side source between a bore and a closed cavity.
With only the main load, the other face receives zero incoming characteristic.
Changing reservoir conditions rebases the force and flux history while retaining
the fluid's stored physical state. Conditions are staged on each private trial;
failed root solves, source refusals and cancellation preserve the accepted ghost,
fluid and both acoustic histories. The normal reservoir mass and energy accounts
include the resulting exchange.

`examples/resolved-dg-transverse-drive.performance` combines a changing side
pressure/temperature, a closed upstream cavity, the bore's computed radiation,
and the original per-tile reservoir observer at the same microphone. The declared
common baffle places the reservoir tiles with origin `[-0.05, 0]` m and zero
rotation; the selected gas-propagation law reaches both sets of paths. This is
one-way observation and a coarse synthetic integration example, not a measured
mouth or a qualified sounding flute. The reservoir has no finite
mass/heat inventory and retains the grid's original composition; changing its
temperature does not add a new species-transport or anatomical airway model.

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

`NsfBoundaryDrive::PressureTemperature { pressure_offset_pa,
inflow_temperature_k }` adds a physical inflow temperature at that same source
face. Pressure and the live returning characteristic still determine flow.
Where flow enters the resolved volume, the face temperature is the supplied
value in kelvin. Where it leaves, the adjacent fluid retains its original
outgoing entropy. The original material source supplies density, caloric state,
transport and conservative energy flux; composition remains fixed. The held
temperature is admitted even during outflow, and failed changes or trials
preserve the previous controls, fields and histories.

The native selection is `boundary pressure-temperature COUNT`, followed by
`hold END_OBSERVATION_SAMPLE PRESSURE_OFFSET_PA INFLOW_TEMPERATURE_K` rows.
Both controls follow the original physical interval clock. An independently
controlled patch remainder accepts `pressure-temperature PA K` as well.
The manufactured `crates/fs-couple/examples/resolved-dg-warm-source.performance`
demonstrates a warm source through the original constriction, live DG bore and
material exterior observer. Its temperatures are numerical inputs, not measured
player data. This feature does not add humidity/species transport, finite airway
thermal capacity or mean enthalpy transport into the acoustic DG bore.

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

### Physical breath ramps

Native histories can vary their physical source controls with linear ramps:

```text
boundary pressure 3
ramp 480 0 20
hold 1440 20
ramp 1920 20 0
```

A scalar ramp is `ramp END_OBSERVATION_SAMPLE START_VALUE END_VALUE`.
Values retain the selected boundary's units: pressure offsets and upstream
series pressure rises in Pa, incoming characteristics in m/s. For
`pressure-temperature`, `reservoir` and indexed `aperture` controls, each
endpoint is a pressure/temperature pair:

```text
boundary pressure-temperature 2
ramp 960 0 300 20 303
ramp 1920 20 303 0 300
```

Each segment begins at the previous segment's exclusive endpoint, or zero for
the first segment. Its final value is the value approached at its exclusive
end; the following segment owns that instant. Both endpoint temperatures must
lie in the original fluid temperature window even if the final endpoint is
never sampled. Reservoir pressure remains absolute; pressure-temperature
source pressure remains relative to the original face reference. An indexed
`boundary aperture INDEX COUNT` retains that same physical opening throughout
every ramp; it never interpolates or switches the opening identity.

The reader evaluates each ramp at the **start of the original physical
NSF/DG interval**, measured on the observation clock, then retains that control
through the interval. It does not change the source clock or interpolate PCM.
The existing exact-clock path retains integer sample positions; different
observation clocks use the same physical interval duration as the source.
As with holds, physical control discretization requires interval refinement.
Hold-only input preserves its existing arithmetic.

`crates/fs-couple/examples/resolved-dg-warm-ramp.performance` exercises warm
pressure ramps through the original resolved volume, live bore and exterior
observer. It is a synthetic integration source, not a measured player envelope
or instrument-fidelity result.

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
[optional sgs exploration MODEL COEFFICIENT FILTER_WIDTH_LU MAX_EDDY_TO_MOLECULAR TURBULENT_PRANDTL]
characteristic AXIS low|high MAX_RELATIVE_PERTURBATION
[optional coupled_patch START0 END0 START1 END1 pressure REMAINDER_OFFSET_PA]
[optional source_patch START0 END0 START1 END1 pressure REMAINDER_OFFSET_PA]
[optional reservoir_apertures APERTURE_COUNT; repeat the following aperture block APERTURE_COUNT times]
[optional reservoir_aperture AXIS low|high reference COUNT]
[or explicit reservoir_aperture AXIS low|high DENSITY_KG_M3 VX_M_S VY_M_S VZ_M_S TEMPERATURE_K COUNT]
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
[repeat exactly COUNT increasing hold or ramp records]
end
```

Alternatives use the same original parsers and owners:

- `sgs exploration smagorinsky|wale|vreman C DELTA RATIO PR_T` attaches the
  original `NsfSgsConfig` after the declared openings. `C` is the model's own
  coefficient (Vreman uses `c`, not its square root); `DELTA` is in original
  lattice cells, so its physical width is `DELTA * DX_M`. `RATIO` bounds modeled
  over molecular viscosity and `PR_T` controls modeled heat diffusion. No
  values are inferred from the instrument name or desired sound.
- `coupled_patch START0 END0 START1 END1 pressure PA` selects an explicit bore
  aperture immediately after `characteristic`; `characteristic M_S` can replace
  `pressure PA` for an incoming-wave condition on its remainder. Cell bounds are
  half-open, nonempty, and inside the original transverse grid. The selected
  aperture area must match the main DG inlet. Omitting the record preserves the
  existing complete-face path. The canonical
  `examples/resolved-dg-split-aperture.performance` and Instrument Lab starter
  exercise the native aperture, surrounding ambient condition and bore-outlet
  microphone together. This coarse manufactured fixture is not flute geometry.
- `source_patch START0 END0 START1 END1 pressure PA` selects an independent
  original patch on the opposite source face. The `boundary` history drives
  those cells; the remainder retains `pressure PA` or `characteristic M_S`.
  The main DG inlet must match the selected bore patch, or its complete face
  when no bore patch is declared. If `upstream_dg` is present, its inlet must
  match the selected source area. One source and one coupled declaration may
  appear in either order; duplicate roles refuse, including mixed indexed/SI
  declarations. Neither record changes the original grid.
- `gas mixture T_K P_PA mole|mass COUNT SPECIES FRACTION ...` and the shared
  `gas pack NAME` records retain the original named species and complete
  composition in `NamedMixtureNsfGas`. The original material owner evaluates
  every evolving cell and thermal face: mixture EOS, molecular viscosity,
  conductivity, sound speed and `cp/R` for SGS heat remain coupled. Its named
  near-ambient law has constant source heat capacities, so conserved sensible
  energy and temperature reconstruction use exactly `(cp - R) T`; no fitted
  single-species `GasSpec` or NASA9 approximation is substituted. One canonical
  source state binds the fluid, reference reservoirs, main duct and optional
  upstream duct. The native record retains the shared reader's canonical mole
  fractions, including for mass-fraction declarations.
  The original 273.15–323.15 K, 10–200 kPa and water-saturation limits apply to
  the reference, declared temperature window and actual evolved/controlled
  states. A saturated reference cannot authorize colder or supersaturated
  states. Composition is frozen: this does not add species advection,
  condensation, chemistry or finite-inventory thermal transport into the DG
  bore. The NSF bulk/shear viscosity ratio remains an explicit source input.
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
  `dg_terminal_closed PORT_ID` on a finite `dg_side_chimney` makes its distal
  end an actual closed DG boundary. Its air column, material and wall memory
  remain in the coupled bore, harmonic response and continuous renderer, while
  its capped endpoint has no exterior-radiation row or microphone flow. The
  shared 32-side-port bound counts internal sites independently of exterior
  openings. A flush-site closure retains the existing terminal-channel law.
- `dg_molecular` selects the existing causal molecular-relaxation law for the
  original DG gas region, using exactly the same grammar as native reed/flue
  sources. Its species composition comes from the original named-mixture or
  thermochemical source record, or the axial field. A dry-air scalar snapshot has no species
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
- `boundary pressure-temperature COUNT` adds an inflow temperature in kelvin
  after each pressure offset. Temperatures outside the declared fluid window,
  missing temperature fields and extra fields refuse. The original material
  and characteristic-domain checks apply to every physical reconstruction.
- `boundary reservoir COUNT` selects the existing transverse aperture instead
  of an axial actuator. Each record is `hold END_OBSERVATION_SAMPLE P_ABSOLUTE_PA
  TEMPERATURE_K`, with positive pressure and temperature in the original gas
  domain. It requires exactly one aperture and may coexist with `upstream_dg`;
  neither axial face then receives an added series-pressure source.
- `boundary aperture INDEX COUNT` uses the same held absolute pressure and
  temperature on one explicitly indexed opening. Unknown indices refuse; a
  selected passive reaction opening cannot also be the prescribed breath source.
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
selects the face coupled to the DG inlet. The `hold` or `ramp` history drives the opposite
face, its explicit upstream series actuator, or the explicitly selected
transverse reservoir. Endpoints are exclusive indices on the **observation** clock and must
cover the declared window exactly. Each physical interval uses the segment's value
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

### Resolved escape flow at the same exterior microphone

`reservoir_observation baffled OX OY ANGLE_RAD`, immediately after `observation`
and before `observer_gas`, includes the original transverse reservoir aperture
in the same emitted pressure. `NsfDgPerformance::with_reservoir_observation`
is the equivalent typed selection for a single declared opening. For multiple
openings, repeat `reservoir_observation aperture INDEX baffled OX OY ANGLE_RAD`
or call `with_reservoir_observation_at` once per selected identity, before
`observer_gas`. The pressure sum follows the declared observation order.
Every actual fluid face tile has its own
signed volume-flow input; the receiver does not substitute bore flow, infer a
uniform mouth piston, or divide compressible mass flux by reference density.
The source owner averages fluid/ghost velocity at that physical face and uses
the original endpoint trapezoid and all physical substeps. The streamed mass
and energy accounts remain separate. Failed trials retain all accepted fields
and tile-flow accounts.

Placement explicitly declares a common infinite baffle with the existing bore
openings. For reservoir normal axis `a`, the local coordinates use grid axes
`(a+1)%3` and `(a+2)%3`; the second coordinate reverses on a low face, preserving
the outward normal. Grid edge zero maps to `[OX, OY]`, then the supplied angle
rotates the actual tiles in the plane. Tile areas and spacing remain `dx²` and
`dx`; neither aperture size nor surface flow is independently authored. The
same microphone observes both surfaces with their separate travel times and
signs. Overlapping openings, underresolved tile phase, and nonquiescent or
mismatched observed-opening/exterior reference gases refuse. A distinct breath
opening may have its own prescribed temperature.

The shared `BaffledTilePressure` adapter partitions larger surfaces into original
receivers with at most 128 inputs each, retaining every tile's signed flow, area,
position and physical delay. It never averages tiles into one piston. The
explicit reference allowance admits at most 120,000 tiles, with aggregate
retained history and kernel limits of 256 MiB each; material preparation budgets
also apply across all selected surfaces, together with the bore's material
paths. Adding an aperture cannot multiply the shared tile, storage or gas-filter
allowance. Reaching a budget refuses construction.
Larger banks reuse the original linear pressure differentiator per piece.
All pieces preview before any observer history accepts a sample.

The shared differentiating Rayleigh kernel observes every tile after the same
conservative physical-to-output clock projection as the bore. A selected
`observer_gas` law reaches every selected path, with only numerical padding
aligned once to their common maximum before the pressure sum. Physical flight
times and gas phase remain distinct. The native sidecar identifies original
cell coordinates, placement, flow convention, path material and latency;
multiple selections appear in `reservoir_observations`, while a single
selection retains the existing `reservoir_observation` record.
`resolved-dg-side-aperture.performance` selects a single escape;
`resolved-dg-multiple-escapes.performance` observes two independently placed
openings with different original tile counts. These are numerical composition
fixtures, not measured instrument geometry.

Without the optional `reservoir_feedback` row below, this is one-way
observation of a prescribed-reservoir opening. Observation alone adds no
reaction to the fluid, no transparent boundary, no acoustic/hydrodynamic
decomposition, no free-field body scattering, and no separate observation of
the complementary region of a characteristic patch. Source and boundary
placement/refinement remain necessary; this is not measured flute validation.

### Joint reaction of the resolved mouth and bore outlets

After the explicit reservoir observation placement, an optional native row
selects the same-interval exterior feedback:

```text
reservoir_feedback PRESSURE_LO_PA PRESSURE_HI_PA QUAD_NX QUAD_NY
```

When more than one opening is observed, name the original reaction owner:

```text
reservoir_feedback aperture INDEX PRESSURE_LO_PA PRESSURE_HI_PA QUAD_NX QUAD_NY
```

This remains one original joint mouth/bore reaction. Other selected escapes
contribute their prescribed-boundary flows to observation without acquiring a
second reaction law or mutual coupling between all reservoirs.

The pressure bounds are ghost-pressure offsets relative to the original source
reference, in Pa. The two positive quadrature counts refine the existing
rectangular radiation source independently of the fluid grid. This selection
requires zero placement rotation and a complete rectangular set of actual
reservoir cells. Its dimensions come from those cells; the reader never
rescales an aperture to match an acoustic section. The original full-band
surface, fit, causality, state and passivity admissions still apply.

The source may use the original axial characteristic, pressure or
pressure/temperature control, or the finite upstream duct's pressure rise.
`boundary reservoir` cannot prescribe the same transverse opening that the
passive exterior must solve; this conflicting pair refuses before construction.
The typed `step_with_reservoir` likewise refuses without advancing any owner.
For separately declared boundaries, `boundary aperture INDEX COUNT` and
`step_with_reservoir_aperture` can prescribe breath while the distinct observed
escape receives this reaction, including with a finite upstream duct attached.

The typed equivalent is `NsfDg::with_reservoir_radiation` with an explicit
`NsfDgReservoirRadiation`. The original joint radiation model receives two
independent coordinates: flow into the bore and flow out through the reservoir.
The same matrix retains every mouth/outlet and outlet/outlet interaction.
No equality between bore flow and escaping flow is imposed on a resolved
volume that can store and redirect gas.

At each physical interval, bounded nested solves match both the bore face and
the reservoir's actual mean acoustic pressure to the original coupled acoustic
response. An attached finite upstream duct adds its existing pressure-rise
equation to this same transaction. Every numerical candidate starts from the
accepted physical state; cancellation, unbracketed roots, exhausted total work
budgets and excessive pressure/work errors publish none of the participants.
The accepted frame exposes the reservoir pressure, independent flow, radiation
work and signed projection discrepancy. Its `resolved_exchange` retains the
complete original fluid interval, including every NSF substep. Pressure matching
and acoustic work admission both use the cached FV predictor quadrature that
updates conservative total energy. This work is the spatial average of
`p_cell * u_cell` and `p_ghost * u_ghost`, integrated with the original predictor
time trapezoid. It is not a product of independently averaged pressure and flow,
nor a reconstruction from the corrected endpoint state.

The per-tile endpoint pressure/work traces remain available for observation
diagnostics. They cannot replace the conservative exchange in the reaction
equation. Original streamed mass and conservative total-energy flux remain
separate; the acoustic projection does not absorb enthalpy, conduction or
viscous work.

The microphone still observes each original reservoir tile and each circular
bore outlet once, at their actual positions, with the selected shared gas law.
The radiation load uses a uniform mouth-pressure projection. Differences from
the spatially and temporally resolved pressure work remain measured covariance;
they are not heat, manufactured damping or an extra acoustic source. This
optional reaction does not establish a transparent reservoir boundary, resolved
external hydrodynamics, free-field body scattering or measured flute accuracy.
Physical interface, geometry, quadrature and clock refinement remain required.

`resolved-dg-radiating-mouth.performance` is a native integration fixture that
selects the joint feedback and material observation paths. It retains twelve
original NSF steps per 8 kHz acoustic interval and is a numerical example,
not a measured instrument or a realism preset. The Instrument Lab includes it
as a starter with the same original data. Its reference-gas aperture follows
source pressure/composition changes without a separately reauthored air density.

The existing `resolved-dg-slit-edge.performance` Instrument Lab source now also
selects this reaction and the original per-tile microphone. Its 27 escape cells
form a 1.5 mm by 4.5 mm rectangle centered at `[-0.05, 0]` m in the shared baffle.
The explicit observation origin is `[-0.05125, -0.00725]` m and the native
`reservoir_feedback -5 5 16 32` row supplies the ghost-pressure search domain
and independent radiation quadrature. The unchanged physical clock requires
the joint radiation source to resolve a band extending to approximately 127 kHz.
The circular surfaces therefore use 48 radial rings and 96 angular points,
with a 128 MiB cold-fitting memory allowance, so they satisfy the original
phase-resolution and pair-distance-storage admissions. These numerical
resources do not relax the fit, passivity or coupling error tolerances.
The original voxel slit/edge, gas, pressure history and physical clock remain
the source of the calculation. The circular receiver observes the bore outlets once;
the tile receiver observes the mouth once. This coarse manufactured fixture
still requires geometry, boundary, clock and acoustic refinement before any
edge-tone, self-oscillation or instrument-fidelity claim.

The shared radiation admission includes all 32 allowed DG side ports, the
axial outlet and one independent mouth (34 coordinates). Every self and mutual
term remains present, with the existing fitting-memory and runtime-state
bounds.

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

## Physical geometry for grid refinement

Native resolved input can keep both solids and openings in metres while the
Cartesian grid is refined. Replace the indexed `solids` block with:

```text
solid_geometry SHAPES MAX_CELL_EVALUATIONS
solid_box add|subtract X0 Y0 Z0 X1 Y1 Z1
solid_wedge add|subtract X0 Y0 Z0 X1 Y1 Z1 NX NY NZ OFFSET_M
```

Rows compose in declaration order, starting with an addition. A wedge is the
original FRep box intersected with its normalized-normal half-space; the offset
is a distance in metres. All bounds must be inside the grid. The complete
shape-count times cell-count work must fit the declared allowance. A primitive
that selects no cell centers refuses, so a thin windway or edge cannot silently
disappear.

The existing source and bore aperture records also accept physical bounds:

```text
source_patch_m A0 A1 B0 B1 pressure PA
coupled_patch_m A0 A1 B0 B1 characteristic M_S
reservoir_aperture_m AXIS low|high RHO VX VY VZ T A0 A1 B0 B1
reservoir_aperture_m AXIS low|high reference A0 A1 B0 B1
```

One characteristic patch per face is allowed. Source and bore patches may
coexist and have different areas. Each remainder can use the same independent
`pressure`, `pressure-temperature`, or `characteristic` controls as the indexed
form. A and B are the two **increasing Cartesian axes other than the normal**:
for normal y they are x,z. Bounds are half-open and select original face-tile
centers. The physical reservoir record generates its original cell list and
does not take subsequent `reservoir_cell` rows. Low/high retains the original
normal orientation and flux sign.

The `reference` form binds the same original resolved gas state as the indexed
reservoir record. Its density and temperature follow the source material data,
while its physical aperture bounds remain unchanged.

The resulting solids use the existing halfway bounce-back walls; openings use
the original face quadrature. Each represented aperture edge is within half a
cell width of its declared position. This is not subcell geometry or a
convergence certificate. Keep the SI geometry fixed, refine the grid and clocks,
and retain the existing collar, gas, cell-work and matching DG inlet-area checks.
`crates/fs-couple/examples/resolved-dg-si-slit-edge.performance` gives the
complete manufactured slit/edge input with physical solids and apertures.

### Steady gas in open bore branches

The shared `DgSourceInput` now accepts `dg_gas_network mole SPECIES_COUNT ID...`
after side/terminal declarations and before `dg_radiation`. Its
`gas_inlet`, `gas_bore`, `gas_cell`, `gas_side` and `gas_solve` records declare
one source-backed stationary mass/species/enthalpy problem. The complete grammar
is documented on `DgSourceInput::read`.

Each side discharge reduces downstream main-bore flow. Finite chimneys own
thermal cells from their actual geometry; flush openings have no invented
thermal volume. The solved regional states feed the original DG viscothermal
and molecular operators. Closed chimneys retain their separately declared gas.
This route uses prescribed steady flow, not an inferred mean acoustic flow,
transient NSF heat/species transport, or a predicted pressure drop.

### Curved SI solids

The same `solid_geometry` block also accepts these original FRep primitives:

```text
solid_cylinder add|subtract X0 Y0 Z0 X1 Y1 Z1 RADIUS_M
solid_sphere add|subtract CX CY CZ RADIUS_M
```

Cylinder endpoints are the centers of its flat caps in metres. The original
infinite-cylinder primitive is intersected with a finite box, then rigidly
rotated and translated to those endpoints. Reversing the endpoints describes
the same solid. Sphere coordinates are its center and radius in metres. A round
passage can be authored by adding an outer cylinder and subtracting a narrower
one with the same endpoints; subsequent shapes retain declaration-order CSG.

The entire physical shape, including its radial extent, must lie in the grid.
Nonpositive radii, coincident cylinder endpoints, unresolved shapes, exhausted
cell-evaluation budgets and cancellation refuse. The projected circular-cap
bounds permit an axial end cap exactly on a grid boundary without permitting
radial overhang. These primitives reach the original stationary NSF mask and
retain its halfway bounce-back and characteristic-collar restrictions; they do
not create cut cells or supply missing measured instrument/player anatomy.

The focused original-owner tests check oblique cylinders, reversed endpoints,
finite caps, hollow-passage subtraction and independent analytic volume changes
on three grid rungs. This is bounded manufactured-geometry evidence, not a
universal convergence or real-flute qualification.

### Independent source and bore openings

`NsfDg::new_on_patches` requires one original characteristic patch on each
axial face. Their rectangles, areas, aperture drives and surrounding boundary
laws remain independent. The existing prescribed-source, finite-upstream and
reciprocal reservoir-radiation paths all use these same two selected traces.
Every selected duct still requires exact source-gas, area and clock admission.

The original grid exposes `characteristic_patch_on(high)`,
`characteristic_patch_drive_on(high)`, `open_patch_traces_on(high)` and
`last_patch_exchange_on(high)`; trials expose `patch_exchange_on(high)`.
`set_boundary_drives_with_patches` and the two plural preview methods carry
one optional patch control per physical low/high face. An omitted control
retains that face's current patch drive; a supplied control requires a patch.
The singular accessors return `None` with two patches, and singular operations
that would ambiguously change one patch refuse.

The original coupled root, NSF update and DG propagators are retained. Each
patch's flow and pressure-work enter only its own acoustic projection. Full-face
mass and total-energy accounts already include both patches and their separate
remainders; do not add selected-patch accounts to those totals again. All owners
and controls publish together after the original pressure, work and cancellation
checks. The characteristic-collar, projection, explicit thermal-backflow and
frozen-composition limits still apply. Two apertures do not establish resolved
lip anatomy, full bore heat/species transport or qualified flute sound.

## Finite thermochemical source chamber

A selected transverse opening can now exchange mass and complete total energy
with a finite well-mixed chamber. Its volume is fixed unless an explicit
thermochemical actuation drive supplies volume rate and heat power. Attach it with
`NsfDg::with_finite_reservoir(aperture, volume_m3, limits, gate)`, then advance
with `NsfDgDrive::FiniteReservoir { aperture }`. The initial inventory comes
from that opening's original density and temperature and the **same frozen
named-mixture or thermochemical source law** owned by the NSF grid. Named gas
uses its retained species/composition and original constant heat capacities;
NASA9 keeps its original variable-capacity law and bounded caloric inverse.
The legacy scalar `GasSpec` medium still refuses this source-owned inventory.

The chamber's accepted pressure and temperature drive the original reservoir
ghost during the next complete NSF/DG interval. Actual streamed mass and the
original conservative FV total-energy exchange update its retained inventory.
Thermochemical inventory remains the original `ClosedGasParcel`. Named inventory
retains `NamedMixtureNsfGas`, fixed volume and an actual source-resolved state;
its sensible energy is exactly `m (cp - R) T`, with zero formation offset. Both
exchange directions and nonzero energy with zero net mass use that same law.
The named temperature increment is `(dE - dm * u_old) / (m_new * cv)`, followed
by the original mixture's actual-temperature/pressure evaluation. No discharge
coefficient, authored jet speed, fitted gas or replacement EOS enters this connection. Incoming kinetic energy thermalizes in the
quiescent well-mixed chamber. Its rigid support carries the unmodeled momentum
reaction; this owner provides mass and energy storage, not resolved chamber
momentum or player anatomy.

### Native records

After all reservoir aperture records and before `initial`, add:

```text
finite_reservoir APERTURE VOLUME_M3 MAX_REL_MASS MAX_REL_ENERGY MASS_ABS_KG ENERGY_ABS_J REL_TOL
```

The two movement fractions must be strictly between zero and one; absolute
mass and energy balance tolerances must be positive; the relative balance
tolerance is in `(0, 1e-3]`. The selected opening must have zero imposed ghost
velocity. Its authored density and temperature initialize the source, including
its initial pressure through the original material EOS.

The boundary history names the same opening and holds only its ownership:

```text
boundary finite-reservoir 0 1
hold 64
end
```

The endpoint is on the existing observation clock. No pressure or temperature
values follow this passive hold. Prescribed source controls, a different
selected opening, or a ramp of passive ownership conflict with this declaration
and refuse. The separate actuation history below supplies physical volume rate
and heat power for a thermochemical chamber. Existing pressure/temperature/
characteristic holds and ramps retain their original grammar for sources
without a finite chamber.

A separate opening can retain `reservoir_feedback` and its material microphone.
The source also composes with the original downstream duct, optional finite
upstream/cork duct, and independently selected axial patches. The finite source
and passive radiation cannot own the same opening, in either attachment order.
Chamber reconstruction and conservation checks occur before publication at all
three fluid/acoustic acceptance sites. Failed inversion, insufficient inventory,
unsupported source states, failed balances or cancellation leave every physical
owner at the preceding accepted interval.

### Energy chart and interval limits

The original NSF total energy uses a fixed sensible chart,
`u_s(T) = u_NASA(T) - u_offset`. With `dm` and `dE_s` positive **out of NSF
and into the chamber**, the original-law inversion uses:

```text
m_next = m_old + dm
rho_next = m_next / chamber_volume
du = (dE_s - dm * u_s(T_old)) / m_next
```

The reconstructed parcel's actual mass is audited against `dm`. Its energy
audit uses the original cancellation-safe NASA internal-energy increment:

```text
stored_change = m_actual_next * delta_u_NASA
              + (m_actual_next - m_old) * u_s(T_old)
energy_residual = stored_change - dE_s
```

The full FV energy already contains pressure work, advection, conduction and
viscous work. Adding pressure work or acoustic work again would debit it twice.
The streamed mass is not replaced by the distinct continuum volume-flow
quadrature. Absolute tolerances and change-scaled relative tolerances bound
both accepted residuals; absolute formation energy never scales those budgets.

Each original substep records the absolute mass and energy transferred through
each aperture. Interval trials sum these variations before opposite signs can
cancel. The chamber bounds those sums relative to its starting inventory and
checks all corners of the resulting mass/energy envelope against the original
caloric temperature bracket. A short withdrawal followed by reverse flow
therefore cannot hide depletion or an unsupported temperature behind a small
net exchange. Nonzero increments that cannot change the represented inventory
refuse instead of silently disappearing. Zero net mass and energy preserve an
exact inventory clone after the movement checks.

For a named source, the outward-rounded mass/energy envelope also bounds
all possible prefix temperatures and pressures. Its most conservative
low-temperature/high-pressure corner is checked through the original mixture
owner, preserving water-saturation and material-pressure refusals even when
the final net exchange returns toward the initial state. Both envelope and
accepted endpoint must remain within the receiving NSF grid's authored
temperature window. The reconstructed endpoint also passes the original NSF
reservoir-ghost admission, including its local viscosity and temperature-change
limits, before it can become the next interval's source. These checks precede
publication; a refused source cannot leave the fluid, either acoustic duct or
its accepted inventory advanced.

Library callers inspect `finite_reservoir_inventory()`, `finite_reservoir_limits()`,
and `NsfDgFrame::finite_reservoir` for accepted mass, pressure, temperature,
source exchange, sensible storage change and conservation residuals. The common
inventory snapshot also identifies its named-mixture or thermochemical source
and retains the sensible energy and energy reference. The existing
`finite_reservoir()` accessor remains the NASA9 `ClosedGasParcel` view and
returns `None` for a named inventory; use the common snapshot or aperture getter
when testing whether any finite source is attached. Native metadata reports both
source kinds through the common snapshot, with the chamber volume, limits and
partition.

This is an **explicit first-order source partition**: pressure and temperature
are held over each coupled physical interval, then updated. Refine that interval
and the grid for source-feedback accuracy. The acoustic DG ducts still omit
mean mass/enthalpy advection; the chamber does not fill that separate gap.
The optional actuation below supplies imposed volume work and heat. It does
not solve player anatomy, moving-wall inertia, source refill, composition
evolution or experimental flute qualification.
The focused integration and native-input regressions are authored numerical
checks; the direct GitHub implementation session did not execute native builds,
runtime tests, convergence campaigns or listening comparisons.

### Supplied-source chamber example

`crates/fs-couple/examples/resolved-dg-finite-reservoir.performance` connects
a 10 mL nitrogen chamber to opening 0, a separate radiating escape to opening 1,
the original downstream bore, a finite closed cork-side duct, and the existing
material microphone. It retains the previous separate-opening fixture's physical
mesh while the original nitrogen properties determine its physical timestep.
The 384 samples at 48 kHz are an 8 ms numerical integration window, not a
playable phrase or a measured flute model. The 48 kHz observation clock is
accepted by the native PCM command without `--decimate`; it does not change
the original physical interval or the twelve NSF substeps per coupled interval.

The resolved gas starts at 300 K and 95000 Pa. The quiescent chamber starts at
300.15 K with density `1.0664144914103022 kg/m3`, corresponding to approximately
95002 Pa under the supplied N2 species EOS. That density follows the original
species molar mass `28.01340 g/mol`; the transport correlation's distinct fit
molar mass is not substituted into thermodynamics. The 280..340 K caloric window
and 95000 Pa reference stay within the supplied dilute transport source's
declared domain. Source admission still applies to every actual state.

Build the three existing source packs and render from the repository root.
Use fresh output paths for repeated runs:

```sh
mkdir -p target/resolved-chamber-packs
for source in nitrogen nitrogen-nasa9-mcbride2002 nitrogen-transport-lemmon2004
do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$source/manifest.tsv" \
    --out "target/resolved-chamber-packs/$source.pack"
done
cargo run -p fs-couple --bin music_render -- resolved \
  crates/fs-couple/examples/resolved-dg-finite-reservoir.performance \
  target/resolved-dg-finite-reservoir.wav --block 8 \
  --gas-component target/resolved-chamber-packs/nitrogen.pack \
    target/resolved-chamber-packs/nitrogen-nasa9-mcbride2002.pack \
    target/resolved-chamber-packs/nitrogen-transport-lemmon2004.pack \
  --gas-water none --gas-phase frozen
```

Instrument Lab includes the same file verbatim as **Resolved flow — finite
breath chamber**. Attach the same three packs with the existing gas-source
controls. The browser still delegates the solve to the original native renderer.
Changing `finite_reservoir` volume changes stored gas capacity; changing the
source aperture's density or temperature changes its original initial state.
Neither edit supplies a desired pitch, artificial turbulence or calibrated
audio. The example and new native regressions require an actual native run
before stability, numerical accuracy or audibility claims.

### Physical chamber motion and heat

`NsfDgDrive::FiniteReservoirActuation` adds an explicit volume rate in
`m3/s`, heat power in `W`, and a maximum relative volume change per
physical interval. Positive volume rate expands the chamber; negative rate
compresses it. Positive heat power supplies energy to the gas. These commands
are physical boundary inputs: the original gas law determines pressure and
temperature, and the original resolved opening determines mass flow.

This attachment requires the retained NASA9 thermochemical parcel. Passive
finite sources retain both the named-mixture and thermochemical paths. The
named source is not converted to an invented NASA model to enable motion.

The native boundary header and rows are:

```text
boundary finite-reservoir-actuation APERTURE MAX_REL_VOLUME SEGMENTS
hold END_SAMPLE VOLUME_RATE_M3_S HEAT_POWER_W
ramp END_SAMPLE START_VOLUME_RATE START_HEAT_POWER END_VOLUME_RATE END_HEAT_POWER
```

`MAX_REL_VOLUME` must be strictly between zero and one. Holds and ramps use
the existing observation-clock endpoints. The current rate and power are sampled
at each original physical interval's start; the original physical interval
duration supplies the volume increment and heat:

```text
delta_volume = volume_rate * physical_dt
supplied_heat = heat_power * physical_dt
```

Callback size, observation resampling and PCM decimation do not retime those
physical increments. A ramp's endpoint remains an exclusive control endpoint,
as in the other native drive histories. Its within-interval hold approximation
needs physical-timestep refinement.

The source transaction first reconstructs the original fluid exchange at fixed
volume. Before publication, the resulting parcel calls its existing
`ClosedGasParcel::midpoint_heat_volume` with the requested final volume and
whole-interval heat. That owner solves the coupled heat/volume update in
conservative coordinates, retains the original caloric and transport sources,
and supplies one common midpoint pressure for the work:

```text
work_by_gas = midpoint_pressure * delta_volume
source_energy_change = original_NSF_energy_out + supplied_heat - work_by_gas
```

The chamber update checks this complete balance using original-law energy
increments and the actual final mass. The original NSF total energy already
contains the opening's pressure and viscous work; those terms are not added
again. Volume work belongs to the chamber motion and remains a separate
original receipt.

The physical volume increment must fit its declared relative bound and retain
a positive, representable final volume. The existing source energy-movement
bound now includes the sum of absolute original fluid substep energy transfers,
absolute supplied heat and absolute chamber work. Opposing heat, work and flow
cannot hide excessive movement through cancellation. Material-domain, original
caloric solve, source-ghost, mass/energy and cancellation checks all precede the
same publication of fluid, ducts and source.

`NsfFiniteReservoirFrame::volume_m3` remains the starting volume;
`final_volume_m3` is the accepted final volume. Its optional `actuation`
contains the original midpoint pressure and `GasHeatVolumeStep`, including
supplied heat, work, storage residual, caloric iterations and entropy change.
That entropy change is a diagnostic, not an entropy-production or exact
isentrope certificate.

Midpoint integration applies to the closed heat/volume stage only. The complete
flow-then-actuation composition remains first order, and the accepted pressure
and temperature affect the next resolved interval. This is imposed motion of a
well-mixed source, without a solved piston, diaphragm, respiratory tissue or
resolved internal chamber field. The acoustic ducts' separate mean-flow
transport limitation remains unchanged.

### Actuated chamber starter

`crates/fs-couple/examples/resolved-dg-actuated-reservoir.performance` retains
the supplied-source chamber example's geometry, source packs, separate escape,
both ducts and material microphone. It adds the following 8 ms program:

```text
boundary finite-reservoir-actuation 0 0.01 4
ramp 96 0 0 -1e-7 0.0005
hold 192 -1e-7 0.0005
ramp 288 -1e-7 0.0005 1e-7 0
hold 384 1e-7 0
```

The source compresses while heat power rises to 0.5 mW, holds those rates,
transitions to expansion as heat falls, then expands without supplied heat.
Use the same three nitrogen packs and native render command above, replacing
the input filename and choosing a fresh output path. Instrument Lab imports
the canonical file as **Resolved flow — chamber work and heat**.

The small imposed motion is an authored numerical exercise. It is not measured
lung motion, a playable phrase, or evidence of flute sound quality. Regression
sources compare the original material candidate, complete energy balance,
feedback, rollback and native control timing; native execution and physical
refinement remain required before stronger claims.
