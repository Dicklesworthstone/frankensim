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
complete face; its outward normal flow is positive into the DG inlet. The other
face receives either a caller-supplied incoming characteristic in metres per
second or a physical pressure relative to its original reference in pascals.
For a zero-velocity reference, the original boundary uses

`J− = u_normal − (p − p_reference)/(rho_reference c_reference)`.

Only the incoming characteristic on the coupled face is solved. For every trial,
the original NSF owner advances from the same accepted state while the original
DG owner supplies its current step-mean law `p_DG = H + Z Q`. The scalar solve
matches the mean reference-relative NSF face pressure to that law, using the
face's actual integrated continuum volume flow.

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

The complete reference gas state, including transport properties, must equal the
DG inlet gas. The actual complete face area must match inlet area. The DG interval
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
`boundary_history()` exposes those exact original controls; pressure samples are
never returned as characteristic velocities. The original scalar constructors
select the characteristic variant and retain their existing physical arithmetic.

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
observation inlet-pressure
[optional existing observer_gas record, only with a baffled observer]
boundary COUNT
hold END_OBSERVATION_SAMPLE INCOMING_CHARACTERISTIC_M_S
[repeat exactly COUNT increasing hold records]
end
```

Alternatives use the same original parsers and owners:

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
- `boundary pressure COUNT` selects pressure offsets in Pa for each `hold`
  value. `boundary characteristic COUNT` explicitly names the original m/s
  condition. The existing bare `boundary COUNT` has the same characteristic
  meaning and arithmetic. Unknown kinds and trailing fields refuse.

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
path.

Additional regressions live in `fs-lbm/tests/nsf3_open_boundaries.rs`,
`fs-couple/src/nsf_dg/performance/clock.rs`,
`fs-couple/src/nsf_dg/performance/file/tests.rs`,
`fs-couple/src/dg_port/radiation/tests.rs` and
`fs-couple/tests/music_render_wind/resolved.rs`. These cover finite-wall corner
selection, original interval overlap, source-backed temperature/composition
substitutions at fixed geometry, exact-step radiation, native schema refusal,
PCM partition identity and provenance clocks. These are numerical/interface
regressions; they do not establish measured instrument fidelity.
