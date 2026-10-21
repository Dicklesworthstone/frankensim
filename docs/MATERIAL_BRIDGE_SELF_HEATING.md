# Spatial material self-heating and gravity sag

`fs_couple::strip_enthalpy::bridge` couples the existing deposited bridge to the
shared enthalpy network. Unlike one uniform section, each axial cell and layer
has its own temperature and Maxwell history. The original beam integrates its
compatible deflection, reactions and independent gravity/dead-load work.

`Bridge::advance_with_cell_energy` retains the physical accounts from the actual
section updates, keyed by `(axial_cell_id, layer_id)`. It does not advance a
material point twice or distribute total loss by mass. The original unrecorded
bridge/support path keeps its existing arithmetic and checks.

## One coupled interval

Supply a `BridgeThermalBinding` for each born axial-cell/layer, with its original
`LayerThermalBinding`. Thermal CellIds need not equal either mechanical ID. All
receiving cells must be distinct; extra passive thermal cells are allowed.
The original caloric object, mechanical law, material source, reference volume
and layer mass must agree. Each axial cell's thermal masses must also sum to the
reference mass used by gravity. The caloric owner admits the original histories;
no material point is reinitialized to make an incompatible binding fit.

`BridgeEnthalpyConfig` reuses `StripEnthalpyConfig` for section/nonlinear controls
and adds the original bridge's integrated deflection limit. Loads are endpoint
values: nonnegative gravity [m/s2] and one signed downward line load [N/m] per
axial cell in `Bridge::cells()` order. Initial energy guesses [J] instead use
sorted thermal CellId order, including passive cells.

```rust,ignore
let transport = network.prepare(&thermal, cancelled)?;
let candidate = fs_couple::strip_enthalpy::bridge::step_transport_candidate(
    &thermal, &beam, &bindings, &transport, target_s, gravity_m_s2,
    &additional_downward_n_m, config, &energy_guesses_j, cancelled,
)?;
let (thermal_next, beam_next, report) = candidate.into_parts();
```

Every Newton and finite-difference trial starts from identical old beam and heat
histories. The original thermal residual moves/stores each layer's physical
Maxwell heat; numerical dissipation remains separate. Pinned material transport
is evaluated at that actual proposed endpoint. The original beam then integrates
the resulting spatial strains/curvatures. No second Maxwell, beam, heat-rate or
nonlinear solver is introduced. The generic `step_candidate` accepts an explicit
coefficient callback for independent heaters/boundaries; it must not include
Maxwell heat again or mutate external physical histories.

Fresh final admission checks local thermal residuals, total thermal balance and
`delta elastic + delta thermal + numerical loss - load work - other heat`.
Gravity is part of the independently integrated load work; do not count it again
as heat or also subtract gravitational potential. The returned thermal audit
already includes physical material heat. A late source, budget, mechanical or
cancellation failure returns neither candidate history.

## Verification and limits

Two bridge regressions compare recorded/unrecorded states under both built-in
supports, an independent 49:1 axial heat ratio, and late rollback. Five coupled
regressions cover actual endpoint replay, nonuniform heating, conductivity-only
changes to sag, lossless numerical damping, unloading memory/whole-history energy,
reference-mass and source-context refusals, and final cancellation.

```
rch exec -- cargo test -p fs-solid --test bridge_layer_energy
rch exec -- cargo test -p fs-couple --test bridge_enthalpy
```

Native compilation/tests were unavailable in the editing environment. Independent
NumPy/SciPy equations passed 54 parameter cases plus named histories, with maximum
combined-energy discrepancy 1.55e-15 J. In the synthetic two-cell cantilever,
changing only conductivity from 0.0002 to 0.02 W/(m K) changes downward tip
displacement from 5.8818604 to 5.8778790 mm under unchanged load/reference mass.
These are equation checks, not native Rust or empirical material validation.

This is the existing quasistatic, slender, small-strain, perfectly bonded beam
with left clamp or pin/roller supports, not drawing-tension or moving-contact
feedback. Thermal geometry is fixed, and each axial-cell/layer has one temperature.
The caloric model still requires constant density, temperature-independent spring
energy/free thermal strain and a single solid phase; conductivity and relaxation
may vary within source support. No inertia, melting, yielding, buckling, adhesion
or automatic mesh/time-error estimate is supplied. Split at scheduled events and
load changes: endpoint events remain pending until both histories are accepted.
The base interval leaves target-time events pending. The `growth` composition
below accepts those events together with the solved interval.

## Growing spatial material on one accepted history

`bridge::growth::GrowingBridge` binds an already born single-layer clamped seed,
its complete `BridgeThermalBinding` list, a selected thermal body ID and bounded
future `TipBirth` descriptors. Each future descriptor supplies its source binding
and unstrained projected length; the actual scheduled thermal cell supplies mass
and enthalpy. The shared admission checks require density, reference volume,
caloric object, original mechanical law and current temperatures to agree.

`step_candidate` and `step_transport_candidate` first use the original spatial
coupled interval. The thermal timeline then applies its target-time events in
the authored order, with zero thermal substeps. The existing bridge append owner
initializes each new section stress-free along the accepted tip tangent. Old
positions, branch histories and sources are preserved. Appends poll cancellation
inside their original per-cell checks; none run inside Newton/Jacobian trials.

The result contains interval accounts and a separate `BirthAudit`, including the
original thermal incoming mass/enthalpy and ordered mechanical tip reports. A
final combined storage/work gate includes incoming enthalpy once. Failure in the
last birth rejects the entire interval and all earlier births. Other thermal
bodies may arrive without acquiring invented beam sections, but a registered tip
cannot change physical body. A missing selected-body binding is an error.

Callers still split at events and load changes. Loads and energy guesses describe
the old topology. After acceptance, use `PreparedMaterialNetwork::extend_thermal`
with the complete new cell/contact description to retain existing source pins.
The thermal contacts remain explicitly authored fixed-reference geometry, not
contacts detected from the deflected beam. Birth is not an instantaneous impact
or equilibrium solve: `birth_equilibrium_pending` remains true until a subsequent
positive-time coupled advance. Multi-layer tip extension, drawing tension,
inertia, phase change and a moving nozzle constraint remain outside this model.

Six additional Rust regressions cover manual-owner equivalence followed by
source-preserving conduction and sag, schedule-order growth, final invalid-birth
rollback, cancellation/replay, current-time initialization with passive enthalpy,
and bounded complete registration. Run with:

```
rch exec -- cargo test -p fs-couple --test bridge_enthalpy --test bridge_enthalpy_growth
```

These Rust tests remain unexecuted here because DSR/RCH/Cargo are unavailable.
Independent Python/SciPy equations checked 27 one-to-two-cell growth histories;
maximum energy discrepancy was 4.41e-16 J. These checks confirm the reference
equations and tip-tangent continuity, not execution of the Rust event composition.
