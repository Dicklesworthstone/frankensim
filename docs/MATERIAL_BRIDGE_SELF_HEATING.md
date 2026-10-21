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
Existing bridge growth may provide the next accepted specimen, but this module
does not itself perform atomic spatial birth/contact synchronization.
