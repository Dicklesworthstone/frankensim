# Electrical heating on the conservative enthalpy timeline

Enable `fs-lbm`'s `enthalpy-coupling` feature to use
`fs_lbm::conjugate::enthalpy::advance_joule_enthalpy`. It connects the existing
voxel current/contact solver to the existing `fs-conduction` enthalpy network.
It is not another thermal solver or an apparent-heat-capacity approximation.
Default LBM consumers do not enable the additional conduction dependency.

## Sources use the actual thermal substeps

The new optional `DepositionObserver::sources` hook supplies signed watts at
active cells, after start-time events and before each original thermal substep.
Those powers are added once to authored boundary sources, then removed from the
boundary description before observations and subsequent events. Omitted entries
mean zero additional power. Existing observers have an empty default hook.
The original conduction CFL, event ordering, reference masses, enthalpy updates,
thermal laws, conservation checks, substep budgets and atomic failures remain.
The source owner receives the actual duration, not a separately chosen clock.

The electrical adapter uses that hook to compute cellwise bulk AND contact Joule
heat. One unique active thermal `CellId` receives each electrical voxel's watts.
The explicit correspondence must have matching volumes and centreline midpoints
within floating-point geometry roundoff. Extra unmapped thermal cells can exchange
heat through existing links without acquiring an electrical source. Mapped cells
must be initialized before the advance. Electrical contacts do not implicitly
create thermal links, thermal resistance or a mechanical boundary.

The callback receives `(voxel_index, thermal_cell, temperature_kelvin)`.
`thermal_cell` retains material ID, reference mass and **specific enthalpy**.
Use the original material chart to resolve phase fraction or another supported
caloric coordinate. Temperature alone cannot recover position on an isothermal
latent-heat plateau. `Some(resistivity)` declares finite positive Ohmic response;
`None` explicitly declares insulation. Unsupported data must be an error.
The ideal electrical electrode cells do not need bulk resistivity queries, but
remain finite thermal receiving cells: declared contact heat is not discarded.

## Accounting and failure behavior

The original current solver owns the current distribution and contact heat
partition. Its heat is deposited once; never add a second I-squared-R source.
The adapter independently checks summed electrode power, the electrical heat
account and mapped cell powers. Their energies are integrated over the actual
thermal substeps. The original `EnergyAudit::external_heat_j` already includes
Joule heating together with other boundary sources/exchanges; adding the new
Joule account again would double-count it. Incoming mass/enthalpy from unrelated
thermal birth events retains its original separate account.

The result distinguishes the last **applied start-state** electrical solution
from a fresh **endpoint diagnostic/support** solution. The latter supplies no
heat during the preceding step. Its work counts against the explicit whole-call
current-solve budget. Invalid geometry, missing material support, an exhausted
budget, failed energy checks or cancellation returns no candidate; the original
thermal state remains unchanged. Source-query receipts and any interior callback
side effects remain caller-owned. The callback should evaluate an immutable law,
not debit an external reservoir or advance kinetic history during a proposal.

## Model boundaries

This is stationary, fixed-geometry Ohmic/caloric coupling. The thermal law and
network must be justified for the represented control volumes. Spatial heat
exchange follows the original network links. The conduction CFL does not bound
source stiffness: choose and refine `StepConfig::max_dt_s`. Currents and contact
coefficients are held over each substep; this is not the magnetic midpoint API.
Drive changes must occur at separate advance endpoints. Existing wired electrical
admission and insulated exterior faces remain; no artificial reference electrode
or implied current leakage is introduced.

Latent and sensible heat live in the SAME enthalpy coordinate. The adapter does
not add latent heat to a second capacity account. This does not implement melt
flow, a moving free surface, phase-density volume changes, pressure/expansion
work, magnetic induction, kinetic phase transitions, or experimental validation.
A chart with significant phase-volume changes needs an appropriate volume/work
model; it cannot be silently treated as fixed-volume physics here.

## Data-only example

```sh
cargo run -p fs-lbm --features enthalpy-coupling --example joule_phase_change -- \
  crates/fs-lbm/examples/joule_phase_change/latent_20k.txt 10 60 0.05
cargo run -p fs-lbm --features enthalpy-coupling --example joule_phase_change -- \
  crates/fs-lbm/examples/joule_phase_change/latent_40k.txt 10 60 0.05
```

The example has three fixed 1 cm caloric cells, including thermally evolving
ideal-electrode cells. It uses 300 K initial/ambient temperature, explicit
0.1 W/(m2 K) exterior heat transfer, geometry-derived thermal link conductances,
and a voltage pulse switched to zero exactly halfway through. It reports actual
centre temperature, liquid mass fraction, enthalpy, endpoint current, stored
energy, electrical work and the original thermal energy balance. Source shutoff
does not reset thermal or phase history. Zero voltage shorts the ideal electrodes;
this stationary model has no motional source to sustain current after shutoff.

The bounded example-only table begins with density and constant conductivity in
SI, then 2..256 `h[J/kg] T[K] liquid_mass_fraction resistivity[ohm m]` rows. It
rejects inputs above 64 KiB and extrapolation. The existing phase-curve owner
validates the caloric chart; resistance is explicitly linearly interpolated in h.
The two synthetic files change the latent interval from 20 to 40 kJ/kg while
retaining density, conductivity, sensible slopes and the same resistance versus
phase fraction. They are not source-backed measurements of a named substance.
Their equal phase densities make the neglected density-change work explicit.

```sh
cargo test -p fs-conduction --test generated_enthalpy_sources
cargo test -p fs-lbm --features enthalpy-coupling --test joule_enthalpy
cargo test -p fs-lbm --features enthalpy-coupling --example joule_phase_change
```

Five source-hook, seven electrical/enthalpy and two example tests are authored.
Cargo/rustc were absent in the implementation environment, so native compilation,
Rust tests, executable runs and lockfile regeneration remain unverified/outstanding.
Independent scalar/circuit and three-cell heat equations check the described
accounts and temporal refinement; those are not execution of repository Rust.
