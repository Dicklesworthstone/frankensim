# Shared material transport on the enthalpy timeline

`fs_lbm::conjugate::enthalpy::transport::advance_joule_enthalpy_with_material_transport`
connects phase/state-dependent thermal transport to the existing electrical and
conservative enthalpy solvers. Enable the existing `enthalpy-coupling` feature.
No additional solver, source database or caloric law is introduced.

## One material evaluation, two physical uses

The bulk callback receives `(electrical_voxel_index, thermal_cell, kelvin)` and
returns `EnthalpyTransport<F>`: the original concrete material/phase/source object,
Ohmic resistivity and isotropic thermal conductivity. The cell carries invariant
reference mass, material ID, geometry and specific enthalpy. Use that enthalpy to
query the existing phase chart: equal Kelvin temperatures on a latent plateau do
not imply equal phase fraction or transport properties. This response never
replaces the original storage law with a new density or heat capacity.

Every mapped cell is sampled once per epoch, INCLUDING ideal electrode cells
which still transport and store heat. The electrical solver ignores electrode
bulk resistivity. The thermal-link callback receives the exact two ordered bulk
samples plus the original `ContactPatch`. It supplies total conductance [W/K]
and an optional concrete source/law payload. There is no automatic interface
resistance or mixing rule. Source-backed consumers can retain the original
state-point resolution and selected claims in F; synthetic scalars do not
acquire measurement authority by being wrapped here.

Explicit `link_ids` select the thermal links to refresh. They must be unique,
active and have two mapped receiving cells. Other links retain their original
constant conductance, including links to additional unmapped thermal cells.
Missing selected links or unsupported source states refuse; they do not silently
fall back to the authored value. Resolved zero thermal conductance insulates but
does not uncover either surface. Electrical insulation is explicitly None.

## Actual substeps, unchanged history and energy ownership

The default-empty `DepositionObserver::conductances` hook follows `sources` at
the same unchanged substep start, after all events due at that time. The thermal
owner uses the resulting values in its original pairwise heat exchange AND its
CFL calculation. Authored link values are restored before observation/events,
so later legacy calls cannot retain stale transport. Per-side covered areas,
cell masses, enthalpy, topology and the event cursor are not reinitialized.

The existing electrical solver deposits bulk/contact Joule heat once. Its work
account uses the thermal owner's actual selected duration, not a separate clock.
Thermal `external_heat_j` already includes this source. Internal thermal links
exchange equal opposite power and are not another external heat contribution.

The report retains the last applied material/link epoch, its matching electrical
report, and separately evaluated endpoint material/link/current support. Endpoints
supply no heat retroactively. No Clone is required on material payloads; retained
Arcs share immutable evaluations, not an advancing physical history. Failed source
queries, link admission, current budgets, heat balances or cancellation return no
candidate. Interior callback side effects remain caller-owned and must not debit
reservoirs or advance irreversible material history during evaluation.

The whole-call current-solve budget also limits sampling epochs: each substep
and the final endpoint consume one solve. Memory adds only the current/last/endpoint
material and link samples to existing thermal/current scratch, not a trajectory
of full fields. Cell and link counts remain subject to the existing budgets.

## Data-only executable

`joule_phase_change` retains its original arguments and four-column tables.
Adding a fifth column to EVERY enthalpy row declares local conductivity in
W/(m K) and selects the shared transport API. Electrical resistivity and
conductivity are interpolated with the SAME enthalpy bracket; no extrapolation
is permitted. The original 64 KiB and 256-knot bounds remain. The header's
conductivity is the authored link reference, not a fallback for missing rows.
The extended output appends endpoint centre conductivity; legacy CSV is unchanged.

```sh
cargo run -p fs-lbm --features enthalpy-coupling --example joule_phase_change -- \
  crates/fs-lbm/examples/joule_phase_change/transport_constant.txt 10 60 0.05
cargo run -p fs-lbm --features enthalpy-coupling --example joule_phase_change -- \
  crates/fs-lbm/examples/joule_phase_change/transport_phase_dependent.txt 10 60 0.05
```

Only conductivity differs between these synthetic files. Density, caloric/latent
curve, resistivity, geometry and voltage pulse are identical. The example declares
a two-half-cell series-resistance law with no extra thermal contact resistance.
The same electrical solver and enthalpy kernel run both cases. Voltage shutoff
preserves heat/phase history and continuing internal heat redistribution.

## Limits and verification

This is explicit stationary, fixed-volume caloric coupling. Rates are frozen
within each actual substep. The frozen-conductance CFL is not a bound on arbitrary
nonlinear material/source stiffness; max_dt refinement remains necessary. There
is no melt flow, moving front, phase-density volume work, induction, implicit
feedback, mechanical dissipation heating or experimental validation. The separate
LBM momentum paths are unchanged. Significant phase-volume changes need another
admitted physical model rather than an in-place density edit here.

Five thermal-hook, seven electrical/transport and two example regressions were
added. They exercise the actual Rust APIs, but compilation and execution remain
unverified: RCH/Cargo/rustc were unavailable in the implementation environment.
Independent nonlinear two-body equations show first-order time refinement; four
three-cell data-file reference runs close heat accounting below 1.4e-13 J. These
are mathematical references, not native executable or measured material results.

```sh
rch exec -- cargo test -p fs-conduction --test material_substep_transport
rch exec -- cargo test -p fs-lbm --features enthalpy-coupling --test joule_enthalpy --test joule_enthalpy_transport
rch exec -- cargo test -p fs-lbm --features enthalpy-coupling --example joule_phase_change
```
