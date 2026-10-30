# Finite gas cooling into radiating walls

This executable cools two finite gas layers into two finite wall profiles. It
composes the original NASA gas inventory, two-zone network, stored wall solver
and participating enclosure radiation. Wall storage and radiative exchange are
solved together at the wall endpoint temperatures. Accepted gas and wall states
advance together, with identical internal heat exchanges booked to each owner.

## Run

From the repository root:

```bash
cargo run --release -p fs-airflow --features compartment-fire \
  --example radiating_compartment -- gray 10 0.1
cargo run --release -p fs-airflow --features compartment-fire \
  --example radiating_compartment -- spectral 10 0.1
cargo run --release -p fs-airflow --features compartment-fire \
  --example radiating_compartment -- spectral 10 0.05
```

Arguments are optical mode (`gray` or `spectral`), end time [s] and maximum
timestep [s]. Defaults are `spectral 10 0.1`. End time must be 0.001..300 s,
maximum timestep 0.001..0.5 s, and requested and accepted work is bounded by
50,000 steps. CSV goes to stdout; model attribution and accumulated numerical
diagnostics go to stderr. `--help` prints the invocation.

The [executable](../radiating_compartment.rs) is a numerical coupling example.
Its optical and wall properties are explicitly authored inputs. It provides no
experimental validation or molecular absorption data.

## Original owners and coupling API

The [wall/radiation adapter](../../src/compartment/radiation/walls.rs) exposes
`ParticipatingEnclosure::wall_workspace`, `candidate_room_walls` and
`candidate_walls`, with `WallOptics::Gray` or `WallOptics::Spectral`. The room
entry uses the actual original gas temperatures and clock, and rejects optical
or convective contact through an absent layer. The direct-temperature entry
allows a caller to supply its own gas endpoint temperatures explicitly.

For trial surface powers q, each original
[`Wall1d`](../../src/compartment/wall.rs) candidate supplies its endpoint face
temperature T_w(q). The original [gray enclosure](../../src/compartment/radiation.rs)
or [spectral enclosure](../../src/compartment/radiation/spectral.rs) then supplies
R(T_w(q), T_g). A bounded Newton solve with backtracking solves q - R = 0.
Conduction, enthalpy storage, half-cell surface resistance, radiosity and Planck
integration remain in their original modules. Each Jacobian column perturbs
only its own wall; prepared wall, optical and dense-LU workspaces are reused.
Original thermal trials and source-bearing reports still allocate.

The runnable example demonstrates the complete acceptance sequence:

1. Prepare one workspace for the original walls and enclosure. Supply one wall
   per enclosure surface, matching areas and accepted clocks. Supply explicit
   layer absorption, optical paths, convective contacts and exterior boundaries.
2. Inside `ZoneNetwork::step_with_sources`, call `candidate_room_walls` using
   the attempted duration. Resolve each original material-card emissivity at
   every requested trial temperature. Previous accepted surface powers are an
   optional numerical initial guess, not another prescribed heat source.
3. Call `candidate.report.add_gas_heat(&mut rates[room_index])` once. This adds
   **all wall convection plus enclosure gas radiation**. Adding the inner
   radiation report's heat or individual walls' convective heat again would
   double-count an internal exchange.
4. Return the wall candidate as the network's source payload. Install its walls
   only after the complete network step accepts. A refused attempt or
   cancellation leaves the accepted gas and wall owners unchanged; a shorter
   attempt recomputes the duration-dependent candidates.

At acceptance the same original radiative powers feed the returned walls and
gas account. A final optical evaluation at the actual returned wall faces gates
their temperature and heat-power consistency. No imbalance is redistributed
to a gas layer. Exterior convection and prescribed external irradiation remain
separate from internal enclosure exchange.

`EnclosureWallSolve::default()` permits 32 Newton iterations, 24 backtracking
probes per update, 4,096 radiation evaluations and 16,384 wall evaluations. Its
endpoint gates are 1e-4 W and 1e-7 K, with a 1e-4 J per-candidate combined-energy
gate. Reports retain original optical sources, actual queried temperatures,
residuals and work counts. Spectral quadrature and omitted wavelength coverage
remain separate diagnostics from wall/radiation nonlinear residuals.

## Explicit physical inputs

| Input | Value |
| --- | --- |
| Room floor area / height / volume | 2 m² / 2 m / 4 m³ |
| Initial lower / upper temperature | 550 K / 900 K |
| Initial interface / pressure | 1 m / 101,325 Pa |
| Gas composition | 0.79 N2 / 0.21 O2 by mole |
| Initial floor / ceiling wall temperature | 300 K / 350 K |
| Each wall area / thickness / cells | 2 m² / 0.02 m / 4 |
| Wall density / heat capacity / conductivity | 500 kg/m³ / 1,000 J/(kg K) / 1 W/(m K) |
| Floor / ceiling emissivity | 0.7 / 0.9 |
| Interior convection coefficient | 8 W/(m² K) |
| Lower / upper absorption coefficient | 0.2 /m / 0.7 /m |
| External convection / irradiation | Zero |

The floor contacts only the lower gas and the ceiling only the upper gas.
Wall calorics use T = 300 + h/1000 [K], admitted over 200..1200 K; the authored
emissivity cards share that temperature domain and state no uncertainty.
The original [NASA9 air source](../compartment_fire/air.rs) retains the original
200..1000 K N2/O2 coefficients from McBride, Zehe and Gordon,
NASA TP-2002-211556, as distributed in
[Cantera v3.2.0 airNASA9.yaml](https://github.com/Cantera/cantera/blob/v3.2.0/data/airNASA9.yaml)
(Git blob `5975d1303819332c1bd77a783d9873b66f65c952`). No caloric extrapolation
or refitting is used.

Optical geometry is one reciprocal normal chord between infinite parallel
plates. Segment lengths are regenerated from the accepted gas interface and
room height before each interval, then held during that interval and its
retries. This is an authored angular approximation, not exact hemispherical
slab transport or an automatically generated finite-room enclosure.

The absorption coefficients are numerical inputs, **not measured N2/O2
spectroscopy**. Spectral mode uses one constant-opacity band with those same
coefficients over 1e-7..1e-3 m and the same gray surface emissivities. It therefore
checks the finite-window gray limit; its small difference from gray mode is
not evidence of improved nongray physics. The separate
[coupling tests](../../src/compartment/radiation/walls/tests.rs) also exercise
selective-band optics. No soot population, combustion source or chemistry is
included in this example.

## Actual retained runs

All three native runs reached 10 s, returned success and required zero whole-step
halvings. The normal Cargo release build independently reproduced the coarse
spectral CSV exactly. Values below are final states or the maximum absolute
residual over the run, as labelled.

| Quantity | Gray, dt ≤ 0.1 s | Spectral, dt ≤ 0.1 s | Spectral, dt ≤ 0.05 s |
| --- | ---: | ---: | ---: |
| Lower temperature [K] | 489.336353128 | 489.336360861 | 489.393112232 |
| Upper temperature [K] | 553.012435004 | 553.012461790 | 553.680766594 |
| Pressure [Pa] | 76204.542545326 | 76204.544765415 | 76247.392331867 |
| Interface height [m] | 1.182988556483 | 1.182988540712 | 1.182460876374 |
| Floor inner face [K] | 323.601962511 | 323.601959182 | 323.609350739 |
| Ceiling inner face [K] | 380.199205218 | 380.199204469 | 380.226252746 |
| Maximum combined energy residual [J] | 2.984051e-5 | 3.388710e-5 | 7.009879e-5 |
| Maximum endpoint power residual [W] | 7.338125e-5 | 7.337270e-5 | 8.082434e-5 |
| Maximum endpoint temperature residual [K] | 9.063035e-8 | 9.061989e-8 | 9.982307e-8 |
| Accepted steps | 101 | 101 | 200 |
| Wall / radiation evaluations | 1616 / 1010 | 1616 / 1010 | 3196 / 1994 |
| Wavelength evaluations | 0 | 332160 | 657792 |

With zero external heat, the checked whole-run balance is the change in
original gas physical internal energy plus the change in both stored wall
energies. The executable refuses an absolute defect above 0.05 J, per-layer
species changes above 1e-10 kg, or different accepted gas/wall clocks. Its gas
heat columns include both convection and radiation. The negative values record
gas cooling, with opposite heat entering the walls.

At dt ≤ 0.1 s, spectral minus gray changes the final upper gas temperature by
2.6786e-5 K and pressure by 0.002220089 Pa. Halving the spectral timestep changes
the lower and upper temperatures by +0.056751371 K and +0.668304804 K, pressure
by +42.847566452 Pa, interface by -0.000527664338 m, and floor/ceiling faces by
+0.007391557 K / +0.027048277 K. This is one timestep refinement on the same
four-cell walls; it does not establish time order or mesh independence.

Gas temperatures, contact geometry and optical coefficients are frozen during
each source interval. Implicit wall/radiation feedback therefore leaves the
complete gas/source integration first order. Small conservation residuals do
not establish time accuracy, a fully implicit gas or chemistry solution, an
entropy theorem, or unconditional convergence. The authored optics, angular
quadrature, wall properties and finite spectral window remain separate model
approximations.

[Machine-readable inputs, commands and results](../../../../data/reference/compartment_fire/radiating_compartment_example.json)
retain final CSV values, accumulated work and both differences above.
