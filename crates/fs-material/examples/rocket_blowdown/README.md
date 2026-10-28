# Finite-reservoir rocket nozzle

This example runs the actual variable-heat-capacity rocket nozzle and finite
adiabatic reservoir APIs. Gas mass, temperature and pressure fall as propellant
leaves. The same evolving source model determines choking, supersonic exit
conditions, momentum thrust, pressure thrust and specific impulse. Thrust is
integrated with the reservoir's midpoint mass-flow step.

## Run the included manufactured case

From the repository root:

```sh
cargo run -p fs-material --example rocket_blowdown -- \
  --synthetic crates/fs-material/examples/material_gas_spring/mixture-a.tsv \
  500 300000 > rocket.csv
```

This reuses the existing gas-spring TSV loader and data. Its two species and
caloric/transport laws are **synthetic numerical examples**, not measured rocket
propellants. The initial tank is 500 K and 300,000 Pa. Defaults are 100 steps of
0.02 s, a 0.002 m³ tank, a 10⁻⁶ m² throat, exit/throat area ratio 1.5 and vacuum
ambient. These defaults keep the complete example trajectory within the supplied
200–2,000 K and 1,000–10⁷ Pa source coverage.

Optional arguments are `--steps`, `--dt`, `--volume`, `--throat-area`,
`--area-ratio`, `--ambient`, `--temperature-min` and `--temperature-max`, all in
coherent SI units. `--help` prints the defaults. The two temperature options
bound numerical root searches; they do not enlarge material-source validity.
At most 10,000 steps and 64 arguments are accepted. Halve `--dt` and double
`--steps` to assess transient convergence over the same physical duration.

## Run supplied canonical source packs

```sh
cargo run -p fs-material --example rocket_blowdown -- \
  --packs reservoir.gas species.pack nasa9.pack transport-material-card.pack \
  --steps 100 --dt 0.02 --volume 0.002 --throat-area 0.000001 \
  --area-ratio 1.5 --ambient 0 --temperature-min 200 --temperature-max 1800
```

This route uses the existing `read_source_parcel` owner. Supply one canonical
species pack, its NASA-9 model pack, its canonical **material-card** transport
pack, and an existing gas record such as `thermochemical 500 300000 mole 1 N2 1`.
The record's species identity must match the supplied pack; names do not select
or synthesize material data. Records are bounded to 8 KiB, each pack to 16 MiB,
and synthetic input to 64 KiB. Frozen-gas phase policy is explicit; no reaction
or condensation model is implied.

Every accepted tank, throat, shock and nozzle-exit query must fit the actual source
temperature **and pressure** domains. The included nitrogen and helium
near-atmospheric transport packs cover only 75,000–101,325 Pa, so they do not
cover the pressure drop across a choked nozzle. This example refuses those
unsupported states; it does not extrapolate or substitute an air law. Supply
appropriately supported source packs for a physical propellant study.

## CSV interpretation

Rows contain reservoir time/mass/temperature/pressure, exhaust mass flow,
exit Mach and pressure, thrust, specific impulse, accumulated impulse,
expelled mass, exported stagnation enthalpy and mass/energy residuals. The
first flow is evaluated at time zero. Thereafter reservoir columns describe
the step endpoint and flow columns describe the mass-midpoint state used by
the integrator; `flow_time_s` identifies the associated midpoint time.

`energy_balance_residual_j` is the tank's internal-energy change plus exported
stagnation enthalpy. It is not normalized by formation energy. Exported
stagnation enthalpy already includes nozzle kinetic energy and pressure work;
neither is an additional tank-energy sink. The two nozzle residual columns and
step mass-quadrature residual expose the underlying numerical closures.

The nozzle is ideal, inviscid, adiabatic and frozen in composition. This
example has no combustion, real-gas EOS, nozzle heat transfer, external vehicle
motion or drag. The nozzle supports source-valid unchoked discharge, internal
normal shocks and the attached ideal external-compression branch. It does not
predict viscous separation, side loads or transient shock motion. A failed step prints its time and
reason and exits with failure; earlier CSV rows remain the accepted prefix.

## Advance a vehicle with the same finite gas

The coupled example consumes the exact mass and impulse accepted by this
reservoir/nozzle step, including pressure thrust once:

```sh
cargo run -p fs-couple --example rocket_flight -- \
  crates/fs-material/examples/material_gas_spring/mixture-a.tsv 500 300000
```

It uses the same manufactured gas and tank/nozzle geometry, a 0.1 kg dry
vehicle, vacuum, +x thrust, and zero gravity or other external force. Its CSV
adds position, velocity and gas/vehicle mass and impulse residuals. Optional
positional arguments after pressure are step count, timestep and dry mass.
The public coupling API also accepts explicit gravity, external force and
thrust direction in a matching inertial frame. This is point-vehicle
translation; attitude, atmospheric drag, slosh and combustion are not modeled.
