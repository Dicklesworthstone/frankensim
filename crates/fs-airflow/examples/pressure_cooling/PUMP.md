# Pump-driven material cooling

The existing `pressure_cooling` command now accepts a fixed-speed pump head
curve instead of a prescribed driving pressure. It uses the SAME material
importer, source-backed caloric/transport laws and actual tetrahedral conduction
callback. One network branch represents series flow; two represent parallel
flow. The existing fixed-pressure forms remain unchanged.

```sh
cargo run -p fs-airflow --example pressure_cooling -- --parallel --pump \
  crates/fs-airflow/examples/pressure_cooling/pump.txt 20 90 \
  crates/fs-airflow/examples/pressure_cooling/temperature_dependent.txt \
  20 0.0002 0.012
```

Omit `--parallel` for the original two-region series circuit. After the pump
file, supply minimum/maximum NETWORK differential pressure [Pa], material file,
solid heating [W], and each branch's minimum/maximum mass flow [kg/s]. Negative
heat extracts energy. Pressure and flow brackets must support every trial;
there is no silent expansion, extrapolation or fallback material/source.

The pump table contains required `source`, `identifier`, `license`,
`gravity_m_s2`, `minimum_nu_m2_s` and `maximum_nu_m2_s` headers, followed by
`Q_m3_s H_m` samples (put `#` before an explanatory column header). The source
identifies a fixed-speed curve. Increasing flows and nonincreasing nonnegative
heads are required. The example loader retains the license; the library curve
retains its original samples, gravity, viscosity domain and source identity.
The supplied curve, gravity of 10 m/s2 and coolant tables are synthetic fixtures,
not measured equipment or experimentally validated fluids.

For each pressure trial, the library recomputes branch flow and thermal
transport, solves the common solid, and evaluates the pump at the resulting
TOTAL INLET volume flow. Pump pressure is `rho_in * g * H(Q_in)`. It uses actual
inlet density, not a water constant or a heated-branch density. The inlet
kinematic viscosity must remain in the explicitly supported curve range. This
range is a caller declaration, not a computed viscosity correction. A head
curve from one fluid is not automatically valid for every other liquid.

Both the pump pressure residual and all original wall/reference/coefficient/
branch-flow/regional-power gates must pass. The output adds the actual pump
head, available pressure, pressure residual, inlet volume flow, hydraulic power,
source/license and total network/flow evaluations. The returned library value
retains original solid fields, sent boundary states and recomputed fluid states.
The full command has shared caps of 250000 branch thermal evaluations and 4096
solid callbacks across every pump-pressure trial; 48 pump candidates are allowed.

For the constant table and illustrative pump, independent model equations
predict 50 Pa / 0.01 kg/s total in parallel, or 80 Pa / 0.004 kg/s in series.
At 20 W those imply mixed outlets of 301 K and 302.5 K respectively. These are
analytic checks, NOT claimed output from executing the Rust/FEM program.

## Physical boundary

This is steady once-through cooling with ideal manifolds and a stable,
low-compressibility liquid at the original thermodynamic source pressure.
No recirculating reservoir, pump speed dynamics, affinity-law scaling, efficiency,
shaft/electrical power, viscosity correction, cavitation/NPSH, boiling or pump
motor heat is inferred. There is no uniqueness/stability certificate. Pump
hydraulic output and branch pressure-loss dissipation are reported separately,
not asserted equal when regional density changes. Neither is automatically
deposited as sensible heat. Every branch and the whole network retain the
original 0.01 W neglected-hydraulic-dissipation allowance.

Head, pressure and power conventions follow the pump-specific quantities in
KSB's centrifugal pump lexicon (Pump affinity laws). The code deliberately does
not apply its speed scaling or assume its applicability conditions:
https://www.ksb.com/en-global/centrifugal-pump-lexicon/article/pump-affinity-laws-1117584

```sh
cargo test -p fs-airflow --test liquid_pump
cargo test -p fs-airflow --example pressure_cooling
```

Nine new library-integration tests and four new example tests cover pump and
material substitutions, original result retention, budgets, domain refusals,
signed heat and the actual FEM connection. They were authored but not executed
in the authoring environment: Cargo/rustc are unavailable. Independent Python
pump/network reference equations are not Rust, IQN or FEM execution evidence.
