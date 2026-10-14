# Material-data-driven pressure cooling

This executable connects an authored liquid property table to the production
material resolvers, pressure-selected mass flow, temperature-dependent
convection, sensible-enthalpy transport and tetrahedral solid conduction.
It adds no material equation or separate solid/coupling solver.

From the workspace root, after materializing the pinned workspace dependencies:

```sh
cargo run -p fs-airflow --example pressure_cooling -- \
  crates/fs-airflow/examples/pressure_cooling/constant.txt 50 20 0.0002 0.008
cargo run -p fs-airflow --example pressure_cooling -- \
  crates/fs-airflow/examples/pressure_cooling/temperature_dependent.txt 50 20 0.0002 0.008
```

Arguments after the table path are differential driving pressure [Pa], total
solid heat input [W], and the two positive mass-flow bracket endpoints [kg/s].
Negative heat input extracts heat. The inlet is 300 K. Every trial must remain
within the source and correlation domains; a missing root is a refusal, not a
request to extrapolate or silently choose a different coolant.

The two supplied tables are **synthetic illustrations, not measured liquids**.
Their four material properties agree at 300 K, but their temperature laws differ.
Only change the table path to compare them; geometry, pressure, heating, solver
and numerical controls remain the same. Output includes conserved mass flow,
fluid outlet, both mean solid wall temperatures, independently accumulated
solid and coolant heat flows, pressure loss, omitted hydraulic power and the
largest regional heat mismatch. No expected output is substituted for a solve.

## Table format

Required headers are `material`, `source`, `license`, and `pressure_pa`.
Preserve the actual source citation and license when supplying real data;
these strings do not authenticate measurements or establish their accuracy.
The source pressure is **absolute material-query pressure**, not the differential
driving pressure in the command line. Other source conditions are unsupported by
this small example-only importer; the library supports complete query points.

Each following row has five positive, finite, coherent-SI numbers:

```text
T_K rho_kg_m3 Cp_J_kg_K k_W_m_K mu_Pa_s
```

The column-name line is explanatory; put a leading `#` on it in a real file.
There must be 2–256 strictly temperature-ordered rows. Each property becomes an
ordinary `fs-matdb` curve. The shared liquid law pins those claims, and the
existing caloric owner derives its enthalpy chart from the same Cp and density
sources with an explicit 0.5 J/kg chord-approximation allowance. An invalid or
out-of-support material is never replaced with water or air defaults.

## Model and verification boundaries

The fixed solid is a 0.06 x 0.10 x 0.01 m block with declared conductivity
5 W/(m K). Its two top-face regions each measure 0.003 m2; other faces are
adiabatic. Both regions share one series coolant stream. Each has a declared
0.01 m correlation diameter and 2 m developed-flow length, independently of
its actual wetted patch area. The declared hydraulic factor is 5e9 /m3 per
region, not inferred from the heat-transfer geometry. This is a reduced
conjugate model, not a fully resolved coolant-channel mesh.

Thermal/hydraulic properties use regional bulk-mean temperature. The existing
flow model treats thermodynamic pressure as fixed, assumes a stable liquid
phase and neglects acceleration pressure. It does not model pumps, cavitation,
boiling or fluid/solid deformation. Hydraulic dissipation is reported, not
added to sensible heat, and every trial must fit the explicit 0.01 W omission
allowance. Numerical and model approximations are not material validation.

For the constant table, the model equations predict 0.0025 kg/s at 50 Pa;
20 W heat uptake then predicts an outlet of 304 K. These are analytic checks,
not a claim that this executable has run in the authoring environment.

```sh
cargo test -p fs-airflow --example pressure_cooling
cargo test -p fs-airflow --test liquid_pressure
cargo test -p fs-airflow hydraulic --lib
```

The three example regressions exercise malformed table admission, data-only
material substitution through the actual coupled FEM path, signed cooling,
and missing-pressure-bracket refusal. They were authored but **not executed
locally**: Cargo and rustc are unavailable in the authoring environment.
