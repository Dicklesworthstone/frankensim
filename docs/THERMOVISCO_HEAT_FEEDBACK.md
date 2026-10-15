# Load, dissipation, heat and material response in one step

`fs_material::thermovisco::thermal::CaloricPoint` binds an existing Maxwell
material history to the existing caloric chart and uniform thermal-owner
protocol. `advance` preserves prescribed-strain behavior; `advance_controlled`
also accepts `MechanicalControl::Stress { stress_pa, tolerance_pa }`.

At every trial temperature, stress control solves the existing affine law
`sigma = C epsilon + b` for the new total engineering strain. It then executes
the original `Point::advance`, checks the actual stress, and supplies ONLY
physical dashpot dissipation to thermal transport. The returned temperature
changes the relaxation time and therefore creep, stress and heat generation.
The material history, caloric state, mass and duration are held fixed at their
accepted old values throughout all trials. There is no repeated time advancement
inside the nonlinear iteration and no frozen-temperature mechanical endpoint.

The thermal callback is the existing `UniformEnthalpyStep` protocol, for example
`|input| environment.advance(cx, input)` with
`fs_conduction::lumped::LumpedThermalEnvironment`. Its real boundary heat, energy
allowance and concrete report are retained. A finite-reservoir callback must
return tentative reservoir state, not mutate that reservoir during trials.
Only the returned `FeedbackStep` may be published with those other states.

## Distinct energy and control accounts

With temperature-independent spring energies and free strain, the existing
backward-Euler law supplies, per reference volume,

`delta_elastic + physical_dissipation + numerical_dissipation = sigma_new delta_epsilon`.

The thermal owner solves
`mass delta_h = volume physical_dissipation + external_heat`.
Numerical dissipation is not a physical heat source. For stress control, imposed
boundary work is `volume prescribed_stress delta_epsilon`, independent of the
actual constitutive stress. The combined gate uses this imposed work, not a
reaction chosen to make the account close. Its residual, the independent stress
residual, physical heat defect and endpoint-temperature defect are all exposed.
The scalar material tangent remains isothermal; it is not the derivative of the
coupled thermal solve.

All admissible dashpot dissipation is assumed to thermalize locally in this
rung. The chart must be strictly increasing, single-solid-phase and
constant-density; mass must match reference volume times that density. Spring
modulus and free-strain tables must be temperature-independent. Relaxation and
heat capacity may vary with temperature. Variable elastic free energy, thermal
eigenstrain, phase changes and volume work require additional thermodynamic
terms and are refused. Independent model inputs are not automatically qualified
as the same measured material. A material name never selects missing properties.

## Runnable coupled example

The `fs-couple` example `self_heating_creep` uses the actual convection/radiation
owner, not the test suite's analytic adiabatic callback. Its command is:

```sh
cargo run --release -p fs-couple --example self_heating_creep -- --demo 5 2 2000000 4.8 48
cargo run --release -p fs-couple --example self_heating_creep -- --demo 10 2 2000000 4.8 48
```

Arguments are Cp [J/(kg K)], film coefficient [W/(m2 K)], imposed stress [Pa],
duration [s], and an even step count. These inputs change ONLY Cp. Both
experiments apply the same stress for the first half and remove it for the
second, retaining heat and branch memory. A zero film coefficient explicitly
selects insulation. Changing that coefficient changes cooling and the subsequent
mechanical response. The original Biot gate still refuses unsupported lumping.

The example data are deliberately synthetic: E_inf=100 MPa, one 300 MPa branch
with a 1 s relaxation time at 300 K and a multiplier decreasing to 0.05 at
320 K. Cp=5 or 10 is a sensitivity demonstration, not a realistic material
identification. Density is 1000 kg/m3, reference length 0.1 m, square section
area 1e-5 m2; thermal volume/surface are fixed. Conductivity is 10 W/(m K),
ambient is 300 K, and radiation is explicitly disabled. All choices are printed.
CSV separates extension, stress, physical heat, algorithmic loss, imposed work,
boundary heat, caloric/elastic storage and cumulative energy defect. It prints
only after the complete bounded trajectory succeeds.

## Evidence and limits

The twelve `thermovisco_heat_feedback` regressions include the seven previously
pending tests, signed constant-load analytical creep, data-only caloric feedback,
strain/stress endpoint equivalence, unloading, and transactional refusals. Three
example tests exercise the actual thermal solver, changed cooling and Biot/late
cancellation failures. Run:

```sh
rch exec -- cargo test -p fs-material --test thermovisco_heat_feedback
rch exec -- cargo test -p fs-couple --example self_heating_creep
```

Native Rust execution was unavailable during implementation. Independent
Python/SciPy equations checked 288 strain-controlled and 648 stress-controlled
thermal cases. Four 48-step example reference histories close cumulative energy
to 1.48e-14 J; they are not native output or experimental evidence. Isothermal
creep time refinement approaches first order. Finite tolerances, fixed-point
iteration and one successful step do not imply unconditional stability or a
physical error bound. Stress discontinuities use endpoint work and can produce
substantial numerical loss; no unloading impulse or inertial response is claimed.

This is uniform, small-strain constitutive feedback, not a spatial structural
solver, dynamics, metal plastic creep, melting, or thermally induced stiffness
softening. It does not retrofit every rod/beam consumer automatically. Existing
spatial rod and deposition solvers remain the owners of their geometry and
compatibility; they must explicitly compose this admitted feedback where valid.

## Source-backed material substitution

`fs_material::thermovisco::source::ResolvedThermoMaxwell` now lowers one actual
material card into the same `ThermoMaxwell`. The library query declares the
roles of equilibrium modulus, ordered branch modulus/time pairs, common modulus
and time multipliers, and free thermal strain. No property is inferred from a
material name. Dimensional SI property keys, scalar or T-only curves, continuous
support over a fixed-context temperature interval, at most 30 branches and 1024
union knots are admitted. Every source knot inside the interval is retained;
independent branch-specific temperature changes cannot masquerade as one common
shift. Interior conflicts, missing data and extrapolation refuse.

`law()` is directly consumable by the existing deposited-layer, bridge and point
implementations. `samples()` retains the actual selected claims and coordinates;
`keys()` retains their ordered physical roles. `caloric_point()` binds the shared
caloric chart only when it comes from the same card, fixed nonthermal context and
coordinate conventions. It derives reference mass from the source density and
specified volume, and keeps all of the original heat-feedback physical limits.
Source agreement is not a guarantee that a measured specimen fits the model.

The existing executable now reads either a canonical material pack or an editable
SI table. Geometry, ambient conditions, load schedule, solver and CSV fields stay
the same as the original demo. These two commands differ only in the table's Cp:

```sh
cargo run --release -p fs-couple --example self_heating_creep -- --table crates/fs-couple/examples/self_heating_creep/cp5.txt 2 2000000 4.8 48
cargo run --release -p fs-couple --example self_heating_creep -- --table crates/fs-couple/examples/self_heating_creep/cp10.txt 2 2000000 4.8 48
cargo run --release -p fs-couple --example self_heating_creep -- --pack material.fsmat 1 2 2000000 4.8 48
```

For `--pack`, the extra integer explicitly selects 0..30 branches; the example
uses the existing `NormalizedMaterialCardPack` decoder and the declared 280..320 K
span with a legacy `T` axis and no other coordinates. The pack must supply
`equilibrium_young_modulus`, `relaxing_young_modulus_i`, `relaxation_time_i`
(zero-based branch indices), `modulus_scale`, `relaxation_scale`,
`free_thermal_strain`, `specific_heat_capacity`, `density`, and
`thermal_conductivity`. These use dimensional Pa, s, dimensionless, J/(kg K),
kg/m3 and W/(m K) keys respectively. Extra required material conditions refuse;
the library supports caller-declared fixed coordinates beyond this small example.

A table starts with `frankensim.thermovisco.v1`, followed by exactly one
`equilibrium E_inf_Pa`, zero or more `branch E_Pa tau_s`, and 2..1024 ordered rows:
`state T_K modulus_scale relaxation_scale free_strain Cp_J_kg_K rho_kg_m3 k_W_m_K`.
Its temperature span comes from those rows and must include the initial 300 K.
Comments begin with `#`; files are limited to 16 MiB. This is declared input,
not independently sourced material evidence. The two supplied tables are synthetic.
They are converted to ordinary material-card claims, then pass through the same
library source, caloric and transport evaluators as canonical packs.

Conductivity and its source receipts now govern the actual lumped/Biot admission.
In the admitted uniform-temperature limit it is not an arbitrary extra cooling
knob: the external film coefficient remains an explicit boundary input. Radiation
is explicitly disabled, not assigned a measured zero emissivity. Missing
conductivity never falls back to the demo's value. File mode prints card,
caloric and thermal-transport identities instead of the demo's material constants.

Six additional library regressions and three additional executable regressions
cover source-table equivalence, data-only Cp/density effects, preserved nonlinear
history, interior source conflicts, context/units, missing inputs and real Biot
refusal. Native tests remain unexecuted in the editing environment (Cargo/RCH/DSR
absent). Four independent parsed-table reference trajectories reproduce the
original 48-step heat/creep equations with cumulative defects below 1.87e-14 J;
these are not Rust execution or experimental validation.

```sh
rch exec -- cargo test -p fs-material --test thermovisco_source
rch exec -- cargo test -p fs-couple --example self_heating_creep
```
