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
