# Load, dissipation, heat and material response in one step

`fs_material::thermovisco::thermal::CaloricPoint` binds an existing Maxwell
material history to the existing caloric chart and uniform thermal-owner
protocol. `advance` preserves prescribed-strain behavior; `advance_controlled`
also accepts `MechanicalControl::Stress { stress_pa, tolerance_pa }`.

At every trial temperature, stress control solves the existing affine law
`sigma = C epsilon + b` for the new total engineering strain. It then executes
the original `Point::advance`, checks the actual stress, and supplies physical
dashpot dissipation plus the admitted signed reversible exchange to thermal
transport. The returned temperature changes spring stiffness and relaxation time,
and therefore creep, stress and heat generation.
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

For temperature-independent spring energies and free strain, the existing
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

With an admitted affine common modulus multiplier `S(T)=a+b*T` and affine free
thermal strain `f(T)` with slope `alpha`, the original retained branch strains
define reference spring potential
`Xi=(E_inf*e^2 + sum(E_j*q_j^2))/2`. Mechanical free energy is `psi=S*Xi`,
mechanical entropy is `s=alpha*sigma-b*Xi`, and mechanical internal energy is
`U=psi+T*s`. Native mechanical strain subtracts the actual free-strain change
from the point's own birth temperature; retained branch strains follow that
same strain increment. Positive, negative and zero expansion slopes are allowed.
The original native mechanical update reports coefficient work
`C_psi=(S_new-S_old)*Xi_old` and thermal-eigenstrain work
`W_th=-sigma_new*(f_new-f_old)`; its branch state and time integration are unchanged.
The SAME caloric callback receives

`Q_internal = volume*(D - C_psi - W_th) - delta(volume*T*s)`.

The coefficient-work and entropy terms are reversible exchange and may cool the specimen.
Physical dashpot dissipation `D` stays nonnegative and numerical loss remains
outside thermal storage. `FeedbackStep::applied_heat_j` is this signed total;
`mechanical_energy.dissipated` remains the independent physical loss density.
The report additionally exposes `reversible_heat_j`,
`mechanical_entropy_energy_change_j`, and `mechanical_internal_energy_change_j`.
`CaloricPoint::mechanical_internal_energy_j()` and `mechanical_entropy_j_k()`
read the actual accepted endpoint contributions. Global balances use internal
energy, whereas the original point's `stored_energy_j_m3()` is spring free energy.

Every source/law knot participates in affine admission. A modulus or free-strain
kink refuses even if endpoint values alone suggest an affine line. Endpoint integration remains first order
and need not conserve adiabatic entropy exactly; its error is not extra material
heating. These quantities do not claim a measured constitutive law.

All admissible dashpot dissipation is assumed to thermalize locally in this
rung. The chart must be strictly increasing, single-solid-phase and
constant-reference-density; mass must match reference volume times that density.
Both the common spring multiplier and free thermal strain must be globally affine.
Relaxation and heat capacity may vary with temperature. Non-affine spring/free-strain
laws, phase changes and finite-deformation volume work require additional terms and
are refused. Independent model inputs are not automatically qualified
as the same measured material. A material name never selects missing properties.

### Thermal expansion and positive heat capacity

Expansion changes the fixed-strain heat capacity. With the viscous strains held
fixed, `e_T=q_j,T=-alpha`, so
`psi_TT=alpha^2*S*K0-2*alpha*b*A0`, where `K0=E_inf+sum(E_j)` and
`A0=E_inf*e+sum(E_j*q_j)`. The total volumetric capacity is
`rho*dh/dT-T*psi_TT`. This follows the free-energy derivative convention in
[Abeyaratne, Continuum Mechanics, Vol. II, §8.4](https://web.mit.edu/abeyaratne/Volumes/RCA_Vol_II.pdf);
the present implementation is a bounded uniaxial hereditary model.

Admission uses the minimum actual caloric-chart segment capacity and a sufficient
full-range bound `T_max*(alpha^2*S_max*K0+2*abs(alpha*b)*A_max)`, with
`A_max=strain_limit*(E_inf+2*sum(E_j))`. For native viscous strain `v=e-q`,
backward Euler gives `v_new=(1-a)*e_new+a*v_old`. Birth at zero and bounded
mechanical strain therefore bound every branch by `abs(q)<=2*strain_limit`;
expansion-bearing endpoint histories check that envelope explicitly. The bound
includes all instantaneous branch stiffness and the mixed expansion/modulus
term even when the current material is stress-free. Insufficient capacity
refuses before publishing a coupled state. Zero expansion retains the prior path.

The existing compatible strip/beam solve turns uniform free heating into
extension, constrained heating into stress, and unequal layer expansion into
curvature. Finite capacitor energy can drive that response through the original
Joule and conduction paths. The oscillator likewise responds through its native
force and inertia; its finite-gas pressure is the original fixed-volume EOS,
not a new acoustic or moving-wall pressure-work model. Focused regressions use
independent free/constrained thermal solutions, cross-section equilibrium,
retained Maxwell histories and complete energy balances. Native Rust execution
of these new cases remains unrun in the editing environment.

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
boundary heat, spring free energy, caloric storage, reversible heat, actual
deposited heat, mechanical internal-energy change and mechanical entropy. Its
cumulative energy defect uses internal energy, including the entropy term. It
prints only after the complete bounded trajectory succeeds.

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

The constitutive law remains small-strain generalized Maxwell behavior; it does
not model metal plastic creep or melting. The existing `thermovisco_dynamics`
owner uses the same internal-energy account for inertial motion and finite-gas
heat exchange. The layered-strip and spatial-beam enthalpy couplers integrate
that account over the original two Gauss points in every layer, retaining the
native geometry, compatibility solve and branch histories. Their electrical
heating path deposits signed material exchange beside the original circuit
heat in the same thermal cells. Contact and radiation transport remain with
their existing owners.

These spatial reports expose `mechanical_internal_energy_change_j`,
`mechanical_entropy_energy_change_j` and `reversible_heat_j`. Their deposited
material heat is the signed sum of physical dashpot loss and reversible exchange;
numerical loss stays outside thermal storage. Scheduled-growth wrappers include
the same entropy contribution through the complete interval and birth account.
New material is born stress-free in the native solid owner, with zero branch
history and zero mechanical entropy contribution, while existing histories are
retained. Affine-source regressions exercise continuation after actual births
and rollback when a later birth is invalid.

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

Four additional regressions cover independent signed affine-elastic heat,
heating at held strain, a hereditary loading/hold/unload history with a closed
scalar temperature reference, and source-knot derivative-jump refusal. The
existing cancellation/budget regression also uses affine stiffness. Native
execution remains unrun in this editing environment; syntax parsing and equation
checks do not establish experimental accuracy.

The existing executable now reads either a canonical material pack or an editable
SI table. Geometry, ambient conditions, load schedule, solver and CSV fields stay
the same as the original demo. These two commands differ only in the table's Cp:

```sh
cargo run --release -p fs-couple --example self_heating_creep -- --table crates/fs-couple/examples/self_heating_creep/cp5.txt 2 2000000 4.8 48
cargo run --release -p fs-couple --example self_heating_creep -- --table crates/fs-couple/examples/self_heating_creep/cp10.txt 2 2000000 4.8 48
cargo run --release -p fs-couple --example self_heating_creep -- --pack material.fsmat 1 2 2000000 4.8 48
```

To change the thermal spring response alone, compare `cp5.txt` with
`cp5_affine.txt` through the same executable:

```sh
cargo run --release -p fs-couple --example self_heating_creep -- --table crates/fs-couple/examples/self_heating_creep/cp5.txt 2 2000000 4.8 48
cargo run --release -p fs-couple --example self_heating_creep -- --table crates/fs-couple/examples/self_heating_creep/cp5_affine.txt 2 2000000 4.8 48
```

The second table changes only the common spring multiplier to
`S(T)=1-0.005*(T-300)`, with the same initial stiffness, heat capacity,
relaxation law, density, conductivity, geometry, load and cooling boundary.
The original material and transport owners then solve the altered temperature,
creep and unloading history. This synthetic sensitivity case is not a measured
metal or polymer. The executable regression independently checks the affine
law's `U=2.5*psi/S` and `s=0.005*psi/S`, loading-induced reversible cooling,
unloading heat, and cumulative internal-energy closure with actual convection.
The new Rust regression remains unrun here; an independent equation calculation
checks its reference trajectory, not native execution or experimental accuracy.

To add expansion through material data, use `cp5_expanding.txt` with the same
command arguments:

```sh
cargo run --release -p fs-couple --example self_heating_creep -- --table crates/fs-couple/examples/self_heating_creep/cp5_expanding.txt 2 2000000 4.8 48
```

This table differs from `cp5_affine.txt` only in
`free_strain=0.00002*(T-300)`. The point's birth datum remains 300 K. The existing
load/transport trajectory now includes free expansion and the entropy contribution
`alpha*stress*volume`; no geometry, material name, spectrum or output filter is
changed to obtain the response. The same CSV channels expose the result.
Its regression checks the first-load identity
`Q_rev=-alpha*T_birth*stress*volume+S_prime*T_new*Xi` and the actual internal-energy
balance through subsequent creep and unloading. This table is also synthetic.
The independent equation reference predicts a first-step temperature of about
296.486636 K, compared with 298.729676 K without expansion; these are calculation
references, not executed Rust output or measured-material validation.

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
not independently sourced material evidence. All supplied tables are synthetic.
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
