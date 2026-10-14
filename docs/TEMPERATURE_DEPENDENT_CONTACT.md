# Source-backed temperature-dependent thermal contacts

The opt-in `fs_conduction::interface::temperature` API closes a missing part
of Material Reality MR20: an interface can change its resistance as the two
bodies heat, instead of using one frozen source sample for the whole run.
Existing fixed-contact and native `.fsim` behavior is unchanged.

`TemperatureContact::new` selects a claim on the original ordered interface
card at an explicit reference query. Every later local evaluation pins that
same claim, retains pressure and all nonthermal coordinates, and changes only
the explicitly named absolute-temperature axis. A source-domain exit refuses;
it cannot silently select a different claim or extrapolate the original one.

Each matching triangle uses the declared closure

```text
T_interface = (mean_face(T_A) + mean_face(T_B)) / 2
R_face = original_selected_claim(T_interface, fixed_context)
q_A_to_B = (T_A - T_B) / R_face
```

The two bulk fields are not required to coincide. For P1 triangles their
arithmetic nodal mean is the exact face mean. Temperature is sampled per face
pair, not once per named surface. The original consistent P1 contact mass
matrix still integrates the jump; this is not a lumped contact flux. Fixed
interfaces may coexist in the same complete contact set.

## Capacity-storage advancement

`TemperatureContacts::advance_capacity` uses the existing backward-Euler
linear or nonlinear-conductivity solver. All contact trials receive the same
physical old temperature. Relaxation changes only the next source-query guess.
A returned temperature must satisfy the full free-node residual

```text
C (T_new - T_old) + dt [K_bulk(T_new) T_new
                       + K_contact(R(T_new)) T_new - b] = 0
```

under a fixed, caller-declared absolute joule tolerance, and the independent
storage-minus-net-input energy tolerance. A small temperature increment or
closed whole-body energy sum is not enough: opposing interface residuals can
cancel in that sum. The accepted field remains the last counted solve; there
is no hidden replacement solve after the physical audit.

The result separates the actually applied contacts from freshly queried
physical endpoint contacts. Both retain source values, pins and local query
receipts. `conduction.contact_fluxes` belongs to the applied operator;
`physical_contacts.frozen_interfaces().fluxes(&conduction.temperature)` gives
the actual source-law flux. Contact transfer is internal and never added to
external heat. Cumulative Krylov, Newton and backtracking work includes every
contact trial, not just the final inner solve.

Admitted exterior data are prescribed Neumann/Robin and endpoint Dirichlet
values. Nonlinear bulk conductivity uses the existing Newton policy. Source
failure, cancellation, exhausted work or failed physical gates returns no new
state. The caller advances its history only after success.

## Checks and limits

`cargo test -p fs-conduction --test temperature_contact` covers G0 source and
quantity admission, G1 eight-node spatial references with unequal capacities
and unequal timesteps, nonlinear bulk conductivity, G3 frozen-law rejection,
G4 cancellation and failure without history mutation, and G5 replay.
The numerical reference uses independently assembled tetrahedral stiffness,
lumped capacity and consistent contact blocks. A frozen-law first endpoint
closes whole-body energy to about 1e-14 J but misses the full physical residual
by about 0.21 J; the new algorithm rejects it. These equation checks ran in
Python; the Rust tests have not run in the implementation environment because
Cargo/rustc are unavailable.

This is an Estimated fixed-geometry model with a declared interface-temperature
closure, not experimental validation or a continuum-error certificate. It
holds pressure, gap, surface condition, medium and contact pairing fixed.
Frozen snapshots do not expose the missing dR/dT term and must not be used as
complete nonlinear-contact adjoints. Use the explicit full-contact derivative
API below for capacity storage. This capacity forward entry point does not
close airflow or radiation feedback. Native project lowering, nonmatching
dynamic contacts, contact mechanics and changing source selection are absent.

## Heterogeneous phase-change storage

`TemperatureContacts::advance_enthalpy` advances the same contact law through
`HeterogeneousEnthalpyBackwardEuler`. Every node retains its assigned chart
and frozen reference mass. The primary balance is

```text
m_i * (h_new_i - h_old_i) + dt * (transport_i(T(h_new), fraction(h_new)) - load_i) = 0
```

All nested trials start from the identical old specific enthalpy. Only the
contact-query temperature guess is relaxed; never invert temperature to obtain
enthalpy on a plateau. Two equal-temperature states can contain different
latent energy, and the full physical residual checks that distinction.
Nonlinear conductivity and declared phase-conductivity laws remain owned by
the existing enthalpy solver. Source pins and nonthermal query coordinates
stay fixed just as on the capacity path.

`advance_enthalpy_with_ambient_radiation` additionally composes the existing
area-mean gray-patch radiation solver. Original convection, signed radiative
heat and internal contact transfer remain separate. Each outer contact trial
contains a complete radiative solve from the same old h. Before acceptance,
the audit re-evaluates both actual contact resistance and nonlinear radiation
at the returned field. Neither the outer residual target nor the energy gate
is replaced by a temporarily frozen inner problem's tolerance.

`ContactEnthalpyStep` retains the final inner conduction state, applied and
physical contacts, and independent `physical` balances. Optional `radiation`
and `combined_boundary` describe the last inner radiative response. Top-level
`radiation_iterations`, `nonlinear_iterations` and `krylov_iterations` count
all nested contact trials. Their checked allowances include the full product
of contact, radiation, Newton and Krylov limits. Failures publish no new h.

### Executable example

The two-tetrahedron example explicitly declares different phase charts,
reference densities 1 and 2 kg/m3, nonlinear conductivity, 10 W of total
heating, a source-backed synthetic R(T) curve and emissivity 0.2. The cold
material melts at 310 K. All inputs are numerical fixtures, not validated
properties of named substances. Both command arguments are required:

```sh
cargo run -p fs-conduction --example temperature_contact_phase -- 300 1
cargo run -p fs-conduction --example temperature_contact_phase -- 340 1
cargo run -p fs-conduction --example temperature_contact_phase -- 340 2
```

The first argument is reservoir temperature in kelvin; the second scales the
whole synthetic contact-resistance source curve. CSV output contains every
accepted node's temperature, specific enthalpy and liquid fraction, plus the
contact resistance, step energy, signed radiation and actual iteration counts.
Three unequal timesteps (0.15, 0.35 and 0.5 s) sum to one second. Internal
contact heat is never counted as another source. The energy reference is
10 J minus the integrated signed radiation, not always 10 J.

### Phase-change verification and boundaries

```sh
cargo test -p fs-conduction --test temperature_contact_enthalpy
```

Four tests cover independent nonuniform P1 roots, contact-driven melting,
unequal reference masses, stationary but distinct latent histories, source
and work refusal, cancellation, replay, radiative heating/cooling and nested
work accounting. An independent Python monolithic-root comparison against
partitioned contact iterations covered 15 nonradiative and six radiative
endpoints; the largest specific-enthalpy difference was 6.87e-8 J/kg. The
radiative reference energy defect was below 3.5e-12 J. These are mathematical
reference checks, not execution of the Rust code. Native tests and example
compilation remain unverified in the implementation environment.

The phase-change path is a thermal, stationary reference-mass model. It does
not move liquid, change geometry/contact mechanics or infer pressure-dependent
phase equilibrium. Enthalpy Dirichlet conditions remain unsupported. Dynamic
contact is not yet lowered from native `.fsim` projects. The complete dynamic
contact derivative is implemented below for capacity storage, not yet for the
enthalpy/radiative trajectory. A frozen-contact adjoint is not its substitute.
Airflow feedback remains separate for this opt-in contact driver.

## Complete contact-law derivatives and capacity adjoints

`TemperatureContacts::linearize_temperature` evaluates the local source laws
and returns `ContactTemperatureLinearization`, with actual contact loads and
both Jacobian actions. For one face, let B be its consistent P1 mass matrix,
D map the full field to the matched A-minus-B trace, and m average the six
trace nodes. With `q = D^T B D T / R(m T)`, the exact temperature Jacobian is

```text
J_contact = D^T B D / R - (R_prime / R) q m
```

The original matching owner supplies D and B. There is no nearest-face search,
new quadrature, finite-difference PDE derivative or mean-jump replacement.
Even when a face's integrated mean jump is zero, its nonuniform nodal jump
can make the second term nonzero. Its transpose is generally different from
its forward action. Both actions preserve internal heat conservation.

Source slope evaluation uses the selected claim's exact piecewise-linear
segments. Queries at adjacent representable temperatures establish the
original pin's two-sided support; they do not estimate a numerical slope.
Continuously supported constant scalars have slope zero. Redundant equal-slope
knots are admitted, but slope-changing knots, support endpoints and unsupported
interpolation refuse derivatives even when the primal remains valid.

`BackwardEuler::linearize_temperature_contact_endpoint` binds an actual
accepted capacity endpoint, its immutable old temperature and the original
absolute physical residual/energy gates. It rechecks that field with freshly
queried resistance, then attaches the second contact term to the existing
`C/dt + J_bulk + K_contact` operator. The complete nonsymmetric transpose uses
the existing bounded FGMRES owner, including for temperature-independent bulk
conductivity. No primal solve or replacement endpoint is introduced.

The returned `ContactCapacityLinearization` exposes the existing source-density,
source-scale, history, capacity and prescribed-Robin derivative operations.
`pullback_with_resistance_controls` adds one independent multiplier for each
interface's entire resistance law, using the same complete physical dual.
Source pins, geometry, pressure and all nonthermal context stay fixed. The
reported contact control is at multiplier one; equal source cards do not merge
named interfaces. These watt-coordinate contractions need no extra dt. Reverse
a trajectory by applying the history pullback to that same dual at each saved
endpoint, summing only controls that precede the selected observation.

```sh
cargo test -p fs-conduction --test temperature_contact_tangent
```

Five new tests cover fresh-query Jacobian differences, transpose duality,
orientation/replay, zero-mean nonuniform jumps, domain/knot/work/cancellation
refusals, and three-step trajectory derivatives against fresh initial states,
sources, capacities and immutable resistance cards. The independent Python
P1 reference gives 0.614056765497 for a combined control direction versus
0.614056766153 from fresh solves; freezing contact instead gives 0.916052407587.
The nonlinear-conductivity counterpart agrees within 6.0e-10. These are
independent mathematical checks, not native Rust test results.

At implementation time GitHub's native job completed dependency bootstrap but
Cargo rejected a pre-existing package cycle through fs-feec, fs-couple,
fs-aeroac, fs-lbm and fs-conduction before compilation. Local Rust tools are
also unavailable. The native tests therefore remain unexecuted, not passed.
These local fixed-grid derivatives do not establish physical validation,
invertibility, global smoothness, geometry sensitivity, or an enthalpy/radiative
contact adjoint. The explicit face allowance includes all fixed and variable
contacts and checks sparse staging sizes; it is not a whole-process RSS limit.
