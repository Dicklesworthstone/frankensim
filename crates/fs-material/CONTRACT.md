# fs-material — CONTRACT

> **Layout (restored 2026-11-09).** This contract is an index plus part files. Once this file grew
> past 1 MiB (1,048,576 bytes), some commits replaced the whole contract with their one new
> section; no edit below that size ever did. The cause has not been identified: the cutting commits
> are ordinary multi-file commits, so a tool that cannot read past 1 MiB and then writes the file
> whole is the likely kind of cause.
> This file was cut to a single section 8 times between 2026-11-05 and 2026-11-09; the last cut (2026-11-09 12:27) removed 1.7 MB.
>
> - All earlier contract text is in `contracts/part-NN.md`, unchanged. A higher part number is newer;
>   inside a part, sections are newest first, as they were in the single file.
> - Add a new section under "New sections" below, newest first.
> - There is no size limit on this file or on the parts (owner ruling, 2026-11-09). If a tool shows
>   you this file empty or cut short, do not write it back whole: add your section with a patch, and
>   check that `git diff --stat` shows no large deletion here before you commit.
> - The ten required sections are in the parts named under their headings below.

## New sections (newest first)

## Identity and state words of a thermochemical gas parcel

For the durable snapshots of fs-gas's cavity, the gas model and its closed
parcel write their own words (`thermochemical_gas.rs` and
`thermochemical_gas/energy.rs`, forwarded by `ClosedGasParcel`), each struct
by a pattern with no rest.

- `ThermochemicalGasModel::identity_words`: the mixture (fs-thermochem's
  words), each species' transport card by its `content_hash()` with its
  selection policy (a pinned plan by each property name and claim hash), the
  water component, and every pack receipt with its four fields.
- `IsochoricGasParcel::identity_words` names a parcel as a reference: the
  model, every axis of its point with its value and its quantity, its mass and
  its volume. The resolved state is not written: `model.evaluate(point)` is
  the one place it comes from. The point is fs-matdb's and is read through
  `axes` and `axis_quantities`, which return its two maps.
- `state_words` writes where a closed parcel stands on its reference: its
  volume, its temperature and its pressure, three words. In a closed parcel
  the model, the mass and the other axes never change, and
  `stored_energy_change_to` refuses a state where they differ.
- `at_volume_temperature(V, T)` is the parcel that stands there. The model,
  the mass, the specific gas constant and the other axes are the reference's;
  the volume and the temperature are the caller's. The pressure is not given:
  it is `(mass / V) r T`, the expression and the order in which
  `discrete_gradient_heat_volume` builds a step's endpoint, and the resolved
  state is evaluated at that point. Nothing is solved. As in that step, a
  resolved density that is not `mass / V` to 64 eps is refused. A caller that
  holds three words writes the rebuilt parcel again and compares, as fs-gas
  does: a pressure that is not the one this computes is then refused.

Tested in `crates/fs-gas/tests/component.rs`, where the parcel builders are.
The identity words are equal for a parcel built again and differ for another
NASA-9 coefficient, formation datum, transport card, selection policy, water
component and state. A parcel two steps from its reference is rebuilt from
its volume and its temperature equal in every field, its pressure among
them. No volume, a temperature that is not one, and a temperature and a
pressure the sources do not cover are refused. Not tested: two models that
differ only in a pack receipt (a model assembled from cards has none), a
pinned selection, a point with other axes than the builders' one, and the
density refusal: the pressure is computed from that density, so no volume
and temperature given here reach it.

## A card's Newtonian liquid in a laminar pipe

`fluid_transport::laminar_pipe::MaterialLaminarPipe` is the Hagen-Poiseuille
resistance of fully developed laminar flow in a circular pipe of fixed radius
`r` and length `L`, `Q = G (p_a - p_b)` with `G = pi r^4 / (8 mu L)`, whose
dynamic viscosity `mu` is a material card's at the temperature it is asked
at. It is prepared from a `ResolvedMaterialStatePoint` holding either a
dynamic viscosity alone, or a kinematic viscosity and a density (then
`mu = rho nu`; `ViscosityBasis` says which). Any other property, a dynamic
viscosity with a density and a kinematic one without are refused, all as
`MaterialStatePointError::InvalidDerived` (no new error type).

Preparation takes the temperature axis and a declared range that must hold
the original point's positive, absolute temperature. Every other coordinate
must be a fixed external condition of the claims, with a dimensionless
provenance flag bounded to one value held at its query value, as the Ohmic
node and the thermal owner hold it. Each claim's pinned line over the range
is read once, its values checked against the original admissibility, and the
card is resolved at every knot of every line, so a range that leaves a claim,
or a viscosity that is not positive inside it, refuses at preparation and
not in a later step.

`viscosity_and_slope(T)` resolves the SAME pinned claims through the card
again at `T` (`ResolvedMaterialStatePoint::at_query_point`) and returns `mu`
with `d mu/dT`. The slope is the pinned line's on the segment `T` lies in:
the right segment at an interior knot, the left one at the last knot, as the
Ohmic node takes its resistance slope. With a kinematic claim it is
`rho' nu + rho nu'`, each on its own line's segment. `conductance_and_slope`
returns `G` and `-G mu'/mu`. A temperature outside the range, or one the card
refuses, is refused, never extrapolated.

`LaminarPipeComponent`, behind the existing `component` feature, puts the
pipe behind the ports of the Carreau pipe: two `FluidPressureFlux` `EffortIn`
ports `a` and `b` and one heat port. It returns `Q`, `-Q` and `-G dp^2`: all
of the viscous power leaves through the heat port. Each proposal resolves the
card again at the heat port's temperature. The class is `Smooth` within a
segment of the card's lines. The JVP and VJP act with the full 3 x 3
Jacobian, including the temperature column `-Q mu'/mu`, `Q mu'/mu`,
`Q dp mu'/mu`. The receipt books `dt P` dissipated and `dt P / T` entropy,
with the kernel's `64 eps` energy tolerance. The component's error is the
pipe's `MaterialStatePointError`. Values that are not finite, a timestep that
is not positive, a wrong slot count and cancellation are refused before any
output is written, each as `InvalidDerived` naming what was refused. The
identity binds the card's content hash, each selected property's receipt
hash, the temperature axis, the reference temperature, the range, the
geometry and the basis. The codec writes no state words.

`tests/laminar_pipe_component.rs` checks, on manufactured tables (no
measured liquid claimed):
- a dynamic-viscosity card and a kinematic one with its density on other
  knots, against hand values within 1e-12: outputs, conductance, slope, and
  the temperature and pressure columns of the tangent, at knots of either
  line too; the conformance kit with the tangent check inside segments;
- a pipe on a fixed 6e4 Pa drop through the kernel into a 1 J/K enthalpy
  body, against a hand backward step bisected at the logarithmic mean the
  body's port publishes, within 2e-9 K over 60 steps across the knot (the
  body ends at 360.02 K and the flow has risen 2.07 times), with the boundary
  work, viscous heat and stored energy within 1e-9 J;
- the preparation refusals, a viscosity that reaches zero at a knot inside
  the range among them;
- the proposal refusals, and a kernel step the card cannot close (3e5 Pa for
  5 s), refused with the pipe named and the body kept bit for bit;
- the identity.

This is a quasisteady, incompressible resistance with one temperature for
the whole liquid. It adds no fluid storage, inertia, entrance length,
transition to turbulence (laminar flow is the caller's assumption; no
Reynolds number is checked), wall slip, pressure dependence of the
viscosity, thermal expansion of the pipe, heat carried by the flow, or
propagated viscosity uncertainty.

## A card's Ohmic conductor behind resistor ports

`conductor::component::OhmicResistor`, behind the existing `component`
feature, puts a `MaterialOhmicNode` behind two electrical terminals (`a`, `b`)
and one heat port. In every proposal the node resolves the card's
`R(T) = rho(T) L / A` again at the temperature the heat port receives, so a
connected thermal owner moves the resistance with no pairing-specific loop.
All of the Joule heat `I^2 R` leaves through the heat port, as `-P`. The
state is `()`.

`Terminals::Current` receives the terminal currents (`FlowIn`, the form of a
resistor in a current loop): the current through the conductor is
`i_a - i_b`, and `a` and `b` return `V = I R` and `-V`. `Terminals::Voltage`
receives the terminal voltages (`EffortIn`, a branch between two nodes): the
current is the node's `current_a_for_drive(Voltage(v_a - v_b))`, and `a` and
`b` return it and its negative. The port kinds and roles are those of the
kernel's constant-resistance `Damper` and `Admittance` in the electrical
domain, so a scene can use either.

The class is `Smooth` within a segment of the card's table. The JVP and VJP
act with the full Jacobian, including the temperature column `dR/dT` from
`MaterialOhmicNode::resolve_and_slope`. That new public method (no new
`*_with_*` name: the enumeration ratchet counts them) returns the resolved
conductor with `dR/dT` at one temperature, refuses a
node bound to a liquid mass fraction, and at an interior knot takes the
right segment's slope (the last knot takes the left one's), as the node's
own graph tangent does. The receipt books `dt P` dissipated and `dt P / T`
entropy, with the kernel's `64 eps` energy tolerance. The component's error
is the node's own `ConductorError` (no new error type). A temperature outside
the node's admitted interval is the node's refusal, unchanged, never an
extrapolation. Values that are not finite, a timestep that is not positive,
a wrong slot count and cancellation are refused before any output is
written, each as `ConductorError::InvalidInput` naming what was refused. A
node with a liquid mass fraction axis is refused at construction: the heat
port carries a temperature only. The identity binds
the card and receipt content hashes, the temperature axis, the reference
temperature and interval, the geometry, any free expansion, and the terminal
form. The codec writes no state words.

`tests/ohmic_component.rs` checks, on a manufactured three-knot table (no
measured material claimed):
- both forms against hand values within 1e-12, and the current form against
  the node's graph evaluation bit for bit;
- the hand temperature column at the interior knot, and the conformance kit
  with the tangent check inside segments;
- a heater on 1 V through the kernel into a 1 J/K enthalpy body, against a
  hand backward step bisected at the logarithmic mean the body's port
  publishes, within 2e-9 K over 60 steps across the knot, with the boundary
  work, Joule heat and stored energy within 1e-9 J;
- a current loop whose inductance's lost energy is the body's gain within
  1e-9 J;
- the refusals, and a kernel step the card cannot close (10 V for 5 s),
  refused with the wire named and the body kept bit for bit;
- the identity.

The first draft of that kernel refusal drove 3 V, which the rising
resistance closes near 370 K; the kernel took it, and the test was
corrected.

**Provenance flags are held.** Before this section, a conductor refused
every seed card that has a resistivity, because each one bounds its
resistivity's validity with provenance flags (`engineering_reference_model =
1`, `source_specimens_matched = 0`, `annealed_material = 1` on Nickel 200),
and `ScalarConductorResponse::admit` and the Ohmic node's own check counted
those flags as unowned coordinates. Both now hold a validity coordinate
that the claim's value does not depend on, that the validity declares
exactly dimensionless and bounds to one value, at the finite value the query
gives it (fs-matdb has already checked it against the bound). The thermal
owner (`ThermalMaterialPoint::admit_static`) already held such flags the same
way. An untyped bound has no known dimension and still refuses; so does a
dimensionless range, such as a phase fraction over [0, 1], and every value
axis that is not owned or a typed fixed external condition. The rule is
`ResolvedScalarProperty::is_held_flag`, private to the crate, and the
crate-private variants `fixed_external_support` and
`owned_coordinate_support` take it as a switch; the public predicates and
every other caller keep the strict rule, and so does the phase-fraction
path. The added test in `tests/ohmic_component.rs` holds a dimensionless flag
with the unflagged card's values to the bit, and refuses an untyped flag and
a dimensionless range.

This is memoryless DC conduction of fixed geometry, at nominal values. It
adds no inductance, capacitance, skin effect, contact resistance, heat
storage of its own, or propagated resistivity uncertainty. The node's other
limits are unchanged.

## A thermal path along a named absolute-temperature axis

`state_point::thermal_path_axis(point)` names the coordinate a thermal path
runs along. That is the legacy untyped `T` when the point has one. Otherwise
it is the point's one coordinate declared a static absolute temperature,
whatever its name; the seed cards name theirs `temperature`. It is `None`
when the point has neither, or has more than one such coordinate.
`integrate_isotropic_thermal_expansion` and the orthotropic path use it both
to read the two temperatures and as the abscissa of the expansion curve.
A point that has `T` behaves exactly as before. On a named axis, a current
point that retypes the coordinate as anything but an absolute temperature is
refused, as is a change of any other coordinate. Before this, the path looked
only for `T`, so it refused every card whose temperature axis has another
name. Among the seed manifests, 463 name it `temperature` and 16 name it `T`.
fs-mbd-material's tests/thermal.rs checks the named axis against `T`, bit for
bit, together with the refusals.

## Temperature-dependent layer conduction on the composition contract

`thermal::conduction::component::LayerConductor`, behind the existing
`component` feature, makes the original `TwoLayerConduction` available on
two heat ports. It owns both conductivity curves and the admitted
`LayerGeometry`, so it can be moved into `fs-system::SystemBuilder` without
leaking material data or borrowing an application-owned chart. A proposal
borrows those immutable curves only for the native layer solve. Original
knots, temperature support, planar/cylindrical/spherical resistance metrics
and area-specific contact resistance are retained. No conductivity is frozen
at the start of a simulation, and no second flux or interface solver is added.

Ports `a` and `b` receive their owners' temperatures [K] and return `Q`, `-Q`
[W], positive into the conductor. The state is `()`. The native solve checks
each material's own bulk-to-interface temperature span; a resistive contact
can bridge disjoint spans, while an unsupported contact face is refused.
The full `Smooth` JVP and VJP use the native analytic endpoint derivatives.
The Jacobian is generally nonsymmetric when conductivity varies with
temperature; the reverse action is its transpose, not the forward action.

Energy storage, physical dissipation and numerical loss are zero. The
component generates `dt Q (1/T_b - 1/T_a)` of entropy, computed without a
subtraction of nearly equal reciprocal temperatures. Its receipt retains
the native heat-continuity residual and tolerance; contact-jump acceptance
remains in the native solve. No heat or entropy is fitted to an expected
trajectory. Proposals and tangent actions poll cancellation, keep the
original bounded scalar iteration, allocate no heap on success, and leave
outputs unchanged on refusal. Errors remain `DepositionError`.

The native Component snapshot codec writes no state words. Its parameter
identity retains both complete ordered conductivity knot sequences, the
original geometry kind and coordinates as well as its resistance metrics,
contact resistance, and all local iteration/tolerance controls. Equivalent
planar resistance coordinates do not erase the original radial geometry.

`tests/thermal_conduction_component.rs` checks the mandatory conformance kit,
constant-conductivity series/contact resistance and reversal, the independent
Kirchhoff integral for affine conductivity, both tangent actions and their
finite differences, disjoint supports, and input/domain/cancellation refusals.
Direct codec tests retain both charts, original geometry and solve controls.
Two owned enthalpy bodies and this conductor run on the actual generic kernel
against an independent scalar caloric/integrated-conductivity oracle; their
energy, entropy and temperature feedback are checked. A source-domain exit
must preserve accepted state, time, ledger and published port values.

This remains source-free, quasisteady conduction in fixed geometry. It adds
no transient layer storage, convection, internal heat generation, evolving
phase distribution, contact-resistance temperature law, uncertainty model or
experimental authority. Conductivity curves cloned from a resolved card keep
their numbers and support; card provenance is retained by the preparing
caller, not inferred by this adapter. The thermal owners choose which
endpoint, midpoint or entropy-conjugate temperatures their heat ports return.

## Thermal Carreau pipes on the composition contract

`fluid_transport::rheology::pipe::component::CarreauPipeComponent`, behind
the existing `component` feature, exposes the original `LaminarCarreauPipe`
through two pressure/volume-flow terminals (`a`, `b`) and one heat port.
Inputs are the two pressures and the uniform material temperature; outputs
are the original signed radial flow, its opposite, and negative original
viscous power. The state is `()`. Its temperature comes from the connected
thermal owner, and `at_temperature` retains the source law, geometry and
authored physical wall-rate ceiling at every trial. The full reduced-rate
support, including the zero-shear core, must still be admitted after shifting.

The `Smooth` port map and its full three-by-three JVP/VJP use the original
pressure conductance and `CarreauPipeProfile::temperature_flow_derivative_m3_s_k`.
For a vertical shift, `Q_T = -(log a_T)' p Q_p`; for a vertical-and-rate
shift, `Q_T = -(log a_T)' Q`. These identities use the already accepted
radial profile and original analytic shift derivative, with no perturbed
source solves or additional quadrature. Both pressure terminals and the
heat row participate in the tangent. Equal pressure offsets change nothing.

Physical dissipation leaves as heat once, and entropy generated is that heat
over the supplied temperature. The native reconstructed pressure residual
remains the receipt's residual and contributes its measured work defect to
the energy allowance; the adapter neither replaces native heat by requested
pressure times flow nor books that difference as numerical loss. Endpoint
flow is `Untold`, since this algebraic response has no midpoint flow history.
Success uses stack storage. The native stress inverse retains its hard work
bound, quadrature polls cancellation at every sample, and any refusal leaves
caller outputs unchanged. The error remains `MaterialError`.

The native Component snapshot codec writes no state words. Its parameter
identity binds the five Carreau-Yasuda parameters and the rate support, the
shift model's kind and coefficients, its reference temperature and validity
band, the shift mode, the radius, length and wall-rate ceiling, and all four
local solve and quadrature controls. The temperature the pipe was prepared at
is not bound: every trial rebinds the law to the heat port's temperature.

`tests/carreau_pipe_component.rs` exercises Poiseuille limits, reversal and
pressure-reference invariance; an independent analytic stress-inverse radial
oracle for a nonlinear Carreau law; both thermal shift conventions and the
complete analytic tangent; source residuals; work/support/cancellation
refusals; and two native pipes, a declared hydraulic compliance and a shared
`UniformEnthalpyBody` on the actual `fs-system` kernel. The R-C-R fixture
compares pressure, thermal storage, delivered work and viscous heat over
60 steps against independent elimination of Poiseuille flow, midpoint
compliance and logarithmic mean temperature. Two codec tests check the empty
state and that a one-ulp change of each bound real parameter, or a change by
one of each integer control, changes the identity. `fs-system` is a development dependency only; production
composition uses `fs-ports`.

This is fully developed, incompressible, quasisteady no-slip resistance at
uniform temperature. It supplies no fluid inertia, entrance or turbulent
correction, cavitation, pressure-dependent constitutive law, advected
enthalpy or species, or thermal spatial resolution. Its wall-rate ceiling
remains an authored envelope, not a Reynolds-number transition certificate.

## Tabulated hardening: a roundoff trial on the yield surface

`TabulatedHardening::solve_increment` returns no increment when the trial
exceeds the yield value by at most 32 machine epsilons of it and the exact
return is below the resolution of the accumulated plastic strain
(`alpha + delta == alpha`): the nearest representable history is the
accepted one. An accepted state re-evaluated at its own strain assembles its
trial stress apart from the yield value and can land that close above it;
the law had refused its own state there ("plastic consistency return is not
representable"). A larger unrepresentable return is still refused, and every
return that solved before keeps its value. `TabulatedJ2` and the tabulated
rate law take a zero increment from a plastic trial as no hardening
increment, the native radial return at the current yield value.

## A uniform enthalpy body on the composition contract

`phase::component::UniformEnthalpyBody`, behind the feature `component`
(optional dependencies `fs-ports` and `fs-exec`; the default build has
neither as a direct dependency and is unchanged), is a fixed mass of the
material of one `EquilibriumEnthalpyPhaseCurve` as a `Component` of the
composition contract. Its state is the chart's own `EquilibriumPhaseState`.
It adds no physics: no heat capacity, no apparent latent heat, no transport
law, no geometry. Its error is `PhaseStateError`.

Its one port, `heat`, receives the heat rate [W] into the body and returns
the chart's entropy-conjugate temperature of the step that heat makes,
`thermodynamic_mean_temperature(h0, h0 + dt Q / m)`: the logarithmic mean
over an interval of constant heat capacity, the transition temperature on a
latent plateau. The heat over that temperature is the body's entropy change
exactly, so the body generates no entropy and a conductance given this
temperature produces `G (theta_a - theta_b)^2 / (theta_a theta_b) dt` and no
more. The state reached is `state_at_specific_enthalpy` of the new enthalpy.
The receipt stores `m (h1 - h0)` as represented and books the remainder of
the heat as numerical loss; heat too small to move the enthalpy is not
refused. Stored energy is `m h` in the chart's datum.

Its class is `Monotone1d`: one input, one output, and the mean never falls
as the heat rises (it is flat across a plateau). The tangent is the chart's
derivative of the mean in the final enthalpy times `dt / m`. That derivative
is continuous in the heat even where a step ends on a corner of the chart;
only at no heat on a corner is it the chart's one-sided choice (half the
right derivative, the left one at the upper end of the chart). The class has
a cost today: fs-system differences a component without a tangent only if
its class is smooth or exactly linear, so in a loop this body must be added
with the tangent it offers.

A state that does not move (no heat, or heat too small to change the
enthalpy) is returned as it is, sign of a zero enthalpy and identity
included, as the chart's own `advance_specific_energy` returns it.

The native Component snapshot codec stores the exact specific-enthalpy bits
and resolves the complete state through the same original curve. Encoding
checks the reconstructed native state identity; decoding refuses malformed
word counts, nonfinite or unsupported enthalpy, and unrepresentable stored
energy. The parameter identity binds mass and the complete curve identity,
including its original card, knots, caloric/density laws and retained source
context. The codec adds no independent thermal history or replacement law.

Refused, with candidate and outputs untouched: a mass that is not positive
and finite; widths other than one; a heat rate that is not a number
(`NonFiniteSpecificEnergyIncrement`); a step that is not positive and
finite; a state of another chart (`StateCurveMismatch`); a step that leaves
the chart (`OutsideEnthalpyDomain`, the chart's own); a mass so large that
`m h` is not a finite number of joules (both the stored energy and the
step); cancellation (`EntropyEvaluationCancelled`). It does not shorten a
step.

Not here: the body's volume, surface and density, a conductance law, heat
generated inside the body, expansion work, and any limit on the step for
accuracy. The enumeration ratchet is unchanged: no new `*_with_*` function,
`*Error` or `*Candidate` type, or trait.

Tests. `tests/phase_temperature_derivative.rs`, module `body`, 4 tests with
the feature (`cargo test -p fs-material --features component --test
phase_temperature_derivative`): the conformance kit and the tangent check
on ten steps of a melting chart (within the solid; ending on and starting
on each of the two corners of the plateau; onto it, along it and across it;
back down across both), with state and temperature equal to the chart's own
answers to the bit, heat over temperature equal to the chart's entropy
change to 1e-13, the slope equal to the chart's derivative to the bit; the
mean rising with nine heats from one state; the one-sided slope at rest on
each corner, confirmed from the right by a difference and exactly zero from
the left; a heat of 1e-30 W taken and booked as lost; a state at an enthalpy
of minus zero kept as it is under no heat; and every refusal above.
Codec tests round-trip every native state field through latent and original
affine-capacity charts, and check malformed/foreign states and parameter identity.
`crates/fs-couple/tests/gas_piston_graph.rs`, 3 tests: the body, the
kernel's `Conductance` and the fs-gas cavity as a declared graph against the
piston module's own finite-body step
(`GasPiston::thermal_discrete_gradient_candidate`) on that module's
fixtures, a step at a time, with every bound computed from the two solves'
tolerances before the run. Per step the module's heat may be off by its
1e-10 J source tolerance (as much again is allowed for the graph) and each
side's gas energy by the owner's 1e-7 J/kg. The exact step does not amplify
a difference between two runs (contact damps the difference of the two
stores' energies and keeps their sum), so after n steps either store's
energy differs by at most n errors of both kinds; with a moving wall the two
first-order couplings between volume and gas temperature are added.
Measured against allowed:

| fixture | steps | body | gas | heat a step | entropy |
|---|---|---|---|---|---|
| fixed walls, body from 300 K, 4000 W/K | 12 | 2.3e-11 K (8.0e-9) | 2.1e-11 K (7.8e-9) | 1.2e-11 J (3.4e-9) | 4.7e-15 J/K (4.9e-10) |
| fixed walls, body from 650 K, 4000 W/K | 12 | 3.3e-11 K (8.0e-9) | 2.9e-11 K (7.8e-9) | 1.7e-11 J (3.4e-9) | 5.8e-16 J/K (6.7e-10) |
| wall thrown at -2 and 2 m/s, 2000 W/K | 12 | 1.1e-10 K (8.0e-9) | 8.0e-11 K (7.6e-9) | 5.2e-11 J (1.8e-9) | 6.9e-15 J/K (5.7e-10) |
| melting through a plateau at 350 K, 2000 W/K | 40 | 4.1e-8 J/kg (1.3e-5) | 6.5e-11 K (2.3e-8) | 4.1e-11 J (5.5e-9) | 3.3e-15 J/K (2.2e-9) |
| freezing through a plateau at 500 K, 2000 W/K | 40 | 1.0e-7 J/kg (1.3e-5) | 1.5e-10 K (2.3e-8) | 7.5e-11 J (6.5e-9) | 1.7e-15 J/K (2.6e-9) |

With the wall moving, its velocity agrees to 9.4e-13 m/s (1.3e-10) and the
gas volume to 3.0e-18 m3 (8.1e-16). The states above do not depend on the
tangents the graph is given; the iterations do, and each run is held to
three a step (measured 27 and 28 in 12 steps, 25 behind the moving wall, 77
and 69 in 40). An independent reader's copies with the body's slope halved
or doubled (and, on the fixed sensible run, zeroed) still met every bound
above, in 107 to 339 iterations.
These are not small steps: in its twelve
the body from 650 K gives 66 J and cools by 132 K, and the bodies on a
plateau leave it, one molten at 351.5 K and one frozen at 464.8 K. With the
kernel's ideal `HeatCapacity` in the body's
place, which returns the arithmetic mean of its end temperatures, the same
comparison is off by 0.065 K and 0.21 K in the body, 0.032 J and 0.10 J in a
step's heat and 0.19 K and 0.61 K in the temperature the contact was given:
more than a million times the bounds, as the test requires.

In those comparisons the conductance is constant (the module's own latent
test makes it depend on the liquid fraction, which the kernel's element does
not). Three of the six runs use the module fixture's heat bracket of 100 J
either way: both moving runs and the fixed run from 650 K. The other three
have a narrower one because the module, alone, refuses a step with its own,
which the test checks and prints before it narrows anything: from 300 K
(60 J used) it refuses its second step and from the plateau at 350 K (20 J)
its tenth, both `OutsideBracket` from the gas's own temperature solve at the
far end of the bracket; from the plateau at 500 K (20 J) it refuses its
sixth, `OutsideEnthalpyDomain` from the body's chart at -4,213 J/kg. (The
first landing said all three kinds of run needed a narrower bracket, and its
correction that the fixed sensible and the latent runs did for the gas's
sake; both generalised from the first failure of a loop.) The module's
answer depends a little on its bracket: with 60 J the moving row read
6.3e-11 K, 4.9e-11 K and 3.1e-11 J, and the 650 K row 5.2e-11 K, 4.6e-11 K
and 2.6e-11 J.

What these comparisons cannot see: an error inside the chart's mean or the
gas owner's step, which both sides call; an error in how the step length or
the body's mass enters that is right at 5e-5 s and 1e-3 kg, the only values
used; heat that rounds away, a step that starts on a corner, a heat that
reverses, and a chart made from heat-capacity data. The module's midpoint
thermal step, a conductance that depends on state, and its real-gas and
reactive thermal variants are not compared.

## Causal pressure-driven molecular compressibility

`RelaxationState::compressibility_relaxation` retains the original evaluated
species/reservoir state and lowers its complete `Cp(s)` to a positive finite
pressure-memory law `beta(s) = beta_inf + sum delta_beta/(1+s*tau)`.
The bounded original symmetric Jacobi owner diagonalizes the rank-one inverse
of `Cp_inf + sum C_j/(1+s*tau_j)`. These pressure-driven times differ from the
caloric and compression-driven times. Equal source times combine exactly;
zero-strength reservoirs remain in the retained source but allocate no arm.
No frequency fit, attenuation conversion or extra bulk viscosity is introduced.

Each arm obeys `tau*q_dot + q = p`, stores `delta_beta*q^2/2` per volume and
dissipates `delta_beta*(p-q)^2/tau`. A spatial pressure consumer must use the
frozen compressibility for its instantaneous mass; the memories supply the
relaxed increment once. Original gas coherence, finite positivity, eigenpairs
and relaxed/frozen limits are checked before accepting a realization.
Harmonic response uses the existing peak `exp(-i*omega*t)` convention.

Four focused regressions compare the independent one-reservoir pressure time
and strength, multi-reservoir response against the reciprocal original caloric
modulus over frequencies from zero to `f64::MAX`, equal/zero/empty reservoirs,
physical pressure scaling and original-state admission. They ran in a bounded
native build of unchanged original gas/caloric/Jacobi numerical owners; the
largest sampled relative compressibility discrepancy was `2.144e-16`.
This is a frozen-state, linear molecular constitutive realization; it supplies
no shear viscosity, Fourier heat conduction, material calibration or full
workspace-build claim.

## Causal Kelvin/Maxwell plane-stress ramps

`PreparedTensorKelvinPlaneStress::advance` now accepts arbitrary affine
engineering-strain ramps, holds and reversals with positive finite duration.
The existing constrained QR/SVD preparation and native `PronyRamp` owner evolve
the original anisotropic Maxwell history. Total transverse stress, including
Kelvin viscosity, remains zero; out-of-plane strains and rates remain physical
state. Zero/rank-one/rank-two transverse Kelvin operators use the existing
algebraic traction constraints without a viscosity floor or branch reduction.
Initial incompatible algebraic spring traction refuses without projection.
Strain and branch memory remain continuous across prescribed-rate changes.

Full original stress reconstruction supplies endpoint and interval-mean
reactions and signed work. Independent modal and residual-Gram heat remains
nonnegative. Original equilibrium/branch coordinates supply a separately
evaluated local storage change, retaining sub-roundoff held-state release.
Conservative sums of monotone modal excursions admit the complete full-tensor
strain path; this is a numerical sufficient bound, not an interval certificate.
All failures and cancellation leave the accepted history unchanged.

`ramp_tangent` returns analytical full-six-component endpoint strain, stress,
mean stress and rate derivatives with respect to the three final engineering
strains at fixed incoming history, coefficients and duration. It includes
normal relaxation and viscous rate dependence without perturbed solves.
Both operations reuse preparation in O(n²) work with no additional dense
factorization, n = 3 + 6 times the original branch count. They do not evolve
temperature, fit branches, extend a finite-band calibration claim to arbitrary
ramps, infer finite-strain behavior or admit instantaneous strain jumps.

Four focused regression groups use the existing independent original ODE/DAE
test seam: ramps/holds/reversals in every ordered material plane, full endpoint
and mean finite-difference tangents, subdivision/tiny-hold/algebraic limits,
and atomic domain/source/cancellation refusals.

## Complexation-controlled mineral solubility

`osmosis::solubility::electrolyte::complexation::ComplexedElectrolyteSolubilitySolution`
binds the original source-qualified `CopperSulfidePrecipitation` and its complete
metal/ligand chemistry. Solvent-volume changes and explicit analytical dosing
re-solve the original coupled Cu/S/CN phase, ligand and charge balances. Complete
metal, ligand, other acid-family and spectator inventories survive independently
of the actual dissolved report and retained Cu2S formula stock.

The shared dissolved material view includes H+, OH-, acid-family forms, charged
metal members and neutral metal complexes exactly once. Charged members already
occur in the native aqueous spectator report; only neutral members are appended.
Bound ligand units and solid formula units are not additional osmotic particles.
Every update preserves the original species axis, source receipts, physical
volume/activity conversion and bounded cancellation semantics.

This is the selected fixed-temperature, fixed-pressure ideal dilute source
model, including its admitted pH domain and pure-water ion background. It does
not supply concentrated-electrolyte activities, omitted complexes, competing
solids, heat or kinetic rates. Neutral-complex particle counting, dilution,
analytical dosing, source refusal and transactional cancellation are covered by
focused native regressions using explicitly authored Gibbs test data.

## Parts

| Part | Bytes | `##` sections | First (newest) section | Last (oldest) section |
|---|---|---|---|---|
| [`contracts/part-04.md`](contracts/part-04.md) | 435615 | 108 | Original strain domains in checked solid consumers | No-claim boundaries |
| [`contracts/part-03.md`](contracts/part-03.md) | 447617 | 94 | Authority-separated multi-case identifiability schema (bead I10.1) | Radial material conduction and finite thermal feedback |
| [`contracts/part-02.md`](contracts/part-02.md) | 450128 | 89 | Implicit spatial freezing with original source tangents | Concentric annular Carreau transport |
| [`contracts/part-01.md`](contracts/part-01.md) | 363756 | 116 | Physical viscous stress and shear traction | Original equilibrium factors for mechanical projection consumers |

## Purpose and layer

Full text: [`contracts/part-04.md`](contracts/part-04.md), under the same heading.

## Public types and semantics

Full text: [`contracts/part-04.md`](contracts/part-04.md), under the same heading.

## Invariants

Full text: [`contracts/part-04.md`](contracts/part-04.md), under the same heading.

## Error model

Full text: [`contracts/part-04.md`](contracts/part-04.md), under the same heading.

## Determinism class

Full text: [`contracts/part-04.md`](contracts/part-04.md), under the same heading.

## Cancellation behavior

Full text: [`contracts/part-04.md`](contracts/part-04.md), under the same heading.

## Unsafe boundary

Full text: [`contracts/part-04.md`](contracts/part-04.md), under the same heading.

## Feature flags

Full text: [`contracts/part-04.md`](contracts/part-04.md), under the same heading.

## Conformance tests

Full text: [`contracts/part-04.md`](contracts/part-04.md), under the same heading.

## No-claim boundaries

Full text: [`contracts/part-04.md`](contracts/part-04.md), under the same heading.
