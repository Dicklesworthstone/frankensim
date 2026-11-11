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
right derivative, the left one at the upper end of the chart).

Refused, with candidate and outputs untouched: a mass that is not positive
and finite; widths other than one; a heat rate that is not a number
(`NonFiniteSpecificEnergyIncrement`); a step that is not positive and
finite; a state of another chart (`StateCurveMismatch`); a step that leaves
the chart (`OutsideEnthalpyDomain`, the chart's own); cancellation
(`EntropyEvaluationCancelled`). It does not shorten a step.

Not here: the body's volume, surface and density, a conductance law, heat
generated inside the body, expansion work, and any limit on the step for
accuracy. The enumeration ratchet is unchanged: no new `*_with_*` function,
`*Error` or `*Candidate` type, or trait.

Tests. `tests/phase_temperature_derivative.rs`, module `body`, 2 tests with
the feature (`cargo test -p fs-material --features component --test
phase_temperature_derivative`): the conformance kit and the tangent check
on nine steps of a melting chart (within the solid; ending on, starting on
and crossing each corner of the plateau; along the plateau; back down across
both), with state and temperature equal to the chart's own answers to the
bit, heat over temperature equal to the chart's entropy change to 1e-13, the
slope equal to the chart's derivative to the bit; the mean rising with nine
heats from one state; the one-sided slope at rest on each corner, confirmed
from the right by a difference and exactly zero from the left; a heat of
1e-30 W taken and booked as lost; and every refusal above.
`crates/fs-couple/tests/gas_piston_graph.rs`, 3 tests: the body, the
kernel's `Conductance` and the fs-gas cavity as a declared graph against the
piston module's own finite-body step
(`GasPiston::thermal_discrete_gradient_candidate`) on that module's
fixtures, a step at a time, with every bound computed from the two solves'
tolerances before the run (the module's 1e-10 J source tolerance and the gas
owner's 1e-7 J/kg, summed over the steps). Measured against allowed:

| fixture | steps | body | gas | heat a step | entropy |
|---|---|---|---|---|---|
| fixed walls, body from 300 K, 4000 W/K | 12 | 2.3e-11 K (4.8e-9) | 2.1e-11 K (7.8e-9) | 1.2e-11 J (2.7e-9) | 4.7e-15 J/K (3.9e-10) |
| fixed walls, body from 650 K, 4000 W/K | 12 | 5.2e-11 K (4.8e-9) | 4.6e-11 K (7.8e-9) | 2.6e-11 J (2.7e-9) | 9.3e-16 J/K (3.9e-10) |
| wall thrown at -2 and 2 m/s, 2000 W/K | 12 | 1.1e-10 K (4.8e-9) | 8.0e-11 K (6.9e-9) | 5.2e-11 J (1.4e-9) | 6.9e-15 J/K (4.3e-10) |
| melting through a plateau at 350 K, 2000 W/K | 40 | 4.1e-8 J/kg (8.0e-6) | 6.5e-11 K (2.3e-8) | 4.1e-11 J (4.3e-9) | 3.3e-15 J/K (1.7e-9) |
| freezing through a plateau at 500 K, 2000 W/K | 40 | 1.0e-7 J/kg (8.0e-6) | 1.5e-10 K (2.3e-8) | 7.5e-11 J (4.9e-9) | 1.7e-15 J/K (1.9e-9) |

With the wall moving, its velocity agrees to 9.4e-13 m/s (1.2e-10) and the
gas volume to 3.0e-18 m3 (7.2e-16). These are not small steps: in its twelve
the body from 650 K gives 66 J and cools by 132 K, and the bodies on a
plateau leave it, one molten at 351.5 K and one frozen at 464.8 K. With the
kernel's ideal `HeatCapacity` in the body's
place, which returns the arithmetic mean of its end temperatures, the same
comparison is off by 0.065 K and 0.21 K in the body, 0.032 J and 0.10 J in a
step's heat and 0.19 K and 0.61 K in the temperature the contact was given:
more than a million times the bounds, as the test requires.

In those comparisons the conductance is constant (the module's own latent
test makes it depend on the liquid fraction, which the kernel's element does
not). The moving run uses the module fixture's heat bracket of 100 J either
way; the fixed sensible runs use 60 J and the latent runs 20 J, because with
100 J the module refuses the second step of the first and the tenth of the
second (observed: `OutsideBracket` from the gas's own temperature solve at
the far end of the bracket, which its solve tries). The module's answer
depends a little on its bracket: with 60 J the moving run's row read
6.3e-11 K, 4.9e-11 K, 3.1e-11 J and 3.3e-14 J/K. The module's midpoint
thermal step, a conductance that
depends on state, and its real-gas and reactive thermal variants are not
compared.

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
