# fs-material — CONTRACT

> **Layout (restored 2026-11-09).** This contract is an index plus part files. Once this file grew
> past 1 MiB (1,048,576 bytes), edits through one write path replaced the whole contract with their
> one new section; no edit below that size ever did. The likely cause is the GitHub contents API,
> which returns empty content for files above 1 MB, so that writer saw an empty file.
> This file was cut to a single section 8 times between 2026-11-05 and 2026-11-09; the last cut (2026-11-09 12:27) removed 1.7 MB.
>
> - All earlier contract text is in `contracts/part-NN.md`, unchanged. A higher part number is newer;
>   inside a part, sections are newest first, as they were in the single file.
> - Add a new section under "New sections" below, newest first.
> - When this file passes 600 KB, move the sections under "New sections" into a new
>   highest-numbered part and add its row to the table. No file here may exceed 900 KB.
> - The ten required sections are in the parts named under their headings below.

## New sections (newest first)

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
