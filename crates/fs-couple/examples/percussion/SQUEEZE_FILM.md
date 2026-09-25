# Pressure-mediated hi-hat motion

`--squeeze-film INPUT.fsf` adds a distributed Reynolds-film load to `hihat`,
`hihat-wav`, and `hihat-mic`. Both actual thickness-offset inner skins supply
the local separation and reciprocal modal force projection. Pressure changes
mechanical motion before the acoustic observer; it is not an output filter,
open/closed crossfade, or synthesized air-noise channel.

The default is a **quasistatic incompressible, no-slip, thin-gap viscous image**.
An explicit `isothermal` record selects the compressible storage image described
below, using the same geometric cells and passages. The default image has
no pressure storage, fluid inertia, turbulent edge flow, or rarefied-gas law.
The caller must establish small slopes, negligible inertia/compressibility,
and continuum/no-slip applicability. Numeric admission alone does not prove
those assumptions. No calibrated instrument, native render, listening result,
or real-time performance is claimed.

## Complete input

The strict UTF-8 format permits comments and blank lines, with an 8 KiB ceiling.
Each required record occurs exactly once; `isothermal` is optional and occurs
at most once. All quantities are SI:

```text
frankensim-squeeze-film-v1
annulus,0.02,0.1,16,4
viscosity_pa_s,0.000018
limits,0.000001,0.002,1000
boundary,open,open
```

These are **synthetic test inputs**, not measured specimen data or a guarantee
that the included hi-hat fits this gap domain. Choose compatible geometry,
physical boundaries and validity limits; the program will not widen a limit
or relocate a sample to force a successful run.

`annulus` supplies inner/outer radii and radial/azimuthal cell counts. Inner
radius must be positive, with at least one radial cell, three azimuths and at
most 64 cells total. This work ceiling is not convergence certification.

`limits` supplies minimum retained volume-cell gap, maximum cell/channel gap,
and maximum absolute **gauge** pressure. In the default incompressible image,
keep the pressure limit small relative to ambient pressure. These are refusal thresholds, never clamps or substitute
gaps.

`boundary` supplies inner then outer conditions: `open` connects the actual
face aperture to zero-gauge ambient pressure; `sealed` supplies zero normal
flow. Without compressible storage, at least one boundary must initially be open. An omitted bell region is
not automatically an ambient vent, and mounting hardware may obstruct a hole.
No vent is inferred from a material name or missing cells.

```sh
cargo run --release -p fs-couple --example percussion -- \
  hihat path/to/pair.fshh 12000 \
  --squeeze-film path/to/film.fsf \
  --analytic-newton --impact-substeps 8 511 --strike-speed-m-s 0
```

The complete pedal program must fit the duration. Existing two-stick, force
playback, stereo microphone and radiation options remain available. Without
this flag the original construction is unchanged. CSV `loss_j` includes fluid
dissipation in the existing mechanical ledger.

## Default resistance image: one mechanical time equation

For cell closure rows B, projected sector areas A and the conductance graph L:

```text
h = h_reference - B q
G = width * max(h_face, 0)^3 / (12 * viscosity * length)
L(h) p = A B v
D = transpose(B) A p
v dot D = sum_faces G * pressure_drop^2 >= 0
```

The same areas and closure rows supply both swept-volume forcing and pressure
reaction. Radial and periodic azimuthal passages sample their own physical
apertures. The lower shell uses the same proper rotation as contact. Missing
or ambiguous surface locations refuse rather than snap. Sticks and pedal have
zero direct fluid-force columns; pressure reaches them through contacts/mounts.
Common rigid translation generates no squeeze load.

`fs-tribo::resistive_film` owns the shared graph and pressure law;
`fs-flux::resistive_film` remains a compatibility reexport. `fs-couple` supplies trial configuration and the
same discrete-gradient effort used by contact. `fs-phs` remains the only time
owner. Analytic Newton includes `L dp = A B dv - (dL) p`, not lagged pressure.
Pressure/tangent evaluation uses bounded stack scratch; throughput is unmeasured.

A nonpositive channel aperture blocks that passage exactly. Collapsed volume
cells and disconnected trapped pockets refuse: incompressible resistance
cannot determine trapped gas pressure. There is no gap floor, diagonal repair,
pressure clipping or invented leakage. Endpoint drainage/gap checks precede
publication of motion and felt history, as do the usual energy/cancellation
gates. Trapped gas requires the compressible image below; moving contact and
topology handoff remain outside these fixed-chart images.

The projected grid/material chart stays fixed. Sliding seals, changing normals,
roughness leakage and fluid acoustic radiation are outside this image. BEM
still observes a stationary reference scene; air loading does not make its
scattering geometry follow the closing gap.

## Focused checks

```sh
cargo test --release -p fs-flux --no-default-features --lib resistive_film
cargo test --release -p fs-couple --lib render::plate::impact::squeeze
cargo test --release -p fs-couple --example percussion hihat -- --test-threads=1
```

Authored tests cover exact single-cell pressure, inverse-cubic resistance,
suction, reciprocal work, analytic tangent, closure/trapping, transactional
refusals, joint momentum/energy, prepared/reference agreement and annular
refinement. Native execution remains required. An independently executed
numerical reference gives annular load errors of 53.6%, 14.2%, 3.63% and 0.914%
with 2, 4, 8 and 16 radial cells at four azimuths in a uniform synthetic gap.
That illustrates refinement, not instrument fidelity or a Rust test result.

## Compressible storage: seal, retain pressure, then reopen

Add one optional record to the same input file:

```text
isothermal,100000,293.15,287.05
```

The entries are absolute ambient pressure [Pa], fixed gas temperature [K], and
specific gas constant [J/(kg K)]. These illustrative air-like inputs are not
measured instrument conditions. All must be positive and finite. Omission
preserves the original resistance-only model; missing data never selects gas
properties by name. Duplicate records refuse. The existing `limits` still
bounds volume-cell gaps and absolute **gauge** pressure; it never clamps them.
The external acoustic medium is not silently retuned by this record.

Both `boundary,sealed,sealed` and open boundaries are legal in this image.
A sealed boundary means an explicitly supplied zero-flux physical boundary,
not an automatic seal at an arbitrarily cropped annulus. Every retained cell
begins at ambient pressure for its actual initial volume. As the skins move,
cell masses evolve with the same fs-phs step as shell/stick/pedal motion and
felt memory. A closed passage has exactly zero mobility. It traps the gas
already present, rather than losing its mass, resetting pressure, adding a
leak, or requiring an incompressible pressure solve. Reopening resumes flow
from those retained masses and the current gaps. Volume-cell collapse and
nonpositive mass/pressure still refuse; gas cannot be squeezed to zero volume.

For cell volume V and actual gas mass m, absolute pressure is `p=m R T/V`.
The storage used with gauge-pressure wall reaction is the relative free energy

```text
H_gas = p0 V [u log(u) - u + 1],  u = p/p0.
```

Its volume derivative supplies the same reciprocal force on the two skins.
At the physical gradient, an open passage transports signed mass at
`G (p_i^2-p_j^2)/(2 R T)`. Positive logarithmic-mean mobility extends that law to
the time owner's discrete-gradient effort while preserving nonnegative
free-energy dissipation and pairwise internal mass conservation. Gas gradient,
Hessian and mobility tangents are analytic, including pressure and gap changes.
There is no second time integrator or endpoint force correction.

CSV appends `gas_min_absolute_pa`, `gas_max_absolute_pa`, `gas_mass_kg`, and
`gas_free_energy_j` only for this selection. These are endpoint physical gas
observations, **not microphone channels**. Gas free energy is already included
in `total_energy_j`, and its dissipation in `loss_j`; do not add it twice.
The accounting is relative to a fixed pressure/temperature reservoir, not the
total internal thermal energy of an isolated gas. No gas heat capacity, thermal
relaxation, shock, turbulent jet or whistle sound has been introduced.

`hihat-wav` and `hihat-mic` use the resulting physical shell motion with the
unchanged acoustic path. Source indices, two flexible sticks, player-force
clock, prepared analytic execution, material/radiation memory and full-state
rollback remain shared. Gas states are appended after the existing histories;
none becomes a structural or radiation source coordinate. Fixed-reference BEM
still does not follow gap-dependent scattering geometry.

The model assumes ideal gas, constant temperature, low-inertia laminar thin-film
flow and no slip. These are caller-validated applicability conditions, not
facts certified by passing the numeric gap/pressure checks. The standard
isothermal compressible Reynolds relation is described in
[COMSOL's modified Reynolds gas-flow reference](https://doc.comsol.com/6.3/doc/com.comsol.help.cfd/cfd_ug_fluidflow_thinfilm.12.21.html).
The discrete free-energy/mobility construction above is explicit in the code.

Focused native tests (authored; execution results must be checked separately):

```sh
cargo test --release -p fs-tribo --lib resistive_film::isothermal
cargo test --release -p fs-couple --lib render::plate::impact::gas_film
cargo test --release -p fs-couple --example percussion hihat::gas_tests
```
