# Pressure-mediated hi-hat motion

`--squeeze-film INPUT.fsf` adds a distributed Reynolds-film load to `hihat`,
`hihat-wav`, and `hihat-mic`. Both actual thickness-offset inner skins supply
the local separation and reciprocal modal force projection. Pressure changes
mechanical motion before the acoustic observer; it is not an output filter,
open/closed crossfade, or synthesized air-noise channel.

This is a **quasistatic incompressible, no-slip, thin-gap viscous image**. It has
no pressure storage, fluid inertia, turbulent edge flow, or rarefied-gas law.
The caller must establish small slopes, negligible inertia/compressibility,
and continuum/no-slip applicability. Numeric admission alone does not prove
those assumptions. No calibrated instrument, native render, listening result,
or real-time performance is claimed.

## Complete input

The strict UTF-8 format permits comments and blank lines, with an 8 KiB ceiling.
Each record occurs exactly once. All quantities are SI:

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
and maximum absolute **gauge** pressure. Keep the pressure limit small relative
to ambient pressure. These are refusal thresholds, never clamps or substitute
gaps.

`boundary` supplies inner then outer conditions: `open` connects the actual
face aperture to zero-gauge ambient pressure; `sealed` supplies zero normal
flow. At least one boundary must initially be open. An omitted bell region is
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

## One pressure law, one mechanical time equation

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

`fs-flux::resistive_film` owns pressure elimination, alongside its existing
compressible gas-film image. `fs-couple` supplies trial configuration and the
same discrete-gradient effort used by contact. `fs-phs` remains the only time
owner. Analytic Newton includes `L dp = A B dv - (dL) p`, not lagged pressure.
Pressure/tangent evaluation uses bounded stack scratch; throughput is unmeasured.

A nonpositive channel aperture blocks that passage exactly. Collapsed volume
cells and disconnected trapped pockets refuse: incompressible resistance
cannot determine trapped gas pressure. There is no gap floor, diagonal repair,
pressure clipping or invented leakage. Endpoint drainage/gap checks precede
publication of motion and felt history, as do the usual energy/cancellation
gates. Compressible storage and contact/topology handoff remain separate work.

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
