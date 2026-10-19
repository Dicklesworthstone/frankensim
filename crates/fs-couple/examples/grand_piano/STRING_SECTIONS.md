# Material-driven bonded speaking strings

The grand-piano example accepts **`--string-sections FILE`**. Both input versions
bind source-declared layer properties through `fs-matdb` claims and the shared
`fs-material::section` owner. The same section supplies speaking-string mass
per length, bending rigidity, and axial rigidity. Version 2 also supplies
internal bending loss. The existing contact, string, bridge, soundboard,
and acoustic output owners consume these coefficients without selecting
physics by a material's name.

## Elastic section (v1)

Input is UTF-8 with these two header records:

```text
frankensim-string-sections-v1
layer,key,outer_radius_m,density_kg_m3,young_pa,source,license
layer,69,0.0005,7800,200000000000,illustrative synthetic material,CC0-1.0
```

This file covers a one-key scale containing MIDI key 69. For a full piano,
declare at least one ordered layer for **every** key, including unplayed
keys. Another row for the same key with a larger radius adds an annular
layer. The first layer starts at radius zero. Hollow strings, missing keys,
invalid properties, reversed radii, and unknown keys refuse atomically.
There may be at most 128 layers per key. Source and license must be
nonempty; importing a declaration does not independently verify it.

The shared section integrates `mu = sum(rho A)`, `EI = sum(E I)`, and
`EA = sum(E A)` over these exact annuli. Tension remains the imposed
prestress. No automatic retuning hides a density or stiffness substitution.

## Elastic section and internal bending loss (v2)

Use the same CLI option with this versioned input:

```text
frankensim-string-sections-v2
layer,key,outer_radius_m,density_kg_m3,young_pa,bending_viscosity_pa_s,omega_min_rad_s,omega_max_rad_s,source,license
layer,69,0.0005,7800,200000000000,100000000,0,150000,illustrative synthetic loss law,CC0-1.0
```

These numbers are a **synthetic example, not measured metal loss data**.
The viscosity must describe the solid's bending strain-rate response
in Pa s. Fluid viscosity with the same units is not interchangeable.
Each layer supplies a finite nonnegative viscosity and an ordered finite
angular-frequency interval in **rad/s**. Zero viscosity is an explicit
lossless internal bending law; an omitted value is not treated as zero.

`ConcentricElasticSection::with_kelvin_voigt_bending` integrates
`C_b = sum(eta I)` and intersects every layer's source-supported band.
The same frozen binding carries mu, EI, EA, C_b and that interval into
the existing string bank. Each mode receives the Kelvin-Voigt rate
`C_b k^4 / mu` multiplying its velocity relative to the moving bridge.
The bridge receives the reciprocal reaction through the existing passive
damping operator. An uncurved moving chord gains no bending-loss traction.

This applies to every physical unison member, both transverse directions
when enabled, and each speaking and duplex segment. Preparation checks
each retained mode's **actual** frequency after tension offsets and any
requested tuning. A short duplex or detuned member outside the supported
band refuses preparation; the renderer does not extrapolate the loss law.

V2 replaces the demonstration Maxwell bending spectrum. The separately
authored common velocity-drag rate remains an approximation, and the render
report states that scope. V1 continues to use the explicitly authored
demonstration bending spectrum because density and modulus do not identify
internal loss. Neither version can use the source-specific
`--rt0425-string-damping` law, which is not identified for a substituted
material. A later v1 substitution clears an earlier v2 loss binding.

## Rendering and nonlinear extension

For a one-key scale using the ordinary SI scale schema:

```bash
cargo run --release -p fs-couple --example grand_piano -- \
  --scale one-key-scale.csv --string-sections a4-bonded.fsss \
  --render output.wav --note 69 --duration 3
```

Change the material file while retaining geometry and tension. Density
changes inertia and frequencies, modulus changes stiffness and
inharmonicity, and v2 viscosity changes internal decay. These changes
reach the same time-domain mechanics used by the pressure WAV; no output
filter or timbre control is inserted. Explicit `--concert-pitch` or
`--stretch-tuning` retunes physical tension after material binding.

Both versions work with `--string-stretching`. Use the token `section`
in place of numeric EA in each stretching or longitudinal-wave row, so
geometric extension, nonlinear tension and axial wave speed use the
same layer integrals. Strain, slope and retention limits remain explicit.
See [STRING_STRETCHING.md](STRING_STRETCHING.md) for those rows.
`--dump-scale` exports the effective elastic coefficients for inspection;
it is not a loss-law or provenance round trip. Preserve the section input
as the material declaration.

## Physical scope

This is a frozen, perfectly bonded concentric rod with a linear elastic
law and optional Kelvin-Voigt bending loss. It does not model axial creep,
temperature evolution, nonlinear damping, fatigue, contact properties, or
winding slip. A loose helical bass-string winding needs its own pitch,
wrapping inertia and slip mechanics; it must not be represented as a bonded
annulus merely because the metal is copper, silver, or gold.

The source records remain declarations, not certified specimen measurements.
The existing energy, mesh and audio-band limits still apply. The added
regressions exercise viscosity-only trajectory changes and energy closure,
moving-bridge reactions, layered integrals across all span/direction groups,
actual modal support checks, and missing/zero/invalid data.
