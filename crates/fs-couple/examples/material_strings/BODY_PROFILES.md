# Directional body materials in the shared string performer

`material_strings PERFORMANCE NEW_OUTPUT.wav --body-profile PROFILE` selects
an explicit rectangular material chart instead of the legacy isotropic body
image. Omitting the flag retains the original construction and material profile.
The performance's `body` record still supplies its pack, dimensions, thickness,
mode count and observer distance. String layers, support prescriptions, torsion,
loads, plucks, bow tracks, shared gas and optional cavity are unchanged.

This is the existing `from_viscous_material_plate_chart` path, not a new solver.
The material tensor is projected on the SAME eigenbasis as bridge displacement,
force and radiation. Full modal viscous cross terms participate in the shared
string/body/friction solve; the body is not a bank of independent decay filters.
No extra scalar body damping or an elastic-modulus-derived viscosity is added.

## Run the synthetic fixture

From the repository root, with the existing toolchain and pinned dependencies:

```bash
demo="$(mktemp -d)"
cp crates/fs-couple/examples/material_strings/*.tsv "$demo/"
cp crates/fs-couple/examples/material_strings/pluck.material-strings "$demo/"
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/wire-manifest.tsv" --out "$demo/wire.fsmatpk"
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/directional-body-manifest.tsv" --out "$demo/body.fsmatpk"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/pluck.material-strings" "$demo/directional.wav" \
  --body-profile crates/fs-couple/examples/material_strings/orthotropic.body
```

The coefficients are synthetic, not measured spruce, maple or any other material.
Change only the body pack or its grain angle to compare results; there is no pitch
normalization, per-mode gain fitting or replacement of the physical string input.
Complete source and model admission still precedes `create_new` output.

## Profile grammar

The fixed-order UTF-8/LF records, with no ignored fields, are:

```text
frankensim-string-body-v1
mesh NX NY
elastic orthotropic12
grain ANGLE_RAD
support simply-supported
pretension FORCE_PER_LENGTH_N_M
window LAMBDA_LOW LAMBDA_HIGH
bridge IX IY
```

`NX`/`NY` are interval counts in `2..=12`. Node `(IX,IY)` is at
`(IX*length/NX,IY*width/NY)` and must be interior. This is a physical nodal point
attachment, not a supplied modal vector; all strings currently use this one
footprint. The separate string records still specify each bridge's stiffness,
damping and station on its complete span including afterlength.

Elastic models are `isotropic`, `orthotropic12` or `orthotropic21`; the ordered
plane and grain angle are passed unchanged to the original material chart.
Support is `simply-supported` or `clamped`. Pretension is nonnegative conservative
plate force per unit length and does NOT acquire material bending loss. The
spectral window contains squared angular frequencies in rad²/s², not note names
or frequencies in Hz. Existing modal/Nyquist/source-band admission still applies.
The finite spatial basis and coupled dashpot time discretization require their
usual mesh, mode-count and time-step refinement.

Every selected body pack supplies `density` and either `young_modulus` plus
`poisson_ratio`, or `young_modulus_1`, `young_modulus_2`, `poisson_ratio_12` and
`shear_modulus_12`. It also supplies ALL six reduced in-plane viscosities in Pa s:
`plane_stress_viscosity_11`, `_22`, `_66`, `_12`, `_16`, `_26`. The engineering
shear convention is the original shared material owner's convention. Signed
couplings and explicit zero entries are preserved; the tensor must be positive
definite or explicitly all zero. Missing entries never mean zero. The legacy
`kelvin_voigt_bending_viscosity` cannot substitute for this complete tensor.
All original state, units, source-band, coefficient and matrix checks apply.

Profile input is bounded to 64 KiB; the original source-pack and complete-network
limits remain. Viscous projection is capped at 16,777,216 operations and 16 body
modes. The dashpot factor uses a 1 N s/m coordinate scale; this changes neither
the projected damping matrix nor material coefficients. The exact profile bytes,
original packs, regional material state, tensor and projected damping remain
attached to the performer. No live source-state mutation is implied.

The current model is a fixed-state, flat, linear rectangular DKT plate with one
uniform material region, not a carved shell, soundpost, bass bar, glue interface
or empirical instrument model. Pressure remains the inherited compact baffled
body/neck observation, not full exterior radiation or radiation back-loading.

Native regression target:

```bash
cargo test -p fs-couple --example material_strings body::tests -- --test-threads=1
```
