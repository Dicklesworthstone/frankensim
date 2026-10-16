# Spatial drumhead tuning through installed tension

A drum specification may add an independent `tension_variation` record for each
head. It changes the real prestress matrix **before** computing modes, contact
projections, cavity coupling and the acoustic boundary. The existing film
material, thickness, mass, fixed rim and time integrators remain unchanged.
No oscillator frequency, audio pitch control or output gain is supplied.

```text
tension_variation,batter,300,-200,80,600,-200,100,400
tension_variation,resonant,-100,150,-40,-300,120,100,-200
```

These are illustrative, unmeasured variations. The full example is
`estimated-varied-tension.fsd`; it includes the required original geometry,
head, mesh and frequency-window records:

```sh
cargo run --release -p fs-couple --example percussion -- \
  drum-modal 4096 \
  --drum-spec crates/fs-couple/examples/percussion/estimated-varied-tension.fsd \
  --strike-position-m 0.06 0.01 > tuned-motion.csv
```

The same specification can be used with `drum-stretch`, `snare`, `snare-off`,
and their existing `-wav`/`-mic` pressure outputs. Numerical-image restrictions,
nonlinear head/wire and material compatibility, force programs, both strikers,
mutes and the existing radiation-feedback admission still apply. The original
head/modal/acoustic bandwidth limits are not extended by this option.

## Resolved multi-lug harmonics

A drum specification may add normalized Airy harmonics such as:

```text
tension_harmonic,batter,10,120,-35
tension_harmonic,resonant,10,60,20
```

The order is an integer 2..12. Orders 6, 8, 10 and 12 are useful bounded images
of common equally spaced lug counts. The final cosine/sine coefficients are N/m
and rotate the stress pattern; they are not individual lug forces or per-mode
detuning controls.

Each term comes from a polynomial Airy stress potential, so its installed
resultant is divergence-free in the film interior. Coefficients are normalized
so their pair magnitude is a conservative tensor-norm bound across the disk.
The sum must stay below the certified tensile margin of the base+affine field;
a request that could create slack/compression refuses rather than being clipped.
Order-12-and-below stress integrals use bounded 8x8 Gauss-Legendre Duffy
quadrature for the polynomial Airy family. With no harmonic rows the original
affine centroid assembly and arithmetic remain unchanged.

```sh
cargo run --release -p fs-couple --example percussion -- \
  drum-stretch 4096 --analytic-newton --impact-substeps 8 511 \
  --drum-spec crates/fs-couple/examples/percussion/illustrative-ten-lug.fsd \
  --strike-position-m 0.06 0.01 > ten-lug-motion.csv
```

The tracked card is illustrative, not a measured torque/force survey. Real lug
measurements can be identified into this basis without inventing oscillator
pitches. Hoop inertia, lug compliance, bearing-edge friction and head seating
remain separate mechanisms.

## SI tensor, not seven independent pitches

After the head name, the first three values are constant additions
`xx,yy,xy` in **newtons per metre**. The next four are `c,d,e,f` in
**newtons per square metre**. With the head record's original base tension `T`
and its centered Cartesian coordinates `x,y` in metres, the full installed
membrane-force resultant is

```text
Nxx = T + xx + c*x + d*y
Nyy = T + yy + e*x + f*y
Nxy =     xy - f*x - c*y
```

`Nxy` is tensor shear, not doubled engineering shear. Both in-plane divergence
equations are zero: this is a continuous equilibrated affine stress field with
rim traction `N*n`. Arbitrarily varying a scalar tension would not generally
satisfy that equilibrium. This affine family resolves directional tuning and a spatial gradient. The optional Airy harmonics resolve circumferential multi-lug stress patterns, but neither image infers individual lug forces, hoop compliance, bearing-edge friction, wrinkling or a measured specimen.

The full tensor must be positive definite at every mesh vertex. Because it is
affine, the vertex tensors also bound the interior of every triangle through
convex combination. Slack, compression, nonfinite values and unresolved positive
margins refuse instead of being clipped. The condition covers the actual
polygonal head, not an unmeshed annulus outside it. Centroid integration is exact
for the affine tensor and the existing constant P1 transverse gradients.

Each head may have at most one complete record; duplicate names, unknown heads,
wrong field counts and nonfinite values refuse. Omission or seven explicit zeros
retains the original uniform image and its assembly arithmetic. Negative
variation coefficients are legal only while the **total** tensor remains tensile.
A specification's uniform base is not counted a second time in the varied pencil.

Changed mode shapes, not just changed eigenvalues, feed every downstream
participant. Geometric head stretching remains the existing incremental
membrane energy about that installed prestress. Runtime does not edit tension,
restart notes, rebuild meshes, drop modes or allocate a second resonator.

Focused native checks:

```sh
cargo test -p fs-plate --lib shell::head::tension
cargo test -p fs-couple --example percussion drum_spec::tension_tests
```

The source includes operator/energy, mode splitting, input admission, actual
contact motion and nonlinear two-hand/cavity retry regressions. Their presence
is not a native passing-test, full-band fidelity or real-time claim.
