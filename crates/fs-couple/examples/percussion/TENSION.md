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
satisfy that equilibrium. This family can resolve directional tuning and a
spatial gradient, but it does **not** reconstruct discrete lug forces, hoop
compliance, bearing-edge friction, wrinkling or a measured specimen.

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
