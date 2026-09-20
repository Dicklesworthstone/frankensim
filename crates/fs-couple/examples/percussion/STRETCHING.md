# Geometric drumhead stretching and physical playing inputs

`drum-stretch`, `drum-stretch-wav` and `drum-stretch-mic` use the same head
geometry, thicknesses, material estimates, installed tensions, eigensolve,
striker, cavity and acoustic boundary as the existing `drum` commands. They
replace only each head's linear potential with the existing fs-plate
`MembraneReduction`: bending and installed tension plus statically relaxed
von Karman stretching. No resonant frequencies are retuned in the renderer.

This matters because stretching creates amplitude-dependent restoring forces
and modal interactions. A harder stroke can change frequency and timbre through
mechanics, not just loudness. Relevant research: Avanzini and Marogna,
“A Modular Physically Based Approach to the Sound Synthesis of Membrane
Percussion Instruments,” IEEE TASLP (2010), DOI 10.1109/TASL.2009.2036903;
Avanzini, Marogna and Bank, “Efficient synthesis of tension modulation in strings
and membranes based on energy estimation,” JASA (2012), DOI 10.1121/1.3651097.
The code does NOT substitute the latter paper's global energy-to-pitch shortcut
for the existing spatially resolved, statically condensed strain model.

## Use

```sh
cargo run --release -p fs-couple --example percussion -- \
  drum-stretch 4096 --strike-speed-m-s 2.0 > drum-stretch.csv

(set -C; cargo run --release -p fs-couple --example percussion -- \
  drum-stretch-mic 48000 20 0.08 0.05 0.35 \
  --strike-speed-m-s 2.0 --strike-position-m 0.06 0.01 > drum-stretch.wav)
```

Named playing options are accepted by the existing drum, prepared-drum, snare
and splash commands too. `--strike-speed-m-s V` sets the actual downward initial
stick-tip velocity, in 0..=20 m/s. The estimated stick effective mass and 0.2 mm
approach gap do not change. Doubling velocity quadruples initial kinetic energy;
it does not multiply the microphone pressure after simulation. Zero launches
no stroke; it does not remove any separately declared initial preload energy.

`--strike-position-m X Y` sets an explicit location in the instrument's reference
XY plane. The existing modal surface fields are interpolated barycentrically;
positions outside the mesh, including a cymbal mounting hole, refuse instead of
snapping to an edge. Contact remains a point idealization with the original
estimated Hertz law and fixed vertical stroke axis, not a resolved finite tip
patch or oblique rolling contact. Omitting both options preserves the old
0.8 m/s stroke and old example contact-location choice.

For an otherwise identical linear comparison replace `drum-stretch-mic` with
`drum-mic`. Do not normalize the two outputs independently. Full-scale pressure
and microphone coordinates retain the meanings documented in AUDIO.md.

## Physical state and numerical limits

The in-plane rim is fixed and the interior is statically relaxed by one cold
fs-la Cholesky factorization shared by every quadratic modal-pair load. The
existing impact/pHS owner advances stick contact, both nonlinear heads and
sealed-air forces together. Membrane data are immutable and shared, not rebuilt
or duplicated each sample. The slope cap is explicitly 0.2. Accepted initial
and endpoint states above it are refused before publishing mechanics or felt
history; internal Newton trials may extrapolate. This cap is a model restriction,
not a discretization error bound or a proof against between-sample overshoot.

Stretching CSV adds batter slope, resonant-head slope and total head stretching
storage to the old columns. That storage is ALREADY part of total energy; do not
add it again. The read-only diagnostics are not synthetic audio channels.

Existing `drum`, `drum-modal` and all snare/splash paths keep their previous
physics unless a playing input is explicitly changed. `drum-stretch` is the
nonlinear reference, not the linear prepared image. The full 160-coordinate
snare-wire bundle is not silently reduced to fit the 64-coordinate nonlinear
host: combining that bundle with stretching awaits a suitable coupled image.

## Verification and limits

Native regressions cover real-film weak/strong periods against independent
continuous Duffing quadrature, two-way air loading, slope refusal/rollback,
playing-input units, unchanged geometric observation maps and contact-driven
nonlinear motion. Native compilation/tests and complete WAV runs were NOT
executed in the implementation environment: cargo/rustc were unavailable.

Executed independent Python checks use a tension-dominated P1 film, NOT the
native DKT pencil. They reconstruct condensed strain energy from nodal motion,
check analytic forces and compare single-mode release with continuous period
quadrature at several time steps. Numerical agreement is not experimental
validation or evidence of realistic sound.

The 80..500 Hz head-mode window and 40..1640 Hz acoustic-fit band remain as in
the original example. Nonlinear harmonics can exceed those bands; neither the
slope cap nor oversampling establishes modal adequacy. In-plane inertia,
wrinkling/slackening, polymer memory, radiation-force feedback, venting and
flexible hardware are absent. BEM still uses the undeformed reference surface.
The general nonlinear time solve still allocates and uses dense finite-difference
Newton work. No native real-time or full-band instrument claim is made.
