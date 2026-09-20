# Distributed enclosed air

`--cavity-modes` enables reciprocal standing-wave air loading on `drum` and
`drum-stretch`, including their `-wav` and `-mic` forms. It works together with
`--prepared-nonlinear`. Without it the original compact-volume model is unchanged.

```sh
(set -C; cargo run --release -p fs-couple --example percussion -- \
  drum-stretch 4096 --cavity-modes --prepared-nonlinear \
  --strike-speed-m-s 4 --strike-position-m 0.06 0.01 > drum-cavity.csv)
(set -C; cargo run --release -p fs-couple --example percussion -- \
  drum-stretch-mic 48000 20 0.08 0.05 0.35 --cavity-modes \
  --prepared-nonlinear --strike-speed-m-s 4 > drum-cavity.wav)
```

The current example still uses the existing reconstructed 14 x 6.5 inch drum:
clear radius 0.1703 m, depth 0.1651 m, the original estimated PET heads/tensions,
Hertz stick contact, and the original declared gas density/sound speed. It does
not introduce a new drum preset or retune its heads. The cylindrical sidewall is
rigid; the two heads supply boundary motion and receive pressure feedback.

## What is computed

The cylindrical basis separates radial P1 Galerkin functions, exact azimuthal
Fourier pairs, and axial cosines. The radial weak Neumann pencil is solved by
`fs-modal`, not by a new eigensolver or an authored resonance list. The axis is
regular. Volume norms and surface overlaps belong to the same pressure basis.
The example uses 32 radial intervals, angular orders through 4, axial orders
through 2, and a 1,100 Hz window with a maximum of eight complete modes. Both
members of every retained angular pair are kept or the construction refuses.

Independent Python/SciPy evaluation of that discretization retains the constant
pressure mode, two transverse pairs near 590.24 and 979.14 Hz, and the first axial
mode near 1,038.76 Hz. These are **uncoupled cavity basis frequencies**, not the
coupled drum's pitches and not measured production-drum resonances. The Rust
implementation must still be executed and its radial/interface/mode windows
refined. At 8/16/32 radial intervals the first-positive-root relative errors
against analytic Neumann Bessel roots decrease approximately fourfold; at 32
intervals they are 0.0387% (m=0), 0.00724% (m=1), and 0.00962% (m=2). This is an
independent equation check, not a native Rust test or full model-convergence proof.

Each actual head triangle contributes three area quadrature points. The existing
`assemble_coupling` integrates head mode shapes against cavity pressure shapes.
Outward signs are opposite on the two heads. The zero mode replaces the old gas
spring exactly at the formulation level; it is **not added on top** of it.
Nonzero modes add acoustic inertia to the same nonlinear energy/storage solve,
so stored acoustic energy changes subsequent head motion and stick interaction.

## Observations and sound

`cavity_internal_pa` remains the uniform compression pressure component. Two
additional CSV columns report pressure at interior points (0.4R,0.2R,0) and
(-0.4R,-0.2R,L), exposing spatially nonuniform pressure. They are not outside
microphones and do not inject forces.

Exterior audio still uses the original solid-head BEM boundary, propagation,
decimation and PCM conversion. Appended gas coordinates have zero external drive
and zero direct solid-radiation projection. They affect sound by changing head
motion, not by being summed into a manufactured microphone signal. The existing
mechanical and audio clocks are unchanged. `drum-modal`, `snare` and cymbal
commands reject the option rather than silently changing their physical model.

## Reusable API and limits

`impact::cavity::CavityCoupling` accepts the existing generic `CavityModes` and
its structural overlap matrix. It is not limited to a cylinder: another cavity
geometry can supply the same carrier. `cylinder::CylindricalCavity` supplies a
parameterized cylindrical basis and point/interface sampling. Head, cavity and
contact bases must all use the declared ordering and physical units.

The new example declares zero acoustic momentum drag; no fabricated
thermoviscous or wall-loss coefficients are fitted. The generic adapter accepts
explicit passive momentum drag but refuses automatic conversion of a nonzero
frequency-domain hysteretic loss factor. This work does not model external
radiation reaction, sidewall elasticity, vents, snare/cavity co-simulation,
nonlinear gas compression or a complete timpani kettle. No calibration,
full-band convergence, native allocation or real-time deadline claim is made.
