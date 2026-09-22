# Distributed enclosed air

`--cavity-modes` enables reciprocal standing-wave air loading on `drum` and
`drum-stretch`, including their `-wav` and `-mic` forms. It works together with
`--prepared-nonlinear`. Without it the original compact-volume model is unchanged.
Prepared `drum-modal` and `snare`/`snare-off` support the same cavity as described
in `SNARE_CAVITY.md`, including explicit momentum drag and vented mechanics.

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
An optional compact sidewall neck supplies additional outward volume flow.

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
With a neck, the uniform pressure includes displaced neck volume; it is not
computed from the old sealed-head contraction alone.

Sealed-cavity exterior audio still uses the original solid-head BEM boundary, propagation,
decimation and PCM conversion. Appended gas coordinates have zero external drive
and zero direct solid-radiation projection. They affect sound by changing head
motion, not by being summed into a manufactured microphone signal. The existing
mechanical and audio clocks are unchanged. Prepared `drum-modal` and `snare`
commands retain their own solver while consuming the same cavity. Cymbal
commands still reject this drum-specific option.

## Vented cavity mechanics

`--cavity-neck radius_m effective_length_m resistance_Pa_s_m3 azimuth_rad z_m`
adds one explicit sidewall opening to `drum` or `drum-stretch` **CSV** mechanics,
with `--cavity-modes`. It also works with `--prepared-nonlinear`, or with the
prepared `drum-modal`, `snare`, and `snare-off` CSV commands. Every physical
parameter is required; there is no hidden vent preset or fitted loss constant.

```sh
# Illustrative declared opening, NOT a measured or calibrated hardware card.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  drum-stretch 4096 --cavity-modes --prepared-nonlinear \
  --cavity-neck 0.005 0.012 1000 0.4 0.08 \
  --strike-speed-m-s 4 --strike-position-m 0.06 0.01 > drum-vent.csv)
```

The opening is a circular patch on the unrolled cylindrical wall: arc offset
`s` maps to azimuth `theta+s/R`. Its area is `S=pi*a^2`. It must fit between the
two heads and have `a<=R/10`. This small-opening chart is not a solid-model cut
of a particular drilled or flared fitting. The example averages all retained
pressure shapes using 8 equal-area radial rings and 32 periodic angular samples
per ring; the reusable `sidewall_averages` API exposes both refinement counts
and a work budget. The average and the pressure-mode norms keep the SAME basis.

The neck has acoustic inertance `L=rho*effective_length/S`. Effective length
includes whatever end corrections the caller declares; the code does not infer
them from an unspecified exterior. In outward displaced volume `v` and flow `Q`,
its equation is `L Q' = average_interior_pressure - resistance*Q`, and `v'=Q`.
Adding the same averaged shapes times `v` to cavity compression makes this a
reciprocal boundary participant, not a one-way filter or an authored oscillator.
The uniform, fixed-wall, lossless limit has `omega_H=c*sqrt(S/(V*effective_length))`.
The actual drum has moving heads and additional pressure modes, so this formula
is not its coupled pitch. The original energy solver includes `L*Q^2/2` and
resistive loss `resistance*Q^2`; there is no second integrator. Initial neck volume
and flow are zero in the example; the reusable `CavityNeck` exposes both explicitly.

Four extra CSV columns give outward displaced volume [m^3], outward flow [m^3/s],
aperture-averaged pressure [Pa], and instantaneous resistive loss power [W]. The
last is an endpoint diagnostic, NOT the step's discrete dissipated energy. The
original `loss_j` and `balance_j` remain the full solver energy accounting.

This model uses a **zero-gauge-pressure reservoir and constant nonnegative
resistance**. It does not supply nonlinear jet separation, flow-dependent losses,
frequency-dependent viscothermal impedance or a distributed duct. Admission
requires `k*max(effective_length,sqrt(S/pi))<=0.3` at each retained acoustic
frequency and fixed-wall neck frequency. That guard is not a certificate for
the nonlinear strike bandwidth, small acoustic pressure, or low neck Mach number;
the caller must assess those limits and refine the cavity/aperture basis.

**Vented WAV/microphone requests refuse before building geometry.** The existing
sealed exterior BEM boundary omits the opening and its radiated flow; reusing it
as though it were the vented instrument would be a false physical model. Sealed
audio is unchanged. Direct vent radiation and its reaction impedance remain
required for a complete vented-instrument audio path.

The generic `CavityCoupling::with_necks` supports up to eight explicit openings,
within the constructor's retained total-coordinate budget, and accepts finite-area
averages from other geometries. Reference execution remains capped at 64;
the prepared snare admits 256 coordinates including all wires, air and necks.
The drum CLI currently exposes one opening.
No new runtime dependency or changed instrument/material/pitch preset is added.

Focused Rust tests cover Helmholtz frequency, physical volume/flow and resistance
work, pressure-basis rescaling, multiple openings, cancellation/exact retry,
finite-aperture Fourier integrals, command admission and struck stretching-head
feedback. These tests are authored; local Rust/Cargo/DSR/RCH were unavailable
during implementation, so native success is not claimed here.

## Reusable API and limits

`impact::cavity::CavityCoupling` accepts the existing generic `CavityModes` and
its structural overlap matrix. It is not limited to a cylinder: another cavity
geometry can supply the same carrier. `cylinder::CylindricalCavity` supplies a
parameterized cylindrical basis and point/interface sampling. Head, cavity and
contact bases must all use the declared ordering and physical units.

The example defaults to zero acoustic momentum drag. `--cavity-drag-per-s D`
explicitly applies a finite nonnegative rate to each nonuniform gas momentum;
the uniform compression mode has no momentum and remains undamped. Both
reference and prepared execution retain this same physical resistance. See
`SNARE_CAVITY.md` for units, bounds, port budgets and sealed-audio examples.
No thermoviscous or wall-loss coefficients are inferred, and the generic
adapter still refuses automatic conversion of a frequency-domain hysteretic
loss factor. This work does not model external
radiation reaction, sidewall elasticity, noncompact vents,
nonlinear gas compression or a complete timpani kettle. No calibration,
full-band convergence, native allocation or real-time deadline claim is made.
