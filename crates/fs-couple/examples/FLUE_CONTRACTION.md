# Separated mouth flow in the original flue simulation

The optional `FlueVoice::with_mouth_contraction` adds an explicit quasi-steady
pressure loss INSIDE the original mouth/bore solve. It acts on acoustic mouth
flow, not injected jet flow, output pressure or a postprocessing gain. Source
density and sound speed are the original voice's coherent `GasState`; the
window is `flue_width_m * cut_up_m`, not bore area or jet thickness.

For area ratio Cc, physical mouth area A and signed volume flow Q:

```
dp_separation = rho Q |Q| / (2 (Cc A)^2)
power_separation = dp_separation Q >= 0
source = L (Q - Q_old)/dt + R Q + dp_separation + p_bore
```

The generic `bernoulli_aperture::contraction::ContractionLoss` uses the existing
continuous Bernoulli owner and its rationalized pressure partition. It also
accepts an arbitrary nonnegative linear port response, so fixed-window flow
loss does not need an instrument-specific solver. Both flow signs are admitted.
The same-sample affine wave boundary includes direct reflection and history;
full reflection constrains flow without inventing a pressure-release fallback.

Supply an area ratio in (0,1], a source/condition label and an explicit maximum
contracted-flow Mach in (0,0.3]. No coefficient is selected from an instrument
name. Unity area ratio still produces Bernoulli separation loss: the original
linear comparison is selected by NOT attaching a contraction. The original
linear mouth resistance is unchanged. Do not encode the same physical loss in
both that resistance and the optional closure when identifying a specimen.

`FlueBlockDiag::separation_loss_j` is separate from radiation and backward-Euler
numerical loss. All three enter the mouth balance. Actual candidate flow must
pass its Mach limit before waves are committed. As on the original voice,
failed physical callbacks poison the run because internal source history may
have advanced; this does not claim rollback of every source filter or RNG.
Late/repeated attachment refuses rather than replacing a running model.

## Data-only comparison through a physical receiver

Use new output filenames and the SAME full-scale pressure for comparisons:

```sh
cargo run --release -p fs-couple --example flue_contraction -- linear.wav none 101325 1
cargo run --release -p fs-couple --example flue_contraction -- contracted.wav 0.6 101325 1
cargo run --release -p fs-couple --example flue_contraction -- dense-gas.wav 0.6 202650 1
cargo run --release -p fs-couple --example flue_contraction -- finite-windway.wav 0.6 101325 1 0.015
```

Arguments are destination, area ratio (or `none`), static ambient pressure [Pa],
PCM full-scale [Pa], and an optional finite windway length [m]. Omitting the
length preserves the original ideal jet. The example is a 0.3 m, 6 mm radius tube and a fixed
10 mm by 4 mm mouth window, with a 1 mm exciting jet thickness. It runs the
original tonal-interim jet at 96 kHz with fractional tube travel and the
existing gas/geometry-derived compact baffled termination. Its canonical
`flue_contraction.gesture` ramps to 300 Pa, releases after 60 ms and retains
traveling waves and mouth motion to the end of a 100 ms window. Distributed
viscothermal losses and material-wall motion are not selected in this example;
those existing source constructors can carry the same mouth attachment.

The shared multirate differentiating Rayleigh receiver consumes every solved
outlet flow and emits exterior pressure at 48 kHz, with a declared 2 kHz band.
Its propagation and numerical filter delays remain in the output. No extra
radiation load, pitch correction, waveform shaping or peak normalization is
introduced. Output is created only after the complete source/receiver run and
PCM encoding succeed; an I/O error can still leave a partial file. Stderr
reports gas, loss work, mouth-balance discrepancy, retained network energy and
clipping. Changing ambient pressure reconstructs the source, mouth, propagation,
termination and receiver from the same changed gas, not an independent density
knob applied to an existing waveform.

## Finite rectangular windway

`FlueVoice::with_viscous_windway(length_m)` uses the existing exit width and
height, plus the actual inlet `GasState`, in the shared `ViscousSlit` law:

```
p_blow = 12 mu L v / h^2 + rho v^2 / 2
Q_jet = w h v
tau = cut_up / (convection_ratio v)
```

The solved bulk speed drives injected volume, Reynolds number, jet transit and
amplification band before the acoustic feedback runs. A viscosity-only material
data change can therefore alter microphone pressure even at unchanged density
and sound speed. Both schedule compilation and applied controls admit the
actual viscous speed against the declared low-Mach/history limits; they never
validate an ideal speed and then silently clamp a slower jet's delay. Attachment
also rechecks an already held pressure and refuses after the first sample.

This adds a quasi-steady parallel-plate windway resistance. The quadratic term
is exit kinetic head; it is not added again to the acoustic mouth loss ledger.
The existing jet-card profile remains fixed. Entrance development, profile
reshaping, windway inertia, turbulent flow and acoustic modulation of the mean
jet require richer models. Supply a physical rectangular channel length for a
recorder/organ windway or another justified channel; the example's optional
15 mm length is authored geometry, not a measured transverse-flute lip length.
The significance of windway profile and length is established experimentally
by Segoufin et al., *Acta Acustica* 86 (2000), pp.649–661
([author record](https://research.tue.nl/en/publications/experimental-study-of-the-influence-of-the-mouth-geometry-on-soun/));
the present reduced law does not reproduce that measured instrument model.

## Scope and evidence

This is a symmetric fixed-window, incompressible quasi-steady separation
closure with no vortex memory or pressure recovery. The vena-contracta pressure
relation is also given in Appendix B of Sierra-Ausin et al., *J. Fluid Mech.*
943 A48 (2022), DOI 10.1017/jfm.2022.417. That relation does not identify the
coefficient or validate this flute's transient acoustics. Source labels in the
example explicitly describe authored inputs, NOT measurements of flute edges.
The existing jet card, receptivity, small-signal bore and infinite-baffle
observation limitations remain. No resolved turbulent noise, empirical flute
timbre, broadband/real-time qualification or piano-level fidelity is claimed.

Focused native checks:

```sh
cargo test -p fs-couple --test contraction_loss --test flue_contraction --test flue_windway
cargo test -p fs-couple --lib flue_loop::mouth_loss::
cargo check -p fs-couple --example flue_contraction
```

Independent reference equations checked 567 manufactured affine boundaries,
180 independently bracketed nonlinear roots and 18,432 complete delay-loop
steps. These are not execution of the Rust implementation. Local Cargo/rustc
were unavailable during implementation; native workflow results must be checked
separately rather than inferred from numerical reference agreement.
