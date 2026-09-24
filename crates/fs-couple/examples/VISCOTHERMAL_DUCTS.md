# Distributed gas-wall loss in played valve ducts

The supplied valve path can now place **viscous and thermal memory inside its
reciprocal acoustic network**, rather than filtering the final waveform. Loss
changes returning pressure waves, valve motion and the flow driving the exterior
receiver. The plate, contact, solid-material memory and pressure phrase still use
the existing owners.

```bash
cargo run --release -p fs-couple --bin music_render -- \
  wind crates/fs-couple/examples/plate-valve-viscothermal.performance \
  /tmp/viscothermal-valve.wav --decimate --block 37
```

The illustrative source retains the preceding 25 x 10 mm plate, asymmetric lay,
Maxwell material history, 250 mm by 7 mm-radius tube and pressure phrase. It runs
at 96 kHz and observes the radiation-loaded outlet 0.2 m away. The explicit loss
use-band is 100–1000 Hz, with 12 spatial cells and eight arms per loss. The outlet
load and receiver upper bands are explicitly 1000 Hz. These are authored source
values and finite-band approximations, **not measured cane or instrument data**.

## Input selection

Append the following optional suffix to a `tube` or `duct_section` record:

```text
viscothermal MIN_FREQUENCY_HZ MAX_FREQUENCY_HZ CELLS ARMS
```

For example:

```text
tube 0.25 0.007 baffled-low-ka 1000 0.002 1048576 viscothermal 100 1000 12 8
duct_section 0 1 0.125 0.007 0.002 viscothermal 100 1000 6 8
```

The first line is a complete tube declaration; the second belongs to an explicit
network. ARMS must be exactly 8, matching the shared owner; other values refuse.
Select each section independently. No suffix means the original
lossless propagation, not an automatically inferred gas-loss model. An unselected
narrow side branch remains lossless; selecting an invalid wide-tube approximation
there refuses instead of silently substituting a different law. The numeric
reflection and geometry-derived radiation alternatives both accept the suffix.
A receiver's upper use-band cannot exceed any selected section's upper band.

`ambient` supplies the actual gas density, sound speed, viscosity, heat-capacity
ratio and Prandtl number through the existing gas model. There is no independently
adjustable attenuation, thermal diffusivity, acoustic radius or output gain.

## Physical and numerical construction

For section radius a, area S and kinematic viscosity nu, add only the first-order
wide-tube Zwikker–Kosten excess terms, with the exp(-i omega t) convention:

```text
Z'_excess = (rho/S) sqrt(2 nu omega)/a (1-i)
Y'_excess = S/(rho c²) (gamma-1) sqrt(2 nu omega/Pr)/a (1-i)
```

The original characteristic delays already contain ideal-gas inertia rho/S and
adiabatic compliance S/(rho c²); neither is added again. Positive Foster arms
from the existing `fs_phs` producer approximate the excess terms. Viscous series
RL states and thermal parallel RC states use the existing `fs_vfit` midpoint
impedance/admittance owners. Every state contributes actual storage and positive
loss to the network's original work balance. Thermal arms are not mislabelled RL
branches, an extra pressure source or a delayed correction to a completed step.

`with_viscothermal_sections` uses symmetric cells: quarter transit, half
viscous load, quarter transit, full thermal load, quarter transit, half viscous
load, quarter transit. Cumulative quarter positions are rounded to the original
mechanical transit grid. All intervening section delays are positive. Their sum is exactly the
original admitted integer transit; original junction and outlet indices remain
unchanged. Geometric length rounding still obeys the original section allowance.
Increasing spatial resolution also requires enough mechanical transit samples.
The returned section records retain original geometry, cell choices and all
actual generated loss-node addresses.

Admission requires 1..128 cells per selected section, at most 1024 cells in
total, and the fixed eight-arm spectrum. Viscous **and thermal** shear numbers
must be at least 10 at the lower band, dt*fmax <= 0.05 and k0*cell_length <= 0.5.
There must be at least four original transit samples per cell. The existing
network payload limit includes the additional states and lowered geometry.

The Foster producer uses fixed poles spanning omega_min/64 to 64*omega_max.
At 65 logarithmically spaced use-band frequencies, both its analogue response
and its implemented bilinear response must pass fixed 5% relative checks on the
complex excess term **and separately its real dissipative part**. A positive
fallback is not accepted merely because it is passive. Failed checks refuse;
coefficients, tolerances and geometry are not retuned to make a source pass.
The actual integer-cell and bilinear-load scattering response is additionally
compared with the homogeneous first-order telegraph section at 65 frequencies;
the fixed maximum absolute scattering discrepancy is 0.03. These are sampled
checks, not a continuum or whole-network error bound. Temporal and spatial
refinement remain necessary for quantitative use.

## Playback and scope

Both `wind` and `ensemble --valve` retain all gas, solid-material, contact,
propagation and receiver state through the existing cancellation/resume path.
Physical travel delay is not cancelled by output-filter latency alignment. PCM
scales remain explicit, and source scales never become hidden ensemble gains.
The sidecar includes original section geometry, bands, loss-node and propagation
section ranges, unchanged transit, and constitutive/scattering discrepancies; expanded nodes retain the actual RL/RC
coefficients. No new waveform writer or test framework is introduced.

This is a **linear, first-order thin-boundary-layer model over a positive declared
frequency band**. It is not DC/Poiseuille flow, the all-regime Bessel law, mean-flow
turbulence, rarefaction, a moving thermal boundary or evolving gas temperature.
The pressure phrase and nonlinear valve are not band-limited automatically; a
finite attack, mean flow or contact transient contains content outside the stated
band. Do not infer full-transient accuracy from the sampled checks or passive
energy accounting. The existing plate, lay, compact radiation and gas-model
qualifications remain in force. Omitted distributed loss and an invalid selected
loss are distinct choices; neither is silently repaired.

The focused native targets remain `plate_aperture` and `music_render_wind`.
The existing core tests cover constitutive response, actual network work and
source motion. The six file/playback additions compare explicit construction,
retained histories, supplied grammar and end-to-end PCM output.
Their execution status must be obtained from the corresponding native run.

This adapter uses the shared `bernoulli_aperture::viscothermal` implementation
introduced by main commit `d77085ad074c6825842020a382298e53b627baa1`; it does
not introduce another loss builder or relaxation-state owner.
