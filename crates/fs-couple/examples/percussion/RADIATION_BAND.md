# Explicit acoustic bandwidth and source-preserving boundary refinement

`--radiation-spec preparation.fra` replaces the hard-coded acoustic preparation
window and fitting order with an explicit bounded request. It applies to every
existing drum, snare and cymbal `-wav`/`-mic` command, including stereo. It does
not retune structural modes, add a second resonator, change contact or material
history, or change the mechanical and PCM clocks. CSV commands reject it.

```sh
# Requested 40..6000 Hz observation of the same nonlinear cymbal trajectory.
# This is an offline preparation request, not a claimed completed render.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  splash-mic 4800 20 --analytic-newton --impact-substeps 8 511 \
  --radiation-spec crates/fs-couple/examples/percussion/radiation-6khz.fra \
  --microphone-right -0.08,0.05,0.35 > wider-band-cymbal.wav)
```

The supplied file is complete and editable:

```text
frankensim-radiation-preparation-v1
band_hz,40,6000
training_intervals,64
max_order,24
subdivisions,1
max_panels,4096
max_dense_work,20000000000000
```

All six records are required exactly once. Whitespace and `#` comments are
allowed; unknown, duplicate, missing or extra fields refuse. Input is limited
to 4096 UTF-8 bytes. The first frequency is positive; the last is larger and
strictly below 21600 Hz, leaving the existing 48 kHz observation Nyquist guard.
The unchanged default request is 40..1640 Hz, 20 training intervals, maximum
order eight, no refinement, 2048 panels and one trillion dense-work units.
Defaults do not fill missing fields in an explicitly supplied file.

## Fitting with a genuinely separate audit

For `N` training intervals, the BEM produces a uniform `4*N+1` frequency lattice.
Indices `4k` are training data; `4k+2` select model order; odd indices are a final
independent audit. The example has 65 training, 64 order-selection and 128 audit
samples. Orders 2, 4, ... up to `max_order` use the existing vector-fitting,
conjugation, prewarping and proper stable digital state-space owners. The first
order passing the selection gate is audited once. An audit failure is terminal:
its data cannot select another order, poles, response scale or retry.

`max_order` must be even in 2..32, and `N` is in `2*max_order..128`. Each fit is
normalized only by its training/selection peak magnitude and restored to its
original signed SI response afterward. Both selection and independent audit
retain the original 15% maximum and 5% RMS relative-error thresholds. These are
source-row peak-relative errors, not pointwise relative errors near a response
zero. Exact-zero source rows remain exactly zero; nonzero audit data against
zero training data refuse. Error limits are not tunable input fields.

One boundary-integral formulation is used across the complete transfer: Plain
CBIE only when the entire band satisfies `k*radius < 0.5`, otherwise
Burton-Miller. Switching discrete formulations partway through a transfer can
introduce a numerical discontinuity. Each frequency's geometry factorization
and source solves are shared between receivers; each microphone retains its
own Green response, filter fit, propagation delay and runtime history.

Successful sampled audits are not a bound between samples or outside the band.
The final regression additionally checks fresh off-grid frequencies, but that
is still not a continuous-frequency or measured-instrument certificate.

## Geometry, spatial resolution and work limits

`subdivisions` in 0..4 reuses the existing conforming aperture refiner, marking
every edge. One level makes four children per source triangle. Midpoints stay
on the same polyhedron; no smoothing, sphere projection or reconstruction of
missing curvature is performed. Every child inherits its original panel's
piecewise-constant normal-velocity row, including signed source values. This
preserves the integrated signed flux of every mechanical input and retains all
state addresses, including the separately supplied prescribed-neck input.

This refines the **acoustic discretization**, not the structural mesh or modal
basis. It cannot restore unresolved structural shape variation or add missing
high-frequency modes. A better measured source geometry, physical resolution
and modal convergence remain independent requirements. The original boundary
construction and mechanical mode limits remain; `max_panels` in 1..4096 is an
explicit separate cold-acoustic ceiling below the BEM owner's dense limit.

Before allocating the refined boundary or starting frequency solves, the host
charges `(4*N+1) * (P^3 + 2*P^2*M + P*M*R)` work units for `P` refined panels,
`M` source fields and `R` receivers. `max_dense_work` is positive and at most
100 trillion. This counter bounds the dominant dense algebraic work; it is
not a measured floating-point operation count, memory receipt, elapsed-time
limit or real-time qualification. Order and sample counts bound fitting
separately. Large admitted requests can still be expensive.

The highest requested frequency is solved first, so the existing BEM
six-panels-per-wavelength refusal is encountered before the lower frequencies.
Nonfinite diagnostics and negative radiation power beyond roundoff still
refuse. No denser fit, larger order or wider band overrides an invalid BEM
solution. Neither panel size admission nor fit quality alone certifies acoustic
boundary convergence, especially near two closely spaced cymbal faces.

For prescribed-vent audio the compact-neck `ka` and `kL` limits are checked at
the **requested upper frequency**, not merely at 1640 Hz. A neck valid in the
old band may refuse a wider one. The uniform-flow neck approximation and
prescribed one-way vent interpretation are not silently broadened.

## Mechanics and failure behavior

Both receivers consume the same mechanical trajectory and decimated generalized
accelerations. Source force files, two sticks, nonlinear heads and wires,
carriers, supplied shell/drum geometry, compliant mutes and material history
remain in their existing owners. The actual sample counts, decimator histories,
PCM full scale and propagation meaning are unchanged. Widening the observer
never changes a force, damping coefficient or mode to fit a desired sound.

A preparation refusal or cancellation occurs before mechanical stepping. The
passed cancellation gate is checked between refinement levels, frequency
solves and individual fits; the existing dense BEM/LU itself is not made
interruptible. As before, this offline WAV front door publishes no partial
file on failure. A late mechanical or observer failure may have advanced the
instrument: discard that failed experiment, rather than treating it as a
transactional audio-block retry.

Radiation is still linear and one-way about an undeformed stationary boundary.
There is no radiation mass/damping feedback, moving-boundary solve, room model,
mute scattering or calibrated microphone electronics. A wider acoustic window
is **not full-band cymbal realism**: nonlinear harmonics, structural truncation,
geometry, radiation convergence and measured validation must also be checked.
The estimated 6 kHz request demonstrates syntax and explicit work limits; no
successful complete render, listening result or calibrated specimen is claimed.

Focused native tests:

```sh
cargo test --release -p fs-couple --example percussion observer_fit -- --test-threads=1
cargo test --release -p fs-couple --example percussion radiation_spec -- --test-threads=1
cargo test --release -p fs-couple --example percussion radiation_band -- --test-threads=1
cargo test --release -p fs-couple --example percussion cold_boundary_refinement -- --test-threads=1
```
