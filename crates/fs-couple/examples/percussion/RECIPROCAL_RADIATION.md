# Reciprocal exterior acoustics for ordinary percussion playback

The existing `--radiation-spec` file accepts an optional `image,reciprocal`
record. Omission or `image,collocation` preserves the original operator.
Selection happens before acoustic preparation, never after a power/fit failure.

```sh
# Same nonlinear splash and stand/mallet mechanics; alternate exterior image.
(set -C; cargo run --release -p fs-couple --example percussion -- --shell-profile crates/fs-couple/examples/percussion/estimated-splash.profile \
  splash-mic 4800 20 --analytic-newton \
  --radiation-spec crates/fs-couple/examples/percussion/reciprocal-shell.radiation \
  > reciprocal-splash.wav)

# Attach the exterior reaction as well as observing the same pressure field.
(set -C; cargo run --release -p fs-couple --example percussion -- --shell-profile crates/fs-couple/examples/percussion/estimated-splash.profile \
  splash-mic 4800 20 --analytic-newton --radiation-feedback \
  --radiation-spec crates/fs-couple/examples/percussion/reciprocal-shell.radiation \
  > reciprocal-loaded-splash.wav)
```

These commands request renders, not observed successful recordings. Native
execution and deadline/listening qualification remain required. The narrow
reference bandwidth and estimated specimen materials are unchanged.

The existing `fs_bem::reciprocal_radiation` owner solves the actual closed
triangles in its continuous-pressure P1 / auxiliary-density P0 representation.
The microphone evaluates that representation, not a collocation reconstruction
from mean panel pressure. Near-field/directional receivers, finite-point pressure,
far-field phase, shared stereo source solves, independent receiver fitting,
decimation and PCM all retain their existing owners and physical units.

Optional feedback projects the SAME boundary pressures against the transpose
of the SAME velocity rows, retaining signed cross-source work. Its complete
Hermitian power check and passive impedance fit are unchanged. Without feedback
there is no load, impedance fit or acoustic-memory attachment; the observer
cannot alter mechanical motion. Coupled vents retain their existing required
reciprocal image even when an old file omits the image record.

This addresses the missing playback route to an existing variational alternative
for the documented negative-power collocation failure (bead
`frankensim-rc-root-q61wp.68.7`). It does NOT disable the collocation power gate,
clip eigenvalues of measured/computed impedance, replace a refused source with
noise or a modal sum, or claim that the native smoke is fixed without execution.

The owner still refuses bands outside its slab-based interior-resonance
exclusion or angular-resolution chart. A large/deep bowl may not admit the same
band as a thin cymbal; it is not a general Burton-Miller replacement. Full
Galerkin geometry, matrix and frequency-work limits remain enforced in addition
to the existing preparation preflight. The supplied larger work allowance is
explicit and does not make preparation real time. Pressure-space, surface,
angular and spectral convergence remain independent requirements.

Four focused authored Rust regressions check backward-compatible selection,
ordinary one-way thin-body pressure with two receivers and off-grid frequencies,
signed reciprocal work, and budget/cancellation/resonance refusal without
fallback. Run `cargo test --release -p fs-couple --example percussion
reciprocal_playback_tests -- --test-threads=1`. They have not run locally.
