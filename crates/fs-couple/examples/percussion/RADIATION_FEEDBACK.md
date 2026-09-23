# Radiation reaction inside the played drum/cymbal mechanics

`--radiation-feedback` adds a passive approximation of the actual exterior
BEM force/velocity matrix to the SAME nonlinear time solve as sticks, heads,
cymbal shells, enclosed air, wires and felt. The microphone then observes the
motion WITH that reaction. This is not additional PCM damping or a one-way
pressure render labelled as a load.

```sh
# Nonlinear drum with sealed cavity and actual exterior reaction.
# These are preparation requests, not claims of completed/calibrated renders.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  drum-stretch-mic 4800 20 --cavity-modes \
  --radiation-feedback --analytic-newton --impact-substeps 8 511 \
  --microphone-right -0.08,0.05,0.35 > loaded-drum.wav)

# One pair of independently vibrating/contacting shells, one loaded BEM scene.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  hihat-mic crates/fs-couple/examples/percussion/estimated-hihat.fshh 4800 20 \
  --radiation-feedback --analytic-newton --impact-substeps 8 511 \
  > loaded-pair.wav)
```

The option accepts the nonlinear-capable `splash`, `drum`, `drum-stretch`,
`snare`, `snare-off`, and `hihat` pressure command families (`-wav`/`-mic`).
An ordinary linear-only `drum-modal` or snare image cannot retain this memory
and refuses. Snare playback must explicitly select its nonlinear-capable image,
for example `--head-stretching`. No strands are removed to make the load fit.
CSV commands and vented drums refuse: the latter need separation of existing
neck end corrections from exterior loading to avoid double-counted inertance.

## One geometry solve for reaction and microphones

For every frequency, the existing BEM batch already solves the complete normal
velocity basis. Those SAME surface pressure fields now supply all receivers
AND `Z_ij = integral(b_i * p_j) dA`, after conversion from the existing unit
acceleration input to unit velocity. Both projections use actual panel areas
and unchanged signed modal rows. Reaction on the solid is `-Z*v`. No modal
averaging, output gain, second mesh, or per-receiver mechanical solve occurs.

The complete Hermitian power matrix is admitted with the existing fs-phs PSD
check; positive diagonal powers alone are insufficient. Raw complex BEM data
are never clipped or made symmetric to force admission. Failed numerical
passivity, spatial resolution or model fit refuses before any mechanical step.

The positive-real pole model and fixed-dictionary PSD-residue fitter are shared
with `grand_piano` rather than copied. The piano's split runtime is unchanged;
percussion instead uses its existing joint Gonzalez discrete-gradient solver.
Each acoustic pole has positive frequency, nonnegative damping, and a signed
coupling row. The same row drives its stored energy and transposes the reaction
onto the original mechanical momenta. This guarantees passive model structure,
not accuracy of the model for an arbitrary supplied boundary.

Even-index BEM samples fit coefficients. Odd-index samples are held out.
Training-only resistive weighting prevents small radiation resistance from
being hidden by a larger reactive load. Existing complex peak/RMS error limits
are 5%/2%; Hermitian-power peak/RMS limits are 10%/5%. They are not user-adjustable.
The observer retains its separate order-selection/final-audit scheme and its
existing 15%/5% peak-relative limits. Neither failure falls back to one-way audio.

## State, work and supported combinations

Acoustic coordinates start at rest and follow all original mechanical,
Kelvin and material-memory coordinates. All original player ports and source
addresses remain fixed. The actual PHS interconnection supplies equal and
opposite internal power exchange; acoustic storage and damping enter the
existing TOTAL work/loss ledger exactly once. End-of-render acoustic energy
and instantaneous loss power on stderr are components, not extra debits.

Analytic storage tangents include the new coordinates. Finite-difference
preparation and bounded impact substeps solve the same loaded model. Mechanical
failure restores acoustic history together with the original contact/material
state, and consumes no staged stick/pedal force. Output clocks, decimators,
pressure scale and receiver propagation remain unchanged. The offline WAV
writer still publishes no partial file; a late audio failure is terminal for
that experiment, not an independently retryable audio-block transaction.

The existing separate material/stand damping is not silently removed. Supplied
modal damping should represent intrinsic structural loss rather than a measured
TOTAL decay already containing radiation. Adding this load to such a total-loss
fit would double-count radiation. No measured material identification is claimed.

## Explicit limits

`--radiation-spec` retains its bandwidth/refinement/work controls. Feedback
additionally requires the COMPLETE source basis to have at most 32 columns,
8..64 training intervals (33..257 total BEM samples), at most 256 factored poles,
and at most 1024 total scalar mechanical/internal/acoustic states. Excess
work refuses, never truncates ranks or retained modes. The dictionary includes
storage poles up to eight times the requested top frequency; every pole must
satisfy the original mechanical-rate guard. The coupling norm times the time
step must not exceed 0.25, as in the existing piano load runtime. These are
resolution/work guards, not temporal error bounds or real-time qualification.

The same original BEM work cap bounds its solves; matrix projection and fitting
are additionally bounded by the above source/sample/pole limits. Cancellation
is checked between frequency solves and fits. Existing dense BEM/fitting kernels
are not made internally interruptible by this composition.

This is linear, stationary-reference exterior loading fitted over a declared
band. There is no moving-boundary radiation, new room/microphone model, full-band
certificate or measured specimen calibration. Hi-hat closure retains one
reference separation: squeeze-film air and gap-dependent scattering remain
absent. Sampled fit gates do not establish interpolation/extrapolation accuracy.

Focused native checks:

```sh
cargo test --release -p fs-phs --lib port_load
cargo test --release -p fs-couple --lib render::plate::impact::radiation
cargo test --release -p fs-couple --example percussion radiation_feedback -- --test-threads=1
```
