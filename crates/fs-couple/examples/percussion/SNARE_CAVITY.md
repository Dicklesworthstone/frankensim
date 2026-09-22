# Snare wires and distributed air in one mechanical solve

`--cavity-modes` now works with `snare`, `snare-off`, and `drum-modal`, including
all their existing `-wav` and `-mic` variants. It extends the original nonlinear
path in `CAVITY.md`: the old restriction to `drum` and `drum-stretch` no longer
applies to sealed prepared cavities. Commands without the option are unchanged.

```sh
# Mechanical onset: 20 distinct strands, two heads, spatial enclosed air.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  snare 4096 --cavity-modes --strike-speed-m-s 4 \
  --strike-position-m 0.06 0.01 > snare-cavity.csv)

# The same physical system observed through the existing exterior BEM path.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  snare-mic 4800 2 0.08 0.05 0.35 --cavity-modes \
  --drum-spec crates/fs-couple/examples/percussion/estimated_drum.fsd \
  --strike-position-m 0.06 0.01 > snare-cavity.wav)
```

Without `--head-stretching`, these commands retain the prepared linear-head
modal/contact image and reject `--prepared-nonlinear`. Selecting
`--head-stretching` explicitly combines both nonlinear heads with the full wire
bank and its existing contact loss; `--prepared-nonlinear`, `--analytic-newton`
and `--impact-substeps` then operate on that joint nonlinear system. See
[NONLINEAR_SNARE.md](NONLINEAR_SNARE.md). Both images admit the same resistive
neck for mechanics CSV. Vented exterior audio still refuses. No nonlinear
potential or declared loss is silently dropped.

## Actual coupling, not a pressure overlay

The original film FEM, installed tensions, damping, striker, wire geometries,
and contact laws are unchanged. The reference snare retains all 20 strands,
160 wire coordinates and 240 distributed wire/head contact points, plus the
stick contact. Each strand still moves independently. The cavity basis and head
surface quadrature are the same ones used by the nonlinear drum. A shared
construction supplies the uniform compression mode and nonuniform pressure
modes; it does not install a second copy of the compact gas spring.

Each acoustic spring acts through its complete signed column: both heads plus
its air inertia. That column cannot be replaced by independent pairwise springs
without changing the Hamiltonian's cross terms. The new reusable
`CavityCoupling::build_linear` compiles a direct sum of the unchanged diagonal
coordinates into the existing prepared connection/contact solver. No eigenbasis
is changed, no wire is homogenized, and no new time integrator is introduced.
All cavity and contact reactions are simultaneous, not previous-sample forcing.

Striker, heads and wires keep their original state addresses; acoustic inertia
is appended last. Wires and striker have zero direct enclosed-volume coupling.
Wires affect air through the resonant head. External audio still observes actual
head motion through the original closed exterior boundary and BEM transfers.
Neither interior pressure nor an extra wire/gas oscillator is summed into the
microphone signal. Direct wire radiation, radiation backreaction and room
scattering remain outside this model.

## Vents and acoustic losses

`--cavity-neck radius_m effective_length_m resistance_Pa_s_m3 azimuth_rad z_m`
now works with **mechanics CSV** for `drum-modal`, `snare`, and `snare-off`, as
well as `drum` and `drum-stretch`. It requires `--cavity-modes`. The existing
finite-aperture pressure averages, compactness guard, gas inertia, initial zero
neck volume/flow, and explicit nonnegative resistance are retained. Pressure
includes outward displaced neck volume; loss includes the actual `R*Q^2`
resistance. The four existing neck CSV columns report volume, flow, pressure
and endpoint resistive power. Endpoint power is not the step's energy loss.

`--cavity-drag-per-s D` explicitly supplies the same momentum-drag rate [1/s]
to each retained **nonuniform** cavity coordinate, in either mechanical image.
It requires `--cavity-modes`, accepts finite `0..=100000`, and defaults to zero.
The uniform pressure mode remains undamped because it has no inertial state.
This is a declared constant-rate model, not a measured wall/thermal loss,
pressure-decay fit, linewidth, or frequency-independent material loss factor.
Drag removes `D*p_j^2` from each mass-normalized gas momentum; it never applies
an output envelope or silently changes head/wire damping. A sealed render may
use it: the existing BEM observer receives the resulting changed head motion.

```sh
# Illustrative inputs, not measured vent geometry or identified acoustic loss.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  snare 20000 --cavity-modes --cavity-drag-per-s 50 \
  --cavity-neck 0.005 0.012 5000 0.4 0.08 \
  --strike-speed-m-s 2 --second-stick-position-m -0.05 0.02 \
  --second-stick-speed-m-s 1.6 --muffler batter 0.07 0.01 0.1 > vented-snare.csv)

# Sealed audio with acoustic momentum loss; no vent is requested.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  snare-mic 48000 20 --cavity-modes --cavity-drag-per-s 50 > lossy-snare.wav)
```

Both stick-force performance files, supplied drum specifications and fixed
mufflers continue to compose with these controls. All resistance, contact and
cavity reactions use one accepted mechanical step. **Any vented WAV/microphone
request still refuses** before geometry: the closed exterior omits aperture
radiation and cannot honestly represent a vented instrument's sound.

## Budgets and evidence

The prepared snare keeps its existing limits of 256 total coordinates, 512
contact points and 50 million contact setup terms. Air coordinates count toward
that same total; overflow refuses rather than truncating the wire bank. The
explicit nonlinear-head path also admits up to 256 total coordinates, retaining
the complete wire bank instead of truncating to the former 64-mode ceiling.
Neck additions preserve the caller's original total-coordinate budget,
including across repeated additions. No wire is removed to make room for a vent.

Prepared commands with a vent or positive acoustic drag declare 32 bilateral
connections and one million connection-setup terms, enough for eight cavity
springs, seven gas drags, one vent drag and sixteen solid mufflers. The original
zero-loss/no-vent commands keep their existing eight-connection envelope.
Library builders never enlarge caller budgets: every nonzero inertial drag
consumes one grounded viscous link in addition to springs and supplied ports.
These links use the existing simultaneous solve, not previous-sample forcing
or a separate decay stage. This is passive second-order coupling, not the exact
full damped exponential; strong drag still requires time-refinement checks.

The original library tests cover uniform-volume equivalence, independent coupled
eigen-dynamics under time refinement, basis-rescaling invariance, many-body
contact and cancellation/retry, and explicit loss/budget admission. Five original example
tests cover real 20-strand assembly, changed head motion and spatial pressure,
prepared/reference onset comparison, supplied-geometry audio construction, and
command/physics restrictions. The audio-construction test does not execute BEM
fitting or certify a waveform. Focused native commands are:

```sh
cargo test --release -p fs-couple --lib render::plate::impact::cavity::prepared::tests
cargo test --release -p fs-couple --lib render::plate::impact::linear::free_drag_tests
cargo test --release -p fs-couple --example percussion cavity::loss_tests -- --test-threads=1
cargo test --release -p fs-couple --example percussion -- --test-threads=1
```

The conservative equations were also evaluated independently in Python: halving
the test step reduced the trajectory error by about fourfold. That is not a
native Rust result. Native builds, full-band/modal convergence, specimen
calibration and real-time deadlines require their own execution and evidence.
Additional native regressions cover damped Helmholtz flow/loss and refinement,
132-coordinate vent admission, lossy pressure-basis rescaling, inertial force
and work scaling, and the actual two-stick/20-strand/muffler/vent composition.
Native tests and audio exports remain unexecuted in the editing environment;
independent arithmetic checks are not evidence that the Rust tests passed.
