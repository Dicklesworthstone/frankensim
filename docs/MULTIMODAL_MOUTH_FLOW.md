# Mouth flow from the original multimodal bore

`fs_duct::modal::outlet` exposes `mouth_response`, `mouth_response_tabulated`,
`mouth_response_load` and `mouth_response_flanged`. Each returns the original
input `ModalResponse` and the complete row-major map from inlet modal volume
flow to mouth modal volume flow. Mode zero is total volume flow; higher modes
retain the same wall-normalized radial basis as the original modal owner.

The calculation runs inside the existing backward reflection recurrence.
At the terminal plane it retains `H = Zc^-1 (I-R)`. Propagating toward the inlet
updates `H <- H E` with the original decaying `E = diag(exp(i k L))`. Junction
pressure and flow continuity give the downstream forward amplitudes from the
upstream ones; `H` composes that map. It never evaluates a growing evanescent
exponential, inverts the projection matrix, or solves a substitute scalar bore.
Only one additional N-by-N map is retained, not a station-by-station field.
The impedance-only entry points disable this extra calculation and preserve
their original arithmetic and public result type.

`MouthResponse::per_incident_plane_pressure(reference)` gives every mouth-mode
flow per unit incident plane pressure at the specified positive real reference.
It imposes zero higher-mode INLET flow, matching the Z00 closure used by
`MmLineBank`. This is different from driving the full reflection matrix with
zero higher-mode incident waves. The positive-frequency table path validates
its original mouth plane, radius and gas; a separate explicit-load evaluation
can use declared static impedance without extending the table below its band.

The flanged entry point retains the complete original exterior radiation matrix
and its interior/exterior-gas admission. Closed and ideal pressure-release ends
impose zero flow and zero pressure on every retained mode. Scalar radiation,
explicit scalar impedance and tabulated loads retain their matched higher-mode
mouth approximation. Cone endpoints are still the
original midpoint-radius staircase image: returned endpoint radii identify the
actual discrete stations, and convergence still requires refining that image.
This change neither alters that geometry nor certifies arbitrary imported loads
as causal. Mouth flow is not a microphone signal or an exterior radiation model.

Rust regressions cover analytic resistive cylinders, long evanescent spans,
scalar-chain degeneration in both junction directions, an independent
simultaneous all-wave system, flanged input/output power, strict-table/DC
separation, and zero full-mode flow through a rigid cap after radius changes. Native Rust compilation and test execution remain unverified in the
implementation environment. An independent Python/SciPy derivation check passed
160 simultaneous-system comparisons (1-4 modes, both radius-step directions,
complex propagation and above/below-cutoff frequencies), with worst relative
flow error 2.03e-15, plus a long-evanescence case. That checks the equations, not
the compiled Rust implementation or experimental instrument accuracy.

Run the native regressions with:

```sh
cargo test -p fs-duct --lib modal::outlet::tests
```

## Scheduled modal mouth-flow runtime

`fs_couple::brass_sequence::mouth::BrassMouthSequence` constructs the original
`BrassSequence` and its output transfers together from the same ducts, gas,
loads, modal resolution and controls. It does not attach an unrelated output
model to an already-running voice. The input impedance and output transfer
therefore share the original modal bore and physical terminal law.

The runtime admits the original voice's analytic terminations and strict
referenced scalar tables. Closed and ideal pressure-release conditions act on
every mode; scalar radiation and tables retain matched higher-mouth modes. The frequency-domain full flanged-matrix API above is not yet
an input-load option for `BrassVoice`. Diagnostic table splices and matched
continuations are refused by the new runtime. A strict table covers every
positive output-ladder bin; DC uses its separately declared static impedance
at the same 1e-3 rad/s probe as the input bank, never table extrapolation.

`block(mouthpiece_pa, outgoing_pa, mouth_modes_m3_s)` returns three aligned
physical streams. The first two buffers have M samples; the third has M*N
sample-major values, N = `mode_count()`. All mode channels are retained. Mode
zero is total outlet flow, not microphone pressure. Every path's filters receive
the actual outgoing history from sample zero, including inactive paths. Hard
and faded switches use the original event indices and fade ramp; a same-path
assignment remains a no-op. Paths must share the same actual discrete aperture
radius before modal flow vectors can be crossfaded. This remains the existing
linear-path valve interpolation, not a simulation of valve travel or aperture
motion. It does not replace the held-bore lip-file microphone adapter.

`MouthFlowConfig` bounds each FFT ceiling (512..=65536) and the sum of all
path/mode ceiling taps before allocation. The original input bank retains its
separate limits. Each mouth mode must pass its own last-quarter energy AND
common-prefix tap-change gates. Scaling precedes squaring, so tiny but nonzero
modes are not erased by underflow. Exact-zero channels pass after a second
rung. The existing time-image owner supplies conjugation, Nyquist shaping and
precursor treatment; there is no clipping, fitted gain, or normalization.
Exhausted ladders refuse rather than silently truncating. Each realized filter
stores its taps and an equally long history; all paths are evaluated during
rendering. These bounds do not imply real-time throughput.

Callback size/horizon errors occur before advancement. Physical/filter failures
poison continuation because the source or an output prefix may have advanced;
no rollback is claimed. Observation never changes the lip/cup/source state.
Four runtime regression groups plus a scale-invariant convergence test cover
manual source parity, direct convolution across valve changes, callback
partitions, DC/zero modes, budget and buffer refusals, and failed continuation.
They remain natively unexecuted. Independent Python streaming arithmetic matched
direct convolution across 503 partitions, three paths and three nonzero modal
channels (maximum difference 5.56e-16), and checked mode scaling from 1e-200 to
1e200. Those checks do not execute the Rust acoustic realizer.

The existing manufactured example now exposes an optional physical CSV trace
alongside its mouthpiece-pressure WAV. Both files are create-only and the WAV
keeps its declared pressure scale; neither is presented as a measured instrument
or exterior microphone:

```sh
cargo test -p fs-couple --lib brass_sequence::mouth::
cargo run -p fs-couple --example brass_valve_phrase -- mouthpiece.wav mouth.csv
```

CSV columns are source sample, physical time, mouthpiece pressure, outgoing
inlet pressure, and every original mouth-mode volume-flow channel. Native example
execution is still unverified in the implementation environment.
