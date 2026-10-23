# Physical sensitivity through actual cooling solves

`cooling-network-uq BASE UQ --sensitivity shapley` attributes the variance of the
cooling response to its original uncertain physical parameters, including their
declared dependence. It uses `fs_uq::ShapleyExecution` conditional Gaussian
pairs for Gaussian dependence, or physical pick-freeze pairs for explicitly
independent Gaussian/uniform mixtures. It never substitutes independent swaps
for correlated sampling or reports surrogate-germ Sobol coefficients.
Every callback still invokes the existing cooling product and its watchdog.

```sh
cargo run -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-shapley-hotspot.json --sensitivity shapley
```

The example declares correlated inlet temperature and component power. Its
24 physical calls supply four base rows of six solves each. This tiny design
illustrates the interface, not a precision or hardware qualification. For d
inputs, samples MUST equal N*2*(2^d-1), with N>=2; nothing is rounded, added or
silently discarded. The existing 10,000-physical-call ceiling applies. Up to
eight varying physical inputs are supported. Independent laws may mix Gaussian
and uniform marginals; joint-gaussian laws require all-Gaussian marginals and
support admitted singular correlations. Keep constants in the base model.
Zero uncertainty, unknown dependence, correlated uniforms and invalid conditional
factors refuse before any physical work or file output.

The original target lowering, names and units are retained. Draws outside the
physical input domain fail the experiment; they are not clipped or replaced.
Steady objectives and explicit transient-sampled-peak observations both use the
same adapter. Each transient sample is one COMPLETE declared trajectory,
including fixed repeated cycles, not an endpoint treated as an independent draw.

## Interpretation

The completed `frankensim.cooling-network-uq.shapley.v1` report includes physical
parameter `effects` with dimensionless `index`, `base_output_std_dev_k`, original
parameter laws, solver counts, and `conditional_variance_fractions`.
The latter array is indexed by a bitmask in the original parameter order:
bit j means physical parameter j is conditioned on. Entry S estimates
E[Var(response | X_S)] / Var(response). The empty-set fraction is one; the full
set is zero under the fixed deterministic-model assumption. Every other entry
comes from actual paired cooling evaluations.

Shapley effects share dependence and interaction contributions. They are not
causal effects, local derivatives, independent-input main/total Sobol indices,
or sensitivity of a probability inferred by thresholding the response. A
configured temperature limit does not change this response into an indicator.
Finite-sample effects can be negative or exceed one and are not clipped.
Their sum being one is an algebraic identity, not evidence of accuracy.
Undefined normalization emits a diagnostic, not invented zero effects.
There is no confidence interval, convergence decision or continuous-time peak
bound. The original `--sensitivity sobol` behavior is unchanged.

## Uniform tolerances mixed with Gaussian scatter

`uq-shapley-mixed-hotspot.json` declares Gaussian inlet scatter and independent
uniform component-power tolerance. Both laws reach the original physical model:

```sh
frankensim --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-shapley-mixed-hotspot.json \
  --sensitivity shapley --workers 3
```

The 24 calls are four complete six-evaluation rows, an interface example rather
than a precision claim. No uniform is approximated by a Gaussian or sampled at
its midpoint. Pair coordinates retain the declared marginal, name and unit;
independence must be explicitly declared, not inferred from bounded supports.
Epistemic bounds do not automatically define a uniform probability distribution.

Mixed and all-uniform designs report `method: conditional-product-philox`.
Existing all-Gaussian designs keep `conditional-gaussian-philox`, their original
samples and their original arithmetic. The effect definition and `shapley.v1`
schema are unchanged. Serial and parallel execution use the same producer;
their existing checkpoints bind each uniform endpoint as well as every Gaussian
parameter. Changing a support or law on restart refuses rather than reweighting
previous observations. The recovery and failure rules below apply to both.

Native checks for this extension:

```sh
cargo test -p fs-uq --test product_shapley --test product_shapley_mixed
cargo test -p fs-cli --test cooling_shapley_mixed
```

## Recover every paid conditional half-pair

```sh
frankensim --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-shapley-hotspot.json --sensitivity shapley \
  --checkpoint first.shapley --max-new-samples 3

frankensim --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-shapley-hotspot.json --sensitivity shapley \
  --resume first.shapley --checkpoint completed.shapley
```

The first command exits BUDGET after three accepted physical solves, halfway
through the pair for subset mask one. Progress reports include next row,
subset mask and half (0=first, 1=second), but NO sensitivity estimate or fractions.
The second resumes the exact unpaid partner without repeating those three solves.
Wall-time cancellation also retains completed work when a destination is supplied.
An unfinished transient retries the same full trajectory from its initial state.
Without a checkpoint, timeout emits only a diagnostic, never shortened-design
statistics. A complete checkpoint performs no additional physics on resume.

The existing bounded atomic transport saves after every accepted call. Output
paths must be new and source checkpoints are never overwritten. Filesystem errors
refuse; a failed physical solve invalidates the current output so its earlier
successful prefix cannot be silently resumed past the failure. `FSSHP001` binds
the original allocation, every input-law field, conditional ordering, model,
base-request bytes, target lowering and executable. Reordered inputs, a changed
model or seed, another sampler's checkpoint and corruption refuse before output
reservation. Only per-invocation wall time and chunk length may change.
Checksums are not authentication: use trusted files and the same deterministic
runtime profile. The library exposes `ShapleyExecution::checkpoint/restore` for
other consumers without requiring cooling or the CLI.

Other sampling modes, mean controls, sequential compliance and candidate
selection are incompatible and refuse before reading inputs. No checkpoint
changes or new dependencies are required for existing MC/QMC/Sobol consumers.

Focused native checks:

```sh
cargo test -p fs-uq --test product_shapley_checkpoint
cargo test -p fs-cli --bin frankensim uq_command::execute::shapley
cargo test -p fs-cli --test cooling_shapley
```

The tests cover native physical reconstruction, complete-transient observations,
chunked continuation, changed identities, failed physics and incompatible modes.
They must run natively before claiming a Rust execution pass or physical accuracy.

## Concurrent correlated sensitivity

Add `--workers N`, with N in 1..64, to evaluate original conditional pairs using
concurrent cooling children. This is the same physical experiment, not latent
Sobol sensitivity or a replacement independent-input law. The complete result
uses the original `shapley.v1` renderer and ordered reductions. Worker identity,
completion order and wave width do not change the conditional samples.

```sh
frankensim --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-shapley-hotspot.json --sensitivity shapley \
  --workers 3 --checkpoint first.shapley-batch --max-new-samples 3

frankensim --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-shapley-hotspot.json --sensitivity shapley \
  --workers 1 --resume first.shapley-batch --checkpoint complete.shapley-batch
```

A wave may split pairs or cross coalition and row boundaries. Completed later
siblings survive interruption of an earlier call, including on disk, and only
missing ordinals are retried. `shapley.parallel.progress.v1` exposes contiguous
and out-of-order counts plus the next missing row/subset/half, but no sensitivity
estimate. A completed checkpoint requires no further physical solves. Actual
model failures, nonfinite outputs and caught worker panics remain terminal.

`FSSHPB01` embeds the original conditional-prefix checkpoint and bounded pending
slots. Repeat both selectors on resume, changing only worker count, wall time
or chunk allowance. Serial `FSSHP001` and batch formats are not interchangeable.
The existing atomic file transport saves AFTER each drained wave; a crash can
lose its unsaved wave but not the previous published checkpoint. The worker cap
is not an aggregate memory/CPU cap: solver resources are still per child. See
`PARALLEL_UQ.md` for shared resource and draining semantics. No measured speedup,
confidence guarantee, causal interpretation or model-validation claim is added.

Other schedulers may use `fs_uq::product_shapley::ShapleyBatchExecution` directly.
Focused native checks, required before claiming a runtime/bit-parity pass:

```sh
cargo test -p fs-uq --test product_shapley_batch
cargo test -p fs-cli --bin frankensim uq_command::execute::parallel
cargo test -p fs-cli --test cooling_shapley_parallel
```
