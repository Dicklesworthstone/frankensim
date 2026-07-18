# CONTRACT: fs-rep-neural

Neural implicit charts: small coordinate MLPs as shapes, Lipschitz-constrained
so certified bounds remain available.

## Purpose and layer

Layer L2 (MORPH / representation). Pure Rust; depends on L1 `fs-ivl` for
outward-rounded interval arithmetic and on `fs-math` for the exact deterministic
elementary-function implementation covered by those ULP budgets. Spectral
diagnostics and upper bounds remain in-house.

## Public types and semantics

- `Layer::new(weights, bias)` — a dense affine layer; `spectral_norm(&weights)`
  (power iteration on `WᵀW`); `spectral_normalize(layer, bound)` (scale so the
  certified spectral-norm upper bound is at most `bound`).
- `MlpSdf::new(layers, bound)` — spectrally normalizes each layer to `bound` and
  records the certified global Lipschitz constant `L = Π Uᵢ`, where each `Uᵢ`
  is an outward-rounded spectral-norm upper bound (tanh is
  1-Lipschitz). `eval` (tanh hidden, linear output), `eval_grad` (finite
  differences), `eval_interval(lo, hi)` (IBP output enclosure), `lipschitz`,
  `input_dim`, `topology_hint`. The corresponding `try_eval`, `try_eval_grad`,
  and `try_eval_interval` entry points return a structured
  `InputDimensionError` instead of panicking at an untrusted boundary.
- `safe_step_radius(value, lipschitz)` — `|value|/L`, the provably safe
  sphere-tracing step.
- `TopologyHint::Unknown` — the only variant; topology is never inferred from
  the fit.

## Invariants

- The certified Lipschitz constant is a valid UPPER bound: no sampled pair
  violates `|f(x) − f(y)| ≤ L·‖x − y‖`.
- IBP is SOUND: `eval_interval(lo, hi)` encloses `f(x)` for every `x` in the box.
  Every affine product/sum is outward-rounded and every hidden `tanh` uses
  `fs-ivl`'s deterministic five-ULP enclosure. Point evaluation uses the same
  `fs_math::det::tanh` primitive, exposed as
  `MLP_ACTIVATION_SEMANTICS=fs-rep-neural-det-tanh-v1`, rather than an
  ungoverned platform `tanh`. A degenerate input box may widen by the accumulated
  rounding budget but must contain the separately evaluated point.
- Point, gradient, and both interval endpoint vectors must contain exactly
  `input_dim()` coordinates. No evaluator truncates, pads, or otherwise
  reinterprets a malformed vector; all four checks share one admission helper.
- The analytic gradient of the continuous real MLP satisfies `‖∇f(x)‖ ≤ L`
  everywhere. `eval_grad` is only a rounded central-finite-difference
  diagnostic: its coordinate secants use different line segments and it has no
  gradient-certificate authority.
- A sphere-trace step of `safe_step_radius(f(x), L)` never tunnels: `f` cannot
  change sign within that radius.
- `topology_hint` is always `Unknown` (honest — never claimed from the loss).

## Error model

Layer construction panics on structural misuse. The convenient `eval`,
`eval_grad`, and `eval_interval` methods also panic on an input-dimension
mismatch; their `try_*` counterparts return deterministic structured errors
containing the evaluation surface plus expected and actual dimensions. A
non-finite or inverted correctly dimensioned interval box returns
`(-infinity, +infinity)` as the fail-closed enclosure.

## Determinism class

Fully deterministic: the spectral norm uses a fixed initial vector; point and
interval activation semantics share versioned `fs-math`/`fs-ivl` arithmetic;
eval, IBP, and the Lipschitz constant are pure functions of the weights.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/neural.rs`: spectral norm vs known values; spectral
normalization to a bound; the Lipschitz certificate is never violated; IBP
soundness, degenerate-point enclosure, malformed-box refusal, and deterministic
endpoint replay; exact point-evaluator binding to the interval certifier's
deterministic `tanh`; short, long, and empty vectors across point, gradient, and
both interval endpoints with exact structured diagnostics; deterministic
finite-difference diagnostics with an explicitly non-authoritative generic
bound; a certified sphere-trace step never tunnels; topology honestly unknown;
determinism.

## No-claim boundaries

- The certificate machinery is complete for ANY spectrally-normalized weights;
  FITTING (DeepSDF-style training from source charts via FrankenTorch, eikonal
  regularization, latent-code conditioning for shape FAMILIES with exact
  autograd Jacobians) is the fuller deliverable, staged — this v0 does not train.
- IBP is the interval evaluator; the tighter CROWN-class linear-relaxation bound
  propagation is a follow-on.
- `eval_grad` is a finite-difference diagnostic, not a certificate. An analytic
  AD gradient or interval-derivative enclosure bound to its arithmetic and
  error model is required before derivative evidence may carry authority.
- Watertightness and Hausdorff agreement vs the source chart come from the
  certificate machinery (fs-rep validity-certificates), never from this crate.
