# CONTRACT: fs-sos

Proof-carrying optimization: executable sum-of-squares decomposition checks and
soundly scoped polynomial lower bounds.

## Purpose and layer

Layer L4 (ASCENT). Production code is safe Rust and depends on `fs-ivl` for
exact expansion arithmetic and `fs-math` for error-free products. The PSD path
uses an in-house Jacobi eigensolver. `fs-obs` is test-only evidence plumbing.

## Public types and semantics

- `Poly` is a univariate polynomial with ascending coefficients. It provides
  construction, Horner evaluation, arithmetic, degree/coefficients, and a
  coefficient infinity norm; `square(q)` returns `q * q`.
- `SosCertificate { squares, lower_bound }` represents the proposed
  decomposition `p - lower_bound = sum(squares_i^2)`.
  - `residual(p)` is the floating coefficient infinity norm of the mismatch.
  - `verify(p, tol)` is only a coefficient-residual diagnostic. A positive
    coefficient tolerance is not a global value theorem.
  - `certified_bound_global(p)` returns a global bound only when every
    non-constant residual expansion is exactly zero. It absorbs the enclosed
    constant residual into the returned bound.
  - `certified_bound_on(p, radius)` encloses every residual term and returns a
    sound bound for `|x| <= radius`; a mismatch degrades the bound instead of
    repeating an unsupported claim.
- `certify_quadratic(a, b, c)` constructs the usual completed-square
  certificate for finite `a > 0` when all derived values remain finite.
- `is_psd(matrix, tol)` checks the symmetric part of a square matrix against
  minimum eigenvalue `-tol`; it is the current SDP-feasibility core.
- `lyapunov_certifies_stability(A, P)` verifies the fixed two-dimensional
  quadratic Lyapunov inequalities for a supplied `P`.

## Invariants

- A value returned by `certified_bound_global` is sound for every real `x`.
- A value returned by `certified_bound_on` is sound on its stated finite
  radius, including for a mismatched or overstated input certificate.
- `verify(p, tol)` claims only bounded coefficient mismatch; callers must not
  promote that diagnostic into a value bound.
- The exact dyadic fixture `x^2 - 2x + 3 = (x - 1)^2 + 2` produces an exact
  global bound of `2` and is covered by a bit-complete replay receipt.
- PSD and Lyapunov decisions are made from the symmetric quadratic form, so an
  asymmetric matrix cannot forge a certificate through ignored entries.

## Error model

- `certify_quadratic` returns `None` for non-finite input, `a <= 0`, or
  non-finite derived certificate values.
- `certified_bound_on` returns `None` unless `radius` is finite and positive;
  `certified_bound_global` returns `None` for any nonzero non-constant exact
  residual.
- `is_psd` expects a square, consistently sized matrix. Ragged input is outside
  this v0 API's admitted domain and may panic; shape-typed admission is staged.
- Polynomial arithmetic follows ordinary `f64` behavior for non-finite values.

## Determinism class

Operations are deterministic for the same inputs, build, and ISA. The G5
fixture binds every returned bit and replays exactly on the current build. This
contract does not claim cross-ISA bit equality for square root or eigensolver
paths.

## Cancellation behavior

None. Current operations are finite, synchronous functions without `Cx`.

## Unsafe boundary

None. Workspace lints deny unsafe code.

## Feature flags

None.

## Conformance tests

- `tests/sos.rs` covers polynomial arithmetic, quadratic and multi-square
  certificates, overstated and bogus certificates, exact-global and
  radius-scoped bounds, the historical tolerance-forgery counterexample,
  invalid quadratic input, symmetric-form PSD behavior, Lyapunov verification,
  and deterministic repetition.
- `tests/quadratic_study_replay.rs` is a G5 exact-dyadic production fixture. It
  binds the complete `certify_quadratic` result and derived public verdicts,
  checks retained schema-v1 fixture/result roots, requires byte-identical
  in-process replay, emits wire-valid `fs-obs` evidence, and catches a disclosed
  seeded one-bit square-coefficient mutation at payload, retained-reference,
  semantic, and merge gates.

## No-claim boundaries

- The replay fixture proves one finite exact-dyadic quadratic and one disclosed
  mutation lane. It is not exhaustive tamper testing or cryptographic
  authentication; the current replay root is a non-cryptographic house digest.
- `certify_quadratic` is not claimed expansion-exact for every floating input,
  and `verify(p, tol)` is never a global theorem merely because it passes.
- General univariate or multivariate certificate search, Lasserre/moment
  relaxations, Burer-Monteiro SDP optimization, and Positivstellensatz
  machinery are staged. The crate currently checks supplied certificates and
  constructs only quadratics.
- `lyapunov_certifies_stability` verifies a supplied two-dimensional `P`; it
  does not search for `P` or certify nonlinear regions of attraction.
- Under `docs/CERTIFICATE_REGIMES.md`, this is only a local-stability route
  inside the stated model, equilibrium, parameter domain, and Lyapunov
  assumptions. It cannot be widened into global attraction, long-horizon
  predictive accuracy, broadband validation, or duty-cycle reliability.
- No cross-ISA replay, cancellation/concurrency, persisted or authenticated
  ledger, external-oracle, broad-input, or performance claim is made here.
