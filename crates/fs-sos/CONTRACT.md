# CONTRACT: fs-sos

Proof-carrying optimization: executable sum-of-squares decomposition checks and
soundly scoped polynomial lower bounds.

## Purpose and layer

Layer L4 (ASCENT). Production code is safe Rust and depends on `fs-ivl` for
exact expansion arithmetic and outward-rounded intervals and `fs-math` for
error-free products. The PSD path uses an in-house Jacobi eigensolver; the SDP
engine is in-house (no BLAS/LAPACK). `fs-obs` is test-only evidence plumbing.

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

### Multivariate proof-carrying layer (plan §9.8)

- `MPoly` / `IPoly`: multivariate polynomials with `f64` / interval
  coefficients (BTreeMap-ordered, exact zeros dropped, tiny coefficients never
  trimmed). `lie_derivative(V, f)` encloses `∇V·f` from exact float inputs.
- `sdp::solve`: primal–dual interior point (HKM direction, Mehrotra
  predictor–corrector) for block SDPs with exactly-handled free variables.
  Statuses: `Optimal`, `NearOptimal` (stalled within `tol_inaccurate`),
  heuristic `PrimalInfeasible`/`DualInfeasible`, `MaxIterations`, `Stalled`,
  `NumericalFailure`. Its output is UNTRUSTED.
- `SosProgram`: scalar / free-polynomial / SOS decisions and polynomial
  identities `constant + Σ multiplier·decision ≡ 0` (interval constants and
  multipliers), lowered to the block SDP. `solve_centered` pins scalars and
  maximizes a uniform Gram margin `t` (`Q ⪰ tI`, `t ≤ t_cap`).
- `verify(program, values) -> Certificate`: the only path to a claim. One SOS
  term with a constant nonzero point multiplier per identity is the slack; all
  other values are taken exactly; the residual is enclosed with intervals and
  absorbed into designated slack-Gram entries; positive definiteness of every
  Gram matrix is proved by the interval Cholesky method (Alefeld–Mayer 1993:
  feasibility of interval Cholesky ⇒ every symmetric member is PD).
- `minimize(p, constraints, opts) -> GlobalBound`: Lasserre/Putinar
  relaxation; `lower` is PROVED (backed-off γ re-solved centred and verified),
  `upper` is the interval-enclosed `p` at a feasible point extracted from the
  moments (Henrion–Lasserre column-echelon extraction, mean/eigenvector
  fallbacks, Newton polish when unconstrained). Unconstrained bases are pruned
  by coordinate/degree half-Newton-polytope tests plus iterative diagonal
  consistency (both sound).
- `certify_roa(f, opts) / certify_roa_with(f, P, ε, opts) -> RoaCertificate`:
  `V = xᵀPx` (default: Lyapunov equation of `Df(0)` with `Q = I`), `P ≻ 0`
  proved by interval Cholesky, and the S-procedure identity
  `−V̇ − ε‖x‖² − s·(c − V) = σ₀` proved by `verify` at each reported level;
  bisection only reports verified levels.

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
- A `Certificate` from `verify` implies: every identity holds exactly for real
  polynomials whose coefficients are the given floats (non-slack decisions)
  and some member of the enclosed slack Gram matrix, and every Gram matrix is
  positive definite. No SDP output is ever trusted without it.
- `GlobalBound::lower ≤ min p ≤ GlobalBound::upper` whenever both are present.
- `RoaCertificate::level` was verified: `V̇ ≤ −ε‖x‖² < 0` on
  `{xᵀPx ≤ level} \ {0}`, an ellipsoid, hence an inner estimate of the region
  of attraction of the polynomial model.

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

None. Operations are finite, synchronous functions without `Cx`; the SDP has an
iteration cap and ROA bisection a step cap, so every call terminates.

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
- `tests/multivariate.rs`: certified enclosures against a dense-grid oracle
  (nonconvex quartic) and the literature six-hump-camel minimum (G2, both
  minimizers extracted); the Motzkin falsifier (never SOS: must refuse, never
  claim); a Putinar disk bound (−√2); inflated-bound and forged-Gram
  falsifiers; bit-identical replay (G5); ROA soundness on `ẋ = −x + x³`
  (proved level strictly below the exact 0.5 boundary) and on the reversed
  Van der Pol oscillator, where RK4 trajectories from 64 boundary points of
  the certified ellipsoid all converge and `V̇ < 0` there, while the set stays
  inside the limit cycle; unstable/non-equilibrium refusals.
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
- The SDP engine is a dense second-order interior-point method sized for
  low-dimensional programs (Gram blocks up to a few hundred, Schur systems up
  to roughly a thousand rows). Burer–Monteiro low-rank first-order solving for
  scale, sparse Schur assembly, and chordal decomposition are staged.
- SDP infeasibility statuses are divergence heuristics, not verified Farkas
  certificates; a refusal from `minimize`/`certify_roa` asserts nothing about
  the problem.
- `minimize` uses a caller-fixed relaxation order and makes no claim of
  finite convergence; a large gap means raise the order. The upper bound's
  extraction is heuristic (its VALUE is enclosed, its tightness is not).
- `certify_roa` fixes a quadratic `V`; V–s alternation, polynomial Lyapunov
  functions, and rational/exact-arithmetic certificates are staged. The
  theorem is about the supplied polynomial model only.
- `lyapunov_certifies_stability` verifies a supplied two-dimensional `P`; it
  does not search for `P` (use `solve_lyapunov`/`certify_roa`).
- Under `docs/CERTIFICATE_REGIMES.md`, this is only a local-stability route
  inside the stated model, equilibrium, parameter domain, and Lyapunov
  assumptions. It cannot be widened into global attraction, long-horizon
  predictive accuracy, broadband validation, or duty-cycle reliability.
- No cross-ISA replay, cancellation/concurrency, persisted or authenticated
  ledger, external-oracle, broad-input, or performance claim is made here.
