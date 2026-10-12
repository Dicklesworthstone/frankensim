# Signed-DWR checks inside native topology studies

`bracket-dwr-checked-2d.fsim` uses the existing `projected-volume` study and
`ProjectedOptimizer`. It adds a numerical acceptance check, not another
optimizer, study format, checkpoint format, or physics implementation.

```sh
cargo run -p fs-cli -- --json study \
  examples/marquee/bracket-dwr-checked-2d.fsim dwr-study.db --budget 1

# Substitute the run_id returned by study; use the identical executable.
cargo run -p fs-cli -- --json study --resume study-<receipt-hash> dwr-study.db
cargo run -p fs-cli -- --json report study-<receipt-hash> dwr-study.db
```

The optional declaration follows `:cg-poll-iters` (and `:refine-from`, when
present), before prescribed design regions:

```lisp
:goal-check (elasticity-dwr
  :max-error-fraction 0.5
  :max-solves-per-attempt 34)
```

The fraction must be finite and strictly between zero and one. The solve
allowance must cover `2 * (1 + max-candidates)` additional solves; with sixteen
candidates that is 34. Values above 130 are refused. The native geometry level
remains 2 through 5, with one additional enriched level. The study requires at
least 256 MiB admitted memory before starting geometry or solver work.

This declaration is currently for the single-load, hard-area volume mode.
Combining `:goal-check` and `:mesh-check` is explicitly refused rather than
silently bypassing either policy. Stress-constrained and multi-load goal gates,
spatially adaptive refinement, and 3-D topology are not supplied by this change.
The separate final `(assessment ...)` remains available and does not replace
in-loop admission.

## What the check evaluates

The accepted baseline is assessed once per attempt. Each otherwise admissible
coarse candidate is assessed before the existing optimizer publishes its
geometry, iteration ordinal, or augmented-Lagrangian multiplier. Both DWR
solves use the same original bilinear level set and canonical elasticity
problem; enrichment does not add material, re-project, or redistance it.

Let `g0 = J0_baseline - J0_candidate` and
`g1 = J1_baseline - J1_candidate`. In addition to strict requested relative
improvement on **both** grids, admission requires

```text
abs(eta_signed_baseline - eta_signed_candidate) <= fraction * g0
abs(g1 - g0)                                   <= fraction * g0
```

Every retained coarse measurement must reproduce the original accepted or
projected snapshot, material area and compliance bit for bit. The recomputed
relative residuals on both grids must be finite and below `1e-12`. A rejected
candidate contracts the existing bounded search. Exhaustion returns
`no-feasible-descent`, not convergence or success. A baseline numerical failure
retains the feasible design and returns `numerical-failure`.

**These are Estimated numerical screens, not continuum error bounds.** Shared
bias can survive agreement. The material constraint remains coarse-grid cut
quadrature. Absolute DWR indicator mass is not used as an error radius. Passing
does not establish an optimum, stress admissibility, or physical validation.

## Retention, interruption and reporting

JSON and HTML reports include `constraints.goal_resolution`. The last complete
search retains baseline estimates and, when accepted, candidate estimates and
recomputable gain arithmetic. Terminal candidate refusals use the existing
`terminal_refusals`. Resume checks measurement identity, residuals, policy,
gain arithmetic and work counters before recovery physics.

Each attempt's **worst-case additional solve allowance** is persisted before
work. `attempts_admitted` and `total_reserved_solves` are conservative charges,
not measured completed-solve counts. Interrupted or crashed attempts are not
free, and a retry receives a new charge. The last complete numerical check and
last accepted geometry survive interruption. Existing wall-time accounting
continues independently; the original proposal/projection work is additional
to the declared DWR solve allowance.

Cancellation is checked before and after each DWR assessment, as well as at the
existing optimizer checkpoints. DWR assembly, its two solves and residual
integration remain indivisible; this feature does not promise intra-DWR-CG
cancellation or a wall-time bound for those kernels.

## Native validation targets

```sh
cargo test -p fs-topols signed_dwr_gate_ -- --nocapture
cargo test -p fs-cli signed_dwr_native_ -- --nocapture
```

The regression sources cover signed corrections, strict and floating-point
boundaries, load scaling, hostile evidence, full-family work admission, actual
solver cancellation, exact retry, native reporting and retained-data tampering.
They were authored without a Rust toolchain in the editing environment; their
native execution remains pending. This feature is an in-loop slice of
`frankensim-rc-root-q61wp.16`, not closure of its broader acceptance requirements.
