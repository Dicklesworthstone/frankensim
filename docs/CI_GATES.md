# CI Gates — the Decalogue as machine-checkable policy

Owned by the foundations CI/CD bead. This is the mapping from principle to
mechanical gate, with the honest wiring status of each. Phase gates are
Gauntlet states, not dates (plan §16.1); a gate listed as **deferred** below
is wired the moment its owning crate lands, and this file is updated in the
same change.

## Workflows

| Workflow | Trigger | What it proves |
| --- | --- | --- |
| `ci.yml` | PR + push to `main` | fmt, clippy `-D warnings`, workspace unit + conformance tests, xtask policy gates, constellation drift — on BOTH ISA families (plan §13.3: done = green on both) |
| `nightly.yml` | daily 09:17 UTC + manual | full suites in dev AND release profiles, policy gates, and a run record written into an fs-ledger database retained as an artifact (golden-ledger convention, §16.2) |
| `ci-self-test.yml` | manual | meta-tests: injected failures demonstrably turn the gates red |

## Runner honesty

- `macos-14` = Apple Silicon (aarch64, NEON): a true reference-family runner.
- `ubuntu-24.04` = x86-64, but GitHub-hosted runners do **not** guarantee
  AVX-512. fs-substrate resolves the SIMD tier at startup and the tier is in
  the logs; the x86 lane therefore proves *portable correctness plus whatever
  tier the runner exposes* (typically AVX2). An AVX-512-guaranteed lane needs
  a self-hosted runner and stays **deferred** until one exists. Roofline
  *performance* bands (§14) are meaningless on shared virtual runners and are
  explicitly NOT claimed by CI; they belong to the reference machines with
  ledgered fingerprints.

## Constellation

The workspace path-depends on sibling repos. CI materializes them with
`scripts/ci/checkout_constellation.sh` at the exact `constellation.lock`
pins (shallow fetch by SHA), then `xtask check-constellation` verifies zero
drift. Bumping a sibling is a deliberate act: re-run
`cargo run -p xtask -- lock-constellation` locally and commit the new lock.

## Decalogue mapping (plan patch Rev 25)

| Principle | Gate | Status |
| --- | --- | --- |
| P1 pure Rust, Franken-only | `xtask check-deps` + `check-unsafe` capsule registry lint | **live** (`ci.yml`) |
| P2 determinism | G5 audit: bit-identical re-run of seeded suites across thread counts | **deferred** — owning bead: determinism audit harness; unit-level determinism asserted in crate suites today |
| P3 differentiable-or-certifiable | gradient-check merge gate (adjoint vs dual/FD; a solver without a passing check cannot merge, §8.7) | **deferred** — wired when fs-ad adjoint infra + first solver land; fs-vskeleton already self-gates its adjoint in-run |
| P4 budgets first | Five-Explicits lint on ledger ops | **live in fs-ledger API** (structured rejection); IR-level admission lint deferred to fs-ir |
| P5 structure preservation | exact-sequence / conservation G0 law suites | **partial** — review-checklist status; suites land with fs-feec |
| P6 roofline-honest | perf-regression bands vs ledgered baselines (§14.4) | **deferred** — needs reference-machine runners + roofline harness; explicitly not faked on shared CI |
| P7 cancellation-correct | G4 storm gate (kill/cancel batteries, leak accounting) | **partial** — fs-ledger kill -9 battery runs in `ci.yml` today; executor-wide storms land with fs-exec |
| P8 one data model | conformance suites over shared complex/cochain types | **deferred** — with fs-geom/fs-feec |
| P9 provenance-complete | golden-ledger replay + integrity re-hash; constellation drift gate | **partial** — integrity + drift live; replay gate lands with fs-ledger time travel |
| P10 agent-first | JSONL verdicts from every suite retained as artifacts; structured-error schema validation | **partial** — verdict retention live; catalog no-drift gate lands with fs-ir |

## Gate meta-tests (are the gates real?)

`ci-self-test.yml` proves the wiring by injection:

1. `FS_CI_INJECT_FAILURE=1` arms a deliberately failing xtask test → the test
   gate must go red (asserted by inverting the step's exit code).
2. A seeded formatting violation → `cargo fmt --check` must go red.
3. The xtask unit battery seeds layer/dependency/contract/unsafe violations
   and asserts each is caught (`cargo test -p xtask -- seeded`).

A green `ci-self-test` run is the evidence that a red gate blocks; run it
after any change to the workflows or to xtask's checks.
