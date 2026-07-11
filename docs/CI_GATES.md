# CI Gates — the Decalogue as machine-checkable policy

Owned by the foundations CI/CD bead. This is the mapping from principle to
mechanical gate, with the honest wiring status of each. Phase gates are
Gauntlet states, not dates (plan §16.1); a gate listed as **deferred** below
is wired the moment its owning crate lands, and this file is updated in the
same change.

## DSR is the CI source of truth

GitHub Actions is unavailable/throttled for this account, so FrankenSim uses
DSR as the primary CI and release runner. Agents should use DSR instead of
waiting on or citing GitHub Actions:

```bash
dsr repos info frankensim
dsr quality --tool frankensim
dsr build frankensim --target darwin/arm64
```

The GitHub workflow files remain in the repo as manual executable specs for the
gate shape and as historical documentation. They are not automatic push/PR
criteria.

## Workflows

| Workflow | Trigger | What it proves |
| --- | --- | --- |
| `dsr quality --tool frankensim` | manual DSR command | fmt, clippy `-D warnings`, workspace unit + conformance tests, xtask policy gates, constellation drift |
| `dsr build frankensim --target darwin/arm64` | manual DSR command | native release artifact build for the configured Apple Silicon lane |
| `ci.yml` | manual GitHub dispatch only | archived/manual version of the fmt, clippy, test, policy, and constellation gate shape |
| `nightly.yml` | manual GitHub dispatch only | archived/manual version of the full dev/release suites and retained fs-ledger run record |
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

## xtask policy gates (the `check-all` set)

`cargo run -p xtask -- check-all` runs, and CI treats as one gate:

| Check | What it refuses |
| --- | --- |
| `check-layers` | a crate depending upward/sideways against the L0–L6 map |
| `check-deps` | non-Franken runtime dependencies |
| `check-contracts` | crates missing CONTRACT.md canonical sections |
| `check-unsafe` | unsafe outside registered <300-line capsules with SAFETY.md |
| `check-powi` | optimization-level-dependent `f64::powi` on determinism paths |
| `check-goldens` | golden hashes whose upstream semantic surfaces drifted without a deliberate re-freeze (`golden-couplings.json`, docs/GOLDEN_POLICY.md) |
| `check-claims` | claim-state drift in the tracker mirror |

Each is also runnable alone (same names). Golden re-pins follow
docs/GOLDEN_POLICY.md: committed tree, BOTH build modes, plausible root
cause, coupling row updated in the same commit.

## Constellation

The workspace path-depends on sibling repos. The canonical fresh-checkout
path is `cargo run --manifest-path tools/bootstrap/Cargo.toml` (beads
huq.17 and 1t8i; docs/BOOTSTRAP.md). This standalone, zero-dependency package
builds before the root workspace resolves, fetches every repo from
`constellation.lock`'s recorded remotes at the pinned revisions, applies the
same pinned-head and clean-tree verification to existing and newly cloned
siblings, supports `--offline` cache verification and `--from` mirrors, and
writes fetch provenance. `cargo run -p xtask -- bootstrap-constellation`
remains the in-workspace command after the sibling paths already resolve.
`scripts/ci/checkout_constellation.sh` remains the shell-only
equivalent used by the manual workflow specs; `xtask
check-constellation` then verifies zero drift. Bumping a sibling is a
deliberate act: re-run `cargo run -p xtask -- lock-constellation`
locally and commit the new lock (schema v2 records remotes).

## Decalogue mapping (plan patch Rev 25)

| Principle | Gate | Status |
| --- | --- | --- |
| P1 pure Rust, Franken-only | `xtask check-deps` + `check-unsafe` capsule registry lint | **live** (`ci.yml`) |
| P2 determinism | G5 audit: bit-identical re-run of seeded suites across thread counts | **deferred** — owning bead: determinism audit harness; unit-level determinism asserted in crate suites today |
| P3 differentiable-or-certifiable | gradient-check merge gate (adjoint vs dual/FD; a solver without a passing check cannot merge, §8.7) | **deferred** — wired when fs-ad adjoint infra + first solver land; fs-vskeleton already self-gates its adjoint in-run |
| P4 budgets first | Five-Explicits lint on ledger ops | **live in fs-ledger API** (structured rejection); IR-level admission lint deferred to fs-ir |
| P5 structure preservation | exact-sequence / conservation G0 law suites | **partial** — review-checklist status; suites land with fs-feec |
| P6 roofline-honest | perf-regression bands vs ledgered baselines (§14.4) | **partial** — fs-roofline harness + regress statistics live (fz2.4/fz2.4.1); per-crate `#[ignore]`d perf lanes ledger attainment with anti-collapse floors; TARGET certification requires a quiet reference machine (shared-workbench noise swings axes ±15%+, measured), so gates report targets and assert floors until then |
| P7 cancellation-correct | G4 storm gate (kill/cancel batteries, leak accounting) | **partial** — fs-ledger kill -9 battery runs in `ci.yml` today; executor-wide storms land with fs-exec |
| P8 one data model | conformance suites over shared complex/cochain types | **deferred** — with fs-geom/fs-feec |
| P9 provenance-complete | golden-ledger replay + integrity re-hash; constellation drift + golden-coupling gates | **partial** — integrity, drift, `check-goldens` coupling discipline, and bootstrap fetch-provenance live; replay gate lands with fs-ledger time travel |
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
