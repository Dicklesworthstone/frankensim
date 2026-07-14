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
| `dsr quality --tool frankensim` | manual DSR command | configured fmt, clippy `-D warnings`, workspace unit + conformance tests, xtask policy gates, constellation drift, required quality lanes, x86 cross-check; see the external-wrapper limitations below |
| `dsr build frankensim --target darwin/arm64` | manual DSR command | native release artifact build for the configured Apple Silicon lane |
| `ci.yml` | manual GitHub dispatch only | archived/manual version of the fmt, clippy, test, policy, and constellation gate shape |
| `nightly.yml` | manual GitHub dispatch only | archived/manual version of the full dev/release suites and retained fs-ledger run record |
| `ci-self-test.yml` | manual | meta-tests: injected failures demonstrably turn the gates red |

## Runner honesty

- `macos-14` = Apple Silicon (aarch64, NEON): a true reference-family runner.
- The `scripts/ci/x86_cross_check.sh` lane (bead ebro) is the ARM-INVISIBLE
  BREAKAGE firewall: `cargo check --locked --workspace --all-targets --target
  x86_64-unknown-linux-gnu` from the aarch64 dev machines, followed by one
  derived `cargo check --locked` for every test/bench/bin/example target whose
  manifest declares `required-features`. This covers the default target surface and the
  gated targets that `--all-targets` alone skips. No cross-linker is needed
  (`check` does not link). Each invocation owns a collision-proof log directory
  and a separate Cargo target directory; concurrent proofs never share mutable
  build output.
  Because this is part of the required DSR gate, a missing Rust target is a
  named failure with the installation command, not a zero-exit skip.
  Every Cargo command is lockfile-closed. Its retained JSONL verdict binds the
  exact before/after root content, `constellation.lock`, and every clean pinned
  sibling tree; root or sibling movement fails the lane. The sole terminal row
  says whether that provenance is `sealed` or `incomplete`. It cannot prove
  runtime behavior or generated machine-code quality - real x86 execution still
  lives on the Threadripper lanes.
- The `scripts/ci/x86_runtime_sweep.sh` lane (bead yuyy) closes that runtime
  gap for the pre-release flow: on the first quiet, reachable Threadripper
  (ts1/ts2), fast-forwarded to origin/main, it first runs `cargo run --locked
  -p xtask -- check-identities`, which recomputes and byte-compares the
  committed `identity-schemas.json`, then runs the full workspace tests
  (`--no-fail-fast`). Both commands use the same hermetic target directory,
  their exact statuses are retained, and either failure fails the lane. The
  lane is load-gated with bounded retries (the workers double as rch lanes
  whose background traffic can kill long runs), and dirty/diverged clones are
  refused with a named host diagnostic. Local and remote pre/post dirt
  probes disable global excludes, independently enumerate untracked paths using
  only repository `.gitignore` rules, and reject hidden index flags. Exactly one terminal JSONL
  verdict represents the overall run. The lane never stashes, drops, or
  rewrites remote worktree changes; even `Cargo.lock` churn is a refusal.
  Before Cargo, the remote lane retains a
  `scripts/ci/checkout_constellation.sh --snapshot` identity; after Cargo it
  captures the identity again and requires equality. That binds exact root
  content, `constellation.lock`, and every pinned/clean path-dependency sibling.
  Both Cargo commands run with `--locked`, and a pass requires each exact exit
  status, an
  unchanged before/after HEAD and snapshot, clean post-test remote and local
  worktrees, at least one
  completed suite, no failed suite, and a successfully retained full log. The
  retained marker block includes the final porcelain diagnostic when a test or
  build script dirties the clone. The local caller must itself be clean,
  and the fast-forwarded remote HEAD must equal that admitted local HEAD; a
  remote `origin/main` run can never stand in for uncommitted local changes.
  ~20 min;
  named-skip when no quiet host answers - staleness of the last pass row is
  the alarm to watch. First catch of this class: the obq0 fs-cheb per-ISA
  schedule divergence.
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
| `check-identities` | persisted/replayed identity schemas with unclassified source fields, missing exact-field mutations or replay-version refusal guard targets, stale generated coverage, or schema/domain drift without a golden-coupling bump (`identity-schemas.json`) |
| `check-claims` | claim-state drift in the tracker mirror |

Each is also runnable alone (same names). Golden re-pins follow
docs/GOLDEN_POLICY.md: committed tree, BOTH build modes, plausible root
cause, coupling row updated in the same commit.

## Required quality lanes

`scripts/ci/quality_lanes.sh` is a REQUIRED DSR check. It derives every
test, bench, binary, and example target with `required-features` from locked
Cargo metadata and checks/runs each target with its declared features and
`--locked`. The standalone
`crates/fs-wasm` workspace has three required obligations:

1. Native tests against its tracked lock.
2. A locked `wasm32-unknown-unknown` Cargo check.
3. A `wasm-pack build --dev --target web -- --locked` browser build whose
   invocation must not modify the nested lock.

The manifest, `wasm-pack`, and the Rust wasm target are DSR host prerequisites.
Absence is a named FAILURE, not an advisory skip. Every admitted invocation gets
a collision-proof run directory under `target/quality-lanes/` containing full
per-lane logs, an isolated wasm-pack output, `verdicts.jsonl`, and
`inventory.json`; it also owns a separate Cargo target directory. Override paths
are normalized to absolute paths before any nested-workspace command runs. The
inventory SHA-256 hashes canonical input file CONTENTS (manifests, locks, contracts, Rust
sources/tests, CI scripts, and the derived gated-target list), not merely
filenames. The same composite
root + clean-constellation snapshot used by the x86 lane must match before and
after the run.
Per-lane rows and the summary are explicitly `provisional` and carry JSON `null`
for the not-yet-known closing snapshot. The literal final JSONL row is
`quality-proof-seal`; it binds the SHA-256 of every preceding row to the final
snapshot and reports `sealed` or `incomplete`. Consumers MUST verify that seal
rather than treating an earlier row as an independent receipt. Once the run
directory exists, an EXIT finalizer also writes an `incomplete`/`fail` seal for
initial-admission errors and later implicit shell failures; a provisional row
cannot become the terminal record merely because `set -e` fired. SIGKILL,
host/process loss, or a storage failure can still prevent the finalizer from
writing, so a missing or malformed seal is always a failed proof, never a skip.

## Constellation

The workspace path-depends on sibling repos. The canonical fresh-checkout
path is `cargo run --manifest-path tools/bootstrap/Cargo.toml` (beads
huq.17 and 1t8i; docs/BOOTSTRAP.md). This standalone, zero-dependency package
builds before the root workspace resolves, fetches every repo from
`constellation.lock`'s recorded remotes at the pinned revisions, applies the
same pinned-head and clean-tree verification to existing and newly cloned
siblings, and refuses a noncanonical/oversized lock, duplicate or unknown
library, mismatched lock hash, path-unsafe library identity, or missing, wrong,
duplicate, unknown, or type-invalid top-level identity metadata before deriving
or accessing any destination. The canonical v2 lock identity is
`org.frankensim.xtask.constellation-lock.v1` at integer version `1`; its existing
row-only hash still covers `(lib, version, git_head)`. Clean verification
disables global ignores, catches files hidden by repository-local excludes, and
rejects `assume-unchanged` or `skip-worktree` index entries. It supports
`--offline` cache verification and `--from` mirrors, and writes v2 fetch
provenance with identity domain
`org.frankensim.xtask.constellation-bootstrap-provenance.v1` at version `1`; the
standalone and in-workspace producers share the same top-level and row shape,
distinguish the canonical upstream remote from the selected mirror/transport,
and record whether that transport was actually used. `cargo run -p xtask --
bootstrap-constellation` remains the in-workspace command after the sibling
paths already resolve.
`scripts/ci/checkout_constellation.sh` remains the shell-only
equivalent used by the manual workflow specs. It admits the same exact lock
schema/domain/version and refuses malformed identity metadata before any sibling
operation; retained synthetic coverage runs checkout, verify-only, and snapshot
through that parser. Both it and `xtask
check-constellation` now require every sibling to be at the pinned head with a
clean tracked/untracked working tree. The script's `--snapshot` mode is the
canonical repo-local CI content identity. Snapshot v2 length-frames HEAD, index
state, sorted paths, Git-semantic modes, regular-file SHA-256s, symlink target
bytes, explicit missing entries, the exact lock bytes, and every clean pinned
sibling head/tree. It never hashes rendered `git diff` text and fails closed on
unsupported special files. Cleanliness forces executable-bit and full untracked
reporting, disables global excludes, then independently enumerates untracked
paths using only the repository's `.gitignore` rules; ignored build artifacts
are intentionally outside source identity. Stable dirty FrankenSim roots are
therefore hashed exactly, while any dirty sibling is refused.

The standalone Rust, shell, and in-workspace xtask bootstrap paths initialize
new siblings in place with a local incomplete-state marker before fetching. A retry may resume
a clean marked checkout, or a clean unborn checkout with the exact expected
origin, without deleting anything. Ordinary dirty repositories and ordinary
repositories at a mismatched HEAD are never repaired automatically. Bumping a
sibling remains deliberate: re-run `cargo run -p xtask -- lock-constellation`
locally and commit the new lock (schema v2 records remotes).

## External DSR wrapper limitations

The repo-local scripts above fail closed and retain their own complete evidence,
but the current DSR wrapper/configuration lives outside this repository and is
not yet an evidence-complete aggregate. As of this revision, external DSR:

- can return success when its repository configuration is missing, malformed,
  unreadable, or resolves to zero checks;
- captures command output in memory and retains only the first 1,000 bytes in
  its result instead of storing a full per-command log;
- does not bind all configured checks to one before/after root + constellation
  snapshot; and
- reports dry-run commands as passed even though they were not executed.

Its separately configured top-level Cargo checks/build also need `--locked`.
Until those external issues are fixed, a bare "DSR 8/8" line is not sufficient
evidence by itself. Cite the exact repo-local `verdicts.jsonl`, full-log paths,
before/after snapshot, HEAD, and lock identities emitted by these scripts. DSR
remains the required orchestrator; this is an explicit no-claim boundary on its
current aggregate receipt, not permission to substitute GitHub Actions.

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
| P10 agent-first | JSONL verdicts from every suite retained as artifacts; structured-error schema validation | **partial** — repo-local quality/x86 verdict and full-log retention live; external DSR aggregation limitations are named above; catalog no-drift gate lands with fs-ir |

## Gate meta-tests (are the gates real?)

`ci-self-test.yml` proves the wiring by injection:

1. `FS_CI_INJECT_FAILURE=1` arms a deliberately failing xtask test → the test
   gate must go red (asserted by inverting the step's exit code).
2. A seeded formatting violation → `cargo fmt --check` must go red.
3. The xtask unit battery seeds layer/dependency/contract/unsafe violations
   and asserts each is caught (`cargo test -p xtask -- seeded`).

A green `ci-self-test` run is the evidence that a red gate blocks; run it
after any change to the workflows or to xtask's checks.
