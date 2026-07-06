# FrankenSim Workspace Conventions

The mechanical rails for a workspace built by a swarm of agents (plan §16.2).
The Decalogue is the tie-breaker: disputes resolve by principle number, not
seniority. This file states conventions; `cargo run -p xtask -- check-all`
enforces the checkable subset.

## One crate = one contract

Every `crates/fs-*` crate ships a `CONTRACT.md` (template:
`docs/CONTRACT_TEMPLATE.md`) with all ten required sections, plus a
`tests/conformance.rs` suite. Agents integrate against contracts and the
FrankenScript IR — never against another crate's internals. "Done" is defined
as conformance green on both target ISAs, nothing else (plan §13.3).

## Layer discipline (enforced: `check-layers`)

Every crate declares its layer in `[package.metadata.frankensim]`:

| Layer | Meaning | May depend on |
| --- | --- | --- |
| `UTIL` | dependency-free utilities (fs-qty, fs-obs) | UTIL |
| `L0`   | SUBSTRATE | UTIL, L0 |
| `L1`   | BEDROCK | UTIL, L0, L1 |
| `L2`   | MORPH | UTIL, L0–L2 |
| `L3`   | FLUX | UTIL, L0–L3 |
| `L4`   | ASCENT | UTIL, L0–L4 |
| `L5`   | LUMEN | UTIL, L0–L3, L5 (**not** L4 — siblings) |
| `L6`   | HELM | everything |
| `TOOL` | repo tooling (xtask) | anything; nothing depends on it |

Lower layers never know about higher ones. If you need to violate this, stop
and redesign the API (AGENTS.md).

## Dependency policy (enforced: `check-deps`)

Runtime dependencies are exactly: `std`, workspace `fs-*` crates, and the
Franken constellation (asupersync, FrankenSQLite, FrankenNumpy, FrankenTorch,
FrankenScipy, FrankenPandas, FrankenNetworkx). No BLAS/LAPACK/C/C++/FFI, no
general crates.io dependencies in production paths (Decalogue P1).
Dev-dependencies are exempt but reported: dev-only comparison oracles must be
isolated and documented in the crate's CONTRACT.md.

## Unsafe policy

`unsafe_code` is denied workspace-wide via `[workspace.lints]`. Unsafe is
permitted only in registered capsule modules: `<300` lines, a safe façade, a
`SAFETY.md` per the unsafe-capsule bead (invariants, aliasing, alignment,
lifetimes, panic/cancellation/concurrency behavior, Miri/model-check/fuzz
coverage, caller obligations). The capsule registry is the source of truth;
unregistered `#[allow(unsafe_code)]` or `unsafe` fails CI via `cargo run -p xtask -- check-unsafe` (part of check-all).

Panic policy (decided, per plan §5.2): panics are captured inside tile scopes
and converted to structured diagnostics before crossing layers — never
`panic = abort` in shipped kernels, never a process abort mid-campaign.

## Feature flags

- `[S]` solid work: no flag, default path.
- `[F]` frontier: `frontier-<name>`, default **off**.
- `[M]` moonshot: `moonshot-<name>`, default **off**, promoted only after the
  Gauntlet's certifier trials (§13.2). Nothing `[M]` gates anything `[S]`.

## Toolchain

Pinned nightly (`rust-toolchain.toml`), edition 2024. Nightly features are
permitted only in narrow, documented places (currently: const-generic dimension
arithmetic in `fs-qty`; `experimental-portable-simd` tier in `fs-simd`).
Everything else must compile on stable in principle — treat each nightly
feature as a liability with a documented fallback.

## Manifest conventions (what the xtask parser relies on)

Workspace manifests are generated/edited to a known shape: section headers on
their own line, one dependency per line (`name = { path = ".." }` inline
tables on a single line). The xtask TOML parser handles exactly this subset
and fails loudly otherwise — keep manifests boring.

## Output style

Core library code never prints casually. Observability is structured events
(fs-obs schema) or ledger records. Test suites emit JSON-lines verdicts with
seeds, fixture hashes, and case ids — every failure must be reproducible from
its log line alone.

## Compiler checks

After substantive changes:

```bash
cargo fmt --check
cargo clippy --all-targets -- -D warnings
cargo test --workspace
cargo run -p xtask -- check-all
```

Use RCH lanes for compute-heavy builds in shared-agent environments
(AGENTS.md).
