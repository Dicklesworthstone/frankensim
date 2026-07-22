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

## Determinism tiers and the libm doctrine (bead frankensim-lyms)

A crate's `CONTRACT.md` determinism class is the claim surface, and the
claim decides the code rules:

- **Cross-ISA bitwise** ("fully deterministic", "identical reports
  bitwise", golden-bearing lower-layer crates): every transcendental
  (`ln`, `exp`, `sin`, `cos`, `tan`, `atan2`, `hypot`, `powf`, `cbrt`,
  inverse/hyperbolic forms) must route through `fs_math::det`, which is
  built from correctly-rounded IEEE-754 ops and is bit-identical across
  ISAs by construction. Platform libm is NOT correctly rounded and
  differs by ≥1 ULP across ISAs and libm versions. `sqrt` and exact ops
  (`abs`, `floor`, `rem_euclid`, `mul_add`, `to_degrees`) are exempt —
  IEEE-754 requires correct rounding for them. Crates at this tier are
  registered in `LIBM_DOCTRINE_CRATES` (xtask `check-libm`, part of
  `check-all`); dev-only oracle comparisons escape with
  `// det-ok: <reason>` on the same or preceding line.
- **Same-ISA bitwise**: raw libm is permitted (one build of one libm is
  self-consistent). The CONTRACT must say "same-ISA" (or otherwise scope
  the claim) rather than claiming unqualified full determinism.
- **Statistical / fast-mode**: out of scope for this doctrine.

Migrating a crate from same-ISA to cross-ISA shifts last-ULP outputs:
re-check every golden in the crate under the golden-bump protocol
(`docs/GOLDEN_POLICY.md`) in the same change, then add the crate to
`LIBM_DOCTRINE_CRATES` so the doctrine is enforced, not documented.

## Claim-integrity defect class (bead f85xj.2.1)

A **claim-integrity defect** exists when any public surface — API return type,
evidence color, certificate, report line, CONTRACT/README sentence, WASM/CLI
export, ledger row — can assert a **stronger** epistemic state than its actual
evidence establishes. A false certificate is worse than an ordinary wrong
answer, so these are tracked as their own countable, gateable class rather than
scattered among ordinary bugs. The normative definition, decision rules, audit
method, and known-answer set live in [`docs/CLAIM_INTEGRITY.md`](CLAIM_INTEGRITY.md);
`check-claims` (part of `check-all`) lints that it stays present and intact.

Label taxonomy — every such bead carries all three:

| Label | Meaning |
| --- | --- |
| `claim-integrity` | Mandatory class membership. `br list -l claim-integrity` is the live inventory and the gate reads exactly this label. |
| exactly one `severity:*` | `severity:default-path` (P0, reachable on a default/public path), `severity:gated` (P1, non-default feature or opt-in API only), `severity:doc-only` (P2, prose overstates honest code). |
| `crate:<name>` | Every crate whose surface can emit the claim. No crate scope means **global** (fail closed, blocks every promotion), never unscoped-and-harmless. |

Defects are filed as `--type=bug`; the type is what separates *exposure* from
*work*. A `claim-integrity`-labelled `epic`/`task`/`feature` is an E02 program
bead (this doctrine, its sweep, its gate) and is exempt from the severity and
ownership rules — otherwise the gate would count its own epic as an open P0 and
block every promotion forever.

Severity is decided by reachability, escalates rather than averages, and
resolves ambiguity toward the stronger severity. Under-claiming is never a
claim-integrity defect. Bead bodies must carry a minimal repro and the honest
claim the surface should make instead. Run
`scripts/ci/claim_integrity_inventory.sh` for the checked live inventory; it
fails when an open P0 has no owner, when severity labels are missing or
duplicated, or when the beads store cannot be read.

## EXTREAL critical-path triage (bead f85xj.16.2)

The EXTREAL vertical is the program critical path. `.beads/issues.jsonl`
remains authoritative for issue state and dependency edges;
`vertical-capability-graph.json` is a checked projection that binds every
`capability-maturity.json` capability to implementing Beads and names the
owners of the four integration seams. It intentionally does not copy Beads
edges or assert that membership, issue closure, or crate-local green proves a
capability mature.

Every Bead whose id starts with `frankensim-extreal-program-f85xj` carries the
`extreal` label. Use robot modes only:

```bash
bv --robot-plan --label extreal
bv --robot-triage --label extreal
bv --robot-alerts --label extreal
```

The plan result must expose a summary, parallel tracks, status, and a data
hash. The triage result must expose quick reference counts, recommendations,
blockers, and project health. If a local `bv` version rejects `--label` for a
particular robot command, use `bv --recipe actionable --robot-plan` to inspect
the global actionable graph and retain the limitation in the triage record;
never fall back to bare interactive `bv`.

The named seam owners are integration owners, not crate owners:

| Seam | Owner | Accountability |
| --- | --- | --- |
| `schema-scenario` | `jemanuel` | Versioned project inputs preserve persistent scenario identity and the Five Explicits. |
| `physics-evidence` | `jemanuel` | Physics QoIs retain numerical, model-form, experimental, and regime authority. |
| `corpus-scorecard` | `jemanuel` | Corpus identity, domain, uncertainty, and limitations reach the public scorecard. |
| `cli-session` | `jemanuel` | CLI work enters governed sessions with budgets, cancellation, replay, and structured output. |

Crate-local tests do not establish seam health. The owner must demand an e2e
exercise across the seam before reporting it integrated. Work outside the
critical path needs the `f85xj.16.3` displacement-rule justification: either a
named trust-risk reduction or capped research with a falsifier, budget, and no
claim on critical-path capacity.

### Weekly recipe

1. Run `bv --robot-plan --label extreal`; inspect blockers and whether each
   track advances the Cooling 0.1 vertical or a linked `ty23`/`v6dn` gate.
2. Run `bv --robot-triage --label extreal`; challenge the top recommendations
   against seam risk, maturity gaps, and the actual dependency chain.
3. Check `cargo run -p xtask -- check-critical-path`; repair missing labels,
   stale capability bindings, missing owners, or dangling issue ids before
   selecting work.
4. For off-path work, record the `f85xj.16.3` justification before assigning
   capacity. Record the command, UTC timestamp, data hash, result shape, and
   counts from the weekly run; do not silently replace prior receipts.

The first retained exercise ran `bv --robot-plan --label extreal` at
`2026-07-22T21:19:43Z`. It returned the documented summary/tracks/status/hash
shape with data hash `5ba7d5d6e7f40830`, 178 scoped issues, 68 actionable, and
61 blocked. The complete deterministic receipt fields are retained in
`vertical-capability-graph.json`; those counts are historical observations,
not promises about the live graph.

## Moonshot portfolio and displacement rule (bead f85xj.16.3)

`moonshot-portfolio.json` is the checked declaration layer for `[M]` work;
Beads remains authoritative for labels, status, and dependency edges. A Bead
carrying `ambition:M`, `ambition-m`, or a `moonshot*` label consumes WIP only
while its status is `in_progress`. Open entries are unadmitted backlog, not
quietly active research. Closed entries are terminal history.

`fs_govern::MoonshotPortfolio` is the typed policy core. It projects every
declaration into the existing proof-lane ledger, exposes only consumed-value
`liquidate` and one-for-one replacement transitions, retains exact terminal
dispositions, and prevents terminal Bead revival. The xtask command is the
explicit live-world adapter around that pure state machine.

The initial live sweep found 508 tagged Beads: 6 active, 482 open backlog, and
20 closed. The immutable v1 baseline names those six active Beads, so editing
the current snapshot cannot manufacture a seventh slot or erase history. The
status-quo cap is therefore six. Completion or shelving without a replacement
adds a `liquidate` history row and reduces both active WIP and the declared cap
in the same update; it cannot leave a refillable slot. Every later admission
adds a `replace` row that displaces one named active lane and records one of
three terminal dispositions:

- `completed`, with its retained state/evidence locator;
- `falsified`, with the counterevidence locator;
- `shelved-with-state`, with enough state to resume only through a future
  displacement decision.

Every active declaration names its owner, falsifier observation and decision
rule, effort cap, calendar deadline, next quarterly falsifier review, evidence
required at that review, and an explicit no-critical-path-capacity statement.
The latter is not trusted as prose: `cargo run -p xtask -- check-moonshots`
reads blocking edges from `.beads/issues.jsonl`, derives the product-path set
from `vertical-capability-graph.json`, and rejects direct or transitive
moonshot reachability into that set. `check-all` includes the same audit.

At each quarterly falsifier review, evaluate the named observation against
retained evidence. A hit closes or shelves the mechanism under its registered
decision rule; an absent or inconclusive test does not count as survival and
does not widen the budget. Changes to the active inventory must update the
checked snapshot, declaration, and disposition history in the same change as
the Beads status. A terminal Bead cannot be replayed as a later replacement.

The v1 file retains one hypothetical-only admission rehearsal. It takes the
real open rank-deficient TSQR theorem Bead from the unadmitted backlog, names
the active spectral-bridge lane it would displace, preserves that lane as
`shelved-with-state`, and records the seven admission decisions. The expected
result is an admissible replacement at six-of-six active slots. The rehearsal
changes no Bead state and grants no scientific or product authority; it
exercises the policy shape and live dependency graph only.

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
