# Abstraction Consolidation Review

Bead: `frankensim-extreal-program-f85xj.16.8`. Record:
[`consolidation-review.json`](../consolidation-review.json). Gate:
`cargo run -p xtask -- check-consolidation` (also runs inside `check-all`).
Cadence: once per release train (the e13.4 anchor).

The record is not a hand-maintained list. `check-consolidation` re-derives the
usage sweep from the manifests on every run, so a disposition that the tree has
outgrown fails the gate instead of quietly becoming fiction.

Across the current 143 `fs-*` crate directories (142 native root-workspace
members plus the standalone nested `fs-wasm` workspace), an abstraction that no
supported workflow exercises is not neutral. It costs comprehension,
maintenance, and CI time, and it dilutes the signal of what the product actually
is. The moonshot WIP cap
(`docs/CONVENTIONS.md`, bead 16.3) governs **new** speculative work; this review
governs the **existing** inventory's slow accretion.

It is the retirement-side counterpart to the schema freeze
([`docs/SCHEMA_POLICY.md`](SCHEMA_POLICY.md)): that policy decides what is
promised, this one decides what stays.

## The usage sweep

A **supported workflow** is a vertical, campaign, or e2e lane — concretely, every
`fs-*-e2e` crate plus the named verticals and campaigns and the product CLI. The
sweep takes the transitive closure of `fs-*` dependency edges from every
workflow root, counting **runtime and dev-dependency edges**, and reports the
crates outside it.

Dev edges are counted deliberately: a crate used only as a test oracle by a
workflow *is* exercised by that workflow, and excluding dev edges would
overstate the candidate set.

The root list is **part of the reviewed record**, not an implementation detail.
Adding a workflow root changes the candidate set, so a root-list change is a
reviewable event.

### What the sweep proves, and what it does not

It proves the absence of a dependency **edge** from a workflow root. It does not
prove absence of use. A crate exercised only through a doc example, a binary, or
a harness outside the crate graph would still be flagged. Every flagged crate
therefore gets a human-reviewed disposition; the sweep ranks attention, it does
not decide.

### Validation

Each run records two checks, so a silently-broken sweep cannot masquerade as a
clean inventory:

- a **known-exercised control** must be reached and unflagged (currently
  `fs-sparse`, 17 dependent crates);
- **zero-dependent spot checks** on flagged crates, confirming by independent
  grep that nothing depends on them.

## Dispositions

Every candidate gets exactly one, recorded with its rationale.

| Disposition | Meaning |
|---|---|
| **KEEP** | A named consumer, or a named trust-risk rationale that justifies the crate independently of current consumption. The default. |
| **CONSOLIDATE** | Merge into a neighbour, with a migration note. |
| **FREEZE** | Explicitly parked and visible: compiles, tests green, contract kept, **no new investment**. A visible state, not a quiet death. |
| **REPAIR-OR-QUARANTINE** | Parked but **not** green. See below. |
| **RETIRE** | Removal proposal only. Agents **propose**; they never execute. Requires the owner's explicit approval under the repository's no-deletion rule. |

### FREEZE has a green precondition

FREEZE asserts "compiles + tests stay green but no new investment". A candidate
whose tests are red therefore **cannot** be frozen — labelling it FREEZE would
assert a green parked state that does not exist, which is exactly the kind of
flattering-label drift the claim-integrity work exists to prevent.

`REPAIR-OR-QUARANTINE` is the honest disposition for that case: parked, known
red, and named as such. This vocabulary gap surfaced in the first review, from
`fs-wasm`.

### Consolidation must not buy false economy

A CONSOLIDATE may never break the layer rules or merge crates across determinism
classes. `check-layers`, `check-contracts`, and `check-deps` gate the
consolidation itself — the post-consolidation gate run is the acceptance
evidence, not the reviewer's judgement.

## Relationship to the maturity registry

The sweep is not only a cleanup list. It is an **L3-promotion falsifier**.

L3 means "integrated workflow". A capability whose crates no supported workflow
reaches cannot honestly hold an L3 claim, whatever its test quality. The first
review demonstrated this working in both directions:

- `geometry.topology-certificates` is registered **L2** with the boundary note
  "Not L3 — no e2e lane exercises it". The sweep reached that conclusion
  independently, from the dependency graph alone.
- `wasm.browser-flagships` (**L1**) is likewise unreached, consistent with its
  recorded build break.

Every review therefore cross-references the registry and records, per
capability, whether a workflow consumer exists.

## Running a review

1. Run the usage sweep from the workflow roots; record totals and both
   validation checks.
2. Cross-reference `capability-maturity.json`: every capability whose crates are
   all unreached gets an explicit disposition.
3. Give every candidate a disposition plus rationale, with mechanical evidence
   (test-file count, source-file count, last commit, dependent count).
4. Execute any FREEZE/CONSOLIDATE: make the state **visible** where a developer
   will meet it — a notice in the crate's `CONTRACT.md` — and attach the green
   gate as acceptance evidence.
5. Append the review to `consolidation-review.json` so the next review starts
   from state.

## First review: CONS-2026-07-A (2026-07-24)

141 crates, 23 workflow roots. 115 reached through runtime edges, 121 including
dev edges, leaving **20** exercised by no supported workflow.

Dispositions: **18 KEEP**, **1 FREEZE**, **1 REPAIR-OR-QUARANTINE**. No RETIRE
was proposed.

The headline finding is deliberately unglamorous: every candidate had been
committed within the preceding three weeks and all but two ship test files. This
inventory is **"not yet composed into a workflow", not "abandoned"** — which is
why KEEP is the default and why a review that mostly keeps things is a
successful review, not a wasted one. The accretion this cadence exists to catch
is a future condition; the first run establishes the baseline that makes it
detectable.

`fs-dimine` was frozen as the worked example: no workflow consumer and no
dependent crate anywhere, the oldest candidate, a single dependency, and 9/9
conformance tests verified green at review time. `fs-wasm` produced the
vocabulary gap that added `REPAIR-OR-QUARANTINE`.
