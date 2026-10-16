# Capability Maturity Levels

Bead `frankensim-extreal-program-f85xj.16.1`. This file defines the levels;
`capability-maturity.json` records where each capability actually sits; and
`cargo run -p xtask -- check-maturity` enforces the mechanical parts.

"Is this capability experimental, verified, integrated, validated, or
supported?" previously had no queryable answer — status lived in README prose,
CONTRACT no-claim sections, and folklore. This registry makes it a fact you can
query, and makes **promotion an event** rather than a mood.

## What a capability is

A capability is a **product-meaningful unit of work a user could ask for** —
"steady conduction solve", "evidence packaging", "deterministic sparse
assembly". It is deliberately *not* a crate:

- one crate can host several capabilities (`fs-evidence` hosts both colour
  algebra and model-card validity domains);
- one capability can span several crates (evidence packaging spans
  `fs-package`, `fs-checker`, and `fs-crosswalk`).

Crates are an implementation detail; capabilities are what maturity is claimed
about. The `crates` field on each entry is scope metadata for the claim-integrity
gate, not the identity of the capability.

## The levels

Each level requires **everything below it**, plus its own bar. A capability sits
at the highest level whose criteria its recorded evidence actually meets — not
the level someone hopes it reaches.

### L1 — Experimental component

The code exists and is honest about itself.

- Implemented and compiling in the workspace.
- `CONTRACT.md` present for every hosting crate, with its no-claim boundaries
  filled in rather than boilerplate.
- No claim of numerical correctness beyond what the contract states.

L1 is not an insult. Most of this workspace is L1, and saying so is the point.

### L2 — Numerically verified

The mathematics is checked against something independent of itself.

- Named Gauntlet tiers green for the capability (typically G0 algebraic laws,
  plus G1 manufactured solutions or convergence order where a PDE is involved).
- Evidence refs resolve to real, currently-passing tests.
- An independent oracle, closed form, symmetry law, or exact/interval check —
  *not only* a golden hash. A golden alone freezes behavior; it does not verify
  meaning, so a capability whose only evidence is a golden stays L1.

### L3 — Integrated workflow

The capability composes with its neighbours end to end.

- Exercised by an e2e lane, not only unit tests.
- Its inputs and outputs cross at least one crate boundary under the real
  types (no fixtures standing in for the adjacent layer).
- Failure modes are typed and the lane covers at least one of them.

### L4 — Experimentally validated

The model has been compared against reality over a stated domain.

- Validation against an external corpus (`fs-vvreg`-registered datasets), with
  the validity domain recorded.
- Coverage stated: which regions of the domain the corpus actually covers, and
  which it does not.
- Model-form discrepancy quantified and carried as evidence, not assumed zero.

L4 is the first level at which the word "validated" may appear in user-facing
prose about the capability.

### L5 — Supported product

Someone has promised to keep it working.

- A written support policy: what is guaranteed, for how long, and to whom.
- Migration promises for breaking changes.
- A named owner accountable for regressions.

Nothing in this repository is L5 today, and nothing should be marked L5 until
the support policy exists as a document rather than an intention.

## Promotion is an event

Level increases are **promotions** and are treated as governed events:

- `check-maturity` compares the working registry against the last committed
  one and reports every promotion it finds.
- A promotion must be justified by evidence refs that resolve — the check
  refuses a promotion whose new level cites evidence that does not exist.
- Bead `.2.3` (the claim-integrity promotion gate) consumes exactly this
  promotion signal: an open `severity:default-path` claim-integrity defect whose
  crate scope overlaps a promoted capability blocks the promotion.

**Demotions are always allowed** and are never blocked. Lowering a claim is
how the registry stays honest, so the check logs a demotion as a policy note
and passes. A system that makes it procedurally harder to weaken a claim than
to strengthen one will accumulate false claims by construction.

## Staleness

Every entry carries `last_review`. A level is a claim about the present, and an
unreviewed claim decays: `check-maturity` reports entries whose review is older
than the staleness threshold. Staleness is reported, never fatal — an old review
is a prompt to look, not evidence that the capability broke.

## Evidence refs

Each entry records the evidence justifying its level. Refs are typed so the
check can resolve them mechanically:

| kind | meaning | resolved by |
| --- | --- | --- |
| `test` | `<path>::<fn name>` | the file exists and contains `fn <name>` |
| `contract` | a crate `CONTRACT.md` | the file exists |
| `lane` | a `scripts/ci/*.sh` lane | the script exists |
| `corpus` | an external V&V dataset id | recorded; **not** mechanically resolvable yet |
| `doc` | a repo document | the file exists |

A ref that does not resolve is a violation: the registry may not cite evidence
that is not there. `corpus` refs are recorded but not resolved, because the
V&V corpus registry (`e04`) does not exist yet — that gap is the honest reason
no capability here is L4.

## Capability-to-work projection

`vertical-capability-graph.json` cross-references every capability id here to
the Beads that implement or harden it. The projection is checked by
`cargo run -p xtask -- check-critical-path`: every referenced issue must exist,
every registry capability must appear exactly once, and any L2-or-higher row
must have at least one closed implementing Bead. This is a consistency floor,
not promotion evidence. The evidence refs above remain the authority for the
maturity level, while `.beads/issues.jsonl` remains the authority for issue
status and dependency edges.

The same projection names owners for the four EXTREAL integration seams and
retains the shape and identity of a successful robot-plan exercise. See
`docs/CONVENTIONS.md` for the weekly triage recipe and the explicit no-claim
boundary between a closed task, a green seam, and a mature capability.
