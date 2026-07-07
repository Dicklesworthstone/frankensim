# fs-constraint — CONTRACT

## Purpose and layer

L4 (ASCENT). The constraint CALCULUS (plan §9.1, patch Rev F):
constraints with SEMANTICS, not anonymous `g(x) ≤ 0` — a typed kind
taxonomy with per-kind optimizer treatments, evidence-typed
evaluations, and infeasibility DIAGNOSIS (minimal unsat cores, ranked
calibrated repairs) that turns optimizer failures into design
conversations. fs-opt hosts the expression graphs; this crate owns
what the constraints MEAN.

## Public types and semantics

- `ConstraintKind`: `Hard` (never traded → feasibility restoration),
  `Soft(PenaltyLaw::{Quadratic, Hinge})` (→ penalty term), `Chance
  {level, ChanceEstimator::MonteCarlo{samples, delta}}` (→ estimate
  then act on the BOUND), `Robust{half_widths}` and `Certification
  {ProofKind::{Interval, Sos}}` (→ prove or escalate), `Fabrication
  {process}` and `Code{standard}` (→ domain check; semantics named
  for the ledger). `treatment()` is the routing table.
- `evaluate(problem, spec, x, noise) -> ConstraintEvidence`: status
  (`Satisfied`/`Active`/`Violated`/`NeedsProof`/`Proven`/
  `BoundNotCleared`), EXACT violation certificates for algebraic
  graphs, active-set role, penalty per law. Chance kinds compute a
  Hoeffding lower confidence bound and report satisfied ONLY when the
  bound clears the level — `BoundNotCleared` exists precisely for the
  case where the raw empirical rate clears but the bound does not
  (the validity machinery is the feature). Certification kinds refuse
  "satisfied" pointwise REGARDLESS of how good `g(x)` looks.
- `interval_eval` (the in-house prover): rigorous inclusion over
  fs-opt graphs per node; refuses (teaching reasons) on division
  through zero, domain violations, negative powers, PDE/stochastic
  nodes. `prove_interval` turns a provable domain into a `Proven`
  status + `ProofArtifact::IntervalBound`. Robust kinds are proven
  conservatively over their uncertainty boxes the same way, carrying
  ENCLOSURE certificates.
- `diagnose_infeasibility(problem, specs, domain, cx) -> Diagnosis`:
  elastic-relaxation solve (multi-start projected subgradient descent
  on total hinge violation, deterministic LCG starts) classifies
  feasibility and yields a witness or an unsat core seeded from the
  elastic support and refined by the DELETION FILTER — the core is
  MINIMAL (dropping any member restores feasibility). Repairs
  (relax-bound at graded slacks; drop-soft for soft members) come
  RANKED by Monte-Carlo feasible-volume estimates.
- `serialize_specs`/`parse_specs`: canonical line form (floats as bit
  patterns), identical round-trips, line-numbered refusals.
- `ConstraintEvidence::to_ledger_row` and `Diagnosis::to_json`: the
  Rev S ledger row and the agent-facing diagnosis payload.

## Invariants

1. All seven kinds map to their declared treatments; spec sets
   round-trip identically; ledger rows validate through fs-obs
   (fscon-001).
2. Statuses, roles, exact violation certificates, and both penalty
   laws evaluate as declared (fscon-002).
3. The chance BOUND decides: on an analytic uniform-noise fixture the
   raw rate clearing the level while the Hoeffding bound does not
   yields `BoundNotCleared`, never `Satisfied`; the half-width
   travels as a `StatisticalCertificate` (fscon-003).
4. Certification kinds refuse without artifacts; interval proofs
   succeed exactly on provable domains and refuse honestly otherwise;
   interval containment holds over random nonlinear boxes (G0);
   robust kinds carry enclosures (fscon-004).
5. Unsat cores are MINIMAL against enumeration: the core is jointly
   infeasible, every single deletion restores feasibility, bystanders
   are excluded, feasible systems return witnesses, and elastic
   verdicts match grid enumeration on the seeded fixture family
   (fscon-005).
6. Repairs are ranked by feasibility estimate, soft members offer
   drop actions, and estimates are CALIBRATED against exact
   enumeration (worst gap < 0.05 on the worked mass/strength
   example); the diagnosis payload ships through fs-obs (fscon-006).

## Error model

`ConError` teaching errors: `NotScalar`, `Eval` (fs-opt errors carried
through), `NotProvable{why}` (an honest gap with escalation advice,
not a failure), `BadParam`, `Parse{line, what}`. The interval engine's
`IvalError` names each refusal reason.

## Determinism class

Fully deterministic: LCG-seeded multi-starts and Monte-Carlo streams,
canonical constraint ordering, bitwise float serialization. Identical
inputs give identical diagnoses, estimates, and bytes.

## Cancellation behavior

`elastic_solve` (and therefore `diagnose_infeasibility`) polls
`cx.checkpoint()` per restart and returns the carried `Cancelled`
teaching error between solves.

## Unsafe boundary

None. `#![forbid(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs`, cases fscon-001..fscon-006 — JSON-line
verdicts, seeded LCG randomness, fs-obs events for ledger rows and
the full diagnosis payload. Any reimplementation must pass the suite
unchanged.

## No-claim boundaries

- The interval prover rounds to-nearest; outward-rounded arithmetic
  joins with fs-ivl (containment carries an fp-slack caveat until
  then). SOS certificates are REPRESENTED (`ProofKind::Sos`), not
  executable — fs-sos is a later bead.
- The elastic solver is small-fixture machinery (multi-start FD
  subgradient); the production feasibility-restoration solver is a
  later ASCENT bead. Nonconvex fixtures can defeat it — verdicts are
  cross-checked against enumeration only at conformance scale.
- Chance estimation is Monte-Carlo/Hoeffding v1; e-process anytime
  validity and richer estimators join with the UQ beads.
- Fabrication/Code kinds carry semantics to the ledger; process
  models (fs-fab) and code-check rule packs bind in their beads.
- Repair generation covers bound relaxations and soft drops; material
  /topology switches (the patch's richer vocabulary) need fs-xform
  and fs-fab integration.
- Host problems are single-Rn-variable v1; multi-variable and
  manifold-variable domains generalize with the restoration solver.
