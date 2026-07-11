# CONTRACT: fs-evidence

## Purpose and layer
`Evidence<T>` / `Certified<T>` (patch Rev B): numerical, statistical, and
MODEL-FORM certificates carried inside values, with a conservative
composition algebra, model cards + the registration lint, two-fidelity
discrepancy models with out-of-distribution refusal, model bracketing, and
decision-aware escalation. The reason this crate exists: without model
evidence the system can produce beautifully certified WRONG answers (mesh
error 0.7%, closure discrepancy 10%) — being able to SAY that is the
product. Layer: UTIL (usable by every layer; the bead label said L6, but
the bead scope explicitly demands a low-layer home — this is it). Depends
on fs-obs only.

## Public types and semantics
- `Evidence<T> { value, qoi, numerical, statistical, model, sensitivity,
  provenance, adjoint_ref }` — the traveling noun. `Certified<T>` is the
  opaque newtype whose constructor discipline is `Evidence::certified()`: rigorous
  numerics (Exact/Enclosure), valid statistical parameters, and in-domain
  non-negative model discrepancy; pure math
  certifies with `ModelEvidence::none()` (the explicit "no model involved"
  statement); refusals are structured `CertifyError`s.
- `NumericalCertificate { kind: Exact|Enclosure|Estimate|NoClaim, lo, hi }`
  — the plan Appendix B `Certified` fields (value + interval bound +
  provenance + adjoint hook) kept intact as the numerical slice. Severity
  is ordered; float composition never claims Exact; NoClaim absorbs.
- `StatisticalCertificate { None | EValue{e, alpha} | HalfWidth{...} }` —
  finite non-negative e-values and widths with levels/confidences strictly in
  `(0,1)`; v1 composition is conservative-weakest (see no-claims).
- `ModelEvidence { cards, assumptions, validity, discrepancy_rel,
  in_domain }`; `ValidityDomain` — named-parameter boxes with
  intersection/containment; `SensitivitySummary` — d(qoi)/d(param)
  headlines, merged by magnitude.
- `Evidence::combine(op, a, b, value)` — Add/Sub/Mul/Min/Max on the QoI
  with certificates composed conservatively and provenance chained
  (`ProvenanceHash::chain`, order-sensitive, FNV until the ledger hash).
- `Evidence::assess(threshold_rel) -> DecisionStatus` and
  `UncertaintyBreakdown` — per-source relative bands, first-order-sum
  total, dominant source with the declaration-order tie law (ModelForm
  first: ties escalate the model, the band cheap refinement cannot fix).
  `escalation_advice` maps dominance to RefineNumerics /
  GatherMoreSamples / EscalateModelFidelity (the HELM governor hook).
- `ModelCard` (name, version, ambition tag, assumptions, validity, known
  failures, calibration provenance, discrepancy band) and `ModelRegistry`
  — `register_solver` REFUSES without a registered card (the lint).
- `DiscrepancyModel::fit(&[FidelityPair])` — observed parameter box +
  mean/max relative discrepancy; `query`/`evidence_at` refuse
  out-of-distribution points with the violated parameter named
  (`OutOfDomain`).
- `ModelBracket` — N plausible models; evidence = midrange value, an
  enclosure spanning every member, spread as the model band, and a
  bracket-spread sensitivity entry (the vessel flagship's contact-line
  mitigation).
- `to_ledger_row_json` on evidence and cards — the `evidence` /
  `model_cards` table rows (canonical order, no clocks, no addresses).

- `color` module (bead qmao.1): the THREE-COLOR epistemic schema —
  `Color::{Verified{lo,hi}, Validated{regime: ValidityDomain, dataset},
  Estimated{estimator, dispersion}}` with the `ColorRank` lattice
  (verified > validated > estimated), the TOTAL conservative pairwise
  `compose` (result rank = min of operands; verified intervals combine
  per `IntervalOp` with outward-rounded arithmetic; validated regimes
  INTERSECT, and a disjoint intersection demotes with infinite/no-claim
  dispersion; estimated absorbs everything with additive dispersion),
  `check_regime` (validated is a REGIONAL property: exiting, failing to
  report a regime axis, supplying a non-finite state, or declaring an
  empty/non-finite/inverted regime AUTO-DEMOTES to estimated with a
  `Demotion` flag), `verified_from`
  (the only door to a verified color — non-enclosure certificates
  refuse with the laundering teaching error), and `color_of` (the
  honest bridge from existing Evidence receipts). `Color::payload_json`
  escapes caller-controlled strings and represents non-finite floats as
  tagged JSON strings, never invalid bare numeric tokens. The distinct
  `Color::canonical_bytes` identity encoding is versioned (v1), structurally
  length-prefixed, deterministically ordered, and preserves every IEEE-754 bit;
  display rounding therefore never aliases color identity or authorization.
  Write-time enforcement lives HELM-side in fs-ledger over these types.

- `falsify` module (bead qmao.4): FALSIFIER PAIRING — `FalsifierRegistry`
  (a certificate class CANNOT register without ≥1 independent falsifier;
  `standard()` ships the six proposal pairings: watertightness→ray-parity,
  conservation→independent-quadrature flux audit, adjoint→FD spot checks,
  surrogate→held-out points, symmetry-block→occasional full solves,
  validated-color→held-out anchors), `ship_gate` (the no-falsifier-no-ship
  Gauntlet gate), `FalsifierHistory` (per-class-per-regime pass/hit/compute
  rows; `doubt` = 1 − pass rate with COLD-START = max doubt and a
  never-zero floor), `record_hit` → mandatory `(Tombstone, EstimatorBug)`
  canonical-JSON pair, `allocate_budget` (consequence × doubt ×
  rent-share, normalized; consequence floors for dependent-free claims;
  zero claims spend zero), and `rent_review` (zero-yield classes at
  meaningful volume decay toward a floor — every falsifier pays rent, but
  the pairing rule itself is not killable).

## Invariants
1. Conservativeness (G0, evd-001): composed enclosures contain every
   propagation of operand-enclosed true values (300k seeded samples);
   composed validity is exactly the per-parameter intersection;
   assumptions union sorted; discrepancy bands add; `in_domain` is a
   conjunction. Indeterminate IEEE endpoint arithmetic widens to the whole
   real line instead of discarding NaN corners.
2. Severity monotonicity: composition kind = max operand severity, floored
   at Enclosure for float ops; NoClaim absorbs to infinite bounds; an
   Estimate anywhere poisons `certified()` downstream (evd-006).
3. No card, no solver: `ModelRegistry::register_solver` refuses unknown
   cards with teaching text (evd-002).
4. Out-of-distribution discrepancy queries refuse with the violated
   parameter named — never silent extrapolation; non-finite training or query
   coordinates are unusable, not a way to synthesize or enter a trained box
   (evd-004).
5. Dominance ties break in declaration order (ModelForm, Statistical,
   Numerical) — deterministic verdicts.
6. Ledger rows and provenance chains are deterministic (repeat-identical).
   Evidence and model-card rows stay valid JSON under hostile metadata and
   tagged non-finite no-claim values; evidence rows retain model assumptions
   and sensitivity headlines rather than dropping those semantic slices.
   Public set-like card, assumption, and known-failure vectors are sorted and
   deduplicated again at the durable rendering boundary, so caller mutation or
   insertion order cannot change row identity.
   Provenance chaining is order-sensitive.
7. Color regime checks fail closed: no empty, inverted, or non-finite regime
   and no non-finite current state can retain `Validated`; disjoint validated
   composition and regime exit both carry infinite/no-claim dispersion.
8. Certified trust boundary (gp3.2.1, evd-012): `Certified<T>` is an
   OPAQUE newtype over `Evidence<T>` — private inner, no `DerefMut`, no
   field access for writing. The ONLY constructor is
   `Evidence::certified()`, which validates the ACTUAL numbers, not the
   constructor that claimed them: Exact requires a finite QoI with
   bit-identical bounds; Enclosure requires finite ordered bounds that
   CONTAIN the QoI; statistical e-values, levels, widths, and confidence
   parameters must satisfy their finite domains; model discrepancy must be
   non-negative (positive infinity is the explicit unbounded claim); empty,
   inverted, or non-finite model-validity domains refuse even when a public
   literal asserts `in_domain: true`; Estimate/NoClaim and out-of-domain models
   refuse. Decision breakdown applies the same validity check and assigns an
   infinite model band to an impossible domain.
   Reads flow through `Deref<Target = Evidence<T>>`;
   `Certified::into_evidence()` is the explicit downgrade — the mark is
   lost and any reconstruction must re-enter `certified()`
   (re-validated round trip). Escape hatches are plain `Evidence<T>` or
   `NumericalCertificate::no_claim()`, never a `Certified<T>`.
   MIGRATION: `Certified<T>` was a type alias for `Evidence<T>`;
   callers that mutated or moved fields of a certified value now call
   `.into_evidence()` first (one workspace site: fs-geom conformance),
   and callers that only read fields are unchanged via Deref. This
   crate has no serializer; persisted evidence re-enters through
   `certified()` on ingest by construction.
9. Decision assessment fails closed (evd-013): malformed or negative
   uncertainty becomes an infinite band; infinite totals and malformed
   thresholds cannot become `DecisionGrade`.

## Error model
Structured teaching errors throughout: `CertifyError`, `RegistryError`,
`OutOfDomain`, `FitError` — all `core::error::Error` with actionable
Display text. Constructors are total (enclosure bounds normalize by
swapping); no panics cross the boundary.

## Determinism class
Deterministic: every function is a pure computation over its inputs; all
renderings use sorted (BTreeMap) order; no clocks, no addresses, no
randomness. Bit-stable across runs and platforms up to fs-math-class
scalar-arithmetic divergence.

## Cancellation behavior
No compute loops (bounded small algebra per call); nothing to poll. The
crate is used INSIDE cancellable kernels; it adds no blocking.

## Unsafe boundary
None. `unsafe_code` denied workspace-wide.

## Feature flags
None. The mechanisms are `[S]`-grade bookkeeping; the models they DESCRIBE
carry their own ambition tags in their cards.

## Conformance tests
tests/conformance.rs, cases evd-001..evd-013 (JSON-line verdicts; seeded
cases carry seeds): the G0 conservativeness battery, the registration
lint, the worked model-discrepancy-dominates example (10% closure vs 0.7%
mesh at a 5% threshold → NotDecisionGrade{ModelForm} + escalation advice,
flipping to DecisionGrade with a 2% calibrated closure), the
out-of-distribution refusal on a synthetic two-fidelity corpus, bracketing
spread reporting with deterministic schema-valid rows, and certification
poisoning. The color cases cover disjoint-regime demotion with infinite
dispersion, outward-rounded verified arithmetic, non-finite state/regime
demotion, and deterministic escaping/tagging of hostile JSON payloads.
In-module suites cover the certificate algebra, validity laws, tie-breaking,
provenance chaining, and card rendering.
`evd-013` exercises the public `Evidence` layer against indeterminate interval
arithmetic and malformed statistical/model uncertainty.

- Falsifier registration is total: empty falsifier lists refuse at the
  source; the ship gate names every unpaired class.
- Budget allocation is monotone in consequence AND doubt, with honest
  boundaries (cold-start max, perfect-record floor, dependent-free floor,
  empty-job zero) — property-tested.
- Every falsifier hit produces BOTH a tombstone and an estimator bug
  report; neither is optional.

## No-claim boundaries (colors)

- Verified-interval composition here covers outward-rounded Add/Mul and exact
  endpoint Hull; the full ledger operation algebra composes through fs-ivl
  when wired.
- Estimated dispersion combines additively (conservative); calibrated
  dispersion algebra joins the color-probes bead.

## No-claim boundaries
- Statistical composition is CONSERVATIVE-WEAKEST v1 (half-widths add,
  confidences min, mixed kinds keep the width-bearing certificate);
  proper e-value arithmetic (products under independence, e-BH) is
  fs-eproc's contract and will replace this composition rule behind the
  same API.
- Discrepancy models are honest BOOKKEEPING (observed box + mean/max
  band), not learning; trained/learned discrepancy models (FrankenTorch)
  arrive with fs-surrogate and will implement the same query/refusal
  surface.
- The adjoint hook is carried, never composed here — composed tapes are
  fs-ad's contract.
- First-order band addition is conservative for small relative bands; it
  is NOT a rigorous product-form bound for large ones — fs-ivl composition
  should be used for the numerical slice when bands are large (the
  numerical slice already does interval arithmetic; the TOTAL across
  sources is the first-order sum).
- Ledger persistence is row RENDERING only; the `evidence` /
  `model_cards` tables land with fs-ledger (rows are shaped for that
  migration).
- `ProvenanceHash` is FNV-1a until the BLAKE3-class ledger hash supersedes
  it (same upgrade path as fs-obs).

## No-claim boundaries (falsifiers)

- The registry stores falsifier IDENTITIES and stated methods; executing
  a falsifier (running rays, FD probes, full solves) is each consumer
  kernel's code — this module is the schema, gate, and allocator.
- `consequence` is supplied by the caller (ledger-DAG dependent-weight
  traversal is HELM-side); the allocator's contract is what it does with
  the number, floors included.
- Tombstone/bug payloads are canonical JSON for the ledger; the tombstone
  REGISTRY (Proposal E) and falsifier-log mining (Proposal 9) consume
  them downstream.
- Rent decay is per-class and floor-bounded; per-regime decay and
  quarterly cadence enforcement are governance-bead policy (xpck.6).
