# External Evidence Is a Portfolio, Not a Pyramid

Status: ratified for the EXTREAL validation corpus and claim-admission
scorecards. Schema authority:
`fs_vvreg::portfolio::EVIDENCE_PORTFOLIO_SCHEMA_VERSION`.

## Ruling

External evidence is represented on independent categorical axes. The axes do
not form a ladder, are never averaged into a confidence percentage, and do not
admit a `max(level)` operation. More observations on one axis cannot manufacture
a missing axis.

This matters because evidence answers different questions. A carefully designed
laboratory experiment may control a claim's nuisance variables better than a
large field dataset. A blind result tests a frozen prediction protocol but does
not repair weak instrumentation. Agreement between two codes can expose an
implementation defect but can also preserve a shared modeling misconception.
Field monitoring can reveal operational drift while remaining confounded about
its cause.

## The seven coordinates

| Axis | It can support | It cannot support by itself | Canonical failure modes |
| --- | --- | --- | --- |
| Numerical verification | Correct solution of the declared equations, discretization/order claims, algebraic convergence, and implementation consistency | Physical adequacy of the equations or parameters | Shared oracle, inverse crime, unretained meshes, false asymptotic range, unchecked roundoff |
| Cross-code agreement | Detection of implementation-specific mistakes and sensitivity to an independent numerical path | Truth of a model or physical validation | Common authorship, common libraries, copied formulas, matched defaults, correlated bugs |
| Controlled experimental validation | Predictive adequacy for a named QoI in a measured and controlled regime | Blindness, field transfer, population-wide validity, or independent reproduction | Calibration gaps, hidden tuning, incomplete nuisance control, condition mismatch, selection bias |
| Blind predictive validation | Performance after the prediction and split were frozen | Experimental quality or physical authority when the blind source is not itself controlled | Leaked holdout, post-freeze rule changes, weak preregistration, selective reporting |
| Field monitoring | Behavior and drift under operational conditions | A `Validated` prediction color without controlled-condition provenance | Unknown loads, sensor health, intervention effects, latent confounders, survivorship bias |
| Transferability across regimes | Stability of a claim across explicitly different regimes | Validity outside the tested transfer graph | Sparse regime coverage, covariate shift, hidden extrapolation, changing instrumentation |
| Independent reproduction | Resistance to one team's implementation and analysis choices | Automatic correctness if both reproductions share a source or assumption | Shared independence group, unavailable raw data, copied preprocessing, publication selection |

## A-E corpus migration

The retained A-E tags remain compact source-corpus encodings, but schema v3
interprets them as coordinates:

| Legacy tag | Portfolio coordinates |
| --- | --- |
| A | numerical verification |
| B | cross-code agreement |
| C | controlled experimental validation |
| D | controlled experimental validation + blind predictive validation |
| E | field monitoring |

The letters are labels, not ranks. There is deliberately no ordering relation
on `EvidenceLevel`. Transferability and independent reproduction have no
legacy-letter shortcut; they require explicit portfolio observations.

## Claim admission

Claim admission is exact in claim class, QoI, and regime. Every required axis
must have a same-QoI, same-regime observation:

| Claim class | Required coordinates |
| --- | --- |
| Numerically verified | numerical verification |
| Cross-code consistent | numerical verification + cross-code agreement |
| Validated prediction | controlled experimental validation |
| Blind validated prediction | controlled experimental validation + blind predictive validation |
| Field-supported prediction | controlled experimental validation + field monitoring |
| Transferable prediction | controlled experimental validation + transferability across regimes |
| Independently reproduced prediction | controlled experimental validation + independent reproduction under a distinct source and independence group |

The structural admission receipt is not source authentication and does not mint
an evidence color. Authentication, exact-instance policy, and durable authority
remain package/checker/ledger responsibilities. It only proves that the
non-substitutable coordinate rule ran over the attached content identities.

## Anti-laundering rules

1. Field monitoring alone cannot support `Validated`, regardless of dataset
   size.
2. Repeating, splitting, or relabeling one observation does not increase axis
   coverage.
3. Numerical verification and cross-code agreement do not become experimental
   validation when composed.
4. Controlled experiments do not become blind evidence without a frozen split
   and release.
5. Independent reproduction requires a source and independence group both
   distinct from the controlled experiment.
6. Evidence for a different QoI or regime cannot satisfy the requested claim.
7. Process or standards conformance stays separate from every scientific axis.

These rules mirror the weakest-wins and no-authority-by-composition rules in
`fs-evidence`. Composition can preserve already admitted authority; it cannot
create a missing external anchor.

## Scorecard behavior

The deterministic corpus audit renders all seven counts for every tracked
thermal QoI. It never chooses a highest letter or averages categorical axes.
Zeroes remain visible. The current acquisition warning gate specifically names
missing controlled-experiment coverage, while the complete row exposes thin
blind, field, transferability, and reproduction coordinates for planning.

### Worked refusal: field data are not a lab surrogate

Suppose a year of rack telemetry reports component peak temperature under an
unknown mixture of fan maintenance, workload scheduling, recirculation, and
sensor replacement. Even with original bytes, calibrated current sensors, and
complete lineage, this contributes the `field-monitoring` coordinate only.

A request for `ValidatedPrediction(component-peak-temperature,
cooling-vertical-v1)` is refused with:

```text
missing axis controlled-experimental-validation
```

The data remain useful: they can detect drift, constrain an operational
envelope, or support a field-supported claim after a controlled experiment is
attached. The refusal prevents operational volume from laundering unresolved
confounding into physical-validation authority.

## Review and falsifier

Review this doctrine when a retained case demonstrates that two coordinates do
not remain sufficiently distinct in what they support, or when a claim class
needs a materially different requirement set. This doctrine does not assert
probabilistic independence among observations or axes. A proposed change is
falsified if it allows any absent required coordinate to be satisfied by
duplicating or composing other axes, or if equivalent input order changes the
portfolio/admission identity.

Changes require a schema/domain rotation, G0 missing-axis tests, G3
order/mutation tests, and an updated worked refusal. New axes may be added, but
no new inference rule enters the default policy merely because it is
mathematically attractive.

## No-claim boundaries

- Axis presence does not authenticate source bytes, instruments, laboratories,
  independence, or scientific correctness.
- Dataset counts are coverage inventory, not sample size, statistical power,
  population coverage, or confidence.
- A controlled-experiment coordinate is necessary for the defined validated
  prediction classes, not sufficient for an organizational deployment
  decision.
- The scorecard does not price evidence, rank experiments, infer dependence, or
  choose the next action.
- The current A-E migration binds only the exact mappings above. It does not
  infer transferability or independent reproduction from publication labels.
