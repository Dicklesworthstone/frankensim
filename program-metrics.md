# FrankenSim program metrics dashboard

schema: 1
metrics: 19
measured: 8
no_data: 11
trend_basis: NO-DATA (no generation recorded yet; every trend cell reads `no prior generation`)
source_identity: 611f1afc9ba737a492271bc4484726d947876d0d5e5c309175d24e7af0c3c4a2

This dashboard measures OUTCOMES. A `NO-DATA` row means no measurement machinery exists yet, so no number is invented; a measured `0` means the population is enumerable and the answer is genuinely none. The two are never conflated, and a measured zero is deliberately left visible rather than hidden behind `NO-DATA`.

## Outcome metrics

| metric | value | trend | direction | sources |
| --- | --- | --- | --- | --- |
| Prediction error against blind held-out experimental references | NO-DATA (needs f85xj.7.5) | no prior generation | lower-is-better | - |
| Share of claims replayed bitwise on an independent machine | NO-DATA | no prior generation | higher-is-better | - |
| Compliance verdicts flipped by adding a previously omitted budget term | NO-DATA (needs f85xj.8.7) | no prior generation | neutral | - |
| Time from dirty CAD input to a defensible decision artifact | NO-DATA (needs f85xj.6.11) | no prior generation | lower-is-better | - |
| Empirical coverage of predicted uncertainty intervals | NO-DATA (needs f85xj.7.2) | no prior generation | higher-is-better | - |
| Share of user-facing outputs carrying a complete error budget | NO-DATA (needs f85xj.8.7) | no prior generation | higher-is-better | - |
| Rate at which adversarial challenges were wrongly accepted | NO-DATA (needs f85xj.7.5) | no prior generation | lower-is-better | - |
| Supplier CAD import: clean, repaired, and refused rates | NO-DATA (needs f85xj.11.6) | no prior generation | higher-is-better | - |
| Datasets reproduced by an independent team or implementation lineage | 0 of 28 (0.00%) | no prior generation | higher-is-better | fs_vvreg::corpus seeded validation registry |
| Correctness of certify-or-escalate decisions by learned components | NO-DATA (needs f85xj.14.2) | no prior generation | higher-is-better | - |
| Time to explain a surprising result through ledger lineage | NO-DATA | no prior generation | lower-is-better | - |
| Setup time, diagnosis time, and decision quality from real user sessions | NO-DATA (needs f85xj.7.6) | no prior generation | higher-is-better | - |

## Evidence portfolio metrics

| metric | value | trend | direction | sources |
| --- | --- | --- | --- | --- |
| Registered adversarial challenges actually executed | 0 of 8 (0.00%) | no prior generation | higher-is-better | fs_vvreg::adversarial registry; vv-scorecard.json (fs_vvreg::scorecard) |
| Validation datasets on the blind-predictive axis | 0 of 28 (0.00%) | no prior generation | higher-is-better | fs_vvreg::corpus seeded validation registry |
| Validation datasets supplying an external evidence axis | 9 of 28 (32.14%) | no prior generation | higher-is-better | fs_vvreg::corpus seeded validation registry |
| Claim cells carrying at least one external reference | 13 of 25 (52.00%) | no prior generation | higher-is-better | vv-scorecard.json (fs_vvreg::scorecard) |
| Claim cells with no external reference at all | 12 | no prior generation | lower-is-better | vv-scorecard.json (fs_vvreg::scorecard) |

## Governance metrics

| metric | value | trend | direction | sources |
| --- | --- | --- | --- | --- |
| Registered capabilities at L2 (numerically verified) or above | 11 of 15 (73.33%) | no prior generation | higher-is-better | capability-maturity.json |
| Registered capabilities at L3 (integrated workflow) or above | 0 of 15 (0.00%) | no prior generation | higher-is-better | capability-maturity.json |

## What each metric does not capture

- `blind-prediction-error` — error against a reference bounds nothing on its own: the reference's own uncertainty and the regime it was measured in both constrain what the number means
- `cross-machine-reproducibility` — bitwise replay on a second machine proves determinism, not correctness: two machines can reproduce the same wrong number exactly
- `decision-changes-from-omitted-uncertainty` — a flipped verdict shows the term mattered on re-run; it does not show the new verdict is correct, only that the earlier one was underdetermined
- `decision-turnaround` — when live this will time OUR examples on OUR hardware, which is not the same as a real user's setup, data, or interruptions
- `empirical-interval-coverage` — coverage measured on the calibration population says nothing about coverage under the distribution shift a real design study introduces
- `error-budget-completeness` — a complete budget means every term the model KNOWS about is present; it cannot count terms nobody has identified yet
- `false-acceptance-rate` — a false-acceptance rate only covers the failure modes someone thought to write a challenge for; it is silent about unimagined ones
- `import-admission-rate` — import success measures admission, not fidelity: a file that imports cleanly can still carry geometry that means something different downstream
- `independent-reproduction` — this counts datasets DECLARING the independent-reproduction axis; the declaration is a registry fact, and the current value is a genuine zero rather than an unmeasured one
- `surrogate-escalation-correctness` — escalation correctness is measured against cases where the truth is known, which are systematically the easier ones
- `time-to-explain` — even when live this measures the tool's explain path, not whether the explanation actually convinced the engineer reading it
- `user-study-measurements` — proxies measured on people who already know the system are the opposite of the population this metric is supposed to describe
- `adversarial-suite-execution` — execution is not survival: running a challenge says nothing about whether the program passed it, and an unexecuted suite is a registry of good intentions
- `blind-predictive-datasets` — the honesty exam: a prediction made before the reference was unblinded. The current value is a real zero, which is the single most important number on this dashboard
- `external-reference-datasets` — external means cross-code, controlled-experimental, blind-predictive, or field monitoring; our own numerical verification is deliberately excluded from the numerator because agreeing with ourselves is not external evidence
- `externally-anchored-claim-cells` — anchoring counts REFERENCES attached to a cell, not agreement with them; a cell can be anchored and still predict the reference badly
- `unanchored-claim-cells` — an absolute count, so it grows when the program declares new claim cells; a rising number can mean expanding scope rather than decaying evidence
- `capabilities-at-l2-plus` — registry levels are declarations backed by cited evidence, not independent audits; the registry's own maturity is L1
- `capabilities-at-l3-plus` — L3 requires an admitted end-to-end integration claim; the current value is a real zero, and no crate count or test count can move it

## Why metrics are missing

- `blind-prediction-error` — no ledgered run-result store exists and no Level-D blind reference is admitted, so no model-versus-reality error can be computed (tracked: f85xj.7.5)
- `cross-machine-reproducibility` — cross-ISA determinism is proven per-artifact by golden couplings and per-host by perf baselines, but no lane aggregates them into a program-level replay rate
- `decision-changes-from-omitted-uncertainty` — no ledger diff of decision verdicts across budget-term introductions is retained; this is the single best evidence that the error-budget program has decision value, and it is not yet collected (tracked: f85xj.8.7)
- `decision-turnaround` — no end-to-end acceptance lane records stage timings, so there is nothing to measure without inventing it (tracked: f85xj.6.11)
- `empirical-interval-coverage` — the empirical coverage machinery is not live; nominal coverage is never extrapolated into an empirical claim (tracked: f85xj.7.2)
- `error-budget-completeness` — no audit enumerates user-facing outputs against their budget terms (tracked: f85xj.8.7)
- `false-acceptance-rate` — zero registered adversarial challenges have been executed, so the rate's denominator is empty; a rate over zero trials is unrepresentable, not zero (tracked: f85xj.7.5)
- `import-admission-rate` — no retained real supplier CAD corpus exists, and rates measured on fixtures we authored would be self-graded (tracked: f85xj.11.6)
- `surrogate-escalation-correctness` — no decisive-metrics instrumentation exists for learned components (tracked: f85xj.14.2)
- `time-to-explain` — the ledger exposes no explain-query session instrumentation, so no timing surface exists to read
- `user-study-measurements` — no user-study measurement exists; the nearest current proxies are quickstart timings and the external-reproduction friction log, neither of which is a user study (tracked: f85xj.7.6)

## Deliberately excluded

These are legitimate signals that are NOT outcome metrics. They move without the program getting better at predicting reality, and each has its own lane.

- kernel throughput — a performance diagnostic owned by the roofline lane; a faster kernel is not a better prediction
- crate count — inventory, not capability; the capability maturity registry is the outcome measure
- integration-test file count — inventory, not proof; check-docs already pins it and a test file is not an outcome
- open issue counts — the beads store churns on every unrelated issue edit, which would make this checked artifact stale for every agent in the repository

identity: 2ce9f333a951f45151843f39be493b837731f9ed7ac3ce534d198651f460e251
