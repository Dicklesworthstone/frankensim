# FrankenSim program metrics dashboard

schema: 1
metrics: 27
measured: 13
no_data: 14
trend_basis: NO-DATA (no generation recorded yet; every trend cell reads `no prior generation`)
source_identity: eee9a8580ac0a5eb45ebd58c5c3ffbb702b7dbc8fdc234a161d1805bba68df47

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
| Human-reviewed supplier CAD files admitted cleanly | NO-DATA (needs f85xj.11.6) | no prior generation | higher-is-better | - |
| Human-locked supplier import annotations that disagree with current observations | 0 | no prior generation | lower-is-better | fs_io::supplier_corpus + data/cad-import-corpus/{corpus-v1.tsv,scorecard-summary-v1.json} |
| Human-reviewed supplier CAD files refused by the standing import policy | NO-DATA (needs f85xj.11.6) | no prior generation | neutral | - |
| Human-reviewed supplier CAD files admitted after repair | NO-DATA (needs f85xj.11.6) | no prior generation | neutral | - |
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
| Open beads with no open blocker, from the tracker snapshot | 288 | no prior generation | neutral | spine-metrics.json (xtask spine-metrics beads snapshot) |
| Open beads with at least one open blocker, from the tracker snapshot | 1456 of 1744 (83.48%) | no prior generation | lower-is-better | spine-metrics.json (xtask spine-metrics beads snapshot) |
| Registered capabilities at L2 (numerically verified) or above | 11 of 15 (73.33%) | no prior generation | higher-is-better | capability-maturity.json |
| Registered capabilities at L3 (integrated workflow) or above | 0 of 15 (0.00%) | no prior generation | higher-is-better | capability-maturity.json |
| Reality-check spine beads on the certified tropical critical path | 0 of 5 (0.00%) | no prior generation | higher-is-better | tropical-critical-path.json (fs-tropical over the bead graph, xtask tropical-path) |
| Staged-producer e2e lane stages proven green by a retained checked receipt | NO-DATA (needs frankensim-iakds) | no prior generation | higher-is-better | - |
| Solve pipeline stages executing in the checked spine ratchet | 3 of 6 (50.00%) | no prior generation | higher-is-better | spine-ratchet.json (fs-cli SolveStage table, xtask spine-ratchet) |

## What each metric does not capture

- `blind-prediction-error` — error against a reference bounds nothing on its own: the reference's own uncertainty and the regime it was measured in both constrain what the number means
- `cross-machine-reproducibility` — bitwise replay on a second machine proves determinism, not correctness: two machines can reproduce the same wrong number exactly
- `decision-changes-from-omitted-uncertainty` — a flipped verdict shows the term mattered on re-run; it does not show the new verdict is correct, only that the earlier one was underdetermined
- `decision-turnaround` — when live this will time OUR examples on OUR hardware, which is not the same as a real user's setup, data, or interruptions
- `empirical-interval-coverage` — coverage measured on the calibration population says nothing about coverage under the distribution shift a real design study introduces
- `error-budget-completeness` — a complete budget means every term the model KNOWS about is present; it cannot count terms nobody has identified yet
- `false-acceptance-rate` — a false-acceptance rate only covers the failure modes someone thought to write a challenge for; it is silent about unimagined ones
- `import-admission-rate` — this rate covers only the retained, human-reviewed population; clean admission is not geometry fidelity, and a file can import cleanly while meaning something different downstream
- `import-annotation-regressions` — this is an absolute locked-annotation mismatch count; zero while reviewed is zero means no reviewed regression exposure, not a validated importer
- `import-refusal-rate` — refusal depends on both corpus difficulty and policy strictness; lowering it by weakening quarantine would not be an improvement
- `import-repair-rate` — repair is neither inherently good nor bad: it records that the standing structural policy changed an input before promotion, not that the repaired geometry is equivalent to the supplier's intent
- `independent-reproduction` — this counts datasets DECLARING the independent-reproduction axis; the declaration is a registry fact, and the current value is a genuine zero rather than an unmeasured one
- `surrogate-escalation-correctness` — escalation correctness is measured against cases where the truth is known, which are systematically the easier ones
- `time-to-explain` — even when live this measures the tool's explain path, not whether the explanation actually convinced the engineer reading it
- `user-study-measurements` — proxies measured on people who already know the system are the opposite of the population this metric is supposed to describe
- `adversarial-suite-execution` — execution is not survival: running a challenge says nothing about whether the program passed it, and an unexecuted suite is a registry of good intentions
- `blind-predictive-datasets` — the honesty exam: a prediction made before the reference was unblinded. The current value is a real zero, which is the single most important number on this dashboard
- `external-reference-datasets` — external means cross-code, controlled-experimental, blind-predictive, or field monitoring; our own numerical verification is deliberately excluded from the numerator because agreeing with ourselves is not external evidence
- `externally-anchored-claim-cells` — anchoring counts REFERENCES attached to a cell, not agreement with them; a cell can be anchored and still predict the reference badly
- `unanchored-claim-cells` — an absolute count, so it grows when the program declares new claim cells; a rising number can mean expanding scope rather than decaying evidence
- `beads-actionable` — directionless on purpose: actionable count falls when work completes AND when new blocked work is filed, so neither rise nor fall is inherently good
- `beads-blocked-ratio` — a deliberately regenerated snapshot ratio; the live tracker moves on every br op and is not a checked input, so this row trails the tracker by design
- `capabilities-at-l2-plus` — registry levels are declarations backed by cited evidence, not independent audits; the registry's own maturity is L1
- `capabilities-at-l3-plus` — L3 requires an admitted end-to-end integration claim; the current value is a real zero, and no crate count or test count can move it
- `spine-critical-path-positions` — slack is computed over ESTIMATES with a recorded default for unestimated beads; a spine bead off the path buys nothing by being rushed, and one on it sets the makespan
- `spine-e2e-lane-green` — an out-of-band green run is not a measurement this artifact can cite: the dashboard reads checked inputs, so the honest row today is the gap itself
- `spine-stages-executing` — counts the executing stage PREFIX the product admits and the ratchet pins; it proves the stages execute, not that their answers are correct

## Why metrics are missing

- `blind-prediction-error` — no ledgered run-result store exists and no Level-D blind reference is admitted, so no model-versus-reality error can be computed (tracked: f85xj.7.5)
- `cross-machine-reproducibility` — cross-ISA determinism is proven per-artifact by golden couplings and per-host by perf baselines, but no lane aggregates them into a program-level replay rate
- `decision-changes-from-omitted-uncertainty` — no ledger diff of decision verdicts across budget-term introductions is retained; this is the single best evidence that the error-budget program has decision value, and it is not yet collected (tracked: f85xj.8.7)
- `decision-turnaround` — no end-to-end acceptance lane records stage timings, so there is nothing to measure without inventing it (tracked: f85xj.6.11)
- `empirical-interval-coverage` — the empirical coverage machinery is not live; nominal coverage is never extrapolated into an empirical claim (tracked: f85xj.7.2)
- `error-budget-completeness` — no audit enumerates user-facing outputs against their budget terms (tracked: f85xj.8.7)
- `false-acceptance-rate` — zero registered adversarial challenges have been executed, so the rate's denominator is empty; a rate over zero trials is unrepresentable, not zero (tracked: f85xj.7.5)
- `import-admission-rate` — no supplier import annotation is human-locked yet; proposed outcomes are not a rate denominator (tracked: f85xj.11.6)
- `import-refusal-rate` — no supplier import annotation is human-locked yet; proposed outcomes are not a rate denominator (tracked: f85xj.11.6)
- `import-repair-rate` — no supplier import annotation is human-locked yet; proposed outcomes are not a rate denominator (tracked: f85xj.11.6)
- `surrogate-escalation-correctness` — no decisive-metrics instrumentation exists for learned components (tracked: f85xj.14.2)
- `time-to-explain` — the ledger exposes no explain-query session instrumentation, so no timing surface exists to read
- `user-study-measurements` — no user-study measurement exists; the nearest current proxies are quickstart timings and the external-reproduction friction log, neither of which is a user study (tracked: f85xj.7.6)
- `spine-e2e-lane-green` — the staged-producer e2e lane runs green out-of-band (34/34 full profile, 2026-07-29, RCH) but retains no tracked checked receipt for this dashboard to read; 'lane not built', 'lane run with zero passing stages', and 'lane green with no retained receipt' are three different facts and only the last is true (tracked: frankensim-iakds)

## Deliberately excluded

These are legitimate signals that are NOT outcome metrics. They move without the program getting better at predicting reality, and each has its own lane.

- kernel throughput — a performance diagnostic owned by the roofline lane; a faster kernel is not a better prediction
- crate count — inventory, not capability; the capability maturity registry is the outcome measure
- integration-test file count — inventory, not proof; check-docs already pins it and a test file is not an outcome

identity: 62ce0177406efec7e78cbe51792f79e59fe1cbc74460950707ff02b98757a6f8
