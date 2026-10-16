//! The I08 (evidence-budget co-design planner) VerificationManifest draft
//! (bead frankensim-leapfrog-2026-program-i94v.1.4.8.1).
//!
//! Baseline lattice ([S]): typed evidence actions, calibrated
//! response/cost models, constrained portfolio planning, online
//! idempotent replanning, anti-Goodhart and falsification budgets.
//! Maximal lattice ([F]/[M]): anytime-valid globally adaptive evidence
//! allocation with theorem-aware value of information, certified
//! self-improving cost/error models, and a robust multi-horizon
//! optimality falsifier lane. A weaker receipt closes its own element and
//! is never relabeled as the stronger theorem.

use crate::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin, FixtureSource,
    GauntletTier, ManifestDraft, ObligationRow, OracleRoute, Partition, ToleranceSemantics, Waiver,
};

/// Build the I08 draft. Consumers freeze it themselves; the conformance
/// battery proves it freezes.
#[must_use]
pub fn i08_draft() -> ManifestDraft {
    ManifestDraft {
        initiative: "I08",
        title: "Evidence-budget co-design planner gate: typed actions, calibrated \
                 models, constrained portfolios, idempotent replanning, and \
                 anti-Goodhart budgets",
        version: 1,
        explicits: FiveExplicits {
            units: "SI base units; costs in seconds (wall) and joules where metered, \
                    budgets in the same units; probabilities dimensionless '1'; exact \
                    verdicts use unit 'bit'",
            seeds: "Philox 4x32-10 counter streams keyed 'i08/<fixture-id>/<case-index>'; \
                    development cases use seed indices 0..=4095; held-out cases use \
                    65536..=69631 (disjoint by construction, split frozen here)",
            budgets: "smoke tier <= 60 s on one host; core tier <= 30 min; max tier <= 8 h \
                      on a quiet perf host; <= 16 GiB memory per lane; accuracy budgets \
                      are the per-claim tolerance fields",
            versions: "fs-vmanifest schema v2; toolchain pinned by \
                       rust-toolchain.toml; sibling pins by constellation.lock",
            capabilities: "no network; no FFI; deterministic mode mandatory for every G5 \
                           row; frontier/moonshot lanes stay behind feature flags",
        },
        claims: i08_claims(),
        fixtures: i08_fixtures(),
        obligations: i08_obligations(),
        waivers: i08_waivers(),
        amendment_rules: "Any change is a successor version through FrozenManifest::amend; \
                          the amendment record names every invalidated claim/obligation \
                          descendant; an amendment after campaign start invalidates the \
                          affected evidence, which must be re-earned; there is no in-place \
                          edit path in the type system",
    }
}

#[allow(clippy::too_many_lines)]
fn i08_claims() -> Vec<ClaimSpec> {
    vec![
        ClaimSpec {
            id: "i08-typed-evidence-actions",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Evidence actions (run test, refine model, acquire data, \
                        escalate fidelity) are typed values carrying units, budget \
                        effects, preconditions, and produced-evidence declarations; \
                        ill-typed actions (unit mismatch, negative cost, undeclared \
                        evidence) are refused with a structured diagnostic",
            hypotheses: &[
                "action declarations come from the pinned portfolio fixture families",
                "budget arithmetic is the declared cost unit system in explicits",
            ],
            qoi: "action_admission_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent action type-checker walking the fixture \
                           generator's action certificates (emitted during generation, \
                           not by the planner)",
                independent: true,
                tcb_overlap: "shares fs-qty dimension arithmetic only",
            },
            activation: "first implementation bead of the I08 action leaf opens",
            kill: "any admitted action whose certificate says ill-typed (or vice versa) \
                   returns the action algebra to review with the disagreeing action as \
                   receipt",
            fallback: "narrow the action vocabulary and re-freeze via amendment",
            no_claim: "typing only; no claim that admitted actions are worth running",
        },
        ClaimSpec {
            id: "i08-calibrated-response-cost-models",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Response and cost models fitted on the pinned histories are \
                        calibrated: prediction-interval empirical coverage on the \
                        held-out split lies within the stated band around nominal, \
                        including under the censored-cost decks (censoring handled, \
                        not dropped)",
            hypotheses: &[
                "histories from the pinned portfolio and censored-cost families",
                "nominal coverage levels are the deck-declared grid {0.5, 0.8, 0.9}",
            ],
            qoi: "max_coverage_deviation",
            unit: "1",
            tolerance: ToleranceSemantics::Absolute { atol: 0.05 },
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "direct empirical-coverage counting on the held-out split, \
                           computed without the model-fitting code",
                independent: true,
                tcb_overlap: "shares fs-math deterministic scalar kernels",
            },
            activation: "action leaf green at smoke tier",
            kill: "systematic over-confidence (coverage below band) on any deck class \
                   kills the model family for planning use; models demote to \
                   diagnostics",
            fallback: "worst-case (interval) cost bounds instead of calibrated models",
            no_claim: "calibration on the pinned distributions only; no claim under \
                       distribution shift beyond the budget-shock decks",
        },
        ClaimSpec {
            id: "i08-constrained-portfolio-planning",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Planned portfolios satisfy every declared constraint exactly \
                        (budget caps, dependency edges, mutual exclusions, floor \
                        allocations), and on the small enumerable decks the plan's \
                        declared objective value matches the deck's exhaustively \
                        computed optimum",
            hypotheses: &[
                "portfolios from the pinned cold-plate and motor families",
                "objective and constraints are the deck's frozen declarations",
                "enumerable decks have <= 2^16 feasible portfolios by construction",
            ],
            qoi: "feasibility_and_small_deck_optimality_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "independent constraint checker plus exhaustive enumeration \
                           on the small decks",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "model leaf green at smoke tier",
            kill: "any constraint violation in an emitted plan is a release blocker; \
                   suboptimality on enumerable decks kills the optimality wording \
                   (planner demotes to feasible-plan generator)",
            fallback: "feasible plans with reported (unverified) objective values",
            no_claim: "optimality is claimed only where exhaustive enumeration ran; \
                       large decks get feasibility plus the maximal lanes' bounds",
        },
        ClaimSpec {
            id: "i08-online-idempotent-replanning",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Replanning is deterministic and idempotent: identical planner \
                        state yields byte-identical plans; replanning after an \
                        observation is a pure function of (state, observation, \
                        idempotency key); replaying a logged observation stream \
                        reproduces the logged plan sequence exactly",
            hypotheses: &[
                "observation streams from the pinned optional-stopping schedules",
                "idempotency keys follow the ledger discipline named in explicits",
            ],
            qoi: "replay_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G5,
            oracle: OracleRoute {
                identity: "byte comparison of content-addressed plan artifacts across \
                           replays and worker counts",
                independent: true,
                tcb_overlap: "shares fs-blake3 hashing",
            },
            activation: "planning leaf green at smoke tier",
            kill: "any schedule-dependent plan difference in deterministic mode is a \
                   release blocker",
            fallback: "none: replayability is load-bearing for audit",
            no_claim: "no claim across ISAs; no claim for fast mode",
        },
        ClaimSpec {
            id: "i08-anti-goodhart-falsification-budgets",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Every plan reserves at least the declared falsification floor \
                        for holdout/falsification actions that the optimizing lane \
                        cannot raid under any replanning step, and on the adversarial \
                        Goodhart decks the planner flags proxy gaming within the \
                        deck's declared detection window",
            hypotheses: &[
                "adversarial decks from the pinned Goodhart family with \
                 generator-recorded gaming onset and detection windows",
                "the falsification floor is the deck-declared fraction",
            ],
            qoi: "floor_and_detection_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "budget-ledger auditor (floor accounting) plus the deck's \
                           recorded gaming onset for the detection side",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "planning leaf green at core tier",
            kill: "any replanning step that raids the falsification floor kills \
                   adaptive replanning until redesigned; missed detection past the \
                   window on a majority of decks kills the detection claim",
            fallback: "static falsification allocations fixed at plan time",
            no_claim: "detection is claimed for the pinned gaming taxonomy within \
                       declared windows; no claim against novel gaming strategies",
        },
        ClaimSpec {
            id: "i08-anytime-adaptive-allocation-voi",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Adaptive allocation with theorem-aware value of information: \
                        under the pinned optional-stopping schedules the allocation's \
                        stopping decisions remain anytime-valid (e-process discipline; \
                        no inflation from optional stopping), and allocations agree \
                        with the deck-computed VoI ranking within the stated \
                        top-k band",
            hypotheses: &[
                "decks carry generator-computed VoI rankings under the frozen \
                 theorem-obligation weights",
                "stopping uses the declared e-process family; dependence between \
                 tests follows the pinned dependent-test graphs",
            ],
            qoi: "topk_ranking_agreement",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.8, hi: 1.0 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "the fixture generator's VoI computation plus null-simulation \
                           batteries with adversarial stopping for the anytime-validity \
                           side",
                independent: true,
                tcb_overlap: "shares fs-eproc e-process primitives",
            },
            activation: "all five baseline claims green at core tier",
            kill: "anytime-validity violation in null simulations (error inflation \
                   beyond the declared level) kills adaptive stopping; ranking \
                   agreement below band demotes VoI to advisory",
            fallback: "fixed-allocation plans with declared sample sizes",
            no_claim: "validity under the pinned dependence structures only; VoI \
                       agreement is against the deck's theorem weights, not a claim \
                       of universal informativeness",
        },
        ClaimSpec {
            id: "i08-self-improving-cost-models",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Online cost/error-model updates never degrade held-out \
                        calibration: after each pinned update batch, held-out coverage \
                        deviation is no worse than the pre-update deviation plus the \
                        stated slack, with the update ledgered and reversible",
            hypotheses: &[
                "update batches from the pinned histories in deck order",
                "reversibility means the ledger can reconstruct any prior model \
                 version bit-exactly",
            ],
            qoi: "max_calibration_regression",
            unit: "1",
            tolerance: ToleranceSemantics::Absolute { atol: 0.02 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "held-out coverage counting per update step, computed without \
                           the update code; ledger replay for reversibility",
                independent: true,
                tcb_overlap: "shares fs-math deterministic scalar kernels",
            },
            activation: "all five baseline claims green at core tier",
            kill: "a calibration regression beyond slack that persists across two \
                   update batches kills online updates (models freeze at plan time)",
            fallback: "frozen models refit only at campaign boundaries",
            no_claim: "non-degradation on the pinned update sequences; no claim of \
                       improvement, and no claim under adversarial data injection",
        },
        ClaimSpec {
            id: "i08-robust-multihorizon-optimality",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Refutation,
            statement: "Falsifier lane: search the frozen portfolio space, across the \
                        pinned horizons and budget-shock scenarios, for a plan whose \
                        worst-case declared objective beats the planner's by more than \
                        the declared robustness factor; the optimality claim survives \
                        only while this search comes up empty",
            hypotheses: &[
                "portfolio space, horizons, and shock scenarios are the frozen \
                 enumerations in the fixture families",
                "search budget is the max-tier budget in explicits",
            ],
            qoi: "counterexample_count",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent enumerator/branch-and-bound over the frozen \
                           portfolio space with the deck's worst-case evaluator",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "maximal lane opens only after baseline closes at core tier",
            kill: "first verified counterexample kills the robust-optimality claim; \
                   the beating plan is retained as the receipt",
            fallback: "robustness restated as within-factor bound via amendment, or \
                       the planner demotes to heuristic with reported worst cases",
            no_claim: "an empty search is Estimated evidence bounded by the frozen \
                       space and budget; it is never an optimality proof",
        },
    ]
}

fn i08_fixtures() -> Vec<FixturePin> {
    vec![
        FixturePin {
            id: "i08-cold-plate-portfolio",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: cold-plate thermal-design evidence portfolios: \
candidate tests (correlation lookups, reduced-order runs, \
full-fidelity solves, physical-test slots) with declared costs, \
response distributions, dependency edges, and per-deck frozen \
objective/constraint sets; small decks (<= 2^16 feasible \
portfolios) carry exhaustively computed optima.\n\
SEEDS: Philox stream 'i08/coldplate/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i08-motor-portfolio",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: motor electromagnetic/thermal evidence portfolios \
with the same declaration discipline as the cold-plate family \
(costs, responses, dependencies, frozen objectives, small-deck \
optima).\n\
SEEDS: Philox stream 'i08/motor/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i08-portfolio-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "Held-out variants of both portfolio generators; SEEDS: Philox \
streams 'i08/coldplate/<k>' and 'i08/motor/<k>', k in \
65536..=69631; withheld until adjudication.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i08-censored-cost-histories",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: cost/response histories with right-censoring \
(timeouts, cancelled runs) at recorded censoring times, so \
calibration must handle censoring rather than dropping rows.\n\
SEEDS: Philox stream 'i08/censored/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i08-dependent-test-graphs",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: dependence structures between evidence tests \
(shared fixtures, shared solvers, common-mode failures) as \
declared dependence graphs consumed by the anytime-validity \
lane.\n\
SEEDS: Philox stream 'i08/dependent/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i08-optional-stopping-schedules",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: adversarial optional-stopping observation schedules \
(peeking patterns, early-stop temptations, null streams for \
validity batteries) with recorded ground truth per stream.\n\
SEEDS: Philox stream 'i08/stopping/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i08-goodhart-adversaries",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: adversarial Goodhart loops: response processes that \
reward the proxy objective while degrading the true objective \
from a recorded gaming onset, with declared detection windows.\n\
SEEDS: Philox stream 'i08/goodhart/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i08-goodhart-adversaries-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "Same generator as i08-goodhart-adversaries; SEEDS: Philox \
stream 'i08/goodhart/<k>', k in 65536..=69631; withheld until \
adjudication.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i08-budget-shock-scenarios",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: mid-campaign budget shocks (cuts, windfalls, \
deadline moves) as declared scenario schedules for the robust \
multi-horizon lane.\n\
SEEDS: Philox stream 'i08/shocks/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn i08_obligations() -> Vec<ObligationRow> {
    const UNIT_CASES: &[&str] = &[
        "happy",
        "empty",
        "boundary",
        "max",
        "error",
        "unit-dimension",
        "tie-break",
        "cancellation",
        "migration",
    ];
    vec![
        ObligationRow {
            leaf: "i08-action-algebra",
            claims_covered: &["i08-typed-evidence-actions"],
            unit_cases: UNIT_CASES,
            g0: "generators: portfolio families incl. ill-typed action twins; laws: \
                 admission/refusal certificate agreement, budget-effect additivity; \
                 shrinker: action-set reduction preserving the fault; replay seeds per \
                 explicits",
            decks: &["i08-cold-plate-portfolio", "i08-motor-portfolio"],
            g3_relations: &["action relabeling invariance"],
            g4_schedule: "cancel admission at every boundary class once; drain leaves \
                          no partial action registration",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i08_actions.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i08-actions slice)",
            obs_events: &["action.admitted", "action.refused", "action.cancelled"],
            replay_command: "scripts/e2e/leapfrog/i08_actions.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i08-model-calibration",
            claims_covered: &["i08-calibrated-response-cost-models"],
            unit_cases: UNIT_CASES,
            g0: "generators: portfolio + censored-cost histories; laws: coverage \
                 counting on held-out split, censored rows retained; replay seeds per \
                 explicits",
            decks: &["i08-censored-cost-histories", "i08-portfolio-holdout"],
            g3_relations: &["unit-rescaling covariance of cost models"],
            g4_schedule: "cancel mid-fit; drained fits publish no partial model",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i08_models.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i08-models slice)",
            obs_events: &["model.fitted", "model.coverage", "model.cancelled"],
            replay_command: "scripts/e2e/leapfrog/i08_models.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i08-portfolio-planning",
            claims_covered: &["i08-constrained-portfolio-planning"],
            unit_cases: UNIT_CASES,
            g0: "generators: portfolio families; laws: constraint satisfaction on \
                 every emitted plan, exact optimum agreement on enumerable decks; \
                 replay seeds per explicits",
            decks: &["i08-cold-plate-portfolio", "i08-motor-portfolio"],
            g3_relations: &["constraint permutation invariance of feasibility"],
            g4_schedule: "cancel mid-plan; drained planning publishes no partial plan",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i08_planning.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i08-planning slice)",
            obs_events: &["plan.emitted", "plan.checked", "plan.cancelled"],
            replay_command: "scripts/e2e/leapfrog/i08_planning.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i08-replanning-replay",
            claims_covered: &["i08-online-idempotent-replanning"],
            unit_cases: UNIT_CASES,
            g0: "generators: optional-stopping observation streams; laws: idempotence \
                 on identical state, logged-stream replay equality, idempotency-key \
                 discipline; replay seeds per explicits",
            decks: &["i08-optional-stopping-schedules"],
            g3_relations: &["observation-batch regrouping invariance (same content)"],
            g4_schedule: "storm: cancel/resume replanning at every boundary class, \
                          seeded injection order",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i08_replanning.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i08-replanning slice)",
            obs_events: &["replan.step", "replay.artifact", "replan.cancelled"],
            replay_command: "scripts/e2e/leapfrog/i08_replanning.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i08-goodhart-defense",
            claims_covered: &["i08-anti-goodhart-falsification-budgets"],
            unit_cases: UNIT_CASES,
            g0: "generators: Goodhart adversaries; laws: falsification floor never \
                 raided across replans, detection within declared windows; replay \
                 seeds per explicits",
            decks: &["i08-goodhart-adversaries"],
            g3_relations: &["proxy relabeling does not evade detection"],
            g4_schedule: "cancel during detection; resumed runs report identical \
                          detection steps",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i08_goodhart.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i08-goodhart slice)",
            obs_events: &["floor.audited", "gaming.flagged", "defense.cancelled"],
            replay_command: "scripts/e2e/leapfrog/i08_goodhart.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i08-maximal-lanes",
            claims_covered: &[
                "i08-anytime-adaptive-allocation-voi",
                "i08-self-improving-cost-models",
                "i08-robust-multihorizon-optimality",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: dependent-test graphs, stopping schedules, shock \
                 scenarios, held-out families; laws: null-simulation validity under \
                 adversarial stopping, VoI ranking agreement, per-update calibration \
                 non-degradation with ledger reversibility, portfolio enumeration \
                 replay; replay seeds per explicits",
            decks: &[
                "i08-dependent-test-graphs",
                "i08-optional-stopping-schedules",
                "i08-budget-shock-scenarios",
                "i08-portfolio-holdout",
                "i08-goodhart-adversaries-holdout",
                "i08-industrial-evidence-portfolio",
            ],
            g3_relations: &[
                "held-out vs development distribution parity checks",
                "shock-scenario permutation invariance of worst-case evaluation",
            ],
            g4_schedule: "cancel each maximal lane once mid-campaign; a resumed lane \
                          must report identical verdicts",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i08_maximal.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (i08-maximal slice)",
            obs_events: &[
                "allocation.step",
                "voi.ranking",
                "model.update",
                "optimality.counterexample",
            ],
            replay_command: "scripts/e2e/leapfrog/i08_maximal.sh --replay <artifact-id>",
        },
    ]
}

fn i08_waivers() -> Vec<Waiver> {
    vec![Waiver {
        subject: "i08-industrial-evidence-portfolio",
        reason: "the external industrial evidence-portfolio case study (a vendor \
                 test-campaign record) is not yet licensed or content-addressed; only \
                 its slot is preregistered",
        owner: "I08 implementation beads (frankensim-leapfrog-2026-program-i94v.1.4.x)",
        predicate: "deck bytes pinned with exact edition, license, digest, QoIs, and \
                    acceptance envelopes through the fs-vvreg registry discipline",
        expiry: "before the first Max-tier campaign run; review at each Phase-2 close \
                 burst",
        promotion_effect: "maximal claims may run without this deck but cannot cite it; \
                           any result depending on it stays Estimated until the pin \
                           lands",
    }]
}
