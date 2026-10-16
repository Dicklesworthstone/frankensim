//! The I01 (multi-field equation/compiler) VerificationManifest draft
//! (bead frankensim-leapfrog-2026-program-i94v.1.1.8.1).
//!
//! Baseline lattice ([S]): typed law/port semantics, residual/JVP/VJP
//! actions, block lowering, conservation and exact-sequence obligations,
//! deterministic cancellation-correct replay. Maximal lattice ([F]/[M]):
//! exact-discrete adjoints and DWR, representation-independent coupled
//! action equivalence, generated multiphysics law completeness, and
//! solver-hint optimality under certified budgets. A weaker receipt closes
//! its own element and is never relabeled as the stronger theorem.

use crate::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin, FixtureSource,
    GauntletTier, ManifestDraft, ObligationRow, OracleRoute, Partition, ToleranceSemantics, Waiver,
};

/// Build the I01 draft. Freezing it (and every gate that entails) is
/// proven by the conformance battery; consumers freeze it themselves so
/// no panic path hides inside a static initializer.
#[must_use]
pub fn i01_draft() -> ManifestDraft {
    ManifestDraft {
        initiative: "I01",
        title: "Multi-field equation/compiler gate: typed laws to lowered blocks with \
                 conservation, derivative, and replay obligations",
        version: 1,
        explicits: FiveExplicits {
            units: "SI base units throughout; dimensionless QoIs use unit '1'; exact \
                    bitwise/boolean verdicts use unit 'bit'",
            seeds: "Philox 4x32-10 counter streams keyed 'i01/<fixture-id>/<case-index>'; \
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
        claims: i01_claims(),
        fixtures: i01_fixtures(),
        obligations: i01_obligations(),
        waivers: i01_waivers(),
        amendment_rules: "Any change is a successor version through FrozenManifest::amend; \
                          the amendment record names every invalidated claim/obligation \
                          descendant; an amendment after campaign start invalidates the \
                          affected evidence, which must be re-earned; there is no in-place \
                          edit path in the type system",
    }
}

// One function per manifest section keeps the data readable; the claim
// table is long because every claim carries its full preregistration row.
#[allow(clippy::too_many_lines)]
fn i01_claims() -> Vec<ClaimSpec> {
    vec![
        ClaimSpec {
            id: "i01-typed-law-port-semantics",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "The compiler admits exactly the typed constitutive-law and port \
                        declarations whose dimensions, variances, and port directions \
                        check, and refuses every ill-typed declaration with a structured \
                        diagnostic naming the offending binding",
            hypotheses: &[
                "law/port declarations come from the pinned constitutive-graph fixture family",
                "dimension algebra is the fs-qty six/seven-base system named in explicits",
            ],
            qoi: "admission_verdict_agreement",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent declaration type-checker walking the fixture \
                           generator's own typing certificate (generated with, not by, \
                           the production admission path)",
                independent: true,
                tcb_overlap: "shares fs-qty dimension arithmetic only",
            },
            activation: "first implementation bead of the I01 admission leaf opens",
            kill: "if >5% of generated well-typed fixtures are refused for reasons the \
                   certificate cannot explain, the admission design returns to review \
                   with the refusal corpus as the receipt",
            fallback: "narrow the admissible declaration grammar and re-freeze via amendment",
            no_claim: "no claim about physical meaningfulness of admitted laws; typing \
                       only",
        },
        ClaimSpec {
            id: "i01-residual-jvp-vjp-actions",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Compiled residual, JVP, and VJP actions agree with independent \
                        dual-number forward mode and finite-difference probes on every \
                        pinned deck within the stated band, and JVP/VJP satisfy the \
                        transpose identity <v, Ju> = <J^T v, u>",
            hypotheses: &[
                "operands drawn from the pinned coupled MMS deck family",
                "probe step sizes follow the deck's declared curvature-aware schedule",
                "arithmetic is deterministic-mode f64",
            ],
            qoi: "max_relative_disagreement",
            unit: "1",
            tolerance: ToleranceSemantics::AbsRel {
                atol: 1e-12,
                rtol: 1e-9,
            },
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "fs-ad forward-mode duals plus central finite differences with \
                           independent step-size selection",
                independent: true,
                tcb_overlap: "shares fs-math deterministic scalar kernels",
            },
            activation: "admission leaf green at smoke tier",
            kill: "persistent transpose-identity violation beyond band on any pinned deck \
                   kills the fused-action design; receipt is the violating deck id and \
                   probe log",
            fallback: "unfused reference actions remain the shipping path",
            no_claim: "no claim for f32 or fast-mode arithmetic; no claim beyond the \
                       pinned deck family's regularity class",
        },
        ClaimSpec {
            id: "i01-block-lowering",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Block lowering preserves declared equation structure: variable \
                        blocking, sparsity pattern, units, and coupling topology of the \
                        lowered system match the typed source graph exactly",
            hypotheses: &[
                "source graphs from the pinned constitutive-graph and singular-block \
                 fixture families",
                "lowering runs with default (non-reordering) hints",
            ],
            qoi: "structure_preservation_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent structural differ comparing source-graph metadata \
                           against lowered-block metadata field by field",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "admission leaf green at smoke tier",
            kill: "any structure loss that the differ attributes to lowering (not to a \
                   declared reduction) kills silent-normalization behavior",
            fallback: "lowering emits an explicit reduction receipt for every transform",
            no_claim: "no claim about numerical conditioning of the lowered blocks",
        },
        ClaimSpec {
            id: "i01-conservation-exact-sequence",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Discrete conservation identities and the exact-sequence \
                        obligations (grad/curl/div compositions vanish) hold to \
                        round-off for every compiled operator family on the pinned decks",
            hypotheses: &[
                "operators compiled from the pinned coupled MMS deck family",
                "FEEC-compatible discretization as declared per deck",
            ],
            qoi: "max_identity_residual",
            unit: "1",
            tolerance: ToleranceSemantics::Absolute { atol: 1e-13 },
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "fs-feec incidence-operator composition checked in exact \
                           integer arithmetic where the complex permits, else \
                           compensated summation",
                independent: true,
                tcb_overlap: "shares fs-sparse assembly kernels",
            },
            activation: "lowering leaf green at smoke tier",
            kill: "a persistent nonzero composition on a conforming complex kills the \
                   operator-construction route (this is structural, not a tolerance \
                   issue)",
            fallback: "restrict to the operator families whose compositions vanish and \
                       amend the manifest to narrow scope",
            no_claim: "no claim on non-conforming or cut complexes; those belong to \
                       CutFEM lanes with their own manifests",
        },
        ClaimSpec {
            id: "i01-deterministic-cancellation-replay",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Compiled pipelines replay bit-identically across worker counts \
                        1/2/7 and across cancellation storms with request-drain-finalize \
                        semantics: every interrupted run resumes to the same bits as an \
                        uninterrupted run",
            hypotheses: &[
                "deterministic mode on; reductions use fixed-shape trees keyed by \
                 logical tile identity",
                "cancellation is injected at tile boundaries per the G4 schedule",
            ],
            qoi: "bitwise_replay_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G5,
            oracle: OracleRoute {
                identity: "byte comparison of content-addressed result artifacts across \
                           schedules (the artifact store, not the pipeline, decides \
                           equality)",
                independent: true,
                tcb_overlap: "shares fs-blake3 hashing",
            },
            activation: "residual-action leaf green at core tier",
            kill: "any schedule-dependent bit difference in deterministic mode is a \
                   release blocker, not a tolerance case",
            fallback: "none: determinism is load-bearing for every downstream lane",
            no_claim: "no claim across ISAs (cross-ISA replay is a separate two-host \
                       campaign); no claim for fast mode",
        },
        ClaimSpec {
            id: "i01-exact-discrete-adjoint-dwr",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Exact-discrete adjoints of compiled pipelines match transposed \
                        solves to round-off, and DWR error indicators built from them \
                        bound the true QoI error on the pinned decks within the stated \
                        effectivity band",
            hypotheses: &[
                "adjoints assembled by transposed solves, never by differentiating \
                 through Krylov iterations",
                "true QoI error measured against the deck's manufactured solution",
            ],
            qoi: "dwr_effectivity_index",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.5, hi: 2.0 },
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "manufactured-solution QoI error from the deck generator, \
                           evaluated without the production adjoint stack",
                independent: true,
                tcb_overlap: "shares fs-la dense kernels",
            },
            activation: "all five baseline claims green at core tier",
            kill: "effectivity outside band on a majority of decks kills DWR-driven \
                   adaptivity for I01 (indicators remain diagnostics)",
            fallback: "adjoint-free residual indicators, honestly labeled Estimated",
            no_claim: "no claim that DWR bounds hold off the pinned deck family's \
                       regularity class; indicators are Estimated evidence there",
        },
        ClaimSpec {
            id: "i01-rep-independent-coupled-action",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Coupled residual actions are representation-independent: \
                        semantically equivalent source graphs (pinned rewrite pairs) \
                        compile to actions that agree within the stated band on every \
                        probe state",
            hypotheses: &[
                "rewrite pairs generated with the fixture family's equivalence \
                 certificate",
                "probe states drawn from the pinned deck distribution",
            ],
            qoi: "max_pairwise_action_disagreement",
            unit: "1",
            tolerance: ToleranceSemantics::AbsRel {
                atol: 1e-12,
                rtol: 1e-9,
            },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "the fixture generator's equivalence certificate (constructed \
                           symbolically, independent of compilation)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "all five baseline claims green at core tier",
            kill: "a certified-equivalent pair with divergent actions beyond band kills \
                   representation-dependent fusion",
            fallback: "canonical-form pre-normalization before compilation",
            no_claim: "equivalence is claimed only for pairs carrying the generator's \
                       certificate, not for arbitrary user rewrites",
        },
        ClaimSpec {
            id: "i01-generated-law-completeness",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Refutation,
            statement: "Falsifier lane: search for a member of the declared multiphysics \
                        law basis that the generator cannot express or the compiler \
                        cannot admit; the moonshot completeness claim survives only \
                        while this search comes up empty",
            hypotheses: &[
                "the declared basis is the pinned generated-reference-law family",
                "search budget is the max-tier budget in explicits",
            ],
            qoi: "counterexample_count",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent basis enumerator over the pinned law-family \
                           grammar",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "maximal lane opens only after baseline closes at core tier",
            kill: "first admitted counterexample kills the completeness claim; the \
                   counterexample itself is retained as the receipt",
            fallback: "completeness restated over the surviving sub-basis via amendment",
            no_claim: "an empty search is Estimated evidence bounded by the search \
                       budget; it is never a completeness proof",
        },
        ClaimSpec {
            id: "i01-solver-hint-optimality",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "Compiler solver hints are within the stated factor of the best \
                        exhaustively-enumerated hint configuration under the certified \
                        cost budget on the pinned decks",
            hypotheses: &[
                "hint space is the finite enumeration frozen with the deck family",
                "cost is wall-clock under the max-tier budget on the named perf host \
                 class",
            ],
            qoi: "hint_cost_ratio_vs_best",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 1.0, hi: 1.25 },
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "exhaustive enumeration of the frozen hint space with \
                           fs-roofline-recorded runs",
                independent: true,
                tcb_overlap: "shares fs-roofline measurement plumbing",
            },
            activation: "maximal lane opens only after baseline closes at core tier",
            kill: "ratio beyond band on a majority of decks kills hint synthesis for \
                   I01 (hints demote to advisory)",
            fallback: "static default hints with recorded cost tables",
            no_claim: "optimality is relative to the frozen finite hint space and the \
                       named host class; no universal optimality claim",
        },
    ]
}

fn i01_fixtures() -> Vec<FixturePin> {
    vec![
        FixturePin {
            id: "i01-constitutive-graph-family",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: typed cross-domain constitutive graphs (thermal, \
elastic, electric, porous couplings) with per-node dimension \
signatures and port directions; sizes 2..=64 nodes; each graph \
emitted WITH a typing certificate constructed during generation.\n\
SEEDS: Philox stream 'i01/constitutive/<k>', k in 0..=4095.\n\
ILL-TYPED TWINS: each well-typed graph gets one mutated twin \
(dimension, variance, or port-direction fault) with the fault \
site recorded, for refusal-side testing.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i01-constitutive-graph-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "Same generator and certificate discipline as \
i01-constitutive-graph-family; SEEDS: Philox stream \
'i01/constitutive/<k>', k in 65536..=69631; withheld from \
development, opened only at adjudication.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i01-coupled-mms-deck-family",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: coupled field/port manufactured-solution decks on \
structured 1D/2D/3D complexes; smooth trigonometric/polynomial \
manufactured fields per physics; source terms derived \
symbolically at generation time and emitted with the deck; \
curvature-aware finite-difference probe schedules included.\n\
SEEDS: Philox stream 'i01/mms/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i01-coupled-mms-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "Same generator as i01-coupled-mms-deck-family; SEEDS: Philox \
stream 'i01/mms/<k>', k in 65536..=69631; withheld until \
adjudication.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i01-singular-block-patterns",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: block systems with engineered singular/deficient \
patterns (redundant constraints, zero pivots, decoupled islands, \
rank-deficient couplings) plus the exact expected structural \
diagnosis for each; sizes 2..=32 blocks.\n\
SEEDS: Philox stream 'i01/singular/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i01-generated-reference-laws",
            source: FixtureSource::AuthoredSpec {
                spec: "GRAMMAR: the declared multiphysics law basis (flux-gradient \
pairs, storage laws, source couplings) as a finite generative \
grammar with dimension constraints; enumeration order is the \
grammar's canonical derivation order so the basis is replayable.\n\
SEEDS: none (exhaustive enumeration to depth 4).",
            },
            partition: Partition::Development,
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn i01_obligations() -> Vec<ObligationRow> {
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
            leaf: "i01-law-admission",
            claims_covered: &["i01-typed-law-port-semantics"],
            unit_cases: UNIT_CASES,
            g0: "generators: constitutive-graph family incl. ill-typed twins; validity \
                 predicate: certificate agreement; laws: admission idempotence, refusal \
                 stability; shrinker: node-removal preserving fault site; replay seeds \
                 per explicits",
            decks: &["i01-constitutive-graph-family"],
            g3_relations: &[
                "node relabeling invariance",
                "port-order permutation invariance",
            ],
            g4_schedule: "cancel admission at every tile boundary class once; verify \
                          drain leaves no partial admission state",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i01_law_admission.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i01-law-admission slice)",
            obs_events: &[
                "admission.accepted",
                "admission.refused",
                "admission.cancelled",
            ],
            replay_command: "scripts/e2e/leapfrog/i01_law_admission.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i01-residual-derivative-actions",
            claims_covered: &["i01-residual-jvp-vjp-actions"],
            unit_cases: UNIT_CASES,
            g0: "generators: MMS decks; laws: linearity of JVP in tangent, transpose \
                 identity, zero-tangent annihilation; shrinker: deck coarsening; replay \
                 seeds per explicits",
            decks: &["i01-coupled-mms-deck-family"],
            g3_relations: &[
                "rigid transform invariance of residual norms",
                "unit-rescaling covariance of actions",
            ],
            g4_schedule: "cancel mid-action at tile boundaries; resume must reproduce \
                          uninterrupted bits",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i01_actions.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i01-actions slice)",
            obs_events: &[
                "action.residual",
                "action.jvp",
                "action.vjp",
                "action.cancelled",
            ],
            replay_command: "scripts/e2e/leapfrog/i01_actions.sh --replay <artifact-id>",
        },
        ObligationRow {
            // Leaf named after its entry point (i01_lowering.sh), NOT the
            // claim it covers: claim ids and leaf ids share one evidence-id
            // namespace under schema v2, so reusing the claim id here is a
            // DuplicateId freeze refusal.
            leaf: "i01-lowering",
            claims_covered: &["i01-block-lowering"],
            unit_cases: UNIT_CASES,
            g0: "generators: constitutive graphs + singular block patterns; laws: \
                 structure differ finds zero diffs on identity lowering; every engineered \
                 singularity gets its exact expected diagnosis; replay seeds per explicits",
            decks: &[
                "i01-constitutive-graph-family",
                "i01-singular-block-patterns",
            ],
            g3_relations: &["block permutation invariance of the diagnosis set"],
            g4_schedule: "cancel lowering between block passes; drain leaves no partial \
                          lowered artifact",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i01_lowering.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i01-lowering slice)",
            obs_events: &["lowering.block", "lowering.diagnosis", "lowering.cancelled"],
            replay_command: "scripts/e2e/leapfrog/i01_lowering.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i01-conservation-sequences",
            claims_covered: &["i01-conservation-exact-sequence"],
            unit_cases: UNIT_CASES,
            g0: "generators: MMS decks over conforming complexes; laws: dd=0 composition \
                 residuals, discrete conservation balances; replay seeds per explicits",
            decks: &["i01-coupled-mms-deck-family"],
            g3_relations: &["refinement monotonicity of identity residuals"],
            g4_schedule: "cancel assembly at tile boundaries; verify drained partial \
                          assembly is discarded whole",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i01_conservation.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i01-conservation slice)",
            obs_events: &[
                "assembly.operator",
                "identity.residual",
                "assembly.cancelled",
            ],
            replay_command: "scripts/e2e/leapfrog/i01_conservation.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i01-replay-cancellation",
            claims_covered: &["i01-deterministic-cancellation-replay"],
            unit_cases: UNIT_CASES,
            g0: "generators: MMS decks; laws: artifact-digest equality across schedules; \
                 replay seeds per explicits",
            decks: &["i01-coupled-mms-deck-family"],
            g3_relations: &["worker-count invariance is itself the metamorphic relation"],
            g4_schedule: "storm: cancel/resume at every boundary class, randomized \
                          (seeded) injection order, 128 storms per deck",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay; \
                        cross-ISA explicitly OUT of scope for v1",
            entry_point: "scripts/e2e/leapfrog/i01_replay.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i01-replay slice)",
            obs_events: &["replay.artifact", "storm.injected", "storm.drained"],
            replay_command: "scripts/e2e/leapfrog/i01_replay.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i01-maximal-lanes",
            claims_covered: &[
                "i01-exact-discrete-adjoint-dwr",
                "i01-rep-independent-coupled-action",
                "i01-generated-law-completeness",
                "i01-solver-hint-optimality",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: held-out families + generated reference laws + frozen hint \
                 enumeration; laws: adjoint transpose identity, equivalence-certificate \
                 checks, grammar enumeration replay; replay seeds per explicits",
            decks: &[
                "i01-coupled-mms-holdout",
                "i01-constitutive-graph-holdout",
                "i01-generated-reference-laws",
                "i01-industrial-coupled-benchmark",
            ],
            g3_relations: &[
                "certified rewrite-pair equivalence",
                "held-out vs development distribution parity checks",
            ],
            g4_schedule: "cancel each maximal lane once mid-campaign; a resumed lane must \
                          report identical verdicts",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i01_maximal.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (i01-maximal slice)",
            obs_events: &[
                "adjoint.dwr",
                "equivalence.pair",
                "completeness.counterexample",
                "hints.enumeration",
            ],
            replay_command: "scripts/e2e/leapfrog/i01_maximal.sh --replay <artifact-id>",
        },
    ]
}

fn i01_waivers() -> Vec<Waiver> {
    vec![Waiver {
        subject: "i01-industrial-coupled-benchmark",
        reason: "the external industrial coupled benchmark deck is not yet licensed or \
                 content-addressed; only its slot is preregistered",
        owner: "I01 implementation beads (frankensim-leapfrog-2026-program-i94v.1.1.x)",
        predicate: "deck bytes pinned with exact edition, license, digest, QoIs, and \
                    acceptance envelopes through the fs-vvreg registry discipline",
        expiry: "before the first Max-tier campaign run; review at each Phase-2 close \
                 burst",
        promotion_effect: "maximal claims may run without this deck but cannot cite it; \
                           any result depending on it stays Estimated until the pin \
                           lands",
    }]
}
