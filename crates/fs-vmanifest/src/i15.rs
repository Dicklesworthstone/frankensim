//! The I15 (executable standards compiler) VerificationManifest draft
//! (bead frankensim-leapfrog-2026-program-i94v.1.6.8.1).
//!
//! Baseline lattice ([S]): edition-aware clause identity, typed normative
//! rules and human judgments, applicability, obligation execution, scoped
//! conformance receipts with semantic edition diffs. Maximal lattice
//! ([F]/[M]): proof-carrying cross-standard composition,
//! jurisdiction/edition-aware theorem obligations, a transitive-impact
//! completeness falsifier lane, and machine-checkable assurance that
//! structurally cannot erase human authority. A weaker receipt closes its
//! own element and is never relabeled as the stronger theorem.
//!
//! LICENSING BOUNDARY: every fixture in this manifest is a SYNTHETIC
//! standard-shaped pack authored here; no licensed standard text is
//! embedded anywhere. The real editions enter only through the waived
//! external deck slot, pinned later under the fs-vvreg registry
//! discipline.

use crate::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin, FixtureSource,
    GauntletTier, ManifestDraft, ObligationRow, OracleRoute, Partition, ToleranceSemantics, Waiver,
};

/// Build the I15 draft. Consumers freeze it themselves; the conformance
/// battery proves it freezes.
#[must_use]
pub fn i15_draft() -> ManifestDraft {
    ManifestDraft {
        initiative: "I15",
        title: "Executable standards compiler gate: clause identity, typed rules and \
                 judgments, applicability, obligation execution, and scoped \
                 conformance receipts",
        version: 1,
        explicits: FiveExplicits {
            units: "SI base units for physical quantities in worked examples; clause \
                    and edition identities are content-addressed strings; exact \
                    verdicts use unit 'bit'",
            seeds: "Philox 4x32-10 counter streams keyed 'i15/<fixture-id>/<case-index>'; \
                    development cases use seed indices 0..=4095; held-out cases use \
                    65536..=69631 (disjoint by construction, split frozen here)",
            budgets: "smoke tier <= 60 s on one host; core tier <= 30 min; max tier <= 8 h \
                      on a quiet perf host; <= 16 GiB memory per lane; accuracy budgets \
                      are the per-claim tolerance fields",
            versions: "fs-vmanifest schema v1; manifest version 1; toolchain pinned by \
                       rust-toolchain.toml; sibling pins by constellation.lock",
            capabilities: "no network; no FFI; no licensed standard text in any fixture; \
                           deterministic mode mandatory for every G5 row; \
                           frontier/moonshot lanes stay behind feature flags",
        },
        claims: i15_claims(),
        fixtures: i15_fixtures(),
        obligations: i15_obligations(),
        waivers: i15_waivers(),
        amendment_rules: "Any change is a successor version through FrozenManifest::amend; \
                          the amendment record names every invalidated claim/obligation \
                          descendant; an amendment after campaign start invalidates the \
                          affected evidence, which must be re-earned; there is no in-place \
                          edit path in the type system",
    }
}

#[allow(clippy::too_many_lines)]
fn i15_claims() -> Vec<ClaimSpec> {
    vec![
        ClaimSpec {
            id: "i15-edition-aware-clause-identity",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Clause identities bind (standard, edition, clause path, content \
                        digest) inseparably: the same clause text in two editions gets \
                        two identities, any content change moves the identity, and a \
                        bare clause number without an edition is refused as an \
                        identity",
            hypotheses: &[
                "packs from the pinned synthetic standard-shaped families",
                "identity hashing is the workspace domain-separated BLAKE3 discipline",
            ],
            qoi: "identity_law_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent identity recomputation from the pack generator's \
                           clause records (generated with, not by, the compiler)",
                independent: true,
                tcb_overlap: "shares fs-blake3 hashing",
            },
            activation: "first implementation bead of the I15 identity leaf opens",
            kill: "any identity collision across editions, or acceptance of an \
                   edition-less citation, returns the identity design to review",
            fallback: "none: edition-blind clause identity is the exact failure mode \
                       this initiative exists to close",
            no_claim: "identity only; no claim about the normative meaning of clauses",
        },
        ClaimSpec {
            id: "i15-typed-rules-and-judgments",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Normative rules compile to typed executable obligations with \
                        units and applicability conditions, while declared human \
                        judgments compile to typed AUTHORITY nodes that executable \
                        paths can cite but can never synthesize, satisfy, or overwrite; \
                        a machine path claiming judgment authority is refused at \
                        admission",
            hypotheses: &[
                "rule and judgment declarations from the pinned pack families with \
                 generator-recorded kind certificates",
            ],
            qoi: "rule_judgment_admission_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent kind-checker walking the generator's rule/\
                           judgment certificates",
                independent: true,
                tcb_overlap: "shares fs-qty dimension arithmetic only",
            },
            activation: "identity leaf green at smoke tier",
            kill: "any machine-synthesized node admitted with judgment authority is a \
                   release blocker (this is the human-authority boundary)",
            fallback: "narrow the rule vocabulary; judgments stay declaration-only",
            no_claim: "typing only; no claim that declared judgments are correct, only \
                       that their authority is preserved and attributed",
        },
        ClaimSpec {
            id: "i15-applicability",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Applicability predicates (scope, ratings, configuration \
                        conditions) evaluate against pinned case records to exactly \
                        the generator-recorded applicability verdicts, and \
                        out-of-scope cases refuse with the failing scope condition \
                        named rather than executing anyway",
            hypotheses: &[
                "cases from the pinned pack families with recorded applicability \
                 ground truth including boundary values",
            ],
            qoi: "applicability_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "the pack generator's applicability records",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "rules leaf green at smoke tier",
            kill: "any out-of-scope execution (running an obligation whose scope \
                   condition fails) kills silent-scope behavior",
            fallback: "conservative refusal on any unresolved scope condition",
            no_claim: "applicability against the synthetic packs' recorded scopes; no \
                       claim of legally correct scope interpretation",
        },
        ClaimSpec {
            id: "i15-obligation-execution",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Executing compiled obligations on the pinned worked-example \
                        cases reproduces the generator-computed results within band, \
                        with every intermediate carrying its declared unit and every \
                        clause citation resolving to an edition-aware identity",
            hypotheses: &[
                "worked examples are generator-authored (synthetic packs), computed \
                 symbolically at generation time",
                "arithmetic is deterministic-mode f64",
            ],
            qoi: "max_relative_result_deviation",
            unit: "1",
            tolerance: ToleranceSemantics::AbsRel {
                atol: 1e-12,
                rtol: 1e-9,
            },
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "the generator's symbolically computed worked-example \
                           results",
                independent: true,
                tcb_overlap: "shares fs-math deterministic scalar kernels",
            },
            activation: "applicability leaf green at smoke tier",
            kill: "deviation beyond band that traces to formula compilation (not deck \
                   error) returns the obligation compiler to review",
            fallback: "interpretation-mode execution (slow reference evaluator)",
            no_claim: "agreement with the synthetic packs only; agreement with real \
                       licensed editions is gated on the waived external decks",
        },
        ClaimSpec {
            id: "i15-conformance-receipts-edition-diffs",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Conformance receipts are scoped to (edition set, clause set, \
                        input digest, obligation results) and never claim beyond their \
                        scope; semantic edition diffs classify every clause pair across \
                        the pinned edition changes exactly as the generator recorded \
                        (unchanged, renumbered, semantically changed, added, removed)",
            hypotheses: &[
                "edition-change pairs from the pinned family with recorded diff \
                 classifications",
            ],
            qoi: "receipt_scope_and_diff_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "the generator's diff classification records plus an \
                           independent receipt-scope auditor",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "execution leaf green at core tier",
            kill: "a renumbered clause classified as semantically changed (or vice \
                   versa) on a majority of pairs kills the diff engine; any receipt \
                   citing an out-of-scope clause is a release blocker",
            fallback: "textual diffs honestly labeled non-semantic",
            no_claim: "diff classes are the generator's taxonomy; no claim of legal \
                       equivalence between editions",
        },
        ClaimSpec {
            id: "i15-proof-carrying-cross-standard-composition",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Obligations composed across two packs carry a proof object \
                        (clause identities, unit conversions, applicability joins) \
                        that an independent checker replays without the compiler; \
                        composition over the pinned conflicting-clause decks refuses \
                        with both clause identities named instead of silently picking \
                        one",
            hypotheses: &[
                "cross-pack cases from the pinned families with generator-recorded \
                 composition ground truth and conflict inventories",
            ],
            qoi: "composition_proof_replay_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent proof-object replayer (checker without the \
                           composition engine)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "all five baseline claims green at core tier",
            kill: "a composition whose proof fails replay, or a silent conflict \
                   resolution, kills cross-standard composition (compositions demote \
                   to manual with named conflicts)",
            fallback: "single-standard receipts with manual composition notes",
            no_claim: "composition correctness is proof replay over synthetic packs; \
                       no claim that real standards compose legally",
        },
        ClaimSpec {
            id: "i15-jurisdiction-edition-theorem-obligations",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Theorem obligations parameterized by (jurisdiction, edition) \
                        resolve to exactly the deck's recorded obligation matrix: no \
                        cross-edition leakage, no jurisdiction bleed, and unresolved \
                        (jurisdiction, edition) pairs refuse rather than defaulting",
            hypotheses: &[
                "jurisdiction/edition matrices from the pinned families with recorded \
                 resolution ground truth",
            ],
            qoi: "obligation_matrix_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "the generator's obligation matrix records",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "all five baseline claims green at core tier",
            kill: "any cross-edition leakage (an obligation from edition A satisfied \
                   by edition B's rule) kills edition-generic resolution",
            fallback: "single-(jurisdiction, edition) operation with explicit locks",
            no_claim: "resolution against the synthetic matrices; no claim of legal \
                       correctness of any jurisdiction mapping",
        },
        ClaimSpec {
            id: "i15-transitive-impact-completeness",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Refutation,
            statement: "Falsifier lane: search the pinned edition-change decks \
                        (including held-out) for an edition change whose transitive \
                        impact set, as computed by the compiler, misses an affected \
                        clause, obligation, or previously issued receipt that the \
                        generator's ground truth names; the completeness claim \
                        survives only while this search comes up empty",
            hypotheses: &[
                "ground-truth impact sets are generator-recorded over the pinned \
                 dependency graphs",
                "search budget is the max-tier budget in explicits",
            ],
            qoi: "counterexample_count",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "the generators' recorded transitive impact sets",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "maximal lane opens only after baseline closes at core tier",
            kill: "first verified missed impact kills the completeness claim; the \
                   missed clause/receipt is retained as the receipt",
            fallback: "impact sets restated as lower bounds via amendment",
            no_claim: "an empty search is Estimated evidence bounded by the pinned \
                       graphs and budget; it is never a completeness proof",
        },
        ClaimSpec {
            id: "i15-assurance-preserves-human-authority",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "Every machine-checkable assurance chain that transits a \
                        declared human judgment retains that judgment's authority \
                        reference end-to-end: removing the judgment node structurally \
                        invalidates the chain, and no sequence of machine steps can \
                        reconstruct an equivalent chain without it",
            hypotheses: &[
                "assurance chains from the pinned declared-judgment decks with \
                 recorded judgment inventories",
                "'cannot reconstruct' is checked by exhaustive machine-step closure \
                 over the deck's finite step vocabulary",
            ],
            qoi: "authority_preservation_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent chain auditor plus machine-step closure \
                           enumerator over the deck's step vocabulary",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "maximal lane opens only after baseline closes at core tier",
            kill: "any machine-only reconstruction of a judgment-dependent chain is a \
                   release blocker for the assurance engine (this is the erasure the \
                   claim forbids)",
            fallback: "assurance chains refuse to transit judgments at all (manual \
                       handoff points)",
            no_claim: "non-reconstructibility is relative to the deck's finite step \
                       vocabulary; no claim against arbitrary future machine steps",
        },
    ]
}

fn i15_fixtures() -> Vec<FixturePin> {
    vec![
        FixturePin {
            id: "i15-gear-rating-pack",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: a synthetic gear-rating standard-shaped pack \
(clause tree, typed rating rules with units, scope conditions, \
symbolically computed worked examples, declared human-judgment \
slots for application factors); structure inspired by \
ISO-6336-class rating standards but containing NO licensed text.\n\
SEEDS: Philox stream 'i15/gear/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i15-motor-method-pack",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: a synthetic motor test-method pack (efficiency \
determination methods, rating clauses, measurement-uncertainty \
rules, judgment slots for method selection); IEC-60034-class \
structure, NO licensed text.\n\
SEEDS: Philox stream 'i15/motor/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i15-emc-vv-pack",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: a synthetic EMC/V&V pack (test levels, conformance \
criteria, verification/validation split rules) with cross-pack \
references into the gear and motor packs for composition cases; \
NO licensed text.\n\
SEEDS: Philox stream 'i15/emcvv/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i15-edition-change-pairs",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: edition pairs of the synthetic packs with recorded \
per-clause diff classes (unchanged, renumbered, semantically \
changed, added, removed) and recorded transitive impact sets \
over the packs' dependency graphs.\n\
SEEDS: Philox stream 'i15/editions/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i15-edition-change-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "Same generator as i15-edition-change-pairs; SEEDS: Philox \
stream 'i15/editions/<k>', k in 65536..=69631; withheld until \
adjudication.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i15-conflicting-clause-scenarios",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: cross-pack cases engineered so two clauses impose \
incompatible obligations on the same quantity, with the conflict \
inventory recorded (composition must refuse, naming both \
identities).\n\
SEEDS: Philox stream 'i15/conflicts/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i15-declared-judgment-chains",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: assurance chains transiting declared human-judgment \
nodes, with judgment inventories, the finite machine-step \
vocabulary, and recorded (jurisdiction, edition) obligation \
matrices for the resolution claims.\n\
SEEDS: Philox stream 'i15/judgments/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i15-declared-judgment-chains-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "Same generator as i15-declared-judgment-chains; SEEDS: Philox \
stream 'i15/judgments/<k>', k in 65536..=69631; withheld until \
adjudication.",
            },
            partition: Partition::HeldOut,
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn i15_obligations() -> Vec<ObligationRow> {
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
            leaf: "i15-clause-identity",
            claims_covered: &["i15-edition-aware-clause-identity"],
            unit_cases: UNIT_CASES,
            g0: "generators: all pack families; laws: identity binds edition, content \
                 change moves identity, edition-less citation refused; shrinker: \
                 clause-tree pruning preserving the identity fault; replay seeds per \
                 explicits",
            decks: &["i15-gear-rating-pack", "i15-motor-method-pack"],
            g3_relations: &["clause renumbering changes identity, never aliases it"],
            g4_schedule: "cancel identity assembly at boundary classes; drain leaves \
                          no partial identity table",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i15_identity.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i15-identity slice)",
            obs_events: &["clause.identified", "clause.refused", "identity.cancelled"],
            replay_command: "scripts/e2e/leapfrog/i15_identity.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i15-rules-judgments",
            claims_covered: &["i15-typed-rules-and-judgments"],
            unit_cases: UNIT_CASES,
            g0: "generators: pack families incl. machine-claims-judgment twins; laws: \
                 kind-certificate agreement, judgment-synthesis refusal; replay seeds \
                 per explicits",
            decks: &["i15-gear-rating-pack", "i15-motor-method-pack"],
            g3_relations: &["rule relabeling invariance of kinds"],
            g4_schedule: "cancel rule compilation mid-pack; drain leaves no partial \
                          obligation set",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i15_rules.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i15-rules slice)",
            obs_events: &["rule.compiled", "judgment.declared", "synthesis.refused"],
            replay_command: "scripts/e2e/leapfrog/i15_rules.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i15-applicability",
            claims_covered: &["i15-applicability"],
            unit_cases: UNIT_CASES,
            g0: "generators: pack families with applicability ground truth incl. \
                 boundary values; laws: verdict agreement, named-condition refusals; \
                 replay seeds per explicits",
            decks: &[
                "i15-gear-rating-pack",
                "i15-motor-method-pack",
                "i15-emc-vv-pack",
            ],
            g3_relations: &["scope-preserving case permutations leave verdicts fixed"],
            g4_schedule: "cancel mid-evaluation; resumed evaluation reproduces the \
                          verdict",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i15_applicability.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i15-applicability slice)",
            obs_events: &["scope.evaluated", "scope.refused", "scope.cancelled"],
            replay_command: "scripts/e2e/leapfrog/i15_applicability.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i15-obligation-execution",
            claims_covered: &["i15-obligation-execution"],
            unit_cases: UNIT_CASES,
            g0: "generators: worked-example cases; laws: result agreement within band, \
                 every intermediate unit-checked, every citation edition-aware; replay \
                 seeds per explicits",
            decks: &[
                "i15-gear-rating-pack",
                "i15-motor-method-pack",
                "i15-emc-vv-pack",
            ],
            g3_relations: &["unit-rescaling covariance of worked-example results"],
            g4_schedule: "cancel mid-execution; drained runs publish no partial \
                          results",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i15_execution.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i15-execution slice)",
            obs_events: &[
                "obligation.executed",
                "result.checked",
                "execution.cancelled",
            ],
            replay_command: "scripts/e2e/leapfrog/i15_execution.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i15-receipts-diffs",
            claims_covered: &["i15-conformance-receipts-edition-diffs"],
            unit_cases: UNIT_CASES,
            g0: "generators: edition-change pairs; laws: receipt scope audit, exact \
                 diff classification agreement; replay seeds per explicits",
            decks: &["i15-edition-change-pairs"],
            g3_relations: &["diff(A,B) and diff(B,A) are consistent inverses"],
            g4_schedule: "cancel mid-diff; drained diffs publish no partial \
                          classification",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i15_receipts.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i15-receipts slice)",
            obs_events: &["receipt.issued", "diff.classified", "diff.cancelled"],
            replay_command: "scripts/e2e/leapfrog/i15_receipts.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i15-maximal-lanes",
            claims_covered: &[
                "i15-proof-carrying-cross-standard-composition",
                "i15-jurisdiction-edition-theorem-obligations",
                "i15-transitive-impact-completeness",
                "i15-assurance-preserves-human-authority",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: conflict scenarios, judgment chains, edition holdouts; \
                 laws: proof replay, obligation-matrix agreement, impact-set sweep \
                 bookkeeping, machine-step closure enumeration; replay seeds per \
                 explicits",
            decks: &[
                "i15-conflicting-clause-scenarios",
                "i15-declared-judgment-chains",
                "i15-declared-judgment-chains-holdout",
                "i15-edition-change-holdout",
                "i15-licensed-standards-editions",
            ],
            g3_relations: &[
                "held-out vs development distribution parity checks",
                "composition order invariance up to recorded conflict refusals",
            ],
            g4_schedule: "cancel each maximal lane once mid-campaign; a resumed lane \
                          must report identical verdicts",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i15_maximal.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (i15-maximal slice)",
            obs_events: &[
                "composition.proof",
                "matrix.resolved",
                "impact.counterexample",
                "authority.audited",
            ],
            replay_command: "scripts/e2e/leapfrog/i15_maximal.sh --replay <artifact-id>",
        },
    ]
}

fn i15_waivers() -> Vec<Waiver> {
    vec![Waiver {
        subject: "i15-licensed-standards-editions",
        reason: "real licensed standard editions (ISO 6336-class, IEC 60034-class, \
                 NAFEMS-class documents) cannot be embedded as fixtures; only their \
                 slot is preregistered, and all in-repo packs are synthetic \
                 standard-shaped stand-ins",
        owner: "I15 implementation beads (frankensim-leapfrog-2026-program-i94v.1.6.x)",
        predicate: "edition artifacts pinned with exact edition, license, digest, \
                    QoIs, and acceptance envelopes through the fs-vvreg registry \
                    discipline (see the g2-* unpinned targets there)",
        expiry: "before the first Max-tier campaign run; review at each Phase-2 close \
                 burst",
        promotion_effect: "claims about REAL standards stay Estimated and cannot cite \
                           these decks until the pins land; synthetic-pack results \
                           never relabel as real-edition conformance",
    }]
}
