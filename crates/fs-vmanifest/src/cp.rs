//! The CP (compiler portfolio) VerificationManifest draft
//! (bead frankensim-leapfrog-2026-program-i94v.1.7.1).
//!
//! The portfolio-level gate ABOVE the six compiler initiatives: it
//! freezes the cross-initiative conventions that make law compilation,
//! causalization, hybrid execution, conservation attribution,
//! obligation graphs, and evidence planning composable — typed seam
//! schemas between every stage, one content-addressed semantic-source
//! identity chain, exact-version receipt admission, ModeLedger
//! deferral for hybrid semantics, and the hard separation between core
//! (initiative-G4-class) and maximal (initiative-G7-class) claim
//! lattices. A weaker receipt closes its own element and is never
//! relabeled as the stronger theorem; version-1 prose cards mint no
//! proof authority.

use crate::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin, FixtureSource,
    GauntletTier, ManifestDraft, ObligationRow, OracleRoute, Partition, ToleranceSemantics, Waiver,
};

/// Build the CP draft. Consumers freeze it themselves; the conformance
/// battery proves it freezes.
#[must_use]
pub fn cp_draft() -> ManifestDraft {
    ManifestDraft {
        initiative: "CP",
        title: "Compiler portfolio gate: typed seam schemas across the six-stage \
                law-to-runtime path, one content-addressed semantic-source identity, \
                exact-version receipt admission, ModeLedger hybrid deferral, and a \
                hard core/maximal claim lattice separation",
        version: 1,
        explicits: FiveExplicits {
            units: "SI base units throughout; solve-plan costs and residual \
                    attributions carry deck-declared units; relative discrepancies \
                    dimensionless (unit '1'); exact bitwise/boolean verdicts use \
                    unit 'bit'; seam payload units travel with every schema row",
            seeds: "Philox 4x32-10 counter streams keyed 'cp/<fixture-id>/<case-index>'; \
                    development indices 0..=16383; core held-out indices 65536..=81919; \
                    maximal held-out indices 131072..=147455 (disjoint by construction, \
                    split frozen here; per-holdout subranges pinned in the fixture \
                    specs)",
            budgets: "smoke tier <= 60 s on one host; core tier <= 30 min; max tier <= 8 h \
                      on a quiet perf host; <= 16 GiB memory per lane; composition \
                      budgets are per-machine declarations whose exhaustion is a typed \
                      resumable outcome; accuracy budgets are the per-claim tolerance \
                      fields",
            versions: "fs-vmanifest schema v2; toolchain pinned by \
                       rust-toolchain.toml; sibling pins by constellation.lock; \
                       initiative receipt versions pinned in the seam-schema map",
            capabilities: "no network; no FFI; deterministic mode mandatory for every G5 \
                           row; frontier/moonshot lanes stay behind feature flags; \
                           initiative core-G4 and maximal-G7 receipt dependencies are \
                           waivered until their beads land",
        },
        claims: cp_claims(),
        fixtures: cp_fixtures(),
        obligations: cp_obligations(),
        waivers: cp_waivers(),
        amendment_rules: "Any change is a successor version through FrozenManifest::amend; \
                          the amendment record names every invalidated claim/obligation \
                          descendant; an amendment after campaign start invalidates the \
                          affected evidence, which must be re-earned; there is no in-place \
                          edit path in the type system",
    }
}

#[allow(clippy::too_many_lines)]
fn cp_claims() -> Vec<ClaimSpec> {
    vec![
        ClaimSpec {
            id: "cp-typed-seam-schemas",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Every seam on the six-stage path (law corpus -> causalized \
                        equations -> solve plan -> hybrid execution -> conservation \
                        attribution -> obligation/evidence rows) carries a typed, \
                        versioned schema; seam payloads round-trip bit-stably and \
                        every enumerated malformed payload refuses with a typed \
                        diagnostic naming the violated schema field",
            hypotheses: &[
                "seam payloads drawn from the pinned law-corpus and machine \
                 registries through the pinned seam-schema map",
                "malformed payloads are the enumerated ones: unversioned schema, \
                 unmapped stage pair, dangling identity reference, undeclared \
                 regime or rank/index",
            ],
            qoi: "seam_roundtrip_and_refusal_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent decoder plus the schema generator's own \
                           validity certificate (emitted during generation, not by \
                           the portfolio compiler)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "first implementation bead of the CP seam leaf opens",
            kill: "any accepted malformed payload or refused valid payload returns \
                   the schema map to design review with the disagreeing payload as \
                   receipt",
            fallback: "explicit per-seam normal forms with per-field validation",
            no_claim: "no claim that a schema-valid payload is semantically \
                       adequate; schema admission is typing, not review",
        },
        ClaimSpec {
            id: "cp-semantic-source-identity",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "One semantic source governs the whole path: every stage \
                        output carries the content-addressed identity of the law \
                        corpus it was compiled from, cross-stage comparison happens \
                        only through the declared identity chain, and any \
                        hand-maintained seam drift refuses instead of silently \
                        normalizing",
            hypotheses: &[
                "identity chain from the pinned seam-schema map",
                "stage-local identities declare their embedding into the portfolio \
                 identity chain",
            ],
            qoi: "source_identity_enforcement_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent identity walker replaying stage outputs \
                           against the declared chain (separate from the compiler \
                           path)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "seam leaf green at smoke tier",
            kill: "any silent acceptance of seam drift is Sev-0 for the portfolio \
                   lane",
            fallback: "cross-stage comparisons refuse with the missing identity \
                       edge named",
            no_claim: "no claim that the compiled semantics are physically valid; \
                       only that they descend from one pinned source",
        },
        ClaimSpec {
            id: "cp-receipt-version-admission",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Initiative receipts are content-addressed and admitted only \
                        at their exact pinned versions: the same receipt consumed \
                        through different lanes carries the same identity, distinct \
                        receipt versions always carry distinct identities, and \
                        version mismatches, stale receipts, or identity collisions \
                        are typed refusals, never silent",
            hypotheses: &[
                "receipts from the pinned receipt-version holdout with ground-truth \
                 identity maps",
                "identity inputs are the declared canonical receipt fields, not \
                 serialization artifacts",
            ],
            qoi: "receipt_admission_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "ground-truth identity maps carried by the fixture \
                           generator",
                independent: true,
                tcb_overlap: "shares the hash primitive (fs-blake3) only",
            },
            activation: "seam leaf green at smoke tier",
            kill: "any lane-dependent receipt identity or silent stale admission \
                   blocks the receipt lane",
            fallback: "receipts carry lane provenance suffixes and say so",
            no_claim: "admission stability is relative to the declared canonical \
                       fields; a wrong field choice is a design defect, not an \
                       admission defect",
        },
        ClaimSpec {
            id: "cp-mode-ledger-deferral",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Portfolio hybrid machines defer mode/event semantics to the \
                        ModeLedger contract: every hybrid machine names its ledger \
                        and mode budget, and no compiler lane re-implements mode \
                        switching outside that contract",
            hypotheses: &[
                "hybrid machines from the pinned machine registry",
                "the ModeLedger contract is the semantic authority for \
                 modes/events/resets",
            ],
            qoi: "deferral_completeness_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent scanner for mode-switching implementations \
                           outside the declared contract references",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "seam leaf green at smoke tier",
            kill: "any shadow mode implementation found by the scanner blocks the \
                   hybrid machines",
            fallback: "hybrid machines without named ledgers refuse at admission",
            no_claim: "no claim about the ModeLedger's own correctness; deferral \
                       structure only",
        },
        ClaimSpec {
            id: "cp-core-maximal-separation",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Core (initiative-G4-class) and maximal (initiative-G7-class) \
                        claims are disjoint lattice elements: no maximal evidence \
                        requirement gates core promotion, no core receipt closes a \
                        maximal element, and the separation is enforced by the typed \
                        claim registry, not by convention",
            hypotheses: &[
                "claim registries from the pinned families",
                "lattice element assignments are declarations checked at freeze",
            ],
            qoi: "lattice_separation_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent lattice auditor recomputing promotion gates \
                           from the registry declarations",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "seam leaf green at smoke tier",
            kill: "any cross-lattice gating found by the auditor blocks portfolio \
                   promotion",
            fallback: "none: separation is structural, not degradable",
            no_claim: "separation does not rank the lattices; a maximal theorem can \
                       still be discovered before a core benchmark closes",
        },
        ClaimSpec {
            id: "cp-threshold-edit-falsifier",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Refutation,
            statement: "Adversarial post-result mutation families (edited acceptance \
                        arithmetic, relaxed tolerances, swapped receipt versions, \
                        back-dated seeds, re-serialized manifests) attempt to make a \
                        frozen manifest accept evidence it preregistered as failing; \
                        any success refutes the no-post-result-editing discipline",
            hypotheses: &[
                "fixtures constructed with ground-truth accept/reject proofs \
                 in-spec",
                "the freeze/amend path runs at production settings",
            ],
            qoi: "accepted_mutation_count",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "hand-proved accept/reject constructions carried by the \
                           fixture spec",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "runs in every campaign tier from the first portfolio \
                         commit",
            kill: "this lane is never killed; it is the standing tripwire",
            fallback: "none: a nonzero count is a release blocker",
            no_claim: "absence of refutation is not completeness and cannot prove \
                       manifest immutability; the lane only falsifies",
        },
        ClaimSpec {
            id: "cp-attribution-reproducibility",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Balance-microscope conservation attributions are \
                        reproducible across initiatives from the same execution \
                        receipts: two lanes consuming one receipt produce identical \
                        per-seam residual attributions, and every attribution carries \
                        the receipt provenance it was computed from",
            hypotheses: &[
                "execution receipts from the pinned shared family",
                "attribution closure per the balance-microscope semantics \
                 (whole-machine residual never re-derived portfolio-side)",
            ],
            qoi: "attribution_reproducibility_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "independent attribution recomputation from the raw \
                           receipts (separate implementation from the portfolio \
                           path)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "receipt and lattice leaves green at core tier",
            kill: "any receipt-identical attribution divergence between lanes kills \
                   the assurance lane",
            fallback: "per-initiative attributions labeled with their lane \
                       identity, no portfolio aggregation",
            no_claim: "attribution reproducibility does not validate physical laws \
                       or conservation itself; whole-machine qualification is the \
                       waivered lane",
        },
        ClaimSpec {
            id: "cp-solve-plan-consistency",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Solve-plan objectives and costs shared across the portfolio \
                        agree with the reference causalization lane's records on \
                        shared machine decks within the declared band",
            hypotheses: &[
                "shared machine decks from the pinned registries",
                "the reference causalization lane's records are the reference \
                 route",
            ],
            qoi: "relative_plan_cost_discrepancy",
            unit: "1",
            tolerance: ToleranceSemantics::Relative { rtol: 1e-6 },
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "reference causalization records (an independent lane of \
                           this portfolio, not the portfolio aggregation path)",
                independent: true,
                tcb_overlap: "shares machine deck definitions only",
            },
            activation: "receipt leaf green at core tier",
            kill: "gate failure on any shared deck kills portfolio plan-cost \
                   aggregation until root-caused",
            fallback: "per-lane plan costs only, no portfolio-level objective \
                       composition",
            no_claim: "consistency is against the reference records, not against \
                       ground truth; a shared systematic error passes this gate",
        },
        ClaimSpec {
            id: "cp-law-to-runtime-composition",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "The strongest law-to-runtime path composes every \
                        initiative's maximal G7 receipt into one composed lineage \
                        whose record names every constituent receipt identity, and \
                        composition gaps are localized Unknown windows with budget \
                        receipts, never silent joins",
            hypotheses: &[
                "composition cases from the pinned maximal holdout",
                "constituent receipt identities per the [S] receipt-admission lane",
            ],
            qoi: "composition_lineage_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent lineage walker over composed receipts against \
                           ground-truth composition maps",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "behind the cp-moonshot feature flag on a pre-proof \
                         successor version; frontier lanes green",
            kill: "any silent join on a ground-truth composition fixture kills the \
                   lane",
            fallback: "compositions refuse; constituent receipts reported \
                       individually",
            no_claim: "composition lineage is structural; no claim that composed \
                       runtimes are physically valid; version-1 prose mints no \
                       composition authority",
        },
        ClaimSpec {
            id: "cp-maximal-theorem-scope",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "The maximal (initiative-G7-class) theorem-scope card carries \
                        machine-checkable obligations for every global claim \
                        (whole-machine semantic equivalence, hybrid qualification, \
                        cross-compiler assurance), and no global theorem language \
                        enters a portfolio report without its obligation reference",
            hypotheses: &[
                "theorem cards from the pinned maximal holdout",
                "obligation references are typed registry entries, not prose",
            ],
            qoi: "theorem_scope_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent report scanner for theorem language without \
                           obligation references",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "behind the cp-moonshot feature flag on a pre-proof \
                         successor version; requires the composition lane's waiver \
                         status reviewed",
            kill: "any unreferenced theorem language in a portfolio report is Sev-0 \
                   for this lane",
            fallback: "reports carry lane-local claims only",
            no_claim: "the card constrains language, it proves nothing; version-1 \
                       prose mints no theorem authority",
        },
    ]
}

const POLICY_SPEC: &str = "CP_PORTFOLIO_POLICY_V1
SEAM_SCHEMAS= every stage seam carries a typed, versioned schema; schema admission \
is typing, never semantic review
SOURCE_IDENTITY= one content-addressed semantic-source identity chain; cross-stage \
comparison only through the declared chain; silent seam drift is Sev-0
RECEIPT_VERSIONS= initiative receipts admitted only at their exact pinned versions; \
lane-dependent identities and silent stale admissions refuse
MODE_SEMANTICS= hybrid machines defer to the ModeLedger contract; shadow mode \
implementations are forbidden
LATTICE= core (initiative-G4-class) and maximal (initiative-G7-class) claims are \
disjoint lattice elements; no cross-lattice gating in either direction
ASSURANCE_METRICS= attributions reproducible from shared receipts; whole-machine \
residuals never re-derived portfolio-side; adversarial mutation families run in \
every tier
THEOREM_AUTHORITY= version 1 has prose cards only and mints no proof; moonshot lanes \
activate only on a pre-proof successor version with machine-readable obligations
EVIDENCE_STATES= Verified/Validated/Estimated/Failed/Refuted/Unknown/no-claim are \
distinct; one axis never substitutes for another; a weaker receipt closes its own \
lattice element only
HOLDOUT= held-out fixtures adjudicate; development never touches held-out indices; \
each holdout has exactly one stage-local consumer row
LIFECYCLE= request->drain->finalize with checkpoint boundaries; cancellation is \
drained, never dropped; partial success cannot publish normal authority
LOGGING= structured fs-obs events per obligation row incl. the six lifecycle kinds; \
budgets, seeds, versions, capabilities logged at run start
RETENTION= receipts, failure bundles, and adjudication records are content-addressed \
and retained for replay; raw counters are diagnostics only
FAILURE_BUNDLE= every red or Unknown outcome retains a bounded reproducible bundle \
naming fixture, seed, budget, and disposition
PROMOTION= claims promote only through their preregistered obligation rows on frozen \
manifests; amendment invalidates affected descendants which must re-earn evidence
LEAF_REQUIREMENT= every execution leaf maps to exactly one obligation row; there are \
no unnamed skips";

fn cp_fixtures() -> Vec<FixturePin> {
    vec![
        FixturePin {
            id: "cp-portfolio-policy-v1",
            source: FixtureSource::AuthoredSpec { spec: POLICY_SPEC },
            partition: Partition::Development,
        },
        FixturePin {
            id: "cp-law-corpus-registry",
            source: FixtureSource::AuthoredSpec {
                spec: "cp fixture v1: law corpus registry for the representative \
                       electromechanical-thermal-fluid hybrid machine — lumped \
                       electromechanical laws, thermal network branches, \
                       incompressible flow branches, each row naming regime \
                       validity, rank/index declaration, and identity-chain \
                       embedding; generator emits validity certificates; \
                       development indices 0..=16383; seeds key 'cp/laws/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "cp-machine-registry",
            source: FixtureSource::AuthoredSpec {
                spec: "cp fixture v1: machine registry — motor-drive plus \
                       heat-exchanger plus pump hybrid machine variants naming \
                       their ModeLedgers and mode budgets; development indices \
                       0..=16383; seeds key 'cp/machine/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "cp-seam-schema-map",
            source: FixtureSource::AuthoredSpec {
                spec: "cp fixture v1: the versioned seam-schema map for the \
                       six-stage path — per-seam schema versions, the semantic \
                       source identity chain, pinned initiative receipt versions, \
                       and the declared embeddings of every stage-local schema; \
                       development indices 0..=16383; seeds key 'cp/seam/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "cp-shared-receipt-family",
            source: FixtureSource::AuthoredSpec {
                spec: "cp fixture v1: shared execution receipt family — receipts \
                       with known attribution structures consumed by multiple lanes \
                       for reproducibility checks; development indices 0..=16383; \
                       seeds key 'cp/receipt/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "cp-receipt-version-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "cp fixture v1 HOLDOUT: receipt-version adversaries — the \
                       same receipt through different lanes, near-identical stale \
                       versions, re-serialized receipts, with ground-truth identity \
                       maps in-spec; core held-out indices 65536..=69631; one CP.G3 \
                       consumer (cp-receipt-lanes); seeds key 'cp/receipts/<k>'",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "cp-threshold-adversarial-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "cp fixture v1 HOLDOUT: hand-constructed post-result mutation \
                       families — edited acceptance arithmetic, relaxed tolerances, \
                       swapped receipt versions, back-dated seeds, with \
                       accept/reject proofs in-spec; core held-out indices \
                       69632..=73727; one CP.G3 consumer (cp-assurance-metrics); \
                       seeds key 'cp/threshold/<k>'",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "cp-composition-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "cp fixture v1 HOLDOUT: law-to-runtime composition cases — \
                       composed maximal-receipt lineages with ground-truth \
                       composition maps and theorem-scope cards; maximal held-out \
                       indices 131072..=135167; one CP.G3 consumer \
                       (cp-moonshot-composition); seeds key 'cp/composition/<k>'",
            },
            partition: Partition::HeldOut,
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn cp_obligations() -> Vec<ObligationRow> {
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
            leaf: "cp-stage-admission",
            claims_covered: &["cp-semantic-source-identity", "cp-typed-seam-schemas"],
            unit_cases: UNIT_CASES,
            g0: "generators: seam payloads incl. enumerated malformed classes; \
                 validity predicate: generator certificate agreement; laws: \
                 canonical round-trip bit-stability, refusal totality, identity \
                 chain enforcement; shrinker: payload-field removal preserving the \
                 violated rule; replay seeds per explicits",
            decks: &[
                "cp-law-corpus-registry",
                "cp-machine-registry",
                "cp-portfolio-policy-v1",
                "cp-seam-schema-map",
            ],
            g3_relations: &[
                "payload relabeling invariance",
                "registry-order invariance of admission verdicts",
            ],
            g4_schedule: "request->drain->finalize injection between decode, \
                          identity walk, and admission; checkpoint at each phase \
                          boundary; drained cancellation leaves no partially \
                          admitted payload",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/cp_admission.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (cp-admission slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "seam.admitted",
                "seam.refused",
            ],
            replay_command: "scripts/e2e/leapfrog/cp_admission.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "cp-receipt-lanes",
            claims_covered: &["cp-receipt-version-admission"],
            unit_cases: UNIT_CASES,
            g0: "generators: receipt families incl. the receipt-version holdout's \
                 adversaries; laws: lane-independence of receipt identities, \
                 distinctness of distinct versions, typed stale/collision refusals; \
                 replay seeds per explicits",
            decks: &[
                "cp-law-corpus-registry",
                "cp-portfolio-policy-v1",
                "cp-receipt-version-core-holdout",
            ],
            g3_relations: &[
                "receipt re-serialization invariance of identities",
                "lane-route invariance of identities",
            ],
            g4_schedule: "request->drain->finalize injection mid-derivation; \
                          checkpoint at per-receipt boundaries; drained cancellation \
                          reports derived identities only",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        identity derivations",
            entry_point: "scripts/e2e/leapfrog/cp_receipts.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (cp-receipts slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "receipt.admitted",
                "receipt.refused",
            ],
            replay_command: "scripts/e2e/leapfrog/cp_receipts.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "cp-mode-semantics",
            claims_covered: &["cp-mode-ledger-deferral"],
            unit_cases: UNIT_CASES,
            g0: "generators: hybrid machines with and without named ledgers; laws: \
                 deferral completeness (scanner finds no shadow implementations), \
                 refusal of unnamed ledgers, mode-budget presence; replay seeds per \
                 explicits",
            decks: &["cp-machine-registry", "cp-portfolio-policy-v1"],
            g3_relations: &["ledger renaming covariance of deferral verdicts"],
            g4_schedule: "request->drain->finalize injection mid-scan; checkpoint at \
                          per-machine boundaries; drained cancellation reports \
                          scanned machines only",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        deferral verdicts",
            entry_point: "scripts/e2e/leapfrog/cp_modes.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (cp-modes slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "deferral.verified",
                "shadow.detected",
            ],
            replay_command: "scripts/e2e/leapfrog/cp_modes.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "cp-lattice-separation",
            claims_covered: &["cp-core-maximal-separation"],
            unit_cases: UNIT_CASES,
            g0: "generators: claim registries with declared lattice assignments plus \
                 adversarial cross-gating attempts; laws: auditor recomputation \
                 agreement, cross-lattice gating refusal in both directions; replay \
                 seeds per explicits",
            decks: &["cp-law-corpus-registry", "cp-portfolio-policy-v1"],
            g3_relations: &["claim relabeling invariance of lattice audits"],
            g4_schedule: "request->drain->finalize injection mid-audit; checkpoint \
                          at per-claim boundaries; drained cancellation reports \
                          audited claims only",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        audit verdicts",
            entry_point: "scripts/e2e/leapfrog/cp_lattice.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (cp-lattice slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "crossgate.refused",
                "lattice.audited",
            ],
            replay_command: "scripts/e2e/leapfrog/cp_lattice.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "cp-assurance-metrics",
            claims_covered: &[
                "cp-attribution-reproducibility",
                "cp-solve-plan-consistency",
                "cp-threshold-edit-falsifier",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: shared execution receipts and the threshold adversarial \
                 holdout; laws: attribution reproducibility across lanes, plan-cost \
                 agreement within rtol, zero accepted post-result mutations; replay \
                 seeds per explicits",
            decks: &[
                "cp-portfolio-policy-v1",
                "cp-shared-receipt-family",
                "cp-threshold-adversarial-core-holdout",
            ],
            g3_relations: &[
                "receipt-order invariance of attribution verdicts",
                "seam relabeling invariance of mutation verdicts",
            ],
            g4_schedule: "request->drain->finalize injection with deliberate \
                          receipt-budget exhaustion (typed resumable outcomes); \
                          checkpoint between receipt batches; cancel mid-attribution \
                          with rollback verified",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        attributions, plan costs, and falsifier verdicts",
            entry_point: "scripts/e2e/leapfrog/cp_assurance.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (cp-assurance slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "attribution.compared",
                "mutation.verdict",
                "plancost.compared",
            ],
            replay_command: "scripts/e2e/leapfrog/cp_assurance.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "cp-moonshot-composition",
            claims_covered: &["cp-law-to-runtime-composition", "cp-maximal-theorem-scope"],
            unit_cases: UNIT_CASES,
            g0: "generators: composition cases and theorem cards from the maximal \
                 holdout; laws: composed lineage completeness against ground-truth \
                 maps, localized Unknown windows with budget receipts, no \
                 unreferenced theorem language; replay seeds per explicits",
            decks: &[
                "cp-composition-max-holdout",
                "cp-machine-registry",
                "cp-portfolio-policy-v1",
            ],
            g3_relations: &["composition-order invariance of lineage verdicts"],
            g4_schedule: "the core of this lane IS G4: request->drain->finalize \
                          injection with budget exhaustion, timeout, and \
                          cancellation inside composition and scope scanning; \
                          checkpoint at constituent boundaries; each outcome drained \
                          and typed; BudgetExhausted stays Unknown",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        composed receipts and scope verdicts",
            entry_point: "scripts/e2e/leapfrog/cp_moonshot.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (cp-moonshot slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "budget.receipt",
                "composition.lineage",
                "scope.verdict",
            ],
            replay_command: "scripts/e2e/leapfrog/cp_moonshot.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
    ]
}

fn cp_waivers() -> Vec<Waiver> {
    vec![
        Waiver {
            subject: "cp-assurance-metrics",
            reason: "attribution reproducibility and solve-plan consistency consume \
                     the six compiler initiatives' core G4 execution receipts, none \
                     of which have landed; the initiative lanes are preregistered \
                     but unimplemented",
            owner: "portfolio-cp lane owner",
            predicate: "all six compiler-initiative [S] obligations are green at \
                        core tier AND their core G4 receipts are published at the \
                        versions pinned in the seam-schema map",
            expiry: "first CP campaign review after the initiative lanes land; \
                     re-justify or retire at every manifest amendment",
            promotion_effect: "the [F] claims can close only on \
                               receipt-reproducibility semantics; whole-machine \
                               attribution qualification stays Unknown while the \
                               waiver is live",
        },
        Waiver {
            subject: "cp-moonshot-composition",
            reason: "law-to-runtime composition requires every initiative's maximal \
                     G7 receipt and the portfolio [S]/[F] lanes green at core tier; \
                     no initiative has implemented its maximal lane yet",
            owner: "portfolio-cp lane owner",
            predicate: "initiative maximal G7 receipts published at the pinned \
                        versions AND the CP [S]/[F] obligations are green at core \
                        tier",
            expiry: "first CP campaign review after the maximal receipts land; \
                     re-justify or retire at every manifest amendment",
            promotion_effect: "the [M] claims stay Unknown and cannot close; [S]/[F] \
                               promotion is unaffected because their obligations run \
                               on registry and receipt fixtures only",
        },
    ]
}
