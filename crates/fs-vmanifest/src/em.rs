//! The EM (electromechanical exchange portfolio) VerificationManifest
//! draft (bead frankensim-leapfrog-2026-program-i94v.2.5.1).
//!
//! The portfolio-level gate ABOVE the field/harness/material/AP242/
//! manufacturing initiatives: it freezes the cross-boundary conventions
//! that keep terminal, winding, dielectric, harness, geometry,
//! PMI/GD&T, material/process/lot, loss, and evidence meaning intact
//! across every physics, optimization, exchange, and supplier boundary
//! — typed exchange-deck registries, one pinned semantic-loss taxonomy,
//! content-addressed terminal/net/winding identities, AP242
//! profile-contract deferral, and the hard separation between supported
//! core continuity and maximal topology/full-rung continuity. A weaker
//! receipt closes its own element and is never relabeled as the
//! stronger theorem; version-1 prose cards mint no proof authority.

use crate::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin, FixtureSource,
    GauntletTier, ManifestDraft, ObligationRow, OracleRoute, Partition, ToleranceSemantics, Waiver,
};

/// Build the EM draft. Consumers freeze it themselves; the conformance
/// battery proves it freezes.
#[must_use]
pub fn em_draft() -> ManifestDraft {
    ManifestDraft {
        initiative: "EM",
        title: "Electromechanical exchange portfolio gate: typed exchange-deck \
                registries, one pinned semantic-loss taxonomy across every physics/ \
                optimization/exchange/supplier boundary, content-addressed \
                terminal/net/winding identities, AP242 profile-contract deferral, \
                and a hard core/maximal claim lattice separation",
        version: 1,
        explicits: FiveExplicits {
            units: "SI base units throughout; field and loss channels carry \
                    deck-declared units; relative discrepancies dimensionless \
                    (unit '1'); exact bitwise/boolean verdicts use unit 'bit'; \
                    geometry lengths in metres with deck-declared tolerancing \
                    contexts",
            seeds: "Philox 4x32-10 counter streams keyed 'em/<fixture-id>/<case-index>'; \
                    development indices 0..=16383; core held-out indices 65536..=81919; \
                    maximal held-out indices 131072..=147455 (disjoint by construction, \
                    split frozen here; per-holdout subranges pinned in the fixture \
                    specs)",
            budgets: "smoke tier <= 60 s on one host; core tier <= 30 min; max tier <= 8 h \
                      on a quiet perf host; <= 16 GiB memory per lane; composition \
                      budgets are per-deck declarations whose exhaustion is a typed \
                      resumable outcome; accuracy budgets are the per-claim tolerance \
                      fields",
            versions: "fs-vmanifest schema v2; toolchain pinned by \
                       rust-toolchain.toml; sibling pins by constellation.lock; \
                       AP242 profile/edition pins and initiative receipt versions \
                       live in the taxonomy map",
            capabilities: "no network; no FFI; deterministic mode mandatory for every G5 \
                           row; frontier/moonshot lanes stay behind feature flags; \
                           field/circuit solver-receipt and full-rung dependencies \
                           are waivered until their beads land",
        },
        claims: em_claims(),
        fixtures: em_fixtures(),
        obligations: em_obligations(),
        waivers: em_waivers(),
        amendment_rules: "Any change is a successor version through FrozenManifest::amend; \
                          the amendment record names every invalidated claim/obligation \
                          descendant; an amendment after campaign start invalidates the \
                          affected evidence, which must be re-earned; there is no in-place \
                          edit path in the type system",
    }
}

#[allow(clippy::too_many_lines)]
fn em_claims() -> Vec<ClaimSpec> {
    vec![
        ClaimSpec {
            id: "em-typed-exchange-decks",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Every motor-drive and harness deck in the portfolio \
                        registries names its AP242 profile/edition, \
                        terminal/net/winding map, material/process/lot identity, \
                        and fidelity region; the registries round-trip bit-stably \
                        and refuse every enumerated malformed deck with a typed \
                        diagnostic",
            hypotheses: &[
                "decks drawn from the pinned motor-drive and harness registries",
                "malformed decks are the enumerated ones: unnamed profile/edition, \
                 unmapped terminal or net, undeclared material/process/lot, \
                 missing fidelity region",
            ],
            qoi: "registry_roundtrip_and_refusal_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent decoder plus the registry generator's own \
                           validity certificate (emitted during generation, not by \
                           the portfolio exchange path)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "first implementation bead of the EM registry leaf opens",
            kill: "any accepted malformed deck or refused valid deck returns the \
                   registry to design review with the disagreeing deck as receipt",
            fallback: "explicit per-deck normal forms with per-field validation",
            no_claim: "no claim that a registered deck is manufacturable or \
                       EMC-compliant; registration is typing, not review",
        },
        ClaimSpec {
            id: "em-semantic-loss-taxonomy",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Every physics, optimization, exchange, and supplier \
                        boundary crossing is classified by the pinned semantic-loss \
                        taxonomy as lossless, declared-loss under a named taxon, or \
                        refusal; any boundary crossing with unclassified semantic \
                        loss refuses instead of silently normalizing",
            hypotheses: &[
                "taxonomy from the pinned taxonomy map",
                "lane-local loss classifications declare their embedding into the \
                 portfolio taxonomy",
            ],
            qoi: "loss_classification_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent boundary walker replaying crossings against \
                           the declared taxonomy (separate from the exchange code)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "registry leaf green at smoke tier",
            kill: "any silent semantic loss at a boundary is Sev-0 for the \
                   portfolio lane",
            fallback: "boundary crossings refuse with the missing taxon named",
            no_claim: "no claim that the taxonomy is exhaustive over future \
                       boundaries; only that it is single and enforced over the \
                       pinned ones",
        },
        ClaimSpec {
            id: "em-terminal-net-identity",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Terminal, net, and winding identities are content-addressed \
                        from declared canonical invariants and exchange-route \
                        independent: the same net through different exchange routes \
                        carries the same identity, distinct windings always carry \
                        distinct identities, and identity collisions or merges are \
                        typed refusals, never silent",
            hypotheses: &[
                "nets and windings from the pinned identity holdout with \
                 ground-truth identity maps",
                "identity inputs are the declared canonical connectivity \
                 invariants, not serialization artifacts",
            ],
            qoi: "identity_stability_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "ground-truth identity maps carried by the fixture \
                           generator",
                independent: true,
                tcb_overlap: "shares the hash primitive (fs-blake3) only",
            },
            activation: "registry leaf green at smoke tier",
            kill: "any route-dependent identity or silent merge blocks the \
                   identity lane",
            fallback: "identities carry route provenance suffixes and say so",
            no_claim: "identity stability is relative to the declared canonical \
                       invariants; a wrong invariant choice is a design defect, not \
                       an identity defect",
        },
        ClaimSpec {
            id: "em-profile-edition-deferral",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Portfolio exchange decks defer geometry, PMI, and GD&T \
                        equivalence semantics to the pinned AP242 profile contract: \
                        every deck names its profile/edition and equivalence \
                        relation set, and no portfolio lane re-implements geometry \
                        equivalence outside that contract",
            hypotheses: &[
                "exchange decks from the pinned registries",
                "the AP242 profile contract is the semantic authority for \
                 geometry/PMI/GD&T equivalence",
            ],
            qoi: "deferral_completeness_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent scanner for geometry-equivalence \
                           implementations outside the declared contract references",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "registry leaf green at smoke tier",
            kill: "any shadow equivalence implementation found by the scanner \
                   blocks the exchange decks",
            fallback: "decks without named profiles/editions refuse at \
                       registration",
            no_claim: "no claim about the AP242 profile contract's own \
                       correctness; deferral structure only",
        },
        ClaimSpec {
            id: "em-core-maximal-separation",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Supported-core continuity (core-G4-class) and maximal \
                        topology/full-rung continuity (maximal-G7-class) claims are \
                        disjoint lattice elements: no maximal evidence requirement \
                        gates core promotion, no core receipt closes a maximal \
                        element, and the separation is enforced by the typed claim \
                        registry, not by convention",
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
            activation: "registry leaf green at smoke tier",
            kill: "any cross-lattice gating found by the auditor blocks portfolio \
                   promotion",
            fallback: "none: separation is structural, not degradable",
            no_claim: "separation does not rank the lattices; a maximal theorem can \
                       still be discovered before a core benchmark closes",
        },
        ClaimSpec {
            id: "em-visual-match-falsifier",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Refutation,
            statement: "Adversarial look-alike families (visually identical CAD \
                        with swapped PMI datums, attractive field plots from \
                        perturbed windings, renamed nets on identical geometry, \
                        lot-swapped materials) attempt to make semantically \
                        distinct decks pass continuity or collapse distinct \
                        evidence states; any success refutes the \
                        semantic-continuity discipline",
            hypotheses: &[
                "fixtures constructed with ground-truth continuity/distinctness \
                 proofs in-spec",
                "the continuity checks run at production settings",
            ],
            qoi: "false_continuity_count",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "hand-proved continuity/distinctness constructions \
                           carried by the fixture spec",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "runs in every campaign tier from the first continuity \
                         commit",
            kill: "this lane is never killed; it is the standing tripwire",
            fallback: "none: a nonzero count is a release blocker",
            no_claim: "absence of refutation is not completeness and cannot prove \
                       continuity soundness; the lane only falsifies",
        },
        ClaimSpec {
            id: "em-loss-attribution-reproducibility",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Field, circuit, and loss attributions are reproducible \
                        across initiatives from the same solver receipts: two lanes \
                        consuming one receipt produce identical per-channel loss \
                        attributions, and every attribution carries the receipt \
                        provenance it was computed from",
            hypotheses: &[
                "solver receipts from the pinned shared family",
                "attribution closure per the loss-convention semantics \
                 (whole-machine loss never re-derived portfolio-side)",
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
            activation: "identity and lattice leaves green at core tier",
            kill: "any receipt-identical attribution divergence between lanes kills \
                   the assurance lane",
            fallback: "per-initiative attributions labeled with their lane \
                       identity, no portfolio aggregation",
            no_claim: "attribution reproducibility does not validate field physics \
                       and does not grant regulatory EMC approval, supplier \
                       qualification, or manufacturing certification; \
                       whole-machine qualification is the waivered lane",
        },
        ClaimSpec {
            id: "em-objective-sensitivity-consistency",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Topology-optimization objectives and sensitivities shared \
                        across the portfolio agree with the reference adjoint \
                        lane's records on shared decks within the declared band",
            hypotheses: &[
                "shared decks from the pinned registries",
                "the reference adjoint lane's records are the reference route",
            ],
            qoi: "relative_sensitivity_discrepancy",
            unit: "1",
            tolerance: ToleranceSemantics::Relative { rtol: 1e-6 },
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "reference exact-discrete adjoint records (an independent \
                           lane of this portfolio, not the portfolio aggregation \
                           path)",
                independent: true,
                tcb_overlap: "shares deck definitions only",
            },
            activation: "identity leaf green at core tier",
            kill: "gate failure on any shared deck kills portfolio sensitivity \
                   aggregation until root-caused",
            fallback: "per-lane sensitivities only, no portfolio-level objective \
                       composition",
            no_claim: "consistency is against the reference records, not against \
                       ground truth; a shared systematic error passes this gate",
        },
        ClaimSpec {
            id: "em-full-rung-continuity-composition",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "Maximal full-rung continuity composes every boundary's \
                        maximal G7 receipt (physics -> optimization -> exchange -> \
                        supplier -> manufacturing) into one composed lineage whose \
                        record names every constituent receipt identity, and \
                        composition gaps are localized Unknown windows with budget \
                        receipts, never silent joins",
            hypotheses: &[
                "composition cases from the pinned maximal holdout",
                "constituent receipt identities per the [S] identity lane",
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
            activation: "behind the em-moonshot feature flag on a pre-proof \
                         successor version; frontier lanes green",
            kill: "any silent join on a ground-truth composition fixture kills the \
                   lane",
            fallback: "compositions refuse; constituent receipts reported \
                       individually",
            no_claim: "composition lineage is structural; no claim that composed \
                       rungs are physically or contractually valid; version-1 \
                       prose mints no composition authority",
        },
        ClaimSpec {
            id: "em-maximal-theorem-scope",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "The maximal (maximal-G7-class) theorem-scope card carries \
                        machine-checkable obligations for every global claim \
                        (topology-continuity completeness, EMC-margin \
                        qualification, full-rung equivalence), and no global \
                        theorem language enters a portfolio report without its \
                        obligation reference",
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
            activation: "behind the em-moonshot feature flag on a pre-proof \
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

const POLICY_SPEC: &str = "EM_PORTFOLIO_POLICY_V1
DECK_REGISTRY= every deck names AP242 profile/edition, terminal/net/winding map, \
material/process/lot identity, and fidelity region; registration is typing, never \
manufacturability or EMC review
LOSS_TAXONOMY= one pinned semantic-loss taxonomy over every physics/optimization/ \
exchange/supplier boundary; crossings are lossless, declared-loss, or refusal; \
silent semantic loss is Sev-0
IDENTITY= terminal/net/winding identities are content-addressed from declared \
canonical invariants; route-dependent identities and silent merges refuse
PROFILE_SEMANTICS= exchange decks defer to the pinned AP242 profile contract; \
shadow geometry-equivalence implementations are forbidden; a visual CAD match is \
not semantic continuity
LATTICE= supported-core continuity (core-G4-class) and maximal topology/full-rung \
continuity (maximal-G7-class) claims are disjoint lattice elements; no \
cross-lattice gating in either direction
ASSURANCE_METRICS= attributions reproducible from shared receipts; whole-machine \
losses never re-derived portfolio-side; adversarial look-alike families run in \
every tier
THEOREM_AUTHORITY= version 1 has prose cards only and mints no proof; moonshot \
lanes activate only on a pre-proof successor version with machine-readable \
obligations; no regulatory EMC approval, supplier qualification, or manufacturing \
certification is granted at any version
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

fn em_fixtures() -> Vec<FixturePin> {
    vec![
        FixturePin {
            id: "em-portfolio-policy-v1",
            source: FixtureSource::AuthoredSpec { spec: POLICY_SPEC },
            partition: Partition::Development,
        },
        FixturePin {
            id: "em-motor-drive-deck-registry",
            source: FixtureSource::AuthoredSpec {
                spec: "em fixture v1: motor-drive deck registry — electrostatic and \
                       topology-optimized machine variants naming terminal/net/ \
                       winding maps, AP242 profile/edition, material/process/lot \
                       identities, and fidelity regions; generator emits validity \
                       certificates; development indices 0..=16383; seeds key \
                       'em/motor/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "em-harness-deck-registry",
            source: FixtureSource::AuthoredSpec {
                spec: "em fixture v1: harness deck registry — cable harness \
                       routings, shield terminations, connector nets, and EMC \
                       fixture variants naming material/process/lot identities and \
                       fidelity regions; development indices 0..=16383; seeds key \
                       'em/harness/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "em-taxonomy-map",
            source: FixtureSource::AuthoredSpec {
                spec: "em fixture v1: the pinned semantic-loss taxonomy and profile \
                       contract map — taxon definitions per boundary class, AP242 \
                       profile/edition pins, geometry/PMI/GD&T equivalence relation \
                       sets, initiative receipt version pins, and the declared \
                       embeddings of every lane-local classification; development \
                       indices 0..=16383; seeds key 'em/taxonomy/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "em-shared-solver-receipt-family",
            source: FixtureSource::AuthoredSpec {
                spec: "em fixture v1: shared solver receipt family — field/circuit \
                       receipts with known loss-attribution structures consumed by \
                       multiple lanes for reproducibility checks; development \
                       indices 0..=16383; seeds key 'em/receipt/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "em-net-identity-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "em fixture v1 HOLDOUT: terminal/net/winding identity \
                       adversaries — the same net through different exchange \
                       routes, near-identical windings, renamed nets on identical \
                       geometry, with ground-truth identity maps in-spec; core \
                       held-out indices 65536..=69631; one EM.G3 consumer \
                       (em-identity-lanes); seeds key 'em/identity/<k>'",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "em-lookalike-adversarial-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "em fixture v1 HOLDOUT: hand-constructed look-alike families \
                       — visually identical CAD with swapped PMI datums, perturbed \
                       windings behind attractive field plots, lot-swapped \
                       materials, with continuity/distinctness proofs in-spec; \
                       core held-out indices 69632..=73727; one EM.G3 consumer \
                       (em-assurance-metrics); seeds key 'em/lookalike/<k>'",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "em-fullrung-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "em fixture v1 HOLDOUT: full-rung composition cases — \
                       composed physics/optimization/exchange/supplier/ \
                       manufacturing lineages with ground-truth composition maps \
                       and theorem-scope cards; maximal held-out indices \
                       131072..=135167; one EM.G3 consumer \
                       (em-moonshot-composition); seeds key 'em/fullrung/<k>'",
            },
            partition: Partition::HeldOut,
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn em_obligations() -> Vec<ObligationRow> {
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
            leaf: "em-deck-admission",
            claims_covered: &["em-semantic-loss-taxonomy", "em-typed-exchange-decks"],
            unit_cases: UNIT_CASES,
            g0: "generators: deck registries incl. enumerated malformed classes; \
                 validity predicate: generator certificate agreement; laws: \
                 canonical round-trip bit-stability, refusal totality, taxonomy \
                 classification totality over pinned boundaries; shrinker: \
                 deck-row removal preserving the violated rule; replay seeds per \
                 explicits",
            decks: &[
                "em-harness-deck-registry",
                "em-motor-drive-deck-registry",
                "em-portfolio-policy-v1",
                "em-taxonomy-map",
            ],
            g3_relations: &[
                "deck relabeling invariance",
                "registry-order invariance of admission verdicts",
            ],
            g4_schedule: "request->drain->finalize injection between decode, \
                          taxonomy walk, and registration; checkpoint at each \
                          phase boundary; drained cancellation leaves no partially \
                          registered deck",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/em_admission.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (em-admission slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "deck.admitted",
                "deck.refused",
            ],
            replay_command: "scripts/e2e/leapfrog/em_admission.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "em-identity-lanes",
            claims_covered: &["em-terminal-net-identity"],
            unit_cases: UNIT_CASES,
            g0: "generators: net/winding families incl. the identity holdout's \
                 adversaries; laws: route-independence of identities, distinctness \
                 of distinct windings, typed collision refusals; replay seeds per \
                 explicits",
            decks: &[
                "em-motor-drive-deck-registry",
                "em-net-identity-core-holdout",
                "em-portfolio-policy-v1",
            ],
            g3_relations: &[
                "exchange-route invariance of identities",
                "net relabeling covariance of identity verdicts",
            ],
            g4_schedule: "request->drain->finalize injection mid-derivation; \
                          checkpoint at per-net boundaries; drained cancellation \
                          reports derived identities only",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        identity derivations",
            entry_point: "scripts/e2e/leapfrog/em_identity.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (em-identity slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "identity.collision",
                "identity.derived",
            ],
            replay_command: "scripts/e2e/leapfrog/em_identity.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "em-profile-semantics",
            claims_covered: &["em-profile-edition-deferral"],
            unit_cases: UNIT_CASES,
            g0: "generators: exchange decks with and without named \
                 profiles/editions; laws: deferral completeness (scanner finds no \
                 shadow equivalence implementations), refusal of unnamed profiles, \
                 equivalence-relation-set presence; replay seeds per explicits",
            decks: &[
                "em-harness-deck-registry",
                "em-portfolio-policy-v1",
                "em-taxonomy-map",
            ],
            g3_relations: &["profile renaming covariance of deferral verdicts"],
            g4_schedule: "request->drain->finalize injection mid-scan; checkpoint \
                          at per-deck boundaries; drained cancellation reports \
                          scanned decks only",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        deferral verdicts",
            entry_point: "scripts/e2e/leapfrog/em_profiles.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (em-profiles slice)",
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
            replay_command: "scripts/e2e/leapfrog/em_profiles.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "em-lattice-separation",
            claims_covered: &["em-core-maximal-separation"],
            unit_cases: UNIT_CASES,
            g0: "generators: claim registries with declared lattice assignments plus \
                 adversarial cross-gating attempts; laws: auditor recomputation \
                 agreement, cross-lattice gating refusal in both directions; replay \
                 seeds per explicits",
            decks: &["em-motor-drive-deck-registry", "em-portfolio-policy-v1"],
            g3_relations: &["claim relabeling invariance of lattice audits"],
            g4_schedule: "request->drain->finalize injection mid-audit; checkpoint \
                          at per-claim boundaries; drained cancellation reports \
                          audited claims only",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        audit verdicts",
            entry_point: "scripts/e2e/leapfrog/em_lattice.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (em-lattice slice)",
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
            replay_command: "scripts/e2e/leapfrog/em_lattice.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "em-assurance-metrics",
            claims_covered: &[
                "em-loss-attribution-reproducibility",
                "em-objective-sensitivity-consistency",
                "em-visual-match-falsifier",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: shared solver receipts and the look-alike adversarial \
                 holdout; laws: attribution reproducibility across lanes, \
                 sensitivity agreement within rtol, zero false continuities; \
                 replay seeds per explicits",
            decks: &[
                "em-lookalike-adversarial-core-holdout",
                "em-portfolio-policy-v1",
                "em-shared-solver-receipt-family",
            ],
            g3_relations: &[
                "receipt-order invariance of attribution verdicts",
                "rendering invariance of continuity verdicts",
            ],
            g4_schedule: "request->drain->finalize injection with deliberate \
                          receipt-budget exhaustion (typed resumable outcomes); \
                          checkpoint between receipt batches; cancel \
                          mid-attribution with rollback verified",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        attributions, sensitivities, and falsifier verdicts",
            entry_point: "scripts/e2e/leapfrog/em_assurance.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (em-assurance slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "attribution.compared",
                "continuity.verdict",
                "sensitivity.compared",
            ],
            replay_command: "scripts/e2e/leapfrog/em_assurance.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "em-moonshot-composition",
            claims_covered: &[
                "em-full-rung-continuity-composition",
                "em-maximal-theorem-scope",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: composition cases and theorem cards from the maximal \
                 holdout; laws: composed lineage completeness against ground-truth \
                 maps, localized Unknown windows with budget receipts, no \
                 unreferenced theorem language; replay seeds per explicits",
            decks: &[
                "em-fullrung-max-holdout",
                "em-harness-deck-registry",
                "em-portfolio-policy-v1",
            ],
            g3_relations: &["composition-order invariance of lineage verdicts"],
            g4_schedule: "the core of this lane IS G4: request->drain->finalize \
                          injection with budget exhaustion, timeout, and \
                          cancellation inside composition and scope scanning; \
                          checkpoint at constituent boundaries; each outcome drained \
                          and typed; BudgetExhausted stays Unknown",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        composed receipts and scope verdicts",
            entry_point: "scripts/e2e/leapfrog/em_moonshot.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (em-moonshot slice)",
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
            replay_command: "scripts/e2e/leapfrog/em_moonshot.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
    ]
}

fn em_waivers() -> Vec<Waiver> {
    vec![
        Waiver {
            subject: "em-assurance-metrics",
            reason: "loss-attribution reproducibility and sensitivity consistency \
                     consume the field/circuit initiatives' core G4 solver \
                     receipts, none of which have landed; the initiative lanes are \
                     preregistered but unimplemented",
            owner: "portfolio-em lane owner",
            predicate: "the field/circuit initiative [S] obligations are green at \
                        core tier AND their core G4 receipts are published at the \
                        versions pinned in the taxonomy map",
            expiry: "first EM campaign review after the initiative lanes land; \
                     re-justify or retire at every manifest amendment",
            promotion_effect: "the [F] claims can close only on \
                               receipt-reproducibility semantics; whole-machine \
                               EMC-margin qualification stays Unknown while the \
                               waiver is live",
        },
        Waiver {
            subject: "em-moonshot-composition",
            reason: "full-rung continuity requires every boundary's maximal G7 \
                     receipt and the portfolio [S]/[F] lanes green at core tier; \
                     no initiative has implemented its maximal lane yet",
            owner: "portfolio-em lane owner",
            predicate: "boundary maximal G7 receipts published at the pinned \
                        versions AND the EM [S]/[F] obligations are green at core \
                        tier",
            expiry: "first EM campaign review after the maximal receipts land; \
                     re-justify or retire at every manifest amendment",
            promotion_effect: "the [M] claims stay Unknown and cannot close; [S]/[F] \
                               promotion is unaffected because their obligations run \
                               on registry and receipt fixtures only",
        },
    ]
}
