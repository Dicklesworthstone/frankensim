//! The PD (periodic-dynamics portfolio) VerificationManifest draft
//! (bead frankensim-leapfrog-2026-program-i94v.4.2.1).
//!
//! The portfolio-level gate ABOVE the I09 initiative manifest: it
//! freezes the cross-initiative conventions that make orbit, Floquet,
//! adjoint, continuation, and machine-campaign evidence composable —
//! deck registries with declared orbit classes, one portfolio-wide
//! phase/gauge/symmetry convention, content-addressed section/branch
//! identities, deferred event/reset semantics, and the hard separation
//! between core (I09.G4) and maximal (I09.G7) claim lattices. A weaker
//! receipt closes its own element and is never relabeled as the
//! stronger theorem; version-1 prose cards mint no proof authority.

use crate::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin, FixtureSource,
    GauntletTier, ManifestDraft, ObligationRow, OracleRoute, Partition, ToleranceSemantics, Waiver,
};

/// Build the PD draft. Consumers freeze it themselves; the conformance
/// battery proves it freezes.
#[must_use]
pub fn pd_draft() -> ManifestDraft {
    ManifestDraft {
        initiative: "PD",
        title: "Periodic-dynamics portfolio gate: composable deck registries, one \
                phase/gauge/symmetry convention, content-addressed section and branch \
                identities, deferred hybrid semantics, and a hard core/maximal claim \
                lattice separation",
        version: 1,
        explicits: FiveExplicits {
            units: "SI base units throughout; periods in seconds, multipliers and \
                    quotient metrics dimensionless (unit '1'); exact bitwise/boolean \
                    verdicts use unit 'bit'; deck-declared channel units travel with \
                    every deck row",
            seeds: "Philox 4x32-10 counter streams keyed 'pd/<fixture-id>/<case-index>'; \
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
                       rust-toolchain.toml; sibling pins by constellation.lock",
            capabilities: "no network; no FFI; deterministic mode mandatory for every G5 \
                           row; frontier/moonshot lanes stay behind feature flags; \
                           eigensolver-proof and initiative-lane dependencies are \
                           waivered until their beads land",
        },
        claims: pd_claims(),
        fixtures: pd_fixtures(),
        obligations: pd_obligations(),
        waivers: pd_waivers(),
        amendment_rules: "Any change is a successor version through FrozenManifest::amend; \
                          the amendment record names every invalidated claim/obligation \
                          descendant; an amendment after campaign start invalidates the \
                          affected evidence, which must be re-earned; there is no in-place \
                          edit path in the type system",
    }
}

#[allow(clippy::too_many_lines)]
fn pd_claims() -> Vec<ClaimSpec> {
    vec![
        ClaimSpec {
            id: "pd-typed-portfolio-decks",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Every analytic and machine deck in the portfolio registry \
                        names its orbit class (periodic, relative-periodic, hybrid), \
                        phase/gauge/symmetry convention, section identity, and DAE \
                        rank/index declaration; the registry round-trips bit-stably \
                        and refuses every enumerated malformed class with a typed \
                        diagnostic",
            hypotheses: &[
                "decks drawn from the pinned analytic and machine registries",
                "malformed classes are the enumerated ones: unnamed orbit class, \
                 unmapped convention, dangling section reference, undeclared \
                 rank/index",
            ],
            qoi: "registry_roundtrip_and_refusal_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent decoder plus the registry generator's own \
                           validity certificate (emitted during generation, not by \
                           the portfolio compiler)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "first implementation bead of the PD registry leaf opens",
            kill: "any accepted malformed deck or refused valid deck returns the \
                   registry to design review with the disagreeing deck as receipt",
            fallback: "explicit per-deck normal forms with per-field validation",
            no_claim: "no claim that a registered deck is scientifically adequate; \
                       registration is typing, not review",
        },
        ClaimSpec {
            id: "pd-phase-gauge-convention",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "One portfolio-wide phase/gauge/symmetry convention governs \
                        cross-initiative comparison: results are comparable only \
                        through the declared convention map, and any cross-deck \
                        comparison outside the map refuses instead of silently \
                        normalizing",
            hypotheses: &[
                "convention map from the pinned convention fixture",
                "initiative-local conventions declare their embedding into the \
                 portfolio convention",
            ],
            qoi: "convention_enforcement_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent convention walker replaying comparisons \
                           against the declared map (separate from the comparison \
                           code)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "registry leaf green at smoke tier",
            kill: "any silent cross-convention normalization is Sev-0 for the \
                   portfolio lane",
            fallback: "cross-initiative comparisons refuse with the missing \
                       convention edge named",
            no_claim: "no claim that the chosen convention is canonical; only that \
                       it is single and enforced",
        },
        ClaimSpec {
            id: "pd-section-branch-identity",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Section and branch identities are content-addressed and \
                        solver-independent: the same branch produced by different \
                        solvers carries the same identity, distinct sections always \
                        carry distinct identities, and identity collisions or merges \
                        are typed refusals, never silent",
            hypotheses: &[
                "branches and sections from the pinned identity holdout with \
                 ground-truth identity maps",
                "identity inputs are the declared canonical branch invariants, not \
                 solver artifacts",
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
            kill: "any solver-dependent identity or silent merge blocks the identity \
                   lane",
            fallback: "identities carry solver provenance suffixes and say so",
            no_claim: "identity stability is relative to the declared canonical \
                       invariants; a wrong invariant choice is a design defect, not \
                       an identity defect",
        },
        ClaimSpec {
            id: "pd-deferred-event-semantics",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Portfolio hybrid decks defer event/reset semantics to the \
                        I12 automaton contract: every hybrid deck names its automaton \
                        and event budget, and no portfolio lane re-implements event \
                        detection outside that contract",
            hypotheses: &[
                "hybrid decks from the pinned machine registry",
                "the I12 manifest is the semantic authority for events/resets",
            ],
            qoi: "deferral_completeness_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent scanner for event-detection implementations \
                           outside the declared contract references",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "registry leaf green at smoke tier",
            kill: "any shadow event implementation found by the scanner blocks the \
                   hybrid decks",
            fallback: "hybrid decks without named automata refuse at registration",
            no_claim: "no claim about I12's own correctness; deferral structure only",
        },
        ClaimSpec {
            id: "pd-core-maximal-separation",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Core (I09.G4-class) and maximal (I09.G7-class) claims are \
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
            id: "pd-quotient-metric-falsifier",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Refutation,
            statement: "Adversarial symmetry-quotient families (near-symmetric orbit \
                        pairs, symmetry-broken perturbations, relabeled group \
                        actions) attempt to make distinct orbits identify under the \
                        quotient metric or identical orbits separate; any success \
                        refutes the quotient-metric discipline",
            hypotheses: &[
                "fixtures constructed with ground-truth orbit equivalence proofs \
                 in-spec",
                "the quotient metric runs at production settings",
            ],
            qoi: "false_identification_count",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "hand-proved orbit equivalence/inequivalence constructions \
                           carried by the fixture spec",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "runs in every campaign tier from the first quotient-metric \
                         commit",
            kill: "this lane is never killed; it is the standing tripwire",
            fallback: "none: a nonzero count is a release blocker",
            no_claim: "absence of refutation is not completeness and cannot prove \
                       metric soundness; the lane only falsifies",
        },
        ClaimSpec {
            id: "pd-floquet-cluster-metrics",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Portfolio Floquet cluster assignments and quotient metrics \
                        are reproducible across initiatives from the same monodromy \
                        receipts: two lanes consuming one receipt produce identical \
                        cluster assignments, and cluster boundaries carry the \
                        existence-interval provenance they were computed from",
            hypotheses: &[
                "monodromy receipts from the pinned shared family",
                "existence intervals per the I09 waivered semantics (multiplicity \
                 never minted)",
            ],
            qoi: "cluster_reproducibility_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "independent cluster recomputation from the raw receipts \
                           (separate implementation from the portfolio path)",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "identity and lattice leaves green at core tier",
            kill: "any receipt-identical cluster divergence between lanes kills the \
                   metric lane",
            fallback: "per-initiative clusters labeled with their lane identity, no \
                       portfolio aggregation",
            no_claim: "no multiplicity claims; degenerate-cluster qualification is \
                       the waivered lane",
        },
        ClaimSpec {
            id: "pd-objective-gradient-consistency",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Objectives and gradients shared across the portfolio agree \
                        with the I09 adjoint lane's records on shared decks within \
                        the declared band",
            hypotheses: &[
                "shared decks from the pinned registries",
                "the I09 adjoint lane's records are the reference route",
            ],
            qoi: "relative_gradient_discrepancy",
            unit: "1",
            tolerance: ToleranceSemantics::Relative { rtol: 1e-6 },
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "I09 exact-discrete adjoint records (an independent lane \
                           of this portfolio, not the portfolio aggregation path)",
                independent: true,
                tcb_overlap: "shares deck definitions only",
            },
            activation: "identity leaf green at core tier",
            kill: "gate failure on any shared deck kills portfolio gradient \
                   aggregation until root-caused",
            fallback: "per-lane gradients only, no portfolio-level objective \
                       composition",
            no_claim: "consistency is against the I09 records, not against ground \
                       truth; a shared systematic error passes this gate",
        },
        ClaimSpec {
            id: "pd-portfolio-continuation-composition",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "Cross-initiative continuation composition (machine-campaign \
                        branches composed from I09 orbit branches) carries composed \
                        receipts whose lineage names every constituent branch \
                        identity, and composition gaps are localized Unknown windows \
                        with budget receipts, never silent joins",
            hypotheses: &[
                "composition cases from the pinned maximal holdout",
                "constituent branch identities per the [S] identity lane",
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
            activation: "behind the pd-moonshot feature flag on a pre-proof \
                         successor version; frontier lanes green",
            kill: "any silent join on a ground-truth composition fixture kills the \
                   lane",
            fallback: "compositions refuse; constituent branches reported \
                       individually",
            no_claim: "composition lineage is structural; no claim that composed \
                       dynamics are physically valid; version-1 prose mints no \
                       composition authority",
        },
        ClaimSpec {
            id: "pd-maximal-theorem-scope",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "The maximal (I09.G7-class) theorem-scope card carries \
                        machine-checkable obligations for every global claim \
                        (branch-topology completeness, hybrid continuation \
                        qualification), and no global theorem language enters a \
                        portfolio report without its obligation reference",
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
            activation: "behind the pd-moonshot feature flag on a pre-proof \
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

const POLICY_SPEC: &str = "PD_PORTFOLIO_POLICY_V1
DECK_REGISTRY= every deck names orbit class, convention, section identity, and DAE \
rank/index; registration is typing, never scientific review
CONVENTION= one portfolio phase/gauge/symmetry convention; cross-initiative \
comparison only through the declared map; silent normalization is Sev-0
IDENTITY= section/branch identities are content-addressed from declared canonical \
invariants; solver-dependent identities and silent merges refuse
EVENT_SEMANTICS= hybrid decks defer to the I12 automaton contract; shadow event \
implementations are forbidden
LATTICE= core (I09.G4-class) and maximal (I09.G7-class) claims are disjoint lattice \
elements; no cross-lattice gating in either direction
QUOTIENT_METRICS= cluster assignments reproducible from shared receipts; existence \
intervals never mint multiplicity; adversarial quotient families run in every tier
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

fn pd_fixtures() -> Vec<FixturePin> {
    vec![
        FixturePin {
            id: "pd-portfolio-policy-v1",
            source: FixtureSource::AuthoredSpec { spec: POLICY_SPEC },
            partition: Partition::Development,
        },
        FixturePin {
            id: "pd-analytic-deck-registry",
            source: FixtureSource::AuthoredSpec {
                spec: "pd fixture v1: analytic deck registry — Mathieu/Hill tongues, \
                       impact oscillators, torque-free and forced rigid bodies, each \
                       row naming orbit class, convention embedding, section \
                       identity, and rank/index; generator emits validity \
                       certificates; development indices 0..=16383; seeds key \
                       'pd/analytic/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "pd-machine-deck-registry",
            source: FixtureSource::AuthoredSpec {
                spec: "pd fixture v1: machine deck registry — gear-mesh cycles, \
                       Wankel pose campaigns, rotor whirl sweeps, hybrid backlash \
                       variants naming their I12 automata and event budgets; \
                       development indices 0..=16383; seeds key 'pd/machine/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "pd-convention-map",
            source: FixtureSource::AuthoredSpec {
                spec: "pd fixture v1: the portfolio phase/gauge/symmetry convention \
                       map — phase conditions (orthogonality vs Poincare section), \
                       gauge fixings, symmetry group presentations, and the declared \
                       embeddings of every initiative-local convention; development \
                       indices 0..=16383; seeds key 'pd/convention/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "pd-shared-monodromy-family",
            source: FixtureSource::AuthoredSpec {
                spec: "pd fixture v1: shared monodromy receipt family — receipts \
                       with known cluster structures consumed by multiple lanes for \
                       reproducibility checks; development indices 0..=16383; seeds \
                       key 'pd/monodromy/<k>'",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "pd-branch-identity-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "pd fixture v1 HOLDOUT: branch/section identity adversaries — \
                       same branch through different solvers, near-identical \
                       sections, reparameterized branches, with ground-truth \
                       identity maps in-spec; core held-out indices 65536..=69631; \
                       one PD.G3 consumer (pd-identity-lanes); seeds key \
                       'pd/identity/<k>'",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "pd-quotient-adversarial-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "pd fixture v1 HOLDOUT: hand-proved quotient adversaries — \
                       near-symmetric orbit pairs, symmetry-broken perturbations, \
                       relabeled group actions with equivalence proofs in-spec; \
                       core held-out indices 69632..=73727; one PD.G3 consumer \
                       (pd-metrics-quotient); seeds key 'pd/quotient/<k>'",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "pd-composition-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "pd fixture v1 HOLDOUT: cross-initiative composition cases — \
                       machine-campaign branches composed from named orbit branches \
                       with ground-truth composition maps and theorem-scope cards; \
                       maximal held-out indices 131072..=135167; one PD.G3 consumer \
                       (pd-moonshot-composition); seeds key 'pd/composition/<k>'",
            },
            partition: Partition::HeldOut,
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn pd_obligations() -> Vec<ObligationRow> {
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
            leaf: "pd-deck-admission",
            claims_covered: &["pd-phase-gauge-convention", "pd-typed-portfolio-decks"],
            unit_cases: UNIT_CASES,
            g0: "generators: deck registries incl. enumerated malformed classes; \
                 validity predicate: generator certificate agreement; laws: \
                 canonical round-trip bit-stability, refusal totality, convention \
                 map enforcement; shrinker: deck-row removal preserving the violated \
                 rule; replay seeds per explicits",
            decks: &[
                "pd-analytic-deck-registry",
                "pd-convention-map",
                "pd-machine-deck-registry",
                "pd-portfolio-policy-v1",
            ],
            g3_relations: &[
                "deck relabeling invariance",
                "registry-order invariance of admission verdicts",
            ],
            g4_schedule: "request->drain->finalize injection between decode, \
                          convention walk, and registration; checkpoint at each \
                          phase boundary; drained cancellation leaves no partially \
                          registered deck",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/pd_admission.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (pd-admission slice)",
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
            replay_command: "scripts/e2e/leapfrog/pd_admission.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "pd-identity-lanes",
            claims_covered: &["pd-section-branch-identity"],
            unit_cases: UNIT_CASES,
            g0: "generators: branch/section families incl. the identity holdout's \
                 adversaries; laws: solver-independence of identities, distinctness \
                 of distinct sections, typed collision refusals; replay seeds per \
                 explicits",
            decks: &[
                "pd-analytic-deck-registry",
                "pd-branch-identity-core-holdout",
                "pd-portfolio-policy-v1",
            ],
            g3_relations: &[
                "branch reparameterization invariance of identities",
                "solver-route invariance of identities",
            ],
            g4_schedule: "request->drain->finalize injection mid-derivation; \
                          checkpoint at per-branch boundaries; drained cancellation \
                          reports derived identities only",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        identity derivations",
            entry_point: "scripts/e2e/leapfrog/pd_identity.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (pd-identity slice)",
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
            replay_command: "scripts/e2e/leapfrog/pd_identity.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "pd-event-semantics",
            claims_covered: &["pd-deferred-event-semantics"],
            unit_cases: UNIT_CASES,
            g0: "generators: hybrid decks with and without named automata; laws: \
                 deferral completeness (scanner finds no shadow implementations), \
                 refusal of unnamed automata, event-budget presence; replay seeds \
                 per explicits",
            decks: &["pd-machine-deck-registry", "pd-portfolio-policy-v1"],
            g3_relations: &["automaton renaming covariance of deferral verdicts"],
            g4_schedule: "request->drain->finalize injection mid-scan; checkpoint at \
                          per-deck boundaries; drained cancellation reports scanned \
                          decks only",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        deferral verdicts",
            entry_point: "scripts/e2e/leapfrog/pd_events.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (pd-events slice)",
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
            replay_command: "scripts/e2e/leapfrog/pd_events.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "pd-lattice-separation",
            claims_covered: &["pd-core-maximal-separation"],
            unit_cases: UNIT_CASES,
            g0: "generators: claim registries with declared lattice assignments plus \
                 adversarial cross-gating attempts; laws: auditor recomputation \
                 agreement, cross-lattice gating refusal in both directions; replay \
                 seeds per explicits",
            decks: &["pd-analytic-deck-registry", "pd-portfolio-policy-v1"],
            g3_relations: &["claim relabeling invariance of lattice audits"],
            g4_schedule: "request->drain->finalize injection mid-audit; checkpoint \
                          at per-claim boundaries; drained cancellation reports \
                          audited claims only",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        audit verdicts",
            entry_point: "scripts/e2e/leapfrog/pd_lattice.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (pd-lattice slice)",
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
            replay_command: "scripts/e2e/leapfrog/pd_lattice.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "pd-metrics-quotient",
            claims_covered: &[
                "pd-floquet-cluster-metrics",
                "pd-objective-gradient-consistency",
                "pd-quotient-metric-falsifier",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: shared monodromy receipts and the quotient adversarial \
                 holdout; laws: cluster reproducibility across lanes, gradient \
                 agreement within rtol, zero false identifications/separations; \
                 replay seeds per explicits",
            decks: &[
                "pd-portfolio-policy-v1",
                "pd-quotient-adversarial-core-holdout",
                "pd-shared-monodromy-family",
            ],
            g3_relations: &[
                "group-action relabeling invariance of quotient verdicts",
                "receipt-order invariance of cluster assignments",
            ],
            g4_schedule: "request->drain->finalize injection with deliberate \
                          eigensolver-receipt budget exhaustion (typed resumable \
                          outcomes); checkpoint between receipt batches; cancel \
                          mid-clustering with rollback verified",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        clusters, gradients, and falsifier verdicts",
            entry_point: "scripts/e2e/leapfrog/pd_metrics.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (pd-metrics slice)",
            obs_events: &[
                "request.received",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
                "cluster.assigned",
                "gradient.compared",
                "quotient.verdict",
            ],
            replay_command: "scripts/e2e/leapfrog/pd_metrics.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "pd-moonshot-composition",
            claims_covered: &[
                "pd-maximal-theorem-scope",
                "pd-portfolio-continuation-composition",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: composition cases and theorem cards from the maximal \
                 holdout; laws: composed lineage completeness against ground-truth \
                 maps, localized Unknown windows with budget receipts, no \
                 unreferenced theorem language; replay seeds per explicits",
            decks: &[
                "pd-composition-max-holdout",
                "pd-machine-deck-registry",
                "pd-portfolio-policy-v1",
            ],
            g3_relations: &["composition-order invariance of lineage verdicts"],
            g4_schedule: "the core of this lane IS G4: request->drain->finalize \
                          injection with budget exhaustion, timeout, and \
                          cancellation inside composition and scope scanning; \
                          checkpoint at constituent boundaries; each outcome drained \
                          and typed; BudgetExhausted stays Unknown",
            g5_matrix: "threads {1,2,7} x deterministic mode x same-ISA replay of \
                        composed receipts and scope verdicts",
            entry_point: "scripts/e2e/leapfrog/pd_moonshot.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (pd-moonshot slice)",
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
            replay_command: "scripts/e2e/leapfrog/pd_moonshot.sh --manifest \
                             <manifest-id> --replay <artifact-id>",
        },
    ]
}

fn pd_waivers() -> Vec<Waiver> {
    vec![
        Waiver {
            subject: "pd-metrics-quotient",
            reason: "portfolio Floquet cluster metrics consume the fs-spectral \
                     eigensolver service's receipts, whose central proof has not \
                     landed (bead bfid in_progress), and the I09 initiative lanes \
                     that produce shared monodromy receipts are preregistered but \
                     unimplemented",
            owner: "portfolio-pd lane owner",
            predicate: "bead frankensim-ext-spectral-eigensolver-service-bfid closes \
                        with green central proof AND the I09 [S] obligations are \
                        green at core tier",
            expiry: "first PD campaign review after both dependencies land; \
                     re-justify or retire at every manifest amendment",
            promotion_effect: "the [F] cluster-metric claim can close only on \
                               receipt-reproducibility semantics; degenerate-cluster \
                               qualification stays Unknown while the waiver is live",
        },
        Waiver {
            subject: "pd-moonshot-composition",
            reason: "cross-initiative composition requires the I09 continuation \
                     lane and the I12 [S] event lanes green at core tier; neither \
                     initiative has implemented its preregistered lanes yet",
            owner: "portfolio-pd lane owner",
            predicate: "I09 [S]/[F] obligations green at core tier AND I12 [S] \
                        obligations green at core tier",
            expiry: "first PD campaign review after both initiatives land their \
                     core lanes; re-justify or retire at every manifest amendment",
            promotion_effect: "the [M] claims stay Unknown and cannot close; [S]/[F] \
                               promotion is unaffected because their obligations run \
                               on registry and receipt fixtures only",
        },
    ]
}
