//! The I04 (conservation-defect microscope) VerificationManifest draft
//! (bead frankensim-leapfrog-2026-program-i94v.1.3.8.1).
//!
//! Baseline lattice ([S]): local and machine-scale balance attribution,
//! reset/impulse accounting, spatial/temporal/interface defect
//! localization, counterfactual replay with actionable witnesses.
//! Maximal lattice ([F]/[M]): minimal causal explanations across coupled
//! representations, first-divergence guarantees under hybrid topology
//! changes, certified blame stability, and a counterfactual-completeness
//! falsifier lane. A weaker receipt closes its own element and is never
//! relabeled as the stronger theorem.

use crate::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin, FixtureSource,
    GauntletTier, ManifestDraft, ObligationRow, OracleRoute, Partition, ToleranceSemantics, Waiver,
};

/// Build the I04 draft. Consumers freeze it themselves; the conformance
/// battery proves it freezes.
#[must_use]
pub fn i04_draft() -> ManifestDraft {
    ManifestDraft {
        initiative: "I04",
        title: "Conservation-defect microscope gate: balance attribution, impulse \
                 accounting, defect localization, and counterfactual witnesses",
        version: 1,
        explicits: FiveExplicits {
            units: "SI base units throughout; ledger quantities carry their conserved \
                    dimension (J, kg, kg*m/s, C) explicitly; dimensionless QoIs use \
                    unit '1'; exact verdicts use unit 'bit'",
            seeds: "Philox 4x32-10 counter streams keyed 'i04/<fixture-id>/<case-index>'; \
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
        claims: i04_claims(),
        fixtures: i04_fixtures(),
        obligations: i04_obligations(),
        waivers: i04_waivers(),
        amendment_rules: "Any change is a successor version through FrozenManifest::amend; \
                          the amendment record names every invalidated claim/obligation \
                          descendant; an amendment after campaign start invalidates the \
                          affected evidence, which must be re-earned; there is no in-place \
                          edit path in the type system",
    }
}

#[allow(clippy::too_many_lines)]
fn i04_claims() -> Vec<ClaimSpec> {
    vec![
        ClaimSpec {
            id: "i04-balance-attribution",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Local (per-cell/per-port) and machine-scale balance ledgers \
                        attribute every conservation residual to named terms, the \
                        attribution sums reconcile with the total defect under \
                        compensated summation, and manufactured defect magnitudes are \
                        recovered within band",
            hypotheses: &[
                "runs on the pinned manufactured-defect decks with generator-recorded \
                 injection magnitudes",
                "ledger arithmetic is deterministic-mode f64 with compensated sums",
            ],
            qoi: "max_relative_attribution_error",
            unit: "1",
            tolerance: ToleranceSemantics::AbsRel {
                atol: 1e-12,
                rtol: 1e-9,
            },
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "the defect generator's injection record (magnitude, term, \
                           site) — constructed with, not by, the microscope",
                independent: true,
                tcb_overlap: "shares fs-math deterministic scalar kernels",
            },
            activation: "first implementation bead of the I04 ledger leaf opens",
            kill: "an unattributed residual above band on any pinned deck (ledger \
                   leakage) returns the ledger design to review with the leaking deck \
                   as receipt",
            fallback: "total-defect-only reporting, honestly labeled non-attributing",
            no_claim: "attribution is to declared term families; no claim about defects \
                       injected outside the pinned taxonomy",
        },
        ClaimSpec {
            id: "i04-reset-impulse-accounting",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Impulsive events, resets, and mode switches enter the ledger as \
                        explicit accounted entries: no conserved quantity changes across \
                        an event without a ledger row naming the event, and event-free \
                        intervals show zero event-attributed change",
            hypotheses: &[
                "event schedules from the pinned missed-event decks with \
                 generator-recorded event inventories",
                "event detection itself is the upstream solver's job; the microscope \
                 audits the ledger, not the detector",
            ],
            qoi: "unaccounted_event_jump_count",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 },
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "the deck's event inventory replayed against ledger rows by an \
                           independent auditor walk",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "ledger leaf green at smoke tier",
            kill: "any conserved-quantity jump without a ledger row (including \
                   deliberately missed events in the decks) that the auditor cannot \
                   pair kills silent-jump tolerance",
            fallback: "refuse runs whose event inventory cannot be reconciled",
            no_claim: "no claim that the upstream detector finds every physical event; \
                       the claim is that the LEDGER never hides one it was told about \
                       or shows one it was not",
        },
        ClaimSpec {
            id: "i04-defect-localization",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "The microscope localizes each manufactured defect to its exact \
                        injection site: spatial cell/patch, temporal step/window, or \
                        interface id, as recorded by the generator",
            hypotheses: &[
                "one defect per deck in the localization battery (compound defects are \
                 the maximal lane's business)",
                "sites are drawn from the pinned taxonomy: interior cell, boundary \
                 patch, time window, nonconforming interface, scale-transfer joint, \
                 hidden source term",
            ],
            qoi: "localization_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "the defect generator's site record",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "ledger leaf green at smoke tier",
            kill: "systematic mislocalization within any single site class kills that \
                   class's localization claim (the class is named in the receipt)",
            fallback: "region-level (coarse) localization with honest granularity \
                       labels",
            no_claim: "single-defect localization only in this claim; superposed \
                       defects and blame ordering belong to the maximal lanes",
        },
        ClaimSpec {
            id: "i04-counterfactual-replay-witnesses",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Every reported defect carries an actionable witness: a \
                        counterfactual patch (term correction, interface repair, or \
                        event insertion) such that deterministic replay with the patch \
                        removes the attributed defect within band while leaving \
                        unrelated ledger rows unchanged within band",
            hypotheses: &[
                "witness patches come from the declared patch vocabulary frozen with \
                 the fixture families",
                "replay runs under the same seeds and deterministic mode as the \
                 original",
            ],
            qoi: "residual_defect_after_patch",
            unit: "1",
            tolerance: ToleranceSemantics::AbsRel {
                atol: 1e-12,
                rtol: 1e-9,
            },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "deterministic replay harness comparing ledger snapshots \
                           before/after the patch (artifact-level comparison, not \
                           microscope internals)",
                independent: true,
                tcb_overlap: "shares fs-blake3 artifact hashing",
            },
            activation: "localization leaf green at core tier",
            kill: "witnesses that fail their own counterfactual on the pinned decks \
                   demote the microscope to diagnostic-only (no actionable claims)",
            fallback: "defect reports without patches, honestly labeled non-actionable",
            no_claim: "witness patches are counterfactual demonstrations on the pinned \
                       decks, not physical repair recommendations",
        },
        ClaimSpec {
            id: "i04-minimal-causal-explanations",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "For compound defects spanning coupled representations, the \
                        emitted causal explanation is minimal: removing any element \
                        from the explanation leaves a defect the remaining elements \
                        cannot account for, verified by bounded sub-explanation \
                        enumeration",
            hypotheses: &[
                "compound-defect decks with generator-recorded ground-truth cause sets",
                "minimality is with respect to the declared explanation order frozen \
                 with the fixture family",
            ],
            qoi: "minimality_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "bounded sub-explanation enumerator replaying each proper \
                           subset against the ledger",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "all four baseline claims green at core tier",
            kill: "non-minimal explanations on a majority of compound decks demote \
                   explanations to candidate lists",
            fallback: "ranked candidate causes without minimality claims",
            no_claim: "minimality is relative to the frozen explanation order and the \
                       pinned cause taxonomy",
        },
        ClaimSpec {
            id: "i04-first-divergence-hybrid",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Under hybrid topology changes (mode switches that add/remove \
                        equations, ports, or interfaces), the microscope identifies the \
                        FIRST step or interface at which the balance ledger diverges \
                        from the reference, never a later echo",
            hypotheses: &[
                "hybrid decks with generator-recorded first-divergence ground truth",
                "reference trajectories are the decks' recorded clean runs",
            ],
            qoi: "first_divergence_agreement_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "the deck's recorded first-divergence index, generated with \
                           the injected topology change",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "all four baseline claims green at core tier",
            kill: "echo-reporting (any later-than-first divergence) on hybrid decks \
                   kills first-divergence claims across topology changes",
            fallback: "divergence-set reporting without first-ness claims",
            no_claim: "first-ness is with respect to the deck's step/interface \
                       ordering; no wall-clock or causal-time claim",
        },
        ClaimSpec {
            id: "i04-certified-blame-stability",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "Blame assignment is stable across the pinned equivalent chart \
                        representations of the same scene: certified-equivalent pairs \
                        receive the same attributed cause set (up to the pair's \
                        recorded correspondence map)",
            hypotheses: &[
                "chart pairs carry the fixture generator's equivalence certificate \
                 and correspondence map",
            ],
            qoi: "blame_stability_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "the generator's correspondence map applied by an independent \
                           comparator",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "maximal lane opens only after baseline closes at core tier",
            kill: "representation-dependent blame on certified pairs kills chart-generic \
                   blame claims",
            fallback: "per-representation blame with explicit chart labels",
            no_claim: "stability is claimed only for pairs carrying the generator's \
                       certificate, not for arbitrary re-chartings",
        },
        ClaimSpec {
            id: "i04-counterfactual-completeness",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Refutation,
            statement: "Falsifier lane: search for a manufactured defect (within the \
                        pinned taxonomy, including held-out decks) whose true cause is \
                        absent from the microscope's witness set; the completeness \
                        claim survives only while this search comes up empty",
            hypotheses: &[
                "search space is the pinned defect taxonomy across development and \
                 held-out decks",
                "search budget is the max-tier budget in explicits",
            ],
            qoi: "counterexample_count",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "the defect generators' ground-truth cause records over the \
                           full taxonomy sweep",
                independent: true,
                tcb_overlap: "none",
            },
            activation: "maximal lane opens only after baseline closes at core tier",
            kill: "first verified counterexample kills the completeness claim; the \
                   missing-cause deck is retained as the receipt",
            fallback: "completeness restated over the surviving taxonomy subset via \
                       amendment",
            no_claim: "an empty search is Estimated evidence bounded by the taxonomy \
                       and budget; it is never a completeness proof",
        },
    ]
}

fn i04_fixtures() -> Vec<FixturePin> {
    vec![
        FixturePin {
            id: "i04-manufactured-balance-defects",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: clean conservation-ledger scenes (mass/momentum/\
energy/charge) with single injected defects: term scaling, sign \
flip, dropped flux, spurious source; each deck records injection \
site, term, magnitude, and the clean reference trajectory.\n\
SEEDS: Philox stream 'i04/defects/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i04-manufactured-balance-defects-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "Same generator and record discipline as \
i04-manufactured-balance-defects; SEEDS: Philox stream \
'i04/defects/<k>', k in 65536..=69631; withheld until \
adjudication.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i04-missed-event-scenarios",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: event-driven scenes (impacts, resets, switches) \
where the recorded inventory deliberately includes events the \
run must account for AND decks where an event is deliberately \
omitted from the run so the ledger must expose the unaccounted \
jump; inventories and jump magnitudes recorded.\n\
SEEDS: Philox stream 'i04/events/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i04-nonconforming-interfaces",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: coupled scenes with deliberately nonconforming \
interfaces (mismatched quadrature, gap/overlap bands, one-sided \
transfer operators) and recorded interface defect magnitudes.\n\
SEEDS: Philox stream 'i04/interfaces/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i04-scale-transfer-cases",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: scenes with machine-scale transfers (lumped-to-field \
and field-to-lumped couplings) carrying injected transfer-loss \
defects with recorded magnitudes and joints.\n\
SEEDS: Philox stream 'i04/scale/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i04-hidden-source-terms",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: scenes with hidden source/sink terms injected \
outside the declared physics (numerical damping masquerading as \
physics, undeclared stabilization) with recorded term identity \
and magnitude.\n\
SEEDS: Philox stream 'i04/hidden/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i04-equivalent-chart-pairs",
            source: FixtureSource::AuthoredSpec {
                spec: "GENERATOR: the same defective scene expressed in equivalent \
chart representations (re-parameterized fields, permuted \
components, refined-but-equivalent complexes) with an \
equivalence certificate and blame correspondence map per pair.\n\
SEEDS: Philox stream 'i04/charts/<k>', k in 0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i04-equivalent-chart-pairs-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "Same generator as i04-equivalent-chart-pairs; SEEDS: Philox \
stream 'i04/charts/<k>', k in 65536..=69631; withheld until \
adjudication.",
            },
            partition: Partition::HeldOut,
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn i04_obligations() -> Vec<ObligationRow> {
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
            leaf: "i04-ledger-attribution",
            claims_covered: &["i04-balance-attribution"],
            unit_cases: UNIT_CASES,
            g0: "generators: manufactured-defect decks; laws: attribution sums \
                 reconcile with totals under compensated summation, clean decks \
                 attribute zero; shrinker: term-set reduction preserving the injected \
                 defect; replay seeds per explicits",
            decks: &["i04-manufactured-balance-defects"],
            g3_relations: &[
                "unit-rescaling covariance of ledger rows",
                "term relabeling invariance of totals",
            ],
            g4_schedule: "cancel ledger assembly at tile boundaries; a drained ledger \
                          publishes no partial attribution",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i04_ledger.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i04-ledger slice)",
            obs_events: &["ledger.row", "ledger.reconciled", "ledger.cancelled"],
            replay_command: "scripts/e2e/leapfrog/i04_ledger.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i04-event-accounting",
            claims_covered: &["i04-reset-impulse-accounting"],
            unit_cases: UNIT_CASES,
            g0: "generators: missed-event decks; laws: every inventoried event pairs \
                 with a ledger row, event-free intervals show zero event-attributed \
                 change, omitted events surface as unaccounted jumps; replay seeds per \
                 explicits",
            decks: &["i04-missed-event-scenarios"],
            g3_relations: &["event reordering within a step window is ledger-invariant"],
            g4_schedule: "cancel across event boundaries; drained runs leave no \
                          half-accounted event",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i04_events.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i04-events slice)",
            obs_events: &["event.accounted", "event.unaccounted", "event.cancelled"],
            replay_command: "scripts/e2e/leapfrog/i04_events.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i04-localization",
            claims_covered: &["i04-defect-localization"],
            unit_cases: UNIT_CASES,
            g0: "generators: all single-defect families (defects, interfaces, scale, \
                 hidden sources); laws: exact site recovery per taxonomy class; replay \
                 seeds per explicits",
            decks: &[
                "i04-manufactured-balance-defects",
                "i04-nonconforming-interfaces",
                "i04-scale-transfer-cases",
                "i04-hidden-source-terms",
            ],
            g3_relations: &["site-preserving scene permutations leave localization fixed"],
            g4_schedule: "cancel mid-localization; resume reproduces the uninterrupted \
                          site verdict",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i04_localization.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i04-localization slice)",
            obs_events: &["defect.located", "defect.class", "localization.cancelled"],
            replay_command: "scripts/e2e/leapfrog/i04_localization.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i04-counterfactual-replay",
            claims_covered: &["i04-counterfactual-replay-witnesses"],
            unit_cases: UNIT_CASES,
            g0: "generators: single-defect families with the frozen patch vocabulary; \
                 laws: patched replay removes the attributed defect and preserves \
                 unrelated rows within band; replay determinism; seeds per explicits",
            decks: &[
                "i04-manufactured-balance-defects",
                "i04-hidden-source-terms",
            ],
            g3_relations: &["patch idempotence: re-applying a witness patch changes nothing"],
            g4_schedule: "cancel during counterfactual replay; drained replays publish \
                          no partial verdict",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i04_counterfactual.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i04-counterfactual slice)",
            obs_events: &["witness.patch", "replay.verdict", "replay.cancelled"],
            replay_command: "scripts/e2e/leapfrog/i04_counterfactual.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i04-maximal-lanes",
            claims_covered: &[
                "i04-minimal-causal-explanations",
                "i04-first-divergence-hybrid",
                "i04-certified-blame-stability",
                "i04-counterfactual-completeness",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: chart pairs, held-out families, compound-defect and \
                 hybrid-topology sweeps; laws: sub-explanation enumeration replay, \
                 first-divergence ground-truth agreement, correspondence-map \
                 comparison, taxonomy sweep bookkeeping; replay seeds per explicits",
            decks: &[
                "i04-equivalent-chart-pairs",
                "i04-equivalent-chart-pairs-holdout",
                "i04-manufactured-balance-defects-holdout",
                "i04-industrial-defect-benchmark",
            ],
            g3_relations: &[
                "certified chart-pair equivalence",
                "held-out vs development distribution parity checks",
            ],
            g4_schedule: "cancel each maximal lane once mid-campaign; a resumed lane \
                          must report identical verdicts",
            g5_matrix: "threads {1,2,7} x mode {deterministic} x same-ISA replay",
            entry_point: "scripts/e2e/leapfrog/i04_maximal.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (i04-maximal slice)",
            obs_events: &[
                "explanation.minimality",
                "divergence.first",
                "blame.stability",
                "completeness.counterexample",
            ],
            replay_command: "scripts/e2e/leapfrog/i04_maximal.sh --replay <artifact-id>",
        },
    ]
}

fn i04_waivers() -> Vec<Waiver> {
    vec![Waiver {
        subject: "i04-industrial-defect-benchmark",
        reason: "the external industrial defect-hunt benchmark (a vendor or literature \
                 case with a known conservation bug) is not yet licensed or \
                 content-addressed; only its slot is preregistered",
        owner: "I04 implementation beads (frankensim-leapfrog-2026-program-i94v.1.3.x)",
        predicate: "deck bytes pinned with exact edition, license, digest, QoIs, and \
                    acceptance envelopes through the fs-vvreg registry discipline",
        expiry: "before the first Max-tier campaign run; review at each Phase-2 close \
                 burst",
        promotion_effect: "maximal claims may run without this deck but cannot cite it; \
                           any result depending on it stays Estimated until the pin \
                           lands",
    }]
}
