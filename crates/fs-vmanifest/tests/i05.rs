//! Focused I05 deployment-twin/HIL/WCET VerificationManifest conformance.
//!
//! These tests pin the authored campaign authority. They do not execute a
//! target, HIL rig, timing analyzer, proof kernel, or exhaustive search and
//! therefore mint no engineering or theorem evidence.

use fs_vmanifest::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FixtureSource, FreezeRefusal, ManifestDraft,
    Partition, ToleranceSemantics, i05_draft, obligation_digest,
};
use std::collections::{BTreeMap, BTreeSet};

const POLICY: &str = "i05-campaign-policy-v1";
const UNIT_CASES: [&str; 9] = [
    "boundary",
    "cancellation",
    "empty",
    "error",
    "happy",
    "max",
    "migration",
    "tie-break",
    "unit-dimension",
];

fn claim<'a>(draft: &'a ManifestDraft, id: &str) -> &'a ClaimSpec {
    draft
        .claims
        .iter()
        .find(|claim| claim.id == id)
        .unwrap_or_else(|| panic!("missing I05 claim '{id}'"))
}

fn authored_spec<'a>(draft: &'a ManifestDraft, id: &str) -> &'a str {
    let fixture = draft
        .fixtures
        .iter()
        .find(|fixture| fixture.id == id)
        .unwrap_or_else(|| panic!("missing I05 fixture '{id}'"));
    match fixture.source {
        FixtureSource::AuthoredSpec { spec } => spec,
        FixtureSource::External { .. } => panic!("I05 fixture '{id}' must be an authored spec"),
    }
}

fn authority_ids(draft: &ManifestDraft) -> BTreeSet<&'static str> {
    draft
        .claims
        .iter()
        .map(|claim| claim.id)
        .chain(draft.obligations.iter().map(|row| row.leaf))
        .collect()
}

#[test]
fn i05_seed_freezes_with_exact_lattice_and_partition_counts() {
    let draft = i05_draft();
    assert_eq!(draft.initiative, "I05");
    assert_eq!(draft.version, 1);
    assert_eq!(draft.claims.len(), 12);
    assert_eq!(draft.fixtures.len(), 15);
    assert_eq!(draft.obligations.len(), 7);
    assert_eq!(draft.waivers.len(), 1);

    let lattice = draft.claims.iter().fold([0usize; 3], |mut counts, claim| {
        counts[match claim.ambition {
            Ambition::Solid => 0,
            Ambition::Frontier => 1,
            Ambition::Moonshot => 2,
        }] += 1;
        counts
    });
    assert_eq!(lattice, [7, 2, 3]);
    let refutations: Vec<_> = draft
        .claims
        .iter()
        .filter(|claim| claim.polarity == ClaimPolarity::Refutation)
        .map(|claim| claim.id)
        .collect();
    assert_eq!(refutations, ["i05-false-equivalence-certificate-falsifier"]);

    let held_out: BTreeSet<_> = draft
        .fixtures
        .iter()
        .filter(|fixture| fixture.partition == Partition::HeldOut)
        .map(|fixture| fixture.id)
        .collect();
    assert_eq!(
        held_out,
        BTreeSet::from([
            "i05-certificate-mutants-max-holdout",
            "i05-compositional-timing-max-holdout",
            "i05-hil-clock-faults-core-holdout",
            "i05-safe-state-faults-core-holdout",
            "i05-target-refinement-core-holdout",
            "i05-timing-adversaries-core-holdout",
        ])
    );

    let waiver = draft.waivers[0];
    assert_eq!(waiver.subject, "i05-industrial-target-pack");
    assert!(waiver.reason.contains("synthetic target"));
    assert!(waiver.predicate.contains("production-target pack"));
    assert!(waiver.promotion_effect.contains("no silicon WCET"));

    let frozen = draft.freeze().expect("the I05 seed must freeze");
    assert_eq!(frozen.initiative(), "I05");
    assert_eq!(frozen.version(), 1);
    assert_eq!(frozen.claims().len(), 12);
    assert_eq!(frozen.fixtures().len(), 15);
    assert_eq!(frozen.obligations().len(), 7);
    assert_eq!(frozen.waivers().len(), 1);
}

#[test]
#[allow(clippy::too_many_lines)]
fn i05_obligation_map_is_complete_once_only_and_executable() {
    let draft = i05_draft();
    let expected: BTreeMap<&str, (CampaignTier, BTreeSet<&str>)> = BTreeMap::from([
        (
            "i05-ir-numeric-admission",
            (
                CampaignTier::Smoke,
                BTreeSet::from([
                    "i05-deployment-ir-admission",
                    "i05-fixed-point-range-soundness",
                ]),
            ),
        ),
        (
            "i05-reproducible-codegen",
            (
                CampaignTier::Core,
                BTreeSet::from(["i05-reproducible-no-std-binary"]),
            ),
        ),
        (
            "i05-timing-evidence",
            (
                CampaignTier::Core,
                BTreeSet::from(["i05-timing-evidence-separation"]),
            ),
        ),
        (
            "i05-hil-fault-safety",
            (
                CampaignTier::Core,
                BTreeSet::from(["i05-hil-clock-causality", "i05-safe-state-fault-semantics"]),
            ),
        ),
        (
            "i05-target-refinement",
            (
                CampaignTier::Core,
                BTreeSet::from(["i05-bounded-target-refinement"]),
            ),
        ),
        (
            "i05-compositional-timing",
            (
                CampaignTier::Max,
                BTreeSet::from([
                    "i05-compositional-deadline-contracts",
                    "i05-static-measurement-reconciliation",
                ]),
            ),
        ),
        (
            "i05-maximal-certifiers",
            (
                CampaignTier::Max,
                BTreeSet::from([
                    "i05-complete-finite-target-wcet",
                    "i05-false-equivalence-certificate-falsifier",
                    "i05-machine-checked-target-refinement",
                ]),
            ),
        ),
    ]);
    let mut seen = BTreeMap::<&str, usize>::new();
    assert_eq!(draft.obligations.len(), expected.len());

    for row in &draft.obligations {
        let (tier, claims) = expected
            .get(row.leaf)
            .unwrap_or_else(|| panic!("unexpected I05 leaf '{}'", row.leaf));
        assert_eq!(&row.tier, tier, "wrong tier on {}", row.leaf);
        let actual_claims = row.claims_covered.iter().copied().collect::<BTreeSet<_>>();
        assert_eq!(&actual_claims, claims, "wrong claim map on {}", row.leaf);
        for covered in row.claims_covered {
            *seen.entry(covered).or_default() += 1;
        }
        assert_eq!(
            row.unit_cases.iter().copied().collect::<BTreeSet<_>>(),
            UNIT_CASES.into_iter().collect(),
            "all nine unit classes are load-bearing on {}",
            row.leaf
        );
        assert!(row.decks.contains(&POLICY), "{} omits policy", row.leaf);
        assert!(row.entry_point.starts_with("scripts/e2e/leapfrog/i05_"));
        assert!(row.entry_point.ends_with(".sh"));
        assert!(row.replay_command.starts_with(row.entry_point));
        assert!(row.replay_command.contains("--manifest <manifest-id>"));
        assert!(row.replay_command.contains("--replay <artifact-id>"));
        assert!(row.dsr_lane.starts_with("dsr "));
        for event in [
            "request.received",
            "cancel.requested",
            "drain.completed",
            "finalize.completed",
            "failure_bundle.retained",
            "adjudication.receipt",
        ] {
            assert!(
                row.obs_events.contains(&event),
                "{} omits lifecycle event {event}",
                row.leaf
            );
        }
        for token in ["request->drain->finalize", "checkpoint"] {
            assert!(
                row.g4_schedule.contains(token),
                "{} G4 schedule omits {token}",
                row.leaf
            );
        }
        assert!(row.g5_matrix.contains("deterministic mode"));
    }

    assert_eq!(seen.len(), draft.claims.len());
    for claim in &draft.claims {
        assert_eq!(seen.get(claim.id), Some(&1), "{} coverage", claim.id);
    }

    let frozen = draft.freeze().expect("freeze");
    for row in frozen.obligations() {
        assert!(
            row.claims_covered()
                .windows(2)
                .all(|pair| pair[0] < pair[1])
        );
        assert!(row.unit_cases().windows(2).all(|pair| pair[0] < pair[1]));
        assert!(row.decks().windows(2).all(|pair| pair[0] < pair[1]));
        assert!(row.g3_relations().windows(2).all(|pair| pair[0] < pair[1]));
        assert!(row.obs_events().windows(2).all(|pair| pair[0] < pair[1]));
        let authored = i05_draft()
            .obligations
            .into_iter()
            .find(|candidate| candidate.leaf == row.leaf())
            .expect("authored row");
        assert_eq!(row.digest(), obligation_digest(&authored));
    }
}

#[test]
fn i05_numeric_timing_clock_fault_and_refinement_boundaries_are_pinned() {
    let draft = i05_draft();

    let numeric = claim(&draft, "i05-fixed-point-range-soundness");
    for token in [
        "signedness",
        "word lengths",
        "rational scales",
        "rounding",
        "saturation/trap",
        "exact rational fixed-point error",
    ] {
        assert!(
            numeric.statement.contains(token),
            "numeric claim omits {token}"
        );
    }
    assert!(numeric.hypotheses.iter().any(|h| h.contains("wraparound")));
    assert!(numeric.no_claim.contains("closed-loop refinement"));

    let timing = claim(&draft, "i05-timing-evidence-separation");
    for kind in [
        "StaticUpperBound",
        "CompositionalUpperBound",
        "MeasuredSampleMaximum",
        "Unavailable",
    ] {
        assert!(timing.statement.contains(kind));
    }
    assert!(
        timing
            .statement
            .contains("no measured maximum is relabeled as WCET")
    );
    assert!(timing.no_claim.contains("finite testing cannot prove"));

    let clocks = claim(&draft, "i05-hil-clock-causality");
    for token in [
        "simulation",
        "host monotonic",
        "DAQ",
        "bus",
        "device-clock",
        "causality interval",
        "lag/catch-up",
    ] {
        assert!(
            clocks.statement.contains(token),
            "clock claim omits {token}"
        );
    }
    assert!(clocks.no_claim.contains("pinned rig and load"));

    let faults = claim(&draft, "i05-safe-state-fault-semantics");
    assert!(faults.statement.contains("safe output"));
    assert!(faults.kill.contains("unsafe output"));
    assert!(faults.fallback.contains("independent interlock"));
    assert!(
        faults
            .no_claim
            .contains("not IEC/ISO regulatory certification")
    );

    let refinement = claim(&draft, "i05-bounded-target-refinement");
    assert_eq!(
        refinement.tolerance,
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 }
    );
    assert!(refinement.statement.contains("directed interval"));
    assert!(
        refinement
            .hypotheses
            .iter()
            .any(|h| h.contains("never picks the most favorable timestamp"))
    );
    for boundary in [
        "not an unbounded theorem",
        "plant-model validation",
        "stability proof",
        "WCET proof",
    ] {
        assert!(refinement.no_claim.contains(boundary));
    }
}

#[test]
fn i05_policy_is_the_authority_separation_and_retention_spine() {
    let draft = i05_draft();
    let policy = authored_spec(&draft, POLICY);
    assert_eq!(policy.lines().next(), Some("I05_CAMPAIGN_POLICY_V1"));
    for heading in [
        "TARGET_IDENTITY=",
        "ARITHMETIC=",
        "TIMING_KINDS=",
        "CLOCKS=",
        "FAULTS=",
        "REFINEMENT=",
        "THEOREM_AUTHORITY=",
        "EVIDENCE_STATES=",
        "HOLDOUT=",
        "LIFECYCLE=",
        "LOGGING=",
        "RETENTION=",
        "FAILURE_BUNDLE=",
        "PROMOTION=",
        "LEAF_REQUIREMENT=",
    ] {
        assert!(
            policy.lines().any(|line| line.starts_with(heading)),
            "{heading}"
        );
    }
    assert!(policy.contains("measurement is always a lower bound"));
    assert!(policy.contains("never WCET"));
    assert!(policy.contains("normal completion, safe-state completion, containment failure"));
    assert!(policy.contains("one axis never substitutes for another"));
    assert!(policy.contains("version 1 has prose cards only and mints no proof"));
    assert!(policy.contains("request->drain->finalize"));
    assert!(policy.contains("partial success cannot publish normal authority"));

    for heldout in [
        "i05-timing-adversaries-core-holdout",
        "i05-compositional-timing-max-holdout",
        "i05-hil-clock-faults-core-holdout",
        "i05-safe-state-faults-core-holdout",
        "i05-target-refinement-core-holdout",
        "i05-certificate-mutants-max-holdout",
    ] {
        let spec = authored_spec(&draft, heldout);
        assert!(spec.contains("HOLDOUT"));
        assert!(
            spec.contains("one I05.G3 consumer")
                || spec.contains("One I05.G3 consumer")
                || spec.contains("One I05.G6 consumer")
        );
    }
}

#[test]
fn i05_holdout_ranges_are_disjoint_and_each_has_one_stage_local_consumer() {
    let draft = i05_draft();
    for token in [
        "development indices 0..=16383",
        "core held-out indices 65536..=81919",
        "maximal held-out indices 131072..=147455",
    ] {
        assert!(
            draft.explicits.seeds.contains(token),
            "seed policy omits {token}"
        );
    }
    let expected = [
        (
            "i05-timing-adversaries-core-holdout",
            "65536..=69631",
            "i05-timing-evidence",
            CampaignTier::Core,
        ),
        (
            "i05-hil-clock-faults-core-holdout",
            "69632..=73727",
            "i05-hil-fault-safety",
            CampaignTier::Core,
        ),
        (
            "i05-safe-state-faults-core-holdout",
            "73728..=77823",
            "i05-hil-fault-safety",
            CampaignTier::Core,
        ),
        (
            "i05-target-refinement-core-holdout",
            "77824..=81919",
            "i05-target-refinement",
            CampaignTier::Core,
        ),
        (
            "i05-compositional-timing-max-holdout",
            "131072..=135167",
            "i05-compositional-timing",
            CampaignTier::Max,
        ),
        (
            "i05-certificate-mutants-max-holdout",
            "135168..=147455",
            "i05-maximal-certifiers",
            CampaignTier::Max,
        ),
    ];
    for (fixture, range, leaf, tier) in expected {
        assert!(authored_spec(&draft, fixture).contains(range));
        let consumers: Vec<_> = draft
            .obligations
            .iter()
            .filter(|row| row.decks.contains(&fixture))
            .collect();
        assert_eq!(consumers.len(), 1, "{fixture} must have one consumer");
        assert_eq!(consumers[0].leaf, leaf);
        assert_eq!(consumers[0].tier, tier);
    }
}

#[test]
fn i05_maximal_theorem_and_exhaustiveness_ratchets_mint_no_prose_authority() {
    let draft = i05_draft();
    let theorem = claim(&draft, "i05-machine-checked-target-refinement");
    assert!(theorem.activation.contains("pre-proof successor"));
    assert!(
        theorem
            .hypotheses
            .iter()
            .any(|h| h.contains("proposition AST"))
    );
    assert!(
        theorem
            .hypotheses
            .iter()
            .any(|h| h.contains("axiom policy"))
    );
    assert!(
        theorem
            .no_claim
            .contains("version-1 prose mints no theorem authority")
    );

    let wcet = claim(&draft, "i05-complete-finite-target-wcet");
    assert!(
        wcet.activation
            .contains("machine-readable finite semantics")
    );
    assert!(
        wcet.hypotheses
            .iter()
            .any(|h| h.contains("complete executable"))
    );
    assert!(wcet.hypotheses.iter().any(|h| h.contains("rank/unrank")));
    assert!(wcet.no_claim.contains("mints no exhaustive authority"));
    assert!(wcet.no_claim.contains("model-fidelity qualification"));

    let theorem_card = authored_spec(&draft, "i05-refinement-theorem-card");
    for token in [
        "proposition and definition AST bytes/digests",
        "total decoder semantics",
        "premise schema",
        "deterministic Lean translation",
        "exact axiom allowlist {propext,Quot.sound,Classical.choice}",
        "transitive axiom closure",
        "rejection of sorryAx/custom postulates/native or unsafe oracle shortcuts",
        "retained kernel replay",
    ] {
        assert!(theorem_card.contains(token), "theorem card omits {token}");
    }
    assert!(theorem_card.contains("grants no theorem authority"));

    let wcet_card = authored_spec(&draft, "i05-finite-target-wcet-card");
    for token in [
        "complete grammar",
        "validity/initial/transition/cost predicates",
        "total enumeration or sound symbolic coverage",
        "symmetry quotient with proof obligations",
        "rank/unrank/sharding",
        "independent decoder",
        "preflight",
        "Merkle completeness root",
    ] {
        assert!(wcet_card.contains(token), "WCET card omits {token}");
    }
    assert!(wcet_card.contains("grants no exhaustive or silicon-WCET authority"));

    let falsifier = claim(&draft, "i05-false-equivalence-certificate-falsifier");
    assert_eq!(falsifier.polarity, ClaimPolarity::Refutation);
    assert!(falsifier.kill.contains("intended lane succeeds"));
    assert!(
        falsifier
            .no_claim
            .contains("cannot prove certifier soundness")
    );

    let maximal = draft
        .obligations
        .iter()
        .find(|row| row.leaf == "i05-maximal-certifiers")
        .expect("maximal row");
    for deck in [
        POLICY,
        "i05-refinement-theorem-card",
        "i05-finite-target-wcet-card",
        "i05-certificate-mutants-max-holdout",
        "i05-industrial-target-pack",
    ] {
        assert!(maximal.decks.contains(&deck), "maximal row omits {deck}");
    }
    assert!(maximal.g4_schedule.contains("whole-campaign preflight"));
    assert!(
        maximal
            .g4_schedule
            .contains("BudgetExhausted stays Unknown")
    );
    assert!(maximal.g5_matrix.contains("completeness"));
}

#[test]
fn i05_g3_mutations_refuse_or_move_authority() {
    let baseline = i05_draft().freeze().expect("freeze").digest();

    let mut missing_hypotheses = i05_draft();
    missing_hypotheses
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i05-fixed-point-range-soundness")
        .expect("numeric claim")
        .hypotheses = &[];
    assert!(matches!(
        missing_hypotheses.freeze(),
        Err(FreezeRefusal::BlankField {
            field: "claim.hypotheses",
            ..
        })
    ));

    let mut correlated = i05_draft();
    correlated
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i05-timing-evidence-separation")
        .expect("timing claim")
        .oracle
        .independent = false;
    assert!(matches!(
        correlated.freeze(),
        Err(FreezeRefusal::ProductionOracleReuse { .. })
    ));

    let mut relaxed = i05_draft();
    relaxed
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i05-bounded-target-refinement")
        .expect("refinement claim")
        .tolerance = ToleranceSemantics::Interval { lo: 0.0, hi: 2.0 };
    assert_ne!(
        relaxed.freeze().expect("relaxed freezes").digest(),
        baseline
    );

    let mut swapped_holdout = i05_draft();
    swapped_holdout
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i05-target-refinement-core-holdout")
        .expect("heldout")
        .source = FixtureSource::AuthoredSpec {
        spec: "unauthorized post-result replacement",
    };
    assert_ne!(
        swapped_holdout
            .freeze()
            .expect("replacement freezes")
            .digest(),
        baseline
    );

    let mut repartitioned = i05_draft();
    repartitioned
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i05-timing-adversaries-core-holdout")
        .expect("timing holdout")
        .partition = Partition::Development;
    assert_ne!(
        repartitioned
            .freeze()
            .expect("repartition freezes")
            .digest(),
        baseline
    );

    let mut missing_policy = i05_draft();
    missing_policy
        .fixtures
        .retain(|fixture| fixture.id != POLICY);
    assert!(matches!(
        missing_policy.freeze(),
        Err(FreezeRefusal::OrphanDeck { deck, .. }) if deck == POLICY
    ));
}

#[test]
fn i05_g5_top_level_order_is_not_identity() {
    let expected = i05_draft().freeze().expect("freeze");
    let mut permuted = i05_draft();
    permuted.claims.reverse();
    permuted.fixtures.reverse();
    permuted.obligations.reverse();
    permuted.waivers.reverse();
    let actual = permuted.freeze().expect("permuted freeze");
    assert_eq!(actual.digest(), expected.digest());
    assert_eq!(actual, expected);
}

#[test]
fn i05_g4_chunked_in_memory_assembly_is_identity_equivalent() {
    let one_shot = i05_draft();
    let expected = one_shot.clone().freeze().expect("one-shot freeze");
    let mut staged = ManifestDraft {
        initiative: one_shot.initiative,
        title: one_shot.title,
        version: one_shot.version,
        explicits: one_shot.explicits,
        claims: Vec::new(),
        fixtures: Vec::new(),
        obligations: Vec::new(),
        waivers: Vec::new(),
        amendment_rules: one_shot.amendment_rules,
    };
    for chunk in one_shot.claims.chunks(2) {
        staged.claims.extend_from_slice(chunk);
        staged = staged.clone();
    }
    for chunk in one_shot.fixtures.chunks(3) {
        staged.fixtures.extend_from_slice(chunk);
        staged = staged.clone();
    }
    for chunk in one_shot.obligations.chunks(2) {
        staged.obligations.extend_from_slice(chunk);
        staged = staged.clone();
    }
    for chunk in one_shot.waivers.chunks(1) {
        staged.waivers.extend_from_slice(chunk);
        staged = staged.clone();
    }
    let actual = staged.freeze().expect("chunked freeze");
    assert_eq!(actual.digest(), expected.digest());
    assert_eq!(actual, expected);
}

#[test]
fn i05_amendments_invalidate_exact_targeted_or_global_authority() {
    let predecessor_draft = i05_draft();
    let all = authority_ids(&predecessor_draft);
    let frozen = predecessor_draft.freeze().expect("freeze");

    let mut version_only = i05_draft();
    version_only.version = 2;
    let (_, record) = frozen.amend(version_only).expect("version-only amendment");
    assert!(record.invalidated.is_empty());

    let mut numeric = i05_draft();
    numeric.version = 2;
    numeric
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i05-fixed-point-range-soundness")
        .expect("numeric claim")
        .statement = "successor numeric semantics with an intentionally changed authority identity";
    let (_, numeric_record) = frozen.amend(numeric).expect("numeric amendment");
    assert_eq!(
        numeric_record.invalidated,
        vec![
            "i05-fixed-point-range-soundness",
            "i05-ir-numeric-admission",
        ]
    );

    let mut holdout = i05_draft();
    holdout.version = 2;
    holdout
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i05-target-refinement-core-holdout")
        .expect("heldout")
        .source = FixtureSource::AuthoredSpec {
        spec: "successor target-refinement corpus",
    };
    let (_, holdout_record) = frozen.amend(holdout).expect("holdout amendment");
    assert_eq!(
        holdout_record.invalidated,
        vec!["i05-bounded-target-refinement", "i05-target-refinement"]
    );

    let mut policy = i05_draft();
    policy.version = 2;
    policy
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == POLICY)
        .expect("policy")
        .source = FixtureSource::AuthoredSpec {
        spec: "I05_CAMPAIGN_POLICY_V2 intentionally changed global campaign authority",
    };
    let (_, policy_record) = frozen.amend(policy).expect("policy amendment");
    assert_eq!(
        policy_record
            .invalidated
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        all
    );
    assert_eq!(policy_record.invalidated.len(), 19);

    let mut title = i05_draft();
    title.version = 2;
    title.title = "successor global I05 campaign authority";
    let (_, title_record) = frozen.amend(title).expect("title amendment");
    assert_eq!(
        title_record
            .invalidated
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        all
    );
}
