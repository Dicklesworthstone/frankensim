//! Focused I11 experimental-campaign/DAQ VerificationManifest
//! conformance.
//!
//! These tests pin the authored campaign authority. They do not run an
//! acquisition, excite a rig, authenticate a real stream, or adjudicate
//! a blind partition and therefore mint no engineering evidence.

use fs_vmanifest::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FixtureSource, FreezeRefusal, ManifestDraft,
    Partition, ToleranceSemantics, i11_draft, obligation_digest,
};
use std::collections::{BTreeMap, BTreeSet};

const POLICY: &str = "i11-campaign-policy-v1";
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
        .unwrap_or_else(|| panic!("missing I11 claim '{id}'"))
}

fn authored_spec<'a>(draft: &'a ManifestDraft, id: &str) -> &'a str {
    let fixture = draft
        .fixtures
        .iter()
        .find(|fixture| fixture.id == id)
        .unwrap_or_else(|| panic!("missing I11 fixture '{id}'"));
    match fixture.source {
        FixtureSource::AuthoredSpec { spec } => spec,
        FixtureSource::External { .. } => panic!("I11 fixture '{id}' must be an authored spec"),
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
fn i11_seed_freezes_with_exact_lattice_and_partition_counts() {
    let draft = i11_draft();
    assert_eq!(draft.initiative, "I11");
    assert_eq!(draft.version, 1);
    assert_eq!(draft.claims.len(), 10);
    assert_eq!(draft.fixtures.len(), 8);
    assert_eq!(draft.obligations.len(), 6);
    assert_eq!(draft.waivers.len(), 2);

    let lattice = draft.claims.iter().fold([0usize; 3], |mut counts, claim| {
        counts[match claim.ambition {
            Ambition::Solid => 0,
            Ambition::Frontier => 1,
            Ambition::Moonshot => 2,
        }] += 1;
        counts
    });
    assert_eq!(lattice, [6, 2, 2]);
    let refutations: Vec<_> = draft
        .claims
        .iter()
        .filter(|claim| claim.polarity == ClaimPolarity::Refutation)
        .map(|claim| claim.id)
        .collect();
    assert_eq!(refutations, ["i11-false-provenance-falsifier"]);

    let held_out: BTreeSet<_> = draft
        .fixtures
        .iter()
        .filter(|fixture| fixture.partition == Partition::HeldOut)
        .map(|fixture| fixture.id)
        .collect();
    assert_eq!(
        held_out,
        BTreeSet::from([
            "i11-blinded-adjudication-max-holdout",
            "i11-calibration-fault-core-holdout",
            "i11-saturation-dropout-core-holdout",
        ])
    );

    let rig_waiver = draft
        .waivers
        .iter()
        .find(|waiver| waiver.subject == "i11-excitation-safety")
        .expect("physical-rig waiver");
    assert!(rig_waiver.reason.contains("synthetic rig models"));
    assert!(
        rig_waiver
            .promotion_effect
            .contains("no physical-rig safety authority")
    );
    let moonshot_waiver = draft
        .waivers
        .iter()
        .find(|waiver| waiver.subject == "i11-moonshot-proof-adaptive")
        .expect("moonshot waiver");
    assert!(moonshot_waiver.predicate.contains("h61n"));
    assert!(
        moonshot_waiver
            .promotion_effect
            .contains("[M] claims stay Unknown")
    );

    let frozen = draft.freeze().expect("the I11 seed must freeze");
    assert_eq!(frozen.initiative(), "I11");
    assert_eq!(frozen.version(), 1);
    assert_eq!(frozen.claims().len(), 10);
    assert_eq!(frozen.fixtures().len(), 8);
    assert_eq!(frozen.obligations().len(), 6);
    assert_eq!(frozen.waivers().len(), 2);
}

#[test]
#[allow(clippy::too_many_lines)]
fn i11_obligation_map_is_complete_once_only_and_executable() {
    let draft = i11_draft();
    let expected: BTreeMap<&str, (CampaignTier, BTreeSet<&str>)> = BTreeMap::from([
        (
            "i11-campaign-admission",
            (
                CampaignTier::Core,
                BTreeSet::from(["i11-typed-campaign-ir", "i11-typed-signal-chain"]),
            ),
        ),
        (
            "i11-excitation-safety",
            (
                CampaignTier::Core,
                BTreeSet::from(["i11-safe-excitation-interlocks"]),
            ),
        ),
        (
            "i11-daq-packs",
            (
                CampaignTier::Core,
                BTreeSet::from(["i11-bounded-daq-packs"]),
            ),
        ),
        (
            "i11-calibration-chains",
            (
                CampaignTier::Core,
                BTreeSet::from(["i11-calibration-chain-integrity"]),
            ),
        ),
        (
            "i11-stream-authentication",
            (
                CampaignTier::Core,
                BTreeSet::from([
                    "i11-authenticated-raw-streams",
                    "i11-blind-partition-discipline",
                    "i11-false-provenance-falsifier",
                ]),
            ),
        ),
        (
            "i11-moonshot-proof-adaptive",
            (
                CampaignTier::Max,
                BTreeSet::from([
                    "i11-adaptive-safe-excitation",
                    "i11-proof-carrying-execution",
                ]),
            ),
        ),
    ]);
    let mut seen = BTreeMap::<&str, usize>::new();
    assert_eq!(draft.obligations.len(), expected.len());

    for row in &draft.obligations {
        let (tier, claims) = expected
            .get(row.leaf)
            .unwrap_or_else(|| panic!("unexpected I11 leaf '{}'", row.leaf));
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
        assert!(row.entry_point.starts_with("scripts/e2e/leapfrog/i11_"));
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
        let authored = i11_draft()
            .obligations
            .into_iter()
            .find(|candidate| candidate.leaf == row.leaf())
            .expect("authored row");
        assert_eq!(row.digest(), obligation_digest(&authored));
    }
}

#[test]
fn i11_safety_acquisition_and_provenance_boundaries_are_pinned() {
    let draft = i11_draft();

    let excitation = claim(&draft, "i11-safe-excitation-interlocks");
    assert!(excitation.statement.contains(
        "normal completion, safe-state completion, and containment \
                       failure are distinct recorded outcomes"
    ));
    assert!(excitation.statement.contains("request->drain->finalize"));
    assert!(
        excitation
            .no_claim
            .contains("not IEC/ISO regulatory certification")
    );

    let daq = claim(&draft, "i11-bounded-daq-packs");
    assert!(daq.statement.contains(
        "a dropout is \
                       flagged, never interpolated silently"
    ));
    assert!(daq.kill.contains("silent interpolation"));

    let calibration = claim(&draft, "i11-calibration-chain-integrity");
    assert!(calibration.statement.contains("never silently passes"));
    assert!(calibration.fallback.contains("Estimated"));

    let falsifier = claim(&draft, "i11-false-provenance-falsifier");
    assert_eq!(falsifier.polarity, ClaimPolarity::Refutation);
    assert_eq!(
        falsifier.tolerance,
        ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 }
    );
    assert!(falsifier.kill.contains("never killed"));
    assert!(falsifier.no_claim.contains("cannot prove"));

    let streams = claim(&draft, "i11-authenticated-raw-streams");
    assert!(streams.no_claim.contains(
        "not sensor \
                       truth"
    ));
    assert!(streams.oracle.tcb_overlap.contains("fs-blake3"));

    let blind = claim(&draft, "i11-blind-partition-discipline");
    assert!(blind.statement.contains(
        "unblinding is a \
                       receipted one-way event"
    ));
    assert!(blind.no_claim.contains("out-of-band leakage"));
}

#[test]
fn i11_policy_is_the_authority_separation_and_retention_spine() {
    let draft = i11_draft();
    let policy = authored_spec(&draft, POLICY);
    assert_eq!(policy.lines().next(), Some("I11_CAMPAIGN_POLICY_V1"));
    for heading in [
        "CAMPAIGN_IR=",
        "SIGNAL_CHAIN=",
        "EXCITATION_SAFETY=",
        "DAQ_PACKS=",
        "CALIBRATION=",
        "STREAM_AUTHENTICITY=",
        "BLINDING=",
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
    assert!(policy.contains("a dropout is flagged, never interpolated silently"));
    assert!(policy.contains("unblinding is a receipted one-way event"));
    assert!(policy.contains("one axis never substitutes for another"));
    assert!(policy.contains("version 1 has prose cards only and mints no proof"));
    assert!(policy.contains("request->drain->finalize"));
    assert!(policy.contains("partial success cannot publish normal authority"));

    for heldout in [
        "i11-saturation-dropout-core-holdout",
        "i11-calibration-fault-core-holdout",
        "i11-blinded-adjudication-max-holdout",
    ] {
        let spec = authored_spec(&draft, heldout);
        assert!(spec.contains("HOLDOUT"));
        assert!(spec.contains("one I11.G3 consumer"));
    }
}

#[test]
fn i11_holdout_ranges_are_disjoint_and_each_has_one_stage_local_consumer() {
    let draft = i11_draft();
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
            "i11-saturation-dropout-core-holdout",
            "65536..=69631",
            "i11-daq-packs",
            CampaignTier::Core,
        ),
        (
            "i11-calibration-fault-core-holdout",
            "69632..=73727",
            "i11-calibration-chains",
            CampaignTier::Core,
        ),
        (
            "i11-blinded-adjudication-max-holdout",
            "131072..=135167",
            "i11-moonshot-proof-adaptive",
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
fn i11_moonshot_ratchets_mint_no_prose_authority() {
    let draft = i11_draft();
    for id in [
        "i11-proof-carrying-execution",
        "i11-adaptive-safe-excitation",
    ] {
        let moonshot = claim(&draft, id);
        assert_eq!(moonshot.ambition, Ambition::Moonshot);
        assert!(
            moonshot.activation.contains(
                "pre-proof \
                         successor"
            ),
            "{id} activation must require a successor version"
        );
        assert!(
            moonshot.no_claim.contains("version-1 prose mints no"),
            "{id} no-claim must disclaim version-1 prose authority"
        );
    }
    let adaptive = claim(&draft, "i11-adaptive-safe-excitation");
    assert!(
        adaptive
            .no_claim
            .contains("no autonomous experiment authority")
    );
    let maximal = draft
        .obligations
        .iter()
        .find(|row| row.leaf == "i11-moonshot-proof-adaptive")
        .expect("maximal row");
    assert_eq!(maximal.tier, CampaignTier::Max);
    assert!(
        maximal
            .g4_schedule
            .contains("BudgetExhausted stays Unknown")
    );
    assert!(
        maximal
            .decks
            .contains(&"i11-blinded-adjudication-max-holdout")
    );
}

#[test]
fn i11_g3_mutations_refuse_or_move_authority() {
    let baseline = i11_draft().freeze().expect("freeze").digest();

    let mut missing_hypotheses = i11_draft();
    missing_hypotheses
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i11-typed-signal-chain")
        .expect("chain claim")
        .hypotheses = &[];
    assert!(matches!(
        missing_hypotheses.freeze(),
        Err(FreezeRefusal::BlankField {
            field: "claim.hypotheses",
            ..
        })
    ));

    let mut correlated = i11_draft();
    correlated
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i11-authenticated-raw-streams")
        .expect("streams claim")
        .oracle
        .independent = false;
    assert!(matches!(
        correlated.freeze(),
        Err(FreezeRefusal::ProductionOracleReuse { .. })
    ));

    let mut relaxed = i11_draft();
    relaxed
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i11-false-provenance-falsifier")
        .expect("falsifier claim")
        .tolerance = ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 };
    assert_ne!(
        relaxed.freeze().expect("relaxed freezes").digest(),
        baseline
    );

    let mut swapped_holdout = i11_draft();
    swapped_holdout
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i11-calibration-fault-core-holdout")
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

    let mut repartitioned = i11_draft();
    repartitioned
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i11-blinded-adjudication-max-holdout")
        .expect("blind holdout")
        .partition = Partition::Development;
    assert_ne!(
        repartitioned
            .freeze()
            .expect("repartition freezes")
            .digest(),
        baseline
    );

    let mut missing_policy = i11_draft();
    missing_policy
        .fixtures
        .retain(|fixture| fixture.id != POLICY);
    assert!(matches!(
        missing_policy.freeze(),
        Err(FreezeRefusal::OrphanDeck { deck, .. }) if deck == POLICY
    ));
}

#[test]
fn i11_g5_top_level_order_is_not_identity() {
    let expected = i11_draft().freeze().expect("freeze");
    let mut permuted = i11_draft();
    permuted.claims.reverse();
    permuted.fixtures.reverse();
    permuted.obligations.reverse();
    permuted.waivers.reverse();
    let actual = permuted.freeze().expect("permuted freeze");
    assert_eq!(actual.digest(), expected.digest());
    assert_eq!(actual, expected);
}

#[test]
fn i11_g4_chunked_in_memory_assembly_is_identity_equivalent() {
    let one_shot = i11_draft();
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
    for chunk in one_shot.claims.chunks(3) {
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
fn i11_amendments_invalidate_exact_targeted_or_global_authority() {
    let predecessor_draft = i11_draft();
    let all = authority_ids(&predecessor_draft);
    let frozen = predecessor_draft.freeze().expect("freeze");

    let mut version_only = i11_draft();
    version_only.version = 2;
    let (_, record) = frozen.amend(version_only).expect("version-only amendment");
    assert!(record.invalidated.is_empty());

    let mut daq = i11_draft();
    daq.version = 2;
    daq.claims
        .iter_mut()
        .find(|claim| claim.id == "i11-bounded-daq-packs")
        .expect("daq claim")
        .statement = "successor acquisition semantics with an intentionally changed \
                      authority identity";
    let (_, daq_record) = frozen.amend(daq).expect("daq amendment");
    assert_eq!(
        daq_record.invalidated,
        vec!["i11-bounded-daq-packs", "i11-daq-packs"]
    );

    let mut holdout = i11_draft();
    holdout.version = 2;
    holdout
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i11-calibration-fault-core-holdout")
        .expect("heldout")
        .source = FixtureSource::AuthoredSpec {
        spec: "successor calibration-fault corpus",
    };
    let (_, holdout_record) = frozen.amend(holdout).expect("holdout amendment");
    assert_eq!(
        holdout_record.invalidated,
        vec!["i11-calibration-chain-integrity", "i11-calibration-chains",]
    );

    let mut policy = i11_draft();
    policy.version = 2;
    policy
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == POLICY)
        .expect("policy")
        .source = FixtureSource::AuthoredSpec {
        spec: "I11_CAMPAIGN_POLICY_V2 intentionally changed global campaign authority",
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
    assert_eq!(policy_record.invalidated.len(), 16);

    let mut title = i11_draft();
    title.version = 2;
    title.title = "successor global I11 campaign authority";
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
