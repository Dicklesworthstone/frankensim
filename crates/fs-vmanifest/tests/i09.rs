//! Focused I09 periodic-orbit/Floquet VerificationManifest conformance.
//!
//! These tests pin the authored campaign authority. They do not solve an
//! orbit, assemble a monodromy action, run an eigensolver, or continue a
//! branch and therefore mint no engineering or theorem evidence.

use fs_vmanifest::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FixtureSource, FreezeRefusal, ManifestDraft,
    Partition, ToleranceSemantics, i09_draft, obligation_digest,
};
use std::collections::{BTreeMap, BTreeSet};

const POLICY: &str = "i09-campaign-policy-v1";
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
        .unwrap_or_else(|| panic!("missing I09 claim '{id}'"))
}

fn authored_spec<'a>(draft: &'a ManifestDraft, id: &str) -> &'a str {
    let fixture = draft
        .fixtures
        .iter()
        .find(|fixture| fixture.id == id)
        .unwrap_or_else(|| panic!("missing I09 fixture '{id}'"));
    match fixture.source {
        FixtureSource::AuthoredSpec { spec } => spec,
        FixtureSource::External { .. } => panic!("I09 fixture '{id}' must be an authored spec"),
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
fn i09_seed_freezes_with_exact_lattice_and_partition_counts() {
    let draft = i09_draft();
    assert_eq!(draft.initiative, "I09");
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
    assert_eq!(refutations, ["i09-false-stability-falsifier"]);

    let held_out: BTreeSet<_> = draft
        .fixtures
        .iter()
        .filter(|fixture| fixture.partition == Partition::HeldOut)
        .map(|fixture| fixture.id)
        .collect();
    assert_eq!(
        held_out,
        BTreeSet::from([
            "i09-constrained-dae-orbit-core-holdout",
            "i09-defective-monodromy-core-holdout",
            "i09-hybrid-impact-max-holdout",
        ])
    );

    let cluster_waiver = draft
        .waivers
        .iter()
        .find(|waiver| waiver.subject == "i09-floquet-classification")
        .expect("degenerate-cluster waiver");
    assert!(cluster_waiver.predicate.contains("bfid"));
    assert!(cluster_waiver.promotion_effect.contains("existence-only"));
    let hybrid_waiver = draft
        .waivers
        .iter()
        .find(|waiver| waiver.subject == "i09-moonshot-topology-hybrid")
        .expect("hybrid waiver");
    assert!(hybrid_waiver.predicate.contains("ow2o"));
    assert!(
        hybrid_waiver
            .promotion_effect
            .contains("[M] claims stay Unknown")
    );

    let frozen = draft.freeze().expect("the I09 seed must freeze");
    assert_eq!(frozen.initiative(), "I09");
    assert_eq!(frozen.version(), 1);
    assert_eq!(frozen.claims().len(), 10);
    assert_eq!(frozen.fixtures().len(), 8);
    assert_eq!(frozen.obligations().len(), 6);
    assert_eq!(frozen.waivers().len(), 2);
}

#[test]
#[allow(clippy::too_many_lines)]
fn i09_obligation_map_is_complete_once_only_and_executable() {
    let draft = i09_draft();
    let expected: BTreeMap<&str, (CampaignTier, BTreeSet<&str>)> = BTreeMap::from([
        (
            "i09-problem-admission",
            (
                CampaignTier::Core,
                BTreeSet::from([
                    "i09-typed-orbit-problems",
                    "i09-descriptor-infinity-separation",
                ]),
            ),
        ),
        (
            "i09-cross-method-orbits",
            (
                CampaignTier::Core,
                BTreeSet::from(["i09-cross-method-orbit-agreement"]),
            ),
        ),
        (
            "i09-monodromy-adjoint",
            (
                CampaignTier::Core,
                BTreeSet::from([
                    "i09-consistent-tangent-monodromy",
                    "i09-exact-discrete-adjoints",
                ]),
            ),
        ),
        (
            "i09-floquet-classification",
            (
                CampaignTier::Core,
                BTreeSet::from([
                    "i09-certified-multiplier-intervals",
                    "i09-false-stability-falsifier",
                ]),
            ),
        ),
        (
            "i09-continuation",
            (
                CampaignTier::Core,
                BTreeSet::from(["i09-phase-robust-continuation"]),
            ),
        ),
        (
            "i09-moonshot-topology-hybrid",
            (
                CampaignTier::Max,
                BTreeSet::from([
                    "i09-global-branch-topology",
                    "i09-hybrid-grazing-continuation",
                ]),
            ),
        ),
    ]);
    let mut seen = BTreeMap::<&str, usize>::new();
    assert_eq!(draft.obligations.len(), expected.len());

    for row in &draft.obligations {
        let (tier, claims) = expected
            .get(row.leaf)
            .unwrap_or_else(|| panic!("unexpected I09 leaf '{}'", row.leaf));
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
        assert!(row.entry_point.starts_with("scripts/e2e/leapfrog/i09_"));
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
        let authored = i09_draft()
            .obligations
            .into_iter()
            .find(|candidate| candidate.leaf == row.leaf())
            .expect("authored row");
        assert_eq!(row.digest(), obligation_digest(&authored));
    }
}

#[test]
fn i09_stability_multiplier_and_adjoint_boundaries_are_pinned() {
    let draft = i09_draft();

    let intervals = claim(&draft, "i09-certified-multiplier-intervals");
    assert!(
        intervals
            .statement
            .contains("existence intervals never mint multiplicity")
    );
    assert!(intervals.statement.contains("fs-spectral"));
    assert!(
        intervals.no_claim.contains(
            "degenerate \
                       cluster qualification is the waivered lane"
        ) || intervals.no_claim.contains("waivered lane")
    );

    let falsifier = claim(&draft, "i09-false-stability-falsifier");
    assert_eq!(falsifier.polarity, ClaimPolarity::Refutation);
    assert_eq!(
        falsifier.tolerance,
        ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 }
    );
    assert!(falsifier.statement.contains(
        "no measured transient is \
                       ever relabeled as asymptotic stability"
    ));
    assert!(falsifier.kill.contains("never killed"));
    assert!(falsifier.no_claim.contains("cannot prove"));

    let monodromy = claim(&draft, "i09-consistent-tangent-monodromy");
    assert_eq!(
        monodromy.tolerance,
        ToleranceSemantics::Relative { rtol: 1e-6 }
    );
    assert!(monodromy.statement.contains("ACTUAL discrete"));

    let adjoints = claim(&draft, "i09-exact-discrete-adjoints");
    assert_eq!(
        adjoints.tolerance,
        ToleranceSemantics::Relative { rtol: 1e-6 }
    );
    assert!(adjoints.kill.contains("no silent"));
    assert!(adjoints.no_claim.contains("bifurcation"));

    let descriptor = claim(&draft, "i09-descriptor-infinity-separation");
    assert!(descriptor.statement.contains(
        "no infinity mode \
                       is ever reported as a finite stability multiplier"
    ));

    let agreement = claim(&draft, "i09-cross-method-orbit-agreement");
    assert_eq!(
        agreement.tolerance,
        ToleranceSemantics::AbsRel {
            atol: 1e-8,
            rtol: 1e-6,
        }
    );
    assert!(agreement.no_claim.contains("no uniqueness claim"));
}

#[test]
fn i09_policy_is_the_authority_separation_and_retention_spine() {
    let draft = i09_draft();
    let policy = authored_spec(&draft, POLICY);
    assert_eq!(policy.lines().next(), Some("I09_CAMPAIGN_POLICY_V1"));
    for heading in [
        "ORBIT_PROBLEMS=",
        "METHODS=",
        "MONODROMY=",
        "MULTIPLIER_EVIDENCE=",
        "ADJOINTS=",
        "CONTINUATION=",
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
    assert!(policy.contains("existence intervals never mint multiplicity"));
    assert!(policy.contains("no measured transient is ever relabeled"));
    assert!(policy.contains("one axis never substitutes for another"));
    assert!(policy.contains("version 1 has prose cards only and mints no proof"));
    assert!(policy.contains("request->drain->finalize"));
    assert!(policy.contains("partial success cannot publish normal authority"));

    for heldout in [
        "i09-constrained-dae-orbit-core-holdout",
        "i09-defective-monodromy-core-holdout",
        "i09-hybrid-impact-max-holdout",
    ] {
        let spec = authored_spec(&draft, heldout);
        assert!(spec.contains("HOLDOUT"));
        assert!(spec.contains("one I09.G3 consumer"));
    }
}

#[test]
fn i09_holdout_ranges_are_disjoint_and_each_has_one_stage_local_consumer() {
    let draft = i09_draft();
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
            "i09-constrained-dae-orbit-core-holdout",
            "65536..=69631",
            "i09-problem-admission",
            CampaignTier::Core,
        ),
        (
            "i09-defective-monodromy-core-holdout",
            "69632..=73727",
            "i09-floquet-classification",
            CampaignTier::Core,
        ),
        (
            "i09-hybrid-impact-max-holdout",
            "131072..=135167",
            "i09-moonshot-topology-hybrid",
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
fn i09_moonshot_ratchets_mint_no_prose_authority() {
    let draft = i09_draft();
    for id in [
        "i09-global-branch-topology",
        "i09-hybrid-grazing-continuation",
    ] {
        let moonshot = claim(&draft, id);
        assert_eq!(moonshot.ambition, Ambition::Moonshot);
        assert!(
            moonshot.activation.contains("pre-proof successor"),
            "{id} activation must require a successor version"
        );
        assert!(
            moonshot.no_claim.contains(
                "version-1 \
                       prose mints no"
            ),
            "{id} no-claim must disclaim version-1 prose authority"
        );
    }
    let maximal = draft
        .obligations
        .iter()
        .find(|row| row.leaf == "i09-moonshot-topology-hybrid")
        .expect("maximal row");
    assert_eq!(maximal.tier, CampaignTier::Max);
    assert!(
        maximal
            .g4_schedule
            .contains("BudgetExhausted stays Unknown")
    );
    assert!(maximal.decks.contains(&"i09-hybrid-impact-max-holdout"));
}

#[test]
fn i09_g3_mutations_refuse_or_move_authority() {
    let baseline = i09_draft().freeze().expect("freeze").digest();

    let mut missing_hypotheses = i09_draft();
    missing_hypotheses
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i09-consistent-tangent-monodromy")
        .expect("monodromy claim")
        .hypotheses = &[];
    assert!(matches!(
        missing_hypotheses.freeze(),
        Err(FreezeRefusal::BlankField {
            field: "claim.hypotheses",
            ..
        })
    ));

    let mut correlated = i09_draft();
    correlated
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i09-cross-method-orbit-agreement")
        .expect("agreement claim")
        .oracle
        .independent = false;
    assert!(matches!(
        correlated.freeze(),
        Err(FreezeRefusal::ProductionOracleReuse { .. })
    ));

    let mut relaxed = i09_draft();
    relaxed
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i09-exact-discrete-adjoints")
        .expect("adjoint claim")
        .tolerance = ToleranceSemantics::Relative { rtol: 1e-2 };
    assert_ne!(
        relaxed.freeze().expect("relaxed freezes").digest(),
        baseline
    );

    let mut swapped_holdout = i09_draft();
    swapped_holdout
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i09-defective-monodromy-core-holdout")
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

    let mut repartitioned = i09_draft();
    repartitioned
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i09-hybrid-impact-max-holdout")
        .expect("impact holdout")
        .partition = Partition::Development;
    assert_ne!(
        repartitioned
            .freeze()
            .expect("repartition freezes")
            .digest(),
        baseline
    );

    let mut missing_policy = i09_draft();
    missing_policy
        .fixtures
        .retain(|fixture| fixture.id != POLICY);
    assert!(matches!(
        missing_policy.freeze(),
        Err(FreezeRefusal::OrphanDeck { deck, .. }) if deck == POLICY
    ));
}

#[test]
fn i09_g5_top_level_order_is_not_identity() {
    let expected = i09_draft().freeze().expect("freeze");
    let mut permuted = i09_draft();
    permuted.claims.reverse();
    permuted.fixtures.reverse();
    permuted.obligations.reverse();
    permuted.waivers.reverse();
    let actual = permuted.freeze().expect("permuted freeze");
    assert_eq!(actual.digest(), expected.digest());
    assert_eq!(actual, expected);
}

#[test]
fn i09_g4_chunked_in_memory_assembly_is_identity_equivalent() {
    let one_shot = i09_draft();
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
fn i09_amendments_invalidate_exact_targeted_or_global_authority() {
    let predecessor_draft = i09_draft();
    let all = authority_ids(&predecessor_draft);
    let frozen = predecessor_draft.freeze().expect("freeze");

    let mut version_only = i09_draft();
    version_only.version = 2;
    let (_, record) = frozen.amend(version_only).expect("version-only amendment");
    assert!(record.invalidated.is_empty());

    let mut agreement = i09_draft();
    agreement.version = 2;
    agreement
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i09-cross-method-orbit-agreement")
        .expect("agreement claim")
        .statement = "successor cross-method semantics with an intentionally changed \
                      authority identity";
    let (_, agreement_record) = frozen.amend(agreement).expect("agreement amendment");
    assert_eq!(
        agreement_record.invalidated,
        vec![
            "i09-cross-method-orbit-agreement",
            "i09-cross-method-orbits",
        ]
    );

    let mut holdout = i09_draft();
    holdout.version = 2;
    holdout
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i09-defective-monodromy-core-holdout")
        .expect("heldout")
        .source = FixtureSource::AuthoredSpec {
        spec: "successor adversarial monodromy corpus",
    };
    let (_, holdout_record) = frozen.amend(holdout).expect("holdout amendment");
    assert_eq!(
        holdout_record.invalidated,
        vec![
            "i09-certified-multiplier-intervals",
            "i09-false-stability-falsifier",
            "i09-floquet-classification",
        ]
    );

    let mut policy = i09_draft();
    policy.version = 2;
    policy
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == POLICY)
        .expect("policy")
        .source = FixtureSource::AuthoredSpec {
        spec: "I09_CAMPAIGN_POLICY_V2 intentionally changed global campaign authority",
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

    let mut title = i09_draft();
    title.version = 2;
    title.title = "successor global I09 campaign authority";
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
