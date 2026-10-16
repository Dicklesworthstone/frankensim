//! Focused PD periodic-dynamics portfolio VerificationManifest
//! conformance.
//!
//! These tests pin the authored portfolio authority. They do not solve
//! orbits, cluster multipliers, or compose campaigns and therefore mint
//! no engineering or theorem evidence.

use fs_vmanifest::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FixtureSource, FreezeRefusal, ManifestDraft,
    Partition, ToleranceSemantics, obligation_digest, pd_draft,
};
use std::collections::{BTreeMap, BTreeSet};

const POLICY: &str = "pd-portfolio-policy-v1";
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
        .unwrap_or_else(|| panic!("missing PD claim '{id}'"))
}

fn authored_spec<'a>(draft: &'a ManifestDraft, id: &str) -> &'a str {
    let fixture = draft
        .fixtures
        .iter()
        .find(|fixture| fixture.id == id)
        .unwrap_or_else(|| panic!("missing PD fixture '{id}'"));
    match fixture.source {
        FixtureSource::AuthoredSpec { spec } => spec,
        FixtureSource::External { .. } => panic!("PD fixture '{id}' must be an authored spec"),
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
fn pd_seed_freezes_with_exact_lattice_and_partition_counts() {
    let draft = pd_draft();
    assert_eq!(draft.initiative, "PD");
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
    assert_eq!(refutations, ["pd-quotient-metric-falsifier"]);

    let held_out: BTreeSet<_> = draft
        .fixtures
        .iter()
        .filter(|fixture| fixture.partition == Partition::HeldOut)
        .map(|fixture| fixture.id)
        .collect();
    assert_eq!(
        held_out,
        BTreeSet::from([
            "pd-branch-identity-core-holdout",
            "pd-composition-max-holdout",
            "pd-quotient-adversarial-core-holdout",
        ])
    );

    let metrics_waiver = draft
        .waivers
        .iter()
        .find(|waiver| waiver.subject == "pd-metrics-quotient")
        .expect("metrics waiver");
    assert!(metrics_waiver.predicate.contains("bfid"));
    assert!(
        metrics_waiver
            .promotion_effect
            .contains("degenerate-cluster")
    );
    let composition_waiver = draft
        .waivers
        .iter()
        .find(|waiver| waiver.subject == "pd-moonshot-composition")
        .expect("composition waiver");
    assert!(composition_waiver.predicate.contains("I09 [S]/[F]"));
    assert!(
        composition_waiver
            .promotion_effect
            .contains("[M] claims stay Unknown")
    );

    let frozen = draft.freeze().expect("the PD seed must freeze");
    assert_eq!(frozen.initiative(), "PD");
    assert_eq!(frozen.version(), 1);
    assert_eq!(frozen.claims().len(), 10);
    assert_eq!(frozen.fixtures().len(), 8);
    assert_eq!(frozen.obligations().len(), 6);
    assert_eq!(frozen.waivers().len(), 2);
}

#[test]
#[allow(clippy::too_many_lines)]
fn pd_obligation_map_is_complete_once_only_and_executable() {
    let draft = pd_draft();
    let expected: BTreeMap<&str, (CampaignTier, BTreeSet<&str>)> = BTreeMap::from([
        (
            "pd-deck-admission",
            (
                CampaignTier::Core,
                BTreeSet::from(["pd-typed-portfolio-decks", "pd-phase-gauge-convention"]),
            ),
        ),
        (
            "pd-identity-lanes",
            (
                CampaignTier::Core,
                BTreeSet::from(["pd-section-branch-identity"]),
            ),
        ),
        (
            "pd-event-semantics",
            (
                CampaignTier::Core,
                BTreeSet::from(["pd-deferred-event-semantics"]),
            ),
        ),
        (
            "pd-lattice-separation",
            (
                CampaignTier::Core,
                BTreeSet::from(["pd-core-maximal-separation"]),
            ),
        ),
        (
            "pd-metrics-quotient",
            (
                CampaignTier::Core,
                BTreeSet::from([
                    "pd-floquet-cluster-metrics",
                    "pd-objective-gradient-consistency",
                    "pd-quotient-metric-falsifier",
                ]),
            ),
        ),
        (
            "pd-moonshot-composition",
            (
                CampaignTier::Max,
                BTreeSet::from([
                    "pd-maximal-theorem-scope",
                    "pd-portfolio-continuation-composition",
                ]),
            ),
        ),
    ]);
    let mut seen = BTreeMap::<&str, usize>::new();
    assert_eq!(draft.obligations.len(), expected.len());

    for row in &draft.obligations {
        let (tier, claims) = expected
            .get(row.leaf)
            .unwrap_or_else(|| panic!("unexpected PD leaf '{}'", row.leaf));
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
        assert!(row.entry_point.starts_with("scripts/e2e/leapfrog/pd_"));
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
        let authored = pd_draft()
            .obligations
            .into_iter()
            .find(|candidate| candidate.leaf == row.leaf())
            .expect("authored row");
        assert_eq!(row.digest(), obligation_digest(&authored));
    }
}

#[test]
fn pd_convention_identity_and_lattice_boundaries_are_pinned() {
    let draft = pd_draft();

    let convention = claim(&draft, "pd-phase-gauge-convention");
    assert!(convention.statement.contains(
        "refuses instead of silently \
                        normalizing"
    ));
    assert!(convention.kill.contains("Sev-0"));

    let identity = claim(&draft, "pd-section-branch-identity");
    assert!(identity.statement.contains("solver-independent"));
    assert!(identity.statement.contains("never silent"));
    assert!(identity.oracle.tcb_overlap.contains("fs-blake3"));

    let lattice = claim(&draft, "pd-core-maximal-separation");
    assert!(
        lattice
            .statement
            .contains("no maximal evidence requirement gates core promotion")
    );
    assert!(lattice.fallback.contains("none"));

    let deferral = claim(&draft, "pd-deferred-event-semantics");
    assert!(deferral.statement.contains("I12 automaton contract"));
    assert!(deferral.no_claim.contains("deferral structure only"));

    let falsifier = claim(&draft, "pd-quotient-metric-falsifier");
    assert_eq!(falsifier.polarity, ClaimPolarity::Refutation);
    assert_eq!(
        falsifier.tolerance,
        ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 }
    );
    assert!(falsifier.kill.contains("never killed"));

    let clusters = claim(&draft, "pd-floquet-cluster-metrics");
    assert!(clusters.no_claim.contains("no multiplicity claims"));

    let gradients = claim(&draft, "pd-objective-gradient-consistency");
    assert_eq!(
        gradients.tolerance,
        ToleranceSemantics::Relative { rtol: 1e-6 }
    );
    assert!(gradients.no_claim.contains("shared systematic error"));
}

#[test]
fn pd_policy_is_the_authority_separation_and_retention_spine() {
    let draft = pd_draft();
    let policy = authored_spec(&draft, POLICY);
    assert_eq!(policy.lines().next(), Some("PD_PORTFOLIO_POLICY_V1"));
    for heading in [
        "DECK_REGISTRY=",
        "CONVENTION=",
        "IDENTITY=",
        "EVENT_SEMANTICS=",
        "LATTICE=",
        "QUOTIENT_METRICS=",
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
    assert!(policy.contains("silent normalization is Sev-0"));
    assert!(policy.contains(
        "existence \
intervals never mint multiplicity"
    ));
    assert!(policy.contains("one axis never substitutes for another"));
    assert!(policy.contains("version 1 has prose cards only and mints no proof"));
    assert!(policy.contains("request->drain->finalize"));
    assert!(policy.contains("partial success cannot publish normal authority"));

    for heldout in [
        "pd-branch-identity-core-holdout",
        "pd-quotient-adversarial-core-holdout",
        "pd-composition-max-holdout",
    ] {
        let spec = authored_spec(&draft, heldout);
        assert!(spec.contains("HOLDOUT"));
        assert!(spec.contains("one PD.G3 consumer"));
    }
}

#[test]
fn pd_holdout_ranges_are_disjoint_and_each_has_one_stage_local_consumer() {
    let draft = pd_draft();
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
            "pd-branch-identity-core-holdout",
            "65536..=69631",
            "pd-identity-lanes",
            CampaignTier::Core,
        ),
        (
            "pd-quotient-adversarial-core-holdout",
            "69632..=73727",
            "pd-metrics-quotient",
            CampaignTier::Core,
        ),
        (
            "pd-composition-max-holdout",
            "131072..=135167",
            "pd-moonshot-composition",
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
fn pd_moonshot_ratchets_mint_no_prose_authority() {
    let draft = pd_draft();
    for id in [
        "pd-portfolio-continuation-composition",
        "pd-maximal-theorem-scope",
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
        .find(|row| row.leaf == "pd-moonshot-composition")
        .expect("maximal row");
    assert_eq!(maximal.tier, CampaignTier::Max);
    assert!(
        maximal
            .g4_schedule
            .contains("BudgetExhausted stays Unknown")
    );
    assert!(maximal.decks.contains(&"pd-composition-max-holdout"));
}

#[test]
fn pd_g3_mutations_refuse_or_move_authority() {
    let baseline = pd_draft().freeze().expect("freeze").digest();

    let mut missing_hypotheses = pd_draft();
    missing_hypotheses
        .claims
        .iter_mut()
        .find(|claim| claim.id == "pd-phase-gauge-convention")
        .expect("convention claim")
        .hypotheses = &[];
    assert!(matches!(
        missing_hypotheses.freeze(),
        Err(FreezeRefusal::BlankField {
            field: "claim.hypotheses",
            ..
        })
    ));

    let mut correlated = pd_draft();
    correlated
        .claims
        .iter_mut()
        .find(|claim| claim.id == "pd-floquet-cluster-metrics")
        .expect("cluster claim")
        .oracle
        .independent = false;
    assert!(matches!(
        correlated.freeze(),
        Err(FreezeRefusal::ProductionOracleReuse { .. })
    ));

    let mut relaxed = pd_draft();
    relaxed
        .claims
        .iter_mut()
        .find(|claim| claim.id == "pd-objective-gradient-consistency")
        .expect("gradient claim")
        .tolerance = ToleranceSemantics::Relative { rtol: 1e-2 };
    assert_ne!(
        relaxed.freeze().expect("relaxed freezes").digest(),
        baseline
    );

    let mut swapped_holdout = pd_draft();
    swapped_holdout
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "pd-quotient-adversarial-core-holdout")
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

    let mut repartitioned = pd_draft();
    repartitioned
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "pd-composition-max-holdout")
        .expect("composition holdout")
        .partition = Partition::Development;
    assert_ne!(
        repartitioned
            .freeze()
            .expect("repartition freezes")
            .digest(),
        baseline
    );

    let mut missing_policy = pd_draft();
    missing_policy
        .fixtures
        .retain(|fixture| fixture.id != POLICY);
    assert!(matches!(
        missing_policy.freeze(),
        Err(FreezeRefusal::OrphanDeck { deck, .. }) if deck == POLICY
    ));
}

#[test]
fn pd_g5_top_level_order_is_not_identity() {
    let expected = pd_draft().freeze().expect("freeze");
    let mut permuted = pd_draft();
    permuted.claims.reverse();
    permuted.fixtures.reverse();
    permuted.obligations.reverse();
    permuted.waivers.reverse();
    let actual = permuted.freeze().expect("permuted freeze");
    assert_eq!(actual.digest(), expected.digest());
    assert_eq!(actual, expected);
}

#[test]
fn pd_g4_chunked_in_memory_assembly_is_identity_equivalent() {
    let one_shot = pd_draft();
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
fn pd_amendments_invalidate_exact_targeted_or_global_authority() {
    let predecessor_draft = pd_draft();
    let all = authority_ids(&predecessor_draft);
    let frozen = predecessor_draft.freeze().expect("freeze");

    let mut version_only = pd_draft();
    version_only.version = 2;
    let (_, record) = frozen.amend(version_only).expect("version-only amendment");
    assert!(record.invalidated.is_empty());

    let mut identity = pd_draft();
    identity.version = 2;
    identity
        .claims
        .iter_mut()
        .find(|claim| claim.id == "pd-section-branch-identity")
        .expect("identity claim")
        .statement = "successor identity semantics with an intentionally changed \
                      authority identity";
    let (_, identity_record) = frozen.amend(identity).expect("identity amendment");
    assert_eq!(
        identity_record.invalidated,
        vec!["pd-identity-lanes", "pd-section-branch-identity"]
    );

    let mut holdout = pd_draft();
    holdout.version = 2;
    holdout
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "pd-branch-identity-core-holdout")
        .expect("heldout")
        .source = FixtureSource::AuthoredSpec {
        spec: "successor identity-adversary corpus",
    };
    let (_, holdout_record) = frozen.amend(holdout).expect("holdout amendment");
    assert_eq!(
        holdout_record.invalidated,
        vec!["pd-identity-lanes", "pd-section-branch-identity"]
    );

    let mut policy = pd_draft();
    policy.version = 2;
    policy
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == POLICY)
        .expect("policy")
        .source = FixtureSource::AuthoredSpec {
        spec: "PD_PORTFOLIO_POLICY_V2 intentionally changed global campaign authority",
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

    let mut title = pd_draft();
    title.version = 2;
    title.title = "successor global PD campaign authority";
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
