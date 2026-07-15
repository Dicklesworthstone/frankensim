//! Focused CP compiler-portfolio VerificationManifest conformance.
//!
//! These tests pin the authored portfolio authority. They do not
//! compile laws, generate solve plans, or compose runtimes and
//! therefore mint no engineering or theorem evidence.

use fs_vmanifest::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FixtureSource, FreezeRefusal, ManifestDraft,
    Partition, ToleranceSemantics, cp_draft, obligation_digest,
};
use std::collections::{BTreeMap, BTreeSet};

const POLICY: &str = "cp-portfolio-policy-v1";
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
        .unwrap_or_else(|| panic!("missing CP claim '{id}'"))
}

fn authored_spec<'a>(draft: &'a ManifestDraft, id: &str) -> &'a str {
    let fixture = draft
        .fixtures
        .iter()
        .find(|fixture| fixture.id == id)
        .unwrap_or_else(|| panic!("missing CP fixture '{id}'"));
    match fixture.source {
        FixtureSource::AuthoredSpec { spec } => spec,
        FixtureSource::External { .. } => panic!("CP fixture '{id}' must be an authored spec"),
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
fn cp_seed_freezes_with_exact_lattice_and_partition_counts() {
    let draft = cp_draft();
    assert_eq!(draft.initiative, "CP");
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
    assert_eq!(refutations, ["cp-threshold-edit-falsifier"]);

    let held_out: BTreeSet<_> = draft
        .fixtures
        .iter()
        .filter(|fixture| fixture.partition == Partition::HeldOut)
        .map(|fixture| fixture.id)
        .collect();
    assert_eq!(
        held_out,
        BTreeSet::from([
            "cp-composition-max-holdout",
            "cp-receipt-version-core-holdout",
            "cp-threshold-adversarial-core-holdout",
        ])
    );

    let assurance_waiver = draft
        .waivers
        .iter()
        .find(|waiver| waiver.subject == "cp-assurance-metrics")
        .expect("assurance waiver");
    assert!(assurance_waiver.predicate.contains("core G4 receipts"));
    assert!(assurance_waiver.promotion_effect.contains("whole-machine"));
    let composition_waiver = draft
        .waivers
        .iter()
        .find(|waiver| waiver.subject == "cp-moonshot-composition")
        .expect("composition waiver");
    assert!(composition_waiver.predicate.contains("maximal G7 receipts"));
    assert!(
        composition_waiver
            .promotion_effect
            .contains("[M] claims stay Unknown")
    );

    let frozen = draft.freeze().expect("the CP seed must freeze");
    assert_eq!(frozen.initiative(), "CP");
    assert_eq!(frozen.version(), 1);
    assert_eq!(frozen.claims().len(), 10);
    assert_eq!(frozen.fixtures().len(), 8);
    assert_eq!(frozen.obligations().len(), 6);
    assert_eq!(frozen.waivers().len(), 2);
}

#[test]
#[allow(clippy::too_many_lines)]
fn cp_obligation_map_is_complete_once_only_and_executable() {
    let draft = cp_draft();
    let expected: BTreeMap<&str, (CampaignTier, BTreeSet<&str>)> = BTreeMap::from([
        (
            "cp-stage-admission",
            (
                CampaignTier::Core,
                BTreeSet::from(["cp-typed-seam-schemas", "cp-semantic-source-identity"]),
            ),
        ),
        (
            "cp-receipt-lanes",
            (
                CampaignTier::Core,
                BTreeSet::from(["cp-receipt-version-admission"]),
            ),
        ),
        (
            "cp-mode-semantics",
            (
                CampaignTier::Core,
                BTreeSet::from(["cp-mode-ledger-deferral"]),
            ),
        ),
        (
            "cp-lattice-separation",
            (
                CampaignTier::Core,
                BTreeSet::from(["cp-core-maximal-separation"]),
            ),
        ),
        (
            "cp-assurance-metrics",
            (
                CampaignTier::Core,
                BTreeSet::from([
                    "cp-attribution-reproducibility",
                    "cp-solve-plan-consistency",
                    "cp-threshold-edit-falsifier",
                ]),
            ),
        ),
        (
            "cp-moonshot-composition",
            (
                CampaignTier::Max,
                BTreeSet::from(["cp-law-to-runtime-composition", "cp-maximal-theorem-scope"]),
            ),
        ),
    ]);
    let mut seen = BTreeMap::<&str, usize>::new();
    assert_eq!(draft.obligations.len(), expected.len());

    for row in &draft.obligations {
        let (tier, claims) = expected
            .get(row.leaf)
            .unwrap_or_else(|| panic!("unexpected CP leaf '{}'", row.leaf));
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
        assert!(row.entry_point.starts_with("scripts/e2e/leapfrog/cp_"));
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
        let authored = cp_draft()
            .obligations
            .into_iter()
            .find(|candidate| candidate.leaf == row.leaf())
            .expect("authored row");
        assert_eq!(row.digest(), obligation_digest(&authored));
    }
}

#[test]
fn cp_identity_admission_and_lattice_boundaries_are_pinned() {
    let draft = cp_draft();

    let source = claim(&draft, "cp-semantic-source-identity");
    assert!(
        source
            .statement
            .contains("refuses instead of silently normalizing")
    );
    assert!(source.kill.contains("Sev-0"));

    let receipts = claim(&draft, "cp-receipt-version-admission");
    assert!(receipts.statement.contains("content-addressed"));
    assert!(receipts.statement.contains("never silent"));
    assert!(receipts.oracle.tcb_overlap.contains("fs-blake3"));

    let lattice = claim(&draft, "cp-core-maximal-separation");
    assert!(
        lattice
            .statement
            .contains("no maximal evidence requirement gates core promotion")
    );
    assert!(lattice.fallback.contains("none"));

    let deferral = claim(&draft, "cp-mode-ledger-deferral");
    assert!(deferral.statement.contains("ModeLedger contract"));
    assert!(deferral.no_claim.contains("deferral structure only"));

    let falsifier = claim(&draft, "cp-threshold-edit-falsifier");
    assert_eq!(falsifier.polarity, ClaimPolarity::Refutation);
    assert_eq!(
        falsifier.tolerance,
        ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 }
    );
    assert!(falsifier.kill.contains("never killed"));

    let attribution = claim(&draft, "cp-attribution-reproducibility");
    assert!(
        attribution
            .no_claim
            .contains("does not validate physical laws")
    );

    let plans = claim(&draft, "cp-solve-plan-consistency");
    assert_eq!(plans.tolerance, ToleranceSemantics::Relative { rtol: 1e-6 });
    assert!(plans.no_claim.contains("shared systematic error"));
}

#[test]
fn cp_policy_is_the_authority_separation_and_retention_spine() {
    let draft = cp_draft();
    let policy = authored_spec(&draft, POLICY);
    assert_eq!(policy.lines().next(), Some("CP_PORTFOLIO_POLICY_V1"));
    for heading in [
        "SEAM_SCHEMAS=",
        "SOURCE_IDENTITY=",
        "RECEIPT_VERSIONS=",
        "MODE_SEMANTICS=",
        "LATTICE=",
        "ASSURANCE_METRICS=",
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
    assert!(policy.contains("silent seam drift is Sev-0"));
    assert!(policy.contains("exact pinned versions"));
    assert!(policy.contains("one axis never substitutes for another"));
    assert!(policy.contains("version 1 has prose cards only and mints no proof"));
    assert!(policy.contains("request->drain->finalize"));
    assert!(policy.contains("partial success cannot publish normal authority"));

    for heldout in [
        "cp-receipt-version-core-holdout",
        "cp-threshold-adversarial-core-holdout",
        "cp-composition-max-holdout",
    ] {
        let spec = authored_spec(&draft, heldout);
        assert!(spec.contains("HOLDOUT"));
        assert!(spec.contains("one CP.G3 consumer"));
    }
}

#[test]
fn cp_holdout_ranges_are_disjoint_and_each_has_one_stage_local_consumer() {
    let draft = cp_draft();
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
            "cp-receipt-version-core-holdout",
            "65536..=69631",
            "cp-receipt-lanes",
            CampaignTier::Core,
        ),
        (
            "cp-threshold-adversarial-core-holdout",
            "69632..=73727",
            "cp-assurance-metrics",
            CampaignTier::Core,
        ),
        (
            "cp-composition-max-holdout",
            "131072..=135167",
            "cp-moonshot-composition",
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
fn cp_moonshot_ratchets_mint_no_prose_authority() {
    let draft = cp_draft();
    for id in ["cp-law-to-runtime-composition", "cp-maximal-theorem-scope"] {
        let moonshot = claim(&draft, id);
        assert_eq!(moonshot.ambition, Ambition::Moonshot);
        assert!(
            moonshot.activation.contains("pre-proof successor"),
            "{id} activation must require a successor version"
        );
        assert!(
            moonshot.no_claim.contains("version-1 prose mints no"),
            "{id} no-claim must disclaim version-1 prose authority"
        );
    }
    let maximal = draft
        .obligations
        .iter()
        .find(|row| row.leaf == "cp-moonshot-composition")
        .expect("maximal row");
    assert_eq!(maximal.tier, CampaignTier::Max);
    assert!(
        maximal
            .g4_schedule
            .contains("BudgetExhausted stays Unknown")
    );
    assert!(maximal.decks.contains(&"cp-composition-max-holdout"));
}

#[test]
fn cp_g3_mutations_refuse_or_move_authority() {
    let baseline = cp_draft().freeze().expect("freeze").digest();

    let mut missing_hypotheses = cp_draft();
    missing_hypotheses
        .claims
        .iter_mut()
        .find(|claim| claim.id == "cp-semantic-source-identity")
        .expect("source-identity claim")
        .hypotheses = &[];
    assert!(matches!(
        missing_hypotheses.freeze(),
        Err(FreezeRefusal::BlankField {
            field: "claim.hypotheses",
            ..
        })
    ));

    let mut correlated = cp_draft();
    correlated
        .claims
        .iter_mut()
        .find(|claim| claim.id == "cp-attribution-reproducibility")
        .expect("attribution claim")
        .oracle
        .independent = false;
    assert!(matches!(
        correlated.freeze(),
        Err(FreezeRefusal::ProductionOracleReuse { .. })
    ));

    let mut relaxed = cp_draft();
    relaxed
        .claims
        .iter_mut()
        .find(|claim| claim.id == "cp-solve-plan-consistency")
        .expect("plan claim")
        .tolerance = ToleranceSemantics::Relative { rtol: 1e-2 };
    assert_ne!(
        relaxed.freeze().expect("relaxed freezes").digest(),
        baseline
    );

    let mut swapped_holdout = cp_draft();
    swapped_holdout
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "cp-threshold-adversarial-core-holdout")
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

    let mut repartitioned = cp_draft();
    repartitioned
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "cp-composition-max-holdout")
        .expect("composition holdout")
        .partition = Partition::Development;
    assert_ne!(
        repartitioned
            .freeze()
            .expect("repartition freezes")
            .digest(),
        baseline
    );

    let mut missing_policy = cp_draft();
    missing_policy
        .fixtures
        .retain(|fixture| fixture.id != POLICY);
    assert!(matches!(
        missing_policy.freeze(),
        Err(FreezeRefusal::OrphanDeck { deck, .. }) if deck == POLICY
    ));
}

#[test]
fn cp_g5_top_level_order_is_not_identity() {
    let expected = cp_draft().freeze().expect("freeze");
    let mut permuted = cp_draft();
    permuted.claims.reverse();
    permuted.fixtures.reverse();
    permuted.obligations.reverse();
    permuted.waivers.reverse();
    let actual = permuted.freeze().expect("permuted freeze");
    assert_eq!(actual.digest(), expected.digest());
    assert_eq!(actual, expected);
}

#[test]
fn cp_g4_chunked_in_memory_assembly_is_identity_equivalent() {
    let one_shot = cp_draft();
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
fn cp_amendments_invalidate_exact_targeted_or_global_authority() {
    let predecessor_draft = cp_draft();
    let all = authority_ids(&predecessor_draft);
    let frozen = predecessor_draft.freeze().expect("freeze");

    let mut version_only = cp_draft();
    version_only.version = 2;
    let (_, record) = frozen.amend(version_only).expect("version-only amendment");
    assert!(record.invalidated.is_empty());

    let mut admission = cp_draft();
    admission.version = 2;
    admission
        .claims
        .iter_mut()
        .find(|claim| claim.id == "cp-receipt-version-admission")
        .expect("admission claim")
        .statement = "successor receipt-admission semantics with an intentionally \
                      changed authority identity";
    let (_, admission_record) = frozen.amend(admission).expect("admission amendment");
    assert_eq!(
        admission_record.invalidated,
        vec!["cp-receipt-lanes", "cp-receipt-version-admission"]
    );

    let mut holdout = cp_draft();
    holdout.version = 2;
    holdout
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "cp-receipt-version-core-holdout")
        .expect("heldout")
        .source = FixtureSource::AuthoredSpec {
        spec: "successor receipt-adversary corpus",
    };
    let (_, holdout_record) = frozen.amend(holdout).expect("holdout amendment");
    assert_eq!(
        holdout_record.invalidated,
        vec!["cp-receipt-lanes", "cp-receipt-version-admission"]
    );

    let mut policy = cp_draft();
    policy.version = 2;
    policy
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == POLICY)
        .expect("policy")
        .source = FixtureSource::AuthoredSpec {
        spec: "CP_PORTFOLIO_POLICY_V2 intentionally changed global campaign authority",
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

    let mut title = cp_draft();
    title.version = 2;
    title.title = "successor global CP campaign authority";
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
