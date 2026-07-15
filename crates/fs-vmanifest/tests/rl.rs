//! Focused RL reality-loop portfolio manifest conformance.
//!
//! The tests lock preregistration authority only. They do not run a loop,
//! validate a lot, calibrate an instrument, replay a campaign, prove the
//! composition theorem, or serve a deployed answer.

use fs_vmanifest::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FixtureSource, FreezeRefusal, ManifestDraft,
    Partition, ToleranceSemantics, obligation_digest, rl_draft,
};
use std::collections::{BTreeMap, BTreeSet};

const POLICY: &str = "rl-campaign-policy-v1";
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
        .unwrap_or_else(|| panic!("missing RL claim '{id}'"))
}

fn authored_spec<'a>(draft: &'a ManifestDraft, id: &str) -> &'a str {
    let fixture = draft
        .fixtures
        .iter()
        .find(|fixture| fixture.id == id)
        .unwrap_or_else(|| panic!("missing RL fixture '{id}'"));
    match fixture.source {
        FixtureSource::AuthoredSpec { spec } => spec,
        FixtureSource::External { .. } => panic!("RL fixture '{id}' must be authored"),
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
fn rl_seed_freezes_with_exact_lattice_corpus_and_waiver() {
    let draft = rl_draft();
    assert_eq!(draft.initiative, "RL");
    assert_eq!(draft.version, 1);
    assert_eq!(draft.claims.len(), 12);
    assert_eq!(draft.fixtures.len(), 12);
    assert_eq!(draft.obligations.len(), 5);
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
    assert_eq!(refutations, ["rl-forged-loop-certificate-falsifier"]);

    let heldout: BTreeSet<_> = draft
        .fixtures
        .iter()
        .filter(|fixture| fixture.partition == Partition::HeldOut)
        .map(|fixture| fixture.id)
        .collect();
    assert_eq!(
        heldout,
        BTreeSet::from([
            "rl-adaptive-hil-max-holdout",
            "rl-blind-loop-core-holdout",
            "rl-identity-custody-core-holdout",
            "rl-loop-certifier-mutants-max-holdout",
            "rl-replay-invalidation-core-holdout",
        ])
    );
    let waiver = draft.waivers[0];
    assert_eq!(waiver.subject, "rl-governed-physical-loop-pack");
    assert!(waiver.reason.contains("synthetic loop worlds"));
    assert!(waiver.predicate.contains("rl-physical-protocol-card"));
    assert!(waiver.promotion_effect.contains("no physical-loop"));

    let frozen = draft.freeze().expect("the RL seed must freeze");
    assert_eq!(frozen.initiative(), "RL");
    assert_eq!(frozen.version(), 1);
    assert_eq!(frozen.claims().len(), 12);
    assert_eq!(frozen.fixtures().len(), 12);
    assert_eq!(frozen.obligations().len(), 5);
    assert_eq!(frozen.waivers().len(), 1);
}

#[test]
#[allow(clippy::too_many_lines)]
fn rl_obligation_map_is_once_only_complete_and_operational() {
    let draft = rl_draft();
    let expected: BTreeMap<&str, (CampaignTier, BTreeSet<&str>)> = BTreeMap::from([
        (
            "rl-identity-calibration-spine",
            (
                CampaignTier::Core,
                BTreeSet::from([
                    "rl-crossstage-identity-continuity",
                    "rl-calibration-measurement-graph",
                ]),
            ),
        ),
        (
            "rl-blind-custody-ownership",
            (
                CampaignTier::Core,
                BTreeSet::from([
                    "rl-blind-partition-custody",
                    "rl-uncertainty-ownership-conservation",
                ]),
            ),
        ),
        (
            "rl-replay-reproof-gate",
            (
                CampaignTier::Core,
                BTreeSet::from([
                    "rl-core-loop-replay",
                    "rl-substitution-selective-reproof",
                    "rl-deployment-gate-composition",
                ]),
            ),
        ),
        (
            "rl-adaptive-hil-loop",
            (
                CampaignTier::Max,
                BTreeSet::from(["rl-adaptive-oed-loop", "rl-hil-timing-composition"]),
            ),
        ),
        (
            "rl-loop-certifier",
            (
                CampaignTier::Max,
                BTreeSet::from([
                    "rl-physical-campaign-parity",
                    "rl-endtoend-composition-theorem",
                    "rl-forged-loop-certificate-falsifier",
                ]),
            ),
        ),
    ]);
    let mut seen = BTreeMap::<&str, usize>::new();
    for row in &draft.obligations {
        let (tier, claims) = expected
            .get(row.leaf)
            .unwrap_or_else(|| panic!("unexpected RL leaf '{}'", row.leaf));
        assert_eq!(&row.tier, tier);
        let actual = row.claims_covered.iter().copied().collect::<BTreeSet<_>>();
        assert_eq!(&actual, claims, "wrong claims on {}", row.leaf);
        for covered in row.claims_covered {
            *seen.entry(covered).or_default() += 1;
        }
        assert_eq!(
            row.unit_cases.iter().copied().collect::<BTreeSet<_>>(),
            UNIT_CASES.into_iter().collect()
        );
        assert!(row.decks.contains(&POLICY));
        assert!(row.entry_point.starts_with("scripts/e2e/leapfrog/rl_"));
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
                "{} omits {event}",
                row.leaf
            );
        }
        for token in ["request->drain->finalize", "checkpoint", "retain"] {
            assert!(
                row.g4_schedule.contains(token),
                "{} omits {token}",
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
        let authored = rl_draft()
            .obligations
            .into_iter()
            .find(|candidate| candidate.leaf == row.leaf())
            .expect("authored row");
        assert_eq!(row.digest(), obligation_digest(&authored));
    }
}

#[test]
fn rl_seam_semantics_bind_the_loop_not_the_stages() {
    let draft = rl_draft();
    let spine = claim(&draft, "rl-crossstage-identity-continuity");
    assert!(spine.statement.contains("content address"));
    assert!(spine.hypotheses.iter().any(|h| h.contains(
        "display names, file paths and \
                 timestamps carry no identity authority"
    ) || h.contains("carry no identity authority")));
    assert!(spine.no_claim.contains("instance manifests"));

    let custody = claim(&draft, "rl-blind-partition-custody");
    assert!(
        custody
            .statement
            .contains("no stage — lot selection, coupon fabrication records, DAQ streams")
            || custody.statement.contains("read another stage")
    );
    assert!(
        custody
            .hypotheses
            .iter()
            .any(|h| h.contains("side channels are in scope"))
    );
    assert!(custody.fallback.contains("burn the contaminated partition"));

    let ownership = claim(&draft, "rl-uncertainty-ownership-conservation");
    assert!(ownership.statement.contains("exactly one owner"));
    assert!(
        ownership
            .statement
            .contains("neither drop nor double-count")
    );
    assert!(ownership.fallback.contains("OwnershipUnknown"));

    let replay = claim(&draft, "rl-core-loop-replay");
    assert!(replay.statement.contains("bit-identical"));
    assert!(replay.no_claim.contains("not correctness"));

    let reproof = claim(&draft, "rl-substitution-selective-reproof");
    assert!(reproof.statement.contains("authenticated lineage"));
    assert!(reproof.statement.contains("invalidated receipt"));
    assert!(
        draft
            .amendment_rules
            .contains("invalidates every RL receipt that bound the predecessor digest")
    );

    let gate = claim(&draft, "rl-deployment-gate-composition");
    for refusal in [
        "StaleEvidence",
        "OutOfDomain",
        "InvalidatedSupport",
        "OwnershipUnknown",
        "Unknown",
    ] {
        assert!(gate.statement.contains(refusal), "gate omits {refusal}");
    }
    assert!(gate.statement.contains("weakest-wins"));
    assert!(gate.no_claim.contains("outside the declared ContextOfUse"));

    let timing = claim(&draft, "rl-hil-timing-composition");
    assert!(timing.statement.contains("target-exact"));
    assert!(timing.statement.contains("never widens"));
    assert!(
        timing
            .hypotheses
            .iter()
            .any(|h| h.contains("never promotes to a WCET bound"))
    );
}

#[test]
fn rl_policy_is_the_authority_and_retention_spine() {
    let draft = rl_draft();
    let policy = authored_spec(&draft, POLICY);
    assert_eq!(policy.lines().next(), Some("RL_CAMPAIGN_POLICY_V1"));
    for heading in [
        "IDENTITY=",
        "STAGES=",
        "BLINDNESS=",
        "UNCERTAINTY=",
        "CALIBRATION=",
        "REPLAY=",
        "REPROOF=",
        "DEPLOYMENT=",
        "ADAPTIVE=",
        "TIMING=",
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
    assert!(policy.contains("a signature authenticates a receipt, not the loop that produced it"));
    assert!(policy.contains("the portfolio owns only the seams"));
    assert!(policy.contains("no stage may read another stage's holdout"));
    assert!(policy.contains("exactly one owner at every stage"));
    assert!(policy.contains("composition is weakest-wins"));
    assert!(policy.contains(
        "a measured sample maximum never promotes to a WCET bound anywhere in the chain"
    ));
    assert!(policy.contains("no automatic physical actuation"));
    assert!(policy.contains("axes never substitute"));
    assert!(policy.contains("partial success cannot publish normal authority"));
}

#[test]
fn rl_holdouts_have_disjoint_ranges_and_one_stage_local_consumer() {
    let draft = rl_draft();
    for token in [
        "development indices 0..=16383",
        "core held-out 65536..=77823",
        "maximal held-out 131072..=143359",
    ] {
        assert!(
            draft.explicits.seeds.contains(token),
            "seed policy omits {token}"
        );
    }
    assert!(draft.explicits.seeds.contains("never reused as loop seeds"));
    let expected = [
        (
            "rl-identity-custody-core-holdout",
            "65536..=69631",
            "rl-identity-calibration-spine",
            CampaignTier::Core,
        ),
        (
            "rl-blind-loop-core-holdout",
            "69632..=73727",
            "rl-blind-custody-ownership",
            CampaignTier::Core,
        ),
        (
            "rl-replay-invalidation-core-holdout",
            "73728..=77823",
            "rl-replay-reproof-gate",
            CampaignTier::Core,
        ),
        (
            "rl-adaptive-hil-max-holdout",
            "131072..=135167",
            "rl-adaptive-hil-loop",
            CampaignTier::Max,
        ),
        (
            "rl-loop-certifier-mutants-max-holdout",
            "135168..=143359",
            "rl-loop-certifier",
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
fn rl_maximal_theorem_parity_and_falsifier_mint_no_prose_authority() {
    let draft = rl_draft();
    let theorem = claim(&draft, "rl-endtoend-composition-theorem");
    assert!(theorem.activation.contains("pre-proof successor"));
    assert!(
        theorem
            .hypotheses
            .iter()
            .any(|h| h.contains("open-world boundaries"))
    );
    assert!(
        theorem
            .no_claim
            .contains("version-1 prose mints no theorem authority")
    );

    let parity = claim(&draft, "rl-physical-campaign-parity");
    assert!(parity.statement.contains("measured"));
    assert!(parity.statement.contains("never assumed"));
    assert!(
        parity
            .activation
            .contains("rl-governed-physical-loop-pack waiver is")
    );
    assert!(parity.no_claim.contains("accreditation"));

    let theorem_card = authored_spec(&draft, "rl-composition-theorem-card");
    for token in [
        "stage-contract/receipt/chain/property ASTs",
        "total premise map",
        "deterministic Lean translation/roundtrip",
        "exact axiom allowlist {propext,Quot.sound,Classical.choice}",
        "transitive closure",
        "nonvacuity witnesses and retained kernel replay",
    ] {
        assert!(theorem_card.contains(token), "theorem card omits {token}");
    }
    assert!(theorem_card.contains("Prose grants no composition authority"));

    let protocol_card = authored_spec(&draft, "rl-physical-protocol-card");
    for token in [
        "receipt-parity requirements",
        "no physical-only bypass",
        "synthetic-to-physical transfer ledgers",
        "third-party replay constraints",
        "grants no physical",
    ] {
        assert!(protocol_card.contains(token), "protocol card omits {token}");
    }

    let falsifier = claim(&draft, "rl-forged-loop-certificate-falsifier");
    assert_eq!(falsifier.polarity, ClaimPolarity::Refutation);
    assert!(falsifier.kill.contains("intended lane succeeds"));
    assert!(
        falsifier
            .no_claim
            .contains("finite adversarial corpus cannot prove")
    );

    let certifier = draft
        .obligations
        .iter()
        .find(|row| row.leaf == "rl-loop-certifier")
        .expect("certifier row");
    for deck in [
        POLICY,
        "rl-composition-theorem-card",
        "rl-physical-protocol-card",
        "rl-loop-certifier-mutants-max-holdout",
        "rl-governed-physical-loop-pack",
    ] {
        assert!(certifier.decks.contains(&deck), "certifier omits {deck}");
    }
    assert!(certifier.g4_schedule.contains("whole-campaign preflight"));
    assert!(
        certifier
            .g4_schedule
            .contains("BudgetExhausted stays Unknown")
    );
}

#[test]
fn rl_g3_mutations_refuse_or_move_authority() {
    let baseline = rl_draft().freeze().expect("freeze").digest();

    let mut missing_hypotheses = rl_draft();
    missing_hypotheses
        .claims
        .iter_mut()
        .find(|claim| claim.id == "rl-core-loop-replay")
        .expect("replay claim")
        .hypotheses = &[];
    assert!(matches!(
        missing_hypotheses.freeze(),
        Err(FreezeRefusal::BlankField {
            field: "claim.hypotheses",
            ..
        })
    ));

    let mut correlated = rl_draft();
    correlated
        .claims
        .iter_mut()
        .find(|claim| claim.id == "rl-endtoend-composition-theorem")
        .expect("theorem claim")
        .oracle
        .independent = false;
    assert!(matches!(
        correlated.freeze(),
        Err(FreezeRefusal::ProductionOracleReuse { .. })
    ));

    let mut relaxed = rl_draft();
    relaxed
        .claims
        .iter_mut()
        .find(|claim| claim.id == "rl-adaptive-oed-loop")
        .expect("adaptive claim")
        .tolerance = ToleranceSemantics::Interval { lo: 0.0, hi: 0.20 };
    assert_ne!(relaxed.freeze().expect("relaxed freeze").digest(), baseline);

    let mut swapped = rl_draft();
    swapped
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "rl-replay-invalidation-core-holdout")
        .expect("replay holdout")
        .source = FixtureSource::AuthoredSpec {
        spec: "unauthorized post-result replacement",
    };
    assert_ne!(
        swapped.freeze().expect("replacement freeze").digest(),
        baseline
    );

    let mut repartitioned = rl_draft();
    repartitioned
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "rl-adaptive-hil-max-holdout")
        .expect("adaptive holdout")
        .partition = Partition::Development;
    assert_ne!(
        repartitioned.freeze().expect("repartition freeze").digest(),
        baseline
    );

    let mut missing_policy = rl_draft();
    missing_policy
        .fixtures
        .retain(|fixture| fixture.id != POLICY);
    assert!(matches!(
        missing_policy.freeze(),
        Err(FreezeRefusal::OrphanDeck { deck, .. }) if deck == POLICY
    ));
}

#[test]
fn rl_g5_top_level_order_is_not_identity() {
    let expected = rl_draft().freeze().expect("freeze");
    let mut permuted = rl_draft();
    permuted.claims.reverse();
    permuted.fixtures.reverse();
    permuted.obligations.reverse();
    permuted.waivers.reverse();
    let actual = permuted.freeze().expect("permuted freeze");
    assert_eq!(actual.digest(), expected.digest());
    assert_eq!(actual, expected);
}

#[test]
fn rl_g4_chunked_in_memory_assembly_is_identity_equivalent() {
    let one_shot = rl_draft();
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
fn rl_amendments_invalidate_exact_targeted_or_global_authority() {
    let predecessor = rl_draft();
    let all = authority_ids(&predecessor);
    let frozen = predecessor.freeze().expect("freeze");

    let mut version_only = rl_draft();
    version_only.version = 2;
    let (_, version_record) = frozen.amend(version_only).expect("version-only");
    assert!(version_record.invalidated.is_empty());

    let mut graph = rl_draft();
    graph.version = 2;
    graph
        .claims
        .iter_mut()
        .find(|claim| claim.id == "rl-calibration-measurement-graph")
        .expect("calibration claim")
        .statement = "successor calibration graph authority";
    let (_, graph_record) = frozen.amend(graph).expect("graph amendment");
    assert_eq!(
        graph_record.invalidated,
        vec![
            "rl-calibration-measurement-graph",
            "rl-identity-calibration-spine",
        ]
    );

    let mut holdout = rl_draft();
    holdout.version = 2;
    holdout
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "rl-blind-loop-core-holdout")
        .expect("blind holdout")
        .source = FixtureSource::AuthoredSpec {
        spec: "successor blind-loop holdout",
    };
    let (_, holdout_record) = frozen.amend(holdout).expect("holdout amendment");
    assert_eq!(
        holdout_record.invalidated,
        vec![
            "rl-blind-custody-ownership",
            "rl-blind-partition-custody",
            "rl-uncertainty-ownership-conservation",
        ]
    );

    let mut policy = rl_draft();
    policy.version = 2;
    policy
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == POLICY)
        .expect("policy")
        .source = FixtureSource::AuthoredSpec {
        spec: "RL_CAMPAIGN_POLICY_V2 changed global authority",
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
    assert_eq!(policy_record.invalidated.len(), 17);

    let mut title = rl_draft();
    title.version = 2;
    title.title = "successor global RL portfolio authority";
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
