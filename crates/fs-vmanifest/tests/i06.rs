//! Focused I06 material-lot passport/substitution-impact manifest conformance.
//!
//! The tests lock preregistration authority only. They do not authenticate a
//! physical lot, calibrate a model, identify a causal effect, approve a
//! substitution, prove impact completeness, or optimize a real decision.

use fs_vmanifest::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FixtureSource, FreezeRefusal, ManifestDraft,
    Partition, ToleranceSemantics, i06_draft, obligation_digest,
};
use std::collections::{BTreeMap, BTreeSet};

const POLICY: &str = "i06-campaign-policy-v1";
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
        .unwrap_or_else(|| panic!("missing I06 claim '{id}'"))
}

fn authored_spec<'a>(draft: &'a ManifestDraft, id: &str) -> &'a str {
    let fixture = draft
        .fixtures
        .iter()
        .find(|fixture| fixture.id == id)
        .unwrap_or_else(|| panic!("missing I06 fixture '{id}'"));
    match fixture.source {
        FixtureSource::AuthoredSpec { spec } => spec,
        FixtureSource::External { .. } => panic!("I06 fixture '{id}' must be authored"),
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
fn i06_seed_freezes_with_exact_lattice_corpus_and_waiver() {
    let draft = i06_draft();
    assert_eq!(draft.initiative, "I06");
    assert_eq!(draft.version, 1);
    assert_eq!(draft.claims.len(), 12);
    assert_eq!(draft.fixtures.len(), 17);
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
    assert_eq!(
        refutations,
        ["i06-false-provenance-impact-certificate-falsifier"]
    );

    let heldout: BTreeSet<_> = draft
        .fixtures
        .iter()
        .filter(|fixture| fixture.partition == Partition::HeldOut)
        .map(|fixture| fixture.id)
        .collect();
    assert_eq!(
        heldout,
        BTreeSet::from([
            "i06-causal-impact-max-holdout",
            "i06-coupled-properties-max-holdout",
            "i06-decision-certifier-mutants-max-holdout",
            "i06-posterior-coverage-core-holdout",
            "i06-provenance-adversaries-core-holdout",
            "i06-substitution-impact-core-holdout",
            "i06-supplier-drift-core-holdout",
        ])
    );
    let waiver = draft.waivers[0];
    assert_eq!(waiver.subject, "i06-governed-industrial-lot-pack");
    assert!(waiver.reason.contains("synthetic material-shaped data"));
    assert!(waiver.predicate.contains("multi-supplier governed pack"));
    assert!(
        waiver
            .promotion_effect
            .contains("no real supplier/lot qualification")
    );

    let frozen = draft.freeze().expect("the I06 seed must freeze");
    assert_eq!(frozen.initiative(), "I06");
    assert_eq!(frozen.version(), 1);
    assert_eq!(frozen.claims().len(), 12);
    assert_eq!(frozen.fixtures().len(), 17);
    assert_eq!(frozen.obligations().len(), 7);
    assert_eq!(frozen.waivers().len(), 1);
}

#[test]
#[allow(clippy::too_many_lines)]
fn i06_obligation_map_is_once_only_complete_and_operational() {
    let draft = i06_draft();
    let expected: BTreeMap<&str, (CampaignTier, BTreeSet<&str>)> = BTreeMap::from([
        (
            "i06-passport-observation-ingestion",
            (
                CampaignTier::Core,
                BTreeSet::from([
                    "i06-passport-identity-custody",
                    "i06-property-observation-semantics",
                ]),
            ),
        ),
        (
            "i06-hierarchical-posterior",
            (
                CampaignTier::Core,
                BTreeSet::from(["i06-hierarchical-lot-posterior"]),
            ),
        ),
        (
            "i06-substitution-impact",
            (
                CampaignTier::Core,
                BTreeSet::from([
                    "i06-contextual-substitution-admission",
                    "i06-selective-impact-invalidation",
                ]),
            ),
        ),
        (
            "i06-drift-decision-replay",
            (
                CampaignTier::Core,
                BTreeSet::from([
                    "i06-anytime-valid-supplier-drift",
                    "i06-replayable-substitution-decision",
                ]),
            ),
        ),
        (
            "i06-coupled-property-model",
            (
                CampaignTier::Max,
                BTreeSet::from(["i06-coupled-property-posterior"]),
            ),
        ),
        (
            "i06-causal-impact-certifier",
            (
                CampaignTier::Max,
                BTreeSet::from([
                    "i06-identifiable-causal-transport",
                    "i06-transitive-impact-completeness",
                ]),
            ),
        ),
        (
            "i06-robust-decision-certifier",
            (
                CampaignTier::Max,
                BTreeSet::from([
                    "i06-false-provenance-impact-certificate-falsifier",
                    "i06-robust-decision-optimality",
                ]),
            ),
        ),
    ]);
    let mut seen = BTreeMap::<&str, usize>::new();
    for row in &draft.obligations {
        let (tier, claims) = expected
            .get(row.leaf)
            .unwrap_or_else(|| panic!("unexpected I06 leaf '{}'", row.leaf));
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
        assert!(row.entry_point.starts_with("scripts/e2e/leapfrog/i06_"));
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
        let authored = i06_draft()
            .obligations
            .into_iter()
            .find(|candidate| candidate.leaf == row.leaf())
            .expect("authored row");
        assert_eq!(row.digest(), obligation_digest(&authored));
    }
}

#[test]
fn i06_property_semantics_prevent_material_quantity_confusion() {
    let draft = i06_draft();
    let property = claim(&draft, "i06-property-observation-semantics");
    for token in [
        "property kind",
        "censoring",
        "detection-limit",
        "uncertainty/covariance",
        "specimen geometry/location/orientation",
        "temperature",
        "pressure",
        "frequency",
        "strain/shear rate",
        "humidity",
        "surface/process/history state",
        "validity domain",
    ] {
        assert!(
            property.statement.contains(token),
            "property claim omits {token}"
        );
    }
    let hypothesis_text = property.hypotheses.join(" ");
    for distinction in [
        "electrical/thermal conductivity",
        "intrinsic permeability/gas permeability/permeance/diffusivity",
        "magnetic moment/magnetization/permeability/coercivity",
        "dynamic/kinematic/non-Newtonian viscosity",
        "latent heat per mass/volume",
        "ductility and hardness scale",
        "advancing/receding/static contact angle",
    ] {
        assert!(
            hypothesis_text.contains(distinction),
            "missing {distinction}"
        );
    }
    assert!(property.fallback.contains("do not guess a unit"));
    assert!(property.no_claim.contains("outside its validity domain"));

    let ontology = authored_spec(&draft, "i06-property-ontology-v1");
    for token in [
        "Rockwell/Vickers/Brinell hardness as distinct",
        "intrinsic porous permeability, gas permeability, permeance, diffusivity",
        "latent heat per mass/volume",
        "advancing/receding/static contact angle",
        "forbidden cross-kind conversion corpus",
    ] {
        assert!(ontology.contains(token), "ontology omits {token}");
    }
}

#[test]
fn i06_identity_posterior_substitution_impact_and_drift_no_claims_are_exact() {
    let draft = i06_draft();
    let identity = claim(&draft, "i06-passport-identity-custody");
    assert!(
        identity
            .hypotheses
            .iter()
            .any(|h| h.contains("signatures authenticate bytes"))
    );
    for boundary in [
        "physical composition",
        "absence of counterfeiting",
        "legal title",
        "supplier qualification",
    ] {
        assert!(identity.no_claim.contains(boundary));
    }

    let posterior = claim(&draft, "i06-hierarchical-lot-posterior");
    assert_eq!(
        posterior.tolerance,
        ToleranceSemantics::Absolute { atol: 0.03 }
    );
    assert!(posterior.statement.contains("without pseudoreplication"));
    for owner in ["aleatory", "epistemic", "measurement", "model discrepancy"] {
        assert!(posterior.hypotheses.iter().any(|h| h.contains(owner)));
    }

    let substitution = claim(&draft, "i06-contextual-substitution-admission");
    for state in ["Compatible", "Incompatible", "Unknown"] {
        assert!(substitution.statement.contains(state));
    }
    assert!(
        substitution
            .hypotheses
            .iter()
            .any(|h| h.contains("grade/spec-name equality"))
    );
    assert!(substitution.no_claim.contains("not supplier approval"));

    let impact = claim(&draft, "i06-selective-impact-invalidation");
    assert!(
        impact
            .hypotheses
            .iter()
            .any(|h| h.contains("relative to this declared graph"))
    );
    assert!(impact.no_claim.contains("not proof that every real causal"));

    let drift = claim(&draft, "i06-anytime-valid-supplier-drift");
    assert_eq!(
        drift.tolerance,
        ToleranceSemantics::Interval { lo: 0.0, hi: 0.05 }
    );
    for state in ["Alarm", "NoAlarm", "DataInvalid", "ModelUnknown"] {
        assert!(drift.statement.contains(state));
    }
    assert!(
        drift
            .no_claim
            .contains("NoAlarm proves neither no physical drift")
    );
}

#[test]
fn i06_policy_is_the_authority_and_retention_spine() {
    let draft = i06_draft();
    let policy = authored_spec(&draft, POLICY);
    assert_eq!(policy.lines().next(), Some("I06_CAMPAIGN_POLICY_V1"));
    for heading in [
        "IDENTITY=",
        "AUTHENTICITY=",
        "PROPERTY=",
        "UNCERTAINTY=",
        "SUBSTITUTION=",
        "IMPACT=",
        "DRIFT=",
        "CAUSAL=",
        "DECISION=",
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
    assert!(policy.contains("valid signature proves byte/key binding only"));
    assert!(policy.contains("grade/spec-name equality or nominal overlap grants no equivalence"));
    assert!(policy.contains("exact selective invalidation is relative to the graph"));
    assert!(policy.contains("NoAlarm never proves no drift"));
    assert!(policy.contains("learned association is not a causal graph"));
    assert!(policy.contains("no automatic procurement/manufacturing action"));
    assert!(policy.contains("axes never substitute"));
    assert!(policy.contains("partial success cannot publish normal authority"));
}

#[test]
fn i06_holdouts_have_disjoint_ranges_and_one_stage_local_consumer() {
    let draft = i06_draft();
    for token in [
        "development indices 0..=16383",
        "core held-out 65536..=81919",
        "maximal held-out 131072..=147455",
    ] {
        assert!(
            draft.explicits.seeds.contains(token),
            "seed policy omits {token}"
        );
    }
    let expected = [
        (
            "i06-provenance-adversaries-core-holdout",
            "65536..=69631",
            "i06-passport-observation-ingestion",
            CampaignTier::Core,
        ),
        (
            "i06-posterior-coverage-core-holdout",
            "69632..=73727",
            "i06-hierarchical-posterior",
            CampaignTier::Core,
        ),
        (
            "i06-substitution-impact-core-holdout",
            "73728..=77823",
            "i06-substitution-impact",
            CampaignTier::Core,
        ),
        (
            "i06-supplier-drift-core-holdout",
            "77824..=81919",
            "i06-drift-decision-replay",
            CampaignTier::Core,
        ),
        (
            "i06-coupled-properties-max-holdout",
            "131072..=135167",
            "i06-coupled-property-model",
            CampaignTier::Max,
        ),
        (
            "i06-causal-impact-max-holdout",
            "135168..=139263",
            "i06-causal-impact-certifier",
            CampaignTier::Max,
        ),
        (
            "i06-decision-certifier-mutants-max-holdout",
            "139264..=147455",
            "i06-robust-decision-certifier",
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
#[allow(clippy::too_many_lines)]
fn i06_maximal_causal_completeness_and_decision_ratchets_mint_no_prose_authority() {
    let draft = i06_draft();
    let causal = claim(&draft, "i06-identifiable-causal-transport");
    assert!(causal.statement.contains("only identifiable"));
    assert!(causal.statement.contains("CausalUnknown"));
    assert!(causal.hypotheses.iter().any(|h| h.contains("SUTVA")));
    assert!(
        causal
            .hypotheses
            .iter()
            .any(|h| h.contains("graph discovery output"))
    );
    assert!(causal.no_claim.contains("does not validate that graph"));

    let completeness = claim(&draft, "i06-transitive-impact-completeness");
    assert!(completeness.activation.contains("pre-proof successor"));
    assert!(
        completeness
            .hypotheses
            .iter()
            .any(|h| h.contains("proposition/definition ASTs"))
    );
    assert!(
        completeness
            .hypotheses
            .iter()
            .any(|h| h.contains("open-world boundaries"))
    );
    assert!(
        completeness
            .no_claim
            .contains("version-1 prose mints no theorem")
    );

    let optimal = claim(&draft, "i06-robust-decision-optimality");
    assert_eq!(
        optimal.tolerance,
        ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 }
    );
    assert!(
        optimal
            .hypotheses
            .iter()
            .any(|h| h.contains("globality is only over"))
    );
    assert!(
        optimal
            .kill
            .contains("independently better admissible policy")
    );
    assert!(optimal.no_claim.contains("misspecified utility"));

    let theorem_card = authored_spec(&draft, "i06-impact-causal-theorem-card");
    for token in [
        "graph/SCM/authority/proposition/definition ASTs",
        "complete adapters/open-world boundaries",
        "total premise map",
        "deterministic Lean translation/roundtrip",
        "exact axiom allowlist {propext,Quot.sound,Classical.choice}",
        "transitive closure",
        "nonvacuity witnesses",
        "retained kernel replay",
    ] {
        assert!(theorem_card.contains(token), "theorem card omits {token}");
    }
    assert!(theorem_card.contains("Prose grants no theorem"));

    let decision_card = authored_spec(&draft, "i06-decision-grammar-card");
    for token in [
        "finite incumbent/candidate",
        "feasibility and interval utility semantics",
        "validity and exclusion order",
        "exact enumeration or verified relaxation bounds",
        "rank/unrank/sharding",
        "independent decoder/checker",
        "preflight and completeness root",
    ] {
        assert!(decision_card.contains(token), "decision card omits {token}");
    }
    assert!(decision_card.contains("no global-optimality or decision authority"));

    let falsifier = claim(&draft, "i06-false-provenance-impact-certificate-falsifier");
    assert_eq!(falsifier.polarity, ClaimPolarity::Refutation);
    assert!(falsifier.kill.contains("intended lane succeeds"));
    assert!(
        falsifier
            .no_claim
            .contains("finite adversarial corpus cannot prove")
    );

    let robust = draft
        .obligations
        .iter()
        .find(|row| row.leaf == "i06-robust-decision-certifier")
        .expect("robust row");
    for deck in [
        POLICY,
        "i06-decision-grammar-card",
        "i06-decision-certifier-mutants-max-holdout",
        "i06-governed-industrial-lot-pack",
    ] {
        assert!(robust.decks.contains(&deck), "robust row omits {deck}");
    }
    assert!(robust.g4_schedule.contains("whole-campaign preflight"));
    assert!(robust.g4_schedule.contains("BudgetExhausted stays Unknown"));
}

#[test]
fn i06_g3_mutations_refuse_or_move_authority() {
    let baseline = i06_draft().freeze().expect("freeze").digest();

    let mut missing_hypotheses = i06_draft();
    missing_hypotheses
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i06-hierarchical-lot-posterior")
        .expect("posterior claim")
        .hypotheses = &[];
    assert!(matches!(
        missing_hypotheses.freeze(),
        Err(FreezeRefusal::BlankField {
            field: "claim.hypotheses",
            ..
        })
    ));

    let mut correlated = i06_draft();
    correlated
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i06-identifiable-causal-transport")
        .expect("causal claim")
        .oracle
        .independent = false;
    assert!(matches!(
        correlated.freeze(),
        Err(FreezeRefusal::ProductionOracleReuse { .. })
    ));

    let mut relaxed = i06_draft();
    relaxed
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i06-hierarchical-lot-posterior")
        .expect("posterior claim")
        .tolerance = ToleranceSemantics::Absolute { atol: 0.10 };
    assert_ne!(relaxed.freeze().expect("relaxed freeze").digest(), baseline);

    let mut swapped = i06_draft();
    swapped
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i06-substitution-impact-core-holdout")
        .expect("substitution holdout")
        .source = FixtureSource::AuthoredSpec {
        spec: "unauthorized post-result replacement",
    };
    assert_ne!(
        swapped.freeze().expect("replacement freeze").digest(),
        baseline
    );

    let mut repartitioned = i06_draft();
    repartitioned
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i06-causal-impact-max-holdout")
        .expect("causal holdout")
        .partition = Partition::Development;
    assert_ne!(
        repartitioned.freeze().expect("repartition freeze").digest(),
        baseline
    );

    let mut missing_policy = i06_draft();
    missing_policy
        .fixtures
        .retain(|fixture| fixture.id != POLICY);
    assert!(matches!(
        missing_policy.freeze(),
        Err(FreezeRefusal::OrphanDeck { deck, .. }) if deck == POLICY
    ));
}

#[test]
fn i06_g5_top_level_order_is_not_identity() {
    let expected = i06_draft().freeze().expect("freeze");
    let mut permuted = i06_draft();
    permuted.claims.reverse();
    permuted.fixtures.reverse();
    permuted.obligations.reverse();
    permuted.waivers.reverse();
    let actual = permuted.freeze().expect("permuted freeze");
    assert_eq!(actual.digest(), expected.digest());
    assert_eq!(actual, expected);
}

#[test]
fn i06_g4_chunked_in_memory_assembly_is_identity_equivalent() {
    let one_shot = i06_draft();
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
fn i06_amendments_invalidate_exact_targeted_or_global_authority() {
    let predecessor = i06_draft();
    let all = authority_ids(&predecessor);
    let frozen = predecessor.freeze().expect("freeze");

    let mut version_only = i06_draft();
    version_only.version = 2;
    let (_, version_record) = frozen.amend(version_only).expect("version-only");
    assert!(version_record.invalidated.is_empty());

    let mut property = i06_draft();
    property.version = 2;
    property
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i06-property-observation-semantics")
        .expect("property claim")
        .statement = "successor property observation authority";
    let (_, property_record) = frozen.amend(property).expect("property amendment");
    assert_eq!(
        property_record.invalidated,
        vec![
            "i06-passport-observation-ingestion",
            "i06-property-observation-semantics",
        ]
    );

    let mut holdout = i06_draft();
    holdout.version = 2;
    holdout
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i06-substitution-impact-core-holdout")
        .expect("holdout")
        .source = FixtureSource::AuthoredSpec {
        spec: "successor substitution/impact holdout",
    };
    let (_, holdout_record) = frozen.amend(holdout).expect("holdout amendment");
    assert_eq!(
        holdout_record.invalidated,
        vec![
            "i06-contextual-substitution-admission",
            "i06-selective-impact-invalidation",
            "i06-substitution-impact",
        ]
    );

    let mut policy = i06_draft();
    policy.version = 2;
    policy
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == POLICY)
        .expect("policy")
        .source = FixtureSource::AuthoredSpec {
        spec: "I06_CAMPAIGN_POLICY_V2 changed global authority",
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

    let mut title = i06_draft();
    title.version = 2;
    title.title = "successor global I06 campaign authority";
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
