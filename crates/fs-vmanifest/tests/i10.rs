//! Focused I10 constitutive-identifiability/coupon-design manifest
//! conformance.
//!
//! The tests lock preregistration authority only. They do not identify a
//! parameter, verify an adjoint, prove a quotient theorem, admit a coupon for
//! fabrication, discriminate rival laws, or optimize a real experiment
//! campaign.

use fs_vmanifest::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FixtureSource, FreezeRefusal, ManifestDraft,
    Partition, ToleranceSemantics, i10_draft, obligation_digest,
};
use std::collections::{BTreeMap, BTreeSet};

const POLICY: &str = "i10-campaign-policy-v1";
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
        .unwrap_or_else(|| panic!("missing I10 claim '{id}'"))
}

fn authored_spec<'a>(draft: &'a ManifestDraft, id: &str) -> &'a str {
    let fixture = draft
        .fixtures
        .iter()
        .find(|fixture| fixture.id == id)
        .unwrap_or_else(|| panic!("missing I10 fixture '{id}'"));
    match fixture.source {
        FixtureSource::AuthoredSpec { spec } => spec,
        FixtureSource::External { .. } => panic!("I10 fixture '{id}' must be authored"),
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
fn i10_seed_freezes_with_exact_lattice_corpus_and_waiver() {
    let draft = i10_draft();
    assert_eq!(draft.initiative, "I10");
    assert_eq!(draft.version, 1);
    assert_eq!(draft.claims.len(), 11);
    assert_eq!(draft.fixtures.len(), 14);
    assert_eq!(draft.obligations.len(), 6);
    assert_eq!(draft.waivers.len(), 1);

    let lattice = draft.claims.iter().fold([0usize; 3], |mut counts, claim| {
        counts[match claim.ambition {
            Ambition::Solid => 0,
            Ambition::Frontier => 1,
            Ambition::Moonshot => 2,
        }] += 1;
        counts
    });
    assert_eq!(lattice, [6, 2, 3]);
    let refutations: Vec<_> = draft
        .claims
        .iter()
        .filter(|claim| claim.polarity == ClaimPolarity::Refutation)
        .map(|claim| claim.id)
        .collect();
    assert_eq!(
        refutations,
        ["i10-false-identifiability-design-certificate-falsifier"]
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
            "i10-blind-experiment-gain-core-holdout",
            "i10-design-certifier-mutants-max-holdout",
            "i10-discrimination-max-holdout",
            "i10-gauge-quotient-max-holdout",
            "i10-identifiability-core-holdout",
            "i10-schema-adjoint-adversaries-core-holdout",
        ])
    );
    let waiver = draft.waivers[0];
    assert_eq!(waiver.subject, "i10-instrumented-lab-campaign-pack");
    assert!(waiver.reason.contains("synthetic coupon-shaped data"));
    assert!(waiver.predicate.contains("governed instrumented-lab pack"));
    assert!(
        waiver
            .promotion_effect
            .contains("no physical-lab identifiability")
    );

    let frozen = draft.freeze().expect("the I10 seed must freeze");
    assert_eq!(frozen.initiative(), "I10");
    assert_eq!(frozen.version(), 1);
    assert_eq!(frozen.claims().len(), 11);
    assert_eq!(frozen.fixtures().len(), 14);
    assert_eq!(frozen.obligations().len(), 6);
    assert_eq!(frozen.waivers().len(), 1);
}

#[test]
#[allow(clippy::too_many_lines)]
fn i10_obligation_map_is_once_only_complete_and_operational() {
    let draft = i10_draft();
    let expected: BTreeMap<&str, (CampaignTier, BTreeSet<&str>)> = BTreeMap::from([
        (
            "i10-law-schema-adjoint",
            (
                CampaignTier::Core,
                BTreeSet::from([
                    "i10-law-experiment-schema-identity",
                    "i10-sensitivity-adjoint-consistency",
                ]),
            ),
        ),
        (
            "i10-identifiability-analysis",
            (
                CampaignTier::Core,
                BTreeSet::from([
                    "i10-structural-identifiability-verdicts",
                    "i10-practical-identifiability-sloppiness",
                ]),
            ),
        ),
        (
            "i10-coupon-design-gain",
            (
                CampaignTier::Core,
                BTreeSet::from([
                    "i10-manufacturable-coupon-admission",
                    "i10-heldout-design-gain",
                ]),
            ),
        ),
        (
            "i10-robust-discrimination-designer",
            (
                CampaignTier::Max,
                BTreeSet::from([
                    "i10-discrepancy-aware-robust-design",
                    "i10-anytime-valid-adaptive-discrimination",
                ]),
            ),
        ),
        (
            "i10-quotient-theorem-certifier",
            (
                CampaignTier::Max,
                BTreeSet::from(["i10-gauge-quotient-completeness-theorem"]),
            ),
        ),
        (
            "i10-global-design-certifier",
            (
                CampaignTier::Max,
                BTreeSet::from([
                    "i10-global-design-optimality",
                    "i10-false-identifiability-design-certificate-falsifier",
                ]),
            ),
        ),
    ]);
    let mut seen = BTreeMap::<&str, usize>::new();
    for row in &draft.obligations {
        let (tier, claims) = expected
            .get(row.leaf)
            .unwrap_or_else(|| panic!("unexpected I10 leaf '{}'", row.leaf));
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
        assert!(row.entry_point.starts_with("scripts/e2e/leapfrog/i10_"));
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
        let authored = i10_draft()
            .obligations
            .into_iter()
            .find(|candidate| candidate.leaf == row.leaf())
            .expect("authored row");
        assert_eq!(row.digest(), obligation_digest(&authored));
    }
}

#[test]
fn i10_gauge_semantics_prevent_false_identifiability() {
    let draft = i10_draft();
    let structural = claim(&draft, "i10-structural-identifiability-verdicts");
    for state in ["Identifiable", "NonIdentifiable", "Unknown"] {
        assert!(
            structural.statement.contains(state),
            "structural claim omits {state}"
        );
    }
    assert!(structural.statement.contains("gauge witness"));
    assert!(
        structural
            .statement
            .contains("never mints Identifiable from numerics alone")
    );
    let hypothesis_text = structural.hypotheses.join(" ");
    for confounding in [
        "kinematic and isotropic hardening under monotone proportional loading",
        "Prony relabeling",
        "modulus-thickness products",
        "gauge-equivalent parameter sets are one physical authority",
        "full-rank numeric FIM alone",
    ] {
        assert!(
            hypothesis_text.contains(confounding),
            "missing {confounding}"
        );
    }
    assert!(structural.kill.contains("known gauge family"));

    let practical = claim(&draft, "i10-practical-identifiability-sloppiness");
    assert_eq!(
        practical.tolerance,
        ToleranceSemantics::Absolute { atol: 0.05 }
    );
    assert!(
        practical
            .hypotheses
            .iter()
            .any(|h| h.contains("not a certificate of structural non-identifiability"))
    );
    assert!(
        practical
            .statement
            .contains("prior/constraint-induced stiffness")
    );

    let counterexamples = authored_spec(&draft, "i10-gauge-symmetry-counterexamples");
    for token in [
        "Ogden exponent/modulus permutations",
        "coincident time constants",
        "modulus-thickness and modulus-area products",
        "split under cyclic or non-proportional paths",
        "unit-rescaling gauges",
        "independent invariance witness",
    ] {
        assert!(counterexamples.contains(token), "corpus omits {token}");
    }

    let catalog = authored_spec(&draft, "i10-constitutive-law-catalog");
    for family in [
        "Ogden hyperelasticity",
        "Prony-series linear viscoelasticity",
        "Voce isotropic",
        "Armstrong-Frederick kinematic",
        "Norton-Bailey creep",
        "kinematic-versus-isotropic monotone confounding",
    ] {
        assert!(catalog.contains(family), "catalog omits {family}");
    }
}

#[test]
fn i10_sensitivity_design_and_gain_no_claims_are_exact() {
    let draft = i10_draft();
    let adjoint = claim(&draft, "i10-sensitivity-adjoint-consistency");
    assert_eq!(adjoint.tolerance, ToleranceSemantics::Exact);
    assert!(adjoint.statement.contains("outward-rounded"));
    assert!(
        adjoint
            .hypotheses
            .iter()
            .any(|h| h.contains("refused, not zero-filled"))
    );
    assert!(adjoint.no_claim.contains("not identifiability"));

    let admission = claim(&draft, "i10-manufacturable-coupon-admission");
    for state in ["Feasible", "Infeasible", "Unknown"] {
        assert!(admission.statement.contains(state));
    }
    assert!(
        admission
            .statement
            .contains("feasibility precedes any information criterion")
    );
    assert!(
        admission
            .hypotheses
            .iter()
            .any(|h| h.contains("however favorable"))
    );
    assert!(admission.no_claim.contains("not fabrication approval"));

    let gain = claim(&draft, "i10-heldout-design-gain");
    assert_eq!(
        gain.tolerance,
        ToleranceSemantics::Interval { lo: 0.0, hi: 0.05 }
    );
    assert!(
        gain.hypotheses
            .iter()
            .any(|h| h.contains("unidentified directions report unchanged-by-construction"))
    );
    assert!(gain.no_claim.contains("not real-laboratory gain"));

    let discrimination = claim(&draft, "i10-anytime-valid-adaptive-discrimination");
    assert_eq!(
        discrimination.tolerance,
        ToleranceSemantics::Interval { lo: 0.0, hi: 0.05 }
    );
    for state in ["Selected", "Undecided", "DataInvalid", "ModelUnknown"] {
        assert!(discrimination.statement.contains(state));
    }
    assert!(
        discrimination
            .no_claim
            .contains("Undecided never proves model equivalence")
    );

    let robust = claim(&draft, "i10-discrepancy-aware-robust-design");
    assert!(
        robust
            .hypotheses
            .iter()
            .any(|h| h.contains("missing physics never contributes a favorable point value"))
    );
    assert!(robust.no_claim.contains("unmodeled physics"));
}

#[test]
fn i10_policy_is_the_authority_and_retention_spine() {
    let draft = i10_draft();
    let policy = authored_spec(&draft, POLICY);
    assert_eq!(policy.lines().next(), Some("I10_CAMPAIGN_POLICY_V1"));
    for heading in [
        "IDENTITY=",
        "LAW=",
        "SENSITIVITY=",
        "IDENTIFIABILITY=",
        "SLOPPINESS=",
        "DESIGN=",
        "GAIN=",
        "DISCRIMINATION=",
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
    assert!(policy.contains("a full-rank numeric FIM alone never mints Identifiable"));
    assert!(policy.contains("gauge-equivalent parameter sets are one physical authority"));
    assert!(policy.contains("design feasibility precedes information optimality"));
    assert!(policy.contains("no automatic lab actuation"));
    assert!(policy.contains("unidentified directions never count toward gain"));
    assert!(policy.contains("Undecided never proves model equivalence"));
    assert!(policy.contains("axes never substitute"));
    assert!(policy.contains("partial success cannot publish normal authority"));
}

#[test]
fn i10_holdouts_have_disjoint_ranges_and_one_stage_local_consumer() {
    let draft = i10_draft();
    for token in [
        "development indices 0..=16383",
        "core held-out 65536..=77823",
        "maximal held-out 131072..=147455",
    ] {
        assert!(
            draft.explicits.seeds.contains(token),
            "seed policy omits {token}"
        );
    }
    let expected = [
        (
            "i10-schema-adjoint-adversaries-core-holdout",
            "65536..=69631",
            "i10-law-schema-adjoint",
            CampaignTier::Core,
        ),
        (
            "i10-identifiability-core-holdout",
            "69632..=73727",
            "i10-identifiability-analysis",
            CampaignTier::Core,
        ),
        (
            "i10-blind-experiment-gain-core-holdout",
            "73728..=77823",
            "i10-coupon-design-gain",
            CampaignTier::Core,
        ),
        (
            "i10-discrimination-max-holdout",
            "131072..=135167",
            "i10-robust-discrimination-designer",
            CampaignTier::Max,
        ),
        (
            "i10-gauge-quotient-max-holdout",
            "135168..=139263",
            "i10-quotient-theorem-certifier",
            CampaignTier::Max,
        ),
        (
            "i10-design-certifier-mutants-max-holdout",
            "139264..=147455",
            "i10-global-design-certifier",
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
fn i10_maximal_theorem_and_design_ratchets_mint_no_prose_authority() {
    let draft = i10_draft();
    let quotient = claim(&draft, "i10-gauge-quotient-completeness-theorem");
    assert!(quotient.statement.contains("complete"));
    assert!(quotient.statement.contains("independent"));
    assert!(quotient.activation.contains("pre-proof successor"));
    assert!(
        quotient
            .hypotheses
            .iter()
            .any(|h| h.contains("law/observable/group ASTs"))
    );
    assert!(
        quotient
            .no_claim
            .contains("version-1 prose mints no theorem authority")
    );

    let optimal = claim(&draft, "i10-global-design-optimality");
    assert_eq!(
        optimal.tolerance,
        ToleranceSemantics::Interval { lo: 0.0, hi: 0.0 }
    );
    assert!(
        optimal
            .hypotheses
            .iter()
            .any(|h| h.contains("pre-candidate successor"))
    );
    assert!(
        optimal
            .hypotheses
            .iter()
            .any(|h| h.contains("cannot hide rounding"))
    );
    assert!(
        optimal
            .kill
            .contains("independently better admissible design")
    );
    assert!(optimal.no_claim.contains("finite encoded grammar"));

    let theorem_card = authored_spec(&draft, "i10-quotient-theorem-card");
    for token in [
        "law/observable/group ASTs and canonical bytes",
        "total premise map",
        "deterministic Lean translation/roundtrip",
        "exact axiom allowlist {propext,Quot.sound,Classical.choice}",
        "transitive closure",
        "nonvacuity witnesses and retained kernel replay",
    ] {
        assert!(theorem_card.contains(token), "theorem card omits {token}");
    }
    assert!(theorem_card.contains("Prose grants no theorem"));

    let grammar_card = authored_spec(&draft, "i10-design-grammar-card");
    for token in [
        "finite",
        "feasibility and interval information-criterion semantics",
        "validity and exclusion order",
        "exact enumeration or verified relaxation bounds",
        "rank/unrank/sharding",
        "independent decoder/checker",
        "preflight and completeness root",
    ] {
        assert!(grammar_card.contains(token), "grammar card omits {token}");
    }
    assert!(grammar_card.contains("no global-optimality or design authority"));

    let falsifier = claim(
        &draft,
        "i10-false-identifiability-design-certificate-falsifier",
    );
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
        .find(|row| row.leaf == "i10-global-design-certifier")
        .expect("certifier row");
    for deck in [
        POLICY,
        "i10-design-grammar-card",
        "i10-design-certifier-mutants-max-holdout",
        "i10-instrumented-lab-campaign-pack",
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
fn i10_g3_mutations_refuse_or_move_authority() {
    let baseline = i10_draft().freeze().expect("freeze").digest();

    let mut missing_hypotheses = i10_draft();
    missing_hypotheses
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i10-practical-identifiability-sloppiness")
        .expect("practical claim")
        .hypotheses = &[];
    assert!(matches!(
        missing_hypotheses.freeze(),
        Err(FreezeRefusal::BlankField {
            field: "claim.hypotheses",
            ..
        })
    ));

    let mut correlated = i10_draft();
    correlated
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i10-anytime-valid-adaptive-discrimination")
        .expect("discrimination claim")
        .oracle
        .independent = false;
    assert!(matches!(
        correlated.freeze(),
        Err(FreezeRefusal::ProductionOracleReuse { .. })
    ));

    let mut relaxed = i10_draft();
    relaxed
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i10-practical-identifiability-sloppiness")
        .expect("practical claim")
        .tolerance = ToleranceSemantics::Absolute { atol: 0.20 };
    assert_ne!(relaxed.freeze().expect("relaxed freeze").digest(), baseline);

    let mut swapped = i10_draft();
    swapped
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i10-identifiability-core-holdout")
        .expect("identifiability holdout")
        .source = FixtureSource::AuthoredSpec {
        spec: "unauthorized post-result replacement",
    };
    assert_ne!(
        swapped.freeze().expect("replacement freeze").digest(),
        baseline
    );

    let mut repartitioned = i10_draft();
    repartitioned
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i10-gauge-quotient-max-holdout")
        .expect("quotient holdout")
        .partition = Partition::Development;
    assert_ne!(
        repartitioned.freeze().expect("repartition freeze").digest(),
        baseline
    );

    let mut missing_policy = i10_draft();
    missing_policy
        .fixtures
        .retain(|fixture| fixture.id != POLICY);
    assert!(matches!(
        missing_policy.freeze(),
        Err(FreezeRefusal::OrphanDeck { deck, .. }) if deck == POLICY
    ));
}

#[test]
fn i10_g5_top_level_order_is_not_identity() {
    let expected = i10_draft().freeze().expect("freeze");
    let mut permuted = i10_draft();
    permuted.claims.reverse();
    permuted.fixtures.reverse();
    permuted.obligations.reverse();
    permuted.waivers.reverse();
    let actual = permuted.freeze().expect("permuted freeze");
    assert_eq!(actual.digest(), expected.digest());
    assert_eq!(actual, expected);
}

#[test]
fn i10_g4_chunked_in_memory_assembly_is_identity_equivalent() {
    let one_shot = i10_draft();
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
fn i10_amendments_invalidate_exact_targeted_or_global_authority() {
    let predecessor = i10_draft();
    let all = authority_ids(&predecessor);
    let frozen = predecessor.freeze().expect("freeze");

    let mut version_only = i10_draft();
    version_only.version = 2;
    let (_, version_record) = frozen.amend(version_only).expect("version-only");
    assert!(version_record.invalidated.is_empty());

    let mut adjoint = i10_draft();
    adjoint.version = 2;
    adjoint
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i10-sensitivity-adjoint-consistency")
        .expect("adjoint claim")
        .statement = "successor adjoint sensitivity authority";
    let (_, adjoint_record) = frozen.amend(adjoint).expect("adjoint amendment");
    assert_eq!(
        adjoint_record.invalidated,
        vec![
            "i10-law-schema-adjoint",
            "i10-sensitivity-adjoint-consistency",
        ]
    );

    let mut holdout = i10_draft();
    holdout.version = 2;
    holdout
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == "i10-identifiability-core-holdout")
        .expect("holdout")
        .source = FixtureSource::AuthoredSpec {
        spec: "successor identifiability holdout",
    };
    let (_, holdout_record) = frozen.amend(holdout).expect("holdout amendment");
    assert_eq!(
        holdout_record.invalidated,
        vec![
            "i10-identifiability-analysis",
            "i10-practical-identifiability-sloppiness",
            "i10-structural-identifiability-verdicts",
        ]
    );

    let mut policy = i10_draft();
    policy.version = 2;
    policy
        .fixtures
        .iter_mut()
        .find(|fixture| fixture.id == POLICY)
        .expect("policy")
        .source = FixtureSource::AuthoredSpec {
        spec: "I10_CAMPAIGN_POLICY_V2 changed global authority",
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

    let mut title = i10_draft();
    title.version = 2;
    title.title = "successor global I10 campaign authority";
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
