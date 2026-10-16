//! Focused I07 AP242 semantic manufacturing round-trip conformance.
//!
//! These tests pin the exact claim lattice, profile and authority boundaries,
//! stage-separated fixtures, execution contracts, mutation sensitivity,
//! canonical assembly, and transitive amendment invalidation. Generic schema
//! gate ordering remains in `tests/conformance.rs`.

use fs_vmanifest::{
    Ambition, CampaignTier, ClaimPolarity, FixturePin, FixtureSource, FreezeRefusal, GauntletTier,
    ManifestDraft, Partition, ToleranceSemantics, i07_draft,
};
use std::collections::{BTreeMap, BTreeSet};

const CAMPAIGN_POLICY_FIXTURE: &str = "i07-campaign-policy-v1";
const PROFILE_REGISTRY_FIXTURE: &str = "i07-profile-registry-v1";
const REAL_CORPUS_SLOT: &str = "i07-governed-real-ap242-corpus";

const CLAIM_AUTHORITY: &[(&str, Ambition, ClaimPolarity, GauntletTier)] = &[
    (
        "i07-part21-profile-bounded-parse",
        Ambition::Solid,
        ClaimPolarity::Affirmative,
        GauntletTier::G0,
    ),
    (
        "i07-product-assembly-configuration-lineage",
        Ambition::Solid,
        ClaimPolarity::Affirmative,
        GauntletTier::G3,
    ),
    (
        "i07-representation-qualified-geometry-roundtrip",
        Ambition::Solid,
        ClaimPolarity::Affirmative,
        GauntletTier::G2,
    ),
    (
        "i07-semantic-pmi-gdt-datum-texture",
        Ambition::Solid,
        ClaimPolarity::Affirmative,
        GauntletTier::G3,
    ),
    (
        "i07-material-process-lot-requirement-evidence",
        Ambition::Solid,
        ClaimPolarity::Affirmative,
        GauntletTier::G3,
    ),
    (
        "i07-harness-ewis-connectivity-semantics",
        Ambition::Solid,
        ClaimPolarity::Affirmative,
        GauntletTier::G3,
    ),
    (
        "i07-exhaustive-semantic-loss-receipt",
        Ambition::Solid,
        ClaimPolarity::Affirmative,
        GauntletTier::G3,
    ),
    (
        "i07-canonical-export-determinism-and-containment",
        Ambition::Solid,
        ClaimPolarity::Affirmative,
        GauntletTier::G5,
    ),
    (
        "i07-safe-opaque-extension-capsules",
        Ambition::Frontier,
        ClaimPolarity::Affirmative,
        GauntletTier::G3,
    ),
    (
        "i07-kinematic-pair-semantic-roundtrip",
        Ambition::Frontier,
        ClaimPolarity::Affirmative,
        GauntletTier::G3,
    ),
    (
        "i07-partial-bidirectional-edition-migration",
        Ambition::Frontier,
        ClaimPolarity::Affirmative,
        GauntletTier::G3,
    ),
    (
        "i07-governed-real-ap242-corpus-validation",
        Ambition::Frontier,
        ClaimPolarity::Affirmative,
        GauntletTier::G2,
    ),
    (
        "i07-proof-carrying-bidirectional-semantic-equivalence",
        Ambition::Moonshot,
        ClaimPolarity::Affirmative,
        GauntletTier::G3,
    ),
    (
        "i07-semantic-sheaf-descent-obstruction-theorem",
        Ambition::Moonshot,
        ClaimPolarity::Affirmative,
        GauntletTier::G3,
    ),
    (
        "i07-semantic-equivalence-counterexample-search",
        Ambition::Moonshot,
        ClaimPolarity::Refutation,
        GauntletTier::G3,
    ),
];

const CLAIM_ORACLES: &[(&str, &str)] = &[
    (
        "i07-part21-profile-bounded-parse",
        "i07.oracle.part21.v1 at fs-vmanifest-oracles/i07/part21.rs::parse_generated_ast_and_refusals",
    ),
    (
        "i07-product-assembly-configuration-lineage",
        "i07.oracle.product_graph.v1 at fs-vmanifest-oracles/i07/product_graph.rs::check_typed_graph_isomorphism",
    ),
    (
        "i07-representation-qualified-geometry-roundtrip",
        "i07.oracle.geometry.v1 composite fs-vmanifest-oracles/i07/geometry.rs::enclose_roundtrip and independent interval subdivision/oriented-intersection/winding checkers",
    ),
    (
        "i07-semantic-pmi-gdt-datum-texture",
        "i07.oracle.pmi.v1 at fs-vmanifest-oracles/i07/pmi.rs::check_semantic_atoms_and_associations",
    ),
    (
        "i07-material-process-lot-requirement-evidence",
        "i07.oracle.material_lineage.v1 at fs-vmanifest-oracles/i07/material_lineage.rs::check_reference_graph",
    ),
    (
        "i07-harness-ewis-connectivity-semantics",
        "i07.oracle.harness.v1 at fs-vmanifest-oracles/i07/harness.rs::check_multigraph_isomorphism",
    ),
    (
        "i07-exhaustive-semantic-loss-receipt",
        "i07.oracle.semantic_loss.v1 at fs-vmanifest-oracles/i07/semantic_loss.rs::reconcile_independent_inventory",
    ),
    (
        "i07-canonical-export-determinism-and-containment",
        "i07.oracle.canonical_writer.v1 at fs-vmanifest-oracles/i07/canonical_writer.rs::independent_reparse_and_compare",
    ),
    (
        "i07-safe-opaque-extension-capsules",
        "i07.oracle.opaque_capsule.v1 at fs-vmanifest-oracles/i07/opaque_capsule.rs::check_bytes_tokens_refs_and_hazards",
    ),
    (
        "i07-kinematic-pair-semantic-roundtrip",
        "i07.oracle.kinematics.v1 at fs-vmanifest-oracles/i07/kinematics.rs::check_relation_and_tangent_graph",
    ),
    (
        "i07-partial-bidirectional-edition-migration",
        "i07.oracle.edition_migration.v1 at fs-vmanifest-oracles/i07/edition_migration.rs::replay_lens_and_crosswalk",
    ),
    (
        "i07-governed-real-ap242-corpus-validation",
        "i07.oracle.governed_real_corpus.v1 at fs-vmanifest-oracles/i07/governed_real_corpus.rs::check_custody_profiles_semantics_and_scoped_axes",
    ),
    (
        "i07-proof-carrying-bidirectional-semantic-equivalence",
        "i07.oracle.semantic_equivalence_lean.v1 at proofs/i07/SemanticEquivalence.lean::supportedRoundTripEquivalence checked by pinned Lean4 kernel receipt",
    ),
    (
        "i07-semantic-sheaf-descent-obstruction-theorem",
        "i07.oracle.semantic_descent_lean.v1 at proofs/i07/SemanticDescent.lean::finiteDescentAndObstruction checked by pinned Lean4 kernel receipt",
    ),
    (
        "i07-semantic-equivalence-counterexample-search",
        "i07.oracle.semantic_falsifier.v1 at fs-vmanifest-oracles/i07/semantic_falsifier.rs::verify_membership_cover_and_minimize",
    ),
];

const CLAIM_ACCEPTANCE: &[(&str, &str, &str, ToleranceSemantics)] = &[
    (
        "i07-part21-profile-bounded-parse",
        "exact_profile_token_entity_reference_and_refusal_verdict",
        "bit",
        ToleranceSemantics::Exact,
    ),
    (
        "i07-product-assembly-configuration-lineage",
        "exact_product_occurrence_configuration_lineage_graph_isomorphism_verdict",
        "bit",
        ToleranceSemantics::Exact,
    ),
    (
        "i07-representation-qualified-geometry-roundtrip",
        "maximum_preregistered_normalized_bidirectional_geometry_validation_and_topology_defect",
        "1",
        ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
    ),
    (
        "i07-semantic-pmi-gdt-datum-texture",
        "exact_semantic_pmi_datum_texture_and_association_graph_verdict",
        "bit",
        ToleranceSemantics::Exact,
    ),
    (
        "i07-material-process-lot-requirement-evidence",
        "exact_material_process_lot_requirement_evidence_lineage_verdict",
        "bit",
        ToleranceSemantics::Exact,
    ),
    (
        "i07-harness-ewis-connectivity-semantics",
        "exact_harness_connectivity_containment_route_and_identity_graph_verdict",
        "bit",
        ToleranceSemantics::Exact,
    ),
    (
        "i07-exhaustive-semantic-loss-receipt",
        "exact_source_target_atom_reconciliation_and_loss_classification_verdict",
        "bit",
        ToleranceSemantics::Exact,
    ),
    (
        "i07-canonical-export-determinism-and-containment",
        "exact_canonical_bytes_reimport_and_no_partial_publication_verdict",
        "bit",
        ToleranceSemantics::Exact,
    ),
    (
        "i07-safe-opaque-extension-capsules",
        "exact_opaque_octet_token_reference_closure_and_safety_verdict",
        "bit",
        ToleranceSemantics::Exact,
    ),
    (
        "i07-kinematic-pair-semantic-roundtrip",
        "exact_kinematic_relation_frame_freedom_limit_and_coupling_equivalence_verdict",
        "bit",
        ToleranceSemantics::Exact,
    ),
    (
        "i07-partial-bidirectional-edition-migration",
        "exact_admitted_lens_laws_lineage_and_migration_loss_verdict",
        "bit",
        ToleranceSemantics::Exact,
    ),
    (
        "i07-governed-real-ap242-corpus-validation",
        "exact_governed_corpus_coverage_semantic_receipt_and_scoped_validation_verdict",
        "bit",
        ToleranceSemantics::Exact,
    ),
    (
        "i07-proof-carrying-bidirectional-semantic-equivalence",
        "independent_bidirectional_semantic_equivalence_theorem_checker_verdict",
        "bit",
        ToleranceSemantics::Exact,
    ),
    (
        "i07-semantic-sheaf-descent-obstruction-theorem",
        "independent_semantic_descent_gluing_and_obstruction_checker_verdict",
        "bit",
        ToleranceSemantics::Exact,
    ),
    (
        "i07-semantic-equivalence-counterexample-search",
        "exact_nonvacuity_completeness_and_zero_verified_in_domain_counterexample_verdict",
        "bit",
        ToleranceSemantics::Exact,
    ),
];

const LEAVES: &[&str] = &[
    "i07-falsifier-max",
    "i07-governed-real-corpus-validation-max",
    "i07-harness-loss-core",
    "i07-kinematics-max",
    "i07-migration-max",
    "i07-opaque-max",
    "i07-parser-writer-core",
    "i07-pmi-material-core",
    "i07-structure-geometry-core",
    "i07-theorem-descent-max",
];

const LEAF_DECK_MAP: &[(&str, &[&str])] = &[
    (
        "i07-parser-writer-core",
        &[
            CAMPAIGN_POLICY_FIXTURE,
            PROFILE_REGISTRY_FIXTURE,
            "i07-part21-synthetic-corpus",
            "i07-part21-core-holdout",
            "i07-adversarial-reference-graphs",
            "i07-semantic-loss-unsupported-corpus",
        ],
    ),
    (
        "i07-structure-geometry-core",
        &[
            CAMPAIGN_POLICY_FIXTURE,
            PROFILE_REGISTRY_FIXTURE,
            "i07-product-assembly-configurations",
            "i07-product-assembly-core-holdout",
            "i07-mixed-exact-tessellated-geometry",
            "i07-mixed-geometry-core-holdout",
            "i07-adversarial-reference-graphs",
        ],
    ),
    (
        "i07-pmi-material-core",
        &[
            CAMPAIGN_POLICY_FIXTURE,
            PROFILE_REGISTRY_FIXTURE,
            "i07-semantic-pmi-gdt-texture",
            "i07-semantic-pmi-core-holdout",
            "i07-material-process-lot-references",
            "i07-material-process-core-holdout",
            "i07-semantic-loss-unsupported-corpus",
        ],
    ),
    (
        "i07-harness-loss-core",
        &[
            CAMPAIGN_POLICY_FIXTURE,
            PROFILE_REGISTRY_FIXTURE,
            "i07-harness-ewis-graphs",
            "i07-harness-ewis-core-holdout",
            "i07-semantic-loss-unsupported-corpus",
            "i07-semantic-loss-core-holdout",
            "i07-adversarial-reference-graphs",
        ],
    ),
    (
        "i07-opaque-max",
        &[
            CAMPAIGN_POLICY_FIXTURE,
            PROFILE_REGISTRY_FIXTURE,
            "i07-safe-opaque-extensions",
            "i07-safe-opaque-max-holdout",
            "i07-adversarial-reference-graphs",
            "i07-semantic-loss-unsupported-corpus",
        ],
    ),
    (
        "i07-kinematics-max",
        &[
            CAMPAIGN_POLICY_FIXTURE,
            PROFILE_REGISTRY_FIXTURE,
            "i07-product-assembly-configurations",
            "i07-kinematic-pair-assemblies",
            "i07-kinematic-pair-max-holdout",
            "i07-semantic-loss-unsupported-corpus",
        ],
    ),
    (
        "i07-migration-max",
        &[
            CAMPAIGN_POLICY_FIXTURE,
            PROFILE_REGISTRY_FIXTURE,
            "i07-edition-profile-migration-pairs",
            "i07-edition-migration-max-holdout",
            "i07-semantic-loss-unsupported-corpus",
        ],
    ),
    (
        "i07-governed-real-corpus-validation-max",
        &[
            CAMPAIGN_POLICY_FIXTURE,
            PROFILE_REGISTRY_FIXTURE,
            REAL_CORPUS_SLOT,
        ],
    ),
    (
        "i07-theorem-descent-max",
        &[
            CAMPAIGN_POLICY_FIXTURE,
            PROFILE_REGISTRY_FIXTURE,
            "i07-semantic-equivalence-theorem-card",
            "i07-semantic-descent-theorem-card",
            "i07-product-assembly-configurations",
            "i07-semantic-loss-unsupported-corpus",
        ],
    ),
    (
        "i07-falsifier-max",
        &[
            CAMPAIGN_POLICY_FIXTURE,
            PROFILE_REGISTRY_FIXTURE,
            "i07-semantic-adversary-microgrammar",
            "i07-semantic-adversary-max-holdout",
            "i07-semantic-equivalence-theorem-card",
            "i07-semantic-descent-theorem-card",
            "i07-edition-profile-migration-pairs",
        ],
    ),
];

const HELD_OUT_FIXTURES: &[&str] = &[
    "i07-edition-migration-max-holdout",
    "i07-harness-ewis-core-holdout",
    "i07-kinematic-pair-max-holdout",
    "i07-material-process-core-holdout",
    "i07-mixed-geometry-core-holdout",
    "i07-part21-core-holdout",
    "i07-product-assembly-core-holdout",
    "i07-safe-opaque-max-holdout",
    "i07-semantic-adversary-max-holdout",
    "i07-semantic-loss-core-holdout",
    "i07-semantic-pmi-core-holdout",
];

const HELD_OUT_AUTHORITY: &[(&str, CampaignTier, &str)] = &[
    (
        "i07-edition-migration-max-holdout",
        CampaignTier::Max,
        "i07-migration-max",
    ),
    (
        "i07-harness-ewis-core-holdout",
        CampaignTier::Core,
        "i07-harness-loss-core",
    ),
    (
        "i07-kinematic-pair-max-holdout",
        CampaignTier::Max,
        "i07-kinematics-max",
    ),
    (
        "i07-material-process-core-holdout",
        CampaignTier::Core,
        "i07-pmi-material-core",
    ),
    (
        "i07-mixed-geometry-core-holdout",
        CampaignTier::Core,
        "i07-structure-geometry-core",
    ),
    (
        "i07-part21-core-holdout",
        CampaignTier::Core,
        "i07-parser-writer-core",
    ),
    (
        "i07-product-assembly-core-holdout",
        CampaignTier::Core,
        "i07-structure-geometry-core",
    ),
    (
        "i07-safe-opaque-max-holdout",
        CampaignTier::Max,
        "i07-opaque-max",
    ),
    (
        "i07-semantic-adversary-max-holdout",
        CampaignTier::Max,
        "i07-falsifier-max",
    ),
    (
        "i07-semantic-loss-core-holdout",
        CampaignTier::Core,
        "i07-harness-loss-core",
    ),
    (
        "i07-semantic-pmi-core-holdout",
        CampaignTier::Core,
        "i07-pmi-material-core",
    ),
];

const SEED_FIXTURES: &[(&str, &str, &str)] = &[
    ("i07-part21-synthetic-corpus", "i07/part21", "0..=4095"),
    ("i07-part21-core-holdout", "i07/part21", "65536..=69631"),
    (
        "i07-product-assembly-configurations",
        "i07/assembly",
        "0..=4095",
    ),
    (
        "i07-product-assembly-core-holdout",
        "i07/assembly",
        "65536..=69631",
    ),
    (
        "i07-mixed-exact-tessellated-geometry",
        "i07/geometry",
        "0..=4095",
    ),
    (
        "i07-mixed-geometry-core-holdout",
        "i07/geometry",
        "65536..=69631",
    ),
    ("i07-semantic-pmi-gdt-texture", "i07/pmi", "0..=4095"),
    ("i07-semantic-pmi-core-holdout", "i07/pmi", "65536..=69631"),
    (
        "i07-material-process-lot-references",
        "i07/material",
        "0..=4095",
    ),
    (
        "i07-material-process-core-holdout",
        "i07/material",
        "65536..=69631",
    ),
    ("i07-harness-ewis-graphs", "i07/harness", "0..=4095"),
    (
        "i07-harness-ewis-core-holdout",
        "i07/harness",
        "65536..=69631",
    ),
    (
        "i07-semantic-loss-unsupported-corpus",
        "i07/loss",
        "0..=4095",
    ),
    (
        "i07-semantic-loss-core-holdout",
        "i07/loss",
        "65536..=69631",
    ),
    (
        "i07-adversarial-reference-graphs",
        "i07/reference-graph",
        "0..=4095",
    ),
    ("i07-safe-opaque-extensions", "i07/opaque", "0..=4095"),
    (
        "i07-safe-opaque-max-holdout",
        "i07/opaque",
        "131072..=135167",
    ),
    (
        "i07-kinematic-pair-assemblies",
        "i07/kinematics",
        "0..=4095",
    ),
    (
        "i07-kinematic-pair-max-holdout",
        "i07/kinematics",
        "131072..=135167",
    ),
    (
        "i07-edition-profile-migration-pairs",
        "i07/migration",
        "0..=4095",
    ),
    (
        "i07-edition-migration-max-holdout",
        "i07/migration",
        "131072..=135167",
    ),
    (
        "i07-semantic-adversary-microgrammar",
        "i07/adversary",
        "0..=4095",
    ),
    (
        "i07-semantic-adversary-max-holdout",
        "i07/adversary",
        "131072..=135167",
    ),
];

const UNIT_CASES: &[&str] = &[
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

const COMMON_OBS_EVENTS: &[&str] = &[
    "artifact.published",
    "campaign.admitted",
    "campaign.requested",
    "case.cancel.requested",
    "case.cancel.observed",
    "case.cancelled",
    "case.drain.started",
    "case.drain.completed",
    "case.failed",
    "case.finalized",
    "case.started",
    "checkpoint.committed",
    "claim.adjudicated",
    "failure_bundle.retained",
    "evidence.atomic_publish_committed",
    "integrity.failed",
    "loss.itemized",
    "resume.completed",
];

fn authored_spec(fixture: &FixturePin) -> &'static str {
    match fixture.source {
        FixtureSource::AuthoredSpec { spec } => spec,
        FixtureSource::External { .. } => {
            panic!(
                "fixture '{}' must be an authored I07 specification",
                fixture.id
            )
        }
    }
}

fn fixture<'a>(fixtures: &'a [FixturePin], id: &str) -> &'a FixturePin {
    fixtures
        .iter()
        .find(|candidate| candidate.id == id)
        .unwrap_or_else(|| panic!("missing I07 fixture '{id}'"))
}

fn all_claim_ids() -> BTreeSet<&'static str> {
    CLAIM_AUTHORITY.iter().map(|(id, _, _, _)| *id).collect()
}

fn all_authority_ids() -> BTreeSet<&'static str> {
    all_claim_ids()
        .into_iter()
        .chain(LEAVES.iter().copied())
        .collect()
}

fn reversed_static(items: &'static [&'static str]) -> &'static [&'static str] {
    let mut reversed = items.to_vec();
    reversed.reverse();
    Box::leak(reversed.into_boxed_slice())
}

fn exact_section<'a>(text: &'a str, start: &str, end: &str) -> &'a str {
    let (_, tail) = text
        .split_once(start)
        .unwrap_or_else(|| panic!("missing exact section start '{start}'"));
    let (section, _) = tail
        .split_once(end)
        .unwrap_or_else(|| panic!("missing exact section end '{end}'"));
    section
}

fn normalized_table_heads(section: &str) -> Vec<String> {
    section
        .split(';')
        .map(str::trim)
        .filter(|entry| !entry.is_empty())
        .map(|entry| {
            entry
                .split_ascii_whitespace()
                .take(3)
                .collect::<Vec<_>>()
                .join(" ")
        })
        .collect()
}

fn normalized_table_entries(section: &str) -> Vec<String> {
    section
        .split(';')
        .map(str::trim)
        .filter(|entry| !entry.is_empty())
        .map(|entry| entry.split_ascii_whitespace().collect::<Vec<_>>().join(" "))
        .collect()
}

#[test]
#[allow(clippy::too_many_lines)]
fn i07_freezes_with_exact_lattice_and_once_only_mapping() {
    let frozen = i07_draft().freeze().expect("the I07 seed must freeze");
    assert_eq!(frozen.initiative(), "I07");
    assert_eq!(frozen.version(), 1);
    assert!(
        i07_draft()
            .title
            .contains("governed real-corpus validation")
    );
    assert!(
        i07_draft()
            .title
            .contains("proof-carrying sheaf/stack descent")
    );
    assert_eq!(frozen.claims().len(), 15);
    assert_eq!(frozen.fixtures().len(), 27);
    assert_eq!(frozen.obligations().len(), 10);
    assert_eq!(frozen.waivers().len(), 1);

    let mut lattice = BTreeMap::new();
    let mut polarity = BTreeMap::new();
    for claim in frozen.claims() {
        *lattice
            .entry(format!("{:?}", claim.ambition))
            .or_insert(0_usize) += 1;
        *polarity
            .entry(format!("{:?}", claim.polarity))
            .or_insert(0_usize) += 1;
    }
    assert_eq!(lattice.get("Solid"), Some(&8));
    assert_eq!(lattice.get("Frontier"), Some(&4));
    assert_eq!(lattice.get("Moonshot"), Some(&3));
    assert_eq!(polarity.get("Affirmative"), Some(&14));
    assert_eq!(polarity.get("Refutation"), Some(&1));

    for (id, ambition, expected_polarity, evidence_tier) in CLAIM_AUTHORITY {
        let claim = frozen
            .claim(id)
            .unwrap_or_else(|| panic!("missing claim '{id}'"));
        assert_eq!(claim.ambition, *ambition, "ambition drift for {id}");
        assert_eq!(
            claim.polarity, *expected_polarity,
            "polarity drift for {id}"
        );
        assert_eq!(claim.evidence_tier, *evidence_tier, "tier drift for {id}");
        assert!(
            claim.oracle.independent,
            "oracle is not independent for {id}"
        );
        assert!(!claim.oracle.identity.trim().is_empty());
        assert!(!claim.oracle.tcb_overlap.trim().is_empty());
        assert!(!claim.fallback.trim().is_empty());
        assert!(!claim.no_claim.trim().is_empty());
    }
    assert_eq!(
        frozen
            .claims()
            .iter()
            .map(|claim| claim.id)
            .collect::<BTreeSet<_>>(),
        all_claim_ids()
    );
    assert_eq!(CLAIM_ORACLES.len(), CLAIM_AUTHORITY.len());
    assert_eq!(
        CLAIM_ORACLES
            .iter()
            .map(|(id, _)| *id)
            .collect::<BTreeSet<_>>(),
        all_claim_ids(),
        "every I07 claim must have one exact independent-oracle route"
    );
    assert_eq!(
        CLAIM_ORACLES
            .iter()
            .map(|(_, identity)| *identity)
            .collect::<BTreeSet<_>>()
            .len(),
        CLAIM_ORACLES.len(),
        "I07 claims must not alias the same oracle identity"
    );
    for (id, identity) in CLAIM_ORACLES {
        assert_eq!(
            frozen.claim(id).expect("oracle claim").oracle.identity,
            *identity,
            "independent-oracle route drift for {id}"
        );
    }
    assert_eq!(CLAIM_ACCEPTANCE.len(), CLAIM_AUTHORITY.len());
    assert_eq!(
        CLAIM_ACCEPTANCE
            .iter()
            .map(|(id, _, _, _)| *id)
            .collect::<BTreeSet<_>>(),
        all_claim_ids(),
        "every I07 claim must have an exact QoI/unit/tolerance regression lock"
    );
    for (id, qoi, unit, tolerance) in CLAIM_ACCEPTANCE {
        let claim = frozen
            .claim(id)
            .unwrap_or_else(|| panic!("missing acceptance claim '{id}'"));
        assert_eq!(claim.qoi, *qoi, "QoI drift for {id}");
        assert_eq!(claim.unit, *unit, "unit drift for {id}");
        assert_eq!(claim.tolerance, *tolerance, "tolerance drift for {id}");
    }

    let mut coverage = BTreeMap::<&str, usize>::new();
    let mut leaves = BTreeSet::new();
    let mut tiers = BTreeMap::new();
    assert_eq!(LEAF_DECK_MAP.len(), LEAVES.len());
    assert_eq!(
        LEAF_DECK_MAP
            .iter()
            .map(|(leaf, _)| *leaf)
            .collect::<BTreeSet<_>>(),
        LEAVES.iter().copied().collect::<BTreeSet<_>>(),
        "the exact deck map must cover every I07 leaf once"
    );
    for row in frozen.obligations() {
        assert!(leaves.insert(row.leaf()), "duplicate leaf '{}'", row.leaf());
        *tiers.entry(format!("{:?}", row.tier())).or_insert(0_usize) += 1;
        let expected_decks = LEAF_DECK_MAP
            .iter()
            .find_map(|(leaf, decks)| (*leaf == row.leaf()).then_some(*decks))
            .unwrap_or_else(|| panic!("missing exact deck map for '{}'", row.leaf()));
        assert_eq!(
            row.decks().len(),
            expected_decks.len(),
            "deck multiplicity drift for {}",
            row.leaf()
        );
        assert_eq!(
            row.decks().iter().copied().collect::<BTreeSet<_>>(),
            expected_decks.iter().copied().collect::<BTreeSet<_>>(),
            "exact deck-set drift for {}",
            row.leaf()
        );
        assert_eq!(
            row.unit_cases().iter().copied().collect::<BTreeSet<_>>(),
            UNIT_CASES.iter().copied().collect::<BTreeSet<_>>(),
            "unit-case taxonomy drift for {}",
            row.leaf()
        );
        for claim in row.claims_covered() {
            *coverage.entry(claim).or_insert(0) += 1;
            let ambition = frozen.claim(claim).expect("covered claim").ambition;
            assert_eq!(
                row.tier(),
                if ambition == Ambition::Solid {
                    CampaignTier::Core
                } else {
                    CampaignTier::Max
                },
                "I07 lattice stage drift for {claim}"
            );
        }
    }
    assert_eq!(leaves, LEAVES.iter().copied().collect::<BTreeSet<_>>());
    assert_eq!(tiers.get("Core"), Some(&4));
    assert_eq!(tiers.get("Max"), Some(&6));
    assert_eq!(
        coverage.keys().copied().collect::<BTreeSet<_>>(),
        all_claim_ids()
    );
    assert!(
        coverage.values().all(|count| *count == 1),
        "every claim must be owned by exactly one execution leaf: {coverage:?}"
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn i07_profile_policy_and_authority_boundaries_are_frozen() {
    let frozen = i07_draft().freeze().expect("freeze");
    let policy = authored_spec(fixture(frozen.fixtures(), CAMPAIGN_POLICY_FIXTURE));
    let profile = authored_spec(fixture(frozen.fixtures(), PROFILE_REGISTRY_FIXTURE));
    let waiver = frozen.waivers()[0];

    for required in [
        "freeze is preregistration only",
        "ORTHOGONAL_STATES",
        "ExecutionState={Succeeded,Failed,Cancelled,TimedOut,BudgetExhausted,InfrastructureFailed,IntegrityFailed}",
        "ClaimState={Pending,Passed,Failed,Refuted,Unknown,Unsupported}",
        "EvidenceCompleteness={CompleteEvidence,PartialEvidence,NoEvidence}",
        "OpaquePreservedNoSemanticAuthority",
        "StandardAxis={NotAssessed,SyntheticProfileOnly,ScopedClauseConformance}",
        "InteroperabilityAxis={NotAssessed,ScopedRequiredMatrixValidated}",
        "IndustrialCorpusAxis={NotAssessed,PublicCorpusOnly,GovernedBlindScopedValidation}",
        "RegulatoryAxis={NotGranted}",
        "AXIS_AGGREGATION",
        "CROSS_AXIS_VALIDITY",
        "case.finalized forbids ClaimState=Pending",
        "ClaimState=Passed and every favorable non-NotAssessed validation axis require ExecutionState=Succeeded",
        "EvidenceCompleteness=CompleteEvidence",
        "ClaimState=Refuted likewise requires CompleteEvidence",
        "ClaimState=Unsupported requires CompleteEvidence of exact non-applicability",
        "PartialEvidence or NoEvidence forces ClaimState=Unknown",
        "ExecutionState=Failed, Cancelled, TimedOut, BudgetExhausted, InfrastructureFailed or IntegrityFailed cannot carry Passed, Refuted",
        "PartialEvidence and NoEvidence have zero promotion authority",
        "Public/previously accessed corpora remain PublicCorpusOnly",
        "PROFILE_BINDING",
        "PROFILE_HASH_ALGEBRA",
        "H(context,payload)=BLAKE3::derive_key",
        "FRAME_BYTES(x)=U64LE(len(x)) || x",
        "FRAME_LIST(xs)=U64LE(count(xs))",
        "TAGGED_FIELD(tag,payload)=U32LE(tag)||FRAME_BYTES(payload)",
        "tag is the one-based field ordinal",
        "exact_case_sensitive_UTF8_without_normalization",
        "raw 32-byte output",
        "Digest payloads are exactly 32 nonzero raw bytes",
        "PROFILE_BOOTSTRAP_CAPS",
        "MAX_PROFILE_UTF8_BYTES=65536",
        "MAX_PROFILE_NESTED_COLLECTION_ITEMS=65536",
        "MAX_PROFILE_NESTED_PAYLOAD_BYTES=16777216",
        "MAX_PROFILE_DEFINITION_SCHEMA_BYTES=4194304",
        "MAX_CANONICAL_PROFILE_DEFINITION_BYTES=16777216",
        "MAX_PROFILE_DEFINITION_FIELDS=4096",
        "MAX_PROFILE_REGISTRY_ENTRIES=4096",
        "MAX_PROFILE_REGISTRY_PROJECTION_BYTES=67108864",
        "MAX_REQUIRED_ASSIGNMENT_KEYS=1048576",
        "MAX_REQUIRED_ASSIGNMENT_UNIVERSE_BYTES=268435456",
        "MAX_PROFILE_ASSIGNMENT_BYTES=536870912",
        "checked arithmetic before allocation, recursion or hash finalization",
        "enforced before the decoder can trust or consult limit_policy_digest",
        "MAX_PROFILE_NESTED_PAYLOAD_BYTES bounds byte_len(x) inside FRAME_BYTES(x), excluding that frame's outer eight-byte U64LE length",
        "complete frame size is checked_add(8,byte_len(x))",
        "Top-level ProfileKeyBytes, ProfileDefinitionSchemaBytes",
        "CanonicalProfileDefinitionBytes, ProfileRegistry, RequiredAssignmentUniverse and ProfileAssignment",
        "projections each have depth 0",
        "depth 64 is admitted",
        "attempt to enter depth 65 is refused before recursion or allocation",
        "stream into their domain-separated hashes instead of requiring whole-projection allocation",
        "exact-payload-cap, payload-cap-plus-one, exact-depth-64, depth-65, U64::MAX length/count",
        "checked-add/multiply overflow, truncated nested payload, duplicate/order and trailing-byte mutation tests",
        "CanonicalProfileDefinitionBytes is FRAME_UTF8('i07-profile-definition-v1')",
        "ProfileDefinitionSchemaBytes is the exact content-addressed canonical schema artifact",
        "FRAME_UTF8('i07-profile-definition-schema-v1')",
        "ProfileDefinitionSchemaRoot=H('org.frankensim.i07.profile-definition-schema.v1'",
        "raw(ProfileDefinitionSchemaRoot)||U32LE(field_count)",
        "field_count must equal the exact schema table cardinality",
        "definition-AST fields, each as one TAGGED_FIELD in exact strictly increasing schema ordinal order",
        "maps are lexicographically raw-key sorted and unique",
        "schema bytes/root and complete definition-field table are mandatory content-addressed inputs",
        "missing, mismatched or unavailable schema root is IntegrityFailed",
        "raw(ProfileKeyDigest)||raw(ProfileDefinitionSchemaRoot)",
        "ProfileRegistryRoot",
        "FRAME_UTF8('i07-profile-registry-v1')",
        "ProfileAssignmentRoot",
        "ProfileAssignmentRoot=H('org.frankensim.i07.profile-assignment.v1',raw(ProfileRegistryRoot)",
        "AssignmentTargetKind={Fixture=1,Case=2,CorpusItem=3,Stratum=4}",
        "ProfileRole={Input=1,Output=2,SourceEndpoint=3,TargetEndpoint=4,ReferenceOracle=5,MigrationIntermediate=6,TheoremInstantiation=7}",
        "AssignmentLogicalKey",
        "RequiredAssignmentUniverseRoot=H('org.frankensim.i07.required-assignment-universe.v1'",
        "decoded AssignmentRow logical-key set must equal, not merely contain",
        "A duplicate ProfileKeyDigest is IntegrityFailed",
        "Duplicate, missing or unexpected logical keys",
        "Ordinals have no profile-authority semantics",
        "FILE_SCHEMA must match the already bound profile",
        "--profile <profile-root> means exactly the assigned ProfileEntryRoot",
        "Every required assignment row runs",
        "profile-definition-schema root, registry root, required-assignment-universe root, assignment root",
        "FrozenManifest::amend successor",
        "REQUEST_SET",
        "never a favorable caller-selected subset",
        "EXIT_CODES: 0=all_requested_claims_passed",
        "30=integrity_failed; 31=infrastructure_failed",
        "fail-closed code in order 30,31,21,22,11,10,20,13,12,0",
        "request_cancel? -> drain -> finalize -> adjudicate -> atomic_publish",
        "successful publication requires one complete atomic AP242/loss-receipt pair",
        "non-success publishes no AP242 artifact and records the explicit absent-artifact disposition",
        "No path publishes one member of a success pair, a partial member, or an unrecorded absence",
        "LIFECYCLE_EVENT_CONTRACT",
        "every submitted attempt emits exactly one campaign.requested",
        "A refused admission emits no campaign.admitted, case.started or work event",
        "Every admitted attempt emits exactly one campaign.admitted and then exactly one case.started",
        "first accepted cancellation for one open scope/epoch is the primary request",
        "exactly one case.cancel.requested binds request id, scope root, request-event identity",
        "primary_request_logical_sequence",
        "supervisor-clock request interval",
        "calibration-artifact root and conservative timestamp-uncertainty enclosure",
        "atomically seals the request-time observer-set root/count, descendant-frontier root/count, registration epoch and exact terminal reservation",
        "worst-case terminal row-and-encoded-byte reservation",
        "a racing registration either commits before the seal",
        "included with deterministic catalog membership and reservation receipts, or is refused",
        "Exact-key byte-identical replay returns the primary receipt without another event",
        "AlreadySealed, AlreadyDraining or AlreadyFinalized response",
        "Every live cancellable logical work unit in the sealed observer set emits exactly one case.cancel.observed",
        "binding the primary request id/event identity/primary_request_logical_sequence",
        "sealed observer root/count, descendant-frontier root/count, registration epoch, admitted unit/tile id",
        "observation_event_logical_sequence",
        "optional last-completed boundary",
        "supervisor-domain calibrated observation interval with an outward uncertainty enclosure",
        "On a cancellation path, exactly one case.drain.started follows the primary seal",
        "On a non-cancellation path, exactly one case.drain.started follows the earliest applicable work-completion, domain-failure, budget, campaign-timeout or infrastructure drain cause",
        "Observation forbids every new domain, candidate, protected-access, adjudication or publication effect in that unit",
        "only already reserved bounded cleanup, descendant join, exit-receipt, drain and finalization work may continue",
        "On a cancellation path, case.drain.completed carries the sealed roots/counts and a complete exit-reconciliation root",
        "every admitted descendant has one matching membership-plus-exit receipt",
        "On a non-cancellation path, it instead carries the exact admitted-frontier root/count and complete exit reconciliation",
        "it neither requires nor fabricates a cancellation observation",
        "empty sealed observer set requires a present",
        "independently reconciled root with cardinality zero and exactly zero observations",
        "request-to-observation maximum is exactly zero at seal",
        "DrainTrigger={CancellationObserved,EmptyObserverSeal,ObservationTimeoutDrain,InfrastructureFailure,NonCancellationDrain}",
        "primary-seal interval for an independently verified empty observer set",
        "causal rank InfrastructureFailure=0,CancellationObserved=1,EmptyObserverSeal=2,ObservationTimeoutDrain=3,NonCancellationDrain=4",
        "inclusive observation deadline is checked_add(request_interval_lower,250000000 ns) for Core",
        "checked_add(request_interval_lower,1000000000 ns) for Max",
        "Conservative trigger-to-drained latency is checked_sub(drained_interval_upper,trigger_interval_lower)",
        "conservative drained-to-finalized latency is checked_sub(finalized_interval_upper,drained_interval_lower)",
        "Core caps are inclusively 2000000000 ns and 2000000000 ns respectively",
        "Max caps are inclusively 8000000000 ns and 8000000000 ns",
        "first nanosecond outside its inclusive bound",
        "observation whose interval upper equals the inclusive observation deadline",
        "drain/finalization conservative latency equal to its inclusive cap, is on time",
        "At the derived timeout onset checked_add(deadline,1 ns), an observation completion loses to and latches TimedOut",
        "conservative latency of checked_add(cap,1 ns) latches TimedOut",
        "TERMINAL_ARBITRATION",
        "IntegrityFailed regardless of an operational candidate",
        "earliest terminal-eligible logical boundary wins",
        "InfrastructureFailed > TimedOut > Failed > Cancelled > BudgetExhausted > Succeeded",
        "Cancelled is eligible only when every required observation, descendant exit, drain and finalization receipt is complete",
        "domain/acceptance effect whose own durable commit point precedes the primary cancellation",
        "no post-cancellation effect may do so",
        "missed observation, trigger-to-drained, drained-to-finalized or campaign deadline supplies TimedOut",
        "ExecutionState=IntegrityFailed and EvidenceCompleteness=NoEvidence",
        "case.cancelled is the exactly-once terminal-selection event",
        "case.failed is exactly once for selected Failed, TimedOut, BudgetExhausted, InfrastructureFailed or IntegrityFailed",
        "Every admitted attempt emits exactly one case.finalized",
        "final execution/input/integrity/domain/support state and immutable evidence payload",
        "exactly one claim.adjudicated per requested claim binds its claim/evidence axes",
        "Successful adjudication requires exactly one evidence.atomic_publish_committed followed by exactly one artifact.published",
        "may publish at most one paired atomic terminal/failure evidence transaction",
        "with zero promotion authority",
        "a partial/unpaired publication is IntegrityFailed",
        "Duplicate, missing, out-of-order or inapplicable lifecycle events",
        "work/cancellation/drain after case.finalized",
        "conservative latency_i is checked_sub(observation_interval_upper_i,request_interval_lower)",
        "in one supervisor monotonic domain",
        "receipt-bound calibration and outward drift enclosure",
        "inclusively capped at 250000000 ns for Core and 1000000000 ns for Max",
        "unknown/duplicate/extra/root-, count- or epoch-mismatched observations is IntegrityFailed",
        "missing or late observation crossing its inclusive deadline is TimedOut",
        "preflight inability to establish a common clock or bounded observer is an Unsupported admission refusal",
        "receipt-verified observer/supervisor loss is InfrastructureFailed",
        "None can succeed",
        "Raw calibrated timestamps, uncertainty and watchdog arrival remain telemetry",
        "request event identities, logical sequences, unit/boundary ids, sealed root/epoch and selected disposition remain canonical",
        "exactly one DispositionRecord per source atom",
        "exactly one OriginRecord per target atom",
        "forward/inverse incidence relations must be exact transposes",
        "Split, merge, synthesis and deletion",
        "legitimately empty admitted source and target inventory",
        "production import/export path performs no network or URI resolution",
        "Development-only theorem lanes",
        "HOLDOUT_REALIZATION",
        "JOINT_REVEAL",
        "Sequential reveal",
        "GOVERNED_EXTERNAL",
        "I07GovernedCorpusDischargeReceiptV1",
        "exact waiver subject/predicate",
        "domain-separated transaction-intent digest",
        "One atomic FrozenManifest::amend authority transaction",
        "I07GovernedCorpusDischargeEnvelopeV1",
        "same deck id",
        "verifies the AmendmentRecord",
        "advances the authority head",
        "TRANSACTION_INTENT_V1 defines TransactionIntentDigest=BLAKE3::derive_key(",
        "'org.frankensim.i07.governed-transaction-intent.v1',P)",
        "canonical successor intent for exactly one fenced authority transaction",
        "predecessor manifest and governance-stage receipt",
        "expected authority-head digest/generation",
        "immutable CandidateFrozen candidate/checker/toolchain/decision/ProfileRegistryRoot/ProfileAssignmentRoot commitment",
        "lexically ordered retired waiver subjects",
        "role-addressed (slot_id,artifact_role,artifact_schema_digest,protected_root) records",
        "role-addressed future envelope slots that bind each output artifact_schema_digest before Pending",
        "Separately sorted slot-id and root lists are forbidden",
        "coupled transaction group",
        "governance stage",
        "governance stage/scope",
        "I07_ENVELOPE_DIGEST_PENDING_V1",
        "I07_AMENDMENT_RECORD_PENDING_V1",
        "encoded union is Pending and cannot alias a digest",
        "initiative_id='I07'",
        "schema_identity='i07-governed-transaction-intent-v1'",
        "checked successor_version=predecessor_version+1",
        "governance_stage=RealizationCommitted",
        "expected-head compare-and-swap",
        "prevents ABA replay",
        "retired-waiver list must equal the exact predecessor-minus-successor waiver-set difference",
        "form an exact role-addressed bijection",
        "Each future output artifact schema is a mandatory content-addressed canonical schema",
        "schema whose bytes are available and independently rehashed before realization",
        "realized envelope must decode under the byte-identical precommitted schema digest",
        "caller-selected or merely nonzero schema digest has no authority",
        "GOVERNED_OUTPUT_SCHEMA_AUTHORITY_V1",
        "protocol authority is closed to exactly one output and has no extension row",
        "GovernedOutputSchemaMatrixBytesV1 is byte-for-byte",
        "FRAME_UTF8('i07-governed-output-schema-matrix-v1')",
        "protocol_kind=GovernedRealAp242CorpusDischarge=1",
        "governance_stage=RealizationCommitted=3",
        "authority_scope='i07.governed-real-ap242-corpus-discharge'",
        "target_slot_id='i07-governed-real-ap242-corpus'",
        "role_prefix='GovernedRealAp242CorpusDischarge/'",
        "min_count=max_count=1",
        "related_waiver_rule=ExactRetiredSubject=1",
        "protected_binding_rule=Required=1",
        "exact matrix length is (8+36)+3*2+(8+8+40)+(8+8+30)+(8+8+33)+2*8+2*2=221 bytes",
        "sole retired waiver subject, protected slot_id, target_slot_id and FutureArtifact related_waiver_subject",
        "artifact_role='GovernedRealAp242CorpusDischarge/'||related_waiver_subject",
        "'GovernedRealAp242CorpusDischarge/i07-governed-real-ap242-corpus'",
        "GovernedOutputSchemaMatrixDigest=BLAKE3::derive_key(",
        "'org.frankensim.i07.governed-output-schema-matrix.v1',GovernedOutputSchemaMatrixBytesV1)",
        "I07GovernedCorpusDischargeEnvelopeSchemaBytesV1 starts with",
        "FRAME_UTF8('i07-governed-corpus-discharge-envelope-schema-v1')||U16LE(34)",
        "SchemaField=U16LE(ordinal)||FRAME_BYTES(Utf8(field_name))||U16LE(field_type)||U16LE(1)||FRAME_BYTES(ConstraintAst)",
        "field_type is Utf8=1,U16LE=2,U64LE=3,Digest=4",
        "EqIntentTag=U16LE(4)||U16LE(tag)",
        "EqCandidateFrozenComponent=U16LE(6)||U16LE(component)",
        "EqTransactionIntentDigest=U16LE(7)",
        "EqPrecommitAuthorizationDigest=U16LE(8)",
        "EqAuthorizationBinding=U16LE(9)||FRAME_BYTES(Utf8(binding_name))",
        "EqAuthorizationBoundU64Range=U16LE(10)||FRAME_BYTES(Utf8(binding_name))||U64LE(min)||U64LE(max)",
        "authorization-binding namespace is closed to exactly the seventeen case-sensitive binding_name values",
        "fields 14..23 and 28..34",
        "missing, extra, duplicate, aliased or differently typed binding is IntegrityFailed",
        "I07GovernedCorpusDischargeEnvelopeBytesV1 is exactly",
        "FRAME_UTF8('i07-governed-corpus-discharge-envelope-v1')||U16LE(34)",
        "concat(U16LE(ordinal)||U64LE(payload_byte_len)||payload in increasing ordinal order)",
        "Every one of the 34 schema constraints is checked",
        "exact envelope length is",
        "(8+41)+2+34*(2+8)+(6*8+41+3+40+30+30+63)+2*2+8+25*32=1458 bytes",
        "I07GovernedCorpusDischargeEnvelopeSchemaDigest=BLAKE3::derive_key(",
        "'org.frankensim.i07.governed-corpus-discharge-envelope-schema.v1'",
        "I07GovernedCorpusDischargeEnvelopeDigest=BLAKE3::derive_key(",
        "'org.frankensim.i07.governed-corpus-discharge-envelope.v1'",
        "realized FutureArtifact digest and successor-installed same-ID External root",
        "GovernedOutputSchemaRecordV1=U16LE(1)||U16LE(3)",
        "GovernedOutputSchemaRecordV1 is exactly 2*2+(8+8+40)+(8+8+63)+32=171 bytes",
        "GovernedOutputSchemaSetRoot=BLAKE3::derive_key(",
        "'org.frankensim.i07.governed-output-schema-set.v1'",
        "GovernedOutputSchemaMembershipProofV1 is exactly",
        "verification recomputes the singleton set root",
        "proof is exactly 32+2*8+8+171=227 bytes",
        "GovernanceCommitted freezes the exact matrix bytes",
        "CandidateFrozen copies the same set root",
        "P binds it at tag 0x0014",
        "precommit authorization carries the exact membership proof for the sole FutureArtifact",
        "all three root equalities",
        "both exact digest constraints, and every authorization binding must byte-equal its already committed source",
        "permissive substitute schema",
        "matrix extension",
        "cross-protocol, cross-stage or cross-scope schema",
        "GovernanceCommitted/CandidateFrozen/P root mismatch",
        "Untrusted matrix bytes are preflight-capped at 4096",
        "exact 221-byte comparison",
        "schema bytes at 4194304",
        "membership-proof bytes at 8192",
        "exact 227-byte comparison",
        "realized-envelope bytes are exactly 1458",
        "schema field count is exactly 34 and schema-set count exactly 1",
        "depth 64 is admitted",
        "attempt to enter depth 65 is refused before recursion or allocation",
        "two independently implemented matrix/schema/set/proof/envelope encoders and decoders",
        "published exact-byte and exact-root KATs",
        "permissive-schema, missing/extra-output",
        "exact-cap/cap-plus-one",
        "exact-envelope-length/length-plus-one",
        "depth-64/depth-65",
        "envelope-slot-id set, final-successor slot id and AmendmentRecord slot id are pairwise disjoint",
        "ReceiptSchemaSetRoot=H('org.frankensim.i07.receipt-schema-set.v1'",
        "PrecommitAuthorizationSchemaBytes and CommitReceiptSchemaBytes are mandatory content-addressed canonical schema artifacts",
        "FRAME_UTF8('i07-governed-corpus-discharge-precommit-authorization-schema-v1')",
        "FRAME_UTF8('i07-governed-corpus-discharge-commit-receipt-schema-v1')",
        "PrecommitAuthorizationSchemaDigest=H(",
        "'org.frankensim.i07.precommit-authorization-schema.v1',PrecommitAuthorizationSchemaBytes)",
        "CommitReceiptSchemaDigest=H('org.frankensim.i07.commit-receipt-schema.v1'",
        "exactly those two case-sensitive role/schema pairs",
        "arbitrary nonzero schema digests and role aliases are forbidden",
        "recomputes both schema digests and the set root from the exact available schema bytes",
        "embeds its role, schema digest and ReceiptSchemaSetRoot membership proof",
        "missing, extra, swapped or mismatched role/schema is IntegrityFailed",
        "may bind the already-existing predecessor-stage receipt digest",
        "never serializes the newly created authorization/commit receipt digest",
        "precommit authorization binds the predecessor/stage receipt, expected old head, ReceiptSchemaSetRoot, GovernedOutputSchemaSetRoot and P",
        "successor embeds the envelope digest",
        "AmendmentRecord binds predecessor and final successor",
        "I07GovernedCorpusDischargeCommitReceiptV1",
        "NewAuthorityHeadDigest=H('org.frankensim.i07.authority-head.v1'",
        "raw(old_head_digest)||U64LE(old_generation)||raw(TransactionIntentDigest)",
        "NewAuthorityHeadDigest||U64LE(checked_add(old_generation,1))",
        "current head byte-equals expected_authority_head",
        "durably writes exactly the derived new head plus the commit receipt in the same atomic commit",
        "CAS failure, generation overflow or a different proposed new head performs no waiver retirement",
        "strict P -> authorization -> envelope -> successor -> AmendmentRecord -> NewAuthorityHead -> commit-receipt DAG",
        "25 exact ASCII bytes I07_TRANSACTION_INTENT_V1 followed by one trailing byte 0x00",
        "then U16LE(20)",
        "Digest is exactly 32 nonzero raw bytes",
        "GovernanceCommitted=1,CandidateFrozen=2,RealizationCommitted=3,RevealedForAdjudication=4,Closed=5",
        "ProtectedBinding=FRAME_BYTES(Utf8(slot_id))||FRAME_BYTES(Utf8(artifact_role))||Digest(artifact_schema)||Digest(protected_root)",
        "ProtectedBindingList=U64LE(count)||concat(FRAME_BYTES(ProtectedBinding))",
        "FutureArtifact=FRAME_BYTES(Utf8(slot_id))||FRAME_BYTES(Utf8(related_waiver_subject))||FRAME_BYTES(Utf8(artifact_role))||Digest(artifact_schema)||digest_union",
        "FutureArtifactList=U64LE(count)||concat(FRAME_BYTES(FutureArtifact))",
        "FutureDigest=FRAME_BYTES(Utf8(slot_id))||digest_union",
        "digest_union is exactly byte 0=Pending with no following bytes or byte 1=Digest followed by 32 raw bytes",
        "tags 0x000e..0x0010 require Pending",
        "P is at most 16777216 bytes",
        "every Utf8 payload is at most 65536 bytes",
        "checked arithmetic against remaining input and these caps before allocation",
        "missing, extra, duplicate, unknown, out-of-order, empty-required, non-minimal, overflowed or trailing bytes",
        "raw or all-zero digest, role or artifact-schema swap",
        "unavailable/mismatched future artifact schema",
        "permissive or unregistered output schema",
        "arbitrary same-ID External fixture",
        "receipt alone",
        "waiver removal alone",
        "split successor",
        "post-reveal amendment",
        "structurally frozen successor",
        "grants no discharge or promotion authority",
        "two independent encoders must match a published 20-field transaction-intent KAT",
        "header/NUL/count/every-tag/type/length/cardinality/order/role-swap/bijection/union/cap/overflow/trailing-byte mutation suite",
        "One event kind is at most 256 UTF-8 bytes",
        "one complete record including envelope and LF is at most 1048576 bytes",
        "durable flush is at most 1024 rows and 4194304 bytes",
        "one run/scope is at most 65536 rows and 67108864 bytes",
        "governor retains at most 268435456 bytes across concurrent runs",
        "worst-case encoded bytes including envelope, redaction and framing overhead",
        "exact non-borrowable reservation",
        "independently serviced bounded writer",
        "Provisional->Released for refused admission",
        "Provisional->Reserved->Consumed|Released for admitted capacity",
        "finalization emits one reservation-settlement root",
        "crash recovery replays the durable ledger to the same settlement without double consumption or leak",
        "checked global earliest-deadline-first demand-bound test",
        "sum of remaining worst-case durable service time of all already admitted and proposed priority segments",
        "Per-run segment_count*worst_case_durable_service_time is necessary but never sufficient",
        "priority traffic from concurrent runs cannot queue cancellation observations",
        "byte, row or service-time overcommit",
        "No created event is silently dropped",
        "ExecutionState=InfrastructureFailed and EvidenceCompleteness=PartialEvidence",
        "durable-service bound compatible with the exact observation, drain and finalization SLOs",
        "ACCESSIBILITY_AGENT_PARITY",
        "DETERMINISM_PROJECTION",
        "Raw wall/host timestamps",
        "never participate in G5 equality",
        "initially Unmeasured and no speed claim",
    ] {
        assert!(
            policy.contains(required),
            "campaign policy lost '{required}'"
        );
    }
    let future_artifact_encoding = "FutureArtifact=FRAME_BYTES(Utf8(slot_id))||FRAME_BYTES(Utf8(related_waiver_subject))||FRAME_BYTES(Utf8(artifact_role))||Digest(artifact_schema)||digest_union";
    assert_eq!(
        policy.match_indices(future_artifact_encoding).count(),
        1,
        "future output encoding must bind exactly one artifact-schema digest before Pending"
    );
    assert!(
        !policy.contains("FutureArtifact=FRAME_BYTES(Utf8(slot_id))||FRAME_BYTES(Utf8(related_waiver_subject))||FRAME_BYTES(Utf8(artifact_role))||digest_union"),
        "schema-free future output encoding must never regain authority"
    );
    let governed_matrix_encoding = "FRAME_UTF8('i07-governed-output-schema-matrix-v1')||U16LE(1)||U16LE(1)||U16LE(3)||FRAME_BYTES(Utf8('i07.governed-real-ap242-corpus-discharge'))||FRAME_BYTES(Utf8('i07-governed-real-ap242-corpus'))||FRAME_BYTES(Utf8('GovernedRealAp242CorpusDischarge/'))||U64LE(1)||U64LE(1)||U16LE(1)||U16LE(1)";
    assert_eq!(
        policy.match_indices(governed_matrix_encoding).count(),
        1,
        "the closed governed-output matrix must have one exact canonical encoding"
    );
    let matrix_header = "i07-governed-output-schema-matrix-v1";
    let authority_scope = "i07.governed-real-ap242-corpus-discharge";
    let governed_slot = "i07-governed-real-ap242-corpus";
    let role_prefix = "GovernedRealAp242CorpusDischarge/";
    let governed_role = "GovernedRealAp242CorpusDischarge/i07-governed-real-ap242-corpus";
    let exact_matrix_bytes = (8 + matrix_header.len())
        + 3 * 2
        + (8 + 8 + authority_scope.len())
        + (8 + 8 + governed_slot.len())
        + (8 + 8 + role_prefix.len())
        + 2 * 8
        + 2 * 2;
    assert_eq!(exact_matrix_bytes, 221);
    let exact_schema_record_bytes =
        2 * 2 + (8 + 8 + authority_scope.len()) + (8 + 8 + governed_role.len()) + 32;
    assert_eq!(exact_schema_record_bytes, 171);
    let exact_membership_proof_bytes = 32 + 2 * 8 + 8 + exact_schema_record_bytes;
    assert_eq!(exact_membership_proof_bytes, 227);
    let governed_schema_set_formula = "GovernedOutputSchemaSetRoot=BLAKE3::derive_key('org.frankensim.i07.governed-output-schema-set.v1',raw32(GovernedOutputSchemaMatrixDigest)||U64LE(1)||FRAME_BYTES(GovernedOutputSchemaRecordV1))";
    assert_eq!(
        policy.match_indices(governed_schema_set_formula).count(),
        1,
        "the singleton governed-output schema-set algebra must remain exact"
    );
    let envelope_schema_identity = "i07-governed-corpus-discharge-envelope-v1";
    let initiative = "I07";
    let exact_envelope_bytes = (8 + envelope_schema_identity.len())
        + 2
        + 34 * (2 + 8)
        + (6 * 8
            + envelope_schema_identity.len()
            + initiative.len()
            + authority_scope.len()
            + governed_slot.len()
            + governed_slot.len()
            + governed_role.len())
        + 2 * 2
        + 8
        + 25 * 32;
    assert_eq!(exact_envelope_bytes, 1458);
    for exact_schema_guard in [
        "role-addressed future envelope slots that bind each output artifact_schema_digest before Pending",
        "schema whose bytes are available and independently rehashed before realization",
        "realized envelope must decode under the byte-identical precommitted schema digest",
        "role or artifact-schema swap",
        "unavailable/mismatched future artifact schema",
    ] {
        assert_eq!(
            policy.match_indices(exact_schema_guard).count(),
            1,
            "future output schema guard must occur exactly once: {exact_schema_guard}"
        );
    }
    let intent_fields = [
        ("0x0001", "initiative_id", "Utf8"),
        ("0x0002", "schema_identity", "Utf8"),
        ("0x0003", "successor_version", "U64LE"),
        ("0x0004", "predecessor_manifest_digest", "Digest"),
        ("0x0005", "expected_authority_head", "AuthorityHead"),
        ("0x0006", "predecessor_stage_receipt_digest", "Digest"),
        ("0x0007", "candidate_freeze_commitment_digest", "Digest"),
        ("0x0008", "retired_waiver_subjects", "Utf8List"),
        ("0x0009", "protected_bindings", "ProtectedBindingList"),
        ("0x000a", "coupled_transaction_group", "Utf8"),
        ("0x000b", "governance_stage", "U16LE"),
        ("0x000c", "authority_scope", "Utf8"),
        ("0x000d", "mutation_fence", "MutationFence"),
        ("0x000e", "envelope_slots", "FutureArtifactList"),
        ("0x000f", "final_successor_slot", "FutureDigest"),
        ("0x0010", "amendment_record_slot", "FutureDigest"),
        ("0x0011", "receipt_schema_set_root", "Digest"),
        ("0x0012", "profile_registry_root", "Digest"),
        ("0x0013", "profile_assignment_root", "Digest"),
        ("0x0014", "governed_output_schema_set_root", "Digest"),
    ];
    assert_eq!(intent_fields.len(), 20);
    let intent_header = b"I07_TRANSACTION_INTENT_V1\0";
    assert_eq!(intent_header.len(), 26);
    assert_eq!(intent_header.last(), Some(&0));
    assert!(!intent_header.ends_with(br"\0"));
    assert!(!policy.contains(r"I07_TRANSACTION_INTENT_V1\0"));
    let mut previous_mapping_end = 0;
    for (tag, name, field_type) in intent_fields {
        let exact_mapping = format!("{tag} {name} {field_type}");
        assert_eq!(
            policy.match_indices(exact_mapping.as_str()).count(),
            1,
            "transaction-intent mapping must occur exactly once: {exact_mapping}"
        );
        let mapping_start = policy
            .find(exact_mapping.as_str())
            .expect("exact transaction-intent mapping");
        assert!(
            mapping_start >= previous_mapping_end,
            "transaction-intent mappings must remain in increasing tag order: {exact_mapping}"
        );
        previous_mapping_end = mapping_start + exact_mapping.len();
    }
    let profile_key_fields = [
        (1, "part21_encoding_edition_id", "Utf8"),
        (2, "part21_encoding_artifact_digest", "Digest"),
        (3, "ap242_application_protocol_id", "Utf8"),
        (4, "ap242_edition_id", "Utf8"),
        (5, "ap242_schema_artifact_digest", "Digest"),
        (6, "conformance_class_id", "Utf8"),
        (7, "implementation_method_id", "Utf8"),
        (8, "character_encoding_id", "Utf8"),
        (9, "signature_policy_id", "Utf8"),
        (10, "external_reference_policy_id", "Utf8"),
        (11, "supported_entity_inventory_digest", "Digest"),
        (12, "supported_semantic_atom_inventory_digest", "Digest"),
        (13, "validation_rule_inventory_digest", "Digest"),
        (14, "allowed_reference_scc_digest", "Digest"),
        (15, "unit_frame_policy_digest", "Digest"),
        (16, "canonical_writer_policy_digest", "Digest"),
        (17, "no_claim_policy_digest", "Digest"),
        (18, "coverage_floor_digest", "Digest"),
        (19, "limit_policy_digest", "Digest"),
    ];
    let mut previous_profile_field_end = 0;
    for (ordinal, name, field_type) in profile_key_fields {
        let exact_mapping = format!("{ordinal} {name} {field_type}");
        assert_eq!(
            policy.match_indices(exact_mapping.as_str()).count(),
            1,
            "ProfileKey mapping must occur exactly once: {exact_mapping}"
        );
        let mapping_start = policy
            .find(exact_mapping.as_str())
            .expect("exact ProfileKey mapping");
        assert!(
            mapping_start >= previous_profile_field_end,
            "ProfileKey mappings must remain in increasing ordinal order: {exact_mapping}"
        );
        previous_profile_field_end = mapping_start + exact_mapping.len();
    }
    let assignment_key_fields = [
        "1 leaf_id Utf8",
        "2 claim_or_domain_id Utf8",
        "3 target_kind U32LE",
        "4 target_id Utf8",
        "5 profile_role U32LE",
    ];
    let mut previous_assignment_field_end = 0;
    for exact_mapping in assignment_key_fields {
        assert_eq!(
            policy.match_indices(exact_mapping).count(),
            1,
            "AssignmentLogicalKey mapping must occur exactly once: {exact_mapping}"
        );
        let mapping_start = policy
            .find(exact_mapping)
            .expect("exact AssignmentLogicalKey mapping");
        assert!(
            mapping_start >= previous_assignment_field_end,
            "AssignmentLogicalKey mappings must remain in increasing ordinal order"
        );
        previous_assignment_field_end = mapping_start + exact_mapping.len();
    }
    let exact_intent_table = exact_section(
        policy,
        "fields in strictly increasing order: ",
        ". Digest is exactly 32 nonzero raw bytes",
    );
    assert_eq!(
        normalized_table_heads(exact_intent_table),
        intent_fields
            .iter()
            .map(|(tag, name, field_type)| format!("{tag} {name} {field_type}"))
            .collect::<Vec<_>>(),
        "the transaction table must contain exactly the twenty declared fields and no hidden extension"
    );
    let mut exact_intent_entries = intent_fields
        .iter()
        .map(|(tag, name, field_type)| format!("{tag} {name} {field_type}"))
        .collect::<Vec<_>>();
    exact_intent_entries[10].push_str(
        " where GovernanceCommitted=1,CandidateFrozen=2,RealizationCommitted=3,RevealedForAdjudication=4,Closed=5",
    );
    assert_eq!(
        normalized_table_entries(exact_intent_table),
        exact_intent_entries,
        "the transaction table must not hide trailing per-field grammar"
    );
    assert!(
        !policy.contains("then U16LE(19), then exactly nineteen"),
        "the attacker-selectable nineteen-field intent must not regain authority"
    );
    let governed_envelope_schema_fields = [
        "1 schema_identity Utf8 EqUtf8('i07-governed-corpus-discharge-envelope-v1')",
        "2 initiative_id Utf8 EqUtf8('I07')",
        "3 protocol_kind U16LE EqU16(1)",
        "4 governance_stage U16LE EqU16(3)",
        "5 authority_scope Utf8 EqUtf8('i07.governed-real-ap242-corpus-discharge')",
        "6 slot_id Utf8 EqUtf8('i07-governed-real-ap242-corpus')",
        "7 related_waiver_subject Utf8 EqUtf8('i07-governed-real-ap242-corpus')",
        "8 artifact_role Utf8 EqUtf8('GovernedRealAp242CorpusDischarge/i07-governed-real-ap242-corpus')",
        "9 transaction_intent_digest Digest EqTransactionIntentDigest",
        "10 precommit_authorization_digest Digest EqPrecommitAuthorizationDigest",
        "11 governed_output_schema_set_root Digest EqIntentTag(0x0014)",
        "12 profile_registry_root Digest EqIntentTag(0x0012)",
        "13 profile_assignment_root Digest EqIntentTag(0x0013)",
        "14 standard_application_protocol_part21_identity_root Digest EqAuthorizationBinding('standard_application_protocol_part21_identity_root')",
        "15 conformance_class_root Digest EqAuthorizationBinding('conformance_class_root')",
        "16 license_access_custody_root Digest EqAuthorizationBinding('license_access_custody_root')",
        "17 corpus_item_root Digest EqAuthorizationBinding('corpus_item_root')",
        "18 corpus_item_count U64LE EqAuthorizationBoundU64Range('corpus_item_count',1,1048576)",
        "19 expectation_root Digest EqAuthorizationBinding('expectation_root')",
        "20 semantic_inventory_root Digest EqAuthorizationBinding('semantic_inventory_root')",
        "21 qoi_band_root Digest EqAuthorizationBinding('qoi_band_root')",
        "22 clause_interoperability_matrix_root Digest EqAuthorizationBinding('clause_interoperability_matrix_root')",
        "23 public_governed_blind_authority_root Digest EqAuthorizationBinding('public_governed_blind_authority_root')",
        "24 candidate_root Digest EqCandidateFrozenComponent(1)",
        "25 checker_root Digest EqCandidateFrozenComponent(2)",
        "26 toolchain_root Digest EqCandidateFrozenComponent(3)",
        "27 decision_root Digest EqCandidateFrozenComponent(4)",
        "28 joint_candidate_checker_expectation_root Digest EqAuthorizationBinding('joint_candidate_checker_expectation_root')",
        "29 oracle_tcb_threat_graph_root Digest EqAuthorizationBinding('oracle_tcb_threat_graph_root')",
        "30 independent_review_receipt_root Digest EqAuthorizationBinding('independent_review_receipt_root')",
        "31 revocation_receipt_root Digest EqAuthorizationBinding('revocation_receipt_root')",
        "32 semantic_loss_and_component_result_root Digest EqAuthorizationBinding('semantic_loss_and_component_result_root')",
        "33 redaction_policy_root Digest EqAuthorizationBinding('redaction_policy_root')",
        "34 validation_axis_root Digest EqAuthorizationBinding('validation_axis_root')",
    ];
    let exact_governed_envelope_schema_table = exact_section(
        policy,
        "The exact field table is: ",
        ". I07GovernedCorpusDischargeEnvelopeBytesV1",
    );
    assert_eq!(
        normalized_table_entries(exact_governed_envelope_schema_table),
        governed_envelope_schema_fields
            .iter()
            .map(|entry| (*entry).to_owned())
            .collect::<Vec<_>>(),
        "the governed AP242 discharge envelope schema must remain closed and exact"
    );
    assert_eq!(
        normalized_table_heads(exact_governed_envelope_schema_table).len(),
        34,
        "the governed output schema cannot add or omit a field"
    );
    assert_eq!(
        governed_envelope_schema_fields
            .iter()
            .filter(|field| {
                field.contains("EqAuthorizationBinding(")
                    || field.contains("EqAuthorizationBoundU64Range(")
            })
            .count(),
        17,
        "the precommit-authorization binding namespace must stay closed to seventeen fields"
    );
    for forbidden in [
        "GovernedRealAp242CorpusDischargeEnvelope/",
        "governedrealap242corpusdischarge/",
        "authority_scope='i07.governed-real-ap242-corpus'",
        "min_count=0",
        "max_count=4096",
        "schema-set count exactly 0",
    ] {
        assert!(
            !policy.contains(forbidden),
            "permissive governed output-schema alias regained authority: {forbidden}"
        );
    }
    let exact_profile_key_table = exact_section(
        policy,
        "followed by exactly these nineteen TAGGED_FIELD entries, once each in increasing ordinal order: ",
        ". Utf8 payloads use FRAME_UTF8",
    );
    assert_eq!(
        normalized_table_heads(exact_profile_key_table),
        profile_key_fields
            .iter()
            .map(|(ordinal, name, field_type)| format!("{ordinal} {name} {field_type}"))
            .collect::<Vec<_>>(),
        "the ProfileKey table must contain exactly nineteen fields and no hidden extension"
    );
    assert_eq!(
        normalized_table_entries(exact_profile_key_table),
        profile_key_fields
            .iter()
            .map(|(ordinal, name, field_type)| format!("{ordinal} {name} {field_type}"))
            .collect::<Vec<_>>(),
        "the ProfileKey table must not hide trailing per-field grammar"
    );
    let exact_assignment_key_table = exact_section(
        policy,
        "AssignmentLogicalKey=(leaf_id,claim_or_domain_id,target_kind,target_id,profile_role), encoded as exactly ",
        ", each as one TAGGED_FIELD in increasing ordinal order",
    );
    assert_eq!(
        normalized_table_heads(exact_assignment_key_table),
        assignment_key_fields
            .iter()
            .map(|entry| (*entry).to_owned())
            .collect::<Vec<_>>(),
        "the AssignmentLogicalKey table must contain exactly five fields and no hidden extension"
    );
    assert_eq!(
        normalized_table_entries(exact_assignment_key_table),
        assignment_key_fields
            .iter()
            .map(|entry| (*entry).to_owned())
            .collect::<Vec<_>>(),
        "the AssignmentLogicalKey table must not hide trailing per-field grammar"
    );
    let parser_activation = frozen
        .claim("i07-part21-profile-bounded-parse")
        .expect("parser claim")
        .activation;
    assert!(parser_activation.contains("PROFILE_INSTANCE_GATE"));
    assert!(parser_activation.contains("KAT_GATE"));
    assert!(parser_activation.contains("ProfileDefinitionSchemaRoot"));
    assert!(parser_activation.contains("ProfileKeyDigest/ProfileEntryRoot"));
    assert!(parser_activation.contains("required-assignment-universe"));
    for required in [
        "schema_id is the exact 23-byte ASCII string i07-profile-registry-v1",
        "canonical decoder refuses a missing, extra, duplicate, unknown, out-of-order, case-folded, normalized, or trailing schema/profile field",
        "part21_encoding_edition_id",
        "part21_encoding_artifact_digest",
        "ap242_application_protocol_id",
        "ap242_edition_id",
        "ap242_schema_artifact_digest",
        "conformance_class_id",
        "implementation_method_id",
        "character_encoding_id",
        "signature_policy_id",
        "external_reference_policy_id",
        "supported_entity_inventory_digest",
        "supported_semantic_atom_inventory_digest",
        "validation_rule_inventory_digest",
        "allowed_reference_scc_digest",
        "unit_frame_policy_digest",
        "canonical_writer_policy_digest",
        "no_claim_policy_digest",
        "coverage_floor_digest",
        "limit_policy_digest",
        "allowed reference SCC classes",
        "PROFILE_INSTANCE_GATE",
        "ProfileKeyDigest",
        "ProfileDefinitionSchemaBytes",
        "ProfileDefinitionSchemaRoot",
        "ProfileEntryRoot",
        "ProfileRegistryRoot",
        "RequiredAssignmentUniverseRoot",
        "ProfileAssignmentRoot",
        "total immutable mapping",
        "claim_or_domain_id",
        "AssignmentTargetKind",
        "CorpusItem",
        "ProfileRole",
        "SourceEndpoint",
        "duplicate/missing/unexpected logical keys",
        "repeated entries hidden by ordinal aliases",
        "favorable replay-time choice",
        "required claim/leaf set",
        "target/role assignment matrix",
        "nonvacuity matrix",
        "KAT_GATE vector",
        "--profile <profile-root> is exactly the assigned ProfileEntryRoot",
        "schema/registry/universe/assignment membership proofs",
        "no code infers edition/profile from a friendly name",
        "Real ISO artifacts are not embedded here",
    ] {
        assert!(
            profile.contains(required),
            "profile registry lost '{required}'"
        );
    }

    assert_eq!(waiver.subject, REAL_CORPUS_SLOT);
    assert_eq!(
        waiver.owner,
        "I07 manifest/corpus owner frankensim-leapfrog-2026-program-i94v.2.2.8.1 and the named G2/G5 independent-reproduction adjudicators"
    );
    assert!(
        waiver
            .reason
            .contains("cannot be invented or embedded here")
    );
    for required in [
        "I07GovernedCorpusDischargeReceiptV1",
        "exact waiver subject and predicate",
        "one atomic FrozenManifest::amend authority transaction",
        "same-ID I07GovernedCorpusDischargeEnvelopeV1 External root",
        "removes this Waiver row",
        "verifies the AmendmentRecord",
        "advances the authority head",
        "exact standard/application-protocol/Part-21 edition",
        "schema/profile/conformance-class artifact digest",
        "license/access/custody",
        "semantic-role profile assignment",
        "corpus item and expectation root",
        "QoI/band",
        "standard clause pack or interoperability matrix",
        "public-versus-governed-blind authority",
        "joint candidate/checker/expectation commitment",
        "oracle/TCB threat graph",
        "independent review and revocation receipts",
        "predecessor manifest",
        "exact transaction intent",
        "verified AmendmentRecord separately binds the final successor",
        "atomic transaction advances the authority head",
        "raw/all-zero digest",
        "fixture-only replacement",
        "receipt-only change",
        "waiver-only removal",
        "split transaction is not discharge",
    ] {
        assert!(
            waiver.predicate.contains(required),
            "governed external waiver predicate lost '{required}'"
        );
    }
    for required in [
        "before any protected corpus byte, aggregate, derived label or expectation is accessed",
        "before activation of i07-governed-real-ap242-corpus-validation",
        "review again for every corpus item or expectation root",
        "standard/application-protocol/Part-21 edition",
        "schema/profile/conformance class",
        "license/access/custody",
        "semantic-role assignment",
        "QoI/band",
        "clause pack",
        "interoperability matrix",
        "candidate/checker/toolchain/decision",
        "oracle/TCB",
        "reviewer/revocation",
        "envelope/receipt",
        "transaction-intent or authority-head change",
    ] {
        assert!(
            waiver.expiry.contains(required),
            "governed external waiver expiry lost '{required}'"
        );
    }
    assert!(
        waiver
            .promotion_effect
            .contains("only i07-governed-real-ap242-corpus-validation")
    );
    assert!(
        waiver
            .promotion_effect
            .contains("no synthetic result may be relabeled as real-corpus")
    );

    let governed_activation = frozen
        .claim("i07-governed-real-ap242-corpus-validation")
        .expect("governed external claim")
        .activation;
    for required in [
        "one atomic fs-vvreg authority transaction",
        "same-ID typed external discharge-envelope root",
        "removes the external-corpus waiver",
        "verifies the AmendmentRecord",
        "advances the authority head",
        "only then grants one-shot capabilities",
    ] {
        assert!(
            governed_activation.contains(required),
            "governed claim activation lost '{required}'"
        );
    }

    let explicits = frozen.explicits();
    assert!(explicits.seeds.contains("never IID, secrecy"));
    assert!(explicits.seeds.contains("M0=0xD2511F53"));
    assert!(explicits.seeds.contains("M1=0xCD9E8D57"));
    assert!(explicits.seeds.contains("W0=0x9E3779B9"));
    assert!(explicits.seeds.contains("W1=0xBB67AE85"));
    assert!(
        explicits
            .seeds
            .contains("k0=LE32(d[0..4]); k1=LE32(d[4..8])")
    );
    assert!(explicits.seeds.contains("d=BLAKE3::derive_key"));
    assert!(explicits.seeds.contains("c0=low32(case_index)"));
    assert!(explicits.seeds.contains("c3=high32(output_block_ordinal)"));
    assert!(
        explicits
            .seeds
            .contains("Lane is the Philox output-word index")
    );
    assert!(
        explicits
            .seeds
            .contains("byte offset n selects block floor(n/16)")
    );
    assert!(explicits.seeds.contains("Native-endian casts"));
    assert!(explicits.seeds.contains("KAT_GATE"));
    for required in [
        "development and heldout-format alias/digest/key/counter/four-output-word vectors",
        "every counter word and output lane",
        "exact LE serialized bytes",
        "nonzero high32 case and block words",
        "block 0-to-1 and byte-offset 15-to-16 boundaries",
        "case-sensitive/UTF-8 alias twins",
        "Cross-alias derived-key/counter collisions fail admission",
    ] {
        assert!(
            explicits.seeds.contains(required),
            "I07 KAT gate lost '{required}'"
        );
    }
    assert!(
        explicits
            .seeds
            .contains("version-1 prose does not invent those computed words")
    );
    assert!(explicits.units.contains("exact receipt/reconciliation"));
    assert!(explicits.units.contains("use unit 'bit'"));
    assert!(
        explicits
            .units
            .contains("numeric geometry/error defect scores use unit '1'")
    );
    assert!(
        explicits
            .budgets
            .contains("proof declaration/checker steps <=256")
    );
    assert!(explicits.budgets.contains("request observation <=1 s"));
    assert!(
        explicits
            .budgets
            .contains("governed-corpus work <=256 items or <=4096 semantic atoms")
    );
    assert!(
        explicits
            .budgets
            .contains("terminate-drain-finalize within the same bound")
    );
    assert!(
        explicits
            .versions
            .contains("material/process/lot/reference envelope schema v1")
    );
    assert!(
        explicits
            .capabilities
            .contains("frontier-ap242-semantic-roundtrip")
    );
    assert!(
        explicits
            .capabilities
            .contains("moonshot-ap242-proof-carrying-semantics")
    );
    assert!(!explicits.capabilities.contains("semantic-lossless"));
    assert!(
        explicits
            .capabilities
            .contains("no network, FFI, foreign code")
    );
    assert!(
        explicits
            .capabilities
            .contains("no licensed-text embedding")
    );
    assert!(
        explicits
            .capabilities
            .contains("drawing approval, or regulatory authority")
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn i07_fixture_partitions_and_seed_ranges_are_stage_separated() {
    let frozen = i07_draft().freeze().expect("freeze");
    let seed_ids: BTreeSet<_> = SEED_FIXTURES.iter().map(|(id, _, _)| *id).collect();
    let unseeded_ids = [
        CAMPAIGN_POLICY_FIXTURE,
        PROFILE_REGISTRY_FIXTURE,
        "i07-semantic-descent-theorem-card",
        "i07-semantic-equivalence-theorem-card",
    ]
    .into_iter()
    .collect::<BTreeSet<_>>();
    assert_eq!(
        seed_ids.len(),
        SEED_FIXTURES.len(),
        "duplicate seed fixture id"
    );
    assert_eq!(
        seed_ids
            .union(&unseeded_ids)
            .copied()
            .collect::<BTreeSet<_>>(),
        frozen
            .fixtures()
            .iter()
            .map(|pin| pin.id)
            .collect::<BTreeSet<_>>(),
        "every fixture must be explicitly classified as seeded or governance/theorem metadata"
    );
    let held_out: BTreeSet<_> = frozen
        .fixtures()
        .iter()
        .filter(|pin| pin.partition == Partition::HeldOut)
        .map(|pin| pin.id)
        .collect();
    assert_eq!(
        held_out,
        HELD_OUT_FIXTURES.iter().copied().collect::<BTreeSet<_>>()
    );
    assert_eq!(
        HELD_OUT_AUTHORITY
            .iter()
            .map(|(id, _, _)| *id)
            .collect::<BTreeSet<_>>(),
        held_out
    );

    for &(held_out_id, expected_tier, expected_consumer) in HELD_OUT_AUTHORITY {
        let consumers: Vec<_> = frozen
            .obligations()
            .iter()
            .filter(|row| row.decks().contains(&held_out_id))
            .collect();
        assert_eq!(
            consumers.len(),
            1,
            "held-out fixture '{held_out_id}' must have one named consumer"
        );
        assert_eq!(consumers[0].leaf(), expected_consumer);
        assert_eq!(consumers[0].tier(), expected_tier);
        let spec = authored_spec(fixture(frozen.fixtures(), held_out_id));
        assert!(spec.contains("STAGE-HOLDOUT"));
        assert!(spec.contains("Public deterministic replay"));
        assert!(
            spec.contains(&format!("Sole consumer {expected_consumer}")),
            "held-out fixture '{held_out_id}' has stale textual consumer authority"
        );
    }

    for (id, alias, range) in SEED_FIXTURES {
        let pin = fixture(frozen.fixtures(), id);
        let spec = authored_spec(pin);
        assert!(
            spec.contains(&format!("alias '{alias}'")),
            "seed alias drift for {id}"
        );
        assert!(
            spec.contains(&format!("k={range}")),
            "seed range drift for {id}"
        );
        match pin.partition {
            Partition::Development => assert_eq!(*range, "0..=4095"),
            Partition::HeldOut => assert!(
                matches!(*range, "65536..=69631" | "131072..=135167"),
                "held-out range overlap for {id}: {range}"
            ),
        }
    }

    for row in frozen.obligations() {
        assert!(row.decks().contains(&CAMPAIGN_POLICY_FIXTURE));
        assert!(row.decks().contains(&PROFILE_REGISTRY_FIXTURE));
    }

    let waiver_consumers: Vec<_> = frozen
        .obligations()
        .iter()
        .filter(|row| row.decks().contains(&REAL_CORPUS_SLOT))
        .map(|row| row.leaf())
        .collect();
    assert_eq!(
        waiver_consumers,
        ["i07-governed-real-corpus-validation-max"]
    );
    for row in frozen.obligations() {
        if row.leaf() != "i07-governed-real-corpus-validation-max" {
            assert!(
                !row.decks().contains(&REAL_CORPUS_SLOT),
                "synthetic leaf '{}' must not inherit external-corpus waiver authority",
                row.leaf()
            );
        }
    }

    let multi_holdout_leaves: BTreeSet<_> = frozen
        .obligations()
        .iter()
        .filter_map(|row| {
            let heldout_count = row
                .decks()
                .iter()
                .filter(|deck| held_out.contains(*deck))
                .count();
            (heldout_count > 1).then_some(row.leaf())
        })
        .collect();
    assert_eq!(
        multi_holdout_leaves,
        BTreeSet::from([
            "i07-harness-loss-core",
            "i07-pmi-material-core",
            "i07-structure-geometry-core",
        ])
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn i07_obligations_have_complete_execution_contracts() {
    let frozen = i07_draft().freeze().expect("freeze");
    assert_eq!(
        COMMON_OBS_EVENTS
            .iter()
            .copied()
            .collect::<BTreeSet<_>>()
            .len(),
        COMMON_OBS_EVENTS.len(),
        "the I07 common event contract itself must not contain aliases"
    );
    for row in frozen.obligations() {
        assert!(row.entry_point().starts_with("scripts/e2e/leapfrog/i07_"));
        assert!(row.entry_point().ends_with(".sh"));
        assert_eq!(
            row.dsr_lane(),
            format!(
                "env FRANKENSIM_VMANIFEST_LEAF={} dsr quality --tool frankensim",
                row.leaf()
            ),
            "{} must name an exact executable DSR slice",
            row.leaf()
        );
        for replay_fragment in [row.entry_point(), "--manifest", "--profile", "--replay"] {
            assert!(
                row.replay_command().contains(replay_fragment),
                "replay contract for '{}' lost '{replay_fragment}'",
                row.leaf()
            );
        }
        for g0_fragment in ["GENERATORS:", "VALIDITY:", "LAWS:", "SHRINKERS:", "REPLAY"] {
            assert!(
                row.g0().contains(g0_fragment),
                "{} lost {g0_fragment}",
                row.leaf()
            );
        }
        for g4_fragment in ["Cancel", "fault", "Drain", "finalize", "publish"] {
            assert!(
                row.g4_schedule()
                    .to_ascii_lowercase()
                    .contains(&g4_fragment.to_ascii_lowercase()),
                "{} lost G4 concept {g4_fragment}",
                row.leaf()
            );
        }
        for g5_fragment in ["aarch64-apple", "x86_64", "deterministic", "Cross"] {
            assert!(
                row.g5_matrix()
                    .to_ascii_lowercase()
                    .contains(&g5_fragment.to_ascii_lowercase()),
                "{} lost G5 axis {g5_fragment}",
                row.leaf()
            );
        }
        assert!(
            row.g5_matrix().contains("threads") || row.g5_matrix().contains("workers"),
            "{} lost its concurrency axis",
            row.leaf()
        );
        for event in COMMON_OBS_EVENTS {
            assert_eq!(
                row.obs_events()
                    .iter()
                    .filter(|observed| *observed == event)
                    .count(),
                1,
                "{} must declare lifecycle/common event {event} exactly once",
                row.leaf()
            );
        }
        assert!(!row.g3_relations().is_empty());
    }

    let governed = frozen
        .obligations()
        .iter()
        .find(|row| row.leaf() == "i07-governed-real-corpus-validation-max")
        .expect("governed corpus leaf");
    for fragment in [
        "semantic-role assignment",
        "item-addressable",
        "duplicate entries under ordinal aliases refuse",
        "ScopedClauseConformance",
        "ScopedRequiredMatrixValidated",
        "GovernedBlindScopedValidation",
        "arbitrary-digest",
        "same-ID-fixture-only",
        "receipt-only",
        "waiver-only",
        "split-transaction",
        "typed I07GovernedCorpusDischargeEnvelopeV1 and receipt",
        "atomic amendment/authority-head transition",
        "generic manifest syntax or partial discharge grants no authority",
        "predecessor/successor manifests",
        "transaction intent",
        "authority head",
    ] {
        assert!(
            governed.g0().contains(fragment),
            "governed G0 lost '{fragment}'"
        );
    }
    for fragment in [
        "arbitrary same-ID External digest",
        "removing its waiver",
        "typed envelope",
        "verified receipt",
        "one atomic authority-head transaction",
        "cannot preserve discharge authority",
    ] {
        assert!(
            governed
                .g3_relations()
                .iter()
                .any(|relation| relation.contains(fragment)),
            "governed G3 relations lost '{fragment}'"
        );
    }
    for fragment in [
        "at most 256 corpus items or 4096 semantic atoms per cancellable batch",
        "discharge-envelope/receipt verification",
        "atomic amendment and authority-head commit",
        "every independent reader/import/export heartbeat",
        "terminate-drain-finalize",
        "request observation <=1 s",
        "After the primary cancellation no new capability issue/use, protected read, reveal, adjudication, publication or authority-head advance may begin",
        "only an operation whose own durable commit point preceded that cancellation remains committed",
        "in-flight precommit reconciles conservatively without performing a new effect",
        "partial discharge",
    ] {
        assert!(
            governed.g4_schedule().contains(fragment),
            "governed G4 lost '{fragment}'"
        );
    }
    assert!(
        governed
            .g5_matrix()
            .contains("canonical logical-event/adjudication projection")
    );
    assert!(
        governed
            .g5_matrix()
            .contains("Raw timestamps, durations, scheduler/process telemetry")
    );
    assert!(governed.g5_matrix().contains("never for bitwise equality"));

    let parser_writer = frozen
        .obligations()
        .iter()
        .find(|row| row.leaf() == "i07-parser-writer-core")
        .expect("parser/writer leaf");
    assert!(parser_writer.g4_schedule().contains("TimedOut"));
    assert!(!parser_writer.g4_schedule().contains("TimeOut"));
    assert!(
        parser_writer
            .g5_matrix()
            .contains("canonical-logical-event-order")
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn i07_load_bearing_semantic_clauses_are_pinned() {
    let frozen = i07_draft().freeze().expect("freeze");
    let claim_text = |id: &str| {
        let claim = frozen
            .claim(id)
            .unwrap_or_else(|| panic!("missing claim '{id}'"));
        let mut text = String::from(claim.statement);
        for hypothesis in claim.hypotheses {
            text.push(' ');
            text.push_str(hypothesis);
        }
        text.push(' ');
        text.push_str(claim.kill);
        text.push(' ');
        text.push_str(claim.no_claim);
        text
    };

    let parser = claim_text("i07-part21-profile-bounded-parse");
    assert!(parser.contains("ProfileKeyDigest/ProfileEntryRoot"));
    assert!(parser.contains("is the complete expanded authority identity"));
    assert!(parser.contains("WHERE/UNIQUE rule inventory"));
    assert!(parser.contains("legal strongly connected reference components"));
    assert!(parser.contains("cycles are never rejected merely because they are cycles"));

    let assembly = claim_text("i07-product-assembly-configuration-lineage");
    assert!(assembly.contains("independent of exchange-file #numbering"));
    assert!(assembly.contains("not traversal order or geometry proximity"));
    assert!(assembly.contains("typed lineage morphisms or refusal"));

    let geometry = claim_text("i07-representation-qualified-geometry-roundtrip");
    assert!(geometry.contains("tessellated input is never relabeled exact"));
    assert!(geometry.contains("sampled distance or visual overlays are falsifiers only"));
    assert!(geometry.contains("absent premises yield Unknown"));
    let geometry_fixture = authored_spec(fixture(
        frozen.fixtures(),
        "i07-mixed-exact-tessellated-geometry",
    ));
    assert!(geometry_fixture.contains("s_H=max(d_source_to_target,d_target_to_source)"));
    assert!(
        geometry_fixture.contains("independently certified upper endpoints, never sample maxima")
    );
    assert!(geometry_fixture.contains("s_J=||J_source-J_target||_F"));
    assert!(geometry_fixture.contains("outward-rounded interval/ball enclosures"));
    assert!(geometry_fixture.contains(
        "Every A,V,J, reference scale, absolute difference, maximum, norm, sum and quotient"
    ));
    assert!(geometry_fixture.contains("one dependency-aware outward-rounded expression"));
    assert!(geometry_fixture.contains("certified finite upper endpoints"));
    assert!(geometry_fixture.contains("never point estimates"));
    assert!(geometry_fixture.contains("no V/J score for an open shell"));
    assert!(geometry_fixture.contains("GEOMETRY_ACCEPTANCE is the conjunction"));
    assert!(geometry_fixture.contains("finite s_geometry<=1"));
    assert!(geometry_fixture.contains("self-intersection-freedom"));
    assert!(geometry_fixture.contains("certificate state blocks promotion"));
    assert!(!geometry_fixture.contains("floors+"));

    let pmi = claim_text("i07-semantic-pmi-gdt-datum-texture");
    assert!(pmi.contains("semantic PMI and presentation-only graphical PMI are disjoint types"));
    assert!(pmi.contains("no drawing approval"));

    let material = claim_text("i07-material-process-lot-requirement-evidence");
    assert!(material.contains("without inventing inaccessible property values"));
    assert!(material.contains("silent URI fetch"));

    let harness = claim_text("i07-harness-ewis-connectivity-semantics");
    assert!(harness.contains("separate edge families"));
    assert!(harness.contains("no inferred connectivity"));
    assert!(harness.contains("no ampacity, voltage-drop, EMC"));

    let loss = claim_text("i07-exhaustive-semantic-loss-receipt");
    assert!(loss.contains("each source atom has exactly one DispositionRecord"));
    assert!(loss.contains("each target atom has exactly one OriginRecord"));
    assert!(loss.contains("typed zero/one/many counterpart set"));
    assert!(
        loss.contains("forward and inverse source-target incidence relations are exact transposes")
    );
    assert!(loss.contains("split, merge, normalization, synthesis, and deletion"));
    assert!(loss.contains("PreservedExact"));
    assert!(loss.contains("DroppedWithRefusal"));
    assert!(loss.contains("orthogonal terminal states"));
    assert!(loss.contains("legitimately empty admitted inventories"));

    let canonical = claim_text("i07-canonical-export-determinism-and-containment");
    assert!(canonical.contains("content root are deterministic"));
    assert!(canonical.contains("performs no network access or foreign execution"));
    assert!(canonical.contains("transactional temporary artifact"));
    assert!(canonical.contains("reimports to the frozen semantic receipt"));

    let opaque = claim_text("i07-safe-opaque-extension-capsules");
    assert!(opaque.contains("SafeOpaqueCapsule"));
    assert!(opaque.contains("explicit NoSemanticAuthority"));
    assert!(opaque.contains("never proves vendor meaning"));

    let kinematics = claim_text("i07-kinematic-pair-semantic-roundtrip");
    assert!(kinematics.contains("configuration relation and tangent/Pfaffian relation"));
    assert!(kinematics.contains("sampled motion"));
    assert!(kinematics.contains("no mobility completeness"));

    let migration = claim_text("i07-partial-bidirectional-edition-migration");
    assert!(migration.contains("partial bidirectional lens"));
    assert!(migration.contains("GetPut, PutGet, and PutPut laws"));
    assert!(migration.contains("rather than inventing an inverse"));

    let governed = claim_text("i07-governed-real-ap242-corpus-validation");
    assert!(governed.contains("ProfileKeyDigest/ProfileEntryRoot"));
    assert!(governed.contains("fixture/case/corpus-item/stratum"));
    assert!(governed.contains("semantic profile role"));
    assert!(governed.contains("ScopedClauseConformance"));
    assert!(governed.contains("ScopedRequiredMatrixValidated"));
    assert!(governed.contains("GovernedBlindScopedValidation"));
    assert!(governed.contains("PublicCorpusOnly"));
    assert!(governed.contains("repeated entry disguised by another ordinal"));
    assert!(governed.contains("public-as-blind relabeling"));
    assert!(governed.contains("no universal AP242 conformance"));

    let equivalence = claim_text("i07-proof-carrying-bidirectional-semantic-equivalence");
    assert!(equivalence.contains("form an equivalence"));
    assert!(equivalence.contains("typed directed-error/loss quantaloid"));
    assert!(equivalence.contains("homs are complete lattices in a frozen authority order"));
    assert!(equivalence.contains("preserves arbitrary joins in each variable"));
    assert!(equivalence.contains("without treating a finite tolerance as equality"));
    assert!(
        equivalence
            .contains("tolerance closeness is never silently made symmetric, transitive, or equal")
    );
    assert!(equivalence.contains("{propext, Quot.sound, Classical.choice}"));
    assert!(equivalence.contains("manifest version 1 prose mints no theorem color"));

    let descent = claim_text("i07-semantic-sheaf-descent-obstruction-theorem");
    assert!(descent.contains("strict local sections glue uniquely"));
    assert!(descent.contains("complete nerve"));
    assert!(descent.contains("every nonempty finite intersection"));
    assert!(descent.contains("pair/triple/higher-overlap"));
    assert!(descent.contains("groupoid/stack descent"));
    assert!(descent.contains("retaining stabilizers/automorphisms"));
    assert!(descent.contains("uniqueness up to gauge is not ordinary sheaf uniqueness"));
    assert!(descent.contains("closed non-exact cochain witness"));
    assert!(descent.contains("a mismatch is not automatically H1"));
    assert!(descent.contains("nonabelian/categorical descent"));
    assert!(descent.contains("independently checked d^2=0"));
    assert!(descent.contains("never laundered into vector-space H1"));
    assert!(descent.contains("pinned-cover Cech H1"));
    assert!(descent.contains("derived sheaf H1 requires"));

    let falsifier = claim_text("i07-semantic-equivalence-counterexample-search");
    assert!(falsifier.contains("N=16,777,216"));
    assert!(falsifier.contains("manifest-version-1 prose has no exhaustive authority"));
    assert!(falsifier.contains("empty bounded search is not a proof"));

    for theorem_fixture in [
        "i07-semantic-equivalence-theorem-card",
        "i07-semantic-descent-theorem-card",
    ] {
        let spec = authored_spec(fixture(frozen.fixtures(), theorem_fixture));
        assert!(
            spec.contains("Version-1 prose mints") || spec.contains("manifest version 1 is prose")
        );
        assert!(spec.contains("FORMAL_PROJECTION_GATE"));
    }
    let descent_card = authored_spec(fixture(
        frozen.fixtures(),
        "i07-semantic-descent-theorem-card",
    ));
    assert!(descent_card.contains("SHEAF_H1_RATCHET"));
    assert!(descent_card.contains("Leray/acyclic-cover comparison"));
    assert!(descent_card.contains("conclusion remains cover-relative"));
    assert!(descent_card.contains("finite-semantic-sheaf-and-stack-descent"));
    assert!(descent_card.contains("complete nerve with every nonempty finite intersection"));
    assert!(descent_card.contains("pair/triple/higher-overlap"));
    assert!(descent_card.contains("effective normalization quotient"));
    assert!(descent_card.contains("groupoid/stack object-and-morphism descent"));
    assert!(descent_card.contains("automorphism structure"));
    let grammar = authored_spec(fixture(
        frozen.fixtures(),
        "i07-semantic-adversary-microgrammar",
    ));
    assert!(grammar.contains("16^6=N=16,777,216"));
    assert!(grammar.contains("M0_FORMALIZATION_GATE"));
    assert!(grammar.contains(
        "version-1 prose has no cardinality, enumeration, quotient, or exhaustive authority"
    ));
}

#[test]
fn i07_identity_is_stable_and_all_set_order_is_invariant() {
    let frozen = i07_draft().freeze().expect("freeze");
    assert_eq!(
        frozen.digest(),
        i07_draft().freeze().expect("refreeze").digest()
    );

    let mut reordered = i07_draft();
    reordered.claims.reverse();
    reordered.fixtures.reverse();
    reordered.obligations.reverse();
    reordered.waivers.reverse();
    for row in &mut reordered.obligations {
        row.claims_covered = reversed_static(row.claims_covered);
        row.unit_cases = reversed_static(row.unit_cases);
        row.decks = reversed_static(row.decks);
        row.g3_relations = reversed_static(row.g3_relations);
        row.obs_events = reversed_static(row.obs_events);
    }
    assert_eq!(
        reordered.freeze().expect("reordered freeze").digest(),
        frozen.digest()
    );
}

#[test]
fn i07_chunked_in_memory_assembly_matches_one_shot_freeze() {
    // I07 G4 precursor only: this checks completed in-memory content identity.
    // It does not claim process restart, serialized checkpoints, corruption
    // detection, cancellation, draining, or runtime fault containment.
    let one_shot = i07_draft();
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
    for chunk in one_shot.fixtures.chunks(4) {
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
    let chunked = staged.freeze().expect("chunked in-memory freeze");
    assert_eq!(chunked.digest(), expected.digest());
    assert_eq!(chunked, expected);
}

#[test]
#[allow(clippy::too_many_lines)]
fn i07_g3_mutations_refuse_or_move_frozen_authority() {
    let baseline = i07_draft().freeze().expect("freeze").digest();

    let mut weakened = i07_draft();
    weakened
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i07-semantic-sheaf-descent-obstruction-theorem")
        .expect("descent claim")
        .hypotheses = &["local labels happen to match"];
    assert_ne!(
        weakened.freeze().expect("weakened successor").digest(),
        baseline
    );

    let mut correlated = i07_draft();
    correlated
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i07-exhaustive-semantic-loss-receipt")
        .expect("loss claim")
        .oracle
        .independent = false;
    assert!(matches!(
        correlated.freeze(),
        Err(FreezeRefusal::ProductionOracleReuse { .. })
    ));

    let mut rerouted_oracle = i07_draft();
    rerouted_oracle
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i07-harness-ewis-connectivity-semantics")
        .expect("harness claim")
        .oracle
        .identity = "i07.oracle.harness.v1 test-only alternate route";
    assert_ne!(
        rerouted_oracle
            .freeze()
            .expect("oracle-rerouted successor")
            .digest(),
        baseline,
        "an oracle-route change must move frozen authority"
    );

    let mut relaxed = i07_draft();
    relaxed
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i07-representation-qualified-geometry-roundtrip")
        .expect("geometry claim")
        .tolerance = ToleranceSemantics::Interval { lo: 0.0, hi: 2.0 };
    assert_ne!(
        relaxed.freeze().expect("relaxed successor").digest(),
        baseline
    );

    let mut swapped_holdout = i07_draft();
    swapped_holdout
        .fixtures
        .iter_mut()
        .find(|pin| pin.id == "i07-semantic-adversary-max-holdout")
        .expect("adversary holdout")
        .partition = Partition::Development;
    assert_ne!(
        swapped_holdout
            .freeze()
            .expect("swapped holdout successor")
            .digest(),
        baseline
    );

    let mut missing_development_deck = i07_draft();
    assert_eq!(
        fixture(
            &missing_development_deck.fixtures,
            "i07-part21-synthetic-corpus",
        )
        .partition,
        Partition::Development,
        "the mutation target must remain a Development deck"
    );
    missing_development_deck
        .obligations
        .iter_mut()
        .find(|row| row.leaf == "i07-parser-writer-core")
        .expect("parser/writer leaf")
        .decks = &[
        CAMPAIGN_POLICY_FIXTURE,
        PROFILE_REGISTRY_FIXTURE,
        "i07-part21-core-holdout",
        "i07-adversarial-reference-graphs",
        "i07-semantic-loss-unsupported-corpus",
    ];
    assert_ne!(
        missing_development_deck
            .freeze()
            .expect("development-deck-removed successor")
            .digest(),
        baseline,
        "removing a leaf's Development deck must move frozen authority"
    );

    let mut missing_policy = i07_draft();
    missing_policy
        .fixtures
        .retain(|pin| pin.id != CAMPAIGN_POLICY_FIXTURE);
    assert!(matches!(
        missing_policy.freeze(),
        Err(FreezeRefusal::OrphanDeck { deck, .. }) if deck == CAMPAIGN_POLICY_FIXTURE
    ));
}

#[test]
#[allow(clippy::too_many_lines)]
fn i07_amendments_invalidate_exact_transitive_authority() {
    let frozen = i07_draft().freeze().expect("freeze");

    let mut version_only = i07_draft();
    version_only.version = 2;
    let (version_successor, version_record) =
        frozen.amend(version_only).expect("version-only amendment");
    assert!(version_record.invalidated.is_empty());
    assert_eq!(
        (version_record.from_version, version_record.to_version),
        (1, 2)
    );
    assert_eq!(version_record.from_digest, frozen.digest());
    assert_eq!(version_record.to_digest, version_successor.digest());
    assert_ne!(version_record.from_digest, version_record.to_digest);

    let mut changed_pmi = i07_draft();
    changed_pmi.version = 2;
    changed_pmi
        .claims
        .iter_mut()
        .find(|claim| claim.id == "i07-semantic-pmi-gdt-datum-texture")
        .expect("PMI claim")
        .statement = "test-only changed PMI semantics";
    let (_, pmi_record) = frozen.amend(changed_pmi).expect("PMI amendment");
    assert_eq!(
        pmi_record.invalidated,
        vec![
            "i07-pmi-material-core",
            "i07-semantic-pmi-gdt-datum-texture"
        ]
    );

    let mut changed_harness_fixture = i07_draft();
    changed_harness_fixture.version = 2;
    changed_harness_fixture
        .fixtures
        .iter_mut()
        .find(|pin| pin.id == "i07-harness-ewis-graphs")
        .expect("harness fixture")
        .source = FixtureSource::AuthoredSpec {
        spec: "test-only changed harness graph generator semantics",
    };
    let (_, harness_record) = frozen
        .amend(changed_harness_fixture)
        .expect("harness fixture amendment");
    assert_eq!(
        harness_record.invalidated,
        vec![
            "i07-exhaustive-semantic-loss-receipt",
            "i07-harness-ewis-connectivity-semantics",
            "i07-harness-loss-core",
        ]
    );

    let mut changed_policy = i07_draft();
    changed_policy.version = 2;
    changed_policy
        .fixtures
        .iter_mut()
        .find(|pin| pin.id == CAMPAIGN_POLICY_FIXTURE)
        .expect("campaign policy")
        .source = FixtureSource::AuthoredSpec {
        spec: "I07_CAMPAIGN_POLICY_V2 test-only unauthorized authority change",
    };
    let (_, policy_record) = frozen.amend(changed_policy).expect("policy amendment");
    assert_eq!(
        policy_record
            .invalidated
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        all_authority_ids()
    );
    assert_eq!(policy_record.invalidated.len(), 25);

    let mut changed_profile_registry = i07_draft();
    changed_profile_registry.version = 2;
    changed_profile_registry
        .fixtures
        .iter_mut()
        .find(|pin| pin.id == PROFILE_REGISTRY_FIXTURE)
        .expect("profile registry")
        .source = FixtureSource::AuthoredSpec {
        spec: "I07_PROFILE_REGISTRY_V2 test-only changed global profile authority",
    };
    let (_, profile_record) = frozen
        .amend(changed_profile_registry)
        .expect("profile-registry amendment");
    assert_eq!(
        profile_record
            .invalidated
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        all_authority_ids(),
        "the profile registry is a global authority dependency"
    );

    let mut changed_waiver = i07_draft();
    changed_waiver.version = 2;
    changed_waiver.waivers[0].expiry = "test-only unauthorized indefinite waiver";
    let (_, waiver_record) = frozen
        .amend(changed_waiver)
        .expect("real-corpus waiver amendment");
    assert_eq!(
        waiver_record
            .invalidated
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "i07-governed-real-ap242-corpus-validation",
            "i07-governed-real-corpus-validation-max",
        ]),
        "the real-corpus waiver must invalidate only its scoped claim and sole leaf"
    );

    // The generic manifest layer can freeze this syntactically valid successor,
    // but the I07 policy deliberately grants it no scientific authority. Only
    // fs-vvreg can verify the typed discharge and atomic authority transition.
    let mut structural_only_successor = i07_draft();
    structural_only_successor.version = 2;
    structural_only_successor.waivers.clear();
    structural_only_successor.fixtures.push(FixturePin {
        id: REAL_CORPUS_SLOT,
        source: FixtureSource::External {
            digest_hex: "0000000000000000000000000000000000000000000000000000000000000000",
        },
        partition: Partition::HeldOut,
    });
    let (slot_successor, slot_record) = frozen
        .amend(structural_only_successor)
        .expect("structural-only same-id external deck-slot replacement");
    assert!(slot_successor.waivers().is_empty());
    let replaced = fixture(slot_successor.fixtures(), REAL_CORPUS_SLOT);
    assert_eq!(replaced.partition, Partition::HeldOut);
    assert!(matches!(
        replaced.source,
        FixtureSource::External { digest_hex }
            if digest_hex == "0000000000000000000000000000000000000000000000000000000000000000"
    ));
    assert_eq!(
        slot_record
            .invalidated
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "i07-governed-real-ap242-corpus-validation",
            "i07-governed-real-corpus-validation-max",
        ]),
        "real-corpus slot replacement must touch no synthetic authority"
    );
    let successor_policy =
        authored_spec(fixture(slot_successor.fixtures(), CAMPAIGN_POLICY_FIXTURE));
    for boundary in [
        "raw or all-zero digest",
        "arbitrary same-ID External fixture",
        "waiver removal alone",
        "structurally frozen successor",
        "without the verified fs-vvreg transaction grants no discharge or promotion authority",
    ] {
        assert!(
            successor_policy.contains(boundary),
            "structural-only successor lost no-authority boundary '{boundary}'"
        );
    }
    let governed_activation = slot_successor
        .claim("i07-governed-real-ap242-corpus-validation")
        .expect("governed external claim")
        .activation;
    assert!(
        governed_activation.contains("one atomic fs-vvreg authority transaction")
            && governed_activation.contains("verifies the AmendmentRecord")
            && governed_activation.contains("advances the authority head"),
        "generic structural freeze must not erase the external authority gate"
    );

    let mut changed_explicits = i07_draft();
    changed_explicits.version = 2;
    changed_explicits.explicits.capabilities =
        "test-only changed capability authority; every descendant must invalidate";
    let (_, explicits_record) = frozen
        .amend(changed_explicits)
        .expect("Five Explicits amendment");
    assert_eq!(
        explicits_record
            .invalidated
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        all_authority_ids()
    );

    let mut changed_title = i07_draft();
    changed_title.version = 2;
    changed_title.title = "test-only changed I07 campaign semantics";
    let (_, title_record) = frozen.amend(changed_title).expect("title amendment");
    assert_eq!(
        title_record
            .invalidated
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        all_authority_ids()
    );
}
