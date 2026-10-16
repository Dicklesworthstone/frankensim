//! Machine-IR graph-bound fit/clearance admission (Gauntlet G0/G3/G5).

use core::num::NonZeroU64;

use fs_blake3::ContentHash;
use fs_ir::machine::manufacturing::ManufacturingArtifactRefV1;
use fs_ir::machine::manufacturing::fit_clearance::{
    FitAllowanceErrorV1, FitAllowanceV1, FitEndpointRoleV1, FitFeatureSelectorV1, FitLengthUnitV1,
    FitPairTargetV1, FitPresentationRefV1, FitRegimeV1, FitRequirementIdV1, FitRequirementV1,
    FitSemanticSourceRefV1, FitSpecificationRefV1, MAX_MACHINE_FIT_REQUIREMENTS_V1,
    MachineFitClearanceAdmissionErrorV1, MachineFitClearanceDraftV1, PositiveFitLengthErrorV1,
    PositiveFitLengthV1, SignedFitLengthErrorV1, SignedFitLengthV1,
};
use fs_ir::machine::{
    AdmittedMachineGraph, BodyId, ContactFeatureId, MachineGraphDraft, MaterialBinding,
    MaterialCardRef, MaterialTarget, ModelRef, SubsystemId, SubsystemSpec,
};

fn nz(value: u64) -> NonZeroU64 {
    NonZeroU64::new(value).expect("fixture version is nonzero")
}

fn body(key: &str) -> BodyId {
    BodyId::new(key).expect("fixture body key is canonical")
}

fn feature(key: &str) -> ContactFeatureId {
    ContactFeatureId::new(key).expect("fixture contact-feature key is canonical")
}

fn requirement_id(key: &str) -> FitRequirementIdV1 {
    FitRequirementIdV1::new(key).expect("fixture requirement key is canonical")
}

fn signed(value: f64, unit: FitLengthUnitV1) -> SignedFitLengthV1 {
    SignedFitLengthV1::try_new(value, unit).expect("fixture signed length is admitted")
}

fn positive(value: f64, unit: FitLengthUnitV1) -> PositiveFitLengthV1 {
    PositiveFitLengthV1::try_new(value, unit).expect("fixture positive length is admitted")
}

fn allowance(minimum: f64, maximum: f64, unit: FitLengthUnitV1) -> FitAllowanceV1 {
    FitAllowanceV1::try_new(signed(minimum, unit), signed(maximum, unit))
        .expect("fixture allowance is admitted")
}

fn artifact(namespace: &str, byte: u8) -> ManufacturingArtifactRefV1 {
    ManufacturingArtifactRefV1::new(namespace, nz(1), ContentHash([byte; 32]))
        .expect("fixture artifact coordinate is canonical")
}

fn specification(byte: u8) -> FitSpecificationRefV1 {
    FitSpecificationRefV1::new(artifact("fit/specification", byte))
}

fn semantic_source(byte: u8) -> FitSemanticSourceRefV1 {
    FitSemanticSourceRefV1::new(artifact("fit/semantic-source", byte))
}

fn presentation(byte: u8) -> FitPresentationRefV1 {
    FitPresentationRefV1::new(artifact("fit/presentation", byte))
}

fn selector(body_key: &str, feature_key: &str) -> FitFeatureSelectorV1 {
    FitFeatureSelectorV1::new(body(body_key), feature(feature_key))
}

fn material(target: BodyId, key: &str, byte: u8) -> MaterialBinding {
    MaterialBinding {
        target: MaterialTarget::Body(target),
        material: MaterialCardRef::new(key, nz(1), [byte; 32])
            .expect("fixture material is canonical"),
    }
}

fn admitted_graph(model_byte: u8) -> AdmittedMachineGraph {
    let internal = body("body/internal");
    let internal_alternate = body("body/internal-alternate");
    let external = body("body/external");
    let external_alternate = body("body/external-alternate");
    MachineGraphDraft {
        clocks: Vec::new(),
        subsystems: vec![
            SubsystemSpec {
                id: SubsystemId::new("subsystem/internal").expect("canonical subsystem"),
                model: ModelRef::new("models/fit-internal", nz(1), [model_byte; 32])
                    .expect("canonical model"),
                bodies: vec![internal.clone(), internal_alternate.clone()],
                surface_patches: Vec::new(),
                contact_features: vec![
                    feature("contact/internal/a"),
                    feature("contact/internal/b"),
                ],
                state_slots: Vec::new(),
            },
            SubsystemSpec {
                id: SubsystemId::new("subsystem/external").expect("canonical subsystem"),
                model: ModelRef::new("models/fit-external", nz(1), [0x52; 32])
                    .expect("canonical model"),
                bodies: vec![external.clone(), external_alternate.clone()],
                surface_patches: Vec::new(),
                contact_features: vec![
                    feature("contact/external/a"),
                    feature("contact/external/b"),
                ],
                state_slots: Vec::new(),
            },
        ],
        terminals: Vec::new(),
        ports: Vec::new(),
        relations: Vec::new(),
        materials: vec![
            material(internal, "materials/internal", 1),
            material(internal_alternate, "materials/internal-alternate", 2),
            material(external, "materials/external", 3),
            material(external_alternate, "materials/external-alternate", 4),
        ],
        interfaces: Vec::new(),
    }
    .admit()
    .expect("fit fixture graph admits")
}

#[allow(clippy::too_many_arguments)]
fn requirement(
    id: &str,
    internal_body: &str,
    internal_feature: &str,
    external_body: &str,
    external_feature: &str,
    basic_size: PositiveFitLengthV1,
    allowance: FitAllowanceV1,
    specification_byte: u8,
    source_byte: u8,
    presentation_byte: Option<u8>,
) -> FitRequirementV1 {
    FitRequirementV1::new(
        requirement_id(id),
        FitPairTargetV1::new(
            selector(internal_body, internal_feature),
            selector(external_body, external_feature),
        ),
        basic_size,
        allowance,
        specification(specification_byte),
        semantic_source(source_byte),
        presentation_byte.map(presentation),
    )
}

fn valid_requirement(id: &str, regime: FitRegimeV1, seed: u8) -> FitRequirementV1 {
    let (gaps, internal_feature, external_feature) = match regime {
        FitRegimeV1::Clearance => (
            allowance(10.0, 30.0, FitLengthUnitV1::Micrometre),
            "contact/internal/a",
            "contact/external/a",
        ),
        FitRegimeV1::Transition => (
            allowance(-5.0, 10.0, FitLengthUnitV1::Micrometre),
            "contact/internal/a",
            "contact/external/b",
        ),
        FitRegimeV1::Interference => (
            allowance(-30.0, -5.0, FitLengthUnitV1::Micrometre),
            "contact/internal/b",
            "contact/external/a",
        ),
    };
    requirement(
        id,
        "body/internal",
        internal_feature,
        "body/external",
        external_feature,
        positive(10.0, FitLengthUnitV1::Millimetre),
        gaps,
        seed,
        seed.wrapping_add(1),
        Some(seed.wrapping_add(2)),
    )
}

#[derive(Clone)]
struct RequirementFixture {
    id: &'static str,
    internal_body: &'static str,
    internal_feature: &'static str,
    external_body: &'static str,
    external_feature: &'static str,
    basic_size: PositiveFitLengthV1,
    allowance: FitAllowanceV1,
    specification_byte: u8,
    source_byte: u8,
    presentation_byte: Option<u8>,
}

impl RequirementFixture {
    fn baseline() -> Self {
        Self {
            id: "fit/single",
            internal_body: "body/internal",
            internal_feature: "contact/internal/a",
            external_body: "body/external",
            external_feature: "contact/external/a",
            basic_size: positive(10.0, FitLengthUnitV1::Millimetre),
            allowance: allowance(10.0, 30.0, FitLengthUnitV1::Micrometre),
            specification_byte: 0x61,
            source_byte: 0x62,
            presentation_byte: Some(0x63),
        }
    }

    fn build(self) -> FitRequirementV1 {
        requirement(
            self.id,
            self.internal_body,
            self.internal_feature,
            self.external_body,
            self.external_feature,
            self.basic_size,
            self.allowance,
            self.specification_byte,
            self.source_byte,
            self.presentation_byte,
        )
    }
}

fn singleton(requirement: FitRequirementV1) -> MachineFitClearanceDraftV1 {
    MachineFitClearanceDraftV1 {
        requirements: vec![requirement],
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn mfc_001_regimes_roles_units_and_semantic_fields_reach_identity() {
    let graph = admitted_graph(0x41);
    let draft = MachineFitClearanceDraftV1 {
        requirements: vec![
            valid_requirement("fit/transition", FitRegimeV1::Transition, 0x41),
            valid_requirement("fit/interference", FitRegimeV1::Interference, 0x44),
            valid_requirement("fit/clearance", FitRegimeV1::Clearance, 0x47),
        ],
    };
    let baseline = draft
        .clone()
        .admit_against(&graph)
        .expect("three fit regimes admit");
    let mut reordered = draft;
    reordered.requirements.reverse();
    let reordered = reordered
        .admit_against(&graph)
        .expect("caller order is non-semantic");
    assert_eq!(baseline.graph(), graph.identity());
    assert_eq!(baseline.identity(), reordered.identity());
    assert_eq!(
        baseline.identity_receipt().canonical_preimage(),
        reordered.identity_receipt().canonical_preimage()
    );
    assert_eq!(
        baseline
            .requirements()
            .iter()
            .map(|entry| entry.id().canonical_key())
            .collect::<Vec<_>>(),
        ["fit/clearance", "fit/interference", "fit/transition"]
    );
    assert_eq!(
        baseline
            .requirements()
            .iter()
            .map(|entry| entry.allowance().regime())
            .collect::<Vec<_>>(),
        [
            FitRegimeV1::Clearance,
            FitRegimeV1::Interference,
            FitRegimeV1::Transition,
        ]
    );

    let metre = signed(1.0, FitLengthUnitV1::Metre);
    let millimetre = signed(1_000.0, FitLengthUnitV1::Millimetre);
    assert_eq!(metre.metres_bits(), millimetre.metres_bits());
    assert_ne!(metre, millimetre, "submitted unit remains identity-bearing");
    assert_eq!(
        signed(-0.0, FitLengthUnitV1::Micrometre),
        signed(0.0, FitLengthUnitV1::Micrometre)
    );
    assert_eq!(FitLengthUnitV1::Inch.symbol(), "in");

    let base_fixture = RequirementFixture::baseline();
    let base = singleton(base_fixture.clone().build())
        .admit_against(&graph)
        .expect("base fit admits");
    let mut variants = Vec::new();
    let mut changed = base_fixture.clone();
    changed.id = "fit/renamed";
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.internal_body = "body/internal-alternate";
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.internal_feature = "contact/internal/b";
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.external_body = "body/external-alternate";
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.external_feature = "contact/external/b";
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.basic_size = positive(0.01, FitLengthUnitV1::Metre);
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.allowance = allowance(10.0, 31.0, FitLengthUnitV1::Micrometre);
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.allowance = allowance(11.0, 30.0, FitLengthUnitV1::Micrometre);
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.allowance = allowance(-1.0, 30.0, FitLengthUnitV1::Micrometre);
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.specification_byte = 0x64;
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.source_byte = 0x65;
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.presentation_byte = None;
    variants.push(changed);
    let mut changed = base_fixture;
    (changed.internal_body, changed.external_body) = (changed.external_body, changed.internal_body);
    (changed.internal_feature, changed.external_feature) =
        (changed.external_feature, changed.internal_feature);
    variants.push(changed);
    for variant in variants {
        let admitted = singleton(variant.build())
            .admit_against(&graph)
            .expect("one-field identity mutation remains structurally admissible");
        assert_ne!(base.identity(), admitted.identity());
    }

    let changed_graph = singleton(RequirementFixture::baseline().build())
        .admit_against(&admitted_graph(0x42))
        .expect("same selectors admit against changed graph");
    assert_ne!(base.identity(), changed_graph.identity());
}

#[test]
#[allow(clippy::too_many_lines)]
fn mfc_002_constructors_and_graph_admission_refuse_invalid_fit_state() {
    assert_eq!(
        SignedFitLengthV1::try_new(f64::NAN, FitLengthUnitV1::Metre),
        Err(SignedFitLengthErrorV1::NonFinite)
    );
    assert_eq!(
        SignedFitLengthV1::try_new(f64::from_bits(1), FitLengthUnitV1::Nanometre),
        Err(SignedFitLengthErrorV1::SiUnderflow)
    );
    assert_eq!(
        PositiveFitLengthV1::try_new(0.0, FitLengthUnitV1::Millimetre),
        Err(PositiveFitLengthErrorV1::NonPositive)
    );
    assert_eq!(
        PositiveFitLengthV1::try_new(-1.0, FitLengthUnitV1::Millimetre),
        Err(PositiveFitLengthErrorV1::NonPositive)
    );
    assert_eq!(
        FitAllowanceV1::try_new(
            signed(1.0, FitLengthUnitV1::Micrometre),
            signed(-1.0, FitLengthUnitV1::Micrometre),
        ),
        Err(FitAllowanceErrorV1::Inverted)
    );
    assert_eq!(
        FitAllowanceV1::try_new(
            signed(0.0, FitLengthUnitV1::Micrometre),
            signed(-0.0, FitLengthUnitV1::Micrometre),
        ),
        Err(FitAllowanceErrorV1::DegenerateZero)
    );
    assert_eq!(
        allowance(0.0, 1.0, FitLengthUnitV1::Micrometre).regime(),
        FitRegimeV1::Clearance
    );
    assert_eq!(
        allowance(-1.0, 0.0, FitLengthUnitV1::Micrometre).regime(),
        FitRegimeV1::Interference
    );
    assert_eq!(
        allowance(-1.0, 1.0, FitLengthUnitV1::Micrometre).regime(),
        FitRegimeV1::Transition
    );

    let graph = admitted_graph(0x41);
    assert_eq!(
        MachineFitClearanceDraftV1 {
            requirements: Vec::new(),
        }
        .admit_against(&graph),
        Err(MachineFitClearanceAdmissionErrorV1::NoRequirements)
    );

    let duplicate = RequirementFixture::baseline().build();
    assert_eq!(
        MachineFitClearanceDraftV1 {
            requirements: vec![duplicate.clone(), duplicate],
        }
        .admit_against(&graph),
        Err(MachineFitClearanceAdmissionErrorV1::DuplicateRequirement {
            requirement: requirement_id("fit/single"),
        })
    );

    let reversed_pair = MachineFitClearanceDraftV1 {
        requirements: vec![
            RequirementFixture::baseline().build(),
            requirement(
                "fit/reversed",
                "body/external-alternate",
                "contact/external/a",
                "body/internal-alternate",
                "contact/internal/a",
                positive(10.0, FitLengthUnitV1::Millimetre),
                allowance(-5.0, 5.0, FitLengthUnitV1::Micrometre),
                0x64,
                0x65,
                None,
            ),
        ],
    };
    assert_eq!(
        reversed_pair.admit_against(&graph),
        Err(MachineFitClearanceAdmissionErrorV1::DuplicateFeaturePair {
            first: requirement_id("fit/reversed"),
            duplicate: requirement_id("fit/single"),
        })
    );

    let same_body = requirement(
        "fit/same-body",
        "body/internal",
        "contact/internal/a",
        "body/internal",
        "contact/internal/b",
        positive(10.0, FitLengthUnitV1::Millimetre),
        allowance(1.0, 2.0, FitLengthUnitV1::Micrometre),
        0x61,
        0x62,
        None,
    );
    assert_eq!(
        singleton(same_body).admit_against(&graph),
        Err(MachineFitClearanceAdmissionErrorV1::SameBody {
            requirement: requirement_id("fit/same-body"),
            body: body("body/internal"),
        })
    );

    let same_feature = requirement(
        "fit/same-feature",
        "body/internal",
        "contact/internal/a",
        "body/internal-alternate",
        "contact/internal/a",
        positive(10.0, FitLengthUnitV1::Millimetre),
        allowance(1.0, 2.0, FitLengthUnitV1::Micrometre),
        0x61,
        0x62,
        None,
    );
    assert_eq!(
        singleton(same_feature).admit_against(&graph),
        Err(MachineFitClearanceAdmissionErrorV1::SameFeature {
            requirement: requirement_id("fit/same-feature"),
            feature: feature("contact/internal/a"),
        })
    );

    let invalid_selectors = [
        (
            requirement(
                "fit/unknown-internal-body",
                "body/missing",
                "contact/internal/a",
                "body/external",
                "contact/external/a",
                positive(10.0, FitLengthUnitV1::Millimetre),
                allowance(1.0, 2.0, FitLengthUnitV1::Micrometre),
                0x61,
                0x62,
                None,
            ),
            MachineFitClearanceAdmissionErrorV1::UnknownBody {
                requirement: requirement_id("fit/unknown-internal-body"),
                role: FitEndpointRoleV1::Internal,
                body: body("body/missing"),
            },
        ),
        (
            requirement(
                "fit/unknown-internal-feature",
                "body/internal",
                "contact/missing-internal",
                "body/external",
                "contact/external/a",
                positive(10.0, FitLengthUnitV1::Millimetre),
                allowance(1.0, 2.0, FitLengthUnitV1::Micrometre),
                0x61,
                0x62,
                None,
            ),
            MachineFitClearanceAdmissionErrorV1::UnknownFeature {
                requirement: requirement_id("fit/unknown-internal-feature"),
                role: FitEndpointRoleV1::Internal,
                feature: feature("contact/missing-internal"),
            },
        ),
        (
            requirement(
                "fit/unknown-external-body",
                "body/internal",
                "contact/internal/a",
                "body/missing-external",
                "contact/external/a",
                positive(10.0, FitLengthUnitV1::Millimetre),
                allowance(1.0, 2.0, FitLengthUnitV1::Micrometre),
                0x61,
                0x62,
                None,
            ),
            MachineFitClearanceAdmissionErrorV1::UnknownBody {
                requirement: requirement_id("fit/unknown-external-body"),
                role: FitEndpointRoleV1::External,
                body: body("body/missing-external"),
            },
        ),
        (
            requirement(
                "fit/unknown-external-feature",
                "body/internal",
                "contact/internal/a",
                "body/external",
                "contact/missing",
                positive(10.0, FitLengthUnitV1::Millimetre),
                allowance(1.0, 2.0, FitLengthUnitV1::Micrometre),
                0x61,
                0x62,
                None,
            ),
            MachineFitClearanceAdmissionErrorV1::UnknownFeature {
                requirement: requirement_id("fit/unknown-external-feature"),
                role: FitEndpointRoleV1::External,
                feature: feature("contact/missing"),
            },
        ),
    ];
    for (candidate, expected) in invalid_selectors {
        assert_eq!(singleton(candidate).admit_against(&graph), Err(expected));
    }

    let internal_owner_mismatch = requirement(
        "fit/internal-owner-mismatch",
        "body/internal",
        "contact/external/b",
        "body/external",
        "contact/external/a",
        positive(10.0, FitLengthUnitV1::Millimetre),
        allowance(1.0, 2.0, FitLengthUnitV1::Micrometre),
        0x61,
        0x62,
        None,
    );
    assert!(matches!(
        singleton(internal_owner_mismatch).admit_against(&graph),
        Err(MachineFitClearanceAdmissionErrorV1::FeatureOwnerMismatch {
            requirement,
            role: FitEndpointRoleV1::Internal,
            body,
            feature,
            body_owner,
            feature_owner,
        }) if requirement == requirement_id("fit/internal-owner-mismatch")
            && body == crate::body("body/internal")
            && feature == crate::feature("contact/external/b")
            && body_owner == SubsystemId::new("subsystem/internal").unwrap()
            && feature_owner == SubsystemId::new("subsystem/external").unwrap()
    ));

    let external_owner_mismatch = requirement(
        "fit/external-owner-mismatch",
        "body/internal",
        "contact/internal/a",
        "body/external",
        "contact/internal/b",
        positive(10.0, FitLengthUnitV1::Millimetre),
        allowance(1.0, 2.0, FitLengthUnitV1::Micrometre),
        0x61,
        0x62,
        None,
    );
    assert!(matches!(
        singleton(external_owner_mismatch).admit_against(&graph),
        Err(MachineFitClearanceAdmissionErrorV1::FeatureOwnerMismatch {
            requirement,
            role: FitEndpointRoleV1::External,
            body,
            feature,
            body_owner,
            feature_owner,
        }) if requirement == requirement_id("fit/external-owner-mismatch")
            && body == crate::body("body/external")
            && feature == crate::feature("contact/internal/b")
            && body_owner == SubsystemId::new("subsystem/external").unwrap()
            && feature_owner == SubsystemId::new("subsystem/internal").unwrap()
    ));

    assert_eq!(
        MachineFitClearanceAdmissionErrorV1::NoRequirements.code(),
        "MachineFitNoRequirements"
    );
}

fn boundary_graph() -> AdmittedMachineGraph {
    let internal = body("body/boundary-internal");
    let external = body("body/boundary-external");
    let mut features = vec![feature("contact/boundary/internal")];
    features.extend(
        (0..MAX_MACHINE_FIT_REQUIREMENTS_V1)
            .map(|index| feature(&format!("contact/boundary/external-{index:04}"))),
    );
    MachineGraphDraft {
        clocks: Vec::new(),
        subsystems: vec![SubsystemSpec {
            id: SubsystemId::new("subsystem/boundary").expect("canonical subsystem"),
            model: ModelRef::new("models/fit-boundary", nz(1), [0x71; 32])
                .expect("canonical model"),
            bodies: vec![internal.clone(), external.clone()],
            surface_patches: Vec::new(),
            contact_features: features,
            state_slots: Vec::new(),
        }],
        terminals: Vec::new(),
        ports: Vec::new(),
        relations: Vec::new(),
        materials: vec![
            material(internal, "materials/boundary-internal", 0x72),
            material(external, "materials/boundary-external", 0x73),
        ],
        interfaces: Vec::new(),
    }
    .admit()
    .expect("exact-cap fit graph admits")
}

fn boundary_draft() -> MachineFitClearanceDraftV1 {
    MachineFitClearanceDraftV1 {
        requirements: (0..MAX_MACHINE_FIT_REQUIREMENTS_V1)
            .map(|index| {
                requirement(
                    &format!("fit/boundary-r{index:04}"),
                    "body/boundary-internal",
                    "contact/boundary/internal",
                    "body/boundary-external",
                    &format!("contact/boundary/external-{index:04}"),
                    positive(10.0, FitLengthUnitV1::Millimetre),
                    allowance(1.0, 2.0, FitLengthUnitV1::Micrometre),
                    0x74,
                    0x75,
                    Some(0x76),
                )
            })
            .collect(),
    }
}

#[test]
fn mfc_003_exact_resource_cap_admits_and_one_over_refuses_before_deduplication() {
    let graph = boundary_graph();
    let exact = boundary_draft();
    let admitted = exact
        .clone()
        .admit_against(&graph)
        .expect("exact fit-requirement cap admits");
    assert_eq!(
        admitted.requirements().len(),
        MAX_MACHINE_FIT_REQUIREMENTS_V1
    );

    let mut too_many = exact;
    let repeated = too_many.requirements[0].clone();
    too_many.requirements.push(repeated);
    assert_eq!(
        too_many.admit_against(&graph),
        Err(MachineFitClearanceAdmissionErrorV1::RequirementLimit {
            actual: MAX_MACHINE_FIT_REQUIREMENTS_V1 + 1,
            max: MAX_MACHINE_FIT_REQUIREMENTS_V1,
        })
    );
}

#[test]
fn mfc_004_identical_input_replays_the_complete_receipt() {
    let graph = admitted_graph(0x41);
    let first = singleton(RequirementFixture::baseline().build())
        .admit_against(&graph)
        .expect("first replay admits");
    let second = singleton(RequirementFixture::baseline().build())
        .admit_against(&graph)
        .expect("second replay admits");
    assert_eq!(first.identity(), second.identity());
    assert_eq!(first.identity_receipt(), second.identity_receipt());
    assert_eq!(
        first.identity_receipt().canonical_preimage(),
        second.identity_receipt().canonical_preimage()
    );
}
