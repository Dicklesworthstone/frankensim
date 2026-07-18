//! Machine-IR surface-texture admission (Gauntlet G0/G3/G5).

use core::num::NonZeroU64;

use fs_blake3::ContentHash;
use fs_ir::machine::manufacturing::ManufacturingArtifactRefV1;
use fs_ir::machine::manufacturing::surface_texture::{
    MAX_MACHINE_SURFACE_TEXTURE_OBSERVATIONS_V1, MAX_MACHINE_SURFACE_TEXTURE_REQUIREMENTS_V1,
    MachineSurfaceTextureAdmissionErrorV1, MachineSurfaceTextureDraftV1, SurfaceLayV1,
    SurfaceProductionRuleV1, SurfaceTextureCalibrationRefV1, SurfaceTextureFilterRefV1,
    SurfaceTextureLengthErrorV1, SurfaceTextureLengthUnitV1, SurfaceTextureLengthV1,
    SurfaceTextureLimitV1, SurfaceTextureMeasurementRefV1, SurfaceTextureMetricV1,
    SurfaceTextureObservationContextIssueV1, SurfaceTextureObservationIdV1,
    SurfaceTextureObservationV1, SurfaceTexturePresentationRefV1, SurfaceTextureRequirementIdV1,
    SurfaceTextureRequirementIssueV1, SurfaceTextureRequirementV1,
    SurfaceTextureSemanticSourceRefV1, SurfaceTextureStandardRefV1,
};
use fs_ir::machine::{
    AdmittedMachineGraph, BodyId, FrameBinding, MachineGraphDraft, MaterialBinding,
    MaterialCardRef, MaterialTarget, ModelRef, OrientationParity, SubsystemId, SubsystemSpec,
    SurfacePatchId,
};

fn nz(value: u64) -> NonZeroU64 {
    NonZeroU64::new(value).expect("fixture version is nonzero")
}

fn body(key: &str) -> BodyId {
    BodyId::new(key).expect("fixture body key is canonical")
}

fn patch(key: &str) -> SurfacePatchId {
    SurfacePatchId::new(key).expect("fixture surface key is canonical")
}

fn requirement_id(key: &str) -> SurfaceTextureRequirementIdV1 {
    SurfaceTextureRequirementIdV1::new(key).expect("fixture requirement key is canonical")
}

fn observation_id(key: &str) -> SurfaceTextureObservationIdV1 {
    SurfaceTextureObservationIdV1::new(key).expect("fixture observation key is canonical")
}

fn length(value: f64, unit: SurfaceTextureLengthUnitV1) -> SurfaceTextureLengthV1 {
    SurfaceTextureLengthV1::try_new(value, unit).expect("fixture length is admitted")
}

fn artifact(namespace: &str, byte: u8) -> ManufacturingArtifactRefV1 {
    ManufacturingArtifactRefV1::new(namespace, nz(1), ContentHash([byte; 32]))
        .expect("fixture artifact coordinate is canonical")
}

fn filter(byte: u8) -> SurfaceTextureFilterRefV1 {
    SurfaceTextureFilterRefV1::new(artifact("surface-texture/filter", byte))
}

fn standard(byte: u8) -> SurfaceTextureStandardRefV1 {
    SurfaceTextureStandardRefV1::new(artifact("surface-texture/standard", byte))
}

fn semantic_source(byte: u8) -> SurfaceTextureSemanticSourceRefV1 {
    SurfaceTextureSemanticSourceRefV1::new(artifact("surface-texture/semantic", byte))
}

fn presentation(byte: u8) -> SurfaceTexturePresentationRefV1 {
    SurfaceTexturePresentationRefV1::new(artifact("surface-texture/presentation", byte))
}

fn measurement(byte: u8) -> SurfaceTextureMeasurementRefV1 {
    SurfaceTextureMeasurementRefV1::new(artifact("surface-texture/measurement", byte))
}

fn calibration(byte: u8) -> SurfaceTextureCalibrationRefV1 {
    SurfaceTextureCalibrationRefV1::new(artifact("surface-texture/calibration", byte))
}

fn frame(orientation: OrientationParity) -> FrameBinding {
    FrameBinding::new("frame/surface-texture", orientation).expect("fixture frame is canonical")
}

fn material(target: BodyId, key: &str, byte: u8) -> MaterialBinding {
    MaterialBinding {
        target: MaterialTarget::Body(target),
        material: MaterialCardRef::new(key, nz(1), [byte; 32])
            .expect("fixture material is canonical"),
    }
}

fn admitted_graph(model_byte: u8) -> AdmittedMachineGraph {
    let part = body("body/part");
    let alternate = body("body/part-alternate");
    let other = body("body/other");
    MachineGraphDraft {
        clocks: Vec::new(),
        subsystems: vec![
            SubsystemSpec {
                id: SubsystemId::new("subsystem/part").expect("canonical subsystem"),
                model: ModelRef::new("models/surface-texture-part", nz(1), [model_byte; 32])
                    .expect("canonical model"),
                bodies: vec![part.clone(), alternate.clone()],
                surface_patches: vec![
                    patch("surface/part/primary"),
                    patch("surface/part/secondary"),
                ],
                contact_features: Vec::new(),
                state_slots: Vec::new(),
            },
            SubsystemSpec {
                id: SubsystemId::new("subsystem/other").expect("canonical subsystem"),
                model: ModelRef::new("models/surface-texture-other", nz(1), [0x52; 32])
                    .expect("canonical model"),
                bodies: vec![other.clone()],
                surface_patches: vec![patch("surface/other/primary")],
                contact_features: Vec::new(),
                state_slots: Vec::new(),
            },
        ],
        terminals: Vec::new(),
        ports: Vec::new(),
        relations: Vec::new(),
        materials: vec![
            material(part, "materials/part", 1),
            material(alternate, "materials/part-alternate", 2),
            material(other, "materials/other", 3),
        ],
        interfaces: Vec::new(),
    }
    .admit()
    .expect("surface-texture fixture graph admits")
}

#[allow(clippy::too_many_arguments)]
fn requirement(
    id: &str,
    declared_body: &str,
    surface: &str,
    metric: SurfaceTextureMetricV1,
    limit: SurfaceTextureLimitV1,
    cutoff: SurfaceTextureLengthV1,
    evaluation: SurfaceTextureLengthV1,
    filter_byte: u8,
    lay: SurfaceLayV1,
    lay_frame: Option<FrameBinding>,
    production_rule: SurfaceProductionRuleV1,
    standard_byte: u8,
    source_byte: u8,
    presentation_byte: Option<u8>,
) -> SurfaceTextureRequirementV1 {
    SurfaceTextureRequirementV1::new(
        requirement_id(id),
        body(declared_body),
        patch(surface),
        metric,
        limit,
        cutoff,
        evaluation,
        filter(filter_byte),
        lay,
        lay_frame,
        production_rule,
        standard(standard_byte),
        semantic_source(source_byte),
        presentation_byte.map(presentation),
    )
}

fn valid_requirement(
    id: &str,
    surface: &str,
    metric: SurfaceTextureMetricV1,
    seed: u8,
) -> SurfaceTextureRequirementV1 {
    requirement(
        id,
        "body/part",
        surface,
        metric,
        SurfaceTextureLimitV1::Maximum(length(0.8, SurfaceTextureLengthUnitV1::Micrometre)),
        length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
        length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
        0x31,
        SurfaceLayV1::Parallel,
        Some(frame(OrientationParity::Preserving)),
        SurfaceProductionRuleV1::MaterialRemovalRequired,
        seed,
        seed.wrapping_add(1),
        Some(seed.wrapping_add(2)),
    )
}

#[allow(clippy::too_many_arguments)]
fn observation(
    id: &str,
    requirement: &str,
    metric: SurfaceTextureMetricV1,
    cutoff: SurfaceTextureLengthV1,
    evaluation: SurfaceTextureLengthV1,
    filter_byte: u8,
    measured: SurfaceTextureLengthV1,
    uncertainty: SurfaceTextureLengthV1,
    measurement_byte: u8,
    calibration_byte: u8,
) -> SurfaceTextureObservationV1 {
    SurfaceTextureObservationV1::new(
        observation_id(id),
        requirement_id(requirement),
        metric,
        cutoff,
        evaluation,
        filter(filter_byte),
        measured,
        uncertainty,
        measurement(measurement_byte),
        calibration(calibration_byte),
    )
}

fn valid_observation(
    id: &str,
    requirement: &str,
    metric: SurfaceTextureMetricV1,
    seed: u8,
) -> SurfaceTextureObservationV1 {
    observation(
        id,
        requirement,
        metric,
        length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
        length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
        0x31,
        length(0.72, SurfaceTextureLengthUnitV1::Micrometre),
        length(0.02, SurfaceTextureLengthUnitV1::Micrometre),
        seed,
        seed.wrapping_add(1),
    )
}

fn valid_draft() -> MachineSurfaceTextureDraftV1 {
    MachineSurfaceTextureDraftV1 {
        requirements: vec![
            valid_requirement(
                "texture/b",
                "surface/part/secondary",
                SurfaceTextureMetricV1::Rz,
                0x41,
            ),
            valid_requirement(
                "texture/a",
                "surface/part/primary",
                SurfaceTextureMetricV1::Ra,
                0x44,
            ),
        ],
        observations: vec![
            valid_observation(
                "observation/b",
                "texture/b",
                SurfaceTextureMetricV1::Rz,
                0x51,
            ),
            valid_observation(
                "observation/a",
                "texture/a",
                SurfaceTextureMetricV1::Ra,
                0x54,
            ),
        ],
    }
}

#[derive(Clone)]
struct RequirementFixture {
    id: &'static str,
    declared_body: &'static str,
    surface: &'static str,
    metric: SurfaceTextureMetricV1,
    limit: SurfaceTextureLimitV1,
    cutoff: SurfaceTextureLengthV1,
    evaluation: SurfaceTextureLengthV1,
    filter_byte: u8,
    lay: SurfaceLayV1,
    lay_frame: Option<FrameBinding>,
    production_rule: SurfaceProductionRuleV1,
    standard_byte: u8,
    source_byte: u8,
    presentation_byte: Option<u8>,
}

impl RequirementFixture {
    fn baseline() -> Self {
        Self {
            id: "texture/single",
            declared_body: "body/part",
            surface: "surface/part/primary",
            metric: SurfaceTextureMetricV1::Ra,
            limit: SurfaceTextureLimitV1::Maximum(length(
                0.8,
                SurfaceTextureLengthUnitV1::Micrometre,
            )),
            cutoff: length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
            evaluation: length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
            filter_byte: 0x31,
            lay: SurfaceLayV1::Parallel,
            lay_frame: Some(frame(OrientationParity::Preserving)),
            production_rule: SurfaceProductionRuleV1::MaterialRemovalRequired,
            standard_byte: 0x61,
            source_byte: 0x62,
            presentation_byte: Some(0x63),
        }
    }

    fn build(self) -> SurfaceTextureRequirementV1 {
        requirement(
            self.id,
            self.declared_body,
            self.surface,
            self.metric,
            self.limit,
            self.cutoff,
            self.evaluation,
            self.filter_byte,
            self.lay,
            self.lay_frame,
            self.production_rule,
            self.standard_byte,
            self.source_byte,
            self.presentation_byte,
        )
    }
}

#[derive(Clone)]
struct ObservationFixture {
    id: &'static str,
    requirement: &'static str,
    metric: SurfaceTextureMetricV1,
    cutoff: SurfaceTextureLengthV1,
    evaluation: SurfaceTextureLengthV1,
    filter_byte: u8,
    measured: SurfaceTextureLengthV1,
    uncertainty: SurfaceTextureLengthV1,
    measurement_byte: u8,
    calibration_byte: u8,
}

impl ObservationFixture {
    fn baseline() -> Self {
        Self {
            id: "observation/observed",
            requirement: "texture/observed",
            metric: SurfaceTextureMetricV1::Ra,
            cutoff: length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
            evaluation: length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
            filter_byte: 0x31,
            measured: length(0.72, SurfaceTextureLengthUnitV1::Micrometre),
            uncertainty: length(0.02, SurfaceTextureLengthUnitV1::Micrometre),
            measurement_byte: 0x74,
            calibration_byte: 0x75,
        }
    }

    fn build(self) -> SurfaceTextureObservationV1 {
        observation(
            self.id,
            self.requirement,
            self.metric,
            self.cutoff,
            self.evaluation,
            self.filter_byte,
            self.measured,
            self.uncertainty,
            self.measurement_byte,
            self.calibration_byte,
        )
    }
}

fn singleton(requirement: SurfaceTextureRequirementV1) -> MachineSurfaceTextureDraftV1 {
    MachineSurfaceTextureDraftV1 {
        requirements: vec![requirement],
        observations: Vec::new(),
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn mst_001_types_units_order_and_every_semantic_role_reach_identity() {
    let graph = admitted_graph(0x41);
    let baseline = valid_draft()
        .admit_against(&graph)
        .expect("valid surface-texture state admits");

    let mut reordered_draft = valid_draft();
    reordered_draft.requirements.reverse();
    reordered_draft.observations.reverse();
    let reordered = reordered_draft
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
        ["texture/a", "texture/b"]
    );
    assert_eq!(
        baseline
            .observations()
            .iter()
            .map(|entry| entry.id().canonical_key())
            .collect::<Vec<_>>(),
        ["observation/a", "observation/b"]
    );

    let metre = length(1.0, SurfaceTextureLengthUnitV1::Metre);
    let millimetre = length(1_000.0, SurfaceTextureLengthUnitV1::Millimetre);
    assert_eq!(metre.metres().to_bits(), millimetre.metres().to_bits());
    assert_ne!(metre, millimetre, "submitted unit remains identity-bearing");
    assert_eq!(SurfaceTextureLengthUnitV1::Micrometre.symbol(), "um");
    assert_eq!(SurfaceTextureMetricV1::Ra.symbol(), "Ra");

    let cross_unit = MachineSurfaceTextureDraftV1 {
        requirements: vec![requirement(
            "texture/cross-unit",
            "body/part",
            "surface/part/primary",
            SurfaceTextureMetricV1::Ra,
            SurfaceTextureLimitV1::Maximum(length(0.8, SurfaceTextureLengthUnitV1::Micrometre)),
            length(1.0, SurfaceTextureLengthUnitV1::Metre),
            length(2.0, SurfaceTextureLengthUnitV1::Metre),
            0x31,
            SurfaceLayV1::Parallel,
            Some(frame(OrientationParity::Preserving)),
            SurfaceProductionRuleV1::MaterialRemovalRequired,
            0x61,
            0x62,
            None,
        )],
        observations: vec![observation(
            "observation/cross-unit",
            "texture/cross-unit",
            SurfaceTextureMetricV1::Ra,
            length(1_000.0, SurfaceTextureLengthUnitV1::Millimetre),
            length(2_000.0, SurfaceTextureLengthUnitV1::Millimetre),
            0x31,
            length(0.72, SurfaceTextureLengthUnitV1::Micrometre),
            length(0.02, SurfaceTextureLengthUnitV1::Micrometre),
            0x64,
            0x65,
        )],
    }
    .admit_against(&graph)
    .expect("equivalent coherent-SI evaluation context admits across source units");
    assert_eq!(
        cross_unit.requirements()[0].filter_cutoff().unit(),
        SurfaceTextureLengthUnitV1::Metre
    );
    assert_eq!(
        cross_unit.observations()[0].filter_cutoff().unit(),
        SurfaceTextureLengthUnitV1::Millimetre
    );

    let base_fixture = RequirementFixture::baseline();
    let base = singleton(base_fixture.clone().build())
        .admit_against(&graph)
        .expect("base singleton admits");
    let mut variants = Vec::new();
    let mut changed = base_fixture.clone();
    changed.id = "texture/renamed";
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.declared_body = "body/part-alternate";
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.surface = "surface/part/secondary";
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.metric = SurfaceTextureMetricV1::Rq;
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.limit = SurfaceTextureLimitV1::Range {
        minimum: length(0.4, SurfaceTextureLengthUnitV1::Micrometre),
        maximum: length(0.8, SurfaceTextureLengthUnitV1::Micrometre),
    };
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.cutoff = length(800.0, SurfaceTextureLengthUnitV1::Micrometre);
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.evaluation = length(4_000.0, SurfaceTextureLengthUnitV1::Micrometre);
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.filter_byte = 0x32;
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.lay = SurfaceLayV1::Perpendicular;
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.lay_frame = Some(frame(OrientationParity::Reversing));
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.production_rule = SurfaceProductionRuleV1::MaterialRemovalProhibited;
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.standard_byte = 0x64;
    variants.push(changed);
    let mut changed = base_fixture.clone();
    changed.source_byte = 0x65;
    variants.push(changed);
    let mut changed = base_fixture;
    changed.presentation_byte = None;
    variants.push(changed);
    for variant in variants {
        let admitted = singleton(variant.build())
            .admit_against(&graph)
            .expect("identity mutation remains structurally admissible");
        assert_ne!(base.identity(), admitted.identity());
    }

    let observation_fixture = ObservationFixture::baseline();
    let observed_base = MachineSurfaceTextureDraftV1 {
        requirements: vec![valid_requirement(
            "texture/observed",
            "surface/part/primary",
            SurfaceTextureMetricV1::Ra,
            0x71,
        )],
        observations: vec![observation_fixture.clone().build()],
    }
    .admit_against(&graph)
    .expect("observed baseline admits");
    let mut observation_variants = Vec::new();
    let mut changed = observation_fixture.clone();
    changed.id = "observation/renamed";
    observation_variants.push(changed);
    let mut changed = observation_fixture.clone();
    changed.measured = length(720.0, SurfaceTextureLengthUnitV1::Nanometre);
    observation_variants.push(changed);
    let mut changed = observation_fixture.clone();
    changed.uncertainty = length(0.03, SurfaceTextureLengthUnitV1::Micrometre);
    observation_variants.push(changed);
    let mut changed = observation_fixture.clone();
    changed.measurement_byte = 0x76;
    observation_variants.push(changed);
    let mut changed = observation_fixture;
    changed.calibration_byte = 0x77;
    observation_variants.push(changed);
    for changed_observation in observation_variants {
        let changed = MachineSurfaceTextureDraftV1 {
            requirements: vec![valid_requirement(
                "texture/observed",
                "surface/part/primary",
                SurfaceTextureMetricV1::Ra,
                0x71,
            )],
            observations: vec![changed_observation.build()],
        }
        .admit_against(&graph)
        .expect("changed observation remains context-compatible");
        assert_ne!(observed_base.identity(), changed.identity());
    }

    let changed_graph = singleton(valid_requirement(
        "texture/single",
        "surface/part/primary",
        SurfaceTextureMetricV1::Ra,
        0x61,
    ))
    .admit_against(&admitted_graph(0x42))
    .expect("selectors admit against changed graph");
    assert_ne!(base.identity(), changed_graph.identity());
}

#[test]
#[allow(clippy::too_many_lines)]
fn mst_002_admission_refuses_aliases_invalid_geometry_and_context_mismatch() {
    let graph = admitted_graph(0x41);
    assert_eq!(
        MachineSurfaceTextureDraftV1 {
            requirements: Vec::new(),
            observations: Vec::new(),
        }
        .admit_against(&graph),
        Err(MachineSurfaceTextureAdmissionErrorV1::NoRequirements)
    );

    let duplicate = valid_requirement(
        "texture/duplicate",
        "surface/part/primary",
        SurfaceTextureMetricV1::Ra,
        0x21,
    );
    assert_eq!(
        MachineSurfaceTextureDraftV1 {
            requirements: vec![duplicate.clone(), duplicate],
            observations: Vec::new(),
        }
        .admit_against(&graph),
        Err(
            MachineSurfaceTextureAdmissionErrorV1::DuplicateRequirement {
                requirement: requirement_id("texture/duplicate"),
            }
        )
    );

    let selector_alias = MachineSurfaceTextureDraftV1 {
        requirements: vec![
            valid_requirement(
                "texture/first",
                "surface/part/primary",
                SurfaceTextureMetricV1::Ra,
                0x21,
            ),
            requirement(
                "texture/second",
                "body/part-alternate",
                "surface/part/primary",
                SurfaceTextureMetricV1::Ra,
                SurfaceTextureLimitV1::Maximum(length(0.9, SurfaceTextureLengthUnitV1::Micrometre)),
                length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
                length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
                0x31,
                SurfaceLayV1::Parallel,
                Some(frame(OrientationParity::Preserving)),
                SurfaceProductionRuleV1::MaterialRemovalRequired,
                0x24,
                0x25,
                Some(0x26),
            ),
        ],
        observations: Vec::new(),
    };
    assert_eq!(
        selector_alias.admit_against(&graph),
        Err(
            MachineSurfaceTextureAdmissionErrorV1::DuplicateRequirementSelector {
                first: requirement_id("texture/first"),
                duplicate: requirement_id("texture/second"),
            }
        )
    );

    for (candidate, expected) in [
        (
            requirement(
                "texture/unknown-body",
                "body/missing",
                "surface/part/primary",
                SurfaceTextureMetricV1::Ra,
                SurfaceTextureLimitV1::Maximum(length(0.8, SurfaceTextureLengthUnitV1::Micrometre)),
                length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
                length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
                0x31,
                SurfaceLayV1::Parallel,
                Some(frame(OrientationParity::Preserving)),
                SurfaceProductionRuleV1::Unspecified,
                0x41,
                0x42,
                None,
            ),
            MachineSurfaceTextureAdmissionErrorV1::UnknownBody {
                requirement: requirement_id("texture/unknown-body"),
                body: body("body/missing"),
            },
        ),
        (
            requirement(
                "texture/unknown-patch",
                "body/part",
                "surface/missing",
                SurfaceTextureMetricV1::Ra,
                SurfaceTextureLimitV1::Maximum(length(0.8, SurfaceTextureLengthUnitV1::Micrometre)),
                length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
                length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
                0x31,
                SurfaceLayV1::Parallel,
                Some(frame(OrientationParity::Preserving)),
                SurfaceProductionRuleV1::Unspecified,
                0x41,
                0x42,
                None,
            ),
            MachineSurfaceTextureAdmissionErrorV1::UnknownSurfacePatch {
                requirement: requirement_id("texture/unknown-patch"),
                patch: patch("surface/missing"),
            },
        ),
    ] {
        assert_eq!(singleton(candidate).admit_against(&graph), Err(expected));
    }

    assert!(matches!(
        singleton(requirement(
            "texture/cross-owner",
            "body/part",
            "surface/other/primary",
            SurfaceTextureMetricV1::Ra,
            SurfaceTextureLimitV1::Maximum(length(
                0.8,
                SurfaceTextureLengthUnitV1::Micrometre,
            )),
            length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
            length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
            0x31,
            SurfaceLayV1::Parallel,
            Some(frame(OrientationParity::Preserving)),
            SurfaceProductionRuleV1::Unspecified,
            0x41,
            0x42,
            None,
        ))
        .admit_against(&graph),
        Err(MachineSurfaceTextureAdmissionErrorV1::SurfaceOwnerMismatch {
            requirement,
            body,
            patch,
            body_owner,
            patch_owner,
        }) if requirement == requirement_id("texture/cross-owner")
            && body == crate::body("body/part")
            && patch == crate::patch("surface/other/primary")
            && body_owner == SubsystemId::new("subsystem/part").unwrap()
            && patch_owner == SubsystemId::new("subsystem/other").unwrap()
    ));

    let invalid_cases = [
        (
            SurfaceTextureLimitV1::Maximum(length(0.0, SurfaceTextureLengthUnitV1::Micrometre)),
            length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
            length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
            SurfaceLayV1::Parallel,
            Some(frame(OrientationParity::Preserving)),
            SurfaceTextureRequirementIssueV1::NonPositiveLimit,
        ),
        (
            SurfaceTextureLimitV1::Range {
                minimum: length(0.8, SurfaceTextureLengthUnitV1::Micrometre),
                maximum: length(0.8, SurfaceTextureLengthUnitV1::Micrometre),
            },
            length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
            length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
            SurfaceLayV1::Parallel,
            Some(frame(OrientationParity::Preserving)),
            SurfaceTextureRequirementIssueV1::RangeNotIncreasing,
        ),
        (
            SurfaceTextureLimitV1::Maximum(length(0.8, SurfaceTextureLengthUnitV1::Micrometre)),
            length(0.0, SurfaceTextureLengthUnitV1::Millimetre),
            length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
            SurfaceLayV1::Parallel,
            Some(frame(OrientationParity::Preserving)),
            SurfaceTextureRequirementIssueV1::NonPositiveFilterCutoff,
        ),
        (
            SurfaceTextureLimitV1::Maximum(length(0.8, SurfaceTextureLengthUnitV1::Micrometre)),
            length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
            length(0.0, SurfaceTextureLengthUnitV1::Millimetre),
            SurfaceLayV1::Parallel,
            Some(frame(OrientationParity::Preserving)),
            SurfaceTextureRequirementIssueV1::NonPositiveEvaluationLength,
        ),
        (
            SurfaceTextureLimitV1::Maximum(length(0.8, SurfaceTextureLengthUnitV1::Micrometre)),
            length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
            length(0.4, SurfaceTextureLengthUnitV1::Millimetre),
            SurfaceLayV1::Parallel,
            Some(frame(OrientationParity::Preserving)),
            SurfaceTextureRequirementIssueV1::EvaluationShorterThanCutoff,
        ),
        (
            SurfaceTextureLimitV1::Maximum(length(0.8, SurfaceTextureLengthUnitV1::Micrometre)),
            length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
            length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
            SurfaceLayV1::Parallel,
            None,
            SurfaceTextureRequirementIssueV1::MissingLayFrame,
        ),
        (
            SurfaceTextureLimitV1::Maximum(length(0.8, SurfaceTextureLengthUnitV1::Micrometre)),
            length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
            length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
            SurfaceLayV1::Unspecified,
            Some(frame(OrientationParity::Preserving)),
            SurfaceTextureRequirementIssueV1::UnexpectedLayFrame,
        ),
        (
            SurfaceTextureLimitV1::Maximum(length(0.8, SurfaceTextureLengthUnitV1::Micrometre)),
            length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
            length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
            SurfaceLayV1::Particulate,
            Some(frame(OrientationParity::Preserving)),
            SurfaceTextureRequirementIssueV1::UnexpectedLayFrame,
        ),
    ];
    for (index, (limit, cutoff, evaluation, lay, lay_frame, issue)) in
        invalid_cases.into_iter().enumerate()
    {
        let id = format!("texture/invalid-{index}");
        let candidate = requirement(
            &id,
            "body/part",
            "surface/part/primary",
            SurfaceTextureMetricV1::Ra,
            limit,
            cutoff,
            evaluation,
            0x31,
            lay,
            lay_frame,
            SurfaceProductionRuleV1::Unspecified,
            0x41,
            0x42,
            None,
        );
        assert_eq!(
            singleton(candidate).admit_against(&graph),
            Err(MachineSurfaceTextureAdmissionErrorV1::InvalidRequirement {
                requirement: requirement_id(&id),
                issue,
            })
        );
    }

    let particulate = requirement(
        "texture/particulate",
        "body/part",
        "surface/part/primary",
        SurfaceTextureMetricV1::Ra,
        SurfaceTextureLimitV1::Maximum(length(0.8, SurfaceTextureLengthUnitV1::Micrometre)),
        length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
        length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
        0x31,
        SurfaceLayV1::Particulate,
        None,
        SurfaceProductionRuleV1::Unspecified,
        0x41,
        0x42,
        None,
    );
    singleton(particulate)
        .admit_against(&graph)
        .expect("non-directional particulate lay admits without a frame");
    assert!(!SurfaceLayV1::Particulate.requires_frame());

    let base_requirement = valid_requirement(
        "texture/context",
        "surface/part/primary",
        SurfaceTextureMetricV1::Ra,
        0x31,
    );
    let context_cases = [
        (
            SurfaceTextureMetricV1::Rq,
            length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
            length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
            0x31,
            SurfaceTextureObservationContextIssueV1::Metric,
        ),
        (
            SurfaceTextureMetricV1::Ra,
            length(0.9, SurfaceTextureLengthUnitV1::Millimetre),
            length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
            0x31,
            SurfaceTextureObservationContextIssueV1::FilterCutoff,
        ),
        (
            SurfaceTextureMetricV1::Ra,
            length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
            length(4.8, SurfaceTextureLengthUnitV1::Millimetre),
            0x31,
            SurfaceTextureObservationContextIssueV1::EvaluationLength,
        ),
        (
            SurfaceTextureMetricV1::Ra,
            length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
            length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
            0x32,
            SurfaceTextureObservationContextIssueV1::FilterSpecification,
        ),
    ];
    for (index, (metric, cutoff, evaluation, filter_byte, issue)) in
        context_cases.into_iter().enumerate()
    {
        let observation_key = format!("observation/mismatch-{index}");
        let draft = MachineSurfaceTextureDraftV1 {
            requirements: vec![base_requirement.clone()],
            observations: vec![observation(
                &observation_key,
                "texture/context",
                metric,
                cutoff,
                evaluation,
                filter_byte,
                length(0.7, SurfaceTextureLengthUnitV1::Micrometre),
                length(0.02, SurfaceTextureLengthUnitV1::Micrometre),
                0x51,
                0x52,
            )],
        };
        assert_eq!(
            draft.admit_against(&graph),
            Err(
                MachineSurfaceTextureAdmissionErrorV1::ObservationContextMismatch {
                    observation: observation_id(&observation_key),
                    requirement: requirement_id("texture/context"),
                    issue,
                }
            )
        );
    }

    let duplicate_observation = valid_observation(
        "observation/duplicate",
        "texture/context",
        SurfaceTextureMetricV1::Ra,
        0x51,
    );
    assert_eq!(
        MachineSurfaceTextureDraftV1 {
            requirements: vec![base_requirement.clone()],
            observations: vec![duplicate_observation.clone(), duplicate_observation],
        }
        .admit_against(&graph),
        Err(
            MachineSurfaceTextureAdmissionErrorV1::DuplicateObservation {
                observation: observation_id("observation/duplicate"),
            }
        )
    );
    assert_eq!(
        MachineSurfaceTextureDraftV1 {
            requirements: vec![base_requirement.clone()],
            observations: vec![valid_observation(
                "observation/missing",
                "texture/missing",
                SurfaceTextureMetricV1::Ra,
                0x51,
            )],
        }
        .admit_against(&graph),
        Err(
            MachineSurfaceTextureAdmissionErrorV1::MissingObservationRequirement {
                observation: observation_id("observation/missing"),
                requirement: requirement_id("texture/missing"),
            }
        )
    );

    let visibly_over_limit = MachineSurfaceTextureDraftV1 {
        requirements: vec![base_requirement],
        observations: vec![observation(
            "observation/over-limit",
            "texture/context",
            SurfaceTextureMetricV1::Ra,
            length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
            length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
            0x31,
            length(8.0, SurfaceTextureLengthUnitV1::Micrometre),
            length(0.02, SurfaceTextureLengthUnitV1::Micrometre),
            0x51,
            0x52,
        )],
    }
    .admit_against(&graph)
    .expect("an observation is retained without automatic pass/fail");
    assert_eq!(visibly_over_limit.observations().len(), 1);

    assert_eq!(
        SurfaceTextureLengthV1::try_new(f64::NAN, SurfaceTextureLengthUnitV1::Metre),
        Err(SurfaceTextureLengthErrorV1::NonFinite)
    );
    assert_eq!(
        SurfaceTextureLengthV1::try_new(-1.0, SurfaceTextureLengthUnitV1::Metre),
        Err(SurfaceTextureLengthErrorV1::Negative)
    );
    assert_eq!(
        SurfaceTextureLengthV1::try_new(f64::from_bits(1), SurfaceTextureLengthUnitV1::Nanometre,),
        Err(SurfaceTextureLengthErrorV1::SiUnderflow)
    );
    assert_eq!(
        SurfaceTextureLengthV1::try_new(-0.0, SurfaceTextureLengthUnitV1::Metre),
        SurfaceTextureLengthV1::try_new(0.0, SurfaceTextureLengthUnitV1::Metre)
    );
    assert_eq!(
        SurfaceTextureRequirementIssueV1::MissingLayFrame.code(),
        "SurfaceTextureMissingLayFrame"
    );
    assert_eq!(
        MachineSurfaceTextureAdmissionErrorV1::NoRequirements.code(),
        "MachineSurfaceTextureNoRequirements"
    );
}

fn boundary_graph() -> AdmittedMachineGraph {
    let target = body("body/boundary");
    MachineGraphDraft {
        clocks: Vec::new(),
        subsystems: vec![SubsystemSpec {
            id: SubsystemId::new("subsystem/boundary").expect("canonical subsystem"),
            model: ModelRef::new("models/surface-texture-boundary", nz(1), [0x61; 32])
                .expect("canonical model"),
            bodies: vec![target.clone()],
            surface_patches: (0..(MAX_MACHINE_SURFACE_TEXTURE_REQUIREMENTS_V1 / 4))
                .map(|index| patch(&format!("surface/boundary/p{index:04}")))
                .collect(),
            contact_features: Vec::new(),
            state_slots: Vec::new(),
        }],
        terminals: Vec::new(),
        ports: Vec::new(),
        relations: Vec::new(),
        materials: vec![material(target, "materials/boundary", 0x62)],
        interfaces: Vec::new(),
    }
    .admit()
    .expect("exact-cap surface-texture graph admits")
}

fn metric_for(index: usize) -> SurfaceTextureMetricV1 {
    match index % 4 {
        0 => SurfaceTextureMetricV1::Ra,
        1 => SurfaceTextureMetricV1::Rq,
        2 => SurfaceTextureMetricV1::Rz,
        _ => SurfaceTextureMetricV1::Rt,
    }
}

fn boundary_draft() -> MachineSurfaceTextureDraftV1 {
    let requirements = (0..MAX_MACHINE_SURFACE_TEXTURE_REQUIREMENTS_V1)
        .map(|index| {
            requirement(
                &format!("texture/boundary-r{index:04}"),
                "body/boundary",
                &format!("surface/boundary/p{:04}", index / 4),
                metric_for(index),
                SurfaceTextureLimitV1::Maximum(length(0.8, SurfaceTextureLengthUnitV1::Micrometre)),
                length(0.8, SurfaceTextureLengthUnitV1::Millimetre),
                length(4.0, SurfaceTextureLengthUnitV1::Millimetre),
                0x31,
                SurfaceLayV1::Parallel,
                Some(frame(OrientationParity::Preserving)),
                SurfaceProductionRuleV1::MaterialRemovalRequired,
                0x71,
                0x72,
                Some(0x73),
            )
        })
        .collect::<Vec<_>>();
    let observations = (0..MAX_MACHINE_SURFACE_TEXTURE_OBSERVATIONS_V1)
        .map(|index| {
            valid_observation(
                &format!("observation/boundary-o{index:04}"),
                &format!("texture/boundary-r{index:04}"),
                metric_for(index),
                0x74,
            )
        })
        .collect();
    MachineSurfaceTextureDraftV1 {
        requirements,
        observations,
    }
}

#[test]
fn mst_003_exact_resource_caps_admit_and_one_over_refuses_before_deduplication() {
    let graph = boundary_graph();
    let exact = boundary_draft();
    let admitted = exact
        .clone()
        .admit_against(&graph)
        .expect("simultaneous exact requirement and observation caps admit");
    assert_eq!(
        admitted.requirements().len(),
        MAX_MACHINE_SURFACE_TEXTURE_REQUIREMENTS_V1
    );
    assert_eq!(
        admitted.observations().len(),
        MAX_MACHINE_SURFACE_TEXTURE_OBSERVATIONS_V1
    );

    let mut too_many_requirements = exact.clone();
    let repeated_requirement = too_many_requirements.requirements[0].clone();
    too_many_requirements
        .requirements
        .push(repeated_requirement);
    assert_eq!(
        too_many_requirements.admit_against(&graph),
        Err(MachineSurfaceTextureAdmissionErrorV1::RequirementLimit {
            actual: MAX_MACHINE_SURFACE_TEXTURE_REQUIREMENTS_V1 + 1,
            max: MAX_MACHINE_SURFACE_TEXTURE_REQUIREMENTS_V1,
        })
    );

    let mut too_many_observations = exact;
    let repeated_observation = too_many_observations.observations[0].clone();
    too_many_observations
        .observations
        .push(repeated_observation);
    assert_eq!(
        too_many_observations.admit_against(&graph),
        Err(MachineSurfaceTextureAdmissionErrorV1::ObservationLimit {
            actual: MAX_MACHINE_SURFACE_TEXTURE_OBSERVATIONS_V1 + 1,
            max: MAX_MACHINE_SURFACE_TEXTURE_OBSERVATIONS_V1,
        })
    );
}

#[test]
fn mst_004_identical_input_replays_the_complete_receipt() {
    let graph = admitted_graph(0x41);
    let first = valid_draft()
        .admit_against(&graph)
        .expect("first replay admits");
    let second = valid_draft()
        .admit_against(&graph)
        .expect("second replay admits");
    assert_eq!(first.identity(), second.identity());
    assert_eq!(first.identity_receipt(), second.identity_receipt());
    assert_eq!(
        first.identity_receipt().canonical_preimage(),
        second.identity_receipt().canonical_preimage()
    );
}
