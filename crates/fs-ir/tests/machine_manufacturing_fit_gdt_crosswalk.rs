//! Machine-IR fit/GD&T applicability crosswalk (Gauntlet G0/G3/G5).

use core::num::NonZeroU64;

use fs_blake3::ContentHash;
use fs_ir::machine::manufacturing::ManufacturingArtifactRefV1;
use fs_ir::machine::manufacturing::datum_system::{
    AdmittedMachineDatumSystemV1, DatumFeatureBindingV1, DatumFeatureIdV1, DatumFeatureTargetV1,
    DatumReferenceFrameIdV1, DatumReferenceFrameV1, MachineDatumSystemDraftV1,
};
use fs_ir::machine::manufacturing::fit_clearance::{
    AdmittedMachineFitClearanceV1, FitAllowanceV1, FitFeatureSelectorV1, FitLengthUnitV1,
    FitPairTargetV1, FitPresentationRefV1, FitRegimeV1, FitRequirementIdV1, FitRequirementV1,
    FitSemanticSourceRefV1, FitSpecificationRefV1, MachineFitClearanceDraftV1, PositiveFitLengthV1,
    SignedFitLengthV1,
};
use fs_ir::machine::manufacturing::fit_gdt_crosswalk::{
    FitGdtEndpointLinkV1, FitGdtEndpointRoleV1, MACHINE_FIT_GDT_CROSSWALK_SCHEMA_VERSION_V1,
    MAX_MACHINE_FIT_GDT_LINKS_V1, MachineFitGdtCrosswalkAdmissionErrorV1,
    MachineFitGdtCrosswalkDraftV1,
};
use fs_ir::machine::manufacturing::geometric_tolerance::{
    AdmittedMachineGeometricToleranceV1, GeometricCharacteristicV1, GeometricToleranceControlIdV1,
    GeometricToleranceControlV1, GeometricToleranceLengthUnitV1, GeometricToleranceLengthV1,
    GeometricTolerancePresentationRefV1, GeometricToleranceSemanticSourceRefV1,
    GeometricToleranceSpecificationRefV1, MachineGeometricToleranceDraftV1,
};
use fs_ir::machine::{
    AdmittedMachineGraph, BodyId, ContactFeatureId, MachineGraphDraft, MaterialBinding,
    MaterialCardRef, MaterialTarget, ModelRef, SubsystemId, SubsystemSpec, SurfacePatchId,
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

fn patch(key: &str) -> SurfacePatchId {
    SurfacePatchId::new(key).expect("fixture surface-patch key is canonical")
}

fn fit_id(key: &str) -> FitRequirementIdV1 {
    FitRequirementIdV1::new(key).expect("fixture fit key is canonical")
}

fn control_id(key: &str) -> GeometricToleranceControlIdV1 {
    GeometricToleranceControlIdV1::new(key).expect("fixture control key is canonical")
}

fn datum_id(key: &str) -> DatumFeatureIdV1 {
    DatumFeatureIdV1::new(key).expect("fixture datum key is canonical")
}

fn frame_id(key: &str) -> DatumReferenceFrameIdV1 {
    DatumReferenceFrameIdV1::new(key).expect("fixture frame key is canonical")
}

fn material(target: BodyId, key: &str, byte: u8) -> MaterialBinding {
    MaterialBinding {
        target: MaterialTarget::Body(target),
        material: MaterialCardRef::new(key, nz(1), [byte; 32])
            .expect("fixture material is canonical"),
    }
}

fn artifact(namespace: &str, byte: u8) -> ManufacturingArtifactRefV1 {
    ManufacturingArtifactRefV1::new(namespace, nz(1), ContentHash([byte; 32]))
        .expect("fixture artifact is canonical")
}

fn graph(model_byte: u8) -> AdmittedMachineGraph {
    let internal = body("body/internal");
    let external = body("body/external");
    MachineGraphDraft {
        clocks: Vec::new(),
        subsystems: vec![SubsystemSpec {
            id: SubsystemId::new("subsystem/fit-gdt").expect("canonical subsystem"),
            model: ModelRef::new("models/fit-gdt", nz(1), [model_byte; 32])
                .expect("canonical model"),
            bodies: vec![internal.clone(), external.clone()],
            surface_patches: vec![
                patch("surface/internal/datum"),
                patch("surface/external/datum"),
                patch("surface/internal/control"),
                patch("surface/internal/control-alternate"),
                patch("surface/external/control"),
                patch("surface/external/control-alternate"),
            ],
            contact_features: vec![feature("contact/internal"), feature("contact/external")],
            state_slots: Vec::new(),
        }],
        terminals: Vec::new(),
        ports: Vec::new(),
        relations: Vec::new(),
        materials: vec![
            material(internal, "materials/internal", 0x21),
            material(external, "materials/external", 0x22),
        ],
        interfaces: Vec::new(),
    }
    .admit()
    .expect("fit/GD&T fixture graph admits")
}

fn datum_feature(id: &str, declared_body: &str, surface: &str) -> DatumFeatureBindingV1 {
    DatumFeatureBindingV1::new(
        datum_id(id),
        body(declared_body),
        DatumFeatureTargetV1::SurfacePatch(patch(surface)),
    )
}

fn datum_catalog(graph: &AdmittedMachineGraph) -> AdmittedMachineDatumSystemV1 {
    MachineDatumSystemDraftV1 {
        datum_features: vec![
            datum_feature("datum/internal", "body/internal", "surface/internal/datum"),
            datum_feature("datum/external", "body/external", "surface/external/datum"),
        ],
        reference_frames: vec![
            DatumReferenceFrameV1::new(
                frame_id("datum-frame/internal"),
                datum_id("datum/internal"),
                None,
                None,
            ),
            DatumReferenceFrameV1::new(
                frame_id("datum-frame/external"),
                datum_id("datum/external"),
                None,
                None,
            ),
        ],
    }
    .admit_against(graph)
    .expect("fit/GD&T datum catalog admits")
}

fn signed_fit(value: f64) -> SignedFitLengthV1 {
    SignedFitLengthV1::try_new(value, FitLengthUnitV1::Micrometre)
        .expect("fixture signed fit length admits")
}

fn fit_requirement(basic_size_mm: f64) -> FitRequirementV1 {
    FitRequirementV1::new(
        fit_id("fit/main"),
        FitPairTargetV1::new(
            FitFeatureSelectorV1::new(body("body/internal"), feature("contact/internal")),
            FitFeatureSelectorV1::new(body("body/external"), feature("contact/external")),
        ),
        PositiveFitLengthV1::try_new(basic_size_mm, FitLengthUnitV1::Millimetre)
            .expect("fixture basic size admits"),
        FitAllowanceV1::try_new(signed_fit(10.0), signed_fit(30.0))
            .expect("fixture clearance allowance admits"),
        FitSpecificationRefV1::new(artifact("fit/specification", 0x31)),
        FitSemanticSourceRefV1::new(artifact("fit/semantic-source", 0x32)),
        Some(FitPresentationRefV1::new(artifact(
            "fit/presentation",
            0x33,
        ))),
    )
}

fn fit_catalog(graph: &AdmittedMachineGraph, basic_size_mm: f64) -> AdmittedMachineFitClearanceV1 {
    MachineFitClearanceDraftV1 {
        requirements: vec![fit_requirement(basic_size_mm)],
    }
    .admit_against(graph)
    .expect("fit catalog admits")
}

fn geometric_length(value_um: f64) -> GeometricToleranceLengthV1 {
    GeometricToleranceLengthV1::try_new(value_um, GeometricToleranceLengthUnitV1::Micrometre)
        .expect("fixture geometric length admits")
}

fn control(
    id: &str,
    declared_body: &str,
    controlled_patch: &str,
    zone_width_um: f64,
    byte: u8,
) -> GeometricToleranceControlV1 {
    GeometricToleranceControlV1::new(
        control_id(id),
        body(declared_body),
        patch(controlled_patch),
        GeometricCharacteristicV1::Flatness,
        geometric_length(zone_width_um),
        None,
        GeometricToleranceSpecificationRefV1::new(artifact("gdt/specification", byte)),
        GeometricToleranceSemanticSourceRefV1::new(artifact(
            "gdt/semantic-source",
            byte.wrapping_add(1),
        )),
        Some(GeometricTolerancePresentationRefV1::new(artifact(
            "gdt/presentation",
            byte.wrapping_add(2),
        ))),
    )
}

fn geometric_catalog(
    graph: &AdmittedMachineGraph,
    datum: &AdmittedMachineDatumSystemV1,
    zone_width_um: f64,
) -> AdmittedMachineGeometricToleranceV1 {
    MachineGeometricToleranceDraftV1 {
        controls: vec![
            control(
                "gdt/internal",
                "body/internal",
                "surface/internal/control",
                zone_width_um,
                0x41,
            ),
            control(
                "gdt/internal-alternate",
                "body/internal",
                "surface/internal/control-alternate",
                zone_width_um,
                0x44,
            ),
            control(
                "gdt/external",
                "body/external",
                "surface/external/control",
                zone_width_um,
                0x47,
            ),
            control(
                "gdt/external-alternate",
                "body/external",
                "surface/external/control-alternate",
                zone_width_um,
                0x4a,
            ),
        ],
    }
    .admit_against(graph, datum)
    .expect("geometric-tolerance catalog admits")
}

struct Catalogs {
    graph: AdmittedMachineGraph,
    fit: AdmittedMachineFitClearanceV1,
    geometric: AdmittedMachineGeometricToleranceV1,
}

fn catalogs(model_byte: u8, basic_size_mm: f64, zone_width_um: f64) -> Catalogs {
    let graph = graph(model_byte);
    let datum = datum_catalog(&graph);
    let fit = fit_catalog(&graph, basic_size_mm);
    let geometric = geometric_catalog(&graph, &datum, zone_width_um);
    Catalogs {
        graph,
        fit,
        geometric,
    }
}

fn link(requirement: &str, role: FitGdtEndpointRoleV1, control: &str) -> FitGdtEndpointLinkV1 {
    FitGdtEndpointLinkV1::new(fit_id(requirement), role, control_id(control))
}

fn valid_links() -> Vec<FitGdtEndpointLinkV1> {
    vec![
        link("fit/main", FitGdtEndpointRoleV1::External, "gdt/external"),
        link("fit/main", FitGdtEndpointRoleV1::Internal, "gdt/internal"),
    ]
}

#[test]
fn mfgc_001_role_complete_crosswalk_is_graph_bound_resolved_and_replayable() {
    let catalogs = catalogs(0x11, 10.0, 5.0);
    let baseline = MachineFitGdtCrosswalkDraftV1 {
        links: valid_links(),
    }
    .admit_against(&catalogs.fit, &catalogs.geometric)
    .expect("valid fit/GD&T crosswalk admits");
    let reordered = MachineFitGdtCrosswalkDraftV1 {
        links: valid_links().into_iter().rev().collect(),
    }
    .admit_against(&catalogs.fit, &catalogs.geometric)
    .expect("caller link order is non-semantic");

    assert_eq!(baseline.graph(), catalogs.graph.identity());
    assert_eq!(baseline.fit_catalog(), catalogs.fit.identity());
    assert_eq!(
        baseline.fit_catalog_receipt(),
        catalogs.fit.identity_receipt()
    );
    assert_eq!(
        baseline.geometric_tolerance_catalog(),
        catalogs.geometric.identity()
    );
    assert_eq!(
        baseline.geometric_tolerance_catalog_receipt(),
        catalogs.geometric.identity_receipt()
    );
    assert_eq!(baseline.identity(), reordered.identity());
    assert_eq!(baseline.identity_receipt(), reordered.identity_receipt());
    assert_eq!(baseline.endpoints().len(), 2);

    let internal = &baseline.endpoints()[0];
    assert_eq!(internal.link().fit_requirement(), &fit_id("fit/main"));
    assert_eq!(internal.link().role(), FitGdtEndpointRoleV1::Internal);
    assert_eq!(
        internal.link().geometric_control(),
        &control_id("gdt/internal")
    );
    assert_eq!(internal.declared_body(), &body("body/internal"));
    assert_eq!(internal.fit_feature(), &feature("contact/internal"));
    assert_eq!(
        internal.controlled_patch(),
        &patch("surface/internal/control")
    );
    assert_eq!(
        internal.characteristic(),
        GeometricCharacteristicV1::Flatness
    );
    assert_eq!(
        internal.zone_width().submitted_value().to_bits(),
        5.0_f64.to_bits()
    );
    assert_eq!(internal.datum_frame(), None);

    let external = &baseline.endpoints()[1];
    assert_eq!(external.link().role(), FitGdtEndpointRoleV1::External);
    assert_eq!(external.declared_body(), &body("body/external"));
    assert_eq!(external.fit_feature(), &feature("contact/external"));
    assert_eq!(
        external.controlled_patch(),
        &patch("surface/external/control")
    );
    assert_eq!(
        catalogs.fit.requirements()[0].allowance().regime(),
        FitRegimeV1::Clearance
    );
    assert_eq!(
        baseline.identity_receipt().canonical_preimage(),
        reordered.identity_receipt().canonical_preimage()
    );
}

#[test]
fn mfgc_002_catalog_content_and_applicability_links_are_identity_semantic() {
    let baseline_catalogs = catalogs(0x11, 10.0, 5.0);
    let baseline = MachineFitGdtCrosswalkDraftV1 {
        links: valid_links(),
    }
    .admit_against(&baseline_catalogs.fit, &baseline_catalogs.geometric)
    .expect("baseline crosswalk admits");

    let fit_changed = catalogs(0x11, 11.0, 5.0);
    let fit_changed_receipt = MachineFitGdtCrosswalkDraftV1 {
        links: valid_links(),
    }
    .admit_against(&fit_changed.fit, &fit_changed.geometric)
    .expect("fit-mutated crosswalk admits");
    assert_ne!(baseline.fit_catalog(), fit_changed_receipt.fit_catalog());
    assert_ne!(baseline.identity(), fit_changed_receipt.identity());

    let gdt_changed = catalogs(0x11, 10.0, 6.0);
    let gdt_changed_receipt = MachineFitGdtCrosswalkDraftV1 {
        links: valid_links(),
    }
    .admit_against(&gdt_changed.fit, &gdt_changed.geometric)
    .expect("GD&T-mutated crosswalk admits");
    assert_ne!(
        baseline.geometric_tolerance_catalog(),
        gdt_changed_receipt.geometric_tolerance_catalog()
    );
    assert_ne!(baseline.identity(), gdt_changed_receipt.identity());

    let alternate_link = MachineFitGdtCrosswalkDraftV1 {
        links: vec![
            link(
                "fit/main",
                FitGdtEndpointRoleV1::Internal,
                "gdt/internal-alternate",
            ),
            link("fit/main", FitGdtEndpointRoleV1::External, "gdt/external"),
        ],
    }
    .admit_against(&baseline_catalogs.fit, &baseline_catalogs.geometric)
    .expect("alternate applicability link admits");
    assert_ne!(baseline.identity(), alternate_link.identity());
    assert_ne!(
        baseline.endpoints()[0].controlled_patch(),
        alternate_link.endpoints()[0].controlled_patch()
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One refusal matrix preserves exact diagnostic fields.
fn mfgc_003_admission_refuses_graph_ids_aliases_unknowns_body_mismatch_and_gaps() {
    let catalogs = catalogs(0x11, 10.0, 5.0);
    assert_eq!(
        MachineFitGdtCrosswalkDraftV1 { links: Vec::new() }
            .admit_against(&catalogs.fit, &catalogs.geometric),
        Err(MachineFitGdtCrosswalkAdmissionErrorV1::NoLinks)
    );
    assert_eq!(
        MachineFitGdtCrosswalkDraftV1 {
            links: vec![valid_links()[0].clone(); MAX_MACHINE_FIT_GDT_LINKS_V1 + 1],
        }
        .admit_against(&catalogs.fit, &catalogs.geometric),
        Err(MachineFitGdtCrosswalkAdmissionErrorV1::LinkLimit {
            actual: MAX_MACHINE_FIT_GDT_LINKS_V1 + 1,
            max: MAX_MACHINE_FIT_GDT_LINKS_V1,
        })
    );

    let other_graph = graph(0x12);
    let other_datum = datum_catalog(&other_graph);
    let other_geometric = geometric_catalog(&other_graph, &other_datum, 5.0);
    assert!(matches!(
        MachineFitGdtCrosswalkDraftV1 {
            links: valid_links(),
        }
        .admit_against(&catalogs.fit, &other_geometric),
        Err(MachineFitGdtCrosswalkAdmissionErrorV1::CatalogGraphMismatch {
            fit_graph,
            geometric_tolerance_graph,
        }) if fit_graph == catalogs.fit.graph()
            && geometric_tolerance_graph == other_geometric.graph()
    ));

    assert_eq!(
        MachineFitGdtCrosswalkDraftV1 {
            links: vec![
                link("fit/main", FitGdtEndpointRoleV1::Internal, "gdt/internal",),
                link(
                    "fit/main",
                    FitGdtEndpointRoleV1::Internal,
                    "gdt/internal-alternate",
                ),
                link("fit/main", FitGdtEndpointRoleV1::External, "gdt/external",),
            ],
        }
        .admit_against(&catalogs.fit, &catalogs.geometric),
        Err(MachineFitGdtCrosswalkAdmissionErrorV1::DuplicateEndpoint {
            requirement: fit_id("fit/main"),
            role: FitGdtEndpointRoleV1::Internal,
            first_control: control_id("gdt/internal"),
            duplicate_control: control_id("gdt/internal-alternate"),
        })
    );
    assert_eq!(
        MachineFitGdtCrosswalkDraftV1 {
            links: vec![link(
                "fit/missing",
                FitGdtEndpointRoleV1::Internal,
                "gdt/internal",
            )],
        }
        .admit_against(&catalogs.fit, &catalogs.geometric),
        Err(
            MachineFitGdtCrosswalkAdmissionErrorV1::UnknownFitRequirement {
                requirement: fit_id("fit/missing"),
                role: FitGdtEndpointRoleV1::Internal,
            }
        )
    );
    assert_eq!(
        MachineFitGdtCrosswalkDraftV1 {
            links: vec![link(
                "fit/main",
                FitGdtEndpointRoleV1::Internal,
                "gdt/missing",
            )],
        }
        .admit_against(&catalogs.fit, &catalogs.geometric),
        Err(
            MachineFitGdtCrosswalkAdmissionErrorV1::UnknownGeometricControl {
                requirement: fit_id("fit/main"),
                role: FitGdtEndpointRoleV1::Internal,
                control: control_id("gdt/missing"),
            }
        )
    );
    assert_eq!(
        MachineFitGdtCrosswalkDraftV1 {
            links: vec![
                link("fit/main", FitGdtEndpointRoleV1::Internal, "gdt/external",),
                link("fit/main", FitGdtEndpointRoleV1::External, "gdt/external",),
            ],
        }
        .admit_against(&catalogs.fit, &catalogs.geometric),
        Err(
            MachineFitGdtCrosswalkAdmissionErrorV1::DeclaredBodyMismatch {
                requirement: fit_id("fit/main"),
                role: FitGdtEndpointRoleV1::Internal,
                control: control_id("gdt/external"),
                fit_body: body("body/internal"),
                geometric_tolerance_body: body("body/external"),
            }
        )
    );
    assert_eq!(
        MachineFitGdtCrosswalkDraftV1 {
            links: vec![link(
                "fit/main",
                FitGdtEndpointRoleV1::Internal,
                "gdt/internal",
            )],
        }
        .admit_against(&catalogs.fit, &catalogs.geometric),
        Err(MachineFitGdtCrosswalkAdmissionErrorV1::MissingEndpoint {
            requirement: fit_id("fit/main"),
            role: FitGdtEndpointRoleV1::External,
        })
    );
    assert_eq!(
        MachineFitGdtCrosswalkAdmissionErrorV1::NoLinks.code(),
        "MachineFitGdtNoLinks"
    );
}

fn boundary_graph() -> AdmittedMachineGraph {
    let internal = body("body/boundary-internal");
    let external = body("body/boundary-external");
    let mut features = vec![feature("contact/boundary/internal")];
    features.extend(
        (0..MAX_MACHINE_FIT_GDT_LINKS_V1 / 2)
            .map(|index| feature(&format!("contact/boundary/external-{index:04}"))),
    );
    MachineGraphDraft {
        clocks: Vec::new(),
        subsystems: vec![SubsystemSpec {
            id: SubsystemId::new("subsystem/boundary").expect("canonical subsystem"),
            model: ModelRef::new("models/fit-gdt-boundary", nz(1), [0x71; 32])
                .expect("canonical model"),
            bodies: vec![internal.clone(), external.clone()],
            surface_patches: vec![
                patch("surface/boundary/internal-datum"),
                patch("surface/boundary/external-datum"),
                patch("surface/boundary/internal-control"),
                patch("surface/boundary/external-control"),
            ],
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
    .expect("fit/GD&T boundary graph admits")
}

fn boundary_datum(graph: &AdmittedMachineGraph) -> AdmittedMachineDatumSystemV1 {
    MachineDatumSystemDraftV1 {
        datum_features: vec![
            datum_feature(
                "datum/boundary-internal",
                "body/boundary-internal",
                "surface/boundary/internal-datum",
            ),
            datum_feature(
                "datum/boundary-external",
                "body/boundary-external",
                "surface/boundary/external-datum",
            ),
        ],
        reference_frames: vec![
            DatumReferenceFrameV1::new(
                frame_id("datum-frame/boundary-internal"),
                datum_id("datum/boundary-internal"),
                None,
                None,
            ),
            DatumReferenceFrameV1::new(
                frame_id("datum-frame/boundary-external"),
                datum_id("datum/boundary-external"),
                None,
                None,
            ),
        ],
    }
    .admit_against(graph)
    .expect("boundary datum catalog admits")
}

fn boundary_fit(graph: &AdmittedMachineGraph) -> AdmittedMachineFitClearanceV1 {
    MachineFitClearanceDraftV1 {
        requirements: (0..MAX_MACHINE_FIT_GDT_LINKS_V1 / 2)
            .map(|index| {
                FitRequirementV1::new(
                    fit_id(&format!("fit/boundary-{index:04}")),
                    FitPairTargetV1::new(
                        FitFeatureSelectorV1::new(
                            body("body/boundary-internal"),
                            feature("contact/boundary/internal"),
                        ),
                        FitFeatureSelectorV1::new(
                            body("body/boundary-external"),
                            feature(&format!("contact/boundary/external-{index:04}")),
                        ),
                    ),
                    PositiveFitLengthV1::try_new(10.0, FitLengthUnitV1::Millimetre)
                        .expect("boundary size admits"),
                    FitAllowanceV1::try_new(signed_fit(1.0), signed_fit(2.0))
                        .expect("boundary allowance admits"),
                    FitSpecificationRefV1::new(artifact("fit/boundary-specification", 0x74)),
                    FitSemanticSourceRefV1::new(artifact("fit/boundary-source", 0x75)),
                    None,
                )
            })
            .collect(),
    }
    .admit_against(graph)
    .expect("boundary fit catalog admits")
}

fn boundary_geometric(
    graph: &AdmittedMachineGraph,
    datum: &AdmittedMachineDatumSystemV1,
) -> AdmittedMachineGeometricToleranceV1 {
    MachineGeometricToleranceDraftV1 {
        controls: vec![
            control(
                "gdt/boundary-internal",
                "body/boundary-internal",
                "surface/boundary/internal-control",
                5.0,
                0x76,
            ),
            control(
                "gdt/boundary-external",
                "body/boundary-external",
                "surface/boundary/external-control",
                5.0,
                0x79,
            ),
        ],
    }
    .admit_against(graph, datum)
    .expect("boundary geometric-tolerance catalog admits")
}

#[test]
fn mfgc_004_exact_endpoint_link_cap_admits() {
    let graph = boundary_graph();
    let datum = boundary_datum(&graph);
    let fit = boundary_fit(&graph);
    let geometric = boundary_geometric(&graph, &datum);
    let links = (0..MAX_MACHINE_FIT_GDT_LINKS_V1 / 2)
        .flat_map(|index| {
            let requirement = format!("fit/boundary-{index:04}");
            [
                link(
                    &requirement,
                    FitGdtEndpointRoleV1::Internal,
                    "gdt/boundary-internal",
                ),
                link(
                    &requirement,
                    FitGdtEndpointRoleV1::External,
                    "gdt/boundary-external",
                ),
            ]
        })
        .collect::<Vec<_>>();
    assert_eq!(links.len(), MAX_MACHINE_FIT_GDT_LINKS_V1);
    let admitted = MachineFitGdtCrosswalkDraftV1 { links }
        .admit_against(&fit, &geometric)
        .expect("exact endpoint-link cap admits");
    assert_eq!(admitted.endpoints().len(), MAX_MACHINE_FIT_GDT_LINKS_V1);
    assert_eq!(MACHINE_FIT_GDT_CROSSWALK_SCHEMA_VERSION_V1, 1);
}
