//! Machine-IR V2 assembly chronology/topology admission (Gauntlet G0/G3/G5).

use core::num::NonZeroU64;

use std::collections::BTreeSet;

use fs_blake3::identity::{ContentId, StrongIdentity};
use fs_blake3::{ContentHash, DomainHasher};
use fs_ir::IR_VERSION;
use fs_ir::machine::manufacturing::ManufacturingArtifactRefV1;
use fs_ir::machine::manufacturing::assembly::{
    AssemblyExecutionEvidenceRefV2, AssemblyJointFamilyV2, AssemblyLifecycleV2, AssemblyPathRefV2,
    AssemblyPreloadErrorV2, AssemblyPreloadUnitV2, AssemblyPreloadV2, AssemblyProcedureRefV2,
    AssemblyStepIdV2, AssemblyStepV2, AssemblyTopologyIssueV2, BoltStackParticipantV2,
    BoltStackRoleV2, JointFeatureUseIdV2, JointFeatureUseV2, JointOccurrenceIdV2,
    JointOccurrenceV2, JointTopologyV2, MACHINE_ASSEMBLY_AVAILABILITY_COMMITMENT_VERSION_V2,
    MACHINE_ASSEMBLY_IDENTITY_LIMITS_V2, MACHINE_ASSEMBLY_SCHEMA_VERSION_V2,
    MAX_MACHINE_ASSEMBLY_INITIAL_BODIES_V2, MAX_MACHINE_ASSEMBLY_INTRODUCTIONS_PER_STEP_V2,
    MAX_MACHINE_ASSEMBLY_OCCURRENCES_PER_STEP_V2, MAX_MACHINE_ASSEMBLY_OCCURRENCES_V2,
    MAX_MACHINE_ASSEMBLY_PARTICIPANTS_PER_OCCURRENCE_V2, MAX_MACHINE_ASSEMBLY_STEPS_V2,
    MachineAssemblyAdmissionErrorV2, MachineAssemblyDraftV2, MachineAssemblyIdV2,
    PhysicalFeatureUsePolicyV2,
};
use fs_ir::machine::{
    AdmittedMachineGraph, BodyId, ContactFeatureId, MachineGraphDraft, MaterialBinding,
    MaterialCardRef, MaterialTarget, ModelRef, SubsystemId, SubsystemSpec,
};

fn nz_v2(value: u64) -> NonZeroU64 {
    NonZeroU64::new(value).expect("fixture version must be nonzero")
}

fn body_v2(key: &str) -> BodyId {
    BodyId::new(key).unwrap_or_else(|error| panic!("fixture body key {key:?} must admit: {error}"))
}

fn feature_v2(key: &str) -> ContactFeatureId {
    ContactFeatureId::new(key)
        .unwrap_or_else(|error| panic!("fixture feature key {key:?} must admit: {error}"))
}

fn step_id_v2(key: &str) -> AssemblyStepIdV2 {
    AssemblyStepIdV2::new(key)
        .unwrap_or_else(|error| panic!("fixture step key {key:?} must admit: {error}"))
}

fn occurrence_id_v2(key: &str) -> JointOccurrenceIdV2 {
    JointOccurrenceIdV2::new(key)
        .unwrap_or_else(|error| panic!("fixture occurrence key {key:?} must admit: {error}"))
}

fn use_id_v2(key: &str) -> JointFeatureUseIdV2 {
    JointFeatureUseIdV2::new(key)
        .unwrap_or_else(|error| panic!("fixture feature-use key {key:?} must admit: {error}"))
}

fn artifact_v2(namespace: &str, byte: u8) -> ManufacturingArtifactRefV1 {
    artifact_coordinate_v2(namespace, 1, byte)
}

fn artifact_coordinate_v2(
    namespace: &str,
    schema_version: u64,
    byte: u8,
) -> ManufacturingArtifactRefV1 {
    ManufacturingArtifactRefV1::new(namespace, nz_v2(schema_version), ContentHash([byte; 32]))
        .unwrap_or_else(|error| panic!("fixture artifact {namespace:?} must admit: {error}"))
}

fn planned_v2(procedure_byte: u8, path_byte: u8) -> AssemblyLifecycleV2 {
    planned_coordinates_v2(
        "assembly-v2/procedure",
        1,
        procedure_byte,
        "assembly-v2/path",
        1,
        path_byte,
    )
}

fn planned_coordinates_v2(
    procedure_namespace: &str,
    procedure_schema_version: u64,
    procedure_byte: u8,
    path_namespace: &str,
    path_schema_version: u64,
    path_byte: u8,
) -> AssemblyLifecycleV2 {
    AssemblyLifecycleV2::Planned {
        procedure: AssemblyProcedureRefV2::new(artifact_coordinate_v2(
            procedure_namespace,
            procedure_schema_version,
            procedure_byte,
        )),
        path: AssemblyPathRefV2::new(artifact_coordinate_v2(
            path_namespace,
            path_schema_version,
            path_byte,
        )),
    }
}

fn execution_claimed_v2(
    procedure_byte: u8,
    path_byte: u8,
    evidence_byte: u8,
) -> AssemblyLifecycleV2 {
    execution_claimed_coordinates_v2(
        "assembly-v2/procedure",
        1,
        procedure_byte,
        "assembly-v2/path",
        1,
        path_byte,
        "assembly-v2/execution-evidence",
        1,
        evidence_byte,
    )
}

#[allow(clippy::too_many_arguments)]
fn execution_claimed_coordinates_v2(
    procedure_namespace: &str,
    procedure_schema_version: u64,
    procedure_byte: u8,
    path_namespace: &str,
    path_schema_version: u64,
    path_byte: u8,
    evidence_namespace: &str,
    evidence_schema_version: u64,
    evidence_byte: u8,
) -> AssemblyLifecycleV2 {
    AssemblyLifecycleV2::ExecutionClaimed {
        procedure: AssemblyProcedureRefV2::new(artifact_coordinate_v2(
            procedure_namespace,
            procedure_schema_version,
            procedure_byte,
        )),
        path: AssemblyPathRefV2::new(artifact_coordinate_v2(
            path_namespace,
            path_schema_version,
            path_byte,
        )),
        evidence: AssemblyExecutionEvidenceRefV2::new(artifact_coordinate_v2(
            evidence_namespace,
            evidence_schema_version,
            evidence_byte,
        )),
    }
}

fn preload_v2(value: f64, unit: AssemblyPreloadUnitV2) -> AssemblyPreloadV2 {
    AssemblyPreloadV2::try_new(value, unit).unwrap_or_else(|error| {
        panic!(
            "fixture preload {value} {} must admit: {error}",
            unit.symbol()
        )
    })
}

fn feature_use_v2(use_key: &str, body_key: &str, feature_key: &str) -> JointFeatureUseV2 {
    feature_use_with_policy_v2(
        use_key,
        body_key,
        feature_key,
        PhysicalFeatureUsePolicyV2::Reusable,
    )
}

fn feature_use_with_policy_v2(
    use_key: &str,
    body_key: &str,
    feature_key: &str,
    policy: PhysicalFeatureUsePolicyV2,
) -> JointFeatureUseV2 {
    JointFeatureUseV2::new(
        use_id_v2(use_key),
        fs_ir::machine::manufacturing::assembly::AssemblyFeatureSelectorV2::new(
            body_v2(body_key),
            feature_v2(feature_key),
        ),
        policy,
    )
}

fn material_v2(target: BodyId, key: &str, byte: u8) -> MaterialBinding {
    MaterialBinding {
        target: MaterialTarget::Body(target),
        material: MaterialCardRef::new(key, nz_v2(1), [byte; 32])
            .unwrap_or_else(|error| panic!("fixture material {key:?} must admit: {error}")),
    }
}

const ASSEMBLY_BODY_KEYS_V2: &[&str] = &[
    "body/base-a",
    "body/base-b",
    "body/bolt",
    "body/nut",
    "body/washer",
    "body/shaft",
    "body/hub",
    "body/key",
    "body/external",
    "body/internal",
    "body/spare",
];

fn assembly_feature_keys_v2() -> Vec<String> {
    ASSEMBLY_BODY_KEYS_V2
        .iter()
        .flat_map(|body_key| {
            let suffix = body_key.strip_prefix("body/").expect("fixture body prefix");
            [
                format!("contact/{suffix}/main"),
                format!("contact/{suffix}/alternate"),
            ]
        })
        .collect()
}

fn admitted_graph_v2(model_byte: u8) -> AdmittedMachineGraph {
    let assembly_bodies = ASSEMBLY_BODY_KEYS_V2
        .iter()
        .map(|key| body_v2(key))
        .collect::<Vec<_>>();
    let mut materials = assembly_bodies
        .iter()
        .cloned()
        .enumerate()
        .map(|(index, body)| {
            material_v2(
                body,
                &format!("materials/assembly-v2-{index}"),
                u8::try_from(index + 1).expect("fixture material byte fits"),
            )
        })
        .collect::<Vec<_>>();
    let other = body_v2("body/other");
    materials.push(material_v2(other.clone(), "materials/other-v2", 0x7e));

    MachineGraphDraft {
        clocks: Vec::new(),
        subsystems: vec![
            SubsystemSpec {
                id: SubsystemId::new("subsystem/assembly-v2").expect("canonical subsystem"),
                model: ModelRef::new("models/assembly-v2", nz_v2(1), [model_byte; 32])
                    .expect("canonical model"),
                bodies: assembly_bodies,
                surface_patches: Vec::new(),
                contact_features: assembly_feature_keys_v2()
                    .iter()
                    .map(|key| feature_v2(key))
                    .collect(),
                state_slots: Vec::new(),
            },
            SubsystemSpec {
                id: SubsystemId::new("subsystem/other-v2").expect("canonical subsystem"),
                model: ModelRef::new("models/other-v2", nz_v2(1), [0x7f; 32])
                    .expect("canonical model"),
                bodies: vec![other],
                surface_patches: Vec::new(),
                contact_features: vec![feature_v2("contact/other/main")],
                state_slots: Vec::new(),
            },
        ],
        terminals: Vec::new(),
        ports: Vec::new(),
        relations: Vec::new(),
        materials,
        interfaces: Vec::new(),
    }
    .admit()
    .expect("V2 assembly fixture graph must admit")
}

fn bolt_occurrence_v2() -> JointOccurrenceV2 {
    JointOccurrenceV2::new(
        occurrence_id_v2("occurrence/bolt"),
        JointTopologyV2::PreloadedBolt {
            clamped_members: vec![
                feature_use_v2("use/bolt/clamped-b", "body/base-b", "contact/base-b/main"),
                feature_use_v2("use/bolt/clamped-a", "body/base-a", "contact/base-a/main"),
            ],
            fastener_stack: vec![
                BoltStackParticipantV2::new(
                    2,
                    BoltStackRoleV2::Nut,
                    feature_use_v2("use/bolt/nut", "body/nut", "contact/nut/main"),
                ),
                BoltStackParticipantV2::new(
                    0,
                    BoltStackRoleV2::Bolt,
                    feature_use_v2("use/bolt/bolt", "body/bolt", "contact/bolt/main"),
                ),
                BoltStackParticipantV2::new(
                    1,
                    BoltStackRoleV2::Washer,
                    feature_use_v2("use/bolt/washer", "body/washer", "contact/washer/main"),
                ),
            ],
            preload: preload_v2(2.0, AssemblyPreloadUnitV2::Kilonewton),
        },
        planned_v2(0x10, 0x11),
    )
}

fn weld_occurrence_v2() -> JointOccurrenceV2 {
    JointOccurrenceV2::new(
        occurrence_id_v2("occurrence/weld"),
        JointTopologyV2::Weld {
            members: vec![
                feature_use_v2("use/weld/base-a", "body/base-a", "contact/base-a/main"),
                feature_use_v2("use/weld/base-b", "body/base-b", "contact/base-b/main"),
            ],
        },
        execution_claimed_v2(0x20, 0x21, 0x22),
    )
}

fn adhesive_occurrence_v2() -> JointOccurrenceV2 {
    JointOccurrenceV2::new(
        occurrence_id_v2("occurrence/adhesive"),
        JointTopologyV2::AdhesiveBond {
            adherends: vec![
                feature_use_v2("use/adhesive/base-a", "body/base-a", "contact/base-a/main"),
                feature_use_v2("use/adhesive/base-b", "body/base-b", "contact/base-b/main"),
            ],
        },
        planned_v2(0x30, 0x31),
    )
}

fn key_occurrence_v2() -> JointOccurrenceV2 {
    JointOccurrenceV2::new(
        occurrence_id_v2("occurrence/key"),
        JointTopologyV2::Key {
            shaft: feature_use_v2("use/key/shaft", "body/shaft", "contact/shaft/main"),
            hub: feature_use_v2("use/key/hub", "body/hub", "contact/hub/main"),
            key: feature_use_v2("use/key/body", "body/key", "contact/key/main"),
        },
        planned_v2(0x40, 0x41),
    )
}

fn spline_occurrence_v2() -> JointOccurrenceV2 {
    JointOccurrenceV2::new(
        occurrence_id_v2("occurrence/spline"),
        JointTopologyV2::Spline {
            external: feature_use_v2(
                "use/spline/external",
                "body/external",
                "contact/external/main",
            ),
            internal: feature_use_v2(
                "use/spline/internal",
                "body/internal",
                "contact/internal/main",
            ),
        },
        planned_v2(0x50, 0x51),
    )
}

fn interference_occurrence_v2() -> JointOccurrenceV2 {
    JointOccurrenceV2::new(
        occurrence_id_v2("occurrence/interference"),
        JointTopologyV2::InterferenceFit {
            external: feature_use_v2(
                "use/interference/external",
                "body/external",
                "contact/external/main",
            ),
            internal: feature_use_v2(
                "use/interference/internal",
                "body/internal",
                "contact/internal/main",
            ),
        },
        execution_claimed_v2(0x60, 0x61, 0x62),
    )
}

fn valid_draft_v2() -> MachineAssemblyDraftV2 {
    MachineAssemblyDraftV2 {
        initial_available_bodies: vec![
            body_v2("body/internal"),
            body_v2("body/base-b"),
            body_v2("body/shaft"),
            body_v2("body/external"),
            body_v2("body/base-a"),
            body_v2("body/hub"),
        ],
        steps: vec![
            AssemblyStepV2::new(
                step_id_v2("step/key"),
                2,
                vec![body_v2("body/key")],
                vec![occurrence_id_v2("occurrence/key")],
            ),
            AssemblyStepV2::new(
                step_id_v2("step/bolt"),
                0,
                vec![
                    body_v2("body/washer"),
                    body_v2("body/bolt"),
                    body_v2("body/nut"),
                ],
                vec![occurrence_id_v2("occurrence/bolt")],
            ),
            AssemblyStepV2::new(
                step_id_v2("step/directed"),
                3,
                Vec::new(),
                vec![
                    occurrence_id_v2("occurrence/spline"),
                    occurrence_id_v2("occurrence/interference"),
                ],
            ),
            AssemblyStepV2::new(
                step_id_v2("step/hybrid"),
                1,
                Vec::new(),
                vec![
                    occurrence_id_v2("occurrence/weld"),
                    occurrence_id_v2("occurrence/adhesive"),
                ],
            ),
        ],
        occurrences: vec![
            spline_occurrence_v2(),
            bolt_occurrence_v2(),
            interference_occurrence_v2(),
            key_occurrence_v2(),
            adhesive_occurrence_v2(),
            weld_occurrence_v2(),
        ],
    }
}

fn oracle_append_bytes_v2(out: &mut Vec<u8>, bytes: &[u8]) {
    out.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    out.extend_from_slice(bytes);
}

fn oracle_append_rows_v2(out: &mut Vec<u8>, rows: impl IntoIterator<Item = Vec<u8>>) {
    let rows = rows.into_iter().collect::<Vec<_>>();
    out.extend_from_slice(&(rows.len() as u64).to_le_bytes());
    for row in rows {
        oracle_append_bytes_v2(out, &row);
    }
}

const ORACLE_CANONICAL_FRAME_VERSION_V1: u32 = 1;
const ORACLE_PROBLEM_SEMANTIC_ROLE_TAG_V1: u8 = 12;
const ORACLE_FINITE_EXACT_BITS_POLICY_TAG_V1: u8 = 1;
const ORACLE_REQUIRED_PRESENCE_TAG_V1: u8 = 1;
const ORACLE_U64_WIRE_TAG_V1: u8 = 3;
const ORACLE_BYTES_WIRE_TAG_V1: u8 = 2;
const ORACLE_ORDERED_BYTES_WIRE_TAG_V1: u8 = 8;
const ORACLE_ASSEMBLY_SCHEMA_VERSION_V2: u32 = 2;
const ORACLE_AVAILABILITY_COMMITMENT_VERSION_V2: u32 = 1;
const ORACLE_ASSEMBLY_DOMAIN_V2: &str = "org.frankensim.fs-ir.machine.manufacturing-assembly.v2";
const ORACLE_ASSEMBLY_NAME_V2: &str = "admitted-machine-assembly";
const ORACLE_ASSEMBLY_CONTEXT_V2: &str = "one exact Machine graph, initial availability set, family-specific physical occurrences, and versioned authenticated chronological availability-transition chain";
const ORACLE_INITIAL_AVAILABILITY_DOMAIN_V2: &str = "org.frankensim.fs-ir.machine.manufacturing-assembly.v2/availability-transition-chain.v1/initial";
const ORACLE_STEP_AVAILABILITY_DOMAIN_V2: &str =
    "org.frankensim.fs-ir.machine.manufacturing-assembly.v2/availability-transition-chain.v1/step";
const ORACLE_ASSEMBLY_FIELDS_V2: &[(&str, u8)] = &[
    ("assembly-schema-version", ORACLE_U64_WIRE_TAG_V1),
    ("frankenscript-ir-version", ORACLE_U64_WIRE_TAG_V1),
    ("machine-graph", ORACLE_BYTES_WIRE_TAG_V1),
    ("initial-available-bodies", ORACLE_ORDERED_BYTES_WIRE_TAG_V1),
    ("joint-occurrences", ORACLE_ORDERED_BYTES_WIRE_TAG_V1),
    ("assembly-steps", ORACLE_ORDERED_BYTES_WIRE_TAG_V1),
];

#[derive(Debug, Clone)]
struct OracleAssemblyReceiptV2 {
    identity: MachineAssemblyIdV2,
    canonical_preimage: ContentId,
    canonical_frame: Vec<u8>,
    schema_id: [u8; 32],
    collection_items: u64,
    initial_rows: Vec<Vec<u8>>,
    occurrence_rows: Vec<Vec<u8>>,
    step_rows: Vec<Vec<u8>>,
    availability_roots: Vec<(ContentId, ContentId)>,
    initial_availability_preimage_bytes: u64,
    transition_preimage_bytes: u64,
    availability_rows_hashed: u64,
}

fn oracle_append_field_prefix_v2(out: &mut Vec<u8>, ordinal: u32, name: &str, wire_tag: u8) {
    out.push(0xf0);
    out.extend_from_slice(&ordinal.to_le_bytes());
    oracle_append_bytes_v2(out, name.as_bytes());
    out.extend_from_slice(&[wire_tag, ORACLE_REQUIRED_PRESENCE_TAG_V1]);
}

fn oracle_append_u64_field_v2(out: &mut Vec<u8>, ordinal: u32, name: &str, value: u64) {
    oracle_append_field_prefix_v2(out, ordinal, name, ORACLE_U64_WIRE_TAG_V1);
    out.extend_from_slice(&value.to_le_bytes());
}

fn oracle_append_bytes_field_v2(out: &mut Vec<u8>, ordinal: u32, name: &str, value: &[u8]) {
    oracle_append_field_prefix_v2(out, ordinal, name, ORACLE_BYTES_WIRE_TAG_V1);
    oracle_append_bytes_v2(out, value);
}

fn oracle_append_ordered_bytes_field_v2(
    out: &mut Vec<u8>,
    ordinal: u32,
    name: &str,
    rows: &[Vec<u8>],
) {
    oracle_append_field_prefix_v2(out, ordinal, name, ORACLE_ORDERED_BYTES_WIRE_TAG_V1);
    oracle_append_rows_v2(out, rows.iter().cloned());
}

fn oracle_schema_id_v2() -> [u8; 32] {
    let mut descriptor = Vec::new();
    descriptor.extend_from_slice(b"FSSCHEM\x02");
    descriptor.extend_from_slice(&ORACLE_CANONICAL_FRAME_VERSION_V1.to_le_bytes());
    oracle_append_bytes_v2(&mut descriptor, ORACLE_ASSEMBLY_DOMAIN_V2.as_bytes());
    oracle_append_bytes_v2(&mut descriptor, ORACLE_ASSEMBLY_NAME_V2.as_bytes());
    descriptor.extend_from_slice(&ORACLE_ASSEMBLY_SCHEMA_VERSION_V2.to_le_bytes());
    oracle_append_bytes_v2(&mut descriptor, ORACLE_ASSEMBLY_CONTEXT_V2.as_bytes());
    descriptor.extend_from_slice(&(ORACLE_ASSEMBLY_FIELDS_V2.len() as u64).to_le_bytes());
    for (name, wire_tag) in ORACLE_ASSEMBLY_FIELDS_V2 {
        oracle_append_bytes_v2(&mut descriptor, name.as_bytes());
        descriptor.extend_from_slice(&[*wire_tag, ORACLE_REQUIRED_PRESENCE_TAG_V1, 0]);
    }
    let mut hasher = DomainHasher::new("org.frankensim.fs-blake3.schema-id.v1");
    hasher.update(&descriptor);
    *hasher.finalize().as_bytes()
}

fn oracle_frame_v2(
    graph_id: &[u8; 32],
    assembly_schema_version: u32,
    ir_version: u32,
    initial_rows: &[Vec<u8>],
    occurrence_rows: &[Vec<u8>],
    step_rows: &[Vec<u8>],
) -> Vec<u8> {
    let schema_id = oracle_schema_id_v2();
    let mut frame = Vec::new();
    frame.extend_from_slice(b"FSID\0\0\0\x01");
    frame.extend_from_slice(&ORACLE_CANONICAL_FRAME_VERSION_V1.to_le_bytes());
    frame.extend_from_slice(&[
        ORACLE_PROBLEM_SEMANTIC_ROLE_TAG_V1,
        ORACLE_FINITE_EXACT_BITS_POLICY_TAG_V1,
    ]);
    oracle_append_bytes_v2(&mut frame, ORACLE_ASSEMBLY_DOMAIN_V2.as_bytes());
    oracle_append_bytes_v2(&mut frame, ORACLE_ASSEMBLY_NAME_V2.as_bytes());
    frame.extend_from_slice(&schema_id);
    frame.extend_from_slice(&ORACLE_ASSEMBLY_SCHEMA_VERSION_V2.to_le_bytes());
    oracle_append_bytes_v2(&mut frame, ORACLE_ASSEMBLY_CONTEXT_V2.as_bytes());
    frame.extend_from_slice(&(ORACLE_ASSEMBLY_FIELDS_V2.len() as u32).to_le_bytes());
    for (ordinal, (name, wire_tag)) in ORACLE_ASSEMBLY_FIELDS_V2.iter().enumerate() {
        frame.extend_from_slice(&(ordinal as u32).to_le_bytes());
        oracle_append_bytes_v2(&mut frame, name.as_bytes());
        frame.extend_from_slice(&[*wire_tag, ORACLE_REQUIRED_PRESENCE_TAG_V1]);
    }
    oracle_append_u64_field_v2(
        &mut frame,
        0,
        "assembly-schema-version",
        u64::from(assembly_schema_version),
    );
    oracle_append_u64_field_v2(
        &mut frame,
        1,
        "frankenscript-ir-version",
        u64::from(ir_version),
    );
    oracle_append_bytes_field_v2(&mut frame, 2, "machine-graph", graph_id);
    oracle_append_ordered_bytes_field_v2(&mut frame, 3, "initial-available-bodies", initial_rows);
    oracle_append_ordered_bytes_field_v2(&mut frame, 4, "joint-occurrences", occurrence_rows);
    oracle_append_ordered_bytes_field_v2(&mut frame, 5, "assembly-steps", step_rows);
    frame.push(0xff);
    frame.extend_from_slice(&(ORACLE_ASSEMBLY_FIELDS_V2.len() as u32).to_le_bytes());
    frame
}

fn oracle_identity_from_frame_v2(frame: &[u8]) -> MachineAssemblyIdV2 {
    let mut hasher = DomainHasher::new("org.frankensim.fs-blake3.canonical-identity-frame.v1");
    hasher.update(frame);
    MachineAssemblyIdV2::parse_slice(hasher.finalize().as_bytes())
        .expect("a BLAKE3 digest is a typed assembly identity")
}

fn oracle_body_row_v2(body: &BodyId) -> Vec<u8> {
    let mut row = Vec::new();
    oracle_append_bytes_v2(&mut row, body.identity().as_bytes());
    oracle_append_bytes_v2(&mut row, body.canonical_key().as_bytes());
    row
}

fn oracle_initial_availability_root_v2(bodies: &[BodyId]) -> (ContentId, u64) {
    let mut preimage = Vec::new();
    oracle_append_bytes_v2(
        &mut preimage,
        ORACLE_INITIAL_AVAILABILITY_DOMAIN_V2.as_bytes(),
    );
    preimage.extend_from_slice(&ORACLE_AVAILABILITY_COMMITMENT_VERSION_V2.to_le_bytes());
    oracle_append_rows_v2(&mut preimage, bodies.iter().map(oracle_body_row_v2));
    let bytes = preimage.len() as u64;
    (ContentId::of_bytes(&preimage), bytes)
}

fn oracle_feature_row_v2(feature: &ContactFeatureId) -> Vec<u8> {
    let mut row = Vec::new();
    oracle_append_bytes_v2(&mut row, feature.identity().as_bytes());
    oracle_append_bytes_v2(&mut row, feature.canonical_key().as_bytes());
    row
}

fn oracle_artifact_row_v2(artifact: &ManufacturingArtifactRefV1) -> Vec<u8> {
    let mut row = Vec::new();
    oracle_append_bytes_v2(&mut row, artifact.namespace().as_bytes());
    row.extend_from_slice(&artifact.schema_version().get().to_le_bytes());
    row.extend_from_slice(artifact.content_hash().as_bytes());
    row
}

fn oracle_use_row_v2(feature_use: &JointFeatureUseV2) -> Vec<u8> {
    let mut selector = Vec::new();
    oracle_append_bytes_v2(
        &mut selector,
        feature_use.selector().declared_body().identity().as_bytes(),
    );
    oracle_append_bytes_v2(
        &mut selector,
        feature_use
            .selector()
            .declared_body()
            .canonical_key()
            .as_bytes(),
    );
    let feature_row = oracle_feature_row_v2(feature_use.selector().contact_feature());
    selector.extend_from_slice(&feature_row);

    let mut row = Vec::new();
    oracle_append_bytes_v2(&mut row, feature_use.id().canonical_key().as_bytes());
    row.push(match feature_use.policy() {
        PhysicalFeatureUsePolicyV2::Reusable => 1,
        PhysicalFeatureUsePolicyV2::ExclusiveWithinAssembly => 2,
    });
    oracle_append_bytes_v2(&mut row, &selector);
    row
}

fn oracle_lifecycle_row_v2(lifecycle: &AssemblyLifecycleV2) -> Vec<u8> {
    let mut row = Vec::new();
    match lifecycle {
        AssemblyLifecycleV2::Planned { procedure, path } => {
            row.push(1);
            oracle_append_bytes_v2(&mut row, &oracle_artifact_row_v2(procedure.artifact()));
            oracle_append_bytes_v2(&mut row, &oracle_artifact_row_v2(path.artifact()));
        }
        AssemblyLifecycleV2::ExecutionClaimed {
            procedure,
            path,
            evidence,
        } => {
            row.push(2);
            oracle_append_bytes_v2(&mut row, &oracle_artifact_row_v2(procedure.artifact()));
            oracle_append_bytes_v2(&mut row, &oracle_artifact_row_v2(path.artifact()));
            oracle_append_bytes_v2(&mut row, &oracle_artifact_row_v2(evidence.artifact()));
        }
    }
    row
}

fn oracle_topology_row_v2(topology: &JointTopologyV2) -> Vec<u8> {
    let mut row = Vec::new();
    match topology {
        JointTopologyV2::PreloadedBolt {
            clamped_members,
            fastener_stack,
            preload,
        } => {
            row.push(1);
            let mut clamped_members = clamped_members.iter().collect::<Vec<_>>();
            clamped_members.sort_by_cached_key(|member| oracle_use_row_v2(member));
            oracle_append_rows_v2(&mut row, clamped_members.into_iter().map(oracle_use_row_v2));
            let mut fastener_stack = fastener_stack.iter().collect::<Vec<_>>();
            fastener_stack.sort_by(|left, right| {
                left.position().cmp(&right.position()).then_with(|| {
                    oracle_stack_participant_row_v2(left)
                        .cmp(&oracle_stack_participant_row_v2(right))
                })
            });
            oracle_append_rows_v2(
                &mut row,
                fastener_stack
                    .into_iter()
                    .map(oracle_stack_participant_row_v2),
            );
            row.extend_from_slice(&preload.submitted_bits().to_le_bytes());
            row.push(match preload.unit() {
                AssemblyPreloadUnitV2::Newton => 1,
                AssemblyPreloadUnitV2::Kilonewton => 2,
            });
            row.extend_from_slice(&preload.newtons_bits().to_le_bytes());
        }
        JointTopologyV2::Weld { members } => {
            row.push(2);
            let mut members = members.iter().collect::<Vec<_>>();
            members.sort_by_cached_key(|member| oracle_use_row_v2(member));
            oracle_append_rows_v2(&mut row, members.into_iter().map(oracle_use_row_v2));
        }
        JointTopologyV2::AdhesiveBond { adherends } => {
            row.push(3);
            let mut adherends = adherends.iter().collect::<Vec<_>>();
            adherends.sort_by_cached_key(|member| oracle_use_row_v2(member));
            oracle_append_rows_v2(&mut row, adherends.into_iter().map(oracle_use_row_v2));
        }
        JointTopologyV2::Key { shaft, hub, key } => {
            row.push(4);
            oracle_append_bytes_v2(&mut row, &oracle_use_row_v2(shaft));
            oracle_append_bytes_v2(&mut row, &oracle_use_row_v2(hub));
            oracle_append_bytes_v2(&mut row, &oracle_use_row_v2(key));
        }
        JointTopologyV2::Spline { external, internal } => {
            row.push(5);
            oracle_append_bytes_v2(&mut row, &oracle_use_row_v2(external));
            oracle_append_bytes_v2(&mut row, &oracle_use_row_v2(internal));
        }
        JointTopologyV2::InterferenceFit { external, internal } => {
            row.push(6);
            oracle_append_bytes_v2(&mut row, &oracle_use_row_v2(external));
            oracle_append_bytes_v2(&mut row, &oracle_use_row_v2(internal));
        }
    }
    row
}

fn oracle_stack_participant_row_v2(participant: &BoltStackParticipantV2) -> Vec<u8> {
    let mut row = Vec::new();
    row.extend_from_slice(&participant.position().to_le_bytes());
    row.push(match participant.role() {
        BoltStackRoleV2::Bolt => 1,
        BoltStackRoleV2::Nut => 2,
        BoltStackRoleV2::Washer => 3,
        BoltStackRoleV2::Spacer => 4,
        BoltStackRoleV2::LockingElement => 5,
    });
    oracle_append_bytes_v2(&mut row, &oracle_use_row_v2(participant.feature_use()));
    row
}

fn oracle_occurrence_row_v2(occurrence: &JointOccurrenceV2) -> Vec<u8> {
    let mut row = Vec::new();
    oracle_append_bytes_v2(&mut row, occurrence.id().canonical_key().as_bytes());
    oracle_append_bytes_v2(&mut row, &oracle_topology_row_v2(occurrence.topology()));
    oracle_append_bytes_v2(&mut row, &oracle_lifecycle_row_v2(occurrence.lifecycle()));
    row
}

fn oracle_transition_root_v2(
    prior_root: ContentId,
    available_before_count: usize,
    step: &AssemblyStepV2,
    available_after_count: usize,
) -> (ContentId, u64) {
    let mut preimage = Vec::new();
    oracle_append_bytes_v2(&mut preimage, ORACLE_STEP_AVAILABILITY_DOMAIN_V2.as_bytes());
    preimage.extend_from_slice(&ORACLE_AVAILABILITY_COMMITMENT_VERSION_V2.to_le_bytes());
    oracle_append_bytes_v2(&mut preimage, prior_root.as_bytes());
    preimage.extend_from_slice(&(available_before_count as u64).to_le_bytes());
    oracle_append_bytes_v2(&mut preimage, step.id().canonical_key().as_bytes());
    preimage.extend_from_slice(&step.ordinal().to_le_bytes());
    oracle_append_rows_v2(
        &mut preimage,
        step.introduced_bodies().iter().map(oracle_body_row_v2),
    );
    preimage.extend_from_slice(&(available_after_count as u64).to_le_bytes());
    let bytes = preimage.len() as u64;
    (ContentId::of_bytes(&preimage), bytes)
}

fn oracle_step_row_v2(
    step: &AssemblyStepV2,
    available_before_count: usize,
    availability_before_root: ContentId,
    available_after_count: usize,
    availability_after_root: ContentId,
) -> Vec<u8> {
    let mut row = Vec::new();
    oracle_append_bytes_v2(&mut row, step.id().canonical_key().as_bytes());
    row.extend_from_slice(&step.ordinal().to_le_bytes());
    oracle_append_rows_v2(
        &mut row,
        step.introduced_bodies().iter().map(oracle_body_row_v2),
    );
    oracle_append_rows_v2(
        &mut row,
        step.occurrence_ids().iter().map(|id| {
            let mut id_row = Vec::new();
            oracle_append_bytes_v2(&mut id_row, id.canonical_key().as_bytes());
            id_row
        }),
    );
    row.extend_from_slice(&(available_before_count as u64).to_le_bytes());
    oracle_append_bytes_v2(&mut row, availability_before_root.as_bytes());
    row.extend_from_slice(&(available_after_count as u64).to_le_bytes());
    oracle_append_bytes_v2(&mut row, availability_after_root.as_bytes());
    row
}

fn oracle_receipt_v2(
    graph_id: &[u8; 32],
    draft: &MachineAssemblyDraftV2,
    assembly_schema_version: u32,
    ir_version: u32,
) -> OracleAssemblyReceiptV2 {
    let mut initial_bodies = draft.initial_available_bodies.clone();
    initial_bodies
        .sort_by(|left, right| left.identity().as_bytes().cmp(right.identity().as_bytes()));
    let initial_rows = initial_bodies
        .iter()
        .map(oracle_body_row_v2)
        .collect::<Vec<_>>();

    let mut occurrences = draft.occurrences.clone();
    occurrences.sort_by(|left, right| {
        left.id()
            .canonical_key()
            .as_bytes()
            .cmp(right.id().canonical_key().as_bytes())
    });
    let occurrence_rows = occurrences
        .iter()
        .map(oracle_occurrence_row_v2)
        .collect::<Vec<_>>();

    let mut steps = draft
        .steps
        .iter()
        .map(|step| {
            let mut introduced_bodies = step.introduced_bodies().to_vec();
            introduced_bodies
                .sort_by(|left, right| left.identity().as_bytes().cmp(right.identity().as_bytes()));
            let mut occurrence_ids = step.occurrence_ids().to_vec();
            occurrence_ids.sort_by(|left, right| {
                left.canonical_key()
                    .as_bytes()
                    .cmp(right.canonical_key().as_bytes())
            });
            AssemblyStepV2::new(
                step.id().clone(),
                step.ordinal(),
                introduced_bodies,
                occurrence_ids,
            )
        })
        .collect::<Vec<_>>();
    steps.sort_by(|left, right| {
        left.ordinal().cmp(&right.ordinal()).then_with(|| {
            left.id()
                .canonical_key()
                .as_bytes()
                .cmp(right.id().canonical_key().as_bytes())
        })
    });

    let (initial_root, initial_preimage_bytes) =
        oracle_initial_availability_root_v2(&initial_bodies);
    let mut prior_root = initial_root;
    let mut available_count = initial_bodies.len();
    let mut transition_preimage_bytes = 0_u64;
    let mut availability_rows_hashed = initial_bodies.len() as u64;
    let mut availability_roots = Vec::with_capacity(steps.len());
    let mut step_rows = Vec::with_capacity(steps.len());
    for step in &steps {
        let after_count = available_count + step.introduced_bodies().len();
        let (after_root, transition_bytes) =
            oracle_transition_root_v2(prior_root, available_count, step, after_count);
        transition_preimage_bytes += transition_bytes;
        availability_rows_hashed += step.introduced_bodies().len() as u64;
        availability_roots.push((prior_root, after_root));
        step_rows.push(oracle_step_row_v2(
            step,
            available_count,
            prior_root,
            after_count,
            after_root,
        ));
        available_count = after_count;
        prior_root = after_root;
    }

    let canonical_frame = oracle_frame_v2(
        graph_id,
        assembly_schema_version,
        ir_version,
        &initial_rows,
        &occurrence_rows,
        &step_rows,
    );
    let collection_items = (initial_rows.len() + occurrence_rows.len() + step_rows.len()) as u64;
    OracleAssemblyReceiptV2 {
        identity: oracle_identity_from_frame_v2(&canonical_frame),
        canonical_preimage: ContentId::of_bytes(&canonical_frame),
        canonical_frame,
        schema_id: oracle_schema_id_v2(),
        collection_items,
        initial_rows,
        occurrence_rows,
        step_rows,
        availability_roots,
        initial_availability_preimage_bytes: initial_preimage_bytes,
        transition_preimage_bytes,
        availability_rows_hashed,
    }
}

fn singleton_draft_v2(
    occurrence: JointOccurrenceV2,
    initial_available_bodies: Vec<BodyId>,
    introduced_bodies: Vec<BodyId>,
) -> MachineAssemblyDraftV2 {
    let occurrence_id = occurrence.id().clone();
    MachineAssemblyDraftV2 {
        initial_available_bodies,
        steps: vec![AssemblyStepV2::new(
            step_id_v2("step/single"),
            0,
            introduced_bodies,
            vec![occurrence_id],
        )],
        occurrences: vec![occurrence],
    }
}

#[test]
fn mas2_001_all_families_multi_body_chronology_and_independent_oracle_agree() {
    let graph = admitted_graph_v2(0x41);
    let draft = valid_draft_v2();
    let oracle = oracle_receipt_v2(
        graph.identity().as_bytes(),
        &draft,
        MACHINE_ASSEMBLY_SCHEMA_VERSION_V2,
        IR_VERSION,
    );
    let admitted = draft
        .admit_against(&graph)
        .expect("complete V2 family fixture must admit");

    let families = admitted
        .occurrences()
        .iter()
        .map(|occurrence| occurrence.topology().family())
        .collect::<BTreeSet<_>>();
    assert_eq!(
        families,
        BTreeSet::from([
            AssemblyJointFamilyV2::PreloadedBolt,
            AssemblyJointFamilyV2::Weld,
            AssemblyJointFamilyV2::AdhesiveBond,
            AssemblyJointFamilyV2::Key,
            AssemblyJointFamilyV2::Spline,
            AssemblyJointFamilyV2::InterferenceFit,
        ]),
        "every closed family must survive V2 admission"
    );
    assert_eq!(
        admitted.steps()[0].step().introduced_bodies().len(),
        3,
        "bolt/nut/washer must enter atomically in one multi-body step"
    );
    assert!(
        admitted.steps()[1].step().introduced_bodies().is_empty(),
        "hybrid weld/adhesive continuation must require no chronology-derived endpoint order"
    );
    assert_eq!(
        admitted.steps()[2].step().introduced_bodies(),
        &[body_v2("body/key")],
        "the physically distinct key body must enter in its own one-body transition"
    );
    assert_eq!(
        admitted.steps()[0].available_before_count() + 3,
        admitted.steps()[0].available_after_count(),
        "availability may publish only after the complete bolt-stack step validates"
    );

    assert_eq!(
        admitted.identity(),
        oracle.identity,
        "independently serialized canonical rows must reproduce the production identity"
    );
    assert_eq!(
        admitted.identity_receipt().canonical_preimage(),
        oracle.canonical_preimage,
        "independent preimage oracle must match, not merely a second production replay"
    );
    assert_eq!(
        admitted.identity_receipt().canonical_bytes(),
        oracle.canonical_frame.len() as u64,
        "independent frame accounting must pin the exact canonical byte count"
    );
    assert_eq!(
        admitted.identity_receipt().schema_id().as_bytes(),
        &oracle.schema_id,
        "independent schema-descriptor framing must reproduce the production schema ID"
    );
    assert_eq!(
        admitted.identity_receipt().collection_items(),
        oracle.collection_items,
        "independent collection accounting must match the production receipt"
    );
    assert_eq!(
        oracle.identity.to_hex(),
        "cb9f0a76761cdd4ef74c10d74964626741a0a6141f60b4f1874e40cf73ee2f60",
        "reviewed complete V2 semantic-identity golden must remain frozen"
    );
    assert_eq!(
        oracle.canonical_preimage.to_hex(),
        "79fdd746e1860b4ad8d4f807f7b378cc10b26915b308fd32d7e09dc25068bdd0",
        "reviewed complete V2 canonical-preimage golden must remain frozen"
    );
    assert_eq!(
        ContentHash(oracle.schema_id).to_hex(),
        "a5e95139cd763a2e4a77133b28124f0f14f902978769ddcb61326685ab123897",
        "reviewed independently framed assembly schema ID must remain frozen"
    );
    assert_eq!(oracle.canonical_frame.len(), 6_367);
    assert_eq!(oracle.collection_items, 16);
    assert_eq!(oracle.availability_rows_hashed, 10);
    assert_eq!(oracle.initial_availability_preimage_bytes, 517);
    assert_eq!(oracle.transition_preimage_bytes, 1_021);
    assert_eq!(
        admitted.initial_availability_root().to_hex(),
        "b038fd342a606287633b88ecca0fec5555957183da25fbcff91149a3950d81d1"
    );
    assert_eq!(
        admitted
            .steps()
            .last()
            .expect("four steps")
            .availability_after_root()
            .to_hex(),
        "03169c137fa239857ce397eb07dffba8c2cfbe8ddf5304ae1ca7b883e39c4727"
    );

    let mut tag_drift_rows = oracle.occurrence_rows.clone();
    let first_id_len = u64::from_le_bytes(
        tag_drift_rows[0][..8]
            .try_into()
            .expect("occurrence ID length prefix"),
    ) as usize;
    let family_tag_offset = 8 + first_id_len + 8;
    tag_drift_rows[0][family_tag_offset] ^= 0x40;
    let tag_drift_frame = oracle_frame_v2(
        graph.identity().as_bytes(),
        MACHINE_ASSEMBLY_SCHEMA_VERSION_V2,
        IR_VERSION,
        &oracle.initial_rows,
        &tag_drift_rows,
        &oracle.step_rows,
    );
    assert_ne!(
        ContentId::of_bytes(&tag_drift_frame),
        oracle.canonical_preimage,
        "a hard-coded family-tag drift must fail the frozen preimage root"
    );

    let mut order_drift_rows = oracle.occurrence_rows.clone();
    order_drift_rows.reverse();
    let order_drift_frame = oracle_frame_v2(
        graph.identity().as_bytes(),
        MACHINE_ASSEMBLY_SCHEMA_VERSION_V2,
        IR_VERSION,
        &oracle.initial_rows,
        &order_drift_rows,
        &oracle.step_rows,
    );
    assert_ne!(
        oracle_identity_from_frame_v2(&order_drift_frame),
        oracle.identity,
        "ordered-row drift must fail the frozen semantic root"
    );

    let mut framing_drift_rows = oracle.occurrence_rows.clone();
    framing_drift_rows[0][0] ^= 1;
    let framing_drift_frame = oracle_frame_v2(
        graph.identity().as_bytes(),
        MACHINE_ASSEMBLY_SCHEMA_VERSION_V2,
        IR_VERSION,
        &oracle.initial_rows,
        &framing_drift_rows,
        &oracle.step_rows,
    );
    assert_ne!(
        ContentId::of_bytes(&framing_drift_frame),
        oracle.canonical_preimage,
        "row-length framing drift must fail the frozen preimage root"
    );
    assert_eq!(
        admitted.initial_availability_root(),
        oracle.availability_roots[0].0,
        "the retained chain seed must be the independently recomputed canonical initial-set root"
    );
    for (step, (before, after)) in admitted.steps().iter().zip(&oracle.availability_roots) {
        assert_eq!(step.availability_before_root(), *before);
        assert_eq!(step.availability_after_root(), *after);
    }

    let stale_schema = oracle_receipt_v2(
        graph.identity().as_bytes(),
        &valid_draft_v2(),
        MACHINE_ASSEMBLY_SCHEMA_VERSION_V2 - 1,
        IR_VERSION,
    );
    let stale_ir = oracle_receipt_v2(
        graph.identity().as_bytes(),
        &valid_draft_v2(),
        MACHINE_ASSEMBLY_SCHEMA_VERSION_V2,
        IR_VERSION - 1,
    );
    assert_ne!(
        oracle.identity, stale_schema.identity,
        "assembly schema version must be an independent semantic field"
    );
    assert_ne!(
        oracle.identity, stale_ir.identity,
        "FrankenScript IR version must be an independent semantic field"
    );
}

#[test]
fn mas2_002_unordered_sets_are_symmetric_while_directed_roles_remain_distinct() {
    let graph = admitted_graph_v2(0x42);
    let baseline = valid_draft_v2()
        .admit_against(&graph)
        .expect("baseline V2 fixture must admit");

    let mut permuted = valid_draft_v2();
    permuted.initial_available_bodies.reverse();
    permuted.steps = permuted
        .steps
        .into_iter()
        .map(|step| {
            let mut introduced = step.introduced_bodies().to_vec();
            let mut occurrence_ids = step.occurrence_ids().to_vec();
            introduced.reverse();
            occurrence_ids.reverse();
            AssemblyStepV2::new(
                step.id().clone(),
                step.ordinal(),
                introduced,
                occurrence_ids,
            )
        })
        .rev()
        .collect();
    permuted.occurrences = permuted
        .occurrences
        .into_iter()
        .map(|occurrence| {
            let mut topology = occurrence.topology().clone();
            match &mut topology {
                JointTopologyV2::PreloadedBolt {
                    clamped_members,
                    fastener_stack,
                    ..
                } => {
                    clamped_members.reverse();
                    fastener_stack.reverse();
                }
                JointTopologyV2::Weld { members } => members.reverse(),
                JointTopologyV2::AdhesiveBond { adherends } => adherends.reverse(),
                JointTopologyV2::Key { .. }
                | JointTopologyV2::Spline { .. }
                | JointTopologyV2::InterferenceFit { .. } => {}
            }
            JointOccurrenceV2::new(
                occurrence.id().clone(),
                topology,
                occurrence.lifecycle().clone(),
            )
        })
        .rev()
        .collect();
    let permuted = permuted
        .admit_against(&graph)
        .expect("symmetric permutations must admit");
    assert_eq!(
        baseline.identity(),
        permuted.identity(),
        "caller order and genuinely unordered hyperedge-member order must be non-semantic"
    );

    let mut reversed = valid_draft_v2();
    reversed.occurrences = reversed
        .occurrences
        .into_iter()
        .map(|occurrence| match occurrence.topology().clone() {
            JointTopologyV2::Spline { external, internal } => JointOccurrenceV2::new(
                occurrence.id().clone(),
                JointTopologyV2::Spline {
                    external: internal,
                    internal: external,
                },
                occurrence.lifecycle().clone(),
            ),
            topology => JointOccurrenceV2::new(
                occurrence.id().clone(),
                topology,
                occurrence.lifecycle().clone(),
            ),
        })
        .collect();
    let reversed = reversed
        .admit_against(&graph)
        .expect("directed role reversal remains structurally admissible");
    assert_ne!(
        baseline.identity(),
        reversed.identity(),
        "external/internal reversal must move identity even though chronology is unchanged"
    );
}

#[test]
fn mas2_003_lifecycle_is_truthful_and_equal_underlying_artifacts_are_allowed() {
    let graph = admitted_graph_v2(0x43);
    assert_eq!(
        planned_v2(0x01, 0x02).tag(),
        1,
        "the Planned lifecycle tag is a pinned V2 wire discriminant"
    );
    assert_eq!(
        execution_claimed_v2(0x01, 0x02, 0x03).tag(),
        2,
        "the ExecutionClaimed lifecycle tag is a pinned V2 wire discriminant"
    );
    assert_eq!(
        MACHINE_ASSEMBLY_AVAILABILITY_COMMITMENT_VERSION_V2,
        ORACLE_AVAILABILITY_COMMITMENT_VERSION_V2,
        "production and independently frozen availability-chain versions must agree"
    );
    for (family, tag) in [
        (AssemblyJointFamilyV2::PreloadedBolt, 1),
        (AssemblyJointFamilyV2::Weld, 2),
        (AssemblyJointFamilyV2::AdhesiveBond, 3),
        (AssemblyJointFamilyV2::Key, 4),
        (AssemblyJointFamilyV2::Spline, 5),
        (AssemblyJointFamilyV2::InterferenceFit, 6),
    ] {
        assert_eq!(family.tag(), tag, "every closed family tag is frozen");
    }
    for (role, tag) in [
        (BoltStackRoleV2::Bolt, 1),
        (BoltStackRoleV2::Nut, 2),
        (BoltStackRoleV2::Washer, 3),
        (BoltStackRoleV2::Spacer, 4),
        (BoltStackRoleV2::LockingElement, 5),
    ] {
        assert_eq!(role.tag(), tag, "every fastener-stack role tag is frozen");
    }
    assert_eq!(PhysicalFeatureUsePolicyV2::Reusable.tag(), 1);
    assert_eq!(PhysicalFeatureUsePolicyV2::ExclusiveWithinAssembly.tag(), 2);
    assert_eq!(AssemblyPreloadUnitV2::Newton.tag(), 1);
    assert_eq!(AssemblyPreloadUnitV2::Kilonewton.tag(), 2);
    let shared = artifact_v2("assembly-v2/shared-coordinate", 0x71);
    let planned = JointOccurrenceV2::new(
        occurrence_id_v2("occurrence/equal-artifacts"),
        JointTopologyV2::Weld {
            members: vec![
                feature_use_v2(
                    "use/equal/base-a",
                    "body/base-a",
                    "contact/base-a/alternate",
                ),
                feature_use_v2(
                    "use/equal/base-b",
                    "body/base-b",
                    "contact/base-b/alternate",
                ),
            ],
        },
        AssemblyLifecycleV2::Planned {
            procedure: AssemblyProcedureRefV2::new(shared.clone()),
            path: AssemblyPathRefV2::new(shared.clone()),
        },
    );
    let planned = singleton_draft_v2(
        planned,
        vec![body_v2("body/base-a"), body_v2("body/base-b")],
        Vec::new(),
    )
    .admit_against(&graph)
    .expect("equal underlying procedure/path coordinates are valid nominal role reuse");
    assert!(
        planned.occurrences()[0]
            .lifecycle()
            .execution_evidence()
            .is_none(),
        "Planned must not fabricate execution evidence"
    );

    let claimed = JointOccurrenceV2::new(
        occurrence_id_v2("occurrence/equal-artifacts"),
        JointTopologyV2::Weld {
            members: vec![
                feature_use_v2(
                    "use/equal/base-a",
                    "body/base-a",
                    "contact/base-a/alternate",
                ),
                feature_use_v2(
                    "use/equal/base-b",
                    "body/base-b",
                    "contact/base-b/alternate",
                ),
            ],
        },
        AssemblyLifecycleV2::ExecutionClaimed {
            procedure: AssemblyProcedureRefV2::new(shared.clone()),
            path: AssemblyPathRefV2::new(shared.clone()),
            evidence: AssemblyExecutionEvidenceRefV2::new(shared),
        },
    );
    let claimed = singleton_draft_v2(
        claimed,
        vec![body_v2("body/base-a"), body_v2("body/base-b")],
        Vec::new(),
    )
    .admit_against(&graph)
    .expect("equal coordinates across three typed claimed roles are allowed");
    assert!(
        claimed.occurrences()[0]
            .lifecycle()
            .execution_evidence()
            .is_some(),
        "ExecutionClaimed must bind an evidence coordinate"
    );
    assert_ne!(
        planned.identity(),
        claimed.identity(),
        "lifecycle discriminant and evidence payload must move identity"
    );
}

#[test]
fn mas2_004_reusable_features_support_hybrid_and_rework_occurrences() {
    let graph = admitted_graph_v2(0x44);
    let admitted = valid_draft_v2()
        .admit_against(&graph)
        .expect("reusable physical features across bolt/weld/adhesive occurrences must admit");
    let reuse_count = admitted
        .occurrences()
        .iter()
        .flat_map(|occurrence| occurrence.topology().participants())
        .filter(|(_, feature_use)| {
            feature_use.selector().contact_feature() == &feature_v2("contact/base-a/main")
        })
        .count();
    assert_eq!(
        reuse_count, 3,
        "one durable physical feature must support three separately identified uses"
    );
}

#[test]
fn mas2_005_duplicate_occurrence_and_use_identities_refuse_deterministically() {
    let graph = admitted_graph_v2(0x45);

    let mut duplicate_occurrence = valid_draft_v2();
    duplicate_occurrence
        .occurrences
        .push(duplicate_occurrence.occurrences[0].clone());
    let error = duplicate_occurrence
        .admit_against(&graph)
        .expect_err("duplicate physical occurrence identity must refuse");
    assert_eq!(
        error,
        MachineAssemblyAdmissionErrorV2::DuplicateOccurrence {
            occurrence: occurrence_id_v2("occurrence/spline"),
        },
        "sorted duplicate occurrence diagnosis must retain the exact durable ID"
    );
    assert_eq!(error.code(), "MachineAssemblyDuplicateOccurrence");

    let duplicate_use = JointOccurrenceV2::new(
        occurrence_id_v2("occurrence/duplicate-use"),
        JointTopologyV2::Weld {
            members: vec![
                feature_use_v2("use/duplicate", "body/base-a", "contact/base-a/alternate"),
                feature_use_v2("use/duplicate", "body/base-b", "contact/base-b/alternate"),
            ],
        },
        planned_v2(0x72, 0x73),
    );
    let error = singleton_draft_v2(
        duplicate_use,
        vec![body_v2("body/base-a"), body_v2("body/base-b")],
        Vec::new(),
    )
    .admit_against(&graph)
    .expect_err("duplicate feature-use identity must refuse before chronology publication");
    assert_eq!(
        error,
        MachineAssemblyAdmissionErrorV2::DuplicateFeatureUse {
            feature_use: use_id_v2("use/duplicate"),
            first: occurrence_id_v2("occurrence/duplicate-use"),
            duplicate: occurrence_id_v2("occurrence/duplicate-use"),
        },
        "duplicate-use refusal must retain both occurrence coordinates"
    );
}

#[test]
fn mas2_006_exclusive_policy_is_explicit_while_reusable_policy_is_not_one_shot() {
    let graph = admitted_graph_v2(0x46);
    let first = JointOccurrenceV2::new(
        occurrence_id_v2("occurrence/exclusive-a"),
        JointTopologyV2::Weld {
            members: vec![
                feature_use_v2(
                    "use/exclusive-a/base-a",
                    "body/base-a",
                    "contact/base-a/alternate",
                ),
                feature_use_v2(
                    "use/exclusive-a/base-b",
                    "body/base-b",
                    "contact/base-b/alternate",
                ),
            ],
        },
        planned_v2(0x74, 0x75),
    );
    let second = JointOccurrenceV2::new(
        occurrence_id_v2("occurrence/exclusive-b"),
        JointTopologyV2::AdhesiveBond {
            adherends: vec![
                feature_use_with_policy_v2(
                    "use/exclusive-b/base-a",
                    "body/base-a",
                    "contact/base-a/alternate",
                    PhysicalFeatureUsePolicyV2::ExclusiveWithinAssembly,
                ),
                feature_use_v2(
                    "use/exclusive-b/base-b",
                    "body/base-b",
                    "contact/base-b/alternate",
                ),
            ],
        },
        planned_v2(0x76, 0x77),
    );
    let draft = MachineAssemblyDraftV2 {
        initial_available_bodies: vec![body_v2("body/base-a"), body_v2("body/base-b")],
        steps: vec![AssemblyStepV2::new(
            step_id_v2("step/exclusive"),
            0,
            Vec::new(),
            vec![first.id().clone(), second.id().clone()],
        )],
        occurrences: vec![second, first],
    };
    let error = draft
        .admit_against(&graph)
        .expect_err("an explicit exclusive use must reject a competing reusable use");
    assert_eq!(
        error,
        MachineAssemblyAdmissionErrorV2::ExclusiveFeatureReuse {
            feature: feature_v2("contact/base-a/alternate"),
            first: use_id_v2("use/exclusive-a/base-a"),
            duplicate: use_id_v2("use/exclusive-b/base-a"),
        },
        "exclusive reuse diagnosis must be caller-order invariant"
    );
}

#[test]
fn mas2_007_coownership_does_not_fabricate_containment_authority() {
    let graph = admitted_graph_v2(0x47);
    let same_owner_wrong_body = JointOccurrenceV2::new(
        occurrence_id_v2("occurrence/no-containment"),
        JointTopologyV2::Weld {
            members: vec![
                feature_use_v2(
                    "use/no-containment/shaft",
                    "body/shaft",
                    "contact/base-a/main",
                ),
                feature_use_v2(
                    "use/no-containment/base-b",
                    "body/base-b",
                    "contact/base-b/main",
                ),
            ],
        },
        planned_v2(0x78, 0x79),
    );
    let admitted = singleton_draft_v2(
        same_owner_wrong_body,
        vec![body_v2("body/shaft"), body_v2("body/base-b")],
        Vec::new(),
    )
    .admit_against(&graph)
    .expect("same-subsystem body/feature selection cannot be rejected as false containment proof");
    let retained = admitted.occurrences()[0]
        .topology()
        .participants()
        .into_iter()
        .find(|(_, feature_use)| feature_use.id() == &use_id_v2("use/no-containment/shaft"))
        .expect("wrong-body selector use must remain in the receipt");
    assert_eq!(
        retained.1.selector().contact_feature(),
        &feature_v2("contact/base-a/main"),
        "receipt must retain the authority-free caller selector exactly"
    );

    let cross_owner = JointOccurrenceV2::new(
        occurrence_id_v2("occurrence/cross-owner"),
        JointTopologyV2::Weld {
            members: vec![
                feature_use_v2(
                    "use/cross-owner/base-a",
                    "body/base-a",
                    "contact/other/main",
                ),
                feature_use_v2(
                    "use/cross-owner/base-b",
                    "body/base-b",
                    "contact/base-b/main",
                ),
            ],
        },
        planned_v2(0x7a, 0x7b),
    );
    let error = singleton_draft_v2(
        cross_owner,
        vec![body_v2("body/base-a"), body_v2("body/base-b")],
        Vec::new(),
    )
    .admit_against(&graph)
    .expect_err("cross-owner selector must refuse even though containment remains unclaimed");
    assert!(
        matches!(
            error,
            MachineAssemblyAdmissionErrorV2::ParticipantOwnerMismatch {
                ref occurrence,
                ref body,
                ref feature,
                ..
            } if *occurrence == occurrence_id_v2("occurrence/cross-owner")
                && *body == body_v2("body/base-a")
                && *feature == feature_v2("contact/other/main")
        ),
        "cross-owner refusal must retain occurrence, body, and feature; got {error:?}"
    );
}

#[test]
fn mas2_008_closed_topology_invariants_and_preload_constructor_fail_closed() {
    assert_eq!(
        AssemblyPreloadV2::try_new(f64::NAN, AssemblyPreloadUnitV2::Newton),
        Err(AssemblyPreloadErrorV2::NonFinite),
        "NaN preload must refuse at construction"
    );
    assert_eq!(
        AssemblyPreloadV2::try_new(0.0, AssemblyPreloadUnitV2::Newton),
        Err(AssemblyPreloadErrorV2::NonPositive),
        "zero preload must refuse at construction"
    );
    assert_eq!(
        AssemblyPreloadV2::try_new(f64::MAX, AssemblyPreloadUnitV2::Kilonewton),
        Err(AssemblyPreloadErrorV2::SiNonFinite),
        "SI normalization overflow must refuse at construction"
    );

    let graph = admitted_graph_v2(0x48);
    let invalid_bolt = JointOccurrenceV2::new(
        occurrence_id_v2("occurrence/invalid-bolt"),
        JointTopologyV2::PreloadedBolt {
            clamped_members: vec![
                feature_use_v2(
                    "use/invalid-bolt/a",
                    "body/base-a",
                    "contact/base-a/alternate",
                ),
                feature_use_v2(
                    "use/invalid-bolt/b",
                    "body/base-b",
                    "contact/base-b/alternate",
                ),
            ],
            fastener_stack: vec![BoltStackParticipantV2::new(
                0,
                BoltStackRoleV2::Washer,
                feature_use_v2(
                    "use/invalid-bolt/washer",
                    "body/washer",
                    "contact/washer/alternate",
                ),
            )],
            preload: preload_v2(100.0, AssemblyPreloadUnitV2::Newton),
        },
        planned_v2(0x7c, 0x7d),
    );
    let error = singleton_draft_v2(
        invalid_bolt,
        vec![body_v2("body/base-a"), body_v2("body/base-b")],
        vec![body_v2("body/washer")],
    )
    .admit_against(&graph)
    .expect_err("preloaded-bolt topology without exactly one bolt must refuse");
    assert_eq!(
        error,
        MachineAssemblyAdmissionErrorV2::InvalidTopology {
            occurrence: occurrence_id_v2("occurrence/invalid-bolt"),
            issue: AssemblyTopologyIssueV2::BoltCount { actual: 0 },
        }
    );
}

fn mutation_weld_v2(
    occurrence_key: &str,
    use_a_key: &str,
    feature_a_key: &str,
    policy_a: PhysicalFeatureUsePolicyV2,
    topology_is_adhesive: bool,
    lifecycle: AssemblyLifecycleV2,
) -> JointOccurrenceV2 {
    mutation_weld_on_body_v2(
        occurrence_key,
        use_a_key,
        "body/base-a",
        feature_a_key,
        policy_a,
        topology_is_adhesive,
        lifecycle,
    )
}

fn mutation_weld_on_body_v2(
    occurrence_key: &str,
    use_a_key: &str,
    body_a_key: &str,
    feature_a_key: &str,
    policy_a: PhysicalFeatureUsePolicyV2,
    topology_is_adhesive: bool,
    lifecycle: AssemblyLifecycleV2,
) -> JointOccurrenceV2 {
    let left = feature_use_with_policy_v2(use_a_key, body_a_key, feature_a_key, policy_a);
    let right = feature_use_v2(
        "use/mutation/base-b",
        "body/base-b",
        "contact/base-b/alternate",
    );
    let topology = if topology_is_adhesive {
        JointTopologyV2::AdhesiveBond {
            adherends: vec![left, right],
        }
    } else {
        JointTopologyV2::Weld {
            members: vec![left, right],
        }
    };
    JointOccurrenceV2::new(occurrence_id_v2(occurrence_key), topology, lifecycle)
}

fn admit_mutation_weld_v2(
    graph: &AdmittedMachineGraph,
    occurrence: JointOccurrenceV2,
) -> MachineAssemblyIdV2 {
    singleton_draft_v2(
        occurrence,
        vec![
            body_v2("body/base-a"),
            body_v2("body/base-b"),
            body_v2("body/shaft"),
        ],
        Vec::new(),
    )
    .admit_against(graph)
    .expect("isolated semantic mutation fixture must remain admissible")
    .identity()
}

#[test]
fn mas2_009_isolated_semantic_mutations_move_every_occurrence_field() {
    let graph = admitted_graph_v2(0x49);
    let baseline = admit_mutation_weld_v2(
        &graph,
        mutation_weld_v2(
            "occurrence/mutation",
            "use/mutation/base-a",
            "contact/base-a/alternate",
            PhysicalFeatureUsePolicyV2::Reusable,
            false,
            planned_v2(0x80, 0x81),
        ),
    );

    let mutations = [
        (
            "occurrence-id",
            mutation_weld_v2(
                "occurrence/mutation-changed",
                "use/mutation/base-a",
                "contact/base-a/alternate",
                PhysicalFeatureUsePolicyV2::Reusable,
                false,
                planned_v2(0x80, 0x81),
            ),
        ),
        (
            "feature-use-id",
            mutation_weld_v2(
                "occurrence/mutation",
                "use/mutation/base-a-changed",
                "contact/base-a/alternate",
                PhysicalFeatureUsePolicyV2::Reusable,
                false,
                planned_v2(0x80, 0x81),
            ),
        ),
        (
            "feature-use-policy",
            mutation_weld_v2(
                "occurrence/mutation",
                "use/mutation/base-a",
                "contact/base-a/alternate",
                PhysicalFeatureUsePolicyV2::ExclusiveWithinAssembly,
                false,
                planned_v2(0x80, 0x81),
            ),
        ),
        (
            "declared-body-selector",
            mutation_weld_on_body_v2(
                "occurrence/mutation",
                "use/mutation/base-a",
                "body/shaft",
                "contact/base-a/alternate",
                PhysicalFeatureUsePolicyV2::Reusable,
                false,
                planned_v2(0x80, 0x81),
            ),
        ),
        (
            "physical-feature-selector",
            mutation_weld_v2(
                "occurrence/mutation",
                "use/mutation/base-a",
                "contact/base-a/main",
                PhysicalFeatureUsePolicyV2::Reusable,
                false,
                planned_v2(0x80, 0x81),
            ),
        ),
        (
            "family-payload",
            mutation_weld_v2(
                "occurrence/mutation",
                "use/mutation/base-a",
                "contact/base-a/alternate",
                PhysicalFeatureUsePolicyV2::Reusable,
                true,
                planned_v2(0x80, 0x81),
            ),
        ),
        (
            "procedure-artifact-namespace",
            mutation_weld_v2(
                "occurrence/mutation",
                "use/mutation/base-a",
                "contact/base-a/alternate",
                PhysicalFeatureUsePolicyV2::Reusable,
                false,
                planned_coordinates_v2(
                    "assembly-v2/procedure-changed",
                    1,
                    0x80,
                    "assembly-v2/path",
                    1,
                    0x81,
                ),
            ),
        ),
        (
            "procedure-artifact-schema-version",
            mutation_weld_v2(
                "occurrence/mutation",
                "use/mutation/base-a",
                "contact/base-a/alternate",
                PhysicalFeatureUsePolicyV2::Reusable,
                false,
                planned_coordinates_v2(
                    "assembly-v2/procedure",
                    2,
                    0x80,
                    "assembly-v2/path",
                    1,
                    0x81,
                ),
            ),
        ),
        (
            "procedure-artifact-content-hash",
            mutation_weld_v2(
                "occurrence/mutation",
                "use/mutation/base-a",
                "contact/base-a/alternate",
                PhysicalFeatureUsePolicyV2::Reusable,
                false,
                planned_v2(0x82, 0x81),
            ),
        ),
        (
            "path-artifact-namespace",
            mutation_weld_v2(
                "occurrence/mutation",
                "use/mutation/base-a",
                "contact/base-a/alternate",
                PhysicalFeatureUsePolicyV2::Reusable,
                false,
                planned_coordinates_v2(
                    "assembly-v2/procedure",
                    1,
                    0x80,
                    "assembly-v2/path-changed",
                    1,
                    0x81,
                ),
            ),
        ),
        (
            "path-artifact-schema-version",
            mutation_weld_v2(
                "occurrence/mutation",
                "use/mutation/base-a",
                "contact/base-a/alternate",
                PhysicalFeatureUsePolicyV2::Reusable,
                false,
                planned_coordinates_v2(
                    "assembly-v2/procedure",
                    1,
                    0x80,
                    "assembly-v2/path",
                    2,
                    0x81,
                ),
            ),
        ),
        (
            "path-artifact-content-hash",
            mutation_weld_v2(
                "occurrence/mutation",
                "use/mutation/base-a",
                "contact/base-a/alternate",
                PhysicalFeatureUsePolicyV2::Reusable,
                false,
                planned_v2(0x80, 0x83),
            ),
        ),
        (
            "unequal-procedure-path-role-swap",
            mutation_weld_v2(
                "occurrence/mutation",
                "use/mutation/base-a",
                "contact/base-a/alternate",
                PhysicalFeatureUsePolicyV2::Reusable,
                false,
                planned_coordinates_v2(
                    "assembly-v2/path",
                    1,
                    0x81,
                    "assembly-v2/procedure",
                    1,
                    0x80,
                ),
            ),
        ),
        (
            "lifecycle-discriminant-and-evidence",
            mutation_weld_v2(
                "occurrence/mutation",
                "use/mutation/base-a",
                "contact/base-a/alternate",
                PhysicalFeatureUsePolicyV2::Reusable,
                false,
                execution_claimed_v2(0x80, 0x81, 0x84),
            ),
        ),
    ];

    for (field, occurrence) in mutations {
        let mutated = admit_mutation_weld_v2(&graph, occurrence);
        assert_ne!(
            baseline, mutated,
            "isolated {field} mutation must move aggregate identity"
        );
    }

    let claimed_baseline = admit_mutation_weld_v2(
        &graph,
        mutation_weld_v2(
            "occurrence/mutation",
            "use/mutation/base-a",
            "contact/base-a/alternate",
            PhysicalFeatureUsePolicyV2::Reusable,
            false,
            execution_claimed_v2(0x80, 0x81, 0x84),
        ),
    );
    let claimed_evidence_mutations = [
        (
            "namespace",
            execution_claimed_coordinates_v2(
                "assembly-v2/procedure",
                1,
                0x80,
                "assembly-v2/path",
                1,
                0x81,
                "assembly-v2/execution-evidence-changed",
                1,
                0x84,
            ),
        ),
        (
            "schema-version",
            execution_claimed_coordinates_v2(
                "assembly-v2/procedure",
                1,
                0x80,
                "assembly-v2/path",
                1,
                0x81,
                "assembly-v2/execution-evidence",
                2,
                0x84,
            ),
        ),
        ("content-hash", execution_claimed_v2(0x80, 0x81, 0x85)),
    ];
    for (coordinate, lifecycle) in claimed_evidence_mutations {
        let mutated = admit_mutation_weld_v2(
            &graph,
            mutation_weld_v2(
                "occurrence/mutation",
                "use/mutation/base-a",
                "contact/base-a/alternate",
                PhysicalFeatureUsePolicyV2::Reusable,
                false,
                lifecycle,
            ),
        );
        assert_ne!(
            claimed_baseline, mutated,
            "isolated execution-evidence {coordinate} mutation must move identity with the lifecycle tag fixed"
        );
    }

    let graph_mutated = admitted_graph_v2(0x4a);
    let graph_mutated_id = admit_mutation_weld_v2(
        &graph_mutated,
        mutation_weld_v2(
            "occurrence/mutation",
            "use/mutation/base-a",
            "contact/base-a/alternate",
            PhysicalFeatureUsePolicyV2::Reusable,
            false,
            planned_v2(0x80, 0x81),
        ),
    );
    assert_ne!(
        baseline, graph_mutated_id,
        "exact admitted Machine graph identity must move assembly identity"
    );
}

#[test]
fn mas2_010_preload_stack_and_availability_transition_fields_move_identity() {
    let graph = admitted_graph_v2(0x4b);
    let baseline = singleton_draft_v2(
        bolt_occurrence_v2(),
        vec![body_v2("body/base-a"), body_v2("body/base-b")],
        vec![
            body_v2("body/bolt"),
            body_v2("body/nut"),
            body_v2("body/washer"),
        ],
    )
    .admit_against(&graph)
    .expect("baseline bolt transition must admit");

    let mut source_unit_occurrence = bolt_occurrence_v2();
    let JointTopologyV2::PreloadedBolt {
        clamped_members,
        fastener_stack,
        ..
    } = source_unit_occurrence.topology().clone()
    else {
        unreachable!("bolt fixture topology")
    };
    source_unit_occurrence = JointOccurrenceV2::new(
        source_unit_occurrence.id().clone(),
        JointTopologyV2::PreloadedBolt {
            clamped_members,
            fastener_stack,
            preload: preload_v2(2_000.0, AssemblyPreloadUnitV2::Newton),
        },
        source_unit_occurrence.lifecycle().clone(),
    );
    let source_unit = singleton_draft_v2(
        source_unit_occurrence,
        vec![body_v2("body/base-a"), body_v2("body/base-b")],
        vec![
            body_v2("body/bolt"),
            body_v2("body/nut"),
            body_v2("body/washer"),
        ],
    )
    .admit_against(&graph)
    .expect("equal-SI alternative source-unit declaration must admit");
    assert_ne!(
        baseline.identity(),
        source_unit.identity(),
        "submitted preload bits and unit must move identity even when SI force is equal"
    );

    let mut stack_occurrence = bolt_occurrence_v2();
    let JointTopologyV2::PreloadedBolt {
        clamped_members,
        mut fastener_stack,
        preload,
    } = stack_occurrence.topology().clone()
    else {
        unreachable!("bolt fixture topology")
    };
    let washer_index = fastener_stack
        .iter()
        .position(|participant| participant.role() == BoltStackRoleV2::Washer)
        .expect("fixture washer");
    let washer_use = fastener_stack[washer_index].feature_use().clone();
    fastener_stack[washer_index] =
        BoltStackParticipantV2::new(1, BoltStackRoleV2::Spacer, washer_use);
    stack_occurrence = JointOccurrenceV2::new(
        stack_occurrence.id().clone(),
        JointTopologyV2::PreloadedBolt {
            clamped_members,
            fastener_stack,
            preload,
        },
        stack_occurrence.lifecycle().clone(),
    );
    let stack_role = singleton_draft_v2(
        stack_occurrence,
        vec![body_v2("body/base-a"), body_v2("body/base-b")],
        vec![
            body_v2("body/bolt"),
            body_v2("body/nut"),
            body_v2("body/washer"),
        ],
    )
    .admit_against(&graph)
    .expect("typed stack-role mutation remains admissible");
    assert_ne!(
        baseline.identity(),
        stack_role.identity(),
        "typed fastener-stack role must move identity"
    );

    let stack_position_occurrence = bolt_occurrence_v2();
    let JointTopologyV2::PreloadedBolt {
        clamped_members,
        fastener_stack,
        preload,
    } = stack_position_occurrence.topology().clone()
    else {
        unreachable!("bolt fixture topology")
    };
    let fastener_stack = fastener_stack
        .into_iter()
        .map(|participant| {
            let position = match participant.position() {
                0 => 1,
                1 => 0,
                position => position,
            };
            BoltStackParticipantV2::new(
                position,
                participant.role(),
                participant.feature_use().clone(),
            )
        })
        .collect();
    let stack_position_occurrence = JointOccurrenceV2::new(
        stack_position_occurrence.id().clone(),
        JointTopologyV2::PreloadedBolt {
            clamped_members,
            fastener_stack,
            preload,
        },
        stack_position_occurrence.lifecycle().clone(),
    );
    let stack_position = singleton_draft_v2(
        stack_position_occurrence,
        vec![body_v2("body/base-a"), body_v2("body/base-b")],
        vec![
            body_v2("body/bolt"),
            body_v2("body/nut"),
            body_v2("body/washer"),
        ],
    )
    .admit_against(&graph)
    .expect("contiguous stack-position reassignment remains admissible");
    assert_ne!(
        baseline.identity(),
        stack_position.identity(),
        "physical fastener-stack position must move identity independently of role and feature use"
    );

    let key_introduced = singleton_draft_v2(
        key_occurrence_v2(),
        vec![body_v2("body/shaft"), body_v2("body/hub")],
        vec![body_v2("body/key")],
    )
    .admit_against(&graph)
    .expect("key introduction transition must admit");
    let key_initial = singleton_draft_v2(
        key_occurrence_v2(),
        vec![
            body_v2("body/shaft"),
            body_v2("body/hub"),
            body_v2("body/key"),
        ],
        Vec::new(),
    )
    .admit_against(&graph)
    .expect("already-available key topology must admit");
    assert_ne!(
        key_introduced.identity(),
        key_initial.identity(),
        "initial availability plus explicit before/introduction/after transition must be semantic"
    );
}

#[test]
fn mas2_011_step_identity_ordinal_and_occurrence_schedule_are_semantic() {
    let graph = admitted_graph_v2(0x4c);
    let baseline = valid_draft_v2()
        .admit_against(&graph)
        .expect("baseline chronology must admit");

    let mut step_id_mutation = valid_draft_v2();
    step_id_mutation.steps = step_id_mutation
        .steps
        .into_iter()
        .map(|step| {
            let id = if step.id() == &step_id_v2("step/hybrid") {
                step_id_v2("step/hybrid-renamed")
            } else {
                step.id().clone()
            };
            AssemblyStepV2::new(
                id,
                step.ordinal(),
                step.introduced_bodies().to_vec(),
                step.occurrence_ids().to_vec(),
            )
        })
        .collect();
    let step_id_mutation = step_id_mutation
        .admit_against(&graph)
        .expect("step identity mutation must remain admissible");
    assert_ne!(
        baseline.identity(),
        step_id_mutation.identity(),
        "chronological step ID must move identity"
    );

    let mut ordinal_mutation = valid_draft_v2();
    ordinal_mutation.steps = ordinal_mutation
        .steps
        .into_iter()
        .map(|step| {
            let ordinal = match step.id().canonical_key() {
                "step/hybrid" => 2,
                "step/key" => 1,
                _ => step.ordinal(),
            };
            AssemblyStepV2::new(
                step.id().clone(),
                ordinal,
                step.introduced_bodies().to_vec(),
                step.occurrence_ids().to_vec(),
            )
        })
        .collect();
    let ordinal_mutation = ordinal_mutation
        .admit_against(&graph)
        .expect("independent key/hybrid step ordinal swap must remain admissible");
    assert_ne!(
        baseline.identity(),
        ordinal_mutation.identity(),
        "contiguous ordinal assignment and derived before/after rows must move identity"
    );

    let mut schedule_mutation = valid_draft_v2();
    schedule_mutation.steps = schedule_mutation
        .steps
        .into_iter()
        .map(|step| {
            let occurrence_ids = match step.id().canonical_key() {
                "step/hybrid" => vec![
                    occurrence_id_v2("occurrence/weld"),
                    occurrence_id_v2("occurrence/interference"),
                ],
                "step/directed" => vec![
                    occurrence_id_v2("occurrence/spline"),
                    occurrence_id_v2("occurrence/adhesive"),
                ],
                _ => step.occurrence_ids().to_vec(),
            };
            AssemblyStepV2::new(
                step.id().clone(),
                step.ordinal(),
                step.introduced_bodies().to_vec(),
                occurrence_ids,
            )
        })
        .collect();
    let schedule_mutation = schedule_mutation
        .admit_against(&graph)
        .expect("schedule reassignment over already available bodies must admit");
    assert_ne!(
        baseline.identity(),
        schedule_mutation.identity(),
        "step-to-occurrence schedule links must move identity"
    );
}

#[test]
fn mas2_012_atomic_chronology_refuses_unavailable_unused_and_duplicate_links() {
    let graph = admitted_graph_v2(0x4d);

    let unavailable = singleton_draft_v2(
        key_occurrence_v2(),
        vec![body_v2("body/shaft"), body_v2("body/hub")],
        Vec::new(),
    );
    let error = unavailable
        .admit_against(&graph)
        .expect_err("key participant unavailable at the atomic step must refuse");
    assert!(
        matches!(
            error,
            MachineAssemblyAdmissionErrorV2::ParticipantBodyUnavailable {
                ref step,
                ref occurrence,
                ref body,
                ..
            } if *step == step_id_v2("step/single")
                && *occurrence == occurrence_id_v2("occurrence/key")
                && *body == body_v2("body/key")
        ),
        "unavailable-body refusal must identify step, occurrence, and body; got {error:?}"
    );

    let unused = singleton_draft_v2(
        weld_occurrence_v2(),
        vec![body_v2("body/base-a"), body_v2("body/base-b")],
        vec![body_v2("body/spare")],
    );
    assert_eq!(
        unused.admit_against(&graph),
        Err(
            MachineAssemblyAdmissionErrorV2::IntroducedBodyDoesNotParticipate {
                step: step_id_v2("step/single"),
                body: body_v2("body/spare"),
            }
        ),
        "an introduction must not publish unless the new body participates in the complete step"
    );

    let occurrence = weld_occurrence_v2();
    let occurrence_id = occurrence.id().clone();
    let duplicate_reference = MachineAssemblyDraftV2 {
        initial_available_bodies: vec![body_v2("body/base-a"), body_v2("body/base-b")],
        steps: vec![AssemblyStepV2::new(
            step_id_v2("step/duplicate-reference"),
            0,
            Vec::new(),
            vec![occurrence_id.clone(), occurrence_id.clone()],
        )],
        occurrences: vec![occurrence],
    };
    assert_eq!(
        duplicate_reference.admit_against(&graph),
        Err(
            MachineAssemblyAdmissionErrorV2::DuplicateOccurrenceReference {
                step: step_id_v2("step/duplicate-reference"),
                occurrence: occurrence_id,
            }
        ),
        "one step cannot schedule the same physical occurrence twice"
    );

    let occurrence = weld_occurrence_v2();
    let occurrence_id = occurrence.id().clone();
    let scheduled_twice = MachineAssemblyDraftV2 {
        initial_available_bodies: vec![body_v2("body/base-a"), body_v2("body/base-b")],
        steps: vec![
            AssemblyStepV2::new(
                step_id_v2("step/first-schedule"),
                0,
                Vec::new(),
                vec![occurrence_id.clone()],
            ),
            AssemblyStepV2::new(
                step_id_v2("step/second-schedule"),
                1,
                Vec::new(),
                vec![occurrence_id.clone()],
            ),
        ],
        occurrences: vec![occurrence],
    };
    assert_eq!(
        scheduled_twice.admit_against(&graph),
        Err(MachineAssemblyAdmissionErrorV2::OccurrenceScheduledTwice {
            occurrence: occurrence_id,
            first: step_id_v2("step/first-schedule"),
            duplicate: step_id_v2("step/second-schedule"),
        }),
        "a physical occurrence belongs to exactly one chronological step"
    );
}

fn growing_schedule_fixture_v2(
    step_count: usize,
) -> (AdmittedMachineGraph, MachineAssemblyDraftV2) {
    let anchor_body = body_v2("body/growing-anchor");
    let anchor_feature = feature_v2("contact/growing-anchor");
    let introduced = (0..step_count)
        .map(|index| {
            (
                body_v2(&format!("body/growing-{index:04}")),
                feature_v2(&format!("contact/growing-{index:04}")),
            )
        })
        .collect::<Vec<_>>();
    let mut bodies = vec![anchor_body.clone()];
    bodies.extend(introduced.iter().map(|(body, _)| body.clone()));
    let mut features = vec![anchor_feature];
    features.extend(introduced.iter().map(|(_, feature)| feature.clone()));
    let materials = bodies
        .iter()
        .cloned()
        .enumerate()
        .map(|(index, body)| {
            material_v2(
                body,
                &format!("materials/growing-{index:04}"),
                u8::try_from(index % 254 + 1).expect("fixture material byte fits"),
            )
        })
        .collect();
    let graph = MachineGraphDraft {
        clocks: Vec::new(),
        subsystems: vec![SubsystemSpec {
            id: SubsystemId::new("subsystem/growing").expect("canonical subsystem"),
            model: ModelRef::new("models/growing", nz_v2(1), [0xb1; 32]).expect("canonical model"),
            bodies,
            surface_patches: Vec::new(),
            contact_features: features,
            state_slots: Vec::new(),
        }],
        terminals: Vec::new(),
        ports: Vec::new(),
        relations: Vec::new(),
        materials,
        interfaces: Vec::new(),
    }
    .admit()
    .expect("growing availability fixture graph must admit");

    let occurrences = introduced
        .iter()
        .enumerate()
        .map(|(index, (body, feature))| {
            JointOccurrenceV2::new(
                occurrence_id_v2(&format!("occurrence/growing-{index:04}")),
                JointTopologyV2::Weld {
                    members: vec![
                        JointFeatureUseV2::new(
                            use_id_v2(&format!("use/growing-{index:04}/anchor")),
                            fs_ir::machine::manufacturing::assembly::AssemblyFeatureSelectorV2::new(
                                anchor_body.clone(),
                                feature_v2("contact/growing-anchor"),
                            ),
                            PhysicalFeatureUsePolicyV2::Reusable,
                        ),
                        JointFeatureUseV2::new(
                            use_id_v2(&format!("use/growing-{index:04}/introduced")),
                            fs_ir::machine::manufacturing::assembly::AssemblyFeatureSelectorV2::new(
                                body.clone(),
                                feature.clone(),
                            ),
                            PhysicalFeatureUsePolicyV2::Reusable,
                        ),
                    ],
                },
                planned_v2(0xb2, 0xb3),
            )
        })
        .collect::<Vec<_>>();
    let steps = occurrences
        .iter()
        .enumerate()
        .map(|(index, occurrence)| {
            AssemblyStepV2::new(
                step_id_v2(&format!("step/growing-{index:04}")),
                u32::try_from(index).expect("fixture ordinal fits u32"),
                vec![introduced[index].0.clone()],
                vec![occurrence.id().clone()],
            )
        })
        .collect();
    (
        graph,
        MachineAssemblyDraftV2 {
            initial_available_bodies: vec![anchor_body],
            steps,
            occurrences,
        },
    )
}

#[test]
fn mas2_013_availability_transition_chain_is_replayable_history_bound_and_linear() {
    let (small_graph, small_draft) = growing_schedule_fixture_v2(16);
    let small_oracle = oracle_receipt_v2(
        small_graph.identity().as_bytes(),
        &small_draft,
        MACHINE_ASSEMBLY_SCHEMA_VERSION_V2,
        IR_VERSION,
    );
    let small = small_draft
        .admit_against(&small_graph)
        .expect("16-step growing schedule must admit");
    let (large_graph, large_draft) = growing_schedule_fixture_v2(32);
    let large_oracle = oracle_receipt_v2(
        large_graph.identity().as_bytes(),
        &large_draft,
        MACHINE_ASSEMBLY_SCHEMA_VERSION_V2,
        IR_VERSION,
    );
    let large = large_draft
        .admit_against(&large_graph)
        .expect("32-step growing schedule must admit");

    assert_eq!(small_oracle.availability_rows_hashed, 17);
    assert_eq!(large_oracle.availability_rows_hashed, 33);
    assert_eq!(
        small_oracle.initial_availability_preimage_bytes,
        large_oracle.initial_availability_preimage_bytes,
        "the common one-body initial set is hashed exactly once"
    );
    assert_eq!(
        large_oracle.transition_preimage_bytes,
        2 * small_oracle.transition_preimage_bytes,
        "doubling fixed-width one-body transitions must exactly double transition-hash work"
    );
    let root_width = core::mem::size_of::<ContentId>();
    assert_eq!(root_width, 32, "a retained chain root is one fixed digest");
    let small_retained_root_bytes = small.steps().len() * 2 * root_width;
    let large_retained_root_bytes = large.steps().len() * 2 * root_width;
    assert_eq!(small_retained_root_bytes, 1_024);
    assert_eq!(large_retained_root_bytes, 2_048);
    assert_eq!(
        large_retained_root_bytes,
        2 * small_retained_root_bytes,
        "the two retained fixed-size roots per step must scale linearly"
    );
    for (admitted, (before, after)) in small.steps().iter().zip(&small_oracle.availability_roots) {
        assert_eq!(admitted.availability_before_root(), *before);
        assert_eq!(admitted.availability_after_root(), *after);
    }
    for (admitted, (before, after)) in large.steps().iter().zip(&large_oracle.availability_roots) {
        assert_eq!(admitted.availability_before_root(), *before);
        assert_eq!(admitted.availability_after_root(), *after);
    }
    assert_eq!(
        small.steps()[15].availability_after_root(),
        large.steps()[15].availability_after_root(),
        "an independent verifier needs only the common initial set and first 16 transitions"
    );
    for pair in large.steps().windows(2) {
        assert_eq!(
            pair[0].availability_after_root(),
            pair[1].availability_before_root(),
            "each before root must be exactly the preceding after root"
        );
    }

    let (mutation_graph, baseline_draft) = growing_schedule_fixture_v2(3);
    let baseline = baseline_draft
        .clone()
        .admit_against(&mutation_graph)
        .expect("baseline three-step schedule must admit");
    let mut reordered = baseline_draft.clone();
    reordered.steps = reordered
        .steps
        .into_iter()
        .map(|step| {
            let ordinal = match step.ordinal() {
                0 => 1,
                1 => 0,
                value => value,
            };
            AssemblyStepV2::new(
                step.id().clone(),
                ordinal,
                step.introduced_bodies().to_vec(),
                step.occurrence_ids().to_vec(),
            )
        })
        .collect();
    let reordered = reordered
        .admit_against(&mutation_graph)
        .expect("reordered valid introduction history must admit");
    assert_ne!(
        baseline.steps()[0].availability_after_root(),
        reordered.steps()[0].availability_after_root(),
        "changing the first introduced body/order must change the first transition root"
    );
    assert_ne!(
        baseline.steps()[2].availability_after_root(),
        reordered.steps()[2].availability_after_root(),
        "equal final sets reached through different histories must retain different chain roots"
    );

    let mut renamed = baseline_draft;
    renamed.steps = renamed
        .steps
        .into_iter()
        .map(|step| {
            let id = if step.ordinal() == 0 {
                step_id_v2("step/growing-renamed")
            } else {
                step.id().clone()
            };
            AssemblyStepV2::new(
                id,
                step.ordinal(),
                step.introduced_bodies().to_vec(),
                step.occurrence_ids().to_vec(),
            )
        })
        .collect();
    let renamed = renamed
        .admit_against(&mutation_graph)
        .expect("renamed valid step must admit");
    assert_ne!(
        baseline.steps()[0].availability_after_root(),
        renamed.steps()[0].availability_after_root(),
        "moving a transition under a different stable step identity must change its root"
    );
}

fn boundary_occurrence_v2(index: usize) -> JointOccurrenceV2 {
    JointOccurrenceV2::new(
        occurrence_id_v2(&format!("occurrence/boundary-{index:04}")),
        JointTopologyV2::Spline {
            external: feature_use_v2(
                &format!("use/boundary-{index:04}/external"),
                "body/external",
                "contact/external/main",
            ),
            internal: feature_use_v2(
                &format!("use/boundary-{index:04}/internal"),
                "body/internal",
                "contact/internal/main",
            ),
        },
        planned_v2(0x90, 0x91),
    )
}

fn boundary_draft_v2() -> MachineAssemblyDraftV2 {
    let occurrences = (0..MAX_MACHINE_ASSEMBLY_OCCURRENCES_V2)
        .map(boundary_occurrence_v2)
        .collect::<Vec<_>>();
    let steps = occurrences
        .iter()
        .enumerate()
        .map(|(index, occurrence)| {
            AssemblyStepV2::new(
                step_id_v2(&format!("step/boundary-{index:04}")),
                u32::try_from(index).expect("step boundary fits u32"),
                Vec::new(),
                vec![occurrence.id().clone()],
            )
        })
        .collect();
    MachineAssemblyDraftV2 {
        initial_available_bodies: vec![body_v2("body/external"), body_v2("body/internal")],
        steps,
        occurrences,
    }
}

#[test]
fn mas2_014_exact_step_occurrence_and_participant_caps_preflight_n_plus_one() {
    let graph = admitted_graph_v2(0x4e);
    let exact = boundary_draft_v2();
    let admitted = exact
        .clone()
        .admit_against(&graph)
        .expect("exact 4,096-step/occurrence boundary must admit");
    assert_eq!(
        admitted.steps().len(),
        MAX_MACHINE_ASSEMBLY_STEPS_V2,
        "exact chronological step cap must be retained"
    );
    assert_eq!(
        admitted.occurrences().len(),
        MAX_MACHINE_ASSEMBLY_OCCURRENCES_V2,
        "exact physical occurrence cap must be retained"
    );

    let mut occurrence_over = exact.clone();
    occurrence_over.occurrences.push(boundary_occurrence_v2(0));
    assert_eq!(
        occurrence_over.admit_against(&graph),
        Err(MachineAssemblyAdmissionErrorV2::OccurrenceLimit {
            actual: MAX_MACHINE_ASSEMBLY_OCCURRENCES_V2 + 1,
            max: MAX_MACHINE_ASSEMBLY_OCCURRENCES_V2,
        }),
        "raw occurrence N+1 must refuse before duplicate-ID or nested processing"
    );

    let mut step_over = exact;
    step_over.steps.push(AssemblyStepV2::new(
        step_id_v2("step/boundary-over"),
        u32::try_from(MAX_MACHINE_ASSEMBLY_STEPS_V2).expect("step boundary fits u32"),
        Vec::new(),
        vec![occurrence_id_v2("occurrence/boundary-0000")],
    ));
    assert_eq!(
        step_over.admit_against(&graph),
        Err(MachineAssemblyAdmissionErrorV2::StepLimit {
            actual: MAX_MACHINE_ASSEMBLY_STEPS_V2 + 1,
            max: MAX_MACHINE_ASSEMBLY_STEPS_V2,
        }),
        "raw chronological step N+1 must refuse before schedule validation"
    );

    let too_many_members = (0..=MAX_MACHINE_ASSEMBLY_PARTICIPANTS_PER_OCCURRENCE_V2)
        .map(|index| {
            feature_use_v2(
                &format!("use/participant-over-{index}"),
                if index % 2 == 0 {
                    "body/base-a"
                } else {
                    "body/base-b"
                },
                if index % 2 == 0 {
                    "contact/base-a/main"
                } else {
                    "contact/base-b/main"
                },
            )
        })
        .collect();
    let participant_over = JointOccurrenceV2::new(
        occurrence_id_v2("occurrence/participant-over"),
        JointTopologyV2::Weld {
            members: too_many_members,
        },
        planned_v2(0x92, 0x93),
    );
    let error = singleton_draft_v2(
        participant_over,
        vec![body_v2("body/base-a"), body_v2("body/base-b")],
        Vec::new(),
    )
    .admit_against(&graph)
    .expect_err("raw participant N+1 must refuse before sorting or graph traversal");
    assert_eq!(
        error,
        MachineAssemblyAdmissionErrorV2::InvalidTopology {
            occurrence: occurrence_id_v2("occurrence/participant-over"),
            issue: AssemblyTopologyIssueV2::ParticipantLimit {
                actual: MAX_MACHINE_ASSEMBLY_PARTICIPANTS_PER_OCCURRENCE_V2 + 1,
                max: MAX_MACHINE_ASSEMBLY_PARTICIPANTS_PER_OCCURRENCE_V2,
            },
        }
    );
}

fn initial_body_boundary_fixture_v2() -> (AdmittedMachineGraph, MachineAssemblyDraftV2) {
    let bodies = (0..MAX_MACHINE_ASSEMBLY_INITIAL_BODIES_V2)
        .map(|index| body_v2(&format!("body/initial-boundary-{index:04}")))
        .collect::<Vec<_>>();
    let features = (0..2)
        .map(|index| feature_v2(&format!("contact/initial-boundary-{index:04}")))
        .collect::<Vec<_>>();
    let materials = bodies
        .iter()
        .cloned()
        .enumerate()
        .map(|(index, body)| {
            material_v2(
                body,
                &format!("materials/initial-boundary-{index:04}"),
                u8::try_from(index % 254 + 1).expect("fixture material byte fits"),
            )
        })
        .collect();
    let graph = MachineGraphDraft {
        clocks: Vec::new(),
        subsystems: vec![SubsystemSpec {
            id: SubsystemId::new("subsystem/initial-boundary").expect("canonical subsystem"),
            model: ModelRef::new("models/initial-boundary", nz_v2(1), [0xc1; 32])
                .expect("canonical model"),
            bodies: bodies.clone(),
            surface_patches: Vec::new(),
            contact_features: features.clone(),
            state_slots: Vec::new(),
        }],
        terminals: Vec::new(),
        ports: Vec::new(),
        relations: Vec::new(),
        materials,
        interfaces: Vec::new(),
    }
    .admit()
    .expect("initial-body boundary graph must admit");
    let occurrence = JointOccurrenceV2::new(
        occurrence_id_v2("occurrence/initial-boundary"),
        JointTopologyV2::Weld {
            members: vec![
                JointFeatureUseV2::new(
                    use_id_v2("use/initial-boundary-0000"),
                    fs_ir::machine::manufacturing::assembly::AssemblyFeatureSelectorV2::new(
                        bodies[0].clone(),
                        features[0].clone(),
                    ),
                    PhysicalFeatureUsePolicyV2::Reusable,
                ),
                JointFeatureUseV2::new(
                    use_id_v2("use/initial-boundary-0001"),
                    fs_ir::machine::manufacturing::assembly::AssemblyFeatureSelectorV2::new(
                        bodies[1].clone(),
                        features[1].clone(),
                    ),
                    PhysicalFeatureUsePolicyV2::Reusable,
                ),
            ],
        },
        planned_v2(0xc2, 0xc3),
    );
    let occurrence_id = occurrence.id().clone();
    (
        graph,
        MachineAssemblyDraftV2 {
            initial_available_bodies: bodies,
            steps: vec![AssemblyStepV2::new(
                step_id_v2("step/initial-boundary"),
                0,
                Vec::new(),
                vec![occurrence_id],
            )],
            occurrences: vec![occurrence],
        },
    )
}

fn introduction_boundary_fixture_v2() -> (AdmittedMachineGraph, MachineAssemblyDraftV2) {
    let anchor = body_v2("body/introduction-anchor");
    let introduced = (0..MAX_MACHINE_ASSEMBLY_INTRODUCTIONS_PER_STEP_V2)
        .map(|index| {
            (
                body_v2(&format!("body/introduction-{index:04}")),
                feature_v2(&format!("contact/introduction-{index:04}")),
            )
        })
        .collect::<Vec<_>>();
    let mut bodies = vec![anchor.clone()];
    bodies.extend(introduced.iter().map(|(body, _)| body.clone()));
    let materials = bodies
        .iter()
        .cloned()
        .enumerate()
        .map(|(index, body)| {
            material_v2(
                body,
                &format!("materials/introduction-{index:04}"),
                u8::try_from(index % 254 + 1).expect("fixture material byte fits"),
            )
        })
        .collect();
    let graph = MachineGraphDraft {
        clocks: Vec::new(),
        subsystems: vec![SubsystemSpec {
            id: SubsystemId::new("subsystem/introduction-boundary").expect("canonical subsystem"),
            model: ModelRef::new("models/introduction-boundary", nz_v2(1), [0xc4; 32])
                .expect("canonical model"),
            bodies,
            surface_patches: Vec::new(),
            contact_features: introduced
                .iter()
                .map(|(_, feature)| feature.clone())
                .collect(),
            state_slots: Vec::new(),
        }],
        terminals: Vec::new(),
        ports: Vec::new(),
        relations: Vec::new(),
        materials,
        interfaces: Vec::new(),
    }
    .admit()
    .expect("introduction boundary graph must admit");
    let occurrence = JointOccurrenceV2::new(
        occurrence_id_v2("occurrence/introduction-boundary"),
        JointTopologyV2::Weld {
            members: introduced
                .iter()
                .enumerate()
                .map(|(index, (body, feature))| {
                    JointFeatureUseV2::new(
                        use_id_v2(&format!("use/introduction-{index:04}")),
                        fs_ir::machine::manufacturing::assembly::AssemblyFeatureSelectorV2::new(
                            body.clone(),
                            feature.clone(),
                        ),
                        PhysicalFeatureUsePolicyV2::Reusable,
                    )
                })
                .collect(),
        },
        planned_v2(0xc5, 0xc6),
    );
    let occurrence_id = occurrence.id().clone();
    (
        graph,
        MachineAssemblyDraftV2 {
            initial_available_bodies: vec![anchor],
            steps: vec![AssemblyStepV2::new(
                step_id_v2("step/introduction-boundary"),
                0,
                introduced.iter().map(|(body, _)| body.clone()).collect(),
                vec![occurrence_id],
            )],
            occurrences: vec![occurrence],
        },
    )
}

#[test]
fn mas2_015_initial_introduction_and_step_reference_caps_have_exact_n_plus_one_evidence() {
    let (initial_graph, initial_exact) = initial_body_boundary_fixture_v2();
    let admitted = initial_exact
        .clone()
        .admit_against(&initial_graph)
        .expect("exact 4,096 initial-body boundary must admit");
    assert_eq!(
        admitted.initial_available_bodies().len(),
        MAX_MACHINE_ASSEMBLY_INITIAL_BODIES_V2
    );
    let mut initial_over = initial_exact;
    initial_over
        .initial_available_bodies
        .push(body_v2("body/initial-boundary-over"));
    assert_eq!(
        initial_over.admit_against(&initial_graph),
        Err(MachineAssemblyAdmissionErrorV2::InitialBodyLimit {
            actual: MAX_MACHINE_ASSEMBLY_INITIAL_BODIES_V2 + 1,
            max: MAX_MACHINE_ASSEMBLY_INITIAL_BODIES_V2,
        }),
        "raw initial-body N+1 must refuse before graph ownership or duplicate analysis"
    );

    let (introduction_graph, introduction_exact) = introduction_boundary_fixture_v2();
    let admitted = introduction_exact
        .clone()
        .admit_against(&introduction_graph)
        .expect("exact 64-introduction step boundary must admit");
    assert_eq!(
        admitted.steps()[0].step().introduced_bodies().len(),
        MAX_MACHINE_ASSEMBLY_INTRODUCTIONS_PER_STEP_V2
    );
    let mut introduction_over = introduction_exact;
    let step = introduction_over.steps.pop().expect("one boundary step");
    let mut introduced_bodies = step.introduced_bodies().to_vec();
    introduced_bodies.push(body_v2("body/introduction-over"));
    introduction_over.steps.push(AssemblyStepV2::new(
        step.id().clone(),
        step.ordinal(),
        introduced_bodies,
        step.occurrence_ids().to_vec(),
    ));
    assert_eq!(
        introduction_over.admit_against(&introduction_graph),
        Err(MachineAssemblyAdmissionErrorV2::StepIntroductionLimit {
            step: step_id_v2("step/introduction-boundary"),
            actual: MAX_MACHINE_ASSEMBLY_INTRODUCTIONS_PER_STEP_V2 + 1,
            max: MAX_MACHINE_ASSEMBLY_INTRODUCTIONS_PER_STEP_V2,
        }),
        "raw per-step introduction N+1 must refuse before graph traversal"
    );

    let occurrence_graph = admitted_graph_v2(0xc7);
    let occurrences = (0..MAX_MACHINE_ASSEMBLY_OCCURRENCES_PER_STEP_V2)
        .map(boundary_occurrence_v2)
        .collect::<Vec<_>>();
    let occurrence_ids = occurrences
        .iter()
        .map(|occurrence| occurrence.id().clone())
        .collect::<Vec<_>>();
    let exact_refs = MachineAssemblyDraftV2 {
        initial_available_bodies: vec![body_v2("body/external"), body_v2("body/internal")],
        steps: vec![AssemblyStepV2::new(
            step_id_v2("step/reference-boundary"),
            0,
            Vec::new(),
            occurrence_ids,
        )],
        occurrences,
    };
    let admitted = exact_refs
        .clone()
        .admit_against(&occurrence_graph)
        .expect("exact 64-occurrence-reference step boundary must admit");
    assert_eq!(
        admitted.steps()[0].step().occurrence_ids().len(),
        MAX_MACHINE_ASSEMBLY_OCCURRENCES_PER_STEP_V2
    );
    let mut refs_over = exact_refs;
    let step = refs_over.steps.pop().expect("one reference boundary step");
    let mut occurrence_ids = step.occurrence_ids().to_vec();
    occurrence_ids.push(occurrence_id_v2("occurrence/reference-over"));
    refs_over.steps.push(AssemblyStepV2::new(
        step.id().clone(),
        step.ordinal(),
        step.introduced_bodies().to_vec(),
        occurrence_ids,
    ));
    assert_eq!(
        refs_over.admit_against(&occurrence_graph),
        Err(MachineAssemblyAdmissionErrorV2::StepOccurrenceLimit {
            step: step_id_v2("step/reference-boundary"),
            actual: MAX_MACHINE_ASSEMBLY_OCCURRENCES_PER_STEP_V2 + 1,
            max: MAX_MACHINE_ASSEMBLY_OCCURRENCES_PER_STEP_V2,
        }),
        "raw per-step occurrence-reference N+1 must refuse before lookup or duplicate analysis"
    );
}

fn maximum_key_v2(prefix: &str, index: usize) -> String {
    let mut key = format!("{prefix}{index}a");
    assert!(key.len() <= 128, "maximum-key prefix must fit");
    key.extend(std::iter::repeat_n('a', 128 - key.len()));
    assert_eq!(key.len(), 128, "maximum grammar key must be 128 bytes");
    key
}

fn maximum_width_graph_v2() -> (AdmittedMachineGraph, Vec<(String, String)>) {
    let pairs = (0..MAX_MACHINE_ASSEMBLY_PARTICIPANTS_PER_OCCURRENCE_V2)
        .map(|index| {
            (
                maximum_key_v2("body", index),
                maximum_key_v2("feature", index),
            )
        })
        .collect::<Vec<_>>();
    let bodies = pairs
        .iter()
        .map(|(body_key, _)| body_v2(body_key))
        .collect::<Vec<_>>();
    let materials = bodies
        .iter()
        .cloned()
        .enumerate()
        .map(|(index, body)| {
            material_v2(
                body,
                &format!("materials/maximum-{index}"),
                u8::try_from(index + 1).expect("maximum fixture material byte fits"),
            )
        })
        .collect::<Vec<_>>();
    let graph = MachineGraphDraft {
        clocks: Vec::new(),
        subsystems: vec![SubsystemSpec {
            id: SubsystemId::new("subsystem/maximum-width").expect("canonical subsystem"),
            model: ModelRef::new("models/maximum-width", nz_v2(1), [0xa0; 32])
                .expect("canonical model"),
            bodies,
            surface_patches: Vec::new(),
            contact_features: pairs
                .iter()
                .map(|(_, feature_key)| feature_v2(feature_key))
                .collect(),
            state_slots: Vec::new(),
        }],
        terminals: Vec::new(),
        ports: Vec::new(),
        relations: Vec::new(),
        materials,
        interfaces: Vec::new(),
    }
    .admit()
    .expect("maximum-width graph must admit");
    (graph, pairs)
}

#[test]
fn mas2_016_true_maximum_occurrence_field_and_safe_rectangular_envelope_are_pinned() {
    let (graph, pairs) = maximum_width_graph_v2();
    let clamped_members = pairs[..2]
        .iter()
        .enumerate()
        .map(|(index, (body_key, feature_key))| {
            feature_use_v2(&maximum_key_v2("featureuse", index), body_key, feature_key)
        })
        .collect::<Vec<_>>();
    let fastener_stack = pairs[2..]
        .iter()
        .enumerate()
        .map(|(position, (body_key, feature_key))| {
            let role = match position {
                0 => BoltStackRoleV2::Bolt,
                1 => BoltStackRoleV2::Nut,
                _ => BoltStackRoleV2::Washer,
            };
            BoltStackParticipantV2::new(
                u16::try_from(position).expect("maximum stack position fits u16"),
                role,
                feature_use_v2(
                    &maximum_key_v2("featureuse", position + 2),
                    body_key,
                    feature_key,
                ),
            )
        })
        .collect::<Vec<_>>();
    let occurrence = JointOccurrenceV2::new(
        occurrence_id_v2(&maximum_key_v2("occurrence", 0)),
        JointTopologyV2::PreloadedBolt {
            clamped_members,
            fastener_stack,
            preload: preload_v2(f64::MAX / 2.0, AssemblyPreloadUnitV2::Newton),
        },
        AssemblyLifecycleV2::ExecutionClaimed {
            procedure: AssemblyProcedureRefV2::new(artifact_v2(
                &maximum_key_v2("procedure", 0),
                0xa1,
            )),
            path: AssemblyPathRefV2::new(artifact_v2(&maximum_key_v2("path", 0), 0xa2)),
            evidence: AssemblyExecutionEvidenceRefV2::new(artifact_v2(
                &maximum_key_v2("evidence", 0),
                0xa3,
            )),
        },
    );
    let occurrence_id = occurrence.id().clone();
    let draft = MachineAssemblyDraftV2 {
        initial_available_bodies: pairs
            .iter()
            .map(|(body_key, _)| body_v2(body_key))
            .collect(),
        steps: vec![AssemblyStepV2::new(
            step_id_v2(&maximum_key_v2("step", 0)),
            0,
            Vec::new(),
            vec![occurrence_id],
        )],
        occurrences: vec![occurrence],
    };
    let oracle = oracle_receipt_v2(
        graph.identity().as_bytes(),
        &draft,
        MACHINE_ASSEMBLY_SCHEMA_VERSION_V2,
        IR_VERSION,
    );
    let occurrence_row = oracle_occurrence_row_v2(&draft.occurrences[0]);
    let step_row = oracle_step_row_v2(
        &draft.steps[0],
        draft.initial_available_bodies.len(),
        oracle.availability_roots[0].0,
        draft.initial_available_bodies.len(),
        oracle.availability_roots[0].1,
    );
    let admitted = draft
        .admit_against(&graph)
        .expect("true 64-participant maximum-width row must admit");
    assert_eq!(
        occurrence_row.len(),
        33_741,
        "ExecutionClaimed PreloadedBolt with 64 maximum-width participants is the true maximum occurrence row"
    );
    assert_eq!(
        step_row.len(),
        396,
        "maximum-width IDs plus history-bound before/after roots must pin the exact step-row byte count"
    );
    let computed_max_occurrence_field = 8_u64
        + u64::try_from(MAX_MACHINE_ASSEMBLY_OCCURRENCES_V2).expect("occurrence cap fits u64")
            * (8 + u64::try_from(occurrence_row.len()).expect("row length fits u64"));
    assert_eq!(
        computed_max_occurrence_field, 138_235_912,
        "computed 4,096-row true maximum-width occurrence field must remain reviewable and pinned"
    );
    assert!(
        computed_max_occurrence_field <= MACHINE_ASSEMBLY_IDENTITY_LIMITS_V2.max_field_bytes(),
        "computed true maximum-width occurrence field must fit the declared canonical envelope"
    );
    let computed_independent_max_step_row = 136_u64
        + 4
        + 8
        + u64::try_from(MAX_MACHINE_ASSEMBLY_INTRODUCTIONS_PER_STEP_V2)
            .expect("introduction cap fits u64")
            * (8 + 176)
        + 8
        + u64::try_from(MAX_MACHINE_ASSEMBLY_OCCURRENCES_PER_STEP_V2)
            .expect("per-step occurrence cap fits u64")
            * (8 + 136)
        + 2 * (8 + 8 + 32);
    assert_eq!(
        computed_independent_max_step_row, 21_244,
        "one independently maximized step row includes 64 maximum-key introductions and references plus both count/root pairs"
    );
    let rectangular_step_field_envelope = 8_u64
        + u64::try_from(MAX_MACHINE_ASSEMBLY_STEPS_V2).expect("step cap fits u64")
            * (8 + computed_independent_max_step_row);
    assert_eq!(rectangular_step_field_envelope, 87_048_200);
    assert!(
        rectangular_step_field_envelope <= MACHINE_ASSEMBLY_IDENTITY_LIMITS_V2.max_field_bytes(),
        "the rectangular all-rows-at-independent-maximum step envelope must fit the field limit"
    );
    let computed_independent_max_initial_field = 8_u64
        + u64::try_from(MAX_MACHINE_ASSEMBLY_INITIAL_BODIES_V2).expect("initial-body cap fits u64")
            * (8 + 176);
    assert_eq!(computed_independent_max_initial_field, 753_672);
    // This is a deliberately conservative rectangular resource envelope, not
    // a simultaneously realizable assembly. Exactly-once scheduling limits
    // total step references to at most 4,096 occurrences, so 4,096 steps
    // cannot each carry 64 references; the Machine-graph owned-element cap
    // likewise prevents maximum introductions in every step from coexisting
    // with the maximum initial set.
    let safe_rectangular_independent_maxima_envelope = computed_max_occurrence_field
        + rectangular_step_field_envelope
        + computed_independent_max_initial_field;
    assert_eq!(safe_rectangular_independent_maxima_envelope, 226_037_784);
    assert!(
        safe_rectangular_independent_maxima_envelope
            < MACHINE_ASSEMBLY_IDENTITY_LIMITS_V2.max_canonical_bytes(),
        "the safe rectangular sum of independent collection maxima must leave room for fixed schema framing"
    );
    assert_eq!(
        admitted.identity_receipt().canonical_bytes(),
        oracle.canonical_frame.len() as u64,
        "maximum-width production frame must match independent exact byte accounting"
    );
    assert_eq!(
        admitted.identity_receipt().canonical_preimage(),
        oracle.canonical_preimage,
        "maximum-width production frame must match the independent golden preimage oracle"
    );
    assert_eq!(
        admitted.identity(),
        oracle.identity,
        "maximum-width semantic identity must match the independent raw-frame oracle"
    );
}
