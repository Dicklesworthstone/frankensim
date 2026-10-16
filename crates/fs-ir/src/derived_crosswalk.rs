//! Nominal Machine-IR crosswalk candidates for admitted derived geometry (RD.1b).
//!
//! Machine IR does not yet have an admitted semantic model graph. This module
//! therefore binds one exact, already-admitted [`AdmittedDerivedGeometryV1`]
//! to explicit *selectors* on the Machine-IR side and to nominal mapping
//! artifacts. The returned token is a content-addressed structural candidate,
//! not an executable conversion or a statement that either model is true.
//!
//! In particular, this module exposes no model construction, map execution,
//! inverse, composition, equivalence, physical-validity, semantic-preservation,
//! or evidence-transport API. A later Machine-IR admission lane must replace
//! the nominal selectors with strong identities before any such authority can
//! be considered.

#![allow(clippy::too_many_lines)]

use core::fmt;

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, ChildSpec, EvidenceNodeId,
    Field, FieldSpec, IdentityReceipt, StrongIdentity, WireType,
};
use fs_exec::Cx;
use fs_geom::derived::{
    AdmittedDerivedGeometryV1, DerivedFrameIdV1, DerivedGeometryIdV1, DerivedModelVersionIdV1,
    DerivedSubjectIdV1, DerivedUnitSystemIdV1,
};

use crate::IR_VERSION;

/// Current schema for nominal derived-geometry to Machine-IR crosswalk candidates.
pub const DERIVED_MACHINE_MODEL_CROSSWALK_CANDIDATE_SCHEMA_VERSION_V1: u32 = 1;

const DERIVED_MACHINE_MODEL_CROSSWALK_IDENTITY_LIMITS_V1: CanonicalLimits =
    CanonicalLimits::new(1 << 15, 1 << 14, 17, 128, 4096);
static DERIVED_GEOMETRY_CHILD_V1: ChildSpec = ChildSpec::for_identity::<DerivedGeometryIdV1>();

trait DigestBytes {
    fn digest_bytes(&self) -> &[u8; 32];
}

macro_rules! opaque_crosswalk_id {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name([u8; 32]);

        impl $name {
            /// Construct one explicitly nominal selector or artifact identity.
            /// The bytes alone confer no semantic or scientific authority.
            #[must_use]
            pub const fn from_bytes(bytes: [u8; 32]) -> Self {
                Self(bytes)
            }

            /// Borrow the exact retained bytes.
            #[must_use]
            pub const fn as_bytes(&self) -> &[u8; 32] {
                &self.0
            }
        }

        impl DigestBytes for $name {
            fn digest_bytes(&self) -> &[u8; 32] {
                self.as_bytes()
            }
        }
    };
}

opaque_crosswalk_id!(
    /// Nominal selector for a future admitted Machine-IR model.
    MachineIrModelSelectorIdV1
);
opaque_crosswalk_id!(
    /// Nominal selector for a future admitted Machine-IR subject/entity.
    MachineIrSubjectSelectorIdV1
);
opaque_crosswalk_id!(
    /// Nominal selector for an immutable Machine-IR model version.
    MachineIrModelVersionSelectorIdV1
);
opaque_crosswalk_id!(
    /// Nominal selector for a Machine-IR frame convention.
    MachineIrFrameSelectorIdV1
);
opaque_crosswalk_id!(
    /// Nominal selector for a Machine-IR unit-system convention.
    MachineIrUnitSystemSelectorIdV1
);
opaque_crosswalk_id!(
    /// Nominal artifact relating the exact derived and Machine-IR subjects.
    MachineIrSubjectCrosswalkIdV1
);
opaque_crosswalk_id!(
    /// Nominal artifact relating the exact derived and Machine-IR model versions.
    MachineIrModelVersionCrosswalkIdV1
);
opaque_crosswalk_id!(
    /// Nominal artifact relating the exact derived and Machine-IR frames.
    MachineIrFrameCrosswalkIdV1
);
opaque_crosswalk_id!(
    /// Nominal artifact relating the exact derived and Machine-IR unit systems.
    MachineIrUnitSystemCrosswalkIdV1
);
opaque_crosswalk_id!(
    /// Nominal aggregate declaration tying the four crosswalk artifacts together.
    MachineIrDerivedCrosswalkDeclarationIdV1
);
opaque_crosswalk_id!(
    /// Explicit artifact denying conversion, equivalence, and evidence authority.
    MachineIrCrosswalkNoAuthorityIdV1
);

/// Domain-separated identity schema for one nominal Machine-IR crosswalk candidate.
pub enum DerivedMachineModelCrosswalkCandidateIdentitySchemaV1 {}

impl CanonicalSchema for DerivedMachineModelCrosswalkCandidateIdentitySchemaV1 {
    const DOMAIN: &'static str =
        "org.frankensim.fs-ir.derived-machine-model-crosswalk-candidate.v1";
    const NAME: &'static str = "derived-machine-model-crosswalk-candidate";
    const VERSION: u32 = DERIVED_MACHINE_MODEL_CROSSWALK_CANDIDATE_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "one exact admitted derived geometry, nominal Machine-IR model/subject/version/frame/unit selectors, four nominal endpoint-specific crosswalk artifacts, one aggregate declaration, and an explicit no-authority boundary";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("machine-ir-language-version", WireType::U64),
        FieldSpec::required("machine-model-selector", WireType::Bytes),
        FieldSpec::required("machine-subject-selector", WireType::Bytes),
        FieldSpec::required("machine-model-version-selector", WireType::Bytes),
        FieldSpec::required("machine-frame-selector", WireType::Bytes),
        FieldSpec::required("machine-unit-system-selector", WireType::Bytes),
        FieldSpec::child_of("derived-geometry", &DERIVED_GEOMETRY_CHILD_V1),
        FieldSpec::required("derived-subject", WireType::Bytes),
        FieldSpec::required("derived-model-version", WireType::Bytes),
        FieldSpec::required("derived-frame", WireType::Bytes),
        FieldSpec::required("derived-unit-system", WireType::Bytes),
        FieldSpec::required("subject-crosswalk", WireType::Bytes),
        FieldSpec::required("model-version-crosswalk", WireType::Bytes),
        FieldSpec::required("frame-crosswalk", WireType::Bytes),
        FieldSpec::required("unit-system-crosswalk", WireType::Bytes),
        FieldSpec::required("nominal-crosswalk-declaration", WireType::Bytes),
        FieldSpec::required("no-authority", WireType::Bytes),
    ];
}

/// Typed evidence-node identity of one structural crosswalk candidate.
pub type DerivedMachineModelCrosswalkCandidateIdV1 =
    EvidenceNodeId<DerivedMachineModelCrosswalkCandidateIdentitySchemaV1>;

/// Decoded nominal crosswalk candidate.
///
/// The redundant derived selectors are intentional: admission checks them
/// against `derived_geometry`, so a retained mapping artifact cannot silently
/// float to another subject, model version, frame, or unit system.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DerivedMachineModelCrosswalkCandidateIrV1 {
    /// Decoded candidate schema version.
    pub schema_version: u32,
    /// Exact FrankenScript/Machine-IR language version targeted by the selectors.
    pub machine_ir_version: u32,
    /// Nominal Machine-IR model selector.
    pub machine_model: MachineIrModelSelectorIdV1,
    /// Nominal Machine-IR subject selector.
    pub machine_subject: MachineIrSubjectSelectorIdV1,
    /// Nominal immutable Machine-IR model-version selector.
    pub machine_model_version: MachineIrModelVersionSelectorIdV1,
    /// Nominal Machine-IR frame selector.
    pub machine_frame: MachineIrFrameSelectorIdV1,
    /// Nominal Machine-IR unit-system selector.
    pub machine_unit_system: MachineIrUnitSystemSelectorIdV1,
    /// Exact sealed derived-geometry identity supplied at admission.
    pub derived_geometry: DerivedGeometryIdV1,
    /// Exact derived subject redundantly checked against the sealed geometry.
    pub derived_subject: DerivedSubjectIdV1,
    /// Exact derived model version redundantly checked against the sealed geometry.
    pub derived_model_version: DerivedModelVersionIdV1,
    /// Exact derived frame redundantly checked against the sealed geometry.
    pub derived_frame: DerivedFrameIdV1,
    /// Exact derived unit system redundantly checked against the sealed geometry.
    pub derived_unit_system: DerivedUnitSystemIdV1,
    /// Nominal subject crosswalk artifact.
    pub subject_crosswalk: MachineIrSubjectCrosswalkIdV1,
    /// Nominal immutable-model-version crosswalk artifact.
    pub model_version_crosswalk: MachineIrModelVersionCrosswalkIdV1,
    /// Nominal frame crosswalk artifact.
    pub frame_crosswalk: MachineIrFrameCrosswalkIdV1,
    /// Nominal unit-system crosswalk artifact.
    pub unit_system_crosswalk: MachineIrUnitSystemCrosswalkIdV1,
    /// Nominal declaration joining the endpoint-specific artifacts.
    pub nominal_crosswalk: MachineIrDerivedCrosswalkDeclarationIdV1,
    /// Mandatory explicit denial of authority.
    pub no_authority: MachineIrCrosswalkNoAuthorityIdV1,
}

/// Structured refusal from Machine-IR crosswalk-candidate admission.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DerivedMachineModelCrosswalkCandidateErrorV1 {
    /// Unsupported decoded crosswalk schema version.
    UnsupportedSchemaVersion {
        /// Supplied version.
        found: u32,
        /// Sole supported version.
        supported: u32,
    },
    /// The candidate targets an unsupported Machine-IR language version.
    UnsupportedMachineIrVersion {
        /// Supplied language version.
        found: u32,
        /// Exact version owned by this build.
        supported: u32,
    },
    /// A required selector or artifact used the all-zero sentinel.
    MissingIdentity {
        /// Stable identity field.
        field: &'static str,
    },
    /// The raw derived-geometry identity does not name the supplied sealed object.
    DerivedGeometryIdentityMismatch,
    /// A redundant derived selector differs from the supplied sealed object.
    DerivedSelectorMismatch {
        /// Stable selector field.
        field: &'static str,
    },
    /// Cooperative cancellation was observed before publication.
    Cancelled {
        /// Stable admission stage.
        stage: &'static str,
    },
    /// Canonical identity construction failed.
    Identity(CanonicalError),
}

impl fmt::Display for DerivedMachineModelCrosswalkCandidateErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "derived Machine-IR model crosswalk candidate refused: {self:?}"
        )
    }
}

impl core::error::Error for DerivedMachineModelCrosswalkCandidateErrorV1 {}

/// Sealed nominal crosswalk candidate.
///
/// This token deliberately exposes selectors and lineage only. It is not a
/// `VersionedProgram`, an admitted machine model, a derived morphism, an
/// equivalence, or an evidence-transport capability.
#[derive(Debug, PartialEq, Eq)]
pub struct AdmittedDerivedMachineModelCrosswalkCandidateV1 {
    machine_ir_version: u32,
    machine_model: MachineIrModelSelectorIdV1,
    machine_subject: MachineIrSubjectSelectorIdV1,
    machine_model_version: MachineIrModelVersionSelectorIdV1,
    machine_frame: MachineIrFrameSelectorIdV1,
    machine_unit_system: MachineIrUnitSystemSelectorIdV1,
    derived_geometry: DerivedGeometryIdV1,
    derived_subject: DerivedSubjectIdV1,
    derived_model_version: DerivedModelVersionIdV1,
    derived_frame: DerivedFrameIdV1,
    derived_unit_system: DerivedUnitSystemIdV1,
    subject_crosswalk: MachineIrSubjectCrosswalkIdV1,
    model_version_crosswalk: MachineIrModelVersionCrosswalkIdV1,
    frame_crosswalk: MachineIrFrameCrosswalkIdV1,
    unit_system_crosswalk: MachineIrUnitSystemCrosswalkIdV1,
    nominal_crosswalk: MachineIrDerivedCrosswalkDeclarationIdV1,
    no_authority: MachineIrCrosswalkNoAuthorityIdV1,
    receipt: IdentityReceipt<DerivedMachineModelCrosswalkCandidateIdV1>,
}

impl AdmittedDerivedMachineModelCrosswalkCandidateV1 {
    /// Exact Machine-IR language version named by this candidate.
    #[must_use]
    pub const fn machine_ir_version(&self) -> u32 {
        self.machine_ir_version
    }

    /// Nominal Machine-IR model selector.
    #[must_use]
    pub const fn machine_model(&self) -> MachineIrModelSelectorIdV1 {
        self.machine_model
    }

    /// Nominal Machine-IR subject selector.
    #[must_use]
    pub const fn machine_subject(&self) -> MachineIrSubjectSelectorIdV1 {
        self.machine_subject
    }

    /// Nominal Machine-IR model-version selector.
    #[must_use]
    pub const fn machine_model_version(&self) -> MachineIrModelVersionSelectorIdV1 {
        self.machine_model_version
    }

    /// Nominal Machine-IR frame selector.
    #[must_use]
    pub const fn machine_frame(&self) -> MachineIrFrameSelectorIdV1 {
        self.machine_frame
    }

    /// Nominal Machine-IR unit-system selector.
    #[must_use]
    pub const fn machine_unit_system(&self) -> MachineIrUnitSystemSelectorIdV1 {
        self.machine_unit_system
    }

    /// Exact sealed derived geometry.
    #[must_use]
    pub const fn derived_geometry(&self) -> DerivedGeometryIdV1 {
        self.derived_geometry
    }

    /// Exact derived subject.
    #[must_use]
    pub const fn derived_subject(&self) -> DerivedSubjectIdV1 {
        self.derived_subject
    }

    /// Exact derived model version.
    #[must_use]
    pub const fn derived_model_version(&self) -> DerivedModelVersionIdV1 {
        self.derived_model_version
    }

    /// Exact derived frame.
    #[must_use]
    pub const fn derived_frame(&self) -> DerivedFrameIdV1 {
        self.derived_frame
    }

    /// Exact derived unit system.
    #[must_use]
    pub const fn derived_unit_system(&self) -> DerivedUnitSystemIdV1 {
        self.derived_unit_system
    }

    /// Nominal subject crosswalk artifact.
    #[must_use]
    pub const fn subject_crosswalk(&self) -> MachineIrSubjectCrosswalkIdV1 {
        self.subject_crosswalk
    }

    /// Nominal model-version crosswalk artifact.
    #[must_use]
    pub const fn model_version_crosswalk(&self) -> MachineIrModelVersionCrosswalkIdV1 {
        self.model_version_crosswalk
    }

    /// Nominal frame crosswalk artifact.
    #[must_use]
    pub const fn frame_crosswalk(&self) -> MachineIrFrameCrosswalkIdV1 {
        self.frame_crosswalk
    }

    /// Nominal unit-system crosswalk artifact.
    #[must_use]
    pub const fn unit_system_crosswalk(&self) -> MachineIrUnitSystemCrosswalkIdV1 {
        self.unit_system_crosswalk
    }

    /// Nominal aggregate crosswalk declaration.
    #[must_use]
    pub const fn nominal_crosswalk(&self) -> MachineIrDerivedCrosswalkDeclarationIdV1 {
        self.nominal_crosswalk
    }

    /// Explicit no-authority artifact.
    #[must_use]
    pub const fn no_authority(&self) -> MachineIrCrosswalkNoAuthorityIdV1 {
        self.no_authority
    }

    /// Typed candidate identity.
    #[must_use]
    pub const fn id(&self) -> DerivedMachineModelCrosswalkCandidateIdV1 {
        self.receipt.id()
    }

    /// Canonical receipt and construction limits.
    #[must_use]
    pub const fn identity_receipt(
        &self,
    ) -> IdentityReceipt<DerivedMachineModelCrosswalkCandidateIdV1> {
        self.receipt
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SealedDerivedSelectorsV1 {
    geometry: DerivedGeometryIdV1,
    subject: DerivedSubjectIdV1,
    model_version: DerivedModelVersionIdV1,
    frame: DerivedFrameIdV1,
    unit_system: DerivedUnitSystemIdV1,
}

impl SealedDerivedSelectorsV1 {
    fn from_admitted(geometry: &AdmittedDerivedGeometryV1) -> Self {
        Self {
            geometry: geometry.id(),
            subject: geometry.ir().subject,
            model_version: geometry.ir().model_version,
            frame: geometry.ir().frame,
            unit_system: geometry.ir().unit_system,
        }
    }
}

fn is_zero(bytes: &[u8; 32]) -> bool {
    bytes.iter().all(|byte| *byte == 0)
}

fn validate_candidate(
    ir: &DerivedMachineModelCrosswalkCandidateIrV1,
    sealed: SealedDerivedSelectorsV1,
) -> Result<(), DerivedMachineModelCrosswalkCandidateErrorV1> {
    if ir.schema_version != DERIVED_MACHINE_MODEL_CROSSWALK_CANDIDATE_SCHEMA_VERSION_V1 {
        return Err(
            DerivedMachineModelCrosswalkCandidateErrorV1::UnsupportedSchemaVersion {
                found: ir.schema_version,
                supported: DERIVED_MACHINE_MODEL_CROSSWALK_CANDIDATE_SCHEMA_VERSION_V1,
            },
        );
    }
    if ir.machine_ir_version != IR_VERSION {
        return Err(
            DerivedMachineModelCrosswalkCandidateErrorV1::UnsupportedMachineIrVersion {
                found: ir.machine_ir_version,
                supported: IR_VERSION,
            },
        );
    }

    for (field, bytes) in [
        ("machine-model-selector", ir.machine_model.digest_bytes()),
        (
            "machine-subject-selector",
            ir.machine_subject.digest_bytes(),
        ),
        (
            "machine-model-version-selector",
            ir.machine_model_version.digest_bytes(),
        ),
        ("machine-frame-selector", ir.machine_frame.digest_bytes()),
        (
            "machine-unit-system-selector",
            ir.machine_unit_system.digest_bytes(),
        ),
        ("derived-geometry", ir.derived_geometry.as_bytes()),
        ("derived-subject", ir.derived_subject.as_bytes()),
        ("derived-model-version", ir.derived_model_version.as_bytes()),
        ("derived-frame", ir.derived_frame.as_bytes()),
        ("derived-unit-system", ir.derived_unit_system.as_bytes()),
        ("subject-crosswalk", ir.subject_crosswalk.digest_bytes()),
        (
            "model-version-crosswalk",
            ir.model_version_crosswalk.digest_bytes(),
        ),
        ("frame-crosswalk", ir.frame_crosswalk.digest_bytes()),
        (
            "unit-system-crosswalk",
            ir.unit_system_crosswalk.digest_bytes(),
        ),
        (
            "nominal-crosswalk-declaration",
            ir.nominal_crosswalk.digest_bytes(),
        ),
        ("no-authority", ir.no_authority.digest_bytes()),
    ] {
        if is_zero(bytes) {
            return Err(DerivedMachineModelCrosswalkCandidateErrorV1::MissingIdentity { field });
        }
    }

    if ir.derived_geometry != sealed.geometry {
        return Err(DerivedMachineModelCrosswalkCandidateErrorV1::DerivedGeometryIdentityMismatch);
    }
    for (field, matches) in [
        ("derived-subject", ir.derived_subject == sealed.subject),
        (
            "derived-model-version",
            ir.derived_model_version == sealed.model_version,
        ),
        ("derived-frame", ir.derived_frame == sealed.frame),
        (
            "derived-unit-system",
            ir.derived_unit_system == sealed.unit_system,
        ),
    ] {
        if !matches {
            return Err(
                DerivedMachineModelCrosswalkCandidateErrorV1::DerivedSelectorMismatch { field },
            );
        }
    }
    Ok(())
}

fn candidate_receipt<C: FnMut() -> bool>(
    ir: &DerivedMachineModelCrosswalkCandidateIrV1,
    cancelled: &mut C,
) -> Result<
    IdentityReceipt<DerivedMachineModelCrosswalkCandidateIdV1>,
    DerivedMachineModelCrosswalkCandidateErrorV1,
> {
    let map_identity_error = |error| match error {
        CanonicalError::Cancelled { .. } => {
            DerivedMachineModelCrosswalkCandidateErrorV1::Cancelled {
                stage: "crosswalk-identity",
            }
        }
        other => DerivedMachineModelCrosswalkCandidateErrorV1::Identity(other),
    };
    CanonicalEncoder::<DerivedMachineModelCrosswalkCandidateIdV1, _>::new(
        DERIVED_MACHINE_MODEL_CROSSWALK_IDENTITY_LIMITS_V1,
        || cancelled(),
    )
    .map_err(map_identity_error)?
    .u64(
        Field::new(0, "machine-ir-language-version"),
        u64::from(ir.machine_ir_version),
    )
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(1, "machine-model-selector"),
            ir.machine_model.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(2, "machine-subject-selector"),
            ir.machine_subject.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(3, "machine-model-version-selector"),
            ir.machine_model_version.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(4, "machine-frame-selector"),
            ir.machine_frame.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(5, "machine-unit-system-selector"),
            ir.machine_unit_system.as_bytes(),
        )
    })
    .and_then(|encoder| encoder.child(Field::new(6, "derived-geometry"), ir.derived_geometry))
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(7, "derived-subject"),
            ir.derived_subject.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(8, "derived-model-version"),
            ir.derived_model_version.as_bytes(),
        )
    })
    .and_then(|encoder| encoder.bytes(Field::new(9, "derived-frame"), ir.derived_frame.as_bytes()))
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(10, "derived-unit-system"),
            ir.derived_unit_system.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(11, "subject-crosswalk"),
            ir.subject_crosswalk.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(12, "model-version-crosswalk"),
            ir.model_version_crosswalk.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(13, "frame-crosswalk"),
            ir.frame_crosswalk.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(14, "unit-system-crosswalk"),
            ir.unit_system_crosswalk.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(15, "nominal-crosswalk-declaration"),
            ir.nominal_crosswalk.as_bytes(),
        )
    })
    .and_then(|encoder| encoder.bytes(Field::new(16, "no-authority"), ir.no_authority.as_bytes()))
    .and_then(|encoder| encoder.finish())
    .map_err(map_identity_error)
}

fn admit_candidate_with<C: FnMut() -> bool>(
    ir: &DerivedMachineModelCrosswalkCandidateIrV1,
    sealed: SealedDerivedSelectorsV1,
    mut cancelled: C,
) -> Result<
    AdmittedDerivedMachineModelCrosswalkCandidateV1,
    DerivedMachineModelCrosswalkCandidateErrorV1,
> {
    if cancelled() {
        return Err(DerivedMachineModelCrosswalkCandidateErrorV1::Cancelled {
            stage: "crosswalk-entry",
        });
    }
    validate_candidate(ir, sealed)?;
    let receipt = candidate_receipt(ir, &mut cancelled)?;
    if cancelled() {
        return Err(DerivedMachineModelCrosswalkCandidateErrorV1::Cancelled {
            stage: "crosswalk-publication",
        });
    }
    Ok(AdmittedDerivedMachineModelCrosswalkCandidateV1 {
        machine_ir_version: ir.machine_ir_version,
        machine_model: ir.machine_model,
        machine_subject: ir.machine_subject,
        machine_model_version: ir.machine_model_version,
        machine_frame: ir.machine_frame,
        machine_unit_system: ir.machine_unit_system,
        derived_geometry: ir.derived_geometry,
        derived_subject: ir.derived_subject,
        derived_model_version: ir.derived_model_version,
        derived_frame: ir.derived_frame,
        derived_unit_system: ir.derived_unit_system,
        subject_crosswalk: ir.subject_crosswalk,
        model_version_crosswalk: ir.model_version_crosswalk,
        frame_crosswalk: ir.frame_crosswalk,
        unit_system_crosswalk: ir.unit_system_crosswalk,
        nominal_crosswalk: ir.nominal_crosswalk,
        no_authority: ir.no_authority,
        receipt,
    })
}

/// Bind nominal Machine-IR selectors to one exact admitted derived geometry.
///
/// Admission checks schema and Machine-IR language versions, rejects every
/// zero selector/artifact, checks the raw geometry identity against the supplied
/// sealed object, and checks the redundant derived subject/model/frame/unit
/// selectors against that same object. It does not inspect a Machine-IR model
/// graph or execute any mapping artifact because no admitted machine-side model
/// boundary exists yet.
///
/// # Errors
/// Returns a typed refusal for schema/version, missing identity, raw-to-sealed
/// geometry mismatch, derived-selector mismatch, cancellation, or canonical
/// identity failure. No partial token escapes.
#[must_use = "a raw crosswalk declaration has no conversion or equivalence authority"]
pub fn admit_derived_machine_model_crosswalk_candidate_v1(
    ir: &DerivedMachineModelCrosswalkCandidateIrV1,
    derived_geometry: &AdmittedDerivedGeometryV1,
    cx: &Cx<'_>,
) -> Result<
    AdmittedDerivedMachineModelCrosswalkCandidateV1,
    DerivedMachineModelCrosswalkCandidateErrorV1,
> {
    admit_candidate_with(
        ir,
        SealedDerivedSelectorsV1::from_admitted(derived_geometry),
        || cx.checkpoint().is_err(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn geometry(seed: u8) -> DerivedGeometryIdV1 {
        DerivedGeometryIdV1::parse_slice(&[seed; 32]).expect("32-byte geometry selector")
    }

    fn sealed(seed: u8) -> SealedDerivedSelectorsV1 {
        SealedDerivedSelectorsV1 {
            geometry: geometry(seed),
            subject: DerivedSubjectIdV1::from_bytes([seed.wrapping_add(1); 32]),
            model_version: DerivedModelVersionIdV1::from_bytes([seed.wrapping_add(2); 32]),
            frame: DerivedFrameIdV1::from_bytes([seed.wrapping_add(3); 32]),
            unit_system: DerivedUnitSystemIdV1::from_bytes([seed.wrapping_add(4); 32]),
        }
    }

    fn candidate(sealed: SealedDerivedSelectorsV1) -> DerivedMachineModelCrosswalkCandidateIrV1 {
        DerivedMachineModelCrosswalkCandidateIrV1 {
            schema_version: DERIVED_MACHINE_MODEL_CROSSWALK_CANDIDATE_SCHEMA_VERSION_V1,
            machine_ir_version: IR_VERSION,
            machine_model: MachineIrModelSelectorIdV1::from_bytes([20; 32]),
            machine_subject: MachineIrSubjectSelectorIdV1::from_bytes([21; 32]),
            machine_model_version: MachineIrModelVersionSelectorIdV1::from_bytes([22; 32]),
            machine_frame: MachineIrFrameSelectorIdV1::from_bytes([23; 32]),
            machine_unit_system: MachineIrUnitSystemSelectorIdV1::from_bytes([24; 32]),
            derived_geometry: sealed.geometry,
            derived_subject: sealed.subject,
            derived_model_version: sealed.model_version,
            derived_frame: sealed.frame,
            derived_unit_system: sealed.unit_system,
            subject_crosswalk: MachineIrSubjectCrosswalkIdV1::from_bytes([25; 32]),
            model_version_crosswalk: MachineIrModelVersionCrosswalkIdV1::from_bytes([26; 32]),
            frame_crosswalk: MachineIrFrameCrosswalkIdV1::from_bytes([27; 32]),
            unit_system_crosswalk: MachineIrUnitSystemCrosswalkIdV1::from_bytes([28; 32]),
            nominal_crosswalk: MachineIrDerivedCrosswalkDeclarationIdV1::from_bytes([29; 32]),
            no_authority: MachineIrCrosswalkNoAuthorityIdV1::from_bytes([30; 32]),
        }
    }

    fn admit(
        ir: &DerivedMachineModelCrosswalkCandidateIrV1,
        sealed: SealedDerivedSelectorsV1,
    ) -> Result<
        AdmittedDerivedMachineModelCrosswalkCandidateV1,
        DerivedMachineModelCrosswalkCandidateErrorV1,
    > {
        admit_candidate_with(ir, sealed, || false)
    }

    #[test]
    fn g0_schema_replay_and_accessors_are_exact() {
        let sealed = sealed(10);
        let ir = candidate(sealed);
        let first = admit(&ir, sealed).expect("candidate should admit");
        let second = admit(&ir, sealed).expect("replay should admit");

        assert_eq!(first, second);
        assert_eq!(first.machine_ir_version(), IR_VERSION);
        assert_eq!(first.machine_model(), ir.machine_model);
        assert_eq!(first.machine_subject(), ir.machine_subject);
        assert_eq!(first.machine_model_version(), ir.machine_model_version);
        assert_eq!(first.machine_frame(), ir.machine_frame);
        assert_eq!(first.machine_unit_system(), ir.machine_unit_system);
        assert_eq!(first.derived_geometry(), ir.derived_geometry);
        assert_eq!(first.derived_subject(), ir.derived_subject);
        assert_eq!(first.derived_model_version(), ir.derived_model_version);
        assert_eq!(first.derived_frame(), ir.derived_frame);
        assert_eq!(first.derived_unit_system(), ir.derived_unit_system);
        assert_eq!(first.subject_crosswalk(), ir.subject_crosswalk);
        assert_eq!(first.model_version_crosswalk(), ir.model_version_crosswalk);
        assert_eq!(first.frame_crosswalk(), ir.frame_crosswalk);
        assert_eq!(first.unit_system_crosswalk(), ir.unit_system_crosswalk);
        assert_eq!(first.nominal_crosswalk(), ir.nominal_crosswalk);
        assert_eq!(first.no_authority(), ir.no_authority);
        assert_eq!(first.id(), first.identity_receipt().id());
        assert_eq!(
            DerivedMachineModelCrosswalkCandidateIdentitySchemaV1::FIELDS.len(),
            17
        );
    }

    #[test]
    fn g5_every_retained_identity_moves_the_candidate_identity() {
        let sealed = sealed(40);
        let ir = candidate(sealed);
        let baseline = admit(&ir, sealed).expect("baseline").id();

        let mut changed_version = ir;
        changed_version.machine_ir_version = IR_VERSION.wrapping_add(1);
        let changed_version_id = candidate_receipt(&changed_version, &mut || false)
            .expect("identity encoder accepts an explicitly supplied u32")
            .id();
        assert_ne!(changed_version_id, baseline);

        macro_rules! assert_machine_field_moves {
            ($field:ident, $value:expr) => {{
                let mut changed = ir;
                changed.$field = $value;
                assert_ne!(
                    admit(&changed, sealed).expect(stringify!($field)).id(),
                    baseline
                );
            }};
        }

        assert_machine_field_moves!(
            machine_model,
            MachineIrModelSelectorIdV1::from_bytes([41; 32])
        );
        assert_machine_field_moves!(
            machine_subject,
            MachineIrSubjectSelectorIdV1::from_bytes([42; 32])
        );
        assert_machine_field_moves!(
            machine_model_version,
            MachineIrModelVersionSelectorIdV1::from_bytes([43; 32])
        );
        assert_machine_field_moves!(
            machine_frame,
            MachineIrFrameSelectorIdV1::from_bytes([44; 32])
        );
        assert_machine_field_moves!(
            machine_unit_system,
            MachineIrUnitSystemSelectorIdV1::from_bytes([45; 32])
        );
        assert_machine_field_moves!(
            subject_crosswalk,
            MachineIrSubjectCrosswalkIdV1::from_bytes([46; 32])
        );
        assert_machine_field_moves!(
            model_version_crosswalk,
            MachineIrModelVersionCrosswalkIdV1::from_bytes([47; 32])
        );
        assert_machine_field_moves!(
            frame_crosswalk,
            MachineIrFrameCrosswalkIdV1::from_bytes([48; 32])
        );
        assert_machine_field_moves!(
            unit_system_crosswalk,
            MachineIrUnitSystemCrosswalkIdV1::from_bytes([49; 32])
        );
        assert_machine_field_moves!(
            nominal_crosswalk,
            MachineIrDerivedCrosswalkDeclarationIdV1::from_bytes([50; 32])
        );
        assert_machine_field_moves!(
            no_authority,
            MachineIrCrosswalkNoAuthorityIdV1::from_bytes([51; 32])
        );

        macro_rules! assert_derived_field_moves {
            ($field:ident, $sealed_field:ident, $value:expr) => {{
                let mut changed_ir = ir;
                let mut changed_sealed = sealed;
                changed_ir.$field = $value;
                changed_sealed.$sealed_field = $value;
                assert_ne!(
                    admit(&changed_ir, changed_sealed)
                        .expect(stringify!($field))
                        .id(),
                    baseline
                );
            }};
        }

        assert_derived_field_moves!(derived_geometry, geometry, geometry(52));
        assert_derived_field_moves!(
            derived_subject,
            subject,
            DerivedSubjectIdV1::from_bytes([53; 32])
        );
        assert_derived_field_moves!(
            derived_model_version,
            model_version,
            DerivedModelVersionIdV1::from_bytes([54; 32])
        );
        assert_derived_field_moves!(derived_frame, frame, DerivedFrameIdV1::from_bytes([55; 32]));
        assert_derived_field_moves!(
            derived_unit_system,
            unit_system,
            DerivedUnitSystemIdV1::from_bytes([56; 32])
        );
    }

    #[test]
    fn g3_schema_language_and_sealed_selector_mutations_refuse() {
        let sealed = sealed(60);
        let ir = candidate(sealed);

        let mut changed = ir;
        changed.schema_version += 1;
        assert!(matches!(
            admit(&changed, sealed),
            Err(DerivedMachineModelCrosswalkCandidateErrorV1::UnsupportedSchemaVersion { .. })
        ));

        let mut changed = ir;
        changed.machine_ir_version += 1;
        assert!(matches!(
            admit(&changed, sealed),
            Err(DerivedMachineModelCrosswalkCandidateErrorV1::UnsupportedMachineIrVersion { .. })
        ));

        let mut changed = ir;
        changed.derived_geometry = geometry(61);
        assert_eq!(
            admit(&changed, sealed),
            Err(DerivedMachineModelCrosswalkCandidateErrorV1::DerivedGeometryIdentityMismatch)
        );

        macro_rules! assert_selector_refuses {
            ($field:ident, $value:expr, $name:literal) => {{
                let mut changed = ir;
                changed.$field = $value;
                assert_eq!(
                    admit(&changed, sealed),
                    Err(
                        DerivedMachineModelCrosswalkCandidateErrorV1::DerivedSelectorMismatch {
                            field: $name,
                        }
                    )
                );
            }};
        }

        assert_selector_refuses!(
            derived_subject,
            DerivedSubjectIdV1::from_bytes([62; 32]),
            "derived-subject"
        );
        assert_selector_refuses!(
            derived_model_version,
            DerivedModelVersionIdV1::from_bytes([63; 32]),
            "derived-model-version"
        );
        assert_selector_refuses!(
            derived_frame,
            DerivedFrameIdV1::from_bytes([64; 32]),
            "derived-frame"
        );
        assert_selector_refuses!(
            derived_unit_system,
            DerivedUnitSystemIdV1::from_bytes([65; 32]),
            "derived-unit-system"
        );
    }

    #[test]
    fn g3_zero_selectors_and_artifacts_refuse_before_identity() {
        let sealed = sealed(70);
        let ir = candidate(sealed);

        macro_rules! assert_zero_refuses {
            ($field:ident, $value:expr, $name:literal) => {{
                let mut changed = ir;
                changed.$field = $value;
                assert_eq!(
                    admit(&changed, sealed),
                    Err(
                        DerivedMachineModelCrosswalkCandidateErrorV1::MissingIdentity {
                            field: $name,
                        }
                    )
                );
            }};
        }

        assert_zero_refuses!(
            machine_model,
            MachineIrModelSelectorIdV1::from_bytes([0; 32]),
            "machine-model-selector"
        );
        assert_zero_refuses!(
            machine_subject,
            MachineIrSubjectSelectorIdV1::from_bytes([0; 32]),
            "machine-subject-selector"
        );
        assert_zero_refuses!(
            machine_model_version,
            MachineIrModelVersionSelectorIdV1::from_bytes([0; 32]),
            "machine-model-version-selector"
        );
        assert_zero_refuses!(
            machine_frame,
            MachineIrFrameSelectorIdV1::from_bytes([0; 32]),
            "machine-frame-selector"
        );
        assert_zero_refuses!(
            machine_unit_system,
            MachineIrUnitSystemSelectorIdV1::from_bytes([0; 32]),
            "machine-unit-system-selector"
        );
        assert_zero_refuses!(derived_geometry, geometry(0), "derived-geometry");
        assert_zero_refuses!(
            derived_subject,
            DerivedSubjectIdV1::from_bytes([0; 32]),
            "derived-subject"
        );
        assert_zero_refuses!(
            derived_model_version,
            DerivedModelVersionIdV1::from_bytes([0; 32]),
            "derived-model-version"
        );
        assert_zero_refuses!(
            derived_frame,
            DerivedFrameIdV1::from_bytes([0; 32]),
            "derived-frame"
        );
        assert_zero_refuses!(
            derived_unit_system,
            DerivedUnitSystemIdV1::from_bytes([0; 32]),
            "derived-unit-system"
        );
        assert_zero_refuses!(
            subject_crosswalk,
            MachineIrSubjectCrosswalkIdV1::from_bytes([0; 32]),
            "subject-crosswalk"
        );
        assert_zero_refuses!(
            model_version_crosswalk,
            MachineIrModelVersionCrosswalkIdV1::from_bytes([0; 32]),
            "model-version-crosswalk"
        );
        assert_zero_refuses!(
            frame_crosswalk,
            MachineIrFrameCrosswalkIdV1::from_bytes([0; 32]),
            "frame-crosswalk"
        );
        assert_zero_refuses!(
            unit_system_crosswalk,
            MachineIrUnitSystemCrosswalkIdV1::from_bytes([0; 32]),
            "unit-system-crosswalk"
        );
        assert_zero_refuses!(
            nominal_crosswalk,
            MachineIrDerivedCrosswalkDeclarationIdV1::from_bytes([0; 32]),
            "nominal-crosswalk-declaration"
        );
        assert_zero_refuses!(
            no_authority,
            MachineIrCrosswalkNoAuthorityIdV1::from_bytes([0; 32]),
            "no-authority"
        );
    }

    #[test]
    fn g4_entry_cancellation_publishes_no_candidate() {
        let sealed = sealed(80);
        let ir = candidate(sealed);
        assert_eq!(
            admit_candidate_with(&ir, sealed, || true),
            Err(DerivedMachineModelCrosswalkCandidateErrorV1::Cancelled {
                stage: "crosswalk-entry",
            })
        );
    }

    #[test]
    fn g4_identity_and_publication_cancellation_publish_no_candidate() {
        use core::cell::Cell;

        let sealed = sealed(81);
        let ir = candidate(sealed);

        let identity_probes = Cell::new(0_u32);
        let identity_cancelled = admit_candidate_with(&ir, sealed, || {
            let probe = identity_probes.get().saturating_add(1);
            identity_probes.set(probe);
            probe == 2
        });
        assert_eq!(
            identity_cancelled,
            Err(DerivedMachineModelCrosswalkCandidateErrorV1::Cancelled {
                stage: "crosswalk-identity",
            })
        );

        let baseline_probes = Cell::new(0_u32);
        admit_candidate_with(&ir, sealed, || {
            baseline_probes.set(baseline_probes.get().saturating_add(1));
            false
        })
        .expect("uncancelled probe-count baseline");
        let publication_probe = baseline_probes.get();
        assert!(publication_probe > 2);

        let replay_probes = Cell::new(0_u32);
        let publication_cancelled = admit_candidate_with(&ir, sealed, || {
            let probe = replay_probes.get().saturating_add(1);
            replay_probes.set(probe);
            probe == publication_probe
        });
        assert_eq!(
            publication_cancelled,
            Err(DerivedMachineModelCrosswalkCandidateErrorV1::Cancelled {
                stage: "crosswalk-publication",
            })
        );
    }
}
