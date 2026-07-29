//! Semantic publication selection and pure result-size projections.

use crate::canonical::CanonicalFrameV1;
use crate::catalog::{
    DestinationAdmissionModeV2, PlatformPathProfileV2, PublicationProtocolV2, RunnerCommandV2,
};
use crate::construction::{
    ConstructionErrorKindV2, ConstructionErrorV2, ConstructionFixedObservationV2,
    ConstructionObservedV2,
};
use crate::limits::{
    PublicationStorageProjectionV2, SYSTEM_PUBLICATION_OBJECT_COUNT_V2,
    SystemPublicationObjectRoleV2,
};
use crate::path::{ContentStoreObjectKeyV1, LogicalBundlePathV1};
use fs_blake3::ContentHash;

/// Domain for the local semantic publication-selection projection.
pub const PUBLICATION_SELECTION_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.publication-selection.v1";
/// Canonical identity domain for an abstract whole-publication projection.
pub const PUBLICATION_STORAGE_PROJECTION_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.publication-storage-projection.v1";
/// Maximum abstract artifact rows in one publication projection.
pub const PUBLICATION_STORAGE_ARTIFACT_MAX_V2: usize = 256;
/// Maximum declared non-payload bytes in one abstract ContentStore envelope.
pub const CONTENT_STORE_ENVELOPE_NON_PAYLOAD_MAX_BYTES_V2: u64 = 4096;

/// Maximum complete canonical command-result frame.
pub const COMMAND_RESULT_STDOUT_MAX_BYTES_V2: u64 = 5 * 1024 * 1024;
/// Maximum embedded lifecycle document.
pub const LIFECYCLE_DOCUMENT_MAX_BYTES_V2: u64 = 4 * 1024 * 1024;
/// Maximum embedded catalog.
pub const RUNNER_CATALOG_MAX_BYTES_V2: u64 = 1024 * 1024;
/// Maximum embedded published-bundle receipt.
pub const PUBLISHED_BUNDLE_RECEIPT_MAX_BYTES_V2: u64 = 1024 * 1024;
/// Maximum canonical failure frame written to stderr.
pub const FAILURE_STDERR_MAX_BYTES_V2: u64 = 16 * 1024;

/// Logical destination carried by semantic invocation data.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PublicationTargetV2 {
    /// POSIX descriptor-relative logical bundle path.
    PosixRelative(LogicalBundlePathV1),
    /// Windows handle-relative logical bundle path.
    WindowsRelative(LogicalBundlePathV1),
    /// ContentStore logical object key.
    ContentStoreLogicalKey(ContentStoreObjectKeyV1),
}

impl PublicationTargetV2 {
    /// Path profile implied by this tagged target.
    #[must_use]
    pub const fn path_profile(&self) -> PlatformPathProfileV2 {
        match self {
            Self::PosixRelative(_) => PlatformPathProfileV2::PosixDescriptorRelativeV1,
            Self::WindowsRelative(_) => PlatformPathProfileV2::WindowsHandleRelativeV1,
            Self::ContentStoreLogicalKey(_) => PlatformPathProfileV2::ContentStoreObjectKeyV1,
        }
    }

    /// Exact validated logical bytes.
    #[must_use]
    pub fn as_str(&self) -> &str {
        match self {
            Self::PosixRelative(path) | Self::WindowsRelative(path) => path.as_str(),
            Self::ContentStoreLogicalKey(key) => key.as_str(),
        }
    }
}

/// Frozen semantic publication intent.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct PublicationSelectionV2 {
    path_profile: PlatformPathProfileV2,
    protocol: PublicationProtocolV2,
    destination_mode: DestinationAdmissionModeV2,
    target: PublicationTargetV2,
}

impl PublicationSelectionV2 {
    /// Validate the exact profile/protocol/target compatibility matrix.
    pub fn new(
        path_profile: PlatformPathProfileV2,
        protocol: PublicationProtocolV2,
        destination_mode: DestinationAdmissionModeV2,
        target: PublicationTargetV2,
    ) -> Result<Self, ConstructionErrorV2> {
        let compatible = matches!(
            (&path_profile, &protocol, &target),
            (
                PlatformPathProfileV2::PosixDescriptorRelativeV1,
                PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
                PublicationTargetV2::PosixRelative(_)
            ) | (
                PlatformPathProfileV2::WindowsHandleRelativeV1,
                PublicationProtocolV2::WindowsHandleReplaceAndDirectoryFlushV1,
                PublicationTargetV2::WindowsRelative(_)
            ) | (
                PlatformPathProfileV2::ContentStoreObjectKeyV1,
                PublicationProtocolV2::ContentStoreAtomicCommitV1,
                PublicationTargetV2::ContentStoreLogicalKey(_)
            )
        );
        if !compatible {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "publication.profile_protocol_target",
                "the exact frozen profile/protocol/target cell",
                ConstructionObservedV2::closed_triple(
                    &path_profile,
                    &protocol,
                    &target.path_profile(),
                ),
            ));
        }
        Ok(Self {
            path_profile,
            protocol,
            destination_mode,
            target,
        })
    }

    /// Logical path profile.
    #[must_use]
    pub const fn path_profile(&self) -> PlatformPathProfileV2 {
        self.path_profile
    }

    /// Publication protocol.
    #[must_use]
    pub const fn protocol(&self) -> PublicationProtocolV2 {
        self.protocol
    }

    /// Destination admission mode.
    #[must_use]
    pub const fn destination_mode(&self) -> DestinationAdmissionModeV2 {
        self.destination_mode
    }

    /// Validated logical target.
    #[must_use]
    pub const fn target(&self) -> &PublicationTargetV2 {
        &self.target
    }

    pub(crate) fn semantic_projection_root(&self) -> ContentHash {
        let mut frame = CanonicalFrameV1::new(b"FSPUBSEL\x01", 1024)
            .expect("fixed projection header is within its static bound");
        frame
            .push_u16("publication.path_profile", self.path_profile.code())
            .expect("validated publication selection fits its static bound");
        frame
            .push_u16("publication.protocol", self.protocol.code())
            .expect("validated publication selection fits its static bound");
        frame
            .push_u16("publication.destination_mode", self.destination_mode.code())
            .expect("validated publication selection fits its static bound");
        frame
            .push_u16("publication.target_tag", self.target.path_profile().code())
            .expect("validated publication selection fits its static bound");
        frame
            .push_str("publication.target", self.target.as_str())
            .expect("validated publication target is at most 240 bytes");
        frame.root(PUBLICATION_SELECTION_DOMAIN_V1)
    }
}

impl PublicationStorageProjectionV2<'_> {
    /// Canonically identify every presented storage row and aggregate field.
    ///
    /// This identity method is intentionally non-authoritative and does not
    /// assert that the accounting equations hold. It bounds row counts before
    /// allocation and preserves invalid presented projections so mutation and
    /// refusal evidence can identify the exact offending input.
    pub fn canonical_root_v2(&self) -> Result<ContentHash, ConstructionErrorV2> {
        if self.artifacts.len() > PUBLICATION_STORAGE_ARTIFACT_MAX_V2 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "publication_storage.artifacts",
                "at most 256 abstract artifact rows",
                self.artifacts.len(),
            ));
        }
        if self.system_objects.len() > SYSTEM_PUBLICATION_OBJECT_COUNT_V2 as usize {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "publication_storage.system_objects",
                "at most six abstract system-object rows",
                self.system_objects.len(),
            ));
        }

        let mut frame = CanonicalFrameV1::new(b"FSPSTOR\x01", 64 * 1024)?;
        frame.push_u32(
            "publication_storage.artifact_count",
            u32::try_from(self.artifacts.len()).map_err(|_| {
                ConstructionErrorV2::new(
                    ConstructionErrorKindV2::TooLarge,
                    "publication_storage.artifact_count",
                    "a count representable as u32",
                    self.artifacts.len(),
                )
            })?,
        )?;
        for artifact in self.artifacts {
            encode_storage_row(
                &mut frame,
                "publication_storage.artifact",
                artifact.protocol,
                artifact.encoded_bytes,
                artifact.stored_bytes,
                artifact.envelope_non_payload_bytes,
            )?;
        }
        frame.push_u32(
            "publication_storage.system_object_count",
            u32::try_from(self.system_objects.len()).map_err(|_| {
                ConstructionErrorV2::new(
                    ConstructionErrorKindV2::TooLarge,
                    "publication_storage.system_object_count",
                    "a count representable as u32",
                    self.system_objects.len(),
                )
            })?,
        )?;
        for object in self.system_objects {
            frame.push_u16("publication_storage.system_object.role", object.role as u16)?;
            encode_storage_row(
                &mut frame,
                "publication_storage.system_object",
                object.protocol,
                object.encoded_bytes,
                object.stored_bytes,
                object.envelope_non_payload_bytes,
            )?;
        }
        frame.push_u64(
            "publication_storage.artifact_encoded_bytes",
            self.artifact_encoded_bytes,
        )?;
        frame.push_u64(
            "publication_storage.artifact_stored_bytes",
            self.artifact_stored_bytes,
        )?;
        frame.push_u64(
            "publication_storage.system_publication_stored_bytes",
            self.system_publication_stored_bytes,
        )?;
        frame.push_u64(
            "publication_storage.publication_stored_bytes",
            self.publication_stored_bytes,
        )?;
        Ok(frame.root(PUBLICATION_STORAGE_PROJECTION_DOMAIN_V1))
    }

    /// Validate exact per-row storage relations, all six system roles, checked
    /// aggregate sums, and the whole-publication equation.
    ///
    /// This is abstract length algebra only. It does not validate concrete
    /// bundle paths, envelope bytes, physical store generations, durability,
    /// or publication authority.
    pub fn validate_accounting_v2(&self) -> Result<(), ConstructionErrorV2> {
        if self.artifacts.len() > PUBLICATION_STORAGE_ARTIFACT_MAX_V2 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "publication_storage.artifacts",
                "at most 256 abstract artifact rows",
                self.artifacts.len(),
            ));
        }
        if self.system_objects.len() != SYSTEM_PUBLICATION_OBJECT_COUNT_V2 as usize {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::OutOfRange,
                "publication_storage.system_objects",
                "exactly one row for each of six logical system-object roles",
                self.system_objects.len(),
            ));
        }

        let mut artifact_encoded = 0_u64;
        let mut artifact_stored = 0_u64;
        for artifact in self.artifacts {
            validate_storage_row(
                "publication_storage.artifact",
                artifact.protocol,
                artifact.encoded_bytes,
                artifact.stored_bytes,
                artifact.envelope_non_payload_bytes,
            )?;
            artifact_encoded = checked_storage_sum(
                "publication_storage.artifact_encoded_bytes",
                artifact_encoded,
                artifact.encoded_bytes,
            )?;
            artifact_stored = checked_storage_sum(
                "publication_storage.artifact_stored_bytes",
                artifact_stored,
                artifact.stored_bytes,
            )?;
        }

        let mut system_stored = 0_u64;
        for (object, expected_role) in self
            .system_objects
            .iter()
            .zip(SystemPublicationObjectRoleV2::ALL)
        {
            if object.role != expected_role {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::OutOfOrder,
                    "publication_storage.system_object.role",
                    "all six unique logical roles in canonical order",
                    object.role as u16,
                ));
            }
            validate_storage_row(
                "publication_storage.system_object",
                object.protocol,
                object.encoded_bytes,
                object.stored_bytes,
                object.envelope_non_payload_bytes,
            )?;
            system_stored = checked_storage_sum(
                "publication_storage.system_publication_stored_bytes",
                system_stored,
                object.stored_bytes,
            )?;
        }

        validate_accounting_equality(
            "publication_storage.artifact_encoded_bytes",
            artifact_encoded,
            self.artifact_encoded_bytes,
        )?;
        validate_accounting_equality(
            "publication_storage.artifact_stored_bytes",
            artifact_stored,
            self.artifact_stored_bytes,
        )?;
        validate_accounting_equality(
            "publication_storage.system_publication_stored_bytes",
            system_stored,
            self.system_publication_stored_bytes,
        )?;
        let whole = checked_storage_sum(
            "publication_storage.publication_stored_bytes",
            artifact_stored,
            system_stored,
        )?;
        validate_accounting_equality(
            "publication_storage.publication_stored_bytes",
            whole,
            self.publication_stored_bytes,
        )
    }

    /// Validate accounting and return the same canonical presented-projection
    /// root only when every intrinsic equation holds.
    pub fn validated_canonical_root_v2(&self) -> Result<ContentHash, ConstructionErrorV2> {
        self.validate_accounting_v2()?;
        self.canonical_root_v2()
    }
}

fn encode_storage_row(
    frame: &mut CanonicalFrameV1,
    field: &'static str,
    protocol: PublicationProtocolV2,
    encoded_bytes: u64,
    stored_bytes: u64,
    envelope_non_payload_bytes: u64,
) -> Result<(), ConstructionErrorV2> {
    frame.push_u16(field, protocol.code())?;
    frame.push_u64(field, encoded_bytes)?;
    frame.push_u64(field, stored_bytes)?;
    frame.push_u64(field, envelope_non_payload_bytes)
}

fn validate_storage_row(
    field: &'static str,
    protocol: PublicationProtocolV2,
    encoded_bytes: u64,
    stored_bytes: u64,
    envelope_non_payload_bytes: u64,
) -> Result<(), ConstructionErrorV2> {
    let expected_stored = match protocol {
        PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1
        | PublicationProtocolV2::WindowsHandleReplaceAndDirectoryFlushV1 => {
            if envelope_non_payload_bytes != 0 {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    field,
                    "zero non-payload envelope bytes for Posix or Windows",
                    envelope_non_payload_bytes,
                ));
            }
            encoded_bytes
        }
        PublicationProtocolV2::ContentStoreAtomicCommitV1 => {
            if envelope_non_payload_bytes > CONTENT_STORE_ENVELOPE_NON_PAYLOAD_MAX_BYTES_V2 {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::OutOfRange,
                    field,
                    "at most 4096 ContentStore non-payload envelope bytes",
                    envelope_non_payload_bytes,
                ));
            }
            checked_storage_sum(field, encoded_bytes, envelope_non_payload_bytes)?
        }
    };
    validate_accounting_equality(field, expected_stored, stored_bytes)
}

fn checked_storage_sum(
    field: &'static str,
    left: u64,
    right: u64,
) -> Result<u64, ConstructionErrorV2> {
    left.checked_add(right).ok_or_else(|| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::ArithmeticOverflow,
            field,
            "checked u64 storage accounting",
            right,
        )
    })
}

fn validate_accounting_equality(
    field: &'static str,
    expected: u64,
    observed: u64,
) -> Result<(), ConstructionErrorV2> {
    if expected != observed {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            field,
            "the independently recomputed exact byte count",
            observed,
        ));
    }
    Ok(())
}

/// Pure feasibility projection for one atomic user-visible command result.
///
/// This is length algebra only. It does not construct a result envelope or
/// claim that lifecycle or receipt bytes exist.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SymbolicCommandResultPlanV2 {
    command: RunnerCommandV2,
    framing_bytes: u64,
    lifecycle_bytes: u64,
    catalog_bytes: u64,
    receipt_bytes: u64,
    complete_stdout_bytes: u64,
    resident_staging_bytes: u64,
}

impl SymbolicCommandResultPlanV2 {
    /// Validate exact result-shape presence, checked joint length, stdout
    /// grant, and resident staging feasibility.
    pub fn new(
        command: RunnerCommandV2,
        framing_bytes: u64,
        lifecycle_bytes: u64,
        catalog_bytes: u64,
        receipt_bytes: u64,
        stdout_grant_bytes: u64,
        resident_staging_grant_bytes: u64,
    ) -> Result<Self, ConstructionErrorV2> {
        let (needs_lifecycle, needs_catalog, needs_receipt) = match command {
            RunnerCommandV2::List => (false, true, false),
            RunnerCommandV2::Check | RunnerCommandV2::SelfTest => (true, false, false),
            RunnerCommandV2::Run | RunnerCommandV2::Negative | RunnerCommandV2::Replay => {
                (true, false, true)
            }
        };
        validate_presence("result.lifecycle", needs_lifecycle, lifecycle_bytes)?;
        validate_presence("result.catalog", needs_catalog, catalog_bytes)?;
        validate_presence("result.receipt", needs_receipt, receipt_bytes)?;
        validate_cap(
            "result.lifecycle",
            lifecycle_bytes,
            LIFECYCLE_DOCUMENT_MAX_BYTES_V2,
        )?;
        validate_cap("result.catalog", catalog_bytes, RUNNER_CATALOG_MAX_BYTES_V2)?;
        validate_cap(
            "result.receipt",
            receipt_bytes,
            PUBLISHED_BUNDLE_RECEIPT_MAX_BYTES_V2,
        )?;

        let complete_stdout_bytes = framing_bytes
            .checked_add(lifecycle_bytes)
            .and_then(|value| value.checked_add(catalog_bytes))
            .and_then(|value| value.checked_add(receipt_bytes))
            .ok_or_else(|| {
                ConstructionErrorV2::new(
                    ConstructionErrorKindV2::ArithmeticOverflow,
                    "result.complete_stdout_bytes",
                    "checked u64 sum",
                    ConstructionObservedV2::fixed(ConstructionFixedObservationV2::Overflow),
                )
            })?;
        let admitted_stdout = stdout_grant_bytes.min(COMMAND_RESULT_STDOUT_MAX_BYTES_V2);
        validate_cap(
            "result.complete_stdout_bytes",
            complete_stdout_bytes,
            admitted_stdout,
        )?;

        // Atomic output stages exactly one complete frame. A caller cannot
        // represent an unbudgeted second lifecycle copy through this plan.
        let resident_staging_bytes = complete_stdout_bytes;
        validate_cap(
            "result.resident_staging_bytes",
            resident_staging_bytes,
            resident_staging_grant_bytes,
        )?;

        Ok(Self {
            command,
            framing_bytes,
            lifecycle_bytes,
            catalog_bytes,
            receipt_bytes,
            complete_stdout_bytes,
            resident_staging_bytes,
        })
    }

    /// Command whose result shape was projected.
    #[must_use]
    pub const fn command(&self) -> RunnerCommandV2 {
        self.command
    }

    /// Complete checked stdout-frame length.
    #[must_use]
    pub const fn complete_stdout_bytes(&self) -> u64 {
        self.complete_stdout_bytes
    }

    /// Exact single-copy resident staging requirement.
    #[must_use]
    pub const fn resident_staging_bytes(&self) -> u64 {
        self.resident_staging_bytes
    }

    /// Framing bytes.
    #[must_use]
    pub const fn framing_bytes(&self) -> u64 {
        self.framing_bytes
    }

    /// Embedded lifecycle bytes.
    #[must_use]
    pub const fn lifecycle_bytes(&self) -> u64 {
        self.lifecycle_bytes
    }

    /// Embedded catalog bytes.
    #[must_use]
    pub const fn catalog_bytes(&self) -> u64 {
        self.catalog_bytes
    }

    /// Embedded receipt bytes.
    #[must_use]
    pub const fn receipt_bytes(&self) -> u64 {
        self.receipt_bytes
    }
}

/// Validate the bounded canonical failure-frame length.
pub fn validate_failure_stderr_bytes_v2(length: u64) -> Result<(), ConstructionErrorV2> {
    if length == 0 {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Missing,
            "result.failure_stderr_bytes",
            "one nonempty canonical failure frame",
            length,
        ));
    }
    validate_cap(
        "result.failure_stderr_bytes",
        length,
        FAILURE_STDERR_MAX_BYTES_V2,
    )
}

fn validate_presence(
    field: &'static str,
    required: bool,
    bytes: u64,
) -> Result<(), ConstructionErrorV2> {
    match (required, bytes == 0) {
        (true, true) => Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Missing,
            field,
            "nonzero bytes for the command result variant",
            bytes,
        )),
        (false, false) => Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Unexpected,
            field,
            "zero bytes for the command result variant",
            bytes,
        )),
        _ => Ok(()),
    }
}

fn validate_cap(field: &'static str, observed: u64, cap: u64) -> Result<(), ConstructionErrorV2> {
    if observed > cap {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::OutOfRange,
            field,
            "value at or below its admitted byte cap",
            observed,
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        COMMAND_RESULT_STDOUT_MAX_BYTES_V2, CONTENT_STORE_ENVELOPE_NON_PAYLOAD_MAX_BYTES_V2,
        FAILURE_STDERR_MAX_BYTES_V2, LIFECYCLE_DOCUMENT_MAX_BYTES_V2,
        PUBLICATION_STORAGE_ARTIFACT_MAX_V2, PUBLISHED_BUNDLE_RECEIPT_MAX_BYTES_V2,
        PublicationSelectionV2, PublicationTargetV2, RUNNER_CATALOG_MAX_BYTES_V2,
        SymbolicCommandResultPlanV2, validate_failure_stderr_bytes_v2,
    };
    use crate::catalog::{
        DestinationAdmissionModeV2, PlatformPathProfileV2, PublicationProtocolV2, RunnerCommandV2,
    };
    use crate::path::{ContentStoreObjectKeyV1, LogicalBundlePathV1};
    use crate::{
        ConstructionErrorKindV2,
        limits::{
            ArtifactStorageProjectionV2, PublicationStorageProjectionV2,
            SystemObjectStorageProjectionV2, SystemPublicationObjectRoleV2,
        },
    };

    fn targets() -> [PublicationTargetV2; 3] {
        [
            PublicationTargetV2::PosixRelative(
                LogicalBundlePathV1::new("runner/seal").expect("POSIX fixture"),
            ),
            PublicationTargetV2::WindowsRelative(
                LogicalBundlePathV1::new("runner/seal").expect("Windows fixture"),
            ),
            PublicationTargetV2::ContentStoreLogicalKey(
                ContentStoreObjectKeyV1::new("objects/seal").expect("ContentStore fixture"),
            ),
        ]
    }

    fn system_objects(
        protocol: PublicationProtocolV2,
        encoded_bytes: u64,
        envelope_non_payload_bytes: u64,
    ) -> [SystemObjectStorageProjectionV2; 6] {
        let stored_bytes = encoded_bytes
            .checked_add(envelope_non_payload_bytes)
            .expect("bounded fixture");
        SystemPublicationObjectRoleV2::ALL.map(|role| SystemObjectStorageProjectionV2 {
            role,
            protocol,
            encoded_bytes,
            stored_bytes,
            envelope_non_payload_bytes,
        })
    }

    fn publication_projection<'a>(
        artifacts: &'a [ArtifactStorageProjectionV2],
        system_objects: &'a [SystemObjectStorageProjectionV2],
    ) -> PublicationStorageProjectionV2<'a> {
        let artifact_encoded_bytes = artifacts
            .iter()
            .map(|artifact| artifact.encoded_bytes)
            .sum();
        let artifact_stored_bytes = artifacts.iter().map(|artifact| artifact.stored_bytes).sum();
        let system_publication_stored_bytes = system_objects
            .iter()
            .map(|object| object.stored_bytes)
            .sum();
        PublicationStorageProjectionV2 {
            artifacts,
            system_objects,
            artifact_encoded_bytes,
            artifact_stored_bytes,
            system_publication_stored_bytes,
            publication_stored_bytes: artifact_stored_bytes + system_publication_stored_bytes,
        }
    }

    #[test]
    fn profile_protocol_and_target_are_one_exact_cell() {
        let path = LogicalBundlePathV1::new("runner/seal").expect("valid logical path");
        assert!(
            PublicationSelectionV2::new(
                PlatformPathProfileV2::PosixDescriptorRelativeV1,
                PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
                DestinationAdmissionModeV2::Absent,
                PublicationTargetV2::PosixRelative(path.clone()),
            )
            .is_ok()
        );
        assert!(
            PublicationSelectionV2::new(
                PlatformPathProfileV2::WindowsHandleRelativeV1,
                PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
                DestinationAdmissionModeV2::Absent,
                PublicationTargetV2::WindowsRelative(path),
            )
            .is_err()
        );
    }

    #[test]
    fn compatibility_matrix_accepts_exactly_three_cells_per_destination_mode() {
        let profiles = [
            PlatformPathProfileV2::PosixDescriptorRelativeV1,
            PlatformPathProfileV2::WindowsHandleRelativeV1,
            PlatformPathProfileV2::ContentStoreObjectKeyV1,
        ];
        let protocols = [
            PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
            PublicationProtocolV2::WindowsHandleReplaceAndDirectoryFlushV1,
            PublicationProtocolV2::ContentStoreAtomicCommitV1,
        ];
        let modes = [
            DestinationAdmissionModeV2::Absent,
            DestinationAdmissionModeV2::PreExistingEmpty,
        ];
        let mut accepted = 0;
        for (profile_index, profile) in profiles.into_iter().enumerate() {
            for (protocol_index, protocol) in protocols.into_iter().enumerate() {
                for (target_index, target) in targets().into_iter().enumerate() {
                    for mode in modes {
                        let result =
                            PublicationSelectionV2::new(profile, protocol, mode, target.clone());
                        let exact =
                            profile_index == protocol_index && protocol_index == target_index;
                        assert_eq!(
                            result.is_ok(),
                            exact,
                            "{profile:?}/{protocol:?}/{:?}/{mode:?}",
                            target.path_profile()
                        );
                        accepted += usize::from(result.is_ok());
                    }
                }
            }
        }
        assert_eq!(accepted, 6);
    }

    #[test]
    fn each_semantic_publication_mutation_moves_the_projection_root() {
        let make = |mode, path: &str| {
            PublicationSelectionV2::new(
                PlatformPathProfileV2::PosixDescriptorRelativeV1,
                PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
                mode,
                PublicationTargetV2::PosixRelative(
                    LogicalBundlePathV1::new(path).expect("fixture path"),
                ),
            )
            .expect("fixture selection")
        };
        let base = make(DestinationAdmissionModeV2::Absent, "runner/seal");
        assert_ne!(
            base.semantic_projection_root(),
            make(DestinationAdmissionModeV2::PreExistingEmpty, "runner/seal")
                .semantic_projection_root()
        );
        assert_ne!(
            base.semantic_projection_root(),
            make(DestinationAdmissionModeV2::Absent, "runner/other").semantic_projection_root()
        );
        let windows = PublicationSelectionV2::new(
            PlatformPathProfileV2::WindowsHandleRelativeV1,
            PublicationProtocolV2::WindowsHandleReplaceAndDirectoryFlushV1,
            DestinationAdmissionModeV2::Absent,
            targets()[1].clone(),
        )
        .expect("Windows selection");
        let content_store = PublicationSelectionV2::new(
            PlatformPathProfileV2::ContentStoreObjectKeyV1,
            PublicationProtocolV2::ContentStoreAtomicCommitV1,
            DestinationAdmissionModeV2::Absent,
            targets()[2].clone(),
        )
        .expect("ContentStore selection");
        assert_ne!(
            base.semantic_projection_root(),
            windows.semantic_projection_root()
        );
        assert_ne!(
            base.semantic_projection_root(),
            content_store.semantic_projection_root()
        );
    }

    #[test]
    fn physical_observations_have_no_semantic_selection_field() {
        let selection = PublicationSelectionV2::new(
            PlatformPathProfileV2::ContentStoreObjectKeyV1,
            PublicationProtocolV2::ContentStoreAtomicCommitV1,
            DestinationAdmissionModeV2::PreExistingEmpty,
            PublicationTargetV2::ContentStoreLogicalKey(
                ContentStoreObjectKeyV1::new("objects/seal").expect("valid object key"),
            ),
        )
        .expect("compatible selection");
        let before = selection.semantic_projection_root();
        let opaque_generation_before = 17_u64;
        let opaque_generation_after = 18_u64;
        assert_ne!(opaque_generation_before, opaque_generation_after);
        let after = selection.clone().semantic_projection_root();
        assert_eq!(before, after);
    }

    #[test]
    fn whole_publication_projection_has_one_domain_separated_canonical_root() {
        let artifacts = [
            ArtifactStorageProjectionV2 {
                protocol: PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
                encoded_bytes: 17,
                stored_bytes: 17,
                envelope_non_payload_bytes: 0,
            },
            ArtifactStorageProjectionV2 {
                protocol: PublicationProtocolV2::ContentStoreAtomicCommitV1,
                encoded_bytes: 23,
                stored_bytes: 30,
                envelope_non_payload_bytes: 7,
            },
        ];
        let system_objects =
            system_objects(PublicationProtocolV2::ContentStoreAtomicCommitV1, 5, 3);
        let projection = publication_projection(&artifacts, &system_objects);
        projection
            .validate_accounting_v2()
            .expect("independently recomputed accounting");
        let root = projection
            .validated_canonical_root_v2()
            .expect("validated canonical root");
        assert_eq!(
            root,
            projection
                .canonical_root_v2()
                .expect("presented canonical root")
        );
        assert_eq!(
            root,
            publication_projection(&artifacts, &system_objects)
                .validated_canonical_root_v2()
                .expect("deterministic reconstruction")
        );
    }

    #[test]
    fn every_publication_storage_field_moves_identity_and_bad_accounting_refuses() {
        let artifacts = [ArtifactStorageProjectionV2 {
            protocol: PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
            encoded_bytes: 17,
            stored_bytes: 17,
            envelope_non_payload_bytes: 0,
        }];
        let system_objects = system_objects(
            PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
            5,
            0,
        );
        let base_projection = publication_projection(&artifacts, &system_objects);
        let base = base_projection
            .canonical_root_v2()
            .expect("base projection root");

        let mut artifact_protocol = artifacts;
        artifact_protocol[0].protocol = PublicationProtocolV2::ContentStoreAtomicCommitV1;
        let protocol_projection = publication_projection(&artifact_protocol, &system_objects);
        assert_ne!(
            protocol_projection
                .canonical_root_v2()
                .expect("protocol root"),
            base
        );
        protocol_projection
            .validate_accounting_v2()
            .expect("zero-overhead ContentStore row remains a valid relation");

        for field in 0..3 {
            let mut mutation = artifacts;
            match field {
                0 => mutation[0].encoded_bytes += 1,
                1 => mutation[0].stored_bytes += 1,
                2 => mutation[0].envelope_non_payload_bytes += 1,
                _ => unreachable!("three artifact fields"),
            }
            let mutated = PublicationStorageProjectionV2 {
                artifacts: &mutation,
                ..base_projection
            };
            assert_ne!(
                mutated.canonical_root_v2().expect("artifact mutation root"),
                base
            );
            assert!(mutated.validate_accounting_v2().is_err());
        }

        for field in 0..5 {
            let mut mutation = system_objects;
            match field {
                0 => mutation[0].role = SystemPublicationObjectRoleV2::RunTerminal,
                1 => mutation[0].protocol = PublicationProtocolV2::ContentStoreAtomicCommitV1,
                2 => mutation[0].encoded_bytes += 1,
                3 => mutation[0].stored_bytes += 1,
                4 => mutation[0].envelope_non_payload_bytes += 1,
                _ => unreachable!("five system-object fields"),
            }
            let mutated = PublicationStorageProjectionV2 {
                system_objects: &mutation,
                ..base_projection
            };
            assert_ne!(
                mutated
                    .canonical_root_v2()
                    .expect("system-object mutation root"),
                base
            );
            if field != 1 {
                assert!(mutated.validate_accounting_v2().is_err());
            }
        }

        for field in 0..4 {
            let mut mutated = base_projection;
            match field {
                0 => mutated.artifact_encoded_bytes += 1,
                1 => mutated.artifact_stored_bytes += 1,
                2 => mutated.system_publication_stored_bytes += 1,
                3 => mutated.publication_stored_bytes += 1,
                _ => unreachable!("four aggregate fields"),
            }
            assert_ne!(
                mutated
                    .canonical_root_v2()
                    .expect("aggregate mutation root"),
                base
            );
            assert!(mutated.validate_accounting_v2().is_err());
        }
    }

    #[test]
    fn publication_storage_zero_max_role_envelope_and_overflow_boundaries_are_exact() {
        let zero_system = system_objects(
            PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
            0,
            0,
        );
        publication_projection(&[], &zero_system)
            .validate_accounting_v2()
            .expect("zero artifacts and six zero-byte system objects");

        let maximum_artifacts = vec![
            ArtifactStorageProjectionV2 {
                protocol: PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
                encoded_bytes: 1,
                stored_bytes: 1,
                envelope_non_payload_bytes: 0,
            };
            PUBLICATION_STORAGE_ARTIFACT_MAX_V2
        ];
        publication_projection(&maximum_artifacts, &zero_system)
            .validated_canonical_root_v2()
            .expect("256 artifacts are structurally and arithmetically admitted");
        let one_over_artifacts = vec![
            ArtifactStorageProjectionV2 {
                protocol: PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
                encoded_bytes: 0,
                stored_bytes: 0,
                envelope_non_payload_bytes: 0,
            };
            PUBLICATION_STORAGE_ARTIFACT_MAX_V2 + 1
        ];
        assert_eq!(
            publication_projection(&one_over_artifacts, &zero_system)
                .validate_accounting_v2()
                .expect_err("257 artifacts")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );

        let exact_envelope = [ArtifactStorageProjectionV2 {
            protocol: PublicationProtocolV2::ContentStoreAtomicCommitV1,
            encoded_bytes: 1,
            stored_bytes: 1 + CONTENT_STORE_ENVELOPE_NON_PAYLOAD_MAX_BYTES_V2,
            envelope_non_payload_bytes: CONTENT_STORE_ENVELOPE_NON_PAYLOAD_MAX_BYTES_V2,
        }];
        publication_projection(&exact_envelope, &zero_system)
            .validate_accounting_v2()
            .expect("exact ContentStore envelope boundary");
        let one_over_envelope = [ArtifactStorageProjectionV2 {
            protocol: PublicationProtocolV2::ContentStoreAtomicCommitV1,
            encoded_bytes: 1,
            stored_bytes: 2 + CONTENT_STORE_ENVELOPE_NON_PAYLOAD_MAX_BYTES_V2,
            envelope_non_payload_bytes: CONTENT_STORE_ENVELOPE_NON_PAYLOAD_MAX_BYTES_V2 + 1,
        }];
        assert!(
            publication_projection(&one_over_envelope, &zero_system)
                .validate_accounting_v2()
                .is_err()
        );

        let mut duplicate_role = zero_system;
        duplicate_role[0].role = SystemPublicationObjectRoleV2::RunTerminal;
        assert!(
            publication_projection(&[], &duplicate_role)
                .validate_accounting_v2()
                .is_err()
        );
        assert!(
            publication_projection(&[], &zero_system[..5])
                .validate_accounting_v2()
                .is_err()
        );

        let overflow = [
            ArtifactStorageProjectionV2 {
                protocol: PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
                encoded_bytes: u64::MAX,
                stored_bytes: u64::MAX,
                envelope_non_payload_bytes: 0,
            },
            ArtifactStorageProjectionV2 {
                protocol: PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
                encoded_bytes: 1,
                stored_bytes: 1,
                envelope_non_payload_bytes: 0,
            },
        ];
        let overflow_projection = PublicationStorageProjectionV2 {
            artifacts: &overflow,
            system_objects: &zero_system,
            artifact_encoded_bytes: u64::MAX,
            artifact_stored_bytes: u64::MAX,
            system_publication_stored_bytes: 0,
            publication_stored_bytes: u64::MAX,
        };
        assert_eq!(
            overflow_projection
                .validate_accounting_v2()
                .expect_err("artifact sum overflow")
                .kind(),
            ConstructionErrorKindV2::ArithmeticOverflow
        );
    }

    #[test]
    fn atomic_result_projection_rejects_wrong_presence_and_second_copy_pressure() {
        assert!(
            SymbolicCommandResultPlanV2::new(RunnerCommandV2::List, 32, 0, 128, 0, 1024, 1024,)
                .is_ok()
        );
        assert!(
            SymbolicCommandResultPlanV2::new(RunnerCommandV2::Run, 32, 256, 0, 128, 1024, 300,)
                .is_err()
        );
        assert!(
            SymbolicCommandResultPlanV2::new(RunnerCommandV2::Run, 32, 256, 1, 128, 1024, 1024,)
                .is_err()
        );
    }

    #[test]
    fn result_and_failure_caps_are_inclusive() {
        assert!(
            SymbolicCommandResultPlanV2::new(
                RunnerCommandV2::List,
                COMMAND_RESULT_STDOUT_MAX_BYTES_V2 - 1,
                0,
                1,
                0,
                COMMAND_RESULT_STDOUT_MAX_BYTES_V2,
                COMMAND_RESULT_STDOUT_MAX_BYTES_V2,
            )
            .is_ok()
        );
        assert!(validate_failure_stderr_bytes_v2(FAILURE_STDERR_MAX_BYTES_V2).is_ok());
        assert!(validate_failure_stderr_bytes_v2(FAILURE_STDERR_MAX_BYTES_V2 + 1).is_err());
        assert!(validate_failure_stderr_bytes_v2(0).is_err());
    }

    #[test]
    fn every_command_result_shape_and_nested_boundary_is_exact() {
        let cases = [
            (RunnerCommandV2::List, 0, 1, 0),
            (RunnerCommandV2::Check, 1, 0, 0),
            (RunnerCommandV2::SelfTest, 1, 0, 0),
            (RunnerCommandV2::Run, 1, 0, 1),
            (RunnerCommandV2::Negative, 1, 0, 1),
            (RunnerCommandV2::Replay, 1, 0, 1),
        ];
        for (command, lifecycle, catalog, receipt) in cases {
            let plan =
                SymbolicCommandResultPlanV2::new(command, 1, lifecycle, catalog, receipt, 4, 4)
                    .expect("minimal exact command shape");
            assert_eq!(plan.command(), command);
            assert_eq!(plan.complete_stdout_bytes(), 2 + receipt);
            assert_eq!(plan.resident_staging_bytes(), plan.complete_stdout_bytes());

            assert!(
                SymbolicCommandResultPlanV2::new(
                    command,
                    1,
                    u64::from(lifecycle == 0),
                    catalog,
                    receipt,
                    8,
                    8,
                )
                .is_err(),
                "{command:?} lifecycle presence"
            );
        }

        assert!(
            SymbolicCommandResultPlanV2::new(
                RunnerCommandV2::Check,
                1,
                LIFECYCLE_DOCUMENT_MAX_BYTES_V2,
                0,
                0,
                COMMAND_RESULT_STDOUT_MAX_BYTES_V2,
                COMMAND_RESULT_STDOUT_MAX_BYTES_V2,
            )
            .is_ok()
        );
        assert!(
            SymbolicCommandResultPlanV2::new(
                RunnerCommandV2::Check,
                1,
                LIFECYCLE_DOCUMENT_MAX_BYTES_V2 + 1,
                0,
                0,
                COMMAND_RESULT_STDOUT_MAX_BYTES_V2,
                COMMAND_RESULT_STDOUT_MAX_BYTES_V2,
            )
            .is_err()
        );
        assert!(
            SymbolicCommandResultPlanV2::new(
                RunnerCommandV2::List,
                1,
                0,
                RUNNER_CATALOG_MAX_BYTES_V2 + 1,
                0,
                COMMAND_RESULT_STDOUT_MAX_BYTES_V2,
                COMMAND_RESULT_STDOUT_MAX_BYTES_V2,
            )
            .is_err()
        );
        assert!(
            SymbolicCommandResultPlanV2::new(
                RunnerCommandV2::Run,
                1,
                1,
                0,
                PUBLISHED_BUNDLE_RECEIPT_MAX_BYTES_V2 + 1,
                COMMAND_RESULT_STDOUT_MAX_BYTES_V2,
                COMMAND_RESULT_STDOUT_MAX_BYTES_V2,
            )
            .is_err()
        );
        assert!(
            SymbolicCommandResultPlanV2::new(
                RunnerCommandV2::Run,
                u64::MAX,
                1,
                0,
                1,
                u64::MAX,
                u64::MAX,
            )
            .is_err()
        );
        assert!(
            SymbolicCommandResultPlanV2::new(RunnerCommandV2::Run, 1, 1, 0, 1, 3, 2,).is_err(),
            "the exact one-copy staging requirement cannot exceed its grant"
        );
    }
}
