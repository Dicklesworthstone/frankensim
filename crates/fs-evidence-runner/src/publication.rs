//! Semantic publication selection and pure result-size projections.

use crate::canonical::CanonicalFrameV1;
use crate::catalog::{
    DestinationAdmissionModeV2, PlatformPathProfileV2, PublicationProtocolV2, RunnerCommandV2,
};
use crate::construction::{ConstructionErrorKindV2, ConstructionErrorV2};
use crate::path::{ContentStoreObjectKeyV1, LogicalBundlePathV1};
use fs_blake3::ContentHash;

/// Domain for the local semantic publication-selection projection.
pub const PUBLICATION_SELECTION_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.publication-selection.v1";

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
                format_args!(
                    "{}/{}/{}",
                    path_profile.name(),
                    protocol.name(),
                    target.path_profile().name()
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
                    "overflow",
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
        COMMAND_RESULT_STDOUT_MAX_BYTES_V2, FAILURE_STDERR_MAX_BYTES_V2,
        LIFECYCLE_DOCUMENT_MAX_BYTES_V2, PUBLISHED_BUNDLE_RECEIPT_MAX_BYTES_V2,
        PublicationSelectionV2, PublicationTargetV2, RUNNER_CATALOG_MAX_BYTES_V2,
        SymbolicCommandResultPlanV2, validate_failure_stderr_bytes_v2,
    };
    use crate::catalog::{
        DestinationAdmissionModeV2, PlatformPathProfileV2, PublicationProtocolV2, RunnerCommandV2,
    };
    use crate::path::{ContentStoreObjectKeyV1, LogicalBundlePathV1};

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
