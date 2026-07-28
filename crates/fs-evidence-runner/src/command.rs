//! Frozen command-selection provenance and applicability table.

use crate::budget::AdmittedRunnerBudgetsV2;
use crate::catalog::{ArtifactDispositionV2, RunProfileV2, RunnerCommandV2};
use crate::construction::{ConstructionErrorKindV2, ConstructionErrorV2};
use crate::identity::{CaseManifestRootV2, SourceIdentityRootV2};
use crate::publication::PublicationSelectionV2;
use crate::value::StableTokenV2;

/// Exact source of family, mode, and profile for one command.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CommandSelectionProvenanceV2 {
    /// Fixed sealed preflight manifest for `Check`.
    FixedPreflight(CaseManifestRootV2),
    /// Fixed sealed bounded internal manifest for `SelfTest`.
    FixedSelfTest(CaseManifestRootV2),
    /// Explicit caller selection for `Run`.
    CallerRun,
    /// Immutable sealed negative-case manifest.
    SealedNegative(CaseManifestRootV2),
    /// Immutable sealed source lineage for `Replay`.
    SealedReplay(SourceIdentityRootV2),
}

/// Family/mode/profile selected through an exact provenance route.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CommandSelectionV2 {
    provenance: CommandSelectionProvenanceV2,
    family: StableTokenV2,
    mode: StableTokenV2,
    profile: RunProfileV2,
}

impl CommandSelectionV2 {
    /// Selection supplied only by the fixed `Check` preflight manifest.
    #[must_use]
    pub const fn fixed_preflight(
        manifest: CaseManifestRootV2,
        family: StableTokenV2,
        mode: StableTokenV2,
        profile: RunProfileV2,
    ) -> Self {
        Self {
            provenance: CommandSelectionProvenanceV2::FixedPreflight(manifest),
            family,
            mode,
            profile,
        }
    }

    /// Selection supplied only by the fixed bounded `SelfTest` manifest.
    #[must_use]
    pub const fn fixed_self_test(
        manifest: CaseManifestRootV2,
        family: StableTokenV2,
        mode: StableTokenV2,
        profile: RunProfileV2,
    ) -> Self {
        Self {
            provenance: CommandSelectionProvenanceV2::FixedSelfTest(manifest),
            family,
            mode,
            profile,
        }
    }

    /// Explicit caller family/mode/profile for `Run`.
    #[must_use]
    pub const fn caller_run(
        family: StableTokenV2,
        mode: StableTokenV2,
        profile: RunProfileV2,
    ) -> Self {
        Self {
            provenance: CommandSelectionProvenanceV2::CallerRun,
            family,
            mode,
            profile,
        }
    }

    /// Selection supplied by one immutable negative-case manifest.
    #[must_use]
    pub const fn sealed_negative(
        manifest: CaseManifestRootV2,
        family: StableTokenV2,
        mode: StableTokenV2,
        profile: RunProfileV2,
    ) -> Self {
        Self {
            provenance: CommandSelectionProvenanceV2::SealedNegative(manifest),
            family,
            mode,
            profile,
        }
    }

    /// Selection supplied by one immutable replay source lineage.
    #[must_use]
    pub const fn sealed_replay(
        source: SourceIdentityRootV2,
        family: StableTokenV2,
        mode: StableTokenV2,
        profile: RunProfileV2,
    ) -> Self {
        Self {
            provenance: CommandSelectionProvenanceV2::SealedReplay(source),
            family,
            mode,
            profile,
        }
    }

    /// Exact selection provenance.
    #[must_use]
    pub const fn provenance(&self) -> &CommandSelectionProvenanceV2 {
        &self.provenance
    }

    /// Selected registered family token.
    #[must_use]
    pub const fn family(&self) -> &StableTokenV2 {
        &self.family
    }

    /// Selected registered mode token.
    #[must_use]
    pub const fn mode(&self) -> &StableTokenV2 {
        &self.mode
    }

    /// Selected profile.
    #[must_use]
    pub const fn profile(&self) -> RunProfileV2 {
        self.profile
    }
}

/// Validated command intent after sealed-manifest projection and budget
/// admission. It is not a parser result, lifecycle record, or execution grant.
///
/// Validated fields cannot be mutated after construction:
///
/// ```compile_fail
/// use fs_evidence_runner::CommandIntentV2;
///
/// fn forge_command(intent: &mut CommandIntentV2) {
///     intent.command = fs_evidence_runner::RunnerCommandV2::Replay;
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CommandIntentV2 {
    command: RunnerCommandV2,
    selection: Option<CommandSelectionV2>,
    budgets: Option<AdmittedRunnerBudgetsV2>,
    disposition: Option<ArtifactDispositionV2>,
    publication_selection: Option<PublicationSelectionV2>,
}

impl CommandIntentV2 {
    /// Construct `List`, whose selectors, budgets, disposition, and
    /// publication selection are all typed absence.
    #[must_use]
    pub const fn list() -> Self {
        Self {
            command: RunnerCommandV2::List,
            selection: None,
            budgets: None,
            disposition: None,
            publication_selection: None,
        }
    }

    /// Validate the exact command/provenance/profile/disposition/publication
    /// table for every non-List command.
    pub fn new(
        command: RunnerCommandV2,
        selection: CommandSelectionV2,
        budgets: AdmittedRunnerBudgetsV2,
        publication_selection: Option<PublicationSelectionV2>,
    ) -> Result<Self, ConstructionErrorV2> {
        if command == RunnerCommandV2::List {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Unexpected,
                "command.selection",
                "typed absence for List; use CommandIntentV2::list",
                "present",
            ));
        }
        let provenance_matches = matches!(
            (command, selection.provenance()),
            (
                RunnerCommandV2::Check,
                CommandSelectionProvenanceV2::FixedPreflight(_)
            ) | (
                RunnerCommandV2::SelfTest,
                CommandSelectionProvenanceV2::FixedSelfTest(_)
            ) | (
                RunnerCommandV2::Run,
                CommandSelectionProvenanceV2::CallerRun
            ) | (
                RunnerCommandV2::Negative,
                CommandSelectionProvenanceV2::SealedNegative(_)
            ) | (
                RunnerCommandV2::Replay,
                CommandSelectionProvenanceV2::SealedReplay(_)
            )
        );
        if !provenance_matches {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "command.selection_provenance",
                "the frozen command-specific source",
                format_args!(
                    "{}/{}",
                    command.name(),
                    provenance_name(selection.provenance())
                ),
            ));
        }

        if budgets.profile() != selection.profile {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "command.profile",
                "the same profile in selection and admitted budgets",
                format_args!("{}/{}", selection.profile.name(), budgets.profile().name()),
            ));
        }
        let disposition = match command {
            RunnerCommandV2::Check | RunnerCommandV2::SelfTest => {
                ArtifactDispositionV2::LifecycleOnlyNoBundle
            }
            RunnerCommandV2::Run | RunnerCommandV2::Negative | RunnerCommandV2::Replay => {
                ArtifactDispositionV2::DurableBundleRequired
            }
            RunnerCommandV2::List => unreachable!("List was handled before the table"),
        };
        if budgets.disposition() != disposition {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "command.disposition",
                "the command's fixed disposition",
                budgets.disposition().name(),
            ));
        }
        match (disposition, publication_selection.is_some()) {
            (ArtifactDispositionV2::LifecycleOnlyNoBundle, false)
            | (ArtifactDispositionV2::DurableBundleRequired, true) => {}
            (ArtifactDispositionV2::LifecycleOnlyNoBundle, true) => {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Unexpected,
                    "command.publication_selection",
                    "absence for a lifecycle-only command",
                    "present",
                ));
            }
            (ArtifactDispositionV2::DurableBundleRequired, false) => {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Missing,
                    "command.publication_selection",
                    "one selection for a durable command",
                    "absent",
                ));
            }
        }

        Ok(Self {
            command,
            selection: Some(selection),
            budgets: Some(budgets),
            disposition: Some(disposition),
            publication_selection,
        })
    }

    /// Command.
    #[must_use]
    pub const fn command(&self) -> RunnerCommandV2 {
        self.command
    }

    /// Family/mode/profile projection, absent only for `List`.
    #[must_use]
    pub const fn selection(&self) -> Option<&CommandSelectionV2> {
        self.selection.as_ref()
    }

    /// Contextually admitted budgets, absent only for `List`.
    #[must_use]
    pub const fn budgets(&self) -> Option<&AdmittedRunnerBudgetsV2> {
        self.budgets.as_ref()
    }

    /// Fixed disposition, absent only for `List`.
    #[must_use]
    pub const fn disposition(&self) -> Option<ArtifactDispositionV2> {
        self.disposition
    }

    /// Singular publication selection for durable commands.
    #[must_use]
    pub const fn publication_selection(&self) -> Option<&PublicationSelectionV2> {
        self.publication_selection.as_ref()
    }
}

fn provenance_name(provenance: &CommandSelectionProvenanceV2) -> &'static str {
    match provenance {
        CommandSelectionProvenanceV2::FixedPreflight(_) => "fixed-preflight",
        CommandSelectionProvenanceV2::FixedSelfTest(_) => "fixed-self-test",
        CommandSelectionProvenanceV2::CallerRun => "caller-run",
        CommandSelectionProvenanceV2::SealedNegative(_) => "sealed-negative",
        CommandSelectionProvenanceV2::SealedReplay(_) => "sealed-replay",
    }
}

#[cfg(test)]
mod tests {
    use super::{CommandIntentV2, CommandSelectionProvenanceV2, CommandSelectionV2};
    use crate::budget::{AdmittedRunnerBudgetsV2, RunnerBudgetsCandidateV2, RunnerBudgetsV2};
    use crate::catalog::{
        ArtifactDispositionV2, DestinationAdmissionModeV2, DigestRoleV2, LogicalUnitV2,
        PlatformPathProfileV2, PublicationProtocolV2, RunProfileV2, RunnerCommandV2,
    };
    use crate::identity::{CaseManifestRootV2, SourceIdentityRootV2};
    use crate::limits::RunnerLimitsV2;
    use crate::path::LogicalBundlePathV1;
    use crate::publication::{PublicationSelectionV2, PublicationTargetV2};
    use crate::value::StableTokenV2;

    fn token(value: &str) -> StableTokenV2 {
        StableTokenV2::new(value).expect("fixture token")
    }

    fn manifest(byte: u8) -> CaseManifestRootV2 {
        CaseManifestRootV2::parse_presented(
            DigestRoleV2::CaseManifest,
            CaseManifestRootV2::DESCRIPTOR.domain(),
            &format!("{byte:02x}").repeat(32),
        )
        .expect("fixture manifest")
    }

    fn source(byte: u8) -> SourceIdentityRootV2 {
        SourceIdentityRootV2::parse_presented(
            DigestRoleV2::Source,
            SourceIdentityRootV2::DESCRIPTOR.domain(),
            &format!("{byte:02x}").repeat(32),
        )
        .expect("fixture source")
    }

    fn selection_for(command: RunnerCommandV2, profile: RunProfileV2) -> CommandSelectionV2 {
        let family = token("family.fixture");
        let mode = token("mode.fixture");
        match command {
            RunnerCommandV2::Check => {
                CommandSelectionV2::fixed_preflight(manifest(1), family, mode, profile)
            }
            RunnerCommandV2::SelfTest => {
                CommandSelectionV2::fixed_self_test(manifest(2), family, mode, profile)
            }
            RunnerCommandV2::Run => CommandSelectionV2::caller_run(family, mode, profile),
            RunnerCommandV2::Negative => {
                CommandSelectionV2::sealed_negative(manifest(3), family, mode, profile)
            }
            RunnerCommandV2::Replay => {
                CommandSelectionV2::sealed_replay(source(4), family, mode, profile)
            }
            RunnerCommandV2::List => panic!("List has typed selector absence"),
        }
    }

    fn candidate(disposition: ArtifactDispositionV2) -> RunnerBudgetsCandidateV2 {
        let durable = disposition == ArtifactDispositionV2::DurableBundleRequired;
        RunnerBudgetsCandidateV2 {
            wall_time_ns: 100_000_000_000,
            max_resident_bytes: 1024 * 1024 * 1024,
            max_child_processes: u32::from(durable) * 8,
            max_parallel_children: u32::from(durable) * 4,
            logical_work_limit: 1000,
            logical_work_unit: LogicalUnitV2::Operations,
            lifecycle_encoded_bytes: 1000,
            command_result_stdout_bytes: 4000,
            combined_child_stdout_bytes: u64::from(durable) * 2000,
            combined_child_stderr_bytes: u64::from(durable) * 1000,
            artifact_encoded_bytes: u64::from(durable) * 100,
            artifact_stored_bytes: u64::from(durable) * 104,
            artifact_expanded_bytes: u64::from(durable) * 200,
            system_publication_stored_bytes: u64::from(durable) * 72,
            publication_stored_bytes: u64::from(durable) * 176,
            stop_observation_ns: 1_000_000_000,
            drain_ns: 1_000_000_000,
            finalize_ns: 1_000_000_000,
        }
    }

    fn budgets(
        profile: RunProfileV2,
        disposition: ArtifactDispositionV2,
    ) -> AdmittedRunnerBudgetsV2 {
        RunnerBudgetsV2::try_new(candidate(disposition))
            .expect("intrinsic fixture budgets")
            .admit(profile, disposition, &RunnerLimitsV2::base(profile))
            .expect("contextual fixture budgets")
    }

    fn publication() -> PublicationSelectionV2 {
        PublicationSelectionV2::new(
            PlatformPathProfileV2::PosixDescriptorRelativeV1,
            PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
            DestinationAdmissionModeV2::Absent,
            PublicationTargetV2::PosixRelative(
                LogicalBundlePathV1::new("results/bundle").expect("fixture path"),
            ),
        )
        .expect("fixture publication selection")
    }

    #[test]
    fn list_has_no_selectors_budgets_disposition_or_publication() {
        let list = CommandIntentV2::list();
        assert_eq!(list.command(), RunnerCommandV2::List);
        assert!(list.selection().is_none());
        assert!(list.budgets().is_none());
        assert!(list.disposition().is_none());
        assert!(list.publication_selection().is_none());
    }

    #[test]
    fn run_selection_has_no_default_mode_and_records_caller_provenance() {
        let selection = CommandSelectionV2::caller_run(
            token("family.fixture"),
            token("smoke"),
            RunProfileV2::Smoke,
        );
        assert_eq!(selection.family().as_str(), "family.fixture");
        assert_eq!(selection.mode().as_str(), "smoke");
        assert_eq!(selection.profile(), RunProfileV2::Smoke);
        assert!(matches!(
            selection.provenance(),
            CommandSelectionProvenanceV2::CallerRun
        ));
        assert!(StableTokenV2::new("").is_err());
    }

    #[test]
    fn exact_command_table_accepts_all_and_only_frozen_cells() {
        let commands = [
            RunnerCommandV2::Check,
            RunnerCommandV2::SelfTest,
            RunnerCommandV2::Run,
            RunnerCommandV2::Negative,
            RunnerCommandV2::Replay,
        ];
        for command in commands {
            for profile in [RunProfileV2::Smoke, RunProfileV2::Full] {
                let disposition = match command {
                    RunnerCommandV2::Check | RunnerCommandV2::SelfTest => {
                        ArtifactDispositionV2::LifecycleOnlyNoBundle
                    }
                    RunnerCommandV2::Run | RunnerCommandV2::Negative | RunnerCommandV2::Replay => {
                        ArtifactDispositionV2::DurableBundleRequired
                    }
                    RunnerCommandV2::List => unreachable!(),
                };
                let publication =
                    (disposition == ArtifactDispositionV2::DurableBundleRequired).then(publication);
                let intent = CommandIntentV2::new(
                    command,
                    selection_for(command, profile),
                    budgets(profile, disposition),
                    publication,
                )
                .expect("exact frozen command cell");
                assert_eq!(intent.command(), command);
                assert_eq!(
                    intent.selection().expect("present selection").profile(),
                    profile
                );
                assert_eq!(intent.disposition(), Some(disposition));
                assert_eq!(
                    intent.publication_selection().is_some(),
                    disposition == ArtifactDispositionV2::DurableBundleRequired
                );
            }
        }
    }

    #[test]
    fn every_cross_command_provenance_cell_refuses() {
        let commands = [
            RunnerCommandV2::Check,
            RunnerCommandV2::SelfTest,
            RunnerCommandV2::Run,
            RunnerCommandV2::Negative,
            RunnerCommandV2::Replay,
        ];
        for command in commands {
            let disposition = match command {
                RunnerCommandV2::Check | RunnerCommandV2::SelfTest => {
                    ArtifactDispositionV2::LifecycleOnlyNoBundle
                }
                _ => ArtifactDispositionV2::DurableBundleRequired,
            };
            for provenance_owner in commands {
                let result = CommandIntentV2::new(
                    command,
                    selection_for(provenance_owner, RunProfileV2::Smoke),
                    budgets(RunProfileV2::Smoke, disposition),
                    (disposition == ArtifactDispositionV2::DurableBundleRequired).then(publication),
                );
                assert_eq!(
                    result.is_ok(),
                    command == provenance_owner,
                    "{command:?} with {provenance_owner:?} provenance"
                );
            }
        }
    }

    #[test]
    fn profile_disposition_and_publication_presence_cannot_drift() {
        assert!(
            CommandIntentV2::new(
                RunnerCommandV2::Run,
                selection_for(RunnerCommandV2::Run, RunProfileV2::Smoke),
                budgets(
                    RunProfileV2::Full,
                    ArtifactDispositionV2::DurableBundleRequired
                ),
                Some(publication()),
            )
            .is_err()
        );
        assert!(
            CommandIntentV2::new(
                RunnerCommandV2::Run,
                selection_for(RunnerCommandV2::Run, RunProfileV2::Smoke),
                budgets(
                    RunProfileV2::Smoke,
                    ArtifactDispositionV2::LifecycleOnlyNoBundle
                ),
                Some(publication()),
            )
            .is_err()
        );
        assert!(
            CommandIntentV2::new(
                RunnerCommandV2::Check,
                selection_for(RunnerCommandV2::Check, RunProfileV2::Smoke),
                budgets(
                    RunProfileV2::Smoke,
                    ArtifactDispositionV2::LifecycleOnlyNoBundle
                ),
                Some(publication()),
            )
            .is_err()
        );
        assert!(
            CommandIntentV2::new(
                RunnerCommandV2::Replay,
                selection_for(RunnerCommandV2::Replay, RunProfileV2::Smoke),
                budgets(
                    RunProfileV2::Smoke,
                    ArtifactDispositionV2::DurableBundleRequired
                ),
                None,
            )
            .is_err()
        );
        assert!(
            CommandIntentV2::new(
                RunnerCommandV2::List,
                selection_for(RunnerCommandV2::Run, RunProfileV2::Smoke),
                budgets(
                    RunProfileV2::Smoke,
                    ArtifactDispositionV2::DurableBundleRequired
                ),
                Some(publication()),
            )
            .is_err()
        );
    }
}
