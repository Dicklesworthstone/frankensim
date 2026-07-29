//! Frozen command-selection provenance and applicability table.

use crate::budget::AdmittedRunnerBudgetsV2;
use crate::catalog::{
    ArtifactDispositionV2, DiagnosticCodeV2, ProofExitV2, RunProfileV2, RunnerCommandV2,
};
use crate::construction::{ConstructionErrorKindV2, ConstructionErrorV2};
use crate::identity::{CaseManifestRootV2, SourceIdentityRootV2};
use crate::publication::PublicationSelectionV2;
use crate::value::StableTokenV2;

/// Cardinality of one caller-presented command selector.
///
/// This is a non-wire, pure-validation vocabulary. `Duplicate` means the same
/// selector value was repeated; `Ambiguous` means multiple distinct values
/// were presented for a singular selector.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CommandSelectorCardinalityV2 {
    /// No value was presented.
    Absent,
    /// Exactly one value was presented.
    Singular,
    /// One value was presented more than once.
    Duplicate,
    /// Multiple distinct values were presented.
    Ambiguous,
}

impl CommandSelectorCardinalityV2 {
    /// Exact closed cardinality vocabulary.
    pub const ALL: [Self; 4] = [
        Self::Absent,
        Self::Singular,
        Self::Duplicate,
        Self::Ambiguous,
    ];
}

/// One caller-selectable command dimension in deterministic validation order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CommandSelectorFieldV2 {
    /// Registered family selector.
    Family,
    /// Registered family-mode selector.
    Mode,
    /// Smoke-or-Full profile selector.
    Profile,
    /// Immutable negative-case selector.
    NegativeCase,
    /// Immutable replay-source selector.
    ReplaySource,
}

impl CommandSelectorFieldV2 {
    /// Frozen first-error validation order.
    pub const ALL: [Self; 5] = [
        Self::Family,
        Self::Mode,
        Self::Profile,
        Self::NegativeCase,
        Self::ReplaySource,
    ];

    /// Stable field name used by diagnostics and logs.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Family => "command.selector.family",
            Self::Mode => "command.selector.mode",
            Self::Profile => "command.selector.profile",
            Self::NegativeCase => "command.selector.negative_case",
            Self::ReplaySource => "command.selector.replay_source",
        }
    }
}

/// Exact applicability requirement for one caller selector.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CommandSelectorExpectationV2 {
    /// The command forbids the caller selector.
    Absent,
    /// The command requires exactly one caller selector value.
    Singular,
}

/// Exact Usage-class boundary crossed by caller selector validation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CommandSelectorUsageKindV2 {
    /// A required singular selector was omitted.
    Missing,
    /// One selector value was repeated.
    Duplicate,
    /// Multiple distinct values competed for a singular selector.
    Ambiguous,
    /// A selector was presented to a command that does not accept it.
    Inapplicable,
}

/// Caller selector cardinalities before sealed-manifest or source projection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CommandSelectorPresenceV2 {
    family: CommandSelectorCardinalityV2,
    mode: CommandSelectorCardinalityV2,
    profile: CommandSelectorCardinalityV2,
    negative_case: CommandSelectorCardinalityV2,
    replay_source: CommandSelectorCardinalityV2,
}

impl CommandSelectorPresenceV2 {
    /// Construct one complete caller-selector cardinality vector.
    #[must_use]
    pub const fn new(
        family: CommandSelectorCardinalityV2,
        mode: CommandSelectorCardinalityV2,
        profile: CommandSelectorCardinalityV2,
        negative_case: CommandSelectorCardinalityV2,
        replay_source: CommandSelectorCardinalityV2,
    ) -> Self {
        Self {
            family,
            mode,
            profile,
            negative_case,
            replay_source,
        }
    }

    /// Exact valid caller-selector vector for a command.
    ///
    /// `Check` and `SelfTest` obtain their internal family/mode/profile from
    /// sealed manifests, while `Negative` and `Replay` obtain those values
    /// from their selected immutable case/source. Those derived values are not
    /// caller selectors and therefore remain absent here.
    #[must_use]
    pub const fn exact_for(command: RunnerCommandV2) -> Self {
        match command {
            RunnerCommandV2::List | RunnerCommandV2::Check | RunnerCommandV2::SelfTest => {
                Self::new(
                    CommandSelectorCardinalityV2::Absent,
                    CommandSelectorCardinalityV2::Absent,
                    CommandSelectorCardinalityV2::Absent,
                    CommandSelectorCardinalityV2::Absent,
                    CommandSelectorCardinalityV2::Absent,
                )
            }
            RunnerCommandV2::Run => Self::new(
                CommandSelectorCardinalityV2::Singular,
                CommandSelectorCardinalityV2::Singular,
                CommandSelectorCardinalityV2::Singular,
                CommandSelectorCardinalityV2::Absent,
                CommandSelectorCardinalityV2::Absent,
            ),
            RunnerCommandV2::Negative => Self::new(
                CommandSelectorCardinalityV2::Absent,
                CommandSelectorCardinalityV2::Absent,
                CommandSelectorCardinalityV2::Absent,
                CommandSelectorCardinalityV2::Singular,
                CommandSelectorCardinalityV2::Absent,
            ),
            RunnerCommandV2::Replay => Self::new(
                CommandSelectorCardinalityV2::Absent,
                CommandSelectorCardinalityV2::Absent,
                CommandSelectorCardinalityV2::Absent,
                CommandSelectorCardinalityV2::Absent,
                CommandSelectorCardinalityV2::Singular,
            ),
        }
    }

    /// Cardinality of one named selector.
    #[must_use]
    pub const fn cardinality(&self, field: CommandSelectorFieldV2) -> CommandSelectorCardinalityV2 {
        match field {
            CommandSelectorFieldV2::Family => self.family,
            CommandSelectorFieldV2::Mode => self.mode,
            CommandSelectorFieldV2::Profile => self.profile,
            CommandSelectorFieldV2::NegativeCase => self.negative_case,
            CommandSelectorFieldV2::ReplaySource => self.replay_source,
        }
    }
}

/// Deterministic Usage-class refusal for a caller-selector vector.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CommandSelectorUsageV2 {
    kind: CommandSelectorUsageKindV2,
    field: CommandSelectorFieldV2,
    expected: CommandSelectorExpectationV2,
    observed: CommandSelectorCardinalityV2,
}

impl CommandSelectorUsageV2 {
    /// Exact Usage boundary.
    #[must_use]
    pub const fn proof_exit(&self) -> ProofExitV2 {
        ProofExitV2::Usage
    }

    /// Exact base diagnostic required for this Usage boundary.
    #[must_use]
    pub const fn diagnostic_code(&self) -> DiagnosticCodeV2 {
        DiagnosticCodeV2::RunnerUsage
    }

    /// Stable refusal class.
    #[must_use]
    pub const fn kind(&self) -> CommandSelectorUsageKindV2 {
        self.kind
    }

    /// First invalid selector in the frozen validation order.
    #[must_use]
    pub const fn field(&self) -> CommandSelectorFieldV2 {
        self.field
    }

    /// Exact command-specific applicability requirement.
    #[must_use]
    pub const fn expected(&self) -> CommandSelectorExpectationV2 {
        self.expected
    }

    /// Exact presented cardinality.
    #[must_use]
    pub const fn observed(&self) -> CommandSelectorCardinalityV2 {
        self.observed
    }

    /// Stable diagnostic owner.
    #[must_use]
    pub const fn owner(&self) -> &'static str {
        "fs-evidence-runner.command-selectors"
    }
}

/// Validate the complete caller-selector applicability vector for one command.
///
/// This function performs no token parsing, manifest lookup, defaulting, or
/// I/O. A successful `Negative` or `Replay` presence check still requires the
/// later sealed projection represented by [`CommandSelectionV2`].
pub fn validate_command_selector_presence_v2(
    command: RunnerCommandV2,
    presented: CommandSelectorPresenceV2,
) -> Result<(), CommandSelectorUsageV2> {
    let exact = CommandSelectorPresenceV2::exact_for(command);
    for field in CommandSelectorFieldV2::ALL {
        let expected = match exact.cardinality(field) {
            CommandSelectorCardinalityV2::Absent => CommandSelectorExpectationV2::Absent,
            CommandSelectorCardinalityV2::Singular => CommandSelectorExpectationV2::Singular,
            CommandSelectorCardinalityV2::Duplicate | CommandSelectorCardinalityV2::Ambiguous => {
                unreachable!("the frozen command table contains only absent or singular cells")
            }
        };
        let observed = presented.cardinality(field);
        let valid = matches!(
            (expected, observed),
            (
                CommandSelectorExpectationV2::Absent,
                CommandSelectorCardinalityV2::Absent
            ) | (
                CommandSelectorExpectationV2::Singular,
                CommandSelectorCardinalityV2::Singular
            )
        );
        if valid {
            continue;
        }
        let kind = match observed {
            CommandSelectorCardinalityV2::Absent => CommandSelectorUsageKindV2::Missing,
            CommandSelectorCardinalityV2::Singular => CommandSelectorUsageKindV2::Inapplicable,
            CommandSelectorCardinalityV2::Duplicate => CommandSelectorUsageKindV2::Duplicate,
            CommandSelectorCardinalityV2::Ambiguous => CommandSelectorUsageKindV2::Ambiguous,
        };
        return Err(CommandSelectorUsageV2 {
            kind,
            field,
            expected,
            observed,
        });
    }
    Ok(())
}

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
/// The safe `List` constructor produces typed absence for every
/// command-inapplicable field:
///
/// ```
/// use fs_evidence_runner::{CommandIntentV2, RunnerCommandV2};
///
/// let intent = CommandIntentV2::list();
/// assert_eq!(intent.command(), RunnerCommandV2::List);
/// assert!(intent.selection().is_none());
/// assert!(intent.budgets().is_none());
/// assert!(intent.disposition().is_none());
/// assert!(intent.publication_selection().is_none());
/// ```
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
///
/// A validated command intent cannot be converted into an authority scope:
///
/// ```compile_fail
/// use fs_evidence_runner::CommandIntentV2;
/// use fs_evidence_runner::identity::AuthorityScopeRootV2;
///
/// let intent = CommandIntentV2::list();
/// let _authority: AuthorityScopeRootV2 = intent.into();
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
    use super::{
        CommandIntentV2, CommandSelectionProvenanceV2, CommandSelectionV2,
        CommandSelectorCardinalityV2, CommandSelectorExpectationV2, CommandSelectorFieldV2,
        CommandSelectorPresenceV2, CommandSelectorUsageKindV2,
        validate_command_selector_presence_v2,
    };
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
    #[allow(
        clippy::too_many_lines,
        reason = "the test intentionally preserves the complete command-by-selector-by-publication acceptance table as one literal oracle"
    )]
    fn exact_command_table_accepts_all_and_only_frozen_cells() {
        let selector_commands = [
            RunnerCommandV2::List,
            RunnerCommandV2::Check,
            RunnerCommandV2::SelfTest,
            RunnerCommandV2::Run,
            RunnerCommandV2::Negative,
            RunnerCommandV2::Replay,
        ];
        for command in selector_commands {
            let required = match command {
                RunnerCommandV2::List | RunnerCommandV2::Check | RunnerCommandV2::SelfTest => {
                    [false, false, false, false, false]
                }
                RunnerCommandV2::Run => [true, true, true, false, false],
                RunnerCommandV2::Negative => [false, false, false, true, false],
                RunnerCommandV2::Replay => [false, false, false, false, true],
            };
            let mut accepted = 0_u32;
            for ordinal in 0..4_usize.pow(5) {
                let mut remainder = ordinal;
                let mut cells = [CommandSelectorCardinalityV2::Absent; 5];
                for cell in &mut cells {
                    *cell = CommandSelectorCardinalityV2::ALL[remainder % 4];
                    remainder /= 4;
                }
                let presented = CommandSelectorPresenceV2::new(
                    cells[0], cells[1], cells[2], cells[3], cells[4],
                );
                let expected_error = CommandSelectorFieldV2::ALL
                    .into_iter()
                    .zip(required)
                    .zip(cells)
                    .find_map(|((field, required), observed)| {
                        let expected = if required {
                            CommandSelectorExpectationV2::Singular
                        } else {
                            CommandSelectorExpectationV2::Absent
                        };
                        let exact = if required {
                            CommandSelectorCardinalityV2::Singular
                        } else {
                            CommandSelectorCardinalityV2::Absent
                        };
                        (observed != exact).then(|| {
                            let kind = match observed {
                                CommandSelectorCardinalityV2::Absent => {
                                    CommandSelectorUsageKindV2::Missing
                                }
                                CommandSelectorCardinalityV2::Singular => {
                                    CommandSelectorUsageKindV2::Inapplicable
                                }
                                CommandSelectorCardinalityV2::Duplicate => {
                                    CommandSelectorUsageKindV2::Duplicate
                                }
                                CommandSelectorCardinalityV2::Ambiguous => {
                                    CommandSelectorUsageKindV2::Ambiguous
                                }
                            };
                            (kind, field, expected, observed)
                        })
                    });
                match (
                    validate_command_selector_presence_v2(command, presented),
                    expected_error,
                ) {
                    (Ok(()), None) => accepted += 1,
                    (Err(error), Some((kind, field, expected, observed))) => {
                        assert_eq!(error.proof_exit(), crate::catalog::ProofExitV2::Usage);
                        assert_eq!(
                            error.diagnostic_code(),
                            crate::catalog::DiagnosticCodeV2::RunnerUsage
                        );
                        assert_eq!(error.kind(), kind);
                        assert_eq!(error.field(), field);
                        assert_eq!(error.field().name(), field.name());
                        assert_eq!(error.expected(), expected);
                        assert_eq!(error.observed(), observed);
                        assert_eq!(error.owner(), "fs-evidence-runner.command-selectors");
                    }
                    (actual, expected) => panic!(
                        "selector Cartesian mismatch for {command:?}, ordinal {ordinal}: \
                         actual={actual:?}, expected={expected:?}"
                    ),
                }
            }
            assert_eq!(
                accepted, 1,
                "exactly one of all 4^5 selector-cardinality cells is valid for {command:?}"
            );
        }

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
