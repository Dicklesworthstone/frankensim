//! Zero-dependency constellation admission protocol.
//!
//! This module is path-included by the standalone bootstrap, so it deliberately
//! uses only `std` and remains compatible with Rust 2021. It owns the policy
//! state machine, not operating-system authority: fixed-size identifiers name
//! capabilities established by an entrypoint, but never authenticate or
//! recreate those capabilities. In particular, decoding produces an inert
//! [`RecordedAdmission`], never a live [`AdmissionMachine`].

use std::fmt;

/// Canonical binary envelope magic.
pub const ADMISSION_MAGIC: &[u8; 8] = b"FSADM\0\0\0";
/// Canonical schema name retained by human and ledger adapters.
pub const ADMISSION_SCHEMA: &str = "frankensim-constellation-admission-v3";
/// Domain separating this protocol from every other identity surface.
pub const ADMISSION_DOMAIN: &str = "org.frankensim.xtask.constellation-admission.v3";
/// Canonical schema version.
pub const ADMISSION_VERSION: u16 = 3;
/// Maximum path-capability slots retained by one request.
pub const MAX_PATH_CAPABILITIES: usize = 5;
/// Maximum executable-capability slots retained by one request.
pub const MAX_EXECUTABLE_CAPABILITIES: usize = 2;
/// Maximum accepted transition count, including retries and terminal steps.
pub const MAX_ADMISSION_EVENTS: usize = 1_024;
/// Hard envelope bound checked before decoder allocation.
pub const MAX_ADMISSION_BYTES: usize = 262_144;

const ID_BYTES: usize = 32;

/// Opaque, externally established identity carried without secret material.
///
/// The admission core never derives this value from a path, URL, environment,
/// or weak local hash. Adapters must bind a strong identity in their own trust
/// domain and pass only its redacted bytes here.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AuthorityId([u8; ID_BYTES]);

impl AuthorityId {
    /// Admit a nonzero opaque identity.
    ///
    /// # Errors
    /// Refuses the all-zero sentinel, which cannot name live authority.
    pub fn try_from_bytes(bytes: [u8; ID_BYTES]) -> Result<Self, AdmissionError> {
        if bytes == [0; ID_BYTES] {
            return Err(AdmissionError::new(AdmissionRule::ZeroIdentity));
        }
        Ok(Self(bytes))
    }

    /// Exact redacted identity bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; ID_BYTES] {
        &self.0
    }
}

impl fmt::Debug for AuthorityId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "AuthorityId({:02x}{:02x}{:02x}{:02x}..)",
            self.0[0], self.0[1], self.0[2], self.0[3]
        )
    }
}

/// Closed command classes; no free-form command can inherit authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
pub enum CommandClass {
    /// Diagnostic inspection with no mutation or publication authority.
    Diagnostic = 0,
    /// Read-only verification of an already materialized constellation.
    VerifyOnly = 1,
    /// Read-only stable constellation snapshot.
    Snapshot = 2,
    /// Pinned-source materialization and verification.
    Bootstrap = 3,
    /// Deliberate constellation-lock publication.
    LockConstellation = 4,
    /// Repository quality proof under DSR.
    DsrQuality = 5,
    /// Release build under DSR.
    DsrBuild = 6,
    /// Narrow remote Cargo probe under RCH.
    RchProbe = 7,
    /// Shell adapter for read-only verification.
    ShellVerifyOnly = 8,
    /// Shell adapter for a read-only stable snapshot.
    ShellSnapshot = 9,
    /// Shell adapter for pinned-source bootstrap.
    ShellBootstrap = 10,
}

impl CommandClass {
    fn from_tag(tag: u8) -> Option<Self> {
        Some(match tag {
            0 => Self::Diagnostic,
            1 => Self::VerifyOnly,
            2 => Self::Snapshot,
            3 => Self::Bootstrap,
            4 => Self::LockConstellation,
            5 => Self::DsrQuality,
            6 => Self::DsrBuild,
            7 => Self::RchProbe,
            8 => Self::ShellVerifyOnly,
            9 => Self::ShellSnapshot,
            10 => Self::ShellBootstrap,
            _ => return None,
        })
    }

    fn required_paths(self) -> &'static [PathSlot] {
        match self {
            Self::Diagnostic => &[PathSlot::WorkspaceRoot],
            Self::VerifyOnly | Self::Snapshot | Self::ShellVerifyOnly | Self::ShellSnapshot => {
                &[PathSlot::WorkspaceRoot, PathSlot::ConstellationLock]
            }
            Self::Bootstrap | Self::ShellBootstrap => &[
                PathSlot::WorkspaceRoot,
                PathSlot::ConstellationLock,
                PathSlot::DestinationRoot,
            ],
            Self::LockConstellation => &[PathSlot::WorkspaceRoot, PathSlot::ConstellationLock],
            Self::DsrQuality | Self::DsrBuild | Self::RchProbe => &[
                PathSlot::WorkspaceRoot,
                PathSlot::ConstellationLock,
                PathSlot::OutputRoot,
            ],
        }
    }

    fn required_executables(self) -> &'static [ExecutableSlot] {
        match self {
            Self::Diagnostic => &[],
            Self::VerifyOnly | Self::Snapshot | Self::Bootstrap | Self::LockConstellation => {
                &[ExecutableSlot::Git]
            }
            Self::ShellVerifyOnly | Self::ShellSnapshot | Self::ShellBootstrap => {
                &[ExecutableSlot::Git, ExecutableSlot::Shell]
            }
            Self::DsrQuality | Self::DsrBuild => &[ExecutableSlot::Git, ExecutableSlot::Dsr],
            Self::RchProbe => &[ExecutableSlot::Git, ExecutableSlot::Rch],
        }
    }

    fn requires_offline_fetch(self) -> bool {
        matches!(
            self,
            Self::Diagnostic
                | Self::VerifyOnly
                | Self::Snapshot
                | Self::LockConstellation
                | Self::ShellVerifyOnly
                | Self::ShellSnapshot
        )
    }

    fn admits_publication(self, publication: PublicationAuthority) -> bool {
        matches!(
            (self, publication),
            (
                Self::Diagnostic
                    | Self::VerifyOnly
                    | Self::Snapshot
                    | Self::ShellVerifyOnly
                    | Self::ShellSnapshot,
                PublicationAuthority::Prohibited
            ) | (
                Self::Bootstrap,
                PublicationAuthority::BootstrapReceipt { .. }
            ) | (
                Self::ShellBootstrap,
                PublicationAuthority::BootstrapReceipt { .. }
            ) | (
                Self::LockConstellation,
                PublicationAuthority::ConstellationLock { .. }
            ) | (
                Self::DsrQuality | Self::DsrBuild | Self::RchProbe,
                PublicationAuthority::ProofReceipt { .. }
            )
        )
    }

    fn path_is_allowed(self, slot: PathSlot, publication: PublicationAuthority) -> bool {
        self.required_paths().contains(&slot)
            || matches!(
                (self, slot),
                (
                    Self::Bootstrap | Self::ShellBootstrap,
                    PathSlot::SourceCache
                )
            )
            || (slot == PathSlot::PublicationTarget && !publication.is_prohibited())
    }

    fn executable_is_allowed(self, slot: ExecutableSlot) -> bool {
        self.required_executables().contains(&slot)
    }

    fn admits_read_only_completion(self) -> bool {
        matches!(
            self,
            Self::Diagnostic
                | Self::VerifyOnly
                | Self::Snapshot
                | Self::ShellVerifyOnly
                | Self::ShellSnapshot
        )
    }
}

/// Explicit authority to contact a pinned transport.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FetchAuthority {
    /// No network operation is authorized, independent of the numeric budget.
    Offline,
    /// Only a transport already pinned by admitted lock data may be contacted.
    PinnedTransport {
        /// Identity of the live network capability held by the adapter.
        capability: AuthorityId,
    },
}

/// Explicit publication authority; absence never grants a sink.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PublicationAuthority {
    /// Publication is prohibited.
    Prohibited,
    /// Bootstrap provenance may be published to the named path capability.
    BootstrapReceipt {
        /// Identity of the publication-target capability.
        capability: AuthorityId,
    },
    /// A deliberate lock operation may publish the constellation lock.
    ConstellationLock {
        /// Identity of the publication-target capability.
        capability: AuthorityId,
    },
    /// A proof runner may publish one bounded proof receipt.
    ProofReceipt {
        /// Identity of the publication-target capability.
        capability: AuthorityId,
    },
}

impl PublicationAuthority {
    fn is_prohibited(self) -> bool {
        matches!(self, Self::Prohibited)
    }
}

/// Typed evidence that one authorized publication transaction succeeded.
///
/// The adapter establishes the external observation, but the admission core
/// binds that observation to the exact request, attempt, publication target,
/// conditional fence, and receipt that were authorized before the effect. A
/// payload from another request or publication transaction cannot finalize the
/// current machine. Every identity is nonzero because [`AuthorityId`] has no
/// public zero-valued constructor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PublicationSuccessEvidence {
    request_identity: AuthorityId,
    attempt: u16,
    publication_capability: AuthorityId,
    fence: AuthorityId,
    receipt: AuthorityId,
    observation: AuthorityId,
}

impl PublicationSuccessEvidence {
    /// Bind an external success observation to one authorized transaction.
    #[must_use]
    pub const fn new(
        request_identity: AuthorityId,
        attempt: u16,
        publication_capability: AuthorityId,
        fence: AuthorityId,
        receipt: AuthorityId,
        observation: AuthorityId,
    ) -> Self {
        Self {
            request_identity,
            attempt,
            publication_capability,
            fence,
            receipt,
            observation,
        }
    }

    /// Request whose publication was observed.
    #[must_use]
    pub const fn request_identity(self) -> AuthorityId {
        self.request_identity
    }

    /// Attempt whose publication was observed.
    #[must_use]
    pub const fn attempt(self) -> u16 {
        self.attempt
    }

    /// Publication-target capability consumed by the external commit.
    #[must_use]
    pub const fn publication_capability(self) -> AuthorityId {
        self.publication_capability
    }

    /// Conditional fence consumed by the external commit.
    #[must_use]
    pub const fn fence(self) -> AuthorityId {
        self.fence
    }

    /// Exact authorized content/receipt identity that became durable.
    #[must_use]
    pub const fn receipt(self) -> AuthorityId {
        self.receipt
    }

    /// Strong identity of the external success/durability observation.
    #[must_use]
    pub const fn observation(self) -> AuthorityId {
        self.observation
    }
}

/// Closed path capability roles.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
pub enum PathSlot {
    /// Canonical FrankenSim checkout root.
    WorkspaceRoot = 0,
    /// Exact constellation lock input.
    ConstellationLock = 1,
    /// Parent directory that receives or already contains siblings.
    DestinationRoot = 2,
    /// Optional admitted source-cache root.
    SourceCache = 3,
    /// Directory for proof logs and build artifacts.
    OutputRoot = 4,
    /// Exact eventual publication target.
    PublicationTarget = 5,
}

impl PathSlot {
    fn from_tag(tag: u8) -> Option<Self> {
        Some(match tag {
            0 => Self::WorkspaceRoot,
            1 => Self::ConstellationLock,
            2 => Self::DestinationRoot,
            3 => Self::SourceCache,
            4 => Self::OutputRoot,
            5 => Self::PublicationTarget,
            _ => return None,
        })
    }
}

/// One path role bound to an externally established capability identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PathCapability {
    slot: PathSlot,
    identity: AuthorityId,
}

impl PathCapability {
    /// Bind a closed path role to a redacted identity.
    #[must_use]
    pub const fn new(slot: PathSlot, identity: AuthorityId) -> Self {
        Self { slot, identity }
    }

    /// Bound role.
    #[must_use]
    pub const fn slot(&self) -> PathSlot {
        self.slot
    }

    /// Externally established capability identity.
    #[must_use]
    pub const fn identity(&self) -> AuthorityId {
        self.identity
    }
}

/// Closed executable capability roles.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
pub enum ExecutableSlot {
    /// Scrubbed, pinned Git executable.
    Git = 0,
    /// Hermetic shell adapter executable.
    Shell = 1,
    /// DSR runner executable.
    Dsr = 2,
    /// RCH runner executable.
    Rch = 3,
}

impl ExecutableSlot {
    fn from_tag(tag: u8) -> Option<Self> {
        Some(match tag {
            0 => Self::Git,
            1 => Self::Shell,
            2 => Self::Dsr,
            3 => Self::Rch,
            _ => return None,
        })
    }
}

/// One executable role bound to an externally established capability identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExecutableCapability {
    slot: ExecutableSlot,
    identity: AuthorityId,
}

impl ExecutableCapability {
    /// Bind a closed executable role to a redacted identity.
    #[must_use]
    pub const fn new(slot: ExecutableSlot, identity: AuthorityId) -> Self {
        Self { slot, identity }
    }

    /// Bound role.
    #[must_use]
    pub const fn slot(&self) -> ExecutableSlot {
        self.slot
    }

    /// Externally established capability identity.
    #[must_use]
    pub const fn identity(&self) -> AuthorityId {
        self.identity
    }
}

/// Explicit asupersync/cancellation and monotonic-clock binding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CxBinding {
    cx: AuthorityId,
    cancellation: AuthorityId,
    clock: AuthorityId,
    max_unpolled_work: u64,
}

impl CxBinding {
    /// Bind runtime context, cancellation source, clock, and poll granularity.
    ///
    /// # Errors
    /// Refuses a zero poll interval.
    pub fn try_new(
        cx: AuthorityId,
        cancellation: AuthorityId,
        clock: AuthorityId,
        max_unpolled_work: u64,
    ) -> Result<Self, AdmissionError> {
        if max_unpolled_work == 0 {
            return Err(AdmissionError::new(AdmissionRule::ZeroPollInterval));
        }
        if cx == cancellation || cx == clock || cancellation == clock {
            return Err(AdmissionError::new(AdmissionRule::CxIdentityAliasing));
        }
        Ok(Self {
            cx,
            cancellation,
            clock,
            max_unpolled_work,
        })
    }

    /// Runtime context identity.
    #[must_use]
    pub const fn cx(&self) -> AuthorityId {
        self.cx
    }

    /// Cancellation-source identity.
    #[must_use]
    pub const fn cancellation(&self) -> AuthorityId {
        self.cancellation
    }

    /// Monotonic-clock identity.
    #[must_use]
    pub const fn clock(&self) -> AuthorityId {
        self.clock
    }

    /// Maximum work units permitted between cancellation polls.
    #[must_use]
    pub const fn max_unpolled_work(&self) -> u64 {
        self.max_unpolled_work
    }
}

/// Deadline in caller-supplied ticks from one explicitly bound clock.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DeadlineBudget {
    clock: AuthorityId,
    not_after_tick: u64,
}

impl DeadlineBudget {
    /// Construct an inclusive monotonic deadline.
    #[must_use]
    pub const fn new(clock: AuthorityId, not_after_tick: u64) -> Self {
        Self {
            clock,
            not_after_tick,
        }
    }

    /// Bound clock identity.
    #[must_use]
    pub const fn clock(&self) -> AuthorityId {
        self.clock
    }

    /// Inclusive last admitted tick.
    #[must_use]
    pub const fn not_after_tick(&self) -> u64 {
        self.not_after_tick
    }
}

/// Work and memory caps.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ComputeBudget {
    /// Total logical work units across retries.
    pub work_units: u64,
    /// Total admitted memory-byte reservations across retries.
    pub memory_bytes: u64,
}

/// Process, file, and output caps.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IoBudget {
    /// Total child-process starts.
    pub processes: u32,
    /// Total file opens or creations.
    pub files: u32,
    /// Total external stdout, stderr, provenance, and proof-output bytes.
    ///
    /// The admission protocol's own canonical audit envelope is separately
    /// bounded by [`MAX_ADMISSION_BYTES`] and is not charged as operation output.
    pub output_bytes: u64,
}

/// Network caps, independent of network authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NetworkBudget {
    /// Total network requests.
    pub requests: u32,
    /// Total admitted network bytes.
    pub bytes: u64,
}

/// Complete explicit resource budget.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdmissionBudgets {
    deadline: DeadlineBudget,
    compute: ComputeBudget,
    io: IoBudget,
    network: NetworkBudget,
    retries: u16,
}

impl AdmissionBudgets {
    /// Construct one complete budget; there is deliberately no `Default`.
    #[must_use]
    pub const fn new(
        deadline: DeadlineBudget,
        compute: ComputeBudget,
        io: IoBudget,
        network: NetworkBudget,
        retries: u16,
    ) -> Self {
        Self {
            deadline,
            compute,
            io,
            network,
            retries,
        }
    }

    /// Deadline budget.
    #[must_use]
    pub const fn deadline(&self) -> DeadlineBudget {
        self.deadline
    }

    /// Compute budgets.
    #[must_use]
    pub const fn compute(&self) -> ComputeBudget {
        self.compute
    }

    /// Process/file/output budgets.
    #[must_use]
    pub const fn io(&self) -> IoBudget {
        self.io
    }

    /// Numeric network budgets, which do not grant network authority.
    #[must_use]
    pub const fn network(&self) -> NetworkBudget {
        self.network
    }

    /// Maximum additional attempts after the initial attempt.
    #[must_use]
    pub const fn retries(&self) -> u16 {
        self.retries
    }
}

/// Explicit trust-anchor disposition at request construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrustAnchorState {
    /// No live trust anchor is bound; only an unanchored record can result.
    Unanchored,
    /// Exact expected anchor and generation.
    Anchored {
        /// Strong identity of the expected trust anchor.
        identity: AuthorityId,
        /// Expected immutable generation.
        generation: u64,
    },
}

/// Borrowed construction input checked before the core allocates capability storage.
#[derive(Debug)]
pub struct AdmissionContextSpec<'a> {
    /// Stable idempotency identity retained across retries.
    pub request_identity: AuthorityId,
    /// Closed command class.
    pub command: CommandClass,
    /// Explicit offline or pinned-fetch authority.
    pub fetch: FetchAuthority,
    /// Explicit publication disposition.
    pub publication: PublicationAuthority,
    /// Runtime cancellation and clock binding.
    pub cx: CxBinding,
    /// Complete budgets.
    pub budgets: AdmissionBudgets,
    /// Expected trust-anchor disposition.
    pub trust_anchor: TrustAnchorState,
    /// Path capability bindings.
    pub path_capabilities: &'a [PathCapability],
    /// Executable capability bindings.
    pub executable_capabilities: &'a [ExecutableCapability],
}

/// Validated, canonical request policy. This is not a live OS capability.
#[derive(Debug, PartialEq, Eq)]
pub struct AdmissionContext {
    request_identity: AuthorityId,
    command: CommandClass,
    fetch: FetchAuthority,
    publication: PublicationAuthority,
    cx: CxBinding,
    budgets: AdmissionBudgets,
    trust_anchor: TrustAnchorState,
    path_capabilities: Vec<PathCapability>,
    executable_capabilities: Vec<ExecutableCapability>,
}

impl AdmissionContext {
    /// Validate and canonicalize a complete context.
    ///
    /// Counts and fixed command requirements are checked before capability
    /// vectors are allocated. Slot order is canonicalized and duplicate roles
    /// refuse rather than overwrite.
    ///
    /// # Errors
    /// Returns a deterministic structured refusal for every missing, duplicate,
    /// inconsistent, or over-limit declaration.
    pub fn try_new(spec: AdmissionContextSpec<'_>) -> Result<Self, AdmissionError> {
        if spec.path_capabilities.len() > MAX_PATH_CAPABILITIES {
            return Err(AdmissionError::new(AdmissionRule::TooManyPathCapabilities));
        }
        if spec.executable_capabilities.len() > MAX_EXECUTABLE_CAPABILITIES {
            return Err(AdmissionError::new(
                AdmissionRule::TooManyExecutableCapabilities,
            ));
        }
        if spec.budgets.deadline.clock != spec.cx.clock {
            return Err(AdmissionError::new(AdmissionRule::ClockBindingMismatch));
        }
        match spec.fetch {
            FetchAuthority::Offline => {
                if spec.budgets.network
                    != (NetworkBudget {
                        requests: 0,
                        bytes: 0,
                    })
                {
                    return Err(AdmissionError::new(
                        AdmissionRule::OfflineNetworkBudgetConflict,
                    ));
                }
            }
            FetchAuthority::PinnedTransport { .. } => {
                if spec.budgets.network.requests == 0 || spec.budgets.network.bytes == 0 {
                    return Err(AdmissionError::new(
                        AdmissionRule::FetchAuthorityWithoutBudget,
                    ));
                }
            }
        }
        if spec.command.requires_offline_fetch() && !matches!(spec.fetch, FetchAuthority::Offline) {
            return Err(AdmissionError::new(
                AdmissionRule::NetworkForbiddenForCommand,
            ));
        }
        if !spec.command.admits_publication(spec.publication) {
            let rule = if spec.publication.is_prohibited()
                && !spec.command.admits_read_only_completion()
            {
                AdmissionRule::PublicationRequiredForCommand
            } else {
                AdmissionRule::PublicationForbiddenForCommand
            };
            return Err(AdmissionError::new(rule));
        }
        if !spec.publication.is_prohibited() && spec.budgets.io.output_bytes == 0 {
            return Err(AdmissionError::new(
                AdmissionRule::PublicationWithoutOutputBudget,
            ));
        }

        validate_path_capabilities(spec.command, spec.publication, spec.path_capabilities)?;
        validate_executable_capabilities(spec.command, spec.executable_capabilities)?;
        if let Some(publication_identity) = publication_identity(spec.publication) {
            let target_identity = spec
                .path_capabilities
                .iter()
                .find(|capability| capability.slot == PathSlot::PublicationTarget)
                .map(|capability| capability.identity);
            if target_identity != Some(publication_identity) {
                return Err(AdmissionError::new(
                    AdmissionRule::PublicationCapabilityMismatch,
                ));
            }
        }

        let mut path_capabilities = Vec::new();
        path_capabilities
            .try_reserve_exact(spec.path_capabilities.len())
            .map_err(|_| AdmissionError::new(AdmissionRule::AllocationFailed))?;
        path_capabilities.extend_from_slice(spec.path_capabilities);
        path_capabilities.sort_by_key(|capability| capability.slot);

        let mut executable_capabilities = Vec::new();
        executable_capabilities
            .try_reserve_exact(spec.executable_capabilities.len())
            .map_err(|_| AdmissionError::new(AdmissionRule::AllocationFailed))?;
        executable_capabilities.extend_from_slice(spec.executable_capabilities);
        executable_capabilities.sort_by_key(|capability| capability.slot);

        Ok(Self {
            request_identity: spec.request_identity,
            command: spec.command,
            fetch: spec.fetch,
            publication: spec.publication,
            cx: spec.cx,
            budgets: spec.budgets,
            trust_anchor: spec.trust_anchor,
            path_capabilities,
            executable_capabilities,
        })
    }

    /// Stable idempotency identity.
    #[must_use]
    pub const fn request_identity(&self) -> AuthorityId {
        self.request_identity
    }

    /// Closed command class.
    #[must_use]
    pub const fn command(&self) -> CommandClass {
        self.command
    }

    /// Fetch authority.
    #[must_use]
    pub const fn fetch(&self) -> FetchAuthority {
        self.fetch
    }

    /// Publication disposition.
    #[must_use]
    pub const fn publication(&self) -> PublicationAuthority {
        self.publication
    }

    /// Runtime binding.
    #[must_use]
    pub const fn cx(&self) -> CxBinding {
        self.cx
    }

    /// Complete resource budgets.
    #[must_use]
    pub const fn budgets(&self) -> AdmissionBudgets {
        self.budgets
    }

    /// Trust-anchor state.
    #[must_use]
    pub const fn trust_anchor(&self) -> TrustAnchorState {
        self.trust_anchor
    }

    /// Canonically ordered path capability bindings.
    #[must_use]
    pub fn path_capabilities(&self) -> &[PathCapability] {
        &self.path_capabilities
    }

    /// Canonically ordered executable capability bindings.
    #[must_use]
    pub fn executable_capabilities(&self) -> &[ExecutableCapability] {
        &self.executable_capabilities
    }
}

const fn publication_identity(publication: PublicationAuthority) -> Option<AuthorityId> {
    match publication {
        PublicationAuthority::Prohibited => None,
        PublicationAuthority::BootstrapReceipt { capability }
        | PublicationAuthority::ConstellationLock { capability }
        | PublicationAuthority::ProofReceipt { capability } => Some(capability),
    }
}

fn validate_path_capabilities(
    command: CommandClass,
    publication: PublicationAuthority,
    capabilities: &[PathCapability],
) -> Result<(), AdmissionError> {
    if capabilities.iter().enumerate().any(|(index, left)| {
        capabilities[index + 1..]
            .iter()
            .any(|right| right.slot == left.slot)
    }) {
        return Err(AdmissionError::new(AdmissionRule::DuplicatePathCapability));
    }
    if capabilities
        .iter()
        .any(|capability| !command.path_is_allowed(capability.slot, publication))
    {
        return Err(AdmissionError::new(AdmissionRule::UnexpectedPathCapability));
    }
    for required in command.required_paths() {
        if !capabilities
            .iter()
            .any(|capability| capability.slot == *required)
        {
            return Err(AdmissionError::new(AdmissionRule::MissingPathCapability));
        }
    }
    if !publication.is_prohibited()
        && !capabilities
            .iter()
            .any(|capability| capability.slot == PathSlot::PublicationTarget)
    {
        return Err(AdmissionError::new(AdmissionRule::MissingPathCapability));
    }
    Ok(())
}

fn validate_executable_capabilities(
    command: CommandClass,
    capabilities: &[ExecutableCapability],
) -> Result<(), AdmissionError> {
    if capabilities.iter().enumerate().any(|(index, left)| {
        capabilities[index + 1..]
            .iter()
            .any(|right| right.slot == left.slot)
    }) {
        return Err(AdmissionError::new(
            AdmissionRule::DuplicateExecutableCapability,
        ));
    }
    if capabilities
        .iter()
        .any(|capability| !command.executable_is_allowed(capability.slot))
    {
        return Err(AdmissionError::new(
            AdmissionRule::UnexpectedExecutableCapability,
        ));
    }
    for required in command.required_executables() {
        if !capabilities
            .iter()
            .any(|capability| capability.slot == *required)
        {
            return Err(AdmissionError::new(
                AdmissionRule::MissingExecutableCapability,
            ));
        }
    }
    Ok(())
}

/// Canonical admission states. These are truth values, not display strings.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
pub enum StateKind {
    /// Context exists but has not completed stable preflight.
    Diagnostic = 0,
    /// Preflight had no live trust anchor and minted no authority.
    Unanchored = 1,
    /// Stable preflight admitted bounded work or publication.
    Admitted = 2,
    /// A definitive pre-effect failure established no mutation and no authority.
    Refused = 3,
    /// Cancellation completed request, drain, and explicit finalization.
    Cancelled = 4,
    /// Effects or cleanup could not be proven, so no success/refusal claim exists.
    Indeterminate = 5,
}

/// Diagnostic/preflight phase.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiagnosticPhase {
    /// Context admitted structurally; no external observation is bound.
    Created,
    /// First complete observation is bound and awaits a stability recheck.
    Preflighted {
        /// Exact observed source snapshot.
        snapshot: AuthorityId,
        /// Exact observed trust anchor.
        anchor: AuthorityId,
        /// Exact observed anchor generation.
        generation: u64,
    },
}

/// Admitted phase; publication still requires a second stability barrier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdmittedPhase {
    /// Bounded work may be reserved; no publication permit exists.
    Active {
        /// Snapshot admitted by the preflight recheck.
        snapshot: AuthorityId,
        /// Trust anchor admitted by the preflight recheck.
        anchor: AuthorityId,
        /// Stable admitted generation.
        generation: u64,
    },
    /// Read-only work is quiescent and awaits an exact post-work recheck.
    CompletionPending {
        /// Snapshot admitted before work.
        snapshot: AuthorityId,
        /// Trust anchor admitted before work.
        anchor: AuthorityId,
        /// Stable admitted generation.
        generation: u64,
        /// Evidence that read-only workers and retained outputs are quiescent.
        quiescence: AuthorityId,
    },
    /// Publication-prohibited work completed under a stable post-work recheck.
    Completed {
        /// Inert identity of the finalized read-only result receipt.
        receipt: AuthorityId,
    },
    /// Work is closed and publication awaits a final stability recheck.
    PublicationPending {
        /// Snapshot admitted before work.
        snapshot: AuthorityId,
        /// Trust anchor admitted before work.
        anchor: AuthorityId,
        /// Stable admitted generation.
        generation: u64,
        /// Evidence that effectful workers and staged outputs were quiesced.
        quiescence: AuthorityId,
    },
    /// Final recheck succeeded; one receipt may be authorized for publication.
    PublicationReady {
        /// Rechecked snapshot.
        snapshot: AuthorityId,
        /// Rechecked anchor.
        anchor: AuthorityId,
        /// Rechecked generation.
        generation: u64,
        /// Conditional publication fence established by the recheck adapter.
        fence: AuthorityId,
        /// Evidence that work was quiesced before the recheck.
        quiescence: AuthorityId,
    },
    /// A single-use publication identity is authorized for an external commit.
    PublicationCommitting {
        /// Conditional publication fence that the sink must consume.
        fence: AuthorityId,
        /// Exact content/receipt identity authorized for publication.
        receipt: AuthorityId,
        /// Evidence that work was quiesced before authorization.
        quiescence: AuthorityId,
    },
    /// Publication finalized with one exact receipt and success observation.
    Published {
        /// Inert identity of the finalized publication receipt.
        receipt: AuthorityId,
        /// Transaction-bound evidence of external commit and durability success.
        success: PublicationSuccessEvidence,
    },
}

/// Explicit terminal-record finalization phase.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TerminalPhase {
    /// The terminal fact exists but its audit record is not finalized.
    Pending,
    /// The terminal audit record is finalized.
    Finalized,
}

/// Cancellation phase; only `Cancelled { phase: Finalized, .. }` is a completed
/// cancellation claim.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CancellationPhase {
    /// New work is closed and drain obligations are fixed.
    Requested,
    /// At least one deterministic drain observation occurred.
    Draining,
    /// All obligations drained and finalization completed.
    Finalized,
}

/// Public, inert state view.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdmissionState {
    /// Diagnostic state.
    Diagnostic(DiagnosticPhase),
    /// Unanchored terminal fact.
    Unanchored(TerminalPhase),
    /// Admitted work/publication state.
    Admitted(AdmittedPhase),
    /// Definitive refusal and its record phase.
    Refused {
        /// Rule that established definitive refusal.
        rule: AdmissionRule,
        /// Audit-record finalization phase.
        phase: TerminalPhase,
    },
    /// Cancellation progress and remaining drain obligations.
    Cancelled {
        /// Deterministic cancellation cause.
        cause: CancellationCause,
        /// Request/drain/finalization phase.
        phase: CancellationPhase,
        /// Remaining child/process/file/output obligations.
        remaining: DrainObligations,
    },
    /// Uncertain terminal fact and its record phase.
    Indeterminate {
        /// Rule that prevented a stronger terminal claim.
        rule: AdmissionRule,
        /// Audit-record finalization phase.
        phase: TerminalPhase,
    },
}

impl AdmissionState {
    /// Coarse canonical state kind.
    #[must_use]
    pub const fn kind(self) -> StateKind {
        match self {
            Self::Diagnostic(_) => StateKind::Diagnostic,
            Self::Unanchored(_) => StateKind::Unanchored,
            Self::Admitted(_) => StateKind::Admitted,
            Self::Refused { .. } => StateKind::Refused,
            Self::Cancelled { .. } => StateKind::Cancelled,
            Self::Indeterminate { .. } => StateKind::Indeterminate,
        }
    }

    /// Whether the current attempt has a finalized terminal record.
    ///
    /// A retryable terminal may still begin a new attempt with a fresh child Cx;
    /// this method makes no statement about request-level retry legality.
    #[must_use]
    pub const fn attempt_is_finalized(self) -> bool {
        matches!(
            self,
            Self::Unanchored(TerminalPhase::Finalized)
                | Self::Refused {
                    phase: TerminalPhase::Finalized,
                    ..
                }
                | Self::Cancelled {
                    phase: CancellationPhase::Finalized,
                    ..
                }
                | Self::Indeterminate {
                    phase: TerminalPhase::Finalized,
                    ..
                }
                | Self::Admitted(AdmittedPhase::Completed { .. } | AdmittedPhase::Published { .. })
        )
    }
}

/// Trust-anchor observation made by an adapter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnchorObservation {
    /// No live anchor could be established.
    Unavailable,
    /// Exact observed anchor and immutable generation.
    Observed {
        /// Observed anchor identity.
        identity: AuthorityId,
        /// Observed generation.
        generation: u64,
    },
}

/// Resource dimension charged transactionally before effectful work.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BudgetCharge {
    /// Logical work units.
    Work(u64),
    /// Memory-byte reservation.
    Memory(u64),
    /// Child-process start count.
    Processes(u32),
    /// File-open or file-create count.
    Files(u32),
    /// Retained output bytes.
    Output(u64),
    /// Network request and byte reservation.
    Network {
        /// Request count.
        requests: u32,
        /// Byte count.
        bytes: u64,
    },
}

/// Complete cumulative consumption across every attempt.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BudgetConsumption {
    /// Logical work units consumed.
    pub work_units: u64,
    /// Memory-byte reservations consumed.
    pub memory_bytes: u64,
    /// Child-process starts consumed.
    pub processes: u32,
    /// File opens/creates consumed.
    pub files: u32,
    /// Retained output bytes consumed.
    pub output_bytes: u64,
    /// Network requests consumed.
    pub network_requests: u32,
    /// Network bytes consumed.
    pub network_bytes: u64,
    /// Additional attempts consumed.
    pub retries: u16,
}

/// Explicit drain obligations fixed when cancellation closes new work.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct DrainObligations {
    /// Live child processes that must exit.
    pub processes: u32,
    /// Open or staged files that must be finalized or closed.
    pub files: u32,
    /// Pending bounded output fragments that must be drained.
    pub outputs: u32,
}

impl DrainObligations {
    const fn is_empty(self) -> bool {
        self.processes == 0 && self.files == 0 && self.outputs == 0
    }

    fn checked_sub(self, progress: Self) -> Option<Self> {
        Some(Self {
            processes: self.processes.checked_sub(progress.processes)?,
            files: self.files.checked_sub(progress.files)?,
            outputs: self.outputs.checked_sub(progress.outputs)?,
        })
    }
}

/// Deterministic cancellation cause; no unbounded free-form text is retained.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum CancellationCause {
    /// Caller explicitly requested cancellation.
    Requested = 0,
    /// Inclusive deadline elapsed.
    Deadline = 1,
    /// Parent scope requested cancellation.
    ParentScope = 2,
    /// A bounded fault injector requested cancellation.
    InjectedFault = 3,
}

impl CancellationCause {
    fn from_tag(tag: u8) -> Option<Self> {
        Some(match tag {
            0 => Self::Requested,
            1 => Self::Deadline,
            2 => Self::ParentScope,
            3 => Self::InjectedFault,
            _ => return None,
        })
    }
}

/// Closed reasons for a definitive refusal before any uncertain effect exists.
///
/// These tags are deliberately disjoint from indeterminate and publication
/// failure tags. An adapter cannot inject an arbitrary internal codec or state
/// machine rule into a canonical refusal receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
pub enum RefusalReason {
    /// A bounded external input failed its adapter contract.
    InputRejected = 0,
    /// A source artifact failed validation before mutation.
    SourceRejected = 1,
    /// A required live capability could not be established.
    CapabilityRejected = 2,
    /// A required executable failed pre-effect admission.
    ExecutableRejected = 3,
    /// Declared work cannot fit the explicit budget envelope.
    BudgetInfeasible = 4,
    /// The inclusive deadline elapsed before an effect began.
    DeadlineElapsed = 5,
    /// Closed command policy denied the requested pre-effect operation.
    PolicyDenied = 6,
}

impl RefusalReason {
    const fn from_tag(tag: u8) -> Option<Self> {
        Some(match tag {
            0 => Self::InputRejected,
            1 => Self::SourceRejected,
            2 => Self::CapabilityRejected,
            3 => Self::ExecutableRejected,
            4 => Self::BudgetInfeasible,
            5 => Self::DeadlineElapsed,
            6 => Self::PolicyDenied,
            _ => return None,
        })
    }

    /// Stable terminal rule exposed by states and events.
    #[must_use]
    pub const fn rule(self) -> AdmissionRule {
        match self {
            Self::InputRejected => AdmissionRule::InputRejected,
            Self::SourceRejected => AdmissionRule::SourceRejected,
            Self::CapabilityRejected => AdmissionRule::CapabilityRejected,
            Self::ExecutableRejected => AdmissionRule::ExecutableRejected,
            Self::BudgetInfeasible => AdmissionRule::BudgetInfeasible,
            Self::DeadlineElapsed => AdmissionRule::DeadlineExceeded,
            Self::PolicyDenied => AdmissionRule::PolicyDenied,
        }
    }
}

/// Closed reasons for uncertainty after an effect may have occurred.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
pub enum IndeterminateReason {
    /// The adapter cannot prove whether an external effect occurred.
    EffectOutcomeUnknown = 16,
    /// Known drain obligations did not complete.
    DrainIncomplete = 17,
    /// The adapter could not obtain a trustworthy drain observation.
    DrainObservationUnavailable = 18,
    /// Terminal-record finalization may have partially occurred.
    FinalizationOutcomeUnknown = 19,
    /// Source or anchor moved after admitted work began.
    PostEffectStabilityChanged = 20,
}

impl IndeterminateReason {
    const fn from_tag(tag: u8) -> Option<Self> {
        Some(match tag {
            16 => Self::EffectOutcomeUnknown,
            17 => Self::DrainIncomplete,
            18 => Self::DrainObservationUnavailable,
            19 => Self::FinalizationOutcomeUnknown,
            20 => Self::PostEffectStabilityChanged,
            _ => return None,
        })
    }

    /// Stable terminal rule exposed by states and events.
    #[must_use]
    pub const fn rule(self) -> AdmissionRule {
        match self {
            Self::EffectOutcomeUnknown => AdmissionRule::EffectOutcomeUnknown,
            Self::DrainIncomplete => AdmissionRule::DrainIncomplete,
            Self::DrainObservationUnavailable => AdmissionRule::DrainObservationUnavailable,
            Self::FinalizationOutcomeUnknown => AdmissionRule::FinalizationOutcomeUnknown,
            Self::PostEffectStabilityChanged => AdmissionRule::PostEffectStabilityChanged,
        }
    }
}

/// Closed reasons for uncertainty after a publication receipt was authorized.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
pub enum PublicationFailureReason {
    /// The external commit reported failure after authorization.
    CommitFailed = 32,
    /// The adapter cannot determine whether the commit occurred.
    CommitOutcomeUnknown = 33,
    /// Commit visibility is known but durability is not.
    DurabilityUnknown = 34,
    /// Durable content exists but receipt finalization is uncertain.
    ReceiptFinalizationUnknown = 35,
}

impl PublicationFailureReason {
    const fn from_tag(tag: u8) -> Option<Self> {
        Some(match tag {
            32 => Self::CommitFailed,
            33 => Self::CommitOutcomeUnknown,
            34 => Self::DurabilityUnknown,
            35 => Self::ReceiptFinalizationUnknown,
            _ => return None,
        })
    }

    /// Stable terminal rule exposed by states and events.
    #[must_use]
    pub const fn rule(self) -> AdmissionRule {
        match self {
            Self::CommitFailed => AdmissionRule::PublicationCommitFailed,
            Self::CommitOutcomeUnknown => AdmissionRule::PublicationOutcomeUnknown,
            Self::DurabilityUnknown => AdmissionRule::PublicationDurabilityUnknown,
            Self::ReceiptFinalizationUnknown => AdmissionRule::PublicationFinalizationUnknown,
        }
    }
}

/// Runtime/replay action classes used by the exhaustive transition table.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
pub enum TransitionKind {
    /// Bind the first complete observation.
    Preflight = 0,
    /// Re-observe the exact preflight snapshot before admitting work.
    StabilityRecheck = 1,
    /// Charge one budget dimension before work or allocation.
    Charge = 2,
    /// Establish a definitive pre-effect refusal.
    Refuse = 3,
    /// Close new work and fix drain obligations.
    RequestCancellation = 4,
    /// Record deterministic drain progress.
    Drain = 5,
    /// Finalize an unanchored, refused, cancelled, or indeterminate record.
    FinalizeTerminal = 6,
    /// Declare uncertainty after possible effect or failed cleanup.
    DeclareIndeterminate = 7,
    /// Close work and request publication recheck.
    BeginPublication = 8,
    /// Recheck the stable snapshot after work and before publication.
    PublicationRecheck = 9,
    /// Finalize exactly one publication receipt.
    FinalizePublication = 10,
    /// Start another bounded attempt under the same request identity.
    Retry = 11,
    /// Record one explicit cancellation poll boundary.
    PollCancellation = 12,
    /// Authorize one receipt identity against the rechecked publication fence.
    AuthorizePublication = 13,
    /// Record an external publication failure as indeterminate.
    PublicationFailed = 14,
    /// Close publication-prohibited work before its post-work stability check.
    BeginReadOnlyCompletion = 15,
    /// Recheck and finalize one publication-prohibited result receipt.
    CompleteReadOnly = 16,
}

impl TransitionKind {
    fn from_tag(tag: u8) -> Option<Self> {
        Some(match tag {
            0 => Self::Preflight,
            1 => Self::StabilityRecheck,
            2 => Self::Charge,
            3 => Self::Refuse,
            4 => Self::RequestCancellation,
            5 => Self::Drain,
            6 => Self::FinalizeTerminal,
            7 => Self::DeclareIndeterminate,
            8 => Self::BeginPublication,
            9 => Self::PublicationRecheck,
            10 => Self::FinalizePublication,
            11 => Self::Retry,
            12 => Self::PollCancellation,
            13 => Self::AuthorizePublication,
            14 => Self::PublicationFailed,
            15 => Self::BeginReadOnlyCompletion,
            16 => Self::CompleteReadOnly,
            _ => return None,
        })
    }
}

/// Complete deterministic transition payload.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdmissionTransition {
    /// First complete source/anchor observation.
    Preflight {
        /// Complete source snapshot identity.
        snapshot: AuthorityId,
        /// Anchor observation.
        anchor: AnchorObservation,
        /// Caller-supplied monotonic tick.
        at_tick: u64,
    },
    /// Stability recheck of the preflight observation.
    StabilityRecheck {
        /// Reobserved source snapshot.
        snapshot: AuthorityId,
        /// Reobserved anchor.
        anchor: AnchorObservation,
        /// Caller-supplied monotonic tick.
        at_tick: u64,
    },
    /// Transactional resource charge.
    Charge {
        /// Resource reservation.
        charge: BudgetCharge,
        /// Caller-supplied monotonic tick.
        at_tick: u64,
    },
    /// Definitive refusal before an uncertain effect exists.
    Refuse {
        /// Closed deterministic pre-effect reason.
        reason: RefusalReason,
        /// Identity of the bounded observation establishing the refusal.
        observation: AuthorityId,
    },
    /// Close new work and fix drain obligations.
    RequestCancellation {
        /// Deterministic cause.
        cause: CancellationCause,
        /// Complete outstanding obligations.
        obligations: DrainObligations,
        /// Identity of the complete liveness observation fixing those obligations.
        observation: AuthorityId,
    },
    /// Deterministic drain progress.
    Drain {
        /// Obligations proven complete by this step.
        completed: DrainObligations,
        /// Identity of the drain observation proving this progress.
        observation: AuthorityId,
    },
    /// Finalize a terminal audit record.
    FinalizeTerminal {
        /// Identity of the explicit finalization receipt.
        receipt: AuthorityId,
    },
    /// Declare that mutation or cleanup cannot be proven.
    DeclareIndeterminate {
        /// Closed deterministic post-effect uncertainty reason.
        reason: IndeterminateReason,
        /// Identity of the observation establishing uncertainty.
        observation: AuthorityId,
    },
    /// Close work before the publication stability barrier.
    BeginPublication {
        /// Evidence that workers and staged outputs are quiescent.
        quiescence: AuthorityId,
        /// Caller-supplied monotonic tick.
        at_tick: u64,
    },
    /// Reobserve source and anchor immediately before publication.
    PublicationRecheck {
        /// Reobserved source snapshot.
        snapshot: AuthorityId,
        /// Reobserved anchor.
        anchor: AnchorObservation,
        /// Conditional fence that the publication sink must consume.
        fence: AuthorityId,
        /// Caller-supplied monotonic tick.
        at_tick: u64,
    },
    /// Authorize one exact receipt identity before the external commit begins.
    AuthorizePublication {
        /// Strong identity of the content/receipt to publish.
        receipt: AuthorityId,
        /// Caller-supplied monotonic tick.
        at_tick: u64,
    },
    /// Finalize the already-authorized receipt after a successful external commit.
    FinalizePublication {
        /// Nonzero observation bound to the exact authorized transaction.
        success: PublicationSuccessEvidence,
    },
    /// Preserve uncertainty when the authorized external commit fails or is ambiguous.
    PublicationFailed {
        /// Closed publication failure or uncertainty reason.
        reason: PublicationFailureReason,
        /// Identity of the external commit/durability observation.
        observation: AuthorityId,
    },
    /// Close publication-prohibited work before its final stability check.
    BeginReadOnlyCompletion {
        /// Evidence that read-only workers and retained outputs are quiescent.
        quiescence: AuthorityId,
        /// Caller-supplied monotonic tick.
        at_tick: u64,
    },
    /// Reobserve and finalize a publication-prohibited result.
    CompleteReadOnly {
        /// Reobserved source snapshot.
        snapshot: AuthorityId,
        /// Reobserved trust anchor.
        anchor: AnchorObservation,
        /// Exact result receipt identity.
        receipt: AuthorityId,
        /// Caller-supplied monotonic tick.
        at_tick: u64,
    },
    /// Retry under the same idempotency identity without replenishing budgets.
    Retry {
        /// Fresh attempt-scoped child Cx and cancellation binding.
        next_cx: CxBinding,
        /// Caller-supplied monotonic tick.
        at_tick: u64,
    },
    /// Record one bounded cancellation poll boundary.
    PollCancellation {
        /// Logical work completed since the preceding poll.
        work_since_last_poll: u64,
        /// Caller-supplied monotonic tick.
        at_tick: u64,
    },
}

impl AdmissionTransition {
    /// Closed action kind.
    #[must_use]
    pub const fn kind(self) -> TransitionKind {
        match self {
            Self::Preflight { .. } => TransitionKind::Preflight,
            Self::StabilityRecheck { .. } => TransitionKind::StabilityRecheck,
            Self::Charge { .. } => TransitionKind::Charge,
            Self::Refuse { .. } => TransitionKind::Refuse,
            Self::RequestCancellation { .. } => TransitionKind::RequestCancellation,
            Self::Drain { .. } => TransitionKind::Drain,
            Self::FinalizeTerminal { .. } => TransitionKind::FinalizeTerminal,
            Self::DeclareIndeterminate { .. } => TransitionKind::DeclareIndeterminate,
            Self::BeginPublication { .. } => TransitionKind::BeginPublication,
            Self::PublicationRecheck { .. } => TransitionKind::PublicationRecheck,
            Self::AuthorizePublication { .. } => TransitionKind::AuthorizePublication,
            Self::FinalizePublication { .. } => TransitionKind::FinalizePublication,
            Self::PublicationFailed { .. } => TransitionKind::PublicationFailed,
            Self::BeginReadOnlyCompletion { .. } => TransitionKind::BeginReadOnlyCompletion,
            Self::CompleteReadOnly { .. } => TransitionKind::CompleteReadOnly,
            Self::Retry { .. } => TransitionKind::Retry,
            Self::PollCancellation { .. } => TransitionKind::PollCancellation,
        }
    }

    fn observed_tick(self) -> Option<u64> {
        match self {
            Self::Preflight { at_tick, .. }
            | Self::StabilityRecheck { at_tick, .. }
            | Self::Charge { at_tick, .. }
            | Self::BeginPublication { at_tick, .. }
            | Self::PublicationRecheck { at_tick, .. }
            | Self::AuthorizePublication { at_tick, .. }
            | Self::BeginReadOnlyCompletion { at_tick, .. }
            | Self::CompleteReadOnly { at_tick, .. }
            | Self::Retry { at_tick, .. }
            | Self::PollCancellation { at_tick, .. } => Some(at_tick),
            Self::Refuse { .. }
            | Self::RequestCancellation { .. }
            | Self::Drain { .. }
            | Self::FinalizeTerminal { .. }
            | Self::DeclareIndeterminate { .. }
            | Self::FinalizePublication { .. }
            | Self::PublicationFailed { .. } => None,
        }
    }
}

/// Pure coarse transition-matrix candidate check.
///
/// This deliberately over-approximates phase-specific legality; [`AdmissionMachine::apply`]
/// performs the exact phase guard. The table exists so every coarse state/action
/// permutation has a deterministic answer rather than a missing match arm. An
/// adapter must never treat `true` here as authority to perform an effect.
#[must_use]
pub const fn transition_kind_may_apply(from: StateKind, action: TransitionKind) -> bool {
    match from {
        StateKind::Diagnostic => matches!(
            action,
            TransitionKind::Preflight
                | TransitionKind::StabilityRecheck
                | TransitionKind::Charge
                | TransitionKind::PollCancellation
                | TransitionKind::Refuse
        ),
        StateKind::Unanchored => matches!(
            action,
            TransitionKind::FinalizeTerminal | TransitionKind::Retry
        ),
        StateKind::Admitted => matches!(
            action,
            TransitionKind::Charge
                | TransitionKind::PollCancellation
                | TransitionKind::RequestCancellation
                | TransitionKind::DeclareIndeterminate
                | TransitionKind::BeginPublication
                | TransitionKind::BeginReadOnlyCompletion
                | TransitionKind::CompleteReadOnly
                | TransitionKind::PublicationRecheck
                | TransitionKind::AuthorizePublication
                | TransitionKind::FinalizePublication
                | TransitionKind::PublicationFailed
        ),
        StateKind::Refused => {
            matches!(
                action,
                TransitionKind::FinalizeTerminal | TransitionKind::Retry
            )
        }
        StateKind::Indeterminate => matches!(action, TransitionKind::FinalizeTerminal),
        StateKind::Cancelled => matches!(
            action,
            TransitionKind::Drain
                | TransitionKind::FinalizeTerminal
                | TransitionKind::DeclareIndeterminate
                | TransitionKind::Retry
        ),
    }
}

/// One canonical event emitted only after a transition commits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdmissionEvent {
    sequence: u16,
    attempt: u16,
    from: StateKind,
    action: TransitionKind,
    to: StateKind,
    terminal_rule: Option<AdmissionRule>,
    consumption: BudgetConsumption,
    observed_tick: Option<u64>,
}

impl AdmissionEvent {
    /// Zero-based deterministic sequence.
    #[must_use]
    pub const fn sequence(&self) -> u16 {
        self.sequence
    }

    /// Zero-based attempt ordinal.
    #[must_use]
    pub const fn attempt(&self) -> u16 {
        self.attempt
    }

    /// Source state.
    #[must_use]
    pub const fn from(&self) -> StateKind {
        self.from
    }

    /// Action kind.
    #[must_use]
    pub const fn action(&self) -> TransitionKind {
        self.action
    }

    /// Destination state.
    #[must_use]
    pub const fn to(&self) -> StateKind {
        self.to
    }

    /// Rule when the legal observation deterministically produced a terminal fact.
    #[must_use]
    pub const fn terminal_rule(&self) -> Option<AdmissionRule> {
        self.terminal_rule
    }

    /// Complete consumption after the event.
    #[must_use]
    pub const fn consumption(&self) -> BudgetConsumption {
        self.consumption
    }

    /// Monotonic clock observation carried by this event, when applicable.
    #[must_use]
    pub const fn observed_tick(&self) -> Option<u64> {
        self.observed_tick
    }
}

/// Live deterministic policy machine.
///
/// It is intentionally neither `Clone` nor deserializable. Adapters retain the
/// actual OS/process/cancellation handles and must require this machine's state
/// before performing effects.
#[derive(Debug)]
pub struct AdmissionMachine {
    context: AdmissionContext,
    current_cx: CxBinding,
    state: AdmissionState,
    consumption: BudgetConsumption,
    attempt: u16,
    last_observed_tick: Option<u64>,
    history: Vec<AdmissionTransition>,
    events: Vec<AdmissionEvent>,
}

impl AdmissionMachine {
    /// Create a diagnostic machine and preallocate its bounded audit storage.
    ///
    /// # Errors
    /// Refuses if bounded event/history storage cannot be reserved.
    pub fn try_new(context: AdmissionContext) -> Result<Self, AdmissionError> {
        let mut history = Vec::new();
        history
            .try_reserve_exact(MAX_ADMISSION_EVENTS)
            .map_err(|_| AdmissionError::new(AdmissionRule::AllocationFailed))?;
        let mut events = Vec::new();
        events
            .try_reserve_exact(MAX_ADMISSION_EVENTS)
            .map_err(|_| AdmissionError::new(AdmissionRule::AllocationFailed))?;
        let current_cx = context.cx;
        Ok(Self {
            context,
            current_cx,
            state: AdmissionState::Diagnostic(DiagnosticPhase::Created),
            consumption: BudgetConsumption::default(),
            attempt: 0,
            last_observed_tick: None,
            history,
            events,
        })
    }

    /// Immutable validated request context.
    #[must_use]
    pub const fn context(&self) -> &AdmissionContext {
        &self.context
    }

    /// Attempt-scoped live Cx binding; retries replace it with a fresh child.
    #[must_use]
    pub const fn current_cx(&self) -> CxBinding {
        self.current_cx
    }

    /// Current inert state view.
    #[must_use]
    pub const fn state(&self) -> AdmissionState {
        self.state
    }

    /// Complete cumulative resource consumption.
    #[must_use]
    pub const fn consumption(&self) -> BudgetConsumption {
        self.consumption
    }

    /// Current zero-based attempt ordinal.
    #[must_use]
    pub const fn attempt(&self) -> u16 {
        self.attempt
    }

    /// Last committed monotonic clock observation across every attempt.
    #[must_use]
    pub const fn last_observed_tick(&self) -> Option<u64> {
        self.last_observed_tick
    }

    /// Deterministically ordered committed events.
    #[must_use]
    pub fn events(&self) -> &[AdmissionEvent] {
        &self.events
    }

    /// Exact committed transition history, including evidence identities.
    #[must_use]
    pub fn history(&self) -> &[AdmissionTransition] {
        &self.history
    }

    /// Apply one transition transactionally.
    ///
    /// All guards, checked arithmetic, deadline, and event-cap checks complete
    /// before any field is changed. On error the machine is byte-for-byte
    /// semantically unchanged and no event is appended.
    ///
    /// # Errors
    /// Returns a bounded structured refusal for an illegal transition, exceeded
    /// budget/deadline, denied authority, incomplete drain, or event cap.
    pub fn apply(
        &mut self,
        transition: AdmissionTransition,
    ) -> Result<AdmissionEvent, AdmissionError> {
        if self.history.len() >= MAX_ADMISSION_EVENTS {
            return Err(AdmissionError::new(AdmissionRule::EventLimitExceeded));
        }
        if let Some(tick) = transition.observed_tick() {
            if self
                .last_observed_tick
                .is_some_and(|last_observed_tick| tick < last_observed_tick)
            {
                return Err(AdmissionError::new(AdmissionRule::ClockRegression));
            }
            if tick > self.context.budgets.deadline.not_after_tick {
                return Err(AdmissionError::new(AdmissionRule::DeadlineExceeded));
            }
        }
        if matches!(
            self.state,
            AdmissionState::Admitted(
                AdmittedPhase::Completed { .. } | AdmittedPhase::Published { .. }
            )
        ) {
            return Err(AdmissionError::new(AdmissionRule::IllegalTransition));
        }
        if matches!(
            (self.state, transition),
            (
                AdmissionState::Indeterminate {
                    phase: TerminalPhase::Finalized,
                    ..
                },
                AdmissionTransition::Retry { .. }
            )
        ) {
            return Err(AdmissionError::new(
                AdmissionRule::IndeterminateRetryForbidden,
            ));
        }
        if !transition_kind_may_apply(self.state.kind(), transition.kind()) {
            return Err(AdmissionError::new(AdmissionRule::IllegalTransition));
        }

        let from = self.state.kind();
        let evaluated = self.evaluate(transition)?;
        let events_after = self
            .history
            .len()
            .checked_add(1)
            .ok_or_else(|| AdmissionError::new(AdmissionRule::EventLimitExceeded))?;
        let remaining_capacity = MAX_ADMISSION_EVENTS
            .checked_sub(events_after)
            .ok_or_else(|| AdmissionError::new(AdmissionRule::EventLimitExceeded))?;
        if remaining_capacity < minimum_terminal_events(evaluated.state) {
            return Err(AdmissionError::new(AdmissionRule::TerminalHeadroomRequired));
        }
        let sequence = u16::try_from(self.events.len())
            .map_err(|_| AdmissionError::new(AdmissionRule::EventLimitExceeded))?;
        let event = AdmissionEvent {
            sequence,
            attempt: evaluated.attempt,
            from,
            action: transition.kind(),
            to: evaluated.state.kind(),
            terminal_rule: evaluated.terminal_rule,
            consumption: evaluated.consumption,
            observed_tick: transition.observed_tick(),
        };

        self.state = evaluated.state;
        self.consumption = evaluated.consumption;
        self.attempt = evaluated.attempt;
        self.current_cx = evaluated.cx;
        if let Some(tick) = transition.observed_tick() {
            self.last_observed_tick = Some(tick);
        }
        self.history.push(transition);
        self.events.push(event);
        Ok(event)
    }

    fn evaluate(&self, transition: AdmissionTransition) -> Result<Evaluated, AdmissionError> {
        let mut evaluated = Evaluated {
            state: self.state,
            consumption: self.consumption,
            attempt: self.attempt,
            cx: self.current_cx,
            terminal_rule: None,
        };
        match (self.state, transition) {
            (
                AdmissionState::Diagnostic(DiagnosticPhase::Created),
                AdmissionTransition::Preflight {
                    snapshot, anchor, ..
                },
            ) => match (self.context.trust_anchor, anchor) {
                (TrustAnchorState::Unanchored, _) => {
                    evaluated.state = AdmissionState::Unanchored(TerminalPhase::Pending);
                    evaluated.terminal_rule = Some(AdmissionRule::TrustAnchorUnavailable);
                }
                (TrustAnchorState::Anchored { .. }, AnchorObservation::Unavailable) => {
                    evaluated.state = refused(AdmissionRule::TrustAnchorUnavailable);
                    evaluated.terminal_rule = Some(AdmissionRule::TrustAnchorUnavailable);
                }
                (
                    TrustAnchorState::Anchored {
                        identity: expected_identity,
                        generation: expected_generation,
                    },
                    AnchorObservation::Observed {
                        identity,
                        generation,
                    },
                ) if identity == expected_identity && generation == expected_generation => {
                    evaluated.state = AdmissionState::Diagnostic(DiagnosticPhase::Preflighted {
                        snapshot,
                        anchor: identity,
                        generation,
                    });
                }
                (TrustAnchorState::Anchored { .. }, AnchorObservation::Observed { .. }) => {
                    evaluated.state = refused(AdmissionRule::TrustAnchorMismatch);
                    evaluated.terminal_rule = Some(AdmissionRule::TrustAnchorMismatch);
                }
            },
            (
                AdmissionState::Diagnostic(DiagnosticPhase::Preflighted {
                    snapshot: expected_snapshot,
                    anchor: expected_anchor,
                    generation: expected_generation,
                }),
                AdmissionTransition::StabilityRecheck {
                    snapshot,
                    anchor:
                        AnchorObservation::Observed {
                            identity,
                            generation,
                        },
                    ..
                },
            ) if snapshot == expected_snapshot
                && identity == expected_anchor
                && generation == expected_generation =>
            {
                evaluated.state = AdmissionState::Admitted(AdmittedPhase::Active {
                    snapshot,
                    anchor: identity,
                    generation,
                });
            }
            (
                AdmissionState::Diagnostic(DiagnosticPhase::Preflighted { .. }),
                AdmissionTransition::StabilityRecheck { .. },
            ) => {
                evaluated.state = refused(AdmissionRule::StabilityChanged);
                evaluated.terminal_rule = Some(AdmissionRule::StabilityChanged);
            }
            (
                AdmissionState::Diagnostic(_)
                | AdmissionState::Admitted(AdmittedPhase::Active { .. }),
                AdmissionTransition::Charge { charge, .. },
            ) => {
                evaluated.consumption = charge_budget(
                    self.context.fetch,
                    self.context.budgets,
                    self.consumption,
                    charge,
                )?;
            }
            (
                AdmissionState::Diagnostic(_)
                | AdmissionState::Admitted(AdmittedPhase::Active { .. }),
                AdmissionTransition::PollCancellation {
                    work_since_last_poll,
                    ..
                },
            ) => {
                if work_since_last_poll > self.current_cx.max_unpolled_work {
                    return Err(AdmissionError::new(AdmissionRule::PollIntervalExceeded));
                }
            }
            (
                AdmissionState::Diagnostic(_),
                AdmissionTransition::Refuse {
                    reason,
                    observation: _,
                },
            ) => {
                let rule = reason.rule();
                evaluated.state = refused(rule);
                evaluated.terminal_rule = Some(rule);
            }
            (
                AdmissionState::Admitted(
                    AdmittedPhase::Active { .. }
                    | AdmittedPhase::CompletionPending { .. }
                    | AdmittedPhase::PublicationPending { .. }
                    | AdmittedPhase::PublicationReady { .. },
                ),
                AdmissionTransition::RequestCancellation {
                    cause,
                    obligations,
                    observation: _,
                },
            ) => {
                evaluated.state = AdmissionState::Cancelled {
                    cause,
                    phase: CancellationPhase::Requested,
                    remaining: obligations,
                };
            }
            (
                AdmissionState::Cancelled {
                    cause,
                    phase: CancellationPhase::Requested | CancellationPhase::Draining,
                    remaining,
                },
                AdmissionTransition::Drain {
                    completed,
                    observation: _,
                },
            ) => {
                let Some(remaining) = remaining.checked_sub(completed) else {
                    return Err(AdmissionError::new(AdmissionRule::DrainOverrun));
                };
                evaluated.state = AdmissionState::Cancelled {
                    cause,
                    phase: CancellationPhase::Draining,
                    remaining,
                };
            }
            (
                AdmissionState::Unanchored(TerminalPhase::Pending),
                AdmissionTransition::FinalizeTerminal { .. },
            ) => {
                evaluated.state = AdmissionState::Unanchored(TerminalPhase::Finalized);
            }
            (
                AdmissionState::Refused {
                    rule,
                    phase: TerminalPhase::Pending,
                },
                AdmissionTransition::FinalizeTerminal { .. },
            ) => {
                evaluated.state = AdmissionState::Refused {
                    rule,
                    phase: TerminalPhase::Finalized,
                };
            }
            (
                AdmissionState::Indeterminate {
                    rule,
                    phase: TerminalPhase::Pending,
                },
                AdmissionTransition::FinalizeTerminal { .. },
            ) => {
                evaluated.state = AdmissionState::Indeterminate {
                    rule,
                    phase: TerminalPhase::Finalized,
                };
            }
            (
                AdmissionState::Cancelled {
                    cause,
                    phase: CancellationPhase::Draining,
                    remaining,
                },
                AdmissionTransition::FinalizeTerminal { .. },
            ) => {
                if !remaining.is_empty() {
                    return Err(AdmissionError::new(AdmissionRule::DrainIncomplete));
                }
                evaluated.state = AdmissionState::Cancelled {
                    cause,
                    phase: CancellationPhase::Finalized,
                    remaining,
                };
            }
            (
                AdmissionState::Cancelled {
                    phase: CancellationPhase::Requested,
                    ..
                },
                AdmissionTransition::FinalizeTerminal { .. },
            ) => return Err(AdmissionError::new(AdmissionRule::DrainIncomplete)),
            (
                AdmissionState::Admitted(
                    AdmittedPhase::Active { .. }
                    | AdmittedPhase::CompletionPending { .. }
                    | AdmittedPhase::PublicationPending { .. }
                    | AdmittedPhase::PublicationReady { .. }
                    | AdmittedPhase::PublicationCommitting { .. },
                )
                | AdmissionState::Cancelled {
                    phase: CancellationPhase::Requested | CancellationPhase::Draining,
                    ..
                },
                AdmissionTransition::DeclareIndeterminate {
                    reason,
                    observation: _,
                },
            ) => {
                let rule = reason.rule();
                evaluated.state = AdmissionState::Indeterminate {
                    rule,
                    phase: TerminalPhase::Pending,
                };
                evaluated.terminal_rule = Some(rule);
            }
            (
                AdmissionState::Admitted(AdmittedPhase::Active {
                    snapshot,
                    anchor,
                    generation,
                }),
                AdmissionTransition::BeginReadOnlyCompletion { quiescence, .. },
            ) => {
                if !self.context.command.admits_read_only_completion()
                    || !self.context.publication.is_prohibited()
                {
                    return Err(AdmissionError::new(
                        AdmissionRule::ReadOnlyCompletionForbidden,
                    ));
                }
                evaluated.state = AdmissionState::Admitted(AdmittedPhase::CompletionPending {
                    snapshot,
                    anchor,
                    generation,
                    quiescence,
                });
            }
            (
                AdmissionState::Admitted(AdmittedPhase::CompletionPending {
                    snapshot: expected_snapshot,
                    anchor: expected_anchor,
                    generation: expected_generation,
                    ..
                }),
                AdmissionTransition::CompleteReadOnly {
                    snapshot,
                    anchor:
                        AnchorObservation::Observed {
                            identity,
                            generation,
                        },
                    receipt,
                    ..
                },
            ) if snapshot == expected_snapshot
                && identity == expected_anchor
                && generation == expected_generation =>
            {
                evaluated.state = AdmissionState::Admitted(AdmittedPhase::Completed { receipt });
            }
            (
                AdmissionState::Admitted(AdmittedPhase::CompletionPending { .. }),
                AdmissionTransition::CompleteReadOnly { .. },
            ) => {
                let rule = IndeterminateReason::PostEffectStabilityChanged.rule();
                evaluated.state = AdmissionState::Indeterminate {
                    rule,
                    phase: TerminalPhase::Pending,
                };
                evaluated.terminal_rule = Some(rule);
            }
            (
                AdmissionState::Admitted(AdmittedPhase::Active {
                    snapshot,
                    anchor,
                    generation,
                }),
                AdmissionTransition::BeginPublication { quiescence, .. },
            ) => {
                if self.context.publication.is_prohibited() {
                    return Err(AdmissionError::new(
                        AdmissionRule::PublicationAuthorityDenied,
                    ));
                }
                evaluated.state = AdmissionState::Admitted(AdmittedPhase::PublicationPending {
                    snapshot,
                    anchor,
                    generation,
                    quiescence,
                });
            }
            (
                AdmissionState::Admitted(AdmittedPhase::PublicationPending {
                    snapshot: expected_snapshot,
                    anchor: expected_anchor,
                    generation: expected_generation,
                    quiescence,
                }),
                AdmissionTransition::PublicationRecheck {
                    snapshot,
                    anchor:
                        AnchorObservation::Observed {
                            identity,
                            generation,
                        },
                    fence,
                    ..
                },
            ) if snapshot == expected_snapshot
                && identity == expected_anchor
                && generation == expected_generation =>
            {
                evaluated.state = AdmissionState::Admitted(AdmittedPhase::PublicationReady {
                    snapshot,
                    anchor: identity,
                    generation,
                    fence,
                    quiescence,
                });
            }
            (
                AdmissionState::Admitted(AdmittedPhase::PublicationPending { .. }),
                AdmissionTransition::PublicationRecheck { .. },
            ) => {
                evaluated.state = AdmissionState::Indeterminate {
                    rule: AdmissionRule::PublicationStabilityChanged,
                    phase: TerminalPhase::Pending,
                };
                evaluated.terminal_rule = Some(AdmissionRule::PublicationStabilityChanged);
            }
            (
                AdmissionState::Admitted(AdmittedPhase::PublicationReady {
                    fence, quiescence, ..
                }),
                AdmissionTransition::AuthorizePublication { receipt, .. },
            ) => {
                evaluated.state = AdmissionState::Admitted(AdmittedPhase::PublicationCommitting {
                    fence,
                    receipt,
                    quiescence,
                });
            }
            (
                AdmissionState::Admitted(AdmittedPhase::PublicationCommitting {
                    fence,
                    receipt,
                    ..
                }),
                AdmissionTransition::FinalizePublication { success },
            ) => {
                let Some(publication_capability) = publication_identity(self.context.publication)
                else {
                    return Err(AdmissionError::new(
                        AdmissionRule::PublicationAuthorityDenied,
                    ));
                };
                if success.request_identity != self.context.request_identity
                    || success.attempt != self.attempt
                    || success.publication_capability != publication_capability
                    || success.fence != fence
                    || success.receipt != receipt
                {
                    return Err(AdmissionError::new(
                        AdmissionRule::PublicationSuccessEvidenceMismatch,
                    ));
                }
                evaluated.state =
                    AdmissionState::Admitted(AdmittedPhase::Published { receipt, success });
            }
            (
                AdmissionState::Admitted(AdmittedPhase::PublicationCommitting { .. }),
                AdmissionTransition::PublicationFailed {
                    reason,
                    observation: _,
                },
            ) => {
                let rule = reason.rule();
                evaluated.state = AdmissionState::Indeterminate {
                    rule,
                    phase: TerminalPhase::Pending,
                };
                evaluated.terminal_rule = Some(rule);
            }
            (
                AdmissionState::Unanchored(TerminalPhase::Finalized)
                | AdmissionState::Refused {
                    phase: TerminalPhase::Finalized,
                    ..
                }
                | AdmissionState::Cancelled {
                    phase: CancellationPhase::Finalized,
                    ..
                },
                AdmissionTransition::Retry { next_cx, .. },
            ) => {
                if next_cx.clock != self.context.budgets.deadline.clock {
                    return Err(AdmissionError::new(AdmissionRule::ClockBindingMismatch));
                }
                let reuses_prior_binding = reuses_attempt_authority(next_cx, self.context.cx)
                    || self.history.iter().any(|transition| {
                        matches!(
                            transition,
                            AdmissionTransition::Retry {
                                next_cx: prior_cx,
                                ..
                            } if reuses_attempt_authority(next_cx, *prior_cx)
                        )
                    });
                if reuses_prior_binding {
                    return Err(AdmissionError::new(AdmissionRule::RetryCxNotFresh));
                }
                let retries = self
                    .consumption
                    .retries
                    .checked_add(1)
                    .ok_or_else(|| AdmissionError::new(AdmissionRule::RetryBudgetExceeded))?;
                if retries > self.context.budgets.retries {
                    return Err(AdmissionError::new(AdmissionRule::RetryBudgetExceeded));
                }
                evaluated.attempt = self
                    .attempt
                    .checked_add(1)
                    .ok_or_else(|| AdmissionError::new(AdmissionRule::RetryBudgetExceeded))?;
                evaluated.consumption.retries = retries;
                evaluated.cx = next_cx;
                evaluated.state = AdmissionState::Diagnostic(DiagnosticPhase::Created);
            }
            // Indeterminate work may have escaped. Reusing its old authority would
            // risk duplicate mutation, so an explicit new request is required.
            (
                AdmissionState::Indeterminate {
                    phase: TerminalPhase::Finalized,
                    ..
                },
                AdmissionTransition::Retry { .. },
            ) => {
                return Err(AdmissionError::new(
                    AdmissionRule::IndeterminateRetryForbidden,
                ));
            }
            _ => return Err(AdmissionError::new(AdmissionRule::IllegalTransition)),
        }
        Ok(evaluated)
    }

    /// Encode the complete policy and committed history into canonical bytes.
    ///
    /// The envelope contains only redacted identities, policy, budgets, and
    /// deterministic transitions. It contains no path, URL, environment,
    /// credential, process handle, or live capability.
    ///
    /// # Errors
    /// Refuses bounded output allocation failure.
    pub fn encode_canonical(&self) -> Result<Vec<u8>, AdmissionError> {
        let mut encoder = Encoder::new()?;
        encode_context(&mut encoder, &self.context);
        encoder.u16(
            u16::try_from(self.history.len())
                .map_err(|_| AdmissionError::new(AdmissionRule::EventLimitExceeded))?,
        );
        for transition in &self.history {
            encode_transition(&mut encoder, *transition);
        }
        encoder.finish()
    }

    /// Decode and replay canonical bytes into an inert audit receipt.
    ///
    /// No live machine, path, executable, network, cancellation, or publication
    /// authority is reconstructed. An entrypoint must perform fresh admission.
    ///
    /// # Errors
    /// Refuses oversized, truncated, malformed, noncanonical, unknown-version,
    /// impossible-history, or trailing-byte envelopes.
    pub fn decode_recorded(bytes: &[u8]) -> Result<RecordedAdmission, AdmissionError> {
        if bytes.len() > MAX_ADMISSION_BYTES {
            return Err(AdmissionError::new(AdmissionRule::EncodedInputTooLarge));
        }
        let mut decoder = Decoder::new(bytes);
        let context = decode_context(&mut decoder)?;
        let transition_count = usize::from(decoder.u16()?);
        if transition_count > MAX_ADMISSION_EVENTS {
            return Err(AdmissionError::new(AdmissionRule::EventLimitExceeded));
        }
        let mut transitions = Vec::new();
        transitions
            .try_reserve_exact(transition_count)
            .map_err(|_| AdmissionError::new(AdmissionRule::AllocationFailed))?;
        for _ in 0..transition_count {
            transitions.push(decode_transition(&mut decoder)?);
        }
        if !decoder.is_finished() {
            return Err(AdmissionError::new(AdmissionRule::TrailingBytes));
        }

        let mut replay = Self::try_new(context)?;
        for &transition in &transitions {
            replay
                .apply(transition)
                .map_err(|_| AdmissionError::new(AdmissionRule::ImpossibleHistory))?;
        }
        let canonical = replay.encode_canonical()?;
        if canonical != bytes {
            return Err(AdmissionError::new(AdmissionRule::NonCanonicalEncoding));
        }
        let AdmissionMachine {
            context,
            current_cx,
            state,
            consumption,
            attempt,
            last_observed_tick,
            history,
            events,
        } = replay;
        Ok(RecordedAdmission {
            context,
            state,
            consumption,
            attempt,
            current_cx,
            last_observed_tick,
            history,
            events,
            canonical,
        })
    }
}

fn reuses_attempt_authority(next: CxBinding, prior: CxBinding) -> bool {
    next.cx == prior.cx
        || next.cx == prior.cancellation
        || next.cancellation == prior.cx
        || next.cancellation == prior.cancellation
}

struct Evaluated {
    state: AdmissionState,
    consumption: BudgetConsumption,
    attempt: u16,
    cx: CxBinding,
    terminal_rule: Option<AdmissionRule>,
}

const fn minimum_terminal_events(state: AdmissionState) -> usize {
    match state {
        AdmissionState::Diagnostic(_) | AdmissionState::Admitted(AdmittedPhase::Active { .. }) => 2,
        AdmissionState::Admitted(AdmittedPhase::CompletionPending { .. })
        | AdmissionState::Admitted(AdmittedPhase::PublicationPending { .. })
        | AdmissionState::Admitted(AdmittedPhase::PublicationReady { .. })
        | AdmissionState::Admitted(AdmittedPhase::PublicationCommitting { .. }) => 2,
        AdmissionState::Unanchored(TerminalPhase::Pending)
        | AdmissionState::Refused {
            phase: TerminalPhase::Pending,
            ..
        }
        | AdmissionState::Indeterminate {
            phase: TerminalPhase::Pending,
            ..
        } => 1,
        AdmissionState::Cancelled {
            phase: CancellationPhase::Requested,
            ..
        } => 2,
        AdmissionState::Cancelled {
            phase: CancellationPhase::Draining,
            remaining,
            ..
        } => {
            if remaining.is_empty() {
                1
            } else {
                2
            }
        }
        AdmissionState::Unanchored(TerminalPhase::Finalized)
        | AdmissionState::Refused {
            phase: TerminalPhase::Finalized,
            ..
        }
        | AdmissionState::Cancelled {
            phase: CancellationPhase::Finalized,
            ..
        }
        | AdmissionState::Indeterminate {
            phase: TerminalPhase::Finalized,
            ..
        }
        | AdmissionState::Admitted(
            AdmittedPhase::Completed { .. } | AdmittedPhase::Published { .. },
        ) => 0,
    }
}

const fn refused(rule: AdmissionRule) -> AdmissionState {
    AdmissionState::Refused {
        rule,
        phase: TerminalPhase::Pending,
    }
}

fn charge_budget(
    fetch: FetchAuthority,
    budgets: AdmissionBudgets,
    current: BudgetConsumption,
    charge: BudgetCharge,
) -> Result<BudgetConsumption, AdmissionError> {
    let mut next = current;
    match charge {
        BudgetCharge::Work(value) => {
            next.work_units =
                checked_charge_u64(current.work_units, value, budgets.compute.work_units)?;
        }
        BudgetCharge::Memory(value) => {
            next.memory_bytes =
                checked_charge_u64(current.memory_bytes, value, budgets.compute.memory_bytes)?;
        }
        BudgetCharge::Processes(value) => {
            next.processes = checked_charge_u32(current.processes, value, budgets.io.processes)?;
        }
        BudgetCharge::Files(value) => {
            next.files = checked_charge_u32(current.files, value, budgets.io.files)?;
        }
        BudgetCharge::Output(value) => {
            next.output_bytes =
                checked_charge_u64(current.output_bytes, value, budgets.io.output_bytes)?;
        }
        BudgetCharge::Network { requests, bytes } => {
            if matches!(fetch, FetchAuthority::Offline) {
                return Err(AdmissionError::new(AdmissionRule::NetworkAuthorityDenied));
            }
            let next_requests =
                checked_charge_u32(current.network_requests, requests, budgets.network.requests)?;
            let next_bytes =
                checked_charge_u64(current.network_bytes, bytes, budgets.network.bytes)?;
            next.network_requests = next_requests;
            next.network_bytes = next_bytes;
        }
    }
    Ok(next)
}

fn checked_charge_u64(current: u64, added: u64, cap: u64) -> Result<u64, AdmissionError> {
    let next = current
        .checked_add(added)
        .ok_or_else(|| AdmissionError::new(AdmissionRule::BudgetOverflow))?;
    if next > cap {
        return Err(AdmissionError::new(AdmissionRule::BudgetExceeded));
    }
    Ok(next)
}

fn checked_charge_u32(current: u32, added: u32, cap: u32) -> Result<u32, AdmissionError> {
    let next = current
        .checked_add(added)
        .ok_or_else(|| AdmissionError::new(AdmissionRule::BudgetOverflow))?;
    if next > cap {
        return Err(AdmissionError::new(AdmissionRule::BudgetExceeded));
    }
    Ok(next)
}

/// Inert decoded audit record. It cannot be converted into live authority.
#[derive(Debug, PartialEq, Eq)]
pub struct RecordedAdmission {
    context: AdmissionContext,
    state: AdmissionState,
    consumption: BudgetConsumption,
    attempt: u16,
    current_cx: CxBinding,
    last_observed_tick: Option<u64>,
    history: Vec<AdmissionTransition>,
    events: Vec<AdmissionEvent>,
    canonical: Vec<u8>,
}

impl RecordedAdmission {
    /// Complete inert request policy retained by the record.
    #[must_use]
    pub const fn context(&self) -> &AdmissionContext {
        &self.context
    }

    /// Stable request identity retained for correlation only.
    #[must_use]
    pub const fn request_identity(&self) -> AuthorityId {
        self.context.request_identity
    }

    /// Recorded command class.
    #[must_use]
    pub const fn command(&self) -> CommandClass {
        self.context.command
    }

    /// Recorded terminal or in-flight state; this grants no authority.
    #[must_use]
    pub const fn state(&self) -> AdmissionState {
        self.state
    }

    /// Recorded cumulative consumption.
    #[must_use]
    pub const fn consumption(&self) -> BudgetConsumption {
        self.consumption
    }

    /// Recorded attempt ordinal.
    #[must_use]
    pub const fn attempt(&self) -> u16 {
        self.attempt
    }

    /// Final recorded attempt-scoped Cx identity; this is not a live handle.
    #[must_use]
    pub const fn current_cx(&self) -> CxBinding {
        self.current_cx
    }

    /// Last replayed monotonic tick across every recorded attempt.
    #[must_use]
    pub const fn last_observed_tick(&self) -> Option<u64> {
        self.last_observed_tick
    }

    /// Exact deterministic transition history, including evidence identities.
    #[must_use]
    pub fn history(&self) -> &[AdmissionTransition] {
        &self.history
    }

    /// Recorded deterministic events.
    #[must_use]
    pub fn events(&self) -> &[AdmissionEvent] {
        &self.events
    }

    /// Exact canonical bytes that were decoded and replay-validated.
    #[must_use]
    pub fn canonical_bytes(&self) -> &[u8] {
        &self.canonical
    }
}

/// Deterministic admission/refusal rules.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
pub enum AdmissionRule {
    /// An all-zero authority identity was supplied.
    ZeroIdentity = 0,
    /// Cancellation polling interval was zero.
    ZeroPollInterval = 1,
    /// Deadline and Cx used different clocks.
    ClockBindingMismatch = 2,
    /// Too many path capability rows were supplied.
    TooManyPathCapabilities = 3,
    /// Too many executable capability rows were supplied.
    TooManyExecutableCapabilities = 4,
    /// A path role appeared more than once.
    DuplicatePathCapability = 5,
    /// A command-required path role was absent.
    MissingPathCapability = 6,
    /// An executable role appeared more than once.
    DuplicateExecutableCapability = 7,
    /// A command-required executable role was absent.
    MissingExecutableCapability = 8,
    /// Offline authority conflicted with a positive network budget.
    OfflineNetworkBudgetConflict = 9,
    /// Fetch authority had no positive request and byte budget.
    FetchAuthorityWithoutBudget = 10,
    /// The publication-authority variant is forbidden for the command class.
    PublicationForbiddenForCommand = 11,
    /// Publication authority had no output budget.
    PublicationWithoutOutputBudget = 12,
    /// Transition is illegal from the exact current phase.
    IllegalTransition = 13,
    /// Caller-supplied monotonic tick exceeded the inclusive deadline.
    DeadlineExceeded = 14,
    /// A resource charge exceeded its cap.
    BudgetExceeded = 15,
    /// Checked resource accounting overflowed.
    BudgetOverflow = 16,
    /// Network use was attempted under offline authority.
    NetworkAuthorityDenied = 17,
    /// Publication was attempted without explicit publication authority.
    PublicationAuthorityDenied = 18,
    /// No live trust anchor was available.
    TrustAnchorUnavailable = 19,
    /// Observed trust anchor or generation did not match the request.
    TrustAnchorMismatch = 20,
    /// Source or anchor changed between preflight observations.
    StabilityChanged = 21,
    /// Source or anchor changed at the publication barrier.
    PublicationStabilityChanged = 22,
    /// Drain progress exceeded a fixed outstanding obligation.
    DrainOverrun = 23,
    /// Cancellation finalization was attempted with live obligations.
    DrainIncomplete = 24,
    /// Retry cap or attempt ordinal was exceeded.
    RetryBudgetExceeded = 25,
    /// An indeterminate attempt cannot reuse the same authority.
    IndeterminateRetryForbidden = 26,
    /// Deterministic event bound was exceeded.
    EventLimitExceeded = 27,
    /// Bounded allocation failed.
    AllocationFailed = 28,
    /// Canonical input exceeded the hard envelope limit.
    EncodedInputTooLarge = 29,
    /// Canonical envelope magic was unknown.
    UnknownMagic = 30,
    /// Canonical schema version was unknown.
    UnknownSchema = 31,
    /// Canonical input ended before a complete field.
    TruncatedEncoding = 32,
    /// Canonical input carried an unknown numeric tag.
    UnknownTag = 33,
    /// Canonical input retained trailing bytes.
    TrailingBytes = 34,
    /// Decoded transition history could not have occurred legally.
    ImpossibleHistory = 35,
    /// Input was semantically valid but not in its one canonical encoding.
    NonCanonicalEncoding = 36,
    /// The command class is intrinsically offline/read-only.
    NetworkForbiddenForCommand = 37,
    /// Publication authority identity did not match the publication path slot.
    PublicationCapabilityMismatch = 38,
    /// A path capability was not admissible for the command class.
    UnexpectedPathCapability = 39,
    /// An executable capability was not admissible for the command class.
    UnexpectedExecutableCapability = 40,
    /// Caller-supplied monotonic ticks regressed.
    ClockRegression = 41,
    /// A transition would consume the event slots needed to terminalize safely.
    TerminalHeadroomRequired = 42,
    /// Work between explicit cancellation polls exceeded the Cx bound.
    PollIntervalExceeded = 43,
    /// Retry did not bind a fresh child Cx and cancellation source.
    RetryCxNotFresh = 44,
    /// Cx, cancellation, and clock roles reused one identity.
    CxIdentityAliasing = 45,
    /// A bounded adapter input was definitively rejected before effect.
    InputRejected = 46,
    /// A source artifact was definitively rejected before effect.
    SourceRejected = 47,
    /// A required live capability could not be established.
    CapabilityRejected = 48,
    /// A required executable failed pre-effect admission.
    ExecutableRejected = 49,
    /// Declared work cannot fit the explicit budget envelope.
    BudgetInfeasible = 50,
    /// Closed command policy denied the requested pre-effect operation.
    PolicyDenied = 51,
    /// The adapter cannot prove whether an external effect occurred.
    EffectOutcomeUnknown = 52,
    /// A trustworthy drain observation could not be obtained.
    DrainObservationUnavailable = 53,
    /// Terminal-record finalization may have partially occurred.
    FinalizationOutcomeUnknown = 54,
    /// Source or anchor moved after admitted work began.
    PostEffectStabilityChanged = 55,
    /// An authorized external publication commit reported failure.
    PublicationCommitFailed = 56,
    /// The outcome of an authorized publication commit is unknown.
    PublicationOutcomeUnknown = 57,
    /// Publication visibility is known but durability is not.
    PublicationDurabilityUnknown = 58,
    /// Durable publication exists but receipt finalization is uncertain.
    PublicationFinalizationUnknown = 59,
    /// The command or publication policy cannot complete through read-only success.
    ReadOnlyCompletionForbidden = 60,
    /// Publication-success evidence did not bind the authorized transaction.
    PublicationSuccessEvidenceMismatch = 61,
    /// A publication-producing command omitted its required publication authority.
    PublicationRequiredForCommand = 62,
}

impl AdmissionRule {
    /// Stable machine-readable rule code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::ZeroIdentity => "admission-zero-identity",
            Self::ZeroPollInterval => "admission-zero-poll-interval",
            Self::ClockBindingMismatch => "admission-clock-binding-mismatch",
            Self::TooManyPathCapabilities => "admission-path-capability-limit",
            Self::TooManyExecutableCapabilities => "admission-executable-capability-limit",
            Self::DuplicatePathCapability => "admission-duplicate-path-capability",
            Self::MissingPathCapability => "admission-missing-path-capability",
            Self::DuplicateExecutableCapability => "admission-duplicate-executable-capability",
            Self::MissingExecutableCapability => "admission-missing-executable-capability",
            Self::OfflineNetworkBudgetConflict => "admission-offline-network-budget-conflict",
            Self::FetchAuthorityWithoutBudget => "admission-fetch-authority-without-budget",
            Self::PublicationForbiddenForCommand => "admission-publication-forbidden-command",
            Self::PublicationWithoutOutputBudget => "admission-publication-without-output-budget",
            Self::IllegalTransition => "admission-illegal-transition",
            Self::DeadlineExceeded => "admission-deadline-exceeded",
            Self::BudgetExceeded => "admission-budget-exceeded",
            Self::BudgetOverflow => "admission-budget-overflow",
            Self::NetworkAuthorityDenied => "admission-network-authority-denied",
            Self::PublicationAuthorityDenied => "admission-publication-authority-denied",
            Self::TrustAnchorUnavailable => "admission-trust-anchor-unavailable",
            Self::TrustAnchorMismatch => "admission-trust-anchor-mismatch",
            Self::StabilityChanged => "admission-stability-changed",
            Self::PublicationStabilityChanged => "admission-publication-stability-changed",
            Self::DrainOverrun => "admission-drain-overrun",
            Self::DrainIncomplete => "admission-drain-incomplete",
            Self::RetryBudgetExceeded => "admission-retry-budget-exceeded",
            Self::IndeterminateRetryForbidden => "admission-indeterminate-retry-forbidden",
            Self::EventLimitExceeded => "admission-event-limit",
            Self::AllocationFailed => "admission-allocation-failed",
            Self::EncodedInputTooLarge => "admission-encoded-input-limit",
            Self::UnknownMagic => "admission-unknown-magic",
            Self::UnknownSchema => "admission-unknown-schema",
            Self::TruncatedEncoding => "admission-truncated-encoding",
            Self::UnknownTag => "admission-unknown-tag",
            Self::TrailingBytes => "admission-trailing-bytes",
            Self::ImpossibleHistory => "admission-impossible-history",
            Self::NonCanonicalEncoding => "admission-noncanonical-encoding",
            Self::NetworkForbiddenForCommand => "admission-network-forbidden-command",
            Self::PublicationCapabilityMismatch => "admission-publication-capability-mismatch",
            Self::UnexpectedPathCapability => "admission-unexpected-path-capability",
            Self::UnexpectedExecutableCapability => "admission-unexpected-executable-capability",
            Self::ClockRegression => "admission-clock-regression",
            Self::TerminalHeadroomRequired => "admission-terminal-headroom-required",
            Self::PollIntervalExceeded => "admission-poll-interval-exceeded",
            Self::RetryCxNotFresh => "admission-retry-cx-not-fresh",
            Self::CxIdentityAliasing => "admission-cx-identity-aliasing",
            Self::InputRejected => "admission-input-rejected",
            Self::SourceRejected => "admission-source-rejected",
            Self::CapabilityRejected => "admission-capability-rejected",
            Self::ExecutableRejected => "admission-executable-rejected",
            Self::BudgetInfeasible => "admission-budget-infeasible",
            Self::PolicyDenied => "admission-policy-denied",
            Self::EffectOutcomeUnknown => "admission-effect-outcome-unknown",
            Self::DrainObservationUnavailable => "admission-drain-observation-unavailable",
            Self::FinalizationOutcomeUnknown => "admission-finalization-outcome-unknown",
            Self::PostEffectStabilityChanged => "admission-post-effect-stability-changed",
            Self::PublicationCommitFailed => "admission-publication-commit-failed",
            Self::PublicationOutcomeUnknown => "admission-publication-outcome-unknown",
            Self::PublicationDurabilityUnknown => "admission-publication-durability-unknown",
            Self::PublicationFinalizationUnknown => "admission-publication-finalization-unknown",
            Self::ReadOnlyCompletionForbidden => "admission-read-only-completion-forbidden",
            Self::PublicationSuccessEvidenceMismatch => {
                "admission-publication-success-evidence-mismatch"
            }
            Self::PublicationRequiredForCommand => "admission-publication-required-command",
        }
    }

    /// Ranked, closed remedy codes. No environment or secret text is retained.
    #[must_use]
    pub const fn remedies(self) -> &'static [Remedy] {
        match self {
            Self::MissingPathCapability | Self::MissingExecutableCapability => {
                &[Remedy::BindRequiredCapability, Remedy::UseDiagnosticMode]
            }
            Self::UnexpectedPathCapability | Self::UnexpectedExecutableCapability => &[
                Remedy::RemoveUnexpectedCapability,
                Remedy::UseDiagnosticMode,
            ],
            Self::OfflineNetworkBudgetConflict | Self::NetworkAuthorityDenied => {
                &[Remedy::UseOfflineCache, Remedy::BindPinnedNetworkCapability]
            }
            Self::NetworkForbiddenForCommand => {
                &[Remedy::UseOfflineCache, Remedy::UseDiagnosticMode]
            }
            Self::FetchAuthorityWithoutBudget | Self::BudgetExceeded | Self::BudgetOverflow => {
                &[Remedy::ReduceWork, Remedy::IncreaseExplicitBudget]
            }
            Self::DeadlineExceeded => &[Remedy::IncreaseExplicitDeadline, Remedy::ReduceWork],
            Self::ClockRegression => &[Remedy::RestartPreflight, Remedy::InspectDiagnostic],
            Self::PollIntervalExceeded => &[Remedy::ReduceWork, Remedy::CompleteDrain],
            Self::TerminalHeadroomRequired => &[Remedy::FinalizeCurrentAttempt, Remedy::ReduceWork],
            Self::TrustAnchorUnavailable | Self::TrustAnchorMismatch => {
                &[Remedy::BindTrustAnchor, Remedy::UseDiagnosticMode]
            }
            Self::StabilityChanged | Self::PublicationStabilityChanged => {
                &[Remedy::RestartPreflight, Remedy::QuiesceWriters]
            }
            Self::DrainOverrun | Self::DrainIncomplete => {
                &[Remedy::CompleteDrain, Remedy::MarkIndeterminate]
            }
            Self::PublicationAuthorityDenied
            | Self::PublicationWithoutOutputBudget
            | Self::PublicationCapabilityMismatch => {
                &[Remedy::BindPublicationCapability, Remedy::UseDiagnosticMode]
            }
            Self::PublicationForbiddenForCommand => {
                &[Remedy::UseDiagnosticMode, Remedy::InspectDiagnostic]
            }
            Self::PublicationRequiredForCommand => {
                &[Remedy::BindPublicationCapability, Remedy::UseDiagnosticMode]
            }
            Self::RetryCxNotFresh | Self::CxIdentityAliasing => {
                &[Remedy::BindFreshChildCx, Remedy::StartNewRequest]
            }
            Self::RetryBudgetExceeded | Self::IndeterminateRetryForbidden => {
                &[Remedy::StartNewRequest, Remedy::InspectPriorAttempt]
            }
            Self::InputRejected => &[Remedy::InspectDiagnostic, Remedy::UseSupportedSchema],
            Self::SourceRejected => &[Remedy::RestartPreflight, Remedy::InspectDiagnostic],
            Self::CapabilityRejected | Self::ExecutableRejected => {
                &[Remedy::BindRequiredCapability, Remedy::InspectDiagnostic]
            }
            Self::BudgetInfeasible => &[Remedy::ReduceWork, Remedy::IncreaseExplicitBudget],
            Self::PolicyDenied | Self::ReadOnlyCompletionForbidden => {
                &[Remedy::UseDiagnosticMode, Remedy::InspectDiagnostic]
            }
            Self::PublicationSuccessEvidenceMismatch => {
                &[Remedy::MarkIndeterminate, Remedy::InspectDiagnostic]
            }
            Self::EffectOutcomeUnknown
            | Self::FinalizationOutcomeUnknown
            | Self::PublicationCommitFailed
            | Self::PublicationOutcomeUnknown
            | Self::PublicationDurabilityUnknown
            | Self::PublicationFinalizationUnknown => {
                &[Remedy::MarkIndeterminate, Remedy::InspectDiagnostic]
            }
            Self::DrainObservationUnavailable => {
                &[Remedy::CompleteDrain, Remedy::MarkIndeterminate]
            }
            Self::PostEffectStabilityChanged => &[Remedy::RestartPreflight, Remedy::QuiesceWriters],
            Self::UnknownMagic
            | Self::UnknownSchema
            | Self::TruncatedEncoding
            | Self::UnknownTag
            | Self::TrailingBytes
            | Self::ImpossibleHistory
            | Self::NonCanonicalEncoding
            | Self::EncodedInputTooLarge => {
                &[Remedy::UseSupportedSchema, Remedy::RegenerateReceipt]
            }
            _ => &[Remedy::InspectDiagnostic, Remedy::UseDiagnosticMode],
        }
    }
}

/// Closed ranked remedy codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
pub enum Remedy {
    /// Bind the exact required capability role.
    BindRequiredCapability = 0,
    /// Use a no-authority diagnostic command.
    UseDiagnosticMode = 1,
    /// Use already materialized pinned sources.
    UseOfflineCache = 2,
    /// Bind a pinned network capability explicitly.
    BindPinnedNetworkCapability = 3,
    /// Reduce requested work.
    ReduceWork = 4,
    /// Increase an explicit numeric budget deliberately.
    IncreaseExplicitBudget = 5,
    /// Increase the explicit monotonic deadline deliberately.
    IncreaseExplicitDeadline = 6,
    /// Bind the exact expected trust anchor.
    BindTrustAnchor = 7,
    /// Restart from a complete first observation.
    RestartPreflight = 8,
    /// Quiesce writers that move admitted source state.
    QuiesceWriters = 9,
    /// Complete every fixed drain obligation.
    CompleteDrain = 10,
    /// Preserve uncertainty rather than claiming clean refusal/cancellation.
    MarkIndeterminate = 11,
    /// Bind an explicit publication target capability.
    BindPublicationCapability = 12,
    /// Create a fresh idempotency identity and admission request.
    StartNewRequest = 13,
    /// Inspect the inert prior-attempt receipt before proceeding.
    InspectPriorAttempt = 14,
    /// Use an explicitly supported schema reader.
    UseSupportedSchema = 15,
    /// Regenerate the receipt from a live admitted execution.
    RegenerateReceipt = 16,
    /// Inspect the stable structured rule and phase.
    InspectDiagnostic = 17,
    /// Remove a capability role that the command cannot consume.
    RemoveUnexpectedCapability = 18,
    /// Stop adding work and explicitly terminalize the current attempt.
    FinalizeCurrentAttempt = 19,
    /// Bind a new attempt-scoped child Cx and cancellation source.
    BindFreshChildCx = 20,
}

/// Bounded structured admission error.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdmissionError {
    rule: AdmissionRule,
}

impl AdmissionError {
    const fn new(rule: AdmissionRule) -> Self {
        Self { rule }
    }

    /// Stable refusal rule.
    #[must_use]
    pub const fn rule(self) -> AdmissionRule {
        self.rule
    }

    /// Deterministically ranked closed remedies.
    #[must_use]
    pub const fn remedies(self) -> &'static [Remedy] {
        self.rule.remedies()
    }
}

impl fmt::Display for AdmissionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.rule.code())
    }
}

impl std::error::Error for AdmissionError {}

struct Encoder {
    bytes: Vec<u8>,
}

impl Encoder {
    fn new() -> Result<Self, AdmissionError> {
        let mut bytes = Vec::new();
        bytes
            .try_reserve_exact(MAX_ADMISSION_BYTES)
            .map_err(|_| AdmissionError::new(AdmissionRule::AllocationFailed))?;
        bytes.extend_from_slice(ADMISSION_MAGIC);
        bytes.extend_from_slice(&ADMISSION_VERSION.to_le_bytes());
        bytes.extend_from_slice(&(ADMISSION_SCHEMA.len() as u16).to_le_bytes());
        bytes.extend_from_slice(ADMISSION_SCHEMA.as_bytes());
        bytes.extend_from_slice(&(ADMISSION_DOMAIN.len() as u16).to_le_bytes());
        bytes.extend_from_slice(ADMISSION_DOMAIN.as_bytes());
        Ok(Self { bytes })
    }

    fn finish(self) -> Result<Vec<u8>, AdmissionError> {
        if self.bytes.len() > MAX_ADMISSION_BYTES {
            return Err(AdmissionError::new(AdmissionRule::EncodedInputTooLarge));
        }
        Ok(self.bytes)
    }

    fn u8(&mut self, value: u8) {
        self.bytes.push(value);
    }

    fn u16(&mut self, value: u16) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn u32(&mut self, value: u32) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn u64(&mut self, value: u64) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn id(&mut self, value: AuthorityId) {
        self.bytes.extend_from_slice(value.as_bytes());
    }
}

struct Decoder<'a> {
    bytes: &'a [u8],
    cursor: usize,
}

impl<'a> Decoder<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, cursor: 0 }
    }

    fn header(&mut self) -> Result<(), AdmissionError> {
        let magic = self.take(ADMISSION_MAGIC.len())?;
        if magic != ADMISSION_MAGIC {
            return Err(AdmissionError::new(AdmissionRule::UnknownMagic));
        }
        if self.u16()? != ADMISSION_VERSION {
            return Err(AdmissionError::new(AdmissionRule::UnknownSchema));
        }
        let schema_length = usize::from(self.u16()?);
        if self.take(schema_length)? != ADMISSION_SCHEMA.as_bytes() {
            return Err(AdmissionError::new(AdmissionRule::UnknownSchema));
        }
        let domain_length = usize::from(self.u16()?);
        if self.take(domain_length)? != ADMISSION_DOMAIN.as_bytes() {
            return Err(AdmissionError::new(AdmissionRule::UnknownSchema));
        }
        Ok(())
    }

    fn is_finished(&self) -> bool {
        self.cursor == self.bytes.len()
    }

    fn take(&mut self, length: usize) -> Result<&'a [u8], AdmissionError> {
        let end = self
            .cursor
            .checked_add(length)
            .ok_or_else(|| AdmissionError::new(AdmissionRule::TruncatedEncoding))?;
        let value = self
            .bytes
            .get(self.cursor..end)
            .ok_or_else(|| AdmissionError::new(AdmissionRule::TruncatedEncoding))?;
        self.cursor = end;
        Ok(value)
    }

    fn u8(&mut self) -> Result<u8, AdmissionError> {
        Ok(self.take(1)?[0])
    }

    fn u16(&mut self) -> Result<u16, AdmissionError> {
        let mut bytes = [0; 2];
        bytes.copy_from_slice(self.take(2)?);
        Ok(u16::from_le_bytes(bytes))
    }

    fn u32(&mut self) -> Result<u32, AdmissionError> {
        let mut bytes = [0; 4];
        bytes.copy_from_slice(self.take(4)?);
        Ok(u32::from_le_bytes(bytes))
    }

    fn u64(&mut self) -> Result<u64, AdmissionError> {
        let mut bytes = [0; 8];
        bytes.copy_from_slice(self.take(8)?);
        Ok(u64::from_le_bytes(bytes))
    }

    fn id(&mut self) -> Result<AuthorityId, AdmissionError> {
        let mut bytes = [0; ID_BYTES];
        bytes.copy_from_slice(self.take(ID_BYTES)?);
        AuthorityId::try_from_bytes(bytes)
    }
}

fn encode_context(encoder: &mut Encoder, context: &AdmissionContext) {
    encoder.id(context.request_identity);
    encoder.u8(context.command as u8);
    match context.fetch {
        FetchAuthority::Offline => encoder.u8(0),
        FetchAuthority::PinnedTransport { capability } => {
            encoder.u8(1);
            encoder.id(capability);
        }
    }
    match context.publication {
        PublicationAuthority::Prohibited => encoder.u8(0),
        PublicationAuthority::BootstrapReceipt { capability } => {
            encoder.u8(1);
            encoder.id(capability);
        }
        PublicationAuthority::ConstellationLock { capability } => {
            encoder.u8(2);
            encoder.id(capability);
        }
        PublicationAuthority::ProofReceipt { capability } => {
            encoder.u8(3);
            encoder.id(capability);
        }
    }
    encoder.id(context.cx.cx);
    encoder.id(context.cx.cancellation);
    encoder.id(context.cx.clock);
    encoder.u64(context.cx.max_unpolled_work);
    encoder.id(context.budgets.deadline.clock);
    encoder.u64(context.budgets.deadline.not_after_tick);
    encoder.u64(context.budgets.compute.work_units);
    encoder.u64(context.budgets.compute.memory_bytes);
    encoder.u32(context.budgets.io.processes);
    encoder.u32(context.budgets.io.files);
    encoder.u64(context.budgets.io.output_bytes);
    encoder.u32(context.budgets.network.requests);
    encoder.u64(context.budgets.network.bytes);
    encoder.u16(context.budgets.retries);
    match context.trust_anchor {
        TrustAnchorState::Unanchored => encoder.u8(0),
        TrustAnchorState::Anchored {
            identity,
            generation,
        } => {
            encoder.u8(1);
            encoder.id(identity);
            encoder.u64(generation);
        }
    }
    encoder.u8(context.path_capabilities.len() as u8);
    for capability in &context.path_capabilities {
        encoder.u8(capability.slot as u8);
        encoder.id(capability.identity);
    }
    encoder.u8(context.executable_capabilities.len() as u8);
    for capability in &context.executable_capabilities {
        encoder.u8(capability.slot as u8);
        encoder.id(capability.identity);
    }
}

fn decode_context(decoder: &mut Decoder<'_>) -> Result<AdmissionContext, AdmissionError> {
    decoder.header()?;
    let request_identity = decoder.id()?;
    let command = CommandClass::from_tag(decoder.u8()?)
        .ok_or_else(|| AdmissionError::new(AdmissionRule::UnknownTag))?;
    let fetch = match decoder.u8()? {
        0 => FetchAuthority::Offline,
        1 => FetchAuthority::PinnedTransport {
            capability: decoder.id()?,
        },
        _ => return Err(AdmissionError::new(AdmissionRule::UnknownTag)),
    };
    let publication = match decoder.u8()? {
        0 => PublicationAuthority::Prohibited,
        1 => PublicationAuthority::BootstrapReceipt {
            capability: decoder.id()?,
        },
        2 => PublicationAuthority::ConstellationLock {
            capability: decoder.id()?,
        },
        3 => PublicationAuthority::ProofReceipt {
            capability: decoder.id()?,
        },
        _ => return Err(AdmissionError::new(AdmissionRule::UnknownTag)),
    };
    let cx = CxBinding::try_new(decoder.id()?, decoder.id()?, decoder.id()?, decoder.u64()?)?;
    let deadline = DeadlineBudget::new(decoder.id()?, decoder.u64()?);
    let budgets = AdmissionBudgets::new(
        deadline,
        ComputeBudget {
            work_units: decoder.u64()?,
            memory_bytes: decoder.u64()?,
        },
        IoBudget {
            processes: decoder.u32()?,
            files: decoder.u32()?,
            output_bytes: decoder.u64()?,
        },
        NetworkBudget {
            requests: decoder.u32()?,
            bytes: decoder.u64()?,
        },
        decoder.u16()?,
    );
    let trust_anchor = match decoder.u8()? {
        0 => TrustAnchorState::Unanchored,
        1 => TrustAnchorState::Anchored {
            identity: decoder.id()?,
            generation: decoder.u64()?,
        },
        _ => return Err(AdmissionError::new(AdmissionRule::UnknownTag)),
    };
    let path_count = usize::from(decoder.u8()?);
    if path_count > MAX_PATH_CAPABILITIES {
        return Err(AdmissionError::new(AdmissionRule::TooManyPathCapabilities));
    }
    let mut paths = Vec::new();
    paths
        .try_reserve_exact(path_count)
        .map_err(|_| AdmissionError::new(AdmissionRule::AllocationFailed))?;
    for _ in 0..path_count {
        let slot = PathSlot::from_tag(decoder.u8()?)
            .ok_or_else(|| AdmissionError::new(AdmissionRule::UnknownTag))?;
        paths.push(PathCapability::new(slot, decoder.id()?));
    }
    let executable_count = usize::from(decoder.u8()?);
    if executable_count > MAX_EXECUTABLE_CAPABILITIES {
        return Err(AdmissionError::new(
            AdmissionRule::TooManyExecutableCapabilities,
        ));
    }
    let mut executables = Vec::new();
    executables
        .try_reserve_exact(executable_count)
        .map_err(|_| AdmissionError::new(AdmissionRule::AllocationFailed))?;
    for _ in 0..executable_count {
        let slot = ExecutableSlot::from_tag(decoder.u8()?)
            .ok_or_else(|| AdmissionError::new(AdmissionRule::UnknownTag))?;
        executables.push(ExecutableCapability::new(slot, decoder.id()?));
    }
    AdmissionContext::try_new(AdmissionContextSpec {
        request_identity,
        command,
        fetch,
        publication,
        cx,
        budgets,
        trust_anchor,
        path_capabilities: &paths,
        executable_capabilities: &executables,
    })
}

fn encode_anchor(encoder: &mut Encoder, anchor: AnchorObservation) {
    match anchor {
        AnchorObservation::Unavailable => encoder.u8(0),
        AnchorObservation::Observed {
            identity,
            generation,
        } => {
            encoder.u8(1);
            encoder.id(identity);
            encoder.u64(generation);
        }
    }
}

fn decode_anchor(decoder: &mut Decoder<'_>) -> Result<AnchorObservation, AdmissionError> {
    match decoder.u8()? {
        0 => Ok(AnchorObservation::Unavailable),
        1 => Ok(AnchorObservation::Observed {
            identity: decoder.id()?,
            generation: decoder.u64()?,
        }),
        _ => Err(AdmissionError::new(AdmissionRule::UnknownTag)),
    }
}

fn encode_obligations(encoder: &mut Encoder, obligations: DrainObligations) {
    encoder.u32(obligations.processes);
    encoder.u32(obligations.files);
    encoder.u32(obligations.outputs);
}

fn encode_cx(encoder: &mut Encoder, cx: CxBinding) {
    encoder.id(cx.cx);
    encoder.id(cx.cancellation);
    encoder.id(cx.clock);
    encoder.u64(cx.max_unpolled_work);
}

fn decode_cx(decoder: &mut Decoder<'_>) -> Result<CxBinding, AdmissionError> {
    CxBinding::try_new(decoder.id()?, decoder.id()?, decoder.id()?, decoder.u64()?)
}

fn decode_obligations(decoder: &mut Decoder<'_>) -> Result<DrainObligations, AdmissionError> {
    Ok(DrainObligations {
        processes: decoder.u32()?,
        files: decoder.u32()?,
        outputs: decoder.u32()?,
    })
}

fn encode_publication_success(encoder: &mut Encoder, success: PublicationSuccessEvidence) {
    encoder.id(success.request_identity);
    encoder.u16(success.attempt);
    encoder.id(success.publication_capability);
    encoder.id(success.fence);
    encoder.id(success.receipt);
    encoder.id(success.observation);
}

fn decode_publication_success(
    decoder: &mut Decoder<'_>,
) -> Result<PublicationSuccessEvidence, AdmissionError> {
    Ok(PublicationSuccessEvidence::new(
        decoder.id()?,
        decoder.u16()?,
        decoder.id()?,
        decoder.id()?,
        decoder.id()?,
        decoder.id()?,
    ))
}

fn encode_transition(encoder: &mut Encoder, transition: AdmissionTransition) {
    encoder.u8(transition.kind() as u8);
    match transition {
        AdmissionTransition::Preflight {
            snapshot,
            anchor,
            at_tick,
        }
        | AdmissionTransition::StabilityRecheck {
            snapshot,
            anchor,
            at_tick,
        } => {
            encoder.id(snapshot);
            encode_anchor(encoder, anchor);
            encoder.u64(at_tick);
        }
        AdmissionTransition::PublicationRecheck {
            snapshot,
            anchor,
            fence,
            at_tick,
        } => {
            encoder.id(snapshot);
            encode_anchor(encoder, anchor);
            encoder.id(fence);
            encoder.u64(at_tick);
        }
        AdmissionTransition::Charge { charge, at_tick } => {
            match charge {
                BudgetCharge::Work(value) => {
                    encoder.u8(0);
                    encoder.u64(value);
                }
                BudgetCharge::Memory(value) => {
                    encoder.u8(1);
                    encoder.u64(value);
                }
                BudgetCharge::Processes(value) => {
                    encoder.u8(2);
                    encoder.u32(value);
                }
                BudgetCharge::Files(value) => {
                    encoder.u8(3);
                    encoder.u32(value);
                }
                BudgetCharge::Output(value) => {
                    encoder.u8(4);
                    encoder.u64(value);
                }
                BudgetCharge::Network { requests, bytes } => {
                    encoder.u8(5);
                    encoder.u32(requests);
                    encoder.u64(bytes);
                }
            }
            encoder.u64(at_tick);
        }
        AdmissionTransition::Refuse {
            reason,
            observation,
        } => {
            encoder.u8(reason as u8);
            encoder.id(observation);
        }
        AdmissionTransition::DeclareIndeterminate {
            reason,
            observation,
        } => {
            encoder.u8(reason as u8);
            encoder.id(observation);
        }
        AdmissionTransition::RequestCancellation {
            cause,
            obligations,
            observation,
        } => {
            encoder.u8(cause as u8);
            encode_obligations(encoder, obligations);
            encoder.id(observation);
        }
        AdmissionTransition::Drain {
            completed,
            observation,
        } => {
            encode_obligations(encoder, completed);
            encoder.id(observation);
        }
        AdmissionTransition::FinalizeTerminal { receipt } => encoder.id(receipt),
        AdmissionTransition::BeginPublication {
            quiescence,
            at_tick,
        } => {
            encoder.id(quiescence);
            encoder.u64(at_tick);
        }
        AdmissionTransition::AuthorizePublication { receipt, at_tick } => {
            encoder.id(receipt);
            encoder.u64(at_tick);
        }
        AdmissionTransition::BeginReadOnlyCompletion {
            quiescence,
            at_tick,
        } => {
            encoder.id(quiescence);
            encoder.u64(at_tick);
        }
        AdmissionTransition::CompleteReadOnly {
            snapshot,
            anchor,
            receipt,
            at_tick,
        } => {
            encoder.id(snapshot);
            encode_anchor(encoder, anchor);
            encoder.id(receipt);
            encoder.u64(at_tick);
        }
        AdmissionTransition::FinalizePublication { success } => {
            encode_publication_success(encoder, success);
        }
        AdmissionTransition::PublicationFailed {
            reason,
            observation,
        } => {
            encoder.u8(reason as u8);
            encoder.id(observation);
        }
        AdmissionTransition::Retry { next_cx, at_tick } => {
            encode_cx(encoder, next_cx);
            encoder.u64(at_tick);
        }
        AdmissionTransition::PollCancellation {
            work_since_last_poll,
            at_tick,
        } => {
            encoder.u64(work_since_last_poll);
            encoder.u64(at_tick);
        }
    }
}

fn decode_transition(decoder: &mut Decoder<'_>) -> Result<AdmissionTransition, AdmissionError> {
    let kind = TransitionKind::from_tag(decoder.u8()?)
        .ok_or_else(|| AdmissionError::new(AdmissionRule::UnknownTag))?;
    Ok(match kind {
        TransitionKind::Preflight => AdmissionTransition::Preflight {
            snapshot: decoder.id()?,
            anchor: decode_anchor(decoder)?,
            at_tick: decoder.u64()?,
        },
        TransitionKind::StabilityRecheck => AdmissionTransition::StabilityRecheck {
            snapshot: decoder.id()?,
            anchor: decode_anchor(decoder)?,
            at_tick: decoder.u64()?,
        },
        TransitionKind::Charge => {
            let charge = match decoder.u8()? {
                0 => BudgetCharge::Work(decoder.u64()?),
                1 => BudgetCharge::Memory(decoder.u64()?),
                2 => BudgetCharge::Processes(decoder.u32()?),
                3 => BudgetCharge::Files(decoder.u32()?),
                4 => BudgetCharge::Output(decoder.u64()?),
                5 => BudgetCharge::Network {
                    requests: decoder.u32()?,
                    bytes: decoder.u64()?,
                },
                _ => return Err(AdmissionError::new(AdmissionRule::UnknownTag)),
            };
            AdmissionTransition::Charge {
                charge,
                at_tick: decoder.u64()?,
            }
        }
        TransitionKind::Refuse => AdmissionTransition::Refuse {
            reason: RefusalReason::from_tag(decoder.u8()?)
                .ok_or_else(|| AdmissionError::new(AdmissionRule::UnknownTag))?,
            observation: decoder.id()?,
        },
        TransitionKind::RequestCancellation => AdmissionTransition::RequestCancellation {
            cause: CancellationCause::from_tag(decoder.u8()?)
                .ok_or_else(|| AdmissionError::new(AdmissionRule::UnknownTag))?,
            obligations: decode_obligations(decoder)?,
            observation: decoder.id()?,
        },
        TransitionKind::Drain => AdmissionTransition::Drain {
            completed: decode_obligations(decoder)?,
            observation: decoder.id()?,
        },
        TransitionKind::FinalizeTerminal => AdmissionTransition::FinalizeTerminal {
            receipt: decoder.id()?,
        },
        TransitionKind::DeclareIndeterminate => AdmissionTransition::DeclareIndeterminate {
            reason: IndeterminateReason::from_tag(decoder.u8()?)
                .ok_or_else(|| AdmissionError::new(AdmissionRule::UnknownTag))?,
            observation: decoder.id()?,
        },
        TransitionKind::BeginPublication => AdmissionTransition::BeginPublication {
            quiescence: decoder.id()?,
            at_tick: decoder.u64()?,
        },
        TransitionKind::PublicationRecheck => AdmissionTransition::PublicationRecheck {
            snapshot: decoder.id()?,
            anchor: decode_anchor(decoder)?,
            fence: decoder.id()?,
            at_tick: decoder.u64()?,
        },
        TransitionKind::AuthorizePublication => AdmissionTransition::AuthorizePublication {
            receipt: decoder.id()?,
            at_tick: decoder.u64()?,
        },
        TransitionKind::BeginReadOnlyCompletion => AdmissionTransition::BeginReadOnlyCompletion {
            quiescence: decoder.id()?,
            at_tick: decoder.u64()?,
        },
        TransitionKind::CompleteReadOnly => AdmissionTransition::CompleteReadOnly {
            snapshot: decoder.id()?,
            anchor: decode_anchor(decoder)?,
            receipt: decoder.id()?,
            at_tick: decoder.u64()?,
        },
        TransitionKind::FinalizePublication => AdmissionTransition::FinalizePublication {
            success: decode_publication_success(decoder)?,
        },
        TransitionKind::PublicationFailed => AdmissionTransition::PublicationFailed {
            reason: PublicationFailureReason::from_tag(decoder.u8()?)
                .ok_or_else(|| AdmissionError::new(AdmissionRule::UnknownTag))?,
            observation: decoder.id()?,
        },
        TransitionKind::Retry => AdmissionTransition::Retry {
            next_cx: decode_cx(decoder)?,
            at_tick: decoder.u64()?,
        },
        TransitionKind::PollCancellation => AdmissionTransition::PollCancellation {
            work_since_last_poll: decoder.u64()?,
            at_tick: decoder.u64()?,
        },
    })
}
