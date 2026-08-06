//! User-facing static admission for the Euler-disc cinematic pipeline.
//!
//! This module deliberately stops at a real, write-free plan. The durable job
//! graph, whole-bundle verifier, and quarantined mux executor have separate
//! owner Beads; non-dry invocations fail closed until those implementations
//! land. A successful inspect/plan admission proves bounded configuration and
//! asset integrity. It proves trajectory integrity only for `--trajectory`,
//! and resource admission only when a complete host-resource tuple was
//! supplied; `--run-reduced` remains explicitly unverified. No invocation
//! here proves render convergence, physical validation, authenticated host
//! capabilities, or artifact production.

use std::fmt::Write as _;
use std::fs::File;
use std::io::Read;
use std::path::{Component, Path, PathBuf};

use fs_euler_disc_e2e::{
    EULER_RENDER_TRAJECTORY_SCHEMA_VERSION, EulerRenderTrajectoryArtifact,
    RenderTrajectoryCodecBudget,
};
use fs_evidence::{
    cinematic_budget::{
        CinematicBudgetError, CinematicBudgetRepair, CinematicQualityProfile, CinematicQualityTier,
        CinematicResourceAvailability, CinematicResourceDeficit, CinematicResourceEstimate,
        CinematicResourceKind, ResourceLimitSource, admit_cinematic_budget,
    },
    cinematic_config::{CinematicConfig, CinematicMuxRequest},
    cinematic_config_codec::{
        CINEMATIC_CONFIG_DOCUMENT_SCHEMA, CinematicAssetAccessError, CinematicAssetAdmissionBudget,
        CinematicAssetDeclaration, CinematicConfigDocument, CinematicConfigDocumentError,
        MAX_CINEMATIC_CONFIG_DOCUMENT_BYTES, MAX_CINEMATIC_RESOLVED_ASSET_BYTES,
        MAX_CINEMATIC_RESOLVED_ASSET_TOTAL_BYTES,
    },
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

use crate::{CommandOutput, exit};

/// Stable result schema for cinematic CLI records.
pub const CINEMATIC_CLI_RESULT_SCHEMA: &str = "frankensim.cinematic.cli-result.v1";
/// Stable diagnostic schema for cinematic CLI records.
pub const CINEMATIC_CLI_DIAGNOSTIC_SCHEMA: &str = "frankensim.cinematic.cli-diagnostic.v1";
/// Full reconstructible config document schema consumed by this CLI.
pub const CINEMATIC_CLI_CONFIG_SCHEMA: &str = CINEMATIC_CONFIG_DOCUMENT_SCHEMA;
/// Maximum config bytes accepted by filesystem-backed CLI admission.
pub const MAX_CINEMATIC_CONFIG_BYTES: u64 = MAX_CINEMATIC_CONFIG_DOCUMENT_BYTES as u64;
/// Maximum bytes read for one external cinematic asset during static admission.
pub const MAX_CINEMATIC_ASSET_BYTES: u64 = MAX_CINEMATIC_RESOLVED_ASSET_BYTES as u64;
/// Maximum bytes read across all external assets during static admission.
pub const MAX_CINEMATIC_TOTAL_ASSET_BYTES: u64 = MAX_CINEMATIC_RESOLVED_ASSET_TOTAL_BYTES as u64;
/// Maximum encoded trajectory bytes admitted by the CLI's inspection seam.
pub const MAX_CINEMATIC_TRAJECTORY_BYTES: u64 = 1024 * 1024 * 1024;

const READ_TILE_BYTES: usize = 64 * 1024;
const CINEMATIC_USAGE: &str = "frankensim [--json] cinematic <inspect|storyboard|daily|representative-4k-frame|final|resume|verify|mux> <config.fscine> <trajectory-source> [--dry-run] [--memory-bytes <n> --free-storage-bytes <n> --wall-time-s <n> --workers <n> --paths-per-second <n>]; trajectory-source is --trajectory <artifact> for verify/mux, and either --trajectory <artifact> or --run-reduced for every other mode";
const STATIC_AUTHORITY: &str = "bounded-static-composition-admission";
const STATIC_NO_CLAIM: &str = "does not execute or verify a render, produce audio/media artifacts, authenticate declared capabilities, prove render convergence, or validate Euler-disc physics/acoustics";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum CinematicMode {
    Unknown,
    Inspect,
    Storyboard,
    Daily,
    Representative4kFrame,
    Final,
    Resume,
    Verify,
    Mux,
}

impl CinematicMode {
    fn parse(value: &str) -> Option<Self> {
        match value {
            "inspect" => Some(Self::Inspect),
            "storyboard" => Some(Self::Storyboard),
            "daily" => Some(Self::Daily),
            "representative-4k-frame" => Some(Self::Representative4kFrame),
            "final" => Some(Self::Final),
            "resume" => Some(Self::Resume),
            "verify" => Some(Self::Verify),
            "mux" => Some(Self::Mux),
            _ => None,
        }
    }

    const fn name(self) -> &'static str {
        match self {
            Self::Unknown => "unknown",
            Self::Inspect => "inspect",
            Self::Storyboard => "storyboard",
            Self::Daily => "daily",
            Self::Representative4kFrame => "representative-4k-frame",
            Self::Final => "final",
            Self::Resume => "resume",
            Self::Verify => "verify",
            Self::Mux => "mux",
        }
    }

    const fn required_tier(self) -> Option<CinematicQualityTier> {
        match self {
            Self::Storyboard => Some(CinematicQualityTier::StoryboardSmoke),
            Self::Daily => Some(CinematicQualityTier::Daily1080p),
            Self::Representative4kFrame => Some(CinematicQualityTier::Qualification4kFrame),
            Self::Final => Some(CinematicQualityTier::Final4k),
            Self::Unknown | Self::Inspect | Self::Resume | Self::Verify | Self::Mux => None,
        }
    }

    const fn needs_resources(self) -> bool {
        matches!(
            self,
            Self::Storyboard
                | Self::Daily
                | Self::Representative4kFrame
                | Self::Final
                | Self::Resume
        )
    }

    const fn downstream_owner(self) -> Option<&'static str> {
        match self {
            Self::Unknown | Self::Inspect => None,
            Self::Storyboard | Self::Daily => Some("frankensim-h7xu5.8.3"),
            Self::Representative4kFrame | Self::Final | Self::Resume => {
                Some("frankensim-h7xu5.8.2")
            }
            Self::Verify => Some("frankensim-h7xu5.8.4"),
            Self::Mux => Some("frankensim-h7xu5.8.5"),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum TrajectorySource {
    Existing(PathBuf),
    ReducedCampaign,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CinematicRequest {
    mode: CinematicMode,
    config_path: PathBuf,
    trajectory_source: TrajectorySource,
    dry_run: bool,
    availability: Option<CinematicResourceAvailability>,
}

#[derive(Default)]
struct AvailabilityBuilder {
    memory_bytes: Option<u64>,
    free_storage_bytes: Option<u64>,
    wall_time_available_s: Option<u64>,
    worker_capacity: Option<u16>,
    measured_camera_paths_per_second: Option<u64>,
}

impl AvailabilityBuilder {
    fn any(&self) -> bool {
        self.memory_bytes.is_some()
            || self.free_storage_bytes.is_some()
            || self.wall_time_available_s.is_some()
            || self.worker_capacity.is_some()
            || self.measured_camera_paths_per_second.is_some()
    }

    fn finish(
        self,
        required: bool,
        mode: CinematicMode,
    ) -> Result<Option<CinematicResourceAvailability>, Failure> {
        if !required && !self.any() {
            return Ok(None);
        }
        let missing = [
            ("--memory-bytes", self.memory_bytes.is_none()),
            ("--free-storage-bytes", self.free_storage_bytes.is_none()),
            ("--wall-time-s", self.wall_time_available_s.is_none()),
            ("--workers", self.worker_capacity.is_none()),
            (
                "--paths-per-second",
                self.measured_camera_paths_per_second.is_none(),
            ),
        ]
        .into_iter()
        .find_map(|(flag, missing)| missing.then_some(flag));
        if let Some(flag) = missing {
            return Err(Failure::one(
                exit::USAGE,
                Diagnostic::new(
                    mode,
                    "cinematic-resource-flag-missing",
                    "arguments.resources",
                    format!("resource admission is missing `{flag}`"),
                    vec![format!(
                        "supply `{flag}` together with all five resource facts"
                    )],
                ),
            ));
        }
        Ok(Some(CinematicResourceAvailability {
            memory_bytes: self.memory_bytes.unwrap_or(0),
            free_storage_bytes: self.free_storage_bytes.unwrap_or(0),
            wall_time_available_s: self.wall_time_available_s.unwrap_or(0),
            worker_capacity: self.worker_capacity.unwrap_or(0),
            measured_camera_paths_per_second: self.measured_camera_paths_per_second.unwrap_or(0),
        }))
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Diagnostic {
    mode: CinematicMode,
    code: &'static str,
    field_path: String,
    message: String,
    unit: Option<&'static str>,
    required: Option<u64>,
    available: Option<u64>,
    repairs: Vec<String>,
}

impl Diagnostic {
    fn new(
        mode: CinematicMode,
        code: &'static str,
        field_path: impl Into<String>,
        message: impl Into<String>,
        repairs: Vec<String>,
    ) -> Self {
        Self {
            mode,
            code,
            field_path: field_path.into(),
            message: message.into(),
            unit: None,
            required: None,
            available: None,
            repairs,
        }
    }

    fn with_resource(mut self, unit: &'static str, required: u64, available: u64) -> Self {
        self.unit = Some(unit);
        self.required = Some(required);
        self.available = Some(available);
        self
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Failure {
    exit_code: u8,
    diagnostics: Vec<Diagnostic>,
}

impl Failure {
    fn one(exit_code: u8, diagnostic: Diagnostic) -> Self {
        Self {
            exit_code,
            diagnostics: vec![diagnostic],
        }
    }
}

#[derive(Clone, Debug)]
struct TrajectoryFacts {
    source: &'static str,
    verified: bool,
    byte_len: Option<u64>,
    sample_count: Option<u32>,
    transition_count: Option<u32>,
    chunk_count: Option<u32>,
}

#[derive(Clone, Debug)]
struct StaticPlan {
    mode: CinematicMode,
    profile: CinematicQualityProfile,
    config: CinematicConfig,
    namespace: String,
    trajectory: TrajectoryFacts,
    estimate: Option<CinematicResourceEstimate>,
    dry_run: bool,
}

pub(super) fn run(arguments: Vec<String>, json: bool) -> CommandOutput {
    if matches!(arguments.as_slice(), [flag] if flag == "help" || flag == "--help") {
        return help(json);
    }
    let mode = arguments
        .first()
        .and_then(|value| CinematicMode::parse(value))
        .unwrap_or(CinematicMode::Unknown);
    let request = match parse_request(&arguments) {
        Ok(request) => request,
        Err(failure) => return format_failure(mode, json, failure),
    };
    let gate = CancelGate::new_clock_free();
    run_cinematic_with_gate_request(request, json, &gate)
}

/// Execute cinematic parsing and static admission with a caller-owned
/// cancellation gate. This is the bounded G4 seam used before the future job
/// graph owns process-level signal propagation.
#[must_use]
pub fn run_cinematic_with_gate(
    arguments: impl IntoIterator<Item = String>,
    json: bool,
    gate: &CancelGate,
) -> CommandOutput {
    let arguments = arguments.into_iter().collect::<Vec<_>>();
    if matches!(arguments.as_slice(), [flag] if flag == "help" || flag == "--help") {
        return help(json);
    }
    let mode = arguments
        .first()
        .and_then(|value| CinematicMode::parse(value))
        .unwrap_or(CinematicMode::Unknown);
    match parse_request(&arguments) {
        Ok(request) => run_cinematic_with_gate_request(request, json, gate),
        Err(failure) => format_failure(mode, json, failure),
    }
}

fn run_cinematic_with_gate_request(
    request: CinematicRequest,
    json: bool,
    gate: &CancelGate,
) -> CommandOutput {
    if gate.is_requested() {
        return format_failure(request.mode, json, cancelled(request.mode));
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x63_69_6e_65_6d_61_74_69,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        build_static_plan(&request, &cx)
    });
    finish_static_build(request.mode, json, gate, result)
}

fn finish_static_build(
    mode: CinematicMode,
    json: bool,
    gate: &CancelGate,
    result: Result<StaticPlan, Failure>,
) -> CommandOutput {
    if gate.is_requested() {
        return format_failure(mode, json, cancelled(mode));
    }
    match result {
        Ok(plan) => finish_plan(json, plan),
        Err(failure) => format_failure(mode, json, failure),
    }
}

fn parse_request(arguments: &[String]) -> Result<CinematicRequest, Failure> {
    let Some(mode_name) = arguments.first() else {
        return Err(usage_failure(
            CinematicMode::Unknown,
            "a cinematic mode is required",
        ));
    };
    let mode = CinematicMode::parse(mode_name)
        .ok_or_else(|| usage_failure(CinematicMode::Unknown, "the cinematic mode is unknown"))?;
    let Some(config) = arguments.get(1) else {
        return Err(usage_failure(mode, "a cinematic config path is required"));
    };
    if config.is_empty() || config.starts_with('-') {
        return Err(usage_failure(
            mode,
            "the cinematic config operand is invalid",
        ));
    }

    let mut trajectory = None;
    let mut dry_run = false;
    let mut availability = AvailabilityBuilder::default();
    let mut index = 2usize;
    while index < arguments.len() {
        match arguments[index].as_str() {
            "--dry-run" if !dry_run => {
                dry_run = true;
                index += 1;
            }
            "--run-reduced" if trajectory.is_none() => {
                trajectory = Some(TrajectorySource::ReducedCampaign);
                index += 1;
            }
            "--trajectory" if trajectory.is_none() => {
                let value = flag_value(arguments, index, mode, "--trajectory")?;
                trajectory = Some(TrajectorySource::Existing(PathBuf::from(value)));
                index += 2;
            }
            "--memory-bytes" if availability.memory_bytes.is_none() => {
                availability.memory_bytes =
                    Some(parse_u64_flag(arguments, index, mode, "--memory-bytes")?);
                index += 2;
            }
            "--free-storage-bytes" if availability.free_storage_bytes.is_none() => {
                availability.free_storage_bytes = Some(parse_u64_flag(
                    arguments,
                    index,
                    mode,
                    "--free-storage-bytes",
                )?);
                index += 2;
            }
            "--wall-time-s" if availability.wall_time_available_s.is_none() => {
                availability.wall_time_available_s =
                    Some(parse_u64_flag(arguments, index, mode, "--wall-time-s")?);
                index += 2;
            }
            "--workers" if availability.worker_capacity.is_none() => {
                availability.worker_capacity =
                    Some(parse_u16_flag(arguments, index, mode, "--workers")?);
                index += 2;
            }
            "--paths-per-second" if availability.measured_camera_paths_per_second.is_none() => {
                availability.measured_camera_paths_per_second = Some(parse_u64_flag(
                    arguments,
                    index,
                    mode,
                    "--paths-per-second",
                )?);
                index += 2;
            }
            _ => {
                return Err(usage_failure(
                    mode,
                    "a flag is unknown, duplicated, missing its value, or conflicts with the trajectory selector",
                ));
            }
        }
    }
    let trajectory_source = trajectory.ok_or_else(|| {
        let message = if matches!(mode, CinematicMode::Verify | CinematicMode::Mux) {
            "verify and mux require `--trajectory <artifact>`"
        } else {
            "exactly one of `--trajectory <artifact>` or `--run-reduced` is required"
        };
        usage_failure(mode, message)
    })?;
    if matches!(mode, CinematicMode::Verify | CinematicMode::Mux)
        && matches!(&trajectory_source, TrajectorySource::ReducedCampaign)
    {
        return Err(usage_failure(
            mode,
            "verify and mux consume an existing trajectory artifact; `--run-reduced` is not a valid source for these modes",
        ));
    }
    let availability = availability.finish(mode.needs_resources(), mode)?;
    Ok(CinematicRequest {
        mode,
        config_path: PathBuf::from(config),
        trajectory_source,
        dry_run,
        availability,
    })
}

fn flag_value<'a>(
    arguments: &'a [String],
    index: usize,
    mode: CinematicMode,
    flag: &'static str,
) -> Result<&'a str, Failure> {
    let value = arguments
        .get(index + 1)
        .filter(|value| !value.is_empty() && !value.starts_with('-'))
        .ok_or_else(|| usage_failure(mode, &format!("`{flag}` requires one operand")))?;
    Ok(value)
}

fn parse_u64_flag(
    arguments: &[String],
    index: usize,
    mode: CinematicMode,
    flag: &'static str,
) -> Result<u64, Failure> {
    let value = flag_value(arguments, index, mode, flag)?;
    value
        .parse::<u64>()
        .map_err(|_| usage_failure(mode, &format!("`{flag}` requires a non-negative integer")))
}

fn parse_u16_flag(
    arguments: &[String],
    index: usize,
    mode: CinematicMode,
    flag: &'static str,
) -> Result<u16, Failure> {
    let value = flag_value(arguments, index, mode, flag)?;
    value
        .parse::<u16>()
        .map_err(|_| usage_failure(mode, &format!("`{flag}` requires a non-negative integer")))
}

fn usage_failure(mode: CinematicMode, message: &str) -> Failure {
    Failure::one(
        exit::USAGE,
        Diagnostic::new(
            mode,
            "cinematic-cli-usage",
            "arguments",
            message,
            vec![CINEMATIC_USAGE.to_owned()],
        ),
    )
}

fn build_static_plan(request: &CinematicRequest, cx: &Cx<'_>) -> Result<StaticPlan, Failure> {
    checkpoint(request.mode, cx)?;
    let config_bytes = read_bounded_file(
        &request.config_path,
        MAX_CINEMATIC_CONFIG_BYTES,
        request.mode,
        "config",
        cx,
    )?;
    let document = CinematicConfigDocument::from_bytes(&config_bytes)
        .map_err(|error| document_failure(request.mode, error))?;
    if !is_safe_relative(Path::new(document.artifact_locator_hint())) {
        return Err(Failure::one(
            exit::REFUSED,
            Diagnostic::new(
                request.mode,
                "cinematic-artifact-root-not-relative",
                "config.artifact_root",
                "the artifact root must be a normalized repository-relative path",
                vec!["replace absolute, parent, or dot components with a normalized relative output root".to_owned()],
            ),
        ));
    }
    let profile = document
        .quality_profile()
        .map_err(|error| document_failure(request.mode, error))?;
    validate_trajectory_reference(document.trajectory(), request.mode)?;
    if request.mode == CinematicMode::Mux
        && matches!(document.mux_request(), CinematicMuxRequest::None)
    {
        return Err(Failure::one(
            exit::REFUSED,
            Diagnostic::new(
                request.mode,
                "cinematic-mux-not-requested",
                "config.mux",
                "mux mode requires an explicit quarantined-adapter request",
                vec![
                    "version the configuration with a supported mux request and declared quarantined-mux capability"
                        .to_owned(),
                ],
            ),
        ));
    }
    if let Some(required) = request.mode.required_tier()
        && profile.input().tier != required
    {
        return Err(Failure::one(
            exit::REFUSED,
            Diagnostic::new(
                request.mode,
                "cinematic-profile-mode-conflict",
                "config.quality_profile",
                format!(
                    "mode `{}` requires `{}`, but the document selects `{}`",
                    request.mode.name(),
                    tier_name(required),
                    tier_name(profile.input().tier),
                ),
                vec![format!(
                    "create or select a versioned `{}` configuration; quality is never changed implicitly",
                    tier_name(required)
                )],
            ),
        ));
    }
    let estimate = request
        .availability
        .map(|available| admit_cinematic_budget(&profile, available))
        .transpose()
        .map_err(|error| budget_failure(request.mode, error))?
        .map(|admitted| admitted.estimate());

    let base = request
        .config_path
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let base = base.canonicalize().map_err(|_| {
        Failure::one(
            exit::INPUT,
            Diagnostic::new(
                request.mode,
                "cinematic-config-base-unavailable",
                "config",
                "the configuration directory could not be resolved",
                vec!["provide a readable configuration in an existing directory".to_owned()],
            ),
        )
    })?;
    let mut aggregate_asset_bytes = 0u64;
    let config = document
        .admit_with_asset_resolver(
            CinematicAssetAdmissionBudget::DEFAULT,
            || {
                cx.checkpoint()
                    .map_err(|_| CinematicAssetAccessError::Cancelled)
            },
            |_, _, declaration| resolve_asset(&base, declaration, &mut aggregate_asset_bytes, cx),
        )
        .map_err(|error| document_failure(request.mode, error))?;
    checkpoint(request.mode, cx)?;

    let trajectory = match &request.trajectory_source {
        TrajectorySource::Existing(path) => {
            verify_trajectory(path, document.trajectory(), request.mode, cx)?
        }
        TrajectorySource::ReducedCampaign => TrajectoryFacts {
            source: "requested-reduced-campaign",
            verified: false,
            byte_len: None,
            sample_count: None,
            transition_count: None,
            chunk_count: None,
        },
    };
    Ok(StaticPlan {
        mode: request.mode,
        profile,
        config,
        namespace: document.artifact_namespace().to_owned(),
        trajectory,
        estimate,
        dry_run: request.dry_run,
    })
}

fn resolve_asset(
    base: &Path,
    declaration: &CinematicAssetDeclaration,
    aggregate: &mut u64,
    cx: &Cx<'_>,
) -> Result<Vec<u8>, CinematicAssetAccessError> {
    if cx.checkpoint().is_err() {
        return Err(CinematicAssetAccessError::Cancelled);
    }
    let relative = Path::new(declaration.locator_hint());
    if !is_safe_relative(relative) {
        return Err(CinematicAssetAccessError::Unavailable);
    }
    let path = base
        .join(relative)
        .canonicalize()
        .map_err(|_| CinematicAssetAccessError::Unavailable)?;
    if !path.starts_with(base) {
        return Err(CinematicAssetAccessError::Unavailable);
    }
    let metadata = path
        .metadata()
        .map_err(|_| CinematicAssetAccessError::Unavailable)?;
    if !metadata.is_file() {
        return Err(CinematicAssetAccessError::Unavailable);
    }
    let bytes = metadata.len();
    if bytes > MAX_CINEMATIC_ASSET_BYTES {
        return Err(CinematicAssetAccessError::TooLarge);
    }
    let next = aggregate
        .checked_add(bytes)
        .ok_or(CinematicAssetAccessError::TooLarge)?;
    if next > MAX_CINEMATIC_TOTAL_ASSET_BYTES {
        return Err(CinematicAssetAccessError::TooLarge);
    }
    let result =
        read_file_bytes(&path, MAX_CINEMATIC_ASSET_BYTES, cx).map_err(|kind| match kind {
            ReadFailure::Cancelled => CinematicAssetAccessError::Cancelled,
            ReadFailure::Unavailable => CinematicAssetAccessError::Unavailable,
            ReadFailure::TooLarge => CinematicAssetAccessError::TooLarge,
            ReadFailure::Capacity => CinematicAssetAccessError::Capacity,
        })?;
    *aggregate = aggregate
        .checked_add(result.len() as u64)
        .ok_or(CinematicAssetAccessError::TooLarge)?;
    if *aggregate > MAX_CINEMATIC_TOTAL_ASSET_BYTES {
        return Err(CinematicAssetAccessError::TooLarge);
    }
    Ok(result)
}

fn is_safe_relative(path: &Path) -> bool {
    !path.as_os_str().is_empty()
        && !path.is_absolute()
        && path
            .components()
            .all(|component| matches!(component, Component::Normal(_)))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ReadFailure {
    Cancelled,
    Unavailable,
    TooLarge,
    Capacity,
}

fn read_file_bytes(path: &Path, maximum: u64, cx: &Cx<'_>) -> Result<Vec<u8>, ReadFailure> {
    if cx.checkpoint().is_err() {
        return Err(ReadFailure::Cancelled);
    }
    let metadata = path.metadata().map_err(|_| ReadFailure::Unavailable)?;
    if !metadata.is_file() {
        return Err(ReadFailure::Unavailable);
    }
    let declared = metadata.len();
    if declared > maximum {
        return Err(ReadFailure::TooLarge);
    }
    let mut file = File::open(path).map_err(|_| ReadFailure::Unavailable)?;
    if !file
        .metadata()
        .map_err(|_| ReadFailure::Unavailable)?
        .is_file()
    {
        return Err(ReadFailure::Unavailable);
    }
    let capacity = usize::try_from(declared).map_err(|_| ReadFailure::TooLarge)?;
    let mut output = Vec::new();
    output
        .try_reserve_exact(capacity)
        .map_err(|_| ReadFailure::Capacity)?;
    let mut tile = [0u8; READ_TILE_BYTES];
    loop {
        if cx.checkpoint().is_err() {
            return Err(ReadFailure::Cancelled);
        }
        let count = file.read(&mut tile).map_err(|_| ReadFailure::Unavailable)?;
        if count == 0 {
            break;
        }
        let next = (output.len() as u64)
            .checked_add(count as u64)
            .ok_or(ReadFailure::TooLarge)?;
        if next > maximum {
            return Err(ReadFailure::TooLarge);
        }
        output
            .try_reserve(count)
            .map_err(|_| ReadFailure::Capacity)?;
        output.extend_from_slice(&tile[..count]);
    }
    Ok(output)
}

fn read_bounded_file(
    path: &Path,
    maximum: u64,
    mode: CinematicMode,
    field: &'static str,
    cx: &Cx<'_>,
) -> Result<Vec<u8>, Failure> {
    read_file_bytes(path, maximum, cx).map_err(|kind| {
        let (code, message, repair) = match kind {
            ReadFailure::Cancelled => return cancelled(mode),
            ReadFailure::Unavailable => (
                "cinematic-input-unavailable",
                "the input could not be opened or read",
                "provide a readable regular file",
            ),
            ReadFailure::TooLarge => (
                "cinematic-input-too-large",
                "the input exceeds its bounded byte ceiling",
                "provide a smaller admitted input",
            ),
            ReadFailure::Capacity => (
                "cinematic-input-capacity",
                "bounded input storage could not be reserved",
                "reduce input size or make more memory available",
            ),
        };
        Failure::one(
            exit::INPUT,
            Diagnostic::new(mode, code, field, message, vec![repair.to_owned()]),
        )
    })
}

fn verify_trajectory(
    path: &Path,
    expected: fs_evidence::cinematic_config::CinematicComponentRef,
    mode: CinematicMode,
    cx: &Cx<'_>,
) -> Result<TrajectoryFacts, Failure> {
    let metadata = path.metadata().map_err(|_| {
        Failure::one(
            exit::INPUT,
            Diagnostic::new(
                mode,
                "cinematic-trajectory-unavailable",
                "trajectory",
                "the trajectory artifact could not be opened",
                vec!["provide a readable canonical Euler trajectory artifact".to_owned()],
            ),
        )
    })?;
    if !metadata.is_file() {
        return Err(Failure::one(
            exit::INPUT,
            Diagnostic::new(
                mode,
                "cinematic-trajectory-unavailable",
                "trajectory",
                "the trajectory input is not a regular file",
                vec!["provide a readable canonical Euler trajectory artifact".to_owned()],
            ),
        ));
    }
    let byte_len = metadata.len();
    if byte_len > MAX_CINEMATIC_TRAJECTORY_BYTES {
        return Err(Failure::one(
            exit::INPUT,
            Diagnostic::new(
                mode,
                "cinematic-trajectory-too-large",
                "trajectory",
                "the trajectory artifact exceeds the CLI inspection ceiling",
                vec![
                    "provide a bounded trajectory artifact or raise the versioned CLI envelope"
                        .to_owned(),
                ],
            )
            .with_resource("bytes", byte_len, MAX_CINEMATIC_TRAJECTORY_BYTES),
        ));
    }
    let mut file = File::open(path).map_err(|_| {
        Failure::one(
            exit::INPUT,
            Diagnostic::new(
                mode,
                "cinematic-trajectory-unavailable",
                "trajectory",
                "the trajectory artifact could not be opened",
                vec!["provide a readable canonical Euler trajectory artifact".to_owned()],
            ),
        )
    })?;
    if !file
        .metadata()
        .map_err(|_| {
            Failure::one(
                exit::INPUT,
                Diagnostic::new(
                    mode,
                    "cinematic-trajectory-unavailable",
                    "trajectory",
                    "the trajectory artifact metadata could not be read",
                    vec!["provide a readable canonical Euler trajectory artifact".to_owned()],
                ),
            )
        })?
        .is_file()
    {
        return Err(Failure::one(
            exit::INPUT,
            Diagnostic::new(
                mode,
                "cinematic-trajectory-unavailable",
                "trajectory",
                "the opened trajectory input is not a regular file",
                vec!["provide a readable canonical Euler trajectory artifact".to_owned()],
            ),
        ));
    }
    let budget = RenderTrajectoryCodecBudget {
        max_artifact_bytes: MAX_CINEMATIC_TRAJECTORY_BYTES,
        ..RenderTrajectoryCodecBudget::DEFAULT
    };
    let artifact =
        EulerRenderTrajectoryArtifact::read_from(&mut file, budget, cx).map_err(|_| {
            if cx.is_cancel_requested() {
                cancelled(mode)
            } else {
                Failure::one(
                    exit::INPUT,
                    Diagnostic::new(
                        mode,
                        "cinematic-trajectory-refused",
                        "trajectory",
                        "the trajectory codec refused the artifact",
                        vec![
                            "reproduce or repair the canonical trajectory artifact and retry"
                                .to_owned(),
                        ],
                    ),
                )
            }
        })?;
    let receipt = artifact.receipt();
    if receipt.artifact_identity() != expected.identity() {
        return Err(Failure::one(
            exit::REFUSED,
            Diagnostic::new(
                mode,
                "cinematic-trajectory-identity-mismatch",
                "config.trajectory",
                "the verified trajectory bytes do not match the configured identity",
                vec!["select the exact configured artifact or version the configuration with the new identity".to_owned()],
            ),
        ));
    }
    Ok(TrajectoryFacts {
        source: "verified-artifact",
        verified: true,
        byte_len: Some(receipt.byte_len()),
        sample_count: Some(receipt.sample_count()),
        transition_count: Some(receipt.transition_count()),
        chunk_count: Some(receipt.chunk_count()),
    })
}

fn validate_trajectory_reference(
    reference: fs_evidence::cinematic_config::CinematicComponentRef,
    mode: CinematicMode,
) -> Result<(), Failure> {
    if reference.version() == u32::from(EULER_RENDER_TRAJECTORY_SCHEMA_VERSION) {
        return Ok(());
    }
    Err(Failure::one(
        exit::REFUSED,
        Diagnostic::new(
            mode,
            "cinematic-trajectory-version-mismatch",
            "config.trajectory",
            "the trajectory reference version is not the supported trajectory schema version",
            vec![format!(
                "set the trajectory reference version to {EULER_RENDER_TRAJECTORY_SCHEMA_VERSION} after producing that exact schema"
            )],
        ),
    ))
}

fn checkpoint(mode: CinematicMode, cx: &Cx<'_>) -> Result<(), Failure> {
    if cx.checkpoint().is_err() {
        Err(cancelled(mode))
    } else {
        Ok(())
    }
}

fn cancelled(mode: CinematicMode) -> Failure {
    Failure::one(
        exit::CANCELLED,
        Diagnostic::new(
            mode,
            "cinematic-cancelled",
            "execution.cancel_gate",
            "cinematic admission observed cancellation before publication",
            vec!["retry the same immutable request when cancellation is cleared".to_owned()],
        ),
    )
}

fn document_failure(mode: CinematicMode, error: CinematicConfigDocumentError) -> Failure {
    if matches!(
        error,
        CinematicConfigDocumentError::AssetAccess {
            kind: CinematicAssetAccessError::Cancelled,
            ..
        }
    ) {
        return cancelled(mode);
    }
    let exit_code = match error {
        CinematicConfigDocumentError::DocumentTooLarge { .. }
        | CinematicConfigDocumentError::InvalidUtf8
        | CinematicConfigDocumentError::AssetAccess { .. }
        | CinematicConfigDocumentError::Capacity => exit::INPUT,
        _ => exit::REFUSED,
    };
    let mut diagnostic = Diagnostic::new(
        mode,
        error.code(),
        error.field_path(),
        "the versioned cinematic configuration was not admitted",
        vec![document_repair(&error).to_owned()],
    );
    if let CinematicConfigDocumentError::DocumentTooLarge { bytes, maximum } = error {
        diagnostic = diagnostic.with_resource("bytes", bytes as u64, maximum as u64);
    }
    Failure::one(exit_code, diagnostic)
}

fn document_repair(error: &CinematicConfigDocumentError) -> &'static str {
    match error {
        CinematicConfigDocumentError::UnsupportedSchema => {
            "use the exact supported cinematic config document schema"
        }
        CinematicConfigDocumentError::BudgetProfileMismatch { .. } => {
            "bind both budget references to the selected named profile identity and identity version"
        }
        CinematicConfigDocumentError::AssetAccess { .. } => {
            "provide a readable, repository-relative asset within the admitted byte envelope"
        }
        CinematicConfigDocumentError::AssetIdentityMismatch { .. } => {
            "restore the expected asset bytes or version the configuration with the changed asset"
        }
        CinematicConfigDocumentError::MissingField { .. } => {
            "add the required field without introducing hidden defaults"
        }
        _ => "repair the named logical field and retry strict admission",
    }
}

fn budget_failure(mode: CinematicMode, error: CinematicBudgetError) -> Failure {
    match error {
        CinematicBudgetError::InsufficientResources {
            deficits, repairs, ..
        } => Failure {
            exit_code: exit::REFUSED,
            diagnostics: deficits
                .iter()
                .map(|deficit| budget_deficit_diagnostic(mode, *deficit, &repairs))
                .collect(),
        },
        CinematicBudgetError::MissingThroughputMeasurement => Failure::one(
            exit::REFUSED,
            Diagnostic::new(
                mode,
                "cinematic-budget-missing-throughput",
                "resources.paths_per_second",
                "static wall-time admission requires a positive measured camera-path throughput",
                vec!["measure this renderer/profile on the target host and supply --paths-per-second".to_owned()],
            )
            .with_resource("camera_paths_per_second", 1, 0),
        ),
        other => Failure::one(
            exit::REFUSED,
            Diagnostic::new(
                mode,
                other.code(),
                "config.quality_profile",
                "the named cinematic quality profile was structurally refused",
                vec!["repair and version the quality profile; do not silently substitute another tier".to_owned()],
            ),
        ),
    }
}

fn budget_deficit_diagnostic(
    mode: CinematicMode,
    deficit: CinematicResourceDeficit,
    repairs: &[CinematicBudgetRepair],
) -> Diagnostic {
    let (suffix, field, unit) = match (deficit.kind, deficit.source) {
        (CinematicResourceKind::LiveMemoryBytes, ResourceLimitSource::ProfileEnvelope) => (
            "profile-memory",
            "config.quality_profile.memory_ceiling_bytes",
            "bytes",
        ),
        (CinematicResourceKind::LiveMemoryBytes, ResourceLimitSource::HostAvailability) => {
            ("host-memory", "resources.memory_bytes", "bytes")
        }
        (CinematicResourceKind::StorageBytes, ResourceLimitSource::ProfileEnvelope) => (
            "profile-storage",
            "config.quality_profile.output_ceiling_bytes",
            "bytes",
        ),
        (CinematicResourceKind::StorageBytes, ResourceLimitSource::HostAvailability) => (
            "host-storage",
            "resources.usable_storage_bytes_after_reserve",
            "bytes",
        ),
        (CinematicResourceKind::PerFrameWallTimeSeconds, _) => (
            "per-frame-time",
            "config.quality_profile.per_frame_wall_time_ceiling_s",
            "seconds",
        ),
        (CinematicResourceKind::SequenceWallTimeSeconds, ResourceLimitSource::ProfileEnvelope) => (
            "profile-sequence-time",
            "config.quality_profile.sequence_wall_time_ceiling_s",
            "seconds",
        ),
        (CinematicResourceKind::SequenceWallTimeSeconds, ResourceLimitSource::HostAvailability) => {
            ("host-sequence-time", "resources.wall_time_s", "seconds")
        }
        (CinematicResourceKind::Workers, _) => ("host-workers", "resources.workers", "workers"),
    };
    Diagnostic::new(
        mode,
        match suffix {
            "profile-memory" => "cinematic-budget-profile-memory",
            "host-memory" => "cinematic-budget-host-memory",
            "profile-storage" => "cinematic-budget-profile-storage",
            "host-storage" => "cinematic-budget-host-storage",
            "per-frame-time" => "cinematic-budget-per-frame-time",
            "profile-sequence-time" => "cinematic-budget-profile-sequence-time",
            "host-sequence-time" => "cinematic-budget-host-sequence-time",
            _ => "cinematic-budget-host-workers",
        },
        field,
        "the conservative cinematic estimate exceeds an explicit limit",
        repairs
            .iter()
            .map(|repair| repair_name(*repair).to_owned())
            .collect(),
    )
    .with_resource(unit, deficit.required, deficit.available)
}

const fn repair_name(repair: CinematicBudgetRepair) -> &'static str {
    match repair {
        CinematicBudgetRepair::IncreaseHostMemory => "increase-host-memory",
        CinematicBudgetRepair::IncreaseFreeStorage => "increase-free-storage",
        CinematicBudgetRepair::ExtendWallTime => "extend-wall-time",
        CinematicBudgetRepair::IncreaseWorkerCapacity => "increase-worker-capacity",
        CinematicBudgetRepair::LowerPreviewSppWithNewConfiguration => {
            "lower-preview-spp-with-new-configuration"
        }
        CinematicBudgetRepair::ReducePreviewAovsWithNewConfiguration => {
            "reduce-preview-aovs-with-new-configuration"
        }
        CinematicBudgetRepair::ShortenRangeWithNewConfiguration => {
            "shorten-range-with-new-configuration"
        }
        CinematicBudgetRepair::RaiseProfileEnvelopeWithNewConfiguration => {
            "raise-profile-envelope-with-new-configuration"
        }
    }
}

fn finish_plan(json: bool, plan: StaticPlan) -> CommandOutput {
    if plan.mode == CinematicMode::Inspect || plan.dry_run {
        return CommandOutput {
            exit_code: exit::SUCCESS,
            stdout: format_plan(
                json,
                &plan,
                if plan.mode == CinematicMode::Inspect {
                    "inspected"
                } else {
                    "planned"
                },
                None,
            ),
            stderr: String::new(),
        };
    }
    let Some(dependency) = plan.mode.downstream_owner() else {
        return format_failure(
            plan.mode,
            json,
            Failure::one(
                exit::UNAVAILABLE,
                Diagnostic::new(
                    plan.mode,
                    "cinematic-stage-owner-missing",
                    "execution.stage",
                    "the admitted mode has no authoritative execution-stage owner",
                    vec![
                        "bind the mode to a concrete producer Bead before enabling execution"
                            .to_owned(),
                    ],
                ),
            ),
        );
    };
    let diagnostic = Diagnostic::new(
        plan.mode,
        "cinematic-stage-unavailable",
        "execution.stage",
        format!(
            "static admission passed, but `{}` cannot execute until `{dependency}` supplies its authoritative stage",
            plan.mode.name(),
        ),
        vec![format!(
            "complete and verify `{dependency}`; do not substitute placeholder media or a mock success"
        )],
    );
    CommandOutput {
        exit_code: exit::UNAVAILABLE,
        stdout: format_plan(json, &plan, "unavailable", Some(dependency)),
        stderr: format_diagnostic(json, &diagnostic),
    }
}

#[allow(clippy::too_many_lines)]
fn format_plan(json: bool, plan: &StaticPlan, status: &str, dependency: Option<&str>) -> String {
    let profile = plan.profile.input();
    if !json {
        let mut out = format!(
            "status={status}\ncommand=cinematic\nmode={}\nquality_profile={}\nresolution={}x{}\nfirst_frame={}\nframe_count={}\nspp_floor={}\nspp_ceiling={}\nprofile_identity={}\ncomposition_identity={}\nconfigured_trajectory_artifact_identity={}\ntrajectory_partition_identity={}\nimage_identity={}\naudio_identity={}\nmux_identity={}\nartifact_namespace={}\ntrajectory_source={}\ntrajectory_verified={}\nresource_admission={}\ncapability_authority=caller-declared-unverified\nwould_write=false\nauthority={STATIC_AUTHORITY}\nno_claim={STATIC_NO_CLAIM}\n",
            plan.mode.name(),
            tier_name(profile.tier),
            profile.width_pixels,
            profile.height_pixels,
            profile.first_frame,
            profile.frame_count,
            profile.spp_floor,
            profile.spp_ceiling,
            plan.profile.identity().to_hex(),
            plan.config.composition_identity().to_hex(),
            plan.config.input().trajectory.identity().to_hex(),
            plan.config.trajectory_identity().to_hex(),
            plan.config.image_identity().to_hex(),
            plan.config.audio_identity().to_hex(),
            plan.config.mux_identity().to_hex(),
            escape_text(&plan.namespace),
            plan.trajectory.source,
            plan.trajectory.verified,
            if plan.estimate.is_some() {
                "admitted"
            } else {
                "not-requested"
            },
        );
        if let Some(estimate) = plan.estimate {
            let _ = write!(
                out,
                "estimated_camera_paths={}\nestimated_live_memory_bytes={}\nestimated_total_storage_bytes={}\nestimated_per_frame_wall_time_s={}\nestimated_sequence_wall_time_s={}\n",
                estimate.camera_paths,
                estimate.live_memory_bytes,
                estimate.total_storage_bytes,
                estimate.per_frame_wall_time_s,
                estimate.sequence_wall_time_s,
            );
        }
        if let Some(dependency) = dependency {
            let _ = writeln!(out, "dependency={dependency}");
        }
        return out;
    }

    let mut out = String::from("{\"schema\":");
    push_json_string(&mut out, CINEMATIC_CLI_RESULT_SCHEMA);
    out.push_str(",\"command\":\"cinematic\",\"mode\":");
    push_json_string(&mut out, plan.mode.name());
    out.push_str(",\"status\":");
    push_json_string(&mut out, status);
    out.push_str(",\"quality_profile\":");
    push_json_string(&mut out, tier_name(profile.tier));
    let _ = write!(
        out,
        ",\"width_pixels\":{},\"height_pixels\":{},\"first_frame\":{},\"frame_count\":{},\"spp_floor\":{},\"spp_ceiling\":{}",
        profile.width_pixels,
        profile.height_pixels,
        profile.first_frame,
        profile.frame_count,
        profile.spp_floor,
        profile.spp_ceiling,
    );
    for (field, identity) in [
        ("profile_identity", plan.profile.identity()),
        ("composition_identity", plan.config.composition_identity()),
        (
            "configured_trajectory_artifact_identity",
            plan.config.input().trajectory.identity(),
        ),
        (
            "trajectory_partition_identity",
            plan.config.trajectory_identity(),
        ),
        ("image_identity", plan.config.image_identity()),
        ("audio_identity", plan.config.audio_identity()),
        ("mux_identity", plan.config.mux_identity()),
    ] {
        out.push(',');
        push_json_string(&mut out, field);
        out.push(':');
        push_json_string(&mut out, &identity.to_hex());
    }
    out.push_str(",\"artifact_namespace\":");
    push_json_string(&mut out, &plan.namespace);
    out.push_str(",\"trajectory\":{\"source\":");
    push_json_string(&mut out, plan.trajectory.source);
    let _ = write!(
        out,
        ",\"verified\":{},\"byte_len\":{},\"sample_count\":{},\"transition_count\":{},\"chunk_count\":{}}}",
        plan.trajectory.verified,
        optional_u64(plan.trajectory.byte_len),
        optional_u32(plan.trajectory.sample_count),
        optional_u32(plan.trajectory.transition_count),
        optional_u32(plan.trajectory.chunk_count),
    );
    out.push_str(",\"resource_estimate\":");
    if let Some(estimate) = plan.estimate {
        let _ = write!(
            out,
            "{{\"camera_paths\":{},\"live_memory_bytes\":{},\"total_storage_bytes\":{},\"per_frame_wall_time_s\":{},\"sequence_wall_time_s\":{}}}",
            estimate.camera_paths,
            estimate.live_memory_bytes,
            estimate.total_storage_bytes,
            estimate.per_frame_wall_time_s,
            estimate.sequence_wall_time_s,
        );
    } else {
        out.push_str("null");
    }
    out.push_str(",\"planned_stages\":[");
    for (index, stage) in planned_stages(plan.mode).iter().enumerate() {
        if index != 0 {
            out.push(',');
        }
        push_json_string(&mut out, stage);
    }
    out.push_str("],\"resource_admission\":");
    push_json_string(
        &mut out,
        if plan.estimate.is_some() {
            "admitted"
        } else {
            "not-requested"
        },
    );
    out.push_str(",\"capability_authority\":\"caller-declared-unverified\",\"would_write\":false");
    if let Some(dependency) = dependency {
        out.push_str(",\"dependency\":");
        push_json_string(&mut out, dependency);
    }
    out.push_str(",\"authority\":");
    push_json_string(&mut out, STATIC_AUTHORITY);
    out.push_str(",\"no_claim\":");
    push_json_string(&mut out, STATIC_NO_CLAIM);
    out.push_str("}\n");
    out
}

fn planned_stages(mode: CinematicMode) -> &'static [&'static str] {
    match mode {
        CinematicMode::Unknown => &[],
        CinematicMode::Inspect => &["config", "assets", "trajectory"],
        CinematicMode::Storyboard
        | CinematicMode::Daily
        | CinematicMode::Representative4kFrame
        | CinematicMode::Final
        | CinematicMode::Resume => &[
            "trajectory",
            "raw-frames",
            "image-finishing",
            "audio",
            "sequence-verification",
        ],
        CinematicMode::Verify => &["sequence-verification"],
        CinematicMode::Mux => &["mux-derivative"],
    }
}

fn format_failure(mode: CinematicMode, json: bool, failure: Failure) -> CommandOutput {
    let mut stderr = String::new();
    for diagnostic in &failure.diagnostics {
        stderr.push_str(&format_diagnostic(json, diagnostic));
    }
    let status = if failure.exit_code == exit::CANCELLED {
        "cancelled"
    } else {
        "refused"
    };
    let stdout = if json {
        format!(
            "{{\"schema\":\"{CINEMATIC_CLI_RESULT_SCHEMA}\",\"command\":\"cinematic\",\"mode\":\"{}\",\"status\":\"{status}\",\"finding_count\":{},\"would_write\":false}}\n",
            mode.name(),
            failure.diagnostics.len(),
        )
    } else {
        format!(
            "status={status}\ncommand=cinematic\nmode={}\nfinding_count={}\nwould_write=false\n",
            mode.name(),
            failure.diagnostics.len(),
        )
    };
    CommandOutput {
        exit_code: failure.exit_code,
        stdout,
        stderr,
    }
}

fn format_diagnostic(json: bool, diagnostic: &Diagnostic) -> String {
    if !json {
        let mut out = format!(
            "ERROR {} [{}]: {}\n",
            diagnostic.code,
            escape_text(&diagnostic.field_path),
            escape_text(&diagnostic.message),
        );
        if let (Some(unit), Some(required), Some(available)) =
            (diagnostic.unit, diagnostic.required, diagnostic.available)
        {
            let _ = writeln!(
                out,
                "RESOURCE unit={unit} required={required} available={available}"
            );
        }
        for (index, repair) in diagnostic.repairs.iter().enumerate() {
            let _ = writeln!(out, "FIX {}: {}", index + 1, escape_text(repair));
        }
        return out;
    }
    let mut out = String::from("{\"schema\":");
    push_json_string(&mut out, CINEMATIC_CLI_DIAGNOSTIC_SCHEMA);
    out.push_str(",\"command\":\"cinematic\",\"mode\":");
    push_json_string(&mut out, diagnostic.mode.name());
    out.push_str(",\"severity\":\"error\",\"code\":");
    push_json_string(&mut out, diagnostic.code);
    out.push_str(",\"field_path\":");
    push_json_string(&mut out, &diagnostic.field_path);
    out.push_str(",\"message\":");
    push_json_string(&mut out, &diagnostic.message);
    if let (Some(unit), Some(required), Some(available)) =
        (diagnostic.unit, diagnostic.required, diagnostic.available)
    {
        out.push_str(",\"unit\":");
        push_json_string(&mut out, unit);
        let _ = write!(out, ",\"required\":{required},\"available\":{available}");
    }
    out.push_str(",\"ranked_fixes\":[");
    for (index, repair) in diagnostic.repairs.iter().enumerate() {
        if index != 0 {
            out.push(',');
        }
        push_json_string(&mut out, repair);
    }
    out.push_str("]}\n");
    out
}

fn help(json: bool) -> CommandOutput {
    let stdout = if json {
        let mut out = format!(
            "{{\"schema\":\"{CINEMATIC_CLI_RESULT_SCHEMA}\",\"command\":\"cinematic-help\",\"status\":\"ok\",\"usage\":"
        );
        push_json_string(&mut out, CINEMATIC_USAGE);
        out.push_str(",\"config_schema\":");
        push_json_string(&mut out, CINEMATIC_CLI_CONFIG_SCHEMA);
        out.push_str("}\n");
        out
    } else {
        format!(
            "{CINEMATIC_USAGE}\nconfig schema: {CINEMATIC_CLI_CONFIG_SCHEMA}\nexample: frankensim --json cinematic final configs/euler.fscine --trajectory artifacts/euler.trajectory --dry-run --memory-bytes 17179869184 --free-storage-bytes 549755813888 --wall-time-s 5184000 --workers 64 --paths-per-second 10000000\n"
        )
    };
    CommandOutput {
        exit_code: exit::SUCCESS,
        stdout,
        stderr: String::new(),
    }
}

const fn tier_name(tier: CinematicQualityTier) -> &'static str {
    match tier {
        CinematicQualityTier::StoryboardSmoke => "storyboard-smoke",
        CinematicQualityTier::Daily1080p => "daily-1080p",
        CinematicQualityTier::Qualification4kFrame => "qualification-4k-frame",
        CinematicQualityTier::Final4k => "final-4k",
    }
}

fn optional_u64(value: Option<u64>) -> String {
    value.map_or_else(|| "null".to_owned(), |value| value.to_string())
}

fn optional_u32(value: Option<u32>) -> String {
    value.map_or_else(|| "null".to_owned(), |value| value.to_string())
}

fn push_json_string(output: &mut String, value: &str) {
    output.push('"');
    for character in value.chars() {
        match character {
            '"' => output.push_str("\\\""),
            '\\' => output.push_str("\\\\"),
            '\n' => output.push_str("\\n"),
            '\r' => output.push_str("\\r"),
            '\t' => output.push_str("\\t"),
            character if character <= '\u{1f}' => {
                let _ = write!(output, "\\u{:04x}", u32::from(character));
            }
            character => output.push(character),
        }
    }
    output.push('"');
}

fn escape_text(value: &str) -> String {
    let mut output = String::new();
    for character in value.chars() {
        match character {
            '\\' => output.push_str("\\\\"),
            '\n' => output.push_str("\\n"),
            '\r' => output.push_str("\\r"),
            '\t' => output.push_str("\\t"),
            character if character.is_control() => {
                let _ = write!(output, "\\u{{{:x}}}", u32::from(character));
            }
            character => output.push(character),
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn post_build_cancellation_has_priority_over_an_ordinary_refusal() {
        let gate = CancelGate::new_clock_free();
        gate.request();
        let ordinary = Failure::one(
            exit::INPUT,
            Diagnostic::new(
                CinematicMode::Inspect,
                "ordinary-input-refusal",
                "config",
                "ordinary refusal",
                vec!["ordinary repair".to_owned()],
            ),
        );
        let output = finish_static_build(CinematicMode::Inspect, true, &gate, Err(ordinary));
        assert_eq!(output.exit_code, exit::CANCELLED);
        assert!(output.stdout.contains("\"status\":\"cancelled\""));
        assert!(output.stderr.contains("cinematic-cancelled"));
        assert!(!output.stderr.contains("ordinary-input-refusal"));
    }
}
