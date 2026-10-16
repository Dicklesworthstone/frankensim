//! Deterministic, film-specific orchestration for the Euler cinematic.
//!
//! This module owns the dependency graph and its transactional state machine;
//! it does not implement a second scheduler or manufacture missing renderer,
//! temporal-finishing, audio, or mux algorithms. Callers supply those stage
//! implementations through [`CinematicJobBackend`]. A node becomes reusable
//! only after the backend's owner-specific checker succeeds and publication
//! returns an exact [`CinematicPublishedArtifact`].
//!
//! Node identities exclude orchestration resource ceilings: changing a
//! scheduling allowance does not invalidate scientifically identical bytes.
//! The complete plan identity includes those ceilings so admission and replay
//! still retain the exact resource contract.

use core::{fmt, panic::AssertUnwindSafe};
use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::{ContentHash, DomainHasher, hash_domain};
use fs_evidence::cinematic_config::{CinematicConfig, CinematicMuxRequest};
use fs_exec::{AdmittedBudget, BudgetRefusal, Cx};

use crate::render_checkpoint::euler_render_checkpoint_frame_identity;
use crate::render_scene_bridge::EulerCinematicScene;
use crate::render_sharding::{EulerRenderFrameInput, EulerUniformRenderPlan};

/// Semantic version of the immutable cinematic DAG.
pub const CINEMATIC_JOB_SCHEMA_VERSION: u16 = 1;
/// Domain of one node identity.
pub const CINEMATIC_JOB_NODE_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.cinematic-job-node.v1";
/// Domain of one node's expected output identity.
pub const CINEMATIC_JOB_OUTPUT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.cinematic-job-output.v1";
/// Domain of the complete plan, including its resource contract.
pub const CINEMATIC_JOB_PLAN_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.cinematic-job-plan.v1";
/// Domain of a bounded resume snapshot.
pub const CINEMATIC_JOB_SNAPSHOT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.cinematic-job-snapshot.v1";
/// Domain of the canonical camera-shot/frame partition.
pub const CINEMATIC_SHOT_PLAN_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.cinematic-shot-plan.v1";

const ZERO_HASH: ContentHash = ContentHash([0; 32]);
const MAX_EVENTS_PER_NODE: u64 = 5;

/// One domain stage in the fixed Euler film graph.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CinematicJobKind {
    /// Verify or produce the canonical trajectory artifact.
    Trajectory,
    /// Render one immutable tile/sample rectangle.
    RenderShard { shard_ordinal: u64 },
    /// Merge all shards of one event-delimited raw segment.
    MergeRawSegment {
        frame_ordinal: u64,
        segment_index: u64,
    },
    /// Finish one segment using its complete temporal neighborhood.
    FinishSegment {
        frame_ordinal: u64,
        segment_index: u64,
    },
    /// Seal the checked image-sequence inventory.
    ImageSequence,
    /// Derive the sample-aligned control stream.
    AudioControls,
    /// Map physical/artistic controls to excitation.
    AudioExcitation,
    /// Resample excitation onto the audio clock.
    AudioResampling,
    /// Run the strictly ordered modal synthesizer.
    ModalSynthesis,
    /// Mix, meter, and encode the WAV master.
    AudioMaster,
    /// Independently verify the image/audio bundle.
    BundleVerification,
    /// Produce an explicitly non-authoritative delivery derivative.
    MuxDerivative,
}

impl CinematicJobKind {
    fn tag(self) -> u8 {
        match self {
            Self::Trajectory => 1,
            Self::RenderShard { .. } => 2,
            Self::MergeRawSegment { .. } => 3,
            Self::FinishSegment { .. } => 4,
            Self::ImageSequence => 5,
            Self::AudioControls => 6,
            Self::AudioExcitation => 7,
            Self::AudioResampling => 8,
            Self::ModalSynthesis => 9,
            Self::AudioMaster => 10,
            Self::BundleVerification => 11,
            Self::MuxDerivative => 12,
        }
    }

    fn encode(self, output: &mut Vec<u8>) {
        output.push(self.tag());
        match self {
            Self::RenderShard { shard_ordinal } => {
                output.extend_from_slice(&shard_ordinal.to_le_bytes());
            }
            Self::MergeRawSegment {
                frame_ordinal,
                segment_index,
            }
            | Self::FinishSegment {
                frame_ordinal,
                segment_index,
            } => {
                output.extend_from_slice(&frame_ordinal.to_le_bytes());
                output.extend_from_slice(&segment_index.to_le_bytes());
            }
            _ => {}
        }
    }
}

/// Artifact family expected from one node.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum CinematicArtifactKind {
    Trajectory = 1,
    RawRenderShard,
    RawSegment,
    FinishedSegment,
    ImageSequence,
    AudioControls,
    AudioExcitation,
    ResampledAudio,
    ModalAudio,
    WavMaster,
    VerifiedBundle,
    MuxDerivative,
}

impl CinematicArtifactKind {
    fn for_job(kind: CinematicJobKind) -> Self {
        match kind {
            CinematicJobKind::Trajectory => Self::Trajectory,
            CinematicJobKind::RenderShard { .. } => Self::RawRenderShard,
            CinematicJobKind::MergeRawSegment { .. } => Self::RawSegment,
            CinematicJobKind::FinishSegment { .. } => Self::FinishedSegment,
            CinematicJobKind::ImageSequence => Self::ImageSequence,
            CinematicJobKind::AudioControls => Self::AudioControls,
            CinematicJobKind::AudioExcitation => Self::AudioExcitation,
            CinematicJobKind::AudioResampling => Self::ResampledAudio,
            CinematicJobKind::ModalSynthesis => Self::ModalAudio,
            CinematicJobKind::AudioMaster => Self::WavMaster,
            CinematicJobKind::BundleVerification => Self::VerifiedBundle,
            CinematicJobKind::MuxDerivative => Self::MuxDerivative,
        }
    }

    fn from_tag(tag: u8) -> Option<Self> {
        Some(match tag {
            1 => Self::Trajectory,
            2 => Self::RawRenderShard,
            3 => Self::RawSegment,
            4 => Self::FinishedSegment,
            5 => Self::ImageSequence,
            6 => Self::AudioControls,
            7 => Self::AudioExcitation,
            8 => Self::ResampledAudio,
            9 => Self::ModalAudio,
            10 => Self::WavMaster,
            11 => Self::VerifiedBundle,
            12 => Self::MuxDerivative,
            _ => return None,
        })
    }
}

/// Nonzero per-node execution and publication ceiling.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CinematicNodeBudget {
    work_units: u64,
    max_output_bytes: u64,
}

/// Exact domain work represented by a node. These counters describe logical
/// render cells and paths, not elapsed time or a completion-time promise.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct CinematicNodeWork {
    render_frame_ordinal: Option<u64>,
    render_segment_index: Option<u64>,
    render_tile_start: u64,
    render_tiles: u64,
    samples_per_tile: u64,
    render_paths: u64,
    shot_ordinal: Option<u64>,
}

impl CinematicNodeWork {
    const NONE: Self = Self {
        render_frame_ordinal: None,
        render_segment_index: None,
        render_tile_start: 0,
        render_tiles: 0,
        samples_per_tile: 0,
        render_paths: 0,
        shot_ordinal: None,
    };

    fn render_shard(
        render_frame_ordinal: u64,
        render_segment_index: u64,
        render_tile_start: u64,
        render_tiles: u64,
        samples_per_tile: u64,
        render_paths: u64,
    ) -> Result<Self, CinematicJobPlanError> {
        if render_tiles == 0 || samples_per_tile == 0 || render_paths == 0 {
            return Err(CinematicJobPlanError::Incompatible(
                "empty render-shard work",
            ));
        }
        let _ = checked_add(render_tile_start, render_tiles, "render tile range")?;
        let _ = checked_mul(render_tiles, samples_per_tile, "render tile-sample work")?;
        Ok(Self {
            render_frame_ordinal: Some(render_frame_ordinal),
            render_segment_index: Some(render_segment_index),
            render_tile_start,
            render_tiles,
            samples_per_tile,
            render_paths,
            shot_ordinal: None,
        })
    }

    const fn finished_segment(shot_ordinal: u64) -> Self {
        Self {
            render_frame_ordinal: None,
            render_segment_index: None,
            render_tile_start: 0,
            render_tiles: 0,
            samples_per_tile: 0,
            render_paths: 0,
            shot_ordinal: Some(shot_ordinal),
        }
    }

    /// Parent frame ordinal for render-shard work.
    #[must_use]
    pub const fn render_frame_ordinal(self) -> Option<u64> {
        self.render_frame_ordinal
    }

    /// Event-delimited segment index for render-shard work.
    #[must_use]
    pub const fn render_segment_index(self) -> Option<u64> {
        self.render_segment_index
    }

    /// Inclusive row-major start of this shard's tile range.
    #[must_use]
    pub const fn render_tile_start(self) -> u64 {
        self.render_tile_start
    }

    /// Number of logical render tiles in this node.
    #[must_use]
    pub const fn render_tiles(self) -> u64 {
        self.render_tiles
    }

    /// Number of samples represented for each logical tile.
    #[must_use]
    pub const fn samples_per_tile(self) -> u64 {
        self.samples_per_tile
    }

    /// Exact path count retained by the render sharding plan.
    #[must_use]
    pub const fn render_paths(self) -> u64 {
        self.render_paths
    }

    /// Canonical shot ordinal for a temporal finishing node.
    #[must_use]
    pub const fn shot_ordinal(self) -> Option<u64> {
        self.shot_ordinal
    }

    fn tile_samples(self) -> u64 {
        self.render_tiles
            .checked_mul(self.samples_per_tile)
            .expect("node work is checked during planning")
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct DerivedShotRange {
    identity: ContentHash,
    shot_id: u64,
    first_frame_position: u64,
    frame_count: u64,
}

impl CinematicNodeBudget {
    /// Construct a dimensioned node allowance.
    pub fn try_new(work_units: u64, max_output_bytes: u64) -> Result<Self, CinematicJobPlanError> {
        if work_units == 0 {
            return Err(CinematicJobPlanError::InvalidLimit("node work units"));
        }
        if max_output_bytes == 0 {
            return Err(CinematicJobPlanError::InvalidLimit("node output bytes"));
        }
        Ok(Self {
            work_units,
            max_output_bytes,
        })
    }

    #[must_use]
    pub const fn work_units(self) -> u64 {
        self.work_units
    }

    #[must_use]
    pub const fn max_output_bytes(self) -> u64 {
        self.max_output_bytes
    }
}

/// Explicit allowances for every stage family.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CinematicJobBudgets {
    pub trajectory: CinematicNodeBudget,
    pub render_shard: CinematicNodeBudget,
    pub raw_merge: CinematicNodeBudget,
    pub temporal_finish: CinematicNodeBudget,
    pub image_sequence: CinematicNodeBudget,
    pub audio_controls: CinematicNodeBudget,
    pub audio_excitation: CinematicNodeBudget,
    pub audio_resampling: CinematicNodeBudget,
    pub modal_synthesis: CinematicNodeBudget,
    pub audio_master: CinematicNodeBudget,
    pub bundle_verification: CinematicNodeBudget,
    pub mux_derivative: CinematicNodeBudget,
}

impl CinematicJobBudgets {
    /// Apply one already-validated allowance to every family. Useful for
    /// bounded smoke jobs; production callers normally declare each field.
    #[must_use]
    pub const fn uniform(value: CinematicNodeBudget) -> Self {
        Self {
            trajectory: value,
            render_shard: value,
            raw_merge: value,
            temporal_finish: value,
            image_sequence: value,
            audio_controls: value,
            audio_excitation: value,
            audio_resampling: value,
            modal_synthesis: value,
            audio_master: value,
            bundle_verification: value,
            mux_derivative: value,
        }
    }

    fn for_kind(self, kind: CinematicJobKind) -> CinematicNodeBudget {
        match kind {
            CinematicJobKind::Trajectory => self.trajectory,
            CinematicJobKind::RenderShard { .. } => self.render_shard,
            CinematicJobKind::MergeRawSegment { .. } => self.raw_merge,
            CinematicJobKind::FinishSegment { .. } => self.temporal_finish,
            CinematicJobKind::ImageSequence => self.image_sequence,
            CinematicJobKind::AudioControls => self.audio_controls,
            CinematicJobKind::AudioExcitation => self.audio_excitation,
            CinematicJobKind::AudioResampling => self.audio_resampling,
            CinematicJobKind::ModalSynthesis => self.modal_synthesis,
            CinematicJobKind::AudioMaster => self.audio_master,
            CinematicJobKind::BundleVerification => self.bundle_verification,
            CinematicJobKind::MuxDerivative => self.mux_derivative,
        }
    }
}

/// Content identities of the concrete stage implementations/checkers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CinematicStageIdentities {
    pub trajectory: ContentHash,
    pub render_shard: ContentHash,
    pub raw_merge: ContentHash,
    pub temporal_finish: ContentHash,
    pub image_sequence: ContentHash,
    pub audio_controls: ContentHash,
    pub audio_excitation: ContentHash,
    pub audio_resampling: ContentHash,
    pub modal_synthesis: ContentHash,
    pub audio_master: ContentHash,
    pub bundle_verifier: ContentHash,
    pub mux_adapter: ContentHash,
}

impl CinematicStageIdentities {
    fn for_kind(self, kind: CinematicJobKind) -> ContentHash {
        match kind {
            CinematicJobKind::Trajectory => self.trajectory,
            CinematicJobKind::RenderShard { .. } => self.render_shard,
            CinematicJobKind::MergeRawSegment { .. } => self.raw_merge,
            CinematicJobKind::FinishSegment { .. } => self.temporal_finish,
            CinematicJobKind::ImageSequence => self.image_sequence,
            CinematicJobKind::AudioControls => self.audio_controls,
            CinematicJobKind::AudioExcitation => self.audio_excitation,
            CinematicJobKind::AudioResampling => self.audio_resampling,
            CinematicJobKind::ModalSynthesis => self.modal_synthesis,
            CinematicJobKind::AudioMaster => self.audio_master,
            CinematicJobKind::BundleVerification => self.bundle_verifier,
            CinematicJobKind::MuxDerivative => self.mux_adapter,
        }
    }

    fn validate(self, include_mux: bool) -> Result<(), CinematicJobPlanError> {
        for (name, identity) in [
            ("trajectory implementation", self.trajectory),
            ("render-shard implementation", self.render_shard),
            ("raw-merge implementation", self.raw_merge),
            ("temporal-finisher implementation", self.temporal_finish),
            ("image-sequence implementation", self.image_sequence),
            ("audio-controls implementation", self.audio_controls),
            ("audio-excitation implementation", self.audio_excitation),
            ("audio-resampling implementation", self.audio_resampling),
            ("modal-synthesis implementation", self.modal_synthesis),
            ("audio-master implementation", self.audio_master),
            ("bundle-verifier implementation", self.bundle_verifier),
        ] {
            require_nonzero(name, identity)?;
        }
        if include_mux {
            require_nonzero("mux-adapter implementation", self.mux_adapter)?;
        }
        Ok(())
    }
}

/// Bounds checked before retaining the graph or a resume snapshot.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CinematicJobLimits {
    pub max_nodes: u64,
    pub max_dependencies_per_node: u64,
    pub max_total_dependencies: u64,
    pub max_total_output_bytes: u64,
    pub max_snapshot_records: u64,
    pub max_snapshot_bytes: u64,
    pub max_events: u64,
}

impl CinematicJobLimits {
    fn validate(self) -> Result<(), CinematicJobPlanError> {
        for (name, value) in [
            ("maximum nodes", self.max_nodes),
            (
                "maximum dependencies per node",
                self.max_dependencies_per_node,
            ),
            ("maximum total dependencies", self.max_total_dependencies),
            ("maximum total output bytes", self.max_total_output_bytes),
            ("maximum snapshot records", self.max_snapshot_records),
            ("maximum snapshot bytes", self.max_snapshot_bytes),
            ("maximum events", self.max_events),
        ] {
            if value == 0 {
                return Err(CinematicJobPlanError::InvalidLimit(name));
            }
        }
        Ok(())
    }
}

/// One immutable node. Dependency indices are canonical, sorted, unique, and
/// always point to earlier nodes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CinematicJobNode {
    identity: ContentHash,
    expected_output_identity: ContentHash,
    kind: CinematicJobKind,
    artifact_kind: CinematicArtifactKind,
    dependencies: Vec<u32>,
    budget: CinematicNodeBudget,
    work: CinematicNodeWork,
}

impl CinematicJobNode {
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    #[must_use]
    pub const fn expected_output_identity(&self) -> ContentHash {
        self.expected_output_identity
    }

    #[must_use]
    pub const fn kind(&self) -> CinematicJobKind {
        self.kind
    }

    #[must_use]
    pub const fn artifact_kind(&self) -> CinematicArtifactKind {
        self.artifact_kind
    }

    #[must_use]
    pub fn dependencies(&self) -> &[u32] {
        &self.dependencies
    }

    #[must_use]
    pub const fn budget(&self) -> CinematicNodeBudget {
        self.budget
    }

    /// Exact logical render work carried by this node. Render counters are
    /// zero for non-render stages; finishing nodes additionally retain their
    /// canonical shot ordinal.
    #[must_use]
    pub const fn work(&self) -> CinematicNodeWork {
        self.work
    }
}

/// Canonical topological plan for one cinematic composition.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CinematicJobPlan {
    identity: ContentHash,
    configuration_identity: ContentHash,
    render_plan_identity: ContentHash,
    shot_plan_identity: ContentHash,
    shot_count: u64,
    nodes: Vec<CinematicJobNode>,
    total_dependencies: u64,
    total_work_units: u64,
    total_output_bytes: u64,
    limits: CinematicJobLimits,
}

impl CinematicJobPlan {
    /// Construct the graph from the admitted configuration, exact scene-bound
    /// prepared frames, and their immutable render topology. Camera-shot
    /// boundaries are derived from the scene's admitted exposures; callers
    /// cannot supply a second shot oracle.
    pub fn try_new(
        configuration: &CinematicConfig,
        render_plan: &EulerUniformRenderPlan,
        scene: &EulerCinematicScene<'_>,
        render_frames: &[EulerRenderFrameInput<'_>],
        bundle_expectation_identity: ContentHash,
        stages: CinematicStageIdentities,
        budgets: CinematicJobBudgets,
        limits: CinematicJobLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, CinematicJobPlanError> {
        checkpoint_plan(cx)?;
        if configuration.input().trajectory.identity() != render_plan.source_trajectory_identity() {
            return Err(CinematicJobPlanError::Incompatible(
                "configuration and render trajectory",
            ));
        }
        if configuration.input().timeline.identity() != render_plan.sequence_identity() {
            return Err(CinematicJobPlanError::Incompatible(
                "configuration timeline and render sequence",
            ));
        }
        limits.validate()?;
        require_nonzero(
            "bundle finalization expectation identity",
            bundle_expectation_identity,
        )?;
        enforce_limit(
            "shot frame inventory",
            u64::try_from(render_plan.frames().len())
                .map_err(|_| CinematicJobPlanError::ArithmeticOverflow("render frame count"))?,
            limits.max_nodes,
        )?;
        let (shots, shot_plan_identity) =
            derive_scene_shots(scene, render_plan, render_frames, cx)?;
        let include_mux = matches!(
            configuration.input().mux_request,
            CinematicMuxRequest::QuarantinedAdapter { .. }
        );
        stages.validate(include_mux)?;
        let shard_count = render_plan.shards().len() as u64;
        let segment_count = render_plan.segments().len() as u64;
        let _ = stage_total(budgets, shard_count, segment_count, include_mux, |budget| {
            budget.work_units
        })?;
        let total_output_bytes =
            stage_total(budgets, shard_count, segment_count, include_mux, |budget| {
                budget.max_output_bytes
            })?;
        enforce_limit(
            "total output bytes",
            total_output_bytes,
            limits.max_total_output_bytes,
        )?;
        let topology =
            RenderTopology::from_render_plan(render_plan, &shots, include_mux, limits, cx)?;
        let shot_count = u64::try_from(shots.len())
            .map_err(|_| CinematicJobPlanError::ArithmeticOverflow("shot count"))?;
        let sources = PlanSources {
            configuration_identity: configuration.composition_identity(),
            trajectory_partition_identity: configuration.trajectory_identity(),
            trajectory_artifact_identity: configuration.input().trajectory.identity(),
            image_identity: configuration.image_identity(),
            audio_identity: configuration.audio_identity(),
            mux_identity: configuration.mux_identity(),
            bundle_expectation_identity,
            shot_plan_identity,
            shot_count,
            include_mux,
        };
        build_plan(sources, topology, stages, budgets, limits, cx)
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    #[must_use]
    pub const fn configuration_identity(&self) -> ContentHash {
        self.configuration_identity
    }

    #[must_use]
    pub const fn render_plan_identity(&self) -> ContentHash {
        self.render_plan_identity
    }

    /// Identity of the exact camera-shot/frame partition used for progress.
    #[must_use]
    pub const fn shot_plan_identity(&self) -> ContentHash {
        self.shot_plan_identity
    }

    /// Number of camera shots in the canonical partition.
    #[must_use]
    pub const fn shot_count(&self) -> u64 {
        self.shot_count
    }

    #[must_use]
    pub fn nodes(&self) -> &[CinematicJobNode] {
        &self.nodes
    }

    #[must_use]
    pub const fn total_dependencies(&self) -> u64 {
        self.total_dependencies
    }

    #[must_use]
    pub const fn total_work_units(&self) -> u64 {
        self.total_work_units
    }

    #[must_use]
    pub const fn total_output_bytes(&self) -> u64 {
        self.total_output_bytes
    }

    #[must_use]
    pub const fn limits(&self) -> CinematicJobLimits {
        self.limits
    }

    /// Deterministic ready frontier for an external dispatcher. The built-in
    /// conductor remains sequential and deterministic; renderer/audio stages
    /// may internally use their existing bounded crews.
    #[must_use]
    pub fn ready_frontier(&self, completed: &BTreeSet<ContentHash>) -> Vec<&CinematicJobNode> {
        self.nodes
            .iter()
            .filter(|node| {
                !completed.contains(&node.identity)
                    && node.dependencies.iter().all(|dependency| {
                        self.nodes
                            .get(*dependency as usize)
                            .is_some_and(|parent| completed.contains(&parent.identity))
                    })
            })
            .collect()
    }
}

#[derive(Clone, Copy)]
struct PlanSources {
    configuration_identity: ContentHash,
    trajectory_partition_identity: ContentHash,
    trajectory_artifact_identity: ContentHash,
    image_identity: ContentHash,
    audio_identity: ContentHash,
    mux_identity: ContentHash,
    bundle_expectation_identity: ContentHash,
    shot_plan_identity: ContentHash,
    shot_count: u64,
    include_mux: bool,
}

#[derive(Clone)]
struct TopologyShard {
    ordinal: u64,
    logical_identity: ContentHash,
    frame_ordinal: u64,
    segment_index: u64,
    tile_start: u64,
    tile_count: u64,
    samples_per_tile: u64,
    path_count: u64,
}

#[derive(Clone)]
struct TopologySegment {
    frame_ordinal: u64,
    frame_position: u64,
    segment_index: u64,
    shot_ordinal: u64,
    shot_identity: ContentHash,
    frame_identity: ContentHash,
    shard_indices: Vec<usize>,
    neighbor_segment_indices: Vec<usize>,
}

#[derive(Clone)]
struct RenderTopology {
    plan_identity: ContentHash,
    sequence_identity: ContentHash,
    source_trajectory_identity: ContentHash,
    shards: Vec<TopologyShard>,
    segments: Vec<TopologySegment>,
}

fn derive_scene_shots(
    scene: &EulerCinematicScene<'_>,
    plan: &EulerUniformRenderPlan,
    render_frames: &[EulerRenderFrameInput<'_>],
    cx: &Cx<'_>,
) -> Result<(Vec<DerivedShotRange>, ContentHash), CinematicJobPlanError> {
    if scene.scene_identity() != plan.scene_identity()
        || scene.source_configuration_identity() != plan.source_configuration_identity()
        || scene.source_trajectory_identity() != plan.source_trajectory_identity()
    {
        return Err(CinematicJobPlanError::Incompatible(
            "scene and render-plan authority",
        ));
    }
    if render_frames.len() != plan.frames().len() {
        return Err(CinematicJobPlanError::Incompatible(
            "render frame inventory",
        ));
    }

    let mut frames_by_ordinal = BTreeMap::new();
    for input in render_frames {
        checkpoint_plan(cx)?;
        if frames_by_ordinal
            .insert(input.frame_ordinal(), *input)
            .is_some()
        {
            return Err(CinematicJobPlanError::Incompatible(
                "duplicate render frame ordinal",
            ));
        }
    }

    let mut shots = Vec::<DerivedShotRange>::new();
    shots
        .try_reserve_exact(plan.frames().len())
        .map_err(|_| CinematicJobPlanError::Capacity("derived shot ranges"))?;
    let mut seen_shot_ids = BTreeSet::new();
    for (frame_position, planned_frame) in plan.frames().iter().enumerate() {
        checkpoint_plan(cx)?;
        let input = frames_by_ordinal
            .get(&planned_frame.frame_ordinal())
            .copied()
            .ok_or(CinematicJobPlanError::Incompatible(
                "missing render frame input",
            ))?;
        let prepared = input.prepared();
        let prepared_segment_count = u64::try_from(prepared.segments().len())
            .map_err(|_| CinematicJobPlanError::ArithmeticOverflow("prepared segment count"))?;
        if prepared_segment_count != planned_frame.segment_count() {
            return Err(CinematicJobPlanError::Incompatible(
                "prepared frame segment inventory",
            ));
        }
        let first_segment = usize::try_from(planned_frame.first_segment())
            .map_err(|_| CinematicJobPlanError::ArithmeticOverflow("first segment index"))?;
        let canonical_frame_position = u64::try_from(frame_position)
            .map_err(|_| CinematicJobPlanError::ArithmeticOverflow("frame position"))?;
        let mut frame_shot_id = None;
        for segment_index in 0..prepared.segments().len() {
            if segment_index.is_multiple_of(256) {
                checkpoint_plan(cx)?;
            }
            let planned_segment_index = first_segment.checked_add(segment_index).ok_or(
                CinematicJobPlanError::ArithmeticOverflow("planned segment index"),
            )?;
            let planned_segment = plan.segments().get(planned_segment_index).ok_or(
                CinematicJobPlanError::Incompatible("planned frame segment inventory"),
            )?;
            let canonical_segment_index = u64::try_from(segment_index)
                .map_err(|_| CinematicJobPlanError::ArithmeticOverflow("segment index"))?;
            let expected_frame_identity =
                euler_render_checkpoint_frame_identity(prepared, segment_index)
                    .map_err(|_| CinematicJobPlanError::Incompatible("prepared frame identity"))?;
            if planned_segment.frame_ordinal() != planned_frame.frame_ordinal()
                || planned_segment.frame_position() != canonical_frame_position
                || planned_segment.segment_index() != canonical_segment_index
                || planned_segment.frame_identity() != expected_frame_identity
            {
                return Err(CinematicJobPlanError::Incompatible(
                    "render plan and prepared frame binding",
                ));
            }
            let (shutter, cut_side) = scene
                .prepared_segment_shard_binding(prepared, segment_index)
                .map_err(|_| CinematicJobPlanError::Incompatible("prepared scene binding"))?;
            let shot_id = scene
                .camera()
                .admit_shutter(cx, shutter, cut_side)
                .map_err(|_| CinematicJobPlanError::Incompatible("camera shot admission"))?
                .shot_id();
            if shot_id == 0 || frame_shot_id.is_some_and(|prior| prior != shot_id) {
                return Err(CinematicJobPlanError::Incompatible(
                    "frame crosses camera-shot boundary",
                ));
            }
            frame_shot_id = Some(shot_id);
        }
        let shot_id =
            frame_shot_id.ok_or(CinematicJobPlanError::Incompatible("empty prepared frame"))?;
        if let Some(last) = shots.last_mut()
            && last.shot_id == shot_id
        {
            last.frame_count = checked_add(last.frame_count, 1, "shot frame count")?;
            continue;
        }
        if !seen_shot_ids.insert(shot_id) {
            return Err(CinematicJobPlanError::Incompatible(
                "non-contiguous camera-shot re-entry",
            ));
        }
        let mut shot_hasher =
            DomainHasher::new("org.frankensim.fs-euler-disc-e2e.cinematic-scene-shot.v1");
        shot_hasher.update(&CINEMATIC_JOB_SCHEMA_VERSION.to_le_bytes());
        shot_hasher.update(scene.source_configuration_identity().as_bytes());
        shot_hasher.update(&shot_id.to_le_bytes());
        shots.push(DerivedShotRange {
            identity: shot_hasher.finalize(),
            shot_id,
            first_frame_position: canonical_frame_position,
            frame_count: 1,
        });
    }
    let identity =
        validate_and_identify_shots(&shots, plan.frames().len(), || checkpoint_plan(cx))?;
    Ok((shots, identity))
}

fn validate_and_identify_shots(
    shots: &[DerivedShotRange],
    frame_count: usize,
    mut checkpoint: impl FnMut() -> Result<(), CinematicJobPlanError>,
) -> Result<ContentHash, CinematicJobPlanError> {
    if shots.is_empty() {
        return Err(CinematicJobPlanError::Incompatible("empty shot plan"));
    }
    let frame_count = u64::try_from(frame_count)
        .map_err(|_| CinematicJobPlanError::ArithmeticOverflow("render frame count"))?;
    let shot_count = u64::try_from(shots.len())
        .map_err(|_| CinematicJobPlanError::ArithmeticOverflow("shot count"))?;
    if shot_count > frame_count {
        return Err(CinematicJobPlanError::Incompatible(
            "more shots than render frames",
        ));
    }
    let mut expected_start = 0_u64;
    let mut hasher = DomainHasher::new(CINEMATIC_SHOT_PLAN_IDENTITY_DOMAIN);
    hasher.update(&CINEMATIC_JOB_SCHEMA_VERSION.to_le_bytes());
    hasher.update(&shot_count.to_le_bytes());
    for (shot_index, shot) in shots.iter().enumerate() {
        if shot_index.is_multiple_of(256) {
            checkpoint()?;
        }
        require_nonzero("shot identity", shot.identity)?;
        if shot.frame_count == 0 || shot.first_frame_position != expected_start {
            return Err(CinematicJobPlanError::Incompatible(
                "non-canonical shot frame partition",
            ));
        }
        expected_start = checked_add(expected_start, shot.frame_count, "shot frame partition")?;
        hasher.update(shot.identity.as_bytes());
        hasher.update(&shot.first_frame_position.to_le_bytes());
        hasher.update(&shot.frame_count.to_le_bytes());
    }
    if expected_start != frame_count {
        return Err(CinematicJobPlanError::Incompatible(
            "shot plan does not cover render frames",
        ));
    }
    Ok(hasher.finalize())
}

fn shot_for_frame_position(
    shots: &[DerivedShotRange],
    frame_position: u64,
) -> Result<(usize, DerivedShotRange), CinematicJobPlanError> {
    let insertion = shots.partition_point(|shot| shot.first_frame_position <= frame_position);
    let shot_index = insertion
        .checked_sub(1)
        .ok_or(CinematicJobPlanError::Incompatible(
            "segment frame position outside shot plan",
        ))?;
    let shot = shots[shot_index];
    let shot_end = checked_add(
        shot.first_frame_position,
        shot.frame_count,
        "shot frame-position end",
    )?;
    if frame_position >= shot_end {
        return Err(CinematicJobPlanError::Incompatible(
            "segment frame position outside shot plan",
        ));
    }
    Ok((shot_index, shot))
}

impl RenderTopology {
    fn from_render_plan(
        plan: &EulerUniformRenderPlan,
        shots: &[DerivedShotRange],
        include_mux: bool,
        limits: CinematicJobLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, CinematicJobPlanError> {
        limits.validate()?;
        let shard_count = plan.shards().len() as u64;
        let segment_count = plan.segments().len() as u64;
        let node_count = checked_add(
            checked_add(
                shard_count,
                checked_mul(segment_count, 2, "segment stages")?,
                "render nodes",
            )?,
            8 + u64::from(include_mux),
            "maximum topology node count",
        )?;
        enforce_limit("nodes", node_count, limits.max_nodes)?;
        enforce_limit("snapshot records", node_count, limits.max_snapshot_records)?;
        let snapshot_bytes = checked_add(
            SNAPSHOT_HEADER_BYTES,
            checked_mul(SNAPSHOT_RECORD_BYTES, node_count, "snapshot record bytes")?,
            "snapshot bytes",
        )?;
        enforce_limit("snapshot bytes", snapshot_bytes, limits.max_snapshot_bytes)?;
        enforce_limit(
            "events",
            checked_mul(node_count, MAX_EVENTS_PER_NODE, "event bound")?,
            limits.max_events,
        )?;
        let _ = node_capacity(node_count)?;
        enforce_limit(
            "dependencies per node",
            segment_count,
            limits.max_dependencies_per_node,
        )?;
        let mut segment_shots = Vec::new();
        segment_shots
            .try_reserve_exact(plan.segments().len())
            .map_err(|_| CinematicJobPlanError::Capacity("segment shot assignments"))?;
        for segment in plan.segments() {
            checkpoint_plan(cx)?;
            let (shot_index, shot) = shot_for_frame_position(shots, segment.frame_position())?;
            let shot_ordinal = u64::try_from(shot_index)
                .map_err(|_| CinematicJobPlanError::ArithmeticOverflow("shot ordinal"))?;
            segment_shots.push((shot_ordinal, shot.identity));
        }
        let mut neighbor_windows = Vec::new();
        neighbor_windows
            .try_reserve_exact(plan.segments().len())
            .map_err(|_| CinematicJobPlanError::Capacity("finishing neighbor windows"))?;
        let radius = u64::from(plan.finishing_neighbor_radius());
        let mut finishing_dependencies = 0_u64;
        let mut largest_merge = 0_u64;
        let nonfinishing_dependencies = checked_add(
            checked_mul(shard_count, 2, "render dependency total")?,
            checked_add(
                segment_count,
                7 + u64::from(include_mux),
                "sequence and A/V dependency total",
            )?,
            "non-finishing dependency total",
        )?;
        for (target_index, target) in plan.segments().iter().enumerate() {
            checkpoint_plan(cx)?;
            let shot_index = usize::try_from(segment_shots[target_index].0)
                .map_err(|_| CinematicJobPlanError::ArithmeticOverflow("shot index"))?;
            let shot = shots[shot_index];
            let shot_end = checked_add(
                shot.first_frame_position,
                shot.frame_count,
                "shot frame-position end",
            )?;
            let first_neighbor_frame = target
                .frame_position()
                .saturating_sub(radius)
                .max(shot.first_frame_position);
            let last_neighbor_frame = target
                .frame_position()
                .saturating_add(radius)
                .min(shot_end - 1);
            let mut neighbor_count = 0_u64;
            let mut frame_position = first_neighbor_frame;
            loop {
                if (frame_position - first_neighbor_frame).is_multiple_of(256) {
                    checkpoint_plan(cx)?;
                }
                if frame_position != target.frame_position() {
                    let frame_index = usize::try_from(frame_position).map_err(|_| {
                        CinematicJobPlanError::ArithmeticOverflow("neighbor frame index")
                    })?;
                    let frame = plan.frames().get(frame_index).ok_or(
                        CinematicJobPlanError::Incompatible("neighbor frame position"),
                    )?;
                    neighbor_count = checked_add(
                        neighbor_count,
                        frame.segment_count(),
                        "finishing neighbor count",
                    )?;
                    let prospective_finish =
                        checked_add(neighbor_count, 1, "finishing dependencies")?;
                    enforce_limit(
                        "dependencies per node",
                        prospective_finish,
                        limits.max_dependencies_per_node,
                    )?;
                    enforce_limit(
                        "total dependencies",
                        checked_add(
                            nonfinishing_dependencies,
                            checked_add(
                                finishing_dependencies,
                                prospective_finish,
                                "prospective finishing dependency total",
                            )?,
                            "prospective topology dependency total",
                        )?,
                        limits.max_total_dependencies,
                    )?;
                }
                if frame_position == last_neighbor_frame {
                    break;
                }
                frame_position += 1;
            }
            let finish_count = checked_add(neighbor_count, 1, "finishing dependencies")?;
            enforce_limit(
                "dependencies per node",
                finish_count,
                limits.max_dependencies_per_node,
            )?;
            finishing_dependencies = checked_add(
                finishing_dependencies,
                finish_count,
                "finishing dependency total",
            )?;
            enforce_limit(
                "total dependencies",
                checked_add(
                    nonfinishing_dependencies,
                    finishing_dependencies,
                    "topology dependency total",
                )?,
                limits.max_total_dependencies,
            )?;
            largest_merge = largest_merge.max(target.shard_count());
            neighbor_windows.push((first_neighbor_frame, last_neighbor_frame, neighbor_count));
        }
        enforce_limit(
            "dependencies per node",
            largest_merge,
            limits.max_dependencies_per_node,
        )?;
        let total_dependencies = checked_add(
            nonfinishing_dependencies,
            finishing_dependencies,
            "topology dependency total",
        )?;
        enforce_limit(
            "total dependencies",
            total_dependencies,
            limits.max_total_dependencies,
        )?;

        let mut shards = Vec::new();
        shards
            .try_reserve_exact(plan.shards().len())
            .map_err(|_| CinematicJobPlanError::Capacity("render topology shards"))?;
        for shard in plan.shards() {
            checkpoint_plan(cx)?;
            shards.push(TopologyShard {
                ordinal: shard.shard_ordinal(),
                logical_identity: shard.logical_shard_identity(),
                frame_ordinal: shard.frame_ordinal(),
                segment_index: shard.segment_index(),
                tile_start: shard.tile_start(),
                tile_count: shard
                    .tile_end()
                    .checked_sub(shard.tile_start())
                    .ok_or(CinematicJobPlanError::Incompatible("shard tile range"))?,
                samples_per_tile: u64::from(
                    shard
                        .sample_end()
                        .checked_sub(shard.sample_start())
                        .ok_or(CinematicJobPlanError::Incompatible("shard sample range"))?,
                ),
                path_count: shard.path_count(),
            });
        }
        let mut segments = Vec::new();
        segments
            .try_reserve_exact(plan.segments().len())
            .map_err(|_| CinematicJobPlanError::Capacity("render topology segments"))?;
        for (segment_index, segment) in plan.segments().iter().enumerate() {
            checkpoint_plan(cx)?;
            let first = usize::try_from(segment.first_shard())
                .map_err(|_| CinematicJobPlanError::ArithmeticOverflow("first shard index"))?;
            let count = usize::try_from(segment.shard_count())
                .map_err(|_| CinematicJobPlanError::ArithmeticOverflow("shard count"))?;
            let end = first
                .checked_add(count)
                .ok_or(CinematicJobPlanError::ArithmeticOverflow(
                    "segment shard range",
                ))?;
            if end > shards.len() {
                return Err(CinematicJobPlanError::Incompatible("segment shard range"));
            }
            let mut shard_indices = Vec::new();
            shard_indices
                .try_reserve_exact(count)
                .map_err(|_| CinematicJobPlanError::Capacity("segment shard indices"))?;
            shard_indices.extend(first..end);
            let (first_neighbor_frame, last_neighbor_frame, neighbor_count) =
                neighbor_windows[segment_index];
            let neighbor_capacity = usize::try_from(neighbor_count).map_err(|_| {
                CinematicJobPlanError::ArithmeticOverflow("finishing neighbor capacity")
            })?;
            let mut neighbor_segment_indices = Vec::new();
            neighbor_segment_indices
                .try_reserve_exact(neighbor_capacity)
                .map_err(|_| CinematicJobPlanError::Capacity("finishing neighbor indices"))?;
            let mut frame_position = first_neighbor_frame;
            loop {
                if (frame_position - first_neighbor_frame).is_multiple_of(256) {
                    checkpoint_plan(cx)?;
                }
                if frame_position != segment.frame_position() {
                    let frame_index = usize::try_from(frame_position).map_err(|_| {
                        CinematicJobPlanError::ArithmeticOverflow("neighbor frame index")
                    })?;
                    let frame = plan.frames().get(frame_index).ok_or(
                        CinematicJobPlanError::Incompatible("neighbor frame position"),
                    )?;
                    let first = usize::try_from(frame.first_segment()).map_err(|_| {
                        CinematicJobPlanError::ArithmeticOverflow("neighbor segment start")
                    })?;
                    let end = frame
                        .first_segment()
                        .checked_add(frame.segment_count())
                        .ok_or(CinematicJobPlanError::ArithmeticOverflow(
                            "neighbor segment end",
                        ))?;
                    let end = usize::try_from(end).map_err(|_| {
                        CinematicJobPlanError::ArithmeticOverflow("neighbor segment end")
                    })?;
                    for candidate_index in first..end {
                        if candidate_index.is_multiple_of(256) {
                            checkpoint_plan(cx)?;
                        }
                        neighbor_segment_indices.push(candidate_index);
                    }
                }
                if frame_position == last_neighbor_frame {
                    break;
                }
                frame_position += 1;
            }
            debug_assert_eq!(neighbor_segment_indices.len(), neighbor_capacity);
            let (shot_ordinal, shot_identity) = segment_shots[segment_index];
            segments.push(TopologySegment {
                frame_ordinal: segment.frame_ordinal(),
                frame_position: segment.frame_position(),
                segment_index: segment.segment_index(),
                shot_ordinal,
                shot_identity,
                frame_identity: segment.frame_identity(),
                shard_indices,
                neighbor_segment_indices,
            });
        }
        Ok(Self {
            plan_identity: plan.plan_identity(),
            sequence_identity: plan.sequence_identity(),
            source_trajectory_identity: plan.source_trajectory_identity(),
            shards,
            segments,
        })
    }
}

fn build_plan(
    sources: PlanSources,
    topology: RenderTopology,
    stages: CinematicStageIdentities,
    budgets: CinematicJobBudgets,
    limits: CinematicJobLimits,
    cx: &Cx<'_>,
) -> Result<CinematicJobPlan, CinematicJobPlanError> {
    checkpoint_plan(cx)?;
    limits.validate()?;
    stages.validate(sources.include_mux)?;
    for (name, identity) in [
        ("configuration identity", sources.configuration_identity),
        (
            "trajectory partition identity",
            sources.trajectory_partition_identity,
        ),
        (
            "trajectory artifact identity",
            sources.trajectory_artifact_identity,
        ),
        ("image identity", sources.image_identity),
        ("audio identity", sources.audio_identity),
        (
            "bundle finalization expectation identity",
            sources.bundle_expectation_identity,
        ),
        ("shot plan identity", sources.shot_plan_identity),
        ("render plan identity", topology.plan_identity),
        ("render sequence identity", topology.sequence_identity),
        (
            "render trajectory identity",
            topology.source_trajectory_identity,
        ),
    ] {
        require_nonzero(name, identity)?;
    }
    if sources.include_mux {
        require_nonzero("mux identity", sources.mux_identity)?;
    }
    if sources.shot_count == 0 {
        return Err(CinematicJobPlanError::Incompatible("empty shot plan"));
    }
    if sources.trajectory_artifact_identity != topology.source_trajectory_identity {
        return Err(CinematicJobPlanError::Incompatible(
            "trajectory source identity",
        ));
    }
    validate_topology(&topology)?;
    validate_topology_shots(&topology, sources.shot_count)?;

    let node_count = checked_add(
        checked_add(
            u64::try_from(topology.shards.len())
                .map_err(|_| CinematicJobPlanError::ArithmeticOverflow("shard node count"))?,
            checked_mul(
                u64::try_from(topology.segments.len())
                    .map_err(|_| CinematicJobPlanError::ArithmeticOverflow("segment node count"))?,
                2,
                "segment stages",
            )?,
            "render nodes",
        )?,
        8 + u64::from(sources.include_mux),
        "total node count",
    )?;
    enforce_limit("nodes", node_count, limits.max_nodes)?;
    enforce_limit("snapshot records", node_count, limits.max_snapshot_records)?;
    let snapshot_bytes = checked_add(
        SNAPSHOT_HEADER_BYTES,
        checked_mul(SNAPSHOT_RECORD_BYTES, node_count, "snapshot record bytes")?,
        "snapshot bytes",
    )?;
    enforce_limit("snapshot bytes", snapshot_bytes, limits.max_snapshot_bytes)?;
    let worst_events = checked_mul(node_count, MAX_EVENTS_PER_NODE, "event bound")?;
    enforce_limit("events", worst_events, limits.max_events)?;

    let shard_count = topology.shards.len() as u64;
    let segment_count = topology.segments.len() as u64;
    let mut finishing_dependencies = 0_u64;
    let mut largest_merge = 0_u64;
    let mut largest_finish = 0_u64;
    for segment in &topology.segments {
        let merge_count = segment.shard_indices.len() as u64;
        let finish_count = checked_add(
            segment.neighbor_segment_indices.len() as u64,
            1,
            "finishing dependency count",
        )?;
        largest_merge = largest_merge.max(merge_count);
        largest_finish = largest_finish.max(finish_count);
        finishing_dependencies = checked_add(
            finishing_dependencies,
            finish_count,
            "finishing dependency total",
        )?;
    }
    let total_dependencies_preflight = checked_add(
        checked_add(
            checked_mul(shard_count, 2, "render dependency total")?,
            finishing_dependencies,
            "render and finishing dependency total",
        )?,
        checked_add(
            segment_count,
            7 + u64::from(sources.include_mux),
            "sequence and A/V dependency total",
        )?,
        "total dependency preflight",
    )?;
    enforce_limit(
        "total dependencies",
        total_dependencies_preflight,
        limits.max_total_dependencies,
    )?;
    for observed in [largest_merge, largest_finish, segment_count, 2] {
        enforce_limit(
            "dependencies per node",
            observed,
            limits.max_dependencies_per_node,
        )?;
    }

    let total_work_preflight = stage_total(
        budgets,
        shard_count,
        segment_count,
        sources.include_mux,
        |budget| budget.work_units,
    )?;
    let total_output_preflight = stage_total(
        budgets,
        shard_count,
        segment_count,
        sources.include_mux,
        |budget| budget.max_output_bytes,
    )?;
    enforce_limit(
        "total output bytes",
        total_output_preflight,
        limits.max_total_output_bytes,
    )?;

    let node_capacity = node_capacity(node_count)?;
    let mut nodes = Vec::new();
    nodes
        .try_reserve_exact(node_capacity)
        .map_err(|_| CinematicJobPlanError::Capacity("cinematic nodes"))?;

    let trajectory_index = push_node(
        &mut nodes,
        CinematicJobKind::Trajectory,
        sources.trajectory_partition_identity,
        stages.trajectory,
        &[sources.trajectory_artifact_identity],
        &[],
        budgets.trajectory,
        CinematicNodeWork::NONE,
    )?;

    let mut shard_nodes = Vec::new();
    shard_nodes
        .try_reserve_exact(topology.shards.len())
        .map_err(|_| CinematicJobPlanError::Capacity("shard node index"))?;
    for shard in &topology.shards {
        checkpoint_plan(cx)?;
        let mut local = Vec::with_capacity(32 + 8 * 3);
        local.extend_from_slice(topology.plan_identity.as_bytes());
        local.extend_from_slice(&shard.ordinal.to_le_bytes());
        local.extend_from_slice(&shard.frame_ordinal.to_le_bytes());
        local.extend_from_slice(&shard.segment_index.to_le_bytes());
        let local_identity = hash_domain(
            "org.frankensim.fs-euler-disc-e2e.cinematic-render-shard-local.v1",
            &local,
        );
        shard_nodes.push(push_node(
            &mut nodes,
            CinematicJobKind::RenderShard {
                shard_ordinal: shard.ordinal,
            },
            sources.image_identity,
            stages.render_shard,
            &[shard.logical_identity, local_identity],
            &[trajectory_index],
            budgets.render_shard,
            CinematicNodeWork::render_shard(
                shard.frame_ordinal,
                shard.segment_index,
                shard.tile_start,
                shard.tile_count,
                shard.samples_per_tile,
                shard.path_count,
            )?,
        )?);
    }

    let mut merge_nodes = Vec::new();
    merge_nodes
        .try_reserve_exact(topology.segments.len())
        .map_err(|_| CinematicJobPlanError::Capacity("merge node index"))?;
    for segment in &topology.segments {
        checkpoint_plan(cx)?;
        let dependencies = segment
            .shard_indices
            .iter()
            .map(|index| shard_nodes[*index])
            .collect::<Vec<_>>();
        merge_nodes.push(push_node(
            &mut nodes,
            CinematicJobKind::MergeRawSegment {
                frame_ordinal: segment.frame_ordinal,
                segment_index: segment.segment_index,
            },
            sources.image_identity,
            stages.raw_merge,
            &[topology.plan_identity, segment.frame_identity],
            &dependencies,
            budgets.raw_merge,
            CinematicNodeWork::NONE,
        )?);
    }

    let mut finish_nodes = Vec::new();
    finish_nodes
        .try_reserve_exact(topology.segments.len())
        .map_err(|_| CinematicJobPlanError::Capacity("finish node index"))?;
    for (index, segment) in topology.segments.iter().enumerate() {
        checkpoint_plan(cx)?;
        let mut dependencies = Vec::new();
        dependencies
            .try_reserve_exact(segment.neighbor_segment_indices.len() + 1)
            .map_err(|_| CinematicJobPlanError::Capacity("finishing dependencies"))?;
        dependencies.push(merge_nodes[index]);
        dependencies.extend(
            segment
                .neighbor_segment_indices
                .iter()
                .map(|neighbor| merge_nodes[*neighbor]),
        );
        dependencies.sort_unstable();
        dependencies.dedup();
        finish_nodes.push(push_node(
            &mut nodes,
            CinematicJobKind::FinishSegment {
                frame_ordinal: segment.frame_ordinal,
                segment_index: segment.segment_index,
            },
            sources.image_identity,
            stages.temporal_finish,
            &[
                topology.plan_identity,
                segment.frame_identity,
                hash_u64("cinematic-frame-position-v1", segment.frame_position),
                segment.shot_identity,
            ],
            &dependencies,
            budgets.temporal_finish,
            CinematicNodeWork::finished_segment(segment.shot_ordinal),
        )?);
    }

    let image_sequence_index = push_node(
        &mut nodes,
        CinematicJobKind::ImageSequence,
        sources.image_identity,
        stages.image_sequence,
        &[topology.plan_identity, topology.sequence_identity],
        &finish_nodes,
        budgets.image_sequence,
        CinematicNodeWork::NONE,
    )?;
    let audio_controls_index = push_node(
        &mut nodes,
        CinematicJobKind::AudioControls,
        sources.audio_identity,
        stages.audio_controls,
        &[topology.sequence_identity],
        &[trajectory_index],
        budgets.audio_controls,
        CinematicNodeWork::NONE,
    )?;
    let audio_excitation_index = push_node(
        &mut nodes,
        CinematicJobKind::AudioExcitation,
        sources.audio_identity,
        stages.audio_excitation,
        &[],
        &[audio_controls_index],
        budgets.audio_excitation,
        CinematicNodeWork::NONE,
    )?;
    let audio_resampling_index = push_node(
        &mut nodes,
        CinematicJobKind::AudioResampling,
        sources.audio_identity,
        stages.audio_resampling,
        &[],
        &[audio_excitation_index],
        budgets.audio_resampling,
        CinematicNodeWork::NONE,
    )?;
    let modal_synthesis_index = push_node(
        &mut nodes,
        CinematicJobKind::ModalSynthesis,
        sources.audio_identity,
        stages.modal_synthesis,
        &[],
        &[audio_resampling_index],
        budgets.modal_synthesis,
        CinematicNodeWork::NONE,
    )?;
    let audio_master_index = push_node(
        &mut nodes,
        CinematicJobKind::AudioMaster,
        sources.audio_identity,
        stages.audio_master,
        &[],
        &[modal_synthesis_index],
        budgets.audio_master,
        CinematicNodeWork::NONE,
    )?;
    let bundle_index = push_node(
        &mut nodes,
        CinematicJobKind::BundleVerification,
        hash_pair(
            "org.frankensim.fs-euler-disc-e2e.cinematic-av-partitions.v1",
            sources.image_identity,
            sources.audio_identity,
        ),
        stages.bundle_verifier,
        &[
            topology.plan_identity,
            topology.sequence_identity,
            sources.bundle_expectation_identity,
        ],
        &[image_sequence_index, audio_master_index],
        budgets.bundle_verification,
        CinematicNodeWork::NONE,
    )?;
    if sources.include_mux {
        let _ = push_node(
            &mut nodes,
            CinematicJobKind::MuxDerivative,
            sources.mux_identity,
            stages.mux_adapter,
            &[],
            &[bundle_index],
            budgets.mux_derivative,
            CinematicNodeWork::NONE,
        )?;
    }

    debug_assert_eq!(nodes.len(), node_capacity);
    validate_node_work_totals(&nodes)?;
    let (total_dependencies, total_work_units, total_output_bytes) = plan_totals(&nodes)?;
    debug_assert_eq!(total_dependencies, total_dependencies_preflight);
    debug_assert_eq!(total_work_units, total_work_preflight);
    debug_assert_eq!(total_output_bytes, total_output_preflight);
    enforce_limit(
        "total dependencies",
        total_dependencies,
        limits.max_total_dependencies,
    )?;
    enforce_limit(
        "total output bytes",
        total_output_bytes,
        limits.max_total_output_bytes,
    )?;
    for node in &nodes {
        enforce_limit(
            "dependencies per node",
            node.dependencies.len() as u64,
            limits.max_dependencies_per_node,
        )?;
    }
    let identity = plan_identity(
        sources.configuration_identity,
        topology.plan_identity,
        sources.shot_plan_identity,
        &nodes,
        limits,
        total_dependencies,
        total_work_units,
        total_output_bytes,
    );
    Ok(CinematicJobPlan {
        identity,
        configuration_identity: sources.configuration_identity,
        render_plan_identity: topology.plan_identity,
        shot_plan_identity: sources.shot_plan_identity,
        shot_count: sources.shot_count,
        nodes,
        total_dependencies,
        total_work_units,
        total_output_bytes,
        limits,
    })
}

fn stage_total(
    budgets: CinematicJobBudgets,
    shard_count: u64,
    segment_count: u64,
    include_mux: bool,
    value: impl Fn(CinematicNodeBudget) -> u64,
) -> Result<u64, CinematicJobPlanError> {
    let mut total = checked_mul(
        value(budgets.render_shard),
        shard_count,
        "render stage total",
    )?;
    total = checked_add(
        total,
        checked_mul(value(budgets.raw_merge), segment_count, "merge stage total")?,
        "render and merge stage total",
    )?;
    total = checked_add(
        total,
        checked_mul(
            value(budgets.temporal_finish),
            segment_count,
            "finishing stage total",
        )?,
        "render stage total",
    )?;
    for budget in [
        budgets.trajectory,
        budgets.image_sequence,
        budgets.audio_controls,
        budgets.audio_excitation,
        budgets.audio_resampling,
        budgets.modal_synthesis,
        budgets.audio_master,
        budgets.bundle_verification,
    ] {
        total = checked_add(total, value(budget), "fixed stage total")?;
    }
    if include_mux {
        total = checked_add(total, value(budgets.mux_derivative), "mux stage total")?;
    }
    Ok(total)
}

fn validate_topology(topology: &RenderTopology) -> Result<(), CinematicJobPlanError> {
    if topology.shards.is_empty() || topology.segments.is_empty() {
        return Err(CinematicJobPlanError::Incompatible("empty render topology"));
    }
    for (index, shard) in topology.shards.iter().enumerate() {
        if shard.ordinal != index as u64 {
            return Err(CinematicJobPlanError::Incompatible(
                "non-canonical shard ordinal",
            ));
        }
        require_nonzero("logical shard identity", shard.logical_identity)?;
    }
    let mut previous = None;
    for (index, segment) in topology.segments.iter().enumerate() {
        let key = (segment.frame_ordinal, segment.segment_index);
        if previous.is_some_and(|value| value >= key) {
            return Err(CinematicJobPlanError::Incompatible(
                "non-canonical segment order",
            ));
        }
        previous = Some(key);
        require_nonzero("frame identity", segment.frame_identity)?;
        require_nonzero("shot identity", segment.shot_identity)?;
        if segment.shard_indices.is_empty()
            || segment
                .shard_indices
                .iter()
                .any(|shard| *shard >= topology.shards.len())
            || segment.neighbor_segment_indices.iter().any(|neighbor| {
                *neighbor >= topology.segments.len()
                    || *neighbor == index
                    || topology.segments[*neighbor].shot_ordinal != segment.shot_ordinal
            })
        {
            return Err(CinematicJobPlanError::Incompatible(
                "segment dependency topology",
            ));
        }
        let mut canonical_shards = segment.shard_indices.clone();
        canonical_shards.sort_unstable();
        canonical_shards.dedup();
        if canonical_shards != segment.shard_indices
            || segment.shard_indices.iter().any(|shard_index| {
                let shard = &topology.shards[*shard_index];
                shard.frame_ordinal != segment.frame_ordinal
                    || shard.segment_index != segment.segment_index
            })
        {
            return Err(CinematicJobPlanError::Incompatible(
                "non-canonical segment shards",
            ));
        }
        let mut canonical = segment.neighbor_segment_indices.clone();
        canonical.sort_unstable();
        canonical.dedup();
        if canonical != segment.neighbor_segment_indices {
            return Err(CinematicJobPlanError::Incompatible(
                "non-canonical finishing neighbors",
            ));
        }
    }
    Ok(())
}

fn validate_topology_shots(
    topology: &RenderTopology,
    expected_shot_count: u64,
) -> Result<(), CinematicJobPlanError> {
    let mut identities = BTreeMap::<u64, ContentHash>::new();
    for segment in &topology.segments {
        if segment.shot_ordinal >= expected_shot_count {
            return Err(CinematicJobPlanError::Incompatible(
                "shot ordinal outside shot plan",
            ));
        }
        match identities.insert(segment.shot_ordinal, segment.shot_identity) {
            Some(previous) if previous != segment.shot_identity => {
                return Err(CinematicJobPlanError::Incompatible(
                    "inconsistent shot identity",
                ));
            }
            _ => {}
        }
    }
    if identities.len() as u64 != expected_shot_count
        || identities
            .keys()
            .copied()
            .enumerate()
            .any(|(expected, observed)| expected as u64 != observed)
    {
        return Err(CinematicJobPlanError::Incompatible(
            "non-canonical topology shot coverage",
        ));
    }
    Ok(())
}

fn push_node(
    nodes: &mut Vec<CinematicJobNode>,
    kind: CinematicJobKind,
    partition_identity: ContentHash,
    implementation_identity: ContentHash,
    local_inputs: &[ContentHash],
    dependencies: &[u32],
    budget: CinematicNodeBudget,
    work: CinematicNodeWork,
) -> Result<u32, CinematicJobPlanError> {
    require_nonzero("node partition identity", partition_identity)?;
    require_nonzero("node implementation identity", implementation_identity)?;
    for input in local_inputs {
        require_nonzero("node local input", *input)?;
    }
    let work_matches_kind = match kind {
        CinematicJobKind::RenderShard { .. } => {
            work.render_frame_ordinal.is_some()
                && work.render_segment_index.is_some()
                && work.render_tiles > 0
                && work.samples_per_tile > 0
                && work.render_paths > 0
                && work.shot_ordinal.is_none()
        }
        CinematicJobKind::FinishSegment { .. } => {
            work.render_frame_ordinal.is_none()
                && work.render_segment_index.is_none()
                && work.render_tile_start == 0
                && work.render_tiles == 0
                && work.samples_per_tile == 0
                && work.render_paths == 0
                && work.shot_ordinal.is_some()
        }
        _ => work == CinematicNodeWork::NONE,
    };
    if !work_matches_kind {
        return Err(CinematicJobPlanError::Incompatible(
            "node domain-work classification",
        ));
    }
    let index = u32::try_from(nodes.len())
        .map_err(|_| CinematicJobPlanError::ArithmeticOverflow("node index"))?;
    let mut previous = None;
    for dependency in dependencies {
        if *dependency >= index || previous.is_some_and(|value| value >= *dependency) {
            return Err(CinematicJobPlanError::Incompatible("node dependency order"));
        }
        previous = Some(*dependency);
    }
    let mut preimage = Vec::new();
    preimage.extend_from_slice(&CINEMATIC_JOB_SCHEMA_VERSION.to_le_bytes());
    kind.encode(&mut preimage);
    preimage.extend_from_slice(partition_identity.as_bytes());
    preimage.extend_from_slice(implementation_identity.as_bytes());
    preimage.extend_from_slice(&(local_inputs.len() as u64).to_le_bytes());
    for input in local_inputs {
        preimage.extend_from_slice(input.as_bytes());
    }
    preimage.extend_from_slice(&(dependencies.len() as u64).to_le_bytes());
    for dependency in dependencies {
        preimage.extend_from_slice(nodes[*dependency as usize].identity.as_bytes());
    }
    match (work.render_frame_ordinal, work.render_segment_index) {
        (Some(frame_ordinal), Some(segment_index)) => {
            preimage.push(1);
            preimage.extend_from_slice(&frame_ordinal.to_le_bytes());
            preimage.extend_from_slice(&segment_index.to_le_bytes());
        }
        (None, None) => preimage.push(0),
        _ => {
            return Err(CinematicJobPlanError::Incompatible(
                "partial render-work coordinate",
            ));
        }
    }
    preimage.extend_from_slice(&work.render_tile_start.to_le_bytes());
    preimage.extend_from_slice(&work.render_tiles.to_le_bytes());
    preimage.extend_from_slice(&work.samples_per_tile.to_le_bytes());
    preimage.extend_from_slice(&work.render_paths.to_le_bytes());
    let identity = hash_domain(CINEMATIC_JOB_NODE_IDENTITY_DOMAIN, &preimage);
    let artifact_kind = CinematicArtifactKind::for_job(kind);
    let mut output_preimage = Vec::with_capacity(35);
    output_preimage.extend_from_slice(&CINEMATIC_JOB_SCHEMA_VERSION.to_le_bytes());
    output_preimage.push(artifact_kind as u8);
    output_preimage.extend_from_slice(identity.as_bytes());
    let expected_output_identity =
        hash_domain(CINEMATIC_JOB_OUTPUT_IDENTITY_DOMAIN, &output_preimage);
    nodes.push(CinematicJobNode {
        identity,
        expected_output_identity,
        kind,
        artifact_kind,
        dependencies: dependencies.to_vec(),
        budget,
        work,
    });
    Ok(index)
}

fn validate_node_work_totals(nodes: &[CinematicJobNode]) -> Result<(), CinematicJobPlanError> {
    let mut tiles = 0_u64;
    let mut tile_samples = 0_u64;
    let mut paths = 0_u64;
    for node in nodes {
        tiles = checked_add(tiles, node.work.render_tiles, "total render tiles")?;
        tile_samples = checked_add(
            tile_samples,
            checked_mul(
                node.work.render_tiles,
                node.work.samples_per_tile,
                "node render tile samples",
            )?,
            "total render tile samples",
        )?;
        paths = checked_add(paths, node.work.render_paths, "total render paths")?;
    }
    let _ = (tiles, tile_samples, paths);
    Ok(())
}

fn plan_totals(nodes: &[CinematicJobNode]) -> Result<(u64, u64, u64), CinematicJobPlanError> {
    let mut dependencies = 0_u64;
    let mut work = 0_u64;
    let mut output = 0_u64;
    for node in nodes {
        dependencies = checked_add(
            dependencies,
            node.dependencies.len() as u64,
            "dependency total",
        )?;
        work = checked_add(work, node.budget.work_units, "work total")?;
        output = checked_add(output, node.budget.max_output_bytes, "output total")?;
    }
    Ok((dependencies, work, output))
}

#[allow(clippy::too_many_arguments)]
fn plan_identity(
    configuration_identity: ContentHash,
    render_plan_identity: ContentHash,
    shot_plan_identity: ContentHash,
    nodes: &[CinematicJobNode],
    limits: CinematicJobLimits,
    total_dependencies: u64,
    total_work_units: u64,
    total_output_bytes: u64,
) -> ContentHash {
    let mut hasher = DomainHasher::new(CINEMATIC_JOB_PLAN_IDENTITY_DOMAIN);
    hasher.update(&CINEMATIC_JOB_SCHEMA_VERSION.to_le_bytes());
    hasher.update(configuration_identity.as_bytes());
    hasher.update(render_plan_identity.as_bytes());
    hasher.update(shot_plan_identity.as_bytes());
    hasher.update(&(nodes.len() as u64).to_le_bytes());
    hasher.update(&total_dependencies.to_le_bytes());
    hasher.update(&total_work_units.to_le_bytes());
    hasher.update(&total_output_bytes.to_le_bytes());
    for limit in [
        limits.max_nodes,
        limits.max_dependencies_per_node,
        limits.max_total_dependencies,
        limits.max_total_output_bytes,
        limits.max_snapshot_records,
        limits.max_snapshot_bytes,
        limits.max_events,
    ] {
        hasher.update(&limit.to_le_bytes());
    }
    for node in nodes {
        hasher.update(node.identity.as_bytes());
        hasher.update(node.expected_output_identity.as_bytes());
        hasher.update(&node.budget.work_units.to_le_bytes());
        hasher.update(&node.budget.max_output_bytes.to_le_bytes());
    }
    hasher.finalize()
}

fn hash_pair(domain: &str, first: ContentHash, second: ContentHash) -> ContentHash {
    let mut bytes = [0_u8; 64];
    bytes[..32].copy_from_slice(first.as_bytes());
    bytes[32..].copy_from_slice(second.as_bytes());
    hash_domain(domain, &bytes)
}

fn hash_u64(domain: &str, value: u64) -> ContentHash {
    hash_domain(domain, &value.to_le_bytes())
}

fn checked_add(left: u64, right: u64, what: &'static str) -> Result<u64, CinematicJobPlanError> {
    left.checked_add(right)
        .ok_or(CinematicJobPlanError::ArithmeticOverflow(what))
}

fn checked_mul(left: u64, right: u64, what: &'static str) -> Result<u64, CinematicJobPlanError> {
    left.checked_mul(right)
        .ok_or(CinematicJobPlanError::ArithmeticOverflow(what))
}

fn enforce_limit(
    resource: &'static str,
    observed: u64,
    limit: u64,
) -> Result<(), CinematicJobPlanError> {
    if observed > limit {
        Err(CinematicJobPlanError::Limit {
            resource,
            observed,
            limit,
        })
    } else {
        Ok(())
    }
}

fn node_capacity(node_count: u64) -> Result<usize, CinematicJobPlanError> {
    let index_capacity = checked_add(u64::from(u32::MAX), 1, "node index capacity")?;
    enforce_limit("u32 node-index capacity", node_count, index_capacity)?;
    usize::try_from(node_count)
        .map_err(|_| CinematicJobPlanError::ArithmeticOverflow("node capacity"))
}

fn require_nonzero(name: &'static str, identity: ContentHash) -> Result<(), CinematicJobPlanError> {
    if identity == ZERO_HASH {
        Err(CinematicJobPlanError::ZeroIdentity(name))
    } else {
        Ok(())
    }
}

fn checkpoint_plan(cx: &Cx<'_>) -> Result<(), CinematicJobPlanError> {
    cx.checkpoint()
        .map_err(|_| CinematicJobPlanError::Cancelled)
}

/// Fail-closed graph construction error.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CinematicJobPlanError {
    Cancelled,
    ZeroIdentity(&'static str),
    InvalidLimit(&'static str),
    Incompatible(&'static str),
    ArithmeticOverflow(&'static str),
    Capacity(&'static str),
    Limit {
        resource: &'static str,
        observed: u64,
        limit: u64,
    },
}

impl fmt::Display for CinematicJobPlanError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Cancelled => formatter.write_str("cinematic job planning cancelled"),
            Self::ZeroIdentity(name) => write!(formatter, "zero cinematic identity: {name}"),
            Self::InvalidLimit(name) => write!(formatter, "invalid cinematic limit: {name}"),
            Self::Incompatible(name) => write!(formatter, "incompatible cinematic input: {name}"),
            Self::ArithmeticOverflow(name) => {
                write!(formatter, "cinematic arithmetic overflow: {name}")
            }
            Self::Capacity(name) => write!(formatter, "cinematic allocation refused: {name}"),
            Self::Limit {
                resource,
                observed,
                limit,
            } => write!(
                formatter,
                "cinematic {resource} limit exceeded: observed {observed}, limit {limit}"
            ),
        }
    }
}

impl std::error::Error for CinematicJobPlanError {}

/// Candidate bytes described before checked publication.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CinematicArtifactDescriptor {
    artifact_kind: CinematicArtifactKind,
    output_identity: ContentHash,
    content_identity: ContentHash,
    encoded_bytes_hash: ContentHash,
    encoded_bytes: u64,
}

impl CinematicArtifactDescriptor {
    /// Construct a bounded, nonzero artifact descriptor. The conductor still
    /// checks it against the selected node before invoking the checker.
    pub fn try_new(
        artifact_kind: CinematicArtifactKind,
        output_identity: ContentHash,
        content_identity: ContentHash,
        encoded_bytes_hash: ContentHash,
        encoded_bytes: u64,
    ) -> Result<Self, CinematicArtifactError> {
        for (name, identity) in [
            ("output identity", output_identity),
            ("content identity", content_identity),
            ("encoded bytes hash", encoded_bytes_hash),
        ] {
            if identity == ZERO_HASH {
                return Err(CinematicArtifactError::ZeroIdentity(name));
            }
        }
        if encoded_bytes == 0 {
            return Err(CinematicArtifactError::InvalidLength);
        }
        Ok(Self {
            artifact_kind,
            output_identity,
            content_identity,
            encoded_bytes_hash,
            encoded_bytes,
        })
    }

    #[must_use]
    pub const fn artifact_kind(self) -> CinematicArtifactKind {
        self.artifact_kind
    }

    #[must_use]
    pub const fn output_identity(self) -> ContentHash {
        self.output_identity
    }

    #[must_use]
    pub const fn content_identity(self) -> ContentHash {
        self.content_identity
    }

    #[must_use]
    pub const fn encoded_bytes_hash(self) -> ContentHash {
        self.encoded_bytes_hash
    }

    #[must_use]
    pub const fn encoded_bytes(self) -> u64 {
        self.encoded_bytes
    }
}

/// Checked, atomically discoverable completion record.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CinematicPublishedArtifact {
    node_identity: ContentHash,
    descriptor: CinematicArtifactDescriptor,
}

impl CinematicPublishedArtifact {
    /// Bind a published descriptor to the node whose checker accepted it.
    pub fn try_new(
        node_identity: ContentHash,
        descriptor: CinematicArtifactDescriptor,
    ) -> Result<Self, CinematicArtifactError> {
        if node_identity == ZERO_HASH {
            return Err(CinematicArtifactError::ZeroIdentity("node identity"));
        }
        Ok(Self {
            node_identity,
            descriptor,
        })
    }

    #[must_use]
    pub const fn node_identity(self) -> ContentHash {
        self.node_identity
    }

    #[must_use]
    pub const fn descriptor(self) -> CinematicArtifactDescriptor {
        self.descriptor
    }
}

/// Descriptor/snapshot validation failure.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CinematicArtifactError {
    ZeroIdentity(&'static str),
    InvalidLength,
    DuplicateNode,
    NonCanonicalOrder,
    InvalidSnapshot(&'static str),
    SnapshotLimit {
        resource: &'static str,
        observed: u64,
        limit: u64,
    },
}

impl fmt::Display for CinematicArtifactError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroIdentity(name) => write!(formatter, "zero artifact identity: {name}"),
            Self::InvalidLength => formatter.write_str("artifact byte length must be positive"),
            Self::DuplicateNode => formatter.write_str("duplicate cinematic snapshot node"),
            Self::NonCanonicalOrder => {
                formatter.write_str("cinematic snapshot records are not canonical")
            }
            Self::InvalidSnapshot(reason) => {
                write!(formatter, "invalid cinematic snapshot: {reason}")
            }
            Self::SnapshotLimit {
                resource,
                observed,
                limit,
            } => write!(
                formatter,
                "cinematic snapshot {resource} limit exceeded: observed {observed}, limit {limit}"
            ),
        }
    }
}

impl std::error::Error for CinematicArtifactError {}

const SNAPSHOT_MAGIC: &[u8; 8] = b"FSCJOB01";
const SNAPSHOT_HEADER_BYTES: u64 = 92;
const SNAPSHOT_RECORD_BYTES: u64 = 144;

/// Bounded canonical node-to-artifact index used for crash resume. A snapshot
/// may come from an older plan: reconciliation reuses only exact node IDs and
/// independently rechecks every candidate.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CinematicJobSnapshot {
    identity: ContentHash,
    source_plan_identity: ContentHash,
    records: Vec<CinematicPublishedArtifact>,
    encoded_bytes: u64,
}

impl CinematicJobSnapshot {
    /// Canonicalize a set of published records under caller limits.
    pub fn try_new(
        source_plan_identity: ContentHash,
        mut records: Vec<CinematicPublishedArtifact>,
        limits: CinematicJobLimits,
    ) -> Result<Self, CinematicArtifactError> {
        if source_plan_identity == ZERO_HASH {
            return Err(CinematicArtifactError::ZeroIdentity(
                "snapshot source plan identity",
            ));
        }
        let count = u64::try_from(records.len()).map_err(|_| {
            CinematicArtifactError::InvalidSnapshot("record count conversion overflow")
        })?;
        if count > limits.max_snapshot_records {
            return Err(CinematicArtifactError::SnapshotLimit {
                resource: "record count",
                observed: count,
                limit: limits.max_snapshot_records,
            });
        }
        let encoded_bytes = SNAPSHOT_RECORD_BYTES
            .checked_mul(count)
            .and_then(|records| SNAPSHOT_HEADER_BYTES.checked_add(records))
            .ok_or(CinematicArtifactError::InvalidSnapshot(
                "encoded length overflow",
            ))?;
        if encoded_bytes > limits.max_snapshot_bytes {
            return Err(CinematicArtifactError::SnapshotLimit {
                resource: "encoded bytes",
                observed: encoded_bytes,
                limit: limits.max_snapshot_bytes,
            });
        }
        records.sort_unstable_by_key(|record| record.node_identity);
        if records
            .windows(2)
            .any(|pair| pair[0].node_identity == pair[1].node_identity)
        {
            return Err(CinematicArtifactError::DuplicateNode);
        }
        let identity = snapshot_identity(source_plan_identity, &records);
        Ok(Self {
            identity,
            source_plan_identity,
            records,
            encoded_bytes,
        })
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    #[must_use]
    pub const fn source_plan_identity(&self) -> ContentHash {
        self.source_plan_identity
    }

    #[must_use]
    pub fn records(&self) -> &[CinematicPublishedArtifact] {
        &self.records
    }

    #[must_use]
    pub const fn encoded_bytes(&self) -> u64 {
        self.encoded_bytes
    }

    /// Fixed-width canonical bytes suitable for content-addressed storage.
    pub fn encode_canonical(&self) -> Result<Vec<u8>, CinematicArtifactError> {
        let capacity = usize::try_from(self.encoded_bytes).map_err(|_| {
            CinematicArtifactError::InvalidSnapshot("encoded length exceeds address space")
        })?;
        let mut bytes = Vec::new();
        bytes
            .try_reserve_exact(capacity)
            .map_err(|_| CinematicArtifactError::InvalidSnapshot("snapshot allocation refused"))?;
        bytes.extend_from_slice(SNAPSHOT_MAGIC);
        bytes.extend_from_slice(&CINEMATIC_JOB_SCHEMA_VERSION.to_le_bytes());
        bytes.extend_from_slice(&0_u16.to_le_bytes());
        bytes.extend_from_slice(&self.encoded_bytes.to_le_bytes());
        bytes.extend_from_slice(self.identity.as_bytes());
        bytes.extend_from_slice(self.source_plan_identity.as_bytes());
        bytes.extend_from_slice(&(self.records.len() as u64).to_le_bytes());
        for record in &self.records {
            encode_snapshot_record(&mut bytes, *record);
        }
        if bytes.len() != capacity {
            return Err(CinematicArtifactError::InvalidSnapshot(
                "internal encoded length mismatch",
            ));
        }
        Ok(bytes)
    }

    /// Decode and self-verify a snapshot before any record is considered for
    /// reuse. Owner-specific artifact checks still run during reconciliation.
    pub fn decode_canonical(
        bytes: &[u8],
        limits: CinematicJobLimits,
    ) -> Result<Self, CinematicArtifactError> {
        if bytes.len() as u64 > limits.max_snapshot_bytes {
            return Err(CinematicArtifactError::SnapshotLimit {
                resource: "encoded bytes",
                observed: bytes.len() as u64,
                limit: limits.max_snapshot_bytes,
            });
        }
        let mut reader = SnapshotReader::new(bytes);
        if reader.take(8)? != SNAPSHOT_MAGIC {
            return Err(CinematicArtifactError::InvalidSnapshot("magic"));
        }
        if reader.u16()? != CINEMATIC_JOB_SCHEMA_VERSION || reader.u16()? != 0 {
            return Err(CinematicArtifactError::InvalidSnapshot(
                "schema or reserved field",
            ));
        }
        let declared_bytes = reader.u64()?;
        if declared_bytes != bytes.len() as u64 {
            return Err(CinematicArtifactError::InvalidSnapshot("declared length"));
        }
        let declared_identity = reader.hash()?;
        let source_plan_identity = reader.hash()?;
        let count = reader.u64()?;
        if count > limits.max_snapshot_records {
            return Err(CinematicArtifactError::SnapshotLimit {
                resource: "record count",
                observed: count,
                limit: limits.max_snapshot_records,
            });
        }
        let expected_bytes = SNAPSHOT_RECORD_BYTES
            .checked_mul(count)
            .and_then(|records| SNAPSHOT_HEADER_BYTES.checked_add(records))
            .ok_or(CinematicArtifactError::InvalidSnapshot(
                "encoded length overflow",
            ))?;
        if expected_bytes != declared_bytes {
            return Err(CinematicArtifactError::InvalidSnapshot("record length"));
        }
        let capacity = usize::try_from(count)
            .map_err(|_| CinematicArtifactError::InvalidSnapshot("record capacity"))?;
        let mut records = Vec::new();
        records
            .try_reserve_exact(capacity)
            .map_err(|_| CinematicArtifactError::InvalidSnapshot("record allocation refused"))?;
        for _ in 0..count {
            records.push(decode_snapshot_record(&mut reader)?);
        }
        if !reader.is_finished() {
            return Err(CinematicArtifactError::InvalidSnapshot("trailing bytes"));
        }
        if records
            .windows(2)
            .any(|pair| pair[0].node_identity >= pair[1].node_identity)
        {
            return Err(CinematicArtifactError::NonCanonicalOrder);
        }
        let snapshot = Self::try_new(source_plan_identity, records, limits)?;
        if snapshot.identity != declared_identity || snapshot.encoded_bytes != declared_bytes {
            return Err(CinematicArtifactError::InvalidSnapshot("identity"));
        }
        Ok(snapshot)
    }
}

fn snapshot_identity(
    source_plan_identity: ContentHash,
    records: &[CinematicPublishedArtifact],
) -> ContentHash {
    let mut hasher = DomainHasher::new(CINEMATIC_JOB_SNAPSHOT_IDENTITY_DOMAIN);
    hasher.update(&CINEMATIC_JOB_SCHEMA_VERSION.to_le_bytes());
    hasher.update(source_plan_identity.as_bytes());
    hasher.update(&(records.len() as u64).to_le_bytes());
    for record in records {
        hash_snapshot_record(&mut hasher, *record);
    }
    hasher.finalize()
}

fn hash_snapshot_record(hasher: &mut DomainHasher, record: CinematicPublishedArtifact) {
    let descriptor = record.descriptor;
    hasher.update(record.node_identity.as_bytes());
    hasher.update(descriptor.output_identity.as_bytes());
    hasher.update(&[descriptor.artifact_kind as u8]);
    hasher.update(descriptor.content_identity.as_bytes());
    hasher.update(descriptor.encoded_bytes_hash.as_bytes());
    hasher.update(&descriptor.encoded_bytes.to_le_bytes());
}

fn encode_snapshot_record(output: &mut Vec<u8>, record: CinematicPublishedArtifact) {
    let descriptor = record.descriptor;
    output.extend_from_slice(record.node_identity.as_bytes());
    output.extend_from_slice(descriptor.output_identity.as_bytes());
    output.push(descriptor.artifact_kind as u8);
    output.extend_from_slice(&[0_u8; 7]);
    output.extend_from_slice(descriptor.content_identity.as_bytes());
    output.extend_from_slice(descriptor.encoded_bytes_hash.as_bytes());
    output.extend_from_slice(&descriptor.encoded_bytes.to_le_bytes());
}

fn decode_snapshot_record(
    reader: &mut SnapshotReader<'_>,
) -> Result<CinematicPublishedArtifact, CinematicArtifactError> {
    let node_identity = reader.hash()?;
    let output_identity = reader.hash()?;
    let artifact_kind = CinematicArtifactKind::from_tag(reader.u8()?)
        .ok_or(CinematicArtifactError::InvalidSnapshot("artifact kind"))?;
    if reader.take(7)? != [0_u8; 7] {
        return Err(CinematicArtifactError::InvalidSnapshot(
            "record reserved field",
        ));
    }
    let content_identity = reader.hash()?;
    let encoded_bytes_hash = reader.hash()?;
    let encoded_bytes = reader.u64()?;
    let descriptor = CinematicArtifactDescriptor::try_new(
        artifact_kind,
        output_identity,
        content_identity,
        encoded_bytes_hash,
        encoded_bytes,
    )?;
    CinematicPublishedArtifact::try_new(node_identity, descriptor)
}

struct SnapshotReader<'a> {
    bytes: &'a [u8],
    cursor: usize,
}

impl<'a> SnapshotReader<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, cursor: 0 }
    }

    fn take(&mut self, count: usize) -> Result<&'a [u8], CinematicArtifactError> {
        let end = self
            .cursor
            .checked_add(count)
            .ok_or(CinematicArtifactError::InvalidSnapshot("reader overflow"))?;
        let value = self
            .bytes
            .get(self.cursor..end)
            .ok_or(CinematicArtifactError::InvalidSnapshot("truncated"))?;
        self.cursor = end;
        Ok(value)
    }

    fn u8(&mut self) -> Result<u8, CinematicArtifactError> {
        Ok(self.take(1)?[0])
    }

    fn u16(&mut self) -> Result<u16, CinematicArtifactError> {
        let bytes: [u8; 2] = self
            .take(2)?
            .try_into()
            .map_err(|_| CinematicArtifactError::InvalidSnapshot("u16"))?;
        Ok(u16::from_le_bytes(bytes))
    }

    fn u64(&mut self) -> Result<u64, CinematicArtifactError> {
        let bytes: [u8; 8] = self
            .take(8)?
            .try_into()
            .map_err(|_| CinematicArtifactError::InvalidSnapshot("u64"))?;
        Ok(u64::from_le_bytes(bytes))
    }

    fn hash(&mut self) -> Result<ContentHash, CinematicArtifactError> {
        let identity = ContentHash::from_slice(self.take(32)?)
            .ok_or(CinematicArtifactError::InvalidSnapshot("hash"))?;
        if identity == ZERO_HASH {
            return Err(CinematicArtifactError::InvalidSnapshot("zero hash"));
        }
        Ok(identity)
    }

    fn is_finished(&self) -> bool {
        self.cursor == self.bytes.len()
    }
}

/// Stable backend failure code. Messages remain bounded and deterministic;
/// detailed diagnostics belong in the backend's own structured log.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CinematicNodeFailure {
    code: &'static str,
    retryable: bool,
}

impl CinematicNodeFailure {
    /// Validate an ASCII machine code (1..=64 bytes).
    pub fn try_new(code: &'static str, retryable: bool) -> Result<Self, CinematicArtifactError> {
        if code.is_empty()
            || code.len() > 64
            || !code
                .bytes()
                .all(|byte| byte.is_ascii_lowercase() || byte.is_ascii_digit() || byte == b'_')
        {
            return Err(CinematicArtifactError::InvalidSnapshot(
                "backend failure code",
            ));
        }
        Ok(Self { code, retryable })
    }

    #[must_use]
    pub const fn code(self) -> &'static str {
        self.code
    }

    #[must_use]
    pub const fn retryable(self) -> bool {
        self.retryable
    }
}

/// Result of independently checking a previously published artifact.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CinematicReuseVerdict {
    Valid,
    Invalid,
}

/// Synchronous stage boundary consumed by the deterministic conductor.
///
/// Implementations may use existing scoped renderer/audio crews internally,
/// but every method must return only after those children have drained. In
/// particular, an unwind may not leave detached work behind. `publish` must
/// atomically make the exact returned record discoverable by `discover`.
pub trait CinematicJobBackend {
    type Staged;

    /// Look up a previously checked publication by exact node identity. The
    /// snapshot record is a non-authoritative locator hint; returning it still
    /// requires proving that the atomic publication is currently discoverable.
    fn discover(
        &mut self,
        node: &CinematicJobNode,
        snapshot_hint: Option<CinematicPublishedArtifact>,
        cx: &Cx<'_>,
    ) -> Result<Option<CinematicPublishedArtifact>, CinematicNodeFailure>;

    /// Run the stage-specific decoder/checker over an existing publication.
    fn verify_existing(
        &mut self,
        node: &CinematicJobNode,
        artifact: CinematicPublishedArtifact,
        cx: &Cx<'_>,
    ) -> Result<CinematicReuseVerdict, CinematicNodeFailure>;

    /// Produce an unpublished candidate.
    fn stage(
        &mut self,
        node: &CinematicJobNode,
        cx: &Cx<'_>,
    ) -> Result<Self::Staged, CinematicNodeFailure>;

    /// Describe the exact staged bytes without publishing them.
    fn describe_staged(&self, staged: &Self::Staged) -> CinematicArtifactDescriptor;

    /// Run the owner-specific checker. Success is necessary but not itself a
    /// completion record.
    fn check_staged(
        &mut self,
        node: &CinematicJobNode,
        staged: &Self::Staged,
        cx: &Cx<'_>,
    ) -> Result<(), CinematicNodeFailure>;

    /// Atomically publish a checked candidate and return the discoverable
    /// record. The conductor rejects any changed descriptor.
    fn publish(
        &mut self,
        node: &CinematicJobNode,
        staged: Self::Staged,
        cx: &Cx<'_>,
    ) -> Result<CinematicPublishedArtifact, CinematicNodeFailure>;
}

/// Stable transition phase for actionable event logs and fault injection.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CinematicJobPhase {
    Reconcile,
    VerifyExisting,
    Stage,
    Check,
    Publish,
    Dependency,
}

/// One bounded, deterministic event. Logical counters are useful progress;
/// they are not a wall-clock completion promise.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CinematicJobEventKind {
    Reused,
    Recovered,
    ReuseRejected,
    StageStarted,
    StageFinished,
    CheckPassed,
    Published,
    Failed,
    Panicked,
    Blocked,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CinematicJobEvent {
    ordinal: u64,
    node_identity: ContentHash,
    job_kind: CinematicJobKind,
    phase: CinematicJobPhase,
    kind: CinematicJobEventKind,
    completed_nodes: u64,
    remaining_nodes: u64,
}

impl CinematicJobEvent {
    #[must_use]
    pub const fn ordinal(self) -> u64 {
        self.ordinal
    }

    #[must_use]
    pub const fn node_identity(self) -> ContentHash {
        self.node_identity
    }

    #[must_use]
    pub const fn job_kind(self) -> CinematicJobKind {
        self.job_kind
    }

    #[must_use]
    pub const fn phase(self) -> CinematicJobPhase {
        self.phase
    }

    #[must_use]
    pub const fn kind(self) -> CinematicJobEventKind {
        self.kind
    }

    #[must_use]
    pub const fn completed_nodes(self) -> u64 {
        self.completed_nodes
    }

    #[must_use]
    pub const fn remaining_nodes(self) -> u64 {
        self.remaining_nodes
    }
}

/// Terminal failure attached to its exact node and transition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CinematicJobFailureRecord {
    pub node_identity: ContentHash,
    pub job_kind: CinematicJobKind,
    pub phase: CinematicJobPhase,
    pub failure: CinematicNodeFailure,
    pub panicked: bool,
}

/// Monotone logical progress summary.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct CinematicJobProgress {
    /// Nodes in the admitted plan.
    pub total_nodes: u64,
    /// Reused or newly executed nodes with checked publications.
    pub completed_nodes: u64,
    /// Completed nodes discovered and independently rechecked this pass.
    pub reused_nodes: u64,
    /// Checked publications recovered without an exact snapshot hint. These
    /// nodes remain complete but force every descendant to rebuild.
    pub recovered_nodes: u64,
    /// Completed nodes staged, checked, and published this pass.
    pub executed_nodes: u64,
    /// Nodes whose backend transition failed this pass.
    pub failed_nodes: u64,
    /// Nodes skipped because a dependency was incomplete.
    pub blocked_nodes: u64,
    /// Nodes not complete at the end of this pass.
    pub remaining_nodes: u64,
    /// Declared work units charged for newly executed nodes.
    pub completed_work_units: u64,
    /// Declared work units of nodes that still lack checked publications.
    pub estimated_remaining_work_units: u64,
    /// Immutable render shards in the plan.
    pub total_render_shards: u64,
    /// Render shards with checked publications.
    pub completed_render_shards: u64,
    /// Distinct logical frame-segment tiles after coalescing sample shards.
    pub total_render_tiles: u64,
    /// Distinct tiles whose every sample shard has a checked publication.
    pub completed_render_tiles: u64,
    /// Logical `(tile, sample)` cells covered by all shards.
    pub total_render_tile_samples: u64,
    /// Logical `(tile, sample)` cells covered by completed shards.
    pub completed_render_tile_samples: u64,
    /// Exact renderer path count retained by all shards.
    pub total_render_paths: u64,
    /// Exact renderer path count retained by completed shards.
    pub completed_render_paths: u64,
    /// Event-delimited segments requiring temporal finishing.
    pub total_finished_segments: u64,
    /// Event-delimited segments with finished publications.
    pub completed_finished_segments: u64,
    /// Distinct logical frames represented by finishing nodes.
    pub total_frames: u64,
    /// Frames whose every event-delimited segment is complete.
    pub completed_frames: u64,
    /// Continuous camera shots derived from scene-admitted exposures.
    pub total_shots: u64,
    /// Camera shots whose every event-delimited frame segment is complete.
    pub completed_shots: u64,
    /// Audio stages from controls through the WAV master.
    pub total_audio_stages: u64,
    /// Audio stages with checked publications.
    pub completed_audio_stages: u64,
}

/// Honest terminal state of one conductor pass.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CinematicRunDisposition {
    Complete,
    Failed,
    Cancelled,
    Refused,
}

/// Complete pass report. `Complete` is possible only when every current-plan
/// node has an independently checked publication.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CinematicRunReport {
    disposition: CinematicRunDisposition,
    plan_identity: ContentHash,
    progress: CinematicJobProgress,
    events: Vec<CinematicJobEvent>,
    failures: Vec<CinematicJobFailureRecord>,
    budget_refusal: Option<BudgetRefusal>,
    budget_consumption: Option<fs_exec::BudgetConsumption>,
    snapshot: CinematicJobSnapshot,
}

impl CinematicRunReport {
    #[must_use]
    pub const fn disposition(&self) -> CinematicRunDisposition {
        self.disposition
    }

    #[must_use]
    pub const fn plan_identity(&self) -> ContentHash {
        self.plan_identity
    }

    #[must_use]
    pub const fn progress(&self) -> CinematicJobProgress {
        self.progress
    }

    #[must_use]
    pub fn events(&self) -> &[CinematicJobEvent] {
        &self.events
    }

    #[must_use]
    pub fn failures(&self) -> &[CinematicJobFailureRecord] {
        &self.failures
    }

    #[must_use]
    pub const fn budget_refusal(&self) -> Option<BudgetRefusal> {
        self.budget_refusal
    }

    #[must_use]
    pub const fn budget_consumption(&self) -> Option<fs_exec::BudgetConsumption> {
        self.budget_consumption
    }

    #[must_use]
    pub const fn snapshot(&self) -> &CinematicJobSnapshot {
        &self.snapshot
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum NodeState {
    Pending,
    Reused,
    Recovered,
    Executed,
    Failed,
    Blocked,
}

impl NodeState {
    fn is_complete(self) -> bool {
        matches!(self, Self::Reused | Self::Recovered | Self::Executed)
    }

    fn rebuilt(self) -> bool {
        matches!(self, Self::Recovered | Self::Executed)
    }
}

/// Reconcile and execute one plan in canonical topological order.
///
/// The prior snapshot is only a bounded lookup hint. A backend discovery pass
/// also catches the crash window after atomic publication but before snapshot
/// persistence. Any missing/corrupt/rebuilt dependency forces descendants to
/// rerun, preventing a stale artifact from entering a new final manifest.
pub fn run_cinematic_job_plan<B: CinematicJobBackend>(
    plan: &CinematicJobPlan,
    prior: Option<&CinematicJobSnapshot>,
    backend: &mut B,
    cx: &Cx<'_>,
) -> CinematicRunReport {
    let node_count = plan.nodes.len();
    let mut states = vec![NodeState::Pending; node_count];
    let mut records = BTreeMap::new();
    let event_capacity = node_count.saturating_mul(MAX_EVENTS_PER_NODE as usize);
    let mut events = Vec::with_capacity(event_capacity);
    let mut failures = Vec::new();
    let prior = prior.filter(|snapshot| {
        snapshot.records.len() as u64 <= plan.limits.max_snapshot_records
            && snapshot.encoded_bytes <= plan.limits.max_snapshot_bytes
    });

    let mut admitted = match AdmittedBudget::admit_ambient(cx, plan.total_work_units) {
        Ok(admitted) => admitted,
        Err(refusal) => {
            return finish_report(
                plan,
                CinematicRunDisposition::Refused,
                &states,
                &records,
                events,
                failures,
                Some(refusal),
                None,
            );
        }
    };

    macro_rules! return_if_stopped_after_backend {
        ($phase:literal) => {
            if let Err(refusal) = admitted.observe_deadline($phase, cx) {
                let disposition = disposition_for_refusal(refusal);
                return finish_report(
                    plan,
                    disposition,
                    &states,
                    &records,
                    events,
                    failures,
                    Some(refusal),
                    Some(admitted.consumption()),
                );
            }
        };
    }

    for (index, node) in plan.nodes.iter().enumerate() {
        if let Err(refusal) = admitted.checkpoint("cinematic-node-boundary", cx) {
            let disposition = disposition_for_refusal(refusal);
            return finish_report(
                plan,
                disposition,
                &states,
                &records,
                events,
                failures,
                Some(refusal),
                Some(admitted.consumption()),
            );
        }

        if node
            .dependencies
            .iter()
            .any(|dependency| !states[*dependency as usize].is_complete())
        {
            states[index] = NodeState::Blocked;
            push_event(
                &mut events,
                plan,
                &states,
                node,
                CinematicJobPhase::Dependency,
                CinematicJobEventKind::Blocked,
            );
            continue;
        }
        let dependency_rebuilt = node
            .dependencies
            .iter()
            .any(|dependency| states[*dependency as usize].rebuilt());

        if !dependency_rebuilt {
            let hinted = prior.and_then(|snapshot| {
                snapshot
                    .records
                    .binary_search_by_key(&node.identity, |record| record.node_identity)
                    .ok()
                    .map(|position| snapshot.records[position])
            });
            let discovered = call_backend(node, CinematicJobPhase::Reconcile, || {
                backend.discover(node, hinted, cx)
            });
            match discovered {
                Ok(Some(record)) if record_matches(node, record) => {
                    if let Err(refusal) = admitted.checkpoint("cinematic-reuse-check", cx) {
                        let disposition = disposition_for_refusal(refusal);
                        return finish_report(
                            plan,
                            disposition,
                            &states,
                            &records,
                            events,
                            failures,
                            Some(refusal),
                            Some(admitted.consumption()),
                        );
                    }
                    let verified = call_backend(node, CinematicJobPhase::VerifyExisting, || {
                        backend.verify_existing(node, record, cx)
                    });
                    match verified {
                        Ok(CinematicReuseVerdict::Valid) => {
                            let exact_snapshot_match = hinted == Some(record);
                            states[index] = if exact_snapshot_match {
                                NodeState::Reused
                            } else {
                                NodeState::Recovered
                            };
                            records.insert(node.identity, record);
                            push_event(
                                &mut events,
                                plan,
                                &states,
                                node,
                                CinematicJobPhase::VerifyExisting,
                                if exact_snapshot_match {
                                    CinematicJobEventKind::Reused
                                } else {
                                    CinematicJobEventKind::Recovered
                                },
                            );
                            if let Err(refusal) =
                                admitted.checkpoint("cinematic-after-reuse-check", cx)
                            {
                                let disposition = disposition_for_refusal(refusal);
                                return finish_report(
                                    plan,
                                    disposition,
                                    &states,
                                    &records,
                                    events,
                                    failures,
                                    Some(refusal),
                                    Some(admitted.consumption()),
                                );
                            }
                            continue;
                        }
                        Ok(CinematicReuseVerdict::Invalid) => push_event(
                            &mut events,
                            plan,
                            &states,
                            node,
                            CinematicJobPhase::VerifyExisting,
                            CinematicJobEventKind::ReuseRejected,
                        ),
                        Err(BackendCallError::Failure(failure)) => {
                            return_if_stopped_after_backend!("cinematic-after-verify-error");
                            record_failure(
                                &mut states,
                                index,
                                &mut events,
                                &mut failures,
                                plan,
                                node,
                                CinematicJobPhase::VerifyExisting,
                                failure,
                                false,
                            );
                            continue;
                        }
                        Err(BackendCallError::Panicked) => {
                            return_if_stopped_after_backend!("cinematic-after-verify-panic");
                            record_failure(
                                &mut states,
                                index,
                                &mut events,
                                &mut failures,
                                plan,
                                node,
                                CinematicJobPhase::VerifyExisting,
                                panic_failure(),
                                true,
                            );
                            continue;
                        }
                    }
                }
                Ok(Some(_)) => push_event(
                    &mut events,
                    plan,
                    &states,
                    node,
                    CinematicJobPhase::Reconcile,
                    CinematicJobEventKind::ReuseRejected,
                ),
                Ok(None) => {}
                Err(BackendCallError::Failure(failure)) => {
                    return_if_stopped_after_backend!("cinematic-after-discover-error");
                    record_failure(
                        &mut states,
                        index,
                        &mut events,
                        &mut failures,
                        plan,
                        node,
                        CinematicJobPhase::Reconcile,
                        failure,
                        false,
                    );
                    continue;
                }
                Err(BackendCallError::Panicked) => {
                    return_if_stopped_after_backend!("cinematic-after-discover-panic");
                    record_failure(
                        &mut states,
                        index,
                        &mut events,
                        &mut failures,
                        plan,
                        node,
                        CinematicJobPhase::Reconcile,
                        panic_failure(),
                        true,
                    );
                    continue;
                }
            }
        } else {
            push_event(
                &mut events,
                plan,
                &states,
                node,
                CinematicJobPhase::Dependency,
                CinematicJobEventKind::ReuseRejected,
            );
        }

        if let Err(refusal) = admitted.checkpoint("cinematic-before-stage", cx) {
            let disposition = disposition_for_refusal(refusal);
            return finish_report(
                plan,
                disposition,
                &states,
                &records,
                events,
                failures,
                Some(refusal),
                Some(admitted.consumption()),
            );
        }

        push_event(
            &mut events,
            plan,
            &states,
            node,
            CinematicJobPhase::Stage,
            CinematicJobEventKind::StageStarted,
        );
        if let Err(refusal) = admitted.charge_cost("cinematic-stage", node.budget.work_units) {
            let disposition = disposition_for_refusal(refusal);
            return finish_report(
                plan,
                disposition,
                &states,
                &records,
                events,
                failures,
                Some(refusal),
                Some(admitted.consumption()),
            );
        }
        let staged = call_backend(node, CinematicJobPhase::Stage, || backend.stage(node, cx));
        let staged = match staged {
            Ok(staged) => staged,
            Err(BackendCallError::Failure(failure)) => {
                return_if_stopped_after_backend!("cinematic-after-stage-error");
                record_failure(
                    &mut states,
                    index,
                    &mut events,
                    &mut failures,
                    plan,
                    node,
                    CinematicJobPhase::Stage,
                    failure,
                    false,
                );
                continue;
            }
            Err(BackendCallError::Panicked) => {
                return_if_stopped_after_backend!("cinematic-after-stage-panic");
                record_failure(
                    &mut states,
                    index,
                    &mut events,
                    &mut failures,
                    plan,
                    node,
                    CinematicJobPhase::Stage,
                    panic_failure(),
                    true,
                );
                continue;
            }
        };
        push_event(
            &mut events,
            plan,
            &states,
            node,
            CinematicJobPhase::Stage,
            CinematicJobEventKind::StageFinished,
        );
        if let Err(refusal) = admitted.checkpoint("cinematic-before-check", cx) {
            let disposition = disposition_for_refusal(refusal);
            return finish_report(
                plan,
                disposition,
                &states,
                &records,
                events,
                failures,
                Some(refusal),
                Some(admitted.consumption()),
            );
        }
        let descriptor =
            std::panic::catch_unwind(AssertUnwindSafe(|| backend.describe_staged(&staged)));
        let descriptor = match descriptor {
            Ok(descriptor) => descriptor,
            Err(_) => {
                return_if_stopped_after_backend!("cinematic-after-describe-panic");
                record_failure(
                    &mut states,
                    index,
                    &mut events,
                    &mut failures,
                    plan,
                    node,
                    CinematicJobPhase::Check,
                    panic_failure(),
                    true,
                );
                continue;
            }
        };
        if let Err(refusal) = admitted.checkpoint("cinematic-after-describe", cx) {
            let disposition = disposition_for_refusal(refusal);
            return finish_report(
                plan,
                disposition,
                &states,
                &records,
                events,
                failures,
                Some(refusal),
                Some(admitted.consumption()),
            );
        }
        if !descriptor_matches(node, descriptor) {
            record_failure(
                &mut states,
                index,
                &mut events,
                &mut failures,
                plan,
                node,
                CinematicJobPhase::Check,
                contract_failure(),
                false,
            );
            continue;
        }
        let checked = call_backend(node, CinematicJobPhase::Check, || {
            backend.check_staged(node, &staged, cx)
        });
        match checked {
            Ok(()) => {}
            Err(BackendCallError::Failure(failure)) => {
                return_if_stopped_after_backend!("cinematic-after-check-error");
                record_failure(
                    &mut states,
                    index,
                    &mut events,
                    &mut failures,
                    plan,
                    node,
                    CinematicJobPhase::Check,
                    failure,
                    false,
                );
                continue;
            }
            Err(BackendCallError::Panicked) => {
                return_if_stopped_after_backend!("cinematic-after-check-panic");
                record_failure(
                    &mut states,
                    index,
                    &mut events,
                    &mut failures,
                    plan,
                    node,
                    CinematicJobPhase::Check,
                    panic_failure(),
                    true,
                );
                continue;
            }
        }
        push_event(
            &mut events,
            plan,
            &states,
            node,
            CinematicJobPhase::Check,
            CinematicJobEventKind::CheckPassed,
        );
        if let Err(refusal) = admitted.checkpoint("cinematic-before-publish", cx) {
            let disposition = disposition_for_refusal(refusal);
            return finish_report(
                plan,
                disposition,
                &states,
                &records,
                events,
                failures,
                Some(refusal),
                Some(admitted.consumption()),
            );
        }
        let published = call_backend(node, CinematicJobPhase::Publish, || {
            backend.publish(node, staged, cx)
        });
        let published = match published {
            Ok(published) => published,
            Err(BackendCallError::Failure(failure)) => {
                return_if_stopped_after_backend!("cinematic-after-publish-error");
                record_failure(
                    &mut states,
                    index,
                    &mut events,
                    &mut failures,
                    plan,
                    node,
                    CinematicJobPhase::Publish,
                    failure,
                    false,
                );
                continue;
            }
            Err(BackendCallError::Panicked) => {
                return_if_stopped_after_backend!("cinematic-after-publish-panic");
                record_failure(
                    &mut states,
                    index,
                    &mut events,
                    &mut failures,
                    plan,
                    node,
                    CinematicJobPhase::Publish,
                    panic_failure(),
                    true,
                );
                continue;
            }
        };
        return_if_stopped_after_backend!("cinematic-after-publish-result");
        if published.descriptor != descriptor || !record_matches(node, published) {
            record_failure(
                &mut states,
                index,
                &mut events,
                &mut failures,
                plan,
                node,
                CinematicJobPhase::Publish,
                contract_failure(),
                false,
            );
            continue;
        }
        states[index] = NodeState::Executed;
        records.insert(node.identity, published);
        push_event(
            &mut events,
            plan,
            &states,
            node,
            CinematicJobPhase::Publish,
            CinematicJobEventKind::Published,
        );
        if let Err(refusal) = admitted.checkpoint("cinematic-after-publish", cx) {
            let disposition = disposition_for_refusal(refusal);
            return finish_report(
                plan,
                disposition,
                &states,
                &records,
                events,
                failures,
                Some(refusal),
                Some(admitted.consumption()),
            );
        }
    }

    let disposition = if states.iter().all(|state| state.is_complete()) {
        CinematicRunDisposition::Complete
    } else {
        CinematicRunDisposition::Failed
    };
    finish_report(
        plan,
        disposition,
        &states,
        &records,
        events,
        failures,
        None,
        Some(admitted.consumption()),
    )
}

#[derive(Clone, Copy)]
enum BackendCallError {
    Failure(CinematicNodeFailure),
    Panicked,
}

fn call_backend<T>(
    _node: &CinematicJobNode,
    _phase: CinematicJobPhase,
    call: impl FnOnce() -> Result<T, CinematicNodeFailure>,
) -> Result<T, BackendCallError> {
    match std::panic::catch_unwind(AssertUnwindSafe(call)) {
        Ok(Ok(value)) => Ok(value),
        Ok(Err(failure)) => Err(BackendCallError::Failure(failure)),
        Err(_) => Err(BackendCallError::Panicked),
    }
}

fn descriptor_matches(node: &CinematicJobNode, descriptor: CinematicArtifactDescriptor) -> bool {
    descriptor.artifact_kind == node.artifact_kind
        && descriptor.output_identity == node.expected_output_identity
        && descriptor.content_identity != ZERO_HASH
        && descriptor.encoded_bytes_hash != ZERO_HASH
        && descriptor.encoded_bytes > 0
        && descriptor.encoded_bytes <= node.budget.max_output_bytes
}

fn record_matches(node: &CinematicJobNode, record: CinematicPublishedArtifact) -> bool {
    record.node_identity == node.identity && descriptor_matches(node, record.descriptor)
}

fn panic_failure() -> CinematicNodeFailure {
    CinematicNodeFailure {
        code: "backend_panicked",
        retryable: true,
    }
}

fn contract_failure() -> CinematicNodeFailure {
    CinematicNodeFailure {
        code: "artifact_contract_mismatch",
        retryable: false,
    }
}

#[allow(clippy::too_many_arguments)]
fn record_failure(
    states: &mut [NodeState],
    index: usize,
    events: &mut Vec<CinematicJobEvent>,
    failures: &mut Vec<CinematicJobFailureRecord>,
    plan: &CinematicJobPlan,
    node: &CinematicJobNode,
    phase: CinematicJobPhase,
    failure: CinematicNodeFailure,
    panicked: bool,
) {
    states[index] = NodeState::Failed;
    failures.push(CinematicJobFailureRecord {
        node_identity: node.identity,
        job_kind: node.kind,
        phase,
        failure,
        panicked,
    });
    push_event(
        events,
        plan,
        states,
        node,
        phase,
        if panicked {
            CinematicJobEventKind::Panicked
        } else {
            CinematicJobEventKind::Failed
        },
    );
}

fn push_event(
    events: &mut Vec<CinematicJobEvent>,
    plan: &CinematicJobPlan,
    states: &[NodeState],
    node: &CinematicJobNode,
    phase: CinematicJobPhase,
    kind: CinematicJobEventKind,
) {
    let completed_nodes = states.iter().filter(|state| state.is_complete()).count() as u64;
    events.push(CinematicJobEvent {
        ordinal: events.len() as u64,
        node_identity: node.identity,
        job_kind: node.kind,
        phase,
        kind,
        completed_nodes,
        remaining_nodes: plan.nodes.len() as u64 - completed_nodes,
    });
}

fn disposition_for_refusal(refusal: BudgetRefusal) -> CinematicRunDisposition {
    if matches!(refusal, BudgetRefusal::Cancelled { .. }) {
        CinematicRunDisposition::Cancelled
    } else {
        CinematicRunDisposition::Refused
    }
}

#[allow(clippy::too_many_arguments)]
fn finish_report(
    plan: &CinematicJobPlan,
    disposition: CinematicRunDisposition,
    states: &[NodeState],
    records: &BTreeMap<ContentHash, CinematicPublishedArtifact>,
    events: Vec<CinematicJobEvent>,
    failures: Vec<CinematicJobFailureRecord>,
    budget_refusal: Option<BudgetRefusal>,
    budget_consumption: Option<fs_exec::BudgetConsumption>,
) -> CinematicRunReport {
    let snapshot = CinematicJobSnapshot::try_new(
        plan.identity,
        records.values().copied().collect(),
        plan.limits,
    )
    .expect("current-plan records were preflighted and validated");
    let progress = progress(plan, states);
    CinematicRunReport {
        disposition,
        plan_identity: plan.identity,
        progress,
        events,
        failures,
        budget_refusal,
        budget_consumption,
        snapshot,
    }
}

fn progress(plan: &CinematicJobPlan, states: &[NodeState]) -> CinematicJobProgress {
    let mut value = CinematicJobProgress {
        total_nodes: states.len() as u64,
        ..CinematicJobProgress::default()
    };
    let mut frame_segments = BTreeMap::<u64, (u64, u64)>::new();
    let mut shot_segments = BTreeMap::<u64, (u64, u64)>::new();
    let mut render_tile_ranges = BTreeMap::<(u64, u64, u64, u64), bool>::new();
    for (node, state) in plan.nodes.iter().zip(states) {
        let complete = state.is_complete();
        match state {
            NodeState::Pending => {}
            NodeState::Reused => {
                value.completed_nodes += 1;
                value.reused_nodes += 1;
            }
            NodeState::Recovered => {
                value.completed_nodes += 1;
                value.recovered_nodes += 1;
            }
            NodeState::Executed => {
                value.completed_nodes += 1;
                value.executed_nodes += 1;
                value.completed_work_units = value
                    .completed_work_units
                    .saturating_add(node.budget.work_units);
            }
            NodeState::Failed => value.failed_nodes += 1,
            NodeState::Blocked => value.blocked_nodes += 1,
        }
        match node.kind {
            CinematicJobKind::RenderShard { .. } => {
                value.total_render_shards += 1;
                value.total_render_tile_samples += node.work.tile_samples();
                value.total_render_paths += node.work.render_paths;
                let tile_range = (
                    node.work
                        .render_frame_ordinal
                        .expect("render shards retain a frame ordinal"),
                    node.work
                        .render_segment_index
                        .expect("render shards retain a segment index"),
                    node.work.render_tile_start,
                    node.work.render_tiles,
                );
                render_tile_ranges
                    .entry(tile_range)
                    .and_modify(|all_sample_shards_complete| {
                        *all_sample_shards_complete &= complete;
                    })
                    .or_insert(complete);
                if complete {
                    value.completed_render_shards += 1;
                    value.completed_render_tile_samples += node.work.tile_samples();
                    value.completed_render_paths += node.work.render_paths;
                }
            }
            CinematicJobKind::FinishSegment { frame_ordinal, .. } => {
                value.total_finished_segments += 1;
                value.completed_finished_segments += u64::from(complete);
                let frame = frame_segments.entry(frame_ordinal).or_default();
                frame.0 += 1;
                frame.1 += u64::from(complete);
                let shot_ordinal = node
                    .work
                    .shot_ordinal
                    .expect("finishing nodes retain a validated shot ordinal");
                let shot = shot_segments.entry(shot_ordinal).or_default();
                shot.0 += 1;
                shot.1 += u64::from(complete);
            }
            CinematicJobKind::AudioControls
            | CinematicJobKind::AudioExcitation
            | CinematicJobKind::AudioResampling
            | CinematicJobKind::ModalSynthesis
            | CinematicJobKind::AudioMaster => {
                value.total_audio_stages += 1;
                value.completed_audio_stages += u64::from(complete);
            }
            _ => {}
        }
    }
    value.total_frames = frame_segments.len() as u64;
    value.completed_frames = frame_segments
        .values()
        .filter(|(total, completed)| total == completed)
        .count() as u64;
    value.total_shots = shot_segments.len() as u64;
    value.completed_shots = shot_segments
        .values()
        .filter(|(total, completed)| total == completed)
        .count() as u64;
    for ((_, _, _, tile_count), all_sample_shards_complete) in render_tile_ranges {
        value.total_render_tiles += tile_count;
        if all_sample_shards_complete {
            value.completed_render_tiles += tile_count;
        }
    }
    debug_assert_eq!(value.total_shots, plan.shot_count);
    value.remaining_nodes = value.total_nodes - value.completed_nodes;
    value.estimated_remaining_work_units = plan
        .nodes
        .iter()
        .zip(states)
        .filter(|(_, state)| !state.is_complete())
        .map(|(node, _)| node.budget.work_units)
        .fold(0_u64, u64::saturating_add);
    value
}

#[cfg(test)]
mod tests {
    use std::collections::{BTreeMap, BTreeSet};
    use std::sync::Arc;

    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey, VirtualClock};

    use super::*;

    const STREAM: StreamKey = StreamKey {
        seed: 0x4349_4e45,
        kernel_id: 0x4a4f_4253,
        tile: 0,
        iteration: 0,
    };

    fn identity(label: &str) -> ContentHash {
        hash_domain("org.frankensim.test.cinematic-job.v1", label.as_bytes())
    }

    fn stages() -> CinematicStageIdentities {
        CinematicStageIdentities {
            trajectory: identity("stage-trajectory"),
            render_shard: identity("stage-render-shard"),
            raw_merge: identity("stage-raw-merge"),
            temporal_finish: identity("stage-temporal-finish"),
            image_sequence: identity("stage-image-sequence"),
            audio_controls: identity("stage-audio-controls"),
            audio_excitation: identity("stage-audio-excitation"),
            audio_resampling: identity("stage-audio-resampling"),
            modal_synthesis: identity("stage-modal-synthesis"),
            audio_master: identity("stage-audio-master"),
            bundle_verifier: identity("stage-bundle-verifier"),
            mux_adapter: identity("stage-mux-adapter"),
        }
    }

    fn limits() -> CinematicJobLimits {
        CinematicJobLimits {
            max_nodes: 24,
            max_dependencies_per_node: 64,
            max_total_dependencies: 1_024,
            max_total_output_bytes: 1_000_000,
            max_snapshot_records: 24,
            max_snapshot_bytes: 1_000_000,
            max_events: 1_024,
        }
    }

    fn topology(render_label: &str, trajectory: ContentHash) -> RenderTopology {
        let mut shards = Vec::new();
        let mut segments = Vec::new();
        for index in 0..5_usize {
            shards.push(TopologyShard {
                ordinal: index as u64,
                logical_identity: identity(&format!("{render_label}-shard-{index}")),
                frame_ordinal: index as u64,
                segment_index: 0,
                tile_start: 0,
                tile_count: 1,
                samples_per_tile: 1,
                path_count: 1,
            });
            let mut neighbors = Vec::new();
            let in_second_shot = index >= 2;
            if index > 0 && ((index - 1) >= 2) == in_second_shot {
                neighbors.push(index - 1);
            }
            if index + 1 < 5 && ((index + 1) >= 2) == in_second_shot {
                neighbors.push(index + 1);
            }
            segments.push(TopologySegment {
                frame_ordinal: index as u64,
                frame_position: index as u64,
                segment_index: 0,
                shot_ordinal: u64::from(in_second_shot),
                shot_identity: if index < 2 {
                    identity("shot-a")
                } else {
                    identity("shot-b")
                },
                frame_identity: identity(&format!("{render_label}-frame-{index}")),
                shard_indices: vec![index],
                neighbor_segment_indices: neighbors,
            });
        }
        RenderTopology {
            plan_identity: identity(&format!("{render_label}-plan")),
            sequence_identity: identity("sequence"),
            source_trajectory_identity: trajectory,
            shards,
            segments,
        }
    }

    fn sources(
        trajectory_label: &str,
        image_label: &str,
        audio_label: &str,
        mux_label: &str,
    ) -> PlanSources {
        let trajectory = identity(trajectory_label);
        let image_component = identity(image_label);
        let audio_component = identity(audio_label);
        let mux_component = identity(mux_label);
        let image = hash_pair(
            "org.frankensim.test.cinematic-image-partition.v1",
            trajectory,
            image_component,
        );
        let audio = hash_pair(
            "org.frankensim.test.cinematic-audio-partition.v1",
            trajectory,
            audio_component,
        );
        let mux = hash_pair(
            "org.frankensim.test.cinematic-mux-partition.v1",
            hash_pair(
                "org.frankensim.test.cinematic-av-partition.v1",
                image,
                audio,
            ),
            mux_component,
        );
        let mut configuration = Vec::new();
        for value in [trajectory, image_component, audio_component, mux_component] {
            configuration.extend_from_slice(value.as_bytes());
        }
        let shots = [
            DerivedShotRange {
                identity: identity("shot-a"),
                shot_id: 1,
                first_frame_position: 0,
                frame_count: 2,
            },
            DerivedShotRange {
                identity: identity("shot-b"),
                shot_id: 2,
                first_frame_position: 2,
                frame_count: 3,
            },
        ];
        PlanSources {
            configuration_identity: hash_domain(
                "org.frankensim.test.cinematic-configuration.v1",
                &configuration,
            ),
            trajectory_partition_identity: hash_domain(
                "org.frankensim.test.cinematic-trajectory-partition.v1",
                trajectory.as_bytes(),
            ),
            trajectory_artifact_identity: trajectory,
            image_identity: image,
            audio_identity: audio,
            mux_identity: mux,
            bundle_expectation_identity: identity("bundle-expectation-a"),
            shot_plan_identity: validate_and_identify_shots(&shots, 5, || Ok(()))
                .expect("shot plan"),
            shot_count: shots.len() as u64,
            include_mux: true,
        }
    }

    fn plan(
        trajectory_label: &str,
        image_label: &str,
        audio_label: &str,
        mux_label: &str,
    ) -> CinematicJobPlan {
        plan_with_expectation(
            trajectory_label,
            image_label,
            audio_label,
            mux_label,
            "bundle-expectation-a",
        )
    }

    fn plan_with_expectation(
        trajectory_label: &str,
        image_label: &str,
        audio_label: &str,
        mux_label: &str,
        expectation_label: &str,
    ) -> CinematicJobPlan {
        let mut sources = sources(trajectory_label, image_label, audio_label, mux_label);
        sources.bundle_expectation_identity = identity(expectation_label);
        let render_label = format!("render-{trajectory_label}-{image_label}");
        let topology = topology(&render_label, sources.trajectory_artifact_identity);
        let node_budget = CinematicNodeBudget::try_new(1, 1_024).expect("budget");
        let gate = CancelGate::new_clock_free();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                STREAM,
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            build_plan(
                sources,
                topology,
                stages(),
                CinematicJobBudgets::uniform(node_budget),
                limits(),
                &cx,
            )
            .expect("test plan")
        })
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    enum BackendOp {
        Discover,
        Verify,
        Stage,
        Check,
        Publish,
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct BackendCall {
        operation: BackendOp,
        node: ContentHash,
        kind: CinematicJobKind,
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct CancelHook {
        operation: BackendOp,
        node: ContentHash,
    }

    #[derive(Default)]
    struct MemoryBackend {
        stored: BTreeMap<ContentHash, CinematicPublishedArtifact>,
        calls: Vec<BackendCall>,
        discover_hints: Vec<Option<ContentHash>>,
        corrupt: BTreeSet<ContentHash>,
        fail_stage: BTreeSet<ContentHash>,
        fail_check: BTreeSet<ContentHash>,
        panic_stage: BTreeSet<ContentHash>,
        cancel_hook: Option<CancelHook>,
        cancel_gate: Option<Arc<CancelGate>>,
        deadline_hook: Option<CancelHook>,
        deadline_clock: Option<Arc<VirtualClock>>,
    }

    impl MemoryBackend {
        fn record(&mut self, operation: BackendOp, node: &CinematicJobNode) {
            self.calls.push(BackendCall {
                operation,
                node: node.identity,
                kind: node.kind,
            });
            if self.cancel_hook
                == Some(CancelHook {
                    operation,
                    node: node.identity,
                })
            {
                self.cancel_gate
                    .as_ref()
                    .expect("cancel gate attached")
                    .request();
            }
            if self.deadline_hook
                == Some(CancelHook {
                    operation,
                    node: node.identity,
                })
            {
                self.deadline_clock
                    .as_ref()
                    .expect("deadline clock attached")
                    .advance(2);
            }
        }

        fn stage_calls(&self) -> Vec<CinematicJobKind> {
            self.calls
                .iter()
                .filter(|call| call.operation == BackendOp::Stage)
                .map(|call| call.kind)
                .collect()
        }

        fn staged_set(&self) -> BTreeSet<CinematicJobKind> {
            self.stage_calls().into_iter().collect()
        }

        fn clear_calls(&mut self) {
            self.calls.clear();
            self.discover_hints.clear();
        }
    }

    impl CinematicJobBackend for MemoryBackend {
        type Staged = CinematicArtifactDescriptor;

        fn discover(
            &mut self,
            node: &CinematicJobNode,
            snapshot_hint: Option<CinematicPublishedArtifact>,
            _cx: &Cx<'_>,
        ) -> Result<Option<CinematicPublishedArtifact>, CinematicNodeFailure> {
            self.record(BackendOp::Discover, node);
            self.discover_hints
                .push(snapshot_hint.map(CinematicPublishedArtifact::node_identity));
            Ok(self.stored.get(&node.identity).copied())
        }

        fn verify_existing(
            &mut self,
            node: &CinematicJobNode,
            artifact: CinematicPublishedArtifact,
            _cx: &Cx<'_>,
        ) -> Result<CinematicReuseVerdict, CinematicNodeFailure> {
            self.record(BackendOp::Verify, node);
            Ok(
                if self.corrupt.contains(&node.identity)
                    || self.stored.get(&node.identity) != Some(&artifact)
                {
                    CinematicReuseVerdict::Invalid
                } else {
                    CinematicReuseVerdict::Valid
                },
            )
        }

        fn stage(
            &mut self,
            node: &CinematicJobNode,
            _cx: &Cx<'_>,
        ) -> Result<Self::Staged, CinematicNodeFailure> {
            self.record(BackendOp::Stage, node);
            if self.panic_stage.contains(&node.identity) {
                panic!("injected synchronous backend panic");
            }
            if self.fail_stage.contains(&node.identity) {
                return Err(
                    CinematicNodeFailure::try_new("injected_stage_failure", true)
                        .expect("valid test failure"),
                );
            }
            let content_identity = hash_pair(
                "org.frankensim.test.cinematic-content.v1",
                node.identity,
                node.expected_output_identity,
            );
            let encoded_bytes_hash = hash_pair(
                "org.frankensim.test.cinematic-bytes.v1",
                content_identity,
                node.identity,
            );
            CinematicArtifactDescriptor::try_new(
                node.artifact_kind,
                node.expected_output_identity,
                content_identity,
                encoded_bytes_hash,
                64,
            )
            .map_err(|_| {
                CinematicNodeFailure::try_new("test_descriptor_failure", false)
                    .expect("valid test failure")
            })
        }

        fn describe_staged(&self, staged: &Self::Staged) -> CinematicArtifactDescriptor {
            *staged
        }

        fn check_staged(
            &mut self,
            node: &CinematicJobNode,
            _staged: &Self::Staged,
            _cx: &Cx<'_>,
        ) -> Result<(), CinematicNodeFailure> {
            self.record(BackendOp::Check, node);
            if self.fail_check.contains(&node.identity) {
                Err(
                    CinematicNodeFailure::try_new("injected_check_failure", true)
                        .expect("valid test failure"),
                )
            } else {
                Ok(())
            }
        }

        fn publish(
            &mut self,
            node: &CinematicJobNode,
            staged: Self::Staged,
            _cx: &Cx<'_>,
        ) -> Result<CinematicPublishedArtifact, CinematicNodeFailure> {
            self.record(BackendOp::Publish, node);
            let published =
                CinematicPublishedArtifact::try_new(node.identity, staged).map_err(|_| {
                    CinematicNodeFailure::try_new("test_publish_failure", false)
                        .expect("valid test failure")
                })?;
            self.corrupt.remove(&node.identity);
            self.stored.insert(node.identity, published);
            Ok(published)
        }
    }

    fn run(
        plan: &CinematicJobPlan,
        prior: Option<&CinematicJobSnapshot>,
        backend: &mut MemoryBackend,
        gate: Arc<CancelGate>,
        budget: Budget,
    ) -> CinematicRunReport {
        backend.cancel_gate = Some(Arc::clone(&gate));
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                gate.as_ref(),
                arena,
                STREAM,
                budget,
                ExecMode::Deterministic,
            );
            run_cinematic_job_plan(plan, prior, backend, &cx)
        })
    }

    fn run_with_clock(
        plan: &CinematicJobPlan,
        backend: &mut MemoryBackend,
        gate: Arc<CancelGate>,
        budget: Budget,
        clock: &VirtualClock,
    ) -> CinematicRunReport {
        backend.cancel_gate = Some(Arc::clone(&gate));
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                gate.as_ref(),
                arena,
                STREAM,
                budget,
                ExecMode::Deterministic,
            )
            .with_time_source(clock);
            run_cinematic_job_plan(plan, None, backend, &cx)
        })
    }

    fn node(plan: &CinematicJobPlan, kind: CinematicJobKind) -> &CinematicJobNode {
        plan.nodes
            .iter()
            .find(|node| node.kind == kind)
            .expect("node exists")
    }

    fn image_rerun_set() -> BTreeSet<CinematicJobKind> {
        let mut expected = BTreeSet::new();
        for ordinal in 0..5 {
            expected.insert(CinematicJobKind::RenderShard {
                shard_ordinal: ordinal,
            });
            expected.insert(CinematicJobKind::MergeRawSegment {
                frame_ordinal: ordinal,
                segment_index: 0,
            });
            expected.insert(CinematicJobKind::FinishSegment {
                frame_ordinal: ordinal,
                segment_index: 0,
            });
        }
        expected.extend([
            CinematicJobKind::ImageSequence,
            CinematicJobKind::BundleVerification,
            CinematicJobKind::MuxDerivative,
        ]);
        expected
    }

    fn sound_rerun_set() -> BTreeSet<CinematicJobKind> {
        BTreeSet::from([
            CinematicJobKind::AudioControls,
            CinematicJobKind::AudioExcitation,
            CinematicJobKind::AudioResampling,
            CinematicJobKind::ModalSynthesis,
            CinematicJobKind::AudioMaster,
            CinematicJobKind::BundleVerification,
            CinematicJobKind::MuxDerivative,
        ])
    }

    fn middle_shard_rerun_set() -> BTreeSet<CinematicJobKind> {
        BTreeSet::from([
            CinematicJobKind::RenderShard { shard_ordinal: 2 },
            CinematicJobKind::MergeRawSegment {
                frame_ordinal: 2,
                segment_index: 0,
            },
            CinematicJobKind::FinishSegment {
                frame_ordinal: 2,
                segment_index: 0,
            },
            CinematicJobKind::FinishSegment {
                frame_ordinal: 3,
                segment_index: 0,
            },
            CinematicJobKind::ImageSequence,
            CinematicJobKind::BundleVerification,
            CinematicJobKind::MuxDerivative,
        ])
    }

    #[test]
    fn canonical_graph_exposes_parallel_frontier_and_snapshot_round_trip() {
        let job_plan = plan("trajectory-a", "image-a", "audio-a", "mux-a");
        assert_eq!(job_plan.nodes.len(), 24);
        assert_eq!(
            job_plan,
            plan("trajectory-a", "image-a", "audio-a", "mux-a")
        );

        let mut completed = BTreeSet::new();
        assert_eq!(
            job_plan
                .ready_frontier(&completed)
                .iter()
                .map(|node| node.kind)
                .collect::<Vec<_>>(),
            vec![CinematicJobKind::Trajectory]
        );
        completed.insert(node(&job_plan, CinematicJobKind::Trajectory).identity);
        let frontier = job_plan
            .ready_frontier(&completed)
            .iter()
            .map(|node| node.kind)
            .collect::<Vec<_>>();
        assert_eq!(
            frontier
                .iter()
                .filter(|kind| matches!(kind, CinematicJobKind::RenderShard { .. }))
                .count(),
            5
        );
        assert!(frontier.contains(&CinematicJobKind::AudioControls));

        let mut backend = MemoryBackend::default();
        let report = run(
            &job_plan,
            None,
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        assert_eq!(report.disposition, CinematicRunDisposition::Complete);
        let bytes = report.snapshot.encode_canonical().expect("encode snapshot");
        let decoded = CinematicJobSnapshot::decode_canonical(&bytes, limits()).expect("decode");
        assert_eq!(decoded, report.snapshot);
        let mut corrupt = bytes;
        *corrupt.last_mut().expect("snapshot byte") ^= 1;
        assert!(CinematicJobSnapshot::decode_canonical(&corrupt, limits()).is_err());
    }

    #[test]
    fn fresh_run_checks_before_publish_and_exact_resume_executes_nothing() {
        let plan = plan("trajectory-a", "image-a", "audio-a", "mux-a");
        let mut backend = MemoryBackend::default();
        let first = run(
            &plan,
            None,
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        assert_eq!(first.disposition, CinematicRunDisposition::Complete);
        assert_eq!(first.progress.executed_nodes, first.progress.total_nodes);
        assert_eq!(first.progress.total_render_shards, 5);
        assert_eq!(first.progress.completed_render_shards, 5);
        assert_eq!(first.progress.total_render_tiles, 5);
        assert_eq!(first.progress.completed_render_tiles, 5);
        assert_eq!(first.progress.total_render_tile_samples, 5);
        assert_eq!(first.progress.completed_render_tile_samples, 5);
        assert_eq!(first.progress.total_render_paths, 5);
        assert_eq!(first.progress.completed_render_paths, 5);
        assert_eq!(first.progress.total_finished_segments, 5);
        assert_eq!(first.progress.completed_finished_segments, 5);
        assert_eq!(first.progress.total_frames, 5);
        assert_eq!(first.progress.completed_frames, 5);
        assert_eq!(first.progress.total_shots, 2);
        assert_eq!(first.progress.completed_shots, 2);
        assert_eq!(first.progress.total_audio_stages, 5);
        assert_eq!(first.progress.completed_audio_stages, 5);
        for node in &plan.nodes {
            let transitions = backend
                .calls
                .iter()
                .filter(|call| call.node == node.identity)
                .map(|call| call.operation)
                .collect::<Vec<_>>();
            assert!(
                transitions.windows(3).any(
                    |window| window == [BackendOp::Stage, BackendOp::Check, BackendOp::Publish]
                ),
                "missing stage/check/publish gate for {:?}",
                node.kind
            );
        }

        backend.clear_calls();
        let second = run(
            &plan,
            Some(&first.snapshot),
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        assert_eq!(second.disposition, CinematicRunDisposition::Complete);
        assert_eq!(second.progress.reused_nodes, second.progress.total_nodes);
        assert_eq!(second.progress.recovered_nodes, 0);
        assert!(backend.stage_calls().is_empty());
        assert_eq!(
            backend
                .calls
                .iter()
                .filter(|call| call.operation == BackendOp::Discover)
                .count(),
            plan.nodes().len()
        );
        assert_eq!(
            backend
                .calls
                .iter()
                .filter(|call| call.operation == BackendOp::Verify)
                .count(),
            plan.nodes().len()
        );
        assert_eq!(second.snapshot.identity, first.snapshot.identity);
    }

    #[test]
    fn snapshot_over_current_limits_is_ignored_without_bypassing_live_discovery() {
        let plan = plan("trajectory-a", "image-a", "audio-a", "mux-a");
        let mut backend = MemoryBackend::default();
        let first = run(
            &plan,
            None,
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        let mut records = first.snapshot.records.clone();
        records.push(
            CinematicPublishedArtifact::try_new(
                identity("foreign-snapshot-node"),
                records[0].descriptor,
            )
            .expect("synthetic foreign record"),
        );
        let mut loose_limits = limits();
        loose_limits.max_snapshot_records = 25;
        let oversized = CinematicJobSnapshot::try_new(plan.identity, records, loose_limits)
            .expect("snapshot valid under prior loose limits");

        backend.clear_calls();
        let resumed = run(
            &plan,
            Some(&oversized),
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        assert_eq!(resumed.disposition, CinematicRunDisposition::Complete);
        assert_eq!(resumed.progress.reused_nodes, 0);
        assert_eq!(resumed.progress.recovered_nodes, 1);
        assert_eq!(
            resumed.progress.executed_nodes + 1,
            resumed.progress.total_nodes
        );
        assert_eq!(
            backend.stage_calls().len() as u64,
            resumed.progress.executed_nodes
        );
        assert_eq!(backend.discover_hints.len(), 1);
        assert!(backend.discover_hints.iter().all(Option::is_none));
        assert_eq!(resumed.snapshot.identity, first.snapshot.identity);
    }

    #[test]
    fn changed_recovered_publication_taints_every_descendant() {
        let plan = plan("trajectory-a", "image-a", "audio-a", "mux-a");
        let mut backend = MemoryBackend::default();
        let baseline = run(
            &plan,
            None,
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        let trajectory = node(&plan, CinematicJobKind::Trajectory);
        let prior = backend.stored[&trajectory.identity];
        let replacement_content = identity("replacement-trajectory-publication");
        let replacement_descriptor = CinematicArtifactDescriptor::try_new(
            prior.descriptor.artifact_kind(),
            prior.descriptor.output_identity(),
            replacement_content,
            hash_pair(
                "org.frankensim.test.replacement-cinematic-bytes.v1",
                replacement_content,
                trajectory.identity,
            ),
            prior.descriptor.encoded_bytes(),
        )
        .expect("valid replacement descriptor");
        backend.stored.insert(
            trajectory.identity,
            CinematicPublishedArtifact::try_new(trajectory.identity, replacement_descriptor)
                .expect("valid replacement publication"),
        );

        backend.clear_calls();
        let resumed = run(
            &plan,
            Some(&baseline.snapshot),
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        assert_eq!(resumed.disposition, CinematicRunDisposition::Complete);
        assert_eq!(resumed.progress.recovered_nodes, 1);
        assert_eq!(resumed.progress.reused_nodes, 0);
        assert_eq!(
            resumed.progress.executed_nodes + 1,
            resumed.progress.total_nodes
        );
        assert!(
            !backend
                .stage_calls()
                .contains(&CinematicJobKind::Trajectory)
        );
        assert_eq!(backend.stage_calls().len() + 1, plan.nodes().len());
        assert_ne!(resumed.snapshot.identity, baseline.snapshot.identity);
    }

    #[test]
    fn configuration_partitions_invalidate_only_their_descendants() {
        let baseline = plan("trajectory-a", "image-a", "audio-a", "mux-a");

        for image_label in ["camera-b", "material-b"] {
            let mut backend = MemoryBackend::default();
            let old = run(
                &baseline,
                None,
                &mut backend,
                Arc::new(CancelGate::new_clock_free()),
                Budget::INFINITE,
            );
            backend.clear_calls();
            let changed = plan("trajectory-a", image_label, "audio-a", "mux-a");
            let report = run(
                &changed,
                Some(&old.snapshot),
                &mut backend,
                Arc::new(CancelGate::new_clock_free()),
                Budget::INFINITE,
            );
            assert_eq!(report.disposition, CinematicRunDisposition::Complete);
            assert_eq!(backend.staged_set(), image_rerun_set());
        }

        let mut backend = MemoryBackend::default();
        let old = run(
            &baseline,
            None,
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        backend.clear_calls();
        let sound_changed = plan("trajectory-a", "image-a", "audio-b", "mux-a");
        let sound_report = run(
            &sound_changed,
            Some(&old.snapshot),
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        assert_eq!(sound_report.disposition, CinematicRunDisposition::Complete);
        assert_eq!(backend.staged_set(), sound_rerun_set());

        backend.clear_calls();
        let mux_changed = plan("trajectory-a", "image-a", "audio-b", "mux-b");
        let _ = run(
            &mux_changed,
            Some(&sound_report.snapshot),
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        assert_eq!(
            backend.staged_set(),
            BTreeSet::from([CinematicJobKind::MuxDerivative])
        );

        let mut trajectory_backend = MemoryBackend::default();
        let trajectory_baseline = run(
            &baseline,
            None,
            &mut trajectory_backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        trajectory_backend.clear_calls();
        let trajectory_changed = plan("trajectory-b", "image-a", "audio-a", "mux-a");
        let trajectory_report = run(
            &trajectory_changed,
            Some(&trajectory_baseline.snapshot),
            &mut trajectory_backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        assert_eq!(
            trajectory_report.disposition,
            CinematicRunDisposition::Complete
        );
        assert_eq!(
            trajectory_backend.staged_set(),
            trajectory_changed
                .nodes
                .iter()
                .map(|node| node.kind)
                .collect()
        );

        let mut expectation_backend = MemoryBackend::default();
        let expectation_baseline = run(
            &baseline,
            None,
            &mut expectation_backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        expectation_backend.clear_calls();
        let expectation_changed = plan_with_expectation(
            "trajectory-a",
            "image-a",
            "audio-a",
            "mux-a",
            "bundle-expectation-b",
        );
        let expectation_report = run(
            &expectation_changed,
            Some(&expectation_baseline.snapshot),
            &mut expectation_backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        assert_eq!(
            expectation_report.disposition,
            CinematicRunDisposition::Complete
        );
        assert_eq!(
            expectation_backend.staged_set(),
            BTreeSet::from([
                CinematicJobKind::BundleVerification,
                CinematicJobKind::MuxDerivative,
            ])
        );
    }

    #[test]
    fn corrupt_middle_shard_rebuilds_only_temporal_dependents_and_final_seals() {
        let plan = plan("trajectory-a", "image-a", "audio-a", "mux-a");
        let mut backend = MemoryBackend::default();
        let first = run(
            &plan,
            None,
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        let shard = node(&plan, CinematicJobKind::RenderShard { shard_ordinal: 2 });
        backend.corrupt.insert(shard.identity);
        backend.clear_calls();
        let report = run(
            &plan,
            Some(&first.snapshot),
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        assert_eq!(report.disposition, CinematicRunDisposition::Complete);
        assert_eq!(backend.staged_set(), middle_shard_rerun_set());
    }

    #[test]
    fn missing_publication_rebuilds_only_its_descendants() {
        let plan = plan("trajectory-a", "image-a", "audio-a", "mux-a");
        let mut backend = MemoryBackend::default();
        let first = run(
            &plan,
            None,
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        let shard_identity =
            node(&plan, CinematicJobKind::RenderShard { shard_ordinal: 2 }).identity;
        assert!(backend.stored.remove(&shard_identity).is_some());
        backend.clear_calls();

        let resumed = run(
            &plan,
            Some(&first.snapshot),
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        assert_eq!(resumed.disposition, CinematicRunDisposition::Complete);
        assert_eq!(backend.staged_set(), middle_shard_rerun_set());
    }

    #[test]
    fn failure_and_panic_block_descendants_while_retry_reuses_independent_work() {
        let plan = plan("trajectory-a", "image-a", "audio-a", "mux-a");
        let failed_shard = node(&plan, CinematicJobKind::RenderShard { shard_ordinal: 2 }).identity;
        let mut backend = MemoryBackend::default();
        backend.fail_stage.insert(failed_shard);
        let failed = run(
            &plan,
            None,
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        assert_eq!(failed.disposition, CinematicRunDisposition::Failed);
        assert_eq!(failed.progress.failed_nodes, 1);
        assert!(
            failed
                .snapshot
                .records
                .iter()
                .any(|record| record.node_identity
                    == node(&plan, CinematicJobKind::AudioMaster).identity)
        );
        assert!(
            failed
                .snapshot
                .records
                .iter()
                .all(|record| record.node_identity
                    != node(&plan, CinematicJobKind::BundleVerification).identity)
        );

        backend.fail_stage.clear();
        backend.clear_calls();
        let retried = run(
            &plan,
            Some(&failed.snapshot),
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        assert_eq!(retried.disposition, CinematicRunDisposition::Complete);
        assert!(
            !backend
                .stage_calls()
                .contains(&CinematicJobKind::AudioMaster)
        );
        assert!(
            backend
                .stage_calls()
                .contains(&CinematicJobKind::BundleVerification)
        );

        let panic_node = node(&plan, CinematicJobKind::AudioResampling).identity;
        let mut panic_backend = MemoryBackend::default();
        panic_backend.panic_stage.insert(panic_node);
        let panicked = run(
            &plan,
            None,
            &mut panic_backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        assert_eq!(panicked.disposition, CinematicRunDisposition::Failed);
        assert!(panicked.failures.iter().any(|failure| {
            failure.node_identity == panic_node
                && failure.phase == CinematicJobPhase::Stage
                && failure.panicked
        }));
        assert!(!panic_backend.stored.contains_key(&panic_node));
    }

    #[test]
    fn cancellation_at_each_transaction_boundary_resumes_to_exact_terminal_snapshot() {
        let plan = plan("trajectory-a", "image-a", "audio-a", "mux-a");
        let mut baseline_backend = MemoryBackend::default();
        let baseline = run(
            &plan,
            None,
            &mut baseline_backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        let trajectory = node(&plan, CinematicJobKind::Trajectory).identity;

        let finishing_gate = Arc::new(CancelGate::new_clock_free());
        let finishing_node = node(
            &plan,
            CinematicJobKind::FinishSegment {
                frame_ordinal: 2,
                segment_index: 0,
            },
        )
        .identity;
        let mut finishing_backend = MemoryBackend {
            cancel_hook: Some(CancelHook {
                operation: BackendOp::Stage,
                node: finishing_node,
            }),
            ..MemoryBackend::default()
        };
        let finishing_cancelled = run(
            &plan,
            None,
            &mut finishing_backend,
            Arc::clone(&finishing_gate),
            Budget::INFINITE,
        );
        assert_eq!(
            finishing_cancelled.disposition,
            CinematicRunDisposition::Cancelled
        );
        assert_eq!(finishing_cancelled.progress.total_shots, 2);
        assert_eq!(finishing_cancelled.progress.completed_shots, 1);
        assert_eq!(finishing_cancelled.progress.completed_frames, 2);
        assert_eq!(finishing_cancelled.progress.completed_finished_segments, 2);
        finishing_backend.cancel_hook = None;
        finishing_backend.clear_calls();
        let finishing_resumed = run(
            &plan,
            Some(&finishing_cancelled.snapshot),
            &mut finishing_backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        assert_eq!(
            finishing_resumed.disposition,
            CinematicRunDisposition::Complete
        );
        assert_eq!(
            finishing_resumed.snapshot.identity,
            baseline.snapshot.identity
        );
        assert!(finishing_resumed.progress.reused_nodes > 1);
        assert!(finishing_resumed.progress.executed_nodes > 1);

        for operation in [BackendOp::Stage, BackendOp::Check, BackendOp::Publish] {
            let gate = Arc::new(CancelGate::new_clock_free());
            let mut backend = MemoryBackend {
                cancel_hook: Some(CancelHook {
                    operation,
                    node: trajectory,
                }),
                ..MemoryBackend::default()
            };
            let cancelled = run(
                &plan,
                None,
                &mut backend,
                Arc::clone(&gate),
                Budget::INFINITE,
            );
            assert_eq!(cancelled.disposition, CinematicRunDisposition::Cancelled);
            backend.cancel_hook = None;
            backend.clear_calls();
            let resumed = run(
                &plan,
                Some(&cancelled.snapshot),
                &mut backend,
                Arc::new(CancelGate::new_clock_free()),
                Budget::INFINITE,
            );
            assert_eq!(resumed.disposition, CinematicRunDisposition::Complete);
            assert_eq!(resumed.snapshot.identity, baseline.snapshot.identity);
            backend.clear_calls();
            let replay = run(
                &plan,
                Some(&resumed.snapshot),
                &mut backend,
                Arc::new(CancelGate::new_clock_free()),
                Budget::INFINITE,
            );
            assert_eq!(replay.progress.executed_nodes, 0);
            assert!(backend.stage_calls().is_empty());
        }

        let verify_gate = Arc::new(CancelGate::new_clock_free());
        let mux = node(&plan, CinematicJobKind::MuxDerivative).identity;
        baseline_backend.cancel_hook = Some(CancelHook {
            operation: BackendOp::Verify,
            node: mux,
        });
        baseline_backend.clear_calls();
        let cancelled = run(
            &plan,
            Some(&baseline.snapshot),
            &mut baseline_backend,
            Arc::clone(&verify_gate),
            Budget::INFINITE,
        );
        assert_eq!(cancelled.disposition, CinematicRunDisposition::Cancelled);
        assert_eq!(
            cancelled.progress.reused_nodes,
            cancelled.progress.total_nodes
        );
        baseline_backend.cancel_hook = None;
        baseline_backend.clear_calls();
        let resumed = run(
            &plan,
            Some(&cancelled.snapshot),
            &mut baseline_backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        assert_eq!(resumed.disposition, CinematicRunDisposition::Complete);
        assert_eq!(resumed.snapshot.identity, baseline.snapshot.identity);

        for panic_instead_of_error in [false, true] {
            let gate = Arc::new(CancelGate::new_clock_free());
            let mut backend = MemoryBackend {
                cancel_hook: Some(CancelHook {
                    operation: BackendOp::Stage,
                    node: mux,
                }),
                ..MemoryBackend::default()
            };
            if panic_instead_of_error {
                backend.panic_stage.insert(mux);
            } else {
                backend.fail_stage.insert(mux);
            }
            let cancelled = run(
                &plan,
                None,
                &mut backend,
                Arc::clone(&gate),
                Budget::INFINITE,
            );
            assert_eq!(cancelled.disposition, CinematicRunDisposition::Cancelled);
            assert!(matches!(
                cancelled.budget_refusal,
                Some(BudgetRefusal::Cancelled { .. })
            ));
            assert!(cancelled.failures.is_empty());
            assert_eq!(cancelled.progress.failed_nodes, 0);
            assert!(!backend.stored.contains_key(&mux));

            backend.cancel_hook = None;
            backend.fail_stage.clear();
            backend.panic_stage.clear();
            backend.clear_calls();
            let resumed = run(
                &plan,
                Some(&cancelled.snapshot),
                &mut backend,
                Arc::new(CancelGate::new_clock_free()),
                Budget::INFINITE,
            );
            assert_eq!(resumed.disposition, CinematicRunDisposition::Complete);
            assert_eq!(resumed.snapshot.identity, baseline.snapshot.identity);
        }

        let deadline_clock = Arc::new(VirtualClock::new());
        let mut deadline_backend = MemoryBackend {
            fail_stage: BTreeSet::from([mux]),
            deadline_hook: Some(CancelHook {
                operation: BackendOp::Stage,
                node: mux,
            }),
            deadline_clock: Some(Arc::clone(&deadline_clock)),
            ..MemoryBackend::default()
        };
        let deadline_report = run_with_clock(
            &plan,
            &mut deadline_backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::with_deadline_at_ns(1).with_cost_quota(plan.total_work_units),
            deadline_clock.as_ref(),
        );
        assert_eq!(
            deadline_report.disposition,
            CinematicRunDisposition::Refused
        );
        assert!(matches!(
            deadline_report.budget_refusal,
            Some(BudgetRefusal::DeadlineExpired {
                phase: "cinematic-after-stage-error",
                ..
            })
        ));
        assert!(deadline_report.failures.is_empty());
        assert_eq!(deadline_report.progress.failed_nodes, 0);

        let gate = Arc::new(CancelGate::new_clock_free());
        gate.request();
        let mut backend = MemoryBackend::default();
        let cancelled = run(&plan, None, &mut backend, gate, Budget::INFINITE);
        assert_eq!(cancelled.disposition, CinematicRunDisposition::Cancelled);
        assert!(backend.calls.is_empty());
    }

    #[test]
    fn admission_refuses_before_backend_work_and_progress_is_monotone() {
        let plan = plan("trajectory-a", "image-a", "audio-a", "mux-a");
        let mut backend = MemoryBackend::default();
        let one_short = Budget::new().with_cost_quota(plan.total_work_units - 1);
        let refused = run(
            &plan,
            None,
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            one_short,
        );
        assert_eq!(refused.disposition, CinematicRunDisposition::Refused);
        assert!(matches!(
            refused.budget_refusal,
            Some(BudgetRefusal::CostPlanExceedsQuota { .. })
        ));
        assert!(backend.calls.is_empty());

        let mut backend = MemoryBackend::default();
        let clock = VirtualClock::new();
        let expired = run_with_clock(
            &plan,
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::ZERO,
            &clock,
        );
        assert_eq!(expired.disposition, CinematicRunDisposition::Refused);
        assert!(matches!(
            expired.budget_refusal,
            Some(BudgetRefusal::DeadlineExpiredAtAdmission { .. })
        ));
        assert!(backend.calls.is_empty());

        let exact = Budget::new().with_cost_quota(plan.total_work_units);
        let completed = run(
            &plan,
            None,
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            exact,
        );
        assert_eq!(completed.disposition, CinematicRunDisposition::Complete);
        assert_eq!(completed.progress.remaining_nodes, 0);
        assert_eq!(completed.progress.estimated_remaining_work_units, 0);
        assert_eq!(
            completed.progress.completed_render_shards,
            completed.progress.total_render_shards
        );
        assert_eq!(
            completed.progress.completed_render_tiles,
            completed.progress.total_render_tiles
        );
        assert_eq!(
            completed.progress.completed_render_tile_samples,
            completed.progress.total_render_tile_samples
        );
        assert_eq!(
            completed.progress.completed_render_paths,
            completed.progress.total_render_paths
        );
        assert_eq!(
            completed.progress.completed_finished_segments,
            completed.progress.total_finished_segments
        );
        assert_eq!(
            completed.progress.completed_frames,
            completed.progress.total_frames
        );
        assert_eq!(
            completed.progress.completed_shots,
            completed.progress.total_shots
        );
        assert_eq!(
            completed.progress.completed_audio_stages,
            completed.progress.total_audio_stages
        );
        assert!(completed.events.windows(2).all(|events| {
            events[0].ordinal + 1 == events[1].ordinal
                && events[0].completed_nodes <= events[1].completed_nodes
                && events[0].remaining_nodes >= events[1].remaining_nodes
        }));
    }

    #[test]
    fn checker_failure_never_publishes_or_enters_snapshot() {
        let plan = plan("trajectory-a", "image-a", "audio-a", "mux-a");
        let trajectory = node(&plan, CinematicJobKind::Trajectory).identity;
        let mut backend = MemoryBackend::default();
        backend.fail_check.insert(trajectory);
        let report = run(
            &plan,
            None,
            &mut backend,
            Arc::new(CancelGate::new_clock_free()),
            Budget::INFINITE,
        );
        assert_eq!(report.disposition, CinematicRunDisposition::Failed);
        assert!(!backend.stored.contains_key(&trajectory));
        assert!(
            report
                .snapshot
                .records
                .iter()
                .all(|record| record.node_identity != trajectory)
        );
        let operations = backend
            .calls
            .iter()
            .filter(|call| call.node == trajectory)
            .map(|call| call.operation)
            .collect::<Vec<_>>();
        assert_eq!(
            operations,
            vec![BackendOp::Discover, BackendOp::Stage, BackendOp::Check]
        );
    }
}
