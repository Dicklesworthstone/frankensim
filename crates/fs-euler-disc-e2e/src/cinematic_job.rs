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

use crate::render_sharding::EulerUniformRenderPlan;

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

const ZERO_HASH: ContentHash = ContentHash([0; 32]);

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
}

/// Canonical topological plan for one cinematic composition.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CinematicJobPlan {
    identity: ContentHash,
    configuration_identity: ContentHash,
    render_plan_identity: ContentHash,
    nodes: Vec<CinematicJobNode>,
    total_dependencies: u64,
    total_work_units: u64,
    total_output_bytes: u64,
    limits: CinematicJobLimits,
}

impl CinematicJobPlan {
    /// Construct the graph from already-admitted configuration and render
    /// topology. Component bytes are intentionally not reconstructed from
    /// `.fscine` hashes; stage backends receive the real typed inputs.
    pub fn try_new(
        configuration: &CinematicConfig,
        render_plan: &EulerUniformRenderPlan,
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
        let topology = RenderTopology::from_render_plan(render_plan, cx)?;
        let sources = PlanSources {
            configuration_identity: configuration.composition_identity(),
            trajectory_partition_identity: configuration.trajectory_identity(),
            trajectory_artifact_identity: configuration.input().trajectory.identity(),
            image_identity: configuration.image_identity(),
            audio_identity: configuration.audio_identity(),
            mux_identity: configuration.mux_identity(),
            include_mux: matches!(
                configuration.input().mux_request,
                CinematicMuxRequest::QuarantinedAdapter { .. }
            ),
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
    include_mux: bool,
}

#[derive(Clone)]
struct TopologyShard {
    ordinal: u64,
    logical_identity: ContentHash,
    frame_ordinal: u64,
    segment_index: u64,
}

#[derive(Clone)]
struct TopologySegment {
    frame_ordinal: u64,
    frame_position: u64,
    segment_index: u64,
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

impl RenderTopology {
    fn from_render_plan(
        plan: &EulerUniformRenderPlan,
        cx: &Cx<'_>,
    ) -> Result<Self, CinematicJobPlanError> {
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
            let shard_indices = (first..end).collect();
            let neighbor_segment_indices = plan
                .finishing_neighbors(segment_index)
                .map_err(|_| CinematicJobPlanError::Incompatible("finishing neighborhood"))?
                .map(|neighbor| {
                    plan.segments()
                        .iter()
                        .position(|candidate| core::ptr::eq(candidate, neighbor))
                        .expect("neighbor originates from plan segment slice")
                })
                .collect();
            segments.push(TopologySegment {
                frame_ordinal: segment.frame_ordinal(),
                frame_position: segment.frame_position(),
                segment_index: segment.segment_index(),
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
    if sources.trajectory_artifact_identity != topology.source_trajectory_identity {
        return Err(CinematicJobPlanError::Incompatible(
            "trajectory source identity",
        ));
    }
    validate_topology(&topology)?;

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
    let worst_events = checked_mul(node_count, 8, "event bound")?;
    enforce_limit("events", worst_events, limits.max_events)?;

    let node_capacity = usize::try_from(node_count)
        .map_err(|_| CinematicJobPlanError::ArithmeticOverflow("node capacity"))?;
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
            ],
            &dependencies,
            budgets.temporal_finish,
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
    )?;
    let audio_controls_index = push_node(
        &mut nodes,
        CinematicJobKind::AudioControls,
        sources.audio_identity,
        stages.audio_controls,
        &[topology.sequence_identity],
        &[trajectory_index],
        budgets.audio_controls,
    )?;
    let audio_excitation_index = push_node(
        &mut nodes,
        CinematicJobKind::AudioExcitation,
        sources.audio_identity,
        stages.audio_excitation,
        &[],
        &[audio_controls_index],
        budgets.audio_excitation,
    )?;
    let audio_resampling_index = push_node(
        &mut nodes,
        CinematicJobKind::AudioResampling,
        sources.audio_identity,
        stages.audio_resampling,
        &[],
        &[audio_excitation_index],
        budgets.audio_resampling,
    )?;
    let modal_synthesis_index = push_node(
        &mut nodes,
        CinematicJobKind::ModalSynthesis,
        sources.audio_identity,
        stages.modal_synthesis,
        &[],
        &[audio_resampling_index],
        budgets.modal_synthesis,
    )?;
    let audio_master_index = push_node(
        &mut nodes,
        CinematicJobKind::AudioMaster,
        sources.audio_identity,
        stages.audio_master,
        &[],
        &[modal_synthesis_index],
        budgets.audio_master,
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
        &[topology.plan_identity, topology.sequence_identity],
        &[image_sequence_index, audio_master_index],
        budgets.bundle_verification,
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
        )?;
    }

    debug_assert_eq!(nodes.len(), node_capacity);
    let (total_dependencies, total_work_units, total_output_bytes) = plan_totals(&nodes)?;
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
        nodes,
        total_dependencies,
        total_work_units,
        total_output_bytes,
        limits,
    })
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
        if segment.shard_indices.is_empty()
            || segment
                .shard_indices
                .iter()
                .any(|shard| *shard >= topology.shards.len())
            || segment
                .neighbor_segment_indices
                .iter()
                .any(|neighbor| *neighbor >= topology.segments.len() || *neighbor == index)
        {
            return Err(CinematicJobPlanError::Incompatible(
                "segment dependency topology",
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

fn push_node(
    nodes: &mut Vec<CinematicJobNode>,
    kind: CinematicJobKind,
    partition_identity: ContentHash,
    implementation_identity: ContentHash,
    local_inputs: &[ContentHash],
    dependencies: &[u32],
    budget: CinematicNodeBudget,
) -> Result<u32, CinematicJobPlanError> {
    require_nonzero("node partition identity", partition_identity)?;
    require_nonzero("node implementation identity", implementation_identity)?;
    for input in local_inputs {
        require_nonzero("node local input", *input)?;
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
    });
    Ok(index)
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
        records.sort_unstable_by_key(|record| record.node_identity);
        if records
            .windows(2)
            .any(|pair| pair[0].node_identity == pair[1].node_identity)
        {
            return Err(CinematicArtifactError::DuplicateNode);
        }
        let count = records.len() as u64;
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

    /// Look up a previously checked publication by exact node identity.
    fn discover(
        &mut self,
        node: &CinematicJobNode,
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
    pub total_nodes: u64,
    pub completed_nodes: u64,
    pub reused_nodes: u64,
    pub executed_nodes: u64,
    pub failed_nodes: u64,
    pub blocked_nodes: u64,
    pub remaining_nodes: u64,
    pub completed_work_units: u64,
    pub estimated_remaining_work_units: u64,
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
    Executed,
    Failed,
    Blocked,
}

impl NodeState {
    fn is_complete(self) -> bool {
        matches!(self, Self::Reused | Self::Executed)
    }

    fn rebuilt(self) -> bool {
        matches!(self, Self::Executed)
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
    let mut events = Vec::with_capacity(node_count.saturating_mul(8));
    let mut failures = Vec::new();
    let prior_records = prior
        .map(|snapshot| {
            snapshot
                .records
                .iter()
                .map(|record| (record.node_identity, *record))
                .collect::<BTreeMap<_, _>>()
        })
        .unwrap_or_default();

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
            let hinted = prior_records.get(&node.identity).copied();
            let discovered = if hinted.is_some() {
                Ok(hinted)
            } else {
                call_backend(node, CinematicJobPhase::Reconcile, || {
                    backend.discover(node, cx)
                })
            };
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
                    match call_backend(node, CinematicJobPhase::VerifyExisting, || {
                        backend.verify_existing(node, record, cx)
                    }) {
                        Ok(CinematicReuseVerdict::Valid) => {
                            states[index] = NodeState::Reused;
                            records.insert(node.identity, record);
                            push_event(
                                &mut events,
                                plan,
                                &states,
                                node,
                                CinematicJobPhase::VerifyExisting,
                                CinematicJobEventKind::Reused,
                            );
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
        let staged = match call_backend(node, CinematicJobPhase::Stage, || backend.stage(node, cx))
        {
            Ok(staged) => staged,
            Err(BackendCallError::Failure(failure)) => {
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
            match std::panic::catch_unwind(AssertUnwindSafe(|| backend.describe_staged(&staged))) {
                Ok(descriptor) => descriptor,
                Err(_) => {
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
        match call_backend(node, CinematicJobPhase::Check, || {
            backend.check_staged(node, &staged, cx)
        }) {
            Ok(()) => {}
            Err(BackendCallError::Failure(failure)) => {
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
        let published = match call_backend(node, CinematicJobPhase::Publish, || {
            backend.publish(node, staged, cx)
        }) {
            Ok(published) => published,
            Err(BackendCallError::Failure(failure)) => {
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
    for (node, state) in plan.nodes.iter().zip(states) {
        match state {
            NodeState::Pending => {}
            NodeState::Reused => {
                value.completed_nodes += 1;
                value.reused_nodes += 1;
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
    }
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
