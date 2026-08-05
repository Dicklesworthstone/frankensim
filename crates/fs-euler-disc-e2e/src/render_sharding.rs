//! Bounded L6 planning and immutable-artifact coordination for uniform Euler
//! cinematic renders.
//!
//! The plan is canonical and content-addressed, but it is not a scheduler.
//! Workers receive a pinned plan plus one selected shard, render that shard in
//! isolation, and return an immutable result artifact. A single coordinator
//! later preflights and reads those artifacts from its own [`Ledger`]. The
//! Design Ledger is deliberately not used as a shared multi-process writer,
//! claim table, lease service, or arbitration mechanism.

use core::fmt;
use std::collections::BTreeMap;

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_ledger::{Ledger, LedgerError, PutReceipt};
use fs_render::tracer::{
    DirectStrategy, Film, FilmTimeMode, RenderShardError, RenderShardLimits,
    RenderShardMergeLimits, RenderTileLayout, Sampler, Settings, UniformRenderShardResult,
    UniformRenderShardSpec, merge_uniform_shards, render_cinematic_shard,
};

use crate::render_checkpoint::{
    EulerRenderCheckpointError, euler_render_checkpoint_frame_identity,
};
use crate::render_scene_bridge::{EulerCinematicScene, EulerPreparedFrame, EulerSceneError};

/// Canonical plan wire schema and identity semantics.
pub const EULER_RENDER_SHARD_PLAN_SCHEMA_VERSION: u16 = 1;
/// Domain-separated identity of the complete canonical plan semantics.
pub const EULER_RENDER_SHARD_PLAN_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.uniform-render-shard-plan.v1";
/// Domain-separated identity of one logical plan shard.
pub const EULER_RENDER_LOGICAL_SHARD_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.uniform-render-logical-shard.v1";
/// Dedicated immutable ledger artifact kind for one renderer shard result.
pub const EULER_RENDER_SHARD_RESULT_ARTIFACT_KIND: &str =
    "fs-euler-disc-e2e.uniform-render-shard-result.v1";

const PLAN_MAGIC: &[u8; 8] = b"FSEURS01";
const PLAN_HEADER_BYTES: u64 = 312;
const PLAN_FRAME_BYTES: u64 = 24;
const PLAN_SEGMENT_BYTES: u64 = 72;
const PLAN_SHARD_BYTES: u64 = 120;

/// Explicit resource limits applied before plan allocation, worker execution,
/// artifact publication, or coordinator reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EulerRenderShardLimits {
    max_frames: u64,
    max_shards: u64,
    max_plan_bytes: u64,
    max_paths_per_shard: u64,
    max_result_bytes_per_shard: u64,
    max_aggregate_result_bytes: u64,
}

impl EulerRenderShardLimits {
    /// Admit six positive, caller-visible caps. No implicit unbounded default
    /// exists because decoded plans and foreign artifact sets are hostile
    /// inputs at this boundary.
    pub fn try_new(
        max_frames: u64,
        max_shards: u64,
        max_plan_bytes: u64,
        max_paths_per_shard: u64,
        max_result_bytes_per_shard: u64,
        max_aggregate_result_bytes: u64,
    ) -> Result<Self, EulerRenderShardingError> {
        if max_frames == 0 {
            return Err(EulerRenderShardingError::InvalidLimit("max_frames"));
        }
        if max_shards == 0 {
            return Err(EulerRenderShardingError::InvalidLimit("max_shards"));
        }
        if max_plan_bytes == 0 {
            return Err(EulerRenderShardingError::InvalidLimit("max_plan_bytes"));
        }
        if max_paths_per_shard == 0 {
            return Err(EulerRenderShardingError::InvalidLimit(
                "max_paths_per_shard",
            ));
        }
        if max_result_bytes_per_shard == 0 {
            return Err(EulerRenderShardingError::InvalidLimit(
                "max_result_bytes_per_shard",
            ));
        }
        if max_aggregate_result_bytes == 0 {
            return Err(EulerRenderShardingError::InvalidLimit(
                "max_aggregate_result_bytes",
            ));
        }
        Ok(Self {
            max_frames,
            max_shards,
            max_plan_bytes,
            max_paths_per_shard,
            max_result_bytes_per_shard,
            max_aggregate_result_bytes,
        })
    }

    /// Maximum logical animation frames.
    #[must_use]
    pub const fn max_frames(self) -> u64 {
        self.max_frames
    }

    /// Maximum total render shards.
    #[must_use]
    pub const fn max_shards(self) -> u64 {
        self.max_shards
    }

    /// Maximum canonical plan byte length.
    #[must_use]
    pub const fn max_plan_bytes(self) -> u64 {
        self.max_plan_bytes
    }

    /// Maximum traced paths in one worker shard.
    #[must_use]
    pub const fn max_paths_per_shard(self) -> u64 {
        self.max_paths_per_shard
    }

    /// Maximum encoded bytes returned by one worker shard.
    #[must_use]
    pub const fn max_result_bytes_per_shard(self) -> u64 {
        self.max_result_bytes_per_shard
    }

    /// Maximum sum of stored result lengths admitted by one merge.
    #[must_use]
    pub const fn max_aggregate_result_bytes(self) -> u64 {
        self.max_aggregate_result_bytes
    }
}

/// One caller frame. All event-delimited segments in `prepared` are retained
/// as independent output films.
#[derive(Clone, Copy)]
pub struct EulerRenderFrameInput<'a> {
    frame_ordinal: u64,
    prepared: &'a EulerPreparedFrame,
}

impl<'a> EulerRenderFrameInput<'a> {
    /// Bind a stable animation-frame ordinal to one scene-prepared exposure.
    #[must_use]
    pub const fn new(frame_ordinal: u64, prepared: &'a EulerPreparedFrame) -> Self {
        Self {
            frame_ordinal,
            prepared,
        }
    }

    /// Stable frame ordinal within the caller's nonzero sequence identity.
    #[must_use]
    pub const fn frame_ordinal(self) -> u64 {
        self.frame_ordinal
    }

    /// Scene-bound prepared exposure.
    #[must_use]
    pub const fn prepared(self) -> &'a EulerPreparedFrame {
        self.prepared
    }
}

/// Canonical frame table entry.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EulerRenderPlannedFrame {
    frame_ordinal: u64,
    first_segment: u64,
    segment_count: u64,
}

impl EulerRenderPlannedFrame {
    /// Stable sequence ordinal.
    #[must_use]
    pub const fn frame_ordinal(self) -> u64 {
        self.frame_ordinal
    }

    /// First entry in [`EulerUniformRenderPlan::segments`].
    #[must_use]
    pub const fn first_segment(self) -> u64 {
        self.first_segment
    }

    /// Number of independent event-delimited output films for this frame.
    #[must_use]
    pub const fn segment_count(self) -> u64 {
        self.segment_count
    }
}

/// Canonical event-delimited output entry.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EulerRenderPlannedSegment {
    frame_ordinal: u64,
    frame_position: u64,
    segment_index: u64,
    frame_identity: ContentHash,
    first_shard: u64,
    shard_count: u64,
}

impl EulerRenderPlannedSegment {
    /// Stable sequence ordinal of the parent frame.
    #[must_use]
    pub const fn frame_ordinal(self) -> u64 {
        self.frame_ordinal
    }

    /// Zero-based canonical frame-table position used by finishing-neighbor
    /// radius semantics.
    #[must_use]
    pub const fn frame_position(self) -> u64 {
        self.frame_position
    }

    /// Exact prepared-frame segment index.
    #[must_use]
    pub const fn segment_index(self) -> u64 {
        self.segment_index
    }

    /// Exact event-delimited frame identity from the durable checkpoint
    /// boundary.
    #[must_use]
    pub const fn frame_identity(self) -> ContentHash {
        self.frame_identity
    }

    /// First entry in [`EulerUniformRenderPlan::shards`].
    #[must_use]
    pub const fn first_shard(self) -> u64 {
        self.first_shard
    }

    /// Number of tile/sample rectangles for this independent film.
    #[must_use]
    pub const fn shard_count(self) -> u64 {
        self.shard_count
    }
}

/// One canonical rectangular block of the logical
/// `(frame-segment, tile, sample)` cell space.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EulerRenderPlannedShard {
    shard_ordinal: u64,
    logical_shard_identity: ContentHash,
    frame_ordinal: u64,
    segment_index: u64,
    frame_identity: ContentHash,
    tile_start: u64,
    tile_end: u64,
    sample_start: u32,
    sample_end: u32,
    path_count: u64,
}

impl EulerRenderPlannedShard {
    /// Canonical plan-wide shard ordinal.
    #[must_use]
    pub const fn shard_ordinal(self) -> u64 {
        self.shard_ordinal
    }

    /// Stable coordinator identity derived from the plan identity and exact
    /// logical cell rectangle.
    #[must_use]
    pub const fn logical_shard_identity(self) -> ContentHash {
        self.logical_shard_identity
    }

    /// Parent frame ordinal.
    #[must_use]
    pub const fn frame_ordinal(self) -> u64 {
        self.frame_ordinal
    }

    /// Parent event-delimited segment index.
    #[must_use]
    pub const fn segment_index(self) -> u64 {
        self.segment_index
    }

    /// Exact parent frame-segment identity.
    #[must_use]
    pub const fn frame_identity(self) -> ContentHash {
        self.frame_identity
    }

    /// Inclusive row-major tile start.
    #[must_use]
    pub const fn tile_start(self) -> u64 {
        self.tile_start
    }

    /// Exclusive row-major tile end.
    #[must_use]
    pub const fn tile_end(self) -> u64 {
        self.tile_end
    }

    /// Inclusive absolute sample start.
    #[must_use]
    pub const fn sample_start(self) -> u32 {
        self.sample_start
    }

    /// Exclusive absolute sample end.
    #[must_use]
    pub const fn sample_end(self) -> u32 {
        self.sample_end
    }

    /// Exact traced paths across all pixels covered by this tile block.
    #[must_use]
    pub const fn path_count(self) -> u64 {
        self.path_count
    }
}

/// Small plan-wide audit summary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EulerRenderPlanSummary {
    /// Logical input frames.
    pub frame_count: u64,
    /// Event-delimited independent films.
    pub segment_count: u64,
    /// Rectangular uniform shards.
    pub shard_count: u64,
    /// Exact total paths across all shards.
    pub total_paths: u64,
    /// Canonical encoded plan length.
    pub encoded_plan_bytes: u64,
}

/// Owned, canonical and replay-pinned uniform render plan.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EulerUniformRenderPlan {
    plan_identity: ContentHash,
    sequence_identity: ContentHash,
    source_trajectory_identity: ContentHash,
    source_configuration_identity: ContentHash,
    scene_identity: ContentHash,
    settings: Settings,
    tile_width: u32,
    tile_height: u32,
    tiles_per_shard: u64,
    samples_per_shard: u32,
    finishing_neighbor_radius: u32,
    limits: EulerRenderShardLimits,
    frames: Vec<EulerRenderPlannedFrame>,
    segments: Vec<EulerRenderPlannedSegment>,
    shards: Vec<EulerRenderPlannedShard>,
    total_paths: u64,
    encoded_plan_bytes: u64,
}

impl EulerUniformRenderPlan {
    /// Canonically sort caller frames, reject duplicates, retain every exact
    /// event-delimited segment, and partition each logical cell exactly once.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        scene: &EulerCinematicScene<'_>,
        sequence_identity: ContentHash,
        inputs: &[EulerRenderFrameInput<'_>],
        settings: Settings,
        tile_width: u32,
        tile_height: u32,
        tiles_per_shard: u64,
        samples_per_shard: u32,
        finishing_neighbor_radius: u32,
        limits: EulerRenderShardLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, EulerRenderShardingError> {
        checkpoint(cx)?;
        require_nonzero_identity("sequence_identity", sequence_identity)?;
        require_nonzero_identity(
            "source_trajectory_identity",
            scene.source_trajectory_identity(),
        )?;
        require_nonzero_identity(
            "source_configuration_identity",
            scene.source_configuration_identity(),
        )?;
        require_nonzero_identity("scene_identity", scene.scene_identity())?;
        validate_settings(settings)?;
        if tiles_per_shard == 0 {
            return Err(EulerRenderShardingError::InvalidPartition(
                "tiles_per_shard",
            ));
        }
        if samples_per_shard == 0 {
            return Err(EulerRenderShardingError::InvalidPartition(
                "samples_per_shard",
            ));
        }
        let frame_count = u64::try_from(inputs.len())
            .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("frame_count"))?;
        if frame_count == 0 {
            return Err(EulerRenderShardingError::EmptyPlan);
        }
        if frame_count > limits.max_frames {
            return Err(EulerRenderShardingError::FrameLimit {
                limit: limits.max_frames,
                observed: frame_count,
            });
        }

        let layout =
            RenderTileLayout::try_new(settings.width, settings.height, tile_width, tile_height)
                .map_err(|_| EulerRenderShardingError::InvalidPartition("tile_layout"))?;

        let mut sorted: Vec<EulerRenderFrameInput<'_>> = Vec::new();
        sorted
            .try_reserve_exact(inputs.len())
            .map_err(|_| EulerRenderShardingError::Capacity("frame inputs"))?;
        sorted.extend_from_slice(inputs);
        sorted.sort_by_key(|input| input.frame_ordinal);
        for pair in sorted.windows(2) {
            if pair[0].frame_ordinal == pair[1].frame_ordinal {
                return Err(EulerRenderShardingError::DuplicateFrameOrdinal(
                    pair[0].frame_ordinal,
                ));
            }
        }

        let mut segment_count = 0_u64;
        for input in &sorted {
            checkpoint(cx)?;
            let count = u64::try_from(input.prepared.segments().len())
                .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("segment_count"))?;
            if count == 0 {
                return Err(EulerRenderShardingError::EmptyPreparedFrame(
                    input.frame_ordinal,
                ));
            }
            segment_count = segment_count.checked_add(count).ok_or(
                EulerRenderShardingError::ArithmeticOverflow("segment_count"),
            )?;
        }

        let tile_blocks = div_ceil_u64(layout.tile_count(), tiles_per_shard)?;
        let sample_blocks = div_ceil_u64(u64::from(settings.spp), u64::from(samples_per_shard))?;
        let shards_per_segment = tile_blocks.checked_mul(sample_blocks).ok_or(
            EulerRenderShardingError::ArithmeticOverflow("shards_per_segment"),
        )?;
        let shard_count = segment_count
            .checked_mul(shards_per_segment)
            .ok_or(EulerRenderShardingError::ArithmeticOverflow("shard_count"))?;
        if shard_count > limits.max_shards {
            return Err(EulerRenderShardingError::ShardLimit {
                limit: limits.max_shards,
                observed: shard_count,
            });
        }
        let encoded_plan_bytes = measured_plan_bytes(frame_count, segment_count, shard_count)?;
        if encoded_plan_bytes > limits.max_plan_bytes {
            return Err(EulerRenderShardingError::PlanByteLimit {
                limit: limits.max_plan_bytes,
                observed: encoded_plan_bytes,
            });
        }

        let mut frames = try_vec_capacity(frame_count, "plan frames")?;
        let mut raw_segments = try_vec_capacity(segment_count, "plan segments")?;
        for (frame_position, input) in sorted.iter().enumerate() {
            checkpoint(cx)?;
            let first_segment = u64::try_from(raw_segments.len())
                .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("first_segment"))?;
            for segment_index in 0..input.prepared.segments().len() {
                scene.prepared_segment_shard_binding(input.prepared, segment_index)?;
                let frame_identity =
                    euler_render_checkpoint_frame_identity(input.prepared, segment_index)?;
                raw_segments.push(RawSegment {
                    frame_ordinal: input.frame_ordinal,
                    frame_position: u64::try_from(frame_position).map_err(|_| {
                        EulerRenderShardingError::ArithmeticOverflow("frame_position")
                    })?,
                    segment_index: u64::try_from(segment_index).map_err(|_| {
                        EulerRenderShardingError::ArithmeticOverflow("segment_index")
                    })?,
                    frame_identity,
                });
            }
            frames.push(EulerRenderPlannedFrame {
                frame_ordinal: input.frame_ordinal,
                first_segment,
                segment_count: u64::try_from(input.prepared.segments().len()).map_err(|_| {
                    EulerRenderShardingError::ArithmeticOverflow("frame segment count")
                })?,
            });
        }

        let canonical = CanonicalPlanInputs {
            sequence_identity,
            source_trajectory_identity: scene.source_trajectory_identity(),
            source_configuration_identity: scene.source_configuration_identity(),
            scene_identity: scene.scene_identity(),
            settings,
            tile_width,
            tile_height,
            tiles_per_shard,
            samples_per_shard,
            finishing_neighbor_radius,
            limits,
        };
        let plan_identity = plan_identity(&canonical, &frames, &raw_segments, cx)?;
        let (segments, shards, total_paths) = partition_segments(
            plan_identity,
            &raw_segments,
            layout,
            settings.spp,
            tiles_per_shard,
            samples_per_shard,
            limits,
            cx,
        )?;
        debug_assert_eq!(u64::try_from(shards.len()).ok(), Some(shard_count));
        Ok(Self {
            plan_identity,
            sequence_identity,
            source_trajectory_identity: canonical.source_trajectory_identity,
            source_configuration_identity: canonical.source_configuration_identity,
            scene_identity: canonical.scene_identity,
            settings,
            tile_width,
            tile_height,
            tiles_per_shard,
            samples_per_shard,
            finishing_neighbor_radius,
            limits,
            frames,
            segments,
            shards,
            total_paths,
            encoded_plan_bytes,
        })
    }

    /// Externally pinnable canonical plan identity.
    #[must_use]
    pub const fn plan_identity(&self) -> ContentHash {
        self.plan_identity
    }

    /// Nonzero caller sequence identity.
    #[must_use]
    pub const fn sequence_identity(&self) -> ContentHash {
        self.sequence_identity
    }

    /// Source trajectory artifact identity.
    #[must_use]
    pub const fn source_trajectory_identity(&self) -> ContentHash {
        self.source_trajectory_identity
    }

    /// Scene-builder configuration identity.
    #[must_use]
    pub const fn source_configuration_identity(&self) -> ContentHash {
        self.source_configuration_identity
    }

    /// Complete rendered scene identity.
    #[must_use]
    pub const fn scene_identity(&self) -> ContentHash {
        self.scene_identity
    }

    /// Exact immutable tracer settings.
    #[must_use]
    pub const fn settings(&self) -> Settings {
        self.settings
    }

    /// Exact logical tile layout.
    pub fn tile_layout(&self) -> Result<RenderTileLayout, EulerRenderShardingError> {
        RenderTileLayout::try_new(
            self.settings.width,
            self.settings.height,
            self.tile_width,
            self.tile_height,
        )
        .map_err(|_| EulerRenderShardingError::InvalidPartition("tile_layout"))
    }

    /// Positive number of consecutive logical tiles in a shard block.
    #[must_use]
    pub const fn tiles_per_shard(&self) -> u64 {
        self.tiles_per_shard
    }

    /// Positive number of consecutive uniform samples in a shard block.
    #[must_use]
    pub const fn samples_per_shard(&self) -> u32 {
        self.samples_per_shard
    }

    /// Canonical frame-position radius used for finishing dependencies.
    #[must_use]
    pub const fn finishing_neighbor_radius(&self) -> u32 {
        self.finishing_neighbor_radius
    }

    /// Exact resource caps frozen into the plan identity.
    #[must_use]
    pub const fn limits(&self) -> EulerRenderShardLimits {
        self.limits
    }

    /// Canonical frames sorted by ordinal.
    #[must_use]
    pub fn frames(&self) -> &[EulerRenderPlannedFrame] {
        &self.frames
    }

    /// Canonical event-delimited films sorted by frame then segment.
    #[must_use]
    pub fn segments(&self) -> &[EulerRenderPlannedSegment] {
        &self.segments
    }

    /// Canonical shards sorted by frame, segment, tile block, then sample
    /// block. Every logical cell appears exactly once.
    #[must_use]
    pub fn shards(&self) -> &[EulerRenderPlannedShard] {
        &self.shards
    }

    /// Neighbor segments required by a temporal finishing pass. The iterator
    /// excludes the target frame itself and includes every independently
    /// rendered segment of frames within the configured canonical-position
    /// radius. This is dependency metadata only; no beauty-film blending is
    /// performed by this module.
    pub fn finishing_neighbors(
        &self,
        segment_index: usize,
    ) -> Result<impl Iterator<Item = &EulerRenderPlannedSegment>, EulerRenderShardingError> {
        let target = self
            .segments
            .get(segment_index)
            .ok_or(EulerRenderShardingError::UnknownSegmentIndex(segment_index))?;
        let target_position = target.frame_position;
        let radius = u64::from(self.finishing_neighbor_radius);
        Ok(self.segments.iter().filter(move |candidate| {
            candidate.frame_position != target_position
                && candidate.frame_position.abs_diff(target_position) <= radius
        }))
    }

    /// Compact audit totals.
    #[must_use]
    pub fn summary(&self) -> EulerRenderPlanSummary {
        EulerRenderPlanSummary {
            frame_count: self.frames.len() as u64,
            segment_count: self.segments.len() as u64,
            shard_count: self.shards.len() as u64,
            total_paths: self.total_paths,
            encoded_plan_bytes: self.encoded_plan_bytes,
        }
    }

    /// Encode the complete canonical plan for local file/process exchange.
    /// The caller cap and the plan's own frozen cap are both enforced before
    /// allocation.
    pub fn encode_canonical(
        &self,
        max_bytes: u64,
        cx: &Cx<'_>,
    ) -> Result<Vec<u8>, EulerRenderShardingError> {
        checkpoint(cx)?;
        let limit = max_bytes.min(self.limits.max_plan_bytes);
        if self.encoded_plan_bytes > limit {
            return Err(EulerRenderShardingError::PlanByteLimit {
                limit,
                observed: self.encoded_plan_bytes,
            });
        }
        let capacity = usize::try_from(self.encoded_plan_bytes)
            .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("encoded plan capacity"))?;
        let mut bytes = Vec::new();
        bytes
            .try_reserve_exact(capacity)
            .map_err(|_| EulerRenderShardingError::Capacity("encoded plan"))?;
        bytes.extend_from_slice(PLAN_MAGIC);
        put_u16(&mut bytes, EULER_RENDER_SHARD_PLAN_SCHEMA_VERSION);
        put_u16(&mut bytes, 0);
        put_u64(&mut bytes, self.encoded_plan_bytes);
        put_hash(&mut bytes, self.plan_identity);
        put_hash(&mut bytes, self.sequence_identity);
        put_hash(&mut bytes, self.source_trajectory_identity);
        put_hash(&mut bytes, self.source_configuration_identity);
        put_hash(&mut bytes, self.scene_identity);
        put_settings(&mut bytes, self.settings);
        put_u32(&mut bytes, self.tile_width);
        put_u32(&mut bytes, self.tile_height);
        put_u64(&mut bytes, self.tiles_per_shard);
        put_u32(&mut bytes, self.samples_per_shard);
        put_u32(&mut bytes, self.finishing_neighbor_radius);
        put_limits(&mut bytes, self.limits);
        put_u64(&mut bytes, self.frames.len() as u64);
        put_u64(&mut bytes, self.segments.len() as u64);
        put_u64(&mut bytes, self.shards.len() as u64);
        put_u64(&mut bytes, self.total_paths);
        for frame in &self.frames {
            checkpoint(cx)?;
            put_u64(&mut bytes, frame.frame_ordinal);
            put_u64(&mut bytes, frame.first_segment);
            put_u64(&mut bytes, frame.segment_count);
        }
        for segment in &self.segments {
            checkpoint(cx)?;
            put_u64(&mut bytes, segment.frame_ordinal);
            put_u64(&mut bytes, segment.frame_position);
            put_u64(&mut bytes, segment.segment_index);
            put_hash(&mut bytes, segment.frame_identity);
            put_u64(&mut bytes, segment.first_shard);
            put_u64(&mut bytes, segment.shard_count);
        }
        for shard in &self.shards {
            checkpoint(cx)?;
            put_u64(&mut bytes, shard.shard_ordinal);
            put_hash(&mut bytes, shard.logical_shard_identity);
            put_u64(&mut bytes, shard.frame_ordinal);
            put_u64(&mut bytes, shard.segment_index);
            put_hash(&mut bytes, shard.frame_identity);
            put_u64(&mut bytes, shard.tile_start);
            put_u64(&mut bytes, shard.tile_end);
            put_u32(&mut bytes, shard.sample_start);
            put_u32(&mut bytes, shard.sample_end);
            put_u64(&mut bytes, shard.path_count);
        }
        if bytes.len() != capacity {
            return Err(EulerRenderShardingError::Codec(
                "internal plan length mismatch",
            ));
        }
        Ok(bytes)
    }

    /// Strictly decode canonical plan bytes under both an external expected
    /// identity pin and a caller byte cap. Truncation, trailing bytes,
    /// reserved-bit changes, noncanonical order, and recomputation mismatch
    /// all refuse before a plan is returned.
    pub fn decode_canonical(
        bytes: &[u8],
        max_bytes: u64,
        expected_plan_identity: ContentHash,
        cx: &Cx<'_>,
    ) -> Result<Self, EulerRenderShardingError> {
        checkpoint(cx)?;
        require_nonzero_identity("expected_plan_identity", expected_plan_identity)?;
        let observed_len = u64::try_from(bytes.len())
            .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("plan input length"))?;
        if observed_len > max_bytes {
            return Err(EulerRenderShardingError::PlanByteLimit {
                limit: max_bytes,
                observed: observed_len,
            });
        }
        let mut reader = PlanReader::new(bytes);
        if reader.take(PLAN_MAGIC.len())? != PLAN_MAGIC {
            return Err(EulerRenderShardingError::Codec("invalid plan magic"));
        }
        if reader.u16()? != EULER_RENDER_SHARD_PLAN_SCHEMA_VERSION {
            return Err(EulerRenderShardingError::Codec(
                "unsupported plan schema version",
            ));
        }
        if reader.u16()? != 0 {
            return Err(EulerRenderShardingError::Codec(
                "nonzero reserved plan header",
            ));
        }
        let declared_len = reader.u64()?;
        if declared_len != observed_len {
            return Err(EulerRenderShardingError::Codec(
                "declared plan length mismatch",
            ));
        }
        let declared_plan_identity = reader.hash()?;
        if declared_plan_identity != expected_plan_identity {
            return Err(EulerRenderShardingError::PlanIdentityMismatch {
                expected: expected_plan_identity,
                actual: declared_plan_identity,
            });
        }
        let canonical = CanonicalPlanInputs {
            sequence_identity: reader.hash()?,
            source_trajectory_identity: reader.hash()?,
            source_configuration_identity: reader.hash()?,
            scene_identity: reader.hash()?,
            settings: reader.settings()?,
            tile_width: reader.u32()?,
            tile_height: reader.u32()?,
            tiles_per_shard: reader.u64()?,
            samples_per_shard: reader.u32()?,
            finishing_neighbor_radius: reader.u32()?,
            limits: reader.limits()?,
        };
        require_nonzero_identity("sequence_identity", canonical.sequence_identity)?;
        require_nonzero_identity(
            "source_trajectory_identity",
            canonical.source_trajectory_identity,
        )?;
        require_nonzero_identity(
            "source_configuration_identity",
            canonical.source_configuration_identity,
        )?;
        require_nonzero_identity("scene_identity", canonical.scene_identity)?;
        validate_settings(canonical.settings)?;
        if canonical.tiles_per_shard == 0 {
            return Err(EulerRenderShardingError::InvalidPartition(
                "tiles_per_shard",
            ));
        }
        if canonical.samples_per_shard == 0 {
            return Err(EulerRenderShardingError::InvalidPartition(
                "samples_per_shard",
            ));
        }
        let layout = RenderTileLayout::try_new(
            canonical.settings.width,
            canonical.settings.height,
            canonical.tile_width,
            canonical.tile_height,
        )
        .map_err(|_| EulerRenderShardingError::InvalidPartition("tile_layout"))?;
        let frame_count = reader.u64()?;
        let segment_count = reader.u64()?;
        let shard_count = reader.u64()?;
        let declared_total_paths = reader.u64()?;
        if frame_count == 0 {
            return Err(EulerRenderShardingError::EmptyPlan);
        }
        if frame_count > canonical.limits.max_frames {
            return Err(EulerRenderShardingError::FrameLimit {
                limit: canonical.limits.max_frames,
                observed: frame_count,
            });
        }
        if shard_count > canonical.limits.max_shards {
            return Err(EulerRenderShardingError::ShardLimit {
                limit: canonical.limits.max_shards,
                observed: shard_count,
            });
        }
        let expected_tile_blocks = div_ceil_u64(layout.tile_count(), canonical.tiles_per_shard)?;
        let expected_sample_blocks = div_ceil_u64(
            u64::from(canonical.settings.spp),
            u64::from(canonical.samples_per_shard),
        )?;
        let expected_shards = segment_count
            .checked_mul(expected_tile_blocks)
            .and_then(|count| count.checked_mul(expected_sample_blocks))
            .ok_or(EulerRenderShardingError::ArithmeticOverflow(
                "decoded shard count",
            ))?;
        if shard_count != expected_shards {
            return Err(EulerRenderShardingError::Codec(
                "declared shard count does not match partition",
            ));
        }
        let measured = measured_plan_bytes(frame_count, segment_count, shard_count)?;
        if measured != observed_len {
            return Err(EulerRenderShardingError::Codec(
                "noncanonical plan table lengths",
            ));
        }
        if measured > canonical.limits.max_plan_bytes {
            return Err(EulerRenderShardingError::PlanByteLimit {
                limit: canonical.limits.max_plan_bytes,
                observed: measured,
            });
        }

        let mut frames = try_vec_capacity(frame_count, "decoded frames")?;
        for _ in 0..frame_count {
            checkpoint(cx)?;
            frames.push(EulerRenderPlannedFrame {
                frame_ordinal: reader.u64()?,
                first_segment: reader.u64()?,
                segment_count: reader.u64()?,
            });
        }
        let mut raw_segments = try_vec_capacity(segment_count, "decoded segments")?;
        let mut declared_segments = try_vec_capacity(segment_count, "decoded segment table")?;
        for _ in 0..segment_count {
            checkpoint(cx)?;
            let raw = RawSegment {
                frame_ordinal: reader.u64()?,
                frame_position: reader.u64()?,
                segment_index: reader.u64()?,
                frame_identity: reader.hash()?,
            };
            let planned = EulerRenderPlannedSegment {
                frame_ordinal: raw.frame_ordinal,
                frame_position: raw.frame_position,
                segment_index: raw.segment_index,
                frame_identity: raw.frame_identity,
                first_shard: reader.u64()?,
                shard_count: reader.u64()?,
            };
            raw_segments.push(raw);
            declared_segments.push(planned);
        }
        let mut declared_shards = try_vec_capacity(shard_count, "decoded shards")?;
        for _ in 0..shard_count {
            checkpoint(cx)?;
            declared_shards.push(EulerRenderPlannedShard {
                shard_ordinal: reader.u64()?,
                logical_shard_identity: reader.hash()?,
                frame_ordinal: reader.u64()?,
                segment_index: reader.u64()?,
                frame_identity: reader.hash()?,
                tile_start: reader.u64()?,
                tile_end: reader.u64()?,
                sample_start: reader.u32()?,
                sample_end: reader.u32()?,
                path_count: reader.u64()?,
            });
        }
        if !reader.is_finished() {
            return Err(EulerRenderShardingError::Codec("trailing plan bytes"));
        }
        validate_canonical_tables(&frames, &raw_segments)?;
        let actual_plan_identity = plan_identity(&canonical, &frames, &raw_segments, cx)?;
        if actual_plan_identity != declared_plan_identity {
            return Err(EulerRenderShardingError::PlanIdentityMismatch {
                expected: declared_plan_identity,
                actual: actual_plan_identity,
            });
        }
        let (segments, shards, total_paths) = partition_segments(
            actual_plan_identity,
            &raw_segments,
            layout,
            canonical.settings.spp,
            canonical.tiles_per_shard,
            canonical.samples_per_shard,
            canonical.limits,
            cx,
        )?;
        if declared_segments != segments || declared_shards != shards {
            return Err(EulerRenderShardingError::Codec(
                "noncanonical segment or shard table",
            ));
        }
        if total_paths != declared_total_paths {
            return Err(EulerRenderShardingError::Codec(
                "declared total path count mismatch",
            ));
        }
        let plan = Self {
            plan_identity: actual_plan_identity,
            sequence_identity: canonical.sequence_identity,
            source_trajectory_identity: canonical.source_trajectory_identity,
            source_configuration_identity: canonical.source_configuration_identity,
            scene_identity: canonical.scene_identity,
            settings: canonical.settings,
            tile_width: canonical.tile_width,
            tile_height: canonical.tile_height,
            tiles_per_shard: canonical.tiles_per_shard,
            samples_per_shard: canonical.samples_per_shard,
            finishing_neighbor_radius: canonical.finishing_neighbor_radius,
            limits: canonical.limits,
            frames,
            segments,
            shards,
            total_paths,
            encoded_plan_bytes: measured,
        };
        let canonical_bytes = plan.encode_canonical(max_bytes, cx)?;
        if canonical_bytes != bytes {
            return Err(EulerRenderShardingError::Codec(
                "noncanonical plan encoding",
            ));
        }
        Ok(plan)
    }
}

/// Immutable pointer returned after the single coordinator publishes one
/// canonical renderer result to its local ledger.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct EulerRenderShardArtifactRef {
    logical_shard_identity: ContentHash,
    artifact_hash: ContentHash,
}

impl EulerRenderShardArtifactRef {
    /// Construct a transport reference. Trust is established only when merge
    /// preflights the ledger envelope and strictly decodes the bound result.
    #[must_use]
    pub const fn new(logical_shard_identity: ContentHash, artifact_hash: ContentHash) -> Self {
        Self {
            logical_shard_identity,
            artifact_hash,
        }
    }

    /// Logical coordinator shard identity.
    #[must_use]
    pub const fn logical_shard_identity(self) -> ContentHash {
        self.logical_shard_identity
    }

    /// Ledger content hash of the canonical result bytes.
    #[must_use]
    pub const fn artifact_hash(self) -> ContentHash {
        self.artifact_hash
    }
}

/// Coordinator publication receipt. `deduped` is true for an exact canonical
/// duplicate and is therefore an idempotent success, not a second result.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EulerRenderShardArtifactReceipt {
    /// Reusable immutable reference.
    pub artifact: EulerRenderShardArtifactRef,
    /// Renderer-internal exact spec identity.
    pub renderer_shard_identity: ContentHash,
    /// Canonical renderer result identity.
    pub result_identity: ContentHash,
    /// Stored canonical byte length.
    pub len: u64,
    /// Whether identical bytes and envelope were already present.
    pub deduped: bool,
}

/// Execute exactly one selected uniform shard after recomputing the complete
/// scene/frame/segment binding. There is intentionally no adaptive or
/// work-stealing sample repartition at this boundary.
pub fn execute_uniform_render_shard(
    plan: &EulerUniformRenderPlan,
    scene: &EulerCinematicScene<'_>,
    inputs: &[EulerRenderFrameInput<'_>],
    shard_index: usize,
    cx: &Cx<'_>,
) -> Result<UniformRenderShardResult, EulerRenderShardingError> {
    let indexed_inputs = index_frame_inputs(plan, scene, inputs, cx)?;
    let bound = bind_shard_spec_from_index(plan, scene, &indexed_inputs, shard_index, cx)?;
    Ok(render_cinematic_shard(
        scene.scene(),
        scene.camera(),
        bound.cut_side,
        cx,
        &bound.spec,
        bound.shutter,
    )?)
}

/// Publish a completed shard result into the coordinator's immutable local
/// artifact store. All identities and the complete canonical payload are
/// validated before `Ledger::put_artifact`; workers must return bytes to the
/// coordinator rather than concurrently opening this ledger for writes.
pub fn store_uniform_render_shard_artifact(
    ledger: &Ledger,
    plan: &EulerUniformRenderPlan,
    scene: &EulerCinematicScene<'_>,
    inputs: &[EulerRenderFrameInput<'_>],
    shard_index: usize,
    result: &UniformRenderShardResult,
    cx: &Cx<'_>,
) -> Result<EulerRenderShardArtifactReceipt, EulerRenderShardingError> {
    checkpoint(cx)?;
    let indexed_inputs = index_frame_inputs(plan, scene, inputs, cx)?;
    let shard = *plan
        .shards
        .get(shard_index)
        .ok_or(EulerRenderShardingError::UnknownShardIndex(shard_index))?;
    let bound = bind_shard_spec_from_index(plan, scene, &indexed_inputs, shard_index, cx)?;
    store_bound_shard_result(ledger, plan, shard, &bound, result, cx)
}

fn store_bound_shard_result(
    ledger: &Ledger,
    plan: &EulerUniformRenderPlan,
    shard: EulerRenderPlannedShard,
    bound: &BoundShardSpec,
    result: &UniformRenderShardResult,
    cx: &Cx<'_>,
) -> Result<EulerRenderShardArtifactReceipt, EulerRenderShardingError> {
    if result.spec() != &bound.spec
        || result.shard_identity() != bound.spec.shard_identity()
        || result.spec().plan_identity() != plan.plan_identity
    {
        return Err(EulerRenderShardingError::ShardResultBindingMismatch);
    }
    let bytes = result.encode_canonical(plan.limits.max_result_bytes_per_shard, cx)?;
    let len = u64::try_from(bytes.len())
        .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("result byte length"))?;
    if len > plan.limits.max_result_bytes_per_shard {
        return Err(EulerRenderShardingError::ResultByteLimit {
            limit: plan.limits.max_result_bytes_per_shard,
            observed: len,
        });
    }
    if len != result.encoded_result_bytes() {
        return Err(EulerRenderShardingError::ShardResultBindingMismatch);
    }
    let meta = format!(
        "{{\"logical_shard_identity\":\"{}\",\"plan_identity\":\"{}\",\"renderer_shard_identity\":\"{}\",\"result_identity\":\"{}\"}}",
        shard.logical_shard_identity,
        plan.plan_identity,
        result.shard_identity(),
        result.result_identity(),
    );
    checkpoint(cx)?;
    let PutReceipt {
        hash,
        len: stored_len,
        deduped,
        ..
    } = ledger.put_artifact(EULER_RENDER_SHARD_RESULT_ARTIFACT_KIND, &bytes, Some(&meta))?;
    debug_assert_eq!(stored_len, len, "ledger receipt length must match input");
    Ok(EulerRenderShardArtifactReceipt {
        artifact: EulerRenderShardArtifactRef::new(shard.logical_shard_identity, hash),
        renderer_shard_identity: result.shard_identity(),
        result_identity: result.result_identity(),
        len,
        deduped,
    })
}

/// Strict cross-process ingestion seam. A worker returns only canonical bytes;
/// the coordinator reconstructs the externally pinned spec from the original
/// scene/prepared-frame binding, decodes privately, and only then delegates to
/// [`store_uniform_render_shard_artifact`]. No worker opens the coordinator's
/// ledger for writing.
pub fn store_uniform_render_shard_artifact_bytes(
    ledger: &Ledger,
    plan: &EulerUniformRenderPlan,
    scene: &EulerCinematicScene<'_>,
    inputs: &[EulerRenderFrameInput<'_>],
    shard_index: usize,
    bytes: &[u8],
    cx: &Cx<'_>,
) -> Result<EulerRenderShardArtifactReceipt, EulerRenderShardingError> {
    checkpoint(cx)?;
    let observed = u64::try_from(bytes.len())
        .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("worker result length"))?;
    if observed > plan.limits.max_result_bytes_per_shard {
        return Err(EulerRenderShardingError::ResultByteLimit {
            limit: plan.limits.max_result_bytes_per_shard,
            observed,
        });
    }
    let indexed_inputs = index_frame_inputs(plan, scene, inputs, cx)?;
    let bound = bind_shard_spec_from_index(plan, scene, &indexed_inputs, shard_index, cx)?;
    let result = UniformRenderShardResult::decode_canonical(
        bytes,
        plan.limits.max_result_bytes_per_shard,
        &bound.spec,
        plan.plan_identity,
        bound.spec.shard_identity(),
        cx,
    )?;
    let shard = *plan
        .shards
        .get(shard_index)
        .ok_or(EulerRenderShardingError::UnknownShardIndex(shard_index))?;
    store_bound_shard_result(ledger, plan, shard, &bound, &result, cx)
}

/// Load and merge exactly one event-delimited film. All unique artifact
/// envelope lengths and kinds are preflighted together under the aggregate cap
/// before any payload byte is read. Exact duplicate references collapse;
/// differing valid results for one logical shard conflict; missing, foreign,
/// corrupt, or mismatched inputs return no film.
#[allow(clippy::too_many_arguments)]
pub fn merge_uniform_render_segment_artifacts(
    ledger: &Ledger,
    plan: &EulerUniformRenderPlan,
    scene: &EulerCinematicScene<'_>,
    inputs: &[EulerRenderFrameInput<'_>],
    segment_plan_index: usize,
    artifacts: &[EulerRenderShardArtifactRef],
    cx: &Cx<'_>,
) -> Result<Film, EulerRenderShardingError> {
    checkpoint(cx)?;
    let indexed_inputs = index_frame_inputs(plan, scene, inputs, cx)?;
    let segment = *plan.segments.get(segment_plan_index).ok_or(
        EulerRenderShardingError::UnknownSegmentIndex(segment_plan_index),
    )?;
    let start = usize::try_from(segment.first_shard)
        .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("segment shard start"))?;
    let end_u64 = segment.first_shard.checked_add(segment.shard_count).ok_or(
        EulerRenderShardingError::ArithmeticOverflow("segment shard end"),
    )?;
    let end = usize::try_from(end_u64)
        .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("segment shard end"))?;
    let selected = plan
        .shards
        .get(start..end)
        .ok_or(EulerRenderShardingError::Codec(
            "segment shard range outside plan",
        ))?;

    let mut expected_by_logical = BTreeMap::new();
    for (local_index, shard) in selected.iter().enumerate() {
        if expected_by_logical
            .insert(shard.logical_shard_identity, local_index)
            .is_some()
        {
            return Err(EulerRenderShardingError::Codec(
                "duplicate logical shard identity in plan",
            ));
        }
    }

    let mut unique_by_logical = BTreeMap::new();
    for artifact in artifacts {
        checkpoint(cx)?;
        if !expected_by_logical.contains_key(&artifact.logical_shard_identity) {
            return Err(EulerRenderShardingError::ShardResultBindingMismatch);
        }
        match unique_by_logical.get(&artifact.logical_shard_identity) {
            Some(previous) if previous != artifact => {
                return Err(EulerRenderShardingError::ConflictingShardResult(
                    artifact.logical_shard_identity,
                ));
            }
            Some(_) => {}
            None => {
                unique_by_logical.insert(artifact.logical_shard_identity, *artifact);
            }
        }
    }
    for shard in selected {
        if !unique_by_logical.contains_key(&shard.logical_shard_identity) {
            return Err(EulerRenderShardingError::MissingArtifact(
                shard.logical_shard_identity,
            ));
        }
    }

    // One canonical result byte stream can only bind one exact renderer spec.
    // Refuse an untrusted hash aliased to multiple logical shards before read.
    let mut logical_by_hash = BTreeMap::new();
    for artifact in unique_by_logical.values() {
        match logical_by_hash.insert(artifact.artifact_hash, artifact.logical_shard_identity) {
            Some(previous) if previous != artifact.logical_shard_identity => {
                return Err(EulerRenderShardingError::ShardResultBindingMismatch);
            }
            _ => {}
        }
    }

    let mut aggregate_bytes = 0_u64;
    for artifact_hash in logical_by_hash.keys() {
        checkpoint(cx)?;
        let info = ledger
            .artifact_info(artifact_hash)?
            .ok_or(EulerRenderShardingError::MissingArtifact(*artifact_hash))?;
        if info.kind != EULER_RENDER_SHARD_RESULT_ARTIFACT_KIND {
            return Err(EulerRenderShardingError::ForeignArtifactKind {
                hash: *artifact_hash,
                actual: info.kind,
            });
        }
        if info.len > plan.limits.max_result_bytes_per_shard {
            return Err(EulerRenderShardingError::ResultByteLimit {
                limit: plan.limits.max_result_bytes_per_shard,
                observed: info.len,
            });
        }
        aggregate_bytes = aggregate_bytes.checked_add(info.len).ok_or(
            EulerRenderShardingError::ArithmeticOverflow("aggregate result bytes"),
        )?;
        if aggregate_bytes > plan.limits.max_aggregate_result_bytes {
            return Err(EulerRenderShardingError::AggregateResultByteLimit {
                limit: plan.limits.max_aggregate_result_bytes,
                observed: aggregate_bytes,
            });
        }
    }

    let first_selected = *selected.first().ok_or(EulerRenderShardingError::Codec(
        "segment has no planned shards",
    ))?;
    let semantics = bind_segment_semantics(plan, scene, &indexed_inputs, first_selected, cx)?;
    let mut expected_specs = try_vec_capacity(segment.shard_count, "expected shard specs")?;
    for shard in selected {
        checkpoint(cx)?;
        if shard.frame_ordinal != first_selected.frame_ordinal
            || shard.segment_index != first_selected.segment_index
            || shard.frame_identity != first_selected.frame_identity
        {
            return Err(EulerRenderShardingError::Codec(
                "segment shard range crosses a frame binding",
            ));
        }
        expected_specs.push(renderer_spec_for_shard(plan, *shard, semantics)?);
    }
    let mut results: Vec<Option<UniformRenderShardResult>> =
        try_vec_capacity(segment.shard_count, "decoded shard results")?;
    results.resize_with(selected.len(), || None);
    for artifact in unique_by_logical.values() {
        checkpoint(cx)?;
        let local_index = expected_by_logical[&artifact.logical_shard_identity];
        let expected_spec = &expected_specs[local_index];
        let bytes = ledger
            .get_artifact_bounded(
                &artifact.artifact_hash,
                plan.limits.max_result_bytes_per_shard,
            )?
            .ok_or(EulerRenderShardingError::MissingArtifact(
                artifact.artifact_hash,
            ))?;
        let decoded = UniformRenderShardResult::decode_canonical(
            &bytes,
            plan.limits.max_result_bytes_per_shard,
            expected_spec,
            plan.plan_identity,
            expected_spec.shard_identity(),
            cx,
        )?;
        match &results[local_index] {
            Some(previous) if previous.result_identity() != decoded.result_identity() => {
                return Err(EulerRenderShardingError::ConflictingShardResult(
                    artifact.logical_shard_identity,
                ));
            }
            Some(_) => {}
            None => results[local_index] = Some(decoded),
        }
    }
    let mut complete_results = try_vec_capacity(segment.shard_count, "complete shard results")?;
    for (local_index, result) in results.into_iter().enumerate() {
        complete_results.push(result.ok_or(EulerRenderShardingError::MissingArtifact(
            selected[local_index].logical_shard_identity,
        ))?);
    }
    let merge_limits = RenderShardMergeLimits::try_new(
        plan.limits.max_aggregate_result_bytes,
        plan.limits.max_aggregate_result_bytes,
    )?;
    Ok(merge_uniform_shards(
        &expected_specs,
        &complete_results,
        merge_limits,
        cx,
    )?)
}

struct BoundShardSpec {
    spec: UniformRenderShardSpec,
    shutter: fs_render::motion::ShutterInterval,
    cut_side: fs_render::camera::CutSide,
}

#[derive(Clone, Copy)]
struct BoundSegmentSemantics {
    shutter: fs_render::motion::ShutterInterval,
    cut_side: fs_render::camera::CutSide,
    time_mode: FilmTimeMode,
}

fn index_frame_inputs<'prepared>(
    plan: &EulerUniformRenderPlan,
    scene: &EulerCinematicScene<'_>,
    inputs: &[EulerRenderFrameInput<'prepared>],
    cx: &Cx<'_>,
) -> Result<BTreeMap<u64, &'prepared EulerPreparedFrame>, EulerRenderShardingError> {
    checkpoint(cx)?;
    validate_scene_binding(plan, scene)?;
    let observed = u64::try_from(inputs.len())
        .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("frame input count"))?;
    if observed > plan.limits.max_frames {
        return Err(EulerRenderShardingError::FrameLimit {
            limit: plan.limits.max_frames,
            observed,
        });
    }
    if inputs.len() != plan.frames.len() {
        return Err(EulerRenderShardingError::FrameBindingMismatch);
    }
    let mut indexed = BTreeMap::new();
    for input in inputs {
        checkpoint(cx)?;
        if indexed
            .insert(input.frame_ordinal, input.prepared)
            .is_some()
        {
            return Err(EulerRenderShardingError::DuplicateFrameOrdinal(
                input.frame_ordinal,
            ));
        }
    }
    for frame in &plan.frames {
        checkpoint(cx)?;
        if !indexed.contains_key(&frame.frame_ordinal) {
            return Err(EulerRenderShardingError::FrameBindingMismatch);
        }
    }
    Ok(indexed)
}

fn bind_shard_spec_from_index(
    plan: &EulerUniformRenderPlan,
    scene: &EulerCinematicScene<'_>,
    inputs: &BTreeMap<u64, &EulerPreparedFrame>,
    shard_index: usize,
    cx: &Cx<'_>,
) -> Result<BoundShardSpec, EulerRenderShardingError> {
    checkpoint(cx)?;
    validate_scene_binding(plan, scene)?;
    let shard = *plan
        .shards
        .get(shard_index)
        .ok_or(EulerRenderShardingError::UnknownShardIndex(shard_index))?;
    let semantics = bind_segment_semantics(plan, scene, inputs, shard, cx)?;
    let spec = renderer_spec_for_shard(plan, shard, semantics)?;
    Ok(BoundShardSpec {
        spec,
        shutter: semantics.shutter,
        cut_side: semantics.cut_side,
    })
}

fn bind_segment_semantics(
    plan: &EulerUniformRenderPlan,
    scene: &EulerCinematicScene<'_>,
    inputs: &BTreeMap<u64, &EulerPreparedFrame>,
    shard: EulerRenderPlannedShard,
    cx: &Cx<'_>,
) -> Result<BoundSegmentSemantics, EulerRenderShardingError> {
    checkpoint(cx)?;
    let prepared = inputs
        .get(&shard.frame_ordinal)
        .copied()
        .ok_or(EulerRenderShardingError::FrameBindingMismatch)?;
    let segment_index = usize::try_from(shard.segment_index)
        .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("segment_index"))?;
    let (shutter, cut_side) = scene.prepared_segment_shard_binding(prepared, segment_index)?;
    let actual_frame_identity = euler_render_checkpoint_frame_identity(prepared, segment_index)?;
    if actual_frame_identity != shard.frame_identity {
        return Err(EulerRenderShardingError::FrameBindingMismatch);
    }
    let exposure = scene
        .camera()
        .admit_shutter(cx, shutter, cut_side)
        .map_err(EulerSceneError::from)?;
    let time_mode = FilmTimeMode::Cinematic {
        shutter,
        stream_identity: plan.settings.seed,
        shot_id: exposure.shot_id(),
    };
    Ok(BoundSegmentSemantics {
        shutter,
        cut_side,
        time_mode,
    })
}

fn renderer_spec_for_shard(
    plan: &EulerUniformRenderPlan,
    shard: EulerRenderPlannedShard,
    semantics: BoundSegmentSemantics,
) -> Result<UniformRenderShardSpec, EulerRenderShardingError> {
    let renderer_limits = RenderShardLimits::try_new(
        plan.limits.max_paths_per_shard,
        plan.limits.max_result_bytes_per_shard,
    )?;
    let spec = UniformRenderShardSpec::try_new(
        plan.plan_identity,
        shard.frame_identity,
        shard.frame_ordinal,
        plan.settings,
        semantics.time_mode,
        plan.tile_layout()?,
        shard.tile_start,
        shard.tile_end,
        shard.sample_start,
        shard.sample_end,
        renderer_limits,
    )?;
    if spec.path_count() != shard.path_count {
        return Err(EulerRenderShardingError::ShardResultBindingMismatch);
    }
    Ok(spec)
}

fn validate_scene_binding(
    plan: &EulerUniformRenderPlan,
    scene: &EulerCinematicScene<'_>,
) -> Result<(), EulerRenderShardingError> {
    if scene.source_trajectory_identity() != plan.source_trajectory_identity {
        return Err(EulerRenderShardingError::SceneBindingMismatch(
            "source trajectory identity",
        ));
    }
    if scene.source_configuration_identity() != plan.source_configuration_identity {
        return Err(EulerRenderShardingError::SceneBindingMismatch(
            "source configuration identity",
        ));
    }
    if scene.scene_identity() != plan.scene_identity {
        return Err(EulerRenderShardingError::SceneBindingMismatch(
            "scene identity",
        ));
    }
    Ok(())
}

#[derive(Clone, Copy)]
struct CanonicalPlanInputs {
    sequence_identity: ContentHash,
    source_trajectory_identity: ContentHash,
    source_configuration_identity: ContentHash,
    scene_identity: ContentHash,
    settings: Settings,
    tile_width: u32,
    tile_height: u32,
    tiles_per_shard: u64,
    samples_per_shard: u32,
    finishing_neighbor_radius: u32,
    limits: EulerRenderShardLimits,
}

#[derive(Clone, Copy)]
struct RawSegment {
    frame_ordinal: u64,
    frame_position: u64,
    segment_index: u64,
    frame_identity: ContentHash,
}

fn validate_canonical_tables(
    frames: &[EulerRenderPlannedFrame],
    segments: &[RawSegment],
) -> Result<(), EulerRenderShardingError> {
    let mut segment_cursor = 0_u64;
    for (frame_position, frame) in frames.iter().enumerate() {
        if frame.segment_count == 0 || frame.first_segment != segment_cursor {
            return Err(EulerRenderShardingError::Codec(
                "noncanonical frame segment range",
            ));
        }
        if frame_position > 0 && frames[frame_position - 1].frame_ordinal >= frame.frame_ordinal {
            return Err(EulerRenderShardingError::Codec(
                "frames not strictly sorted by ordinal",
            ));
        }
        let end = segment_cursor.checked_add(frame.segment_count).ok_or(
            EulerRenderShardingError::ArithmeticOverflow("frame segment end"),
        )?;
        let start = usize::try_from(segment_cursor)
            .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("segment cursor"))?;
        let end_usize = usize::try_from(end)
            .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("segment end"))?;
        let frame_segments =
            segments
                .get(start..end_usize)
                .ok_or(EulerRenderShardingError::Codec(
                    "frame segment range outside table",
                ))?;
        for (segment_index, segment) in frame_segments.iter().enumerate() {
            require_nonzero_identity("frame_identity", segment.frame_identity)?;
            if segment.frame_ordinal != frame.frame_ordinal
                || segment.frame_position != frame_position as u64
                || segment.segment_index != segment_index as u64
            {
                return Err(EulerRenderShardingError::Codec(
                    "noncanonical frame-segment ordering",
                ));
            }
        }
        segment_cursor = end;
    }
    if usize::try_from(segment_cursor).ok() != Some(segments.len()) {
        return Err(EulerRenderShardingError::Codec(
            "unowned segment table suffix",
        ));
    }
    Ok(())
}

fn put_u16(bytes: &mut Vec<u8>, value: u16) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn put_u32(bytes: &mut Vec<u8>, value: u32) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn put_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn put_hash(bytes: &mut Vec<u8>, value: ContentHash) {
    bytes.extend_from_slice(value.as_bytes());
}

fn put_settings(bytes: &mut Vec<u8>, settings: Settings) {
    put_u32(bytes, settings.width);
    put_u32(bytes, settings.height);
    put_u32(bytes, settings.spp);
    put_u32(bytes, settings.max_depth);
    bytes.push(sampler_tag(settings.sampler));
    bytes.push(strategy_tag(settings.strategy));
    put_u16(bytes, 0);
    put_u64(bytes, settings.seed);
}

fn put_limits(bytes: &mut Vec<u8>, limits: EulerRenderShardLimits) {
    put_u64(bytes, limits.max_frames);
    put_u64(bytes, limits.max_shards);
    put_u64(bytes, limits.max_plan_bytes);
    put_u64(bytes, limits.max_paths_per_shard);
    put_u64(bytes, limits.max_result_bytes_per_shard);
    put_u64(bytes, limits.max_aggregate_result_bytes);
}

struct PlanReader<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> PlanReader<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    fn take(&mut self, len: usize) -> Result<&'a [u8], EulerRenderShardingError> {
        let end = self
            .offset
            .checked_add(len)
            .ok_or(EulerRenderShardingError::Codec("plan offset overflow"))?;
        let value = self
            .bytes
            .get(self.offset..end)
            .ok_or(EulerRenderShardingError::Codec("truncated plan"))?;
        self.offset = end;
        Ok(value)
    }

    fn u16(&mut self) -> Result<u16, EulerRenderShardingError> {
        let bytes: [u8; 2] = self
            .take(2)?
            .try_into()
            .map_err(|_| EulerRenderShardingError::Codec("invalid u16"))?;
        Ok(u16::from_le_bytes(bytes))
    }

    fn u32(&mut self) -> Result<u32, EulerRenderShardingError> {
        let bytes: [u8; 4] = self
            .take(4)?
            .try_into()
            .map_err(|_| EulerRenderShardingError::Codec("invalid u32"))?;
        Ok(u32::from_le_bytes(bytes))
    }

    fn u64(&mut self) -> Result<u64, EulerRenderShardingError> {
        let bytes: [u8; 8] = self
            .take(8)?
            .try_into()
            .map_err(|_| EulerRenderShardingError::Codec("invalid u64"))?;
        Ok(u64::from_le_bytes(bytes))
    }

    fn hash(&mut self) -> Result<ContentHash, EulerRenderShardingError> {
        ContentHash::from_slice(self.take(32)?)
            .ok_or(EulerRenderShardingError::Codec("invalid content hash"))
    }

    fn settings(&mut self) -> Result<Settings, EulerRenderShardingError> {
        let width = self.u32()?;
        let height = self.u32()?;
        let spp = self.u32()?;
        let max_depth = self.u32()?;
        let sampler = match self.take(1)?[0] {
            0 => Sampler::Iid,
            1 => Sampler::OwenSobol,
            _ => return Err(EulerRenderShardingError::Codec("unknown sampler tag")),
        };
        let strategy = match self.take(1)?[0] {
            0 => DirectStrategy::NeeOnly,
            1 => DirectStrategy::BsdfOnly,
            2 => DirectStrategy::Mis,
            _ => {
                return Err(EulerRenderShardingError::Codec(
                    "unknown direct-strategy tag",
                ));
            }
        };
        if self.u16()? != 0 {
            return Err(EulerRenderShardingError::Codec(
                "nonzero reserved settings bytes",
            ));
        }
        let seed = self.u64()?;
        Ok(Settings {
            width,
            height,
            spp,
            max_depth,
            sampler,
            strategy,
            seed,
        })
    }

    fn limits(&mut self) -> Result<EulerRenderShardLimits, EulerRenderShardingError> {
        EulerRenderShardLimits::try_new(
            self.u64()?,
            self.u64()?,
            self.u64()?,
            self.u64()?,
            self.u64()?,
            self.u64()?,
        )
    }

    fn is_finished(&self) -> bool {
        self.offset == self.bytes.len()
    }
}

fn partition_segments(
    plan_identity: ContentHash,
    raw_segments: &[RawSegment],
    layout: RenderTileLayout,
    spp: u32,
    tiles_per_shard: u64,
    samples_per_shard: u32,
    limits: EulerRenderShardLimits,
    cx: &Cx<'_>,
) -> Result<
    (
        Vec<EulerRenderPlannedSegment>,
        Vec<EulerRenderPlannedShard>,
        u64,
    ),
    EulerRenderShardingError,
> {
    let renderer_limits = RenderShardLimits::try_new(
        limits.max_paths_per_shard,
        limits.max_result_bytes_per_shard,
    )?;
    let tile_blocks = div_ceil_u64(layout.tile_count(), tiles_per_shard)?;
    let sample_blocks = div_ceil_u64(u64::from(spp), u64::from(samples_per_shard))?;
    let shards_per_segment = tile_blocks.checked_mul(sample_blocks).ok_or(
        EulerRenderShardingError::ArithmeticOverflow("shards_per_segment"),
    )?;
    let shard_count = u64::try_from(raw_segments.len())
        .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("segment_count"))?
        .checked_mul(shards_per_segment)
        .ok_or(EulerRenderShardingError::ArithmeticOverflow("shard_count"))?;
    if shard_count > limits.max_shards {
        return Err(EulerRenderShardingError::ShardLimit {
            limit: limits.max_shards,
            observed: shard_count,
        });
    }
    // Prove every worker and per-segment coordinator cap before reserving the
    // full shard table. The layout and partition are identical for each
    // segment, so one cancellable pass is sufficient.
    let mut paths_per_segment = 0_u64;
    let mut result_bytes_per_segment = 0_u64;
    let mut preflight_tile_start = 0_u64;
    while preflight_tile_start < layout.tile_count() {
        checkpoint(cx)?;
        let preflight_tile_end = preflight_tile_start
            .checked_add(tiles_per_shard)
            .unwrap_or(u64::MAX)
            .min(layout.tile_count());
        let block_pixels =
            tile_block_pixel_count(layout, preflight_tile_start, preflight_tile_end)?;
        let result_bytes = renderer_limits.admit_result_pixels(block_pixels)?;
        let mut preflight_sample_start = 0_u32;
        while preflight_sample_start < spp {
            checkpoint(cx)?;
            let preflight_sample_end = preflight_sample_start
                .saturating_add(samples_per_shard)
                .min(spp);
            let path_count = block_pixels
                .checked_mul(u64::from(preflight_sample_end - preflight_sample_start))
                .ok_or(EulerRenderShardingError::ArithmeticOverflow(
                    "shard path_count",
                ))?;
            if path_count > limits.max_paths_per_shard {
                return Err(EulerRenderShardingError::PathLimit {
                    limit: limits.max_paths_per_shard,
                    observed: path_count,
                });
            }
            paths_per_segment = paths_per_segment.checked_add(path_count).ok_or(
                EulerRenderShardingError::ArithmeticOverflow("segment paths"),
            )?;
            result_bytes_per_segment = result_bytes_per_segment.checked_add(result_bytes).ok_or(
                EulerRenderShardingError::ArithmeticOverflow("segment result bytes"),
            )?;
            if result_bytes_per_segment > limits.max_aggregate_result_bytes {
                return Err(EulerRenderShardingError::AggregateResultByteLimit {
                    limit: limits.max_aggregate_result_bytes,
                    observed: result_bytes_per_segment,
                });
            }
            preflight_sample_start = preflight_sample_end;
        }
        preflight_tile_start = preflight_tile_end;
    }
    let total_paths = paths_per_segment
        .checked_mul(
            u64::try_from(raw_segments.len())
                .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("segment_count"))?,
        )
        .ok_or(EulerRenderShardingError::ArithmeticOverflow("total_paths"))?;

    let mut segments = try_vec_capacity(
        u64::try_from(raw_segments.len())
            .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("segment_count"))?,
        "planned segments",
    )?;
    let mut shards = try_vec_capacity(shard_count, "planned shards")?;

    for raw in raw_segments {
        checkpoint(cx)?;
        let first_shard = u64::try_from(shards.len())
            .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("first_shard"))?;
        let mut tile_start = 0_u64;
        while tile_start < layout.tile_count() {
            checkpoint(cx)?;
            let tile_end = tile_start
                .checked_add(tiles_per_shard)
                .unwrap_or(u64::MAX)
                .min(layout.tile_count());
            let block_pixels = tile_block_pixel_count(layout, tile_start, tile_end)?;
            let mut sample_start = 0_u32;
            while sample_start < spp {
                checkpoint(cx)?;
                let sample_end = sample_start.saturating_add(samples_per_shard).min(spp);
                let sample_count = u64::from(sample_end - sample_start);
                let path_count = block_pixels.checked_mul(sample_count).ok_or(
                    EulerRenderShardingError::ArithmeticOverflow("shard path_count"),
                )?;
                let shard_ordinal = u64::try_from(shards.len())
                    .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("shard_ordinal"))?;
                let logical_shard_identity = logical_shard_identity(
                    plan_identity,
                    shard_ordinal,
                    *raw,
                    tile_start,
                    tile_end,
                    sample_start,
                    sample_end,
                );
                shards.push(EulerRenderPlannedShard {
                    shard_ordinal,
                    logical_shard_identity,
                    frame_ordinal: raw.frame_ordinal,
                    segment_index: raw.segment_index,
                    frame_identity: raw.frame_identity,
                    tile_start,
                    tile_end,
                    sample_start,
                    sample_end,
                    path_count,
                });
                sample_start = sample_end;
            }
            tile_start = tile_end;
        }
        segments.push(EulerRenderPlannedSegment {
            frame_ordinal: raw.frame_ordinal,
            frame_position: raw.frame_position,
            segment_index: raw.segment_index,
            frame_identity: raw.frame_identity,
            first_shard,
            shard_count: shards_per_segment,
        });
    }
    Ok((segments, shards, total_paths))
}

fn plan_identity(
    canonical: &CanonicalPlanInputs,
    frames: &[EulerRenderPlannedFrame],
    segments: &[RawSegment],
    cx: &Cx<'_>,
) -> Result<ContentHash, EulerRenderShardingError> {
    let mut hasher = DomainHasher::new(EULER_RENDER_SHARD_PLAN_IDENTITY_DOMAIN);
    hasher.update(&EULER_RENDER_SHARD_PLAN_SCHEMA_VERSION.to_le_bytes());
    hash_canonical_inputs(&mut hasher, canonical);
    hasher.update(
        &u64::try_from(frames.len())
            .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("frame_count"))?
            .to_le_bytes(),
    );
    hasher.update(
        &u64::try_from(segments.len())
            .map_err(|_| EulerRenderShardingError::ArithmeticOverflow("segment_count"))?
            .to_le_bytes(),
    );
    for frame in frames {
        checkpoint(cx)?;
        hasher.update(&frame.frame_ordinal.to_le_bytes());
        hasher.update(&frame.first_segment.to_le_bytes());
        hasher.update(&frame.segment_count.to_le_bytes());
    }
    for segment in segments {
        checkpoint(cx)?;
        hasher.update(&segment.frame_ordinal.to_le_bytes());
        hasher.update(&segment.frame_position.to_le_bytes());
        hasher.update(&segment.segment_index.to_le_bytes());
        hasher.update(segment.frame_identity.as_bytes());
    }
    checkpoint(cx)?;
    Ok(hasher.finalize())
}

fn hash_canonical_inputs(hasher: &mut DomainHasher, canonical: &CanonicalPlanInputs) {
    hasher.update(canonical.sequence_identity.as_bytes());
    hasher.update(canonical.source_trajectory_identity.as_bytes());
    hasher.update(canonical.source_configuration_identity.as_bytes());
    hasher.update(canonical.scene_identity.as_bytes());
    hash_settings(hasher, canonical.settings);
    hasher.update(&canonical.tile_width.to_le_bytes());
    hasher.update(&canonical.tile_height.to_le_bytes());
    hasher.update(&canonical.tiles_per_shard.to_le_bytes());
    hasher.update(&canonical.samples_per_shard.to_le_bytes());
    hasher.update(&canonical.finishing_neighbor_radius.to_le_bytes());
    hasher.update(&canonical.limits.max_frames.to_le_bytes());
    hasher.update(&canonical.limits.max_shards.to_le_bytes());
    hasher.update(&canonical.limits.max_plan_bytes.to_le_bytes());
    hasher.update(&canonical.limits.max_paths_per_shard.to_le_bytes());
    hasher.update(&canonical.limits.max_result_bytes_per_shard.to_le_bytes());
    hasher.update(&canonical.limits.max_aggregate_result_bytes.to_le_bytes());
}

fn logical_shard_identity(
    plan_identity: ContentHash,
    shard_ordinal: u64,
    segment: RawSegment,
    tile_start: u64,
    tile_end: u64,
    sample_start: u32,
    sample_end: u32,
) -> ContentHash {
    let mut hasher = DomainHasher::new(EULER_RENDER_LOGICAL_SHARD_IDENTITY_DOMAIN);
    hasher.update(&EULER_RENDER_SHARD_PLAN_SCHEMA_VERSION.to_le_bytes());
    hasher.update(plan_identity.as_bytes());
    hasher.update(&shard_ordinal.to_le_bytes());
    hasher.update(&segment.frame_ordinal.to_le_bytes());
    hasher.update(&segment.segment_index.to_le_bytes());
    hasher.update(segment.frame_identity.as_bytes());
    hasher.update(&tile_start.to_le_bytes());
    hasher.update(&tile_end.to_le_bytes());
    hasher.update(&sample_start.to_le_bytes());
    hasher.update(&sample_end.to_le_bytes());
    hasher.finalize()
}

fn tile_block_pixel_count(
    layout: RenderTileLayout,
    tile_start: u64,
    tile_end: u64,
) -> Result<u64, EulerRenderShardingError> {
    layout
        .pixel_count_in_range(tile_start, tile_end)
        .ok_or(EulerRenderShardingError::InvalidPartition("tile range"))
}

fn measured_plan_bytes(
    frames: u64,
    segments: u64,
    shards: u64,
) -> Result<u64, EulerRenderShardingError> {
    PLAN_HEADER_BYTES
        .checked_add(frames.checked_mul(PLAN_FRAME_BYTES).ok_or(
            EulerRenderShardingError::ArithmeticOverflow("plan frame bytes"),
        )?)
        .and_then(|value| value.checked_add(segments.checked_mul(PLAN_SEGMENT_BYTES)?))
        .and_then(|value| value.checked_add(shards.checked_mul(PLAN_SHARD_BYTES)?))
        .ok_or(EulerRenderShardingError::ArithmeticOverflow(
            "plan byte length",
        ))
}

fn div_ceil_u64(value: u64, divisor: u64) -> Result<u64, EulerRenderShardingError> {
    if divisor == 0 {
        return Err(EulerRenderShardingError::InvalidPartition("divisor"));
    }
    Ok(value / divisor + u64::from(value % divisor != 0))
}

fn try_vec_capacity<T>(
    count: u64,
    field: &'static str,
) -> Result<Vec<T>, EulerRenderShardingError> {
    let count =
        usize::try_from(count).map_err(|_| EulerRenderShardingError::ArithmeticOverflow(field))?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| EulerRenderShardingError::Capacity(field))?;
    Ok(values)
}

fn require_nonzero_identity(
    field: &'static str,
    identity: ContentHash,
) -> Result<(), EulerRenderShardingError> {
    if identity.as_bytes().iter().all(|byte| *byte == 0) {
        Err(EulerRenderShardingError::ZeroIdentity(field))
    } else {
        Ok(())
    }
}

fn validate_settings(settings: Settings) -> Result<(), EulerRenderShardingError> {
    if settings.width == 0 {
        return Err(EulerRenderShardingError::InvalidSettings("width"));
    }
    if settings.height == 0 {
        return Err(EulerRenderShardingError::InvalidSettings("height"));
    }
    if settings.spp == 0 {
        return Err(EulerRenderShardingError::InvalidSettings("spp"));
    }
    if settings.max_depth == 0 {
        return Err(EulerRenderShardingError::InvalidSettings("max_depth"));
    }
    Ok(())
}

fn hash_settings(hasher: &mut DomainHasher, settings: Settings) {
    hasher.update(&settings.width.to_le_bytes());
    hasher.update(&settings.height.to_le_bytes());
    hasher.update(&settings.spp.to_le_bytes());
    hasher.update(&settings.max_depth.to_le_bytes());
    hasher.update(&[sampler_tag(settings.sampler)]);
    hasher.update(&[strategy_tag(settings.strategy)]);
    hasher.update(&settings.seed.to_le_bytes());
}

const fn sampler_tag(sampler: Sampler) -> u8 {
    match sampler {
        Sampler::Iid => 0,
        Sampler::OwenSobol => 1,
    }
}

const fn strategy_tag(strategy: DirectStrategy) -> u8 {
    match strategy {
        DirectStrategy::NeeOnly => 0,
        DirectStrategy::BsdfOnly => 1,
        DirectStrategy::Mis => 2,
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), EulerRenderShardingError> {
    cx.checkpoint()
        .map_err(|_| EulerRenderShardingError::Cancelled)
}

/// Structured fail-closed refusal from plan, worker, codec, artifact, or merge
/// coordination. Errors never contain a partial film.
#[derive(Debug)]
pub enum EulerRenderShardingError {
    /// Execution scope requested cancellation.
    Cancelled,
    /// A mandatory identity was the all-zero sentinel.
    ZeroIdentity(&'static str),
    /// A caller resource cap was zero.
    InvalidLimit(&'static str),
    /// Uniform settings were outside the admitted fixed-SPP domain.
    InvalidSettings(&'static str),
    /// Tile/sample partition metadata was invalid.
    InvalidPartition(&'static str),
    /// No logical input frame was supplied.
    EmptyPlan,
    /// One prepared frame contained no segment.
    EmptyPreparedFrame(u64),
    /// Two caller entries named the same stable ordinal.
    DuplicateFrameOrdinal(u64),
    /// Frame count exceeded the explicit cap.
    FrameLimit {
        /// Admitted maximum frames.
        limit: u64,
        /// Required frames.
        observed: u64,
    },
    /// Shard count exceeded the explicit cap.
    ShardLimit {
        /// Admitted maximum shards.
        limit: u64,
        /// Required shards.
        observed: u64,
    },
    /// Canonical plan length exceeded the explicit cap.
    PlanByteLimit {
        /// Admitted maximum bytes.
        limit: u64,
        /// Required or supplied bytes.
        observed: u64,
    },
    /// One shard exceeded its traced-path cap.
    PathLimit {
        /// Admitted maximum paths.
        limit: u64,
        /// Required paths.
        observed: u64,
    },
    /// One encoded shard result exceeded its cap.
    ResultByteLimit {
        /// Admitted maximum bytes.
        limit: u64,
        /// Required or supplied bytes.
        observed: u64,
    },
    /// Stored result metadata exceeded the aggregate read cap.
    AggregateResultByteLimit {
        /// Admitted maximum bytes.
        limit: u64,
        /// Required bytes.
        observed: u64,
    },
    /// Checked integer arithmetic refused an unrepresentable plan.
    ArithmeticOverflow(&'static str),
    /// Fallible retained allocation refused the bounded request.
    Capacity(&'static str),
    /// Caller selected no such plan segment.
    UnknownSegmentIndex(usize),
    /// Caller selected no such plan shard.
    UnknownShardIndex(usize),
    /// A supplied prepared frame did not match the pinned shard.
    FrameBindingMismatch,
    /// Source trajectory/configuration/scene did not match the plan.
    SceneBindingMismatch(&'static str),
    /// Renderer sharding refused execution, codec, or merge.
    Renderer(RenderShardError),
    /// Euler scene admission refused the binding.
    Scene(EulerSceneError),
    /// Checkpoint frame identity construction refused the binding.
    Checkpoint(EulerRenderCheckpointError),
    /// Design Ledger storage or integrity check failed.
    Ledger(LedgerError),
    /// Canonical plan bytes were malformed.
    Codec(&'static str),
    /// Canonical plan pin did not match.
    PlanIdentityMismatch {
        /// Trusted external plan identity.
        expected: ContentHash,
        /// Identity carried or recomputed from the plan.
        actual: ContentHash,
    },
    /// One immutable artifact was missing.
    MissingArtifact(ContentHash),
    /// One stored artifact used another artifact kind.
    ForeignArtifactKind {
        /// Content hash of the rejected artifact.
        hash: ContentHash,
        /// Stored artifact kind.
        actual: String,
    },
    /// The same logical shard was submitted with differing valid result bytes.
    ConflictingShardResult(ContentHash),
    /// A result artifact did not bind the selected logical shard.
    ShardResultBindingMismatch,
}

impl fmt::Display for EulerRenderShardingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for EulerRenderShardingError {}

impl From<RenderShardError> for EulerRenderShardingError {
    fn from(error: RenderShardError) -> Self {
        Self::Renderer(error)
    }
}

impl From<EulerSceneError> for EulerRenderShardingError {
    fn from(error: EulerSceneError) -> Self {
        Self::Scene(error)
    }
}

impl From<EulerRenderCheckpointError> for EulerRenderShardingError {
    fn from(error: EulerRenderCheckpointError) -> Self {
        Self::Checkpoint(error)
    }
}

impl From<LedgerError> for EulerRenderShardingError {
    fn from(error: LedgerError) -> Self {
        Self::Ledger(error)
    }
}
