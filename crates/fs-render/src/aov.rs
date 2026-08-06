//! Aligned cinematic diagnostic and denoising AOV accumulation.
//!
//! This module is deliberately opt-in and parallel to [`crate::tracer::Film`].
//! The legacy RGB film, EXR bytes, checkpoints, and shard formats are not
//! extended in place. Every diagnostic sample is instead taken from the exact
//! accepted primary intersection that produced its beauty contribution.
//!
//! Categorical identities use deterministic per-frame palettes. EXR stores the
//! small exact palette index as `FLOAT` and stores the lossless `u64`/BLAKE3
//! mapping in header attributes; a hash is never rounded through `f32`.

use core::fmt;
use core::mem::size_of;
use fs_blake3::{ContentHash, DomainHasher};
use fs_img::{Channel, ExrAttribute, PixelType};

use crate::motion::ShutterInterval;
use crate::motion_vectors::MOTION_VECTOR_SEMANTICS_VERSION;
use crate::spectral::{xyz_e_to_d65, xyz_to_linear_srgb};
use crate::tracer::{
    ADAPTIVE_SAMPLING_SEMANTICS_VERSION, AdaptiveDecision, AdaptiveFilm, AdaptiveSamplingConfig,
    CINEMATIC_CAMERA_TRACER_BIT_SEMANTICS_VERSION, DIELECTRIC_TRACER_BIT_SEMANTICS_VERSION,
    DirectStrategy, Film, LIGHTING_TRACER_BIT_SEMANTICS_VERSION, MATERIAL_CONTENT_IDENTITY_DOMAIN,
    MOTION_TRACER_BIT_SEMANTICS_VERSION, Sampler, Scene, Settings, Shape,
    TRACER_BIT_SEMANTICS_VERSION,
};

/// Bit-affecting channel, accumulation, invalid-value, and palette semantics.
pub const CINEMATIC_AOV_SEMANTICS_VERSION: u32 = 1;
/// Deterministic nearest-primary categorical selection semantics.
pub const CINEMATIC_AOV_CATEGORY_SEMANTICS_VERSION: u32 = 1;
/// Linear-sRGB material-albedo extraction semantics.
pub const CINEMATIC_AOV_ALBEDO_SEMANTICS_VERSION: u32 = 1;
/// Domain for complete cinematic AOV configuration identities.
pub const CINEMATIC_AOV_CONFIG_IDENTITY_DOMAIN: &str =
    "org.frankensim.render.cinematic-aov-config.v1";

const MAX_EXACT_F32_INTEGER: u32 = 1 << 24;
const STRING_ATTRIBUTE_TYPE: &str = "string";

/// Frozen AOV channel bundles.
///
/// The channel counts include `R`, `G`, and `B`. `FinalDiagnostic` retains the
/// direct/indirect/emission split promised by the cinematic quality profile in
/// addition to the denoising core and stable-ID diagnostics.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum CinematicAovProfile {
    /// Linear-sRGB beauty only; useful for exercising the new artifact binding
    /// while producing no diagnostic planes.
    BeautyOnly = 0,
    /// Beauty, albedo, shading normal, axial depth, primary coverage, raw-Y
    /// variance, and current-to-previous raster motion.
    DailyCore = 1,
    /// Daily core plus geometric normal, direct/indirect/emission beauty
    /// decomposition, stable object/material palette indices, sample count,
    /// and an explicit validity bit mask.
    FinalDiagnostic = 2,
}

impl CinematicAovProfile {
    /// Number of single-precision EXR channels in this exact profile.
    #[must_use]
    pub const fn float_channel_count(self) -> u32 {
        match self {
            Self::BeautyOnly => 3,
            Self::DailyCore => 14,
            Self::FinalDiagnostic => 30,
        }
    }

    const fn has_common(self) -> bool {
        !matches!(self, Self::BeautyOnly)
    }

    const fn has_final(self) -> bool {
        matches!(self, Self::FinalDiagnostic)
    }

    const fn name(self) -> &'static str {
        match self {
            Self::BeautyOnly => "beauty-only-v1",
            Self::DailyCore => "daily-core-v1",
            Self::FinalDiagnostic => "final-diagnostic-v1",
        }
    }
}

/// Explicit retained/export-memory and palette admission limits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CinematicAovLimits {
    max_pixels: u64,
    max_retained_bytes: u64,
    max_export_plane_bytes: u64,
    max_palette_entries: u32,
}

impl CinematicAovLimits {
    /// Construct nonzero limits. Palette entries must fit exactly in an EXR
    /// `FLOAT` index after reserving index zero for unavailable/background.
    pub fn try_new(
        max_pixels: u64,
        max_retained_bytes: u64,
        max_export_plane_bytes: u64,
        max_palette_entries: u32,
    ) -> Result<Self, CinematicAovError> {
        if max_pixels == 0
            || max_retained_bytes == 0
            || max_export_plane_bytes == 0
            || max_palette_entries == 0
            || max_palette_entries >= MAX_EXACT_F32_INTEGER
        {
            return Err(CinematicAovError::InvalidLimits);
        }
        Ok(Self {
            max_pixels,
            max_retained_bytes,
            max_export_plane_bytes,
            max_palette_entries,
        })
    }

    /// Maximum admitted raster pixels.
    #[must_use]
    pub const fn max_pixels(self) -> u64 {
        self.max_pixels
    }

    /// Maximum retained beauty plus AOV accumulator bytes.
    #[must_use]
    pub const fn max_retained_bytes(self) -> u64 {
        self.max_retained_bytes
    }

    /// Maximum aggregate bytes in staged single-precision EXR planes.
    #[must_use]
    pub const fn max_export_plane_bytes(self) -> u64 {
        self.max_export_plane_bytes
    }

    /// Maximum nonzero entries in either exact identity palette.
    #[must_use]
    pub const fn max_palette_entries(self) -> u32 {
        self.max_palette_entries
    }
}

impl Default for CinematicAovLimits {
    fn default() -> Self {
        // Admits one 4K final-diagnostic frame and its fallible full-frame
        // transactional staging copy. This accounts for owned accumulator
        // payloads, not allocator overhead, EXR encoder scratch, metadata, or
        // process RSS.
        Self {
            max_pixels: 3_840 * 2_160,
            max_retained_bytes: 6 * 1024 * 1024 * 1024,
            max_export_plane_bytes: 1024 * 1024 * 1024,
            max_palette_entries: 65_535,
        }
    }
}

/// L6-supplied identities and exact frame-time context retained in the image
/// artifact. These digests are assertions supplied by the composition layer;
/// `fs-render` does not mint trajectory or scene authority from borrowed data.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CinematicAovProvenance {
    frame_index: u64,
    frame_time_s: f64,
    previous_frame_time_s: f64,
    next_frame_time_s: f64,
    source_trajectory_identity: ContentHash,
    scene_identity: ContentHash,
    composition_identity: ContentHash,
}

impl CinematicAovProvenance {
    /// Admit finite, nondecreasing frame-reference times and nonzero content
    /// identities.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        frame_index: u64,
        frame_time_s: f64,
        previous_frame_time_s: f64,
        next_frame_time_s: f64,
        source_trajectory_identity: ContentHash,
        scene_identity: ContentHash,
        composition_identity: ContentHash,
    ) -> Result<Self, CinematicAovError> {
        if !previous_frame_time_s.is_finite()
            || !frame_time_s.is_finite()
            || !next_frame_time_s.is_finite()
            || previous_frame_time_s > frame_time_s
            || frame_time_s > next_frame_time_s
        {
            return Err(CinematicAovError::InvalidFrameReferences);
        }
        for (field, identity) in [
            ("source_trajectory_identity", source_trajectory_identity),
            ("scene_identity", scene_identity),
            ("composition_identity", composition_identity),
        ] {
            if identity.as_bytes().iter().all(|byte| *byte == 0) {
                return Err(CinematicAovError::MissingIdentity { field });
            }
        }
        Ok(Self {
            frame_index,
            frame_time_s: canonical_zero(frame_time_s),
            previous_frame_time_s: canonical_zero(previous_frame_time_s),
            next_frame_time_s: canonical_zero(next_frame_time_s),
            source_trajectory_identity,
            scene_identity,
            composition_identity,
        })
    }

    /// Master-frame ordinal.
    #[must_use]
    pub const fn frame_index(self) -> u64 {
        self.frame_index
    }

    /// Current frame presentation time in seconds.
    #[must_use]
    pub const fn frame_time_s(self) -> f64 {
        self.frame_time_s
    }

    /// Previous-frame reference time in seconds.
    #[must_use]
    pub const fn previous_frame_time_s(self) -> f64 {
        self.previous_frame_time_s
    }

    /// Next-frame reference time in seconds.
    #[must_use]
    pub const fn next_frame_time_s(self) -> f64 {
        self.next_frame_time_s
    }

    /// Accepted source trajectory artifact identity.
    #[must_use]
    pub const fn source_trajectory_identity(self) -> ContentHash {
        self.source_trajectory_identity
    }

    /// Caller-supplied complete scene identity.
    #[must_use]
    pub const fn scene_identity(self) -> ContentHash {
        self.scene_identity
    }

    /// L6 cinematic composition identity.
    #[must_use]
    pub const fn composition_identity(self) -> ContentHash {
        self.composition_identity
    }
}

/// Complete typed AOV configuration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CinematicAovConfig {
    profile: CinematicAovProfile,
    provenance: CinematicAovProvenance,
    limits: CinematicAovLimits,
    identity: ContentHash,
}

impl CinematicAovConfig {
    /// Bind a frozen profile, external provenance, and explicit resource limits.
    #[must_use]
    pub fn new(
        profile: CinematicAovProfile,
        provenance: CinematicAovProvenance,
        limits: CinematicAovLimits,
    ) -> Self {
        let mut hasher = DomainHasher::new(CINEMATIC_AOV_CONFIG_IDENTITY_DOMAIN);
        hasher.update(&CINEMATIC_AOV_SEMANTICS_VERSION.to_le_bytes());
        hasher.update(&CINEMATIC_AOV_CATEGORY_SEMANTICS_VERSION.to_le_bytes());
        hasher.update(&CINEMATIC_AOV_ALBEDO_SEMANTICS_VERSION.to_le_bytes());
        hasher.update(&[profile as u8]);
        hasher.update(&provenance.frame_index.to_le_bytes());
        for time in [
            provenance.frame_time_s,
            provenance.previous_frame_time_s,
            provenance.next_frame_time_s,
        ] {
            hasher.update(&canonical_zero(time).to_bits().to_le_bytes());
        }
        hasher.update(provenance.source_trajectory_identity.as_bytes());
        hasher.update(provenance.scene_identity.as_bytes());
        hasher.update(provenance.composition_identity.as_bytes());
        hasher.update(&limits.max_pixels.to_le_bytes());
        hasher.update(&limits.max_retained_bytes.to_le_bytes());
        hasher.update(&limits.max_export_plane_bytes.to_le_bytes());
        hasher.update(&limits.max_palette_entries.to_le_bytes());
        Self {
            profile,
            provenance,
            limits,
            identity: hasher.finalize(),
        }
    }

    /// Frozen channel bundle.
    #[must_use]
    pub const fn profile(self) -> CinematicAovProfile {
        self.profile
    }

    /// External source/frame provenance.
    #[must_use]
    pub const fn provenance(self) -> CinematicAovProvenance {
        self.provenance
    }

    /// Explicit resource envelope.
    #[must_use]
    pub const fn limits(self) -> CinematicAovLimits {
        self.limits
    }

    /// Deterministic identity of every configuration field.
    #[must_use]
    pub const fn identity(self) -> ContentHash {
        self.identity
    }

    pub(crate) const fn captures_primary(self) -> bool {
        self.profile.has_common()
    }

    pub(crate) const fn captures_ids(self) -> bool {
        self.profile.has_final()
    }
}

/// Validity flags stored exactly in `diagnostic.validity`.
pub mod validity {
    /// At least one primary surface contributed to surface-filtered planes.
    pub const PRIMARY: u32 = 1 << 0;
    /// At least one material had an admitted linear-sRGB albedo.
    pub const ALBEDO: u32 = 1 << 1;
    /// At least one primary supplied an authored shading normal. The normal AOV
    /// otherwise deliberately falls back to the geometric normal.
    pub const AUTHORED_SHADING_NORMAL: u32 = 1 << 2;
    /// At least one exact current-to-previous projected motion vector existed.
    pub const PREVIOUS_MOTION: u32 = 1 << 3;
    /// Nearest categorical sample had a stable instance object identity.
    pub const OBJECT_ID: u32 = 1 << 4;
    /// Nearest categorical sample had a material identity.
    pub const MATERIAL_ID: u32 = 1 << 5;
    /// Direct/indirect/emission decomposition consumed at least one sample.
    pub const CONTRIBUTION_SPLIT: u32 = 1 << 6;
}

/// One fully prepared AOV observation from the same accepted path sample.
/// Constructed only by the tracer integration seam.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct AlignedAovSample {
    pub beauty_xyz: [f64; 3],
    pub direct_xyz: [f64; 3],
    pub indirect_xyz: [f64; 3],
    pub emission_xyz: [f64; 3],
    pub pixel_jitter: [f64; 2],
    pub absolute_sample: u32,
    pub primary: Option<AlignedAovPrimary>,
}

/// Primary-surface fields after material, camera-depth, motion, and palette
/// evaluation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct AlignedAovPrimary {
    pub primitive_index: usize,
    pub object_palette_index: u32,
    pub material_palette_index: u32,
    pub albedo_linear_rgb: Option<[f64; 3]>,
    pub geometric_normal_world: [f64; 3],
    pub shading_normal_world: [f64; 3],
    pub has_authored_shading_normal: bool,
    pub depth_m: f64,
    pub previous_motion_pixels: Option<[f64; 2]>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct CommonPixel {
    albedo_sum: [f64; 3],
    shading_normal_sum: [f64; 3],
    depth_sum_m: f64,
    previous_motion_sum_pixels: [f64; 2],
    mean_y: f64,
    m2_y: f64,
    accepted_count: u32,
    primary_count: u32,
    albedo_count: u32,
    authored_shading_normal_count: u32,
    previous_motion_count: u32,
}

impl CommonPixel {
    const EMPTY: Self = Self {
        albedo_sum: [0.0; 3],
        shading_normal_sum: [0.0; 3],
        depth_sum_m: 0.0,
        previous_motion_sum_pixels: [0.0; 2],
        mean_y: 0.0,
        m2_y: 0.0,
        accepted_count: 0,
        primary_count: 0,
        albedo_count: 0,
        authored_shading_normal_count: 0,
        previous_motion_count: 0,
    };

    fn push(&mut self, sample: AlignedAovSample) -> Result<(), CinematicAovError> {
        validate_vec("beauty_xyz", sample.beauty_xyz)?;
        let mut next = *self;
        let next_count = checked_increment(next.accepted_count)?;
        let y = sample.beauty_xyz[1];
        let delta = y - next.mean_y;
        let mean_y = next.mean_y + delta / f64::from(next_count);
        let m2_y = delta.mul_add(y - mean_y, next.m2_y);
        if !mean_y.is_finite() || !m2_y.is_finite() || m2_y < 0.0 {
            return Err(CinematicAovError::NonFiniteChannel {
                channel: "variance.Y",
            });
        }
        next.mean_y = canonical_zero(mean_y);
        next.m2_y = canonical_zero(m2_y);
        next.accepted_count = next_count;

        if let Some(primary) = sample.primary {
            validate_primary(primary)?;
            next.primary_count = checked_increment(next.primary_count)?;
            add3(
                &mut next.shading_normal_sum,
                primary.shading_normal_world,
                "normal",
            )?;
            next.depth_sum_m = checked_add(next.depth_sum_m, primary.depth_m, "depth.Z")?;
            if primary.has_authored_shading_normal {
                next.authored_shading_normal_count =
                    checked_increment(next.authored_shading_normal_count)?;
            }
            if let Some(albedo) = primary.albedo_linear_rgb {
                validate_vec("albedo", albedo)?;
                add3(&mut next.albedo_sum, albedo, "albedo")?;
                next.albedo_count = checked_increment(next.albedo_count)?;
            }
            if let Some(motion) = primary.previous_motion_pixels {
                validate_vec("motion.prev", motion)?;
                add2(&mut next.previous_motion_sum_pixels, motion, "motion.prev")?;
                next.previous_motion_count = checked_increment(next.previous_motion_count)?;
            }
        }
        *self = next;
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct CategoricalPrimary {
    present: bool,
    distance_squared: f64,
    absolute_sample: u32,
    primitive_index: u64,
    object_palette_index: u32,
    material_palette_index: u32,
}

impl CategoricalPrimary {
    const NONE: Self = Self {
        present: false,
        distance_squared: 0.0,
        absolute_sample: 0,
        primitive_index: 0,
        object_palette_index: 0,
        material_palette_index: 0,
    };

    fn from_sample(sample: AlignedAovSample) -> Result<Self, CinematicAovError> {
        let Some(primary) = sample.primary else {
            return Ok(Self::NONE);
        };
        let dx = sample.pixel_jitter[0] - 0.5;
        let dy = sample.pixel_jitter[1] - 0.5;
        let distance_squared = dx.mul_add(dx, dy * dy);
        if !distance_squared.is_finite() || distance_squared < 0.0 {
            return Err(CinematicAovError::NonFiniteChannel {
                channel: "categorical sample rank",
            });
        }
        Ok(Self {
            present: true,
            distance_squared: canonical_zero(distance_squared),
            absolute_sample: sample.absolute_sample,
            primitive_index: u64::try_from(primary.primitive_index)
                .map_err(|_| CinematicAovError::InvalidPrimary)?,
            object_palette_index: primary.object_palette_index,
            material_palette_index: primary.material_palette_index,
        })
    }

    fn tie_rank(self) -> (u32, u64, u32, u32) {
        (
            self.absolute_sample,
            self.primitive_index,
            self.object_palette_index,
            self.material_palette_index,
        )
    }

    fn is_better_than(self, current: Self) -> bool {
        if !self.present {
            return false;
        }
        if !current.present {
            return true;
        }
        match self.distance_squared.total_cmp(&current.distance_squared) {
            core::cmp::Ordering::Less => true,
            core::cmp::Ordering::Greater => false,
            core::cmp::Ordering::Equal => self.tie_rank() < current.tie_rank(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct FinalPixel {
    geometric_normal_sum: [f64; 3],
    direct_xyz_sum: [f64; 3],
    indirect_xyz_sum: [f64; 3],
    emission_xyz_sum: [f64; 3],
    nearest_primary: CategoricalPrimary,
}

impl FinalPixel {
    const EMPTY: Self = Self {
        geometric_normal_sum: [0.0; 3],
        direct_xyz_sum: [0.0; 3],
        indirect_xyz_sum: [0.0; 3],
        emission_xyz_sum: [0.0; 3],
        nearest_primary: CategoricalPrimary::NONE,
    };

    fn push(&mut self, sample: AlignedAovSample) -> Result<(), CinematicAovError> {
        validate_vec("direct", sample.direct_xyz)?;
        validate_vec("indirect", sample.indirect_xyz)?;
        validate_vec("emission", sample.emission_xyz)?;
        let mut next = *self;
        add3(&mut next.direct_xyz_sum, sample.direct_xyz, "direct")?;
        add3(&mut next.indirect_xyz_sum, sample.indirect_xyz, "indirect")?;
        add3(&mut next.emission_xyz_sum, sample.emission_xyz, "emission")?;
        if let Some(primary) = sample.primary {
            if primary.material_palette_index == 0 {
                return Err(CinematicAovError::InvalidPrimary);
            }
            add3(
                &mut next.geometric_normal_sum,
                primary.geometric_normal_world,
                "normal_geom",
            )?;
        }
        let candidate = CategoricalPrimary::from_sample(sample)?;
        if candidate.is_better_than(next.nearest_primary) {
            next.nearest_primary = candidate;
        }
        *self = next;
        Ok(())
    }
}

/// Deterministic lossless mappings for `FLOAT` identity channels.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CinematicAovPalette {
    object_ids: Vec<u64>,
    material_identities: Vec<ContentHash>,
}

impl CinematicAovPalette {
    pub(crate) fn from_scene(
        scene: &Scene,
        limits: CinematicAovLimits,
        include_id_channels: bool,
    ) -> Result<Self, CinematicAovError> {
        if !include_id_channels {
            return Ok(Self {
                object_ids: Vec::new(),
                material_identities: Vec::new(),
            });
        }
        let primitive_count = scene.primitives.len();
        let mut object_ids = Vec::new();
        let mut material_identities = Vec::new();
        object_ids
            .try_reserve_exact(primitive_count)
            .map_err(|_| CinematicAovError::AllocationRefused)?;
        material_identities
            .try_reserve_exact(primitive_count)
            .map_err(|_| CinematicAovError::AllocationRefused)?;
        for primitive in &scene.primitives {
            match &primitive.shape {
                Shape::Instance(instance) => object_ids.push(instance.object_id()),
                Shape::AnimatedInstance(instance) => object_ids.push(instance.object_id()),
                Shape::Mesh(_) | Shape::Chart(_) => {}
            }
            material_identities.push(primitive.material.content_identity());
        }
        object_ids.sort_unstable();
        object_ids.dedup();
        material_identities.sort_unstable();
        material_identities.dedup();
        for (kind, count) in [
            ("object", object_ids.len()),
            ("material", material_identities.len()),
        ] {
            let count = u32::try_from(count).map_err(|_| CinematicAovError::PaletteLimit {
                kind,
                requested: u64::MAX,
                limit: u64::from(limits.max_palette_entries),
            })?;
            if count > limits.max_palette_entries {
                return Err(CinematicAovError::PaletteLimit {
                    kind,
                    requested: u64::from(count),
                    limit: u64::from(limits.max_palette_entries),
                });
            }
        }
        Ok(Self {
            object_ids,
            material_identities,
        })
    }

    pub(crate) fn object_index(&self, object_id: u64) -> Result<u32, CinematicAovError> {
        let zero_based = self
            .object_ids
            .binary_search(&object_id)
            .map_err(|_| CinematicAovError::PaletteMismatch)?;
        palette_index(zero_based)
    }

    pub(crate) fn material_index(&self, identity: ContentHash) -> Result<u32, CinematicAovError> {
        let zero_based = self
            .material_identities
            .binary_search(&identity)
            .map_err(|_| CinematicAovError::PaletteMismatch)?;
        palette_index(zero_based)
    }

    /// Sorted nonzero object IDs. EXR index `i + 1` maps to entry `i`.
    #[must_use]
    pub fn object_ids(&self) -> &[u64] {
        &self.object_ids
    }

    /// Sorted material hashes. EXR index `i + 1` maps to entry `i`.
    #[must_use]
    pub fn material_identities(&self) -> &[ContentHash] {
        &self.material_identities
    }

    fn try_clone(&self) -> Result<Self, CinematicAovError> {
        Ok(Self {
            object_ids: fallible_copy(&self.object_ids)?,
            material_identities: fallible_copy(&self.material_identities)?,
        })
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct CinematicAovRenderBinding {
    pub settings: Settings,
    pub shutter: ShutterInterval,
    pub shot_id: u64,
    pub palette: CinematicAovPalette,
    /// Process-local guard against progressive mutation of borrowed scene,
    /// lighting, geometry, or camera values. This is not exported authority.
    pub continuity_fingerprint: ContentHash,
    pub adaptive_policy: Option<AdaptiveSamplingConfig>,
}

impl CinematicAovRenderBinding {
    fn try_clone(&self) -> Result<Self, CinematicAovError> {
        Ok(Self {
            settings: self.settings,
            shutter: self.shutter,
            shot_id: self.shot_id,
            palette: self.palette.try_clone()?,
            continuity_fingerprint: self.continuity_fingerprint,
            adaptive_policy: self.adaptive_policy,
        })
    }
}

/// Raw beauty plus aligned, unnormalised AOV accumulator state.
///
/// This is an estimate artifact. It does not certify physical material values,
/// trajectory truth, target-frame visibility, or Euler-disc mechanics.
#[derive(Debug, PartialEq)]
pub struct CinematicAovFilm {
    beauty: Film,
    config: CinematicAovConfig,
    common: Option<Vec<CommonPixel>>,
    final_diagnostic: Option<Vec<FinalPixel>>,
    binding: Option<CinematicAovRenderBinding>,
    retained_bytes: u64,
}

impl CinematicAovFilm {
    /// Allocate an empty shape-aligned film under the declared retained-memory
    /// limit. No render binding is committed until a nonempty range succeeds.
    pub fn try_new(
        width: u32,
        height: u32,
        config: CinematicAovConfig,
    ) -> Result<Self, CinematicAovError> {
        let pixel_count = checked_pixel_count(width, height)?;
        let pixel_count_u64 = u64::try_from(pixel_count)
            .map_err(|_| CinematicAovError::InvalidDimensions { width, height })?;
        if pixel_count_u64 > config.limits.max_pixels {
            return Err(CinematicAovError::PixelLimit {
                requested: pixel_count_u64,
                limit: config.limits.max_pixels,
            });
        }
        let retained_bytes = retained_bytes(pixel_count, config.profile)?;
        if retained_bytes > config.limits.max_retained_bytes {
            return Err(CinematicAovError::RetainedMemoryLimit {
                requested: retained_bytes,
                limit: config.limits.max_retained_bytes,
            });
        }
        let common = config
            .profile
            .has_common()
            .then(|| fallible_filled(pixel_count, CommonPixel::EMPTY))
            .transpose()?;
        let final_diagnostic = config
            .profile
            .has_final()
            .then(|| fallible_filled(pixel_count, FinalPixel::EMPTY))
            .transpose()?;
        Ok(Self {
            beauty: Film::try_new(width, height).map_err(CinematicAovError::Tracer)?,
            config,
            common,
            final_diagnostic,
            binding: None,
            retained_bytes,
        })
    }

    /// Raw, unmodified beauty film.
    #[must_use]
    pub const fn beauty(&self) -> &Film {
        &self.beauty
    }

    /// Complete AOV configuration.
    #[must_use]
    pub const fn config(&self) -> CinematicAovConfig {
        self.config
    }

    /// Exact retained payload bytes charged at admission.
    #[must_use]
    pub const fn retained_bytes(&self) -> u64 {
        self.retained_bytes
    }

    /// Scene-derived palette after the first nonempty successful range.
    #[must_use]
    pub fn palette(&self) -> Option<&CinematicAovPalette> {
        self.binding.as_ref().map(|binding| &binding.palette)
    }

    /// Per-pixel accepted sample count, when the selected profile has AOVs.
    #[must_use]
    pub fn sample_count(&self, pixel: usize) -> Option<u32> {
        self.common
            .as_ref()
            .and_then(|plane| plane.get(pixel))
            .map(|sample| sample.accepted_count)
    }

    /// Per-pixel primary-hit count, when the selected profile has AOVs.
    #[must_use]
    pub fn primary_count(&self, pixel: usize) -> Option<u32> {
        self.common
            .as_ref()
            .and_then(|plane| plane.get(pixel))
            .map(|sample| sample.primary_count)
    }

    pub(crate) fn beauty_mut(&mut self) -> &mut Film {
        &mut self.beauty
    }

    pub(crate) fn binding(&self) -> Option<&CinematicAovRenderBinding> {
        self.binding.as_ref()
    }

    pub(crate) fn bind(&mut self, binding: CinematicAovRenderBinding) {
        self.binding = Some(binding);
    }

    pub(crate) fn push(
        &mut self,
        pixel: usize,
        sample: AlignedAovSample,
    ) -> Result<(), CinematicAovError> {
        if self.config.profile.has_final()
            && self
                .common
                .as_ref()
                .and_then(|plane| plane.get(pixel))
                .is_some_and(|state| state.accepted_count == MAX_EXACT_F32_INTEGER)
        {
            return Err(CinematicAovError::InexactSampleCount {
                samples: MAX_EXACT_F32_INTEGER.saturating_add(1),
            });
        }
        if let Some(common) = &mut self.common {
            let target = common
                .get_mut(pixel)
                .ok_or(CinematicAovError::ShapeMismatch)?;
            target.push(sample)?;
        }
        if let Some(final_diagnostic) = &mut self.final_diagnostic {
            let target = final_diagnostic
                .get_mut(pixel)
                .ok_or(CinematicAovError::ShapeMismatch)?;
            target.push(sample)?;
        }
        Ok(())
    }

    /// Fallibly clone all accumulator state for an all-or-nothing progressive
    /// append. Publication remains one final assignment in the caller.
    pub(crate) fn try_clone_for_stage(&self) -> Result<Self, CinematicAovError> {
        let staging_peak_bytes = self
            .retained_bytes
            .checked_mul(2)
            .ok_or(CinematicAovError::SizeOverflow)?;
        if staging_peak_bytes > self.config.limits.max_retained_bytes {
            return Err(CinematicAovError::RetainedMemoryLimit {
                requested: staging_peak_bytes,
                limit: self.config.limits.max_retained_bytes,
            });
        }
        Ok(Self {
            beauty: Film {
                width: self.beauty.width,
                height: self.beauty.height,
                xyz: fallible_copy(&self.beauty.xyz)?,
                spp_done: self.beauty.spp_done,
                time_mode: self.beauty.time_mode,
            },
            config: self.config,
            common: self.common.as_deref().map(fallible_copy).transpose()?,
            final_diagnostic: self
                .final_diagnostic
                .as_deref()
                .map(fallible_copy)
                .transpose()?,
            binding: self
                .binding
                .as_ref()
                .map(|value| value.try_clone())
                .transpose()?,
            retained_bytes: self.retained_bytes,
        })
    }

    /// Encode the frozen profile as deterministic linear-sRGB `FLOAT` EXR
    /// channels plus lossless authority, frame, sampler, and palette metadata.
    pub fn to_exr(&self) -> Result<Vec<u8>, CinematicAovError> {
        let binding = self
            .binding
            .as_ref()
            .ok_or(CinematicAovError::UnboundFilm)?;
        self.validate_complete()?;
        let pixel_count = self.beauty.xyz.len();
        let channel_count = u64::from(self.config.profile.float_channel_count());
        let plane_bytes = u64::try_from(pixel_count)
            .ok()
            .and_then(|pixels| pixels.checked_mul(channel_count))
            .and_then(|samples| samples.checked_mul(size_of::<f32>() as u64))
            .ok_or(CinematicAovError::SizeOverflow)?;
        if plane_bytes > self.config.limits.max_export_plane_bytes {
            return Err(CinematicAovError::ExportMemoryLimit {
                requested: plane_bytes,
                limit: self.config.limits.max_export_plane_bytes,
            });
        }

        let spp = f64::from(self.beauty.spp_done);
        let channels = build_exr_channels(
            self.config.profile,
            pixel_count,
            |pixel| self.beauty.xyz[pixel].map(|value| value / spp),
            |_| self.beauty.spp_done,
            self.common.as_deref(),
            self.final_diagnostic.as_deref(),
        )?;
        let attributes = self.exr_attributes(binding)?;
        fs_img::write_exr_with_attributes(
            self.beauty.width,
            self.beauty.height,
            &channels,
            &attributes,
        )
        .map_err(CinematicAovError::Image)
    }

    fn validate_complete(&self) -> Result<(), CinematicAovError> {
        if self.beauty.spp_done == 0
            || (self.config.profile.has_final() && self.beauty.spp_done > MAX_EXACT_F32_INTEGER)
        {
            return Err(CinematicAovError::InexactSampleCount {
                samples: self.beauty.spp_done,
            });
        }
        let expected = self.beauty.xyz.len();
        if self
            .common
            .as_ref()
            .is_some_and(|plane| plane.len() != expected)
            || self
                .final_diagnostic
                .as_ref()
                .is_some_and(|plane| plane.len() != expected)
        {
            return Err(CinematicAovError::ShapeMismatch);
        }
        if let Some(common) = &self.common
            && common
                .iter()
                .any(|pixel| pixel.accepted_count != self.beauty.spp_done)
        {
            return Err(CinematicAovError::SampleAlignmentMismatch);
        }
        Ok(())
    }

    fn exr_attributes(
        &self,
        binding: &CinematicAovRenderBinding,
    ) -> Result<Vec<ExrAttribute>, CinematicAovError> {
        cinematic_exr_attributes(
            self.config,
            binding,
            "uniform",
            self.beauty.spp_done.to_string(),
        )
    }
}

fn cinematic_exr_attributes(
    config: CinematicAovConfig,
    binding: &CinematicAovRenderBinding,
    sample_mode: &'static str,
    rendered_spp: String,
) -> Result<Vec<ExrAttribute>, CinematicAovError> {
    let provenance = config.provenance;
    let mut attributes = Vec::new();
    attributes
        .try_reserve_exact(31)
        .map_err(|_| CinematicAovError::AllocationRefused)?;
    let mut push = |name: &str, value: String| {
        attributes.push(ExrAttribute {
            name: name.to_string(),
            ty: STRING_ATTRIBUTE_TYPE.to_string(),
            value: value.into_bytes(),
        });
    };
    push("frankensim.aov.authority", "raw-estimate".to_string());
    push(
        "frankensim.aov.schemaVersion",
        CINEMATIC_AOV_SEMANTICS_VERSION.to_string(),
    );
    push("frankensim.aov.profile", config.profile.name().to_string());
    push("frankensim.aov.configHash", config.identity.to_hex());
    push(
        "frankensim.aov.channelSemantics",
        "R,G,B=linear-sRGB;albedo=linear-sRGB-reflectance;normal,normal_geom=world-unit-vector;depth.Z=axial-metres;primary.coverage=primary-hit-fraction;variance.Y=unbiased-raw-CIE-Y-sample-variance;motion.prev=target-minus-current-raster-pixels;IDs=exact-palette-index;samples=count;invalid=zero"
            .to_string(),
    );
    push(
        "frankensim.aov.invalidSemantics",
        "zero-with-diagnostic.validity-bitmask-final-profile;daily-core-surface-validity-is-primary.coverage-greater-than-zero"
            .to_string(),
    );
    push(
        "frankensim.aov.materialDomain",
        MATERIAL_CONTENT_IDENTITY_DOMAIN.to_string(),
    );
    push("frankensim.frame.index", provenance.frame_index.to_string());
    push(
        "frankensim.frame.timeSeconds",
        f64_bits_string(provenance.frame_time_s),
    );
    push(
        "frankensim.frame.previousTimeS",
        f64_bits_string(provenance.previous_frame_time_s),
    );
    push(
        "frankensim.frame.nextTimeS",
        f64_bits_string(provenance.next_frame_time_s),
    );
    push(
        "frankensim.source.trajectory",
        provenance.source_trajectory_identity.to_hex(),
    );
    push(
        "frankensim.source.sceneHash",
        provenance.scene_identity.to_hex(),
    );
    push(
        "frankensim.source.composition",
        provenance.composition_identity.to_hex(),
    );
    push("frankensim.render.seed", binding.settings.seed.to_string());
    push(
        "frankensim.render.sampler",
        sampler_name(binding.settings.sampler).to_string(),
    );
    push(
        "frankensim.render.strategy",
        direct_strategy_name(binding.settings.strategy).to_string(),
    );
    push(
        "frankensim.render.maxDepth",
        binding.settings.max_depth.to_string(),
    );
    push("frankensim.render.sampleMode", sample_mode.to_string());
    push("frankensim.render.spp", rendered_spp);
    push(
        "frankensim.render.sppCeiling",
        binding.settings.spp.to_string(),
    );
    if let Some(policy) = binding.adaptive_policy {
        push(
            // OpenEXR attribute names are limited to 31 UTF-8 bytes.
            "frankensim.render.adaptive",
            format!(
                "version={ADAPTIVE_SAMPLING_SEMANTICS_VERSION};minimum={};batch={};absolute={};relative={};darkFloor={}",
                policy.minimum_samples(),
                policy.batch_samples(),
                f64_bits_string(policy.absolute_error()),
                f64_bits_string(policy.relative_error()),
                f64_bits_string(policy.dark_floor())
            ),
        );
    }
    push("frankensim.render.shotId", binding.shot_id.to_string());
    push(
        "frankensim.render.shutterOpenS",
        f64_bits_string(binding.shutter.open_s()),
    );
    push(
        "frankensim.render.shutterCloseS",
        f64_bits_string(binding.shutter.close_s()),
    );
    push(
        "frankensim.render.versions",
        format!(
            "tracer={TRACER_BIT_SEMANTICS_VERSION};motionTracer={MOTION_TRACER_BIT_SEMANTICS_VERSION};cinematicCamera={CINEMATIC_CAMERA_TRACER_BIT_SEMANTICS_VERSION};dielectric={DIELECTRIC_TRACER_BIT_SEMANTICS_VERSION};lighting={LIGHTING_TRACER_BIT_SEMANTICS_VERSION};motionVector={MOTION_VECTOR_SEMANTICS_VERSION};category={CINEMATIC_AOV_CATEGORY_SEMANTICS_VERSION};albedo={CINEMATIC_AOV_ALBEDO_SEMANTICS_VERSION}"
        ),
    );
    push(
        "frankensim.aov.objectPalette",
        encode_object_palette(&binding.palette.object_ids),
    );
    push(
        "frankensim.aov.materialPalette",
        encode_material_palette(&binding.palette.material_identities),
    );
    push(
        "frankensim.aov.paletteZero",
        "0=background-or-unavailable".to_string(),
    );
    Ok(attributes)
}

/// Raw adaptive beauty plus aligned cinematic denoising and diagnostic AOVs.
///
/// Every pixel owns its exact terminal sample count from [`AdaptiveFilm`]. AOV
/// samples never participate in the stopping decision; they merely observe the
/// same accepted path prefix. This remains a raw estimate artifact, not an
/// image-error certificate or a claim about physical scene truth.
#[derive(Debug, PartialEq)]
pub struct AdaptiveCinematicAovFilm {
    pub(crate) beauty: AdaptiveFilm,
    pub(crate) config: CinematicAovConfig,
    common: Option<Vec<CommonPixel>>,
    final_diagnostic: Option<Vec<FinalPixel>>,
    pub(crate) binding: CinematicAovRenderBinding,
    retained_bytes: u64,
}

impl AdaptiveCinematicAovFilm {
    /// Exact adaptive beauty estimator and stopping state.
    #[must_use]
    pub const fn beauty(&self) -> &AdaptiveFilm {
        &self.beauty
    }

    /// Complete AOV configuration.
    #[must_use]
    pub const fn config(&self) -> CinematicAovConfig {
        self.config
    }

    /// Scene-derived palette. It is empty for profiles without ID channels.
    #[must_use]
    pub const fn palette(&self) -> &CinematicAovPalette {
        &self.binding.palette
    }

    /// Exact retained payload bytes charged at admission. Allocator overhead,
    /// metadata strings, EXR encoder scratch, and encoded output are excluded.
    #[must_use]
    pub const fn retained_bytes(&self) -> u64 {
        self.retained_bytes
    }

    /// Encode the adaptive raw estimate with the same frozen AOV schema as the
    /// uniform film. `samples` carries each pixel's actual terminal count.
    pub fn to_exr(&self) -> Result<Vec<u8>, CinematicAovError> {
        self.validate_complete()?;
        let pixel_count = self.beauty.xyz_sums().len();
        validate_export_plane_budget(pixel_count, self.config)?;
        let channels = build_exr_channels(
            self.config.profile,
            pixel_count,
            |pixel| {
                self.beauty
                    .beauty_mean_xyz(pixel)
                    .expect("validated adaptive film owns shape-matched buffers")
            },
            |pixel| self.beauty.sample_counts()[pixel],
            self.common.as_deref(),
            self.final_diagnostic.as_deref(),
        )?;
        let attributes = cinematic_exr_attributes(
            self.config,
            &self.binding,
            "adaptive",
            "per-pixel-channel".to_string(),
        )?;
        fs_img::write_exr_with_attributes(
            self.beauty.width(),
            self.beauty.height(),
            &channels,
            &attributes,
        )
        .map_err(CinematicAovError::Image)
    }

    fn validate_complete(&self) -> Result<(), CinematicAovError> {
        let expected = checked_pixel_count(self.beauty.width(), self.beauty.height())?;
        if self.beauty.xyz_sums().len() != expected
            || self.beauty.running_means_xyz().len() != expected
            || self.beauty.m2_xyz().len() != expected
            || self.beauty.sample_counts().len() != expected
            || self.beauty.decisions().len() != expected
            || self
                .common
                .as_ref()
                .is_some_and(|plane| plane.len() != expected)
            || self
                .final_diagnostic
                .as_ref()
                .is_some_and(|plane| plane.len() != expected)
        {
            return Err(CinematicAovError::ShapeMismatch);
        }
        if self.binding.adaptive_policy != Some(self.beauty.policy())
            || self.binding.settings.width != self.beauty.width()
            || self.binding.settings.height != self.beauty.height()
            || self.binding.settings.spp != self.beauty.maximum_samples()
            || self.binding.settings.sampler != self.beauty.sampler()
            || self.binding.settings.seed != self.beauty.stream_seed()
            || self.beauty.semantics_version() != ADAPTIVE_SAMPLING_SEMANTICS_VERSION
        {
            return Err(CinematicAovError::ProgressiveBindingMismatch);
        }
        if self.beauty.sample_counts().contains(&0)
            || (self.config.profile.has_final()
                && self
                    .beauty
                    .sample_counts()
                    .iter()
                    .any(|samples| *samples > MAX_EXACT_F32_INTEGER))
        {
            return Err(CinematicAovError::SampleAlignmentMismatch);
        }
        if let Some(common) = &self.common {
            for (pixel, state) in common.iter().enumerate() {
                if state.accepted_count != self.beauty.sample_counts()[pixel]
                    || state.mean_y.to_bits() != self.beauty.running_means_xyz()[pixel][1].to_bits()
                    || state.m2_y.to_bits() != self.beauty.m2_xyz()[pixel][1].to_bits()
                {
                    return Err(CinematicAovError::SampleAlignmentMismatch);
                }
            }
        }
        Ok(())
    }
}

/// Private per-pixel AOV construction state used while adaptive beauty remains
/// unpublished. It cannot influence adaptive stopping decisions.
pub(crate) struct AdaptiveAovAccumulator {
    width: u32,
    height: u32,
    config: CinematicAovConfig,
    common: Option<Vec<CommonPixel>>,
    final_diagnostic: Option<Vec<FinalPixel>>,
    retained_bytes: u64,
}

impl AdaptiveAovAccumulator {
    pub(crate) fn try_new(
        width: u32,
        height: u32,
        maximum_samples: u32,
        config: CinematicAovConfig,
    ) -> Result<Self, CinematicAovError> {
        let pixel_count = checked_pixel_count(width, height)?;
        let requested_pixels = u64::try_from(pixel_count)
            .map_err(|_| CinematicAovError::InvalidDimensions { width, height })?;
        if requested_pixels > config.limits.max_pixels {
            return Err(CinematicAovError::PixelLimit {
                requested: requested_pixels,
                limit: config.limits.max_pixels,
            });
        }
        if config.profile.has_final() && maximum_samples > MAX_EXACT_F32_INTEGER {
            return Err(CinematicAovError::InexactSampleCount {
                samples: maximum_samples,
            });
        }
        let retained_bytes = adaptive_retained_bytes(pixel_count, config.profile)?;
        if retained_bytes > config.limits.max_retained_bytes {
            return Err(CinematicAovError::RetainedMemoryLimit {
                requested: retained_bytes,
                limit: config.limits.max_retained_bytes,
            });
        }
        let common = config
            .profile
            .has_common()
            .then(|| fallible_filled(pixel_count, CommonPixel::EMPTY))
            .transpose()?;
        let final_diagnostic = config
            .profile
            .has_final()
            .then(|| fallible_filled(pixel_count, FinalPixel::EMPTY))
            .transpose()?;
        Ok(Self {
            width,
            height,
            config,
            common,
            final_diagnostic,
            retained_bytes,
        })
    }

    pub(crate) fn push(
        &mut self,
        pixel: usize,
        sample: AlignedAovSample,
    ) -> Result<(), CinematicAovError> {
        if let Some(common) = &mut self.common {
            common
                .get_mut(pixel)
                .ok_or(CinematicAovError::ShapeMismatch)?
                .push(sample)?;
        }
        if let Some(final_diagnostic) = &mut self.final_diagnostic {
            final_diagnostic
                .get_mut(pixel)
                .ok_or(CinematicAovError::ShapeMismatch)?
                .push(sample)?;
        }
        Ok(())
    }

    pub(crate) fn publish(
        self,
        beauty: AdaptiveFilm,
        binding: CinematicAovRenderBinding,
    ) -> Result<AdaptiveCinematicAovFilm, CinematicAovError> {
        if (beauty.width(), beauty.height()) != (self.width, self.height) {
            return Err(CinematicAovError::ShapeMismatch);
        }
        let film = AdaptiveCinematicAovFilm {
            beauty,
            config: self.config,
            common: self.common,
            final_diagnostic: self.final_diagnostic,
            binding,
            retained_bytes: self.retained_bytes,
        };
        film.validate_complete()?;
        Ok(film)
    }
}

/// Fail-closed AOV admission, accumulation, and artifact errors.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CinematicAovError {
    /// Underlying tracer/camera/geometry path refused.
    Tracer(crate::tracer::TracerError),
    /// Adaptive policy or beauty-estimator state refused.
    Adaptive(crate::tracer::AdaptiveSamplingError),
    /// EXR writer refused a shape or metadata value.
    Image(fs_img::ImgError),
    /// One or more resource limits were zero or could not encode exact indices.
    InvalidLimits,
    /// Raster dimensions were zero or overflowed.
    InvalidDimensions {
        /// Rejected raster width.
        width: u32,
        /// Rejected raster height.
        height: u32,
    },
    /// Frame reference times were non-finite or unordered.
    InvalidFrameReferences,
    /// Required external identity was the all-zero sentinel.
    MissingIdentity {
        /// Name of the absent identity field.
        field: &'static str,
    },
    /// Raster exceeded the admitted pixel count.
    PixelLimit {
        /// Requested raster pixel count.
        requested: u64,
        /// Configured pixel ceiling.
        limit: u64,
    },
    /// Retained beauty/AOV state exceeded the admitted byte count.
    RetainedMemoryLimit {
        /// Required owned accumulator bytes, including a staging copy when one
        /// is about to be created.
        requested: u64,
        /// Configured owned-accumulator byte ceiling.
        limit: u64,
    },
    /// Staged EXR planes exceeded the admitted byte count.
    ExportMemoryLimit {
        /// Required uncompressed channel-plane bytes.
        requested: u64,
        /// Configured uncompressed channel-plane byte ceiling.
        limit: u64,
    },
    /// Scene palette exceeded its entry limit.
    PaletteLimit {
        /// Palette family (`object` or `material`).
        kind: &'static str,
        /// Required number of nonzero entries.
        requested: u64,
        /// Configured entry ceiling.
        limit: u64,
    },
    /// Scene and retained palette disagreed.
    PaletteMismatch,
    /// Progressive append changed settings, shutter, shot, profile, or palette.
    ProgressiveBindingMismatch,
    /// Frame reference times do not enclose the complete beauty shutter.
    ReferenceTimesDoNotCoverShutter,
    /// A public or restored plane shape disagreed with the raster.
    ShapeMismatch,
    /// AOV sample count diverged from the aligned beauty prefix.
    SampleAlignmentMismatch,
    /// A channel input or accumulated value became non-finite.
    NonFiniteChannel {
        /// Channel or accumulator family that became invalid.
        channel: &'static str,
    },
    /// Primary AOV geometry or palette data was invalid.
    InvalidPrimary,
    /// Per-pixel sample count overflowed.
    SampleCountOverflow,
    /// Sample count cannot be represented exactly by the EXR float channel.
    InexactSampleCount {
        /// Count that cannot be represented exactly in a `FLOAT` channel.
        samples: u32,
    },
    /// No successful nonempty render range has bound the film.
    UnboundFilm,
    /// A fallible vector reservation failed.
    AllocationRefused,
    /// Checked byte or element arithmetic overflowed.
    SizeOverflow,
    /// Internal profile/channel definition disagreed with its frozen count.
    InternalChannelCount,
}

impl fmt::Display for CinematicAovError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Tracer(error) => write!(formatter, "cinematic AOV trace refused: {error}"),
            Self::Adaptive(error) => {
                write!(formatter, "cinematic adaptive AOV refused: {error}")
            }
            Self::Image(error) => write!(formatter, "cinematic AOV EXR refused: {error}"),
            Self::InvalidLimits => formatter.write_str("cinematic AOV limits are invalid"),
            Self::InvalidDimensions { width, height } => {
                write!(
                    formatter,
                    "cinematic AOV raster {width}x{height} is invalid"
                )
            }
            Self::InvalidFrameReferences => {
                formatter.write_str("cinematic AOV frame references must be finite and ordered")
            }
            Self::MissingIdentity { field } => {
                write!(formatter, "cinematic AOV provenance is missing {field}")
            }
            Self::PixelLimit { requested, limit } => {
                write!(
                    formatter,
                    "cinematic AOV raster needs {requested} pixels above limit {limit}"
                )
            }
            Self::RetainedMemoryLimit { requested, limit } => write!(
                formatter,
                "cinematic AOV retained state needs {requested} bytes above limit {limit}"
            ),
            Self::ExportMemoryLimit { requested, limit } => write!(
                formatter,
                "cinematic AOV EXR planes need {requested} bytes above limit {limit}"
            ),
            Self::PaletteLimit {
                kind,
                requested,
                limit,
            } => write!(
                formatter,
                "cinematic AOV {kind} palette needs {requested} entries above limit {limit}"
            ),
            Self::PaletteMismatch => {
                formatter.write_str("cinematic AOV scene identity palette mismatch")
            }
            Self::ProgressiveBindingMismatch => formatter.write_str(
                "cinematic AOV progressive append changed its render or palette binding",
            ),
            Self::ReferenceTimesDoNotCoverShutter => formatter.write_str(
                "cinematic AOV previous/next reference times do not cover the beauty shutter",
            ),
            Self::ShapeMismatch => formatter.write_str("cinematic AOV plane shape mismatch"),
            Self::SampleAlignmentMismatch => {
                formatter.write_str("cinematic AOV and beauty sample prefixes are not aligned")
            }
            Self::NonFiniteChannel { channel } => {
                write!(
                    formatter,
                    "cinematic AOV channel {channel} became non-finite"
                )
            }
            Self::InvalidPrimary => formatter.write_str("cinematic AOV primary sample is invalid"),
            Self::SampleCountOverflow => {
                formatter.write_str("cinematic AOV sample count overflowed")
            }
            Self::InexactSampleCount { samples } => write!(
                formatter,
                "cinematic AOV sample count {samples} is not a positive exact EXR FLOAT integer"
            ),
            Self::UnboundFilm => {
                formatter.write_str("cinematic AOV film has no successful render binding")
            }
            Self::AllocationRefused => formatter.write_str("cinematic AOV allocation was refused"),
            Self::SizeOverflow => formatter.write_str("cinematic AOV size arithmetic overflowed"),
            Self::InternalChannelCount => {
                formatter.write_str("cinematic AOV profile channel count is inconsistent")
            }
        }
    }
}

impl core::error::Error for CinematicAovError {}

impl From<crate::tracer::TracerError> for CinematicAovError {
    fn from(error: crate::tracer::TracerError) -> Self {
        Self::Tracer(error)
    }
}

impl From<crate::tracer::AdaptiveSamplingError> for CinematicAovError {
    fn from(error: crate::tracer::AdaptiveSamplingError) -> Self {
        Self::Adaptive(error)
    }
}

pub(crate) fn validate_reference_times(
    config: CinematicAovConfig,
    shutter: ShutterInterval,
) -> Result<(), CinematicAovError> {
    let provenance = config.provenance;
    if provenance.previous_frame_time_s > shutter.open_s()
        || provenance.frame_time_s < shutter.open_s()
        || provenance.frame_time_s > shutter.close_s()
        || provenance.next_frame_time_s < shutter.close_s()
    {
        Err(CinematicAovError::ReferenceTimesDoNotCoverShutter)
    } else {
        Ok(())
    }
}

pub(crate) fn validate_binding(
    film: &CinematicAovFilm,
    settings: Settings,
    shutter: ShutterInterval,
    shot_id: u64,
    palette: &CinematicAovPalette,
    continuity_fingerprint: ContentHash,
) -> Result<(), CinematicAovError> {
    validate_reference_times(film.config, shutter)?;
    let candidate = CinematicAovRenderBinding {
        settings,
        shutter,
        shot_id,
        palette: palette.try_clone()?,
        continuity_fingerprint,
        adaptive_policy: None,
    };
    if let Some(binding) = &film.binding
        && binding != &candidate
    {
        return Err(CinematicAovError::ProgressiveBindingMismatch);
    }
    Ok(())
}

pub(crate) fn render_binding(
    settings: Settings,
    shutter: ShutterInterval,
    shot_id: u64,
    palette: CinematicAovPalette,
    continuity_fingerprint: ContentHash,
) -> CinematicAovRenderBinding {
    CinematicAovRenderBinding {
        settings,
        shutter,
        shot_id,
        palette,
        continuity_fingerprint,
        adaptive_policy: None,
    }
}

pub(crate) fn adaptive_render_binding(
    settings: Settings,
    shutter: ShutterInterval,
    shot_id: u64,
    palette: CinematicAovPalette,
    continuity_fingerprint: ContentHash,
    adaptive_policy: AdaptiveSamplingConfig,
) -> CinematicAovRenderBinding {
    CinematicAovRenderBinding {
        settings,
        shutter,
        shot_id,
        palette,
        continuity_fingerprint,
        adaptive_policy: Some(adaptive_policy),
    }
}

fn checked_pixel_count(width: u32, height: u32) -> Result<usize, CinematicAovError> {
    let pixels = width
        .checked_mul(height)
        .filter(|pixels| *pixels != 0)
        .ok_or(CinematicAovError::InvalidDimensions { width, height })?;
    usize::try_from(pixels).map_err(|_| CinematicAovError::InvalidDimensions { width, height })
}

fn retained_bytes(
    pixel_count: usize,
    profile: CinematicAovProfile,
) -> Result<u64, CinematicAovError> {
    let bytes_per_pixel = size_of::<[f64; 3]>()
        .checked_add(if profile.has_common() {
            size_of::<CommonPixel>()
        } else {
            0
        })
        .and_then(|bytes| {
            bytes.checked_add(if profile.has_final() {
                size_of::<FinalPixel>()
            } else {
                0
            })
        })
        .ok_or(CinematicAovError::SizeOverflow)?;
    u64::try_from(pixel_count)
        .ok()
        .and_then(|pixels| pixels.checked_mul(bytes_per_pixel as u64))
        .ok_or(CinematicAovError::SizeOverflow)
}

fn adaptive_retained_bytes(
    pixel_count: usize,
    profile: CinematicAovProfile,
) -> Result<u64, CinematicAovError> {
    let beauty_bytes_per_pixel = size_of::<[f64; 3]>()
        .checked_mul(3)
        .and_then(|bytes| bytes.checked_add(size_of::<u32>()))
        .and_then(|bytes| bytes.checked_add(size_of::<AdaptiveDecision>()))
        .ok_or(CinematicAovError::SizeOverflow)?;
    let bytes_per_pixel = beauty_bytes_per_pixel
        .checked_add(if profile.has_common() {
            size_of::<CommonPixel>()
        } else {
            0
        })
        .and_then(|bytes| {
            bytes.checked_add(if profile.has_final() {
                size_of::<FinalPixel>()
            } else {
                0
            })
        })
        .ok_or(CinematicAovError::SizeOverflow)?;
    u64::try_from(pixel_count)
        .ok()
        .and_then(|pixels| pixels.checked_mul(bytes_per_pixel as u64))
        .ok_or(CinematicAovError::SizeOverflow)
}

fn validate_export_plane_budget(
    pixel_count: usize,
    config: CinematicAovConfig,
) -> Result<u64, CinematicAovError> {
    let channel_count = u64::from(config.profile.float_channel_count());
    let plane_bytes = u64::try_from(pixel_count)
        .ok()
        .and_then(|pixels| pixels.checked_mul(channel_count))
        .and_then(|samples| samples.checked_mul(size_of::<f32>() as u64))
        .ok_or(CinematicAovError::SizeOverflow)?;
    if plane_bytes > config.limits.max_export_plane_bytes {
        return Err(CinematicAovError::ExportMemoryLimit {
            requested: plane_bytes,
            limit: config.limits.max_export_plane_bytes,
        });
    }
    Ok(plane_bytes)
}

fn fallible_filled<T: Clone>(len: usize, value: T) -> Result<Vec<T>, CinematicAovError> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(len)
        .map_err(|_| CinematicAovError::AllocationRefused)?;
    output.resize(len, value);
    Ok(output)
}

fn fallible_copy<T: Copy>(source: &[T]) -> Result<Vec<T>, CinematicAovError> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(source.len())
        .map_err(|_| CinematicAovError::AllocationRefused)?;
    output.extend_from_slice(source);
    Ok(output)
}

fn palette_index(zero_based: usize) -> Result<u32, CinematicAovError> {
    let one_based = zero_based
        .checked_add(1)
        .ok_or(CinematicAovError::SizeOverflow)?;
    let index = u32::try_from(one_based).map_err(|_| CinematicAovError::SizeOverflow)?;
    if index >= MAX_EXACT_F32_INTEGER {
        return Err(CinematicAovError::SizeOverflow);
    }
    Ok(index)
}

fn checked_increment(value: u32) -> Result<u32, CinematicAovError> {
    value
        .checked_add(1)
        .ok_or(CinematicAovError::SampleCountOverflow)
}

fn validate_primary(primary: AlignedAovPrimary) -> Result<(), CinematicAovError> {
    validate_vec("normal_geom", primary.geometric_normal_world)?;
    validate_vec("normal", primary.shading_normal_world)?;
    if !primary.depth_m.is_finite() || primary.depth_m <= 0.0 {
        return Err(CinematicAovError::InvalidPrimary);
    }
    Ok(())
}

fn validate_vec<const N: usize>(
    channel: &'static str,
    values: [f64; N],
) -> Result<(), CinematicAovError> {
    if values.iter().any(|value| !value.is_finite()) {
        Err(CinematicAovError::NonFiniteChannel { channel })
    } else {
        Ok(())
    }
}

fn add3(
    target: &mut [f64; 3],
    values: [f64; 3],
    channel: &'static str,
) -> Result<(), CinematicAovError> {
    for index in 0..3 {
        target[index] = checked_add(target[index], values[index], channel)?;
    }
    Ok(())
}

fn add2(
    target: &mut [f64; 2],
    values: [f64; 2],
    channel: &'static str,
) -> Result<(), CinematicAovError> {
    for index in 0..2 {
        target[index] = checked_add(target[index], values[index], channel)?;
    }
    Ok(())
}

fn checked_add(left: f64, right: f64, channel: &'static str) -> Result<f64, CinematicAovError> {
    let sum = left + right;
    if sum.is_finite() {
        Ok(canonical_zero(sum))
    } else {
        Err(CinematicAovError::NonFiniteChannel { channel })
    }
}

fn average(sum: f64, count: u32) -> f64 {
    if count == 0 {
        0.0
    } else {
        sum / f64::from(count)
    }
}

fn normalized(sum: [f64; 3]) -> [f64; 3] {
    let scale = sum
        .iter()
        .fold(0.0_f64, |largest, value| largest.max(value.abs()));
    if scale == 0.0 || !scale.is_finite() {
        return [0.0; 3];
    }
    let scaled = [sum[0] / scale, sum[1] / scale, sum[2] / scale];
    let norm =
        scale * (scaled[0] * scaled[0] + scaled[1] * scaled[1] + scaled[2] * scaled[2]).sqrt();
    if norm == 0.0 || !norm.is_finite() {
        [0.0; 3]
    } else {
        [sum[0] / norm, sum[1] / norm, sum[2] / norm]
    }
}

fn sample_variance(pixel: CommonPixel) -> f64 {
    if pixel.accepted_count < 2 {
        0.0
    } else {
        pixel.m2_y / f64::from(pixel.accepted_count - 1)
    }
}

fn validity_mask(common: CommonPixel, final_diagnostic: FinalPixel) -> u32 {
    let mut mask = 0;
    if common.primary_count != 0 {
        mask |= validity::PRIMARY;
    }
    if common.albedo_count != 0 {
        mask |= validity::ALBEDO;
    }
    if common.authored_shading_normal_count != 0 {
        mask |= validity::AUTHORED_SHADING_NORMAL;
    }
    if common.previous_motion_count != 0 {
        mask |= validity::PREVIOUS_MOTION;
    }
    if final_diagnostic.nearest_primary.object_palette_index != 0 {
        mask |= validity::OBJECT_ID;
    }
    if final_diagnostic.nearest_primary.material_palette_index != 0 {
        mask |= validity::MATERIAL_ID;
    }
    if common.accepted_count != 0 {
        mask |= validity::CONTRIBUTION_SPLIT;
    }
    mask
}

fn build_exr_channels(
    profile: CinematicAovProfile,
    pixel_count: usize,
    beauty_mean_xyz: impl Fn(usize) -> [f64; 3],
    sample_count: impl Fn(usize) -> u32,
    common: Option<&[CommonPixel]>,
    final_diagnostic: Option<&[FinalPixel]>,
) -> Result<Vec<Channel>, CinematicAovError> {
    let mut channels = Vec::new();
    channels
        .try_reserve_exact(profile.float_channel_count() as usize)
        .map_err(|_| CinematicAovError::AllocationRefused)?;
    for component in 0..3 {
        let name = ["R", "G", "B"][component];
        channels.push(float_channel(name, pixel_count, |pixel| {
            let rgb = xyz_to_linear_srgb(xyz_e_to_d65(beauty_mean_xyz(pixel)));
            rgb[component]
        })?);
    }

    if let Some(common) = common {
        for component in 0..3 {
            let name = ["albedo.R", "albedo.G", "albedo.B"][component];
            channels.push(float_channel(name, pixel_count, |pixel| {
                average(
                    common[pixel].albedo_sum[component],
                    common[pixel].albedo_count,
                )
            })?);
        }
        for component in 0..3 {
            let name = ["normal.X", "normal.Y", "normal.Z"][component];
            channels.push(float_channel(name, pixel_count, |pixel| {
                normalized(common[pixel].shading_normal_sum)[component]
            })?);
        }
        channels.push(float_channel("depth.Z", pixel_count, |pixel| {
            average(common[pixel].depth_sum_m, common[pixel].primary_count)
        })?);
        channels.push(float_channel("primary.coverage", pixel_count, |pixel| {
            average(
                f64::from(common[pixel].primary_count),
                common[pixel].accepted_count,
            )
        })?);
        channels.push(float_channel("variance.Y", pixel_count, |pixel| {
            sample_variance(common[pixel])
        })?);
        for component in 0..2 {
            let name = ["motion.prev.X", "motion.prev.Y"][component];
            channels.push(float_channel(name, pixel_count, |pixel| {
                average(
                    common[pixel].previous_motion_sum_pixels[component],
                    common[pixel].previous_motion_count,
                )
            })?);
        }
    }

    if let (Some(common), Some(final_diagnostic)) = (common, final_diagnostic) {
        for component in 0..3 {
            let name = ["normal_geom.X", "normal_geom.Y", "normal_geom.Z"][component];
            channels.push(float_channel(name, pixel_count, |pixel| {
                normalized(final_diagnostic[pixel].geometric_normal_sum)[component]
            })?);
        }
        for (prefix, selector) in [("direct", 0_u8), ("indirect", 1_u8), ("emission", 2_u8)] {
            for component in 0..3 {
                let suffix = ["R", "G", "B"][component];
                let name = format!("{prefix}.{suffix}");
                channels.push(float_channel(&name, pixel_count, |pixel| {
                    let xyz = match selector {
                        0 => final_diagnostic[pixel].direct_xyz_sum,
                        1 => final_diagnostic[pixel].indirect_xyz_sum,
                        _ => final_diagnostic[pixel].emission_xyz_sum,
                    };
                    let divisor = f64::from(sample_count(pixel));
                    let rgb = xyz_to_linear_srgb(xyz_e_to_d65(xyz.map(|value| value / divisor)));
                    rgb[component]
                })?);
            }
        }
        channels.push(float_channel("id.object", pixel_count, |pixel| {
            f64::from(final_diagnostic[pixel].nearest_primary.object_palette_index)
        })?);
        channels.push(float_channel("id.material", pixel_count, |pixel| {
            f64::from(
                final_diagnostic[pixel]
                    .nearest_primary
                    .material_palette_index,
            )
        })?);
        channels.push(float_channel("samples", pixel_count, |pixel| {
            f64::from(common[pixel].accepted_count)
        })?);
        channels.push(float_channel(
            "diagnostic.validity",
            pixel_count,
            |pixel| f64::from(validity_mask(common[pixel], final_diagnostic[pixel])),
        )?);
    }
    if channels.len() != profile.float_channel_count() as usize {
        return Err(CinematicAovError::InternalChannelCount);
    }
    Ok(channels)
}

fn float_channel(
    name: &str,
    len: usize,
    mut sample: impl FnMut(usize) -> f64,
) -> Result<Channel, CinematicAovError> {
    let mut data = Vec::new();
    data.try_reserve_exact(len)
        .map_err(|_| CinematicAovError::AllocationRefused)?;
    for pixel in 0..len {
        let value = sample(pixel);
        if !value.is_finite() {
            return Err(CinematicAovError::NonFiniteChannel {
                channel: "EXR finalization",
            });
        }
        #[allow(clippy::cast_possible_truncation)]
        let value = value as f32;
        if !value.is_finite() {
            return Err(CinematicAovError::NonFiniteChannel {
                channel: "EXR finalization",
            });
        }
        data.push(if value == 0.0 { 0.0 } else { value });
    }
    Ok(Channel {
        name: name.to_string(),
        ty: PixelType::Float,
        data,
    })
}

fn encode_object_palette(values: &[u64]) -> String {
    let mut encoded = String::from("0=unavailable");
    for (zero_based, value) in values.iter().enumerate() {
        encoded.push(';');
        encoded.push_str(&(zero_based + 1).to_string());
        encoded.push('=');
        encoded.push_str(&value.to_string());
    }
    encoded
}

fn encode_material_palette(values: &[ContentHash]) -> String {
    let mut encoded = String::from("0=unavailable");
    for (zero_based, value) in values.iter().enumerate() {
        encoded.push(';');
        encoded.push_str(&(zero_based + 1).to_string());
        encoded.push('=');
        encoded.push_str(&value.to_hex());
    }
    encoded
}

fn f64_bits_string(value: f64) -> String {
    format!(
        "{}@0x{:016x}",
        canonical_zero(value),
        canonical_zero(value).to_bits()
    )
}

const fn sampler_name(sampler: Sampler) -> &'static str {
    match sampler {
        Sampler::Iid => "iid-philox",
        Sampler::OwenSobol => "owen-sobol",
    }
}

const fn direct_strategy_name(strategy: DirectStrategy) -> &'static str {
    match strategy {
        DirectStrategy::NeeOnly => "nee-only",
        DirectStrategy::BsdfOnly => "bsdf-only",
        DirectStrategy::Mis => "mis",
    }
}

const fn canonical_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_blake3::hash_domain;

    fn hash(label: &str) -> ContentHash {
        hash_domain("org.frankensim.test.aov", label.as_bytes())
    }

    fn config(profile: CinematicAovProfile) -> CinematicAovConfig {
        CinematicAovConfig::new(
            profile,
            CinematicAovProvenance::try_new(
                7,
                1.0,
                0.5,
                1.5,
                hash("trajectory"),
                hash("scene"),
                hash("composition"),
            )
            .unwrap(),
            CinematicAovLimits::default(),
        )
    }

    #[test]
    fn g0_profile_counts_and_4k_retained_admission_are_explicit() {
        assert_eq!(CinematicAovProfile::BeautyOnly.float_channel_count(), 3);
        assert_eq!(CinematicAovProfile::DailyCore.float_channel_count(), 14);
        assert_eq!(
            CinematicAovProfile::FinalDiagnostic.float_channel_count(),
            30
        );
        let bytes = retained_bytes(3_840 * 2_160, CinematicAovProfile::FinalDiagnostic).unwrap();
        assert!(
            bytes.checked_mul(2).unwrap() <= CinematicAovLimits::default().max_retained_bytes()
        );
    }

    #[test]
    fn g0_configuration_identity_binds_profile_provenance_and_limits() {
        let daily = config(CinematicAovProfile::DailyCore);
        let final_profile = config(CinematicAovProfile::FinalDiagnostic);
        assert_eq!(daily.identity(), daily.identity());
        assert_ne!(daily.identity(), final_profile.identity());
    }

    #[test]
    fn g0_nearest_category_uses_distance_then_absolute_sample() {
        let primary = |object_palette_index| AlignedAovPrimary {
            primitive_index: object_palette_index as usize,
            object_palette_index,
            material_palette_index: 1,
            albedo_linear_rgb: Some([0.1, 0.2, 0.3]),
            geometric_normal_world: [0.0, 0.0, 1.0],
            shading_normal_world: [0.0, 0.0, 1.0],
            has_authored_shading_normal: false,
            depth_m: 1.0,
            previous_motion_pixels: Some([1.0, 2.0]),
        };
        let sample = |absolute_sample, jitter, object_palette_index| AlignedAovSample {
            beauty_xyz: [1.0, 1.0, 1.0],
            direct_xyz: [1.0, 0.0, 0.0],
            indirect_xyz: [0.0, 1.0, 0.0],
            emission_xyz: [0.0, 0.0, 1.0],
            pixel_jitter: jitter,
            absolute_sample,
            primary: Some(primary(object_palette_index)),
        };
        let mut pixel = FinalPixel::EMPTY;
        pixel.push(sample(9, [0.9, 0.9], 1)).unwrap();
        pixel.push(sample(12, [0.5, 0.5], 2)).unwrap();
        pixel.push(sample(4, [0.5, 0.5], 3)).unwrap();
        assert_eq!(pixel.nearest_primary.object_palette_index, 3);
        assert_eq!(pixel.nearest_primary.absolute_sample, 4);
    }

    #[test]
    fn g0_background_is_finite_zero_with_only_split_validity() {
        let sample = AlignedAovSample {
            beauty_xyz: [0.0; 3],
            direct_xyz: [0.0; 3],
            indirect_xyz: [0.0; 3],
            emission_xyz: [0.0; 3],
            pixel_jitter: [0.5, 0.5],
            absolute_sample: 0,
            primary: None,
        };
        let mut common = CommonPixel::EMPTY;
        let mut final_pixel = FinalPixel::EMPTY;
        common.push(sample).unwrap();
        final_pixel.push(sample).unwrap();
        assert_eq!(
            validity_mask(common, final_pixel),
            validity::CONTRIBUTION_SPLIT
        );
        assert_eq!(average(common.depth_sum_m, common.primary_count), 0.0);
        assert_eq!(normalized(common.shading_normal_sum), [0.0; 3]);
    }

    #[test]
    fn g0_float_channel_refuses_finite_f64_that_overflows_exr_float() {
        assert_eq!(
            float_channel("R", 1, |_| f64::MAX),
            Err(CinematicAovError::NonFiniteChannel {
                channel: "EXR finalization"
            })
        );
        assert_eq!(
            float_channel("R", 1, |_| -f64::MAX),
            Err(CinematicAovError::NonFiniteChannel {
                channel: "EXR finalization"
            })
        );
        let boundary = float_channel("R", 1, |_| f64::from(f32::MAX)).unwrap();
        assert_eq!(boundary.data, [f32::MAX]);
    }
}
