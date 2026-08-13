//! Deterministic animation-aware RGB denoising.
//!
//! The filter consumes raw RGB plus aligned motion/visibility guides and an
//! optional result from the immediately preceding frame. Spatial weights use
//! scene-linear luminance contrast, matching the measured Monte Carlo
//! luminance variance, so isolated chromatic energy is not mistaken for a
//! stable scene edge. Its public result
//! type has private fields and is always biased: there is deliberately no
//! constructor or conversion that can relabel filtered pixels as a raw
//! estimator.  Motion is `previous - current` in raster pixels, matching the
//! cinematic AOV convention used by `fs-render` without depending on that
//! crate.

use core::fmt;
use std::mem::size_of;

/// Frozen implementation/version identity for temporal denoising.
pub const TEMPORAL_DENOISE_PIPELINE_VERSION: &str = "fs-img-temporal-rgb-nearest-reprojection-v3";

/// Exact byte length of a canonical temporal-denoiser configuration.
pub const TEMPORAL_DENOISE_CONFIG_CANONICAL_BYTES: usize = 48;

/// Largest supported number of spatial à-trous passes.
pub const MAX_TEMPORAL_DENOISE_SPATIAL_ITERATIONS: u8 = 8;

const CONFIG_MAGIC: [u8; 8] = *b"FSTDNV3\0";
const CONFIG_VERSION: u16 = 3;
const B3: [f64; 5] = [1.0 / 16.0, 4.0 / 16.0, 6.0 / 16.0, 4.0 / 16.0, 1.0 / 16.0];
// Inverse of fs-render's frozen E-XYZ -> Bradford D65 -> linear-sRGB transform,
// selecting the source CIE-Y row carried by `variance_luminance`.
const RAW_CIE_Y_FROM_LINEAR_SRGB: [f64; 3] = [
    0.222_853_687_870_429_03,
    0.708_672_666_023_707_7,
    0.068_473_658_307_786_68,
];

/// Frozen mapping from current pixels to the previous raster.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TemporalReprojection {
    /// Add `motion.prev` to the integer current-pixel centre and select the
    /// nearest previous integer-pixel centre. Half-pixel ties round toward the
    /// larger coordinate. Samples whose rounded centre is outside the raster
    /// are rejected.
    NearestPixelCenterV1,
}

impl TemporalReprojection {
    const fn tag(self) -> u8 {
        match self {
            Self::NearestPixelCenterV1 => 1,
        }
    }
}

/// Whether a frame may consume the supplied immediately preceding history.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TemporalFrameBoundary {
    /// The frame is continuous with the preceding frame.
    Continuous,
    /// A camera/shot discontinuity resets every pixel to history length one.
    Cut,
}

/// Versioned temporal and joint-RGB spatial filtering parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TemporalDenoiseConfig {
    /// Frozen previous-frame coordinate reconstruction.
    pub reprojection: TemporalReprojection,
    /// Maximum accepted history length per pixel.
    pub max_history_frames: u16,
    /// Hard cap on the contribution of reprojected history, in `[0, 1)`.
    pub max_history_weight: f32,
    /// Maximum absolute change in primary coverage for temporal/spatial guide
    /// compatibility.
    pub coverage_tolerance: f32,
    /// Absolute axial-depth agreement band in metres.
    pub depth_absolute_tolerance_m: f32,
    /// Relative axial-depth agreement band.
    pub depth_relative_tolerance: f32,
    /// Minimum cosine between compatible surface normals.
    pub normal_cosine_threshold: f32,
    /// Number of current-frame standard deviations added to the 3x3 RGB
    /// history clamp interval.
    pub neighborhood_clamp_stddev: f32,
    /// Small nonnegative variance threshold used when two variance estimates
    /// are numerically indistinguishable from zero.
    pub variance_floor: f32,
    /// Number of joint-RGB 5x5 à-trous refinement passes. At least one pass is
    /// required so reset frames remain explicitly filtered derivatives.
    pub spatial_iterations: u8,
    /// Absolute scene-linear luminance edge-stopping sigma for spatial
    /// refinement. The field name is retained for source compatibility with
    /// version two configurations.
    pub spatial_sigma_rgb: f32,
}

impl Default for TemporalDenoiseConfig {
    fn default() -> Self {
        Self {
            reprojection: TemporalReprojection::NearestPixelCenterV1,
            max_history_frames: 32,
            max_history_weight: 0.95,
            coverage_tolerance: 0.2,
            depth_absolute_tolerance_m: 0.002,
            depth_relative_tolerance: 0.01,
            normal_cosine_threshold: 0.8,
            neighborhood_clamp_stddev: 1.0,
            variance_floor: 1.0e-8,
            spatial_iterations: 2,
            spatial_sigma_rgb: 0.15,
        }
    }
}

impl TemporalDenoiseConfig {
    /// Validate and return the exact canonical configuration identity retained
    /// by every result.
    ///
    /// # Errors
    /// [`TemporalDenoiseError::InvalidConfig`] for an invalid field.
    pub fn identity(self) -> Result<TemporalDenoiseConfigIdentity, TemporalDenoiseError> {
        validate_config(self)?;
        let mut bytes = [0_u8; TEMPORAL_DENOISE_CONFIG_CANONICAL_BYTES];
        bytes[..8].copy_from_slice(&CONFIG_MAGIC);
        bytes[8..10].copy_from_slice(&CONFIG_VERSION.to_le_bytes());
        bytes[10] = self.reprojection.tag();
        bytes[11] = self.spatial_iterations;
        bytes[12..14].copy_from_slice(&self.max_history_frames.to_le_bytes());
        // bytes 14..16 are reserved zero bytes.
        let floats = [
            self.max_history_weight,
            self.coverage_tolerance,
            self.depth_absolute_tolerance_m,
            self.depth_relative_tolerance,
            self.normal_cosine_threshold,
            self.neighborhood_clamp_stddev,
            self.variance_floor,
            self.spatial_sigma_rgb,
        ];
        for (index, value) in floats.into_iter().enumerate() {
            let offset = 16 + index * 4;
            bytes[offset..offset + 4].copy_from_slice(&canonical_f32_bits(value).to_le_bytes());
        }
        Ok(TemporalDenoiseConfigIdentity { bytes })
    }
}

/// Exact canonical identity of one validated temporal-denoiser configuration.
/// The field is private so arbitrary bytes cannot be presented as an admitted
/// configuration identity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TemporalDenoiseConfigIdentity {
    bytes: [u8; TEMPORAL_DENOISE_CONFIG_CANONICAL_BYTES],
}

impl TemporalDenoiseConfigIdentity {
    /// Canonical versioned configuration bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; TEMPORAL_DENOISE_CONFIG_CANONICAL_BYTES] {
        &self.bytes
    }
}

/// Caller-owned ceilings for allocations made by one filtering call.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TemporalDenoiseLimits {
    /// Maximum admitted pixels.
    pub max_pixels: u64,
    /// Maximum bytes newly allocated for the returned history plus one spatial
    /// scratch plane. Borrowed raw inputs and borrowed previous history remain
    /// the caller's separately budgeted memory.
    pub max_new_bytes: u64,
}

impl TemporalDenoiseLimits {
    /// Envelope for one native 3840x2160 frame, including both optional u64 ID
    /// planes and one spatial scratch RGB plane.
    #[must_use]
    pub const fn reference_4k() -> Self {
        Self {
            max_pixels: 3_840 * 2_160,
            max_new_bytes: 640 * 1024 * 1024,
        }
    }
}

impl Default for TemporalDenoiseLimits {
    fn default() -> Self {
        Self::reference_4k()
    }
}

/// Explicit, aligned inputs for one raw frame.
///
/// Every plane is row-major with exactly `width * height` elements. Surface
/// guides are valid when coverage is positive. Background pixels use zero
/// depth, zero normal, zero motion, and zero optional IDs. In an optional ID
/// plane, zero on a covered sample means that identity is unavailable; only
/// pairs of nonzero values act as exact categorical equality evidence.
#[derive(Clone, Copy, Debug)]
pub struct TemporalDenoiseInput<'a> {
    /// Zero-based or segment-local monotonically increasing frame index.
    pub frame_index: u64,
    /// Uniform raw Monte Carlo samples contributing to every pixel estimate,
    /// or the hard sample ceiling when `sample_counts_per_pixel` is present.
    /// `variance_luminance` is a sample variance, so the denoiser divides it
    /// by the exact count before using it as variance of the pixel mean.
    pub samples_per_pixel: u32,
    /// Optional exact row-major per-pixel raw sample counts from an adaptive
    /// renderer. Every count must be in `1..=samples_per_pixel`. When absent,
    /// the uniform `samples_per_pixel` count applies to every pixel.
    pub sample_counts_per_pixel: Option<&'a [u32]>,
    /// Raster width.
    pub width: usize,
    /// Raster height.
    pub height: usize,
    /// Scene-linear red estimator samples.
    pub red: &'a [f32],
    /// Scene-linear green estimator samples.
    pub green: &'a [f32],
    /// Scene-linear blue estimator samples.
    pub blue: &'a [f32],
    /// Previous-minus-current raster displacement X, in pixels.
    pub motion_prev_x: &'a [f32],
    /// Previous-minus-current raster displacement Y, in pixels.
    pub motion_prev_y: &'a [f32],
    /// Optional per-pixel proof that the previous-motion vector exists.
    /// `false` rejects temporal history even when the stored zero sentinel
    /// would otherwise look like a valid static correspondence.
    pub previous_motion_valid: Option<&'a [bool]>,
    /// Positive surface axial depth in metres; zero for background.
    pub axial_depth_m: &'a [f32],
    /// World-space unit shading-normal X; zero for background.
    pub normal_x: &'a [f32],
    /// World-space unit shading-normal Y; zero for background.
    pub normal_y: &'a [f32],
    /// World-space unit shading-normal Z; zero for background.
    pub normal_z: &'a [f32],
    /// Primary-hit fraction in `[0, 1]`.
    pub primary_coverage: &'a [f32],
    /// Nonnegative raw luminance variance proxy.
    pub variance_luminance: &'a [f32],
    /// Optional stable object IDs. Zero means background or unavailable;
    /// nonzero values are exact equality labels.
    pub object_ids: Option<&'a [u64]>,
    /// Optional stable material IDs. Zero means background or unavailable;
    /// nonzero values are exact equality labels.
    pub material_ids: Option<&'a [u64]>,
}

/// Mandatory honesty tag returned by [`TemporalDenoisedFrame::provenance`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TemporalDenoiseProvenance {
    /// Animation-aware temporal and spatial filtering was applied. The pixels
    /// are biased and must never serve as a raw/converged estimator.
    BiasedTemporalDenoisedV3 {
        /// Exact versioned configuration identity.
        config_identity: TemporalDenoiseConfigIdentity,
    },
}

/// A private-field animation history and biased display derivative.
///
/// Only [`temporal_denoise_rgb`] constructs this type. It deliberately exposes
/// no raw-provenance constructor, setter, or conversion.
#[derive(Debug, PartialEq)]
pub struct TemporalDenoisedFrame {
    frame_index: u64,
    width: usize,
    height: usize,
    linear_rgb: [Vec<f32>; 3],
    axial_depth_m: Vec<f32>,
    normal: Vec<[f32; 3]>,
    primary_coverage: Vec<f32>,
    estimate_variance: Vec<f32>,
    history_length: Vec<u16>,
    object_ids: Option<Vec<u64>>,
    material_ids: Option<Vec<u64>>,
    config_identity: TemporalDenoiseConfigIdentity,
    retained_bytes: u64,
}

impl TemporalDenoisedFrame {
    /// Source frame index.
    #[must_use]
    pub const fn frame_index(&self) -> u64 {
        self.frame_index
    }

    /// Raster width.
    #[must_use]
    pub const fn width(&self) -> usize {
        self.width
    }

    /// Raster height.
    #[must_use]
    pub const fn height(&self) -> usize {
        self.height
    }

    /// Planar scene-linear biased RGB samples in red, green, blue order.
    #[must_use]
    pub fn linear_rgb(&self) -> [&[f32]; 3] {
        [
            &self.linear_rgb[0],
            &self.linear_rgb[1],
            &self.linear_rgb[2],
        ]
    }

    /// Accepted temporal history length per pixel. Rejected/reset pixels have
    /// length one.
    #[must_use]
    pub fn history_length(&self) -> &[u16] {
        &self.history_length
    }

    /// Exact retained payload bytes charged at admission.
    #[must_use]
    pub const fn retained_bytes(&self) -> u64 {
        self.retained_bytes
    }

    /// Exact versioned configuration identity.
    #[must_use]
    pub const fn config_identity(&self) -> TemporalDenoiseConfigIdentity {
        self.config_identity
    }

    /// Permanent biased provenance. There is no raw relabeling surface.
    #[must_use]
    pub const fn provenance(&self) -> TemporalDenoiseProvenance {
        TemporalDenoiseProvenance::BiasedTemporalDenoisedV3 {
            config_identity: self.config_identity,
        }
    }
}

/// Fail-closed validation, sequence, or resource error.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TemporalDenoiseError {
    /// Width and height must both be nonzero and multiplication must fit.
    InvalidDimensions {
        /// Supplied width.
        width: usize,
        /// Supplied height.
        height: usize,
    },
    /// One input plane has the wrong number of elements.
    Shape {
        /// Stable plane name.
        field: &'static str,
        /// Required elements.
        expected: usize,
        /// Supplied elements.
        got: usize,
    },
    /// One configuration field is outside its admitted range.
    InvalidConfig {
        /// Stable field name.
        field: &'static str,
    },
    /// One sample is nonfinite or violates the guide grammar.
    InvalidSample {
        /// Stable plane/guide name.
        field: &'static str,
        /// Row-major pixel index.
        index: usize,
        /// Stable refusal reason.
        reason: &'static str,
    },
    /// A nonzero first/segment-local frame lacked an explicit cut reset.
    MissingInitialReset {
        /// Supplied frame index.
        frame_index: u64,
    },
    /// Supplied history is not from the immediately preceding frame.
    FrameOrder {
        /// Required next frame index.
        expected: u64,
        /// Supplied current frame index.
        got: u64,
    },
    /// Previous history ended at `u64::MAX`, so no contiguous successor exists.
    FrameIndexOverflow,
    /// Continuous history has another raster shape.
    HistoryShapeMismatch,
    /// Continuous history was produced by another configuration.
    HistoryConfigMismatch,
    /// Optional guide presence changed without a cut.
    HistoryGuideLayoutMismatch {
        /// `object_ids` or `material_ids`.
        field: &'static str,
    },
    /// Pixel admission exceeded the caller ceiling.
    PixelLimit {
        /// Required pixels.
        requested: u64,
        /// Caller ceiling.
        limit: u64,
    },
    /// Newly allocated result-plus-scratch bytes exceeded the caller ceiling.
    WorkingMemoryLimit {
        /// Exact required bytes.
        requested: u64,
        /// Caller ceiling.
        limit: u64,
    },
    /// Checked size arithmetic overflowed.
    SizeOverflow {
        /// Stable quantity name.
        context: &'static str,
    },
    /// The allocator refused a previously admitted vector.
    AllocationRefused {
        /// Stable allocation name.
        resource: &'static str,
        /// Exact vector payload bytes requested.
        requested: u64,
    },
}

impl fmt::Display for TemporalDenoiseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDimensions { width, height } => {
                write!(f, "invalid temporal-denoiser raster {width}x{height}")
            }
            Self::Shape {
                field,
                expected,
                got,
            } => write!(f, "{field}: expected {expected} elements, got {got}"),
            Self::InvalidConfig { field } => {
                write!(f, "invalid temporal-denoiser configuration field {field}")
            }
            Self::InvalidSample {
                field,
                index,
                reason,
            } => write!(f, "invalid {field} at pixel {index}: {reason}"),
            Self::MissingInitialReset { frame_index } => write!(
                f,
                "frame {frame_index} has no history and is not an explicit cut reset"
            ),
            Self::FrameOrder { expected, got } => {
                write!(f, "expected frame {expected} after history, got {got}")
            }
            Self::FrameIndexOverflow => write!(f, "history frame index has no successor"),
            Self::HistoryShapeMismatch => write!(f, "continuous history raster shape changed"),
            Self::HistoryConfigMismatch => {
                write!(f, "continuous history denoiser configuration changed")
            }
            Self::HistoryGuideLayoutMismatch { field } => {
                write!(f, "continuous history changed optional guide {field}")
            }
            Self::PixelLimit { requested, limit } => {
                write!(
                    f,
                    "temporal denoiser needs {requested} pixels above limit {limit}"
                )
            }
            Self::WorkingMemoryLimit { requested, limit } => write!(
                f,
                "temporal denoiser needs {requested} new bytes above limit {limit}"
            ),
            Self::SizeOverflow { context } => {
                write!(f, "temporal denoiser size overflowed for {context}")
            }
            Self::AllocationRefused {
                resource,
                requested,
            } => write!(
                f,
                "allocator refused {requested} bytes for temporal denoiser {resource}"
            ),
        }
    }
}

impl std::error::Error for TemporalDenoiseError {}

/// Filter one frame using an optional immediately preceding biased history.
///
/// A missing history is admitted only for frame zero or an explicit cut. A cut
/// validates but ignores history. Continuous history must match frame order,
/// shape, configuration, and optional-ID layout. Reprojection and every guide
/// test fail closed per pixel; rejected pixels restart at history length one.
/// The final joint-RGB spatial pass uses one edge weight for all three channels
/// and therefore cannot independently skew channel edges.
///
/// # Errors
/// Returns [`TemporalDenoiseError`] before publication for invalid dimensions,
/// configuration, plane shape/content, history, or resource admission.
#[allow(clippy::too_many_lines)]
pub fn temporal_denoise_rgb(
    input: TemporalDenoiseInput<'_>,
    previous: Option<&TemporalDenoisedFrame>,
    boundary: TemporalFrameBoundary,
    config: TemporalDenoiseConfig,
    limits: TemporalDenoiseLimits,
) -> Result<TemporalDenoisedFrame, TemporalDenoiseError> {
    let config_identity = config.identity()?;
    let pixel_count = validate_input(input)?;
    let pixel_count_u64 =
        u64::try_from(pixel_count).map_err(|_| TemporalDenoiseError::SizeOverflow {
            context: "pixel count",
        })?;
    if pixel_count_u64 > limits.max_pixels {
        return Err(TemporalDenoiseError::PixelLimit {
            requested: pixel_count_u64,
            limit: limits.max_pixels,
        });
    }
    validate_history(input, previous, boundary, config_identity)?;

    let retained_bytes = retained_bytes(
        pixel_count,
        input.object_ids.is_some(),
        input.material_ids.is_some(),
    )?;
    let scratch_bytes = vector_bytes::<[f32; 3]>(pixel_count, "spatial RGB scratch")?;
    let required_new_bytes =
        retained_bytes
            .checked_add(scratch_bytes)
            .ok_or(TemporalDenoiseError::SizeOverflow {
                context: "result plus scratch bytes",
            })?;
    if required_new_bytes > limits.max_new_bytes {
        return Err(TemporalDenoiseError::WorkingMemoryLimit {
            requested: required_new_bytes,
            limit: limits.max_new_bytes,
        });
    }

    let mut rgb = [
        try_filled(pixel_count, 0.0_f32, "red result")?,
        try_filled(pixel_count, 0.0_f32, "green result")?,
        try_filled(pixel_count, 0.0_f32, "blue result")?,
    ];
    let mut history_length = try_filled(pixel_count, 1_u16, "history length")?;
    let mut estimate_variance = try_filled(pixel_count, 0.0_f32, "estimate variance")?;

    let use_history = boundary == TemporalFrameBoundary::Continuous && previous.is_some();
    let previous = previous.filter(|_| use_history);
    for index in 0..pixel_count {
        let current_rgb = [input.red[index], input.green[index], input.blue[index]];
        let current_variance = input.variance_luminance[index] / sample_count(input, index) as f32;
        let Some(history) = previous else {
            for channel in 0..3 {
                rgb[channel][index] = current_rgb[channel];
            }
            estimate_variance[index] = current_variance;
            continue;
        };
        let Some(previous_index) = reproject_nearest(input, index) else {
            for channel in 0..3 {
                rgb[channel][index] = current_rgb[channel];
            }
            estimate_variance[index] = current_variance;
            continue;
        };
        if !guides_match_history(input, index, history, previous_index, config) {
            for channel in 0..3 {
                rgb[channel][index] = current_rgb[channel];
            }
            estimate_variance[index] = current_variance;
            continue;
        }

        let old_length = history.history_length[previous_index];
        let count_weight = f64::from(old_length) / (f64::from(old_length) + 1.0);
        let old_variance = f64::from(history.estimate_variance[previous_index]);
        let current_variance_f64 = f64::from(current_variance);
        let variance_sum = old_variance + current_variance_f64;
        let variance_weight = if variance_sum > f64::from(config.variance_floor) {
            current_variance_f64 / variance_sum
        } else {
            count_weight
        };
        let history_weight = count_weight
            .min(variance_weight)
            .min(f64::from(config.max_history_weight));
        let history_rgb =
            core::array::from_fn(|channel| history.linear_rgb[channel][previous_index]);
        let clamped_history = clamp_history_rgb(input, index, history_rgb, config);
        let current_weight = 1.0 - history_weight;
        for channel in 0..3 {
            rgb[channel][index] = (history_weight * f64::from(clamped_history[channel])
                + current_weight * f64::from(current_rgb[channel]))
                as f32;
        }
        let combined_variance = history_weight * history_weight * old_variance
            + current_weight * current_weight * current_variance_f64;
        estimate_variance[index] = combined_variance as f32;
        history_length[index] = old_length.saturating_add(1).min(config.max_history_frames);
    }

    let mut scratch = [
        try_filled(pixel_count, 0.0_f32, "spatial red scratch")?,
        try_filled(pixel_count, 0.0_f32, "spatial green scratch")?,
        try_filled(pixel_count, 0.0_f32, "spatial blue scratch")?,
    ];
    for iteration in 0..config.spatial_iterations {
        spatial_atrous_pass(
            input,
            &rgb,
            &estimate_variance,
            &mut scratch,
            config,
            iteration,
        );
        core::mem::swap(&mut rgb, &mut scratch);
    }

    let axial_depth_m = try_copy(input.axial_depth_m, "depth history")?;
    let primary_coverage = try_copy(input.primary_coverage, "coverage history")?;
    let mut normal = try_filled(pixel_count, [0.0; 3], "normal history")?;
    for (index, sample) in normal.iter_mut().enumerate() {
        *sample = [
            input.normal_x[index],
            input.normal_y[index],
            input.normal_z[index],
        ];
    }
    let object_ids = input
        .object_ids
        .map(|values| try_copy(values, "object-ID history"))
        .transpose()?;
    let material_ids = input
        .material_ids
        .map(|values| try_copy(values, "material-ID history"))
        .transpose()?;

    Ok(TemporalDenoisedFrame {
        frame_index: input.frame_index,
        width: input.width,
        height: input.height,
        linear_rgb: rgb,
        axial_depth_m,
        normal,
        primary_coverage,
        estimate_variance,
        history_length,
        object_ids,
        material_ids,
        config_identity,
        retained_bytes,
    })
}

fn validate_config(config: TemporalDenoiseConfig) -> Result<(), TemporalDenoiseError> {
    if config.max_history_frames == 0 {
        return invalid_config("max_history_frames");
    }
    if !config.max_history_weight.is_finite()
        || config.max_history_weight < 0.0
        || config.max_history_weight >= 1.0
    {
        return invalid_config("max_history_weight");
    }
    if !config.coverage_tolerance.is_finite() || !(0.0..=1.0).contains(&config.coverage_tolerance) {
        return invalid_config("coverage_tolerance");
    }
    for (field, value) in [
        (
            "depth_absolute_tolerance_m",
            config.depth_absolute_tolerance_m,
        ),
        ("depth_relative_tolerance", config.depth_relative_tolerance),
        (
            "neighborhood_clamp_stddev",
            config.neighborhood_clamp_stddev,
        ),
        ("variance_floor", config.variance_floor),
    ] {
        if !value.is_finite() || value < 0.0 {
            return invalid_config(field);
        }
    }
    if !config.normal_cosine_threshold.is_finite()
        || !(-1.0..=1.0).contains(&config.normal_cosine_threshold)
    {
        return invalid_config("normal_cosine_threshold");
    }
    if config.spatial_iterations == 0
        || config.spatial_iterations > MAX_TEMPORAL_DENOISE_SPATIAL_ITERATIONS
    {
        return invalid_config("spatial_iterations");
    }
    if !config.spatial_sigma_rgb.is_finite() || config.spatial_sigma_rgb <= 0.0 {
        return invalid_config("spatial_sigma_rgb");
    }
    Ok(())
}

fn invalid_config<T>(field: &'static str) -> Result<T, TemporalDenoiseError> {
    Err(TemporalDenoiseError::InvalidConfig { field })
}

#[allow(clippy::too_many_lines)] // one fixed plane/guide grammar, validated in refusal order
#[allow(clippy::float_cmp)] // zero is the exact documented invalid/background sentinel
fn validate_input(input: TemporalDenoiseInput<'_>) -> Result<usize, TemporalDenoiseError> {
    if input.samples_per_pixel == 0 {
        return Err(TemporalDenoiseError::InvalidConfig {
            field: "samples_per_pixel",
        });
    }
    if input.width == 0 || input.height == 0 {
        return Err(TemporalDenoiseError::InvalidDimensions {
            width: input.width,
            height: input.height,
        });
    }
    let pixel_count =
        input
            .width
            .checked_mul(input.height)
            .ok_or(TemporalDenoiseError::InvalidDimensions {
                width: input.width,
                height: input.height,
            })?;
    for (field, len) in [
        ("red", input.red.len()),
        ("green", input.green.len()),
        ("blue", input.blue.len()),
        ("motion_prev_x", input.motion_prev_x.len()),
        ("motion_prev_y", input.motion_prev_y.len()),
        ("axial_depth_m", input.axial_depth_m.len()),
        ("normal_x", input.normal_x.len()),
        ("normal_y", input.normal_y.len()),
        ("normal_z", input.normal_z.len()),
        ("primary_coverage", input.primary_coverage.len()),
        ("variance_luminance", input.variance_luminance.len()),
    ] {
        if len != pixel_count {
            return Err(TemporalDenoiseError::Shape {
                field,
                expected: pixel_count,
                got: len,
            });
        }
    }
    if let Some(sample_counts) = input.sample_counts_per_pixel {
        if sample_counts.len() != pixel_count {
            return Err(TemporalDenoiseError::Shape {
                field: "sample_counts_per_pixel",
                expected: pixel_count,
                got: sample_counts.len(),
            });
        }
        if let Some((index, &samples)) = sample_counts
            .iter()
            .enumerate()
            .find(|(_, samples)| **samples == 0 || **samples > input.samples_per_pixel)
        {
            return Err(TemporalDenoiseError::InvalidSample {
                field: "sample_counts_per_pixel",
                index,
                reason: if samples == 0 {
                    "zero"
                } else {
                    "above-sample-ceiling"
                },
            });
        }
    }
    if let Some(previous_motion_valid) = input.previous_motion_valid
        && previous_motion_valid.len() != pixel_count
    {
        return Err(TemporalDenoiseError::Shape {
            field: "previous_motion_valid",
            expected: pixel_count,
            got: previous_motion_valid.len(),
        });
    }
    for (field, values) in [
        ("object_ids", input.object_ids),
        ("material_ids", input.material_ids),
    ] {
        if let Some(values) = values
            && values.len() != pixel_count
        {
            return Err(TemporalDenoiseError::Shape {
                field,
                expected: pixel_count,
                got: values.len(),
            });
        }
    }
    for index in 0..pixel_count {
        for (field, value) in [
            ("red", input.red[index]),
            ("green", input.green[index]),
            ("blue", input.blue[index]),
            ("motion_prev_x", input.motion_prev_x[index]),
            ("motion_prev_y", input.motion_prev_y[index]),
            ("axial_depth_m", input.axial_depth_m[index]),
            ("normal_x", input.normal_x[index]),
            ("normal_y", input.normal_y[index]),
            ("normal_z", input.normal_z[index]),
            ("primary_coverage", input.primary_coverage[index]),
            ("variance_luminance", input.variance_luminance[index]),
        ] {
            if !value.is_finite() {
                return Err(TemporalDenoiseError::InvalidSample {
                    field,
                    index,
                    reason: "nonfinite",
                });
            }
        }
        let coverage = input.primary_coverage[index];
        if !(0.0..=1.0).contains(&coverage) {
            return Err(TemporalDenoiseError::InvalidSample {
                field: "primary_coverage",
                index,
                reason: "outside [0, 1]",
            });
        }
        if input.variance_luminance[index] < 0.0 {
            return Err(TemporalDenoiseError::InvalidSample {
                field: "variance_luminance",
                index,
                reason: "negative",
            });
        }
        let depth = input.axial_depth_m[index];
        let normal = [
            input.normal_x[index],
            input.normal_y[index],
            input.normal_z[index],
        ];
        let motion = [input.motion_prev_x[index], input.motion_prev_y[index]];
        if coverage == 0.0 {
            if depth != 0.0 {
                return invalid_sample("axial_depth_m", index, "background must be zero");
            }
            if normal != [0.0; 3] {
                return invalid_sample("normal", index, "background must be zero");
            }
            if motion != [0.0; 2] {
                return invalid_sample("motion_prev", index, "background must be zero");
            }
            if input.object_ids.is_some_and(|ids| ids[index] != 0) {
                return invalid_sample("object_ids", index, "background must be zero");
            }
            if input.material_ids.is_some_and(|ids| ids[index] != 0) {
                return invalid_sample("material_ids", index, "background must be zero");
            }
        } else {
            if depth <= 0.0 {
                return invalid_sample("axial_depth_m", index, "covered sample must be positive");
            }
            let normal_length_squared = normal
                .into_iter()
                .map(f64::from)
                .map(|value| value * value)
                .sum::<f64>();
            if !(0.9801..=1.0201).contains(&normal_length_squared) {
                return invalid_sample("normal", index, "covered sample must be unit length");
            }
        }
    }
    Ok(pixel_count)
}

fn invalid_sample<T>(
    field: &'static str,
    index: usize,
    reason: &'static str,
) -> Result<T, TemporalDenoiseError> {
    Err(TemporalDenoiseError::InvalidSample {
        field,
        index,
        reason,
    })
}

fn validate_history(
    input: TemporalDenoiseInput<'_>,
    previous: Option<&TemporalDenoisedFrame>,
    boundary: TemporalFrameBoundary,
    config_identity: TemporalDenoiseConfigIdentity,
) -> Result<(), TemporalDenoiseError> {
    let Some(previous) = previous else {
        if input.frame_index != 0 && boundary != TemporalFrameBoundary::Cut {
            return Err(TemporalDenoiseError::MissingInitialReset {
                frame_index: input.frame_index,
            });
        }
        return Ok(());
    };
    let expected = previous
        .frame_index
        .checked_add(1)
        .ok_or(TemporalDenoiseError::FrameIndexOverflow)?;
    if input.frame_index != expected {
        return Err(TemporalDenoiseError::FrameOrder {
            expected,
            got: input.frame_index,
        });
    }
    if boundary == TemporalFrameBoundary::Cut {
        return Ok(());
    }
    if input.width != previous.width || input.height != previous.height {
        return Err(TemporalDenoiseError::HistoryShapeMismatch);
    }
    if config_identity != previous.config_identity {
        return Err(TemporalDenoiseError::HistoryConfigMismatch);
    }
    for (field, present, history_present) in [
        (
            "object_ids",
            input.object_ids.is_some(),
            previous.object_ids.is_some(),
        ),
        (
            "material_ids",
            input.material_ids.is_some(),
            previous.material_ids.is_some(),
        ),
    ] {
        if present != history_present {
            return Err(TemporalDenoiseError::HistoryGuideLayoutMismatch { field });
        }
    }
    Ok(())
}

fn reproject_nearest(input: TemporalDenoiseInput<'_>, index: usize) -> Option<usize> {
    if input
        .previous_motion_valid
        .is_some_and(|valid| !valid[index])
    {
        return None;
    }
    let x = index % input.width;
    let y = index / input.width;
    let previous_x = (x as f64 + f64::from(input.motion_prev_x[index]) + 0.5).floor();
    let previous_y = (y as f64 + f64::from(input.motion_prev_y[index]) + 0.5).floor();
    if previous_x < 0.0
        || previous_y < 0.0
        || previous_x >= input.width as f64
        || previous_y >= input.height as f64
    {
        return None;
    }
    Some(previous_y as usize * input.width + previous_x as usize)
}

fn guides_match_history(
    input: TemporalDenoiseInput<'_>,
    current: usize,
    previous: &TemporalDenoisedFrame,
    old: usize,
    config: TemporalDenoiseConfig,
) -> bool {
    if !coverage_matches(
        input.primary_coverage[current],
        previous.primary_coverage[old],
        config.coverage_tolerance,
    ) {
        return false;
    }
    let current_surface = input.primary_coverage[current] > 0.0;
    let old_surface = previous.primary_coverage[old] > 0.0;
    if current_surface != old_surface {
        return false;
    }
    if current_surface
        && (!depth_matches(
            input.axial_depth_m[current],
            previous.axial_depth_m[old],
            config,
        ) || !normal_matches(
            [
                input.normal_x[current],
                input.normal_y[current],
                input.normal_z[current],
            ],
            previous.normal[old],
            config.normal_cosine_threshold,
        ))
    {
        return false;
    }
    if let (Some(current_ids), Some(previous_ids)) = (input.object_ids, &previous.object_ids)
        && categorical_ids_conflict(current_ids[current], previous_ids[old])
    {
        return false;
    }
    if let (Some(current_ids), Some(previous_ids)) = (input.material_ids, &previous.material_ids)
        && categorical_ids_conflict(current_ids[current], previous_ids[old])
    {
        return false;
    }
    true
}

fn guides_match_current(
    input: TemporalDenoiseInput<'_>,
    center: usize,
    sample: usize,
    config: TemporalDenoiseConfig,
) -> bool {
    if !coverage_matches(
        input.primary_coverage[center],
        input.primary_coverage[sample],
        config.coverage_tolerance,
    ) {
        return false;
    }
    let center_surface = input.primary_coverage[center] > 0.0;
    let sample_surface = input.primary_coverage[sample] > 0.0;
    if center_surface != sample_surface {
        return false;
    }
    if center_surface
        && (!depth_matches(
            input.axial_depth_m[center],
            input.axial_depth_m[sample],
            config,
        ) || !normal_matches(
            [
                input.normal_x[center],
                input.normal_y[center],
                input.normal_z[center],
            ],
            [
                input.normal_x[sample],
                input.normal_y[sample],
                input.normal_z[sample],
            ],
            config.normal_cosine_threshold,
        ))
    {
        return false;
    }
    if input
        .object_ids
        .is_some_and(|ids| categorical_ids_conflict(ids[center], ids[sample]))
        || input
            .material_ids
            .is_some_and(|ids| categorical_ids_conflict(ids[center], ids[sample]))
    {
        return false;
    }
    true
}

fn categorical_ids_conflict(left: u64, right: u64) -> bool {
    left != 0 && right != 0 && left != right
}

fn coverage_matches(left: f32, right: f32, tolerance: f32) -> bool {
    f64::from(left - right).abs() <= f64::from(tolerance)
}

fn depth_matches(left: f32, right: f32, config: TemporalDenoiseConfig) -> bool {
    let left = f64::from(left);
    let right = f64::from(right);
    let tolerance = f64::from(config.depth_absolute_tolerance_m)
        + f64::from(config.depth_relative_tolerance) * left.abs().max(right.abs());
    (left - right).abs() <= tolerance
}

fn normal_matches(left: [f32; 3], right: [f32; 3], threshold: f32) -> bool {
    let dot = left
        .into_iter()
        .zip(right)
        .map(|(a, b)| f64::from(a) * f64::from(b))
        .sum::<f64>();
    dot >= f64::from(threshold)
}

fn clamp_history_rgb(
    input: TemporalDenoiseInput<'_>,
    center: usize,
    history_rgb: [f32; 3],
    config: TemporalDenoiseConfig,
) -> [f32; 3] {
    let x = center % input.width;
    let y = center / input.width;
    let center_rgb = [input.red[center], input.green[center], input.blue[center]];
    let mut minimum = center_rgb;
    let mut maximum = center_rgb;
    for dy in -1_isize..=1 {
        let sy = y.saturating_add_signed(dy).min(input.height - 1);
        for dx in -1_isize..=1 {
            let sx = x.saturating_add_signed(dx).min(input.width - 1);
            let sample = sy * input.width + sx;
            if !guides_match_current(input, center, sample, config) {
                continue;
            }
            let rgb = [input.red[sample], input.green[sample], input.blue[sample]];
            for channel in 0..3 {
                minimum[channel] = minimum[channel].min(rgb[channel]);
                maximum[channel] = maximum[channel].max(rgb[channel]);
            }
        }
    }
    let mean_variance =
        f64::from(input.variance_luminance[center]) / f64::from(sample_count(input, center));
    let expansion = f64::from(config.neighborhood_clamp_stddev) * mean_variance.sqrt();
    // Clamp the history vector along the line from the current sample instead
    // of clipping channels independently. One shared factor preserves gray
    // neutrality and constant-hue lines while still entering every component's
    // admitted neighborhood interval.
    let mut factor = 1.0_f64;
    for channel in 0..3 {
        let center = f64::from(center_rgb[channel]);
        let difference = f64::from(history_rgb[channel]) - center;
        let bound = if difference > 0.0 {
            (f64::from(maximum[channel]) + expansion - center) / difference
        } else if difference < 0.0 {
            (f64::from(minimum[channel]) - expansion - center) / difference
        } else {
            1.0
        };
        factor = factor.min(bound.clamp(0.0, 1.0));
    }
    core::array::from_fn(|channel| {
        (f64::from(center_rgb[channel])
            + factor * (f64::from(history_rgb[channel]) - f64::from(center_rgb[channel])))
            as f32
    })
}

fn sample_count(input: TemporalDenoiseInput<'_>, index: usize) -> u32 {
    input
        .sample_counts_per_pixel
        .map_or(input.samples_per_pixel, |counts| counts[index])
}

fn spatial_atrous_pass(
    input: TemporalDenoiseInput<'_>,
    current: &[Vec<f32>; 3],
    estimate_variance: &[f32],
    next: &mut [Vec<f32>; 3],
    config: TemporalDenoiseConfig,
    iteration: u8,
) {
    let step = 1_isize << iteration;
    let sigma_squared = f64::from(config.spatial_sigma_rgb).powi(2);
    for center in 0..(input.width * input.height) {
        let x = center % input.width;
        let y = center / input.width;
        let center_rgb = core::array::from_fn(|channel| current[channel][center]);
        let mut accumulated = [0.0_f64; 3];
        let mut weight_sum = 0.0_f64;
        for (kernel_y, weight_y) in B3.into_iter().enumerate() {
            let dy = (kernel_y.cast_signed() - 2).saturating_mul(step);
            let sy = y.saturating_add_signed(dy).min(input.height - 1);
            for (kernel_x, weight_x) in B3.into_iter().enumerate() {
                let dx = (kernel_x.cast_signed() - 2).saturating_mul(step);
                let sx = x.saturating_add_signed(dx).min(input.width - 1);
                let sample = sy * input.width + sx;
                if !guides_match_current(input, center, sample, config) {
                    continue;
                }
                let differences = core::array::from_fn::<_, 3, _>(|channel| {
                    f64::from(current[channel][sample] - center_rgb[channel])
                });
                let luminance_difference = differences
                    .into_iter()
                    .zip(RAW_CIE_Y_FROM_LINEAR_SRGB)
                    .map(|(difference, weight)| difference * weight)
                    .sum::<f64>();
                let luminance_difference_squared = luminance_difference * luminance_difference;
                let rgb_difference_squared = differences
                    .into_iter()
                    .map(|difference| difference * difference)
                    .sum::<f64>();
                // `estimate_variance` is the sample-count-adjusted variance of
                // the luminance mean, so it is compared only with luminance
                // contrast. Comparing it with an RGB norm incorrectly protects
                // high-chroma, low-luminance Monte Carlo outliers as color
                // edges. Geometry and identity guides above remain hard stops.
                let noise_scale_squared =
                    f64::from(estimate_variance[center]) + f64::from(estimate_variance[sample]);
                // Retain full RGB edge stopping for converged color detail, but
                // progressively distrust chroma contrast when the only aligned
                // variance evidence says this pixel is noisy. This is an
                // explicit display bias, not a raw-energy clamp.
                let chroma_confidence = sigma_squared / (sigma_squared + noise_scale_squared);
                let difference_squared = luminance_difference_squared
                    + chroma_confidence
                        * (rgb_difference_squared - luminance_difference_squared).max(0.0);
                let color_scale_squared = sigma_squared + noise_scale_squared;
                let weight =
                    weight_x * weight_y * (-difference_squared / color_scale_squared).exp();
                for (channel, sum) in accumulated.iter_mut().enumerate() {
                    *sum += weight * f64::from(current[channel][sample]);
                }
                weight_sum += weight;
            }
        }
        let output = if weight_sum > 0.0 {
            accumulated.map(|sum| (sum / weight_sum) as f32)
        } else {
            center_rgb
        };
        for channel in 0..3 {
            next[channel][center] = output[channel];
        }
    }
}

fn retained_bytes(
    pixel_count: usize,
    has_object_ids: bool,
    has_material_ids: bool,
) -> Result<u64, TemporalDenoiseError> {
    let mut bytes = 0_u64;
    for (context, component_bytes) in [
        ("RGB result", size_of::<[f32; 3]>()),
        ("depth history", size_of::<f32>()),
        ("normal history", size_of::<[f32; 3]>()),
        ("coverage history", size_of::<f32>()),
        ("variance history", size_of::<f32>()),
        ("history length", size_of::<u16>()),
    ] {
        bytes = bytes
            .checked_add(vector_bytes_raw(pixel_count, component_bytes, context)?)
            .ok_or(TemporalDenoiseError::SizeOverflow {
                context: "retained bytes",
            })?;
    }
    for (present, context) in [
        (has_object_ids, "object-ID history"),
        (has_material_ids, "material-ID history"),
    ] {
        if present {
            bytes = bytes
                .checked_add(vector_bytes::<u64>(pixel_count, context)?)
                .ok_or(TemporalDenoiseError::SizeOverflow {
                    context: "retained bytes",
                })?;
        }
    }
    Ok(bytes)
}

fn vector_bytes<T>(count: usize, context: &'static str) -> Result<u64, TemporalDenoiseError> {
    vector_bytes_raw(count, size_of::<T>(), context)
}

fn vector_bytes_raw(
    count: usize,
    element_bytes: usize,
    context: &'static str,
) -> Result<u64, TemporalDenoiseError> {
    let bytes = count
        .checked_mul(element_bytes)
        .ok_or(TemporalDenoiseError::SizeOverflow { context })?;
    u64::try_from(bytes).map_err(|_| TemporalDenoiseError::SizeOverflow { context })
}

fn try_filled<T: Copy>(
    count: usize,
    value: T,
    resource: &'static str,
) -> Result<Vec<T>, TemporalDenoiseError> {
    let requested = vector_bytes::<T>(count, resource)?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| TemporalDenoiseError::AllocationRefused {
            resource,
            requested,
        })?;
    values.resize(count, value);
    Ok(values)
}

fn try_copy<T: Copy>(source: &[T], resource: &'static str) -> Result<Vec<T>, TemporalDenoiseError> {
    let requested = vector_bytes::<T>(source.len(), resource)?;
    let mut values = Vec::new();
    values.try_reserve_exact(source.len()).map_err(|_| {
        TemporalDenoiseError::AllocationRefused {
            resource,
            requested,
        }
    })?;
    values.extend_from_slice(source);
    Ok(values)
}

fn canonical_f32_bits(value: f32) -> u32 {
    if value == 0.0 { 0 } else { value.to_bits() }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_categorical_id_is_unavailable_not_an_equality_or_conflict_claim() {
        for (left, right) in [(0, 0), (0, 1), (1, 0), (1, 1)] {
            assert!(
                !categorical_ids_conflict(left, right),
                "{left} versus {right}"
            );
        }
        assert!(categorical_ids_conflict(1, 2));
    }

    #[derive(Clone)]
    struct Fixture {
        frame_index: u64,
        samples_per_pixel: u32,
        sample_counts_per_pixel: Option<Vec<u32>>,
        width: usize,
        height: usize,
        red: Vec<f32>,
        green: Vec<f32>,
        blue: Vec<f32>,
        motion_x: Vec<f32>,
        motion_y: Vec<f32>,
        motion_valid: Option<Vec<bool>>,
        depth: Vec<f32>,
        normal_x: Vec<f32>,
        normal_y: Vec<f32>,
        normal_z: Vec<f32>,
        coverage: Vec<f32>,
        variance: Vec<f32>,
        object_ids: Option<Vec<u64>>,
        material_ids: Option<Vec<u64>>,
    }

    impl Fixture {
        fn surface(width: usize, height: usize, frame_index: u64, rgb: [f32; 3]) -> Self {
            let count = width * height;
            Self {
                frame_index,
                samples_per_pixel: 1,
                sample_counts_per_pixel: None,
                width,
                height,
                red: vec![rgb[0]; count],
                green: vec![rgb[1]; count],
                blue: vec![rgb[2]; count],
                motion_x: vec![0.0; count],
                motion_y: vec![0.0; count],
                motion_valid: Some(vec![true; count]),
                depth: vec![1.0; count],
                normal_x: vec![0.0; count],
                normal_y: vec![0.0; count],
                normal_z: vec![1.0; count],
                coverage: vec![1.0; count],
                variance: vec![0.01; count],
                object_ids: Some(vec![1; count]),
                material_ids: Some(vec![2; count]),
            }
        }

        fn input(&self) -> TemporalDenoiseInput<'_> {
            TemporalDenoiseInput {
                frame_index: self.frame_index,
                samples_per_pixel: self.samples_per_pixel,
                sample_counts_per_pixel: self.sample_counts_per_pixel.as_deref(),
                width: self.width,
                height: self.height,
                red: &self.red,
                green: &self.green,
                blue: &self.blue,
                motion_prev_x: &self.motion_x,
                motion_prev_y: &self.motion_y,
                previous_motion_valid: self.motion_valid.as_deref(),
                axial_depth_m: &self.depth,
                normal_x: &self.normal_x,
                normal_y: &self.normal_y,
                normal_z: &self.normal_z,
                primary_coverage: &self.coverage,
                variance_luminance: &self.variance,
                object_ids: self.object_ids.as_deref(),
                material_ids: self.material_ids.as_deref(),
            }
        }
    }

    fn run(
        fixture: &Fixture,
        previous: Option<&TemporalDenoisedFrame>,
        boundary: TemporalFrameBoundary,
        config: TemporalDenoiseConfig,
    ) -> Result<TemporalDenoisedFrame, TemporalDenoiseError> {
        temporal_denoise_rgb(
            fixture.input(),
            previous,
            boundary,
            config,
            TemporalDenoiseLimits::default(),
        )
    }

    fn rgb_mse(actual: [&[f32]; 3], expected: [f32; 3]) -> f64 {
        let mut squared_error = 0.0_f64;
        for channel in 0..3 {
            squared_error += actual[channel]
                .iter()
                .map(|&value| {
                    let difference = f64::from(value - expected[channel]);
                    difference * difference
                })
                .sum::<f64>();
        }
        squared_error / (actual[0].len() * 3) as f64
    }

    #[test]
    fn static_history_and_joint_spatial_filter_reduce_variance() {
        let mut first = Fixture::surface(8, 8, 0, [0.5, 0.4, 0.3]);
        let mut second = Fixture::surface(8, 8, 1, [0.5, 0.4, 0.3]);
        for index in 0..64 {
            let sign = if (index + index / 8) % 2 == 0 {
                1.0
            } else {
                -1.0
            };
            first.red[index] += 0.12 * sign;
            first.green[index] += 0.10 * sign;
            first.blue[index] += 0.08 * sign;
            second.red[index] -= 0.12 * sign;
            second.green[index] -= 0.10 * sign;
            second.blue[index] -= 0.08 * sign;
        }
        let config = TemporalDenoiseConfig {
            spatial_sigma_rgb: 1.0,
            neighborhood_clamp_stddev: 2.0,
            ..TemporalDenoiseConfig::default()
        };
        let first_output = run(&first, None, TemporalFrameBoundary::Continuous, config).unwrap();
        let second_output = run(
            &second,
            Some(&first_output),
            TemporalFrameBoundary::Continuous,
            config,
        )
        .unwrap();
        assert!(
            rgb_mse(second_output.linear_rgb(), [0.5, 0.4, 0.3])
                < rgb_mse([&second.red, &second.green, &second.blue], [0.5, 0.4, 0.3]) * 0.1
        );
        assert!(
            second_output
                .history_length()
                .iter()
                .all(|&length| length == 2)
        );
    }

    #[test]
    fn measured_variance_prevents_a_firefly_from_becoming_a_color_edge() {
        let mut noisy = Fixture::surface(9, 9, 0, [0.2, 0.2, 0.2]);
        let center = 4 * 9 + 4;
        noisy.red[center] = 4.0;
        noisy.green[center] = 4.0;
        noisy.blue[center] = 4.0;
        noisy.variance[center] = 16.0;
        let config = TemporalDenoiseConfig {
            spatial_iterations: 2,
            ..TemporalDenoiseConfig::default()
        };

        let filtered = run(&noisy, None, TemporalFrameBoundary::Cut, config).unwrap();
        let filtered_center = filtered.linear_rgb()[0][center];
        assert!(
            filtered_center < 1.0,
            "high measured variance must not preserve a {filtered_center} firefly as a color edge"
        );
        assert!(
            filtered.linear_rgb()[0]
                .iter()
                .all(|&value| value.is_finite() && value >= 0.0)
        );

        noisy.variance.fill(0.0);
        let stable_edge = run(&noisy, None, TemporalFrameBoundary::Cut, config).unwrap();
        assert!(
            stable_edge.linear_rgb()[0][center] > 3.9,
            "without variance evidence the same contrast remains a protected color feature"
        );
    }

    #[test]
    fn luminance_variance_suppresses_a_high_chroma_blue_firefly() {
        let mut noisy = Fixture::surface(9, 9, 0, [0.2, 0.2, 0.2]);
        let center = 4 * 9 + 4;
        noisy.blue[center] = 4.0;
        noisy.variance[center] = 0.5;
        let config = TemporalDenoiseConfig {
            spatial_iterations: 2,
            ..TemporalDenoiseConfig::default()
        };

        let filtered = run(&noisy, None, TemporalFrameBoundary::Cut, config).unwrap();
        let filtered_blue = filtered.linear_rgb()[2][center];
        assert!(
            filtered_blue < 1.0,
            "measured luminance variance must not preserve a {filtered_blue} blue firefly"
        );

        noisy.variance.fill(0.0);
        let stable_edge = run(&noisy, None, TemporalFrameBoundary::Cut, config).unwrap();
        assert!(
            stable_edge.linear_rgb()[2][center] > 2.5,
            "without variance evidence the same chromatic contrast remains protected"
        );
    }

    #[test]
    fn moving_identity_edge_and_disocclusion_restart_without_trails() {
        let mut first = Fixture::surface(5, 1, 0, [1.0, 0.0, 0.0]);
        let mut second = Fixture::surface(5, 1, 1, [1.0, 0.0, 0.0]);
        for index in 3..5 {
            for fixture in [&mut first, &mut second] {
                fixture.red[index] = 0.0;
                fixture.blue[index] = 1.0;
                fixture.depth[index] = 2.0;
                fixture.object_ids.as_mut().unwrap()[index] = 9;
                fixture.material_ids.as_mut().unwrap()[index] = 10;
            }
        }
        // The blue object moves left by one pixel. Its newly covered pixel 2
        // points at previous blue pixel 3; a disoccluded red pixel 3 is forced
        // to point at the incompatible previous red pixel 2.
        second.red[2] = 0.0;
        second.blue[2] = 1.0;
        second.depth[2] = 2.0;
        second.object_ids.as_mut().unwrap()[2] = 9;
        second.material_ids.as_mut().unwrap()[2] = 10;
        second.motion_x[2] = 1.0;
        second.motion_x[3] = -1.0;
        let config = TemporalDenoiseConfig {
            spatial_sigma_rgb: 1.0,
            ..TemporalDenoiseConfig::default()
        };
        let first_output = run(&first, None, TemporalFrameBoundary::Continuous, config).unwrap();
        let output = run(
            &second,
            Some(&first_output),
            TemporalFrameBoundary::Continuous,
            config,
        )
        .unwrap();
        assert!(output.linear_rgb()[2][2] > 0.99 && output.linear_rgb()[0][2] < 0.01);
        assert_eq!(output.history_length()[2], 2);
        assert!(output.linear_rgb()[2][3] > 0.99 && output.linear_rgb()[0][3] < 0.01);
        assert_eq!(output.history_length()[3], 1);
    }

    #[test]
    fn unavailable_zero_motion_restarts_instead_of_reusing_static_history() {
        let first = Fixture::surface(2, 1, 0, [1.0, 0.0, 0.0]);
        let mut second = Fixture::surface(2, 1, 1, [0.0, 1.0, 0.0]);
        second.motion_valid.as_mut().unwrap()[0] = false;
        let config = TemporalDenoiseConfig {
            neighborhood_clamp_stddev: 0.0,
            ..TemporalDenoiseConfig::default()
        };
        let first_output = run(&first, None, TemporalFrameBoundary::Continuous, config).unwrap();
        let output = run(
            &second,
            Some(&first_output),
            TemporalFrameBoundary::Continuous,
            config,
        )
        .unwrap();
        assert_eq!(output.history_length(), &[1, 2]);
        assert_eq!(output.linear_rgb()[0][0].to_bits(), 0.0_f32.to_bits());
        assert_eq!(output.linear_rgb()[1][0].to_bits(), 1.0_f32.to_bits());
    }

    #[test]
    fn depth_normal_coverage_and_bounds_rejections_reset_history() {
        let first = Fixture::surface(5, 1, 0, [0.25; 3]);
        let mut second = Fixture::surface(5, 1, 1, [0.75; 3]);
        second.depth[0] = 2.0;
        second.normal_x[1] = 1.0;
        second.normal_z[1] = 0.0;
        second.coverage[2] = 0.5;
        second.motion_x[3] = 100.0;
        // Pixel 4 becomes background: an explicit surface disocclusion.
        second.red[4] = 0.05;
        second.green[4] = 0.05;
        second.blue[4] = 0.05;
        second.depth[4] = 0.0;
        second.normal_z[4] = 0.0;
        second.coverage[4] = 0.0;
        second.object_ids.as_mut().unwrap()[4] = 0;
        second.material_ids.as_mut().unwrap()[4] = 0;
        let config = TemporalDenoiseConfig {
            coverage_tolerance: 0.1,
            neighborhood_clamp_stddev: 0.0,
            ..TemporalDenoiseConfig::default()
        };
        let first_output = run(&first, None, TemporalFrameBoundary::Continuous, config).unwrap();
        let output = run(
            &second,
            Some(&first_output),
            TemporalFrameBoundary::Continuous,
            config,
        )
        .unwrap();
        assert_eq!(output.history_length(), &[1, 1, 1, 1, 1]);
        assert!(output.linear_rgb().iter().all(|plane| plane[4] < 0.06));
    }

    #[test]
    fn cut_reset_equals_fresh_reset_and_changes_may_cross_cut() {
        let first = Fixture::surface(3, 2, 0, [1.0, 0.0, 0.0]);
        let mut cut = Fixture::surface(3, 2, 1, [0.0, 1.0, 0.0]);
        cut.object_ids = None;
        cut.material_ids = None;
        let first_output = run(
            &first,
            None,
            TemporalFrameBoundary::Continuous,
            TemporalDenoiseConfig::default(),
        )
        .unwrap();
        let changed_config = TemporalDenoiseConfig {
            spatial_sigma_rgb: 0.5,
            ..TemporalDenoiseConfig::default()
        };
        let with_old_history = run(
            &cut,
            Some(&first_output),
            TemporalFrameBoundary::Cut,
            changed_config,
        )
        .unwrap();
        let fresh = run(&cut, None, TemporalFrameBoundary::Cut, changed_config).unwrap();
        assert_eq!(with_old_history, fresh);
        assert!(fresh.history_length().iter().all(|&length| length == 1));
    }

    #[test]
    fn malformed_guides_shapes_and_frame_order_fail_closed() {
        let mut malformed = Fixture::surface(2, 2, 0, [0.2; 3]);
        malformed.motion_x.pop();
        assert_eq!(
            run(
                &malformed,
                None,
                TemporalFrameBoundary::Continuous,
                TemporalDenoiseConfig::default()
            ),
            Err(TemporalDenoiseError::Shape {
                field: "motion_prev_x",
                expected: 4,
                got: 3,
            })
        );

        let valid = Fixture::surface(2, 2, 0, [0.2; 3]);
        let mut zero_samples = valid.clone();
        zero_samples.samples_per_pixel = 0;
        assert_eq!(
            run(
                &zero_samples,
                None,
                TemporalFrameBoundary::Continuous,
                TemporalDenoiseConfig::default()
            ),
            Err(TemporalDenoiseError::InvalidConfig {
                field: "samples_per_pixel"
            })
        );
        let mut malformed_counts = valid.clone();
        malformed_counts.samples_per_pixel = 8;
        malformed_counts.sample_counts_per_pixel = Some(vec![2; 3]);
        assert_eq!(
            run(
                &malformed_counts,
                None,
                TemporalFrameBoundary::Continuous,
                TemporalDenoiseConfig::default()
            ),
            Err(TemporalDenoiseError::Shape {
                field: "sample_counts_per_pixel",
                expected: 4,
                got: 3,
            })
        );
        for (samples, reason) in [(0, "zero"), (9, "above-sample-ceiling")] {
            let mut invalid_counts = valid.clone();
            invalid_counts.samples_per_pixel = 8;
            invalid_counts.sample_counts_per_pixel = Some(vec![2, samples, 2, 2]);
            assert_eq!(
                run(
                    &invalid_counts,
                    None,
                    TemporalFrameBoundary::Continuous,
                    TemporalDenoiseConfig::default()
                ),
                Err(TemporalDenoiseError::InvalidSample {
                    field: "sample_counts_per_pixel",
                    index: 1,
                    reason,
                })
            );
        }
        assert_eq!(
            run(
                &valid,
                None,
                TemporalFrameBoundary::Continuous,
                TemporalDenoiseConfig {
                    spatial_iterations: 0,
                    ..TemporalDenoiseConfig::default()
                }
            ),
            Err(TemporalDenoiseError::InvalidConfig {
                field: "spatial_iterations"
            })
        );

        let mut invalid = Fixture::surface(2, 2, 0, [0.2; 3]);
        invalid.normal_z[1] = f32::NAN;
        assert!(matches!(
            run(
                &invalid,
                None,
                TemporalFrameBoundary::Continuous,
                TemporalDenoiseConfig::default()
            ),
            Err(TemporalDenoiseError::InvalidSample {
                field: "normal_z",
                index: 1,
                reason: "nonfinite"
            })
        ));

        let missing_reset = Fixture::surface(2, 2, 5, [0.2; 3]);
        assert_eq!(
            run(
                &missing_reset,
                None,
                TemporalFrameBoundary::Continuous,
                TemporalDenoiseConfig::default()
            ),
            Err(TemporalDenoiseError::MissingInitialReset { frame_index: 5 })
        );

        let first = Fixture::surface(2, 2, 0, [0.2; 3]);
        let first_output = run(
            &first,
            None,
            TemporalFrameBoundary::Continuous,
            TemporalDenoiseConfig::default(),
        )
        .unwrap();
        let gap = Fixture::surface(2, 2, 2, [0.2; 3]);
        assert_eq!(
            run(
                &gap,
                Some(&first_output),
                TemporalFrameBoundary::Continuous,
                TemporalDenoiseConfig::default()
            ),
            Err(TemporalDenoiseError::FrameOrder {
                expected: 1,
                got: 2
            })
        );

        let mut changed_layout = Fixture::surface(2, 2, 1, [0.2; 3]);
        assert_eq!(
            run(
                &changed_layout,
                Some(&first_output),
                TemporalFrameBoundary::Continuous,
                TemporalDenoiseConfig {
                    spatial_sigma_rgb: 0.5,
                    ..TemporalDenoiseConfig::default()
                }
            ),
            Err(TemporalDenoiseError::HistoryConfigMismatch)
        );
        changed_layout.object_ids = None;
        assert_eq!(
            run(
                &changed_layout,
                Some(&first_output),
                TemporalFrameBoundary::Continuous,
                TemporalDenoiseConfig::default()
            ),
            Err(TemporalDenoiseError::HistoryGuideLayoutMismatch {
                field: "object_ids"
            })
        );
    }

    #[test]
    fn adaptive_sample_counts_scale_each_pixels_estimator_variance_exactly() {
        let mut fixture = Fixture::surface(2, 1, 0, [0.2; 3]);
        fixture.samples_per_pixel = 8;
        fixture.sample_counts_per_pixel = Some(vec![2, 8]);
        fixture.variance = vec![8.0, 8.0];

        let result = run(
            &fixture,
            None,
            TemporalFrameBoundary::Cut,
            TemporalDenoiseConfig::default(),
        )
        .expect("valid heterogeneous adaptive counts");

        assert_eq!(result.estimate_variance, [4.0, 1.0]);
    }

    #[test]
    fn deterministic_replay_is_bit_exact() {
        fn sequence() -> TemporalDenoisedFrame {
            let mut first = Fixture::surface(7, 3, 0, [0.3, 0.4, 0.5]);
            let mut second = Fixture::surface(7, 3, 1, [0.31, 0.39, 0.52]);
            for index in 0..21 {
                first.red[index] += index as f32 * 0.001;
                second.blue[index] -= index as f32 * 0.0005;
            }
            let config = TemporalDenoiseConfig::default();
            let history = run(&first, None, TemporalFrameBoundary::Continuous, config).unwrap();
            run(
                &second,
                Some(&history),
                TemporalFrameBoundary::Continuous,
                config,
            )
            .unwrap()
        }
        assert_eq!(sequence(), sequence());
    }

    #[test]
    fn joint_rgb_filter_preserves_gray_and_constant_hue_lines() {
        let mut gray = Fixture::surface(6, 2, 0, [0.5; 3]);
        for index in 0..12 {
            let delta = (index as f32 - 5.5) * 0.01;
            gray.red[index] += delta;
            gray.green[index] += delta;
            gray.blue[index] += delta;
        }
        let gray_output = run(
            &gray,
            None,
            TemporalFrameBoundary::Continuous,
            TemporalDenoiseConfig::default(),
        )
        .unwrap();
        let gray_rgb = gray_output.linear_rgb();
        for index in 0..gray_rgb[0].len() {
            assert_eq!(gray_rgb[0][index].to_bits(), gray_rgb[1][index].to_bits());
            assert_eq!(gray_rgb[1][index].to_bits(), gray_rgb[2][index].to_bits());
        }

        let mut hue = Fixture::surface(6, 2, 0, [0.2, 0.4, 0.8]);
        for index in 0..12 {
            let scale = 0.5 + index as f32 * 0.05;
            hue.red[index] = 0.2 * scale;
            hue.green[index] = 0.4 * scale;
            hue.blue[index] = 0.8 * scale;
        }
        let hue_config = TemporalDenoiseConfig {
            spatial_sigma_rgb: 1.0,
            ..TemporalDenoiseConfig::default()
        };
        let hue_output = run(&hue, None, TemporalFrameBoundary::Continuous, hue_config).unwrap();
        let hue_rgb = hue_output.linear_rgb();
        for index in 0..hue_rgb[0].len() {
            assert!((hue_rgb[1][index] - 2.0 * hue_rgb[0][index]).abs() < 2.0e-7);
            assert!((hue_rgb[2][index] - 4.0 * hue_rgb[0][index]).abs() < 4.0e-7);
        }

        let mut hue_next = Fixture::surface(6, 2, 1, [0.2, 0.4, 0.8]);
        for index in 0..12 {
            let scale = 1.1 - index as f32 * 0.04;
            hue_next.red[index] = 0.2 * scale;
            hue_next.green[index] = 0.4 * scale;
            hue_next.blue[index] = 0.8 * scale;
        }
        let temporal_hue = run(
            &hue_next,
            Some(&hue_output),
            TemporalFrameBoundary::Continuous,
            hue_config,
        )
        .unwrap();
        let temporal_rgb = temporal_hue.linear_rgb();
        for index in 0..temporal_rgb[0].len() {
            assert!((temporal_rgb[1][index] - 2.0 * temporal_rgb[0][index]).abs() < 3.0e-7);
            assert!((temporal_rgb[2][index] - 4.0 * temporal_rgb[0][index]).abs() < 6.0e-7);
        }
    }

    #[test]
    fn result_is_permanently_biased_and_identity_is_canonical() {
        let fixture = Fixture::surface(2, 2, 0, [0.4; 3]);
        let config = TemporalDenoiseConfig::default();
        let output = run(&fixture, None, TemporalFrameBoundary::Continuous, config).unwrap();
        assert_eq!(
            output.provenance(),
            TemporalDenoiseProvenance::BiasedTemporalDenoisedV3 {
                config_identity: config.identity().unwrap()
            }
        );
        assert_eq!(&output.config_identity().as_bytes()[..8], &CONFIG_MAGIC);

        let positive_zero = TemporalDenoiseConfig {
            depth_absolute_tolerance_m: 0.0,
            ..config
        };
        let mut signed_zero = positive_zero;
        signed_zero.depth_absolute_tolerance_m = -0.0;
        assert_eq!(
            signed_zero.identity().unwrap(),
            positive_zero.identity().unwrap()
        );
    }

    #[test]
    fn exact_memory_and_pixel_limits_refuse_before_work() {
        let fixture = Fixture::surface(2, 2, 0, [0.4; 3]);
        let config = TemporalDenoiseConfig::default();
        let exact_retained = retained_bytes(4, true, true).unwrap();
        let exact_scratch = vector_bytes::<[f32; 3]>(4, "test").unwrap();
        let exact = exact_retained + exact_scratch;
        let admitted = temporal_denoise_rgb(
            fixture.input(),
            None,
            TemporalFrameBoundary::Continuous,
            config,
            TemporalDenoiseLimits {
                max_pixels: 4,
                max_new_bytes: exact,
            },
        )
        .unwrap();
        assert_eq!(admitted.retained_bytes(), exact_retained);
        assert_eq!(
            temporal_denoise_rgb(
                fixture.input(),
                None,
                TemporalFrameBoundary::Continuous,
                config,
                TemporalDenoiseLimits {
                    max_pixels: 3,
                    max_new_bytes: exact,
                }
            ),
            Err(TemporalDenoiseError::PixelLimit {
                requested: 4,
                limit: 3
            })
        );
        assert_eq!(
            temporal_denoise_rgb(
                fixture.input(),
                None,
                TemporalFrameBoundary::Continuous,
                config,
                TemporalDenoiseLimits {
                    max_pixels: 4,
                    max_new_bytes: exact - 1,
                }
            ),
            Err(TemporalDenoiseError::WorkingMemoryLimit {
                requested: exact,
                limit: exact - 1
            })
        );
    }
}
