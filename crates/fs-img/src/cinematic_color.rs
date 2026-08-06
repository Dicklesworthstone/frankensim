//! Deterministic cinematic color and preview derivatives.
//!
//! The renderer's linear-sRGB planes are immutable inputs. This module emits
//! explicitly display-referred derivatives and retains the exact versioned
//! parameters needed to replay them. It is not an ICC/OCIO implementation and
//! its optional box bloom is a visualization effect, not an optical lens model.

use core::{fmt, mem::size_of};

use crate::film::{hable_filmic_unclamped, srgb_encode};

/// Version of the frozen cinematic color configuration and algorithms.
pub const CINEMATIC_COLOR_PIPELINE_VERSION: u16 = 1;
/// Exact byte length returned by [`CinematicColorConfig::canonical_bytes`].
pub const CINEMATIC_COLOR_CONFIG_CANONICAL_BYTES: usize = 64;
/// Hard ceiling for one preview: 8K UHD.
pub const MAX_CINEMATIC_PREVIEW_PIXELS: usize = 7_680 * 4_320;
/// Hard ceiling for the radius of the bounded box-bloom derivative.
pub const MAX_CINEMATIC_GLARE_RADIUS_PX: u16 = 256;

const MAX_ABS_EXPOSURE_EV: i32 = 32;
const MAX_WHITE_BALANCE_GAIN: f64 = 16.0;
const MAX_LINEAR_COMPONENT: f64 = 1.0e12;
const MAX_GLARE_STRENGTH: f64 = 8.0;

/// Scene-linear working-space interpretation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum CinematicWorkingSpace {
    /// Linear sRGB primaries with a D65 white point.
    LinearSrgbD65 = 1,
}

/// Display target of the encoded preview samples.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum CinematicDisplayTarget {
    /// IEC sRGB transfer function and D65 primaries.
    SrgbD65 = 1,
}

/// Frozen display tone-curve formulation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum CinematicToneCurve {
    /// Historical Hable/Uncharted-2 curve normalized at linear value 11.2.
    HableV1 = 1,
    /// Krzysztof Narkowicz's five-coefficient ACES-filmic fit:
    /// `x(2.51x + 0.03) / (x(2.43x + 0.59) + 0.14)`.
    ///
    /// This is deliberately named as a fit: it is not an ACES reference
    /// rendering transform and does not claim ACES/OCIO conformance.
    AcesFittedNarkowiczV1 = 2,
}

/// Explicit out-of-display-gamut handling after the tone curve.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum CinematicGamutMap {
    /// Clamp each display-linear channel independently to `[0, 1]`.
    ClipV1 = 1,
    /// If a channel exceeds one, scale the whole RGB triple by its maximum.
    /// This preserves non-negative RGB channel ratios but is not perceptual.
    RgbRatioScaleV1 = 2,
}

/// Explicit handling of negative scene-linear channels.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum CinematicNegativePolicy {
    /// Count every negative adjusted channel, then clamp it to zero before
    /// bloom and tone mapping.
    ClampToZeroCountedV1 = 1,
}

/// Integer sample precision of the preview derivative.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum PreviewBitDepth {
    /// Eight bits per encoded channel.
    Eight = 8,
    /// Sixteen bits per encoded channel.
    Sixteen = 16,
}

impl PreviewBitDepth {
    const fn bytes_per_channel(self) -> usize {
        match self {
            Self::Eight => 1,
            Self::Sixteen => 2,
        }
    }
}

/// Deterministic quantization-dither policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PreviewDither {
    /// Round to nearest without dither.
    Disabled,
    /// Add one keyed uniform `[-0.5, 0.5)`-LSB variate before rounding.
    UniformHalfLsbV1 {
        /// Explicit replay seed; it is keyed again by pixel and channel.
        seed: u64,
    },
}

/// Optional scene-linear visualization glare.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PreviewGlare {
    /// Do not add a glare derivative.
    Disabled,
    /// Add a separable, zero-boundary box bloom to values above a threshold.
    BoxBloomV1 {
        /// Radius of each one-dimensional box pass in pixels.
        radius_px: u16,
        /// Threshold in exposed, white-balanced linear-sRGB channel units.
        threshold_linear: f64,
        /// Multiplier applied to the normalized bright-pass convolution.
        strength: f64,
    },
}

/// Complete version-one cinematic color configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CinematicColorConfig {
    /// Interpretation of the three input planes.
    pub working_space: CinematicWorkingSpace,
    /// Interpretation of encoded output samples.
    pub display_target: CinematicDisplayTarget,
    /// Frozen tone-curve implementation.
    pub tone_curve: CinematicToneCurve,
    /// Explicit post-tone-map gamut policy.
    pub gamut_map: CinematicGamutMap,
    /// Explicit pre-tone-map negative-channel policy.
    pub negative_policy: CinematicNegativePolicy,
    /// Artistic exposure in exact powers of two.
    pub exposure_ev: i32,
    /// Positive per-channel gains in linear RGB.
    pub white_balance_gains: [f64; 3],
    /// Integer precision of the output samples.
    pub bit_depth: PreviewBitDepth,
    /// Quantization-dither algorithm and seed.
    pub dither: PreviewDither,
    /// Optional labeled visualization glare.
    pub glare: PreviewGlare,
}

impl CinematicColorConfig {
    /// Frozen reference display configuration for 16-bit sRGB previews.
    #[must_use]
    pub const fn reference_srgb_16() -> Self {
        Self {
            working_space: CinematicWorkingSpace::LinearSrgbD65,
            display_target: CinematicDisplayTarget::SrgbD65,
            tone_curve: CinematicToneCurve::AcesFittedNarkowiczV1,
            gamut_map: CinematicGamutMap::RgbRatioScaleV1,
            negative_policy: CinematicNegativePolicy::ClampToZeroCountedV1,
            exposure_ev: 0,
            white_balance_gains: [1.0; 3],
            bit_depth: PreviewBitDepth::Sixteen,
            dither: PreviewDither::UniformHalfLsbV1 { seed: 0 },
            glare: PreviewGlare::Disabled,
        }
    }

    /// Validate numeric fields and hard algorithm bounds before allocation.
    pub fn validate(self) -> Result<Self, CinematicColorError> {
        if !(-MAX_ABS_EXPOSURE_EV..=MAX_ABS_EXPOSURE_EV).contains(&self.exposure_ev) {
            return Err(CinematicColorError::InvalidConfig {
                field: "exposure_ev",
                reason: "must lie in [-32, 32]",
            });
        }
        for gain in self.white_balance_gains {
            if !gain.is_finite() || !(0.0..=MAX_WHITE_BALANCE_GAIN).contains(&gain) || gain == 0.0 {
                return Err(CinematicColorError::InvalidConfig {
                    field: "white_balance_gains",
                    reason: "each gain must be finite and in (0, 16]",
                });
            }
        }
        if let PreviewGlare::BoxBloomV1 {
            radius_px,
            threshold_linear,
            strength,
        } = self.glare
        {
            if radius_px == 0 || radius_px > MAX_CINEMATIC_GLARE_RADIUS_PX {
                return Err(CinematicColorError::InvalidConfig {
                    field: "glare.radius_px",
                    reason: "must lie in [1, 256]",
                });
            }
            if !threshold_linear.is_finite()
                || !(0.0..=MAX_LINEAR_COMPONENT).contains(&threshold_linear)
            {
                return Err(CinematicColorError::InvalidConfig {
                    field: "glare.threshold_linear",
                    reason: "must be finite and in [0, 1e12]",
                });
            }
            if !strength.is_finite()
                || !(0.0..=MAX_GLARE_STRENGTH).contains(&strength)
                || strength == 0.0
            {
                return Err(CinematicColorError::InvalidConfig {
                    field: "glare.strength",
                    reason: "must be finite and in (0, 8]; use Disabled for zero",
                });
            }
        }
        Ok(self)
    }

    /// Canonical fixed-width encoding of every semantic parameter.
    ///
    /// Disabled variant payloads are encoded as zero. All floating-point
    /// zeros are normalized to positive zero, so equivalent configurations
    /// have byte-identical encodings.
    pub fn canonical_bytes(self) -> Result<Vec<u8>, CinematicColorError> {
        let config = self.validate()?;
        let mut bytes = Vec::with_capacity(CINEMATIC_COLOR_CONFIG_CANONICAL_BYTES);
        bytes.extend_from_slice(&CINEMATIC_COLOR_PIPELINE_VERSION.to_le_bytes());
        bytes.push(config.working_space as u8);
        bytes.push(config.display_target as u8);
        bytes.push(config.tone_curve as u8);
        bytes.push(config.gamut_map as u8);
        bytes.push(config.negative_policy as u8);
        bytes.extend_from_slice(&config.exposure_ev.to_le_bytes());
        for gain in config.white_balance_gains {
            push_canonical_f64(&mut bytes, gain);
        }
        bytes.push(config.bit_depth as u8);
        match config.dither {
            PreviewDither::Disabled => {
                bytes.push(0);
                bytes.extend_from_slice(&0_u64.to_le_bytes());
            }
            PreviewDither::UniformHalfLsbV1 { seed } => {
                bytes.push(1);
                bytes.extend_from_slice(&seed.to_le_bytes());
            }
        }
        match config.glare {
            PreviewGlare::Disabled => {
                bytes.push(0);
                bytes.extend_from_slice(&0_u16.to_le_bytes());
                bytes.extend_from_slice(&0_u64.to_le_bytes());
                bytes.extend_from_slice(&0_u64.to_le_bytes());
            }
            PreviewGlare::BoxBloomV1 {
                radius_px,
                threshold_linear,
                strength,
            } => {
                bytes.push(1);
                bytes.extend_from_slice(&radius_px.to_le_bytes());
                push_canonical_f64(&mut bytes, threshold_linear);
                push_canonical_f64(&mut bytes, strength);
            }
        }
        debug_assert_eq!(bytes.len(), CINEMATIC_COLOR_CONFIG_CANONICAL_BYTES);
        Ok(bytes)
    }

    /// Decode and validate the exact canonical version-one representation.
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, CinematicColorError> {
        if bytes.len() != CINEMATIC_COLOR_CONFIG_CANONICAL_BYTES {
            return Err(CinematicColorError::InvalidCanonicalConfig {
                field: "config",
                reason: "wrong fixed-width byte length",
            });
        }
        let mut reader = CanonicalReader::new(bytes);
        if reader.u16()? != CINEMATIC_COLOR_PIPELINE_VERSION {
            return Err(CinematicColorError::InvalidCanonicalConfig {
                field: "version",
                reason: "unsupported pipeline version",
            });
        }
        let working_space = match reader.u8()? {
            1 => CinematicWorkingSpace::LinearSrgbD65,
            _ => return Err(invalid_tag("working_space")),
        };
        let display_target = match reader.u8()? {
            1 => CinematicDisplayTarget::SrgbD65,
            _ => return Err(invalid_tag("display_target")),
        };
        let tone_curve = match reader.u8()? {
            1 => CinematicToneCurve::HableV1,
            2 => CinematicToneCurve::AcesFittedNarkowiczV1,
            _ => return Err(invalid_tag("tone_curve")),
        };
        let gamut_map = match reader.u8()? {
            1 => CinematicGamutMap::ClipV1,
            2 => CinematicGamutMap::RgbRatioScaleV1,
            _ => return Err(invalid_tag("gamut_map")),
        };
        let negative_policy = match reader.u8()? {
            1 => CinematicNegativePolicy::ClampToZeroCountedV1,
            _ => return Err(invalid_tag("negative_policy")),
        };
        let exposure_ev = reader.i32()?;
        let white_balance_gains = [reader.f64()?, reader.f64()?, reader.f64()?];
        let bit_depth = match reader.u8()? {
            8 => PreviewBitDepth::Eight,
            16 => PreviewBitDepth::Sixteen,
            _ => return Err(invalid_tag("bit_depth")),
        };
        let dither_tag = reader.u8()?;
        let dither_seed = reader.u64()?;
        let dither = match dither_tag {
            0 => PreviewDither::Disabled,
            1 => PreviewDither::UniformHalfLsbV1 { seed: dither_seed },
            _ => return Err(invalid_tag("dither")),
        };
        let glare_tag = reader.u8()?;
        let glare_radius = reader.u16()?;
        let glare_threshold = reader.f64()?;
        let glare_strength = reader.f64()?;
        let glare = match glare_tag {
            0 => PreviewGlare::Disabled,
            1 => PreviewGlare::BoxBloomV1 {
                radius_px: glare_radius,
                threshold_linear: glare_threshold,
                strength: glare_strength,
            },
            _ => return Err(invalid_tag("glare")),
        };
        let config = Self {
            working_space,
            display_target,
            tone_curve,
            gamut_map,
            negative_policy,
            exposure_ev,
            white_balance_gains,
            bit_depth,
            dither,
            glare,
        }
        .validate()?;
        if config.canonical_bytes()?.as_slice() != bytes {
            return Err(CinematicColorError::InvalidCanonicalConfig {
                field: "config",
                reason: "non-canonical reserved payload or floating-point zero",
            });
        }
        Ok(config)
    }
}

impl Default for CinematicColorConfig {
    fn default() -> Self {
        Self::reference_srgb_16()
    }
}

/// Caller-supplied limits checked before preview allocations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CinematicColorLimits {
    max_pixels: usize,
    max_working_bytes: usize,
}

impl CinematicColorLimits {
    /// Construct a bounded preview envelope.
    pub fn try_new(
        max_pixels: usize,
        max_working_bytes: usize,
    ) -> Result<Self, CinematicColorError> {
        if max_pixels == 0 || max_pixels > MAX_CINEMATIC_PREVIEW_PIXELS {
            return Err(CinematicColorError::InvalidLimits {
                field: "max_pixels",
            });
        }
        if max_working_bytes == 0 {
            return Err(CinematicColorError::InvalidLimits {
                field: "max_working_bytes",
            });
        }
        Ok(Self {
            max_pixels,
            max_working_bytes,
        })
    }

    /// Envelope sized for one 3840×2160 16-bit preview plus bloom scratch.
    #[must_use]
    pub const fn reference_4k() -> Self {
        Self {
            max_pixels: 3_840 * 2_160,
            max_working_bytes: 320 * 1024 * 1024,
        }
    }

    /// Admitted pixel ceiling.
    #[must_use]
    pub const fn max_pixels(self) -> usize {
        self.max_pixels
    }

    /// Admitted live working-byte ceiling.
    #[must_use]
    pub const fn max_working_bytes(self) -> usize {
        self.max_working_bytes
    }
}

/// Authority class of a cinematic color output.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum CinematicPreviewAuthority {
    /// Display-referred visualization derived from an immutable raw estimate.
    DisplayReferredDerivativeV1 = 1,
}

/// Interleaved RGB samples ready for the matching fs-img PNG writer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CinematicPreviewSamples {
    /// Eight-bit RGB samples.
    U8(Vec<u8>),
    /// Sixteen-bit RGB samples in host integers; `write_png16` emits big-endian.
    U16(Vec<u16>),
}

impl CinematicPreviewSamples {
    /// Number of interleaved channel samples.
    #[must_use]
    pub fn len(&self) -> usize {
        match self {
            Self::U8(samples) => samples.len(),
            Self::U16(samples) => samples.len(),
        }
    }

    /// Whether no channel samples are present.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Borrow eight-bit samples when this preview uses eight-bit precision.
    #[must_use]
    pub fn as_u8(&self) -> Option<&[u8]> {
        match self {
            Self::U8(samples) => Some(samples),
            Self::U16(_) => None,
        }
    }

    /// Borrow sixteen-bit samples when this preview uses sixteen-bit precision.
    #[must_use]
    pub fn as_u16(&self) -> Option<&[u16]> {
        match self {
            Self::U8(_) => None,
            Self::U16(samples) => Some(samples),
        }
    }

    /// Exact sample precision carried by this payload.
    #[must_use]
    pub const fn bit_depth(&self) -> PreviewBitDepth {
        match self {
            Self::U8(_) => PreviewBitDepth::Eight,
            Self::U16(_) => PreviewBitDepth::Sixteen,
        }
    }
}

/// Replay and diagnostic metadata for one display derivative.
#[derive(Debug, Clone, PartialEq)]
pub struct CinematicColorMetadata {
    canonical_config: Vec<u8>,
    negative_linear_channels: u64,
    over_range_linear_channels: u64,
    gamut_mapped_pixels: u64,
    glare_source_rgb_sum: f64,
    glare_added_rgb_sum: f64,
    required_working_bytes: usize,
}

impl CinematicColorMetadata {
    /// Exact versioned configuration bytes needed to replay the transform.
    #[must_use]
    pub fn canonical_config(&self) -> &[u8] {
        &self.canonical_config
    }

    /// Exposed/white-balanced channels below zero, which were handled visibly.
    #[must_use]
    pub const fn negative_linear_channels(&self) -> u64 {
        self.negative_linear_channels
    }

    /// Exposed/white-balanced channels above nominal display white (`1.0`).
    #[must_use]
    pub const fn over_range_linear_channels(&self) -> u64 {
        self.over_range_linear_channels
    }

    /// Pixels on which the chosen gamut operation changed a tone-mapped value.
    #[must_use]
    pub const fn gamut_mapped_pixels(&self) -> u64 {
        self.gamut_mapped_pixels
    }

    /// Sum of thresholded RGB bright-pass components before convolution.
    #[must_use]
    pub const fn glare_source_rgb_sum(&self) -> f64 {
        self.glare_source_rgb_sum
    }

    /// Sum of RGB components added by bloom after zero-boundary convolution.
    #[must_use]
    pub const fn glare_added_rgb_sum(&self) -> f64 {
        self.glare_added_rgb_sum
    }

    /// Exact output-plus-scratch bytes admitted before allocation.
    #[must_use]
    pub const fn required_working_bytes(&self) -> usize {
        self.required_working_bytes
    }
}

/// One complete display-referred cinematic preview.
#[derive(Debug, Clone, PartialEq)]
pub struct CinematicPreview {
    width: u32,
    height: u32,
    samples: CinematicPreviewSamples,
    metadata: CinematicColorMetadata,
}

impl CinematicPreview {
    /// Raster width.
    #[must_use]
    pub const fn width(&self) -> u32 {
        self.width
    }

    /// Raster height.
    #[must_use]
    pub const fn height(&self) -> u32 {
        self.height
    }

    /// Interleaved RGB preview samples.
    #[must_use]
    pub const fn samples(&self) -> &CinematicPreviewSamples {
        &self.samples
    }

    /// Replay and handling metadata.
    #[must_use]
    pub const fn metadata(&self) -> &CinematicColorMetadata {
        &self.metadata
    }

    /// This output is always a display derivative, never a raw estimate.
    #[must_use]
    pub const fn authority(&self) -> CinematicPreviewAuthority {
        CinematicPreviewAuthority::DisplayReferredDerivativeV1
    }
}

/// Stable classification of an invalid numeric input.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CinematicInputValueClass {
    /// IEEE NaN.
    Nan,
    /// Positive infinity.
    PositiveInfinity,
    /// Negative infinity.
    NegativeInfinity,
    /// A finite value exceeded the admitted scene-linear magnitude.
    MagnitudeTooLarge,
}

/// Fail-closed cinematic color admission or processing error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CinematicColorError {
    /// A configuration field is outside its frozen v1 domain.
    InvalidConfig {
        /// Stable field path.
        field: &'static str,
        /// Concise repair constraint.
        reason: &'static str,
    },
    /// Canonical bytes are malformed, non-canonical, or from another version.
    InvalidCanonicalConfig {
        /// Stable field path when one field can be identified.
        field: &'static str,
        /// Stable reason.
        reason: &'static str,
    },
    /// A caller-supplied processing ceiling is invalid.
    InvalidLimits {
        /// Stable field path.
        field: &'static str,
    },
    /// Raster dimensions are zero or overflow the address space.
    InvalidDimensions,
    /// One RGB plane disagrees with the raster dimensions.
    PlaneShape {
        /// Zero-based R/G/B plane index.
        channel: usize,
        /// Required samples.
        expected: usize,
        /// Supplied samples.
        got: usize,
    },
    /// Raster exceeds the admitted pixel envelope.
    PixelLimit {
        /// Required pixels.
        required: usize,
        /// Admitted pixels.
        available: usize,
    },
    /// Output plus algorithm scratch exceeds the admitted byte envelope.
    WorkingMemoryLimit {
        /// Required bytes.
        required: usize,
        /// Admitted bytes.
        available: usize,
    },
    /// First invalid component in canonical raster/channel order.
    InvalidLinearInput {
        /// Zero-based row-major pixel index.
        pixel: usize,
        /// Zero-based R/G/B channel index.
        channel: usize,
        /// Whether the raw or adjusted value failed.
        stage: &'static str,
        /// Stable numeric classification.
        class: CinematicInputValueClass,
    },
    /// A checked size calculation overflowed.
    ArithmeticOverflow {
        /// Stable calculation name.
        context: &'static str,
    },
    /// The allocator refused a pre-admitted output or scratch buffer.
    AllocationFailed {
        /// Stable buffer name.
        context: &'static str,
        /// Exact requested payload bytes.
        requested_bytes: usize,
    },
}

impl fmt::Display for CinematicColorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfig { field, reason } => {
                write!(f, "invalid cinematic color field `{field}`: {reason}")
            }
            Self::InvalidCanonicalConfig { field, reason } => {
                write!(
                    f,
                    "invalid cinematic color canonical field `{field}`: {reason}"
                )
            }
            Self::InvalidLimits { field } => {
                write!(f, "invalid cinematic color limit `{field}`")
            }
            Self::InvalidDimensions => write!(f, "cinematic preview dimensions are invalid"),
            Self::PlaneShape {
                channel,
                expected,
                got,
            } => write!(
                f,
                "cinematic RGB plane {channel}: expected {expected} samples, got {got}"
            ),
            Self::PixelLimit {
                required,
                available,
            } => write!(
                f,
                "cinematic preview needs {required} pixels, limit is {available}"
            ),
            Self::WorkingMemoryLimit {
                required,
                available,
            } => write!(
                f,
                "cinematic preview needs {required} working bytes, limit is {available}"
            ),
            Self::InvalidLinearInput {
                pixel,
                channel,
                stage,
                class,
            } => write!(
                f,
                "invalid {stage} linear input at pixel {pixel}, channel {channel}: {class:?}"
            ),
            Self::ArithmeticOverflow { context } => {
                write!(
                    f,
                    "cinematic color arithmetic overflow while computing {context}"
                )
            }
            Self::AllocationFailed {
                context,
                requested_bytes,
            } => write!(
                f,
                "cinematic color allocator refused {requested_bytes} bytes for {context}"
            ),
        }
    }
}

impl std::error::Error for CinematicColorError {}

/// Produce a deterministic display derivative from immutable planar linear RGB.
///
/// All validation, shape, pixel, and working-memory checks complete before the
/// output or bloom scratch is allocated. The three input planes are never
/// modified. Non-finite and excessively large scene-linear values refuse; the
/// returned metadata counts every finite negative and over-range channel.
pub fn transform_cinematic_preview(
    width: u32,
    height: u32,
    linear_rgb: [&[f32]; 3],
    config: CinematicColorConfig,
    limits: CinematicColorLimits,
) -> Result<CinematicPreview, CinematicColorError> {
    let config = config.validate()?;
    let canonical_config = config.canonical_bytes()?;
    let admission = admit_preview(width, height, linear_rgb, config, limits)?;
    let mut samples = allocate_samples(
        config.bit_depth,
        admission.sample_count,
        admission.output_bytes,
    )?;
    let (gamut_mapped_pixels, glare_added_rgb_sum) = match config.glare {
        PreviewGlare::Disabled => (
            render_without_glare(
                linear_rgb,
                config,
                admission.exposure_scale,
                admission.pixel_count,
                &mut samples,
            )?,
            0.0,
        ),
        PreviewGlare::BoxBloomV1 {
            radius_px,
            threshold_linear,
            strength,
        } => render_box_bloom(
            width,
            height,
            linear_rgb,
            config,
            admission,
            radius_px,
            threshold_linear,
            strength,
            &mut samples,
        )?,
    };

    Ok(CinematicPreview {
        width,
        height,
        samples,
        metadata: CinematicColorMetadata {
            canonical_config,
            negative_linear_channels: admission.stats.negative_linear_channels,
            over_range_linear_channels: admission.stats.over_range_linear_channels,
            gamut_mapped_pixels,
            glare_source_rgb_sum: admission.stats.glare_source_rgb_sum,
            glare_added_rgb_sum,
            required_working_bytes: admission.required_working_bytes,
        },
    })
}

#[derive(Debug, Clone, Copy)]
struct PreviewAdmission {
    pixel_count: usize,
    sample_count: usize,
    output_bytes: usize,
    glare_bytes: usize,
    required_working_bytes: usize,
    exposure_scale: f64,
    stats: InputStats,
}

#[derive(Debug, Clone, Copy, Default)]
struct InputStats {
    negative_linear_channels: u64,
    over_range_linear_channels: u64,
    glare_source_rgb_sum: f64,
}

fn admit_preview(
    width: u32,
    height: u32,
    linear_rgb: [&[f32]; 3],
    config: CinematicColorConfig,
    limits: CinematicColorLimits,
) -> Result<PreviewAdmission, CinematicColorError> {
    let pixel_count = checked_pixel_count(width, height)?;
    if pixel_count > limits.max_pixels || pixel_count > MAX_CINEMATIC_PREVIEW_PIXELS {
        return Err(CinematicColorError::PixelLimit {
            required: pixel_count,
            available: limits.max_pixels.min(MAX_CINEMATIC_PREVIEW_PIXELS),
        });
    }
    for (channel, plane) in linear_rgb.iter().enumerate() {
        if plane.len() != pixel_count {
            return Err(CinematicColorError::PlaneShape {
                channel,
                expected: pixel_count,
                got: plane.len(),
            });
        }
    }
    let sample_count = checked_mul(pixel_count, 3, "output sample count")?;
    let output_bytes = checked_mul(
        sample_count,
        config.bit_depth.bytes_per_channel(),
        "output byte count",
    )?;
    let glare_bytes = if matches!(config.glare, PreviewGlare::Disabled) {
        0
    } else {
        checked_mul(
            pixel_count,
            size_of::<[f64; 3]>(),
            "glare scratch byte count",
        )?
    };
    let required_working_bytes =
        checked_add(output_bytes, glare_bytes, "aggregate working byte count")?;
    if required_working_bytes > limits.max_working_bytes {
        return Err(CinematicColorError::WorkingMemoryLimit {
            required: required_working_bytes,
            available: limits.max_working_bytes,
        });
    }
    let exposure_scale = fs_math::det::powi(2.0, config.exposure_ev);
    let stats = scan_input_stats(linear_rgb, config, exposure_scale, pixel_count)?;
    Ok(PreviewAdmission {
        pixel_count,
        sample_count,
        output_bytes,
        glare_bytes,
        required_working_bytes,
        exposure_scale,
        stats,
    })
}

fn scan_input_stats(
    linear_rgb: [&[f32]; 3],
    config: CinematicColorConfig,
    exposure_scale: f64,
    pixel_count: usize,
) -> Result<InputStats, CinematicColorError> {
    let mut stats = InputStats::default();
    for pixel in 0..pixel_count {
        for value in adjusted_pixel(linear_rgb, pixel, config, exposure_scale)? {
            stats.negative_linear_channels += u64::from(value < 0.0);
            stats.over_range_linear_channels += u64::from(value > 1.0);
            if let PreviewGlare::BoxBloomV1 {
                threshold_linear, ..
            } = config.glare
            {
                stats.glare_source_rgb_sum +=
                    (handle_negative(value, config.negative_policy) - threshold_linear).max(0.0);
            }
        }
    }
    Ok(stats)
}

fn render_without_glare(
    linear_rgb: [&[f32]; 3],
    config: CinematicColorConfig,
    exposure_scale: f64,
    pixel_count: usize,
    samples: &mut CinematicPreviewSamples,
) -> Result<u64, CinematicColorError> {
    let mut gamut_mapped_pixels = 0_u64;
    for pixel in 0..pixel_count {
        let mut linear = adjusted_pixel(linear_rgb, pixel, config, exposure_scale)?;
        for value in &mut linear {
            *value = handle_negative(*value, config.negative_policy);
        }
        gamut_mapped_pixels += u64::from(write_encoded_pixel(samples, pixel, linear, config));
    }
    Ok(gamut_mapped_pixels)
}

#[allow(clippy::too_many_arguments)]
fn render_box_bloom(
    width: u32,
    height: u32,
    linear_rgb: [&[f32]; 3],
    config: CinematicColorConfig,
    admission: PreviewAdmission,
    radius_px: u16,
    threshold_linear: f64,
    strength: f64,
    samples: &mut CinematicPreviewSamples,
) -> Result<(u64, f64), CinematicColorError> {
    let width = usize::try_from(width).map_err(|_| CinematicColorError::ArithmeticOverflow {
        context: "width conversion",
    })?;
    let height = usize::try_from(height).map_err(|_| CinematicColorError::ArithmeticOverflow {
        context: "height conversion",
    })?;
    let radius = usize::from(radius_px);
    let diameter = radius
        .checked_mul(2)
        .and_then(|value| value.checked_add(1))
        .ok_or(CinematicColorError::ArithmeticOverflow {
            context: "glare kernel diameter",
        })?;
    let normalization = 1.0 / diameter as f64;
    let horizontal = horizontal_bright_pass(
        width,
        height,
        linear_rgb,
        config,
        admission,
        radius,
        threshold_linear,
        normalization,
    )?;
    render_vertical_bloom(
        width,
        height,
        linear_rgb,
        config,
        admission.exposure_scale,
        radius,
        normalization,
        strength,
        &horizontal,
        samples,
    )
}

#[allow(clippy::too_many_arguments)]
fn horizontal_bright_pass(
    width: usize,
    height: usize,
    linear_rgb: [&[f32]; 3],
    config: CinematicColorConfig,
    admission: PreviewAdmission,
    radius: usize,
    threshold: f64,
    normalization: f64,
) -> Result<Vec<[f64; 3]>, CinematicColorError> {
    let mut horizontal = Vec::new();
    horizontal
        .try_reserve_exact(admission.pixel_count)
        .map_err(|_| CinematicColorError::AllocationFailed {
            context: "glare scratch",
            requested_bytes: admission.glare_bytes,
        })?;
    horizontal.resize(admission.pixel_count, [0.0; 3]);
    for y in 0..height {
        let row = y * width;
        let mut sum = [0.0; 3];
        for x in 0..=radius.min(width - 1) {
            add_rgb(
                &mut sum,
                bright_pixel(
                    linear_rgb,
                    row + x,
                    config,
                    admission.exposure_scale,
                    threshold,
                ),
            );
        }
        for x in 0..width {
            horizontal[row + x] = nonnegative_rgb(scale_rgb(sum, normalization));
            if x >= radius {
                sub_rgb(
                    &mut sum,
                    bright_pixel(
                        linear_rgb,
                        row + x - radius,
                        config,
                        admission.exposure_scale,
                        threshold,
                    ),
                );
            }
            let add_x = x.saturating_add(radius).saturating_add(1);
            if add_x < width {
                add_rgb(
                    &mut sum,
                    bright_pixel(
                        linear_rgb,
                        row + add_x,
                        config,
                        admission.exposure_scale,
                        threshold,
                    ),
                );
            }
        }
    }
    Ok(horizontal)
}

#[allow(clippy::too_many_arguments)]
fn render_vertical_bloom(
    width: usize,
    height: usize,
    linear_rgb: [&[f32]; 3],
    config: CinematicColorConfig,
    exposure_scale: f64,
    radius: usize,
    normalization: f64,
    strength: f64,
    horizontal: &[[f64; 3]],
    samples: &mut CinematicPreviewSamples,
) -> Result<(u64, f64), CinematicColorError> {
    let mut gamut_mapped_pixels = 0_u64;
    let mut glare_added_rgb_sum = 0.0;
    for x in 0..width {
        let mut sum = [0.0; 3];
        for y in 0..=radius.min(height - 1) {
            add_rgb(&mut sum, horizontal[y * width + x]);
        }
        for y in 0..height {
            let pixel = y * width + x;
            let blurred = nonnegative_rgb(scale_rgb(sum, normalization * strength));
            glare_added_rgb_sum += blurred[0] + blurred[1] + blurred[2];
            let mut linear = adjusted_pixel(linear_rgb, pixel, config, exposure_scale)?;
            for channel in 0..3 {
                linear[channel] =
                    handle_negative(linear[channel], config.negative_policy) + blurred[channel];
            }
            gamut_mapped_pixels += u64::from(write_encoded_pixel(samples, pixel, linear, config));
            if y >= radius {
                sub_rgb(&mut sum, horizontal[(y - radius) * width + x]);
            }
            let add_y = y.saturating_add(radius).saturating_add(1);
            if add_y < height {
                add_rgb(&mut sum, horizontal[add_y * width + x]);
            }
        }
    }
    Ok((gamut_mapped_pixels, glare_added_rgb_sum))
}

fn adjusted_pixel(
    linear_rgb: [&[f32]; 3],
    pixel: usize,
    config: CinematicColorConfig,
    exposure_scale: f64,
) -> Result<[f64; 3], CinematicColorError> {
    let mut adjusted = [0.0; 3];
    for channel in 0..3 {
        let raw = f64::from(linear_rgb[channel][pixel]);
        if !raw.is_finite() {
            return Err(invalid_value(pixel, channel, "raw", raw));
        }
        let value = (raw * config.white_balance_gains[channel]) * exposure_scale;
        if !value.is_finite() {
            return Err(invalid_value(pixel, channel, "adjusted", value));
        }
        if value.abs() > MAX_LINEAR_COMPONENT {
            return Err(CinematicColorError::InvalidLinearInput {
                pixel,
                channel,
                stage: "adjusted",
                class: CinematicInputValueClass::MagnitudeTooLarge,
            });
        }
        adjusted[channel] = value;
    }
    Ok(adjusted)
}

fn allocate_samples(
    bit_depth: PreviewBitDepth,
    sample_count: usize,
    output_bytes: usize,
) -> Result<CinematicPreviewSamples, CinematicColorError> {
    match bit_depth {
        PreviewBitDepth::Eight => {
            let mut values = Vec::new();
            values.try_reserve_exact(sample_count).map_err(|_| {
                CinematicColorError::AllocationFailed {
                    context: "8-bit preview output",
                    requested_bytes: output_bytes,
                }
            })?;
            values.resize(sample_count, 0_u8);
            Ok(CinematicPreviewSamples::U8(values))
        }
        PreviewBitDepth::Sixteen => {
            let mut values = Vec::new();
            values.try_reserve_exact(sample_count).map_err(|_| {
                CinematicColorError::AllocationFailed {
                    context: "16-bit preview output",
                    requested_bytes: output_bytes,
                }
            })?;
            values.resize(sample_count, 0_u16);
            Ok(CinematicPreviewSamples::U16(values))
        }
    }
}

fn bright_pixel(
    linear_rgb: [&[f32]; 3],
    pixel: usize,
    config: CinematicColorConfig,
    exposure_scale: f64,
    threshold: f64,
) -> [f64; 3] {
    let mut bright = [0.0; 3];
    for channel in 0..3 {
        let value = (f64::from(linear_rgb[channel][pixel]) * config.white_balance_gains[channel])
            * exposure_scale;
        bright[channel] = (handle_negative(value, config.negative_policy) - threshold).max(0.0);
    }
    bright
}

fn write_encoded_pixel(
    samples: &mut CinematicPreviewSamples,
    pixel: usize,
    linear: [f64; 3],
    config: CinematicColorConfig,
) -> bool {
    let (display_linear, gamut_changed) = tone_and_gamut(linear, config);
    let base = pixel * 3;
    match samples {
        CinematicPreviewSamples::U8(values) => {
            for channel in 0..3 {
                values[base + channel] = quantize(
                    srgb_encode(display_linear[channel]),
                    u64::from(u8::MAX),
                    config.dither,
                    pixel,
                    channel,
                ) as u8;
            }
        }
        CinematicPreviewSamples::U16(values) => {
            for channel in 0..3 {
                values[base + channel] = quantize(
                    srgb_encode(display_linear[channel]),
                    u64::from(u16::MAX),
                    config.dither,
                    pixel,
                    channel,
                ) as u16;
            }
        }
    }
    gamut_changed
}

fn tone_and_gamut(linear: [f64; 3], config: CinematicColorConfig) -> ([f64; 3], bool) {
    let mut mapped = linear.map(|value| match config.tone_curve {
        CinematicToneCurve::HableV1 => {
            hable_filmic_unclamped(handle_negative(value, config.negative_policy))
        }
        CinematicToneCurve::AcesFittedNarkowiczV1 => {
            aces_fitted(handle_negative(value, config.negative_policy))
        }
    });
    let changed = mapped.iter().any(|value| !(0.0..=1.0).contains(value));
    match config.gamut_map {
        CinematicGamutMap::ClipV1 => {
            for value in &mut mapped {
                *value = value.clamp(0.0, 1.0);
            }
        }
        CinematicGamutMap::RgbRatioScaleV1 => {
            for value in &mut mapped {
                *value = value.max(0.0);
            }
            let maximum = mapped[0].max(mapped[1]).max(mapped[2]);
            if maximum > 1.0 {
                for value in &mut mapped {
                    *value /= maximum;
                }
            }
        }
    }
    (mapped, changed)
}

fn handle_negative(value: f64, policy: CinematicNegativePolicy) -> f64 {
    match policy {
        CinematicNegativePolicy::ClampToZeroCountedV1 => value.max(0.0),
    }
}

fn aces_fitted(x: f64) -> f64 {
    const A: f64 = 2.51;
    const B: f64 = 0.03;
    const C: f64 = 2.43;
    const D: f64 = 0.59;
    const E: f64 = 0.14;
    (x * (A * x + B)) / (x * (C * x + D) + E)
}

fn quantize(
    encoded: f64,
    maximum: u64,
    dither: PreviewDither,
    pixel: usize,
    channel: usize,
) -> u64 {
    let encoded = encoded.clamp(0.0, 1.0);
    if encoded <= 0.0 {
        return 0;
    }
    if encoded >= 1.0 {
        return maximum;
    }
    let noise = match dither {
        PreviewDither::Disabled => 0.0,
        PreviewDither::UniformHalfLsbV1 { seed } => dither_lsb(seed, pixel, channel),
    };
    (encoded * maximum as f64 + noise + 0.5).clamp(0.0, maximum as f64) as u64
}

fn dither_lsb(seed: u64, pixel: usize, channel: usize) -> f64 {
    let pixel_key = u64::try_from(pixel).unwrap_or(u64::MAX);
    let channel_key = u64::try_from(channel).unwrap_or(u64::MAX);
    let mut value = seed
        ^ pixel_key.wrapping_mul(0x9e37_79b9_7f4a_7c15)
        ^ channel_key.wrapping_mul(0xd1b5_4a32_d192_ed03);
    value = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^= value >> 31;
    (value >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0) - 0.5
}

fn checked_pixel_count(width: u32, height: u32) -> Result<usize, CinematicColorError> {
    if width == 0 || height == 0 {
        return Err(CinematicColorError::InvalidDimensions);
    }
    let pixels = u64::from(width).checked_mul(u64::from(height)).ok_or(
        CinematicColorError::ArithmeticOverflow {
            context: "pixel count",
        },
    )?;
    usize::try_from(pixels).map_err(|_| CinematicColorError::ArithmeticOverflow {
        context: "pixel count conversion",
    })
}

fn checked_mul(
    left: usize,
    right: usize,
    context: &'static str,
) -> Result<usize, CinematicColorError> {
    left.checked_mul(right)
        .ok_or(CinematicColorError::ArithmeticOverflow { context })
}

fn checked_add(
    left: usize,
    right: usize,
    context: &'static str,
) -> Result<usize, CinematicColorError> {
    left.checked_add(right)
        .ok_or(CinematicColorError::ArithmeticOverflow { context })
}

fn invalid_value(
    pixel: usize,
    channel: usize,
    stage: &'static str,
    value: f64,
) -> CinematicColorError {
    let class = if value.is_nan() {
        CinematicInputValueClass::Nan
    } else if value.is_sign_positive() {
        CinematicInputValueClass::PositiveInfinity
    } else {
        CinematicInputValueClass::NegativeInfinity
    };
    CinematicColorError::InvalidLinearInput {
        pixel,
        channel,
        stage,
        class,
    }
}

fn add_rgb(sum: &mut [f64; 3], value: [f64; 3]) {
    for channel in 0..3 {
        sum[channel] += value[channel];
    }
}

fn sub_rgb(sum: &mut [f64; 3], value: [f64; 3]) {
    for channel in 0..3 {
        sum[channel] -= value[channel];
    }
}

fn scale_rgb(value: [f64; 3], scale: f64) -> [f64; 3] {
    [value[0] * scale, value[1] * scale, value[2] * scale]
}

fn nonnegative_rgb(value: [f64; 3]) -> [f64; 3] {
    [value[0].max(0.0), value[1].max(0.0), value[2].max(0.0)]
}

fn push_canonical_f64(bytes: &mut Vec<u8>, value: f64) {
    let bits = if value == 0.0 { 0 } else { value.to_bits() };
    bytes.extend_from_slice(&bits.to_le_bytes());
}

fn invalid_tag(field: &'static str) -> CinematicColorError {
    CinematicColorError::InvalidCanonicalConfig {
        field,
        reason: "unknown canonical enum tag",
    }
}

struct CanonicalReader<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> CanonicalReader<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    fn take<const N: usize>(&mut self) -> Result<[u8; N], CinematicColorError> {
        let end =
            self.offset
                .checked_add(N)
                .ok_or(CinematicColorError::InvalidCanonicalConfig {
                    field: "config",
                    reason: "field offset overflow",
                })?;
        let source = self.bytes.get(self.offset..end).ok_or(
            CinematicColorError::InvalidCanonicalConfig {
                field: "config",
                reason: "truncated field",
            },
        )?;
        let mut value = [0_u8; N];
        value.copy_from_slice(source);
        self.offset = end;
        Ok(value)
    }

    fn u8(&mut self) -> Result<u8, CinematicColorError> {
        Ok(self.take::<1>()?[0])
    }

    fn u16(&mut self) -> Result<u16, CinematicColorError> {
        Ok(u16::from_le_bytes(self.take()?))
    }

    fn u64(&mut self) -> Result<u64, CinematicColorError> {
        Ok(u64::from_le_bytes(self.take()?))
    }

    fn i32(&mut self) -> Result<i32, CinematicColorError> {
        Ok(i32::from_le_bytes(self.take()?))
    }

    fn f64(&mut self) -> Result<f64, CinematicColorError> {
        Ok(f64::from_bits(self.u64()?))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn limits(pixels: usize) -> CinematicColorLimits {
        CinematicColorLimits::try_new(pixels, 16 * 1024 * 1024).unwrap()
    }

    fn planes(values: &[[f32; 3]]) -> [Vec<f32>; 3] {
        let mut out = [Vec::new(), Vec::new(), Vec::new()];
        for value in values {
            for channel in 0..3 {
                out[channel].push(value[channel]);
            }
        }
        out
    }

    fn slices(planes: &[Vec<f32>; 3]) -> [&[f32]; 3] {
        [&planes[0], &planes[1], &planes[2]]
    }

    #[test]
    fn g0_config_canonical_round_trip_is_exact_and_closed() {
        let config = CinematicColorConfig {
            glare: PreviewGlare::BoxBloomV1 {
                radius_px: 17,
                threshold_linear: 1.25,
                strength: 0.2,
            },
            ..CinematicColorConfig::reference_srgb_16()
        };
        let bytes = config.canonical_bytes().unwrap();
        assert_eq!(bytes.len(), CINEMATIC_COLOR_CONFIG_CANONICAL_BYTES);
        assert_eq!(
            CinematicColorConfig::from_canonical_bytes(&bytes).unwrap(),
            config
        );

        let mut bad_tag = bytes.clone();
        bad_tag[4] = 99;
        assert!(CinematicColorConfig::from_canonical_bytes(&bad_tag).is_err());
        let mut trailing = bytes.clone();
        trailing.push(0);
        assert!(CinematicColorConfig::from_canonical_bytes(&trailing).is_err());
        for prefix in 0..bytes.len() {
            assert!(CinematicColorConfig::from_canonical_bytes(&bytes[..prefix]).is_err());
        }

        let disabled = CinematicColorConfig {
            dither: PreviewDither::Disabled,
            glare: PreviewGlare::Disabled,
            ..config
        };
        let mut noncanonical = disabled.canonical_bytes().unwrap();
        noncanonical[37] = 1;
        assert!(CinematicColorConfig::from_canonical_bytes(&noncanonical).is_err());
    }

    #[test]
    fn g3_every_semantic_config_field_moves_canonical_bytes() {
        let base = CinematicColorConfig::reference_srgb_16();
        let baseline = base.canonical_bytes().unwrap();
        let mutations = [
            CinematicColorConfig {
                tone_curve: CinematicToneCurve::HableV1,
                ..base
            },
            CinematicColorConfig {
                gamut_map: CinematicGamutMap::ClipV1,
                ..base
            },
            CinematicColorConfig {
                exposure_ev: 1,
                ..base
            },
            CinematicColorConfig {
                white_balance_gains: [1.1, 1.0, 1.0],
                ..base
            },
            CinematicColorConfig {
                bit_depth: PreviewBitDepth::Eight,
                ..base
            },
            CinematicColorConfig {
                dither: PreviewDither::Disabled,
                ..base
            },
            CinematicColorConfig {
                dither: PreviewDither::UniformHalfLsbV1 { seed: 1 },
                ..base
            },
            CinematicColorConfig {
                glare: PreviewGlare::BoxBloomV1 {
                    radius_px: 1,
                    threshold_linear: 1.0,
                    strength: 0.1,
                },
                ..base
            },
        ];
        for mutation in mutations {
            assert_ne!(mutation.canonical_bytes().unwrap(), baseline);
        }
    }

    #[test]
    fn g0_frozen_curves_are_monotone_and_have_known_anchors() {
        assert_eq!(aces_fitted(0.0).to_bits(), 0.0_f64.to_bits());
        assert!((aces_fitted(1.0) - 0.803_797_468_354_430_2).abs() < 1.0e-15);
        let mut last_hable = -1.0;
        let mut last_aces = -1.0;
        for sample in 0..=10_000 {
            let x = f64::from(sample) * 0.001;
            let hable = hable_filmic_unclamped(x);
            let aces = aces_fitted(x);
            assert!(hable >= last_hable);
            assert!(aces >= last_aces);
            last_hable = hable;
            last_aces = aces;
        }
    }

    #[test]
    fn g0_gray_neutrality_exposure_and_integer_endpoints_hold() {
        let input = planes(&[[0.0; 3], [0.18; 3], [1.0e6_f32; 3]]);
        let config = CinematicColorConfig {
            tone_curve: CinematicToneCurve::HableV1,
            gamut_map: CinematicGamutMap::ClipV1,
            bit_depth: PreviewBitDepth::Eight,
            dither: PreviewDither::Disabled,
            ..CinematicColorConfig::reference_srgb_16()
        };
        let base = transform_cinematic_preview(3, 1, slices(&input), config, limits(3)).unwrap();
        let CinematicPreviewSamples::U8(base_samples) = base.samples() else {
            panic!("expected u8 preview");
        };
        assert_eq!(&base_samples[0..3], &[0, 0, 0]);
        assert_eq!(base_samples[3], base_samples[4]);
        assert_eq!(base_samples[4], base_samples[5]);
        assert_eq!(&base_samples[6..9], &[255, 255, 255]);

        let brighter = transform_cinematic_preview(
            3,
            1,
            slices(&input),
            CinematicColorConfig {
                exposure_ev: 1,
                ..config
            },
            limits(3),
        )
        .unwrap();
        let CinematicPreviewSamples::U8(brighter_samples) = brighter.samples() else {
            panic!("expected u8 preview");
        };
        assert!(brighter_samples[3] > base_samples[3]);
    }

    #[test]
    fn g0_white_balance_and_gamut_policies_have_visible_frozen_effects() {
        let gray = planes(&[[0.18_f32; 3]]);
        let balanced = transform_cinematic_preview(
            1,
            1,
            slices(&gray),
            CinematicColorConfig {
                white_balance_gains: [2.0, 1.0, 0.5],
                bit_depth: PreviewBitDepth::Sixteen,
                dither: PreviewDither::Disabled,
                ..CinematicColorConfig::default()
            },
            limits(1),
        )
        .unwrap();
        let CinematicPreviewSamples::U16(samples) = balanced.samples() else {
            panic!("expected u16 preview");
        };
        assert!(samples[0] > samples[1] && samples[1] > samples[2]);

        let extreme = planes(&[[100.0, 1.0, 0.1]]);
        let base = CinematicColorConfig {
            tone_curve: CinematicToneCurve::HableV1,
            bit_depth: PreviewBitDepth::Sixteen,
            dither: PreviewDither::Disabled,
            ..CinematicColorConfig::default()
        };
        let clipped = transform_cinematic_preview(
            1,
            1,
            slices(&extreme),
            CinematicColorConfig {
                gamut_map: CinematicGamutMap::ClipV1,
                ..base
            },
            limits(1),
        )
        .unwrap();
        let ratio_scaled = transform_cinematic_preview(
            1,
            1,
            slices(&extreme),
            CinematicColorConfig {
                gamut_map: CinematicGamutMap::RgbRatioScaleV1,
                ..base
            },
            limits(1),
        )
        .unwrap();
        assert_ne!(clipped.samples(), ratio_scaled.samples());
        assert_eq!(clipped.metadata().gamut_mapped_pixels(), 1);
        assert_eq!(ratio_scaled.metadata().gamut_mapped_pixels(), 1);
    }

    #[test]
    fn g0_hdr_highlight_steps_remain_ordered_before_the_display_ceiling() {
        let input = planes(&[[1.0_f32; 3], [2.0_f32; 3], [4.0_f32; 3]]);
        let preview = transform_cinematic_preview(
            3,
            1,
            slices(&input),
            CinematicColorConfig {
                bit_depth: PreviewBitDepth::Sixteen,
                dither: PreviewDither::Disabled,
                ..CinematicColorConfig::default()
            },
            limits(3),
        )
        .unwrap();
        let samples = preview.samples().as_u16().unwrap();
        assert!(samples[0] < samples[3]);
        assert!(samples[3] < samples[6]);
        assert!(samples[6] < u16::MAX);
    }

    #[test]
    fn g0_nonfinite_shape_pixel_and_memory_fail_before_work() {
        let invalid_glare = CinematicColorConfig {
            glare: PreviewGlare::BoxBloomV1 {
                radius_px: 1,
                threshold_linear: 1.0,
                strength: 0.0,
            },
            ..CinematicColorConfig::default()
        };
        assert!(matches!(
            invalid_glare.validate(),
            Err(CinematicColorError::InvalidConfig {
                field: "glare.strength",
                ..
            })
        ));

        let malformed = [vec![0.0_f32], Vec::new(), vec![0.0_f32]];
        assert_eq!(
            transform_cinematic_preview(
                1,
                1,
                slices(&malformed),
                CinematicColorConfig::default(),
                limits(1),
            ),
            Err(CinematicColorError::PlaneShape {
                channel: 1,
                expected: 1,
                got: 0,
            })
        );

        let finite = planes(&[[0.0; 3], [0.0; 3]]);
        assert!(matches!(
            transform_cinematic_preview(
                2,
                1,
                slices(&finite),
                CinematicColorConfig::default(),
                limits(1),
            ),
            Err(CinematicColorError::PixelLimit { .. })
        ));
        assert!(matches!(
            transform_cinematic_preview(
                2,
                1,
                slices(&finite),
                CinematicColorConfig::default(),
                CinematicColorLimits::try_new(2, 1).unwrap(),
            ),
            Err(CinematicColorError::WorkingMemoryLimit { .. })
        ));

        for (value, class) in [
            (f32::NAN, CinematicInputValueClass::Nan),
            (f32::INFINITY, CinematicInputValueClass::PositiveInfinity),
            (
                f32::NEG_INFINITY,
                CinematicInputValueClass::NegativeInfinity,
            ),
        ] {
            let invalid = [vec![0.0_f32], vec![value], vec![0.0_f32]];
            assert_eq!(
                transform_cinematic_preview(
                    1,
                    1,
                    slices(&invalid),
                    CinematicColorConfig::default(),
                    limits(1),
                ),
                Err(CinematicColorError::InvalidLinearInput {
                    pixel: 0,
                    channel: 1,
                    stage: "raw",
                    class,
                })
            );
        }
    }

    #[test]
    fn g3_negative_and_over_range_handling_is_visible_and_input_is_immutable() {
        let input = planes(&[[-0.25, 0.5, 2.0], [4.0, -1.0, 0.25]]);
        let before = input.clone();
        let preview = transform_cinematic_preview(
            2,
            1,
            slices(&input),
            CinematicColorConfig::default(),
            limits(2),
        )
        .unwrap();
        assert_eq!(input, before);
        assert_eq!(preview.metadata().negative_linear_channels(), 2);
        assert_eq!(preview.metadata().over_range_linear_channels(), 2);
        assert_eq!(
            preview.authority(),
            CinematicPreviewAuthority::DisplayReferredDerivativeV1
        );
        assert_eq!(
            CinematicColorConfig::from_canonical_bytes(preview.metadata().canonical_config())
                .unwrap(),
            CinematicColorConfig::default()
        );
    }

    #[test]
    fn g5_dither_is_replay_exact_seed_sensitive_and_preserves_endpoints() {
        let mut values = vec![[0.18_f32; 3]; 512];
        values[0] = [0.0; 3];
        values[1] = [1.0e6_f32; 3];
        let input = planes(&values);
        let base = CinematicColorConfig {
            bit_depth: PreviewBitDepth::Eight,
            dither: PreviewDither::UniformHalfLsbV1 { seed: 17 },
            ..CinematicColorConfig::default()
        };
        let first = transform_cinematic_preview(512, 1, slices(&input), base, limits(512)).unwrap();
        let replay =
            transform_cinematic_preview(512, 1, slices(&input), base, limits(512)).unwrap();
        assert_eq!(first, replay);
        let changed = transform_cinematic_preview(
            512,
            1,
            slices(&input),
            CinematicColorConfig {
                dither: PreviewDither::UniformHalfLsbV1 { seed: 18 },
                ..base
            },
            limits(512),
        )
        .unwrap();
        assert_ne!(first.samples(), changed.samples());
        let CinematicPreviewSamples::U8(samples) = first.samples() else {
            panic!("expected u8 preview");
        };
        assert_eq!(&samples[0..3], &[0, 0, 0]);
        assert_eq!(&samples[3..6], &[255, 255, 255]);
    }

    #[test]
    fn g0_half_lsb_dither_is_centered_over_a_keyed_sequence() {
        let mut total = 0_u64;
        let samples = 16_384_usize;
        for pixel in 0..samples {
            total += quantize(
                0.5,
                u64::from(u8::MAX),
                PreviewDither::UniformHalfLsbV1 { seed: 0x5eed },
                pixel,
                0,
            );
        }
        let mean = total as f64 / samples as f64;
        assert!((mean - 127.5).abs() < 0.02, "mean quantized value {mean}");
    }

    #[test]
    fn g3_box_bloom_is_local_thresholded_and_explicitly_accounted() {
        let mut values = vec![[0.0_f32; 3]; 49];
        values[3 * 7] = [8.0, 8.0, 8.0];
        let input = planes(&values);
        let config = CinematicColorConfig {
            bit_depth: PreviewBitDepth::Eight,
            dither: PreviewDither::Disabled,
            glare: PreviewGlare::BoxBloomV1 {
                radius_px: 1,
                threshold_linear: 1.0,
                strength: 0.5,
            },
            ..CinematicColorConfig::default()
        };
        let preview =
            transform_cinematic_preview(7, 7, slices(&input), config, limits(49)).unwrap();
        let CinematicPreviewSamples::U8(samples) = preview.samples() else {
            panic!("expected u8 preview");
        };
        let neighbor = (3 * 7 + 1) * 3;
        let far_edge = (3 * 7 + 6) * 3;
        assert!(
            samples[neighbor] > 0,
            "bloom must reach its admitted radius"
        );
        assert_eq!(samples[far_edge], 0, "zero-boundary bloom must not wrap");
        assert!(preview.metadata().glare_source_rgb_sum() > 0.0);
        assert!(preview.metadata().glare_added_rgb_sum() > 0.0);

        let below_threshold = planes(&vec![[0.5_f32; 3]; 49]);
        let with_bloom =
            transform_cinematic_preview(7, 7, slices(&below_threshold), config, limits(49))
                .unwrap();
        let without_bloom = transform_cinematic_preview(
            7,
            7,
            slices(&below_threshold),
            CinematicColorConfig {
                glare: PreviewGlare::Disabled,
                ..config
            },
            limits(49),
        )
        .unwrap();
        assert_eq!(with_bloom.samples(), without_bloom.samples());
        assert_eq!(with_bloom.metadata().glare_source_rgb_sum(), 0.0);
        assert_eq!(with_bloom.metadata().glare_added_rgb_sum(), 0.0);
    }

    #[test]
    fn g0_constant_bright_field_stays_constant_away_from_zero_boundaries() {
        let input = planes(&vec![[2.0_f32; 3]; 49]);
        let config = CinematicColorConfig {
            bit_depth: PreviewBitDepth::Sixteen,
            dither: PreviewDither::Disabled,
            glare: PreviewGlare::BoxBloomV1 {
                radius_px: 1,
                threshold_linear: 1.0,
                strength: 0.25,
            },
            ..CinematicColorConfig::default()
        };
        let preview =
            transform_cinematic_preview(7, 7, slices(&input), config, limits(49)).unwrap();
        let CinematicPreviewSamples::U16(samples) = preview.samples() else {
            panic!("expected u16 preview");
        };
        let center = (3 * 7 + 3) * 3;
        let adjacent = (3 * 7 + 4) * 3;
        assert_eq!(
            &samples[center..center + 3],
            &samples[adjacent..adjacent + 3]
        );
        assert_eq!(
            preview.metadata().required_working_bytes(),
            49 * (3 * size_of::<u16>() + size_of::<[f64; 3]>())
        );
    }

    #[test]
    fn g0_preview_samples_round_trip_through_matching_png_depths() {
        use crate::{PngColor, read_png, write_png8, write_png16};

        let input = planes(&[[0.0, 0.18, 1.0], [2.0, 0.5, 0.01]]);
        for depth in [PreviewBitDepth::Eight, PreviewBitDepth::Sixteen] {
            let preview = transform_cinematic_preview(
                2,
                1,
                slices(&input),
                CinematicColorConfig {
                    bit_depth: depth,
                    dither: PreviewDither::Disabled,
                    ..CinematicColorConfig::default()
                },
                limits(2),
            )
            .unwrap();
            assert_eq!(preview.samples().bit_depth(), depth);
            match preview.samples() {
                CinematicPreviewSamples::U8(samples) => {
                    let encoded = write_png8(2, 1, PngColor::Rgb, samples).unwrap();
                    let decoded = read_png(&encoded).unwrap();
                    assert_eq!(decoded.depth, 8);
                    assert_eq!(decoded.bytes.as_slice(), samples.as_slice());
                }
                CinematicPreviewSamples::U16(samples) => {
                    let encoded = write_png16(2, 1, PngColor::Rgb, samples).unwrap();
                    let decoded = read_png(&encoded).unwrap();
                    assert_eq!(decoded.depth, 16);
                    assert_eq!(decoded.samples16().as_slice(), samples.as_slice());
                }
            }
        }
    }

    #[test]
    fn g0_reference_4k_envelope_covers_worst_preview_payload_and_scratch() {
        let pixels = 3_840_usize * 2_160;
        let required = pixels * (3 * size_of::<u16>() + size_of::<[f64; 3]>());
        let limits = CinematicColorLimits::reference_4k();
        assert_eq!(limits.max_pixels(), pixels);
        assert!(required <= limits.max_working_bytes());
    }
}
