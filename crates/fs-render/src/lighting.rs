//! Deterministic multi-emitter and environment-light sampling.
//!
//! Rectangular emitters are ordered by content rather than caller insertion
//! order. Their selection weights are emitted luminance times the exact solid
//! angle subtended at the shading point. A canonical Y-up latitude/longitude
//! environment uses texel luminance times exact texel solid angle. The same
//! mixture probabilities are exposed for NEE samples and BSDF-hit MIS.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::OnceLock;

use fs_blake3::{ContentHash, DomainHasher, hash_domain};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::det;

use crate::spectral::{LiftedSpectrum, lift_rgb};

const PI: f64 = core::f64::consts::PI;
const TAU: f64 = core::f64::consts::TAU;
const RECTANGLE_IDENTITY_DOMAIN: &str = "org.frankensim.fs-render.rect-light.v1";
const ENVIRONMENT_SEMANTIC_DOMAIN: &str = "org.frankensim.fs-render.environment-map.v1";
const ENVIRONMENT_NATIVE_SOURCE_DOMAIN: &str =
    "org.frankensim.fs-render.environment-native-source.v1";
const ENVIRONMENT_EXR_SOURCE_DOMAIN: &str = "org.frankensim.fs-render.environment-exr-source.v1";
const ENVIRONMENT_PROVENANCE_DOMAIN: &str = "org.frankensim.fs-render.environment-provenance.v1";
const ENVIRONMENT_SEMANTICS_VERSION: u32 = 1;
const ENVIRONMENT_IMPORTER_VERSION: u32 = 1;
const RECTANGULARITY_TOLERANCE: f64 = 1.0e-10;
const PLANE_DISTANCE_TOLERANCE: f64 = 1.0e-12;
const LEGACY_RECTANGLE_COSINE_CUTOFF: f64 = 1.0e-9;
const MULTI_LIGHT_RECTANGLE_COSINE_CUTOFF: f64 = 1.0e-12;
const SPECTRAL_LUT_EDGE: usize = 9;

/// Admission failures for rectangular and environment lighting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LightingError {
    /// The caller's execution authority requested cancellation while scene
    /// emitters were being admitted.
    Cancelled,
    /// A rectangular emitter had non-finite, zero, non-orthogonal, or
    /// non-positive-radiance metadata.
    InvalidRectangle {
        /// Input light index.
        light_index: usize,
    },
    /// Two rectangles had the same content identity and would double-count one
    /// physical source.
    DuplicateRectangle {
        /// First input index.
        first: usize,
        /// Second input index.
        second: usize,
    },
    /// Two NEE rectangles named the same emissive primitive.
    DuplicatePrimitive {
        /// Repeated primitive index.
        primitive_index: usize,
    },
    /// Environment dimensions were zero, overflowed, or disagreed with the
    /// supplied pixel count.
    InvalidEnvironmentDimensions,
    /// A linear-radiance texel was negative or non-finite.
    InvalidEnvironmentPixel {
        /// Row-major pixel index.
        pixel_index: usize,
    },
    /// Environment rotation was non-finite.
    InvalidEnvironmentRotation,
    /// The EXR was outside the in-house subset or did not contain exactly one
    /// usable R, G, and B plane.
    UnsupportedEnvironmentExr,
    /// Allocation for deterministic sampling metadata failed.
    EnvironmentTooLarge,
    /// The admitted scene had neither a positive rectangle nor a non-black
    /// environment.
    NoFiniteEmitter,
}

impl core::fmt::Display for LightingError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled => formatter.write_str("lighting admission cancelled"),
            Self::InvalidRectangle { light_index } => {
                write!(
                    formatter,
                    "invalid rectangular emitter at index {light_index}"
                )
            }
            Self::DuplicateRectangle { first, second } => write!(
                formatter,
                "rectangular emitters {first} and {second} have one content identity"
            ),
            Self::DuplicatePrimitive { primitive_index } => write!(
                formatter,
                "multiple rectangular emitters name primitive {primitive_index}"
            ),
            Self::InvalidEnvironmentDimensions => {
                formatter.write_str("invalid environment-map dimensions or pixel count")
            }
            Self::InvalidEnvironmentPixel { pixel_index } => write!(
                formatter,
                "environment pixel {pixel_index} is negative or non-finite"
            ),
            Self::InvalidEnvironmentRotation => {
                formatter.write_str("environment rotation must be finite")
            }
            Self::UnsupportedEnvironmentExr => formatter.write_str(
                "environment EXR must use the supported scanline subset with R/G/B planes",
            ),
            Self::EnvironmentTooLarge => {
                formatter.write_str("environment sampling metadata allocation failed")
            }
            Self::NoFiniteEmitter => formatter.write_str("scene has no positive finite emitter"),
        }
    }
}

impl core::error::Error for LightingError {}

/// One static rectangular area emitter. The same rectangle must be present as
/// emissive scene geometry at `prim` so BSDF-hit and NEE paths share a target.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RectLight {
    /// One world-space corner.
    pub corner: Point3,
    /// First edge vector.
    pub edge_u: Vec3,
    /// Second, orthogonal edge vector.
    pub edge_v: Vec3,
    /// Index of the matching emissive primitive.
    pub prim: usize,
    /// Spectral emitted radiance and positive scale.
    pub emission: (LiftedSpectrum, f64),
}

impl RectLight {
    /// Rectangle area in square world units.
    #[must_use]
    pub fn area(&self) -> f64 {
        cross(self.edge_u, self.edge_v).norm()
    }

    /// Unit normal from `edge_u × edge_v`. Call only after scene admission.
    #[must_use]
    pub fn normal(&self) -> Vec3 {
        let normal = cross(self.edge_u, self.edge_v);
        normal.scale(1.0 / normal.norm())
    }

    /// Content identity of geometry and emission. The primitive index is
    /// intentionally excluded because insertion order is not physical content.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        let mut hasher = DomainHasher::new(RECTANGLE_IDENTITY_DOMAIN);
        for value in [
            self.corner.x,
            self.corner.y,
            self.corner.z,
            self.edge_u.x,
            self.edge_u.y,
            self.edge_u.z,
            self.edge_v.x,
            self.edge_v.y,
            self.edge_v.z,
            self.emission.0.c[0],
            self.emission.0.c[1],
            self.emission.0.c[2],
            self.emission.1,
        ] {
            hasher.update(&canonical_f64_bits(value).to_le_bytes());
        }
        hasher.finalize()
    }

    fn is_valid(&self) -> bool {
        let values = [
            self.corner.x,
            self.corner.y,
            self.corner.z,
            self.edge_u.x,
            self.edge_u.y,
            self.edge_u.z,
            self.edge_v.x,
            self.edge_v.y,
            self.edge_v.z,
            self.emission.0.c[0],
            self.emission.0.c[1],
            self.emission.0.c[2],
            self.emission.1,
        ];
        if values.iter().any(|value| !value.is_finite()) || self.emission.1 <= 0.0 {
            return false;
        }
        let u_norm = self.edge_u.norm();
        let v_norm = self.edge_v.norm();
        let area = self.area();
        u_norm > 0.0
            && v_norm > 0.0
            && area > 0.0
            && u_norm.is_finite()
            && v_norm.is_finite()
            && area.is_finite()
            && self.edge_u.dot(self.edge_v).abs() <= RECTANGULARITY_TOLERANCE * u_norm * v_norm
    }

    fn luminance(&self) -> f64 {
        let rgb = self.emission.0.rgb();
        linear_srgb_luminance(rgb) * self.emission.1
    }
}

/// Declared spherical layout of a canonical environment artifact.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnvironmentLayout {
    /// Rows run north-to-south (`+Y` to `-Y`); columns run around `+Y`
    /// starting at `+X` and increasing toward `+Z`.
    LatitudeLongitudeYUp,
}

/// Color interpretation of environment samples.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnvironmentColorInterpretation {
    /// Scene-linear sRGB radiance. Values may exceed one but must be finite and
    /// nonnegative.
    LinearSrgbRadiance,
}

/// Origin of the canonical environment pixels.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnvironmentSourceKind {
    /// Constructed directly from row-major linear-sRGB pixels.
    NativeLinearSrgb,
    /// Decoded by `fs-img` from its supported single-part scanline EXR subset.
    FrankenExrSubset,
}

/// Immutable canonical environment-map artifact.
#[derive(Debug, Clone)]
pub struct EnvironmentMap {
    width: u32,
    height: u32,
    pixels: Vec<[f32; 3]>,
    row_cdf: Vec<f64>,
    column_cdf: Vec<f64>,
    total_weight: f64,
    rotation_y_radians: f64,
    source_kind: EnvironmentSourceKind,
    semantic_hash: ContentHash,
    source_hash: ContentHash,
    provenance_hash: ContentHash,
}

impl EnvironmentMap {
    /// Admit a row-major, scene-linear sRGB radiance map.
    pub fn try_from_linear_srgb(
        width: u32,
        height: u32,
        pixels: Vec<[f32; 3]>,
        rotation_y_radians: f64,
    ) -> Result<Self, LightingError> {
        let source_hash = hash_environment_source(width, height, &pixels);
        Self::try_from_pixels(
            width,
            height,
            pixels,
            rotation_y_radians,
            EnvironmentSourceKind::NativeLinearSrgb,
            source_hash,
        )
    }

    /// Decode the in-house single-part scanline EXR subset and admit its R/G/B
    /// planes as scene-linear radiance. The raw source bytes remain bound in
    /// the provenance hash.
    pub fn try_from_exr(bytes: &[u8], rotation_y_radians: f64) -> Result<Self, LightingError> {
        let decoded =
            fs_img::read_exr(bytes).map_err(|_| LightingError::UnsupportedEnvironmentExr)?;
        let (Some(red), Some(green), Some(blue)) = (
            unique_exr_channel(&decoded.channels, "R"),
            unique_exr_channel(&decoded.channels, "G"),
            unique_exr_channel(&decoded.channels, "B"),
        ) else {
            return Err(LightingError::UnsupportedEnvironmentExr);
        };
        let pixel_count = checked_pixel_count(decoded.width, decoded.height)?;
        if red.len() != pixel_count || green.len() != pixel_count || blue.len() != pixel_count {
            return Err(LightingError::UnsupportedEnvironmentExr);
        }
        let mut pixels = Vec::new();
        pixels
            .try_reserve_exact(pixel_count)
            .map_err(|_| LightingError::EnvironmentTooLarge)?;
        for index in 0..pixel_count {
            pixels.push([red[index], green[index], blue[index]]);
        }
        Self::try_from_pixels(
            decoded.width,
            decoded.height,
            pixels,
            rotation_y_radians,
            EnvironmentSourceKind::FrankenExrSubset,
            hash_domain(ENVIRONMENT_EXR_SOURCE_DOMAIN, bytes),
        )
    }

    fn try_from_pixels(
        width: u32,
        height: u32,
        mut pixels: Vec<[f32; 3]>,
        rotation_y_radians: f64,
        source_kind: EnvironmentSourceKind,
        source_hash: ContentHash,
    ) -> Result<Self, LightingError> {
        let pixel_count = checked_pixel_count(width, height)?;
        if pixels.len() != pixel_count {
            return Err(LightingError::InvalidEnvironmentDimensions);
        }
        if !rotation_y_radians.is_finite() {
            return Err(LightingError::InvalidEnvironmentRotation);
        }
        for (pixel_index, pixel) in pixels.iter_mut().enumerate() {
            if pixel
                .iter()
                .any(|channel| !channel.is_finite() || *channel < 0.0)
            {
                return Err(LightingError::InvalidEnvironmentPixel { pixel_index });
            }
            for channel in pixel {
                if *channel == 0.0 {
                    *channel = 0.0;
                }
            }
        }
        let mut rotation_y_radians = rotation_y_radians.rem_euclid(TAU);
        if rotation_y_radians == 0.0 {
            rotation_y_radians = 0.0;
        }

        let mut row_cdf = Vec::new();
        row_cdf
            .try_reserve_exact(height as usize)
            .map_err(|_| LightingError::EnvironmentTooLarge)?;
        let mut column_cdf = Vec::new();
        column_cdf
            .try_reserve_exact(pixel_count)
            .map_err(|_| LightingError::EnvironmentTooLarge)?;
        let texel_solid_angle_scale = TAU / f64::from(width);
        let mut total_weight = 0.0;
        for row in 0..height {
            let theta_min = PI * f64::from(row) / f64::from(height);
            let theta_max = PI * f64::from(row + 1) / f64::from(height);
            let texel_solid_angle =
                texel_solid_angle_scale * (det::cos(theta_min) - det::cos(theta_max));
            let row_start = row as usize * width as usize;
            let mut row_luminance = 0.0;
            for pixel in &pixels[row_start..row_start + width as usize] {
                row_luminance += linear_srgb_luminance(pixel.map(f64::from));
                column_cdf.push(row_luminance);
            }
            let weight = row_luminance * texel_solid_angle;
            total_weight += weight;
            row_cdf.push(total_weight);
        }
        if !total_weight.is_finite() {
            return Err(LightingError::InvalidEnvironmentPixel { pixel_index: 0 });
        }

        let semantic_hash = hash_environment_semantics(width, height, &pixels, rotation_y_radians);
        let mut provenance = DomainHasher::new(ENVIRONMENT_PROVENANCE_DOMAIN);
        provenance.update(&ENVIRONMENT_IMPORTER_VERSION.to_le_bytes());
        provenance.update(&[match source_kind {
            EnvironmentSourceKind::NativeLinearSrgb => 0,
            EnvironmentSourceKind::FrankenExrSubset => 1,
        }]);
        provenance.update(source_hash.as_bytes());
        provenance.update(semantic_hash.as_bytes());
        let provenance_hash = provenance.finalize();
        Ok(Self {
            width,
            height,
            pixels,
            row_cdf,
            column_cdf,
            total_weight,
            rotation_y_radians,
            source_kind,
            semantic_hash,
            source_hash,
            provenance_hash,
        })
    }

    /// Pixel width.
    #[must_use]
    pub const fn width(&self) -> u32 {
        self.width
    }

    /// Pixel height.
    #[must_use]
    pub const fn height(&self) -> u32 {
        self.height
    }

    /// Canonical Y-axis rotation in `[0, 2π)`.
    #[must_use]
    pub const fn rotation_y_radians(&self) -> f64 {
        self.rotation_y_radians
    }

    /// Declared spherical layout.
    #[must_use]
    pub const fn layout(&self) -> EnvironmentLayout {
        EnvironmentLayout::LatitudeLongitudeYUp
    }

    /// Declared color interpretation.
    #[must_use]
    pub const fn color_interpretation(&self) -> EnvironmentColorInterpretation {
        EnvironmentColorInterpretation::LinearSrgbRadiance
    }

    /// Source/import kind.
    #[must_use]
    pub const fn source_kind(&self) -> EnvironmentSourceKind {
        self.source_kind
    }

    /// Hash of pixels, layout, interpretation, dimensions, and rotation.
    #[must_use]
    pub const fn semantic_hash(&self) -> ContentHash {
        self.semantic_hash
    }

    /// Hash of the native pixel source or raw supported-EXR bytes.
    #[must_use]
    pub const fn source_hash(&self) -> ContentHash {
        self.source_hash
    }

    /// Hash binding semantic content to importer version and source lineage.
    #[must_use]
    pub const fn provenance_hash(&self) -> ContentHash {
        self.provenance_hash
    }

    /// Whether every texel is black.
    #[must_use]
    pub fn is_black(&self) -> bool {
        self.total_weight == 0.0
    }

    /// Importance-sample a direction. The returned PDF is conditional on
    /// selecting this environment, not the scene's outer light mixture.
    #[must_use]
    pub fn sample(&self, u1: f64, u2: f64) -> Option<EnvironmentLightSample> {
        if self.total_weight <= 0.0 || !unit_sample(u1) || !unit_sample(u2) {
            return None;
        }
        let (row, within_row) = select_cdf(&self.row_cdf, self.total_weight, u1)?;
        let row_start = row * self.width as usize;
        let row_cdf = &self.column_cdf[row_start..row_start + self.width as usize];
        let row_total = *row_cdf.last()?;
        let (column, within_column) = select_cdf(row_cdf, row_total, u2)?;

        let theta_min = PI * row as f64 / f64::from(self.height);
        let theta_max = PI * (row + 1) as f64 / f64::from(self.height);
        let cos_theta =
            det::cos(theta_min) + within_row * (det::cos(theta_max) - det::cos(theta_min));
        let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
        let phi =
            TAU * (column as f64 + within_column) / f64::from(self.width) + self.rotation_y_radians;
        let direction = Vec3::new(
            sin_theta * det::cos(phi),
            cos_theta,
            sin_theta * det::sin(phi),
        );
        let pixel_index = row_start + column;
        Some(EnvironmentLightSample {
            direction,
            emission: environment_emission(self.pixels[pixel_index]),
            pdf_solid_angle: linear_srgb_luminance(self.pixels[pixel_index].map(f64::from))
                / self.total_weight,
        })
    }

    /// Evaluate piecewise-constant radiance and its conditional directional
    /// PDF for a world-space direction.
    #[must_use]
    pub fn evaluate(&self, direction: Vec3) -> Option<EnvironmentEvaluation> {
        let norm = direction.norm();
        if !norm.is_finite() || norm <= 0.0 {
            return None;
        }
        let direction = direction.scale(1.0 / norm);
        let theta = det::acos(direction.y.clamp(-1.0, 1.0));
        let phi = (det::atan2(direction.z, direction.x) - self.rotation_y_radians).rem_euclid(TAU);
        let row = ((theta / PI) * f64::from(self.height)).floor() as usize;
        let column = ((phi / TAU) * f64::from(self.width)).floor() as usize;
        let row = row.min(self.height as usize - 1);
        let column = column.min(self.width as usize - 1);
        let pixel = self.pixels[row * self.width as usize + column];
        let conditional_pdf = if self.total_weight > 0.0 {
            linear_srgb_luminance(pixel.map(f64::from)) / self.total_weight
        } else {
            0.0
        };
        Some(EnvironmentEvaluation {
            emission: environment_emission(pixel),
            pdf_solid_angle: conditional_pdf,
        })
    }
}

/// Sample of one rectangular emitter, including the complete scene-mixture
/// solid-angle PDF.
#[derive(Debug, Clone, Copy)]
pub struct RectLightSample {
    /// Original input-light index.
    pub light_index: usize,
    /// Matching emissive primitive.
    pub primitive_index: usize,
    /// Uniformly sampled point on the rectangle.
    pub point: Point3,
    /// Rectangle unit normal.
    pub normal: Vec3,
    /// Emitted spectral radiance.
    pub emission: (LiftedSpectrum, f64),
    /// Selection probability times uniform-area density converted to solid
    /// angle at the shading point.
    pub pdf_solid_angle: f64,
}

/// Sample of the environment, including the complete scene-mixture
/// solid-angle PDF.
#[derive(Debug, Clone, Copy)]
pub struct EnvironmentLightSample {
    /// Sampled world direction.
    pub direction: Vec3,
    /// Piecewise-constant texel radiance.
    pub emission: (LiftedSpectrum, f64),
    /// Directional PDF, including environment selection probability when
    /// returned by [`AdmittedLighting::sample`].
    pub pdf_solid_angle: f64,
}

/// Environment radiance and directional PDF for a queried direction.
#[derive(Debug, Clone, Copy)]
pub struct EnvironmentEvaluation {
    /// Piecewise-constant texel radiance.
    pub emission: (LiftedSpectrum, f64),
    /// Directional PDF, including environment selection probability when
    /// returned by [`AdmittedLighting::environment_evaluation`].
    pub pdf_solid_angle: f64,
}

/// One direct-light sample.
#[derive(Debug, Clone, Copy)]
pub enum LightSample {
    /// Finite rectangular source.
    Rectangle(RectLightSample),
    /// Infinite environment source.
    Environment(EnvironmentLightSample),
}

/// Kind of a candidate in selection diagnostics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LightKind {
    /// Rectangular source; carries the original input index.
    Rectangle(usize),
    /// Environment source.
    Environment,
}

/// One candidate's point-dependent selection diagnostics.
#[derive(Debug, Clone, Copy)]
pub struct LightSelectionEntry {
    /// Candidate kind.
    pub kind: LightKind,
    /// Stable content identity used for ordering.
    pub identity: ContentHash,
    /// Unnormalized incident-luminance weight.
    pub weight: f64,
    /// Normalized selection probability.
    pub probability: f64,
}

/// Complete point-dependent selection diagnostics.
#[derive(Debug, Clone)]
pub struct LightSelectionDiagnostics {
    /// Candidates in deterministic content-identity order.
    pub entries: Vec<LightSelectionEntry>,
    /// Sum of all unnormalized weights.
    pub total_weight: f64,
}

#[derive(Debug, Clone, Copy)]
struct Candidate {
    kind: LightKind,
    identity: ContentHash,
}

/// Validated, construction-order-independent view of a scene's emitters.
#[derive(Debug)]
pub struct AdmittedLighting<'a> {
    rectangles: &'a [RectLight],
    environment: Option<&'a EnvironmentMap>,
    candidates: Vec<Candidate>,
}

impl<'a> AdmittedLighting<'a> {
    /// Validate rectangle geometry/emission, reject duplicate identities and
    /// primitive bindings, and build deterministic candidate order.
    pub fn try_new(
        rectangles: &'a [RectLight],
        environment: Option<&'a EnvironmentMap>,
    ) -> Result<Self, LightingError> {
        Self::try_new_controlled(rectangles, environment, || Ok(()))
    }

    /// Cancellable form of [`Self::try_new`] used by production render
    /// admission. It polls before every emitter and every ordered-candidate
    /// publication step, so an arbitrarily large caller-owned light list does
    /// not become one uninterruptible preflight operation.
    pub fn try_new_cancellable(
        cx: &Cx<'_>,
        rectangles: &'a [RectLight],
        environment: Option<&'a EnvironmentMap>,
    ) -> Result<Self, LightingError> {
        Self::try_new_controlled(rectangles, environment, || {
            cx.checkpoint().map_err(|_| LightingError::Cancelled)
        })
    }

    fn try_new_controlled(
        rectangles: &'a [RectLight],
        environment: Option<&'a EnvironmentMap>,
        mut checkpoint: impl FnMut() -> Result<(), LightingError>,
    ) -> Result<Self, LightingError> {
        checkpoint()?;
        let candidate_capacity = rectangles
            .len()
            .checked_add(usize::from(environment.is_some()))
            .ok_or(LightingError::EnvironmentTooLarge)?;
        // Ordered insertion avoids handing an unbounded caller-owned light
        // list to an uncancellable sort. The tuple preserves the former stable
        // sort semantics for the cryptographic-collision case: rectangles in
        // input order, followed by the environment.
        let mut ordered_candidates = BTreeMap::new();
        let mut rectangle_identity_owners = BTreeMap::new();
        let mut duplicate_rectangle = None;
        let mut primitives = BTreeSet::new();
        for (light_index, rectangle) in rectangles.iter().enumerate() {
            checkpoint()?;
            let luminance = rectangle.luminance();
            if !rectangle.is_valid() || !luminance.is_finite() || luminance <= 0.0 {
                return Err(LightingError::InvalidRectangle { light_index });
            }
            if !primitives.insert(rectangle.prim) {
                return Err(LightingError::DuplicatePrimitive {
                    primitive_index: rectangle.prim,
                });
            }
            let identity = rectangle.identity();
            if let Some(first) = rectangle_identity_owners.get(&identity).copied() {
                if duplicate_rectangle
                    .as_ref()
                    .is_none_or(|(current, _, _)| identity < *current)
                {
                    duplicate_rectangle = Some((identity, first, light_index));
                }
            } else {
                rectangle_identity_owners.insert(identity, light_index);
            }
            ordered_candidates.insert(
                (identity, 0_u8, light_index),
                Candidate {
                    kind: LightKind::Rectangle(light_index),
                    identity,
                },
            );
        }
        checkpoint()?;
        // A black map still supplies explicit black miss radiance, but it is
        // not a sampling candidate. In particular, adding a black map must not
        // perturb the frozen one-rectangle random-number mapping.
        if let Some(environment) = environment.filter(|map| !map.is_black()) {
            let identity = environment.semantic_hash();
            ordered_candidates.insert(
                (identity, 1_u8, 0),
                Candidate {
                    kind: LightKind::Environment,
                    // Sampling order is a property of rendered content. Two
                    // byte-distinct containers with identical canonical pixels
                    // therefore retain one stream mapping; source lineage remains
                    // available separately through `provenance_hash`.
                    identity,
                },
            );
        }
        if let Some((_, first, second)) = duplicate_rectangle {
            return Err(LightingError::DuplicateRectangle { first, second });
        }
        if rectangles.is_empty() && environment.is_none_or(EnvironmentMap::is_black) {
            return Err(LightingError::NoFiniteEmitter);
        }
        let mut candidates = Vec::new();
        candidates
            .try_reserve_exact(candidate_capacity)
            .map_err(|_| LightingError::EnvironmentTooLarge)?;
        for (_, candidate) in ordered_candidates {
            checkpoint()?;
            candidates.push(candidate);
        }
        checkpoint()?;
        Ok(Self {
            rectangles,
            environment,
            candidates,
        })
    }

    /// Sample one candidate with a deterministic residual remap. The exact
    /// single-rectangle/no-environment case passes `u1,u2` straight through.
    #[must_use]
    pub fn sample(&self, origin: Point3, u1: f64, u2: f64) -> Option<LightSample> {
        if !unit_sample(u1) || !unit_sample(u2) {
            return None;
        }
        if self.candidates.len() == 1
            && let LightKind::Rectangle(light_index) = self.candidates[0].kind
        {
            return self.sample_rectangle(
                light_index,
                origin,
                u1,
                u2,
                1.0,
                LEGACY_RECTANGLE_COSINE_CUTOFF,
            );
        }
        let total_weight = self.total_weight(origin);
        if total_weight <= 0.0 || !total_weight.is_finite() {
            return None;
        }
        let target = u1 * total_weight;
        let mut cumulative = 0.0;
        let mut last_positive = None;
        for candidate in &self.candidates {
            let weight = self.candidate_weight(*candidate, origin);
            if weight <= 0.0 {
                continue;
            }
            last_positive = Some((*candidate, cumulative, weight));
            let upper = cumulative + weight;
            if target < upper {
                let residual = clamp_unit_open((target - cumulative) / weight);
                return self.sample_candidate(
                    *candidate,
                    origin,
                    residual,
                    u2,
                    weight / total_weight,
                );
            }
            cumulative = upper;
        }
        let (candidate, lower, weight) = last_positive?;
        let residual = clamp_unit_open((target - lower) / weight);
        self.sample_candidate(candidate, origin, residual, u2, weight / total_weight)
    }

    fn sample_candidate(
        &self,
        candidate: Candidate,
        origin: Point3,
        u1: f64,
        u2: f64,
        selection_probability: f64,
    ) -> Option<LightSample> {
        match candidate.kind {
            LightKind::Rectangle(light_index) => self.sample_rectangle(
                light_index,
                origin,
                u1,
                u2,
                selection_probability,
                MULTI_LIGHT_RECTANGLE_COSINE_CUTOFF,
            ),
            LightKind::Environment => {
                let mut sample = self.environment?.sample(u1, u2)?;
                sample.pdf_solid_angle *= selection_probability;
                Some(LightSample::Environment(sample))
            }
        }
    }

    fn sample_rectangle(
        &self,
        light_index: usize,
        origin: Point3,
        u1: f64,
        u2: f64,
        selection_probability: f64,
        cosine_cutoff: f64,
    ) -> Option<LightSample> {
        let light = self.rectangles.get(light_index)?;
        let point = light
            .corner
            .offset(light.edge_u.scale(u1))
            .offset(light.edge_v.scale(u2));
        let direction = point.delta_from(origin);
        let distance_squared = direction.dot(direction);
        if !(distance_squared > 0.0 && distance_squared.is_finite()) {
            return None;
        }
        let direction_unit = direction.scale(1.0 / distance_squared.sqrt());
        let cosine = light.normal().dot(direction_unit).abs();
        if cosine <= cosine_cutoff {
            return None;
        }
        let conditional_pdf = distance_squared / (cosine * light.area());
        Some(LightSample::Rectangle(RectLightSample {
            light_index,
            primitive_index: light.prim,
            point,
            normal: light.normal(),
            emission: light.emission,
            pdf_solid_angle: selection_probability * conditional_pdf,
        }))
    }

    /// Full mixture PDF for a BSDF-sampled point on an admitted rectangle.
    #[must_use]
    pub fn rect_mixture_pdf(&self, light_index: usize, origin: Point3, hit_point: Point3) -> f64 {
        let Some(light) = self.rectangles.get(light_index) else {
            return 0.0;
        };
        let total = self.total_weight(origin);
        let weight = self.rectangle_weight(light, origin);
        if total <= 0.0 || weight <= 0.0 {
            return 0.0;
        }
        let direction = hit_point.delta_from(origin);
        let distance_squared = direction.dot(direction);
        if !(distance_squared > 0.0 && distance_squared.is_finite()) {
            return 0.0;
        }
        let direction = direction.scale(1.0 / distance_squared.sqrt());
        let cosine = light.normal().dot(direction).abs();
        if cosine <= 1.0e-12 {
            0.0
        } else {
            let conditional_pdf = distance_squared / (cosine * light.area());
            (weight / total) * conditional_pdf
        }
    }

    /// Original rectangle index for an emissive primitive, if NEE can sample
    /// that primitive.
    #[must_use]
    pub fn rect_index_for_primitive(&self, primitive_index: usize) -> Option<usize> {
        self.rectangles
            .iter()
            .position(|light| light.prim == primitive_index)
    }

    /// Environment radiance and full mixture PDF for a BSDF-sampled miss.
    #[must_use]
    pub fn environment_evaluation(
        &self,
        origin: Point3,
        direction: Vec3,
    ) -> Option<EnvironmentEvaluation> {
        let environment = self.environment?;
        let mut evaluation = environment.evaluate(direction)?;
        let total = self.total_weight(origin);
        if total > 0.0 {
            evaluation.pdf_solid_angle *= environment.total_weight / total;
        } else {
            evaluation.pdf_solid_angle = 0.0;
        }
        Some(evaluation)
    }

    /// Point-dependent ordered candidate probabilities for diagnostics and
    /// statistical tests.
    #[must_use]
    pub fn diagnostics(&self, origin: Point3) -> LightSelectionDiagnostics {
        let total_weight = self.total_weight(origin);
        let entries = self
            .candidates
            .iter()
            .map(|candidate| {
                let weight = self.candidate_weight(*candidate, origin);
                LightSelectionEntry {
                    kind: candidate.kind,
                    identity: candidate.identity,
                    weight,
                    probability: if total_weight > 0.0 {
                        weight / total_weight
                    } else {
                        0.0
                    },
                }
            })
            .collect();
        LightSelectionDiagnostics {
            entries,
            total_weight,
        }
    }

    /// Whether this admission uses the frozen one-rectangle/no-environment
    /// estimator rather than the separately versioned lighting extension.
    #[must_use]
    pub(crate) fn is_legacy_compatibility_path(&self) -> bool {
        self.rectangles.len() == 1 && self.environment.is_none_or(EnvironmentMap::is_black)
    }

    fn total_weight(&self, origin: Point3) -> f64 {
        self.candidates
            .iter()
            .map(|candidate| self.candidate_weight(*candidate, origin))
            .sum()
    }

    fn candidate_weight(&self, candidate: Candidate, origin: Point3) -> f64 {
        match candidate.kind {
            LightKind::Rectangle(light_index) => self
                .rectangles
                .get(light_index)
                .map_or(0.0, |light| self.rectangle_weight(light, origin)),
            LightKind::Environment => self.environment.map_or(0.0, |map| map.total_weight),
        }
    }

    fn rectangle_weight(&self, light: &RectLight, origin: Point3) -> f64 {
        light.luminance() * rectangle_solid_angle(light, origin)
    }
}

fn unique_exr_channel<'a>(channels: &'a [fs_img::Channel], name: &str) -> Option<&'a [f32]> {
    let mut matching = channels
        .iter()
        .filter(|channel| channel.name == name)
        .map(|channel| channel.data.as_slice());
    let channel = matching.next()?;
    matching.next().is_none().then_some(channel)
}

fn checked_pixel_count(width: u32, height: u32) -> Result<usize, LightingError> {
    if width == 0 || height == 0 {
        return Err(LightingError::InvalidEnvironmentDimensions);
    }
    (width as usize)
        .checked_mul(height as usize)
        .ok_or(LightingError::InvalidEnvironmentDimensions)
}

fn canonical_f64_bits(value: f64) -> u64 {
    if value == 0.0 {
        0.0_f64.to_bits()
    } else {
        value.to_bits()
    }
}

fn canonical_f32_bits(value: f32) -> u32 {
    if value == 0.0 {
        0.0_f32.to_bits()
    } else {
        value.to_bits()
    }
}

fn hash_environment_source(width: u32, height: u32, pixels: &[[f32; 3]]) -> ContentHash {
    let mut hasher = DomainHasher::new(ENVIRONMENT_NATIVE_SOURCE_DOMAIN);
    hasher.update(&width.to_le_bytes());
    hasher.update(&height.to_le_bytes());
    for pixel in pixels {
        for channel in pixel {
            hasher.update(&canonical_f32_bits(*channel).to_le_bytes());
        }
    }
    hasher.finalize()
}

fn hash_environment_semantics(
    width: u32,
    height: u32,
    pixels: &[[f32; 3]],
    rotation_y_radians: f64,
) -> ContentHash {
    let mut hasher = DomainHasher::new(ENVIRONMENT_SEMANTIC_DOMAIN);
    hasher.update(&ENVIRONMENT_SEMANTICS_VERSION.to_le_bytes());
    hasher.update(&[0]); // LatitudeLongitudeYUp.
    hasher.update(&[0]); // LinearSrgbRadiance.
    hasher.update(&width.to_le_bytes());
    hasher.update(&height.to_le_bytes());
    hasher.update(&canonical_f64_bits(rotation_y_radians).to_le_bytes());
    for pixel in pixels {
        for channel in pixel {
            hasher.update(&canonical_f32_bits(*channel).to_le_bytes());
        }
    }
    hasher.finalize()
}

fn linear_srgb_luminance(rgb: [f64; 3]) -> f64 {
    0.212_672_9 * rgb[0] + 0.715_152_2 * rgb[1] + 0.072_175_0 * rgb[2]
}

fn environment_emission(pixel: [f32; 3]) -> (LiftedSpectrum, f64) {
    let scale = pixel.into_iter().fold(0.0_f32, f32::max);
    if scale == 0.0 {
        // The spectral shape is unobservable at zero scale. Avoid paying the
        // one-time LUT construction cost for an explicitly black miss.
        return (LiftedSpectrum { c: [0.0; 3] }, 0.0);
    }
    let normalized = pixel.map(|channel| f64::from(channel / scale));
    (interpolate_spectral_lut(normalized), f64::from(scale))
}

fn spectral_lut() -> &'static [LiftedSpectrum] {
    static LUT: OnceLock<Vec<LiftedSpectrum>> = OnceLock::new();
    LUT.get_or_init(|| {
        let mut table = Vec::with_capacity(SPECTRAL_LUT_EDGE.pow(3));
        let denominator = (SPECTRAL_LUT_EDGE - 1) as f64;
        for red in 0..SPECTRAL_LUT_EDGE {
            for green in 0..SPECTRAL_LUT_EDGE {
                for blue in 0..SPECTRAL_LUT_EDGE {
                    table.push(lift_rgb([
                        red as f64 / denominator,
                        green as f64 / denominator,
                        blue as f64 / denominator,
                    ]));
                }
            }
        }
        table
    })
}

fn spectral_lut_index(red: usize, green: usize, blue: usize) -> usize {
    (red * SPECTRAL_LUT_EDGE + green) * SPECTRAL_LUT_EDGE + blue
}

fn interpolate_spectral_lut(rgb: [f64; 3]) -> LiftedSpectrum {
    let scale = (SPECTRAL_LUT_EDGE - 1) as f64;
    let coordinate = rgb.map(|channel| channel.clamp(0.0, 1.0) * scale);
    let low = coordinate.map(|value| value.floor() as usize);
    let high = low.map(|value| (value + 1).min(SPECTRAL_LUT_EDGE - 1));
    let fraction = [
        coordinate[0] - low[0] as f64,
        coordinate[1] - low[1] as f64,
        coordinate[2] - low[2] as f64,
    ];
    let table = spectral_lut();
    let mut coefficients = [0.0; 3];
    for red_side in 0..2 {
        for green_side in 0..2 {
            for blue_side in 0..2 {
                let red = if red_side == 0 { low[0] } else { high[0] };
                let green = if green_side == 0 { low[1] } else { high[1] };
                let blue = if blue_side == 0 { low[2] } else { high[2] };
                let weight = if red_side == 0 {
                    1.0 - fraction[0]
                } else {
                    fraction[0]
                } * if green_side == 0 {
                    1.0 - fraction[1]
                } else {
                    fraction[1]
                } * if blue_side == 0 {
                    1.0 - fraction[2]
                } else {
                    fraction[2]
                };
                let spectrum = table[spectral_lut_index(red, green, blue)];
                for (coefficient, source) in coefficients.iter_mut().zip(spectrum.c) {
                    *coefficient += weight * source;
                }
            }
        }
    }
    LiftedSpectrum { c: coefficients }
}

fn unit_sample(value: f64) -> bool {
    value.is_finite() && (0.0..1.0).contains(&value)
}

fn clamp_unit_open(value: f64) -> f64 {
    value.clamp(0.0, 1.0_f64.next_down())
}

fn select_cdf(cdf: &[f64], total: f64, sample: f64) -> Option<(usize, f64)> {
    if cdf.is_empty() || total <= 0.0 || !total.is_finite() || !unit_sample(sample) {
        return None;
    }
    let target = sample * total;
    let mut low = 0usize;
    let mut high = cdf.len();
    while low < high {
        let middle = low + (high - low) / 2;
        if target < cdf[middle] {
            high = middle;
        } else {
            low = middle + 1;
        }
    }
    // A rounded target may equal the final CDF even though `sample < 1`.
    // Fall back to the final positive bucket without changing its residual
    // convention.
    let mut index = low.min(cdf.len() - 1);
    while index > 0 && cdf[index] == cdf[index - 1] {
        index -= 1;
    }
    let lower = if index == 0 { 0.0 } else { cdf[index - 1] };
    let weight = cdf[index] - lower;
    if weight <= 0.0 {
        return None;
    }
    Some((index, clamp_unit_open((target - lower) / weight)))
}

fn rectangle_solid_angle(light: &RectLight, origin: Point3) -> f64 {
    let normal = light.normal();
    let plane_distance = light.corner.delta_from(origin).dot(normal).abs();
    let length_scale = light
        .edge_u
        .norm()
        .max(light.edge_v.norm())
        .max(f64::MIN_POSITIVE);
    if plane_distance <= PLANE_DISTANCE_TOLERANCE * length_scale {
        return 0.0;
    }
    let p0 = light.corner.delta_from(origin);
    let p1 = light.corner.offset(light.edge_u).delta_from(origin);
    let p2 = light
        .corner
        .offset(light.edge_u)
        .offset(light.edge_v)
        .delta_from(origin);
    let p3 = light.corner.offset(light.edge_v).delta_from(origin);
    triangle_solid_angle(p0, p1, p2) + triangle_solid_angle(p0, p2, p3)
}

fn triangle_solid_angle(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    let (a_norm, b_norm, c_norm) = (a.norm(), b.norm(), c.norm());
    if a_norm <= 0.0 || b_norm <= 0.0 || c_norm <= 0.0 {
        return 0.0;
    }
    let numerator = a.dot(cross(b, c)).abs();
    let denominator =
        a_norm * b_norm * c_norm + a.dot(b) * c_norm + b.dot(c) * a_norm + c.dot(a) * b_norm;
    2.0 * det::atan2(numerator, denominator)
}

fn cross(left: Vec3, right: Vec3) -> Vec3 {
    Vec3::new(
        left.y * right.z - left.z * right.y,
        left.z * right.x - left.x * right.z,
        left.x * right.y - left.y * right.x,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_img::{Channel, PixelType, write_exr};

    fn rectangle(x: f64, primitive: usize, scale: f64) -> RectLight {
        RectLight {
            corner: Point3::new(x, 1.0, -0.5),
            edge_u: Vec3::new(1.0, 0.0, 0.0),
            edge_v: Vec3::new(0.0, 0.0, 1.0),
            prim: primitive,
            emission: (lift_rgb([1.0, 0.8, 0.6]), scale),
        }
    }

    fn assert_close(actual: f64, expected: f64, tolerance: f64, context: &str) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{context}: actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
        );
    }

    #[test]
    fn g0_rectangle_admission_rejects_geometry_power_and_duplicates() {
        let valid = rectangle(-0.5, 7, 2.0);
        AdmittedLighting::try_new(&[valid], None).expect("valid rectangle");

        let mut skew = valid;
        skew.edge_v = Vec3::new(0.25, 0.0, 1.0);
        assert_eq!(
            AdmittedLighting::try_new(&[skew], None).unwrap_err(),
            LightingError::InvalidRectangle { light_index: 0 }
        );
        let mut zero_power = valid;
        zero_power.emission.1 = 0.0;
        assert_eq!(
            AdmittedLighting::try_new(&[zero_power], None).unwrap_err(),
            LightingError::InvalidRectangle { light_index: 0 }
        );
        let mut overflowing_spectrum = valid;
        overflowing_spectrum.emission.0.c = [f64::MAX; 3];
        assert_eq!(
            AdmittedLighting::try_new(&[overflowing_spectrum], None).unwrap_err(),
            LightingError::InvalidRectangle { light_index: 0 },
            "a finite coefficient vector whose evaluated luminance overflows must refuse at admission"
        );
        let mut duplicate_primitive = rectangle(1.0, 7, 2.0);
        duplicate_primitive.emission.1 = 3.0;
        assert_eq!(
            AdmittedLighting::try_new(&[valid, duplicate_primitive], None).unwrap_err(),
            LightingError::DuplicatePrimitive { primitive_index: 7 }
        );
        let mut duplicate_identity = valid;
        duplicate_identity.prim = 8;
        assert!(matches!(
            AdmittedLighting::try_new(&[valid, duplicate_identity], None),
            Err(LightingError::DuplicateRectangle { .. })
        ));
    }

    #[test]
    fn g4_light_admission_polls_during_validation_and_ordered_publication() {
        let lights = [
            rectangle(-2.0, 10, 1.0),
            rectangle(0.0, 11, 2.0),
            rectangle(2.0, 12, 3.0),
        ];
        // Poll four interrupts before the third input is inspected. Poll six
        // interrupts when canonical candidates begin moving into published
        // storage. These pins keep both phases bounded without timing races.
        for cancellation_poll in [4, 6] {
            let mut polls = 0;
            let result = AdmittedLighting::try_new_controlled(&lights, None, || {
                polls += 1;
                if polls == cancellation_poll {
                    Err(LightingError::Cancelled)
                } else {
                    Ok(())
                }
            });
            assert_eq!(
                result.unwrap_err(),
                LightingError::Cancelled,
                "cancellation_poll={cancellation_poll} polls={polls}: admission continued past its bounded checkpoint"
            );
            assert_eq!(polls, cancellation_poll);
        }
    }

    #[test]
    fn g0_black_environment_is_valid_but_cannot_be_the_only_emitter() {
        let black = EnvironmentMap::try_from_linear_srgb(2, 2, vec![[0.0; 3]; 4], 0.0)
            .expect("finite black map");
        assert!(black.is_black());
        assert_eq!(
            AdmittedLighting::try_new(&[], Some(&black)).unwrap_err(),
            LightingError::NoFiniteEmitter
        );
        let rectangle = rectangle(0.0, 1, 1.0);
        let rectangles = [rectangle];
        let without = AdmittedLighting::try_new(&rectangles, None).unwrap();
        let with_black = AdmittedLighting::try_new(&rectangles, Some(&black))
            .expect("black environment may accompany a finite rectangle");
        let origin = Point3::new(0.0, 0.0, 0.0);
        let LightSample::Rectangle(without) = without.sample(origin, 0.371, 0.829).unwrap() else {
            panic!("rectangle-only rig selected an environment")
        };
        let LightSample::Rectangle(with_black) = with_black.sample(origin, 0.371, 0.829).unwrap()
        else {
            panic!("black environment became a sampling candidate")
        };
        assert_eq!(
            without.point.x.to_bits(),
            with_black.point.x.to_bits(),
            "a black environment must not perturb the one-light stream"
        );
        assert_eq!(without.point.y.to_bits(), with_black.point.y.to_bits());
        assert_eq!(without.point.z.to_bits(), with_black.point.z.to_bits());
        assert_eq!(
            without.pdf_solid_angle.to_bits(),
            with_black.pdf_solid_angle.to_bits()
        );
    }

    #[test]
    fn g0_single_rectangle_preserves_surface_sample_arithmetic() {
        let light = rectangle(-0.5, 3, 4.0);
        let lights = [light];
        let admitted = AdmittedLighting::try_new(&lights, None).unwrap();
        let origin = Point3::new(0.0, 0.0, 0.0);
        let LightSample::Rectangle(sample) = admitted.sample(origin, 0.25, 0.75).unwrap() else {
            panic!("single rectangle selected environment");
        };
        let expected = light
            .corner
            .offset(light.edge_u.scale(0.25))
            .offset(light.edge_v.scale(0.75));
        assert_eq!(sample.point.x.to_bits(), expected.x.to_bits());
        assert_eq!(sample.point.y.to_bits(), expected.y.to_bits());
        assert_eq!(sample.point.z.to_bits(), expected.z.to_bits());
        assert_eq!(sample.primitive_index, 3);

        let displacement = expected.delta_from(origin);
        let distance_squared = displacement.dot(displacement);
        let direction = displacement.scale(1.0 / distance_squared.sqrt());
        let expected_pdf = distance_squared / (light.normal().dot(direction).abs() * light.area());
        assert_eq!(sample.pdf_solid_angle.to_bits(), expected_pdf.to_bits());
    }

    #[test]
    fn g5_single_rectangle_preserves_the_legacy_grazing_nee_cutoff() {
        let light = rectangle(-0.5, 3, 4.0);
        let lights = [light];
        let admitted = AdmittedLighting::try_new(&lights, None).unwrap();
        let origin = Point3::new(1.0e10, 0.0, 0.0);
        assert!(
            admitted.sample(origin, 0.5, 0.5).is_none(),
            "the frozen tracer-v1 NEE path discarded rectangle samples with cosine <= 1e-9"
        );
        let hit = light
            .corner
            .offset(light.edge_u.scale(0.5))
            .offset(light.edge_v.scale(0.5));
        assert!(
            admitted.rect_mixture_pdf(0, origin, hit) > 0.0,
            "the independent tracer-v1 BSDF-hit PDF used its original 1e-12 cutoff"
        );
    }

    #[test]
    fn g5_rectangle_construction_order_does_not_change_stream_mapping() {
        let left = rectangle(-2.0, 10, 2.0);
        let right = rectangle(1.0, 11, 5.0);
        let forward = [left, right];
        let reversed = [right, left];
        let forward = AdmittedLighting::try_new(&forward, None).unwrap();
        let reversed = AdmittedLighting::try_new(&reversed, None).unwrap();
        let origin = Point3::new(0.0, 0.0, 0.0);
        for index in 0..4096 {
            let u1 = (f64::from(index) + 0.5) / 4096.0;
            let u2 = (f64::from((index * 1543) % 4096) + 0.5) / 4096.0;
            let LightSample::Rectangle(a) = forward.sample(origin, u1, u2).unwrap() else {
                panic!("forward selected environment");
            };
            let LightSample::Rectangle(b) = reversed.sample(origin, u1, u2).unwrap() else {
                panic!("reversed selected environment");
            };
            assert_eq!(a.primitive_index, b.primitive_index, "sample={index}");
            assert_eq!(a.point.x.to_bits(), b.point.x.to_bits(), "sample={index}");
            assert_eq!(a.point.y.to_bits(), b.point.y.to_bits(), "sample={index}");
            assert_eq!(a.point.z.to_bits(), b.point.z.to_bits(), "sample={index}");
            assert_eq!(
                a.pdf_solid_angle.to_bits(),
                b.pdf_solid_angle.to_bits(),
                "sample={index}"
            );
        }
    }

    #[test]
    fn g0_constant_environment_has_uniform_sphere_pdf_and_normalized_mass() {
        let map =
            EnvironmentMap::try_from_linear_srgb(8, 4, vec![[2.0, 2.0, 2.0]; 32], 0.0).unwrap();
        let expected_pdf = 1.0 / (4.0 * PI);
        for index in 0..1024 {
            let u1 = (f64::from(index) + 0.5) / 1024.0;
            let u2 = (f64::from((index * 613) % 1024) + 0.5) / 1024.0;
            let sample = map.sample(u1, u2).unwrap();
            assert_close(sample.direction.norm(), 1.0, 2.0e-15, "sample norm");
            assert_close(
                sample.pdf_solid_angle,
                expected_pdf,
                2.0e-15,
                "constant environment PDF",
            );
        }
        let mut integrated = 0.0;
        for row in 0..map.height as usize {
            let theta_min = PI * row as f64 / f64::from(map.height);
            let theta_max = PI * (row + 1) as f64 / f64::from(map.height);
            let solid_angle =
                TAU / f64::from(map.width) * (det::cos(theta_min) - det::cos(theta_max));
            for column in 0..map.width as usize {
                let pixel = map.pixels[row * map.width as usize + column];
                integrated +=
                    linear_srgb_luminance(pixel.map(f64::from)) / map.total_weight * solid_angle;
            }
        }
        assert_close(integrated, 1.0, 3.0e-15, "environment PDF integral");
    }

    #[test]
    fn g3_bright_texel_importance_sampling_never_selects_black_texels() {
        let mut pixels = vec![[0.0; 3]; 8];
        pixels[6] = [1000.0, 50.0, 1.0];
        let map = EnvironmentMap::try_from_linear_srgb(4, 2, pixels, 0.0).unwrap();
        for index in 0..512 {
            let sample = map
                .sample(
                    (f64::from(index) + 0.5) / 512.0,
                    (f64::from((index * 197) % 512) + 0.5) / 512.0,
                )
                .unwrap();
            assert_eq!(sample.emission.1.to_bits(), 1000.0_f64.to_bits());
            assert!(sample.pdf_solid_angle.is_finite());
            assert!(sample.pdf_solid_angle > 0.0);
        }
    }

    #[test]
    fn g3_rotation_seam_and_poles_are_finite_and_directional() {
        let pixels = vec![
            [8.0, 0.1, 0.1],
            [0.1, 8.0, 0.1],
            [0.1, 0.1, 8.0],
            [1.0, 1.0, 1.0],
            [8.0, 0.1, 0.1],
            [0.1, 8.0, 0.1],
            [0.1, 0.1, 8.0],
            [1.0, 1.0, 1.0],
        ];
        let unrotated = EnvironmentMap::try_from_linear_srgb(4, 2, pixels.clone(), 0.0).unwrap();
        let rotated = EnvironmentMap::try_from_linear_srgb(4, 2, pixels, PI / 2.0).unwrap();
        let direction = Vec3::new(1.0, 0.0, 0.0);
        let a = unrotated.evaluate(direction).unwrap();
        let b = rotated.evaluate(direction).unwrap();
        assert_ne!(
            a.emission.0.c.map(f64::to_bits),
            b.emission.0.c.map(f64::to_bits)
        );
        for pole in [Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, -1.0, 0.0)] {
            let evaluation = rotated.evaluate(pole).unwrap();
            assert!(evaluation.pdf_solid_angle.is_finite());
            assert!(evaluation.emission.1.is_finite());
        }
        let left = rotated.evaluate(Vec3::new(1.0, 0.0, -0.0)).unwrap();
        let right = rotated.evaluate(Vec3::new(1.0, 0.0, 0.0)).unwrap();
        assert_eq!(
            left.emission.0.c.map(f64::to_bits),
            right.emission.0.c.map(f64::to_bits),
            "signed zero must not split the longitude seam"
        );
    }

    #[test]
    fn g0_environment_validation_and_exr_provenance_are_explicit() {
        assert_eq!(
            EnvironmentMap::try_from_linear_srgb(0, 1, Vec::new(), 0.0).unwrap_err(),
            LightingError::InvalidEnvironmentDimensions
        );
        assert_eq!(
            EnvironmentMap::try_from_linear_srgb(1, 1, vec![[-0.1, 0.0, 0.0]], 0.0).unwrap_err(),
            LightingError::InvalidEnvironmentPixel { pixel_index: 0 }
        );
        assert_eq!(
            EnvironmentMap::try_from_linear_srgb(1, 1, vec![[0.0; 3]], f64::NAN).unwrap_err(),
            LightingError::InvalidEnvironmentRotation
        );

        let pixels = [[0.25_f32, 1.5, 3.0], [2.0, 0.5, 0.125]];
        let channel = |name: &str, component: usize| Channel {
            name: name.to_string(),
            ty: PixelType::Float,
            data: pixels.iter().map(|pixel| pixel[component]).collect(),
        };
        let bytes = write_exr(2, 1, &[channel("R", 0), channel("G", 1), channel("B", 2)]).unwrap();
        let duplicate_channels = [
            channel("R", 0),
            channel("R", 0),
            channel("G", 1),
            channel("B", 2),
        ];
        assert!(
            unique_exr_channel(&duplicate_channels, "R").is_none(),
            "EXR admission must refuse duplicate required channel names"
        );
        let native = EnvironmentMap::try_from_linear_srgb(2, 1, pixels.to_vec(), 0.25).unwrap();
        let imported = EnvironmentMap::try_from_exr(&bytes, 0.25).unwrap();
        assert_eq!(native.semantic_hash(), imported.semantic_hash());
        assert_ne!(native.source_hash(), imported.source_hash());
        assert_ne!(native.provenance_hash(), imported.provenance_hash());
        assert_eq!(
            imported.source_kind(),
            EnvironmentSourceKind::FrankenExrSubset
        );

        let light = rectangle(-0.5, 12, 3.0);
        let lights = [light];
        let native_lighting = AdmittedLighting::try_new(&lights, Some(&native)).unwrap();
        let imported_lighting = AdmittedLighting::try_new(&lights, Some(&imported)).unwrap();
        let origin = Point3::new(0.0, 0.0, 0.0);
        let native_diagnostics = native_lighting.diagnostics(origin);
        let imported_diagnostics = imported_lighting.diagnostics(origin);
        assert_eq!(
            native_diagnostics
                .entries
                .iter()
                .map(|entry| entry.identity)
                .collect::<Vec<_>>(),
            imported_diagnostics
                .entries
                .iter()
                .map(|entry| entry.identity)
                .collect::<Vec<_>>(),
            "container provenance must not alter sampling order for identical canonical radiance"
        );
        assert_eq!(
            EnvironmentMap::try_from_exr(b"not an exr", 0.0).unwrap_err(),
            LightingError::UnsupportedEnvironmentExr
        );
    }

    #[test]
    fn g0_environment_spectral_lut_is_bounded_finite_and_color_resolving() {
        let red = environment_emission([4.0, 0.1, 0.1]);
        let blue = environment_emission([0.1, 0.1, 4.0]);
        for wavelength in [380.0, 450.0, 550.0, 650.0, 780.0] {
            for emission in [red, blue] {
                let value = emission.0.eval(wavelength) * emission.1;
                assert!(value.is_finite());
                assert!((0.0..=emission.1).contains(&value));
            }
        }
        assert_ne!(red.0.c.map(f64::to_bits), blue.0.c.map(f64::to_bits));
    }

    #[test]
    fn g3_mixture_frequencies_match_reported_probabilities() {
        let lights = [rectangle(-2.0, 1, 1.0), rectangle(1.0, 2, 4.0)];
        let environment =
            EnvironmentMap::try_from_linear_srgb(2, 1, vec![[0.3, 0.4, 0.5]; 2], 0.0).unwrap();
        let admitted = AdmittedLighting::try_new(&lights, Some(&environment)).unwrap();
        let origin = Point3::new(0.0, 0.0, 0.0);
        let diagnostics = admitted.diagnostics(origin);
        assert_close(
            diagnostics
                .entries
                .iter()
                .map(|entry| entry.probability)
                .sum(),
            1.0,
            2.0e-15,
            "selection probability sum",
        );
        let sample_count = 20_000usize;
        let mut counts = [0usize; 3];
        for index in 0..sample_count {
            let u1 = (index as f64 + 0.5) / sample_count as f64;
            let u2 = ((index * 7919) % sample_count) as f64 / sample_count as f64;
            let sample = admitted.sample(origin, u1, u2).unwrap();
            let kind = match sample {
                LightSample::Rectangle(sample) if sample.primitive_index == 1 => 0,
                LightSample::Rectangle(sample) if sample.primitive_index == 2 => 1,
                LightSample::Environment(_) => 2,
                LightSample::Rectangle(sample) => panic!(
                    "unexpected primitive {} at deterministic sample {index}",
                    sample.primitive_index
                ),
            };
            counts[kind] += 1;
        }
        for (kind, observed) in [
            (LightKind::Rectangle(0), counts[0]),
            (LightKind::Rectangle(1), counts[1]),
            (LightKind::Environment, counts[2]),
        ] {
            let expected = diagnostics
                .entries
                .iter()
                .find(|entry| entry.kind == kind)
                .unwrap()
                .probability;
            let observed = observed as f64 / sample_count as f64;
            assert_close(
                observed,
                expected,
                1.0 / sample_count as f64,
                "selection frequency",
            );
        }
    }

    #[test]
    fn g0_sampled_and_bsdf_hit_paths_use_the_same_complete_mixture_pdf() {
        let lights = [rectangle(-2.0, 17, 1.5), rectangle(1.0, 23, 4.0)];
        let environment = EnvironmentMap::try_from_linear_srgb(
            4,
            2,
            vec![
                [0.2, 0.4, 0.8],
                [1.0, 0.3, 0.1],
                [0.1, 0.2, 0.6],
                [0.7, 0.5, 0.2],
                [0.2, 0.4, 0.8],
                [1.0, 0.3, 0.1],
                [0.1, 0.2, 0.6],
                [0.7, 0.5, 0.2],
            ],
            0.37,
        )
        .unwrap();
        let admitted = AdmittedLighting::try_new(&lights, Some(&environment)).unwrap();
        let origin = Point3::new(0.0, 0.0, 0.0);
        let mut saw = [false; 3];
        for index in 0..4096 {
            let u1 = (index as f64 + 0.5) / 4096.0;
            let u2 = ((index * 1543) % 4096) as f64 / 4096.0;
            match admitted.sample(origin, u1, u2).unwrap() {
                LightSample::Rectangle(sample) => {
                    let hit_pdf =
                        admitted.rect_mixture_pdf(sample.light_index, origin, sample.point);
                    assert_eq!(
                        sample.pdf_solid_angle.to_bits(),
                        hit_pdf.to_bits(),
                        "rectangle {} NEE/BSDF-hit PDF mismatch at sample {index}: nee={:.17e}, hit={hit_pdf:.17e}",
                        sample.light_index,
                        sample.pdf_solid_angle
                    );
                    saw[sample.light_index] = true;
                }
                LightSample::Environment(sample) => {
                    let hit = admitted
                        .environment_evaluation(origin, sample.direction)
                        .unwrap();
                    assert_eq!(
                        sample.pdf_solid_angle.to_bits(),
                        hit.pdf_solid_angle.to_bits(),
                        "environment NEE/BSDF-miss PDF mismatch at sample {index}: nee={:.17e}, miss={:.17e}",
                        sample.pdf_solid_angle,
                        hit.pdf_solid_angle
                    );
                    saw[2] = true;
                }
            }
        }
        assert!(
            saw.into_iter().all(core::convert::identity),
            "deterministic PDF partition fixture failed to exercise every candidate: {saw:?}"
        );
    }
}
