//! Dense tiled SDF grids (plan §7.2): f32 STORAGE / f64 EVALUATION on
//! fs-substrate's Morton/tile-major fields, with C¹ triquadratic B-spline
//! reconstruction — continuous gradients matter because shape optimization
//! DIFFERENTIATES through samples (plan §7.6).
//!
//! Error honesty: the chart separately carries a CONSTRUCTED nominal-field
//! reconstruction bound (source Lipschitz times the farthest actual spline
//! support node, f32 quantization, and interpolation roundoff) and the weakest
//! source authority observed relative to the abstract region's signed
//! distance. A source `NoClaim` is never laundered into an enclosure by
//! sampling. MEASURED eikonal statistics ("how much is this NOT a distance
//! field") remain separate, clearly-labeled evidence.

use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::Cx;
use fs_geom::{
    Aabb, Chart, ChartSample, Differentiability, Point3, SamplingDomain, SamplingDomainError,
    TraceStepClaim, Vec3,
};

/// A dense signed-field grid over a box. Finite sampled Lipschitz values admit
/// a nominal reconstruction; rigorous abstract-distance authority additionally
/// requires the source's global `ExactDistance` theorem.
#[derive(Debug)]
pub struct TiledSdf {
    field: fs_substrate::field::TiledField<f32>,
    box_: Aabb,
    /// Strictly increasing representable sample coordinates per axis.
    axis_nodes: [Vec<f64>; 3],
    /// Largest actual grid gap per axis.
    h: [f64; 3],
    /// Downward bound on the smallest actual grid gap per axis.
    h_min: [f64; 3],
    /// Samples per axis.
    n: [u32; 3],
    /// Constructed reconstruction bound relative to the sampled source field.
    nominal_field_bound: f64,
    /// Weakest source authority observed across the grid, after spline
    /// reconstruction demotes `Exact` to `Enclosure`.
    abstract_distance_kind: NumericalKind,
    /// Total reconstruction-plus-source bound relative to abstract region
    /// signed distance. `None` means honest `NoClaim`.
    abstract_distance_bound: Option<f64>,
    /// Maximum source Lipschitz value observed during sampling.
    source_lipschitz: f64,
    /// Preflighted finite reconstructed-field Lipschitz candidate. It is
    /// exposed as a claim only with `has_global_lipschitz_theorem`.
    chart_lipschitz: f64,
    /// Whether the source supplied a global ExactDistance theorem rather than
    /// only sampled local Lipschitz values.
    has_global_lipschitz_theorem: bool,
    /// Monotone maximum absolute sampled value, retained so transactional
    /// incremental refreshes can conservatively update quantization slack.
    max_sample_abs: f64,
    /// Measured max |,∇φ| − 1, over the construction probe set (an
    /// ESTIMATE, clearly labeled — see [`TiledSdf::eikonal_stats`]).
    eikonal_dev: f64,
}

/// Construction failure (Decalogue P10).
#[derive(Debug, Clone, PartialEq)]
pub enum SdfBuildError {
    /// The source support or explicit clip is not an admissible finite
    /// three-dimensional sampling domain.
    SamplingDomain(SamplingDomainError),
    /// A sampling spacing was not finite and strictly positive.
    InvalidSpacing {
        /// Which constructor argument was invalid.
        field: &'static str,
        /// The offending value.
        value: f64,
    },
    /// The source chart certifies no finite Lipschitz bound.
    NoLipschitzBound,
    /// The requested step would need more samples per axis than the cap.
    ResolutionTooFine {
        /// Samples/axis the request needs.
        need: u64,
        /// The cap.
        cap: u64,
    },
    /// Checked multiplication of per-axis sample counts overflowed the
    /// addressable allocation domain.
    SampleCountOverflow {
        /// Samples per axis that could not be multiplied safely.
        dims: [u64; 3],
    },
    /// A nominal dense lattice could not form consecutive coordinates with a
    /// strictly positive representable gap bound.
    DenseLatticeUnrepresentable {
        /// Axis containing the collapsed coordinate.
        axis: usize,
        /// Index of the coordinate that failed to advance.
        index: u32,
        /// Raw bits of the preceding coordinate.
        previous_bits: u64,
        /// Raw bits of the attempted coordinate.
        coordinate_bits: u64,
    },
    /// An adaptive build's worst-case octree would exceed its deterministic
    /// work cap.
    AdaptiveWorkLimit {
        /// Maximum nodes implied by the requested depth.
        need: u128,
        /// Deterministic node cap.
        cap: u64,
    },
    /// An adaptive cell cannot be split into strictly smaller representable
    /// children on one axis.
    AdaptiveSubdivisionUnrepresentable {
        /// Axis whose midpoint rounded to an endpoint or became non-finite.
        axis: usize,
        /// Raw bits of the cell's minimum endpoint.
        min_bits: u64,
        /// Raw bits of the cell's maximum endpoint.
        max_bits: u64,
        /// Raw bits of the attempted midpoint.
        midpoint_bits: u64,
    },
    /// A narrow-band build would scan too many dense lattice points before
    /// sparsification.
    BandScanLimit {
        /// Per-axis lattice dimensions.
        dims: [u64; 3],
        /// Total lattice points required.
        need: u128,
        /// Deterministic scan cap.
        cap: u64,
    },
    /// A narrow-band lattice would leave the signed VDB coordinate domain.
    CoordinateRange {
        /// Axis whose coordinate range is not representable.
        axis: usize,
        /// Samples needed on that axis.
        need: u64,
        /// Largest supported sample count.
        cap: u64,
    },
    /// Construction observed cancellation at a bounded polling point.
    Cancelled,
    /// A source sample was not finite or could not be represented in the
    /// representation's storage field without becoming infinite.
    InvalidSample {
        /// Point at which the invalid sample was observed.
        point: Point3,
        /// Raw f64 bits of the source's nominal field value.
        value_bits: u64,
    },
    /// The constructed nominal-field reconstruction bound overflowed.
    InvalidReconstructionBound {
        /// Non-finite constructed bound.
        value: f64,
    },
    /// The Lipschitz claim derived for the reconstructed chart overflowed or
    /// became invalid.
    InvalidDerivedLipschitz {
        /// Invalid derived claim.
        value: f64,
    },
}

impl core::fmt::Display for SdfBuildError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            SdfBuildError::SamplingDomain(error) => write!(f, "{error}"),
            SdfBuildError::InvalidSpacing { field, value } => write!(
                f,
                "SDF build refused: `{field}` must be finite and strictly positive, got {value}"
            ),
            SdfBuildError::NoLipschitzBound => write!(
                f,
                "dense SDF build refused: the source chart certifies no Lipschitz bound, so \
                 no rigorous sampling error exists; use a certified source"
            ),
            SdfBuildError::ResolutionTooFine { need, cap } => write!(
                f,
                "dense SDF build refused: {need} samples/axis exceed the {cap} cap; coarsen \
                 the step, shrink the box, or use the sparse VDB/adaptive charts"
            ),
            SdfBuildError::SampleCountOverflow { dims } => write!(
                f,
                "SDF build refused: sample dimensions {dims:?} overflow the addressable grid"
            ),
            SdfBuildError::DenseLatticeUnrepresentable {
                axis,
                index,
                previous_bits,
                coordinate_bits,
            } => write!(
                f,
                "dense SDF build refused: axis {axis} node {index} has no usable representable \
                 gap beyond {} (attempted {}); coarsen the node count or rescale coordinates",
                f64::from_bits(*previous_bits),
                f64::from_bits(*coordinate_bits)
            ),
            SdfBuildError::AdaptiveWorkLimit { need, cap } => write!(
                f,
                "adaptive SDF build refused: worst-case octree work {need} nodes exceeds the {cap} node cap; reduce max_depth"
            ),
            SdfBuildError::AdaptiveSubdivisionUnrepresentable {
                axis,
                min_bits,
                max_bits,
                midpoint_bits,
            } => write!(
                f,
                "adaptive SDF build refused: axis {axis} cell [{}, {}] has no finite strictly \
                 interior f64 midpoint (attempted {}); reduce max_depth or rescale coordinates",
                f64::from_bits(*min_bits),
                f64::from_bits(*max_bits),
                f64::from_bits(*midpoint_bits)
            ),
            SdfBuildError::BandScanLimit { dims, need, cap } => write!(
                f,
                "narrow-band build refused: lattice {dims:?} requires {need} samples, exceeding the {cap} sample scan cap; coarsen h or shrink the clip"
            ),
            SdfBuildError::CoordinateRange { axis, need, cap } => write!(
                f,
                "narrow-band build refused: axis {axis} needs {need} samples but signed VDB coordinates allow at most {cap}; coarsen h or shrink the clip"
            ),
            SdfBuildError::Cancelled => write!(f, "SDF build cancelled at a bounded polling point"),
            SdfBuildError::InvalidSample { point, value_bits } => write!(
                f,
                "SDF build refused: source sample at {point:?} is non-finite or outside the target storage range (f64 bits {value_bits:#018x})"
            ),
            SdfBuildError::InvalidReconstructionBound { value } => write!(
                f,
                "SDF build refused: constructed nominal-field reconstruction bound is non-finite ({value})"
            ),
            SdfBuildError::InvalidDerivedLipschitz { value } => write!(
                f,
                "SDF build refused: derived reconstructed-field Lipschitz claim is non-finite or negative ({value})"
            ),
        }
    }
}

impl core::error::Error for SdfBuildError {}

impl From<SamplingDomainError> for SdfBuildError {
    fn from(error: SamplingDomainError) -> Self {
        Self::SamplingDomain(error)
    }
}

/// Per-axis sample cap (beyond this, dense is the wrong tool — the
/// refusal says which chart to use instead).
pub const DENSE_MAX_SAMPLES_PER_AXIS: u64 = 512;

pub(crate) fn finite_positive(value: f64, field: &'static str) -> Result<f64, SdfBuildError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(SdfBuildError::InvalidSpacing { field, value })
    }
}

fn checked_axis_samples(
    span: f64,
    step: f64,
    minimum: u64,
    cap: u64,
) -> Result<u32, SdfBuildError> {
    let cells = (span / step).ceil();
    if !cells.is_finite() || cells < 0.0 || cells > (cap.saturating_sub(1)) as f64 {
        let need = if cells.is_finite() {
            (cells as u64).saturating_add(1)
        } else {
            u64::MAX
        };
        return Err(SdfBuildError::ResolutionTooFine { need, cap });
    }
    let nodes = (cells as u64)
        .checked_add(1)
        .ok_or(SdfBuildError::SampleCountOverflow { dims: [cap; 3] })?
        .max(minimum);
    u32::try_from(nodes).map_err(|_| SdfBuildError::SampleCountOverflow { dims: [nodes; 3] })
}

fn checked_sample_product(dims: [u32; 3]) -> Result<usize, SdfBuildError> {
    let dims_u64 = dims.map(u64::from);
    let total = dims_u64
        .iter()
        .try_fold(1u128, |product, &dim| product.checked_mul(u128::from(dim)))
        .ok_or(SdfBuildError::SampleCountOverflow { dims: dims_u64 })?;
    usize::try_from(total).map_err(|_| SdfBuildError::SampleCountOverflow { dims: dims_u64 })
}

fn stable_lerp(min: f64, max: f64, t: f64) -> f64 {
    if t <= 0.0 {
        return min;
    }
    if t >= 1.0 {
        return max;
    }
    if min.is_sign_negative() == max.is_sign_negative() {
        min + (max - min) * t
    } else {
        min * (1.0 - t) + max * t
    }
}

fn build_axis_nodes(
    min: f64,
    max: f64,
    count: u32,
    axis: usize,
) -> Result<(Vec<f64>, f64, f64), SdfBuildError> {
    let mut nodes = Vec::with_capacity(count as usize);
    let denominator = f64::from(count - 1);
    for index in 0..count {
        let coordinate = if index == count - 1 {
            max
        } else {
            stable_lerp(min, max, f64::from(index) / denominator).clamp(min, max)
        };
        if let Some(&previous) = nodes.last()
            && (!coordinate.is_finite() || coordinate <= previous)
        {
            return Err(SdfBuildError::DenseLatticeUnrepresentable {
                axis,
                index,
                previous_bits: previous.to_bits(),
                coordinate_bits: coordinate.to_bits(),
            });
        }
        nodes.push(coordinate);
    }
    let mut min_gap = f64::INFINITY;
    let mut max_gap = 0.0f64;
    for (lower_index, pair) in nodes.windows(2).enumerate() {
        let gap = pair[1] - pair[0];
        let gap_lower = gap.next_down();
        if !gap.is_finite() || gap <= 0.0 || gap_lower <= 0.0 {
            return Err(SdfBuildError::DenseLatticeUnrepresentable {
                axis,
                index: u32::try_from(lower_index + 1).unwrap_or(u32::MAX),
                previous_bits: pair[0].to_bits(),
                coordinate_bits: pair[1].to_bits(),
            });
        }
        min_gap = min_gap.min(gap_lower);
        max_gap = max_gap.max(gap);
    }
    Ok((nodes, min_gap, max_gap))
}

fn locate_axis(nodes: &[f64], coordinate: f64) -> (usize, usize, f64) {
    let last = nodes.len() - 1;
    if coordinate <= nodes[0] {
        return (0, 1, 0.0);
    }
    if coordinate >= nodes[last] {
        return (last - 1, last, 1.0);
    }
    let upper = nodes.partition_point(|node| *node <= coordinate);
    let lower = upper - 1;
    let width = nodes[upper] - nodes[lower];
    let t = ((coordinate - nodes[lower]) / width).clamp(0.0, 1.0);
    (lower, upper, t)
}

/// Validate how strongly one source sample relates its nominal field value to
/// abstract region signed distance. Certificate fields are public, so a
/// malformed finite claim must fail closed rather than be laundered by the
/// sampler.
pub(crate) fn sample_abstract_distance_authority(
    sample: &ChartSample,
) -> (NumericalKind, Option<f64>) {
    let certificate = sample.error;
    if certificate.kind == NumericalKind::NoClaim {
        return (NumericalKind::NoClaim, None);
    }
    let finite_ordered = certificate.lo.is_finite()
        && certificate.hi.is_finite()
        && certificate.lo <= certificate.hi;
    let consistent = match certificate.kind {
        NumericalKind::Exact => {
            certificate.lo.to_bits() == sample.signed_distance.to_bits()
                && certificate.hi.to_bits() == sample.signed_distance.to_bits()
        }
        NumericalKind::Enclosure | NumericalKind::Estimate => {
            finite_ordered
                && certificate.lo <= sample.signed_distance
                && sample.signed_distance <= certificate.hi
        }
        NumericalKind::NoClaim => false,
    };
    if !finite_ordered || !consistent {
        return (NumericalKind::NoClaim, None);
    }
    let radius = nonnegative_difference_upper(sample.signed_distance, certificate.lo).max(
        nonnegative_difference_upper(certificate.hi, sample.signed_distance),
    );
    if radius.is_finite() {
        (certificate.kind, Some(radius))
    } else {
        (NumericalKind::NoClaim, None)
    }
}

fn checked_stored_sample(sample: &ChartSample, point: Point3) -> Result<f32, SdfBuildError> {
    let stored = sample.signed_distance as f32;
    if sample.signed_distance.is_finite() && stored.is_finite() {
        Ok(stored)
    } else {
        Err(SdfBuildError::InvalidSample {
            point,
            value_bits: sample.signed_distance.to_bits(),
        })
    }
}

fn nonnegative_mul_upper(lhs: f64, rhs: f64) -> f64 {
    if lhs == 0.0 || rhs == 0.0 {
        return 0.0;
    }
    let product = lhs * rhs;
    if product.is_finite() {
        product.next_up()
    } else {
        product
    }
}

fn nonnegative_add_upper(lhs: f64, rhs: f64) -> f64 {
    if lhs == 0.0 {
        return rhs;
    }
    if rhs == 0.0 {
        return lhs;
    }
    let sum = lhs + rhs;
    if sum.is_finite() { sum.next_up() } else { sum }
}

fn nonnegative_div_upper(numerator: f64, denominator: f64) -> f64 {
    if numerator == 0.0 {
        return 0.0;
    }
    let quotient = numerator / denominator;
    if quotient.is_finite() {
        quotient.next_up()
    } else {
        quotient
    }
}

fn nonnegative_difference_upper(lhs: f64, rhs: f64) -> f64 {
    let difference = (lhs - rhs).abs();
    if difference == 0.0 {
        0.0
    } else if difference.is_finite() {
        difference.next_up()
    } else {
        difference
    }
}

fn outward_norm3(x: f64, y: f64, z: f64) -> Option<f64> {
    let scale = x.max(y).max(z);
    if !scale.is_finite() || scale < 0.0 {
        return None;
    }
    if scale == 0.0 {
        return Some(0.0);
    }
    let square = |value: f64| {
        let ratio = nonnegative_div_upper(value, scale);
        nonnegative_mul_upper(ratio, ratio)
    };
    let sum = nonnegative_add_upper(nonnegative_add_upper(square(x), square(y)), square(z));
    let norm = nonnegative_mul_upper(scale, sum.sqrt().next_up());
    norm.is_finite().then_some(norm)
}

fn spline_support_radius(max_gaps: [f64; 3]) -> Result<f64, SdfBuildError> {
    let radius = outward_norm3(
        nonnegative_mul_upper(2.0, max_gaps[0].next_up()),
        nonnegative_mul_upper(2.0, max_gaps[1].next_up()),
        nonnegative_mul_upper(2.0, max_gaps[2].next_up()),
    )
    .unwrap_or(f64::INFINITY);
    if radius.is_finite() {
        Ok(radius)
    } else {
        Err(SdfBuildError::InvalidReconstructionBound { value: radius })
    }
}

fn interpolation_roundoff_bound(max_sample_abs: f64) -> f64 {
    // The tensor-product value path forms 27 weighted contributions. This
    // deliberately loose gamma-style allowance also covers basis arithmetic,
    // accumulation, and subnormal absolute error.
    let relative = (2048.0 * f64::EPSILON).next_up();
    let scaled = nonnegative_mul_upper(max_sample_abs, relative);
    let subnormal_floor = 2048.0 * f64::from_bits(1);
    scaled.max(subnormal_floor).next_up()
}

fn finite_add_upper(lhs: f64, rhs: f64) -> Option<f64> {
    let sum = lhs + rhs;
    sum.is_finite().then(|| sum.next_up())
}

fn finite_subtract_lower(lhs: f64, rhs: f64) -> Option<f64> {
    let difference = lhs - rhs;
    difference.is_finite().then(|| difference.next_down())
}

fn centered_interval_hull(center: f64, radius: f64, nominal: f64) -> Option<(f64, f64)> {
    if !center.is_finite() || !radius.is_finite() || radius < 0.0 || !nominal.is_finite() {
        return None;
    }
    let lo = finite_subtract_lower(center, radius)?.min(nominal);
    let hi = finite_add_upper(center, radius)?.max(nominal);
    (lo.is_finite() && hi.is_finite() && lo <= hi).then_some((lo, hi))
}

fn finite_norm_with_upper(x: f64, y: f64, z: f64) -> Option<(f64, f64)> {
    let (x, y, z) = (x.abs(), y.abs(), z.abs());
    let scale = x.max(y).max(z);
    if !scale.is_finite() {
        return None;
    }
    if scale == 0.0 {
        return Some((0.0, 0.0));
    }
    let nominal = scale
        * ((x / scale) * (x / scale) + (y / scale) * (y / scale) + (z / scale) * (z / scale))
            .sqrt();
    let component_upper = |value: f64| if value == 0.0 { 0.0 } else { value.next_up() };
    let upper = outward_norm3(component_upper(x), component_upper(y), component_upper(z))?;
    (nominal.is_finite() && upper.is_finite()).then_some((nominal, upper))
}

fn normalized_direction(direction: Vec3) -> Option<Vec3> {
    let scale = direction
        .x
        .abs()
        .max(direction.y.abs())
        .max(direction.z.abs());
    if !scale.is_finite() || scale == 0.0 {
        return None;
    }
    let scaled = Vec3::new(
        direction.x / scale,
        direction.y / scale,
        direction.z / scale,
    );
    let norm = (scaled.x * scaled.x + scaled.y * scaled.y + scaled.z * scaled.z).sqrt();
    if !norm.is_finite() || norm == 0.0 {
        return None;
    }
    let unit = Vec3::new(scaled.x / norm, scaled.y / norm, scaled.z / norm);
    (unit.x.is_finite() && unit.y.is_finite() && unit.z.is_finite()).then_some(unit)
}

fn difference_quotient_interval(bound: f64, origin: f64, direction: f64) -> Option<(f64, f64)> {
    let difference = bound - origin;
    if !difference.is_finite() || direction == 0.0 || !direction.is_finite() {
        return None;
    }
    let difference_lo = difference.next_down();
    let difference_hi = difference.next_up();
    if !difference_lo.is_finite() || !difference_hi.is_finite() {
        return None;
    }
    let quotient_a = difference_lo / direction;
    let quotient_b = difference_hi / direction;
    if !quotient_a.is_finite() || !quotient_b.is_finite() {
        return None;
    }
    let lo = quotient_a.min(quotient_b).next_down();
    let hi = quotient_a.max(quotient_b).next_up();
    (lo.is_finite() && hi.is_finite()).then_some((lo, hi))
}

fn ray_box_interval(origin: Point3, direction: Vec3, box_: Aabb) -> Option<(f64, f64)> {
    let origins = [origin.x, origin.y, origin.z];
    let directions = [direction.x, direction.y, direction.z];
    let minima = [box_.min.x, box_.min.y, box_.min.z];
    let maxima = [box_.max.x, box_.max.y, box_.max.z];
    let mut enter = f64::NEG_INFINITY;
    let mut exit = f64::INFINITY;
    for axis in 0..3 {
        if directions[axis] == 0.0 {
            if origins[axis] < minima[axis] || origins[axis] > maxima[axis] {
                return None;
            }
            continue;
        }
        let lower_t = difference_quotient_interval(minima[axis], origins[axis], directions[axis])?;
        let upper_t = difference_quotient_interval(maxima[axis], origins[axis], directions[axis])?;
        let near = lower_t.0.min(upper_t.0);
        let far = lower_t.1.max(upper_t.1);
        enter = enter.max(near);
        exit = exit.min(far);
        if enter > exit {
            return None;
        }
    }
    let enter = enter.max(0.0);
    (enter.is_finite() && exit.is_finite() && enter <= exit && exit >= 0.0).then_some((enter, exit))
}

fn f32_storage_rounding_bound(max_sample_abs: f64) -> f64 {
    // Round-to-nearest f64→f32 conversion has at most half a minimum f32
    // subnormal of absolute error near zero. A relative-only term misses
    // underflow to zero and subnormal spacing entirely.
    let half_min_subnormal = 0.5 * f64::from(f32::from_bits(1));
    nonnegative_mul_upper(max_sample_abs, f64::from(f32::EPSILON)).max(half_min_subnormal)
}

fn constructed_nominal_bound(
    lipschitz: f64,
    support_radius: f64,
    max_sample_abs: f64,
) -> Result<f64, SdfBuildError> {
    let sampling = nonnegative_mul_upper(lipschitz, support_radius);
    let quantization = f32_storage_rounding_bound(max_sample_abs);
    let stored_magnitude = nonnegative_add_upper(max_sample_abs, quantization);
    let interpolation = interpolation_roundoff_bound(stored_magnitude);
    let bound = nonnegative_add_upper(nonnegative_add_upper(sampling, quantization), interpolation);
    if bound.is_finite() && bound >= 0.0 {
        Ok(bound)
    } else {
        Err(SdfBuildError::InvalidReconstructionBound { value: bound })
    }
}

fn constructed_chart_lipschitz(
    source_lipschitz: f64,
    nominal_field_bound: f64,
    h_min: f64,
) -> Result<f64, SdfBuildError> {
    let reconstruction_slope =
        nonnegative_div_upper(nonnegative_mul_upper(2.0, nominal_field_bound), h_min);
    let value = nonnegative_add_upper(source_lipschitz, reconstruction_slope);
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(SdfBuildError::InvalidDerivedLipschitz { value })
    }
}

fn composed_abstract_bound(
    kind: &mut NumericalKind,
    nominal_field_bound: f64,
    source_error_radius: f64,
) -> Option<f64> {
    // Spline reconstruction is approximate even if every stored source
    // sample was exact.
    *kind = (*kind).max(NumericalKind::Enclosure);
    if *kind == NumericalKind::NoClaim {
        return None;
    }
    let bound = nonnegative_add_upper(nominal_field_bound, source_error_radius);
    if bound.is_finite() {
        Some(bound)
    } else {
        *kind = NumericalKind::NoClaim;
        None
    }
}

/// Measured eikonal statistics (`|∇φ| − 1` over a seeded probe set):
/// evidence for sphere-tracing safety, LABELED as measurement.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EikonalStats {
    /// Mean absolute deviation.
    pub mean_abs_dev: f64,
    /// Maximum absolute deviation observed.
    pub max_abs_dev: f64,
    /// Probe count.
    pub probes: u64,
}

impl TiledSdf {
    /// Sample `source` over its inflated support at step `target_h`.
    /// The nominal field-fit band uses the outward distance to the farthest
    /// actual spline-support node, plus f32 storage and interpolation
    /// roundoff. It becomes a rigorous
    /// abstract-distance enclosure only when `source` supplies the global
    /// `ExactDistance` theorem; sampled local Lipschitz maxima alone yield at
    /// best Estimate authority.
    ///
    /// # Errors
    /// [`SdfBuildError`] (teaching refusals; nothing runs before checks).
    pub fn build(
        source: &dyn Chart,
        target_h: f64,
        cx: &Cx<'_>,
    ) -> Result<TiledSdf, SdfBuildError> {
        let target_h = finite_positive(target_h, "target_h")?;
        let padding = finite_positive(3.0 * target_h, "3 * target_h")?;
        let support = SamplingDomain::admit(source.support(), None)?.bounds();
        let box_ = SamplingDomain::admit(support.inflate(padding), None)?.bounds();
        Self::build_in_domain(source, box_, target_h, cx)
    }

    /// Sample the geometric intersection `source ∩ clip` at step
    /// `target_h`. The explicit clip is validated before any source
    /// evaluation or grid allocation.
    ///
    /// # Errors
    /// [`SdfBuildError`] when the clip, spacing, or requested grid is not
    /// admissible.
    pub fn build_clipped(
        source: &dyn Chart,
        target_h: f64,
        clip: Aabb,
        cx: &Cx<'_>,
    ) -> Result<TiledSdf, SdfBuildError> {
        let clipped = fs_geom::ClippedChart::new(source, clip)?;
        Self::build(&clipped, target_h, cx)
    }

    #[allow(clippy::too_many_lines)] // One ordered admission, sampling, authority, and publication transaction.
    fn build_in_domain(
        source: &dyn Chart,
        box_: Aabb,
        target_h: f64,
        cx: &Cx<'_>,
    ) -> Result<TiledSdf, SdfBuildError> {
        let domain = SamplingDomain::admit(box_, None)?;
        let box_ = domain.bounds();
        let spans = domain.spans();
        let edges = [spans.x, spans.y, spans.z];
        let endpoints = [
            (box_.min.x, box_.max.x),
            (box_.min.y, box_.max.y),
            (box_.min.z, box_.max.z),
        ];
        let mut n = [0u32; 3];
        let mut h = [0.0f64; 3];
        let mut h_min = [0.0f64; 3];
        let mut axis_nodes = [Vec::new(), Vec::new(), Vec::new()];
        for a in 0..3 {
            n[a] = checked_axis_samples(edges[a], target_h, 4, DENSE_MAX_SAMPLES_PER_AXIS)?;
            let (nodes, min_gap, max_gap) =
                build_axis_nodes(endpoints[a].0, endpoints[a].1, n[a], a)?;
            h_min[a] = min_gap;
            h[a] = max_gap;
            axis_nodes[a] = nodes;
        }
        let _total_samples = checked_sample_product(n)?;
        let has_global_lipschitz_theorem =
            source.trace_step_claim() == TraceStepClaim::ExactDistance;
        cx.checkpoint().map_err(|_| SdfBuildError::Cancelled)?;
        let probe = source.eval(
            Point3::new(
                f64::midpoint(box_.min.x, box_.max.x),
                f64::midpoint(box_.min.y, box_.max.y),
                f64::midpoint(box_.min.z, box_.max.z),
            ),
            cx,
        );
        let mut lipschitz = match probe.lipschitz {
            Some(l) if l.is_finite() && l >= 0.0 => l,
            _ if has_global_lipschitz_theorem => 1.0,
            _ => return Err(SdfBuildError::NoLipschitzBound),
        };
        if has_global_lipschitz_theorem {
            // An exact signed-distance theorem implies global L = 1. A local
            // sample hint may tighten other chart classes, but cannot weaken
            // this theorem or the rigorous reconstruction band.
            lipschitz = lipschitz.max(1.0);
        }
        let (mut abstract_distance_kind, probe_radius) = sample_abstract_distance_authority(&probe);
        let mut source_error_radius = probe_radius.unwrap_or(0.0);
        cx.checkpoint().map_err(|_| SdfBuildError::Cancelled)?;
        let grid = fs_substrate::tile::TileGrid::new(n, fs_substrate::tile::TileEdge::E8)
            .expect("caps keep the grid within Morton bounds");
        let mut field = fs_substrate::field::TiledField::new(grid, 0.0f32);
        let mut max_abs = 0.0f64;
        for k in 0..n[2] {
            for j in 0..n[1] {
                cx.checkpoint().map_err(|_| SdfBuildError::Cancelled)?;
                for i in 0..n[0] {
                    let p = Point3::new(
                        axis_nodes[0][i as usize],
                        axis_nodes[1][j as usize],
                        axis_nodes[2][k as usize],
                    );
                    let sample = source.eval(p, cx);
                    let stored = checked_stored_sample(&sample, p)?;
                    match sample.lipschitz {
                        Some(local) if local.is_finite() && local >= 0.0 => {
                            lipschitz = lipschitz.max(local);
                        }
                        _ if has_global_lipschitz_theorem => {}
                        _ => return Err(SdfBuildError::NoLipschitzBound),
                    }
                    let (sample_kind, sample_radius) = sample_abstract_distance_authority(&sample);
                    abstract_distance_kind = abstract_distance_kind.max(sample_kind);
                    if let Some(radius) = sample_radius {
                        source_error_radius = source_error_radius.max(radius);
                    }
                    max_abs = max_abs.max(sample.signed_distance.abs());
                    field.set([i, j, k], stored);
                }
            }
        }
        cx.checkpoint().map_err(|_| SdfBuildError::Cancelled)?;
        let support_radius = spline_support_radius(h)?;
        let nominal_field_bound = constructed_nominal_bound(lipschitz, support_radius, max_abs)?;
        let chart_lipschitz = constructed_chart_lipschitz(
            lipschitz,
            nominal_field_bound,
            h_min[0].min(h_min[1]).min(h_min[2]),
        )?;
        if !has_global_lipschitz_theorem {
            // ChartSample::lipschitz is local. A finite maximum observed on a
            // grid is useful for a nominal fit but is not a global theorem
            // over unsampled cell interiors, so it cannot mint an enclosure.
            abstract_distance_kind = abstract_distance_kind.max(NumericalKind::Estimate);
        }
        let abstract_distance_bound = composed_abstract_bound(
            &mut abstract_distance_kind,
            nominal_field_bound,
            source_error_radius,
        );
        let mut sdf = TiledSdf {
            field,
            box_,
            axis_nodes,
            h,
            h_min,
            n,
            nominal_field_bound,
            abstract_distance_kind,
            abstract_distance_bound,
            source_lipschitz: lipschitz,
            chart_lipschitz,
            has_global_lipschitz_theorem,
            max_sample_abs: max_abs,
            eikonal_dev: 0.0,
        };
        sdf.eikonal_dev = sdf.measure_eikonal(0x51DF_0001, 2_000, cx)?.max_abs_dev;
        cx.checkpoint().map_err(|_| SdfBuildError::Cancelled)?;
        Ok(sdf)
    }

    /// Reconstruction error relative to the sampled source field.
    ///
    /// This compatibility accessor is not abstract-region signed-distance
    /// authority. Consumers making such a claim must inspect
    /// [`Self::abstract_distance_kind`] and [`Self::abstract_distance_bound`].
    #[must_use]
    pub fn bound(&self) -> f64 {
        self.nominal_field_bound
    }

    /// Reconstruction error relative to the sampled source field.
    #[must_use]
    pub fn nominal_field_bound(&self) -> f64 {
        self.nominal_field_bound
    }

    /// Weakest abstract signed-distance authority observed during sampling.
    #[must_use]
    pub fn abstract_distance_kind(&self) -> NumericalKind {
        self.abstract_distance_kind
    }

    /// Total error relative to abstract region signed distance when every
    /// sampled source value carried a finite enclosure or estimate.
    #[must_use]
    pub fn abstract_distance_bound(&self) -> Option<f64> {
        self.abstract_distance_bound
    }

    /// Monotonically weaken this field's authority relative to abstract region
    /// signed distance.
    ///
    /// Numerical kinds are severity ordered, so requests stronger than the
    /// current kind are ignored. Downgrading to `NoClaim` also clears the
    /// abstract-distance bound; the finite nominal-field bound is retained.
    /// This infallible operation lets a representation-specific validator add
    /// knowledge such as a heuristic mesh sign without risking authority
    /// laundering or a partially applied transaction.
    #[must_use]
    pub fn downgrade_abstract_distance_authority(
        &mut self,
        requested: NumericalKind,
    ) -> NumericalKind {
        self.abstract_distance_kind = self.abstract_distance_kind.max(requested);
        if self.abstract_distance_kind == NumericalKind::NoClaim {
            self.abstract_distance_bound = None;
        }
        self.abstract_distance_kind
    }

    /// Largest actual representable grid gaps per axis.
    #[must_use]
    pub fn steps(&self) -> [f64; 3] {
        self.h
    }

    /// Re-sample every grid sample whose position lies inside `region`
    /// (inflated by one spline-support cell so reconstruction around the
    /// edit seam sees fresh data) — the incremental path for locally
    /// edited sources. Samples are recomputed at EXACTLY the original
    /// positions, so an incremental update is bit-identical to a full
    /// rebuild of the same source (the converter beads' G5 law). Returns
    /// the number of samples refreshed.
    ///
    /// # Errors
    /// [`SdfBuildError`] when the dirty region is not a finite admissible box,
    /// does not intersect the field, a source sample is invalid, or
    /// cancellation is observed mid-update. Refreshes are transactional:
    /// cancellation or refusal leaves every existing sample and both authority
    /// summaries unchanged.
    pub fn resample_box(
        &mut self,
        source: &dyn Chart,
        region: fs_geom::Aabb,
        cx: &Cx<'_>,
    ) -> Result<u64, SdfBuildError> {
        let dirty = SamplingDomain::admit(self.box_, Some(region))?.bounds();
        let max_gap_upper = self.h[0].max(self.h[1]).max(self.h[2]).next_up();
        let pad = finite_positive(
            nonnegative_mul_upper(2.0, max_gap_upper),
            "dirty-region spline padding",
        )?;
        let r = SamplingDomain::admit(self.box_, Some(dirty.inflate(pad)))?.bounds();
        let node_range = |axis: usize, min: f64, max: f64| -> (u32, u32) {
            let nodes = &self.axis_nodes[axis];
            let last = nodes.len() - 1;
            let first_not_less = nodes.partition_point(|node| *node < min);
            let first_greater = nodes.partition_point(|node| *node <= max);
            let lo = first_not_less.saturating_sub(1).min(last);
            let hi = first_greater.min(last).max(lo);
            (lo as u32, hi as u32)
        };
        let (i0, i1) = node_range(0, r.min.x, r.max.x);
        let (j0, j1) = node_range(1, r.min.y, r.max.y);
        let (k0, k1) = node_range(2, r.min.z, r.max.z);
        let update_dims = [i1 - i0 + 1, j1 - j0 + 1, k1 - k0 + 1];
        let update_count = checked_sample_product(update_dims)?;
        let source_has_global_lipschitz_theorem =
            source.trace_step_claim() == TraceStepClaim::ExactDistance;
        let has_global_lipschitz_theorem =
            self.has_global_lipschitz_theorem && source_has_global_lipschitz_theorem;
        let mut samples = Vec::with_capacity(update_count);
        let mut source_lipschitz = self.source_lipschitz;
        if source_has_global_lipschitz_theorem {
            source_lipschitz = source_lipschitz.max(1.0);
        }
        let mut max_sample_abs = self.max_sample_abs;
        let mut abstract_distance_kind = self.abstract_distance_kind;
        let mut source_error_radius = self.abstract_distance_bound.map_or(0.0, |bound| {
            nonnegative_difference_upper(bound, self.nominal_field_bound)
        });
        for k in k0..=k1 {
            for j in j0..=j1 {
                cx.checkpoint().map_err(|_| SdfBuildError::Cancelled)?;
                for i in i0..=i1 {
                    let p = Point3::new(
                        self.axis_nodes[0][i as usize],
                        self.axis_nodes[1][j as usize],
                        self.axis_nodes[2][k as usize],
                    );
                    let sample = source.eval(p, cx);
                    let stored = checked_stored_sample(&sample, p)?;
                    match sample.lipschitz {
                        Some(local) if local.is_finite() && local >= 0.0 => {
                            source_lipschitz = source_lipschitz.max(local);
                        }
                        _ if source_has_global_lipschitz_theorem => {}
                        _ => return Err(SdfBuildError::NoLipschitzBound),
                    }
                    let (sample_kind, sample_radius) = sample_abstract_distance_authority(&sample);
                    abstract_distance_kind = abstract_distance_kind.max(sample_kind);
                    if let Some(radius) = sample_radius {
                        source_error_radius = source_error_radius.max(radius);
                    }
                    max_sample_abs = max_sample_abs.max(sample.signed_distance.abs());
                    samples.push(stored);
                }
            }
        }
        let support_radius = spline_support_radius(self.h)?;
        let nominal_field_bound =
            constructed_nominal_bound(source_lipschitz, support_radius, max_sample_abs)?;
        let chart_lipschitz = constructed_chart_lipschitz(
            source_lipschitz,
            nominal_field_bound,
            self.h_min[0].min(self.h_min[1]).min(self.h_min[2]),
        )?;
        if !has_global_lipschitz_theorem {
            abstract_distance_kind = abstract_distance_kind.max(NumericalKind::Estimate);
        }
        let abstract_distance_bound = composed_abstract_bound(
            &mut abstract_distance_kind,
            nominal_field_bound,
            source_error_radius,
        );
        // Do not expose a partial refresh. In particular, observe a request
        // raised by the final source evaluation before mutating the live field.
        cx.checkpoint().map_err(|_| SdfBuildError::Cancelled)?;
        debug_assert_eq!(samples.len(), update_count);
        let cells =
            (k0..=k1).flat_map(|k| (j0..=j1).flat_map(move |j| (i0..=i1).map(move |i| [i, j, k])));
        for (cell, sample) in cells.zip(samples) {
            self.field.set(cell, sample);
        }
        self.source_lipschitz = source_lipschitz;
        self.chart_lipschitz = chart_lipschitz;
        self.has_global_lipschitz_theorem = has_global_lipschitz_theorem;
        self.max_sample_abs = max_sample_abs;
        self.nominal_field_bound = nominal_field_bound;
        self.abstract_distance_kind = abstract_distance_kind;
        self.abstract_distance_bound = abstract_distance_bound;
        Ok(u64::from(update_dims[0]) * u64::from(update_dims[1]) * u64::from(update_dims[2]))
    }

    /// Quadratic B-spline basis at fractional offset `t ∈ [0,1)` for the
    /// three support samples (standard cardinal quadratic B-spline in the
    /// C¹ warped index coordinate).
    fn bspline_w(t: f64) -> [f64; 3] {
        [
            0.5 * (1.0 - t) * (1.0 - t),
            0.5 + t * (1.0 - t),
            0.5 * t * t,
        ]
    }

    /// Derivative of [`Self::bspline_w`] with respect to t.
    fn bspline_dw(t: f64) -> [f64; 3] {
        [t - 1.0, 1.0 - 2.0 * t, t]
    }

    fn sample_raw(&self, i: i64, j: i64, k: i64) -> f64 {
        let c = [
            i.clamp(0, i64::from(self.n[0]) - 1) as u32,
            j.clamp(0, i64::from(self.n[1]) - 1) as u32,
            k.clamp(0, i64::from(self.n[2]) - 1) as u32,
        ];
        f64::from(self.field.get(c))
    }

    /// Map one physical axis onto the cardinal spline index while honoring
    /// the actual representable nodes. The monotone cubic has one shared
    /// endpoint slope, `1 / max_gap`, so adjacent physical cells meet C¹ even
    /// when their floating-point gaps differ by one or more ulps.
    fn spline_axis_coordinate(&self, axis: usize, coordinate: f64) -> (f64, f64) {
        let nodes = &self.axis_nodes[axis];
        let (lower, upper, t) = locate_axis(nodes, coordinate);
        let width = nodes[upper] - nodes[lower];
        let alpha = (width / self.h[axis]).clamp(0.0, 1.0);
        let complement = 1.0 - alpha;
        let smooth = t * t * (3.0 - 2.0 * t);
        let warped = alpha * t + complement * smooth;
        let du_dt = alpha + 6.0 * complement * t * (1.0 - t);
        (lower as f64 + warped, du_dt / width)
    }

    /// Evaluate value and gradient of the C¹ triquadratic reconstruction.
    fn spline_eval(&self, x: Point3) -> (f64, Vec3) {
        // Cell-centered quadratic B-spline: base sample nearest the point,
        // fractional offset in [-0.5, 0.5) mapped to t = frac + 0.5.
        let mut base = [0i64; 3];
        let mut t = [0.0f64; 3];
        let mapped = [
            self.spline_axis_coordinate(0, x.x),
            self.spline_axis_coordinate(1, x.y),
            self.spline_axis_coordinate(2, x.z),
        ];
        for a in 0..3 {
            let c = mapped[a].0.clamp(0.0, f64::from(self.n[a] - 1));
            let b = c.round();
            base[a] = b as i64;
            t[a] = (c - b) + 0.5;
        }
        let (wx, wy, wz) = (
            Self::bspline_w(t[0]),
            Self::bspline_w(t[1]),
            Self::bspline_w(t[2]),
        );
        let (dwx, dwy, dwz) = (
            Self::bspline_dw(t[0]),
            Self::bspline_dw(t[1]),
            Self::bspline_dw(t[2]),
        );
        let mut v = 0.0f64;
        let mut g = [0.0f64; 3];
        for (dk, wzk) in wz.iter().enumerate() {
            for (dj, wyj) in wy.iter().enumerate() {
                for (di, wxi) in wx.iter().enumerate() {
                    let s = self.sample_raw(
                        base[0] + i64::try_from(di).expect("di<3") - 1,
                        base[1] + i64::try_from(dj).expect("dj<3") - 1,
                        base[2] + i64::try_from(dk).expect("dk<3") - 1,
                    );
                    v += s * wxi * wyj * wzk;
                    g[0] += s * dwx[di] * wyj * wzk;
                    g[1] += s * wxi * dwy[dj] * wzk;
                    g[2] += s * wxi * wyj * dwz[dk];
                }
            }
        }
        (
            v,
            Vec3::new(g[0] * mapped[0].1, g[1] * mapped[1].1, g[2] * mapped[2].1),
        )
    }

    /// Measured eikonal statistics over a seeded probe set (evidence for
    /// sphere-tracing safety; a MEASUREMENT, not a certificate).
    ///
    /// # Errors
    /// [`SdfBuildError::Cancelled`] when cancellation is observed at a bounded
    /// 256-probe stride or before publishing the statistics.
    pub fn measure_eikonal(
        &self,
        seed: u64,
        probes: u64,
        cx: &Cx<'_>,
    ) -> Result<EikonalStats, SdfBuildError> {
        let mut state = seed | 1;
        let mut unit = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 11) as f64) / (1u64 << 53) as f64
        };
        let (mut sum, mut max) = (0.0f64, 0.0f64);
        for probe in 0..probes {
            if probe.is_multiple_of(256) {
                cx.checkpoint().map_err(|_| SdfBuildError::Cancelled)?;
            }
            let p = Point3::new(
                stable_lerp(self.box_.min.x, self.box_.max.x, unit()),
                stable_lerp(self.box_.min.y, self.box_.max.y, unit()),
                stable_lerp(self.box_.min.z, self.box_.max.z, unit()),
            );
            let (_, g) = self.spline_eval(p);
            let dev = (g.norm() - 1.0).abs();
            sum += dev;
            max = max.max(dev);
        }
        cx.checkpoint().map_err(|_| SdfBuildError::Cancelled)?;
        Ok(EikonalStats {
            mean_abs_dev: sum / probes.max(1) as f64,
            max_abs_dev: max,
            probes,
        })
    }

    /// Mean curvature via central second differences at step `h` (an
    /// ESTIMATE — certified curvature stencils are fs-ivl-integration
    /// follow-up work; see CONTRACT no-claims).
    #[must_use]
    pub fn mean_curvature_estimate(&self, x: Point3) -> f64 {
        let h = self.h[0].min(self.h[1]).min(self.h[2]);
        let f = |p: Point3| self.spline_eval(p).0;
        let lap = (f(x.offset(Vec3::new(h, 0.0, 0.0)))
            + f(x.offset(Vec3::new(-h, 0.0, 0.0)))
            + f(x.offset(Vec3::new(0.0, h, 0.0)))
            + f(x.offset(Vec3::new(0.0, -h, 0.0)))
            + f(x.offset(Vec3::new(0.0, 0.0, h)))
            + f(x.offset(Vec3::new(0.0, 0.0, -h)))
            - 6.0 * f(x))
            / (h * h);
        0.5 * lap
    }

    /// Sphere-trace a ray only when this field carries rigorous abstract
    /// signed-distance authority. Estimate/NoClaim fields conservatively
    /// return no hit through this compatibility API; callers must not treat
    /// the nominal-field reconstruction bound as a true-SDF certificate.
    /// Rigorous steps shrink by the total abstract-distance bound and
    /// Lipschitz safety factor (plan §10.2).
    #[must_use]
    pub fn raycast(&self, origin: Point3, dir: Vec3, t_max: f64, cx: &Cx<'_>) -> Option<f64> {
        if !self.has_global_lipschitz_theorem
            || !matches!(
                self.abstract_distance_kind,
                NumericalKind::Exact | NumericalKind::Enclosure
            )
            || !origin.x.is_finite()
            || !origin.y.is_finite()
            || !origin.z.is_finite()
            || !t_max.is_finite()
            || t_max < 0.0
        {
            return None;
        }
        let bound = self.abstract_distance_bound?;
        let d = normalized_direction(dir)?;
        let (mut t, box_exit) = ray_box_interval(origin, d, self.box_)?;
        let trace_exit = box_exit.min(t_max);
        if t > trace_exit {
            return None;
        }
        let lip = self.chart_lipschitz();
        for _ in 0..10_000 {
            if cx.is_cancel_requested() || t > trace_exit {
                return None;
            }
            let p = origin.offset(d.scale(t));
            let sample = self.eval(p, cx);
            if !sample.signed_distance.is_finite()
                || !matches!(
                    sample.error.kind,
                    NumericalKind::Exact | NumericalKind::Enclosure
                )
            {
                return None;
            }
            let sd = sample.signed_distance;
            if sd <= bound {
                return Some(t);
            }
            // Safe step: the field moves at most `lip` per unit distance,
            // and lies within `bound` of the true SDF.
            let step = ((sd - bound) / lip).max(1e-7);
            let next = t + step;
            if !step.is_finite() || !next.is_finite() || next <= t {
                return None;
            }
            t = next;
        }
        None
    }

    /// The preflighted Lipschitz candidate for the reconstruction. Callers
    /// expose it only when the source supplied a global theorem.
    fn chart_lipschitz(&self) -> f64 {
        self.chart_lipschitz
    }
}

impl Chart for TiledSdf {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let no_claim = |signed_distance: f64| ChartSample {
            signed_distance,
            gradient: None,
            lipschitz: None,
            error: NumericalCertificate::no_claim(),
        };
        if !x.x.is_finite() || !x.y.is_finite() || !x.z.is_finite() {
            return no_claim(f64::NAN);
        }
        let clamped = Point3::new(
            x.x.clamp(self.box_.min.x, self.box_.max.x),
            x.y.clamp(self.box_.min.y, self.box_.max.y),
            x.z.clamp(self.box_.min.z, self.box_.max.z),
        );
        let Some((dist_out, dist_out_upper)) =
            finite_norm_with_upper(x.x - clamped.x, x.y - clamped.y, x.z - clamped.z)
        else {
            return no_claim(f64::NAN);
        };
        let (v, g) = self.spline_eval(clamped);
        if !v.is_finite() {
            return no_claim(v);
        }
        let bound = self.abstract_distance_bound.unwrap_or(0.0);
        let abstract_error = |nominal: f64, radius: f64| {
            let Some((lo, hi)) = self
                .abstract_distance_bound
                .and_then(|_| centered_interval_hull(v, radius, nominal))
            else {
                return NumericalCertificate::no_claim();
            };
            match self.abstract_distance_kind {
                NumericalKind::Exact | NumericalKind::Enclosure => {
                    NumericalCertificate::enclosure(lo, hi)
                }
                NumericalKind::Estimate => NumericalCertificate::estimate(lo, hi),
                NumericalKind::NoClaim => NumericalCertificate::no_claim(),
            }
        };
        if dist_out == 0.0 {
            let gradient = (g.x.is_finite() && g.y.is_finite() && g.z.is_finite()).then_some(g);
            ChartSample {
                signed_distance: v,
                gradient,
                lipschitz: self
                    .has_global_lipschitz_theorem
                    .then_some(self.chart_lipschitz()),
                error: abstract_error(v, bound),
            }
        } else {
            let sd = v + dist_out;
            if !sd.is_finite() {
                return no_claim(sd);
            }
            // ExactDistance is a global unit-Lipschitz theorem. A source's
            // local `Some(L > 1)` hint is a valid loose hint, but it must not
            // replace the theorem used to propagate a rigorous enclosure.
            let propagation_lipschitz = if self.has_global_lipschitz_theorem {
                1.0
            } else {
                self.source_lipschitz
            };
            let propagation = nonnegative_mul_upper(propagation_lipschitz, dist_out_upper);
            let radius = nonnegative_add_upper(bound, propagation);
            ChartSample {
                signed_distance: sd,
                gradient: None,
                lipschitz: self
                    .has_global_lipschitz_theorem
                    .then_some(self.chart_lipschitz()),
                error: abstract_error(sd, radius),
            }
        }
    }

    fn support(&self) -> Aabb {
        self.box_
    }

    fn name(&self) -> &'static str {
        "rep-sdf/dense"
    }

    fn differentiability(&self) -> Differentiability {
        Differentiability::C1
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use asupersync::types::Budget;
    use fs_exec::{CancelGate, ExecMode, StreamKey};
    use std::sync::{
        Mutex,
        atomic::{AtomicU64, Ordering},
    };

    struct CountingPlane {
        offset: f64,
        evals: AtomicU64,
    }

    impl CountingPlane {
        fn new(offset: f64) -> Self {
            Self {
                offset,
                evals: AtomicU64::new(0),
            }
        }

        fn eval_count(&self) -> u64 {
            self.evals.load(Ordering::Relaxed)
        }
    }

    impl Chart for CountingPlane {
        fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
            self.evals.fetch_add(1, Ordering::Relaxed);
            let signed_distance = x.x - self.offset;
            ChartSample {
                signed_distance,
                gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
                lipschitz: Some(1.0),
                error: NumericalCertificate::enclosure(signed_distance, signed_distance),
            }
        }

        fn support(&self) -> Aabb {
            Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
        }

        fn name(&self) -> &'static str {
            "test/counting-plane"
        }

        fn trace_step_claim(&self) -> TraceStepClaim {
            TraceStepClaim::ExactDistance
        }
    }

    struct UnderstatedExactPlane;

    impl Chart for UnderstatedExactPlane {
        fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
            ChartSample {
                signed_distance: x.x,
                gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
                lipschitz: Some(0.0),
                error: NumericalCertificate::exact(x.x),
            }
        }

        fn support(&self) -> Aabb {
            Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
        }

        fn name(&self) -> &'static str {
            "test/understated-exact-plane"
        }

        fn trace_step_claim(&self) -> TraceStepClaim {
            TraceStepClaim::ExactDistance
        }
    }

    struct OverstatedExactPlane;

    impl Chart for OverstatedExactPlane {
        fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
            ChartSample {
                signed_distance: x.x,
                gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
                lipschitz: Some(2.0),
                error: NumericalCertificate::exact(x.x),
            }
        }

        fn support(&self) -> Aabb {
            Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
        }

        fn name(&self) -> &'static str {
            "test/overstated-exact-plane"
        }

        fn trace_step_claim(&self) -> TraceStepClaim {
            TraceStepClaim::ExactDistance
        }
    }

    #[derive(Default)]
    struct RecordingExactPlane {
        points: Mutex<Vec<Point3>>,
    }

    impl Chart for RecordingExactPlane {
        fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
            self.points.lock().expect("recording lock").push(x);
            ChartSample {
                signed_distance: x.x,
                gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
                lipschitz: Some(1.0),
                error: NumericalCertificate::exact(x.x),
            }
        }

        fn support(&self) -> Aabb {
            Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
        }

        fn name(&self) -> &'static str {
            "test/recording-exact-plane"
        }

        fn trace_step_claim(&self) -> TraceStepClaim {
            TraceStepClaim::ExactDistance
        }
    }

    struct TinyTranslatedExactPlane;

    impl TinyTranslatedExactPlane {
        const OFFSET: f64 = 1.0e-40;
        const EXTENT: f64 = 4.0e-48;
    }

    impl Chart for TinyTranslatedExactPlane {
        fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
            let signed_distance = x.x - Self::OFFSET;
            ChartSample {
                signed_distance,
                gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
                lipschitz: Some(1.0),
                error: NumericalCertificate::exact(signed_distance),
            }
        }

        fn support(&self) -> Aabb {
            Aabb::new(
                Point3::new(Self::OFFSET - Self::EXTENT, -Self::EXTENT, -Self::EXTENT),
                Point3::new(Self::OFFSET + Self::EXTENT, Self::EXTENT, Self::EXTENT),
            )
        }

        fn name(&self) -> &'static str {
            "test/tiny-translated-exact-plane"
        }

        fn trace_step_claim(&self) -> TraceStepClaim {
            TraceStepClaim::ExactDistance
        }
    }

    struct CancellingPlane<'a> {
        gate: &'a CancelGate,
        evals: AtomicU64,
    }

    impl Chart for CancellingPlane<'_> {
        fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
            if self.evals.fetch_add(1, Ordering::Relaxed) == 0 {
                self.gate.request();
            }
            let signed_distance = x.x - 0.25;
            ChartSample {
                signed_distance,
                gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
                lipschitz: Some(1.0),
                error: NumericalCertificate::enclosure(signed_distance, signed_distance),
            }
        }

        fn support(&self) -> Aabb {
            Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
        }

        fn name(&self) -> &'static str {
            "test/cancelling-plane"
        }
    }

    fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                gate,
                arena,
                StreamKey {
                    seed: 0xD3E5,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn stored_samples(sdf: &TiledSdf) -> Vec<f32> {
        let mut samples = Vec::new();
        for k in 0..sdf.n[2] {
            for j in 0..sdf.n[1] {
                for i in 0..sdf.n[0] {
                    samples.push(sdf.field.get([i, j, k]));
                }
            }
        }
        samples
    }

    #[test]
    fn dense_build_observes_cancellation_without_source_cooperation() {
        let gate = CancelGate::new();
        gate.request();
        let source = CountingPlane::new(0.0);
        with_cx(&gate, |cx| {
            let error = TiledSdf::build(&source, 0.5, cx).expect_err("cancelled build");
            assert_eq!(error, SdfBuildError::Cancelled);
        });
        assert_eq!(
            source.eval_count(),
            0,
            "checkpoint precedes chart evaluation"
        );
    }

    #[test]
    fn dense_lattice_is_representable_before_evaluation() {
        let gate = CancelGate::new();
        let source = CountingPlane::new(0.0);
        let translated = 1.0_f64;
        let two_ulps = translated.next_up().next_up();
        let box_ = Aabb::new(
            Point3::new(translated, translated, translated),
            Point3::new(two_ulps, two_ulps, two_ulps),
        );
        with_cx(&gate, |cx| {
            let error = TiledSdf::build_in_domain(&source, box_, 1.0, cx)
                .expect_err("four nominal nodes collapse onto adjacent floats");
            match error {
                SdfBuildError::DenseLatticeUnrepresentable {
                    axis,
                    index,
                    previous_bits,
                    coordinate_bits,
                } => {
                    assert_eq!(axis, 0);
                    assert!(index > 0);
                    assert_eq!(previous_bits, coordinate_bits);
                }
                other => panic!("unexpected error: {other:?}"),
            }
        });
        assert_eq!(source.eval_count(), 0, "invalid step refuses before eval");
        assert!(matches!(
            constructed_chart_lipschitz(f64::MAX, f64::MAX, f64::MIN_POSITIVE),
            Err(SdfBuildError::InvalidDerivedLipschitz { .. })
        ));
        let min = -0.5 * f64::MAX;
        let max = 0.5 * f64::MAX;
        let (nodes, min_gap, max_gap) = build_axis_nodes(min, max, 4, 0).expect("axis");
        assert!(nodes.iter().all(|node| node.is_finite()));
        assert!(nodes.windows(2).all(|pair| pair[0] < pair[1]));
        assert_eq!(nodes[3].to_bits(), max.to_bits());
        assert!(min_gap > 0.0 && max_gap >= min_gap);
    }

    #[test]
    fn translated_few_ulp_lattice_is_actual_strict_and_does_not_launder_authority() {
        let gate = CancelGate::new();
        with_cx(&gate, |cx| {
            let min = 1.0_f64;
            let mut max = min;
            for _ in 0..8 {
                max = max.next_up();
            }
            let box_ = Aabb::new(Point3::new(min, min, min), Point3::new(max, max, max));
            let mut sdf = TiledSdf::build_in_domain(&OverstatedExactPlane, box_, max - min, cx)
                .expect("eight ulps admit four strict actual nodes");
            assert_eq!(sdf.abstract_distance_kind(), NumericalKind::Enclosure);
            for nodes in &sdf.axis_nodes {
                assert_eq!(nodes.len(), 4);
                assert!(nodes.windows(2).all(|pair| pair[0] < pair[1]));
                assert!(
                    nodes
                        .windows(2)
                        .any(|pair| (pair[1] - pair[0]).to_bits() != sdf.h[0].to_bits()),
                    "fixture must exercise a nonuniform representable lattice"
                );
            }
            let mut coordinate = min;
            loop {
                let point = Point3::new(coordinate, coordinate, coordinate);
                let sample = sdf.eval(point, cx);
                let truth = OverstatedExactPlane.eval(point, cx).signed_distance;
                assert_eq!(sample.error.kind, NumericalKind::Enclosure);
                assert!(
                    sample.error.lo <= truth && truth <= sample.error.hi,
                    "truth {truth:?} escaped {:?} at {coordinate:?}",
                    sample.error
                );
                if coordinate.to_bits() == max.to_bits() {
                    break;
                }
                coordinate = coordinate.next_up();
            }

            let recorder = RecordingExactPlane::default();
            let refreshed = sdf
                .resample_box(&recorder, box_, cx)
                .expect("full translated refresh");
            assert_eq!(refreshed, 64);
            let recorded = recorder.points.lock().expect("recording lock");
            assert_eq!(recorded.len(), 64);
            for point in recorded.iter() {
                for (axis, coordinate) in [point.x, point.y, point.z].into_iter().enumerate() {
                    assert!(
                        sdf.axis_nodes[axis]
                            .iter()
                            .any(|node| node.to_bits() == coordinate.to_bits()),
                        "resample must reuse stored axis node {coordinate:?} on axis {axis}"
                    );
                }
            }
        });
    }

    #[test]
    fn exact_distance_theorem_has_unit_lipschitz_floor() {
        let gate = CancelGate::new();
        with_cx(&gate, |cx| {
            let sdf = TiledSdf::build(&UnderstatedExactPlane, 0.5, cx).expect("build");
            assert_eq!(sdf.abstract_distance_kind(), NumericalKind::Enclosure);
            assert!(sdf.source_lipschitz >= 1.0);
            assert!(sdf.chart_lipschitz >= 1.0);
            assert!(
                sdf.nominal_field_bound >= 2.0 * sdf.h[0].max(sdf.h[1]).max(sdf.h[2]),
                "Some(0) must not understate the exact-distance reconstruction band"
            );
            assert!(
                sdf.eval(Point3::new(0.0, 0.0, 0.0), cx)
                    .lipschitz
                    .is_some_and(|bound| bound >= 1.0)
            );
        });
    }

    #[test]
    fn f32_underflow_error_is_inside_the_advertised_enclosure() {
        let gate = CancelGate::new();
        with_cx(&gate, |cx| {
            let source = TinyTranslatedExactPlane;
            let sdf = TiledSdf::build(&source, 1.0e-48, cx).expect("build");
            let half_min_subnormal = 0.5 * f64::from(f32::from_bits(1));
            assert!(sdf.nominal_field_bound() >= half_min_subnormal);
            assert_eq!(sdf.abstract_distance_kind(), NumericalKind::Enclosure);

            let mut worst_storage_error = 0.0f64;
            for k in 0..sdf.n[2] {
                for j in 0..sdf.n[1] {
                    for i in 0..sdf.n[0] {
                        let point = Point3::new(
                            sdf.axis_nodes[0][i as usize],
                            sdf.axis_nodes[1][j as usize],
                            sdf.axis_nodes[2][k as usize],
                        );
                        let truth = source.eval(point, cx).signed_distance;
                        let stored = f64::from(sdf.field.get([i, j, k]));
                        worst_storage_error = worst_storage_error.max((stored - truth).abs());
                        let reconstructed = sdf.eval(point, cx);
                        assert_eq!(reconstructed.error.kind, NumericalKind::Enclosure);
                        assert!(
                            reconstructed.error.lo <= truth && truth <= reconstructed.error.hi,
                            "truth {truth:e} escaped [{:e}, {:e}]",
                            reconstructed.error.lo,
                            reconstructed.error.hi
                        );
                    }
                }
            }
            let h_max = sdf.h[0].max(sdf.h[1]).max(sdf.h[2]);
            let legacy_relative_only = 2.0 * h_max + sdf.max_sample_abs * f64::from(f32::EPSILON);
            assert!(
                worst_storage_error > legacy_relative_only,
                "fixture must exercise the old subnormal underbound"
            );
            assert!(worst_storage_error <= sdf.nominal_field_bound());
        });
    }

    #[test]
    fn abstract_distance_authority_can_only_be_weakened() {
        let gate = CancelGate::new();
        with_cx(&gate, |cx| {
            let mut sdf = TiledSdf::build(&CountingPlane::new(0.0), 0.5, cx).expect("build");
            assert_eq!(sdf.abstract_distance_kind(), NumericalKind::Enclosure);
            assert_eq!(
                sdf.eval(Point3::new(f64::MAX, 0.0, 0.0), cx).error.kind,
                NumericalKind::NoClaim,
                "overflowed outside-support interval fails closed"
            );
            let rigorous_bound = sdf.abstract_distance_bound();
            assert_eq!(
                sdf.downgrade_abstract_distance_authority(NumericalKind::Exact),
                NumericalKind::Enclosure,
                "an Exact request cannot strengthen an enclosure"
            );
            assert_eq!(sdf.abstract_distance_bound(), rigorous_bound);
            assert_eq!(
                sdf.downgrade_abstract_distance_authority(NumericalKind::Estimate),
                NumericalKind::Estimate
            );
            assert_eq!(sdf.abstract_distance_bound(), rigorous_bound);
            assert_eq!(
                sdf.downgrade_abstract_distance_authority(NumericalKind::Enclosure),
                NumericalKind::Estimate,
                "a later enclosure request cannot undo a downgrade"
            );
            assert_eq!(
                sdf.downgrade_abstract_distance_authority(NumericalKind::NoClaim),
                NumericalKind::NoClaim
            );
            assert!(sdf.abstract_distance_bound().is_none());
            assert_eq!(
                sdf.downgrade_abstract_distance_authority(NumericalKind::Exact),
                NumericalKind::NoClaim
            );
            assert_eq!(
                sdf.eval(Point3::new(0.0, 0.0, 0.0), cx).error.kind,
                NumericalKind::NoClaim
            );
        });
    }

    #[test]
    fn exact_distance_outside_enclosure_uses_unit_theorem_and_outward_arithmetic() {
        let gate = CancelGate::new();
        with_cx(&gate, |cx| {
            let sdf = TiledSdf::build(&OverstatedExactPlane, 0.5, cx).expect("build");
            assert!(
                sdf.source_lipschitz >= 2.0,
                "fixture retains the loose hint"
            );
            for point in [
                Point3::new(3.0, 0.0, 0.0),
                Point3::new(0.25, 4.0, -3.0),
                Point3::new(-3.0, 2.0, 2.0),
            ] {
                let sample = sdf.eval(point, cx);
                let truth = OverstatedExactPlane.eval(point, cx).signed_distance;
                assert_eq!(sample.error.kind, NumericalKind::Enclosure);
                assert!(sample.error.lo <= truth && truth <= sample.error.hi);
                assert!(
                    sample.error.lo <= sample.signed_distance
                        && sample.signed_distance <= sample.error.hi,
                    "published nominal must remain inside its own certificate hull"
                );
            }
            for point in [
                Point3::new(f64::NAN, 0.0, 0.0),
                Point3::new(f64::INFINITY, 0.0, 0.0),
                Point3::new(f64::MAX, f64::MAX, 0.0),
            ] {
                assert_eq!(
                    sdf.eval(point, cx).error.kind,
                    NumericalKind::NoClaim,
                    "nonfinite query or overflowed distance arithmetic must fail closed"
                );
            }
            let entry = sdf
                .raycast(
                    Point3::new(-100.0, 0.0, 0.0),
                    Vec3::new(1.0, 0.0, 0.0),
                    200.0,
                    cx,
                )
                .expect("ray reaches the finite field domain");
            assert!(
                entry > 90.0,
                "outside rays must enter the stored AABB before sampling, got {entry}"
            );
        });
    }

    #[test]
    fn eikonal_measurement_observes_cancellation() {
        let gate = CancelGate::new();
        with_cx(&gate, |cx| {
            let sdf = TiledSdf::build(&CountingPlane::new(0.0), 0.5, cx).expect("build");
            gate.request();
            assert_eq!(
                sdf.measure_eikonal(1, 10_000, cx)
                    .expect_err("cancelled probes"),
                SdfBuildError::Cancelled
            );
        });
    }

    #[test]
    fn cancelled_resample_leaves_field_bitwise_unchanged() {
        let gate = CancelGate::new();
        with_cx(&gate, |cx| {
            let mut sdf = TiledSdf::build(&CountingPlane::new(0.0), 0.5, cx).expect("build");
            let before = stored_samples(&sdf);
            let before_nominal = sdf.nominal_field_bound().to_bits();
            let before_kind = sdf.abstract_distance_kind();
            let before_abstract = sdf.abstract_distance_bound().map(f64::to_bits);
            let region = sdf.support();
            let source = CancellingPlane {
                gate: &gate,
                evals: AtomicU64::new(0),
            };
            let error = sdf
                .resample_box(&source, region, cx)
                .expect_err("source requests cancellation during staging");
            assert_eq!(error, SdfBuildError::Cancelled);
            assert_eq!(stored_samples(&sdf), before);
            assert_eq!(sdf.nominal_field_bound().to_bits(), before_nominal);
            assert_eq!(sdf.abstract_distance_kind(), before_kind);
            assert_eq!(
                sdf.abstract_distance_bound().map(f64::to_bits),
                before_abstract
            );
        });
    }
}
