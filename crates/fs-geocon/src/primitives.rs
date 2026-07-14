//! The value/derivative primitives: thickness aggregation, draft
//! angles, envelopes, and certified volume. Smooth forms exist for the
//! optimizer; certified forms exist only where an explicit theorem supports
//! them. Authority travels with every returned value and is never inferred
//! from a hard sampled reduction.

use fs_evidence::NumericalKind;
use fs_exec::Cx;
use fs_geom::{
    Aabb, Chart, ChartSample, Point3, SamplingDomain, SamplingDomainError, TraceStepClaim, Vec3,
};
use fs_query::{QueryError, thickness_at, thickness_at_clipped};

/// Thickness aggregation over boundary samples.
#[derive(Debug, Clone)]
pub struct ThicknessReport {
    /// Smooth soft-minimum: the C¹ optimizer value (`≥ hard_min`,
    /// converging down as p grows).
    pub soft_min: f64,
    /// The hard sampled minimum of local estimates. This is not a certified
    /// global thickness bound.
    pub hard_min: f64,
    /// Numerical authority inherited from the generic local thickness oracle;
    /// currently always [`NumericalKind::Estimate`].
    pub authority: NumericalKind,
    /// Indices of samples violating the requirement (LOCALIZATION).
    pub violating: Vec<usize>,
    /// Samples the oracle skipped (medial degeneracies), counted.
    pub skipped: u32,
}

/// Smooth minimum thickness over `samples` with mean p-norm
/// aggregation: `T_soft = (mean(t_i^{-p}))^{-1/p}` — a smooth
/// OVER-approximation of the minimum that converges DOWN to it as `p`
/// grows (exact when samples are uniform, which keeps lever
/// derivatives clean). The optimizer differentiates `soft_min`; the
/// sampled comparison value is `hard_min` and the localized estimated
/// violation list — the two are reported side by side and both retain
/// [`NumericalKind::Estimate`] authority.
///
/// # Errors
/// [`fs_query::QueryError`] teaching errors carried through.
pub fn min_thickness_soft(
    chart: &dyn Chart,
    samples: &[Point3],
    required: f64,
    p: f64,
    cx: &Cx<'_>,
) -> Result<ThicknessReport, QueryError> {
    min_thickness_soft_impl(chart, samples, required, p, None, cx)
}

/// Smooth minimum thickness over `samples`, restricted to an explicit finite
/// `clip`. This is the local counterpart of [`min_thickness_soft`]; it does
/// not turn the clipped result into a global minimum-thickness claim.
///
/// # Errors
/// [`QueryError::SamplingDomain`] when the clip cannot resolve the chart's
/// support into a usable finite volume, plus the local query errors propagated
/// by [`thickness_at_clipped`].
pub fn min_thickness_soft_clipped(
    chart: &dyn Chart,
    samples: &[Point3],
    required: f64,
    p: f64,
    clip: Aabb,
    cx: &Cx<'_>,
) -> Result<ThicknessReport, QueryError> {
    min_thickness_soft_impl(chart, samples, required, p, Some(clip), cx)
}

fn min_thickness_soft_impl(
    chart: &dyn Chart,
    samples: &[Point3],
    required: f64,
    p: f64,
    clip: Option<Aabb>,
    cx: &Cx<'_>,
) -> Result<ThicknessReport, QueryError> {
    SamplingDomain::admit(chart.support(), clip)?;
    cx.checkpoint().map_err(|_| QueryError::Cancelled)?;
    if !required.is_finite() || required < 0.0 {
        return Err(QueryError::InvalidThicknessArithmetic {
            reason: "the required thickness must be finite and non-negative",
        });
    }
    if !p.is_finite() || p <= 0.0 {
        return Err(QueryError::InvalidThicknessArithmetic {
            reason: "the soft-min exponent must be finite and positive",
        });
    }
    let mut inv_sum = 0.0;
    let mut count = 0u32;
    let mut hard_min = f64::INFINITY;
    let mut violating = Vec::new();
    let mut skipped = 0u32;
    for (i, &s) in samples.iter().enumerate() {
        cx.checkpoint().map_err(|_| QueryError::Cancelled)?;
        if !(s.x.is_finite() && s.y.is_finite() && s.z.is_finite()) {
            return Err(QueryError::InvalidThicknessSample {
                at: [s.x, s.y, s.z],
            });
        }
        let thickness = match clip {
            Some(domain) => thickness_at_clipped(chart, s, domain, cx),
            None => thickness_at(chart, s, cx),
        };
        match thickness {
            Ok(t) => {
                if !t.value.is_finite() || t.value <= 0.0 || t.authority != NumericalKind::Estimate
                {
                    return Err(QueryError::InvalidThicknessArithmetic {
                        reason: "the local thickness oracle returned an invalid or unexpectedly authoritative estimate",
                    });
                }
                let inverse_power = t.value.powf(-p);
                if !inverse_power.is_finite() || inverse_power <= 0.0 {
                    return Err(QueryError::InvalidThicknessArithmetic {
                        reason: "a finite local thickness overflowed the soft-min power",
                    });
                }
                inv_sum += inverse_power;
                if !inv_sum.is_finite() || inv_sum <= 0.0 {
                    return Err(QueryError::InvalidThicknessArithmetic {
                        reason: "the soft-min accumulator is not finite and positive",
                    });
                }
                count = count
                    .checked_add(1)
                    .ok_or(QueryError::InvalidThicknessArithmetic {
                        reason: "the successful thickness-sample count overflowed",
                    })?;
                hard_min = hard_min.min(t.value);
                if t.value < required {
                    violating.push(i);
                }
            }
            Err(
                QueryError::NoGradient { .. }
                | QueryError::NotOnBoundary { .. }
                | QueryError::NoOppositeWall,
            ) => {
                skipped = skipped
                    .checked_add(1)
                    .ok_or(QueryError::InvalidThicknessArithmetic {
                        reason: "the skipped thickness-sample count overflowed",
                    })?;
            }
            Err(error) => return Err(error),
        }
    }
    if count == 0 {
        return Err(QueryError::NoThicknessSamples { skipped });
    }
    let soft_min = (inv_sum / f64::from(count)).powf(-1.0 / p);
    if !soft_min.is_finite() || soft_min <= 0.0 || !hard_min.is_finite() || hard_min <= 0.0 {
        return Err(QueryError::InvalidThicknessArithmetic {
            reason: "the aggregated thickness estimate is not finite and positive",
        });
    }
    cx.checkpoint().map_err(|_| QueryError::Cancelled)?;
    Ok(ThicknessReport {
        soft_min,
        hard_min,
        authority: NumericalKind::Estimate,
        violating,
        skipped,
    })
}

/// Draft-angle assessment against a pull direction.
#[derive(Debug, Clone)]
pub struct DraftReport {
    /// Smooth penalty: mean of squared hinges `max(sinα − n·d, 0)²`
    /// over the assessed samples (C¹ in the normals).
    pub penalty: f64,
    /// EXACT violating regions: sample indices with insufficient draft.
    pub violating: Vec<usize>,
    /// Undercuts (normals pointing AGAINST the pull): worse than mere
    /// low draft; flagged separately.
    pub undercuts: Vec<usize>,
    /// Worst deficit `sinα − n·d` observed.
    pub worst_deficit: f64,
}

/// Assess draft for the mold half pulled along `pull` (unit): surface
/// normals must satisfy `n·pull ≥ sin(min_draft)`. Samples whose
/// normals oppose the pull are undercuts. Samples nearly perpendicular
/// to the pull's mirror-half (`n·pull < −cos_tolerance`) belong to the
/// other mold half and are skipped — the v1 parting model is the plane
/// perpendicular to the pull.
///
/// # Errors
/// [`fs_query::QueryError::NoGradient`] where the chart has no normal.
pub fn draft_violations(
    chart: &dyn Chart,
    samples: &[Point3],
    pull: Vec3,
    min_draft: f64,
    cx: &Cx<'_>,
) -> Result<DraftReport, fs_query::QueryError> {
    let pn = pull.norm().max(1e-300);
    let d = pull.scale(1.0 / pn);
    let sin_a = min_draft.sin();
    let mut penalty = 0.0;
    let mut violating = Vec::new();
    let mut undercuts = Vec::new();
    let mut worst = 0.0f64;
    let mut assessed = 0u32;
    for (i, &s) in samples.iter().enumerate() {
        let sample = chart.eval(s, cx);
        let Some(g) = sample.gradient else {
            return Err(fs_query::QueryError::NoGradient {
                at: [s.x, s.y, s.z],
            });
        };
        let n = g.scale(1.0 / g.norm().max(1e-300));
        let nd = n.dot(d);
        if nd < -0.5 {
            continue; // the other mold half's face
        }
        assessed += 1;
        let deficit = sin_a - nd;
        if deficit > 0.0 {
            if nd < -1e-9 {
                undercuts.push(i);
            } else {
                violating.push(i);
            }
            penalty += deficit * deficit;
            worst = worst.max(deficit);
        }
    }
    penalty /= f64::from(assessed.max(1));
    Ok(DraftReport {
        penalty,
        violating,
        undercuts,
        worst_deficit: worst,
    })
}

/// Envelope containment assessment.
#[derive(Debug, Clone)]
pub struct EnvelopeReport {
    /// Sampled worst signed distance of the design boundary into the
    /// forbidden side (`> 0` means violation).
    pub worst: f64,
    /// Smooth log-sum-exp aggregate: `≥ worst` (conservative), within
    /// `ln(n)/β` of it (the C¹ value the optimizer differentiates).
    pub soft_worst: f64,
    /// Violating sample indices.
    pub violating: Vec<usize>,
}

/// Containment: every design-boundary sample must satisfy
/// `φ_allowed ≤ 0` (inside the allowed region). For keep-outs, pass
/// the keep-out's COMPLEMENT semantics by supplying `flip = true`
/// (violation when the sample is INSIDE the keep-out).
pub fn envelope_violation(
    allowed: &dyn Chart,
    design_boundary: &[Point3],
    beta: f64,
    flip: bool,
    cx: &Cx<'_>,
) -> EnvelopeReport {
    let mut worst = f64::NEG_INFINITY;
    let mut violating = Vec::new();
    // Sum-form log-sum-exp: (1/β)·ln(Σ exp(β·g_i)) ≥ max(g_i) — a
    // CONSERVATIVE smooth upper bound, so driving the soft value to 0
    // drives the true worst to 0 (never stops short).
    let mut acc = 0.0;
    let mut max_g = f64::NEG_INFINITY;
    let gs: Vec<f64> = design_boundary
        .iter()
        .map(|&p| {
            let sd = allowed.eval(p, cx).signed_distance;
            if flip { -sd } else { sd }
        })
        .collect();
    for g in &gs {
        max_g = max_g.max(*g);
    }
    for (i, &g) in gs.iter().enumerate() {
        worst = worst.max(g);
        if g > 0.0 {
            violating.push(i);
        }
        acc += ((g - max_g) * beta).exp();
    }
    let soft_worst = if gs.is_empty() {
        0.0
    } else {
        max_g + acc.ln() / beta
    };
    EnvelopeReport {
        worst,
        soft_worst,
        violating,
    }
}

/// A rigorous volume enclosure.
#[derive(Debug, Clone, Copy)]
pub struct VolumeEnclosure {
    /// Certain lower bound (sure-inside cells).
    pub lo: f64,
    /// Certain upper bound (lower + the uncertainty band).
    pub hi: f64,
    /// Requested maximum grid step.
    pub h: f64,
}

/// Deterministic upper bound on cells admitted by either volume sampler.
pub const VOLUME_MAX_CELLS: u64 = 16_777_216;

/// Structured volume-integration failure. Every grid-preflight variant is
/// returned before the first chart evaluation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VolumeError {
    /// The explicit integration box was not a usable finite 3-D domain.
    SamplingDomain(SamplingDomainError),
    /// A caller-supplied spacing or smoothing width was not finite and
    /// strictly positive.
    InvalidSpacing {
        /// Which argument failed validation.
        field: &'static str,
        /// The rejected value.
        value: f64,
    },
    /// An axis would need more cells than a `u64` can represent.
    CellCountOverflow {
        /// Axis index in x/y/z order.
        axis: usize,
        /// Finite admitted span on that axis.
        span: f64,
        /// Requested maximum cell spacing.
        h: f64,
    },
    /// Checked multiplication of the three axis counts overflowed.
    CellProductOverflow {
        /// Per-axis cell counts involved in the overflow.
        dims: [u64; 3],
    },
    /// The requested grid exceeds the deterministic integration-work cap.
    WorkLimit {
        /// Per-axis cell counts.
        dims: [u64; 3],
        /// Total cells requested.
        need: u128,
        /// Deterministic cell cap.
        cap: u64,
    },
    /// The admitted partition has no finite positive representable cell or
    /// total volume.
    InvalidCellMeasure {
        /// Actual cell widths after partitioning each span.
        widths: [f64; 3],
    },
    /// A normalized cell center or its conservative cell radius could not be
    /// represented.
    NonRepresentableCellGeometry {
        /// Cell index in x/y/z order.
        index: [u64; 3],
        /// Best available center coordinates.
        point: Point3,
    },
    /// A chart returned a non-finite nominal field value.
    InvalidSample {
        /// Point at which evaluation failed.
        point: Point3,
        /// Raw bits of the rejected value.
        value_bits: u64,
    },
    /// A rigorous enclosure requires the chart's global exact-distance
    /// theorem; a sample-local Lipschitz number is insufficient.
    UncertifiedChart {
        /// The weaker theorem actually advertised by the chart.
        claim: TraceStepClaim,
    },
    /// A purported exact-distance chart returned weak or malformed evidence
    /// at one integration cell.
    InvalidCertificate {
        /// Point at which the certificate failed closed.
        point: Point3,
        /// Certificate authority actually returned.
        kind: NumericalKind,
        /// Raw lower-endpoint bits.
        lo_bits: u64,
        /// Raw upper-endpoint bits.
        hi_bits: u64,
    },
    /// Integration observed cancellation at a bounded polling point.
    Cancelled {
        /// Cells fully evaluated before cancellation was observed.
        completed_cells: u64,
    },
}

impl core::fmt::Display for VolumeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::SamplingDomain(error) => write!(f, "{error}"),
            Self::InvalidSpacing { field, value } => write!(
                f,
                "volume integration refused: `{field}` must be finite and strictly positive, got {value}"
            ),
            Self::CellCountOverflow { axis, span, h } => write!(
                f,
                "volume integration refused: axis {axis} span {span} at spacing {h} has no representable cell count"
            ),
            Self::CellProductOverflow { dims } => write!(
                f,
                "volume integration refused: cell dimensions {dims:?} overflow the checked work count"
            ),
            Self::WorkLimit { dims, need, cap } => write!(
                f,
                "volume integration refused: grid {dims:?} requires {need} cells, exceeding the {cap} cell cap; coarsen h or shrink the domain"
            ),
            Self::InvalidCellMeasure { widths } => write!(
                f,
                "volume integration refused: cell widths {widths:?} have no finite positive representable volume"
            ),
            Self::NonRepresentableCellGeometry { index, point } => write!(
                f,
                "volume integration refused: cell {index:?} has non-representable center/radius geometry at {point:?}"
            ),
            Self::InvalidSample { point, value_bits } => write!(
                f,
                "volume integration refused: chart sample at {point:?} is non-finite (f64 bits {value_bits:#018x})"
            ),
            Self::UncertifiedChart { claim } => write!(
                f,
                "certified volume refused: chart advertises {claim:?}, but v1 requires the global ExactDistance theorem"
            ),
            Self::InvalidCertificate {
                point,
                kind,
                lo_bits,
                hi_bits,
            } => write!(
                f,
                "certified volume refused: chart returned malformed or non-rigorous {kind:?} evidence at {point:?} (bounds {lo_bits:#018x}..={hi_bits:#018x})"
            ),
            Self::Cancelled { completed_cells } => write!(
                f,
                "volume integration cancelled after {completed_cells} completed cells"
            ),
        }
    }
}

impl core::error::Error for VolumeError {}

impl From<SamplingDomainError> for VolumeError {
    fn from(error: SamplingDomainError) -> Self {
        Self::SamplingDomain(error)
    }
}

#[derive(Debug, Clone, Copy)]
struct VolumeGrid {
    bounds: Aabb,
    /// Outward enclosures of each exact endpoint difference.
    span_bounds: [[f64; 2]; 3],
    dims: [u64; 3],
    /// Outward enclosures of the exact uniform cell widths.
    width_bounds: [[f64; 2]; 3],
    /// Nearest-rounded cell volume used only by the explicitly non-certifying
    /// smooth estimator.
    cell_volume_estimate: f64,
    /// Outward enclosure of one exact uniform cell volume.
    cell_volume_bounds: [f64; 2],
}

fn finite_positive(value: f64, field: &'static str) -> Result<f64, VolumeError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(VolumeError::InvalidSpacing { field, value })
    }
}

fn checked_axis_cells(span_upper: f64, h: f64, axis: usize) -> Result<u64, VolumeError> {
    let cells = (span_upper / h).next_up().ceil().max(1.0);
    if !cells.is_finite() || cells >= u64::MAX as f64 {
        return Err(VolumeError::CellCountOverflow {
            axis,
            span: span_upper,
            h,
        });
    }
    Ok(cells as u64)
}

fn positive_product_down(lhs: f64, rhs: f64) -> f64 {
    (lhs * rhs).next_down().max(0.0)
}

fn positive_product_up(lhs: f64, rhs: f64) -> f64 {
    (lhs * rhs).next_up()
}

fn positive_quotient_down(numerator: f64, denominator: f64) -> f64 {
    (numerator / denominator).next_down().max(0.0)
}

fn positive_quotient_up(numerator: f64, denominator: f64) -> f64 {
    (numerator / denominator).next_up()
}

fn positive_sum_up(lhs: f64, rhs: f64) -> f64 {
    (lhs + rhs).next_up()
}

fn volume_grid(domain: Aabb, h: f64) -> Result<VolumeGrid, VolumeError> {
    let h = finite_positive(h, "h")?;
    let admitted = SamplingDomain::admit(domain, None)?;
    let bounds = admitted.bounds();
    let nearest_spans = admitted.spans();
    let nearest_spans = [nearest_spans.x, nearest_spans.y, nearest_spans.z];
    let span_bounds = nearest_spans.map(|span| [span.next_down().max(0.0), span.next_up()]);
    if span_bounds
        .iter()
        .any(|bounds| !bounds[1].is_finite() || bounds[1] <= 0.0)
    {
        return Err(VolumeError::InvalidCellMeasure {
            widths: nearest_spans,
        });
    }
    let dims = [
        checked_axis_cells(span_bounds[0][1], h, 0)?,
        checked_axis_cells(span_bounds[1][1], h, 1)?,
        checked_axis_cells(span_bounds[2][1], h, 2)?,
    ];
    let need = dims
        .iter()
        .try_fold(1u128, |product, &dim| product.checked_mul(u128::from(dim)))
        .ok_or(VolumeError::CellProductOverflow { dims })?;
    if need > u128::from(VOLUME_MAX_CELLS) {
        return Err(VolumeError::WorkLimit {
            dims,
            need,
            cap: VOLUME_MAX_CELLS,
        });
    }

    let width_bounds = core::array::from_fn(|axis| {
        let count = dims[axis] as f64;
        [
            positive_quotient_down(span_bounds[axis][0], count),
            positive_quotient_up(span_bounds[axis][1], count),
        ]
    });
    let widths = core::array::from_fn(|axis| nearest_spans[axis] / dims[axis] as f64);
    let cell_volume_estimate = widths[0] * widths[1] * widths[2];
    let cell_volume_lo = positive_product_down(
        positive_product_down(width_bounds[0][0], width_bounds[1][0]),
        width_bounds[2][0],
    );
    let cell_volume_hi = positive_product_up(
        positive_product_up(width_bounds[0][1], width_bounds[1][1]),
        width_bounds[2][1],
    );
    let total_volume_hi = positive_product_up(cell_volume_hi, need as f64);
    if width_bounds
        .iter()
        .any(|width| !width[1].is_finite() || width[1] <= 0.0)
        || !cell_volume_estimate.is_finite()
        || cell_volume_estimate <= 0.0
        || !cell_volume_hi.is_finite()
        || cell_volume_hi <= 0.0
        || !total_volume_hi.is_finite()
        || total_volume_hi <= 0.0
    {
        return Err(VolumeError::InvalidCellMeasure { widths });
    }

    Ok(VolumeGrid {
        bounds,
        span_bounds,
        dims,
        width_bounds,
        cell_volume_estimate,
        cell_volume_bounds: [cell_volume_lo, cell_volume_hi],
    })
}

fn cell_center_and_radius(
    grid: &VolumeGrid,
    index: [u64; 3],
) -> Result<(Point3, f64), VolumeError> {
    let mins = [grid.bounds.min.x, grid.bounds.min.y, grid.bounds.min.z];
    let maxs = [grid.bounds.max.x, grid.bounds.max.y, grid.bounds.max.z];
    let mut centers = [0.0; 3];
    let mut axis_radii = [0.0; 3];
    for axis in 0..3 {
        // VOLUME_MAX_CELLS is below 2^53, so both integer conversions and
        // `index + 1/2` are exact before the outward-rounded division.
        let numerator = index[axis] as f64 + 0.5;
        let denominator = grid.dims[axis] as f64;
        let fraction = [
            (numerator / denominator).next_down().max(0.0),
            (numerator / denominator).next_up().min(1.0),
        ];
        let offset_lo = positive_product_down(grid.span_bounds[axis][0], fraction[0]);
        let offset_hi = positive_product_up(grid.span_bounds[axis][1], fraction[1]);
        let center_lo = (mins[axis] + offset_lo)
            .next_down()
            .clamp(mins[axis], maxs[axis]);
        let center_hi = (mins[axis] + offset_hi)
            .next_up()
            .clamp(mins[axis], maxs[axis]);
        let center = f64::midpoint(center_lo, center_hi);
        if !center.is_finite() || center < center_lo || center > center_hi {
            let point = Point3::new(centers[0], centers[1], centers[2]);
            return Err(VolumeError::NonRepresentableCellGeometry { index, point });
        }
        centers[axis] = center;
        let center_error = (center - center_lo)
            .abs()
            .max((center_hi - center).abs())
            .next_up();
        let half_width = positive_product_up(0.5, grid.width_bounds[axis][1]);
        axis_radii[axis] = positive_sum_up(center_error, half_width);
        if !axis_radii[axis].is_finite() || axis_radii[axis] <= 0.0 {
            let point = Point3::new(centers[0], centers[1], centers[2]);
            return Err(VolumeError::NonRepresentableCellGeometry { index, point });
        }
    }
    let point = Point3::new(centers[0], centers[1], centers[2]);
    // L1 dominates Euclidean distance and avoids under-rounded sqrt/square
    // arithmetic. Each addition is directed upward.
    let radius = positive_sum_up(positive_sum_up(axis_radii[0], axis_radii[1]), axis_radii[2]);
    if !radius.is_finite() || radius <= 0.0 {
        return Err(VolumeError::NonRepresentableCellGeometry { index, point });
    }
    Ok((point, radius))
}

fn rigorous_distance_bounds(
    chart: &dyn Chart,
    point: Point3,
    sample: &ChartSample,
    completed_cells: u64,
    cx: &Cx<'_>,
) -> Result<(f64, f64), VolumeError> {
    if !sample.signed_distance.is_finite() {
        return Err(VolumeError::InvalidSample {
            point,
            value_bits: sample.signed_distance.to_bits(),
        });
    }
    volume_checkpoint(cx, completed_cells)?;
    let certificate = chart.trace_value_enclosure(point, sample, cx);
    volume_checkpoint(cx, completed_cells)?;
    let valid = certificate.lo.is_finite()
        && certificate.hi.is_finite()
        && certificate.lo <= sample.signed_distance
        && sample.signed_distance <= certificate.hi
        && match certificate.kind {
            NumericalKind::Exact => {
                certificate.lo.to_bits() == sample.signed_distance.to_bits()
                    && certificate.hi.to_bits() == sample.signed_distance.to_bits()
            }
            NumericalKind::Enclosure => true,
            NumericalKind::Estimate | NumericalKind::NoClaim => false,
        };
    if !valid {
        return Err(VolumeError::InvalidCertificate {
            point,
            kind: certificate.kind,
            lo_bits: certificate.lo.to_bits(),
            hi_bits: certificate.hi.to_bits(),
        });
    }
    Ok((certificate.lo, certificate.hi))
}

fn volume_checkpoint(cx: &Cx<'_>, completed_cells: u64) -> Result<(), VolumeError> {
    cx.checkpoint()
        .map_err(|_| VolumeError::Cancelled { completed_cells })
}

/// Certified volume over an EXPLICIT integration domain (fixed independently
/// of design levers, so lever derivatives see the shape change, not grid
/// realignment). `h` is the requested maximum cell width; each axis is
/// partitioned into `ceil(span / h)` cells, sampled at normalized centers.
/// Cells farther inside than an outward-rounded L1 radius enclosing their
/// exact partition cell are SURELY inside; cells within that distance of the
/// zero set form the uncertainty band. The
/// true volume lies in `[lo, hi]` for charts advertising the global
/// [`TraceStepClaim::ExactDistance`] theorem and returning a rigorous finite
/// enclosure at every cell center. Local Lipschitz samples and weak distance
/// estimates cannot authorize this result.
///
/// # Errors
/// [`VolumeError`] reports invalid/unbounded domains, invalid spacing,
/// excessive or unrepresentable work, invalid chart samples, or cancellation.
pub fn volume_certified(
    chart: &dyn Chart,
    domain: &fs_geom::Aabb,
    h: f64,
    cx: &Cx<'_>,
) -> Result<VolumeEnclosure, VolumeError> {
    let grid = volume_grid(*domain, h)?;
    let claim = chart.trace_step_claim();
    if claim != TraceStepClaim::ExactDistance {
        return Err(VolumeError::UncertifiedChart { claim });
    }
    let mut sure = 0u64;
    let mut band = 0u64;
    let mut completed_cells = 0u64;
    for i in 0..grid.dims[0] {
        for j in 0..grid.dims[1] {
            for k in 0..grid.dims[2] {
                if completed_cells.is_multiple_of(256) {
                    volume_checkpoint(cx, completed_cells)?;
                }
                let index = [i, j, k];
                let (p, radius) = cell_center_and_radius(&grid, index)?;
                let sample = chart.eval(p, cx);
                volume_checkpoint(cx, completed_cells)?;
                let (lo, hi) = rigorous_distance_bounds(chart, p, &sample, completed_cells, cx)?;
                completed_cells += 1;
                let cell_hi = positive_sum_up(hi, radius);
                let cell_lo = (lo - radius).next_down();
                if cell_hi <= 0.0 {
                    sure += 1;
                } else if cell_lo < 0.0 {
                    band += 1;
                }
            }
        }
    }
    volume_checkpoint(cx, completed_cells)?;
    Ok(VolumeEnclosure {
        lo: positive_product_down(sure as f64, grid.cell_volume_bounds[0]),
        hi: positive_product_up((sure + band) as f64, grid.cell_volume_bounds[1]),
        h,
    })
}

/// Smoothed volume: `Σ cell_volume·σ(−φ/ε)` over the same normalized
/// maximum-step grid as [`volume_certified`], using the logistic mollifier.
/// This is the C¹ value whose lever derivative matches the Hadamard shape
/// derivative on fixtures (the battery's validation).
///
/// # Errors
/// [`VolumeError`] reports invalid/unbounded domains, invalid spacing or
/// smoothing width, excessive or unrepresentable work, invalid chart samples,
/// or cancellation.
pub fn volume_smooth(
    chart: &dyn Chart,
    domain: &fs_geom::Aabb,
    h: f64,
    epsilon: f64,
    cx: &Cx<'_>,
) -> Result<f64, VolumeError> {
    finite_positive(epsilon, "epsilon")?;
    let grid = volume_grid(*domain, h)?;
    let mut acc = 0.0;
    let mut completed_cells = 0u64;
    for i in 0..grid.dims[0] {
        for j in 0..grid.dims[1] {
            for k in 0..grid.dims[2] {
                if completed_cells.is_multiple_of(256) {
                    volume_checkpoint(cx, completed_cells)?;
                }
                let (p, _) = cell_center_and_radius(&grid, [i, j, k])?;
                let sd = chart.eval(p, cx).signed_distance;
                volume_checkpoint(cx, completed_cells)?;
                if !sd.is_finite() {
                    return Err(VolumeError::InvalidSample {
                        point: p,
                        value_bits: sd.to_bits(),
                    });
                }
                acc += grid.cell_volume_estimate / (1.0 + (sd / epsilon).exp());
                completed_cells += 1;
            }
        }
    }
    volume_checkpoint(cx, completed_cells)?;
    Ok(acc)
}

#[cfg(test)]
mod volume_rounding_tests {
    use super::*;

    #[test]
    fn nextafter_cell_geometry_and_measure_are_outward_enclosed() {
        let x_min = f64::from_bits(0x4330_0000_0000_0000); // Exactly 2^52.
        let x_max = x_min.next_up();
        let domain = Aabb::new(Point3::new(x_min, -1.0, 0.0), Point3::new(x_max, 1.0, 1.0));
        let grid = volume_grid(domain, 2.0).expect("finite nextafter domain");
        assert_eq!(grid.dims, [1, 2, 1]);
        let (center, radius) = cell_center_and_radius(&grid, [0, 0, 0])
            .expect("unrepresentable ideal midpoint gets an enclosing stored center");
        assert!(center.x == x_min || center.x == x_max);
        assert!(center.y.is_finite() && center.z.is_finite() && radius.is_finite());

        let nearest_measure = grid.cell_volume_estimate;
        assert!(grid.cell_volume_bounds[0] <= nearest_measure);
        assert!(nearest_measure <= grid.cell_volume_bounds[1]);
        assert!(grid.cell_volume_bounds[0] < grid.cell_volume_bounds[1]);
        assert!(
            radius >= positive_product_up(0.5, grid.width_bounds[0][1]),
            "the radius must include at least the outward half-width before center-rounding error"
        );
    }
}
