//! CONVERTER SDF → NURBS (plan §7.3 edge 4, bead wqd.12; [F] — behind
//! the `nurbs-refit` feature until its Gauntlet tier is green): spline
//! RE-FITTING with thin-plate smoothing and an honest sampled error
//! report. THE STRATEGIC ROLE (§7.2): Booleans route through F-rep/SDF
//! and re-fit to splines when a spline chart is required — this edge is
//! what makes the honest NURBS Boolean policy work.
//!
//! v1 pipeline (star-shaped domains): radial projection through a
//! (u, v) direction grid finds retained sign-bracket targets by BISECTION ON
//! THE CALLER'S SCALAR CLOSURE. Without an admitted continuity/root witness a
//! target is not authoritative evidence of a zero set or source surface;
//! tensor-product B-spline least squares with discrete thin-plate
//! (control-lattice Laplacian) regularization; exact G⁰ seam closure by
//! control-column tying, G¹ measured.
//!
//! Error honesty: the report keeps `max sampled |f(S(u,v))|` separate from the
//! geometric probe-spacing estimate `L_S·h_probe`, where `h_probe` is the
//! retained probe grid covering radius in parameter space. Adding them is dimensionally and
//! analytically justified only when the caller separately proves a compatible
//! unit-Lipschitz field model. This generic closure API carries neither that
//! authority, a metric-error bound, nor directed-rounding evidence, and `|f|`
//! is not generically an upper geometric distance. Retained projection-target
//! coverage is likewise sampled. A future admitted-field API may promote these
//! with the required units, metric regularity, and interval evidence. Large
//! localized paired-parameter fit residuals produce structured warnings with
//! locations. Thin geometry is one possible cause, but smoothing, inadequate
//! patch density, conditioning, or caller-field behavior can produce the same
//! signal; features missed by every ray remain outside this API's visibility.

use crate::NurbsError;
use crate::basis::{AdmittedKnotVector, KnotVector};
use crate::closest::norm3;
use crate::surface::{AdmittedNurbsSurface, NurbsSurface};
use core::mem::size_of;
use fs_exec::Cx;
use fs_math::det;

/// The fitting knobs (the ErrBudget-style trade: patch density vs
/// fidelity, priced by the router cost model).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RefitConfig {
    /// Control-net size along u (the seam direction).
    pub nu: usize,
    /// Control-net size along v.
    pub nv: usize,
    /// B-spline degree (both directions).
    pub degree: usize,
    /// Thin-plate (bending-energy) weight.
    pub lambda: f64,
    /// Sample-grid resolution along u.
    pub samples_u: usize,
    /// Sample-grid resolution along v.
    pub samples_v: usize,
    /// Residuals above this trigger localized fit-residual warnings.
    pub warn_residual: f64,
    /// Dense-probe resolution per axis for sampled field residuals and the
    /// separate geometric probe-spacing estimate.
    pub probe: usize,
}

impl Default for RefitConfig {
    fn default() -> Self {
        RefitConfig {
            nu: 12,
            nv: 12,
            degree: 3,
            lambda: 1e-4,
            samples_u: 36,
            samples_v: 36,
            warn_residual: 5e-2,
            probe: 96,
        }
    }
}

/// Static safety envelope for the legacy closure-based refit API. A successor
/// budgeted/cancellable API will make these caller-visible ledger values.
const REFIT_MAX_ALLOC_BYTES: usize = 256 * 1024 * 1024;
const REFIT_MAX_PROBE_POINTS: usize = 4 * 1024 * 1024;
const REFIT_MAX_WORK_UNITS: u128 = 1_000_000_000;
const REFIT_CANCELLATION_STRIDE: usize = 64;

fn refit_structure_error(what: impl Into<String>) -> NurbsError {
    NurbsError::Structure { what: what.into() }
}

fn refit_allocation_error(what: &'static str) -> NurbsError {
    NurbsError::Domain {
        what: format!("{what} allocation was refused"),
    }
}

fn checked_refit_work_product(values: &[u128], stage: &str) -> Result<u128, NurbsError> {
    values.iter().try_fold(1u128, |work, value| {
        work.checked_mul(*value).ok_or_else(|| {
            refit_structure_error(format!("refit {stage} work accounting overflows u128"))
        })
    })
}

fn checked_refit_work_sum(values: &[u128], stage: &str) -> Result<u128, NurbsError> {
    values.iter().try_fold(0u128, |work, value| {
        work.checked_add(*value).ok_or_else(|| {
            refit_structure_error(format!("refit {stage} work accounting overflows u128"))
        })
    })
}

fn try_vec_with_capacity<T>(capacity: usize, what: &'static str) -> Result<Vec<T>, NurbsError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(capacity)
        .map_err(|_| refit_allocation_error(what))?;
    Ok(values)
}

fn try_filled_vec<T: Copy>(len: usize, value: T, what: &'static str) -> Result<Vec<T>, NurbsError> {
    let mut values = try_vec_with_capacity(len, what)?;
    for _ in 0..len {
        values.push(value);
    }
    Ok(values)
}

fn try_filled_matrix<T: Copy>(
    rows: usize,
    cols: usize,
    value: T,
    what: &'static str,
) -> Result<Vec<Vec<T>>, NurbsError> {
    rows.checked_mul(cols)
        .ok_or_else(|| refit_allocation_error(what))?;
    let mut matrix = try_vec_with_capacity(rows, what)?;
    for _ in 0..rows {
        matrix.push(try_filled_vec(cols, value, what)?);
    }
    Ok(matrix)
}

/// Validate dimensions and derive all allocation/work sizes before the first
/// field evaluation or allocation.
fn validate_refit_request(
    center: [f64; 3],
    r_max: f64,
    config: &RefitConfig,
) -> Result<(usize, usize, usize), NurbsError> {
    if center.iter().any(|coordinate| !coordinate.is_finite()) {
        return Err(refit_structure_error("refit center must be finite"));
    }
    if !r_max.is_finite() || r_max <= 0.0 {
        return Err(refit_structure_error(
            "refit radial extent must be finite and positive",
        ));
    }
    if config.degree == 0
        || config.nu < 2
        || config.nv < 2
        || config.degree >= config.nu
        || config.degree >= config.nv
    {
        return Err(refit_structure_error(
            "refit needs degree >= 1 and at least degree+1 control points on each axis",
        ));
    }
    if config.samples_u == 0 || config.samples_v == 0 || config.probe == 0 {
        return Err(refit_structure_error(
            "refit sample and probe resolutions must be positive",
        ));
    }
    if !config.lambda.is_finite()
        || config.lambda < 0.0
        || !config.warn_residual.is_finite()
        || config.warn_residual < 0.0
    {
        return Err(refit_structure_error(
            "refit lambda and warning threshold must be finite and non-negative",
        ));
    }

    let control_points = config
        .nu
        .checked_mul(config.nv)
        .ok_or_else(|| refit_structure_error("refit control-grid size overflow"))?;
    let sample_points = config
        .samples_u
        .checked_mul(config.samples_v)
        .ok_or_else(|| refit_structure_error("refit sample-grid size overflow"))?;
    let minimum_probe = config
        .samples_u
        .max(config.samples_v)
        .checked_mul(2)
        .ok_or_else(|| refit_structure_error("refit probe-axis size overflow"))?;
    let probe = config.probe.max(minimum_probe);
    let probe_points = probe
        .checked_mul(probe)
        .ok_or_else(|| refit_structure_error("refit probe-grid size overflow"))?;
    if probe_points > REFIT_MAX_PROBE_POINTS {
        return Err(refit_structure_error(format!(
            "refit probe grid {probe_points} exceeds static cap {REFIT_MAX_PROBE_POINTS}"
        )));
    }

    let row_scalars = sample_points
        .checked_mul(control_points)
        .ok_or_else(|| refit_structure_error("refit sample-matrix size overflow"))?;
    let dense_scalars = control_points
        .checked_mul(control_points)
        .ok_or_else(|| refit_structure_error("refit normal-matrix size overflow"))?;
    let bytes_for = |count: usize, element_bytes: usize| {
        count
            .checked_mul(element_bytes)
            .ok_or_else(|| refit_structure_error("refit allocation-byte estimate overflow"))
    };
    let mut allocation_bytes = 0usize;
    let mut add_bytes = |bytes: usize| -> Result<(), NurbsError> {
        allocation_bytes = allocation_bytes
            .checked_add(bytes)
            .ok_or_else(|| refit_structure_error("refit aggregate allocation size overflow"))?;
        Ok(())
    };
    // Conservative simultaneously-live envelope. Include nested `Vec` headers
    // and every sample-sized side buffer rather than counting only f64 matrix
    // payloads. Allocator metadata is implementation-defined and remains outside
    // this deterministic process cap.
    add_bytes(bytes_for(row_scalars, size_of::<f64>())?)?;
    add_bytes(bytes_for(sample_points, size_of::<Vec<f64>>())?)?;
    add_bytes(bytes_for(sample_points, size_of::<[f64; 3]>())?)?; // targets
    add_bytes(bytes_for(sample_points, size_of::<[f64; 2]>())?)?; // uvs
    add_bytes(bytes_for(
        sample_points,
        size_of::<LocalizedFitResidualWarning>(),
    )?)?;
    add_bytes(bytes_for(dense_scalars, size_of::<f64>())?)?;
    add_bytes(bytes_for(control_points, size_of::<Vec<f64>>())?)?;
    add_bytes(bytes_for(control_points, size_of::<f64>())?)?; // rhs
    add_bytes(bytes_for(control_points, size_of::<[f64; 3]>())?)?; // net
    add_bytes(bytes_for(config.nu, size_of::<Vec<[f64; 3]>>())?)?;
    add_bytes(bytes_for(control_points, size_of::<f64>())?)?; // weights
    add_bytes(bytes_for(config.nu, size_of::<Vec<f64>>())?)?;
    add_bytes(bytes_for(control_points, size_of::<[f64; 4]>())?)?; // surface cpw
    add_bytes(bytes_for(config.nu, size_of::<Vec<[f64; 4]>>())?)?;
    let knot_overhead = config
        .degree
        .checked_add(1)
        .and_then(|value| value.checked_mul(2))
        .ok_or_else(|| refit_structure_error("refit knot allocation size overflow"))?;
    let knot_scalars = config
        .nu
        .checked_add(config.nv)
        .and_then(|value| value.checked_add(knot_overhead))
        .ok_or_else(|| refit_structure_error("refit knot allocation size overflow"))?;
    add_bytes(bytes_for(knot_scalars, size_of::<f64>())?)?;
    if allocation_bytes > REFIT_MAX_ALLOC_BYTES {
        return Err(refit_structure_error(format!(
            "refit allocation estimate {allocation_bytes} bytes exceeds static cap {REFIT_MAX_ALLOC_BYTES}"
        )));
    }

    let active_basis = config
        .degree
        .checked_add(1)
        .and_then(|value| value.checked_mul(value))
        .ok_or_else(|| refit_structure_error("refit active-basis size overflow"))?;
    let basis_triangle = config
        .degree
        .checked_mul(config.degree + 1)
        .map(|product| product / 2)
        .ok_or_else(|| refit_structure_error("refit basis-work size overflow"))?;
    let knots_per_axis = config
        .degree
        .checked_add(1)
        .and_then(|order| config.nu.checked_add(order).map(|u_knots| (order, u_knots)))
        .and_then(|(order, u_knots)| {
            config
                .nv
                .checked_add(order)
                .and_then(|v_knots| u_knots.checked_add(v_knots))
        })
        .ok_or_else(|| refit_structure_error("refit probe knot-scan size overflow"))?;
    let assembly_work = checked_refit_work_product(
        &[
            sample_points as u128,
            active_basis as u128,
            control_points as u128,
        ],
        "sample assembly",
    )?;
    let factor_work = checked_refit_work_product(
        &[
            control_points as u128,
            control_points as u128,
            control_points as u128,
        ],
        "normal factorization",
    )?;
    let rhs_and_report_work = checked_refit_work_product(
        &[sample_points as u128, control_points as u128, 6],
        "right-hand-side and report",
    )?;
    let triangular_solve_work = checked_refit_work_product(
        &[control_points as u128, control_points as u128, 3],
        "triangular solve",
    )?;
    let projection_evaluations =
        checked_refit_work_product(&[sample_points as u128, 42], "radial projection")?;
    // The report binds one admitted immutable surface across all probes. Keep
    // charging the former owning-evaluator scan on every probe as conservative
    // headroom: this validate-once migration must not silently loosen the
    // legacy process ceiling. Tightening that ceiling belongs to the successor
    // caller-budgeted API. The arbitrary closure's own cost remains outside
    // this legacy static model.
    let basis_triangle_work =
        checked_refit_work_product(&[basis_triangle as u128, 2], "probe basis triangle")?;
    let active_basis_work =
        checked_refit_work_product(&[active_basis as u128, 8], "probe active basis")?;
    let per_probe_work = checked_refit_work_sum(
        &[
            control_points as u128,
            knots_per_axis as u128,
            basis_triangle_work,
            active_basis_work,
            16,
        ],
        "per-probe",
    )?;
    let probe_work =
        checked_refit_work_product(&[probe_points as u128, per_probe_work], "probe grid")?;
    let total_work = checked_refit_work_sum(
        &[
            assembly_work,
            factor_work,
            rhs_and_report_work,
            triangular_solve_work,
            projection_evaluations,
            probe_work,
        ],
        "aggregate",
    )?;
    if total_work > REFIT_MAX_WORK_UNITS {
        return Err(refit_structure_error(format!(
            "refit work estimate {total_work} exceeds static cap {REFIT_MAX_WORK_UNITS}"
        )));
    }
    Ok((control_points, sample_points, probe))
}

/// A localized paired-parameter fit residual above the configured threshold.
/// This diagnoses a retained mismatch, not its cause: thin geometry,
/// smoothing, global under-resolution, conditioning, or caller-field behavior
/// can all produce the same observation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LocalizedFitResidualWarning {
    /// Parameter location of the offending sample.
    pub uv: [f64; 2],
    /// World-space location.
    pub point: [f64; 3],
    /// The residual left behind.
    pub residual: f64,
}

/// The fit report: sampled and analytic-model estimates with no continuum or
/// geometric-distance certificate.
#[derive(Debug, Clone, PartialEq)]
pub struct RefitReport {
    /// RMS fit residual over the sample grid.
    pub rms_residual: f64,
    /// Worst fit residual.
    pub max_residual: f64,
    /// Worst paired-parameter residual from a retained radial sign-bracket
    /// target to the fitted spline point at the same `(u,v)`. This is an upper
    /// witness for target-to-spline point-set distance, not that distance
    /// itself. A generic closure also does not prove that targets lie on a
    /// source surface or even on a continuous field's zero set.
    pub projected_target_to_spline_sampled: f64,
    /// Sampled maximum of `|field(S(u,v))|` (spline → source-field direction).
    pub spline_to_field_sampled: f64,
    /// Geometric probe-spacing estimate `(L_u + L_v) · h_probe` from the fitted
    /// surface to the nearest retained probe in parameter space. This has
    /// position units, not arbitrary field-value units, and therefore is not
    /// added to `spline_to_field_sampled` by this generic API. Ordinary f64
    /// arithmetic makes it an estimate rather than an outward enclosure.
    pub spline_probe_spacing_estimate: f64,
    /// Numerically evaluated hodograph Lipschitz estimate used above. The
    /// analytic formula is sound for this non-rational unit-weight surface,
    /// but the f64 result is not outward-rounded.
    pub spline_lipschitz_estimate: f64,
    /// Max G¹ seam deviation (angle proxy: 1 − cos between u-tangents
    /// across the seam); G⁰ is exact by construction. Positive infinity is the
    /// explicit no-claim value when either tangent direction is undefined.
    pub seam_g1_max: f64,
    /// Retained seam samples at which at least one tangent direction was
    /// exactly degenerate. Such samples force `seam_g1_max = +∞` rather than
    /// silently looking perfect.
    pub seam_g1_degenerate_samples: usize,
    /// Whether the current seam-direction diagnostic excludes `v=0` and
    /// `v=1`. Pole tangent directions require a separate chart-aware audit;
    /// the exclusion is machine-visible rather than implied coverage.
    pub seam_g1_excludes_v_endpoints: bool,
    /// Localized fit-residual warnings. Empty means only that no retained
    /// paired-parameter sample exceeded the configured threshold; it neither
    /// diagnoses fit quality between samples nor proves a geometric feature
    /// absent.
    pub warnings: Vec<LocalizedFitResidualWarning>,
}

/// The refit result.
#[derive(Debug)]
pub struct Refit {
    /// The fitted surface (u closed by control tying, v open).
    pub surface: NurbsSurface<f64>,
    /// The honest report.
    pub report: RefitReport,
}

/// Direction of the (u, v) spherical parameterization. `v` runs
/// SOUTH → NORTH (φ = π(1 − v)) so the fitted surface's du × dv
/// normals point OUTWARD — the orientation the signed chart
/// presentation (and the sheaf comparison against source fields)
/// relies on.
fn direction(u: f64, v: f64) -> [f64; 3] {
    let theta = 2.0 * std::f64::consts::PI * u;
    let phi = std::f64::consts::PI * (1.0 - v);
    let sin_phi = det::sin(phi);
    [
        sin_phi * det::cos(theta),
        sin_phi * det::sin(theta),
        det::cos(phi),
    ]
}

/// Transactional outcome of a cancellation-aware radial projection.
///
/// Cancellation never exposes the partially narrowed sign bracket or a
/// provisional radius.
#[must_use]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RadialProjectionRun {
    /// All 40 fixed bisection steps and the final publication checkpoint
    /// completed.
    Complete {
        /// Scalar `r` in `center + r * direction`.
        radius: f64,
    },
    /// Cancellation was observed before a field evaluation or publication.
    Cancelled,
}

/// Transactional outcome of cancellation-aware normal-matrix assembly.
///
/// Cancellation drops the partially initialized or accumulated dense matrix;
/// callers receive either the complete matrix or no matrix at all.
#[must_use]
#[derive(Debug, Clone, PartialEq)]
pub enum RefitNormalAssemblyRun {
    /// The complete `B^T B + lambda L^T L` matrix is safe to publish.
    Complete {
        /// Dense row-major normal matrix.
        matrix: Vec<Vec<f64>>,
    },
    /// Cancellation was observed before publication.
    Cancelled,
}

/// Transactional outcome of cancellation-aware dense Cholesky factorization.
///
/// The primitive consumes its input matrix. Cancellation or error drops that
/// storage; only a complete factor can be returned to the caller.
#[must_use]
#[derive(Debug, Clone, PartialEq)]
pub enum RefitNormalFactorRun {
    /// The complete lower-triangular factor is safe to publish. Entries above
    /// the diagonal retain their input values and are not part of the factor.
    Complete {
        /// Dense storage containing the Cholesky factor in its lower triangle.
        factor: Vec<Vec<f64>>,
    },
    /// Cancellation was observed before publication.
    Cancelled,
}

/// Transactional outcome of a cancellation-aware solve against a completed
/// refit Cholesky factor.
///
/// The primitive consumes its right-hand side. Cancellation or error drops the
/// partial forward/back substitution state; only a complete solution escapes.
#[must_use]
#[derive(Debug, Clone, PartialEq)]
pub enum RefitNormalSolveRun {
    /// Both triangular substitutions and the publication checkpoint completed.
    Complete {
        /// Complete solution vector.
        solution: Vec<f64>,
    },
    /// Cancellation was observed before publication.
    Cancelled,
}

/// Transactional outcome of a cancellation-aware measured surface Lipschitz
/// estimate.
///
/// The returned directional values use ordinary `f64` arithmetic. They are
/// measured estimates, not outward enclosures or certified Lipschitz bounds.
#[must_use]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RefitLipschitzEstimateRun {
    /// The complete U/V control-hodograph traversal is safe to publish.
    Complete {
        /// Estimated U-direction derivative magnitude ceiling.
        u_estimate: f64,
        /// Estimated V-direction derivative magnitude ceiling.
        v_estimate: f64,
    },
    /// Cancellation was observed before publication.
    Cancelled,
}

fn validate_radial_projection_request(
    center: [f64; 3],
    dir: [f64; 3],
    r_max: f64,
) -> Result<(), NurbsError> {
    if center.iter().any(|coordinate| !coordinate.is_finite()) {
        return Err(refit_structure_error(
            "radial projection center must be finite",
        ));
    }
    if dir.iter().any(|coordinate| !coordinate.is_finite())
        || dir.iter().all(|coordinate| *coordinate == 0.0)
    {
        return Err(refit_structure_error(
            "radial projection direction must be finite and nonzero",
        ));
    }
    if !r_max.is_finite() || r_max <= 0.0 {
        return Err(refit_structure_error(
            "radial projection extent must be finite and positive",
        ));
    }
    Ok(())
}

fn radial_field_value_with_poll(
    field: &dyn Fn([f64; 3]) -> f64,
    center: [f64; 3],
    dir: [f64; 3],
    r: f64,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<Option<f64>, NurbsError> {
    if should_cancel() {
        return Ok(None);
    }
    let point = [
        center[0] + r * dir[0],
        center[1] + r * dir[1],
        center[2] + r * dir[2],
    ];
    if point.iter().any(|coordinate| !coordinate.is_finite()) {
        return Err(refit_structure_error(
            "radial field sample point is not representable",
        ));
    }
    let value = field(point);
    if !value.is_finite() {
        return Err(refit_structure_error(format!(
            "implicit field returned non-finite value at {point:?}"
        )));
    }
    Ok(Some(value))
}

fn radial_bracket_error(dir: [f64; 3]) -> NurbsError {
    NurbsError::Structure {
        what: format!(
            "radial bracket failed along {dir:?}: refit v1 needs a star-shaped \
             domain around the given center (field(center) < 0 < field(center + r_max·dir))"
        ),
    }
}

fn project_radial_with_poll(
    field: &dyn Fn([f64; 3]) -> f64,
    center: [f64; 3],
    dir: [f64; 3],
    r_max: f64,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<RadialProjectionRun, NurbsError> {
    validate_radial_projection_request(center, dir, r_max)?;
    let (mut lo, mut hi) = (0.0f64, r_max);
    let Some(at_lo) = radial_field_value_with_poll(field, center, dir, lo, should_cancel)? else {
        return Ok(RadialProjectionRun::Cancelled);
    };
    if at_lo >= 0.0 {
        return Err(radial_bracket_error(dir));
    }
    let Some(at_hi) = radial_field_value_with_poll(field, center, dir, hi, should_cancel)? else {
        return Ok(RadialProjectionRun::Cancelled);
    };
    if at_hi <= 0.0 {
        return Err(radial_bracket_error(dir));
    }
    for _ in 0..40 {
        let mid = f64::midpoint(lo, hi);
        let Some(value) = radial_field_value_with_poll(field, center, dir, mid, should_cancel)?
        else {
            return Ok(RadialProjectionRun::Cancelled);
        };
        if value < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    if should_cancel() {
        return Ok(RadialProjectionRun::Cancelled);
    }
    Ok(RadialProjectionRun::Complete {
        radius: f64::midpoint(lo, hi),
    })
}

/// Bisect the implicit field along `center + r * direction` for a sign
/// crossing, with cancellation before every field call and final publication.
///
/// Constant-time request refusals precede the first checkpoint. The caller's
/// direction supplies the ray parameterization; this operation does not grant
/// continuity, root-existence, unit-direction, or geometric-distance
/// authority to the closure or returned scalar.
///
/// # Errors
///
/// Returns a structured [`NurbsError`] for malformed inputs, a non-finite field
/// value or sample point, or a missing strict negative-to-positive bracket.
pub fn project_radial_with_cx(
    field: &dyn Fn([f64; 3]) -> f64,
    center: [f64; 3],
    direction: [f64; 3],
    r_max: f64,
    cx: &Cx<'_>,
) -> Result<RadialProjectionRun, NurbsError> {
    let mut should_cancel = || cx.checkpoint().is_err();
    project_radial_with_poll(field, center, direction, r_max, &mut should_cancel)
}

/// Bisect the implicit field along `center + r·dir` for a sign crossing.
fn project_radial(
    field: &dyn Fn([f64; 3]) -> f64,
    center: [f64; 3],
    dir: [f64; 3],
    r_max: f64,
) -> Result<f64, NurbsError> {
    let mut never_cancel = || false;
    match project_radial_with_poll(field, center, dir, r_max, &mut never_cancel)? {
        RadialProjectionRun::Complete { radius } => Ok(radius),
        RadialProjectionRun::Cancelled => Err(refit_structure_error(
            "non-cancellable radial projection observed cancellation",
        )),
    }
}

fn refit_poll_due(
    operations_since_poll: &mut usize,
    should_cancel: &mut impl FnMut() -> bool,
) -> bool {
    *operations_since_poll += 1;
    if *operations_since_poll < REFIT_CANCELLATION_STRIDE {
        return false;
    }
    *operations_since_poll = 0;
    should_cancel()
}

fn preflight_refit_normal_assembly(
    row_count: usize,
    nu: usize,
    nv: usize,
    lambda: f64,
) -> Result<usize, NurbsError> {
    if nu == 0 || nv == 0 {
        return Err(refit_structure_error(
            "refit normal-matrix axes must be nonzero",
        ));
    }
    if !lambda.is_finite() || lambda < 0.0 {
        return Err(refit_structure_error(
            "refit normal-matrix lambda must be finite and non-negative",
        ));
    }
    let n = nu
        .checked_mul(nv)
        .ok_or_else(|| refit_structure_error("refit normal-matrix dimension overflow"))?;
    let matrix_scalars = n
        .checked_mul(n)
        .ok_or_else(|| refit_structure_error("refit normal-matrix size overflow"))?;
    let matrix_bytes = matrix_scalars
        .checked_mul(size_of::<f64>())
        .and_then(|bytes| {
            n.checked_mul(size_of::<Vec<f64>>())
                .and_then(|headers| bytes.checked_add(headers))
        })
        .ok_or_else(|| refit_structure_error("refit normal-matrix byte estimate overflow"))?;
    if matrix_bytes > REFIT_MAX_ALLOC_BYTES {
        return Err(refit_structure_error(format!(
            "refit normal-matrix allocation estimate {matrix_bytes} bytes exceeds static cap {REFIT_MAX_ALLOC_BYTES}"
        )));
    }

    let row_count = row_count as u128;
    let n_work = n as u128;
    let row_validation_work =
        checked_refit_work_product(&[row_count, n_work], "normal row validation")?;
    let matrix_initialization_work =
        checked_refit_work_product(&[n_work, n_work], "normal matrix initialization")?;
    let gram_work =
        checked_refit_work_product(&[row_count, n_work, n_work], "normal Gram assembly")?;
    // A five-point lattice stencil contributes at most 5 * 5 outer products.
    let regularization_work = checked_refit_work_product(&[n_work, 25], "normal regularization")?;
    let total_work = checked_refit_work_sum(
        &[
            row_validation_work,
            matrix_initialization_work,
            gram_work,
            regularization_work,
            2,
        ],
        "normal assembly",
    )?;
    if total_work > REFIT_MAX_WORK_UNITS {
        return Err(refit_structure_error(format!(
            "refit normal-matrix work estimate {total_work} exceeds static cap {REFIT_MAX_WORK_UNITS}"
        )));
    }
    Ok(n)
}

fn preflight_refit_normal_factor(n: usize) -> Result<(), NurbsError> {
    if n == 0 {
        return Err(refit_structure_error(
            "refit normal factorization requires a nonempty square matrix",
        ));
    }
    let n_work = n as u128;
    let shape_work = n_work;
    let value_and_symmetry_work =
        checked_refit_work_product(&[n_work, n_work, 2], "normal factor validation")?;
    let factor_work =
        checked_refit_work_product(&[n_work, n_work, n_work], "normal factorization")?;
    let total_work = checked_refit_work_sum(
        &[shape_work, value_and_symmetry_work, factor_work, 1],
        "normal factorization",
    )?;
    if total_work > REFIT_MAX_WORK_UNITS {
        return Err(refit_structure_error(format!(
            "refit normal-factor work estimate {total_work} exceeds static cap {REFIT_MAX_WORK_UNITS}"
        )));
    }
    Ok(())
}

fn preflight_refit_normal_solve(factor_rows: usize, rhs_len: usize) -> Result<(), NurbsError> {
    if rhs_len == 0 {
        return Err(refit_structure_error(
            "refit normal solve requires a nonempty right-hand side",
        ));
    }
    if factor_rows != rhs_len {
        return Err(refit_structure_error(format!(
            "refit normal solve factor has {factor_rows} rows for right-hand side length {rhs_len}"
        )));
    }
    let n_work = rhs_len as u128;
    let shape_and_rhs_work = checked_refit_work_product(&[n_work, 2], "normal solve inputs")?;
    let factor_and_solve_work =
        checked_refit_work_product(&[n_work, n_work, 3], "normal triangular solve")?;
    let total_work = checked_refit_work_sum(
        &[shape_and_rhs_work, factor_and_solve_work, 1],
        "normal solve",
    )?;
    if total_work > REFIT_MAX_WORK_UNITS {
        return Err(refit_structure_error(format!(
            "refit normal-solve work estimate {total_work} exceeds static cap {REFIT_MAX_WORK_UNITS}"
        )));
    }
    Ok(())
}

/// Factor a dense symmetric-positive-definite matrix with bounded
/// cancellation polling.
///
/// Count-derived worst-case work refusal precedes the first checkpoint. One
/// gate then spans square-shape, finite-value and exact-symmetry validation;
/// deterministic lower-triangle factorization; and final publication. The
/// primitive allocates no derived numerical payload and does not consume the
/// `Cx` budget, solve a right-hand side, or make the full refit pipeline
/// cancellation-aware.
///
/// # Errors
/// Returns a structured [`NurbsError`] for an empty or non-square matrix,
/// non-finite or asymmetric input, checked work refusal, or a non-positive or
/// non-finite pivot.
pub fn factor_refit_normal_with_cx(
    matrix: Vec<Vec<f64>>,
    cx: &Cx<'_>,
) -> Result<RefitNormalFactorRun, NurbsError> {
    let mut should_cancel = || cx.checkpoint().is_err();
    factor_refit_normal_with_poll(matrix, &mut should_cancel)
}

#[allow(clippy::needless_range_loop)]
fn factor_refit_normal_with_poll(
    mut matrix: Vec<Vec<f64>>,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<RefitNormalFactorRun, NurbsError> {
    let n = matrix.len();
    preflight_refit_normal_factor(n)?;
    if should_cancel() {
        return Ok(RefitNormalFactorRun::Cancelled);
    }

    let mut operations_since_poll = 0usize;
    for (row_index, row) in matrix.iter().enumerate() {
        if row.len() != n {
            return Err(refit_structure_error(format!(
                "refit normal-factor row {row_index} has length {}, expected {n}",
                row.len()
            )));
        }
        if refit_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(RefitNormalFactorRun::Cancelled);
        }
    }
    for (row_index, row) in matrix.iter().enumerate() {
        for value in row {
            if !value.is_finite() {
                return Err(refit_structure_error(format!(
                    "refit normal-factor row {row_index} contains a non-finite value"
                )));
            }
            if refit_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(RefitNormalFactorRun::Cancelled);
            }
        }
    }
    for i in 0..n {
        for j in 0..i {
            if matrix[i][j] != matrix[j][i] {
                return Err(refit_structure_error(format!(
                    "refit normal-factor input is asymmetric at ({i}, {j})"
                )));
            }
            if refit_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(RefitNormalFactorRun::Cancelled);
            }
        }
    }

    for i in 0..n {
        for j in 0..=i {
            let mut sum = matrix[i][j];
            let (ri, rj) = (&matrix[i], &matrix[j]);
            for (x, y) in ri[..j].iter().zip(&rj[..j]) {
                sum -= x * y;
                if refit_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(RefitNormalFactorRun::Cancelled);
                }
            }
            if i == j {
                if !sum.is_finite() || sum <= 0.0 {
                    return Err(NurbsError::Structure {
                        what: "normal equations not SPD (raise lambda or sample count)".to_string(),
                    });
                }
                matrix[i][i] = det::sqrt(sum);
            } else {
                matrix[i][j] = sum / matrix[j][j];
                if !matrix[i][j].is_finite() {
                    return Err(refit_structure_error(
                        "normal-equation factorization became non-finite",
                    ));
                }
            }
            if refit_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(RefitNormalFactorRun::Cancelled);
            }
        }
    }
    if should_cancel() {
        return Ok(RefitNormalFactorRun::Cancelled);
    }
    Ok(RefitNormalFactorRun::Complete { factor: matrix })
}

/// Dense symmetric-positive-definite Cholesky factorization. The factor is
/// shared across all three coordinate right-hand sides.
fn cholesky_factor(matrix: Vec<Vec<f64>>) -> Result<Vec<Vec<f64>>, NurbsError> {
    let mut never_cancel = || false;
    match factor_refit_normal_with_poll(matrix, &mut never_cancel)? {
        RefitNormalFactorRun::Complete { factor } => Ok(factor),
        RefitNormalFactorRun::Cancelled => Err(refit_structure_error(
            "non-cancellable refit normal factorization observed cancellation",
        )),
    }
}

/// Solve one right-hand side using a completed lower-triangular refit factor
/// with bounded cancellation polling.
///
/// Count-derived dimensions and worst-case work are refused before the first
/// checkpoint. One gate then spans factor shape/lower-triangle validation,
/// right-hand-side validation, deterministic forward/back substitution, finite
/// arithmetic checks, and final publication. The borrowed factor is never
/// modified; cancellation drops the consumed right-hand side. This primitive
/// does not consume the `Cx` budget, prove conditioning, or make the full refit
/// pipeline cancellation-aware.
///
/// # Errors
/// Returns a structured [`NurbsError`] for dimension/shape mismatch,
/// non-finite factor or right-hand-side values, non-positive diagonal entries,
/// checked work refusal, or non-finite substitution arithmetic.
pub fn solve_refit_normal_with_cx(
    factor: &[Vec<f64>],
    rhs: Vec<f64>,
    cx: &Cx<'_>,
) -> Result<RefitNormalSolveRun, NurbsError> {
    let mut should_cancel = || cx.checkpoint().is_err();
    solve_refit_normal_with_poll(factor, rhs, &mut should_cancel)
}

#[allow(clippy::needless_range_loop)]
fn solve_refit_normal_with_poll(
    factor: &[Vec<f64>],
    mut rhs: Vec<f64>,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<RefitNormalSolveRun, NurbsError> {
    let n = rhs.len();
    preflight_refit_normal_solve(factor.len(), n)?;
    if should_cancel() {
        return Ok(RefitNormalSolveRun::Cancelled);
    }

    let mut operations_since_poll = 0usize;
    for (row_index, row) in factor.iter().enumerate() {
        if row.len() != n {
            return Err(refit_structure_error(format!(
                "refit normal-solve factor row {row_index} has length {}, expected {n}",
                row.len()
            )));
        }
        if refit_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(RefitNormalSolveRun::Cancelled);
        }
    }
    for i in 0..n {
        for j in 0..=i {
            let value = factor[i][j];
            if !value.is_finite() {
                return Err(refit_structure_error(format!(
                    "refit normal-solve factor contains a non-finite lower entry at ({i}, {j})"
                )));
            }
            if i == j && value <= 0.0 {
                return Err(refit_structure_error(format!(
                    "refit normal-solve factor has a non-positive diagonal at {i}"
                )));
            }
            if refit_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(RefitNormalSolveRun::Cancelled);
            }
        }
    }
    for (index, value) in rhs.iter().enumerate() {
        if !value.is_finite() {
            return Err(refit_structure_error(format!(
                "refit normal-solve right-hand side contains a non-finite value at {index}"
            )));
        }
        if refit_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(RefitNormalSolveRun::Cancelled);
        }
    }

    for i in 0..n {
        let mut sum = rhs[i];
        for k in 0..i {
            sum -= factor[i][k] * rhs[k];
            if refit_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(RefitNormalSolveRun::Cancelled);
            }
        }
        rhs[i] = sum / factor[i][i];
        if !rhs[i].is_finite() {
            return Err(refit_structure_error(
                "normal-equation forward solve became non-finite",
            ));
        }
        if refit_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(RefitNormalSolveRun::Cancelled);
        }
    }
    for i in (0..n).rev() {
        let mut sum = rhs[i];
        for k in (i + 1)..n {
            sum -= factor[k][i] * rhs[k];
            if refit_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(RefitNormalSolveRun::Cancelled);
            }
        }
        rhs[i] = sum / factor[i][i];
        if !rhs[i].is_finite() {
            return Err(refit_structure_error(
                "normal-equation back solve became non-finite",
            ));
        }
        if refit_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(RefitNormalSolveRun::Cancelled);
        }
    }
    if should_cancel() {
        return Ok(RefitNormalSolveRun::Cancelled);
    }
    Ok(RefitNormalSolveRun::Complete { solution: rhs })
}

fn cholesky_solve_factored(factor: &[Vec<f64>], rhs: Vec<f64>) -> Result<Vec<f64>, NurbsError> {
    let mut never_cancel = || false;
    match solve_refit_normal_with_poll(factor, rhs, &mut never_cancel)? {
        RefitNormalSolveRun::Complete { solution } => Ok(solution),
        RefitNormalSolveRun::Cancelled => Err(refit_structure_error(
            "non-cancellable refit normal solve observed cancellation",
        )),
    }
}

fn open_uniform_knots(n: usize, degree: usize) -> Result<KnotVector<f64>, NurbsError> {
    let order = degree
        .checked_add(1)
        .ok_or_else(|| refit_structure_error("refit knot order overflow"))?;
    let inner = n
        .checked_sub(degree)
        .ok_or_else(|| refit_structure_error("refit degree exceeds control count"))?;
    if inner == 0 {
        return Err(refit_structure_error(
            "refit degree must be less than the control count",
        ));
    }
    let knot_count = n
        .checked_add(order)
        .ok_or_else(|| refit_structure_error("refit knot count overflow"))?;
    let mut knots = try_vec_with_capacity(knot_count, "refit knot vector")?;
    for _ in 0..order {
        knots.push(0.0);
    }
    #[allow(clippy::cast_precision_loss)]
    for k in 1..inner {
        knots.push(k as f64 / inner as f64);
    }
    knots.extend(std::iter::repeat_n(1.0, order));
    KnotVector::new(knots, degree)
}

/// Row of basis values over the whole control axis (dense, small).
fn basis_row(kv: AdmittedKnotVector<'_, f64>, t: f64) -> Result<Vec<f64>, NurbsError> {
    let (span, vals) = kv.basis(t)?;
    let n = kv.control_count();
    let mut row = try_filled_vec(n, 0.0f64, "refit dense basis row")?;
    let p = kv.degree();
    for (c, &b) in vals.iter().enumerate() {
        row[span - p + c] = b;
    }
    Ok(row)
}

fn preflight_refit_lipschitz_estimate(rows: usize, cols: usize) -> Result<(), NurbsError> {
    if rows == 0 || cols == 0 {
        return Err(refit_structure_error(
            "refit Lipschitz estimate requires a nonempty control net",
        ));
    }
    let controls =
        checked_refit_work_product(&[rows as u128, cols as u128], "Lipschitz control traversal")?;
    // Each visit performs at most two three-coordinate projections,
    // differences/norms, span divisions, and maxima. This intentionally
    // overprices the fixed f64 scalar work.
    let traversal_work = checked_refit_work_product(&[controls, 32], "Lipschitz estimate")?;
    let total_work = checked_refit_work_sum(&[traversal_work, 2], "Lipschitz estimate")?;
    if total_work > REFIT_MAX_WORK_UNITS {
        return Err(refit_structure_error(format!(
            "refit Lipschitz work estimate {total_work} exceeds static cap {REFIT_MAX_WORK_UNITS}"
        )));
    }
    Ok(())
}

/// Analytic spline Lipschitz formula from the hodograph hull. The derivative
/// curve `S'(u) = Σ Dᵢ Nᵢ,ₚ₋₁(u)` has control points
/// `Dᵢ = p·ΔCᵢ / (u_{i+p+1} − u_{i+1})`, and B-spline bases are a nonnegative
/// partition of unity, so `|S'(u)| ≤ maxᵢ‖Dᵢ‖ = L`.
///
/// The per-difference knot span `u_{i+p+1} − u_{i+1}` MUST be used: the closed
/// form `L ≤ max‖ΔC‖·(n−p)` only holds for the uniform interior span
/// `p/(n−p)`. On a clamped open-uniform knot vector the END spans collapse
/// (for ΔC₀, `u_{p+1} − u₁ = 1/(n−p)`, one interval), so `p/span = p·(n−p)` —
/// the closed form UNDER-bounds by up to a factor `p` when the largest control
/// difference sits near the clamp, which would make the estimate too tight.
/// The implementation uses ordinary f64 arithmetic and therefore returns an
/// estimate, not an outward-rounded enclosure.
///
/// Count-derived traversal-work refusal precedes the first checkpoint. One
/// gate then spans the complete admitted control-hodograph traversal,
/// finite-arithmetic checks, and final publication, polling after at most 64
/// control visits. The primitive allocates no numerical payload and does not
/// consume the `Cx` budget or grant distance, field, or no-tunneling authority.
///
/// # Errors
/// Returns a structured [`NurbsError`] for checked work refusal or non-finite
/// projection, norm, or coefficient arithmetic.
pub fn estimate_refit_surface_lipschitz_with_cx(
    surface: AdmittedNurbsSurface<'_, f64>,
    cx: &Cx<'_>,
) -> Result<RefitLipschitzEstimateRun, NurbsError> {
    let mut should_cancel = || cx.checkpoint().is_err();
    estimate_refit_surface_lipschitz_with_poll(surface, &mut should_cancel)
}

fn estimate_refit_surface_lipschitz_with_poll(
    admitted: AdmittedNurbsSurface<'_, f64>,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<RefitLipschitzEstimateRun, NurbsError> {
    let surface = admitted.source();
    let rows = surface.cpw.len();
    let cols = surface.cpw[0].len();
    preflight_refit_lipschitz_estimate(rows, cols)?;
    if should_cancel() {
        return Ok(RefitLipschitzEstimateRun::Cancelled);
    }

    let p_u = surface.knots_u.degree;
    let p_v = surface.knots_v.degree;
    let ku = &surface.knots_u.knots;
    let kv = &surface.knots_v.knots;
    let cart = |h: &[f64; 4]| [h[0] / h[3], h[1] / h[3], h[2] / h[3]];
    let dist = |a: [f64; 3], b: [f64; 3]| -> f64 { norm3([a[0] - b[0], a[1] - b[1], a[2] - b[2]]) };
    let mut lu = 0.0f64;
    let mut lv = 0.0f64;
    let mut controls_since_poll = 0usize;
    for i in 0..rows {
        for j in 0..cols {
            let c = cart(&surface.cpw[i][j]);
            if c.iter().any(|coordinate| !coordinate.is_finite()) {
                return Err(refit_structure_error(
                    "refit Lipschitz control projection became non-finite",
                ));
            }
            if i + 1 < rows {
                let dc = dist(cart(&surface.cpw[i + 1][j]), c);
                let span = ku[i + p_u + 1] - ku[i + 1];
                if span > 0.0 {
                    #[allow(clippy::cast_precision_loss)]
                    let coef = p_u as f64 * dc / span;
                    if !coef.is_finite() {
                        return Err(refit_structure_error(
                            "refit U Lipschitz estimate became non-finite",
                        ));
                    }
                    lu = lu.max(coef);
                }
            }
            if j + 1 < cols {
                let dc = dist(cart(&surface.cpw[i][j + 1]), c);
                let span = kv[j + p_v + 1] - kv[j + 1];
                if span > 0.0 {
                    #[allow(clippy::cast_precision_loss)]
                    let coef = p_v as f64 * dc / span;
                    if !coef.is_finite() {
                        return Err(refit_structure_error(
                            "refit V Lipschitz estimate became non-finite",
                        ));
                    }
                    lv = lv.max(coef);
                }
            }
            if refit_poll_due(&mut controls_since_poll, should_cancel) {
                return Ok(RefitLipschitzEstimateRun::Cancelled);
            }
        }
    }
    if should_cancel() {
        return Ok(RefitLipschitzEstimateRun::Cancelled);
    }
    Ok(RefitLipschitzEstimateRun::Complete {
        u_estimate: lu,
        v_estimate: lv,
    })
}

fn lipschitz_estimate(surface: AdmittedNurbsSurface<'_, f64>) -> Result<(f64, f64), NurbsError> {
    let mut never_cancel = || false;
    match estimate_refit_surface_lipschitz_with_poll(surface, &mut never_cancel)? {
        RefitLipschitzEstimateRun::Complete {
            u_estimate,
            v_estimate,
        } => Ok((u_estimate, v_estimate)),
        RefitLipschitzEstimateRun::Cancelled => Err(refit_structure_error(
            "non-cancellable refit Lipschitz estimate observed cancellation",
        )),
    }
}

/// Fit one scalar/vector LSQ system: `(BᵀB + λ LᵀL) c = Bᵀy` where `L`
/// is the discrete control-lattice Laplacian (thin-plate proxy).
///
/// Checked dimension, requested-output, and worst-case work refusals precede
/// the first checkpoint. One gate then covers borrowed row validation,
/// fallible matrix initialization, deterministic Gram and regularization
/// accumulation, finite-arithmetic checks, and final publication. A cancelled
/// run exposes no partial matrix. This primitive does not consume the `Cx`
/// budget or make the full refit pipeline cancellation-aware.
///
/// # Errors
/// Returns a structured [`NurbsError`] for invalid dimensions, row shapes,
/// non-finite inputs or arithmetic, checked work/retained-memory refusal, or
/// allocation refusal.
pub fn assemble_refit_normal_with_cx(
    rows_b: &[Vec<f64>],
    nu: usize,
    nv: usize,
    lambda: f64,
    cx: &Cx<'_>,
) -> Result<RefitNormalAssemblyRun, NurbsError> {
    let mut should_cancel = || cx.checkpoint().is_err();
    assemble_refit_normal_with_poll(rows_b, nu, nv, lambda, &mut should_cancel)
}

#[allow(clippy::needless_range_loop, clippy::too_many_lines)]
fn assemble_refit_normal_with_poll(
    rows_b: &[Vec<f64>],
    nu: usize,
    nv: usize,
    lambda: f64,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<RefitNormalAssemblyRun, NurbsError> {
    let n = preflight_refit_normal_assembly(rows_b.len(), nu, nv, lambda)?;
    if should_cancel() {
        return Ok(RefitNormalAssemblyRun::Cancelled);
    }

    let mut operations_since_poll = 0usize;
    for (row_index, row) in rows_b.iter().enumerate() {
        if row.len() != n {
            return Err(refit_structure_error(format!(
                "refit normal-matrix row {row_index} has length {}, expected {n}",
                row.len()
            )));
        }
        for value in row {
            if !value.is_finite() {
                return Err(refit_structure_error(format!(
                    "refit normal-matrix row {row_index} contains a non-finite value"
                )));
            }
            if refit_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(RefitNormalAssemblyRun::Cancelled);
            }
        }
    }

    let mut a = try_vec_with_capacity(n, "refit normal matrix")?;
    if should_cancel() {
        return Ok(RefitNormalAssemblyRun::Cancelled);
    }
    for _ in 0..n {
        let mut row = try_vec_with_capacity(n, "refit normal matrix")?;
        if should_cancel() {
            return Ok(RefitNormalAssemblyRun::Cancelled);
        }
        for _ in 0..n {
            row.push(0.0);
            if refit_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(RefitNormalAssemblyRun::Cancelled);
            }
        }
        a.push(row);
    }

    for row in rows_b {
        for i in 0..n {
            if row[i] == 0.0 {
                if refit_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(RefitNormalAssemblyRun::Cancelled);
                }
                continue;
            }
            for j in 0..n {
                if row[j] != 0.0 {
                    let next = a[i][j] + row[i] * row[j];
                    if !next.is_finite() {
                        return Err(refit_structure_error(
                            "refit normal-matrix Gram assembly became non-finite",
                        ));
                    }
                    a[i][j] = next;
                }
                if refit_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(RefitNormalAssemblyRun::Cancelled);
                }
            }
        }
    }
    // Thin-plate: Laplacian rows (4-neighbor) on the control lattice.
    let idx = |i: usize, j: usize| i * nv + j;
    for i in 0..nu {
        for j in 0..nv {
            let mut stencil = [(0usize, 0.0f64); 5];
            let mut stencil_len = 1usize;
            stencil[0] = (idx(i, j), 0.0);
            let mut degree = 0.0f64;
            if i > 0 {
                stencil[stencil_len] = (idx(i - 1, j), -1.0);
                stencil_len += 1;
                degree += 1.0;
            }
            if i + 1 < nu {
                stencil[stencil_len] = (idx(i + 1, j), -1.0);
                stencil_len += 1;
                degree += 1.0;
            }
            if j > 0 {
                stencil[stencil_len] = (idx(i, j - 1), -1.0);
                stencil_len += 1;
                degree += 1.0;
            }
            if j + 1 < nv {
                stencil[stencil_len] = (idx(i, j + 1), -1.0);
                stencil_len += 1;
                degree += 1.0;
            }
            stencil[0].1 = degree;
            for &(p, wp) in &stencil[..stencil_len] {
                for &(q, wq) in &stencil[..stencil_len] {
                    let next = a[p][q] + lambda * wp * wq;
                    if !next.is_finite() {
                        return Err(refit_structure_error(
                            "refit normal-matrix regularization became non-finite",
                        ));
                    }
                    a[p][q] = next;
                    if refit_poll_due(&mut operations_since_poll, should_cancel) {
                        return Ok(RefitNormalAssemblyRun::Cancelled);
                    }
                }
            }
        }
    }
    if should_cancel() {
        return Ok(RefitNormalAssemblyRun::Cancelled);
    }
    Ok(RefitNormalAssemblyRun::Complete { matrix: a })
}

fn assemble_normal(
    rows_b: &[Vec<f64>],
    nu: usize,
    nv: usize,
    lambda: f64,
) -> Result<Vec<Vec<f64>>, NurbsError> {
    let mut never_cancel = || false;
    match assemble_refit_normal_with_poll(rows_b, nu, nv, lambda, &mut never_cancel)? {
        RefitNormalAssemblyRun::Complete { matrix } => Ok(matrix),
        RefitNormalAssemblyRun::Cancelled => Err(refit_structure_error(
            "non-cancellable refit normal assembly observed cancellation",
        )),
    }
}

/// The implicit-field → NURBS refit (radial pipeline; star-shaped domains).
///
/// # Errors
/// Invalid configuration, static allocation/work-cap refusal, radial bracket
/// failure, non-finite field/evaluation/report arithmetic, and degenerate
/// systems are returned as structured [`NurbsError`] values.
#[allow(clippy::too_many_lines)]
pub fn refit_radial(
    field: &dyn Fn([f64; 3]) -> f64,
    center: [f64; 3],
    r_max: f64,
    config: &RefitConfig,
) -> Result<Refit, NurbsError> {
    let (control_points, sample_points, probe) = validate_refit_request(center, r_max, config)?;
    let (nu, nv) = (config.nu, config.nv);
    let ku = open_uniform_knots(nu, config.degree)?;
    let kv = open_uniform_knots(nv, config.degree)?;
    // Construction already validated both sealed owners. Bind one admitted
    // view per axis and reuse it for every dense sample row instead of
    // rescanning the same immutable knots for every (u, v) pair.
    let admitted_ku = ku.admitted_after_validation();
    let admitted_kv = kv.admitted_after_validation();
    // Sample the field: radial projections on a (u, v) grid.
    let (mu, mv) = (config.samples_u, config.samples_v);
    let mut rows_b = try_vec_with_capacity(sample_points, "refit sample rows")?;
    let mut targets = try_vec_with_capacity(sample_points, "refit targets")?;
    let mut uvs = try_vec_with_capacity(sample_points, "refit parameters")?;
    for a in 0..mu {
        for b in 0..mv {
            #[allow(clippy::cast_precision_loss)]
            let (u, v) = ((a as f64 + 0.5) / mu as f64, (b as f64 + 0.5) / mv as f64);
            let dir = direction(u, v);
            let r = project_radial(field, center, dir, r_max)?;
            let target = [
                center[0] + r * dir[0],
                center[1] + r * dir[1],
                center[2] + r * dir[2],
            ];
            if target.iter().any(|coordinate| !coordinate.is_finite()) {
                return Err(refit_structure_error(
                    "projected refit target is not representable",
                ));
            }
            targets.push(target);
            let bu = basis_row(admitted_ku, u)?;
            let bv = basis_row(admitted_kv, v)?;
            let mut row = try_filled_vec(control_points, 0.0f64, "refit sample matrix row")?;
            for (i, &wu) in bu.iter().enumerate() {
                if wu == 0.0 {
                    continue;
                }
                for (j, &wv) in bv.iter().enumerate() {
                    if wv != 0.0 {
                        row[i * nv + j] = wu * wv;
                    }
                }
            }
            rows_b.push(row);
            uvs.push([u, v]);
        }
    }
    // Assemble and factor once, then solve the three coordinate right-hand
    // sides against the same deterministic factor.
    let mut net = try_filled_matrix(nu, nv, [0.0f64; 3], "refit control net")?;
    let factor = cholesky_factor(assemble_normal(&rows_b, nu, nv, config.lambda)?)?;
    for axis in 0..3 {
        let mut rhs = try_filled_vec(control_points, 0.0f64, "refit right-hand side")?;
        for (row, t) in rows_b.iter().zip(&targets) {
            for (k, &w) in row.iter().enumerate() {
                if w != 0.0 {
                    rhs[k] += w * t[axis];
                }
            }
        }
        let rhs = cholesky_solve_factored(&factor, rhs)?;
        for i in 0..nu {
            for j in 0..nv {
                net[i][j][axis] = rhs[i * nv + j];
            }
        }
    }
    // EXACT G0 seam closure: tie the u-boundary control columns.
    let (first_row, rest) = net.split_first_mut().expect("nu >= 2");
    let last_row = rest.last_mut().expect("nu >= 2");
    for (c0, c1) in first_row.iter_mut().zip(last_row.iter_mut()) {
        let avg = [
            f64::midpoint(c0[0], c1[0]),
            f64::midpoint(c0[1], c1[1]),
            f64::midpoint(c0[2], c1[2]),
        ];
        *c0 = avg;
        *c1 = avg;
    }
    let weights = try_filled_matrix(nu, nv, 1.0f64, "refit unit weights")?;
    let surface = NurbsSurface::new(ku, kv, &net, &weights)?;
    let report_surface = surface.admit()?;
    // ---- The honest report -------------------------------------------
    let mut rms = 0.0f64;
    let mut max_res = 0.0f64;
    let mut warnings = try_vec_with_capacity(sample_points, "refit warnings")?;
    for ((row, t), uv) in rows_b.iter().zip(&targets).zip(&uvs) {
        let mut p = [0.0f64; 3];
        for (k, &w) in row.iter().enumerate() {
            if w != 0.0 {
                let (i, j) = (k / nv, k % nv);
                for axis in 0..3 {
                    p[axis] += w * net[i][j][axis];
                }
            }
        }
        let r = norm3([p[0] - t[0], p[1] - t[1], p[2] - t[2]]);
        if !r.is_finite() {
            return Err(refit_structure_error(
                "refit residual arithmetic is non-finite",
            ));
        }
        rms += r * r;
        if !rms.is_finite() {
            return Err(refit_structure_error(
                "refit RMS accumulation is non-finite",
            ));
        }
        max_res = max_res.max(r);
        if r > config.warn_residual {
            warnings.push(LocalizedFitResidualWarning {
                uv: *uv,
                point: *t,
                residual: r,
            });
        }
    }
    #[allow(clippy::cast_precision_loss)]
    let rms_residual = det::sqrt(rms / sample_points as f64);
    // Spline → field: dense probe plus an analytic-model Lipschitz estimate;
    // the other reported direction stays the sampled fit-target worst case and
    // does not claim that a generic closure's targets belong to a source set.
    let mut sampled = 0.0f64;
    for a in 0..probe {
        for b in 0..probe {
            #[allow(clippy::cast_precision_loss)]
            let (u, v) = (
                (a as f64 + 0.5) / probe as f64,
                (b as f64 + 0.5) / probe as f64,
            );
            let p = report_surface.eval(u, v)?;
            let point = [p[0], p[1], p[2]];
            if point.iter().any(|coordinate| !coordinate.is_finite()) {
                return Err(refit_structure_error(format!(
                    "fitted surface returned a non-finite probe point at ({u}, {v})"
                )));
            }
            let field_value = field(point);
            if !field_value.is_finite() {
                return Err(refit_structure_error(format!(
                    "implicit field returned non-finite value at probe {point:?}"
                )));
            }
            sampled = sampled.max(field_value.abs());
        }
    }
    let coverage = max_res;
    let (lip_u, lip_v) = lipschitz_estimate(report_surface)?;
    let lip = lip_u + lip_v;
    #[allow(clippy::cast_precision_loss)]
    let probe_param_radius = 0.5 / probe as f64;
    let probe_spacing_estimate = (lip_u + lip_v) * probe_param_radius;
    if !rms_residual.is_finite()
        || !sampled.is_finite()
        || !lip.is_finite()
        || !probe_spacing_estimate.is_finite()
    {
        return Err(refit_structure_error(
            "refit report arithmetic is non-finite",
        ));
    }
    let (seam_g1, seam_g1_degenerate_samples) = seam_g1_diagnostic_admitted(report_surface)?;
    Ok(Refit {
        surface,
        report: RefitReport {
            rms_residual,
            max_residual: max_res,
            projected_target_to_spline_sampled: coverage,
            spline_to_field_sampled: sampled,
            spline_probe_spacing_estimate: probe_spacing_estimate,
            spline_lipschitz_estimate: lip,
            seam_g1_max: seam_g1,
            seam_g1_degenerate_samples,
            seam_g1_excludes_v_endpoints: true,
            warnings,
        },
    })
}

#[cfg(test)]
fn seam_g1_diagnostic(surface: &NurbsSurface<f64>) -> Result<(f64, usize), NurbsError> {
    seam_g1_diagnostic_admitted(surface.admit()?)
}

fn seam_g1_diagnostic_admitted(
    surface: AdmittedNurbsSurface<'_, f64>,
) -> Result<(f64, usize), NurbsError> {
    // Compare u-tangents across the exactly closed seam. Normalize each tangent
    // separately to avoid overflow/underflow in n0*n1, and clamp the rounded dot
    // product to the mathematical cosine range.
    let mut seam_g1 = 0.0f64;
    let mut degenerate = 0usize;
    // `v=0` and `v=1` are deliberately outside this diagnostic's typed scope:
    // radial/revolution fits commonly collapse them to poles where a tangent
    // direction is undefined. The report makes that exclusion machine-visible
    // pending a chart-aware pole audit.
    for b in 1..24 {
        let v = f64::from(b) / 24.0;
        let (_, du0, _) = surface.partials(0.0, v)?;
        let (_, du1, _) = surface.partials(1.0, v)?;
        let n0 = norm3(du0);
        let n1 = norm3(du1);
        if !n0.is_finite() || !n1.is_finite() {
            return Err(refit_structure_error(
                "refit seam-derivative arithmetic is non-finite",
            ));
        }
        if n0 == 0.0 || n1 == 0.0 {
            degenerate = degenerate.saturating_add(1);
            seam_g1 = f64::INFINITY;
            continue;
        }
        let unit0 = du0.map(|value| value / n0);
        let unit1 = du1.map(|value| value / n1);
        let cosang =
            (unit0[0] * unit1[0] + unit0[1] * unit1[1] + unit0[2] * unit1[2]).clamp(-1.0, 1.0);
        if !cosang.is_finite() {
            return Err(refit_structure_error(
                "refit seam-angle arithmetic is non-finite",
            ));
        }
        seam_g1 = seam_g1.max(1.0 - cosang);
    }
    Ok((seam_g1, degenerate))
}

#[cfg(test)]
mod tests {
    use super::*;
    use asupersync::types::Budget;
    use fs_exec::{CancelGate, ExecMode, StreamKey};
    use std::cell::Cell;

    fn with_refit_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        if cancelled {
            gate.request();
        }
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0xAEF1_7001,
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

    // G0/G5: the live transactional entry point preserves the fixed search.
    #[test]
    fn radial_projection_with_cx_matches_the_legacy_fixed_search() {
        let field = |point: [f64; 3]| point[0] - 1.0;
        let expected = project_radial(&field, [0.0; 3], [1.0, 0.0, 0.0], 2.0)
            .expect("legacy radial projection");
        with_refit_cx(false, |cx| {
            assert_eq!(
                project_radial_with_cx(&field, [0.0; 3], [1.0, 0.0, 0.0], 2.0, cx)
                    .expect("live cancellable projection"),
                RadialProjectionRun::Complete { radius: expected }
            );
        });
    }

    // G4: cancellation before work is observable and side-effect-free.
    #[test]
    fn radial_projection_pre_cancel_does_not_evaluate_the_field() {
        let calls = Cell::new(0usize);
        let field = |point: [f64; 3]| {
            calls.set(calls.get() + 1);
            point[0] - 1.0
        };
        with_refit_cx(true, |cx| {
            assert_eq!(
                project_radial_with_cx(&field, [0.0; 3], [1.0, 0.0, 0.0], 2.0, cx)
                    .expect("pre-cancellation is a terminal state"),
                RadialProjectionRun::Cancelled
            );
        });
        assert_eq!(calls.get(), 0);
    }

    // G4: constant-time and sampled refusals keep deterministic precedence.
    #[test]
    fn radial_projection_refusals_precede_or_dominate_cancellation() {
        let calls = Cell::new(0usize);
        let polls = Cell::new(0usize);
        let field = |_: [f64; 3]| {
            calls.set(calls.get() + 1);
            f64::NAN
        };
        let mut cancel_immediately = || {
            polls.set(polls.get() + 1);
            true
        };
        let malformed = project_radial_with_poll(
            &field,
            [f64::NAN, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            2.0,
            &mut cancel_immediately,
        )
        .expect_err("malformed constant-time input must refuse before polling");
        assert!(malformed.to_string().contains("center must be finite"));
        assert_eq!(polls.get(), 0);
        assert_eq!(calls.get(), 0);

        let mut cancel_after_field_call = || {
            polls.set(polls.get() + 1);
            calls.get() > 0
        };
        let nonfinite = project_radial_with_poll(
            &field,
            [0.0; 3],
            [1.0, 0.0, 0.0],
            2.0,
            &mut cancel_after_field_call,
        )
        .expect_err("a sampled non-finite field value must not become cancellation");
        assert!(nonfinite.to_string().contains("non-finite"));
        assert_eq!(polls.get(), 1);
        assert_eq!(calls.get(), 1);
    }

    // G4/G5: the legacy lower-bracket refusal remains a one-sample boundary.
    #[test]
    fn radial_projection_lower_bracket_refusal_short_circuits_the_upper_sample() {
        let calls = Cell::new(0usize);
        let polls = Cell::new(0usize);
        let field = |_: [f64; 3]| {
            calls.set(calls.get() + 1);
            1.0
        };
        let mut cancel_at_second_poll = || {
            polls.set(polls.get() + 1);
            polls.get() == 2
        };
        let error = project_radial_with_poll(
            &field,
            [0.0; 3],
            [1.0, 0.0, 0.0],
            2.0,
            &mut cancel_at_second_poll,
        )
        .expect_err("the failed lower bracket must refuse before the upper sample");
        assert!(error.to_string().contains("radial bracket failed"));
        assert_eq!(polls.get(), 1);
        assert_eq!(calls.get(), 1);
    }

    // G4/G5: a replayed checkpoint ordinal cancels at the same boundary.
    #[test]
    fn radial_projection_cancels_at_a_deterministic_bisection_boundary() {
        const CANCEL_AT_POLL: usize = 10;
        let calls = Cell::new(0usize);
        let polls = Cell::new(0usize);
        let field = |point: [f64; 3]| {
            calls.set(calls.get() + 1);
            point[0] - 1.0
        };
        let mut cancel_at_poll = || {
            polls.set(polls.get() + 1);
            polls.get() == CANCEL_AT_POLL
        };
        assert_eq!(
            project_radial_with_poll(&field, [0.0; 3], [1.0, 0.0, 0.0], 2.0, &mut cancel_at_poll,)
                .expect("mid-search cancellation"),
            RadialProjectionRun::Cancelled
        );
        assert_eq!(polls.get(), CANCEL_AT_POLL);
        assert_eq!(calls.get(), CANCEL_AT_POLL - 1);
    }

    // G4: completed local work is not published after final cancellation.
    #[test]
    fn radial_projection_has_a_final_publication_checkpoint() {
        let healthy_calls = Cell::new(0usize);
        let healthy_polls = Cell::new(0usize);
        let healthy_field = |point: [f64; 3]| {
            healthy_calls.set(healthy_calls.get() + 1);
            point[0] - 1.0
        };
        let mut count_without_cancelling = || {
            healthy_polls.set(healthy_polls.get() + 1);
            false
        };
        assert!(matches!(
            project_radial_with_poll(
                &healthy_field,
                [0.0; 3],
                [1.0, 0.0, 0.0],
                2.0,
                &mut count_without_cancelling,
            )
            .expect("healthy fixed search"),
            RadialProjectionRun::Complete { .. }
        ));
        assert_eq!(healthy_calls.get(), 42);
        assert_eq!(healthy_polls.get(), 43);

        let replay_calls = Cell::new(0usize);
        let replay_polls = Cell::new(0usize);
        let replay_field = |point: [f64; 3]| {
            replay_calls.set(replay_calls.get() + 1);
            point[0] - 1.0
        };
        let mut cancel_at_publication = || {
            replay_polls.set(replay_polls.get() + 1);
            replay_polls.get() == healthy_polls.get()
        };
        assert_eq!(
            project_radial_with_poll(
                &replay_field,
                [0.0; 3],
                [1.0, 0.0, 0.0],
                2.0,
                &mut cancel_at_publication,
            )
            .expect("publication cancellation"),
            RadialProjectionRun::Cancelled
        );
        assert_eq!(replay_calls.get(), healthy_calls.get());
        assert_eq!(replay_polls.get(), healthy_polls.get());
    }

    // G0/G5: the cancellable primitive preserves the deterministic thin-plate
    // matrix, including the four-neighbor regularizer's exact update order.
    #[test]
    fn normal_assembly_with_cx_matches_the_known_two_by_two_lattice() {
        let rows = vec![
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0, 0.0],
            vec![0.0, 0.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0],
        ];
        let expected = vec![
            vec![2.5, -1.0, -1.0, 0.5],
            vec![-1.0, 2.5, 0.5, -1.0],
            vec![-1.0, 0.5, 2.5, -1.0],
            vec![0.5, -1.0, -1.0, 2.5],
        ];
        with_refit_cx(false, |cx| {
            assert_eq!(
                assemble_refit_normal_with_cx(&rows, 2, 2, 0.25, cx)
                    .expect("cancellable normal assembly"),
                RefitNormalAssemblyRun::Complete {
                    matrix: expected.clone(),
                }
            );
        });
        assert_eq!(
            assemble_normal(&rows, 2, 2, 0.25).expect("legacy normal assembly"),
            expected
        );
    }

    // G4/G5: cancellation is replayable during validation/assembly and at the
    // final publication checkpoint; every cancelled path drops the matrix.
    #[test]
    fn normal_assembly_polling_is_bounded_and_transactional() {
        let rows = vec![vec![1.0; 8]; 16];
        let run = |cancel_at: Option<usize>| {
            let polls = Cell::new(0usize);
            let outcome = assemble_refit_normal_with_poll(&rows, 2, 4, 0.125, &mut || {
                polls.set(polls.get() + 1);
                cancel_at == Some(polls.get())
            })
            .expect("valid normal assembly");
            (outcome, polls.get())
        };

        let (complete, healthy_polls) = run(None);
        assert!(matches!(complete, RefitNormalAssemblyRun::Complete { .. }));
        assert!(healthy_polls > 4);
        assert_eq!(run(Some(1)), (RefitNormalAssemblyRun::Cancelled, 1));
        let middle_poll = healthy_polls / 2;
        assert_eq!(
            run(Some(middle_poll)),
            (RefitNormalAssemblyRun::Cancelled, middle_poll)
        );
        assert_eq!(
            run(Some(healthy_polls)),
            (RefitNormalAssemblyRun::Cancelled, healthy_polls)
        );
    }

    // G4: count-derived refusals precede cancellation, while traversed source
    // refusals remain inside the cancellable gate and allocate no output.
    #[test]
    fn normal_assembly_refusals_are_typed_and_preflighted() {
        let polls = Cell::new(0usize);
        let dimension_error = assemble_refit_normal_with_poll(&[], usize::MAX, 2, 0.0, &mut || {
            polls.set(polls.get() + 1);
            true
        })
        .expect_err("dimension overflow must precede cancellation");
        assert!(dimension_error.to_string().contains("dimension overflow"));
        assert_eq!(polls.get(), 0);

        let work_error = preflight_refit_normal_assembly(usize::MAX, 1, 1, 0.0)
            .expect_err("unbounded borrowed-row work must be refused");
        assert!(work_error.to_string().contains("work estimate"));

        let malformed = assemble_refit_normal_with_poll(&[vec![1.0; 3]], 2, 2, 0.0, &mut || {
            polls.set(polls.get() + 1);
            false
        })
        .expect_err("malformed borrowed row must be refused");
        assert!(malformed.to_string().contains("length 3, expected 4"));
        assert_eq!(polls.get(), 1);

        with_refit_cx(true, |cx| {
            assert_eq!(
                assemble_refit_normal_with_cx(&[vec![1.0; 4]], 2, 2, 0.0, cx)
                    .expect("pre-cancellation is a terminal state"),
                RefitNormalAssemblyRun::Cancelled
            );
        });
    }

    // G0/G5: the transactional factor preserves the legacy deterministic
    // lower-triangle Cholesky arithmetic and leaves the unused upper triangle.
    #[test]
    fn normal_factor_with_cx_matches_the_known_spd_factor() {
        let matrix = vec![vec![4.0, 2.0], vec![2.0, 5.0]];
        let expected = vec![vec![2.0, 2.0], vec![1.0, 2.0]];
        with_refit_cx(false, |cx| {
            assert_eq!(
                factor_refit_normal_with_cx(matrix.clone(), cx)
                    .expect("cancellable normal factorization"),
                RefitNormalFactorRun::Complete {
                    factor: expected.clone(),
                }
            );
        });
        assert_eq!(
            cholesky_factor(matrix).expect("legacy normal factorization"),
            expected
        );
    }

    // G4/G5: the same checkpoint ordinal cancels factorization during source
    // validation, arithmetic, or final publication without exposing mutation.
    #[test]
    fn normal_factor_polling_is_bounded_and_transactional() {
        let mut matrix = vec![vec![0.0; 16]; 16];
        for (index, row) in matrix.iter_mut().enumerate() {
            row[index] = 1.0;
        }
        let run = |cancel_at: Option<usize>| {
            let polls = Cell::new(0usize);
            let outcome = factor_refit_normal_with_poll(matrix.clone(), &mut || {
                polls.set(polls.get() + 1);
                cancel_at == Some(polls.get())
            })
            .expect("valid normal factorization");
            (outcome, polls.get())
        };

        let (complete, healthy_polls) = run(None);
        assert!(matches!(complete, RefitNormalFactorRun::Complete { .. }));
        assert!(healthy_polls > 4);
        assert_eq!(run(Some(1)), (RefitNormalFactorRun::Cancelled, 1));
        let middle_poll = healthy_polls / 2;
        assert_eq!(
            run(Some(middle_poll)),
            (RefitNormalFactorRun::Cancelled, middle_poll)
        );
        assert_eq!(
            run(Some(healthy_polls)),
            (RefitNormalFactorRun::Cancelled, healthy_polls)
        );
    }

    // G4: constant/count refusals beat cancellation, and traversed structural
    // refusals remain typed inside the gate.
    #[test]
    fn normal_factor_refusals_are_typed_and_preflighted() {
        let work_error = preflight_refit_normal_factor(1_000)
            .expect_err("cubic work above the process cap must be refused");
        assert!(work_error.to_string().contains("work estimate"));

        let polls = Cell::new(0usize);
        let empty = factor_refit_normal_with_poll(Vec::new(), &mut || {
            polls.set(polls.get() + 1);
            true
        })
        .expect_err("empty input must refuse before cancellation");
        assert!(empty.to_string().contains("nonempty square matrix"));
        assert_eq!(polls.get(), 0);

        let malformed = factor_refit_normal_with_poll(vec![vec![1.0], vec![0.0, 1.0]], &mut || {
            polls.set(polls.get() + 1);
            false
        })
        .expect_err("non-square input must be refused");
        assert!(malformed.to_string().contains("length 1, expected 2"));
        assert_eq!(polls.get(), 1);

        let asymmetric =
            factor_refit_normal_with_poll(vec![vec![1.0, 0.0], vec![1.0, 1.0]], &mut || false)
                .expect_err("asymmetric input must not be factored as SPD");
        assert!(asymmetric.to_string().contains("asymmetric at (1, 0)"));
    }

    // G0/G5: the transactional solve preserves the legacy forward/back
    // substitution order for a factor with an exactly known solution.
    #[test]
    fn normal_solve_with_cx_matches_the_known_solution() {
        let factor = vec![vec![2.0, 2.0], vec![1.0, 2.0]];
        let rhs = vec![8.0, 12.0];
        let expected = vec![1.0, 2.0];
        with_refit_cx(false, |cx| {
            assert_eq!(
                solve_refit_normal_with_cx(&factor, rhs.clone(), cx)
                    .expect("cancellable normal solve"),
                RefitNormalSolveRun::Complete {
                    solution: expected.clone(),
                }
            );
        });
        assert_eq!(
            cholesky_solve_factored(&factor, rhs).expect("legacy normal solve"),
            expected
        );
    }

    // G4/G5: validation, both substitutions, and solution publication share
    // replayable checkpoints, and cancellation never exposes a partial RHS.
    #[test]
    fn normal_solve_polling_is_bounded_and_transactional() {
        let mut factor = vec![vec![0.0; 32]; 32];
        for (index, row) in factor.iter_mut().enumerate() {
            row[index] = 1.0;
        }
        let rhs = vec![1.0; 32];
        let run = |cancel_at: Option<usize>| {
            let polls = Cell::new(0usize);
            let outcome = solve_refit_normal_with_poll(&factor, rhs.clone(), &mut || {
                polls.set(polls.get() + 1);
                cancel_at == Some(polls.get())
            })
            .expect("valid normal solve");
            (outcome, polls.get())
        };

        let (complete, healthy_polls) = run(None);
        assert!(matches!(complete, RefitNormalSolveRun::Complete { .. }));
        assert!(healthy_polls > 4);
        assert_eq!(run(Some(1)), (RefitNormalSolveRun::Cancelled, 1));
        let middle_poll = healthy_polls / 2;
        assert_eq!(
            run(Some(middle_poll)),
            (RefitNormalSolveRun::Cancelled, middle_poll)
        );
        assert_eq!(
            run(Some(healthy_polls)),
            (RefitNormalSolveRun::Cancelled, healthy_polls)
        );
    }

    // G4: constant/count refusals precede cancellation, while traversed factor
    // and RHS failures remain typed inside the gate.
    #[test]
    fn normal_solve_refusals_are_typed_and_preflighted() {
        let work_error = preflight_refit_normal_solve(20_000, 20_000)
            .expect_err("quadratic solve work above the cap must be refused");
        assert!(work_error.to_string().contains("work estimate"));

        let polls = Cell::new(0usize);
        let dimension_error = solve_refit_normal_with_poll(&[], vec![1.0], &mut || {
            polls.set(polls.get() + 1);
            true
        })
        .expect_err("dimension mismatch must precede cancellation");
        assert!(dimension_error.to_string().contains("0 rows"));
        assert_eq!(polls.get(), 0);

        let malformed =
            solve_refit_normal_with_poll(&[vec![1.0], vec![0.0, 1.0]], vec![1.0, 1.0], &mut || {
                polls.set(polls.get() + 1);
                false
            })
            .expect_err("non-square factor must be refused");
        assert!(malformed.to_string().contains("length 1, expected 2"));
        assert_eq!(polls.get(), 1);

        let nonpositive = solve_refit_normal_with_poll(&[vec![0.0]], vec![1.0], &mut || false)
            .expect_err("non-positive factor diagonal must be refused");
        assert!(nonpositive.to_string().contains("non-positive diagonal"));

        let nonfinite_rhs =
            solve_refit_normal_with_poll(&[vec![1.0]], vec![f64::NAN], &mut || false)
                .expect_err("non-finite RHS must be refused");
        assert!(nonfinite_rhs.to_string().contains("right-hand side"));
    }

    #[test]
    fn refit_admission_refuses_invalid_or_unbounded_work_before_field_evaluation() {
        let calls = Cell::new(0usize);
        let field = |point: [f64; 3]| {
            calls.set(calls.get() + 1);
            point[0]
        };
        let zero_samples = RefitConfig {
            samples_u: 0,
            ..RefitConfig::default()
        };
        assert!(refit_radial(&field, [0.0; 3], 1.0, &zero_samples).is_err());

        let overflowing_grid = RefitConfig {
            nu: usize::MAX,
            nv: 2,
            degree: 1,
            ..RefitConfig::default()
        };
        assert!(refit_radial(&field, [0.0; 3], 1.0, &overflowing_grid).is_err());

        let excessive_grid = RefitConfig {
            nu: 4096,
            nv: 4096,
            degree: 1,
            ..RefitConfig::default()
        };
        assert!(refit_radial(&field, [0.0; 3], 1.0, &excessive_grid).is_err());

        let excessive_work = RefitConfig {
            nu: 32,
            nv: 32,
            degree: 3,
            ..RefitConfig::default()
        };
        let error = refit_radial(&field, [0.0; 3], 1.0, &excessive_work)
            .expect_err("work above the process cap must refuse before sampling");
        assert!(
            error.to_string().contains("work estimate"),
            "the reachable work-only refusal must remain distinct: {error}"
        );
        assert_eq!(
            calls.get(),
            0,
            "all shape/cap refusals must precede field evaluation"
        );
    }

    #[test]
    fn refit_work_accounting_is_exact_and_overflow_fallible() {
        assert_eq!(
            checked_refit_work_product(&[7, 11, 13], "test product").expect("small product"),
            1_001
        );
        assert_eq!(
            checked_refit_work_sum(&[7, 11, 13], "test sum").expect("small sum"),
            31
        );
        assert!(matches!(
            checked_refit_work_product(&[u128::MAX, 2], "test product"),
            Err(NurbsError::Structure { ref what })
                if what == "refit test product work accounting overflows u128"
        ));
        assert!(matches!(
            checked_refit_work_sum(&[u128::MAX, 1], "test sum"),
            Err(NurbsError::Structure { ref what })
                if what == "refit test sum work accounting overflows u128"
        ));
    }

    #[test]
    fn admitted_refit_basis_rows_match_owning_basis_without_rescans() {
        let knots = open_uniform_knots(4, 2).expect("quadratic refit knots");
        let parameter = 0.375;
        let (span, values) = knots.basis(parameter).expect("owning basis oracle");
        let mut expected = vec![0.0; knots.control_count()];
        for (offset, value) in values.into_iter().enumerate() {
            expected[span - knots.degree() + offset] = value;
        }
        assert_eq!(
            basis_row(knots.admitted_after_validation(), parameter)
                .expect("admitted dense basis row"),
            expected
        );

        let owning_error = knots
            .basis(-0.25)
            .expect_err("owning out-of-domain refusal");
        let admitted_error = basis_row(knots.admitted_after_validation(), -0.25)
            .expect_err("admitted out-of-domain refusal");
        assert_eq!(admitted_error, owning_error);
    }

    #[test]
    fn refit_refuses_nonfinite_field_samples() {
        let config = RefitConfig {
            nu: 2,
            nv: 2,
            degree: 1,
            samples_u: 2,
            samples_v: 2,
            probe: 2,
            ..RefitConfig::default()
        };
        let error = refit_radial(&|_| f64::NAN, [0.0; 3], 1.0, &config)
            .expect_err("non-finite fields must refuse");
        assert!(error.to_string().contains("non-finite"));
    }

    #[test]
    fn refit_owned_buffers_are_fallible_and_shape_preserving() {
        let matrix = try_filled_matrix(2, 3, 7.0f64, "test matrix").expect("small matrix");
        assert_eq!(matrix, vec![vec![7.0; 3]; 2]);

        let error = try_filled_vec(usize::MAX, 0.0f64, "test vector")
            .expect_err("unrepresentable capacity must be a typed refusal");
        assert!(
            matches!(error, NurbsError::Domain { ref what } if what == "test vector allocation was refused")
        );

        let error = try_filled_matrix(1, usize::MAX, 0.0f64, "test matrix")
            .expect_err("unrepresentable inner capacity must be a typed refusal");
        assert!(
            matches!(error, NurbsError::Domain { ref what } if what == "test matrix allocation was refused")
        );
    }

    #[test]
    fn lipschitz_estimate_uses_the_collapsed_clamp_span() {
        // Regression: `lipschitz_estimate` must divide by the ACTUAL hodograph knot
        // span, not the closed form (n−p). On a clamped open-uniform knot vector
        // the END span collapses from p/(n−p) (interior) to 1/(n−p), so a control
        // difference at the clamp has a true hodograph coefficient p·(n−p), not
        // (n−p). Under-bounding there makes the analytic estimate too tight.
        let (n, p) = (8usize, 3usize);
        let ku = open_uniform_knots(n, p).expect("u knots");
        let kv = open_uniform_knots(2, 1).expect("v knots"); // linear in v
        // Large jump ONLY between the first two u-rows (the clamped end); every
        // other u-difference is zero, so max‖ΔC_u‖ lives in the collapsed span.
        let jump = 2.0;
        let net: Vec<Vec<[f64; 3]>> = (0..n)
            .map(|i| {
                let x = if i == 0 { 0.0 } else { jump };
                vec![[x, 0.0, 0.0], [x, 1.0, 0.0]]
            })
            .collect();
        let weights = vec![vec![1.0, 1.0]; n];
        let surface = NurbsSurface::new(ku, kv, &net, &weights).expect("surface");
        let admitted = surface.admit().expect("admitted surface");
        let (lu, lv) = lipschitz_estimate(admitted).expect("measured Lipschitz estimate");
        with_refit_cx(false, |cx| {
            assert_eq!(
                estimate_refit_surface_lipschitz_with_cx(admitted, cx)
                    .expect("cancellable measured estimate"),
                RefitLipschitzEstimateRun::Complete {
                    u_estimate: lu,
                    v_estimate: lv,
                }
            );
        });
        // Analytic formula: p · jump / (1/(n−p)) = p·(n−p)·jump.
        let expected = p as f64 * (n - p) as f64 * jump;
        let closed_form = (n - p) as f64 * jump; // the old factor-p under-estimate
        assert!(
            (lu - expected).abs() < 1e-9,
            "L_u must use the collapsed clamp span: got {lu}, expected {expected}, \
             old closed-form under-estimate {closed_form}"
        );
        assert!(
            lu > closed_form + 1e-9,
            "the per-span estimate must exceed the closed-form under-estimate ({lu} vs {closed_form})"
        );
    }

    // G4/G5: the control traversal and tuple publication share replayable
    // checkpoints, and cancellation never publishes one directional estimate.
    #[test]
    fn lipschitz_estimate_polling_is_bounded_and_transactional() {
        let knots = open_uniform_knots(9, 1).expect("linear knots");
        let net: Vec<Vec<[f64; 3]>> = (0_u32..9)
            .map(|i| {
                (0_u32..9)
                    .map(|j| [f64::from(i), f64::from(j), f64::from(i + j) * 0.125])
                    .collect()
            })
            .collect();
        let weights = vec![vec![1.0; 9]; 9];
        let surface = NurbsSurface::new(knots.clone(), knots, &net, &weights).expect("surface");
        let admitted = surface.admit().expect("admitted surface");
        let run = |cancel_at: Option<usize>| {
            let polls = Cell::new(0usize);
            let outcome = estimate_refit_surface_lipschitz_with_poll(admitted, &mut || {
                polls.set(polls.get() + 1);
                cancel_at == Some(polls.get())
            })
            .expect("finite measured estimate");
            (outcome, polls.get())
        };

        let (complete, healthy_polls) = run(None);
        assert!(matches!(
            complete,
            RefitLipschitzEstimateRun::Complete { .. }
        ));
        assert_eq!(healthy_polls, 3, "entry, 64-control stride, publication");
        assert_eq!(run(Some(1)), (RefitLipschitzEstimateRun::Cancelled, 1));
        assert_eq!(run(Some(2)), (RefitLipschitzEstimateRun::Cancelled, 2));
        assert_eq!(run(Some(3)), (RefitLipschitzEstimateRun::Cancelled, 3));
        with_refit_cx(true, |cx| {
            assert_eq!(
                estimate_refit_surface_lipschitz_with_cx(admitted, cx)
                    .expect("pre-cancellation is terminal"),
                RefitLipschitzEstimateRun::Cancelled
            );
        });
    }

    // G0/G4: hostile counts refuse before traversal, and finite admitted
    // controls cannot silently publish an infinite measured coefficient.
    #[test]
    fn lipschitz_estimate_refuses_unbounded_or_nonfinite_arithmetic() {
        let work_error = preflight_refit_lipschitz_estimate(usize::MAX, 2)
            .expect_err("unbounded traversal work must be refused");
        assert!(work_error.to_string().contains("work estimate"));

        let knots = open_uniform_knots(2, 1).expect("linear knots");
        let net = vec![
            vec![[-f64::MAX, 0.0, 0.0]; 2],
            vec![[f64::MAX, 0.0, 0.0]; 2],
        ];
        let weights = vec![vec![1.0; 2]; 2];
        let surface = NurbsSurface::new(knots.clone(), knots, &net, &weights)
            .expect("finite extreme surface");
        let error = estimate_refit_surface_lipschitz_with_poll(
            surface.admit().expect("admitted extreme surface"),
            &mut || false,
        )
        .expect_err("overflowed coefficient must not publish");
        assert!(
            error
                .to_string()
                .contains("Lipschitz estimate became non-finite")
        );
    }

    #[test]
    fn degenerate_seam_tangents_are_explicitly_no_claim() {
        let line = open_uniform_knots(2, 1).expect("linear knots");
        let points = vec![vec![[0.0; 3]; 2]; 2];
        let weights = vec![vec![1.0; 2]; 2];
        let surface =
            NurbsSurface::new(line.clone(), line, &points, &weights).expect("degenerate surface");
        let (g1, degenerate) = seam_g1_diagnostic(&surface).expect("bounded diagnostic");
        let admitted = surface.admit().expect("admitted surface");
        let admitted_result =
            seam_g1_diagnostic_admitted(admitted).expect("admitted bounded diagnostic");
        assert_eq!(admitted_result, (g1, degenerate));
        assert!(g1.is_infinite(), "undefined tangent direction is no-claim");
        assert_eq!(degenerate, 23, "every retained seam sample is named");
    }
}
