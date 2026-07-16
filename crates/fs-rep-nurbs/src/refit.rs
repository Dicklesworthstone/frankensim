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
use crate::basis::{AdmittedKnotVector, BASIS_MAX_WORK_UNITS, BasisRun, KnotVector};
use crate::closest::norm3;
use crate::surface::{
    AdmittedNurbsSurface, NurbsSurface, SurfaceEvaluationRun, SurfacePartialsRun,
};
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
const REFIT_SEAM_G1_SAMPLE_COUNT: usize = 23;

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
    add_bytes(bytes_for(control_points, size_of::<[f64; 3]>())?)?; // x/y/z RHS vectors
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

/// Transactional outcome of the fixed paired-side G1 seam diagnostic.
///
/// Positive infinity is a complete measured no-claim result when one or more
/// retained tangent pairs are degenerate; cancellation publishes no tuple.
#[must_use]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RefitSeamG1DiagnosticRun {
    /// All 23 interior parameter pairs were measured.
    Complete {
        /// Largest retained `1 - cos(angle)` deviation, or positive infinity
        /// when a sampled tangent direction is undefined.
        max_deviation: f64,
        /// Number of retained parameter pairs with a zero tangent direction.
        degenerate_samples: usize,
    },
    /// Cancellation was observed before tuple publication.
    Cancelled,
}

/// Transactional outcome of cancellation-aware dense basis-row expansion for
/// the refit sample matrix.
#[must_use]
#[derive(Debug, Clone, PartialEq)]
pub enum RefitDenseBasisRun {
    /// The complete control-axis row is safe to publish.
    Complete {
        /// One value per control, with inactive basis entries set to zero.
        values: Vec<f64>,
    },
    /// Cancellation was observed before publication.
    Cancelled,
}

/// Transactional outcome of cancellation-aware tensor sample-basis assembly.
#[must_use]
#[derive(Debug, Clone, PartialEq)]
pub enum RefitSampleBasisRowRun {
    /// Both dense axis rows and their complete tensor product were assembled.
    Complete {
        /// Row-major U-by-V basis values for one refit sample.
        values: Vec<f64>,
    },
    /// Cancellation was observed before publication.
    Cancelled,
}

/// Transactional outcome of cancellation-aware refit right-hand-side assembly.
#[must_use]
#[derive(Debug, Clone, PartialEq)]
pub enum RefitRightHandSidesRun {
    /// Borrowed inputs were validated and all three RHS vectors were accumulated.
    Complete {
        /// X/Y/Z right-hand sides, each with one entry per control.
        right_hand_sides: [Vec<f64>; 3],
    },
    /// Cancellation was observed before publication.
    Cancelled,
}

/// Transactional outcome of cancellation-aware paired-parameter residual
/// summarization.
///
/// The returned values describe only the retained sample pairs supplied to the
/// primitive. They do not certify continuum fit quality or geometric distance.
#[must_use]
#[derive(Debug, Clone, PartialEq)]
pub enum RefitResidualSummaryRun {
    /// Every borrowed input was validated and the complete summary is safe to
    /// publish.
    Complete {
        /// Root-mean-square residual over the retained paired samples.
        rms_residual: f64,
        /// Largest residual over the retained paired samples.
        max_residual: f64,
        /// Samples whose residual strictly exceeded the requested threshold.
        warnings: Vec<LocalizedFitResidualWarning>,
    },
    /// Cancellation was observed before publication.
    Cancelled,
}

/// Transactional outcome of a dense measured refit field probe.
///
/// The scalar is only the largest retained `abs(field(surface(u, v)))` sample.
/// It carries no field, zero-set, distance, or between-sample authority.
#[must_use]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RefitFieldProbeRun {
    /// Every grid point and opaque field call completed with finite arithmetic.
    Complete {
        /// Maximum absolute field value over the retained cell-center grid.
        max_abs_field: f64,
    },
    /// Cancellation was observed before scalar publication.
    Cancelled,
}

/// Transactional outcome of cancellation-aware open-uniform refit knot
/// construction.
#[must_use]
#[derive(Debug, PartialEq)]
pub enum RefitOpenUniformKnotRun {
    /// Filling and sealed structural validation completed.
    Complete {
        /// Complete checked open-uniform knot owner.
        knots: KnotVector<f64>,
    },
    /// Cancellation was observed before publication.
    Cancelled,
}

/// Transactional outcome of exact periodic-u refit seam closure.
#[must_use]
#[derive(Debug, PartialEq)]
pub enum RefitSeamClosureRun {
    /// The complete rectangular net has finite controls and tied boundary rows.
    Complete {
        /// Exclusively owned control net after exact G0 closure.
        net: Vec<Vec<[f64; 3]>>,
    },
    /// Cancellation was observed; the consumed net was dropped unpublished.
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RefitOpenUniformKnotPlan {
    order: usize,
    inner: usize,
    knot_count: usize,
}

fn preflight_open_uniform_refit_knots(
    n: usize,
    degree: usize,
) -> Result<RefitOpenUniformKnotPlan, NurbsError> {
    if degree == 0 {
        return Err(refit_structure_error(
            "refit open-uniform knot degree must be positive",
        ));
    }
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
    let validation_work = KnotVector::<f64>::validation_work_for(knot_count, degree)?;
    let build_work =
        checked_refit_work_sum(&[knot_count as u128, 2], "open-uniform knot construction")?;
    let total_work = checked_refit_work_sum(
        &[build_work, validation_work],
        "open-uniform knot aggregate",
    )?;
    KnotVector::<f64>::enforce_work(total_work, "refit open-uniform knot construction")?;
    let retained_bytes = knot_count
        .checked_mul(size_of::<f64>())
        .ok_or_else(|| refit_structure_error("refit knot byte estimate overflow"))?;
    if retained_bytes > REFIT_MAX_ALLOC_BYTES {
        return Err(refit_structure_error(format!(
            "refit knot vector retains {retained_bytes} requested bytes above static cap {REFIT_MAX_ALLOC_BYTES}"
        )));
    }
    Ok(RefitOpenUniformKnotPlan {
        order,
        inner,
        knot_count,
    })
}

/// Construct one checked open-uniform refit knot vector with bounded
/// cancellation polling.
///
/// Aggregate fill-plus-validation work and retained knot payload are admitted
/// before allocation. One gate then spans fallible reservation, ordered leading
/// clamp/interior/trailing-clamp fill, complete sealed knot validation, and
/// final owner publication. Cancellation drops the partial or fully filled but
/// unpublished knot storage. The primitive does not consume the `Cx` budget or
/// make the full refit pipeline cancellation-aware.
///
/// # Errors
/// Returns a structured dimension, checked work/retained-byte, allocation,
/// finite-arithmetic, or sealed knot-validation refusal.
pub fn build_open_uniform_refit_knots_with_cx(
    control_count: usize,
    degree: usize,
    cx: &Cx<'_>,
) -> Result<RefitOpenUniformKnotRun, NurbsError> {
    let mut should_cancel = || cx.checkpoint().is_err();
    build_open_uniform_refit_knots_with_poll(control_count, degree, &mut should_cancel)
}

fn build_open_uniform_refit_knots_with_poll(
    control_count: usize,
    degree: usize,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<RefitOpenUniformKnotRun, NurbsError> {
    let plan = preflight_open_uniform_refit_knots(control_count, degree)?;
    if should_cancel() {
        return Ok(RefitOpenUniformKnotRun::Cancelled);
    }
    let mut knots = try_vec_with_capacity(plan.knot_count, "refit knot vector")?;
    if should_cancel() {
        return Ok(RefitOpenUniformKnotRun::Cancelled);
    }
    let mut operations_since_poll = 0usize;
    for _ in 0..plan.order {
        knots.push(0.0);
        if refit_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(RefitOpenUniformKnotRun::Cancelled);
        }
    }
    #[allow(clippy::cast_precision_loss)]
    for k in 1..plan.inner {
        let value = k as f64 / plan.inner as f64;
        if !value.is_finite() {
            return Err(refit_structure_error(
                "refit open-uniform knot arithmetic became non-finite",
            ));
        }
        knots.push(value);
        if refit_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(RefitOpenUniformKnotRun::Cancelled);
        }
    }
    for _ in 0..plan.order {
        knots.push(1.0);
        if refit_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(RefitOpenUniformKnotRun::Cancelled);
        }
    }
    debug_assert_eq!(knots.len(), plan.knot_count);
    match KnotVector::new_with_poll(knots, degree, should_cancel)? {
        Some(knots) => Ok(RefitOpenUniformKnotRun::Complete { knots }),
        None => Ok(RefitOpenUniformKnotRun::Cancelled),
    }
}

fn open_uniform_knots(n: usize, degree: usize) -> Result<KnotVector<f64>, NurbsError> {
    let mut never_cancel = || false;
    match build_open_uniform_refit_knots_with_poll(n, degree, &mut never_cancel)? {
        RefitOpenUniformKnotRun::Complete { knots } => Ok(knots),
        RefitOpenUniformKnotRun::Cancelled => Err(refit_structure_error(
            "non-cancellable open-uniform refit knot construction observed cancellation",
        )),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RefitSeamClosurePlan {
    rows: usize,
    cols: usize,
    control_count: usize,
}

fn preflight_refit_u_seam_closure(
    rows: usize,
    cols: usize,
) -> Result<RefitSeamClosurePlan, NurbsError> {
    if rows < 2 {
        return Err(refit_structure_error(
            "refit U-seam closure requires at least two control rows",
        ));
    }
    if cols == 0 {
        return Err(refit_structure_error(
            "refit U-seam closure requires a nonempty control row",
        ));
    }
    let control_count = rows
        .checked_mul(cols)
        .ok_or_else(|| refit_structure_error("refit U-seam control-grid size overflow"))?;
    // Charge rectangular-shape checks, three finite coordinates plus one
    // checkpoint unit per control, conservative midpoint arithmetic, entry
    // gates, and final publication.
    let validation_work =
        checked_refit_work_product(&[control_count as u128, 4], "U-seam control validation")?;
    let closure_work = checked_refit_work_product(&[cols as u128, 16], "U-seam midpoint closure")?;
    let total_work = checked_refit_work_sum(
        &[rows as u128, validation_work, closure_work, 3],
        "U-seam closure aggregate",
    )?;
    if total_work > REFIT_MAX_WORK_UNITS {
        return Err(refit_structure_error(format!(
            "refit U-seam closure requests {total_work} work units above static cap {REFIT_MAX_WORK_UNITS}"
        )));
    }

    let row_header_bytes = (rows as u128)
        .checked_mul(size_of::<Vec<[f64; 3]>>() as u128)
        .ok_or_else(|| refit_structure_error("refit U-seam row-byte accounting overflow"))?;
    let control_bytes = (control_count as u128)
        .checked_mul(size_of::<[f64; 3]>() as u128)
        .ok_or_else(|| refit_structure_error("refit U-seam control-byte accounting overflow"))?;
    let retained_bytes = row_header_bytes
        .checked_add(control_bytes)
        .ok_or_else(|| refit_structure_error("refit U-seam retained-byte accounting overflow"))?;
    if retained_bytes > REFIT_MAX_ALLOC_BYTES as u128 {
        return Err(refit_structure_error(format!(
            "refit U-seam net retains {retained_bytes} requested bytes above static cap {REFIT_MAX_ALLOC_BYTES}"
        )));
    }
    Ok(RefitSeamClosurePlan {
        rows,
        cols,
        control_count,
    })
}

/// Tie the first and last U rows of one owned refit control net exactly.
///
/// Count-derived aggregate work and the valid rectangular net's requested
/// payload are admitted before the first checkpoint or scan. One gate then
/// spans row-shape validation, finite-control validation, deterministic midpoint
/// tying, and final publication, polling after at most 64 rows or controls.
/// Taking exclusive ownership makes cancellation transactional: a partially
/// tied net is dropped and never returned. The primitive allocates no numerical
/// payload, does not consume the `Cx` budget, and grants no G1 or surface-validity
/// authority.
///
/// # Errors
/// Returns a structured shape, checked work/retained-byte, or non-finite-control
/// refusal.
pub fn close_refit_u_seam_with_cx(
    net: Vec<Vec<[f64; 3]>>,
    cx: &Cx<'_>,
) -> Result<RefitSeamClosureRun, NurbsError> {
    let mut should_cancel = || cx.checkpoint().is_err();
    close_refit_u_seam_with_poll(net, &mut should_cancel)
}

fn close_refit_u_seam_with_poll(
    mut net: Vec<Vec<[f64; 3]>>,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<RefitSeamClosureRun, NurbsError> {
    let rows = net.len();
    let cols = net.first().map_or(0, Vec::len);
    let plan = preflight_refit_u_seam_closure(rows, cols)?;
    if should_cancel() {
        return Ok(RefitSeamClosureRun::Cancelled);
    }

    let mut operations_since_poll = 0usize;
    for (row_index, row) in net.iter().enumerate() {
        if row.len() != plan.cols {
            return Err(refit_structure_error(format!(
                "refit U-seam row {row_index} has length {}, expected {}",
                row.len(),
                plan.cols
            )));
        }
        if refit_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(RefitSeamClosureRun::Cancelled);
        }
    }
    if should_cancel() {
        return Ok(RefitSeamClosureRun::Cancelled);
    }

    operations_since_poll = 0;
    for row in &net {
        for control in row {
            if control.iter().any(|coordinate| !coordinate.is_finite()) {
                return Err(refit_structure_error(
                    "refit U-seam control coordinates must be finite",
                ));
            }
            if refit_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(RefitSeamClosureRun::Cancelled);
            }
        }
    }
    if should_cancel() {
        return Ok(RefitSeamClosureRun::Cancelled);
    }

    let (first_row, rest) = net
        .split_first_mut()
        .ok_or_else(|| refit_structure_error("refit U-seam control net is empty"))?;
    let last_row = rest
        .last_mut()
        .ok_or_else(|| refit_structure_error("refit U-seam control net has only one row"))?;
    operations_since_poll = 0;
    for (first, last) in first_row.iter_mut().zip(last_row.iter_mut()) {
        let midpoint = [
            f64::midpoint(first[0], last[0]),
            f64::midpoint(first[1], last[1]),
            f64::midpoint(first[2], last[2]),
        ];
        if midpoint.iter().any(|coordinate| !coordinate.is_finite()) {
            return Err(refit_structure_error(
                "refit U-seam midpoint arithmetic became non-finite",
            ));
        }
        *first = midpoint;
        *last = midpoint;
        if refit_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(RefitSeamClosureRun::Cancelled);
        }
    }
    debug_assert_eq!(net.len(), plan.rows);
    debug_assert_eq!(net.iter().map(Vec::len).sum::<usize>(), plan.control_count);
    if should_cancel() {
        return Ok(RefitSeamClosureRun::Cancelled);
    }
    Ok(RefitSeamClosureRun::Complete { net })
}

fn close_refit_u_seam(net: Vec<Vec<[f64; 3]>>) -> Result<Vec<Vec<[f64; 3]>>, NurbsError> {
    let mut never_cancel = || false;
    match close_refit_u_seam_with_poll(net, &mut never_cancel)? {
        RefitSeamClosureRun::Complete { net } => Ok(net),
        RefitSeamClosureRun::Cancelled => Err(refit_structure_error(
            "non-cancellable refit U-seam closure observed cancellation",
        )),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RefitDenseBasisPlan {
    total_work: u128,
    peak_scalars: usize,
}

fn preflight_refit_dense_basis(
    control_count: usize,
    degree: usize,
) -> Result<RefitDenseBasisPlan, NurbsError> {
    if control_count == 0 {
        return Err(refit_structure_error(
            "refit dense basis requires at least one control",
        ));
    }
    let order = degree
        .checked_add(1)
        .ok_or_else(|| refit_structure_error("refit dense basis order overflows usize"))?;
    if order > control_count {
        return Err(refit_structure_error(
            "refit dense basis order exceeds the control count",
        ));
    }
    let triangle_work =
        checked_refit_work_product(&[degree as u128, order as u128], "dense basis triangle")? / 2;
    let basis_work = checked_refit_work_sum(
        &[triangle_work, order as u128, control_count as u128],
        "dense basis evaluation",
    )?;
    let dense_work = checked_refit_work_sum(
        &[control_count as u128, order as u128, 2],
        "dense basis expansion",
    )?;
    let total_work = checked_refit_work_sum(&[basis_work, dense_work], "dense basis aggregate")?;
    if total_work > BASIS_MAX_WORK_UNITS {
        return Err(refit_structure_error(format!(
            "refit dense basis requests {total_work} work units above defensive ceiling {BASIS_MAX_WORK_UNITS}"
        )));
    }

    let basis_scratch_scalars = order
        .checked_mul(3)
        .ok_or_else(|| refit_structure_error("refit dense basis scratch size overflows usize"))?;
    let expansion_scalars = control_count
        .checked_add(order)
        .ok_or_else(|| refit_structure_error("refit dense basis live size overflows usize"))?;
    let peak_scalars = basis_scratch_scalars.max(expansion_scalars);
    let retained_bytes = peak_scalars
        .checked_mul(size_of::<f64>())
        .ok_or_else(|| refit_structure_error("refit dense basis byte estimate overflows usize"))?;
    if retained_bytes > REFIT_MAX_ALLOC_BYTES {
        return Err(refit_structure_error(format!(
            "refit dense basis retains {retained_bytes} requested bytes above static cap {REFIT_MAX_ALLOC_BYTES}"
        )));
    }
    Ok(RefitDenseBasisPlan {
        total_work,
        peak_scalars,
    })
}

/// Evaluate and expand one admitted nonzero basis row into the complete refit
/// control axis with bounded cancellation polling.
///
/// Aggregate basis-plus-expansion work and simultaneously live scalar payload
/// are admitted before the nested basis scratch allocates. One gate then spans
/// admitted span/basis evaluation, fallible dense-output allocation, zero fill,
/// active-value placement, and final publication. Cancellation exposes neither
/// a partial nonzero basis nor a partial dense row. The primitive does not
/// consume the `Cx` budget or make the full refit pipeline cancellation-aware.
///
/// # Errors
/// Returns the admitted basis parameter/arithmetic/allocation refusals or a
/// structured aggregate work, retained-byte, or dense-allocation refusal.
pub fn evaluate_refit_dense_basis_with_cx(
    knots: AdmittedKnotVector<'_, f64>,
    parameter: f64,
    cx: &Cx<'_>,
) -> Result<RefitDenseBasisRun, NurbsError> {
    let mut should_cancel = || cx.checkpoint().is_err();
    evaluate_refit_dense_basis_with_poll(knots, parameter, &mut should_cancel)
}

fn evaluate_refit_dense_basis_with_poll(
    knots: AdmittedKnotVector<'_, f64>,
    parameter: f64,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<RefitDenseBasisRun, NurbsError> {
    let control_count = knots.control_count();
    let degree = knots.degree();
    preflight_refit_dense_basis(control_count, degree)?;
    let (span, active) = match knots.basis_with_poll(parameter, should_cancel)? {
        BasisRun::Complete { span, values } => (span, values),
        BasisRun::Cancelled => return Ok(RefitDenseBasisRun::Cancelled),
    };

    if should_cancel() {
        return Ok(RefitDenseBasisRun::Cancelled);
    }
    let mut values = try_vec_with_capacity(control_count, "refit dense basis row")?;
    let mut operations_since_poll = 0usize;
    for _ in 0..control_count {
        values.push(0.0);
        if refit_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(RefitDenseBasisRun::Cancelled);
        }
    }
    for (offset, &value) in active.iter().enumerate() {
        values[span - degree + offset] = value;
        if refit_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(RefitDenseBasisRun::Cancelled);
        }
    }
    if should_cancel() {
        return Ok(RefitDenseBasisRun::Cancelled);
    }
    Ok(RefitDenseBasisRun::Complete { values })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RefitSampleBasisRowPlan {
    control_count: usize,
}

fn preflight_refit_sample_basis_row_counts(
    control_count_u: usize,
    degree_u: usize,
    control_count_v: usize,
    degree_v: usize,
) -> Result<RefitSampleBasisRowPlan, NurbsError> {
    let control_count = control_count_u
        .checked_mul(control_count_v)
        .ok_or_else(|| refit_structure_error("refit sample basis row size overflows usize"))?;
    let plan_u = preflight_refit_dense_basis(control_count_u, degree_u)?;
    let plan_v = preflight_refit_dense_basis(control_count_v, degree_v)?;
    let tensor_work =
        checked_refit_work_product(&[control_count as u128, 4], "sample basis tensor product")?;
    let total_work = checked_refit_work_sum(
        &[plan_u.total_work, plan_v.total_work, tensor_work, 3],
        "sample basis row aggregate",
    )?;
    if total_work > REFIT_MAX_WORK_UNITS {
        return Err(refit_structure_error(format!(
            "refit sample basis row requests {total_work} work units above static cap {REFIT_MAX_WORK_UNITS}"
        )));
    }

    let v_phase_scalars = control_count_u
        .checked_add(plan_v.peak_scalars)
        .ok_or_else(|| refit_structure_error("refit sample V-basis live size overflows usize"))?;
    let output_phase_scalars = control_count_u
        .checked_add(control_count_v)
        .and_then(|scalars| scalars.checked_add(control_count))
        .ok_or_else(|| refit_structure_error("refit sample row live size overflows usize"))?;
    let peak_scalars = plan_u
        .peak_scalars
        .max(v_phase_scalars)
        .max(output_phase_scalars);
    let retained_bytes = peak_scalars
        .checked_mul(size_of::<f64>())
        .ok_or_else(|| refit_structure_error("refit sample row byte estimate overflows usize"))?;
    if retained_bytes > REFIT_MAX_ALLOC_BYTES {
        return Err(refit_structure_error(format!(
            "refit sample basis row retains {retained_bytes} requested bytes above static cap {REFIT_MAX_ALLOC_BYTES}"
        )));
    }
    Ok(RefitSampleBasisRowPlan { control_count })
}

fn preflight_refit_sample_basis_row(
    knots_u: AdmittedKnotVector<'_, f64>,
    knots_v: AdmittedKnotVector<'_, f64>,
    u: f64,
    v: f64,
) -> Result<RefitSampleBasisRowPlan, NurbsError> {
    knots_u
        .source()
        .preflight_parameter(u, "refit U sample basis")?;
    knots_v
        .source()
        .preflight_parameter(v, "refit V sample basis")?;
    preflight_refit_sample_basis_row_counts(
        knots_u.control_count(),
        knots_u.degree(),
        knots_v.control_count(),
        knots_v.degree(),
    )
}

/// Assemble one complete tensor-product sample row from admitted U/V knots.
///
/// U-then-V parameter checks, aggregate nested-basis-plus-tensor work, and peak
/// simultaneously live U-row/V-row/output payload are admitted before the first
/// nested basis allocation. One gate spans both admitted basis evaluations,
/// fallible output reservation, deterministic row-major tensor fill, finite
/// arithmetic checks, and final publication. Cancellation drops every local row
/// and exposes no partial sample. The primitive does not consume the `Cx` budget
/// or grant fit, conditioning, or geometric-error authority.
///
/// # Errors
/// Returns the admitted basis parameter/arithmetic/allocation refusals or a
/// structured aggregate work, retained-byte, tensor-size, or output-allocation
/// refusal.
pub fn build_refit_sample_basis_row_with_cx(
    knots_u: AdmittedKnotVector<'_, f64>,
    knots_v: AdmittedKnotVector<'_, f64>,
    u: f64,
    v: f64,
    cx: &Cx<'_>,
) -> Result<RefitSampleBasisRowRun, NurbsError> {
    let mut should_cancel = || cx.checkpoint().is_err();
    build_refit_sample_basis_row_with_poll(knots_u, knots_v, u, v, &mut should_cancel)
}

fn build_refit_sample_basis_row_with_poll(
    knots_u: AdmittedKnotVector<'_, f64>,
    knots_v: AdmittedKnotVector<'_, f64>,
    u: f64,
    v: f64,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<RefitSampleBasisRowRun, NurbsError> {
    let plan = preflight_refit_sample_basis_row(knots_u, knots_v, u, v)?;
    let values_u = match evaluate_refit_dense_basis_with_poll(knots_u, u, should_cancel)? {
        RefitDenseBasisRun::Complete { values } => values,
        RefitDenseBasisRun::Cancelled => return Ok(RefitSampleBasisRowRun::Cancelled),
    };
    let values_v = match evaluate_refit_dense_basis_with_poll(knots_v, v, should_cancel)? {
        RefitDenseBasisRun::Complete { values } => values,
        RefitDenseBasisRun::Cancelled => return Ok(RefitSampleBasisRowRun::Cancelled),
    };
    if should_cancel() {
        return Ok(RefitSampleBasisRowRun::Cancelled);
    }

    let mut values = try_vec_with_capacity(plan.control_count, "refit sample basis row")?;
    if should_cancel() {
        return Ok(RefitSampleBasisRowRun::Cancelled);
    }
    let mut operations_since_poll = 0usize;
    for &value_u in &values_u {
        for &value_v in &values_v {
            let value = if value_u == 0.0 || value_v == 0.0 {
                0.0
            } else {
                value_u * value_v
            };
            if !value.is_finite() {
                return Err(refit_structure_error(
                    "refit sample basis tensor arithmetic became non-finite",
                ));
            }
            values.push(value);
            if refit_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(RefitSampleBasisRowRun::Cancelled);
            }
        }
    }
    debug_assert_eq!(values.len(), plan.control_count);
    if should_cancel() {
        return Ok(RefitSampleBasisRowRun::Cancelled);
    }
    Ok(RefitSampleBasisRowRun::Complete { values })
}

fn refit_sample_basis_row(
    knots_u: AdmittedKnotVector<'_, f64>,
    knots_v: AdmittedKnotVector<'_, f64>,
    u: f64,
    v: f64,
) -> Result<Vec<f64>, NurbsError> {
    let mut never_cancel = || false;
    match build_refit_sample_basis_row_with_poll(knots_u, knots_v, u, v, &mut never_cancel)? {
        RefitSampleBasisRowRun::Complete { values } => Ok(values),
        RefitSampleBasisRowRun::Cancelled => Err(refit_structure_error(
            "non-cancellable refit sample basis assembly observed cancellation",
        )),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RefitRightHandSidesPlan {
    row_count: usize,
    control_count: usize,
}

fn preflight_refit_right_hand_sides_counts(
    row_count: usize,
    target_count: usize,
    control_count: usize,
) -> Result<RefitRightHandSidesPlan, NurbsError> {
    if row_count == 0 {
        return Err(refit_structure_error(
            "refit right-hand-side assembly requires at least one sample row",
        ));
    }
    if row_count != target_count {
        return Err(refit_structure_error(format!(
            "refit right-hand-side has {row_count} sample rows but {target_count} targets"
        )));
    }
    if control_count == 0 {
        return Err(refit_structure_error(
            "refit right-hand-side assembly requires at least one control",
        ));
    }
    let entries = checked_refit_work_product(
        &[row_count as u128, control_count as u128],
        "right-hand-side matrix traversal",
    )?;
    let row_validation =
        checked_refit_work_product(&[row_count as u128, 4], "right-hand-side row validation")?;
    let accumulation = checked_refit_work_product(&[entries, 12], "right-hand-side accumulation")?;
    let output_fill =
        checked_refit_work_product(&[control_count as u128, 3], "right-hand-side output fill")?;
    let total_work = checked_refit_work_sum(
        &[row_validation, entries, output_fill, accumulation, 5],
        "right-hand-side aggregate",
    )?;
    if total_work > REFIT_MAX_WORK_UNITS {
        return Err(refit_structure_error(format!(
            "refit right-hand-side requests {total_work} work units above static cap {REFIT_MAX_WORK_UNITS}"
        )));
    }
    let retained_scalars = control_count.checked_mul(3).ok_or_else(|| {
        refit_structure_error("refit right-hand-side retained-scalar estimate overflows usize")
    })?;
    let retained_bytes = retained_scalars
        .checked_mul(size_of::<f64>())
        .ok_or_else(|| {
            refit_structure_error("refit right-hand-side byte estimate overflows usize")
        })?;
    if retained_bytes > REFIT_MAX_ALLOC_BYTES {
        return Err(refit_structure_error(format!(
            "refit right-hand-sides retain {retained_bytes} requested bytes above static cap {REFIT_MAX_ALLOC_BYTES}"
        )));
    }
    Ok(RefitRightHandSidesPlan {
        row_count,
        control_count,
    })
}

/// Assemble all three coordinate right-hand sides `B^T targets` in one pass.
///
/// Count compatibility, aggregate borrowed validation/three-axis accumulation
/// work, and requested output payload are admitted before the first checkpoint
/// or scan. One gate then spans deterministic row/target/value validation,
/// three fallible zero-output allocations, finite row-major accumulation, and
/// final publication, polling after at most 64 rows, values, or matrix entries.
/// Cancellation drops all partial RHS vectors. Borrowed inputs are excluded from
/// the retained envelope; the primitive does not consume the `Cx` budget or
/// grant solve, conditioning, fit, or geometric authority.
///
/// # Errors
/// Returns a structured dimension, checked work/retained-byte, allocation,
/// borrowed-value, or finite-arithmetic refusal.
pub fn assemble_refit_right_hand_sides_with_cx(
    rows_b: &[Vec<f64>],
    targets: &[[f64; 3]],
    cx: &Cx<'_>,
) -> Result<RefitRightHandSidesRun, NurbsError> {
    let mut should_cancel = || cx.checkpoint().is_err();
    assemble_refit_right_hand_sides_with_poll(rows_b, targets, &mut should_cancel)
}

fn assemble_refit_right_hand_sides_with_poll(
    rows_b: &[Vec<f64>],
    targets: &[[f64; 3]],
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<RefitRightHandSidesRun, NurbsError> {
    let control_count = rows_b.first().map_or(0, Vec::len);
    let plan = preflight_refit_right_hand_sides_counts(rows_b.len(), targets.len(), control_count)?;
    if should_cancel() {
        return Ok(RefitRightHandSidesRun::Cancelled);
    }

    let mut operations_since_poll = 0usize;
    for (row_index, (row, target)) in rows_b.iter().zip(targets).enumerate() {
        if row.len() != plan.control_count {
            return Err(refit_structure_error(format!(
                "refit right-hand-side row {row_index} has length {}, expected {}",
                row.len(),
                plan.control_count
            )));
        }
        if target.iter().any(|coordinate| !coordinate.is_finite()) {
            return Err(refit_structure_error(format!(
                "refit right-hand-side target {row_index} must be finite"
            )));
        }
        if refit_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(RefitRightHandSidesRun::Cancelled);
        }
    }
    if should_cancel() {
        return Ok(RefitRightHandSidesRun::Cancelled);
    }

    operations_since_poll = 0;
    for row in rows_b {
        for value in row {
            if !value.is_finite() {
                return Err(refit_structure_error(
                    "refit right-hand-side sample matrix values must be finite",
                ));
            }
            if refit_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(RefitRightHandSidesRun::Cancelled);
            }
        }
    }
    if should_cancel() {
        return Ok(RefitRightHandSidesRun::Cancelled);
    }

    let rhs_x = try_vec_with_capacity(plan.control_count, "refit X right-hand side")?;
    if should_cancel() {
        return Ok(RefitRightHandSidesRun::Cancelled);
    }
    let rhs_y = try_vec_with_capacity(plan.control_count, "refit Y right-hand side")?;
    if should_cancel() {
        return Ok(RefitRightHandSidesRun::Cancelled);
    }
    let rhs_z = try_vec_with_capacity(plan.control_count, "refit Z right-hand side")?;
    if should_cancel() {
        return Ok(RefitRightHandSidesRun::Cancelled);
    }
    let mut right_hand_sides = [rhs_x, rhs_y, rhs_z];
    operations_since_poll = 0;
    for right_hand_side in &mut right_hand_sides {
        for _ in 0..plan.control_count {
            right_hand_side.push(0.0);
            if refit_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(RefitRightHandSidesRun::Cancelled);
            }
        }
    }

    operations_since_poll = 0;
    for (row, target) in rows_b.iter().zip(targets) {
        for (control_index, &weight) in row.iter().enumerate() {
            if weight != 0.0 {
                for axis in 0..3 {
                    let contribution = weight * target[axis];
                    let accumulated = right_hand_sides[axis][control_index] + contribution;
                    if !contribution.is_finite() || !accumulated.is_finite() {
                        return Err(refit_structure_error(
                            "refit right-hand-side accumulation became non-finite",
                        ));
                    }
                    right_hand_sides[axis][control_index] = accumulated;
                }
            }
            if refit_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(RefitRightHandSidesRun::Cancelled);
            }
        }
    }
    debug_assert_eq!(rows_b.len(), plan.row_count);
    if should_cancel() {
        return Ok(RefitRightHandSidesRun::Cancelled);
    }
    Ok(RefitRightHandSidesRun::Complete { right_hand_sides })
}

fn assemble_refit_right_hand_sides(
    rows_b: &[Vec<f64>],
    targets: &[[f64; 3]],
) -> Result<[Vec<f64>; 3], NurbsError> {
    let mut never_cancel = || false;
    match assemble_refit_right_hand_sides_with_poll(rows_b, targets, &mut never_cancel)? {
        RefitRightHandSidesRun::Complete { right_hand_sides } => Ok(right_hand_sides),
        RefitRightHandSidesRun::Cancelled => Err(refit_structure_error(
            "non-cancellable refit right-hand-side assembly observed cancellation",
        )),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RefitResidualSummaryPlan {
    row_count: usize,
    control_count: usize,
    net_rows: usize,
    net_cols: usize,
}

fn preflight_refit_residual_summary_counts(
    row_count: usize,
    target_count: usize,
    uv_count: usize,
    control_count: usize,
    net_rows: usize,
    net_cols: usize,
    warn_residual: f64,
) -> Result<RefitResidualSummaryPlan, NurbsError> {
    if !warn_residual.is_finite() || warn_residual < 0.0 {
        return Err(refit_structure_error(
            "refit residual warning threshold must be finite and non-negative",
        ));
    }
    if row_count == 0 {
        return Err(refit_structure_error(
            "refit residual summary requires at least one sample row",
        ));
    }
    if row_count != target_count || row_count != uv_count {
        return Err(refit_structure_error(format!(
            "refit residual summary has {row_count} sample rows, {target_count} targets, and {uv_count} parameter pairs"
        )));
    }
    if control_count == 0 || net_rows == 0 || net_cols == 0 {
        return Err(refit_structure_error(
            "refit residual summary requires a nonempty sample row and control net",
        ));
    }
    let net_control_count = net_rows
        .checked_mul(net_cols)
        .ok_or_else(|| refit_structure_error("refit residual control-net size overflows usize"))?;
    if net_control_count != control_count {
        return Err(refit_structure_error(format!(
            "refit residual sample rows have {control_count} controls but the control net has {net_control_count}"
        )));
    }

    let entries = checked_refit_work_product(
        &[row_count as u128, control_count as u128],
        "residual sample-matrix traversal",
    )?;
    let row_validation =
        checked_refit_work_product(&[row_count as u128, 8], "residual row validation")?;
    let net_validation = checked_refit_work_product(
        &[net_control_count as u128, 4],
        "residual control-net validation",
    )?;
    // Each matrix entry performs deterministic index traversal and at most
    // three multiply-add pairs. The coefficient intentionally overprices that
    // fixed scalar work.
    let evaluation = checked_refit_work_product(&[entries, 8], "residual point evaluation")?;
    let report_work = checked_refit_work_product(&[row_count as u128, 24], "residual report")?;
    let total_work = checked_refit_work_sum(
        &[
            row_validation,
            entries,
            net_validation,
            evaluation,
            report_work,
            4,
        ],
        "residual summary aggregate",
    )?;
    if total_work > REFIT_MAX_WORK_UNITS {
        return Err(refit_structure_error(format!(
            "refit residual summary requests {total_work} work units above static cap {REFIT_MAX_WORK_UNITS}"
        )));
    }

    let warning_bytes = row_count
        .checked_mul(size_of::<LocalizedFitResidualWarning>())
        .ok_or_else(|| refit_structure_error("refit residual warning bytes overflow usize"))?;
    if warning_bytes > REFIT_MAX_ALLOC_BYTES {
        return Err(refit_structure_error(format!(
            "refit residual warnings retain {warning_bytes} requested bytes above static cap {REFIT_MAX_ALLOC_BYTES}"
        )));
    }

    Ok(RefitResidualSummaryPlan {
        row_count,
        control_count,
        net_rows,
        net_cols,
    })
}

/// Summarize paired sample-target residuals against a fitted control net.
///
/// Count compatibility, aggregate borrowed validation/evaluation work, and the
/// worst-case warning payload are admitted before the first checkpoint or
/// scan. One caller gate then spans rectangular and finite borrowed-input
/// validation, fallible warning reservation, deterministic U-major point
/// evaluation, finite RMS/max accumulation, warning construction, and final
/// publication. Polls occur after at most 64 rows, controls, or matrix entries;
/// cancellation drops all provisional statistics and warnings.
///
/// Borrowed matrix, target, parameter, and control-net storage is excluded from
/// the retained envelope. This measured primitive does not consume the `Cx`
/// budget or grant field, zero-set, Hausdorff, topology, or continuum-fit
/// authority.
///
/// # Errors
/// Returns a structured threshold, dimension, checked work/retained-byte,
/// allocation, borrowed-value, or finite-arithmetic refusal.
pub fn summarize_refit_residuals_with_cx(
    rows_b: &[Vec<f64>],
    targets: &[[f64; 3]],
    uvs: &[[f64; 2]],
    net: &[Vec<[f64; 3]>],
    warn_residual: f64,
    cx: &Cx<'_>,
) -> Result<RefitResidualSummaryRun, NurbsError> {
    let mut should_cancel = || cx.checkpoint().is_err();
    summarize_refit_residuals_with_poll(
        rows_b,
        targets,
        uvs,
        net,
        warn_residual,
        &mut should_cancel,
    )
}

fn summarize_refit_residuals_with_poll(
    rows_b: &[Vec<f64>],
    targets: &[[f64; 3]],
    uvs: &[[f64; 2]],
    net: &[Vec<[f64; 3]>],
    warn_residual: f64,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<RefitResidualSummaryRun, NurbsError> {
    let control_count = rows_b.first().map_or(0, Vec::len);
    let net_cols = net.first().map_or(0, Vec::len);
    let plan = preflight_refit_residual_summary_counts(
        rows_b.len(),
        targets.len(),
        uvs.len(),
        control_count,
        net.len(),
        net_cols,
        warn_residual,
    )?;
    if should_cancel() {
        return Ok(RefitResidualSummaryRun::Cancelled);
    }

    let mut operations_since_poll = 0usize;
    for (row_index, ((row, target), uv)) in rows_b.iter().zip(targets).zip(uvs).enumerate() {
        if row.len() != plan.control_count {
            return Err(refit_structure_error(format!(
                "refit residual row {row_index} has length {}, expected {}",
                row.len(),
                plan.control_count
            )));
        }
        if target.iter().any(|coordinate| !coordinate.is_finite()) {
            return Err(refit_structure_error(format!(
                "refit residual target {row_index} must be finite"
            )));
        }
        if uv.iter().any(|coordinate| !coordinate.is_finite()) {
            return Err(refit_structure_error(format!(
                "refit residual parameter pair {row_index} must be finite"
            )));
        }
        if refit_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(RefitResidualSummaryRun::Cancelled);
        }
    }
    if should_cancel() {
        return Ok(RefitResidualSummaryRun::Cancelled);
    }

    operations_since_poll = 0;
    for (row_index, row) in net.iter().enumerate() {
        if row.len() != plan.net_cols {
            return Err(refit_structure_error(format!(
                "refit residual control-net row {row_index} has length {}, expected {}",
                row.len(),
                plan.net_cols
            )));
        }
        for control in row {
            if control.iter().any(|coordinate| !coordinate.is_finite()) {
                return Err(refit_structure_error(format!(
                    "refit residual control-net row {row_index} must be finite"
                )));
            }
            if refit_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(RefitResidualSummaryRun::Cancelled);
            }
        }
    }
    if should_cancel() {
        return Ok(RefitResidualSummaryRun::Cancelled);
    }

    operations_since_poll = 0;
    for row in rows_b {
        for value in row {
            if !value.is_finite() {
                return Err(refit_structure_error(
                    "refit residual sample matrix values must be finite",
                ));
            }
            if refit_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(RefitResidualSummaryRun::Cancelled);
            }
        }
    }
    if should_cancel() {
        return Ok(RefitResidualSummaryRun::Cancelled);
    }

    let mut warnings = try_vec_with_capacity(plan.row_count, "refit residual warnings")?;
    if should_cancel() {
        return Ok(RefitResidualSummaryRun::Cancelled);
    }

    let mut rms_square_sum = 0.0f64;
    let mut max_residual = 0.0f64;
    operations_since_poll = 0;
    for ((row, target), uv) in rows_b.iter().zip(targets).zip(uvs) {
        let mut point = [0.0f64; 3];
        for (control_index, &weight) in row.iter().enumerate() {
            if weight != 0.0 {
                let i = control_index / plan.net_cols;
                let j = control_index % plan.net_cols;
                for axis in 0..3 {
                    let contribution = weight * net[i][j][axis];
                    let accumulated = point[axis] + contribution;
                    if !contribution.is_finite() || !accumulated.is_finite() {
                        return Err(refit_structure_error(
                            "refit residual point evaluation became non-finite",
                        ));
                    }
                    point[axis] = accumulated;
                }
            }
            if refit_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(RefitResidualSummaryRun::Cancelled);
            }
        }

        let residual = norm3([
            point[0] - target[0],
            point[1] - target[1],
            point[2] - target[2],
        ]);
        if !residual.is_finite() {
            return Err(refit_structure_error(
                "refit residual arithmetic is non-finite",
            ));
        }
        let squared = residual * residual;
        let accumulated = rms_square_sum + squared;
        if !squared.is_finite() || !accumulated.is_finite() {
            return Err(refit_structure_error(
                "refit RMS accumulation is non-finite",
            ));
        }
        rms_square_sum = accumulated;
        max_residual = max_residual.max(residual);
        if residual > warn_residual {
            warnings.push(LocalizedFitResidualWarning {
                uv: *uv,
                point: *target,
                residual,
            });
        }
        if refit_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(RefitResidualSummaryRun::Cancelled);
        }
    }
    debug_assert_eq!(net.len(), plan.net_rows);
    debug_assert!(warnings.len() <= plan.row_count);
    #[allow(clippy::cast_precision_loss)]
    let rms_residual = det::sqrt(rms_square_sum / plan.row_count as f64);
    if !rms_residual.is_finite() {
        return Err(refit_structure_error("refit RMS result is non-finite"));
    }
    if should_cancel() {
        return Ok(RefitResidualSummaryRun::Cancelled);
    }
    Ok(RefitResidualSummaryRun::Complete {
        rms_residual,
        max_residual,
        warnings,
    })
}

fn summarize_refit_residuals(
    rows_b: &[Vec<f64>],
    targets: &[[f64; 3]],
    uvs: &[[f64; 2]],
    net: &[Vec<[f64; 3]>],
    warn_residual: f64,
) -> Result<(f64, f64, Vec<LocalizedFitResidualWarning>), NurbsError> {
    let mut never_cancel = || false;
    match summarize_refit_residuals_with_poll(
        rows_b,
        targets,
        uvs,
        net,
        warn_residual,
        &mut never_cancel,
    )? {
        RefitResidualSummaryRun::Complete {
            rms_residual,
            max_residual,
            warnings,
        } => Ok((rms_residual, max_residual, warnings)),
        RefitResidualSummaryRun::Cancelled => Err(refit_structure_error(
            "non-cancellable refit residual summary observed cancellation",
        )),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RefitFieldProbePlan {
    probe: usize,
}

fn preflight_refit_field_probe_envelope(
    probe: usize,
    evaluation_work: u128,
    evaluation_retained_bytes: u128,
) -> Result<RefitFieldProbePlan, NurbsError> {
    if probe == 0 {
        return Err(refit_structure_error(
            "refit field probe resolution must be positive",
        ));
    }
    let probe_points = probe
        .checked_mul(probe)
        .ok_or_else(|| refit_structure_error("refit field probe grid size overflows usize"))?;
    if probe_points > REFIT_MAX_PROBE_POINTS {
        return Err(refit_structure_error(format!(
            "refit field probe grid {probe_points} exceeds static cap {REFIT_MAX_PROBE_POINTS}"
        )));
    }
    let nested_work = checked_refit_work_product(
        &[probe_points as u128, evaluation_work],
        "field-probe surface evaluations",
    )?;
    // Grid parameter construction, pre/post opaque-call gates, finite/absolute
    // reduction, and loop bookkeeping. The opaque closure's own cost cannot be
    // inferred from this API and remains outside the measured envelope.
    let reduction_work =
        checked_refit_work_product(&[probe_points as u128, 16], "field-probe reduction")?;
    let total_work =
        checked_refit_work_sum(&[nested_work, reduction_work, 2], "field-probe aggregate")?;
    if total_work > REFIT_MAX_WORK_UNITS {
        return Err(refit_structure_error(format!(
            "refit field probe requests {total_work} work units above static cap {REFIT_MAX_WORK_UNITS}"
        )));
    }
    if evaluation_retained_bytes > REFIT_MAX_ALLOC_BYTES as u128 {
        return Err(refit_structure_error(format!(
            "refit field probe evaluation retains {evaluation_retained_bytes} requested bytes above static cap {REFIT_MAX_ALLOC_BYTES}"
        )));
    }
    Ok(RefitFieldProbePlan { probe })
}

fn preflight_refit_field_probe(
    surface: AdmittedNurbsSurface<'_, f64>,
    probe: usize,
) -> Result<RefitFieldProbePlan, NurbsError> {
    let (evaluation_work, evaluation_retained_bytes) = surface.evaluation_envelope()?;
    preflight_refit_field_probe_envelope(probe, evaluation_work, evaluation_retained_bytes)
}

/// Measure `max(abs(field(surface(u, v))))` on a dense cell-center grid.
///
/// Grid size, all repeated admitted-surface evaluation work, and the peak
/// requested nested basis payload are composed before the first checkpoint.
/// One caller gate then spans deterministic row-major parameters, every nested
/// surface evaluation, cancellation immediately before each opaque field call,
/// finite absolute-value reduction, and final scalar publication. Cancellation
/// exposes no provisional maximum.
///
/// The field closure itself is opaque: one in-flight call is non-preemptible
/// and its cost, side effects, continuity, units, Lipschitz behavior, zero-set
/// meaning, and distance authority are outside this measured primitive. The
/// method does not consume the `Cx` budget or own request drain/finalize.
///
/// # Errors
/// Returns a structured grid/count/work/retained-byte refusal, a nested surface
/// evaluation error, or a non-finite field-value refusal.
pub fn measure_refit_field_probe_with_cx(
    surface: AdmittedNurbsSurface<'_, f64>,
    field: &dyn Fn([f64; 3]) -> f64,
    probe: usize,
    cx: &Cx<'_>,
) -> Result<RefitFieldProbeRun, NurbsError> {
    let mut should_cancel = || cx.checkpoint().is_err();
    measure_refit_field_probe_with_poll(surface, field, probe, &mut should_cancel)
}

fn measure_refit_field_probe_with_poll(
    surface: AdmittedNurbsSurface<'_, f64>,
    field: &dyn Fn([f64; 3]) -> f64,
    probe: usize,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<RefitFieldProbeRun, NurbsError> {
    let plan = preflight_refit_field_probe(surface, probe)?;
    if should_cancel() {
        return Ok(RefitFieldProbeRun::Cancelled);
    }

    let mut max_abs_field = 0.0f64;
    for a in 0..plan.probe {
        for b in 0..plan.probe {
            #[allow(clippy::cast_precision_loss)]
            let (u, v) = (
                (a as f64 + 0.5) / plan.probe as f64,
                (b as f64 + 0.5) / plan.probe as f64,
            );
            let point = match surface.eval_with_poll(u, v, should_cancel)? {
                SurfaceEvaluationRun::Complete { point } => point,
                SurfaceEvaluationRun::Cancelled => return Ok(RefitFieldProbeRun::Cancelled),
            };
            if point.iter().any(|coordinate| !coordinate.is_finite()) {
                return Err(refit_structure_error(format!(
                    "fitted surface returned a non-finite probe point at ({u}, {v})"
                )));
            }
            if should_cancel() {
                return Ok(RefitFieldProbeRun::Cancelled);
            }
            let field_value = field(point);
            if !field_value.is_finite() {
                return Err(refit_structure_error(format!(
                    "implicit field returned non-finite value at probe {point:?}"
                )));
            }
            max_abs_field = max_abs_field.max(field_value.abs());
            if should_cancel() {
                return Ok(RefitFieldProbeRun::Cancelled);
            }
        }
    }
    if should_cancel() {
        return Ok(RefitFieldProbeRun::Cancelled);
    }
    Ok(RefitFieldProbeRun::Complete { max_abs_field })
}

fn measure_refit_field_probe(
    surface: AdmittedNurbsSurface<'_, f64>,
    field: &dyn Fn([f64; 3]) -> f64,
    probe: usize,
) -> Result<f64, NurbsError> {
    let mut never_cancel = || false;
    match measure_refit_field_probe_with_poll(surface, field, probe, &mut never_cancel)? {
        RefitFieldProbeRun::Complete { max_abs_field } => Ok(max_abs_field),
        RefitFieldProbeRun::Cancelled => Err(refit_structure_error(
            "non-cancellable refit field probe observed cancellation",
        )),
    }
}

/// Row of basis values over the whole control axis (dense, small).
fn basis_row(kv: AdmittedKnotVector<'_, f64>, t: f64) -> Result<Vec<f64>, NurbsError> {
    let mut never_cancel = || false;
    match evaluate_refit_dense_basis_with_poll(kv, t, &mut never_cancel)? {
        RefitDenseBasisRun::Complete { values } => Ok(values),
        RefitDenseBasisRun::Cancelled => Err(refit_structure_error(
            "non-cancellable refit dense basis observed cancellation",
        )),
    }
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
            let row = refit_sample_basis_row(admitted_ku, admitted_kv, u, v)?;
            debug_assert_eq!(row.len(), control_points);
            rows_b.push(row);
            uvs.push([u, v]);
        }
    }
    // Assemble and factor once, then solve the three coordinate right-hand
    // sides against the same deterministic factor.
    let mut net = try_filled_matrix(nu, nv, [0.0f64; 3], "refit control net")?;
    let factor = cholesky_factor(assemble_normal(&rows_b, nu, nv, config.lambda)?)?;
    let right_hand_sides = assemble_refit_right_hand_sides(&rows_b, &targets)?;
    for (axis, rhs) in right_hand_sides.into_iter().enumerate() {
        debug_assert_eq!(rhs.len(), control_points);
        let rhs = cholesky_solve_factored(&factor, rhs)?;
        for i in 0..nu {
            for j in 0..nv {
                net[i][j][axis] = rhs[i * nv + j];
            }
        }
    }
    // EXACT G0 seam closure: tie the u-boundary control rows through the same
    // aggregate-admitted primitive exposed to successor orchestration.
    let net = close_refit_u_seam(net)?;
    let weights = try_filled_matrix(nu, nv, 1.0f64, "refit unit weights")?;
    let surface = NurbsSurface::new(ku, kv, &net, &weights)?;
    let report_surface = surface.admit()?;
    // ---- The honest report -------------------------------------------
    let (rms_residual, max_res, warnings) =
        summarize_refit_residuals(&rows_b, &targets, &uvs, &net, config.warn_residual)?;
    // Spline → field: dense probe plus an analytic-model Lipschitz estimate;
    // the other reported direction stays the sampled fit-target worst case and
    // does not claim that a generic closure's targets belong to a source set.
    let sampled = measure_refit_field_probe(report_surface, field, probe)?;
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

fn preflight_refit_seam_g1_envelope(
    partial_work: u128,
    partial_retained_bytes: u128,
) -> Result<(), NurbsError> {
    let evaluation_count = checked_refit_work_product(
        &[REFIT_SEAM_G1_SAMPLE_COUNT as u128, 2],
        "seam G1 partial evaluations",
    )?;
    let partials_work = checked_refit_work_product(
        &[evaluation_count, partial_work],
        "seam G1 partial evaluations",
    )?;
    // Per sample: two nested preflight seams plus stable norms,
    // normalization, dot product, clamp, degeneracy bookkeeping, and polling.
    let diagnostic_work = checked_refit_work_product(
        &[REFIT_SEAM_G1_SAMPLE_COUNT as u128, 48],
        "seam G1 diagnostic",
    )?;
    let total_work =
        checked_refit_work_sum(&[partials_work, diagnostic_work, 2], "seam G1 aggregate")?;
    if total_work > REFIT_MAX_WORK_UNITS {
        return Err(refit_structure_error(format!(
            "refit seam G1 diagnostic requests {total_work} work units above static cap {REFIT_MAX_WORK_UNITS}"
        )));
    }
    if partial_retained_bytes > REFIT_MAX_ALLOC_BYTES as u128 {
        return Err(refit_structure_error(format!(
            "refit seam G1 partials retain {partial_retained_bytes} requested bytes above static cap {REFIT_MAX_ALLOC_BYTES}"
        )));
    }
    Ok(())
}

fn preflight_refit_seam_g1(surface: AdmittedNurbsSurface<'_, f64>) -> Result<(), NurbsError> {
    let (partial_work, partial_retained_bytes) = surface.partials_envelope()?;
    preflight_refit_seam_g1_envelope(partial_work, partial_retained_bytes)
}

/// Measure the fitted surface's fixed interior paired-side G1 seam diagnostic.
///
/// The exact envelope for one admitted surface-partials evaluation is composed
/// across all 46 nested evaluations before the first checkpoint. One caller
/// gate then spans U=0/U=1 partials at 23 deterministic interior V parameters,
/// stable tangent normalization, degeneracy accounting, and final tuple
/// publication. Nested evaluators retain their own parameter, continuity,
/// allocation, and finite-arithmetic refusals. Cancellation publishes neither
/// the running maximum nor a partial degenerate count.
///
/// The primitive allocates no owned numerical output. Its ordinary-f64 result
/// is a measured diagnostic that deliberately excludes V endpoints and grants
/// no G1 continuity, pole, topology, geometric-error, or certification
/// authority. It does not consume the `Cx` budget or own drain/finalize.
///
/// # Errors
/// Returns the admitted partial evaluator's structured refusals or an aggregate
/// checked work/retained-byte or finite seam-arithmetic refusal.
pub fn measure_refit_seam_g1_with_cx(
    surface: AdmittedNurbsSurface<'_, f64>,
    cx: &Cx<'_>,
) -> Result<RefitSeamG1DiagnosticRun, NurbsError> {
    let mut should_cancel = || cx.checkpoint().is_err();
    measure_refit_seam_g1_with_poll(surface, &mut should_cancel)
}

fn measure_refit_seam_g1_with_poll(
    surface: AdmittedNurbsSurface<'_, f64>,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<RefitSeamG1DiagnosticRun, NurbsError> {
    preflight_refit_seam_g1(surface)?;
    if should_cancel() {
        return Ok(RefitSeamG1DiagnosticRun::Cancelled);
    }

    // Compare u-tangents across the exactly closed seam. Normalize each tangent
    // separately to avoid overflow/underflow in n0*n1, and clamp the rounded dot
    // product to the mathematical cosine range.
    let mut seam_g1 = 0.0f64;
    let mut degenerate = 0usize;
    // `v=0` and `v=1` are deliberately outside this diagnostic's typed scope:
    // radial/revolution fits commonly collapse them to poles where a tangent
    // direction is undefined. The report makes that exclusion machine-visible
    // pending a chart-aware pole audit.
    for b in 1..=REFIT_SEAM_G1_SAMPLE_COUNT {
        #[allow(clippy::cast_precision_loss)]
        let v = b as f64 / (REFIT_SEAM_G1_SAMPLE_COUNT + 1) as f64;
        let (_, du0, _) = match surface.partials_with_poll(0.0, v, should_cancel)? {
            SurfacePartialsRun::Complete { partials } => partials,
            SurfacePartialsRun::Cancelled => return Ok(RefitSeamG1DiagnosticRun::Cancelled),
        };
        let (_, du1, _) = match surface.partials_with_poll(1.0, v, should_cancel)? {
            SurfacePartialsRun::Complete { partials } => partials,
            SurfacePartialsRun::Cancelled => return Ok(RefitSeamG1DiagnosticRun::Cancelled),
        };
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
            if should_cancel() {
                return Ok(RefitSeamG1DiagnosticRun::Cancelled);
            }
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
        if should_cancel() {
            return Ok(RefitSeamG1DiagnosticRun::Cancelled);
        }
    }
    if should_cancel() {
        return Ok(RefitSeamG1DiagnosticRun::Cancelled);
    }
    Ok(RefitSeamG1DiagnosticRun::Complete {
        max_deviation: seam_g1,
        degenerate_samples: degenerate,
    })
}

fn seam_g1_diagnostic_admitted(
    surface: AdmittedNurbsSurface<'_, f64>,
) -> Result<(f64, usize), NurbsError> {
    let mut never_cancel = || false;
    match measure_refit_seam_g1_with_poll(surface, &mut never_cancel)? {
        RefitSeamG1DiagnosticRun::Complete {
            max_deviation,
            degenerate_samples,
        } => Ok((max_deviation, degenerate_samples)),
        RefitSeamG1DiagnosticRun::Cancelled => Err(refit_structure_error(
            "non-cancellable refit seam G1 diagnostic observed cancellation",
        )),
    }
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

    // G0/G5: the cancellable primitive preserves the exact open-uniform
    // construction, including deterministic interior-knot placement.
    #[test]
    fn open_uniform_refit_knots_with_cx_match_the_expected_sequence() {
        let expected = [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0];
        with_refit_cx(false, |cx| {
            let RefitOpenUniformKnotRun::Complete { knots } =
                build_open_uniform_refit_knots_with_cx(4, 2, cx)
                    .expect("cancellable open-uniform knots")
            else {
                panic!("healthy construction must complete");
            };
            assert_eq!(knots.degree(), 2);
            assert_eq!(knots.knots(), expected);
        });

        let legacy = open_uniform_knots(4, 2).expect("legacy open-uniform knots");
        assert_eq!(legacy.degree(), 2);
        assert_eq!(legacy.knots(), expected);
    }

    // G4/G5: one callback spans allocation, fixed-stride fill, sealed
    // validation, and final publication without leaking partial ownership.
    #[test]
    fn open_uniform_refit_knot_polling_is_bounded_and_transactional() {
        let run = |cancel_at: Option<usize>| {
            let polls = Cell::new(0usize);
            let outcome = build_open_uniform_refit_knots_with_poll(130, 1, &mut || {
                polls.set(polls.get() + 1);
                cancel_at == Some(polls.get())
            })
            .expect("valid open-uniform knot construction");
            (outcome, polls.get())
        };

        let (complete, healthy_polls) = run(None);
        assert!(matches!(complete, RefitOpenUniformKnotRun::Complete { .. }));
        assert!(healthy_polls > 4);
        for cancel_at in [1, healthy_polls / 2, healthy_polls] {
            assert_eq!(
                run(Some(cancel_at)),
                (RefitOpenUniformKnotRun::Cancelled, cancel_at)
            );
        }
        with_refit_cx(true, |cx| {
            assert_eq!(
                build_open_uniform_refit_knots_with_cx(130, 1, cx)
                    .expect("pre-cancellation is terminal"),
                RefitOpenUniformKnotRun::Cancelled
            );
        });
    }

    // G0/G4: shape and checked-arithmetic refusals win before the first
    // cancellation observation and therefore before allocation.
    #[test]
    fn open_uniform_refit_knot_refusals_are_preflighted() {
        for (control_count, degree, expected) in [
            (4, 0, "degree must be positive"),
            (usize::MAX, 1, "knot count overflow"),
        ] {
            let polls = Cell::new(0usize);
            let error =
                build_open_uniform_refit_knots_with_poll(control_count, degree, &mut || {
                    polls.set(polls.get() + 1);
                    true
                })
                .expect_err("invalid knot request must be refused");
            assert!(
                error.to_string().contains(expected),
                "unexpected refusal for ({control_count}, {degree}): {error}"
            );
            assert_eq!(polls.get(), 0);
        }
    }

    // G0/G5: the consuming primitive preserves deterministic componentwise
    // midpoint tying and leaves every interior row unchanged.
    #[test]
    fn refit_u_seam_closure_matches_exact_midpoints() {
        let source = vec![
            vec![[0.0, 0.0, 0.0], [2.0, 4.0, 6.0]],
            vec![[9.0, 8.0, 7.0], [6.0, 5.0, 4.0]],
            vec![[2.0, 2.0, 2.0], [4.0, 6.0, 8.0]],
        ];
        let expected = vec![
            vec![[1.0, 1.0, 1.0], [3.0, 5.0, 7.0]],
            source[1].clone(),
            vec![[1.0, 1.0, 1.0], [3.0, 5.0, 7.0]],
        ];
        with_refit_cx(false, |cx| {
            assert_eq!(
                close_refit_u_seam_with_cx(source.clone(), cx)
                    .expect("cancellable exact seam closure"),
                RefitSeamClosureRun::Complete {
                    net: expected.clone()
                }
            );
        });
        assert_eq!(
            close_refit_u_seam(source).expect("synchronous exact seam closure"),
            expected
        );
    }

    // G4/G5: shape validation, finite scanning, mutation, and owner publication
    // share replayable checkpoints; a partially tied consumed net never escapes.
    #[test]
    fn refit_u_seam_closure_polling_is_bounded_and_transactional() {
        let source = vec![vec![[0.0; 3]; 130], vec![[2.0; 3]; 130]];
        let run = |cancel_at: Option<usize>| {
            let polls = Cell::new(0usize);
            let outcome = close_refit_u_seam_with_poll(source.clone(), &mut || {
                polls.set(polls.get() + 1);
                cancel_at == Some(polls.get())
            })
            .expect("valid seam closure");
            (outcome, polls.get())
        };

        let (complete, healthy_polls) = run(None);
        assert!(matches!(complete, RefitSeamClosureRun::Complete { .. }));
        assert!(healthy_polls > 6);
        for cancel_at in [1, healthy_polls / 2, healthy_polls] {
            assert_eq!(
                run(Some(cancel_at)),
                (RefitSeamClosureRun::Cancelled, cancel_at)
            );
        }
        with_refit_cx(true, |cx| {
            assert_eq!(
                close_refit_u_seam_with_cx(source, cx).expect("pre-cancellation is terminal"),
                RefitSeamClosureRun::Cancelled
            );
        });
    }

    // G0/G4: count refusals precede polling; malformed rectangular structure
    // and non-finite controls refuse before any boundary mutation.
    #[test]
    fn refit_u_seam_closure_refusals_are_typed_and_preflighted() {
        let overflow = preflight_refit_u_seam_closure(usize::MAX, 2)
            .expect_err("unrepresentable grid must be refused");
        assert!(overflow.to_string().contains("control-grid size overflow"));

        let polls = Cell::new(0usize);
        let too_short = close_refit_u_seam_with_poll(vec![vec![[0.0; 3]]], &mut || {
            polls.set(polls.get() + 1);
            true
        })
        .expect_err("single-row grid must be refused");
        assert!(too_short.to_string().contains("at least two control rows"));
        assert_eq!(polls.get(), 0);

        let malformed =
            close_refit_u_seam_with_poll(vec![vec![[0.0; 3]; 2], vec![[1.0; 3]]], &mut || false)
                .expect_err("ragged grid must be refused");
        assert!(malformed.to_string().contains("row 1 has length 1"));

        let nonfinite = close_refit_u_seam_with_poll(
            vec![vec![[0.0; 3]], vec![[f64::NAN, 0.0, 0.0]]],
            &mut || false,
        )
        .expect_err("non-finite grid must be refused");
        assert!(nonfinite.to_string().contains("must be finite"));
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
        with_refit_cx(false, |cx| {
            assert_eq!(
                evaluate_refit_dense_basis_with_cx(
                    knots.admitted_after_validation(),
                    parameter,
                    cx,
                )
                .expect("cancellable dense basis row"),
                RefitDenseBasisRun::Complete {
                    values: expected.clone(),
                }
            );
        });

        let owning_error = knots
            .basis(-0.25)
            .expect_err("owning out-of-domain refusal");
        let admitted_error = basis_row(knots.admitted_after_validation(), -0.25)
            .expect_err("admitted out-of-domain refusal");
        assert_eq!(admitted_error, owning_error);
    }

    // G4/G5: the same callback spans the admitted basis core, dense allocation,
    // fixed-stride zero fill, active placement, and final publication.
    #[test]
    fn refit_dense_basis_polling_is_bounded_and_transactional() {
        let knots = open_uniform_knots(130, 1).expect("long linear knots");
        let admitted = knots.admitted_after_validation();
        let parameter = 0.75;
        let basis_polls = Cell::new(0usize);
        assert!(matches!(
            admitted
                .basis_with_poll(parameter, &mut || {
                    basis_polls.set(basis_polls.get() + 1);
                    false
                })
                .expect("admitted basis"),
            BasisRun::Complete { .. }
        ));

        let run = |cancel_at: Option<usize>| {
            let polls = Cell::new(0usize);
            let outcome = evaluate_refit_dense_basis_with_poll(admitted, parameter, &mut || {
                polls.set(polls.get() + 1);
                cancel_at == Some(polls.get())
            })
            .expect("valid dense basis expansion");
            (outcome, polls.get())
        };
        let (complete, healthy_polls) = run(None);
        assert!(matches!(complete, RefitDenseBasisRun::Complete { .. }));
        let first_dense_poll = basis_polls.get() + 1;
        assert!(healthy_polls > first_dense_poll);
        assert_eq!(
            run(Some(first_dense_poll)),
            (RefitDenseBasisRun::Cancelled, first_dense_poll)
        );
        assert_eq!(
            run(Some(healthy_polls)),
            (RefitDenseBasisRun::Cancelled, healthy_polls)
        );
        with_refit_cx(true, |cx| {
            assert_eq!(
                evaluate_refit_dense_basis_with_cx(admitted, parameter, cx)
                    .expect("pre-cancellation is terminal"),
                RefitDenseBasisRun::Cancelled
            );
        });
    }

    // G0/G4: aggregate hostile work is refused before the nested basis core can
    // allocate any scratch storage.
    #[test]
    fn refit_dense_basis_refuses_unbounded_aggregate_work() {
        let error = preflight_refit_dense_basis(usize::MAX, 1)
            .expect_err("unbounded dense basis work must be refused");
        assert!(error.to_string().contains("work units"));
    }

    // G0/G5: tensor composition preserves the exact U-major row produced by
    // independently evaluated admitted dense basis rows.
    #[test]
    fn refit_sample_basis_row_matches_axis_tensor_product() {
        let knots_u = open_uniform_knots(4, 2).expect("quadratic U knots");
        let knots_v = open_uniform_knots(3, 1).expect("linear V knots");
        let admitted_u = knots_u.admitted_after_validation();
        let admitted_v = knots_v.admitted_after_validation();
        let (u, v) = (0.375, 0.25);
        let values_u = basis_row(admitted_u, u).expect("dense U basis");
        let values_v = basis_row(admitted_v, v).expect("dense V basis");
        let mut expected = Vec::with_capacity(values_u.len() * values_v.len());
        for &value_u in &values_u {
            for &value_v in &values_v {
                expected.push(if value_u == 0.0 || value_v == 0.0 {
                    0.0
                } else {
                    value_u * value_v
                });
            }
        }

        assert_eq!(
            refit_sample_basis_row(admitted_u, admitted_v, u, v).expect("synchronous sample row"),
            expected
        );
        with_refit_cx(false, |cx| {
            assert_eq!(
                build_refit_sample_basis_row_with_cx(admitted_u, admitted_v, u, v, cx)
                    .expect("cancellable sample row"),
                RefitSampleBasisRowRun::Complete {
                    values: expected.clone()
                }
            );
        });
    }

    // G4/G5: one callback crosses both admitted basis phases, output
    // reservation/fill, and final publication without exposing a partial row.
    #[test]
    fn refit_sample_basis_row_polling_is_bounded_and_transactional() {
        let knots = open_uniform_knots(130, 1).expect("long linear knots");
        let admitted = knots.admitted_after_validation();
        let run = |cancel_at: Option<usize>| {
            let polls = Cell::new(0usize);
            let outcome =
                build_refit_sample_basis_row_with_poll(admitted, admitted, 0.25, 0.75, &mut || {
                    polls.set(polls.get() + 1);
                    cancel_at == Some(polls.get())
                })
                .expect("valid sample basis row");
            (outcome, polls.get())
        };

        let (complete, healthy_polls) = run(None);
        assert!(matches!(complete, RefitSampleBasisRowRun::Complete { .. }));
        assert!(healthy_polls > 16);
        for cancel_at in [1, healthy_polls / 2, healthy_polls] {
            assert_eq!(
                run(Some(cancel_at)),
                (RefitSampleBasisRowRun::Cancelled, cancel_at)
            );
        }
        with_refit_cx(true, |cx| {
            assert_eq!(
                build_refit_sample_basis_row_with_cx(admitted, admitted, 0.25, 0.75, cx)
                    .expect("pre-cancellation is terminal"),
                RefitSampleBasisRowRun::Cancelled
            );
        });
    }

    // G0/G4: U/V parameter and hostile aggregate count/work/payload refusals
    // precede the first nested basis allocation.
    #[test]
    fn refit_sample_basis_row_refusals_are_preflighted() {
        let knots = open_uniform_knots(4, 2).expect("quadratic knots");
        let admitted = knots.admitted_after_validation();
        let polls = Cell::new(0usize);
        let parameter_error = build_refit_sample_basis_row_with_poll(
            admitted,
            admitted,
            -0.25,
            f64::NAN,
            &mut || {
                polls.set(polls.get() + 1);
                true
            },
        )
        .expect_err("invalid U parameter must win before cancellation");
        assert!(parameter_error.to_string().contains("refit U sample basis"));
        assert_eq!(polls.get(), 0);

        let overflow = preflight_refit_sample_basis_row_counts(usize::MAX, 1, 2, 1)
            .expect_err("unrepresentable tensor size must be refused");
        assert!(overflow.to_string().contains("row size overflows"));

        let work = preflight_refit_sample_basis_row_counts(16_000, 1, 16_000, 1)
            .expect_err("unbounded tensor work must be refused");
        assert!(work.to_string().contains("work units"));

        let retained = preflight_refit_sample_basis_row_counts(6_000, 1, 6_000, 1)
            .expect_err("unbounded simultaneously live payload must be refused");
        assert!(retained.to_string().contains("retains"));
    }

    // G0/G5: the cancellable RHS primitive preserves deterministic row-major
    // accumulation for an exactly known small system.
    #[test]
    fn refit_right_hand_sides_with_cx_match_known_accumulation() {
        let rows = vec![vec![1.0, 0.0, 2.0], vec![0.5, -1.0, 0.0]];
        let targets = [[2.0, 3.0, 4.0], [6.0, 7.0, 8.0]];
        let expected = [
            vec![5.0, -6.0, 4.0],
            vec![6.5, -7.0, 6.0],
            vec![8.0, -8.0, 8.0],
        ];
        assert_eq!(
            assemble_refit_right_hand_sides(&rows, &targets).expect("synchronous right-hand sides"),
            expected
        );
        with_refit_cx(false, |cx| {
            assert_eq!(
                assemble_refit_right_hand_sides_with_cx(&rows, &targets, cx)
                    .expect("cancellable right-hand sides"),
                RefitRightHandSidesRun::Complete {
                    right_hand_sides: expected.clone()
                }
            );
        });
    }

    // G4/G5: validation, fallible output initialization, accumulation, and
    // publication share replayable checkpoints and expose no partial RHS.
    #[test]
    fn refit_right_hand_sides_polling_is_bounded_and_transactional() {
        let rows = vec![vec![1.0; 130]; 130];
        let targets = vec![[1.0, 2.0, 3.0]; 130];
        let run = |cancel_at: Option<usize>| {
            let polls = Cell::new(0usize);
            let outcome = assemble_refit_right_hand_sides_with_poll(&rows, &targets, &mut || {
                polls.set(polls.get() + 1);
                cancel_at == Some(polls.get())
            })
            .expect("valid right-hand sides");
            (outcome, polls.get())
        };

        let (complete, healthy_polls) = run(None);
        assert!(matches!(complete, RefitRightHandSidesRun::Complete { .. }));
        assert!(healthy_polls > 32);
        for cancel_at in [1, healthy_polls / 2, healthy_polls] {
            assert_eq!(
                run(Some(cancel_at)),
                (RefitRightHandSidesRun::Cancelled, cancel_at)
            );
        }
        with_refit_cx(true, |cx| {
            assert_eq!(
                assemble_refit_right_hand_sides_with_cx(&rows, &targets, cx)
                    .expect("pre-cancellation is terminal"),
                RefitRightHandSidesRun::Cancelled
            );
        });
    }

    // G0/G4: constant/count refusals win before polling and malformed borrowed
    // values refuse before output allocation or accumulation.
    #[test]
    fn refit_right_hand_sides_refusals_are_typed_and_preflighted() {
        let polls = Cell::new(0usize);
        let empty = assemble_refit_right_hand_sides_with_poll(&[], &[], &mut || {
            polls.set(polls.get() + 1);
            true
        })
        .expect_err("empty system must be refused");
        assert!(empty.to_string().contains("at least one sample row"));
        assert_eq!(polls.get(), 0);

        let mismatch = preflight_refit_right_hand_sides_counts(2, 1, 3)
            .expect_err("sample/target mismatch must be refused");
        assert!(mismatch.to_string().contains("2 sample rows but 1 targets"));

        let work = preflight_refit_right_hand_sides_counts(1_000_000, 1_000_000, 1_000_000)
            .expect_err("unbounded traversal work must be refused");
        assert!(work.to_string().contains("work units"));

        let retained_controls = REFIT_MAX_ALLOC_BYTES / (3 * size_of::<f64>()) + 1;
        let retained = preflight_refit_right_hand_sides_counts(1, 1, retained_controls)
            .expect_err("unbounded three-vector output payload must be refused");
        assert!(retained.to_string().contains("retains"));

        let malformed = assemble_refit_right_hand_sides_with_poll(
            &[vec![1.0, 2.0], vec![3.0]],
            &[[1.0; 3], [2.0; 3]],
            &mut || false,
        )
        .expect_err("ragged sample matrix must be refused");
        assert!(malformed.to_string().contains("row 1 has length 1"));

        let nonfinite =
            assemble_refit_right_hand_sides_with_poll(&[vec![f64::NAN]], &[[1.0; 3]], &mut || {
                false
            })
            .expect_err("non-finite matrix value must be refused");
        assert!(
            nonfinite
                .to_string()
                .contains("matrix values must be finite")
        );

        let overflowed = assemble_refit_right_hand_sides_with_poll(
            &[vec![f64::MAX]],
            &[[2.0, 1.0, 1.0]],
            &mut || false,
        )
        .expect_err("finite inputs with overflowed accumulation must be refused");
        assert!(
            overflowed
                .to_string()
                .contains("accumulation became non-finite")
        );
    }

    // G0/G5: the cancellable primitive preserves the legacy U-major point
    // evaluation and strict warning threshold for an exactly known sample set.
    #[test]
    fn refit_residual_summary_matches_known_paired_samples() {
        let rows = vec![
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.5, 0.5, 0.0],
            vec![0.0, 0.0, 0.0, 1.0],
        ];
        let targets = [[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 2.0, 2.0]];
        let uvs = [[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]];
        let net = vec![
            vec![[0.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
            vec![[2.0, 0.0, 0.0], [2.0, 2.0, 0.0]],
        ];
        let rms_residual = det::sqrt(5.0 / 3.0);
        let warnings = vec![LocalizedFitResidualWarning {
            uv: uvs[2],
            point: targets[2],
            residual: 2.0,
        }];

        assert_eq!(
            summarize_refit_residuals(&rows, &targets, &uvs, &net, 1.0)
                .expect("synchronous residual summary"),
            (rms_residual, 2.0, warnings.clone())
        );
        with_refit_cx(false, |cx| {
            assert_eq!(
                summarize_refit_residuals_with_cx(&rows, &targets, &uvs, &net, 1.0, cx)
                    .expect("cancellable residual summary"),
                RefitResidualSummaryRun::Complete {
                    rms_residual,
                    max_residual: 2.0,
                    warnings: warnings.clone(),
                }
            );
        });
    }

    // G4/G5: validation, warning reservation, point evaluation, accumulation,
    // and final publication share replayable checkpoints without partial output.
    #[test]
    fn refit_residual_summary_polling_is_bounded_and_transactional() {
        let rows = vec![vec![1.0 / 130.0; 130]; 130];
        let targets = vec![[0.0; 3]; 130];
        let uvs = vec![[0.25, 0.75]; 130];
        let net = vec![vec![[1.0, 2.0, 3.0]; 13]; 10];
        let run = |cancel_at: Option<usize>| {
            let polls = Cell::new(0usize);
            let outcome = summarize_refit_residuals_with_poll(
                &rows,
                &targets,
                &uvs,
                &net,
                f64::MAX,
                &mut || {
                    polls.set(polls.get() + 1);
                    cancel_at == Some(polls.get())
                },
            )
            .expect("valid residual summary");
            (outcome, polls.get())
        };

        let (complete, healthy_polls) = run(None);
        assert!(matches!(complete, RefitResidualSummaryRun::Complete { .. }));
        assert!(healthy_polls > 512);
        for cancel_at in [1, healthy_polls / 2, healthy_polls] {
            assert_eq!(
                run(Some(cancel_at)),
                (RefitResidualSummaryRun::Cancelled, cancel_at)
            );
        }
        with_refit_cx(true, |cx| {
            assert_eq!(
                summarize_refit_residuals_with_cx(&rows, &targets, &uvs, &net, f64::MAX, cx,)
                    .expect("pre-cancellation is terminal"),
                RefitResidualSummaryRun::Cancelled
            );
        });
    }

    // G0/G4: count-only refusals precede polling and all borrowed values are
    // validated before warning storage or point accumulation is published.
    #[test]
    fn refit_residual_summary_refusals_are_preflighted_and_finite() {
        let polls = Cell::new(0usize);
        let threshold =
            summarize_refit_residuals_with_poll(&[], &[], &[], &[], f64::NAN, &mut || {
                polls.set(polls.get() + 1);
                true
            })
            .expect_err("invalid threshold must be refused before polling");
        assert!(threshold.to_string().contains("threshold must be finite"));
        assert_eq!(polls.get(), 0);

        let empty = preflight_refit_residual_summary_counts(0, 0, 0, 0, 0, 0, 0.0)
            .expect_err("empty samples must be refused");
        assert!(empty.to_string().contains("at least one sample row"));

        let mismatch = preflight_refit_residual_summary_counts(2, 1, 2, 1, 1, 1, 0.0)
            .expect_err("sample/target mismatch must be refused");
        assert!(mismatch.to_string().contains("2 sample rows, 1 targets"));

        let size_overflow = preflight_refit_residual_summary_counts(1, 1, 1, 1, usize::MAX, 2, 0.0)
            .expect_err("unrepresentable control-net size must be refused");
        assert!(size_overflow.to_string().contains("overflows usize"));

        let work = preflight_refit_residual_summary_counts(
            1_000_000, 1_000_000, 1_000_000, 1_000_000, 1_000, 1_000, 0.0,
        )
        .expect_err("unbounded traversal work must be refused");
        assert!(work.to_string().contains("work units"));

        let retained_rows = REFIT_MAX_ALLOC_BYTES / size_of::<LocalizedFitResidualWarning>() + 1;
        let retained = preflight_refit_residual_summary_counts(
            retained_rows,
            retained_rows,
            retained_rows,
            1,
            1,
            1,
            0.0,
        )
        .expect_err("unbounded warning payload must be refused");
        assert!(retained.to_string().contains("warnings retain"));

        let ragged_rows = summarize_refit_residuals_with_poll(
            &[vec![1.0, 0.0], vec![1.0]],
            &[[0.0; 3]; 2],
            &[[0.0; 2]; 2],
            &[vec![[0.0; 3], [0.0; 3]]],
            0.0,
            &mut || false,
        )
        .expect_err("ragged sample rows must be refused");
        assert!(ragged_rows.to_string().contains("row 1 has length 1"));

        let ragged_net = summarize_refit_residuals_with_poll(
            &[vec![1.0, 0.0]],
            &[[0.0; 3]],
            &[[0.0; 2]],
            &[vec![[0.0; 3]], vec![]],
            0.0,
            &mut || false,
        )
        .expect_err("ragged control net must be refused");
        assert!(ragged_net.to_string().contains("control-net row 1"));

        let nonfinite_value = summarize_refit_residuals_with_poll(
            &[vec![f64::NAN]],
            &[[0.0; 3]],
            &[[0.0; 2]],
            &[vec![[0.0; 3]]],
            0.0,
            &mut || false,
        )
        .expect_err("non-finite matrix value must be refused");
        assert!(
            nonfinite_value
                .to_string()
                .contains("sample matrix values must be finite")
        );

        let overflowed = summarize_refit_residuals_with_poll(
            &[vec![f64::MAX]],
            &[[0.0; 3]],
            &[[0.0; 2]],
            &[vec![[2.0, 0.0, 0.0]]],
            0.0,
            &mut || false,
        )
        .expect_err("finite inputs with overflowed evaluation must be refused");
        assert!(
            overflowed
                .to_string()
                .contains("point evaluation became non-finite")
        );
    }

    // G0/G5: the Cx-aware probe preserves the legacy U-major cell-center grid
    // and measured absolute-field reduction on an exactly known surface.
    #[test]
    fn refit_field_probe_with_cx_matches_known_cell_centers() {
        let knots = open_uniform_knots(2, 1).expect("linear knots");
        let net = vec![
            vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
        ];
        let weights = vec![vec![1.0; 2]; 2];
        let surface = NurbsSurface::new(knots.clone(), knots, &net, &weights).expect("surface");
        let admitted = surface.admit().expect("admitted surface");
        let field = |point: [f64; 3]| point[0] + point[1];

        assert_eq!(
            measure_refit_field_probe(admitted, &field, 2).expect("synchronous field probe"),
            1.5
        );
        with_refit_cx(false, |cx| {
            assert_eq!(
                measure_refit_field_probe_with_cx(admitted, &field, 2, cx)
                    .expect("cancellable field probe"),
                RefitFieldProbeRun::Complete { max_abs_field: 1.5 }
            );
        });
    }

    // G4/G5: all nested evaluations, pre-field boundaries, reductions, and
    // final scalar publication share one deterministic cancellation callback.
    #[test]
    fn refit_field_probe_polling_is_bounded_and_transactional() {
        let knots = open_uniform_knots(2, 1).expect("linear knots");
        let net = vec![
            vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
        ];
        let weights = vec![vec![1.0; 2]; 2];
        let surface = NurbsSurface::new(knots.clone(), knots, &net, &weights).expect("surface");
        let admitted = surface.admit().expect("admitted surface");
        let run = |cancel_at: Option<usize>| {
            let polls = Cell::new(0usize);
            let field_calls = Cell::new(0usize);
            let field = |point: [f64; 3]| {
                field_calls.set(field_calls.get() + 1);
                point[0] - point[1]
            };
            let outcome = measure_refit_field_probe_with_poll(admitted, &field, 17, &mut || {
                polls.set(polls.get() + 1);
                cancel_at == Some(polls.get())
            })
            .expect("finite field probe");
            (outcome, polls.get(), field_calls.get())
        };

        let (complete, healthy_polls, healthy_field_calls) = run(None);
        assert!(matches!(complete, RefitFieldProbeRun::Complete { .. }));
        assert_eq!(healthy_field_calls, 17 * 17);
        assert!(healthy_polls > healthy_field_calls * 2);
        for cancel_at in [1, healthy_polls / 2, healthy_polls] {
            let (outcome, polls, field_calls) = run(Some(cancel_at));
            assert_eq!(outcome, RefitFieldProbeRun::Cancelled);
            assert_eq!(polls, cancel_at);
            assert!(field_calls <= healthy_field_calls);
        }

        let public_calls = Cell::new(0usize);
        let public_field = |_: [f64; 3]| {
            public_calls.set(public_calls.get() + 1);
            0.0
        };
        with_refit_cx(true, |cx| {
            assert_eq!(
                measure_refit_field_probe_with_cx(admitted, &public_field, 17, cx)
                    .expect("pre-cancellation is terminal"),
                RefitFieldProbeRun::Cancelled
            );
        });
        assert_eq!(public_calls.get(), 0);
    }

    // G0/G4: hostile grid/work/payload counts refuse before polling, and a
    // non-finite opaque result cannot be published as a measured residual.
    #[test]
    fn refit_field_probe_refusals_are_preflighted_and_finite() {
        let empty = preflight_refit_field_probe_envelope(0, 1, 0)
            .expect_err("zero probe resolution must be refused");
        assert!(empty.to_string().contains("must be positive"));

        let overflow = preflight_refit_field_probe_envelope(usize::MAX, 1, 0)
            .expect_err("unrepresentable probe grid must be refused");
        assert!(overflow.to_string().contains("overflows usize"));

        let capped = preflight_refit_field_probe_envelope(2_049, 1, 0)
            .expect_err("grid above the static point cap must be refused");
        assert!(capped.to_string().contains("exceeds static cap"));

        let work = preflight_refit_field_probe_envelope(2, REFIT_MAX_WORK_UNITS, 0)
            .expect_err("repeated nested evaluation work must be aggregate-bounded");
        assert!(work.to_string().contains("work units"));

        let work_overflow = preflight_refit_field_probe_envelope(2, u128::MAX, 0)
            .expect_err("nested evaluation multiplication must be checked");
        assert!(work_overflow.to_string().contains("overflows u128"));

        let retained =
            preflight_refit_field_probe_envelope(1, 0, REFIT_MAX_ALLOC_BYTES as u128 + 1)
                .expect_err("peak nested basis payload must honor the refit cap");
        assert!(retained.to_string().contains("evaluation retains"));

        let knots = open_uniform_knots(2, 1).expect("linear knots");
        let net = vec![vec![[0.0; 3]; 2]; 2];
        let weights = vec![vec![1.0; 2]; 2];
        let surface = NurbsSurface::new(knots.clone(), knots, &net, &weights).expect("surface");
        let error = measure_refit_field_probe_with_poll(
            surface.admit().expect("admitted surface"),
            &|_| f64::NAN,
            1,
            &mut || false,
        )
        .expect_err("non-finite field samples must be refused");
        assert!(error.to_string().contains("non-finite value at probe"));
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

    // G0/G5: the compound Cx entry preserves the fixed 23-sample diagnostic
    // and exact paired-side angle reduction.
    #[test]
    fn seam_g1_with_cx_matches_the_fixed_interior_diagnostic() {
        let knots_u = open_uniform_knots(3, 2).expect("quadratic U knots");
        let knots_v = open_uniform_knots(2, 1).expect("linear V knots");
        let net = vec![
            vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        ];
        let weights = vec![vec![1.0; 2]; 3];
        let surface =
            NurbsSurface::new(knots_u, knots_v, &net, &weights).expect("closed test surface");
        let admitted = surface.admit().expect("admitted test surface");
        let expected = seam_g1_diagnostic_admitted(admitted).expect("fixed seam diagnostic");
        assert_eq!(expected, (2.0, 0));

        with_refit_cx(false, |cx| {
            assert_eq!(
                measure_refit_seam_g1_with_cx(admitted, cx).expect("cancellable seam diagnostic"),
                RefitSeamG1DiagnosticRun::Complete {
                    max_deviation: expected.0,
                    degenerate_samples: expected.1,
                }
            );
        });
    }

    // G4/G5: all 46 nested partial evaluations, per-sample reductions, and
    // final tuple publication share one replayable cancellation callback.
    #[test]
    fn seam_g1_polling_is_transactional_across_nested_partials() {
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
            let outcome = measure_refit_seam_g1_with_poll(admitted, &mut || {
                polls.set(polls.get() + 1);
                cancel_at == Some(polls.get())
            })
            .expect("finite seam diagnostic");
            (outcome, polls.get())
        };

        let (complete, healthy_polls) = run(None);
        assert!(matches!(
            complete,
            RefitSeamG1DiagnosticRun::Complete { .. }
        ));
        assert!(healthy_polls > REFIT_SEAM_G1_SAMPLE_COUNT * 2);
        for cancel_at in [1, healthy_polls / 2, healthy_polls] {
            assert_eq!(
                run(Some(cancel_at)),
                (RefitSeamG1DiagnosticRun::Cancelled, cancel_at)
            );
        }
        with_refit_cx(true, |cx| {
            assert_eq!(
                measure_refit_seam_g1_with_cx(admitted, cx).expect("pre-cancellation is terminal"),
                RefitSeamG1DiagnosticRun::Cancelled
            );
        });
    }

    // G0/G4: aggregate repeated-partials work and peak temporary payload are
    // rejected synchronously before the first nested evaluation or checkpoint.
    #[test]
    fn seam_g1_aggregate_envelope_refuses_hostile_requests() {
        let work = preflight_refit_seam_g1_envelope(REFIT_MAX_WORK_UNITS, 0)
            .expect_err("46 individually plausible calls must still be aggregate-bounded");
        assert!(work.to_string().contains("work units"));

        let work_overflow = preflight_refit_seam_g1_envelope(u128::MAX, 0)
            .expect_err("aggregate repeated-partials multiplication must be checked");
        assert!(work_overflow.to_string().contains("overflows u128"));

        let retained = preflight_refit_seam_g1_envelope(0, REFIT_MAX_ALLOC_BYTES as u128 + 1)
            .expect_err("peak nested temporary payload must honor the refit cap");
        assert!(retained.to_string().contains("partials retain"));
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
        with_refit_cx(false, |cx| {
            assert_eq!(
                measure_refit_seam_g1_with_cx(admitted, cx)
                    .expect("degenerate diagnostic is a complete no-claim"),
                RefitSeamG1DiagnosticRun::Complete {
                    max_deviation: f64::INFINITY,
                    degenerate_samples: REFIT_SEAM_G1_SAMPLE_COUNT,
                }
            );
        });
    }
}
