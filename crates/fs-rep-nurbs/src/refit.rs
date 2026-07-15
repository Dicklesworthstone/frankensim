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
use crate::basis::KnotVector;
use crate::closest::norm3;
use crate::surface::NurbsSurface;
use core::mem::size_of;
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

fn refit_structure_error(what: impl Into<String>) -> NurbsError {
    NurbsError::Structure { what: what.into() }
}

fn refit_allocation_error(what: &'static str) -> NurbsError {
    NurbsError::Domain {
        what: format!("{what} allocation was refused"),
    }
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
    let assembly_work = (sample_points as u128)
        .saturating_mul(active_basis as u128)
        .saturating_mul(control_points as u128);
    let factor_work = (control_points as u128)
        .saturating_mul(control_points as u128)
        .saturating_mul(control_points as u128);
    let rhs_and_report_work = (sample_points as u128)
        .saturating_mul(control_points as u128)
        .saturating_mul(6);
    let triangular_solve_work = (control_points as u128)
        .saturating_mul(control_points as u128)
        .saturating_mul(3);
    let projection_evaluations = (sample_points as u128).saturating_mul(42);
    // Every retained probe calls the public surface evaluator. That path
    // revalidates the whole mutable control net and knot vectors, builds two
    // Cox-de Boor triangles, and accumulates the active tensor basis. Charge
    // those stages rather than pretending one probe is one work unit. The
    // arbitrary closure's own cost remains outside this legacy static model.
    let per_probe_work = (control_points as u128)
        .saturating_add(knots_per_axis as u128)
        .saturating_add((basis_triangle as u128).saturating_mul(2))
        .saturating_add((active_basis as u128).saturating_mul(8))
        .saturating_add(16);
    let probe_work = (probe_points as u128).saturating_mul(per_probe_work);
    let total_work = assembly_work
        .saturating_add(factor_work)
        .saturating_add(rhs_and_report_work)
        .saturating_add(triangular_solve_work)
        .saturating_add(projection_evaluations)
        .saturating_add(probe_work);
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

/// Bisect the implicit field along `center + r·dir` for a sign crossing.
fn project_radial(
    field: &dyn Fn([f64; 3]) -> f64,
    center: [f64; 3],
    dir: [f64; 3],
    r_max: f64,
) -> Result<f64, NurbsError> {
    let at = |r: f64| -> Result<f64, NurbsError> {
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
        Ok(value)
    };
    let (mut lo, mut hi) = (0.0f64, r_max);
    if at(lo)? >= 0.0 || at(hi)? <= 0.0 {
        return Err(NurbsError::Structure {
            what: format!(
                "radial bracket failed along {dir:?}: refit v1 needs a star-shaped \
                 domain around the given center (field(center) < 0 < field(center + r_max·dir))"
            ),
        });
    }
    for _ in 0..40 {
        let mid = f64::midpoint(lo, hi);
        if at(mid)? < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(f64::midpoint(lo, hi))
}

/// Dense symmetric-positive-definite Cholesky factorization in place. The
/// factor is shared across all three coordinate right-hand sides.
fn cholesky_factor(a: &mut [Vec<f64>]) -> Result<(), NurbsError> {
    let n = a.len();
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i][j];
            let (ri, rj) = (&a[i], &a[j]);
            for (x, y) in ri[..j].iter().zip(&rj[..j]) {
                sum -= x * y;
            }
            if i == j {
                if !sum.is_finite() || sum <= 0.0 {
                    return Err(NurbsError::Structure {
                        what: "normal equations not SPD (raise lambda or sample count)".to_string(),
                    });
                }
                a[i][i] = det::sqrt(sum);
            } else {
                a[i][j] = sum / a[j][j];
                if !a[i][j].is_finite() {
                    return Err(refit_structure_error(
                        "normal-equation factorization became non-finite",
                    ));
                }
            }
        }
    }
    Ok(())
}

/// Solve one right-hand side using a factor produced by
/// [`cholesky_factor`].
fn cholesky_solve_factored(a: &[Vec<f64>], b: &mut [f64]) -> Result<(), NurbsError> {
    let n = b.len();
    for i in 0..n {
        let mut sum = b[i];
        for k in 0..i {
            sum -= a[i][k] * b[k];
        }
        b[i] = sum / a[i][i];
        if !b[i].is_finite() {
            return Err(refit_structure_error(
                "normal-equation forward solve became non-finite",
            ));
        }
    }
    for i in (0..n).rev() {
        let mut sum = b[i];
        for k in (i + 1)..n {
            sum -= a[k][i] * b[k];
        }
        b[i] = sum / a[i][i];
        if !b[i].is_finite() {
            return Err(refit_structure_error(
                "normal-equation back solve became non-finite",
            ));
        }
    }
    Ok(())
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
fn basis_row(kv: &KnotVector<f64>, n: usize, t: f64) -> Result<Vec<f64>, NurbsError> {
    let (span, vals) = kv.basis(t)?;
    let mut row = try_filled_vec(n, 0.0f64, "refit dense basis row")?;
    let p = kv.degree;
    for (c, &b) in vals.iter().enumerate() {
        row[span - p + c] = b;
    }
    Ok(row)
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
/// estimate, not an outward-rounded enclosure. Returns (L_u, L_v).
fn lipschitz_estimate(surface: &NurbsSurface<f64>) -> (f64, f64) {
    let p_u = surface.knots_u.degree;
    let p_v = surface.knots_v.degree;
    let ku = &surface.knots_u.knots;
    let kv = &surface.knots_v.knots;
    let cart = |h: &[f64; 4]| [h[0] / h[3], h[1] / h[3], h[2] / h[3]];
    let dist = |a: [f64; 3], b: [f64; 3]| -> f64 { norm3([a[0] - b[0], a[1] - b[1], a[2] - b[2]]) };
    let rows = surface.cpw.len();
    let cols = surface.cpw[0].len();
    let mut lu = 0.0f64;
    let mut lv = 0.0f64;
    for i in 0..rows {
        for j in 0..cols {
            let c = cart(&surface.cpw[i][j]);
            if i + 1 < rows {
                let dc = dist(cart(&surface.cpw[i + 1][j]), c);
                let span = ku[i + p_u + 1] - ku[i + 1];
                if span > 0.0 {
                    #[allow(clippy::cast_precision_loss)]
                    let coef = p_u as f64 * dc / span;
                    lu = lu.max(coef);
                }
            }
            if j + 1 < cols {
                let dc = dist(cart(&surface.cpw[i][j + 1]), c);
                let span = kv[j + p_v + 1] - kv[j + 1];
                if span > 0.0 {
                    #[allow(clippy::cast_precision_loss)]
                    let coef = p_v as f64 * dc / span;
                    lv = lv.max(coef);
                }
            }
        }
    }
    (lu, lv)
}

/// Fit one scalar/vector LSQ system: `(BᵀB + λ LᵀL) c = Bᵀy` where `L`
/// is the discrete control-lattice Laplacian (thin-plate proxy).
#[allow(clippy::needless_range_loop)]
fn assemble_normal(
    rows_b: &[Vec<f64>],
    nu: usize,
    nv: usize,
    lambda: f64,
) -> Result<Vec<Vec<f64>>, NurbsError> {
    let n = nu
        .checked_mul(nv)
        .ok_or_else(|| refit_structure_error("refit normal-matrix dimension overflow"))?;
    let mut a = try_filled_matrix(n, n, 0.0f64, "refit normal matrix")?;
    for row in rows_b {
        for i in 0..n {
            if row[i] == 0.0 {
                continue;
            }
            for j in 0..n {
                if row[j] != 0.0 {
                    a[i][j] += row[i] * row[j];
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
                    a[p][q] += lambda * wp * wq;
                }
            }
        }
    }
    Ok(a)
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
            let bu = basis_row(&ku, nu, u)?;
            let bv = basis_row(&kv, nv, v)?;
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
    let mut factor = assemble_normal(&rows_b, nu, nv, config.lambda)?;
    cholesky_factor(&mut factor)?;
    for axis in 0..3 {
        let mut rhs = try_filled_vec(control_points, 0.0f64, "refit right-hand side")?;
        for (row, t) in rows_b.iter().zip(&targets) {
            for (k, &w) in row.iter().enumerate() {
                if w != 0.0 {
                    rhs[k] += w * t[axis];
                }
            }
        }
        cholesky_solve_factored(&factor, &mut rhs)?;
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
            let p = surface.eval(u, v)?;
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
    let (lip_u, lip_v) = lipschitz_estimate(&surface);
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
    let (seam_g1, seam_g1_degenerate_samples) = seam_g1_diagnostic(&surface)?;
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

fn seam_g1_diagnostic(surface: &NurbsSurface<f64>) -> Result<(f64, usize), NurbsError> {
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
    use std::cell::Cell;

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
        assert_eq!(
            calls.get(),
            0,
            "all shape/cap refusals must precede field evaluation"
        );
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
        let (lu, _lv) = lipschitz_estimate(&surface);
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

    #[test]
    fn degenerate_seam_tangents_are_explicitly_no_claim() {
        let line = open_uniform_knots(2, 1).expect("linear knots");
        let points = vec![vec![[0.0; 3]; 2]; 2];
        let weights = vec![vec![1.0; 2]; 2];
        let surface =
            NurbsSurface::new(line.clone(), line, &points, &weights).expect("degenerate surface");
        let (g1, degenerate) = seam_g1_diagnostic(&surface).expect("bounded diagnostic");
        assert!(g1.is_infinite(), "undefined tangent direction is no-claim");
        assert_eq!(degenerate, 23, "every retained seam sample is named");
    }
}
