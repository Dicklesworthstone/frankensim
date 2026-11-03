//! Whole affine paths of a sphere and a six-node quadratic triangle.
//!
//! Each original Lagrange node and the sphere centre move affinely in time.
//! Bernstein convex hulls and the shared bounded global distance query prove
//! clearance over complete time windows. A midpoint query is never itself a
//! clearance proof: its lower distance must cover both rounded reconstruction
//! error and the Hausdorff motion of the relative Bernstein controls.

use super::super::{
    TriangleQueryError,
    sweep::{LinearSweepControls, LinearSweepReport, LinearSweepVerdict},
};
use super::{
    QuadraticQueryControls, QuadraticQueryError, closest_quadratic_triangle_point,
    enclosed_norm_upper, relative_bernstein_controls,
};
use fs_exec::Cx;
use fs_ivl::Interval;
use fs_ivl::expansion::{diff_expansion, expansion_sign, fast_expansion_sum_zeroelim};

/// Separate limits for time subdivision, each distance solve, and total cells.
#[derive(Clone, Copy, Debug)]
pub struct QuadraticSweepControls {
    pub time: LinearSweepControls,
    /// Only the distance tolerance is needed by the sweep. The local solver's
    /// witness tolerance is widened because clearance needs no unique witness.
    pub distance: QuadraticQueryControls,
    /// Aggregate distance cells, including a final solve that exhausts its cap.
    /// Must be positive and at most 2^24; per-query limits still apply.
    pub maximum_distance_cells: usize,
}

impl Default for QuadraticSweepControls {
    fn default() -> Self {
        Self {
            time: LinearSweepControls::default(),
            distance: QuadraticQueryControls::default(),
            maximum_distance_cells: 1 << 20,
        }
    }
}

impl QuadraticSweepControls {
    /// Validate even when a caller's broad phase finds no candidate patches.
    pub fn validate(self) -> Result<(), QuadraticQueryError> {
        validate(self)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadraticSweepReport {
    pub path: LinearSweepReport,
    pub examined_distance_cells: usize,
}

/// Prove nonpenetration throughout the declared six-node affine motion.
///
/// Nodes are corners 0, 1, 2 followed by midside nodes 01, 12, 20. Endpoints
/// are original binary64 coordinates interpreted as exact real input. All
/// admission bounds use fs-ivl directed arithmetic. Relative controls cancel
/// common translation before bounding motion; all six nodes affect the proof.
///
/// `Clear` covers every fraction in [0, 1]. `PossibleContact` identifies the
/// earliest unresolved window, including grazing, and is never an impact-time
/// root. Work exhaustion and cancellation are typed refusals. This proves only
/// the declared coordinate interpolation, not a continuum trajectory, patch
/// injectivity, material validity, or absence of interior solid containment.
pub fn linear_sphere_quadratic_triangle(
    center_start: [f64; 3],
    center_end: [f64; 3],
    radius: f64,
    nodes_start: [[f64; 3]; 6],
    nodes_end: [[f64; 3]; 6],
    controls: QuadraticSweepControls,
    cx: &Cx<'_>,
) -> Result<QuadraticSweepReport, QuadraticQueryError> {
    poll(cx)?;
    controls.validate()?;
    if !radius.is_finite()
        || radius <= 0.0
        || center_start
            .into_iter()
            .chain(center_end)
            .chain(nodes_start.into_iter().flatten())
            .chain(nodes_end.into_iter().flatten())
            .any(|value| !value.is_finite())
    {
        return Err(QuadraticQueryError::Invalid(
            "finite affine geometry and positive sphere radius required",
        ));
    }
    let start = relative_bernstein_controls(nodes_start, center_start)?;
    let end = relative_bernstein_controls(nodes_end, center_end)?;
    let stationary = center_start == center_end && nodes_start == nodes_end;
    let mut pending = vec![(0.0_f64, 1.0_f64)];
    let mut examined_windows = 0;
    let mut examined_distance_cells = 0;
    let mut stationary_distance = None;
    let mut minimum = f64::INFINITY;
    while let Some((lo, hi)) = pending.pop() {
        poll(cx)?;
        if examined_windows == controls.time.maximum_windows {
            return Err(QuadraticQueryError::Geometry(
                TriangleQueryError::WorkBudgetExhausted {
                    maximum_visits: controls.time.maximum_windows,
                },
            ));
        }
        examined_windows += 1;
        let mid = f64::midpoint(lo, hi);
        let sample_fraction = if stationary { 0.0 } else { mid };
        let sample = relative_sample(
            nodes_start,
            center_start,
            nodes_end,
            center_end,
            sample_fraction,
        )?;
        let direction = sample[0].map(|value| -value);
        let mut clearance = plane_clearance(start, end, radius, direction, lo, hi);
        if clearance.is_none()
            && lo == 0.0
            && hi == 1.0
            && exact_axis_separation(center_start, center_end, radius, nodes_start, nodes_end)
        {
            clearance = Some(0.0);
        }
        if clearance.is_none() {
            let closest = if let Some(cached) = stationary_distance {
                cached
            } else {
                let remaining = controls.maximum_distance_cells - examined_distance_cells;
                if remaining == 0 {
                    return Err(QuadraticQueryError::CellBudgetExhausted {
                        maximum_cells: controls.maximum_distance_cells,
                    });
                }
                let mut distance_controls = controls.distance;
                distance_controls.maximum_cells = distance_controls.maximum_cells.min(remaining);
                distance_controls.maximum_active_cells = distance_controls
                    .maximum_active_cells
                    .min(distance_controls.maximum_cells);
                distance_controls.witness_tolerance = f64::MAX;
                let result =
                    closest_quadratic_triangle_point([0.0; 3], sample, distance_controls, cx)?;
                examined_distance_cells += result.examined_cells;
                // A static unresolved touch must not rerun the global solve
                // at every time subdivision. The surrogate is identical.
                if stationary {
                    stationary_distance = Some(result);
                }
                result
            };
            let direction = closest.point.map(|value| -value);
            clearance = plane_clearance(start, end, radius, direction, lo, hi);
            if clearance.is_none() {
                let sign = if direction
                    .into_iter()
                    .zip(closest.unit_normal)
                    .map(|(a, b)| a * b)
                    .sum::<f64>()
                    < 0.0
                {
                    -1.0
                } else {
                    1.0
                };
                clearance = plane_clearance(
                    start,
                    end,
                    radius,
                    closest.unit_normal.map(|value| sign * value),
                    lo,
                    hi,
                );
            }
            if clearance.is_none() {
                let surrogate = relative_bernstein_controls(sample, [0.0; 3])?;
                let mut reconstruction = 0.0_f64;
                let mut motion = 0.0_f64;
                let half_width = Interval::point(
                    (Interval::point(mid) - Interval::point(lo))
                        .hi()
                        .max((Interval::point(hi) - Interval::point(mid)).hi()),
                );
                for i in 0..6 {
                    reconstruction =
                        reconstruction.max(enclosed_norm_upper(core::array::from_fn(|axis| {
                            affine_value(start[i][axis], end[i][axis], sample_fraction)
                                - surrogate[i][axis]
                        }))?);
                    if !stationary {
                        // Scale the control velocity before its norm to avoid
                        // overflowing a speed whose window displacement fits.
                        motion = motion.max(enclosed_norm_upper(core::array::from_fn(|axis| {
                            (end[i][axis] - start[i][axis]) * half_width
                        }))?);
                    }
                }
                let lower = (Interval::point(closest.distance_lower)
                    - Interval::point(reconstruction)
                    - Interval::point(motion)
                    - Interval::point(radius))
                .lo();
                if lower.is_finite() && lower >= 0.0 {
                    clearance = Some(lower);
                }
            }
        }
        if let Some(gap) = clearance {
            minimum = minimum.min(gap);
            continue;
        }
        if hi - lo <= controls.time.fraction_tolerance || mid == lo || mid == hi {
            poll(cx)?;
            return Ok(QuadraticSweepReport {
                path: LinearSweepReport {
                    verdict: LinearSweepVerdict::PossibleContact {
                        fraction_start: lo,
                        fraction_end: hi,
                    },
                    examined_windows,
                },
                examined_distance_cells,
            });
        }
        // Visit earlier time first, preserving the original sweep semantics.
        pending.push((mid, hi));
        pending.push((lo, mid));
    }
    poll(cx)?;
    Ok(QuadraticSweepReport {
        path: LinearSweepReport {
            verdict: LinearSweepVerdict::Clear {
                minimum_clearance: minimum,
            },
            examined_windows,
        },
        examined_distance_cells,
    })
}

fn validate(controls: QuadraticSweepControls) -> Result<(), QuadraticQueryError> {
    super::validate_controls(controls.distance)?;
    if controls.time.maximum_windows == 0
        || controls.time.maximum_windows > 1 << 20
        || !controls.time.fraction_tolerance.is_finite()
        || controls.time.fraction_tolerance <= 0.0
        || controls.time.fraction_tolerance >= 1.0
        || controls.maximum_distance_cells == 0
        || controls.maximum_distance_cells > 1 << 24
    {
        return Err(QuadraticQueryError::Invalid(
            "positive bounded time/distance work and fraction tolerance in (0,1) required",
        ));
    }
    Ok(())
}

fn poll(cx: &Cx<'_>) -> Result<(), QuadraticQueryError> {
    super::super::poll(cx).map_err(QuadraticQueryError::Geometry)
}

fn affine_value(start: Interval, end: Interval, fraction: f64) -> Interval {
    if fraction == 0.0 {
        start
    } else if fraction == 1.0 {
        end
    } else {
        let fraction = Interval::point(fraction);
        (Interval::point(1.0) - fraction) * start + fraction * end
    }
}

// A rounded surrogate is only a query input. Its entire Bernstein-control
// error is enclosed separately before a sampled distance can admit a window.
fn relative_sample(
    nodes_start: [[f64; 3]; 6],
    center_start: [f64; 3],
    nodes_end: [[f64; 3]; 6],
    center_end: [f64; 3],
    fraction: f64,
) -> Result<[[f64; 3]; 6], QuadraticQueryError> {
    let mut result = [[0.0; 3]; 6];
    for i in 0..6 {
        for axis in 0..3 {
            let value = affine_value(
                Interval::point(nodes_start[i][axis]) - Interval::point(center_start[axis]),
                Interval::point(nodes_end[i][axis]) - Interval::point(center_end[axis]),
                fraction,
            );
            if !value.lo().is_finite() || !value.hi().is_finite() {
                return Err(QuadraticQueryError::NumericalRange);
            }
            result[i][axis] = f64::midpoint(value.lo(), value.hi());
        }
    }
    Ok(result)
}

fn plane_clearance(
    start: [[Interval; 3]; 6],
    end: [[Interval; 3]; 6],
    radius: f64,
    direction: [f64; 3],
    lo: f64,
    hi: f64,
) -> Option<f64> {
    let scale = direction.into_iter().map(f64::abs).fold(0.0, f64::max);
    if scale == 0.0 || !scale.is_finite() {
        return None;
    }
    let normal = direction.map(|value| Interval::point(value / scale));
    let length = normal
        .into_iter()
        .fold(Interval::point(0.0), |sum, value| sum + value * value)
        .sqrt();
    if length.lo() <= 0.0 || !length.hi().is_finite() {
        return None;
    }
    let support = Interval::point(radius) * length;
    let project = |point: [Interval; 3]| {
        point
            .into_iter()
            .zip(normal)
            .fold(Interval::point(0.0), |sum, (value, n)| sum - value * n)
    };
    let mut lower = f64::INFINITY;
    for i in 0..6 {
        let a = project(start[i]);
        let b = project(end[i]);
        let projection = affine_value(a, b, lo).hull(affine_value(a, b, hi));
        lower = lower.min(((projection - support) / length).lo());
    }
    (lower.is_finite() && lower >= 0.0).then_some(lower)
}

// Exact original-endpoint support inequalities preserve axis-aligned touching.
// For an edge control B = 2 M - (A + D)/2, test twice the support gap:
// sign * [2(C-M) + (A-M) + (D-M)] - r - r >= 0. Expansion additions avoid
// rounding the control or halving a subnormal coordinate. This only proves the
// complete original path; rounded subwindow endpoints are never substituted.
fn exact_axis_separation(
    center_start: [f64; 3],
    center_end: [f64; 3],
    radius: f64,
    nodes_start: [[f64; 3]; 6],
    nodes_end: [[f64; 3]; 6],
) -> bool {
    if radius > f64::MAX / 64.0
        || center_start
            .into_iter()
            .chain(center_end)
            .chain(nodes_start.into_iter().flatten())
            .chain(nodes_end.into_iter().flatten())
            .any(|value| value.abs() > f64::MAX / 64.0)
    {
        return false;
    }
    let nonnegative = |mut gap: Vec<f64>, sign: f64, doubled: bool| {
        if sign < 0.0 {
            for value in &mut gap {
                *value = -*value;
            }
        }
        gap = fast_expansion_sum_zeroelim(&gap, &[-radius]);
        if doubled {
            gap = fast_expansion_sum_zeroelim(&gap, &[-radius]);
        }
        expansion_sign(&gap) >= 0
    };
    for axis in 0..3 {
        for sign in [-1.0, 1.0] {
            let separated = [(center_start, nodes_start), (center_end, nodes_end)]
                .into_iter()
                .all(|(center, nodes)| {
                    (0..3).all(|i| {
                        nonnegative(diff_expansion(center[axis], nodes[i][axis]), sign, false)
                    }) && [(0, 1, 3), (1, 2, 4), (2, 0, 5)]
                        .into_iter()
                        .all(|(a, b, m)| {
                            let cm = diff_expansion(center[axis], nodes[m][axis]);
                            let gap = fast_expansion_sum_zeroelim(&cm, &cm);
                            let gap = fast_expansion_sum_zeroelim(
                                &gap,
                                &diff_expansion(nodes[a][axis], nodes[m][axis]),
                            );
                            let gap = fast_expansion_sum_zeroelim(
                                &gap,
                                &diff_expansion(nodes[b][axis], nodes[m][axis]),
                            );
                            nonnegative(gap, sign, true)
                        })
                });
            if separated {
                return true;
            }
        }
    }
    false
}

#[cfg(test)]
mod tests;
