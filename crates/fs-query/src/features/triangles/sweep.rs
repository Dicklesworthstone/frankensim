//! Complete straight sphere/triangle paths, using the original interval owner.
//!
//! A fixed separating plane is checked at BOTH endpoints of a subwindow.
//! Every vertex and the sphere centre are affine in time, so this proves
//! separation throughout that window. Closest-point samples only suggest plane
//! directions; they never admit a path. Windows without a proof are subdivided
//! and eventually return PossibleContact, including unresolved grazing. This
//! certifies the declared linear coordinate interpolation, not a curved flow.

use super::{TriangleQueryError, closest_triangle_point};
use fs_exec::Cx;
use fs_ivl::Interval;
use fs_ivl::expansion::{diff_expansion, expansion_sign, fast_expansion_sum_zeroelim};

/// Maximum examined windows and unresolved fractional-time width.
#[derive(Clone, Copy, Debug)]
pub struct LinearSweepControls {
    pub maximum_windows: usize,
    pub fraction_tolerance: f64,
}
impl Default for LinearSweepControls {
    fn default() -> Self {
        Self {
            maximum_windows: 4096,
            fraction_tolerance: 1e-8,
        }
    }
}

/// A complete proof or the first window that cannot be proved clear.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum LinearSweepVerdict {
    /// Every fraction in [0,1] is nonpenetrating. Exact axis-aligned touch may
    /// have zero clearance; a positive value is a rigorous separation bound.
    Clear { minimum_clearance: f64 },
    /// This is a possible encounter, never a claimed impact or a time root.
    PossibleContact {
        fraction_start: f64,
        fraction_end: f64,
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearSweepReport {
    pub verdict: LinearSweepVerdict,
    pub examined_windows: usize,
}

/// Conservative upper clearance over an affine sphere/deforming-triangle path.
///
/// Each endpoint's barycentric weights are nonnegative and normalized by their
/// enclosed exact real sum. Interpolating those weights produces a convex point
/// of the SAME moving triangle at every time; the proof never substitutes
/// a rounded reconstructed point that might lie outside it. Relative coordinates
/// use triangle edges before accumulation, retaining small separations far from
/// the world origin. All arithmetic is directed through the original fs-ivl owner.
///
/// The relative point path is quadratic. Its Bernstein controls are the two
/// endpoint relative vectors and (C0+C1-sum(w0*T1)-sum(w1*T0))/2. Their maximum
/// enclosed norm bounds the complete path by convex-hull containment and norm
/// convexity. The middle control is evaluated as half the endpoint-vector sum
/// plus half sum(delta(w_i)*delta(edge_i)), removing common translation and
/// avoiding cancellation between large, opposite sliding displacements.
/// On a stationary plane under steady sliding all three controls
/// equal the normal separation, avoiding a false timestep restriction from
/// holding the starting closest point fixed while the sphere moves tangentially.
/// Distance to the closest triangle point cannot be larger than distance to
/// this convex point. The returned upper clearance subtracts the ORIGINAL
/// sphere radius. This establishes nearness only: nonpenetration still needs
/// the separate whole-path lower-bound admission. Work and storage are constant.
#[allow(clippy::too_many_arguments)]
pub fn linear_sphere_triangle_nearness(
    center_start: [f64; 3],
    center_end: [f64; 3],
    radius: f64,
    triangle_start: [[f64; 3]; 3],
    triangle_end: [[f64; 3]; 3],
    weights_start: [f64; 3],
    weights_end: [f64; 3],
    cx: &Cx<'_>,
) -> Result<f64, TriangleQueryError> {
    super::poll(cx)?;
    if !radius.is_finite()
        || radius <= 0.0
        || center_start
            .into_iter()
            .chain(center_end)
            .chain(triangle_start.into_iter().flatten())
            .chain(triangle_end.into_iter().flatten())
            .any(|x| !x.is_finite())
    {
        return Err(TriangleQueryError::Invalid(
            "finite affine geometry and positive sphere radius required",
        ));
    }
    let normalized_start = enclosed_barycentric(weights_start)?;
    let normalized_end = enclosed_barycentric(weights_end)?;
    let relative = |center: [f64; 3], triangle: [[f64; 3]; 3], normalized: [Interval; 3]| {
        core::array::from_fn(|axis| {
            let anchor = Interval::point(triangle[0][axis]);
            Interval::point(center[axis])
                - anchor
                - normalized[1] * (Interval::point(triangle[1][axis]) - anchor)
                - normalized[2] * (Interval::point(triangle[2][axis]) - anchor)
        })
    };
    let start_relative: [Interval; 3] = relative(center_start, triangle_start, normalized_start);
    let end_relative: [Interval; 3] = relative(center_end, triangle_end, normalized_end);
    let half = Interval::point(0.5);
    let middle_relative = core::array::from_fn(|axis| {
        let mut middle = half * start_relative[axis] + half * end_relative[axis];
        // Identical original endpoint data prove the corresponding exact
        // change is zero, without subtracting two dependent enclosures.
        if weights_start != weights_end && triangle_start != triangle_end {
            for node in 1..3 {
                let weight_change = normalized_end[node] - normalized_start[node];
                let start_edge = Interval::point(triangle_start[node][axis])
                    - Interval::point(triangle_start[0][axis]);
                let end_edge = Interval::point(triangle_end[node][axis])
                    - Interval::point(triangle_end[0][axis]);
                middle = middle + half * weight_change * (end_edge - start_edge);
            }
        }
        middle
    });
    let start = enclosed_norm_upper(start_relative)?;
    let end = enclosed_norm_upper(end_relative)?;
    let middle = enclosed_norm_upper(middle_relative)?;
    let upper = (Interval::point(start.max(middle).max(end)) - Interval::point(radius)).hi();
    super::poll(cx)?;
    super::finite(upper)
}

fn enclosed_barycentric(weights: [f64; 3]) -> Result<[Interval; 3], TriangleQueryError> {
    let sum: f64 = weights.iter().sum();
    if weights
        .iter()
        .any(|w| !w.is_finite() || !(0.0..=1.0).contains(w))
        || (sum - 1.0).abs() > 64.0 * f64::EPSILON
    {
        return Err(TriangleQueryError::Invalid(
            "convex near-unit barycentric weights required",
        ));
    }
    let denominator = weights
        .into_iter()
        .fold(Interval::point(0.0), |sum, w| sum + Interval::point(w));
    if !denominator.lo().is_finite() || denominator.lo() <= 0.0 || !denominator.hi().is_finite() {
        return Err(TriangleQueryError::NumericalRange);
    }
    Ok(weights.map(|w| Interval::point(w) / denominator))
}

/// Upper clearance over the complete affine path of two sphere centres.
///
/// Convexity bounds the centre distance by its endpoint maximum; each ORIGINAL
/// radius is subtracted separately with directed rounding. This is a nearness
/// upper bound, not the independent nonpenetration lower-bound certificate.
pub fn linear_spheres_nearness(
    first_start: [f64; 3],
    first_end: [f64; 3],
    first_radius: f64,
    second_start: [f64; 3],
    second_end: [f64; 3],
    second_radius: f64,
    cx: &Cx<'_>,
) -> Result<f64, TriangleQueryError> {
    super::poll(cx)?;
    if [first_radius, second_radius]
        .iter()
        .any(|r| !r.is_finite() || *r <= 0.0)
        || first_start
            .into_iter()
            .chain(first_end)
            .chain(second_start)
            .chain(second_end)
            .any(|x| !x.is_finite())
    {
        return Err(TriangleQueryError::Invalid(
            "finite affine centres and positive original radii required",
        ));
    }
    let relative = |first: [f64; 3], second: [f64; 3]| {
        core::array::from_fn(|axis| Interval::point(first[axis]) - Interval::point(second[axis]))
    };
    let start = enclosed_norm_upper(relative(first_start, second_start))?;
    let end = enclosed_norm_upper(relative(first_end, second_end))?;
    let upper = (Interval::point(start.max(end))
        - Interval::point(first_radius)
        - Interval::point(second_radius))
    .hi();
    super::poll(cx)?;
    super::finite(upper)
}

fn enclosed_norm_upper(vector: [Interval; 3]) -> Result<f64, TriangleQueryError> {
    if vector
        .iter()
        .any(|v| !v.lo().is_finite() || !v.hi().is_finite())
    {
        return Err(TriangleQueryError::NumericalRange);
    }
    let scale = vector.iter().map(|v| v.abs().hi()).fold(0.0, f64::max);
    if scale == 0.0 {
        return Ok(0.0);
    }
    let squared = vector.into_iter().fold(Interval::point(0.0), |sum, v| {
        let x = (v / Interval::point(scale)).abs();
        sum + x * x
    });
    if !squared.hi().is_finite() || squared.hi() < 0.0 {
        return Err(TriangleQueryError::NumericalRange);
    }
    super::finite((squared.sqrt() * Interval::point(scale)).hi())
}

/// Certify nonpenetration of a sphere and one linearly deforming triangle.
///
/// Endpoints are original floating-point coordinates interpreted as exact real
/// input. All admission arithmetic uses fs-ivl directed enclosures. An exact
/// expansion fallback admits axis-aligned touching/receding paths without
/// turning an epsilon into permitted penetration. No direction normalization
/// accuracy is assumed by the proof: sphere support uses an enclosed |n|.
pub fn linear_sphere_triangle(
    center_start: [f64; 3],
    center_end: [f64; 3],
    radius: f64,
    triangle_start: [[f64; 3]; 3],
    triangle_end: [[f64; 3]; 3],
    controls: LinearSweepControls,
    cx: &Cx<'_>,
) -> Result<LinearSweepReport, TriangleQueryError> {
    linear_sphere_triangle_with_radii(
        center_start,
        center_end,
        [radius, 0.0],
        triangle_start,
        triangle_end,
        controls,
        cx,
    )
}

// The interval fast path uses an enclosing support radius. The exact cold
// path retains each ORIGINAL radius, so upward rounding of a radius sum cannot
// turn an actual touching sphere pair into an artificial overlap.
fn linear_sphere_triangle_with_radii(
    center_start: [f64; 3],
    center_end: [f64; 3],
    radii: [f64; 2],
    triangle_start: [[f64; 3]; 3],
    triangle_end: [[f64; 3]; 3],
    controls: LinearSweepControls,
    cx: &Cx<'_>,
) -> Result<LinearSweepReport, TriangleQueryError> {
    super::poll(cx)?;
    validate(controls)?;
    if !radii[0].is_finite()
        || radii[0] <= 0.0
        || !radii[1].is_finite()
        || radii[1] < 0.0
        || center_start
            .into_iter()
            .chain(center_end)
            .chain(triangle_start.into_iter().flatten())
            .chain(triangle_end.into_iter().flatten())
            .any(|x| !x.is_finite())
    {
        return Err(TriangleQueryError::Invalid(
            "finite coordinates and positive sphere radius required",
        ));
    }
    let radius = if radii[1] == 0.0 {
        radii[0]
    } else {
        (Interval::point(radii[0]) + Interval::point(radii[1])).hi()
    };
    if !radius.is_finite() {
        return Err(TriangleQueryError::NumericalRange);
    }
    let mut pending = vec![(0.0_f64, 1.0_f64)];
    let mut examined = 0;
    let mut minimum = f64::INFINITY;
    while let Some((lo, hi)) = pending.pop() {
        super::poll(cx)?;
        if examined == controls.maximum_windows {
            return Err(TriangleQueryError::WorkBudgetExhausted {
                maximum_visits: controls.maximum_windows,
            });
        }
        examined += 1;
        let mid = f64::midpoint(lo, hi);
        let center = lerp(center_start, center_end, mid);
        let triangle = core::array::from_fn(|i| lerp(triangle_start[i], triangle_end[i], mid));
        // A vertex direction also works for a degenerate triangle representing
        // one moving point (the shared sphere/sphere path owner below).
        let direction = core::array::from_fn(|i| center[i] - triangle[0][i]);
        let mut clearance = plane_clearance(
            center_start,
            center_end,
            radius,
            triangle_start,
            triangle_end,
            direction,
            lo,
            hi,
        );
        // Approximate geometry only selects directions. Even a poorly rounded
        // or nearly collapsed sample cannot make the interval proof unsound.
        if clearance.is_none()
            && let Ok(closest) = closest_triangle_point(center, triangle)
        {
            let direction = core::array::from_fn(|i| center[i] - closest.point[i]);
            clearance = plane_clearance(
                center_start,
                center_end,
                radius,
                triangle_start,
                triangle_end,
                direction,
                lo,
                hi,
            );
            if clearance.is_none() {
                let sign = if dot(direction, closest.unit_normal) < 0.0 {
                    -1.0
                } else {
                    1.0
                };
                clearance = plane_clearance(
                    center_start,
                    center_end,
                    radius,
                    triangle_start,
                    triangle_end,
                    closest.unit_normal.map(|x| sign * x),
                    lo,
                    hi,
                );
            }
        }
        if let Some(gap) = clearance {
            minimum = minimum.min(gap);
            continue;
        }
        // Exact arithmetic is a cold equality fallback after the inexpensive
        // interval proof. Its whole-path result remains valid for this window.
        if lo == 0.0
            && hi == 1.0
            && exact_axis_separation(
                center_start,
                center_end,
                radii,
                triangle_start,
                triangle_end,
            )
        {
            return Ok(LinearSweepReport {
                verdict: LinearSweepVerdict::Clear {
                    minimum_clearance: 0.0,
                },
                examined_windows: examined,
            });
        }
        if hi - lo <= controls.fraction_tolerance || !(mid > lo && mid < hi) {
            return Ok(LinearSweepReport {
                verdict: LinearSweepVerdict::PossibleContact {
                    fraction_start: lo,
                    fraction_end: hi,
                },
                examined_windows: examined,
            });
        }
        // First unresolved time is deterministic; later intervals are private.
        pending.push((mid, hi));
        pending.push((lo, mid));
    }
    Ok(LinearSweepReport {
        verdict: LinearSweepVerdict::Clear {
            minimum_clearance: minimum,
        },
        examined_windows: examined,
    })
}

/// Prove separation of two declared straight sphere-centre trajectories.
/// The fast radius sum is rounded upward; exact axis-aligned separation uses
/// both ORIGINAL radii, admitting real touching/receding paths without allowing
/// penetration. This uses the same separating-plane authority as sphere/triangle
/// paths, with a point as the second convex set.
#[allow(clippy::too_many_arguments)]
pub fn linear_spheres(
    first_start: [f64; 3],
    first_end: [f64; 3],
    first_radius: f64,
    second_start: [f64; 3],
    second_end: [f64; 3],
    second_radius: f64,
    controls: LinearSweepControls,
    cx: &Cx<'_>,
) -> Result<LinearSweepReport, TriangleQueryError> {
    if !first_radius.is_finite()
        || first_radius <= 0.0
        || !second_radius.is_finite()
        || second_radius <= 0.0
    {
        return Err(TriangleQueryError::Invalid(
            "two positive finite sphere radii required",
        ));
    }
    linear_sphere_triangle_with_radii(
        first_start,
        first_end,
        [first_radius, second_radius],
        [second_start; 3],
        [second_end; 3],
        controls,
        cx,
    )
}

/// Admit orientation throughout one linear tetrahedron interpolation.
///
/// Endpoint sign is supplied by the immutable reference owner (+1 or -1).
/// This only excludes element inversion/collapse, not global mesh intersection
/// or constitutive-domain violations between native integration endpoints.
pub fn linear_tetrahedron_orientation(
    start: [[f64; 3]; 4],
    end: [[f64; 3]; 4],
    orientation: f64,
    controls: LinearSweepControls,
    cx: &Cx<'_>,
) -> Result<LinearSweepReport, TriangleQueryError> {
    super::poll(cx)?;
    validate(controls)?;
    if ![1.0, -1.0].contains(&orientation)
        || start
            .into_iter()
            .flatten()
            .chain(end.into_iter().flatten())
            .any(|x| !x.is_finite())
    {
        return Err(TriangleQueryError::Invalid(
            "finite tetrahedra and reference orientation +/-1 required",
        ));
    }
    let mut pending = vec![(0.0_f64, 1.0_f64)];
    let mut examined = 0;
    let mut minimum = f64::INFINITY;
    while let Some((lo, hi)) = pending.pop() {
        super::poll(cx)?;
        if examined == controls.maximum_windows {
            return Err(TriangleQueryError::WorkBudgetExhausted {
                maximum_visits: controls.maximum_windows,
            });
        }
        examined += 1;
        let t = Interval::new(lo, hi);
        // Differences first retain translation correlation. A relative edge is
        // affine; its hull is bounded by the two endpoint evaluations.
        let edges: [[Interval; 3]; 3] = core::array::from_fn(|i| {
            core::array::from_fn(|axis| {
                let a = Interval::point(start[i + 1][axis]) - Interval::point(start[0][axis]);
                let b = Interval::point(end[i + 1][axis]) - Interval::point(end[0][axis]);
                affine_hull(a, b, t.lo(), t.hi())
            })
        });
        let determinant = (edges[0][0] * (edges[1][1] * edges[2][2] - edges[1][2] * edges[2][1])
            - edges[0][1] * (edges[1][0] * edges[2][2] - edges[1][2] * edges[2][0])
            + edges[0][2] * (edges[1][0] * edges[2][1] - edges[1][1] * edges[2][0]))
            * Interval::point(orientation);
        if determinant.lo().is_finite() && determinant.lo() > 0.0 {
            minimum = minimum.min(determinant.lo());
            continue;
        }
        let mid = f64::midpoint(lo, hi);
        if hi - lo <= controls.fraction_tolerance || !(mid > lo && mid < hi) {
            return Ok(LinearSweepReport {
                verdict: LinearSweepVerdict::PossibleContact {
                    fraction_start: lo,
                    fraction_end: hi,
                },
                examined_windows: examined,
            });
        }
        pending.push((mid, hi));
        pending.push((lo, mid));
    }
    // Here the lower bound is signed six-volume, in coordinate units cubed.
    Ok(LinearSweepReport {
        verdict: LinearSweepVerdict::Clear {
            minimum_clearance: minimum,
        },
        examined_windows: examined,
    })
}

fn validate(controls: LinearSweepControls) -> Result<(), TriangleQueryError> {
    if controls.maximum_windows == 0
        || controls.maximum_windows > 1 << 20
        || !controls.fraction_tolerance.is_finite()
        || controls.fraction_tolerance <= 0.0
        || controls.fraction_tolerance >= 1.0
    {
        return Err(TriangleQueryError::Invalid(
            "bounded positive sweep work and fractional tolerance in (0,1) required",
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn plane_clearance(
    c0: [f64; 3],
    c1: [f64; 3],
    radius: f64,
    a: [[f64; 3]; 3],
    b: [[f64; 3]; 3],
    direction: [f64; 3],
    lo: f64,
    hi: f64,
) -> Option<f64> {
    let scale = direction.into_iter().map(f64::abs).fold(0.0, f64::max);
    if !scale.is_finite() || scale == 0.0 {
        return None;
    }
    let n = direction.map(|x| x / scale);
    let length = n
        .into_iter()
        .fold(Interval::point(0.0), |sum, x| {
            sum + Interval::point(x) * Interval::point(x)
        })
        .sqrt();
    if !length.hi().is_finite() || length.lo() <= 0.0 {
        return None;
    }
    let support = Interval::point(radius) * length;
    let mut lower = f64::INFINITY;
    for i in 0..3 {
        let project = |c: [f64; 3], vertex: [f64; 3]| {
            (0..3).fold(Interval::point(0.0), |sum, axis| {
                sum + (Interval::point(c[axis]) - Interval::point(vertex[axis]))
                    * Interval::point(n[axis])
            })
        };
        let projected = affine_hull(project(c0, a[i]), project(c1, b[i]), lo, hi);
        lower = lower.min(((projected - support) / length).lo());
    }
    (lower.is_finite() && lower >= 0.0).then_some(lower)
}

fn affine_hull(a: Interval, b: Interval, lo: f64, hi: f64) -> Interval {
    let at = |t: f64| {
        if t == 0.0 {
            a
        } else if t == 1.0 {
            b
        } else {
            (Interval::point(1.0) - Interval::point(t)) * a + Interval::point(t) * b
        }
    };
    let x = at(lo);
    let y = at(hi);
    Interval::new(x.lo().min(y.lo()), x.hi().max(y.hi()))
}

fn exact_axis_separation(
    c0: [f64; 3],
    c1: [f64; 3],
    radii: [f64; 2],
    a: [[f64; 3]; 3],
    b: [[f64; 3]; 3],
) -> bool {
    // The cold linear expansion tail never multiplies and cannot underflow.
    // This finite range preflight also keeps every intermediate sum finite.
    if radii.into_iter().any(|radius| radius > f64::MAX / 16.0)
        || c0
            .into_iter()
            .chain(c1)
            .chain(a.into_iter().flatten())
            .chain(b.into_iter().flatten())
            .any(|x| x.abs() > f64::MAX / 16.0)
    {
        return false;
    }
    for axis in 0..3 {
        for sign in [1.0, -1.0] {
            let separated = [(c0, a), (c1, b)].into_iter().all(|(c, triangle)| {
                triangle.into_iter().all(|vertex| {
                    let mut difference = diff_expansion(c[axis], vertex[axis]);
                    if sign < 0.0 {
                        for x in &mut difference {
                            *x = -*x;
                        }
                    }
                    for radius in radii {
                        if radius != 0.0 {
                            difference = fast_expansion_sum_zeroelim(&difference, &[-radius]);
                        }
                    }
                    expansion_sign(&difference) >= 0
                })
            });
            if separated {
                return true;
            }
        }
    }
    false
}

fn lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    core::array::from_fn(|i| (1.0 - t) * a[i] + t * b[i])
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    fn run(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
        let gate = CancelGate::new();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            f(
                &gate,
                &Cx::new(
                    &gate,
                    arena,
                    StreamKey {
                        seed: 7,
                        kernel_id: 1,
                        tile: 0,
                        iteration: 0,
                    },
                    Budget::INFINITE,
                    ExecMode::Deterministic,
                ),
            );
        });
    }
    const TRI: [[f64; 3]; 3] = [[-4.0, -4.0, 0.0], [4.0, -4.0, 0.0], [0.0, 4.0, 0.0]];
    #[test]
    fn affine_nearness_encloses_a_material_point_on_a_deforming_triangle() {
        run(|_, cx| {
            let end = [[-8.0, -4.0, 1.0], [8.0, -4.0, 1.0], [0.0, 12.0, 1.0]];
            let weights = [0.25, 0.25, 0.5];
            let gap = 1.0 / 1_048_576.0;
            // Fixed barycentric P moves from (0,0,0) to (0,4,1), while
            // every intermediate relative centre vector is (0,0,.5+gap).
            let upper = linear_sphere_triangle_nearness(
                [0.0, 0.0, 0.5 + gap],
                [0.0, 4.0, 1.5 + gap],
                0.5,
                TRI,
                end,
                weights,
                weights,
                cx,
            )
            .unwrap();
            assert!(upper >= gap && upper < gap + 1e-12);
            let receding = linear_sphere_triangle_nearness(
                [0.0, 0.0, 0.5 + gap],
                [0.0, 4.0, 2.5 + gap],
                0.5,
                TRI,
                end,
                weights,
                weights,
                cx,
            )
            .unwrap();
            assert!(receding >= 1.0 + gap);
            // A source history may not count that interval as sustained
            // contact simply because its initial point was near the surface.
            assert!(receding > 2.0 * gap);
            // Moving endpoint witnesses follow one metre of flat sliding.
            // Holding the start point fixed would give a false O(d^2/R) gap;
            // all three actual Bernstein controls have the same normal offset.
            let sliding = linear_sphere_triangle_nearness(
                [0.0, 0.0, 0.5 + gap],
                [1.0, 0.0, 0.5 + gap],
                0.5,
                TRI,
                TRI,
                weights,
                [0.125, 0.375, 0.5],
                cx,
            )
            .unwrap();
            assert!(sliding >= gap && sliding < gap + 1e-12);
        });
    }

    #[test]
    fn nearness_uses_original_material_coordinates_and_is_not_nonpenetration() {
        run(|_, cx| {
            let large = 4_503_599_627_370_496.0; // 2^52: half an ULP is lost.
            let triangle = [
                [large, 0.0, 0.0],
                [large + 1.0, 0.0, 0.0],
                [large, 1.0, 0.0],
            ];
            // The exact convex material point is (large+.5,0,0). Rounding
            // it first to a world-space float would erase its .5 separation.
            let upper = linear_sphere_triangle_nearness(
                [large, 0.0, 1.0],
                [large, 0.0, 1.0],
                1.0,
                triangle,
                triangle,
                [0.5, 0.5, 0.0],
                [0.5, 0.5, 0.0],
                cx,
            )
            .unwrap();
            assert!(upper >= 1.25_f64.sqrt() - 1.0 && upper < 0.119);
            let near = linear_sphere_triangle_nearness(
                [0.0, 0.0, 0.5],
                [0.0, 0.0, -0.5],
                0.5,
                TRI,
                TRI,
                [0.25, 0.25, 0.5],
                [0.25, 0.25, 0.5],
                cx,
            )
            .unwrap();
            assert!(near < 1e-12);
            assert!(matches!(
                linear_sphere_triangle(
                    [0.0, 0.0, 0.5],
                    [0.0, 0.0, -0.5],
                    0.5,
                    TRI,
                    TRI,
                    Default::default(),
                    cx,
                )
                .unwrap()
                .verdict,
                LinearSweepVerdict::PossibleContact { .. }
            ));
        });
    }

    #[test]
    fn changing_weights_need_the_quadratic_middle_control() {
        run(|_, cx| {
            let start = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
            let end = [[0.0; 3], [3.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
            // P_x(t)=(1-t)(1+2t), C_x(t)=1-t. Thus the relative x
            // control is -1 at its quadratic middle, despite zero endpoints.
            let upper = linear_sphere_triangle_nearness(
                [1.0, 0.0, 0.5],
                [0.0, 0.0, 0.5],
                0.5,
                start,
                end,
                [0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
                cx,
            )
            .unwrap();
            assert!(upper >= 1.25_f64.sqrt() - 0.5 && upper < 0.619);
        });
    }

    #[test]
    fn affine_pair_nearness_keeps_original_radii_and_bounded_range_refusals() {
        run(|gate, cx| {
            let gap = 1.0 / 1_048_576.0;
            let upper = linear_spheres_nearness(
                [0.0; 3],
                [1.0, 2.0, 3.0],
                0.125,
                [0.375 + gap, 0.0, 0.0],
                [1.375 + gap, 2.0, 3.0],
                0.25,
                cx,
            )
            .unwrap();
            assert!(upper >= gap && upper < gap + 1e-12);
            let large = linear_spheres_nearness(
                [1e200, 0.0, 0.0],
                [1e200, 0.0, 0.0],
                0.25e200,
                [0.0; 3],
                [0.0; 3],
                0.25e200,
                cx,
            )
            .unwrap();
            assert!(large >= 0.5e200 && large < 0.50000000000001e200);
            assert!(matches!(
                linear_spheres_nearness(
                    [f64::MAX, 0.0, 0.0],
                    [f64::MAX, 0.0, 0.0],
                    0.5,
                    [-f64::MAX, 0.0, 0.0],
                    [-f64::MAX, 0.0, 0.0],
                    0.5,
                    cx,
                ),
                Err(TriangleQueryError::NumericalRange)
            ));
            assert!(matches!(
                linear_sphere_triangle_nearness(
                    [0.0, 0.0, 1.0],
                    [0.0, 0.0, 1.0],
                    0.5,
                    TRI,
                    TRI,
                    [0.5; 3],
                    [0.5; 3],
                    cx,
                ),
                Err(TriangleQueryError::Invalid(_))
            ));
            gate.request();
            assert!(matches!(
                linear_spheres_nearness(
                    [0.0; 3],
                    [0.0; 3],
                    0.5,
                    [1.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    0.5,
                    cx,
                ),
                Err(TriangleQueryError::Cancelled)
            ));
        });
    }

    #[test]
    fn clear_endpoints_cannot_hide_a_pass_through() {
        run(|_, cx| {
            let report = linear_sphere_triangle(
                [0.0, 0.0, 2.0],
                [0.0, 0.0, -2.0],
                0.5,
                TRI,
                TRI,
                Default::default(),
                cx,
            )
            .unwrap();
            let LinearSweepVerdict::PossibleContact {
                fraction_start,
                fraction_end,
            } = report.verdict
            else {
                panic!("missed crossing")
            };
            // Exact first touch: 2 - 4t = 1/2, hence t = 3/8.
            assert!(fraction_start <= 0.375 && fraction_end >= 0.375 - 2e-8);
        });
    }
    #[test]
    fn exact_touch_receding_and_moving_surface_share_one_path_proof() {
        run(|_, cx| {
            let end = TRI.map(|p| [p[0] + 3.0, p[1], p[2] + 1.0]);
            let report = linear_sphere_triangle(
                [0.0, 0.0, 0.5],
                [3.0, 0.0, 2.0],
                0.5,
                TRI,
                end,
                Default::default(),
                cx,
            )
            .unwrap();
            assert!(matches!(report.verdict, LinearSweepVerdict::Clear { .. }));
            // One ULP of real penetration is never reclassified as touch.
            let report = linear_sphere_triangle(
                [0.0, 0.0, 0.5_f64.next_down()],
                [3.0, 0.0, 2.0],
                0.5,
                TRI,
                end,
                Default::default(),
                cx,
            )
            .unwrap();
            assert!(matches!(
                report.verdict,
                LinearSweepVerdict::PossibleContact {
                    fraction_start: 0.0,
                    ..
                }
            ));
        });
    }
    #[test]
    fn endpoint_positive_tetrahedra_can_collapse_twice_between_endpoints() {
        run(|_, cx| {
            let a = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
            let b = [
                [0.0; 3],
                [-2.0, 0.0, 0.0],
                [0.0, -0.5, 0.0],
                [0.0, 0.0, 1.0],
            ];
            // det = (1-3t)(1-1.5t): endpoints positive, interior inverted.
            let report = linear_tetrahedron_orientation(a, b, 1.0, Default::default(), cx).unwrap();
            assert!(
                matches!(report.verdict, LinearSweepVerdict::PossibleContact { fraction_start, .. } if fraction_start <= 1.0/3.0)
            );
            let translated = a.map(|p| [p[0] + 10.0, p[1] - 7.0, p[2] + 3.0]);
            assert!(matches!(
                linear_tetrahedron_orientation(a, translated, 1.0, Default::default(), cx)
                    .unwrap()
                    .verdict,
                LinearSweepVerdict::Clear { .. }
            ));
        });
    }
    #[test]
    fn sphere_pairs_need_a_whole_path_proof_including_grazing() {
        run(|_, cx| {
            let crossing = linear_spheres(
                [-1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                0.125,
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                0.125,
                Default::default(),
                cx,
            )
            .unwrap();
            assert!(
                matches!(crossing.verdict, LinearSweepVerdict::PossibleContact {
                fraction_start, .. } if fraction_start <= 0.4375)
            );
            let grazing = linear_spheres(
                [-1.0, 0.25, 0.0],
                [1.0, 0.25, 0.0],
                0.125,
                [0.0; 3],
                [0.0; 3],
                0.125,
                Default::default(),
                cx,
            )
            .unwrap();
            // Exact constant y separation proves nonpenetration even though
            // the trajectories touch at their interior closest approach.
            assert_eq!(
                grazing.verdict,
                LinearSweepVerdict::Clear {
                    minimum_clearance: 0.0
                }
            );
            let clear = linear_spheres(
                [-1.0, 0.5, 0.0],
                [1.0, 0.5, 0.0],
                0.125,
                [0.0; 3],
                [0.0; 3],
                0.125,
                Default::default(),
                cx,
            )
            .unwrap();
            assert!(matches!(clear.verdict, LinearSweepVerdict::Clear { .. }));
        });
    }

    #[test]
    fn sphere_pair_exact_touch_and_recession_keep_original_radii() {
        run(|_, cx| {
            for second_end in [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]] {
                let touching = linear_spheres(
                    [0.0; 3],
                    [0.0; 3],
                    0.5,
                    [1.0, 0.0, 0.0],
                    second_end,
                    0.5,
                    LinearSweepControls {
                        maximum_windows: 1,
                        ..Default::default()
                    },
                    cx,
                )
                .unwrap();
                assert_eq!(touching.examined_windows, 1);
                assert_eq!(
                    touching.verdict,
                    LinearSweepVerdict::Clear {
                        minimum_clearance: 0.0
                    }
                );
            }
            let penetrating = linear_spheres(
                [0.0; 3],
                [0.0; 3],
                0.5,
                [1.0_f64.next_down(), 0.0, 0.0],
                [2.0, 0.0, 0.0],
                0.5,
                Default::default(),
                cx,
            )
            .unwrap();
            assert!(matches!(
                penetrating.verdict,
                LinearSweepVerdict::PossibleContact {
                    fraction_start: 0.0,
                    ..
                }
            ));
        });
    }

    #[test]
    fn rounded_radius_sum_cannot_hide_real_pair_overlap() {
        run(|_, cx| {
            // As exact binary64 inputs, 0.1 + 0.2 is greater than 0.3.
            // Compare both original radii through expansions, never an epsilon
            // or an assumed exactness of the rounded radius sum.
            for (separation, clear) in [(0.3, false), (0.1 + 0.2, true)] {
                let report = linear_spheres(
                    [0.0; 3],
                    [0.0; 3],
                    0.1,
                    [separation, 0.0, 0.0],
                    [separation, 0.0, 0.0],
                    0.2,
                    Default::default(),
                    cx,
                )
                .unwrap();
                assert_eq!(
                    matches!(report.verdict, LinearSweepVerdict::Clear { .. }),
                    clear,
                    "separation={separation:.17}: {:?}",
                    report.verdict,
                );
            }
        });
    }
    #[test]
    fn exhausted_or_cancelled_queries_publish_no_clear_path() {
        run(|gate, cx| {
            assert!(matches!(
                linear_sphere_triangle(
                    [0.0, 0.0, 2.0],
                    [0.0, 0.0, -2.0],
                    0.5,
                    TRI,
                    TRI,
                    LinearSweepControls {
                        maximum_windows: 1,
                        ..Default::default()
                    },
                    cx
                ),
                Err(TriangleQueryError::WorkBudgetExhausted { .. })
            ));
            gate.request();
            assert!(matches!(
                linear_sphere_triangle(
                    [0.0, 0.0, 2.0],
                    [0.0, 0.0, 3.0],
                    0.5,
                    TRI,
                    TRI,
                    Default::default(),
                    cx
                ),
                Err(TriangleQueryError::Cancelled)
            ));
        });
    }
}
