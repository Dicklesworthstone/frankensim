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
