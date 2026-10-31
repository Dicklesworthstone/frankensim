//! Source-distance scales and scattering angular-cell crossings of line segments.
//!
//! These floating-point roots guide quadrature; they are not certified
//! intersections or numerical error enclosures. The opposite cone nappe and
//! azimuth half-plane are removed by their signs, retaining ambiguous zeros.

/// Add a strictly interior normalized parameter without merging nearby roots.
pub(crate) fn insert_knot(knots: &mut Vec<f64>, value: f64, limit: usize) -> Result<(), KnotError> {
    if !value.is_finite() || !(0.0 < value && value < 1.0) {
        return Ok(());
    }
    let index = match knots.binary_search_by(|existing| existing.total_cmp(&value)) {
        Ok(_) => return Ok(()),
        Err(index) => index,
    };
    if knots.len() >= limit {
        return Err(KnotError::Budget);
    }
    if knots.len() == knots.capacity() {
        // Grow in small bounded blocks, never reserving the entire caller's
        // potentially very large panel budget for a single beam.
        knots
            .try_reserve_exact((limit - knots.len()).min(16))
            .map_err(|_| KnotError::Allocation)?;
    }
    knots.insert(index, value);
    Ok(())
}

#[derive(Debug, PartialEq, Eq)]
pub(crate) enum KnotError {
    Budget,
    Allocation,
}

/// Distance-based camera parameters: the closest segment point, followed by
/// offsets of one, two, four, ... closest-source distances along the ray.
/// The caller merges these with angular knots under one shared work budget.
#[derive(Clone, Copy, Debug)]
pub(crate) struct DistanceBreakpoints {
    center: f64,
    radius: f64,
    emit_center: bool,
    emit_left: bool,
}

/// Construct normalized parameters without allocating. The physical spacing
/// is d_min / |direction| in the ray's parameter, including a collinear source
/// beyond an endpoint. A source on the segment refuses; scene admission owns
/// the corresponding public singularity error before calling this helper.
pub(crate) fn distance_breakpoints(
    camera_to_source: [f64; 3],
    direction: [f64; 3],
    distance: f64,
) -> Option<DistanceBreakpoints> {
    if !distance.is_finite()
        || distance <= 0.0
        || camera_to_source
            .iter()
            .chain(&direction)
            .any(|value| !value.is_finite())
    {
        return None;
    }
    let scale = camera_to_source
        .iter()
        .fold(distance, |largest, value| largest.max(value.abs()));
    let q = camera_to_source.map(|value| value / scale);
    let distance = distance / scale;
    let norm = scaled_norm(direction);
    let extent = norm * distance;
    if !extent.is_finite() || extent <= 0.0 {
        return None;
    }
    let axis = direction.map(|value| value / norm);
    let projection = dot(q, axis).clamp(0.0, extent);
    let closest = std::array::from_fn(|i| axis[i].mul_add(-projection, q[i]));
    let radius = scaled_norm(closest) / extent;
    if radius.is_nan() || radius <= 0.0 {
        return None;
    }
    Some(DistanceBreakpoints {
        center: projection / extent,
        radius,
        // A segment shorter than one natural distance scale has no narrow
        // inverse-square feature to partition, including at its center.
        emit_center: radius < 1.0,
        emit_left: true,
    })
}

impl Iterator for DistanceBreakpoints {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.emit_center {
            self.emit_center = false;
            if let Some(center) = interior(self.center) {
                return Some(center);
            }
        }
        while self.radius < self.center.max(1.0 - self.center) {
            let root = if self.emit_left {
                self.center - self.radius
            } else {
                self.center + self.radius
            };
            if !self.emit_left {
                self.radius *= 2.0;
            }
            self.emit_left = !self.emit_left;
            if let Some(root) = interior(root) {
                return Some(root);
            }
        }
        None
    }
}

fn scaled_norm(value: [f64; 3]) -> f64 {
    let scale = value
        .iter()
        .fold(0.0_f64, |largest, component| largest.max(component.abs()));
    if scale == 0.0 {
        return 0.0;
    }
    let scaled = value.map(|component| component / scale);
    scale * dot(scaled, scaled).sqrt()
}

/// Crossings of `z/|r| = mu` for `r(u)=origin+u*delta`.
/// Coordinates have already been scaled to avoid squared-distance overflow.
pub(crate) fn cone_crossings(origin: [f64; 3], delta: [f64; 3], mu: f64) -> [Option<f64>; 2] {
    if mu == 0.0 {
        return [plane_crossing(origin, delta, [0.0, 0.0, 1.0]), None];
    }
    let magnitude = mu.abs();
    let sine_squared = (1.0 - magnitude) * (1.0 + magnitude);
    let cosine_squared = magnitude * magnitude;
    let radial_delta = delta[0].mul_add(delta[0], delta[1] * delta[1]);
    let a = difference_of_products(
        sine_squared,
        delta[2] * delta[2],
        cosine_squared,
        radial_delta,
    );
    let half_b = mixed_coefficient(origin, delta, sine_squared, cosine_squared);
    // Recompute the constant from the translated 3D point. Merely using a
    // stable quadratic formula on unshifted coefficients still loses narrow
    // beams when C subtracts two much larger nearly equal terms.
    let shift = if a == 0.0 {
        0.0
    } else {
        (-half_b / a).clamp(0.0, 1.0)
    };
    let shifted = std::array::from_fn(|i| delta[i].mul_add(shift, origin[i]));
    let half_b = mixed_coefficient(shifted, delta, sine_squared, cosine_squared);
    let radial_origin = shifted[0].mul_add(shifted[0], shifted[1] * shifted[1]);
    let c = difference_of_products(
        sine_squared,
        shifted[2] * shifted[2],
        cosine_squared,
        radial_origin,
    );
    quadratic_roots(a, half_b, c).map(|root| {
        root.and_then(|root| interior(root + shift)).filter(|u| {
            let z = delta[2].mul_add(*u, origin[2]);
            if mu < 0.0 { z <= 0.0 } else { z >= 0.0 }
        })
    })
}

/// Camera-segment crossing of a point source's propagation phase cosine.
/// `camera_to_source` and `distance` share one length scale; the caller scales
/// them before forming these products. The view direction is finite and near
/// unit length. Its actual norm is retained to match the incident-field dot
/// product, including the ray admission's permitted unit-vector roundoff.
pub(crate) fn phase_crossing(
    camera_to_source: [f64; 3],
    direction: [f64; 3],
    distance: f64,
    cosine: f64,
) -> Option<f64> {
    let norm = dot(direction, direction).sqrt();
    let mu = cosine / norm;
    if !(mu > -1.0 && mu < 1.0) {
        return None;
    }
    let axis = direction.map(|component| component / norm);
    let projection = dot(camera_to_source, axis);
    let perpendicular = [
        difference_of_products(camera_to_source[1], axis[2], camera_to_source[2], axis[1]),
        difference_of_products(camera_to_source[2], axis[0], camera_to_source[0], axis[2]),
        difference_of_products(camera_to_source[0], axis[1], camera_to_source[1], axis[0]),
    ];
    let impact = dot(perpendicular, perpendicular).sqrt();
    if impact == 0.0 {
        // A collinear source outside the admitted segment has constant phase.
        // Sources on the segment are rejected by scene singularity admission.
        return None;
    }
    // In the camera-axis frame mu=z/sqrt(b²+z²), z=c-|direction|*t.
    // Its inverse is monotone, with one root and no opposite cone nappe.
    let sine = ((1.0 - mu.abs()) * (1.0 + mu.abs())).sqrt();
    interior((-impact).mul_add(mu / sine, projection) / (norm * distance))
}

/// Crossings of one azimuth ray, excluding the opposite half of its plane.
pub(crate) fn azimuth_crossing(
    origin: [f64; 3],
    delta: [f64; 3],
    sine: f64,
    cosine: f64,
) -> Option<f64> {
    plane_crossing(origin, delta, [-sine, cosine, 0.0]).filter(|u| {
        let x = delta[0].mul_add(*u, origin[0]);
        let y = delta[1].mul_add(*u, origin[1]);
        cosine.mul_add(x, sine * y) >= 0.0
    })
}

/// Crossing of a plane through the fixture origin, including its opposite ray.
pub(crate) fn plane_crossing(origin: [f64; 3], delta: [f64; 3], normal: [f64; 3]) -> Option<f64> {
    let offset = dot(origin, normal);
    let slope = dot(delta, normal);
    if slope == 0.0 {
        None
    } else {
        interior(-offset / slope)
    }
}

fn mixed_coefficient(origin: [f64; 3], delta: [f64; 3], sine: f64, cosine: f64) -> f64 {
    difference_of_products(
        sine,
        origin[2] * delta[2],
        cosine,
        origin[0].mul_add(delta[0], origin[1] * delta[1]),
    )
}

fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left[0].mul_add(right[0], left[1].mul_add(right[1], left[2] * right[2]))
}

fn difference_of_products(a: f64, b: f64, c: f64, d: f64) -> f64 {
    let product = c * d;
    a.mul_add(b, -product) + (-c).mul_add(d, product)
}

fn interior(value: f64) -> Option<f64> {
    (value.is_finite() && 0.0 < value && value < 1.0).then_some(value)
}

// Solve A*x² + 2*half_B*x + C = 0. The roots can lie outside the segment.
fn quadratic_roots(a: f64, half_b: f64, c: f64) -> [Option<f64>; 2] {
    let scale = a.abs().max(half_b.abs()).max(c.abs());
    if scale == 0.0 || !scale.is_finite() {
        return [None, None];
    }
    let (a, half_b, c) = (a / scale, half_b / scale, c / scale);
    if a == 0.0 {
        return if half_b == 0.0 {
            [None, None]
        } else {
            [Some(-c / (2.0 * half_b)), None]
        };
    }
    let discriminant = difference_of_products(half_b, half_b, a, c);
    if discriminant < 0.0 {
        return [None, None];
    }
    if discriminant == 0.0 {
        return [Some(-half_b / a), None];
    }
    let q = -half_b - discriminant.sqrt().copysign(half_b);
    [Some(q / a), Some(c / q)]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn distance_knots(source: [f64; 3], direction: [f64; 3], length: f64) -> Vec<f64> {
        let mut knots = Vec::new();
        for root in distance_breakpoints(source, direction, length).unwrap() {
            insert_knot(&mut knots, root, 256).unwrap();
        }
        knots
    }

    #[test]
    fn g1_distance_scales_resolve_an_inverse_square_peak_missed_by_refinement() {
        let (b, c, length, sigma) = (0.01_f64, 5.110918_f64, 10.0_f64, 0.01_f64);
        let knots = distance_knots([b, 0.0, c], [0.0, 0.0, 1.0], length);
        assert_eq!(knots.len(), 19);
        // Actual thin-medium midpoint/Beer mathematics, using the production
        // distance knots. No support-independent stopping heuristic is changed.
        let integrate = |knots: &[f64], count: u32| {
            let (mut start, mut sum) = (0.0, 0.0);
            for boundary in knots.iter().copied().chain([1.0]) {
                let end = length * boundary;
                let width = (end - start) / f64::from(count);
                let weight = -(-sigma * width).exp_m1();
                for index in 0..count {
                    let t = start + (f64::from(index) + 0.5) * width;
                    let r = b.hypot(t - c);
                    let incident =
                        100.0 * (-sigma * r).exp() / (4.0 * core::f64::consts::PI * r * r);
                    sum += incident * weight * (-sigma * (start + f64::from(index) * width)).exp();
                }
                start = end;
            }
            sum
        };
        // Independent substitution theta=atan((t-c)/b), dt/r²=dtheta/b.
        // Standard transcendental functions and Simpson weights supply a
        // reference without production source, phase, or panel evaluators.
        let (lo, hi) = ((-c / b).atan(), ((length - c) / b).atan());
        let count = 32_768;
        let step = (hi - lo) / f64::from(count);
        let mut reference = 0.0;
        for index in 0..=count {
            let theta = lo + f64::from(index) * step;
            let t = c + b * theta.tan();
            let r = b / theta.cos();
            let weight = if index == 0 || index == count {
                1.0
            } else if index % 2 == 0 {
                2.0
            } else {
                4.0
            };
            reference += weight * (-sigma * (t + r)).exp();
        }
        reference *= step / 3.0 * (100.0 * sigma / (4.0 * core::f64::consts::PI * b));
        let coarse = integrate(&[], 32);
        let refined = integrate(&[], 64);
        assert!((coarse - refined).abs() / refined < 1e-4);
        assert!((refined / reference - 1.0).abs() > 0.45);
        let resolved = integrate(&knots, 32);
        assert!(
            (resolved / reference - 1.0).abs() < 5e-5,
            "{resolved} versus {reference}"
        );
        assert!((integrate(&knots, 64) - reference).abs() < (resolved - reference).abs());
    }

    #[test]
    fn g3_distance_knots_preserve_si_scale_rotation_and_actual_ray_speed() {
        let q = [0.01, 0.0, 5.110918];
        let baseline = distance_knots(q, [0.0, 0.0, 1.0], 10.0);
        for scale in [1e-100, 0.125, 128.0, 1e100] {
            let rotated = [q[2] * scale, q[1] * scale, -q[0] * scale];
            let actual = distance_knots(rotated, [1.0, 0.0, 0.0], 10.0 * scale);
            assert_eq!(actual.len(), baseline.len());
            for (actual, expected) in actual.iter().zip(&baseline) {
                assert!((actual - expected).abs() < 3e-14);
            }
        }
        let faster = distance_knots(q, [0.0, 0.0, 2.0], 5.0);
        assert_eq!(faster.len(), baseline.len());
        for (actual, expected) in faster.iter().zip(&baseline) {
            assert!((actual - expected).abs() < 3e-14);
        }
    }

    #[test]
    fn g0_distance_knots_handle_endpoints_and_share_the_existing_budget() {
        for (q, expected) in [([0.0, 0.0, -0.5], 0.5), ([0.0, 0.0, 1.5], 0.5)] {
            let knots = distance_knots(q, [0.0, 0.0, 1.0], 1.0);
            assert_eq!(knots.len(), 1);
            assert!((knots[0] - expected).abs() < 3e-15);
        }
        assert!(distance_knots([2.0, 0.0, 0.5], [0.0, 0.0, 1.0], 1.0).is_empty());
        for q in [[0.0; 3], [0.0, 0.0, 0.5], [0.0, 0.0, 1.0]] {
            assert!(distance_breakpoints(q, [0.0, 0.0, 1.0], 1.0).is_none());
        }
        assert!(distance_breakpoints([1.0; 3], [0.0; 3], 1.0).is_none());
        assert!(distance_breakpoints([1.0; 3], [0.0, 0.0, 1.0], 0.0).is_none());
        assert!(distance_breakpoints([f64::NAN; 3], [0.0, 0.0, 1.0], 1.0).is_none());
        let q = [0.01, 0.0, 5.110918];
        let mut knots = distance_knots(q, [0.0, 0.0, 1.0], 10.0);
        let limit = knots.len();
        let snapshot = knots.clone();
        for repeated in distance_breakpoints(q, [0.0, 0.0, 1.0], 10.0).unwrap() {
            insert_knot(&mut knots, repeated, limit).unwrap();
        }
        assert_eq!(knots, snapshot);
        let other = distance_breakpoints([0.01, 0.0, 4.9], [0.0, 0.0, 1.0], 10.0)
            .unwrap()
            .next()
            .unwrap();
        assert_eq!(
            insert_knot(&mut knots, other, limit),
            Err(KnotError::Budget)
        );
        assert_eq!(knots, snapshot);
    }

    fn sorted(roots: [Option<f64>; 2]) -> Vec<f64> {
        let mut roots: Vec<_> = roots.into_iter().flatten().collect();
        roots.sort_by(f64::total_cmp);
        roots
    }

    #[test]
    fn g1_phase_crossings_match_forward_and_backward_right_triangles() {
        let q = [0.5, 0.25, 0.0];
        let direction = [1.0, 0.0, 0.0];
        let diagonal_cosine = core::f64::consts::FRAC_1_SQRT_2;
        for (cosine, expected) in [
            (diagonal_cosine, 0.25),
            (0.0, 0.5),
            (-diagonal_cosine, 0.75),
        ] {
            let actual = phase_crossing(q, direction, 1.0, cosine).unwrap();
            assert!((actual - expected).abs() < 2.0e-16);
        }
        let norm = 1.0 + 4.0 * f64::EPSILON;
        let actual = phase_crossing(q, [norm, 0.0, 0.0], 1.0, norm * diagonal_cosine).unwrap();
        assert!((actual - 0.25 / norm).abs() < 2.0e-16);
    }

    #[test]
    fn g1_narrow_phase_cells_keep_distinct_boundaries_between_quadrature_nodes() {
        let q = [0.5, 0.3, 0.0];
        for center in [0.2_f64, -0.2] {
            let mut knots = Vec::new();
            for cosine in [center - 1.0e-5, center, center + 1.0e-5] {
                let actual = phase_crossing(q, [1.0, 0.0, 0.0], 1.0, cosine).unwrap();
                let expected = 0.5 - 0.3 / cosine.acos().tan();
                assert!((actual - expected).abs() < 3.0e-16);
                insert_knot(&mut knots, actual, 3).unwrap();
            }
            assert_eq!(knots.len(), 3);
            assert!(knots[0] < knots[1] && knots[1] < knots[2]);
            // The unpartitioned two- and four-panel midpoint evaluations
            // never encounter either of these nonzero phase lobes.
            for node in [0.125, 0.25, 0.375, 0.625, 0.75, 0.875] {
                assert!(node < knots[0] || node > knots[2]);
            }
        }
    }

    #[test]
    fn g3_phase_crossings_preserve_rigid_rotation_and_common_length_scale() {
        let rotate = |v: [f64; 3]| {
            [
                (v[0] + 2.0 * v[1] - 2.0 * v[2]) / 3.0,
                (2.0 * v[0] + v[1] + 2.0 * v[2]) / 3.0,
                (2.0 * v[0] - 2.0 * v[1] - v[2]) / 3.0,
            ]
        };
        let q = [0.55, 0.25, -0.125];
        for cosine in [-0.6, 0.0, 0.6] {
            let baseline = phase_crossing(q, [1.0, 0.0, 0.0], 1.0, cosine).unwrap();
            for scale in [1.0e-100, 0.125, 128.0, 1.0e100] {
                let actual = phase_crossing(
                    rotate(q.map(|component| component * scale)),
                    rotate([1.0, 0.0, 0.0]),
                    scale,
                    cosine,
                )
                .unwrap();
                assert!((actual - baseline).abs() < 8.0e-16);
            }
        }
    }

    #[test]
    fn g0_constant_phase_and_endpoint_crossings_add_no_interior_partition() {
        for q in [[-1.0, 0.0, 0.0], [2.0, 0.0, 0.0]] {
            for cosine in [-0.95, 0.0, 0.95] {
                assert_eq!(phase_crossing(q, [1.0, 0.0, 0.0], 1.0, cosine), None);
            }
        }
        for cosine in [-1.0, 1.0, f64::NAN] {
            assert_eq!(
                phase_crossing([0.5, 0.3, 0.0], [1.0, 0.0, 0.0], 1.0, cosine),
                None
            );
        }
        for projection in [0.0, 1.0] {
            assert_eq!(
                phase_crossing([projection, 0.3, 0.0], [1.0, 0.0, 0.0], 1.0, 0.0),
                None
            );
        }
    }

    #[test]
    fn g1_narrow_beam_crossings_match_independent_tangent_geometry() {
        let angle = 0.01_f64;
        let center = 1.0 / 16.0;
        let roots = sorted(cone_crossings(
            [-center, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            angle.cos(),
        ));
        assert_eq!(roots.len(), 2);
        for (actual, expected) in roots
            .into_iter()
            .zip([center - angle.tan(), center + angle.tan()])
        {
            assert!((actual - expected).abs() < 2.0e-14);
        }
    }

    #[test]
    fn g1_nearest_representable_polar_cutoff_keeps_two_distinct_crossings() {
        let mu = f64::from_bits(1.0_f64.to_bits() - 1);
        let center = 0.9;
        let radius = ((1.0 - mu) * (1.0 + mu)).sqrt() / mu;
        let roots = sorted(cone_crossings([-center, 0.0, 1.0], [1.0, 0.0, 0.0], mu));
        assert_eq!(roots, [center - radius, center + radius]);
    }

    #[test]
    fn g1_azimuth_and_hemisphere_planes_match_linear_intersections() {
        let angle = 0.2_f64;
        let (sine, cosine) = angle.sin_cos();
        let root = azimuth_crossing([1.0, -0.25, 0.5], [0.0, 1.0, 0.0], sine, cosine).unwrap();
        assert!((root - (0.25 + angle.tan())).abs() < 2.0e-16);
        assert_eq!(
            azimuth_crossing([-1.0, -0.25, 0.5], [0.0, 1.0, 0.0], sine, cosine),
            None
        );
        assert_eq!(
            cone_crossings([1.0, 0.0, -0.25], [0.0, 0.0, 1.0], 0.0),
            [Some(0.25), None]
        );
        assert_eq!(
            plane_crossing([0.0; 3], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
            None
        );
    }

    #[test]
    fn g0_tangent_linear_degenerate_and_endpoint_cases_are_finite() {
        assert_eq!(
            cone_crossings([-0.5, 0.0, 1.0], [1.0, 0.0, 0.0], 1.0),
            [Some(0.5), None]
        );
        assert_eq!(quadratic_roots(0.0, 1.0, -0.5), [Some(0.25), None]);
        assert_eq!(quadratic_roots(0.0, 0.0, 0.0), [None, None]);
        assert_eq!(
            cone_crossings([0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 0.5),
            [None, None]
        );
        assert_eq!(
            plane_crossing([0.0; 3], [0.0, 0.0, 1.0], [0.0, 0.0, 1.0]),
            None
        );
        assert_eq!(
            plane_crossing([0.0, 0.0, -1.0], [0.0, 0.0, 1.0], [0.0, 0.0, 1.0]),
            None
        );
    }

    #[test]
    fn g3_reflecting_the_cone_axis_and_scaling_coordinates_preserve_boundaries() {
        let baseline = sorted(cone_crossings([-0.2, 0.0, 0.5], [1.0, 0.0, 0.0], 0.9));
        for scale in [1.0e-100, 0.125, 128.0, 1.0e100] {
            for mu in [-0.9, 0.9] {
                let roots = sorted(cone_crossings(
                    [-0.2 * scale, 0.0, (0.5 * scale).copysign(mu)],
                    [scale, 0.0, 0.0],
                    mu,
                ));
                assert_eq!(roots.len(), baseline.len());
                for (actual, expected) in roots.into_iter().zip(&baseline) {
                    assert!((actual - expected).abs() < 5.0e-16);
                }
            }
        }
    }

    #[test]
    fn g0_equatorial_table_boundary_does_not_keep_the_opposite_nappe() {
        let mu = core::f64::consts::FRAC_PI_2.cos();
        let roots = sorted(cone_crossings([1.0, 0.0, -0.5], [0.0, 0.0, 1.0], mu));
        assert_eq!(roots.len(), 1);
        assert!(roots[0] >= 0.5);
        let negative = sorted(cone_crossings([1.0, 0.0, -0.5], [0.0, 0.0, 1.0], -mu));
        assert_eq!(negative.len(), 1);
        assert!(negative[0] <= 0.5);
    }

    #[test]
    fn g0_exact_dedup_preserves_thin_intervals_and_global_budget() {
        let first = 0.5;
        let adjacent = f64::from_bits(0.5_f64.to_bits() + 1);
        let mut knots = Vec::new();
        insert_knot(&mut knots, adjacent, 2).unwrap();
        insert_knot(&mut knots, first, 2).unwrap();
        insert_knot(&mut knots, first, 2).unwrap();
        assert_eq!(knots, [first, adjacent]);
        assert_eq!(insert_knot(&mut knots, 0.75, 2), Err(KnotError::Budget));
        assert_eq!(knots, [first, adjacent]);
        for outside in [0.0, 1.0, f64::NAN, f64::INFINITY] {
            insert_knot(&mut knots, outside, 2).unwrap();
        }
        assert_eq!(knots, [first, adjacent]);
    }
}
