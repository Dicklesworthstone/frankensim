//! Angular-cell crossings of a scaled fixture-local line segment.
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

    fn sorted(roots: [Option<f64>; 2]) -> Vec<f64> {
        let mut roots: Vec<_> = roots.into_iter().flatten().collect();
        roots.sort_by(f64::total_cmp);
        roots
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
