use super::*;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

fn unit(v: [f64; 3]) -> [f64; 3] {
    let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    v.map(|x| x / norm)
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}

fn independent_sag(c: f64, k: f64, a: [f64; 7], x: f64, y: f64) -> f64 {
    let r2 = x * x + y * y;
    c * r2 / (1.0 + (1.0 - (1.0 + k) * c * c * r2).sqrt())
        + a.into_iter()
            .enumerate()
            .map(|(i, a)| a * r2.powi((i + 2) as i32))
            .sum::<f64>()
}

fn coefficients() -> [f64; 7] {
    [-2e5, 3e8, -5e11, 2e14, -5e16, 2e19, -3e21]
}

#[test]
fn analytic_conics_keep_sphere_parabola_and_plane_limits() {
    let c = 20.0;
    let r = 0.012;
    for k in [0.0, -1.0, -2.0, 0.4] {
        let surface = AsphericSurface::try_new(c, k, [0.0; 7], 0.02).unwrap();
        let (point, distance) = surface
            .intersect(0.08, [r, 0.0, 0.0], [0.0, 0.0, 1.0])
            .unwrap()
            .unwrap();
        let sag = independent_sag(c, k, [0.0; 7], r, 0.0);
        close(point[2], 0.08 + sag, 2e-16);
        close(distance, point[2], 0.0);
        let normal = surface.normal([r, 0.0]).unwrap();
        close(
            normal[0],
            c * r / (1.0 - (1.0 + k) * c * c * r * r).sqrt(),
            3e-16,
        );
        assert_eq!(normal[1..], [0.0, -1.0]);
    }
    let plane = AsphericSurface::try_new(0.0, -2.0, [0.0; 7], 0.02).unwrap();
    assert_eq!(
        plane
            .intersect(0.08, [0.01, 0.0, 0.0], [0.0, 0.0, 1.0])
            .unwrap(),
        Some(([0.01, 0.0, 0.08], 0.08))
    );
}

#[test]
fn analytic_conic_focus_removes_single_interface_spherical_aberration() {
    let eta = 1.0 / 1.5;
    let focus = 0.25;
    let c = 1.0 / (focus * (1.0 - eta));
    // ni*z + nt*sqrt(r²+(f-z)²)=nt*f gives k=-(ni/nt)² and
    // vertex radius f*(1-ni/nt), independently of the ray solver.
    let asphere = AsphericSurface::try_new(c, -eta * eta, [0.0; 7], 0.035).unwrap();
    let sphere = AsphericSurface::try_new(c, 0.0, [0.0; 7], 0.035).unwrap();
    let mut spherical_error = 0.0_f64;
    for radius in [0.005, 0.015, 0.03] {
        for (surface, corrected) in [(asphere, true), (sphere, false)] {
            let (point, _) = surface
                .intersect(0.08, [radius, 0.0, 0.0], [0.0, 0.0, 1.0])
                .unwrap()
                .unwrap();
            let normal = unit(surface.normal([radius, 0.0]).unwrap());
            let cosine = -normal[2];
            let transmitted_cosine = (1.0 - eta * eta * (1.0 - cosine * cosine)).sqrt();
            let ray = core::array::from_fn::<_, 3, _>(|j| {
                (if j == 2 { eta } else { 0.0 }) + (eta * cosine - transmitted_cosine) * normal[j]
            });
            let intercept = point[2] - radius * ray[2] / ray[0];
            if corrected {
                close(intercept, 0.08 + focus, 1e-14);
            } else {
                spherical_error = spherical_error.max((intercept - 0.08 - focus).abs());
            }
        }
    }
    assert!(spherical_error > 1e-3);
}

#[test]
fn polynomial_roots_and_normals_match_independent_known_geometry() {
    let (c, k, a, radius, vertex) = (20.0, -0.7, coefficients(), 0.015, 0.1);
    let surface = AsphericSurface::try_new(c, k, a, radius).unwrap();
    let direction = unit([0.02, 0.03, 1.0]);
    for (x, y) in [(0.0, 0.0), (0.007, -0.004), (-0.004, 0.009)] {
        let z = vertex + independent_sag(c, k, a, x, y);
        let expected = z / direction[2];
        let origin = [
            x - direction[0] * expected,
            y - direction[1] * expected,
            0.0,
        ];
        let (point, distance) = surface
            .intersect(vertex, origin, direction)
            .unwrap()
            .unwrap();
        close(distance, expected, 2e-15);
        for (actual, expected) in point.into_iter().zip([x, y, z]) {
            close(actual, expected, 2e-15);
        }
        let normal = surface.normal([x, y]).unwrap();
        let h = 1e-7;
        for j in 0..2 {
            let (dx, dy) = if j == 0 { (h, 0.0) } else { (0.0, h) };
            let slope = (independent_sag(c, k, a, x + dx, y + dy)
                - independent_sag(c, k, a, x - dx, y - dy))
                / (2.0 * h);
            close(normal[j], slope, 4e-10);
        }
    }
}

#[test]
fn analytic_bounds_include_nonmonotone_sag_and_keep_coefficient_units() {
    let (c, k, a, radius) = (20.0, -0.7, coefficients(), 0.015);
    let profile = AsphericSurface::try_new(c, k, a, radius).unwrap();
    let [low, high] = profile.sag_bounds();
    assert!(low < 0.0 && high > 0.0);
    // Samples are an independent control; production admission uses the
    // analytic full-disk term bounds and does not infer extrema from samples.
    for i in 0..=127 {
        let r = radius * i as f64 / 127.0;
        let expected = independent_sag(c, k, a, r, 0.0);
        assert!(expected >= low && expected <= high);
        close(profile.sag(r * r), expected, 8e-18);
    }
    for scale in [1e-3_f64, 1e3] {
        let scaled = AsphericSurface::try_new(
            c / scale,
            k,
            core::array::from_fn(|j| a[j] / scale.powi((2 * j + 3) as i32)),
            radius * scale,
        )
        .unwrap();
        let x = 0.007;
        close(
            scaled.sag(x * x * scale * scale) / scale,
            profile.sag(x * x),
            8e-18,
        );
        let normal = scaled.normal([x * scale, 0.0]).unwrap();
        close(normal[0], profile.normal([x, 0.0]).unwrap()[0], 2e-15);
        let original = profile
            .intersect(0.1, [x, 0.0, 0.0], [0.0, 0.0, 1.0])
            .unwrap()
            .unwrap();
        let changed = scaled
            .intersect(0.1 * scale, [x * scale, 0.0, 0.0], [0.0, 0.0, 1.0])
            .unwrap()
            .unwrap();
        close(changed.1 / scale, original.1, 2e-15);
    }
}

#[test]
fn rear_surface_solid_angle_jacobian_matches_direction_derivatives() {
    let c = 20.0;
    let k = -0.7;
    let a = coefficients();
    let profile = AsphericSurface::try_new(c, k, a, 0.015).unwrap();
    let sensor = [-0.01, 0.007, 0.0];
    let vertex = 0.08;
    let ray = |x, y| {
        unit([
            x - sensor[0],
            y - sensor[1],
            vertex + independent_sag(c, k, a, x, y),
        ])
    };
    for (x, y) in [(0.001, 0.002), (0.008, -0.004)] {
        let d = ray(x, y);
        let h = 1e-7;
        let dx =
            core::array::from_fn::<_, 3, _>(|j| (ray(x + h, y)[j] - ray(x - h, y)[j]) / (2.0 * h));
        let dy =
            core::array::from_fn::<_, 3, _>(|j| (ray(x, y + h)[j] - ray(x, y - h)[j]) / (2.0 * h));
        let cross = [
            dx[1] * dy[2] - dx[2] * dy[1],
            dx[2] * dy[0] - dx[0] * dy[2],
            dx[0] * dy[1] - dx[1] * dy[0],
        ];
        let numerical = dot(d, cross).abs();
        let normal = profile.normal([x, y]).unwrap();
        let v = [
            x - sensor[0],
            y - sensor[1],
            vertex + profile.sag(x * x + y * y),
        ];
        let analytic = dot(normal, d).abs() / normal[2].abs() / dot(v, v);
        close(numerical / analytic, 1.0, 3e-9);
        // Multiplying this dOmega/dA by the sensor cosine and vertex² gives
        // the unchanged native aperture-normalized geometry response.
        assert!(d[2] * analytic * vertex * vertex > 0.0);
    }
}

#[test]
fn aperture_misses_and_unqualified_crossings_remain_distinct() {
    let profile = AsphericSurface::try_new(20.0, -0.7, coefficients(), 0.015).unwrap();
    assert!(
        profile
            .intersect(0.1, [0.05, 0.0, 0.0], [0.0, 0.0, 1.0])
            .unwrap()
            .is_none()
    );
    assert!(
        profile
            .intersect(0.1, [0.0; 3], unit([1.0, 0.0, 0.1]))
            .is_err()
    );
    assert!(
        profile
            .intersect(0.1, [0.0; 3], [0.0, 0.0, -1.0])
            .unwrap()
            .is_none()
    );
    assert!(
        profile
            .intersect(0.1, [f64::NAN, 0.0, 0.0], [0.0, 0.0, 1.0])
            .is_err()
    );
}

#[test]
fn reverse_conic_and_polynomial_intersections_recover_known_geometry() {
    let vertex = 0.1;
    let direction = unit([0.02, -0.03, 1.0]);
    let reverse = direction.map(|component| -component);
    for (curvature, conic, coefficients) in [
        (20.0, -0.7, [0.0; 7]),
        (0.0, -1.0, [0.0; 7]),
        (-12.0, 0.3, coefficients()),
        (20.0, -0.7, coefficients()),
    ] {
        let profile = AsphericSurface::try_new(curvature, conic, coefficients, 0.015).unwrap();
        for (x, y) in [(0.0, 0.0), (0.007, -0.004), (-0.004, 0.009)] {
            let expected = [
                x,
                y,
                vertex + independent_sag(curvature, conic, coefficients, x, y),
            ];
            let from_sensor = core::array::from_fn(|j| expected[j] - 0.07 * direction[j]);
            let from_scene = core::array::from_fn(|j| expected[j] - 0.13 * reverse[j]);
            let (forward_hit, forward_distance) = profile
                .intersect(vertex, from_sensor, direction)
                .unwrap()
                .unwrap();
            let (backward_hit, backward_distance) = profile
                .intersect(vertex, from_scene, reverse)
                .unwrap()
                .unwrap();
            close(forward_distance, 0.07, 3e-15);
            close(backward_distance, 0.13, 3e-15);
            for coordinate in 0..3 {
                close(forward_hit[coordinate], expected[coordinate], 3e-15);
                close(backward_hit[coordinate], expected[coordinate], 3e-15);
            }
        }
    }
}

#[test]
fn reverse_aperture_miss_and_unproved_crossing_remain_distinct() {
    let profile = AsphericSurface::try_new(20.0, -0.7, coefficients(), 0.015).unwrap();
    assert!(
        profile
            .intersect(0.1, [0.05, 0.0, 0.2], [0.0, 0.0, -1.0])
            .unwrap()
            .is_none()
    );
    assert!(
        profile
            .intersect(0.1, [0.0, 0.0, 0.2], unit([1.0, 0.0, -0.1]))
            .is_err()
    );
}

#[test]
fn conic_domain_and_unrepresentable_prescriptions_refuse() {
    assert!(AsphericSurface::try_new(100.0, 0.0, [0.0; 7], 0.01).is_none());
    assert!(AsphericSurface::try_new(20.0, f64::NAN, [0.0; 7], 0.01).is_none());
    assert!(AsphericSurface::try_new(20.0, -1.0, [f64::INFINITY; 7], 0.01).is_none());
    assert!(AsphericSurface::try_new(20.0, -1.0, [0.0; 7], 1e-200).is_none());
    assert!(AsphericSurface::try_new(20.0, -1.0, [1e-300; 7], 1e-10).is_none());
    // A hyperboloid need not obey the spherical |c*R|<1 restriction.
    let hyperbola = AsphericSurface::try_new(100.0, -2.0, [0.0; 7], 0.02).unwrap();
    assert!(
        hyperbola
            .intersect(0.1, [0.015, 0.0, 0.0], [0.0, 0.0, 1.0])
            .unwrap()
            .is_some()
    );
}
