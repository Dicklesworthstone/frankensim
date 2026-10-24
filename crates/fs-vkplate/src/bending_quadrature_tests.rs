//! Bending discretization checks independent of a recorded instrument.
use super::*;

fn uniform_plate(h: f64, nu: f64) -> Body {
    Body::plate_without_surface_loss(
        "bending test",
        &PlateSpec {
            radius: 0.16,
            hole: 0.0,
            grid: h,
            edge: Edge::Free,
            thickness: &|_| 0.001,
            rest_height: &|_| 0.0,
            // Hold D=1 N m; the constant-D bulk equation is independent of nu.
            young: 12.0 * (1.0 - nu * nu) / 0.001_f64.powi(3),
            poisson: nu,
            density: 8000.0,
            tension: 0.0,
            bending: true,
            nonlinear: Nonlinear::None,
            t60: [(100.0, 1.0), (1000.0, 1.0)],
            support: None,
            muffle: None,
        },
    )
    .unwrap()
}

#[test]
fn g0_uniform_bending_satisfies_discrete_integration_by_parts() {
    let h = 0.02_f64;
    for nu in [-0.2, 0.0, 0.33, 0.49] {
        let mut body = uniform_plate(h, nu);
        let centre = body.xy.iter().position(|p| *p == [0.0, 0.0]).unwrap();
        for (ax, ay) in [(0.8_f64, 0.0_f64), (0.6, 0.8), (1.4, 1.1)] {
            let w: Vec<_> = body
                .xy
                .iter()
                .map(|[x, y]| (ax * x / h + ay * y / h).cos())
                .collect();
            let actual = body.linear_stiffness_apply(&w).unwrap()[centre];
            // Independent Fourier symbol of the five-point Laplacian squared.
            // Centered mixed derivatives incorrectly suppress diagonal waves.
            let lap = 4.0 * ((ax / 2.0).sin().powi(2) + (ay / 2.0).sin().powi(2)) / h.powi(2);
            let expected = h.powi(2) * lap.powi(2);
            assert!(
                (actual / expected - 1.0).abs() < 2e-12,
                "nu={nu}, wave=({ax},{ay}), actual={actual}, expected={expected}"
            );
        }
    }
}

#[test]
fn g1_bending_wave_force_converges_to_kirchhoff_love_at_second_order() {
    let (kx, ky) = (31.0_f64, 23.0_f64);
    let continuum = (kx * kx + ky * ky).powi(2);
    let mut previous = None;
    for h in [0.02_f64, 0.01, 0.005] {
        let mut body = uniform_plate(h, 0.33);
        let centre = body.xy.iter().position(|p| *p == [0.0, 0.0]).unwrap();
        let w: Vec<_> = body
            .xy
            .iter()
            .map(|[x, y]| (kx * x + ky * y).cos())
            .collect();
        let force_density = body.linear_stiffness_apply(&w).unwrap()[centre] / h.powi(2);
        let error = (force_density / continuum - 1.0).abs();
        if let Some(old) = previous {
            assert!(
                error < 0.27 * old && error > 0.23 * old,
                "h={h}, relative error={error}, previous={old}"
            );
        }
        previous = Some(error);
    }
}

#[test]
fn g0_free_bending_keeps_rigid_motion_and_sparse_evaluation_order() {
    let body = uniform_plate(0.02, 0.33);
    for field in [
        body.xy.iter().map(|_| 1.0).collect::<Vec<_>>(),
        body.xy
            .iter()
            .map(|[x, y]| 0.3 + 2.0 * x - 0.7 * y)
            .collect(),
    ] {
        let mut force = vec![0.0; body.n];
        body.stiff.mul_add(&field, 1.0, &mut force);
        for (i, f) in force.iter().enumerate() {
            let rounding_scale: f64 = (body.stiff.ptr[i]..body.stiff.ptr[i + 1])
                .map(|q| (body.stiff.val[q] * field[body.stiff.col[q]]).abs())
                .sum();
            assert!(f.abs() <= 32.0 * f64::EPSILON * rounding_scale);
        }
    }
    let field: Vec<_> = body
        .xy
        .iter()
        .enumerate()
        .map(|(i, [x, y])| (17.0 * x + 31.0 * y).sin() + if i % 2 == 0 { 0.7 } else { -0.4 })
        .collect();
    let mut optimized = vec![0.0; body.n];
    body.stiff.mul_add(&field, 1.0, &mut optimized);
    for (i, actual) in optimized.iter().enumerate() {
        let mut expected = 0.0;
        for q in body.stiff.ptr[i]..body.stiff.ptr[i + 1] {
            expected += body.stiff.val[q] * field[body.stiff.col[q]];
        }
        assert_eq!(
            actual.to_bits(),
            expected.to_bits(),
            "CSR order changed at row {i}"
        );
    }
    assert!(dot(&field, &optimized) > 0.0);
}
