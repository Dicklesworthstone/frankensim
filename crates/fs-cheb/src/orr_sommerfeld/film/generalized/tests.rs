use super::*;
use crate::orr_sommerfeld::film::FilmSpectrumBudget;

fn budget() -> FilmSpectrumBudget {
    FilmSpectrumBudget {
        coarse_order: 16,
        fine_order: 24,
        modes: 1,
        max_qr_sweeps: 256,
        max_svd_sweeps: 64,
        relative_residual_tolerance: 1.0e-8,
        absolute_resolution_tolerance: 1.0e-7,
        relative_resolution_tolerance: 1.0e-6,
    }
}
fn parameters() -> FilmParameters {
    FilmParameters {
        reynolds: 1.0,
        kh: 0.5,
        normal_gravity: 0.0,
        capillarity: 10.0,
    }
}
fn samples(
    solver: &PreparedFilmSpectrum,
    source: impl Fn(f64) -> FilmProfilePoint,
) -> (Vec<FilmProfilePoint>, Vec<FilmProfilePoint>) {
    let (coarse, fine) = solver.profile_coordinates();
    (
        coarse.iter().map(|&y| source(y)).collect(),
        fine.iter().map(|&y| source(y)).collect(),
    )
}

#[test]
fn original_constant_profile_recovers_the_native_nusselt_spectrum() {
    let solver = PreparedFilmSpectrum::new(budget(), || false).unwrap();
    let (coarse, fine) = samples(&solver, |y| FilmProfilePoint {
        velocity: 3.0 * y - 1.5 * y * y,
        velocity_gradient: 3.0 * (1.0 - y),
        velocity_curvature: -3.0,
        apparent_viscosity: 1.0,
        stress_tangent: 1.0,
    });
    let prepared = solver
        .prepare_profile(&coarse, &fine, 0.0, || false)
        .unwrap();
    let actual = prepared.solve(parameters(), || false).unwrap()[0].eigenvalue;
    let native = solver.solve(parameters(), || false).unwrap()[0].eigenvalue;
    assert!((actual - native).abs() < 2.0e-9, "{actual:?} != {native:?}");
    assert!(matches!(
        solver.prepare_profile(&coarse, &fine, 0.0, || true),
        Err(FilmSpectrumError::Cancelled)
    ));
    assert!(matches!(
        solver.prepare_profile(&coarse[1..], &fine, 0.0, || false),
        Err(FilmSpectrumError::Invalid(_))
    ));
}

#[test]
fn resolved_shear_thinning_matches_independent_qz_and_exact_flux_wave_speed() {
    let solver = PreparedFilmSpectrum::new(budget(), || false).unwrap();
    let beta: f64 = 0.5;
    let t = 2.0 * beta / (beta.sqrt().asin() / beta.sqrt() - (1.0 - beta).sqrt());
    // Exact monotone stress law mu=(1+lambda²*gamma²)^(-1/2), below
    // its finite stress ceiling. mu=sqrt(1-beta*z²), M=mu³, z=1-y.
    // T is chosen analytically so integral U dy=1. No fitted base is used.
    let (coarse, fine) = samples(&solver, |y| {
        let z = 1.0 - y;
        let root = (1.0 - beta * z * z).sqrt();
        FilmProfilePoint {
            velocity: t / beta * (root - (1.0 - beta).sqrt()),
            velocity_gradient: t * z / root,
            velocity_curvature: -t / (root * root * root),
            apparent_viscosity: root,
            stress_tangent: root * root * root,
        }
    });
    let prepared = solver
        .prepare_profile(&coarse, &fine, 0.0, || false)
        .unwrap();
    // Independent generalized-QZ solve retaining eta as an additional unknown;
    // its bulk viscosity derivatives were analytic, not nodally differentiated.
    let expected = C64::new(-0.024_616_882_32, -1.304_593_661_62);
    let actual = prepared.solve(parameters(), || false).unwrap()[0].eigenvalue;
    assert!(
        (actual - expected).abs() < 2.0e-6,
        "{actual:?} != {expected:?}"
    );
    let alpha = 0.0001;
    let long_wave = prepared
        .solve(
            FilmParameters {
                kh: alpha,
                capillarity: 0.0,
                ..parameters()
            },
            || false,
        )
        .unwrap()[0]
        .eigenvalue;
    // Original nonlinear flux obeys dq/dh=h*gamma_wall. This wave speed
    // is 3.503876..., not the Newtonian value 3 on mean-speed scales.
    let exact_speed = t / (1.0 - beta).sqrt();
    assert!((-long_wave.im / alpha - exact_speed).abs() < 2.0e-6);
}

#[test]
fn nonzero_surface_tangent_gradient_is_retained_in_normal_stress() {
    let solver = PreparedFilmSpectrum::new(budget(), || false).unwrap();
    let beta: f64 = 0.3;
    let t = 1.0 / (-1.0 / (2.0 * beta) - 1.0 / beta.powi(2) - (1.0 - beta).ln() / beta.powi(3));
    // Another explicit monotone law: mu=1/(1+a*gamma), on the finite
    // admitted range. Its base has mu=1-beta*z, M=mu² and M'_surface=2beta.
    // All one-sided coefficient derivatives are finite in this fixture.
    let (coarse, fine) = samples(&solver, |y| {
        let z = 1.0 - y;
        let mu = 1.0 - beta * z;
        FilmProfilePoint {
            velocity: t / beta.powi(2) * ((mu / (1.0 - beta)).ln() - beta * y),
            velocity_gradient: t * z / mu,
            velocity_curvature: -t / (mu * mu),
            apparent_viscosity: mu,
            stress_tangent: mu * mu,
        }
    });
    let prepared = solver
        .prepare_profile(&coarse, &fine, 2.0 * beta, || false)
        .unwrap();
    let actual = prepared.solve(parameters(), || false).unwrap()[0].eigenvalue;
    let expected = C64::new(-0.038_448_923_26, -1.280_104_790_63);
    assert!(
        (actual - expected).abs() < 2.0e-6,
        "{actual:?} != {expected:?}"
    );
    let omitted = solver
        .prepare_profile(&coarse, &fine, 0.0, || false)
        .unwrap()
        .solve(parameters(), || false)
        .unwrap()[0]
        .eigenvalue;
    assert!((actual - omitted).abs() > 0.01);
}
