use super::*;
use super::super::tests::{close, limits, model};
use fs_math::det;

#[test]
fn receiver_pressure_and_all_gradient_slopes_match_complete_independent_solves() {
    let surface = SpherePanels::icosphere(0.2, 0).unwrap();
    let n = surface.areas().len();
    let velocity: Vec<_> = surface.centroids().iter()
        .map(|x| C64::new(0.4 + x[2], 0.1 + 0.2 * x[0])).collect();
    let slope: Vec<_> = surface.centroids().iter()
        .map(|x| C64::new(-0.0001, 0.0002 * (1.0 + x[1]))).collect();
    let points = [[0.8, 0.1, 0.3], [-0.4, 0.7, 0.2]];
    let omega = 700.0;
    let h = 0.3;
    for relaxing in [false, true] {
        let material = model(relaxing);
        let response = material.solve_batch_with_frequency_derivative(
            &surface, omega, &[&velocity], Some(&[&slope]),
            Formulation::BurtonMiller, &limits(n), || false,
        ).unwrap().remove(0);
        let actual = response.observe(&surface, &points, 0.1, || false).unwrap();
        let original = response.response.observe(&surface, &points, 0.1, || false).unwrap();
        let samples: Vec<_> = [-2.0, -1.0, 1.0, 2.0].into_iter().map(|offset| {
            let drive: Vec<_> = velocity.iter().zip(&slope)
                .map(|(&v, &s)| v + s.scale(offset * h)).collect();
            material.solve_batch(
                &surface, omega + offset * h, &[&drive],
                Formulation::BurtonMiller, &limits(n), || false,
            ).unwrap().remove(0).observe(&surface, &points, 0.1, || false).unwrap()
        }).collect();
        let difference = |values: [C64; 4]| {
            (values[0] - values[1].scale(8.0) + values[2].scale(8.0) - values[3])
                .scale(1.0 / (12.0 * h))
        };
        for i in 0..points.len() {
            close(actual[i].field.pressure_pa, original[i].pressure_pa, 2e-13);
            let expected = difference(std::array::from_fn(|j| samples[j][i].pressure_pa));
            close(actual[i].pressure_derivative_pa_per_rad_s, expected, 3e-7);
            close(actual[i].pressure_log_frequency_derivative_s().unwrap(),
                expected / original[i].pressure_pa, 3e-7);
            assert_eq!(actual[i].pressure_group_delay_s().unwrap(),
                actual[i].pressure_log_frequency_derivative_s().unwrap().im);
            for axis in 0..3 {
                close(actual[i].field.gradient_pa_per_m[axis], original[i].gradient_pa_per_m[axis], 2e-13);
                let expected = difference(std::array::from_fn(|j| samples[j][i].gradient_pa_per_m[axis]));
                close(actual[i].gradient_derivative_pa_per_m_per_rad_s[axis], expected, 3e-7);
            }
        }
    }
}

#[test]
fn sphere_receiver_slope_and_group_delay_converge_to_radial_wave() {
    let material = model(true);
    let (omega, radius, distance) = (600.0, 0.2, 0.8);
    let wave = material.source().relaxing_frequency_response(omega).unwrap();
    let (k, dk) = (wave.wavenumber_per_m, wave.wavenumber_derivative_s_m);
    let rho = material.source().gas().density;
    let denominator = C64::new(-1.0 / radius - k.im, k.re);
    let q = C64::new(0.0, rho * omega);
    let z = q / denominator;
    let dz = C64::new(0.0, rho) / denominator
        - q * C64::new(-dk.im, dk.re) / (denominator * denominator);
    let phase = k.scale(distance - radius);
    let propagation = C64::new(det::cos(phase.re), det::sin(phase.re))
        .scale(det::exp(-phase.im) * radius / distance);
    let expected_p = z * propagation;
    let expected_dp = (dz + z * C64::new(-dk.im, dk.re).scale(distance - radius)) * propagation;
    let expected_gradient_slope = expected_dp * C64::new(-1.0 / distance - k.im, k.re)
        + expected_p * C64::new(-dk.im, dk.re);
    let expected_delay = (expected_dp / expected_p).im;
    let mut errors = Vec::new();
    for level in 0..=2 {
        let surface = SpherePanels::icosphere(radius, level).unwrap();
        let n = surface.areas().len();
        let velocity = vec![C64::ONE; n];
        let response = material.solve_batch_with_frequency_derivative(
            &surface, omega, &[&velocity], None,
            Formulation::BurtonMiller, &limits(n), || false,
        ).unwrap().remove(0);
        let observed = response.observe(&surface, &[[0.0, 0.0, distance]], 0.1, || false).unwrap()[0];
        errors.push((observed.pressure_derivative_pa_per_rad_s - expected_dp).abs() / expected_dp.abs());
        if level == 2 {
            close(observed.gradient_derivative_pa_per_m_per_rad_s[2], expected_gradient_slope, 0.12);
            assert!((observed.pressure_group_delay_s().unwrap() - expected_delay).abs()
                < 0.12 * expected_delay.abs());
        }
    }
    assert!(errors[1] < errors[0] && errors[2] < errors[1], "{errors:?}");
    assert!(errors[2] < 0.10, "{errors:?}");
}

#[test]
fn receiver_refuses_cross_wiring_invalid_tangents_undefined_phase_and_late_cancellation() {
    let material = model(true);
    let surface = SpherePanels::icosphere(0.2, 0).unwrap();
    let n = surface.areas().len();
    let velocity = vec![C64::ONE; n];
    let response = material.solve_batch_with_frequency_derivative(
        &surface, 700.0, &[&velocity], None,
        Formulation::BurtonMiller, &limits(n), || false,
    ).unwrap().remove(0);
    let points = [[0.8, 0.1, 0.3]];
    let different = SpherePanels::icosphere(0.21, 0).unwrap();
    assert!(response.observe(&different, &points, 0.1, || false).is_err());
    assert!(response.observe(&surface, &[[0.0; 3]], 0.1, || false).is_err());
    let mut incomplete = response.clone();
    incomplete.pressure_derivative_pa_per_rad_s.truncate(n - 1);
    assert!(incomplete.observe(&surface, &points, 0.1, || false).is_err());
    let mut invalid = response.clone();
    invalid.wavenumber_derivative_s_m.im = f64::NAN;
    assert!(invalid.observe(&surface, &points, 0.1, || false).is_err());
    let mut calls = 0;
    let mut sample = response.observe(&surface, &points, 0.1, || { calls += 1; false }).unwrap()[0];
    let mut seen = 0;
    assert!(matches!(response.observe(&surface, &points, 0.1, || { seen += 1; seen == calls }),
        Err(MolecularError::Boundary(HelmholtzError::Cancelled))));
    sample.field.pressure_pa = C64::ZERO;
    assert!(sample.pressure_group_delay_s().is_err());
    sample.field.pressure_pa = C64::from_re(f64::INFINITY);
    assert!(sample.pressure_log_frequency_derivative_s().is_err());
}
