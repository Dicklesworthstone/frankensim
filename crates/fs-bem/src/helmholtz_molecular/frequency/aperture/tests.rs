use super::*;
use super::super::tests::{close, limits, model};

fn openings(surface: &SpherePanels) -> Vec<Aperture> {
    [-1.0, 1.0].into_iter().map(|side| Aperture {
        shape: surface.centroids().iter()
            .map(|p| if side * p[2] > 0.0 { 1.0 } else { 0.0 }).collect(),
        reference_point: [0.0, 0.0, side * 0.2],
        normal: [0.0, 0.0, side],
    }).collect()
}

#[test]
fn complete_aperture_load_and_flow_receiver_derivatives_match_original_frequency_solves() {
    let material = model(true);
    let surface = SpherePanels::icosphere(0.2, 0).unwrap();
    let n = surface.areas().len();
    let apertures = openings(&surface);
    let omega = 700.0;
    let h = 0.3;
    let flows = [C64::new(2e-4, -5e-5), C64::new(-8e-5, 4e-5)];
    let slopes = [C64::new(1e-7, 2e-8), C64::new(-2e-8, -4e-8)];
    let points = [[0.8, 0.1, 0.3], [-0.4, 0.7, 0.2]];
    let actual = material.aperture_response_with_frequency_derivative(
        &surface, &apertures, omega, Formulation::BurtonMiller, &limits(n), || false,
    ).unwrap();
    let original = material.aperture_response(
        &surface, &apertures, omega, Formulation::BurtonMiller, &limits(n), || false,
    ).unwrap();
    assert_eq!(actual.response().impedance().declaration, original.impedance().declaration);
    assert_eq!(actual.response().impedance().rows[0].impedance, original.impedance().rows[0].impedance);
    let observed = actual.observe_volume_flows(
        &surface, &flows, Some(&slopes), &points, 0.1, || false,
    ).unwrap();
    let original_fields = original.observe_volume_flows(&surface, &flows, &points, 0.1, || false).unwrap();
    let samples: Vec<_> = [-2.0, -1.0, 1.0, 2.0].into_iter().map(|offset| {
        let response = material.aperture_response(
            &surface, &apertures, omega + offset * h,
            Formulation::BurtonMiller, &limits(n), || false,
        ).unwrap();
        let flow: Vec<_> = flows.iter().zip(&slopes).map(|(&v, &s)| v + s.scale(offset * h)).collect();
        let fields = response.observe_volume_flows(&surface, &flow, &points, 0.1, || false).unwrap();
        let acoustic = response.impedance().acoustic(0).unwrap();
        (response, fields, acoustic)
    }).collect();
    let difference = |values: [C64; 4]| {
        (values[0] - values[1].scale(8.0) + values[2].scale(8.0) - values[3])
            .scale(1.0 / (12.0 * h))
    };
    let acoustic = actual.acoustic_derivative().unwrap();
    for i in 0..apertures.len() * apertures.len() {
        let expected = difference(std::array::from_fn(|j| samples[j].0.impedance().rows[0].impedance[i]));
        close(actual.impedance_derivative()[i], expected, 3e-8);
        let expected_acoustic = difference(std::array::from_fn(|j| samples[j].2[i]));
        close(acoustic[i], expected_acoustic, 3e-8);
    }
    for i in 0..points.len() {
        close(observed[i].field.pressure_pa, original_fields[i].pressure_pa, 2e-13);
        let expected = difference(std::array::from_fn(|j| samples[j].1[i].pressure_pa));
        close(observed[i].pressure_derivative_pa_per_rad_s, expected, 3e-7);
        for axis in 0..3 {
            let expected = difference(std::array::from_fn(|j| samples[j].1[i].gradient_pa_per_m[axis]));
            close(observed[i].gradient_derivative_pa_per_m_per_rad_s[axis], expected, 3e-7);
        }
    }
}

#[test]
fn zero_crossing_volume_flow_retains_the_full_pressure_slope_and_fixed_flows_agree() {
    let material = model(true);
    let surface = SpherePanels::icosphere(0.2, 0).unwrap();
    let n = surface.areas().len();
    let apertures = openings(&surface);
    let response = material.aperture_response_with_frequency_derivative(
        &surface, &apertures, 700.0, Formulation::BurtonMiller, &limits(n), || false,
    ).unwrap();
    let flows = [C64::new(2e-4, 5e-5), C64::new(-8e-5, 4e-5)];
    let zero = [C64::ZERO; 2];
    let points = [[0.8, 0.1, 0.3]];
    let crossing = response.observe_volume_flows(
        &surface, &zero, Some(&flows), &points, 0.1, || false,
    ).unwrap()[0];
    let reference = response.response().observe_volume_flows(
        &surface, &flows, &points, 0.1, || false,
    ).unwrap()[0];
    assert_eq!(crossing.field.pressure_pa, C64::ZERO);
    close(crossing.pressure_derivative_pa_per_rad_s, reference.pressure_pa, 2e-13);
    for axis in 0..3 {
        close(crossing.gradient_derivative_pa_per_m_per_rad_s[axis], reference.gradient_pa_per_m[axis], 2e-13);
    }
    assert!(crossing.pressure_group_delay_s().is_err());
    let held = response.observe_volume_flows(&surface, &flows, None, &points, 0.1, || false).unwrap()[0];
    let explicit = response.observe_volume_flows(&surface, &flows, Some(&zero), &points, 0.1, || false).unwrap()[0];
    assert_eq!(held.pressure_derivative_pa_per_rad_s, explicit.pressure_derivative_pa_per_rad_s);
    assert!(response.observe_volume_flows(&surface, &flows, Some(&zero[..1]), &points, 0.1, || false).is_err());
    assert!(response.observe_volume_flows(&surface, &flows, None, &points, 0.1, || true).is_err());
}

#[test]
fn zero_net_flow_keeps_generalized_load_but_refuses_plane_port_derivatives() {
    let material = model(true);
    let surface = SpherePanels::icosphere(0.2, 0).unwrap();
    let n = surface.areas().len();
    let mut shape = vec![0.0; n];
    // The two products cancel exactly, independent of unequal panel areas.
    shape[0] = surface.areas()[1];
    shape[1] = -surface.areas()[0];
    let apertures = [Aperture { shape, reference_point: [0.0; 3], normal: [0.0, 0.0, 1.0] }];
    let response = material.aperture_response_with_frequency_derivative(
        &surface, &apertures, 700.0, Formulation::BurtonMiller, &limits(n), || false,
    ).unwrap();
    assert!(response.impedance_derivative().iter().all(|&z| finite(z)));
    assert_eq!(response.response().impedance().volume_gains[0], 0.0);
    assert!(response.acoustic_derivative().is_err());
    assert!(response.observe_volume_flows(
        &surface, &[C64::ONE], None, &[[0.8, 0.1, 0.3]], 0.1, || false,
    ).is_err());
    assert!(material.aperture_response_with_frequency_derivative(
        &surface, &[], 700.0, Formulation::BurtonMiller, &limits(n), || false,
    ).is_err());
}
