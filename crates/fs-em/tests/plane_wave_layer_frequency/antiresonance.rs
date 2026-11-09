//! Public coherent-notch contracts checked against independent film formulas.
use super::*;
use fs_em::plane_wave::layers::AntiresonanceControls;

fn notch_controls(recovery_power: Option<f64>) -> AntiresonanceControls {
    AntiresonanceControls {
        absolute_frequency_tolerance_rad_s: 0.1,
        relative_frequency_tolerance: 0.0,
        maximum_evaluations: 160,
        recovery_power,
    }
}

#[test]
fn symmetric_lossless_film_reflection_notch_and_recovery_follow_airy_law() {
    let ambient = material(E, M);
    let exit = DielectricResponse::Lossless(&ambient);
    let beta = (E * M).sqrt();
    let level = 0.1;
    for index in [2.0_f64, 3.0] {
        let film = material(E * index * index, M);
        let thickness = core::f64::consts::PI / (index * W * beta);
        let films = [DielectricLayer {
            response: DielectricResponse::Lossless(&film),
            thickness_m: thickness,
        }];
        let stack = PlanarDielectricLayers::try_new(&ambient, exit, W, &films).unwrap();
        // Independent Airy reflectance for equal lossless exterior media:
        // R = F sin²(delta)/(1 + F sin²(delta)), delta = omega*n*sqrt(E*M)*d.
        // The first nontrivial reflection zero is delta=pi, with unit flux
        // transmission rather than a vanishing physical field everywhere.
        let airy = (index * index - 1.0).powi(2) / (4.0 * index * index);
        let phase_slope = index * beta * thickness;
        let phase_edge = (level / (airy * (1.0 - level))).sqrt().asin();
        let edge_offset = phase_edge / phase_slope;
        for polarization in [Polarization::Te, Polarization::Tm] {
            for ratio in [0.8, 0.9, 1.05, 1.2] {
                let omega = ratio * W;
                let sine_squared = (omega * phase_slope).sin().powi(2);
                let expected = airy * sine_squared / (1.0 + airy * sine_squared);
                let actual = PlanarDielectricLayers::try_new(&ambient, exit, omega, &films)
                    .unwrap()
                    .response(1.0, polarization)
                    .unwrap();
                assert!((actual.reflected_power_fraction - expected).abs() < 2e-13);
            }
            let result = stack
                .driven_antiresonance(
                    1.0,
                    polarization,
                    Constraint::FixedAngle,
                    Observation::ReflectedElectricAmplitude,
                    (0.65 * W, 1.35 * W),
                    notch_controls(Some(level)),
                    || false,
                )
                .unwrap();
            assert!((result.minimum.angular_frequency_rad_s - W).abs() < 0.2);
            assert!(result.observation.abs() < 2e-9);
            // At the coherent zero |dr/domega| = sqrt(F)*d(delta)/domega.
            // This tests a nonzero analytic field derivative through zero;
            // neither sqrt(reflectance) nor a logarithmic phase supplies it.
            let expected_derivative = airy.sqrt() * phase_slope;
            assert!((result.observation_derivative.abs() - expected_derivative).abs() < 5e-21);
            assert!(result.evaluations <= 160);
            assert!(result.maximum_relative_residual <= 512.0 * f64::EPSILON);
            let recovery = result.recovery.unwrap();
            assert_eq!(recovery.power_level, level);
            assert!((recovery.lower.angular_frequency_rad_s - (W - edge_offset)).abs() < 0.2);
            assert!((recovery.upper.angular_frequency_rad_s - (W + edge_offset)).abs() < 0.2);
            assert!((recovery.bandwidth_rad_s - 2.0 * edge_offset).abs() < 0.3);
            for edge in [recovery.lower, recovery.upper] {
                let phase = edge.angular_frequency_rad_s * phase_slope;
                let expected_power =
                    airy * phase.sin().powi(2) / (1.0 + airy * phase.sin().powi(2));
                assert!((expected_power - level).abs() < 2e-9);
            }
            let transparent = stack.response(1.0, polarization).unwrap();
            near(transparent.reflected_tangential_e, C64::ZERO, 2e-13);
            near(transparent.transmitted_tangential_e, -C64::ONE, 2e-13);
            near(transparent.transmitted_conjugate_h, -C64::ONE, 2e-13);
            assert!(transparent.transmitted_power_fraction.is_finite());
            assert!((transparent.transmitted_power_fraction - 1.0).abs() < 2e-13);
            assert_eq!(transparent.absorbed_power_fraction, 0.0);
        }
    }
}

#[test]
fn power_cusps_missing_recovery_edges_and_cancellation_refuse() {
    let ambient = material(E, M);
    let film = material(4.0 * E, M);
    let films = [DielectricLayer {
        response: DielectricResponse::Lossless(&film),
        thickness_m: core::f64::consts::PI / (2.0 * W * (E * M).sqrt()),
    }];
    let stack = PlanarDielectricLayers::try_new(
        &ambient,
        DielectricResponse::Lossless(&ambient),
        W,
        &films,
    )
    .unwrap();
    for observation in [
        Observation::TransmittedPowerFraction,
        Observation::AbsorbedPowerFraction,
        Observation::LayerAbsorbedPowerFraction { layer_index: 0 },
        Observation::AbsorbedPowerDensity,
        Observation::LayerAbsorbedPowerDensity { layer_index: 0 },
    ] {
        assert!(matches!(
            stack.driven_antiresonance(
                1.0,
                Polarization::Te,
                Constraint::FixedAngle,
                observation,
                (0.65 * W, 1.35 * W),
                notch_controls(None),
                || false,
            ),
            Err(ResonanceSearchError::Source(_))
        ));
    }
    // The reflection zero is bracketed, but R=0.1 is outside both endpoints.
    assert!(matches!(
        stack.driven_antiresonance(
            1.0,
            Polarization::Te,
            Constraint::FixedAngle,
            Observation::ReflectedElectricAmplitude,
            (0.98 * W, 1.02 * W),
            notch_controls(Some(0.1)),
            || false,
        ),
        Err(ResonanceSearchError::NotBracketed(_))
    ));
    // Both endpoint slopes lie on the decreasing side of the same notch.
    assert!(matches!(
        stack.driven_antiresonance(
            1.0,
            Polarization::Tm,
            Constraint::FixedAngle,
            Observation::ReflectedElectricAmplitude,
            (0.6 * W, 0.8 * W),
            notch_controls(None),
            || false,
        ),
        Err(ResonanceSearchError::NotBracketed(_))
    ));
    assert!(matches!(
        stack.driven_antiresonance(
            1.0,
            Polarization::Te,
            Constraint::FixedAngle,
            Observation::ReflectedElectricAmplitude,
            (0.65 * W, 1.35 * W),
            notch_controls(None),
            || true,
        ),
        Err(ResonanceSearchError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        stack.driven_antiresonance(
            1.0,
            Polarization::Te,
            Constraint::FixedAngle,
            Observation::ReflectedElectricAmplitude,
            (0.65 * W, 1.35 * W),
            notch_controls(None),
            || {
                polls += 1;
                polls >= 5
            },
        ),
        Err(ResonanceSearchError::Cancelled)
    ));
    assert!(polls >= 5);
}

#[test]
fn oblique_notch_keeps_the_original_fixed_transverse_wavenumber() {
    let ambient = material(E, M);
    let index = 2.0_f64;
    let film = material(index * index * E, M);
    let exit = DielectricResponse::Lossless(&ambient);
    let beta = (E * M).sqrt();
    let prepared_cosine = 0.8_f64;
    let transverse = W * beta * (1.0 - prepared_cosine * prepared_cosine).sqrt();
    // Deliberately put the physical notch away from the prepared carrier.
    // Holding its angle instead of k_parallel gives a different frequency.
    let notch = 1.15 * W;
    let film_normal = ((index * beta * notch).powi(2) - transverse * transverse).sqrt();
    let thickness = core::f64::consts::PI / film_normal;
    let films = [DielectricLayer {
        response: DielectricResponse::Lossless(&film),
        thickness_m: thickness,
    }];
    let stack = PlanarDielectricLayers::try_new(&ambient, exit, W, &films).unwrap();
    let incident_normal = ((beta * notch).powi(2) - transverse * transverse).sqrt();
    let notch_cosine = incident_normal / (beta * notch);
    // delta = d sqrt(n² E M omega² - k_parallel²), with derivative in seconds.
    let phase_slope = thickness * (index * beta).powi(2) * notch / film_normal;
    for polarization in [Polarization::Te, Polarization::Tm] {
        let result = stack
            .driven_antiresonance(
                prepared_cosine,
                polarization,
                Constraint::FixedTransverseWavenumber,
                Observation::ReflectedElectricAmplitude,
                (0.95 * notch, 1.05 * notch),
                notch_controls(None),
                || false,
            )
            .unwrap();
        assert!((result.minimum.angular_frequency_rad_s - notch).abs() < 0.2);
        assert!(result.observation.abs() < 2e-9);
        assert!(result.recovery.is_none());
        let admittance_ratio = match polarization {
            Polarization::Te => film_normal / incident_normal,
            Polarization::Tm => index * index * incident_normal / film_normal,
        };
        let airy_amplitude = 0.5 * (admittance_ratio - 1.0 / admittance_ratio).abs();
        let expected_derivative = airy_amplitude * phase_slope;
        assert!(
            (result.observation_derivative.abs() - expected_derivative).abs()
                < 5e-9 * expected_derivative
        );
        let transparent = PlanarDielectricLayers::try_new(&ambient, exit, notch, &films)
            .unwrap()
            .response(notch_cosine, polarization)
            .unwrap();
        near(transparent.reflected_tangential_e, C64::ZERO, 2e-13);
        assert!(transparent.transmitted_power_fraction.is_finite());
        assert!((transparent.transmitted_power_fraction - 1.0).abs() < 2e-13);
    }
}
