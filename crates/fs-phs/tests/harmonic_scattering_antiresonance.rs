//! Coherent reflection notches retain physical finite-load work, energy and heat.
use fs_math::c64::C64;
use fs_phs::{
    PortHamiltonian, QuadraticStorage,
    harmonic::{
        HarmonicBudget, HarmonicError, HarmonicModel,
        resonance::{AntiresonanceControls, DrivenResonanceError},
        scattering::{ScatteringBudget, ScatteringModel, ScatteringObservation},
    },
};

fn controls(level: Option<f64>) -> AntiresonanceControls {
    AntiresonanceControls {
        absolute_frequency_tolerance_rad_s: 1e-10,
        relative_frequency_tolerance: 0.0,
        maximum_evaluations: 256,
        recovery_power: level,
    }
}

fn prepare(system: &PortHamiltonian, ports: &[(usize, f64)]) -> ScatteringModel {
    ScatteringModel::new(
        HarmonicModel::from_ode(
            system,
            HarmonicBudget {
                max_dimension: 8,
                max_dense_bytes: 100_000,
                max_dense_work: 100_000,
                relative_residual_tolerance: 2e-12,
            },
            || false,
        )
        .unwrap(),
        ports,
        ScatteringBudget {
            maximum_ports: 2,
            maximum_payload_bytes: 100_000,
            maximum_work: 100_000,
        },
        || false,
    )
    .unwrap()
}

#[test]
fn critically_matched_reflection_zero_keeps_absorbed_work_heat_and_stored_energy() {
    let source = fs_phs::mass_spring_damper(1.0, 100.0, 1.0).unwrap();
    let model = prepare(&source, &[(0, 1.0)]);
    let input = [C64::new(0.6, -0.8)];
    let result = model
        .driven_antiresonance(
            &input,
            ScatteringObservation::Outgoing(&[1.0]),
            (7.0, 13.0),
            controls(Some(0.25)),
            || false,
        )
        .unwrap();
    assert!((result.minimum.angular_frequency_rad_s - 10.0).abs() < 1e-10);
    assert!(result.observation.abs() < 1e-10);
    assert!((result.observation_derivative - C64::new(0.0, 1.0) * input[0]).abs() < 1e-9);
    // b/a=-i*x/(2+i*x), x=100/omega-omega. Recovery at |b|^2=1/4
    // occurs at |x|=2/sqrt(3), with unequal frequency offsets around 10.
    let x = 2.0 / 3.0_f64.sqrt();
    let root = (x * x + 400.0).sqrt();
    let band = result.recovery.unwrap();
    assert!((band.lower.angular_frequency_rad_s - (root - x) / 2.0).abs() < 1e-10);
    assert!((band.upper.angular_frequency_rad_s - (root + x) / 2.0).abs() < 1e-10);
    assert!((band.bandwidth_rad_s - x).abs() < 2e-10);
    let response = model
        .factor(result.minimum.angular_frequency_rad_s, || false)
        .unwrap()
        .response(&input, || false)
        .unwrap();
    assert!((response.physical.output[0] - input[0]).abs() < 1e-10);
    assert!((response.physical_inputs[0] - input[0]).abs() < 1e-10);
    assert!((response.physical.mean_stored_energy_j - 0.5).abs() < 1e-10);
    assert!((response.physical.dissipated_power_w - 0.5).abs() < 1e-10);
    assert!((response.physical.supplied_power_w - 0.5).abs() < 1e-10);
    assert!(response.waves[0].outgoing_power_w < 1e-20);
    assert!((response.waves[0].incident_power_w - 0.5).abs() < 1e-10);
    assert!(response.wave_power_residual_w.abs() < 1e-12);
    assert!(result.maximum_relative_residual < 2e-12);
}

#[test]
fn unmatched_dissipation_gives_a_finite_notch_with_its_own_absolute_recovery_width() {
    for damping in [0.2_f64, 3.0_f64] {
        let model = prepare(
            &fs_phs::mass_spring_damper(1.0, 100.0, damping).unwrap(),
            &[(0, 1.0)],
        );
        let input = [C64::new(-0.8, -0.6)];
        let result = model
            .driven_antiresonance(
                &input,
                ScatteringObservation::Outgoing(&[1.0]),
                (5.0, 16.0),
                controls(Some(0.75)),
                || false,
            )
            .unwrap();
        assert!((result.minimum.angular_frequency_rad_s - 10.0).abs() < 1e-10);
        let expected = input[0].scale((1.0 - damping) / (1.0 + damping));
        assert!((result.observation - expected).abs() < 1e-10);
        assert!(result.observation.abs() > 0.49);
        let x = ((0.75 * (1.0 + damping).powi(2) - (1.0 - damping).powi(2)) / 0.25).sqrt();
        let root = (x * x + 400.0).sqrt();
        let band = result.recovery.unwrap();
        assert!((band.lower.angular_frequency_rad_s - (root - x) / 2.0).abs() < 1e-10);
        assert!((band.upper.angular_frequency_rad_s - (root + x) / 2.0).abs() < 1e-10);
    }
}

#[test]
fn coherent_multiport_reflection_and_actual_reaction_minima_retain_source_order() {
    let source = PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, 0.2],
        vec![0.0, 0.0, 1.0, 1.0],
        Box::new(QuadraticStorage::new(vec![100.0, 0.0, 0.0, 1.0], 2).unwrap()),
    )
    .unwrap();
    // Outgoing uses declared wave order; reaction inputs retain source order.
    let model = prepare(&source, &[(1, 1.0), (0, 1.0)]);
    let phase = C64::new(0.3, -0.7);
    let input = [phase, phase * C64::new(0.0, 1.0)];
    let gamma = 2.2_f64;
    for (observation, coefficient, reaction) in [
        (ScatteringObservation::Outgoing(&[0.0, 1.0]), 0.4, false),
        (ScatteringObservation::Input(&[1.0, 0.0]), 2.4, true),
    ] {
        let x = (-coefficient + (coefficient * coefficient + 4.0 * gamma * gamma).sqrt()) / 2.0;
        let expected = ((x * x + 400.0).sqrt() - x) / 2.0;
        let result = model
            .driven_antiresonance(&input, observation, (6.0, 10.0), controls(None), || false)
            .unwrap();
        assert!((result.minimum.angular_frequency_rad_s - expected).abs() < 1e-10);
        assert!(result.minimum.angular_frequency_rad_s < 9.5);
        let velocity = phase * (C64::new(2.0, 2.0) / C64::new(gamma, x));
        let physical = if reaction {
            phase.scale(2.0) - velocity
        } else {
            velocity - phase
        };
        assert!((result.observation - physical).abs() < 1e-10);
        assert!(result.maximum_relative_residual < 2e-12);
    }
}

#[test]
fn invalid_order_missing_recovery_and_cancellation_refuse_without_changing_replay() {
    let model = prepare(
        &fs_phs::mass_spring_damper(1.0, 100.0, 1.0).unwrap(),
        &[(0, 1.0)],
    );
    let input = [C64::ONE];
    assert!(matches!(
        model.driven_antiresonance(
            &input,
            ScatteringObservation::Outgoing(&[1.0, 0.0]),
            (7.0, 13.0),
            controls(None),
            || false,
        ),
        Err(DrivenResonanceError::Invalid(_))
    ));
    assert!(matches!(
        model.driven_antiresonance(
            &input,
            ScatteringObservation::Outgoing(&[1.0]),
            (9.9, 10.1),
            controls(Some(0.25)),
            || false,
        ),
        Err(DrivenResonanceError::NotBracketed(_))
    ));
    let mut checks = 0;
    assert_eq!(
        model.driven_antiresonance(
            &input,
            ScatteringObservation::Outgoing(&[1.0]),
            (7.0, 13.0),
            controls(Some(0.25)),
            || {
                checks += 1;
                checks > 120
            },
        ),
        Err(DrivenResonanceError::Harmonic(HarmonicError::Cancelled))
    );
    let run = || {
        model
            .driven_antiresonance(
                &input,
                ScatteringObservation::Outgoing(&[1.0]),
                (7.0, 13.0),
                controls(Some(0.25)),
                || false,
            )
            .unwrap()
    };
    assert_eq!(run(), run());
}

#[test]
fn a_lossless_single_port_allpass_has_no_reflected_amplitude_notch() {
    let model = prepare(
        &fs_phs::mass_spring_damper(1.0, 100.0, 0.0).unwrap(),
        &[(0, 1.0)],
    );
    for input in [C64::ONE, C64::new(0.6, -0.8), C64::new(-0.3, 0.7)] {
        for interval in [(7.0, 13.0), (8.0, 12.0), (5.0, 16.0), (1.0, 100.0)] {
            let result = model.driven_antiresonance(
                &[input],
                ScatteringObservation::Outgoing(&[1.0]),
                interval,
                controls(None),
                || false,
            );
            assert!(
                result.is_err(),
                "a lossless allpass became a notch: {result:?}"
            );
        }
    }
}
