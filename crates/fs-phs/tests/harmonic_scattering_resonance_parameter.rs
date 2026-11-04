//! Actual loaded peak tuning, including boundary feedback and wave interference.
use fs_math::c64::C64;
use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
use fs_phs::harmonic::resonance::{DrivenResonanceError, ResonanceControls};
use fs_phs::harmonic::scattering::{ScatteringBudget, ScatteringModel, ScatteringObservation};
use fs_phs::harmonic::{
    HarmonicBudget, HarmonicError, HarmonicModel, HarmonicParameterDirection as Direction,
};
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch};

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 32,
        max_dense_bytes: 2_000_000,
        max_dense_work: 4_000_000,
        relative_residual_tolerance: 2e-12,
    }
}
fn wave_budget() -> ScatteringBudget {
    ScatteringBudget {
        maximum_ports: 8,
        maximum_payload_bytes: 500_000,
        maximum_work: 1_000_000,
    }
}
fn controls(width: bool) -> ResonanceControls {
    ResonanceControls {
        absolute_frequency_tolerance_rad_s: 1e-11,
        relative_frequency_tolerance: 0.0,
        maximum_evaluations: 256,
        measure_half_power_bandwidth: width,
    }
}
fn loaded(source: &PortHamiltonian, ports: &[(usize, f64)]) -> ScatteringModel {
    ScatteringModel::new(
        HarmonicModel::from_ode(source, budget(), || false).unwrap(),
        ports,
        wave_budget(),
        || false,
    )
    .unwrap()
}
fn oscillator(m: f64, k: f64, c: f64, g: f64) -> PortHamiltonian {
    PortHamiltonian::new(
        2,
        1,
        vec![0., 1., -1., 0.],
        vec![0., 0., 0., c],
        vec![0., g],
        Box::new(QuadraticStorage::new(vec![k, 0., 0., 1. / m], 2).unwrap()),
    )
    .unwrap()
}
fn near(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + 4e-8 * expected.abs(),
        "{actual:e} != {expected:e}"
    );
}
fn complex_near(actual: C64, expected: C64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + 4e-8 * expected.abs(),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn loaded_mass_stiffness_loss_coupling_and_reference_move_the_exact_peak_and_width() {
    let (m, k, g, z) = (2.0_f64, 200.0_f64, 0.8_f64, 0.4_f64);
    let a = C64::new(0.7, -0.3);
    for c in [0.0, 4.0] {
        let model = loaded(&oscillator(m, k, c, g), &[(0, z)]);
        let gamma = c + g * g / z;
        let omega = (k / m).sqrt();
        let width = gamma / m;
        let quality = omega / width;
        let root = (gamma * gamma + 4. * m * k).sqrt();
        let lower = (root - gamma) / (2. * m);
        let upper = (root + gamma) / (2. * m);
        let y = a.scale(2. * g * g / (z * c + g * g));
        for (dm, dk, dc, dg, dz) in [
            (1., 0., 0., 0., 0.),
            (0., 1., 0., 0., 0.),
            (0., 0., 1., 0., 0.),
            (0., 0., 0., 1., 0.),
            (0., 0., 0., 0., 1.),
            (0.2, -0.3, 0.1, -0.4, 0.2),
        ] {
            let dq = [dk, 0., 0., -dm / (m * m)];
            let dr = [0., 0., 0., dc];
            let dp = [0., dg];
            let result = model
                .driven_resonance_with_parameter_derivative(
                    &[a],
                    ScatteringObservation::Output(&[1.]),
                    Direction {
                        storage: Some(&dq),
                        dissipation: Some(&dr),
                        ports: Some(&dp),
                        ..Default::default()
                    },
                    None,
                    Some(&[dz]),
                    (5., 16.),
                    controls(true),
                    || false,
                )
                .unwrap();
            let d_gamma = dc + 2. * g * dg / z - g * g * dz / (z * z);
            let d_omega = omega * 0.5 * (dk / k - dm / m);
            let d_width = d_gamma / m - gamma * dm / (m * m);
            let d_root = (gamma * d_gamma + 2. * (k * dm + m * dk)) / root;
            let d_log = 2. * dg / g - (c * dz + z * dc + 2. * g * dg) / (z * c + g * g);
            near(result.resonance.peak.angular_frequency_rad_s, omega, 1e-11);
            near(result.peak_frequency_derivative_rad_s, d_omega, 2e-10);
            complex_near(result.resonance.observation, y, 1e-10);
            complex_near(result.observation_derivative, y.scale(d_log), 2e-10);
            near(result.log_amplitude_derivative, d_log, 2e-10);
            near(result.phase_derivative, 0., 2e-10);
            let band = result.half_power_derivative.unwrap();
            near(
                band.lower_frequency_derivative_rad_s,
                (d_root - d_gamma) / (2. * m) - lower * dm / m,
                2e-10,
            );
            near(
                band.upper_frequency_derivative_rad_s,
                (d_root + d_gamma) / (2. * m) - upper * dm / m,
                2e-10,
            );
            near(band.bandwidth_derivative_rad_s, d_width, 2e-10);
            near(
                band.response_quality_factor_derivative,
                quality * (d_omega / omega - d_width / width),
                2e-9,
            );
            assert_eq!(result.evaluations, result.resonance.evaluations + 3);
            assert!(result.evaluations <= controls(true).maximum_evaluations);
            assert!(result.maximum_relative_residual <= 2e-12);
        }
    }
}

#[test]
fn displacement_and_storage_effort_follow_the_damped_loaded_peak() {
    let (m, k, c, g, z) = (2.0_f64, 200.0_f64, 4.0_f64, 0.8_f64, 0.4_f64);
    let (dm, dk, dc, dg, dz) = (0.2, 0.7, -0.1, 0.3, -0.2);
    let model = loaded(&oscillator(m, k, c, g), &[(0, z)]);
    let a = C64::new(0.7, -0.3);
    let gamma = c + g * g / z;
    let d_gamma = dc + 2. * g * dg / z - g * g * dz / (z * z);
    let omega = (k / m - gamma * gamma / (2. * m * m)).sqrt();
    let d_omega = (dk / m - k * dm / (m * m) - gamma * d_gamma / (m * m)
        + gamma * gamma * dm / (m * m * m))
        / (2. * omega);
    let dynamic = C64::new(k - m * omega * omega, -gamma * omega);
    let q = a.scale(2. * g / z) / dynamic;
    let d_dynamic = C64::new(
        dk - dm * omega * omega - 2. * m * omega * d_omega,
        -d_gamma * omega - gamma * d_omega,
    );
    let dq = (a.scale(2. * dg / z - 2. * g * dz / (z * z)) - d_dynamic * q) / dynamic;
    let dq_matrix = [dk, 0., 0., -dm / (m * m)];
    let dr = [0., 0., 0., dc];
    let dp = [0., dg];
    for (observation, expected, slope) in [
        (ScatteringObservation::State(&[1., 0.]), q, dq),
        (
            ScatteringObservation::Effort(&[1., 0.]),
            q.scale(k),
            dq.scale(k) + q.scale(dk),
        ),
    ] {
        let result = model
            .driven_resonance_with_parameter_derivative(
                &[a],
                observation,
                Direction {
                    storage: Some(&dq_matrix),
                    dissipation: Some(&dr),
                    ports: Some(&dp),
                    ..Default::default()
                },
                None,
                Some(&[dz]),
                (5., 16.),
                controls(false),
                || false,
            )
            .unwrap();
        near(result.resonance.peak.angular_frequency_rad_s, omega, 1e-11);
        near(result.peak_frequency_derivative_rad_s, d_omega, 2e-10);
        complex_near(result.resonance.observation, expected, 2e-10);
        complex_near(result.observation_derivative, slope, 2e-10);
        near(result.phase_derivative, (slope / expected).im, 2e-10);
    }
}

#[test]
fn changing_both_port_maps_and_references_differentiates_actual_force_and_reflection() {
    let a = [C64::ONE, C64::new(0., 1.)];
    let da = [C64::new(0.03, 0.04), C64::new(-0.02, 0.06)];
    let prepare = |theta: f64| {
        let source = PortHamiltonian::new(
            2,
            2,
            vec![0., 1., -1., 0.],
            vec![0., 0., 0., 0.2 + 0.08 * theta],
            vec![0., 0., 1. + 0.15 * theta, 1. - 0.11 * theta],
            Box::new(
                QuadraticStorage::new(
                    vec![100. - 0.6 * theta, 0., 0., 1. / (1. + 0.12 * theta)],
                    2,
                )
                .unwrap(),
            ),
        )
        .unwrap();
        // Wave order deliberately reverses source order. Both references change.
        loaded(&source, &[(1, 1. + 0.17 * theta), (0, 1. - 0.23 * theta)])
    };
    let model = prepare(0.);
    for observation in [
        ScatteringObservation::Input(&[1., 0.]),
        ScatteringObservation::Outgoing(&[0., 1.]),
    ] {
        let result = model
            .driven_resonance_with_parameter_derivative(
                &a,
                observation,
                Direction {
                    storage: Some(&[-0.6, 0., 0., -0.12]),
                    dissipation: Some(&[0., 0., 0., 0.08]),
                    ports: Some(&[0., 0., 0.15, -0.11]),
                    ..Default::default()
                },
                Some(&da),
                Some(&[0.17, -0.23]),
                (10., 16.),
                controls(false),
                || false,
            )
            .unwrap();
        // Independent physical reassembly and peak search are only a test
        // oracle. Both observations include incident feedthrough and therefore
        // peak away from the internal velocity resonance at omega=10.
        let step = 1e-4;
        let shifted = |theta| {
            prepare(theta)
                .driven_resonance(
                    &[a[0] + da[0].scale(theta), a[1] + da[1].scale(theta)],
                    observation,
                    (10., 16.),
                    controls(false),
                    || false,
                )
                .unwrap()
        };
        let lower = shifted(-step);
        let upper = shifted(step);
        let slope = (upper.observation - lower.observation).scale(0.5 / step);
        near(
            result.peak_frequency_derivative_rad_s,
            (upper.peak.angular_frequency_rad_s - lower.peak.angular_frequency_rad_s) / (2. * step),
            2e-6,
        );
        complex_near(result.observation_derivative, slope, 2e-6);
        near(
            result.log_amplitude_derivative,
            (slope / result.resonance.observation).re,
            2e-6,
        );
        near(
            result.phase_derivative,
            (slope / result.resonance.observation).im,
            2e-6,
        );
        assert!(result.resonance.peak.angular_frequency_rad_s > 11.);
        assert_eq!(result.evaluations, result.resonance.evaluations + 1);
        assert!(result.maximum_relative_residual <= 2e-12);
    }
}

#[test]
fn incident_gain_and_phase_only_move_the_peak_amplitude_at_extreme_source_scales() {
    let model = loaded(&oscillator(2., 200., 4., 0.8), &[(0, 0.4)]);
    let alpha = C64::new(-0.2, 0.3);
    for scale in [1e-80, 1., 1e80] {
        let a = C64::new(0.7, -0.3).scale(scale);
        let result = model
            .driven_resonance_with_parameter_derivative(
                &[a],
                ScatteringObservation::Output(&[1.]),
                Direction::default(),
                Some(&[alpha * a]),
                None,
                (5., 16.),
                controls(true),
                || false,
            )
            .unwrap();
        near(result.peak_frequency_derivative_rad_s, 0., 1e-10);
        complex_near(
            result.observation_derivative.scale(1. / scale),
            (result.resonance.observation * alpha).scale(1. / scale),
            1e-10,
        );
        near(result.log_amplitude_derivative, alpha.re, 1e-10);
        near(result.phase_derivative, alpha.im, 1e-10);
        let band = result.half_power_derivative.unwrap();
        near(band.lower_frequency_derivative_rad_s, 0., 1e-10);
        near(band.upper_frequency_derivative_rad_s, 0., 1e-10);
        near(band.bandwidth_derivative_rad_s, 0., 1e-10);
        near(band.response_quality_factor_derivative, 0., 1e-9);
    }
}

#[test]
fn loaded_descriptor_circuit_tunes_width_with_original_resistor_and_wave_reference() {
    let (r, l, c, z, dz) = (0.2_f64, 0.5_f64, 0.02_f64, 2.0_f64, 0.3_f64);
    let circuit = assemble_circuit(&CircuitGraph {
        node_count: 4,
        transformers: vec![],
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Inductor { henries: l }),
            (2, 3, Branch::Capacitor { farads: c }),
            (3, 0, Branch::Resistor { ohms: r }),
        ],
    })
    .unwrap();
    let dr: Vec<_> = circuit
        .system
        .structure()
        .1
        .iter()
        .map(|v| -v / r)
        .collect();
    let model = ScatteringModel::new(
        HarmonicModel::from_descriptor(&circuit.system, budget(), || false).unwrap(),
        &[(0, z)],
        wave_budget(),
        || false,
    )
    .unwrap();
    let a = C64::new(0.8, -0.2);
    let result = model
        .driven_resonance_with_parameter_derivative(
            &[a],
            ScatteringObservation::Output(&[1.]),
            Direction {
                dissipation: Some(&dr),
                ..Default::default()
            },
            None,
            Some(&[dz]),
            (6., 15.),
            controls(true),
            || false,
        )
        .unwrap();
    let width = (r + 1. / z) / l;
    let d_width = (1. - dz / (z * z)) / l;
    near(
        result.resonance.peak.angular_frequency_rad_s,
        1. / (l * c).sqrt(),
        1e-10,
    );
    near(result.peak_frequency_derivative_rad_s, 0., 1e-9);
    complex_near(
        result.observation_derivative,
        -a.scale(2. * (z + r * dz) / (1. + z * r).powi(2)),
        1e-9,
    );
    let band = result.half_power_derivative.unwrap();
    near(band.bandwidth_derivative_rad_s, d_width, 1e-9);
    near(
        band.response_quality_factor_derivative,
        -10. * d_width / (width * width),
        1e-8,
    );
    assert!(result.maximum_relative_residual <= 2e-12);
}

#[test]
fn loaded_maxwell_relaxation_keeps_the_memory_state_in_peak_and_bandwidth_tuning() {
    let tau = 0.05_f64;
    let prepare = |time| {
        let source = fs_phs::mass_spring_damper(1., 100., 1.)
            .unwrap()
            .with_relaxation_branches(vec![RelaxationBranch {
                projection: vec![1., 0.],
                stiffness: 50.,
                relaxation_time_s: time,
            }])
            .unwrap();
        loaded(&source, &[(0, 2.)])
    };
    let model = prepare(tau);
    let mut dr = [0.; 9];
    dr[8] = -1. / (tau * tau);
    let input = [C64::new(0.8, -0.2)];
    let result = model
        .driven_resonance_with_parameter_derivative(
            &input,
            ScatteringObservation::Output(&[1.]),
            Direction {
                dissipation: Some(&dr),
                ..Default::default()
            },
            None,
            None,
            (5., 20.),
            controls(true),
            || false,
        )
        .unwrap();
    let step = 1e-5;
    let shifted = |time| {
        prepare(time)
            .driven_resonance(
                &input,
                ScatteringObservation::Output(&[1.]),
                (5., 20.),
                controls(true),
                || false,
            )
            .unwrap()
    };
    let lower = shifted(tau - step);
    let upper = shifted(tau + step);
    near(
        result.peak_frequency_derivative_rad_s,
        (upper.peak.angular_frequency_rad_s - lower.peak.angular_frequency_rad_s) / (2. * step),
        3e-6,
    );
    complex_near(
        result.observation_derivative,
        (upper.observation - lower.observation).scale(0.5 / step),
        3e-6,
    );
    let band = result.half_power_derivative.unwrap();
    near(
        band.bandwidth_derivative_rad_s,
        (upper.half_power.unwrap().bandwidth_rad_s - lower.half_power.unwrap().bandwidth_rad_s)
            / (2. * step),
        5e-6,
    );
    near(
        band.response_quality_factor_derivative,
        (upper.half_power.unwrap().response_quality_factor
            - lower.half_power.unwrap().response_quality_factor)
            / (2. * step),
        5e-6,
    );
    assert!(result.maximum_relative_residual <= 2e-12);
}

#[test]
fn loaded_peak_tuning_refuses_incomplete_directions_caps_and_cancelled_work() {
    let source = oscillator(1., 100., 0.2, 1.);
    let model = loaded(&source, &[(0, 1.)]);
    let compute = |controls, source_slope: Option<&[C64]>, dz: Option<&[f64]>| {
        model.driven_resonance_with_parameter_derivative(
            &[C64::ONE],
            ScatteringObservation::Output(&[1.]),
            Direction::default(),
            source_slope,
            dz,
            (7., 13.),
            controls,
            || false,
        )
    };
    assert!(matches!(
        compute(controls(true), Some(&[]), None),
        Err(DrivenResonanceError::Invalid(_))
    ));
    assert!(matches!(
        compute(controls(true), None, Some(&[])),
        Err(DrivenResonanceError::Invalid(_))
    ));
    assert!(matches!(
        compute(controls(true), None, Some(&[f64::NAN])),
        Err(DrivenResonanceError::Invalid(_))
    ));
    let mut cap = controls(true);
    cap.maximum_evaluations = 5;
    assert!(matches!(
        compute(cap, None, None),
        Err(DrivenResonanceError::Invalid(_))
    ));
    // The constructor and individual tangent queries fit these wave budgets;
    // the additional simultaneous peak arrays/work must still be admitted.
    for (bytes, work) in [(2500, 1_000_000), (500_000, 600)] {
        let limited = ScatteringModel::new(
            HarmonicModel::from_ode(&source, budget(), || false).unwrap(),
            &[(0, 1.)],
            ScatteringBudget {
                maximum_ports: 1,
                maximum_payload_bytes: bytes,
                maximum_work: work,
            },
            || false,
        )
        .unwrap();
        assert!(matches!(
            limited.driven_resonance_with_parameter_derivative(
                &[C64::ONE],
                ScatteringObservation::Output(&[1.]),
                Direction::default(),
                None,
                None,
                (7., 13.),
                controls(true),
                || false
            ),
            Err(DrivenResonanceError::Invalid(_))
        ));
    }
    let mut polls = 0;
    assert_eq!(
        model.driven_resonance_with_parameter_derivative(
            &[C64::ONE],
            ScatteringObservation::Output(&[1.]),
            Direction::default(),
            None,
            None,
            (7., 13.),
            controls(true),
            || {
                polls += 1;
                polls > 120
            }
        ),
        Err(DrivenResonanceError::Harmonic(HarmonicError::Cancelled))
    );
    let first = compute(controls(true), None, None).unwrap();
    assert_eq!(first, compute(controls(true), None, None).unwrap());
    assert_eq!(first.peak_frequency_derivative_rad_s, 0.);
    assert_eq!(first.observation_derivative, C64::ZERO);
}
