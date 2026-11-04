//! Local observed-peak/linewidth parameter derivatives against physical oracles.
use fs_math::c64::C64;
use fs_phs::RelaxationBranch;
use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
use fs_phs::harmonic::resonance::{DrivenResonanceError, HarmonicObservation, ResonanceControls};
use fs_phs::harmonic::{
    HarmonicBudget, HarmonicError, HarmonicModel, HarmonicParameterDirection as Direction,
};

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 32,
        max_dense_bytes: 2_000_000,
        max_dense_work: 4_000_000,
        relative_residual_tolerance: 2e-12,
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
fn close(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() < absolute + 3e-8 * expected.abs(),
        "{actual} != {expected}"
    );
}
fn complex_close(actual: C64, expected: C64, absolute: f64) {
    assert!(
        (actual - expected).abs() < absolute + 3e-8 * expected.abs(),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn velocity_peak_and_linewidth_differentiate_mass_stiffness_and_damping() {
    let (m, k, c) = (2.0_f64, 200.0_f64, 4.0_f64);
    let system = fs_phs::mass_spring_damper(m, k, c).unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    let dk = [1.0, 0.0, 0.0, 0.0];
    let dm = [0.0, 0.0, 0.0, -1.0 / (m * m)];
    let dc = [0.0, 0.0, 0.0, 1.0];
    let force = C64::new(0.8, -0.3);
    let peak = (k / m).sqrt();
    let width = c / m;
    let quality = peak / width;
    for (direction, dw, db, dq, dy) in [
        (
            Direction {
                storage: Some(&dk),
                ..Default::default()
            },
            1.0 / (2.0 * m * peak),
            0.0,
            quality / (2.0 * k),
            C64::ZERO,
        ),
        (
            Direction {
                storage: Some(&dm),
                ..Default::default()
            },
            -peak / (2.0 * m),
            -c / (m * m),
            quality / (2.0 * m),
            C64::ZERO,
        ),
        (
            Direction {
                dissipation: Some(&dc),
                ..Default::default()
            },
            0.0,
            1.0 / m,
            -quality / c,
            -force.scale(1.0 / (c * c)),
        ),
    ] {
        let result = model
            .driven_resonance_with_parameter_derivative(
                &[force],
                HarmonicObservation::Output(&[1.0]),
                direction,
                None,
                (5.0, 16.0),
                controls(true),
                || false,
            )
            .unwrap();
        close(result.resonance.peak.angular_frequency_rad_s, peak, 1e-11);
        close(result.peak_frequency_derivative_rad_s, dw, 2e-11);
        complex_close(result.observation_derivative, dy, 3e-11);
        let band = result.half_power_derivative.unwrap();
        close(band.bandwidth_derivative_rad_s, db, 1e-10);
        close(band.response_quality_factor_derivative, dq, 1e-9);
        assert!(result.log_amplitude_curvature < 0.0);
        assert!(result.maximum_relative_residual <= 2e-12);
        assert_eq!(result.evaluations, result.resonance.evaluations + 3);
        assert!(result.evaluations <= controls(true).maximum_evaluations);
    }
}

#[test]
fn narrow_resonance_curvature_preserves_a_high_q_linewidth_sensitivity() {
    let (mass, stiffness, damping) = (2.0_f64, 200.0_f64, 1e-4_f64);
    let source = fs_phs::mass_spring_damper(mass, stiffness, damping).unwrap();
    let model = HarmonicModel::from_ode(&source, budget(), || false).unwrap();
    let dc = [0.0, 0.0, 0.0, 1.0];
    let result = model
        .driven_resonance_with_parameter_derivative(
            &[C64::ONE],
            HarmonicObservation::Output(&[1.0]),
            Direction {
                dissipation: Some(&dc),
                ..Default::default()
            },
            None,
            (8.0, 12.0),
            controls(true),
            || false,
        )
        .unwrap();
    close(result.peak_frequency_derivative_rad_s, 0.0, 2e-10);
    close(result.log_amplitude_derivative, -1.0 / damping, 1e-6);
    let band = result.half_power_derivative.unwrap();
    close(band.bandwidth_derivative_rad_s, 1.0 / mass, 3e-7);
    let expected_quality_slope = -10.0 * mass / (damping * damping);
    assert!((band.response_quality_factor_derivative / expected_quality_slope - 1.0).abs() < 3e-6);
    assert!(result.evaluations < 190);
}

#[test]
fn displacement_peak_motion_is_not_the_natural_or_free_decay_frequency_motion() {
    let (m, k, c) = (2.0_f64, 200.0_f64, 8.0_f64);
    let system = fs_phs::mass_spring_damper(m, k, c).unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    let dc = [0.0, 0.0, 0.0, 1.0];
    let force = C64::new(0.8, -0.3);
    let result = model
        .driven_resonance_with_parameter_derivative(
            &[force],
            HarmonicObservation::State(&[1.0, 0.0]),
            Direction {
                dissipation: Some(&dc),
                ..Default::default()
            },
            None,
            (4.0, 15.0),
            controls(true),
            || false,
        )
        .unwrap();
    let peak = (k / m - c * c / (2.0 * m * m)).sqrt();
    let dp = -c / (2.0 * m * m * peak);
    close(result.peak_frequency_derivative_rad_s, dp, 3e-11);
    let free_decay = (k / m - c * c / (4.0 * m * m)).sqrt();
    assert!((dp - (-c / (4.0 * m * m * free_decay))).abs() > 0.04);
    let dynamic = C64::new(k - m * peak * peak, -c * peak);
    let dynamic_derivative = C64::new(-2.0 * m * peak, -c).scale(dp) + C64::new(0.0, -peak);
    complex_close(
        result.observation_derivative,
        -force * dynamic_derivative / (dynamic * dynamic),
        1e-12,
    );
    close(
        result.log_amplitude_derivative,
        -1.0 / c + c / (4.0 * m * k - c * c),
        1e-11,
    );
    let minimum = c * c * k / m - c.powi(4) / (4.0 * m * m);
    let d_minimum = 2.0 * c * k / m - c.powi(3) / (m * m);
    let lower = (peak * peak - minimum.sqrt() / m).sqrt();
    let upper = (peak * peak + minimum.sqrt() / m).sqrt();
    let dl = (-c / (m * m) - d_minimum / (2.0 * m * minimum.sqrt())) / (2.0 * lower);
    let du = (-c / (m * m) + d_minimum / (2.0 * m * minimum.sqrt())) / (2.0 * upper);
    let band = result.half_power_derivative.unwrap();
    close(band.lower_frequency_derivative_rad_s, dl, 2e-10);
    close(band.upper_frequency_derivative_rad_s, du, 2e-10);
    close(band.bandwidth_derivative_rad_s, du - dl, 3e-10);
}

#[test]
fn source_gain_phase_and_port_map_change_amplitude_without_detuning_the_peak() {
    let system = fs_phs::mass_spring_damper(1.0, 100.0, 0.5).unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    let dg = [0.0, 1.0];
    let alpha = C64::new(0.3, -0.4);
    for force in [
        C64::new(0.8, -0.2),
        C64::new(1e-60, -2e-60),
        C64::new(-1e60, 2e60),
    ] {
        let result = model
            .driven_resonance_with_parameter_derivative(
                &[force],
                HarmonicObservation::Output(&[-3.0]),
                Direction {
                    ports: Some(&dg),
                    ..Default::default()
                },
                Some(&[alpha * force]),
                (7.0, 13.0),
                controls(true),
                || false,
            )
            .unwrap();
        close(result.peak_frequency_derivative_rad_s, 0.0, 2e-10);
        // G changes both source and observation: d log y=2+alpha.
        close(result.log_amplitude_derivative, 2.0 + alpha.re, 2e-10);
        close(result.phase_derivative, alpha.im, 2e-10);
        complex_close(
            result.observation_derivative,
            result.resonance.observation * (C64::from_re(2.0) + alpha),
            force.abs() * 1e-9,
        );
        close(
            result
                .half_power_derivative
                .unwrap()
                .bandwidth_derivative_rad_s,
            0.0,
            1e-9,
        );
    }
}

#[test]
fn descriptor_rlc_resistance_moves_linewidth_without_detuning_the_current_peak() {
    let (r, l, c) = (0.02_f64, 0.25_f64, 0.04_f64);
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
    let model = HarmonicModel::from_descriptor(&circuit.system, budget(), || false).unwrap();
    // Circuit resistor contributes conductance 1/R in original algebraic rows.
    let dr: Vec<_> = circuit
        .system
        .structure()
        .1
        .iter()
        .map(|x| -x / r)
        .collect();
    let voltage = C64::new(2.0, -0.4);
    let result = model
        .driven_resonance_with_parameter_derivative(
            &[voltage],
            HarmonicObservation::Output(&[1.0]),
            Direction {
                dissipation: Some(&dr),
                ..Default::default()
            },
            None,
            (7.0, 14.0),
            controls(true),
            || false,
        )
        .unwrap();
    close(result.peak_frequency_derivative_rad_s, 0.0, 2e-9);
    complex_close(
        result.observation_derivative,
        -voltage.scale(1.0 / (r * r)),
        1e-6,
    );
    let band = result.half_power_derivative.unwrap();
    close(band.bandwidth_derivative_rad_s, 1.0 / l, 3e-8);
    close(
        band.response_quality_factor_derivative,
        -l / (r * r * (l * c).sqrt()),
        1e-5,
    );
}

#[test]
fn maxwell_relaxation_peak_derivative_retains_the_causal_internal_state() {
    let tau = 0.05_f64;
    let prepare = |time| {
        let source = fs_phs::mass_spring_damper(1.0, 100.0, 1.0)
            .unwrap()
            .with_relaxation_branches(vec![RelaxationBranch {
                projection: vec![1.0, 0.0],
                stiffness: 50.0,
                relaxation_time_s: time,
            }])
            .unwrap();
        HarmonicModel::from_ode(&source, budget(), || false).unwrap()
    };
    let model = prepare(tau);
    let mut dr = [0.0; 9];
    dr[8] = -1.0 / (tau * tau);
    let input = [C64::new(0.8, -0.2)];
    let result = model
        .driven_resonance_with_parameter_derivative(
            &input,
            HarmonicObservation::Output(&[1.0]),
            Direction {
                dissipation: Some(&dr),
                ..Default::default()
            },
            None,
            (6.0, 18.0),
            controls(false),
            || false,
        )
        .unwrap();
    // Independent parameter perturbation of the physical relaxation time checks
    // the mixed tangent and implicit-root chain rule; production uses neither.
    let step = 1e-5;
    let a = prepare(tau - step)
        .driven_resonance(
            &input,
            HarmonicObservation::Output(&[1.0]),
            (6.0, 18.0),
            controls(false),
            || false,
        )
        .unwrap();
    let b = prepare(tau + step)
        .driven_resonance(
            &input,
            HarmonicObservation::Output(&[1.0]),
            (6.0, 18.0),
            controls(false),
            || false,
        )
        .unwrap();
    close(
        result.peak_frequency_derivative_rad_s,
        (b.peak.angular_frequency_rad_s - a.peak.angular_frequency_rad_s) / (2.0 * step),
        2e-6,
    );
    complex_close(
        result.observation_derivative,
        (b.observation - a.observation).scale(0.5 / step),
        2e-6,
    );
    assert_eq!(result.half_power_derivative, None);
    assert_eq!(result.evaluations, result.resonance.evaluations + 1);
}

#[test]
fn refusal_and_cancellation_keep_resonance_parameter_queries_transactional() {
    let system = fs_phs::mass_spring_damper(1.0, 100.0, 0.5).unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    let mut cap = controls(true);
    cap.maximum_evaluations = 5;
    assert!(matches!(
        model.driven_resonance_with_parameter_derivative(
            &[C64::ONE],
            HarmonicObservation::Output(&[1.0]),
            Direction::default(),
            None,
            (7.0, 13.0),
            cap,
            || false,
        ),
        Err(DrivenResonanceError::Invalid(_))
    ));
    let mut polls = 0;
    assert_eq!(
        model.driven_resonance_with_parameter_derivative(
            &[C64::ONE],
            HarmonicObservation::Output(&[1.0]),
            Direction::default(),
            None,
            (7.0, 13.0),
            controls(true),
            || {
                polls += 1;
                polls > 120
            },
        ),
        Err(DrivenResonanceError::Harmonic(HarmonicError::Cancelled))
    );
    let compute = || {
        model
            .driven_resonance_with_parameter_derivative(
                &[C64::ONE],
                HarmonicObservation::Output(&[1.0]),
                Direction::default(),
                None,
                (7.0, 13.0),
                controls(true),
                || false,
            )
            .unwrap()
    };
    let first = compute();
    assert_eq!(first, compute());
    assert_eq!(first.peak_frequency_derivative_rad_s, 0.0);
    assert_eq!(first.observation_derivative, C64::ZERO);
}
