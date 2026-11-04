//! G1/G3: actual driven peaks and linewidths, rather than free-decay poles.
use fs_math::c64::C64;
use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
use fs_phs::harmonic::resonance::{
    DrivenResonanceError, HarmonicObservation, ResonanceControls, ResonanceFrequency,
};
use fs_phs::harmonic::{HarmonicBudget, HarmonicError, HarmonicModel};

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 32,
        max_dense_bytes: 2_000_000,
        max_dense_work: 4_000_000,
        relative_residual_tolerance: 1e-12,
    }
}

fn controls() -> ResonanceControls {
    ResonanceControls {
        absolute_frequency_tolerance_rad_s: 1e-10,
        relative_frequency_tolerance: 0.0,
        maximum_evaluations: 256,
        measure_half_power_bandwidth: true,
    }
}

fn frequency(actual: ResonanceFrequency, expected: f64, tolerance: f64) {
    assert!(
        (actual.angular_frequency_rad_s - expected).abs() <= tolerance,
        "{actual:?} does not resolve {expected}"
    );
    assert!(
        actual.bracket_rad_s.1 - actual.bracket_rad_s.0 <= tolerance,
        "frequency bracket was not refined: {actual:?}"
    );
    assert!(
        actual.angular_frequency_rad_s >= actual.bracket_rad_s.0
            && actual.angular_frequency_rad_s <= actual.bracket_rad_s.1
    );
}

#[test]
fn displacement_peak_and_half_height_differ_from_natural_and_free_decay_frequencies() {
    let (m, k, c) = (2.0_f64, 200.0_f64, 8.0_f64);
    let system = fs_phs::mass_spring_damper(m, k, c).unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    let result = model
        .driven_resonance(
            &[C64::new(0.8, -0.3)],
            HarmonicObservation::State(&[1.0, 0.0]),
            (4.0, 15.0),
            controls(),
            || false,
        )
        .unwrap();
    // |q/F|^-2 = m^2(omega^2-omega_peak^2)^2 + D_min.
    let peak_squared = k / m - c * c / (2.0 * m * m);
    let minimum = c * c * k / m - c.powi(4) / (4.0 * m * m);
    frequency(result.peak, peak_squared.sqrt(), 1e-10);
    assert!((result.peak.angular_frequency_rad_s - (k / m).sqrt()).abs() > 0.3);
    assert!(
        (result.peak.angular_frequency_rad_s - (k / m - (c / (2.0 * m)).powi(2)).sqrt()).abs()
            > 0.1
    );
    let band = result.half_power.unwrap();
    frequency(
        band.lower,
        (peak_squared - minimum.sqrt() / m).sqrt(),
        1e-10,
    );
    frequency(
        band.upper,
        (peak_squared + minimum.sqrt() / m).sqrt(),
        1e-10,
    );
    let omega = result.peak.angular_frequency_rad_s;
    let exact = C64::new(0.8, -0.3) / C64::new(k - m * omega * omega, -c * omega);
    assert!((result.observation - exact).abs() < 1e-14);
    assert!(result.maximum_relative_residual < 1e-12);
    assert!(result.evaluations <= controls().maximum_evaluations);
}

#[test]
fn very_narrow_velocity_peak_resolves_physical_linewidth_without_a_frequency_grid() {
    let (m, k, c) = (2.0_f64, 200.0_f64, 1e-4_f64);
    let system = fs_phs::mass_spring_damper(m, k, c).unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    let result = model
        .driven_resonance(
            &[C64::ONE],
            HarmonicObservation::Output(&[1.0]),
            (8.0, 12.0),
            controls(),
            || false,
        )
        .unwrap();
    frequency(result.peak, 10.0, 1e-10);
    let band = result.half_power.unwrap();
    let root = (c * c + 4.0 * m * k).sqrt();
    frequency(band.lower, (root - c) / (2.0 * m), 1e-10);
    frequency(band.upper, (root + c) / (2.0 * m), 1e-10);
    assert!((band.bandwidth_rad_s - c / m).abs() < 1e-10);
    assert!((band.response_quality_factor / (10.0 * m / c) - 1.0).abs() < 3e-6);
    assert!((result.observation.re * c - 1.0).abs() < 1e-10);
    assert!(result.evaluations < 180);
}

#[test]
fn constrained_rlc_current_uses_the_same_peak_and_bandwidth_owner() {
    let (resistance, inductance, capacitance) = (0.02_f64, 0.25_f64, 0.04_f64);
    let graph = CircuitGraph {
        node_count: 4,
        transformers: vec![],
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (
                1,
                2,
                Branch::Inductor {
                    henries: inductance,
                },
            ),
            (
                2,
                3,
                Branch::Capacitor {
                    farads: capacitance,
                },
            ),
            (3, 0, Branch::Resistor { ohms: resistance }),
        ],
    };
    let circuit = assemble_circuit(&graph).unwrap();
    let model = HarmonicModel::from_descriptor(&circuit.system, budget(), || false).unwrap();
    let result = model
        .driven_resonance(
            &[C64::new(2.0, -0.4)],
            HarmonicObservation::Output(&[1.0]),
            (7.0, 14.0),
            controls(),
            || false,
        )
        .unwrap();
    frequency(result.peak, 1.0 / (inductance * capacitance).sqrt(), 1e-10);
    let band = result.half_power.unwrap();
    assert!((band.bandwidth_rad_s - resistance / inductance).abs() < 1e-10);
    assert!((result.observation - C64::new(2.0, -0.4).scale(1.0 / resistance)).abs() < 2e-7);
}

#[test]
fn excitation_and_observation_scaling_preserve_frequency_and_bandwidth() {
    let system = fs_phs::mass_spring_damper(1.0, 100.0, 0.5).unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    let evaluate = |drive, weight| {
        model
            .driven_resonance(
                &[drive],
                HarmonicObservation::Effort(&[weight, 0.0]),
                (6.0, 14.0),
                controls(),
                || false,
            )
            .unwrap()
    };
    let reference = evaluate(C64::ONE, 1.0);
    assert_eq!(reference, evaluate(C64::ONE, 1.0));
    for (drive, weight) in [
        (C64::new(1e-60, -2e-60), -3.0),
        (C64::new(-1e60, 2e60), 0.25),
    ] {
        let result = evaluate(drive, weight);
        assert!(
            (result.peak.angular_frequency_rad_s - reference.peak.angular_frequency_rad_s).abs()
                < 1e-10
        );
        assert!(
            (result.half_power.unwrap().bandwidth_rad_s
                - reference.half_power.unwrap().bandwidth_rad_s)
                .abs()
                < 2e-10
        );
    }
}

#[test]
fn missing_brackets_silence_and_unresolved_bandwidth_do_not_manufacture_resonances() {
    let system = fs_phs::mass_spring_damper(1.0, 100.0, 1.0).unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    for range in [(1.0, 2.0), (20.0, 30.0), (9.9, 10.1)] {
        assert!(matches!(
            model.driven_resonance(
                &[C64::ONE],
                HarmonicObservation::Output(&[1.0]),
                range,
                controls(),
                || false
            ),
            Err(DrivenResonanceError::NotBracketed(_))
        ));
    }
    let peak_only = ResonanceControls {
        measure_half_power_bandwidth: false,
        ..controls()
    };
    assert!(
        model
            .driven_resonance(
                &[C64::ONE],
                HarmonicObservation::Output(&[1.0]),
                (9.9, 10.1),
                peak_only,
                || false
            )
            .unwrap()
            .half_power
            .is_none()
    );
    assert!(matches!(
        model.driven_resonance(
            &[C64::ZERO],
            HarmonicObservation::Output(&[1.0]),
            (8.0, 12.0),
            controls(),
            || false
        ),
        Err(DrivenResonanceError::NotBracketed(_))
    ));
    assert!(matches!(
        model.driven_resonance(
            &[C64::ONE],
            HarmonicObservation::Output(&[0.0]),
            (8.0, 12.0),
            controls(),
            || false
        ),
        Err(DrivenResonanceError::Invalid(_))
    ));
}

#[test]
fn cancellation_and_evaluation_limit_leave_the_prepared_model_reusable() {
    let system = fs_phs::mass_spring_damper(1.0, 100.0, 0.1).unwrap();
    let model = HarmonicModel::from_ode(&system, budget(), || false).unwrap();
    let limited = ResonanceControls {
        maximum_evaluations: 3,
        ..controls()
    };
    assert!(matches!(
        model.driven_resonance(
            &[C64::ONE],
            HarmonicObservation::Output(&[1.0]),
            (8.0, 12.0),
            limited,
            || false
        ),
        Err(DrivenResonanceError::EvaluationLimit { maximum: 3 })
    ));
    let mut polls = 0;
    assert!(matches!(
        model.driven_resonance(
            &[C64::ONE],
            HarmonicObservation::Output(&[1.0]),
            (8.0, 12.0),
            controls(),
            || {
                polls += 1;
                polls > 100
            }
        ),
        Err(DrivenResonanceError::Harmonic(HarmonicError::Cancelled))
    ));
    let retried = model
        .driven_resonance(
            &[C64::ONE],
            HarmonicObservation::Output(&[1.0]),
            (8.0, 12.0),
            controls(),
            || false,
        )
        .unwrap();
    frequency(retried.peak, 10.0, 1e-10);
}
