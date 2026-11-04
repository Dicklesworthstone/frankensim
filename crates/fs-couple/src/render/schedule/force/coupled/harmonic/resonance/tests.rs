use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};

const OMEGA: f64 = 1000.0;
const ZETA: f64 = 0.08;

fn oscillator(modes: usize) -> CoupledModalSystem {
    let gate = CancelGate::new_clock_free();
    let model = ModalAcousticTimeModel::try_new(
        48_000,
        vec![
            ModalAcousticMode {
                angular_frequency_rad_s: OMEGA,
                damping_ratio: ZETA,
                pressure_per_modal_velocity: C64::ZERO,
            };
            modes
        ],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    CoupledModalSystem::new_with_columns(
        vec![model],
        vec![],
        ModalCouplingConfig {
            max_modes: 8,
            max_connections: 8,
            max_setup_terms: 1_000_000,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 100.0,
            maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 1e6,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-9,
        },
        &gate,
    )
    .unwrap()
}

fn controls() -> ResonanceControls {
    ResonanceControls {
        absolute_frequency_tolerance_rad_s: 1e-6,
        relative_frequency_tolerance: 1e-12,
        maximum_evaluations: 256,
        measure_half_power_bandwidth: true,
    }
}

fn frequency(actual: ResonanceFrequency, expected: f64) {
    assert!(actual.bracket_rad_s.0 <= expected + 1e-9);
    assert!(actual.bracket_rad_s.1 >= expected - 1e-9);
    assert!(actual.bracket_rad_s.1 - actual.bracket_rad_s.0 <= 1.001e-6);
    assert!((actual.angular_frequency_rad_s - expected).abs() <= 1e-6);
}

#[test]
fn displacement_and_velocity_peaks_and_widths_match_original_oscillator_equations() {
    let gate = CancelGate::new_clock_free();
    let system = oscillator(1);
    let force = C64::new(0.7, -0.2);
    let drag = 2.0 * ZETA * OMEGA;
    let displacement_peak_sq = OMEGA * OMEGA - 0.5 * drag * drag;
    let minimum_denominator = drag * drag * (OMEGA * OMEGA - 0.25 * drag * drag);
    let displacement_edges = (
        det::sqrt(displacement_peak_sq - det::sqrt(minimum_denominator)),
        det::sqrt(displacement_peak_sq + det::sqrt(minimum_denominator)),
    );
    let velocity_edges = (
        0.5 * (det::sqrt(drag * drag + 4.0 * OMEGA * OMEGA) - drag),
        0.5 * (det::sqrt(drag * drag + 4.0 * OMEGA * OMEGA) + drag),
    );
    let weight = [2.3];
    let mut results = Vec::new();
    for (observation, peak, edges, velocity) in [
        (
            ModalHarmonicObservation::Displacement(&weight),
            det::sqrt(displacement_peak_sq),
            displacement_edges,
            false,
        ),
        (
            ModalHarmonicObservation::Velocity(&weight),
            OMEGA,
            velocity_edges,
            true,
        ),
    ] {
        let result = system
            .driven_resonance(
                &[force],
                observation,
                (600.0, 1400.0),
                controls(),
                100_000,
                &gate,
            )
            .unwrap();
        frequency(result.peak, peak);
        let band = result.half_power.unwrap();
        frequency(band.lower, edges.0);
        frequency(band.upper, edges.1);
        assert!((band.bandwidth_rad_s - (edges.1 - edges.0)).abs() < 2e-6);
        assert!((band.response_quality_factor - peak / (edges.1 - edges.0)).abs() < 2e-7);
        let omega = result.peak.angular_frequency_rad_s;
        let displacement =
            force.scale(weight[0]) / C64::new((OMEGA - omega) * (OMEGA + omega), -drag * omega);
        let expected = if velocity {
            C64::new(0.0, -omega) * displacement
        } else {
            displacement
        };
        assert!((result.observation - expected).abs() < 1e-12 * expected.abs());
        assert!(
            (result.observation_derivative / result.observation)
                .re
                .abs()
                < 1e-9
        );
        assert!(result.maximum_relative_residual < 1e-11);
        assert!(result.evaluations <= controls().maximum_evaluations);
        results.push(result);
    }
    assert!(
        results[1].peak.angular_frequency_rad_s - results[0].peak.angular_frequency_rad_s > 6.0
    );
    assert_eq!(system.samples_rendered(), 0);
}

#[test]
fn acceleration_peak_and_bandwidth_follow_original_inertial_motion() {
    let gate = CancelGate::new_clock_free();
    let system = oscillator(1);
    let force = C64::new(0.7, -0.2);
    let weight = [2.3];
    let drag = 2.0 * ZETA * OMEGA;
    let displacement_peak_sq = OMEGA * OMEGA - 0.5 * drag * drag;
    let minimum_denominator = drag * drag * (OMEGA * OMEGA - 0.25 * drag * drag);
    // Accelerance at omega has the same normalized magnitude as receptance
    // at OMEGA²/omega. This gives independent closed-form peak and edge roots.
    let peak = OMEGA * OMEGA / det::sqrt(displacement_peak_sq);
    let edges = (
        OMEGA * OMEGA / det::sqrt(displacement_peak_sq + det::sqrt(minimum_denominator)),
        OMEGA * OMEGA / det::sqrt(displacement_peak_sq - det::sqrt(minimum_denominator)),
    );
    let result = system
        .driven_resonance(
            &[force],
            ModalHarmonicObservation::Acceleration(&weight),
            (600.0, 1400.0),
            controls(),
            100_000,
            &gate,
        )
        .unwrap();
    frequency(result.peak, peak);
    let band = result.half_power.unwrap();
    frequency(band.lower, edges.0);
    frequency(band.upper, edges.1);
    assert!((band.bandwidth_rad_s - (edges.1 - edges.0)).abs() < 2e-6);
    let omega = result.peak.angular_frequency_rad_s;
    let denominator = C64::new((OMEGA - omega) * (OMEGA + omega), -drag * omega);
    let q = force.scale(weight[0]) / denominator;
    let q_slope = -(C64::new(-2.0 * omega, -drag) * q) / denominator;
    let acceleration = q.scale(-omega * omega);
    let acceleration_slope = q_slope.scale(-omega * omega) - q.scale(2.0 * omega);
    assert!((result.observation - acceleration).abs() < 1e-12 * acceleration.abs());
    assert!((result.observation_derivative - acceleration_slope).abs() < 1e-12);
    assert!(
        (result.observation_derivative / result.observation)
            .re
            .abs()
            < 1e-9
    );
    assert!(
        peak - OMEGA > 6.0,
        "acceleration must not reuse the velocity peak"
    );
    assert_eq!(system.samples_rendered(), 0);
}

#[test]
fn quiet_signed_modal_interference_retains_the_observed_resonance() {
    let gate = CancelGate::new_clock_free();
    let system = oscillator(3);
    // Equal and opposite original modal responses surround a much quieter
    // coherent contribution. Ordinary summation loses that contribution.
    let result = system
        .driven_resonance(
            &[C64::ONE, C64::from_re(1e-16), C64::from_re(-1.0)],
            ModalHarmonicObservation::Velocity(&[1e16, 1e16, 1e16]),
            (600.0, 1400.0),
            controls(),
            100_000,
            &gate,
        )
        .unwrap();
    frequency(result.peak, OMEGA);
    assert!((result.observation.abs() - 1.0 / (2.0 * ZETA * OMEGA)).abs() < 1e-13);
    assert!((result.half_power.unwrap().bandwidth_rad_s - 2.0 * ZETA * OMEGA).abs() < 2e-6);
    assert_eq!(system.samples_rendered(), 0);
}

#[test]
fn modal_peak_refusals_preserve_source_errors_and_state() {
    let gate = CancelGate::new_clock_free();
    let system = oscillator(1);
    let observation = ModalHarmonicObservation::Velocity(&[1.0]);
    assert!(matches!(
        system.driven_resonance(
            &[C64::ZERO],
            observation,
            (600.0, 1400.0),
            controls(),
            100_000,
            &gate,
        ),
        Err(ResonanceSearchError::NotBracketed(_))
    ));
    assert!(matches!(
        system.driven_resonance(
            &[C64::ONE],
            observation,
            (200.0, 400.0),
            controls(),
            100_000,
            &gate,
        ),
        Err(ResonanceSearchError::NotBracketed(_))
    ));
    assert!(matches!(
        system.driven_resonance(
            &[C64::ONE],
            observation,
            (950.0, 1050.0),
            controls(),
            100_000,
            &gate,
        ),
        Err(ResonanceSearchError::NotBracketed(_))
    ));
    for weights in [&[][..], &[0.0][..], &[f64::NAN][..], &[1.0, 1.0][..]] {
        assert!(matches!(
            system.driven_resonance(
                &[C64::ONE],
                ModalHarmonicObservation::Displacement(weights),
                (600.0, 1400.0),
                controls(),
                100_000,
                &gate,
            ),
            Err(ResonanceSearchError::Source(
                ModalCouplingError::Invalid { .. }
            ))
        ));
    }
    assert!(matches!(
        system.driven_resonance(
            &[C64::ONE],
            observation,
            (600.0, 1400.0),
            controls(),
            1,
            &gate,
        ),
        Err(ResonanceSearchError::Source(
            ModalCouplingError::Invalid { .. }
        ))
    ));
    let mut limited = controls();
    limited.maximum_evaluations = 3;
    assert!(matches!(
        system.driven_resonance(
            &[C64::ONE],
            observation,
            (600.0, 1400.0),
            limited,
            100_000,
            &gate,
        ),
        Err(ResonanceSearchError::EvaluationLimit { maximum: 3 })
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        system.driven_resonance(
            &[C64::ONE],
            observation,
            (600.0, 1400.0),
            controls(),
            100_000,
            &cancelled,
        ),
        Err(ResonanceSearchError::Source(ModalCouplingError::Cancelled))
            | Err(ResonanceSearchError::Cancelled)
    ));
    assert_eq!(system.samples_rendered(), 0);
    assert!(system.material_memory().is_empty());
    assert!(system.radiation_state().is_empty());
}
