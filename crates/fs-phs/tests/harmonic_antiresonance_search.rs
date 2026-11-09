//! Observed notches retain coherent zeros, native load energy, and one search owner.
use core::convert::Infallible;
use std::cell::Cell;

use fs_math::c64::C64;
use fs_phs::{
    PortHamiltonian, QuadraticStorage,
    harmonic::{
        HarmonicBudget, HarmonicModel,
        resonance::{
            AntiresonanceControls, DrivenResonanceError, HarmonicObservation, ResonanceSearchError,
            search_antiresonance,
        },
    },
};

fn controls() -> AntiresonanceControls {
    AntiresonanceControls {
        absolute_frequency_tolerance_rad_s: 1e-10,
        relative_frequency_tolerance: 0.0,
        maximum_evaluations: 256,
        recovery_power: None,
    }
}

#[test]
fn exact_coherent_zero_keeps_its_tangent_and_absolute_recovery_edges_at_extreme_units() {
    let center = 13.25;
    for scale in [1e-120, 1.0, 1e120] {
        let phase = C64::new(0.6, -0.8).scale(scale);
        let queries = Cell::new(0_usize);
        let center_queries = Cell::new(0_usize);
        let tangent_queries = Cell::new(0_usize);
        let level = 0.25 * scale * scale;
        let result = search_antiresonance(
            (center - 2.0, center + 2.0),
            AntiresonanceControls {
                recovery_power: Some(level),
                ..controls()
            },
            || false,
            |omega, need_derivative, cancelled| {
                assert!(!cancelled());
                queries.set(queries.get() + 1);
                center_queries.set(center_queries.get() + usize::from(omega == center));
                tangent_queries.set(tangent_queries.get() + usize::from(need_derivative));
                Ok::<_, Infallible>((
                    phase.scale(omega - center),
                    if need_derivative {
                        phase
                    } else {
                        // Recovery queries must not require or inspect tangents.
                        C64::new(f64::NAN, f64::NAN)
                    },
                    queries.get() as f64 * 1e-17,
                ))
            },
        )
        .unwrap();
        assert_eq!(result.minimum.angular_frequency_rad_s, center);
        assert_eq!(result.observation, C64::ZERO);
        assert_eq!(result.observation_derivative, phase);
        assert!(result.minimum.bracket_rad_s.0 < center);
        assert!(result.minimum.bracket_rad_s.1 > center);
        assert!(result.minimum.bracket_rad_s.1 - result.minimum.bracket_rad_s.0 <= 1e-10);
        assert_eq!(
            center_queries.get(),
            1,
            "the stationary solve must be reused"
        );
        let band = result.recovery.unwrap();
        assert_eq!(band.power_level, level);
        assert!((band.lower.angular_frequency_rad_s - (center - 0.5)).abs() <= 1e-10);
        assert!((band.upper.angular_frequency_rad_s - (center + 0.5)).abs() <= 1e-10);
        assert!((band.bandwidth_rad_s - 1.0).abs() <= 2e-10);
        assert_eq!(result.evaluations, queries.get());
        assert!(tangent_queries.get() > 2 && tangent_queries.get() < queries.get());
        assert_eq!(
            result.maximum_relative_residual,
            queries.get() as f64 * 1e-17
        );
    }
}

#[test]
fn damped_nonzero_notch_resolves_off_grid_minimum_and_analytic_recovery_width() {
    // Complete coherent direct-plus-resonant response. Finite internal loss
    // leaves a/b at the notch; its measured width is not a pole's linewidth.
    let (center, a, b, level) = (13.234567891_f64, 0.03_f64, 0.4_f64, 0.25_f64);
    let phase = C64::new(-0.8, -0.6);
    let result = search_antiresonance(
        (10.0, 16.0),
        AntiresonanceControls {
            recovery_power: Some(level),
            ..controls()
        },
        || false,
        |omega, need_derivative, _| {
            let delta = omega - center;
            let denominator = C64::new(delta, b);
            let value = phase * (C64::new(delta, a) / denominator);
            let derivative = if need_derivative {
                phase * (C64::new(0.0, b - a) / (denominator * denominator))
            } else {
                C64::ZERO
            };
            Ok::<_, Infallible>((value, derivative, 3e-16))
        },
    )
    .unwrap();
    assert!((result.minimum.angular_frequency_rad_s - center).abs() <= 1e-10);
    assert!(result.minimum.bracket_rad_s.0 <= center);
    assert!(result.minimum.bracket_rad_s.1 >= center);
    assert!((result.observation - phase.scale(a / b)).abs() < 3e-10);
    assert!(result.observation.abs() > 0.07);
    let band = result.recovery.unwrap();
    let offset = ((level * b * b - a * a) / (1.0 - level)).sqrt();
    assert!((band.lower.angular_frequency_rad_s - center + offset).abs() <= 1e-10);
    assert!((band.upper.angular_frequency_rad_s - center - offset).abs() <= 1e-10);
    assert!((band.bandwidth_rad_s - 2.0 * offset).abs() <= 2e-10);
    assert_eq!(result.maximum_relative_residual, 3e-16);
}

#[test]
fn source_amplitude_and_phase_remain_in_the_total_notch_derivative() {
    let (center, floor, source_slope) = (20.0_f64, 0.2_f64, 0.7_f64);
    let phase = C64::new(0.6, -0.8);
    let result = search_antiresonance(
        (19.8, 20.2),
        controls(),
        || false,
        |omega, _, _| {
            let delta = omega - center;
            let source = 1.0 + source_slope * delta;
            let response = C64::new(delta, floor);
            Ok::<_, Infallible>((
                phase * response.scale(source),
                phase * (C64::from_re(source) + response.scale(source_slope)),
                0.0,
            ))
        },
    )
    .unwrap();
    let discriminant = (1.0 - 8.0 * source_slope.powi(2) * floor.powi(2)).sqrt();
    let offset = -2.0 * source_slope * floor * floor / (1.0 + discriminant);
    assert!((result.minimum.angular_frequency_rad_s - center - offset).abs() <= 1e-10);
    assert!(result.minimum.angular_frequency_rad_s < center - 0.02);
    let expected = phase * C64::new(offset, floor).scale(1.0 + source_slope * offset);
    assert!((result.observation - expected).abs() <= 1e-10);
    assert!(result.recovery.is_none());
}

#[test]
fn silent_monotone_maximum_and_flat_observations_do_not_become_notches() {
    for case in 0..4 {
        let result = search_antiresonance(
            (9.0, 11.0),
            controls(),
            || false,
            |omega, _, _| {
                let delta = omega - 10.0;
                let (value, derivative) = match case {
                    0 => (0.0, 0.0),
                    1 => (omega, 1.0),
                    2 => (2.0 - delta * delta, -2.0 * delta),
                    _ => {
                        let outside = (delta.abs() - 0.25).max(0.0);
                        (outside * outside, 2.0 * delta.signum() * outside)
                    }
                };
                Ok::<_, Infallible>((C64::from_re(value), C64::from_re(derivative), 0.0))
            },
        );
        assert!(
            matches!(result, Err(ResonanceSearchError::NotBracketed(_))),
            "case {case}: {result:?}"
        );
    }
    for recovery_power in [0.01, 9.0] {
        let result = search_antiresonance(
            (9.0, 11.0),
            AntiresonanceControls {
                recovery_power: Some(recovery_power),
                ..controls()
            },
            || false,
            |omega, _, _| Ok::<_, Infallible>((C64::new(omega - 10.0, 0.2), C64::ONE, 0.0)),
        );
        assert!(matches!(result, Err(ResonanceSearchError::NotBracketed(_))));
    }
}

#[derive(Debug)]
struct OriginalRefusal {
    message: String,
}

#[test]
fn recovery_source_refusal_and_one_total_evaluation_cap_are_preserved() {
    let failure = search_antiresonance(
        (9.0, 11.0),
        AntiresonanceControls {
            recovery_power: Some(0.25),
            ..controls()
        },
        || false,
        |omega, derivative, _| {
            if !derivative {
                return Err(OriginalRefusal {
                    message: "original coupled load refused its recovery solve".to_owned(),
                });
            }
            Ok((C64::from_re(omega - 10.0), C64::ONE, 0.0))
        },
    );
    match failure {
        Err(ResonanceSearchError::Source(error)) => assert!(error.message.contains("coupled load")),
        other => panic!("native source refusal was lost: {other:?}"),
    }
    let queries = Cell::new(0);
    let capped = search_antiresonance(
        (9.0, 11.0),
        AntiresonanceControls {
            maximum_evaluations: 3,
            ..controls()
        },
        || false,
        |omega, _, _| {
            queries.set(queries.get() + 1);
            Ok::<_, Infallible>((C64::from_re(omega - 10.0), C64::ONE, 0.0))
        },
    );
    assert_eq!(
        capped,
        Err(ResonanceSearchError::EvaluationLimit { maximum: 3 })
    );
    assert_eq!(queries.get(), 3);
}

#[test]
fn cancellation_before_after_and_at_publication_never_returns_a_partial_notch() {
    let checks = Cell::new(0);
    let trial = |stop_at: usize| {
        checks.set(0);
        search_antiresonance(
            (9.0, 11.0),
            controls(),
            || {
                checks.set(checks.get() + 1);
                checks.get() == stop_at
            },
            |omega, _, _| Ok::<_, Infallible>((C64::from_re(omega - 10.0), C64::ONE, 0.0)),
        )
    };
    trial(usize::MAX).unwrap();
    let complete_checks = checks.get();
    for stop_at in [1, 3, complete_checks] {
        assert_eq!(trial(stop_at), Err(ResonanceSearchError::Cancelled));
    }
    let cancel = Cell::new(false);
    let queries = Cell::new(0);
    let after_source = search_antiresonance(
        (9.0, 11.0),
        controls(),
        || cancel.get(),
        |omega, _, _| {
            queries.set(queries.get() + 1);
            cancel.set(queries.get() == 3);
            Ok::<_, Infallible>((C64::from_re(omega - 10.0), C64::ONE, 0.0))
        },
    );
    assert_eq!(queries.get(), 3);
    assert_eq!(after_source, Err(ResonanceSearchError::Cancelled));
}

#[test]
fn invalid_controls_nonfinite_phasors_and_invalid_residuals_refuse() {
    for level in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        let result = search_antiresonance(
            (9.0, 11.0),
            AntiresonanceControls {
                recovery_power: Some(level),
                ..controls()
            },
            || false,
            |_, _, _| -> Result<(C64, C64, f64), Infallible> {
                panic!("invalid controls must refuse before physical work")
            },
        );
        assert!(matches!(result, Err(ResonanceSearchError::Invalid(_))));
    }
    for sample in [
        (C64::new(f64::NAN, 0.0), C64::ONE, 0.0),
        (C64::ONE, C64::new(0.0, f64::INFINITY), 0.0),
        (C64::ONE, C64::ONE, -1.0),
        (C64::ONE, C64::ONE, f64::NAN),
    ] {
        let result = search_antiresonance(
            (9.0, 11.0),
            controls(),
            || false,
            |_, _, _| Ok::<_, Infallible>(sample),
        );
        assert!(matches!(result, Err(ResonanceSearchError::Invalid(_))));
    }
}

// Original host/absorber displacements and momenta. At omega=sqrt(ka/ma)=10
// the host is silent while the lossless absorber carries finite energy and an
// opposing force. Host damping remains present in the original dynamic pencil.
fn absorber() -> HarmonicModel {
    let (host_mass, absorber_mass, host_stiffness, absorber_stiffness) = (2.0, 0.5, 200.0, 50.0);
    let q = vec![
        host_stiffness + absorber_stiffness,
        0.0,
        -absorber_stiffness,
        0.0,
        0.0,
        1.0 / host_mass,
        0.0,
        0.0,
        -absorber_stiffness,
        0.0,
        absorber_stiffness,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0 / absorber_mass,
    ];
    let j = vec![
        0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0,
    ];
    let mut r = vec![0.0; 16];
    r[5] = 0.4;
    let system = PortHamiltonian::new(
        4,
        1,
        j,
        r,
        vec![0.0, 1.0, 0.0, 0.0],
        Box::new(QuadraticStorage::new(q, 4).unwrap()),
    )
    .unwrap();
    HarmonicModel::from_ode(
        &system,
        HarmonicBudget {
            max_dimension: 4,
            max_dense_bytes: 16_384,
            max_dense_work: 10_000,
            relative_residual_tolerance: 1e-10,
        },
        || false,
    )
    .unwrap()
}

#[test]
fn native_complete_model_finds_cancellation_without_erasing_absorber_energy() {
    let model = absorber();
    let input = [C64::new(0.6, -0.8)];
    for observation in [
        HarmonicObservation::State(&[1.0, 0.0, 0.0, 0.0]),
        HarmonicObservation::Output(&[1.0]),
    ] {
        let result = model
            .driven_antiresonance(&input, observation, (9.0, 11.0), controls(), || false)
            .unwrap();
        assert!((result.minimum.angular_frequency_rad_s - 10.0).abs() <= 1e-10);
        assert!(result.observation.abs() < 1e-11);
        assert!(result.observation_derivative.abs() > 0.0039);
        assert!(result.maximum_relative_residual < 1e-10);
        let actual = model
            .factor(result.minimum.angular_frequency_rad_s, || false)
            .unwrap()
            .response(&input, || false)
            .unwrap();
        let projected = match observation {
            HarmonicObservation::State(_) => actual.state[0],
            _ => actual.output[0],
        };
        assert_eq!(result.observation, projected);
        assert!((actual.state[2] + input[0].scale(1.0 / 50.0)).abs() < 1e-11);
        assert!((actual.mean_stored_energy_j - 0.01).abs() < 1e-10);
    }
    assert!(matches!(
        model.driven_antiresonance(
            &input,
            HarmonicObservation::State(&[0.0; 4]),
            (9.0, 11.0),
            controls(),
            || false,
        ),
        Err(DrivenResonanceError::Invalid(_))
    ));
    assert!(matches!(
        model.driven_antiresonance(
            &[C64::ZERO],
            HarmonicObservation::Output(&[1.0]),
            (9.0, 11.0),
            controls(),
            || false,
        ),
        Err(DrivenResonanceError::NotBracketed(_))
    ));
}
