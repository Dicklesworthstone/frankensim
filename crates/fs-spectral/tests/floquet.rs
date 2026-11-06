//! Closed-form monodromies keep physical stability separate from phase freedom.

use fs_math::c64::C64;
use fs_spectral::{
    dense::{DenseSpectrumError, DenseSpectrumOptions},
    floquet::{
        FloquetClassification, FloquetControls, FloquetError, analyze_monodromy,
        assess_monodromy_with_multipliers,
    },
};

fn controls() -> FloquetControls {
    FloquetControls {
        spectrum: DenseSpectrumOptions {
            maximum_dimension: 16,
            maximum_dense_bytes: 1 << 20,
            maximum_work: 1_000_000_000,
            qr_max_sweeps_per_window: 160,
            svd_max_sweeps: 100,
            svd_tolerance: 2e-14,
            maximum_relative_residual: 1e-11,
            minimum_eigenbasis_reciprocal_condition: 1e-8,
        },
        unit_circle_guard: 1e-4,
        phase_multiplier_tolerance: 1e-8,
        maximum_phase_residual: 1e-10,
    }
}

fn real(values: &[f64]) -> Vec<C64> {
    values.iter().copied().map(C64::from_re).collect()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * (1.0 + expected.abs()),
        "{actual:.16e} != {expected:.16e}"
    );
}

fn near_complex(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * (1.0 + expected.abs()),
        "{actual:?} != {expected:?}"
    );
}

// S * (rho R) * S^-1, S=diag(4,1), cos(theta)=3/5, sin(theta)=4/5.
// This is a coupled, nonnormal real period map with known complex multipliers.
fn rotating_map(radius: f64) -> Vec<C64> {
    real(&[0.6 * radius, -3.2 * radius, 0.2 * radius, 0.6 * radius])
}

fn rotating_multipliers(radius: f64) -> Vec<C64> {
    vec![
        C64::new(0.6 * radius, -0.8 * radius),
        C64::new(0.6 * radius, 0.8 * radius),
    ]
}

#[test]
fn forced_nonnormal_modes_preserve_complex_multipliers_and_growth_sign() {
    let period = 0.4;
    for (radius, expected_class) in [
        (0.7, FloquetClassification::Stable),
        (1.25, FloquetClassification::Unstable),
    ] {
        let matrix = rotating_map(radius);
        let expected = rotating_multipliers(radius);
        let discovered = analyze_monodromy(&matrix, 2, period, None, controls(), || false).unwrap();
        let supplied = assess_monodromy_with_multipliers(
            &matrix,
            2,
            period,
            &expected,
            None,
            controls(),
            || false,
        )
        .unwrap();
        for result in [&discovered, &supplied] {
            assert_eq!(result.classification, expected_class);
            assert_eq!(result.spectrum.values.len(), 2);
            assert_eq!(result.spectrum.right_vectors.len(), 2);
            for (&actual, &known) in result.spectrum.values.iter().zip(&expected) {
                near_complex(actual, known);
            }
            assert!(
                result
                    .spectrum
                    .relative_residuals
                    .iter()
                    .all(|&r| r < 1e-11)
            );
            assert!(result.spectrum.eigenbasis_reciprocal_condition > 0.1);
            near(result.spectral_radius, radius);
            near(result.effective_spectral_radius.unwrap(), radius);
            near(
                result.dominant_growth_rate_per_s.unwrap(),
                radius.ln() / period,
            );
            assert_eq!(result.phase_multiplier_index, None);
            assert_eq!(result.phase_residual, None);
        }
    }
}

#[test]
fn autonomous_full_flow_identifies_one_phase_mode_and_keeps_transverse_memory() {
    let period = 1.25;
    // The full original phase direction is (1,2,3), not a coordinate vector.
    // A common complex phase cannot change its physical eigendirection.
    let phase_factor = C64::new(0.6, 0.8);
    let flow: Vec<_> = [1.0, 2.0, 3.0]
        .into_iter()
        .map(|x| phase_factor * C64::from_re(x))
        .collect();
    for (transverse, expected_class, phase_index) in [
        (0.5, FloquetClassification::Stable, 2),
        (1.2, FloquetClassification::Unstable, 1),
    ] {
        // Eigenvectors e0, e1, (1,2,3); the first memory multiplier is 1/4.
        let matrix = real(&[
            0.25,
            0.0,
            0.25,
            0.0,
            transverse,
            2.0 * (1.0 - transverse) / 3.0,
            0.0,
            0.0,
            1.0,
        ]);
        let result =
            analyze_monodromy(&matrix, 3, period, Some(&flow), controls(), || false).unwrap();
        assert_eq!(result.classification, expected_class);
        assert_eq!(result.spectrum.values.len(), 3);
        assert_eq!(result.spectrum.right_vectors.len(), 3);
        assert_eq!(result.phase_multiplier_index, Some(phase_index));
        near_complex(result.spectrum.values[phase_index], C64::ONE);
        near_complex(result.spectrum.values[0], C64::from_re(0.25));
        assert!(result.phase_residual.unwrap() <= 1e-12);
        near(result.spectral_radius, transverse.max(1.0));
        near(result.effective_spectral_radius.unwrap(), transverse);
        near(
            result.dominant_growth_rate_per_s.unwrap(),
            transverse.ln() / period,
        );
    }
}

#[test]
fn guarded_unit_modes_and_a_false_phase_direction_remain_unresolved() {
    let matrix = real(&[0.5, 0.0, 0.0, 1.0]);
    let forced = analyze_monodromy(&matrix, 2, 1.0, None, controls(), || false).unwrap();
    assert_eq!(forced.classification, FloquetClassification::Unresolved);
    assert_eq!(forced.phase_multiplier_index, None);
    near(forced.effective_spectral_radius.unwrap(), 1.0);
    near(forced.dominant_growth_rate_per_s.unwrap(), 0.0);

    // Merely finding an eigenvalue +1 does not validate the supplied flow.
    let false_flow = real(&[1.0, 1.0]);
    let false_phase =
        analyze_monodromy(&matrix, 2, 1.0, Some(&false_flow), controls(), || false).unwrap();
    assert_eq!(
        false_phase.classification,
        FloquetClassification::Unresolved
    );
    assert!(false_phase.phase_residual.unwrap() > controls().maximum_phase_residual);
    assert_eq!(false_phase.effective_spectral_radius, None);
    assert_eq!(false_phase.dominant_growth_rate_per_s, None);
    assert_eq!(false_phase.spectrum.values.len(), 2);

    let repeated = real(&[0.5, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
    let flow = real(&[0.0, 1.0, 0.0]);
    let ambiguous =
        analyze_monodromy(&repeated, 3, 1.0, Some(&flow), controls(), || false).unwrap();
    assert_eq!(ambiguous.classification, FloquetClassification::Unresolved);
    assert_eq!(ambiguous.phase_multiplier_index, None);
    assert_eq!(ambiguous.effective_spectral_radius, None);
    assert_eq!(ambiguous.dominant_growth_rate_per_s, None);
    assert_eq!(ambiguous.spectrum.values, real(&[0.5, 1.0, 1.0]));

    // One root is inside the tight phase tolerance; another is outside it but
    // inside the wider +1 isolation guard. It must not become a stable quotient.
    let near_neutral = real(&[1.0 - controls().unit_circle_guard / 2.0, 0.0, 0.0, 1.0]);
    let flow = real(&[0.0, 1.0]);
    let near_cluster =
        analyze_monodromy(&near_neutral, 2, 1.0, Some(&flow), controls(), || false).unwrap();
    assert_eq!(
        near_cluster.classification,
        FloquetClassification::Unresolved
    );
    assert_eq!(near_cluster.effective_spectral_radius, None);
    assert_eq!(near_cluster.spectrum.values.len(), 2);
}

#[test]
fn supplied_multipliers_require_complete_correct_independent_modes() {
    let matrix = rotating_map(0.75);
    let exact = rotating_multipliers(0.75);
    assert!(
        assess_monodromy_with_multipliers(&matrix, 2, 1.0, &exact, None, controls(), || false,)
            .is_ok()
    );
    assert!(matches!(
        assess_monodromy_with_multipliers(&matrix, 2, 1.0, &exact[..1], None, controls(), || false,),
        Err(FloquetError::Invalid(_))
            | Err(FloquetError::Spectrum(DenseSpectrumError::Invalid { .. }))
    ));
    for wrong in [real(&[0.4, 0.8]), vec![exact[0], exact[0]]] {
        assert!(matches!(
            assess_monodromy_with_multipliers(&matrix, 2, 1.0, &wrong, None, controls(), || false,),
            Err(FloquetError::Spectrum(_))
        ));
    }

    // Correct repeated roots of a defective map still do not supply an
    // independent complete modal basis, so no pruned stable answer is allowed.
    let jordan = real(&[0.5, 1.0, 0.0, 0.5]);
    let repeated = real(&[0.5, 0.5]);
    assert!(matches!(
        assess_monodromy_with_multipliers(&jordan, 2, 1.0, &repeated, None, controls(), || false,),
        Err(FloquetError::Spectrum(
            DenseSpectrumError::EigenpairResidual { .. }
        )) | Err(FloquetError::Spectrum(
            DenseSpectrumError::UnresolvedEigenbasis { .. }
        ))
    ));
}

#[test]
fn malformed_phase_limits_resource_caps_and_cancellation_refuse_atomically() {
    let matrix = rotating_map(0.75);
    let exact = rotating_multipliers(0.75);
    for period in [0.0, -1.0, f64::NAN] {
        assert!(matches!(
            analyze_monodromy(&matrix, 2, period, None, controls(), || false),
            Err(FloquetError::Invalid(_))
        ));
    }
    for flow in [
        real(&[0.0, 0.0]),
        real(&[1.0]),
        vec![C64::ONE, C64::new(f64::NAN, 0.0)],
    ] {
        assert!(matches!(
            analyze_monodromy(&matrix, 2, 1.0, Some(&flow), controls(), || false),
            Err(FloquetError::Invalid(_))
        ));
    }
    for invalid in [
        FloquetControls {
            unit_circle_guard: 0.0,
            ..controls()
        },
        FloquetControls {
            phase_multiplier_tolerance: controls().unit_circle_guard,
            ..controls()
        },
        FloquetControls {
            maximum_phase_residual: 1.0,
            ..controls()
        },
    ] {
        assert!(matches!(
            analyze_monodromy(&matrix, 2, 1.0, None, invalid, || false),
            Err(FloquetError::Invalid(_))
        ));
    }
    for spectrum in [
        DenseSpectrumOptions {
            maximum_dimension: 1,
            ..controls().spectrum
        },
        DenseSpectrumOptions {
            maximum_dense_bytes: 1,
            ..controls().spectrum
        },
        DenseSpectrumOptions {
            maximum_work: 0,
            ..controls().spectrum
        },
    ] {
        assert!(matches!(
            analyze_monodromy(
                &matrix,
                2,
                1.0,
                None,
                FloquetControls {
                    spectrum,
                    ..controls()
                },
                || false,
            ),
            Err(FloquetError::Spectrum(DenseSpectrumError::Invalid { .. }))
        ));
    }
    for stop in [1, 4, 12] {
        let mut polls = 0;
        assert!(matches!(
            analyze_monodromy(&matrix, 2, 1.0, None, controls(), || {
                polls += 1;
                polls >= stop
            }),
            Err(FloquetError::Cancelled)
        ));
        let mut polls = 0;
        assert!(matches!(
            assess_monodromy_with_multipliers(&matrix, 2, 1.0, &exact, None, controls(), || {
                polls += 1;
                polls >= stop
            },),
            Err(FloquetError::Cancelled)
        ));
    }
}
