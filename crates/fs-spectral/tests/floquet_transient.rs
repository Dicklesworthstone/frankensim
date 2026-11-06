//! Independent finite-horizon oracles, including defective stable maps.

use fs_la::eigen_complex::ComplexSvdError;
use fs_math::c64::C64;
use fs_spectral::floquet::transient::{
    StroboscopicGrowthControls, StroboscopicGrowthError, StroboscopicGrowthSample,
    stroboscopic_growth,
};

fn controls() -> StroboscopicGrowthControls {
    StroboscopicGrowthControls {
        maximum_dimension: 16,
        maximum_horizons: 64,
        maximum_periods: u64::MAX,
        maximum_payload_bytes: 8 << 20,
        maximum_work: 5_000_000_000,
        svd_max_sweeps: 100,
        svd_tolerance: 2e-14,
        maximum_relative_residual: 2e-11,
    }
}

fn real(values: &[f64]) -> Vec<C64> {
    values.iter().copied().map(C64::from_re).collect()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-10 * expected.abs().max(1e-300),
        "{actual:.16e} != {expected:.16e}"
    );
}

// Check the claimed directions against the independently authored ACTUAL power,
// not the implementation's normalized matrix or its own residual fields.
fn directions(sample: &StroboscopicGrowthSample, power: &[C64], n: usize, gain: f64) {
    let initial = &sample.initial_direction;
    let final_direction = sample.final_direction.as_ref().unwrap();
    assert_eq!(initial.len(), n);
    assert_eq!(final_direction.len(), n);
    let norm = |v: &[C64]| v.iter().fold(0.0_f64, |sum, z| sum.hypot(z.abs()));
    near(norm(initial), 1.0);
    near(norm(final_direction), 1.0);
    let tolerance = 5e-10 * gain.max(1e-300);
    for row in 0..n {
        let image = (0..n).fold(C64::ZERO, |sum, j| sum + power[row * n + j] * initial[j]);
        assert!(
            (image - final_direction[row].scale(gain)).abs() < tolerance,
            "forward singular direction failed at horizon {} row {row}",
            sample.periods,
        );
        let adjoint = (0..n).fold(C64::ZERO, |sum, j| {
            sum + power[j * n + row].conj() * final_direction[j]
        });
        assert!(
            (adjoint - initial[row].scale(gain)).abs() < tolerance,
            "adjoint singular direction failed at horizon {} row {row}",
            sample.periods,
        );
    }
    assert!(sample.forward_relative_residual <= controls().maximum_relative_residual);
    assert!(sample.adjoint_relative_residual <= controls().maximum_relative_residual);
}

#[test]
fn normal_damped_complex_rotation_reduces_to_modal_gain_and_full_directions() {
    let radius = 0.8;
    let matrix = vec![
        C64::from_re(0.6 * radius),
        C64::new(0.0, 0.8 * radius),
        C64::new(0.0, 0.8 * radius),
        C64::from_re(0.6 * radius),
    ];
    let horizons = [0, 1, 2, 7, 31];
    let result = stroboscopic_growth(&matrix, 2, &horizons, controls(), || false).unwrap();
    assert_eq!(result.dimension, 2);
    assert_eq!(result.samples.len(), horizons.len());
    let theta = 0.8_f64.atan2(0.6);
    for (sample, &k) in result.samples.iter().zip(&horizons) {
        assert_eq!(sample.periods, k);
        let gain = radius.powi(k as i32);
        near(sample.gain.unwrap(), gain);
        if k == 0 {
            assert_eq!(sample.log_gain, Some(0.0));
        } else {
            near(sample.log_gain.unwrap(), k as f64 * radius.ln());
        }
        let (sine, cosine) = (k as f64 * theta).sin_cos();
        let exact = vec![
            C64::from_re(gain * cosine),
            C64::new(0.0, gain * sine),
            C64::new(0.0, gain * sine),
            C64::from_re(gain * cosine),
        ];
        directions(sample, &exact, 2, gain);
    }
}

#[test]
fn stable_jordan_map_exposes_transient_amplification_then_decay_without_eigenvectors() {
    let radius = 0.75_f64;
    let coupling = 4.0;
    let matrix = real(&[radius, coupling, 0.0, radius]);
    let horizons = [0, 1, 2, 4, 8, 16, 32, 64];
    let result = stroboscopic_growth(&matrix, 2, &horizons, controls(), || false).unwrap();
    for (sample, &k) in result.samples.iter().zip(&horizons) {
        // M^k = [rho^k, k*b*rho^(k-1); 0, rho^k]. Its top singular value
        // is (hypot(2*a,b_k)+|b_k|)/2, independent of any numerical SVD.
        let diagonal = radius.powi(k as i32);
        let off = if k == 0 {
            0.0
        } else {
            k as f64 * coupling * radius.powi(k as i32 - 1)
        };
        let gain = ((2.0 * diagonal).hypot(off) + off.abs()) / 2.0;
        near(sample.gain.unwrap(), gain);
        if k == 0 {
            assert_eq!(sample.log_gain, Some(0.0));
        } else {
            near(sample.log_gain.unwrap(), gain.ln());
        }
        directions(sample, &real(&[diagonal, off, 0.0, diagonal]), 2, gain);
    }
    assert!(result.samples[1].gain.unwrap() > 4.0);
    assert!(result.samples[3].gain.unwrap() > result.samples[1].gain.unwrap());
    assert!(result.samples.last().unwrap().gain.unwrap() < 1e-4);
    // Every asymptotic multiplier is still exactly 0.75. No complete
    // eigenvector basis exists, and none is needed for these finite powers.
}

#[test]
fn sparse_huge_horizons_keep_logarithmic_gain_when_ordinary_gain_is_out_of_range() {
    let bounded = StroboscopicGrowthControls {
        maximum_work: 1_000_000,
        ..controls()
    };
    let horizons = [0, 2, 1_000, 1_000_000, u64::MAX];
    let growing = stroboscopic_growth(&real(&[1.25]), 1, &horizons, bounded, || false).unwrap();
    assert!(growing.matrix_products < 400);
    assert!(growing.admitted_work <= bounded.maximum_work);
    for (sample, &k) in growing.samples.iter().zip(&horizons) {
        if k == 0 {
            assert_eq!(sample.log_gain, Some(0.0));
        } else {
            near(sample.log_gain.unwrap(), k as f64 * 1.25_f64.ln());
        }
        near(sample.initial_direction[0].abs(), 1.0);
        near(sample.final_direction.as_ref().unwrap()[0].abs(), 1.0);
        if k <= 1_000 {
            near(sample.gain.unwrap(), 1.25_f64.powi(k as i32));
        } else {
            assert_eq!(sample.gain, None);
            assert!(sample.log_gain.unwrap().is_finite());
        }
    }
    let decaying =
        stroboscopic_growth(&real(&[0.5]), 1, &[0, 1_000_000, u64::MAX], bounded, || {
            false
        })
        .unwrap();
    for sample in decaying.samples.iter().skip(1) {
        assert_eq!(sample.gain, None);
        near(
            sample.log_gain.unwrap(),
            sample.periods as f64 * 0.5_f64.ln(),
        );
        assert!(sample.final_direction.is_some());
    }
    // Scaling must happen before squaring: these ordinary powers cannot be
    // represented even though their logarithms and optimal directions can.
    for value in [1e300_f64, 1e-300] {
        let result = stroboscopic_growth(&real(&[value]), 1, &[2], bounded, || false).unwrap();
        let sample = &result.samples[0];
        assert_eq!(sample.gain, None);
        near(sample.log_gain.unwrap(), 2.0 * value.ln());
        assert!(sample.final_direction.is_some());
    }
}

#[test]
fn true_nilpotent_zero_powers_are_explicit_and_lost_nonzero_components_refuse() {
    let nilpotent = vec![C64::ZERO, C64::new(2.0, 1.0), C64::ZERO, C64::ZERO];
    let result = stroboscopic_growth(&nilpotent, 2, &[0, 1, 2, 37], controls(), || false).unwrap();
    assert_eq!(result.samples[0].gain, Some(1.0));
    let first = &result.samples[1];
    near(first.gain.unwrap(), 5.0_f64.sqrt());
    directions(first, &nilpotent, 2, 5.0_f64.sqrt());
    for sample in result.samples.iter().skip(2) {
        assert_eq!(sample.gain, Some(0.0));
        assert_eq!(sample.log_gain, None);
        assert_eq!(sample.final_direction, None);
        assert_eq!(sample.initial_direction.len(), 2);
        near(
            sample
                .initial_direction
                .iter()
                .fold(0.0_f64, |s, z| s.hypot(z.abs())),
            1.0,
        );
        assert_eq!(sample.forward_relative_residual, 0.0);
        assert_eq!(sample.adjoint_relative_residual, 0.0);
    }
    let incompatible_range = real(&[1e300, 0.0, 0.0, 1e-300]);
    assert!(matches!(
        stroboscopic_growth(&incompatible_range, 2, &[1], controls(), || false),
        Err(StroboscopicGrowthError::Arithmetic(_))
    ));
    // M^2 has one NONZERO entry eps^2; losing that product must not fabricate
    // a zero gain. A single global matrix scale cannot resolve this range.
    let eps = 1e-200;
    let tiny_product = real(&[0.0, 1.0, eps, 0.0, 0.0, 0.0, 0.0, eps, 0.0]);
    assert!(matches!(
        stroboscopic_growth(&tiny_product, 3, &[2], controls(), || false),
        Err(StroboscopicGrowthError::Arithmetic(_))
    ));
}

#[test]
fn malformed_queries_total_limits_native_iteration_limits_and_cancellation_refuse() {
    let matrix = real(&[0.75, 4.0, 0.0, 0.75]);
    for horizons in [&[][..], &[1, 1][..], &[3, 2][..]] {
        assert!(matches!(
            stroboscopic_growth(&matrix, 2, horizons, controls(), || false),
            Err(StroboscopicGrowthError::Invalid(_))
        ));
    }
    assert!(matches!(
        stroboscopic_growth(&matrix[..3], 2, &[1], controls(), || false),
        Err(StroboscopicGrowthError::Invalid(_))
    ));
    assert!(matches!(
        stroboscopic_growth(&[C64::new(f64::INFINITY, 0.0)], 1, &[1], controls(), || {
            false
        }),
        Err(StroboscopicGrowthError::Invalid(_))
    ));
    for limited in [
        StroboscopicGrowthControls {
            maximum_dimension: 1,
            ..controls()
        },
        StroboscopicGrowthControls {
            maximum_horizons: 1,
            ..controls()
        },
        StroboscopicGrowthControls {
            maximum_periods: 0,
            ..controls()
        },
        StroboscopicGrowthControls {
            maximum_work: 0,
            ..controls()
        },
        StroboscopicGrowthControls {
            maximum_payload_bytes: 1,
            ..controls()
        },
    ] {
        assert!(matches!(
            stroboscopic_growth(&matrix, 2, &[0, 1], limited, || false),
            Err(StroboscopicGrowthError::Budget { .. })
        ));
    }
    let one = stroboscopic_growth(&matrix, 2, &[3], controls(), || false).unwrap();
    // Each individual horizon fits; the complete three-SVD campaign does not.
    // Limits cannot reset for each returned sample or omit retained directions.
    for limited in [
        StroboscopicGrowthControls {
            maximum_work: one.admitted_work,
            ..controls()
        },
        StroboscopicGrowthControls {
            maximum_payload_bytes: one.admitted_payload_bytes,
            ..controls()
        },
    ] {
        assert!(matches!(
            stroboscopic_growth(&matrix, 2, &[1, 2, 3], limited, || false),
            Err(StroboscopicGrowthError::Budget { .. })
        ));
    }
    let dense = real(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0]);
    assert!(matches!(
        stroboscopic_growth(
            &dense,
            3,
            &[1],
            StroboscopicGrowthControls {
                svd_max_sweeps: 1,
                ..controls()
            },
            || false,
        ),
        Err(StroboscopicGrowthError::SingularVectors(
            ComplexSvdError::NotConverged { .. }
        ))
    ));
    for stop in [1, 5, 30, 100] {
        let mut polls = 0;
        assert!(matches!(
            stroboscopic_growth(&matrix, 2, &[0, 3, 17], controls(), || {
                polls += 1;
                polls >= stop
            }),
            Err(StroboscopicGrowthError::Cancelled)
        ));
    }
}
