use super::*;

fn options() -> DenseSpectrumOptions {
    DenseSpectrumOptions {
        maximum_dimension: 16,
        maximum_dense_bytes: 1 << 20,
        maximum_work: 1_000_000_000,
        qr_max_sweeps_per_window: 160,
        svd_max_sweeps: 100,
        svd_tolerance: 2e-14,
        maximum_relative_residual: 1e-11,
        minimum_eigenbasis_reciprocal_condition: 1e-8,
    }
}
fn real(a: &[f64]) -> Vec<C64> {
    a.iter().copied().map(C64::from_re).collect()
}
fn near(a: C64, b: C64) {
    assert!((a - b).abs() < 2e-12 * (1. + b.abs()), "{a:?} != {b:?}");
}

#[test]
fn g1_nonsymmetric_complex_modes_keep_original_residuals_and_deterministic_vectors() {
    let a = real(&[-1., 6., -2., -3.]);
    let first = nonsymmetric_spectrum(&a, 2, options(), || false).unwrap();
    let second = nonsymmetric_spectrum(&a, 2, options(), || false).unwrap();
    near(first.values[0], C64::new(-2., -11_f64.sqrt()));
    near(first.values[1], C64::new(-2., 11_f64.sqrt()));
    assert!(first.relative_residuals.iter().all(|&v| v < 2e-13));
    assert!(first.eigenbasis_reciprocal_condition > 0.1);
    assert_eq!(first.values, second.values);
    assert_eq!(first.right_vectors, second.right_vectors);
    let bad =
        eigenpair_row_residuals(&a, 2, first.values[0], &[C64::ONE, C64::ONE], || false).unwrap();
    assert!(bad.iter().any(|&v| v > 0.1));
}

#[test]
fn g1_exact_semisimple_multiplicity_zero_and_unstable_roots_are_not_dropped() {
    let a = real(&[
        -2., 0., 0., 0., 0., -2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.,
    ]);
    let spectrum = nonsymmetric_spectrum(&a, 4, options(), || false).unwrap();
    assert_eq!(spectrum.values, real(&[-2., -2., 0., 1.]));
    assert_eq!(spectrum.right_vectors.len(), 4);
    assert_eq!(spectrum.relative_residuals, vec![0.; 4]);
    near(
        C64::from_re(spectrum.eigenbasis_reciprocal_condition),
        C64::ONE,
    );
    let inner: C64 = spectrum.right_vectors[0]
        .iter()
        .zip(&spectrum.right_vectors[1])
        .fold(C64::ZERO, |sum, (&a, &b)| sum + a.conj() * b);
    assert_eq!(inner, C64::ZERO);
}

#[test]
fn g0_defective_and_unresolved_bases_refuse_without_a_pruned_answer() {
    let jordan = real(&[-1., 1., 0., -1.]);
    assert!(matches!(
        nonsymmetric_spectrum(&jordan, 2, options(), || false),
        Err(DenseSpectrumError::EigenpairResidual { .. })
            | Err(DenseSpectrumError::UnresolvedEigenbasis { .. })
    ));
    let near_jordan = real(&[0., 1., 1e-10, 0.]);
    let strict = DenseSpectrumOptions {
        minimum_eigenbasis_reciprocal_condition: 1e-3,
        ..options()
    };
    assert!(matches!(
        nonsymmetric_spectrum(&near_jordan, 2, strict, || false),
        Err(DenseSpectrumError::UnresolvedEigenbasis { .. })
    ));
    assert!(nonsymmetric_spectrum(&near_jordan, 2, options(), || false).is_ok());
}

#[test]
fn g4_complete_capacity_and_native_cancellation_do_not_publish_partial_modes() {
    let a = real(&[-1., 6., -2., -3.]);
    for budget in [
        DenseSpectrumOptions {
            maximum_dimension: 1,
            ..options()
        },
        DenseSpectrumOptions {
            maximum_dense_bytes: 1,
            ..options()
        },
        DenseSpectrumOptions {
            maximum_work: 0,
            ..options()
        },
    ] {
        assert!(matches!(
            nonsymmetric_spectrum(&a, 2, budget, || false),
            Err(DenseSpectrumError::Invalid { .. })
        ));
    }
    for stop in [1, 3, 10] {
        let mut calls = 0;
        assert!(matches!(
            nonsymmetric_spectrum(&a, 2, options(), || {
                calls += 1;
                calls >= stop
            }),
            Err(DenseSpectrumError::Cancelled)
        ));
    }
}
