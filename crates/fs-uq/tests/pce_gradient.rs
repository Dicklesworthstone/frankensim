//! Analytic surrogate derivatives, including zero tensor factors.
use fs_uq::pce::{PceModel, hermite_orthonormal};

fn model() -> PceModel {
    PceModel {
        dim: 3,
        indices: vec![vec![0, 0, 0], vec![1, 0, 0], vec![0, 2, 0],
                      vec![1, 1, 0], vec![1, 0, 3]],
        coefficients: vec![7.0, 2.0, -3.0, 4.0, 0.5],
    }
}

#[test]
fn closed_form_gradients_and_finite_differences_agree() {
    let model = model();
    for x in [[0.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.2, -0.7, 1.3],
              [1.0, 0.0, fs_math::det::sqrt(3.0)]] {
        let (value, gradient) = model.eval_with_gradient(&x);
        assert!((value - model.eval(&x)).abs() < 1e-13);
        let expected = [
            2.0 + 4.0 * x[1] + 0.5 * hermite_orthonormal(3, x[2]),
            -3.0 * fs_math::det::sqrt(2.0) * x[1] + 4.0 * x[0],
            0.5 * x[0] * fs_math::det::sqrt(3.0) * hermite_orthonormal(2, x[2]),
        ];
        for j in 0..3 {
            assert!((gradient[j] - expected[j]).abs() < 2e-13);
            let mut plus = x;
            let mut minus = x;
            plus[j] += 1e-5;
            minus[j] -= 1e-5;
            let fd = (model.eval(&plus) - model.eval(&minus)) / 2e-5;
            assert!((gradient[j] - fd).abs() < 2e-8);
        }
    }
}

#[test]
fn zero_basis_factors_do_not_destroy_nonzero_derivatives() {
    let model = PceModel { dim: 3, indices: vec![vec![1, 1, 0]], coefficients: vec![2.0] };
    let (value, gradient) = model.eval_with_gradient(&[0.0, 3.0, 17.0]);
    assert_eq!(value, 0.0);
    assert!((gradient[0] - 6.0).abs() < 1e-14);
    assert_eq!(gradient[1..], [0.0, 0.0]);
    assert_eq!(model.eval_with_gradient(&[0.0; 3]).1, vec![0.0; 3]);
}

#[test]
fn constant_surrogates_have_well_defined_zero_gradients() {
    for dim in [0, 3] {
        let model = PceModel {
            dim, indices: vec![vec![0; dim]], coefficients: vec![7.0],
        };
        let (value, gradient) = model.eval_with_gradient(&vec![0.0; dim]);
        assert!((value - 7.0).abs() < 1e-14);
        assert_eq!(gradient, vec![0.0; dim]);
    }
}

#[test]
fn malformed_evaluation_refuses_instead_of_silently_truncating() {
    assert!(std::panic::catch_unwind(|| model().eval_with_gradient(&[1.0, 2.0])).is_err());
    assert!(std::panic::catch_unwind(|| model().eval_with_gradient(&[f64::NAN; 3])).is_err());
    let mut malformed = model();
    malformed.indices[0].pop();
    assert!(std::panic::catch_unwind(|| malformed.eval_with_gradient(&[0.0; 3])).is_err());
    let mut malformed = model();
    malformed.coefficients.pop();
    assert!(std::panic::catch_unwind(|| malformed.eval_with_gradient(&[0.0; 3])).is_err());
}
