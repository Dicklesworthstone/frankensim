//! High-dimensional PCE regression, basis preflight, and input refusals.
use fs_uq::pce::{fit_pce, hermite_orthonormal, total_degree_basis_size};

fn response(x: &[f64]) -> f64 {
    3.0 + 0.7 * x[0] - 1.2 * x[19]
        + 0.4 * x[2] * x[11] + 0.3 * hermite_orthonormal(2, x[5])
}

#[test]
fn twenty_germ_quadratic_fit_recovers_a_known_surrogate() {
    let dim = 20;
    let mut xi = vec![vec![0.0; dim]];
    for j in 0..dim {
        for sign in [-1.0, 1.0] {
            let mut row = vec![0.0; dim];
            row[j] = sign;
            xi.push(row);
        }
        for k in j + 1..dim {
            let mut row = vec![0.0; dim];
            row[j] = 1.0;
            row[k] = 1.0;
            xi.push(row);
        }
    }
    // This unisolvent design has one row per quadratic basis term. Repeat
    // it to obey the existing 2*basis regression admission requirement.
    assert_eq!(Some(xi.len()), total_degree_basis_size(dim, 2));
    xi.extend_from_within(..);
    let y: Vec<_> = xi.iter().map(|row| response(row)).collect();
    let model = fit_pce(&xi, &y, 2);
    assert_eq!(model.indices.len(), 231);
    for sample in 0..12 {
        let x: Vec<_> = (0..dim).map(|j| {
            ((sample * 17 + j * 7) % 31) as f64 / 10.0 - 1.5
        }).collect();
        assert!((model.eval(&x) - response(&x)).abs() < 1e-7);
        let (value, gradient) = model.eval_with_gradient(&x);
        assert!((value - response(&x)).abs() < 1e-7);
        for j in 0..dim {
            let expected = match j {
                0 => 0.7,
                2 => 0.4 * x[11],
                5 => 0.3 * fs_math::det::sqrt(2.0) * x[5],
                11 => 0.4 * x[2],
                19 => -1.2,
                _ => 0.0,
            };
            assert!((gradient[j] - expected).abs() < 1e-7);
        }
    }
    // The fitted surrogate feeds the global-sensitivity API directly.
    let sensitivity = model.sobol_indices().unwrap();
    let variance = 0.7 * 0.7 + 1.2 * 1.2 + 0.4 * 0.4 + 0.3 * 0.3;
    for j in 0..dim {
        let main = match j { 0 => 0.49, 5 => 0.09, 19 => 1.44, _ => 0.0 };
        let total = main + if j == 2 || j == 11 { 0.16 } else { 0.0 };
        assert!((sensitivity.first_order[j] - main / variance).abs() < 1e-7);
        assert!((sensitivity.total_order[j] - total / variance).abs() < 1e-7);
    }
    let replay = fit_pce(&xi, &y, 2);
    assert_eq!(model.indices, replay.indices);
    assert!(model.coefficients.iter().zip(&replay.coefficients)
        .all(|(a, b)| a.to_bits() == b.to_bits()));
}

#[test]
fn malformed_samples_refuse_instead_of_changing_the_regression() {
    let refuses = |xi: Vec<Vec<f64>>, y: Vec<f64>, p| {
        assert!(std::panic::catch_unwind(|| fit_pce(&xi, &y, p)).is_err());
    };
    refuses(vec![], vec![], 1);
    refuses(vec![vec![0.0]; 2], vec![0.0], 0);
    refuses(vec![vec![0.0], vec![0.0, 1.0]], vec![0.0; 2], 0);
    refuses(vec![vec![f64::NAN]; 2], vec![0.0; 2], 0);
    refuses(vec![vec![0.0]; 2], vec![f64::INFINITY; 2], 0);
    // Refuse before the old multi-billion-candidate enumeration.
    refuses(vec![vec![0.0; 20]; 2], vec![0.0; 2], 2);
    refuses(vec![vec![0.0]; 2], vec![0.0; 2], usize::MAX);
    refuses(vec![vec![f64::MAX]; 6], vec![0.0; 6], 2);
}

#[test]
fn zero_germ_regression_retains_the_constant_model() {
    let model = fit_pce(&[vec![], vec![]], &[4.0, 4.0], usize::MAX);
    assert_eq!(model.indices, vec![Vec::<usize>::new()]);
    assert!((model.mean() - 4.0).abs() < 1e-10);
    assert!((model.eval(&[]) - 4.0).abs() < 1e-10);
}
