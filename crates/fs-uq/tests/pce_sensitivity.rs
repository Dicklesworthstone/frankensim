//! Known-answer Sobol decomposition, scale invariance and public-field errors.
use fs_uq::pce::{PceModel, PceSensitivityError};

fn example() -> PceModel {
    PceModel {
        dim: 4,
        indices: vec![vec![0, 0, 0, 0], vec![1, 0, 0, 0], vec![0, 1, 0, 0],
                      vec![1, 1, 0, 0], vec![0, 0, 2, 0]],
        coefficients: vec![10.0, 2.0, 3.0, 4.0, 5.0],
    }
}

fn close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (&a, &b) in actual.iter().zip(expected) {
        assert!((a - b).abs() <= 2e-14, "{a} != {b}");
    }
}

#[test]
fn mixed_main_effects_interaction_and_inactive_germ() {
    let report = example().sobol_indices().unwrap();
    close(&report.first_order, &[4.0 / 54.0, 9.0 / 54.0, 25.0 / 54.0, 0.0]);
    close(&report.total_order, &[20.0 / 54.0, 25.0 / 54.0, 25.0 / 54.0, 0.0]);
    assert_eq!(report.components.iter().map(|c| c.variables.clone()).collect::<Vec<_>>(),
               vec![vec![0], vec![0, 1], vec![1], vec![2]]);
    close(&report.components.iter().map(|c| c.index).collect::<Vec<_>>(),
          &[4.0 / 54.0, 16.0 / 54.0, 9.0 / 54.0, 25.0 / 54.0]);
    assert!((report.components.iter().map(|c| c.index).sum::<f64>() - 1.0).abs() < 2e-14);
}

#[test]
fn polynomial_degree_is_not_interaction_order() {
    let model = PceModel {
        dim: 3,
        indices: vec![vec![1, 0, 0], vec![5, 0, 0], vec![1, 2, 3]],
        coefficients: vec![3.0, 4.0, 5.0],
    };
    let report = model.sobol_indices().unwrap();
    close(&report.first_order, &[0.5, 0.0, 0.0]);
    close(&report.total_order, &[1.0, 0.5, 0.5]);
    assert_eq!(report.components.len(), 2);
    assert_eq!(report.components[1].variables, vec![0, 1, 2]);
}

#[test]
fn output_units_constant_offset_and_term_permutation_do_not_change_indices() {
    let expected = example().sobol_indices().unwrap();
    for scale in [1e300, 1e-300, -1e300, -1e-300] {
        let mut model = example();
        for c in &mut model.coefficients {
            *c *= scale;
        }
        // The constant is not included in the normalization scale.
        model.coefficients[0] = f64::MAX;
        let report = model.sobol_indices().unwrap();
        close(&report.first_order, &expected.first_order);
        close(&report.total_order, &expected.total_order);
    }
    let mut permuted = example();
    permuted.indices.reverse();
    permuted.coefficients.reverse();
    assert_eq!(permuted.sobol_indices().unwrap(), expected);
}

#[test]
fn high_dimensional_supports_are_not_limited_by_integer_bit_masks() {
    let mut alpha = vec![0; 130];
    alpha[0] = 2;
    alpha[129] = 1;
    let model = PceModel { dim: 130, indices: vec![alpha], coefficients: vec![1.0] };
    let report = model.sobol_indices().unwrap();
    assert!(report.first_order.iter().all(|&v| v == 0.0));
    close(&[report.total_order[0], report.total_order[129]], &[1.0, 1.0]);
    assert_eq!(report.components[0].variables, vec![0, 129]);
}

#[test]
fn undefined_or_malformed_models_return_errors() {
    use PceSensitivityError::{CoefficientCount, DuplicateBasis, EmptyModel,
                             NonFiniteCoefficient, TermDimension, ZeroVariance};
    let mut model = example();
    model.coefficients.pop();
    assert!(matches!(model.sobol_indices(), Err(CoefficientCount { .. })));
    let mut model = example();
    model.indices[1].pop();
    assert!(matches!(model.sobol_indices(), Err(TermDimension { term: 1, .. })));
    let mut model = example();
    model.indices[2] = model.indices[1].clone();
    assert_eq!(model.sobol_indices().unwrap_err(), DuplicateBasis { first: 1, second: 2 });
    for non_finite in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut model = example();
        model.coefficients[0] = non_finite;
        assert_eq!(model.sobol_indices().unwrap_err(), NonFiniteCoefficient { term: 0 });
    }
    let mut model = example();
    model.coefficients[1..].fill(0.0);
    assert_eq!(model.sobol_indices().unwrap_err(), ZeroVariance);
    let model = PceModel { dim: 0, indices: vec![vec![]], coefficients: vec![7.0] };
    assert_eq!(model.sobol_indices().unwrap_err(), ZeroVariance);
    let model = PceModel { dim: 1, indices: vec![], coefficients: vec![] };
    assert_eq!(model.sobol_indices().unwrap_err(), EmptyModel);
}
