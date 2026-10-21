use fs_uq::{PceBasis, ProductPceModel};
use fs_uq::pce::{PceModel, PceSensitivityError};

fn close(a: f64, b: f64) {
    assert!((a-b).abs() < 2e-11 * (1.0+a.abs().max(b.abs())), "{a} != {b}");
}
fn legendre(degree: usize) -> ProductPceModel {
    if degree == 0 {
        ProductPceModel::new(vec![PceBasis::Uniform01], vec![vec![0]], vec![1.0]).unwrap()
    } else {
        ProductPceModel::new(vec![PceBasis::Uniform01], vec![vec![0], vec![degree]], vec![0.0, 1.0]).unwrap()
    }
}

#[test]
fn shifted_legendre_is_orthonormal_and_has_regular_endpoint_derivatives() {
    // Nine-point Gauss-Legendre integrates all products through degree sixteen.
    let nodes = [
        -0.9681602395076261, -0.8360311073266358, -0.6133714327005904,
        -0.3242534234038089, 0.0, 0.3242534234038089, 0.6133714327005904,
        0.8360311073266358, 0.9681602395076261,
    ];
    let weights = [
        0.0812743883615744, 0.1806481606948574, 0.2606106964029354,
        0.3123470770400029, 0.3302393550012598, 0.3123470770400029,
        0.2606106964029354, 0.1806481606948574, 0.0812743883615744,
    ];
    for n in 0..=8 {
        let a = legendre(n);
        for m in 0..=8 {
            let b = legendre(m);
            let integral = nodes.iter().zip(weights).map(|(&x, w)| {
                0.5*w*a.eval(&[(x+1.0)/2.0])*b.eval(&[(x+1.0)/2.0])
            }).sum();
            close(integral, if n == m { 1.0 } else { 0.0 });
        }
        let scale = fs_math::det::sqrt((2*n+1) as f64);
        for u in [0.0, 1.0] {
            let (value, gradient) = a.eval_with_gradient(&[u]);
            let sign = if u == 0.0 && n % 2 == 1 { -1.0 } else { 1.0 };
            close(value, sign*scale);
            let derivative_sign = if u == 0.0 { -sign } else { sign };
            close(gradient[0], derivative_sign*scale*(n*(n+1)) as f64);
        }
    }
}

#[test]
fn mixed_tensor_has_correct_moments_sensitivity_and_gradients_at_basis_zeros() {
    let model = ProductPceModel::new(
        vec![PceBasis::StandardNormal, PceBasis::Uniform01],
        vec![vec![0,0], vec![1,0], vec![0,1], vec![1,1]],
        vec![10.0, 2.0, 3.0, 4.0],
    ).unwrap();
    close(model.mean(), 10.0);
    close(model.variance(), 29.0);
    let indices = model.sobol_indices().unwrap();
    close(indices.first_order[0], 4.0/29.0);
    close(indices.first_order[1], 9.0/29.0);
    close(indices.total_order[0], 20.0/29.0);
    close(indices.total_order[1], 25.0/29.0);
    for z in [0.0, -0.75, 2.0] {
        for u in [0.0, 0.5, 1.0] {
            let l = fs_math::det::sqrt(3.0)*(2.0*u-1.0);
            let (value, gradient) = model.eval_with_gradient(&[z,u]);
            close(value, 10.0+2.0*z+3.0*l+4.0*z*l);
            close(gradient[0], 2.0+4.0*l);
            close(gradient[1], 2.0*fs_math::det::sqrt(3.0)*(3.0+4.0*z));
            close(value, model.eval(&[z,u]));
        }
    }
}

#[test]
fn gaussian_only_model_preserves_legacy_values_gradients_and_sensitivity() {
    let legacy = PceModel {
        dim: 2,
        indices: vec![vec![0,0], vec![0,1], vec![1,0], vec![2,0], vec![1,2]],
        coefficients: vec![2.0, -3.0, 0.5, 0.125, 0.25],
    };
    let model = ProductPceModel::new(
        vec![PceBasis::StandardNormal; 2], legacy.indices.clone(), legacy.coefficients.clone(),
    ).unwrap();
    for x in [[0.0,0.0], [0.5,-1.0], [-1.25,2.0], [3.0,-0.25]] {
        assert_eq!(model.eval(&x).to_bits(), legacy.eval(&x).to_bits());
        assert_eq!(model.eval_with_gradient(&x), legacy.eval_with_gradient(&x));
    }
    assert_eq!(model.sobol_indices(), legacy.sobol_indices());
    assert_eq!(model.indices(), legacy.indices.as_slice());
    assert_eq!(model.coefficients(), legacy.coefficients.as_slice());
}

#[test]
fn invalid_models_germs_and_unrepresentable_evaluations_refuse() {
    let uniform = vec![PceBasis::Uniform01];
    for (indices, coefficients) in [
        (vec![], vec![]), (vec![vec![0]], vec![]),
        (vec![vec![0]], vec![f64::NAN]), (vec![vec![1]], vec![1.0]),
        (vec![vec![0],vec![0]], vec![1.0,2.0]),
        (vec![vec![0],vec![9]], vec![1.0,2.0]),
        (vec![vec![0,0]], vec![1.0]),
    ] {
        assert!(ProductPceModel::new(uniform.clone(), indices, coefficients).is_err());
    }
    assert!(ProductPceModel::new(vec![PceBasis::Uniform01;17], vec![vec![0;17]], vec![1.0]).is_err());
    let constant = legendre(0);
    assert_eq!(constant.sobol_indices(), Err(PceSensitivityError::ZeroVariance));
    for x in [-0.01, 1.01, f64::NAN, f64::INFINITY] {
        assert!(constant.try_eval(&[x]).is_err());
        assert!(constant.try_eval_with_gradient(&[x]).is_err());
    }
    assert!(constant.try_eval(&[]).is_err());
    let gaussian = ProductPceModel::new(
        vec![PceBasis::StandardNormal], vec![vec![0],vec![2]], vec![0.0,1.0],
    ).unwrap();
    assert!(gaussian.try_eval(&[f64::MAX]).is_err());
    assert!(gaussian.try_eval_with_gradient(&[f64::MAX]).is_err());
}
