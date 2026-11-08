//! Physical residuals must qualify the returned, unit-restored load factor.
use fs_solid::SolidError;
use fs_solid::stability::{
    self, BucklingResult,
    sparse::{self, SparseBucklingLimits},
};
use fs_sparse::{Coo, Csr};

fn diagonal(values: &[f64]) -> Csr {
    let mut coo = Coo::new(values.len(), values.len());
    for (i, &value) in values.iter().enumerate() {
        coo.push(i, i, value);
    }
    coo.assemble()
}

fn solve(k: &Csr, g: &Csr, sparse_route: bool, count: usize) -> Result<BucklingResult, SolidError> {
    let map: Vec<_> = (0..k.nrows()).map(Some).collect();
    if sparse_route {
        sparse::buckling_loads(k, g, &map, count, 150, &SparseBucklingLimits::default())
    } else {
        stability::buckling_loads(k, g, &map, count, 150)
    }
}

fn original_residual(k: &Csr, g: &Csr, load: f64, mode: &[f64]) -> f64 {
    let mut kp = vec![0.0; mode.len()];
    let mut gp = vec![0.0; mode.len()];
    k.spmv(mode, &mut kp);
    g.spmv(mode, &mut gp);
    let (mut residual, mut knorm, mut gnorm) = (0.0_f64, 0.0_f64, 0.0_f64);
    for (&a, &b) in kp.iter().zip(&gp) {
        residual = residual.hypot(load.mul_add(b, a));
        knorm = knorm.hypot(a);
        gnorm = gnorm.hypot(load * b);
    }
    residual / (knorm + gnorm)
}

#[test]
fn coarse_subnormal_load_refuses_even_with_an_exact_transformed_root() {
    let k = diagonal(&[1e-300; 9]);
    let g = diagonal(&[-1e20; 9]);
    // The transformed operator is exactly I, but the physical root 1e-320
    // has too few significant bits: its rounded value misses the original
    // scalar equation by ~5.6e-6 in the declared relative residual metric.
    let rounded = 1e-300 / 1e20;
    let scalar_residual = (1e-300_f64 - rounded * 1e20).abs() / (1e-300_f64 + rounded * 1e20);
    assert!(scalar_residual > 1e-6);
    for sparse_route in [false, true] {
        assert!(
            matches!(
                solve(&k, &g, sparse_route, 1),
                Err(SolidError::SolveFailed { rel_residual, .. })
                    if rel_residual.is_finite() && rel_residual > 1e-8
            ),
            "route sparse={sparse_route}"
        );
    }
}

#[test]
fn exactly_representable_subnormal_loads_pass_the_actual_physical_gate() {
    let stiffness = f64::from_bits(23_u64 << 52); // 2^-1000
    let reference = f64::from_bits(1073_u64 << 52); // 2^50
    let expected = f64::from_bits(1_u64 << 24); // 2^-1050, exactly subnormal
    let k = diagonal(&[stiffness; 9]);
    let g = diagonal(&[-reference; 9]);
    for sparse_route in [false, true] {
        let result = solve(&k, &g, sparse_route, 1).unwrap();
        assert_eq!(result.loads, vec![expected]);
        assert!(result.relative_residuals[0] <= 1e-8);
        assert!(original_residual(&k, &g, result.loads[0], &result.modes[0]) <= 1e-8);
    }
}

#[test]
fn subnormal_intermediate_ratio_does_not_corrupt_normal_physical_loads() {
    let k = diagonal(&[1e-300; 3]);
    let g = diagonal(&[-1e-80, -5e-81, 1e20]);
    // Request the complete 3x3 dense native pencil, avoiding an iterative
    // relative-resolution claim for tiny positive roots beside strong tension.
    // K/G_scale is subnormal, but dividing by theta=1e-100 restores a normal
    // physical load. Both scalar equations have independent known roots.
    let result = solve(&k, &g, false, 3).unwrap();
    assert_eq!(result.loads.len(), 2);
    for (i, expected) in [1e-220, 2e-220].into_iter().enumerate() {
        assert!((result.loads[i] / expected - 1.0).abs() < 2e-13);
        assert!(result.relative_residuals[i] <= 1e-8);
        assert!(original_residual(&k, &g, result.loads[i], &result.modes[i]) <= 1e-8);
    }
}
