//! CPU differential, equivalence-envelope, determinism, and fault test suite for accelerator pilot
//! (bead `frankensim-extreal-program-f85xj.15.3.3`).

use fs_sparse::{Coo, Csr, run_accelerator_spmv_pilot};

#[test]
fn test_differential_empty_matrix() {
    let coo = Coo::new(0, 0);
    let csr = coo.assemble();
    let x: Vec<f64> = vec![];
    let (y, receipt) = run_accelerator_spmv_pilot(&csr, &x, "test_empty");
    assert_eq!(y.len(), 0);
    assert!(receipt.envelope.passed);
    assert!(receipt.cancellation_drain_verified);
}

#[test]
fn test_differential_diagonal_matrix() {
    let mut coo = Coo::new(100, 100);
    for i in 0..100 {
        coo.push(i, i, (i + 1) as f64 * 1.5);
    }
    let csr = coo.assemble();
    let x: Vec<f64> = (0..100).map(|i| i as f64).collect();
    let (y, receipt) = run_accelerator_spmv_pilot(&csr, &x, "test_diag");

    assert_eq!(y.len(), 100);
    for i in 0..100 {
        let expected = (i + 1) as f64 * 1.5 * (i as f64);
        assert!((y[i] - expected).abs() < 1e-12);
    }
    assert!(receipt.envelope.passed);
    assert_eq!(receipt.matrix_shape, (100, 100, 100));
}

#[test]
fn test_differential_laplacian_1d_determinism() {
    // 1D tridiagonal Laplacian: [-1, 2, -1]
    let n = 200;
    let mut coo = Coo::new(n, n);
    for i in 0..n {
        coo.push(i, i, 2.0);
        if i > 0 {
            coo.push(i, i - 1, -1.0);
        }
        if i + 1 < n {
            coo.push(i, i + 1, -1.0);
        }
    }
    let csr = coo.assemble();
    let x: Vec<f64> = (0..n).map(|i| (i as f64).sin()).collect();

    // Run 1
    let (y1, r1) = run_accelerator_spmv_pilot(&csr, &x, "run_laplacian_1");
    // Run 2
    let (y2, r2) = run_accelerator_spmv_pilot(&csr, &x, "run_laplacian_2");

    assert_eq!(y1, y2, "Fixed-order reduction must be bit-identical");
    assert!(r1.envelope.passed);
    assert!(r2.envelope.passed);
    assert_eq!(r1.matrix_shape.2, 3 * n - 2);
}

#[test]
fn test_differential_extreme_dynamic_range() {
    let mut coo = Coo::new(3, 3);
    coo.push(0, 0, 1e20);
    coo.push(0, 1, 1e-20);
    coo.push(1, 1, -1e15);
    coo.push(2, 2, 4.2);
    let csr = coo.assemble();

    let x = vec![1.0, 1e20, 2.0];
    let (y, receipt) = run_accelerator_spmv_pilot(&csr, &x, "test_extreme");

    assert_eq!(y[0], 1e20 + 1.0);
    assert_eq!(y[1], -1e35);
    assert!((y[2] - 8.4).abs() < 1e-12);
    assert!(receipt.envelope.passed);
}
