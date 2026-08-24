//! Wasm SIMD128 Tier-1w capsule conformance and scalar referee equivalence (bead `frankensim-wf-root-guzez.1.5`).

use fs_simd::scalar;
use fs_simd::wasm;

#[test]
fn wasm_simd128_axpy_matches_scalar() {
    let a = 2.5;
    let x = [1.0, 2.0, 3.0, 4.0, 5.0];
    let mut y_wasm = [10.0, 20.0, 30.0, 40.0, 50.0];
    let mut y_scalar = y_wasm;

    wasm::axpy(a, &x, &mut y_wasm);
    scalar::axpy(a, &x, &mut y_scalar);

    for (w, s) in y_wasm.iter().zip(&y_scalar) {
        assert_eq!(w.to_bits(), s.to_bits(), "axpy must match scalar bitwise");
    }
}

#[test]
fn wasm_simd128_scale_matches_scalar() {
    let a = -1.25;
    let mut x_wasm = [1.0, -2.0, 3.5, 4.0, -5.5];
    let mut x_scalar = x_wasm;

    wasm::scale(a, &mut x_wasm);
    scalar::scale(a, &mut x_scalar);

    for (w, s) in x_wasm.iter().zip(&x_scalar) {
        assert_eq!(w.to_bits(), s.to_bits(), "scale must match scalar bitwise");
    }
}

#[test]
fn wasm_simd128_mul_and_fma_match_scalar() {
    let a = [1.5, 2.5, 3.5, 4.5, 5.5];
    let b = [2.0, 3.0, 4.0, 5.0, 6.0];
    let c = [10.0, 20.0, 30.0, 40.0, 50.0];

    let mut out_mul_w = [0.0; 5];
    let mut out_mul_s = [0.0; 5];
    wasm::mul_elem(&a, &b, &mut out_mul_w);
    scalar::mul_elem(&a, &b, &mut out_mul_s);
    for (w, s) in out_mul_w.iter().zip(&out_mul_s) {
        assert_eq!(w.to_bits(), s.to_bits());
    }

    let mut out_fma_w = [0.0; 5];
    let mut out_fma_s = [0.0; 5];
    wasm::fma3(&a, &b, &c, &mut out_fma_w);
    scalar::fma3(&a, &b, &c, &mut out_fma_s);
    for (w, s) in out_fma_w.iter().zip(&out_fma_s) {
        assert_eq!(w.to_bits(), s.to_bits());
    }
}

#[test]
fn wasm_simd128_dot_and_sum_match_scalar() {
    let x = [1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7];
    let y = [7.7, 6.6, 5.5, 4.4, 3.3, 2.2, 1.1];

    let dot_w = wasm::dot(&x, &y);
    let dot_s = scalar::dot(&x, &y);
    assert!((dot_w - dot_s).abs() < 1e-12);

    let sum_w = wasm::sum(&x);
    let sum_s = scalar::sum(&x);
    assert!((sum_w - sum_s).abs() < 1e-12);
}
