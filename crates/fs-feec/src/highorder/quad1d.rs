//! 1D Gauss–Legendre quadrature and Legendre/Lobatto polynomial
//! evaluation — the scalar bedrock of the tensor-product high-order
//! families (tfz.6). Deterministic: nodes come from Newton iteration
//! with fixed starting guesses and a fixed iteration count, all
//! arithmetic through fs-math strict kernels.

use fs_math::det;

/// Legendre polynomial L_n(x) and its derivative, by the three-term
/// recurrence (numerically stable on [-1, 1]).
#[must_use]
pub fn legendre(n: usize, x: f64) -> (f64, f64) {
    if n == 0 {
        return (1.0, 0.0);
    }
    let (mut pm1, mut p) = (1.0f64, x);
    for k in 1..n {
        let kf = k as f64;
        let next = ((2.0 * kf + 1.0) * x * p - kf * pm1) / (kf + 1.0);
        pm1 = p;
        p = next;
    }
    // L'_n from the standard identity (guarded at |x| = 1 where the
    // denominator vanishes; endpoint derivative is n(n+1)/2 · (±1)ⁿ⁺¹).
    let denom = x.mul_add(-x, 1.0);
    let dp = if denom.abs() < 1e-14 {
        let nf = n as f64;
        let sign = if x > 0.0 || (n + 1).is_multiple_of(2) {
            1.0
        } else {
            -1.0
        };
        sign * nf * (nf + 1.0) / 2.0
    } else {
        (n as f64) * x.mul_add(p, -pm1) / -denom
    };
    (p, dp)
}

/// Gauss–Legendre nodes and weights on [-1, 1]: `n` points integrate
/// polynomials of degree ≤ 2n−1 exactly. Chebyshev starting guesses +
/// 25 Newton steps (converged far past f64 by then; FIXED count keeps
/// the bit pattern independent of convergence-test vagaries).
#[must_use]
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    assert!(n >= 1, "quadrature needs at least one point");
    let mut nodes = vec![0.0f64; n];
    let mut weights = vec![0.0f64; n];
    for i in 0..n {
        // Chebyshev-angle initial guess for root i (descending order).
        let theta = std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5);
        let mut x = det::cos(theta);
        for _ in 0..25 {
            let (p, dp) = legendre(n, x);
            x -= p / dp;
        }
        let (_, dp) = legendre(n, x);
        nodes[i] = x;
        weights[i] = 2.0 / (x.mul_add(-x, 1.0) * dp * dp);
    }
    // Ascending order (Newton from descending guesses gives descending
    // roots); fixed deterministic reorder.
    nodes.reverse();
    weights.reverse();
    (nodes, weights)
}

/// Lobatto hierarchical shape functions on [-1, 1] and derivatives:
/// N_0 = (1−x)/2, N_1 = (1+x)/2 (the vertex pair), and for k ≥ 2 the
/// integrated-Legendre bubbles
/// N_k = (L_k(x) − L_{k−2}(x)) / √(2(2k−1)), which vanish at ±1.
/// Returns (values, derivatives), length `order + 1`.
///
/// # Panics
/// If `order < 1` or the basis extent cannot be represented by `usize`.
#[must_use]
pub fn lobatto_shapes(order: usize, x: f64) -> (Vec<f64>, Vec<f64>) {
    assert!(order >= 1, "continuous H1 basis needs order >= 1");
    let basis_len = order.checked_add(1).expect("Lobatto basis extent overflow");
    let mut vals = Vec::with_capacity(basis_len);
    let mut ders = Vec::with_capacity(basis_len);
    vals.push(f64::midpoint(1.0, -x));
    ders.push(-0.5);
    vals.push(f64::midpoint(1.0, x));
    ders.push(0.5);
    for k in 2..=order {
        let (lk, dlk) = legendre(k, x);
        let (lk2, dlk2) = legendre(k - 2, x);
        let scale = 1.0 / det::sqrt(2.0 * (2.0 * k as f64 - 1.0));
        vals.push((lk - lk2) * scale);
        ders.push((dlk - dlk2) * scale);
    }
    (vals, ders)
}

/// 1D element mass and stiffness matrices for the order-r Lobatto
/// basis on a physical cell of width `h` (reference [-1, 1] scaled by
/// the affine map): M_e = (h/2)·∫N_i N_j, K_e = (2/h)·∫N'_i N'_j.
/// Quadrature with r+2 points (exact: integrands have degree ≤ 2r).
///
/// # Panics
/// If the basis, quadrature, or dense matrix extent cannot be represented by
/// `usize`, or if `order < 1` (as required by the continuous H¹ basis).
#[must_use]
pub fn element_matrices(order: usize, h: f64) -> (Vec<f64>, Vec<f64>) {
    assert!(order >= 1, "element matrices need order >= 1");
    let n = order.checked_add(1).expect("element basis extent overflow");
    let nq = order
        .checked_add(2)
        .expect("element quadrature extent overflow");
    let matrix_len = n.checked_mul(n).expect("element matrix extent overflow");
    let (qx, qw) = gauss_legendre(nq);
    let mut mass = vec![0.0f64; matrix_len];
    let mut stiff = vec![0.0f64; matrix_len];
    for (&x, &w) in qx.iter().zip(&qw) {
        let (vals, ders) = lobatto_shapes(order, x);
        for i in 0..n {
            for j in 0..n {
                mass[i * n + j] = (w * vals[i]).mul_add(vals[j], mass[i * n + j]);
                stiff[i * n + j] = (w * ders[i]).mul_add(ders[j], stiff[i * n + j]);
            }
        }
    }
    for v in &mut mass {
        *v *= h / 2.0;
    }
    for v in &mut stiff {
        *v *= 2.0 / h;
    }
    (mass, stiff)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn panic_text(payload: &(dyn std::any::Any + Send)) -> &str {
        payload
            .downcast_ref::<&str>()
            .copied()
            .or_else(|| payload.downcast_ref::<String>().map(String::as_str))
            .unwrap_or("non-string panic")
    }

    #[test]
    fn polynomial_extents_fail_before_wrapping_or_allocating() {
        let basis = std::panic::catch_unwind(|| lobatto_shapes(usize::MAX, 0.0))
            .expect_err("basis extent must not wrap");
        assert!(panic_text(basis.as_ref()).contains("basis extent overflow"));

        let element_basis = std::panic::catch_unwind(|| element_matrices(usize::MAX, 1.0))
            .expect_err("element basis extent must not wrap");
        assert!(panic_text(element_basis.as_ref()).contains("element basis extent overflow"));

        let matrix = std::panic::catch_unwind(|| element_matrices(usize::MAX / 2, 1.0))
            .expect_err("element matrix extent must not wrap");
        assert!(panic_text(matrix.as_ref()).contains("matrix extent overflow"));

        let quadrature = std::panic::catch_unwind(|| element_matrices(usize::MAX - 1, 1.0))
            .expect_err("element quadrature extent must not wrap");
        assert!(panic_text(quadrature.as_ref()).contains("quadrature extent overflow"));
    }
}
