//! fs-iga — isogeometric analysis (1D core). Layer: L3.
//!
//! Galerkin directly on B-spline spaces: the GEOMETRY BASIS IS THE ANALYSIS
//! BASIS, so the CAD→mesh information massacre never happens. This v0 is a 1D
//! B-spline Poisson solver on a clamped, uniform knot vector with k-refinement
//! (the IGA superpower — smooth high-order bases with few DOFs).
//!
//! It solves `−u''(x) = g(x)` on `[0, 1]` with homogeneous Dirichlet boundary
//! conditions by assembling the stiffness matrix `Kᵢⱼ = ∫ Nᵢ' Nⱼ'` and load
//! `fᵢ = ∫ Nᵢ g` over each knot span with Gauss–Legendre quadrature. Because a
//! degree-`p` B-spline space contains every polynomial up to degree `p`, a
//! polynomial solution is reproduced EXACTLY (to roundoff); on a smooth
//! non-polynomial solution the error falls sharply with degree.
//!
//! Deterministic; Gauss nodes use fs-math's strict cosine so raising the
//! spline degree does not introduce a platform-libm dependency.

/// A clamped, uniform B-spline space on `[0, 1]`.
#[derive(Debug, Clone, PartialEq)]
pub struct BsplineSpace {
    degree: usize,
    knots: Vec<f64>,
}

impl BsplineSpace {
    /// A clamped uniform space of the given `degree` with `elements` knot spans.
    ///
    /// # Panics
    /// If `degree == 0` or `elements == 0`.
    #[must_use]
    pub fn clamped_uniform(degree: usize, elements: usize) -> BsplineSpace {
        assert!(
            degree >= 1 && elements >= 1,
            "degree and elements must be >= 1"
        );
        // p+1 repeated 0s (clamp), interior uniform knots, p+1 repeated 1s.
        let mut knots = vec![0.0; degree + 1];
        for e in 1..elements {
            knots.push(e as f64 / elements as f64);
        }
        knots.resize(knots.len() + degree + 1, 1.0);
        BsplineSpace { degree, knots }
    }

    /// The number of basis functions (DOFs).
    #[must_use]
    pub fn num_basis(&self) -> usize {
        self.knots.len() - self.degree - 1
    }

    /// The `i`-th B-spline basis function value at `x` (Cox–de Boor).
    #[must_use]
    pub fn basis(&self, i: usize, x: f64) -> f64 {
        cox_de_boor(&self.knots, self.degree, i, clamp01(x))
    }

    /// The derivative of the `i`-th basis function at `x`.
    #[must_use]
    pub fn basis_deriv(&self, i: usize, x: f64) -> f64 {
        cox_de_boor_deriv(&self.knots, self.degree, i, clamp01(x))
    }
}

fn clamp01(x: f64) -> f64 {
    // the half-open Cox–de Boor convention gives 0 at x = 1; nudge inside.
    if x >= 1.0 { 1.0 - f64::EPSILON } else { x }
}

fn cox_de_boor(knots: &[f64], p: usize, i: usize, x: f64) -> f64 {
    if p == 0 {
        return f64::from(u8::from(knots[i] <= x && x < knots[i + 1]));
    }
    let mut val = 0.0;
    let d1 = knots[i + p] - knots[i];
    if d1 > 0.0 {
        val += (x - knots[i]) / d1 * cox_de_boor(knots, p - 1, i, x);
    }
    let d2 = knots[i + p + 1] - knots[i + 1];
    if d2 > 0.0 {
        val += (knots[i + p + 1] - x) / d2 * cox_de_boor(knots, p - 1, i + 1, x);
    }
    val
}

fn cox_de_boor_deriv(knots: &[f64], p: usize, i: usize, x: f64) -> f64 {
    if p == 0 {
        return 0.0;
    }
    let mut val = 0.0;
    let d1 = knots[i + p] - knots[i];
    if d1 > 0.0 {
        val += p as f64 / d1 * cox_de_boor(knots, p - 1, i, x);
    }
    let d2 = knots[i + p + 1] - knots[i + 1];
    if d2 > 0.0 {
        val -= p as f64 / d2 * cox_de_boor(knots, p - 1, i + 1, x);
    }
    val
}

/// The IGA solution: the B-spline coefficients + the space.
#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    space: BsplineSpace,
    coeffs: Vec<f64>,
}

impl Solution {
    /// The solution coefficients.
    #[must_use]
    pub fn coeffs(&self) -> &[f64] {
        &self.coeffs
    }

    /// Evaluate `uₕ(x) = Σ cᵢ Nᵢ(x)`.
    #[must_use]
    pub fn eval(&self, x: f64) -> f64 {
        (0..self.space.num_basis())
            .map(|i| self.coeffs[i] * self.space.basis(i, x))
            .sum()
    }

    /// The L2 error against an exact solution, by Gauss quadrature per span.
    #[must_use]
    pub fn l2_error(&self, exact: impl Fn(f64) -> f64) -> f64 {
        let mut sq = 0.0;
        for_each_quadrature(&self.space.knots, self.space.degree + 1, |x, w| {
            let e = self.eval(x) - exact(x);
            sq += w * e * e;
        });
        sq.sqrt()
    }
}

/// A structured failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IgaError {
    /// The space is too small for a Dirichlet problem (fewer than 3 DOFs).
    TooFewDofs,
    /// The assembled system was singular.
    Singular,
}

/// Deterministic Gauss-Legendre nodes and weights on [-1, 1]. The Newton
/// iteration is symmetric by construction and uses a fixed iteration ceiling.
fn gauss_legendre_rule(order: usize) -> Vec<(f64, f64)> {
    assert!(order > 0, "quadrature order must be positive");
    let mut rule = vec![(0.0, 0.0); order];
    let half = order.div_ceil(2);
    for i in 0..half {
        let angle = core::f64::consts::PI * (i as f64 + 0.75) / (order as f64 + 0.5);
        let mut root = fs_math::det::cos(angle);
        for _ in 0..64 {
            let (value, derivative) = legendre_value_derivative(order, root);
            let next = root - value / derivative;
            if (next - root).abs() <= 4.0 * f64::EPSILON * next.abs().max(1.0) {
                root = next;
                break;
            }
            root = next;
        }
        let (_, derivative) = legendre_value_derivative(order, root);
        let weight = 2.0 / ((1.0 - root * root) * derivative * derivative);
        let opposite = order - 1 - i;
        if i == opposite {
            rule[i] = (0.0, weight);
        } else {
            rule[i] = (-root, weight);
            rule[opposite] = (root, weight);
        }
    }
    rule
}

fn legendre_value_derivative(order: usize, x: f64) -> (f64, f64) {
    let mut previous = 0.0;
    let mut current = 1.0;
    for degree in 1..=order {
        let next = ((2.0 * degree as f64 - 1.0) * x * current
            - (degree - 1) as f64 * previous)
            / degree as f64;
        previous = current;
        current = next;
    }
    let derivative = order as f64 * (x * current - previous) / (x * x - 1.0);
    (current, derivative)
}

fn for_each_quadrature(knots: &[f64], order: usize, mut f: impl FnMut(f64, f64)) {
    let rule = gauss_legendre_rule(order);
    for w in knots.windows(2) {
        let (a, b) = (w[0], w[1]);
        if b <= a {
            continue;
        }
        let (mid, half) = (f64::midpoint(a, b), (b - a) / 2.0);
        for &(node, weight) in &rule {
            let x = mid + half * node;
            f(x, weight * half);
        }
    }
}

/// Solve `−u'' = g` on `[0, 1]` with `u(0) = u(1) = 0` on the B-spline space.
///
/// # Errors
/// [`IgaError::TooFewDofs`] if the space is too small; [`IgaError::Singular`]
/// if the reduced system is singular.
pub fn solve_poisson(space: &BsplineSpace, g: impl Fn(f64) -> f64) -> Result<Solution, IgaError> {
    let n = space.num_basis();
    if n < 3 {
        return Err(IgaError::TooFewDofs);
    }
    // full stiffness + load.
    let mut k = vec![vec![0.0; n]; n];
    let mut load = vec![0.0; n];
    // p+1 Gauss points integrate stiffness products (degree 2p-2) exactly
    // and avoid the rank ceiling caused by a fixed quadrature rule.
    for_each_quadrature(&space.knots, space.degree + 1, |x, w| {
        let dphi: Vec<f64> = (0..n).map(|i| space.basis_deriv(i, x)).collect();
        let phi: Vec<f64> = (0..n).map(|i| space.basis(i, x)).collect();
        let gx = g(x);
        for i in 0..n {
            load[i] += w * phi[i] * gx;
            for j in 0..n {
                k[i][j] += w * dphi[i] * dphi[j];
            }
        }
    });
    // homogeneous Dirichlet: clamp basis 0 and n-1 (interpolatory at the ends).
    // solve the interior block 1..n-1.
    let interior: Vec<usize> = (1..n - 1).collect();
    let m = interior.len();
    let mut a = vec![vec![0.0; m]; m];
    let mut b = vec![0.0; m];
    for (ri, &i) in interior.iter().enumerate() {
        b[ri] = load[i];
        for (cj, &j) in interior.iter().enumerate() {
            a[ri][cj] = k[i][j];
        }
    }
    let sol = gauss_solve(a, b).ok_or(IgaError::Singular)?;
    let mut coeffs = vec![0.0; n];
    for (ri, &i) in interior.iter().enumerate() {
        coeffs[i] = sol[ri];
    }
    Ok(Solution {
        space: space.clone(),
        coeffs,
    })
}

/// Gaussian elimination with partial pivoting.
fn gauss_solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        let piv = (col..n).max_by(|&r1, &r2| {
            a[r1][col]
                .abs()
                .partial_cmp(&a[r2][col].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;
        if a[piv][col].abs() <= 1e-14 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        let pivot = a[col].clone();
        let bcol = b[col];
        for r in (col + 1)..n {
            let f = a[r][col] / pivot[col];
            for (arc, pc) in a[r].iter_mut().zip(&pivot).skip(col) {
                *arc -= f * pc;
            }
            b[r] -= f * bcol;
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let s: f64 = a[i][(i + 1)..n]
            .iter()
            .zip(&x[(i + 1)..n])
            .map(|(aij, xj)| aij * xj)
            .sum();
        x[i] = (b[i] - s) / a[i][i];
    }
    Some(x)
}
