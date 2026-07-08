//! fs-symmetry — symmetry harvesting (plan addendum, Proposal 13). Layer: L1.
//!
//! Declared (or detected) symmetry is harvested as BOTH correctness AND speed:
//! a `k`-fold cyclic symmetry BLOCK-DIAGONALIZES a solve by irreducible
//! representation, turning one `n×n` solve into `n` independent scalar solves —
//! the win condition "`k`-fold symmetry ≈ `k`× speedup". Concretely, a
//! circulant operator (the algebraic signature of cyclic symmetry) is
//! diagonalized by the DISCRETE FOURIER TRANSFORM (the isotypic projection for
//! the cyclic group), so [`solve_circulant`] costs `O(n²)` (naive DFT) instead
//! of `O(n³)`.
//!
//! Symmetry is rarely exact, so this is not binary. [`cyclic_residual`]
//! measures how far a field is from `k`-fold symmetric and CERTIFIES the
//! asymmetry residual; [`symmetrize`] splits a field into its symmetric part
//! and an asymmetric remainder; and [`symmetrized_solve`] solves the symmetric
//! part cheaply while returning a CERTIFIED PERTURBATION BOUND that provably
//! contains the correction due to the asymmetry (the remainder enters as a
//! right-hand-side perturbation with an interval bound).
//!
//! v1 covers the CYCLIC group `Cₙ` (turbomachinery, fasteners, most
//! architectural regularity); dihedral + geometric-hashing detection are noted
//! as follow-ons. Everything is deterministic; no dependencies.

use std::f64::consts::TAU;

/// A minimal complex number (no external complex crate needed).
#[derive(Debug, Clone, Copy, PartialEq)]
struct Cx {
    re: f64,
    im: f64,
}

impl Cx {
    const ZERO: Cx = Cx { re: 0.0, im: 0.0 };
    fn new(re: f64, im: f64) -> Cx {
        Cx { re, im }
    }
    fn from_angle(theta: f64) -> Cx {
        Cx {
            re: theta.cos(),
            im: theta.sin(),
        }
    }
    fn add(self, o: Cx) -> Cx {
        Cx::new(self.re + o.re, self.im + o.im)
    }
    fn mul(self, o: Cx) -> Cx {
        Cx::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
    fn div(self, o: Cx) -> Cx {
        let d = o.re * o.re + o.im * o.im;
        Cx::new(
            (self.re * o.re + self.im * o.im) / d,
            (self.im * o.re - self.re * o.im) / d,
        )
    }
    fn modulus(self) -> f64 {
        self.re.hypot(self.im)
    }
}

/// A structured symmetry failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SymmetryError {
    /// An input vector is empty.
    EmptyInput,
    /// Two vectors have mismatched lengths.
    LengthMismatch {
        /// Expected length.
        expected: usize,
        /// Actual length.
        found: usize,
    },
    /// The fold count is zero.
    ZeroFold,
    /// The vector length is not divisible by the fold count (no clean sectors).
    NotDivisible {
        /// The vector length.
        len: usize,
        /// The requested fold.
        k_fold: usize,
    },
    /// The circulant operator is singular (a zero eigenvalue).
    Singular,
}

/// The cyclic group `Cₙ`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CyclicGroup {
    /// The order.
    pub n: usize,
}

impl CyclicGroup {
    /// `Cₙ`.
    #[must_use]
    pub fn new(n: usize) -> CyclicGroup {
        CyclicGroup { n }
    }

    /// The character of irreducible representation `irrep` evaluated at group
    /// element `element`: `e^{2πi · irrep · element / n}` as `(re, im)`. The
    /// `Cₙ` character table is exactly the roots of unity.
    #[must_use]
    pub fn character(&self, irrep: usize, element: usize) -> (f64, f64) {
        if self.n == 0 {
            return (1.0, 0.0);
        }
        let theta = TAU * ((irrep * element) % self.n) as f64 / self.n as f64;
        (theta.cos(), theta.sin())
    }
}

/// Naive DFT (`O(n²)`); `sign = -1.0` forward, `+1.0` inverse (unnormalized).
fn dft(x: &[Cx], sign: f64) -> Vec<Cx> {
    let n = x.len();
    (0..n)
        .map(|k| {
            let mut acc = Cx::ZERO;
            for (j, &xj) in x.iter().enumerate() {
                let theta = sign * TAU * (k * j) as f64 / n as f64;
                acc = acc.add(xj.mul(Cx::from_angle(theta)));
            }
            acc
        })
        .collect()
}

/// The eigenvalues of the circulant `C[i][j] = first_row[(j-i) mod n]`. Its
/// eigenvector `(v_k)_j = e^{2πi jk/n}` has eigenvalue `Σ_m c_m e^{2πi mk/n}`,
/// i.e. the inverse-sign DFT of the first row.
fn circulant_eigenvalues(first_row: &[f64]) -> Vec<Cx> {
    let c: Vec<Cx> = first_row.iter().map(|&v| Cx::new(v, 0.0)).collect();
    dft(&c, 1.0)
}

/// Multiply the circulant matrix defined by `first_row` (row 0; row `i` is
/// `first_row` rotated right by `i`) by `x`.
///
/// # Errors
/// [`SymmetryError`] on empty or length-mismatched input.
pub fn circulant_matvec(first_row: &[f64], x: &[f64]) -> Result<Vec<f64>, SymmetryError> {
    let n = first_row.len();
    if n == 0 {
        return Err(SymmetryError::EmptyInput);
    }
    if x.len() != n {
        return Err(SymmetryError::LengthMismatch {
            expected: n,
            found: x.len(),
        });
    }
    // C[i][j] = first_row[(j - i) mod n].
    let y = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| first_row[(j + n - i % n) % n] * x[j])
                .sum::<f64>()
        })
        .collect();
    Ok(y)
}

/// Solve `C x = rhs` for the circulant `C` (given by `first_row`) by isotypic
/// (DFT) block-diagonalization: `x = IDFT(DFT(rhs) / eigenvalues)`. This is the
/// `O(n²)` symmetry-adapted solve; the result solves the full system.
///
/// # Errors
/// [`SymmetryError`] on empty / mismatched input or a singular operator.
pub fn solve_circulant(first_row: &[f64], rhs: &[f64]) -> Result<Vec<f64>, SymmetryError> {
    let n = first_row.len();
    if n == 0 {
        return Err(SymmetryError::EmptyInput);
    }
    if rhs.len() != n {
        return Err(SymmetryError::LengthMismatch {
            expected: n,
            found: rhs.len(),
        });
    }
    let lambda = circulant_eigenvalues(first_row);
    if lambda.iter().any(|l| l.modulus() <= 1e-12) {
        return Err(SymmetryError::Singular);
    }
    let b: Vec<Cx> = rhs.iter().map(|&v| Cx::new(v, 0.0)).collect();
    let bhat = dft(&b, -1.0);
    let xhat: Vec<Cx> = bhat
        .iter()
        .zip(&lambda)
        .map(|(bk, lk)| bk.div(*lk))
        .collect();
    let x = dft(&xhat, 1.0);
    // inverse DFT normalization + take the real part.
    Ok(x.iter().map(|c| c.re / n as f64).collect())
}

/// The certified asymmetry residual of a field under `k`-fold cyclic symmetry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SymmetryResidual {
    /// `||v - rotate(v, len/k)||₂` — the un-normalized residual.
    pub residual: f64,
    /// The residual relative to `||v||₂` (0 if `v` is zero).
    pub relative: f64,
    /// Is the field exactly `k`-fold symmetric (within tolerance)?
    pub is_exact: bool,
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// Certify how far `v` is from `k`-fold cyclic symmetry (invariance under a
/// rotation by `len / k_fold` samples).
///
/// # Errors
/// [`SymmetryError`] on empty input, zero fold, or a length not divisible by
/// `k_fold`.
pub fn cyclic_residual(v: &[f64], k_fold: usize) -> Result<SymmetryResidual, SymmetryError> {
    let n = v.len();
    if n == 0 {
        return Err(SymmetryError::EmptyInput);
    }
    if k_fold == 0 {
        return Err(SymmetryError::ZeroFold);
    }
    if !n.is_multiple_of(k_fold) {
        return Err(SymmetryError::NotDivisible { len: n, k_fold });
    }
    let shift = n / k_fold;
    let diff: Vec<f64> = (0..n).map(|i| v[i] - v[(i + shift) % n]).collect();
    let residual = norm(&diff);
    let base = norm(v);
    let relative = if base <= f64::EPSILON {
        0.0
    } else {
        residual / base
    };
    Ok(SymmetryResidual {
        residual,
        relative,
        is_exact: residual <= 1e-12,
    })
}

/// Split `v` into its `k`-fold-symmetric part (the average over the `k`
/// rotations — the projection onto the trivial isotypic component) and the
/// asymmetric remainder (`v - symmetric`).
///
/// # Errors
/// [`SymmetryError`] on empty input, zero fold, or an indivisible length.
pub fn symmetrize(v: &[f64], k_fold: usize) -> Result<(Vec<f64>, Vec<f64>), SymmetryError> {
    let n = v.len();
    if n == 0 {
        return Err(SymmetryError::EmptyInput);
    }
    if k_fold == 0 {
        return Err(SymmetryError::ZeroFold);
    }
    if !n.is_multiple_of(k_fold) {
        return Err(SymmetryError::NotDivisible { len: n, k_fold });
    }
    let shift = n / k_fold;
    let sym: Vec<f64> = (0..n)
        .map(|i| {
            let s: f64 = (0..k_fold).map(|j| v[(i + j * shift) % n]).sum();
            s / k_fold as f64
        })
        .collect();
    let asym: Vec<f64> = (0..n).map(|i| v[i] - sym[i]).collect();
    Ok((sym, asym))
}

/// A symmetrized solve with a certified perturbation bound.
#[derive(Debug, Clone, PartialEq)]
pub struct PerturbationBound {
    /// The solution of the symmetric part of the right-hand side.
    pub symmetric_solution: Vec<f64>,
    /// `||asymmetric remainder||₂` — the certified asymmetry residual.
    pub asymmetry_residual: f64,
    /// A certified bound on `||x_full - symmetric_solution||₂`, namely
    /// `asymmetry_residual / λ_min` (`λ_min` = smallest `|eigenvalue|`).
    pub correction_bound: f64,
}

/// Solve the `k`-fold-symmetric part of `rhs` against the circulant `first_row`
/// and certify a bound on the correction the asymmetric remainder would add.
/// The bound provably CONTAINS `||x_full − symmetric_solution||`.
///
/// # Errors
/// [`SymmetryError`] on bad input or a singular operator.
pub fn symmetrized_solve(
    first_row: &[f64],
    rhs: &[f64],
    k_fold: usize,
) -> Result<PerturbationBound, SymmetryError> {
    if rhs.len() != first_row.len() {
        return Err(SymmetryError::LengthMismatch {
            expected: first_row.len(),
            found: rhs.len(),
        });
    }
    let (b_sym, b_asym) = symmetrize(rhs, k_fold)?;
    let symmetric_solution = solve_circulant(first_row, &b_sym)?;
    let lambda_min = circulant_eigenvalues(first_row)
        .iter()
        .map(|l| l.modulus())
        .fold(f64::INFINITY, f64::min);
    let asymmetry_residual = norm(&b_asym);
    Ok(PerturbationBound {
        symmetric_solution,
        asymmetry_residual,
        correction_bound: asymmetry_residual / lambda_min,
    })
}
