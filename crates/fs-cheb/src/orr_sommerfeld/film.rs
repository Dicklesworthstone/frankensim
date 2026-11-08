//! Free-surface temporal stability of a Newtonian Nusselt film.
//!
//! Length and velocity scales are thickness h and mean speed Ubar. With
//! y in [0,1], U=3y-3y²/2, alpha=kh, Re=Ubar*h/nu, G=g_n*h/Ubar² and
//! S=sigma/(rho*Ubar²*h), disturbances are exp(lambda*t + i*alpha*x).
//! The original equations are
//! (D²-alpha²)² phi/Re - i*alpha*[U*(D²-alpha²)phi + 3phi]
//!     = lambda*(D²-alpha²)phi.
//! At the wall phi=Dphi=0. At the free surface:
//! (D²+alpha²)phi=3eta;
//! lambda*eta=-i*alpha*(phi + 3eta/2);
//! lambda*Dphi=(D³-3alpha²D)phi/Re - i*alpha*3Dphi/2
//!     - i*alpha*(G+S*alpha²)*eta.
//!
//! These are the Newtonian limit of Chimetta & Franklin (2023), equations
//! (20), (22)-(29), with the normal coordinate directed FROM the wall and
//! scales changed to h and Ubar: DOI 10.1016/j.jnnfm.2023.105153,
//! <https://arxiv.org/pdf/2311.11986>.
//!
//! The coefficient basis phi=y²*p enforces both wall conditions identically.
//! The tangential condition eliminates eta exactly; it does not impose a rigid
//! lid. Native Chebyshev differentiation builds this basis. An equilibrated
//! shift-invert pencil and the existing complex LU/QR retain the surface and
//! shear modes. A second order and native SVD original-pencil residuals admit
//! each requested mode. These are numerical checks, not certified enclosures.

use crate::{Cheb1, lobatto_points};
use fs_la::eigen_complex::{
    ComplexSvdError, EigError, LuComplexError, eig_with_cancel, lu_complex_with_cancel,
    svd_complex_with_cancel,
};
use fs_math::c64::C64;

/// Explicitly dense film lane; p has at most this Chebyshev degree.
pub const MAX_FILM_ORDER: usize = 64;

/// Dimensionless data for the same shear-free Nusselt base profile.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilmParameters {
    /// Mean-speed Reynolds number Ubar*h/nu = q/nu, strictly positive.
    pub reynolds: f64,
    /// Positive streamwise angular wavenumber times thickness, k*h.
    pub kh: f64,
    /// Wall-directed normal acceleration g_n*h/Ubar², nonnegative.
    pub normal_gravity: f64,
    /// Surface tension sigma/(rho*Ubar²*h), nonnegative.
    pub capillarity: f64,
}

/// Bounded spectral resolution and explicit acceptance tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilmSpectrumBudget {
    /// Degree of p in phi=y²*p, at least eight.
    pub coarse_order: usize,
    /// At least coarse_order+4 and at most MAX_FILM_ORDER.
    pub fine_order: usize,
    /// Number of rightmost modes; at most coarse_order+1.
    pub modes: usize,
    /// Native QR sweep cap per active eigenvalue window, in 1..=4096.
    pub max_qr_sweeps: usize,
    /// Native Jacobi SVD sweep cap per residual query, in 1..=4096.
    pub max_svd_sweeps: usize,
    /// Maximum componentwise original-pencil backward residual, in (0,1).
    pub relative_residual_tolerance: f64,
    /// Absolute tolerance on the complex eigenvalue change between orders.
    pub absolute_resolution_tolerance: f64,
    /// Relative tolerance on that change, using the larger eigenvalue modulus.
    pub relative_resolution_tolerance: f64,
}

/// A temporal mode admitted at both authored collocation orders.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilmEigenmode {
    /// lambda on time scale h/Ubar; positive real part grows.
    pub eigenvalue: C64,
    /// Complex-plane distance to the matched mode at the coarser order.
    /// A resolution indicator, not a rigorous eigenvalue error bound.
    pub order_difference: f64,
    /// Largest original-pencil backward residual at the two orders.
    pub relative_residual: f64,
}

/// No partial spectrum is returned on numerical refusal or cancellation.
#[derive(Clone, Debug, PartialEq)]
pub enum FilmSpectrumError {
    Invalid(&'static str),
    Cancelled,
    Singular {
        index: usize,
    },
    Eigensolve {
        window_hi: usize,
    },
    SingularVector(ComplexSvdError),
    Residual {
        relative: f64,
        tolerance: f64,
    },
    Resolution {
        mode: usize,
        difference: f64,
        tolerance: f64,
    },
}

impl core::fmt::Display for FilmSpectrumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "free-surface film spectrum refused: {self:?}")
    }
}
impl std::error::Error for FilmSpectrumError {}

/// Reusable native coefficient/derivative grids at both authored orders.
///
/// Preparation is O(N³) at these bounded degrees. Each solve retains the full
/// discretized pencil, with bounded LU/QR work and one residual SVD per requested
/// mode and order. No parameter fitting, mode removal or long-wave expansion
/// enters the operator. Surface growth is not a turbulent-transition criterion.
pub struct PreparedFilmSpectrum {
    budget: FilmSpectrumBudget,
    coarse: Grid,
    fine: Grid,
}

impl PreparedFilmSpectrum {
    pub fn new(
        budget: FilmSpectrumBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, FilmSpectrumError> {
        poll(&mut cancelled)?;
        if budget.coarse_order < 8
            || budget.coarse_order > MAX_FILM_ORDER - 4
            || budget.fine_order < budget.coarse_order + 4
            || budget.fine_order > MAX_FILM_ORDER
            || budget.modes == 0
            || budget.modes > budget.coarse_order + 1
            || !(1..=4096).contains(&budget.max_qr_sweeps)
            || !(1..=4096).contains(&budget.max_svd_sweeps)
            || !positive(budget.relative_residual_tolerance)
            || budget.relative_residual_tolerance >= 1.0
            || !positive(budget.absolute_resolution_tolerance)
            || !budget.relative_resolution_tolerance.is_finite()
            || budget.relative_resolution_tolerance < 0.0
            || budget.relative_resolution_tolerance >= 1.0
        {
            return Err(FilmSpectrumError::Invalid("invalid film spectral budget"));
        }
        let coarse = Grid::new(budget.coarse_order, &mut cancelled)?;
        let fine = Grid::new(budget.fine_order, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(Self {
            budget,
            coarse,
            fine,
        })
    }

    /// Rightmost modes, descending growth with deterministic imaginary tie-break.
    /// Both orders' rightmost sets must match one-to-one within the authored
    /// tolerance. Every returned mode passes the original-equation residual.
    pub fn solve(
        &self,
        parameters: FilmParameters,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<FilmEigenmode>, FilmSpectrumError> {
        poll(&mut cancelled)?;
        if !positive(parameters.reynolds)
            || !positive(parameters.kh)
            || !parameters.normal_gravity.is_finite()
            || parameters.normal_gravity < 0.0
            || !parameters.capillarity.is_finite()
            || parameters.capillarity < 0.0
        {
            return Err(FilmSpectrumError::Invalid(
                "finite positive Re/kh and nonnegative gravity/capillarity required",
            ));
        }
        let coarse = self.coarse.pencil(parameters, &mut cancelled)?;
        let fine = self.fine.pencil(parameters, &mut cancelled)?;
        let coarse_values = coarse.eigenvalues(self.budget, &mut cancelled)?;
        let fine_values = fine.eigenvalues(self.budget, &mut cancelled)?;
        let mut used = vec![false; self.budget.modes];
        let mut modes = Vec::with_capacity(self.budget.modes);
        for (mode, &value) in fine_values.iter().take(self.budget.modes).enumerate() {
            poll(&mut cancelled)?;
            let index = (0..self.budget.modes)
                .filter(|&index| !used[index])
                .min_by(|&a, &b| {
                    (value - coarse_values[a])
                        .abs()
                        .total_cmp(&(value - coarse_values[b]).abs())
                        .then(a.cmp(&b))
                })
                .ok_or(FilmSpectrumError::Invalid(
                    "spectral mode matching exhausted",
                ))?;
            used[index] = true;
            let difference = (value - coarse_values[index]).abs();
            let tolerance = self.budget.absolute_resolution_tolerance
                + self.budget.relative_resolution_tolerance
                    * value.abs().max(coarse_values[index].abs());
            if !difference.is_finite() || !tolerance.is_finite() || difference > tolerance {
                return Err(FilmSpectrumError::Resolution {
                    mode,
                    difference,
                    tolerance,
                });
            }
            let residual = fine
                .residual(value, self.budget, &mut cancelled)?
                .max(coarse.residual(coarse_values[index], self.budget, &mut cancelled)?);
            modes.push(FilmEigenmode {
                eigenvalue: value,
                order_difference: difference,
                relative_residual: residual,
            });
        }
        poll(&mut cancelled)?;
        Ok(modes)
    }
}

struct Grid {
    nodes: Vec<f64>,
    // Derivatives D^j[y²*T_k(2y-1)], j=0..4, in native Chebyshev coefficient order.
    derivatives: [Vec<f64>; 5],
}

impl Grid {
    fn new(order: usize, cancelled: &mut impl FnMut() -> bool) -> Result<Self, FilmSpectrumError> {
        let n = order + 1;
        let nodes: Vec<_> = lobatto_points(order)
            .into_iter()
            .map(|x| 0.5 * (x + 1.0))
            .collect();
        let mut derivatives: [Vec<f64>; 5] = core::array::from_fn(|_| vec![0.0; n * n]);
        for column in 0..n {
            poll(cancelled)?;
            let mut coefficients = vec![0.0; column + 1];
            coefficients[column] = if column == 0 { 2.0 } else { 1.0 };
            let mut polynomial = Cheb1::from_coeffs(0.0, 1.0, coefficients);
            for values in &mut derivatives {
                for (row, &y) in nodes.iter().enumerate() {
                    values[row * n + column] = polynomial.eval(y);
                }
                polynomial = polynomial.differentiate();
            }
        }
        for (row, &y) in nodes.iter().enumerate() {
            poll(cancelled)?;
            for column in 0..n {
                let index = row * n + column;
                let [p, dp, d2p, d3p, d4p] = core::array::from_fn(|j| derivatives[j][index]);
                derivatives[0][index] = y * y * p;
                derivatives[1][index] = y * y * dp + 2.0 * y * p;
                derivatives[2][index] = y * y * d2p + 4.0 * y * dp + 2.0 * p;
                derivatives[3][index] = y * y * d3p + 6.0 * y * d2p + 6.0 * dp;
                derivatives[4][index] = y * y * d4p + 8.0 * y * d3p + 12.0 * d2p;
            }
        }
        Ok(Self { nodes, derivatives })
    }

    fn pencil(
        &self,
        data: FilmParameters,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Pencil, FilmSpectrumError> {
        let n = self.nodes.len();
        let mut a = vec![C64::ZERO; n * n];
        let mut b = vec![C64::ZERO; n * n];
        let [d0, d1, d2, d3, d4] = &self.derivatives;
        let alpha = data.kh;
        let alpha2 = alpha * alpha;
        let alpha4 = alpha2 * alpha2;
        let restoring = data.normal_gravity + data.capillarity * alpha2;
        for i in 1..(n - 1) {
            poll(cancelled)?;
            let y = self.nodes[i];
            let u = 3.0 * y - 1.5 * y * y;
            for j in 0..n {
                let source = i * n + j;
                let target = (i - 1) * n + j;
                let laplace = d2[source] - alpha2 * d0[source];
                a[target] = C64::new(
                    (d4[source] - 2.0 * alpha2 * d2[source] + alpha4 * d0[source]) / data.reynolds,
                    -alpha * (u * laplace + 3.0 * d0[source]),
                );
                b[target] = C64::from_re(laplace);
            }
        }
        for j in 0..n {
            let eta = (d2[j] + alpha2 * d0[j]) / 3.0;
            let normal = (n - 2) * n + j;
            a[normal] = C64::new(
                (d3[j] - 3.0 * alpha2 * d1[j]) / data.reynolds,
                -alpha * (1.5 * d1[j] + restoring * eta),
            );
            b[normal] = C64::from_re(d1[j]);
            let kinematic = (n - 1) * n + j;
            a[kinematic] = C64::new(0.0, -alpha * (d0[j] + 1.5 * eta));
            b[kinematic] = C64::from_re(eta);
        }
        if a.iter().chain(&b).any(|&z| !finite(z)) {
            return Err(FilmSpectrumError::Invalid(
                "film spectral pencil is not representable",
            ));
        }
        // Equivalent row/coordinate scaling, shared by the eigen solve AND its
        // original-pencil residual. No physical coefficient is changed.
        for i in 0..n {
            poll(cancelled)?;
            let scale = (0..n)
                .map(|j| a[i * n + j].abs().max(b[i * n + j].abs()))
                .fold(0.0, f64::max);
            if !positive(scale) {
                return Err(FilmSpectrumError::Invalid(
                    "zero or unrepresentable film pencil row",
                ));
            }
            for j in 0..n {
                a[i * n + j] = divide(a[i * n + j], scale);
                b[i * n + j] = divide(b[i * n + j], scale);
            }
        }
        for j in 0..n {
            poll(cancelled)?;
            let scale = (0..n)
                .map(|i| a[i * n + j].abs().max(b[i * n + j].abs()))
                .fold(0.0, f64::max);
            if !positive(scale) {
                return Err(FilmSpectrumError::Invalid(
                    "zero or unrepresentable film pencil column",
                ));
            }
            for i in 0..n {
                a[i * n + j] = divide(a[i * n + j], scale);
                b[i * n + j] = divide(b[i * n + j], scale);
            }
        }
        Ok(Pencil { a, b, n })
    }
}

struct Pencil {
    a: Vec<C64>,
    b: Vec<C64>,
    n: usize,
}

impl Pencil {
    fn eigenvalues(
        &self,
        budget: FilmSpectrumBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, FilmSpectrumError> {
        // Algebraic shift-invert: (A-B)^-1 B has mu=1/(lambda-1).
        // Direct B^-1 A magnifies the fourth-derivative dynamic range and loses
        // near-neutral growth. The shift is not a physical damping term.
        let shifted: Vec<_> = self.a.iter().zip(&self.b).map(|(&a, &b)| a - b).collect();
        let factor = lu_complex_with_cancel(&shifted, self.n, &mut *cancelled).map_err(lu_error)?;
        let mut matrix = vec![C64::ZERO; self.n * self.n];
        let mut column = vec![C64::ZERO; self.n];
        for j in 0..self.n {
            poll(cancelled)?;
            for (i, value) in column.iter_mut().enumerate() {
                *value = self.b[i * self.n + j];
            }
            factor
                .solve_with_cancel(&mut column, &mut *cancelled)
                .map_err(lu_error)?;
            for (i, &value) in column.iter().enumerate() {
                matrix[i * self.n + j] = value;
            }
        }
        let values = eig_with_cancel(&matrix, self.n, budget.max_qr_sweeps, &mut *cancelled)
            .map_err(|error| match error {
                EigError::Cancelled => FilmSpectrumError::Cancelled,
                EigError::Invalid { what } => FilmSpectrumError::Invalid(what),
                EigError::NotConverged(failure) => FilmSpectrumError::Eigensolve {
                    window_hi: failure.window_hi,
                },
            })?;
        let mut values: Vec<_> = values
            .into_iter()
            .map(|mu| C64::ONE + C64::ONE / mu)
            .collect();
        if values.iter().any(|&value| !finite(value)) {
            return Err(FilmSpectrumError::Invalid(
                "unresolved infinite film eigenvalue",
            ));
        }
        values.sort_by(|a, b| b.re.total_cmp(&a.re).then_with(|| b.im.total_cmp(&a.im)));
        Ok(values)
    }

    fn residual(
        &self,
        lambda: C64,
        budget: FilmSpectrumBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, FilmSpectrumError> {
        let mut residual_matrix = vec![C64::ZERO; self.n * self.n];
        for i in 0..self.n {
            poll(cancelled)?;
            let scale: f64 = (0..self.n)
                .map(|j| self.a[i * self.n + j].abs() + lambda.abs() * self.b[i * self.n + j].abs())
                .sum();
            if !positive(scale) {
                return Err(FilmSpectrumError::Invalid(
                    "film residual row is not representable",
                ));
            }
            for j in 0..self.n {
                let index = i * self.n + j;
                residual_matrix[index] = divide(self.a[index] - lambda * self.b[index], scale);
            }
        }
        let svd = svd_complex_with_cancel(
            &residual_matrix,
            self.n,
            self.n,
            64.0 * f64::EPSILON,
            budget.max_svd_sweeps,
            &mut *cancelled,
        )
        .map_err(|error| match error {
            ComplexSvdError::Cancelled => FilmSpectrumError::Cancelled,
            other => FilmSpectrumError::SingularVector(other),
        })?;
        let mut maximum = 0.0_f64;
        for i in 0..self.n {
            poll(cancelled)?;
            let mut value = C64::ZERO;
            let mut scale = 0.0;
            for j in 0..self.n {
                let index = i * self.n + j;
                let vector = svd.v[j * self.n + self.n - 1];
                value = value + (self.a[index] - lambda * self.b[index]) * vector;
                scale += (self.a[index].abs() + lambda.abs() * self.b[index].abs()) * vector.abs();
            }
            let relative = if scale > 0.0 {
                value.abs() / scale
            } else {
                value.abs()
            };
            if !relative.is_finite() || relative > budget.relative_residual_tolerance {
                return Err(FilmSpectrumError::Residual {
                    relative,
                    tolerance: budget.relative_residual_tolerance,
                });
            }
            maximum = maximum.max(relative);
        }
        Ok(maximum)
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}
fn divide(value: C64, scale: f64) -> C64 {
    C64::new(value.re / scale, value.im / scale)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), FilmSpectrumError> {
    if cancelled() {
        Err(FilmSpectrumError::Cancelled)
    } else {
        Ok(())
    }
}
fn lu_error(error: LuComplexError) -> FilmSpectrumError {
    match error {
        LuComplexError::Cancelled => FilmSpectrumError::Cancelled,
        LuComplexError::Singular { index } => FilmSpectrumError::Singular { index },
    }
}

#[cfg(test)]
mod tests;
