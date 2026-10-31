//! A bounded dense maximizer of `a^T G a + 2 h^T a`, `||a||_2 <= 1`.
//!
//! For numerically admitted PSD G, the global trust-region conditions are
//! `(lambda I-G)a=h`, `lambda I-G >= 0`, feasibility and complementarity.
//! The regular branch solves the secular equation in `delta=lambda-lambda_max`;
//! the hard branch fills the remaining norm in the top eigenspace. This is not
//! Steihaug iteration: zero linear forcing must still find the adverse mode.
//! See Adachi et al., "Solving the trust-region subproblem by a generalized
//! eigenvalue problem", Theorem 1.1, for the underlying optimality conditions.
//!
//! Reports contain numerical residuals and a numerical primal/dual comparison,
//! not interval enclosures or an exact global-optimality certificate. No matrix
//! eigenvalues are clipped or diagonally shifted to hide indefiniteness.
use std::ops::ControlFlow;

/// Bounds the existing non-interruptible cyclic-Jacobi kernel.
pub const MAX_QUADRATIC_BALL_DIMENSION: usize = 16;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QuadraticBallOptions {
    /// Secular iterations, in 1..=1024. Jacobi has its own fixed sweep bound.
    pub max_iterations: usize,
    /// Dimensionless symmetry, PSD, eigenpair and final numerical gates.
    pub tolerance: f64,
    /// Minimum relative positive shift for a differentiable-envelope consumer.
    /// Must be at least tolerance. This is a numerical conditioning gate.
    pub separation_tolerance: f64,
}
impl Default for QuadraticBallOptions {
    fn default() -> Self {
        Self {
            max_iterations: 256,
            tolerance: 1e-10,
            separation_tolerance: 1e-7,
        }
    }
}
impl QuadraticBallOptions {
    pub fn validate(self) -> Result<(), QuadraticBallError> {
        if !(1..=1024).contains(&self.max_iterations)
            || !self.tolerance.is_finite()
            || self.tolerance < 128.0 * f64::EPSILON
            || self.tolerance > 1e-3
            || !self.separation_tolerance.is_finite()
            || self.separation_tolerance < self.tolerance
            || self.separation_tolerance >= 1.0
        {
            return Err(QuadraticBallError::Invalid(
                "invalid quadratic-ball numerical policy",
            ));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QuadraticBallBranch {
    Regular,
    /// The linear term is orthogonal to the numerical top eigenspace and the
    /// singular shifted solve leaves room for a top-eigenspace component.
    HardCase,
    /// G and h are both exactly zero; every feasible point is a maximizer.
    Flat,
}

#[derive(Debug, Clone, PartialEq)]
pub enum QuadraticBallError {
    Invalid(&'static str),
    Cancelled,
    Spectrum(&'static str),
    Unrepresentable(&'static str),
    NoConvergence { iterations: usize },
    Evidence(&'static str),
}
impl std::fmt::Display for QuadraticBallError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "quadratic-ball maximum refused: {self:?}")
    }
}
impl std::error::Error for QuadraticBallError {}

#[derive(Debug, Clone, PartialEq)]
pub struct QuadraticBallMaximum {
    /// One deterministic numerical witness in the caller's coefficient basis.
    pub coefficients: Vec<f64>,
    /// Actual admitted symmetric quadratic at the returned witness; no constant.
    pub value: f64,
    /// Multiplier for `1-a^T a`, in the supplied objective units.
    pub multiplier: f64,
    /// Numerical shifted-system dual `lambda + h^T (lambda I-G)^dagger h`.
    /// Its floating-point evaluation is not a certified upper enclosure.
    pub dual_value: f64,
    /// Signed dual_value-value, without clipping small negative roundoff.
    pub duality_gap: f64,
    /// Absolute gap evaluated before undoing objective normalization.
    pub relative_duality_gap: f64,
    pub branch: QuadraticBallBranch,
    /// Largest absolute entry of the original G and h (zero for Flat).
    pub normalization: f64,
    /// `(lambda-lambda_max)/max(lambda, ||h||_2)` in normalized arithmetic.
    pub relative_separation: f64,
    /// Actual `(lambda I-G)a-h` norm relative to the reported evidence scale.
    pub stationarity_residual: f64,
    /// Positive part of `a^T a-1` at the returned witness.
    pub feasibility_residual: f64,
    pub complementarity_residual: f64,
    /// Max entry residual for G normalized by its own largest absolute entry.
    pub eigenpair_residual: f64,
    pub orthogonality_defect: f64,
    pub symmetry_defect: f64,
    /// Original-unit minimum numerical eigenvalue, never clipped to zero.
    pub minimum_eigenvalue: f64,
    pub iterations: usize,
    gradient_admissible: bool,
}
impl QuadraticBallMaximum {
    /// True only for an admitted regular branch with a resolved positive shift.
    /// Hard/flat/near-hard witnesses remain usable values, not smooth gradients.
    #[must_use]
    pub const fn gradient_admissible(&self) -> bool {
        self.gradient_admissible
    }
}

fn poll(checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), QuadraticBallError> {
    if checkpoint().is_break() {
        Err(QuadraticBallError::Cancelled)
    } else {
        Ok(())
    }
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn norm(a: &[f64]) -> f64 {
    let scale = a.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
    if scale == 0.0 || !scale.is_finite() {
        return scale;
    }
    scale * fs_math::det::sqrt(a.iter().map(|x| (x / scale) * (x / scale)).sum())
}
fn rescale(value: f64, scale: f64) -> Result<f64, QuadraticBallError> {
    let result = value * scale;
    if !result.is_finite() || (value != 0.0 && result == 0.0) {
        Err(QuadraticBallError::Unrepresentable(
            "objective or multiplier range",
        ))
    } else {
        Ok(result)
    }
}

// Stable norm of the secular solution. An overflowing trial is simply on the
// infeasible side of the bracket, never a usable field or objective.
fn coordinates(g: &[f64], gaps: &[f64], delta: f64) -> Vec<f64> {
    g.iter()
        .zip(gaps)
        .map(|(&g, &gap)| if g == 0.0 { 0.0 } else { g / (gap + delta) })
        .collect()
}
fn secular(
    g: &[f64],
    gaps: &[f64],
    options: QuadraticBallOptions,
    checkpoint: &mut impl FnMut() -> ControlFlow<()>,
) -> Result<(Vec<f64>, f64, QuadraticBallBranch, usize), QuadraticBallError> {
    let top: Vec<_> = g
        .iter()
        .zip(gaps)
        .map(|(&g, &gap)| if gap == 0.0 { g } else { 0.0 })
        .collect();
    let top_norm = norm(&top);
    if top_norm == 0.0 {
        let mut y = coordinates(g, gaps, 0.0);
        let length = norm(&y);
        if length <= 1.0 {
            let index = gaps
                .iter()
                .position(|gap| *gap == 0.0)
                .expect("largest eigenvalue has zero gap");
            y[index] = fs_math::det::sqrt((1.0 - length) * (1.0 + length));
            return Ok((y, 0.0, QuadraticBallBranch::HardCase, 0));
        }
    }
    let mut low = top_norm;
    let mut high = norm(g);
    let target = (0.01 * options.tolerance).min(1e-12);
    let upper = coordinates(g, gaps, high);
    if (norm(&upper) - 1.0).abs() <= target {
        return Ok((upper, high, QuadraticBallBranch::Regular, 0));
    }
    for iteration in 1..=options.max_iterations {
        poll(checkpoint)?;
        // A nonzero top projection can be arbitrarily smaller than the other
        // components. Geometric bracketing resolves its exponent range first.
        let mut middle = if low > 0.0 && high > 4.0 * low {
            fs_math::det::sqrt(low) * fs_math::det::sqrt(high)
        } else {
            0.5 * low + 0.5 * high
        };
        if middle <= low || middle >= high {
            middle = 0.5 * low + 0.5 * high;
        }
        if middle <= low || middle >= high {
            return Ok((
                coordinates(g, gaps, high),
                high,
                QuadraticBallBranch::Regular,
                iteration,
            ));
        }
        let y = coordinates(g, gaps, middle);
        let length = norm(&y);
        if (length - 1.0).abs() <= target {
            return Ok((y, middle, QuadraticBallBranch::Regular, iteration));
        }
        if length > 1.0 {
            low = middle;
        } else {
            high = middle;
        }
    }
    Err(QuadraticBallError::NoConvergence {
        iterations: options.max_iterations,
    })
}

/// Numerically solve the global PSD quadratic maximization on the unit ball.
///
/// G is symmetric row-major, h has 1..=16 entries, and all inputs are finite.
/// Symmetry and PSD are admitted relative to G's own scale, independently of
/// h. Tiny negative numerical eigenvalues within the policy are reported and
/// used unchanged; larger negative eigenvalues refuse. Symmetry is checked
/// before using the symmetric part, which defines the same real quadratic.
/// Common objective normalization supports uniform extreme scales without
/// forming an overflowing reciprocal. Unrepresentable range loss refuses.
///
/// The existing fs-la Jacobi kernel is polled immediately before and after its
/// dimension-bounded call. Secular iterations and other loops poll separately.
/// There is no callback inside Jacobi and no arbitrary-conditioning guarantee.
/// A hard-case sign is fixed by the first largest-magnitude entry of its native
/// eigenvector; repeated eigenvalues do not acquire a unique direction claim.
pub fn maximize_quadratic_ball(
    gram: &[f64],
    linear: &[f64],
    options: QuadraticBallOptions,
    mut checkpoint: impl FnMut() -> ControlFlow<()>,
) -> Result<QuadraticBallMaximum, QuadraticBallError> {
    poll(&mut checkpoint)?;
    options.validate()?;
    let n = linear.len();
    if n == 0
        || n > MAX_QUADRATIC_BALL_DIMENSION
        || gram.len() != n * n
        || gram.iter().chain(linear).any(|v| !v.is_finite())
    {
        return Err(QuadraticBallError::Invalid(
            "finite square G and 1..=16 linear entries required",
        ));
    }
    let gram_scale = gram.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    let normalization = linear.iter().map(|v| v.abs()).fold(gram_scale, f64::max);
    if normalization == 0.0 {
        let mut coefficients = vec![0.0; n];
        coefficients[0] = 1.0;
        poll(&mut checkpoint)?;
        return Ok(QuadraticBallMaximum {
            coefficients,
            value: 0.0,
            multiplier: 0.0,
            dual_value: 0.0,
            duality_gap: 0.0,
            relative_duality_gap: 0.0,
            branch: QuadraticBallBranch::Flat,
            normalization,
            relative_separation: 0.0,
            stationarity_residual: 0.0,
            feasibility_residual: 0.0,
            complementarity_residual: 0.0,
            eigenpair_residual: 0.0,
            orthogonality_defect: 0.0,
            symmetry_defect: 0.0,
            minimum_eigenvalue: 0.0,
            iterations: 0,
            gradient_admissible: false,
        });
    }
    let mut matrix = Vec::with_capacity(n * n);
    for &value in gram {
        let scaled = if gram_scale == 0.0 {
            0.0
        } else {
            value / gram_scale
        };
        if value != 0.0 && scaled == 0.0 {
            return Err(QuadraticBallError::Unrepresentable(
                "Gram normalization underflow",
            ));
        }
        matrix.push(scaled);
    }
    let mut symmetry_defect = 0.0_f64;
    for i in 0..n {
        poll(&mut checkpoint)?;
        for j in i + 1..n {
            symmetry_defect = symmetry_defect.max((matrix[i * n + j] - matrix[j * n + i]).abs());
            let value = 0.5 * matrix[i * n + j] + 0.5 * matrix[j * n + i];
            matrix[i * n + j] = value;
            matrix[j * n + i] = value;
        }
    }
    if symmetry_defect > options.tolerance {
        return Err(QuadraticBallError::Spectrum("nonsymmetric Gram matrix"));
    }
    poll(&mut checkpoint)?;
    let (eigenvalues, mut vectors) = fs_la::eigen::jacobi_eigh(&matrix, n);
    poll(&mut checkpoint)?;
    if eigenvalues.iter().chain(&vectors).any(|v| !v.is_finite())
        || eigenvalues[0] < -options.tolerance
        || eigenvalues[n - 1] < 0.0
    {
        return Err(QuadraticBallError::Spectrum(
            "Gram matrix is not numerically positive semidefinite",
        ));
    }
    for k in 0..n {
        let mut pivot = 0;
        for i in 1..n {
            if vectors[i * n + k].abs() > vectors[pivot * n + k].abs() {
                pivot = i;
            }
        }
        if vectors[pivot * n + k] < 0.0 {
            for i in 0..n {
                vectors[i * n + k] = -vectors[i * n + k];
            }
        }
    }
    let mut eigenpair_residual = 0.0_f64;
    let mut orthogonality_defect = 0.0_f64;
    for i in 0..n {
        poll(&mut checkpoint)?;
        for k in 0..n {
            let applied: f64 = (0..n).map(|j| matrix[i * n + j] * vectors[j * n + k]).sum();
            eigenpair_residual =
                eigenpair_residual.max((applied - eigenvalues[k] * vectors[i * n + k]).abs());
            let inner: f64 = (0..n)
                .map(|j| vectors[j * n + i] * vectors[j * n + k])
                .sum();
            orthogonality_defect =
                orthogonality_defect.max((inner - if i == k { 1.0 } else { 0.0 }).abs());
        }
    }
    if eigenpair_residual > options.tolerance || orthogonality_defect > options.tolerance {
        return Err(QuadraticBallError::Spectrum(
            "unresolved eigenpairs or orthogonality",
        ));
    }
    let ratio = gram_scale / normalization;
    if gram_scale != 0.0 && ratio == 0.0 {
        return Err(QuadraticBallError::Unrepresentable(
            "relative Gram scale underflow",
        ));
    }
    let mut h = Vec::with_capacity(n);
    for &value in linear {
        let scaled = value / normalization;
        if value != 0.0 && scaled == 0.0 {
            return Err(QuadraticBallError::Unrepresentable(
                "linear normalization underflow",
            ));
        }
        h.push(scaled);
    }
    for value in &mut matrix {
        let scaled = *value * ratio;
        if *value != 0.0 && scaled == 0.0 {
            return Err(QuadraticBallError::Unrepresentable(
                "relative Gram entry underflow",
            ));
        }
        *value = scaled;
    }
    let largest = eigenvalues[n - 1] * ratio;
    let gaps: Vec<_> = eigenvalues
        .iter()
        .map(|d| (eigenvalues[n - 1] - d) * ratio)
        .collect();
    let g: Vec<_> = (0..n)
        .map(|k| (0..n).map(|i| vectors[i * n + k] * h[i]).sum())
        .collect();
    if norm(&g) == 0.0 && h.iter().any(|h| *h != 0.0) {
        return Err(QuadraticBallError::Unrepresentable(
            "linear eigenspace projection collapsed",
        ));
    }
    let (y, delta, branch, iterations) = secular(&g, &gaps, options, &mut checkpoint)?;
    let multiplier_scaled = largest + delta;
    let mut coefficients: Vec<_> = (0..n)
        .map(|i| (0..n).map(|k| vectors[i * n + k] * y[k]).sum())
        .collect();
    let length = norm(&coefficients);
    if !length.is_finite() || length == 0.0 {
        return Err(QuadraticBallError::Evidence("invalid secular witness"));
    }
    for a in &mut coefficients {
        *a /= length;
    }
    // Make an overshoot from the final basis rotation strictly feasible rather
    // than publishing an out-of-ball force. The actual KKT/gap gates follow.
    for _ in 0..3 {
        let square = dot(&coefficients, &coefficients);
        if square <= 1.0 {
            break;
        }
        let length = fs_math::det::sqrt(square);
        let upper = f64::from_bits(length.to_bits() + 1);
        for a in &mut coefficients {
            *a /= upper;
        }
    }
    let square = dot(&coefficients, &coefficients);
    if square > 1.0 {
        return Err(QuadraticBallError::Evidence("witness exceeds unit ball"));
    }
    let applied: Vec<_> = (0..n)
        .map(|i| dot(&matrix[i * n..(i + 1) * n], &coefficients))
        .collect();
    let value_scaled = dot(&coefficients, &applied) + 2.0 * dot(&h, &coefficients);
    let dual_scaled = multiplier_scaled + dot(&g, &y);
    let hnorm = norm(&h);
    let matrix_norm = matrix
        .chunks_exact(n)
        .map(|row| row.iter().map(|v| v.abs()).sum::<f64>())
        .fold(0.0_f64, f64::max);
    let evidence_scale = value_scaled
        .abs()
        .max(dual_scaled.abs())
        .max(multiplier_scaled.abs())
        .max(matrix_norm)
        .max(hnorm)
        .max(f64::MIN_POSITIVE);
    let residual: Vec<_> = coefficients
        .iter()
        .zip(&applied)
        .zip(&h)
        .map(|((&a, &ga), &h)| multiplier_scaled * a - ga - h)
        .collect();
    let stationarity_residual = norm(&residual) / evidence_scale;
    let feasibility_residual = (square - 1.0).max(0.0);
    let complementarity_residual = (multiplier_scaled * (1.0 - square)).abs() / evidence_scale;
    let relative_duality_gap = (dual_scaled - value_scaled).abs() / evidence_scale;
    if [
        value_scaled,
        dual_scaled,
        stationarity_residual,
        feasibility_residual,
        complementarity_residual,
        relative_duality_gap,
    ]
    .iter()
    .any(|v| !v.is_finite())
        || [
            stationarity_residual,
            feasibility_residual,
            complementarity_residual,
            relative_duality_gap,
        ]
        .iter()
        .any(|v| *v > options.tolerance)
    {
        return Err(QuadraticBallError::Evidence(
            "quadratic-ball KKT or primal-dual gate failed",
        ));
    }
    let shift = multiplier_scaled - largest;
    let relative_separation = shift / multiplier_scaled.max(hnorm).max(f64::MIN_POSITIVE);
    let spectral_uncertainty =
        8.0 * n as f64 * (eigenpair_residual * ratio + orthogonality_defect * largest);
    let gradient_admissible = branch == QuadraticBallBranch::Regular
        && relative_separation > options.separation_tolerance
        && shift > spectral_uncertainty;
    let value = rescale(value_scaled, normalization)?;
    let dual_value = rescale(dual_scaled, normalization)?;
    let multiplier = rescale(multiplier_scaled, normalization)?;
    let minimum_eigenvalue = eigenvalues[0] * gram_scale;
    if !minimum_eigenvalue.is_finite() {
        return Err(QuadraticBallError::Unrepresentable(
            "minimum eigenvalue range",
        ));
    }
    poll(&mut checkpoint)?;
    Ok(QuadraticBallMaximum {
        coefficients,
        value,
        multiplier,
        dual_value,
        duality_gap: dual_value - value,
        relative_duality_gap,
        branch,
        normalization,
        relative_separation,
        stationarity_residual,
        feasibility_residual,
        complementarity_residual,
        eigenpair_residual,
        orthogonality_defect,
        symmetry_defect,
        minimum_eigenvalue,
        iterations,
        gradient_admissible,
    })
}
