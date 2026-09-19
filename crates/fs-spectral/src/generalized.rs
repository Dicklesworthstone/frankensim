//! Matrix-free low modes of an SPD pencil A x = lambda B x.
//!
//! Block inverse iteration uses caller-owned A solves, twice-reorthogonalized
//! B-inner products, and fs-la's small Rayleigh--Ritz kernel. It never treats
//! B^-1 A as Euclidean symmetric, factors a global mass matrix, or silently
//! substitutes lumped mass. Storage is O(n * block_size), not O(n^2).
//!
//! Residuals are recomputed with the original operators. They are numerical
//! estimates, not certified spectral enclosures or proof that no lower mode
//! was missed. A and B must be fixed, symmetric and positive definite on the
//! caller's admitted space. Remove constraints/nullspaces before adapting it.

use fs_la::eigen::{MAX_EIGEN_UNPOLLED_SCALAR_WORK, MAX_EIGEN_WORK_ELEMENTS, jacobi_eigh};

/// Domain-neutral generalized operator. Each call must bound its own work;
/// the domain owns linear-solve tolerances, budgets, and cancellation.
pub trait GeneralizedOp {
    /// Domain error, preserved without replacing it by a partial eigenresult.
    type Error;
    /// Number of unconstrained unknowns.
    fn dim(&self) -> usize;
    /// Write A x into y.
    fn apply_a(&mut self, x: &[f64], y: &mut [f64]) -> Result<(), Self::Error>;
    /// Write B x into y (B is the actual positive-definite metric).
    fn apply_b(&mut self, x: &[f64], y: &mut [f64]) -> Result<(), Self::Error>;
    /// Solve A y = rhs under the caller's accuracy and work budgets.
    fn solve_a(&mut self, rhs: &[f64], y: &mut [f64]) -> Result<(), Self::Error>;
    /// Poll before each bounded algebraic stage and before publishing state.
    fn checkpoint(&mut self) -> Result<(), Self::Error>;
}

/// Explicit block, accuracy, and deterministic initialization settings.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeneralizedOptions {
    /// Number of low Ritz pairs requested.
    pub count: usize,
    /// Search-block width; oversampling helps resolve clustered low modes.
    pub block_size: usize,
    /// Acceptance threshold for ||A x - lambda B x|| / (||A x|| + |lambda| ||B x||).
    pub relative_tolerance: f64,
    /// Initialization seed, independent of thread scheduling.
    pub seed: u64,
}

/// Numerical failure or a caller-owned interruption.
#[derive(Debug, Clone, PartialEq)]
pub enum GeneralizedError<E> {
    /// Invalid dimensions, tolerance, or impractical aggregate work.
    InvalidInput(&'static str),
    /// Operator, linear-budget, or cancellation failure from the domain.
    Operator(E),
    /// Non-finite arithmetic, nonpositive metric, or dependent search block.
    Breakdown(&'static str),
    /// No complete, residual-qualified answer within the supplied step budget.
    Unconverged {
        /// Accepted outer steps, including earlier resumed runs.
        iterations: usize,
        /// Worst relative residual among the requested modes.
        worst_residual: f64,
    },
}

impl<E: std::fmt::Display> std::fmt::Display for GeneralizedError<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidInput(message) | Self::Breakdown(message) => f.write_str(message),
            Self::Operator(error) => std::fmt::Display::fmt(error, f),
            Self::Unconverged { iterations, worst_residual } => write!(f,
                "generalized eigen solve unconverged after {iterations} steps (residual {worst_residual:e})"),
        }
    }
}

impl<E: std::error::Error + 'static> std::error::Error for GeneralizedError<E> {}

/// A complete set of requested numerical Ritz estimates, in ascending order.
#[derive(Debug, Clone, PartialEq)]
pub struct GeneralizedReport {
    /// Eigenvalues (squared angular frequencies for an elasticity pencil).
    pub values: Vec<f64>,
    /// B-normalized vectors in the caller's unconstrained coordinate space.
    pub vectors: Vec<Vec<f64>>,
    /// Recomputed, scale-relative pencil residuals in the same order.
    pub relative_residuals: Vec<f64>,
    /// Accepted block inverse iterations.
    pub iterations: usize,
}

/// Resumable accepted state. Clone is a checkpoint. A failed step does not
/// alter the last accepted block; operator-side work accounting is NOT rolled
/// back. Do not resume with different A or B: start a new state with warm seeds.
#[derive(Debug, Clone, PartialEq)]
pub struct GeneralizedState {
    n: usize,
    options: GeneralizedOptions,
    vectors: Vec<Vec<f64>>,
    metric_vectors: Vec<Vec<f64>>,
    values: Vec<f64>,
    residuals: Vec<f64>,
    iterations: usize,
}

impl GeneralizedState {
    /// Initialize a deterministic full search block. Optional warm seeds must
    /// provide exactly `block_size` finite, independent vectors of length n.
    pub fn new<O: GeneralizedOp>(
        op: &mut O,
        options: GeneralizedOptions,
        seeds: Option<&[Vec<f64>]>,
    ) -> Result<Self, GeneralizedError<O::Error>> {
        let n = op.dim();
        let b = options.block_size;
        if n == 0 || options.count == 0 || options.count > b || b > n || b > 64 {
            return Err(GeneralizedError::InvalidInput("require 1 <= count <= block_size <= min(n, 64)"));
        }
        if !options.relative_tolerance.is_finite()
            || options.relative_tolerance <= 0.0 || options.relative_tolerance >= 1.0 {
            return Err(GeneralizedError::InvalidInput("relative tolerance must lie in (0, 1)"));
        }
        // Retained state, candidate/rotated blocks, small dense work and scratch.
        let storage = n.checked_mul(b).and_then(|v| v.checked_mul(10))
            .and_then(|v| v.checked_add(8 * b * b));
        let work = n.checked_mul(b).and_then(|v| v.checked_mul(b))
            .and_then(|v| v.checked_mul(16));
        if storage.is_none_or(|v| v > MAX_EIGEN_WORK_ELEMENTS)
            || work.is_none_or(|v| v > MAX_EIGEN_UNPOLLED_SCALAR_WORK) {
            return Err(GeneralizedError::InvalidInput("generalized search block exceeds work/memory cap"));
        }
        op.checkpoint().map_err(GeneralizedError::Operator)?;
        let mut vectors = if let Some(seeds) = seeds {
            if seeds.len() != b || seeds.iter().any(|v| v.len() != n || !finite(v)) {
                return Err(GeneralizedError::InvalidInput("warm-start block shape or values are invalid"));
            }
            seeds.to_vec()
        } else {
            let mut key = options.seed;
            (0..b).map(|_| (0..n).map(|_| {
                // SplitMix64: fixed integer arithmetic and a 53-bit f64 map.
                key = key.wrapping_add(0x9e37_79b9_7f4a_7c15);
                let mut z = key;
                z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
                z ^= z >> 31;
                (z >> 11) as f64 * (1.0 / 4_503_599_627_370_496.0) - 1.0
            }).collect()).collect()
        };
        let (metric_vectors, values, residuals) = extract(op, &mut vectors)?;
        op.checkpoint().map_err(GeneralizedError::Operator)?;
        Ok(Self { n, options, vectors, metric_vectors, values, residuals, iterations: 0 })
    }

    /// Perform one transactional block step. Returns true only when all
    /// requested modes satisfy the recomputed residual threshold.
    pub fn step<O: GeneralizedOp>(&mut self, op: &mut O) -> Result<bool, GeneralizedError<O::Error>> {
        if op.dim() != self.n {
            return Err(GeneralizedError::InvalidInput("operator dimension changed during continuation"));
        }
        op.checkpoint().map_err(GeneralizedError::Operator)?;
        let iterations = self.iterations.checked_add(1)
            .ok_or(GeneralizedError::InvalidInput("iteration count overflow"))?;
        let mut vectors = Vec::with_capacity(self.options.block_size);
        for rhs in &self.metric_vectors {
            op.checkpoint().map_err(GeneralizedError::Operator)?;
            let mut vector = vec![f64::NAN; self.n];
            op.solve_a(rhs, &mut vector).map_err(GeneralizedError::Operator)?;
            if !finite(&vector) {
                return Err(GeneralizedError::Breakdown("inverse application is non-finite"));
            }
            vectors.push(vector);
        }
        let (metric_vectors, values, residuals) = extract(op, &mut vectors)?;
        op.checkpoint().map_err(GeneralizedError::Operator)?;
        self.vectors = vectors;
        self.metric_vectors = metric_vectors;
        self.values = values;
        self.residuals = residuals;
        self.iterations = iterations;
        Ok(self.converged())
    }

    /// Run at most `additional_steps`, preserving resumability on exhaustion.
    pub fn run<O: GeneralizedOp>(
        &mut self, op: &mut O, additional_steps: usize,
    ) -> Result<GeneralizedReport, GeneralizedError<O::Error>> {
        if op.dim() != self.n {
            return Err(GeneralizedError::InvalidInput("operator dimension changed during continuation"));
        }
        op.checkpoint().map_err(GeneralizedError::Operator)?;
        for _ in 0..additional_steps {
            if self.converged() { break; }
            self.step(op)?;
        }
        op.checkpoint().map_err(GeneralizedError::Operator)?;
        if !self.converged() {
            return Err(GeneralizedError::Unconverged {
                iterations: self.iterations,
                worst_residual: self.residuals[..self.options.count].iter().copied().fold(0.0, f64::max),
            });
        }
        let k = self.options.count;
        Ok(GeneralizedReport {
            values: self.values[..k].to_vec(), vectors: self.vectors[..k].to_vec(),
            relative_residuals: self.residuals[..k].to_vec(), iterations: self.iterations,
        })
    }

    /// Accepted full search block, suitable as explicit seeds for a new pencil.
    #[must_use]
    pub fn warm_start(&self) -> &[Vec<f64>] { &self.vectors }

    fn converged(&self) -> bool {
        self.residuals[..self.options.count].iter().all(|r| *r <= self.options.relative_tolerance)
    }
}

fn finite(x: &[f64]) -> bool { x.iter().all(|v| v.is_finite()) }
fn dot(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y).fold(0.0, |sum, (&x, &y)| x.mul_add(y, sum))
}
fn norm(x: &[f64]) -> f64 { x.iter().fold(0.0_f64, |sum, &v| sum.hypot(v)) }

type Extraction = (Vec<Vec<f64>>, Vec<f64>, Vec<f64>);

fn extract<O: GeneralizedOp>(
    op: &mut O, vectors: &mut Vec<Vec<f64>>,
) -> Result<Extraction, GeneralizedError<O::Error>> {
    let n = op.dim();
    let b = vectors.len();
    let mut metric: Vec<Vec<f64>> = Vec::with_capacity(b);
    for i in 0..b {
        op.checkpoint().map_err(GeneralizedError::Operator)?;
        let mut bv = vec![f64::NAN; n];
        op.apply_b(&vectors[i], &mut bv).map_err(GeneralizedError::Operator)?;
        let original = dot(&vectors[i], &bv);
        if !finite(&bv) || !original.is_finite() || original <= 0.0 {
            return Err(GeneralizedError::Breakdown("mass metric is not positive on the search block"));
        }
        for _ in 0..2 {
            for j in 0..i {
                let coefficient = dot(&vectors[j], &bv);
                for d in 0..n {
                    let correction = coefficient * vectors[j][d];
                    vectors[i][d] -= correction;
                    bv[d] -= coefficient * metric[j][d];
                }
            }
        }
        // Reapply B instead of trusting the orthogonalization recurrence.
        bv.fill(f64::NAN);
        op.apply_b(&vectors[i], &mut bv).map_err(GeneralizedError::Operator)?;
        let length_squared = dot(&vectors[i], &bv);
        if !finite(&bv) || !length_squared.is_finite()
            || length_squared <= 128.0 * f64::EPSILON * f64::EPSILON * original {
            return Err(GeneralizedError::Breakdown("mass-orthogonal search block lost rank"));
        }
        let length = length_squared.sqrt();
        for d in 0..n { vectors[i][d] /= length; bv[d] /= length; }
        metric.push(bv);
    }
    let mut av = Vec::with_capacity(b);
    for vector in vectors.iter() {
        op.checkpoint().map_err(GeneralizedError::Operator)?;
        let mut output = vec![f64::NAN; n];
        op.apply_a(vector, &mut output).map_err(GeneralizedError::Operator)?;
        if !finite(&output) { return Err(GeneralizedError::Breakdown("stiffness application is non-finite")); }
        av.push(output);
    }
    let mut projected = vec![0.0; b * b];
    for i in 0..b {
        for j in i..b {
            let value = 0.5 * dot(&vectors[i], &av[j]) + 0.5 * dot(&vectors[j], &av[i]);
            projected[i * b + j] = value;
            projected[j * b + i] = value;
        }
    }
    let scale = projected.iter().map(|v| v.abs()).fold(0.0, f64::max);
    if !finite(&projected) || scale <= 0.0 {
        return Err(GeneralizedError::Breakdown("invalid projected stiffness"));
    }
    for value in &mut projected { *value /= scale; }
    op.checkpoint().map_err(GeneralizedError::Operator)?;
    let (mut values, rotation) = jacobi_eigh(&projected, b);
    let mut rotated = vec![vec![0.0; n]; b];
    let mut residuals = Vec::with_capacity(b);
    for k in 0..b {
        values[k] *= scale;
        if !values[k].is_finite() || values[k] <= 0.0 {
            return Err(GeneralizedError::Breakdown("nonpositive Ritz stiffness: check supports/nullspace"));
        }
        for j in 0..b {
            for d in 0..n {
                rotated[k][d] = rotation[j * b + k].mul_add(vectors[j][d], rotated[k][d]);
            }
        }
        op.checkpoint().map_err(GeneralizedError::Operator)?;
        metric[k].fill(f64::NAN);
        op.apply_b(&rotated[k], &mut metric[k]).map_err(GeneralizedError::Operator)?;
        let length_squared = dot(&rotated[k], &metric[k]);
        if !length_squared.is_finite() || length_squared <= 0.0 {
            return Err(GeneralizedError::Breakdown("invalid Ritz mass normalization"));
        }
        let length = length_squared.sqrt();
        for d in 0..n { rotated[k][d] /= length; metric[k][d] /= length; }
        av[k].fill(f64::NAN);
        op.apply_a(&rotated[k], &mut av[k]).map_err(GeneralizedError::Operator)?;
        if !finite(&av[k]) || !finite(&metric[k]) || !finite(&rotated[k]) {
            return Err(GeneralizedError::Breakdown("non-finite Ritz vector or operator image"));
        }
        // Scale the vectors before forming either the residual or denominator.
        let magnitude = av[k].iter().map(|v| v.abs())
            .chain(metric[k].iter().map(|v| (values[k] * v).abs())).fold(0.0, f64::max);
        if !magnitude.is_finite() || magnitude <= 0.0 {
            return Err(GeneralizedError::Breakdown("invalid residual scale"));
        }
        let ka: Vec<f64> = av[k].iter().map(|v| v / magnitude).collect();
        let mb: Vec<f64> = metric[k].iter().map(|v| values[k] * v / magnitude).collect();
        let residual: Vec<f64> = ka.iter().zip(&mb).map(|(a, b)| a - b).collect();
        let relative = norm(&residual) / (norm(&ka) + norm(&mb));
        if !relative.is_finite() { return Err(GeneralizedError::Breakdown("non-finite pencil residual")); }
        residuals.push(relative);
    }
    *vectors = rotated;
    Ok((metric, values, residuals))
}

#[cfg(test)]
mod tests {
    use super::*;

    // A = S^T diag(lambda) S, B = S^T S, where S is upper bidiagonal.
    // These noncommuting SPD matrices have the prescribed generalized roots.
    struct Pencil { lambda: Vec<f64>, fail_solve: bool }
    impl Pencil {
        fn s(x: &[f64]) -> Vec<f64> {
            (0..x.len()).map(|i| x[i] + if i + 1 < x.len() { 0.25 * x[i + 1] } else { 0.0 }).collect()
        }
        fn st(x: &[f64], y: &mut [f64]) {
            for i in 0..x.len() { y[i] = x[i] + if i > 0 { 0.25 * x[i - 1] } else { 0.0 }; }
        }
    }
    impl GeneralizedOp for Pencil {
        type Error = &'static str;
        fn dim(&self) -> usize { self.lambda.len() }
        fn apply_a(&mut self, x: &[f64], y: &mut [f64]) -> Result<(), Self::Error> {
            let mut sx = Self::s(x);
            for (v, lambda) in sx.iter_mut().zip(&self.lambda) { *v *= lambda; }
            Self::st(&sx, y); Ok(())
        }
        fn apply_b(&mut self, x: &[f64], y: &mut [f64]) -> Result<(), Self::Error> {
            Self::st(&Self::s(x), y); Ok(())
        }
        fn solve_a(&mut self, rhs: &[f64], y: &mut [f64]) -> Result<(), Self::Error> {
            if self.fail_solve { return Err("cancelled linear solve"); }
            y.copy_from_slice(rhs);
            for i in 1..y.len() { let correction = 0.25 * y[i - 1]; y[i] -= correction; }
            for (v, lambda) in y.iter_mut().zip(&self.lambda) { *v /= lambda; }
            for i in (0..y.len() - 1).rev() { let correction = 0.25 * y[i + 1]; y[i] -= correction; }
            Ok(())
        }
        fn checkpoint(&mut self) -> Result<(), Self::Error> { Ok(()) }
    }
    fn options() -> GeneralizedOptions {
        GeneralizedOptions { count: 3, block_size: 5, relative_tolerance: 1e-10, seed: 42 }
    }
    fn pencil() -> Pencil {
        Pencil { lambda: vec![1.0, 1.0, 2.0, 4.0, 7.0, 12.0, 16.0, 20.0], fail_solve: false }
    }

    #[test]
    fn noncommuting_metric_and_repeated_low_modes() {
        let mut op = pencil();
        let mut state = GeneralizedState::new(&mut op, options(), None).unwrap();
        let report = state.run(&mut op, 100).unwrap();
        for (actual, expected) in report.values.iter().zip([1.0, 1.0, 2.0]) {
            assert!((actual - expected).abs() < 1e-9);
        }
        for i in 0..3 {
            let mut bv = vec![0.0; op.dim()];
            op.apply_b(&report.vectors[i], &mut bv).unwrap();
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((dot(&report.vectors[j], &bv) - expected).abs() < 1e-10);
            }
        }
        assert!(report.relative_residuals.iter().all(|r| *r < 1e-10));
    }

    #[test]
    fn split_steps_replay_and_failures_preserve_accepted_state() {
        let mut op = pencil();
        let mut whole = GeneralizedState::new(&mut op, options(), None).unwrap();
        let mut split = whole.clone();
        for _ in 0..6 { whole.step(&mut op).unwrap(); }
        for _ in 0..2 { split.step(&mut op).unwrap(); }
        let saved = split.clone();
        op.fail_solve = true;
        assert!(matches!(split.step(&mut op), Err(GeneralizedError::Operator(_))));
        assert_eq!(split, saved);
        op.fail_solve = false;
        for _ in 0..4 { split.step(&mut op).unwrap(); }
        assert_eq!(whole, split);
    }

    #[test]
    fn no_partial_success_on_zero_budget_or_rank_deficient_seeds() {
        let mut op = pencil();
        let mut state = GeneralizedState::new(&mut op, options(), None).unwrap();
        assert!(matches!(state.run(&mut op, 0), Err(GeneralizedError::Unconverged { iterations: 0, .. })));
        let seeds = vec![vec![1.0; op.dim()]; options().block_size];
        assert!(matches!(GeneralizedState::new(&mut op, options(), Some(&seeds)),
            Err(GeneralizedError::Breakdown(_))));
    }

    #[test]
    fn full_block_and_common_physical_scaling() {
        for scale in [1e-12, 1.0, 1e12] {
            let mut op = pencil();
            for lambda in &mut op.lambda { *lambda *= scale; }
            let mut settings = options(); settings.block_size = op.dim();
            let mut state = GeneralizedState::new(&mut op, settings, None).unwrap();
            let report = state.run(&mut op, 2).unwrap();
            for (actual, expected) in report.values.iter().zip([1.0, 1.0, 2.0]) {
                assert!((actual / scale - expected).abs() < 1e-9);
            }
        }
    }
}
