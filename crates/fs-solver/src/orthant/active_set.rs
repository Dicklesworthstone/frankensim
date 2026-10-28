//! Deterministic active-set minimization of a strictly convex orthant quadratic.
//!
//! Solve `min 0.5*x^T*H*x + b^T*x` with `x >= 0`. The supplied row-major
//! Hessian must admit positive normalized Cholesky pivots. Singular PSD and
//! indefinite systems refuse; no regularization or arbitrary null-space choice
//! converts them into a uniqueness claim. The separate Gram-coordinate solver
//! remains the PSD baseline. This numerical owner knows nothing about contact.

use fs_la::factor::{FactorError, cholesky};

/// Bounded dense-solve work and acceptance tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpdOrthantConfig {
    /// Maximum dimension, in `1..=64`; empty problems are also admitted.
    pub max_variables: usize,
    /// Maximum principal-subproblem solves, in `1..=4096`.
    pub max_iterations: usize,
    /// Positive finite absolute KKT tolerance in the original gradient units.
    pub gradient_tolerance: f64,
    /// Positive symmetry tolerance after diagonal normalization, at most `1e-6`.
    pub symmetry_tolerance: f64,
    /// Minimum normalized squared Cholesky pivot, in `[64*EPSILON, 1e-4]`.
    pub relative_pivot_tolerance: f64,
}

impl Default for SpdOrthantConfig {
    fn default() -> Self {
        Self {
            max_variables: 64,
            max_iterations: 1024,
            gradient_tolerance: 1e-10,
            symmetry_tolerance: 1e-12,
            relative_pivot_tolerance: 1e-12,
        }
    }
}

/// A nonnegative point with independently reconstructed original-system KKT
/// diagnostics. The SPD admission is a numerical test, not an interval proof.
#[derive(Clone, Debug, PartialEq)]
pub struct SpdOrthantSolution {
    /// Accepted primal coordinates.
    pub point: Vec<f64>,
    /// Original supplied `H*x+b`, recomputed independently of the factor solve.
    pub gradient: Vec<f64>,
    /// Quadratic objective, evaluated as `x dot (gradient+b)/2`.
    pub objective: f64,
    /// Maximum positive-coordinate stationarity or zero-coordinate dual error.
    pub residual: f64,
    /// Maximum negative primal coordinate magnitude.
    pub primal_residual: f64,
    /// Maximum negative gradient magnitude.
    pub dual_residual: f64,
    /// Maximum absolute `point[i]*gradient[i]`, in objective units.
    pub complementarity_residual: f64,
    /// Completed principal-subproblem solves.
    pub iterations: usize,
    /// Minimum squared Cholesky pivot after normalization; one for an empty problem.
    pub min_normalized_pivot: f64,
}

/// Refusal without a partially accepted point or an infeasibility claim.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SpdOrthantError {
    /// Invalid dimensions, limits, nonfinite input, or symmetry.
    Invalid(&'static str),
    /// Finite input led to unrepresentable arithmetic.
    NonFinite,
    /// A nonpositive normalized pivot; singular and indefinite are not conflated
    /// with a proof that the constrained quadratic is infeasible or unbounded.
    NotPositiveDefinite {
        /// Pivot index.
        index: usize,
    },
    /// A positive pivot is below the admitted numerical independence threshold.
    Degenerate {
        /// Pivot index.
        index: usize,
        /// Squared normalized Cholesky pivot.
        pivot: f64,
    },
    /// The active-set solve budget was exhausted.
    IterationLimit {
        /// Principal subproblems solved.
        iterations: usize,
        /// Last independently reconstructed KKT residual.
        residual: f64,
    },
    /// No entering variable remains, but original-system stationarity does not
    /// meet the requested tolerance; a finite solve is not called converged.
    ResidualNotSatisfied {
        /// Independently reconstructed KKT residual.
        residual: f64,
    },
    /// Caller cancellation, including the final publication checkpoint.
    Cancelled,
}

impl core::fmt::Display for SpdOrthantError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "SPD orthant solve: {self:?}")
    }
}
impl std::error::Error for SpdOrthantError {}

/// Lawson--Hanson primal active-set solve using the existing dense Cholesky
/// owner. The most negative scaled inactive gradient enters; ties retain the
/// lowest index. A principal solve that crosses a bound takes the first
/// feasible line-search step and releases its blocking coordinate.
///
/// Diagonal equilibration changes variables, not the physical quadratic.
/// Admitted small asymmetry is averaged for factorization; publication still
/// requires the KKT residual of the original supplied matrix. No diagonal
/// shift or truncation of negative multipliers hides a failed solve.
/// Cancellation is checked at rows, before/after each at-most-64-dimensional
/// factorization and solve, and immediately before publication. All state is
/// private scratch. This bounded candidate does not provide an adjoint.
pub fn solve_spd_orthant(
    hessian: &[f64],
    linear: &[f64],
    config: SpdOrthantConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SpdOrthantSolution, SpdOrthantError> {
    let problem = Problem::new(hessian, linear, config, &mut cancelled)?;
    let n = linear.len();
    let mut scaled_point = vec![0.0; n];
    let mut point = vec![0.0; n];
    let mut gradient = linear.to_vec();
    let mut passive = vec![false; n];
    let mut iterations = 0;
    loop {
        let residual =
            problem.evaluate(&scaled_point, &mut point, &mut gradient, &mut cancelled)?;
        if residual <= config.gradient_tolerance {
            return problem.finish(point, gradient, residual, iterations, &mut cancelled);
        }
        let mut entering = None;
        let mut most_negative = f64::INFINITY;
        for i in 0..n {
            if !passive[i] && gradient[i] < -config.gradient_tolerance {
                let scaled_gradient = finite(gradient[i] / problem.scales[i])?;
                if scaled_gradient < most_negative {
                    entering = Some(i);
                    most_negative = scaled_gradient;
                }
            }
        }
        let Some(entering) = entering else {
            return Err(SpdOrthantError::ResidualNotSatisfied { residual });
        };
        passive[entering] = true;
        loop {
            if iterations == config.max_iterations {
                let residual =
                    problem.evaluate(&scaled_point, &mut point, &mut gradient, &mut cancelled)?;
                return Err(SpdOrthantError::IterationLimit {
                    iterations,
                    residual,
                });
            }
            let candidate = problem.solve_active(&passive, &mut cancelled)?;
            iterations += 1;
            let mut fraction = 1.0;
            let mut blocker = None;
            for i in 0..n {
                if passive[i] && candidate[i] < 0.0 {
                    let scale = scaled_point[i].max(-candidate[i]);
                    let numerator = scaled_point[i] / scale;
                    let ratio = numerator / (numerator + (-candidate[i] / scale));
                    if blocker.is_none() || ratio < fraction {
                        fraction = ratio;
                        blocker = Some(i);
                    }
                }
            }
            let Some(blocker) = blocker else {
                scaled_point = candidate;
                break;
            };
            for i in 0..n {
                if passive[i] {
                    // This is the feasible bound step, not clipping a solved
                    // multiplier. Roundoff at a reached bound is set to zero.
                    scaled_point[i] =
                        finite((1.0 - fraction) * scaled_point[i] + fraction * candidate[i])?
                            .max(0.0);
                    if scaled_point[i] == 0.0 {
                        passive[i] = false;
                    }
                }
            }
            scaled_point[blocker] = 0.0;
            passive[blocker] = false;
        }
    }
}

struct Problem<'a> {
    raw: &'a [f64],
    linear: &'a [f64],
    matrix: Vec<f64>,
    scaled_linear: Vec<f64>,
    scales: Vec<f64>,
    min_pivot: f64,
}

impl<'a> Problem<'a> {
    fn new(
        raw: &'a [f64],
        linear: &'a [f64],
        config: SpdOrthantConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, SpdOrthantError> {
        poll(cancelled)?;
        if !(1..=64).contains(&config.max_variables)
            || !(1..=4096).contains(&config.max_iterations)
            || !config.gradient_tolerance.is_finite()
            || config.gradient_tolerance <= 0.0
            || !config.symmetry_tolerance.is_finite()
            || config.symmetry_tolerance <= 0.0
            || config.symmetry_tolerance > 1e-6
            || !(64.0 * f64::EPSILON..=1e-4).contains(&config.relative_pivot_tolerance)
        {
            return Err(SpdOrthantError::Invalid("configuration"));
        }
        let n = linear.len();
        if n > config.max_variables || raw.len() != n * n || linear.iter().any(|v| !v.is_finite()) {
            return Err(SpdOrthantError::Invalid("dimensions or linear term"));
        }
        let mut scales = Vec::with_capacity(n);
        let mut scaled_linear = Vec::with_capacity(n);
        for i in 0..n {
            poll(cancelled)?;
            if raw[i * n..(i + 1) * n].iter().any(|v| !v.is_finite()) {
                return Err(SpdOrthantError::Invalid("nonfinite Hessian"));
            }
            if raw[i * n + i] <= 0.0 {
                return Err(SpdOrthantError::NotPositiveDefinite { index: i });
            }
            let scale = raw[i * n + i].sqrt();
            scales.push(scale);
            scaled_linear.push(finite(linear[i] / scale)?);
        }
        let mut matrix = vec![0.0; n * n];
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..=i {
                let a = finite(raw[i * n + j] / scales[i] / scales[j])?;
                let b = finite(raw[j * n + i] / scales[j] / scales[i])?;
                if (a - b).abs() > config.symmetry_tolerance {
                    return Err(SpdOrthantError::Invalid("Hessian symmetry"));
                }
                let entry = 0.5 * a + 0.5 * b;
                matrix[i * n + j] = entry;
                matrix[j * n + i] = entry;
            }
        }
        poll(cancelled)?;
        let factor = cholesky(&matrix, n).map_err(factor_error)?;
        poll(cancelled)?;
        let mut min_pivot = if n == 0 { 1.0 } else { f64::MAX };
        for i in 0..n {
            let pivot = finite(factor.l(i, i) * factor.l(i, i))?;
            if pivot < config.relative_pivot_tolerance {
                return Err(SpdOrthantError::Degenerate { index: i, pivot });
            }
            min_pivot = min_pivot.min(pivot);
        }
        Ok(Self {
            raw,
            linear,
            matrix,
            scaled_linear,
            scales,
            min_pivot,
        })
    }

    fn solve_active(
        &self,
        passive: &[bool],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, SpdOrthantError> {
        poll(cancelled)?;
        let n = self.linear.len();
        let indices: Vec<_> = (0..n).filter(|&i| passive[i]).collect();
        let m = indices.len();
        let mut matrix = Vec::with_capacity(m * m);
        let mut rhs = Vec::with_capacity(m);
        for &i in &indices {
            poll(cancelled)?;
            for &j in &indices {
                matrix.push(self.matrix[i * n + j]);
            }
            rhs.push(-self.scaled_linear[i]);
        }
        let factor = cholesky(&matrix, m).map_err(factor_error)?;
        poll(cancelled)?;
        factor.solve(&mut rhs);
        poll(cancelled)?;
        let mut candidate = vec![0.0; n];
        for (&i, value) in indices.iter().zip(rhs) {
            candidate[i] = finite(value)?;
        }
        Ok(candidate)
    }

    fn evaluate(
        &self,
        scaled: &[f64],
        point: &mut [f64],
        gradient: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, SpdOrthantError> {
        let n = self.linear.len();
        for i in 0..n {
            point[i] = finite(scaled[i] / self.scales[i])?;
        }
        let mut residual = 0.0_f64;
        for i in 0..n {
            poll(cancelled)?;
            gradient[i] =
                finite(crate::dot(&self.raw[i * n..(i + 1) * n], point) + self.linear[i])?;
            residual = residual.max(if point[i] > 0.0 {
                gradient[i].abs()
            } else {
                (-gradient[i]).max(0.0)
            });
        }
        Ok(residual)
    }

    fn finish(
        &self,
        point: Vec<f64>,
        gradient: Vec<f64>,
        residual: f64,
        iterations: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpdOrthantSolution, SpdOrthantError> {
        let mut primal_residual = 0.0_f64;
        let mut dual_residual = 0.0_f64;
        let mut complementarity_residual = 0.0_f64;
        let mut objective_terms = Vec::with_capacity(point.len());
        for i in 0..point.len() {
            primal_residual = primal_residual.max(-point[i]);
            dual_residual = dual_residual.max(-gradient[i]);
            complementarity_residual =
                complementarity_residual.max(finite(point[i] * gradient[i])?.abs());
            objective_terms.push(0.5 * gradient[i] + 0.5 * self.linear[i]);
        }
        let objective = finite(crate::dot(&point, &objective_terms))?;
        poll(cancelled)?;
        Ok(SpdOrthantSolution {
            point,
            gradient,
            objective,
            residual,
            primal_residual,
            dual_residual,
            complementarity_residual,
            iterations,
            min_normalized_pivot: self.min_pivot,
        })
    }
}

fn finite(value: f64) -> Result<f64, SpdOrthantError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SpdOrthantError::NonFinite)
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SpdOrthantError> {
    if cancelled() {
        Err(SpdOrthantError::Cancelled)
    } else {
        Ok(())
    }
}
fn factor_error(error: FactorError) -> SpdOrthantError {
    match error {
        FactorError::NotSpd { index } | FactorError::Singular { index } => {
            SpdOrthantError::NotPositiveDefinite { index }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config() -> SpdOrthantConfig {
        SpdOrthantConfig {
            gradient_tolerance: 1e-12,
            ..SpdOrthantConfig::default()
        }
    }

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
            "{actual} != {expected}"
        );
    }

    fn release_problem() -> ([f64; 9], [f64; 3]) {
        // The first entering coordinate has a negative unconstrained value
        // once all three are considered. Its bound must be released, not
        // repaired by clipping the final unconstrained solution.
        (
            [1.0, 0.625, 0.625, 0.625, 1.0, 0.0, 0.625, 0.0, 1.0],
            [-1.125, -1.0, -1.0],
        )
    }

    #[test]
    fn coupled_quadratic_matches_independent_two_contact_solution() {
        let h = [2.5, -1.5, -1.5, 2.5];
        let b = [-1.8, -1.2];
        let solution = solve_spd_orthant(&h, &b, config(), || false).unwrap();
        close(solution.point[0], 1.575);
        close(solution.point[1], 1.425);
        close(solution.objective, -2.2725);
        for i in 0..2 {
            let gradient = h[2 * i] * solution.point[0] + h[2 * i + 1] * solution.point[1] + b[i];
            close(solution.gradient[i], gradient);
        }
        assert!(solution.residual <= config().gradient_tolerance);
        assert!(solution.primal_residual == 0.0);
        assert!(solution.dual_residual < 1e-12);
        assert!(solution.complementarity_residual < 1e-12);
        assert_eq!(solution.iterations, 2);
    }

    #[test]
    fn bound_crossing_releases_a_previously_active_coordinate() {
        let (h, b) = release_problem();
        let solution = solve_spd_orthant(&h, &b, config(), || false).unwrap();
        for (actual, expected) in solution.point.iter().zip([0.0, 1.0, 1.0]) {
            close(*actual, expected);
        }
        for (actual, expected) in solution.gradient.iter().zip([0.125, 0.0, 0.0]) {
            close(*actual, expected);
        }
        close(solution.objective, -1.0);
        assert_eq!(solution.iterations, 4);
        assert!(matches!(
            solve_spd_orthant(
                &h,
                &b,
                SpdOrthantConfig {
                    max_iterations: 3,
                    ..config()
                },
                || false
            ),
            Err(SpdOrthantError::IterationLimit { iterations: 3, .. })
        ));
    }

    #[test]
    fn permutation_and_diagonal_unit_scaling_preserve_the_minimizer() {
        let (h, b) = release_problem();
        let baseline = solve_spd_orthant(&h, &b, config(), || false).unwrap();
        for order in [[2, 0, 1], [1, 2, 0]] {
            let permuted_h: Vec<_> = order
                .iter()
                .flat_map(|&i| order.iter().map(move |&j| h[i * 3 + j]))
                .collect();
            let permuted_b = order.map(|i| b[i]);
            let solution = solve_spd_orthant(&permuted_h, &permuted_b, config(), || false).unwrap();
            for i in 0..3 {
                close(solution.point[i], baseline.point[order[i]]);
            }
        }
        let scales = [2.0_f64.powi(-20), 4.0, 2.0_f64.powi(20)];
        let scaled_h: Vec<_> = (0..3)
            .flat_map(|i| (0..3).map(move |j| h[i * 3 + j] * scales[i] * scales[j]))
            .collect();
        let scaled_b: Vec<_> = (0..3).map(|i| b[i] * scales[i]).collect();
        let solution = solve_spd_orthant(&scaled_h, &scaled_b, config(), || false).unwrap();
        for i in 0..3 {
            close(solution.point[i] * scales[i], baseline.point[i]);
        }
        close(solution.objective, baseline.objective);
        let tiny_units = solve_spd_orthant(&[1e-308], &[-1.0], config(), || false).unwrap();
        assert!((tiny_units.point[0] / 1e308 - 1.0).abs() < 2e-14);
        assert!(tiny_units.min_normalized_pivot > 0.99);
    }

    #[test]
    fn empty_inactive_and_maximum_dimension_problems_are_supported() {
        let empty = solve_spd_orthant(&[], &[], config(), || false).unwrap();
        assert!(empty.point.is_empty() && empty.gradient.is_empty());
        assert_eq!(empty.iterations, 0);
        let inactive =
            solve_spd_orthant(&[2.0, -1.0, -1.0, 2.0], &[1.0, 0.0], config(), || false).unwrap();
        assert_eq!(inactive.point, [0.0, 0.0]);
        assert_eq!(inactive.iterations, 0);
        let mut identity = vec![0.0; 64 * 64];
        for i in 0..64 {
            identity[i * 64 + i] = 1.0;
        }
        let maximum = solve_spd_orthant(&identity, &[-1.0; 64], config(), || false).unwrap();
        assert_eq!(maximum.point, [1.0; 64]);
        assert_eq!(maximum.iterations, 64);
        assert!(matches!(
            solve_spd_orthant(&vec![0.0; 65 * 65], &[0.0; 65], config(), || false),
            Err(SpdOrthantError::Invalid(_))
        ));
    }

    #[test]
    fn singular_indefinite_and_near_dependent_hessians_refuse_without_regularization() {
        for h in [
            [1.0, 1.0, 1.0, 1.0],
            [1.0, 2.0, 2.0, 1.0],
            [1.0, 0.0, 0.0, 0.0],
        ] {
            assert!(matches!(
                solve_spd_orthant(&h, &[-1.0, -1.0], config(), || false),
                Err(SpdOrthantError::NotPositiveDefinite { .. })
            ));
            // A zero optimizer does not waive the advertised SPD admission.
            assert!(solve_spd_orthant(&h, &[1.0, 1.0], config(), || false).is_err());
        }
        let almost_parallel = 1.0 - 4.0 * f64::EPSILON;
        assert!(matches!(
            solve_spd_orthant(
                &[1.0, almost_parallel, almost_parallel, 1.0],
                &[-1.0, -1.0],
                config(),
                || false
            ),
            Err(SpdOrthantError::Degenerate { .. })
        ));
    }

    #[test]
    fn malformed_inputs_unrepresentable_points_and_false_residual_claims_refuse() {
        for bad in [
            SpdOrthantConfig {
                gradient_tolerance: 0.0,
                ..config()
            },
            SpdOrthantConfig {
                gradient_tolerance: f64::NAN,
                ..config()
            },
            SpdOrthantConfig {
                symmetry_tolerance: 0.0,
                ..config()
            },
            SpdOrthantConfig {
                relative_pivot_tolerance: 1e-30,
                ..config()
            },
            SpdOrthantConfig {
                max_iterations: 0,
                ..config()
            },
        ] {
            assert!(matches!(
                solve_spd_orthant(&[1.0], &[-1.0], bad, || false),
                Err(SpdOrthantError::Invalid(_))
            ));
        }
        assert!(solve_spd_orthant(&[1.0], &[], config(), || false).is_err());
        assert!(solve_spd_orthant(&[f64::NAN], &[-1.0], config(), || false).is_err());
        assert!(solve_spd_orthant(&[1.0], &[f64::INFINITY], config(), || false).is_err());
        assert!(
            solve_spd_orthant(&[1.0, 0.0, 0.1, 1.0], &[-1.0, -1.0], config(), || false).is_err()
        );
        assert_eq!(
            solve_spd_orthant(&[1e-308], &[-2.0], config(), || false),
            Err(SpdOrthantError::NonFinite)
        );
        // Tiny admitted asymmetry is not allowed to disappear from the
        // original-system residual when the requested accuracy resolves it.
        assert!(matches!(
            solve_spd_orthant(
                &[1.0, 0.5 + 1e-13, 0.5, 1.0],
                &[-1.0, -1.0],
                SpdOrthantConfig {
                    gradient_tolerance: 1e-16,
                    ..config()
                },
                || false
            ),
            Err(SpdOrthantError::ResidualNotSatisfied { .. })
        ));
    }

    #[test]
    fn cancellation_is_latched_at_every_boundary_and_retry_replays() {
        let (h, b) = release_problem();
        let mut calls = 0;
        let expected = solve_spd_orthant(&h, &b, config(), || {
            calls += 1;
            false
        })
        .unwrap();
        for stop in 1..=calls {
            let mut observed = 0;
            let result = solve_spd_orthant(&h, &b, config(), || {
                observed += 1;
                observed == stop
            });
            assert_eq!(result, Err(SpdOrthantError::Cancelled));
            assert_eq!(observed, stop);
        }
        assert_eq!(
            solve_spd_orthant(&h, &b, config(), || false).unwrap(),
            expected
        );
    }
}
