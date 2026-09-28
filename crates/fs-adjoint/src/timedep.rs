//! Time-dependent adjoints for backward-Euler heat conduction:
//! (M + h K) u_{i+1} = M u_i + h b_i.
//!
//! Initial-state and source derivatives use transposed linear solves, not
//! differentiated Krylov iterations. The source adjoint does not need primal
//! checkpoints: its linear step Jacobian is independent of the temperature.
//! The legacy initial-state interface retains its treeverse recompute counter.

use fs_ad::revolve::{checkpointed_adjoint, min_budget};
use fs_solver::{CgState, CsrOp};
use fs_sparse::Csr;
use fs_sparse::precond::IdentityPrecond;

/// Backward-Euler heat problem on caller-assembled symmetric matrices.
///
/// The caller must supply positive-definite mass and an operator M + h K
/// suitable for conjugate gradients. Dimension and finite-value checks do
/// not prove symmetry or positive definiteness. This is an assembled discrete
/// model, not a conjugate-heat-transfer assembler or a physical certificate.
pub struct HeatAdjoint {
    sys: CsrOp,
    mass: Csr,
    assembled_h: f64,
    /// Time step used to assemble the operator. Changing this field is refused;
    /// construct a new problem to change the time step.
    pub h: f64,
    /// Step count. A zero-step trajectory is admitted.
    pub steps: usize,
}

fn validate_vector(values: &[f64], n: usize, name: &str) {
    assert_eq!(values.len(), n, "{name}: dimension mismatch");
    assert!(
        values.iter().all(|v| v.is_finite()),
        "{name}: nonfinite value"
    );
}

fn validate_matrix(matrix: &Csr, n: usize, name: &str) {
    assert_eq!(matrix.nrows(), n, "{name}: row dimension mismatch");
    assert_eq!(matrix.ncols(), n, "{name}: column dimension mismatch");
    for r in 0..n {
        assert!(
            matrix.row(r).1.iter().all(|v| v.is_finite()),
            "{name}: nonfinite entry in row {r}"
        );
    }
}

impl HeatAdjoint {
    /// Build from mass and stiffness (interior-reduced, symmetric).
    ///
    /// Panics on malformed dimensions, nonfinite entries or a nonpositive or
    /// nonfinite step. Symmetry and positive definiteness remain caller duties.
    #[must_use]
    pub fn new(mass: Csr, stiffness: &Csr, h: f64, steps: usize) -> HeatAdjoint {
        let n = mass.nrows();
        assert!(
            n > 0,
            "heat problem must have at least one degree of freedom"
        );
        assert!(
            h.is_finite() && h > 0.0,
            "heat time step must be finite and positive"
        );
        validate_matrix(&mass, n, "mass");
        validate_matrix(stiffness, n, "stiffness");
        let mut coo = fs_sparse::Coo::new(n, n);
        for r in 0..n {
            let (cols, vals) = mass.row(r);
            for (&c, &v) in cols.iter().zip(vals) {
                coo.push(r, c, v);
            }
            let (cols, vals) = stiffness.row(r);
            for (&c, &v) in cols.iter().zip(vals) {
                coo.push(r, c, h * v);
            }
        }
        let system = coo.assemble();
        validate_matrix(&system, n, "assembled heat operator");
        HeatAdjoint {
            sys: CsrOp::symmetric(system),
            mass,
            assembled_h: h,
            h,
            steps,
        }
    }

    fn validate_initial(&self, u0: &[f64]) {
        assert_eq!(
            self.h.to_bits(),
            self.assembled_h.to_bits(),
            "heat time step changed after assembly; construct a new HeatAdjoint"
        );
        validate_vector(u0, self.mass.nrows(), "initial temperature");
    }

    fn validate_sources(&self, sources: &[Vec<f64>]) {
        assert_eq!(
            sources.len(),
            self.steps,
            "one source vector is required per step"
        );
        for source in sources {
            validate_vector(source, self.mass.nrows(), "heat source");
        }
    }

    fn solve(&self, rhs: &[f64]) -> Vec<f64> {
        validate_vector(rhs, self.mass.nrows(), "heat solve right-hand side");
        let mut state = CgState::new(&self.sys, &IdentityPrecond, rhs);
        let report = state.run(&self.sys, &IdentityPrecond, 1e-13, 10_000);
        assert!(report.converged, "heat solve failed: {report:?}");
        validate_vector(&state.x, self.mass.nrows(), "heat solve result");
        state.x
    }

    fn step_forward_source(&self, u: &[f64], source: Option<&[f64]>) -> Vec<f64> {
        let mut rhs = vec![0.0; self.mass.nrows()];
        self.mass.spmv(u, &mut rhs);
        if let Some(source) = source {
            for (value, &load) in rhs.iter_mut().zip(source) {
                *value += self.h * load;
            }
        }
        self.solve(&rhs)
    }

    /// One unforced step: solve (M + h K) u⁺ = M u.
    fn step_forward(&self, u: &[f64]) -> Vec<f64> {
        self.step_forward_source(u, None)
    }

    /// Pull a state cotangent through Mᵀ (M + h K)⁻ᵀ. Both supplied
    /// matrices are symmetric under the constructor's mathematical contract.
    fn step_reverse(&self, lambda: &[f64]) -> Vec<f64> {
        let multiplier = self.solve(lambda);
        let mut out = vec![0.0; self.mass.nrows()];
        self.mass.spmv(&multiplier, &mut out);
        validate_vector(&out, self.mass.nrows(), "initial-state cotangent");
        out
    }

    /// Run the unforced forward sweep to the terminal state.
    #[must_use]
    pub fn forward(&self, u0: &[f64]) -> Vec<f64> {
        self.validate_initial(u0);
        let mut u = u0.to_vec();
        for _ in 0..self.steps {
            u = self.step_forward(&u);
        }
        u
    }

    /// Run a prescribed load schedule. `sources[i]` is the assembled load
    /// vector b_i for interval i, BEFORE multiplication by h. It is not a
    /// nodal source density: apply any required mass weighting before calling.
    /// The complete schedule is checked before the first solve.
    #[must_use]
    pub fn forward_with_sources(&self, u0: &[f64], sources: &[Vec<f64>]) -> Vec<f64> {
        self.validate_initial(u0);
        self.validate_sources(sources);
        let mut u = u0.to_vec();
        for source in sources {
            u = self.step_forward_source(&u, Some(source));
        }
        u
    }
}

fn terminal_misfit(terminal: &[f64], target: &[f64]) -> (f64, Vec<f64>) {
    validate_vector(target, terminal.len(), "target temperature");
    let seed: Vec<f64> = terminal.iter().zip(target).map(|(a, b)| a - b).collect();
    validate_vector(&seed, terminal.len(), "terminal cotangent");
    let objective = seed.iter().map(|v| 0.5 * v * v).sum::<f64>();
    assert!(objective.is_finite(), "terminal objective overflow");
    (objective, seed)
}

/// Gradient of J = ½‖u_N − target‖² with respect to u₀.
///
/// The second return value is the legacy treeverse RECOMPUTATION count; it
/// excludes the `problem.steps` forward solves used to obtain the terminal
/// seed. No primal trajectory is stored in full.
#[must_use]
pub fn heat_initial_gradient(problem: &HeatAdjoint, u0: &[f64], target: &[f64]) -> (Vec<f64>, u64) {
    problem.validate_initial(u0);
    validate_vector(target, problem.mass.nrows(), "target temperature");
    let forward = |_i: usize, u: &Vec<f64>| -> Vec<f64> { problem.step_forward(u) };
    let reverse =
        |_i: usize, _u: &Vec<f64>, bar: Vec<f64>| -> Vec<f64> { problem.step_reverse(&bar) };
    let terminal = problem.forward(u0);
    let (_, seed) = terminal_misfit(&terminal, target);
    let (bar, stats) = checkpointed_adjoint(
        &u0.to_vec(),
        problem.steps,
        min_budget(problem.steps),
        &forward,
        &reverse,
        seed,
    );
    (bar, stats.forward_steps)
}

/// Terminal-misfit derivatives for a prescribed heat-source schedule.
#[derive(Debug, Clone, PartialEq)]
pub struct HeatSourceGradient {
    /// Value of ½‖u_N − target‖².
    pub objective: f64,
    /// Terminal temperatures from the same forward solve as the gradient.
    pub terminal: Vec<f64>,
    /// Derivative with respect to each initial temperature.
    pub initial: Vec<f64>,
    /// `sources[i][j]` is dJ/db_i[j], in chronological order. Includes h.
    pub sources: Vec<Vec<f64>>,
    /// Actual forward linear solves, including terminal-seed construction.
    pub forward_evaluations: u64,
}

/// Differentiate the actual forced transient solve with respect to u₀ and
/// every b_i. For z_i = (M + h K)⁻ᵀ λ_{i+1}, the pullbacks are
/// λ_i = Mᵀ z_i and dJ/db_i = h z_i.
///
/// Requires N forward and N transposed solves, independent of the number of
/// source controls. Live primal storage is O(n); the requested gradient output
/// itself is O(N n). No primal recomputation or finite differencing is needed.
#[must_use]
pub fn heat_source_gradient(
    problem: &HeatAdjoint,
    u0: &[f64],
    target: &[f64],
    sources: &[Vec<f64>],
) -> HeatSourceGradient {
    validate_vector(target, problem.mass.nrows(), "target temperature");
    let terminal = problem.forward_with_sources(u0, sources);
    let (objective, mut bar) = terminal_misfit(&terminal, target);
    let n = problem.mass.nrows();
    let mut gradients = vec![Vec::new(); problem.steps];
    for i in (0..problem.steps).rev() {
        let multiplier = problem.solve(&bar);
        gradients[i] = multiplier.iter().map(|value| problem.h * value).collect();
        validate_vector(&gradients[i], n, "source cotangent");
        problem.mass.spmv(&multiplier, &mut bar);
        validate_vector(&bar, n, "initial-state cotangent");
    }
    HeatSourceGradient {
        objective,
        terminal,
        initial: bar,
        sources: gradients,
        forward_evaluations: problem.steps as u64,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix(a: [[f64; 2]; 2]) -> Csr {
        let mut coo = fs_sparse::Coo::new(2, 2);
        for (i, row) in a.iter().enumerate() {
            for (j, &value) in row.iter().enumerate() {
                coo.push(i, j, value);
            }
        }
        coo.assemble()
    }

    fn problem(steps: usize) -> HeatAdjoint {
        // SPD but NONCOMMUTING matrices catch reversed pullback ordering.
        HeatAdjoint::new(
            matrix([[2.0, 0.3], [0.3, 1.0]]),
            &matrix([[3.0, -0.8], [-0.8, 2.0]]),
            0.17,
            steps,
        )
    }

    fn loss(problem: &HeatAdjoint, u0: &[f64], sources: &[Vec<f64>]) -> f64 {
        let terminal = problem.forward_with_sources(u0, sources);
        terminal_misfit(&terminal, &[0.2, -0.1]).0
    }

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 2e-8 * expected.abs().max(1.0),
            "{actual} != {expected}"
        );
    }

    #[test]
    fn scheduled_source_and_initial_gradients_match_resolve_differences() {
        let heat = problem(4);
        let u0 = vec![1.2, -0.4];
        let sources = vec![
            vec![0.7, 0.2],
            vec![-0.3, 0.9],
            vec![1.1, -0.2],
            vec![0.4, 0.6],
        ];
        let result = heat_source_gradient(&heat, &u0, &[0.2, -0.1], &sources);
        let eps = 1e-5;
        close(result.objective, loss(&heat, &u0, &sources));
        for j in 0..2 {
            let mut plus = u0.clone();
            let mut minus = u0.clone();
            plus[j] += eps;
            minus[j] -= eps;
            close(
                result.initial[j],
                (loss(&heat, &plus, &sources) - loss(&heat, &minus, &sources)) / (2.0 * eps),
            );
        }
        for i in 0..heat.steps {
            for j in 0..2 {
                let mut plus = sources.clone();
                let mut minus = sources.clone();
                plus[i][j] += eps;
                minus[i][j] -= eps;
                close(
                    result.sources[i][j],
                    (loss(&heat, &u0, &plus) - loss(&heat, &u0, &minus)) / (2.0 * eps),
                );
            }
        }
        assert_eq!(result.forward_evaluations, 4);
        assert_eq!(
            result,
            heat_source_gradient(&heat, &u0, &[0.2, -0.1], &sources)
        );
    }

    #[test]
    fn scalar_forcing_matches_backward_euler_closed_form() {
        let mut mass = fs_sparse::Coo::new(1, 1);
        let mut stiffness = fs_sparse::Coo::new(1, 1);
        mass.push(0, 0, 2.0);
        stiffness.push(0, 0, 3.0);
        let heat = HeatAdjoint::new(mass.assemble(), &stiffness.assemble(), 0.25, 3);
        let sources = vec![vec![4.0]; 3];
        let result = heat_source_gradient(&heat, &[1.0], &[0.5], &sources);
        let a: f64 = 2.0 / 2.75;
        let terminal = a.powi(3) + (1.0 / 2.75) * (1.0 + a + a * a);
        close(result.terminal[0], terminal);
        close(result.initial[0], (terminal - 0.5) * a.powi(3));
        for i in 0..3 {
            close(
                result.sources[i][0],
                (terminal - 0.5) * (0.25 / 2.75) * a.powi((2 - i) as i32),
            );
        }
    }

    #[test]
    fn zero_schedule_matches_legacy_unforced_gradient() {
        let heat = problem(7);
        let u0 = [0.8, -0.5];
        let target = [0.2, -0.1];
        let sources = vec![vec![0.0; 2]; 7];
        let result = heat_source_gradient(&heat, &u0, &target, &sources);
        assert_eq!(result.terminal, heat.forward(&u0));
        assert_eq!(result.initial, heat_initial_gradient(&heat, &u0, &target).0);
    }

    #[test]
    fn empty_horizon_needs_no_sources_or_solves() {
        let result = heat_source_gradient(&problem(0), &[0.8, -0.5], &[0.2, -0.1], &[]);
        assert_eq!(result.terminal, vec![0.8, -0.5]);
        close(result.initial[0], 0.6);
        close(result.initial[1], -0.4);
        assert!(result.sources.is_empty());
        assert_eq!(result.forward_evaluations, 0);
    }

    #[test]
    #[should_panic(expected = "target temperature: dimension mismatch")]
    fn refuses_truncated_target_before_solving() {
        heat_initial_gradient(&problem(2), &[1.0, 2.0], &[0.0]);
    }

    #[test]
    #[should_panic(expected = "one source vector is required per step")]
    fn refuses_incomplete_schedule() {
        problem(2).forward_with_sources(&[1.0, 2.0], &[vec![0.0; 2]]);
    }

    #[test]
    #[should_panic(expected = "heat source: dimension mismatch")]
    fn refuses_truncated_source() {
        problem(1).forward_with_sources(&[1.0, 2.0], &[vec![0.0]]);
    }

    #[test]
    #[should_panic(expected = "heat source: nonfinite value")]
    fn refuses_nonfinite_source() {
        problem(1).forward_with_sources(&[1.0, 2.0], &[vec![0.0, f64::NAN]]);
    }

    #[test]
    #[should_panic(expected = "heat time step changed after assembly")]
    fn refuses_stale_time_step_operator() {
        let mut heat = problem(1);
        heat.h = 0.5;
        heat.forward_with_sources(&[1.0, 2.0], &[vec![0.0; 2]]);
    }
}
