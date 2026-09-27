//! Bound-controlled corrections of the complete stored solid/reference system.
//! FGMRES proposes fields for A-B*C; only the existing outward coupled maximum
//! enclosure admits them. A frozen-solid residual cannot terminate this solve.

use std::cell::RefCell;

use fs_solver::{FgmresState, LinearOp};
use fs_sparse::Csr;

use super::{ConductionError, Cx, LinearRobinFeedbackAnalyzer, LinearRobinMaximumAnalysis,
    invalid, poll, zeros};
use super::super::{LinearGoalSolve, LinearGoalSolveConfig, LinearGoalStop};

impl LinearRobinFeedbackAnalyzer<'_> {
    /// Area-weighted wall temperatures in this analyzer's bound port order.
    /// Used to rebuild air states after adopting a corrected nodal field.
    /// The field must preserve all prescribed values and material support.
    /// This evaluation does not claim energy balance or solve convergence.
    ///
    /// # Errors
    /// Inadmissible field, nonfinite arithmetic/allocation failure or cancellation.
    pub fn wall_mean_temperatures(
        &self, cx: &Cx<'_>, temperature: &[f64],
    ) -> Result<Vec<f64>, ConductionError> {
        super::super::validate_field(cx, self.solid.problem, self.solid.dofs(), temperature)?;
        let mut walls = zeros(cx, self.ports.len())?;
        let mut work = Poll { cx, visits: 0 };
        for (out, port) in walls.iter_mut().zip(&self.ports) {
            for (vertices, area) in &port.faces {
                work.tick()?;
                let weight = finite((area / 3.0) / port.area_m2)?;
                for &vertex in vertices {
                    work.tick()?;
                    *out = finite(weight.mul_add(temperature[vertex], *out))?;
                }
            }
        }
        poll(cx, 0)?;
        Ok(walls)
    }

    /// Correct a regional maximum with the complete affine Robin feedback.
    ///
    /// Solves `(A-B*C) delta = b+B*(d+C*T)-A*T` by the existing restarted
    /// FGMRES, not CG: the coupled operator need not be symmetric or a
    /// contractive fixed-point map. Neither the air/reference variables nor
    /// prescribed-node contributions are frozen. Prepared inverse and response
    /// columns are reused and independently checked by `analyze_maximum`.
    ///
    /// `check_every` caps each restart cycle (also capped by the free dimension
    /// and remaining iterations). Every restart/recomputed defect spends the
    /// same `max_primal_iterations` allowance. The inner residual only decides
    /// when to recompute a defect; it NEVER admits goal success. The outward
    /// error bound on the actual rounded full temperature field does that.
    /// A missing bound returns `BoundUnavailable` without speculative work.
    ///
    /// The result is a temperature candidate and its stored-system bound, NOT
    /// a `ConductionSolution` or an air/energy report. Consumers must recompute
    /// their physical wall/air/flux outputs before publishing a changed field.
    /// Nonlinear conductivity, radiation and moving geometry remain unsupported.
    ///
    /// # Errors
    /// Invalid controls, inadmissible fields, arithmetic/allocation failure or
    /// cancellation. Budget/stagnation outcomes retain the best checked field.
    pub fn solve_maximum_to_goal(
        &self,
        cx: &Cx<'_>,
        initial_temperature: &[f64],
        region_vertices: &[usize],
        config: LinearGoalSolveConfig,
    ) -> Result<LinearGoalSolve<LinearRobinMaximumAnalysis>, ConductionError> {
        self.solve_maximum_to_goal_observed(
            cx, initial_temperature, region_vertices, config, |_, _| {},
        )
    }

    /// Coupled maximum control with progress after each completed enclosure.
    /// The callback sees all checked candidates, including rejected ones.
    /// Cancellation is checked AFTER it, including a passing initial/candidate
    /// check. Neither this analyzer nor the caller's initial field is mutated.
    ///
    /// # Errors
    /// The same refusals as `solve_maximum_to_goal`.
    #[allow(clippy::too_many_lines)] // One loop owns the shared work and publication boundary.
    pub fn solve_maximum_to_goal_observed(
        &self,
        cx: &Cx<'_>,
        initial_temperature: &[f64],
        region_vertices: &[usize],
        config: LinearGoalSolveConfig,
        mut observe: impl FnMut(usize, &LinearRobinMaximumAnalysis),
    ) -> Result<LinearGoalSolve<LinearRobinMaximumAnalysis>, ConductionError> {
        poll(cx, 0)?;
        if !(config.absolute_tolerance.is_finite() && config.absolute_tolerance > 0.0)
            || !(1..=32).contains(&config.check_every)
        {
            return Err(invalid("coupled goal solve needs a finite positive tolerance and check_every in 1..=32"));
        }
        let analysis = self.analyze_maximum(cx, initial_temperature, region_vertices)?;
        observe(0, &analysis);
        poll(cx, 0)?;
        let mut result = LinearGoalSolve {
            temperature: initial_temperature.to_vec(), analysis,
            stop: LinearGoalStop::IterationBudget, primal_iterations: 0,
            defect_corrections: 0, goal_checks: 1,
        };
        let Some(mut best_bound) = result.analysis.algebraic_half_width_k() else {
            result.stop = LinearGoalStop::BoundUnavailable;
            poll(cx, 0)?;
            return Ok(result);
        };
        if best_bound <= config.absolute_tolerance {
            result.stop = LinearGoalStop::GoalTolerance;
            poll(cx, 0)?;
            return Ok(result);
        }
        if config.max_primal_iterations == 0 { poll(cx, 0)?; return Ok(result); }

        let operator = FeedbackOp {
            matrix: &self.solid.response.matrix,
            injection: &self.injection,
            feedback: &self.feedback,
            cx,
            ports: RefCell::new(zeros(cx, self.offset.len())?),
            fault: RefCell::new(None),
        };
        let preconditioner = crate::solve::spd_preconditioner(&self.solid.response.matrix);
        poll(cx, 0)?;
        'attempt: loop {
            poll(cx, result.primal_iterations)?;
            let base = self.solid.dofs().gather(&result.temperature);
            let mut defect = operator.residual(&self.solid.rhs, &self.offset, &base)?;
            let mut scale = 0.0_f64;
            for (i, &value) in defect.iter().enumerate() {
                if i % 512 == 0 { poll(cx, result.primal_iterations)?; }
                scale = scale.max(value.abs());
            }
            if scale == 0.0 {
                // A zero rounded defect is not a zero outward residual.
                result.stop = LinearGoalStop::NoProgress;
                break;
            }
            for (i, value) in defect.iter_mut().enumerate() {
                if i % 512 == 0 { poll(cx, result.primal_iterations)?; }
                *value /= scale;
            }
            let start_bound = best_bound;
            let mut state = FgmresState::new(&defect, config.check_every.min(operator.n()));
            loop {
                poll(cx, result.primal_iterations)?;
                let before = state.iters;
                let cycles_before = state.history.len();
                state.restart = config.check_every.min(operator.n()).min(
                    config.max_primal_iterations - result.primal_iterations,
                );
                // A recurrence cutoff, not an acceptance tolerance. Further
                // work can recompute the PHYSICAL defect under the same cap.
                state.run(&operator, &preconditioner, &defect, f64::EPSILON, 1);
                result.primal_iterations += state.iters - before;
                operator.check()?;
                poll(cx, result.primal_iterations)?;
                if !state.x.iter().all(|value| value.is_finite()) { break; }
                let mut candidate = self.solid.dofs().prescribed().to_vec();
                for (i, &vertex) in self.solid.dofs().free().iter().enumerate() {
                    if i % 512 == 0 { poll(cx, result.primal_iterations)?; }
                    candidate[vertex] = finite(scale.mul_add(state.x[i], base[i]))?;
                }
                let checked = self.analyze_maximum(cx, &candidate, region_vertices)?;
                result.goal_checks = result.goal_checks.saturating_add(1);
                observe(result.primal_iterations, &checked);
                poll(cx, result.primal_iterations)?;
                if let Some(bound) = checked.algebraic_half_width_k()
                    && bound < best_bound
                {
                    best_bound = bound;
                    result.temperature = candidate;
                    result.analysis = checked;
                }
                if best_bound <= config.absolute_tolerance {
                    result.stop = LinearGoalStop::GoalTolerance;
                    break 'attempt;
                }
                if result.primal_iterations == config.max_primal_iterations {
                    result.stop = LinearGoalStop::IterationBudget;
                    break 'attempt;
                }
                if state.iters == before || state.history.len() == cycles_before
                    || !state.rel_residual().is_finite()
                    || state.rel_residual() <= f64::EPSILON
                {
                    break;
                }
            }
            if result.primal_iterations == config.max_primal_iterations {
                result.stop = LinearGoalStop::IterationBudget;
                break;
            }
            if best_bound >= start_bound {
                result.stop = LinearGoalStop::NoProgress;
                break;
            }
            if result.defect_corrections == config.max_defect_corrections {
                result.stop = LinearGoalStop::DefectCorrectionBudget;
                break;
            }
            result.defect_corrections += 1;
        }
        poll(cx, result.primal_iterations)?;
        Ok(result)
    }
}

/// Borrowed low-rank operator. A trait apply cannot return an error, so a
/// first fault is latched, output is poisoned, and the driver MUST check the
/// latch after the Krylov call. Cancellation never masquerades as stagnation.
struct FeedbackOp<'a, 'cx> {
    matrix: &'a Csr,
    injection: &'a Csr,
    feedback: &'a Csr,
    cx: &'a Cx<'cx>,
    ports: RefCell<Vec<f64>>,
    fault: RefCell<Option<ConductionError>>,
}

impl FeedbackOp<'_, '_> {
    fn check(&self) -> Result<(), ConductionError> {
        match self.fault.borrow_mut().take() {
            Some(error) => Err(error),
            None => poll(self.cx, 0),
        }
    }

    fn apply_checked(&self, x: &[f64], y: &mut [f64], transpose: bool)
        -> Result<(), ConductionError>
    {
        poll(self.cx, 0)?;
        let mut work = Poll { cx: self.cx, visits: 0 };
        let mut ports = self.ports.borrow_mut();
        if transpose {
            // A^T x - C^T (B^T x), using the ACTUAL stored transpose.
            for value in ports.iter_mut() { work.tick()?; *value = 0.0; }
            for value in y.iter_mut() { work.tick()?; *value = 0.0; }
            for (i, &xi) in x.iter().enumerate() {
                work.tick()?;
                let (indices, values) = self.injection.row(i);
                for (&j, &a) in indices.iter().zip(values) {
                    work.tick()?; ports[j] = finite(a.mul_add(xi, ports[j]))?;
                }
                let (indices, values) = self.matrix.row(i);
                for (&j, &a) in indices.iter().zip(values) {
                    work.tick()?; y[j] = finite(a.mul_add(xi, y[j]))?;
                }
            }
            for (i, &value) in ports.iter().enumerate() {
                work.tick()?;
                let (indices, values) = self.feedback.row(i);
                for (&j, &a) in indices.iter().zip(values) {
                    work.tick()?; y[j] = finite((-a).mul_add(value, y[j]))?;
                }
            }
        } else {
            for (j, out) in ports.iter_mut().enumerate() {
                work.tick()?;
                *out = row_dot(self.feedback, j, x, 0.0, 1.0, &mut work)?;
            }
            for (i, out) in y.iter_mut().enumerate() {
                work.tick()?;
                let ax = row_dot(self.matrix, i, x, 0.0, 1.0, &mut work)?;
                *out = row_dot(self.injection, i, &ports, ax, -1.0, &mut work)?;
            }
        }
        poll(self.cx, 0)
    }

    fn apply_or_latch(&self, x: &[f64], y: &mut [f64], transpose: bool) {
        if self.fault.borrow().is_some() {
            y.fill(f64::NAN);
            return;
        }
        if let Err(error) = self.apply_checked(x, y, transpose) {
            *self.fault.borrow_mut() = Some(error);
            y.fill(f64::NAN);
        }
    }

    fn residual(&self, rhs: &[f64], offset: &[f64], x: &[f64])
        -> Result<Vec<f64>, ConductionError>
    {
        poll(self.cx, 0)?;
        let mut work = Poll { cx: self.cx, visits: 0 };
        let mut ports = self.ports.borrow_mut();
        for (j, out) in ports.iter_mut().enumerate() {
            work.tick()?;
            *out = row_dot(self.feedback, j, x, offset[j], 1.0, &mut work)?;
        }
        let mut residual = zeros(self.cx, rhs.len())?;
        for (i, out) in residual.iter_mut().enumerate() {
            work.tick()?;
            let r = row_dot(self.matrix, i, x, rhs[i], -1.0, &mut work)?;
            *out = row_dot(self.injection, i, &ports, r, 1.0, &mut work)?;
        }
        poll(self.cx, 0)?;
        Ok(residual)
    }
}

impl LinearOp for FeedbackOp<'_, '_> {
    fn n(&self) -> usize { self.matrix.nrows() }
    fn apply(&self, x: &[f64], y: &mut [f64]) { self.apply_or_latch(x, y, false); }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) { self.apply_or_latch(x, y, true); }
}

struct Poll<'a, 'cx> { cx: &'a Cx<'cx>, visits: usize }
impl Poll<'_, '_> {
    fn tick(&mut self) -> Result<(), ConductionError> {
        self.visits += 1;
        if self.visits == 512 { self.visits = 0; poll(self.cx, 0)?; }
        Ok(())
    }
}

fn row_dot(
    matrix: &Csr, row: usize, x: &[f64], mut value: f64, sign: f64, work: &mut Poll<'_, '_>,
) -> Result<f64, ConductionError> {
    let (indices, values) = matrix.row(row);
    for (&j, &a) in indices.iter().zip(values) {
        work.tick()?;
        value = finite((sign * a).mul_add(x[j], value))?;
    }
    Ok(value)
}

fn finite(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() { Ok(value) }
    else { Err(ConductionError::NonFinite { field: "coupled goal correction", bits: value.to_bits() }) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    fn matrix(rows: &[&[f64]]) -> Csr {
        let mut coo = fs_sparse::Coo::new(rows.len(), rows[0].len());
        for (i, row) in rows.iter().enumerate() {
            for (j, &value) in row.iter().enumerate() { coo.push(i, j, value); }
        }
        coo.assemble()
    }
    fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
        let gate = CancelGate::new_clock_free();
        fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
            let cx = Cx::new(&gate, arena, StreamKey { seed: 73, kernel_id: 740,
                tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic);
            f(&gate, &cx);
        });
    }

    #[test]
    fn low_rank_operator_and_actual_transpose_match_dense_actions() {
        with_cx(|_, cx| {
            let a = matrix(&[&[4.0, 1.0], &[2.0, 3.0]]);
            let b = matrix(&[&[1.0], &[2.0]]);
            let c = matrix(&[&[0.25, -0.5]]);
            let op = FeedbackOp { matrix: &a, injection: &b, feedback: &c, cx,
                ports: RefCell::new(vec![0.0]), fault: RefCell::new(None) };
            // A-B*C = [[3.75, 1.5], [1.5, 4]], except use asymmetric A
            // below to ensure a copied forward action cannot pass as transpose.
            let mut x = [0.0; 2];
            op.apply(&[2.0, -1.0], &mut x); op.check().unwrap();
            assert_eq!(x, [6.0, -1.0]);
            let a = matrix(&[&[4.0, 1.0], &[3.0, 3.0]]);
            let op = FeedbackOp { matrix: &a, injection: &b, feedback: &c, cx,
                ports: RefCell::new(vec![0.0]), fault: RefCell::new(None) };
            op.apply_transpose(&[2.0, -1.0], &mut x); op.check().unwrap();
            assert_eq!(x, [5.0, -1.0]);
            assert_eq!(op.residual(&[1.0, 2.0], &[3.0], &[2.0, -1.0]).unwrap(), vec![-2.0, 7.0]);
        });
    }

    #[test]
    fn cancellation_inside_operator_is_an_error_not_numerical_stagnation() {
        with_cx(|gate, cx| {
            let a = matrix(&[&[2.0]]); let b = matrix(&[&[1.0]]); let c = matrix(&[&[0.5]]);
            let op = FeedbackOp { matrix: &a, injection: &b, feedback: &c, cx,
                ports: RefCell::new(vec![0.0]), fault: RefCell::new(None) };
            gate.request();
            let mut y = [0.0]; op.apply(&[1.0], &mut y);
            assert!(y[0].is_nan());
            assert!(matches!(op.check(), Err(ConductionError::Cancelled { .. })));
        });
    }
}
