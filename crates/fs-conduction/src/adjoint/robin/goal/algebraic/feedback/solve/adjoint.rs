//! Full affine-feedback pullbacks: (A - B C)^T lambda = weights.
//! Reuse the primal's actual low-rank operator, including its transpose and
//! cancellation latch. No solid-only adjoint or differentiated Krylov trace.

use super::{ConductionError, Cx, FeedbackOp, FgmresState, LinearOp,
    LinearRobinFeedbackAnalyzer, Poll, RefCell, finite, invalid, poll, row_dot, zeros};
use crate::LinearConfig;
use crate::adjoint::robin::{RobinGradient, admit_linear, failed, vector};

struct Transpose<'a, 'op, 'cx>(&'a FeedbackOp<'op, 'cx>);
impl LinearOp for Transpose<'_, '_, '_> {
    fn n(&self) -> usize { self.0.n() }
    fn apply(&self, x: &[f64], y: &mut [f64]) { self.0.apply_or_latch(x, y, true); }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.0.apply_or_latch(x, y, false);
    }
}

impl LinearRobinFeedbackAnalyzer<'_> {
    /// Differentiate a nodal linear functional through the COMPLETE feedback.
    ///
    /// With r = d + D W T, solve (A - B C)^T lambda = weights once.
    /// `references` contains derivatives with respect to additive changes in d,
    /// NOT changes to frozen references; `nodal_load` contains derivatives with
    /// respect to assembled nodal watts (zero on prescribed nodes). `log_htc`
    /// changes a Robin coefficient while holding D, geometry and all other
    /// parameters fixed. A physical air h/flow derivative ALSO differentiates
    /// D and d; the log_htc result alone must not be used for that derivative.
    ///
    /// A linear functional is smooth even at a tie in nodal temperature. A
    /// maximum consumer must separately establish its active vertex or identify
    /// this as the derivative of a selected nodal functional, not of the max.
    /// The log_htc contraction is evaluated at the supplied material-admissible
    /// field; its primal error is not hidden inside the adjoint residual.
    ///
    /// Require the existing complete-system inverse evidence before a nonzero
    /// solve. Reuse prepared inverse/response data, but independently check the
    /// final rescaled adjoint's true relative residual. All restart cycles share
    /// `linear.max_iterations`; restart storage is capped at 32 vectors. This is
    /// an Estimated stored-system derivative, not a gradient-error enclosure,
    /// physical validation, continuum/shape or nonlinear-material derivative.
    ///
    /// # Errors
    /// Invalid controls, weights or field, missing inverse evidence, exhausted
    /// Krylov work, nonfinite arithmetic/allocation, or cancellation. No partial
    /// gradient is returned. Neither the analyzer nor caller inputs are mutated.
    pub fn pullback_affine_controls(
        &self, cx: &Cx<'_>, temperature: &[f64], full_weights: &[f64],
        linear: LinearConfig,
    ) -> Result<RobinGradient, ConductionError> {
        poll(cx, 0)?;
        admit_linear(linear)?;
        if !(1..=32).contains(&linear.restart) {
            return Err(invalid("coupled adjoint restart must be in 1..=32"));
        }
        let dofs = self.solid.dofs();
        super::super::super::validate_field(cx, self.solid.problem, dofs, temperature)?;
        vector(cx, full_weights, temperature.len())?;
        let mut load = zeros(cx, temperature.len())?;
        let mut offsets = zeros(cx, self.ports.len())?;
        let mut log_htc = zeros(cx, self.ports.len())?;
        let mut work = Poll { cx, visits: 0 };
        let mut scale = 0.0_f64;
        for &vertex in dofs.free() {
            work.tick()?;
            scale = scale.max(full_weights[vertex].abs());
        }
        if scale == 0.0 {
            poll(cx, 0)?;
            return Ok(RobinGradient { references: offsets, log_htc, nodal_load: load,
                relative_residual: 0.0, iterations: 0 });
        }
        // The immutable analyzer verifies its OWN inverse. Do not accept a
        // caller-supplied stability scalar or a compatible singular solve.
        let assessment = self.analyze_maximum(cx, temperature, &dofs.free()[..1])?;
        if assessment.coupled().coupled_inverse_infinity_upper().is_none() {
            return Err(invalid("coupled adjoint requires a verified complete-system inverse"));
        }
        let rhs = dofs.free().iter().map(|&v| full_weights[v] / scale).collect::<Vec<_>>();
        let operator = FeedbackOp { matrix: &self.solid.response.matrix,
            injection: &self.injection, feedback: &self.feedback, cx,
            ports: RefCell::new(zeros(cx, self.ports.len())?), fault: RefCell::new(None) };
        let transposed = Transpose(&operator);
        let pre = crate::solve::spd_preconditioner(&self.solid.response.matrix);
        let mut state = FgmresState::new(&rhs, linear.restart.min(dofs.n()));
        let mut residual = 1.0;
        while state.iters < linear.max_iterations {
            poll(cx, state.iters)?;
            let before = state.iters;
            state.restart = linear.restart.min(dofs.n()).min(linear.max_iterations - before);
            state.run(&transposed, &pre, &rhs, linear.tolerance, 1);
            operator.check()?;
            residual = actual_residual(cx, &operator, &state.x, &rhs)?;
            if residual < linear.tolerance || state.iters == before { break; }
        }
        poll(cx, state.iters)?;
        if residual >= linear.tolerance { return Err(failed(state.iters, residual, linear)); }
        let mut rounded = zeros(cx, dofs.n())?;
        for (i, &vertex) in dofs.free().iter().enumerate() {
            work.tick()?;
            load[vertex] = finite(state.x[i] * scale)?;
            // Rescaling can round/underflow. Check the gradient we RETURN,
            // rather than only the normalized iterative candidate.
            rounded[i] = finite(load[vertex] / scale)?;
        }
        residual = actual_residual(cx, &operator, &rounded, &rhs)?;
        if residual >= linear.tolerance { return Err(failed(state.iters, residual, linear)); }
        for (i, &vertex) in dofs.free().iter().enumerate() {
            let (columns, values) = self.injection.row(i);
            for (&j, &value) in columns.iter().zip(values) {
                work.tick()?;
                offsets[j] = finite(value.mul_add(load[vertex], offsets[j]))?;
            }
        }
        let free = dofs.gather(temperature);
        for (j, port) in self.ports.iter().enumerate() {
            let change = row_dot(&self.feedback, j, &free, self.offset[j], 1.0, &mut work)?;
            let reference = finite(port.reference_k + change)?;
            for (vertices, area) in &port.faces {
                let mut sum_lambda = 0.0;
                let mut sum_temperature = 0.0;
                let mut diagonal = 0.0;
                for &v in vertices {
                    work.tick()?;
                    let delta = finite(temperature[v] - reference)?;
                    sum_lambda = finite(sum_lambda + load[v])?;
                    sum_temperature = finite(sum_temperature + delta)?;
                    diagonal = finite(load[v].mul_add(delta, diagonal))?;
                }
                // Consistent P1 face mass: area/12 * (ones + identity).
                let contraction = finite(sum_lambda.mul_add(sum_temperature, diagonal))?;
                log_htc[j] = finite(log_htc[j]
                    - finite(port.htc_w_m2_k * (area / 12.0) * contraction)?)?;
            }
        }
        poll(cx, state.iters)?;
        Ok(RobinGradient { references: offsets, log_htc, nodal_load: load,
            relative_residual: residual, iterations: state.iters })
    }
}

fn actual_residual(
    cx: &Cx<'_>, operator: &FeedbackOp<'_, '_>, x: &[f64], rhs: &[f64],
) -> Result<f64, ConductionError> {
    vector(cx, x, rhs.len())?;
    let mut ax = zeros(cx, rhs.len())?;
    operator.apply_checked(x, &mut ax, true)?;
    let mut scale = 0.0_f64;
    let mut work = Poll { cx, visits: 0 };
    for (&a, &b) in ax.iter().zip(rhs) {
        work.tick()?;
        scale = scale.max(a.abs()).max(b.abs());
    }
    if scale == 0.0 { poll(cx, 0)?; return Ok(0.0); }
    let mut residual = zeros(cx, rhs.len())?;
    let mut normalized = zeros(cx, rhs.len())?;
    for (i, (&a, &b)) in ax.iter().zip(rhs).enumerate() {
        work.tick()?;
        residual[i] = finite(b / scale - a / scale)?;
        normalized[i] = b / scale;
    }
    let value = finite(fs_solver::norm2(&residual)
        / fs_solver::norm2(&normalized).max(f64::MIN_POSITIVE))?;
    poll(cx, 0)?;
    Ok(value)
}
