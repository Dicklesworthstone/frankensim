//! One bounded original-law Gauss/Newton iteration for both spatial dimensions.
//! Geometry and constitutive kernels own their evaluations and checked tangent;
//! the electrode quotient, source work, line search and work limits are shared.
use super::electrodes::Dofs;
use super::*;

pub(super) trait Problem {
    type Evaluation;
    fn new_evaluation(&self) -> Self::Evaluation;
    fn evaluate(
        &self,
        potential: &[f64],
        out: &mut Self::Evaluation,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error>;
    fn nodal_reactions(evaluation: &Self::Evaluation) -> &[f64];
    fn coenergy_change(
        &self,
        before: &Self::Evaluation,
        after: &Self::Evaluation,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error>;
    fn direction(
        &self,
        evaluation: &Self::Evaluation,
        dofs: &Dofs,
        rhs: &[f64],
        linear: CheckedCgConfig,
        maximum_applications: usize,
        applications: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, SolveReport), Error>;
}

impl Problem for PlanarNonlinearElectrostatics {
    type Evaluation = Evaluation;
    fn new_evaluation(&self) -> Evaluation {
        Evaluation::new(self)
    }
    fn evaluate(
        &self,
        potential: &[f64],
        out: &mut Evaluation,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        PlanarNonlinearElectrostatics::evaluate(self, potential, out, cancelled)
    }
    fn nodal_reactions(evaluation: &Evaluation) -> &[f64] {
        &evaluation.nodal
    }
    fn coenergy_change(
        &self,
        before: &Evaluation,
        after: &Evaluation,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        PlanarNonlinearElectrostatics::coenergy_change(self, before, after, cancelled)
    }
    fn direction(
        &self,
        evaluation: &Evaluation,
        dofs: &Dofs,
        rhs: &[f64],
        linear: CheckedCgConfig,
        maximum_applications: usize,
        applications: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, SolveReport), Error> {
        tangent::direction(
            self,
            &evaluation.tangent,
            dofs,
            rhs,
            linear,
            maximum_applications,
            applications,
            cancelled,
        )
    }
}

fn objective_change<P: Problem>(
    problem: &P,
    dofs: &Dofs,
    rhs: &[f64],
    before_v: &[f64],
    before: &P::Evaluation,
    after_v: &[f64],
    after: &P::Evaluation,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let (mut sum, mut correction) = (problem.coenergy_change(before, after, cancelled)?, 0.0);
    for (v, &q) in rhs.iter().enumerate() {
        poll(cancelled)?;
        add(&mut sum, &mut correction, -q * (after_v[v] - before_v[v]))?;
    }
    for (i, &v) in dofs.representatives.iter().enumerate() {
        add(
            &mut sum,
            &mut correction,
            -dofs.charge[i] * (after_v[v] - before_v[v]),
        )?;
    }
    Ok(finite(sum + correction)?)
}

pub(super) fn iterate<P: Problem>(
    problem: &P,
    dofs: &Dofs,
    rhs: &[f64],
    mut potential: Vec<f64>,
    settings: NonlinearElectrostaticControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<f64>, P::Evaluation, NonlinearElectrostaticReport), Error> {
    poll(cancelled)?;
    solve::controls(settings)?;
    let mut evaluation = problem.new_evaluation();
    let mut trial = problem.new_evaluation();
    problem.evaluate(&potential, &mut evaluation, cancelled)?;
    let mut residual = vec![0.0; dofs.len()];
    let mut correction = vec![0.0; dofs.len()];
    dofs.residual(
        P::nodal_reactions(&evaluation),
        rhs,
        &mut residual,
        &mut correction,
    )?;
    let initial_norm = finite(norm(&residual))?;
    let target = finite(
        settings
            .absolute_residual_tolerance_c
            .max(settings.relative_residual_tolerance * initial_norm)
            .min(settings.charge_tolerance_c / (2.0 * (dofs.len().max(1) as f64).sqrt())),
    )?;
    if !positive(target) {
        return Err(invalid("physical charge target is unrepresentable"));
    }
    let mut report = NonlinearElectrostaticReport {
        newton_iterations: 0,
        linear_iterations: 0,
        backtracks: 0,
        operator_applications: 0,
        residual_norm_c: initial_norm,
        residual_target_c: target,
    };
    let mut proposed = potential.clone();
    while report.residual_norm_c > target {
        poll(cancelled)?;
        if report.newton_iterations == settings.max_newton_iterations {
            return Err(Error::Nonconvergence(report));
        }
        // The compensated-residual scratch is idle during the linear solve.
        // Reuse it for -R, without allocating another quotient-sized vector.
        for (out, &r) in correction.iter_mut().zip(&residual) {
            *out = -r;
        }
        let (direction, linear) = problem.direction(
            &evaluation,
            dofs,
            &correction,
            settings.linear,
            settings.max_operator_applications,
            &mut report.operator_applications,
            cancelled,
        )?;
        report.linear_iterations = report
            .linear_iterations
            .checked_add(linear.iters)
            .ok_or_else(|| invalid("nonlinear work count overflow"))?;
        let (mut derivative, mut cd) = (0.0, 0.0);
        for (&r, &d) in residual.iter().zip(&direction) {
            add(&mut derivative, &mut cd, r * d)?;
        }
        let derivative = finite(derivative + cd)?;
        if derivative >= 0.0 {
            return Err(Error::Nonconvergence(report));
        }
        let mut alpha = 1.0_f64;
        let mut accepted = false;
        for backtrack in 0..=settings.max_backtracks {
            poll(cancelled)?;
            for v in 0..proposed.len() {
                proposed[v] = if let Some(i) = dofs.index[v] {
                    finite(alpha.mul_add(direction[i], potential[v]))?
                } else {
                    potential[v]
                };
            }
            match problem.evaluate(&proposed, &mut trial, cancelled) {
                Ok(()) => {
                    let delta = objective_change(
                        problem,
                        dofs,
                        rhs,
                        &potential,
                        &evaluation,
                        &proposed,
                        &trial,
                        cancelled,
                    )?;
                    if delta <= 1e-4 * alpha * derivative {
                        accepted = true;
                        break;
                    }
                    // A trial already satisfying the original Gauss target
                    // needs no numerically lost Armijo decrement. Independent
                    // original field/charge/energy audits still precede return.
                    dofs.residual(
                        P::nodal_reactions(&trial),
                        rhs,
                        &mut residual,
                        &mut correction,
                    )?;
                    if finite(norm(&residual))? <= target {
                        accepted = true;
                        break;
                    }
                }
                Err(Error::FieldRange { .. }) => {}
                Err(error) => return Err(error),
            }
            if backtrack == settings.max_backtracks {
                break;
            }
            alpha *= 0.5;
            report.backtracks = report
                .backtracks
                .checked_add(1)
                .ok_or_else(|| invalid("line search work count overflow"))?;
        }
        if !accepted {
            return Err(Error::Nonconvergence(report));
        }
        std::mem::swap(&mut potential, &mut proposed);
        std::mem::swap(&mut evaluation, &mut trial);
        dofs.residual(
            P::nodal_reactions(&evaluation),
            rhs,
            &mut residual,
            &mut correction,
        )?;
        report.residual_norm_c = finite(norm(&residual))?;
        report.newton_iterations += 1;
    }
    // Recompute the accepted ORIGINAL material response before publication.
    problem.evaluate(&potential, &mut evaluation, cancelled)?;
    dofs.residual(
        P::nodal_reactions(&evaluation),
        rhs,
        &mut residual,
        &mut correction,
    )?;
    report.residual_norm_c = finite(norm(&residual))?;
    if report.residual_norm_c > target {
        return Err(Error::Nonconvergence(report));
    }
    poll(cancelled)?;
    Ok((potential, evaluation, report))
}
