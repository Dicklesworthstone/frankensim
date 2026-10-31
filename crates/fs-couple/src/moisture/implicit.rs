//! Endpoint-dependent Fick transport using the existing material and LU owners.

use super::*;
use fs_la::tridiagonal::{TridiagonalError, TridiagonalWorkspace};
use fs_solver::{CheckedCgSolution, ResidualClaim, StallDiagnosis};

/// Bounded Picard controls for fully implicit moisture-dependent diffusivity.
/// These are numerical controls, not fitted material or physical rate constants.
#[derive(Clone, Copy, Debug)]
pub struct MoistureNonlinearConfig {
    /// True Euclidean residual of the complete endpoint-dependent discrete
    /// equations, divided by the endpoint system's RHS norm, in (0,1).
    pub relative_tolerance: f64,
    /// Maximum source-refresh/linear-solve iterations per proposed interval.
    pub maximum_iterations: usize,
    /// Convex relaxation of the Picard iterate, in (0,1]. A failed iteration
    /// budget refuses the whole interval; neither time nor tolerance is changed.
    pub relaxation: f64,
}
impl MoistureNonlinearConfig {
    pub(super) fn validate(self) -> Result<(), MoistureError> {
        if !positive(self.relative_tolerance)
            || self.relative_tolerance >= 1.0
            || self.maximum_iterations == 0
            || !positive(self.relaxation)
            || self.relaxation > 1.0
        {
            return Err(MoistureError::Invalid(
                "implicit diffusion needs relative tolerance in (0,1), positive iteration budget and relaxation in (0,1]",
            ));
        }
        Ok(())
    }
}

/// Independent residual of endpoint-dependent diffusion and its actual work.
#[derive(Clone, Copy, Debug)]
pub struct MoistureNonlinearReport {
    pub iterations: usize,
    /// All tridiagonal substitutions, including linear defect corrections.
    pub tridiagonal_solves: usize,
    pub relative_residual: f64,
    /// Uses endpoint diffusivity for both surface transfers; never the lagged
    /// linear solve's flux to make nonlinear conservation appear satisfied.
    pub water_balance_residual_kg: f64,
}

fn linear_error(error: TridiagonalError) -> MoistureError {
    match error {
        TridiagonalError::Cancelled => MoistureError::Cancelled,
        error => MoistureError::Tridiagonal(error),
    }
}

impl MoistureSlab {
    pub(super) fn solve_implicit(
        &self,
        dt: f64,
        faces: &[ResolvedMoistureFace; 2],
        start: &[ResolvedMaterialStatePoint],
        step: MoistureStepConfig,
        controls: MoistureNonlinearConfig,
        cx: &Cx<'_>,
    ) -> Result<
        (
            CheckedCgSolution,
            Vec<ResolvedMaterialStatePoint>,
            [(f64, f64); 2],
            Option<MoistureNonlinearReport>,
        ),
        MoistureError,
    > {
        let n = self.moisture.len();
        let mut workspace = TridiagonalWorkspace::new(n).map_err(linear_error)?;
        let mut iterate = self.moisture.clone();
        let mut system = self.transport_system(dt, faces, start, cx)?;
        let mut total_solves = 0;
        for iteration in 1..=controls.maximum_iterations {
            poll(cx)?;
            let solution =
                system
                    .0
                    .solve_tridiagonal(&system.1, step.linear, &mut workspace, cx)?;
            total_solves += solution.report.iters;
            if !solution.report.converged_euclidean() {
                return Err(MoistureError::NotConverged(solution.report));
            }
            // Re-query the ORIGINAL pinned source at every proposed endpoint.
            // Unsupported intermediate or final states refuse transactionally.
            let endpoint = self.transport_states(&solution.x, cx)?;
            let endpoint_system = self.transport_system(dt, faces, &endpoint, cx)?;
            let (_, relative_residual) =
                endpoint_system
                    .0
                    .defect(&solution.x, &endpoint_system.1, cx)?;
            let transfer = self.boundary_transfer(&solution.x, faces, &endpoint_system.2);
            let water_after: f64 = solution
                .x
                .iter()
                .map(|m| self.dry_mass_per_cell_kg * m)
                .sum();
            let water_balance = water_after - self.water_mass_kg() - transfer[0] - transfer[1];
            let report = MoistureNonlinearReport {
                iterations: iteration,
                tridiagonal_solves: total_solves,
                relative_residual,
                water_balance_residual_kg: water_balance,
            };
            if relative_residual < controls.relative_tolerance
                && water_balance.is_finite()
                && water_balance.abs() <= step.water_balance_tolerance_kg
            {
                return Ok((solution, endpoint, endpoint_system.2, Some(report)));
            }
            if iteration == controls.maximum_iterations {
                return Err(MoistureError::NonlinearNotConverged(report));
            }
            if controls.relaxation == 1.0 {
                iterate = solution.x;
                system = endpoint_system;
            } else {
                for (old, candidate) in iterate.iter_mut().zip(&solution.x) {
                    *old = controls.relaxation.mul_add(candidate - *old, *old);
                }
                let states = self.transport_states(&iterate, cx)?;
                system = self.transport_system(dt, faces, &states, cx)?;
            }
        }
        unreachable!("validated positive nonlinear iteration budget")
    }
}

impl DiffusionStep {
    /// Residual against the original conservative operator. Scaled Euclidean
    /// norms avoid squaring the physical load magnitude directly.
    fn defect(
        &self,
        x: &[f64],
        rhs: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, f64), MoistureError> {
        poll(cx)?;
        let mut residual = vec![0.0; rhs.len()];
        LinearOp::apply(self, x, &mut residual);
        let mut rscale: f64 = 0.0;
        let mut bscale: f64 = 0.0;
        for (i, (r, b)) in residual.iter_mut().zip(rhs).enumerate() {
            if i % 128 == 0 {
                poll(cx)?;
            }
            *r = *b - *r;
            if !r.is_finite() {
                return Err(MoistureError::Invalid("nonfinite diffusion residual"));
            }
            rscale = rscale.max(r.abs());
            bscale = bscale.max(b.abs());
        }
        let relative = if rscale == 0.0 {
            0.0
        } else if bscale == 0.0 {
            f64::INFINITY
        } else {
            let rn = fs_solver::norm2(&residual.iter().map(|r| r / rscale).collect::<Vec<_>>());
            let bn = fs_solver::norm2(&rhs.iter().map(|b| b / bscale).collect::<Vec<_>>());
            (rscale / bscale) * (rn / bn)
        };
        poll(cx)?;
        Ok((residual, relative))
    }

    /// Reuse native O(n) adjacent-pivot LU, then check the original equation.
    /// Correction attempts reuse the same factors and stay inside both work
    /// limits. A direct solve supplies no forward-error or continuum certificate.
    fn solve_tridiagonal(
        &self,
        rhs: &[f64],
        controls: CheckedCgConfig,
        workspace: &mut TridiagonalWorkspace,
        cx: &Cx<'_>,
    ) -> Result<CheckedCgSolution, MoistureError> {
        if !positive(controls.tolerance) || controls.tolerance >= 1.0 {
            return Err(MoistureError::Invalid("linear tolerance must lie in (0,1)"));
        }
        let off: Vec<_> = self.off.iter().map(|a| -*a).collect();
        let mut x = vec![0.0; rhs.len()];
        let (mut residual, mut relative) = self.defect(&x, rhs, cx)?;
        let mut history = vec![relative];
        let mut solves = 0;
        let maximum_solves = controls
            .max_iterations
            .min(controls.max_corrections.saturating_add(1));
        if relative >= controls.tolerance && maximum_solves > 0 {
            workspace
                .factor(&off, &self.diagonal, &off, || cx.checkpoint().is_err())
                .map_err(linear_error)?;
        }
        let mut correction = vec![0.0; rhs.len()];
        let mut diagnosis = StallDiagnosis::BudgetExhausted;
        while relative >= controls.tolerance && solves < maximum_solves {
            workspace
                .solve_factored_into(&residual, &mut correction, || cx.checkpoint().is_err())
                .map_err(linear_error)?;
            solves += 1;
            let candidate: Vec<_> = x.iter().zip(&correction).map(|(x, dx)| x + dx).collect();
            let (defect, measured) = self.defect(&candidate, rhs, cx)?;
            if !measured.is_finite() {
                diagnosis = StallDiagnosis::Breakdown;
                break;
            }
            if measured >= relative {
                diagnosis = StallDiagnosis::Plateau;
                break;
            }
            x = candidate;
            residual = defect;
            relative = measured;
            history.push(relative);
        }
        poll(cx)?;
        Ok(CheckedCgSolution {
            x,
            report: SolveReport::from_claim_with_diagnosis(
                solves,
                ResidualClaim::TrueEuclidean(relative),
                controls.tolerance,
                history,
                diagnosis,
            ),
            corrections: solves.saturating_sub(1),
        })
    }
}
