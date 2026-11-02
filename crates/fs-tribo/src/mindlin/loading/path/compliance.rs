//! Work-conjugate scalar body response around the original analytic path.

use super::{
    HertzMindlinContact, HertzMindlinState, MindlinPathError, MindlinPathStep, PreparedPath,
    TangentialStepInput, TriboError, checked_add, checked_mul, finite, poll,
};

impl HertzMindlinContact {
    /// Resolve simultaneous normal/shear loading against scalar body compliance.
    ///
    /// The Hertz indentation is affine from the accepted state to the endpoint
    /// `input.normal_n`. The first input slip is RELEASED relative motion, and
    /// the solved actual shear obeys `du = du_free - C * average_traction`.
    /// Normal motion remains prescribed; there is no normal/body feedback or
    /// elastic cross-compliance. Shear is affine along each candidate path.
    ///
    /// The returned average traction drives both the body response and contact
    /// work. For a midpoint inertial step, the mechanical owner can supply
    /// `C = dt^2/(2*m_effective)` and its released relative displacement.
    /// The owner still accepts the trajectory and refines curved body motion.
    ///
    /// Zero compliance delegates exactly to `advance_prescribed_path`, and
    /// fixed normal load delegates to the original `TangentialLaw` solve.
    /// Varying-load trials cache normal geometry and allocate no history. A
    /// nonnegative analytic tangent supports safeguarded Newton iterations within
    /// the caller's force/iteration budgets. Cancellation or refusal leaves the
    /// accepted state untouched. Only the final response materializes memory.
    pub fn advance_coupled_path(
        &self,
        state: &HertzMindlinState,
        input: &TangentialStepInput,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MindlinPathStep, MindlinPathError> {
        self.advance_path(state, input, false, cancelled)
    }
}

struct Trial {
    residual_n: f64,
    tolerance_n: f64,
    increment_m: f64,
    derivative: f64,
    returned_n: f64,
}

impl PreparedPath {
    fn trial(
        &self,
        state: &HertzMindlinState,
        input: &TangentialStepInput,
        force: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Trial, MindlinPathError> {
        poll(cancelled)?;
        let increment_m = input.slip_increment_m[0]
            - checked_mul(input.compliance_m_per_n, force, "mindlin_body_response_m")?;
        finite(increment_m, "mindlin_resolved_increment_m")?;
        let response = self.response(state, increment_m, None, cancelled)?;
        let residual_n = force - response.average_n;
        finite(residual_n, "mindlin_path_force_residual_n")?;
        Ok(Trial {
            residual_n,
            tolerance_n: input.solve.tolerance(force, response.average_n)?,
            increment_m,
            derivative: checked_add(
                1.0,
                checked_mul(
                    input.compliance_m_per_n,
                    response.derivative_n_per_m,
                    "mindlin_path_body_tangent",
                )?,
                "mindlin_path_root_tangent",
            )?,
            returned_n: response.average_n,
        })
    }

    fn admitted(
        &self,
        state: &HertzMindlinState,
        input: &TangentialStepInput,
        trial: &Trial,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, MindlinPathError> {
        if trial.residual_n.abs() > trial.tolerance_n {
            return Ok(false);
        }
        // The body uses the constitutive average, not the root's trial force.
        // Check that returned effort on its own load line; a large C*dQ/du
        // otherwise amplifies an apparently converged force residual.
        let returned = self.trial(state, input, trial.returned_n, cancelled)?;
        Ok(returned.residual_n.abs() <= returned.tolerance_n)
    }

    pub(super) fn solve_increment(
        &self,
        state: &HertzMindlinState,
        input: &TangentialStepInput,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, usize), MindlinPathError> {
        if input.compliance_m_per_n == 0.0 {
            return Ok((input.slip_increment_m[0], 0));
        }
        // Each local spring's force is nondecreasing in imposed shear. Hence
        // F - Qbar(free-C*F) is strictly increasing for C>=0, even through
        // entry, yield, detachment and reversal. The largest endpoint normal
        // load bounds every local capacity along the affine indentation path.
        let bound = self.capacity_n.min(input.solve.maximum_force_n);
        let mut lo = -bound;
        let mut hi = bound;
        for (force, lower) in [(lo, true), (hi, false)] {
            let trial = self.trial(state, input, force, cancelled)?;
            if self.admitted(state, input, &trial, cancelled)? {
                return Ok((trial.increment_m, 0));
            }
            if (lower && trial.residual_n > 0.0) || (!lower && trial.residual_n < 0.0) {
                return Err(TriboError::OutsideApplicability {
                    field: "maximum_force_n",
                    value: trial.returned_n.abs(),
                    minimum: 0.0,
                    maximum: input.solve.maximum_force_n,
                }
                .into());
            }
        }
        let mut force = 0.0;
        let mut final_residual = 0.0;
        let mut final_tolerance = 0.0;
        for iteration in 1..=input.solve.max_iterations {
            let trial = self.trial(state, input, force, cancelled)?;
            if self.admitted(state, input, &trial, cancelled)? {
                return Ok((trial.increment_m, iteration));
            }
            final_residual = trial.residual_n;
            final_tolerance = trial.tolerance_n;
            if trial.residual_n < 0.0 {
                lo = force;
            } else {
                hi = force;
            }
            let proposed = force - trial.residual_n / trial.derivative;
            force = if proposed.is_finite() && proposed > lo && proposed < hi {
                proposed
            } else {
                f64::midpoint(lo, hi)
            };
        }
        Err(TriboError::TangentialSolve {
            residual_n: final_residual,
            tolerance_n: final_tolerance,
            iterations: input.solve.max_iterations,
        }
        .into())
    }
}
