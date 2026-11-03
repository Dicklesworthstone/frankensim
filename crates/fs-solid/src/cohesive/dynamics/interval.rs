//! Bounded intervals of the original energy-audited cohesive midpoint owner.
//!
//! A rejected Coulomb reversal brackets admissible shorter candidates from the
//! SAME predecessor. No rejected damage, interpolated velocity, artificial
//! impulse or numerical heat enters the trajectory. This localizes a boundary
//! of the discrete straight-path model; it does not certify the first continuous
//! event, monotonicity of a guard, or a unique softening branch.

use super::{
    CohesiveEnergy, CohesiveMotion, DynamicsError, DynamicsPatchDissipation, DynamicsStep, Error,
    finite, positive,
};
use fs_exec::Cx;

/// Numerical subdivision and cumulative work limits, independent of each
/// original native step's unchanged physical acceptance tolerances.
#[derive(Clone, Copy, Debug)]
pub struct IntervalConfig {
    pub minimum_step_s: f64,
    /// Width of the admitted/refused time bracket, not a certified event error.
    pub reversal_time_tolerance_s: f64,
    pub maximum_trials: usize,
    pub maximum_substeps: usize,
    /// Each attempt reserves its remaining native residual allowance. Successful
    /// attempts refund unused evaluations; failures charge their full allowance.
    /// Each evaluation visits every original material patch. Native final energy
    /// checks additionally use four mass/elastic operator applications per trial.
    pub maximum_total_residual_evaluations: usize,
    /// The same conservative reservation policy applies to nonlinear iterations.
    pub maximum_total_newton_iterations: usize,
}

/// Complete requested interval. All physical losses are sums of accepted native
/// paths, including every original patch; they are not endpoint-slip estimates.
pub struct DynamicsInterval<'m, 'a> {
    pub motion: CohesiveMotion<'m, 'a>,
    pub accepted_substeps: usize,
    pub attempted_substeps: usize,
    pub charged_residual_evaluations: usize,
    pub charged_newton_iterations: usize,
    pub maximum_force_residual_n: f64,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    pub bulk_stored_before_j: f64,
    pub bulk_stored_after_j: f64,
    pub cohesive_stored_before_j: f64,
    pub cohesive_stored_after_j: f64,
    pub cohesive_energy_j: CohesiveEnergy,
    pub viscous_dissipated_j: f64,
    pub friction_dissipated_j: f64,
    pub patch_dissipation: Vec<DynamicsPatchDissipation>,
    /// Original patch order. Maximum over the incoming jump and every accepted
    /// straight local material path; rejected stopping trials do not contribute.
    pub maximum_normal_jump_m: Vec<f64>,
    pub external_work_j: f64,
    pub prescribed_jump_work_j: f64,
    /// Independently reconstructed from interval endpoints and accepted work.
    pub energy_residual_j: f64,
    /// Sum of the original accepted native first-law allowances.
    pub energy_tolerance_j: f64,
}

impl<'m, 'a> CohesiveMotion<'m, 'a> {
    /// Complete a requested interval with a CONSTANT applied force and affine
    /// prescribed offsets. Every accepted substep uses the original operators,
    /// material histories, source coefficients and physical tolerances. The
    /// ordinary maximum step applies to each substep; the interval can be longer.
    ///
    /// Only `UnresolvedFrictionReversal` triggers refinement. All other failures,
    /// cancellation and outer-budget exhaustion return no partially adopted
    /// interval. Stopping/capture retains the original velocity tolerance and
    /// actual residual kinetic energy; no velocity is projected to zero.
    pub fn advance_interval(
        &self,
        constant_force_n: &[f64],
        endpoint_prescribed_jumps_m: &[[f64; 3]],
        duration_s: f64,
        controls: IntervalConfig,
        cx: &Cx<'_>,
    ) -> Result<DynamicsInterval<'m, 'a>, DynamicsError> {
        self.advance_interval_controlled(
            constant_force_n,
            endpoint_prescribed_jumps_m,
            duration_s,
            controls,
            &mut || cx.checkpoint().is_err(),
        )
    }

    pub fn advance_interval_controlled(
        &self,
        constant_force_n: &[f64],
        endpoint_prescribed_jumps_m: &[[f64; 3]],
        duration_s: f64,
        controls: IntervalConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DynamicsInterval<'m, 'a>, DynamicsError> {
        poll(cancelled)?;
        positive(duration_s)?;
        positive(controls.minimum_step_s)?;
        positive(controls.reversal_time_tolerance_s)?;
        if controls.maximum_trials == 0
            || controls.maximum_substeps == 0
            || controls.maximum_total_residual_evaluations == 0
            || controls.maximum_total_newton_iterations == 0
        {
            return Err(Error::Budget("nonzero cohesive interval work limits required").into());
        }
        let end_time = finite(self.time_s + duration_s)?;
        if end_time <= self.time_s
            || constant_force_n.len() != self.state.system.n
            || endpoint_prescribed_jumps_m.len() != self.state.points.len()
        {
            return Err(Error::Invalid("cohesive interval clock or input dimensions").into());
        }
        let mut initial_offsets = Vec::with_capacity(self.state.points.len());
        for (patch, point) in self.state.system.patches.iter().zip(&self.state.points) {
            poll(cancelled)?;
            let free = patch.jump(&self.state.displacements, [0.0; 3])?;
            let jump = point.jump_m();
            initial_offsets.push([jump[0] - free[0], jump[1] - free[1], jump[2] - free[2]]);
        }
        let mut offsets = initial_offsets.clone();
        let mut work = Work {
            controls,
            trials: 0,
            residuals: 0,
            iterations: 0,
        };
        let mut current = self.clone();
        let mut result: Option<DynamicsInterval<'m, 'a>> = None;
        while current.time_s < end_time {
            poll(cancelled)?;
            if result.as_ref().map_or(0, |r| r.accepted_substeps) >= controls.maximum_substeps {
                return Err(Error::Budget("cohesive interval substeps").into());
            }
            let remaining = end_time - current.time_s;
            let h = remaining.min(self.config.maximum_step_s);
            let mut trial = |h: f64, work: &mut Work| {
                // Re-adding a rounded remainder can overshoot by one ULP.
                // Keep the exact requested endpoint without changing native h.
                let at = if h == remaining {
                    end_time
                } else {
                    finite(current.time_s + h)?
                };
                let fraction = if at == end_time {
                    1.0
                } else {
                    (at - self.time_s) / duration_s
                };
                for (p, offset) in offsets.iter_mut().enumerate() {
                    poll(cancelled)?;
                    for d in 0..3 {
                        offset[d] = if fraction == 1.0 {
                            endpoint_prescribed_jumps_m[p][d]
                        } else {
                            finite(
                                initial_offsets[p][d]
                                    + fraction
                                        * (endpoint_prescribed_jumps_m[p][d]
                                            - initial_offsets[p][d]),
                            )?
                        };
                    }
                }
                work.attempt(&current, constant_force_n, &offsets, h, at, cancelled)
            };
            let accepted = match trial(h, &mut work) {
                Ok(step) => step,
                Err(DynamicsError::Model(Error::UnresolvedFrictionReversal { mut patch })) => {
                    let mut low = 0.0;
                    let mut high = h;
                    let mut best = None;
                    loop {
                        if low > 0.0 && high - low <= controls.reversal_time_tolerance_s {
                            break best
                                .ok_or(Error::Invalid("missing admitted reversal bracket"))?;
                        }
                        let mid = low + (high - low) * 0.5;
                        if mid < controls.minimum_step_s || mid <= low || mid >= high {
                            return Err(Error::Budget("cohesive reversal time resolution").into());
                        }
                        match trial(mid, &mut work) {
                            Ok(step) => {
                                if let Some(binding) = &current.coulomb {
                                    if binding.resolved_stop(&current, &step.motion, mid, patch)? {
                                        break step;
                                    }
                                }
                                low = mid;
                                best = Some(step);
                            }
                            Err(DynamicsError::Model(Error::UnresolvedFrictionReversal {
                                patch: next_patch,
                            })) => {
                                high = mid;
                                patch = next_patch;
                            }
                            Err(error) => return Err(error),
                        }
                    }
                }
                Err(error) => return Err(error),
            };
            current = accepted.motion.clone();
            if let Some(total) = &mut result {
                total.append(accepted)?;
            } else {
                result = Some(DynamicsInterval::first(accepted, self));
            }
        }
        let mut result = result.ok_or(Error::Invalid("empty cohesive interval"))?;
        result.attempted_substeps = work.trials;
        result.charged_residual_evaluations = work.residuals;
        result.charged_newton_iterations = work.iterations;
        result.energy_residual_j = finite(
            result.external_work_j + result.prescribed_jump_work_j
                - (result.kinetic_after_j - result.kinetic_before_j)
                - (result.bulk_stored_after_j - result.bulk_stored_before_j)
                - (result.cohesive_stored_after_j - result.cohesive_stored_before_j)
                - result.cohesive_energy_j.dissipated
                - result.viscous_dissipated_j
                - result.friction_dissipated_j,
        )?;
        if result.energy_residual_j.abs() > result.energy_tolerance_j {
            return Err(DynamicsError::Energy {
                residual_j: result.energy_residual_j,
                tolerance_j: result.energy_tolerance_j,
            });
        }
        poll(cancelled)?;
        Ok(result)
    }
}

struct Work {
    controls: IntervalConfig,
    trials: usize,
    residuals: usize,
    iterations: usize,
}
impl Work {
    fn attempt<'m, 'a>(
        &mut self,
        motion: &CohesiveMotion<'m, 'a>,
        force: &[f64],
        offsets: &[[f64; 3]],
        h: f64,
        at: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DynamicsStep<'m, 'a>, DynamicsError> {
        poll(cancelled)?;
        if self.trials >= self.controls.maximum_trials {
            return Err(Error::Budget("cohesive interval trials").into());
        }
        let residuals = motion
            .config
            .solve
            .max_operator_evaluations
            .min(self.controls.maximum_total_residual_evaluations - self.residuals);
        let iterations = motion
            .config
            .solve
            .max_newton_iterations
            .min(self.controls.maximum_total_newton_iterations - self.iterations);
        if residuals == 0 || iterations == 0 {
            return Err(Error::Budget("cohesive interval cumulative work").into());
        }
        self.trials += 1;
        self.residuals += residuals;
        self.iterations += iterations;
        let mut limited = motion.clone();
        limited.config.solve.max_operator_evaluations = residuals;
        limited.config.solve.max_newton_iterations = iterations;
        let mut step = limited.advance_at_controlled(force, offsets, h, at, cancelled)?;
        self.residuals -= residuals - step.residual_evaluations;
        self.iterations -= iterations - step.nonlinear.iterations;
        step.motion.config = motion.config;
        Ok(step)
    }
}

impl<'m, 'a> DynamicsInterval<'m, 'a> {
    fn first(s: DynamicsStep<'m, 'a>, incoming: &CohesiveMotion<'m, 'a>) -> Self {
        let maximum_normal_jump_m = incoming
            .state
            .points
            .iter()
            .zip(&s.motion.state.points)
            .map(|(before, after)| before.jump_m()[0].max(after.jump_m()[0]))
            .collect();
        Self {
            motion: s.motion,
            accepted_substeps: 1,
            attempted_substeps: 0,
            charged_residual_evaluations: 0,
            charged_newton_iterations: 0,
            maximum_force_residual_n: s.maximum_force_residual_n,
            kinetic_before_j: s.kinetic_before_j,
            kinetic_after_j: s.kinetic_after_j,
            bulk_stored_before_j: s.bulk_stored_before_j,
            bulk_stored_after_j: s.bulk_stored_after_j,
            cohesive_stored_before_j: s.cohesive_stored_before_j,
            cohesive_stored_after_j: s.cohesive_stored_after_j,
            cohesive_energy_j: s.cohesive_energy_j,
            viscous_dissipated_j: s.viscous_dissipated_j,
            friction_dissipated_j: s.friction_dissipated_j,
            patch_dissipation: s.patch_dissipation,
            maximum_normal_jump_m,
            external_work_j: s.external_work_j,
            prescribed_jump_work_j: s.prescribed_jump_work_j,
            energy_residual_j: s.energy_residual_j,
            energy_tolerance_j: s.energy_tolerance_j,
        }
    }
    fn append(&mut self, s: DynamicsStep<'m, 'a>) -> Result<(), Error> {
        self.accepted_substeps += 1;
        self.maximum_force_residual_n = self
            .maximum_force_residual_n
            .max(s.maximum_force_residual_n);
        self.kinetic_after_j = s.kinetic_after_j;
        self.bulk_stored_after_j = s.bulk_stored_after_j;
        self.cohesive_stored_after_j = s.cohesive_stored_after_j;
        self.cohesive_energy_j.work =
            finite(self.cohesive_energy_j.work + s.cohesive_energy_j.work)?;
        self.cohesive_energy_j.stored_change =
            finite(self.cohesive_energy_j.stored_change + s.cohesive_energy_j.stored_change)?;
        self.cohesive_energy_j.dissipated =
            finite(self.cohesive_energy_j.dissipated + s.cohesive_energy_j.dissipated)?;
        self.cohesive_energy_j.defect =
            finite(self.cohesive_energy_j.defect + s.cohesive_energy_j.defect)?;
        self.viscous_dissipated_j = finite(self.viscous_dissipated_j + s.viscous_dissipated_j)?;
        self.friction_dissipated_j = finite(self.friction_dissipated_j + s.friction_dissipated_j)?;
        self.external_work_j = finite(self.external_work_j + s.external_work_j)?;
        self.prescribed_jump_work_j =
            finite(self.prescribed_jump_work_j + s.prescribed_jump_work_j)?;
        self.energy_tolerance_j = finite(self.energy_tolerance_j + s.energy_tolerance_j)?;
        for (a, b) in self.patch_dissipation.iter_mut().zip(s.patch_dissipation) {
            a.fracture_j = finite(a.fracture_j + b.fracture_j)?;
            a.viscous_j = finite(a.viscous_j + b.viscous_j)?;
            a.friction_j = finite(a.friction_j + b.friction_j)?;
        }
        for (maximum, point) in self
            .maximum_normal_jump_m
            .iter_mut()
            .zip(&s.motion.state.points)
        {
            *maximum = maximum.max(point.jump_m()[0]);
        }
        self.motion = s.motion;
        Ok(())
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
