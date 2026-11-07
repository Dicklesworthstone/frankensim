//! Allocation-free coupling probes of the original analytic shear path.
use super::{
    HertzMindlinContact, HertzMindlinState, MindlinPathError, PreparedPath, TangentialLaw,
    TangentialStepInput, TriboError, checked_add, checked_mul, finite, poll,
};

/// Immutable scalar path response for an outer mechanical solve.
///
/// This has no candidate history. The final `advance_coupled_path` must still
/// admit the retained-memory ceiling before any mechanical state is published.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MindlinPathTrial {
    /// Actual affine shear travel after the original scalar recoil solve [m].
    pub slip_increment_m: f64,
    /// Work-conjugate path-average effort entering the contact [N].
    pub average_force_n: f64,
    /// Final constitutive effort, distinct from the path average [N].
    pub endpoint_force_n: f64,
    /// Derivative of average effort with respect to ACTUAL shear travel [N/m].
    /// It is not the derivative with respect to released body displacement.
    pub average_tangent_n_per_m: f64,
    pub stored_change_j: f64,
    pub dissipated_j: f64,
    pub closure_residual_j: f64,
    pub compatibility_residual_m: f64,
    pub iterations: usize,
}

impl HertzMindlinContact {
    /// Probe the same simultaneous indentation/shear solve without allocating,
    /// cloning, retaining or publishing any reversal memory.
    ///
    /// This shares the original input/force/iteration admission and analytic
    /// response with `advance_coupled_path`. Normal motion is prescribed; the
    /// optional scalar compliance resolves the actual shear travel. Only the
    /// final history-capacity check is deferred to the publishing operation.
    pub fn trial_coupled_path(
        &self,
        state: &HertzMindlinState,
        input: &TangentialStepInput,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MindlinPathTrial, MindlinPathError> {
        self.admit_path(state, input, false, cancelled)?;
        let (increment, iterations, average, endpoint, derivative, change, heat) =
            if input.normal_n == state.normal_force_n() {
                if let Some(fixed) = &state.fixed {
                    fixed.owner.admit_step(fixed, input)?;
                    let (increment, iterations) = if fixed.owner.capacity_n == 0.0 {
                        (input.slip_increment_m[0], 0)
                    } else {
                        fixed.owner.solve_increment(fixed, input)?
                    };
                    let response = fixed.owner.response(fixed, increment, None)?;
                    (
                        increment,
                        iterations,
                        response.average_n,
                        response.endpoint_n,
                        response.derivative_n_per_m,
                        response.stored_change_j,
                        response.heat_j,
                    )
                } else {
                    (input.slip_increment_m[0], 0, 0.0, 0.0, 0.0, 0.0, 0.0)
                }
            } else if self.parameters.friction_coefficient == 0.0 {
                // Preserve the original geometry/representability checks even
                // though this explicit zero-friction law carries no shear.
                if input.normal_n > 0.0 {
                    self.normal_geometry(input.normal_n)?;
                }
                (input.slip_increment_m[0], 0, 0.0, 0.0, 0.0, 0.0, 0.0)
            } else {
                let prepared = PreparedPath::new(self, state, input.normal_n)?;
                let (increment, iterations) = prepared.solve_increment(state, input, cancelled)?;
                let response = prepared.response(state, increment, None, cancelled)?;
                (
                    increment,
                    iterations,
                    response.average_n,
                    response.endpoint_n,
                    response.derivative_n_per_m,
                    response.energy_j - self.stored_energy_j(state),
                    response.heat_j,
                )
            };
        if average.abs() > input.solve.maximum_force_n {
            return Err(TriboError::OutsideApplicability {
                field: "maximum_force_n",
                value: average.abs(),
                minimum: 0.0,
                maximum: input.solve.maximum_force_n,
            }
            .into());
        }
        checked_add(state.slip_m, increment, "mindlin_accumulated_slip_m")?;
        let work = checked_mul(average, increment, "mindlin_path_shear_work_j")?;
        let closure_residual_j = change + heat - work;
        let compatibility_residual_m = increment
            - (input.slip_increment_m[0]
                - checked_mul(input.compliance_m_per_n, average, "mindlin_body_response_m")?);
        finite(closure_residual_j, "mindlin_path_closure_j")?;
        finite(compatibility_residual_m, "mindlin_compatibility_residual_m")?;
        poll(cancelled)?;
        Ok(MindlinPathTrial {
            slip_increment_m: increment,
            average_force_n: average,
            endpoint_force_n: endpoint,
            average_tangent_n_per_m: derivative,
            stored_change_j: change,
            dissipated_j: heat,
            closure_residual_j,
            compatibility_residual_m,
            iterations,
        })
    }
}
