//! Rebuild the original operator without resetting physical material memory.
use super::*;

/// Explicit allowance for a source-controlled change at an accepted instant.
#[derive(Clone, Copy, Debug)]
pub struct MaterialChangeBudget {
    /// Maximum magnitude of signed external coefficient work [J]. This is a
    /// source allowance, not a tolerance or a license to deposit it as heat.
    pub max_absolute_parameter_work_j: f64,
    /// Relative factor-memory reconstruction error, strictly positive and no
    /// greater than 1e-8. Ill-conditioned or rank-deficient recovery refuses.
    pub relative_memory_tolerance: f64,
}

impl TensorMaterialHistory<'_> {
    // Source provenance and corresponding branch identity are checked by the
    // mesh/source binding owner before this numerical continuation is called.
    pub(super) fn rebind_body<'b>(
        &self,
        body: &'b TensorMaterialBody,
        budget: MaterialChangeBudget,
        gate: &CancelGate,
    ) -> Result<(TensorMaterialHistory<'b>, f64)> {
        poll(gate)?;
        if !budget.max_absolute_parameter_work_j.is_finite()
            || budget.max_absolute_parameter_work_j < 0.0
            || !budget.relative_memory_tolerance.is_finite()
            || budget.relative_memory_tolerance <= 0.0
            || budget.relative_memory_tolerance > 1e-8
            || !self.body.plane_points.is_empty()
            || !body.plane_points.is_empty()
            || self.body.coordinates() != body.coordinates()
            || self.body.points.len() != body.points.len()
            || self.body.inertia.matrix() != body.inertia.matrix()
        {
            return Err(bad(
                "material continuation requires fixed full-3D kinematics/inertia and bounded work/memory controls",
            ));
        }
        let n = body.coordinates();
        let mut state = self.state[..2 * n].to_vec();
        let mut offset = 2 * n;
        for (old, next) in self.body.points.iter().zip(&body.points) {
            poll(gate)?;
            let count = 6 * old.branches().len();
            let memory = old
                .transfer_factor_memory(
                    next,
                    &self.state[offset..offset + count],
                    budget.relative_memory_tolerance,
                    || gate.is_requested(),
                )
                .map_err(TensorMaterialError::Material)?;
            state.extend(memory);
            offset += count;
        }
        body.admit_state(&state)?;
        let change = body.system.hamiltonian(&state) - self.stored_energy_j();
        let parameter_work = self.parameter_work + change;
        if !change.is_finite()
            || change.abs() > budget.max_absolute_parameter_work_j
            || !parameter_work.is_finite()
        {
            return Err(bad(
                "source coefficient work exceeds its finite external allowance",
            ));
        }
        poll(gate)?;
        Ok((
            TensorMaterialHistory {
                body,
                state,
                time: self.time,
                steps: self.steps,
                work: self.work,
                parameter_work,
                thermal_coefficient_work: self.thermal_coefficient_work,
                heat: self.heat,
                initial_energy: self.initial_energy,
            },
            change,
        ))
    }
}
