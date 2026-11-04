//! Coefficient continuation at the same instant of the original trajectory.

use super::*;
use crate::linear3::QuadraticThermalStrainField;
use crate::linear3::visco::source::ViscoSourceError;

/// Candidate motion with the same time, velocity, memories and lifetime budget.
/// The caller publishes it only after this complete continuation succeeds.
pub struct ViscoDynamicMaterialChange<'p, 'a> {
    pub trajectory: ViscoDynamicState<'p, 'a>,
    /// Exact signed material storage change at this instant [J]. No time passes,
    /// no impulse is applied, and neither Kelvin nor Maxwell heat is produced.
    pub parameter_work_j: f64,
}

fn same_controls(a: ViscoDynamicConfig, b: ViscoDynamicConfig) -> bool {
    a.maximum_steps == b.maximum_steps
        && a.maximum_free_dofs == b.maximum_free_dofs
        && a.solve.maximum_iterations == b.solve.maximum_iterations
        && a.solve.relative_residual_tolerance.to_bits()
            == b.solve.relative_residual_tolerance.to_bits()
        && a.absolute_energy_tolerance_j.to_bits() == b.absolute_energy_tolerance_j.to_bits()
        && a.relative_energy_tolerance.to_bits() == b.relative_energy_tolerance.to_bits()
}

impl ViscoDynamicState<'_, '_> {
    /// Continue motion using a rebuilt original-source material operator.
    ///
    /// Prepare `next` from this solid's `at_query_points` result with the SAME
    /// duration, trajectory limits and solve/energy controls. The existing
    /// source owner proves the original geometry, supports, material frames,
    /// branch ordering and reference mass. Physical Maxwell spring strains,
    /// displacement, velocity, time and accepted step count all survive. The
    /// remaining trajectory budget therefore cannot be restarted or enlarged.
    ///
    /// This is instantaneous prescribed parameter work at fixed kinematics.
    /// Coefficients remain frozen within the following native time intervals;
    /// refinement of that source/time splitting remains the caller's task.
    /// An accepted thermal field needs explicit unchanged replacement receipts
    /// through `rebind_materials_with_thermal_strain` instead.
    pub fn rebind_materials<'p, 'a>(
        &self,
        next: &'p ViscoDynamicOperator<'a>,
        cx: &Cx<'_>,
    ) -> Result<ViscoDynamicMaterialChange<'p, 'a>, ViscoSourceError> {
        self.rebind_dynamic_materials(next, None, cx)
    }

    /// Continue source coefficients while retaining every accepted thermal
    /// sample and its B-bar projection. Only the elastic receipt identity may
    /// change. The history owner requires each current thermal temperature to
    /// equal its new element coefficient-query temperature; it never averages
    /// a current temperature gradient into an elementwise source coefficient.
    /// Original thermal reference temperatures, eigenstrain work and all other
    /// trajectory ledgers are retained. Refusal leaves this motion untouched.
    pub fn rebind_materials_with_thermal_strain<'p, 'a, 't>(
        &self,
        next: &'p ViscoDynamicOperator<'a>,
        thermal: impl Into<QuadraticThermalStrainField<'t>>,
        cx: &Cx<'_>,
    ) -> Result<ViscoDynamicMaterialChange<'p, 'a>, ViscoSourceError> {
        self.rebind_dynamic_materials(next, Some(thermal.into()), cx)
    }

    fn rebind_dynamic_materials<'p, 'a>(
        &self,
        next: &'p ViscoDynamicOperator<'a>,
        thermal: Option<QuadraticThermalStrainField<'_>>,
        cx: &Cx<'_>,
    ) -> Result<ViscoDynamicMaterialChange<'p, 'a>, ViscoSourceError> {
        checkpoint(cx)?;
        if next.dt_s.to_bits() != self.prepared.dt_s.to_bits()
            || !same_controls(next.config, self.prepared.config)
        {
            return Err(ViscoSourceError::Invalid(
                "material continuation requires the original dynamic duration and all trajectory controls",
            ));
        }
        let changed = self
            .material
            .rebind_materials_impl(next.system, thermal, cx)?;
        let parameter_work_j = finite(self.parameter_work_j + changed.parameter_work_j)?;
        let trajectory = ViscoDynamicState {
            prepared: next,
            material: changed.history,
            velocity_mps: self.velocity_mps.clone(),
            // Source continuation retains exactly the same reference mass and
            // velocity, so both momentum and this computed kinetic value survive.
            kinetic_j: self.kinetic_j,
            external_work_j: self.external_work_j,
            eigenstrain_work_j: self.eigenstrain_work_j,
            parameter_work_j,
            time_s: self.time_s,
            steps: self.steps,
        };
        checkpoint(cx)?;
        Ok(ViscoDynamicMaterialChange {
            trajectory,
            parameter_work_j: changed.parameter_work_j,
        })
    }
}
