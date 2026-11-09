//! Reuse the original hanging-node map for both sides of the thermal block.

use std::ops::ControlFlow;

use crate::elastic3::ElasticityError3;
use crate::elastic3::adaptive::AdaptiveElasticity3;
use crate::elastic3::thermal::{
    PreparedThermalResidualJacobian3, ThermalCoordinateMap3, ThermalJacobianLimits3,
};

impl AdaptiveElasticity3 {
    /// Freeze the native mixed block in independent master coordinates.
    /// Uses the original T for physical admission and T^T for subsequent force
    /// actions; no expanded global field or copied constraint map is prepared.
    /// Limits and physical/source admission match the Cartesian preparation.
    pub fn prepare_thermal_residual_jacobian(
        &self,
        u: &[f64],
        maximum_kinematic_strain: f64,
        limits: ThermalJacobianLimits3,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<PreparedThermalResidualJacobian3, ElasticityError3> {
        PreparedThermalResidualJacobian3::prepare(
            &self.raw,
            ThermalCoordinateMap3 { rows: Some(&self.rows), fixed: &self.fixed },
            u, maximum_kinematic_strain, limits, None, None, checkpoint,
        )
    }

    /// Include actual embedded motion and its local temperature slope with the
    /// original constraint map. Pure callback, frozen-geometry, homogeneous box
    /// support and preflight rules match the Cartesian with-motion preparation.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_thermal_residual_jacobian_with_motion(
        &self,
        u: &[f64],
        prescribed: &dyn Fn([f64; 3], [f64; 3]) -> [f64; 3],
        prescribed_temperature_derivative: Option<&dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
        maximum_kinematic_strain: f64,
        limits: ThermalJacobianLimits3,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<PreparedThermalResidualJacobian3, ElasticityError3> {
        PreparedThermalResidualJacobian3::prepare(
            &self.raw,
            ThermalCoordinateMap3 { rows: Some(&self.rows), fixed: &self.fixed },
            u, maximum_kinematic_strain, limits, Some(prescribed),
            prescribed_temperature_derivative, checkpoint,
        )
    }

    /// Allocation-free T^T J_T delta_temperature in independent master rows.
    /// Owner, scales and finite dimensions are admitted before output writes.
    /// Cancellation/overflow during accumulation can leave partial output.
    pub fn thermal_residual_temperature_apply_into(
        &self,
        prepared: &PreparedThermalResidualJacobian3,
        temperature_direction: &[f64],
        output: &mut [f64],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        prepared.apply_mapped(
            &self.raw,
            ThermalCoordinateMap3 { rows: Some(&self.rows), fixed: &self.fixed },
            temperature_direction, output, false, checkpoint,
        )
    }

    /// Allocation-free J_T^T T dual using the same cached forward coefficients.
    /// Neither algebraic directions nor duals consume physical strain budgets.
    pub fn thermal_residual_temperature_transpose_into(
        &self,
        prepared: &PreparedThermalResidualJacobian3,
        dual: &[f64],
        output: &mut [f64],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        prepared.apply_mapped(
            &self.raw,
            ThermalCoordinateMap3 { rows: Some(&self.rows), fixed: &self.fixed },
            dual, output, true, checkpoint,
        )
    }
}
