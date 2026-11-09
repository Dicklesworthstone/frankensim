//! Temperature transpose through the original homogeneous hanging-node map.

use std::ops::ControlFlow;

use fs_solver::op::LinearOp;

use crate::elastic3::ElasticityError3;
use crate::elastic3::adaptive::AdaptiveElasticity3;

impl AdaptiveElasticity3 {
    /// `dual^T d(K(T)u-b_free(T)-f)/dT_cell` on the original adaptive space.
    ///
    /// Temperatures are independent native source states per retained leaf;
    /// geometry, reference mass, scales and external force stay fixed. The same
    /// homogeneous-clamp-aware hanging-node reconstruction acts on the physical
    /// field and the algebraic dual. No raw global displacement copy is made.
    /// Bulk, selected Nitsche and both ghost-side temperature derivatives all
    /// use the original retained geometry. Nonzero imposed displacement is
    /// outside this contract. Physical-domain, smooth-stabilization-branch,
    /// point-budget and cancellation rules are those of
    /// [`crate::elastic3::CutElasticity3::thermal_residual_temperature_pullback`].
    pub fn thermal_residual_temperature_pullback(
        &self,
        u: &[f64],
        dual: &[f64],
        maximum_kinematic_strain: f64,
        max_points: usize,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        self.raw.thermal_residual_temperature_pullback_mapped(
            u,
            dual,
            self.n(),
            maximum_kinematic_strain,
            max_points,
            |node, component, field| self.physical_coefficient(field, node, component),
            checkpoint,
        )
    }
}
