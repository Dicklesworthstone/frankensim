//! Native thermal-stress objective partials through the original hanging map.

use std::ops::ControlFlow;

use fs_solver::op::LinearOp;

use super::super::AdaptiveElasticity3;
use crate::elastic3::ElasticityError3;
use crate::elastic3::thermal::ThermoelasticStressPullback3;

impl AdaptiveElasticity3 {
    /// Native thermal-stress partials in independent master coordinates.
    ///
    /// Cotangents pair directly with physical stress components in original
    /// bulk-point order and already include any desired integration weights.
    /// Total and elastic strain admission use cell-local reconstruction with
    /// the original hanging constraints and homogeneous master clamps. The
    /// algebraic cotangents have no material-strain ceiling. Cell-scale and
    /// temperature partials follow `leaves()` order; source C'(T) and integrated,
    /// independently oriented expansion are retained exactly.
    ///
    /// The displacement cotangent uses the existing reference stress transpose
    /// and original constraint transpose. No global primal expansion, derivative
    /// matrix, geometry regeneration or finite-difference material action occurs.
    /// These are explicit partials: subtract the full residual temperature VJP
    /// with the equilibrium adjoint to differentiate a solved stress objective.
    pub fn thermoelastic_stress_pullback(
        &self,
        u: &[f64],
        maximum_kinematic_strain: f64,
        derivatives: &[[f64; 6]],
        mut poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<ThermoelasticStressPullback3, ElasticityError3> {
        let direct = self.raw.thermoelastic_stress_partials_mapped(
            u,
            self.n(),
            maximum_kinematic_strain,
            derivatives,
            |cell, field| {
                std::array::from_fn(|i| self.physical_coefficient(field, cell.nodes[i / 3], i % 3))
            },
            &mut poll,
        )?;
        let displacement =
            self.reference_bulk_stress_pullback(&direct.scaled_cotangents, &mut poll)?;
        if poll().is_break() {
            return Err(ElasticityError3::Cancelled);
        }
        Ok(ThermoelasticStressPullback3 {
            displacement,
            scales: direct.scales,
            temperatures: direct.temperatures,
        })
    }
}
