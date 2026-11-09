//! Physical entropy coupling through the original adaptive displacement map.

use std::ops::ControlFlow;

use crate::elastic3::ElasticityError3;
use crate::elastic3::adaptive::AdaptiveElasticity3;
use crate::elastic3::thermal::{
    PreparedThermoelasticEntropy3, ThermalCoordinateMap3, ThermoelasticEntropyLimits3,
};

impl AdaptiveElasticity3 {
    /// Integrate native bulk source entropy and its displacement derivative on
    /// original physical leaves. No global field or second hanging map is built.
    /// Physical admission, explicit limits and cache semantics match the raw
    /// Cartesian preparation.
    pub fn prepare_thermoelastic_entropy(
        &self,
        u: &[f64],
        maximum_kinematic_strain: f64,
        limits: ThermoelasticEntropyLimits3,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<PreparedThermoelasticEntropy3, ElasticityError3> {
        PreparedThermoelasticEntropy3::prepare(
            &self.raw,
            ThermalCoordinateMap3 { rows: Some(&self.rows), fixed: &self.fixed },
            u, maximum_kinematic_strain, limits, checkpoint,
        )
    }

    /// Allocation-free S_u T direction in original cell entropy rows.
    /// Algebraic directions have no physical strain limit. Owner/scales/shape/
    /// finite admission precedes writes; cancellation or overflow can leave
    /// partial output, which must be discarded on error.
    pub fn thermoelastic_entropy_displacement_apply_into(
        &self,
        prepared: &PreparedThermoelasticEntropy3,
        direction: &[f64],
        output: &mut [f64],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        prepared.apply_mapped(
            &self.raw,
            ThermalCoordinateMap3 { rows: Some(&self.rows), fixed: &self.fixed },
            direction, output, false, false, checkpoint,
        )
    }

    /// Allocation-free T^T S_u^T cell cotangent, retaining the exact same
    /// physical column and original support mask as the forward action.
    pub fn thermoelastic_entropy_displacement_transpose_into(
        &self,
        prepared: &PreparedThermoelasticEntropy3,
        entropy_cotangent: &[f64],
        output: &mut [f64],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        prepared.apply_mapped(
            &self.raw,
            ThermalCoordinateMap3 { rows: Some(&self.rows), fixed: &self.fixed },
            entropy_cotangent, output, true, false, checkpoint,
        )
    }

    /// Reversible mechanical heat rate from reconstructed velocity, in watts per cell.
    /// Semantics match the Cartesian heat-rate action: physical bulk source
    /// only, with signed elastic capacity supplied separately by cell_potential.
    /// It adds no damping, ghost/support storage, or thermal integrator.
    pub fn thermoelastic_reversible_heat_rate_into(
        &self,
        prepared: &PreparedThermoelasticEntropy3,
        velocity: &[f64],
        output: &mut [f64],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        prepared.apply_mapped(
            &self.raw,
            ThermalCoordinateMap3 { rows: Some(&self.rows), fixed: &self.fixed },
            velocity, output, false, true, checkpoint,
        )
    }
}
