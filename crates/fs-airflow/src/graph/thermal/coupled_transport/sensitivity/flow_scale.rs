//! Common flow/capacity-scale derivative of the complete coupled problem.
//!
//! Scale ALL branch and external heat-capacity rates by s, holding their
//! proportions, signs, topology, inlet temperatures and hA fixed. Scaling hA
//! and capacity together leaves every transported temperature unchanged and
//! multiplies every transported heat rate by s. Euler's identity therefore
//! gives dJ/dln(s) = J_heat - sum_i dJ/dln(hA_i) for the AIR pullback.
//!
//! The weights must include the converged interface multiplier; using only
//! the explicit air objective would miss the shared solid's response. The
//! solid's log(h) gradient must NOT be subtracted in this identity: its
//! boundary operator does not scale with air capacity.

use super::*;

/// Existing thermal controls plus a common signed-flow scaling control.
#[derive(Debug, Clone)]
pub struct CoupledFlowScaleGradient {
    /// Total coupled derivatives with respect to the original controls.
    pub thermal: CoupledGradient,
    /// Objective units per ln(s) when all branch/external capacity rates scale
    /// by s at fixed hA. This is not an individual branch-flow derivative.
    pub log_flow_scale: f64,
}

impl CoupledLinearization<'_, '_> {
    /// Extend the stationary coupled adjoint with a common flow-scale control.
    /// One additional air-only reverse sweep is performed, not a solid solve.
    ///
    /// # Errors
    /// Preserves the original adjoint's admission, residual and budget checks;
    /// nonfinite contraction arithmetic and cancellation also refuse.
    pub fn pullback_flow_scale(
        &self,
        cx: &Cx<'_>,
        objective: &CoupledObjective,
        config: InterfaceSolveConfig,
    ) -> Result<CoupledFlowScaleGradient> {
        let thermal = self.pullback(cx, objective, config)?;
        self.with_flow_scale(cx, objective, thermal)
    }

    /// Extend the accelerated coupled adjoint with a common flow-scale control.
    /// Each solve starts fresh IQN history and checks the true transpose
    /// equation before this contraction. Geometry, material laws, loads and
    /// fixed-temperature boundaries retain the original binding.
    ///
    /// # Errors
    /// The refusals of [`Self::pullback_iqn`] and nonfinite/cancelled contraction.
    pub fn pullback_flow_scale_iqn(
        &self,
        cx: &Cx<'_>,
        objective: &CoupledObjective,
        config: InterfaceSolveConfig,
        acceleration: IqnIlsConfig,
    ) -> Result<CoupledFlowScaleGradient> {
        let thermal = self.pullback_iqn(cx, objective, config, acceleration)?;
        self.with_flow_scale(cx, objective, thermal)
    }

    fn with_flow_scale(
        &self,
        cx: &Cx<'_>,
        objective: &CoupledObjective,
        thermal: CoupledGradient,
    ) -> Result<CoupledFlowScaleGradient> {
        poll(cx)?;
        let mut weights = objective.air.clone();
        for (weight, multiplier) in weights.references.iter_mut().zip(&thermal.interface_adjoint) {
            *weight = finite(*weight + multiplier)?;
        }
        let air = self.air.pullback(cx, &weights)?;
        let primal = self.air.primal();
        // Temperature functionals have degree zero; watt functionals have
        // degree one. Keep even the raw imbalance/defect terms rather than
        // silently treating an admitted numerical residual as exactly zero.
        let mut log_flow_scale = 0.0;
        for (weight, value) in [
            (weights.wall_heat_rate, primal.wall_heat_rate_w),
            (weights.external_heat_gain, primal.external_heat_gain_w),
            (weights.heat_imbalance, primal.heat_imbalance_w),
            (weights.hydraulic_energy_defect, primal.hydraulic_energy_defect_w),
        ] {
            log_flow_scale = finite(log_flow_scale + finite(weight * value)?)?;
        }
        for value in air.log_conductances {
            poll(cx)?;
            log_flow_scale = finite(log_flow_scale - value)?;
        }
        poll(cx)?;
        Ok(CoupledFlowScaleGradient { thermal, log_flow_scale })
    }
}
