//! Individual flow controls after the shared-solid interface adjoint closes.
use super::*;
use crate::graph::GraphSolution;

impl CoupledLinearization<'_, '_> {
    /// Immutable hydraulic state used by this bound solid/air fixed point.
    #[must_use]
    pub fn hydraulics(&self) -> &GraphSolution { self.air.hydraulics() }

    /// Return the original thermal gradient and total signed-flow derivatives.
    ///
    /// The second tuple entry is in original graph branch order, in objective
    /// units/(m^3/s). It includes solid feedback, changing junction mixing,
    /// exchanger capacities and external supply/exhaust exchange at fixed h.
    /// Flow-derived convection must be chained through the returned log_htc.
    /// This does not include hydraulic resistance/pressure sensitivities yet;
    /// contract the returned weights with an accepted hydraulic adjoint.
    ///
    /// One coupled adjoint plus one air-only reverse sweep is used, independent
    /// of the number of branches. The individual-flow producer refuses stagnant
    /// branches, reservoir sign switches and imposed return feedback.
    pub fn pullback_branch_flows(
        &self, cx: &Cx<'_>, objective: &CoupledObjective, config: InterfaceSolveConfig,
    ) -> Result<(CoupledGradient, Vec<f64>)> {
        let thermal = self.pullback(cx, objective, config)?;
        self.with_branch_flows(cx, objective, thermal)
    }

    /// The same total flow gradient with checked IQN-ILS interface acceleration.
    /// No primal iteration history or per-control solid solves are differentiated.
    pub fn pullback_branch_flows_iqn(
        &self, cx: &Cx<'_>, objective: &CoupledObjective, config: InterfaceSolveConfig,
        acceleration: IqnIlsConfig,
    ) -> Result<(CoupledGradient, Vec<f64>)> {
        let thermal = self.pullback_iqn(cx, objective, config, acceleration)?;
        self.with_branch_flows(cx, objective, thermal)
    }

    fn with_branch_flows(
        &self, cx: &Cx<'_>, objective: &CoupledObjective, thermal: CoupledGradient,
    ) -> Result<(CoupledGradient, Vec<f64>)> {
        poll(cx)?;
        let mut weights = objective.air.clone();
        if weights.references.len() != thermal.interface_adjoint.len() {
            return Err(CoupledSensitivityError::InvalidInput("flow pullback interface order mismatch"));
        }
        for (weight, multiplier) in weights.references.iter_mut().zip(&thermal.interface_adjoint) {
            poll(cx)?;
            *weight = finite(*weight + multiplier)?;
        }
        let flows = self.air.pullback_branch_flows(cx, &weights)?;
        poll(cx)?;
        Ok((thermal, flows))
    }
}

#[cfg(test)]
mod tests;
