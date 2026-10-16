//! Reverse derivatives of individual signed flows, at fixed flow directions.
//!
//! Capacity changes alter both NTU and junction mixing fractions. Reservoir
//! exchange is recomputed from incident flows, as in TransportNetwork::new;
//! free-node imbalance is retained rather than differentiated as exactly zero.
//! No hydraulic solve is performed here. Contract these weights with an
//! accepted hydraulic tangent/adjoint to obtain physical network controls.

use super::*;
use crate::graph::GraphSolution;

impl TransportLinearization<'_, '_> {
    /// Exact immutable hydraulic state used by this transport linearization.
    #[must_use]
    pub fn hydraulics(&self) -> &GraphSolution { self.network.hydraulics() }

    /// Scalar objective derivatives per signed branch flow (objective/(m^3/s)).
    ///
    /// Output order and signs match the original hydraulic graph. This includes
    /// exchanger capacity, changing junction fractions, external enthalpy and
    /// the disclosed free-node energy defect. Wall/supply temperatures, fluid
    /// properties and hA stay fixed. Add any flow-derived h response separately.
    ///
    /// This differentiates the smooth extension off the conservation manifold;
    /// the physical direction must come from a hydraulic solve. A single reverse
    /// DAG sweep suffices for all branches. Zero flows, zero reservoir exchange
    /// and imposed return feedback refuse, not silently frozen derivatives.
    pub fn pullback_branch_flows(
        &self, cx: &Cx<'_>, objective: &TransportObjective,
    ) -> Result<Vec<f64>, TransportError> {
        self.check_objective(cx, objective)?;
        if self.network.feedback.is_some() {
            return Err(TransportError::InvalidInput("individual flow derivatives do not freeze imposed return feedback"));
        }
        let n = self.mixers.len();
        let mut fixed = vec![false; n];
        for boundary in &self.network.flow.boundaries {
            poll(cx)?;
            fixed[boundary.node] = true;
            if self.network.external_capacity[boundary.node] == 0.0 {
                return Err(TransportError::Node { node: boundary.node,
                    reason: "zero reservoir exchange has no admitted two-sided flow derivative" });
            }
        }
        for (branch, state) in self.network.flow.branches.iter().enumerate() {
            poll(cx)?;
            if state.flow.value() == 0.0 {
                return Err(TransportError::Branch { branch,
                    reason: "zero flow has no admitted individual flow derivative" });
            }
        }
        let heat_weight = finite(objective.wall_heat_rate - objective.heat_imbalance, "flow heat weight")?;
        let external_weight = finite(objective.external_heat_gain + objective.heat_imbalance, "flow external weight")?;
        let mut nodes = objective.node_temperatures.clone();
        // Derivatives w.r.t. outgoing-minus-incoming capacity at each node.
        let mut balances = vec![0.0; n];
        let mut logarithmic = vec![0.0; self.upstream.len()];
        let reference = self.primal.enthalpy_reference_k;
        for &node in &self.network.order {
            poll(cx)?;
            let Some(temperature) = self.primal.node_temperatures_k[node] else { continue; };
            let external = self.network.external_capacity[node];
            if external <= 0.0 { add(&mut nodes[node], -external * external_weight)?; }
            add(&mut nodes[node], self.network.capacity_residual[node] * objective.hydraulic_energy_defect)?;
            balances[node] = if fixed[node] {
                let entering = self.network.inlet_temperatures[node].unwrap_or(temperature);
                finite(-external_weight * (entering - reference), "reservoir capacity derivative")?
            } else {
                finite(objective.hydraulic_energy_defect * (temperature - reference), "residual capacity derivative")?
            };
        }
        for &node in self.network.order.iter().rev() {
            poll(cx)?;
            let Some(temperature) = self.primal.node_temperatures_k[node] else { continue; };
            let mixer = &self.mixers[node];
            let mut total = self.network.external_capacity[node].max(0.0);
            for &(branch, _) in &mixer.incoming {
                poll(cx)?;
                add(&mut total, self.network.capacities[branch])?;
            }
            positive(total, "flow derivative incoming mixing capacity")?;
            if let Some(supply) = self.network.inlet_temperatures[node] {
                add(&mut balances[node], (nodes[node] / total) * (supply - temperature))?;
            }
            for &(branch, fraction) in mixer.incoming.iter().rev() {
                poll(cx)?;
                let outlet = self.primal.branches[branch].outlet_temperature_k
                    .ok_or(TransportError::InvalidInput("missing flow derivative branch outlet"))?;
                // dT_mix/dln(C_e) = fraction_e * (T_e - T_mix).
                add(&mut logarithmic[branch], nodes[node] * fraction * (outlet - temperature))?;
                let mut weight = finite(objective.branch_outlets[branch] + fraction * nodes[node], "flow outlet adjoint")?;
                if let BranchThermalModel::Exchange(rows) = &self.network.models[branch] {
                    for (row, segment) in self.network.row_range(branch).zip(rows).rev() {
                        poll(cx)?;
                        let kernel = self.segments[row];
                        let conductance = positive(segment.area_m2() * segment.htc_w_per_m2_k(), "flow segment conductance")?;
                        // Temperature rows depend on U/C. For heat, dQ/dln(C)
                        // = U*(g-a)*D, using the existing small-NTU stable row.
                        add(&mut logarithmic[branch], dot(
                            [-kernel.outlet[2], -kernel.reference[2], conductance * kernel.reference[2]],
                            [weight, objective.references[row], heat_weight],
                        )?)?;
                        weight = dot([kernel.outlet[0], kernel.reference[0], kernel.heat[0]],
                            [weight, objective.references[row], heat_weight])?;
                    }
                }
                add(&mut nodes[self.upstream[branch]], weight)?;
            }
        }
        let mut result = Vec::with_capacity(logarithmic.len());
        for (branch, value) in logarithmic.iter_mut().enumerate() {
            poll(cx)?;
            add(value, self.network.capacities[branch]
                * (balances[self.upstream[branch]] - balances[self.network.downstream[branch]]))?;
            result.push(finite(*value / self.network.flow.branches[branch].flow.value(), "signed branch flow gradient")?);
        }
        poll(cx)?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
