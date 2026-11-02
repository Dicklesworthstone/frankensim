//! Simultaneous mixing for CLOSED, thermally anchored hydraulic circulation.
//!
//! Each existing branch is affine: Tout = a Tin + b(Twall), 0 <= a <= 1.
//! At each node i, T_i - sum_{e:j->i}(w_e a_e T_j) = sum_{e:j->i}(w_e b_e),
//! with incoming-capacity mixing weights. A heat-exchanging branch anchors its downstream node;
//! reachability then establishes uniqueness in each active component. These
//! are frozen-flow sensible-heat equations, not buoyancy/temperature feedback.

use super::*;
use super::recirculation::solve_feedback;

/// Hard bound on simultaneous native circulation temperatures.
pub const MAX_CIRCULATION_NODES: usize = 64;

#[derive(Debug)]
pub(super) struct Circulation {
    matrix: Vec<Vec<f64>>,
    pub(super) incoming_capacities: Vec<f64>,
}

impl Circulation {
    pub(super) fn new(network: &TransportNetwork<'_>, cx: &Cx<'_>)
        -> Result<Self, TransportError> {
        let n = network.order.len();
        let mut incoming_capacities = vec![0.0; n];
        for &node in &network.order {
            for &branch in &network.outgoing[node] {
                poll(cx)?;
                let down = network.downstream[branch];
                incoming_capacities[down] = finite(incoming_capacities[down]
                    + network.capacities[branch], "circulation incoming capacity")?;
            }
        }
        let mut matrix = vec![vec![0.0; n]; n];
        let mut anchored = vec![false; n];
        let mut ready = BTreeSet::new();
        for node in 0..n {
            poll(cx)?;
            matrix[node][node] = 1.0;
            if incoming_capacities[node] == 0.0 && !network.outgoing[node].is_empty() {
                return Err(TransportError::Node { node,
                    reason: "closed circulation has outgoing flow without an incoming stream" });
            }
            for &branch in &network.outgoing[node] {
                poll(cx)?;
                let down = network.downstream[branch];
                let weight = network.capacities[branch] / incoming_capacities[down];
                let mut attenuation = 1.0;
                if let BranchThermalModel::Exchange(rows) = &network.models[branch] {
                    for row in rows {
                        poll(cx)?;
                        attenuation *= fs_math::det::exp(-row.ntu(network.capacities[branch]));
                    }
                }
                matrix[down][node] = finite(matrix[down][node] - weight * attenuation,
                    "circulation mixing coefficient")?;
                if attenuation < 1.0 && !anchored[down] {
                    anchored[down] = true;
                    ready.insert(down);
                }
            }
        }
        while let Some(node) = ready.pop_first() {
            poll(cx)?;
            for &branch in &network.outgoing[node] {
                poll(cx)?;
                let down = network.downstream[branch];
                if !anchored[down] {
                    anchored[down] = true;
                    ready.insert(down);
                }
            }
        }
        for node in 0..n {
            if incoming_capacities[node] > 0.0 && !anchored[node] {
                return Err(TransportError::Node { node,
                    reason: "closed circulation needs a numerically resolved wall-heat anchor in every active component" });
            }
        }
        // Reuse the existing bounded elimination and backward-residual gate.
        // This also refuses an anchor too weak to resolve the linear system.
        solve_feedback(cx, &matrix, &vec![0.0; n], false)?;
        Ok(Self { matrix, incoming_capacities })
    }

    pub(super) fn solve(&self, cx: &Cx<'_>, rhs: &[f64], transpose: bool)
        -> Result<Vec<f64>, TransportError> {
        solve_feedback(cx, &self.matrix, rhs, transpose)
    }

    pub(super) fn march(&self, network: &TransportNetwork<'_>, cx: &Cx<'_>, walls: &[f64])
        -> Result<TransportMarch, TransportError> {
        // Center at an actual wall. A uniform-wall equilibrium then has an
        // exactly zero RHS even when all heat anchors are weak. The declared
        // network reference still owns the separate energy-accounting chart.
        let reference = walls.first().copied().unwrap_or(network.reference_k);
        let mut rhs = vec![0.0; network.order.len()];
        for (branch, model) in network.models.iter().enumerate() {
            poll(cx)?;
            if let BranchThermalModel::Exchange(_) = model {
                // The original owner computes every NTU/effectiveness and
                // wall contribution, including its small-NTU arithmetic.
                let local_walls = &walls[network.row_range(branch)];
                let baseline = network.path(branch, reference)?.march(local_walls)?;
                let mut forcing = 0.0;
                for (segment, &wall) in baseline.segments.iter().zip(local_walls) {
                    poll(cx)?;
                    forcing = finite(forcing + ((wall - reference) - forcing)
                        * segment.effectiveness, "circulation branch forcing")?;
                }
                let down = network.downstream[branch];
                let weight = network.capacities[branch] / self.incoming_capacities[down];
                rhs[down] = finite(rhs[down] + weight * forcing,
                    "circulation wall-temperature forcing")?;
            }
        }
        let mut temperatures = self.solve(cx, &rhs, false)?;
        for temperature in &mut temperatures {
            poll(cx)?;
            *temperature = positive(reference + *temperature, "circulating node temperature")?;
        }
        // Re-march all physical branches and independently recompute every
        // mixer and the enclosing heat balance at the accepted node values.
        network.march_once_with_nodes(cx, Some(walls), &network.inlet_temperatures,
            Some(&temperatures))
    }
}

#[cfg(test)]
mod tests;
