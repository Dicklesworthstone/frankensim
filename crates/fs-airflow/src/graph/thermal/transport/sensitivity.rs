//! Matrix-free tangents and adjoints of the admitted thermal transport DAG.
//!
//! The controls are wall temperatures, external supply temperatures, and
//! `ln(U)` for each segment conductance `U = h A`. One forward sweep applies
//! the Jacobian; one reverse sweep applies its transpose, independent of the
//! number of controls. Neither sweep reruns the primal march or a solid solve.
//!
//! With `x = U/C`, `a = exp(-x)`, `e = 1-a`, `g = e/x`, and `D = Tw-Tin`,
//! the segment derivatives with respect to `(Tin, Tw, ln(U))` are
//! `Tout: (a, e, x*a*D)`, `Tref: (g, 1-g, (g-a)*D)`, and
//! `Q: (-C*e, C*e, U*a*D)`. The two small differences in the reference row
//! use series at small NTU instead of subtracting nearly equal numbers.
//!
//! These differentiate the smooth nominal model, not floating-point rounding,
//! admission predicates, or the heat-budget test. Flows, fluid properties and
//! flow directions stay fixed. A geometry change that also changes hydraulic
//! resistance needs a hydraulic derivative; shared-solid feedback needs the
//! solid Jacobian and an implicit coupled solve. Neither is silently included.

use fs_exec::Cx;
use fs_math::det;

use super::{BranchThermalModel, TransportError, TransportMarch, TransportNetwork, finite, poll, positive};

/// An input perturbation. Temperature entries are differences, not absolute K.
#[derive(Debug, Clone, PartialEq)]
pub struct TransportDirection {
    /// One wall-temperature perturbation per `TransportNetwork::regions()` row.
    pub walls_k: Vec<f64>,
    /// One supply-temperature perturbation per hydraulic node; zero off supplies.
    pub inlets_k: Vec<f64>,
    /// One relative conductance perturbation `dU/U` per region, in row order.
    pub log_conductances: Vec<f64>,
}

/// Weights defining a scalar linear functional of the transported outputs.
/// Units belong to the functional; for a temperature objective its temperature
/// weights are dimensionless and its watt weights are K/W.
#[derive(Debug, Clone, PartialEq)]
pub struct TransportObjective {
    /// One weight per node; stagnant/unknown nodes require zero.
    pub node_temperatures: Vec<f64>,
    /// One weight per branch outlet in hydraulic graph order; stagnant = zero.
    pub branch_outlets: Vec<f64>,
    /// One weight per effective Robin reference, in region order.
    pub references: Vec<f64>,
    /// Weight on total heat leaving the walls.
    pub wall_heat_rate: f64,
    /// Weight on external sensible-enthalpy gain.
    pub external_heat_gain: f64,
    /// Weight on external gain minus wall heat; never a corrected balance.
    pub heat_imbalance: f64,
    /// Weight on the separately disclosed hydraulic energy defect.
    pub hydraulic_energy_defect: f64,
}

/// Jacobian-vector product. `None` keeps the primal's unknown stagnant state.
#[derive(Debug, Clone, PartialEq)]
pub struct TransportDifferential {
    /// Node-temperature changes, K.
    pub node_temperatures_k: Vec<Option<f64>>,
    /// Branch-outlet temperature changes, K.
    pub branch_outlets_k: Vec<Option<f64>>,
    /// Effective Robin-reference changes, K.
    pub reference_temperatures_k: Vec<f64>,
    /// Change of total heat leaving the walls, W.
    pub wall_heat_rate_w: f64,
    /// Change of external sensible-enthalpy gain, W.
    pub external_heat_gain_w: f64,
    /// Change of the raw heat imbalance, W.
    pub heat_imbalance_w: f64,
    /// Change of the disclosed hydraulic energy defect, W.
    pub hydraulic_energy_defect_w: f64,
}

/// Transposed-Jacobian product for one scalar objective.
#[derive(Debug, Clone, PartialEq)]
pub struct TransportGradient {
    /// Objective units per kelvin of wall temperature, in region order.
    pub walls: Vec<f64>,
    /// Objective units per kelvin of supply temperature; zero off supplies.
    pub inlets: Vec<f64>,
    /// Objective units per unit change in `ln(h A)`, in region order.
    /// Divide by the nominal conductance for a derivative with respect to U.
    pub log_conductances: Vec<f64>,
}

#[derive(Debug, Clone, Copy)]
struct SegmentJacobian {
    outlet: [f64; 3],
    reference: [f64; 3],
    heat: [f64; 3],
}

#[derive(Debug, Default)]
struct MixingRow {
    supply: f64,
    incoming: Vec<(usize, f64)>,
}

/// A linearization bound to an immutable network and an actually admitted march.
/// Storage and every tangent/adjoint sweep are linear in nodes, edges and
/// segments. Publicly assembled or mismatched `TransportMarch` values cannot
/// be substituted for the private primal used to construct these derivatives.
#[derive(Debug)]
pub struct TransportLinearization<'network, 'flow> {
    network: &'network TransportNetwork<'flow>,
    primal: TransportMarch,
    segments: Vec<SegmentJacobian>,
    mixers: Vec<MixingRow>,
    upstream: Vec<usize>,
    reference_node: Option<usize>,
}

impl<'flow> TransportNetwork<'flow> {
    /// Run and admit one primal march, then prepare reusable tangent/adjoint
    /// sweeps. The derivatives hold hydraulics fixed and do not include the
    /// response of a shared solid. Nonfinite derivative arithmetic refuses.
    pub fn linearize<'network>(
        &'network self, cx: &Cx<'_>, walls_k: &[f64],
    ) -> Result<TransportLinearization<'network, 'flow>, TransportError> {
        let primal = self.march(cx, walls_k)?;
        let mut segments = Vec::with_capacity(walls_k.len());
        for (branch, model) in self.models.iter().enumerate() {
            poll(cx)?;
            if let BranchThermalModel::Exchange(rows) = model {
                let march = primal.branches[branch].march.as_ref()
                    .ok_or(TransportError::InvalidInput("missing primal exchanger march"))?;
                for (index, (row, state)) in rows.iter().zip(&march.segments).enumerate() {
                    poll(cx)?;
                    let driving = finite(walls_k[self.offsets[branch] + index]
                        - state.inlet_temperature_k, "segment linearization temperature difference")?;
                    segments.push(segment_jacobian(state.ntu, state.effectiveness,
                        row.area_m2() * row.htc_w_per_m2_k(), self.capacities[branch], driving)?);
                }
            }
        }
        let mut mixers: Vec<MixingRow> = (0..self.order.len()).map(|_| MixingRow::default()).collect();
        let mut totals: Vec<f64> = self.external_capacity.iter().map(|c| c.max(0.0)).collect();
        let mut upstream = vec![0; self.models.len()];
        // This is the primal's arrival order, including the external supply
        // before any arriving branch. Normalization uses INCOMING capacity:
        // an admitted hydraulic residual is not silently repaired by outgoing C.
        for &node in &self.order {
            poll(cx)?;
            for &branch in &self.outgoing[node] {
                poll(cx)?;
                upstream[branch] = node;
                let down = self.downstream[branch];
                totals[down] = positive(totals[down] + self.capacities[branch], "mixing derivative capacity")?;
                mixers[down].incoming.push((branch, self.capacities[branch]));
            }
        }
        for (node, row) in mixers.iter_mut().enumerate() {
            poll(cx)?;
            if totals[node] > 0.0 {
                row.supply = self.external_capacity[node].max(0.0) / totals[node];
                for (_, weight) in &mut row.incoming { poll(cx)?; *weight /= totals[node]; }
            }
        }
        poll(cx)?;
        Ok(TransportLinearization { network: self, primal, segments, mixers, upstream,
            reference_node: self.inlet_temperatures.iter().position(Option::is_some) })
    }
}

impl TransportLinearization<'_, '_> {
    /// The admitted primal to which every derivative is bound.
    #[must_use]
    pub const fn primal(&self) -> &TransportMarch { &self.primal }

    /// Correctly sized zero perturbations; only declared supply slots may vary.
    #[must_use]
    pub fn zero_direction(&self) -> TransportDirection {
        TransportDirection { walls_k: vec![0.0; self.segments.len()],
            inlets_k: vec![0.0; self.mixers.len()], log_conductances: vec![0.0; self.segments.len()] }
    }

    /// Correctly sized zero objective weights.
    #[must_use]
    pub fn zero_objective(&self) -> TransportObjective {
        TransportObjective { node_temperatures: vec![0.0; self.mixers.len()],
            branch_outlets: vec![0.0; self.upstream.len()], references: vec![0.0; self.segments.len()],
            wall_heat_rate: 0.0, external_heat_gain: 0.0, heat_imbalance: 0.0, hydraulic_energy_defect: 0.0 }
    }

    /// Apply the full transport Jacobian in one topological pass.
    pub fn apply(&self, cx: &Cx<'_>, direction: &TransportDirection) -> Result<TransportDifferential, TransportError> {
        self.check_direction(cx, direction)?;
        let mut result = TransportDifferential {
            node_temperatures_k: vec![None; self.mixers.len()], branch_outlets_k: vec![None; self.upstream.len()],
            reference_temperatures_k: vec![0.0; self.segments.len()], wall_heat_rate_w: 0.0,
            external_heat_gain_w: 0.0, heat_imbalance_w: 0.0, hydraulic_energy_defect_w: 0.0,
        };
        let reference = self.reference_node.map_or(0.0, |node| direction.inlets_k[node]);
        for &node in &self.network.order {
            poll(cx)?;
            if self.primal.node_temperatures_k[node].is_none() { continue; }
            let mixer = &self.mixers[node];
            let mut temperature = finite(mixer.supply * direction.inlets_k[node], "supply tangent")?;
            for &(branch, weight) in &mixer.incoming {
                poll(cx)?;
                let upstream = result.branch_outlets_k[branch]
                    .ok_or(TransportError::InvalidInput("missing upstream tangent"))?;
                add(&mut temperature, weight * upstream)?;
            }
            result.node_temperatures_k[node] = Some(temperature);
            let external = self.network.external_capacity[node];
            let entering = if external > 0.0 { direction.inlets_k[node] } else { temperature };
            add(&mut result.external_heat_gain_w, -external * (entering - reference))?;
            add(&mut result.hydraulic_energy_defect_w,
                self.network.capacity_residual[node] * (temperature - reference))?;
            for &branch in &self.network.outgoing[node] {
                poll(cx)?;
                let mut inlet = temperature;
                for row in self.network.row_range(branch) {
                    poll(cx)?;
                    let input = [inlet, direction.walls_k[row], direction.log_conductances[row]];
                    let kernel = self.segments[row];
                    result.reference_temperatures_k[row] = dot(kernel.reference, input)?;
                    add(&mut result.wall_heat_rate_w, dot(kernel.heat, input)?)?;
                    inlet = dot(kernel.outlet, input)?;
                }
                result.branch_outlets_k[branch] = Some(inlet);
            }
        }
        result.heat_imbalance_w = finite(result.external_heat_gain_w - result.wall_heat_rate_w, "heat imbalance tangent")?;
        poll(cx)?;
        Ok(result)
    }

    /// Differentiate one weighted output with respect to ALL controls in one
    /// reverse pass. No dense Jacobian or per-control primal runs are needed.
    pub fn pullback(&self, cx: &Cx<'_>, objective: &TransportObjective) -> Result<TransportGradient, TransportError> {
        self.check_objective(cx, objective)?;
        let mut gradient = TransportGradient { walls: vec![0.0; self.segments.len()],
            inlets: vec![0.0; self.mixers.len()], log_conductances: vec![0.0; self.segments.len()] };
        let mut nodes = objective.node_temperatures.clone();
        let heat_weight = finite(objective.wall_heat_rate - objective.heat_imbalance, "wall heat adjoint weight")?;
        let external_weight = finite(objective.external_heat_gain + objective.heat_imbalance, "external heat adjoint weight")?;
        let mut reference_weight = 0.0;
        for &node in &self.network.order {
            poll(cx)?;
            if self.primal.node_temperatures_k[node].is_none() { continue; }
            let external = self.network.external_capacity[node];
            let residual = self.network.capacity_residual[node];
            if external > 0.0 { add(&mut gradient.inlets[node], -external * external_weight)?; }
            else { add(&mut nodes[node], -external * external_weight)?; }
            add(&mut nodes[node], residual * objective.hydraulic_energy_defect)?;
            add(&mut reference_weight, external * external_weight - residual * objective.hydraulic_energy_defect)?;
        }
        if let Some(node) = self.reference_node { add(&mut gradient.inlets[node], reference_weight)?; }
        for &node in self.network.order.iter().rev() {
            poll(cx)?;
            if self.primal.node_temperatures_k[node].is_none() { continue; }
            let mixer = &self.mixers[node];
            add(&mut gradient.inlets[node], mixer.supply * nodes[node])?;
            for &(branch, fraction) in mixer.incoming.iter().rev() {
                poll(cx)?;
                let mut outlet_weight = finite(objective.branch_outlets[branch] + fraction * nodes[node], "branch adjoint")?;
                for row in self.network.row_range(branch).rev() {
                    poll(cx)?;
                    let kernel = self.segments[row];
                    let weights = [outlet_weight, objective.references[row], heat_weight];
                    gradient.walls[row] = dot([kernel.outlet[1], kernel.reference[1], kernel.heat[1]], weights)?;
                    gradient.log_conductances[row] = dot([kernel.outlet[2], kernel.reference[2], kernel.heat[2]], weights)?;
                    outlet_weight = dot([kernel.outlet[0], kernel.reference[0], kernel.heat[0]], weights)?;
                }
                add(&mut nodes[self.upstream[branch]], outlet_weight)?;
            }
        }
        poll(cx)?;
        Ok(gradient)
    }

    fn check_direction(&self, cx: &Cx<'_>, direction: &TransportDirection) -> Result<(), TransportError> {
        check_vector(cx, &direction.walls_k, self.segments.len())?;
        check_vector(cx, &direction.log_conductances, self.segments.len())?;
        check_vector(cx, &direction.inlets_k, self.mixers.len())?;
        for (node, &delta) in direction.inlets_k.iter().enumerate() {
            poll(cx)?;
            if self.network.inlet_temperatures[node].is_none() && delta != 0.0 {
                return Err(TransportError::Node { node, reason: "only declared external supplies have inlet derivatives" });
            }
        }
        Ok(())
    }

    fn check_objective(&self, cx: &Cx<'_>, objective: &TransportObjective) -> Result<(), TransportError> {
        check_vector(cx, &objective.node_temperatures, self.mixers.len())?;
        check_vector(cx, &objective.branch_outlets, self.upstream.len())?;
        check_vector(cx, &objective.references, self.segments.len())?;
        check_vector(cx, &[objective.wall_heat_rate, objective.external_heat_gain,
            objective.heat_imbalance, objective.hydraulic_energy_defect], 4)?;
        for (node, &weight) in objective.node_temperatures.iter().enumerate() {
            poll(cx)?;
            if self.primal.node_temperatures_k[node].is_none() && weight != 0.0 {
                return Err(TransportError::Node { node, reason: "unknown stagnant temperature has no objective derivative" });
            }
        }
        for (branch, &weight) in objective.branch_outlets.iter().enumerate() {
            poll(cx)?;
            if self.primal.branches[branch].outlet_temperature_k.is_none() && weight != 0.0 {
                return Err(TransportError::Branch { branch, reason: "unknown stagnant outlet has no objective derivative" });
            }
        }
        Ok(())
    }
}

fn segment_jacobian(x: f64, effectiveness: f64, conductance: f64, capacity: f64, driving: f64)
    -> Result<SegmentJacobian, TransportError>
{
    let attenuation = det::exp(-x);
    let g = effectiveness / x;
    let (one_minus_g, g_minus_a) = if x < 1.0e-3 {
        (x * (0.5 + x * (-1.0/6.0 + x * (1.0/24.0 + x * (-1.0/120.0 + x * (1.0/720.0 - x/5040.0))))),
         x * (0.5 + x * (-1.0/3.0 + x * (1.0/8.0 + x * (-1.0/30.0 + x * (1.0/144.0 - x/840.0))))))
    } else { (1.0 - g, g - attenuation) };
    let heat_wall = capacity * effectiveness;
    let kernel = SegmentJacobian {
        outlet: [attenuation, effectiveness, (x * attenuation) * driving],
        reference: [g, one_minus_g, g_minus_a * driving],
        heat: [-heat_wall, heat_wall, (conductance * attenuation) * driving],
    };
    for value in kernel.outlet.into_iter().chain(kernel.reference).chain(kernel.heat) {
        finite(value, "nonfinite segment derivative")?;
    }
    Ok(kernel)
}

fn check_vector(cx: &Cx<'_>, values: &[f64], expected: usize) -> Result<(), TransportError> {
    poll(cx)?;
    if values.len() != expected { return Err(TransportError::InvalidInput("transport derivative shape mismatch")); }
    for &value in values { poll(cx)?; finite(value, "nonfinite transport derivative input")?; }
    Ok(())
}

fn dot(coefficients: [f64; 3], values: [f64; 3]) -> Result<f64, TransportError> {
    finite(coefficients[0] * values[0] + coefficients[1] * values[1]
        + coefficients[2] * values[2], "nonfinite transport derivative product")
}

fn add(total: &mut f64, value: f64) -> Result<(), TransportError> {
    *total = finite(*total + value, "nonfinite transport derivative accumulation")?;
    Ok(())
}

#[cfg(test)]
mod tests;
