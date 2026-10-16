//! Prescribed return-air fractions around a passive pressure network.
//!
//! A return link models an external, adiabatic pressure-reset path, NOT an
//! extra passive hydraulic edge. Fractions refer to the receiving supply's
//! capacity flow. Fresh makeup and the undrawn exhaust close the outer mass
//! and sensible-heat balance. No fan heat, return-duct loss, residence time,
//! humidity, or solved recirculation hydraulics are inferred.
//!
//! At frozen walls and hA, the existing transport is affine in inlet
//! temperatures. Eliminate its DAG to the bounded supply system (I - R H)x=b;
//! re-march the actual solved supplies and gate every mixer and outer heat
//! balance. Strictly positive fresh makeup avoids an unanchored pure-return
//! loop. This is a nominal floating-point solve, not an interval certificate.

use super::*;

/// Maximum independently mixed external supplies in one feedback solve.
pub const MAX_RECIRCULATION_SUPPLIES: usize = 64;
/// Maximum return links, admitted before allocating a feedback matrix.
pub const MAX_RECIRCULATION_LINKS: usize = 256;

/// An imposed adiabatic return from an external sink to an external supply.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RecirculationLink {
    /// Existing pressure-boundary node with positive external supply flow.
    pub supply_node: usize,
    /// Existing pressure-boundary node with net external exhaust flow.
    pub return_node: usize,
    /// Fraction of the receiving supply, not a fraction of the return exhaust.
    /// Fractions into each supply must sum to less than one. Zero is allowed.
    pub fraction: f64,
}

/// One accepted fresh/return mixture. Capacities are mass flow times cp.
#[derive(Debug, Clone, PartialEq)]
pub struct RecirculatedSupply {
    /// External supply node.
    pub node: usize,
    /// Temperature of the caller-declared fresh makeup, K.
    pub fresh_temperature_k: f64,
    /// Actually marched mixture of fresh and returned air, K.
    pub mixed_temperature_k: f64,
    /// Fresh makeup capacity rate, W/K.
    pub fresh_capacity_w_per_k: f64,
    /// Returned capacity rate entering this supply, W/K.
    pub return_capacity_w_per_k: f64,
}

/// An accepted return draw; the source exhaust is not consumed twice.
#[derive(Debug, Clone, PartialEq)]
pub struct RecirculatedStream {
    /// Receiving supply node.
    pub supply_node: usize,
    /// Donating exhaust node.
    pub return_node: usize,
    /// Declared receiving-supply fraction.
    pub fraction: f64,
    /// Capacity flow removed from this exhaust, W/K.
    pub capacity_w_per_k: f64,
    /// Temperature at the actual solved return node, K.
    pub temperature_k: f64,
}

/// Independently accumulated outer heat balance after removing return streams.
#[derive(Debug, Clone, PartialEq)]
pub struct RecirculationReport {
    /// Affected supplies, sorted by node.
    pub supplies: Vec<RecirculatedSupply>,
    /// Positive return draws, sorted by (supply, return).
    pub streams: Vec<RecirculatedStream>,
    /// Undrawn exhaust enthalpy minus fresh makeup enthalpy, W.
    /// Distinct from the pressure-network boundary balance in TransportMarch.
    pub external_heat_gain_w: f64,
    /// Outer gain minus heat actually received from the walls, W.
    pub heat_imbalance_w: f64,
    /// Largest independently recomputed mixture-temperature residual, K.
    pub max_mixing_residual_k: f64,
}

#[derive(Debug)]
pub(super) struct Feedback {
    pub(super) links: Vec<RecirculationLink>,
    pub(super) sources: Vec<usize>,
    pub(super) fractions: Vec<f64>,
    pub(super) matrix: Vec<Vec<f64>>,
    zero_heat_matrix: Vec<Vec<f64>>,
    remaining_exhaust: Vec<f64>,
    pub(super) tolerance_k: f64,
}

impl<'a> TransportNetwork<'a> {
    /// Bind imposed return fractions to this exact hydraulic/thermal network.
    /// Inlet declarations now denote FRESH makeup temperatures. The mixed
    /// supply temperatures are unknowns solved on every march. Fixed-flow
    /// wall, fresh-inlet and hA derivatives include the return feedback.
    ///
    /// The return path is outside the passive graph: it must supply whatever
    /// pressure reset is needed, with zero modeled heat addition. Its pressure
    /// loss, fan work and delay are NOT solved. Fully closed (100% return)
    /// supplies refuse. Resource caps and cancellation apply before solid work.
    pub fn with_recirculation(
        mut self, cx: &Cx<'_>, mut links: Vec<RecirculationLink>,
        temperature_tolerance_k: f64,
    ) -> Result<Self, TransportError> {
        poll(cx)?;
        positive(temperature_tolerance_k, "recirculation temperature tolerance")?;
        if self.feedback.is_some() || links.len() > MAX_RECIRCULATION_LINKS {
            return Err(TransportError::InvalidInput("recirculation already bound or link cap exceeded"));
        }
        links.sort_by_key(|link| (link.supply_node, link.return_node));
        let mut previous = None;
        for link in &links {
            poll(cx)?;
            if link.supply_node >= self.order.len() || link.return_node >= self.order.len()
                || self.external_capacity[link.supply_node] <= 0.0
                || self.external_capacity[link.return_node] >= 0.0
                || !(0.0..1.0).contains(&link.fraction)
                || previous == Some((link.supply_node, link.return_node))
            {
                return Err(TransportError::InvalidInput("return links need unique sink-to-supply endpoints and fractions in [0,1)"));
            }
            previous = Some((link.supply_node, link.return_node));
        }
        links.retain(|link| link.fraction > 0.0);
        if links.is_empty() { return Ok(self); }
        let sources: Vec<_> = links.iter().map(|l| l.supply_node).collect::<BTreeSet<_>>().into_iter().collect();
        if sources.len() > MAX_RECIRCULATION_SUPPLIES {
            return Err(TransportError::InvalidInput("recirculation supply cap exceeded"));
        }
        let mut fractions = vec![0.0; sources.len()];
        let mut drawn = vec![0.0; self.order.len()];
        for link in &links {
            poll(cx)?;
            let i = sources.binary_search(&link.supply_node).expect("collected supply");
            fractions[i] = finite(fractions[i] + link.fraction, "total return fraction")?;
            let capacity = positive(link.fraction * self.external_capacity[link.supply_node], "return capacity")?;
            drawn[link.return_node] = finite(drawn[link.return_node] + capacity, "exhaust draw")?;
        }
        if fractions.iter().any(|&f| f >= 1.0) {
            return Err(TransportError::InvalidInput("every recycled supply needs positive fresh makeup"));
        }
        let mut remaining_exhaust = vec![0.0; self.order.len()];
        for node in 0..self.order.len() {
            poll(cx)?;
            let available = (-self.external_capacity[node]).max(0.0);
            if drawn[node] > available {
                return Err(TransportError::Node { node, reason: "return draws exceed available exhaust capacity" });
            }
            remaining_exhaust[node] = available - drawn[node];
        }
        let matrix = feedback_matrix(&self, cx, &sources, &links, true)?;
        let zero_heat_matrix = feedback_matrix(&self, cx, &sources, &links, false)?;
        self.feedback = Some(Feedback { links, sources, fractions, matrix, zero_heat_matrix,
            remaining_exhaust, tolerance_k: temperature_tolerance_k });
        Ok(self)
    }

    /// The imposed return topology; empty for the original once-through model.
    #[must_use]
    pub fn recirculation_links(&self) -> &[RecirculationLink] {
        self.feedback.as_ref().map_or(&[], |f| f.links.as_slice())
    }
}

impl Feedback {
    pub(super) fn march(&self, network: &TransportNetwork<'_>, cx: &Cx<'_>,
        walls: Option<&[f64]>) -> Result<TransportMarch, TransportError> {
        poll(cx)?;
        let reference = network.reference_k;
        let mut inlets = network.inlet_temperatures.clone();
        for &node in &self.sources { inlets[node] = Some(reference); }
        let baseline = network.march_once(cx, walls, &inlets)?;
        // Solve in deviations from a common K reference; do not form products
        // of huge capacity rates and absolute temperatures in the linear solve.
        let mut rhs = Vec::with_capacity(self.sources.len());
        for (i, &node) in self.sources.iter().enumerate() {
            rhs.push((1.0 - self.fractions[i]) *
                (network.inlet_temperatures[node].expect("admitted fresh source") - reference));
        }
        for link in &self.links {
            poll(cx)?;
            let i = self.sources.binary_search(&link.supply_node).expect("admitted source");
            let t = baseline.node_temperatures_k[link.return_node]
                .ok_or(TransportError::InvalidInput("return node has no transported temperature"))?;
            rhs[i] = finite(rhs[i] + link.fraction * (t - reference), "return temperature right-hand side")?;
        }
        let matrix = if walls.is_some() { &self.matrix } else { &self.zero_heat_matrix };
        let solved = solve_feedback(cx, matrix, &rhs, false)?;
        for (&node, &delta) in self.sources.iter().zip(&solved) {
            inlets[node] = Some(positive(reference + delta, "solved mixed inlet")?);
        }
        let mut result = network.march_once(cx, walls, &inlets)?;
        result.recirculation = Some(self.audit(network, cx, &inlets, &result)?);
        poll(cx)?;
        Ok(result)
    }

    fn audit(&self, network: &TransportNetwork<'_>, cx: &Cx<'_>, inlets: &[Option<f64>],
        march: &TransportMarch) -> Result<RecirculationReport, TransportError> {
        let reference = march.enthalpy_reference_k;
        let mut report = RecirculationReport { supplies: Vec::new(), streams: Vec::new(),
            external_heat_gain_w: 0.0, heat_imbalance_w: 0.0, max_mixing_residual_k: 0.0 };
        let mut proposed = Vec::with_capacity(self.sources.len());
        for (i, &node) in self.sources.iter().enumerate() {
            proposed.push((1.0 - self.fractions[i]) *
                (network.inlet_temperatures[node].expect("fresh source") - reference));
        }
        for link in &self.links {
            poll(cx)?;
            let temperature = march.node_temperatures_k[link.return_node]
                .ok_or(TransportError::InvalidInput("accepted return node has no temperature"))?;
            let i = self.sources.binary_search(&link.supply_node).expect("admitted source");
            proposed[i] = finite(proposed[i] + link.fraction * (temperature - reference), "recomputed mixture")?;
            report.streams.push(RecirculatedStream { supply_node: link.supply_node,
                return_node: link.return_node, fraction: link.fraction,
                capacity_w_per_k: link.fraction * network.external_capacity[link.supply_node],
                temperature_k: temperature });
        }
        let mut mixer_defect = 0.0;
        for (i, &node) in self.sources.iter().enumerate() {
            poll(cx)?;
            let mixed = inlets[node].expect("solved supply");
            let residual = finite((mixed - reference) - proposed[i], "mixture residual")?.abs();
            report.max_mixing_residual_k = report.max_mixing_residual_k.max(residual);
            mixer_defect = finite(mixer_defect + residual * network.external_capacity[node], "mixing heat defect")?;
            report.supplies.push(RecirculatedSupply { node,
                fresh_temperature_k: network.inlet_temperatures[node].expect("fresh source"),
                mixed_temperature_k: mixed,
                fresh_capacity_w_per_k: (1.0 - self.fractions[i]) * network.external_capacity[node],
                return_capacity_w_per_k: self.fractions[i] * network.external_capacity[node] });
        }
        if report.max_mixing_residual_k > self.tolerance_k
            || mixer_defect > network.config.absolute_heat_tolerance_w {
            return Err(TransportError::InvalidInput("recirculation mixture temperature or watt residual exceeded budget"));
        }
        let mut gross = 0.0;
        for node in 0..network.order.len() {
            poll(cx)?;
            let heat = if network.external_capacity[node] > 0.0 {
                let fraction = self.sources.binary_search(&node).map_or(0.0, |i| self.fractions[i]);
                -(network.external_capacity[node] * (1.0 - fraction)) *
                    (network.inlet_temperatures[node].expect("fresh source") - reference)
            } else if self.remaining_exhaust[node] > 0.0 {
                self.remaining_exhaust[node] * (march.node_temperatures_k[node]
                    .ok_or(TransportError::InvalidInput("exhaust has no temperature"))? - reference)
            } else { 0.0 };
            report.external_heat_gain_w = finite(report.external_heat_gain_w + heat, "fresh/exhaust heat gain")?;
            gross = finite(gross + heat.abs(), "gross fresh/exhaust heat")?;
        }
        // Keep gross wall exchange, not just its possibly cancelling net sum.
        let mut gross_wall = 0.0;
        for branch in &march.branches {
            if let Some(path) = &branch.march {
                for segment in &path.segments {
                    poll(cx)?;
                    gross_wall = finite(gross_wall + segment.heat_rate_w.abs(), "gross wall heat")?;
                }
            }
        }
        report.heat_imbalance_w = finite(report.external_heat_gain_w - march.wall_heat_rate_w, "outer return-loop heat imbalance")?;
        let tolerance = finite(network.config.absolute_heat_tolerance_w
            + network.config.relative_heat_tolerance * gross.max(gross_wall), "outer heat tolerance")?;
        if report.heat_imbalance_w.abs() > tolerance {
            return Err(TransportError::HeatImbalance { residual_w: report.heat_imbalance_w,
                tolerance_w: tolerance, hydraulic_defect_w: march.hydraulic_energy_defect_w });
        }
        Ok(report)
    }
}

/// Eliminate passive DAG mixing/attenuation without differencing temperatures.
fn feedback_matrix(network: &TransportNetwork<'_>, cx: &Cx<'_>, sources: &[usize],
    links: &[RecirculationLink], exchange: bool) -> Result<Vec<Vec<f64>>, TransportError> {
    let n = sources.len();
    let mut matrix = vec![vec![0.0; n]; n];
    for (i, row) in matrix.iter_mut().enumerate() { row[i] = 1.0; }
    let mut attenuation = vec![1.0; network.models.len()];
    if exchange {
        for (branch, value) in attenuation.iter_mut().enumerate() {
            poll(cx)?;
            if let BranchThermalModel::Exchange(rows) = &network.models[branch] {
                let path = network.path(branch, network.reference_k)?;
                let state = path.march(&vec![network.reference_k; rows.len()])?;
                for segment in &state.segments {
                    poll(cx)?;
                    *value *= fs_math::det::exp(-segment.ntu);
                }
            }
        }
    }
    for (column, &source) in sources.iter().enumerate() {
        poll(cx)?;
        let mut capacities: Vec<f64> = network.external_capacity.iter().map(|c| c.max(0.0)).collect();
        let mut values = vec![0.0; network.order.len()];
        values[source] = 1.0;
        for &node in &network.order {
            poll(cx)?;
            for &branch in &network.outgoing[node] {
                poll(cx)?;
                let down = network.downstream[branch];
                let total = positive(capacities[down] + network.capacities[branch], "feedback mixing capacity")?;
                let arriving = values[node] * attenuation[branch];
                let previous = values[down];
                values[down] = previous + (arriving - previous) * (network.capacities[branch] / total);
                capacities[down] = total;
            }
        }
        for link in links {
            poll(cx)?;
            let row = sources.binary_search(&link.supply_node).expect("collected source");
            matrix[row][column] = finite(matrix[row][column] - link.fraction * values[link.return_node], "feedback coefficient")?;
        }
    }
    Ok(matrix)
}

/// Bounded partial-pivot elimination for the reduced (at most 64 supply) system.
/// Both primal and transpose solves recompute a backward residual. This check
/// is numerical evidence only; no conditioning or interval claim is made.
pub(super) fn solve_feedback(cx: &Cx<'_>, matrix: &[Vec<f64>], rhs: &[f64], transpose: bool)
    -> Result<Vec<f64>, TransportError> {
    poll(cx)?;
    let n = rhs.len();
    if n == 0 || n > MAX_RECIRCULATION_SUPPLIES || matrix.len() != n || matrix.iter().any(|row| row.len() != n) {
        return Err(TransportError::InvalidInput("feedback solve shape or cap"));
    }
    let at = |i: usize, j: usize| if transpose { matrix[j][i] } else { matrix[i][j] };
    let mut a: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| at(i,j)).collect()).collect();
    let mut b = rhs.to_vec();
    for value in a.iter().flatten().chain(&b) { finite(*value, "feedback solve input")?; }
    for k in 0..n {
        poll(cx)?;
        let mut pivot = k;
        for i in k+1..n { if a[i][k].abs() > a[pivot][k].abs() { pivot = i; } }
        if a[pivot][k] == 0.0 { return Err(TransportError::InvalidInput("singular return-temperature system")); }
        a.swap(k, pivot); b.swap(k, pivot);
        for i in k+1..n {
            poll(cx)?;
            let factor = finite(a[i][k] / a[k][k], "feedback elimination multiplier")?;
            a[i][k] = 0.0;
            for j in k+1..n { a[i][j] = finite(a[i][j] - factor * a[k][j], "feedback elimination")?; }
            b[i] = finite(b[i] - factor * b[k], "feedback right-hand elimination")?;
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        poll(cx)?;
        let mut value = b[i];
        for j in i+1..n { value = finite(value - a[i][j] * x[j], "feedback back substitution")?; }
        x[i] = finite(value / a[i][i], "feedback solution")?;
    }
    for i in 0..n {
        poll(cx)?;
        let mut residual = -rhs[i];
        let mut scale = rhs[i].abs();
        for (j, &value) in x.iter().enumerate() {
            let term = finite(at(i,j) * value, "feedback residual term")?;
            residual = finite(residual + term, "feedback residual")?;
            scale = finite(scale + term.abs(), "feedback residual scale")?;
        }
        if residual.abs() > 256.0 * f64::EPSILON * scale.max(f64::MIN_POSITIVE) {
            return Err(TransportError::InvalidInput("return-temperature linear residual did not close"));
        }
    }
    Ok(x)
}

#[cfg(test)]
mod tests;
