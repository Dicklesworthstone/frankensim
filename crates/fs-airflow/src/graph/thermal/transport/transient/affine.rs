//! Conservative direction-resolved transport over an affine flow interval.
//!
//! The original momentum owner supplies initial, mean and final signed rates.
//! Each reversing connection becomes two native donor-cell transfers, carrying
//! the independently resolved enthalpy of each donor. This exactly integrates
//! the positive/negative swept volume of the affine interpolant, while retaining
//! the native held-endpoint thermal approximation. It is not exact time-ordered
//! transport; time/spatial refinement remains necessary for residence and mixing.
//!
//! External streams are explicit. Their separate inlet/outlet rates must equal
//! the directional integrals of the node's net affine flow, rather than the
//! positive/negative parts of its signed mean. Each inlet keeps its authored
//! source enthalpy; no hydraulic pressure boundary invents a reservoir state.

use super::*;
use fs_math::dd::Dd;

pub mod rates;
pub use rates::AffineVolumeRates;

/// Actual directional external mass rates in original hydraulic node order.
/// Both can be positive during a reversal. These are rate requirements only;
/// the caller must still author each entering stream's original caloric state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AffineNodeMassStreams {
    pub inlet_kg_s: f64,
    pub outlet_kg_s: f64,
}

/// Integrate each node's NET external flow before separating its directions.
/// Internal branch exchanges are summed with compensation in their original
/// graph order. Splitting each incident flow first would invent reservoir
/// mixing at an ordinary through-flow junction.
pub fn node_mass_streams(
    cx: &Cx<'_>,
    graph: &LossGraph,
    reference_density: Density,
    initial_m3_s: &[f64],
    mean_m3_s: &[f64],
    final_m3_s: &[f64],
) -> Result<Vec<AffineNodeMassStreams>, TransientTransportError> {
    let rates = admit_rates(cx, graph, initial_m3_s, mean_m3_s, final_m3_s)?;
    let rho = reference_density.value();
    if !rho.is_finite() || rho <= 0.0 {
        return Err(TransientTransportError::Invalid(
            "positive affine reference density required",
        ));
    }
    let zero = Dd::from_f64(0.0);
    let mut nodes = vec![[zero; 3]; graph.node_count()];
    let mut scales = vec![zero; graph.node_count()];
    for (branch, edge) in graph.branches().iter().enumerate() {
        poll(cx)?;
        let values = [
            initial_m3_s[branch],
            rates[branch].mean_m3_s,
            final_m3_s[branch],
        ];
        for (node, sign) in [(edge.from, 1.0), (edge.to, -1.0)] {
            for (sum, rate) in nodes[node].iter_mut().zip(values) {
                *sum = *sum + Dd::from_f64(sign * rate);
            }
            scales[node] =
                scales[node] + Dd::from_f64(values.into_iter().map(f64::abs).fold(0.0, f64::max));
        }
    }
    let mut result = Vec::with_capacity(nodes.len());
    for (node, sums) in nodes.into_iter().enumerate() {
        poll(cx)?;
        let [initial, mean, final_flow] = sums.map(Dd::to_f64);
        let scale = scales[node].to_f64() * rho;
        // Summation and the native port mean have only arithmetic accuracy.
        // Integrate the actual endpoint incidence; never floor a small physical
        // supply merely because a much larger through-flow crosses this node.
        let allowance = 128.0 * f64::EPSILON * scales[node].to_f64();
        let endpoint_mean = 0.5 * initial + 0.5 * final_flow;
        if !allowance.is_finite() || (mean - endpoint_mean).abs() > allowance {
            return Err(TransientTransportError::Invalid(
                "node incidence lost the affine mean",
            ));
        }
        let rates = AffineVolumeRates::try_new(initial, endpoint_mean, final_flow)
            .map_err(TransientTransportError::Invalid)?;
        let inlet = rho * rates.forward_m3_s;
        let outlet = rho * rates.reverse_m3_s;
        if ![inlet, outlet, scale].into_iter().all(f64::is_finite)
            || (rates.forward_m3_s > 0.0 && inlet == 0.0)
            || (rates.reverse_m3_s > 0.0 && outlet == 0.0)
        {
            return Err(TransientTransportError::Invalid(
                "affine external mass rate is unrepresentable",
            ));
        }
        result.push(AffineNodeMassStreams {
            inlet_kg_s: inlet,
            outlet_kg_s: outlet,
        });
    }
    Ok(result)
}

pub(super) fn admit_rates(
    cx: &Cx<'_>,
    graph: &LossGraph,
    initial: &[f64],
    mean: &[f64],
    final_flow: &[f64],
) -> Result<Vec<AffineVolumeRates>, TransientTransportError> {
    poll(cx)?;
    if [initial.len(), mean.len(), final_flow.len()]
        .into_iter()
        .any(|n| n != graph.branches().len())
    {
        return Err(TransientTransportError::Invalid(
            "one original affine history per branch required",
        ));
    }
    initial
        .iter()
        .zip(mean)
        .zip(final_flow)
        .map(|((&a, &m), &b)| {
            poll(cx)?;
            AffineVolumeRates::try_new(a, m, b).map_err(TransientTransportError::Invalid)
        })
        .collect()
}

pub(super) fn admit_boundaries(
    cx: &Cx<'_>,
    node_cells: &[CellId],
    expected: &[AffineNodeMassStreams],
    boundaries: &[MassBoundary],
) -> Result<(), TransientTransportError> {
    if node_cells.len() != expected.len() {
        return Err(TransientTransportError::Invalid(
            "complete affine hydraulic node map required",
        ));
    }
    let mut actual: BTreeMap<CellId, [Dd; 2]> = node_cells
        .iter()
        .map(|&id| (id, [Dd::from_f64(0.0); 2]))
        .collect();
    if actual.len() != node_cells.len() {
        return Err(TransientTransportError::Invalid(
            "distinct affine hydraulic node cells required",
        ));
    }
    for boundary in boundaries {
        poll(cx)?;
        let (cell, kg_s, side) = match *boundary {
            MassBoundary::Inlet { cell, kg_s, .. } => (cell, kg_s, 0),
            MassBoundary::Outlet { cell, kg_s, .. } => (cell, kg_s, 1),
        };
        let sums = actual
            .get_mut(&cell)
            .ok_or(TransientTransportError::Invalid(
                "affine external stream must bind an original hydraulic junction",
            ))?;
        if !kg_s.is_finite() || kg_s <= 0.0 {
            return Err(TransientTransportError::Invalid(
                "positive finite affine boundary mass rate required",
            ));
        }
        sums[side] = sums[side] + Dd::from_f64(kg_s);
    }
    for (&cell, expected) in node_cells.iter().zip(expected) {
        poll(cx)?;
        for (actual, required) in actual[&cell]
            .iter()
            .zip([expected.inlet_kg_s, expected.outlet_kg_s])
        {
            let actual = actual.to_f64();
            // Only the actual directional supply sets this physical-transfer
            // scale. A much larger unrelated through-flow cannot hide a
            // missing entering enthalpy or authorize a fictitious reservoir.
            let allowance = 128.0 * f64::EPSILON * required.max(actual);
            if !actual.is_finite()
                || (actual == 0.0) != (required == 0.0)
                || (actual - required).abs() > allowance
            {
                return Err(TransientTransportError::Invalid(
                    "external inlet and outlet must each match the affine swept mass, including reversal",
                ));
            }
        }
    }
    Ok(())
}

impl GraphEnthalpyTransport {
    /// Bind one affine branch history using both original upwind donors.
    ///
    /// The mean must be the original owner's mean, consistent with its initial
    /// and final rates. Base connection IDs stay in original branch order;
    /// a reversing branch additionally uses `branch_count + branch_index`.
    /// External IDs must be disjoint. Every actual transfer consumes the
    /// existing connection budget and native outgoing-mass/CFL accounting.
    #[allow(clippy::too_many_arguments)]
    pub fn from_affine_volume_flows(
        cx: &Cx<'_>,
        graph: &LossGraph,
        initial_m3_s: &[f64],
        mean_m3_s: &[f64],
        final_m3_s: &[f64],
        reference_density: Density,
        node_cells: &[CellId],
        thermal: &DepositionSimulation<'_>,
        boundaries: Vec<MassBoundary>,
        maximum_connections: usize,
    ) -> Result<Self, TransientTransportError> {
        let rates = admit_rates(cx, graph, initial_m3_s, mean_m3_s, final_m3_s)?;
        let streams = node_mass_streams(
            cx,
            graph,
            reference_density,
            initial_m3_s,
            mean_m3_s,
            final_m3_s,
        )?;
        admit_boundaries(cx, node_cells, &streams, &boundaries)?;
        let base = Self::from_volume_flows(
            cx,
            graph,
            mean_m3_s,
            reference_density,
            node_cells,
            thermal,
            boundaries,
            maximum_connections,
        )?;
        base.with_affine_rates(cx, thermal, reference_density, &rates)
    }

    pub(super) fn with_affine_rates(
        mut self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        density: Density,
        rates: &[AffineVolumeRates],
    ) -> Result<Self, TransientTransportError> {
        let count = self.advection.flows().len();
        if rates.len() != count {
            return Err(TransientTransportError::Invalid(
                "complete expanded affine transport history required",
            ));
        }
        let additional = rates.iter().filter(|r| r.counterflow_m3_s() > 0.0).count();
        if count
            .checked_add(additional)
            .and_then(|n| n.checked_add(self.advection.boundaries().len()))
            .is_none_or(|n| n > self.max_connections)
        {
            return Err(TransientTransportError::Invalid(
                "affine swept transport exceeds the native connection budget",
            ));
        }
        if additional == 0 {
            return Ok(self);
        }
        let mut flows = self.advection.flows().to_vec();
        flows.reserve(additional);
        for (index, rate) in rates.iter().enumerate() {
            poll(cx)?;
            if rate.counterflow_m3_s() == 0.0 {
                continue;
            }
            let original = &mut flows[index];
            let (primary, opposite) = if rate.mean_m3_s >= 0.0 {
                (rate.forward_m3_s, -rate.reverse_m3_s)
            } else {
                (-rate.reverse_m3_s, rate.forward_m3_s)
            };
            let primary = density.value() * primary;
            let opposite = density.value() * opposite;
            if !primary.is_finite() || !opposite.is_finite() || primary == 0.0 || opposite == 0.0 {
                return Err(TransientTransportError::Invalid(
                    "affine internal swept mass is unrepresentable",
                ));
            }
            original.kg_s = primary;
            let additional = MassFlow {
                id: u64::try_from(count.checked_add(index).ok_or(
                    TransientTransportError::Invalid("affine transfer ID overflow"),
                )?)
                .map_err(|_| TransientTransportError::Invalid("affine transfer ID overflow"))?,
                a: original.a,
                b: original.b,
                kg_s: opposite,
            };
            flows.push(additional);
        }
        self.advection = Advected::new((), flows, self.max_connections)?
            .with_boundaries(self.advection.boundaries().to_vec(), self.max_connections)?;
        self.endpoint_transport(cx, thermal)?;
        Ok(self)
    }
}
