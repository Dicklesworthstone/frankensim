//! The original fixed-inventory mass-flow graph on the SAME joint heat step.
use super::*;
use fs_thermal::enthalpy::{CellId, DepositionError, advection::MassBoundary};

/// Original flow powers evaluated at the accepted thermal endpoint. These
/// observations are already included in the native energy account exactly once.
/// Reference-subtracted edge powers need not have the sign of the mass flux.
#[derive(Debug)]
pub struct RelaxationEnthalpyFlow {
    /// Inlet/outlet source powers in the native flow's common-datum coordinates.
    pub endpoint_sources_w: BTreeMap<CellId, f64>,
    /// Original pairwise internal transfers; never external heat or Joule loss.
    pub endpoint_transfers: Vec<HeatTransfer>,
    /// Signed net inlet minus outlet enthalpy, INCLUDED in external_thermal_heat_j.
    pub external_enthalpy_j: f64,
    /// Signed internal Q/T at cells, INCLUDED in unclassified entropy. This is
    /// neither entropy carried by mass nor a mixing/second-law certificate.
    pub internal_at_cell_entropy_j_k: f64,
}

impl<F> FrameRelaxationThermalSystem<'_, '_, '_, '_, F> {
    /// Compose prescribed balanced coolant flow with the existing complete heat
    /// source callback. Every original circuit, field, material-memory and
    /// caloric candidate sees the SAME endpoint. Retain the graph/source pins
    /// across calls; changing inlet data or rates is an explicit intervention.
    ///
    /// The callback has the same obligations as `step_heat_network_with` and
    /// must NOT add this graph's powers itself. Internal flow enters the native
    /// pairwise residual; only actual inlet/outlet enthalpy is external supply.
    /// Flow through a bound structural layer refuses: its mechanical material
    /// history is a fixed inventory, not an Eulerian mixing control volume.
    ///
    /// Native advection admits caloric identities, enthalpy storage and balanced
    /// incoming/outgoing mass at EVERY flow cell. This does not solve momentum,
    /// pumps, pressure work, changing mass or evolving species. Upwind spatial
    /// diffusion and first-order time splitting require independent refinement.
    /// Failure publishes no state, including the borrowed flow observer's audit.
    #[allow(clippy::too_many_arguments)]
    pub fn step_heat_network_with_advection<S, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        external: &[f64],
        heat: RelaxationHeatInput<'_, '_, '_, '_>,
        electrical_input: Option<RelaxationElectricalInput<'_>>,
        solve: FrameThermalSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
        flow: &Advected<()>,
        resolve: impl Fn(&DepositionSimulation<'_>)
            -> core::result::Result<RelaxationHeatCoefficients<S>, E>,
    ) -> core::result::Result<RelaxationHeatReport<S>, RelaxationTransportError<E>> {
        self.step_heat_network_with_optional_flow(thermal, external, heat,
            electrical_input, solve, energy_guess_j, gate, Some(flow), resolve)
    }
}

pub(super) fn admit_inventory(flow: &Advected<()>, solid: &BTreeSet<CellId>,
    mut checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>) -> Result<()> {
    checkpoint()?;
    for edge in flow.flows() {
        checkpoint()?;
        if edge.kg_s != 0.0 && (solid.contains(&edge.a) || solid.contains(&edge.b)) {
            return Err(invalid("mass flow cannot cross a retained structural material inventory"));
        }
    }
    for boundary in flow.boundaries() {
        checkpoint()?;
        let cell = match boundary {
            MassBoundary::Inlet { cell, .. } | MassBoundary::Outlet { cell, .. } => cell,
        };
        if solid.contains(cell) {
            return Err(invalid("mass inlet/outlet cannot replace a structural layer's material history"));
        }
    }
    checkpoint()?;
    Ok(())
}

pub(super) fn evaluate(flow: Option<&Advected<()>>, thermal: &DepositionSimulation<'_>,
    dt: f64, sources: &mut BTreeMap<CellId, f64>,
    mut checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>)
    -> Result<(Option<RelaxationEnthalpyFlow>, Vec<HeatTransfer>)> {
    let Some(flow) = flow else { return Ok((None, Vec::new())); };
    let resolved = flow.endpoint_transport(thermal, &mut checkpoint)?;
    let source = resolved.heat.bulk.sources_w;
    let mut external = 0.0;
    for (&cell, &power) in &source {
        checkpoint()?;
        let energy = finite(dt * power)?;
        if power != 0.0 && energy == 0.0 {
            return Err(invalid("advected boundary energy underflow"));
        }
        external = finite(external + energy)?;
        let old = sources.get(&cell).copied().unwrap_or(0.0);
        let next = finite(old + power)?;
        if power != 0.0 && next == old {
            return Err(invalid("advected boundary power disappears beside independent heat"));
        }
        if power != 0.0 { sources.insert(cell, next); }
    }
    let mut internal_entropy = 0.0;
    for transfer in &resolved.transfers {
        checkpoint()?;
        let a = thermal.temperature_k(transfer.a)?;
        let b = thermal.temperature_k(transfer.b)?;
        // Evaluate p*(1/Tb-1/Ta) without subtracting nearly equal reciprocals.
        internal_entropy = finite(internal_entropy
            + finite(dt * (transfer.power_w / a) * ((a-b) / b))?)?;
    }
    checkpoint()?;
    let report = RelaxationEnthalpyFlow { endpoint_sources_w: source,
        endpoint_transfers: resolved.transfers.clone(), external_enthalpy_j: external,
        internal_at_cell_entropy_j_k: internal_entropy };
    Ok((Some(report), resolved.transfers))
}

/// One coherent endpoint evaluation from an existing physical source owner.
/// The heat callback has NOT already added this graph's transported enthalpy.
/// Keep pressure/material/geometry receipts in `heat.bulk.source.1`; the accepted
/// payload is returned without a post-publication source re-evaluation.
pub struct RelaxationFlowCoefficients<S> {
    pub heat: RelaxationHeatCoefficients<S>,
    pub flow: Advected<()>,
}

impl<F> FrameRelaxationThermalSystem<'_, '_, '_, '_, F> {
    /// Resolve balanced coolant flow and heat from the SAME proposed endpoint.
    /// This composes an existing pressure/flow evaluator with the original
    /// coupled thermal residual; it adds no momentum or flow time integrator.
    ///
    /// `resolve` is repeatable and returns both the graph and heat coefficients
    /// from one trial state. It runs during initial admission and every native
    /// thermal evaluation. Every returned graph is checked for fixed-solid
    /// exclusion, native caloric identity, enthalpy storage, fixed mass balance
    /// and bounded transport. The accepted graph's powers and source payload
    /// accompany the report; discarded evaluations cannot advance flow audits.
    ///
    /// The callback must account for pressure work/hydraulic dissipation rather
    /// than hiding it in convective heat, and preserve original source support.
    /// There is no coefficient lag, cached last-good flow, mass-flow regularizer,
    /// or separately accepted hydraulic state. Source refusal and cancellation
    /// preserve ALL accepted histories through the existing transaction.
    #[allow(clippy::too_many_arguments)]
    pub fn step_heat_network_with_resolved_advection<S, E>(
        &mut self, thermal: &mut DepositionSimulation<'_>, external: &[f64],
        heat: RelaxationHeatInput<'_, '_, '_, '_>,
        electrical_input: Option<RelaxationElectricalInput<'_>>,
        solve: FrameThermalSolve, energy_guess_j: &[f64], gate: &CancelGate,
        resolve: impl Fn(&DepositionSimulation<'_>)
            -> core::result::Result<RelaxationFlowCoefficients<S>, E>,
    ) -> core::result::Result<RelaxationHeatReport<S>, RelaxationTransportError<E>> {
        self.step_heat_network_with_flow_source(thermal, external, heat, electrical_input,
            solve, energy_guess_j, gate, None, |state| {
                resolve(state).map(|resolved| Evaluation {
                    coefficients: resolved.heat, flow: Some(resolved.flow),
                })
            })
    }
}
