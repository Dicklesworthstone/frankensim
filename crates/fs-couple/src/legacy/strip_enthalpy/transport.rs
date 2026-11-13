//! Source-pinned heat transport in the original layered self-heating solve.
//!
//! Prepare the existing MaterialNetwork once per declared source/topology epoch.
//! Its original resolver evaluates conductivity at every proposed endpoint; this
//! adapter supplies those coefficients to the original strip/enthalpy residual.
//! No electrical state is advanced and resistor-bearing networks are refused.

use std::cell::RefCell;
use std::collections::BTreeSet;
use crate::circuit_enthalpy::material::{MaterialNetworkError, MaterialNetworkSample,
    PreparedMaterialNetwork};
use fs_material::thermovisco::source::SourceError;
use fs_solid::deposition::{Control, LayeredStrip};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use super::{LayerThermalBinding, StripEnthalpyCandidate, StripEnthalpyConfig,
    StripEnthalpyError, ThermalCoefficients};

/// Original source refusals remain available to the caller.
#[derive(Debug)]
pub enum TransportError {
    Network(MaterialNetworkError),
    Context { layer: u64, error: SourceError },
    Invalid(&'static str),
}
impl core::fmt::Display for TransportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "layered material transport: {self:?}")
    }
}
impl std::error::Error for TransportError {}

/// Resolve the SAME pinned transport claims at the actual thermal endpoint.
///
/// All thermal cells/links, including passive bodies, remain subject to the
/// material network's complete-binding checks. Conducting mechanical layers
/// additionally require the mechanical card, fixed coordinates and temperature
/// convention. Missing conductivity is not interpreted as insulation. An
/// isolated layer may omit it, as admitted by the original network.
///
/// The returned source retains complete original material samples. Independent
/// heat sources can be added to `sources_w` before returning it from a custom
/// callback; never add Maxwell dissipation, which the strip owner deposits once.
/// Resistor branches refuse rather than silently discarding electrical work.
pub fn coefficients(
    transport: &PreparedMaterialNetwork<'_, '_>, thermal: &DepositionSimulation<'_>,
    bindings: &[LayerThermalBinding<'_>], mut cancelled: impl FnMut() -> bool,
) -> Result<ThermalCoefficients<MaterialNetworkSample>, TransportError> {
    let cancelled_error = || TransportError::Network(
        MaterialNetworkError::Thermal(DepositionError::Cancelled));
    if cancelled() { return Err(cancelled_error()); }
    if bindings.len() > thermal.cells().len() {
        return Err(TransportError::Invalid("mechanical bindings exceed active thermal cells"));
    }
    let resolved = transport.resolve(thermal, &mut cancelled).map_err(TransportError::Network)?;
    if !resolved.resistors.is_empty() {
        return Err(TransportError::Invalid("electrical branches require the circuit energy owner"));
    }
    let mut layers = BTreeSet::new();
    let mut cells = BTreeSet::new();
    for binding in bindings {
        if cancelled() { return Err(cancelled_error()); }
        if !layers.insert(binding.layer) || !cells.insert(binding.cell) {
            return Err(TransportError::Invalid("duplicate mechanical layer or heat destination"));
        }
        // MaterialNetworkSample has stable sorted CellId order.
        let index = resolved.source.cells.binary_search_by_key(&binding.cell, |s| s.cell)
            .map_err(|_| TransportError::Invalid("mechanical layer lacks a transport cell"))?;
        let sample = &resolved.source.cells[index];
        if let Some(state) = &sample.conductivity {
            binding.material.admit_state_point(state).map_err(|error|
                TransportError::Context { layer: binding.layer, error })?;
        }
    }
    if cancelled() { return Err(cancelled_error()); }
    Ok(ThermalCoefficients { sources_w: resolved.sources_w,
        conductances_w_k: resolved.conductances_w_k, source: resolved.source })
}

/// Advance compatible deformation, physical self-heating and material conduction.
///
/// Reuse `transport` across accepted steps to retain source selection. Each
/// Newton/Jacobian trial re-evaluates k(T) rather than preparing a new selection.
/// Final reports contain the exact material samples used by the accepted
/// residual. Unsupported source domains reject both histories without fallback.
/// Geometry and caloric/Maxwell limits remain those of the original owners.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'t, 'm>(
    thermal: &DepositionSimulation<'t>, strip: &LayeredStrip<'m>,
    bindings: &[LayerThermalBinding<'_>], transport: &PreparedMaterialNetwork<'_, '_>,
    target_s: f64, control: Control, config: StripEnthalpyConfig,
    initial_energy_increments_j: &[f64], cancelled: impl FnMut() -> bool,
) -> Result<StripEnthalpyCandidate<'t, 'm, MaterialNetworkSample>, StripEnthalpyError<TransportError>> {
    let cancel = RefCell::new(cancelled);
    super::step_candidate(thermal, strip, bindings, target_s, control, config,
        initial_energy_increments_j,
        |state| coefficients(transport, state, bindings, || (*cancel.borrow_mut())()),
        || (*cancel.borrow_mut())())
}

impl<'t, 'm> super::growth::GrowingStrip<'t, 'm> {
    /// Material transport over one interval, followed by the original atomic
    /// endpoint births. The prepared network describes the interval's OLD
    /// topology; newborn cells are not aged or queried inside Newton trials.
    ///
    /// After a topology change the caller must explicitly prepare a complete
    /// network for the next interval. This method never silently reselects
    /// sources or infers newborn contacts. The report's source describes the
    /// interval's left limit, while events account for right-limit births.
    #[allow(clippy::too_many_arguments)]
    pub fn step_transport_candidate(
        &self, transport: &PreparedMaterialNetwork<'_, '_>, target_s: f64,
        control: Control, config: StripEnthalpyConfig, initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, super::growth::GrowthReport<MaterialNetworkSample>),
        super::growth::GrowthError<TransportError>> {
        let cancel = RefCell::new(cancelled);
        self.step_candidate(target_s, control, config, initial_energy_increments_j,
            |state| coefficients(transport, state, self.bindings(), || (*cancel.borrow_mut())()),
            || (*cancel.borrow_mut())())
    }
}
