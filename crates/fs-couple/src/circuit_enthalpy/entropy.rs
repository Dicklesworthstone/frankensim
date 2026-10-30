//! Entropy on the original electrical/thermal candidate, not a second solver.
//!
//! Physical electrical loss is thermalized at its explicit receiving cells.
//! Stored charge/flux, supplied electrical work and numerical damping do not
//! themselves supply heat entropy. Independent source watts remain unclassified.
//! All integrals and endpoint fluxes use the original thermal material laws.

use super::{
    CircuitEnthalpyCandidate, CircuitEnthalpyError, CircuitEnthalpyInput,
    CircuitEnthalpySolve, CircuitThermalCoefficients, DrivenCoefficients,
    coupled_driver_candidate_entropy_option, finite,
};
use fs_thermal::enthalpy::entropy::EntropyAudit;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::cell::RefCell;
use std::collections::BTreeMap;

/// Thermal entropy over one accepted electrical/thermal interval, in J/K.
///
/// The sum of reservoir input, conduction, thermalized electrical loss, other
/// source-at-cell contribution, known reversible transfer and the three
/// numerical terms equals storage.
/// This is an accounting identity, not a timestep-accuracy certificate. A
/// receiving cell must represent the physical thermalization temperature;
/// unresolved contact/junction temperatures are not inferred from a heat split.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CircuitEntropyAudit {
    /// Exact original caloric integral, including latent storage once per cell.
    pub thermal_storage_change_j_k: f64,
    /// Heat received from each actual ambient/support reservoir divided by its
    /// own Kelvin temperature, not the receiving cell's temperature.
    pub reservoir_entropy_input_j_k: f64,
    /// Original finite-temperature links and reservoir contacts, counted once.
    pub conduction_generation_j_k: f64,
    /// Sum Q_loss/T_cell for the ORIGINAL physical loss deposited in each cell.
    /// No numerical driver loss, stored field energy or electrical work is added.
    pub electrical_loss_generation_j_k: f64,
    /// Independent authored/callback energy divided by receiving-cell Kelvin.
    /// This is not an assigned reservoir entropy or irreversible production.
    pub unclassified_source_at_cell_j_k: f64,
    /// Signed sum of Q/T from a specifically resolved reversible source owner.
    /// Zero for the ordinary resistor path. For a Kelvin thermoelectric leg the
    /// complete junction/bulk sum is zero in exact arithmetic; its signed
    /// floating-point residue remains visible, never irreversible production.
    pub reversible_transfer_entropy_j_k: f64,
    /// Any independent authored or callback source was nonzero, even if signed
    /// sources cancel in one cell or in the total account.
    pub has_unclassified_sources: bool,
    /// Exact caloric storage minus sum m*(h_end-h_start)/T_end. The signed
    /// backward-Euler time-quadrature defect is not physical heating.
    pub temporal_defect_j_k: f64,
    /// Sum of the original per-cell ENERGY residual divided by endpoint Kelvin.
    /// A converged nonlinear solve need not have an exactly zero residual.
    pub algebraic_residual_j_k: f64,
    /// Remaining signed arithmetic/reduction discrepancy. In particular, a
    /// small independent source is not lost by subtracting large Joule totals.
    pub accounting_residual_j_k: f64,
}

/// Solve the SAME resistor/enthalpy interval and return its entropy account.
///
/// Newton, Jacobian stencils and electrical evaluation budgets are unchanged.
/// Entropy is evaluated only during the existing final fresh acceptance, with
/// no additional source callback or circuit solve. Every active thermal law
/// must supply its exact dh/T integral, including passive cells. Refusal or
/// cancellation returns neither physical endpoint. The ordinary energy-only
/// entry point does not perform this work and retains its original return type.
///
/// Only physical resistor loss is classified here. Other source watts, including
/// reversible heat authored by another owner, remain explicitly unclassified.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'m, S, E>(
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    input: CircuitEnthalpyInput<'_>,
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<CircuitThermalCoefficients<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<(CircuitEnthalpyCandidate<'m, S>, CircuitEntropyAudit), CircuitEnthalpyError<E>> {
    let cancel = RefCell::new(cancelled);
    let (next, circuit, report, entropy) = coupled_driver_candidate_entropy_option(
        thermal, target_s, input, &[], controls, initial_energy_increments_j,
        |state, _remaining| {
            let coefficients = resolve(state).map_err(CircuitEnthalpyError::Source)?;
            let circuit = input.circuit.step_candidate(coefficients.resistors, input.ports,
                target_s - thermal.time_s(), controls.circuit, || (*cancel.borrow_mut())())?;
            Ok(DrivenCoefficients {
                circuit, sources_w: coefficients.sources_w,
                conductances_w_k: coefficients.conductances_w_k,
                source: coefficients.source, circuit_trials: 1,
            })
        }, true, || (*cancel.borrow_mut())(),
    )?;
    let entropy = entropy.ok_or(CircuitEnthalpyError::Invalid("missing requested circuit entropy"))?;
    Ok((CircuitEnthalpyCandidate { thermal: next, circuit, report }, entropy))
}

/// The same native entropy account with total ambient/reservoir conductances
/// supplied by the physical boundary owner. Actual reservoir temperatures and
/// accepted endpoint heat supply entropy directly; no geometry correction is
/// misclassified as applied work or counted alongside the original boundary.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate_with_boundary_conductances<'m, S, E>(
    thermal: &DepositionSimulation<'m>, target_s: f64,
    input: CircuitEnthalpyInput<'_>, controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<super::BoundaryCircuitThermalCoefficients<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<(CircuitEnthalpyCandidate<'m, S>, CircuitEntropyAudit), CircuitEnthalpyError<E>> {
    let (candidate, entropy) = super::boundary_candidate(thermal, target_s, input, controls,
        initial_energy_increments_j, resolve, true, cancelled)?;
    Ok((candidate, entropy.ok_or(CircuitEnthalpyError::Invalid("missing requested circuit entropy"))?))
}

#[derive(Default)]
pub(super) struct Sum { value: f64, correction: f64 }
impl Sum {
    pub(super) fn add(&mut self, term: f64) -> Result<(), DepositionError> {
        let term = finite(term)?;
        let next = finite(self.value + term)?;
        let lost = if self.value.abs() >= term.abs() {
            (self.value - next) + term
        } else { (term - next) + self.value };
        self.correction = finite(self.correction + lost)?;
        self.value = next;
        Ok(())
    }
    pub(super) fn total(self) -> Result<f64, DepositionError> { finite(self.value + self.correction) }
}
pub(super) fn power_entropy(power: f64, dt: f64, temperature: f64) -> Result<f64, DepositionError> {
    let heat = finite(dt * finite(power)?)?;
    let entropy = finite(heat / temperature)?;
    if power != 0.0 && (heat == 0.0 || entropy == 0.0) {
        return Err(DepositionError::Arithmetic);
    }
    Ok(entropy)
}

// Capture independent sources BEFORE mixing them with physical electrical heat.
// O(cells + sparse sources), only on the final candidate. Temperatures are in
// the exact stable CellId order used by the native energy residual.
pub(super) struct SourceSample {
    temperatures_k: Vec<f64>,
    energy_at_endpoint_temperature_j_k: f64,
    electrical_loss_j_k: f64,
    other_source_j_k: f64,
    has_other_sources: bool,
}
impl SourceSample {
    pub(super) fn new(
        initial: &DepositionSimulation<'_>,
        state: &DepositionSimulation<'_>,
        losses_w: &BTreeMap<CellId, f64>,
        sources_w: &BTreeMap<CellId, f64>,
        dt: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if !dt.is_finite() || dt <= 0.0 || initial.cells().len() != state.cells().len() {
            return Err(DepositionError::Arithmetic);
        }
        for (index, (&id, &power)) in losses_w.iter().enumerate() {
            if index % 256 == 0 { checkpoint()?; }
            if !state.cells().contains_key(&id) || !power.is_finite() || power < 0.0 {
                return Err(DepositionError::InvalidInput("entropy needs nonnegative physical loss at active cells"));
            }
        }
        for (index, (&id, &power)) in sources_w.iter().enumerate() {
            if index % 256 == 0 { checkpoint()?; }
            if !state.cells().contains_key(&id) || !power.is_finite() {
                return Err(DepositionError::InvalidInput("entropy source must name an active cell"));
            }
        }
        let mut temperatures_k = Vec::with_capacity(state.cells().len());
        let mut energy = Sum::default();
        let mut losses = Sum::default();
        let mut other = Sum::default();
        let mut has_other_sources = false;
        for (index, (&id, cell)) in state.cells().iter().enumerate() {
            if index % 256 == 0 { checkpoint()?; }
            let old = initial.cells().get(&id).ok_or(DepositionError::Arithmetic)?;
            let temperature = state.temperature_k(id)?;
            if !temperature.is_finite() || temperature <= 0.0 || old.mass_kg != cell.mass_kg {
                return Err(DepositionError::MaterialDomain);
            }
            temperatures_k.push(temperature);
            let stored = finite(cell.mass_kg * finite(cell.specific_h_j_kg - old.specific_h_j_kg)?)?;
            energy.add(finite(stored / temperature)?)?;
            let loss = losses_w.get(&id).copied().unwrap_or(0.0);
            if loss != 0.0 { losses.add(power_entropy(loss, dt, temperature)?)?; }
            for source in [cell.boundary.source_w, sources_w.get(&id).copied().unwrap_or(0.0)] {
                if source != 0.0 {
                    has_other_sources = true;
                    other.add(power_entropy(source, dt, temperature)?)?;
                }
            }
        }
        checkpoint()?;
        Ok(Self { temperatures_k,
            energy_at_endpoint_temperature_j_k: energy.total()?,
            electrical_loss_j_k: losses.total()?, other_source_j_k: other.total()?,
            has_other_sources })
    }

    pub(super) fn finish(
        &self,
        thermal: &EntropyAudit,
        residuals_j: &[f64],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<CircuitEntropyAudit, DepositionError> {
        checkpoint()?;
        if residuals_j.len() != self.temperatures_k.len() { return Err(DepositionError::Arithmetic); }
        let mut algebraic = Sum::default();
        for (index, (&residual, &temperature)) in residuals_j.iter().zip(&self.temperatures_k).enumerate() {
            if index % 256 == 0 { checkpoint()?; }
            algebraic.add(finite(residual / temperature)?)?;
        }
        let algebraic = algebraic.total()?;
        let temporal = finite(thermal.thermal_storage_change_j_k - self.energy_at_endpoint_temperature_j_k)?;
        checkpoint()?;
        let mut audit = CircuitEntropyAudit {
            thermal_storage_change_j_k: thermal.thermal_storage_change_j_k,
            reservoir_entropy_input_j_k: thermal.reservoir_entropy_input_j_k,
            conduction_generation_j_k: thermal.conduction_generation_j_k,
            electrical_loss_generation_j_k: self.electrical_loss_j_k,
            unclassified_source_at_cell_j_k: self.other_source_j_k,
            has_unclassified_sources: self.has_other_sources,
            reversible_transfer_entropy_j_k: 0.0,
            temporal_defect_j_k: temporal, algebraic_residual_j_k: algebraic,
            accounting_residual_j_k: 0.0,
        };
        audit.classify_sources(self.other_source_j_k, self.has_other_sources, 0.0)?;
        Ok(audit)
    }
}

impl CircuitEntropyAudit {
    // A specific physical owner may identify previously unclassified sources.
    // Rebuild the account from directly accumulated original terms, never by
    // subtracting large known transfer from a small unknown source remainder.
    pub(super) fn classify_sources(
        &mut self, independent_j_k: f64, has_independent: bool, reversible_j_k: f64,
    ) -> Result<(), DepositionError> {
        self.unclassified_source_at_cell_j_k = finite(independent_j_k)?;
        self.has_unclassified_sources = has_independent;
        self.reversible_transfer_entropy_j_k = finite(reversible_j_k)?;
        let mut accounting = Sum::default();
        for value in [self.thermal_storage_change_j_k,
            -self.reservoir_entropy_input_j_k, -self.conduction_generation_j_k,
            -self.electrical_loss_generation_j_k, -self.unclassified_source_at_cell_j_k,
            -self.reversible_transfer_entropy_j_k,
            -self.temporal_defect_j_k, -self.algebraic_residual_j_k] {
            accounting.add(value)?;
        }
        self.accounting_residual_j_k = accounting.total()?;
        Ok(())
    }
}
