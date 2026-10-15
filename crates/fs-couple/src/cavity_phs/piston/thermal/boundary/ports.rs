//! Multiple thermal contacts on ONE gas, evaluated at ONE energy midpoint.
//! Distinct hot and cold contacts remain distinct even when their net heat is
//! zero. No averaged temperature, repeated gas step, or inferred mixture law.
use super::*;
use std::cell::RefCell;
use std::collections::BTreeMap;
pub use super::heat::{BoundaryHeatCandidate, BoundaryHeatStep};

/// A stable physical contact and its other-side temperature for this trial.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalPort {
    /// Caller-owned contact identity, unique within this gas boundary.
    pub id: u64,
    /// Positive absolute temperature [K], not a second stored thermal state.
    pub temperature_k: f64,
}
/// Actual heat and original material evaluation for one contact.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalPortSample<S> {
    /// Exact input contact identity and trial temperature.
    pub port: ThermalPort,
    /// Actual conductance [W/K] at the shared gas midpoint.
    pub conductance_w_k: f64,
    /// dt*G*(T_contact-T_gas_mid) [J], positive INTO gas. A finite counterpart
    /// debits THIS signed value once. Never split net heat by conductance: hot
    /// and cold contacts can carry nonzero opposite heats with zero net input.
    pub heat_j: f64,
    /// Original source result at this contact; no type-erased material preset.
    pub source: S,
}
/// Source records are in stable contact-ID order, independent of input order.
pub type ThermalPortsCandidate<S> = BoundaryHeatCandidate<Vec<ThermalPortSample<S>>>;

impl GasPiston {
    /// Simultaneous heat transfer through multiple boundaries on the same gas.
    /// One scalar solve uses the SUM of contact heats; every material callback
    /// sees the same actual gas midpoint and every contact retains its own heat.
    /// Empty contacts reproduce the isolated gas/momentum step. Explicit zero
    /// conductance does not bypass source evaluation or state admission.
    ///
    /// `maximum_ports` bounds allocation and per-trial callbacks before work.
    /// IDs must be unique. Callbacks must be repeatable, without changes to live
    /// physical history. Cancellation is polled between contact evaluations.
    /// The existing feedback and piston controls bound all nested solves.
    ///
    /// A finite network must debit EACH evaluated contact heat and accept all
    /// endpoints together. The difference between their sum and the heat used
    /// by the gas is explicitly bounded by the source and total-energy gates;
    /// contact heat is never rescaled to conceal that residual. Temperatures
    /// supplied by a network must come from its same nonlinear trial. This is a
    /// uniform gas with thermal ports, not gas mixing, mass flow or convection.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_ports_candidate<S, E>(
        &self, input: PistonInput, controls: PistonSolve, exchange: HeatExchangeSolve,
        ports: &[ThermalPort], maximum_ports: usize,
        mut resolve: impl FnMut(ThermalPort, &ResolvedThermochemicalGasState)
            -> Result<ConductanceSample<S>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ThermalPortsCandidate<S>, ThermalPistonError<E>> {
        let cancel = RefCell::new(cancelled);
        poll(&mut || (*cancel.borrow_mut())())?;
        if ports.len() > maximum_ports {
            return Err(ThermalPistonError::Invalid("thermal contact budget exceeded"));
        }
        let mut ordered = BTreeMap::new();
        for port in ports {
            poll(&mut || (*cancel.borrow_mut())())?;
            if !positive(port.temperature_k) || ordered.insert(port.id, *port).is_some() {
                return Err(ThermalPistonError::Invalid("unique contact IDs and positive finite temperatures required"));
            }
        }
        heat::solve(self, input, controls, exchange, |gas| {
            let mut total = 0.0;
            let mut sources = Vec::with_capacity(ordered.len());
            for &port in ordered.values() {
                poll(&mut || (*cancel.borrow_mut())())?;
                let sample = resolve(port, gas).map_err(ThermalPistonError::Interface)?;
                let q = heat::conductive_heat(input.duration_s, sample.conductance_w_k,
                    port.temperature_k, gas.state().temperature)?;
                total += q;
                if !total.is_finite() {
                    return Err(ThermalPistonError::Invalid("summed contact heat is unrepresentable"));
                }
                sources.push(ThermalPortSample { port, conductance_w_k: sample.conductance_w_k,
                    heat_j: q, source: sample.source });
            }
            Ok((total, sources))
        }, || (*cancel.borrow_mut())())
    }
}

#[cfg(test)]
mod tests;
