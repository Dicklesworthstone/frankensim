//! Reversible Kelvin heat is transported entropy, not electrical dissipation.
//!
//! Reuse the exact accepted Seebeck/Peltier/Thomson samples and original thermal
//! states. Independent authored sources are read directly before classification;
//! subtracting a large reversible heat aggregate is never used to recover them.

use super::*;
use crate::circuit_enthalpy::entropy::{Sum, power_entropy};

/// Advance the original coupled thermoelectric owner with final-only entropy.
///
/// Joule loss and Fourier conduction produce entropy. Signed Kelvin exchange
/// transports entropy between the actual junction and bulk temperatures, with
/// zero sum for EACH leg in exact arithmetic. Its small signed numerical residue
/// is reported as reversible transfer, not clamped or relabeled as production.
/// Other authored sources stay unclassified. Cell storage is integrated once,
/// even when several legs share a junction or a bulk cell.
///
/// This uses the same source pins, Newton solve, evaluation budgets and energy
/// gates as `super::step_candidate`. Entropy runs on its existing final fresh
/// trial; it adds no transport/Seebeck source resolution, circuit solve, or time
/// advancement. The Kelvin
/// check is a floating-point identity check, not a continuum accuracy bound.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'m>(
    material: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
    circuit: &MaterialCircuit, ports: &[f64], legs: &[ThermoelectricLeg],
    electrical_bases: &[(usize, ElectricalPropertyBasis)], controls: CircuitEnthalpySolve,
    guess_j: &[f64], mut cancelled: impl FnMut() -> bool,
) -> Result<(ThermoelectricCandidate<'m>, CircuitEntropyAudit), Error> {
    let dt = target_s - thermal.time_s();
    let (thermal, next_circuit, coupled, audit) = coupled_candidate_entropy_option(
        material, thermal, target_s, circuit, ports, legs, electrical_bases,
        controls, guess_j,
        |sources, inputs, dt, cancelled| circuit.step_candidate(sources, inputs, dt,
            controls.circuit, cancelled), true, &mut cancelled,
    )?;
    let mut audit = audit.ok_or(invalid("missing requested thermoelectric entropy"))?;
    classify(&thermal, &coupled.source.legs, dt, &mut audit, &mut cancelled)?;
    let candidate = finish_candidate(thermal, next_circuit, coupled, controls, cancelled)?;
    Ok((candidate, audit))
}

// A Kelvin leg has two entropy edges and three thermal stores. Check each leg
// independently so errors in distinct legs cannot cancel. The bound covers a
// small fixed number of rounded products, quotients and compensated additions;
// it is explicitly not a thermodynamic modeling or time-discretization budget.
fn kelvin_residual(powers_w: [f64; 3], temperatures_k: [f64; 3], dt: f64)
    -> Result<f64, Error> {
    let mut sum = Sum::default();
    let mut magnitude = Sum::default();
    for (power, temperature) in powers_w.into_iter().zip(temperatures_k) {
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(invalid("Kelvin entropy requires positive finite original temperatures"));
        }
        let term = power_entropy(power, dt, temperature)?;
        sum.add(term)?;
        magnitude.add(term.abs())?;
    }
    let residual = sum.total()?;
    let allowance = finite(128.0 * f64::EPSILON * magnitude.total()?)?;
    if residual.abs() > allowance {
        return Err(invalid("thermoelectric junction and bulk entropy transfers do not close"));
    }
    Ok(residual)
}

fn classify(
    thermal: &DepositionSimulation<'_>, legs: &[ThermoelectricLegSample], dt: f64,
    audit: &mut CircuitEntropyAudit, mut cancelled: impl FnMut() -> bool,
) -> Result<(), Error> {
    let mut other = Sum::default();
    let mut reversible = Sum::default();
    let mut has_other = false;
    for (index, (&id, cell)) in thermal.cells().iter().enumerate() {
        if index % 256 == 0 && cancelled() { return Err(CircuitEnthalpyError::Cancelled); }
        if cell.boundary.source_w != 0.0 {
            has_other = true;
            other.add(power_entropy(cell.boundary.source_w, dt, thermal.temperature_k(id)?)?)?;
        }
    }
    for leg in legs {
        if cancelled() { return Err(CircuitEnthalpyError::Cancelled); }
        let [a, b] = leg.material.temperature_k();
        let bulk = leg.material.bulk_temperature_k();
        if a != thermal.temperature_k(leg.terminals[0])?
            || b != thermal.temperature_k(leg.terminals[1])?
            || bulk != thermal.temperature_k(leg.bulk)? {
            return Err(invalid("Kelvin entropy must use the same accepted material temperatures"));
        }
        reversible.add(kelvin_residual(
            [leg.transfer.terminal_heat_w[0], leg.transfer.bulk_heat_w, leg.transfer.terminal_heat_w[1]],
            [a, bulk, b], dt,
        )?)?;
    }
    audit.classify_sources(other.total()?, has_other, reversible.total()?)?;
    if cancelled() { return Err(CircuitEnthalpyError::Cancelled); }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn original_kelvin_edges_reverse_without_producing_entropy() {
        for (left, right) in [(0.02, -0.01), (-0.03, 0.04), (0.0, 0.0), (0.02, 0.02)] {
            for current in [-2.0, 0.0, 2.0] {
                let temperatures = [400.0, 430.0, 300.0];
                let heat = [-left*temperatures[0]*current,
                    (left-right)*temperatures[1]*current, right*temperatures[2]*current];
                let value = kelvin_residual(heat, temperatures, 0.25).unwrap();
                assert!(value.abs() < 1e-16);
            }
        }
    }

    #[test]
    fn omitting_variable_seebeck_bulk_heat_or_using_the_wrong_temperature_refuses() {
        let temperatures = [400.0, 430.0, 300.0];
        let heat = [-8.0, 12.9, -3.0];
        assert!(kelvin_residual(heat, temperatures, 0.1).is_ok());
        assert!(kelvin_residual([heat[0], 0.0, heat[2]], temperatures, 0.1).is_err());
        assert!(kelvin_residual(heat, [400.0, 300.0, 300.0], 0.1).is_err());
        assert!(kelvin_residual(heat, [400.0, 0.0, 300.0], 0.1).is_err());
    }
}
