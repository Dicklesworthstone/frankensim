//! Explicit nonlinear source replacement at the original native coordinates.

use super::*;

// Each original polarization potential is quadratic in its OWN charge.
// Scale its larger endpoint energy by relative capacitance change, avoiding
// subtraction of large energies or nearly equal inverse capacitances.
fn pole_work(charge: f64, before: f64, after: f64) -> Result<f64, MaterialCircuitError> {
    if charge == 0.0 || before == after {
        return Ok(0.0);
    }
    let larger_energy = finite((0.5 * charge) * (charge / before.min(after)))?;
    let work = finite(larger_energy * ((before - after) / before.max(after)))?;
    if work == 0.0 {
        return Err(invalid("nonzero polarization parameter work underflowed"));
    }
    Ok(work)
}

pub(super) fn parameter_work(
    circuit: &DebyeCircuit,
    source: &ResolvedDebyeCapacitor,
    before: &MaterialDielectricCapacitor,
    after: &MaterialDielectricCapacitor,
    polarization: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, MaterialCircuitError> {
    let original = circuit.source();
    if before.branch() != after.branch()
        || before.area_m2().to_bits() != after.area_m2().to_bits()
        || before.gap_m().to_bits() != after.gap_m().to_bits()
        || original.area_m2().to_bits() != source.area_m2().to_bits()
        || original.gap_m().to_bits() != source.gap_m().to_bits()
        || original.high_frequency_property() != source.high_frequency_property()
        || original.pole_properties() != source.pole_properties()
        || original.poles().len() != source.poles().len()
    {
        return Err(invalid(
            "nonlinear dielectric replacement requires fixed branch, geometry and ordered pole roles",
        ));
    }
    let charge = circuit.instantaneous_charge_c();
    // The original native law supplies support and representability checks.
    // These endpoints are never subtracted to infer parameter work.
    before.stored_energy_j(charge)?;
    after.stored_energy_j(charge)?;
    original.stored_energy_j(0.0, polarization)?;
    source.stored_energy_j(0.0, polarization)?;
    let mut work = crate::material_circuit::capacitance::parameter_work(
        Some(before.curve()),
        Some(after.curve()),
        original.instantaneous_capacitance_f(),
        charge,
        cancelled,
    )?;
    let mut correction = 0.0;
    for ((old, new), &charge) in original
        .poles()
        .iter()
        .zip(source.poles())
        .zip(polarization)
    {
        poll(cancelled)?;
        accumulate(
            &mut work,
            &mut correction,
            pole_work(charge, old.capacitance_f, new.capacitance_f)?,
        )?;
    }
    finite(work + correction)
}

impl DebyeCircuit {
    /// Explicitly replace BOTH original nonlinear instantaneous and Debye
    /// sources at fixed geometry, pole ordering, charge history and fluxes.
    /// Each proposed source pair must describe the same card/context and pass
    /// the original instantaneous/polarization-role admission. A different
    /// material card is an explicit source change, never an implicit refresh.
    ///
    /// Signed work is integral_0^q_inf (V_new-V_old) dq plus each original
    /// pole's fixed-charge energy change. Relaxation-time changes do no stored
    /// work, and no contribution is counted as irreversible heat. The physical
    /// parameter owner must supply/receive this work when accepting the result.
    /// Incompatible ideal clamps and missing source support refuse atomically.
    pub fn rebind_reversible_sources_candidate(
        &self,
        source: ResolvedDebyeCapacitor,
        instantaneous: MaterialDielectricCapacitor,
        inputs: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<DebyeRebindCandidate, MaterialCircuitError> {
        self.rebind(source, Some(instantaneous), None, inputs, cancelled)
    }

    /// The same conservative nonlinear replacement with one complete original
    /// source binding per non-polarization resistor. Existing conductor owners
    /// retain their complete claims, geometry and temperature-only context, as
    /// in the ordinary Debye refresh. Original source/transformer topology and
    /// all unrelated storage coordinates remain unchanged.
    pub fn rebind_reversible_sources_and_resistors_candidate(
        &self,
        source: ResolvedDebyeCapacitor,
        instantaneous: MaterialDielectricCapacitor,
        resistors: Vec<MaterialResistor>,
        inputs: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<DebyeRebindCandidate, MaterialCircuitError> {
        self.rebind(
            source,
            Some(instantaneous),
            Some(resistors),
            inputs,
            cancelled,
        )
    }
}
