//! Local energy transport and molecular heat from the original acoustic field.
//!
//! For exp(-i omega t), Euler momentum and the original caloric compliance give
//! `I = Re(p conj(v))/2` and `Q = omega Im(1/K) |p|²/2`, hence `div(I) = -Q`
//! away from sources. Q is physical irreversible molecular heating per volume,
//! not outgoing radiation, wall absorption, or an amplitude-attenuation fit.

use super::*;
use crate::near_field::VectorObservation;

/// Cycle-mean transport and irreversible gas heating at one receiver.
/// The source pressure and velocity are peak phasors, not RMS values.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MolecularPowerDensity {
    /// Signed Cartesian active intensity [W/m²]. Reactive near-field velocity
    /// remains in the original field and is not replaced by a plane-wave value.
    pub active_intensity_w_m2: [f64; 3],
    /// Irreversible molecular bulk heat [W/m³]. This excludes wall losses,
    /// shear viscosity and thermal conduction, which this exterior does not solve.
    pub absorbed_power_density_w_m3: f64,
}

/// Co-located acoustic fields and their physical power observations.
pub struct MolecularPowerObservation {
    /// Original adaptive pressure/velocity and their quadrature discrepancies.
    /// Those estimates retain their existing meaning; they are not a bound on
    /// the nonlinear power values, BEM error or constitutive model error.
    pub acoustic: VectorObservation,
    /// One row per original receiver, in exactly the same order.
    pub power: Vec<MolecularPowerDensity>,
}

fn finish(
    source: &RelaxationState,
    omega: f64,
    acoustic: VectorObservation,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<MolecularPowerObservation, MolecularError> {
    poll(cancelled)?;
    let compliance = source.caloric_response(omega)?.bulk_modulus_pa.recip();
    // Scale before squaring pressure. Forming |p|² first overflows for large
    // amplitudes even when weak absorption leaves a representable heat density.
    let root_loss = det::sqrt(0.5) * det::sqrt(omega) * det::sqrt(compliance.im);
    if !root_loss.is_finite() || root_loss < 0.0 {
        return Err(bad("molecular bulk heat coefficient is unrepresentable"));
    }
    let mut power = Vec::with_capacity(acoustic.scalar.pressure.len());
    for (&pressure, &velocity) in acoustic
        .scalar
        .pressure
        .iter()
        .zip(&acoustic.particle_velocity_m_s)
    {
        poll(cancelled)?;
        power.push(density(pressure, velocity, root_loss)?);
    }
    poll(cancelled)?;
    Ok(MolecularPowerObservation { acoustic, power })
}

fn density(
    pressure: C64,
    velocity: [C64; 3],
    root_loss: f64,
) -> Result<MolecularPowerDensity, MolecularError> {
    let active_intensity_w_m2 = velocity.map(|v| (pressure.scale(0.5) * v.conj()).re);
    let absorbed_power_density_w_m3 = pressure.scale(root_loss).norm_sq();
    if !finite(pressure)
        || !velocity.into_iter().all(finite)
        || active_intensity_w_m2.iter().any(|v| !v.is_finite())
        || !absorbed_power_density_w_m3.is_finite()
    {
        return Err(bad("molecular local power is unrepresentable"));
    }
    Ok(MolecularPowerDensity {
        active_intensity_w_m2,
        absorbed_power_density_w_m3,
    })
}

impl MolecularRadiation {
    /// Resolve molecular heating and active intensity from the same adaptive
    /// outgoing pressure/velocity field. No extra propagation loss is applied.
    /// To integrate gas heating, the caller supplies physical volume quadrature
    /// and convergence independently of the BEM/receiver quadrature.
    pub fn observe_power(
        &self,
        surface: &SpherePanels,
        points_m: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: MolecularNearFieldOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularPowerObservation, MolecularError> {
        let acoustic = self.observe_near_field(
            surface,
            points_m,
            minimum_clearance_m,
            options,
            &mut cancelled,
        )?;
        finish(&self.source, self.omega, acoustic, &mut cancelled)
    }
}

impl PreparedMolecularNearField<'_> {
    /// Reuse the original source-bound receiver rows, including all geometry,
    /// frequency and material-identity checks. Power adds O(receivers) work;
    /// material evaluation occurs once per field, with no new solve/quadrature.
    pub fn evaluate_power(
        &self,
        solution: &MolecularRadiation,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularPowerObservation, MolecularError> {
        let acoustic = self.evaluate(solution, &mut cancelled)?;
        finish(&solution.source, solution.omega, acoustic, &mut cancelled)
    }
}

impl MolecularScatteredField<'_, '_> {
    /// Molecular heat and intensity of the TOTAL incident plus scattered field.
    /// Pressure and velocity are summed before quadratic work is evaluated, so
    /// incident/scattered interference is retained. Summing their individual
    /// absorption densities would lose these physical cross terms.
    pub fn total_power_at(
        &self,
        points_m: &[[f64; 3]],
        minimum_receiver_clearance_m: f64,
        options: MolecularNearFieldOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularPowerObservation, MolecularError> {
        let acoustic = self.total_pressure_velocity_adaptive_at(
            points_m,
            minimum_receiver_clearance_m,
            options,
            &mut cancelled,
        )?;
        let source = &self.solution().scattered;
        finish(&source.source, source.omega, acoustic, &mut cancelled)
    }
}

#[cfg(test)]
#[path = "power/tests.rs"]
mod tests;
