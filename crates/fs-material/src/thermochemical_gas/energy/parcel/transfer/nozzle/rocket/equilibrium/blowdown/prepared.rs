//! One original birth-isentrope workspace for external finite-outflow consumers.

use super::{EquilibriumBlowdownCandidate, EquilibriumBlowdownError, EquilibriumBlowdownSolve,
    EquilibriumBlowdownStep,
    EquilibriumGasReservoir, reservoir::{check_cancel, positive}};
use fs_thermochem::equilibrium::composition::isentrope::{
    EquilibriumIsentropeV1, density::EquilibriumDensityQueryV1,
};

/// Reusable original reservoir thermodynamics; trials never mutate the tank.
///
/// No alternate energy equation is introduced. Density samples and the final
/// finite mass/enthalpy exchange use the same private owners as rocket blowdown.
/// Access to the SP path lets a higher-layer nozzle reuse its actual chemistry
/// workspace rather than prepare a different entropy reference for each trial.
pub struct EquilibriumReservoirDepletion<'a> {
    reservoir: &'a EquilibriumGasReservoir,
    path: EquilibriumIsentropeV1<'a>,
    controls: EquilibriumBlowdownSolve,
}

impl EquilibriumGasReservoir {
    /// Prepare one source-bound depletion path for a finite physical interval.
    pub fn prepare_depletion(
        &self,
        controls: EquilibriumBlowdownSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumReservoirDepletion<'_>, EquilibriumBlowdownError> {
        check_cancel(&mut cancelled)?;
        controls.validate(self.isentropic_sample().pressure.value())?;
        let path = self.path(&mut cancelled)?;
        check_cancel(&mut cancelled)?;
        Ok(EquilibriumReservoirDepletion { reservoir: self, path, controls })
    }
}

impl<'a> EquilibriumReservoirDepletion<'a> {
    /// Original reusable chemistry. Queries may update numerical warm starts;
    /// the borrowed physical tank does not move. Replacing this path is allowed
    /// by Rust, but a changed birth reference, feed or policy is rejected before
    /// every subsequent density sample and finite endpoint publication.
    pub fn isentrope(&mut self) -> &mut EquilibriumIsentropeV1<'a> {
        &mut self.path
    }

    /// Reconstruct actual density m/V on the original chemical isentrope.
    /// Only positive masses no larger than the accepted reservoir are admitted.
    pub fn sample_remaining_mass(
        &mut self,
        mass_kg: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumDensityQueryV1, EquilibriumBlowdownError> {
        check_cancel(&mut cancelled)?;
        self.reservoir.admit_prepared_path(&self.path)?;
        if !positive(mass_kg) || mass_kg > self.reservoir.mass_kg() {
            return Err(EquilibriumBlowdownError::Invalid("depletion sample must retain a positive original tank inventory"));
        }
        self.reservoir.density_at_mass(&mut self.path, mass_kg, self.controls, &mut cancelled)
    }

    /// Evaluate the original finite U-loss account without allocating source receipts.
    /// The accepted native current sample is the remaining tank, allowing a
    /// pressure-driven consumer to reuse it for the following nozzle query.
    pub fn sample_withdrawal(
        &mut self,
        mass_kg: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumBlowdownStep, EquilibriumBlowdownError> {
        check_cancel(&mut cancelled)?;
        self.reservoir.admit_prepared_path(&self.path)?;
        self.reservoir.step_at_mass(
            &mut self.path, mass_kg, self.controls, &mut cancelled,
        )
    }

    /// Admit the finite endpoint and its integrated stagnation enthalpy once.
    /// The endpoint is freshly reconstructed, so an intervening nozzle query
    /// cannot leave a throat state masquerading as the remaining reservoir.
    pub fn finish_remaining_mass(
        &mut self,
        mass_kg: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumBlowdownCandidate, EquilibriumBlowdownError> {
        check_cancel(&mut cancelled)?;
        self.reservoir.admit_prepared_path(&self.path)?;
        self.reservoir.candidate_at_mass(&mut self.path, mass_kg, self.controls, &mut cancelled)
    }
}
