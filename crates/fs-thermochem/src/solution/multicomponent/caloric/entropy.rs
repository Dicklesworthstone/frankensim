//! Original species-reduced entropy for homogeneous and separated liquid flow.

use super::equilibrium::MulticomponentEquilibriumCaloricStateV1;
use super::{MulticomponentCaloricStateV1, Result, Sum, checked, scaled};
use crate::solution::caloric::PureLiquidHeatCapacityV1;
use fs_qty::Temperature;

impl MulticomponentCaloricStateV1 {
    /// Homogeneous entropy after the shared pure-species constants at the
    /// original common datum are removed [J/K]. Original Cp/T integrals and
    /// ideal/nonideal mixing entropy remain. Conserved component balances
    /// cancel the omitted constants analytically across different compositions.
    ///
    /// This is not an absolute entropy or a pressure-dependent equation of
    /// state. For small same-inventory thermal changes use the model's direct
    /// `entropy_change_j_per_k` rather than subtracting anchored values.
    pub fn reduced_entropy_j_per_k(&self) -> Result<f64> {
        let mut entropy = sensible_entropy(
            &self.component_moles,
            &self.owner.heat_capacity,
            self.temperature(),
        )?;
        entropy.add(scaled(
            self.total_moles,
            self.mixing.mixing_entropy_j_per_mol_k(),
            "source-relative homogeneous mixing entropy",
        )?);
        checked(entropy.value(), "source-relative homogeneous entropy")
    }
}

impl MulticomponentEquilibriumCaloricStateV1 {
    /// Entropy of the actual accepted phases after removing only the shared
    /// pure-species constants at the common datum [J/K]. Pure Cp/T integrals
    /// use the conserved inventory; each phase contributes its own mixing
    /// entropy with its actual amount and composition. A homogeneous feed's
    /// mixing entropy cannot substitute for these separated phases.
    ///
    /// This shares the homogeneous owner's entropy convention, so stationary
    /// species conservation cancels the omitted constants in an open flow
    /// balance. It supplies no absolute entropy, pressure response, phase
    /// selection or kinetics. For small same-inventory increments use the
    /// model's direct `equilibrium_entropy_change_j_per_k` instead.
    pub fn reduced_entropy_j_per_k(&self) -> Result<f64> {
        let mut entropy = sensible_entropy(
            self.component_moles(),
            self.heat_capacities(),
            self.temperature(),
        )?;
        for phase in self.candidate().phases() {
            entropy.add(scaled(
                phase.total_moles(),
                phase.thermodynamics().mixing_entropy_j_per_mol_k(),
                "source-relative equilibrium phase mixing entropy",
            )?);
        }
        checked(entropy.value(), "source-relative equilibrium entropy")
    }
}

fn sensible_entropy(
    amounts: &[f64],
    sources: &[PureLiquidHeatCapacityV1],
    temperature: Temperature,
) -> Result<Sum> {
    let mut entropy = Sum::default();
    for (&amount, source) in amounts.iter().zip(sources) {
        if amount != 0.0 {
            entropy.add(scaled(
                amount,
                source.entropy_change_j_per_mol_k(source.reference_temperature(), temperature)?,
                "source-relative liquid sensible entropy",
            )?);
        }
    }
    Ok(entropy)
}
