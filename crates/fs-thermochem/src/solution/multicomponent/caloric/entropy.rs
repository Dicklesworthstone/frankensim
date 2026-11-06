//! Original species-reduced entropy for homogeneous liquid transport.

use super::{MulticomponentCaloricStateV1, Result, Sum, checked, scaled};

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
        let mut entropy = Sum::default();
        for (&amount, source) in self.component_moles.iter().zip(&self.owner.heat_capacity) {
            if amount != 0.0 {
                entropy.add(scaled(
                    amount,
                    source.entropy_change_j_per_mol_k(
                        source.reference_temperature(),
                        self.temperature(),
                    )?,
                    "source-relative homogeneous sensible entropy",
                )?);
            }
        }
        entropy.add(scaled(
            self.total_moles,
            self.mixing.mixing_entropy_j_per_mol_k(),
            "source-relative homogeneous mixing entropy",
        )?);
        checked(entropy.value(), "source-relative homogeneous entropy")
    }
}
