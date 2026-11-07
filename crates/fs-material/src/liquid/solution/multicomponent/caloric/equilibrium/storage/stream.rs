//! Original equilibrium storage on the existing stationary species flow owner.
//!
//! Whole-mixture streams share the homogeneous liquid owner's original pure
//! enthalpy and species-reduced entropy references. Local phase redistribution
//! changes calorics inside each cell; phase-selective transport, transient
//! composition, pressure evolution and phase morphology require other owners.

use super::{
    MulticomponentEquilibriumCaloricSolveV1, MulticomponentLiquidCaloricError,
    MulticomponentLiquidCaloricLaw, MulticomponentLiquidEquilibriumThermalStorage,
    ResolvedMulticomponentLiquidEquilibriumCaloricState, Result, poll,
};
use crate::liquid::solution::multicomponent::caloric::storage::stream::StreamBinding;

impl MulticomponentLiquidCaloricLaw<'_, '_> {
    /// Bind an actual separated inventory to stationary whole-mixture flow.
    /// Authored canonical MASS fractions must sum exactly to one and agree
    /// with the original phase inventory to relative arithmetic roundoff.
    /// Positive traces and exact absences must agree without an absolute floor.
    /// The fractions remain unchanged for the native species balance.
    ///
    /// Each cell retains the initial local phase count/support and one common
    /// temperature. Existing bounded continuation and the physical equilibrium
    /// Cp determine its changing temperature, phase amounts and compositions.
    /// Incoming homogeneous and separated streams can share this same source
    /// convention: their actual mixing enthalpy and entropy enter exactly once.
    /// No extra heat-of-mixing source or separately integrated temperature is
    /// introduced. The original stationary species owner admits the flow graph.
    ///
    /// Complete source compatibility is independent of composition, prepared
    /// temperature bracket and phase count. A separate local identity includes
    /// the actual cold phase allocations. A source, branch, datum or arithmetic
    /// refusal, or cancellation, leaves the supplied phase state unchanged.
    pub fn equilibrium_stream_storage_from_mass_fractions(
        &self,
        initial: &ResolvedMulticomponentLiquidEquilibriumCaloricState,
        mass_fractions: &[f64],
        controls: MulticomponentEquilibriumCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentLiquidEquilibriumThermalStorage> {
        poll(&mut cancelled)?;
        self.admit_equilibrium_caloric(initial)?;
        if mass_fractions.len() != self.liquid.components().len()
            || mass_fractions
                .iter()
                .any(|&y| !y.is_finite() || !(0.0..=1.0).contains(&y))
            || mass_fractions.iter().sum::<f64>() != 1.0
        {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "equilibrium streams need exact unit mass fractions in original species order",
            ));
        }
        let mass = initial.equilibrium().total_mass_kg();
        for ((&fraction, &amount), component) in mass_fractions
            .iter()
            .zip(initial.caloric_state().component_moles())
            .zip(self.liquid.components())
        {
            poll(&mut cancelled)?;
            let original = (amount * component.molar_mass.value()) / mass;
            if !original.is_finite()
                || (fraction == 0.0) != (amount == 0.0)
                || (amount > 0.0 && original <= 0.0)
                || (fraction - original).abs()
                    > 128.0 * f64::EPSILON * fraction.abs().max(original.abs())
            {
                return Err(MulticomponentLiquidCaloricError::Invalid(
                    "stream mass fractions differ from the original equilibrium inventory",
                ));
            }
        }
        let mut storage = self.equilibrium_thermal_storage(initial, controls, &mut cancelled)?;
        let cold = &storage.endpoints[0];
        let lower = cold.enthalpy_j() / storage.mass_kg;
        let upper = lower + storage.upper_specific_h_j_kg;
        let entropy = cold.reduced_entropy_j_per_k()? / storage.mass_kg;
        if !lower.is_finite() || !upper.is_finite() || upper <= lower || !entropy.is_finite() {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "equilibrium stream enthalpy or entropy is unrepresentable in the common datum",
            ));
        }
        let mut branch = b"equilibrium".to_vec();
        branch.extend_from_slice(&(cold.candidate().phases().len() as u64).to_le_bytes());
        for phase in cold.candidate().phases() {
            for amount in phase.component_moles() {
                branch.extend_from_slice(&amount.to_bits().to_le_bytes());
            }
        }
        storage.lower_specific_h_j_kg = lower;
        storage.upper_specific_h_j_kg = upper;
        storage.stream = Some(StreamBinding::new(
            &self.owner,
            mass_fractions,
            storage.temperature_domain(),
            [lower, entropy],
            &branch,
        ));
        poll(&mut cancelled)?;
        Ok(storage)
    }
}
