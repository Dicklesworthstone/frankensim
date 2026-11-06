//! Nonideal liquid enthalpy on the existing stationary-species flow graph.
//!
//! Each cell keeps its admitted stationary composition while the native thermal
//! owner evolves enthalpy. Inlets carry THEIR original mixture enthalpy in the
//! same pure-species datum; receiving cells recover temperature from their own
//! equilibrium calorics. Heat of mixing and temperature-dependent Cp therefore
//! enter exactly once, without temperature averaging or a separate mixer solver.
//!
//! The same storage also binds the existing implicit conduction adapter. This
//! is prescribed mass-balanced bulk flow at common fixed pressure, not transient
//! species storage, phase-selective flow, a momentum solve or an inferred fluid
//! geometry. Original pure Cp/T integrals and accepted phase mixing entropy
//! supply the native flow entropy account without inferred absolute species
//! entropies or a complete pressure/chemical-work model.

use super::{
    BinaryCaloricSolveV1, BinaryLiquidCaloricError, BinaryLiquidCaloricLaw,
    BinaryLiquidThermalStorage, SourceOwner, poll,
};
use fs_blake3::ContentHash;
use fs_thermal::enthalpy::advection::species::{CaloricComposition, CaloricEntropyAnchor};

#[derive(Clone, Debug, PartialEq)]
pub(super) struct StreamBinding {
    species: [String; 2],
    mass_fractions: [f64; 2],
    enthalpy_reference: String,
    entropy_reference: String,
    entropy_anchor_h_j_kg: f64,
    reduced_entropy_anchor_j_kg_k: f64,
    caloric_identity: [u8; 32],
}

impl StreamBinding {
    pub(super) fn composition(&self) -> CaloricComposition<'_> {
        CaloricComposition {
            species: &self.species,
            mass_fractions: &self.mass_fractions,
            enthalpy_reference: &self.enthalpy_reference,
        }
    }

    pub(super) const fn caloric_identity(&self) -> [u8; 32] {
        self.caloric_identity
    }

    pub(super) fn entropy_anchor(&self) -> CaloricEntropyAnchor<'_> {
        CaloricEntropyAnchor {
            specific_h_j_kg: self.entropy_anchor_h_j_kg,
            reduced_specific_entropy_j_kg_k: self.reduced_entropy_anchor_j_kg_k,
            entropy_reference: &self.entropy_reference,
        }
    }
}

impl BinaryLiquidCaloricLaw<'_, '_> {
    /// Prepare a frozen-composition liquid stream on the original common
    /// pure-species enthalpy datum. Both original pure species have zero sensible
    /// enthalpy at their admitted common datum; the original mixing enthalpy is
    /// retained there, including any two-liquid equilibrium.
    ///
    /// Fractions are MASS fractions in the original pair order and must sum to
    /// exactly one. They are retained without a mole-to-mass round trip so an
    /// existing `StationarySpecies` result can bind directly. Component amounts
    /// use the original molar masses; a lost positive trace is refused.
    ///
    /// The common reference binds both pure Cp source identities, species,
    /// molar masses, static pressure and the complete interaction law. It is
    /// independent of mixture composition and numerical solver budgets. A
    /// separate storage identity includes this composition and temperature
    /// domain. An unchanged source may therefore supply every inlet and mixer
    /// cell while retaining different physical Cp(T) and enthalpy curves.
    ///
    /// The existing implicit thermal graph admits these laws through
    /// `Advected::with_stationary_species` and `with_implicit_advection`.
    /// `with_implicit_advection_entropy` additionally retains original mixing
    /// and Cp/T contributions on the same accepted flow stages. Its anchor
    /// removes only shared linear pure-species entropy constants at the common
    /// source datum. The entropy convention is independent of composition and
    /// prepared temperature bracket; it does not infer absolute entropy or a
    /// pressure/chemical-work closure from heat capacity.
    pub fn stream_storage_from_mass_fractions(
        &self,
        mass_fractions: [f64; 2],
        controls: BinaryCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryLiquidThermalStorage, BinaryLiquidCaloricError> {
        poll(&mut cancelled)?;
        if mass_fractions
            .iter()
            .any(|&y| !y.is_finite() || !(0.0..=1.0).contains(&y))
            || mass_fractions[0] + mass_fractions[1] != 1.0
        {
            return Err(BinaryLiquidCaloricError::Invalid(
                "stream storage requires exact unit mass fractions in original species order",
            ));
        }
        let mut component_moles = [0.0; 2];
        for i in 0..2 {
            component_moles[i] = mass_fractions[i] / self.owner.liquid.molar_masses[i].value();
            if !component_moles[i].is_finite()
                || (mass_fractions[i] > 0.0 && component_moles[i] <= 0.0)
            {
                return Err(BinaryLiquidCaloricError::Invalid(
                    "stream mass-to-mole conversion loses a positive component",
                ));
            }
        }
        let mut storage = self.thermal_storage(component_moles, controls, &mut cancelled)?;
        let lower = storage.endpoints[0].enthalpy_j() / storage.mass_kg;
        let upper = lower + storage.upper_specific_h_j_kg;
        if !lower.is_finite() || !upper.is_finite() || upper <= lower {
            return Err(BinaryLiquidCaloricError::Invalid(
                "stream enthalpy domain is unrepresentable in the original common datum",
            ));
        }
        let reduced_entropy_anchor_j_kg_k =
            storage.endpoints[0].reduced_entropy_j_per_k()? / storage.mass_kg;
        if !reduced_entropy_anchor_j_kg_k.is_finite() {
            return Err(BinaryLiquidCaloricError::Invalid(
                "source-relative stream entropy is not representable",
            ));
        }
        let reference = common_reference(&self.owner);
        let mut identity = Vec::with_capacity(64);
        identity.extend_from_slice(reference.as_bytes());
        for value in mass_fractions
            .into_iter()
            .chain(controls.temperature_bracket.map(|t| t.value()))
        {
            identity.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        storage.lower_specific_h_j_kg = lower;
        storage.upper_specific_h_j_kg = upper;
        storage.stream = Some(StreamBinding {
            species: self
                .owner
                .model
                .thermodynamics()
                .species()
                .each_ref()
                .map(|species| species.to_string()),
            mass_fractions,
            enthalpy_reference: format!("binary-liquid-isobaric-v1:{}", reference.to_hex()),
            entropy_reference: format!(
                "binary-liquid-species-reduced-entropy-v1:{}",
                reference.to_hex()
            ),
            entropy_anchor_h_j_kg: lower,
            reduced_entropy_anchor_j_kg_k,
            caloric_identity: *fs_blake3::hash_domain(
                "org.frankensim.fs-material.binary-liquid-stream-storage.v1",
                &identity,
            )
            .as_bytes(),
        });
        poll(&mut cancelled)?;
        Ok(storage)
    }
}

// Encode the original caloric data, not a material label, sampled answer or
// current composition. Pure source identities already retain their datum,
// complete knots, original context and receipts. Extra transport properties
// and a mixture card's reference composition do not change the energy datum.
fn common_reference(owner: &SourceOwner) -> ContentHash {
    let pair = owner.model.thermodynamics();
    let mut bytes = Vec::new();
    for i in 0..2 {
        let name = pair.species()[i].as_str().as_bytes();
        bytes.extend_from_slice(&(name.len() as u64).to_le_bytes());
        bytes.extend_from_slice(name);
        bytes.extend_from_slice(&owner.liquid.molar_masses[i].value().to_bits().to_le_bytes());
        bytes.extend_from_slice(owner.sources[i].identity().as_bytes());
    }
    for value in [
        owner.sources[0].pressure_pa(),
        pair.chi_constant(),
        pair.chi_inverse_temperature_k(),
        pair.excess_heat_capacity_coefficient_j_per_mol_k(),
        pair.excess_heat_capacity_reference_temperature().value(),
        pair.temperature_domain()[0].value(),
        pair.temperature_domain()[1].value(),
    ] {
        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    fs_blake3::hash_domain(
        "org.frankensim.fs-material.binary-liquid-caloric-reference.v1",
        &bytes,
    )
}
