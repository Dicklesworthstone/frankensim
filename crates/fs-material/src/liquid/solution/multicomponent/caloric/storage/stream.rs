//! Original multicomponent liquid enthalpy and entropy on native mass flow.
//!
//! Stationary species fix each cell's composition. Inlets carry their own
//! original enthalpy, while receiving cells invert their own calorics. The
//! shared pure-species datum retains heat of mixing without a second mixer.
//! Original Cp/T and nonideal mixing entropy anchor the native entropy account.
//! No absolute species entropy, phase-selective flow or transient composition
//! storage is inferred. Pressure and reference geometry remain fixed.

use super::{
    MulticomponentCaloricSolveV1, MulticomponentLiquidCaloricError, MulticomponentLiquidCaloricLaw,
    MulticomponentLiquidThermalStorage, Result, SourceOwner, poll,
};
use fs_blake3::ContentHash;
use fs_thermal::enthalpy::advection::species::{CaloricComposition, CaloricEntropyAnchor};

#[derive(Clone, Debug, PartialEq)]
pub(super) struct StreamBinding {
    species: Box<[String]>,
    mass_fractions: Box<[f64]>,
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

impl MulticomponentLiquidCaloricLaw<'_, '_> {
    /// Prepare an explicitly homogeneous stationary stream in the original
    /// canonical species order. Authored MASS fractions must sum exactly to
    /// one; they are retained unchanged for the native stationary species owner.
    /// Original species molar masses perform the mass-to-mole conversion, and
    /// an unrepresentable positive trace refuses instead of disappearing.
    ///
    /// Source compatibility binds complete original pure Cp receipts, common
    /// datum, species, molar masses, static pressure and interaction matrices.
    /// Different compositions and thermal brackets share this reference.
    /// A separate identity binds the actual composition and prepared domain.
    /// The original species-reduced entropy removes only common linear species
    /// constants; ideal/nonideal mixing entropy and all Cp/T contributions stay.
    ///
    /// The resulting law supplies `Advected::with_stationary_species`, native
    /// implicit advection and `with_implicit_advection_entropy` through their
    /// existing interfaces. Homogeneous metastable states remain constrained;
    /// this constructor does not certify spontaneous miscibility.
    pub fn homogeneous_stream_storage_from_mass_fractions(
        &self,
        mass_fractions: &[f64],
        controls: MulticomponentCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentLiquidThermalStorage> {
        poll(&mut cancelled)?;
        if mass_fractions.len() != self.liquid.components().len()
            || mass_fractions
                .iter()
                .any(|&y| !y.is_finite() || !(0.0..=1.0).contains(&y))
            || mass_fractions.iter().sum::<f64>() != 1.0
        {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "stream storage needs exact unit mass fractions in original species order",
            ));
        }
        let mut amounts = Vec::with_capacity(mass_fractions.len());
        for (&fraction, component) in mass_fractions.iter().zip(self.liquid.components()) {
            let amount = fraction / component.molar_mass.value();
            if !amount.is_finite() || (fraction > 0.0 && amount <= 0.0) {
                return Err(MulticomponentLiquidCaloricError::Invalid(
                    "stream mass-to-mole conversion loses a positive component",
                ));
            }
            amounts.push(amount);
        }
        let mut storage = self.homogeneous_thermal_storage(&amounts, controls, &mut cancelled)?;
        let lower = storage.endpoints[0].enthalpy_j() / storage.mass_kg;
        let upper = lower + storage.upper_specific_h_j_kg;
        let reduced_entropy_anchor_j_kg_k =
            storage.endpoints[0].reduced_entropy_j_per_k()? / storage.mass_kg;
        if !lower.is_finite()
            || !upper.is_finite()
            || upper <= lower
            || !reduced_entropy_anchor_j_kg_k.is_finite()
        {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "stream enthalpy or entropy is unrepresentable in the original datum",
            ));
        }
        let reference = common_reference(&self.owner);
        let mut identity = Vec::with_capacity(32 + 8 * (mass_fractions.len() + 2));
        identity.extend_from_slice(reference.as_bytes());
        for value in mass_fractions
            .iter()
            .copied()
            .chain(controls.temperature_bracket.map(|t| t.value()))
        {
            identity.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        storage.lower_specific_h_j_kg = lower;
        storage.upper_specific_h_j_kg = upper;
        storage.stream = Some(StreamBinding {
            species: self
                .model()
                .thermodynamics()
                .species()
                .iter()
                .map(ToString::to_string)
                .collect(),
            mass_fractions: mass_fractions.into(),
            enthalpy_reference: format!(
                "homogeneous-multicomponent-liquid-isobaric-v1:{}",
                reference.to_hex()
            ),
            entropy_reference: format!(
                "homogeneous-multicomponent-liquid-reduced-entropy-v1:{}",
                reference.to_hex()
            ),
            entropy_anchor_h_j_kg: lower,
            reduced_entropy_anchor_j_kg_k,
            caloric_identity: *fs_blake3::hash_domain(
                "org.frankensim.fs-material.homogeneous-multicomponent-liquid-stream.v1",
                &identity,
            )
            .as_bytes(),
        });
        poll(&mut cancelled)?;
        Ok(storage)
    }
}

// Canonical source data determine compatibility, independently of composition,
// numerical solver budgets or the unrelated mixture transport properties.
fn common_reference(owner: &SourceOwner) -> ContentHash {
    let model = owner.model.thermodynamics();
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&(model.species().len() as u64).to_le_bytes());
    for ((species, component), source) in model
        .species()
        .iter()
        .zip(&owner.liquid.components)
        .zip(&owner.sources)
    {
        let name = species.as_str().as_bytes();
        bytes.extend_from_slice(&(name.len() as u64).to_le_bytes());
        bytes.extend_from_slice(name);
        bytes.extend_from_slice(&component.molar_mass.value().to_bits().to_le_bytes());
        bytes.extend_from_slice(source.identity().as_bytes());
    }
    for value in [
        owner.sources[0].pressure_pa(),
        model.temperature_domain()[0].value(),
        model.temperature_domain()[1].value(),
    ]
    .into_iter()
    .chain(model.chi_constants().iter().copied())
    .chain(model.chi_inverse_temperature_k().iter().copied())
    {
        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    fs_blake3::hash_domain(
        "org.frankensim.fs-material.homogeneous-multicomponent-liquid-reference.v1",
        &bytes,
    )
}
