//! Durable open-parcel coordinates on an immutable native source catalog.
//!
//! Store accepted mass and the declared composition directly. Reconstructing
//! either from rounded species masses would change subsequent transport. The
//! native material and mixture owners retain all source data and evaluate the
//! restored point; this module supplies no caloric, mixing or flow law.

use fs_matdb::QueryPoint;
use fs_thermochem::{Composition, CompositionBasis};

use super::super::{GasEnergyError, IsochoricGasParcel, positive, updated_point};
use super::ClosedGasParcel;
use super::transfer::species::boundary::admit_catalog;

const VERSION: u64 = u64::from_le_bytes(*b"gasopen1");
const HEADER_WORDS: usize = 6;

fn auxiliary(point: &QueryPoint) -> impl Iterator<Item = (&str, u64)> {
    point
        .axes()
        .iter()
        .filter(|(axis, _)| axis.as_str() != "temperature" && axis.as_str() != "pressure")
        .map(|(axis, value)| (axis.as_str(), value.to_bits()))
}

fn invalid(message: &'static str) -> GasEnergyError {
    GasEnergyError::Invalid(message)
}

impl ClosedGasParcel {
    /// Admit a variable inventory on this parcel's original source catalog.
    ///
    /// Active species may disappear or return, and mass, volume, temperature
    /// and pressure may change. Every retained caloric law, transport source,
    /// source-pack receipt, explicit water identity and auxiliary query axis
    /// must remain the reference's. This is source admission, not a claim that
    /// any particular transport or energy balance has closed.
    pub fn admit_inventory_state(&self, state: &Self) -> Result<(), GasEnergyError> {
        admit_catalog(state, self.model())?;
        if self.point().axis_quantities() != state.point().axis_quantities()
            || !auxiliary(self.point()).eq(auxiliary(state.point()))
        {
            return Err(invalid(
                "open gas state changes auxiliary material conditions",
            ));
        }
        Ok(())
    }

    /// Append a source-bound open state, leaving `words` unchanged on refusal.
    ///
    /// The version, accepted mass, volume, temperature, pressure and declared
    /// fraction basis precede one fraction per reference species. Exact zero
    /// slots preserve depleted species. Keeping the original declared basis
    /// and fraction bits also preserves the native model's evaluation receipts.
    /// The reference's complete identity must be bound by the enclosing codec.
    pub fn encode_inventory_state(
        &self,
        state: &Self,
        words: &mut Vec<u64>,
    ) -> Result<(), GasEnergyError> {
        self.admit_inventory_state(state)?;
        let catalog = self.model().thermodynamics().components();
        let active = state.model().thermodynamics();
        let declared = active.declared_composition();
        let basis = match declared.basis() {
            CompositionBasis::MassFraction => 0,
            CompositionBasis::MoleFraction => 1,
            CompositionBasis::VolumeFraction => {
                return Err(invalid("open gas snapshot requires a mass or mole basis"));
            }
        };
        let fractions: Vec<f64> = catalog
            .iter()
            .map(|source| {
                active
                    .components()
                    .iter()
                    .position(|held| held.species() == source.species())
                    .map_or(0.0, |index| declared.fractions()[index])
            })
            .collect();
        let composition = Composition::new(declared.basis(), fractions.clone())
            .map_err(|error| crate::thermochemical_gas::failure("inventory snapshot", error))?;
        // A codec must reproduce the complete model, including receipt order,
        // rather than silently canonicalizing a different accepted object.
        if self.model().with_composition(composition)? != *state.model() {
            return Err(invalid(
                "open gas snapshot cannot reproduce its complete source model",
            ));
        }
        let resolved = state.state().state();
        words.extend([
            VERSION,
            state.mass_kg().to_bits(),
            state.volume_m3().to_bits(),
            resolved.temperature.to_bits(),
            resolved.pressure.to_bits(),
            basis,
        ]);
        words.extend(fractions.into_iter().map(f64::to_bits));
        Ok(())
    }

    /// Rebuild an open state through the original mixture and material owners.
    ///
    /// The stored pressure is evaluated exactly, preserving an accepted
    /// endpoint's rounding. Its native density must agree with stored mass / V
    /// within the existing parcel tolerance of 64 eps. Unsupported source
    /// states, invalid fractions, vacuum, nonfinite coordinates, changed word
    /// counts/versions and noncanonical zero slots refuse. Nothing is published.
    pub fn decode_inventory_state(&self, words: &[u64]) -> Result<Self, GasEnergyError> {
        let count = self.model().thermodynamics().components().len();
        if words.len() != HEADER_WORDS + count || words[0] != VERSION {
            return Err(invalid("open gas snapshot version or word count differs"));
        }
        let [mass, volume, temperature, pressure] =
            [words[1], words[2], words[3], words[4]].map(f64::from_bits);
        let density = mass / volume;
        if ![mass, volume, temperature, pressure, density]
            .into_iter()
            .all(positive)
        {
            return Err(invalid(
                "open gas snapshot needs positive finite physical coordinates",
            ));
        }
        let basis = match words[5] {
            0 => CompositionBasis::MassFraction,
            1 => CompositionBasis::MoleFraction,
            _ => {
                return Err(invalid(
                    "open gas snapshot has an unsupported fraction basis",
                ));
            }
        };
        let composition = Composition::new(
            basis,
            words[HEADER_WORDS..]
                .iter()
                .copied()
                .map(f64::from_bits)
                .collect(),
        )
        .map_err(|error| crate::thermochemical_gas::failure("inventory snapshot", error))?;
        let model = self.model().with_composition(composition)?;
        let point = updated_point(self.point(), temperature, pressure)?;
        let resolved = model.evaluate(&point)?;
        if !((resolved.state().density / density - 1.0).abs() <= 64.0 * f64::EPSILON) {
            return Err(invalid(
                "open gas snapshot density disagrees with its mass and volume",
            ));
        }
        let state = Self::from_isochoric(IsochoricGasParcel {
            model,
            point,
            state: resolved,
            mass_kg: mass,
            volume_m3: volume,
        });
        let mut rebuilt = Vec::with_capacity(words.len());
        self.encode_inventory_state(&state, &mut rebuilt)?;
        if rebuilt != words {
            return Err(invalid("open gas snapshot coordinates are not canonical"));
        }
        Ok(state)
    }
}
