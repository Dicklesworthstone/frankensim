//! Sensible heat and system entropy for a prescribed uniform temperature change.
//!
//! The original regular-solution law has `chi = A + B/T`: at fixed composition
//! its mixing enthalpy and mixing entropy do not change with temperature. The
//! temperature candidate also keeps concentration, geometry and kappa fixed,
//! so its interfacial energy is unchanged. Only the original pure-liquid Cp
//! and Cp/T increments belong in this receipt. The reduced mixing free-energy
//! change remains a separate external-parameter account on the candidate.
//!
//! No phase flash, thermal root or adiabatic feedback occurs. This is the heat
//! required to realize the prescribed temperature while retaining the actual
//! spatial compositions. Pure-liquid source qualification and the maintained
//! pressure/fixed-volume approximation remain the caller's responsibility.
//! System entropy change does not specify reservoir entropy or production.

use super::super::{BinaryMiscibilityError, Sum, finite, poll};
use super::BinaryMiscibilityTemperatureCandidate;
use fs_qty::Temperature;
use fs_thermochem::solution::caloric::{BinaryCaloricErrorV1, BinaryLiquidCaloricModelV1};

/// An incompatible original pair, unsupported pure-liquid source, or an
/// incomplete caloric receipt. No error alters the temperature candidate.
#[derive(Clone, Debug, PartialEq)]
pub enum BinaryMiscibilityCaloricError {
    Invalid(&'static str),
    Solver(BinaryMiscibilityError),
    Caloric(BinaryCaloricErrorV1),
}
impl core::fmt::Display for BinaryMiscibilityCaloricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "prescribed miscibility temperature calorics: {self:?}")
    }
}
impl std::error::Error for BinaryMiscibilityCaloricError {}
impl From<BinaryMiscibilityError> for BinaryMiscibilityCaloricError {
    fn from(error: BinaryMiscibilityError) -> Self {
        Self::Solver(error)
    }
}
impl From<BinaryCaloricErrorV1> for BinaryMiscibilityCaloricError {
    fn from(error: BinaryCaloricErrorV1) -> Self {
        Self::Caloric(error)
    }
}
type Result<T> = core::result::Result<T, BinaryMiscibilityCaloricError>;

/// Original pure-liquid caloric increments for unchanged field inventories.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryMiscibilityCaloricChange {
    pub temperature_before: Temperature,
    pub temperature_after: Temperature,
    /// Signed external heat [J], positive into the liquid. This is not the
    /// candidate's reduced mixing free-energy change.
    pub heat_input_j: f64,
    /// Signed system entropy change [J/K], without an assumed heat reservoir.
    pub entropy_change_j_per_k: f64,
}

impl BinaryMiscibilityTemperatureCandidate {
    /// Account for the prescribed temperature change using the exact original
    /// pair and ordered pure-liquid Cp sources. Both endpoints must lie inside
    /// their admitted temperature domains. Spatial compositions, including
    /// metastable and unstable cells, retain their homogeneous branches.
    ///
    /// Direct Cp and Cp/T integrals avoid subtracting large reference energies.
    /// Each is evaluated once per species. Contributions retain both words of
    /// the paired cell inventories; no species is reconstructed by subtraction.
    pub fn caloric_change(
        &self,
        caloric: &BinaryLiquidCaloricModelV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityCaloricChange> {
        poll(&mut cancelled)?;
        if caloric.thermodynamics() != self.state().model().thermodynamics() {
            return Err(BinaryMiscibilityCaloricError::Invalid(
                "caloric receipt requires the original full regular-solution pair",
            ));
        }
        let change = self.change();
        let from = change.temperature_before;
        let to = change.temperature_after;
        let mut increments = [[0.0; 2]; 2];
        for (increment, source) in increments.iter_mut().zip(caloric.heat_capacities()) {
            poll(&mut cancelled)?;
            increment[0] = source.enthalpy_change_j_per_mol(from, to)?;
            increment[1] = source.entropy_change_j_per_mol_k(from, to)?;
        }
        let mut heat = Sum::default();
        let mut entropy = Sum::default();
        let state = self.state();
        for (amounts, remainders) in state.amounts.iter().zip(&state.amount_remainders) {
            poll(&mut cancelled)?;
            for species in 0..2 {
                for amount in [amounts[species], remainders[species]] {
                    heat.add(finite(
                        amount * increments[species][0],
                        "field sensible enthalpy increment",
                    )?)?;
                    entropy.add(finite(
                        amount * increments[species][1],
                        "field sensible entropy increment",
                    )?)?;
                }
            }
        }
        let result = BinaryMiscibilityCaloricChange {
            temperature_before: from,
            temperature_after: to,
            heat_input_j: heat.value()?,
            entropy_change_j_per_k: entropy.value()?,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }
}
