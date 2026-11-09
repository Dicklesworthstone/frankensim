//! Neutral-solute saturation with finite solid and liquid inventories.
//!
//! The explicitly supplied concentration is mol per m³ of incompressible
//! solvent, with negligible solute volume. A constant dissolution enthalpy gives
//! `ln(c_sat(T)/c_ref) = -ΔH/R (1/T - 1/T_ref)`. Optional original solid and
//! dissolved Cp(T) curves instead give `ΔH(T) = ΔH_ref + ∫ ΔCp dT` and integrate
//! the SAME van't Hoff equation `d ln(c_sat)/dT = ΔH(T)/(R T²)`. This is a
//! caller-qualified ideal-dilute pair law on a bounded temperature interval, not a solubility
//! inferred from a material name. Endothermic dissolution increases solubility
//! with temperature. No heat is added to a finite thermal owner by this law.
//!
//! At equilibrium `n_liquid = min(n_total, c_sat V_solvent)`; the remaining
//! material is a pure solid. An absent solid admits undersaturation, while a
//! present solid fixes the dissolved concentration. Retaining the total lets a
//! consumer precipitate on solvent loss and dissolve that SAME stock on dilution.
//! The finite-phase inequality follows the USGS PHREEQC equilibrium-phase
//! convention: <https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqc3-html/phreeqc3-13.htm>.
//!
//! This law is for independent, nonreacting neutral solutes. For ionic solids,
//! van't Hoff applies to the activity product Ksp, not to salt concentration:
//! <https://goldbook.iupac.org/terms/view/S05742>. Common ions, pH, complexes,
//! hydrates, competing solids and nonideal solutions require their original
//! reaction/speciation owners. No nucleation barrier, metastable supersaturation,
//! dissolution rate, solid volume, or caloric evolution is supplied here.

pub mod caloric;

/// Finite-rate dissolution and precipitation with fixed or evolving particle area.
pub mod kinetics;

mod heat_capacity;
mod potential;

use crate::solution::caloric::MolarHeatCapacityCurveV1;
use crate::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use fs_math::det::{exp, ln};
use fs_qty::{Amount, AmountConcentration, Temperature, Volume};
use heat_capacity::PhaseHeatCapacities;
use std::sync::Arc;

/// Invalid pair/model inputs, unsupported temperature, or lost numerical range.
#[derive(Clone, Debug, PartialEq)]
pub enum SolubilityErrorV1 {
    /// A required physical quantity or pair identity is invalid.
    Invalid(&'static str),
    /// Temperature lies outside the explicitly supplied inclusive interval.
    OutsideTemperatureDomain,
    /// A positive physical quantity or a finite balance cannot be represented.
    Unrepresentable(&'static str),
}
impl core::fmt::Display for SolubilityErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite neutral-solute solubility: {self:?}")
    }
}
impl std::error::Error for SolubilityErrorV1 {}

/// One declared solvent/neutral-solute pair and its bounded saturation fit.
///
/// The caller owns experimental applicability and the selected pure solid.
/// A numerical constructor is not a measured-data or material-source receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct IdealSolubilityV1 {
    solvent: SpeciesId,
    solute: SpeciesId,
    reference_temperature: Temperature,
    reference_solubility: AmountConcentration,
    dissolution_enthalpy_j_mol: f64,
    temperature_domain: [Temperature; 2],
    phase_heat_capacities: Option<Arc<PhaseHeatCapacities>>,
}
impl IdealSolubilityV1 {
    /// Bind an ideal-dilute saturation fit with an explicit concentration basis.
    ///
    /// `reference_solubility` is mol/m³ of solvent. `dissolution_enthalpy_j_mol`
    /// is the constant van't Hoff slope parameter [J/mol] for this convention;
    /// no molarity/molality conversion or solvent-density correction is inferred.
    /// The positive reference temperature must lie in the inclusive domain.
    pub fn new(
        solvent: SpeciesId,
        solute: SpeciesId,
        reference_temperature: Temperature,
        reference_solubility: AmountConcentration,
        dissolution_enthalpy_j_mol: f64,
        temperature_domain: [Temperature; 2],
    ) -> Result<Self, SolubilityErrorV1> {
        let [low, high] = temperature_domain.map(Temperature::value);
        let reference = reference_temperature.value();
        if solvent == solute
            || !positive(low)
            || !positive(high)
            || high < low
            || !positive(reference)
            || reference < low
            || reference > high
            || !positive(reference_solubility.value())
            || !dissolution_enthalpy_j_mol.is_finite()
        {
            return Err(SolubilityErrorV1::Invalid(
                "distinct pair, positive saturation and bounded temperature fit required",
            ));
        }
        Ok(Self {
            solvent,
            solute,
            reference_temperature,
            reference_solubility,
            dissolution_enthalpy_j_mol,
            temperature_domain,
            phase_heat_capacities: None,
        })
    }

    /// Bind explicitly phase-qualified solid and dissolved molar Cp(T) curves.
    ///
    /// Both must name this solute and cover the FULL original temperature
    /// domain. Their independent sensible datums do not replace this law's
    /// reference dissolution enthalpy. The caller supplies phase/source
    /// qualification; a numerical curve does not establish measured authority.
    /// Preparation merges their knots once. Subsequent evaluation allocates
    /// nothing, and cloning retains immutable shared source data.
    pub fn with_phase_heat_capacities(
        mut self,
        solid: MolarHeatCapacityCurveV1,
        dissolved: MolarHeatCapacityCurveV1,
    ) -> Result<Self, SolubilityErrorV1> {
        let phases =
            PhaseHeatCapacities::new(&self.solute, self.temperature_domain, solid, dissolved)?;
        self.phase_heat_capacities = Some(Arc::new(phases));
        Ok(self)
    }

    /// Original phase curves in fixed `[solid, dissolved]` order, when selected.
    pub fn phase_heat_capacities(&self) -> Option<&[MolarHeatCapacityCurveV1; 2]> {
        self.phase_heat_capacities
            .as_ref()
            .map(|phases| phases.curves())
    }

    /// Original solvent identity, never guessed from a phase or density.
    pub const fn solvent(&self) -> &SpeciesId {
        &self.solvent
    }
    /// Original neutral dissolved component and its pure-solid equivalent.
    pub const fn solute(&self) -> &SpeciesId {
        &self.solute
    }
    /// Inclusive source/model temperature limits [K].
    pub const fn temperature_domain(&self) -> [Temperature; 2] {
        self.temperature_domain
    }
    /// Original concentration reference temperature [K].
    pub const fn reference_temperature(&self) -> Temperature {
        self.reference_temperature
    }
    /// Original saturation concentration [mol/m³ solvent].
    pub const fn reference_solubility(&self) -> AmountConcentration {
        self.reference_solubility
    }
    /// Declared dissolution enthalpy at the original reference temperature [J/mol].
    /// It is constant only when no unequal phase heat capacities are selected.
    pub const fn dissolution_enthalpy_j_mol(&self) -> f64 {
        self.dissolution_enthalpy_j_mol
    }

    /// Dissolved-minus-solid molar enthalpy at the actual temperature [J/mol].
    /// The same original Cp difference controls saturation and caloric feedback.
    pub fn dissolution_enthalpy_at(
        &self,
        temperature: Temperature,
    ) -> Result<f64, SolubilityErrorV1> {
        self.admit_temperature(temperature)?;
        let increment = match &self.phase_heat_capacities {
            Some(phases) => phases
                .enthalpy_increment(self.reference_temperature.value(), temperature.value())?,
            None => 0.0,
        };
        if increment == 0.0 {
            return Ok(self.dissolution_enthalpy_j_mol);
        }
        let enthalpy = self.dissolution_enthalpy_j_mol + increment;
        if !enthalpy.is_finite() {
            return Err(SolubilityErrorV1::Unrepresentable(
                "temperature-dependent dissolution enthalpy",
            ));
        }
        Ok(enthalpy)
    }

    /// Direct change in the original dissolved-minus-solid molar enthalpy
    /// [J/mol]. Integrates the SAME phase-Cp difference, retaining small thermal
    /// increments that subtraction of large absolute enthalpies can erase.
    pub fn dissolution_enthalpy_change_j_mol(
        &self,
        from: Temperature,
        to: Temperature,
    ) -> Result<f64, SolubilityErrorV1> {
        self.admit_temperature(from)?;
        self.admit_temperature(to)?;
        match &self.phase_heat_capacities {
            Some(phases) => phases.enthalpy_increment(from.value(), to.value()),
            None => Ok(0.0),
        }
    }

    /// Direct `ln(c_sat(to)/c_sat(from))` from the original integrated van't
    /// Hoff/Kirchhoff law. No ratio or difference of rounded concentrations is
    /// formed, so neighboring temperatures retain their small response.
    pub fn log_saturation_ratio(
        &self,
        from: Temperature,
        to: Temperature,
    ) -> Result<f64, SolubilityErrorV1> {
        self.admit_temperature(from)?;
        self.admit_temperature(to)?;
        if from == to {
            return Ok(0.0);
        }
        let enthalpy = self.dissolution_enthalpy_at(from)?;
        let mut change =
            (enthalpy / R) * ((to.value() - from.value()) / from.value() / to.value());
        if let Some(phases) = &self.phase_heat_capacities {
            change += phases.log_saturation_correction(from.value(), to.value())? / R;
        }
        if !change.is_finite() {
            return Err(SolubilityErrorV1::Unrepresentable(
                "original saturation logarithmic increment",
            ));
        }
        Ok(change)
    }

    /// Retained identities and phase-Cp source/preparation storage, excluding
    /// inline law fields. Shared immutable allocations are counted in full.
    pub fn heap_bytes(&self) -> usize {
        self.solvent.heap_bytes()
            + self.solute.heap_bytes()
            + self.phase_heat_capacities.as_ref().map_or(0, |phases| {
                core::mem::size_of::<PhaseHeatCapacities>()
                    + 2 * core::mem::size_of::<usize>()
                    + core::mem::align_of::<PhaseHeatCapacities>()
                    + phases.heap_bytes()
            })
    }

    /// Saturation at the actual temperature, with no unqualified extrapolation.
    ///
    /// This scalar evaluation allocates nothing. Logarithmic evaluation avoids
    /// overflowing the exponential when the final concentration is representable.
    /// The original reference value and constant/equal-Cp limits remain bit-exact.
    pub fn saturation_concentration(
        &self,
        temperature: Temperature,
    ) -> Result<AmountConcentration, SolubilityErrorV1> {
        self.admit_temperature(temperature)?;
        let t = temperature.value();
        let reference_t = self.reference_temperature.value();
        if t == reference_t {
            return Ok(self.reference_solubility);
        }
        // (T-Tref)/Tref/T avoids both subtracting nearly equal reciprocals
        // and forming a potentially overflowing product T*Tref.
        let mut exponent = if self.dissolution_enthalpy_j_mol == 0.0 {
            0.0
        } else {
            (self.dissolution_enthalpy_j_mol / R) * ((t - reference_t) / reference_t / t)
        };
        if let Some(phases) = &self.phase_heat_capacities {
            exponent += phases.log_saturation_correction(reference_t, t)? / R;
        }
        if !exponent.is_finite() {
            return Err(SolubilityErrorV1::Unrepresentable(
                "temperature-dependent saturation exponent",
            ));
        }
        if exponent == 0.0 {
            return Ok(self.reference_solubility);
        }
        let concentration = exp(ln(self.reference_solubility.value()) + exponent);
        if !positive(concentration) {
            return Err(SolubilityErrorV1::Unrepresentable(
                "positive saturation concentration",
            ));
        }
        Ok(AmountConcentration::new(concentration))
    }

    fn admit_temperature(&self, temperature: Temperature) -> Result<(), SolubilityErrorV1> {
        let t = temperature.value();
        if !t.is_finite()
            || t < self.temperature_domain[0].value()
            || t > self.temperature_domain[1].value()
        {
            return Err(SolubilityErrorV1::OutsideTemperatureDomain);
        }
        Ok(())
    }

    /// Equilibrate the complete component stock at the actual solvent volume.
    /// `total_solute` includes both dissolved and solid moles; a solid amount
    /// must never be dropped when a subsequent state is constructed.
    pub fn equilibrate(
        &self,
        temperature: Temperature,
        solvent_volume: Volume,
        total_solute: Amount,
    ) -> Result<SolubilityPartitionV1, SolubilityErrorV1> {
        partition_at_saturation(
            self.saturation_concentration(temperature)?,
            solvent_volume,
            total_solute,
        )
    }
}

/// Finite equilibrium on one explicitly supplied saturation/concentration basis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolubilityPartitionV1 {
    /// Complete dissolved-plus-solid component stock [mol].
    pub total_solute: Amount,
    /// Actual incompressible solvent volume [m³], excluding solid/solute volume.
    pub solvent_volume: Volume,
    /// Saturation concentration at the owner's actual temperature [mol/m³].
    pub saturation_concentration: AmountConcentration,
    /// Dissolved capacity of this solvent inventory [mol].
    pub capacity: Amount,
    /// Actual dissolved amount [mol], retained even against a much larger solid.
    pub dissolved: Amount,
    /// Remaining pure solid [mol], available for subsequent redissolution.
    pub solid: Amount,
    /// Compensated dissolved-plus-solid-minus-total residual [mol].
    /// This is a floating-point balance, not an experimental uncertainty bound.
    pub balance_residual: Amount,
}

/// Partition at an already evaluated saturation concentration.
///
/// Fixed-temperature material consumers can cache the original pair-law
/// evaluation and call this allocation-free, noniterative helper at each volume
/// trial. It preserves the finite component total, positive dissolved traces and
/// the absent-solid undersaturation inequality without a concentration floor.
pub fn partition_at_saturation(
    saturation: AmountConcentration,
    solvent_volume: Volume,
    total_solute: Amount,
) -> Result<SolubilityPartitionV1, SolubilityErrorV1> {
    let total = total_solute.value();
    if !positive(saturation.value())
        || !positive(solvent_volume.value())
        || !total.is_finite()
        || total < 0.0
    {
        return Err(SolubilityErrorV1::Invalid(
            "positive saturation/solvent volume and nonnegative total required",
        ));
    }
    let capacity = saturation.value() * solvent_volume.value();
    if !positive(capacity) {
        return Err(SolubilityErrorV1::Unrepresentable(
            "positive dissolved capacity",
        ));
    }
    let dissolved = total.min(capacity);
    let solid = total - dissolved;
    // For a dominant solid stock the rounded sum may equal total exactly while
    // the separately retained dissolved trace is real. Recover the addition's
    // low part before subtracting the original total (Neumaier/two-sum).
    let sum = dissolved + solid;
    let low = if dissolved >= solid {
        (dissolved - sum) + solid
    } else {
        (solid - sum) + dissolved
    };
    let residual = (sum - total) + low;
    if !sum.is_finite() || !residual.is_finite() || residual.abs() > 8.0 * f64::EPSILON * total {
        return Err(SolubilityErrorV1::Unrepresentable(
            "dissolved and solid component balance",
        ));
    }
    Ok(SolubilityPartitionV1 {
        total_solute,
        solvent_volume,
        saturation_concentration: saturation,
        capacity: Amount::new(capacity),
        dissolved: Amount::new(dissolved),
        solid: Amount::new(solid),
        balance_residual: Amount::new(residual),
    })
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

#[cfg(test)]
mod tests;
