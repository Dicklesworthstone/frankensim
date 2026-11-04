//! One original mechanical representative releases heat into one original cell.
//!
//! The representative covers the whole fixed reference volume and mass of the
//! thermal cell. Nonuniform quadrature histories require their own aggregation;
//! one point's release is not implicitly an average over a heterogeneous cell.
//! Both caloric sources, the accepted release T/X, and fixed V0/m0 must match.
//! The original mixture specific volume must also remain unchanged, within
//! arithmetic rounding only. At the common reference pressure this makes the
//! thermal enthalpy increment equal the internal-energy increment. A genuine
//! density change needs a pressure/volume-work owner and refuses here.
//!
//! The original thermal owner applies the heat impulse without advancing time,
//! aging crystals, consuming scheduled events or rebuilding contacts. Selected
//! equilibrium remelting determines any latent response. This bridge neither
//! chooses when support is lost nor supplies fluid stress or moving geometry.
//! The unlocked state and thermal candidate must be adopted together. Failed
//! source admission, energy closure or cancellation publishes neither candidate.

use super::{CrystallizationUnlock, CrystallizationUnlockEnergy, UnlockedCrystallization};
use crate::phase::crystallization::CrystallizationCaloricPair;
use crate::phase::crystallization::mechanics::CrystallizationMechanicsError;
use fs_thermal::enthalpy::crystallization::{
    CrystallizationAudit, CrystallizationImpulseControls, CrystallizationSimulation,
    CrystallizationTimelineError,
};
use fs_thermal::enthalpy::{CellId, DepositionError};
use std::collections::BTreeMap;

/// Source/state mismatch and original thermal refusal retain separate causes.
#[derive(Clone, Debug, PartialEq)]
pub enum CrystallizationUnlockThermalError {
    Mechanics(CrystallizationMechanicsError),
    Thermal(CrystallizationTimelineError),
}

impl From<CrystallizationMechanicsError> for CrystallizationUnlockThermalError {
    fn from(error: CrystallizationMechanicsError) -> Self {
        Self::Mechanics(error)
    }
}
impl From<CrystallizationTimelineError> for CrystallizationUnlockThermalError {
    fn from(error: CrystallizationTimelineError) -> Self {
        Self::Thermal(error)
    }
}
impl From<DepositionError> for CrystallizationUnlockThermalError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(CrystallizationTimelineError::Thermal(error))
    }
}
impl core::fmt::Display for CrystallizationUnlockThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "crystallization unlock thermal transfer: {self:?}")
    }
}
impl std::error::Error for CrystallizationUnlockThermalError {}

/// A complete combined mechanical/thermal energy transfer at unchanged time.
/// The original release is consumed; this result cannot thermalize it again.
#[derive(Debug)]
pub struct ThermalizedCrystallizationUnlock {
    pub cell_id: CellId,
    pub state: UnlockedCrystallization,
    /// Original intensive storage debit and heat/export split [J/kg].
    pub energy: CrystallizationUnlockEnergy,
    /// Thermal-subsystem audit: its `external_heat_j` is the heat arriving from
    /// mechanics, an INTERNAL transfer in this combined accounting boundary.
    pub thermal_audit: CrystallizationAudit,
    /// Original elastic energy density times the fixed cell reference volume.
    pub removed_elastic_energy_j: f64,
    /// Heat requested from the thermal owner, mass times the original heat/kg.
    pub thermalized_energy_j: f64,
    /// Energy that the caller must credit to its explicitly selected other owner.
    pub exported_energy_j: f64,
    /// Thermal storage change plus export minus removed elastic storage [J].
    pub balance_residual_j: f64,
}

impl CrystallizationUnlock {
    /// Apply this authentic original-energy release to its accepted thermal cell.
    ///
    /// The same caloric pair must own the release and BOTH registered kinetic
    /// branches. Release T and crystal MASS fraction must equal the cell's
    /// current accepted state exactly. `V0/m0` must equal the actual fixed cell
    /// geometry volume divided by transported mass, with no current-free-volume
    /// substitution or tolerance large enough to hide a different reference.
    /// The same original pair must give unchanged mixture specific volume
    /// before/after the impulse, within 32 machine epsilons of the larger value.
    /// That rounding allowance is not a physical compressibility tolerance.
    ///
    /// Heat enters the original enthalpy once; the mechanical debit, actual
    /// represented thermal storage and explicit export close the combined
    /// energy account. This is not extra external heating or another latent
    /// heat source. The caller separately owns any exported energy destination.
    /// Original phase-domain, event-order and cancellation refusals propagate.
    /// The caller must adopt both returned states in its existing transaction.
    pub fn thermalize<'a>(
        self,
        pair: &CrystallizationCaloricPair<'_>,
        thermal: &CrystallizationSimulation<'a>,
        cell_id: CellId,
        controls: CrystallizationImpulseControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (
            CrystallizationSimulation<'a>,
            ThermalizedCrystallizationUnlock,
        ),
        CrystallizationUnlockThermalError,
    > {
        checkpoint()?;
        controls.validate()?;
        let invalid = CrystallizationMechanicsError::Invalid;
        let cell = thermal
            .cells()
            .get(&cell_id)
            .ok_or(invalid("unlocking requires an existing thermal cell"))?;
        let current = thermal
            .states()
            .get(&cell_id)
            .ok_or(invalid("unlocking requires an accepted kinetic cell state"))?;
        let law = thermal.material_law(cell.material_id).ok_or(invalid(
            "unlocking requires the original crystallization material",
        ))?;
        let release = self.state.release_volume_state();
        if pair.identity() != release.caloric_pair_identity
            || law.parent_caloric_identity() != Some(*pair.parent().identity().as_bytes())
            || law.crystalline_caloric_identity() != Some(*pair.crystal().identity().as_bytes())
        {
            return Err(CrystallizationMechanicsError::PairMismatch.into());
        }
        if release.temperature_k != current.temperature_k()
            || release.crystalline_mass_fraction != current.crystalline_mass_fraction()
        {
            return Err(
                invalid("mechanical release does not match the accepted thermal T/X").into(),
            );
        }
        let initial_volume = pair
            .volume_state(current.temperature_k(), current.crystalline_mass_fraction())
            .map_err(CrystallizationMechanicsError::from)?;
        let reference_volume_per_mass = cell.geometry.volume_m3 / cell.mass_kg;
        if !reference_volume_per_mass.is_finite()
            || reference_volume_per_mass <= 0.0
            || reference_volume_per_mass != self.state.reference_specific_volume_m3_kg()
        {
            return Err(invalid("mechanical and thermal fixed reference V0/m0 differ").into());
        }
        let removed_elastic_energy_j =
            finite(self.energy.removed_elastic_energy_density_j_m3 * cell.geometry.volume_m3)?;
        let thermalized_energy_j = finite(self.energy.thermalized_energy_j_kg * cell.mass_kg)?;
        let exported_energy_j = finite(self.energy.exported_energy_j_kg * cell.mass_kg)?;
        let partition_residual_j =
            finite(finite(thermalized_energy_j + exported_energy_j)? - removed_elastic_energy_j)?;
        checkpoint()?;
        let (next, thermal_audit) = thermal.apply_heat_impulses(
            &BTreeMap::from([(cell_id, thermalized_energy_j)]),
            controls,
            &mut checkpoint,
        )?;
        // The original pair's declared common domain remains a boundary even
        // when the same branch charts cover a wider kinetic-source interval.
        let final_state = next.states()[&cell_id];
        let final_volume = pair
            .volume_state(
                final_state.temperature_k(),
                final_state.crystalline_mass_fraction(),
            )
            .map_err(CrystallizationMechanicsError::from)?;
        let volume_roundoff = 32.0
            * f64::EPSILON
            * initial_volume
                .specific_volume_m3_kg
                .max(final_volume.specific_volume_m3_kg);
        if (final_volume.specific_volume_m3_kg - initial_volume.specific_volume_m3_kg).abs()
            > volume_roundoff
        {
            return Err(invalid(
                "density-changing unlocking requires a pressure/volume-work owner",
            )
            .into());
        }
        let balance_residual_j = finite(
            finite(thermal_audit.thermal.stored_energy_change_j + exported_energy_j)?
                - removed_elastic_energy_j,
        )?;
        // Bound conversion/partition rounding and the thermal owner's absolute
        // caloric/representation defects separately, without cancellation.
        let absolute_defect_j =
            finite(partition_residual_j.abs() + thermal_audit.absolute_energy_residual_j)?
                .max(balance_residual_j.abs());
        if absolute_defect_j > controls.cell_energy_tolerance_j
            || absolute_defect_j > controls.total_energy_tolerance_j
        {
            return Err(DepositionError::EnergyBalance(absolute_defect_j).into());
        }
        checkpoint()?;
        Ok((
            next,
            ThermalizedCrystallizationUnlock {
                cell_id,
                state: self.state,
                energy: self.energy,
                thermal_audit,
                removed_elastic_energy_j,
                thermalized_energy_j,
                exported_energy_j,
                balance_residual_j,
            },
        ))
    }
}

fn finite(value: f64) -> Result<f64, CrystallizationMechanicsError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CrystallizationMechanicsError::Invalid(
            "extensive unlocking energy is unrepresentable",
        ))
    }
}
