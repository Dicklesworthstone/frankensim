//! Released elastic heat in the original constant-pressure material parcel.
//!
//! One uniform mechanical representative covers the parcel's entire immutable
//! reference volume and mass. The original pair supplies both caloric branches
//! and the actual free volume before and after heating. The pressure reservoir
//! receives the signed work p*(V_after-V_before); the parcel retains delta_U =
//! delta_H - p*delta_V. Its latent heat remains inside the original enthalpy.
//!
//! The initial elastic volume strain must be arithmetic zero: this handoff
//! cannot erase a compressed or dilated mechanical volume without an owner for
//! its relaxation work. The admitted small-strain representative therefore has
//! the original free volume at release, within its infinitesimal kinematics.
//! Finite-strain volume corrections require a separate mechanical owner.
//!
//! This is a quasistatic constant-pressure endpoint transition. It advances no
//! kinetic clock and supplies no inertia, pressure evolution, moving mesh,
//! momentum, viscosity or contact update. Actual parcel volume and fixed
//! mechanical reference volume remain distinct. Source-pressure applicability
//! is the original parcel caller's explicit material admission.
//!
//! The authentic release is consumed. Adopt the returned parcel and unlocked
//! state together; a refusal or cancellation publishes neither candidate.

use super::{CrystallizationUnlock, CrystallizationUnlockEnergy, UnlockedCrystallization};
use crate::phase::crystallization::CrystallizationCaloricPair;
use crate::phase::crystallization::mechanics::CrystallizationMechanicsError;
use fs_thermal::crystallization::CrystallizationError;
use fs_thermal::crystallization::pressure::parcel::IsobaricCrystallizationParcel;
use fs_thermal::crystallization::pressure::{
    IsobaricCrystallizationControls, IsobaricCrystallizationStep,
};

/// Original material admission and caloric/pressure refusals retain their causes.
#[derive(Clone, Debug, PartialEq)]
pub enum CrystallizationUnlockIsobaricError {
    Mechanics(CrystallizationMechanicsError),
    Thermal(CrystallizationError),
}

impl From<CrystallizationMechanicsError> for CrystallizationUnlockIsobaricError {
    fn from(error: CrystallizationMechanicsError) -> Self {
        Self::Mechanics(error)
    }
}
impl From<CrystallizationError> for CrystallizationUnlockIsobaricError {
    fn from(error: CrystallizationError) -> Self {
        Self::Thermal(error)
    }
}
impl core::fmt::Display for CrystallizationUnlockIsobaricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "crystallization unlock isobaric transfer: {self:?}")
    }
}
impl std::error::Error for CrystallizationUnlockIsobaricError {}

/// Complete elastic, internal-energy and constant-pressure work account.
/// This report cannot reconstruct or thermalize the consumed release.
#[derive(Debug)]
pub struct ThermalizedIsobaricCrystallizationUnlock {
    pub state: UnlockedCrystallization,
    /// Original reference storage debit and heat/export partition [J/kg].
    pub energy: CrystallizationUnlockEnergy,
    /// The parcel's actual enthalpy, internal-energy and pressure-work changes.
    /// Its supplied heat is an INTERNAL transfer from the released network.
    pub thermal: IsobaricCrystallizationStep,
    /// Original elastic energy density times fixed mechanical reference volume.
    pub removed_elastic_energy_j: f64,
    /// Original heat per mass times the parcel's conserved mass.
    pub thermalized_energy_j: f64,
    /// Credit owed to the caller's separately selected nonthermal destination.
    /// The pressure reservoir's work is separately reported in `thermal`.
    pub exported_energy_j: f64,
    /// delta_U + actual p*delta_V + export - removed elastic storage [J].
    pub balance_residual_j: f64,
}

impl CrystallizationUnlock {
    /// Transfer this release once to its source-bound isobaric parcel.
    ///
    /// Both original caloric identities, accepted release T/X and immutable
    /// reference V0/m0 must match. The original elastic normal-strain trace
    /// must have magnitude at most 8*EPSILON*sum(abs(normal elastic strains)).
    /// This is a summation-roundoff bound, not a compressibility allowance;
    /// unresolved elastic volume relaxation refuses.
    ///
    /// The initial actual parcel volume is the original pair's free volume
    /// times conserved mass. The same source must also admit the final state
    /// and volume. These source-volume comparisons allow at most 32 machine
    /// epsilons relative rounding, independently of the caller's solve budget.
    ///
    /// The existing parcel adds only the released heat to its original h
    /// coordinate. Pressure work is not subtracted from that heat a second
    /// time: delta_U = delta_H - p*delta_V closes the combined account. Any
    /// selected equilibrium remelting belongs to the original caloric owner.
    /// Neither the parcel clock nor the fixed reference volume changes.
    ///
    /// The caller must adopt both returned states in the same accepted-state
    /// transaction and credit any explicit nonthermal export to its own owner.
    pub fn thermalize_isobaric<'law, 'source>(
        self,
        pair: &CrystallizationCaloricPair<'_>,
        parcel: &IsobaricCrystallizationParcel<'law, 'source>,
        controls: IsobaricCrystallizationControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        (
            IsobaricCrystallizationParcel<'law, 'source>,
            ThermalizedIsobaricCrystallizationUnlock,
        ),
        CrystallizationUnlockIsobaricError,
    > {
        checkpoint(&mut cancelled)?;
        controls.validate()?;
        let invalid = CrystallizationMechanicsError::Invalid;
        let current = parcel.state();
        let law = parcel.law();
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
                invalid("mechanical release does not match the accepted isobaric T/X").into(),
            );
        }
        let reference_volume_per_mass = parcel.reference_volume_m3() / parcel.mass_kg();
        if !reference_volume_per_mass.is_finite()
            || reference_volume_per_mass <= 0.0
            || reference_volume_per_mass != self.state.reference_specific_volume_m3_kg()
        {
            return Err(
                invalid("mechanical and isobaric parcel fixed reference V0/m0 differ").into(),
            );
        }
        if self.state.release_elastic_volumetric_strain.abs()
            > self.state.release_elastic_trace_roundoff
        {
            return Err(invalid(
                "isobaric unlocking requires an owner for elastic volume relaxation",
            )
            .into());
        }
        let initial_volume = pair
            .volume_state(current.temperature_k(), current.crystalline_mass_fraction())
            .map_err(CrystallizationMechanicsError::from)?;
        admit_source_volume(
            parcel.volume_m3(),
            extensive(initial_volume.specific_volume_m3_kg, parcel.mass_kg())?,
        )?;
        let removed_elastic_energy_j = extensive(
            self.energy.removed_elastic_energy_density_j_m3,
            parcel.reference_volume_m3(),
        )?;
        let thermalized_energy_j =
            extensive(self.energy.thermalized_energy_j_kg, parcel.mass_kg())?;
        let exported_energy_j = extensive(self.energy.exported_energy_j_kg, parcel.mass_kg())?;
        let partition_residual_j =
            finite(finite(thermalized_energy_j + exported_energy_j)? - removed_elastic_energy_j)?;
        checkpoint(&mut cancelled)?;
        let (next, thermal) = parcel.apply_heat(thermalized_energy_j, controls, &mut cancelled)?;
        // The original pair's common domain remains authoritative even when
        // the individual kinetic branches have wider temperature support.
        let final_state = next.state();
        let final_volume = pair
            .volume_state(
                final_state.temperature_k(),
                final_state.crystalline_mass_fraction(),
            )
            .map_err(CrystallizationMechanicsError::from)?;
        admit_source_volume(
            next.volume_m3(),
            extensive(final_volume.specific_volume_m3_kg, next.mass_kg())?,
        )?;
        let balance_residual_j = finite(
            finite(
                finite(thermal.internal_energy_change_j + thermal.boundary_work_j)?
                    + exported_energy_j,
            )? - removed_elastic_energy_j,
        )?;
        // Do not allow opposite caloric, pressure-work or partition defects to
        // cancel each other in the combined signed balance.
        let absolute_defect_j =
            finite(partition_residual_j.abs() + thermal.absolute_energy_residual_j)?
                .max(balance_residual_j.abs());
        if absolute_defect_j > controls.energy_tolerance_j {
            return Err(CrystallizationError::Invalid(
                "isobaric unlocking energy residual exceeds its tolerance",
            )
            .into());
        }
        checkpoint(&mut cancelled)?;
        Ok((
            next,
            ThermalizedIsobaricCrystallizationUnlock {
                state: self.state,
                energy: self.energy,
                thermal,
                removed_elastic_energy_j,
                thermalized_energy_j,
                exported_energy_j,
                balance_residual_j,
            },
        ))
    }
}

fn checkpoint(
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), CrystallizationUnlockIsobaricError> {
    if cancelled() {
        Err(CrystallizationError::Cancelled.into())
    } else {
        Ok(())
    }
}

fn finite(value: f64) -> Result<f64, CrystallizationMechanicsError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CrystallizationMechanicsError::Invalid(
            "extensive isobaric unlocking account is unrepresentable",
        ))
    }
}

fn extensive(value: f64, scale: f64) -> Result<f64, CrystallizationMechanicsError> {
    let product = finite(value * scale)?;
    if value != 0.0 && scale != 0.0 && product == 0.0 {
        return Err(CrystallizationMechanicsError::Invalid(
            "nonzero isobaric unlocking account underflows",
        ));
    }
    Ok(product)
}

fn admit_source_volume(
    actual_m3: f64,
    source_m3: f64,
) -> Result<(), CrystallizationMechanicsError> {
    if !actual_m3.is_finite()
        || actual_m3 <= 0.0
        || !source_m3.is_finite()
        || source_m3 <= 0.0
        || (actual_m3 - source_m3).abs() > 32.0 * f64::EPSILON * actual_m3.max(source_m3)
    {
        return Err(CrystallizationMechanicsError::Invalid(
            "isobaric parcel volume differs from its original caloric pair",
        ));
    }
    Ok(())
}
