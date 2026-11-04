//! Fixed-strain release of an explicitly unlocked crystalline elastic network.
//!
//! The caller selects the loss-of-support event. No percolation threshold,
//! stiffness/fraction interpolation or melting kinetics is inferred here. The
//! original full elastic tensor owns recoverable energy; an arbitrary
//! `SmallStrainLaw` declaration is insufficient to infer that energy.
//!
//! Unlocking removes this elastic contribution at unchanged total strain and
//! unchanged caloric state. Its lost energy is partitioned into internal heat
//! and an explicitly exported nonthermal transfer. The thermal owner may add
//! ONLY that heat to its original enthalpy once. Latent heat already belongs to
//! the original phase law. Changes of eigenstrain or coefficients BEFORE this
//! event retain their separate parameter-work account.
//!
//! This is an instantaneous constitutive transition at fixed displacement, not
//! a fluid, momentum, viscosity, pressure or contact-topology solver. A native
//! assembly may correctly refuse the resulting unconstrained zero-stiffness
//! mechanism. Subsequent locking requires another explicit stress-free choice.

use super::{
    CrystallizationCaloricPair, CrystallizationLock, CrystallizationMechanicsError,
    CrystallizationStrainLimits, CrystallizationVolumeState, FrozenCrystallizationElastic, within,
};
use crate::elastic::acoustic::ElasticAcousticState;
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_blake3::ContentHash;
use fs_evidence::{Ambition, ModelCard, ValidityDomain};

/// Atomic transfer of a released network's heat to its original thermal cell.
pub mod thermal;

/// A debit from recoverable elastic storage. Positive heat is deposited into
/// the material; positive export must be credited to the caller's other owner.
/// Neither term is external heating of a boundary enclosing both participants.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrystallizationUnlockEnergy {
    /// Original energy per fixed mechanical reference volume [J/m^3].
    pub removed_elastic_energy_density_j_m3: f64,
    /// Original storage debit, expressed as a positive amount [J/kg].
    pub removed_elastic_energy_j_kg: f64,
    pub thermalized_energy_j_kg: f64,
    pub exported_energy_j_kg: f64,
    /// `-removed + thermalized + exported`, in J/kg.
    pub balance_residual_j_kg: f64,
}

/// One proposal containing the released elastic state and its complete debit.
/// Publish both together with the destination credits in the caller's existing
/// accepted-state transaction. No thermal cell or mechanical field is mutated.
#[derive(Debug)]
pub struct CrystallizationUnlock {
    state: UnlockedCrystallization,
    energy: CrystallizationUnlockEnergy,
}

impl CrystallizationUnlock {
    pub const fn state(&self) -> &UnlockedCrystallization {
        &self.state
    }
    pub const fn energy(&self) -> CrystallizationUnlockEnergy {
        self.energy
    }

    /// Consume the authentic release for an explicitly caller-owned transfer.
    /// Neither returned component can reconstruct another consumable release.
    /// Use [`Self::thermalize`] for an atomic handoff to the original timeline.
    pub fn into_parts(self) -> (UnlockedCrystallization, CrystallizationUnlockEnergy) {
        (self.state, self.energy)
    }
}

/// Zero recoverable elastic contribution after explicit network release.
///
/// This type cannot release energy again. The original frozen snapshot is
/// consumed by unlocking; reconstructing it to repeat the debit is a caller
/// history error. This state supplies no liquid compressibility or viscous
/// stress. Those remain separate physical owners if that regime is needed.
#[derive(Debug)]
pub struct UnlockedCrystallization {
    previous_lock: CrystallizationLock,
    release_volume: CrystallizationVolumeState,
    release_total_strain: Voigt,
    elastic_state_identity: ContentHash,
    reference_specific_volume_m3_kg: f64,
    limits: CrystallizationStrainLimits,
}

impl FrozenCrystallizationElastic<'_, ElasticAcousticState> {
    /// Release the original tensor's entire stored elastic energy at fixed
    /// strain, converting it through the caller's immutable reference `V0/m0`.
    ///
    /// `reference_specific_volume_m3_kg` is the SAME mechanical quadrature
    /// reference volume divided by the SAME transported reference mass. It is
    /// not the evolving free mixture volume or reciprocal current density.
    /// This explicit convention supports the original thermal `mass * h`
    /// account without inventing mass or pressure work during unlocking.
    ///
    /// `heat_fraction` in [0,1] is an explicit dissipation/transfer assumption.
    /// Its complement is exported energy requiring a named caller-owned
    /// destination; it is not silently discarded. A value of one thermalizes
    /// the entire release. Caloric evolution and event selection remain with
    /// their existing owners. All checks precede publishing the proposal.
    pub fn unlock_at_fixed_strain(
        self,
        total_strain: Voigt,
        reference_specific_volume_m3_kg: f64,
        heat_fraction: f64,
    ) -> Result<CrystallizationUnlock, CrystallizationMechanicsError> {
        let invalid = CrystallizationMechanicsError::Invalid;
        if !reference_specific_volume_m3_kg.is_finite()
            || reference_specific_volume_m3_kg <= 0.0
            || !heat_fraction.is_finite()
            || !(0.0..=1.0).contains(&heat_fraction)
        {
            return Err(invalid(
                "positive fixed reference V0/m0 and heat fraction in [0,1] required",
            ));
        }
        let strain = self.elastic_strain(&total_strain)?;
        let response = self
            .elastic
            .strain_response(crate::tensor::to_mandel(&strain))
            .map_err(|_| {
                invalid("original elastic source refuses release strain or stored energy")
            })?;
        let density = response.stored_energy_density_j_m3;
        let removed = density * reference_specific_volume_m3_kg;
        if !removed.is_finite() || removed < 0.0 || (density > 0.0 && removed == 0.0) {
            return Err(invalid(
                "released specific elastic energy is unrepresentable",
            ));
        }
        let thermalized = heat_fraction * removed;
        // Use a remainder so both destination credits account for one debit.
        let exported = removed - thermalized;
        let residual = (thermalized + exported) - removed;
        if !thermalized.is_finite() || !exported.is_finite() || !residual.is_finite() {
            return Err(invalid("unlocking energy partition is unrepresentable"));
        }
        Ok(CrystallizationUnlock {
            state: UnlockedCrystallization {
                previous_lock: self.lock,
                release_volume: self.volume,
                release_total_strain: total_strain,
                elastic_state_identity: response.material_state_identity,
                reference_specific_volume_m3_kg,
                limits: self.limits,
            },
            energy: CrystallizationUnlockEnergy {
                removed_elastic_energy_density_j_m3: density,
                removed_elastic_energy_j_kg: removed,
                thermalized_energy_j_kg: thermalized,
                exported_energy_j_kg: exported,
                balance_residual_j_kg: residual,
            },
        })
    }
}

impl UnlockedCrystallization {
    pub const fn previous_lock(&self) -> &CrystallizationLock {
        &self.previous_lock
    }
    /// Caloric/density state AT release, not a prediction of later cooling.
    pub const fn release_volume_state(&self) -> &CrystallizationVolumeState {
        &self.release_volume
    }
    pub const fn release_total_strain(&self) -> &Voigt {
        &self.release_total_strain
    }
    pub const fn elastic_state_identity(&self) -> ContentHash {
        self.elastic_state_identity
    }
    pub const fn reference_specific_volume_m3_kg(&self) -> f64 {
        self.reference_specific_volume_m3_kg
    }

    /// Propose a new caller-selected stress-free lock using the SAME original
    /// caloric pair and a separately admitted current strain and phase state.
    /// No threshold or crystallization history is reconstructed from stiffness.
    /// The unlocked checkpoint remains unchanged if preparation refuses.
    ///
    /// The new free-volume reference describes subsequent eigenstrain. It
    /// does NOT change the caller's fixed quadrature `V0/m0` or transported mass.
    /// Preparing a new elastic snapshot must independently admit its original
    /// elastic source at that temperature/fraction. Locking itself stores zero
    /// elastic energy at the selected total strain and creates no heat credit.
    pub fn relock_at(
        &self,
        pair: &CrystallizationCaloricPair<'_>,
        temperature_k: f64,
        crystalline_mass_fraction: f64,
        total_strain: Voigt,
    ) -> Result<CrystallizationLock, CrystallizationMechanicsError> {
        if pair.identity() != self.release_volume.caloric_pair_identity {
            return Err(CrystallizationMechanicsError::PairMismatch);
        }
        if !within(
            &total_strain,
            self.limits.maximum_absolute_total_mandel_strain,
        ) {
            return Err(CrystallizationMechanicsError::Invalid(
                "new locking strain leaves the declared total-strain domain",
            ));
        }
        CrystallizationLock::new(pair, temperature_k, crystalline_mass_fraction, total_strain)
    }
}

impl SmallStrainLaw for UnlockedCrystallization {
    type State = ();
    fn initial_state(&self) {}
    fn stress(&self, strain: &Voigt, (): &()) -> Voigt {
        if within(strain, self.limits.maximum_absolute_total_mandel_strain) {
            [0.0; 6]
        } else {
            [f64::NAN; 6]
        }
    }
    fn tangent(&self, strain: &Voigt, (): &()) -> Tangent6 {
        if within(strain, self.limits.maximum_absolute_total_mandel_strain) {
            [[0.0; 6]; 6]
        } else {
            [[f64::NAN; 6]; 6]
        }
    }
    fn update_state(&self, _strain: &Voigt, (): &()) {}
    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true, // identically zero after the release
            dissipation_nonnegative: true,
            polyconvex: Some(true),
            tangent_symmetric: true,
            failure_envelope: "released elastic contribution; declared total-strain bounds; explicit relock only",
        }
    }
    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.unlocked-crystallization-elastic",
            "0.1.0",
            Ambition::Frontier,
            vec![
                "caller-selected complete network release at fixed total strain".into(),
                "zero elastic storage; no automatic recovery of load-bearing support".into(),
                format!("released original elastic state {:?}; caloric pair {:?}",
                    self.elastic_state_identity, self.release_volume.caloric_pair_identity),
            ],
            ValidityDomain::unconstrained().with(
                "strain-magnitude",
                0.0,
                core::f64::consts::FRAC_1_SQRT_2 * self.limits.maximum_absolute_total_mandel_strain,
            ),
            vec![
                "no universal melting criterion, liquid pressure/viscosity, momentum or contact solve".into(),
                "energy destinations and thermal feedback require the caller's accepted-state transaction".into(),
            ],
            f64::INFINITY,
        )
    }
}

#[cfg(test)]
mod tests;
