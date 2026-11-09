//! UV inversion relative to an accepted source-bound SP caloric snapshot.
//!
//! A changing-mass consumer can scale its original feed basis and volume, then
//! request a direct U increment from the accepted upstream/receiver state.
//! The existing TV/Gibbs and UV root owners do all physical closure. This does
//! not change feed elements, rebase the prepared path, or supply a flow law.

use super::*;
use crate::equilibrium::composition::isentrope::EquilibriumCaloricAnchorV1;

/// One UV endpoint and its direct balances against the caller's accepted anchor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumAnchoredEnergyVolumeSampleV1 {
    /// Same accepted native state, with the prepared path's own reference intact.
    pub sample: EquilibriumEnergyVolumeSampleV1,
    /// Requested U-U_anchor [J], on the original feed basis.
    pub target_energy_change_j: f64,
    /// Original source-integrated U-U_anchor [J].
    pub internal_energy_change_j: f64,
    /// Actual minus requested anchored internal energy [J].
    pub energy_residual_j: f64,
    /// Direct source and mixing S-S_anchor [J/K].
    pub entropy_change_j_k: f64,
}

impl EquilibriumEnergyVolumePathV1<'_> {
    /// Close volume and a source-relative energy increment from an accepted SP state.
    ///
    /// The anchor must retain exactly this original chemistry, element/charge
    /// inventory and feed normalization. Its birth entropy need not equal this
    /// path's reference entropy: receiving gas is generally irreversible.
    /// Neither reference is altered, and a failed/cancelled call invalidates the
    /// current publishable trial. No absolute formation-energy subtraction or
    /// unreacted-feed heat-release term enters the target.
    pub fn sample_from_caloric_anchor(
        &mut self,
        anchor: &EquilibriumCaloricAnchorV1,
        volume: Volume,
        energy_change: Energy,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumAnchoredEnergyVolumeSampleV1, EquilibriumVolumeErrorV1> {
        let previous = self.current.take();
        poll(&mut cancelled)?;
        if !volume.value().is_finite() || volume.value() <= 0.0
            || !energy_change.value().is_finite() {
            return Err(EquilibriumVolumeErrorV1::Invalid(
                "positive finite volume and finite anchored energy increment required",
            ));
        }
        let reference = anchor.reference_for(self.owner, &self.work,
            self.controls.minimum_temperature_k, self.controls.maximum_temperature_k)?;
        let preferred = previous.map(|sample| sample.temperature.value());
        let owner = self.owner;
        let (temperature, counts, iterations, volume_residual) =
            self.work.solve_internal_energy(
                owner, volume, energy_change.value(), self.controls, preferred,
                &mut cancelled,
                |work, t, cancel| Ok(reference.energy_change(owner, work, t, cancel)?),
            )?;
        // The ordinary sample still belongs to the immutable UV birth reference.
        // The caller's independent target is admitted below before publication.
        let sample = self.finish_sample(
            temperature, volume, None, counts, iterations, volume_residual,
            &mut cancelled,
        )?;
        self.current = None;
        let du = reference.energy_change(owner, &self.work, temperature, &mut cancelled)?;
        let ds = reference.entropy_change(owner, &self.work, temperature, &mut cancelled)?;
        let residual = du - energy_change.value();
        if !residual.is_finite() || residual.abs() > self.controls.energy_tolerance_j {
            return Err(EquilibriumVolumeErrorV1::InternalEnergyNoConvergence {
                iterations, residual_j: residual,
            });
        }
        poll(&mut cancelled)?;
        self.current = Some(sample);
        Ok(EquilibriumAnchoredEnergyVolumeSampleV1 {
            sample, target_energy_change_j: energy_change.value(),
            internal_energy_change_j: du, energy_residual_j: residual,
            entropy_change_j_k: ds,
        })
    }
}
