//! Volume queries on the existing prepared reacting isentrope.
//!
//! Density inversion retains its original workspace and source-relative
//! enthalpy/entropy increments. Only the accepted endpoint is materialized.
use super::super::IdealGasIsochoricEquilibriumV1;
use crate::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use crate::equilibrium::composition::{
    CompositionEquilibriumErrorV1, EquilibriumExpansionErrorV1,
    isentrope::{
        EquilibriumIsentropeV1,
        density::{
            EquilibriumDensityErrorV1, EquilibriumDensityQueryV1, EquilibriumDensitySolveV1,
        },
    },
};
use fs_qty::{Energy, Volume};

/// Accepted fixed-volume endpoint and the work of its prepared density query.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumVolumeQueryV1 {
    state: IdealGasIsochoricEquilibriumV1,
    query: EquilibriumDensityQueryV1,
}
impl EquilibriumVolumeQueryV1 {
    /// Complete original source-bearing equilibrium at the requested volume.
    pub const fn state(&self) -> &IdealGasIsochoricEquilibriumV1 {
        &self.state
    }
    /// Source-relative increments and accumulated pressure/composition work.
    pub const fn query(&self) -> EquilibriumDensityQueryV1 {
        self.query
    }
    /// Transfer the accepted endpoint without another chemistry solve.
    pub fn into_parts(self) -> (IdealGasIsochoricEquilibriumV1, EquilibriumDensityQueryV1) {
        (self.state, self.query)
    }
}

impl EquilibriumIsentropeV1<'_> {
    /// Change volume on the SAME reference entropy and elemental inventory.
    ///
    /// The conserved feed mass gives the density target m/V. In addition to
    /// density closure, the actual ideal-gas volume NRT/p must meet the same
    /// relative tolerance. Composition mass error therefore cannot hide a
    /// volume defect; callers should admit chemistry more tightly than volume.
    /// The source-relative enthalpy increment remains available for accurate
    /// mechanical work when absolute formation energies are large.
    ///
    /// This changes numerical warm starts, not a physical accepted parcel.
    /// Only a returned result admits a volume endpoint. Every nested query
    /// retains the original temperature, pressure, iteration and source limits.
    pub fn sample_volume(
        &mut self,
        volume: Volume,
        controls: EquilibriumDensitySolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumVolumeQueryV1, EquilibriumDensityErrorV1> {
        // Invalid/zero/nonfinite volumes produce invalid density targets and
        // pass through the original query's cancellation/invalidation owner.
        let target = self.feed_mass_kg() / volume.value();
        let query = self.sample_density(target, controls, &mut cancelled)?;
        let thermodynamic = self.publish_current(&mut cancelled)?;
        let chemical = thermodynamic.equilibrium();
        let nrt = chemical.total_amount().value() * R * chemical.temperature().value();
        let actual_volume = nrt / chemical.pressure().value();
        let relative_volume_residual = (actual_volume - volume.value()) / volume.value();
        let energy = chemical.enthalpy().value() - nrt;
        if !volume.value().is_finite()
            || volume.value() <= 0.0
            || !energy.is_finite()
            || !actual_volume.is_finite()
            || !relative_volume_residual.is_finite()
            || relative_volume_residual.abs() > controls.relative_density_tolerance
        {
            return Err(EquilibriumDensityErrorV1::Invalid(
                "prepared equilibrium failed original volume or internal-energy admission",
            ));
        }
        if cancelled() {
            return Err(EquilibriumExpansionErrorV1::from(
                CompositionEquilibriumErrorV1::Cancelled,
            )
            .into());
        }
        Ok(EquilibriumVolumeQueryV1 {
            state: IdealGasIsochoricEquilibriumV1 {
                thermodynamic,
                volume,
                internal_energy: Energy::new(energy),
                target_internal_energy: None,
                relative_volume_residual,
                pressure_iterations: query.pressure_queries,
            },
            query,
        })
    }
}
