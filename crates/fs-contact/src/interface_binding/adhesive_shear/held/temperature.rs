//! Scalar junction continuation through the original thermal normal owner.
//!
//! The external normal owner supplies its sealed held-indentation or held-load
//! update and keeps all normal mechanical and parameter work. This binding
//! checks its accepted solid area/repulsion, refreshes the same selected shear
//! claims, and binds the successor's actual solid junction. It adds no heat,
//! tangential displacement, elastic shear store, or slip/peel transition.

use super::{AdhesiveShearError, SourceAdhesiveContactUpdate, SourceAdhesiveTangential};
use crate::normal_patch::{
    adhesive::{
        SourceJkrSphereContact,
        temperature::{
            SourceDmtThermalRefresh, SourceJkrThermalRefresh, SourceMaugisThermalRefresh,
        },
    },
    dmt::SourceDmtSphereContact,
    maugis::SourceMaugisSphereContact,
};

impl SourceAdhesiveTangential {
    /// Rebind the original scalar junction to an admitted JKR thermal successor.
    /// The receipt's predecessor must match the currently accepted real area
    /// and repulsion, as well as its original ordered normal/shear source state.
    ///
    /// The caller retains the original normal control and work account. For a
    /// held-load receipt, normal displacement and actuator work are external to
    /// this scalar tangential binding. For held indentation, normal temperature
    /// can still change the solid area and both junction capacities.
    pub fn rebind_jkr_temperature(
        &self,
        original: &SourceJkrSphereContact,
        refreshed: &SourceJkrThermalRefresh,
    ) -> Result<Self, AdhesiveShearError> {
        self.check_thermal_predecessor(SourceAdhesiveContactUpdate::from_jkr(refreshed.before())?)?;
        let shear = self.shear.rebind_jkr_temperature(original, refreshed)?;
        shear.bind_jkr(refreshed.source().normal_state(), refreshed.response())
    }

    /// Rebind the original DMT solid junction through its sealed thermal change.
    /// The finite-range attractive gap has zero solid area before and after a
    /// held-indentation refresh, so it cannot acquire friction or heat. Actual
    /// held-load normal motion remains with the external normal owner.
    pub fn rebind_dmt_temperature(
        &self,
        original: &SourceDmtSphereContact,
        refreshed: &SourceDmtThermalRefresh,
    ) -> Result<Self, AdhesiveShearError> {
        self.check_thermal_predecessor(SourceAdhesiveContactUpdate::from_dmt(refreshed.before())?)?;
        let shear = self.shear.rebind_dmt_temperature(original, refreshed)?;
        shear.bind_dmt(refreshed.source().normal_state(), refreshed.response())
    }

    /// Rebind the full original Maugis equilibrium, including its noncontact or
    /// detached phase. The cohesive annulus never becomes a shear-bearing disk.
    /// The normal owner retains its branch/control and absolute-potential work;
    /// only the actual successor solid junction enters this scalar consumer.
    pub fn rebind_maugis_temperature(
        &self,
        original: &SourceMaugisSphereContact,
        refreshed: &SourceMaugisThermalRefresh,
    ) -> Result<Self, AdhesiveShearError> {
        self.check_thermal_predecessor(SourceAdhesiveContactUpdate::from_maugis_equilibrium(
            original,
            refreshed.before(),
        )?)?;
        let shear = self.shear.rebind_maugis_temperature(original, refreshed)?;
        shear.bind_maugis_equilibrium(refreshed.source(), refreshed.response())
    }

    fn check_thermal_predecessor(
        &self,
        before: SourceAdhesiveContactUpdate,
    ) -> Result<(), AdhesiveShearError> {
        // Use the same identity gate as ordinary normal-geometry updates. An
        // unchanged material identity alone cannot admit an older patch state.
        let previous = self.tangential_for_update(before)?;
        if previous.contact() != self.tangential.contact() {
            return Err(AdhesiveShearError::Invalid {
                reason: "thermal junction refresh must match its accepted solid area and repulsion",
            });
        }
        Ok(())
    }
}
