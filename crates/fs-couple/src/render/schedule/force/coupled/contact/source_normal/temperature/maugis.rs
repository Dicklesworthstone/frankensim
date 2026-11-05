//! Original finite-zone temperature continuation on the complete native network.

use super::super::{
    BoundNormalContact, CancelGate, SourceNormalModalError, SourceNormalModalSystem, finite, limit,
    poll, solve,
};
use super::prepare_maugis_contact_refresh;
use fs_contact::normal_patch::adhesive::temperature::SourceMaugisThermalRefresh;

/// Accepted finite-zone source change at unchanged native sample and geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceNormalMaugisRefreshReport {
    pub sample: u64,
    pub approach_m: f64,
    pub energy_before_j: f64,
    pub energy_after_j: f64,
    pub contact_energy_before_j: f64,
    pub contact_energy_after_j: f64,
    /// Signed work supplied by the prescribed temperature driver, not heat.
    pub parameter_work_j: f64,
    /// Actual contact-potential change minus original signed parameter work.
    pub energy_residual_j: f64,
    /// Contact-sized tolerance; native stored energy does not enlarge it.
    pub energy_tolerance_j: f64,
}

impl SourceNormalModalSystem {
    /// Apply a sealed original Maugis refresh at actual accepted indentation.
    /// Retains the selected contact/noncontact policy and original work budgets.
    /// Native bodies, coefficients, Maxwell history, losses and clock do not
    /// change; the explicit temperature driver owns the returned signed work.
    pub fn apply_maugis_source_refresh(
        &mut self,
        refreshed: &SourceMaugisThermalRefresh,
        gate: &CancelGate,
    ) -> Result<SourceNormalMaugisRefreshReport, SourceNormalModalError> {
        self.apply_maugis_source_refresh_with_admission(refreshed, gate, |report, _| Ok(*report))
    }

    /// Stage a receiving history before publishing the new original contact.
    /// The callback returns a candidate whose later installation is infallible.
    /// Any refusal or late cancellation preserves all accepted native state.
    pub fn apply_maugis_source_refresh_with_admission<R, E>(
        &mut self,
        refreshed: &SourceMaugisThermalRefresh,
        gate: &CancelGate,
        admit: impl FnOnce(&SourceNormalMaugisRefreshReport, &BoundNormalContact) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<SourceNormalModalError>,
    {
        let (candidate, report) = self
            .prepare_maugis_refresh(refreshed, gate)
            .map_err(E::from)?;
        let receiving = admit(&report, &candidate)?;
        poll(Some(gate))
            .map_err(SourceNormalModalError::from)
            .map_err(E::from)?;
        self.contact = candidate;
        Ok(receiving)
    }

    fn prepare_maugis_refresh(
        &self,
        refreshed: &SourceMaugisThermalRefresh,
        gate: &CancelGate,
    ) -> Result<(BoundNormalContact, SourceNormalMaugisRefreshReport), SourceNormalModalError> {
        poll(Some(gate))?;
        let approach = self.approach(self.components())?;
        self.check_approach(approach)?;
        let (candidate, change) =
            prepare_maugis_contact_refresh(&self.contact, self.config, approach, refreshed, gate)?;
        let native_energy = self.network.total_energy_j()?;
        let energy_before = finite(native_energy + change.contact_energy_before_j)?;
        let energy_after = finite(native_energy + change.contact_energy_after_j)?;
        limit(
            "refreshed native/source energy magnitude",
            energy_after.abs(),
            self.network.config.maximum_total_energy_j,
        )?;

        solve::positive_slope(
            candidate.source(),
            approach,
            approach,
            self.compliance_m_per_n,
        )?;
        let report = SourceNormalMaugisRefreshReport {
            sample: self.samples_rendered(),
            approach_m: approach,
            energy_before_j: energy_before,
            energy_after_j: energy_after,
            contact_energy_before_j: change.contact_energy_before_j,
            contact_energy_after_j: change.contact_energy_after_j,
            parameter_work_j: change.parameter_work_j,
            energy_residual_j: change.energy_residual_j,
            energy_tolerance_j: change.energy_tolerance_j,
        };
        poll(Some(gate))?;
        Ok((candidate, report))
    }
}
