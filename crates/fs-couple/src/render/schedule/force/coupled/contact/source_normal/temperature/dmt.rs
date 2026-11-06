//! Original finite-range DMT temperature continuation at held native motion.
//!
//! Contact, attractive gap, separated state and both piecewise corners retain
//! their original source semantics. The update supplies signed parameter work
//! at fixed geometry; it supplies no impulse or inferred caloric heat.

use super::super::{
    BoundNormalContact, CancelGate, SourceNormalModalError, SourceNormalModalSystem, finite, limit,
    poll, solve,
};
use super::prepare_dmt_contact_refresh;
use fs_contact::normal_patch::adhesive::temperature::SourceDmtThermalRefresh;

/// Accepted DMT source change, with the complete native state held fixed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceNormalDmtRefreshReport {
    pub sample: u64,
    pub approach_m: f64,
    pub energy_before_j: f64,
    pub energy_after_j: f64,
    pub contact_energy_before_j: f64,
    pub contact_energy_after_j: f64,
    /// Signed work of the explicit temperature driver, not frictional heat.
    pub parameter_work_j: f64,
    pub energy_residual_j: f64,
    /// Original contact-sized tolerance, independent of unrelated native energy.
    pub energy_tolerance_j: f64,
}

impl SourceNormalModalSystem {
    /// Continue an original sealed DMT source at actual accepted signed approach.
    /// The attraction range, geometry, work column and native histories remain
    /// fixed. The returned signed work belongs to the prescribed thermal driver.
    pub fn apply_dmt_source_refresh(
        &mut self,
        refreshed: &SourceDmtThermalRefresh,
        gate: &CancelGate,
    ) -> Result<SourceNormalDmtRefreshReport, SourceNormalModalError> {
        self.apply_dmt_source_refresh_with_admission(refreshed, gate, |report, _| Ok(*report))
    }

    /// Stage downstream state before publishing the complete source replacement.
    /// The callback returns an infallibly installable candidate. A receiver
    /// refusal or late cancellation leaves all accepted state unchanged.
    pub fn apply_dmt_source_refresh_with_admission<R, E>(
        &mut self,
        refreshed: &SourceDmtThermalRefresh,
        gate: &CancelGate,
        admit: impl FnOnce(&SourceNormalDmtRefreshReport, &BoundNormalContact) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<SourceNormalModalError>,
    {
        let (candidate, report) = self.prepare_dmt_refresh(refreshed, gate).map_err(E::from)?;
        let receiving = admit(&report, &candidate)?;
        poll(Some(gate))
            .map_err(SourceNormalModalError::from)
            .map_err(E::from)?;
        self.contact = candidate;
        Ok(receiving)
    }

    fn prepare_dmt_refresh(
        &self,
        refreshed: &SourceDmtThermalRefresh,
        gate: &CancelGate,
    ) -> Result<(BoundNormalContact, SourceNormalDmtRefreshReport), SourceNormalModalError> {
        poll(Some(gate))?;
        let approach = self.approach(self.components())?;
        self.check_approach(approach)?;
        let (candidate, change) =
            prepare_dmt_contact_refresh(&self.contact, self.config, approach, refreshed, gate)?;
        // A thermal change can strengthen the negative tail while the held
        // point is compressive or exactly at a corner. Preserve the original
        // mechanical envelope screen rather than checking only that point.
        solve::check_dmt_tail_stability(candidate.source(), self.compliance_m_per_n, self.config)?;
        solve::positive_slope(
            candidate.source(),
            approach,
            approach,
            self.compliance_m_per_n,
        )?;
        let native = self.network.total_energy_j()?;
        let before = finite(native + change.contact_energy_before_j)?;
        let after = finite(native + change.contact_energy_after_j)?;
        limit(
            "refreshed native/source energy magnitude",
            after.abs(),
            self.network.config.maximum_total_energy_j,
        )?;
        let report = SourceNormalDmtRefreshReport {
            sample: self.samples_rendered(),
            approach_m: approach,
            energy_before_j: before,
            energy_after_j: after,
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
