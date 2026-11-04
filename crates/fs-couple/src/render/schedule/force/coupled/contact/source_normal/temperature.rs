//! Prescribed original-source temperature continuation at fixed native state.
//!
//! Only the local JKR contact source changes. Native modal/connection/material
//! coefficients, Maxwell history, bodies and clock remain as supplied. Signed
//! contact-potential change is parameter work supplied by the explicit driver;
//! it is not inferred caloric heat, thermal expansion or autonomous feedback.

use super::{
    BoundNormalContact, CancelGate, NormalLawSource, NormalLawSourceError, SourceNormalModalError,
    SourceNormalModalSystem, check_contact_work, finite, invalid, limit, poll, solve,
};
use fs_contact::normal_patch::adhesive::temperature::{
    AdhesiveThermalControl, SourceJkrThermalRefresh,
};
use fs_math::dd::Dd;

/// One accepted source change, at unchanged native sample and coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceNormalJkrRefreshReport {
    pub sample: u64,
    pub approach_m: f64,
    pub energy_before_j: f64,
    pub energy_after_j: f64,
    pub contact_energy_before_j: f64,
    pub contact_energy_after_j: f64,
    pub parameter_work_j: f64,
    /// Contact-potential change minus the original signed parameter work.
    pub energy_residual_j: f64,
    /// Contact-sized tolerance; native stored energy does not enlarge it.
    pub energy_tolerance_j: f64,
}

impl SourceNormalModalSystem {
    /// Apply a sealed original JKR refresh at the actual accepted indentation.
    /// Old-source identity and canonical response prevent stale/replayed work.
    /// The returned signed work must be accounted by the prescribed driver.
    pub fn apply_jkr_source_refresh(
        &mut self,
        refreshed: &SourceJkrThermalRefresh,
        gate: &CancelGate,
    ) -> Result<SourceNormalJkrRefreshReport, SourceNormalModalError> {
        self.apply_jkr_source_refresh_with_admission(refreshed, gate, |report, _| Ok(*report))
    }

    /// Stage a receiving history against the new bound contact before either
    /// owner publishes. The callback returns a candidate, not an accepted-state
    /// mutation. Its installation after success must be infallible. Refusal or
    /// late cancellation preserves the old source and every native state.
    pub fn apply_jkr_source_refresh_with_admission<R, E>(
        &mut self,
        refreshed: &SourceJkrThermalRefresh,
        gate: &CancelGate,
        admit: impl FnOnce(&SourceNormalJkrRefreshReport, &BoundNormalContact) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<SourceNormalModalError>,
    {
        let (candidate, report) = self.prepare_jkr_refresh(refreshed, gate).map_err(E::from)?;
        let receiving = admit(&report, &candidate)?;
        poll(Some(gate))
            .map_err(SourceNormalModalError::from)
            .map_err(E::from)?;
        self.contact = candidate;
        Ok(receiving)
    }

    fn prepare_jkr_refresh(
        &self,
        refreshed: &SourceJkrThermalRefresh,
        gate: &CancelGate,
    ) -> Result<(BoundNormalContact, SourceNormalJkrRefreshReport), SourceNormalModalError> {
        poll(Some(gate))?;
        let NormalLawSource::JkrAttached(previous) = self.contact.source() else {
            return Err(invalid(
                "JKR temperature continuation requires the original attached JKR source",
            )
            .into());
        };
        let approach = self.approach(self.components())?;
        self.check_approach(approach)?;
        let before = previous
            .evaluate_attached(approach)
            .map_err(NormalLawSourceError::Jkr)?;
        let next = refreshed.source();
        let after = next
            .evaluate_attached(approach)
            .map_err(NormalLawSourceError::Jkr)?;
        if refreshed.original_model_identity() != previous.identity()
            || refreshed.control() != AdhesiveThermalControl::HeldIndentation
            || refreshed.mechanical_work_j() != 0.0
            || *refreshed.before() != before
            || *refreshed.response() != after
            || previous.effective_radius_m() != next.effective_radius_m()
            || previous.interaction_range_m() != next.interaction_range_m()
            || previous.envelope() != next.envelope()
        {
            return Err(invalid("JKR refresh must retain original source, actual held approach, attached branch and geometry").into());
        }
        limit(
            "refreshed native JKR force",
            after.normal_force_n.abs(),
            self.config.maximum_force_n,
        )?;
        let work = finite(refreshed.parameter_work_j())?;
        let residual = finite(
            (Dd::from_f64(after.reversible_energy_j)
                - Dd::from_f64(before.reversible_energy_j)
                - Dd::from_f64(work))
            .to_f64(),
        )?;
        let scale = before
            .reversible_energy_j
            .abs()
            .max(after.reversible_energy_j.abs())
            .max(work.abs());
        let tolerance = finite(
            self.config.contact_energy_absolute_tolerance_j
                + self.config.contact_energy_relative_tolerance * scale,
        )?;
        check_contact_work(residual, tolerance)?;
        check_contact_work(finite(refreshed.potential_change_j() - work)?, tolerance)?;
        let native_energy = self.network.total_energy_j()?;
        let energy_before = finite(native_energy + before.reversible_energy_j)?;
        let energy_after = finite(native_energy + after.reversible_energy_j)?;
        limit(
            "refreshed native/source energy magnitude",
            energy_after.abs(),
            self.network.config.maximum_total_energy_j,
        )?;

        let candidate = BoundNormalContact::new(
            NormalLawSource::JkrAttached(Box::new(next.clone())),
            self.contact.work_column().to_vec(),
            self.contact.gap_m(),
        )
        .map_err(super::ModalCouplingError::ContactLaw)?;
        solve::positive_slope(
            candidate.source(),
            approach,
            approach,
            self.compliance_m_per_n,
        )?;
        let report = SourceNormalJkrRefreshReport {
            sample: self.samples_rendered(),
            approach_m: approach,
            energy_before_j: energy_before,
            energy_after_j: energy_after,
            contact_energy_before_j: before.reversible_energy_j,
            contact_energy_after_j: after.reversible_energy_j,
            parameter_work_j: work,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
        };
        poll(Some(gate))?;
        Ok((candidate, report))
    }
}
