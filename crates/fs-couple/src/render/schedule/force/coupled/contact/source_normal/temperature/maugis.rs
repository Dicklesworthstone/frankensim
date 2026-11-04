//! Original finite-zone temperature continuation on the complete native network.

use super::super::{
    BoundNormalContact, CancelGate, ModalCouplingError, NormalLawSource, NormalLawSourceError,
    SourceNormalModalError, SourceNormalModalSystem, check_contact_work, finite, invalid,
    limit, poll, solve,
};
use fs_contact::normal_patch::adhesive::temperature::SourceMaugisThermalRefresh;
use fs_math::dd::Dd;

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
        let (candidate, report) = self.prepare_maugis_refresh(refreshed, gate).map_err(E::from)?;
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
        let NormalLawSource::MaugisAttached(previous) = self.contact.source() else {
            return Err(invalid(
                "Maugis temperature continuation requires the original finite-zone potential",
            )
            .into());
        };
        let approach = self.approach(self.components())?;
        self.check_approach(approach)?;
        let before = previous
            .evaluate(approach)
            .map_err(NormalLawSourceError::Maugis)?;
        let next = previous.with_held_indentation_refresh(approach, refreshed)?;
        let after = next
            .evaluate(approach)
            .map_err(NormalLawSourceError::Maugis)?;
        limit(
            "refreshed native Maugis force",
            after.normal_force_n().abs(),
            self.config.maximum_force_n,
        )?;
        let work = finite(refreshed.parameter_work_j())?;
        let residual = finite(
            (Dd::from_f64(after.potential_j())
                - Dd::from_f64(before.potential_j())
                - Dd::from_f64(work))
            .to_f64(),
        )?;
        let scale = before
            .potential_j()
            .abs()
            .max(after.potential_j().abs())
            .max(work.abs());
        let tolerance = finite(
            self.config.contact_energy_absolute_tolerance_j
                + self.config.contact_energy_relative_tolerance * scale,
        )?;
        check_contact_work(residual, tolerance)?;
        check_contact_work(finite(refreshed.potential_change_j() - work)?, tolerance)?;
        let native_energy = self.network.total_energy_j()?;
        let energy_before = finite(native_energy + before.potential_j())?;
        let energy_after = finite(native_energy + after.potential_j())?;
        limit(
            "refreshed native/source energy magnitude",
            energy_after.abs(),
            self.network.config.maximum_total_energy_j,
        )?;

        let candidate = BoundNormalContact::new(
            next.into(),
            self.contact.work_column().to_vec(),
            self.contact.gap_m(),
        )
        .map_err(ModalCouplingError::ContactLaw)?;
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
            contact_energy_before_j: before.potential_j(),
            contact_energy_after_j: after.potential_j(),
            parameter_work_j: work,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
        };
        poll(Some(gate))?;
        Ok((candidate, report))
    }
}
