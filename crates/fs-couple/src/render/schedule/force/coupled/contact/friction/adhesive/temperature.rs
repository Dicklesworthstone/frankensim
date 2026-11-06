//! Original normal and scalar-junction thermal sources on native friction.
//!
//! These changes hold every native tangential coordinate, material/radiation
//! memory, clock and mechanical observation fixed. The external normal owner
//! must accept the supplied normal update and account for its signed normal
//! parameter/actuator work. That work is never debited again as friction heat.
//! This junction law has no reversible tangential store; its changed capacities
//! affect the next actual native mechanical solve without an imposed impulse.

use super::{
    AdhesiveContactState, AdhesiveModalBinding, CancelGate, ContactBranch, FrictionModalSystem,
    ModalFrictionError, invalid, poll,
};
use fs_blake3::ContentHash;
use fs_contact::normal_patch::{
    adhesive::{
        SourceJkrSphereContact,
        temperature::{
            AdhesiveThermalControl, SourceDmtThermalRefresh, SourceJkrThermalRefresh,
            SourceMaugisThermalRefresh,
        },
    },
    dmt::SourceDmtSphereContact,
    maugis::SourceMaugisSphereContact,
};

/// One sealed original normal-source transition. Normal geometry/control and
/// normal work remain with that owner; this consumer owns scalar junction shear.
#[derive(Clone, Copy)]
pub enum AdhesiveModalSourceRefresh<'a> {
    Jkr {
        original: &'a SourceJkrSphereContact,
        refreshed: &'a SourceJkrThermalRefresh,
    },
    Dmt {
        original: &'a SourceDmtSphereContact,
        refreshed: &'a SourceDmtThermalRefresh,
    },
    Maugis {
        original: &'a SourceMaugisSphereContact,
        refreshed: &'a SourceMaugisThermalRefresh,
    },
}
impl AdhesiveModalSourceRefresh<'_> {
    fn normal_control(self) -> AdhesiveThermalControl {
        match self {
            Self::Jkr { refreshed, .. } => refreshed.control(),
            Self::Dmt { refreshed, .. } => refreshed.control(),
            Self::Maugis { refreshed, .. } => refreshed.control(),
        }
    }
}

/// Material/geometry handoff at an unchanged native sample. This report is not
/// a second normal-work account and introduces no thermal-energy observation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveModalSourceRefreshReport {
    pub sample: u64,
    /// Control declared by the external original normal owner.
    pub normal_control: AdhesiveThermalControl,
    pub original_normal_state_identity: ContentHash,
    pub normal_state_identity: ContentHash,
    pub contact_before: AdhesiveContactState,
    pub contact_after: AdhesiveContactState,
    pub temperature_before_k: Option<f64>,
    pub temperature_after_k: Option<f64>,
}

impl FrictionModalSystem {
    /// Continue the same source-bound scalar junction through an original
    /// normal thermal update. Its predecessor area and repulsion must match the
    /// accepted patch exactly; stale geometry and foreign sources refuse.
    ///
    /// A still-attached junction retains its occupied friction branch. Changed
    /// coefficients/geometry invalidate the load-line cache, while an actual
    /// zero-area endpoint starts a fresh junction. Every native mechanical and
    /// material state, sample count and last frame remains unchanged.
    pub fn apply_adhesive_source_refresh(
        &mut self,
        refresh: AdhesiveModalSourceRefresh<'_>,
        gate: &CancelGate,
    ) -> Result<AdhesiveModalSourceRefreshReport, ModalFrictionError> {
        self.apply_adhesive_source_refresh_with_admission(refresh, gate, |report, _| Ok(*report))
    }

    /// Stage a receiving candidate against both original successor source
    /// states before this owner publishes. The callback must not mutate an
    /// accepted external owner: its returned candidate is installed infallibly
    /// after success. Refusal or late cancellation preserves source, branch,
    /// cache, native history, clock and observations together.
    pub fn apply_adhesive_source_refresh_with_admission<R, E>(
        &mut self,
        refresh: AdhesiveModalSourceRefresh<'_>,
        gate: &CancelGate,
        admit: impl FnOnce(&AdhesiveModalSourceRefreshReport, &AdhesiveModalBinding) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<ModalFrictionError>,
    {
        let (candidate, report) = self
            .prepare_adhesive_source_refresh(refresh, gate)
            .map_err(E::from)?;
        let receiving = admit(&report, &candidate)?;
        poll(Some(gate))
            .map_err(ModalFrictionError::from)
            .map_err(E::from)?;
        // Every fallible source, receiver and cancellation gate has passed.
        if self
            .adhesive
            .as_ref()
            .is_none_or(|previous| previous.law != candidate.law)
        {
            self.load_line = None;
        }
        if report.contact_before.area_m2() == 0.0 || report.contact_after.area_m2() == 0.0 {
            self.branch = ContactBranch::Fresh;
            self.staged_branch = ContactBranch::Fresh;
        }
        self.adhesive = Some(candidate);
        Ok(receiving)
    }

    fn prepare_adhesive_source_refresh(
        &self,
        refresh: AdhesiveModalSourceRefresh<'_>,
        gate: &CancelGate,
    ) -> Result<(AdhesiveModalBinding, AdhesiveModalSourceRefreshReport), ModalFrictionError> {
        poll(Some(gate))?;
        let previous = self
            .adhesive
            .as_ref()
            .and_then(AdhesiveModalBinding::source)
            .ok_or_else(|| {
                invalid("thermal junction continuation requires a retained original source")
            })?;
        let source = match refresh {
            AdhesiveModalSourceRefresh::Jkr {
                original,
                refreshed,
            } => previous.rebind_jkr_temperature(original, refreshed)?,
            AdhesiveModalSourceRefresh::Dmt {
                original,
                refreshed,
            } => previous.rebind_dmt_temperature(original, refreshed)?,
            AdhesiveModalSourceRefresh::Maugis {
                original,
                refreshed,
            } => previous.rebind_maugis_temperature(original, refreshed)?,
        };
        let report = AdhesiveModalSourceRefreshReport {
            sample: self.samples_rendered(),
            normal_control: refresh.normal_control(),
            original_normal_state_identity: previous.normal_state().identity(),
            normal_state_identity: source.normal_state().identity(),
            contact_before: previous.tangential().contact(),
            contact_after: source.tangential().contact(),
            temperature_before_k: previous.temperature_kelvin(),
            temperature_after_k: source.temperature_kelvin(),
        };
        let candidate = AdhesiveModalBinding::from_source(source)?;
        if self.port.law != *candidate.law.law().pressure_law() {
            return Err(
                invalid("thermal junction refresh must retain the original pressure law").into(),
            );
        }
        poll(Some(gate))?;
        Ok((candidate, report))
    }
}
