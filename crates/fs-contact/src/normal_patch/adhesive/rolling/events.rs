//! Original JKR attachment and separation with retained rolling storage.
//!
//! Normal travel is prescribed at fixed rolling angle. The original normal
//! owner selects its branch and integrates its smooth work. Before separation,
//! the original rolling owner reaches the actual displacement spinodal radius;
//! its remaining spring energy is then outgoing mechanical release, separately
//! from the normal snap and from irreversible rolling loss. Attachment starts
//! a stress-free junction. Detached angular travel cannot preload that junction.
//! This extends the existing one-way normal-area approximation, not reciprocal
//! mixed-mode equilibrium, a heat partition, a body impulse or contact discovery.

use super::{
    AdhesiveRollingError, JkrRollingControls, SourceJkrRollingContact, SourceJkrRollingError,
    SourceJkrRollingState, SourceJkrRollingStep,
};
use crate::normal_patch::adhesive::{
    JkrBranch, JkrResponse, JkrTransition, SourceJkrSphereContact,
};
use fs_blake3::ContentHash;
use fs_math::dd::Dd;

/// Original temperature continuation of the complete attached lifecycle.
pub mod temperature;

/// One original normal source and calibration shared by all junction lifetimes.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceJkrRollingEvents {
    // An admitted seed supplies unstrained initial memory at later junctions;
    // immutable original source receipts stay shared across radius changes.
    seed: SourceJkrRollingContact,
}

/// Sealed branch, attached spring and separate accumulated outgoing energies.
#[derive(Clone, Debug, PartialEq)]
pub struct JkrRollingEventState {
    identity: ContentHash,
    normal: JkrResponse,
    attached: Option<(SourceJkrRollingContact, SourceJkrRollingState)>,
    rolling_storage_j: f64,
    normal_release_j: Dd,
    rolling_release_j: Dd,
    rolling_loss_j: Dd,
}
impl JkrRollingEventState {
    pub const fn normal_response(&self) -> &JkrResponse {
        &self.normal
    }
    pub fn attached(&self) -> Option<(&SourceJkrRollingContact, SourceJkrRollingState)> {
        self.attached
            .as_ref()
            .map(|(source, state)| (source, *state))
    }
    pub const fn rolling_storage_j(&self) -> f64 {
        self.rolling_storage_j
    }
    pub const fn accumulated_normal_release_j(&self) -> f64 {
        self.normal_release_j.to_f64()
    }
    pub const fn accumulated_rolling_release_j(&self) -> f64 {
        self.rolling_release_j.to_f64()
    }
    pub fn accumulated_released_energy_j(&self) -> f64 {
        (self.normal_release_j + self.rolling_release_j).to_f64()
    }
    /// Plastic rolling loss stays distinct from branch-event release.
    pub const fn accumulated_rolling_loss_j(&self) -> f64 {
        self.rolling_loss_j.to_f64()
    }
}

/// Original local results and their distinct physical work channels.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrRollingEventReport {
    /// Present only for a prescribed normal segment, including its original jump.
    pub normal_transition: Option<JkrTransition>,
    /// Present only for angular travel on an attached junction.
    pub rolling: Option<SourceJkrRollingStep>,
    pub normal_contact_work_j: f64,
    pub rolling_parameter_work_j: f64,
    pub normal_released_energy_j: f64,
    pub rolling_released_energy_j: f64,
    pub irreversible_rolling_loss_j: f64,
    /// Factored contact-plus-rolling change, preserving small smooth work.
    pub storage_change_j: f64,
    pub applied_work_j: f64,
    /// Normal snap plus released spring storage; no thermal destination inferred.
    pub released_energy_j: f64,
    /// Delta storage + event release + plastic loss - applied work.
    pub energy_residual_j: f64,
}

/// Immutable proposal; diagnostics cannot be rewritten into accepted history.
#[derive(Clone, Debug, PartialEq)]
pub struct JkrRollingEventStep {
    previous: JkrRollingEventState,
    state: JkrRollingEventState,
    report: JkrRollingEventReport,
}
impl JkrRollingEventStep {
    pub const fn previous_state(&self) -> &JkrRollingEventState {
        &self.previous
    }
    pub const fn state(&self) -> &JkrRollingEventState {
        &self.state
    }
    pub const fn report(&self) -> &JkrRollingEventReport {
        &self.report
    }
    pub fn into_parts(self) -> (JkrRollingEventState, JkrRollingEventReport) {
        (self.state, self.report)
    }
}

impl SourceJkrRollingEvents {
    /// Admit the original source and rolling calibration at attachment.
    pub fn new(
        normal: SourceJkrSphereContact,
        controls: JkrRollingControls,
    ) -> Result<Self, SourceJkrRollingError> {
        Ok(Self {
            seed: SourceJkrRollingContact::new(normal, 0.0, controls)?,
        })
    }
    pub fn normal_source(&self) -> &SourceJkrSphereContact {
        self.seed.normal_source()
    }
    pub const fn controls(&self) -> JkrRollingControls {
        self.seed.controls()
    }
    /// Start an admitted attached contact with no rolling preload or prior loss.
    pub fn initial_attached_state(
        &self,
        indentation_m: f64,
    ) -> Result<JkrRollingEventState, SourceJkrRollingError> {
        let attached = self.fresh_attached(indentation_m)?;
        let mut state = self.initial(*attached.0.normal_response());
        state.attached = Some(attached);
        Ok(state)
    }
    /// Zero approach may be initialized detached; the next normal call attaches.
    pub fn initial_detached_state(
        &self,
        indentation_m: f64,
    ) -> Result<JkrRollingEventState, SourceJkrRollingError> {
        self.normal_source().detached_state(indentation_m)?;
        Ok(self.initial(self.normal_source().detached_response(indentation_m)))
    }

    /// One monotone normal segment at held rolling angle. Split reversals into
    /// separate calls. Original event thresholds own work even on an overshoot.
    pub fn advance_normal(
        &self,
        previous: &JkrRollingEventState,
        indentation_m: f64,
    ) -> Result<JkrRollingEventStep, SourceJkrRollingError> {
        self.admit(previous)?;
        let normal = self
            .normal_source()
            .advance_monotone(previous.normal.state(), indentation_m)?;
        let mut state = previous.clone();
        state.normal = normal.response;
        let mut report = empty_report();
        report.normal_transition = Some(normal);
        report.normal_contact_work_j = normal.prescribed_work_j;
        report.normal_released_energy_j = normal.irreversible_work_j;
        if let Some((source, memory)) = previous.attached() {
            let target = if normal.response.branch == JkrBranch::Attached {
                indentation_m
            } else {
                normal
                    .jump
                    .ok_or(SourceJkrRollingError::StateMismatch)?
                    .indentation_m
            };
            let refresh = source.advance_normal(&memory, target)?;
            report.rolling_parameter_work_j = refresh.rolling_parameter_work_j();
            if normal.response.branch == JkrBranch::Attached {
                state.rolling_storage_j = refresh.rolling_storage_after_j();
                state.attached = Some(refresh.into_parts());
            } else {
                report.rolling_released_energy_j = refresh.rolling_storage_after_j();
                state.rolling_storage_j = 0.0;
                state.attached = None;
            }
        } else if normal.response.branch == JkrBranch::Attached {
            state.attached = Some(self.fresh_attached(indentation_m)?);
        }
        Self::finish(previous, state, report)
    }

    /// Prescribe scalar rolling on the current patch using the original return
    /// map. Detached travel has zero contact work and no retained angular memory.
    /// This operation supplies no body motion or normal evolution.
    pub fn advance_rolling(
        &self,
        previous: &JkrRollingEventState,
        angle_increment_rad: f64,
    ) -> Result<JkrRollingEventStep, SourceJkrRollingError> {
        self.admit(previous)?;
        if !angle_increment_rad.is_finite() {
            return Err(AdhesiveRollingError::InvalidInput {
                field: "rolling angle increment",
            }
            .into());
        }
        let mut state = previous.clone();
        let mut report = empty_report();
        if let Some((source, memory)) = previous.attached() {
            let step = source.advance(&memory, angle_increment_rad)?;
            state.rolling_storage_j = step.rolling().storage_after_j;
            state.attached = Some((source.clone(), step.state()));
            report.rolling = Some(step);
            report.irreversible_rolling_loss_j = step.rolling().irreversible_loss_j;
        }
        Self::finish(previous, state, report)
    }

    fn fresh_attached(
        &self,
        indentation_m: f64,
    ) -> Result<(SourceJkrRollingContact, SourceJkrRollingState), SourceJkrRollingError> {
        Ok(self
            .seed
            .advance_normal(&self.seed.initial_state(), indentation_m)?
            .into_parts())
    }
    fn initial(&self, normal: JkrResponse) -> JkrRollingEventState {
        JkrRollingEventState {
            identity: self.seed.identity(),
            normal,
            attached: None,
            rolling_storage_j: 0.0,
            normal_release_j: Dd::ZERO,
            rolling_release_j: Dd::ZERO,
            rolling_loss_j: Dd::ZERO,
        }
    }
    fn admit(&self, state: &JkrRollingEventState) -> Result<(), SourceJkrRollingError> {
        if state.identity != self.seed.identity() {
            return Err(SourceJkrRollingError::StateMismatch);
        }
        Ok(())
    }
    fn finish(
        previous: &JkrRollingEventState,
        mut state: JkrRollingEventState,
        mut report: JkrRollingEventReport,
    ) -> Result<JkrRollingEventStep, SourceJkrRollingError> {
        let rolling = report.rolling.as_ref().map(SourceJkrRollingStep::rolling);
        let storage = dd(report.normal_contact_work_j) - dd(report.normal_released_energy_j)
            + dd(report.rolling_parameter_work_j)
            - dd(report.rolling_released_energy_j)
            + dd(rolling.map_or(0.0, |step| step.storage_change_j));
        let work = dd(report.normal_contact_work_j)
            + dd(report.rolling_parameter_work_j)
            + dd(rolling.map_or(0.0, |step| step.applied_work_j));
        let release = dd(report.normal_released_energy_j) + dd(report.rolling_released_energy_j);
        report.storage_change_j = storage.to_f64();
        report.applied_work_j = work.to_f64();
        report.released_energy_j = release.to_f64();
        report.energy_residual_j =
            (storage + release + dd(report.irreversible_rolling_loss_j) - work).to_f64();
        state.normal_release_j = previous.normal_release_j + dd(report.normal_released_energy_j);
        state.rolling_release_j = previous.rolling_release_j + dd(report.rolling_released_energy_j);
        state.rolling_loss_j = previous.rolling_loss_j + dd(report.irreversible_rolling_loss_j);
        if [
            report.storage_change_j,
            report.applied_work_j,
            report.released_energy_j,
            report.energy_residual_j,
            state.accumulated_normal_release_j(),
            state.accumulated_rolling_release_j(),
            state.accumulated_released_energy_j(),
            state.accumulated_rolling_loss_j(),
        ]
        .iter()
        .any(|value| !value.is_finite())
        {
            return Err(AdhesiveRollingError::Unrepresentable {
                field: "rolling lifecycle energy",
            }
            .into());
        }
        Ok(JkrRollingEventStep {
            previous: previous.clone(),
            state,
            report,
        })
    }
}

fn empty_report() -> JkrRollingEventReport {
    JkrRollingEventReport {
        normal_transition: None,
        rolling: None,
        normal_contact_work_j: 0.0,
        rolling_parameter_work_j: 0.0,
        normal_released_energy_j: 0.0,
        rolling_released_energy_j: 0.0,
        irreversible_rolling_loss_j: 0.0,
        storage_change_j: 0.0,
        applied_work_j: 0.0,
        released_energy_j: 0.0,
        energy_residual_j: 0.0,
    }
}
fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}
