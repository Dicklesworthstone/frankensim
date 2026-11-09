//! Prescribed JKR attachment/separation with retained bonded-twist energy.
//!
//! Normal travel holds the original homogeneous elastic twist. The normal
//! owner supplies its mode-I branch, threshold and smooth work. Before a
//! separation, the original torsional spring reaches that prescribed fold
//! radius; its remaining storage becomes outgoing mechanical release,
//! separately from normal snap and irreversible twisting loss. A new junction
//! starts unstrained, with no preload from detached angular travel.
//!
//! This retains the parent's prescribed-area approximation and its explicit
//! G_III/W screen. It does not locate a coupled I/III instability, recruit
//! stress-free annuli, assign a heat destination or produce a body impulse.

use super::{
    AdhesiveTorsionError, JkrTorsionControls, SourceJkrTorsionContact, SourceJkrTorsionError,
    SourceJkrTorsionState, SourceJkrTorsionStep,
};
use crate::normal_patch::adhesive::{
    JkrBranch, JkrResponse, JkrTransition, SourceJkrSphereContact,
};
use fs_blake3::ContentHash;
use fs_math::dd::Dd;

/// Original temperature continuation of the complete attached lifecycle.
pub mod temperature;

/// One original normal source and calibrated torsional law across junctions.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceJkrTorsionEvents {
    // An admitted seed supplies unstrained initial memory at later junctions;
    // immutable original source receipts stay shared across radius changes.
    seed: SourceJkrTorsionContact,
}

/// Sealed branch, optional loaded spring and separate accumulated energy losses.
#[derive(Clone, Debug, PartialEq)]
pub struct JkrTorsionEventState {
    identity: ContentHash,
    normal: JkrResponse,
    attached: Option<(SourceJkrTorsionContact, SourceJkrTorsionState)>,
    torsion_storage_j: f64,
    normal_release_j: Dd,
    torsion_release_j: Dd,
    torsion_loss_j: Dd,
}
impl JkrTorsionEventState {
    pub const fn normal_response(&self) -> &JkrResponse {
        &self.normal
    }
    pub fn attached(&self) -> Option<(&SourceJkrTorsionContact, SourceJkrTorsionState)> {
        self.attached
            .as_ref()
            .map(|(source, state)| (source, *state))
    }
    pub const fn torsion_storage_j(&self) -> f64 {
        self.torsion_storage_j
    }
    pub const fn accumulated_normal_release_j(&self) -> f64 {
        self.normal_release_j.to_f64()
    }
    pub const fn accumulated_torsion_release_j(&self) -> f64 {
        self.torsion_release_j.to_f64()
    }
    pub fn accumulated_released_energy_j(&self) -> f64 {
        (self.normal_release_j + self.torsion_release_j).to_f64()
    }
    /// Plastic twisting loss, distinct from normal and remaining-spring release.
    pub const fn accumulated_torsion_loss_j(&self) -> f64 {
        self.torsion_loss_j.to_f64()
    }
}

/// Original local result and its distinct physical work/storage/release channels.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrTorsionEventReport {
    /// Present only for prescribed normal travel, including its original jump.
    pub normal_transition: Option<JkrTransition>,
    /// Present only for twisting an attached junction.
    pub torsion: Option<SourceJkrTorsionStep>,
    pub normal_contact_work_j: f64,
    pub torsion_parameter_work_j: f64,
    pub normal_released_energy_j: f64,
    pub torsion_released_energy_j: f64,
    pub irreversible_torsion_loss_j: f64,
    /// Original incremental normal-plus-torsional storage change [J].
    pub storage_change_j: f64,
    pub applied_work_j: f64,
    /// Normal snap plus remaining spring energy, with no destination inferred.
    pub released_energy_j: f64,
    /// Delta storage + event release + plastic loss - applied work [J].
    pub energy_residual_j: f64,
}

/// Immutable source proposal; writable diagnostics cannot become accepted state.
#[derive(Clone, Debug, PartialEq)]
pub struct JkrTorsionEventStep {
    previous: JkrTorsionEventState,
    state: JkrTorsionEventState,
    report: JkrTorsionEventReport,
}
impl JkrTorsionEventStep {
    pub const fn previous_state(&self) -> &JkrTorsionEventState {
        &self.previous
    }
    pub const fn state(&self) -> &JkrTorsionEventState {
        &self.state
    }
    pub const fn report(&self) -> &JkrTorsionEventReport {
        &self.report
    }
    pub fn into_parts(self) -> (JkrTorsionEventState, JkrTorsionEventReport) {
        (self.state, self.report)
    }
}

impl SourceJkrTorsionEvents {
    /// Admit the original material, twist calibration and G_III/W ceiling at
    /// attachment. Every later radius is checked again by the original owner.
    pub fn new(
        normal: SourceJkrSphereContact,
        controls: JkrTorsionControls,
    ) -> Result<Self, SourceJkrTorsionError> {
        Ok(Self {
            seed: SourceJkrTorsionContact::new(normal, 0.0, controls)?,
        })
    }
    pub fn normal_source(&self) -> &SourceJkrSphereContact {
        self.seed.normal_source()
    }
    pub const fn controls(&self) -> JkrTorsionControls {
        self.seed.controls()
    }

    /// Start an admitted attached contact without preload or prior release/loss.
    pub fn initial_attached_state(
        &self,
        indentation_m: f64,
    ) -> Result<JkrTorsionEventState, SourceJkrTorsionError> {
        let attached = self.fresh_attached(indentation_m)?;
        let mut state = self.initial(*attached.0.normal_response());
        state.attached = Some(attached);
        Ok(state)
    }
    /// A detached state at zero approach attaches on its next normal advance.
    pub fn initial_detached_state(
        &self,
        indentation_m: f64,
    ) -> Result<JkrTorsionEventState, SourceJkrTorsionError> {
        self.normal_source().detached_state(indentation_m)?;
        Ok(self.initial(self.normal_source().detached_response(indentation_m)))
    }

    /// One monotone prescribed normal segment at held elastic twist. Split
    /// reversals into separate calls. Original JKR thresholds, including exact
    /// equality, determine the event even when the endpoint overshoots it.
    /// A departing spring reaches the admitted mode-I fold before release.
    pub fn advance_normal(
        &self,
        previous: &JkrTorsionEventState,
        indentation_m: f64,
    ) -> Result<JkrTorsionEventStep, SourceJkrTorsionError> {
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
                    .ok_or(SourceJkrTorsionError::StateMismatch)?
                    .indentation_m
            };
            let refresh = source.advance_normal(&memory, target)?;
            report.torsion_parameter_work_j = refresh.torsion_parameter_work_j();
            if normal.response.branch == JkrBranch::Attached {
                state.torsion_storage_j = refresh.torsion_storage_after_j();
                state.attached = Some(refresh.into_parts());
            } else {
                report.torsion_released_energy_j = refresh.torsion_storage_after_j();
                state.torsion_storage_j = 0.0;
                state.attached = None;
            }
        } else if normal.response.branch == JkrBranch::Attached {
            state.attached = Some(self.fresh_attached(indentation_m)?);
        }
        Self::finish(previous, state, report)
    }

    /// One monotone signed twist on the current patch, through the original
    /// exact elastic/plastic return map. Detached travel produces no contact
    /// work and cannot preload a later junction. Zero angular travel leaves a
    /// marginal attached state attached; a normal call owns its branch event.
    pub fn advance_torsion(
        &self,
        previous: &JkrTorsionEventState,
        angle_increment_rad: f64,
    ) -> Result<JkrTorsionEventStep, SourceJkrTorsionError> {
        self.admit(previous)?;
        if !angle_increment_rad.is_finite() {
            return Err(AdhesiveTorsionError::InvalidInput {
                field: "torsion angle increment",
            }
            .into());
        }
        let mut state = previous.clone();
        let mut report = empty_report();
        if let Some((source, memory)) = previous.attached() {
            let step = source.advance(&memory, angle_increment_rad)?;
            state.torsion_storage_j = step.torsion().storage_after_j;
            state.attached = Some((source.clone(), step.state()));
            report.torsion = Some(step);
            report.irreversible_torsion_loss_j = step.torsion().irreversible_loss_j;
        }
        Self::finish(previous, state, report)
    }

    fn fresh_attached(
        &self,
        indentation_m: f64,
    ) -> Result<(SourceJkrTorsionContact, SourceJkrTorsionState), SourceJkrTorsionError> {
        Ok(self
            .seed
            .advance_normal(&self.seed.initial_state(), indentation_m)?
            .into_parts())
    }
    fn initial(&self, normal: JkrResponse) -> JkrTorsionEventState {
        JkrTorsionEventState {
            identity: self.seed.identity(),
            normal,
            attached: None,
            torsion_storage_j: 0.0,
            normal_release_j: Dd::ZERO,
            torsion_release_j: Dd::ZERO,
            torsion_loss_j: Dd::ZERO,
        }
    }
    fn admit(&self, state: &JkrTorsionEventState) -> Result<(), SourceJkrTorsionError> {
        if state.identity != self.seed.identity() {
            return Err(SourceJkrTorsionError::StateMismatch);
        }
        Ok(())
    }
    fn finish(
        previous: &JkrTorsionEventState,
        mut state: JkrTorsionEventState,
        mut report: JkrTorsionEventReport,
    ) -> Result<JkrTorsionEventStep, SourceJkrTorsionError> {
        let torsion = report.torsion.as_ref().map(SourceJkrTorsionStep::torsion);
        // Keep the original incremental work primitives; subtraction of two
        // rounded absolute stores would lose weak work beside a large preload.
        let storage = dd(report.normal_contact_work_j) - dd(report.normal_released_energy_j)
            + dd(report.torsion_parameter_work_j)
            - dd(report.torsion_released_energy_j)
            + dd(torsion.map_or(0.0, |step| step.storage_change_j));
        let work = dd(report.normal_contact_work_j)
            + dd(report.torsion_parameter_work_j)
            + dd(torsion.map_or(0.0, |step| step.applied_work_j));
        let release = dd(report.normal_released_energy_j) + dd(report.torsion_released_energy_j);
        report.storage_change_j = storage.to_f64();
        report.applied_work_j = work.to_f64();
        report.released_energy_j = release.to_f64();
        report.energy_residual_j =
            (storage + release + dd(report.irreversible_torsion_loss_j) - work).to_f64();
        state.normal_release_j = previous.normal_release_j + dd(report.normal_released_energy_j);
        state.torsion_release_j = previous.torsion_release_j + dd(report.torsion_released_energy_j);
        state.torsion_loss_j = previous.torsion_loss_j + dd(report.irreversible_torsion_loss_j);
        if [
            report.storage_change_j,
            report.applied_work_j,
            report.released_energy_j,
            report.energy_residual_j,
            state.accumulated_normal_release_j(),
            state.accumulated_torsion_release_j(),
            state.accumulated_released_energy_j(),
            state.accumulated_torsion_loss_j(),
        ]
        .iter()
        .any(|value| !value.is_finite())
        {
            return Err(AdhesiveTorsionError::Unrepresentable {
                field: "torsion lifecycle energy",
            }
            .into());
        }
        Ok(JkrTorsionEventStep {
            previous: previous.clone(),
            state,
            report,
        })
    }
}

fn empty_report() -> JkrTorsionEventReport {
    JkrTorsionEventReport {
        normal_transition: None,
        torsion: None,
        normal_contact_work_j: 0.0,
        torsion_parameter_work_j: 0.0,
        normal_released_energy_j: 0.0,
        torsion_released_energy_j: 0.0,
        irreversible_torsion_loss_j: 0.0,
        storage_change_j: 0.0,
        applied_work_j: 0.0,
        released_energy_j: 0.0,
        energy_residual_j: 0.0,
    }
}
fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}
