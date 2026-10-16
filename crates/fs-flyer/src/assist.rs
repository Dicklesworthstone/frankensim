//! TrainingAssist + PilotSAS (bead wf-root-guzez.5.16.3, E4.6c-iii).
//! Plan §2.1/§5.1.4: accessibility aids with HARD, VISIBLE authority —
//! a zero-latency SAS rate damper plus attitude command shaping, both
//! clamped to a declared fraction of the canard travel, with an
//! always-set active flag (the HUD renders it; the model exposes it —
//! assisted flight is never silently presented as authentic).
//!
//! ISOLATION DOCTRINE (plan law): assist parameters are calibration-
//! subset material ONLY. The spec carries a `CalibrationSubsetTag`, and
//! the historical-calibration admission REFUSES anything tagged — a
//! typed refusal, executed by the battery, so assist tuning can never
//! leak into historical claims.

use crate::Refusal;

/// Authority cap: |assist command| ≤ frac × canard travel.
pub const MAX_AUTHORITY_FRAC: f64 = 0.5;

/// Marker: this object belongs to the calibration subset and is barred
/// from historical calibration.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CalibrationSubsetTag;

/// The assist configuration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AssistSpec {
    /// SAS pitch-rate damping [rad of canard per rad/s of q].
    pub sas_rate_gain: f64,
    /// TrainingAssist attitude shaping [rad of canard per rad of θ err].
    pub assist_attitude_gain: f64,
    /// Authority as a fraction of the canard travel (0, MAX].
    pub authority_frac: f64,
    /// Calibration-subset marker (the historical path refuses it).
    pub tag: CalibrationSubsetTag,
}

/// The registered v1 assist (accessible-but-honest levels).
pub const ASSIST_V1: AssistSpec = AssistSpec {
    sas_rate_gain: 1.5,
    assist_attitude_gain: 2.0,
    authority_frac: 0.3,
    tag: CalibrationSubsetTag,
};

/// One tick's assist output.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AssistOutput {
    /// Canard-command increment [rad] (adds to the pilot's surface).
    pub dc_assist_rad: f64,
    /// ALWAYS true while the assist system is engaged — the HUD flag.
    pub active: bool,
    /// True when the authority clamp bounded this tick's output.
    pub clamped: bool,
}

impl AssistSpec {
    /// Admit the spec.
    ///
    /// # Errors
    /// `assist-spec-invalid` (non-finite/negative gains; authority
    /// outside (0, [`MAX_AUTHORITY_FRAC`]] — tested at the cap AND one
    /// ulp past).
    pub fn admit(&self) -> Result<(), Refusal> {
        let ok = self.sas_rate_gain.is_finite()
            && self.sas_rate_gain >= 0.0
            && self.assist_attitude_gain.is_finite()
            && self.assist_attitude_gain >= 0.0
            && self.authority_frac.is_finite()
            && self.authority_frac > 0.0
            && self.authority_frac <= MAX_AUTHORITY_FRAC;
        if !ok {
            return Err(Refusal {
                code: "assist-spec-invalid",
                message: format!("{self:?}"),
                ranked_repairs: vec![format!(
                    "non-negative finite gains; authority in (0, {MAX_AUTHORITY_FRAC}]"
                )],
            });
        }
        Ok(())
    }

    /// One tick: assist command from pitch rate + attitude error,
    /// clamped to the authority window. Zero-latency by design (this is
    /// a sim-side aid, not a human model — documented tier).
    ///
    /// # Errors
    /// Admission refusals; `assist-input-invalid` (non-finite inputs or
    /// non-positive travel).
    pub fn apply(
        &self,
        q_rad_s: f64,
        theta_err_rad: f64,
        canard_travel_rad: f64,
    ) -> Result<AssistOutput, Refusal> {
        self.admit()?;
        if !(q_rad_s.is_finite() && theta_err_rad.is_finite())
            || !(canard_travel_rad.is_finite() && canard_travel_rad > 0.0)
        {
            return Err(Refusal {
                code: "assist-input-invalid",
                message: format!(
                    "q {q_rad_s:?}, theta_err {theta_err_rad:?}, travel {canard_travel_rad:?}"
                ),
                ranked_repairs: vec!["check the state plumbing".into()],
            });
        }
        let raw = -self.sas_rate_gain * q_rad_s - self.assist_attitude_gain * theta_err_rad;
        let authority = self.authority_frac * canard_travel_rad;
        let dc = raw.clamp(-authority, authority);
        Ok(AssistOutput {
            dc_assist_rad: dc,
            active: true,
            clamped: raw.abs() > authority,
        })
    }
}

/// The historical-calibration admission for assist-tagged objects:
/// ALWAYS a typed refusal. Any code path assembling historical
/// calibration inputs must route tagged specs through here — assist
/// tuning is calibration-subset material and can never enter a
/// historical claim (plan isolation law).
///
/// # Errors
/// `assist-in-historical-calibration`, unconditionally.
pub fn historical_calibration_admit(
    spec: &AssistSpec,
) -> Result<core::convert::Infallible, Refusal> {
    Err(Refusal {
        code: "assist-in-historical-calibration",
        message: format!(
            "assist spec (authority {}) carries CalibrationSubsetTag — barred from \
             historical calibration",
            spec.authority_frac
        ),
        ranked_repairs: vec![
            "run historical calibration without assist; assist belongs to the calibration \
             subset only"
                .into(),
        ],
    })
}
