//! Native driven Hertz--Mindlin dynamics and original ordered fretting wear.
//!
//! Native bodies determine the changing normal load and actual compliant shear.
//! The existing mechanical phase runs once. Its irreversible work supplies the
//! pressure/sliding integral for positive mu; the original exact normal mean
//! and accepted shear travel supply the explicit frictionless limit. Both
//! original Archard inventories are proposed only after native work closes.
//! This remains a small-removal volume diagnostic, without wear geometry,
//! changing material coefficients, debris transport or an extra heat sink.

use super::{
    RigidMindlinContact, RigidMindlinControls, RigidMindlinError, RigidMindlinReport,
    advance_rigid_mindlin_driven,
};
use crate::{MassProperties, RigidBodyState, Vec3, Wrench};
use fs_contact::interface_binding::mindlin::fretting::{
    SourceFrettingError,
    loading::{SourceHertzMindlinWear, SourceHertzMindlinWearAdvance, SourceHertzMindlinWearState},
};
use fs_tribo::{TangentialStep, partial_slip::TangentFrame};

/// Original ordered wear/contact source and its two native material attachments.
#[derive(Clone, Copy, Debug)]
pub struct RigidMindlinWearContact<'a> {
    pub law: &'a SourceHertzMindlinWear,
    pub arms_body: [Vec3; 2],
    pub frame: TangentFrame,
}

/// All native momenta, retained history and two wear inventories publish together.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidMindlinWearCandidate {
    pub bodies_after: [RigidBodyState; 2],
    pub contact: SourceHertzMindlinWearAdvance,
    /// Unchanged original native mechanics, including independent external work.
    pub report: RigidMindlinReport,
}

/// No failed mechanical or either-surface wear candidate changes accepted state.
#[derive(Debug)]
pub enum RigidMindlinWearError {
    Mechanical(RigidMindlinError),
    Wear(SourceFrettingError),
    Cancelled,
}
impl From<RigidMindlinError> for RigidMindlinWearError {
    fn from(value: RigidMindlinError) -> Self {
        match value {
            RigidMindlinError::Cancelled => Self::Cancelled,
            value => Self::Mechanical(value),
        }
    }
}
impl From<SourceFrettingError> for RigidMindlinWearError {
    fn from(value: SourceFrettingError) -> Self {
        match value {
            SourceFrettingError::Cancelled => Self::Cancelled,
            value => Self::Wear(value),
        }
    }
}
impl core::fmt::Display for RigidMindlinWearError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native Hertz--Mindlin fretting: {self:?}")
    }
}
impl std::error::Error for RigidMindlinWearError {}

/// Advance unforced native contact and its actual local-sliding wear once.
pub fn advance_rigid_mindlin_wear(
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    contact: RigidMindlinWearContact<'_>,
    state: &SourceHertzMindlinWearState,
    dt_s: f64,
    controls: RigidMindlinControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RigidMindlinWearCandidate, RigidMindlinWearError> {
    advance_rigid_mindlin_driven_wear(
        bodies,
        properties,
        contact,
        state,
        [Wrench::ZERO; 2],
        dt_s,
        controls,
        cancelled,
    )
}

/// Advance held external wrenches, changing normal/shear history and both volumes.
///
/// The original native phase owns body response, source rate limits, normal and
/// shear compatibility, storage, dissipation and world momentum. The added wear
/// calculation does not change its force, heat or work. Separation, unsupported
/// motion and changing thermal/material conditions retain the native/source
/// refusals. A late overflow or cancellation publishes neither bodies nor wear.
#[allow(clippy::too_many_arguments)]
pub fn advance_rigid_mindlin_driven_wear(
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    contact: RigidMindlinWearContact<'_>,
    state: &SourceHertzMindlinWearState,
    external: [Wrench; 2],
    dt_s: f64,
    controls: RigidMindlinControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RigidMindlinWearCandidate, RigidMindlinWearError> {
    let native = advance_rigid_mindlin_driven(
        bodies,
        properties,
        RigidMindlinContact {
            law: contact.law.contact(),
            arms_body: contact.arms_body,
            frame: contact.frame,
        },
        state.contact_state(),
        external,
        dt_s,
        controls,
        cancelled,
    )?;
    if cancelled() {
        return Err(RigidMindlinWearError::Cancelled);
    }
    let path = native.contact;
    let mu = contact
        .law
        .contact()
        .law()
        .parameters()
        .friction_coefficient;
    let measure = if mu > 0.0 {
        path.tangential.dissipated_j / mu
    } else {
        native.report.normal_average_force_n * path.tangential.slip_increment_m[0].abs()
    };
    // The original Archard leaf checks finite/nonnegative measure and both
    // cumulative sums before either body or inventory can be accepted.
    let (volumes, wear) = contact.law.advance_volume_candidates(state, measure)?;
    let shear = path.tangential;
    let next = contact.law.state_from_contact(shear.state, volumes)?;
    if cancelled() {
        return Err(RigidMindlinWearError::Cancelled);
    }
    Ok(RigidMindlinWearCandidate {
        bodies_after: native.bodies_after,
        report: native.report,
        contact: SourceHertzMindlinWearAdvance {
            tangential: TangentialStep {
                traction_n: shear.traction_n,
                slip_increment_m: shear.slip_increment_m,
                state: next,
                regime: shear.regime,
                stored_change_j: shear.stored_change_j,
                dissipated_j: shear.dissipated_j,
                closure_residual_j: shear.closure_residual_j,
                compatibility_residual_m: shear.compatibility_residual_m,
                iterations: shear.iterations,
            },
            normal_response: path.normal_response,
            normal_work_j: path.normal_work_j,
            normal_stored_change_j: path.normal_stored_change_j,
            wear,
            pressure_weighted_sliding_j: measure,
        },
    })
}
