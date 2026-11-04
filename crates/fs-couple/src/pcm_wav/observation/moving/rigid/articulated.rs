//! Original articulated material points and link-frame sources on one clock.
//!
//! Call `prepare_point_motion` on the ORIGINAL inverse/forward dynamics receipt
//! and its matching kinematics before observing a link. Prescribed-base receipts
//! require their original gravity vector; free-floating receipts already carry
//! physical acceleration. This adapter must neither add gravity again nor use
//! the body's spatial-acceleration linear part as Cartesian point acceleration.
//!
//! The attachment is relative to the native link origin, NOT its centre of mass.
//! Joint transport, Coriolis, tangential and centripetal motion belong to
//! `PreparedLinkPointMotion`. All sources/microphones on a link can reuse that
//! one prepared snapshot. No model traversal, dynamics solve, finite difference
//! or allocation occurs here. Source force remains explicitly force ON FLUID;
//! actuator efforts and support reactions are not converted into radiation.
//!
//! Histories retain the same physical world frame, clock and bandwidth contract
//! as standalone rigid histories. They feed the existing fixed/moving-receiver
//! trajectory and pressure owners directly. Laboratory histories in uniform
//! flow use the existing explicit `PreparedKinematicTrajectory::in_uniform_flow`
//! conversion. Observation must not interpolate through an unresolved impact.

use super::{
    BodyFrameCompactSourceStrength, CompactSourceStrength, RigidAcousticAttachment,
    RigidAcousticError, RigidAcousticHistory, RigidAcousticSample, SourceKinematics,
    TimedCompactSource, TimedSourceKinematics, nonfinite, poll, validate_strength,
};
use fs_ga::Vec3;
use fs_mbd::articulated::{ArticulatedError, ArticulatedPointMotion, PreparedLinkPointMotion};

/// Coherent, allocation-free publication of a complete accepted link frame.
pub mod frame;
pub use frame::{LinkAcousticUpdate, append_link_frame, append_link_frame_with_relative_motion};

/// Preserve the original articulated refusal without changing the existing
/// standalone-rigid error API or replacing native failures with generic text.
#[derive(Debug, Clone, PartialEq)]
pub enum ArticulatedAcousticError {
    /// The original articulated point or spatial rotation refused.
    Dynamics(ArticulatedError),
    /// Original source, clock, capacity or cancellation admission refused.
    Observation(RigidAcousticError),
}

impl From<ArticulatedError> for ArticulatedAcousticError {
    fn from(error: ArticulatedError) -> Self {
        Self::Dynamics(error)
    }
}
impl From<RigidAcousticError> for ArticulatedAcousticError {
    fn from(error: RigidAcousticError) -> Self {
        Self::Observation(error)
    }
}
impl core::fmt::Display for ArticulatedAcousticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Dynamics(error) => write!(f, "articulated acoustic motion: {error}"),
            Self::Observation(error) => write!(f, "articulated acoustic observation: {error}"),
        }
    }
}
impl std::error::Error for ArticulatedAcousticError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Dynamics(error) => Some(error),
            Self::Observation(error) => Some(error),
        }
    }
}

impl RigidAcousticAttachment {
    /// Observe a source or ideal microphone fixed to an articulated link.
    ///
    /// `arm_body_m` is measured from that LINK's original coordinate origin.
    /// The supplied snapshot must describe the same link and accepted physical
    /// event on every call. Native mechanical receipts do not own timestamps.
    pub fn kinematics_link(
        self,
        time_s: f64,
        motion: &PreparedLinkPointMotion,
    ) -> Result<TimedSourceKinematics, ArticulatedAcousticError> {
        Ok(timed(time_s, self.link_point(time_s, motion)?))
    }

    /// Transport an explicit link-frame source using the native link rotation.
    ///
    /// The derivative of link-frame COMPONENTS is transported as
    /// `R Fdot_link + omega_world cross (R F_link)`. Even a force constant in
    /// link axes therefore has a nonzero world derivative when the link rotates.
    /// Point acceleration and force derivative remain separate physical inputs.
    pub fn sample_link(
        self,
        time_s: f64,
        motion: &PreparedLinkPointMotion,
        strength: BodyFrameCompactSourceStrength,
    ) -> Result<RigidAcousticSample, ArticulatedAcousticError> {
        validate_strength(strength)?;
        let point = self.link_point(time_s, motion)?;
        sample_at_point(time_s, point, motion, strength)
    }

    fn link_point(
        self,
        time_s: f64,
        motion: &PreparedLinkPointMotion,
    ) -> Result<ArticulatedPointMotion, ArticulatedAcousticError> {
        if !time_s.is_finite() {
            return Err(nonfinite("articulated acoustic accepted timestamp").into());
        }
        let arm = self.arm_body_m();
        Ok(motion.point_motion(Vec3::new(arm.x, arm.y, arm.z))?)
    }
}

// Both fixed and migrating points use this one force/derivative transport.
// The native point query has already run; no duplicate motion evaluation.
pub(super) fn sample_at_point(
    time_s: f64,
    point: ArticulatedPointMotion,
    motion: &PreparedLinkPointMotion,
    strength: BodyFrameCompactSourceStrength,
) -> Result<RigidAcousticSample, ArticulatedAcousticError> {
    let rotation = motion.world_from_link().rotation();
    let force = rotation
        .rotate(vector(strength.force_on_fluid_body_n))
        .map_err(ArticulatedError::from)?;
    let derivative = rotation
        .rotate(vector(strength.force_rate_on_fluid_body_n_s))
        .map_err(ArticulatedError::from)?
        + point.angular_velocity_world.cross(force);
    let force = components(force);
    let derivative = components(derivative);
    if force.iter().chain(&derivative).any(|v| !v.is_finite()) {
        return Err(nonfinite("articulated acoustic force transport derivative").into());
    }
    Ok(RigidAcousticSample {
        kinematics: timed(time_s, point),
        source: TimedCompactSource {
            source_time_s: time_s,
            strength: CompactSourceStrength {
                volume_flow_m3_s: strength.volume_flow_m3_s,
                volume_flow_rate_m3_s2: strength.volume_flow_rate_m3_s2,
                force_on_fluid_n: force,
                force_rate_on_fluid_n_s: derivative,
            },
        },
    })
}

impl RigidAcousticHistory {
    /// Record one accepted articulated event into the ORIGINAL paired history.
    ///
    /// Use one fixed link identity and link-relative attachment throughout the
    /// history. Prepare native dynamics once per link/event and reuse it across
    /// attachments. This method performs no solve or append-time allocation.
    /// Invalid source/motion, a full history, a non-increasing clock or even late
    /// cancellation publishes neither endpoint and never mutates mechanics.
    pub fn append_link(
        &mut self,
        time_s: f64,
        motion: &PreparedLinkPointMotion,
        strength: BodyFrameCompactSourceStrength,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ArticulatedAcousticError> {
        poll(cancelled)?;
        self.admit_append(time_s)?;
        let sample = self.attachment.sample_link(time_s, motion, strength)?;
        poll(cancelled)?;
        self.publish(sample);
        Ok(())
    }
}

pub(super) fn timed(time_s: f64, point: ArticulatedPointMotion) -> TimedSourceKinematics {
    TimedSourceKinematics {
        source_time_s: time_s,
        kinematics: SourceKinematics {
            position_m: components(point.point_world),
            velocity_m_s: components(point.point_velocity_world),
            acceleration_m_s2: components(point.point_acceleration_world),
        },
    }
}
fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
