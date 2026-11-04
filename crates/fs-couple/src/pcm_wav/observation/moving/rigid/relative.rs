//! Explicit source/probe migration in native carrier coordinates.
//!
//! The relative coordinates and both derivatives describe one smooth physical
//! point path. They are not inferred from changes between accepted samples.
//! Force derivatives are TOTAL derivatives along that path in the carrier axes;
//! a spatially varying source law must already include its convective derivative.
//! This module adds no flow/source law or rotation relative to the carrier.

use super::{
    BodyFrameCompactSourceStrength, PreparedRigidPointMotion, RigidAcousticAttachment,
    RigidAcousticError, RigidAcousticHistory, RigidAcousticSample, SourceKinematics,
    TimedSourceKinematics, Vec3, components, nonfinite, poll, sample_at_kinematics,
    validate_strength,
};
use super::articulated::{ArticulatedAcousticError, sample_at_point};
use fs_mbd::articulated::PreparedLinkPointMotion;
use fs_ga::Vec3 as LinkVec3;

/// Motion relative to a history's fixed body-coordinate reference attachment.
///
/// All derivatives are taken in the rotating carrier coordinates at the SAME
/// physical timestamp. The standalone carrier uses its COM; an articulated
/// carrier uses its link origin. Values are validated by each sampling method.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RelativeAcousticMotion {
    /// Displacement from the original reference attachment [m].
    pub displacement_body_m: Vec3,
    /// Time derivative of that displacement in body/link components [m/s].
    pub velocity_body_m_s: Vec3,
    /// Second derivative of those same coordinate components [m/s^2].
    pub acceleration_body_m_s2: Vec3,
}

impl RelativeAcousticMotion {
    fn arm(self, attachment: RigidAcousticAttachment) -> Result<Vec3, RigidAcousticError> {
        if !self.displacement_body_m.is_finite()
            || !self.velocity_body_m_s.is_finite()
            || !self.acceleration_body_m_s2.is_finite()
        {
            return Err(nonfinite("relative acoustic coordinates and derivatives"));
        }
        let arm = attachment.arm_body_m().add(self.displacement_body_m);
        if !arm.is_finite() {
            return Err(nonfinite("relative acoustic attachment position"));
        }
        Ok(arm)
    }
}

impl RigidAcousticAttachment {
    /// Observe a migrating source or ideal microphone on a standalone carrier.
    /// Native motion supplies all carrier and Coriolis contributions exactly once.
    pub fn kinematics_relative(
        self,
        time_s: f64,
        motion: &PreparedRigidPointMotion,
        relative: RelativeAcousticMotion,
    ) -> Result<TimedSourceKinematics, RigidAcousticError> {
        if relative == RelativeAcousticMotion::default() {
            return self.kinematics(time_s, motion);
        }
        if !time_s.is_finite() {
            return Err(nonfinite("relative acoustic accepted timestamp"));
        }
        let point = motion.relative_point_motion(
            relative.arm(self)?, relative.velocity_body_m_s, relative.acceleration_body_m_s2,
        )?;
        Ok(TimedSourceKinematics {
            source_time_s: time_s,
            kinematics: SourceKinematics {
                position_m: components(point.kinematics.point_world),
                velocity_m_s: components(point.kinematics.point_velocity_world),
                acceleration_m_s2: components(point.point_acceleration_world),
            },
        })
    }

    /// Sample explicit migrating source strength in the carrier's rotating axes.
    /// The supplied force derivative includes any source change along the local
    /// path. Only carrier-axis rotation is added by the original transport law.
    pub fn sample_relative(
        self,
        time_s: f64,
        motion: &PreparedRigidPointMotion,
        relative: RelativeAcousticMotion,
        strength: BodyFrameCompactSourceStrength,
    ) -> Result<RigidAcousticSample, RigidAcousticError> {
        validate_strength(strength)?;
        let kinematics = self.kinematics_relative(time_s, motion, relative)?;
        sample_at_kinematics(kinematics, motion, strength)
    }

    /// Observe a point migrating relative to an articulated link's fixed anchor.
    /// Do not supply a child joint's motion again when querying that child link.
    pub fn kinematics_link_relative(
        self,
        time_s: f64,
        motion: &PreparedLinkPointMotion,
        relative: RelativeAcousticMotion,
    ) -> Result<TimedSourceKinematics, ArticulatedAcousticError> {
        if relative == RelativeAcousticMotion::default() {
            return self.kinematics_link(time_s, motion);
        }
        let point = self.relative_link_point(time_s, motion, relative)?;
        Ok(super::articulated::timed(time_s, point))
    }

    /// Sample a migrating link-frame source through the original force transport.
    pub fn sample_link_relative(
        self,
        time_s: f64,
        motion: &PreparedLinkPointMotion,
        relative: RelativeAcousticMotion,
        strength: BodyFrameCompactSourceStrength,
    ) -> Result<RigidAcousticSample, ArticulatedAcousticError> {
        if relative == RelativeAcousticMotion::default() {
            return self.sample_link(time_s, motion, strength);
        }
        validate_strength(strength)?;
        let point = self.relative_link_point(time_s, motion, relative)?;
        sample_at_point(time_s, point, motion, strength)
    }

    fn relative_link_point(
        self,
        time_s: f64,
        motion: &PreparedLinkPointMotion,
        relative: RelativeAcousticMotion,
    ) -> Result<fs_mbd::articulated::ArticulatedPointMotion, ArticulatedAcousticError> {
        if !time_s.is_finite() {
            return Err(nonfinite("relative link acoustic accepted timestamp").into());
        }
        Ok(motion.relative_point_motion(
            link_vector(relative.arm(self)?), link_vector(relative.velocity_body_m_s),
            link_vector(relative.acceleration_body_m_s2),
        )?)
    }
}

impl RigidAcousticHistory {
    /// Append smooth relative motion using the same paired, preallocated history.
    /// The retained attachment is the reference anchor, not the moving endpoint.
    pub fn append_relative(
        &mut self,
        time_s: f64,
        motion: &PreparedRigidPointMotion,
        relative: RelativeAcousticMotion,
        strength: BodyFrameCompactSourceStrength,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RigidAcousticError> {
        poll(cancelled)?;
        self.admit_append(time_s)?;
        let sample = self.attachment.sample_relative(time_s, motion, relative, strength)?;
        poll(cancelled)?;
        self.publish(sample);
        Ok(())
    }

    /// Append native link motion plus explicit local migration, atomically.
    /// Returned slices feed the existing finite-trajectory and pressure pipeline.
    pub fn append_link_relative(
        &mut self,
        time_s: f64,
        motion: &PreparedLinkPointMotion,
        relative: RelativeAcousticMotion,
        strength: BodyFrameCompactSourceStrength,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ArticulatedAcousticError> {
        poll(cancelled)?;
        self.admit_append(time_s)?;
        let sample = self.attachment.sample_link_relative(time_s, motion, relative, strength)?;
        poll(cancelled)?;
        self.publish(sample);
        Ok(())
    }
}

fn link_vector(v: Vec3) -> LinkVec3 { LinkVec3::new(v.x, v.y, v.z) }
