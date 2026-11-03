//! Native rigid motion and body-frame sources on the acoustic emission clock.
//!
//! The original `fs-mbd` owner supplies point position, velocity and acceleration
//! from accepted momenta, full inertia, instantaneous wrench and its gravity.
//! Its world frame must be the acoustic medium's rest frame. Observation does
//! not advance mechanics, estimate a missing load, or turn total body force
//! into an acoustic source. Source force here is explicitly force ON THE FLUID.
//!
//! A rotating force needs the transport term:
//! `F_world = R F_body`, `Fdot_world = R(Fdot_body + omega_body cross F_body)`.
//! `Fdot_body` is the derivative of body-coordinate components, not a world
//! derivative merely expressed in body axes. Volume flow and its derivative
//! are scalars. The point's centripetal/tangential acceleration stays with the
//! trajectory, independently of this force-component rotation.
//!
//! The history collector records finite smooth accepted endpoints. It neither
//! supplies unresolved impacts nor makes a C2 interpolant valid across an
//! impulse, detachment or change of mechanical model. Source bandwidth,
//! finite-source compactness, interpolation convergence and radiation reaction
//! retain their original owners and admission requirements.
//! The declared strength bandwidth must cover the exported medium-frame force,
//! including rotation: a constant body-frame force can still radiate periodically.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::trajectory::TimedSourceKinematics;
use fs_aeroac::curle3d::moving::{CompactSourceStrength, SourceKinematics, TimedCompactSource};
use fs_mbd::{DynamicsError, PreparedRigidPointMotion, Vec3};

/// A source law evaluated in the original body's instantaneous coordinate frame.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BodyFrameCompactSourceStrength {
    /// Outward volume injection [m^3/s], not mass injection or body velocity.
    pub volume_flow_m3_s: f64,
    /// Physical time derivative of volume injection [m^3/s^2].
    pub volume_flow_rate_m3_s2: f64,
    /// Total force ON the fluid [N], in body axes, not force per span.
    pub force_on_fluid_body_n: [f64; 3],
    /// Time derivative of the BODY COMPONENTS of that force [N/s].
    pub force_rate_on_fluid_body_n_s: [f64; 3],
}

/// Geometry and source strengths evaluated at exactly one accepted timestamp.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidAcousticSample {
    /// Native material-point kinematics on the physical emission clock.
    pub kinematics: TimedSourceKinematics,
    /// Fixed-medium-frame source values and derivatives at that SAME time.
    pub source: TimedCompactSource,
}

/// A material point fixed in the original centre-of-mass body coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidAcousticAttachment {
    arm_body_m: Vec3,
}

impl RigidAcousticAttachment {
    /// Admit the actual acoustic source or microphone position in body axes.
    /// No geometric radius, effective area or radiation law is inferred.
    pub fn new(arm_body_m: Vec3) -> Result<Self, RigidAcousticError> {
        if !arm_body_m.is_finite() {
            return Err(nonfinite("rigid acoustic body attachment"));
        }
        Ok(Self { arm_body_m })
    }

    /// Original centre-of-mass-relative attachment [m].
    pub const fn arm_body_m(self) -> Vec3 {
        self.arm_body_m
    }

    /// Observe a native point without inventing an acoustic source strength.
    /// This same path supplies an attached ideal microphone's accepted motion.
    pub fn kinematics(
        self,
        time_s: f64,
        motion: &PreparedRigidPointMotion,
    ) -> Result<TimedSourceKinematics, RigidAcousticError> {
        if !time_s.is_finite() {
            return Err(nonfinite("rigid acoustic accepted timestamp"));
        }
        let point = motion.point_motion(self.arm_body_m)?;
        Ok(TimedSourceKinematics {
            source_time_s: time_s,
            kinematics: SourceKinematics {
                position_m: components(point.kinematics.point_world),
                velocity_m_s: components(point.kinematics.point_velocity_world),
                acceleration_m_s2: components(point.point_acceleration_world),
            },
        })
    }

    /// Transport one explicitly supplied body-frame source into the medium frame.
    ///
    /// Mechanical wrench and acoustic strength remain separate inputs: gravity,
    /// support/contact loads and a net body torque are not automatically fluid
    /// dipoles. If mechanics reports fluid force ON THE BODY, its source owner
    /// must select the physical acoustic contribution and reverse its sign.
    pub fn sample(
        self,
        time_s: f64,
        motion: &PreparedRigidPointMotion,
        strength: BodyFrameCompactSourceStrength,
    ) -> Result<RigidAcousticSample, RigidAcousticError> {
        if !strength.volume_flow_m3_s.is_finite()
            || !strength.volume_flow_rate_m3_s2.is_finite()
            || strength
                .force_on_fluid_body_n
                .iter()
                .chain(&strength.force_rate_on_fluid_body_n_s)
                .any(|value| !value.is_finite())
        {
            return Err(nonfinite("rigid acoustic body-frame strength"));
        }
        let kinematics = self.kinematics(time_s, motion)?;
        let force_body = vector(strength.force_on_fluid_body_n);
        let derivative_body = vector(strength.force_rate_on_fluid_body_n_s)
            .add(motion.angular_velocity_body().cross(force_body));
        if !derivative_body.is_finite() {
            return Err(nonfinite("rigid acoustic force transport derivative"));
        }
        let orientation = motion.pose().orientation();
        let force_world = orientation.rotate_body_to_world(force_body);
        let derivative_world = orientation.rotate_body_to_world(derivative_body);
        if !force_world.is_finite() || !derivative_world.is_finite() {
            return Err(nonfinite("rigid acoustic medium-frame strength"));
        }
        Ok(RigidAcousticSample {
            kinematics,
            source: TimedCompactSource {
                source_time_s: time_s,
                strength: CompactSourceStrength {
                    volume_flow_m3_s: strength.volume_flow_m3_s,
                    volume_flow_rate_m3_s2: strength.volume_flow_rate_m3_s2,
                    force_on_fluid_n: components(force_world),
                    force_rate_on_fluid_n_s: components(derivative_world),
                },
            },
        })
    }
}

/// Explicit retained endpoint and logical payload limits for a finite history.
#[derive(Clone, Copy, Debug)]
pub struct RigidHistoryBudget {
    /// Maximum accepted endpoints, at least two for trajectory interpolation.
    pub maximum_samples: usize,
    /// This owner plus both full-capacity endpoint arrays. Allocator metadata
    /// and allocator rounding are excluded; there is no append-time allocation.
    pub maximum_payload_bytes: usize,
}

/// A bounded source/trajectory history from the same accepted mechanical events.
///
/// Source and microphone trajectories may share each native prepared snapshot.
/// This owner keeps one fixed body attachment and records source strengths as
/// well as its motion. Callers retain the mechanical identity, physical clock
/// and smoothness of this sequence; the native state contains no time field.
/// Both borrowed slices are immutable during observation, enforced by Rust.
pub struct RigidAcousticHistory {
    attachment: RigidAcousticAttachment,
    kinematics: Vec<TimedSourceKinematics>,
    sources: Vec<TimedCompactSource>,
    maximum_samples: usize,
}

impl RigidAcousticHistory {
    /// Reserve the complete history once, before observing mechanical events.
    pub fn new(
        attachment: RigidAcousticAttachment,
        budget: RigidHistoryBudget,
    ) -> Result<Self, RigidAcousticError> {
        let bytes = size_of::<TimedSourceKinematics>()
            .checked_add(size_of::<TimedCompactSource>())
            .and_then(|bytes| bytes.checked_mul(budget.maximum_samples))
            .and_then(|bytes| bytes.checked_add(size_of::<Self>()))
            .ok_or_else(|| invalid("rigid acoustic history payload overflow"))?;
        if budget.maximum_samples < 2 || bytes > budget.maximum_payload_bytes {
            return Err(invalid("rigid acoustic history endpoint/payload budget"));
        }
        let mut kinematics = Vec::new();
        let mut sources = Vec::new();
        kinematics
            .try_reserve_exact(budget.maximum_samples)
            .map_err(|_| invalid("cannot reserve rigid acoustic motion history"))?;
        sources
            .try_reserve_exact(budget.maximum_samples)
            .map_err(|_| invalid("cannot reserve rigid acoustic strength history"))?;
        Ok(Self {
            attachment,
            kinematics,
            sources,
            maximum_samples: budget.maximum_samples,
        })
    }

    /// Stage one complete accepted event; publish both slices only on success.
    /// Failure/cancellation changes neither history nor any mechanical state.
    pub fn append(
        &mut self,
        time_s: f64,
        motion: &PreparedRigidPointMotion,
        strength: BodyFrameCompactSourceStrength,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RigidAcousticError> {
        poll(cancelled)?;
        if self.sources.len() >= self.maximum_samples {
            return Err(invalid(
                "rigid acoustic history endpoint capacity exhausted",
            ));
        }
        if !time_s.is_finite()
            || self
                .sources
                .last()
                .is_some_and(|last| time_s <= last.source_time_s)
        {
            return Err(invalid(
                "rigid acoustic times must be finite and strictly increasing",
            ));
        }
        let sample = self.attachment.sample(time_s, motion, strength)?;
        poll(cancelled)?;
        self.kinematics.push(sample.kinematics);
        self.sources.push(sample.source);
        Ok(())
    }

    /// Accepted trajectory endpoints for the original C2 path preparers.
    pub fn kinematics(&self) -> &[TimedSourceKinematics] {
        &self.kinematics
    }
    /// Accepted fixed-frame source strengths for the shared pressure renderer.
    pub fn sources(&self) -> &[TimedCompactSource] {
        &self.sources
    }
    /// The fixed original body attachment retained by this history.
    pub const fn attachment(&self) -> RigidAcousticAttachment {
        self.attachment
    }
    /// Number of complete paired events.
    pub fn len(&self) -> usize {
        self.sources.len()
    }
    /// Whether no mechanical event has yet been observed.
    pub fn is_empty(&self) -> bool {
        self.sources.is_empty()
    }
    /// Further events that fit the original retained payload.
    pub fn remaining_capacity(&self) -> usize {
        self.maximum_samples - self.sources.len()
    }
}

/// Preserve native dynamics refusals separately from acoustic input admission.
#[derive(Debug, Clone, PartialEq)]
pub enum RigidAcousticError {
    /// The original rigid-motion owner rejected an input or derived quantity.
    Dynamics(DynamicsError),
    /// Source, clock, history budget or cancellation refusal.
    Acoustic(AeroacError),
}

impl From<DynamicsError> for RigidAcousticError {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}
impl core::fmt::Display for RigidAcousticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Dynamics(error) => write!(f, "rigid acoustic motion: {error}"),
            Self::Acoustic(error) => write!(f, "rigid acoustic observation: {error}"),
        }
    }
}
impl std::error::Error for RigidAcousticError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Dynamics(error) => Some(error),
            Self::Acoustic(error) => Some(error),
        }
    }
}

fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn invalid(what: &'static str) -> RigidAcousticError {
    RigidAcousticError::Acoustic(AeroacError::InvalidParameter { what })
}
fn nonfinite(what: &'static str) -> RigidAcousticError {
    RigidAcousticError::Acoustic(AeroacError::NonFinite { what })
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), RigidAcousticError> {
    if cancelled() {
        Err(RigidAcousticError::Acoustic(AeroacError::Cancelled {
            what: "accepted rigid acoustic history",
        }))
    } else {
        Ok(())
    }
}
