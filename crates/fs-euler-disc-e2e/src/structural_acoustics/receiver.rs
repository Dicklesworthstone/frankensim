//! Identified receiver motion on the original native acoustic trajectory owner.

use fs_aeroac::curle3d::moving::trajectory::{
    PreparedKinematicTrajectory, TimedSourceKinematics, TrajectoryBudget,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_material::fluid_acoustics::FluidAcousticState;

use super::retarded_clock::clip_reception_window;
use super::{
    AcousticWorldObserver, PhysicalPressureObserver, RetardedFarFieldObserverControls,
    StructuralModalBasisError, arrival_time, observer_vector, retarded_clock_error,
};
use crate::timeline_resampling::TimelineResampler;

/// Immutable finite microphone motion in the same quiescent world frame as
/// the Euler source. The native C2 quintics and continuous subsonic admission
/// are reused; accepted mechanics and source coefficients are never changed.
#[derive(Debug)]
pub struct AcousticReceiverTrajectory {
    path: PreparedKinematicTrajectory,
    identity: ContentHash,
}

impl AcousticReceiverTrajectory {
    /// Admit a finite receiver history with explicit storage, work and Mach caps.
    ///
    /// The native knot field `source_time_s` denotes the receiver's physical
    /// event time here. It is not emission time or a propagation-shifted clock.
    /// Acceleration and velocity are derivatives in the medium's rest frame.
    ///
    /// # Errors
    /// Refuses invalid knots, non-subsonic interpolants, exhausted budgets or
    /// cancellation through the original trajectory owner.
    pub fn new(
        medium: FluidAcousticState,
        knots: &[TimedSourceKinematics],
        maximum_mach: f64,
        budget: TrajectoryBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, StructuralModalBasisError> {
        let path = PreparedKinematicTrajectory::new(
            medium,
            knots,
            maximum_mach,
            budget,
            &mut || cx.checkpoint().is_err(),
        )?;
        let mut identity =
            DomainHasher::new("org.frankensim.euler-disc.acoustic-receiver-trajectory.v1");
        identity.update(b"native-c2-quintic-quiescent-world");
        for value in [medium.density_kg_m3(), medium.sound_speed_m_s(), maximum_mach] {
            identity.update(&value.to_bits().to_le_bytes());
        }
        identity.update(&(knots.len() as u64).to_le_bytes());
        for (index, knot) in knots.iter().enumerate() {
            if index % 64 == 0 {
                cx.checkpoint().map_err(|_| StructuralModalBasisError::Cancelled)?;
            }
            identity.update(&knot.source_time_s.to_bits().to_le_bytes());
            for value in knot.kinematics.position_m.into_iter()
                .chain(knot.kinematics.velocity_m_s)
                .chain(knot.kinematics.acceleration_m_s2)
            {
                identity.update(&value.to_bits().to_le_bytes());
            }
        }
        cx.checkpoint().map_err(|_| StructuralModalBasisError::Cancelled)?;
        Ok(Self { path, identity: identity.finalize() })
    }

    /// Identity of the original knots, medium and declared motion admission.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Exact finite reception-time domain; no extrapolation is permitted.
    #[must_use]
    pub fn time_interval_s(&self) -> [f64; 2] {
        self.path.time_interval_s()
    }

    /// Whole-path outward speed bound divided by the original sound speed.
    #[must_use]
    pub fn maximum_mach(&self) -> f64 {
        self.path.maximum_mach()
    }
}

/// One receiver of the shared retarded structural/rigid-disc field.
///
/// A moving receiver is evaluated at reception time, once per output query.
/// That position remains fixed while the source emission time is solved.
#[derive(Clone, Copy, Debug)]
pub enum AcousticWorldReceiver<'a> {
    /// Stationary microphone with no finite receiver-history restriction.
    Fixed(AcousticWorldObserver),
    /// Immutable admitted microphone history, sampled without extrapolation.
    Moving(&'a AcousticReceiverTrajectory),
}

impl AcousticWorldReceiver<'_> {
    pub(super) fn validate(self, sound_speed_m_s: f64) -> Result<(), StructuralModalBasisError> {
        match self {
            Self::Fixed(observer) => {
                if observer.position_world_m.iter().any(|v| !v.is_finite()) {
                    return Err(StructuralModalBasisError::InvalidRequest {
                        what: "retarded receiver position must be finite",
                    });
                }
            }
            Self::Moving(receiver) => {
                if receiver.path.medium().sound_speed_m_s().to_bits() != sound_speed_m_s.to_bits() {
                    return Err(StructuralModalBasisError::IdentityMismatch {
                        what: "receiver motion and source radiation must use the same sound speed",
                    });
                }
            }
        }
        Ok(())
    }

    pub(super) fn maximum_mach(self) -> f64 {
        match self {
            Self::Fixed(_) => 0.0,
            Self::Moving(receiver) => receiver.maximum_mach(),
        }
    }

    pub(super) fn observer(self) -> PhysicalPressureObserver {
        match self {
            Self::Fixed(observer) => PhysicalPressureObserver::WorldFixed(observer),
            Self::Moving(receiver) => PhysicalPressureObserver::WorldMoving {
                trajectory_identity: receiver.identity(),
            },
        }
    }

    pub(super) fn position_at(
        self,
        reception_s: f64,
    ) -> Result<AcousticWorldObserver, StructuralModalBasisError> {
        match self {
            Self::Fixed(observer) => Ok(observer),
            Self::Moving(receiver) => Ok(AcousticWorldObserver {
                position_world_m: receiver.path.kinematics_at(reception_s)?.position_m,
            }),
        }
    }

    pub(super) fn reception_window(
        self,
        timeline: &TimelineResampler<'_>,
        emission_interval_s: [f64; 2],
        sound_speed_m_s: f64,
        controls: RetardedFarFieldObserverControls,
        cx: &Cx<'_>,
    ) -> Result<[f64; 2], StructuralModalBasisError> {
        match self {
            Self::Fixed(observer) => Ok([
                arrival_time(timeline, observer, emission_interval_s[0], sound_speed_m_s)?,
                arrival_time(timeline, observer, emission_interval_s[1], sound_speed_m_s)?,
            ]),
            Self::Moving(receiver) => {
                let source_positions = [
                    timeline.position(emission_interval_s[0])?,
                    timeline.position(emission_interval_s[1])?,
                ];
                clip_reception_window(
                    receiver.time_interval_s(),
                    controls.root_time_tolerance_s,
                    controls.maximum_root_iterations,
                    |endpoint, reception_s| {
                        cx.checkpoint().map_err(|_| StructuralModalBasisError::Cancelled)?;
                        let observer = self.position_at(reception_s)?;
                        let (_, radius) = observer_vector(observer, source_positions[endpoint]);
                        let arrival = emission_interval_s[endpoint] + radius / sound_speed_m_s;
                        // Zero distance is valid for an ENDPOINT wavefront
                        // residual at a different reception time. The actual
                        // solved emission pair still must pass the far-field gate.
                        if !(radius >= 0.0 && radius.is_finite() && arrival.is_finite()) {
                            return Err(StructuralModalBasisError::InvalidRequest {
                                what: "receiver history has an invalid wavefront separation",
                            });
                        }
                        // This is exactly the negative of the source-root
                        // residual, including its floating-point operation order.
                        Ok(reception_s - arrival)
                    },
                ).map_err(retarded_clock_error)
            }
        }
    }
}
