//! Reuse admitted compact-source motion without duplicating its clock solver.

use super::check_cancel;
use crate::pcm_wav::observation::moving::{CompactMotion, CompactSourceMotion};
use crate::render::RenderError;
use fs_aeroac::curle3d::moving::trajectory::PreparedSourceTrajectory;
use fs_aeroac::curle3d::moving::{MovingCompactObserver, UniformSourceMotion};
use fs_aeroac::curle3d::quadrupole::transient::moving::PreparedMovingQuadrupoleSample;

/// Fixed-receiver direct paths supplying complete stress-potential derivatives.
///
/// The sealed original `CompactMotion` owners supply continuous admission and
/// Doppler bounds. This subset additionally supplies consistent acceleration
/// and jerk for the quadrupole Hessian. Uniform paths retain their analytic
/// inverse; borrowed quintics retain their original safeguarded inverse and
/// work cap. Mixed ensembles use the existing `CompactSourceMotion` enum.
pub trait QuadrupoleMotion: CompactMotion<Observer = MovingCompactObserver> {
    /// Prepare the original stress-potential Hessian at one receiver event.
    fn prepare_quadrupole_at<Cancel: FnMut() -> bool>(
        self,
        observer: &MovingCompactObserver,
        receiver_time_s: f64,
        maximum_time_residual_s: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingQuadrupoleSample, RenderError>;
}

impl QuadrupoleMotion for UniformSourceMotion {
    fn prepare_quadrupole_at<Cancel: FnMut() -> bool>(
        self,
        observer: &MovingCompactObserver,
        time: f64,
        tolerance: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingQuadrupoleSample, RenderError> {
        check_cancel(cancelled)?;
        observer
            .prepare_uniform_quadrupole(time, &self, tolerance)
            .map_err(RenderError::Acoustic)
    }
}

impl QuadrupoleMotion for &PreparedSourceTrajectory {
    fn prepare_quadrupole_at<Cancel: FnMut() -> bool>(
        self,
        observer: &MovingCompactObserver,
        time: f64,
        tolerance: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingQuadrupoleSample, RenderError> {
        self.admit_observer(observer)?;
        self.prepare_quadrupole_receiver(time, tolerance, cancelled)
            .map_err(RenderError::Acoustic)
    }
}

impl QuadrupoleMotion for CompactSourceMotion<'_> {
    fn prepare_quadrupole_at<Cancel: FnMut() -> bool>(
        self,
        observer: &MovingCompactObserver,
        time: f64,
        tolerance: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingQuadrupoleSample, RenderError> {
        match self {
            Self::Uniform(motion) => {
                motion.prepare_quadrupole_at(observer, time, tolerance, cancelled)
            }
            Self::Trajectory(path) => {
                path.prepare_quadrupole_at(observer, time, tolerance, cancelled)
            }
        }
    }
}
