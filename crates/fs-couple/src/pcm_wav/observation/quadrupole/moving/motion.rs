//! Reuse admitted compact-source motion without duplicating its clock solver.

use super::{check_cancel, sizing};
use crate::pcm_wav::observation::moving::{CompactMotion, CompactSourceMotion};
use crate::render::RenderError;
use fs_aeroac::curle3d::moving::trajectory::{
    PreparedMovingReceiverSource, PreparedReflectedMovingReceiverSource,
    PreparedReflectedSourceTrajectory, PreparedSourceTrajectory,
};
use fs_aeroac::curle3d::moving::{MovingCompactObserver, UniformSourceMotion};
use fs_aeroac::curle3d::quadrupole::transient::moving::{
    PreparedMovingQuadrupoleFieldSample, PreparedMovingQuadrupoleSample,
};

/// Original direct/image paths supplying complete stress-potential derivatives.
///
/// The sealed original `CompactMotion` owners supply continuous admission and
/// Doppler bounds. This subset additionally supplies consistent acceleration
/// and jerk for the quadrupole Hessian. Uniform paths retain their analytic
/// inverse; borrowed quintics retain their original safeguarded inverse and
/// work cap. Native moving microphones and ideal boundaries keep their own
/// receiver identity, complete path group and whole-window Doppler admission.
/// Mixed fixed-receiver ensembles use the existing `CompactSourceMotion` enum.
pub trait QuadrupoleMotion: CompactMotion {
    /// Complete local pressure/velocity operator for one original image path.
    /// Its explicit history argument uses the corresponding virtual receiver;
    /// the prepared field applies image parity and vector transformation once.
    fn prepare_quadrupole_field_path_at<Cancel: FnMut() -> bool>(
        self,
        path: usize,
        observer: &Self::Observer,
        receiver_time_s: f64,
        maximum_time_residual_s: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingQuadrupoleFieldSample, RenderError>;

    /// Prepare the original stress-potential Hessian at one receiver event.
    fn prepare_quadrupole_at<Cancel: FnMut() -> bool>(
        self,
        observer: &Self::Observer,
        receiver_time_s: f64,
        maximum_time_residual_s: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingQuadrupoleSample, RenderError>;

    /// One original path, with image parity already in its pressure operator.
    /// The source reader still receives the PHYSICAL source index and clock.
    fn prepare_quadrupole_path_at<Cancel: FnMut() -> bool>(
        self,
        path: usize,
        observer: &Self::Observer,
        receiver_time_s: f64,
        maximum_time_residual_s: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingQuadrupoleSample, RenderError> {
        if path != 0 {
            return Err(sizing(
                "free-space quadrupole motion has only its direct path",
            ));
        }
        self.prepare_quadrupole_at(
            observer,
            receiver_time_s,
            maximum_time_residual_s,
            cancelled,
        )
    }
}

impl QuadrupoleMotion for UniformSourceMotion {
    fn prepare_quadrupole_field_path_at<Cancel: FnMut() -> bool>(
        self,
        path: usize,
        observer: &MovingCompactObserver,
        time: f64,
        tolerance: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingQuadrupoleFieldSample, RenderError> {
        check_cancel(cancelled)?;
        if path != 0 {
            return Err(sizing("uniform quadrupole has only its direct path"));
        }
        observer
            .prepare_uniform_quadrupole_field(time, &self, tolerance)
            .map_err(RenderError::Acoustic)
    }

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
    fn prepare_quadrupole_field_path_at<Cancel: FnMut() -> bool>(
        self,
        path: usize,
        observer: &MovingCompactObserver,
        time: f64,
        tolerance: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingQuadrupoleFieldSample, RenderError> {
        check_cancel(cancelled)?;
        if path != 0 {
            return Err(sizing("free-space quadrupole has only its direct path"));
        }
        self.admit_observer(observer)?;
        self.prepare_quadrupole_field_receiver(time, tolerance, cancelled)
            .map_err(RenderError::Acoustic)
    }

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
    fn prepare_quadrupole_field_path_at<Cancel: FnMut() -> bool>(
        self,
        path: usize,
        observer: &MovingCompactObserver,
        time: f64,
        tolerance: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingQuadrupoleFieldSample, RenderError> {
        match self {
            Self::Uniform(motion) => {
                motion.prepare_quadrupole_field_path_at(path, observer, time, tolerance, cancelled)
            }
            Self::Trajectory(motion) => {
                motion.prepare_quadrupole_field_path_at(path, observer, time, tolerance, cancelled)
            }
        }
    }

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

impl QuadrupoleMotion for &PreparedReflectedSourceTrajectory<'_> {
    fn prepare_quadrupole_field_path_at<Cancel: FnMut() -> bool>(
        self,
        path: usize,
        observer: &Self::Observer,
        time: f64,
        tolerance: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingQuadrupoleFieldSample, RenderError> {
        check_cancel(cancelled)?;
        self.admit_observer(observer)?;
        self.admit_clock(time, time, tolerance)?;
        self.prepare_quadrupole_field_path(path, time, tolerance, cancelled)
            .map_err(RenderError::Acoustic)
    }

    fn prepare_quadrupole_at<Cancel: FnMut() -> bool>(
        self,
        observer: &Self::Observer,
        time: f64,
        tolerance: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingQuadrupoleSample, RenderError> {
        self.prepare_quadrupole_path_at(0, observer, time, tolerance, cancelled)
    }

    fn prepare_quadrupole_path_at<Cancel: FnMut() -> bool>(
        self,
        path: usize,
        observer: &Self::Observer,
        time: f64,
        tolerance: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingQuadrupoleSample, RenderError> {
        check_cancel(cancelled)?;
        self.admit_observer(observer)?;
        self.admit_clock(time, time, tolerance)?;
        self.prepare_quadrupole_path(path, time, tolerance, cancelled)
            .map_err(RenderError::Acoustic)
    }
}

impl QuadrupoleMotion for &PreparedMovingReceiverSource<'_> {
    fn prepare_quadrupole_field_path_at<Cancel: FnMut() -> bool>(
        self,
        path: usize,
        observer: &Self::Observer,
        time: f64,
        tolerance: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingQuadrupoleFieldSample, RenderError> {
        check_cancel(cancelled)?;
        if path != 0 {
            return Err(sizing("free-space quadrupole has only its direct path"));
        }
        self.admit_observer(observer)?;
        self.admit_clock(time, time, tolerance)?;
        self.prepare_quadrupole_field_receiver(time, cancelled)
            .map_err(RenderError::Acoustic)
    }

    fn prepare_quadrupole_at<Cancel: FnMut() -> bool>(
        self,
        observer: &Self::Observer,
        time: f64,
        tolerance: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingQuadrupoleSample, RenderError> {
        check_cancel(cancelled)?;
        self.admit_observer(observer)?;
        self.admit_clock(time, time, tolerance)?;
        self.prepare_quadrupole_receiver(time, cancelled)
            .map_err(RenderError::Acoustic)
    }
}

impl QuadrupoleMotion for &PreparedReflectedMovingReceiverSource<'_> {
    fn prepare_quadrupole_field_path_at<Cancel: FnMut() -> bool>(
        self,
        path: usize,
        observer: &Self::Observer,
        time: f64,
        tolerance: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingQuadrupoleFieldSample, RenderError> {
        check_cancel(cancelled)?;
        self.admit_observer(observer)?;
        self.admit_clock(time, time, tolerance)?;
        self.prepare_quadrupole_field_path(path, time, cancelled)
            .map_err(RenderError::Acoustic)
    }

    fn prepare_quadrupole_at<Cancel: FnMut() -> bool>(
        self,
        observer: &Self::Observer,
        time: f64,
        tolerance: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingQuadrupoleSample, RenderError> {
        self.prepare_quadrupole_path_at(0, observer, time, tolerance, cancelled)
    }

    fn prepare_quadrupole_path_at<Cancel: FnMut() -> bool>(
        self,
        path: usize,
        observer: &Self::Observer,
        time: f64,
        tolerance: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingQuadrupoleSample, RenderError> {
        check_cancel(cancelled)?;
        self.admit_observer(observer)?;
        self.admit_clock(time, time, tolerance)?;
        self.prepare_quadrupole_path(path, time, cancelled)
            .map_err(RenderError::Acoustic)
    }
}
