//! Admitted motion policies for the shared compact-source history renderer.

use super::{acoustic, admit_clearance, check_cancel, sizing};
use crate::render::RenderError;
use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::trajectory::PreparedSourceTrajectory;
use fs_aeroac::curle3d::moving::{
    MovingCompactObserver, PreparedMovingCompactSample, SourceKinematics, UniformSourceMotion,
};

pub(super) mod sealed {
    pub trait Sealed {}
}

/// Motion accepted by the shared physical-pressure history renderer.
///
/// Implemented for uniform motion, a borrowed prepared curved trajectory,
/// [`CompactSourceMotion`] for heterogeneous fixed-receiver ensembles, and a
/// borrowed prepared moving source/receiver pair. Sealing keeps
/// whole-domain speed/clearance guarantees with the original motion owners.
/// Dispatch is static for homogeneous histories and one enum match for mixed
/// histories. No trajectory is cloned, evolved, or allocated during callbacks.
pub trait CompactMotion: sealed::Sealed + Copy {
    /// The original receiver owner: a fixed microphone or a borrowed path.
    /// Every source in one render must admit this SAME receiver.
    type Observer: Copy;
    /// Verify the original microphone, material state and geometric limits.
    fn admit_observer(self, observer: &Self::Observer) -> Result<(), RenderError>;
    /// Sample the original receiver position at its physical arrival time.
    /// The returned pressure owner always uses source motion in the medium frame.
    fn observer_at(
        self,
        observer: &Self::Observer,
        receiver_time_s: f64,
    ) -> Result<MovingCompactObserver, RenderError>;
    /// Admit the complete requested clock and residual policy before any
    /// endpoint reconciliation. Fixed receivers have no separate clock owner.
    fn admit_clock(
        self,
        _first_time_s: f64,
        _last_time_s: f64,
        _maximum_time_residual_s: f64,
    ) -> Result<(), RenderError> {
        Ok(())
    }
    /// Maximum iterative root steps per receiver sample; uniform motion uses 0.
    fn maximum_root_iterations(self) -> usize;
    /// Evaluate position and its consistent derivatives on the source clock.
    fn kinematics_at(self, source_time_s: f64) -> Result<SourceKinematics, AeroacError>;
    /// Prepare the same physical pressure operator at one receiver timestamp.
    fn prepare_at<Cancel: FnMut() -> bool>(
        self,
        observer: &Self::Observer,
        receiver_time_s: f64,
        maximum_time_residual_s: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingCompactSample, RenderError>;
    /// Admit continuous geometry and return a whole-window Doppler-rate bound.
    fn admit_window<Cancel: FnMut() -> bool>(
        self,
        observer: &Self::Observer,
        first: &PreparedMovingCompactSample,
        last: &PreparedMovingCompactSample,
        cancelled: &mut Cancel,
    ) -> Result<f64, RenderError>;
}

impl sealed::Sealed for UniformSourceMotion {}
impl CompactMotion for UniformSourceMotion {
    type Observer = MovingCompactObserver;

    fn admit_observer(self, _observer: &MovingCompactObserver) -> Result<(), RenderError> {
        Ok(())
    }
    fn observer_at(
        self,
        observer: &MovingCompactObserver,
        _time: f64,
    ) -> Result<MovingCompactObserver, RenderError> {
        Ok(*observer)
    }
    fn maximum_root_iterations(self) -> usize {
        0
    }
    fn kinematics_at(self, time: f64) -> Result<SourceKinematics, AeroacError> {
        UniformSourceMotion::kinematics_at(&self, time)
    }
    fn prepare_at<Cancel: FnMut() -> bool>(
        self,
        observer: &MovingCompactObserver,
        time: f64,
        tolerance: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingCompactSample, RenderError> {
        check_cancel(cancelled)?;
        observer
            .prepare_uniform(time, &self, tolerance)
            .map_err(acoustic)
    }
    fn admit_window<Cancel: FnMut() -> bool>(
        self,
        observer: &MovingCompactObserver,
        first: &PreparedMovingCompactSample,
        last: &PreparedMovingCompactSample,
        cancelled: &mut Cancel,
    ) -> Result<f64, RenderError> {
        check_cancel(cancelled)?;
        admit_clearance(observer, &self, first.source_time_s(), last.source_time_s())?;
        // For uniform motion Mr decreases with source time, so the first
        // emission bounds this whole finite window, including receding paths.
        Ok(first.source_time_rate())
    }
}

impl sealed::Sealed for &PreparedSourceTrajectory {}
impl CompactMotion for &PreparedSourceTrajectory {
    type Observer = MovingCompactObserver;

    fn admit_observer(self, observer: &MovingCompactObserver) -> Result<(), RenderError> {
        if self.observer().medium() != observer.medium()
            || self.observer().observer_position_m() != observer.observer_position_m()
            || self.observer().limits() != observer.limits()
        {
            return Err(sizing(
                "prepared trajectory belongs to different observer, fluid or limits",
            ));
        }
        Ok(())
    }
    fn observer_at(
        self,
        observer: &MovingCompactObserver,
        _time: f64,
    ) -> Result<MovingCompactObserver, RenderError> {
        self.admit_observer(observer)?;
        Ok(*observer)
    }
    fn maximum_root_iterations(self) -> usize {
        PreparedSourceTrajectory::maximum_root_iterations(self)
    }
    fn kinematics_at(self, time: f64) -> Result<SourceKinematics, AeroacError> {
        PreparedSourceTrajectory::kinematics_at(self, time)
    }
    fn prepare_at<Cancel: FnMut() -> bool>(
        self,
        observer: &MovingCompactObserver,
        time: f64,
        tolerance: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingCompactSample, RenderError> {
        self.admit_observer(observer)?;
        self.prepare_receiver(time, tolerance, cancelled)
            .map_err(acoustic)
    }
    fn admit_window<Cancel: FnMut() -> bool>(
        self,
        observer: &MovingCompactObserver,
        first: &PreparedMovingCompactSample,
        last: &PreparedMovingCompactSample,
        cancelled: &mut Cancel,
    ) -> Result<f64, RenderError> {
        check_cancel(cancelled)?;
        self.admit_observer(observer)?;
        let [start, end] = self.source_time_interval_s();
        if first.source_time_s() < start
            || last.source_time_s() > end
            || first.source_time_s() > last.source_time_s()
        {
            return Err(sizing(
                "receiver window exceeds the retained curved trajectory",
            ));
        }
        // Curved paths can turn toward the observer later in the window.
        // Their original whole-domain speed bound must replace the uniform
        // path's first-sample radial-Mach shortcut.
        Ok(self.maximum_source_time_rate())
    }
}

/// Mix analytic uniform paths and prepared curved paths at the same microphone.
/// All curved paths retain their original finite trajectory and observer owner.
#[derive(Clone, Copy)]
pub enum CompactSourceMotion<'a> {
    /// Globally uniform source translation with an analytic retarded root.
    Uniform(UniformSourceMotion),
    /// Original C2 trajectory, speed/clearance bounds and bounded root policy.
    Trajectory(&'a PreparedSourceTrajectory),
}

impl core::fmt::Debug for CompactSourceMotion<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Uniform(motion) => f.debug_tuple("Uniform").field(motion).finish(),
            Self::Trajectory(path) => f
                .debug_struct("Trajectory")
                .field("segments", &path.segment_count())
                .field("source_times_s", &path.source_time_interval_s())
                .finish(),
        }
    }
}

impl sealed::Sealed for CompactSourceMotion<'_> {}
impl CompactMotion for CompactSourceMotion<'_> {
    type Observer = MovingCompactObserver;

    fn admit_observer(self, observer: &MovingCompactObserver) -> Result<(), RenderError> {
        match self {
            Self::Uniform(m) => m.admit_observer(observer),
            Self::Trajectory(m) => m.admit_observer(observer),
        }
    }
    fn observer_at(
        self,
        observer: &MovingCompactObserver,
        time: f64,
    ) -> Result<MovingCompactObserver, RenderError> {
        match self {
            Self::Uniform(m) => m.observer_at(observer, time),
            Self::Trajectory(m) => m.observer_at(observer, time),
        }
    }
    fn maximum_root_iterations(self) -> usize {
        match self {
            Self::Uniform(m) => CompactMotion::maximum_root_iterations(m),
            Self::Trajectory(m) => m.maximum_root_iterations(),
        }
    }
    fn kinematics_at(self, time: f64) -> Result<SourceKinematics, AeroacError> {
        match self {
            Self::Uniform(m) => m.kinematics_at(time),
            Self::Trajectory(m) => m.kinematics_at(time),
        }
    }
    fn prepare_at<Cancel: FnMut() -> bool>(
        self,
        observer: &MovingCompactObserver,
        time: f64,
        tolerance: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingCompactSample, RenderError> {
        match self {
            Self::Uniform(m) => m.prepare_at(observer, time, tolerance, cancelled),
            Self::Trajectory(m) => m.prepare_at(observer, time, tolerance, cancelled),
        }
    }
    fn admit_window<Cancel: FnMut() -> bool>(
        self,
        observer: &MovingCompactObserver,
        first: &PreparedMovingCompactSample,
        last: &PreparedMovingCompactSample,
        cancelled: &mut Cancel,
    ) -> Result<f64, RenderError> {
        match self {
            Self::Uniform(m) => m.admit_window(observer, first, last, cancelled),
            Self::Trajectory(m) => m.admit_window(observer, first, last, cancelled),
        }
    }
}
