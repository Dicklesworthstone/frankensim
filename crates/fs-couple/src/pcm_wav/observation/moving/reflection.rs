//! Ideal reflected paths through the original accepted-history renderer.

use super::motion::{CompactMotion, sealed};
use super::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
    acoustic, check_cancel, sizing,
};
use crate::render::RenderError;
use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::trajectory::{
    IdealReflectionObserver, PreparedReflectedSourceTrajectory,
};
use fs_aeroac::curle3d::moving::{
    MovingCompactObserver, PreparedMovingCompactSample, SourceKinematics,
};

impl sealed::Sealed for &PreparedReflectedSourceTrajectory<'_> {}

impl<'path> CompactMotion for &PreparedReflectedSourceTrajectory<'path> {
    type Observer = &'path IdealReflectionObserver;

    fn admit_observer(self, observer: &Self::Observer) -> Result<(), RenderError> {
        if !core::ptr::eq(self.reflection_observer(), *observer) {
            return Err(sizing(
                "reflected sources must share the original microphone and boundaries",
            ));
        }
        Ok(())
    }

    fn observer_at(
        self,
        observer: &Self::Observer,
        time: f64,
    ) -> Result<MovingCompactObserver, RenderError> {
        self.admit_observer(observer)?;
        self.admit_clock(time, time, 0.0)?;
        Ok(*observer.observer())
    }

    fn admit_clock(self, first: f64, last: f64, tolerance: f64) -> Result<(), RenderError> {
        let [start, end] = self.receiver_time_interval_s();
        if !first.is_finite()
            || !last.is_finite()
            || first < start
            || last > end
            || first > last
            || !tolerance.is_finite()
            || tolerance < 0.0
        {
            return Err(sizing(
                "pressure clock lacks complete direct and reflected source history",
            ));
        }
        Ok(())
    }

    fn maximum_root_iterations(self) -> usize {
        PreparedReflectedSourceTrajectory::maximum_root_iterations(self)
    }
    fn path_count(self) -> usize {
        PreparedReflectedSourceTrajectory::path_count(self)
    }
    fn kinematics_at(self, time: f64) -> Result<SourceKinematics, AeroacError> {
        self.trajectory().kinematics_at(time)
    }

    fn prepare_at<Cancel: FnMut() -> bool>(
        self,
        observer: &Self::Observer,
        time: f64,
        tolerance: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingCompactSample, RenderError> {
        self.prepare_path_at(0, observer, time, tolerance, cancelled)
    }

    fn prepare_path_at<Cancel: FnMut() -> bool>(
        self,
        path: usize,
        observer: &Self::Observer,
        time: f64,
        tolerance: f64,
        cancelled: &mut Cancel,
    ) -> Result<PreparedMovingCompactSample, RenderError> {
        check_cancel(cancelled)?;
        self.admit_observer(observer)?;
        self.admit_clock(time, time, tolerance)?;
        self.prepare_path(path, time, tolerance, cancelled)
            .map_err(acoustic)
    }

    fn prepare_path_endpoint(
        self,
        path: usize,
        observer: &Self::Observer,
        source_time: f64,
        receiver_time: f64,
    ) -> Result<PreparedMovingCompactSample, RenderError> {
        self.admit_observer(observer)?;
        self.admit_clock(receiver_time, receiver_time, 0.0)?;
        self.prepare_path_emission(path, source_time)
            .map_err(acoustic)
    }

    fn admit_window<Cancel: FnMut() -> bool>(
        self,
        observer: &Self::Observer,
        first: &PreparedMovingCompactSample,
        last: &PreparedMovingCompactSample,
        cancelled: &mut Cancel,
    ) -> Result<f64, RenderError> {
        check_cancel(cancelled)?;
        self.admit_observer(observer)?;
        self.admit_clock(first.receiver_time_s(), last.receiver_time_s(), 0.0)?;
        let [start, end] = self.trajectory().source_time_interval_s();
        if first.source_time_s() < start
            || last.source_time_s() > end
            || first.source_time_s() > last.source_time_s()
        {
            return Err(sizing(
                "reflected source window exceeds the original trajectory",
            ));
        }
        Ok(self.maximum_source_time_rate())
    }
}

impl<'history, 'motion, 'path>
    MovingCompactRenderer<'history, &'motion PreparedReflectedSourceTrajectory<'path>>
{
    /// Render actual reflected pressure from original physical source histories.
    ///
    /// All sources share one original fixed microphone and ideal-plane group.
    /// Each strength interpolant is compiled ONCE per physical source, then
    /// queried at every path's own retarded time. Complete path work is charged
    /// to the original source/receiver-pair budget and retained root caps.
    /// Existing coherent sums, finite history, bandwidth admission, transactional
    /// blocks, physical clock, decimation and PCM sinks remain unchanged.
    pub fn with_ideal_reflections<Cancel: FnMut() -> bool>(
        sources: &[CompactSourceHistory<
            'history,
            &'motion PreparedReflectedSourceTrajectory<'path>,
        >],
        config: MovingPressureConfig,
        budget: MovingPressureBudget,
        cancelled: &mut Cancel,
    ) -> Result<Self, RenderError> {
        check_cancel(cancelled)?;
        let first = sources
            .first()
            .ok_or_else(|| sizing("reflected observation needs a physical source"))?;
        Self::new(
            first.motion.reflection_observer(),
            sources,
            config,
            budget,
            cancelled,
        )
    }
}
