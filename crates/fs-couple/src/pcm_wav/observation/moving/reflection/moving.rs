//! Original moving-microphone image paths in the shared pressure renderer.

use super::super::motion::{CompactMotion, sealed};
use super::super::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
    acoustic, check_cancel, sizing,
};
use crate::render::RenderError;
use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::trajectory::{
    MovingReflectionObserver, PreparedReflectedMovingReceiverSource,
};
use fs_aeroac::curle3d::moving::{
    MovingCompactObserver, PreparedMovingCompactSample, SourceKinematics,
};

impl sealed::Sealed for &PreparedReflectedMovingReceiverSource<'_> {}

impl<'path> CompactMotion for &PreparedReflectedMovingReceiverSource<'path> {
    type Observer = &'path MovingReflectionObserver<'path>;

    fn admit_observer(self, observer: &Self::Observer) -> Result<(), RenderError> {
        if !core::ptr::eq(self.reflection_observer(), *observer) {
            return Err(sizing(
                "reflected moving sources must share the original microphone path and walls",
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
        PreparedReflectedMovingReceiverSource::observer_at(self, time).map_err(acoustic)
    }

    fn admit_clock(self, first: f64, last: f64, tolerance: f64) -> Result<(), RenderError> {
        let config = self.config();
        let [start, end] = config.receiver_time_interval_s;
        if !first.is_finite() || !last.is_finite() || first < start || last > end || first > last {
            return Err(sizing(
                "pressure clock exceeds the reflected moving receiver window",
            ));
        }
        if !tolerance.is_finite() || tolerance < 0.0 || config.maximum_time_residual_s > tolerance {
            return Err(sizing(
                "reflected moving receiver root tolerance exceeds the pressure clock budget",
            ));
        }
        Ok(())
    }

    fn maximum_root_iterations(self) -> usize {
        PreparedReflectedMovingReceiverSource::maximum_root_iterations(self)
    }
    fn path_count(self) -> usize {
        PreparedReflectedMovingReceiverSource::path_count(self)
    }
    fn kinematics_at(self, time: f64) -> Result<SourceKinematics, AeroacError> {
        self.source().kinematics_at(time)
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
        self.prepare_path(path, time, cancelled)
            .map(|sample| sample.pressure)
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
        self.prepare_path_emission(path, source_time, receiver_time)
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
        let [start, end] = self.config().receiver_time_interval_s;
        if first.receiver_time_s() < start
            || last.receiver_time_s() > end
            || first.receiver_time_s() > last.receiver_time_s()
            || first.source_time_s() > last.source_time_s()
        {
            return Err(sizing(
                "pressure clock exceeds complete reflected moving history",
            ));
        }
        // All image directions and reflected receiver velocities participate.
        // Reusing only the free-space rate can understate reflected bandwidth.
        Ok(self.maximum_source_time_rate())
    }
}

impl<'history, 'motion, 'path>
    MovingCompactRenderer<'history, &'motion PreparedReflectedMovingReceiverSource<'path>>
{
    /// Render coherent ideal-wall pressure at the SAME original moving
    /// microphone. Each physical strength interpolant is compiled once and
    /// read at every path's own source time. Existing source/path work budgets,
    /// global sample clock, compensated sums and transactional blocks apply.
    pub fn with_moving_receiver_reflections<Cancel: FnMut() -> bool>(
        sources: &[CompactSourceHistory<
            'history,
            &'motion PreparedReflectedMovingReceiverSource<'path>,
        >],
        config: MovingPressureConfig,
        budget: MovingPressureBudget,
        cancelled: &mut Cancel,
    ) -> Result<Self, RenderError> {
        check_cancel(cancelled)?;
        let first = sources
            .first()
            .ok_or_else(|| sizing("reflected moving observation needs a physical source"))?;
        Self::new(
            first.motion.reflection_observer(),
            sources,
            config,
            budget,
            cancelled,
        )
    }
}
