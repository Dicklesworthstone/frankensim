//! The original moving microphone owner in the shared pressure renderer.

use super::motion::{CompactMotion, sealed};
use super::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
    acoustic, check_cancel, sizing,
};
use crate::render::RenderError;
use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::trajectory::{
    PreparedKinematicTrajectory, PreparedMovingReceiverSource,
};
use fs_aeroac::curle3d::moving::{
    MovingCompactLimits, MovingCompactObserver, PreparedMovingCompactSample, SourceKinematics,
};

impl sealed::Sealed for &PreparedMovingReceiverSource<'_> {}

impl<'path> CompactMotion for &PreparedMovingReceiverSource<'path> {
    type Observer = &'path PreparedKinematicTrajectory;

    fn admit_observer(self, observer: &Self::Observer) -> Result<(), RenderError> {
        // Equal endpoints are not sufficient: they could have different interior
        // trajectories, accepted clocks, interpolation or preparation policies.
        if !core::ptr::eq(self.receiver(), *observer) {
            return Err(sizing(
                "coherent moving sources must share the original microphone trajectory",
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
        let config = self.config();
        let [first, last] = config.receiver_time_interval_s;
        if !time.is_finite() || time < first || time > last {
            return Err(sizing("moving microphone time exceeds its prepared window"));
        }
        let receiver = self.receiver().kinematics_at(time).map_err(acoustic)?;
        MovingCompactObserver::new(
            self.source().medium(),
            receiver.position_m,
            MovingCompactLimits {
                maximum_mach: self.source().maximum_mach(),
                minimum_distance_m: config.minimum_distance_m,
            },
        )
        .map_err(acoustic)
    }

    fn maximum_root_iterations(self) -> usize {
        PreparedMovingReceiverSource::maximum_root_iterations(self)
    }

    fn admit_clock(self, first: f64, last: f64, tolerance: f64) -> Result<(), RenderError> {
        let config = self.config();
        let [start, end] = config.receiver_time_interval_s;
        if !first.is_finite() || !last.is_finite() || first < start || last > end || first > last {
            return Err(sizing(
                "pressure clock exceeds the prepared moving receiver window",
            ));
        }
        if !tolerance.is_finite() || tolerance < 0.0 || config.maximum_time_residual_s > tolerance {
            return Err(sizing(
                "moving receiver root tolerance exceeds the pressure clock's residual budget",
            ));
        }
        Ok(())
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
        check_cancel(cancelled)?;
        self.admit_observer(observer)?;
        self.admit_clock(time, time, tolerance)?;
        // The physical pair checks the retarded root and both trajectories.
        // Its clock rate is used for bandwidth admission, never pressure gain.
        self.prepare_receiver(time, cancelled)
            .map(|sample| sample.pressure)
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
                "pressure clock exceeds the prepared moving source/receiver window",
            ));
        }
        // Both paths can turn: the first instantaneous rate is insufficient.
        // The prepared pair retains (1 + M_receiver)/(1 - M_source).
        Ok(self.maximum_source_time_rate())
    }
}

impl<'history, 'motion, 'path>
    MovingCompactRenderer<'history, &'motion PreparedMovingReceiverSource<'path>>
{
    /// Observe accepted sources at their SAME original moving microphone.
    ///
    /// The source/receiver pairs already own continuous motion, clearance and
    /// causal-coverage admission. This constructor borrows their receiver,
    /// checks each pair against it, and reuses the original source interpolator,
    /// coherent pressure sum, global sample clock and transactional callbacks.
    /// No mechanics state, trajectory or strength history is cloned or advanced.
    pub fn with_moving_receiver<Cancel: FnMut() -> bool>(
        sources: &[CompactSourceHistory<'history, &'motion PreparedMovingReceiverSource<'path>>],
        config: MovingPressureConfig,
        budget: MovingPressureBudget,
        cancelled: &mut Cancel,
    ) -> Result<Self, RenderError> {
        check_cancel(cancelled)?;
        let first = sources
            .first()
            .ok_or_else(|| sizing("moving microphone observation needs a physical source"))?;
        Self::new(first.motion.receiver(), sources, config, budget, cancelled)
    }
}
