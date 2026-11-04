//! Fixed laboratory histories in one original uniform-flow receiver.

use super::motion::{CompactMotion, sealed};
use super::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
    acoustic, check_cancel, sizing,
};
use crate::render::RenderError;
use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::flow::{UniformFlowCompactObserver, UniformFlowReceiver};
use fs_aeroac::curle3d::moving::{
    MovingCompactLimits, MovingCompactObserver, PreparedMovingCompactSample, SourceKinematics,
    UniformSourceMotion,
};

impl sealed::Sealed for UniformFlowCompactObserver {}
impl CompactMotion for UniformFlowCompactObserver {
    type Observer = UniformFlowReceiver;

    fn admit_observer(self, observer: &Self::Observer) -> Result<(), RenderError> {
        if self.receiver() != *observer {
            return Err(sizing(
                "uniform-flow sources must share their receiver, fluid and velocity",
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
        let receiver = UniformSourceMotion::new(
            0.0,
            observer.observer_position_m(),
            observer.velocity_m_s().map(|v| -v),
        )
        .and_then(|motion| motion.kinematics_at(time))
        .map_err(acoustic)?;
        // This public query returns fluid-frame geometry. Rendering and endpoint
        // reconciliation use the cached relative path, never these potentially
        // large absolute Galilean positions.
        let distance = self
            .prepare_emission(0.0)
            .map_err(acoustic)?
            .pressure_in_fluid_frame()
            .distance_m();
        MovingCompactObserver::new(
            observer.medium(),
            receiver.position_m,
            MovingCompactLimits {
                minimum_distance_m: distance,
                ..observer.limits()
            },
        )
        .map_err(acoustic)
    }
    fn maximum_root_iterations(self) -> usize {
        0
    }
    fn kinematics_at(self, time: f64) -> Result<SourceKinematics, AeroacError> {
        UniformSourceMotion::new(
            0.0,
            self.source_position_m(),
            self.velocity_m_s().map(|v| -v),
        )?
        .kinematics_at(time)
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
        self.prepare_receiver(time, tolerance)
            .map(|sample| sample.pressure_in_fluid_frame())
            .map_err(acoustic)
    }
    fn prepare_endpoint(
        self,
        observer: &Self::Observer,
        source_time: f64,
        _receiver_time: f64,
    ) -> Result<PreparedMovingCompactSample, RenderError> {
        self.admit_observer(observer)?;
        self.prepare_emission(source_time)
            .map(|sample| sample.pressure_in_fluid_frame())
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
        if first.source_time_s() > last.source_time_s()
            || first.receiver_time_s() > last.receiver_time_s()
        {
            return Err(sizing(
                "uniform-flow observation clocks lost their physical order",
            ));
        }
        // Source and receiver are both fixed in the laboratory: wind changes
        // their path but causes neither Doppler shift nor amplitude modulation.
        Ok(1.0)
    }
}

impl<'history> MovingCompactRenderer<'history, UniformFlowCompactObserver> {
    /// Observe accepted fixed-source histories at their common flow receiver.
    /// Reuses the original Hermite interpolation, finite history/bandwidth
    /// admission, coherent pressure sum, work limits and transactional callbacks.
    pub fn with_uniform_flow<Cancel: FnMut() -> bool>(
        sources: &[CompactSourceHistory<'history, UniformFlowCompactObserver>],
        config: MovingPressureConfig,
        budget: MovingPressureBudget,
        cancelled: &mut Cancel,
    ) -> Result<Self, RenderError> {
        check_cancel(cancelled)?;
        let first = sources
            .first()
            .ok_or_else(|| sizing("uniform-flow observation needs a source"))?;
        Self::new(first.motion.receiver(), sources, config, budget, cancelled)
    }
}
