//! Complete past-force integration at the original moving microphone event.

use super::super::{
    CompactMotion, CompactSourceHistory, MovingPressureBudget, MovingPressureConfig, acoustic,
    check_cancel, sizing,
};
use super::motion::{ForceHistoryPath, ForceMotion};
use super::{ForceHistoryQuadrature, MovingForceRenderer};
use crate::render::RenderError;
use fs_aeroac::curle3d::moving::MovingCompactObserver;
use fs_aeroac::curle3d::moving::trajectory::{
    PreparedKinematicTrajectory, PreparedMovingReceiverSource,
    PreparedReflectedMovingReceiverSource,
};
use fs_bem::near_field::FirstOrder;

fn path_geometry(
    source: &PreparedKinematicTrajectory,
    observer: MovingCompactObserver,
    pressure_sign: f64,
    velocity_transform: Option<[[f64; 3]; 3]>,
    interval: [f64; 2],
) -> Result<ForceHistoryPath, RenderError> {
    let [first, last] = source.time_interval_s();
    let end = interval[1].min(last);
    if interval.iter().any(|t| !t.is_finite()) || interval[0] < first || interval[0] > end {
        return Err(sizing(
            "moving microphone force prehistory lacks original source geometry",
        ));
    }
    let speed = source.maximum_mach() * source.medium().sound_speed_m_s();
    Ok(ForceHistoryPath {
        observer,
        pressure_sign,
        velocity_transform,
        source_end_s: end,
        maximum_source_speed_m_s: if speed == 0.0 { 0.0 } else { speed.next_up() },
        // Causal pair clearance does not bound earlier source positions at
        // today's microphone point. The quadrature admits every history cell.
        minimum_history_distance_m: observer.limits().minimum_distance_m,
    })
}

impl ForceMotion for &PreparedMovingReceiverSource<'_> {
    const MOVING_HISTORY_KERNEL: bool = true;

    fn force_history_path_at(
        self,
        observer: &Self::Observer,
        path: usize,
        receiver_time_s: f64,
        source_interval_s: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ForceHistoryPath, RenderError> {
        check_cancel(cancelled)?;
        if path != 0 {
            return Err(sizing("free-space force history has only its direct path"));
        }
        path_geometry(
            self.source(),
            self.observer_at(observer, receiver_time_s)?,
            1.0,
            None,
            source_interval_s,
        )
    }

    fn next_force_history_knot_after(self, time: f64) -> Result<Option<f64>, RenderError> {
        self.source().next_source_knot_after(time).map_err(acoustic)
    }
}

impl ForceMotion for &PreparedReflectedMovingReceiverSource<'_> {
    const MOVING_HISTORY_KERNEL: bool = true;

    fn force_history_path_at(
        self,
        observer: &Self::Observer,
        path: usize,
        receiver_time_s: f64,
        source_interval_s: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ForceHistoryPath, RenderError> {
        check_cancel(cancelled)?;
        self.admit_observer(observer)?;
        let (virtual_receiver, parity, transform) =
            PreparedReflectedMovingReceiverSource::force_history_path_at(
                self,
                path,
                receiver_time_s,
            )
            .map_err(acoustic)?;
        path_geometry(
            self.source(),
            virtual_receiver,
            parity,
            (path != 0).then_some(transform),
            source_interval_s,
        )
    }

    fn next_force_history_knot_after(self, time: f64) -> Result<Option<f64>, RenderError> {
        self.source().next_source_knot_after(time).map_err(acoustic)
    }
}

impl<'history, 'motion, 'path>
    MovingForceRenderer<'history, &'motion PreparedMovingReceiverSource<'path>>
{
    /// Render complete Q/F fields and directional audio from accepted histories
    /// at their SAME original moving microphone, with explicit zero force
    /// prehistory. No caller integration, new motion, or second strength
    /// interpolant is required. The original fluid frame is retained.
    pub fn with_moving_receiver(
        sources: &[CompactSourceHistory<'history, &'motion PreparedMovingReceiverSource<'path>>],
        microphone: FirstOrder,
        config: MovingPressureConfig,
        budget: MovingPressureBudget,
        quadrature: ForceHistoryQuadrature,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RenderError> {
        check_cancel(cancelled)?;
        let first = sources
            .first()
            .ok_or_else(|| sizing("moving force observation needs a physical source"))?;
        Self::new_with_zero_force_prehistory(
            first.motion.receiver(),
            sources,
            microphone,
            config,
            budget,
            quadrature,
            cancelled,
        )
    }
}

impl<'history, 'motion, 'path>
    MovingForceRenderer<'history, &'motion PreparedReflectedMovingReceiverSource<'path>>
{
    /// Render complete ideal-wall Q/F interference at the SAME original moving
    /// microphone. Each source/image integrates its own causal past at its
    /// current virtual receiver. Pressure parity and polar-vector transforms
    /// enter once, before coherent intensity and directional pickup.
    pub fn with_moving_receiver_reflections(
        sources: &[CompactSourceHistory<
            'history,
            &'motion PreparedReflectedMovingReceiverSource<'path>,
        >],
        microphone: FirstOrder,
        config: MovingPressureConfig,
        budget: MovingPressureBudget,
        quadrature: ForceHistoryQuadrature,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RenderError> {
        check_cancel(cancelled)?;
        let first = sources
            .first()
            .ok_or_else(|| sizing("reflected moving force observation needs a physical source"))?;
        Self::new_with_zero_force_prehistory(
            first.motion.reflection_observer(),
            sources,
            microphone,
            config,
            budget,
            quadrature,
            cancelled,
        )
    }
}
