//! Fixed fluid-frame receiver geometry for force-history prefixes.

use super::super::{
    CompactMotion, CompactSourceMotion, acoustic, admit_clearance, check_cancel, sizing,
};
use crate::render::RenderError;
use fs_aeroac::curle3d::moving::trajectory::{
    PreparedReflectedSourceTrajectory, PreparedSourceTrajectory,
};
use fs_aeroac::curle3d::moving::{MovingCompactObserver, UniformSourceMotion};

/// One original direct or image receiver for a fixed-observer history integral.
#[doc(hidden)]
#[derive(Clone, Copy)]
pub struct ForceHistoryPath {
    pub(in crate::pcm_wav::observation) observer: MovingCompactObserver,
    pub(super) pressure_sign: f64,
    // Includes pressure parity and the inverse image Jacobian exactly once.
    pub(super) velocity_transform: Option<[[f64; 3]; 3]>,
    pub(in crate::pcm_wav::observation) source_end_s: f64,
    pub(in crate::pcm_wav::observation) maximum_source_speed_m_s: f64,
    // Zero for fixed paths whose complete past geometry was admitted earlier.
    pub(in crate::pcm_wav::observation) minimum_history_distance_m: f64,
}

/// Source motion whose every propagation path has one fixed fluid-frame receiver.
///
/// Implemented for uniform sources, accepted curved trajectories, their mixed
/// ensemble, and ideal reflections of accepted curved trajectories. Moving
/// receivers and laboratory-fixed receivers in a mean flow do not satisfy this
/// contract: their history kernel changes its spatial evaluation point.
/// The original sealed motion owners retain continuous speed and clearance.
pub trait FixedForceMotion: CompactMotion {
    /// Admit the ENTIRE retained force prefix, not only current retarded events.
    #[doc(hidden)]
    fn admit_force_history(
        self,
        observer: &Self::Observer,
        path: usize,
        source_interval_s: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ForceHistoryPath, RenderError>;

    /// Original geometry knot strictly after time; uniform motion has none.
    #[doc(hidden)]
    fn next_source_knot_after(self, time: f64) -> Result<Option<f64>, RenderError>;
}

/// Original motion accepted by the full-field force-history renderer.
///
/// Fixed microphones retain their immutable prefix integrals. Native moving
/// microphones retain only the original knot partition and integrate at each
/// current virtual receiver position. This trait is sealed by `CompactMotion`.
pub trait ForceMotion: CompactMotion {
    /// Whether the spatial history kernel changes between receiver events.
    #[doc(hidden)]
    const MOVING_HISTORY_KERNEL: bool;

    /// Bind the current image geometry and complete retained source-time domain.
    #[doc(hidden)]
    fn force_history_path_at(
        self,
        observer: &Self::Observer,
        path: usize,
        receiver_time_s: f64,
        source_interval_s: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ForceHistoryPath, RenderError>;

    /// Next original geometry knot, never an adaptive quadrature node.
    #[doc(hidden)]
    fn next_force_history_knot_after(self, time: f64) -> Result<Option<f64>, RenderError>;
}

// Preserve every existing fixed-receiver implementation and its fast path.
impl<M: FixedForceMotion> ForceMotion for M {
    const MOVING_HISTORY_KERNEL: bool = false;

    fn force_history_path_at(
        self,
        observer: &Self::Observer,
        path: usize,
        _receiver_time_s: f64,
        source_interval_s: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ForceHistoryPath, RenderError> {
        self.admit_force_history(observer, path, source_interval_s, cancelled)
    }

    fn next_force_history_knot_after(self, time: f64) -> Result<Option<f64>, RenderError> {
        self.next_source_knot_after(time)
    }
}

fn trajectory_speed_bound(trajectory: &PreparedSourceTrajectory) -> f64 {
    (trajectory.maximum_source_mach() * trajectory.observer().medium().sound_speed_m_s()).next_up()
}

fn interval(interval: [f64; 2]) -> Result<(), RenderError> {
    if interval.iter().any(|time| !time.is_finite()) || interval[0] > interval[1] {
        return Err(sizing(
            "force history requires a finite ordered source interval",
        ));
    }
    Ok(())
}

fn retained_end(
    trajectory: &PreparedSourceTrajectory,
    requested: [f64; 2],
) -> Result<f64, RenderError> {
    interval(requested)?;
    let [first, last] = trajectory.source_time_interval_s();
    let end = requested[1].min(last);
    if requested[0] < first || requested[0] > end {
        return Err(sizing(
            "force prehistory lacks its original source geometry",
        ));
    }
    // The prepared owner admitted continuous clearance throughout this domain.
    Ok(end)
}

impl FixedForceMotion for UniformSourceMotion {
    fn admit_force_history(
        self,
        observer: &MovingCompactObserver,
        path: usize,
        source_interval_s: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ForceHistoryPath, RenderError> {
        check_cancel(cancelled)?;
        interval(source_interval_s)?;
        if path != 0 {
            return Err(sizing("uniform force history has only its direct path"));
        }
        for time in source_interval_s {
            observer
                .prepare_emission(time, self.kinematics_at(time).map_err(acoustic)?)
                .map_err(acoustic)?;
        }
        admit_clearance(observer, &self, source_interval_s[0], source_interval_s[1])?;
        Ok(ForceHistoryPath {
            observer: *observer,
            pressure_sign: 1.0,
            velocity_transform: None,
            source_end_s: source_interval_s[1],
            minimum_history_distance_m: 0.0,
            maximum_source_speed_m_s: self
                .velocity_m_s()
                .iter()
                .fold(0.0_f64, |sum, v| (sum + v.abs()).next_up()),
        })
    }

    fn next_source_knot_after(self, _time: f64) -> Result<Option<f64>, RenderError> {
        Ok(None)
    }
}

impl FixedForceMotion for &PreparedSourceTrajectory {
    fn admit_force_history(
        self,
        observer: &MovingCompactObserver,
        path: usize,
        source_interval_s: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ForceHistoryPath, RenderError> {
        check_cancel(cancelled)?;
        self.admit_observer(observer)?;
        if path != 0 {
            return Err(sizing("free-space force history has only its direct path"));
        }
        Ok(ForceHistoryPath {
            observer: *observer,
            pressure_sign: 1.0,
            velocity_transform: None,
            source_end_s: retained_end(self, source_interval_s)?,
            maximum_source_speed_m_s: trajectory_speed_bound(self),
            minimum_history_distance_m: 0.0,
        })
    }

    fn next_source_knot_after(self, time: f64) -> Result<Option<f64>, RenderError> {
        PreparedSourceTrajectory::next_source_knot_after(self, time).map_err(acoustic)
    }
}

impl FixedForceMotion for CompactSourceMotion<'_> {
    fn admit_force_history(
        self,
        observer: &MovingCompactObserver,
        path: usize,
        source_interval_s: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ForceHistoryPath, RenderError> {
        match self {
            Self::Uniform(m) => m.admit_force_history(observer, path, source_interval_s, cancelled),
            Self::Trajectory(m) => {
                m.admit_force_history(observer, path, source_interval_s, cancelled)
            }
        }
    }

    fn next_source_knot_after(self, time: f64) -> Result<Option<f64>, RenderError> {
        match self {
            Self::Uniform(m) => m.next_source_knot_after(time),
            Self::Trajectory(m) => FixedForceMotion::next_source_knot_after(m, time),
        }
    }
}

impl FixedForceMotion for &PreparedReflectedSourceTrajectory<'_> {
    fn admit_force_history(
        self,
        observer: &Self::Observer,
        path: usize,
        source_interval_s: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ForceHistoryPath, RenderError> {
        check_cancel(cancelled)?;
        self.admit_observer(observer)?;
        let source_end_s = retained_end(self.trajectory(), source_interval_s)?;
        let source = self
            .trajectory()
            .kinematics_at(source_interval_s[0])
            .map_err(acoustic)?;
        let (virtual_receiver, pressure_sign, velocity_transform) = observer
            .stationary_field_path(path, source.position_m)
            .map_err(acoustic)?;
        // The reflected trajectory owner admitted every image throughout the
        // original curve, including halfspaces and continuous image clearance.
        Ok(ForceHistoryPath {
            observer: virtual_receiver,
            pressure_sign,
            velocity_transform: (path != 0).then_some(velocity_transform),
            source_end_s,
            maximum_source_speed_m_s: trajectory_speed_bound(self.trajectory()),
            minimum_history_distance_m: 0.0,
        })
    }

    fn next_source_knot_after(self, time: f64) -> Result<Option<f64>, RenderError> {
        self.trajectory()
            .next_source_knot_after(time)
            .map_err(acoustic)
    }
}
