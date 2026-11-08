//! Coherent full fields from independently moving sources at one microphone.
//!
//! Every source keeps its original causal propagation and force-history owner.
//! Pressure and velocity from ALL sources and images stay in double-double
//! arithmetic until the final field and signed instantaneous p*u. Adding the
//! rounded public single-source fields can lose a weak force contribution next
//! to a volume source before a different source cancels that volume field.

use super::reflection::ReflectedFieldSum;
use super::{
    PreparedKinematicTrajectory, PreparedMovingReceiverSource,
    PreparedReflectedMovingReceiverSource, invalid, poll,
};
use crate::AeroacError;
use crate::curle3d::moving::{CompactSourceStrength, MovingCompactObserver};

mod multipole;
pub use multipole::{MovingMultipoleHistory, MovingMultipoleSource};

/// Aggregate work admission for one complete receiver event.
/// Native source/receiver preparation budgets remain in force as well.
#[derive(Clone, Copy, Debug)]
pub struct MovingSourceArrayBudget {
    /// Maximum original physical source count.
    pub maximum_sources: usize,
    /// Maximum total direct AND reflected paths across all original sources.
    pub maximum_source_paths: usize,
    /// Sum of each native path's maximum interior causal-root refinements.
    /// Original bounded knot searches and constant endpoint work are additional.
    /// User strength/history callbacks retain their own work and accuracy limits.
    pub maximum_retarded_iterations: usize,
}

/// Simultaneous coherent physical fields at the original microphone event.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingSourceArrayField {
    pub receiver_time_s: f64,
    /// Original source count; image paths are never independent sources.
    pub source_count: usize,
    pub source_path_count: usize,
    /// Total instantaneous pressure [Pa].
    pub pressure_pa: f64,
    /// Total Cartesian acoustic particle velocity in the fluid-rest frame [m/s].
    pub particle_velocity_m_s: [f64; 3],
    /// Signed total instantaneous pressure-work flux [W/m²], formed after sums.
    /// No harmonic half-product, plane-wave conversion or mean-flow advection.
    pub instantaneous_intensity_w_m2: [f64; 3],
    /// Largest absolute native retarded-time residual over every source/path.
    pub maximum_time_residual_s: f64,
}

#[derive(Clone, Copy, Debug)]
enum Sources<'a> {
    Direct(&'a [PreparedMovingReceiverSource<'a>]),
    Reflected(&'a [PreparedReflectedMovingReceiverSource<'a>]),
}

#[derive(Clone, Copy)]
enum Source<'a> {
    Direct(&'a PreparedMovingReceiverSource<'a>),
    Reflected(&'a PreparedReflectedMovingReceiverSource<'a>),
}

impl<'a> Sources<'a> {
    fn len(self) -> usize {
        match self {
            Self::Direct(sources) => sources.len(),
            Self::Reflected(sources) => sources.len(),
        }
    }

    fn at(self, index: usize) -> Source<'a> {
        match self {
            Self::Direct(sources) => Source::Direct(&sources[index]),
            Self::Reflected(sources) => Source::Reflected(&sources[index]),
        }
    }
}

impl<'a> Source<'a> {
    fn receiver(self) -> &'a PreparedKinematicTrajectory {
        match self {
            Self::Direct(source) => source.receiver(),
            Self::Reflected(source) => source.receiver(),
        }
    }

    fn window(self) -> [f64; 2] {
        match self {
            Self::Direct(source) => source.config().receiver_time_interval_s,
            Self::Reflected(source) => source.config().receiver_time_interval_s,
        }
    }

    fn path_count(self) -> usize {
        match self {
            Self::Direct(_) => 1,
            Self::Reflected(source) => source.path_count(),
        }
    }

    fn root_iterations(self) -> usize {
        match self {
            Self::Direct(source) => source.maximum_root_iterations(),
            Self::Reflected(source) => source.maximum_root_iterations(),
        }
    }

    fn same_boundary(self, other: Self) -> bool {
        match (self, other) {
            (Self::Direct(_), Self::Direct(_)) => true,
            (Self::Reflected(left), Self::Reflected(right)) => {
                core::ptr::eq(left.reflection_observer(), right.reflection_observer())
            }
            _ => false,
        }
    }
}

/// Borrowed coherent source array on ONE original moving microphone trajectory.
///
/// Sources may have independent motion, accepted knots and strengths. Every
/// source must bind the SAME receiver object, so the original medium, rest frame
/// and microphone clock cannot silently differ. Reflected sources additionally
/// share the same admitted image owner and its complete path ordering. A fixed
/// microphone is represented by its original stationary kinematic trajectory.
///
/// Construction and queries allocate nothing and never recompile paths. The
/// common observation interval is the intersection of native source windows;
/// no missing prehistory or tail is treated as silence. Source strengths and
/// force integrals remain caller-owned physical histories, not new interpolants.
/// This is deterministic superposition, not an inferred source-noise law.
#[derive(Clone, Copy, Debug)]
pub struct PreparedMovingSourceArray<'a> {
    sources: Sources<'a>,
    receiver: &'a PreparedKinematicTrajectory,
    receiver_time_interval_s: [f64; 2],
    source_path_count: usize,
    maximum_retarded_iterations: usize,
}

impl<'a> PreparedMovingSourceArray<'a> {
    /// Bind already admitted free-space source/receiver pairs in source order.
    /// At least one source is required to establish the physical receiver.
    pub fn new(
        sources: &'a [PreparedMovingReceiverSource<'a>],
        budget: MovingSourceArrayBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::prepare(Sources::Direct(sources), budget, cancelled)
    }

    /// Bind complete orthogonal or pi/m wedge groups in original source order.
    /// All sources borrow the same MovingReflectionObserver; equal path
    /// counts alone never establish common boundary geometry or image identity.
    pub fn new_reflected(
        sources: &'a [PreparedReflectedMovingReceiverSource<'a>],
        budget: MovingSourceArrayBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::prepare(Sources::Reflected(sources), budget, cancelled)
    }

    fn prepare(
        sources: Sources<'a>,
        budget: MovingSourceArrayBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        if sources.len() == 0 || sources.len() > budget.maximum_sources {
            return Err(invalid(
                "moving source array needs a nonempty bounded source list",
            ));
        }
        let first = sources.at(0);
        let receiver = first.receiver();
        let mut window = first.window();
        let mut path_count = 0usize;
        let mut iterations = 0usize;
        for index in 0..sources.len() {
            poll(cancelled)?;
            let source = sources.at(index);
            if !core::ptr::eq(receiver, source.receiver()) || !first.same_boundary(source) {
                return Err(invalid(
                    "moving source array requires one original receiver and boundary group",
                ));
            }
            let interval = source.window();
            window[0] = window[0].max(interval[0]);
            window[1] = window[1].min(interval[1]);
            path_count = path_count
                .checked_add(source.path_count())
                .ok_or_else(|| invalid("moving source array path count overflows"))?;
            iterations = source
                .path_count()
                .checked_mul(source.root_iterations())
                .and_then(|work| iterations.checked_add(work))
                .ok_or_else(|| invalid("moving source array root work overflows"))?;
            if path_count > budget.maximum_source_paths
                || iterations > budget.maximum_retarded_iterations
            {
                return Err(invalid(
                    "moving source array exceeds aggregate path or root work",
                ));
            }
        }
        if window[0] > window[1] {
            return Err(invalid(
                "moving sources have no complete common receiver window",
            ));
        }
        poll(cancelled)?;
        Ok(Self {
            sources,
            receiver,
            receiver_time_interval_s: window,
            source_path_count: path_count,
            maximum_retarded_iterations: iterations,
        })
    }

    pub fn source_count(&self) -> usize {
        self.sources.len()
    }
    pub const fn source_path_count(&self) -> usize {
        self.source_path_count
    }
    pub const fn maximum_retarded_iterations(&self) -> usize {
        self.maximum_retarded_iterations
    }
    pub const fn receiver(&self) -> &'a PreparedKinematicTrajectory {
        self.receiver
    }
    /// Finite physical interval covered by every original source and image.
    pub const fn receiver_time_interval_s(&self) -> [f64; 2] {
        self.receiver_time_interval_s
    }

    /// Observe all original Q/F sources at one common microphone event.
    ///
    /// Source indices exactly follow the constructor's slice. Each callback is
    /// queried at that source/path's OWN retarded emission time on the original
    /// physical clock. The force-history callback additionally receives its
    /// image index and actual virtual observer at this receiver event; direct
    /// propagation uses image index zero. Supply J [N s/m³] from the original
    /// moving-force owner's kernel, without image parity or vector rotation.
    /// The native propagation owner applies both once after evaluation.
    ///
    /// A moving microphone changes the spatial history kernel at every event;
    /// a previous microphone position's integral cannot be reused. Continuous
    /// past-worldline clearance, integration accuracy and explicit initial
    /// forcing remain caller obligations. No source or image can fail silently:
    /// any callback, numeric or cancellation failure refuses the whole result.
    pub fn force_field_at(
        &self,
        receiver_time: f64,
        source_at: &mut impl FnMut(usize, f64) -> Result<CompactSourceStrength, AeroacError>,
        force_history_at: &mut impl FnMut(
            usize,
            usize,
            &MovingCompactObserver,
            f64,
        ) -> Result<[f64; 3], AeroacError>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingSourceArrayField, AeroacError> {
        poll(cancelled)?;
        let [first, last] = self.receiver_time_interval_s;
        if !receiver_time.is_finite() || receiver_time < first || receiver_time > last {
            return Err(invalid(
                "receiver time lacks complete moving-source array history",
            ));
        }
        let mut sum = ReflectedFieldSum::default();
        for index in 0..self.source_count() {
            poll(cancelled)?;
            match self.sources.at(index) {
                Source::Direct(source) => {
                    let original = source.prepare_receiver(receiver_time, cancelled)?.pressure;
                    let time = original.source_time_s();
                    let residual = original.retarded_time_residual_s();
                    let observer = source.observer_at(receiver_time)?;
                    let prepared = observer
                        .prepare_force_sample(original, source.source().kinematics_at(time)?)?;
                    poll(cancelled)?;
                    let strength = source_at(index, time)?;
                    poll(cancelled)?;
                    let history = force_history_at(index, 0, &observer, time)?;
                    let (_, pressure, velocity) =
                        prepared.evaluate_with_fields_dd(strength, history)?;
                    sum.add(pressure, velocity, residual);
                }
                Source::Reflected(source) => source.add_force_field(
                    receiver_time,
                    &mut |time| source_at(index, time),
                    &mut |path, observer, time| force_history_at(index, path, observer, time),
                    &mut sum,
                    cancelled,
                )?,
            }
        }
        poll(cancelled)?;
        let field = sum.finish(receiver_time)?;
        Ok(MovingSourceArrayField {
            receiver_time_s: receiver_time,
            source_count: self.source_count(),
            source_path_count: self.source_path_count,
            pressure_pa: field.pressure_pa,
            particle_velocity_m_s: field.particle_velocity_m_s,
            instantaneous_intensity_w_m2: field.instantaneous_intensity_w_m2,
            maximum_time_residual_s: field.maximum_time_residual_s,
        })
    }
}
