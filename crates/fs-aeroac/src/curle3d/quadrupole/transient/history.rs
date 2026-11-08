//! Original finite stress records as consistent transient acoustic histories.
//!
//! The six physical tensor entries are reconstructed by the same clamped C2
//! cubic spline on the original uniform source clock. Endpoint rates and the
//! integral at the first sample are explicit physical inputs. Q, Qdot, Qddot
//! and the integral are analytic evaluations of that ONE reconstruction.
//! No noisy sample-by-sample differentiation is performed during observation.
//!
//! This is an offline interpolation model, not an antialias filter or a
//! causal streaming predictor: its spline uses the complete supplied record.
//! Samples alone do not establish source bandwidth, derivative accuracy or
//! turbulence/source validity. Refine the original record for those claims.
//! There is no extrapolation, implicit periodic continuation, fabricated
//! silence before the record, or resetting of the near-field integral.

use super::{
    PreparedTransientQuadrupoleObserver, RetardedQuadrupoleField, TransientQuadrupoleSource,
    invalid, nonfinite, poll,
};
use crate::AeroacError;
use crate::curle3d::quadrupole::recorded::QuadrupoleSampleClock;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::dd::Dd;

/// The physical conditions completing one source's finite sampled history.
#[derive(Clone, Copy, Debug)]
pub struct QuadrupoleHistoryEnds {
    /// Qdot at the first recorded instant [N m/s], in original tensor order.
    pub initial_rate_n_m_s: [f64; 6],
    /// Qdot at the last recorded instant [N m/s]. This does not define a tail.
    pub final_rate_n_m_s: [f64; 6],
    /// Integral of Q through the first instant [N m s]. Zero declares that
    /// initial integral; it does not infer unrecorded earlier source values.
    pub initial_integral_n_m_s: [f64; 6],
}

/// Explicit source, time and requested backing-buffer bounds.
#[derive(Clone, Copy, Debug)]
pub struct QuadrupoleHistoryBudget {
    pub maximum_sources: usize,
    pub maximum_samples_per_source: usize,
    /// Includes retained nodes/positions, common tridiagonal factors and the
    /// reused six-channel preparation workspace. Excludes caller inputs,
    /// allocator metadata and the separately budgeted propagation observer.
    pub maximum_workspace_bytes: usize,
}

/// End conditions carry each channel's own value/time and value*time units.
#[derive(Clone, Copy)]
pub(super) struct CubicHistoryEnds<const N: usize> {
    pub(super) initial_rate: [f64; N],
    pub(super) final_rate: [f64; N],
    pub(super) initial_integral: [f64; N],
}

pub(super) struct CubicSourceSample<const N: usize> {
    pub(super) value: [f64; N],
    pub(super) rate: [f64; N],
    pub(super) acceleration: [f64; N],
    pub(super) integral: [f64; N],
}

#[derive(Clone, Copy, Debug)]
struct Node<const N: usize> {
    moment: [f64; N],
    // h² Qddot; no tiny timestep is squared in the interpolation solve.
    curvature: [f64; N],
    integral: [Dd; N],
}

/// Shared clamped C2 reconstruction; physical source units belong to its caller.
///
/// Preparation uses one uniform, positive tridiagonal factor for every source
/// and channel: O(sources*samples) work and storage, no dense solve.
/// Per-source queries use constant work and no allocations. Nodes retain the
/// original values and scaled curvatures; cubic coefficients are recovered
/// locally instead of storing a second full coefficient bank.
#[derive(Debug)]
pub(super) struct CubicHistory<const N: usize> {
    source_positions_m: Vec<[f64; 3]>,
    clock: QuadrupoleSampleClock,
    samples: usize,
    // Source-major nodes. Caller input is sample-major, like native records.
    nodes: Vec<Node<N>>,
    // Authored endpoint conditions remain distinct from solve roundoff in
    // reconstructed polynomial coefficients (needed for causal admission).
    initial_rates: Vec<[f64; N]>,
    workspace_bytes: usize,
    inverse_timestep_s: Dd,
    sixth: Dd,
    third: Dd,
}

impl<const N: usize> CubicHistory<N> {
    pub(super) fn required_workspace_bytes(
        sources: usize,
        samples: usize,
    ) -> Result<usize, AeroacError> {
        sources
            .checked_mul(samples)
            .and_then(|n| n.checked_mul(core::mem::size_of::<Node<N>>()))
            .and_then(|n| n.checked_add(sources.checked_mul(core::mem::size_of::<[f64; 3]>())?))
            .and_then(|n| n.checked_add(sources.checked_mul(core::mem::size_of::<[f64; N]>())?))
            .and_then(|n| {
                n.checked_add(
                    samples.checked_mul(
                        core::mem::size_of::<Dd>() + core::mem::size_of::<[Dd; N]>(),
                    )?,
                )
            })
            .ok_or_else(|| invalid("source history workspace extent overflow"))
    }

    /// The caller admits its physical layout before these indexed callbacks.
    /// All channels use the original stress owner's identical C2 algebra.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        source_positions_m: &[[f64; 3]],
        clock: QuadrupoleSampleClock,
        sample_times_s: &[f64],
        mut value_at: impl FnMut(usize, usize) -> [f64; N],
        mut end_at: impl FnMut(usize) -> CubicHistoryEnds<N>,
        budget: QuadrupoleHistoryBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        let sources = source_positions_m.len();
        let samples = sample_times_s.len();
        let count = sources.checked_mul(samples);
        if sources == 0
            || sources > budget.maximum_sources
            || samples < 2
            || samples > budget.maximum_samples_per_source
        {
            return Err(invalid("source history source/sample shape or budget"));
        }
        clock.validate_last_index(samples - 1)?;
        let h = clock.timestep_s();
        let inverse_h = 1.0 / h;
        if !inverse_h.is_finite()
            || !(inverse_h * inverse_h).is_finite()
            || inverse_h * inverse_h == 0.0
        {
            return Err(invalid(
                "source history derivative scale is unrepresentable",
            ));
        }
        let workspace_bytes = Self::required_workspace_bytes(sources, samples)?;
        if workspace_bytes > budget.maximum_workspace_bytes {
            return Err(invalid("source history exceeds its workspace byte budget"));
        }
        for (index, &time) in sample_times_s.iter().enumerate() {
            if index.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            clock.admit(index, time)?;
        }
        for (index, position) in source_positions_m.iter().enumerate() {
            let conditions = end_at(index);
            if index.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            if position
                .iter()
                .chain(&conditions.initial_rate)
                .chain(&conditions.final_rate)
                .chain(&conditions.initial_integral)
                .any(|v| !v.is_finite())
            {
                return Err(nonfinite("source history geometry or endpoint data"));
            }
        }
        let count = count.ok_or_else(|| invalid("source history node count overflow"))?;
        let mut nodes = reserve(count)?;
        let mut positions = reserve(sources)?;
        let mut initial_rates = reserve(sources)?;
        let mut inverse_pivots = reserve(samples)?;
        let mut work = reserve(samples)?;
        // Scaled-curvature equations have off-diagonals one, diagonal four
        // internally and two at either clamped endpoint. Their common SPD
        // tridiagonal factor is reused for all physical source channels.
        for i in 0..samples {
            if i.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            let diagonal = if i == 0 || i + 1 == samples { 2.0 } else { 4.0 };
            let previous = if i == 0 {
                Dd::ZERO
            } else {
                inverse_pivots[i - 1]
            };
            inverse_pivots.push(Dd::ONE / (Dd::from_f64(diagonal) - previous));
            work.push([Dd::ZERO; N]);
        }
        let h_dd = Dd::from_f64(h);
        let inverse_timestep_s = Dd::ONE / h_dd;
        let sixth = Dd::ONE / Dd::from_f64(6.0);
        let third = sixth * Dd::from_f64(2.0);
        for source in 0..sources {
            poll(&mut cancelled)?;
            positions.push(source_positions_m[source]);
            let ends = end_at(source);
            initial_rates.push(ends.initial_rate);
            let start = nodes.len();
            for i in 0..samples {
                if i.is_multiple_of(256) {
                    poll(&mut cancelled)?;
                }
                let moment = value_at(i, source);
                if moment.iter().any(|v| !v.is_finite()) {
                    return Err(nonfinite("source recorded stress moment"));
                }
                nodes.push(Node {
                    moment,
                    curvature: [0.0; N],
                    integral: [Dd::ZERO; N],
                });
            }
            for i in 0..samples {
                if i.is_multiple_of(256) {
                    poll(&mut cancelled)?;
                }
                for component in 0..N {
                    let q = Dd::from_f64(nodes[start + i].moment[component]);
                    let right = if i == 0 {
                        Dd::from_f64(nodes[start + 1].moment[component])
                            - q
                            - h_dd * Dd::from_f64(ends.initial_rate[component])
                    } else if i + 1 == samples {
                        h_dd * Dd::from_f64(ends.final_rate[component]) - q
                            + Dd::from_f64(nodes[start + i - 1].moment[component])
                    } else {
                        Dd::from_f64(nodes[start + i - 1].moment[component])
                            + Dd::from_f64(nodes[start + i + 1].moment[component])
                            - q * Dd::from_f64(2.0)
                    } * Dd::from_f64(6.0);
                    let previous = if i == 0 {
                        Dd::ZERO
                    } else {
                        work[i - 1][component]
                    };
                    work[i][component] = (right - previous) * inverse_pivots[i];
                }
            }
            for i in (0..samples).rev() {
                if i.is_multiple_of(256) {
                    poll(&mut cancelled)?;
                }
                for component in 0..N {
                    if i + 1 < samples {
                        work[i][component] =
                            work[i][component] - inverse_pivots[i] * work[i + 1][component];
                    }
                    if !finite_dd(work[i][component]) {
                        return Err(nonfinite("source history curvature solve"));
                    }
                    nodes[start + i].curvature[component] = work[i][component].to_f64();
                }
            }
            let mut integral = ends.initial_integral.map(Dd::from_f64);
            for i in 0..samples {
                if i.is_multiple_of(256) {
                    poll(&mut cancelled)?;
                }
                if integral.iter().any(|&v| !finite_dd(v)) {
                    return Err(nonfinite("source history integrated stress"));
                }
                nodes[start + i].integral = integral;
                if i + 1 < samples {
                    let coefficients =
                        coefficients(&nodes[start + i], &nodes[start + i + 1], sixth);
                    for component in 0..N {
                        integral[component] = integral[component]
                            + h_dd * primitive(coefficients[component], Dd::ONE, third);
                    }
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            source_positions_m: positions,
            clock,
            samples,
            nodes,
            initial_rates,
            workspace_bytes,
            inverse_timestep_s,
            sixth,
            third,
        })
    }

    pub(super) const fn clock(&self) -> QuadrupoleSampleClock {
        self.clock
    }
    pub(super) fn source_positions_m(&self) -> &[[f64; 3]] {
        &self.source_positions_m
    }
    pub(super) fn source_count(&self) -> usize {
        self.source_positions_m.len()
    }
    pub(super) const fn sample_count(&self) -> usize {
        self.samples
    }

    /// Original authored initial value, rate and integral. The acceleration
    /// remains that of the reconstructed spline. This preserves exact
    /// physical onset declarations separately from polynomial-solve roundoff;
    /// ordinary interpolation and its derivatives remain unchanged.
    pub(super) fn initial_source(
        &self,
        source: usize,
    ) -> Result<CubicSourceSample<N>, AeroacError> {
        if source >= self.source_count() {
            return Err(invalid("initial source index"));
        }
        let mut initial = self.source_at_elapsed(source, Dd::ZERO)?;
        initial.value = self.nodes[source * self.samples].moment;
        initial.rate = self.initial_rates[source];
        initial.integral = self.nodes[source * self.samples].integral.map(Dd::to_f64);
        Ok(initial)
    }
    /// Requested peak preparation payload, including temporary factors/work.
    pub(super) const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    /// Inclusive original source interval in double-double seconds.
    pub(super) fn source_time_bounds_s(&self) -> [Dd; 2] {
        let first = Dd::from_f64(self.clock.epoch_s());
        [first, first + self.duration_s()]
    }

    pub(super) fn duration_s(&self) -> Dd {
        Dd::from_f64(self.clock.timestep_s()) * Dd::from_f64((self.samples - 1) as f64)
    }

    /// Original source value, matching derivatives and cumulative integral.
    /// Keep a large absolute timestamp in DD until this method subtracts the
    /// original epoch. Query order and callback partitioning change no state.
    pub(super) fn source_at_clock_dd(
        &self,
        source: usize,
        time_s: Dd,
    ) -> Result<CubicSourceSample<N>, AeroacError> {
        let time_s = normalize_time(time_s)?;
        let [first, last] = self.source_time_bounds_s();
        if source >= self.source_count() || time_s.lt(first) || last.lt(time_s) {
            return Err(invalid("source index or uncovered finite history"));
        }
        self.source_at_elapsed(
            source,
            interval_offset(time_s, [first, last], self.duration_s())?,
        )
    }

    pub(super) fn source_at_elapsed(
        &self,
        source: usize,
        elapsed: Dd,
    ) -> Result<CubicSourceSample<N>, AeroacError> {
        if !finite_dd(elapsed) || elapsed.lt(Dd::ZERO) || self.duration_s().lt(elapsed) {
            return Err(invalid("source offset is outside recorded history"));
        }
        let h = Dd::from_f64(self.clock.timestep_s());
        let mut interval =
            ((elapsed * self.inverse_timestep_s).to_f64().floor() as usize).min(self.samples - 2);
        // A DD value just below a knot can round up in the integer estimate.
        if interval != 0 && elapsed.lt(h * Dd::from_f64(interval as f64)) {
            interval -= 1;
        }
        let u = (elapsed - h * Dd::from_f64(interval as f64)) * self.inverse_timestep_s;
        let left = &self.nodes[source * self.samples + interval];
        let right = &self.nodes[source * self.samples + interval + 1];
        let coefficients = coefficients(left, right, self.sixth);
        let inverse_h2 = self.inverse_timestep_s * self.inverse_timestep_s;
        let mut result = CubicSourceSample {
            value: [0.0; N],
            rate: [0.0; N],
            acceleration: [0.0; N],
            integral: [0.0; N],
        };
        for (j, [a, b, c, d]) in coefficients.into_iter().enumerate() {
            result.value[j] = (a + u * (b + u * (c + u * d))).to_f64();
            result.rate[j] = ((b + u * (c * Dd::from_f64(2.0) + u * d * Dd::from_f64(3.0)))
                * self.inverse_timestep_s)
                .to_f64();
            result.acceleration[j] =
                ((c * Dd::from_f64(2.0) + u * d * Dd::from_f64(6.0)) * inverse_h2).to_f64();
            result.integral[j] =
                (left.integral[j] + h * primitive([a, b, c, d], u, self.third)).to_f64();
        }
        if result
            .value
            .iter()
            .chain(&result.rate)
            .chain(&result.acceleration)
            .chain(&result.integral)
            .any(|v| !v.is_finite())
        {
            return Err(nonfinite("source reconstructed history"));
        }
        Ok(result)
    }
}

/// Immutable physical stress histories using the common clamped C2 core.
/// Original values, matching derivatives and cumulative integrals retain one
/// clock and geometry; preparation is O(sources*samples), evaluation O(1).
#[derive(Debug)]
pub struct PreparedQuadrupoleHistory {
    spline: CubicHistory<6>,
}

impl PreparedQuadrupoleHistory {
    /// Reconstruct the original sample-major physical stress record.
    /// The six entries are xx, yy, zz, xy, xz, yz [N m], with physical shear.
    /// Every endpoint rate and initial integral is explicitly supplied.
    /// This remains an offline interpolation model, not an antialias filter,
    /// extrapolator or a claim about source bandwidth or derivative accuracy.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        source_positions_m: &[[f64; 3]],
        clock: QuadrupoleSampleClock,
        sample_times_s: &[f64],
        stress_moments_n_m: &[[f64; 6]],
        ends: &[QuadrupoleHistoryEnds],
        budget: QuadrupoleHistoryBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        if ends.len() != source_positions_m.len()
            || source_positions_m.len().checked_mul(sample_times_s.len())
                != Some(stress_moments_n_m.len())
        {
            return Err(invalid("quadrupole history source/sample shape or budget"));
        }
        let sources = source_positions_m.len();
        let spline = CubicHistory::new(
            source_positions_m,
            clock,
            sample_times_s,
            |sample, source| stress_moments_n_m[sample * sources + source],
            |source| CubicHistoryEnds {
                initial_rate: ends[source].initial_rate_n_m_s,
                final_rate: ends[source].final_rate_n_m_s,
                initial_integral: ends[source].initial_integral_n_m_s,
            },
            budget,
            cancelled,
        )?;
        Ok(Self { spline })
    }
    pub const fn clock(&self) -> QuadrupoleSampleClock {
        self.spline.clock()
    }
    pub fn source_positions_m(&self) -> &[[f64; 3]] {
        self.spline.source_positions_m()
    }
    pub fn source_count(&self) -> usize {
        self.spline.source_count()
    }
    pub const fn sample_count(&self) -> usize {
        self.spline.sample_count()
    }
    /// Authored initial value/rate/integral, with reconstructed acceleration.
    pub fn initial_source(&self, source: usize) -> Result<TransientQuadrupoleSource, AeroacError> {
        self.spline.initial_source(source).map(stress_sample)
    }
    /// Peak preparation payload, including temporary factors and work.
    pub const fn workspace_bytes(&self) -> usize {
        self.spline.workspace_bytes()
    }
    pub fn source_time_bounds_s(&self) -> [Dd; 2] {
        self.spline.source_time_bounds_s()
    }
    fn duration_s(&self) -> Dd {
        self.spline.duration_s()
    }
    pub fn source_at_clock_dd(
        &self,
        source: usize,
        time_s: Dd,
    ) -> Result<TransientQuadrupoleSource, AeroacError> {
        self.spline
            .source_at_clock_dd(source, time_s)
            .map(stress_sample)
    }
    fn source_at_elapsed(
        &self,
        source: usize,
        elapsed: Dd,
    ) -> Result<TransientQuadrupoleSource, AeroacError> {
        self.spline
            .source_at_elapsed(source, elapsed)
            .map(stress_sample)
    }

    /// Prepare a free-space observer with this exact recorded source layout.
    pub fn prepare_observer(
        &self,
        acoustic: FluidAcousticState,
        receiver_m: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RecordedTransientQuadrupoleObserver<'_>, AeroacError> {
        let observer = PreparedTransientQuadrupoleObserver::new_cancellable(
            acoustic,
            receiver_m,
            self.source_positions_m(),
            self.source_count(),
            &mut cancelled,
        )?;
        self.bind_observer(observer, cancelled)
    }

    /// Bind an existing original-layout observer, including its ideal images.
    /// Every path needs covered history at the same receiver time. A record
    /// shorter than its complete echo-delay spread refuses; no echo is omitted.
    pub fn bind_observer(
        &self,
        observer: PreparedTransientQuadrupoleObserver,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RecordedTransientQuadrupoleObserver<'_>, AeroacError> {
        poll(&mut cancelled)?;
        if observer.source_count() != self.source_count() {
            return Err(invalid("transient observer changed recorded source count"));
        }
        let [first, last] = self.source_time_bounds_s();
        let mut earliest = first;
        let mut latest = Dd::from_f64(f64::INFINITY);
        for (i, &position) in self.source_positions_m().iter().enumerate() {
            if i.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            if observer.source_position_m(i) != Some(position) {
                return Err(invalid(
                    "transient observer changed recorded source geometry/order",
                ));
            }
        }
        for i in 0..observer.path_count() {
            if i.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            let delay = observer
                .path_delay_s(i)
                .ok_or_else(|| invalid("transient observer path is unavailable"))?;
            let lo = first + Dd::from_f64(delay);
            let hi = last + Dd::from_f64(delay);
            if !finite_dd(lo) || !finite_dd(hi) || !first.lt(lo) || !last.lt(hi) {
                return Err(invalid("recorded receiver interval is unresolved"));
            }
            if earliest.lt(lo) {
                earliest = lo;
            }
            if hi.lt(latest) {
                latest = hi;
            }
        }
        if latest.lt(earliest) {
            return Err(invalid(
                "finite stress record cannot cover every acoustic path",
            ));
        }
        poll(&mut cancelled)?;
        Ok(RecordedTransientQuadrupoleObserver {
            history: self,
            observer,
            receiver_bounds_s: [earliest, latest],
        })
    }
}

pub(super) fn stress_sample(sample: CubicSourceSample<6>) -> TransientQuadrupoleSource {
    TransientQuadrupoleSource {
        stress_moment_n_m: sample.value,
        stress_rate_n_m_s: sample.rate,
        stress_acceleration_n_m_s2: sample.acceleration,
        stress_integral_n_m_s: sample.integral,
    }
}

/// A full pressure/vector-velocity observer bound to one immutable record.
#[derive(Debug)]
pub struct RecordedTransientQuadrupoleObserver<'h> {
    history: &'h PreparedQuadrupoleHistory,
    observer: PreparedTransientQuadrupoleObserver,
    receiver_bounds_s: [Dd; 2],
}

impl RecordedTransientQuadrupoleObserver<'_> {
    pub fn history(&self) -> &PreparedQuadrupoleHistory {
        self.history
    }
    pub const fn propagation(&self) -> &PreparedTransientQuadrupoleObserver {
        &self.observer
    }
    /// Inclusive intersection of all original direct and reflected histories.
    pub const fn receiver_time_bounds_s(&self) -> [Dd; 2] {
        self.receiver_bounds_s
    }

    /// Query every original source at its own retarded clock before summing
    /// pressure and velocity and forming p*u. O(paths), allocation-free, with
    /// no material/source advancement, extrapolation or history reset.
    /// Clock subtraction is anchored at the nearer declared path endpoint;
    /// the reported inclusive endpoints map to the original source endpoints
    /// without adding a timing tolerance or clipping an uncovered query.
    pub fn evaluate_at_clock_dd(
        &self,
        receiver_time_s: Dd,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RetardedQuadrupoleField, AeroacError> {
        poll(&mut cancelled)?;
        let time = normalize_time(receiver_time_s)?;
        let [first, last] = self.receiver_bounds_s;
        if time.lt(first) || last.lt(time) {
            return Err(invalid("receiver time requires uncovered source history"));
        }
        let source_bounds = self.history.source_time_bounds_s();
        let duration = self.history.duration_s();
        self.observer.evaluate_sources(
            |_, source, delay| {
                let delay = Dd::from_f64(delay);
                let emission = time - delay;
                if !finite_dd(emission) || !emission.lt(time) {
                    return Err(invalid(
                        "transient quadrupole retarded time is unresolved on the DD clock",
                    ));
                }
                let path_bounds = [source_bounds[0] + delay, source_bounds[1] + delay];
                let elapsed = interval_offset(time, path_bounds, duration)?;
                self.history.source_at_elapsed(source, elapsed)
            },
            cancelled,
        )
    }
}

fn reserve<T>(count: usize) -> Result<Vec<T>, AeroacError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| invalid("quadrupole history allocation exceeds available memory"))?;
    Ok(values)
}

fn finite_dd(v: Dd) -> bool {
    v.hi.is_finite() && v.lo.is_finite()
}

fn normalize_time(t: Dd) -> Result<Dd, AeroacError> {
    if !finite_dd(t) {
        return Err(nonfinite("quadrupole history query clock"));
    }
    let t = Dd::from_pair(t.hi, t.lo);
    if !finite_dd(t) {
        return Err(nonfinite("quadrupole history query clock"));
    }
    Ok(t)
}

// Equivalent affine coordinates anchored at either end. Near a large epoch,
// (last + delay) - delay need not round back inside the original DD interval.
// Subtracting from the nearer declared arrival instead maps each endpoint
// exactly and retains a short interior offset. No out-of-range time is clamped.
fn interval_offset(time: Dd, [first, last]: [Dd; 2], duration: Dd) -> Result<Dd, AeroacError> {
    if time.lt(first) || last.lt(time) {
        return Err(invalid("receiver/source clock requires uncovered history"));
    }
    let left = time - first;
    let right = last - time;
    Ok(if left.lt(right) {
        left
    } else {
        duration - right
    })
}

fn coefficients<const N: usize>(left: &Node<N>, right: &Node<N>, sixth: Dd) -> [[Dd; 4]; N] {
    core::array::from_fn(|j| {
        let a = Dd::from_f64(left.moment[j]);
        let m0 = Dd::from_f64(left.curvature[j]);
        let m1 = Dd::from_f64(right.curvature[j]);
        [
            a,
            Dd::from_f64(right.moment[j]) - a - (m0 * Dd::from_f64(2.0) + m1) * sixth,
            m0 * Dd::from_f64(0.5),
            (m1 - m0) * sixth,
        ]
    })
}

fn primitive([a, b, c, d]: [Dd; 4], u: Dd, third: Dd) -> Dd {
    u * (a + u * (b * Dd::from_f64(0.5) + u * (c * third + u * d * Dd::from_f64(0.25))))
}
