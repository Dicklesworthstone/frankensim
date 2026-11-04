//! Finite accepted compact-source histories observed on one microphone clock.
//!
//! Mechanics supplies Q, F and their physical time derivatives; it is never
//! resampled or advanced here. A C1 cubic Hermite interpolant represents each
//! source between accepted endpoints. Its analytic derivative supplies Qdot
//! and Fdot to the original moving-source pressure operator, at the SAME
//! retarded time. Uniform motion retains its analytic inverse; prepared curved
//! motion uses bounded inversion of the original C2 trajectory interpolant.
//! A moving microphone reads that same Eulerian pressure at its accepted
//! arrival-time position. Its velocity changes the observed source clock, not
//! the source's radiation coefficients.
//! Fixed laboratory sources in a uniform flow reuse this same history owner;
//! their original convected pressure paths retain a unit observed clock rate.
//!
//! Input bandwidth, interpolation convergence, finite-source compactness and
//! pressure linearity remain producer responsibilities. Declared receiver
//! bandwidth must include Doppler shift AND geometric amplitude modulation.
//! Hermite interpolation is not an antialias filter. Render adequately sampled
//! pressure first, then use the existing `DecimatedRenderer` when reducing its
//! rate. Neither missing prehistory nor a tail is synthesized.

use super::PressureRenderer;
use crate::render::RenderError;
use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::{
    CompactSourceStrength, MovingCompactObserver, PreparedMovingCompactSample, TimedCompactSource,
    UniformSourceMotion,
};
use fs_math::{dd::Dd, det};

mod flow;
mod motion;
mod receiver;
/// Native rigid-body attachments and bounded accepted acoustic histories.
pub mod rigid;
/// Full stationary compact fields and directional microphones on these histories.
pub mod stationary;
pub use motion::{CompactMotion, CompactSourceMotion};

/// Immutable accepted history and a prescribed trajectory in the medium frame.
#[derive(Clone, Copy, Debug)]
pub struct CompactSourceHistory<'a, M = UniformSourceMotion> {
    /// Source translation, independent of the mechanics' time integration.
    pub motion: M,
    /// Strictly increasing finite physical times with values AND derivatives.
    pub samples: &'a [TimedCompactSource],
    /// Declared source bandwidth [Hz], including physical source harmonics.
    /// Zero admits constant sources; it is not inferred from sparse endpoints.
    pub source_bandwidth_hz: f64,
}

/// One finite uniform receiver clock, before optional pressure decimation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingPressureConfig {
    /// Physical arrival time of output sample zero [s].
    pub first_receiver_time_s: f64,
    /// Pressure observation rate [Hz], independent of source solver steps.
    pub sample_rate_hz: u32,
    /// Complete finite sample window, with no implied ring-down after it.
    pub samples: u64,
    /// Largest output callback, preallocated at construction.
    pub max_block: usize,
    /// Declared received bandwidth [Hz], strictly below receiver Nyquist.
    /// Must include motion-induced amplitude modulation as well as Doppler.
    pub receiver_bandwidth_hz: f64,
    /// Absolute light-cone residual tolerance [s], at most a quarter sample.
    pub maximum_time_residual_s: f64,
}

/// Admission bounds for input processing, retained storage and lifetime work.
#[derive(Clone, Copy, Debug)]
pub struct MovingPressureBudget {
    /// Maximum number of coherently summed physical sources.
    pub maximum_sources: usize,
    /// Total accepted endpoints across all sources, including shared histories.
    pub maximum_source_samples: usize,
    /// Maximum complete receiver samples over this renderer's lifetime.
    pub maximum_output_samples: u64,
    /// Logical owned payload: track records, cubic coefficients and scratch.
    /// Borrowed input, allocator metadata and allocator rounding are excluded.
    pub maximum_payload_bytes: usize,
    /// Sources times output samples. Each uniform pair has one analytic root;
    /// each curved pair additionally retains its prepared trajectory's explicit
    /// root-iteration cap. All pairs have a bounded source-history search.
    pub maximum_source_receiver_pairs: u64,
}

#[derive(Clone, Copy)]
struct CubicInterval {
    // Linear, quadratic and cubic terms in normalized source time, for Q/Fxyz.
    // The constant term and exact endpoint derivatives stay in the source.
    terms: [[f64; 3]; 4],
}

struct PreparedTrack<'a, M> {
    source: CompactSourceHistory<'a, M>,
    first_interval: usize,
}

/// Coherent physical pressure from retained accepted histories, in pascals.
///
/// Preparation compiles all cubic coefficients once. Callbacks allocate
/// nothing and use global integer sample indices, so block partitioning does
/// not change timestamps or arithmetic. A failed/cancelled block leaves the
/// output and sample counter untouched; only private scratch may have changed.
/// Keep the borrowed source histories alive and immutable while observing.
pub struct MovingCompactRenderer<'a, M: CompactMotion = UniformSourceMotion> {
    observer: M::Observer,
    config: MovingPressureConfig,
    tracks: Vec<PreparedTrack<'a, M>>,
    intervals: Vec<CubicInterval>,
    scratch: Vec<f64>,
    completed: u64,
    maximum_retarded_iterations: u64,
}

impl<'a, M: CompactMotion> MovingCompactRenderer<'a, M> {
    /// Admit the complete receiver window and compile source interpolation.
    /// Cancellation is checked while traversing retained input, before output.
    pub fn new<Cancel: FnMut() -> bool>(
        observer: M::Observer,
        sources: &[CompactSourceHistory<'a, M>],
        config: MovingPressureConfig,
        budget: MovingPressureBudget,
        cancelled: &mut Cancel,
    ) -> Result<Self, RenderError> {
        check_cancel(cancelled)?;
        let rate = f64::from(config.sample_rate_hz);
        if sources.is_empty()
            || sources.len() > budget.maximum_sources
            || config.sample_rate_hz == 0
            || config.samples == 0
            || config.samples > budget.maximum_output_samples
            || config.samples > (1_u64 << 53)
            || config.max_block == 0
            || !config.first_receiver_time_s.is_finite()
            || !config.receiver_bandwidth_hz.is_finite()
            || config.receiver_bandwidth_hz < 0.0
            || config.receiver_bandwidth_hz >= 0.5 * rate
            || !config.maximum_time_residual_s.is_finite()
            || !(0.0..=0.25 / rate).contains(&config.maximum_time_residual_s)
        {
            return Err(sizing(
                "moving pressure requires finite clocks, bandwidth and explicit capacities",
            ));
        }
        let source_count =
            u64::try_from(sources.len()).map_err(|_| sizing("source count exceeds u64"))?;
        if source_count
            .checked_mul(config.samples)
            .is_none_or(|pairs| pairs > budget.maximum_source_receiver_pairs)
        {
            return Err(sizing(
                "moving pressure lifetime work exceeds its source/receiver budget",
            ));
        }
        let mut endpoint_count = 0_usize;
        let mut root_iterations_per_sample = 0_u64;
        for (index, source) in sources.iter().enumerate() {
            if index % 64 == 0 {
                check_cancel(cancelled)?;
            }
            if source.samples.len() < 2 {
                return Err(sizing(
                    "compact observation needs two accepted source endpoints",
                ));
            }
            source.motion.admit_observer(&observer)?;
            root_iterations_per_sample = root_iterations_per_sample
                .checked_add(
                    u64::try_from(source.motion.maximum_root_iterations())
                        .map_err(|_| sizing("trajectory root-iteration limit exceeds u64"))?,
                )
                .ok_or_else(|| sizing("combined trajectory root work overflows"))?;
            endpoint_count = endpoint_count
                .checked_add(source.samples.len())
                .ok_or_else(|| sizing("source endpoint count overflow"))?;
        }
        if endpoint_count > budget.maximum_source_samples {
            return Err(sizing("compact source histories exceed the input budget"));
        }
        let maximum_retarded_iterations = config
            .samples
            .checked_mul(root_iterations_per_sample)
            .ok_or_else(|| sizing("finite receiver-window root work overflows"))?;
        let interval_count = endpoint_count - sources.len();
        let payload = sources
            .len()
            .checked_mul(size_of::<PreparedTrack<'a, M>>())
            .and_then(|n| {
                interval_count
                    .checked_mul(size_of::<CubicInterval>())
                    .and_then(|m| n.checked_add(m))
            })
            .and_then(|n| {
                config
                    .max_block
                    .checked_mul(size_of::<f64>())
                    .and_then(|m| n.checked_add(m))
            })
            .ok_or_else(|| sizing("moving pressure payload overflow"))?;
        if payload > budget.maximum_payload_bytes {
            return Err(sizing(
                "moving pressure exceeds its retained payload budget",
            ));
        }
        let last_receiver_time = clock_time(config, config.samples - 1)?;
        let largest_time = config
            .first_receiver_time_s
            .abs()
            .max(last_receiver_time.abs());
        if config.samples > 1 && 1.0 / rate <= largest_time.next_up() - largest_time {
            return Err(sizing(
                "receiver sample spacing is unrepresentable at this time origin",
            ));
        }
        let mut tracks = Vec::new();
        let mut intervals = Vec::new();
        let mut scratch = Vec::new();
        tracks
            .try_reserve_exact(sources.len())
            .map_err(|_| sizing("cannot reserve moving source tracks"))?;
        intervals
            .try_reserve_exact(interval_count)
            .map_err(|_| sizing("cannot reserve source interpolation"))?;
        scratch
            .try_reserve_exact(config.max_block)
            .map_err(|_| sizing("cannot reserve pressure scratch"))?;
        scratch.resize(config.max_block, 0.0);
        for source in sources {
            check_cancel(cancelled)?;
            source.motion.admit_clock(
                config.first_receiver_time_s,
                last_receiver_time,
                config.maximum_time_residual_s,
            )?;
            let first_interval = intervals.len();
            if !source.source_bandwidth_hz.is_finite() || source.source_bandwidth_hz < 0.0 {
                return Err(sizing("source bandwidth must be finite and nonnegative"));
            }
            for (index, endpoint) in source.samples.iter().enumerate() {
                if index % 64 == 0 {
                    check_cancel(cancelled)?;
                }
                if !endpoint.source_time_s.is_finite() || !finite_strength(endpoint.strength) {
                    return Err(acoustic(AeroacError::NonFinite {
                        what: "accepted compact source endpoint",
                    }));
                }
            }
            for pair in source.samples.windows(2) {
                check_cancel(cancelled)?;
                let h = pair[1].source_time_s - pair[0].source_time_s;
                if !h.is_finite() || h <= 0.0 || source.source_bandwidth_hz * h >= 0.5 {
                    return Err(sizing(
                        "source times must increase and resolve the declared source bandwidth",
                    ));
                }
                intervals.push(compile_interval(pair[0].strength, pair[1].strength, h)?);
            }
            let first = prepare_history(
                &observer,
                source,
                config.first_receiver_time_s,
                config.maximum_time_residual_s,
                cancelled,
            )?;
            let last = prepare_history(
                &observer,
                source,
                last_receiver_time,
                config.maximum_time_residual_s,
                cancelled,
            )?;
            // Each motion owner supplies the continuous geometry admission and
            // a Doppler bound valid for the whole requested observation window.
            let maximum_doppler = source
                .motion
                .admit_window(&observer, &first, &last, cancelled)?;
            if source.source_bandwidth_hz * maximum_doppler > config.receiver_bandwidth_hz {
                return Err(sizing(
                    "declared receiver bandwidth does not include source Doppler expansion",
                ));
            }
            tracks.push(PreparedTrack {
                source: *source,
                first_interval,
            });
        }
        check_cancel(cancelled)?;
        Ok(Self {
            observer,
            config,
            tracks,
            intervals,
            scratch,
            completed: 0,
            maximum_retarded_iterations,
        })
    }

    /// Original receiver owner: fixed microphone or borrowed accepted path.
    pub const fn observer(&self) -> &M::Observer {
        &self.observer
    }
    /// Microphone position in the fluid rest frame at one admitted time.
    /// For uniform flow, the common receiver retains its fixed laboratory
    /// position separately; an unrepresentable frame transform may refuse.
    /// Limits belong to the first source; each other source retains its own
    /// prepared motion and clearance bounds.
    pub fn observer_at(&self, time: f64) -> Result<MovingCompactObserver, RenderError> {
        let last = clock_time(self.config, self.config.samples - 1)?;
        if !time.is_finite() || time < self.config.first_receiver_time_s || time > last {
            return Err(sizing(
                "microphone query exceeds the admitted receiver window",
            ));
        }
        self.tracks[0]
            .source
            .motion
            .observer_at(&self.observer, time)
    }
    /// Admitted physical receiver clock, horizon and declared bandwidth.
    pub const fn config(&self) -> MovingPressureConfig {
        self.config
    }
    /// Number of retained sources, in deterministic summation order.
    pub fn source_count(&self) -> usize {
        self.tracks.len()
    }
    /// Completed samples on the receiver clock.
    pub const fn samples_rendered(&self) -> u64 {
        self.completed
    }
    /// Unobserved samples inside the finite admitted window.
    pub fn remaining_samples(&self) -> u64 {
        self.config.samples - self.completed
    }
    /// Largest complete callback that fits the retained scratch buffer.
    pub const fn max_block_len(&self) -> usize {
        self.config.max_block
    }

    /// Checked upper bound on iterative root steps over the complete window.
    /// Uniform sources contribute zero; constant endpoint/search work is
    /// separately bounded by source count and retained history sizes.
    pub const fn maximum_retarded_iterations(&self) -> u64 {
        self.maximum_retarded_iterations
    }

    /// Exact uniform clock expression used by callbacks; no accumulated dt.
    pub fn receiver_time_at(&self, sample: u64) -> Result<f64, RenderError> {
        if sample >= self.config.samples {
            return Err(sizing(
                "receiver timestamp lies outside the admitted window",
            ));
        }
        clock_time(self.config, sample)
    }

    /// Render a complete transactional callback. The existing stream writer
    /// supplies cancellation between callbacks; direct callers can additionally
    /// cancel between samples/source groups with `block_controlled`.
    pub fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.block_controlled(output, &mut || false)
    }

    /// Observe a block with bounded cancellation polling and atomic publication.
    /// A refused block can be retried because accepted source data is immutable.
    pub fn block_controlled<Cancel: FnMut() -> bool>(
        &mut self,
        output: &mut [f64],
        cancelled: &mut Cancel,
    ) -> Result<(), RenderError> {
        if output.is_empty() {
            return Err(RenderError::EmptyBlock);
        }
        if output.len() > self.config.max_block {
            return Err(sizing(
                "moving pressure block exceeds its preallocated capacity",
            ));
        }
        let count =
            u64::try_from(output.len()).map_err(|_| sizing("pressure block exceeds u64"))?;
        self.validate_sample_count(count)?;
        for i in 0..output.len() {
            check_cancel(cancelled)?;
            let time = clock_time(self.config, self.completed + i as u64)?;
            let mut pressure = Dd::ZERO;
            for (j, track) in self.tracks.iter().enumerate() {
                if j % 64 == 0 {
                    check_cancel(cancelled)?;
                }
                let prepared = prepare_history(
                    &self.observer,
                    &track.source,
                    time,
                    self.config.maximum_time_residual_s,
                    cancelled,
                )?;
                let strength = interpolate(track, &self.intervals, prepared.source_time_s())?;
                pressure = pressure
                    + Dd::from_f64(prepared.evaluate(strength).map_err(acoustic)?.pressure_pa);
            }
            let value = pressure.to_f64();
            if !value.is_finite() {
                return Err(acoustic(AeroacError::NonFinite {
                    what: "coherent moving source pressure",
                }));
            }
            self.scratch[i] = value;
        }
        check_cancel(cancelled)?;
        output.copy_from_slice(&self.scratch[..output.len()]);
        self.completed += count;
        Ok(())
    }
}

impl<M: CompactMotion> PressureRenderer for MovingCompactRenderer<'_, M> {
    fn samples_rendered(&self) -> u64 {
        Self::samples_rendered(self)
    }
    fn max_block_len(&self) -> usize {
        Self::max_block_len(self)
    }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        if rate != self.config.sample_rate_hz {
            return Err(sizing("sink rate must match the moving receiver clock"));
        }
        Ok(())
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if samples > self.remaining_samples() {
            return Err(sizing(
                "pressure request exceeds retained source-history coverage",
            ));
        }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        Self::block(self, output)
    }
}

fn compile_interval(
    a: CompactSourceStrength,
    b: CompactSourceStrength,
    h: f64,
) -> Result<CubicInterval, RenderError> {
    let (a, da) = components(a);
    let (b, db) = components(b);
    let terms = core::array::from_fn(|i| {
        let delta = Dd::from_f64(b[i]) - Dd::from_f64(a[i]);
        let left = Dd::from_f64(h) * Dd::from_f64(da[i]);
        let right = Dd::from_f64(h) * Dd::from_f64(db[i]);
        [
            left.to_f64(),
            (Dd::from_f64(3.0) * delta - Dd::from_f64(2.0) * left - right).to_f64(),
            (left + right - Dd::from_f64(2.0) * delta).to_f64(),
        ]
    });
    if terms.iter().flatten().any(|v| !v.is_finite())
        || (0..4).any(|i| (da[i] != 0.0 && h * da[i] == 0.0) || (db[i] != 0.0 && h * db[i] == 0.0))
    {
        return Err(acoustic(AeroacError::InvalidParameter {
            what: "compact source Hermite coefficients are unrepresentable",
        }));
    }
    Ok(CubicInterval { terms })
}

fn interpolate<M>(
    track: &PreparedTrack<'_, M>,
    intervals: &[CubicInterval],
    time: f64,
) -> Result<CompactSourceStrength, RenderError> {
    let samples = track.source.samples;
    let upper = samples.partition_point(|sample| sample.source_time_s <= time);
    if upper == 0 || time > samples[samples.len() - 1].source_time_s {
        return Err(sizing(
            "retarded strength would extrapolate accepted source history",
        ));
    }
    let left = upper - 1;
    if time == samples[left].source_time_s {
        return Ok(samples[left].strength);
    }
    let h = samples[left + 1].source_time_s - samples[left].source_time_s;
    let u = Dd::from_f64((time - samples[left].source_time_s) / h);
    let (start, _) = components(samples[left].strength);
    let terms = intervals[track.first_interval + left].terms;
    let mut values = [0.0; 4];
    let mut slopes = [0.0; 4];
    for i in 0..4 {
        let [b, c, d] = terms[i].map(Dd::from_f64);
        values[i] = (Dd::from_f64(start[i]) + u * (b + u * (c + u * d))).to_f64();
        slopes[i] = ((b + u * (Dd::from_f64(2.0) * c + u * Dd::from_f64(3.0) * d))
            / Dd::from_f64(h))
        .to_f64();
    }
    let strength = CompactSourceStrength {
        volume_flow_m3_s: values[0],
        volume_flow_rate_m3_s2: slopes[0],
        force_on_fluid_n: [values[1], values[2], values[3]],
        force_rate_on_fluid_n_s: [slopes[1], slopes[2], slopes[3]],
    };
    if !finite_strength(strength) {
        return Err(acoustic(AeroacError::NonFinite {
            what: "interpolated compact source strength",
        }));
    }
    Ok(strength)
}

fn prepare_history<M: CompactMotion, Cancel: FnMut() -> bool>(
    observer: &M::Observer,
    source: &CompactSourceHistory<'_, M>,
    time: f64,
    tolerance: f64,
    cancelled: &mut Cancel,
) -> Result<PreparedMovingCompactSample, RenderError> {
    let first = source.samples[0].source_time_s;
    let last = source.samples[source.samples.len() - 1].source_time_s;
    let original_error = match source
        .motion
        .prepare_at(observer, time, tolerance, cancelled)
    {
        Ok(prepared) if (first..=last).contains(&prepared.source_time_s()) => return Ok(prepared),
        Ok(_) => sizing("receiver window is outside retained physical source history"),
        Err(error @ RenderError::Acoustic(AeroacError::Cancelled { .. })) => return Err(error),
        Err(error) => error,
    };
    // Analytic inversion can round just outside an exact retained endpoint.
    // Reuse that physical event ONLY when its rounded arrival is the requested
    // clock value and its measured light-cone residual meets the same budget.
    // This is no extrapolation, tolerance-sized time clamp, or invented history.
    for endpoint_time in [first, last] {
        check_cancel(cancelled)?;
        if let Ok(endpoint) = source
            .motion
            .prepare_endpoint(observer, endpoint_time, time)
            && endpoint.receiver_time_s() == time
            && endpoint.retarded_time_residual_s() <= tolerance
        {
            check_cancel(cancelled)?;
            return Ok(endpoint);
        }
    }
    Err(original_error)
}

fn admit_clearance(
    observer: &MovingCompactObserver,
    motion: &UniformSourceMotion,
    first: f64,
    last: f64,
) -> Result<(), RenderError> {
    let velocity = motion.velocity_m_s();
    let scale = velocity.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return Ok(());
    }
    if last < first {
        return Err(sizing(
            "uniform retarded source times lost their physical order",
        ));
    }
    let position = motion.reference_position_m();
    let elapsed = Dd::from_f64(first) - Dd::from_f64(motion.reference_time_s());
    let receiver = observer.observer_position_m();
    let mut offset = [Dd::ZERO; 3];
    let mut direction = [Dd::ZERO; 3];
    let mut projected = Dd::ZERO;
    let mut squared = Dd::ZERO;
    for i in 0..3 {
        offset[i] = Dd::from_f64(receiver[i])
            - Dd::from_f64(position[i])
            - Dd::from_f64(velocity[i]) * elapsed;
        direction[i] = Dd::from_f64(velocity[i]) / Dd::from_f64(scale);
        projected = projected + offset[i] * direction[i];
        squared = squared + direction[i] * direction[i];
    }
    // Minimize distance on the CONTINUOUS line segment. Rounding the closest
    // event back to a large f64 epoch can miss a collision between clock values.
    // Scaling velocity avoids squaring physical speeds; Dd retains the offset
    // from the original trajectory rather than a rounded source position.
    let extent = Dd::from_f64(scale) * (Dd::from_f64(last) - Dd::from_f64(first));
    let projected = projected / squared;
    let along = if projected.lt(Dd::ZERO) {
        Dd::ZERO
    } else if extent.lt(projected) {
        extent
    } else {
        projected
    };
    let displacement = core::array::from_fn(|i| (offset[i] - direction[i] * along).to_f64());
    let distance = norm(displacement);
    if !distance.is_finite() || distance < observer.limits().minimum_distance_m {
        return Err(sizing(
            "continuous source trajectory violates microphone clearance",
        ));
    }
    Ok(())
}

fn components(s: CompactSourceStrength) -> ([f64; 4], [f64; 4]) {
    (
        [
            s.volume_flow_m3_s,
            s.force_on_fluid_n[0],
            s.force_on_fluid_n[1],
            s.force_on_fluid_n[2],
        ],
        [
            s.volume_flow_rate_m3_s2,
            s.force_rate_on_fluid_n_s[0],
            s.force_rate_on_fluid_n_s[1],
            s.force_rate_on_fluid_n_s[2],
        ],
    )
}
fn finite_strength(s: CompactSourceStrength) -> bool {
    let (values, slopes) = components(s);
    values.iter().chain(&slopes).all(|v| v.is_finite())
}
fn clock_time(config: MovingPressureConfig, sample: u64) -> Result<f64, RenderError> {
    let time = (Dd::from_f64(config.first_receiver_time_s)
        + Dd::from_f64(sample as f64) / Dd::from_f64(f64::from(config.sample_rate_hz)))
    .to_f64();
    if !time.is_finite() {
        return Err(sizing("receiver clock overflow"));
    }
    Ok(time)
}
fn norm(v: [f64; 3]) -> f64 {
    det::hypot(det::hypot(v[0], v[1]), v[2])
}
fn sizing(what: &'static str) -> RenderError {
    RenderError::Sizing { what }
}
fn acoustic(error: AeroacError) -> RenderError {
    RenderError::Acoustic(error)
}
fn check_cancel(cancelled: &mut impl FnMut() -> bool) -> Result<(), RenderError> {
    if cancelled() {
        return Err(acoustic(AeroacError::Cancelled {
            what: "moving compact pressure observation",
        }));
    }
    Ok(())
}
