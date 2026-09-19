//! Sample-accurate control scheduling (music bead 3ez8g.2.3).
//!
//! Events use the render context's integer sample clock and apply BEFORE the
//! named sample. A callback covers a half-open interval: an event at its end
//! belongs to the next nonempty callback. Simultaneous events retain input
//! order, so repeated assignments to one input are explicitly last-write-wins.
//! Only the loop partition changes; all integration remains in the existing
//! voices. No musical-note-to-force calibration or new contact law is implied.
//!
//! Admission sorts and validates the COMPLETE schedule and reserves its log
//! storage. Scheduling itself allocates nothing during rendering. Allocations
//! inside a voice retain that voice's disclosed allocation boundary.

use super::{ControlDelta, GatedRenderOutcome, RenderContext, RenderError};
use fs_exec::CancelGate;
use fs_scenario::gesture::{GestureSchedule, GestureTarget};

/// Lower one pressure track onto an audio sample clock before rendering.
///
/// Each control tick is sampled using the gesture's existing ramp semantics and
/// held until the next tick. Tick `k` applies at `ceil(k * sample_rate / control_rate)`,
/// never before its physical time. `samples` is a half-open horizon from sample
/// zero. The caller must use the same sample rate for the destination voice.
/// Unchanged values are coalesced; the initial value is always emitted for a
/// nonempty horizon. This does not reset vibration or infer pressure from notes.
///
/// `max_ticks` bounds evaluation work and storage, including unchanged ticks.
/// Admission may allocate; rendering the resulting schedule does not.
///
/// # Errors
/// Refuses absent/non-pressure tracks, zero rates, control rates above the audio
/// rate, horizons exceeding the tick budget, and allocation failure.
pub fn pressure_gesture_controls(
    schedule: &GestureSchedule,
    track_id: &str,
    voice: usize,
    sample_rate_hz: u32,
    samples: u64,
    max_ticks: usize,
) -> Result<Vec<ScheduledControl>, RenderError> {
    let control_rate = schedule.control_rate_hz;
    if sample_rate_hz == 0 || control_rate == 0 || control_rate > sample_rate_hz {
        return Err(RenderError::Control {
            what: "pressure gesture requires 0 < control rate <= audio sample rate",
        });
    }
    if !schedule
        .tracks()
        .iter()
        .any(|track| track.id == track_id && track.target == GestureTarget::BlowingPressure)
    {
        return Err(RenderError::Control {
            what: "pressure gesture binding requires an existing blowing-pressure track",
        });
    }
    // Wide integer arithmetic keeps long horizons and nonintegral rate ratios
    // independent of floating-point clock rounding.
    let rate = u128::from(sample_rate_hz);
    let control = u128::from(control_rate);
    let ticks = if samples == 0 {
        0
    } else {
        (u128::from(samples - 1) * control) / rate + 1
    };
    if ticks > max_ticks as u128 {
        return Err(RenderError::Sizing {
            what: "pressure gesture horizon exceeds the control-tick budget",
        });
    }
    let mut events = Vec::new();
    events
        .try_reserve(ticks as usize)
        .map_err(|_| RenderError::Sizing {
            what: "cannot reserve pressure gesture controls",
        })?;
    let mut previous = None;
    for tick in 0..ticks as u64 {
        let pressure_pa = schedule
            .sample(track_id, tick)
            .map_err(|_| RenderError::Control {
                what: "pressure gesture could not be sampled",
            })?;
        if previous != Some(pressure_pa) {
            events.push(ScheduledControl {
                sample: (u128::from(tick) * rate).div_ceil(control) as u64,
                delta: ControlDelta::SetBlowingPressure { voice, pressure_pa },
            });
            previous = Some(pressure_pa);
        }
    }
    Ok(events)
}

/// One input assignment on the context-relative, absolute sample clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScheduledControl {
    /// Apply immediately before rendering this sample (zero is the first).
    pub sample: u64,
    /// Existing voice input operation; stored vibration is not reset.
    pub delta: ControlDelta,
}

/// An admitted finite performance with exclusive ownership of its voice context.
///
/// Exclusive ownership prevents a caller from changing voice topology or
/// advancing the clock behind the schedule. Recover it with [`Self::into_context`].
pub struct ScheduledRenderer {
    context: RenderContext,
    events: Vec<ScheduledControl>,
    next_event: usize,
}

impl ScheduledRenderer {
    /// Admit a complete schedule before any controls or physical state move.
    ///
    /// Input may be unsorted. Equal-time events keep their original order.
    /// `max_events` is the caller's explicit storage/work admission budget.
    /// A previously advanced context is supported, but past events are refused.
    ///
    /// # Errors
    /// Invalid times, event budget, voice inputs, a poisoned context, zero block
    /// capacity, or failure to reserve the complete applied-control log.
    pub fn new(
        mut context: RenderContext,
        mut events: Vec<ScheduledControl>,
        max_events: usize,
    ) -> Result<Self, RenderError> {
        context.validate_controls(&[])?;
        if context.max_block_len() == 0 || events.len() > max_events {
            return Err(RenderError::Sizing {
                what: "schedule requires nonzero block capacity and events within max_events",
            });
        }
        for event in &events {
            if event.sample < context.samples_rendered() || event.sample == u64::MAX {
                return Err(RenderError::Control {
                    what: "scheduled sample is in the past or cannot be rendered before clock overflow",
                });
            }
            context.validate_controls(core::slice::from_ref(&event.delta))?;
        }
        // Stable sorting is deliberate: same-time assignments are ordered data.
        events.sort_by_key(|event| event.sample);
        context
            .controls_applied
            .try_reserve(events.len())
            .map_err(|_| RenderError::Sizing {
                what: "cannot reserve the admitted schedule's control log",
            })?;
        Ok(Self {
            context,
            events,
            next_event: 0,
        })
    }

    /// Completed sample clock, independent of callback size.
    #[must_use]
    pub const fn samples_rendered(&self) -> u64 {
        self.context.samples_rendered()
    }

    /// Read-only access to the hosted context and its internal-block log.
    #[must_use]
    pub const fn context(&self) -> &RenderContext {
        &self.context
    }

    /// Sample-addressed applied events, invariant under callback repartitioning.
    /// Unlike the context's block log, these timestamps name physical samples.
    #[must_use]
    pub fn applied_controls(&self) -> &[ScheduledControl] {
        &self.events[..self.next_event]
    }

    /// The remaining admitted events, in execution order.
    #[must_use]
    pub fn pending_controls(&self) -> &[ScheduledControl] {
        &self.events[self.next_event..]
    }

    /// Recover the context at its current physical state, dropping future events.
    #[must_use]
    pub fn into_context(self) -> RenderContext {
        self.context
    }

    /// Render one host callback, splitting only at actual event boundaries.
    ///
    /// The WHOLE request is size-checked before even a sample-zero control is
    /// applied. Events exactly at the callback end remain pending.
    ///
    /// # Errors
    /// Empty/oversized output or clock overflow refuses without mutation. A voice
    /// refusal poisons the context: discard the entire callback output; completed
    /// internal segments and controls are not rolled back or resumable.
    pub fn block(&mut self, out: &mut [f64]) -> Result<(), RenderError> {
        self.context.validate_block_len(out.len())?;
        let end = self.samples_rendered() + out.len() as u64;
        self.validate_segment_budget(end)?;
        let mut offset = 0;
        while offset < out.len() {
            let now = self.samples_rendered();
            while let Some(event) = self.events.get(self.next_event) {
                if event.sample != now {
                    break;
                }
                // Every delta was admitted; the owned voice topology cannot change.
                self.context
                    .apply_controls(core::slice::from_ref(&event.delta))?;
                self.next_event += 1;
            }
            let until = self
                .events
                .get(self.next_event)
                .map_or(end, |e| e.sample.min(end));
            // The difference is bounded by this already-admitted output slice.
            let len = (until - now) as usize;
            self.context.block(&mut out[offset..offset + len])?;
            offset += len;
        }
        Ok(())
    }

    /// Render fixed-size host callbacks, polling cancellation only BEFORE each
    /// host callback, not at the scheduler's internal event splits. An in-flight
    /// callback drains completely; an event at a cancelled boundary stays pending.
    /// The same object resumes with a fresh/unrequested gate and the next output
    /// slice, without reapplying any already-consumed event.
    ///
    /// # Errors
    /// Invalid output shape is refused up front, then the same errors as [`Self::block`].
    pub fn render_under_gate(
        &mut self,
        gate: &CancelGate,
        out: &mut [f64],
        block_len: usize,
        blocks: usize,
    ) -> Result<GatedRenderOutcome, RenderError> {
        self.context.validate_block_len(block_len)?;
        let required = block_len.checked_mul(blocks).ok_or(RenderError::Sizing {
            what: "scheduled block_len * blocks overflows usize",
        })?;
        let samples = u64::try_from(required).map_err(|_| RenderError::Sizing {
            what: "scheduled render length exceeds the sample clock",
        })?;
        let end = self
            .samples_rendered()
            .checked_add(samples)
            .ok_or(RenderError::Sizing {
                what: "scheduled render would overflow the sample clock",
            })?;
        if out.len() < required {
            return Err(RenderError::Sizing {
                what: "output slice must hold all requested scheduled callbacks",
            });
        }
        // At most one internal segment per sample. This up-front conservative
        // clock check excludes overflow after partially rendering a request.
        if self
            .context
            .blocks_rendered()
            .checked_add(samples)
            .is_none()
        {
            return Err(RenderError::Sizing {
                what: "scheduled render could overflow the internal block clock",
            });
        }
        self.validate_segment_budget(end)?;
        for index in 0..blocks {
            if gate.is_requested() {
                return Ok(GatedRenderOutcome::Cancelled {
                    blocks: index as u64,
                });
            }
            let start = index * block_len;
            self.block(&mut out[start..start + block_len])?;
        }
        Ok(GatedRenderOutcome::Completed {
            blocks: blocks as u64,
        })
    }

    fn validate_segment_budget(&self, end: u64) -> Result<(), RenderError> {
        let mut previous = self.samples_rendered();
        let mut segments = 1_u64;
        for event in self.pending_controls() {
            if event.sample >= end {
                break;
            }
            if event.sample != previous {
                segments = segments.checked_add(1).ok_or(RenderError::Sizing {
                    what: "scheduled segment count overflows",
                })?;
                previous = event.sample;
            }
        }
        if self
            .context
            .blocks_rendered()
            .checked_add(segments)
            .is_none()
        {
            return Err(RenderError::Sizing {
                what: "scheduled callback would overflow the internal block clock",
            });
        }
        Ok(())
    }
}

/// Explicit binding from a typed pressure track to a hosted reed voice.
///
/// No target is guessed from a track name. Every track must have exactly one
/// binding and each voice may have at most one pressure track in a compilation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PressureGestureBinding {
    /// Exact id in the source `GestureSchedule`.
    pub track: String,
    /// Voice slot validated later by [`ScheduledRenderer::new`].
    pub voice: usize,
}

/// Admission failures while lowering typed gestures to audio-sample controls.
#[derive(Debug)]
pub enum GestureCompileError {
    /// The source schedule refused a query.
    Gesture(fs_scenario::gesture::GestureError),
    /// A named track cannot be unambiguously rendered by this adapter.
    Track {
        /// Offending source track id.
        track: String,
        /// The unmet admission condition.
        what: &'static str,
    },
    /// A clock or resource request is invalid.
    Invalid {
        /// The unmet admission condition.
        what: &'static str,
    },
    /// Worst-case source-sampling work exceeds the caller's explicit budget.
    WorkBudget {
        /// Upper bound on track/event visits during sampling.
        required: u128,
        /// Caller-supplied maximum.
        allowed: u64,
    },
}

impl core::fmt::Display for GestureCompileError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Gesture(error) => write!(f, "gesture schedule: {error}"),
            Self::Track { track, what } => write!(f, "gesture track {track:?}: {what}"),
            Self::Invalid { what } => write!(f, "gesture compilation: {what}"),
            Self::WorkBudget { required, allowed } => write!(
                f,
                "gesture compilation needs at most {required} sampling visits; budget is {allowed}"
            ),
        }
    }
}

impl core::error::Error for GestureCompileError {}

/// Lower complete pressure performances through the existing gesture sampler.
///
/// Control tick `k` applies at `ceil(k * audio_rate_hz / control_rate_hz)`:
/// never earlier than its control-clock time, including non-divisor clocks.
/// Each value is held until the next tick. Only bitwise changes are emitted,
/// preserving the initial value at sample zero and avoiding repeated held-input
/// assignments. The half-open render interval is `[0, samples)`.
///
/// This is an OFFLINE adapter, not a second gesture interpolator or physical
/// solver. All tracks must be pressure tracks explicitly bound to distinct
/// voices; unsupported/unbound tracks refuse rather than silently disappear.
/// Overlapping ramps currently refuse because the source sampler does not yet
/// define interruption correctly. Completed voice state is never involved in
/// compilation; pass the result to [`ScheduledRenderer::new`] to validate voice
/// indices and kinds before rendering. The audio rate must match the voices.
///
/// `max_work` caps the worst-case track/event visits by the existing stateless
/// sampler, not just emitted changes (a long held track still costs work).
///
/// # Errors
/// Invalid clocks, bindings, unsupported tracks, overlapping ramps, non-finite
/// samples, work-budget exhaustion, or inability to allocate compiled controls.
pub fn compile_pressure_gestures(
    schedule: &fs_scenario::gesture::GestureSchedule,
    bindings: &[PressureGestureBinding],
    audio_rate_hz: u32,
    samples: u64,
    max_work: u64,
) -> Result<Vec<ScheduledControl>, GestureCompileError> {
    use fs_scenario::gesture::GestureTarget;
    use std::collections::BTreeSet;

    let rate = schedule.control_rate_hz;
    if rate == 0 || audio_rate_hz == 0 || rate > audio_rate_hz {
        return Err(GestureCompileError::Invalid {
            what: "control clock must be positive and no faster than the positive audio clock",
        });
    }
    let track_error = |id: &str, what| GestureCompileError::Track {
        track: id.to_string(),
        what,
    };
    let mut bound = BTreeSet::new();
    let mut voices = BTreeSet::new();
    let mut work_per_tick = 0_u128;
    for binding in bindings {
        let track = schedule
            .tracks()
            .iter()
            .find(|track| track.id == binding.track)
            .ok_or_else(|| track_error(&binding.track, "no source track with this id"))?;
        if !bound.insert(binding.track.as_str()) {
            return Err(track_error(&binding.track, "track is bound more than once"));
        }
        if !voices.insert(binding.voice) {
            return Err(track_error(&binding.track, "voice has more than one pressure track"));
        }
        if !matches!(track.target, GestureTarget::BlowingPressure) {
            return Err(track_error(&binding.track, "target is not a blowing-pressure input"));
        }
        for pair in track.events.windows(2) {
            // Subtraction avoids overflowing time + duration for finite inputs.
            if pair[0].transition_s > pair[1].time_s - pair[0].time_s {
                return Err(track_error(
                    &binding.track,
                    "overlapping ramps need source-sampler interruption support",
                ));
            }
        }
        // Each query may scan all track ids and all events on the selected track.
        work_per_tick += schedule.tracks().len() as u128 + track.events.len() as u128 + 1;
    }
    for track in schedule.tracks() {
        if !bound.contains(track.id.as_str()) {
            return Err(track_error(&track.id, "source track has no voice binding"));
        }
    }
    if samples == 0 || bindings.is_empty() {
        return Ok(Vec::new());
    }
    // ceil(k*a/c) < samples iff k*a <= (samples-1)*c. u128 keeps
    // u64 sample counts times u32 rates exact without a float clock.
    let ticks = ((u128::from(samples) - 1) * u128::from(rate))
        / u128::from(audio_rate_hz)
        + 1;
    let required = ticks.checked_mul(work_per_tick).ok_or(GestureCompileError::Invalid {
        what: "gesture sampling work bound overflows",
    })?;
    if required > u128::from(max_work) {
        return Err(GestureCompileError::WorkBudget { required, allowed: max_work });
    }
    let ticks = u64::try_from(ticks).map_err(|_| GestureCompileError::Invalid {
        what: "control tick count exceeds u64",
    })?;
    let mut previous = vec![None; bindings.len()];
    let mut controls = Vec::new();
    for tick in 0..ticks {
        let sample = (u128::from(tick) * u128::from(audio_rate_hz))
            .div_ceil(u128::from(rate));
        // By the tick bound this sample is strictly less than `samples` (u64).
        let sample = u64::try_from(sample).map_err(|_| GestureCompileError::Invalid {
            what: "compiled control sample exceeds u64",
        })?;
        for (index, binding) in bindings.iter().enumerate() {
            let pressure_pa = schedule
                .sample(&binding.track, tick)
                .map_err(GestureCompileError::Gesture)?;
            if !pressure_pa.is_finite() || pressure_pa < 0.0 {
                return Err(track_error(&binding.track, "sample is not finite nonnegative pressure"));
            }
            let bits = pressure_pa.to_bits();
            if previous[index] != Some(bits) {
                controls.try_reserve(1).map_err(|_| GestureCompileError::Invalid {
                    what: "cannot allocate compiled gesture controls",
                })?;
                controls.push(ScheduledControl {
                    sample,
                    delta: ControlDelta::SetBlowingPressure {
                        voice: binding.voice,
                        pressure_pa,
                    },
                });
                previous[index] = Some(bits);
            }
        }
    }
    Ok(controls)
}

#[cfg(test)]
mod pressure_gesture_tests {
    use super::*;
    use fs_scenario::gesture::{GestureEvent, GestureSchedule, GestureTarget, GestureTrack, GestureValue};

    fn track(id: &str) -> GestureTrack {
        GestureTrack {
            id: id.to_string(),
            target: GestureTarget::BlowingPressure,
            initial: GestureValue::PressurePa(0.0),
            events: vec![GestureEvent {
                time_s: 0.0,
                transition_s: 1.0,
                value: GestureValue::PressurePa(4.0),
            }],
        }
    }

    fn binding(id: &str, voice: usize) -> PressureGestureBinding {
        PressureGestureBinding { track: id.to_string(), voice }
    }

    #[test]
    fn ramps_use_the_source_clock_not_the_audio_buffer_partition() {
        let schedule = GestureSchedule::try_new(4, vec![track("blow")]).unwrap();
        let controls = compile_pressure_gestures(&schedule, &[binding("blow", 2)], 10, 10, 100).unwrap();
        assert_eq!(controls, vec![
            ScheduledControl { sample: 0, delta: ControlDelta::SetBlowingPressure { voice: 2, pressure_pa: 0.0 } },
            ScheduledControl { sample: 3, delta: ControlDelta::SetBlowingPressure { voice: 2, pressure_pa: 1.0 } },
            ScheduledControl { sample: 5, delta: ControlDelta::SetBlowingPressure { voice: 2, pressure_pa: 2.0 } },
            ScheduledControl { sample: 8, delta: ControlDelta::SetBlowingPressure { voice: 2, pressure_pa: 3.0 } },
        ]);
    }

    #[test]
    fn nondivisor_clock_and_half_open_end_never_apply_early() {
        let mut blow = track("blow");
        blow.events = vec![
            GestureEvent { time_s: 1.0 / 7.0, transition_s: 0.0, value: GestureValue::PressurePa(100.0) },
            GestureEvent { time_s: 2.0 / 7.0, transition_s: 0.0, value: GestureValue::PressurePa(0.0) },
        ];
        let schedule = GestureSchedule::try_new(7, vec![blow]).unwrap();
        let compile = |samples| compile_pressure_gestures(&schedule, &[binding("blow", 0)], 48_000, samples, 100).unwrap();
        let before = compile(13_715);
        assert_eq!(before.iter().map(|e| e.sample).collect::<Vec<_>>(), vec![0, 6_858]);
        let after = compile(13_716);
        assert_eq!(after.iter().map(|e| e.sample).collect::<Vec<_>>(), vec![0, 6_858, 13_715]);
    }

    #[test]
    fn holds_are_coalesced_and_every_voice_keeps_its_binding() {
        let mut a = track("a");
        a.events.clear();
        a.initial = GestureValue::PressurePa(1500.0);
        let mut b = a.clone();
        b.id = "b".to_string();
        let schedule = GestureSchedule::try_new(100, vec![a, b]).unwrap();
        let controls = compile_pressure_gestures(&schedule, &[binding("b", 7), binding("a", 3)], 48_000, 48_000, 1000).unwrap();
        assert_eq!(controls.len(), 2);
        assert_eq!(controls[0].delta, ControlDelta::SetBlowingPressure { voice: 7, pressure_pa: 1500.0 });
        assert_eq!(controls[1].delta, ControlDelta::SetBlowingPressure { voice: 3, pressure_pa: 1500.0 });
        assert!(controls.iter().all(|e| e.sample == 0));
    }

    #[test]
    fn complete_binding_and_target_admission_precedes_compilation() {
        let schedule = GestureSchedule::try_new(4, vec![track("a"), track("b")]).unwrap();
        for bindings in [
            vec![binding("a", 0)],
            vec![binding("a", 0), binding("missing", 1)],
            vec![binding("a", 0), binding("a", 1)],
            vec![binding("a", 0), binding("b", 0)],
        ] {
            assert!(matches!(compile_pressure_gestures(&schedule, &bindings, 10, 10, 1000), Err(GestureCompileError::Track { .. })));
        }
        let unsupported = GestureSchedule::try_new(4, vec![GestureTrack {
            id: "pedal".to_string(), target: GestureTarget::SustainPedal,
            initial: GestureValue::Fraction(0.0), events: Vec::new(),
        }]).unwrap();
        assert!(matches!(compile_pressure_gestures(&unsupported, &[binding("pedal", 0)], 10, 10, 1000), Err(GestureCompileError::Track { .. })));
    }

    #[test]
    fn work_budget_counts_sampling_not_only_emitted_changes() {
        let schedule = GestureSchedule::try_new(4, vec![track("blow")]).unwrap();
        // Four ticks * (one track lookup + one event + one initial value).
        assert!(matches!(compile_pressure_gestures(&schedule, &[binding("blow", 0)], 10, 10, 11),
            Err(GestureCompileError::WorkBudget { required: 12, allowed: 11 })));
        assert!(compile_pressure_gestures(&schedule, &[binding("blow", 0)], 10, 10, 12).is_ok());
        assert!(compile_pressure_gestures(&schedule, &[binding("blow", 0)], 10, u64::MAX, 12).is_err());
        assert!(compile_pressure_gestures(&schedule, &[binding("blow", 0)], 10, 0, 0).unwrap().is_empty());
    }

    #[test]
    fn invalid_clocks_and_unsupported_ramp_interruptions_refuse() {
        let mut schedule = GestureSchedule::try_new(4, vec![track("blow")]).unwrap();
        assert!(compile_pressure_gestures(&schedule, &[binding("blow", 0)], 0, 10, 100).is_err());
        assert!(compile_pressure_gestures(&schedule, &[binding("blow", 0)], 3, 10, 100).is_err());
        schedule.control_rate_hz = 0;
        assert!(compile_pressure_gestures(&schedule, &[binding("blow", 0)], 10, 10, 100).is_err());
        let mut blow = track("blow");
        blow.events.push(GestureEvent { time_s: 0.5, transition_s: 0.0, value: GestureValue::PressurePa(0.0) });
        let overlap = GestureSchedule::try_new(4, vec![blow]).unwrap();
        assert!(matches!(compile_pressure_gestures(&overlap, &[binding("blow", 0)], 10, 10, 100), Err(GestureCompileError::Track { .. })));
    }
}

/// Physical actuator histories projected through mass-normalized modal ports.
pub mod force;
