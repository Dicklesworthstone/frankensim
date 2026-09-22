//! Finite mixed performances assembled from the existing scheduled renderer.
//!
//! Each part contributes its current physical state and ALL pending controls.
//! The join is sample zero of a new output timeline: source clocks and modal /
//! contact / body state are not rewound. Pending outer controls are reindexed
//! and shifted by that part's completed sample count. A nested bow scheduler
//! retains its own absolute voice clock. No resampling or gain is introduced.
//!
//! Voices sum in part order, then original slot order. This is acoustic
//! superposition of independent parts, not a new mechanical interconnection.
//! Use an existing coupled-modal voice when parts exchange physical forces.

use crate::render::{ControlDelta, RenderContext, RenderError, RenderVoice};
use crate::render::schedule::{ScheduledControl, ScheduledRenderer};
use fs_exec::CancelGate;

/// Explicit finite output window and construction budgets.
#[derive(Clone, Copy, Debug)]
pub struct EnsembleConfig {
    /// Every part must already run at this exact audio rate [Hz].
    pub sample_rate_hz: u32,
    /// Largest output callback. May not exceed any source context's capacity.
    pub max_block: usize,
    /// Additional samples to render after joining the parts' current states.
    pub samples: u64,
    /// Maximum total voice slots; every supplied part must contain a voice.
    pub max_voices: usize,
    /// Maximum total pending OUTER controls. Nested bow schedules are already
    /// compiled under their own budgets and are moved, not recompiled or copied.
    pub max_events: usize,
}

/// Mapping from one input performance to its slots on the joined timeline.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EnsemblePartOrigin {
    /// The part's completed outer sample count at the join.
    pub source_sample: u64,
    /// First voice slot in the joined renderer.
    pub first_voice: usize,
    /// Number of consecutive slots contributed by this part.
    pub voices: usize,
}

/// Pressure playback over retained parts, bounded by an explicit finite window.
///
/// Prior outer control logs are not replayed or copied: their effects are in
/// the retained state. `origins` records the clock/slot translation; pending
/// controls, including those beyond the requested window, remain inspectable.
/// A bow's independent control history remains in its owned scheduler.
/// Construction may allocate; playback delegates to existing voice steppers
/// with their existing allocation and physical-model qualifications.
pub struct EnsembleRenderer {
    renderer: ScheduledRenderer,
    origins: Vec<EnsemblePartOrigin>,
    samples: u64,
}

fn sizing(what: &'static str) -> RenderError { RenderError::Sizing { what } }
fn control(what: &'static str) -> RenderError { RenderError::Control { what } }

impl EnsembleRenderer {
    /// Join independent, possibly already vibrating scheduled performances.
    ///
    /// The complete source topology, clock, pending-control and finite-window
    /// budgets are checked before moving inputs. A short requested window does
    /// not silently discard later controls. No part is advanced by admission.
    pub fn from_parts(parts: Vec<ScheduledRenderer>, config: EnsembleConfig)
        -> Result<Self, RenderError>
    {
        if parts.is_empty() || config.max_block == 0 || parts.len() > config.max_voices {
            return Err(sizing("ensemble needs nonempty parts, callback capacity and a voice budget"));
        }
        let mut voice_count = 0_usize;
        let mut event_count = 0_usize;
        for part in &parts {
            part.validate_sample_rate(config.sample_rate_hz)?;
            let count = part.context.voices.len();
            if count == 0 || config.max_block > part.context.max_block_len() {
                return Err(sizing("each ensemble part needs voices and the requested callback capacity"));
            }
            voice_count = voice_count.checked_add(count)
                .ok_or_else(|| sizing("ensemble voice count overflow"))?;
            event_count = event_count.checked_add(part.pending_controls().len())
                .ok_or_else(|| sizing("ensemble control count overflow"))?;
            if voice_count > config.max_voices || event_count > config.max_events {
                return Err(sizing("ensemble exceeds its voice or pending-control budget"));
            }
            part.samples_rendered().checked_add(config.samples)
                .ok_or_else(|| sizing("ensemble window would overflow a source sample clock"))?;
            for voice in &part.context.voices {
                if let RenderVoice::BowedString(bow) = voice {
                    if config.samples > bow.remaining_samples()
                        || config.max_block > bow.state().max_block_len()
                    {
                        return Err(sizing("ensemble window or callback exceeds a bowed part's admission"));
                    }
                }
            }
        }
        let mut voices = Vec::new();
        let mut events = Vec::new();
        let mut origins = Vec::new();
        let mut scratch = Vec::new();
        voices.try_reserve_exact(voice_count).map_err(|_| sizing("cannot reserve ensemble voices"))?;
        events.try_reserve_exact(event_count).map_err(|_| sizing("cannot reserve ensemble controls"))?;
        origins.try_reserve_exact(parts.len()).map_err(|_| sizing("cannot reserve ensemble origins"))?;
        scratch.try_reserve_exact(config.max_block).map_err(|_| sizing("cannot reserve ensemble scratch"))?;
        scratch.resize(config.max_block, 0.0);
        for part in parts {
            let source_sample = part.samples_rendered();
            let first_voice = voices.len();
            origins.push(EnsemblePartOrigin {
                source_sample, first_voice, voices: part.context.voices.len(),
            });
            for event in part.events.into_iter().skip(part.next_event) {
                let sample = event.sample.checked_sub(source_sample)
                    .ok_or_else(|| control("ensemble contains a pending event before its source clock"))?;
                let mut delta = event.delta;
                let slot = match &mut delta {
                    ControlDelta::SetModalForce { voice, .. }
                    | ControlDelta::SetPlateForce { voice, .. }
                    | ControlDelta::SetBlowingPressure { voice, .. } => voice,
                };
                *slot = (*slot).checked_add(first_voice)
                    .ok_or_else(|| sizing("ensemble voice index overflow"))?;
                events.push(ScheduledControl { sample, delta });
            }
            voices.extend(part.context.voices);
        }
        let context = RenderContext {
            voices, scratch, blocks_rendered: 0, samples_rendered: 0,
            controls_applied: Vec::new(), poisoned: false,
        };
        let renderer = ScheduledRenderer::new(context, events, config.max_events)?;
        Ok(Self { renderer, origins, samples: config.samples })
    }

    /// Source-clock and voice-slot mapping in original part order.
    #[must_use]
    pub fn origins(&self) -> &[EnsemblePartOrigin] { &self.origins }

    /// Completed samples on the joined timeline, independent of host callbacks.
    #[must_use]
    pub fn samples_rendered(&self) -> u64 { self.renderer.samples_rendered() }

    /// Unrendered samples in the admitted output window.
    #[must_use]
    pub fn remaining_samples(&self) -> u64 { self.samples - self.samples_rendered() }

    /// Inspect physical voices and pending/applied controls without advancing them.
    #[must_use]
    pub fn renderer(&self) -> &ScheduledRenderer { &self.renderer }

    /// Explicitly remove this output-window cap and recover the shared renderer.
    /// Bow horizons and poisoned-state restrictions still apply. This does not
    /// recover discarded pre-join outer logs or reset any physical voice state.
    #[must_use]
    pub fn into_renderer(self) -> ScheduledRenderer { self.renderer }

    fn validate_window(&self, samples: u64) -> Result<(), RenderError> {
        self.renderer.context.validate_controls(&[])?;
        if samples > self.remaining_samples() {
            return Err(sizing("output exceeds the admitted ensemble window"));
        }
        // Conservative: event/host splits cannot create more than one nonempty
        // internal segment per sample. Refuse overflow before any partial work.
        self.renderer.context.blocks_rendered().checked_add(samples)
            .ok_or_else(|| sizing("ensemble window could overflow the internal block clock"))?;
        Ok(())
    }

    /// Render physical observer pressure [Pa], preserving existing event splits.
    /// Empty/oversized or out-of-window requests refuse before state changes.
    /// Physics failures retain the shared renderer's discard/poison semantics.
    pub fn block(&mut self, out: &mut [f64]) -> Result<(), RenderError> {
        let samples = u64::try_from(out.len()).map_err(|_| sizing("ensemble output exceeds u64"))?;
        self.validate_window(samples)?;
        self.renderer.block(out)
    }

    /// Render a complete requested window with a possible short final callback.
    /// The whole window is admitted BEFORE any callback. Cancellation is polled
    /// only between host callbacks, returns the exact written prefix, and leaves
    /// every suffix slot and boundary event untouched for an exact resume.
    pub fn render_under_gate(&mut self, gate: &CancelGate, out: &mut [f64], block_len: usize)
        -> Result<EnsembleRenderOutcome, RenderError>
    {
        if block_len == 0 || block_len > self.renderer.context.max_block_len() {
            return Err(sizing("ensemble callback length is outside its admitted capacity"));
        }
        let requested = u64::try_from(out.len()).map_err(|_| sizing("ensemble output exceeds u64"))?;
        self.validate_window(requested)?;
        let mut samples = 0;
        for block in out.chunks_mut(block_len) {
            if gate.is_requested() { return Ok(EnsembleRenderOutcome::Cancelled { samples }); }
            self.block(block)?;
            samples += block.len();
        }
        Ok(EnsembleRenderOutcome::Completed { samples })
    }
}

/// Exact number of pressure samples written by one gated request.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EnsembleRenderOutcome {
    /// The complete requested window was rendered.
    Completed { /// Written prefix length.
        samples: usize },
    /// Cancellation stopped before the next callback; the same object resumes.
    Cancelled { /// Written prefix length.
        samples: usize },
}

