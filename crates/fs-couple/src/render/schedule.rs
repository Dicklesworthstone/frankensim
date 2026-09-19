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
            let until = self.events.get(self.next_event).map_or(end, |e| e.sample.min(end));
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
        let end = self.samples_rendered().checked_add(samples).ok_or(RenderError::Sizing {
            what: "scheduled render would overflow the sample clock",
        })?;
        if out.len() < required {
            return Err(RenderError::Sizing {
                what: "output slice must hold all requested scheduled callbacks",
            });
        }
        // At most one internal segment per sample. This up-front conservative
        // clock check excludes overflow after partially rendering a request.
        if self.context.blocks_rendered().checked_add(samples).is_none() {
            return Err(RenderError::Sizing {
                what: "scheduled render could overflow the internal block clock",
            });
        }
        self.validate_segment_budget(end)?;
        for index in 0..blocks {
            if gate.is_requested() {
                return Ok(GatedRenderOutcome::Cancelled { blocks: index as u64 });
            }
            let start = index * block_len;
            self.block(&mut out[start..start + block_len])?;
        }
        Ok(GatedRenderOutcome::Completed { blocks: blocks as u64 })
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
        if self.context.blocks_rendered().checked_add(segments).is_none() {
            return Err(RenderError::Sizing {
                what: "scheduled callback would overflow the internal block clock",
            });
        }
        Ok(())
    }
}
