//! Sample-timed momentum inputs over the existing exact-ZOH modal renderer.
//!
//! A generalized impulse J [N s / sqrt(kg)] changes mass-normalized velocity
//! by J, with no displacement change and no elapsed physical time. It is not
//! a held force whose strength depends on callback length. Same-time inputs
//! add before admission; all affected voices commit together. The existing
//! scheduler still owns held controls, sample stepping, mixing and failures.
//!
//! This is an explicitly prescribed ideal impulse, not a hammer/contact law,
//! broadband radiation model, or inference from a musical strike velocity.
//! Only independent `ModalString` slots (generic mass-normalized models) are
//! admitted; coupled-network state must remain owned by its coupling solver.

use super::{ScheduledRenderer, control, sizing};
use crate::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeModel};
use crate::render::{RenderError, RenderVoice};

/// One additive mass-normalized momentum input on the absolute sample clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModalImpulse {
    /// Apply before this sample; an event at a callback's end stays pending.
    pub sample: u64,
    /// Independent modal voice slot in the hosted context.
    pub voice: usize,
    /// Coordinate in that voice's admitted modal basis.
    pub mode: usize,
    /// Signed generalized impulse [N s / sqrt(kg)].
    pub impulse_n_s_per_sqrt_kg: f64,
}

/// Offline compilation limits; callback storage never grows.
#[derive(Clone, Copy, Debug)]
pub struct ImpulseBudget {
    /// Maximum authored scalar modal inputs, before same-time aggregation.
    pub max_events: usize,
    /// Sum of mode counts over distinct (sample, voice) groups. Each admitted
    /// mode reserves a delta, a candidate state and a cloned model's mode/state
    /// records; exact-ZOH coefficients remain shared with the original model.
    pub max_staged_modes: usize,
}

/// Work at an accepted, instantaneous momentum transfer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImpulseReceipt {
    /// Sample boundary where the impulse was applied.
    pub sample: u64,
    /// Voice receiving the summed modal impulse.
    pub voice: usize,
    /// Kinetic-energy change of the represented velocity jump [J]. Position
    /// and elastic energy do not change; this is external work, not damping.
    pub kinetic_energy_change_j: f64,
}

struct StagedImpulse {
    sample: u64,
    voice: usize,
    delta: Vec<f64>,
    states: Vec<ModalAcousticState>,
    model: ModalAcousticTimeModel,
    work_j: f64,
}

/// Prescribed momentum events composed with an existing force performance.
///
/// Construction reserves all staging and receipts. Applying an impulse does
/// not allocate, reset ringing, change held forces or advance the sample clock.
/// The underlying voices retain their own disclosed allocation boundaries.
/// State/energy budgets are checked transactionally at each impulse; pressure
/// and ordinary time-step budgets remain checked by the existing voice step.
pub struct ImpulseRenderer {
    renderer: ScheduledRenderer,
    events: Vec<StagedImpulse>,
    next_event: usize,
    receipts: Vec<ImpulseReceipt>,
}

impl ScheduledRenderer {
    /// Add a complete, explicitly bounded momentum-input history.
    ///
    /// Inputs may be unsorted. Equal-time impulses on a coordinate are summed
    /// in input order, not last-write-wins. A partially rendered performance is
    /// supported, but past events and the unrenderable `u64::MAX` sample refuse.
    /// No state or held input changes during admission.
    ///
    /// # Errors
    /// Invalid event times, voice kinds, coordinates, nonfinite impulses or
    /// aggregate impulses, a poisoned renderer, or either compilation budget.
    pub fn with_modal_impulses(
        self,
        mut inputs: Vec<ModalImpulse>,
        budget: ImpulseBudget,
    ) -> Result<ImpulseRenderer, RenderError> {
        self.context.validate_controls(&[])?;
        if inputs.len() > budget.max_events {
            return Err(sizing("modal impulse event budget exceeded"));
        }
        for input in &inputs {
            if input.sample < self.samples_rendered() || input.sample == u64::MAX {
                return Err(control("modal impulse sample is past or unrenderable"));
            }
            let model = modal_model(&self, input.voice)?;
            if input.mode >= model.modes().len()
                || !input.impulse_n_s_per_sqrt_kg.is_finite()
            {
                return Err(control("modal impulse requires an existing coordinate and finite momentum"));
            }
        }
        // Stable ordering preserves the authored sum order within each group.
        inputs.sort_by_key(|input| (input.sample, input.voice));
        let mut required = 0_usize;
        let mut last = None;
        for input in &inputs {
            let key = (input.sample, input.voice);
            if last != Some(key) {
                required = required.checked_add(modal_model(&self, input.voice)?.modes().len())
                    .ok_or_else(|| sizing("modal impulse staging size overflow"))?;
                if required > budget.max_staged_modes {
                    return Err(sizing("modal impulse staging budget exceeded"));
                }
                last = Some(key);
            }
        }
        let mut events = Vec::new();
        events.try_reserve(inputs.len()).map_err(|_| sizing("cannot reserve modal impulse groups"))?;
        let mut index = 0;
        while index < inputs.len() {
            let first = inputs[index];
            let source = modal_model(&self, first.voice)?;
            let mut delta = vec![0.0; source.modes().len()];
            while index < inputs.len()
                && (inputs[index].sample, inputs[index].voice) == (first.sample, first.voice)
            {
                let input = inputs[index];
                delta[input.mode] += input.impulse_n_s_per_sqrt_kg;
                if !delta[input.mode].is_finite() {
                    return Err(control("summed modal impulse is nonfinite"));
                }
                index += 1;
            }
            events.push(StagedImpulse {
                sample: first.sample,
                voice: first.voice,
                delta,
                states: source.states().to_vec(),
                model: source.clone(),
                work_j: 0.0,
            });
        }
        let mut receipts = Vec::new();
        receipts.try_reserve(events.len()).map_err(|_| sizing("cannot reserve modal impulse receipts"))?;
        Ok(ImpulseRenderer { renderer: self, events, next_event: 0, receipts })
    }
}

fn modal_model(renderer: &ScheduledRenderer, voice: usize) -> Result<&ModalAcousticTimeModel, RenderError> {
    match renderer.context.voices.get(voice) {
        Some(RenderVoice::ModalString(voice)) => Ok(&voice.model),
        Some(_) => Err(control("prescribed modal impulses require an independent modal voice")),
        None => Err(RenderError::UnknownVoice { index: voice }),
    }
}

impl ImpulseRenderer {
    /// Physical sample clock, shared with the original performance.
    #[must_use]
    pub fn samples_rendered(&self) -> u64 { self.renderer.samples_rendered() }

    /// Read the original force/control performance without bypassing impulses.
    #[must_use]
    pub fn renderer(&self) -> &ScheduledRenderer { &self.renderer }

    /// Accepted impulse groups in sample/voice order, independent of callbacks.
    #[must_use]
    pub fn receipts(&self) -> &[ImpulseReceipt] { &self.receipts }

    /// Inspect an independent modal voice's accepted state without advancing it.
    pub fn modal_states(&self, voice: usize) -> Result<&[ModalAcousticState], RenderError> {
        Ok(modal_model(&self.renderer, voice)?.states())
    }

    /// Recover physical state and ordinary pending controls; future impulses
    /// are intentionally dropped. Use this only when ending/replacing a history.
    #[must_use]
    pub fn into_renderer(self) -> ScheduledRenderer { self.renderer }

    /// Render a host callback using the existing sample integrator and mixer.
    ///
    /// Empty/oversized output and clock exhaustion refuse before mutation.
    /// State-dependent impulse refusals poison the performance just like a
    /// voice refusal: discard that callback's output; it cannot be resumed.
    /// A failed same-time group commits none of its affected modal states.
    /// Held assignments at that sample still precede the first ensuing step;
    /// assignments do not themselves change the instantaneous modal state.
    pub fn block(&mut self, out: &mut [f64]) -> Result<(), RenderError> {
        self.renderer.context.validate_block_len(out.len())?;
        let end = self.samples_rendered().checked_add(out.len() as u64)
            .ok_or_else(|| sizing("impulse callback would overflow the sample clock"))?;
        self.validate_segment_budget(end, 1)?;
        let mut written = 0;
        while self.samples_rendered() < end {
            if self.next_event < self.events.len()
                && self.events[self.next_event].sample == self.samples_rendered()
            {
                if let Err(error) = self.apply_current_group() {
                    self.renderer.context.poisoned = true;
                    return Err(error);
                }
            }
            let next = self.events.get(self.next_event).map_or(end, |event| event.sample.min(end));
            let count = (next - self.samples_rendered()) as usize;
            self.renderer.block(&mut out[written..written + count])?;
            written += count;
        }
        Ok(())
    }

    fn validate_segment_budget(&self, end: u64, host_blocks: u64) -> Result<(), RenderError> {
        // A conservative bound over the UNION of both event streams. Counting
        // coincident events more than once is safe and avoids callback-time
        // allocation or a second timeline. Only internal block overflow uses it.
        let mut segments = host_blocks;
        for _ in self.renderer.pending_controls().iter().take_while(|event| event.sample < end) {
            segments = segments.checked_add(1).ok_or_else(|| sizing("impulse segment count overflow"))?;
        }
        for _ in self.events[self.next_event..].iter().take_while(|event| event.sample < end) {
            segments = segments.checked_add(1).ok_or_else(|| sizing("impulse segment count overflow"))?;
        }
        self.renderer.context.blocks_rendered().checked_add(segments)
            .ok_or_else(|| sizing("impulse performance would overflow the internal block clock"))?;
        Ok(())
    }

    fn apply_current_group(&mut self) -> Result<(), RenderError> {
        let begin = self.next_event;
        let sample = self.events[begin].sample;
        let mut end = begin;
        while end < self.events.len() && self.events[end].sample == sample { end += 1; }
        // All candidates are formed from CURRENT accepted state, not the state
        // at compilation. Nothing is published until every voice passes.
        for event in &mut self.events[begin..end] {
            let current = modal_model(&self.renderer, event.voice)?;
            event.states.clone_from_slice(current.states());
            event.work_j = 0.0;
            for ((candidate, old), &impulse) in event.states.iter_mut().zip(current.states()).zip(&event.delta) {
                let before = old.velocity_m_sqrt_kg_per_s;
                let after = before + impulse;
                candidate.velocity_m_sqrt_kg_per_s = after;
                event.work_j += (after - before) * f64::midpoint(before, after);
            }
            if !event.work_j.is_finite() {
                return Err(control("modal impulse work left the finite set"));
            }
            event.model.restore_states(&event.states).map_err(RenderError::Modal)?;
        }
        // Swaps cannot fail or allocate. Workspace and held forces stay with
        // each runtime; the cloned model has the identical admitted basis,
        // clock, limits and shared transition kernel.
        for event in &mut self.events[begin..end] {
            let RenderVoice::ModalString(voice) = &mut self.renderer.context.voices[event.voice] else {
                unreachable!("impulse voice kind was admitted under exclusive ownership")
            };
            core::mem::swap(&mut voice.model, &mut event.model);
            self.receipts.push(ImpulseReceipt {
                sample,
                voice: event.voice,
                kinetic_energy_change_j: event.work_j,
            });
        }
        self.next_event = end;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
    use crate::render::{ControlDelta, ModalStringVoice, RenderContext};
    use crate::render::schedule::ScheduledControl;
    use fs_math::c64::C64;

    fn model() -> ModalAcousticTimeModel {
        ModalAcousticTimeModel::try_new(48_000, vec![ModalAcousticMode {
            angular_frequency_rad_s: 700.0,
            damping_ratio: 0.02,
            pressure_per_modal_velocity: C64 { re: 2.0, im: 0.1 },
        }], ModalAcousticTimeBudget::audible_reference()).unwrap()
    }

    fn renderer(models: Vec<ModalAcousticTimeModel>, controls: Vec<ScheduledControl>) -> ScheduledRenderer {
        let voices = models.into_iter().map(|model| {
            let force = vec![0.0; model.modes().len()];
            RenderVoice::ModalString(ModalStringVoice::new(model, force).unwrap())
        }).collect();
        ScheduledRenderer::new(RenderContext::new(voices, 2048), controls, 20).unwrap()
    }

    fn impulse(sample: u64, voice: usize, value: f64) -> ModalImpulse {
        ModalImpulse { sample, voice, mode: 0, impulse_n_s_per_sqrt_kg: value }
    }

    fn budget() -> ImpulseBudget { ImpulseBudget { max_events: 20, max_staged_modes: 20 } }

    #[test]
    fn impulses_and_held_forces_match_direct_stepping_bitwise() {
        let inputs = vec![impulse(37, 0, 0.75), impulse(0, 0, 0.5), impulse(37, 0, -0.25), impulse(127, 0, -0.3)];
        let controls = vec![ScheduledControl { sample: 37, delta: ControlDelta::SetModalForce {
            voice: 0, mode: 0, force_n_per_sqrt_kg: 0.125,
        }}];
        let mut reference = model();
        let mut force = 0.0;
        let mut expected = Vec::new();
        for sample in 0..512 {
            let delta = match sample { 0 | 37 => 0.5, 127 => -0.3, _ => 0.0 };
            if delta != 0.0 {
                let mut state = reference.states().to_vec();
                state[0].velocity_m_sqrt_kg_per_s += delta;
                reference.restore_states(&state).unwrap();
            }
            if sample == 37 { force = 0.125; }
            expected.push(reference.step(&[force]).unwrap().observer_pressure_pa);
        }
        let mut first_receipts: Option<Vec<ImpulseReceipt>> = None;
        for block in [1, 37, 128, 512] {
            let mut performance = renderer(vec![model()], controls.clone())
                .with_modal_impulses(inputs.clone(), budget()).unwrap();
            let mut actual = vec![0.0; 512];
            for chunk in actual.chunks_mut(block) { performance.block(chunk).unwrap(); }
            assert_eq!(actual, expected, "callback size {block}");
            assert_eq!(performance.samples_rendered(), 512);
            assert_eq!(performance.receipts().len(), 3);
            if let Some(ref prior) = first_receipts { assert_eq!(performance.receipts(), prior.as_slice()); }
            else { first_receipts = Some(performance.receipts().to_vec()); }
        }
    }

    #[test]
    fn impulse_work_and_boundary_ownership_are_physical() {
        let free = ModalAcousticTimeModel::try_free_mass(48_000, 1.0, 0.0, 2.0,
            ModalAcousticTimeBudget::audible_reference()).unwrap();
        let mut performance = renderer(vec![free], Vec::new())
            .with_modal_impulses(vec![impulse(2, 0, -0.5)], budget()).unwrap();
        performance.block(&mut [0.0; 2]).unwrap();
        assert!(performance.receipts().is_empty());
        assert_eq!(performance.modal_states(0).unwrap()[0].velocity_m_sqrt_kg_per_s, 2.0);
        performance.block(&mut [0.0]).unwrap();
        assert_eq!(performance.modal_states(0).unwrap()[0].velocity_m_sqrt_kg_per_s, 1.5);
        assert_eq!(performance.receipts()[0].kinetic_energy_change_j, -0.875);
        assert_eq!(performance.receipts()[0].sample, 2);
    }

    #[test]
    fn simultaneous_impulses_are_additive_before_state_budget_checks() {
        let mut limits = ModalAcousticTimeBudget::audible_reference();
        limits.maximum_abs_velocity_m_sqrt_kg_per_s = 1.0;
        let free = ModalAcousticTimeModel::try_free_mass(48_000, 1.0, 0.0, 0.0, limits).unwrap();
        let mut performance = renderer(vec![free], Vec::new())
            .with_modal_impulses(vec![impulse(0, 0, 2.0), impulse(0, 0, -1.5)], budget()).unwrap();
        performance.block(&mut [0.0]).unwrap();
        assert_eq!(performance.modal_states(0).unwrap()[0].velocity_m_sqrt_kg_per_s, 0.5);
        assert_eq!(performance.receipts()[0].kinetic_energy_change_j, 0.125);
    }

    #[test]
    fn failed_same_time_batch_changes_no_voice_and_poisons_the_callback() {
        let mut limits = ModalAcousticTimeBudget::audible_reference();
        limits.maximum_abs_velocity_m_sqrt_kg_per_s = 1.0;
        let bounded = ModalAcousticTimeModel::try_free_mass(48_000, 1.0, 0.0, 0.0, limits).unwrap();
        let mut performance = renderer(vec![model(), bounded], Vec::new())
            .with_modal_impulses(vec![impulse(0, 0, 0.25), impulse(0, 1, 2.0)], budget()).unwrap();
        let before = performance.modal_states(0).unwrap().to_vec();
        let mut out = [123.0];
        assert!(matches!(performance.block(&mut out), Err(RenderError::Modal(_))));
        assert_eq!(performance.modal_states(0).unwrap(), before.as_slice());
        assert_eq!(performance.samples_rendered(), 0);
        assert!(performance.receipts().is_empty());
        assert_eq!(out, [123.0]);
        assert!(matches!(performance.block(&mut out), Err(RenderError::Poisoned)));
    }

    #[test]
    fn malformed_future_events_and_budgets_refuse_before_rendering() {
        for event in [impulse(u64::MAX, 0, 1.0), impulse(1000, 9, 1.0), impulse(1000, 0, f64::NAN),
            ModalImpulse { mode: 1, ..impulse(1000, 0, 1.0) }]
        {
            assert!(renderer(vec![model()], Vec::new()).with_modal_impulses(vec![event], budget()).is_err());
        }
        assert!(renderer(vec![model()], Vec::new()).with_modal_impulses(vec![impulse(0, 0, 1.0)],
            ImpulseBudget { max_events: 0, max_staged_modes: 20 }).is_err());
        assert!(renderer(vec![model()], Vec::new()).with_modal_impulses(vec![impulse(0, 0, 1.0)],
            ImpulseBudget { max_events: 20, max_staged_modes: 0 }).is_err());
        assert!(renderer(vec![model()], Vec::new()).with_modal_impulses(
            vec![impulse(0, 0, f64::MAX), impulse(0, 0, f64::MAX)], budget()).is_err());
        let mut advanced = renderer(vec![model()], Vec::new());
        advanced.block(&mut [0.0; 2]).unwrap();
        assert!(advanced.with_modal_impulses(vec![impulse(1, 0, 1.0)], budget()).is_err());
    }

    #[test]
    fn oversized_requests_leave_sample_zero_impulses_pending() {
        let mut performance = renderer(vec![model()], Vec::new())
            .with_modal_impulses(vec![impulse(0, 0, 0.5)], budget()).unwrap();
        assert!(matches!(performance.block(&mut []), Err(RenderError::EmptyBlock)));
        assert!(matches!(performance.block(&mut [0.0; 2049]), Err(RenderError::Sizing { .. })));
        assert_eq!(performance.samples_rendered(), 0);
        assert!(performance.receipts().is_empty());
        performance.block(&mut [0.0]).unwrap();
        assert_eq!(performance.receipts()[0].sample, 0);
    }
}
