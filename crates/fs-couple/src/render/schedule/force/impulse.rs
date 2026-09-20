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

/// Physical newton-second inputs projected through the performance's force ports.
pub mod ports {
    //! Physical momentum inputs using the SAME ports as the held-force performance.
    //!
    //! J_p [N s] projects as B J into the mass-normalized basis. The conjugate
    //! velocity is B^T v, so J_p^T B^T (v_before + v_after)/2 is precisely the
    //! kinetic-energy work of the jump, up to rounding. No mass, contact law or
    //! strike-velocity conversion is inferred. The caller supplies the impulse.

    use super::{ImpulseBudget, ImpulseRenderer, ModalImpulse, control, sizing};
    use super::super::{ForceRenderConfig, ModalForceEvent, ModalForceVoice, ScheduledRenderer, project};
    use crate::render::RenderError;

    /// A signed physical impulse at one declared force port.
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub struct PortImpulse {
        /// Apply before this absolute audio sample.
        pub sample: u64,
        /// Index in the supplied modal force voices.
        pub voice: usize,
        /// Index in that voice's admitted physical port matrix.
        pub port: usize,
        /// Applied momentum [N s], independent of callback length and sample rate.
        pub impulse_n_s: f64,
    }

    /// Bounds for offline physical-impulse lowering, separate from force budgets.
    #[derive(Clone, Copy, Debug)]
    pub struct PhysicalImpulseBudget {
        /// Maximum authored physical impulse events, including repeats and zeros.
        pub max_events: usize,
        /// Sum of modes times ports over all distinct (sample, voice) groups.
        pub max_projection_terms: usize,
        /// Sum of mode counts over those groups. Also bounds the scalar modal
        /// impulse history passed to the runtime before any zero coalescing.
        pub max_staged_modes: usize,
    }

    impl ScheduledRenderer {
        /// Build a physical performance with held loads AND instantaneous impulses.
        ///
        /// The same admitted port columns drive both inputs. At each timestamp,
        /// physical impulses sum per port in authored order, then project in mode
        /// order and port order. They never overwrite independently held forces.
        /// Static preload, when explicitly selected, precedes sample-zero impulses.
        /// All event/clock/budget admission is offline; no new stepping scheme or
        /// per-sample physical projection is introduced.
        pub fn from_modal_forces_and_impulses(
            voices: Vec<ModalForceVoice>,
            forces: Vec<ModalForceEvent>,
            mut impulses: Vec<PortImpulse>,
            force_config: ForceRenderConfig,
            impulse_budget: PhysicalImpulseBudget,
        ) -> Result<ImpulseRenderer, RenderError> {
            if impulses.len() > impulse_budget.max_events {
                return Err(sizing("physical impulse event budget exceeded"));
            }
            for input in &impulses {
                let voice = voices.get(input.voice)
                    .ok_or(RenderError::UnknownVoice { index: input.voice })?;
                if input.sample == u64::MAX || input.port >= voice.port_shapes.len()
                    || !input.impulse_n_s.is_finite()
                {
                    return Err(control("physical impulse needs a renderable sample, an existing port and finite newton-seconds"));
                }
            }
            impulses.sort_by_key(|input| (input.sample, input.voice));
            let mut terms = 0_usize;
            let mut modes = 0_usize;
            let mut previous = None;
            for input in &impulses {
                let key = (input.sample, input.voice);
                if previous == Some(key) { continue; }
                let voice = &voices[input.voice];
                let count = voice.model.modes().len();
                let work = count.checked_mul(voice.port_shapes.len())
                    .ok_or_else(|| sizing("physical impulse projection size overflow"))?;
                terms = terms.checked_add(work)
                    .ok_or_else(|| sizing("physical impulse projection work overflow"))?;
                modes = modes.checked_add(count)
                    .ok_or_else(|| sizing("physical impulse staging size overflow"))?;
                if terms > impulse_budget.max_projection_terms || modes > impulse_budget.max_staged_modes {
                    return Err(sizing("physical impulse projection or staging budget exceeded"));
                }
                previous = Some(key);
            }
            let mut lowered = Vec::new();
            lowered.try_reserve(modes).map_err(|_| sizing("cannot reserve projected physical impulses"))?;
            let mut index = 0;
            while index < impulses.len() {
                let first = impulses[index];
                let voice = &voices[first.voice];
                let mut port_impulses = vec![0.0; voice.port_shapes.len()];
                while index < impulses.len()
                    && (impulses[index].sample, impulses[index].voice) == (first.sample, first.voice)
                {
                    let input = impulses[index];
                    port_impulses[input.port] += input.impulse_n_s;
                    if !port_impulses[input.port].is_finite() {
                        return Err(control("summed physical impulse is nonfinite"));
                    }
                    index += 1;
                }
                // This is the existing physical port projection, not a second basis
                // convention. Its linear map applies equally to force or momentum.
                let generalized = project(&voice.port_shapes, &port_impulses)?;
                for (mode, impulse_n_s_per_sqrt_kg) in generalized.into_iter().enumerate() {
                    lowered.push(ModalImpulse {
                        sample: first.sample,
                        voice: first.voice,
                        mode,
                        impulse_n_s_per_sqrt_kg,
                    });
                }
            }
            let renderer = Self::from_modal_forces(voices, forces, force_config)?;
            renderer.with_modal_impulses(lowered, ImpulseBudget {
                max_events: impulse_budget.max_staged_modes,
                max_staged_modes: impulse_budget.max_staged_modes,
            })
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use super::super::super::ForceInitialization;
        use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
        use fs_math::c64::C64;

        fn config() -> ForceRenderConfig {
            ForceRenderConfig { sample_rate_hz: 48_000, max_block: 512,
                max_events: 20, max_controls: 100, max_projection_terms: 100 }
        }
        fn budget() -> PhysicalImpulseBudget {
            PhysicalImpulseBudget { max_events: 20, max_projection_terms: 100, max_staged_modes: 100 }
        }
        fn input(sample: u64, port: usize, impulse_n_s: f64) -> PortImpulse {
            PortImpulse { sample, voice: 0, port, impulse_n_s }
        }
        fn mass() -> ModalForceVoice {
            let model = ModalAcousticTimeModel::try_free_mass(48_000, 4.0, 0.0, 1.0,
                ModalAcousticTimeBudget::audible_reference()).unwrap();
            ModalForceVoice::new(model, vec![vec![0.5]], vec![0.0], ForceInitialization::RetainState).unwrap()
        }
        fn oscillator() -> ModalAcousticTimeModel {
            ModalAcousticTimeModel::try_new(48_000, vec![
                ModalAcousticMode { angular_frequency_rad_s: 300.0, damping_ratio: 0.02,
                    pressure_per_modal_velocity: C64 { re: 1.0, im: 0.1 } },
                ModalAcousticMode { angular_frequency_rad_s: 700.0, damping_ratio: 0.04,
                    pressure_per_modal_velocity: C64 { re: 0.5, im: 0.2 } },
            ], ModalAcousticTimeBudget::audible_reference()).unwrap()
        }

        #[test]
        fn newton_seconds_produce_the_declared_mass_momentum_and_work() {
            let mut performance = ScheduledRenderer::from_modal_forces_and_impulses(
                vec![mass()], Vec::new(), vec![input(0, 0, 2.0)], config(), budget()).unwrap();
            performance.block(&mut [0.0]).unwrap();
            let velocity = performance.modal_states(0).unwrap()[0].velocity_m_sqrt_kg_per_s / 2.0;
            assert_eq!(velocity, 1.5);
            assert_eq!(4.0 * (velocity - 1.0), 2.0);
            assert_eq!(performance.receipts()[0].kinetic_energy_change_j, 2.5);
        }

        #[test]
        fn port_impulses_and_held_loads_match_direct_modal_steps() {
            let columns = vec![vec![2.0, -1.0], vec![0.5, 3.0]];
            let impulses = vec![input(37, 0, -0.01), input(0, 0, 0.02), input(0, 1, 0.01),
                input(37, 0, 0.005), input(127, 1, 0.02)];
            let forces = vec![
                ModalForceEvent { sample: 53, voice: 0, port: 0, force_n: 0.0 },
                ModalForceEvent { sample: 53, voice: 0, port: 1, force_n: 0.1 },
            ];
            let mut reference = oscillator();
            let mut expected = Vec::new();
            for sample in 0..512 {
                let mut physical = [0.0; 2];
                for event in impulses.iter().filter(|event| event.sample == sample) {
                    physical[event.port] += event.impulse_n_s;
                }
                if physical != [0.0; 2] {
                    let mut states = reference.states().to_vec();
                    for mode in 0..2 {
                        states[mode].velocity_m_sqrt_kg_per_s += columns[0][mode] * physical[0]
                            + columns[1][mode] * physical[1];
                    }
                    reference.restore_states(&states).unwrap();
                }
                let physical_force = if sample < 53 { [0.03, -0.02] } else { [0.0, 0.1] };
                let force = [columns[0][0] * physical_force[0] + columns[1][0] * physical_force[1],
                    columns[0][1] * physical_force[0] + columns[1][1] * physical_force[1]];
                expected.push(reference.step(&force).unwrap().observer_pressure_pa);
            }
            for block in [1, 37, 128, 512] {
                let voice = ModalForceVoice::new(oscillator(), columns.clone(), vec![0.03, -0.02],
                    ForceInitialization::RetainState).unwrap();
                let mut performance = ScheduledRenderer::from_modal_forces_and_impulses(
                    vec![voice], forces.clone(), impulses.clone(), config(), budget()).unwrap();
                let mut actual = vec![0.0; 512];
                for chunk in actual.chunks_mut(block) { performance.block(chunk).unwrap(); }
                assert_eq!(actual, expected, "block {block}");
                assert_eq!(performance.receipts().len(), 3);
            }
        }

        #[test]
        fn same_time_physical_events_share_projection_and_staging_work() {
            let exact = PhysicalImpulseBudget { max_events: 2, max_projection_terms: 1, max_staged_modes: 1 };
            let mut performance = ScheduledRenderer::from_modal_forces_and_impulses(vec![mass()], Vec::new(),
                vec![input(0, 0, 6.0), input(0, 0, -4.0)], config(), exact).unwrap();
            performance.block(&mut [0.0]).unwrap();
            assert_eq!(performance.receipts()[0].kinetic_energy_change_j, 2.5);
            assert!(ScheduledRenderer::from_modal_forces_and_impulses(vec![mass()], Vec::new(),
                vec![input(0, 0, 1.0), input(1, 0, 1.0)], config(), exact).is_err());
        }

        #[test]
        fn physical_input_and_resource_refusals_are_explicit() {
            for event in [input(u64::MAX, 0, 1.0), input(0, 1, 1.0), input(0, 0, f64::NAN),
                PortImpulse { voice: 7, ..input(0, 0, 1.0) }]
            {
                assert!(ScheduledRenderer::from_modal_forces_and_impulses(vec![mass()], Vec::new(),
                    vec![event], config(), budget()).is_err());
            }
            for limits in [PhysicalImpulseBudget { max_events: 0, ..budget() },
                PhysicalImpulseBudget { max_projection_terms: 0, ..budget() },
                PhysicalImpulseBudget { max_staged_modes: 0, ..budget() }]
            {
                assert!(ScheduledRenderer::from_modal_forces_and_impulses(vec![mass()], Vec::new(),
                    vec![input(0, 0, 1.0)], config(), limits).is_err());
            }
            assert!(ScheduledRenderer::from_modal_forces_and_impulses(vec![mass()], Vec::new(),
                vec![input(0, 0, f64::MAX), input(0, 0, f64::MAX)], config(), budget()).is_err());
            let wrong_clock = ForceRenderConfig { sample_rate_hz: 44_100, ..config() };
            assert!(ScheduledRenderer::from_modal_forces_and_impulses(vec![mass()], Vec::new(),
                vec![input(0, 0, 1.0)], wrong_clock, budget()).is_err());
        }

        #[test]
        fn empty_impulses_preserve_the_existing_force_performance() {
            let mut plain = ScheduledRenderer::from_modal_forces(vec![mass()], Vec::new(), config()).unwrap();
            let mut wrapped = ScheduledRenderer::from_modal_forces_and_impulses(vec![mass()], Vec::new(),
                Vec::new(), config(), PhysicalImpulseBudget { max_events: 0,
                    max_projection_terms: 0, max_staged_modes: 0 }).unwrap();
            let mut a = [0.0; 64];
            let mut b = [0.0; 64];
            plain.block(&mut a).unwrap();
            wrapped.block(&mut b).unwrap();
            assert_eq!(a, b);
            assert!(wrapped.receipts().is_empty());
            assert_eq!(wrapped.samples_rendered(), plain.samples_rendered());
        }
    }
}

mod gated {
    //! Whole-callback cancellation without consuming boundary impulses.

    use super::{ImpulseRenderer, sizing};
    use crate::render::{GatedRenderOutcome, RenderError};
    use fs_exec::CancelGate;

    impl ImpulseRenderer {
        /// Render bounded host callbacks, polling cancellation before each callback.
        ///
        /// Completed callbacks and their impulse receipts remain accepted. Pending
        /// impulses and held controls at a cancelled boundary are untouched. Resume
        /// this same object with a fresh/unrequested gate; no event is reapplied.
        /// Internal event splits do not introduce additional cancellation points.
        /// The unused tail of `out` is never changed.
        ///
        /// # Errors
        /// Invalid output sizing or clock limits refuse before any state changes,
        /// including when cancellation was already requested. Physics refusals have
        /// the poisoning/discard semantics of `block`, not resumable cancellation.
        pub fn render_under_gate(
            &mut self,
            gate: &CancelGate,
            out: &mut [f64],
            block_len: usize,
            blocks: usize,
        ) -> Result<GatedRenderOutcome, RenderError> {
            self.renderer.context.validate_block_len(block_len)?;
            let required = block_len.checked_mul(blocks)
                .ok_or_else(|| sizing("impulse block_len * blocks overflows usize"))?;
            let samples = u64::try_from(required)
                .map_err(|_| sizing("impulse render length exceeds the sample clock"))?;
            let end = self.samples_rendered().checked_add(samples)
                .ok_or_else(|| sizing("impulse render would overflow the sample clock"))?;
            if out.len() < required {
                return Err(sizing("output must hold every requested impulse callback"));
            }
            let host_blocks = u64::try_from(blocks)
                .map_err(|_| sizing("impulse callback count exceeds the block clock"))?;
            self.validate_segment_budget(end, host_blocks)?;
            for index in 0..blocks {
                if gate.is_requested() {
                    return Ok(GatedRenderOutcome::Cancelled { blocks: index as u64 });
                }
                let start = index * block_len;
                self.block(&mut out[start..start + block_len])?;
            }
            Ok(GatedRenderOutcome::Completed { blocks: host_blocks })
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use super::super::{ImpulseBudget, ModalImpulse};
        use crate::render::{ControlDelta, ModalStringVoice, RenderContext, RenderVoice};
        use crate::render::schedule::{ScheduledControl, ScheduledRenderer};
        use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
        use fs_math::c64::C64;

        fn performance() -> ImpulseRenderer {
            let model = ModalAcousticTimeModel::try_new(48_000, vec![ModalAcousticMode {
                angular_frequency_rad_s: 500.0, damping_ratio: 0.02,
                pressure_per_modal_velocity: C64 { re: 1.0, im: 0.1 },
            }], ModalAcousticTimeBudget::audible_reference()).unwrap();
            let voice = RenderVoice::ModalString(ModalStringVoice::new(model, vec![0.0]).unwrap());
            ScheduledRenderer::new(RenderContext::new(vec![voice], 512), vec![ScheduledControl {
                sample: 37, delta: ControlDelta::SetModalForce { voice: 0, mode: 0, force_n_per_sqrt_kg: 0.1 },
            }], 1).unwrap().with_modal_impulses(vec![
                ModalImpulse { sample: 0, voice: 0, mode: 0, impulse_n_s_per_sqrt_kg: 0.5 },
                ModalImpulse { sample: 37, voice: 0, mode: 0, impulse_n_s_per_sqrt_kg: -0.1 },
            ], ImpulseBudget { max_events: 2, max_staged_modes: 2 }).unwrap()
        }

        #[test]
        fn cancellation_preserves_boundary_impulses_controls_and_exact_resume() {
            let mut reference = performance();
            let mut expected = [0.0; 111];
            reference.block(&mut expected).unwrap();
            let mut resumed = performance();
            let mut actual = [0.0; 111];
            resumed.block(&mut actual[..37]).unwrap();
            assert_eq!(resumed.receipts().len(), 1);
            assert_eq!(resumed.renderer().pending_controls().len(), 1);
            let gate = CancelGate::new();
            gate.request();
            let before = resumed.modal_states(0).unwrap().to_vec();
            let mut untouched = [123.0; 74];
            assert!(matches!(resumed.render_under_gate(&gate, &mut untouched, 37, 2).unwrap(),
                GatedRenderOutcome::Cancelled { blocks: 0 }));
            assert_eq!(untouched, [123.0; 74]);
            assert_eq!(resumed.modal_states(0).unwrap(), before.as_slice());
            assert_eq!(resumed.samples_rendered(), 37);
            assert_eq!(resumed.receipts().len(), 1);
            assert_eq!(resumed.renderer().pending_controls().len(), 1);
            assert!(matches!(resumed.render_under_gate(&CancelGate::new(), &mut actual[37..], 37, 2).unwrap(),
                GatedRenderOutcome::Completed { blocks: 2 }));
            assert_eq!(actual, expected);
            assert_eq!(resumed.receipts(), reference.receipts());
            assert_eq!(resumed.renderer().applied_controls(), reference.renderer().applied_controls());
        }

        #[test]
        fn whole_request_admission_precedes_cancellation_and_mutation() {
            let mut run = performance();
            let gate = CancelGate::new();
            gate.request();
            let mut out = [123.0; 20];
            assert!(matches!(run.render_under_gate(&gate, &mut out, 10, 3), Err(RenderError::Sizing { .. })));
            assert!(matches!(run.render_under_gate(&gate, &mut out, 2, usize::MAX), Err(RenderError::Sizing { .. })));
            assert!(matches!(run.render_under_gate(&gate, &mut out, 0, 1), Err(RenderError::EmptyBlock)));
            assert_eq!(out, [123.0; 20]);
            assert!(run.receipts().is_empty());
            assert_eq!(run.samples_rendered(), 0);
            assert!(matches!(run.render_under_gate(&gate, &mut out, 10, 0).unwrap(),
                GatedRenderOutcome::Completed { blocks: 0 }));
        }
    }
}
