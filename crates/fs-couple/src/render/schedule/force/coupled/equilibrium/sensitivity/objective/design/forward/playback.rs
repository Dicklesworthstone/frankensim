//! Fitted physical parameters -> settled load case -> scheduled transient motion.
//! Reuses candidate binding, the original preload, physical force compilation,
//! contact dynamics and observer transfers. No fitted transient response is implied.
use super::super::{DesignControl, DesignError, EquilibriumDesign, checkpoint};
use crate::render::{ControlDelta, RenderContext, RenderError, RenderVoice};
use crate::render::schedule::ScheduledRenderer;
use crate::render::schedule::force::{
    ForceInitialization, ForceRenderConfig, ModalForceEvent, ModalForceVoice,
};
use crate::render::schedule::force::coupled::render::{
    CoupledModalVoice, contact::multiple::MultiContactModalVoice,
};
use crate::render::schedule::force::coupled::contact::multiple::MultiContactModalSystem;
use fs_exec::CancelGate;

/// A held-force assignment to one load in the selected experiment. Two loads
/// at the same attachment remain independent; releasing one leaves the other.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaseForceEvent {
    /// Apply before this mechanics sample; zero is the first sample.
    pub sample: u64,
    /// Zero-based load index in the selected DesignLoadCase.
    pub load: usize,
    /// New signed physical force [N], not a modal force or an impulse.
    pub force_n: f64,
}

/// Explicit duration and scheduling budgets, separate from a static template's
/// unused audio settings. Mechanical and contact budgets remain unchanged.
#[derive(Clone, Copy, Debug)]
pub struct DesignPlaybackConfig {
    /// Requested mechanics samples. No padding, retiming or tail is inferred.
    pub samples: u64,
    /// Existing force compiler's clock, callback, event and projection budgets.
    pub force: ForceRenderConfig,
}

/// Original physics failures stay inspectable; no substitute sound is produced.
#[derive(Debug)]
pub enum DesignPlaybackError {
    /// Incomplete event, case, duration or scheduling admission.
    Invalid(&'static str),
    /// Original parameter/preload/work/cancellation refusal.
    Design(DesignError),
    /// Original force compiler, contact host or scheduler refusal.
    Render(RenderError),
}
impl core::fmt::Display for DesignPlaybackError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "design playback: {what}"),
            Self::Design(error) => write!(f, "design playback: {error}"),
            Self::Render(error) => write!(f, "design playback: {error}"),
        }
    }
}
impl std::error::Error for DesignPlaybackError {}
impl From<DesignError> for DesignPlaybackError {
    fn from(error: DesignError) -> Self { Self::Design(error) }
}
impl From<RenderError> for DesignPlaybackError {
    fn from(error: RenderError) -> Self { Self::Render(error) }
}

/// Frozen origin of the transient, distinct from transient validation.
#[derive(Clone, Debug, PartialEq)]
pub struct DesignPlaybackInfo {
    /// Selected independent experiment, in source order.
    pub case: usize,
    /// Actual bound parameters, including shared and case-specific assignments.
    pub physical_parameters: Vec<f64>,
    /// Mechanics rate. The wrapper does not resample acoustic observations.
    pub sample_rate_hz: u32,
    /// Requested horizon, including the last short callback.
    pub samples: u64,
    /// Loaded storage, including all spring/contact potentials [J].
    pub initial_energy_j: f64,
}

/// Owns the complete fitted mechanical state, forces and pending schedule.
/// The returned scheduler keeps its existing error/poisoning and cancellation
/// behavior. Callers own the output horizon and may use the shared WAV writer.
pub struct DesignPlayback {
    info: DesignPlaybackInfo,
    renderer: ScheduledRenderer,
}
impl DesignPlayback {
    /// Inspect the settled origin without advancing it.
    #[must_use]
    pub fn info(&self) -> &DesignPlaybackInfo { &self.info }
    /// Transfer the complete state to the existing block/audio API.
    #[must_use]
    pub fn into_renderer(self) -> ScheduledRenderer { self.renderer }
}

impl EquilibriumDesign {
    /// Settle one case at the supplied decision, then continue its actual state.
    ///
    /// Initial forces come from the PARAMETERIZED case, never the unmodified
    /// source load. Events replace individual physical loads without moving Q/V.
    /// Every body, its numerical ceilings, damping and pressure transfer remains
    /// in the original basis. An unobserved source stays acoustically silent.
    ///
    /// Candidate preparation uses the shared whole-family budget admission;
    /// one evaluation and ONE selected-case preload are charged. The caller must
    /// therefore have capacity for a complete family even though only one case
    /// is solved here. Failed work is not refunded; source templates never mutate.
    /// Playback is not another objective/adjoint evaluation and does not certify
    /// static design requirements throughout a transient or identify damping.
    ///
    /// All event and compiler admission precedes the preload. No partial renderer
    /// escapes on any error. Cancellation is polled by the original preload and
    /// at construction boundaries; rendering uses the existing callback gate.
    pub fn playback_case(&self, point: &[f64], case: usize, events: Vec<CaseForceEvent>,
        config: DesignPlaybackConfig, control: &mut DesignControl, gate: &CancelGate)
        -> Result<DesignPlayback, DesignPlaybackError>
    {
        checkpoint(gate)?;
        let Some(source_case) = self.cases.get(case) else {
            return Err(DesignPlaybackError::Invalid("unknown load case"));
        };
        if config.samples == 0 || config.samples > 28_800_000
            || !(1..=65_536).contains(&config.force.max_block)
            || config.force.sample_rate_hz == 0 || events.len() > config.force.max_events {
            return Err(DesignPlaybackError::Invalid("duration, callback, rate or event budget is invalid"));
        }
        for event in &events {
            if event.sample >= config.samples || event.load >= source_case.loads.len() || !event.force_n.is_finite() {
                return Err(DesignPlaybackError::Invalid("force event needs an in-window sample, an existing load and finite newtons"));
            }
        }
        let prepared = self.prepare(point, control, gate)?;
        let selected = &prepared.cases[case];
        let mut maps = vec![(0, 0); selected.loads.len()];
        let mut voices = Vec::with_capacity(self.models.len());
        // This existing compiler is used only to lower physical events. The
        // transient below hosts the ORIGINAL solved network, not these copies.
        for (component, model) in self.models.iter().enumerate() {
            checkpoint(gate)?;
            let mut columns = Vec::new();
            let mut initial = Vec::new();
            for (load, input) in selected.loads.iter().enumerate() {
                if input.attachment.component == component {
                    maps[load] = (component, columns.len());
                    columns.push(input.attachment.shapes.clone());
                    initial.push(input.force_n);
                }
            }
            // A component can have no actuator without disappearing from the
            // mechanical network. Its explicit zero port cannot receive events.
            if columns.is_empty() { columns.push(vec![0.0; model.modes().len()]); initial.push(0.0); }
            voices.push(ModalForceVoice::new(model.clone(), columns, initial, ForceInitialization::RetainState)?);
        }
        let projected = events.into_iter().map(|event| {
            let (voice, port) = maps[event.load];
            ModalForceEvent { sample: event.sample, voice, port, force_n: event.force_n }
        }).collect();
        let compiled = ScheduledRenderer::from_modal_forces(voices, projected, config.force)?;
        let mut controls = compiled.pending_controls().to_vec();
        drop(compiled);
        for event in &mut controls {
            let ControlDelta::SetModalForce { voice, mode, .. } = &mut event.delta else {
                return Err(DesignPlaybackError::Invalid("physical compiler returned a non-modal event"));
            };
            *mode += self.offsets[*voice];
            *voice = 0;
        }
        let solved = self.solve_case(&prepared, case, control, gate)?;
        let initial_energy_j = solved.stored_energy_j;
        let voice = if prepared.contacts.is_empty() {
            RenderVoice::CoupledModal(Box::new(CoupledModalVoice::new(solved.network, solved.external)?))
        } else {
            let system = MultiContactModalSystem::new(solved.network, prepared.contacts, self.budget.contact, gate)
                .map_err(RenderError::Coupled)?;
            RenderVoice::MultiContactModal(Box::new(MultiContactModalVoice::new(system, solved.external)?))
        };
        let context = RenderContext::new(vec![voice], config.force.max_block);
        let renderer = ScheduledRenderer::new(context, controls, config.force.max_controls)?;
        renderer.validate_sample_rate(config.force.sample_rate_hz)?;
        checkpoint(gate)?;
        Ok(DesignPlayback { renderer, info: DesignPlaybackInfo { case, physical_parameters: prepared.parameters,
            sample_rate_hz: config.force.sample_rate_hz, samples: config.samples, initial_energy_j } })
    }
}
