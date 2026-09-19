//! Existing physical-force compiler and event scheduler hosting coupled mechanics.
use super::{CoupledModalSystem, ModalConnection, ModalCouplingConfig};
use super::super::{ForceInitialization, ForceRenderConfig, ModalForceEvent, ModalForceVoice};
use crate::render::{ControlDelta, RenderContext, RenderError, RenderVoice};
use crate::render::schedule::ScheduledRenderer;
use fs_exec::CancelGate;

/// One network in a render slot. Its mode controls use component order, then
/// original mode order; connection reactions remain internal and are never
/// overwritten by external actuator assignments.
pub struct CoupledModalVoice {
    system: CoupledModalSystem,
    held_force: Vec<f64>,
    poisoned: bool,
}
impl CoupledModalVoice {
    /// Wrap an admitted network without changing its accepted physical state.
    pub fn new(system: CoupledModalSystem, held_force: Vec<f64>) -> Result<Self, RenderError> {
        if held_force.len() != system.mode_count() || held_force.iter().any(|f| !f.is_finite()) {
            return Err(invalid("coupled voice requires finite external forces for every retained mode"));
        }
        Ok(Self { system, held_force, poisoned: false })
    }
    /// Original components and the last complete coupled energy/reaction record.
    #[must_use]
    pub const fn system(&self) -> &CoupledModalSystem { &self.system }
    /// Shared mechanical period [s], never resampled by the wrapper.
    #[must_use]
    pub fn sample_period_s(&self) -> f64 { self.system.sample_period_s() }
    /// Validate a held external-mode control without touching any input or state.
    pub fn validate_force(&self, mode: usize, force: f64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if mode >= self.held_force.len() || !force.is_finite() {
            return Err(invalid("coupled mode control is outside the basis or not finite"));
        }
        Ok(())
    }
    pub(in crate::render) fn set_force_admitted(&mut self, mode: usize, force: f64) {
        self.held_force[mode] = force;
    }
    /// Advance complete network samples. A physics error may leave an earlier
    /// callback prefix, so this wrapper poisons itself just like RenderContext.
    /// Direct CoupledModalSystem trials remain independently transactional.
    pub fn step_block(&mut self, out: &mut [f64]) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if out.is_empty() { return Err(RenderError::EmptyBlock); }
        let count = u64::try_from(out.len()).map_err(|_| invalid("coupled block length exceeds the clock"))?;
        self.system.samples_rendered().checked_add(count).ok_or_else(|| invalid("coupled sample clock overflow"))?;
        for sample in out {
            match self.system.step(&self.held_force) {
                Ok(frame) => *sample = frame.observer_pressure_pa,
                Err(error) => { self.poisoned = true; return Err(RenderError::Coupled(error)); }
            }
        }
        Ok(())
    }
}

impl ScheduledRenderer {
    /// Compile existing physical actuator schedules, then join their actual
    /// mechanical components rather than summing independent pressure tracks.
    /// The force projector, same-time ordering, event/log budgets and callback
    /// cancellation are exactly the existing scheduler's. The new operation is
    /// only the two-way spring/damper connection at each mechanical sample.
    ///
    /// Either every component retains its supplied state, or every component
    /// requests static preload with zero initial Q/V. Preload solves the full
    /// network equilibrium; mixed preload/state requests refuse instead of
    /// inventing an inconsistent partly settled assembly.
    pub fn from_coupled_modal_forces(
        mut voices: Vec<ModalForceVoice>,
        events: Vec<ModalForceEvent>,
        force_config: ForceRenderConfig,
        connections: Vec<ModalConnection>,
        coupling: ModalCouplingConfig,
        gate: &CancelGate,
    ) -> Result<Self, RenderError> {
        super::poll(Some(gate)).map_err(RenderError::Coupled)?;
        super::validate_config(coupling).map_err(RenderError::Coupled)?;
        let count = voices.iter().try_fold(0_usize, |n,v| n.checked_add(v.model.modes().len()))
            .ok_or_else(|| invalid("coupled mode count overflow"))?;
        if count > coupling.max_modes || connections.len() > coupling.max_connections {
            return Err(invalid("coupled component/connection counts exceed declared limits"));
        }
        let preload = voices.iter().any(|v| v.initialization == ForceInitialization::StaticPreload);
        if preload {
            for voice in &mut voices {
                if voice.initialization != ForceInitialization::StaticPreload
                    || voice.model.states().iter().any(|s| s.displacement_m_sqrt_kg != 0.0 || s.velocity_m_sqrt_kg_per_s != 0.0) {
                    return Err(invalid("coupled preload must cover all components with zero input vibration"));
                }
                // Defer initialization until the full spring system exists.
                voice.initialization = ForceInitialization::RetainState;
            }
        }
        let prepared = Self::from_modal_forces(voices, events, force_config)?;
        super::poll(Some(gate)).map_err(RenderError::Coupled)?;
        let mut components = Vec::new();
        let mut offsets = Vec::new();
        let mut held = Vec::with_capacity(count);
        for voice in prepared.context.voices {
            let RenderVoice::ModalString(voice) = voice else {
                return Err(invalid("modal force compiler returned a non-modal component"));
            };
            offsets.push(held.len());
            held.extend_from_slice(&voice.held_force);
            components.push(voice.model);
        }
        let mut system = CoupledModalSystem::new(components, connections, coupling, gate)
            .map_err(RenderError::Coupled)?;
        if preload { system.initialize_static_equilibrium(&held, gate).map_err(RenderError::Coupled)?; }
        let mut controls = prepared.events;
        for event in &mut controls {
            let ControlDelta::SetModalForce { voice, mode, .. } = &mut event.delta else {
                return Err(invalid("modal force compiler returned a non-modal assignment"));
            };
            *mode += offsets[*voice];
            *voice = 0;
        }
        let voice = CoupledModalVoice::new(system, held)?;
        let context = RenderContext::new(vec![RenderVoice::CoupledModal(Box::new(voice))], force_config.max_block);
        super::poll(Some(gate)).map_err(RenderError::Coupled)?;
        Self::new(context, controls, force_config.max_controls)
    }
}
fn invalid(what: &'static str) -> RenderError { RenderError::Control { what } }
/// Scheduled normal contact on top of the existing coupled network.
pub mod contact;
