//! Existing scheduler and force compiler hosting simultaneous normal/tangential contacts.
use super::*;
use super::super::super::contact::multiple::{MultiContactConfig, MultiContactModalSystem};
use super::super::super::contact::multiple::friction::ModalFriction;

/// A contact set in one render slot. External controls never replace solved
/// contact forces. The same network state survives all callback boundaries.
pub struct MultiContactModalVoice {
    system: MultiContactModalSystem,
    held_force: Vec<f64>,
    poisoned: bool,
}
impl MultiContactModalVoice {
    /// Host accepted physical states without resetting any component.
    pub fn new(system: MultiContactModalSystem, held_force: Vec<f64>) -> Result<Self, RenderError> {
        if held_force.len() != system.mode_count() || held_force.iter().any(|x| !x.is_finite()) {
            return Err(invalid("multi-contact voice requires finite external forces for every mode"));
        }
        Ok(Self { system, held_force, poisoned: false })
    }
    /// Original components and last complete contact-set diagnostic.
    #[must_use]
    pub const fn system(&self) -> &MultiContactModalSystem { &self.system }
    /// Shared mechanical period [s].
    #[must_use]
    pub fn sample_period_s(&self) -> f64 { self.system.sample_period_s() }
    /// Validate a complete mode assignment before a batch changes any input.
    pub fn validate_force(&self, mode: usize, force: f64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if mode >= self.held_force.len() || !force.is_finite() {
            return Err(invalid("multi-contact force control is outside the finite modal basis"));
        }
        Ok(())
    }
    pub(in crate::render) fn set_force_admitted(&mut self, mode: usize, force: f64) { self.held_force[mode] = force; }
    /// Advance complete joint samples. A failed callback may have an accepted
    /// prefix, so this wrapper refuses continuation; direct system trials remain
    /// retryable. The enclosing scheduler owns callback cancellation boundaries.
    pub fn step_block(&mut self, out: &mut [f64]) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if out.is_empty() { return Err(RenderError::EmptyBlock); }
        let count = u64::try_from(out.len()).map_err(|_| invalid("multi-contact callback length exceeds the clock"))?;
        self.system.samples_rendered().checked_add(count).ok_or_else(|| invalid("multi-contact sample clock overflow"))?;
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
    /// Compile physical actuators with the existing coupled-force compiler, then
    /// attach a jointly solved set of normal contacts to the complete network.
    /// No event timing, modal projection, constitutive law or encoder is copied.
    /// Explicit all-component static preload settles the complete contact set
    /// under the authored contact_set setup/sweep caps and per-contact root caps.
    #[allow(clippy::too_many_arguments)]
    pub fn from_multi_contact_modal_forces(
        voices: Vec<ModalForceVoice>, events: Vec<ModalForceEvent>, force_config: ForceRenderConfig,
        connections: Vec<ModalConnection>, coupling: ModalCouplingConfig,
        contacts: Vec<(ModalContact, ModalContactConfig)>, contact_set: MultiContactConfig, gate: &CancelGate,
    ) -> Result<Self, RenderError> {
        Self::compile_contact_modal_forces(voices, events, force_config, connections, coupling,
            contacts, contact_set, None, gate)
    }

    /// Host explicit 1-D regularized Coulomb friction in the same shared-body
    /// contact solve, actuator compiler, sample clock and callback transaction.
    /// Each normal contact requires Some(authored law) or explicit None.
    /// No static sticking, rigid impact or friction-loaded preload is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn from_frictional_modal_forces(
        voices: Vec<ModalForceVoice>, events: Vec<ModalForceEvent>, force_config: ForceRenderConfig,
        connections: Vec<ModalConnection>, coupling: ModalCouplingConfig,
        contacts: Vec<(ModalContact, ModalContactConfig)>, contact_set: MultiContactConfig,
        friction: Vec<Option<ModalFriction>>, gate: &CancelGate,
    ) -> Result<Self, RenderError> {
        Self::compile_contact_modal_forces(voices, events, force_config, connections, coupling,
            contacts, contact_set, Some(friction), gate)
    }

    #[allow(clippy::too_many_arguments)]
    fn compile_contact_modal_forces(
        mut voices: Vec<ModalForceVoice>, events: Vec<ModalForceEvent>, force_config: ForceRenderConfig,
        connections: Vec<ModalConnection>, coupling: ModalCouplingConfig,
        contacts: Vec<(ModalContact, ModalContactConfig)>, contact_set: MultiContactConfig,
        friction: Option<Vec<Option<ModalFriction>>>, gate: &CancelGate,
    ) -> Result<Self, RenderError> {
        super::super::super::poll(Some(gate)).map_err(RenderError::Coupled)?;
        if friction.is_some() && voices.iter().any(|v| v.initialization != ForceInitialization::RetainState) {
            return Err(invalid("frictional performances require retained states; tangential preload is not inferred"));
        }
        let preload = prepare_preload(&mut voices)?;
        let prepared = Self::from_coupled_modal_forces(voices, events, force_config, connections, coupling, gate)?;
        let mut slots = prepared.context.voices.into_iter();
        let Some(RenderVoice::CoupledModal(voice)) = slots.next() else {
            return Err(invalid("coupled compiler returned an incompatible multi-contact host"));
        };
        if slots.next().is_some() { return Err(invalid("multi-contact host requires one complete network")); }
        let mut voice = *voice;
        if preload {
            voice.system.initialize_contact_equilibrium(&voice.held_force, &contacts, contact_set, gate)
                .map_err(RenderError::Coupled)?;
        }
        let mut system = MultiContactModalSystem::new(voice.system, contacts, contact_set, gate).map_err(RenderError::Coupled)?;
        if let Some(friction) = friction { system = system.with_friction(friction, gate).map_err(RenderError::Coupled)?; }
        let hosted = MultiContactModalVoice::new(system, voice.held_force)?;
        let context = RenderContext::new(vec![RenderVoice::MultiContactModal(Box::new(hosted))], force_config.max_block);
        super::super::super::poll(Some(gate)).map_err(RenderError::Coupled)?;
        Self::new(context, prepared.events, force_config.max_controls)
    }
}
