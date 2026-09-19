//! Existing scheduler hosting an implicitly contacting modal network.
use super::super::contact::{ContactModalSystem, ModalContact, ModalContactConfig};
use super::super::{ModalConnection, ModalCouplingConfig};
use super::super::super::{ForceInitialization, ForceRenderConfig, ModalForceEvent, ModalForceVoice};
use crate::render::{RenderContext, RenderError, RenderVoice};
use crate::render::schedule::ScheduledRenderer;
use fs_exec::CancelGate;

/// Contact-capable network in one existing render slot. External mode controls
/// never overwrite the contact reaction, which is solved again each sample.
/// The contact owner allocates; this is not an allocation-free callback claim.
pub struct ContactModalVoice {
    system: ContactModalSystem,
    held_force: Vec<f64>,
    poisoned: bool,
}
impl ContactModalVoice {
    /// Wrap an admitted contact system without resetting vibration or contact.
    pub fn new(system: ContactModalSystem, held_force: Vec<f64>) -> Result<Self,RenderError> {
        if held_force.len()!=system.mode_count() || held_force.iter().any(|x|!x.is_finite()) {
            return Err(invalid("contact voice requires finite external forces for all modes"));
        }
        Ok(Self {system,held_force,poisoned:false})
    }
    /// Accepted physical states and contact-inclusive diagnostics.
    #[must_use]
    pub const fn system(&self) -> &ContactModalSystem { &self.system }
    /// Shared mechanical period [s].
    #[must_use]
    pub fn sample_period_s(&self) -> f64 { self.system.sample_period_s() }
    /// Validate a held external-mode control before the entire batch applies.
    pub fn validate_force(&self, mode:usize, force:f64) -> Result<(),RenderError> {
        if self.poisoned {return Err(RenderError::Poisoned);}
        if mode>=self.held_force.len() || !force.is_finite() {return Err(invalid("contact force control is outside the finite modal basis"));}
        Ok(())
    }
    pub(in crate::render) fn set_force_admitted(&mut self, mode:usize, force:f64) {self.held_force[mode]=force;}
    /// Advance complete contact transactions. An error poisons this wrapper
    /// because a prefix of the host callback may have completed; direct system
    /// trials remain retryable. The enclosing scheduler owns cancellation.
    pub fn step_block(&mut self,out:&mut [f64]) -> Result<(),RenderError> {
        if self.poisoned {return Err(RenderError::Poisoned);}
        if out.is_empty() {return Err(RenderError::EmptyBlock);}
        let count=u64::try_from(out.len()).map_err(|_|invalid("contact callback length exceeds the clock"))?;
        self.system.samples_rendered().checked_add(count).ok_or_else(||invalid("contact sample clock overflow"))?;
        for sample in out {
            match self.system.step(&self.held_force) {
                Ok(frame)=>*sample=frame.observer_pressure_pa,
                Err(error)=>{self.poisoned=true;return Err(RenderError::Coupled(error));}
            }
        }
        Ok(())
    }
}

impl ScheduledRenderer {
    /// Reuse the existing physical force compiler and coupled network scheduler,
    /// adding one implicit two-body contact. No force-history synthesis or event
    /// timing code is duplicated. Either retain supplied states, or explicitly settle
    /// ALL zero-state components against the contact before sample zero. The
    /// preload uses the existing root cap, one joint sweep (one contact), and
    /// coupling.max_setup_terms as its static setup-work ceiling.
    #[allow(clippy::too_many_arguments)] // all physical owners and budgets are explicit
    pub fn from_contact_modal_forces(
        mut voices:Vec<ModalForceVoice>, events:Vec<ModalForceEvent>, force_config:ForceRenderConfig,
        connections:Vec<ModalConnection>, coupling:ModalCouplingConfig,
        contact:ModalContact, contact_config:ModalContactConfig, gate:&CancelGate,
    ) -> Result<Self,RenderError> {
        super::super::poll(Some(gate)).map_err(RenderError::Coupled)?;
        let preload = prepare_preload(&mut voices)?;
        let prepared=Self::from_coupled_modal_forces(voices,events,force_config,connections,coupling,gate)?;
        let mut slots=prepared.context.voices.into_iter();
        let Some(RenderVoice::CoupledModal(voice))=slots.next() else {
            return Err(invalid("coupled compiler returned an incompatible contact host"));
        };
        if slots.next().is_some() {return Err(invalid("contact host requires exactly one complete network"));}
        let mut voice=*voice;
        let pair = (contact, contact_config);
        if preload {
            let config = super::super::contact::multiple::MultiContactConfig {
                max_contacts: 1, max_sweeps: 1, max_setup_terms: coupling.max_setup_terms,
            };
            voice.system.initialize_contact_equilibrium(&voice.held_force, std::slice::from_ref(&pair), config, gate)
                .map_err(RenderError::Coupled)?;
        }
        let system=ContactModalSystem::new(voice.system,pair.0,pair.1,gate).map_err(RenderError::Coupled)?;
        let hosted=ContactModalVoice::new(system,voice.held_force)?;
        let context=RenderContext::new(vec![RenderVoice::ContactModal(Box::new(hosted))],force_config.max_block);
        super::super::poll(Some(gate)).map_err(RenderError::Coupled)?;
        Self::new(context,prepared.events,force_config.max_controls)
    }
}
fn invalid(what:&'static str)->RenderError {RenderError::Control {what}}

// Defer settling until every contact is present. Do not let the bilateral-only
// compiler initialize a physically different equilibrium first.
fn prepare_preload(voices: &mut [ModalForceVoice]) -> Result<bool, RenderError> {
    let preload = voices.iter().any(|v| v.initialization == ForceInitialization::StaticPreload);
    if preload {
        if voices.iter().any(|v| v.initialization != ForceInitialization::StaticPreload
            || v.model.states().iter().any(|s| s.displacement_m_sqrt_kg != 0.0 || s.velocity_m_sqrt_kg_per_s != 0.0)) {
            return Err(invalid("contact preload requires every component to request static-preload with zero input Q/V"));
        }
        for voice in voices { voice.initialization = ForceInitialization::RetainState; }
    }
    Ok(preload)
}

/// Scheduled rendering of a simultaneously solved normal-contact set.
pub mod multiple;
