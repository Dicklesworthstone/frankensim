//! Physical force performances over admitted, mass-normalized modal images.
//!
//! A port column B_p maps force in N to modal force in N/sqrt(kg).
//! The SAME column maps modal velocity back to port velocity: v_p = B_p^T v.
//! Thus F^T v_port = (B F)^T v_modal. Columns must come from the caller's
//! structural basis (for example FEM mode shapes at an actuator); this module
//! does not invent shapes, pitches, radiation transfers or contact laws.
//!
//! Independent ports are summed, not overwritten. Same-sample assignments are
//! applied together, last assignment to a port wins, and projection always uses
//! mode order then port order. Compilation is bounded and entirely off-callback.

use crate::modal_acoustic_time::ModalAcousticTimeModel;
use crate::render::{ControlDelta, ModalStringVoice, RenderContext, RenderError, RenderVoice};
use super::{ScheduledControl, ScheduledRenderer};

/// How a performance's initial held forces relate to its initial modal state.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ForceInitialization {
    /// Preserve the supplied state, including any existing vibration. A new
    /// force applied to a resting model starts a real load-on transient.
    RetainState,
    /// Explicit static preload: initialize q_k = (B F)_k / omega_k^2 and v=0.
    /// Releasing that force later produces a pluck without resetting vibration.
    /// This assumes the load settled BEFORE the window; it is not a strike law.
    StaticPreload,
}

/// One admitted structural image and physical force ports in its exact basis.
pub struct ModalForceVoice {
    model: ModalAcousticTimeModel,
    port_shapes: Vec<Vec<f64>>,
    initial_force_n: Vec<f64>,
    initialization: ForceInitialization,
}

impl ModalForceVoice {
    /// Bind one column of mass-normalized shapes [1/sqrt(kg)] per physical port.
    /// Each column has one value per retained mode, in the model's mode order.
    /// Signed weights and forces are legal; all values must be finite. Zero
    /// weights are legal (a retained basis need not couple to every actuator).
    /// No physical validity is inferred merely from cardinality checks.
    pub fn new(
        model: ModalAcousticTimeModel,
        port_shapes: Vec<Vec<f64>>,
        initial_force_n: Vec<f64>,
        initialization: ForceInitialization,
    ) -> Result<Self, RenderError> {
        if port_shapes.is_empty() || port_shapes.len() != initial_force_n.len() {
            return Err(control("one initial physical force is required per nonempty port set"));
        }
        if port_shapes.iter().any(|column| {
            column.len() != model.modes().len() || column.iter().any(|x| !x.is_finite())
        }) || initial_force_n.iter().any(|x| !x.is_finite()) {
            return Err(control("force-port columns must match the modal basis and all inputs must be finite"));
        }
        Ok(Self { model, port_shapes, initial_force_n, initialization })
    }
}

/// A signed, zero-order-held actuator force assignment on the audio sample clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModalForceEvent {
    /// Apply before this sample. Sample zero is the first output sample.
    pub sample: u64,
    /// Index of the force-driven image in the performance.
    pub voice: usize,
    /// Column index in that image's admitted force-port matrix.
    pub port: usize,
    /// Physical force [N]. Zero releases THIS port, not the other ports.
    pub force_n: f64,
}

/// Explicit compilation and callback budgets; none are enlarged silently.
#[derive(Clone, Copy, Debug)]
pub struct ForceRenderConfig {
    /// Audio rate. Every supplied model must have this exact sample period.
    pub sample_rate_hz: u32,
    /// Largest host callback, in samples.
    pub max_block: usize,
    /// Maximum number of authored physical force events.
    pub max_events: usize,
    /// Maximum emitted per-mode assignments after projection.
    pub max_controls: usize,
    /// Total multiply-add terms admitted across initial projections and all
    /// changed-voice, same-sample groups. Does not bound the later audio render.
    pub max_projection_terms: usize,
}

struct PreparedVoice {
    voice: RenderVoice,
    port_shapes: Vec<Vec<f64>>,
    force_n: Vec<f64>,
    generalized: Vec<f64>,
}

impl ScheduledRenderer {
    /// Check the output clock against every hosted voice before connecting an
    /// audio sink. This never retunes a model or resamples its pressure history.
    pub fn validate_sample_rate(&self, sample_rate_hz: u32) -> Result<(), RenderError> {
        self.context.validate_controls(&[])?;
        if sample_rate_hz == 0 {
            return Err(control("output sample rate must be positive"));
        }
        let expected = f64::from(sample_rate_hz).recip().to_bits();
        for voice in &self.context.voices {
            let period = match voice {
                RenderVoice::ReedBore(reed) => reed.dt,
                RenderVoice::ModalString(string) => string.model.sample_period_s(),
            };
            if period.to_bits() != expected {
                return Err(control("audio sink sample rate must match every hosted voice"));
            }
        }
        Ok(())
    }

    /// Compile independent physical actuator histories to the existing renderer.
    ///
    /// Validation includes ALL authored events before projection or preload.
    /// Changed generalized inputs lower to ordinary `SetModalForce` controls;
    /// the exact-ZOH runtime, work/energy checks, sample scheduling and callback
    /// cancellation are unchanged. There is no per-sample force projection.
    ///
    /// Same-time changes are atomic at the physical-input level: all port values
    /// change before their sum is projected. A released port cannot silence an
    /// independently held one. No action resets a resonating state.
    pub fn from_modal_forces(
        voices: Vec<ModalForceVoice>,
        mut events: Vec<ModalForceEvent>,
        config: ForceRenderConfig,
    ) -> Result<Self, RenderError> {
        if voices.is_empty() || config.sample_rate_hz == 0 || config.max_block == 0 {
            return Err(control("physical performance needs voices, a positive sample rate and block capacity"));
        }
        if events.len() > config.max_events {
            return Err(sizing("physical force event budget exceeded"));
        }
        let period = f64::from(config.sample_rate_hz).recip();
        for voice in &voices {
            if voice.model.sample_period_s().to_bits() != period.to_bits() {
                return Err(control("all modal images must match the declared audio sample rate"));
            }
        }
        for event in &events {
            let voice = voices.get(event.voice)
                .ok_or(RenderError::UnknownVoice { index: event.voice })?;
            if event.sample == u64::MAX || event.port >= voice.port_shapes.len()
                || !event.force_n.is_finite() {
                return Err(control("force event needs a renderable sample, an existing port and finite newtons"));
            }
        }
        events.sort_by_key(|event| event.sample);
        let mut remaining = config.max_projection_terms;
        let mut prepared = Vec::new();
        prepared.try_reserve(voices.len()).map_err(|_| sizing("cannot reserve force voices"))?;
        for mut voice in voices {
            charge(&mut remaining, &voice.port_shapes)?;
            let generalized = project(&voice.port_shapes, &voice.initial_force_n)?;
            if voice.initialization == ForceInitialization::StaticPreload {
                voice.model.initialize_static_equilibrium(&generalized).map_err(RenderError::Modal)?;
            }
            let runtime = ModalStringVoice::new(voice.model, generalized.clone())
                .map_err(RenderError::Modal)?;
            prepared.push(PreparedVoice {
                voice: RenderVoice::ModalString(runtime),
                port_shapes: voice.port_shapes,
                force_n: voice.initial_force_n,
                generalized,
            });
        }
        let mut controls = Vec::new();
        let mut changed = Vec::new();
        let mut at = 0;
        while at < events.len() {
            let sample = events[at].sample;
            changed.clear();
            while at < events.len() && events[at].sample == sample {
                let event = events[at];
                prepared[event.voice].force_n[event.port] = event.force_n;
                changed.push(event.voice);
                at += 1;
            }
            changed.sort_unstable();
            changed.dedup();
            for &index in &changed {
                let voice = &mut prepared[index];
                charge(&mut remaining, &voice.port_shapes)?;
                let projected = project(&voice.port_shapes, &voice.force_n)?;
                let count = projected.iter().zip(&voice.generalized)
                    .filter(|(a, b)| a.to_bits() != b.to_bits()).count();
                let total = controls.len().checked_add(count)
                    .ok_or_else(|| sizing("projected control count overflow"))?;
                if total > config.max_controls {
                    return Err(sizing("projected modal control budget exceeded"));
                }
                controls.try_reserve(count).map_err(|_| sizing("cannot reserve projected controls"))?;
                for (mode, (&next, &previous)) in projected.iter().zip(&voice.generalized).enumerate() {
                    if next.to_bits() != previous.to_bits() {
                        controls.push(ScheduledControl { sample, delta: ControlDelta::SetModalForce {
                            voice: index, mode, force_n_per_sqrt_kg: next,
                        }});
                    }
                }
                voice.generalized = projected;
            }
        }
        let context = RenderContext::new(prepared.into_iter().map(|v| v.voice).collect(), config.max_block);
        Self::new(context, controls, config.max_controls)
    }
}

fn project(columns: &[Vec<f64>], forces: &[f64]) -> Result<Vec<f64>, RenderError> {
    let mut generalized = vec![0.0; columns[0].len()];
    for (mode, value) in generalized.iter_mut().enumerate() {
        for (column, &force) in columns.iter().zip(forces) {
            *value += column[mode] * force;
            if !value.is_finite() {
                return Err(control("physical force projection produced a non-finite modal input"));
            }
        }
    }
    Ok(generalized)
}

fn charge(remaining: &mut usize, columns: &[Vec<f64>]) -> Result<(), RenderError> {
    let terms = columns.len().checked_mul(columns[0].len())
        .ok_or_else(|| sizing("force projection work overflow"))?;
    *remaining = remaining.checked_sub(terms)
        .ok_or_else(|| sizing("force projection work budget exceeded"))?;
    Ok(())
}

fn control(what: &'static str) -> RenderError { RenderError::Control { what } }
fn sizing(what: &'static str) -> RenderError { RenderError::Sizing { what } }

/// Bounded loading of authored modal images and physical force performances.
pub mod file;
