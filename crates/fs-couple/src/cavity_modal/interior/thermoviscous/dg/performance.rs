//! Continuing physical-pressure playback of the complete resolved field and bore.
//!
//! Every control is a held physical traction or generalized structural force on
//! the original clock. The field's molecular temperatures, finite conductive
//! walls and moving structure continue through the same reciprocal DG step.
//! The shared pressure, decimation, ensemble and PCM owners consume this stream;
//! observation never substitutes an independent sound source or wall-loss law.

use super::{ThermoviscousDg, ThermoviscousDgFrame};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::performance::ApertureObservation;
use crate::pcm_wav::baffled::BaffledPressure;
use crate::pcm_wav::baffled::gas::{GasPropagation, GasPropagationReport};
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use fs_exec::CancelGate;
use std::sync::Arc;

/// One original physical interval's complete external drives.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermoviscousDgControl {
    /// Held pressures [Pa], in original FE port-mode order. The mode attached to
    /// the bore MUST be zero here: its pressure comes from the reciprocal solve.
    pub traction_pa: Vec<f64>,
    /// Held forces [N/sqrt(kg)], in original structural-coordinate order.
    pub wall_forces_n_per_sqrt_kg: Vec<f64>,
}

/// Explicit finite playback storage and callback bounds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ThermoviscousDgPerformanceConfig {
    /// Must match the original physical interval exactly: dt = 1 / rate.
    /// Use the shared pressure decimator for a lower PCM rate.
    pub sample_rate_hz: u32,
    /// Largest callback, in 1..=65_536 samples.
    pub maximum_block: usize,
    /// Total stored scalar traction and force values across the whole history.
    /// This bounds authored input storage, separately from FE preparation limits.
    pub maximum_control_values: usize,
}

/// A finite pressure stream retaining both original physical owners.
///
/// There is no no-allocation claim: the complete field and DG steppers retain
/// their original candidate allocations. Callback boundaries do no arithmetic.
/// A failed physical step leaves both owners at their preceding accepted state;
/// an observation failure retains the accepted physical step. Either failure
/// poisons the stream, preserving the completed output prefix for inspection.
pub struct ThermoviscousDgPerformance<'a, 'c> {
    system: ThermoviscousDg<'a, 'c>,
    controls: Vec<ThermoviscousDgControl>,
    config: ThermoviscousDgPerformanceConfig,
    observation: ApertureObservation,
    receiver: Option<BaffledPressure>,
    observer_gas: Option<GasPropagationReport>,
    terminal_flows: Vec<f64>,
    gate: Arc<CancelGate>,
    rendered: u64,
    poisoned: bool,
}

fn sizing(what: &'static str) -> RenderError {
    RenderError::Sizing { what }
}

impl<'a, 'c> ThermoviscousDgPerformance<'a, 'c> {
    /// Bind complete physical controls to an unadvanced coupled system.
    /// Nonzero admitted initial field or wall motion is retained. The input
    /// window is finite and shares the original coupling step budget.
    ///
    /// Observation selects either the solved mean interface pressure or the
    /// existing DG computed-baffle microphone, using every actual terminal's
    /// outward flow. Other aperture/network observation variants refuse.
    pub fn new(
        system: ThermoviscousDg<'a, 'c>,
        controls: Vec<ThermoviscousDgControl>,
        observation: ApertureObservation,
        config: ThermoviscousDgPerformanceConfig,
        gate: Arc<CancelGate>,
    ) -> Result<Self, RenderError> {
        if gate.is_requested() {
            return Err(RenderError::Voice(AcousticRealizeError::Cancelled));
        }
        if system.accepted_steps() != 0
            || config.sample_rate_hz == 0
            || system.prepared().duration_s().to_bits()
                != (1.0 / f64::from(config.sample_rate_hz)).to_bits()
            || !(1..=65_536).contains(&config.maximum_block)
            || controls.is_empty()
            || controls.len() as u64 > system.limits().maximum_steps
        {
            return Err(sizing(
                "resolved field playback requires a cold exact clock and bounded finite history",
            ));
        }
        let pressures = system.source().inner().port_mode_count();
        let forces = system.state().wall_velocity_m_sqrt_kg_per_s.len();
        let values = pressures
            .checked_add(forces)
            .and_then(|per_step| per_step.checked_mul(controls.len()))
            .ok_or_else(|| sizing("resolved field control storage overflow"))?;
        if values > config.maximum_control_values {
            return Err(sizing("resolved field controls exceed their declared storage bound"));
        }
        for control in &controls {
            if gate.is_requested() {
                return Err(RenderError::Voice(AcousticRealizeError::Cancelled));
            }
            if control.traction_pa.len() != pressures
                || control.wall_forces_n_per_sqrt_kg.len() != forces
                || control.traction_pa[system.port()] != 0.0
                || control
                    .traction_pa
                    .iter()
                    .chain(&control.wall_forces_n_per_sqrt_kg)
                    .any(|value| !value.is_finite())
            {
                return Err(sizing(
                    "resolved field controls need every finite physical drive and zero authored bore reaction",
                ));
            }
        }
        let receiver = match observation {
            ApertureObservation::Inlet => None,
            ApertureObservation::TubeBaffled(location) => Some(
                system
                    .load()
                    .baffled_receiver(config.sample_rate_hz, location)
                    .map_err(RenderError::Voice)?,
            ),
            _ => {
                return Err(sizing(
                    "resolved field observation must name its DG inlet or computed baffled outlet",
                ));
            }
        };
        let terminal_flows = if receiver.is_some() {
            vec![0.0; system.load().terminal_ports().len()]
        } else {
            Vec::new()
        };
        if gate.is_requested() {
            return Err(RenderError::Voice(AcousticRealizeError::Cancelled));
        }
        Ok(Self {
            system,
            controls,
            config,
            observation,
            receiver,
            observer_gas: None,
            terminal_flows,
            gate,
            rendered: 0,
            poisoned: false,
        })
    }

    /// Attach original exterior-gas propagation to the computed microphone.
    /// This changes the observer only; the same original DG exterior still owns
    /// its reaction. The shared receiver retains material attenuation, phase
    /// and causal delay on each original path before summing pressure.
    pub fn with_observer_gas(mut self, propagation: &GasPropagation) -> Result<Self, RenderError> {
        if self.poisoned || self.system.accepted_steps() != 0 || self.observer_gas.is_some() {
            return Err(sizing("resolved field gas observation requires one cold attachment"));
        }
        let ApertureObservation::TubeBaffled(location) = self.observation else {
            return Err(sizing("resolved field gas propagation requires an exterior microphone"));
        };
        let (receiver, report) = self
            .system
            .load()
            .baffled_receiver_with_gas(
                self.config.sample_rate_hz,
                location,
                propagation,
                &self.gate,
            )
            .map_err(RenderError::Voice)?;
        if self.gate.is_requested() {
            return Err(RenderError::Voice(AcousticRealizeError::Cancelled));
        }
        self.receiver = Some(receiver);
        self.observer_gas = Some(report);
        Ok(self)
    }

    pub const fn system(&self) -> &ThermoviscousDg<'a, 'c> {
        &self.system
    }

    pub fn controls(&self) -> &[ThermoviscousDgControl] {
        &self.controls
    }

    pub const fn config(&self) -> ThermoviscousDgPerformanceConfig {
        self.config
    }

    pub const fn observation(&self) -> ApertureObservation {
        self.observation
    }

    pub fn last_frame(&self) -> Option<&ThermoviscousDgFrame> {
        self.system.last_frame()
    }

    pub fn baffled_receiver(&self) -> Option<&BaffledPressure> {
        self.receiver.as_ref()
    }

    pub fn observer_gas(&self) -> Option<&GasPropagationReport> {
        self.observer_gas.as_ref()
    }

    pub fn cancellation_gate(&self) -> &CancelGate {
        &self.gate
    }

    pub fn remaining_samples(&self) -> u64 {
        self.controls.len() as u64 - self.rendered
    }

    fn advance(&mut self) -> Result<f64, RenderError> {
        let control = self
            .controls
            .get(self.system.accepted_steps() as usize)
            .ok_or_else(|| sizing("resolved field physical history is exhausted"))?;
        let pressure = self
            .system
            .step(
                &control.traction_pa,
                &control.wall_forces_n_per_sqrt_kg,
                &self.gate,
            )
            .map_err(|error| RenderError::Voice(AcousticRealizeError::Nonlinear(error.to_string())))?
            .mean_pressure_pa;
        let pressure = if let Some(receiver) = &mut self.receiver {
            let terminals = self.system.load().last_terminal_frames();
            if terminals.is_empty() || terminals.len() != self.terminal_flows.len() {
                return Err(sizing("resolved field microphone requires every original terminal flow"));
            }
            for (flow, terminal) in self.terminal_flows.iter_mut().zip(terminals) {
                *flow = terminal.outward_flow_m3_s;
            }
            receiver
                .step_trace(&self.terminal_flows)
                .map_err(|error| RenderError::Voice(AcousticRealizeError::Nonlinear(error)))?
        } else {
            pressure
        };
        if !pressure.is_finite() {
            return Err(sizing("resolved field observer produced nonfinite physical pressure"));
        }
        self.rendered += 1;
        Ok(pressure)
    }
}

impl PressureRenderer for ThermoviscousDgPerformance<'_, '_> {
    fn samples_rendered(&self) -> u64 {
        self.rendered
    }

    fn max_block_len(&self) -> usize {
        self.config.maximum_block
    }

    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.validate_sample_count(0)?;
        if rate != self.config.sample_rate_hz {
            return Err(sizing("resolved field pressure uses its original physical clock"));
        }
        Ok(())
    }

    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned {
            return Err(RenderError::Poisoned);
        }
        if samples > self.remaining_samples() {
            return Err(sizing("resolved field pressure exceeds its finite physical history"));
        }
        Ok(())
    }

    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.validate_sample_count(output.len() as u64)?;
        if output.is_empty() {
            return Err(RenderError::EmptyBlock);
        }
        if output.len() > self.config.maximum_block {
            return Err(sizing("resolved field callback exceeds its declared capacity"));
        }
        for slot in output {
            match self.advance() {
                Ok(pressure) => *slot = pressure,
                Err(error) => {
                    self.poisoned = true;
                    return Err(error);
                }
            }
        }
        Ok(())
    }
}
