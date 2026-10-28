//! Finite source-plane diagnostics using the canonical pressure gesture compiler.
//! No gain, pitch preset, synthetic microphone or independent source clock exists.

use super::{LoadedApertureNetworkLine, NetworkLineFrame};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::performance::{
    AperturePerformanceConfig, compile_schedule, pressure_candidate,
};
use crate::bernoulli_aperture::tube::TubeDrive;
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use crate::render::schedule::{GestureCompileError, ScheduledControl};
use fs_scenario::gesture::GestureSchedule;

/// Finite source-plane pressure. This internal port signal is useful for
/// diagnostics; it does not describe exterior radiation at a receiver.
pub struct NetworkLinePerformance {
    system: LoadedApertureNetworkLine,
    config: AperturePerformanceConfig,
    schedule: GestureSchedule,
    controls: Vec<ScheduledControl>,
    next: usize,
    held_pressure_pa: f64,
    poisoned: bool,
}

impl NetworkLinePerformance {
    /// Compile the existing pressure schedule at the unchanged physical clock.
    /// A failed physical sample accepts neither its gesture nor either owner;
    /// failed callbacks poison rendering and leave earlier accepted steps visible.
    pub fn new(
        system: LoadedApertureNetworkLine,
        schedule: GestureSchedule,
        config: AperturePerformanceConfig,
    ) -> Result<Self, GestureCompileError> {
        if config.sample_rate_hz == 0
            || config.sample_rate_hz != system.load().sample_rate_hz()
            || config.samples == 0
            || config.samples > (1_u64 << 53)
            || config.samples > system.aperture().max_steps()
            || system.aperture().accepted_steps() != 0
            || !(1..=65_536).contains(&config.max_block)
        {
            return Err(GestureCompileError::Invalid {
                what: "network pressure phrase requires a cold matching clock and finite window/capacity",
            });
        }
        let controls = compile_schedule(&schedule, &config)?;
        Ok(Self {
            system,
            config,
            schedule,
            controls,
            next: 0,
            held_pressure_pa: 0.0,
            poisoned: false,
        })
    }
    /// Complete physical source/network and retained material evidence.
    #[must_use]
    pub fn system(&self) -> &LoadedApertureNetworkLine {
        &self.system
    }
    /// Original canonical gesture schedule.
    #[must_use]
    pub fn schedule(&self) -> &GestureSchedule {
        &self.schedule
    }
    /// Original finite clock and resource declaration.
    #[must_use]
    pub fn config(&self) -> AperturePerformanceConfig {
        self.config
    }
    /// Accepted source-plane observation, including material/contact work.
    #[must_use]
    pub fn last_frame(&self) -> Option<&NetworkLineFrame> {
        self.system.last_frame()
    }
    /// Number of compiled controls committed with their physical samples.
    #[must_use]
    pub fn applied_controls(&self) -> usize {
        self.next
    }
    /// Number of compiled controls still awaiting physical acceptance.
    #[must_use]
    pub fn pending_controls(&self) -> usize {
        self.controls.len() - self.next
    }
    /// Unconsumed physical samples in the explicitly declared phrase.
    #[must_use]
    pub fn remaining_samples(&self) -> u64 {
        self.config.samples - self.system.aperture().accepted_steps()
    }

    fn advance_sample(&mut self) -> Result<f64, RenderError> {
        let sample = self.system.aperture().accepted_steps();
        let (next, pressure) =
            pressure_candidate(&self.controls, self.next, self.held_pressure_pa, sample);
        let frame = self
            .system
            .step(TubeDrive {
                upstream_pressure_pa: pressure,
                body_flow_m3_s: 0.0,
            })
            .map_err(|error| {
                RenderError::Voice(AcousticRealizeError::Nonlinear(error.to_string()))
            })?;
        self.next = next;
        self.held_pressure_pa = pressure;
        Ok(frame.aperture.bore_pressure_pa)
    }
}

impl PressureRenderer for NetworkLinePerformance {
    fn samples_rendered(&self) -> u64 {
        self.system.aperture().accepted_steps()
    }
    fn max_block_len(&self) -> usize {
        self.config.max_block
    }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.validate_sample_count(0)?;
        if rate != self.config.sample_rate_hz {
            return Err(RenderError::Sizing {
                what: "network pressure uses its admitted mechanical clock",
            });
        }
        Ok(())
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned {
            return Err(RenderError::Poisoned);
        }
        if samples > self.remaining_samples() {
            return Err(RenderError::Sizing {
                what: "network pressure exceeds its finite source window",
            });
        }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.validate_sample_count(output.len() as u64)?;
        if output.is_empty() {
            return Err(RenderError::EmptyBlock);
        }
        if output.len() > self.config.max_block {
            return Err(RenderError::Sizing {
                what: "network callback exceeds its declared capacity",
            });
        }
        for slot in output {
            match self.advance_sample() {
                Ok(value) => *slot = value,
                Err(error) => {
                    self.poisoned = true;
                    return Err(error);
                }
            }
        }
        Ok(())
    }
}
