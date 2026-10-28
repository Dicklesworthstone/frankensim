//! Finite material-valve phrases on the original exact DG composition.
//!
//! Scheduling selects pressure and source-mesh loads before each mechanical
//! solve. It neither reduces a paired/modal valve to a scalar reed nor rebuilds
//! its material memory. Exterior observation ingests each accepted terminal
//! flow once, using the same aperture and gas as the radiation reaction.

use super::LoadedApertureDg;
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::dg::DgApertureFrame;
use crate::bernoulli_aperture::performance::{
    ApertureObservation, AperturePerformanceConfig, compile_schedule, pressure_candidate,
};
use crate::bernoulli_aperture::performance::force::ForceSchedule;
use crate::bernoulli_aperture::tube::TubeDrive;
use crate::bernoulli_aperture::two_sided::observer::schedule::spatial_pressures::PressureSchedule;
use crate::pcm_wav::baffled::BaffledPressure;
use crate::pcm_wav::baffled::gas::{GasPropagation, GasPropagationReport};
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use crate::render::schedule::{GestureCompileError, ScheduledControl};
use fs_scenario::gesture::GestureSchedule;

mod controls;

fn invalid(what: &'static str) -> GestureCompileError {
    GestureCompileError::Invalid { what }
}
fn sizing(what: &'static str) -> RenderError { RenderError::Sizing { what } }

/// One finite physical source clock, with immutable access to all retained owners.
/// Initial motion/material/contact energy is allowed. A failed source step leaves
/// its owners unchanged; any failed callback poisons further rendering. Successful
/// physical samples and their controls remain inspectable even if a later
/// observation fails. There is no note envelope, pressure gain or implicit reset.
pub struct LoadedDgPerformance {
    system: LoadedApertureDg,
    observation: ApertureObservation,
    receiver: Option<BaffledPressure>,
    observer_gas: Option<GasPropagationReport>,
    terminal_flows: Vec<f64>,
    config: AperturePerformanceConfig,
    schedule: GestureSchedule,
    controls: Vec<ScheduledControl>,
    next: usize,
    held_pressure_pa: f64,
    plate_forces: [Option<ForceSchedule>; 2],
    spatial_pressures: [Option<PressureSchedule>; 2],
    last: Option<DgApertureFrame>,
    poisoned: bool,
}

impl LoadedDgPerformance {
    /// Compile the canonical blowing-pressure track against a fresh DG source.
    /// The finite window and output rate are mechanical samples. Explicit rate
    /// conversion belongs to the existing pressure/PCM adapter, not the physics.
    /// Only the actual mean inlet pressure or computed baffled outlet is admitted.
    pub fn new(
        system: LoadedApertureDg,
        observation: ApertureObservation,
        schedule: GestureSchedule,
        config: AperturePerformanceConfig,
    ) -> Result<Self, GestureCompileError> {
        let valve = system.aperture();
        if config.sample_rate_hz == 0 || config.samples == 0
            || config.samples > (1_u64 << 53)
            || config.samples > valve.max_steps()
            || !(1..=65_536).contains(&config.max_block)
            || valve.accepted_steps() != 0 || system.load().state().step_count() != 0
            || valve.time_step_s().to_bits()
                != (1.0 / f64::from(config.sample_rate_hz)).to_bits()
        {
            return Err(invalid("loaded DG performance requires a cold matching clock and finite window/capacity"));
        }
        let receiver = match observation {
            ApertureObservation::Inlet => None,
            ApertureObservation::TubeBaffled(location) => Some(
                system.load().baffled_receiver(config.sample_rate_hz, location)
                    .map_err(|error| GestureCompileError::Render(RenderError::Voice(error)))?,
            ),
            _ => return Err(invalid("loaded DG observation must name its inlet or computed baffled outlet")),
        };
        let controls = compile_schedule(&schedule, &config)?;
        let terminal_flows = if receiver.is_some() {
            vec![0.0; system.load().terminal_ports().len()]
        } else { Vec::new() };
        Ok(Self {
            system, observation, receiver, observer_gas: None, terminal_flows, config, schedule, controls, next: 0,
            held_pressure_pa: 0.0, plate_forces: [None, None],
            spatial_pressures: [None, None], last: None, poisoned: false,
        })
    }

    /// Apply an explicit physical exterior propagation law before the first
    /// sample. Material modes, duct state, reactions and accepted flow are
    /// unchanged; the shared surface receiver preserves every path and its
    /// original exterior material, including gas dispersion and causal delay.
    pub fn with_observer_gas(
        mut self, propagation: &GasPropagation, gate: &fs_exec::CancelGate,
    ) -> Result<Self, GestureCompileError> {
        if self.poisoned || self.system.aperture().accepted_steps() != 0
            || self.observer_gas.is_some()
        {
            return Err(invalid("material gas observation requires one cold attachment"));
        }
        let ApertureObservation::TubeBaffled(receiver) = self.observation else {
            return Err(invalid("material gas propagation requires an exterior receiver"));
        };
        let (receiver, report) = self.system.load().baffled_receiver_with_gas(
            self.config.sample_rate_hz, receiver, propagation, gate,
        ).map_err(|error| GestureCompileError::Render(RenderError::Voice(error)))?;
        self.receiver = Some(receiver);
        self.observer_gas = Some(report);
        Ok(self)
    }

    /// Applied exterior law and the source-derived individual path responses.
    pub fn observer_gas(&self) -> Option<&GasPropagationReport> { self.observer_gas.as_ref() }

    /// Complete original mechanical, material and acoustic owners.
    #[must_use]
    pub const fn system(&self) -> &LoadedApertureDg { &self.system }
    /// Original canonical gestures, unaffected by output block partitioning.
    #[must_use]
    pub const fn schedule(&self) -> &GestureSchedule { &self.schedule }
    /// Declared physical window, clock and compilation budgets.
    #[must_use]
    pub const fn config(&self) -> &AperturePerformanceConfig { &self.config }
    /// Physical location of this pressure trace.
    #[must_use]
    pub const fn observation(&self) -> ApertureObservation { self.observation }
    /// Last accepted complete physical balance, including source pressure work.
    #[must_use]
    pub const fn last_frame(&self) -> Option<&DgApertureFrame> { self.last.as_ref() }
    /// Original retarded receiver, including its physical and numerical delays.
    #[must_use]
    pub fn baffled_receiver(&self) -> Option<&BaffledPressure> { self.receiver.as_ref() }
    /// Controls whose physical sample has been accepted.
    #[must_use]
    pub fn applied_controls(&self) -> &[ScheduledControl] { &self.controls[..self.next] }
    /// Controls still pending on the mechanical clock.
    #[must_use]
    pub fn pending_controls(&self) -> &[ScheduledControl] { &self.controls[self.next..] }
    /// Remaining physical samples, before any downstream rate conversion.
    #[must_use]
    pub fn remaining_samples(&self) -> u64 {
        self.config.samples - self.system.aperture().accepted_steps()
    }

    // Pressure compilation retains the existing independent work allowance.
    // Both blades' auxiliary loads share ONE remaining storage/work allowance.
    fn remaining_control_budgets(&self) -> Result<(usize, u64), GestureCompileError> {
        let assignments = self.controls.len() as u128
            + self.plate_forces.iter().flatten().map(|f| f.events.len() as u128).sum::<u128>()
            + self.spatial_pressures.iter().flatten().map(|p| p.events.len() as u128).sum::<u128>();
        let work = self.plate_forces.iter().flatten().map(|f| u128::from(f.work)).sum::<u128>()
            + self.spatial_pressures.iter().flatten().map(|p| u128::from(p.work)).sum::<u128>();
        if assignments > self.config.max_controls as u128
            || work > u128::from(self.config.max_compile_work)
        {
            return Err(invalid("loaded DG source controls exceed their shared compilation allowance"));
        }
        Ok((self.config.max_controls - assignments as usize,
            self.config.max_compile_work - work as u64))
    }

    fn advance_sample(&mut self) -> Result<f64, RenderError> {
        let sample = self.system.aperture().accepted_steps();
        let (next, pressure) = pressure_candidate(
            &self.controls, self.next, self.held_pressure_pa, sample,
        );
        let forces = std::array::from_fn(|i|
            self.plate_forces[i].as_ref().map_or(0.0, |f| f.candidate(sample)));
        let pressures = std::array::from_fn(|i|
            self.spatial_pressures[i].as_ref().map_or(0.0, |p| p.candidate(sample)));
        let frame = self.system.step_with_forces_and_pressures(
            TubeDrive { upstream_pressure_pa: pressure, body_flow_m3_s: 0.0 },
            forces, pressures,
        ).map_err(RenderError::Voice)?;

        // Publish each accepted physical sample before its fallible observer.
        // A later refusal must not make accepted source histories appear pending.
        self.next = next;
        self.held_pressure_pa = pressure;
        for force in self.plate_forces.iter_mut().flatten() { force.accept(sample); }
        for pressure in self.spatial_pressures.iter_mut().flatten() { pressure.accept(sample); }
        self.last = Some(frame);
        let frame = self.last.as_ref().expect("accepted physical frame");
        let output = match self.observation {
            ApertureObservation::Inlet => frame.aperture.bore_pressure_pa,
            ApertureObservation::TubeBaffled(_) => {
                if frame.terminals.len() != self.terminal_flows.len() || frame.terminals.is_empty() {
                    return Err(sizing("loaded DG microphone requires every accepted terminal flow"));
                }
                for (flow, terminal) in self.terminal_flows.iter_mut().zip(&frame.terminals) {
                    *flow = terminal.outward_flow_m3_s;
                }
                self.receiver.as_mut().expect("admitted original DG receiver")
                    .step_trace(&self.terminal_flows)
                    .map_err(|error| RenderError::Voice(AcousticRealizeError::Nonlinear(error)))?
            }
            _ => unreachable!("physical observation admitted at construction"),
        };
        if !output.is_finite() {
            return Err(sizing("loaded DG observation produced a nonfinite physical pressure"));
        }
        Ok(output)
    }
}

impl PressureRenderer for LoadedDgPerformance {
    fn samples_rendered(&self) -> u64 { self.system.aperture().accepted_steps() }
    fn max_block_len(&self) -> usize { self.config.max_block }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.validate_sample_count(0)?;
        if rate != self.config.sample_rate_hz {
            return Err(sizing("loaded DG pressure uses its admitted mechanical clock"));
        }
        Ok(())
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if samples > self.remaining_samples() {
            return Err(sizing("loaded DG pressure exceeds its finite source window"));
        }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.validate_sample_count(output.len() as u64)?;
        if output.is_empty() { return Err(RenderError::EmptyBlock); }
        if output.len() > self.config.max_block {
            return Err(sizing("loaded DG callback exceeds its declared capacity"));
        }
        for slot in output {
            match self.advance_sample() {
                Ok(value) => *slot = value,
                Err(error) => { self.poisoned = true; return Err(error); }
            }
        }
        Ok(())
    }
}
