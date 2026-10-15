//! Scheduled exterior pressure from the solved flow of a physical flue outlet.
//!
//! Jet, resonator, gesture interpolation and Rayleigh propagation stay with
//! their existing owners. Only clock/event splits and physical port wiring live
//! here. This observes ONE circular opening in an infinite rigid baffle, with
//! homogeneous exterior gas equal to its adjacent tube region. It is NOT the
//! complete free-field radiation of a flute, its mouth, all holes or a room.
//! The original tonal/estimated jet remains tonal/estimated, not turbulent CFD.

use super::{FlueError, FlueResonator, FlueVoice};
use crate::acoustic_realize::AcousticRealizeError;
use crate::flue_loop::FlueControl;
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::{ControlDelta, RenderError};
use crate::render::schedule::{
    GestureCompileError, PressureGestureBinding, ScheduledControl, compile_pressure_gestures,
};
use fs_scenario::gesture::GestureSchedule;
use fs_vfit::waveguide::network::NetworkNode;

/// Finite physical output and explicit observation/compilation limits.
#[derive(Clone, Copy, Debug)]
pub struct FlueOutletConfig {
    /// Must equal the source's mechanical clock [Hz]. No resampling or retuning.
    pub sample_rate_hz: u32,
    /// Complete physical time window, including any desired propagation tail.
    pub samples: u64,
    /// Largest host callback and allocated launched-wave scratch [samples].
    pub max_block: usize,
    /// Upper bound on control ticks, even if equal inputs are coalesced.
    pub max_controls: usize,
    /// Existing gesture compiler's source-track/event work allowance.
    pub max_compile_work: u64,
    /// Original physical node address of a degree-one passive terminal.
    pub outlet_node: usize,
    /// Position relative to that opening and bounded surface/band resolution.
    pub receiver: CircularOutletReceiver,
    /// Positive finite microphone-pressure ceiling [Pa], never a normalizer.
    pub maximum_abs_pressure_pa: f64,
}

/// One owned jet/resonator/observer/schedule with a common integer sample clock.
/// Internal p+ remains a diagnostic, not the output. Source histories persist
/// across callbacks; source diagnostics allocate per control segment, not per
/// sample. No general allocation-free or real-time qualification is implied.
pub struct FlueOutletPerformance {
    voice: FlueVoice,
    receiver: BaffledPressure,
    config: FlueOutletConfig,
    schedule: GestureSchedule,
    controls: Vec<ScheduledControl>,
    next: usize,
    launched: Vec<f64>,
    poisoned: bool,
}

fn sizing(what: &'static str) -> RenderError { RenderError::Sizing { what } }
fn invalid(what: &'static str) -> GestureCompileError {
    GestureCompileError::Render(RenderError::Control { what })
}
fn source_error(error: FlueError) -> RenderError {
    match error {
        FlueError::Physics(error) => RenderError::Voice(error),
        error => RenderError::Voice(AcousticRealizeError::Nonlinear(error.to_string())),
    }
}

impl FlueOutletPerformance {
    /// Bind one canonical blowing-pressure track to a quiescent physical graph.
    /// Uses the existing pressure compiler's ramp, clock-rounding and event rules.
    /// Every emitted input is checked against the actual jet's declared domain
    /// before any source or observer advances. No source track is silently ignored.
    ///
    /// The receiver uses the ACTUAL outlet radius and LOCAL outlet gas, not the
    /// inlet state or a separately authored amplitude. The existing load still
    /// owns reaction/storage; observation introduces no extra dissipative load.
    /// Above-band nonlinear content is not certified or silently filtered by the
    /// declared Rayleigh band. Output sampling/refinement remains the caller's
    /// responsibility. Zero blowing releases the source without clearing waves.
    ///
    /// # Errors
    /// Invalid source/topology/observer, conflicting clocks, input-domain or
    /// compilation/allocation limits refuse before performance begins.
    pub fn new(voice: FlueVoice, schedule: GestureSchedule, config: FlueOutletConfig)
        -> Result<Self, GestureCompileError>
    {
        if config.sample_rate_hz == 0 || config.samples == 0 || config.samples > (1_u64 << 53)
            || config.max_block == 0 || config.max_block > 1_048_576 || voice.samples_rendered() != 0
            || voice.sample_period_s().to_bits() != (1.0 / f64::from(config.sample_rate_hz)).to_bits()
            || !config.maximum_abs_pressure_pa.is_finite() || config.maximum_abs_pressure_pa <= 0.0 {
            return Err(invalid("flue pressure needs a fresh source, matching clock and finite positive window/output budgets"));
        }
        voice.validate_blowing_pressure(0.0).map_err(|e| GestureCompileError::Render(source_error(e)))?;
        let network = voice.network().ok_or_else(|| invalid("exterior flue pressure requires a physical tube graph, not an internal TMM pressure trace"))?;
        let graph = network.graph();
        let node = config.outlet_node;
        if !matches!(graph.nodes.get(node), Some(NetworkNode::Termination { .. }
            | NetworkNode::Impedance { .. } | NetworkNode::Relaxation { .. })) {
            return Err(invalid("flue outlet observation requires a passive degree-one terminal"));
        }
        let (index, section) = graph.sections.iter().enumerate().find(|(_, s)| s.nodes.contains(&node))
            .ok_or_else(|| invalid("flue outlet has no physical adjacent section"))?;
        let gas = &network.section_gases()[index];
        let receiver = BaffledPressure::circular_outlet(section.radius_m, config.sample_rate_hz,
            config.receiver, RayleighMedium { density: gas.density, sound_speed: gas.sound_speed })
            .map_err(|e| GestureCompileError::Render(RenderError::Voice(AcousticRealizeError::Nonlinear(e))))?;
        let [track] = schedule.tracks() else {
            return Err(invalid("flue performance needs one explicitly typed blowing-pressure track"));
        };
        let rate = schedule.control_rate_hz;
        if rate == 0 || rate > config.sample_rate_hz {
            return Err(invalid("flue controls must have a positive rate no faster than the source clock"));
        }
        let last_tick = u128::from(config.samples - 1) * u128::from(rate) / u128::from(config.sample_rate_hz);
        if last_tick + 1 > config.max_controls as u128
            || track.events.iter().any(|e| e.time_s > last_tick as f64 / f64::from(rate)) {
            return Err(invalid("flue commands exceed the observed window or control-tick budget"));
        }
        let controls = compile_pressure_gestures(&schedule, &[PressureGestureBinding {
            track: track.id.clone(), voice: 0,
        }], config.sample_rate_hz, config.samples, config.max_compile_work)?;
        for control in &controls {
            let ControlDelta::SetBlowingPressure { voice: 0, pressure_pa } = control.delta else {
                return Err(invalid("flue pressure compiler produced a non-pressure or differently addressed control"));
            };
            voice.validate_blowing_pressure(pressure_pa)
                .map_err(|e| GestureCompileError::Render(source_error(e)))?;
        }
        let mut launched = Vec::new();
        launched.try_reserve_exact(config.max_block)
            .map_err(|_| invalid("flue callback scratch exceeds available allocation"))?;
        launched.resize(config.max_block, 0.0);
        Ok(Self { voice, receiver, config, schedule, controls, next: 0, launched, poisoned: false })
    }

    /// Source mechanics, material gas mapping and internal diagnostics, read-only.
    #[must_use]
    pub const fn voice(&self) -> &FlueVoice { &self.voice }
    /// Exact original physical gesture schedule.
    #[must_use]
    pub const fn schedule(&self) -> &GestureSchedule { &self.schedule }
    /// Admitted physical receiver and finite output limits.
    #[must_use]
    pub const fn config(&self) -> &FlueOutletConfig { &self.config }
    /// Original shared receiver, including physical and numerical delay metadata.
    #[must_use]
    pub const fn receiver(&self) -> &BaffledPressure { &self.receiver }
    /// Samples remaining on the source/observer clock.
    #[must_use]
    pub fn remaining_samples(&self) -> u64 { self.config.samples - self.voice.samples_rendered() }

    fn advance_block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        let end = self.voice.samples_rendered() + output.len() as u64;
        let mut written = 0;
        while written < output.len() {
            let now = self.voice.samples_rendered();
            while let Some(control) = self.controls.get(self.next) {
                if control.sample != now { break; }
                let ControlDelta::SetBlowingPressure { pressure_pa, .. } = control.delta else {
                    unreachable!("pressure-only controls were admitted")
                };
                self.voice.apply(FlueControl::SetBlowingPressure(pressure_pa)).map_err(source_error)?;
                self.next += 1;
            }
            let until = self.controls.get(self.next).map_or(end, |c| c.sample.min(end));
            let count = (until - now) as usize;
            let receiver = &mut self.receiver;
            let node = self.config.outlet_node;
            let ceiling = self.config.maximum_abs_pressure_pa;
            let mut destinations = output[written..written + count].iter_mut();
            self.voice.step_block_observed(&mut self.launched[..count], |line| {
                let FlueResonator::Network(network) = line else { unreachable!("owned physical graph") };
                let flow = network.node_frame(node).expect("owned admitted outlet").net_flow_into_node_m3_s;
                let pressure = receiver.step(&[flow])
                    .map_err(|e| FlueError::Physics(AcousticRealizeError::Nonlinear(e)))?;
                if !pressure.is_finite() || pressure.abs() > ceiling {
                    return Err(FlueError::Invalid { what: "flue microphone pressure exceeds its declared finite ceiling" });
                }
                *destinations.next().expect("one observation per source sample") = pressure;
                Ok(())
            }).map_err(source_error)?;
            written += count;
        }
        Ok(())
    }
}

impl PressureRenderer for FlueOutletPerformance {
    fn samples_rendered(&self) -> u64 { self.voice.samples_rendered() }
    fn max_block_len(&self) -> usize { self.config.max_block }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.validate_sample_count(0)?;
        if rate != self.config.sample_rate_hz { return Err(sizing("flue pressure output must use its actual physical sample clock")); }
        Ok(())
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if samples > self.remaining_samples() { return Err(sizing("flue output exceeds its finite source/observer window")); }
        Ok(())
    }
    /// Render complete observer samples. Cancellation belongs at shared host
    /// callback boundaries. A source/observer refusal poisons the performance;
    /// discard the partial callback. Size/rate refusals never mutate anything.
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.validate_sample_count(output.len() as u64)?;
        if output.is_empty() { return Err(RenderError::EmptyBlock); }
        if output.len() > self.config.max_block { return Err(sizing("flue callback exceeds its construction-time capacity")); }
        let result = self.advance_block(output);
        if result.is_err() { self.poisoned = true; }
        result
    }
}

#[cfg(test)]
mod tests;
