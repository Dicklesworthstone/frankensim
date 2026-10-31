//! Scheduled exterior pressure from the solved flow of a physical flue outlet.
//!
//! Jet, resonator, gesture interpolation and Rayleigh propagation stay with
//! their existing owners. Only clock/event splits and physical port wiring live
//! here. The legacy constructor observes one opening in its adjacent gas;
//! `with_radiation` binds several coplanar openings and an explicit exterior gas.
//! Observation is an infinite-baffle image. A separately installed MIMO source
//! load supplies mutual reaction and retains its medium/port binding here.
//! An explicit rectangular mouth uses its actual signed inlet flow and the SAME
//! exterior gas, clock and Rayleigh kernel. No complete free-field flute/body
//! radiation, turbulent jet noise or room scattering follows.
//! The original tonal/estimated jet remains tonal/estimated, not turbulent CFD.

/// Versioned data-driven flue construction through this same pressure renderer.
pub mod file;
pub(in crate::flue_loop) mod controls;
pub mod aperture_flows;

use super::{FlueError, FlueResonator, FlueVoice};
use crate::acoustic_realize::AcousticRealizeError;
#[cfg(test)]
use crate::flue_loop::FlueControl;
#[cfg(test)]
use crate::pcm_wav::baffled::RayleighMedium;
use crate::pcm_wav::baffled::gas::{GasPropagation, GasPropagationReport};
use crate::pcm_wav::baffled::openings::RectangularOutletReceiver;
use crate::pcm_wav::baffled::openings::duct::{
    DuctOpening, DuctPressure, DuctRadiationSpec, DuctReconstruction, DuctView,
};
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver};
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use crate::render::schedule::GestureCompileError;
use controls::ScheduledFlueControl;
use fs_material::gas::GasState;
use fs_scenario::gesture::GestureSchedule;
#[cfg(test)]
use fs_vfit::waveguide::network::NetworkNode;

/// Finite physical output and explicit observation/compilation limits.
#[derive(Clone, Copy, Debug)]
pub struct FlueOutletConfig {
    /// PCM/observation clock [Hz]. Legacy constructors require the same source
    /// clock; `with_multirate_radiation` explicitly admits a faster source.
    pub sample_rate_hz: u32,
    /// Complete physical time window, including any desired propagation tail.
    pub samples: u64,
    /// Largest host callback [OUTPUT samples]. The launched-wave scratch has
    /// this many times the mechanics/output ratio, at most 1,048,576 samples.
    pub max_block: usize,
    /// Upper bound on all track/tick assignments before equal values coalesce.
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

/// Additional physical opening in the same stationary z=0 rigid baffle.
/// Its radius is taken from the actual adjacent duct, never independently set.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlueRadiatingOutlet {
    /// Original passive terminal node address, distinct from every other outlet.
    pub node: usize,
    /// Center [m] in the baffle plane; the primary outlet is the origin.
    pub center_m: [f64; 2],
}

/// Explicit common exterior environment and additional resolved radiation ports.
#[derive(Clone, Debug, PartialEq)]
pub struct FlueRadiationInput {
    /// Complete homogeneous exterior gas, independently chosen from bore regions.
    /// Static pressure must match; this does not solve mixing or interface flow.
    pub exterior: GasState,
    /// Additional openings; `config.outlet_node` is always included at the origin.
    /// All share the config's receiver position, band and surface resolution.
    pub additional_outlets: Vec<FlueRadiatingOutlet>,
}

/// Optional physical mouth in the common stationary z=0 baffle.
/// An explicit aperture can differ from the jet's width and cut-up. Positive
/// outward flow is the negative of flow into the
/// bore; this is an observation of the solved source, not another jet model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlueMouthRadiation {
    /// Rectangle center [m] relative to the primary outlet at the origin.
    pub center_m: [f64; 2],
    /// Independent acoustic rectangle width and length [m]. Omission uses
    /// the source's flue width and cut-up, appropriate to a rectangular window.
    /// This changes observation and selected reaction, never jet geometry.
    /// Both characteristic and DG routes retain this same physical aperture.
    pub aperture_size_m: Option<[f64; 2]>,
    /// Midpoint surface cells along the acoustic rectangle's two dimensions.
    pub cells: [usize; 2],
}

impl FlueMouthRadiation {
    /// Resolve the acoustic aperture independently of the source jet.
    #[must_use]
    pub fn size_m(&self, geometry: &crate::flue_loop::FlueGeometry) -> [f64; 2] {
        self.aperture_size_m.unwrap_or([geometry.flue_width_m, geometry.cut_up_m])
    }
}

struct MouthPressure {
    input: FlueMouthRadiation,
    inlet_node: usize,
    receiver: BaffledPressure,
    observer_gas: Option<GasPropagationReport>,
    trace: Vec<f64>,
    filled: usize,
}
impl MouthPressure {
    fn push(&mut self, network: &super::FlueNetwork) -> Result<Option<f64>, FlueError> {
        // The network defines flow into a node from its adjacent tube. At the
        // source node this is (p_minus-p_plus)/Z, hence already outward -q.
        let flow = network.node_frame(self.inlet_node)
            .ok_or(FlueError::Invalid { what: "accepted mouth port has no flow" })?
            .net_flow_into_node_m3_s;
        if !flow.is_finite() {
            return Err(FlueError::Invalid { what: "mouth volume flow must be finite" });
        }
        self.trace[self.filled] = flow;
        let pressure = if self.filled + 1 == self.trace.len() {
            Some(self.receiver.step_trace(&self.trace)
                .map_err(|e| FlueError::Physics(AcousticRealizeError::Nonlinear(e)))?)
        } else { None };
        self.filled = if pressure.is_some() { 0 } else { self.filled + 1 };
        Ok(pressure)
    }
}

/// One owned jet/resonator/observer/schedule with a common integer sample clock.
/// Internal p+ remains a diagnostic, not the output. Source histories persist
/// across callbacks; source diagnostics allocate per control segment, not per
/// sample. No general allocation-free or real-time qualification is implied.
pub struct FlueOutletPerformance {
    voice: FlueVoice,
    receiver: DuctPressure,
    config: FlueOutletConfig,
    radiation: Option<FlueRadiationInput>,
    outlet_nodes: Vec<usize>,
    mechanics_rate_hz: u32,
    mechanics_per_output: usize,
    schedule: GestureSchedule,
    controls: Vec<ScheduledFlueControl>,
    next: usize,
    launched: Vec<f64>,
    mouth: Option<MouthPressure>,
    flow_tap: Option<aperture_flows::ApertureFlowTap>,
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
    /// Bind a canonical blowing-pressure track and optional jet-angle track to
    /// a quiescent physical graph. Both use the existing gesture ramp semantics
    /// and apply on the same mechanical clock before their named sample.
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
        Self::build(voice, schedule, config, None, None)
    }

    /// Observe the signed solved flows of multiple openings in ONE exterior gas.
    /// All ports feed one shared Rayleigh kernel at each accepted source sample;
    /// interference and propagation delays are retained before PCM conversion.
    /// The original scheduler, source, wall histories and pressure ceiling apply.
    /// `receiver().mode_pressures` exposes each opening's delayed contribution.
    ///
    /// Openings must be nonoverlapping circles in the same infinite rigid baffle.
    /// Prescribed terminal loads remain unchanged: this observation neither adds
    /// radiation damping nor solves inter-opening radiation reaction, gas mixing
    /// or the pressure transmission at a material interface. Choosing an exterior
    /// gas does not silently rewrite the bore or its authored termination laws.
    ///
    /// # Errors
    /// Original performance refusals, incoherent/different-pressure exterior gas,
    /// duplicate or nonterminal nodes, overlapping or over-budget openings.
    pub fn with_radiation(voice: FlueVoice, schedule: GestureSchedule,
        config: FlueOutletConfig, radiation: FlueRadiationInput)
        -> Result<Self, GestureCompileError>
    {
        Self::build(voice, schedule, config, Some(radiation), None)
    }

    /// Run an already-constructed high-rate jet/bore through the same broadband
    /// differentiating and anti-alias receiver used by structural sound sources.
    /// `mechanics_rate_hz` must match that source exactly; this does not retime,
    /// rebuild or retune it. The integer source/output ratio must be in 1..=16.
    /// Ratio one selects broadband reconstruction without changing the source.
    ///
    /// Gestures compile on the MECHANICAL clock. A control can therefore land
    /// inside an audio frame without throwing away preceding outlet-flow samples.
    /// Every accepted substep of every selected outlet enters the receiver, not
    /// just the last flow of an audio frame. Callback boundaries preserve both
    /// the physical state and every anti-alias/derivative/propagation history.
    /// Filter delay is explicit in OUTPUT samples via `receiver()`; include its
    /// tail in the requested window. The physical source is never delay-shifted.
    ///
    /// The selected source/material/jet domains and one-way infinite-baffle
    /// radiation limits are unchanged. This is improved numerical observation,
    /// not additional jet physics, modes, empirical data or acoustic reaction.
    ///
    /// # Errors
    /// Original physical refusals, conflicting/noninteger clocks, bandwidth or
    /// scratch/window budgets. Admission does not advance a source or receiver.
    pub fn with_multirate_radiation(voice: FlueVoice, schedule: GestureSchedule,
        config: FlueOutletConfig, radiation: FlueRadiationInput, mechanics_rate_hz: u32)
        -> Result<Self, GestureCompileError>
    {
        Self::build(voice, schedule, config, Some(radiation), Some(mechanics_rate_hz))
    }

    fn build(voice: FlueVoice, schedule: GestureSchedule, config: FlueOutletConfig,
        radiation: Option<FlueRadiationInput>, reconstruction_rate: Option<u32>)
        -> Result<Self, GestureCompileError>
    {
        let mechanics_rate_hz = reconstruction_rate.unwrap_or(config.sample_rate_hz);
        if config.sample_rate_hz == 0 || mechanics_rate_hz == 0
            || mechanics_rate_hz % config.sample_rate_hz != 0
            || !(1..=16).contains(&(mechanics_rate_hz / config.sample_rate_hz)) {
            return Err(invalid("flue needs an explicit integer mechanics/output clock ratio in 1..=16"));
        }
        let mechanics_per_output = (mechanics_rate_hz / config.sample_rate_hz) as usize;
        let mechanical_samples = config.samples.checked_mul(mechanics_per_output as u64)
            .filter(|&n| n <= (1_u64 << 53))
            .ok_or_else(|| invalid("flue mechanical time window exceeds its exact sample budget"))?;
        let scratch_samples = config.max_block.checked_mul(mechanics_per_output)
            .filter(|&n| n <= 1_048_576)
            .ok_or_else(|| invalid("flue mechanical callback scratch exceeds its sample budget"))?;
        if config.sample_rate_hz == 0 || config.samples == 0 || config.samples > (1_u64 << 53)
            || config.max_block == 0 || config.max_block > 1_048_576 || voice.samples_rendered() != 0
            || voice.sample_period_s().to_bits() != (1.0 / f64::from(mechanics_rate_hz)).to_bits()
            || !config.maximum_abs_pressure_pa.is_finite() || config.maximum_abs_pressure_pa <= 0.0 {
            return Err(invalid("flue pressure needs a fresh source, matching clock and finite positive window/output budgets"));
        }
        voice.validate_blowing_pressure(0.0).map_err(|e| GestureCompileError::Render(source_error(e)))?;
        let network = voice.network().ok_or_else(|| invalid("exterior flue pressure requires a physical tube graph, not an internal TMM pressure trace"))?;
        if radiation.as_ref().is_some_and(|input|
            input.additional_outlets.len() >= crate::pcm_wav::baffled::MAX_RAYLEIGH_INPUTS) {
            return Err(invalid("flue outlet count exceeds the shared receiver input budget"));
        }
        let mut openings = vec![DuctOpening { network: 0, node: config.outlet_node, center_m: [0.0; 2] }];
        if let Some(input) = &radiation {
            openings.extend(input.additional_outlets.iter().map(|outlet| DuctOpening {
                network: 0, node: outlet.node, center_m: outlet.center_m,
            }));
        }
        let outlet_nodes = openings.iter().map(|p| p.node).collect();
        let reconstruction = match reconstruction_rate {
            None => DuctReconstruction::FiniteDifference { sample_rate_hz: config.sample_rate_hz },
            Some(rate) => DuctReconstruction::Differentiating {
                mechanics_rate_hz: rate, output_rate_hz: config.sample_rate_hz,
            },
        };
        let receiver = DuctPressure::new(&[DuctView {
            graph: network.graph(), gases: network.section_gases(),
        }], DuctRadiationSpec { openings, exterior: radiation.as_ref().map(|r| r.exterior),
            receiver: config.receiver, reconstruction,
        }).map_err(|e| GestureCompileError::Render(RenderError::Voice(AcousticRealizeError::Nonlinear(e))))?;
        if let Some(environment) = network.mimo_environment() {
            environment.validate_receiver(network.graph(), receiver.source(), *receiver.exterior_gas())
                .map_err(|e| GestureCompileError::Render(RenderError::Voice(e)))?;
        }
        let controls = controls::compile(&voice, &schedule,
            controls::FlueControlLimits { max_controls: config.max_controls,
                max_compile_work: config.max_compile_work },
            mechanics_rate_hz, mechanical_samples)?;
        let mut launched = Vec::new();
        launched.try_reserve_exact(scratch_samples)
            .map_err(|_| invalid("flue callback scratch exceeds available allocation"))?;
        launched.resize(scratch_samples, 0.0);
        Ok(Self { voice, receiver, config, radiation, outlet_nodes, mechanics_rate_hz,
            mechanics_per_output, schedule, controls, next: 0, launched, mouth: None, flow_tap: None, poisoned: false })
    }

    /// Include actual mouth volume flow in the same external microphone signal.
    /// Mouth and terminal pressures superpose with their signs and distinct
    /// travel times; each retains all mechanical samples and antialias history.
    /// This selection alone retains the existing source load. Calling
    /// `with_mouth_reaction` explicitly replaces it with the same rectangle's
    /// compact exterior-gas reaction before any sample is rendered.
    ///
    /// This explicitly describes a coplanar baffled assembly. It does not claim
    /// the free-field geometry of an orchestral flute. Mouth and outlet openings
    /// must not overlap. Geometry is frozen; the admitted pressure/angle controls
    /// change excitation without changing either observed surface.
    pub fn with_mouth_radiation(
        mut self,
        input: FlueMouthRadiation,
    ) -> Result<Self, GestureCompileError> {
        if self.poisoned
            || self.voice.samples_rendered() != 0
            || self.mouth.is_some()
            || self.receiver.observer_gas().is_some()
        {
            return Err(invalid("mouth observation requires a fresh performance and exactly one source mouth"));
        }
        let network = self.voice.network().expect("admitted physical flue graph");
        let size = input.size_m(&self.voice.geometry);
        if size.iter().any(|x| !x.is_finite() || *x <= 0.0) {
            return Err(invalid("mouth aperture dimensions must be positive finite"));
        }
        let exterior = self.receiver.exterior_gas();
        if self.voice.gas().pressure.to_bits() != exterior.pressure.to_bits() {
            return Err(invalid("mouth and receiver static pressure differ; mean interface flow is not modeled"));
        }
        let inlet_node = network.graph().nodes.iter()
            .position(|node| matches!(node, fs_vfit::waveguide::network::NetworkNode::Inlet))
            .ok_or_else(|| invalid("mouth observation requires the original physical inlet"))?;
        for port in &self.receiver.source().openings {
            let section = network.graph().sections.iter().find(|s| s.nodes.contains(&port.node))
                .expect("admitted degree-one terminal");
            let dx = ((port.center_m[0] - input.center_m[0]).abs() - 0.5 * size[0]).max(0.0);
            let dy = ((port.center_m[1] - input.center_m[1]).abs() - 0.5 * size[1]).max(0.0);
            if !dx.is_finite() || !dy.is_finite() || dx.hypot(dy) < section.radius_m {
                return Err(invalid("mouth rectangle and radiating terminal disks must not overlap"));
            }
        }
        let geometry = RectangularOutletReceiver { center_m: input.center_m,
            position_m: self.config.receiver.position_m, cells: input.cells,
            maximum_frequency_hz: self.config.receiver.maximum_frequency_hz };
        let receiver = match self.receiver.source().reconstruction {
            DuctReconstruction::FiniteDifference { sample_rate_hz } =>
                BaffledPressure::rectangular_outlet(size, sample_rate_hz, geometry, exterior.into()),
            DuctReconstruction::Differentiating { mechanics_rate_hz, output_rate_hz } =>
                BaffledPressure::rectangular_outlet_multirate_differentiating(size,
                    mechanics_rate_hz, output_rate_hz, geometry, exterior.into()),
        }.map_err(|e| GestureCompileError::Render(RenderError::Voice(AcousticRealizeError::Nonlinear(e))))?;
        self.mouth = Some(MouthPressure { input, inlet_node, receiver, observer_gas: None,
            trace: vec![0.0; self.mechanics_per_output], filled: 0 });
        Ok(self)
    }

    /// Carry the explicit exterior gas law through the mouth and every terminal
    /// path. Attach all openings first. Each original surface retains its own
    /// distance, phase and absorption; a shared numerical padding keeps the
    /// separately prepared mouth and terminal pressures on the same time origin.
    /// This changes observation only, leaving the jet and bore state untouched.
    /// Omission retains the original lossless propagation and arithmetic.
    pub fn with_observer_gas(
        mut self,
        propagation: &GasPropagation,
        gate: &fs_exec::CancelGate,
    ) -> Result<Self, GestureCompileError> {
        if self.poisoned
            || self.voice.samples_rendered() != 0
            || self.receiver.observer_gas().is_some()
        {
            return Err(invalid("flue gas observation requires one cold attachment"));
        }
        let observation_error = |error| {
            GestureCompileError::Render(RenderError::Voice(AcousticRealizeError::Nonlinear(error)))
        };
        let mut receiver = self
            .receiver
            .with_observer_gas(propagation, &mut || gate.is_requested())
            .map_err(observation_error)?;
        let mut mouth_receiver = if let Some(mouth) = &self.mouth {
            Some(
                BaffledPressure::rectangular_outlet_multirate_gas_differentiating(
                    mouth.input.size_m(&self.voice.geometry),
                    self.mechanics_rate_hz,
                    self.config.sample_rate_hz,
                    RectangularOutletReceiver {
                        center_m: mouth.input.center_m,
                        position_m: self.config.receiver.position_m,
                        cells: mouth.input.cells,
                        maximum_frequency_hz: self.config.receiver.maximum_frequency_hz,
                    },
                    receiver.exterior_gas(),
                    propagation,
                    &mut || gate.is_requested(),
                )
                .map_err(observation_error)?,
            )
        } else {
            None
        };
        if let Some((mouth, report)) = &mut mouth_receiver {
            let latency = receiver
                .observer_gas()
                .expect("prepared terminal gas report")
                .additional_latency_samples
                .max(report.additional_latency_samples);
            receiver
                .align_gas_latency(latency)
                .map_err(observation_error)?;
            mouth
                .align_gas_latency(report, latency)
                .map_err(observation_error)?;
        }
        if gate.is_requested() {
            return Err(invalid("flue material observation cancelled"));
        }
        self.receiver = receiver;
        if let Some((receiver, report)) = mouth_receiver {
            let mouth = self.mouth.as_mut().expect("original physical mouth");
            mouth.receiver = receiver;
            mouth.observer_gas = Some(report);
        }
        Ok(self)
    }

    /// Replace the legacy mouth end correction and authored linear loss with
    /// compact self reaction from the SAME rectangle, exterior and source clock
    /// as the attached receiver. The load participates in the current jet/duct
    /// solve; it adds no signal or downstream pressure gain.
    ///
    /// Requires `mouth_radiation` and zero legacy `mouth_loss_rel`. The shared
    /// passive relaxation owner retains near-field storage separately from
    /// escaped energy; an independent separation law can remain attached.
    /// The checked compact band and quadrature limits belong to
    /// [`crate::flue_loop::mouth_radiation::RectangularMouthLoad`]. Mutual loading
    /// between mouth and terminal openings and free-field body scattering remain
    /// outside this independent self-load construction.
    ///
    /// # Errors
    /// Missing mouth observation, late/repeated attachment, nonzero legacy loss,
    /// mismatched source pressure or any original physical load refusal.
    pub fn with_mouth_reaction(mut self, gate: &fs_exec::CancelGate)
        -> Result<Self, GestureCompileError>
    {
        if self.poisoned || self.voice.samples_rendered() != 0 {
            return Err(invalid("mouth reaction requires a fresh performance"));
        }
        let mouth = self.mouth.as_ref()
            .ok_or_else(|| invalid("mouth reaction requires the actual rectangular mouth receiver"))?;
        let load = crate::flue_loop::mouth_radiation::RectangularMouthLoad::new(
            mouth.input.size_m(&self.voice.geometry), mouth.input.cells,
            *self.receiver.exterior_gas(), self.voice.sample_period_s(),
            self.config.receiver.maximum_frequency_hz, gate)
            .map_err(|e| GestureCompileError::Render(source_error(e)))?;
        self.voice.attach_mouth_radiation(load, mouth.input.size_m(&self.voice.geometry))
            .map_err(|e| GestureCompileError::Render(source_error(e)))?;
        Ok(self)
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
    /// Terminal receiver, including physical and numerical delay metadata.
    /// An attached mouth's separate contribution is exposed by `mouth_receiver`.
    #[must_use]
    pub const fn receiver(&self) -> &BaffledPressure { self.receiver.receiver() }
    /// Actual rectangular mouth receiver, absent when only terminals are selected.
    #[must_use]
    pub fn mouth_receiver(&self) -> Option<&BaffledPressure> { self.mouth.as_ref().map(|m| &m.receiver) }
    /// Source gas propagation and path reports for all terminal/tone-hole mouths.
    pub fn observer_gas(&self) -> Option<&GasPropagationReport> {
        self.receiver.observer_gas()
    }
    /// The same exterior law and common numerical padding on the source mouth.
    pub fn mouth_observer_gas(&self) -> Option<&GasPropagationReport> {
        self.mouth
            .as_ref()
            .and_then(|mouth| mouth.observer_gas.as_ref())
    }
    /// Original mouth geometry placement and surface resolution.
    #[must_use]
    pub fn mouth_radiation(&self) -> Option<FlueMouthRadiation> { self.mouth.as_ref().map(|m| m.input) }
    /// Complete propagation FIR support for every selected physical source.
    #[must_use]
    pub fn propagation_delay_samples(&self) -> (usize, usize) {
        let terminal = self.receiver().delay_samples;
        self.mouth_receiver().map_or(terminal, |m|
            (terminal.0.min(m.delay_samples.0), terminal.1.max(m.delay_samples.1)))
    }
    /// Original explicit exterior/environment input, absent for the legacy path.
    #[must_use]
    pub const fn radiation_input(&self) -> Option<&FlueRadiationInput> { self.radiation.as_ref() }
    /// Input order for `receiver().mode_pressures`: primary then additional ports.
    #[must_use]
    pub fn outlet_nodes(&self) -> &[usize] { &self.outlet_nodes }
    /// Actual source clock; independent of the PCM clock only when requested.
    #[must_use]
    pub const fn mechanics_rate_hz(&self) -> u32 { self.mechanics_rate_hz }
    /// Number of complete physical frames consumed per emitted audio sample.
    #[must_use]
    pub const fn mechanics_per_output(&self) -> usize { self.mechanics_per_output }
    /// OUTPUT samples remaining. Source sample counts remain available on `voice`.
    #[must_use]
    pub fn remaining_samples(&self) -> u64 { self.config.samples - self.samples_rendered() }

    fn advance_block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        let ratio = self.mechanics_per_output;
        let end = self.voice.samples_rendered() + output.len() as u64 * ratio as u64;
        let mut written = 0;
        // The shared receiver retains partial packets across control segments.
        // Each successful callback ends on a complete output frame.
        while self.voice.samples_rendered() < end {
            let now = self.voice.samples_rendered();
            while let Some(control) = self.controls.get(self.next) {
                if control.sample != now { break; }
                self.voice.apply(control.delta).map_err(source_error)?;
                self.next += 1;
            }
            let until = self.controls.get(self.next).map_or(end, |c| c.sample.min(end));
            let count = (until - now) as usize;
            let receiver = &mut self.receiver;
            let mouth = &mut self.mouth;
            let flow_tap = &mut self.flow_tap;
            let ceiling = self.config.maximum_abs_pressure_pa;
            let mut destinations = output[written..].iter_mut();
            self.voice.step_block_labium_observed(&mut self.launched[..count], |line, labium| {
                let FlueResonator::Network(network) = line else { unreachable!("owned physical graph") };
                if let Some(tap) = flow_tap.as_mut() { tap.record(network, labium)?; }
                let pressure = receiver.push_frame(|_, node| network.node_frame(node)
                    .map(|frame| frame.net_flow_into_node_m3_s))
                    .map_err(|e| FlueError::Physics(AcousticRealizeError::Nonlinear(e)))?;
                let mouth_pressure = match mouth {
                    Some(mouth) => {
                        let contribution = mouth.push(network)?;
                        if contribution.is_some() != pressure.is_some() {
                            return Err(FlueError::Invalid { what: "mouth and terminal observation clocks disagree" });
                        }
                        contribution
                    }
                    None => None,
                };
                let Some(pressure) = pressure else { return Ok(()); };
                let pressure = pressure + mouth_pressure.unwrap_or(0.0);
                if !pressure.is_finite() || pressure.abs() > ceiling {
                    return Err(FlueError::Invalid { what: "flue microphone pressure exceeds its declared finite ceiling" });
                }
                *destinations.next().expect("one observation per complete mechanical trace") = pressure;
                written += 1;
                Ok(())
            }).map_err(source_error)?;
        }
        debug_assert_eq!(self.receiver.pending_frames(), 0);
        debug_assert_eq!(written, output.len());
        Ok(())
    }
}

impl PressureRenderer for FlueOutletPerformance {
    fn samples_rendered(&self) -> u64 {
        self.voice.samples_rendered() / self.mechanics_per_output as u64
    }
    fn max_block_len(&self) -> usize { self.config.max_block }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.validate_sample_count(0)?;
        if rate != self.config.sample_rate_hz { return Err(sizing("flue pressure output must use its admitted observation clock")); }
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
