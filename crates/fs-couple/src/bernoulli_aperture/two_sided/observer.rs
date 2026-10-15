//! Existing retarded baffled pressure and streaming PCM interface over BOTH loads.
//! Observation consumes terminal volume flow once per accepted source sample.
//! No new microphone law, output normalization, pressure-to-note mapping or
//! radiation reaction is inferred. Declare radiative mechanical loads in the ducts.
use super::{LoadedAperture, NetworkNode, TwoSidedAperture, TwoSidedDrive, TwoSidedFrame, ValveSide};
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use crate::acoustic_realize::AcousticRealizeError;
use fs_exec::CancelGate;

/// Bounded pressure callbacks over one retained valve and two acoustic networks.
/// A physical/observer failure poisons the renderer: a partly completed callback
/// cannot be resumed as though its mechanics and pressure histories still agree.
/// The existing stream writer owns cancellation BETWEEN callbacks. Use the
/// mechanical API directly for per-sample cancellation and detailed diagnostics.
pub struct TwoSidedPressure {
    system: TwoSidedAperture,
    side: ValveSide,
    node: usize,
    receiver: BaffledPressure,
    rate: u32,
    max_block: usize,
    held: TwoSidedDrive,
    gate: CancelGate,
    completed: u64,
    last: Option<TwoSidedFrame>,
    poisoned: bool,
}
fn sizing(what: &'static str) -> RenderError { RenderError::Sizing { what } }
fn observation(error: String) -> RenderError {
    RenderError::Voice(AcousticRealizeError::Nonlinear(error))
}
impl TwoSidedPressure {
    /// Bind an explicitly placed circular-baffle receiver to a physical outlet.
    /// Its radius and exterior gas come from that outlet's ACTUAL section, not
    /// the source inlet. Its prehistory is zero, so an already-running mechanical
    /// source refuses. All source and receiver state then stay exclusively owned.
    pub fn new(system: TwoSidedAperture, side: ValveSide, node: usize,
        receiver: CircularOutletReceiver, sample_rate_hz: u32, max_block: usize)
        -> Result<Self, RenderError>
    {
        if system.valve.accepted_steps() != 0 || sample_rate_hz == 0 || max_block == 0
            || (1.0/f64::from(sample_rate_hz)).to_bits() != system.valve.time_step_s().to_bits() {
            return Err(sizing("two-sided receiver needs an unadvanced source, matching rate and positive callback bound"));
        }
        let s = &system.sides[side.index()];
        if !matches!(s.graph.nodes.get(node), Some(NetworkNode::Termination { .. }
            | NetworkNode::Impedance { .. } | NetworkNode::Relaxation { .. })) {
            return Err(sizing("two-sided exterior observation requires a passive degree-one terminal"));
        }
        let mut connected = s.graph.sections.iter().enumerate().filter(|(_, section)| section.nodes.contains(&node));
        let (index, section) = connected.next().ok_or_else(|| sizing("outlet has no physical section"))?;
        if connected.next().is_some() { return Err(sizing("outlet has multiple physical sections")); }
        let gas = s.gases[index];
        let receiver = BaffledPressure::circular_outlet(section.radius_m, sample_rate_hz, receiver,
            RayleighMedium { density: gas.density, sound_speed: gas.sound_speed }).map_err(observation)?;
        Ok(Self { system, side, node, receiver, rate: sample_rate_hz, max_block,
            held: TwoSidedDrive::default(), gate: CancelGate::new_clock_free(), completed: 0,
            last: None, poisoned: false })
    }
    /// Set physical supplies between callbacks; zero releases the source/actuator,
    /// not the retained vibration. A refused input changes neither held control
    /// nor any physical state. There is no frequency or output-amplitude control.
    pub fn set_drive(&mut self, drive: TwoSidedDrive) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if !drive.source_pressure_pa.is_finite()
            || drive.blade_forces_n.iter().any(|x| !x.is_finite())
            || (matches!(&self.system.valve, LoadedAperture::Single(_)) && drive.blade_forces_n[1] != 0.0) {
            return Err(RenderError::Control { what: "two-sided supplies must be finite and address existing force ports" });
        }
        self.held = drive;
        Ok(())
    }
    /// Read physical states and all original material/geometry/gas declarations.
    #[must_use]
    pub const fn system(&self) -> &TwoSidedAperture { &self.system }
    /// Actual spatial/travel-time observer, including uncompensated delay metadata.
    #[must_use]
    pub const fn receiver(&self) -> &BaffledPressure { &self.receiver }
    /// Last completely observed mechanical frame, absent before first output.
    #[must_use]
    pub const fn last_frame(&self) -> Option<&TwoSidedFrame> { self.last.as_ref() }
}
impl PressureRenderer for TwoSidedPressure {
    fn samples_rendered(&self) -> u64 { self.completed }
    fn max_block_len(&self) -> usize { self.max_block }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if rate != self.rate { return Err(sizing("two-sided pressure sink must use the admitted source rate")); }
        Ok(())
    }
    fn validate_sample_count(&self, count: u64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if self.completed.checked_add(count).is_none_or(|end| end > self.system.valve.max_steps()) {
            return Err(sizing("two-sided pressure request exceeds the source lifetime budget"));
        }
        Ok(())
    }
    fn block(&mut self, out: &mut [f64]) -> Result<(), RenderError> {
        if out.is_empty() || out.len() > self.max_block { return Err(sizing("two-sided callback must fit its positive capacity")); }
        let count = u64::try_from(out.len()).map_err(|_| sizing("pressure sample count overflow"))?;
        self.validate_sample_count(count)?;
        for slot in out {
            let result = (|| {
                let frame = self.system.step_under_gate(self.held, &self.gate).map_err(RenderError::Voice)?;
                let flow = self.system.node_frame(self.side, self.node)
                    .expect("admitted outlet cannot disappear").net_flow_into_node_m3_s;
                let pressure = self.receiver.step(&[flow]).map_err(observation)?;
                Ok::<_, RenderError>((frame, pressure))
            })();
            match result {
                Ok((frame, pressure)) => {
                    *slot = pressure;
                    self.last = Some(frame);
                    self.completed += 1;
                }
                Err(error) => { self.poisoned = true; return Err(error); }
            }
        }
        Ok(())
    }
}
