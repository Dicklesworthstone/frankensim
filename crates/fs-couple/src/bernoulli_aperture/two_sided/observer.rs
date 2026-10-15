//! Existing retarded baffled pressure and streaming PCM interface over BOTH loads.
//! Observation consumes terminal volume flow once per accepted source sample.
//! No new microphone law, output normalization, pressure-to-note mapping or
//! radiation reaction is inferred. Declare radiative mechanical loads in the ducts.
use super::{LoadedAperture, TwoSidedAperture, TwoSidedDrive, TwoSidedFrame, ValveSide};
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver};
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use crate::acoustic_realize::AcousticRealizeError;
use fs_exec::CancelGate;
use crate::pcm_wav::baffled::openings::duct::{
    DuctOpening, DuctPressure, DuctRadiationSpec, DuctReconstruction, DuctView,
};

/// Canonical source-pressure phrases on the retained two-sided physics.
pub mod schedule;

/// Bounded pressure callbacks over one retained valve and two acoustic networks.
/// A physical/observer failure poisons the renderer: a partly completed callback
/// cannot be resumed as though its mechanics and pressure histories still agree.
/// The existing stream writer owns cancellation BETWEEN callbacks. Use the
/// mechanical API directly for per-sample cancellation and detailed diagnostics.
pub struct TwoSidedPressure {
    system: TwoSidedAperture,
    receiver: DuctPressure,
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
        Self::with_radiation(system, DuctRadiationSpec {
            openings: vec![DuctOpening { network: side.index(), node, center_m: [0.0; 2] }],
            exterior: None, receiver,
            reconstruction: DuctReconstruction::FiniteDifference { sample_rate_hz },
        }, max_block)
    }

    /// Bind multiple physical outlets from either acoustic side to the same
    /// shared gas-aware receiver as jet-driven ducts. Network 0 is upstream;
    /// network 1 is downstream. Radius and gas come from each actual graph.
    /// `Differentiating` consumes every source substep before emitting audio;
    /// the actual valve clock MUST match its declared mechanics rate exactly.
    ///
    /// Geometry, blade materials/contact histories, both gas fields, fractional
    /// delays and original terminal loads remain with the existing source. No
    /// valve is rebuilt, retuned or replaced by a sampled pressure track. The
    /// shared infinite-baffle observer does not introduce radiation backreaction.
    ///
    /// # Errors
    /// Late attachment, clock mismatch, zero callback bound, invalid ports/gas,
    /// or original receiver geometry, bandwidth and history-budget admission.
    pub fn with_radiation(system: TwoSidedAperture, radiation: DuctRadiationSpec,
        max_block: usize) -> Result<Self, RenderError>
    {
        let mechanics = radiation.reconstruction.mechanics_rate_hz();
        let rate = radiation.reconstruction.output_rate_hz();
        if system.valve.accepted_steps() != 0 || mechanics == 0 || max_block == 0
            || (1.0/f64::from(mechanics)).to_bits() != system.valve.time_step_s().to_bits() {
            return Err(sizing("two-sided receiver needs an unadvanced source, matching mechanical rate and positive callback bound"));
        }
        let receiver = DuctPressure::new(&[
            DuctView { graph: &system.sides[0].graph, gases: &system.sides[0].gases },
            DuctView { graph: &system.sides[1].graph, gases: &system.sides[1].gases },
        ], radiation).map_err(observation)?;
        Ok(Self { system, receiver, rate, max_block,
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
    pub const fn receiver(&self) -> &BaffledPressure { self.receiver.receiver() }
    /// Exact source-independent outlet/environment/reconstruction declaration.
    pub const fn radiation(&self) -> &DuctRadiationSpec { self.receiver.source() }
    /// Original physical clock, which may be higher than the audio clock.
    pub const fn mechanics_rate_hz(&self) -> u32 {
        self.receiver.source().reconstruction.mechanics_rate_hz()
    }
    /// Required physical steps per completely emitted microphone sample.
    pub const fn mechanics_per_output(&self) -> usize { self.receiver.mechanics_per_output() }
    /// Accepted physical steps; independent of completed output sample count.
    pub fn mechanical_samples_rendered(&self) -> u64 { self.system.valve.accepted_steps() }
    /// Last successfully ingested mechanical frame, absent before stepping.
    /// After a failed callback this is diagnostic, not a resumable output boundary.
    #[must_use]
    pub const fn last_frame(&self) -> Option<&TwoSidedFrame> { self.last.as_ref() }

    // The ONLY mechanical/observer step for held and scheduled playback. The
    // scheduler can change pressure between these calls inside an audio packet.
    fn advance_mechanical(&mut self) -> Result<Option<f64>, RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        let result = (|| {
            let frame = self.system.step_under_gate(self.held, &self.gate).map_err(RenderError::Voice)?;
            let system = &self.system;
            let pressure = self.receiver.push_frame(|network, node| {
                let side = match network { 0 => ValveSide::Upstream, 1 => ValveSide::Downstream, _ => return None };
                system.node_frame(side, node).map(|f| f.net_flow_into_node_m3_s)
            }).map_err(observation)?;
            Ok::<_, RenderError>((frame, pressure))
        })();
        match result {
            Ok((frame, pressure)) => {
                self.last = Some(frame);
                if pressure.is_some() { self.completed += 1; }
                Ok(pressure)
            }
            Err(error) => { self.poisoned = true; Err(error) }
        }
    }
}
impl PressureRenderer for TwoSidedPressure {
    fn samples_rendered(&self) -> u64 { self.completed }
    fn max_block_len(&self) -> usize { self.max_block }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if rate != self.rate { return Err(sizing("two-sided pressure sink must use the admitted output rate")); }
        Ok(())
    }
    fn validate_sample_count(&self, count: u64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if count.checked_mul(self.mechanics_per_output() as u64)
            .and_then(|steps| self.system.valve.accepted_steps().checked_add(steps))
            .is_none_or(|end| end > self.system.valve.max_steps()) {
            return Err(sizing("two-sided pressure request exceeds the source lifetime budget"));
        }
        Ok(())
    }
    fn block(&mut self, out: &mut [f64]) -> Result<(), RenderError> {
        if out.is_empty() || out.len() > self.max_block { return Err(sizing("two-sided callback must fit its positive capacity")); }
        let count = u64::try_from(out.len()).map_err(|_| sizing("pressure sample count overflow"))?;
        self.validate_sample_count(count)?;
        for slot in out {
            for _ in 0..self.mechanics_per_output() {
                if let Some(pressure) = self.advance_mechanical()? { *slot = pressure; }
            }
            debug_assert_eq!(self.receiver.pending_frames(), 0);
        }
        Ok(())
    }
}
