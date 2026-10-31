//! Existing retarded baffled pressure and streaming PCM interface over BOTH loads.
//! Observation consumes terminal volume flow once per accepted source sample.
//! No new microphone law, output normalization, pressure-to-note mapping or
//! radiation reaction is inferred. Declare radiative mechanical loads in the ducts.
use super::{LoadedAperture, TwoSidedAperture, TwoSidedDrive, TwoSidedFrame, ValveSide};
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver};
use crate::pcm_wav::baffled::gas::{GasPropagation, GasPropagationReport};
use crate::pcm_wav::observation::{ObservationError, ObservationPolicy, PressureRenderer};
use crate::bernoulli_aperture::radiation::environment::BaffledOutletBinding;
use crate::render::RenderError;
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::network::boundaries::{BoundaryAperture, BoundaryDrive};
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
    held_spatial_pressures_pa: [f64; 2],
    gate: CancelGate,
    completed: u64,
    last: Option<TwoSidedFrame>,
    poisoned: bool,
    observation_policy: ObservationPolicy,
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
    /// If a side already owns a MIMO exterior, every outlet in that ordered group
    /// must appear with the same geometry, gas and mechanics clock. Validation
    /// remaps global side indices to the group's local graph without relabeling
    /// its atmosphere or operator. This also applies to diagnostic playback.
    ///
    /// # Errors
    /// Late attachment, clock mismatch, zero callback bound, invalid ports/gas,
    /// or original receiver geometry, bandwidth and history-budget admission.
    pub fn with_radiation(
        system: TwoSidedAperture,
        radiation: DuctRadiationSpec,
        max_block: usize,
    ) -> Result<Self, RenderError> {
        Self::with_observation_policy(
            system,
            radiation,
            max_block,
            ObservationPolicy::Diagnostic,
            &[],
        )
    }

    /// Acceptance requires one retained, actual reaction binding per opening.
    /// The raw one-way constructor remains an explicitly diagnostic observer.
    pub fn with_observation_policy(
        system: TwoSidedAperture,
        radiation: DuctRadiationSpec,
        max_block: usize,
        policy: ObservationPolicy,
        bindings: &[BaffledOutletBinding],
    ) -> Result<Self, RenderError> {
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
        system.validate_mimo_receiver(receiver.source(), *receiver.exterior_gas())
            .map_err(RenderError::Voice)?;
        if policy == ObservationPolicy::Acceptance || !bindings.is_empty() {
            let source = receiver.source();
            if bindings.len() != source.openings.len() {
                return Err(observation(
                    ObservationError::UnmatchedReactionLoad.to_string(),
                ));
            }
            for opening in &source.openings {
                let mut found = bindings.iter().filter(|b| b.opening() == *opening);
                let binding = found.next().ok_or_else(|| {
                    observation(ObservationError::UnmatchedReactionLoad.to_string())
                })?;
                if found.next().is_some() {
                    return Err(observation(
                        ObservationError::UnmatchedReactionLoad.to_string(),
                    ));
                }
                binding
                    .validate_receiver(
                        &system.sides[opening.network].graph,
                        source,
                        *receiver.exterior_gas(),
                    )
                    .map_err(|e| observation(e.to_string()))?;
            }
        }
        Ok(Self {
            system,
            receiver,
            rate,
            max_block,
            held: TwoSidedDrive::default(),
            held_spatial_pressures_pa: [0.0; 2],
            gate: CancelGate::new_clock_free(),
            completed: 0,
            last: None,
            poisoned: false,
            observation_policy: policy,
        })
    }
    /// Admission scope, retained independently of the source's mechanical clock.
    pub const fn observation_policy(&self) -> ObservationPolicy {
        self.observation_policy
    }
    /// Bind the shared material propagation law to every original aperture path
    /// before either source or receiver advances. The same exterior gas, signed
    /// port ordering, numerical clocks and mechanical radiation loads remain.
    /// Absorption and dispersion act before pressure summation; the receiver's
    /// residual travel and common causal padding count their delays exactly once.
    ///
    /// This consumes the cold source atomically. Original gas-law, work-budget,
    /// cancellation and differentiating-reconstruction refusals propagate.
    pub fn with_observer_gas(
        mut self,
        propagation: &GasPropagation,
        gate: &CancelGate,
    ) -> Result<Self, RenderError> {
        if self.poisoned || self.mechanical_samples_rendered() != 0
            || self.completed != 0 || self.receiver.mechanical_samples() != 0 {
            return Err(sizing("material gas propagation requires an unadvanced two-sided source and receiver"));
        }
        self.receiver = self.receiver.with_observer_gas(propagation,
            &mut || gate.is_requested()).map_err(observation)?;
        Ok(self)
    }
    /// Original exterior material and individual path-filter accounts. Numerical
    /// padding also appears in `receiver().filter_delay_samples`; physical flight
    /// and molecular dispersion are never removed by ensemble alignment.
    pub fn observer_gas(&self) -> Option<&GasPropagationReport> {
        self.receiver.observer_gas()
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
    /// Set independently prescribed spatial blade pressure amplitudes [Pa].
    /// The material port retains its complete signed nodal trace. This does not
    /// change the uniform series source or either acoustic inlet's swept volume.
    pub fn set_spatial_pressures(&mut self, pressures_pa: [f64; 2]) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        self.system.valve.validate_pressures(pressures_pa).map_err(RenderError::Voice)?;
        self.held_spatial_pressures_pa = pressures_pa;
        Ok(())
    }
    /// Last selected spatial pressure amplitudes; scheduled controls publish
    /// their accepted values separately after successful physical ingestion.
    pub const fn spatial_pressures_pa(&self) -> [f64; 2] { self.held_spatial_pressures_pa }
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

    /// Attach material moving boundaries before either mechanics or reception
    /// advances. The original receiver and both graph histories stay intact.
    /// A scheduled performer must then bind a complete explicit drive program.
    pub fn with_boundary_apertures(mut self, side: ValveSide,
        inputs: Vec<BoundaryAperture>, max_apertures: usize) -> Result<Self, RenderError>
    {
        if self.poisoned || self.mechanical_samples_rendered() != 0 || self.completed != 0 {
            return Err(sizing("moving boundaries require an unadvanced two-sided receiver"));
        }
        self.system = self.system.with_boundary_apertures(side, inputs, max_apertures, &self.gate)
            .map_err(RenderError::Voice)?;
        Ok(self)
    }

    // The ONLY mechanical/observer step for held and scheduled playback. The
    // scheduler can change pressure between these calls inside an audio packet.
    fn advance_mechanical(&mut self) -> Result<Option<f64>, RenderError> {
        self.advance_mechanical_with_boundaries([&[], &[]])
    }

    fn advance_mechanical_with_boundaries(&mut self, drives: [&[BoundaryDrive]; 2])
        -> Result<Option<f64>, RenderError>
    {
        if self.poisoned { return Err(RenderError::Poisoned); }
        let result = (|| {
            let frame = self.system.step_with_boundaries_and_pressures_under_gate(
                self.held, drives, self.held_spatial_pressures_pa, &self.gate)
                .map_err(RenderError::Voice)?;
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
