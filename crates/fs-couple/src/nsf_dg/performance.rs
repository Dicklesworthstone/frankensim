//! Finite physical-pressure playback of the original resolved zone and DG bore.
//! The caller supplies typed incoming-characteristic or pressure samples on the
//! OTHER resolved face. Pressure controls the original face reconstruction;
//! it is not converted to an authored velocity or calibrated source waveform.
//! The coupled face is still solved against the live DG history.
//! Initial resolved motion and thermochemical material state are retained.

use super::{NsfDg, NsfDgDrive, NsfDgFrame, same_clock};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::performance::ApertureObservation;
use crate::pcm_wav::baffled::BaffledPressure;
use crate::pcm_wav::baffled::gas::{GasPropagation, GasPropagationReport};
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use fs_exec::CancelGate;
use fs_lbm::lattice3::Lattice3;
use fs_lbm::nsf3::NsfBoundaryDrive;
use std::sync::Arc;

/// Generic, source-bound resolved geometry and finite boundary history.
pub mod file;
mod clock;
use clock::IntervalClock;
mod reservoir;
pub use reservoir::ReservoirObservation;

fn sizing(what: &'static str) -> RenderError { RenderError::Sizing { what } }

/// Number of original physical intervals needed to cover an output window.
/// Endpoint products determine coverage; neither geometry nor physical dt is
/// rounded to the observer clock. A derived same-clock pair retains its exact
/// one-interval-per-sample path and the existing clock-roundoff admission.
pub fn source_intervals_for(physical_dt_s: f64, output_rate_hz: u32, output_samples: u64)
    -> Result<usize, RenderError>
{
    if same_clock(physical_dt_s, 1.0 / f64::from(output_rate_hz)) {
        if output_samples == 0 || output_samples > clock::EXACT_INTEGER_LIMIT {
            return Err(sizing("resolved observer requires positive exactly representable sample indices"));
        }
        usize::try_from(output_samples).map_err(|_| sizing("resolved output sample count exceeds addressable source data"))
    } else { Ok(IntervalClock::new(physical_dt_s, output_rate_hz, output_samples)
        .map_err(sizing)?.source_intervals()) }
}

/// A finite source on the existing pressure, decimation, ensemble and PCM path.
/// Failed individual physical steps publish neither owner. A failed callback
/// poisons playback, retaining any already accepted prefix and observer state.
pub struct NsfDgPerformance<L: Lattice3<27>> {
    system: NsfDg<L>,
    boundary: Vec<NsfDgDrive>,
    observation: ApertureObservation,
    receiver: Option<BaffledPressure>,
    observer_gas: Option<GasPropagationReport>,
    reservoir: Option<reservoir::ReservoirPressure>,
    terminal_flows: Vec<f64>,
    sample_rate_hz: u32,
    maximum_block: usize,
    gate: Arc<CancelGate>,
    last: Option<NsfDgFrame>,
    clock: Option<IntervalClock>,
    projected_observer: bool,
    output_samples: u64,
    output_cursor: u64,
    poisoned: bool,
}

impl<L: Lattice3<27>> NsfDgPerformance<L> {
    /// Bind one held incoming-characteristic value per coupled interval.
    /// The sample rate matches the existing physical clock; explicit output
    /// decimation belongs to the shared PressureRenderer adapter.
    pub fn new(
        system: NsfDg<L>, observation: ApertureObservation, sample_rate_hz: u32,
        incoming_characteristics_m_s: Vec<f64>, maximum_block: usize,
        gate: Arc<CancelGate>,
    ) -> Result<Self, RenderError> {
        Self::new_with_boundary(system, observation, sample_rate_hz,
            incoming_characteristics_m_s.into_iter().map(NsfBoundaryDrive::IncomingCharacteristic).collect(),
            maximum_block, gate)
    }

    /// Bind typed physical controls on the original matching clock.
    /// Every value names its units and boundary law explicitly. Switching the
    /// drive kind preserves accepted physical fields through the original
    /// force-chart rebase; pressure and incoming-wave values are never mixed.
    pub fn new_with_boundary(
        system: NsfDg<L>, observation: ApertureObservation, sample_rate_hz: u32,
        boundary: Vec<NsfBoundaryDrive>, maximum_block: usize, gate: Arc<CancelGate>,
    ) -> Result<Self, RenderError> {
        Self::new_with_drives(system, observation, sample_rate_hz,
            boundary.into_iter().map(NsfDgDrive::Boundary).collect(), maximum_block, gate)
    }

    /// Bind boundary controls or explicit series pressure rises through a
    /// retained finite upstream load on the original physical clock.
    pub fn new_with_drives(
        system: NsfDg<L>, observation: ApertureObservation, sample_rate_hz: u32,
        boundary: Vec<NsfDgDrive>, maximum_block: usize, gate: Arc<CancelGate>,
    ) -> Result<Self, RenderError> {
        if !same_clock(system.load().propagator().dt(), 1.0 / f64::from(sample_rate_hz)) {
            return Err(sizing("resolved performance requires its original matching physical clock"));
        }
        let output_samples = boundary.len() as u64;
        Self::build(system, observation, sample_rate_hz, boundary,
            output_samples, maximum_block, gate, None, false)
    }

    /// Observe the unchanged physical solver on an independent regular clock.
    /// Incoming data has one value per ORIGINAL coupled interval; its length
    /// must equal `source_intervals_for`. Each original interval-mean pressure
    /// and terminal flow is held only for conservative overlap integration onto
    /// output intervals, then the original Rayleigh observer runs on that clock.
    ///
    /// This is a piecewise-constant temporal projection, with a box-average sinc
    /// response; it is not a band-limited anti-alias filter or an accuracy claim.
    /// Refine the physical and observation clocks for the intended bandwidth,
    /// and use the existing output decimator when needed. Physical state may be
    /// ahead of output time by less than one original interval; it is never
    /// slowed, retimed or reconstructed from observer values. Matching clocks
    /// retain the original one-frame arithmetic and callback behavior; a
    /// bit-identical clock delegates to `new` itself.
    #[allow(clippy::too_many_arguments)]
    pub fn new_resampled(
        system: NsfDg<L>, observation: ApertureObservation, output_rate_hz: u32,
        incoming_characteristics_m_s: Vec<f64>, output_samples: u64,
        maximum_block: usize, gate: Arc<CancelGate>,
    ) -> Result<Self, RenderError> {
        Self::new_resampled_with_boundary(system, observation, output_rate_hz,
            incoming_characteristics_m_s.into_iter().map(NsfBoundaryDrive::IncomingCharacteristic).collect(),
            output_samples, maximum_block, gate)
    }

    /// Observe typed physical boundary controls on an independent regular clock.
    /// The original NSF interval and conservative overlap projection are the
    /// same as in the incoming-characteristic constructor. Source histories
    /// have one value per physical interval and never follow callback timing.
    #[allow(clippy::too_many_arguments)]
    pub fn new_resampled_with_boundary(
        system: NsfDg<L>, observation: ApertureObservation, output_rate_hz: u32,
        boundary: Vec<NsfBoundaryDrive>, output_samples: u64,
        maximum_block: usize, gate: Arc<CancelGate>,
    ) -> Result<Self, RenderError> {
        Self::new_resampled_with_drives(system, observation, output_rate_hz,
            boundary.into_iter().map(NsfDgDrive::Boundary).collect(),
            output_samples, maximum_block, gate)
    }

    /// Observe the same original two-sided or one-sided physical system without
    /// retiming source histories or replacing its declared drive location.
    #[allow(clippy::too_many_arguments)]
    pub fn new_resampled_with_drives(
        system: NsfDg<L>, observation: ApertureObservation, output_rate_hz: u32,
        boundary: Vec<NsfDgDrive>, output_samples: u64,
        maximum_block: usize, gate: Arc<CancelGate>,
    ) -> Result<Self, RenderError> {
        let dt = system.load().propagator().dt();
        if boundary.len() != source_intervals_for(dt, output_rate_hz, output_samples)? {
            return Err(sizing("resolved boundary data must cover exactly the original physical intervals of its output window"));
        }
        if same_clock(dt, 1.0 / f64::from(output_rate_hz)) {
            if dt.to_bits() == (1.0 / f64::from(output_rate_hz)).to_bits() {
                return Self::new_with_drives(system, observation, output_rate_hz,
                    boundary, maximum_block, gate);
            }
            // Existing derived-clock roundoff admission still uses exactly
            // one physical frame per output interval. Only the observer's
            // source-band check differs from the bit-exact legacy API.
            return Self::build(system, observation, output_rate_hz,
                boundary, output_samples, maximum_block, gate, None, true);
        }
        let clock = IntervalClock::new(dt, output_rate_hz, output_samples).map_err(sizing)?;
        Self::build(system, observation, output_rate_hz, boundary,
            output_samples, maximum_block, gate, Some(clock), true)
    }

    #[allow(clippy::too_many_arguments)]
    fn build(
        system: NsfDg<L>, observation: ApertureObservation, sample_rate_hz: u32,
        boundary: Vec<NsfDgDrive>, output_samples: u64,
        maximum_block: usize, gate: Arc<CancelGate>, clock: Option<IntervalClock>,
        projected_observer: bool,
    ) -> Result<Self, RenderError> {
        if gate.is_requested() || system.accepted_steps() != 0
            || sample_rate_hz == 0 || !(1..=65_536).contains(&maximum_block)
            || boundary.is_empty()
            || boundary.len() as u64 > system.limits().maximum_steps
            || !boundary.iter().all(|x| x.is_finite())
            || boundary.iter().any(|drive| system.upstream_load().is_some()
                != matches!(drive, NsfDgDrive::UpstreamPressureRise(_)))
            || output_samples == 0
        {
            return Err(sizing("resolved performance requires a cold matching clock and finite boundary samples/capacity"));
        }
        let receiver = match observation {
            ApertureObservation::Inlet => None,
            ApertureObservation::TubeBaffled(location) => Some(if projected_observer {
                system.load().baffled_receiver_projected(sample_rate_hz, location)
            } else { system.load().baffled_receiver(sample_rate_hz, location) }
                .map_err(RenderError::Voice)?),
            _ => return Err(sizing("resolved pressure observation must name its DG inlet or computed baffled outlet")),
        };
        let terminal_flows = if receiver.is_some() {
            vec![0.0; system.load().terminal_ports().len()]
        } else { Vec::new() };
        Ok(Self {
            system, boundary, observation, receiver,
            observer_gas: None, reservoir: None, terminal_flows, sample_rate_hz, maximum_block, gate,
            last: None, clock, projected_observer, output_samples, output_cursor: 0, poisoned: false,
        })
    }
    /// Bind the shared physical atmosphere law to every original outlet-to-microphone
    /// path before playback. The exterior gas comes from the actual DG radiation
    /// source, including its transport and caloric properties. The original fluid,
    /// bore, terminal reactions and held physical boundary samples remain unchanged.
    ///
    /// Source absorption, dispersion and causal delay use the same Rayleigh owner
    /// as native flute and reed playback. No second filter is applied after the
    /// surface pressure is summed. The retained cancellation gate also owns this
    /// bounded cold preparation; failed preparation publishes no receiver.
    pub fn with_observer_gas(mut self, propagation: &GasPropagation) -> Result<Self, RenderError> {
        if self.poisoned || self.system.accepted_steps() != 0 || self.observer_gas.is_some() {
            return Err(sizing("resolved material gas observation requires one cold attachment"));
        }
        let ApertureObservation::TubeBaffled(location) = self.observation else {
            return Err(sizing("resolved material gas propagation requires an exterior receiver"));
        };
        let (mut receiver, mut report) = if self.projected_observer {
            self.system.load().baffled_receiver_projected_with_gas(
                self.sample_rate_hz, location, propagation, &self.gate)
        } else { self.system.load().baffled_receiver_with_gas(
            self.sample_rate_hz, location, propagation, &self.gate)
        }.map_err(RenderError::Voice)?;
        let mut reservoir_receiver = self.reservoir.as_ref().map(|reservoir| {
            let gas = self.system.load().radiation_record().expect("original exterior").gas_record().state();
            BaffledPressure::new_multirate_gas_differentiating(&reservoir.surface,
                reservoir.cells.len(), self.sample_rate_hz, self.sample_rate_hz,
                location.position_m, gas, propagation, &mut || self.gate.is_requested())
                .map_err(|error| RenderError::Voice(AcousticRealizeError::Nonlinear(error)))
        }).transpose()?;
        if let Some((escape, escape_report)) = &mut reservoir_receiver {
            let latency = report.additional_latency_samples.max(escape_report.additional_latency_samples);
            receiver.align_gas_latency(&mut report, latency)
                .map_err(|error| RenderError::Voice(AcousticRealizeError::Nonlinear(error)))?;
            escape.align_gas_latency(escape_report, latency)
                .map_err(|error| RenderError::Voice(AcousticRealizeError::Nonlinear(error)))?;
        }
        if self.gate.is_requested() { return Err(sizing("resolved material observation cancelled")); }
        self.receiver = Some(receiver);
        self.observer_gas = Some(report);
        if let Some((receiver, report)) = reservoir_receiver {
            let reservoir = self.reservoir.as_mut().expect("original escape aperture");
            reservoir.receiver = receiver;
            reservoir.gas = Some(report);
        }
        Ok(self)
    }
    /// Original exterior material, selected physical law and individual path responses.
    pub fn observer_gas(&self) -> Option<&GasPropagationReport> { self.observer_gas.as_ref() }

    /// Original resolved source, material state and acoustic load.
    pub const fn system(&self) -> &NsfDg<L> { &self.system }
    /// Original physical controls, with explicit location, kinds and units.
    pub fn drive_history(&self) -> &[NsfDgDrive] { &self.boundary }
    /// Last accepted physical interval, even if its later observation failed.
    pub const fn last_frame(&self) -> Option<&NsfDgFrame> { self.last.as_ref() }
    /// Physical location of the emitted pressure.
    pub const fn observation(&self) -> ApertureObservation { self.observation }
    /// Original reaction-bound retarded receiver and its explicit delay.
    pub fn baffled_receiver(&self) -> Option<&BaffledPressure> { self.receiver.as_ref() }
    /// Shared cancellation owner, also polled during each nonlinear interval.
    pub fn cancellation_gate(&self) -> &CancelGate { &self.gate }
    /// Remaining observer intervals. Original physical progress remains on
    /// `system().accepted_steps()` and can differ under temporal projection.
    pub fn remaining_samples(&self) -> u64 {
        self.output_samples - self.output_cursor
    }

    fn advance_source(&mut self) -> Result<(), RenderError> {
        let boundary = *self.boundary.get(self.system.accepted_steps() as usize)
            .ok_or_else(|| sizing("resolved physical boundary history is exhausted"))?;
        let frame = self.system.step_with_drive(boundary, &self.gate)
            .map_err(|e| RenderError::Voice(AcousticRealizeError::Nonlinear(e.to_string())))?;
        self.last = Some(frame);
        Ok(())
    }

    fn advance(&mut self) -> Result<f64, RenderError> {
        if let Some(clock) = self.clock { return self.advance_projected(clock); }
        self.advance_source()?;
        // Preserve the strict path's accepted-prefix semantics even if the
        // subsequent observer fails. Source and output indices coincide here.
        self.output_cursor += 1;
        let frame = self.last.as_ref().expect("accepted original resolved interval");
        if let Some(reservoir) = &mut self.reservoir {
            reservoir.flows.fill(0.0);
            reservoir.accumulate(&frame.reservoir_flows, 1.0)?;
        }
        if self.receiver.is_some() {
            if frame.terminals.len() != self.terminal_flows.len() || frame.terminals.is_empty() {
                return Err(sizing("resolved DG microphone requires every accepted terminal flow"));
            }
            for (flow, terminal) in self.terminal_flows.iter_mut().zip(&frame.terminals) {
                *flow = terminal.outward_flow_m3_s;
            }
        }
        self.observe(frame.inlet_pressure_pa)
    }

    fn advance_projected(&mut self, clock: IntervalClock) -> Result<f64, RenderError> {
        if self.gate.is_requested() {
            return Err(RenderError::Voice(AcousticRealizeError::Nonlinear("resolved observation cancelled".into())));
        }
        let (mut cursor, end) = clock.output_bounds(self.output_cursor).map_err(sizing)?;
        let mut pressure = 0.0;
        self.terminal_flows.fill(0.0);
        if let Some(reservoir) = &mut self.reservoir { reservoir.flows.fill(0.0); }
        while cursor < end {
            while self.last.is_none() || clock.source_end(self.system.accepted_steps()) <= cursor {
                self.advance_source()?;
            }
            let (weight, next) = clock.overlap(cursor, end, self.system.accepted_steps()).map_err(sizing)?;
            let frame = self.last.as_ref().expect("accepted covering physical interval");
            if let Some(reservoir) = &mut self.reservoir {
                reservoir.accumulate(&frame.reservoir_flows, weight)?;
            }
            pressure += weight * frame.inlet_pressure_pa;
            if self.receiver.is_some() {
                if frame.terminals.len() != self.terminal_flows.len() || frame.terminals.is_empty() {
                    return Err(sizing("resolved projected microphone requires every original terminal flow"));
                }
                for (flow, terminal) in self.terminal_flows.iter_mut().zip(&frame.terminals) {
                    *flow += weight * terminal.outward_flow_m3_s;
                }
            }
            cursor = next;
        }
        self.output_cursor += 1;
        self.observe(pressure)
    }

    fn observe(&mut self, inlet_pressure: f64) -> Result<f64, RenderError> {
        let mut pressure = if let Some(receiver) = &mut self.receiver {
            receiver.step_trace(&self.terminal_flows)
                .map_err(|e| RenderError::Voice(AcousticRealizeError::Nonlinear(e)))?
        } else { inlet_pressure };
        if let Some(reservoir) = &mut self.reservoir {
            pressure += reservoir.receiver.step_trace(&reservoir.flows)
                .map_err(|error| RenderError::Voice(AcousticRealizeError::Nonlinear(error)))?;
        }
        if !pressure.is_finite() {
            return Err(sizing("resolved observer produced nonfinite physical pressure"));
        }
        Ok(pressure)
    }
}

impl<L: Lattice3<27>> PressureRenderer for NsfDgPerformance<L> {
    fn samples_rendered(&self) -> u64 { self.output_cursor }
    fn max_block_len(&self) -> usize { self.maximum_block }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.validate_sample_count(0)?;
        if rate != self.sample_rate_hz {
            return Err(sizing("resolved pressure uses its admitted observer clock"));
        }
        Ok(())
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if samples > self.remaining_samples() {
            return Err(sizing("resolved pressure exceeds its finite source window"));
        }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.validate_sample_count(output.len() as u64)?;
        if output.is_empty() { return Err(RenderError::EmptyBlock); }
        if output.len() > self.maximum_block {
            return Err(sizing("resolved pressure callback exceeds its declared capacity"));
        }
        for slot in output {
            match self.advance() {
                Ok(value) => *slot = value,
                Err(error) => { self.poisoned = true; return Err(error); }
            }
        }
        Ok(())
    }
}
