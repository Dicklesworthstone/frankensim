//! Source-independent binding of physical duct ports to the shared receiver.
//! Jets and single/paired valves provide accepted node flows; this owns only
//! geometry admission and complete mechanics-to-audio observation packets.
use super::{
    BaffledPressure, CircularOpening, CircularOutletReceiver, MAX_RAYLEIGH_INPUTS, RayleighMedium,
};
use crate::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec};
use crate::pcm_wav::baffled::gas::{GasPropagation, GasPropagationReport};
use fs_material::gas::GasState;

/// Borrow the realized geometry and its matching ordered regional gas states.
/// Node addresses refer to this graph, including any already lowered losses.
#[derive(Clone, Copy, Debug)]
pub struct DuctView<'a> {
    /// Physical cylindrical sections and their passive node laws.
    pub graph: &'a TubeNetworkSpec,
    /// Exactly one coherent gas per section, not a second authored density.
    pub gases: &'a [GasState],
}

/// A physical terminal, not an oscillator or a pressure-to-audio gain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctOpening {
    /// Index in the ordered graph views supplied to `DuctPressure::new`.
    pub network: usize,
    /// Terminal node in that graph. Equal node numbers in different graphs differ.
    pub node: usize,
    /// Center in a common stationary z=0 baffle [m]; all openings face +z.
    pub center_m: [f64; 2],
}

/// Existing numerical observation images; neither changes the mechanical source.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DuctReconstruction {
    /// Original one-clock finite-difference acceleration and linear travel delay.
    FiniteDifference { sample_rate_hz: u32 },
    /// Existing bandlimited derivative, anti-alias stages and causal travel delay.
    /// Every physical substep is consumed; an integer ratio 1..=16 is required.
    Differentiating { mechanics_rate_hz: u32, output_rate_hz: u32 },
}
impl DuctReconstruction {
    /// Actual source clock, not a rate inferred by rounding a floating dt.
    pub const fn mechanics_rate_hz(self) -> u32 {
        match self { Self::FiniteDifference { sample_rate_hz } => sample_rate_hz,
            Self::Differentiating { mechanics_rate_hz, .. } => mechanics_rate_hz }
    }
    /// Pressure output clock. Numerical filter latency is expressed at this rate.
    pub const fn output_rate_hz(self) -> u32 {
        match self { Self::FiniteDifference { sample_rate_hz } => sample_rate_hz,
            Self::Differentiating { output_rate_hz, .. } => output_rate_hz }
    }
}

/// One original geometric/environmental observation declaration.
#[derive(Clone, Debug, PartialEq)]
pub struct DuctRadiationSpec {
    /// Complete input ordering, with no duplicate physical ports.
    pub openings: Vec<DuctOpening>,
    /// Common exterior gas. None is allowed only for a single opening and uses
    /// its actual local gas; multiple openings must explicitly choose one medium.
    pub exterior: Option<GasState>,
    /// Position in the common baffle frame, resolution and admitted frequency band.
    pub receiver: CircularOutletReceiver,
    /// Numerical observation image and its exact clocks.
    pub reconstruction: DuctReconstruction,
}

/// One clocked receiver shared across wind exciters. Copying retains all history,
/// including a partially filled mechanical packet; no source is copied or stepped.
/// This remains prescribed-flow infinite-baffle radiation, not reaction, interface
/// transmission, gas mixing, arbitrary wall normals, or free-field body scattering.
#[derive(Clone)]
pub struct DuctPressure {
    source: DuctRadiationSpec,
    exterior: GasState,
    geometry: Vec<CircularOpening>,
    receiver: BaffledPressure,
    observer_gas: Option<GasPropagationReport>,
    trace: Vec<f64>,
    ratio: usize,
    filled: usize,
    observed: u64,
}
impl DuctPressure {
    /// Derive each opening radius from its actual section and use the existing
    /// signed Rayleigh receiver. No caller-specified alternative radius is accepted.
    /// All allocation is cold. The original mechanical terminal laws remain intact.
    ///
    /// # Errors
    /// Invalid/duplicate ports, missing regional state, incoherent gas, different
    /// static pressures, absent common exterior, or original receiver admission.
    pub fn new(networks: &[DuctView<'_>], source: DuctRadiationSpec) -> Result<Self, String> {
        if !(1..=MAX_RAYLEIGH_INPUTS).contains(&source.openings.len()) {
            return Err("duct radiation requires a bounded nonempty opening list".into());
        }
        if source.exterior.is_none() && source.openings.len() != 1 {
            return Err("multiple duct openings require an explicit common exterior gas".into());
        }
        let clock = source.reconstruction;
        let ratio = BaffledPressure::opening_clock_ratio(clock.mechanics_rate_hz(), clock.output_rate_hz())?;
        if let Some(gas) = source.exterior { gas.validate_coherence().map_err(|e| e.to_string())?; }
        let mut geometry = Vec::with_capacity(source.openings.len());
        let mut exterior = source.exterior;
        for (i, opening) in source.openings.iter().enumerate() {
            if source.openings[..i].iter().any(|p| p.network == opening.network && p.node == opening.node) {
                return Err("one physical duct port cannot radiate twice into the same receiver".into());
            }
            let view = networks.get(opening.network).ok_or("radiating duct network is absent")?;
            if view.gases.len() != view.graph.sections.len() {
                return Err("radiating duct requires one original gas per physical section".into());
            }
            if !matches!(view.graph.nodes.get(opening.node), Some(NetworkNode::Termination { .. }
                | NetworkNode::Impedance { .. } | NetworkNode::Relaxation { .. })) {
                return Err("radiation requires a passive degree-one physical terminal".into());
            }
            let mut connected = view.graph.sections.iter().enumerate()
                .filter(|(_, s)| s.nodes.contains(&opening.node));
            let (index, section) = connected.next().ok_or("radiating outlet has no section")?;
            if connected.next().is_some() || section.nodes[0] == section.nodes[1]
                || section.nodes.iter().any(|&n| n >= view.graph.nodes.len()) {
                return Err("radiating outlet requires exactly one nondegenerate adjacent section".into());
            }
            let local = view.gases[index];
            local.validate_coherence().map_err(|e| e.to_string())?;
            let gas = exterior.get_or_insert(local);
            if local.pressure.to_bits() != gas.pressure.to_bits() {
                return Err("duct and exterior static pressure differ; mean interface flow is not modeled".into());
            }
            geometry.push(CircularOpening { center_m: opening.center_m, radius_m: section.radius_m });
        }
        let exterior = exterior.expect("nonempty admitted ports");
        let medium = RayleighMedium { density: exterior.density, sound_speed: exterior.sound_speed };
        let receiver = match clock {
            DuctReconstruction::FiniteDifference { sample_rate_hz } =>
                BaffledPressure::circular_openings(&geometry, sample_rate_hz, source.receiver, medium),
            DuctReconstruction::Differentiating { mechanics_rate_hz, output_rate_hz } =>
                BaffledPressure::circular_openings_multirate_differentiating(
                    &geometry, mechanics_rate_hz, output_rate_hz, source.receiver, medium),
        }?;
        let trace = vec![0.0; source.openings.len() * ratio]; // <= 128 * 16
        Ok(Self {
            source,
            exterior,
            geometry,
            receiver,
            observer_gas: None,
            trace,
            ratio,
            filled: 0,
            observed: 0,
        })
    }

    /// Apply one explicit material propagation law to every admitted opening.
    /// Original geometry, signed port ordering, exterior and source/output clocks
    /// are retained. Absorption and dispersion act on each surface path before
    /// pressure is summed. The mechanical source and its terminal loads do not
    /// change. A cold differentiating receiver is required; a finite-difference
    /// receiver is never silently replaced by a different reconstruction.
    pub fn with_observer_gas(
        mut self,
        propagation: &GasPropagation,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, String> {
        if self.observed != 0 || self.observer_gas.is_some() {
            return Err("duct gas observation requires one cold attachment".into());
        }
        let DuctReconstruction::Differentiating {
            mechanics_rate_hz,
            output_rate_hz,
        } = self.source.reconstruction
        else {
            return Err(
                "duct gas observation requires explicit differentiating reconstruction".into(),
            );
        };
        let (receiver, report) = BaffledPressure::circular_openings_multirate_gas_differentiating(
            &self.geometry,
            mechanics_rate_hz,
            output_rate_hz,
            self.source.receiver,
            &self.exterior,
            propagation,
            cancelled,
        )?;
        if cancelled() {
            return Err("duct material observation cancelled".into());
        }
        self.receiver = receiver;
        self.observer_gas = Some(report);
        Ok(self)
    }

    /// Align independently prepared apertures without changing their physical
    /// relative flight times or dispersion. Only common numerical padding moves.
    pub(crate) fn align_gas_latency(&mut self, latency: usize) -> Result<(), String> {
        if self.observed != 0 {
            return Err("duct gas latency alignment requires a cold receiver".into());
        }
        let report = self
            .observer_gas
            .as_mut()
            .ok_or("duct gas latency alignment requires material propagation")?;
        self.receiver.align_gas_latency(report, latency)
    }

    /// Original port ordering, chosen environment and reconstruction controls.
    pub const fn source(&self) -> &DuctRadiationSpec { &self.source }
    /// Actual coherent exterior, including the explicit single-port local choice.
    pub const fn exterior_gas(&self) -> &GasState { &self.exterior }
    /// Shared kernel/history diagnostics; reading does not advance a clock.
    pub const fn receiver(&self) -> &BaffledPressure { &self.receiver }
    /// Original material and path-filter accounts, absent for default propagation.
    pub fn observer_gas(&self) -> Option<&GasPropagationReport> {
        self.observer_gas.as_ref()
    }
    /// Complete physical frames needed for one output frame.
    pub const fn mechanics_per_output(&self) -> usize { self.ratio }
    /// Successfully ingested physical samples, including the pending packet.
    pub const fn mechanical_samples(&self) -> u64 { self.observed }
    /// Completely emitted microphone samples.
    pub fn samples_rendered(&self) -> u64 { self.observed / self.ratio as u64 }
    /// Accepted physical frames awaiting the rest of their audio packet.
    pub const fn pending_frames(&self) -> usize { self.filled }

    /// Ingest all outward volume flows from ONE already accepted mechanics frame.
    /// The owner must call exactly once per source step. A control change between
    /// calls never resets the unfinished trace. None means the packet is incomplete,
    /// not silence. The returned pressure already includes every selected outlet.
    ///
    /// # Errors
    /// Missing/nonfinite flow, counter overflow or the existing receiver refusal.
    /// Refusal does not advance observation counters or any committed FIR history;
    /// scratch in the not-yet-accepted frame is overwritten on retry. A caller that
    /// has advanced its physical source must poison its enclosing performance.
    pub fn push_frame(&mut self, mut flow: impl FnMut(usize, usize) -> Option<f64>)
        -> Result<Option<f64>, String>
    {
        let next = self.observed.checked_add(1).ok_or("duct observation clock overflow")?;
        let ports = self.source.openings.len();
        let start = self.filled * ports;
        for (slot, port) in self.trace[start..start + ports].iter_mut().zip(&self.source.openings) {
            let value = flow(port.network, port.node).ok_or("accepted radiating port has no flow")?;
            if !value.is_finite() { return Err("radiating port flow must be finite".into()); }
            *slot = value;
        }
        let pressure = if self.filled + 1 == self.ratio {
            Some(self.receiver.step_trace(&self.trace)?)
        } else { None };
        self.observed = next;
        self.filled = if pressure.is_some() { 0 } else { self.filled + 1 };
        Ok(pressure)
    }
}
