//! Matching material-diffraction observation of accepted terminal flows.

use super::DgPortLoad;
use crate::bernoulli_aperture::radiation::environment::multiport::bake::MAX_BAFFLED_MIMO_PORTS;
use fs_bem::observation_transfers::{
    IterativeBoundaryObservationBake, LowerBandWindow, ReceiverSpec, TransferDeclaration,
    TransferError,
};
use fs_blake3::ContentHash;
use fs_duct::dg::PortId;
use fs_vfit::observation::ObservationFir;

/// One original causal pressure/gradient FIR bound to the complete aperture
/// order of a material-diffraction DG load. It observes accepted outward mean
/// flows; it never participates in reaction or changes the mechanical state.
///
/// Output order is receiver-major `[p, dp/dx, dp/dy, dp/dz]`, in Pa and Pa/m.
/// Original propagation delays, declared reconstruction latency, band tapers,
/// numerical fit error and unresolved spatial uncertainty remain those of the
/// supplied observation bake. Feeding interval-mean flows does not assert an
/// exact continuous reconstruction of a varying within-step boundary trace.
#[derive(Clone)]
pub struct DgBoundaryRadiationObserver {
    fir: ObservationFir,
    ports: Vec<PortId>,
    declaration_identity: [u8; 32],
    realization_identity: [u8; 32],
    transfer_declaration_hash: ContentHash,
    time_step_s: f64,
    observed_steps: u64,
    receivers: Vec<ReceiverSpec>,
    lower_window: Option<LowerBandWindow>,
    pass_band_upper_hz: f64,
    filter_latency_samples: usize,
}

impl DgPortLoad {
    /// Prepare the matching pressure observer BEFORE the first accepted step.
    /// The original immutable BEM bake and its complete declaration must agree
    /// with each other and with this exact geometry/material/source operator.
    /// A changed receiver needs its own physical bake; an unrelated or reordered
    /// aperture FIR cannot acquire this load's source identity.
    ///
    /// Observe each accepted step once through `DgBoundaryRadiationObserver::step`.
    /// Starting from an already moving load refuses because zero FIR prehistory
    /// would discard real past radiation. No BEM solve or refit is performed.
    pub fn boundary_radiation_observer(
        &self,
        declaration: &TransferDeclaration,
        bake: &IterativeBoundaryObservationBake,
    ) -> Result<DgBoundaryRadiationObserver, TransferError> {
        let source = self
            .boundary_radiation_record()
            .ok_or(TransferError::Invalid(
                "matching boundary observation requires an original material-diffraction source",
            ))?;
        declaration.admit_work()?;
        let transfer = bake.observation().transfer();
        let dt = self.propagator().dt();
        if self.state().step_count() != 0
            || declaration.producer_hash != source.operator_hash()
            || declaration.declaration_hash() != transfer.declaration_hash
            || declaration.sources != source.sources()
            || transfer.sources != source.sources()
            || declaration.receivers != transfer.receivers
            || self.terminal_ports().len() != source.sources().len()
            || self.terminal_ports().len() > MAX_BAFFLED_MIMO_PORTS
            || (1.0 / declaration.fs_hz).to_bits() != dt.to_bits()
            || transfer.fs_hz.to_bits() != declaration.fs_hz.to_bits()
            || source.time_step_s().to_bits() != dt.to_bits()
        {
            return Err(TransferError::Invalid(
                "boundary observation needs the same quiescent load, original operator, aperture order, declaration and clock",
            ));
        }
        Ok(DgBoundaryRadiationObserver {
            fir: bake.observer()?,
            ports: self.terminal_ports().to_vec(),
            declaration_identity: source.declaration_identity(),
            realization_identity: source.realization_identity(),
            transfer_declaration_hash: transfer.declaration_hash,
            time_step_s: dt,
            observed_steps: 0,
            receivers: transfer.receivers.clone(),
            lower_window: bake.observation().lower_band_window(),
            pass_band_upper_hz: transfer.pass_band_upper_hz,
            filter_latency_samples: transfer.filter_latency_samples,
        })
    }
}

impl DgBoundaryRadiationObserver {
    #[must_use]
    pub fn receivers(&self) -> &[ReceiverSpec] {
        &self.receivers
    }

    #[must_use]
    pub const fn transfer_declaration_hash(&self) -> ContentHash {
        self.transfer_declaration_hash
    }

    #[must_use]
    pub const fn lower_band_window(&self) -> Option<LowerBandWindow> {
        self.lower_window
    }

    #[must_use]
    pub const fn pass_band_upper_hz(&self) -> f64 {
        self.pass_band_upper_hz
    }

    #[must_use]
    pub const fn filter_latency_samples(&self) -> usize {
        self.filter_latency_samples
    }

    #[must_use]
    pub const fn observed_steps(&self) -> u64 {
        self.observed_steps
    }

    /// Observe the next accepted mechanical interval through the original FIR.
    /// Duplicate/skipped steps, changed sources or port order refuse without
    /// advancing observation history. Abandoned mechanical candidates have no
    /// observation. Finite-branch source flows come from their OUTER terminals.
    /// The original FIR owns its output allocation and transactional overflow
    /// refusal; the bounded input gather uses only stack storage.
    pub fn step(&mut self, load: &DgPortLoad) -> Result<Vec<f64>, TransferError> {
        let source = load
            .boundary_radiation_record()
            .ok_or(TransferError::Invalid(
                "boundary observer lost its material-diffraction source",
            ))?;
        let next = self
            .observed_steps
            .checked_add(1)
            .ok_or(TransferError::Invalid(
                "boundary observation step count overflow",
            ))?;
        let frames = load.last_terminal_frames();
        if load.state().step_count() != next
            || source.declaration_identity() != self.declaration_identity
            || source.realization_identity() != self.realization_identity
            || load.terminal_ports() != self.ports
            || load.propagator().dt().to_bits() != self.time_step_s.to_bits()
            || frames.len() != self.ports.len()
            || frames
                .iter()
                .zip(&self.ports)
                .any(|(frame, port)| frame.port != *port)
        {
            return Err(TransferError::Invalid(
                "boundary observation must consume each accepted interval exactly once from the same ordered source",
            ));
        }
        let mut flows = [0.0; MAX_BAFFLED_MIMO_PORTS];
        for (flow, frame) in flows.iter_mut().zip(frames) {
            *flow = frame.outward_flow_m3_s;
        }
        let output = self.fir.step(&flows[..frames.len()])?;
        self.observed_steps = next;
        Ok(output)
    }
}

#[cfg(test)]
mod tests;
