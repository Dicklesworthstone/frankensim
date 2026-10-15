//! Fractional sections on the original connected wave/load owner.
use super::{NetworkNode, NetworkSegment, WaveguideError, WaveguideNetwork};
use crate::waveguide::fractional::{FractionalPair, FractionalTransit};

impl WaveguideNetwork {
    /// Build a zero-state network with one explicitly selected residual transit
    /// per section. Every load, pressure junction, paired-network transaction and
    /// buffer remains in the original owner. Fractional storage is included in
    /// wave_stored_energy_j; it is not terminal loss or a pressure gain.
    ///
    /// Allpass coordinates are inline in PassiveWaveguide, so the original
    /// required_memory_bytes includes their complete payload. No per-step
    /// allocation, extra junction delay or second propagation graph is added.
    /// The original constructor remains the exact integer-transit image.
    pub fn with_fractional_transits(nodes: &[NetworkNode], segments: &[NetworkSegment],
        transits: &[FractionalTransit], time_step_s: f64, max_memory_bytes: usize)
        -> Result<Self, WaveguideError>
    {
        if transits.len()!=segments.len() {
            return Err(WaveguideError("one fractional transit is required for every network section"));
        }
        let mut network=Self::new(nodes,segments,time_step_s,max_memory_bytes)?;
        for (line,transit) in network.lines.iter_mut().zip(transits) {
            line.fractional=FractionalPair::new(*transit);
        }
        Ok(network)
    }
    /// Residual transit of an existing section; zero identifies integer travel.
    /// The complete DC transit is its integer delay plus this residual.
    #[must_use]
    pub fn section_fractional_transit(&self, section: usize) -> Option<FractionalTransit> {
        self.lines.get(section).map(|line|line.fractional_transit()
            .unwrap_or_else(||FractionalTransit::new(0.0).expect("exact zero transit")))
    }
    /// Actual normalized memories in the original section and direction order.
    #[must_use]
    pub fn section_fractional_state_pa(&self, section: usize) -> Option<[f64;2]> {
        self.lines.get(section).map(|line|line.fractional_state_pa())
    }
}
