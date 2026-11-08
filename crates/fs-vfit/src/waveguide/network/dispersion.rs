//! Authored lossless phase dispersion in the original network propagation.

use super::{NetworkNode, NetworkSegment, WaveguideError, WaveguideNetwork};
use crate::waveguide::{
    dispersion::{AllpassSection, DispersionStage, TransitResponse, stages},
    fractional::FractionalTransit,
};

impl WaveguideNetwork {
    /// Complete requested network payload including every bidirectional
    /// cascade. Coefficients are one borrowed slice per original section;
    /// caller-owned inputs and allocator overhead are excluded.
    pub fn required_dispersive_memory_bytes(
        node_count: usize, segments: &[NetworkSegment], sections: &[&[AllpassSection]],
    ) -> Result<usize, WaveguideError> {
        if sections.len() != segments.len() {
            return Err(WaveguideError("one dispersion cascade is required for every network section"));
        }
        let mut bytes = Self::required_memory_bytes(node_count, segments)?;
        for cascade in sections {
            bytes = cascade.len().checked_mul(core::mem::size_of::<DispersionStage>())
                .and_then(|extra| bytes.checked_add(extra))
                .ok_or(WaveguideError("network dispersion payload overflows"))?;
        }
        Ok(bytes)
    }

    /// Build fixed cascades in both directions of every original section.
    /// The same junctions, reactive loads, driven/MIMO ports and paired
    /// transaction own all subsequent propagation. Dispersion storage is wave
    /// storage, never dissipation. No live reset, retuning or extra graph exists.
    pub fn with_dispersive_transits(
        nodes: &[NetworkNode], segments: &[NetworkSegment], fractions: &[FractionalTransit],
        sections: &[&[AllpassSection]], time_step_s: f64, max_memory_bytes: usize,
    ) -> Result<Self, WaveguideError> {
        if Self::required_dispersive_memory_bytes(nodes.len(), segments, sections)? > max_memory_bytes {
            return Err(WaveguideError("network dispersion payload exceeds memory budget"));
        }
        let mut network = Self::with_fractional_transits(nodes, segments, fractions, time_step_s, max_memory_bytes)?;
        for (line, cascade) in network.lines.iter_mut().zip(sections) {
            line.dispersion = stages(cascade)?;
        }
        Ok(network)
    }

    /// Immutable actual section states and coefficients, in cascade order.
    #[must_use]
    pub fn section_dispersion_stages(&self, section: usize) -> Option<&[DispersionStage]> {
        self.lines.get(section).map(|line| line.dispersion_stages())
    }

    /// Complete one-way discrete phase/group transit of an original section.
    pub fn section_propagation_response(
        &self, section: usize, omega: f64,
    ) -> Result<TransitResponse, WaveguideError> {
        self.lines.get(section).ok_or(WaveguideError("unknown waveguide section"))?
            .propagation_response(omega)
    }
}
