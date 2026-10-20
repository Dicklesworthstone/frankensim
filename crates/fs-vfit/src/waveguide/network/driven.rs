//! Fixed external ports, with no topology or history replacement at closure.
use super::{NetworkFrame, NetworkNode, WaveguideError, WaveguideNetwork};

/// Old characteristic state exposed to the owner of a degree-one boundary.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DrivenBoundary {
    /// Physical node index in the unchanged original graph.
    pub node: usize,
    /// Original adjacent section index in the network's ordered section input.
    pub section: usize,
    /// Complete integer transit of that same retained section [samples].
    pub one_way_samples: usize,
    /// Original optional sub-sample transit; absent means integer travel.
    pub fractional_transit: Option<crate::waveguide::fractional::FractionalTransit>,
    /// Wave arriving from this endpoint's own adjacent section [Pa].
    pub incoming_pressure_pa: f64,
    /// That section's positive pressure/volume-flow impedance [Pa s/m³].
    pub impedance_pa_s_m3: f64,
}

/// One outgoing wave solved by an external physical boundary owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DrivenWave {
    /// Original degree-one `NetworkNode::Driven` index.
    pub node: usize,
    /// Wave departing from the boundary into the network [Pa].
    pub outgoing_pressure_pa: f64,
}

impl WaveguideNetwork {
    /// Current shared physical clock [s].
    #[must_use]
    pub const fn time_step_s(&self) -> f64 { self.time_step_s }

    /// Number of accepted wave/load updates, excluding previews and refusals.
    #[must_use]
    pub const fn accepted_steps(&self) -> u64 { self.accepted_steps }

    /// Additional external endpoints in their deterministic physical node order.
    pub fn driven_nodes(&self) -> impl Iterator<Item = usize> + '_ {
        self.nodes.iter().enumerate().filter_map(|(n, kind)|
            matches!(kind, NetworkNode::Driven).then_some(n))
    }

    /// Read this endpoint's old wave before solving any boundary. Positive
    /// transit on every section makes all such local solves simultaneous.
    pub fn driven_boundary(&self, node: usize) -> Result<DrivenBoundary, WaveguideError> {
        if self.nodes.get(node) != Some(&NetworkNode::Driven) {
            return Err(WaveguideError("driven boundary requires an explicitly driven node"));
        }
        let port = self.ports[self.offsets[node]];
        let line = &self.lines[port / 2];
        Ok(DrivenBoundary { node, section: port / 2, incoming_pressure_pa: self.arriving(port),
            one_way_samples: line.spec.one_way_samples, fractional_transit: line.fractional_transit(),
            impedance_pa_s_m3: line.spec.impedance_pa_s_m3 })
    }

    pub(super) fn validate_driven(&self, waves: &[DrivenWave]) -> Result<(), WaveguideError> {
        let mut supplied = waves.iter();
        for node in self.driven_nodes() {
            if supplied.next().is_none_or(|wave|
                wave.node != node || !wave.outgoing_pressure_pa.is_finite()) {
                return Err(WaveguideError("every driven endpoint needs one finite outgoing wave in node order"));
            }
        }
        if supplied.next().is_some() {
            return Err(WaveguideError("extra or repeated driven endpoint assignment"));
        }
        Ok(())
    }

    /// Validate simultaneous scattering and every load, then publish all old
    /// endpoints together. Refusal retains waves, reactive state and clock.
    pub fn step_with_driven(&mut self, outgoing: f64, waves: &[DrivenWave])
        -> Result<NetworkFrame, WaveguideError>
    {
        let frame = self.preview_step_with_driven(outgoing, waves)?;
        self.accept_preview();
        Ok(frame)
    }

    /// Publish two coupled networks with independently solved moving endpoints.
    /// Both original previews must pass before either wave, reactive/MIMO load,
    /// observation or clock advances. No accepted history is cloned or reset.
    /// Outgoing inlet waves and endpoint lists use self/other order.
    pub fn step_pair_with_driven(
        &mut self,
        other: &mut Self,
        outgoing: [f64; 2],
        waves: [&[DrivenWave]; 2],
    ) -> Result<[NetworkFrame; 2], WaveguideError> {
        if self.time_step_s.to_bits() != other.time_step_s.to_bits() {
            return Err(WaveguideError("coupled wave networks require identical clocks"));
        }
        let frames = [
            self.preview_step_with_driven(outgoing[0], waves[0])?,
            other.preview_step_with_driven(outgoing[1], waves[1])?,
        ];
        self.accept_preview();
        other.accept_preview();
        Ok(frames)
    }
}

#[cfg(test)]
mod tests;
