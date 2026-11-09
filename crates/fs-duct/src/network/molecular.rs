//! Original molecular storage, absorption and wall diffusion in a duct graph.
//!
//! Each exact edge identity retains its own evaluated molecular reservoirs and
//! optional circular thermal/shear spectra and mechanical material wall.
//! Frequency queries use those original
//! two-ports in the existing fully loaded network solve. No equilibrium sound
//! speed, fitted attenuation or separately added thermal shunt replaces them.
//! The returned prepared source operator feeds the existing coherent-drive,
//! source-statistics and finite-window observation APIs without another solve.
//! [`MolecularNetwork::prepare_fields`] also observes original branch pressure,
//! flow and velocity using those same molecular waves and loaded equations.
//!
//! Geometry is frozen, with one gas state per edge, a common static pressure
//! and an independently declared exterior. Cylinders, cones, profiles, parallel
//! paths and loops retain the existing plane-wave geometry and graph owners.
//! Chimneys must be explicit graph branches with their own mouth boundaries;
//! embedded compact tone holes need additional exterior-port data and refuse.
//! JunctionCorrections changes only those compact holes; both policies are
//! equivalent on the admitted hole-free edge geometry.

use super::{
    CoupledSourceBank, CoupledTermination, DuctNetwork, MAX_NETWORK_CONDITION, NetworkMedia,
    PortSource, PortTermination, PreparedNetworkSources, bad, loaded,
};
use crate::compose::{ComposeOptions, MAX_COMPOSITION_BLOCKS, compose_scattering};
use crate::molecular::{
    MolecularCompliantWall, MolecularDuct, MolecularDuctError, MolecularWallResolution,
};
use crate::{ConeSlices, DuctError, Segment};
use fs_material::relaxation::RelaxationState;
use std::collections::BTreeMap;
use std::sync::Arc;

mod field;
pub use field::{MolecularNetworkFieldRequest, PreparedMolecularNetworkFields};
mod absorption;
pub use absorption::{
    MolecularNetworkCylinderAbsorption, PreparedMolecularNetworkCylinderAbsorption,
};

/// Original molecular law and explicit diffusion/geometry resolution for one
/// authored graph edge. The state's gas must equal that edge's original media
/// record; matching only density or equilibrium sound speed is insufficient.
#[derive(Debug, Clone)]
pub struct MolecularEdge {
    /// Exact authored graph edge identity, independent of vector position.
    pub edge_id: String,
    /// Original evaluated capacities, collision times and equilibrium gas.
    pub state: RelaxationState,
    /// None is an adiabatic slip wall; Some retains joint gas/wall diffusion.
    pub wall: Option<MolecularWallResolution>,
    /// Existing spherical-cone integration resolution, without a silent cap.
    pub cones: ConeSlices,
}

/// Immutable graph and original edge operators. Radial diffusion eigensolves
/// occur at binding; frequency sweeps reuse their positive spectra. A prepared
/// source solve can then be reused for arbitrary coherent source amplitudes.
#[derive(Debug, Clone)]
pub struct MolecularNetwork {
    network: DuctNetwork,
    media: NetworkMedia,
    edges: BTreeMap<String, MolecularDuct>,
    blocks: usize,
}

fn poll(cancelled: &mut impl FnMut() -> bool, completed: usize) -> Result<(), DuctError> {
    if cancelled() {
        Err(DuctError::Cancelled {
            completed_steps: completed,
        })
    } else {
        Ok(())
    }
}

impl MolecularNetwork {
    /// Bind every edge exactly once, retaining the complete original media.
    /// `maximum_blocks` caps the SUM of all actual edge propagation blocks.
    /// It never reduces the caller's cone, profile or diffusion resolution.
    /// Original mechanical pins on `NetworkMedia` enter the same molecular
    /// propagation through their native locally reacting surface law.
    pub fn new(
        network: DuctNetwork,
        media: NetworkMedia,
        edges: Vec<MolecularEdge>,
        maximum_blocks: usize,
    ) -> Result<Self, MolecularDuctError> {
        Self::new_controlled(network, media, edges, maximum_blocks, &mut || false)
    }

    /// Cancellable binding; no partial graph or fallback constitutive law is
    /// returned when an original material, resolution or work bound refuses.
    pub fn new_controlled(
        network: DuctNetwork,
        media: NetworkMedia,
        edges: Vec<MolecularEdge>,
        maximum_blocks: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MolecularDuctError> {
        Self::new_with_compliant_walls_controlled(
            network,
            media,
            edges,
            Vec::new(),
            maximum_blocks,
            cancelled,
        )
    }

    /// Bind source-owned local mechanical walls by exact edge identity.
    /// Unlisted edges use their original media pin, or remain rigid if absent.
    /// Gas reservoirs, gas/wall diffusion and geometry keep their native owner;
    /// no independent attenuation or extra equilibrium thermal shunt is added.
    ///
    /// A wall retains its original Arc, material state, memory response and
    /// frequency/radius admission. Its mechanical work may include structural
    /// loss or radiation and is never automatically classified as gas heat.
    /// The boundary is locally reacting: axial shell dynamics, mutual wall
    /// coupling, large motion and thermoelastic feedback are outside this image.
    ///
    /// # Errors
    /// Duplicate or unknown wall identity, a custom wall overlapping an
    /// existing media pin, and all original gas/geometry/work refusals.
    pub fn new_with_compliant_walls(
        network: DuctNetwork,
        media: NetworkMedia,
        edges: Vec<MolecularEdge>,
        walls: Vec<(String, MolecularCompliantWall)>,
        maximum_blocks: usize,
    ) -> Result<Self, MolecularDuctError> {
        Self::new_with_compliant_walls_controlled(
            network,
            media,
            edges,
            walls,
            maximum_blocks,
            &mut || false,
        )
    }

    /// Cancellable binding of the same original gas and mechanical owners.
    /// All edge identities and conflicting declarations are admitted before
    /// preparing radial spectra. Frequency-dependent walls admit at querying.
    pub fn new_with_compliant_walls_controlled(
        network: DuctNetwork,
        media: NetworkMedia,
        edges: Vec<MolecularEdge>,
        walls: Vec<(String, MolecularCompliantWall)>,
        maximum_blocks: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MolecularDuctError> {
        poll(cancelled, 0)?;
        media.validate_for(&network)?;
        if maximum_blocks == 0 || maximum_blocks > MAX_COMPOSITION_BLOCKS {
            return Err(bad("molecular network composition block cap is invalid").into());
        }
        if edges.len() != network.edges.len() {
            return Err(bad("molecular network needs exactly one law per edge identity").into());
        }
        let mut bindings = BTreeMap::new();
        for edge in edges {
            if bindings.insert(edge.edge_id.clone(), edge).is_some() {
                return Err(bad("molecular network edge identity is duplicated").into());
            }
        }
        if walls.len() > network.edges.len() {
            return Err(bad("molecular network mechanical wall count exceeds edge count").into());
        }
        let mut mechanical: BTreeMap<String, MolecularCompliantWall> = BTreeMap::new();
        for (id, wall) in walls {
            if !network.edges.iter().any(|edge| edge.id == id) {
                return Err(bad("molecular network wall names an unknown edge").into());
            }
            if media.edge_wall(&id).is_some() {
                return Err(bad("molecular network custom wall overlaps an original media pin").into());
            }
            if mechanical.insert(id, wall).is_some() {
                return Err(bad("molecular network mechanical wall identity is duplicated").into());
            }
        }
        // Admit identity and gas before preparing ANY expensive wall spectra.
        for edge in &network.edges {
            if let Some(pin) = media.edge_wall(&edge.id) {
                // WallPin's WallLaw implementation delegates to its original
                // fs-phs response. The complete media still retains the pin.
                mechanical.insert(edge.id.clone(), Arc::new(*pin));
            }
            let binding = bindings
                .get(&edge.id)
                .ok_or_else(|| bad("molecular network omits an authored edge identity"))?;
            if Some(binding.state.gas()) != media.edge_gas(&edge.id) {
                return Err(
                    bad("molecular network state differs from its original edge gas").into(),
                );
            }
            if edge
                .duct
                .segments
                .iter()
                .any(|s| matches!(s, Segment::ToneHole { .. }))
            {
                return Err(
                    bad("molecular network chimneys require explicit graph branches").into(),
                );
            }
        }
        let mut prepared = BTreeMap::new();
        let mut blocks = 0;
        // Stable identities make preparation/cancellation independent of the
        // authored edge-vector order, just like the original graph assembly.
        for (id, binding) in bindings {
            poll(cancelled, blocks)?;
            let edge = network
                .edges
                .iter()
                .find(|edge| edge.id == id)
                .expect("complete exact edge identities admitted");
            let duct = match mechanical.remove(&id) {
                Some(wall) => MolecularDuct::new_with_compliant_wall_controlled(
                    &edge.duct,
                    binding.state,
                    binding.wall,
                    wall,
                    binding.cones,
                    maximum_blocks - blocks,
                    cancelled,
                )?,
                None => MolecularDuct::new_controlled(
                    &edge.duct,
                    binding.state,
                    binding.wall,
                    binding.cones,
                    maximum_blocks - blocks,
                    cancelled,
                )?,
            };
            blocks += duct.blocks();
            prepared.insert(id, duct);
        }
        poll(cancelled, blocks)?;
        Ok(Self {
            network,
            media,
            edges: prepared,
            blocks,
        })
    }

    /// Original frozen topology and exact edge geometry.
    #[must_use]
    pub fn network(&self) -> &DuctNetwork {
        &self.network
    }

    /// Original source-bearing regional and exterior records.
    #[must_use]
    pub fn media(&self) -> &NetworkMedia {
        &self.media
    }

    /// Original prepared constitutive owner at an exact edge identity.
    #[must_use]
    pub fn edge(&self, id: &str) -> Option<&MolecularDuct> {
        self.edges.get(id)
    }

    /// Aggregate actual propagation resolution, without any truncation.
    #[must_use]
    pub fn blocks(&self) -> usize {
        self.blocks
    }

    /// Prepare the complete loaded graph at physical angular frequency [rad/s].
    /// Sources are complex peak amplitudes on the same exp(-i omega t) clock.
    /// Every external port needs an explicit source or scalar/coupled load.
    /// Full source-bank coupling, signed work and retained-wave pole handling
    /// belong to the existing network solver. No open-circuit Z is required.
    pub fn prepare_sources(
        &self,
        omega: f64,
        sources: &[PortSource],
        bank: Option<&CoupledSourceBank>,
        loads: &[PortTermination],
        coupled: Option<&CoupledTermination>,
    ) -> Result<PreparedNetworkSources, MolecularDuctError> {
        self.prepare_sources_controlled(omega, sources, bank, loads, coupled, &mut || false)
    }

    /// Poll at original local wave evaluations, scattering composition and
    /// around the existing bounded dense graph solve. Cancellation never
    /// publishes a partially prepared source operator.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_sources_controlled(
        &self,
        omega: f64,
        sources: &[PortSource],
        bank: Option<&CoupledSourceBank>,
        loads: &[PortTermination],
        coupled: Option<&CoupledTermination>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedNetworkSources, MolecularDuctError> {
        poll(cancelled, 0)?;
        super::validate_medium(self.media.exterior_gas(), omega)?;
        let options = ComposeOptions {
            maximum_condition: MAX_NETWORK_CONDITION,
            maximum_blocks: self.blocks,
            ..ComposeOptions::default()
        };
        let mut edges = BTreeMap::new();
        for (id, duct) in &self.edges {
            let chain = duct.prepared_chain_controlled(omega, cancelled)?;
            edges.insert(
                id.clone(),
                compose_scattering(&chain.blocks, options, cancelled)?,
            );
        }
        poll(cancelled, self.blocks)?;
        let prepared = loaded::prepare_sources_from_edges(
            &self.network,
            &self.media,
            omega,
            &edges,
            sources,
            bank,
            loads,
            coupled,
        )?;
        poll(cancelled, self.blocks)?;
        Ok(prepared)
    }
}

#[cfg(test)]
mod tests;
