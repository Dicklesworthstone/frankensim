//! Original molecular edge laws closed against the existing full BEM exterior.
//!
//! Binding derives the graph and gas records from the same immutable radiation
//! owner. It cannot attach a similarly labelled but geometrically different
//! graph. Molecular/wall spectra are prepared once; each physical frequency
//! then reaches the shared loaded-source and BEM observation machinery.

use super::prepared::Propagation;
use super::*;
use fs_bem::observation_transfers::PreparedRadiationOptions;
use fs_duct::network::molecular::{MolecularEdge, MolecularNetwork};
use fs_duct::network::{CoupledSourceBank, PortSource};

/// Frequency-dependent interior storage, absorption and optional wall diffusion
/// attached to one original network/exterior owner. Neither binding nor a
/// frequency sweep changes its geometry, source records or exterior medium.
///
/// The exterior retains the original homogeneous, nondispersive BEM law. Only
/// the interior edges use the explicitly supplied molecular laws. This is a
/// settled harmonic model, without mean flow, evolving materials, switch-on
/// history or a causal time realization. All original plane-wave, opening-basis,
/// mesh, molecular support and diffusion-resolution limits remain binding.
pub struct MolecularRadiatingNetwork<'a> {
    model: &'a RadiatingNetwork,
    interior: MolecularNetwork,
}

impl RadiatingNetwork {
    /// Retain the original molecular reservoirs on every exact graph edge.
    /// This owner must have been built with [`Self::new_regional`], supplying
    /// explicit records even for a uniform gas. Each edge's relaxation state
    /// must match that original gas; labels alone do not establish identity.
    /// No source record, equilibrium replacement or exterior relaxation law is
    /// inferred from a numerical density/speed pair.
    ///
    /// Wall eigensolves occur here once, within the caller's original radial
    /// resolution and aggregate propagation-block limit. Later source-amplitude
    /// sweeps reuse the prepared transfer, and frequency sweeps reuse these same
    /// immutable edge spectra. The existing BEM budget separately bounds the
    /// complete exterior source basis at each frequency.
    ///
    /// # Errors
    /// Missing original regional records, incomplete/mismatched edge laws,
    /// unsupported embedded tone holes, exceeded original budgets or cancellation.
    pub fn bind_molecular(
        &self,
        edges: Vec<MolecularEdge>,
        maximum_blocks: usize,
        gate: &CancelGate,
    ) -> Result<MolecularRadiatingNetwork<'_>, ExteriorNetworkError> {
        poll(gate)?;
        let media = self.media.as_ref().ok_or_else(|| {
            bad("molecular radiation requires the original explicit regional gas records")
        })?;
        let interior = MolecularNetwork::new_controlled(
            self.driving.network().clone(),
            media.clone(),
            edges,
            maximum_blocks,
            &mut || gate.is_requested(),
        )?;
        poll(gate)?;
        Ok(MolecularRadiatingNetwork {
            model: self,
            interior,
        })
    }
}

impl<'a> MolecularRadiatingNetwork<'a> {
    /// Exact original graph, exterior mesh, openings and receiver binding.
    pub fn model(&self) -> &'a RadiatingNetwork {
        self.model
    }

    /// Exact original edge laws and already prepared molecular/wall spectra.
    pub fn molecular_network(&self) -> &MolecularNetwork {
        &self.interior
    }

    /// Prepare all declared scalar/coupled generators and loads at the actual
    /// angular frequency [rad/s]. The shared full BEM radiation matrix changes
    /// the opening flows before any microphone or field observation is formed.
    ///
    /// The returned existing prepared owner supports coherent drives, complete
    /// source statistics and directional microphones without additional solves.
    /// Source units, exp(-i omega t) phase, signed work and cancellation retain
    /// their existing semantics. The source/load partition must cover every
    /// non-opening port exactly once, as in [`RadiatingNetwork::prepare_sources`].
    pub fn prepare_sources<'p>(
        &'p self,
        omega: f64,
        sources: &'p [PortSource],
        bank: Option<&'p CoupledSourceBank>,
        scalar_loads: &'p [PortTermination],
        gate: &CancelGate,
    ) -> Result<PreparedRadiatingNetwork<'p>, ExteriorNetworkError> {
        self.model.prepare_sources_inner(
            omega,
            Propagation::Molecular(&self.interior),
            sources,
            bank,
            scalar_loads,
            None,
            gate,
        )
    }

    /// Additionally retain the existing positive angular quadrature and
    /// requested directions for the same BEM solutions and molecularly loaded
    /// opening flows. Far-field power is an observation of the existing load
    /// work; it is never added again to the network's dissipation.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_sources_with_radiation<'p>(
        &'p self,
        omega: f64,
        sources: &'p [PortSource],
        bank: Option<&'p CoupledSourceBank>,
        scalar_loads: &'p [PortTermination],
        directions: &[[f64; 3]],
        radiation: PreparedRadiationOptions,
        gate: &CancelGate,
    ) -> Result<PreparedRadiatingNetwork<'p>, ExteriorNetworkError> {
        self.model.prepare_sources_inner(
            omega,
            Propagation::Molecular(&self.interior),
            sources,
            bank,
            scalar_loads,
            Some((directions, radiation)),
            gate,
        )
    }
}
