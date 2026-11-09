//! Original molecular branch fields from the same loaded graph inverse.

use super::*;
use crate::molecular::local::absorption::MolecularCylinder;
use crate::molecular::spatial::PreparedMolecularSpatial;
use crate::network::field::{
    NetworkFieldBudget, NetworkFieldProbe, NetworkFieldRequest, PreparedNetworkFields,
};
use std::cell::OnceCell;

/// Fixed coherent sources/loads and original geometry probes. The molecular
/// network owns the constitutive choice; no independent equilibrium LossModel
/// can replace its molecular reservoirs or original wall diffusion spectra.
pub struct MolecularNetworkFieldRequest<'a> {
    pub omega: f64,
    pub source_clock: &'a str,
    pub sources: &'a [PortSource],
    pub source_bank: Option<&'a CoupledSourceBank>,
    pub loads: &'a [PortTermination],
    pub coupled_load: Option<&'a CoupledTermination>,
    pub probes: &'a [NetworkFieldProbe],
    pub budget: NetworkFieldBudget,
}

/// Source-retaining interior transfer bank on one physical frequency/clock.
/// Its native fields expose coherent samples, finite-window observations and
/// measured Welch cross spectra without another graph or constitutive solve.
/// Sound-speed metadata describes the retained equilibrium reference only;
/// actual phase, storage and loss use the original complex molecular law.
#[derive(Debug)]
pub struct PreparedMolecularNetworkFields<'network> {
    network: &'network MolecularNetwork,
    fields: PreparedNetworkFields,
    // Only cylinders with retained probes, bounded by MAX_FIELD_PROBES. The
    // exact Cp/wave are cached; material and radial spectra remain borrowed.
    cylinders: BTreeMap<String, BTreeMap<usize, MolecularCylinder<'network>>>,
}

impl PreparedMolecularNetworkFields<'_> {
    /// Original edge reservoirs, diffusion spectra, geometry and exterior.
    pub fn network(&self) -> &MolecularNetwork {
        self.network
    }

    /// Existing physical p/U observation bank. The pressure is a peak phasor;
    /// volume flow and mean particle velocity follow each original edge's
    /// from-to direction. Source columns, coherence and temporal quadrature
    /// retain the existing owners and their exact source-coordinate ordering.
    ///
    /// `prepare_cylinder_absorption` explicitly refuses: equilibrium Bessel
    /// loss forms are not the solved joint molecular/thermal constitutive law.
    /// Use this molecular wrapper's own `prepare_cylinder_absorption` instead.
    /// Signed power differences alone do not assign caloric heating.
    pub fn fields(&self) -> &PreparedNetworkFields {
        &self.fields
    }

    pub(super) fn cylinder(
        &self,
        edge: &str,
        segment: usize,
    ) -> Result<MolecularCylinder<'_>, DuctError> {
        self.cylinders
            .get(edge)
            .and_then(|edge| edge.get(&segment))
            .copied()
            .ok_or_else(|| bad("molecular absorption needs its original retained cylinder"))
    }
}

impl MolecularNetwork {
    /// Prepare the actual pressure, directed volume flow and cross-section
    /// mean particle velocity inside coherently loaded molecular branches.
    /// Original cylinders/cones support interior points; profiles support
    /// their original Magnus-panel boundaries. Branches, loops, source-bank
    /// coupling, external loads and open-circuit poles use one retained solve.
    ///
    /// Each edge evaluates its caloric response once, then reuses its original
    /// local waves in BOTH the loaded graph and observation reconstruction.
    /// Prefix/suffix scattering remains stable on opaque edges; cone probes
    /// split the original slice without resampling its material or radius.
    /// Preparation caches at most the network's admitted aggregate block cap;
    /// selected-edge transfer/probe limits retain NetworkFieldBudget. All
    /// later observations are contractions only. The synchronous dense solve
    /// remains bounded and is polled immediately before and after, not within.
    ///
    /// # Errors
    /// The original graph/material/source/load/probe/work/conditioning refusals
    /// and cancellation. No partial transfer bank or fallback law is returned.
    pub fn prepare_fields(
        &self,
        request: &MolecularNetworkFieldRequest<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMolecularNetworkFields<'_>, MolecularDuctError> {
        poll(&mut cancelled, 0)?;
        self.media.validate_for(&self.network)?;
        super::super::validate_medium(self.media.exterior_gas(), request.omega)?;
        let native_request = NetworkFieldRequest {
            omega: request.omega,
            // Only the common geometry/source carrier needs this field. All
            // propagation is supplied by the original molecular operators and
            // None below disables equilibrium volume-loss reconstruction.
            loss: crate::LossModel::Lossless,
            source_clock: request.source_clock,
            sources: request.sources,
            source_bank: request.source_bank,
            loads: request.loads,
            coupled_load: request.coupled_load,
            probes: request.probes,
            budget: request.budget,
        };
        // Filled only after the shared owner has admitted all probes, source
        // counts and transfer budgets. Subsequent closures borrow this one
        // frequency preparation; they neither re-query nor approximate it.
        let cached: OnceCell<BTreeMap<String, PreparedMolecularSpatial<'_>>> = OnceCell::new();
        let fields = PreparedNetworkFields::prepare_with_operators(
            &self.network,
            &native_request,
            None,
            |edge| {
                let gas = self.edges[&edge.id].state().gas();
                (gas.density, gas.sound_speed)
            },
            |_edge| None,
            |selected, cancelled| {
                let options = ComposeOptions {
                    maximum_condition: MAX_NETWORK_CONDITION,
                    maximum_blocks: self.blocks,
                    ..ComposeOptions::default()
                };
                let mut operators = BTreeMap::new();
                let mut original = BTreeMap::new();
                for (id, duct) in &self.edges {
                    poll(&mut || cancelled(), original.len())?;
                    let prepared =
                        duct.prepare_spatial_controlled(request.omega, &mut || cancelled())?;
                    operators.insert(
                        id.clone(),
                        compose_scattering(&prepared.chain.blocks, options, &mut || cancelled())?,
                    );
                    original.insert(id.clone(), prepared);
                }
                cached.set(original).map_err(|_| {
                    bad("molecular interior frequency preparation was already initialized")
                })?;
                poll(&mut || cancelled(), self.blocks)?;
                let result = loaded::prepare_spatial_with_operators(
                    &self.network,
                    request.omega,
                    request.sources,
                    request.source_bank,
                    request.loads,
                    request.coupled_load,
                    selected,
                    |edge| Ok::<_, MolecularDuctError>(operators[&edge.id]),
                    |load, node| {
                        Ok(super::super::load_impedance(
                            load,
                            super::super::node_gas(
                                &self.network,
                                self.media.exterior_gas(),
                                Some(&self.media),
                                node,
                            )
                            .as_ref(),
                            self.media.exterior_gas(),
                            request.omega,
                        )?)
                    },
                )?;
                poll(&mut || cancelled(), self.blocks)?;
                Ok::<_, MolecularDuctError>(result)
            },
            |edge, options, boundaries, cancelled| {
                poll(&mut || cancelled(), 0)?;
                let prepared = &cached.get().expect("loaded preparation precedes probes")[&edge.id];
                if prepared.chain.blocks.len() > options.maximum_blocks {
                    return Err(bad("molecular spatial propagation block budget exceeded").into());
                }
                boundaries.extend_from_slice(&prepared.segment_boundaries);
                Ok(crate::compose::PreparedChain {
                    blocks: prepared.chain.blocks.clone(),
                    reference_impedance: prepared.chain.reference_impedance,
                    min_shear_number: prepared.chain.min_shear_number,
                })
            },
            |edge, block, _radius| {
                cached.get().expect("loaded preparation precedes probes")[&edge.id].wave(block)
            },
            &mut cancelled,
        )?;
        let cached = cached
            .into_inner()
            .expect("loaded frequency preparation completed");
        let mut cylinders: BTreeMap<String, BTreeMap<usize, MolecularCylinder<'_>>> =
            BTreeMap::new();
        for (index, location) in fields.locations().iter().enumerate() {
            poll(&mut cancelled, index)?;
            let segment = match location.probe.section {
                crate::network::field::NetworkFieldSection::Inlet { segment }
                | crate::network::field::NetworkFieldSection::Outlet { segment }
                | crate::network::field::NetworkFieldSection::Cylinder { segment, .. } => segment,
                _ => continue,
            };
            if let Some(cylinder) = cached[&location.probe.edge_id].cylinder(segment) {
                cylinders
                    .entry(location.probe.edge_id.clone())
                    .or_default()
                    .entry(segment)
                    .or_insert(cylinder);
            }
        }
        poll(&mut cancelled, fields.locations().len())?;
        Ok(PreparedMolecularNetworkFields {
            network: self,
            fields,
            cylinders,
        })
    }
}
