//! Interior fields and exterior observations of one BEM-loaded molecular solve.

use super::*;
use fs_duct::network::field::{NetworkFieldBudget, NetworkFieldProbe};
use fs_duct::network::molecular::{
    MolecularNetworkFieldRequest, PreparedMolecularNetworkFields,
};

/// Physical source/load coordinates and original interior probes.
/// The retained exterior owner supplies the complete opening load; no caller
/// can replace it with a second boundary or select an equilibrium loss model.
pub struct MolecularRadiatingFieldRequest<'a> {
    pub omega: f64,
    pub source_clock: &'a str,
    pub sources: &'a [PortSource],
    pub source_bank: Option<&'a CoupledSourceBank>,
    pub scalar_loads: &'a [PortTermination],
    pub probes: &'a [NetworkFieldProbe],
    pub budget: NetworkFieldBudget,
}

/// One native loaded solution supplies both sides of the original openings.
/// Interior fields retain their original material network and exact port bank;
/// exterior observation owns a bounded clone of those SAME solved transfers.
/// No second molecular evaluation, duct factorization or BEM solve is used.
pub struct PreparedMolecularRadiatingFields<'a> {
    interior: PreparedMolecularNetworkFields<'a>,
    radiating: PreparedRadiatingNetwork<'a>,
}

impl<'a> PreparedMolecularRadiatingFields<'a> {
    /// Original pressure/directed-flow fields, coherence, finite windows and
    /// source-owned cylinder losses. Gas absorption and mechanical wall work
    /// keep their native meanings and are not automatically caloric sources.
    pub fn interior(&self) -> &PreparedMolecularNetworkFields<'a> {
        &self.interior
    }

    /// Identically loaded opening flows, physical microphones and optional
    /// far field. Source order, units, clock and full mutual radiation loading
    /// are the same as the interior port bank.
    pub fn radiating(&self) -> &PreparedRadiatingNetwork<'a> {
        &self.radiating
    }
}

impl MolecularRadiatingNetwork<'_> {
    /// Prepare original interior fields and the same loaded exterior together.
    /// The native molecular owner resolves all sources, scalar loads and the
    /// full BEM opening matrix in one loaded factorization. The same solved
    /// opening flows feed the existing microphone/coherence consumers.
    ///
    /// NetworkFieldBudget bounds original probes and spatial work. The existing
    /// exterior field budget additionally admits four port-by-source arrays and
    /// one default-source vector for the retained native port-bank clone.
    /// All geometry, source support, conditioning and cancellation refusals
    /// preserve both original owners; no partial paired result is published.
    pub fn prepare_fields<'p>(
        &'p self,
        request: &MolecularRadiatingFieldRequest<'p>,
        gate: &CancelGate,
    ) -> Result<PreparedMolecularRadiatingFields<'p>, ExteriorNetworkError> {
        self.prepare_fields_inner(request, None, gate)
    }

    /// Also prepare the existing angular observation from the SAME unit-opening
    /// BEM solutions and native loaded transfer bank. Far-field power observes
    /// exterior load work and is never added again to an interior loss channel.
    pub fn prepare_fields_with_radiation<'p>(
        &'p self,
        request: &MolecularRadiatingFieldRequest<'p>,
        directions: &[[f64; 3]],
        radiation: PreparedRadiationOptions,
        gate: &CancelGate,
    ) -> Result<PreparedMolecularRadiatingFields<'p>, ExteriorNetworkError> {
        self.prepare_fields_inner(request, Some((directions, radiation)), gate)
    }

    fn prepare_fields_inner<'p>(
        &'p self,
        request: &MolecularRadiatingFieldRequest<'p>,
        angular: Option<(&[[f64; 3]], PreparedRadiationOptions)>,
        gate: &CancelGate,
    ) -> Result<PreparedMolecularRadiatingFields<'p>, ExteriorNetworkError> {
        let (radiating, interior) = self.model.prepare_sources_with_network(
            request.omega,
            Propagation::Molecular(&self.interior),
            request.sources,
            request.source_bank,
            request.scalar_loads,
            angular,
            true,
            |boundary| {
                let native = MolecularNetworkFieldRequest {
                    omega: request.omega,
                    source_clock: request.source_clock,
                    sources: request.sources,
                    source_bank: request.source_bank,
                    loads: request.scalar_loads,
                    coupled_load: Some(boundary),
                    probes: request.probes,
                    budget: request.budget,
                };
                let fields = self.interior.prepare_fields(&native, || gate.is_requested())?;
                let ports = fields.fields().port_response().clone();
                Ok((ports, fields))
            },
            gate,
        )?;
        Ok(PreparedMolecularRadiatingFields { interior, radiating })
    }
}
