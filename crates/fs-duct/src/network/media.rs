//! Material records on the existing pressure/volume-flow graph.
//!
//! Each edge has one frozen gas. Every exposed boundary sees the explicitly
//! supplied homogeneous exterior, including chimneys embedded in an edge.
//! Junction pressure and volume flow stay in the original nodal solver. No
//! mixture, temperature, density or wave impedance is averaged at a junction.

use super::{DuctError, DuctNetwork, bad};
use fs_material::gas::GasState;
use fs_material::gas_species::input::GasRecord;
use std::collections::BTreeMap;
use std::sync::Arc;

/// Complete immutable gas binding for an existing duct graph.
///
/// Edge identities, rather than array positions, bind the records, so graph
/// reindexing cannot move a material onto another branch. The original records
/// retain thermochemical source identities, query states and composition where
/// provided. A numerical/legacy record receives no additional source authority.
///
/// This is prescribed, no-mean-flow acoustics at one common static pressure.
/// It does not evolve mixing, diffusion, heat or a moving gas interface.
#[derive(Clone, Debug, PartialEq)]
pub struct NetworkMedia {
    // Frequency sweeps retain the original records without copying all source
    // receipts into every NetworkResponse. Only immutable handles are shared.
    edges: Arc<BTreeMap<String, GasRecord>>,
    exterior: Arc<GasRecord>,
}

impl NetworkMedia {
    /// Bind every edge exactly once and declare one exterior atmosphere.
    /// Exterior temperature/composition may differ from every interior edge.
    ///
    /// # Errors
    /// Invalid graph, missing/extra/duplicate edge identity, incoherent gas or
    /// unequal static pressure. No default gas or inferred composition exists.
    pub fn new(
        net: &DuctNetwork,
        edges: Vec<(String, GasRecord)>,
        exterior: GasRecord,
    ) -> Result<Self, DuctError> {
        net.validate()?;
        if edges.len() != net.edges.len() {
            return Err(bad(
                "network media require exactly one record per edge identity",
            ));
        }
        let mut indexed = BTreeMap::new();
        for (id, record) in edges {
            if indexed.insert(id, record).is_some() {
                return Err(bad("network gas edge identity is duplicated"));
            }
        }
        let result = Self {
            edges: Arc::new(indexed),
            exterior: Arc::new(exterior),
        };
        result.validate_for(net)?;
        Ok(result)
    }

    /// Verify complete identity coverage before using the binding with a graph.
    /// The topology may be reordered; every actual edge still uses its own id.
    pub fn validate_for(&self, net: &DuctNetwork) -> Result<(), DuctError> {
        net.validate()?;
        if self.edges.len() != net.edges.len() {
            return Err(bad(
                "network media do not cover the complete edge identity set",
            ));
        }
        crate::validate_gas(self.exterior_gas())?;
        for edge in &net.edges {
            let gas = self
                .edge_gas(&edge.id)
                .ok_or_else(|| bad("network media omit an authored edge identity"))?;
            crate::exterior::admit(gas, self.exterior_gas())?;
        }
        Ok(())
    }

    /// Complete original record at an exact edge identity.
    #[must_use]
    pub fn edge_record(&self, id: &str) -> Option<&GasRecord> {
        self.edges.get(id)
    }

    /// Exact numerical view of the original edge record, without reevaluation.
    #[must_use]
    pub fn edge_gas(&self, id: &str) -> Option<&GasState> {
        self.edge_record(id).map(GasRecord::state)
    }

    /// Original exterior record, including its source/composition authority.
    #[must_use]
    pub fn exterior_record(&self) -> &GasRecord {
        &self.exterior
    }

    /// Homogeneous exterior used by radiation, porous boundaries and observers.
    #[must_use]
    pub fn exterior_gas(&self) -> &GasState {
        self.exterior.state()
    }
}
